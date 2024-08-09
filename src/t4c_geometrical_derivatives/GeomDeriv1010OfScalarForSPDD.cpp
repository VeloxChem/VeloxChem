#include "GeomDeriv1010OfScalarForSPDD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1010_spdd_0(CSimdArray<double>& buffer_1010_spdd,
                     const CSimdArray<double>& buffer_pppd,
                     const CSimdArray<double>& buffer_ppfd,
                     const double a_exp,
                     const double* c_exps) -> void
{
    const auto ndims = buffer_1010_spdd.number_of_columns();

    /// Set up components of auxilary buffer : buffer_pppd

    auto g_x_x_x_xx = buffer_pppd[0];

    auto g_x_x_x_xy = buffer_pppd[1];

    auto g_x_x_x_xz = buffer_pppd[2];

    auto g_x_x_x_yy = buffer_pppd[3];

    auto g_x_x_x_yz = buffer_pppd[4];

    auto g_x_x_x_zz = buffer_pppd[5];

    auto g_x_x_y_xx = buffer_pppd[6];

    auto g_x_x_y_xy = buffer_pppd[7];

    auto g_x_x_y_xz = buffer_pppd[8];

    auto g_x_x_y_yy = buffer_pppd[9];

    auto g_x_x_y_yz = buffer_pppd[10];

    auto g_x_x_y_zz = buffer_pppd[11];

    auto g_x_x_z_xx = buffer_pppd[12];

    auto g_x_x_z_xy = buffer_pppd[13];

    auto g_x_x_z_xz = buffer_pppd[14];

    auto g_x_x_z_yy = buffer_pppd[15];

    auto g_x_x_z_yz = buffer_pppd[16];

    auto g_x_x_z_zz = buffer_pppd[17];

    auto g_x_y_x_xx = buffer_pppd[18];

    auto g_x_y_x_xy = buffer_pppd[19];

    auto g_x_y_x_xz = buffer_pppd[20];

    auto g_x_y_x_yy = buffer_pppd[21];

    auto g_x_y_x_yz = buffer_pppd[22];

    auto g_x_y_x_zz = buffer_pppd[23];

    auto g_x_y_y_xx = buffer_pppd[24];

    auto g_x_y_y_xy = buffer_pppd[25];

    auto g_x_y_y_xz = buffer_pppd[26];

    auto g_x_y_y_yy = buffer_pppd[27];

    auto g_x_y_y_yz = buffer_pppd[28];

    auto g_x_y_y_zz = buffer_pppd[29];

    auto g_x_y_z_xx = buffer_pppd[30];

    auto g_x_y_z_xy = buffer_pppd[31];

    auto g_x_y_z_xz = buffer_pppd[32];

    auto g_x_y_z_yy = buffer_pppd[33];

    auto g_x_y_z_yz = buffer_pppd[34];

    auto g_x_y_z_zz = buffer_pppd[35];

    auto g_x_z_x_xx = buffer_pppd[36];

    auto g_x_z_x_xy = buffer_pppd[37];

    auto g_x_z_x_xz = buffer_pppd[38];

    auto g_x_z_x_yy = buffer_pppd[39];

    auto g_x_z_x_yz = buffer_pppd[40];

    auto g_x_z_x_zz = buffer_pppd[41];

    auto g_x_z_y_xx = buffer_pppd[42];

    auto g_x_z_y_xy = buffer_pppd[43];

    auto g_x_z_y_xz = buffer_pppd[44];

    auto g_x_z_y_yy = buffer_pppd[45];

    auto g_x_z_y_yz = buffer_pppd[46];

    auto g_x_z_y_zz = buffer_pppd[47];

    auto g_x_z_z_xx = buffer_pppd[48];

    auto g_x_z_z_xy = buffer_pppd[49];

    auto g_x_z_z_xz = buffer_pppd[50];

    auto g_x_z_z_yy = buffer_pppd[51];

    auto g_x_z_z_yz = buffer_pppd[52];

    auto g_x_z_z_zz = buffer_pppd[53];

    auto g_y_x_x_xx = buffer_pppd[54];

    auto g_y_x_x_xy = buffer_pppd[55];

    auto g_y_x_x_xz = buffer_pppd[56];

    auto g_y_x_x_yy = buffer_pppd[57];

    auto g_y_x_x_yz = buffer_pppd[58];

    auto g_y_x_x_zz = buffer_pppd[59];

    auto g_y_x_y_xx = buffer_pppd[60];

    auto g_y_x_y_xy = buffer_pppd[61];

    auto g_y_x_y_xz = buffer_pppd[62];

    auto g_y_x_y_yy = buffer_pppd[63];

    auto g_y_x_y_yz = buffer_pppd[64];

    auto g_y_x_y_zz = buffer_pppd[65];

    auto g_y_x_z_xx = buffer_pppd[66];

    auto g_y_x_z_xy = buffer_pppd[67];

    auto g_y_x_z_xz = buffer_pppd[68];

    auto g_y_x_z_yy = buffer_pppd[69];

    auto g_y_x_z_yz = buffer_pppd[70];

    auto g_y_x_z_zz = buffer_pppd[71];

    auto g_y_y_x_xx = buffer_pppd[72];

    auto g_y_y_x_xy = buffer_pppd[73];

    auto g_y_y_x_xz = buffer_pppd[74];

    auto g_y_y_x_yy = buffer_pppd[75];

    auto g_y_y_x_yz = buffer_pppd[76];

    auto g_y_y_x_zz = buffer_pppd[77];

    auto g_y_y_y_xx = buffer_pppd[78];

    auto g_y_y_y_xy = buffer_pppd[79];

    auto g_y_y_y_xz = buffer_pppd[80];

    auto g_y_y_y_yy = buffer_pppd[81];

    auto g_y_y_y_yz = buffer_pppd[82];

    auto g_y_y_y_zz = buffer_pppd[83];

    auto g_y_y_z_xx = buffer_pppd[84];

    auto g_y_y_z_xy = buffer_pppd[85];

    auto g_y_y_z_xz = buffer_pppd[86];

    auto g_y_y_z_yy = buffer_pppd[87];

    auto g_y_y_z_yz = buffer_pppd[88];

    auto g_y_y_z_zz = buffer_pppd[89];

    auto g_y_z_x_xx = buffer_pppd[90];

    auto g_y_z_x_xy = buffer_pppd[91];

    auto g_y_z_x_xz = buffer_pppd[92];

    auto g_y_z_x_yy = buffer_pppd[93];

    auto g_y_z_x_yz = buffer_pppd[94];

    auto g_y_z_x_zz = buffer_pppd[95];

    auto g_y_z_y_xx = buffer_pppd[96];

    auto g_y_z_y_xy = buffer_pppd[97];

    auto g_y_z_y_xz = buffer_pppd[98];

    auto g_y_z_y_yy = buffer_pppd[99];

    auto g_y_z_y_yz = buffer_pppd[100];

    auto g_y_z_y_zz = buffer_pppd[101];

    auto g_y_z_z_xx = buffer_pppd[102];

    auto g_y_z_z_xy = buffer_pppd[103];

    auto g_y_z_z_xz = buffer_pppd[104];

    auto g_y_z_z_yy = buffer_pppd[105];

    auto g_y_z_z_yz = buffer_pppd[106];

    auto g_y_z_z_zz = buffer_pppd[107];

    auto g_z_x_x_xx = buffer_pppd[108];

    auto g_z_x_x_xy = buffer_pppd[109];

    auto g_z_x_x_xz = buffer_pppd[110];

    auto g_z_x_x_yy = buffer_pppd[111];

    auto g_z_x_x_yz = buffer_pppd[112];

    auto g_z_x_x_zz = buffer_pppd[113];

    auto g_z_x_y_xx = buffer_pppd[114];

    auto g_z_x_y_xy = buffer_pppd[115];

    auto g_z_x_y_xz = buffer_pppd[116];

    auto g_z_x_y_yy = buffer_pppd[117];

    auto g_z_x_y_yz = buffer_pppd[118];

    auto g_z_x_y_zz = buffer_pppd[119];

    auto g_z_x_z_xx = buffer_pppd[120];

    auto g_z_x_z_xy = buffer_pppd[121];

    auto g_z_x_z_xz = buffer_pppd[122];

    auto g_z_x_z_yy = buffer_pppd[123];

    auto g_z_x_z_yz = buffer_pppd[124];

    auto g_z_x_z_zz = buffer_pppd[125];

    auto g_z_y_x_xx = buffer_pppd[126];

    auto g_z_y_x_xy = buffer_pppd[127];

    auto g_z_y_x_xz = buffer_pppd[128];

    auto g_z_y_x_yy = buffer_pppd[129];

    auto g_z_y_x_yz = buffer_pppd[130];

    auto g_z_y_x_zz = buffer_pppd[131];

    auto g_z_y_y_xx = buffer_pppd[132];

    auto g_z_y_y_xy = buffer_pppd[133];

    auto g_z_y_y_xz = buffer_pppd[134];

    auto g_z_y_y_yy = buffer_pppd[135];

    auto g_z_y_y_yz = buffer_pppd[136];

    auto g_z_y_y_zz = buffer_pppd[137];

    auto g_z_y_z_xx = buffer_pppd[138];

    auto g_z_y_z_xy = buffer_pppd[139];

    auto g_z_y_z_xz = buffer_pppd[140];

    auto g_z_y_z_yy = buffer_pppd[141];

    auto g_z_y_z_yz = buffer_pppd[142];

    auto g_z_y_z_zz = buffer_pppd[143];

    auto g_z_z_x_xx = buffer_pppd[144];

    auto g_z_z_x_xy = buffer_pppd[145];

    auto g_z_z_x_xz = buffer_pppd[146];

    auto g_z_z_x_yy = buffer_pppd[147];

    auto g_z_z_x_yz = buffer_pppd[148];

    auto g_z_z_x_zz = buffer_pppd[149];

    auto g_z_z_y_xx = buffer_pppd[150];

    auto g_z_z_y_xy = buffer_pppd[151];

    auto g_z_z_y_xz = buffer_pppd[152];

    auto g_z_z_y_yy = buffer_pppd[153];

    auto g_z_z_y_yz = buffer_pppd[154];

    auto g_z_z_y_zz = buffer_pppd[155];

    auto g_z_z_z_xx = buffer_pppd[156];

    auto g_z_z_z_xy = buffer_pppd[157];

    auto g_z_z_z_xz = buffer_pppd[158];

    auto g_z_z_z_yy = buffer_pppd[159];

    auto g_z_z_z_yz = buffer_pppd[160];

    auto g_z_z_z_zz = buffer_pppd[161];

    /// Set up components of auxilary buffer : buffer_ppfd

    auto g_x_x_xxx_xx = buffer_ppfd[0];

    auto g_x_x_xxx_xy = buffer_ppfd[1];

    auto g_x_x_xxx_xz = buffer_ppfd[2];

    auto g_x_x_xxx_yy = buffer_ppfd[3];

    auto g_x_x_xxx_yz = buffer_ppfd[4];

    auto g_x_x_xxx_zz = buffer_ppfd[5];

    auto g_x_x_xxy_xx = buffer_ppfd[6];

    auto g_x_x_xxy_xy = buffer_ppfd[7];

    auto g_x_x_xxy_xz = buffer_ppfd[8];

    auto g_x_x_xxy_yy = buffer_ppfd[9];

    auto g_x_x_xxy_yz = buffer_ppfd[10];

    auto g_x_x_xxy_zz = buffer_ppfd[11];

    auto g_x_x_xxz_xx = buffer_ppfd[12];

    auto g_x_x_xxz_xy = buffer_ppfd[13];

    auto g_x_x_xxz_xz = buffer_ppfd[14];

    auto g_x_x_xxz_yy = buffer_ppfd[15];

    auto g_x_x_xxz_yz = buffer_ppfd[16];

    auto g_x_x_xxz_zz = buffer_ppfd[17];

    auto g_x_x_xyy_xx = buffer_ppfd[18];

    auto g_x_x_xyy_xy = buffer_ppfd[19];

    auto g_x_x_xyy_xz = buffer_ppfd[20];

    auto g_x_x_xyy_yy = buffer_ppfd[21];

    auto g_x_x_xyy_yz = buffer_ppfd[22];

    auto g_x_x_xyy_zz = buffer_ppfd[23];

    auto g_x_x_xyz_xx = buffer_ppfd[24];

    auto g_x_x_xyz_xy = buffer_ppfd[25];

    auto g_x_x_xyz_xz = buffer_ppfd[26];

    auto g_x_x_xyz_yy = buffer_ppfd[27];

    auto g_x_x_xyz_yz = buffer_ppfd[28];

    auto g_x_x_xyz_zz = buffer_ppfd[29];

    auto g_x_x_xzz_xx = buffer_ppfd[30];

    auto g_x_x_xzz_xy = buffer_ppfd[31];

    auto g_x_x_xzz_xz = buffer_ppfd[32];

    auto g_x_x_xzz_yy = buffer_ppfd[33];

    auto g_x_x_xzz_yz = buffer_ppfd[34];

    auto g_x_x_xzz_zz = buffer_ppfd[35];

    auto g_x_x_yyy_xx = buffer_ppfd[36];

    auto g_x_x_yyy_xy = buffer_ppfd[37];

    auto g_x_x_yyy_xz = buffer_ppfd[38];

    auto g_x_x_yyy_yy = buffer_ppfd[39];

    auto g_x_x_yyy_yz = buffer_ppfd[40];

    auto g_x_x_yyy_zz = buffer_ppfd[41];

    auto g_x_x_yyz_xx = buffer_ppfd[42];

    auto g_x_x_yyz_xy = buffer_ppfd[43];

    auto g_x_x_yyz_xz = buffer_ppfd[44];

    auto g_x_x_yyz_yy = buffer_ppfd[45];

    auto g_x_x_yyz_yz = buffer_ppfd[46];

    auto g_x_x_yyz_zz = buffer_ppfd[47];

    auto g_x_x_yzz_xx = buffer_ppfd[48];

    auto g_x_x_yzz_xy = buffer_ppfd[49];

    auto g_x_x_yzz_xz = buffer_ppfd[50];

    auto g_x_x_yzz_yy = buffer_ppfd[51];

    auto g_x_x_yzz_yz = buffer_ppfd[52];

    auto g_x_x_yzz_zz = buffer_ppfd[53];

    auto g_x_x_zzz_xx = buffer_ppfd[54];

    auto g_x_x_zzz_xy = buffer_ppfd[55];

    auto g_x_x_zzz_xz = buffer_ppfd[56];

    auto g_x_x_zzz_yy = buffer_ppfd[57];

    auto g_x_x_zzz_yz = buffer_ppfd[58];

    auto g_x_x_zzz_zz = buffer_ppfd[59];

    auto g_x_y_xxx_xx = buffer_ppfd[60];

    auto g_x_y_xxx_xy = buffer_ppfd[61];

    auto g_x_y_xxx_xz = buffer_ppfd[62];

    auto g_x_y_xxx_yy = buffer_ppfd[63];

    auto g_x_y_xxx_yz = buffer_ppfd[64];

    auto g_x_y_xxx_zz = buffer_ppfd[65];

    auto g_x_y_xxy_xx = buffer_ppfd[66];

    auto g_x_y_xxy_xy = buffer_ppfd[67];

    auto g_x_y_xxy_xz = buffer_ppfd[68];

    auto g_x_y_xxy_yy = buffer_ppfd[69];

    auto g_x_y_xxy_yz = buffer_ppfd[70];

    auto g_x_y_xxy_zz = buffer_ppfd[71];

    auto g_x_y_xxz_xx = buffer_ppfd[72];

    auto g_x_y_xxz_xy = buffer_ppfd[73];

    auto g_x_y_xxz_xz = buffer_ppfd[74];

    auto g_x_y_xxz_yy = buffer_ppfd[75];

    auto g_x_y_xxz_yz = buffer_ppfd[76];

    auto g_x_y_xxz_zz = buffer_ppfd[77];

    auto g_x_y_xyy_xx = buffer_ppfd[78];

    auto g_x_y_xyy_xy = buffer_ppfd[79];

    auto g_x_y_xyy_xz = buffer_ppfd[80];

    auto g_x_y_xyy_yy = buffer_ppfd[81];

    auto g_x_y_xyy_yz = buffer_ppfd[82];

    auto g_x_y_xyy_zz = buffer_ppfd[83];

    auto g_x_y_xyz_xx = buffer_ppfd[84];

    auto g_x_y_xyz_xy = buffer_ppfd[85];

    auto g_x_y_xyz_xz = buffer_ppfd[86];

    auto g_x_y_xyz_yy = buffer_ppfd[87];

    auto g_x_y_xyz_yz = buffer_ppfd[88];

    auto g_x_y_xyz_zz = buffer_ppfd[89];

    auto g_x_y_xzz_xx = buffer_ppfd[90];

    auto g_x_y_xzz_xy = buffer_ppfd[91];

    auto g_x_y_xzz_xz = buffer_ppfd[92];

    auto g_x_y_xzz_yy = buffer_ppfd[93];

    auto g_x_y_xzz_yz = buffer_ppfd[94];

    auto g_x_y_xzz_zz = buffer_ppfd[95];

    auto g_x_y_yyy_xx = buffer_ppfd[96];

    auto g_x_y_yyy_xy = buffer_ppfd[97];

    auto g_x_y_yyy_xz = buffer_ppfd[98];

    auto g_x_y_yyy_yy = buffer_ppfd[99];

    auto g_x_y_yyy_yz = buffer_ppfd[100];

    auto g_x_y_yyy_zz = buffer_ppfd[101];

    auto g_x_y_yyz_xx = buffer_ppfd[102];

    auto g_x_y_yyz_xy = buffer_ppfd[103];

    auto g_x_y_yyz_xz = buffer_ppfd[104];

    auto g_x_y_yyz_yy = buffer_ppfd[105];

    auto g_x_y_yyz_yz = buffer_ppfd[106];

    auto g_x_y_yyz_zz = buffer_ppfd[107];

    auto g_x_y_yzz_xx = buffer_ppfd[108];

    auto g_x_y_yzz_xy = buffer_ppfd[109];

    auto g_x_y_yzz_xz = buffer_ppfd[110];

    auto g_x_y_yzz_yy = buffer_ppfd[111];

    auto g_x_y_yzz_yz = buffer_ppfd[112];

    auto g_x_y_yzz_zz = buffer_ppfd[113];

    auto g_x_y_zzz_xx = buffer_ppfd[114];

    auto g_x_y_zzz_xy = buffer_ppfd[115];

    auto g_x_y_zzz_xz = buffer_ppfd[116];

    auto g_x_y_zzz_yy = buffer_ppfd[117];

    auto g_x_y_zzz_yz = buffer_ppfd[118];

    auto g_x_y_zzz_zz = buffer_ppfd[119];

    auto g_x_z_xxx_xx = buffer_ppfd[120];

    auto g_x_z_xxx_xy = buffer_ppfd[121];

    auto g_x_z_xxx_xz = buffer_ppfd[122];

    auto g_x_z_xxx_yy = buffer_ppfd[123];

    auto g_x_z_xxx_yz = buffer_ppfd[124];

    auto g_x_z_xxx_zz = buffer_ppfd[125];

    auto g_x_z_xxy_xx = buffer_ppfd[126];

    auto g_x_z_xxy_xy = buffer_ppfd[127];

    auto g_x_z_xxy_xz = buffer_ppfd[128];

    auto g_x_z_xxy_yy = buffer_ppfd[129];

    auto g_x_z_xxy_yz = buffer_ppfd[130];

    auto g_x_z_xxy_zz = buffer_ppfd[131];

    auto g_x_z_xxz_xx = buffer_ppfd[132];

    auto g_x_z_xxz_xy = buffer_ppfd[133];

    auto g_x_z_xxz_xz = buffer_ppfd[134];

    auto g_x_z_xxz_yy = buffer_ppfd[135];

    auto g_x_z_xxz_yz = buffer_ppfd[136];

    auto g_x_z_xxz_zz = buffer_ppfd[137];

    auto g_x_z_xyy_xx = buffer_ppfd[138];

    auto g_x_z_xyy_xy = buffer_ppfd[139];

    auto g_x_z_xyy_xz = buffer_ppfd[140];

    auto g_x_z_xyy_yy = buffer_ppfd[141];

    auto g_x_z_xyy_yz = buffer_ppfd[142];

    auto g_x_z_xyy_zz = buffer_ppfd[143];

    auto g_x_z_xyz_xx = buffer_ppfd[144];

    auto g_x_z_xyz_xy = buffer_ppfd[145];

    auto g_x_z_xyz_xz = buffer_ppfd[146];

    auto g_x_z_xyz_yy = buffer_ppfd[147];

    auto g_x_z_xyz_yz = buffer_ppfd[148];

    auto g_x_z_xyz_zz = buffer_ppfd[149];

    auto g_x_z_xzz_xx = buffer_ppfd[150];

    auto g_x_z_xzz_xy = buffer_ppfd[151];

    auto g_x_z_xzz_xz = buffer_ppfd[152];

    auto g_x_z_xzz_yy = buffer_ppfd[153];

    auto g_x_z_xzz_yz = buffer_ppfd[154];

    auto g_x_z_xzz_zz = buffer_ppfd[155];

    auto g_x_z_yyy_xx = buffer_ppfd[156];

    auto g_x_z_yyy_xy = buffer_ppfd[157];

    auto g_x_z_yyy_xz = buffer_ppfd[158];

    auto g_x_z_yyy_yy = buffer_ppfd[159];

    auto g_x_z_yyy_yz = buffer_ppfd[160];

    auto g_x_z_yyy_zz = buffer_ppfd[161];

    auto g_x_z_yyz_xx = buffer_ppfd[162];

    auto g_x_z_yyz_xy = buffer_ppfd[163];

    auto g_x_z_yyz_xz = buffer_ppfd[164];

    auto g_x_z_yyz_yy = buffer_ppfd[165];

    auto g_x_z_yyz_yz = buffer_ppfd[166];

    auto g_x_z_yyz_zz = buffer_ppfd[167];

    auto g_x_z_yzz_xx = buffer_ppfd[168];

    auto g_x_z_yzz_xy = buffer_ppfd[169];

    auto g_x_z_yzz_xz = buffer_ppfd[170];

    auto g_x_z_yzz_yy = buffer_ppfd[171];

    auto g_x_z_yzz_yz = buffer_ppfd[172];

    auto g_x_z_yzz_zz = buffer_ppfd[173];

    auto g_x_z_zzz_xx = buffer_ppfd[174];

    auto g_x_z_zzz_xy = buffer_ppfd[175];

    auto g_x_z_zzz_xz = buffer_ppfd[176];

    auto g_x_z_zzz_yy = buffer_ppfd[177];

    auto g_x_z_zzz_yz = buffer_ppfd[178];

    auto g_x_z_zzz_zz = buffer_ppfd[179];

    auto g_y_x_xxx_xx = buffer_ppfd[180];

    auto g_y_x_xxx_xy = buffer_ppfd[181];

    auto g_y_x_xxx_xz = buffer_ppfd[182];

    auto g_y_x_xxx_yy = buffer_ppfd[183];

    auto g_y_x_xxx_yz = buffer_ppfd[184];

    auto g_y_x_xxx_zz = buffer_ppfd[185];

    auto g_y_x_xxy_xx = buffer_ppfd[186];

    auto g_y_x_xxy_xy = buffer_ppfd[187];

    auto g_y_x_xxy_xz = buffer_ppfd[188];

    auto g_y_x_xxy_yy = buffer_ppfd[189];

    auto g_y_x_xxy_yz = buffer_ppfd[190];

    auto g_y_x_xxy_zz = buffer_ppfd[191];

    auto g_y_x_xxz_xx = buffer_ppfd[192];

    auto g_y_x_xxz_xy = buffer_ppfd[193];

    auto g_y_x_xxz_xz = buffer_ppfd[194];

    auto g_y_x_xxz_yy = buffer_ppfd[195];

    auto g_y_x_xxz_yz = buffer_ppfd[196];

    auto g_y_x_xxz_zz = buffer_ppfd[197];

    auto g_y_x_xyy_xx = buffer_ppfd[198];

    auto g_y_x_xyy_xy = buffer_ppfd[199];

    auto g_y_x_xyy_xz = buffer_ppfd[200];

    auto g_y_x_xyy_yy = buffer_ppfd[201];

    auto g_y_x_xyy_yz = buffer_ppfd[202];

    auto g_y_x_xyy_zz = buffer_ppfd[203];

    auto g_y_x_xyz_xx = buffer_ppfd[204];

    auto g_y_x_xyz_xy = buffer_ppfd[205];

    auto g_y_x_xyz_xz = buffer_ppfd[206];

    auto g_y_x_xyz_yy = buffer_ppfd[207];

    auto g_y_x_xyz_yz = buffer_ppfd[208];

    auto g_y_x_xyz_zz = buffer_ppfd[209];

    auto g_y_x_xzz_xx = buffer_ppfd[210];

    auto g_y_x_xzz_xy = buffer_ppfd[211];

    auto g_y_x_xzz_xz = buffer_ppfd[212];

    auto g_y_x_xzz_yy = buffer_ppfd[213];

    auto g_y_x_xzz_yz = buffer_ppfd[214];

    auto g_y_x_xzz_zz = buffer_ppfd[215];

    auto g_y_x_yyy_xx = buffer_ppfd[216];

    auto g_y_x_yyy_xy = buffer_ppfd[217];

    auto g_y_x_yyy_xz = buffer_ppfd[218];

    auto g_y_x_yyy_yy = buffer_ppfd[219];

    auto g_y_x_yyy_yz = buffer_ppfd[220];

    auto g_y_x_yyy_zz = buffer_ppfd[221];

    auto g_y_x_yyz_xx = buffer_ppfd[222];

    auto g_y_x_yyz_xy = buffer_ppfd[223];

    auto g_y_x_yyz_xz = buffer_ppfd[224];

    auto g_y_x_yyz_yy = buffer_ppfd[225];

    auto g_y_x_yyz_yz = buffer_ppfd[226];

    auto g_y_x_yyz_zz = buffer_ppfd[227];

    auto g_y_x_yzz_xx = buffer_ppfd[228];

    auto g_y_x_yzz_xy = buffer_ppfd[229];

    auto g_y_x_yzz_xz = buffer_ppfd[230];

    auto g_y_x_yzz_yy = buffer_ppfd[231];

    auto g_y_x_yzz_yz = buffer_ppfd[232];

    auto g_y_x_yzz_zz = buffer_ppfd[233];

    auto g_y_x_zzz_xx = buffer_ppfd[234];

    auto g_y_x_zzz_xy = buffer_ppfd[235];

    auto g_y_x_zzz_xz = buffer_ppfd[236];

    auto g_y_x_zzz_yy = buffer_ppfd[237];

    auto g_y_x_zzz_yz = buffer_ppfd[238];

    auto g_y_x_zzz_zz = buffer_ppfd[239];

    auto g_y_y_xxx_xx = buffer_ppfd[240];

    auto g_y_y_xxx_xy = buffer_ppfd[241];

    auto g_y_y_xxx_xz = buffer_ppfd[242];

    auto g_y_y_xxx_yy = buffer_ppfd[243];

    auto g_y_y_xxx_yz = buffer_ppfd[244];

    auto g_y_y_xxx_zz = buffer_ppfd[245];

    auto g_y_y_xxy_xx = buffer_ppfd[246];

    auto g_y_y_xxy_xy = buffer_ppfd[247];

    auto g_y_y_xxy_xz = buffer_ppfd[248];

    auto g_y_y_xxy_yy = buffer_ppfd[249];

    auto g_y_y_xxy_yz = buffer_ppfd[250];

    auto g_y_y_xxy_zz = buffer_ppfd[251];

    auto g_y_y_xxz_xx = buffer_ppfd[252];

    auto g_y_y_xxz_xy = buffer_ppfd[253];

    auto g_y_y_xxz_xz = buffer_ppfd[254];

    auto g_y_y_xxz_yy = buffer_ppfd[255];

    auto g_y_y_xxz_yz = buffer_ppfd[256];

    auto g_y_y_xxz_zz = buffer_ppfd[257];

    auto g_y_y_xyy_xx = buffer_ppfd[258];

    auto g_y_y_xyy_xy = buffer_ppfd[259];

    auto g_y_y_xyy_xz = buffer_ppfd[260];

    auto g_y_y_xyy_yy = buffer_ppfd[261];

    auto g_y_y_xyy_yz = buffer_ppfd[262];

    auto g_y_y_xyy_zz = buffer_ppfd[263];

    auto g_y_y_xyz_xx = buffer_ppfd[264];

    auto g_y_y_xyz_xy = buffer_ppfd[265];

    auto g_y_y_xyz_xz = buffer_ppfd[266];

    auto g_y_y_xyz_yy = buffer_ppfd[267];

    auto g_y_y_xyz_yz = buffer_ppfd[268];

    auto g_y_y_xyz_zz = buffer_ppfd[269];

    auto g_y_y_xzz_xx = buffer_ppfd[270];

    auto g_y_y_xzz_xy = buffer_ppfd[271];

    auto g_y_y_xzz_xz = buffer_ppfd[272];

    auto g_y_y_xzz_yy = buffer_ppfd[273];

    auto g_y_y_xzz_yz = buffer_ppfd[274];

    auto g_y_y_xzz_zz = buffer_ppfd[275];

    auto g_y_y_yyy_xx = buffer_ppfd[276];

    auto g_y_y_yyy_xy = buffer_ppfd[277];

    auto g_y_y_yyy_xz = buffer_ppfd[278];

    auto g_y_y_yyy_yy = buffer_ppfd[279];

    auto g_y_y_yyy_yz = buffer_ppfd[280];

    auto g_y_y_yyy_zz = buffer_ppfd[281];

    auto g_y_y_yyz_xx = buffer_ppfd[282];

    auto g_y_y_yyz_xy = buffer_ppfd[283];

    auto g_y_y_yyz_xz = buffer_ppfd[284];

    auto g_y_y_yyz_yy = buffer_ppfd[285];

    auto g_y_y_yyz_yz = buffer_ppfd[286];

    auto g_y_y_yyz_zz = buffer_ppfd[287];

    auto g_y_y_yzz_xx = buffer_ppfd[288];

    auto g_y_y_yzz_xy = buffer_ppfd[289];

    auto g_y_y_yzz_xz = buffer_ppfd[290];

    auto g_y_y_yzz_yy = buffer_ppfd[291];

    auto g_y_y_yzz_yz = buffer_ppfd[292];

    auto g_y_y_yzz_zz = buffer_ppfd[293];

    auto g_y_y_zzz_xx = buffer_ppfd[294];

    auto g_y_y_zzz_xy = buffer_ppfd[295];

    auto g_y_y_zzz_xz = buffer_ppfd[296];

    auto g_y_y_zzz_yy = buffer_ppfd[297];

    auto g_y_y_zzz_yz = buffer_ppfd[298];

    auto g_y_y_zzz_zz = buffer_ppfd[299];

    auto g_y_z_xxx_xx = buffer_ppfd[300];

    auto g_y_z_xxx_xy = buffer_ppfd[301];

    auto g_y_z_xxx_xz = buffer_ppfd[302];

    auto g_y_z_xxx_yy = buffer_ppfd[303];

    auto g_y_z_xxx_yz = buffer_ppfd[304];

    auto g_y_z_xxx_zz = buffer_ppfd[305];

    auto g_y_z_xxy_xx = buffer_ppfd[306];

    auto g_y_z_xxy_xy = buffer_ppfd[307];

    auto g_y_z_xxy_xz = buffer_ppfd[308];

    auto g_y_z_xxy_yy = buffer_ppfd[309];

    auto g_y_z_xxy_yz = buffer_ppfd[310];

    auto g_y_z_xxy_zz = buffer_ppfd[311];

    auto g_y_z_xxz_xx = buffer_ppfd[312];

    auto g_y_z_xxz_xy = buffer_ppfd[313];

    auto g_y_z_xxz_xz = buffer_ppfd[314];

    auto g_y_z_xxz_yy = buffer_ppfd[315];

    auto g_y_z_xxz_yz = buffer_ppfd[316];

    auto g_y_z_xxz_zz = buffer_ppfd[317];

    auto g_y_z_xyy_xx = buffer_ppfd[318];

    auto g_y_z_xyy_xy = buffer_ppfd[319];

    auto g_y_z_xyy_xz = buffer_ppfd[320];

    auto g_y_z_xyy_yy = buffer_ppfd[321];

    auto g_y_z_xyy_yz = buffer_ppfd[322];

    auto g_y_z_xyy_zz = buffer_ppfd[323];

    auto g_y_z_xyz_xx = buffer_ppfd[324];

    auto g_y_z_xyz_xy = buffer_ppfd[325];

    auto g_y_z_xyz_xz = buffer_ppfd[326];

    auto g_y_z_xyz_yy = buffer_ppfd[327];

    auto g_y_z_xyz_yz = buffer_ppfd[328];

    auto g_y_z_xyz_zz = buffer_ppfd[329];

    auto g_y_z_xzz_xx = buffer_ppfd[330];

    auto g_y_z_xzz_xy = buffer_ppfd[331];

    auto g_y_z_xzz_xz = buffer_ppfd[332];

    auto g_y_z_xzz_yy = buffer_ppfd[333];

    auto g_y_z_xzz_yz = buffer_ppfd[334];

    auto g_y_z_xzz_zz = buffer_ppfd[335];

    auto g_y_z_yyy_xx = buffer_ppfd[336];

    auto g_y_z_yyy_xy = buffer_ppfd[337];

    auto g_y_z_yyy_xz = buffer_ppfd[338];

    auto g_y_z_yyy_yy = buffer_ppfd[339];

    auto g_y_z_yyy_yz = buffer_ppfd[340];

    auto g_y_z_yyy_zz = buffer_ppfd[341];

    auto g_y_z_yyz_xx = buffer_ppfd[342];

    auto g_y_z_yyz_xy = buffer_ppfd[343];

    auto g_y_z_yyz_xz = buffer_ppfd[344];

    auto g_y_z_yyz_yy = buffer_ppfd[345];

    auto g_y_z_yyz_yz = buffer_ppfd[346];

    auto g_y_z_yyz_zz = buffer_ppfd[347];

    auto g_y_z_yzz_xx = buffer_ppfd[348];

    auto g_y_z_yzz_xy = buffer_ppfd[349];

    auto g_y_z_yzz_xz = buffer_ppfd[350];

    auto g_y_z_yzz_yy = buffer_ppfd[351];

    auto g_y_z_yzz_yz = buffer_ppfd[352];

    auto g_y_z_yzz_zz = buffer_ppfd[353];

    auto g_y_z_zzz_xx = buffer_ppfd[354];

    auto g_y_z_zzz_xy = buffer_ppfd[355];

    auto g_y_z_zzz_xz = buffer_ppfd[356];

    auto g_y_z_zzz_yy = buffer_ppfd[357];

    auto g_y_z_zzz_yz = buffer_ppfd[358];

    auto g_y_z_zzz_zz = buffer_ppfd[359];

    auto g_z_x_xxx_xx = buffer_ppfd[360];

    auto g_z_x_xxx_xy = buffer_ppfd[361];

    auto g_z_x_xxx_xz = buffer_ppfd[362];

    auto g_z_x_xxx_yy = buffer_ppfd[363];

    auto g_z_x_xxx_yz = buffer_ppfd[364];

    auto g_z_x_xxx_zz = buffer_ppfd[365];

    auto g_z_x_xxy_xx = buffer_ppfd[366];

    auto g_z_x_xxy_xy = buffer_ppfd[367];

    auto g_z_x_xxy_xz = buffer_ppfd[368];

    auto g_z_x_xxy_yy = buffer_ppfd[369];

    auto g_z_x_xxy_yz = buffer_ppfd[370];

    auto g_z_x_xxy_zz = buffer_ppfd[371];

    auto g_z_x_xxz_xx = buffer_ppfd[372];

    auto g_z_x_xxz_xy = buffer_ppfd[373];

    auto g_z_x_xxz_xz = buffer_ppfd[374];

    auto g_z_x_xxz_yy = buffer_ppfd[375];

    auto g_z_x_xxz_yz = buffer_ppfd[376];

    auto g_z_x_xxz_zz = buffer_ppfd[377];

    auto g_z_x_xyy_xx = buffer_ppfd[378];

    auto g_z_x_xyy_xy = buffer_ppfd[379];

    auto g_z_x_xyy_xz = buffer_ppfd[380];

    auto g_z_x_xyy_yy = buffer_ppfd[381];

    auto g_z_x_xyy_yz = buffer_ppfd[382];

    auto g_z_x_xyy_zz = buffer_ppfd[383];

    auto g_z_x_xyz_xx = buffer_ppfd[384];

    auto g_z_x_xyz_xy = buffer_ppfd[385];

    auto g_z_x_xyz_xz = buffer_ppfd[386];

    auto g_z_x_xyz_yy = buffer_ppfd[387];

    auto g_z_x_xyz_yz = buffer_ppfd[388];

    auto g_z_x_xyz_zz = buffer_ppfd[389];

    auto g_z_x_xzz_xx = buffer_ppfd[390];

    auto g_z_x_xzz_xy = buffer_ppfd[391];

    auto g_z_x_xzz_xz = buffer_ppfd[392];

    auto g_z_x_xzz_yy = buffer_ppfd[393];

    auto g_z_x_xzz_yz = buffer_ppfd[394];

    auto g_z_x_xzz_zz = buffer_ppfd[395];

    auto g_z_x_yyy_xx = buffer_ppfd[396];

    auto g_z_x_yyy_xy = buffer_ppfd[397];

    auto g_z_x_yyy_xz = buffer_ppfd[398];

    auto g_z_x_yyy_yy = buffer_ppfd[399];

    auto g_z_x_yyy_yz = buffer_ppfd[400];

    auto g_z_x_yyy_zz = buffer_ppfd[401];

    auto g_z_x_yyz_xx = buffer_ppfd[402];

    auto g_z_x_yyz_xy = buffer_ppfd[403];

    auto g_z_x_yyz_xz = buffer_ppfd[404];

    auto g_z_x_yyz_yy = buffer_ppfd[405];

    auto g_z_x_yyz_yz = buffer_ppfd[406];

    auto g_z_x_yyz_zz = buffer_ppfd[407];

    auto g_z_x_yzz_xx = buffer_ppfd[408];

    auto g_z_x_yzz_xy = buffer_ppfd[409];

    auto g_z_x_yzz_xz = buffer_ppfd[410];

    auto g_z_x_yzz_yy = buffer_ppfd[411];

    auto g_z_x_yzz_yz = buffer_ppfd[412];

    auto g_z_x_yzz_zz = buffer_ppfd[413];

    auto g_z_x_zzz_xx = buffer_ppfd[414];

    auto g_z_x_zzz_xy = buffer_ppfd[415];

    auto g_z_x_zzz_xz = buffer_ppfd[416];

    auto g_z_x_zzz_yy = buffer_ppfd[417];

    auto g_z_x_zzz_yz = buffer_ppfd[418];

    auto g_z_x_zzz_zz = buffer_ppfd[419];

    auto g_z_y_xxx_xx = buffer_ppfd[420];

    auto g_z_y_xxx_xy = buffer_ppfd[421];

    auto g_z_y_xxx_xz = buffer_ppfd[422];

    auto g_z_y_xxx_yy = buffer_ppfd[423];

    auto g_z_y_xxx_yz = buffer_ppfd[424];

    auto g_z_y_xxx_zz = buffer_ppfd[425];

    auto g_z_y_xxy_xx = buffer_ppfd[426];

    auto g_z_y_xxy_xy = buffer_ppfd[427];

    auto g_z_y_xxy_xz = buffer_ppfd[428];

    auto g_z_y_xxy_yy = buffer_ppfd[429];

    auto g_z_y_xxy_yz = buffer_ppfd[430];

    auto g_z_y_xxy_zz = buffer_ppfd[431];

    auto g_z_y_xxz_xx = buffer_ppfd[432];

    auto g_z_y_xxz_xy = buffer_ppfd[433];

    auto g_z_y_xxz_xz = buffer_ppfd[434];

    auto g_z_y_xxz_yy = buffer_ppfd[435];

    auto g_z_y_xxz_yz = buffer_ppfd[436];

    auto g_z_y_xxz_zz = buffer_ppfd[437];

    auto g_z_y_xyy_xx = buffer_ppfd[438];

    auto g_z_y_xyy_xy = buffer_ppfd[439];

    auto g_z_y_xyy_xz = buffer_ppfd[440];

    auto g_z_y_xyy_yy = buffer_ppfd[441];

    auto g_z_y_xyy_yz = buffer_ppfd[442];

    auto g_z_y_xyy_zz = buffer_ppfd[443];

    auto g_z_y_xyz_xx = buffer_ppfd[444];

    auto g_z_y_xyz_xy = buffer_ppfd[445];

    auto g_z_y_xyz_xz = buffer_ppfd[446];

    auto g_z_y_xyz_yy = buffer_ppfd[447];

    auto g_z_y_xyz_yz = buffer_ppfd[448];

    auto g_z_y_xyz_zz = buffer_ppfd[449];

    auto g_z_y_xzz_xx = buffer_ppfd[450];

    auto g_z_y_xzz_xy = buffer_ppfd[451];

    auto g_z_y_xzz_xz = buffer_ppfd[452];

    auto g_z_y_xzz_yy = buffer_ppfd[453];

    auto g_z_y_xzz_yz = buffer_ppfd[454];

    auto g_z_y_xzz_zz = buffer_ppfd[455];

    auto g_z_y_yyy_xx = buffer_ppfd[456];

    auto g_z_y_yyy_xy = buffer_ppfd[457];

    auto g_z_y_yyy_xz = buffer_ppfd[458];

    auto g_z_y_yyy_yy = buffer_ppfd[459];

    auto g_z_y_yyy_yz = buffer_ppfd[460];

    auto g_z_y_yyy_zz = buffer_ppfd[461];

    auto g_z_y_yyz_xx = buffer_ppfd[462];

    auto g_z_y_yyz_xy = buffer_ppfd[463];

    auto g_z_y_yyz_xz = buffer_ppfd[464];

    auto g_z_y_yyz_yy = buffer_ppfd[465];

    auto g_z_y_yyz_yz = buffer_ppfd[466];

    auto g_z_y_yyz_zz = buffer_ppfd[467];

    auto g_z_y_yzz_xx = buffer_ppfd[468];

    auto g_z_y_yzz_xy = buffer_ppfd[469];

    auto g_z_y_yzz_xz = buffer_ppfd[470];

    auto g_z_y_yzz_yy = buffer_ppfd[471];

    auto g_z_y_yzz_yz = buffer_ppfd[472];

    auto g_z_y_yzz_zz = buffer_ppfd[473];

    auto g_z_y_zzz_xx = buffer_ppfd[474];

    auto g_z_y_zzz_xy = buffer_ppfd[475];

    auto g_z_y_zzz_xz = buffer_ppfd[476];

    auto g_z_y_zzz_yy = buffer_ppfd[477];

    auto g_z_y_zzz_yz = buffer_ppfd[478];

    auto g_z_y_zzz_zz = buffer_ppfd[479];

    auto g_z_z_xxx_xx = buffer_ppfd[480];

    auto g_z_z_xxx_xy = buffer_ppfd[481];

    auto g_z_z_xxx_xz = buffer_ppfd[482];

    auto g_z_z_xxx_yy = buffer_ppfd[483];

    auto g_z_z_xxx_yz = buffer_ppfd[484];

    auto g_z_z_xxx_zz = buffer_ppfd[485];

    auto g_z_z_xxy_xx = buffer_ppfd[486];

    auto g_z_z_xxy_xy = buffer_ppfd[487];

    auto g_z_z_xxy_xz = buffer_ppfd[488];

    auto g_z_z_xxy_yy = buffer_ppfd[489];

    auto g_z_z_xxy_yz = buffer_ppfd[490];

    auto g_z_z_xxy_zz = buffer_ppfd[491];

    auto g_z_z_xxz_xx = buffer_ppfd[492];

    auto g_z_z_xxz_xy = buffer_ppfd[493];

    auto g_z_z_xxz_xz = buffer_ppfd[494];

    auto g_z_z_xxz_yy = buffer_ppfd[495];

    auto g_z_z_xxz_yz = buffer_ppfd[496];

    auto g_z_z_xxz_zz = buffer_ppfd[497];

    auto g_z_z_xyy_xx = buffer_ppfd[498];

    auto g_z_z_xyy_xy = buffer_ppfd[499];

    auto g_z_z_xyy_xz = buffer_ppfd[500];

    auto g_z_z_xyy_yy = buffer_ppfd[501];

    auto g_z_z_xyy_yz = buffer_ppfd[502];

    auto g_z_z_xyy_zz = buffer_ppfd[503];

    auto g_z_z_xyz_xx = buffer_ppfd[504];

    auto g_z_z_xyz_xy = buffer_ppfd[505];

    auto g_z_z_xyz_xz = buffer_ppfd[506];

    auto g_z_z_xyz_yy = buffer_ppfd[507];

    auto g_z_z_xyz_yz = buffer_ppfd[508];

    auto g_z_z_xyz_zz = buffer_ppfd[509];

    auto g_z_z_xzz_xx = buffer_ppfd[510];

    auto g_z_z_xzz_xy = buffer_ppfd[511];

    auto g_z_z_xzz_xz = buffer_ppfd[512];

    auto g_z_z_xzz_yy = buffer_ppfd[513];

    auto g_z_z_xzz_yz = buffer_ppfd[514];

    auto g_z_z_xzz_zz = buffer_ppfd[515];

    auto g_z_z_yyy_xx = buffer_ppfd[516];

    auto g_z_z_yyy_xy = buffer_ppfd[517];

    auto g_z_z_yyy_xz = buffer_ppfd[518];

    auto g_z_z_yyy_yy = buffer_ppfd[519];

    auto g_z_z_yyy_yz = buffer_ppfd[520];

    auto g_z_z_yyy_zz = buffer_ppfd[521];

    auto g_z_z_yyz_xx = buffer_ppfd[522];

    auto g_z_z_yyz_xy = buffer_ppfd[523];

    auto g_z_z_yyz_xz = buffer_ppfd[524];

    auto g_z_z_yyz_yy = buffer_ppfd[525];

    auto g_z_z_yyz_yz = buffer_ppfd[526];

    auto g_z_z_yyz_zz = buffer_ppfd[527];

    auto g_z_z_yzz_xx = buffer_ppfd[528];

    auto g_z_z_yzz_xy = buffer_ppfd[529];

    auto g_z_z_yzz_xz = buffer_ppfd[530];

    auto g_z_z_yzz_yy = buffer_ppfd[531];

    auto g_z_z_yzz_yz = buffer_ppfd[532];

    auto g_z_z_yzz_zz = buffer_ppfd[533];

    auto g_z_z_zzz_xx = buffer_ppfd[534];

    auto g_z_z_zzz_xy = buffer_ppfd[535];

    auto g_z_z_zzz_xz = buffer_ppfd[536];

    auto g_z_z_zzz_yy = buffer_ppfd[537];

    auto g_z_z_zzz_yz = buffer_ppfd[538];

    auto g_z_z_zzz_zz = buffer_ppfd[539];

    /// Set up components of integrals buffer : buffer_1010_spdd

    auto g_x_0_x_0_0_x_xx_xx = buffer_1010_spdd[0];

    auto g_x_0_x_0_0_x_xx_xy = buffer_1010_spdd[1];

    auto g_x_0_x_0_0_x_xx_xz = buffer_1010_spdd[2];

    auto g_x_0_x_0_0_x_xx_yy = buffer_1010_spdd[3];

    auto g_x_0_x_0_0_x_xx_yz = buffer_1010_spdd[4];

    auto g_x_0_x_0_0_x_xx_zz = buffer_1010_spdd[5];

    auto g_x_0_x_0_0_x_xy_xx = buffer_1010_spdd[6];

    auto g_x_0_x_0_0_x_xy_xy = buffer_1010_spdd[7];

    auto g_x_0_x_0_0_x_xy_xz = buffer_1010_spdd[8];

    auto g_x_0_x_0_0_x_xy_yy = buffer_1010_spdd[9];

    auto g_x_0_x_0_0_x_xy_yz = buffer_1010_spdd[10];

    auto g_x_0_x_0_0_x_xy_zz = buffer_1010_spdd[11];

    auto g_x_0_x_0_0_x_xz_xx = buffer_1010_spdd[12];

    auto g_x_0_x_0_0_x_xz_xy = buffer_1010_spdd[13];

    auto g_x_0_x_0_0_x_xz_xz = buffer_1010_spdd[14];

    auto g_x_0_x_0_0_x_xz_yy = buffer_1010_spdd[15];

    auto g_x_0_x_0_0_x_xz_yz = buffer_1010_spdd[16];

    auto g_x_0_x_0_0_x_xz_zz = buffer_1010_spdd[17];

    auto g_x_0_x_0_0_x_yy_xx = buffer_1010_spdd[18];

    auto g_x_0_x_0_0_x_yy_xy = buffer_1010_spdd[19];

    auto g_x_0_x_0_0_x_yy_xz = buffer_1010_spdd[20];

    auto g_x_0_x_0_0_x_yy_yy = buffer_1010_spdd[21];

    auto g_x_0_x_0_0_x_yy_yz = buffer_1010_spdd[22];

    auto g_x_0_x_0_0_x_yy_zz = buffer_1010_spdd[23];

    auto g_x_0_x_0_0_x_yz_xx = buffer_1010_spdd[24];

    auto g_x_0_x_0_0_x_yz_xy = buffer_1010_spdd[25];

    auto g_x_0_x_0_0_x_yz_xz = buffer_1010_spdd[26];

    auto g_x_0_x_0_0_x_yz_yy = buffer_1010_spdd[27];

    auto g_x_0_x_0_0_x_yz_yz = buffer_1010_spdd[28];

    auto g_x_0_x_0_0_x_yz_zz = buffer_1010_spdd[29];

    auto g_x_0_x_0_0_x_zz_xx = buffer_1010_spdd[30];

    auto g_x_0_x_0_0_x_zz_xy = buffer_1010_spdd[31];

    auto g_x_0_x_0_0_x_zz_xz = buffer_1010_spdd[32];

    auto g_x_0_x_0_0_x_zz_yy = buffer_1010_spdd[33];

    auto g_x_0_x_0_0_x_zz_yz = buffer_1010_spdd[34];

    auto g_x_0_x_0_0_x_zz_zz = buffer_1010_spdd[35];

    auto g_x_0_x_0_0_y_xx_xx = buffer_1010_spdd[36];

    auto g_x_0_x_0_0_y_xx_xy = buffer_1010_spdd[37];

    auto g_x_0_x_0_0_y_xx_xz = buffer_1010_spdd[38];

    auto g_x_0_x_0_0_y_xx_yy = buffer_1010_spdd[39];

    auto g_x_0_x_0_0_y_xx_yz = buffer_1010_spdd[40];

    auto g_x_0_x_0_0_y_xx_zz = buffer_1010_spdd[41];

    auto g_x_0_x_0_0_y_xy_xx = buffer_1010_spdd[42];

    auto g_x_0_x_0_0_y_xy_xy = buffer_1010_spdd[43];

    auto g_x_0_x_0_0_y_xy_xz = buffer_1010_spdd[44];

    auto g_x_0_x_0_0_y_xy_yy = buffer_1010_spdd[45];

    auto g_x_0_x_0_0_y_xy_yz = buffer_1010_spdd[46];

    auto g_x_0_x_0_0_y_xy_zz = buffer_1010_spdd[47];

    auto g_x_0_x_0_0_y_xz_xx = buffer_1010_spdd[48];

    auto g_x_0_x_0_0_y_xz_xy = buffer_1010_spdd[49];

    auto g_x_0_x_0_0_y_xz_xz = buffer_1010_spdd[50];

    auto g_x_0_x_0_0_y_xz_yy = buffer_1010_spdd[51];

    auto g_x_0_x_0_0_y_xz_yz = buffer_1010_spdd[52];

    auto g_x_0_x_0_0_y_xz_zz = buffer_1010_spdd[53];

    auto g_x_0_x_0_0_y_yy_xx = buffer_1010_spdd[54];

    auto g_x_0_x_0_0_y_yy_xy = buffer_1010_spdd[55];

    auto g_x_0_x_0_0_y_yy_xz = buffer_1010_spdd[56];

    auto g_x_0_x_0_0_y_yy_yy = buffer_1010_spdd[57];

    auto g_x_0_x_0_0_y_yy_yz = buffer_1010_spdd[58];

    auto g_x_0_x_0_0_y_yy_zz = buffer_1010_spdd[59];

    auto g_x_0_x_0_0_y_yz_xx = buffer_1010_spdd[60];

    auto g_x_0_x_0_0_y_yz_xy = buffer_1010_spdd[61];

    auto g_x_0_x_0_0_y_yz_xz = buffer_1010_spdd[62];

    auto g_x_0_x_0_0_y_yz_yy = buffer_1010_spdd[63];

    auto g_x_0_x_0_0_y_yz_yz = buffer_1010_spdd[64];

    auto g_x_0_x_0_0_y_yz_zz = buffer_1010_spdd[65];

    auto g_x_0_x_0_0_y_zz_xx = buffer_1010_spdd[66];

    auto g_x_0_x_0_0_y_zz_xy = buffer_1010_spdd[67];

    auto g_x_0_x_0_0_y_zz_xz = buffer_1010_spdd[68];

    auto g_x_0_x_0_0_y_zz_yy = buffer_1010_spdd[69];

    auto g_x_0_x_0_0_y_zz_yz = buffer_1010_spdd[70];

    auto g_x_0_x_0_0_y_zz_zz = buffer_1010_spdd[71];

    auto g_x_0_x_0_0_z_xx_xx = buffer_1010_spdd[72];

    auto g_x_0_x_0_0_z_xx_xy = buffer_1010_spdd[73];

    auto g_x_0_x_0_0_z_xx_xz = buffer_1010_spdd[74];

    auto g_x_0_x_0_0_z_xx_yy = buffer_1010_spdd[75];

    auto g_x_0_x_0_0_z_xx_yz = buffer_1010_spdd[76];

    auto g_x_0_x_0_0_z_xx_zz = buffer_1010_spdd[77];

    auto g_x_0_x_0_0_z_xy_xx = buffer_1010_spdd[78];

    auto g_x_0_x_0_0_z_xy_xy = buffer_1010_spdd[79];

    auto g_x_0_x_0_0_z_xy_xz = buffer_1010_spdd[80];

    auto g_x_0_x_0_0_z_xy_yy = buffer_1010_spdd[81];

    auto g_x_0_x_0_0_z_xy_yz = buffer_1010_spdd[82];

    auto g_x_0_x_0_0_z_xy_zz = buffer_1010_spdd[83];

    auto g_x_0_x_0_0_z_xz_xx = buffer_1010_spdd[84];

    auto g_x_0_x_0_0_z_xz_xy = buffer_1010_spdd[85];

    auto g_x_0_x_0_0_z_xz_xz = buffer_1010_spdd[86];

    auto g_x_0_x_0_0_z_xz_yy = buffer_1010_spdd[87];

    auto g_x_0_x_0_0_z_xz_yz = buffer_1010_spdd[88];

    auto g_x_0_x_0_0_z_xz_zz = buffer_1010_spdd[89];

    auto g_x_0_x_0_0_z_yy_xx = buffer_1010_spdd[90];

    auto g_x_0_x_0_0_z_yy_xy = buffer_1010_spdd[91];

    auto g_x_0_x_0_0_z_yy_xz = buffer_1010_spdd[92];

    auto g_x_0_x_0_0_z_yy_yy = buffer_1010_spdd[93];

    auto g_x_0_x_0_0_z_yy_yz = buffer_1010_spdd[94];

    auto g_x_0_x_0_0_z_yy_zz = buffer_1010_spdd[95];

    auto g_x_0_x_0_0_z_yz_xx = buffer_1010_spdd[96];

    auto g_x_0_x_0_0_z_yz_xy = buffer_1010_spdd[97];

    auto g_x_0_x_0_0_z_yz_xz = buffer_1010_spdd[98];

    auto g_x_0_x_0_0_z_yz_yy = buffer_1010_spdd[99];

    auto g_x_0_x_0_0_z_yz_yz = buffer_1010_spdd[100];

    auto g_x_0_x_0_0_z_yz_zz = buffer_1010_spdd[101];

    auto g_x_0_x_0_0_z_zz_xx = buffer_1010_spdd[102];

    auto g_x_0_x_0_0_z_zz_xy = buffer_1010_spdd[103];

    auto g_x_0_x_0_0_z_zz_xz = buffer_1010_spdd[104];

    auto g_x_0_x_0_0_z_zz_yy = buffer_1010_spdd[105];

    auto g_x_0_x_0_0_z_zz_yz = buffer_1010_spdd[106];

    auto g_x_0_x_0_0_z_zz_zz = buffer_1010_spdd[107];

    auto g_x_0_y_0_0_x_xx_xx = buffer_1010_spdd[108];

    auto g_x_0_y_0_0_x_xx_xy = buffer_1010_spdd[109];

    auto g_x_0_y_0_0_x_xx_xz = buffer_1010_spdd[110];

    auto g_x_0_y_0_0_x_xx_yy = buffer_1010_spdd[111];

    auto g_x_0_y_0_0_x_xx_yz = buffer_1010_spdd[112];

    auto g_x_0_y_0_0_x_xx_zz = buffer_1010_spdd[113];

    auto g_x_0_y_0_0_x_xy_xx = buffer_1010_spdd[114];

    auto g_x_0_y_0_0_x_xy_xy = buffer_1010_spdd[115];

    auto g_x_0_y_0_0_x_xy_xz = buffer_1010_spdd[116];

    auto g_x_0_y_0_0_x_xy_yy = buffer_1010_spdd[117];

    auto g_x_0_y_0_0_x_xy_yz = buffer_1010_spdd[118];

    auto g_x_0_y_0_0_x_xy_zz = buffer_1010_spdd[119];

    auto g_x_0_y_0_0_x_xz_xx = buffer_1010_spdd[120];

    auto g_x_0_y_0_0_x_xz_xy = buffer_1010_spdd[121];

    auto g_x_0_y_0_0_x_xz_xz = buffer_1010_spdd[122];

    auto g_x_0_y_0_0_x_xz_yy = buffer_1010_spdd[123];

    auto g_x_0_y_0_0_x_xz_yz = buffer_1010_spdd[124];

    auto g_x_0_y_0_0_x_xz_zz = buffer_1010_spdd[125];

    auto g_x_0_y_0_0_x_yy_xx = buffer_1010_spdd[126];

    auto g_x_0_y_0_0_x_yy_xy = buffer_1010_spdd[127];

    auto g_x_0_y_0_0_x_yy_xz = buffer_1010_spdd[128];

    auto g_x_0_y_0_0_x_yy_yy = buffer_1010_spdd[129];

    auto g_x_0_y_0_0_x_yy_yz = buffer_1010_spdd[130];

    auto g_x_0_y_0_0_x_yy_zz = buffer_1010_spdd[131];

    auto g_x_0_y_0_0_x_yz_xx = buffer_1010_spdd[132];

    auto g_x_0_y_0_0_x_yz_xy = buffer_1010_spdd[133];

    auto g_x_0_y_0_0_x_yz_xz = buffer_1010_spdd[134];

    auto g_x_0_y_0_0_x_yz_yy = buffer_1010_spdd[135];

    auto g_x_0_y_0_0_x_yz_yz = buffer_1010_spdd[136];

    auto g_x_0_y_0_0_x_yz_zz = buffer_1010_spdd[137];

    auto g_x_0_y_0_0_x_zz_xx = buffer_1010_spdd[138];

    auto g_x_0_y_0_0_x_zz_xy = buffer_1010_spdd[139];

    auto g_x_0_y_0_0_x_zz_xz = buffer_1010_spdd[140];

    auto g_x_0_y_0_0_x_zz_yy = buffer_1010_spdd[141];

    auto g_x_0_y_0_0_x_zz_yz = buffer_1010_spdd[142];

    auto g_x_0_y_0_0_x_zz_zz = buffer_1010_spdd[143];

    auto g_x_0_y_0_0_y_xx_xx = buffer_1010_spdd[144];

    auto g_x_0_y_0_0_y_xx_xy = buffer_1010_spdd[145];

    auto g_x_0_y_0_0_y_xx_xz = buffer_1010_spdd[146];

    auto g_x_0_y_0_0_y_xx_yy = buffer_1010_spdd[147];

    auto g_x_0_y_0_0_y_xx_yz = buffer_1010_spdd[148];

    auto g_x_0_y_0_0_y_xx_zz = buffer_1010_spdd[149];

    auto g_x_0_y_0_0_y_xy_xx = buffer_1010_spdd[150];

    auto g_x_0_y_0_0_y_xy_xy = buffer_1010_spdd[151];

    auto g_x_0_y_0_0_y_xy_xz = buffer_1010_spdd[152];

    auto g_x_0_y_0_0_y_xy_yy = buffer_1010_spdd[153];

    auto g_x_0_y_0_0_y_xy_yz = buffer_1010_spdd[154];

    auto g_x_0_y_0_0_y_xy_zz = buffer_1010_spdd[155];

    auto g_x_0_y_0_0_y_xz_xx = buffer_1010_spdd[156];

    auto g_x_0_y_0_0_y_xz_xy = buffer_1010_spdd[157];

    auto g_x_0_y_0_0_y_xz_xz = buffer_1010_spdd[158];

    auto g_x_0_y_0_0_y_xz_yy = buffer_1010_spdd[159];

    auto g_x_0_y_0_0_y_xz_yz = buffer_1010_spdd[160];

    auto g_x_0_y_0_0_y_xz_zz = buffer_1010_spdd[161];

    auto g_x_0_y_0_0_y_yy_xx = buffer_1010_spdd[162];

    auto g_x_0_y_0_0_y_yy_xy = buffer_1010_spdd[163];

    auto g_x_0_y_0_0_y_yy_xz = buffer_1010_spdd[164];

    auto g_x_0_y_0_0_y_yy_yy = buffer_1010_spdd[165];

    auto g_x_0_y_0_0_y_yy_yz = buffer_1010_spdd[166];

    auto g_x_0_y_0_0_y_yy_zz = buffer_1010_spdd[167];

    auto g_x_0_y_0_0_y_yz_xx = buffer_1010_spdd[168];

    auto g_x_0_y_0_0_y_yz_xy = buffer_1010_spdd[169];

    auto g_x_0_y_0_0_y_yz_xz = buffer_1010_spdd[170];

    auto g_x_0_y_0_0_y_yz_yy = buffer_1010_spdd[171];

    auto g_x_0_y_0_0_y_yz_yz = buffer_1010_spdd[172];

    auto g_x_0_y_0_0_y_yz_zz = buffer_1010_spdd[173];

    auto g_x_0_y_0_0_y_zz_xx = buffer_1010_spdd[174];

    auto g_x_0_y_0_0_y_zz_xy = buffer_1010_spdd[175];

    auto g_x_0_y_0_0_y_zz_xz = buffer_1010_spdd[176];

    auto g_x_0_y_0_0_y_zz_yy = buffer_1010_spdd[177];

    auto g_x_0_y_0_0_y_zz_yz = buffer_1010_spdd[178];

    auto g_x_0_y_0_0_y_zz_zz = buffer_1010_spdd[179];

    auto g_x_0_y_0_0_z_xx_xx = buffer_1010_spdd[180];

    auto g_x_0_y_0_0_z_xx_xy = buffer_1010_spdd[181];

    auto g_x_0_y_0_0_z_xx_xz = buffer_1010_spdd[182];

    auto g_x_0_y_0_0_z_xx_yy = buffer_1010_spdd[183];

    auto g_x_0_y_0_0_z_xx_yz = buffer_1010_spdd[184];

    auto g_x_0_y_0_0_z_xx_zz = buffer_1010_spdd[185];

    auto g_x_0_y_0_0_z_xy_xx = buffer_1010_spdd[186];

    auto g_x_0_y_0_0_z_xy_xy = buffer_1010_spdd[187];

    auto g_x_0_y_0_0_z_xy_xz = buffer_1010_spdd[188];

    auto g_x_0_y_0_0_z_xy_yy = buffer_1010_spdd[189];

    auto g_x_0_y_0_0_z_xy_yz = buffer_1010_spdd[190];

    auto g_x_0_y_0_0_z_xy_zz = buffer_1010_spdd[191];

    auto g_x_0_y_0_0_z_xz_xx = buffer_1010_spdd[192];

    auto g_x_0_y_0_0_z_xz_xy = buffer_1010_spdd[193];

    auto g_x_0_y_0_0_z_xz_xz = buffer_1010_spdd[194];

    auto g_x_0_y_0_0_z_xz_yy = buffer_1010_spdd[195];

    auto g_x_0_y_0_0_z_xz_yz = buffer_1010_spdd[196];

    auto g_x_0_y_0_0_z_xz_zz = buffer_1010_spdd[197];

    auto g_x_0_y_0_0_z_yy_xx = buffer_1010_spdd[198];

    auto g_x_0_y_0_0_z_yy_xy = buffer_1010_spdd[199];

    auto g_x_0_y_0_0_z_yy_xz = buffer_1010_spdd[200];

    auto g_x_0_y_0_0_z_yy_yy = buffer_1010_spdd[201];

    auto g_x_0_y_0_0_z_yy_yz = buffer_1010_spdd[202];

    auto g_x_0_y_0_0_z_yy_zz = buffer_1010_spdd[203];

    auto g_x_0_y_0_0_z_yz_xx = buffer_1010_spdd[204];

    auto g_x_0_y_0_0_z_yz_xy = buffer_1010_spdd[205];

    auto g_x_0_y_0_0_z_yz_xz = buffer_1010_spdd[206];

    auto g_x_0_y_0_0_z_yz_yy = buffer_1010_spdd[207];

    auto g_x_0_y_0_0_z_yz_yz = buffer_1010_spdd[208];

    auto g_x_0_y_0_0_z_yz_zz = buffer_1010_spdd[209];

    auto g_x_0_y_0_0_z_zz_xx = buffer_1010_spdd[210];

    auto g_x_0_y_0_0_z_zz_xy = buffer_1010_spdd[211];

    auto g_x_0_y_0_0_z_zz_xz = buffer_1010_spdd[212];

    auto g_x_0_y_0_0_z_zz_yy = buffer_1010_spdd[213];

    auto g_x_0_y_0_0_z_zz_yz = buffer_1010_spdd[214];

    auto g_x_0_y_0_0_z_zz_zz = buffer_1010_spdd[215];

    auto g_x_0_z_0_0_x_xx_xx = buffer_1010_spdd[216];

    auto g_x_0_z_0_0_x_xx_xy = buffer_1010_spdd[217];

    auto g_x_0_z_0_0_x_xx_xz = buffer_1010_spdd[218];

    auto g_x_0_z_0_0_x_xx_yy = buffer_1010_spdd[219];

    auto g_x_0_z_0_0_x_xx_yz = buffer_1010_spdd[220];

    auto g_x_0_z_0_0_x_xx_zz = buffer_1010_spdd[221];

    auto g_x_0_z_0_0_x_xy_xx = buffer_1010_spdd[222];

    auto g_x_0_z_0_0_x_xy_xy = buffer_1010_spdd[223];

    auto g_x_0_z_0_0_x_xy_xz = buffer_1010_spdd[224];

    auto g_x_0_z_0_0_x_xy_yy = buffer_1010_spdd[225];

    auto g_x_0_z_0_0_x_xy_yz = buffer_1010_spdd[226];

    auto g_x_0_z_0_0_x_xy_zz = buffer_1010_spdd[227];

    auto g_x_0_z_0_0_x_xz_xx = buffer_1010_spdd[228];

    auto g_x_0_z_0_0_x_xz_xy = buffer_1010_spdd[229];

    auto g_x_0_z_0_0_x_xz_xz = buffer_1010_spdd[230];

    auto g_x_0_z_0_0_x_xz_yy = buffer_1010_spdd[231];

    auto g_x_0_z_0_0_x_xz_yz = buffer_1010_spdd[232];

    auto g_x_0_z_0_0_x_xz_zz = buffer_1010_spdd[233];

    auto g_x_0_z_0_0_x_yy_xx = buffer_1010_spdd[234];

    auto g_x_0_z_0_0_x_yy_xy = buffer_1010_spdd[235];

    auto g_x_0_z_0_0_x_yy_xz = buffer_1010_spdd[236];

    auto g_x_0_z_0_0_x_yy_yy = buffer_1010_spdd[237];

    auto g_x_0_z_0_0_x_yy_yz = buffer_1010_spdd[238];

    auto g_x_0_z_0_0_x_yy_zz = buffer_1010_spdd[239];

    auto g_x_0_z_0_0_x_yz_xx = buffer_1010_spdd[240];

    auto g_x_0_z_0_0_x_yz_xy = buffer_1010_spdd[241];

    auto g_x_0_z_0_0_x_yz_xz = buffer_1010_spdd[242];

    auto g_x_0_z_0_0_x_yz_yy = buffer_1010_spdd[243];

    auto g_x_0_z_0_0_x_yz_yz = buffer_1010_spdd[244];

    auto g_x_0_z_0_0_x_yz_zz = buffer_1010_spdd[245];

    auto g_x_0_z_0_0_x_zz_xx = buffer_1010_spdd[246];

    auto g_x_0_z_0_0_x_zz_xy = buffer_1010_spdd[247];

    auto g_x_0_z_0_0_x_zz_xz = buffer_1010_spdd[248];

    auto g_x_0_z_0_0_x_zz_yy = buffer_1010_spdd[249];

    auto g_x_0_z_0_0_x_zz_yz = buffer_1010_spdd[250];

    auto g_x_0_z_0_0_x_zz_zz = buffer_1010_spdd[251];

    auto g_x_0_z_0_0_y_xx_xx = buffer_1010_spdd[252];

    auto g_x_0_z_0_0_y_xx_xy = buffer_1010_spdd[253];

    auto g_x_0_z_0_0_y_xx_xz = buffer_1010_spdd[254];

    auto g_x_0_z_0_0_y_xx_yy = buffer_1010_spdd[255];

    auto g_x_0_z_0_0_y_xx_yz = buffer_1010_spdd[256];

    auto g_x_0_z_0_0_y_xx_zz = buffer_1010_spdd[257];

    auto g_x_0_z_0_0_y_xy_xx = buffer_1010_spdd[258];

    auto g_x_0_z_0_0_y_xy_xy = buffer_1010_spdd[259];

    auto g_x_0_z_0_0_y_xy_xz = buffer_1010_spdd[260];

    auto g_x_0_z_0_0_y_xy_yy = buffer_1010_spdd[261];

    auto g_x_0_z_0_0_y_xy_yz = buffer_1010_spdd[262];

    auto g_x_0_z_0_0_y_xy_zz = buffer_1010_spdd[263];

    auto g_x_0_z_0_0_y_xz_xx = buffer_1010_spdd[264];

    auto g_x_0_z_0_0_y_xz_xy = buffer_1010_spdd[265];

    auto g_x_0_z_0_0_y_xz_xz = buffer_1010_spdd[266];

    auto g_x_0_z_0_0_y_xz_yy = buffer_1010_spdd[267];

    auto g_x_0_z_0_0_y_xz_yz = buffer_1010_spdd[268];

    auto g_x_0_z_0_0_y_xz_zz = buffer_1010_spdd[269];

    auto g_x_0_z_0_0_y_yy_xx = buffer_1010_spdd[270];

    auto g_x_0_z_0_0_y_yy_xy = buffer_1010_spdd[271];

    auto g_x_0_z_0_0_y_yy_xz = buffer_1010_spdd[272];

    auto g_x_0_z_0_0_y_yy_yy = buffer_1010_spdd[273];

    auto g_x_0_z_0_0_y_yy_yz = buffer_1010_spdd[274];

    auto g_x_0_z_0_0_y_yy_zz = buffer_1010_spdd[275];

    auto g_x_0_z_0_0_y_yz_xx = buffer_1010_spdd[276];

    auto g_x_0_z_0_0_y_yz_xy = buffer_1010_spdd[277];

    auto g_x_0_z_0_0_y_yz_xz = buffer_1010_spdd[278];

    auto g_x_0_z_0_0_y_yz_yy = buffer_1010_spdd[279];

    auto g_x_0_z_0_0_y_yz_yz = buffer_1010_spdd[280];

    auto g_x_0_z_0_0_y_yz_zz = buffer_1010_spdd[281];

    auto g_x_0_z_0_0_y_zz_xx = buffer_1010_spdd[282];

    auto g_x_0_z_0_0_y_zz_xy = buffer_1010_spdd[283];

    auto g_x_0_z_0_0_y_zz_xz = buffer_1010_spdd[284];

    auto g_x_0_z_0_0_y_zz_yy = buffer_1010_spdd[285];

    auto g_x_0_z_0_0_y_zz_yz = buffer_1010_spdd[286];

    auto g_x_0_z_0_0_y_zz_zz = buffer_1010_spdd[287];

    auto g_x_0_z_0_0_z_xx_xx = buffer_1010_spdd[288];

    auto g_x_0_z_0_0_z_xx_xy = buffer_1010_spdd[289];

    auto g_x_0_z_0_0_z_xx_xz = buffer_1010_spdd[290];

    auto g_x_0_z_0_0_z_xx_yy = buffer_1010_spdd[291];

    auto g_x_0_z_0_0_z_xx_yz = buffer_1010_spdd[292];

    auto g_x_0_z_0_0_z_xx_zz = buffer_1010_spdd[293];

    auto g_x_0_z_0_0_z_xy_xx = buffer_1010_spdd[294];

    auto g_x_0_z_0_0_z_xy_xy = buffer_1010_spdd[295];

    auto g_x_0_z_0_0_z_xy_xz = buffer_1010_spdd[296];

    auto g_x_0_z_0_0_z_xy_yy = buffer_1010_spdd[297];

    auto g_x_0_z_0_0_z_xy_yz = buffer_1010_spdd[298];

    auto g_x_0_z_0_0_z_xy_zz = buffer_1010_spdd[299];

    auto g_x_0_z_0_0_z_xz_xx = buffer_1010_spdd[300];

    auto g_x_0_z_0_0_z_xz_xy = buffer_1010_spdd[301];

    auto g_x_0_z_0_0_z_xz_xz = buffer_1010_spdd[302];

    auto g_x_0_z_0_0_z_xz_yy = buffer_1010_spdd[303];

    auto g_x_0_z_0_0_z_xz_yz = buffer_1010_spdd[304];

    auto g_x_0_z_0_0_z_xz_zz = buffer_1010_spdd[305];

    auto g_x_0_z_0_0_z_yy_xx = buffer_1010_spdd[306];

    auto g_x_0_z_0_0_z_yy_xy = buffer_1010_spdd[307];

    auto g_x_0_z_0_0_z_yy_xz = buffer_1010_spdd[308];

    auto g_x_0_z_0_0_z_yy_yy = buffer_1010_spdd[309];

    auto g_x_0_z_0_0_z_yy_yz = buffer_1010_spdd[310];

    auto g_x_0_z_0_0_z_yy_zz = buffer_1010_spdd[311];

    auto g_x_0_z_0_0_z_yz_xx = buffer_1010_spdd[312];

    auto g_x_0_z_0_0_z_yz_xy = buffer_1010_spdd[313];

    auto g_x_0_z_0_0_z_yz_xz = buffer_1010_spdd[314];

    auto g_x_0_z_0_0_z_yz_yy = buffer_1010_spdd[315];

    auto g_x_0_z_0_0_z_yz_yz = buffer_1010_spdd[316];

    auto g_x_0_z_0_0_z_yz_zz = buffer_1010_spdd[317];

    auto g_x_0_z_0_0_z_zz_xx = buffer_1010_spdd[318];

    auto g_x_0_z_0_0_z_zz_xy = buffer_1010_spdd[319];

    auto g_x_0_z_0_0_z_zz_xz = buffer_1010_spdd[320];

    auto g_x_0_z_0_0_z_zz_yy = buffer_1010_spdd[321];

    auto g_x_0_z_0_0_z_zz_yz = buffer_1010_spdd[322];

    auto g_x_0_z_0_0_z_zz_zz = buffer_1010_spdd[323];

    auto g_y_0_x_0_0_x_xx_xx = buffer_1010_spdd[324];

    auto g_y_0_x_0_0_x_xx_xy = buffer_1010_spdd[325];

    auto g_y_0_x_0_0_x_xx_xz = buffer_1010_spdd[326];

    auto g_y_0_x_0_0_x_xx_yy = buffer_1010_spdd[327];

    auto g_y_0_x_0_0_x_xx_yz = buffer_1010_spdd[328];

    auto g_y_0_x_0_0_x_xx_zz = buffer_1010_spdd[329];

    auto g_y_0_x_0_0_x_xy_xx = buffer_1010_spdd[330];

    auto g_y_0_x_0_0_x_xy_xy = buffer_1010_spdd[331];

    auto g_y_0_x_0_0_x_xy_xz = buffer_1010_spdd[332];

    auto g_y_0_x_0_0_x_xy_yy = buffer_1010_spdd[333];

    auto g_y_0_x_0_0_x_xy_yz = buffer_1010_spdd[334];

    auto g_y_0_x_0_0_x_xy_zz = buffer_1010_spdd[335];

    auto g_y_0_x_0_0_x_xz_xx = buffer_1010_spdd[336];

    auto g_y_0_x_0_0_x_xz_xy = buffer_1010_spdd[337];

    auto g_y_0_x_0_0_x_xz_xz = buffer_1010_spdd[338];

    auto g_y_0_x_0_0_x_xz_yy = buffer_1010_spdd[339];

    auto g_y_0_x_0_0_x_xz_yz = buffer_1010_spdd[340];

    auto g_y_0_x_0_0_x_xz_zz = buffer_1010_spdd[341];

    auto g_y_0_x_0_0_x_yy_xx = buffer_1010_spdd[342];

    auto g_y_0_x_0_0_x_yy_xy = buffer_1010_spdd[343];

    auto g_y_0_x_0_0_x_yy_xz = buffer_1010_spdd[344];

    auto g_y_0_x_0_0_x_yy_yy = buffer_1010_spdd[345];

    auto g_y_0_x_0_0_x_yy_yz = buffer_1010_spdd[346];

    auto g_y_0_x_0_0_x_yy_zz = buffer_1010_spdd[347];

    auto g_y_0_x_0_0_x_yz_xx = buffer_1010_spdd[348];

    auto g_y_0_x_0_0_x_yz_xy = buffer_1010_spdd[349];

    auto g_y_0_x_0_0_x_yz_xz = buffer_1010_spdd[350];

    auto g_y_0_x_0_0_x_yz_yy = buffer_1010_spdd[351];

    auto g_y_0_x_0_0_x_yz_yz = buffer_1010_spdd[352];

    auto g_y_0_x_0_0_x_yz_zz = buffer_1010_spdd[353];

    auto g_y_0_x_0_0_x_zz_xx = buffer_1010_spdd[354];

    auto g_y_0_x_0_0_x_zz_xy = buffer_1010_spdd[355];

    auto g_y_0_x_0_0_x_zz_xz = buffer_1010_spdd[356];

    auto g_y_0_x_0_0_x_zz_yy = buffer_1010_spdd[357];

    auto g_y_0_x_0_0_x_zz_yz = buffer_1010_spdd[358];

    auto g_y_0_x_0_0_x_zz_zz = buffer_1010_spdd[359];

    auto g_y_0_x_0_0_y_xx_xx = buffer_1010_spdd[360];

    auto g_y_0_x_0_0_y_xx_xy = buffer_1010_spdd[361];

    auto g_y_0_x_0_0_y_xx_xz = buffer_1010_spdd[362];

    auto g_y_0_x_0_0_y_xx_yy = buffer_1010_spdd[363];

    auto g_y_0_x_0_0_y_xx_yz = buffer_1010_spdd[364];

    auto g_y_0_x_0_0_y_xx_zz = buffer_1010_spdd[365];

    auto g_y_0_x_0_0_y_xy_xx = buffer_1010_spdd[366];

    auto g_y_0_x_0_0_y_xy_xy = buffer_1010_spdd[367];

    auto g_y_0_x_0_0_y_xy_xz = buffer_1010_spdd[368];

    auto g_y_0_x_0_0_y_xy_yy = buffer_1010_spdd[369];

    auto g_y_0_x_0_0_y_xy_yz = buffer_1010_spdd[370];

    auto g_y_0_x_0_0_y_xy_zz = buffer_1010_spdd[371];

    auto g_y_0_x_0_0_y_xz_xx = buffer_1010_spdd[372];

    auto g_y_0_x_0_0_y_xz_xy = buffer_1010_spdd[373];

    auto g_y_0_x_0_0_y_xz_xz = buffer_1010_spdd[374];

    auto g_y_0_x_0_0_y_xz_yy = buffer_1010_spdd[375];

    auto g_y_0_x_0_0_y_xz_yz = buffer_1010_spdd[376];

    auto g_y_0_x_0_0_y_xz_zz = buffer_1010_spdd[377];

    auto g_y_0_x_0_0_y_yy_xx = buffer_1010_spdd[378];

    auto g_y_0_x_0_0_y_yy_xy = buffer_1010_spdd[379];

    auto g_y_0_x_0_0_y_yy_xz = buffer_1010_spdd[380];

    auto g_y_0_x_0_0_y_yy_yy = buffer_1010_spdd[381];

    auto g_y_0_x_0_0_y_yy_yz = buffer_1010_spdd[382];

    auto g_y_0_x_0_0_y_yy_zz = buffer_1010_spdd[383];

    auto g_y_0_x_0_0_y_yz_xx = buffer_1010_spdd[384];

    auto g_y_0_x_0_0_y_yz_xy = buffer_1010_spdd[385];

    auto g_y_0_x_0_0_y_yz_xz = buffer_1010_spdd[386];

    auto g_y_0_x_0_0_y_yz_yy = buffer_1010_spdd[387];

    auto g_y_0_x_0_0_y_yz_yz = buffer_1010_spdd[388];

    auto g_y_0_x_0_0_y_yz_zz = buffer_1010_spdd[389];

    auto g_y_0_x_0_0_y_zz_xx = buffer_1010_spdd[390];

    auto g_y_0_x_0_0_y_zz_xy = buffer_1010_spdd[391];

    auto g_y_0_x_0_0_y_zz_xz = buffer_1010_spdd[392];

    auto g_y_0_x_0_0_y_zz_yy = buffer_1010_spdd[393];

    auto g_y_0_x_0_0_y_zz_yz = buffer_1010_spdd[394];

    auto g_y_0_x_0_0_y_zz_zz = buffer_1010_spdd[395];

    auto g_y_0_x_0_0_z_xx_xx = buffer_1010_spdd[396];

    auto g_y_0_x_0_0_z_xx_xy = buffer_1010_spdd[397];

    auto g_y_0_x_0_0_z_xx_xz = buffer_1010_spdd[398];

    auto g_y_0_x_0_0_z_xx_yy = buffer_1010_spdd[399];

    auto g_y_0_x_0_0_z_xx_yz = buffer_1010_spdd[400];

    auto g_y_0_x_0_0_z_xx_zz = buffer_1010_spdd[401];

    auto g_y_0_x_0_0_z_xy_xx = buffer_1010_spdd[402];

    auto g_y_0_x_0_0_z_xy_xy = buffer_1010_spdd[403];

    auto g_y_0_x_0_0_z_xy_xz = buffer_1010_spdd[404];

    auto g_y_0_x_0_0_z_xy_yy = buffer_1010_spdd[405];

    auto g_y_0_x_0_0_z_xy_yz = buffer_1010_spdd[406];

    auto g_y_0_x_0_0_z_xy_zz = buffer_1010_spdd[407];

    auto g_y_0_x_0_0_z_xz_xx = buffer_1010_spdd[408];

    auto g_y_0_x_0_0_z_xz_xy = buffer_1010_spdd[409];

    auto g_y_0_x_0_0_z_xz_xz = buffer_1010_spdd[410];

    auto g_y_0_x_0_0_z_xz_yy = buffer_1010_spdd[411];

    auto g_y_0_x_0_0_z_xz_yz = buffer_1010_spdd[412];

    auto g_y_0_x_0_0_z_xz_zz = buffer_1010_spdd[413];

    auto g_y_0_x_0_0_z_yy_xx = buffer_1010_spdd[414];

    auto g_y_0_x_0_0_z_yy_xy = buffer_1010_spdd[415];

    auto g_y_0_x_0_0_z_yy_xz = buffer_1010_spdd[416];

    auto g_y_0_x_0_0_z_yy_yy = buffer_1010_spdd[417];

    auto g_y_0_x_0_0_z_yy_yz = buffer_1010_spdd[418];

    auto g_y_0_x_0_0_z_yy_zz = buffer_1010_spdd[419];

    auto g_y_0_x_0_0_z_yz_xx = buffer_1010_spdd[420];

    auto g_y_0_x_0_0_z_yz_xy = buffer_1010_spdd[421];

    auto g_y_0_x_0_0_z_yz_xz = buffer_1010_spdd[422];

    auto g_y_0_x_0_0_z_yz_yy = buffer_1010_spdd[423];

    auto g_y_0_x_0_0_z_yz_yz = buffer_1010_spdd[424];

    auto g_y_0_x_0_0_z_yz_zz = buffer_1010_spdd[425];

    auto g_y_0_x_0_0_z_zz_xx = buffer_1010_spdd[426];

    auto g_y_0_x_0_0_z_zz_xy = buffer_1010_spdd[427];

    auto g_y_0_x_0_0_z_zz_xz = buffer_1010_spdd[428];

    auto g_y_0_x_0_0_z_zz_yy = buffer_1010_spdd[429];

    auto g_y_0_x_0_0_z_zz_yz = buffer_1010_spdd[430];

    auto g_y_0_x_0_0_z_zz_zz = buffer_1010_spdd[431];

    auto g_y_0_y_0_0_x_xx_xx = buffer_1010_spdd[432];

    auto g_y_0_y_0_0_x_xx_xy = buffer_1010_spdd[433];

    auto g_y_0_y_0_0_x_xx_xz = buffer_1010_spdd[434];

    auto g_y_0_y_0_0_x_xx_yy = buffer_1010_spdd[435];

    auto g_y_0_y_0_0_x_xx_yz = buffer_1010_spdd[436];

    auto g_y_0_y_0_0_x_xx_zz = buffer_1010_spdd[437];

    auto g_y_0_y_0_0_x_xy_xx = buffer_1010_spdd[438];

    auto g_y_0_y_0_0_x_xy_xy = buffer_1010_spdd[439];

    auto g_y_0_y_0_0_x_xy_xz = buffer_1010_spdd[440];

    auto g_y_0_y_0_0_x_xy_yy = buffer_1010_spdd[441];

    auto g_y_0_y_0_0_x_xy_yz = buffer_1010_spdd[442];

    auto g_y_0_y_0_0_x_xy_zz = buffer_1010_spdd[443];

    auto g_y_0_y_0_0_x_xz_xx = buffer_1010_spdd[444];

    auto g_y_0_y_0_0_x_xz_xy = buffer_1010_spdd[445];

    auto g_y_0_y_0_0_x_xz_xz = buffer_1010_spdd[446];

    auto g_y_0_y_0_0_x_xz_yy = buffer_1010_spdd[447];

    auto g_y_0_y_0_0_x_xz_yz = buffer_1010_spdd[448];

    auto g_y_0_y_0_0_x_xz_zz = buffer_1010_spdd[449];

    auto g_y_0_y_0_0_x_yy_xx = buffer_1010_spdd[450];

    auto g_y_0_y_0_0_x_yy_xy = buffer_1010_spdd[451];

    auto g_y_0_y_0_0_x_yy_xz = buffer_1010_spdd[452];

    auto g_y_0_y_0_0_x_yy_yy = buffer_1010_spdd[453];

    auto g_y_0_y_0_0_x_yy_yz = buffer_1010_spdd[454];

    auto g_y_0_y_0_0_x_yy_zz = buffer_1010_spdd[455];

    auto g_y_0_y_0_0_x_yz_xx = buffer_1010_spdd[456];

    auto g_y_0_y_0_0_x_yz_xy = buffer_1010_spdd[457];

    auto g_y_0_y_0_0_x_yz_xz = buffer_1010_spdd[458];

    auto g_y_0_y_0_0_x_yz_yy = buffer_1010_spdd[459];

    auto g_y_0_y_0_0_x_yz_yz = buffer_1010_spdd[460];

    auto g_y_0_y_0_0_x_yz_zz = buffer_1010_spdd[461];

    auto g_y_0_y_0_0_x_zz_xx = buffer_1010_spdd[462];

    auto g_y_0_y_0_0_x_zz_xy = buffer_1010_spdd[463];

    auto g_y_0_y_0_0_x_zz_xz = buffer_1010_spdd[464];

    auto g_y_0_y_0_0_x_zz_yy = buffer_1010_spdd[465];

    auto g_y_0_y_0_0_x_zz_yz = buffer_1010_spdd[466];

    auto g_y_0_y_0_0_x_zz_zz = buffer_1010_spdd[467];

    auto g_y_0_y_0_0_y_xx_xx = buffer_1010_spdd[468];

    auto g_y_0_y_0_0_y_xx_xy = buffer_1010_spdd[469];

    auto g_y_0_y_0_0_y_xx_xz = buffer_1010_spdd[470];

    auto g_y_0_y_0_0_y_xx_yy = buffer_1010_spdd[471];

    auto g_y_0_y_0_0_y_xx_yz = buffer_1010_spdd[472];

    auto g_y_0_y_0_0_y_xx_zz = buffer_1010_spdd[473];

    auto g_y_0_y_0_0_y_xy_xx = buffer_1010_spdd[474];

    auto g_y_0_y_0_0_y_xy_xy = buffer_1010_spdd[475];

    auto g_y_0_y_0_0_y_xy_xz = buffer_1010_spdd[476];

    auto g_y_0_y_0_0_y_xy_yy = buffer_1010_spdd[477];

    auto g_y_0_y_0_0_y_xy_yz = buffer_1010_spdd[478];

    auto g_y_0_y_0_0_y_xy_zz = buffer_1010_spdd[479];

    auto g_y_0_y_0_0_y_xz_xx = buffer_1010_spdd[480];

    auto g_y_0_y_0_0_y_xz_xy = buffer_1010_spdd[481];

    auto g_y_0_y_0_0_y_xz_xz = buffer_1010_spdd[482];

    auto g_y_0_y_0_0_y_xz_yy = buffer_1010_spdd[483];

    auto g_y_0_y_0_0_y_xz_yz = buffer_1010_spdd[484];

    auto g_y_0_y_0_0_y_xz_zz = buffer_1010_spdd[485];

    auto g_y_0_y_0_0_y_yy_xx = buffer_1010_spdd[486];

    auto g_y_0_y_0_0_y_yy_xy = buffer_1010_spdd[487];

    auto g_y_0_y_0_0_y_yy_xz = buffer_1010_spdd[488];

    auto g_y_0_y_0_0_y_yy_yy = buffer_1010_spdd[489];

    auto g_y_0_y_0_0_y_yy_yz = buffer_1010_spdd[490];

    auto g_y_0_y_0_0_y_yy_zz = buffer_1010_spdd[491];

    auto g_y_0_y_0_0_y_yz_xx = buffer_1010_spdd[492];

    auto g_y_0_y_0_0_y_yz_xy = buffer_1010_spdd[493];

    auto g_y_0_y_0_0_y_yz_xz = buffer_1010_spdd[494];

    auto g_y_0_y_0_0_y_yz_yy = buffer_1010_spdd[495];

    auto g_y_0_y_0_0_y_yz_yz = buffer_1010_spdd[496];

    auto g_y_0_y_0_0_y_yz_zz = buffer_1010_spdd[497];

    auto g_y_0_y_0_0_y_zz_xx = buffer_1010_spdd[498];

    auto g_y_0_y_0_0_y_zz_xy = buffer_1010_spdd[499];

    auto g_y_0_y_0_0_y_zz_xz = buffer_1010_spdd[500];

    auto g_y_0_y_0_0_y_zz_yy = buffer_1010_spdd[501];

    auto g_y_0_y_0_0_y_zz_yz = buffer_1010_spdd[502];

    auto g_y_0_y_0_0_y_zz_zz = buffer_1010_spdd[503];

    auto g_y_0_y_0_0_z_xx_xx = buffer_1010_spdd[504];

    auto g_y_0_y_0_0_z_xx_xy = buffer_1010_spdd[505];

    auto g_y_0_y_0_0_z_xx_xz = buffer_1010_spdd[506];

    auto g_y_0_y_0_0_z_xx_yy = buffer_1010_spdd[507];

    auto g_y_0_y_0_0_z_xx_yz = buffer_1010_spdd[508];

    auto g_y_0_y_0_0_z_xx_zz = buffer_1010_spdd[509];

    auto g_y_0_y_0_0_z_xy_xx = buffer_1010_spdd[510];

    auto g_y_0_y_0_0_z_xy_xy = buffer_1010_spdd[511];

    auto g_y_0_y_0_0_z_xy_xz = buffer_1010_spdd[512];

    auto g_y_0_y_0_0_z_xy_yy = buffer_1010_spdd[513];

    auto g_y_0_y_0_0_z_xy_yz = buffer_1010_spdd[514];

    auto g_y_0_y_0_0_z_xy_zz = buffer_1010_spdd[515];

    auto g_y_0_y_0_0_z_xz_xx = buffer_1010_spdd[516];

    auto g_y_0_y_0_0_z_xz_xy = buffer_1010_spdd[517];

    auto g_y_0_y_0_0_z_xz_xz = buffer_1010_spdd[518];

    auto g_y_0_y_0_0_z_xz_yy = buffer_1010_spdd[519];

    auto g_y_0_y_0_0_z_xz_yz = buffer_1010_spdd[520];

    auto g_y_0_y_0_0_z_xz_zz = buffer_1010_spdd[521];

    auto g_y_0_y_0_0_z_yy_xx = buffer_1010_spdd[522];

    auto g_y_0_y_0_0_z_yy_xy = buffer_1010_spdd[523];

    auto g_y_0_y_0_0_z_yy_xz = buffer_1010_spdd[524];

    auto g_y_0_y_0_0_z_yy_yy = buffer_1010_spdd[525];

    auto g_y_0_y_0_0_z_yy_yz = buffer_1010_spdd[526];

    auto g_y_0_y_0_0_z_yy_zz = buffer_1010_spdd[527];

    auto g_y_0_y_0_0_z_yz_xx = buffer_1010_spdd[528];

    auto g_y_0_y_0_0_z_yz_xy = buffer_1010_spdd[529];

    auto g_y_0_y_0_0_z_yz_xz = buffer_1010_spdd[530];

    auto g_y_0_y_0_0_z_yz_yy = buffer_1010_spdd[531];

    auto g_y_0_y_0_0_z_yz_yz = buffer_1010_spdd[532];

    auto g_y_0_y_0_0_z_yz_zz = buffer_1010_spdd[533];

    auto g_y_0_y_0_0_z_zz_xx = buffer_1010_spdd[534];

    auto g_y_0_y_0_0_z_zz_xy = buffer_1010_spdd[535];

    auto g_y_0_y_0_0_z_zz_xz = buffer_1010_spdd[536];

    auto g_y_0_y_0_0_z_zz_yy = buffer_1010_spdd[537];

    auto g_y_0_y_0_0_z_zz_yz = buffer_1010_spdd[538];

    auto g_y_0_y_0_0_z_zz_zz = buffer_1010_spdd[539];

    auto g_y_0_z_0_0_x_xx_xx = buffer_1010_spdd[540];

    auto g_y_0_z_0_0_x_xx_xy = buffer_1010_spdd[541];

    auto g_y_0_z_0_0_x_xx_xz = buffer_1010_spdd[542];

    auto g_y_0_z_0_0_x_xx_yy = buffer_1010_spdd[543];

    auto g_y_0_z_0_0_x_xx_yz = buffer_1010_spdd[544];

    auto g_y_0_z_0_0_x_xx_zz = buffer_1010_spdd[545];

    auto g_y_0_z_0_0_x_xy_xx = buffer_1010_spdd[546];

    auto g_y_0_z_0_0_x_xy_xy = buffer_1010_spdd[547];

    auto g_y_0_z_0_0_x_xy_xz = buffer_1010_spdd[548];

    auto g_y_0_z_0_0_x_xy_yy = buffer_1010_spdd[549];

    auto g_y_0_z_0_0_x_xy_yz = buffer_1010_spdd[550];

    auto g_y_0_z_0_0_x_xy_zz = buffer_1010_spdd[551];

    auto g_y_0_z_0_0_x_xz_xx = buffer_1010_spdd[552];

    auto g_y_0_z_0_0_x_xz_xy = buffer_1010_spdd[553];

    auto g_y_0_z_0_0_x_xz_xz = buffer_1010_spdd[554];

    auto g_y_0_z_0_0_x_xz_yy = buffer_1010_spdd[555];

    auto g_y_0_z_0_0_x_xz_yz = buffer_1010_spdd[556];

    auto g_y_0_z_0_0_x_xz_zz = buffer_1010_spdd[557];

    auto g_y_0_z_0_0_x_yy_xx = buffer_1010_spdd[558];

    auto g_y_0_z_0_0_x_yy_xy = buffer_1010_spdd[559];

    auto g_y_0_z_0_0_x_yy_xz = buffer_1010_spdd[560];

    auto g_y_0_z_0_0_x_yy_yy = buffer_1010_spdd[561];

    auto g_y_0_z_0_0_x_yy_yz = buffer_1010_spdd[562];

    auto g_y_0_z_0_0_x_yy_zz = buffer_1010_spdd[563];

    auto g_y_0_z_0_0_x_yz_xx = buffer_1010_spdd[564];

    auto g_y_0_z_0_0_x_yz_xy = buffer_1010_spdd[565];

    auto g_y_0_z_0_0_x_yz_xz = buffer_1010_spdd[566];

    auto g_y_0_z_0_0_x_yz_yy = buffer_1010_spdd[567];

    auto g_y_0_z_0_0_x_yz_yz = buffer_1010_spdd[568];

    auto g_y_0_z_0_0_x_yz_zz = buffer_1010_spdd[569];

    auto g_y_0_z_0_0_x_zz_xx = buffer_1010_spdd[570];

    auto g_y_0_z_0_0_x_zz_xy = buffer_1010_spdd[571];

    auto g_y_0_z_0_0_x_zz_xz = buffer_1010_spdd[572];

    auto g_y_0_z_0_0_x_zz_yy = buffer_1010_spdd[573];

    auto g_y_0_z_0_0_x_zz_yz = buffer_1010_spdd[574];

    auto g_y_0_z_0_0_x_zz_zz = buffer_1010_spdd[575];

    auto g_y_0_z_0_0_y_xx_xx = buffer_1010_spdd[576];

    auto g_y_0_z_0_0_y_xx_xy = buffer_1010_spdd[577];

    auto g_y_0_z_0_0_y_xx_xz = buffer_1010_spdd[578];

    auto g_y_0_z_0_0_y_xx_yy = buffer_1010_spdd[579];

    auto g_y_0_z_0_0_y_xx_yz = buffer_1010_spdd[580];

    auto g_y_0_z_0_0_y_xx_zz = buffer_1010_spdd[581];

    auto g_y_0_z_0_0_y_xy_xx = buffer_1010_spdd[582];

    auto g_y_0_z_0_0_y_xy_xy = buffer_1010_spdd[583];

    auto g_y_0_z_0_0_y_xy_xz = buffer_1010_spdd[584];

    auto g_y_0_z_0_0_y_xy_yy = buffer_1010_spdd[585];

    auto g_y_0_z_0_0_y_xy_yz = buffer_1010_spdd[586];

    auto g_y_0_z_0_0_y_xy_zz = buffer_1010_spdd[587];

    auto g_y_0_z_0_0_y_xz_xx = buffer_1010_spdd[588];

    auto g_y_0_z_0_0_y_xz_xy = buffer_1010_spdd[589];

    auto g_y_0_z_0_0_y_xz_xz = buffer_1010_spdd[590];

    auto g_y_0_z_0_0_y_xz_yy = buffer_1010_spdd[591];

    auto g_y_0_z_0_0_y_xz_yz = buffer_1010_spdd[592];

    auto g_y_0_z_0_0_y_xz_zz = buffer_1010_spdd[593];

    auto g_y_0_z_0_0_y_yy_xx = buffer_1010_spdd[594];

    auto g_y_0_z_0_0_y_yy_xy = buffer_1010_spdd[595];

    auto g_y_0_z_0_0_y_yy_xz = buffer_1010_spdd[596];

    auto g_y_0_z_0_0_y_yy_yy = buffer_1010_spdd[597];

    auto g_y_0_z_0_0_y_yy_yz = buffer_1010_spdd[598];

    auto g_y_0_z_0_0_y_yy_zz = buffer_1010_spdd[599];

    auto g_y_0_z_0_0_y_yz_xx = buffer_1010_spdd[600];

    auto g_y_0_z_0_0_y_yz_xy = buffer_1010_spdd[601];

    auto g_y_0_z_0_0_y_yz_xz = buffer_1010_spdd[602];

    auto g_y_0_z_0_0_y_yz_yy = buffer_1010_spdd[603];

    auto g_y_0_z_0_0_y_yz_yz = buffer_1010_spdd[604];

    auto g_y_0_z_0_0_y_yz_zz = buffer_1010_spdd[605];

    auto g_y_0_z_0_0_y_zz_xx = buffer_1010_spdd[606];

    auto g_y_0_z_0_0_y_zz_xy = buffer_1010_spdd[607];

    auto g_y_0_z_0_0_y_zz_xz = buffer_1010_spdd[608];

    auto g_y_0_z_0_0_y_zz_yy = buffer_1010_spdd[609];

    auto g_y_0_z_0_0_y_zz_yz = buffer_1010_spdd[610];

    auto g_y_0_z_0_0_y_zz_zz = buffer_1010_spdd[611];

    auto g_y_0_z_0_0_z_xx_xx = buffer_1010_spdd[612];

    auto g_y_0_z_0_0_z_xx_xy = buffer_1010_spdd[613];

    auto g_y_0_z_0_0_z_xx_xz = buffer_1010_spdd[614];

    auto g_y_0_z_0_0_z_xx_yy = buffer_1010_spdd[615];

    auto g_y_0_z_0_0_z_xx_yz = buffer_1010_spdd[616];

    auto g_y_0_z_0_0_z_xx_zz = buffer_1010_spdd[617];

    auto g_y_0_z_0_0_z_xy_xx = buffer_1010_spdd[618];

    auto g_y_0_z_0_0_z_xy_xy = buffer_1010_spdd[619];

    auto g_y_0_z_0_0_z_xy_xz = buffer_1010_spdd[620];

    auto g_y_0_z_0_0_z_xy_yy = buffer_1010_spdd[621];

    auto g_y_0_z_0_0_z_xy_yz = buffer_1010_spdd[622];

    auto g_y_0_z_0_0_z_xy_zz = buffer_1010_spdd[623];

    auto g_y_0_z_0_0_z_xz_xx = buffer_1010_spdd[624];

    auto g_y_0_z_0_0_z_xz_xy = buffer_1010_spdd[625];

    auto g_y_0_z_0_0_z_xz_xz = buffer_1010_spdd[626];

    auto g_y_0_z_0_0_z_xz_yy = buffer_1010_spdd[627];

    auto g_y_0_z_0_0_z_xz_yz = buffer_1010_spdd[628];

    auto g_y_0_z_0_0_z_xz_zz = buffer_1010_spdd[629];

    auto g_y_0_z_0_0_z_yy_xx = buffer_1010_spdd[630];

    auto g_y_0_z_0_0_z_yy_xy = buffer_1010_spdd[631];

    auto g_y_0_z_0_0_z_yy_xz = buffer_1010_spdd[632];

    auto g_y_0_z_0_0_z_yy_yy = buffer_1010_spdd[633];

    auto g_y_0_z_0_0_z_yy_yz = buffer_1010_spdd[634];

    auto g_y_0_z_0_0_z_yy_zz = buffer_1010_spdd[635];

    auto g_y_0_z_0_0_z_yz_xx = buffer_1010_spdd[636];

    auto g_y_0_z_0_0_z_yz_xy = buffer_1010_spdd[637];

    auto g_y_0_z_0_0_z_yz_xz = buffer_1010_spdd[638];

    auto g_y_0_z_0_0_z_yz_yy = buffer_1010_spdd[639];

    auto g_y_0_z_0_0_z_yz_yz = buffer_1010_spdd[640];

    auto g_y_0_z_0_0_z_yz_zz = buffer_1010_spdd[641];

    auto g_y_0_z_0_0_z_zz_xx = buffer_1010_spdd[642];

    auto g_y_0_z_0_0_z_zz_xy = buffer_1010_spdd[643];

    auto g_y_0_z_0_0_z_zz_xz = buffer_1010_spdd[644];

    auto g_y_0_z_0_0_z_zz_yy = buffer_1010_spdd[645];

    auto g_y_0_z_0_0_z_zz_yz = buffer_1010_spdd[646];

    auto g_y_0_z_0_0_z_zz_zz = buffer_1010_spdd[647];

    auto g_z_0_x_0_0_x_xx_xx = buffer_1010_spdd[648];

    auto g_z_0_x_0_0_x_xx_xy = buffer_1010_spdd[649];

    auto g_z_0_x_0_0_x_xx_xz = buffer_1010_spdd[650];

    auto g_z_0_x_0_0_x_xx_yy = buffer_1010_spdd[651];

    auto g_z_0_x_0_0_x_xx_yz = buffer_1010_spdd[652];

    auto g_z_0_x_0_0_x_xx_zz = buffer_1010_spdd[653];

    auto g_z_0_x_0_0_x_xy_xx = buffer_1010_spdd[654];

    auto g_z_0_x_0_0_x_xy_xy = buffer_1010_spdd[655];

    auto g_z_0_x_0_0_x_xy_xz = buffer_1010_spdd[656];

    auto g_z_0_x_0_0_x_xy_yy = buffer_1010_spdd[657];

    auto g_z_0_x_0_0_x_xy_yz = buffer_1010_spdd[658];

    auto g_z_0_x_0_0_x_xy_zz = buffer_1010_spdd[659];

    auto g_z_0_x_0_0_x_xz_xx = buffer_1010_spdd[660];

    auto g_z_0_x_0_0_x_xz_xy = buffer_1010_spdd[661];

    auto g_z_0_x_0_0_x_xz_xz = buffer_1010_spdd[662];

    auto g_z_0_x_0_0_x_xz_yy = buffer_1010_spdd[663];

    auto g_z_0_x_0_0_x_xz_yz = buffer_1010_spdd[664];

    auto g_z_0_x_0_0_x_xz_zz = buffer_1010_spdd[665];

    auto g_z_0_x_0_0_x_yy_xx = buffer_1010_spdd[666];

    auto g_z_0_x_0_0_x_yy_xy = buffer_1010_spdd[667];

    auto g_z_0_x_0_0_x_yy_xz = buffer_1010_spdd[668];

    auto g_z_0_x_0_0_x_yy_yy = buffer_1010_spdd[669];

    auto g_z_0_x_0_0_x_yy_yz = buffer_1010_spdd[670];

    auto g_z_0_x_0_0_x_yy_zz = buffer_1010_spdd[671];

    auto g_z_0_x_0_0_x_yz_xx = buffer_1010_spdd[672];

    auto g_z_0_x_0_0_x_yz_xy = buffer_1010_spdd[673];

    auto g_z_0_x_0_0_x_yz_xz = buffer_1010_spdd[674];

    auto g_z_0_x_0_0_x_yz_yy = buffer_1010_spdd[675];

    auto g_z_0_x_0_0_x_yz_yz = buffer_1010_spdd[676];

    auto g_z_0_x_0_0_x_yz_zz = buffer_1010_spdd[677];

    auto g_z_0_x_0_0_x_zz_xx = buffer_1010_spdd[678];

    auto g_z_0_x_0_0_x_zz_xy = buffer_1010_spdd[679];

    auto g_z_0_x_0_0_x_zz_xz = buffer_1010_spdd[680];

    auto g_z_0_x_0_0_x_zz_yy = buffer_1010_spdd[681];

    auto g_z_0_x_0_0_x_zz_yz = buffer_1010_spdd[682];

    auto g_z_0_x_0_0_x_zz_zz = buffer_1010_spdd[683];

    auto g_z_0_x_0_0_y_xx_xx = buffer_1010_spdd[684];

    auto g_z_0_x_0_0_y_xx_xy = buffer_1010_spdd[685];

    auto g_z_0_x_0_0_y_xx_xz = buffer_1010_spdd[686];

    auto g_z_0_x_0_0_y_xx_yy = buffer_1010_spdd[687];

    auto g_z_0_x_0_0_y_xx_yz = buffer_1010_spdd[688];

    auto g_z_0_x_0_0_y_xx_zz = buffer_1010_spdd[689];

    auto g_z_0_x_0_0_y_xy_xx = buffer_1010_spdd[690];

    auto g_z_0_x_0_0_y_xy_xy = buffer_1010_spdd[691];

    auto g_z_0_x_0_0_y_xy_xz = buffer_1010_spdd[692];

    auto g_z_0_x_0_0_y_xy_yy = buffer_1010_spdd[693];

    auto g_z_0_x_0_0_y_xy_yz = buffer_1010_spdd[694];

    auto g_z_0_x_0_0_y_xy_zz = buffer_1010_spdd[695];

    auto g_z_0_x_0_0_y_xz_xx = buffer_1010_spdd[696];

    auto g_z_0_x_0_0_y_xz_xy = buffer_1010_spdd[697];

    auto g_z_0_x_0_0_y_xz_xz = buffer_1010_spdd[698];

    auto g_z_0_x_0_0_y_xz_yy = buffer_1010_spdd[699];

    auto g_z_0_x_0_0_y_xz_yz = buffer_1010_spdd[700];

    auto g_z_0_x_0_0_y_xz_zz = buffer_1010_spdd[701];

    auto g_z_0_x_0_0_y_yy_xx = buffer_1010_spdd[702];

    auto g_z_0_x_0_0_y_yy_xy = buffer_1010_spdd[703];

    auto g_z_0_x_0_0_y_yy_xz = buffer_1010_spdd[704];

    auto g_z_0_x_0_0_y_yy_yy = buffer_1010_spdd[705];

    auto g_z_0_x_0_0_y_yy_yz = buffer_1010_spdd[706];

    auto g_z_0_x_0_0_y_yy_zz = buffer_1010_spdd[707];

    auto g_z_0_x_0_0_y_yz_xx = buffer_1010_spdd[708];

    auto g_z_0_x_0_0_y_yz_xy = buffer_1010_spdd[709];

    auto g_z_0_x_0_0_y_yz_xz = buffer_1010_spdd[710];

    auto g_z_0_x_0_0_y_yz_yy = buffer_1010_spdd[711];

    auto g_z_0_x_0_0_y_yz_yz = buffer_1010_spdd[712];

    auto g_z_0_x_0_0_y_yz_zz = buffer_1010_spdd[713];

    auto g_z_0_x_0_0_y_zz_xx = buffer_1010_spdd[714];

    auto g_z_0_x_0_0_y_zz_xy = buffer_1010_spdd[715];

    auto g_z_0_x_0_0_y_zz_xz = buffer_1010_spdd[716];

    auto g_z_0_x_0_0_y_zz_yy = buffer_1010_spdd[717];

    auto g_z_0_x_0_0_y_zz_yz = buffer_1010_spdd[718];

    auto g_z_0_x_0_0_y_zz_zz = buffer_1010_spdd[719];

    auto g_z_0_x_0_0_z_xx_xx = buffer_1010_spdd[720];

    auto g_z_0_x_0_0_z_xx_xy = buffer_1010_spdd[721];

    auto g_z_0_x_0_0_z_xx_xz = buffer_1010_spdd[722];

    auto g_z_0_x_0_0_z_xx_yy = buffer_1010_spdd[723];

    auto g_z_0_x_0_0_z_xx_yz = buffer_1010_spdd[724];

    auto g_z_0_x_0_0_z_xx_zz = buffer_1010_spdd[725];

    auto g_z_0_x_0_0_z_xy_xx = buffer_1010_spdd[726];

    auto g_z_0_x_0_0_z_xy_xy = buffer_1010_spdd[727];

    auto g_z_0_x_0_0_z_xy_xz = buffer_1010_spdd[728];

    auto g_z_0_x_0_0_z_xy_yy = buffer_1010_spdd[729];

    auto g_z_0_x_0_0_z_xy_yz = buffer_1010_spdd[730];

    auto g_z_0_x_0_0_z_xy_zz = buffer_1010_spdd[731];

    auto g_z_0_x_0_0_z_xz_xx = buffer_1010_spdd[732];

    auto g_z_0_x_0_0_z_xz_xy = buffer_1010_spdd[733];

    auto g_z_0_x_0_0_z_xz_xz = buffer_1010_spdd[734];

    auto g_z_0_x_0_0_z_xz_yy = buffer_1010_spdd[735];

    auto g_z_0_x_0_0_z_xz_yz = buffer_1010_spdd[736];

    auto g_z_0_x_0_0_z_xz_zz = buffer_1010_spdd[737];

    auto g_z_0_x_0_0_z_yy_xx = buffer_1010_spdd[738];

    auto g_z_0_x_0_0_z_yy_xy = buffer_1010_spdd[739];

    auto g_z_0_x_0_0_z_yy_xz = buffer_1010_spdd[740];

    auto g_z_0_x_0_0_z_yy_yy = buffer_1010_spdd[741];

    auto g_z_0_x_0_0_z_yy_yz = buffer_1010_spdd[742];

    auto g_z_0_x_0_0_z_yy_zz = buffer_1010_spdd[743];

    auto g_z_0_x_0_0_z_yz_xx = buffer_1010_spdd[744];

    auto g_z_0_x_0_0_z_yz_xy = buffer_1010_spdd[745];

    auto g_z_0_x_0_0_z_yz_xz = buffer_1010_spdd[746];

    auto g_z_0_x_0_0_z_yz_yy = buffer_1010_spdd[747];

    auto g_z_0_x_0_0_z_yz_yz = buffer_1010_spdd[748];

    auto g_z_0_x_0_0_z_yz_zz = buffer_1010_spdd[749];

    auto g_z_0_x_0_0_z_zz_xx = buffer_1010_spdd[750];

    auto g_z_0_x_0_0_z_zz_xy = buffer_1010_spdd[751];

    auto g_z_0_x_0_0_z_zz_xz = buffer_1010_spdd[752];

    auto g_z_0_x_0_0_z_zz_yy = buffer_1010_spdd[753];

    auto g_z_0_x_0_0_z_zz_yz = buffer_1010_spdd[754];

    auto g_z_0_x_0_0_z_zz_zz = buffer_1010_spdd[755];

    auto g_z_0_y_0_0_x_xx_xx = buffer_1010_spdd[756];

    auto g_z_0_y_0_0_x_xx_xy = buffer_1010_spdd[757];

    auto g_z_0_y_0_0_x_xx_xz = buffer_1010_spdd[758];

    auto g_z_0_y_0_0_x_xx_yy = buffer_1010_spdd[759];

    auto g_z_0_y_0_0_x_xx_yz = buffer_1010_spdd[760];

    auto g_z_0_y_0_0_x_xx_zz = buffer_1010_spdd[761];

    auto g_z_0_y_0_0_x_xy_xx = buffer_1010_spdd[762];

    auto g_z_0_y_0_0_x_xy_xy = buffer_1010_spdd[763];

    auto g_z_0_y_0_0_x_xy_xz = buffer_1010_spdd[764];

    auto g_z_0_y_0_0_x_xy_yy = buffer_1010_spdd[765];

    auto g_z_0_y_0_0_x_xy_yz = buffer_1010_spdd[766];

    auto g_z_0_y_0_0_x_xy_zz = buffer_1010_spdd[767];

    auto g_z_0_y_0_0_x_xz_xx = buffer_1010_spdd[768];

    auto g_z_0_y_0_0_x_xz_xy = buffer_1010_spdd[769];

    auto g_z_0_y_0_0_x_xz_xz = buffer_1010_spdd[770];

    auto g_z_0_y_0_0_x_xz_yy = buffer_1010_spdd[771];

    auto g_z_0_y_0_0_x_xz_yz = buffer_1010_spdd[772];

    auto g_z_0_y_0_0_x_xz_zz = buffer_1010_spdd[773];

    auto g_z_0_y_0_0_x_yy_xx = buffer_1010_spdd[774];

    auto g_z_0_y_0_0_x_yy_xy = buffer_1010_spdd[775];

    auto g_z_0_y_0_0_x_yy_xz = buffer_1010_spdd[776];

    auto g_z_0_y_0_0_x_yy_yy = buffer_1010_spdd[777];

    auto g_z_0_y_0_0_x_yy_yz = buffer_1010_spdd[778];

    auto g_z_0_y_0_0_x_yy_zz = buffer_1010_spdd[779];

    auto g_z_0_y_0_0_x_yz_xx = buffer_1010_spdd[780];

    auto g_z_0_y_0_0_x_yz_xy = buffer_1010_spdd[781];

    auto g_z_0_y_0_0_x_yz_xz = buffer_1010_spdd[782];

    auto g_z_0_y_0_0_x_yz_yy = buffer_1010_spdd[783];

    auto g_z_0_y_0_0_x_yz_yz = buffer_1010_spdd[784];

    auto g_z_0_y_0_0_x_yz_zz = buffer_1010_spdd[785];

    auto g_z_0_y_0_0_x_zz_xx = buffer_1010_spdd[786];

    auto g_z_0_y_0_0_x_zz_xy = buffer_1010_spdd[787];

    auto g_z_0_y_0_0_x_zz_xz = buffer_1010_spdd[788];

    auto g_z_0_y_0_0_x_zz_yy = buffer_1010_spdd[789];

    auto g_z_0_y_0_0_x_zz_yz = buffer_1010_spdd[790];

    auto g_z_0_y_0_0_x_zz_zz = buffer_1010_spdd[791];

    auto g_z_0_y_0_0_y_xx_xx = buffer_1010_spdd[792];

    auto g_z_0_y_0_0_y_xx_xy = buffer_1010_spdd[793];

    auto g_z_0_y_0_0_y_xx_xz = buffer_1010_spdd[794];

    auto g_z_0_y_0_0_y_xx_yy = buffer_1010_spdd[795];

    auto g_z_0_y_0_0_y_xx_yz = buffer_1010_spdd[796];

    auto g_z_0_y_0_0_y_xx_zz = buffer_1010_spdd[797];

    auto g_z_0_y_0_0_y_xy_xx = buffer_1010_spdd[798];

    auto g_z_0_y_0_0_y_xy_xy = buffer_1010_spdd[799];

    auto g_z_0_y_0_0_y_xy_xz = buffer_1010_spdd[800];

    auto g_z_0_y_0_0_y_xy_yy = buffer_1010_spdd[801];

    auto g_z_0_y_0_0_y_xy_yz = buffer_1010_spdd[802];

    auto g_z_0_y_0_0_y_xy_zz = buffer_1010_spdd[803];

    auto g_z_0_y_0_0_y_xz_xx = buffer_1010_spdd[804];

    auto g_z_0_y_0_0_y_xz_xy = buffer_1010_spdd[805];

    auto g_z_0_y_0_0_y_xz_xz = buffer_1010_spdd[806];

    auto g_z_0_y_0_0_y_xz_yy = buffer_1010_spdd[807];

    auto g_z_0_y_0_0_y_xz_yz = buffer_1010_spdd[808];

    auto g_z_0_y_0_0_y_xz_zz = buffer_1010_spdd[809];

    auto g_z_0_y_0_0_y_yy_xx = buffer_1010_spdd[810];

    auto g_z_0_y_0_0_y_yy_xy = buffer_1010_spdd[811];

    auto g_z_0_y_0_0_y_yy_xz = buffer_1010_spdd[812];

    auto g_z_0_y_0_0_y_yy_yy = buffer_1010_spdd[813];

    auto g_z_0_y_0_0_y_yy_yz = buffer_1010_spdd[814];

    auto g_z_0_y_0_0_y_yy_zz = buffer_1010_spdd[815];

    auto g_z_0_y_0_0_y_yz_xx = buffer_1010_spdd[816];

    auto g_z_0_y_0_0_y_yz_xy = buffer_1010_spdd[817];

    auto g_z_0_y_0_0_y_yz_xz = buffer_1010_spdd[818];

    auto g_z_0_y_0_0_y_yz_yy = buffer_1010_spdd[819];

    auto g_z_0_y_0_0_y_yz_yz = buffer_1010_spdd[820];

    auto g_z_0_y_0_0_y_yz_zz = buffer_1010_spdd[821];

    auto g_z_0_y_0_0_y_zz_xx = buffer_1010_spdd[822];

    auto g_z_0_y_0_0_y_zz_xy = buffer_1010_spdd[823];

    auto g_z_0_y_0_0_y_zz_xz = buffer_1010_spdd[824];

    auto g_z_0_y_0_0_y_zz_yy = buffer_1010_spdd[825];

    auto g_z_0_y_0_0_y_zz_yz = buffer_1010_spdd[826];

    auto g_z_0_y_0_0_y_zz_zz = buffer_1010_spdd[827];

    auto g_z_0_y_0_0_z_xx_xx = buffer_1010_spdd[828];

    auto g_z_0_y_0_0_z_xx_xy = buffer_1010_spdd[829];

    auto g_z_0_y_0_0_z_xx_xz = buffer_1010_spdd[830];

    auto g_z_0_y_0_0_z_xx_yy = buffer_1010_spdd[831];

    auto g_z_0_y_0_0_z_xx_yz = buffer_1010_spdd[832];

    auto g_z_0_y_0_0_z_xx_zz = buffer_1010_spdd[833];

    auto g_z_0_y_0_0_z_xy_xx = buffer_1010_spdd[834];

    auto g_z_0_y_0_0_z_xy_xy = buffer_1010_spdd[835];

    auto g_z_0_y_0_0_z_xy_xz = buffer_1010_spdd[836];

    auto g_z_0_y_0_0_z_xy_yy = buffer_1010_spdd[837];

    auto g_z_0_y_0_0_z_xy_yz = buffer_1010_spdd[838];

    auto g_z_0_y_0_0_z_xy_zz = buffer_1010_spdd[839];

    auto g_z_0_y_0_0_z_xz_xx = buffer_1010_spdd[840];

    auto g_z_0_y_0_0_z_xz_xy = buffer_1010_spdd[841];

    auto g_z_0_y_0_0_z_xz_xz = buffer_1010_spdd[842];

    auto g_z_0_y_0_0_z_xz_yy = buffer_1010_spdd[843];

    auto g_z_0_y_0_0_z_xz_yz = buffer_1010_spdd[844];

    auto g_z_0_y_0_0_z_xz_zz = buffer_1010_spdd[845];

    auto g_z_0_y_0_0_z_yy_xx = buffer_1010_spdd[846];

    auto g_z_0_y_0_0_z_yy_xy = buffer_1010_spdd[847];

    auto g_z_0_y_0_0_z_yy_xz = buffer_1010_spdd[848];

    auto g_z_0_y_0_0_z_yy_yy = buffer_1010_spdd[849];

    auto g_z_0_y_0_0_z_yy_yz = buffer_1010_spdd[850];

    auto g_z_0_y_0_0_z_yy_zz = buffer_1010_spdd[851];

    auto g_z_0_y_0_0_z_yz_xx = buffer_1010_spdd[852];

    auto g_z_0_y_0_0_z_yz_xy = buffer_1010_spdd[853];

    auto g_z_0_y_0_0_z_yz_xz = buffer_1010_spdd[854];

    auto g_z_0_y_0_0_z_yz_yy = buffer_1010_spdd[855];

    auto g_z_0_y_0_0_z_yz_yz = buffer_1010_spdd[856];

    auto g_z_0_y_0_0_z_yz_zz = buffer_1010_spdd[857];

    auto g_z_0_y_0_0_z_zz_xx = buffer_1010_spdd[858];

    auto g_z_0_y_0_0_z_zz_xy = buffer_1010_spdd[859];

    auto g_z_0_y_0_0_z_zz_xz = buffer_1010_spdd[860];

    auto g_z_0_y_0_0_z_zz_yy = buffer_1010_spdd[861];

    auto g_z_0_y_0_0_z_zz_yz = buffer_1010_spdd[862];

    auto g_z_0_y_0_0_z_zz_zz = buffer_1010_spdd[863];

    auto g_z_0_z_0_0_x_xx_xx = buffer_1010_spdd[864];

    auto g_z_0_z_0_0_x_xx_xy = buffer_1010_spdd[865];

    auto g_z_0_z_0_0_x_xx_xz = buffer_1010_spdd[866];

    auto g_z_0_z_0_0_x_xx_yy = buffer_1010_spdd[867];

    auto g_z_0_z_0_0_x_xx_yz = buffer_1010_spdd[868];

    auto g_z_0_z_0_0_x_xx_zz = buffer_1010_spdd[869];

    auto g_z_0_z_0_0_x_xy_xx = buffer_1010_spdd[870];

    auto g_z_0_z_0_0_x_xy_xy = buffer_1010_spdd[871];

    auto g_z_0_z_0_0_x_xy_xz = buffer_1010_spdd[872];

    auto g_z_0_z_0_0_x_xy_yy = buffer_1010_spdd[873];

    auto g_z_0_z_0_0_x_xy_yz = buffer_1010_spdd[874];

    auto g_z_0_z_0_0_x_xy_zz = buffer_1010_spdd[875];

    auto g_z_0_z_0_0_x_xz_xx = buffer_1010_spdd[876];

    auto g_z_0_z_0_0_x_xz_xy = buffer_1010_spdd[877];

    auto g_z_0_z_0_0_x_xz_xz = buffer_1010_spdd[878];

    auto g_z_0_z_0_0_x_xz_yy = buffer_1010_spdd[879];

    auto g_z_0_z_0_0_x_xz_yz = buffer_1010_spdd[880];

    auto g_z_0_z_0_0_x_xz_zz = buffer_1010_spdd[881];

    auto g_z_0_z_0_0_x_yy_xx = buffer_1010_spdd[882];

    auto g_z_0_z_0_0_x_yy_xy = buffer_1010_spdd[883];

    auto g_z_0_z_0_0_x_yy_xz = buffer_1010_spdd[884];

    auto g_z_0_z_0_0_x_yy_yy = buffer_1010_spdd[885];

    auto g_z_0_z_0_0_x_yy_yz = buffer_1010_spdd[886];

    auto g_z_0_z_0_0_x_yy_zz = buffer_1010_spdd[887];

    auto g_z_0_z_0_0_x_yz_xx = buffer_1010_spdd[888];

    auto g_z_0_z_0_0_x_yz_xy = buffer_1010_spdd[889];

    auto g_z_0_z_0_0_x_yz_xz = buffer_1010_spdd[890];

    auto g_z_0_z_0_0_x_yz_yy = buffer_1010_spdd[891];

    auto g_z_0_z_0_0_x_yz_yz = buffer_1010_spdd[892];

    auto g_z_0_z_0_0_x_yz_zz = buffer_1010_spdd[893];

    auto g_z_0_z_0_0_x_zz_xx = buffer_1010_spdd[894];

    auto g_z_0_z_0_0_x_zz_xy = buffer_1010_spdd[895];

    auto g_z_0_z_0_0_x_zz_xz = buffer_1010_spdd[896];

    auto g_z_0_z_0_0_x_zz_yy = buffer_1010_spdd[897];

    auto g_z_0_z_0_0_x_zz_yz = buffer_1010_spdd[898];

    auto g_z_0_z_0_0_x_zz_zz = buffer_1010_spdd[899];

    auto g_z_0_z_0_0_y_xx_xx = buffer_1010_spdd[900];

    auto g_z_0_z_0_0_y_xx_xy = buffer_1010_spdd[901];

    auto g_z_0_z_0_0_y_xx_xz = buffer_1010_spdd[902];

    auto g_z_0_z_0_0_y_xx_yy = buffer_1010_spdd[903];

    auto g_z_0_z_0_0_y_xx_yz = buffer_1010_spdd[904];

    auto g_z_0_z_0_0_y_xx_zz = buffer_1010_spdd[905];

    auto g_z_0_z_0_0_y_xy_xx = buffer_1010_spdd[906];

    auto g_z_0_z_0_0_y_xy_xy = buffer_1010_spdd[907];

    auto g_z_0_z_0_0_y_xy_xz = buffer_1010_spdd[908];

    auto g_z_0_z_0_0_y_xy_yy = buffer_1010_spdd[909];

    auto g_z_0_z_0_0_y_xy_yz = buffer_1010_spdd[910];

    auto g_z_0_z_0_0_y_xy_zz = buffer_1010_spdd[911];

    auto g_z_0_z_0_0_y_xz_xx = buffer_1010_spdd[912];

    auto g_z_0_z_0_0_y_xz_xy = buffer_1010_spdd[913];

    auto g_z_0_z_0_0_y_xz_xz = buffer_1010_spdd[914];

    auto g_z_0_z_0_0_y_xz_yy = buffer_1010_spdd[915];

    auto g_z_0_z_0_0_y_xz_yz = buffer_1010_spdd[916];

    auto g_z_0_z_0_0_y_xz_zz = buffer_1010_spdd[917];

    auto g_z_0_z_0_0_y_yy_xx = buffer_1010_spdd[918];

    auto g_z_0_z_0_0_y_yy_xy = buffer_1010_spdd[919];

    auto g_z_0_z_0_0_y_yy_xz = buffer_1010_spdd[920];

    auto g_z_0_z_0_0_y_yy_yy = buffer_1010_spdd[921];

    auto g_z_0_z_0_0_y_yy_yz = buffer_1010_spdd[922];

    auto g_z_0_z_0_0_y_yy_zz = buffer_1010_spdd[923];

    auto g_z_0_z_0_0_y_yz_xx = buffer_1010_spdd[924];

    auto g_z_0_z_0_0_y_yz_xy = buffer_1010_spdd[925];

    auto g_z_0_z_0_0_y_yz_xz = buffer_1010_spdd[926];

    auto g_z_0_z_0_0_y_yz_yy = buffer_1010_spdd[927];

    auto g_z_0_z_0_0_y_yz_yz = buffer_1010_spdd[928];

    auto g_z_0_z_0_0_y_yz_zz = buffer_1010_spdd[929];

    auto g_z_0_z_0_0_y_zz_xx = buffer_1010_spdd[930];

    auto g_z_0_z_0_0_y_zz_xy = buffer_1010_spdd[931];

    auto g_z_0_z_0_0_y_zz_xz = buffer_1010_spdd[932];

    auto g_z_0_z_0_0_y_zz_yy = buffer_1010_spdd[933];

    auto g_z_0_z_0_0_y_zz_yz = buffer_1010_spdd[934];

    auto g_z_0_z_0_0_y_zz_zz = buffer_1010_spdd[935];

    auto g_z_0_z_0_0_z_xx_xx = buffer_1010_spdd[936];

    auto g_z_0_z_0_0_z_xx_xy = buffer_1010_spdd[937];

    auto g_z_0_z_0_0_z_xx_xz = buffer_1010_spdd[938];

    auto g_z_0_z_0_0_z_xx_yy = buffer_1010_spdd[939];

    auto g_z_0_z_0_0_z_xx_yz = buffer_1010_spdd[940];

    auto g_z_0_z_0_0_z_xx_zz = buffer_1010_spdd[941];

    auto g_z_0_z_0_0_z_xy_xx = buffer_1010_spdd[942];

    auto g_z_0_z_0_0_z_xy_xy = buffer_1010_spdd[943];

    auto g_z_0_z_0_0_z_xy_xz = buffer_1010_spdd[944];

    auto g_z_0_z_0_0_z_xy_yy = buffer_1010_spdd[945];

    auto g_z_0_z_0_0_z_xy_yz = buffer_1010_spdd[946];

    auto g_z_0_z_0_0_z_xy_zz = buffer_1010_spdd[947];

    auto g_z_0_z_0_0_z_xz_xx = buffer_1010_spdd[948];

    auto g_z_0_z_0_0_z_xz_xy = buffer_1010_spdd[949];

    auto g_z_0_z_0_0_z_xz_xz = buffer_1010_spdd[950];

    auto g_z_0_z_0_0_z_xz_yy = buffer_1010_spdd[951];

    auto g_z_0_z_0_0_z_xz_yz = buffer_1010_spdd[952];

    auto g_z_0_z_0_0_z_xz_zz = buffer_1010_spdd[953];

    auto g_z_0_z_0_0_z_yy_xx = buffer_1010_spdd[954];

    auto g_z_0_z_0_0_z_yy_xy = buffer_1010_spdd[955];

    auto g_z_0_z_0_0_z_yy_xz = buffer_1010_spdd[956];

    auto g_z_0_z_0_0_z_yy_yy = buffer_1010_spdd[957];

    auto g_z_0_z_0_0_z_yy_yz = buffer_1010_spdd[958];

    auto g_z_0_z_0_0_z_yy_zz = buffer_1010_spdd[959];

    auto g_z_0_z_0_0_z_yz_xx = buffer_1010_spdd[960];

    auto g_z_0_z_0_0_z_yz_xy = buffer_1010_spdd[961];

    auto g_z_0_z_0_0_z_yz_xz = buffer_1010_spdd[962];

    auto g_z_0_z_0_0_z_yz_yy = buffer_1010_spdd[963];

    auto g_z_0_z_0_0_z_yz_yz = buffer_1010_spdd[964];

    auto g_z_0_z_0_0_z_yz_zz = buffer_1010_spdd[965];

    auto g_z_0_z_0_0_z_zz_xx = buffer_1010_spdd[966];

    auto g_z_0_z_0_0_z_zz_xy = buffer_1010_spdd[967];

    auto g_z_0_z_0_0_z_zz_xz = buffer_1010_spdd[968];

    auto g_z_0_z_0_0_z_zz_yy = buffer_1010_spdd[969];

    auto g_z_0_z_0_0_z_zz_yz = buffer_1010_spdd[970];

    auto g_z_0_z_0_0_z_zz_zz = buffer_1010_spdd[971];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_0_x_0_0_x_xx_xx, g_x_0_x_0_0_x_xx_xy, g_x_0_x_0_0_x_xx_xz, g_x_0_x_0_0_x_xx_yy, g_x_0_x_0_0_x_xx_yz, g_x_0_x_0_0_x_xx_zz, g_x_x_x_xx, g_x_x_x_xy, g_x_x_x_xz, g_x_x_x_yy, g_x_x_x_yz, g_x_x_x_zz, g_x_x_xxx_xx, g_x_x_xxx_xy, g_x_x_xxx_xz, g_x_x_xxx_yy, g_x_x_xxx_yz, g_x_x_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_x_xx_xx[i] = -4.0 * g_x_x_x_xx[i] * a_exp + 4.0 * g_x_x_xxx_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_xx_xy[i] = -4.0 * g_x_x_x_xy[i] * a_exp + 4.0 * g_x_x_xxx_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_xx_xz[i] = -4.0 * g_x_x_x_xz[i] * a_exp + 4.0 * g_x_x_xxx_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_xx_yy[i] = -4.0 * g_x_x_x_yy[i] * a_exp + 4.0 * g_x_x_xxx_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_xx_yz[i] = -4.0 * g_x_x_x_yz[i] * a_exp + 4.0 * g_x_x_xxx_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_xx_zz[i] = -4.0 * g_x_x_x_zz[i] * a_exp + 4.0 * g_x_x_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_0_x_0_0_x_xy_xx, g_x_0_x_0_0_x_xy_xy, g_x_0_x_0_0_x_xy_xz, g_x_0_x_0_0_x_xy_yy, g_x_0_x_0_0_x_xy_yz, g_x_0_x_0_0_x_xy_zz, g_x_x_xxy_xx, g_x_x_xxy_xy, g_x_x_xxy_xz, g_x_x_xxy_yy, g_x_x_xxy_yz, g_x_x_xxy_zz, g_x_x_y_xx, g_x_x_y_xy, g_x_x_y_xz, g_x_x_y_yy, g_x_x_y_yz, g_x_x_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_x_xy_xx[i] = -2.0 * g_x_x_y_xx[i] * a_exp + 4.0 * g_x_x_xxy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_xy_xy[i] = -2.0 * g_x_x_y_xy[i] * a_exp + 4.0 * g_x_x_xxy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_xy_xz[i] = -2.0 * g_x_x_y_xz[i] * a_exp + 4.0 * g_x_x_xxy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_xy_yy[i] = -2.0 * g_x_x_y_yy[i] * a_exp + 4.0 * g_x_x_xxy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_xy_yz[i] = -2.0 * g_x_x_y_yz[i] * a_exp + 4.0 * g_x_x_xxy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_xy_zz[i] = -2.0 * g_x_x_y_zz[i] * a_exp + 4.0 * g_x_x_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_0_x_0_0_x_xz_xx, g_x_0_x_0_0_x_xz_xy, g_x_0_x_0_0_x_xz_xz, g_x_0_x_0_0_x_xz_yy, g_x_0_x_0_0_x_xz_yz, g_x_0_x_0_0_x_xz_zz, g_x_x_xxz_xx, g_x_x_xxz_xy, g_x_x_xxz_xz, g_x_x_xxz_yy, g_x_x_xxz_yz, g_x_x_xxz_zz, g_x_x_z_xx, g_x_x_z_xy, g_x_x_z_xz, g_x_x_z_yy, g_x_x_z_yz, g_x_x_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_x_xz_xx[i] = -2.0 * g_x_x_z_xx[i] * a_exp + 4.0 * g_x_x_xxz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_xz_xy[i] = -2.0 * g_x_x_z_xy[i] * a_exp + 4.0 * g_x_x_xxz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_xz_xz[i] = -2.0 * g_x_x_z_xz[i] * a_exp + 4.0 * g_x_x_xxz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_xz_yy[i] = -2.0 * g_x_x_z_yy[i] * a_exp + 4.0 * g_x_x_xxz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_xz_yz[i] = -2.0 * g_x_x_z_yz[i] * a_exp + 4.0 * g_x_x_xxz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_xz_zz[i] = -2.0 * g_x_x_z_zz[i] * a_exp + 4.0 * g_x_x_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_x_0_x_0_0_x_yy_xx, g_x_0_x_0_0_x_yy_xy, g_x_0_x_0_0_x_yy_xz, g_x_0_x_0_0_x_yy_yy, g_x_0_x_0_0_x_yy_yz, g_x_0_x_0_0_x_yy_zz, g_x_x_xyy_xx, g_x_x_xyy_xy, g_x_x_xyy_xz, g_x_x_xyy_yy, g_x_x_xyy_yz, g_x_x_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_x_yy_xx[i] = 4.0 * g_x_x_xyy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_yy_xy[i] = 4.0 * g_x_x_xyy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_yy_xz[i] = 4.0 * g_x_x_xyy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_yy_yy[i] = 4.0 * g_x_x_xyy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_yy_yz[i] = 4.0 * g_x_x_xyy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_yy_zz[i] = 4.0 * g_x_x_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_x_0_x_0_0_x_yz_xx, g_x_0_x_0_0_x_yz_xy, g_x_0_x_0_0_x_yz_xz, g_x_0_x_0_0_x_yz_yy, g_x_0_x_0_0_x_yz_yz, g_x_0_x_0_0_x_yz_zz, g_x_x_xyz_xx, g_x_x_xyz_xy, g_x_x_xyz_xz, g_x_x_xyz_yy, g_x_x_xyz_yz, g_x_x_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_x_yz_xx[i] = 4.0 * g_x_x_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_yz_xy[i] = 4.0 * g_x_x_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_yz_xz[i] = 4.0 * g_x_x_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_yz_yy[i] = 4.0 * g_x_x_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_yz_yz[i] = 4.0 * g_x_x_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_yz_zz[i] = 4.0 * g_x_x_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_x_0_x_0_0_x_zz_xx, g_x_0_x_0_0_x_zz_xy, g_x_0_x_0_0_x_zz_xz, g_x_0_x_0_0_x_zz_yy, g_x_0_x_0_0_x_zz_yz, g_x_0_x_0_0_x_zz_zz, g_x_x_xzz_xx, g_x_x_xzz_xy, g_x_x_xzz_xz, g_x_x_xzz_yy, g_x_x_xzz_yz, g_x_x_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_x_zz_xx[i] = 4.0 * g_x_x_xzz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_zz_xy[i] = 4.0 * g_x_x_xzz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_zz_xz[i] = 4.0 * g_x_x_xzz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_zz_yy[i] = 4.0 * g_x_x_xzz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_zz_yz[i] = 4.0 * g_x_x_xzz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_zz_zz[i] = 4.0 * g_x_x_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_x_0_x_0_0_y_xx_xx, g_x_0_x_0_0_y_xx_xy, g_x_0_x_0_0_y_xx_xz, g_x_0_x_0_0_y_xx_yy, g_x_0_x_0_0_y_xx_yz, g_x_0_x_0_0_y_xx_zz, g_x_y_x_xx, g_x_y_x_xy, g_x_y_x_xz, g_x_y_x_yy, g_x_y_x_yz, g_x_y_x_zz, g_x_y_xxx_xx, g_x_y_xxx_xy, g_x_y_xxx_xz, g_x_y_xxx_yy, g_x_y_xxx_yz, g_x_y_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_y_xx_xx[i] = -4.0 * g_x_y_x_xx[i] * a_exp + 4.0 * g_x_y_xxx_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_xx_xy[i] = -4.0 * g_x_y_x_xy[i] * a_exp + 4.0 * g_x_y_xxx_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_xx_xz[i] = -4.0 * g_x_y_x_xz[i] * a_exp + 4.0 * g_x_y_xxx_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_xx_yy[i] = -4.0 * g_x_y_x_yy[i] * a_exp + 4.0 * g_x_y_xxx_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_xx_yz[i] = -4.0 * g_x_y_x_yz[i] * a_exp + 4.0 * g_x_y_xxx_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_xx_zz[i] = -4.0 * g_x_y_x_zz[i] * a_exp + 4.0 * g_x_y_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_x_0_x_0_0_y_xy_xx, g_x_0_x_0_0_y_xy_xy, g_x_0_x_0_0_y_xy_xz, g_x_0_x_0_0_y_xy_yy, g_x_0_x_0_0_y_xy_yz, g_x_0_x_0_0_y_xy_zz, g_x_y_xxy_xx, g_x_y_xxy_xy, g_x_y_xxy_xz, g_x_y_xxy_yy, g_x_y_xxy_yz, g_x_y_xxy_zz, g_x_y_y_xx, g_x_y_y_xy, g_x_y_y_xz, g_x_y_y_yy, g_x_y_y_yz, g_x_y_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_y_xy_xx[i] = -2.0 * g_x_y_y_xx[i] * a_exp + 4.0 * g_x_y_xxy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_xy_xy[i] = -2.0 * g_x_y_y_xy[i] * a_exp + 4.0 * g_x_y_xxy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_xy_xz[i] = -2.0 * g_x_y_y_xz[i] * a_exp + 4.0 * g_x_y_xxy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_xy_yy[i] = -2.0 * g_x_y_y_yy[i] * a_exp + 4.0 * g_x_y_xxy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_xy_yz[i] = -2.0 * g_x_y_y_yz[i] * a_exp + 4.0 * g_x_y_xxy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_xy_zz[i] = -2.0 * g_x_y_y_zz[i] * a_exp + 4.0 * g_x_y_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_x_0_x_0_0_y_xz_xx, g_x_0_x_0_0_y_xz_xy, g_x_0_x_0_0_y_xz_xz, g_x_0_x_0_0_y_xz_yy, g_x_0_x_0_0_y_xz_yz, g_x_0_x_0_0_y_xz_zz, g_x_y_xxz_xx, g_x_y_xxz_xy, g_x_y_xxz_xz, g_x_y_xxz_yy, g_x_y_xxz_yz, g_x_y_xxz_zz, g_x_y_z_xx, g_x_y_z_xy, g_x_y_z_xz, g_x_y_z_yy, g_x_y_z_yz, g_x_y_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_y_xz_xx[i] = -2.0 * g_x_y_z_xx[i] * a_exp + 4.0 * g_x_y_xxz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_xz_xy[i] = -2.0 * g_x_y_z_xy[i] * a_exp + 4.0 * g_x_y_xxz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_xz_xz[i] = -2.0 * g_x_y_z_xz[i] * a_exp + 4.0 * g_x_y_xxz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_xz_yy[i] = -2.0 * g_x_y_z_yy[i] * a_exp + 4.0 * g_x_y_xxz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_xz_yz[i] = -2.0 * g_x_y_z_yz[i] * a_exp + 4.0 * g_x_y_xxz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_xz_zz[i] = -2.0 * g_x_y_z_zz[i] * a_exp + 4.0 * g_x_y_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_x_0_x_0_0_y_yy_xx, g_x_0_x_0_0_y_yy_xy, g_x_0_x_0_0_y_yy_xz, g_x_0_x_0_0_y_yy_yy, g_x_0_x_0_0_y_yy_yz, g_x_0_x_0_0_y_yy_zz, g_x_y_xyy_xx, g_x_y_xyy_xy, g_x_y_xyy_xz, g_x_y_xyy_yy, g_x_y_xyy_yz, g_x_y_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_y_yy_xx[i] = 4.0 * g_x_y_xyy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_yy_xy[i] = 4.0 * g_x_y_xyy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_yy_xz[i] = 4.0 * g_x_y_xyy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_yy_yy[i] = 4.0 * g_x_y_xyy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_yy_yz[i] = 4.0 * g_x_y_xyy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_yy_zz[i] = 4.0 * g_x_y_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_x_0_x_0_0_y_yz_xx, g_x_0_x_0_0_y_yz_xy, g_x_0_x_0_0_y_yz_xz, g_x_0_x_0_0_y_yz_yy, g_x_0_x_0_0_y_yz_yz, g_x_0_x_0_0_y_yz_zz, g_x_y_xyz_xx, g_x_y_xyz_xy, g_x_y_xyz_xz, g_x_y_xyz_yy, g_x_y_xyz_yz, g_x_y_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_y_yz_xx[i] = 4.0 * g_x_y_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_yz_xy[i] = 4.0 * g_x_y_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_yz_xz[i] = 4.0 * g_x_y_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_yz_yy[i] = 4.0 * g_x_y_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_yz_yz[i] = 4.0 * g_x_y_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_yz_zz[i] = 4.0 * g_x_y_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_x_0_x_0_0_y_zz_xx, g_x_0_x_0_0_y_zz_xy, g_x_0_x_0_0_y_zz_xz, g_x_0_x_0_0_y_zz_yy, g_x_0_x_0_0_y_zz_yz, g_x_0_x_0_0_y_zz_zz, g_x_y_xzz_xx, g_x_y_xzz_xy, g_x_y_xzz_xz, g_x_y_xzz_yy, g_x_y_xzz_yz, g_x_y_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_y_zz_xx[i] = 4.0 * g_x_y_xzz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_zz_xy[i] = 4.0 * g_x_y_xzz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_zz_xz[i] = 4.0 * g_x_y_xzz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_zz_yy[i] = 4.0 * g_x_y_xzz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_zz_yz[i] = 4.0 * g_x_y_xzz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_zz_zz[i] = 4.0 * g_x_y_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_x_0_x_0_0_z_xx_xx, g_x_0_x_0_0_z_xx_xy, g_x_0_x_0_0_z_xx_xz, g_x_0_x_0_0_z_xx_yy, g_x_0_x_0_0_z_xx_yz, g_x_0_x_0_0_z_xx_zz, g_x_z_x_xx, g_x_z_x_xy, g_x_z_x_xz, g_x_z_x_yy, g_x_z_x_yz, g_x_z_x_zz, g_x_z_xxx_xx, g_x_z_xxx_xy, g_x_z_xxx_xz, g_x_z_xxx_yy, g_x_z_xxx_yz, g_x_z_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_z_xx_xx[i] = -4.0 * g_x_z_x_xx[i] * a_exp + 4.0 * g_x_z_xxx_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_xx_xy[i] = -4.0 * g_x_z_x_xy[i] * a_exp + 4.0 * g_x_z_xxx_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_xx_xz[i] = -4.0 * g_x_z_x_xz[i] * a_exp + 4.0 * g_x_z_xxx_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_xx_yy[i] = -4.0 * g_x_z_x_yy[i] * a_exp + 4.0 * g_x_z_xxx_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_xx_yz[i] = -4.0 * g_x_z_x_yz[i] * a_exp + 4.0 * g_x_z_xxx_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_xx_zz[i] = -4.0 * g_x_z_x_zz[i] * a_exp + 4.0 * g_x_z_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_x_0_x_0_0_z_xy_xx, g_x_0_x_0_0_z_xy_xy, g_x_0_x_0_0_z_xy_xz, g_x_0_x_0_0_z_xy_yy, g_x_0_x_0_0_z_xy_yz, g_x_0_x_0_0_z_xy_zz, g_x_z_xxy_xx, g_x_z_xxy_xy, g_x_z_xxy_xz, g_x_z_xxy_yy, g_x_z_xxy_yz, g_x_z_xxy_zz, g_x_z_y_xx, g_x_z_y_xy, g_x_z_y_xz, g_x_z_y_yy, g_x_z_y_yz, g_x_z_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_z_xy_xx[i] = -2.0 * g_x_z_y_xx[i] * a_exp + 4.0 * g_x_z_xxy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_xy_xy[i] = -2.0 * g_x_z_y_xy[i] * a_exp + 4.0 * g_x_z_xxy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_xy_xz[i] = -2.0 * g_x_z_y_xz[i] * a_exp + 4.0 * g_x_z_xxy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_xy_yy[i] = -2.0 * g_x_z_y_yy[i] * a_exp + 4.0 * g_x_z_xxy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_xy_yz[i] = -2.0 * g_x_z_y_yz[i] * a_exp + 4.0 * g_x_z_xxy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_xy_zz[i] = -2.0 * g_x_z_y_zz[i] * a_exp + 4.0 * g_x_z_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_x_0_x_0_0_z_xz_xx, g_x_0_x_0_0_z_xz_xy, g_x_0_x_0_0_z_xz_xz, g_x_0_x_0_0_z_xz_yy, g_x_0_x_0_0_z_xz_yz, g_x_0_x_0_0_z_xz_zz, g_x_z_xxz_xx, g_x_z_xxz_xy, g_x_z_xxz_xz, g_x_z_xxz_yy, g_x_z_xxz_yz, g_x_z_xxz_zz, g_x_z_z_xx, g_x_z_z_xy, g_x_z_z_xz, g_x_z_z_yy, g_x_z_z_yz, g_x_z_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_z_xz_xx[i] = -2.0 * g_x_z_z_xx[i] * a_exp + 4.0 * g_x_z_xxz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_xz_xy[i] = -2.0 * g_x_z_z_xy[i] * a_exp + 4.0 * g_x_z_xxz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_xz_xz[i] = -2.0 * g_x_z_z_xz[i] * a_exp + 4.0 * g_x_z_xxz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_xz_yy[i] = -2.0 * g_x_z_z_yy[i] * a_exp + 4.0 * g_x_z_xxz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_xz_yz[i] = -2.0 * g_x_z_z_yz[i] * a_exp + 4.0 * g_x_z_xxz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_xz_zz[i] = -2.0 * g_x_z_z_zz[i] * a_exp + 4.0 * g_x_z_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_x_0_x_0_0_z_yy_xx, g_x_0_x_0_0_z_yy_xy, g_x_0_x_0_0_z_yy_xz, g_x_0_x_0_0_z_yy_yy, g_x_0_x_0_0_z_yy_yz, g_x_0_x_0_0_z_yy_zz, g_x_z_xyy_xx, g_x_z_xyy_xy, g_x_z_xyy_xz, g_x_z_xyy_yy, g_x_z_xyy_yz, g_x_z_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_z_yy_xx[i] = 4.0 * g_x_z_xyy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_yy_xy[i] = 4.0 * g_x_z_xyy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_yy_xz[i] = 4.0 * g_x_z_xyy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_yy_yy[i] = 4.0 * g_x_z_xyy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_yy_yz[i] = 4.0 * g_x_z_xyy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_yy_zz[i] = 4.0 * g_x_z_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_x_0_x_0_0_z_yz_xx, g_x_0_x_0_0_z_yz_xy, g_x_0_x_0_0_z_yz_xz, g_x_0_x_0_0_z_yz_yy, g_x_0_x_0_0_z_yz_yz, g_x_0_x_0_0_z_yz_zz, g_x_z_xyz_xx, g_x_z_xyz_xy, g_x_z_xyz_xz, g_x_z_xyz_yy, g_x_z_xyz_yz, g_x_z_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_z_yz_xx[i] = 4.0 * g_x_z_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_yz_xy[i] = 4.0 * g_x_z_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_yz_xz[i] = 4.0 * g_x_z_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_yz_yy[i] = 4.0 * g_x_z_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_yz_yz[i] = 4.0 * g_x_z_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_yz_zz[i] = 4.0 * g_x_z_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_x_0_x_0_0_z_zz_xx, g_x_0_x_0_0_z_zz_xy, g_x_0_x_0_0_z_zz_xz, g_x_0_x_0_0_z_zz_yy, g_x_0_x_0_0_z_zz_yz, g_x_0_x_0_0_z_zz_zz, g_x_z_xzz_xx, g_x_z_xzz_xy, g_x_z_xzz_xz, g_x_z_xzz_yy, g_x_z_xzz_yz, g_x_z_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_z_zz_xx[i] = 4.0 * g_x_z_xzz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_zz_xy[i] = 4.0 * g_x_z_xzz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_zz_xz[i] = 4.0 * g_x_z_xzz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_zz_yy[i] = 4.0 * g_x_z_xzz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_zz_yz[i] = 4.0 * g_x_z_xzz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_zz_zz[i] = 4.0 * g_x_z_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_x_0_y_0_0_x_xx_xx, g_x_0_y_0_0_x_xx_xy, g_x_0_y_0_0_x_xx_xz, g_x_0_y_0_0_x_xx_yy, g_x_0_y_0_0_x_xx_yz, g_x_0_y_0_0_x_xx_zz, g_x_x_xxy_xx, g_x_x_xxy_xy, g_x_x_xxy_xz, g_x_x_xxy_yy, g_x_x_xxy_yz, g_x_x_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_x_xx_xx[i] = 4.0 * g_x_x_xxy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_xx_xy[i] = 4.0 * g_x_x_xxy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_xx_xz[i] = 4.0 * g_x_x_xxy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_xx_yy[i] = 4.0 * g_x_x_xxy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_xx_yz[i] = 4.0 * g_x_x_xxy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_xx_zz[i] = 4.0 * g_x_x_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_x_0_y_0_0_x_xy_xx, g_x_0_y_0_0_x_xy_xy, g_x_0_y_0_0_x_xy_xz, g_x_0_y_0_0_x_xy_yy, g_x_0_y_0_0_x_xy_yz, g_x_0_y_0_0_x_xy_zz, g_x_x_x_xx, g_x_x_x_xy, g_x_x_x_xz, g_x_x_x_yy, g_x_x_x_yz, g_x_x_x_zz, g_x_x_xyy_xx, g_x_x_xyy_xy, g_x_x_xyy_xz, g_x_x_xyy_yy, g_x_x_xyy_yz, g_x_x_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_x_xy_xx[i] = -2.0 * g_x_x_x_xx[i] * a_exp + 4.0 * g_x_x_xyy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_xy_xy[i] = -2.0 * g_x_x_x_xy[i] * a_exp + 4.0 * g_x_x_xyy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_xy_xz[i] = -2.0 * g_x_x_x_xz[i] * a_exp + 4.0 * g_x_x_xyy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_xy_yy[i] = -2.0 * g_x_x_x_yy[i] * a_exp + 4.0 * g_x_x_xyy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_xy_yz[i] = -2.0 * g_x_x_x_yz[i] * a_exp + 4.0 * g_x_x_xyy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_xy_zz[i] = -2.0 * g_x_x_x_zz[i] * a_exp + 4.0 * g_x_x_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_x_0_y_0_0_x_xz_xx, g_x_0_y_0_0_x_xz_xy, g_x_0_y_0_0_x_xz_xz, g_x_0_y_0_0_x_xz_yy, g_x_0_y_0_0_x_xz_yz, g_x_0_y_0_0_x_xz_zz, g_x_x_xyz_xx, g_x_x_xyz_xy, g_x_x_xyz_xz, g_x_x_xyz_yy, g_x_x_xyz_yz, g_x_x_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_x_xz_xx[i] = 4.0 * g_x_x_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_xz_xy[i] = 4.0 * g_x_x_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_xz_xz[i] = 4.0 * g_x_x_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_xz_yy[i] = 4.0 * g_x_x_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_xz_yz[i] = 4.0 * g_x_x_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_xz_zz[i] = 4.0 * g_x_x_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_x_0_y_0_0_x_yy_xx, g_x_0_y_0_0_x_yy_xy, g_x_0_y_0_0_x_yy_xz, g_x_0_y_0_0_x_yy_yy, g_x_0_y_0_0_x_yy_yz, g_x_0_y_0_0_x_yy_zz, g_x_x_y_xx, g_x_x_y_xy, g_x_x_y_xz, g_x_x_y_yy, g_x_x_y_yz, g_x_x_y_zz, g_x_x_yyy_xx, g_x_x_yyy_xy, g_x_x_yyy_xz, g_x_x_yyy_yy, g_x_x_yyy_yz, g_x_x_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_x_yy_xx[i] = -4.0 * g_x_x_y_xx[i] * a_exp + 4.0 * g_x_x_yyy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_yy_xy[i] = -4.0 * g_x_x_y_xy[i] * a_exp + 4.0 * g_x_x_yyy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_yy_xz[i] = -4.0 * g_x_x_y_xz[i] * a_exp + 4.0 * g_x_x_yyy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_yy_yy[i] = -4.0 * g_x_x_y_yy[i] * a_exp + 4.0 * g_x_x_yyy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_yy_yz[i] = -4.0 * g_x_x_y_yz[i] * a_exp + 4.0 * g_x_x_yyy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_yy_zz[i] = -4.0 * g_x_x_y_zz[i] * a_exp + 4.0 * g_x_x_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_x_0_y_0_0_x_yz_xx, g_x_0_y_0_0_x_yz_xy, g_x_0_y_0_0_x_yz_xz, g_x_0_y_0_0_x_yz_yy, g_x_0_y_0_0_x_yz_yz, g_x_0_y_0_0_x_yz_zz, g_x_x_yyz_xx, g_x_x_yyz_xy, g_x_x_yyz_xz, g_x_x_yyz_yy, g_x_x_yyz_yz, g_x_x_yyz_zz, g_x_x_z_xx, g_x_x_z_xy, g_x_x_z_xz, g_x_x_z_yy, g_x_x_z_yz, g_x_x_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_x_yz_xx[i] = -2.0 * g_x_x_z_xx[i] * a_exp + 4.0 * g_x_x_yyz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_yz_xy[i] = -2.0 * g_x_x_z_xy[i] * a_exp + 4.0 * g_x_x_yyz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_yz_xz[i] = -2.0 * g_x_x_z_xz[i] * a_exp + 4.0 * g_x_x_yyz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_yz_yy[i] = -2.0 * g_x_x_z_yy[i] * a_exp + 4.0 * g_x_x_yyz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_yz_yz[i] = -2.0 * g_x_x_z_yz[i] * a_exp + 4.0 * g_x_x_yyz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_yz_zz[i] = -2.0 * g_x_x_z_zz[i] * a_exp + 4.0 * g_x_x_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_x_0_y_0_0_x_zz_xx, g_x_0_y_0_0_x_zz_xy, g_x_0_y_0_0_x_zz_xz, g_x_0_y_0_0_x_zz_yy, g_x_0_y_0_0_x_zz_yz, g_x_0_y_0_0_x_zz_zz, g_x_x_yzz_xx, g_x_x_yzz_xy, g_x_x_yzz_xz, g_x_x_yzz_yy, g_x_x_yzz_yz, g_x_x_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_x_zz_xx[i] = 4.0 * g_x_x_yzz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_zz_xy[i] = 4.0 * g_x_x_yzz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_zz_xz[i] = 4.0 * g_x_x_yzz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_zz_yy[i] = 4.0 * g_x_x_yzz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_zz_yz[i] = 4.0 * g_x_x_yzz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_zz_zz[i] = 4.0 * g_x_x_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_x_0_y_0_0_y_xx_xx, g_x_0_y_0_0_y_xx_xy, g_x_0_y_0_0_y_xx_xz, g_x_0_y_0_0_y_xx_yy, g_x_0_y_0_0_y_xx_yz, g_x_0_y_0_0_y_xx_zz, g_x_y_xxy_xx, g_x_y_xxy_xy, g_x_y_xxy_xz, g_x_y_xxy_yy, g_x_y_xxy_yz, g_x_y_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_y_xx_xx[i] = 4.0 * g_x_y_xxy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_xx_xy[i] = 4.0 * g_x_y_xxy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_xx_xz[i] = 4.0 * g_x_y_xxy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_xx_yy[i] = 4.0 * g_x_y_xxy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_xx_yz[i] = 4.0 * g_x_y_xxy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_xx_zz[i] = 4.0 * g_x_y_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_x_0_y_0_0_y_xy_xx, g_x_0_y_0_0_y_xy_xy, g_x_0_y_0_0_y_xy_xz, g_x_0_y_0_0_y_xy_yy, g_x_0_y_0_0_y_xy_yz, g_x_0_y_0_0_y_xy_zz, g_x_y_x_xx, g_x_y_x_xy, g_x_y_x_xz, g_x_y_x_yy, g_x_y_x_yz, g_x_y_x_zz, g_x_y_xyy_xx, g_x_y_xyy_xy, g_x_y_xyy_xz, g_x_y_xyy_yy, g_x_y_xyy_yz, g_x_y_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_y_xy_xx[i] = -2.0 * g_x_y_x_xx[i] * a_exp + 4.0 * g_x_y_xyy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_xy_xy[i] = -2.0 * g_x_y_x_xy[i] * a_exp + 4.0 * g_x_y_xyy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_xy_xz[i] = -2.0 * g_x_y_x_xz[i] * a_exp + 4.0 * g_x_y_xyy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_xy_yy[i] = -2.0 * g_x_y_x_yy[i] * a_exp + 4.0 * g_x_y_xyy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_xy_yz[i] = -2.0 * g_x_y_x_yz[i] * a_exp + 4.0 * g_x_y_xyy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_xy_zz[i] = -2.0 * g_x_y_x_zz[i] * a_exp + 4.0 * g_x_y_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_x_0_y_0_0_y_xz_xx, g_x_0_y_0_0_y_xz_xy, g_x_0_y_0_0_y_xz_xz, g_x_0_y_0_0_y_xz_yy, g_x_0_y_0_0_y_xz_yz, g_x_0_y_0_0_y_xz_zz, g_x_y_xyz_xx, g_x_y_xyz_xy, g_x_y_xyz_xz, g_x_y_xyz_yy, g_x_y_xyz_yz, g_x_y_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_y_xz_xx[i] = 4.0 * g_x_y_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_xz_xy[i] = 4.0 * g_x_y_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_xz_xz[i] = 4.0 * g_x_y_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_xz_yy[i] = 4.0 * g_x_y_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_xz_yz[i] = 4.0 * g_x_y_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_xz_zz[i] = 4.0 * g_x_y_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_x_0_y_0_0_y_yy_xx, g_x_0_y_0_0_y_yy_xy, g_x_0_y_0_0_y_yy_xz, g_x_0_y_0_0_y_yy_yy, g_x_0_y_0_0_y_yy_yz, g_x_0_y_0_0_y_yy_zz, g_x_y_y_xx, g_x_y_y_xy, g_x_y_y_xz, g_x_y_y_yy, g_x_y_y_yz, g_x_y_y_zz, g_x_y_yyy_xx, g_x_y_yyy_xy, g_x_y_yyy_xz, g_x_y_yyy_yy, g_x_y_yyy_yz, g_x_y_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_y_yy_xx[i] = -4.0 * g_x_y_y_xx[i] * a_exp + 4.0 * g_x_y_yyy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_yy_xy[i] = -4.0 * g_x_y_y_xy[i] * a_exp + 4.0 * g_x_y_yyy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_yy_xz[i] = -4.0 * g_x_y_y_xz[i] * a_exp + 4.0 * g_x_y_yyy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_yy_yy[i] = -4.0 * g_x_y_y_yy[i] * a_exp + 4.0 * g_x_y_yyy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_yy_yz[i] = -4.0 * g_x_y_y_yz[i] * a_exp + 4.0 * g_x_y_yyy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_yy_zz[i] = -4.0 * g_x_y_y_zz[i] * a_exp + 4.0 * g_x_y_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_x_0_y_0_0_y_yz_xx, g_x_0_y_0_0_y_yz_xy, g_x_0_y_0_0_y_yz_xz, g_x_0_y_0_0_y_yz_yy, g_x_0_y_0_0_y_yz_yz, g_x_0_y_0_0_y_yz_zz, g_x_y_yyz_xx, g_x_y_yyz_xy, g_x_y_yyz_xz, g_x_y_yyz_yy, g_x_y_yyz_yz, g_x_y_yyz_zz, g_x_y_z_xx, g_x_y_z_xy, g_x_y_z_xz, g_x_y_z_yy, g_x_y_z_yz, g_x_y_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_y_yz_xx[i] = -2.0 * g_x_y_z_xx[i] * a_exp + 4.0 * g_x_y_yyz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_yz_xy[i] = -2.0 * g_x_y_z_xy[i] * a_exp + 4.0 * g_x_y_yyz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_yz_xz[i] = -2.0 * g_x_y_z_xz[i] * a_exp + 4.0 * g_x_y_yyz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_yz_yy[i] = -2.0 * g_x_y_z_yy[i] * a_exp + 4.0 * g_x_y_yyz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_yz_yz[i] = -2.0 * g_x_y_z_yz[i] * a_exp + 4.0 * g_x_y_yyz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_yz_zz[i] = -2.0 * g_x_y_z_zz[i] * a_exp + 4.0 * g_x_y_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_x_0_y_0_0_y_zz_xx, g_x_0_y_0_0_y_zz_xy, g_x_0_y_0_0_y_zz_xz, g_x_0_y_0_0_y_zz_yy, g_x_0_y_0_0_y_zz_yz, g_x_0_y_0_0_y_zz_zz, g_x_y_yzz_xx, g_x_y_yzz_xy, g_x_y_yzz_xz, g_x_y_yzz_yy, g_x_y_yzz_yz, g_x_y_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_y_zz_xx[i] = 4.0 * g_x_y_yzz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_zz_xy[i] = 4.0 * g_x_y_yzz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_zz_xz[i] = 4.0 * g_x_y_yzz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_zz_yy[i] = 4.0 * g_x_y_yzz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_zz_yz[i] = 4.0 * g_x_y_yzz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_zz_zz[i] = 4.0 * g_x_y_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_x_0_y_0_0_z_xx_xx, g_x_0_y_0_0_z_xx_xy, g_x_0_y_0_0_z_xx_xz, g_x_0_y_0_0_z_xx_yy, g_x_0_y_0_0_z_xx_yz, g_x_0_y_0_0_z_xx_zz, g_x_z_xxy_xx, g_x_z_xxy_xy, g_x_z_xxy_xz, g_x_z_xxy_yy, g_x_z_xxy_yz, g_x_z_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_z_xx_xx[i] = 4.0 * g_x_z_xxy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_xx_xy[i] = 4.0 * g_x_z_xxy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_xx_xz[i] = 4.0 * g_x_z_xxy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_xx_yy[i] = 4.0 * g_x_z_xxy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_xx_yz[i] = 4.0 * g_x_z_xxy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_xx_zz[i] = 4.0 * g_x_z_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_x_0_y_0_0_z_xy_xx, g_x_0_y_0_0_z_xy_xy, g_x_0_y_0_0_z_xy_xz, g_x_0_y_0_0_z_xy_yy, g_x_0_y_0_0_z_xy_yz, g_x_0_y_0_0_z_xy_zz, g_x_z_x_xx, g_x_z_x_xy, g_x_z_x_xz, g_x_z_x_yy, g_x_z_x_yz, g_x_z_x_zz, g_x_z_xyy_xx, g_x_z_xyy_xy, g_x_z_xyy_xz, g_x_z_xyy_yy, g_x_z_xyy_yz, g_x_z_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_z_xy_xx[i] = -2.0 * g_x_z_x_xx[i] * a_exp + 4.0 * g_x_z_xyy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_xy_xy[i] = -2.0 * g_x_z_x_xy[i] * a_exp + 4.0 * g_x_z_xyy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_xy_xz[i] = -2.0 * g_x_z_x_xz[i] * a_exp + 4.0 * g_x_z_xyy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_xy_yy[i] = -2.0 * g_x_z_x_yy[i] * a_exp + 4.0 * g_x_z_xyy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_xy_yz[i] = -2.0 * g_x_z_x_yz[i] * a_exp + 4.0 * g_x_z_xyy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_xy_zz[i] = -2.0 * g_x_z_x_zz[i] * a_exp + 4.0 * g_x_z_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_x_0_y_0_0_z_xz_xx, g_x_0_y_0_0_z_xz_xy, g_x_0_y_0_0_z_xz_xz, g_x_0_y_0_0_z_xz_yy, g_x_0_y_0_0_z_xz_yz, g_x_0_y_0_0_z_xz_zz, g_x_z_xyz_xx, g_x_z_xyz_xy, g_x_z_xyz_xz, g_x_z_xyz_yy, g_x_z_xyz_yz, g_x_z_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_z_xz_xx[i] = 4.0 * g_x_z_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_xz_xy[i] = 4.0 * g_x_z_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_xz_xz[i] = 4.0 * g_x_z_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_xz_yy[i] = 4.0 * g_x_z_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_xz_yz[i] = 4.0 * g_x_z_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_xz_zz[i] = 4.0 * g_x_z_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_x_0_y_0_0_z_yy_xx, g_x_0_y_0_0_z_yy_xy, g_x_0_y_0_0_z_yy_xz, g_x_0_y_0_0_z_yy_yy, g_x_0_y_0_0_z_yy_yz, g_x_0_y_0_0_z_yy_zz, g_x_z_y_xx, g_x_z_y_xy, g_x_z_y_xz, g_x_z_y_yy, g_x_z_y_yz, g_x_z_y_zz, g_x_z_yyy_xx, g_x_z_yyy_xy, g_x_z_yyy_xz, g_x_z_yyy_yy, g_x_z_yyy_yz, g_x_z_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_z_yy_xx[i] = -4.0 * g_x_z_y_xx[i] * a_exp + 4.0 * g_x_z_yyy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_yy_xy[i] = -4.0 * g_x_z_y_xy[i] * a_exp + 4.0 * g_x_z_yyy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_yy_xz[i] = -4.0 * g_x_z_y_xz[i] * a_exp + 4.0 * g_x_z_yyy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_yy_yy[i] = -4.0 * g_x_z_y_yy[i] * a_exp + 4.0 * g_x_z_yyy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_yy_yz[i] = -4.0 * g_x_z_y_yz[i] * a_exp + 4.0 * g_x_z_yyy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_yy_zz[i] = -4.0 * g_x_z_y_zz[i] * a_exp + 4.0 * g_x_z_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_x_0_y_0_0_z_yz_xx, g_x_0_y_0_0_z_yz_xy, g_x_0_y_0_0_z_yz_xz, g_x_0_y_0_0_z_yz_yy, g_x_0_y_0_0_z_yz_yz, g_x_0_y_0_0_z_yz_zz, g_x_z_yyz_xx, g_x_z_yyz_xy, g_x_z_yyz_xz, g_x_z_yyz_yy, g_x_z_yyz_yz, g_x_z_yyz_zz, g_x_z_z_xx, g_x_z_z_xy, g_x_z_z_xz, g_x_z_z_yy, g_x_z_z_yz, g_x_z_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_z_yz_xx[i] = -2.0 * g_x_z_z_xx[i] * a_exp + 4.0 * g_x_z_yyz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_yz_xy[i] = -2.0 * g_x_z_z_xy[i] * a_exp + 4.0 * g_x_z_yyz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_yz_xz[i] = -2.0 * g_x_z_z_xz[i] * a_exp + 4.0 * g_x_z_yyz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_yz_yy[i] = -2.0 * g_x_z_z_yy[i] * a_exp + 4.0 * g_x_z_yyz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_yz_yz[i] = -2.0 * g_x_z_z_yz[i] * a_exp + 4.0 * g_x_z_yyz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_yz_zz[i] = -2.0 * g_x_z_z_zz[i] * a_exp + 4.0 * g_x_z_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_x_0_y_0_0_z_zz_xx, g_x_0_y_0_0_z_zz_xy, g_x_0_y_0_0_z_zz_xz, g_x_0_y_0_0_z_zz_yy, g_x_0_y_0_0_z_zz_yz, g_x_0_y_0_0_z_zz_zz, g_x_z_yzz_xx, g_x_z_yzz_xy, g_x_z_yzz_xz, g_x_z_yzz_yy, g_x_z_yzz_yz, g_x_z_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_z_zz_xx[i] = 4.0 * g_x_z_yzz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_zz_xy[i] = 4.0 * g_x_z_yzz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_zz_xz[i] = 4.0 * g_x_z_yzz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_zz_yy[i] = 4.0 * g_x_z_yzz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_zz_yz[i] = 4.0 * g_x_z_yzz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_zz_zz[i] = 4.0 * g_x_z_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_x_0_z_0_0_x_xx_xx, g_x_0_z_0_0_x_xx_xy, g_x_0_z_0_0_x_xx_xz, g_x_0_z_0_0_x_xx_yy, g_x_0_z_0_0_x_xx_yz, g_x_0_z_0_0_x_xx_zz, g_x_x_xxz_xx, g_x_x_xxz_xy, g_x_x_xxz_xz, g_x_x_xxz_yy, g_x_x_xxz_yz, g_x_x_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_x_xx_xx[i] = 4.0 * g_x_x_xxz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_xx_xy[i] = 4.0 * g_x_x_xxz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_xx_xz[i] = 4.0 * g_x_x_xxz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_xx_yy[i] = 4.0 * g_x_x_xxz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_xx_yz[i] = 4.0 * g_x_x_xxz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_xx_zz[i] = 4.0 * g_x_x_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_x_0_z_0_0_x_xy_xx, g_x_0_z_0_0_x_xy_xy, g_x_0_z_0_0_x_xy_xz, g_x_0_z_0_0_x_xy_yy, g_x_0_z_0_0_x_xy_yz, g_x_0_z_0_0_x_xy_zz, g_x_x_xyz_xx, g_x_x_xyz_xy, g_x_x_xyz_xz, g_x_x_xyz_yy, g_x_x_xyz_yz, g_x_x_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_x_xy_xx[i] = 4.0 * g_x_x_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_xy_xy[i] = 4.0 * g_x_x_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_xy_xz[i] = 4.0 * g_x_x_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_xy_yy[i] = 4.0 * g_x_x_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_xy_yz[i] = 4.0 * g_x_x_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_xy_zz[i] = 4.0 * g_x_x_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_x_0_z_0_0_x_xz_xx, g_x_0_z_0_0_x_xz_xy, g_x_0_z_0_0_x_xz_xz, g_x_0_z_0_0_x_xz_yy, g_x_0_z_0_0_x_xz_yz, g_x_0_z_0_0_x_xz_zz, g_x_x_x_xx, g_x_x_x_xy, g_x_x_x_xz, g_x_x_x_yy, g_x_x_x_yz, g_x_x_x_zz, g_x_x_xzz_xx, g_x_x_xzz_xy, g_x_x_xzz_xz, g_x_x_xzz_yy, g_x_x_xzz_yz, g_x_x_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_x_xz_xx[i] = -2.0 * g_x_x_x_xx[i] * a_exp + 4.0 * g_x_x_xzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_xz_xy[i] = -2.0 * g_x_x_x_xy[i] * a_exp + 4.0 * g_x_x_xzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_xz_xz[i] = -2.0 * g_x_x_x_xz[i] * a_exp + 4.0 * g_x_x_xzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_xz_yy[i] = -2.0 * g_x_x_x_yy[i] * a_exp + 4.0 * g_x_x_xzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_xz_yz[i] = -2.0 * g_x_x_x_yz[i] * a_exp + 4.0 * g_x_x_xzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_xz_zz[i] = -2.0 * g_x_x_x_zz[i] * a_exp + 4.0 * g_x_x_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_x_0_z_0_0_x_yy_xx, g_x_0_z_0_0_x_yy_xy, g_x_0_z_0_0_x_yy_xz, g_x_0_z_0_0_x_yy_yy, g_x_0_z_0_0_x_yy_yz, g_x_0_z_0_0_x_yy_zz, g_x_x_yyz_xx, g_x_x_yyz_xy, g_x_x_yyz_xz, g_x_x_yyz_yy, g_x_x_yyz_yz, g_x_x_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_x_yy_xx[i] = 4.0 * g_x_x_yyz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_yy_xy[i] = 4.0 * g_x_x_yyz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_yy_xz[i] = 4.0 * g_x_x_yyz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_yy_yy[i] = 4.0 * g_x_x_yyz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_yy_yz[i] = 4.0 * g_x_x_yyz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_yy_zz[i] = 4.0 * g_x_x_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_x_0_z_0_0_x_yz_xx, g_x_0_z_0_0_x_yz_xy, g_x_0_z_0_0_x_yz_xz, g_x_0_z_0_0_x_yz_yy, g_x_0_z_0_0_x_yz_yz, g_x_0_z_0_0_x_yz_zz, g_x_x_y_xx, g_x_x_y_xy, g_x_x_y_xz, g_x_x_y_yy, g_x_x_y_yz, g_x_x_y_zz, g_x_x_yzz_xx, g_x_x_yzz_xy, g_x_x_yzz_xz, g_x_x_yzz_yy, g_x_x_yzz_yz, g_x_x_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_x_yz_xx[i] = -2.0 * g_x_x_y_xx[i] * a_exp + 4.0 * g_x_x_yzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_yz_xy[i] = -2.0 * g_x_x_y_xy[i] * a_exp + 4.0 * g_x_x_yzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_yz_xz[i] = -2.0 * g_x_x_y_xz[i] * a_exp + 4.0 * g_x_x_yzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_yz_yy[i] = -2.0 * g_x_x_y_yy[i] * a_exp + 4.0 * g_x_x_yzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_yz_yz[i] = -2.0 * g_x_x_y_yz[i] * a_exp + 4.0 * g_x_x_yzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_yz_zz[i] = -2.0 * g_x_x_y_zz[i] * a_exp + 4.0 * g_x_x_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_x_0_z_0_0_x_zz_xx, g_x_0_z_0_0_x_zz_xy, g_x_0_z_0_0_x_zz_xz, g_x_0_z_0_0_x_zz_yy, g_x_0_z_0_0_x_zz_yz, g_x_0_z_0_0_x_zz_zz, g_x_x_z_xx, g_x_x_z_xy, g_x_x_z_xz, g_x_x_z_yy, g_x_x_z_yz, g_x_x_z_zz, g_x_x_zzz_xx, g_x_x_zzz_xy, g_x_x_zzz_xz, g_x_x_zzz_yy, g_x_x_zzz_yz, g_x_x_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_x_zz_xx[i] = -4.0 * g_x_x_z_xx[i] * a_exp + 4.0 * g_x_x_zzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_zz_xy[i] = -4.0 * g_x_x_z_xy[i] * a_exp + 4.0 * g_x_x_zzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_zz_xz[i] = -4.0 * g_x_x_z_xz[i] * a_exp + 4.0 * g_x_x_zzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_zz_yy[i] = -4.0 * g_x_x_z_yy[i] * a_exp + 4.0 * g_x_x_zzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_zz_yz[i] = -4.0 * g_x_x_z_yz[i] * a_exp + 4.0 * g_x_x_zzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_zz_zz[i] = -4.0 * g_x_x_z_zz[i] * a_exp + 4.0 * g_x_x_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_x_0_z_0_0_y_xx_xx, g_x_0_z_0_0_y_xx_xy, g_x_0_z_0_0_y_xx_xz, g_x_0_z_0_0_y_xx_yy, g_x_0_z_0_0_y_xx_yz, g_x_0_z_0_0_y_xx_zz, g_x_y_xxz_xx, g_x_y_xxz_xy, g_x_y_xxz_xz, g_x_y_xxz_yy, g_x_y_xxz_yz, g_x_y_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_y_xx_xx[i] = 4.0 * g_x_y_xxz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_xx_xy[i] = 4.0 * g_x_y_xxz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_xx_xz[i] = 4.0 * g_x_y_xxz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_xx_yy[i] = 4.0 * g_x_y_xxz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_xx_yz[i] = 4.0 * g_x_y_xxz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_xx_zz[i] = 4.0 * g_x_y_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_x_0_z_0_0_y_xy_xx, g_x_0_z_0_0_y_xy_xy, g_x_0_z_0_0_y_xy_xz, g_x_0_z_0_0_y_xy_yy, g_x_0_z_0_0_y_xy_yz, g_x_0_z_0_0_y_xy_zz, g_x_y_xyz_xx, g_x_y_xyz_xy, g_x_y_xyz_xz, g_x_y_xyz_yy, g_x_y_xyz_yz, g_x_y_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_y_xy_xx[i] = 4.0 * g_x_y_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_xy_xy[i] = 4.0 * g_x_y_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_xy_xz[i] = 4.0 * g_x_y_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_xy_yy[i] = 4.0 * g_x_y_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_xy_yz[i] = 4.0 * g_x_y_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_xy_zz[i] = 4.0 * g_x_y_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_x_0_z_0_0_y_xz_xx, g_x_0_z_0_0_y_xz_xy, g_x_0_z_0_0_y_xz_xz, g_x_0_z_0_0_y_xz_yy, g_x_0_z_0_0_y_xz_yz, g_x_0_z_0_0_y_xz_zz, g_x_y_x_xx, g_x_y_x_xy, g_x_y_x_xz, g_x_y_x_yy, g_x_y_x_yz, g_x_y_x_zz, g_x_y_xzz_xx, g_x_y_xzz_xy, g_x_y_xzz_xz, g_x_y_xzz_yy, g_x_y_xzz_yz, g_x_y_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_y_xz_xx[i] = -2.0 * g_x_y_x_xx[i] * a_exp + 4.0 * g_x_y_xzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_xz_xy[i] = -2.0 * g_x_y_x_xy[i] * a_exp + 4.0 * g_x_y_xzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_xz_xz[i] = -2.0 * g_x_y_x_xz[i] * a_exp + 4.0 * g_x_y_xzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_xz_yy[i] = -2.0 * g_x_y_x_yy[i] * a_exp + 4.0 * g_x_y_xzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_xz_yz[i] = -2.0 * g_x_y_x_yz[i] * a_exp + 4.0 * g_x_y_xzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_xz_zz[i] = -2.0 * g_x_y_x_zz[i] * a_exp + 4.0 * g_x_y_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_x_0_z_0_0_y_yy_xx, g_x_0_z_0_0_y_yy_xy, g_x_0_z_0_0_y_yy_xz, g_x_0_z_0_0_y_yy_yy, g_x_0_z_0_0_y_yy_yz, g_x_0_z_0_0_y_yy_zz, g_x_y_yyz_xx, g_x_y_yyz_xy, g_x_y_yyz_xz, g_x_y_yyz_yy, g_x_y_yyz_yz, g_x_y_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_y_yy_xx[i] = 4.0 * g_x_y_yyz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_yy_xy[i] = 4.0 * g_x_y_yyz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_yy_xz[i] = 4.0 * g_x_y_yyz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_yy_yy[i] = 4.0 * g_x_y_yyz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_yy_yz[i] = 4.0 * g_x_y_yyz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_yy_zz[i] = 4.0 * g_x_y_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_x_0_z_0_0_y_yz_xx, g_x_0_z_0_0_y_yz_xy, g_x_0_z_0_0_y_yz_xz, g_x_0_z_0_0_y_yz_yy, g_x_0_z_0_0_y_yz_yz, g_x_0_z_0_0_y_yz_zz, g_x_y_y_xx, g_x_y_y_xy, g_x_y_y_xz, g_x_y_y_yy, g_x_y_y_yz, g_x_y_y_zz, g_x_y_yzz_xx, g_x_y_yzz_xy, g_x_y_yzz_xz, g_x_y_yzz_yy, g_x_y_yzz_yz, g_x_y_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_y_yz_xx[i] = -2.0 * g_x_y_y_xx[i] * a_exp + 4.0 * g_x_y_yzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_yz_xy[i] = -2.0 * g_x_y_y_xy[i] * a_exp + 4.0 * g_x_y_yzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_yz_xz[i] = -2.0 * g_x_y_y_xz[i] * a_exp + 4.0 * g_x_y_yzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_yz_yy[i] = -2.0 * g_x_y_y_yy[i] * a_exp + 4.0 * g_x_y_yzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_yz_yz[i] = -2.0 * g_x_y_y_yz[i] * a_exp + 4.0 * g_x_y_yzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_yz_zz[i] = -2.0 * g_x_y_y_zz[i] * a_exp + 4.0 * g_x_y_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_x_0_z_0_0_y_zz_xx, g_x_0_z_0_0_y_zz_xy, g_x_0_z_0_0_y_zz_xz, g_x_0_z_0_0_y_zz_yy, g_x_0_z_0_0_y_zz_yz, g_x_0_z_0_0_y_zz_zz, g_x_y_z_xx, g_x_y_z_xy, g_x_y_z_xz, g_x_y_z_yy, g_x_y_z_yz, g_x_y_z_zz, g_x_y_zzz_xx, g_x_y_zzz_xy, g_x_y_zzz_xz, g_x_y_zzz_yy, g_x_y_zzz_yz, g_x_y_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_y_zz_xx[i] = -4.0 * g_x_y_z_xx[i] * a_exp + 4.0 * g_x_y_zzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_zz_xy[i] = -4.0 * g_x_y_z_xy[i] * a_exp + 4.0 * g_x_y_zzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_zz_xz[i] = -4.0 * g_x_y_z_xz[i] * a_exp + 4.0 * g_x_y_zzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_zz_yy[i] = -4.0 * g_x_y_z_yy[i] * a_exp + 4.0 * g_x_y_zzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_zz_yz[i] = -4.0 * g_x_y_z_yz[i] * a_exp + 4.0 * g_x_y_zzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_zz_zz[i] = -4.0 * g_x_y_z_zz[i] * a_exp + 4.0 * g_x_y_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_x_0_z_0_0_z_xx_xx, g_x_0_z_0_0_z_xx_xy, g_x_0_z_0_0_z_xx_xz, g_x_0_z_0_0_z_xx_yy, g_x_0_z_0_0_z_xx_yz, g_x_0_z_0_0_z_xx_zz, g_x_z_xxz_xx, g_x_z_xxz_xy, g_x_z_xxz_xz, g_x_z_xxz_yy, g_x_z_xxz_yz, g_x_z_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_z_xx_xx[i] = 4.0 * g_x_z_xxz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_xx_xy[i] = 4.0 * g_x_z_xxz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_xx_xz[i] = 4.0 * g_x_z_xxz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_xx_yy[i] = 4.0 * g_x_z_xxz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_xx_yz[i] = 4.0 * g_x_z_xxz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_xx_zz[i] = 4.0 * g_x_z_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_x_0_z_0_0_z_xy_xx, g_x_0_z_0_0_z_xy_xy, g_x_0_z_0_0_z_xy_xz, g_x_0_z_0_0_z_xy_yy, g_x_0_z_0_0_z_xy_yz, g_x_0_z_0_0_z_xy_zz, g_x_z_xyz_xx, g_x_z_xyz_xy, g_x_z_xyz_xz, g_x_z_xyz_yy, g_x_z_xyz_yz, g_x_z_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_z_xy_xx[i] = 4.0 * g_x_z_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_xy_xy[i] = 4.0 * g_x_z_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_xy_xz[i] = 4.0 * g_x_z_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_xy_yy[i] = 4.0 * g_x_z_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_xy_yz[i] = 4.0 * g_x_z_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_xy_zz[i] = 4.0 * g_x_z_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_x_0_z_0_0_z_xz_xx, g_x_0_z_0_0_z_xz_xy, g_x_0_z_0_0_z_xz_xz, g_x_0_z_0_0_z_xz_yy, g_x_0_z_0_0_z_xz_yz, g_x_0_z_0_0_z_xz_zz, g_x_z_x_xx, g_x_z_x_xy, g_x_z_x_xz, g_x_z_x_yy, g_x_z_x_yz, g_x_z_x_zz, g_x_z_xzz_xx, g_x_z_xzz_xy, g_x_z_xzz_xz, g_x_z_xzz_yy, g_x_z_xzz_yz, g_x_z_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_z_xz_xx[i] = -2.0 * g_x_z_x_xx[i] * a_exp + 4.0 * g_x_z_xzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_xz_xy[i] = -2.0 * g_x_z_x_xy[i] * a_exp + 4.0 * g_x_z_xzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_xz_xz[i] = -2.0 * g_x_z_x_xz[i] * a_exp + 4.0 * g_x_z_xzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_xz_yy[i] = -2.0 * g_x_z_x_yy[i] * a_exp + 4.0 * g_x_z_xzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_xz_yz[i] = -2.0 * g_x_z_x_yz[i] * a_exp + 4.0 * g_x_z_xzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_xz_zz[i] = -2.0 * g_x_z_x_zz[i] * a_exp + 4.0 * g_x_z_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_x_0_z_0_0_z_yy_xx, g_x_0_z_0_0_z_yy_xy, g_x_0_z_0_0_z_yy_xz, g_x_0_z_0_0_z_yy_yy, g_x_0_z_0_0_z_yy_yz, g_x_0_z_0_0_z_yy_zz, g_x_z_yyz_xx, g_x_z_yyz_xy, g_x_z_yyz_xz, g_x_z_yyz_yy, g_x_z_yyz_yz, g_x_z_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_z_yy_xx[i] = 4.0 * g_x_z_yyz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_yy_xy[i] = 4.0 * g_x_z_yyz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_yy_xz[i] = 4.0 * g_x_z_yyz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_yy_yy[i] = 4.0 * g_x_z_yyz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_yy_yz[i] = 4.0 * g_x_z_yyz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_yy_zz[i] = 4.0 * g_x_z_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_x_0_z_0_0_z_yz_xx, g_x_0_z_0_0_z_yz_xy, g_x_0_z_0_0_z_yz_xz, g_x_0_z_0_0_z_yz_yy, g_x_0_z_0_0_z_yz_yz, g_x_0_z_0_0_z_yz_zz, g_x_z_y_xx, g_x_z_y_xy, g_x_z_y_xz, g_x_z_y_yy, g_x_z_y_yz, g_x_z_y_zz, g_x_z_yzz_xx, g_x_z_yzz_xy, g_x_z_yzz_xz, g_x_z_yzz_yy, g_x_z_yzz_yz, g_x_z_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_z_yz_xx[i] = -2.0 * g_x_z_y_xx[i] * a_exp + 4.0 * g_x_z_yzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_yz_xy[i] = -2.0 * g_x_z_y_xy[i] * a_exp + 4.0 * g_x_z_yzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_yz_xz[i] = -2.0 * g_x_z_y_xz[i] * a_exp + 4.0 * g_x_z_yzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_yz_yy[i] = -2.0 * g_x_z_y_yy[i] * a_exp + 4.0 * g_x_z_yzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_yz_yz[i] = -2.0 * g_x_z_y_yz[i] * a_exp + 4.0 * g_x_z_yzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_yz_zz[i] = -2.0 * g_x_z_y_zz[i] * a_exp + 4.0 * g_x_z_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_x_0_z_0_0_z_zz_xx, g_x_0_z_0_0_z_zz_xy, g_x_0_z_0_0_z_zz_xz, g_x_0_z_0_0_z_zz_yy, g_x_0_z_0_0_z_zz_yz, g_x_0_z_0_0_z_zz_zz, g_x_z_z_xx, g_x_z_z_xy, g_x_z_z_xz, g_x_z_z_yy, g_x_z_z_yz, g_x_z_z_zz, g_x_z_zzz_xx, g_x_z_zzz_xy, g_x_z_zzz_xz, g_x_z_zzz_yy, g_x_z_zzz_yz, g_x_z_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_z_zz_xx[i] = -4.0 * g_x_z_z_xx[i] * a_exp + 4.0 * g_x_z_zzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_zz_xy[i] = -4.0 * g_x_z_z_xy[i] * a_exp + 4.0 * g_x_z_zzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_zz_xz[i] = -4.0 * g_x_z_z_xz[i] * a_exp + 4.0 * g_x_z_zzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_zz_yy[i] = -4.0 * g_x_z_z_yy[i] * a_exp + 4.0 * g_x_z_zzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_zz_yz[i] = -4.0 * g_x_z_z_yz[i] * a_exp + 4.0 * g_x_z_zzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_zz_zz[i] = -4.0 * g_x_z_z_zz[i] * a_exp + 4.0 * g_x_z_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (324-330)

    #pragma omp simd aligned(g_y_0_x_0_0_x_xx_xx, g_y_0_x_0_0_x_xx_xy, g_y_0_x_0_0_x_xx_xz, g_y_0_x_0_0_x_xx_yy, g_y_0_x_0_0_x_xx_yz, g_y_0_x_0_0_x_xx_zz, g_y_x_x_xx, g_y_x_x_xy, g_y_x_x_xz, g_y_x_x_yy, g_y_x_x_yz, g_y_x_x_zz, g_y_x_xxx_xx, g_y_x_xxx_xy, g_y_x_xxx_xz, g_y_x_xxx_yy, g_y_x_xxx_yz, g_y_x_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_x_xx_xx[i] = -4.0 * g_y_x_x_xx[i] * a_exp + 4.0 * g_y_x_xxx_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_xx_xy[i] = -4.0 * g_y_x_x_xy[i] * a_exp + 4.0 * g_y_x_xxx_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_xx_xz[i] = -4.0 * g_y_x_x_xz[i] * a_exp + 4.0 * g_y_x_xxx_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_xx_yy[i] = -4.0 * g_y_x_x_yy[i] * a_exp + 4.0 * g_y_x_xxx_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_xx_yz[i] = -4.0 * g_y_x_x_yz[i] * a_exp + 4.0 * g_y_x_xxx_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_xx_zz[i] = -4.0 * g_y_x_x_zz[i] * a_exp + 4.0 * g_y_x_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (330-336)

    #pragma omp simd aligned(g_y_0_x_0_0_x_xy_xx, g_y_0_x_0_0_x_xy_xy, g_y_0_x_0_0_x_xy_xz, g_y_0_x_0_0_x_xy_yy, g_y_0_x_0_0_x_xy_yz, g_y_0_x_0_0_x_xy_zz, g_y_x_xxy_xx, g_y_x_xxy_xy, g_y_x_xxy_xz, g_y_x_xxy_yy, g_y_x_xxy_yz, g_y_x_xxy_zz, g_y_x_y_xx, g_y_x_y_xy, g_y_x_y_xz, g_y_x_y_yy, g_y_x_y_yz, g_y_x_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_x_xy_xx[i] = -2.0 * g_y_x_y_xx[i] * a_exp + 4.0 * g_y_x_xxy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_xy_xy[i] = -2.0 * g_y_x_y_xy[i] * a_exp + 4.0 * g_y_x_xxy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_xy_xz[i] = -2.0 * g_y_x_y_xz[i] * a_exp + 4.0 * g_y_x_xxy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_xy_yy[i] = -2.0 * g_y_x_y_yy[i] * a_exp + 4.0 * g_y_x_xxy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_xy_yz[i] = -2.0 * g_y_x_y_yz[i] * a_exp + 4.0 * g_y_x_xxy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_xy_zz[i] = -2.0 * g_y_x_y_zz[i] * a_exp + 4.0 * g_y_x_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (336-342)

    #pragma omp simd aligned(g_y_0_x_0_0_x_xz_xx, g_y_0_x_0_0_x_xz_xy, g_y_0_x_0_0_x_xz_xz, g_y_0_x_0_0_x_xz_yy, g_y_0_x_0_0_x_xz_yz, g_y_0_x_0_0_x_xz_zz, g_y_x_xxz_xx, g_y_x_xxz_xy, g_y_x_xxz_xz, g_y_x_xxz_yy, g_y_x_xxz_yz, g_y_x_xxz_zz, g_y_x_z_xx, g_y_x_z_xy, g_y_x_z_xz, g_y_x_z_yy, g_y_x_z_yz, g_y_x_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_x_xz_xx[i] = -2.0 * g_y_x_z_xx[i] * a_exp + 4.0 * g_y_x_xxz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_xz_xy[i] = -2.0 * g_y_x_z_xy[i] * a_exp + 4.0 * g_y_x_xxz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_xz_xz[i] = -2.0 * g_y_x_z_xz[i] * a_exp + 4.0 * g_y_x_xxz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_xz_yy[i] = -2.0 * g_y_x_z_yy[i] * a_exp + 4.0 * g_y_x_xxz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_xz_yz[i] = -2.0 * g_y_x_z_yz[i] * a_exp + 4.0 * g_y_x_xxz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_xz_zz[i] = -2.0 * g_y_x_z_zz[i] * a_exp + 4.0 * g_y_x_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (342-348)

    #pragma omp simd aligned(g_y_0_x_0_0_x_yy_xx, g_y_0_x_0_0_x_yy_xy, g_y_0_x_0_0_x_yy_xz, g_y_0_x_0_0_x_yy_yy, g_y_0_x_0_0_x_yy_yz, g_y_0_x_0_0_x_yy_zz, g_y_x_xyy_xx, g_y_x_xyy_xy, g_y_x_xyy_xz, g_y_x_xyy_yy, g_y_x_xyy_yz, g_y_x_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_x_yy_xx[i] = 4.0 * g_y_x_xyy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_yy_xy[i] = 4.0 * g_y_x_xyy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_yy_xz[i] = 4.0 * g_y_x_xyy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_yy_yy[i] = 4.0 * g_y_x_xyy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_yy_yz[i] = 4.0 * g_y_x_xyy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_yy_zz[i] = 4.0 * g_y_x_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (348-354)

    #pragma omp simd aligned(g_y_0_x_0_0_x_yz_xx, g_y_0_x_0_0_x_yz_xy, g_y_0_x_0_0_x_yz_xz, g_y_0_x_0_0_x_yz_yy, g_y_0_x_0_0_x_yz_yz, g_y_0_x_0_0_x_yz_zz, g_y_x_xyz_xx, g_y_x_xyz_xy, g_y_x_xyz_xz, g_y_x_xyz_yy, g_y_x_xyz_yz, g_y_x_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_x_yz_xx[i] = 4.0 * g_y_x_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_yz_xy[i] = 4.0 * g_y_x_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_yz_xz[i] = 4.0 * g_y_x_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_yz_yy[i] = 4.0 * g_y_x_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_yz_yz[i] = 4.0 * g_y_x_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_yz_zz[i] = 4.0 * g_y_x_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (354-360)

    #pragma omp simd aligned(g_y_0_x_0_0_x_zz_xx, g_y_0_x_0_0_x_zz_xy, g_y_0_x_0_0_x_zz_xz, g_y_0_x_0_0_x_zz_yy, g_y_0_x_0_0_x_zz_yz, g_y_0_x_0_0_x_zz_zz, g_y_x_xzz_xx, g_y_x_xzz_xy, g_y_x_xzz_xz, g_y_x_xzz_yy, g_y_x_xzz_yz, g_y_x_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_x_zz_xx[i] = 4.0 * g_y_x_xzz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_zz_xy[i] = 4.0 * g_y_x_xzz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_zz_xz[i] = 4.0 * g_y_x_xzz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_zz_yy[i] = 4.0 * g_y_x_xzz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_zz_yz[i] = 4.0 * g_y_x_xzz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_zz_zz[i] = 4.0 * g_y_x_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (360-366)

    #pragma omp simd aligned(g_y_0_x_0_0_y_xx_xx, g_y_0_x_0_0_y_xx_xy, g_y_0_x_0_0_y_xx_xz, g_y_0_x_0_0_y_xx_yy, g_y_0_x_0_0_y_xx_yz, g_y_0_x_0_0_y_xx_zz, g_y_y_x_xx, g_y_y_x_xy, g_y_y_x_xz, g_y_y_x_yy, g_y_y_x_yz, g_y_y_x_zz, g_y_y_xxx_xx, g_y_y_xxx_xy, g_y_y_xxx_xz, g_y_y_xxx_yy, g_y_y_xxx_yz, g_y_y_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_y_xx_xx[i] = -4.0 * g_y_y_x_xx[i] * a_exp + 4.0 * g_y_y_xxx_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_xx_xy[i] = -4.0 * g_y_y_x_xy[i] * a_exp + 4.0 * g_y_y_xxx_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_xx_xz[i] = -4.0 * g_y_y_x_xz[i] * a_exp + 4.0 * g_y_y_xxx_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_xx_yy[i] = -4.0 * g_y_y_x_yy[i] * a_exp + 4.0 * g_y_y_xxx_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_xx_yz[i] = -4.0 * g_y_y_x_yz[i] * a_exp + 4.0 * g_y_y_xxx_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_xx_zz[i] = -4.0 * g_y_y_x_zz[i] * a_exp + 4.0 * g_y_y_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (366-372)

    #pragma omp simd aligned(g_y_0_x_0_0_y_xy_xx, g_y_0_x_0_0_y_xy_xy, g_y_0_x_0_0_y_xy_xz, g_y_0_x_0_0_y_xy_yy, g_y_0_x_0_0_y_xy_yz, g_y_0_x_0_0_y_xy_zz, g_y_y_xxy_xx, g_y_y_xxy_xy, g_y_y_xxy_xz, g_y_y_xxy_yy, g_y_y_xxy_yz, g_y_y_xxy_zz, g_y_y_y_xx, g_y_y_y_xy, g_y_y_y_xz, g_y_y_y_yy, g_y_y_y_yz, g_y_y_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_y_xy_xx[i] = -2.0 * g_y_y_y_xx[i] * a_exp + 4.0 * g_y_y_xxy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_xy_xy[i] = -2.0 * g_y_y_y_xy[i] * a_exp + 4.0 * g_y_y_xxy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_xy_xz[i] = -2.0 * g_y_y_y_xz[i] * a_exp + 4.0 * g_y_y_xxy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_xy_yy[i] = -2.0 * g_y_y_y_yy[i] * a_exp + 4.0 * g_y_y_xxy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_xy_yz[i] = -2.0 * g_y_y_y_yz[i] * a_exp + 4.0 * g_y_y_xxy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_xy_zz[i] = -2.0 * g_y_y_y_zz[i] * a_exp + 4.0 * g_y_y_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (372-378)

    #pragma omp simd aligned(g_y_0_x_0_0_y_xz_xx, g_y_0_x_0_0_y_xz_xy, g_y_0_x_0_0_y_xz_xz, g_y_0_x_0_0_y_xz_yy, g_y_0_x_0_0_y_xz_yz, g_y_0_x_0_0_y_xz_zz, g_y_y_xxz_xx, g_y_y_xxz_xy, g_y_y_xxz_xz, g_y_y_xxz_yy, g_y_y_xxz_yz, g_y_y_xxz_zz, g_y_y_z_xx, g_y_y_z_xy, g_y_y_z_xz, g_y_y_z_yy, g_y_y_z_yz, g_y_y_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_y_xz_xx[i] = -2.0 * g_y_y_z_xx[i] * a_exp + 4.0 * g_y_y_xxz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_xz_xy[i] = -2.0 * g_y_y_z_xy[i] * a_exp + 4.0 * g_y_y_xxz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_xz_xz[i] = -2.0 * g_y_y_z_xz[i] * a_exp + 4.0 * g_y_y_xxz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_xz_yy[i] = -2.0 * g_y_y_z_yy[i] * a_exp + 4.0 * g_y_y_xxz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_xz_yz[i] = -2.0 * g_y_y_z_yz[i] * a_exp + 4.0 * g_y_y_xxz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_xz_zz[i] = -2.0 * g_y_y_z_zz[i] * a_exp + 4.0 * g_y_y_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (378-384)

    #pragma omp simd aligned(g_y_0_x_0_0_y_yy_xx, g_y_0_x_0_0_y_yy_xy, g_y_0_x_0_0_y_yy_xz, g_y_0_x_0_0_y_yy_yy, g_y_0_x_0_0_y_yy_yz, g_y_0_x_0_0_y_yy_zz, g_y_y_xyy_xx, g_y_y_xyy_xy, g_y_y_xyy_xz, g_y_y_xyy_yy, g_y_y_xyy_yz, g_y_y_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_y_yy_xx[i] = 4.0 * g_y_y_xyy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_yy_xy[i] = 4.0 * g_y_y_xyy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_yy_xz[i] = 4.0 * g_y_y_xyy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_yy_yy[i] = 4.0 * g_y_y_xyy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_yy_yz[i] = 4.0 * g_y_y_xyy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_yy_zz[i] = 4.0 * g_y_y_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (384-390)

    #pragma omp simd aligned(g_y_0_x_0_0_y_yz_xx, g_y_0_x_0_0_y_yz_xy, g_y_0_x_0_0_y_yz_xz, g_y_0_x_0_0_y_yz_yy, g_y_0_x_0_0_y_yz_yz, g_y_0_x_0_0_y_yz_zz, g_y_y_xyz_xx, g_y_y_xyz_xy, g_y_y_xyz_xz, g_y_y_xyz_yy, g_y_y_xyz_yz, g_y_y_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_y_yz_xx[i] = 4.0 * g_y_y_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_yz_xy[i] = 4.0 * g_y_y_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_yz_xz[i] = 4.0 * g_y_y_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_yz_yy[i] = 4.0 * g_y_y_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_yz_yz[i] = 4.0 * g_y_y_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_yz_zz[i] = 4.0 * g_y_y_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (390-396)

    #pragma omp simd aligned(g_y_0_x_0_0_y_zz_xx, g_y_0_x_0_0_y_zz_xy, g_y_0_x_0_0_y_zz_xz, g_y_0_x_0_0_y_zz_yy, g_y_0_x_0_0_y_zz_yz, g_y_0_x_0_0_y_zz_zz, g_y_y_xzz_xx, g_y_y_xzz_xy, g_y_y_xzz_xz, g_y_y_xzz_yy, g_y_y_xzz_yz, g_y_y_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_y_zz_xx[i] = 4.0 * g_y_y_xzz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_zz_xy[i] = 4.0 * g_y_y_xzz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_zz_xz[i] = 4.0 * g_y_y_xzz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_zz_yy[i] = 4.0 * g_y_y_xzz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_zz_yz[i] = 4.0 * g_y_y_xzz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_zz_zz[i] = 4.0 * g_y_y_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (396-402)

    #pragma omp simd aligned(g_y_0_x_0_0_z_xx_xx, g_y_0_x_0_0_z_xx_xy, g_y_0_x_0_0_z_xx_xz, g_y_0_x_0_0_z_xx_yy, g_y_0_x_0_0_z_xx_yz, g_y_0_x_0_0_z_xx_zz, g_y_z_x_xx, g_y_z_x_xy, g_y_z_x_xz, g_y_z_x_yy, g_y_z_x_yz, g_y_z_x_zz, g_y_z_xxx_xx, g_y_z_xxx_xy, g_y_z_xxx_xz, g_y_z_xxx_yy, g_y_z_xxx_yz, g_y_z_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_z_xx_xx[i] = -4.0 * g_y_z_x_xx[i] * a_exp + 4.0 * g_y_z_xxx_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_xx_xy[i] = -4.0 * g_y_z_x_xy[i] * a_exp + 4.0 * g_y_z_xxx_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_xx_xz[i] = -4.0 * g_y_z_x_xz[i] * a_exp + 4.0 * g_y_z_xxx_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_xx_yy[i] = -4.0 * g_y_z_x_yy[i] * a_exp + 4.0 * g_y_z_xxx_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_xx_yz[i] = -4.0 * g_y_z_x_yz[i] * a_exp + 4.0 * g_y_z_xxx_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_xx_zz[i] = -4.0 * g_y_z_x_zz[i] * a_exp + 4.0 * g_y_z_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (402-408)

    #pragma omp simd aligned(g_y_0_x_0_0_z_xy_xx, g_y_0_x_0_0_z_xy_xy, g_y_0_x_0_0_z_xy_xz, g_y_0_x_0_0_z_xy_yy, g_y_0_x_0_0_z_xy_yz, g_y_0_x_0_0_z_xy_zz, g_y_z_xxy_xx, g_y_z_xxy_xy, g_y_z_xxy_xz, g_y_z_xxy_yy, g_y_z_xxy_yz, g_y_z_xxy_zz, g_y_z_y_xx, g_y_z_y_xy, g_y_z_y_xz, g_y_z_y_yy, g_y_z_y_yz, g_y_z_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_z_xy_xx[i] = -2.0 * g_y_z_y_xx[i] * a_exp + 4.0 * g_y_z_xxy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_xy_xy[i] = -2.0 * g_y_z_y_xy[i] * a_exp + 4.0 * g_y_z_xxy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_xy_xz[i] = -2.0 * g_y_z_y_xz[i] * a_exp + 4.0 * g_y_z_xxy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_xy_yy[i] = -2.0 * g_y_z_y_yy[i] * a_exp + 4.0 * g_y_z_xxy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_xy_yz[i] = -2.0 * g_y_z_y_yz[i] * a_exp + 4.0 * g_y_z_xxy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_xy_zz[i] = -2.0 * g_y_z_y_zz[i] * a_exp + 4.0 * g_y_z_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (408-414)

    #pragma omp simd aligned(g_y_0_x_0_0_z_xz_xx, g_y_0_x_0_0_z_xz_xy, g_y_0_x_0_0_z_xz_xz, g_y_0_x_0_0_z_xz_yy, g_y_0_x_0_0_z_xz_yz, g_y_0_x_0_0_z_xz_zz, g_y_z_xxz_xx, g_y_z_xxz_xy, g_y_z_xxz_xz, g_y_z_xxz_yy, g_y_z_xxz_yz, g_y_z_xxz_zz, g_y_z_z_xx, g_y_z_z_xy, g_y_z_z_xz, g_y_z_z_yy, g_y_z_z_yz, g_y_z_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_z_xz_xx[i] = -2.0 * g_y_z_z_xx[i] * a_exp + 4.0 * g_y_z_xxz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_xz_xy[i] = -2.0 * g_y_z_z_xy[i] * a_exp + 4.0 * g_y_z_xxz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_xz_xz[i] = -2.0 * g_y_z_z_xz[i] * a_exp + 4.0 * g_y_z_xxz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_xz_yy[i] = -2.0 * g_y_z_z_yy[i] * a_exp + 4.0 * g_y_z_xxz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_xz_yz[i] = -2.0 * g_y_z_z_yz[i] * a_exp + 4.0 * g_y_z_xxz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_xz_zz[i] = -2.0 * g_y_z_z_zz[i] * a_exp + 4.0 * g_y_z_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (414-420)

    #pragma omp simd aligned(g_y_0_x_0_0_z_yy_xx, g_y_0_x_0_0_z_yy_xy, g_y_0_x_0_0_z_yy_xz, g_y_0_x_0_0_z_yy_yy, g_y_0_x_0_0_z_yy_yz, g_y_0_x_0_0_z_yy_zz, g_y_z_xyy_xx, g_y_z_xyy_xy, g_y_z_xyy_xz, g_y_z_xyy_yy, g_y_z_xyy_yz, g_y_z_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_z_yy_xx[i] = 4.0 * g_y_z_xyy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_yy_xy[i] = 4.0 * g_y_z_xyy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_yy_xz[i] = 4.0 * g_y_z_xyy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_yy_yy[i] = 4.0 * g_y_z_xyy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_yy_yz[i] = 4.0 * g_y_z_xyy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_yy_zz[i] = 4.0 * g_y_z_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (420-426)

    #pragma omp simd aligned(g_y_0_x_0_0_z_yz_xx, g_y_0_x_0_0_z_yz_xy, g_y_0_x_0_0_z_yz_xz, g_y_0_x_0_0_z_yz_yy, g_y_0_x_0_0_z_yz_yz, g_y_0_x_0_0_z_yz_zz, g_y_z_xyz_xx, g_y_z_xyz_xy, g_y_z_xyz_xz, g_y_z_xyz_yy, g_y_z_xyz_yz, g_y_z_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_z_yz_xx[i] = 4.0 * g_y_z_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_yz_xy[i] = 4.0 * g_y_z_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_yz_xz[i] = 4.0 * g_y_z_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_yz_yy[i] = 4.0 * g_y_z_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_yz_yz[i] = 4.0 * g_y_z_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_yz_zz[i] = 4.0 * g_y_z_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (426-432)

    #pragma omp simd aligned(g_y_0_x_0_0_z_zz_xx, g_y_0_x_0_0_z_zz_xy, g_y_0_x_0_0_z_zz_xz, g_y_0_x_0_0_z_zz_yy, g_y_0_x_0_0_z_zz_yz, g_y_0_x_0_0_z_zz_zz, g_y_z_xzz_xx, g_y_z_xzz_xy, g_y_z_xzz_xz, g_y_z_xzz_yy, g_y_z_xzz_yz, g_y_z_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_z_zz_xx[i] = 4.0 * g_y_z_xzz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_zz_xy[i] = 4.0 * g_y_z_xzz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_zz_xz[i] = 4.0 * g_y_z_xzz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_zz_yy[i] = 4.0 * g_y_z_xzz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_zz_yz[i] = 4.0 * g_y_z_xzz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_zz_zz[i] = 4.0 * g_y_z_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (432-438)

    #pragma omp simd aligned(g_y_0_y_0_0_x_xx_xx, g_y_0_y_0_0_x_xx_xy, g_y_0_y_0_0_x_xx_xz, g_y_0_y_0_0_x_xx_yy, g_y_0_y_0_0_x_xx_yz, g_y_0_y_0_0_x_xx_zz, g_y_x_xxy_xx, g_y_x_xxy_xy, g_y_x_xxy_xz, g_y_x_xxy_yy, g_y_x_xxy_yz, g_y_x_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_x_xx_xx[i] = 4.0 * g_y_x_xxy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_xx_xy[i] = 4.0 * g_y_x_xxy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_xx_xz[i] = 4.0 * g_y_x_xxy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_xx_yy[i] = 4.0 * g_y_x_xxy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_xx_yz[i] = 4.0 * g_y_x_xxy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_xx_zz[i] = 4.0 * g_y_x_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (438-444)

    #pragma omp simd aligned(g_y_0_y_0_0_x_xy_xx, g_y_0_y_0_0_x_xy_xy, g_y_0_y_0_0_x_xy_xz, g_y_0_y_0_0_x_xy_yy, g_y_0_y_0_0_x_xy_yz, g_y_0_y_0_0_x_xy_zz, g_y_x_x_xx, g_y_x_x_xy, g_y_x_x_xz, g_y_x_x_yy, g_y_x_x_yz, g_y_x_x_zz, g_y_x_xyy_xx, g_y_x_xyy_xy, g_y_x_xyy_xz, g_y_x_xyy_yy, g_y_x_xyy_yz, g_y_x_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_x_xy_xx[i] = -2.0 * g_y_x_x_xx[i] * a_exp + 4.0 * g_y_x_xyy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_xy_xy[i] = -2.0 * g_y_x_x_xy[i] * a_exp + 4.0 * g_y_x_xyy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_xy_xz[i] = -2.0 * g_y_x_x_xz[i] * a_exp + 4.0 * g_y_x_xyy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_xy_yy[i] = -2.0 * g_y_x_x_yy[i] * a_exp + 4.0 * g_y_x_xyy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_xy_yz[i] = -2.0 * g_y_x_x_yz[i] * a_exp + 4.0 * g_y_x_xyy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_xy_zz[i] = -2.0 * g_y_x_x_zz[i] * a_exp + 4.0 * g_y_x_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (444-450)

    #pragma omp simd aligned(g_y_0_y_0_0_x_xz_xx, g_y_0_y_0_0_x_xz_xy, g_y_0_y_0_0_x_xz_xz, g_y_0_y_0_0_x_xz_yy, g_y_0_y_0_0_x_xz_yz, g_y_0_y_0_0_x_xz_zz, g_y_x_xyz_xx, g_y_x_xyz_xy, g_y_x_xyz_xz, g_y_x_xyz_yy, g_y_x_xyz_yz, g_y_x_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_x_xz_xx[i] = 4.0 * g_y_x_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_xz_xy[i] = 4.0 * g_y_x_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_xz_xz[i] = 4.0 * g_y_x_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_xz_yy[i] = 4.0 * g_y_x_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_xz_yz[i] = 4.0 * g_y_x_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_xz_zz[i] = 4.0 * g_y_x_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (450-456)

    #pragma omp simd aligned(g_y_0_y_0_0_x_yy_xx, g_y_0_y_0_0_x_yy_xy, g_y_0_y_0_0_x_yy_xz, g_y_0_y_0_0_x_yy_yy, g_y_0_y_0_0_x_yy_yz, g_y_0_y_0_0_x_yy_zz, g_y_x_y_xx, g_y_x_y_xy, g_y_x_y_xz, g_y_x_y_yy, g_y_x_y_yz, g_y_x_y_zz, g_y_x_yyy_xx, g_y_x_yyy_xy, g_y_x_yyy_xz, g_y_x_yyy_yy, g_y_x_yyy_yz, g_y_x_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_x_yy_xx[i] = -4.0 * g_y_x_y_xx[i] * a_exp + 4.0 * g_y_x_yyy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_yy_xy[i] = -4.0 * g_y_x_y_xy[i] * a_exp + 4.0 * g_y_x_yyy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_yy_xz[i] = -4.0 * g_y_x_y_xz[i] * a_exp + 4.0 * g_y_x_yyy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_yy_yy[i] = -4.0 * g_y_x_y_yy[i] * a_exp + 4.0 * g_y_x_yyy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_yy_yz[i] = -4.0 * g_y_x_y_yz[i] * a_exp + 4.0 * g_y_x_yyy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_yy_zz[i] = -4.0 * g_y_x_y_zz[i] * a_exp + 4.0 * g_y_x_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (456-462)

    #pragma omp simd aligned(g_y_0_y_0_0_x_yz_xx, g_y_0_y_0_0_x_yz_xy, g_y_0_y_0_0_x_yz_xz, g_y_0_y_0_0_x_yz_yy, g_y_0_y_0_0_x_yz_yz, g_y_0_y_0_0_x_yz_zz, g_y_x_yyz_xx, g_y_x_yyz_xy, g_y_x_yyz_xz, g_y_x_yyz_yy, g_y_x_yyz_yz, g_y_x_yyz_zz, g_y_x_z_xx, g_y_x_z_xy, g_y_x_z_xz, g_y_x_z_yy, g_y_x_z_yz, g_y_x_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_x_yz_xx[i] = -2.0 * g_y_x_z_xx[i] * a_exp + 4.0 * g_y_x_yyz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_yz_xy[i] = -2.0 * g_y_x_z_xy[i] * a_exp + 4.0 * g_y_x_yyz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_yz_xz[i] = -2.0 * g_y_x_z_xz[i] * a_exp + 4.0 * g_y_x_yyz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_yz_yy[i] = -2.0 * g_y_x_z_yy[i] * a_exp + 4.0 * g_y_x_yyz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_yz_yz[i] = -2.0 * g_y_x_z_yz[i] * a_exp + 4.0 * g_y_x_yyz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_yz_zz[i] = -2.0 * g_y_x_z_zz[i] * a_exp + 4.0 * g_y_x_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (462-468)

    #pragma omp simd aligned(g_y_0_y_0_0_x_zz_xx, g_y_0_y_0_0_x_zz_xy, g_y_0_y_0_0_x_zz_xz, g_y_0_y_0_0_x_zz_yy, g_y_0_y_0_0_x_zz_yz, g_y_0_y_0_0_x_zz_zz, g_y_x_yzz_xx, g_y_x_yzz_xy, g_y_x_yzz_xz, g_y_x_yzz_yy, g_y_x_yzz_yz, g_y_x_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_x_zz_xx[i] = 4.0 * g_y_x_yzz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_zz_xy[i] = 4.0 * g_y_x_yzz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_zz_xz[i] = 4.0 * g_y_x_yzz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_zz_yy[i] = 4.0 * g_y_x_yzz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_zz_yz[i] = 4.0 * g_y_x_yzz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_zz_zz[i] = 4.0 * g_y_x_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (468-474)

    #pragma omp simd aligned(g_y_0_y_0_0_y_xx_xx, g_y_0_y_0_0_y_xx_xy, g_y_0_y_0_0_y_xx_xz, g_y_0_y_0_0_y_xx_yy, g_y_0_y_0_0_y_xx_yz, g_y_0_y_0_0_y_xx_zz, g_y_y_xxy_xx, g_y_y_xxy_xy, g_y_y_xxy_xz, g_y_y_xxy_yy, g_y_y_xxy_yz, g_y_y_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_y_xx_xx[i] = 4.0 * g_y_y_xxy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_xx_xy[i] = 4.0 * g_y_y_xxy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_xx_xz[i] = 4.0 * g_y_y_xxy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_xx_yy[i] = 4.0 * g_y_y_xxy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_xx_yz[i] = 4.0 * g_y_y_xxy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_xx_zz[i] = 4.0 * g_y_y_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (474-480)

    #pragma omp simd aligned(g_y_0_y_0_0_y_xy_xx, g_y_0_y_0_0_y_xy_xy, g_y_0_y_0_0_y_xy_xz, g_y_0_y_0_0_y_xy_yy, g_y_0_y_0_0_y_xy_yz, g_y_0_y_0_0_y_xy_zz, g_y_y_x_xx, g_y_y_x_xy, g_y_y_x_xz, g_y_y_x_yy, g_y_y_x_yz, g_y_y_x_zz, g_y_y_xyy_xx, g_y_y_xyy_xy, g_y_y_xyy_xz, g_y_y_xyy_yy, g_y_y_xyy_yz, g_y_y_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_y_xy_xx[i] = -2.0 * g_y_y_x_xx[i] * a_exp + 4.0 * g_y_y_xyy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_xy_xy[i] = -2.0 * g_y_y_x_xy[i] * a_exp + 4.0 * g_y_y_xyy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_xy_xz[i] = -2.0 * g_y_y_x_xz[i] * a_exp + 4.0 * g_y_y_xyy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_xy_yy[i] = -2.0 * g_y_y_x_yy[i] * a_exp + 4.0 * g_y_y_xyy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_xy_yz[i] = -2.0 * g_y_y_x_yz[i] * a_exp + 4.0 * g_y_y_xyy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_xy_zz[i] = -2.0 * g_y_y_x_zz[i] * a_exp + 4.0 * g_y_y_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (480-486)

    #pragma omp simd aligned(g_y_0_y_0_0_y_xz_xx, g_y_0_y_0_0_y_xz_xy, g_y_0_y_0_0_y_xz_xz, g_y_0_y_0_0_y_xz_yy, g_y_0_y_0_0_y_xz_yz, g_y_0_y_0_0_y_xz_zz, g_y_y_xyz_xx, g_y_y_xyz_xy, g_y_y_xyz_xz, g_y_y_xyz_yy, g_y_y_xyz_yz, g_y_y_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_y_xz_xx[i] = 4.0 * g_y_y_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_xz_xy[i] = 4.0 * g_y_y_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_xz_xz[i] = 4.0 * g_y_y_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_xz_yy[i] = 4.0 * g_y_y_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_xz_yz[i] = 4.0 * g_y_y_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_xz_zz[i] = 4.0 * g_y_y_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (486-492)

    #pragma omp simd aligned(g_y_0_y_0_0_y_yy_xx, g_y_0_y_0_0_y_yy_xy, g_y_0_y_0_0_y_yy_xz, g_y_0_y_0_0_y_yy_yy, g_y_0_y_0_0_y_yy_yz, g_y_0_y_0_0_y_yy_zz, g_y_y_y_xx, g_y_y_y_xy, g_y_y_y_xz, g_y_y_y_yy, g_y_y_y_yz, g_y_y_y_zz, g_y_y_yyy_xx, g_y_y_yyy_xy, g_y_y_yyy_xz, g_y_y_yyy_yy, g_y_y_yyy_yz, g_y_y_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_y_yy_xx[i] = -4.0 * g_y_y_y_xx[i] * a_exp + 4.0 * g_y_y_yyy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_yy_xy[i] = -4.0 * g_y_y_y_xy[i] * a_exp + 4.0 * g_y_y_yyy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_yy_xz[i] = -4.0 * g_y_y_y_xz[i] * a_exp + 4.0 * g_y_y_yyy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_yy_yy[i] = -4.0 * g_y_y_y_yy[i] * a_exp + 4.0 * g_y_y_yyy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_yy_yz[i] = -4.0 * g_y_y_y_yz[i] * a_exp + 4.0 * g_y_y_yyy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_yy_zz[i] = -4.0 * g_y_y_y_zz[i] * a_exp + 4.0 * g_y_y_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (492-498)

    #pragma omp simd aligned(g_y_0_y_0_0_y_yz_xx, g_y_0_y_0_0_y_yz_xy, g_y_0_y_0_0_y_yz_xz, g_y_0_y_0_0_y_yz_yy, g_y_0_y_0_0_y_yz_yz, g_y_0_y_0_0_y_yz_zz, g_y_y_yyz_xx, g_y_y_yyz_xy, g_y_y_yyz_xz, g_y_y_yyz_yy, g_y_y_yyz_yz, g_y_y_yyz_zz, g_y_y_z_xx, g_y_y_z_xy, g_y_y_z_xz, g_y_y_z_yy, g_y_y_z_yz, g_y_y_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_y_yz_xx[i] = -2.0 * g_y_y_z_xx[i] * a_exp + 4.0 * g_y_y_yyz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_yz_xy[i] = -2.0 * g_y_y_z_xy[i] * a_exp + 4.0 * g_y_y_yyz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_yz_xz[i] = -2.0 * g_y_y_z_xz[i] * a_exp + 4.0 * g_y_y_yyz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_yz_yy[i] = -2.0 * g_y_y_z_yy[i] * a_exp + 4.0 * g_y_y_yyz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_yz_yz[i] = -2.0 * g_y_y_z_yz[i] * a_exp + 4.0 * g_y_y_yyz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_yz_zz[i] = -2.0 * g_y_y_z_zz[i] * a_exp + 4.0 * g_y_y_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (498-504)

    #pragma omp simd aligned(g_y_0_y_0_0_y_zz_xx, g_y_0_y_0_0_y_zz_xy, g_y_0_y_0_0_y_zz_xz, g_y_0_y_0_0_y_zz_yy, g_y_0_y_0_0_y_zz_yz, g_y_0_y_0_0_y_zz_zz, g_y_y_yzz_xx, g_y_y_yzz_xy, g_y_y_yzz_xz, g_y_y_yzz_yy, g_y_y_yzz_yz, g_y_y_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_y_zz_xx[i] = 4.0 * g_y_y_yzz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_zz_xy[i] = 4.0 * g_y_y_yzz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_zz_xz[i] = 4.0 * g_y_y_yzz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_zz_yy[i] = 4.0 * g_y_y_yzz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_zz_yz[i] = 4.0 * g_y_y_yzz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_zz_zz[i] = 4.0 * g_y_y_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (504-510)

    #pragma omp simd aligned(g_y_0_y_0_0_z_xx_xx, g_y_0_y_0_0_z_xx_xy, g_y_0_y_0_0_z_xx_xz, g_y_0_y_0_0_z_xx_yy, g_y_0_y_0_0_z_xx_yz, g_y_0_y_0_0_z_xx_zz, g_y_z_xxy_xx, g_y_z_xxy_xy, g_y_z_xxy_xz, g_y_z_xxy_yy, g_y_z_xxy_yz, g_y_z_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_z_xx_xx[i] = 4.0 * g_y_z_xxy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_xx_xy[i] = 4.0 * g_y_z_xxy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_xx_xz[i] = 4.0 * g_y_z_xxy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_xx_yy[i] = 4.0 * g_y_z_xxy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_xx_yz[i] = 4.0 * g_y_z_xxy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_xx_zz[i] = 4.0 * g_y_z_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (510-516)

    #pragma omp simd aligned(g_y_0_y_0_0_z_xy_xx, g_y_0_y_0_0_z_xy_xy, g_y_0_y_0_0_z_xy_xz, g_y_0_y_0_0_z_xy_yy, g_y_0_y_0_0_z_xy_yz, g_y_0_y_0_0_z_xy_zz, g_y_z_x_xx, g_y_z_x_xy, g_y_z_x_xz, g_y_z_x_yy, g_y_z_x_yz, g_y_z_x_zz, g_y_z_xyy_xx, g_y_z_xyy_xy, g_y_z_xyy_xz, g_y_z_xyy_yy, g_y_z_xyy_yz, g_y_z_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_z_xy_xx[i] = -2.0 * g_y_z_x_xx[i] * a_exp + 4.0 * g_y_z_xyy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_xy_xy[i] = -2.0 * g_y_z_x_xy[i] * a_exp + 4.0 * g_y_z_xyy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_xy_xz[i] = -2.0 * g_y_z_x_xz[i] * a_exp + 4.0 * g_y_z_xyy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_xy_yy[i] = -2.0 * g_y_z_x_yy[i] * a_exp + 4.0 * g_y_z_xyy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_xy_yz[i] = -2.0 * g_y_z_x_yz[i] * a_exp + 4.0 * g_y_z_xyy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_xy_zz[i] = -2.0 * g_y_z_x_zz[i] * a_exp + 4.0 * g_y_z_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (516-522)

    #pragma omp simd aligned(g_y_0_y_0_0_z_xz_xx, g_y_0_y_0_0_z_xz_xy, g_y_0_y_0_0_z_xz_xz, g_y_0_y_0_0_z_xz_yy, g_y_0_y_0_0_z_xz_yz, g_y_0_y_0_0_z_xz_zz, g_y_z_xyz_xx, g_y_z_xyz_xy, g_y_z_xyz_xz, g_y_z_xyz_yy, g_y_z_xyz_yz, g_y_z_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_z_xz_xx[i] = 4.0 * g_y_z_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_xz_xy[i] = 4.0 * g_y_z_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_xz_xz[i] = 4.0 * g_y_z_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_xz_yy[i] = 4.0 * g_y_z_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_xz_yz[i] = 4.0 * g_y_z_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_xz_zz[i] = 4.0 * g_y_z_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (522-528)

    #pragma omp simd aligned(g_y_0_y_0_0_z_yy_xx, g_y_0_y_0_0_z_yy_xy, g_y_0_y_0_0_z_yy_xz, g_y_0_y_0_0_z_yy_yy, g_y_0_y_0_0_z_yy_yz, g_y_0_y_0_0_z_yy_zz, g_y_z_y_xx, g_y_z_y_xy, g_y_z_y_xz, g_y_z_y_yy, g_y_z_y_yz, g_y_z_y_zz, g_y_z_yyy_xx, g_y_z_yyy_xy, g_y_z_yyy_xz, g_y_z_yyy_yy, g_y_z_yyy_yz, g_y_z_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_z_yy_xx[i] = -4.0 * g_y_z_y_xx[i] * a_exp + 4.0 * g_y_z_yyy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_yy_xy[i] = -4.0 * g_y_z_y_xy[i] * a_exp + 4.0 * g_y_z_yyy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_yy_xz[i] = -4.0 * g_y_z_y_xz[i] * a_exp + 4.0 * g_y_z_yyy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_yy_yy[i] = -4.0 * g_y_z_y_yy[i] * a_exp + 4.0 * g_y_z_yyy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_yy_yz[i] = -4.0 * g_y_z_y_yz[i] * a_exp + 4.0 * g_y_z_yyy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_yy_zz[i] = -4.0 * g_y_z_y_zz[i] * a_exp + 4.0 * g_y_z_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (528-534)

    #pragma omp simd aligned(g_y_0_y_0_0_z_yz_xx, g_y_0_y_0_0_z_yz_xy, g_y_0_y_0_0_z_yz_xz, g_y_0_y_0_0_z_yz_yy, g_y_0_y_0_0_z_yz_yz, g_y_0_y_0_0_z_yz_zz, g_y_z_yyz_xx, g_y_z_yyz_xy, g_y_z_yyz_xz, g_y_z_yyz_yy, g_y_z_yyz_yz, g_y_z_yyz_zz, g_y_z_z_xx, g_y_z_z_xy, g_y_z_z_xz, g_y_z_z_yy, g_y_z_z_yz, g_y_z_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_z_yz_xx[i] = -2.0 * g_y_z_z_xx[i] * a_exp + 4.0 * g_y_z_yyz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_yz_xy[i] = -2.0 * g_y_z_z_xy[i] * a_exp + 4.0 * g_y_z_yyz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_yz_xz[i] = -2.0 * g_y_z_z_xz[i] * a_exp + 4.0 * g_y_z_yyz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_yz_yy[i] = -2.0 * g_y_z_z_yy[i] * a_exp + 4.0 * g_y_z_yyz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_yz_yz[i] = -2.0 * g_y_z_z_yz[i] * a_exp + 4.0 * g_y_z_yyz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_yz_zz[i] = -2.0 * g_y_z_z_zz[i] * a_exp + 4.0 * g_y_z_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (534-540)

    #pragma omp simd aligned(g_y_0_y_0_0_z_zz_xx, g_y_0_y_0_0_z_zz_xy, g_y_0_y_0_0_z_zz_xz, g_y_0_y_0_0_z_zz_yy, g_y_0_y_0_0_z_zz_yz, g_y_0_y_0_0_z_zz_zz, g_y_z_yzz_xx, g_y_z_yzz_xy, g_y_z_yzz_xz, g_y_z_yzz_yy, g_y_z_yzz_yz, g_y_z_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_z_zz_xx[i] = 4.0 * g_y_z_yzz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_zz_xy[i] = 4.0 * g_y_z_yzz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_zz_xz[i] = 4.0 * g_y_z_yzz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_zz_yy[i] = 4.0 * g_y_z_yzz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_zz_yz[i] = 4.0 * g_y_z_yzz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_zz_zz[i] = 4.0 * g_y_z_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (540-546)

    #pragma omp simd aligned(g_y_0_z_0_0_x_xx_xx, g_y_0_z_0_0_x_xx_xy, g_y_0_z_0_0_x_xx_xz, g_y_0_z_0_0_x_xx_yy, g_y_0_z_0_0_x_xx_yz, g_y_0_z_0_0_x_xx_zz, g_y_x_xxz_xx, g_y_x_xxz_xy, g_y_x_xxz_xz, g_y_x_xxz_yy, g_y_x_xxz_yz, g_y_x_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_x_xx_xx[i] = 4.0 * g_y_x_xxz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_xx_xy[i] = 4.0 * g_y_x_xxz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_xx_xz[i] = 4.0 * g_y_x_xxz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_xx_yy[i] = 4.0 * g_y_x_xxz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_xx_yz[i] = 4.0 * g_y_x_xxz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_xx_zz[i] = 4.0 * g_y_x_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (546-552)

    #pragma omp simd aligned(g_y_0_z_0_0_x_xy_xx, g_y_0_z_0_0_x_xy_xy, g_y_0_z_0_0_x_xy_xz, g_y_0_z_0_0_x_xy_yy, g_y_0_z_0_0_x_xy_yz, g_y_0_z_0_0_x_xy_zz, g_y_x_xyz_xx, g_y_x_xyz_xy, g_y_x_xyz_xz, g_y_x_xyz_yy, g_y_x_xyz_yz, g_y_x_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_x_xy_xx[i] = 4.0 * g_y_x_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_xy_xy[i] = 4.0 * g_y_x_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_xy_xz[i] = 4.0 * g_y_x_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_xy_yy[i] = 4.0 * g_y_x_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_xy_yz[i] = 4.0 * g_y_x_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_xy_zz[i] = 4.0 * g_y_x_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (552-558)

    #pragma omp simd aligned(g_y_0_z_0_0_x_xz_xx, g_y_0_z_0_0_x_xz_xy, g_y_0_z_0_0_x_xz_xz, g_y_0_z_0_0_x_xz_yy, g_y_0_z_0_0_x_xz_yz, g_y_0_z_0_0_x_xz_zz, g_y_x_x_xx, g_y_x_x_xy, g_y_x_x_xz, g_y_x_x_yy, g_y_x_x_yz, g_y_x_x_zz, g_y_x_xzz_xx, g_y_x_xzz_xy, g_y_x_xzz_xz, g_y_x_xzz_yy, g_y_x_xzz_yz, g_y_x_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_x_xz_xx[i] = -2.0 * g_y_x_x_xx[i] * a_exp + 4.0 * g_y_x_xzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_xz_xy[i] = -2.0 * g_y_x_x_xy[i] * a_exp + 4.0 * g_y_x_xzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_xz_xz[i] = -2.0 * g_y_x_x_xz[i] * a_exp + 4.0 * g_y_x_xzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_xz_yy[i] = -2.0 * g_y_x_x_yy[i] * a_exp + 4.0 * g_y_x_xzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_xz_yz[i] = -2.0 * g_y_x_x_yz[i] * a_exp + 4.0 * g_y_x_xzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_xz_zz[i] = -2.0 * g_y_x_x_zz[i] * a_exp + 4.0 * g_y_x_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (558-564)

    #pragma omp simd aligned(g_y_0_z_0_0_x_yy_xx, g_y_0_z_0_0_x_yy_xy, g_y_0_z_0_0_x_yy_xz, g_y_0_z_0_0_x_yy_yy, g_y_0_z_0_0_x_yy_yz, g_y_0_z_0_0_x_yy_zz, g_y_x_yyz_xx, g_y_x_yyz_xy, g_y_x_yyz_xz, g_y_x_yyz_yy, g_y_x_yyz_yz, g_y_x_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_x_yy_xx[i] = 4.0 * g_y_x_yyz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_yy_xy[i] = 4.0 * g_y_x_yyz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_yy_xz[i] = 4.0 * g_y_x_yyz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_yy_yy[i] = 4.0 * g_y_x_yyz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_yy_yz[i] = 4.0 * g_y_x_yyz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_yy_zz[i] = 4.0 * g_y_x_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (564-570)

    #pragma omp simd aligned(g_y_0_z_0_0_x_yz_xx, g_y_0_z_0_0_x_yz_xy, g_y_0_z_0_0_x_yz_xz, g_y_0_z_0_0_x_yz_yy, g_y_0_z_0_0_x_yz_yz, g_y_0_z_0_0_x_yz_zz, g_y_x_y_xx, g_y_x_y_xy, g_y_x_y_xz, g_y_x_y_yy, g_y_x_y_yz, g_y_x_y_zz, g_y_x_yzz_xx, g_y_x_yzz_xy, g_y_x_yzz_xz, g_y_x_yzz_yy, g_y_x_yzz_yz, g_y_x_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_x_yz_xx[i] = -2.0 * g_y_x_y_xx[i] * a_exp + 4.0 * g_y_x_yzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_yz_xy[i] = -2.0 * g_y_x_y_xy[i] * a_exp + 4.0 * g_y_x_yzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_yz_xz[i] = -2.0 * g_y_x_y_xz[i] * a_exp + 4.0 * g_y_x_yzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_yz_yy[i] = -2.0 * g_y_x_y_yy[i] * a_exp + 4.0 * g_y_x_yzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_yz_yz[i] = -2.0 * g_y_x_y_yz[i] * a_exp + 4.0 * g_y_x_yzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_yz_zz[i] = -2.0 * g_y_x_y_zz[i] * a_exp + 4.0 * g_y_x_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (570-576)

    #pragma omp simd aligned(g_y_0_z_0_0_x_zz_xx, g_y_0_z_0_0_x_zz_xy, g_y_0_z_0_0_x_zz_xz, g_y_0_z_0_0_x_zz_yy, g_y_0_z_0_0_x_zz_yz, g_y_0_z_0_0_x_zz_zz, g_y_x_z_xx, g_y_x_z_xy, g_y_x_z_xz, g_y_x_z_yy, g_y_x_z_yz, g_y_x_z_zz, g_y_x_zzz_xx, g_y_x_zzz_xy, g_y_x_zzz_xz, g_y_x_zzz_yy, g_y_x_zzz_yz, g_y_x_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_x_zz_xx[i] = -4.0 * g_y_x_z_xx[i] * a_exp + 4.0 * g_y_x_zzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_zz_xy[i] = -4.0 * g_y_x_z_xy[i] * a_exp + 4.0 * g_y_x_zzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_zz_xz[i] = -4.0 * g_y_x_z_xz[i] * a_exp + 4.0 * g_y_x_zzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_zz_yy[i] = -4.0 * g_y_x_z_yy[i] * a_exp + 4.0 * g_y_x_zzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_zz_yz[i] = -4.0 * g_y_x_z_yz[i] * a_exp + 4.0 * g_y_x_zzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_zz_zz[i] = -4.0 * g_y_x_z_zz[i] * a_exp + 4.0 * g_y_x_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (576-582)

    #pragma omp simd aligned(g_y_0_z_0_0_y_xx_xx, g_y_0_z_0_0_y_xx_xy, g_y_0_z_0_0_y_xx_xz, g_y_0_z_0_0_y_xx_yy, g_y_0_z_0_0_y_xx_yz, g_y_0_z_0_0_y_xx_zz, g_y_y_xxz_xx, g_y_y_xxz_xy, g_y_y_xxz_xz, g_y_y_xxz_yy, g_y_y_xxz_yz, g_y_y_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_y_xx_xx[i] = 4.0 * g_y_y_xxz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_xx_xy[i] = 4.0 * g_y_y_xxz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_xx_xz[i] = 4.0 * g_y_y_xxz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_xx_yy[i] = 4.0 * g_y_y_xxz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_xx_yz[i] = 4.0 * g_y_y_xxz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_xx_zz[i] = 4.0 * g_y_y_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (582-588)

    #pragma omp simd aligned(g_y_0_z_0_0_y_xy_xx, g_y_0_z_0_0_y_xy_xy, g_y_0_z_0_0_y_xy_xz, g_y_0_z_0_0_y_xy_yy, g_y_0_z_0_0_y_xy_yz, g_y_0_z_0_0_y_xy_zz, g_y_y_xyz_xx, g_y_y_xyz_xy, g_y_y_xyz_xz, g_y_y_xyz_yy, g_y_y_xyz_yz, g_y_y_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_y_xy_xx[i] = 4.0 * g_y_y_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_xy_xy[i] = 4.0 * g_y_y_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_xy_xz[i] = 4.0 * g_y_y_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_xy_yy[i] = 4.0 * g_y_y_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_xy_yz[i] = 4.0 * g_y_y_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_xy_zz[i] = 4.0 * g_y_y_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (588-594)

    #pragma omp simd aligned(g_y_0_z_0_0_y_xz_xx, g_y_0_z_0_0_y_xz_xy, g_y_0_z_0_0_y_xz_xz, g_y_0_z_0_0_y_xz_yy, g_y_0_z_0_0_y_xz_yz, g_y_0_z_0_0_y_xz_zz, g_y_y_x_xx, g_y_y_x_xy, g_y_y_x_xz, g_y_y_x_yy, g_y_y_x_yz, g_y_y_x_zz, g_y_y_xzz_xx, g_y_y_xzz_xy, g_y_y_xzz_xz, g_y_y_xzz_yy, g_y_y_xzz_yz, g_y_y_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_y_xz_xx[i] = -2.0 * g_y_y_x_xx[i] * a_exp + 4.0 * g_y_y_xzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_xz_xy[i] = -2.0 * g_y_y_x_xy[i] * a_exp + 4.0 * g_y_y_xzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_xz_xz[i] = -2.0 * g_y_y_x_xz[i] * a_exp + 4.0 * g_y_y_xzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_xz_yy[i] = -2.0 * g_y_y_x_yy[i] * a_exp + 4.0 * g_y_y_xzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_xz_yz[i] = -2.0 * g_y_y_x_yz[i] * a_exp + 4.0 * g_y_y_xzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_xz_zz[i] = -2.0 * g_y_y_x_zz[i] * a_exp + 4.0 * g_y_y_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (594-600)

    #pragma omp simd aligned(g_y_0_z_0_0_y_yy_xx, g_y_0_z_0_0_y_yy_xy, g_y_0_z_0_0_y_yy_xz, g_y_0_z_0_0_y_yy_yy, g_y_0_z_0_0_y_yy_yz, g_y_0_z_0_0_y_yy_zz, g_y_y_yyz_xx, g_y_y_yyz_xy, g_y_y_yyz_xz, g_y_y_yyz_yy, g_y_y_yyz_yz, g_y_y_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_y_yy_xx[i] = 4.0 * g_y_y_yyz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_yy_xy[i] = 4.0 * g_y_y_yyz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_yy_xz[i] = 4.0 * g_y_y_yyz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_yy_yy[i] = 4.0 * g_y_y_yyz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_yy_yz[i] = 4.0 * g_y_y_yyz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_yy_zz[i] = 4.0 * g_y_y_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (600-606)

    #pragma omp simd aligned(g_y_0_z_0_0_y_yz_xx, g_y_0_z_0_0_y_yz_xy, g_y_0_z_0_0_y_yz_xz, g_y_0_z_0_0_y_yz_yy, g_y_0_z_0_0_y_yz_yz, g_y_0_z_0_0_y_yz_zz, g_y_y_y_xx, g_y_y_y_xy, g_y_y_y_xz, g_y_y_y_yy, g_y_y_y_yz, g_y_y_y_zz, g_y_y_yzz_xx, g_y_y_yzz_xy, g_y_y_yzz_xz, g_y_y_yzz_yy, g_y_y_yzz_yz, g_y_y_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_y_yz_xx[i] = -2.0 * g_y_y_y_xx[i] * a_exp + 4.0 * g_y_y_yzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_yz_xy[i] = -2.0 * g_y_y_y_xy[i] * a_exp + 4.0 * g_y_y_yzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_yz_xz[i] = -2.0 * g_y_y_y_xz[i] * a_exp + 4.0 * g_y_y_yzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_yz_yy[i] = -2.0 * g_y_y_y_yy[i] * a_exp + 4.0 * g_y_y_yzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_yz_yz[i] = -2.0 * g_y_y_y_yz[i] * a_exp + 4.0 * g_y_y_yzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_yz_zz[i] = -2.0 * g_y_y_y_zz[i] * a_exp + 4.0 * g_y_y_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (606-612)

    #pragma omp simd aligned(g_y_0_z_0_0_y_zz_xx, g_y_0_z_0_0_y_zz_xy, g_y_0_z_0_0_y_zz_xz, g_y_0_z_0_0_y_zz_yy, g_y_0_z_0_0_y_zz_yz, g_y_0_z_0_0_y_zz_zz, g_y_y_z_xx, g_y_y_z_xy, g_y_y_z_xz, g_y_y_z_yy, g_y_y_z_yz, g_y_y_z_zz, g_y_y_zzz_xx, g_y_y_zzz_xy, g_y_y_zzz_xz, g_y_y_zzz_yy, g_y_y_zzz_yz, g_y_y_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_y_zz_xx[i] = -4.0 * g_y_y_z_xx[i] * a_exp + 4.0 * g_y_y_zzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_zz_xy[i] = -4.0 * g_y_y_z_xy[i] * a_exp + 4.0 * g_y_y_zzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_zz_xz[i] = -4.0 * g_y_y_z_xz[i] * a_exp + 4.0 * g_y_y_zzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_zz_yy[i] = -4.0 * g_y_y_z_yy[i] * a_exp + 4.0 * g_y_y_zzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_zz_yz[i] = -4.0 * g_y_y_z_yz[i] * a_exp + 4.0 * g_y_y_zzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_zz_zz[i] = -4.0 * g_y_y_z_zz[i] * a_exp + 4.0 * g_y_y_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (612-618)

    #pragma omp simd aligned(g_y_0_z_0_0_z_xx_xx, g_y_0_z_0_0_z_xx_xy, g_y_0_z_0_0_z_xx_xz, g_y_0_z_0_0_z_xx_yy, g_y_0_z_0_0_z_xx_yz, g_y_0_z_0_0_z_xx_zz, g_y_z_xxz_xx, g_y_z_xxz_xy, g_y_z_xxz_xz, g_y_z_xxz_yy, g_y_z_xxz_yz, g_y_z_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_z_xx_xx[i] = 4.0 * g_y_z_xxz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_xx_xy[i] = 4.0 * g_y_z_xxz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_xx_xz[i] = 4.0 * g_y_z_xxz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_xx_yy[i] = 4.0 * g_y_z_xxz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_xx_yz[i] = 4.0 * g_y_z_xxz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_xx_zz[i] = 4.0 * g_y_z_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (618-624)

    #pragma omp simd aligned(g_y_0_z_0_0_z_xy_xx, g_y_0_z_0_0_z_xy_xy, g_y_0_z_0_0_z_xy_xz, g_y_0_z_0_0_z_xy_yy, g_y_0_z_0_0_z_xy_yz, g_y_0_z_0_0_z_xy_zz, g_y_z_xyz_xx, g_y_z_xyz_xy, g_y_z_xyz_xz, g_y_z_xyz_yy, g_y_z_xyz_yz, g_y_z_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_z_xy_xx[i] = 4.0 * g_y_z_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_xy_xy[i] = 4.0 * g_y_z_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_xy_xz[i] = 4.0 * g_y_z_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_xy_yy[i] = 4.0 * g_y_z_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_xy_yz[i] = 4.0 * g_y_z_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_xy_zz[i] = 4.0 * g_y_z_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (624-630)

    #pragma omp simd aligned(g_y_0_z_0_0_z_xz_xx, g_y_0_z_0_0_z_xz_xy, g_y_0_z_0_0_z_xz_xz, g_y_0_z_0_0_z_xz_yy, g_y_0_z_0_0_z_xz_yz, g_y_0_z_0_0_z_xz_zz, g_y_z_x_xx, g_y_z_x_xy, g_y_z_x_xz, g_y_z_x_yy, g_y_z_x_yz, g_y_z_x_zz, g_y_z_xzz_xx, g_y_z_xzz_xy, g_y_z_xzz_xz, g_y_z_xzz_yy, g_y_z_xzz_yz, g_y_z_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_z_xz_xx[i] = -2.0 * g_y_z_x_xx[i] * a_exp + 4.0 * g_y_z_xzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_xz_xy[i] = -2.0 * g_y_z_x_xy[i] * a_exp + 4.0 * g_y_z_xzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_xz_xz[i] = -2.0 * g_y_z_x_xz[i] * a_exp + 4.0 * g_y_z_xzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_xz_yy[i] = -2.0 * g_y_z_x_yy[i] * a_exp + 4.0 * g_y_z_xzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_xz_yz[i] = -2.0 * g_y_z_x_yz[i] * a_exp + 4.0 * g_y_z_xzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_xz_zz[i] = -2.0 * g_y_z_x_zz[i] * a_exp + 4.0 * g_y_z_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (630-636)

    #pragma omp simd aligned(g_y_0_z_0_0_z_yy_xx, g_y_0_z_0_0_z_yy_xy, g_y_0_z_0_0_z_yy_xz, g_y_0_z_0_0_z_yy_yy, g_y_0_z_0_0_z_yy_yz, g_y_0_z_0_0_z_yy_zz, g_y_z_yyz_xx, g_y_z_yyz_xy, g_y_z_yyz_xz, g_y_z_yyz_yy, g_y_z_yyz_yz, g_y_z_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_z_yy_xx[i] = 4.0 * g_y_z_yyz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_yy_xy[i] = 4.0 * g_y_z_yyz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_yy_xz[i] = 4.0 * g_y_z_yyz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_yy_yy[i] = 4.0 * g_y_z_yyz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_yy_yz[i] = 4.0 * g_y_z_yyz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_yy_zz[i] = 4.0 * g_y_z_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (636-642)

    #pragma omp simd aligned(g_y_0_z_0_0_z_yz_xx, g_y_0_z_0_0_z_yz_xy, g_y_0_z_0_0_z_yz_xz, g_y_0_z_0_0_z_yz_yy, g_y_0_z_0_0_z_yz_yz, g_y_0_z_0_0_z_yz_zz, g_y_z_y_xx, g_y_z_y_xy, g_y_z_y_xz, g_y_z_y_yy, g_y_z_y_yz, g_y_z_y_zz, g_y_z_yzz_xx, g_y_z_yzz_xy, g_y_z_yzz_xz, g_y_z_yzz_yy, g_y_z_yzz_yz, g_y_z_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_z_yz_xx[i] = -2.0 * g_y_z_y_xx[i] * a_exp + 4.0 * g_y_z_yzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_yz_xy[i] = -2.0 * g_y_z_y_xy[i] * a_exp + 4.0 * g_y_z_yzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_yz_xz[i] = -2.0 * g_y_z_y_xz[i] * a_exp + 4.0 * g_y_z_yzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_yz_yy[i] = -2.0 * g_y_z_y_yy[i] * a_exp + 4.0 * g_y_z_yzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_yz_yz[i] = -2.0 * g_y_z_y_yz[i] * a_exp + 4.0 * g_y_z_yzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_yz_zz[i] = -2.0 * g_y_z_y_zz[i] * a_exp + 4.0 * g_y_z_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (642-648)

    #pragma omp simd aligned(g_y_0_z_0_0_z_zz_xx, g_y_0_z_0_0_z_zz_xy, g_y_0_z_0_0_z_zz_xz, g_y_0_z_0_0_z_zz_yy, g_y_0_z_0_0_z_zz_yz, g_y_0_z_0_0_z_zz_zz, g_y_z_z_xx, g_y_z_z_xy, g_y_z_z_xz, g_y_z_z_yy, g_y_z_z_yz, g_y_z_z_zz, g_y_z_zzz_xx, g_y_z_zzz_xy, g_y_z_zzz_xz, g_y_z_zzz_yy, g_y_z_zzz_yz, g_y_z_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_z_zz_xx[i] = -4.0 * g_y_z_z_xx[i] * a_exp + 4.0 * g_y_z_zzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_zz_xy[i] = -4.0 * g_y_z_z_xy[i] * a_exp + 4.0 * g_y_z_zzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_zz_xz[i] = -4.0 * g_y_z_z_xz[i] * a_exp + 4.0 * g_y_z_zzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_zz_yy[i] = -4.0 * g_y_z_z_yy[i] * a_exp + 4.0 * g_y_z_zzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_zz_yz[i] = -4.0 * g_y_z_z_yz[i] * a_exp + 4.0 * g_y_z_zzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_zz_zz[i] = -4.0 * g_y_z_z_zz[i] * a_exp + 4.0 * g_y_z_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (648-654)

    #pragma omp simd aligned(g_z_0_x_0_0_x_xx_xx, g_z_0_x_0_0_x_xx_xy, g_z_0_x_0_0_x_xx_xz, g_z_0_x_0_0_x_xx_yy, g_z_0_x_0_0_x_xx_yz, g_z_0_x_0_0_x_xx_zz, g_z_x_x_xx, g_z_x_x_xy, g_z_x_x_xz, g_z_x_x_yy, g_z_x_x_yz, g_z_x_x_zz, g_z_x_xxx_xx, g_z_x_xxx_xy, g_z_x_xxx_xz, g_z_x_xxx_yy, g_z_x_xxx_yz, g_z_x_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_x_xx_xx[i] = -4.0 * g_z_x_x_xx[i] * a_exp + 4.0 * g_z_x_xxx_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_xx_xy[i] = -4.0 * g_z_x_x_xy[i] * a_exp + 4.0 * g_z_x_xxx_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_xx_xz[i] = -4.0 * g_z_x_x_xz[i] * a_exp + 4.0 * g_z_x_xxx_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_xx_yy[i] = -4.0 * g_z_x_x_yy[i] * a_exp + 4.0 * g_z_x_xxx_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_xx_yz[i] = -4.0 * g_z_x_x_yz[i] * a_exp + 4.0 * g_z_x_xxx_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_xx_zz[i] = -4.0 * g_z_x_x_zz[i] * a_exp + 4.0 * g_z_x_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (654-660)

    #pragma omp simd aligned(g_z_0_x_0_0_x_xy_xx, g_z_0_x_0_0_x_xy_xy, g_z_0_x_0_0_x_xy_xz, g_z_0_x_0_0_x_xy_yy, g_z_0_x_0_0_x_xy_yz, g_z_0_x_0_0_x_xy_zz, g_z_x_xxy_xx, g_z_x_xxy_xy, g_z_x_xxy_xz, g_z_x_xxy_yy, g_z_x_xxy_yz, g_z_x_xxy_zz, g_z_x_y_xx, g_z_x_y_xy, g_z_x_y_xz, g_z_x_y_yy, g_z_x_y_yz, g_z_x_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_x_xy_xx[i] = -2.0 * g_z_x_y_xx[i] * a_exp + 4.0 * g_z_x_xxy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_xy_xy[i] = -2.0 * g_z_x_y_xy[i] * a_exp + 4.0 * g_z_x_xxy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_xy_xz[i] = -2.0 * g_z_x_y_xz[i] * a_exp + 4.0 * g_z_x_xxy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_xy_yy[i] = -2.0 * g_z_x_y_yy[i] * a_exp + 4.0 * g_z_x_xxy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_xy_yz[i] = -2.0 * g_z_x_y_yz[i] * a_exp + 4.0 * g_z_x_xxy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_xy_zz[i] = -2.0 * g_z_x_y_zz[i] * a_exp + 4.0 * g_z_x_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (660-666)

    #pragma omp simd aligned(g_z_0_x_0_0_x_xz_xx, g_z_0_x_0_0_x_xz_xy, g_z_0_x_0_0_x_xz_xz, g_z_0_x_0_0_x_xz_yy, g_z_0_x_0_0_x_xz_yz, g_z_0_x_0_0_x_xz_zz, g_z_x_xxz_xx, g_z_x_xxz_xy, g_z_x_xxz_xz, g_z_x_xxz_yy, g_z_x_xxz_yz, g_z_x_xxz_zz, g_z_x_z_xx, g_z_x_z_xy, g_z_x_z_xz, g_z_x_z_yy, g_z_x_z_yz, g_z_x_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_x_xz_xx[i] = -2.0 * g_z_x_z_xx[i] * a_exp + 4.0 * g_z_x_xxz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_xz_xy[i] = -2.0 * g_z_x_z_xy[i] * a_exp + 4.0 * g_z_x_xxz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_xz_xz[i] = -2.0 * g_z_x_z_xz[i] * a_exp + 4.0 * g_z_x_xxz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_xz_yy[i] = -2.0 * g_z_x_z_yy[i] * a_exp + 4.0 * g_z_x_xxz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_xz_yz[i] = -2.0 * g_z_x_z_yz[i] * a_exp + 4.0 * g_z_x_xxz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_xz_zz[i] = -2.0 * g_z_x_z_zz[i] * a_exp + 4.0 * g_z_x_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (666-672)

    #pragma omp simd aligned(g_z_0_x_0_0_x_yy_xx, g_z_0_x_0_0_x_yy_xy, g_z_0_x_0_0_x_yy_xz, g_z_0_x_0_0_x_yy_yy, g_z_0_x_0_0_x_yy_yz, g_z_0_x_0_0_x_yy_zz, g_z_x_xyy_xx, g_z_x_xyy_xy, g_z_x_xyy_xz, g_z_x_xyy_yy, g_z_x_xyy_yz, g_z_x_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_x_yy_xx[i] = 4.0 * g_z_x_xyy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_yy_xy[i] = 4.0 * g_z_x_xyy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_yy_xz[i] = 4.0 * g_z_x_xyy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_yy_yy[i] = 4.0 * g_z_x_xyy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_yy_yz[i] = 4.0 * g_z_x_xyy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_yy_zz[i] = 4.0 * g_z_x_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (672-678)

    #pragma omp simd aligned(g_z_0_x_0_0_x_yz_xx, g_z_0_x_0_0_x_yz_xy, g_z_0_x_0_0_x_yz_xz, g_z_0_x_0_0_x_yz_yy, g_z_0_x_0_0_x_yz_yz, g_z_0_x_0_0_x_yz_zz, g_z_x_xyz_xx, g_z_x_xyz_xy, g_z_x_xyz_xz, g_z_x_xyz_yy, g_z_x_xyz_yz, g_z_x_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_x_yz_xx[i] = 4.0 * g_z_x_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_yz_xy[i] = 4.0 * g_z_x_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_yz_xz[i] = 4.0 * g_z_x_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_yz_yy[i] = 4.0 * g_z_x_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_yz_yz[i] = 4.0 * g_z_x_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_yz_zz[i] = 4.0 * g_z_x_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (678-684)

    #pragma omp simd aligned(g_z_0_x_0_0_x_zz_xx, g_z_0_x_0_0_x_zz_xy, g_z_0_x_0_0_x_zz_xz, g_z_0_x_0_0_x_zz_yy, g_z_0_x_0_0_x_zz_yz, g_z_0_x_0_0_x_zz_zz, g_z_x_xzz_xx, g_z_x_xzz_xy, g_z_x_xzz_xz, g_z_x_xzz_yy, g_z_x_xzz_yz, g_z_x_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_x_zz_xx[i] = 4.0 * g_z_x_xzz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_zz_xy[i] = 4.0 * g_z_x_xzz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_zz_xz[i] = 4.0 * g_z_x_xzz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_zz_yy[i] = 4.0 * g_z_x_xzz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_zz_yz[i] = 4.0 * g_z_x_xzz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_zz_zz[i] = 4.0 * g_z_x_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (684-690)

    #pragma omp simd aligned(g_z_0_x_0_0_y_xx_xx, g_z_0_x_0_0_y_xx_xy, g_z_0_x_0_0_y_xx_xz, g_z_0_x_0_0_y_xx_yy, g_z_0_x_0_0_y_xx_yz, g_z_0_x_0_0_y_xx_zz, g_z_y_x_xx, g_z_y_x_xy, g_z_y_x_xz, g_z_y_x_yy, g_z_y_x_yz, g_z_y_x_zz, g_z_y_xxx_xx, g_z_y_xxx_xy, g_z_y_xxx_xz, g_z_y_xxx_yy, g_z_y_xxx_yz, g_z_y_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_y_xx_xx[i] = -4.0 * g_z_y_x_xx[i] * a_exp + 4.0 * g_z_y_xxx_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_xx_xy[i] = -4.0 * g_z_y_x_xy[i] * a_exp + 4.0 * g_z_y_xxx_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_xx_xz[i] = -4.0 * g_z_y_x_xz[i] * a_exp + 4.0 * g_z_y_xxx_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_xx_yy[i] = -4.0 * g_z_y_x_yy[i] * a_exp + 4.0 * g_z_y_xxx_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_xx_yz[i] = -4.0 * g_z_y_x_yz[i] * a_exp + 4.0 * g_z_y_xxx_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_xx_zz[i] = -4.0 * g_z_y_x_zz[i] * a_exp + 4.0 * g_z_y_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (690-696)

    #pragma omp simd aligned(g_z_0_x_0_0_y_xy_xx, g_z_0_x_0_0_y_xy_xy, g_z_0_x_0_0_y_xy_xz, g_z_0_x_0_0_y_xy_yy, g_z_0_x_0_0_y_xy_yz, g_z_0_x_0_0_y_xy_zz, g_z_y_xxy_xx, g_z_y_xxy_xy, g_z_y_xxy_xz, g_z_y_xxy_yy, g_z_y_xxy_yz, g_z_y_xxy_zz, g_z_y_y_xx, g_z_y_y_xy, g_z_y_y_xz, g_z_y_y_yy, g_z_y_y_yz, g_z_y_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_y_xy_xx[i] = -2.0 * g_z_y_y_xx[i] * a_exp + 4.0 * g_z_y_xxy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_xy_xy[i] = -2.0 * g_z_y_y_xy[i] * a_exp + 4.0 * g_z_y_xxy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_xy_xz[i] = -2.0 * g_z_y_y_xz[i] * a_exp + 4.0 * g_z_y_xxy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_xy_yy[i] = -2.0 * g_z_y_y_yy[i] * a_exp + 4.0 * g_z_y_xxy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_xy_yz[i] = -2.0 * g_z_y_y_yz[i] * a_exp + 4.0 * g_z_y_xxy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_xy_zz[i] = -2.0 * g_z_y_y_zz[i] * a_exp + 4.0 * g_z_y_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (696-702)

    #pragma omp simd aligned(g_z_0_x_0_0_y_xz_xx, g_z_0_x_0_0_y_xz_xy, g_z_0_x_0_0_y_xz_xz, g_z_0_x_0_0_y_xz_yy, g_z_0_x_0_0_y_xz_yz, g_z_0_x_0_0_y_xz_zz, g_z_y_xxz_xx, g_z_y_xxz_xy, g_z_y_xxz_xz, g_z_y_xxz_yy, g_z_y_xxz_yz, g_z_y_xxz_zz, g_z_y_z_xx, g_z_y_z_xy, g_z_y_z_xz, g_z_y_z_yy, g_z_y_z_yz, g_z_y_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_y_xz_xx[i] = -2.0 * g_z_y_z_xx[i] * a_exp + 4.0 * g_z_y_xxz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_xz_xy[i] = -2.0 * g_z_y_z_xy[i] * a_exp + 4.0 * g_z_y_xxz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_xz_xz[i] = -2.0 * g_z_y_z_xz[i] * a_exp + 4.0 * g_z_y_xxz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_xz_yy[i] = -2.0 * g_z_y_z_yy[i] * a_exp + 4.0 * g_z_y_xxz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_xz_yz[i] = -2.0 * g_z_y_z_yz[i] * a_exp + 4.0 * g_z_y_xxz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_xz_zz[i] = -2.0 * g_z_y_z_zz[i] * a_exp + 4.0 * g_z_y_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (702-708)

    #pragma omp simd aligned(g_z_0_x_0_0_y_yy_xx, g_z_0_x_0_0_y_yy_xy, g_z_0_x_0_0_y_yy_xz, g_z_0_x_0_0_y_yy_yy, g_z_0_x_0_0_y_yy_yz, g_z_0_x_0_0_y_yy_zz, g_z_y_xyy_xx, g_z_y_xyy_xy, g_z_y_xyy_xz, g_z_y_xyy_yy, g_z_y_xyy_yz, g_z_y_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_y_yy_xx[i] = 4.0 * g_z_y_xyy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_yy_xy[i] = 4.0 * g_z_y_xyy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_yy_xz[i] = 4.0 * g_z_y_xyy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_yy_yy[i] = 4.0 * g_z_y_xyy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_yy_yz[i] = 4.0 * g_z_y_xyy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_yy_zz[i] = 4.0 * g_z_y_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (708-714)

    #pragma omp simd aligned(g_z_0_x_0_0_y_yz_xx, g_z_0_x_0_0_y_yz_xy, g_z_0_x_0_0_y_yz_xz, g_z_0_x_0_0_y_yz_yy, g_z_0_x_0_0_y_yz_yz, g_z_0_x_0_0_y_yz_zz, g_z_y_xyz_xx, g_z_y_xyz_xy, g_z_y_xyz_xz, g_z_y_xyz_yy, g_z_y_xyz_yz, g_z_y_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_y_yz_xx[i] = 4.0 * g_z_y_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_yz_xy[i] = 4.0 * g_z_y_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_yz_xz[i] = 4.0 * g_z_y_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_yz_yy[i] = 4.0 * g_z_y_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_yz_yz[i] = 4.0 * g_z_y_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_yz_zz[i] = 4.0 * g_z_y_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (714-720)

    #pragma omp simd aligned(g_z_0_x_0_0_y_zz_xx, g_z_0_x_0_0_y_zz_xy, g_z_0_x_0_0_y_zz_xz, g_z_0_x_0_0_y_zz_yy, g_z_0_x_0_0_y_zz_yz, g_z_0_x_0_0_y_zz_zz, g_z_y_xzz_xx, g_z_y_xzz_xy, g_z_y_xzz_xz, g_z_y_xzz_yy, g_z_y_xzz_yz, g_z_y_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_y_zz_xx[i] = 4.0 * g_z_y_xzz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_zz_xy[i] = 4.0 * g_z_y_xzz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_zz_xz[i] = 4.0 * g_z_y_xzz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_zz_yy[i] = 4.0 * g_z_y_xzz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_zz_yz[i] = 4.0 * g_z_y_xzz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_zz_zz[i] = 4.0 * g_z_y_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (720-726)

    #pragma omp simd aligned(g_z_0_x_0_0_z_xx_xx, g_z_0_x_0_0_z_xx_xy, g_z_0_x_0_0_z_xx_xz, g_z_0_x_0_0_z_xx_yy, g_z_0_x_0_0_z_xx_yz, g_z_0_x_0_0_z_xx_zz, g_z_z_x_xx, g_z_z_x_xy, g_z_z_x_xz, g_z_z_x_yy, g_z_z_x_yz, g_z_z_x_zz, g_z_z_xxx_xx, g_z_z_xxx_xy, g_z_z_xxx_xz, g_z_z_xxx_yy, g_z_z_xxx_yz, g_z_z_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_z_xx_xx[i] = -4.0 * g_z_z_x_xx[i] * a_exp + 4.0 * g_z_z_xxx_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_xx_xy[i] = -4.0 * g_z_z_x_xy[i] * a_exp + 4.0 * g_z_z_xxx_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_xx_xz[i] = -4.0 * g_z_z_x_xz[i] * a_exp + 4.0 * g_z_z_xxx_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_xx_yy[i] = -4.0 * g_z_z_x_yy[i] * a_exp + 4.0 * g_z_z_xxx_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_xx_yz[i] = -4.0 * g_z_z_x_yz[i] * a_exp + 4.0 * g_z_z_xxx_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_xx_zz[i] = -4.0 * g_z_z_x_zz[i] * a_exp + 4.0 * g_z_z_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (726-732)

    #pragma omp simd aligned(g_z_0_x_0_0_z_xy_xx, g_z_0_x_0_0_z_xy_xy, g_z_0_x_0_0_z_xy_xz, g_z_0_x_0_0_z_xy_yy, g_z_0_x_0_0_z_xy_yz, g_z_0_x_0_0_z_xy_zz, g_z_z_xxy_xx, g_z_z_xxy_xy, g_z_z_xxy_xz, g_z_z_xxy_yy, g_z_z_xxy_yz, g_z_z_xxy_zz, g_z_z_y_xx, g_z_z_y_xy, g_z_z_y_xz, g_z_z_y_yy, g_z_z_y_yz, g_z_z_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_z_xy_xx[i] = -2.0 * g_z_z_y_xx[i] * a_exp + 4.0 * g_z_z_xxy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_xy_xy[i] = -2.0 * g_z_z_y_xy[i] * a_exp + 4.0 * g_z_z_xxy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_xy_xz[i] = -2.0 * g_z_z_y_xz[i] * a_exp + 4.0 * g_z_z_xxy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_xy_yy[i] = -2.0 * g_z_z_y_yy[i] * a_exp + 4.0 * g_z_z_xxy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_xy_yz[i] = -2.0 * g_z_z_y_yz[i] * a_exp + 4.0 * g_z_z_xxy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_xy_zz[i] = -2.0 * g_z_z_y_zz[i] * a_exp + 4.0 * g_z_z_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (732-738)

    #pragma omp simd aligned(g_z_0_x_0_0_z_xz_xx, g_z_0_x_0_0_z_xz_xy, g_z_0_x_0_0_z_xz_xz, g_z_0_x_0_0_z_xz_yy, g_z_0_x_0_0_z_xz_yz, g_z_0_x_0_0_z_xz_zz, g_z_z_xxz_xx, g_z_z_xxz_xy, g_z_z_xxz_xz, g_z_z_xxz_yy, g_z_z_xxz_yz, g_z_z_xxz_zz, g_z_z_z_xx, g_z_z_z_xy, g_z_z_z_xz, g_z_z_z_yy, g_z_z_z_yz, g_z_z_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_z_xz_xx[i] = -2.0 * g_z_z_z_xx[i] * a_exp + 4.0 * g_z_z_xxz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_xz_xy[i] = -2.0 * g_z_z_z_xy[i] * a_exp + 4.0 * g_z_z_xxz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_xz_xz[i] = -2.0 * g_z_z_z_xz[i] * a_exp + 4.0 * g_z_z_xxz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_xz_yy[i] = -2.0 * g_z_z_z_yy[i] * a_exp + 4.0 * g_z_z_xxz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_xz_yz[i] = -2.0 * g_z_z_z_yz[i] * a_exp + 4.0 * g_z_z_xxz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_xz_zz[i] = -2.0 * g_z_z_z_zz[i] * a_exp + 4.0 * g_z_z_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (738-744)

    #pragma omp simd aligned(g_z_0_x_0_0_z_yy_xx, g_z_0_x_0_0_z_yy_xy, g_z_0_x_0_0_z_yy_xz, g_z_0_x_0_0_z_yy_yy, g_z_0_x_0_0_z_yy_yz, g_z_0_x_0_0_z_yy_zz, g_z_z_xyy_xx, g_z_z_xyy_xy, g_z_z_xyy_xz, g_z_z_xyy_yy, g_z_z_xyy_yz, g_z_z_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_z_yy_xx[i] = 4.0 * g_z_z_xyy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_yy_xy[i] = 4.0 * g_z_z_xyy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_yy_xz[i] = 4.0 * g_z_z_xyy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_yy_yy[i] = 4.0 * g_z_z_xyy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_yy_yz[i] = 4.0 * g_z_z_xyy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_yy_zz[i] = 4.0 * g_z_z_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (744-750)

    #pragma omp simd aligned(g_z_0_x_0_0_z_yz_xx, g_z_0_x_0_0_z_yz_xy, g_z_0_x_0_0_z_yz_xz, g_z_0_x_0_0_z_yz_yy, g_z_0_x_0_0_z_yz_yz, g_z_0_x_0_0_z_yz_zz, g_z_z_xyz_xx, g_z_z_xyz_xy, g_z_z_xyz_xz, g_z_z_xyz_yy, g_z_z_xyz_yz, g_z_z_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_z_yz_xx[i] = 4.0 * g_z_z_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_yz_xy[i] = 4.0 * g_z_z_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_yz_xz[i] = 4.0 * g_z_z_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_yz_yy[i] = 4.0 * g_z_z_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_yz_yz[i] = 4.0 * g_z_z_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_yz_zz[i] = 4.0 * g_z_z_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (750-756)

    #pragma omp simd aligned(g_z_0_x_0_0_z_zz_xx, g_z_0_x_0_0_z_zz_xy, g_z_0_x_0_0_z_zz_xz, g_z_0_x_0_0_z_zz_yy, g_z_0_x_0_0_z_zz_yz, g_z_0_x_0_0_z_zz_zz, g_z_z_xzz_xx, g_z_z_xzz_xy, g_z_z_xzz_xz, g_z_z_xzz_yy, g_z_z_xzz_yz, g_z_z_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_z_zz_xx[i] = 4.0 * g_z_z_xzz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_zz_xy[i] = 4.0 * g_z_z_xzz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_zz_xz[i] = 4.0 * g_z_z_xzz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_zz_yy[i] = 4.0 * g_z_z_xzz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_zz_yz[i] = 4.0 * g_z_z_xzz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_zz_zz[i] = 4.0 * g_z_z_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (756-762)

    #pragma omp simd aligned(g_z_0_y_0_0_x_xx_xx, g_z_0_y_0_0_x_xx_xy, g_z_0_y_0_0_x_xx_xz, g_z_0_y_0_0_x_xx_yy, g_z_0_y_0_0_x_xx_yz, g_z_0_y_0_0_x_xx_zz, g_z_x_xxy_xx, g_z_x_xxy_xy, g_z_x_xxy_xz, g_z_x_xxy_yy, g_z_x_xxy_yz, g_z_x_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_x_xx_xx[i] = 4.0 * g_z_x_xxy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_xx_xy[i] = 4.0 * g_z_x_xxy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_xx_xz[i] = 4.0 * g_z_x_xxy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_xx_yy[i] = 4.0 * g_z_x_xxy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_xx_yz[i] = 4.0 * g_z_x_xxy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_xx_zz[i] = 4.0 * g_z_x_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (762-768)

    #pragma omp simd aligned(g_z_0_y_0_0_x_xy_xx, g_z_0_y_0_0_x_xy_xy, g_z_0_y_0_0_x_xy_xz, g_z_0_y_0_0_x_xy_yy, g_z_0_y_0_0_x_xy_yz, g_z_0_y_0_0_x_xy_zz, g_z_x_x_xx, g_z_x_x_xy, g_z_x_x_xz, g_z_x_x_yy, g_z_x_x_yz, g_z_x_x_zz, g_z_x_xyy_xx, g_z_x_xyy_xy, g_z_x_xyy_xz, g_z_x_xyy_yy, g_z_x_xyy_yz, g_z_x_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_x_xy_xx[i] = -2.0 * g_z_x_x_xx[i] * a_exp + 4.0 * g_z_x_xyy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_xy_xy[i] = -2.0 * g_z_x_x_xy[i] * a_exp + 4.0 * g_z_x_xyy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_xy_xz[i] = -2.0 * g_z_x_x_xz[i] * a_exp + 4.0 * g_z_x_xyy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_xy_yy[i] = -2.0 * g_z_x_x_yy[i] * a_exp + 4.0 * g_z_x_xyy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_xy_yz[i] = -2.0 * g_z_x_x_yz[i] * a_exp + 4.0 * g_z_x_xyy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_xy_zz[i] = -2.0 * g_z_x_x_zz[i] * a_exp + 4.0 * g_z_x_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (768-774)

    #pragma omp simd aligned(g_z_0_y_0_0_x_xz_xx, g_z_0_y_0_0_x_xz_xy, g_z_0_y_0_0_x_xz_xz, g_z_0_y_0_0_x_xz_yy, g_z_0_y_0_0_x_xz_yz, g_z_0_y_0_0_x_xz_zz, g_z_x_xyz_xx, g_z_x_xyz_xy, g_z_x_xyz_xz, g_z_x_xyz_yy, g_z_x_xyz_yz, g_z_x_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_x_xz_xx[i] = 4.0 * g_z_x_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_xz_xy[i] = 4.0 * g_z_x_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_xz_xz[i] = 4.0 * g_z_x_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_xz_yy[i] = 4.0 * g_z_x_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_xz_yz[i] = 4.0 * g_z_x_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_xz_zz[i] = 4.0 * g_z_x_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (774-780)

    #pragma omp simd aligned(g_z_0_y_0_0_x_yy_xx, g_z_0_y_0_0_x_yy_xy, g_z_0_y_0_0_x_yy_xz, g_z_0_y_0_0_x_yy_yy, g_z_0_y_0_0_x_yy_yz, g_z_0_y_0_0_x_yy_zz, g_z_x_y_xx, g_z_x_y_xy, g_z_x_y_xz, g_z_x_y_yy, g_z_x_y_yz, g_z_x_y_zz, g_z_x_yyy_xx, g_z_x_yyy_xy, g_z_x_yyy_xz, g_z_x_yyy_yy, g_z_x_yyy_yz, g_z_x_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_x_yy_xx[i] = -4.0 * g_z_x_y_xx[i] * a_exp + 4.0 * g_z_x_yyy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_yy_xy[i] = -4.0 * g_z_x_y_xy[i] * a_exp + 4.0 * g_z_x_yyy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_yy_xz[i] = -4.0 * g_z_x_y_xz[i] * a_exp + 4.0 * g_z_x_yyy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_yy_yy[i] = -4.0 * g_z_x_y_yy[i] * a_exp + 4.0 * g_z_x_yyy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_yy_yz[i] = -4.0 * g_z_x_y_yz[i] * a_exp + 4.0 * g_z_x_yyy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_yy_zz[i] = -4.0 * g_z_x_y_zz[i] * a_exp + 4.0 * g_z_x_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (780-786)

    #pragma omp simd aligned(g_z_0_y_0_0_x_yz_xx, g_z_0_y_0_0_x_yz_xy, g_z_0_y_0_0_x_yz_xz, g_z_0_y_0_0_x_yz_yy, g_z_0_y_0_0_x_yz_yz, g_z_0_y_0_0_x_yz_zz, g_z_x_yyz_xx, g_z_x_yyz_xy, g_z_x_yyz_xz, g_z_x_yyz_yy, g_z_x_yyz_yz, g_z_x_yyz_zz, g_z_x_z_xx, g_z_x_z_xy, g_z_x_z_xz, g_z_x_z_yy, g_z_x_z_yz, g_z_x_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_x_yz_xx[i] = -2.0 * g_z_x_z_xx[i] * a_exp + 4.0 * g_z_x_yyz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_yz_xy[i] = -2.0 * g_z_x_z_xy[i] * a_exp + 4.0 * g_z_x_yyz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_yz_xz[i] = -2.0 * g_z_x_z_xz[i] * a_exp + 4.0 * g_z_x_yyz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_yz_yy[i] = -2.0 * g_z_x_z_yy[i] * a_exp + 4.0 * g_z_x_yyz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_yz_yz[i] = -2.0 * g_z_x_z_yz[i] * a_exp + 4.0 * g_z_x_yyz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_yz_zz[i] = -2.0 * g_z_x_z_zz[i] * a_exp + 4.0 * g_z_x_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (786-792)

    #pragma omp simd aligned(g_z_0_y_0_0_x_zz_xx, g_z_0_y_0_0_x_zz_xy, g_z_0_y_0_0_x_zz_xz, g_z_0_y_0_0_x_zz_yy, g_z_0_y_0_0_x_zz_yz, g_z_0_y_0_0_x_zz_zz, g_z_x_yzz_xx, g_z_x_yzz_xy, g_z_x_yzz_xz, g_z_x_yzz_yy, g_z_x_yzz_yz, g_z_x_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_x_zz_xx[i] = 4.0 * g_z_x_yzz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_zz_xy[i] = 4.0 * g_z_x_yzz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_zz_xz[i] = 4.0 * g_z_x_yzz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_zz_yy[i] = 4.0 * g_z_x_yzz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_zz_yz[i] = 4.0 * g_z_x_yzz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_zz_zz[i] = 4.0 * g_z_x_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (792-798)

    #pragma omp simd aligned(g_z_0_y_0_0_y_xx_xx, g_z_0_y_0_0_y_xx_xy, g_z_0_y_0_0_y_xx_xz, g_z_0_y_0_0_y_xx_yy, g_z_0_y_0_0_y_xx_yz, g_z_0_y_0_0_y_xx_zz, g_z_y_xxy_xx, g_z_y_xxy_xy, g_z_y_xxy_xz, g_z_y_xxy_yy, g_z_y_xxy_yz, g_z_y_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_y_xx_xx[i] = 4.0 * g_z_y_xxy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_xx_xy[i] = 4.0 * g_z_y_xxy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_xx_xz[i] = 4.0 * g_z_y_xxy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_xx_yy[i] = 4.0 * g_z_y_xxy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_xx_yz[i] = 4.0 * g_z_y_xxy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_xx_zz[i] = 4.0 * g_z_y_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (798-804)

    #pragma omp simd aligned(g_z_0_y_0_0_y_xy_xx, g_z_0_y_0_0_y_xy_xy, g_z_0_y_0_0_y_xy_xz, g_z_0_y_0_0_y_xy_yy, g_z_0_y_0_0_y_xy_yz, g_z_0_y_0_0_y_xy_zz, g_z_y_x_xx, g_z_y_x_xy, g_z_y_x_xz, g_z_y_x_yy, g_z_y_x_yz, g_z_y_x_zz, g_z_y_xyy_xx, g_z_y_xyy_xy, g_z_y_xyy_xz, g_z_y_xyy_yy, g_z_y_xyy_yz, g_z_y_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_y_xy_xx[i] = -2.0 * g_z_y_x_xx[i] * a_exp + 4.0 * g_z_y_xyy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_xy_xy[i] = -2.0 * g_z_y_x_xy[i] * a_exp + 4.0 * g_z_y_xyy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_xy_xz[i] = -2.0 * g_z_y_x_xz[i] * a_exp + 4.0 * g_z_y_xyy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_xy_yy[i] = -2.0 * g_z_y_x_yy[i] * a_exp + 4.0 * g_z_y_xyy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_xy_yz[i] = -2.0 * g_z_y_x_yz[i] * a_exp + 4.0 * g_z_y_xyy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_xy_zz[i] = -2.0 * g_z_y_x_zz[i] * a_exp + 4.0 * g_z_y_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (804-810)

    #pragma omp simd aligned(g_z_0_y_0_0_y_xz_xx, g_z_0_y_0_0_y_xz_xy, g_z_0_y_0_0_y_xz_xz, g_z_0_y_0_0_y_xz_yy, g_z_0_y_0_0_y_xz_yz, g_z_0_y_0_0_y_xz_zz, g_z_y_xyz_xx, g_z_y_xyz_xy, g_z_y_xyz_xz, g_z_y_xyz_yy, g_z_y_xyz_yz, g_z_y_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_y_xz_xx[i] = 4.0 * g_z_y_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_xz_xy[i] = 4.0 * g_z_y_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_xz_xz[i] = 4.0 * g_z_y_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_xz_yy[i] = 4.0 * g_z_y_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_xz_yz[i] = 4.0 * g_z_y_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_xz_zz[i] = 4.0 * g_z_y_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (810-816)

    #pragma omp simd aligned(g_z_0_y_0_0_y_yy_xx, g_z_0_y_0_0_y_yy_xy, g_z_0_y_0_0_y_yy_xz, g_z_0_y_0_0_y_yy_yy, g_z_0_y_0_0_y_yy_yz, g_z_0_y_0_0_y_yy_zz, g_z_y_y_xx, g_z_y_y_xy, g_z_y_y_xz, g_z_y_y_yy, g_z_y_y_yz, g_z_y_y_zz, g_z_y_yyy_xx, g_z_y_yyy_xy, g_z_y_yyy_xz, g_z_y_yyy_yy, g_z_y_yyy_yz, g_z_y_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_y_yy_xx[i] = -4.0 * g_z_y_y_xx[i] * a_exp + 4.0 * g_z_y_yyy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_yy_xy[i] = -4.0 * g_z_y_y_xy[i] * a_exp + 4.0 * g_z_y_yyy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_yy_xz[i] = -4.0 * g_z_y_y_xz[i] * a_exp + 4.0 * g_z_y_yyy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_yy_yy[i] = -4.0 * g_z_y_y_yy[i] * a_exp + 4.0 * g_z_y_yyy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_yy_yz[i] = -4.0 * g_z_y_y_yz[i] * a_exp + 4.0 * g_z_y_yyy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_yy_zz[i] = -4.0 * g_z_y_y_zz[i] * a_exp + 4.0 * g_z_y_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (816-822)

    #pragma omp simd aligned(g_z_0_y_0_0_y_yz_xx, g_z_0_y_0_0_y_yz_xy, g_z_0_y_0_0_y_yz_xz, g_z_0_y_0_0_y_yz_yy, g_z_0_y_0_0_y_yz_yz, g_z_0_y_0_0_y_yz_zz, g_z_y_yyz_xx, g_z_y_yyz_xy, g_z_y_yyz_xz, g_z_y_yyz_yy, g_z_y_yyz_yz, g_z_y_yyz_zz, g_z_y_z_xx, g_z_y_z_xy, g_z_y_z_xz, g_z_y_z_yy, g_z_y_z_yz, g_z_y_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_y_yz_xx[i] = -2.0 * g_z_y_z_xx[i] * a_exp + 4.0 * g_z_y_yyz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_yz_xy[i] = -2.0 * g_z_y_z_xy[i] * a_exp + 4.0 * g_z_y_yyz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_yz_xz[i] = -2.0 * g_z_y_z_xz[i] * a_exp + 4.0 * g_z_y_yyz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_yz_yy[i] = -2.0 * g_z_y_z_yy[i] * a_exp + 4.0 * g_z_y_yyz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_yz_yz[i] = -2.0 * g_z_y_z_yz[i] * a_exp + 4.0 * g_z_y_yyz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_yz_zz[i] = -2.0 * g_z_y_z_zz[i] * a_exp + 4.0 * g_z_y_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (822-828)

    #pragma omp simd aligned(g_z_0_y_0_0_y_zz_xx, g_z_0_y_0_0_y_zz_xy, g_z_0_y_0_0_y_zz_xz, g_z_0_y_0_0_y_zz_yy, g_z_0_y_0_0_y_zz_yz, g_z_0_y_0_0_y_zz_zz, g_z_y_yzz_xx, g_z_y_yzz_xy, g_z_y_yzz_xz, g_z_y_yzz_yy, g_z_y_yzz_yz, g_z_y_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_y_zz_xx[i] = 4.0 * g_z_y_yzz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_zz_xy[i] = 4.0 * g_z_y_yzz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_zz_xz[i] = 4.0 * g_z_y_yzz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_zz_yy[i] = 4.0 * g_z_y_yzz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_zz_yz[i] = 4.0 * g_z_y_yzz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_zz_zz[i] = 4.0 * g_z_y_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (828-834)

    #pragma omp simd aligned(g_z_0_y_0_0_z_xx_xx, g_z_0_y_0_0_z_xx_xy, g_z_0_y_0_0_z_xx_xz, g_z_0_y_0_0_z_xx_yy, g_z_0_y_0_0_z_xx_yz, g_z_0_y_0_0_z_xx_zz, g_z_z_xxy_xx, g_z_z_xxy_xy, g_z_z_xxy_xz, g_z_z_xxy_yy, g_z_z_xxy_yz, g_z_z_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_z_xx_xx[i] = 4.0 * g_z_z_xxy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_xx_xy[i] = 4.0 * g_z_z_xxy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_xx_xz[i] = 4.0 * g_z_z_xxy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_xx_yy[i] = 4.0 * g_z_z_xxy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_xx_yz[i] = 4.0 * g_z_z_xxy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_xx_zz[i] = 4.0 * g_z_z_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (834-840)

    #pragma omp simd aligned(g_z_0_y_0_0_z_xy_xx, g_z_0_y_0_0_z_xy_xy, g_z_0_y_0_0_z_xy_xz, g_z_0_y_0_0_z_xy_yy, g_z_0_y_0_0_z_xy_yz, g_z_0_y_0_0_z_xy_zz, g_z_z_x_xx, g_z_z_x_xy, g_z_z_x_xz, g_z_z_x_yy, g_z_z_x_yz, g_z_z_x_zz, g_z_z_xyy_xx, g_z_z_xyy_xy, g_z_z_xyy_xz, g_z_z_xyy_yy, g_z_z_xyy_yz, g_z_z_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_z_xy_xx[i] = -2.0 * g_z_z_x_xx[i] * a_exp + 4.0 * g_z_z_xyy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_xy_xy[i] = -2.0 * g_z_z_x_xy[i] * a_exp + 4.0 * g_z_z_xyy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_xy_xz[i] = -2.0 * g_z_z_x_xz[i] * a_exp + 4.0 * g_z_z_xyy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_xy_yy[i] = -2.0 * g_z_z_x_yy[i] * a_exp + 4.0 * g_z_z_xyy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_xy_yz[i] = -2.0 * g_z_z_x_yz[i] * a_exp + 4.0 * g_z_z_xyy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_xy_zz[i] = -2.0 * g_z_z_x_zz[i] * a_exp + 4.0 * g_z_z_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (840-846)

    #pragma omp simd aligned(g_z_0_y_0_0_z_xz_xx, g_z_0_y_0_0_z_xz_xy, g_z_0_y_0_0_z_xz_xz, g_z_0_y_0_0_z_xz_yy, g_z_0_y_0_0_z_xz_yz, g_z_0_y_0_0_z_xz_zz, g_z_z_xyz_xx, g_z_z_xyz_xy, g_z_z_xyz_xz, g_z_z_xyz_yy, g_z_z_xyz_yz, g_z_z_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_z_xz_xx[i] = 4.0 * g_z_z_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_xz_xy[i] = 4.0 * g_z_z_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_xz_xz[i] = 4.0 * g_z_z_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_xz_yy[i] = 4.0 * g_z_z_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_xz_yz[i] = 4.0 * g_z_z_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_xz_zz[i] = 4.0 * g_z_z_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (846-852)

    #pragma omp simd aligned(g_z_0_y_0_0_z_yy_xx, g_z_0_y_0_0_z_yy_xy, g_z_0_y_0_0_z_yy_xz, g_z_0_y_0_0_z_yy_yy, g_z_0_y_0_0_z_yy_yz, g_z_0_y_0_0_z_yy_zz, g_z_z_y_xx, g_z_z_y_xy, g_z_z_y_xz, g_z_z_y_yy, g_z_z_y_yz, g_z_z_y_zz, g_z_z_yyy_xx, g_z_z_yyy_xy, g_z_z_yyy_xz, g_z_z_yyy_yy, g_z_z_yyy_yz, g_z_z_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_z_yy_xx[i] = -4.0 * g_z_z_y_xx[i] * a_exp + 4.0 * g_z_z_yyy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_yy_xy[i] = -4.0 * g_z_z_y_xy[i] * a_exp + 4.0 * g_z_z_yyy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_yy_xz[i] = -4.0 * g_z_z_y_xz[i] * a_exp + 4.0 * g_z_z_yyy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_yy_yy[i] = -4.0 * g_z_z_y_yy[i] * a_exp + 4.0 * g_z_z_yyy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_yy_yz[i] = -4.0 * g_z_z_y_yz[i] * a_exp + 4.0 * g_z_z_yyy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_yy_zz[i] = -4.0 * g_z_z_y_zz[i] * a_exp + 4.0 * g_z_z_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (852-858)

    #pragma omp simd aligned(g_z_0_y_0_0_z_yz_xx, g_z_0_y_0_0_z_yz_xy, g_z_0_y_0_0_z_yz_xz, g_z_0_y_0_0_z_yz_yy, g_z_0_y_0_0_z_yz_yz, g_z_0_y_0_0_z_yz_zz, g_z_z_yyz_xx, g_z_z_yyz_xy, g_z_z_yyz_xz, g_z_z_yyz_yy, g_z_z_yyz_yz, g_z_z_yyz_zz, g_z_z_z_xx, g_z_z_z_xy, g_z_z_z_xz, g_z_z_z_yy, g_z_z_z_yz, g_z_z_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_z_yz_xx[i] = -2.0 * g_z_z_z_xx[i] * a_exp + 4.0 * g_z_z_yyz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_yz_xy[i] = -2.0 * g_z_z_z_xy[i] * a_exp + 4.0 * g_z_z_yyz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_yz_xz[i] = -2.0 * g_z_z_z_xz[i] * a_exp + 4.0 * g_z_z_yyz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_yz_yy[i] = -2.0 * g_z_z_z_yy[i] * a_exp + 4.0 * g_z_z_yyz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_yz_yz[i] = -2.0 * g_z_z_z_yz[i] * a_exp + 4.0 * g_z_z_yyz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_yz_zz[i] = -2.0 * g_z_z_z_zz[i] * a_exp + 4.0 * g_z_z_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (858-864)

    #pragma omp simd aligned(g_z_0_y_0_0_z_zz_xx, g_z_0_y_0_0_z_zz_xy, g_z_0_y_0_0_z_zz_xz, g_z_0_y_0_0_z_zz_yy, g_z_0_y_0_0_z_zz_yz, g_z_0_y_0_0_z_zz_zz, g_z_z_yzz_xx, g_z_z_yzz_xy, g_z_z_yzz_xz, g_z_z_yzz_yy, g_z_z_yzz_yz, g_z_z_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_z_zz_xx[i] = 4.0 * g_z_z_yzz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_zz_xy[i] = 4.0 * g_z_z_yzz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_zz_xz[i] = 4.0 * g_z_z_yzz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_zz_yy[i] = 4.0 * g_z_z_yzz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_zz_yz[i] = 4.0 * g_z_z_yzz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_zz_zz[i] = 4.0 * g_z_z_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (864-870)

    #pragma omp simd aligned(g_z_0_z_0_0_x_xx_xx, g_z_0_z_0_0_x_xx_xy, g_z_0_z_0_0_x_xx_xz, g_z_0_z_0_0_x_xx_yy, g_z_0_z_0_0_x_xx_yz, g_z_0_z_0_0_x_xx_zz, g_z_x_xxz_xx, g_z_x_xxz_xy, g_z_x_xxz_xz, g_z_x_xxz_yy, g_z_x_xxz_yz, g_z_x_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_x_xx_xx[i] = 4.0 * g_z_x_xxz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_xx_xy[i] = 4.0 * g_z_x_xxz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_xx_xz[i] = 4.0 * g_z_x_xxz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_xx_yy[i] = 4.0 * g_z_x_xxz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_xx_yz[i] = 4.0 * g_z_x_xxz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_xx_zz[i] = 4.0 * g_z_x_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (870-876)

    #pragma omp simd aligned(g_z_0_z_0_0_x_xy_xx, g_z_0_z_0_0_x_xy_xy, g_z_0_z_0_0_x_xy_xz, g_z_0_z_0_0_x_xy_yy, g_z_0_z_0_0_x_xy_yz, g_z_0_z_0_0_x_xy_zz, g_z_x_xyz_xx, g_z_x_xyz_xy, g_z_x_xyz_xz, g_z_x_xyz_yy, g_z_x_xyz_yz, g_z_x_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_x_xy_xx[i] = 4.0 * g_z_x_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_xy_xy[i] = 4.0 * g_z_x_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_xy_xz[i] = 4.0 * g_z_x_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_xy_yy[i] = 4.0 * g_z_x_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_xy_yz[i] = 4.0 * g_z_x_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_xy_zz[i] = 4.0 * g_z_x_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (876-882)

    #pragma omp simd aligned(g_z_0_z_0_0_x_xz_xx, g_z_0_z_0_0_x_xz_xy, g_z_0_z_0_0_x_xz_xz, g_z_0_z_0_0_x_xz_yy, g_z_0_z_0_0_x_xz_yz, g_z_0_z_0_0_x_xz_zz, g_z_x_x_xx, g_z_x_x_xy, g_z_x_x_xz, g_z_x_x_yy, g_z_x_x_yz, g_z_x_x_zz, g_z_x_xzz_xx, g_z_x_xzz_xy, g_z_x_xzz_xz, g_z_x_xzz_yy, g_z_x_xzz_yz, g_z_x_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_x_xz_xx[i] = -2.0 * g_z_x_x_xx[i] * a_exp + 4.0 * g_z_x_xzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_xz_xy[i] = -2.0 * g_z_x_x_xy[i] * a_exp + 4.0 * g_z_x_xzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_xz_xz[i] = -2.0 * g_z_x_x_xz[i] * a_exp + 4.0 * g_z_x_xzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_xz_yy[i] = -2.0 * g_z_x_x_yy[i] * a_exp + 4.0 * g_z_x_xzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_xz_yz[i] = -2.0 * g_z_x_x_yz[i] * a_exp + 4.0 * g_z_x_xzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_xz_zz[i] = -2.0 * g_z_x_x_zz[i] * a_exp + 4.0 * g_z_x_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (882-888)

    #pragma omp simd aligned(g_z_0_z_0_0_x_yy_xx, g_z_0_z_0_0_x_yy_xy, g_z_0_z_0_0_x_yy_xz, g_z_0_z_0_0_x_yy_yy, g_z_0_z_0_0_x_yy_yz, g_z_0_z_0_0_x_yy_zz, g_z_x_yyz_xx, g_z_x_yyz_xy, g_z_x_yyz_xz, g_z_x_yyz_yy, g_z_x_yyz_yz, g_z_x_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_x_yy_xx[i] = 4.0 * g_z_x_yyz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_yy_xy[i] = 4.0 * g_z_x_yyz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_yy_xz[i] = 4.0 * g_z_x_yyz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_yy_yy[i] = 4.0 * g_z_x_yyz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_yy_yz[i] = 4.0 * g_z_x_yyz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_yy_zz[i] = 4.0 * g_z_x_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (888-894)

    #pragma omp simd aligned(g_z_0_z_0_0_x_yz_xx, g_z_0_z_0_0_x_yz_xy, g_z_0_z_0_0_x_yz_xz, g_z_0_z_0_0_x_yz_yy, g_z_0_z_0_0_x_yz_yz, g_z_0_z_0_0_x_yz_zz, g_z_x_y_xx, g_z_x_y_xy, g_z_x_y_xz, g_z_x_y_yy, g_z_x_y_yz, g_z_x_y_zz, g_z_x_yzz_xx, g_z_x_yzz_xy, g_z_x_yzz_xz, g_z_x_yzz_yy, g_z_x_yzz_yz, g_z_x_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_x_yz_xx[i] = -2.0 * g_z_x_y_xx[i] * a_exp + 4.0 * g_z_x_yzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_yz_xy[i] = -2.0 * g_z_x_y_xy[i] * a_exp + 4.0 * g_z_x_yzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_yz_xz[i] = -2.0 * g_z_x_y_xz[i] * a_exp + 4.0 * g_z_x_yzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_yz_yy[i] = -2.0 * g_z_x_y_yy[i] * a_exp + 4.0 * g_z_x_yzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_yz_yz[i] = -2.0 * g_z_x_y_yz[i] * a_exp + 4.0 * g_z_x_yzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_yz_zz[i] = -2.0 * g_z_x_y_zz[i] * a_exp + 4.0 * g_z_x_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (894-900)

    #pragma omp simd aligned(g_z_0_z_0_0_x_zz_xx, g_z_0_z_0_0_x_zz_xy, g_z_0_z_0_0_x_zz_xz, g_z_0_z_0_0_x_zz_yy, g_z_0_z_0_0_x_zz_yz, g_z_0_z_0_0_x_zz_zz, g_z_x_z_xx, g_z_x_z_xy, g_z_x_z_xz, g_z_x_z_yy, g_z_x_z_yz, g_z_x_z_zz, g_z_x_zzz_xx, g_z_x_zzz_xy, g_z_x_zzz_xz, g_z_x_zzz_yy, g_z_x_zzz_yz, g_z_x_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_x_zz_xx[i] = -4.0 * g_z_x_z_xx[i] * a_exp + 4.0 * g_z_x_zzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_zz_xy[i] = -4.0 * g_z_x_z_xy[i] * a_exp + 4.0 * g_z_x_zzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_zz_xz[i] = -4.0 * g_z_x_z_xz[i] * a_exp + 4.0 * g_z_x_zzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_zz_yy[i] = -4.0 * g_z_x_z_yy[i] * a_exp + 4.0 * g_z_x_zzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_zz_yz[i] = -4.0 * g_z_x_z_yz[i] * a_exp + 4.0 * g_z_x_zzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_zz_zz[i] = -4.0 * g_z_x_z_zz[i] * a_exp + 4.0 * g_z_x_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (900-906)

    #pragma omp simd aligned(g_z_0_z_0_0_y_xx_xx, g_z_0_z_0_0_y_xx_xy, g_z_0_z_0_0_y_xx_xz, g_z_0_z_0_0_y_xx_yy, g_z_0_z_0_0_y_xx_yz, g_z_0_z_0_0_y_xx_zz, g_z_y_xxz_xx, g_z_y_xxz_xy, g_z_y_xxz_xz, g_z_y_xxz_yy, g_z_y_xxz_yz, g_z_y_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_y_xx_xx[i] = 4.0 * g_z_y_xxz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_xx_xy[i] = 4.0 * g_z_y_xxz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_xx_xz[i] = 4.0 * g_z_y_xxz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_xx_yy[i] = 4.0 * g_z_y_xxz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_xx_yz[i] = 4.0 * g_z_y_xxz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_xx_zz[i] = 4.0 * g_z_y_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (906-912)

    #pragma omp simd aligned(g_z_0_z_0_0_y_xy_xx, g_z_0_z_0_0_y_xy_xy, g_z_0_z_0_0_y_xy_xz, g_z_0_z_0_0_y_xy_yy, g_z_0_z_0_0_y_xy_yz, g_z_0_z_0_0_y_xy_zz, g_z_y_xyz_xx, g_z_y_xyz_xy, g_z_y_xyz_xz, g_z_y_xyz_yy, g_z_y_xyz_yz, g_z_y_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_y_xy_xx[i] = 4.0 * g_z_y_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_xy_xy[i] = 4.0 * g_z_y_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_xy_xz[i] = 4.0 * g_z_y_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_xy_yy[i] = 4.0 * g_z_y_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_xy_yz[i] = 4.0 * g_z_y_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_xy_zz[i] = 4.0 * g_z_y_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (912-918)

    #pragma omp simd aligned(g_z_0_z_0_0_y_xz_xx, g_z_0_z_0_0_y_xz_xy, g_z_0_z_0_0_y_xz_xz, g_z_0_z_0_0_y_xz_yy, g_z_0_z_0_0_y_xz_yz, g_z_0_z_0_0_y_xz_zz, g_z_y_x_xx, g_z_y_x_xy, g_z_y_x_xz, g_z_y_x_yy, g_z_y_x_yz, g_z_y_x_zz, g_z_y_xzz_xx, g_z_y_xzz_xy, g_z_y_xzz_xz, g_z_y_xzz_yy, g_z_y_xzz_yz, g_z_y_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_y_xz_xx[i] = -2.0 * g_z_y_x_xx[i] * a_exp + 4.0 * g_z_y_xzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_xz_xy[i] = -2.0 * g_z_y_x_xy[i] * a_exp + 4.0 * g_z_y_xzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_xz_xz[i] = -2.0 * g_z_y_x_xz[i] * a_exp + 4.0 * g_z_y_xzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_xz_yy[i] = -2.0 * g_z_y_x_yy[i] * a_exp + 4.0 * g_z_y_xzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_xz_yz[i] = -2.0 * g_z_y_x_yz[i] * a_exp + 4.0 * g_z_y_xzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_xz_zz[i] = -2.0 * g_z_y_x_zz[i] * a_exp + 4.0 * g_z_y_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (918-924)

    #pragma omp simd aligned(g_z_0_z_0_0_y_yy_xx, g_z_0_z_0_0_y_yy_xy, g_z_0_z_0_0_y_yy_xz, g_z_0_z_0_0_y_yy_yy, g_z_0_z_0_0_y_yy_yz, g_z_0_z_0_0_y_yy_zz, g_z_y_yyz_xx, g_z_y_yyz_xy, g_z_y_yyz_xz, g_z_y_yyz_yy, g_z_y_yyz_yz, g_z_y_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_y_yy_xx[i] = 4.0 * g_z_y_yyz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_yy_xy[i] = 4.0 * g_z_y_yyz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_yy_xz[i] = 4.0 * g_z_y_yyz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_yy_yy[i] = 4.0 * g_z_y_yyz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_yy_yz[i] = 4.0 * g_z_y_yyz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_yy_zz[i] = 4.0 * g_z_y_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (924-930)

    #pragma omp simd aligned(g_z_0_z_0_0_y_yz_xx, g_z_0_z_0_0_y_yz_xy, g_z_0_z_0_0_y_yz_xz, g_z_0_z_0_0_y_yz_yy, g_z_0_z_0_0_y_yz_yz, g_z_0_z_0_0_y_yz_zz, g_z_y_y_xx, g_z_y_y_xy, g_z_y_y_xz, g_z_y_y_yy, g_z_y_y_yz, g_z_y_y_zz, g_z_y_yzz_xx, g_z_y_yzz_xy, g_z_y_yzz_xz, g_z_y_yzz_yy, g_z_y_yzz_yz, g_z_y_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_y_yz_xx[i] = -2.0 * g_z_y_y_xx[i] * a_exp + 4.0 * g_z_y_yzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_yz_xy[i] = -2.0 * g_z_y_y_xy[i] * a_exp + 4.0 * g_z_y_yzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_yz_xz[i] = -2.0 * g_z_y_y_xz[i] * a_exp + 4.0 * g_z_y_yzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_yz_yy[i] = -2.0 * g_z_y_y_yy[i] * a_exp + 4.0 * g_z_y_yzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_yz_yz[i] = -2.0 * g_z_y_y_yz[i] * a_exp + 4.0 * g_z_y_yzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_yz_zz[i] = -2.0 * g_z_y_y_zz[i] * a_exp + 4.0 * g_z_y_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (930-936)

    #pragma omp simd aligned(g_z_0_z_0_0_y_zz_xx, g_z_0_z_0_0_y_zz_xy, g_z_0_z_0_0_y_zz_xz, g_z_0_z_0_0_y_zz_yy, g_z_0_z_0_0_y_zz_yz, g_z_0_z_0_0_y_zz_zz, g_z_y_z_xx, g_z_y_z_xy, g_z_y_z_xz, g_z_y_z_yy, g_z_y_z_yz, g_z_y_z_zz, g_z_y_zzz_xx, g_z_y_zzz_xy, g_z_y_zzz_xz, g_z_y_zzz_yy, g_z_y_zzz_yz, g_z_y_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_y_zz_xx[i] = -4.0 * g_z_y_z_xx[i] * a_exp + 4.0 * g_z_y_zzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_zz_xy[i] = -4.0 * g_z_y_z_xy[i] * a_exp + 4.0 * g_z_y_zzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_zz_xz[i] = -4.0 * g_z_y_z_xz[i] * a_exp + 4.0 * g_z_y_zzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_zz_yy[i] = -4.0 * g_z_y_z_yy[i] * a_exp + 4.0 * g_z_y_zzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_zz_yz[i] = -4.0 * g_z_y_z_yz[i] * a_exp + 4.0 * g_z_y_zzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_zz_zz[i] = -4.0 * g_z_y_z_zz[i] * a_exp + 4.0 * g_z_y_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (936-942)

    #pragma omp simd aligned(g_z_0_z_0_0_z_xx_xx, g_z_0_z_0_0_z_xx_xy, g_z_0_z_0_0_z_xx_xz, g_z_0_z_0_0_z_xx_yy, g_z_0_z_0_0_z_xx_yz, g_z_0_z_0_0_z_xx_zz, g_z_z_xxz_xx, g_z_z_xxz_xy, g_z_z_xxz_xz, g_z_z_xxz_yy, g_z_z_xxz_yz, g_z_z_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_z_xx_xx[i] = 4.0 * g_z_z_xxz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_xx_xy[i] = 4.0 * g_z_z_xxz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_xx_xz[i] = 4.0 * g_z_z_xxz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_xx_yy[i] = 4.0 * g_z_z_xxz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_xx_yz[i] = 4.0 * g_z_z_xxz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_xx_zz[i] = 4.0 * g_z_z_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (942-948)

    #pragma omp simd aligned(g_z_0_z_0_0_z_xy_xx, g_z_0_z_0_0_z_xy_xy, g_z_0_z_0_0_z_xy_xz, g_z_0_z_0_0_z_xy_yy, g_z_0_z_0_0_z_xy_yz, g_z_0_z_0_0_z_xy_zz, g_z_z_xyz_xx, g_z_z_xyz_xy, g_z_z_xyz_xz, g_z_z_xyz_yy, g_z_z_xyz_yz, g_z_z_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_z_xy_xx[i] = 4.0 * g_z_z_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_xy_xy[i] = 4.0 * g_z_z_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_xy_xz[i] = 4.0 * g_z_z_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_xy_yy[i] = 4.0 * g_z_z_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_xy_yz[i] = 4.0 * g_z_z_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_xy_zz[i] = 4.0 * g_z_z_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (948-954)

    #pragma omp simd aligned(g_z_0_z_0_0_z_xz_xx, g_z_0_z_0_0_z_xz_xy, g_z_0_z_0_0_z_xz_xz, g_z_0_z_0_0_z_xz_yy, g_z_0_z_0_0_z_xz_yz, g_z_0_z_0_0_z_xz_zz, g_z_z_x_xx, g_z_z_x_xy, g_z_z_x_xz, g_z_z_x_yy, g_z_z_x_yz, g_z_z_x_zz, g_z_z_xzz_xx, g_z_z_xzz_xy, g_z_z_xzz_xz, g_z_z_xzz_yy, g_z_z_xzz_yz, g_z_z_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_z_xz_xx[i] = -2.0 * g_z_z_x_xx[i] * a_exp + 4.0 * g_z_z_xzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_xz_xy[i] = -2.0 * g_z_z_x_xy[i] * a_exp + 4.0 * g_z_z_xzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_xz_xz[i] = -2.0 * g_z_z_x_xz[i] * a_exp + 4.0 * g_z_z_xzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_xz_yy[i] = -2.0 * g_z_z_x_yy[i] * a_exp + 4.0 * g_z_z_xzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_xz_yz[i] = -2.0 * g_z_z_x_yz[i] * a_exp + 4.0 * g_z_z_xzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_xz_zz[i] = -2.0 * g_z_z_x_zz[i] * a_exp + 4.0 * g_z_z_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (954-960)

    #pragma omp simd aligned(g_z_0_z_0_0_z_yy_xx, g_z_0_z_0_0_z_yy_xy, g_z_0_z_0_0_z_yy_xz, g_z_0_z_0_0_z_yy_yy, g_z_0_z_0_0_z_yy_yz, g_z_0_z_0_0_z_yy_zz, g_z_z_yyz_xx, g_z_z_yyz_xy, g_z_z_yyz_xz, g_z_z_yyz_yy, g_z_z_yyz_yz, g_z_z_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_z_yy_xx[i] = 4.0 * g_z_z_yyz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_yy_xy[i] = 4.0 * g_z_z_yyz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_yy_xz[i] = 4.0 * g_z_z_yyz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_yy_yy[i] = 4.0 * g_z_z_yyz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_yy_yz[i] = 4.0 * g_z_z_yyz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_yy_zz[i] = 4.0 * g_z_z_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (960-966)

    #pragma omp simd aligned(g_z_0_z_0_0_z_yz_xx, g_z_0_z_0_0_z_yz_xy, g_z_0_z_0_0_z_yz_xz, g_z_0_z_0_0_z_yz_yy, g_z_0_z_0_0_z_yz_yz, g_z_0_z_0_0_z_yz_zz, g_z_z_y_xx, g_z_z_y_xy, g_z_z_y_xz, g_z_z_y_yy, g_z_z_y_yz, g_z_z_y_zz, g_z_z_yzz_xx, g_z_z_yzz_xy, g_z_z_yzz_xz, g_z_z_yzz_yy, g_z_z_yzz_yz, g_z_z_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_z_yz_xx[i] = -2.0 * g_z_z_y_xx[i] * a_exp + 4.0 * g_z_z_yzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_yz_xy[i] = -2.0 * g_z_z_y_xy[i] * a_exp + 4.0 * g_z_z_yzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_yz_xz[i] = -2.0 * g_z_z_y_xz[i] * a_exp + 4.0 * g_z_z_yzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_yz_yy[i] = -2.0 * g_z_z_y_yy[i] * a_exp + 4.0 * g_z_z_yzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_yz_yz[i] = -2.0 * g_z_z_y_yz[i] * a_exp + 4.0 * g_z_z_yzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_yz_zz[i] = -2.0 * g_z_z_y_zz[i] * a_exp + 4.0 * g_z_z_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (966-972)

    #pragma omp simd aligned(g_z_0_z_0_0_z_zz_xx, g_z_0_z_0_0_z_zz_xy, g_z_0_z_0_0_z_zz_xz, g_z_0_z_0_0_z_zz_yy, g_z_0_z_0_0_z_zz_yz, g_z_0_z_0_0_z_zz_zz, g_z_z_z_xx, g_z_z_z_xy, g_z_z_z_xz, g_z_z_z_yy, g_z_z_z_yz, g_z_z_z_zz, g_z_z_zzz_xx, g_z_z_zzz_xy, g_z_z_zzz_xz, g_z_z_zzz_yy, g_z_z_zzz_yz, g_z_z_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_z_zz_xx[i] = -4.0 * g_z_z_z_xx[i] * a_exp + 4.0 * g_z_z_zzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_zz_xy[i] = -4.0 * g_z_z_z_xy[i] * a_exp + 4.0 * g_z_z_zzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_zz_xz[i] = -4.0 * g_z_z_z_xz[i] * a_exp + 4.0 * g_z_z_zzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_zz_yy[i] = -4.0 * g_z_z_z_yy[i] * a_exp + 4.0 * g_z_z_zzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_zz_yz[i] = -4.0 * g_z_z_z_yz[i] * a_exp + 4.0 * g_z_z_zzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_zz_zz[i] = -4.0 * g_z_z_z_zz[i] * a_exp + 4.0 * g_z_z_zzz_zz[i] * a_exp * c_exps[i];
    }
}

} // t4c_geom namespace


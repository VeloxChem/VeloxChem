#include "GeomDeriv1010OfScalarForSDDD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1010_sddd_0(CSimdArray<double>& buffer_1010_sddd,
                     const CSimdArray<double>& buffer_pdpd,
                     const CSimdArray<double>& buffer_pdfd,
                     const double a_exp,
                     const double* c_exps) -> void
{
    const auto ndims = buffer_1010_sddd.number_of_columns();

    /// Set up components of auxilary buffer : buffer_pdpd

    auto g_x_xx_x_xx = buffer_pdpd[0];

    auto g_x_xx_x_xy = buffer_pdpd[1];

    auto g_x_xx_x_xz = buffer_pdpd[2];

    auto g_x_xx_x_yy = buffer_pdpd[3];

    auto g_x_xx_x_yz = buffer_pdpd[4];

    auto g_x_xx_x_zz = buffer_pdpd[5];

    auto g_x_xx_y_xx = buffer_pdpd[6];

    auto g_x_xx_y_xy = buffer_pdpd[7];

    auto g_x_xx_y_xz = buffer_pdpd[8];

    auto g_x_xx_y_yy = buffer_pdpd[9];

    auto g_x_xx_y_yz = buffer_pdpd[10];

    auto g_x_xx_y_zz = buffer_pdpd[11];

    auto g_x_xx_z_xx = buffer_pdpd[12];

    auto g_x_xx_z_xy = buffer_pdpd[13];

    auto g_x_xx_z_xz = buffer_pdpd[14];

    auto g_x_xx_z_yy = buffer_pdpd[15];

    auto g_x_xx_z_yz = buffer_pdpd[16];

    auto g_x_xx_z_zz = buffer_pdpd[17];

    auto g_x_xy_x_xx = buffer_pdpd[18];

    auto g_x_xy_x_xy = buffer_pdpd[19];

    auto g_x_xy_x_xz = buffer_pdpd[20];

    auto g_x_xy_x_yy = buffer_pdpd[21];

    auto g_x_xy_x_yz = buffer_pdpd[22];

    auto g_x_xy_x_zz = buffer_pdpd[23];

    auto g_x_xy_y_xx = buffer_pdpd[24];

    auto g_x_xy_y_xy = buffer_pdpd[25];

    auto g_x_xy_y_xz = buffer_pdpd[26];

    auto g_x_xy_y_yy = buffer_pdpd[27];

    auto g_x_xy_y_yz = buffer_pdpd[28];

    auto g_x_xy_y_zz = buffer_pdpd[29];

    auto g_x_xy_z_xx = buffer_pdpd[30];

    auto g_x_xy_z_xy = buffer_pdpd[31];

    auto g_x_xy_z_xz = buffer_pdpd[32];

    auto g_x_xy_z_yy = buffer_pdpd[33];

    auto g_x_xy_z_yz = buffer_pdpd[34];

    auto g_x_xy_z_zz = buffer_pdpd[35];

    auto g_x_xz_x_xx = buffer_pdpd[36];

    auto g_x_xz_x_xy = buffer_pdpd[37];

    auto g_x_xz_x_xz = buffer_pdpd[38];

    auto g_x_xz_x_yy = buffer_pdpd[39];

    auto g_x_xz_x_yz = buffer_pdpd[40];

    auto g_x_xz_x_zz = buffer_pdpd[41];

    auto g_x_xz_y_xx = buffer_pdpd[42];

    auto g_x_xz_y_xy = buffer_pdpd[43];

    auto g_x_xz_y_xz = buffer_pdpd[44];

    auto g_x_xz_y_yy = buffer_pdpd[45];

    auto g_x_xz_y_yz = buffer_pdpd[46];

    auto g_x_xz_y_zz = buffer_pdpd[47];

    auto g_x_xz_z_xx = buffer_pdpd[48];

    auto g_x_xz_z_xy = buffer_pdpd[49];

    auto g_x_xz_z_xz = buffer_pdpd[50];

    auto g_x_xz_z_yy = buffer_pdpd[51];

    auto g_x_xz_z_yz = buffer_pdpd[52];

    auto g_x_xz_z_zz = buffer_pdpd[53];

    auto g_x_yy_x_xx = buffer_pdpd[54];

    auto g_x_yy_x_xy = buffer_pdpd[55];

    auto g_x_yy_x_xz = buffer_pdpd[56];

    auto g_x_yy_x_yy = buffer_pdpd[57];

    auto g_x_yy_x_yz = buffer_pdpd[58];

    auto g_x_yy_x_zz = buffer_pdpd[59];

    auto g_x_yy_y_xx = buffer_pdpd[60];

    auto g_x_yy_y_xy = buffer_pdpd[61];

    auto g_x_yy_y_xz = buffer_pdpd[62];

    auto g_x_yy_y_yy = buffer_pdpd[63];

    auto g_x_yy_y_yz = buffer_pdpd[64];

    auto g_x_yy_y_zz = buffer_pdpd[65];

    auto g_x_yy_z_xx = buffer_pdpd[66];

    auto g_x_yy_z_xy = buffer_pdpd[67];

    auto g_x_yy_z_xz = buffer_pdpd[68];

    auto g_x_yy_z_yy = buffer_pdpd[69];

    auto g_x_yy_z_yz = buffer_pdpd[70];

    auto g_x_yy_z_zz = buffer_pdpd[71];

    auto g_x_yz_x_xx = buffer_pdpd[72];

    auto g_x_yz_x_xy = buffer_pdpd[73];

    auto g_x_yz_x_xz = buffer_pdpd[74];

    auto g_x_yz_x_yy = buffer_pdpd[75];

    auto g_x_yz_x_yz = buffer_pdpd[76];

    auto g_x_yz_x_zz = buffer_pdpd[77];

    auto g_x_yz_y_xx = buffer_pdpd[78];

    auto g_x_yz_y_xy = buffer_pdpd[79];

    auto g_x_yz_y_xz = buffer_pdpd[80];

    auto g_x_yz_y_yy = buffer_pdpd[81];

    auto g_x_yz_y_yz = buffer_pdpd[82];

    auto g_x_yz_y_zz = buffer_pdpd[83];

    auto g_x_yz_z_xx = buffer_pdpd[84];

    auto g_x_yz_z_xy = buffer_pdpd[85];

    auto g_x_yz_z_xz = buffer_pdpd[86];

    auto g_x_yz_z_yy = buffer_pdpd[87];

    auto g_x_yz_z_yz = buffer_pdpd[88];

    auto g_x_yz_z_zz = buffer_pdpd[89];

    auto g_x_zz_x_xx = buffer_pdpd[90];

    auto g_x_zz_x_xy = buffer_pdpd[91];

    auto g_x_zz_x_xz = buffer_pdpd[92];

    auto g_x_zz_x_yy = buffer_pdpd[93];

    auto g_x_zz_x_yz = buffer_pdpd[94];

    auto g_x_zz_x_zz = buffer_pdpd[95];

    auto g_x_zz_y_xx = buffer_pdpd[96];

    auto g_x_zz_y_xy = buffer_pdpd[97];

    auto g_x_zz_y_xz = buffer_pdpd[98];

    auto g_x_zz_y_yy = buffer_pdpd[99];

    auto g_x_zz_y_yz = buffer_pdpd[100];

    auto g_x_zz_y_zz = buffer_pdpd[101];

    auto g_x_zz_z_xx = buffer_pdpd[102];

    auto g_x_zz_z_xy = buffer_pdpd[103];

    auto g_x_zz_z_xz = buffer_pdpd[104];

    auto g_x_zz_z_yy = buffer_pdpd[105];

    auto g_x_zz_z_yz = buffer_pdpd[106];

    auto g_x_zz_z_zz = buffer_pdpd[107];

    auto g_y_xx_x_xx = buffer_pdpd[108];

    auto g_y_xx_x_xy = buffer_pdpd[109];

    auto g_y_xx_x_xz = buffer_pdpd[110];

    auto g_y_xx_x_yy = buffer_pdpd[111];

    auto g_y_xx_x_yz = buffer_pdpd[112];

    auto g_y_xx_x_zz = buffer_pdpd[113];

    auto g_y_xx_y_xx = buffer_pdpd[114];

    auto g_y_xx_y_xy = buffer_pdpd[115];

    auto g_y_xx_y_xz = buffer_pdpd[116];

    auto g_y_xx_y_yy = buffer_pdpd[117];

    auto g_y_xx_y_yz = buffer_pdpd[118];

    auto g_y_xx_y_zz = buffer_pdpd[119];

    auto g_y_xx_z_xx = buffer_pdpd[120];

    auto g_y_xx_z_xy = buffer_pdpd[121];

    auto g_y_xx_z_xz = buffer_pdpd[122];

    auto g_y_xx_z_yy = buffer_pdpd[123];

    auto g_y_xx_z_yz = buffer_pdpd[124];

    auto g_y_xx_z_zz = buffer_pdpd[125];

    auto g_y_xy_x_xx = buffer_pdpd[126];

    auto g_y_xy_x_xy = buffer_pdpd[127];

    auto g_y_xy_x_xz = buffer_pdpd[128];

    auto g_y_xy_x_yy = buffer_pdpd[129];

    auto g_y_xy_x_yz = buffer_pdpd[130];

    auto g_y_xy_x_zz = buffer_pdpd[131];

    auto g_y_xy_y_xx = buffer_pdpd[132];

    auto g_y_xy_y_xy = buffer_pdpd[133];

    auto g_y_xy_y_xz = buffer_pdpd[134];

    auto g_y_xy_y_yy = buffer_pdpd[135];

    auto g_y_xy_y_yz = buffer_pdpd[136];

    auto g_y_xy_y_zz = buffer_pdpd[137];

    auto g_y_xy_z_xx = buffer_pdpd[138];

    auto g_y_xy_z_xy = buffer_pdpd[139];

    auto g_y_xy_z_xz = buffer_pdpd[140];

    auto g_y_xy_z_yy = buffer_pdpd[141];

    auto g_y_xy_z_yz = buffer_pdpd[142];

    auto g_y_xy_z_zz = buffer_pdpd[143];

    auto g_y_xz_x_xx = buffer_pdpd[144];

    auto g_y_xz_x_xy = buffer_pdpd[145];

    auto g_y_xz_x_xz = buffer_pdpd[146];

    auto g_y_xz_x_yy = buffer_pdpd[147];

    auto g_y_xz_x_yz = buffer_pdpd[148];

    auto g_y_xz_x_zz = buffer_pdpd[149];

    auto g_y_xz_y_xx = buffer_pdpd[150];

    auto g_y_xz_y_xy = buffer_pdpd[151];

    auto g_y_xz_y_xz = buffer_pdpd[152];

    auto g_y_xz_y_yy = buffer_pdpd[153];

    auto g_y_xz_y_yz = buffer_pdpd[154];

    auto g_y_xz_y_zz = buffer_pdpd[155];

    auto g_y_xz_z_xx = buffer_pdpd[156];

    auto g_y_xz_z_xy = buffer_pdpd[157];

    auto g_y_xz_z_xz = buffer_pdpd[158];

    auto g_y_xz_z_yy = buffer_pdpd[159];

    auto g_y_xz_z_yz = buffer_pdpd[160];

    auto g_y_xz_z_zz = buffer_pdpd[161];

    auto g_y_yy_x_xx = buffer_pdpd[162];

    auto g_y_yy_x_xy = buffer_pdpd[163];

    auto g_y_yy_x_xz = buffer_pdpd[164];

    auto g_y_yy_x_yy = buffer_pdpd[165];

    auto g_y_yy_x_yz = buffer_pdpd[166];

    auto g_y_yy_x_zz = buffer_pdpd[167];

    auto g_y_yy_y_xx = buffer_pdpd[168];

    auto g_y_yy_y_xy = buffer_pdpd[169];

    auto g_y_yy_y_xz = buffer_pdpd[170];

    auto g_y_yy_y_yy = buffer_pdpd[171];

    auto g_y_yy_y_yz = buffer_pdpd[172];

    auto g_y_yy_y_zz = buffer_pdpd[173];

    auto g_y_yy_z_xx = buffer_pdpd[174];

    auto g_y_yy_z_xy = buffer_pdpd[175];

    auto g_y_yy_z_xz = buffer_pdpd[176];

    auto g_y_yy_z_yy = buffer_pdpd[177];

    auto g_y_yy_z_yz = buffer_pdpd[178];

    auto g_y_yy_z_zz = buffer_pdpd[179];

    auto g_y_yz_x_xx = buffer_pdpd[180];

    auto g_y_yz_x_xy = buffer_pdpd[181];

    auto g_y_yz_x_xz = buffer_pdpd[182];

    auto g_y_yz_x_yy = buffer_pdpd[183];

    auto g_y_yz_x_yz = buffer_pdpd[184];

    auto g_y_yz_x_zz = buffer_pdpd[185];

    auto g_y_yz_y_xx = buffer_pdpd[186];

    auto g_y_yz_y_xy = buffer_pdpd[187];

    auto g_y_yz_y_xz = buffer_pdpd[188];

    auto g_y_yz_y_yy = buffer_pdpd[189];

    auto g_y_yz_y_yz = buffer_pdpd[190];

    auto g_y_yz_y_zz = buffer_pdpd[191];

    auto g_y_yz_z_xx = buffer_pdpd[192];

    auto g_y_yz_z_xy = buffer_pdpd[193];

    auto g_y_yz_z_xz = buffer_pdpd[194];

    auto g_y_yz_z_yy = buffer_pdpd[195];

    auto g_y_yz_z_yz = buffer_pdpd[196];

    auto g_y_yz_z_zz = buffer_pdpd[197];

    auto g_y_zz_x_xx = buffer_pdpd[198];

    auto g_y_zz_x_xy = buffer_pdpd[199];

    auto g_y_zz_x_xz = buffer_pdpd[200];

    auto g_y_zz_x_yy = buffer_pdpd[201];

    auto g_y_zz_x_yz = buffer_pdpd[202];

    auto g_y_zz_x_zz = buffer_pdpd[203];

    auto g_y_zz_y_xx = buffer_pdpd[204];

    auto g_y_zz_y_xy = buffer_pdpd[205];

    auto g_y_zz_y_xz = buffer_pdpd[206];

    auto g_y_zz_y_yy = buffer_pdpd[207];

    auto g_y_zz_y_yz = buffer_pdpd[208];

    auto g_y_zz_y_zz = buffer_pdpd[209];

    auto g_y_zz_z_xx = buffer_pdpd[210];

    auto g_y_zz_z_xy = buffer_pdpd[211];

    auto g_y_zz_z_xz = buffer_pdpd[212];

    auto g_y_zz_z_yy = buffer_pdpd[213];

    auto g_y_zz_z_yz = buffer_pdpd[214];

    auto g_y_zz_z_zz = buffer_pdpd[215];

    auto g_z_xx_x_xx = buffer_pdpd[216];

    auto g_z_xx_x_xy = buffer_pdpd[217];

    auto g_z_xx_x_xz = buffer_pdpd[218];

    auto g_z_xx_x_yy = buffer_pdpd[219];

    auto g_z_xx_x_yz = buffer_pdpd[220];

    auto g_z_xx_x_zz = buffer_pdpd[221];

    auto g_z_xx_y_xx = buffer_pdpd[222];

    auto g_z_xx_y_xy = buffer_pdpd[223];

    auto g_z_xx_y_xz = buffer_pdpd[224];

    auto g_z_xx_y_yy = buffer_pdpd[225];

    auto g_z_xx_y_yz = buffer_pdpd[226];

    auto g_z_xx_y_zz = buffer_pdpd[227];

    auto g_z_xx_z_xx = buffer_pdpd[228];

    auto g_z_xx_z_xy = buffer_pdpd[229];

    auto g_z_xx_z_xz = buffer_pdpd[230];

    auto g_z_xx_z_yy = buffer_pdpd[231];

    auto g_z_xx_z_yz = buffer_pdpd[232];

    auto g_z_xx_z_zz = buffer_pdpd[233];

    auto g_z_xy_x_xx = buffer_pdpd[234];

    auto g_z_xy_x_xy = buffer_pdpd[235];

    auto g_z_xy_x_xz = buffer_pdpd[236];

    auto g_z_xy_x_yy = buffer_pdpd[237];

    auto g_z_xy_x_yz = buffer_pdpd[238];

    auto g_z_xy_x_zz = buffer_pdpd[239];

    auto g_z_xy_y_xx = buffer_pdpd[240];

    auto g_z_xy_y_xy = buffer_pdpd[241];

    auto g_z_xy_y_xz = buffer_pdpd[242];

    auto g_z_xy_y_yy = buffer_pdpd[243];

    auto g_z_xy_y_yz = buffer_pdpd[244];

    auto g_z_xy_y_zz = buffer_pdpd[245];

    auto g_z_xy_z_xx = buffer_pdpd[246];

    auto g_z_xy_z_xy = buffer_pdpd[247];

    auto g_z_xy_z_xz = buffer_pdpd[248];

    auto g_z_xy_z_yy = buffer_pdpd[249];

    auto g_z_xy_z_yz = buffer_pdpd[250];

    auto g_z_xy_z_zz = buffer_pdpd[251];

    auto g_z_xz_x_xx = buffer_pdpd[252];

    auto g_z_xz_x_xy = buffer_pdpd[253];

    auto g_z_xz_x_xz = buffer_pdpd[254];

    auto g_z_xz_x_yy = buffer_pdpd[255];

    auto g_z_xz_x_yz = buffer_pdpd[256];

    auto g_z_xz_x_zz = buffer_pdpd[257];

    auto g_z_xz_y_xx = buffer_pdpd[258];

    auto g_z_xz_y_xy = buffer_pdpd[259];

    auto g_z_xz_y_xz = buffer_pdpd[260];

    auto g_z_xz_y_yy = buffer_pdpd[261];

    auto g_z_xz_y_yz = buffer_pdpd[262];

    auto g_z_xz_y_zz = buffer_pdpd[263];

    auto g_z_xz_z_xx = buffer_pdpd[264];

    auto g_z_xz_z_xy = buffer_pdpd[265];

    auto g_z_xz_z_xz = buffer_pdpd[266];

    auto g_z_xz_z_yy = buffer_pdpd[267];

    auto g_z_xz_z_yz = buffer_pdpd[268];

    auto g_z_xz_z_zz = buffer_pdpd[269];

    auto g_z_yy_x_xx = buffer_pdpd[270];

    auto g_z_yy_x_xy = buffer_pdpd[271];

    auto g_z_yy_x_xz = buffer_pdpd[272];

    auto g_z_yy_x_yy = buffer_pdpd[273];

    auto g_z_yy_x_yz = buffer_pdpd[274];

    auto g_z_yy_x_zz = buffer_pdpd[275];

    auto g_z_yy_y_xx = buffer_pdpd[276];

    auto g_z_yy_y_xy = buffer_pdpd[277];

    auto g_z_yy_y_xz = buffer_pdpd[278];

    auto g_z_yy_y_yy = buffer_pdpd[279];

    auto g_z_yy_y_yz = buffer_pdpd[280];

    auto g_z_yy_y_zz = buffer_pdpd[281];

    auto g_z_yy_z_xx = buffer_pdpd[282];

    auto g_z_yy_z_xy = buffer_pdpd[283];

    auto g_z_yy_z_xz = buffer_pdpd[284];

    auto g_z_yy_z_yy = buffer_pdpd[285];

    auto g_z_yy_z_yz = buffer_pdpd[286];

    auto g_z_yy_z_zz = buffer_pdpd[287];

    auto g_z_yz_x_xx = buffer_pdpd[288];

    auto g_z_yz_x_xy = buffer_pdpd[289];

    auto g_z_yz_x_xz = buffer_pdpd[290];

    auto g_z_yz_x_yy = buffer_pdpd[291];

    auto g_z_yz_x_yz = buffer_pdpd[292];

    auto g_z_yz_x_zz = buffer_pdpd[293];

    auto g_z_yz_y_xx = buffer_pdpd[294];

    auto g_z_yz_y_xy = buffer_pdpd[295];

    auto g_z_yz_y_xz = buffer_pdpd[296];

    auto g_z_yz_y_yy = buffer_pdpd[297];

    auto g_z_yz_y_yz = buffer_pdpd[298];

    auto g_z_yz_y_zz = buffer_pdpd[299];

    auto g_z_yz_z_xx = buffer_pdpd[300];

    auto g_z_yz_z_xy = buffer_pdpd[301];

    auto g_z_yz_z_xz = buffer_pdpd[302];

    auto g_z_yz_z_yy = buffer_pdpd[303];

    auto g_z_yz_z_yz = buffer_pdpd[304];

    auto g_z_yz_z_zz = buffer_pdpd[305];

    auto g_z_zz_x_xx = buffer_pdpd[306];

    auto g_z_zz_x_xy = buffer_pdpd[307];

    auto g_z_zz_x_xz = buffer_pdpd[308];

    auto g_z_zz_x_yy = buffer_pdpd[309];

    auto g_z_zz_x_yz = buffer_pdpd[310];

    auto g_z_zz_x_zz = buffer_pdpd[311];

    auto g_z_zz_y_xx = buffer_pdpd[312];

    auto g_z_zz_y_xy = buffer_pdpd[313];

    auto g_z_zz_y_xz = buffer_pdpd[314];

    auto g_z_zz_y_yy = buffer_pdpd[315];

    auto g_z_zz_y_yz = buffer_pdpd[316];

    auto g_z_zz_y_zz = buffer_pdpd[317];

    auto g_z_zz_z_xx = buffer_pdpd[318];

    auto g_z_zz_z_xy = buffer_pdpd[319];

    auto g_z_zz_z_xz = buffer_pdpd[320];

    auto g_z_zz_z_yy = buffer_pdpd[321];

    auto g_z_zz_z_yz = buffer_pdpd[322];

    auto g_z_zz_z_zz = buffer_pdpd[323];

    /// Set up components of auxilary buffer : buffer_pdfd

    auto g_x_xx_xxx_xx = buffer_pdfd[0];

    auto g_x_xx_xxx_xy = buffer_pdfd[1];

    auto g_x_xx_xxx_xz = buffer_pdfd[2];

    auto g_x_xx_xxx_yy = buffer_pdfd[3];

    auto g_x_xx_xxx_yz = buffer_pdfd[4];

    auto g_x_xx_xxx_zz = buffer_pdfd[5];

    auto g_x_xx_xxy_xx = buffer_pdfd[6];

    auto g_x_xx_xxy_xy = buffer_pdfd[7];

    auto g_x_xx_xxy_xz = buffer_pdfd[8];

    auto g_x_xx_xxy_yy = buffer_pdfd[9];

    auto g_x_xx_xxy_yz = buffer_pdfd[10];

    auto g_x_xx_xxy_zz = buffer_pdfd[11];

    auto g_x_xx_xxz_xx = buffer_pdfd[12];

    auto g_x_xx_xxz_xy = buffer_pdfd[13];

    auto g_x_xx_xxz_xz = buffer_pdfd[14];

    auto g_x_xx_xxz_yy = buffer_pdfd[15];

    auto g_x_xx_xxz_yz = buffer_pdfd[16];

    auto g_x_xx_xxz_zz = buffer_pdfd[17];

    auto g_x_xx_xyy_xx = buffer_pdfd[18];

    auto g_x_xx_xyy_xy = buffer_pdfd[19];

    auto g_x_xx_xyy_xz = buffer_pdfd[20];

    auto g_x_xx_xyy_yy = buffer_pdfd[21];

    auto g_x_xx_xyy_yz = buffer_pdfd[22];

    auto g_x_xx_xyy_zz = buffer_pdfd[23];

    auto g_x_xx_xyz_xx = buffer_pdfd[24];

    auto g_x_xx_xyz_xy = buffer_pdfd[25];

    auto g_x_xx_xyz_xz = buffer_pdfd[26];

    auto g_x_xx_xyz_yy = buffer_pdfd[27];

    auto g_x_xx_xyz_yz = buffer_pdfd[28];

    auto g_x_xx_xyz_zz = buffer_pdfd[29];

    auto g_x_xx_xzz_xx = buffer_pdfd[30];

    auto g_x_xx_xzz_xy = buffer_pdfd[31];

    auto g_x_xx_xzz_xz = buffer_pdfd[32];

    auto g_x_xx_xzz_yy = buffer_pdfd[33];

    auto g_x_xx_xzz_yz = buffer_pdfd[34];

    auto g_x_xx_xzz_zz = buffer_pdfd[35];

    auto g_x_xx_yyy_xx = buffer_pdfd[36];

    auto g_x_xx_yyy_xy = buffer_pdfd[37];

    auto g_x_xx_yyy_xz = buffer_pdfd[38];

    auto g_x_xx_yyy_yy = buffer_pdfd[39];

    auto g_x_xx_yyy_yz = buffer_pdfd[40];

    auto g_x_xx_yyy_zz = buffer_pdfd[41];

    auto g_x_xx_yyz_xx = buffer_pdfd[42];

    auto g_x_xx_yyz_xy = buffer_pdfd[43];

    auto g_x_xx_yyz_xz = buffer_pdfd[44];

    auto g_x_xx_yyz_yy = buffer_pdfd[45];

    auto g_x_xx_yyz_yz = buffer_pdfd[46];

    auto g_x_xx_yyz_zz = buffer_pdfd[47];

    auto g_x_xx_yzz_xx = buffer_pdfd[48];

    auto g_x_xx_yzz_xy = buffer_pdfd[49];

    auto g_x_xx_yzz_xz = buffer_pdfd[50];

    auto g_x_xx_yzz_yy = buffer_pdfd[51];

    auto g_x_xx_yzz_yz = buffer_pdfd[52];

    auto g_x_xx_yzz_zz = buffer_pdfd[53];

    auto g_x_xx_zzz_xx = buffer_pdfd[54];

    auto g_x_xx_zzz_xy = buffer_pdfd[55];

    auto g_x_xx_zzz_xz = buffer_pdfd[56];

    auto g_x_xx_zzz_yy = buffer_pdfd[57];

    auto g_x_xx_zzz_yz = buffer_pdfd[58];

    auto g_x_xx_zzz_zz = buffer_pdfd[59];

    auto g_x_xy_xxx_xx = buffer_pdfd[60];

    auto g_x_xy_xxx_xy = buffer_pdfd[61];

    auto g_x_xy_xxx_xz = buffer_pdfd[62];

    auto g_x_xy_xxx_yy = buffer_pdfd[63];

    auto g_x_xy_xxx_yz = buffer_pdfd[64];

    auto g_x_xy_xxx_zz = buffer_pdfd[65];

    auto g_x_xy_xxy_xx = buffer_pdfd[66];

    auto g_x_xy_xxy_xy = buffer_pdfd[67];

    auto g_x_xy_xxy_xz = buffer_pdfd[68];

    auto g_x_xy_xxy_yy = buffer_pdfd[69];

    auto g_x_xy_xxy_yz = buffer_pdfd[70];

    auto g_x_xy_xxy_zz = buffer_pdfd[71];

    auto g_x_xy_xxz_xx = buffer_pdfd[72];

    auto g_x_xy_xxz_xy = buffer_pdfd[73];

    auto g_x_xy_xxz_xz = buffer_pdfd[74];

    auto g_x_xy_xxz_yy = buffer_pdfd[75];

    auto g_x_xy_xxz_yz = buffer_pdfd[76];

    auto g_x_xy_xxz_zz = buffer_pdfd[77];

    auto g_x_xy_xyy_xx = buffer_pdfd[78];

    auto g_x_xy_xyy_xy = buffer_pdfd[79];

    auto g_x_xy_xyy_xz = buffer_pdfd[80];

    auto g_x_xy_xyy_yy = buffer_pdfd[81];

    auto g_x_xy_xyy_yz = buffer_pdfd[82];

    auto g_x_xy_xyy_zz = buffer_pdfd[83];

    auto g_x_xy_xyz_xx = buffer_pdfd[84];

    auto g_x_xy_xyz_xy = buffer_pdfd[85];

    auto g_x_xy_xyz_xz = buffer_pdfd[86];

    auto g_x_xy_xyz_yy = buffer_pdfd[87];

    auto g_x_xy_xyz_yz = buffer_pdfd[88];

    auto g_x_xy_xyz_zz = buffer_pdfd[89];

    auto g_x_xy_xzz_xx = buffer_pdfd[90];

    auto g_x_xy_xzz_xy = buffer_pdfd[91];

    auto g_x_xy_xzz_xz = buffer_pdfd[92];

    auto g_x_xy_xzz_yy = buffer_pdfd[93];

    auto g_x_xy_xzz_yz = buffer_pdfd[94];

    auto g_x_xy_xzz_zz = buffer_pdfd[95];

    auto g_x_xy_yyy_xx = buffer_pdfd[96];

    auto g_x_xy_yyy_xy = buffer_pdfd[97];

    auto g_x_xy_yyy_xz = buffer_pdfd[98];

    auto g_x_xy_yyy_yy = buffer_pdfd[99];

    auto g_x_xy_yyy_yz = buffer_pdfd[100];

    auto g_x_xy_yyy_zz = buffer_pdfd[101];

    auto g_x_xy_yyz_xx = buffer_pdfd[102];

    auto g_x_xy_yyz_xy = buffer_pdfd[103];

    auto g_x_xy_yyz_xz = buffer_pdfd[104];

    auto g_x_xy_yyz_yy = buffer_pdfd[105];

    auto g_x_xy_yyz_yz = buffer_pdfd[106];

    auto g_x_xy_yyz_zz = buffer_pdfd[107];

    auto g_x_xy_yzz_xx = buffer_pdfd[108];

    auto g_x_xy_yzz_xy = buffer_pdfd[109];

    auto g_x_xy_yzz_xz = buffer_pdfd[110];

    auto g_x_xy_yzz_yy = buffer_pdfd[111];

    auto g_x_xy_yzz_yz = buffer_pdfd[112];

    auto g_x_xy_yzz_zz = buffer_pdfd[113];

    auto g_x_xy_zzz_xx = buffer_pdfd[114];

    auto g_x_xy_zzz_xy = buffer_pdfd[115];

    auto g_x_xy_zzz_xz = buffer_pdfd[116];

    auto g_x_xy_zzz_yy = buffer_pdfd[117];

    auto g_x_xy_zzz_yz = buffer_pdfd[118];

    auto g_x_xy_zzz_zz = buffer_pdfd[119];

    auto g_x_xz_xxx_xx = buffer_pdfd[120];

    auto g_x_xz_xxx_xy = buffer_pdfd[121];

    auto g_x_xz_xxx_xz = buffer_pdfd[122];

    auto g_x_xz_xxx_yy = buffer_pdfd[123];

    auto g_x_xz_xxx_yz = buffer_pdfd[124];

    auto g_x_xz_xxx_zz = buffer_pdfd[125];

    auto g_x_xz_xxy_xx = buffer_pdfd[126];

    auto g_x_xz_xxy_xy = buffer_pdfd[127];

    auto g_x_xz_xxy_xz = buffer_pdfd[128];

    auto g_x_xz_xxy_yy = buffer_pdfd[129];

    auto g_x_xz_xxy_yz = buffer_pdfd[130];

    auto g_x_xz_xxy_zz = buffer_pdfd[131];

    auto g_x_xz_xxz_xx = buffer_pdfd[132];

    auto g_x_xz_xxz_xy = buffer_pdfd[133];

    auto g_x_xz_xxz_xz = buffer_pdfd[134];

    auto g_x_xz_xxz_yy = buffer_pdfd[135];

    auto g_x_xz_xxz_yz = buffer_pdfd[136];

    auto g_x_xz_xxz_zz = buffer_pdfd[137];

    auto g_x_xz_xyy_xx = buffer_pdfd[138];

    auto g_x_xz_xyy_xy = buffer_pdfd[139];

    auto g_x_xz_xyy_xz = buffer_pdfd[140];

    auto g_x_xz_xyy_yy = buffer_pdfd[141];

    auto g_x_xz_xyy_yz = buffer_pdfd[142];

    auto g_x_xz_xyy_zz = buffer_pdfd[143];

    auto g_x_xz_xyz_xx = buffer_pdfd[144];

    auto g_x_xz_xyz_xy = buffer_pdfd[145];

    auto g_x_xz_xyz_xz = buffer_pdfd[146];

    auto g_x_xz_xyz_yy = buffer_pdfd[147];

    auto g_x_xz_xyz_yz = buffer_pdfd[148];

    auto g_x_xz_xyz_zz = buffer_pdfd[149];

    auto g_x_xz_xzz_xx = buffer_pdfd[150];

    auto g_x_xz_xzz_xy = buffer_pdfd[151];

    auto g_x_xz_xzz_xz = buffer_pdfd[152];

    auto g_x_xz_xzz_yy = buffer_pdfd[153];

    auto g_x_xz_xzz_yz = buffer_pdfd[154];

    auto g_x_xz_xzz_zz = buffer_pdfd[155];

    auto g_x_xz_yyy_xx = buffer_pdfd[156];

    auto g_x_xz_yyy_xy = buffer_pdfd[157];

    auto g_x_xz_yyy_xz = buffer_pdfd[158];

    auto g_x_xz_yyy_yy = buffer_pdfd[159];

    auto g_x_xz_yyy_yz = buffer_pdfd[160];

    auto g_x_xz_yyy_zz = buffer_pdfd[161];

    auto g_x_xz_yyz_xx = buffer_pdfd[162];

    auto g_x_xz_yyz_xy = buffer_pdfd[163];

    auto g_x_xz_yyz_xz = buffer_pdfd[164];

    auto g_x_xz_yyz_yy = buffer_pdfd[165];

    auto g_x_xz_yyz_yz = buffer_pdfd[166];

    auto g_x_xz_yyz_zz = buffer_pdfd[167];

    auto g_x_xz_yzz_xx = buffer_pdfd[168];

    auto g_x_xz_yzz_xy = buffer_pdfd[169];

    auto g_x_xz_yzz_xz = buffer_pdfd[170];

    auto g_x_xz_yzz_yy = buffer_pdfd[171];

    auto g_x_xz_yzz_yz = buffer_pdfd[172];

    auto g_x_xz_yzz_zz = buffer_pdfd[173];

    auto g_x_xz_zzz_xx = buffer_pdfd[174];

    auto g_x_xz_zzz_xy = buffer_pdfd[175];

    auto g_x_xz_zzz_xz = buffer_pdfd[176];

    auto g_x_xz_zzz_yy = buffer_pdfd[177];

    auto g_x_xz_zzz_yz = buffer_pdfd[178];

    auto g_x_xz_zzz_zz = buffer_pdfd[179];

    auto g_x_yy_xxx_xx = buffer_pdfd[180];

    auto g_x_yy_xxx_xy = buffer_pdfd[181];

    auto g_x_yy_xxx_xz = buffer_pdfd[182];

    auto g_x_yy_xxx_yy = buffer_pdfd[183];

    auto g_x_yy_xxx_yz = buffer_pdfd[184];

    auto g_x_yy_xxx_zz = buffer_pdfd[185];

    auto g_x_yy_xxy_xx = buffer_pdfd[186];

    auto g_x_yy_xxy_xy = buffer_pdfd[187];

    auto g_x_yy_xxy_xz = buffer_pdfd[188];

    auto g_x_yy_xxy_yy = buffer_pdfd[189];

    auto g_x_yy_xxy_yz = buffer_pdfd[190];

    auto g_x_yy_xxy_zz = buffer_pdfd[191];

    auto g_x_yy_xxz_xx = buffer_pdfd[192];

    auto g_x_yy_xxz_xy = buffer_pdfd[193];

    auto g_x_yy_xxz_xz = buffer_pdfd[194];

    auto g_x_yy_xxz_yy = buffer_pdfd[195];

    auto g_x_yy_xxz_yz = buffer_pdfd[196];

    auto g_x_yy_xxz_zz = buffer_pdfd[197];

    auto g_x_yy_xyy_xx = buffer_pdfd[198];

    auto g_x_yy_xyy_xy = buffer_pdfd[199];

    auto g_x_yy_xyy_xz = buffer_pdfd[200];

    auto g_x_yy_xyy_yy = buffer_pdfd[201];

    auto g_x_yy_xyy_yz = buffer_pdfd[202];

    auto g_x_yy_xyy_zz = buffer_pdfd[203];

    auto g_x_yy_xyz_xx = buffer_pdfd[204];

    auto g_x_yy_xyz_xy = buffer_pdfd[205];

    auto g_x_yy_xyz_xz = buffer_pdfd[206];

    auto g_x_yy_xyz_yy = buffer_pdfd[207];

    auto g_x_yy_xyz_yz = buffer_pdfd[208];

    auto g_x_yy_xyz_zz = buffer_pdfd[209];

    auto g_x_yy_xzz_xx = buffer_pdfd[210];

    auto g_x_yy_xzz_xy = buffer_pdfd[211];

    auto g_x_yy_xzz_xz = buffer_pdfd[212];

    auto g_x_yy_xzz_yy = buffer_pdfd[213];

    auto g_x_yy_xzz_yz = buffer_pdfd[214];

    auto g_x_yy_xzz_zz = buffer_pdfd[215];

    auto g_x_yy_yyy_xx = buffer_pdfd[216];

    auto g_x_yy_yyy_xy = buffer_pdfd[217];

    auto g_x_yy_yyy_xz = buffer_pdfd[218];

    auto g_x_yy_yyy_yy = buffer_pdfd[219];

    auto g_x_yy_yyy_yz = buffer_pdfd[220];

    auto g_x_yy_yyy_zz = buffer_pdfd[221];

    auto g_x_yy_yyz_xx = buffer_pdfd[222];

    auto g_x_yy_yyz_xy = buffer_pdfd[223];

    auto g_x_yy_yyz_xz = buffer_pdfd[224];

    auto g_x_yy_yyz_yy = buffer_pdfd[225];

    auto g_x_yy_yyz_yz = buffer_pdfd[226];

    auto g_x_yy_yyz_zz = buffer_pdfd[227];

    auto g_x_yy_yzz_xx = buffer_pdfd[228];

    auto g_x_yy_yzz_xy = buffer_pdfd[229];

    auto g_x_yy_yzz_xz = buffer_pdfd[230];

    auto g_x_yy_yzz_yy = buffer_pdfd[231];

    auto g_x_yy_yzz_yz = buffer_pdfd[232];

    auto g_x_yy_yzz_zz = buffer_pdfd[233];

    auto g_x_yy_zzz_xx = buffer_pdfd[234];

    auto g_x_yy_zzz_xy = buffer_pdfd[235];

    auto g_x_yy_zzz_xz = buffer_pdfd[236];

    auto g_x_yy_zzz_yy = buffer_pdfd[237];

    auto g_x_yy_zzz_yz = buffer_pdfd[238];

    auto g_x_yy_zzz_zz = buffer_pdfd[239];

    auto g_x_yz_xxx_xx = buffer_pdfd[240];

    auto g_x_yz_xxx_xy = buffer_pdfd[241];

    auto g_x_yz_xxx_xz = buffer_pdfd[242];

    auto g_x_yz_xxx_yy = buffer_pdfd[243];

    auto g_x_yz_xxx_yz = buffer_pdfd[244];

    auto g_x_yz_xxx_zz = buffer_pdfd[245];

    auto g_x_yz_xxy_xx = buffer_pdfd[246];

    auto g_x_yz_xxy_xy = buffer_pdfd[247];

    auto g_x_yz_xxy_xz = buffer_pdfd[248];

    auto g_x_yz_xxy_yy = buffer_pdfd[249];

    auto g_x_yz_xxy_yz = buffer_pdfd[250];

    auto g_x_yz_xxy_zz = buffer_pdfd[251];

    auto g_x_yz_xxz_xx = buffer_pdfd[252];

    auto g_x_yz_xxz_xy = buffer_pdfd[253];

    auto g_x_yz_xxz_xz = buffer_pdfd[254];

    auto g_x_yz_xxz_yy = buffer_pdfd[255];

    auto g_x_yz_xxz_yz = buffer_pdfd[256];

    auto g_x_yz_xxz_zz = buffer_pdfd[257];

    auto g_x_yz_xyy_xx = buffer_pdfd[258];

    auto g_x_yz_xyy_xy = buffer_pdfd[259];

    auto g_x_yz_xyy_xz = buffer_pdfd[260];

    auto g_x_yz_xyy_yy = buffer_pdfd[261];

    auto g_x_yz_xyy_yz = buffer_pdfd[262];

    auto g_x_yz_xyy_zz = buffer_pdfd[263];

    auto g_x_yz_xyz_xx = buffer_pdfd[264];

    auto g_x_yz_xyz_xy = buffer_pdfd[265];

    auto g_x_yz_xyz_xz = buffer_pdfd[266];

    auto g_x_yz_xyz_yy = buffer_pdfd[267];

    auto g_x_yz_xyz_yz = buffer_pdfd[268];

    auto g_x_yz_xyz_zz = buffer_pdfd[269];

    auto g_x_yz_xzz_xx = buffer_pdfd[270];

    auto g_x_yz_xzz_xy = buffer_pdfd[271];

    auto g_x_yz_xzz_xz = buffer_pdfd[272];

    auto g_x_yz_xzz_yy = buffer_pdfd[273];

    auto g_x_yz_xzz_yz = buffer_pdfd[274];

    auto g_x_yz_xzz_zz = buffer_pdfd[275];

    auto g_x_yz_yyy_xx = buffer_pdfd[276];

    auto g_x_yz_yyy_xy = buffer_pdfd[277];

    auto g_x_yz_yyy_xz = buffer_pdfd[278];

    auto g_x_yz_yyy_yy = buffer_pdfd[279];

    auto g_x_yz_yyy_yz = buffer_pdfd[280];

    auto g_x_yz_yyy_zz = buffer_pdfd[281];

    auto g_x_yz_yyz_xx = buffer_pdfd[282];

    auto g_x_yz_yyz_xy = buffer_pdfd[283];

    auto g_x_yz_yyz_xz = buffer_pdfd[284];

    auto g_x_yz_yyz_yy = buffer_pdfd[285];

    auto g_x_yz_yyz_yz = buffer_pdfd[286];

    auto g_x_yz_yyz_zz = buffer_pdfd[287];

    auto g_x_yz_yzz_xx = buffer_pdfd[288];

    auto g_x_yz_yzz_xy = buffer_pdfd[289];

    auto g_x_yz_yzz_xz = buffer_pdfd[290];

    auto g_x_yz_yzz_yy = buffer_pdfd[291];

    auto g_x_yz_yzz_yz = buffer_pdfd[292];

    auto g_x_yz_yzz_zz = buffer_pdfd[293];

    auto g_x_yz_zzz_xx = buffer_pdfd[294];

    auto g_x_yz_zzz_xy = buffer_pdfd[295];

    auto g_x_yz_zzz_xz = buffer_pdfd[296];

    auto g_x_yz_zzz_yy = buffer_pdfd[297];

    auto g_x_yz_zzz_yz = buffer_pdfd[298];

    auto g_x_yz_zzz_zz = buffer_pdfd[299];

    auto g_x_zz_xxx_xx = buffer_pdfd[300];

    auto g_x_zz_xxx_xy = buffer_pdfd[301];

    auto g_x_zz_xxx_xz = buffer_pdfd[302];

    auto g_x_zz_xxx_yy = buffer_pdfd[303];

    auto g_x_zz_xxx_yz = buffer_pdfd[304];

    auto g_x_zz_xxx_zz = buffer_pdfd[305];

    auto g_x_zz_xxy_xx = buffer_pdfd[306];

    auto g_x_zz_xxy_xy = buffer_pdfd[307];

    auto g_x_zz_xxy_xz = buffer_pdfd[308];

    auto g_x_zz_xxy_yy = buffer_pdfd[309];

    auto g_x_zz_xxy_yz = buffer_pdfd[310];

    auto g_x_zz_xxy_zz = buffer_pdfd[311];

    auto g_x_zz_xxz_xx = buffer_pdfd[312];

    auto g_x_zz_xxz_xy = buffer_pdfd[313];

    auto g_x_zz_xxz_xz = buffer_pdfd[314];

    auto g_x_zz_xxz_yy = buffer_pdfd[315];

    auto g_x_zz_xxz_yz = buffer_pdfd[316];

    auto g_x_zz_xxz_zz = buffer_pdfd[317];

    auto g_x_zz_xyy_xx = buffer_pdfd[318];

    auto g_x_zz_xyy_xy = buffer_pdfd[319];

    auto g_x_zz_xyy_xz = buffer_pdfd[320];

    auto g_x_zz_xyy_yy = buffer_pdfd[321];

    auto g_x_zz_xyy_yz = buffer_pdfd[322];

    auto g_x_zz_xyy_zz = buffer_pdfd[323];

    auto g_x_zz_xyz_xx = buffer_pdfd[324];

    auto g_x_zz_xyz_xy = buffer_pdfd[325];

    auto g_x_zz_xyz_xz = buffer_pdfd[326];

    auto g_x_zz_xyz_yy = buffer_pdfd[327];

    auto g_x_zz_xyz_yz = buffer_pdfd[328];

    auto g_x_zz_xyz_zz = buffer_pdfd[329];

    auto g_x_zz_xzz_xx = buffer_pdfd[330];

    auto g_x_zz_xzz_xy = buffer_pdfd[331];

    auto g_x_zz_xzz_xz = buffer_pdfd[332];

    auto g_x_zz_xzz_yy = buffer_pdfd[333];

    auto g_x_zz_xzz_yz = buffer_pdfd[334];

    auto g_x_zz_xzz_zz = buffer_pdfd[335];

    auto g_x_zz_yyy_xx = buffer_pdfd[336];

    auto g_x_zz_yyy_xy = buffer_pdfd[337];

    auto g_x_zz_yyy_xz = buffer_pdfd[338];

    auto g_x_zz_yyy_yy = buffer_pdfd[339];

    auto g_x_zz_yyy_yz = buffer_pdfd[340];

    auto g_x_zz_yyy_zz = buffer_pdfd[341];

    auto g_x_zz_yyz_xx = buffer_pdfd[342];

    auto g_x_zz_yyz_xy = buffer_pdfd[343];

    auto g_x_zz_yyz_xz = buffer_pdfd[344];

    auto g_x_zz_yyz_yy = buffer_pdfd[345];

    auto g_x_zz_yyz_yz = buffer_pdfd[346];

    auto g_x_zz_yyz_zz = buffer_pdfd[347];

    auto g_x_zz_yzz_xx = buffer_pdfd[348];

    auto g_x_zz_yzz_xy = buffer_pdfd[349];

    auto g_x_zz_yzz_xz = buffer_pdfd[350];

    auto g_x_zz_yzz_yy = buffer_pdfd[351];

    auto g_x_zz_yzz_yz = buffer_pdfd[352];

    auto g_x_zz_yzz_zz = buffer_pdfd[353];

    auto g_x_zz_zzz_xx = buffer_pdfd[354];

    auto g_x_zz_zzz_xy = buffer_pdfd[355];

    auto g_x_zz_zzz_xz = buffer_pdfd[356];

    auto g_x_zz_zzz_yy = buffer_pdfd[357];

    auto g_x_zz_zzz_yz = buffer_pdfd[358];

    auto g_x_zz_zzz_zz = buffer_pdfd[359];

    auto g_y_xx_xxx_xx = buffer_pdfd[360];

    auto g_y_xx_xxx_xy = buffer_pdfd[361];

    auto g_y_xx_xxx_xz = buffer_pdfd[362];

    auto g_y_xx_xxx_yy = buffer_pdfd[363];

    auto g_y_xx_xxx_yz = buffer_pdfd[364];

    auto g_y_xx_xxx_zz = buffer_pdfd[365];

    auto g_y_xx_xxy_xx = buffer_pdfd[366];

    auto g_y_xx_xxy_xy = buffer_pdfd[367];

    auto g_y_xx_xxy_xz = buffer_pdfd[368];

    auto g_y_xx_xxy_yy = buffer_pdfd[369];

    auto g_y_xx_xxy_yz = buffer_pdfd[370];

    auto g_y_xx_xxy_zz = buffer_pdfd[371];

    auto g_y_xx_xxz_xx = buffer_pdfd[372];

    auto g_y_xx_xxz_xy = buffer_pdfd[373];

    auto g_y_xx_xxz_xz = buffer_pdfd[374];

    auto g_y_xx_xxz_yy = buffer_pdfd[375];

    auto g_y_xx_xxz_yz = buffer_pdfd[376];

    auto g_y_xx_xxz_zz = buffer_pdfd[377];

    auto g_y_xx_xyy_xx = buffer_pdfd[378];

    auto g_y_xx_xyy_xy = buffer_pdfd[379];

    auto g_y_xx_xyy_xz = buffer_pdfd[380];

    auto g_y_xx_xyy_yy = buffer_pdfd[381];

    auto g_y_xx_xyy_yz = buffer_pdfd[382];

    auto g_y_xx_xyy_zz = buffer_pdfd[383];

    auto g_y_xx_xyz_xx = buffer_pdfd[384];

    auto g_y_xx_xyz_xy = buffer_pdfd[385];

    auto g_y_xx_xyz_xz = buffer_pdfd[386];

    auto g_y_xx_xyz_yy = buffer_pdfd[387];

    auto g_y_xx_xyz_yz = buffer_pdfd[388];

    auto g_y_xx_xyz_zz = buffer_pdfd[389];

    auto g_y_xx_xzz_xx = buffer_pdfd[390];

    auto g_y_xx_xzz_xy = buffer_pdfd[391];

    auto g_y_xx_xzz_xz = buffer_pdfd[392];

    auto g_y_xx_xzz_yy = buffer_pdfd[393];

    auto g_y_xx_xzz_yz = buffer_pdfd[394];

    auto g_y_xx_xzz_zz = buffer_pdfd[395];

    auto g_y_xx_yyy_xx = buffer_pdfd[396];

    auto g_y_xx_yyy_xy = buffer_pdfd[397];

    auto g_y_xx_yyy_xz = buffer_pdfd[398];

    auto g_y_xx_yyy_yy = buffer_pdfd[399];

    auto g_y_xx_yyy_yz = buffer_pdfd[400];

    auto g_y_xx_yyy_zz = buffer_pdfd[401];

    auto g_y_xx_yyz_xx = buffer_pdfd[402];

    auto g_y_xx_yyz_xy = buffer_pdfd[403];

    auto g_y_xx_yyz_xz = buffer_pdfd[404];

    auto g_y_xx_yyz_yy = buffer_pdfd[405];

    auto g_y_xx_yyz_yz = buffer_pdfd[406];

    auto g_y_xx_yyz_zz = buffer_pdfd[407];

    auto g_y_xx_yzz_xx = buffer_pdfd[408];

    auto g_y_xx_yzz_xy = buffer_pdfd[409];

    auto g_y_xx_yzz_xz = buffer_pdfd[410];

    auto g_y_xx_yzz_yy = buffer_pdfd[411];

    auto g_y_xx_yzz_yz = buffer_pdfd[412];

    auto g_y_xx_yzz_zz = buffer_pdfd[413];

    auto g_y_xx_zzz_xx = buffer_pdfd[414];

    auto g_y_xx_zzz_xy = buffer_pdfd[415];

    auto g_y_xx_zzz_xz = buffer_pdfd[416];

    auto g_y_xx_zzz_yy = buffer_pdfd[417];

    auto g_y_xx_zzz_yz = buffer_pdfd[418];

    auto g_y_xx_zzz_zz = buffer_pdfd[419];

    auto g_y_xy_xxx_xx = buffer_pdfd[420];

    auto g_y_xy_xxx_xy = buffer_pdfd[421];

    auto g_y_xy_xxx_xz = buffer_pdfd[422];

    auto g_y_xy_xxx_yy = buffer_pdfd[423];

    auto g_y_xy_xxx_yz = buffer_pdfd[424];

    auto g_y_xy_xxx_zz = buffer_pdfd[425];

    auto g_y_xy_xxy_xx = buffer_pdfd[426];

    auto g_y_xy_xxy_xy = buffer_pdfd[427];

    auto g_y_xy_xxy_xz = buffer_pdfd[428];

    auto g_y_xy_xxy_yy = buffer_pdfd[429];

    auto g_y_xy_xxy_yz = buffer_pdfd[430];

    auto g_y_xy_xxy_zz = buffer_pdfd[431];

    auto g_y_xy_xxz_xx = buffer_pdfd[432];

    auto g_y_xy_xxz_xy = buffer_pdfd[433];

    auto g_y_xy_xxz_xz = buffer_pdfd[434];

    auto g_y_xy_xxz_yy = buffer_pdfd[435];

    auto g_y_xy_xxz_yz = buffer_pdfd[436];

    auto g_y_xy_xxz_zz = buffer_pdfd[437];

    auto g_y_xy_xyy_xx = buffer_pdfd[438];

    auto g_y_xy_xyy_xy = buffer_pdfd[439];

    auto g_y_xy_xyy_xz = buffer_pdfd[440];

    auto g_y_xy_xyy_yy = buffer_pdfd[441];

    auto g_y_xy_xyy_yz = buffer_pdfd[442];

    auto g_y_xy_xyy_zz = buffer_pdfd[443];

    auto g_y_xy_xyz_xx = buffer_pdfd[444];

    auto g_y_xy_xyz_xy = buffer_pdfd[445];

    auto g_y_xy_xyz_xz = buffer_pdfd[446];

    auto g_y_xy_xyz_yy = buffer_pdfd[447];

    auto g_y_xy_xyz_yz = buffer_pdfd[448];

    auto g_y_xy_xyz_zz = buffer_pdfd[449];

    auto g_y_xy_xzz_xx = buffer_pdfd[450];

    auto g_y_xy_xzz_xy = buffer_pdfd[451];

    auto g_y_xy_xzz_xz = buffer_pdfd[452];

    auto g_y_xy_xzz_yy = buffer_pdfd[453];

    auto g_y_xy_xzz_yz = buffer_pdfd[454];

    auto g_y_xy_xzz_zz = buffer_pdfd[455];

    auto g_y_xy_yyy_xx = buffer_pdfd[456];

    auto g_y_xy_yyy_xy = buffer_pdfd[457];

    auto g_y_xy_yyy_xz = buffer_pdfd[458];

    auto g_y_xy_yyy_yy = buffer_pdfd[459];

    auto g_y_xy_yyy_yz = buffer_pdfd[460];

    auto g_y_xy_yyy_zz = buffer_pdfd[461];

    auto g_y_xy_yyz_xx = buffer_pdfd[462];

    auto g_y_xy_yyz_xy = buffer_pdfd[463];

    auto g_y_xy_yyz_xz = buffer_pdfd[464];

    auto g_y_xy_yyz_yy = buffer_pdfd[465];

    auto g_y_xy_yyz_yz = buffer_pdfd[466];

    auto g_y_xy_yyz_zz = buffer_pdfd[467];

    auto g_y_xy_yzz_xx = buffer_pdfd[468];

    auto g_y_xy_yzz_xy = buffer_pdfd[469];

    auto g_y_xy_yzz_xz = buffer_pdfd[470];

    auto g_y_xy_yzz_yy = buffer_pdfd[471];

    auto g_y_xy_yzz_yz = buffer_pdfd[472];

    auto g_y_xy_yzz_zz = buffer_pdfd[473];

    auto g_y_xy_zzz_xx = buffer_pdfd[474];

    auto g_y_xy_zzz_xy = buffer_pdfd[475];

    auto g_y_xy_zzz_xz = buffer_pdfd[476];

    auto g_y_xy_zzz_yy = buffer_pdfd[477];

    auto g_y_xy_zzz_yz = buffer_pdfd[478];

    auto g_y_xy_zzz_zz = buffer_pdfd[479];

    auto g_y_xz_xxx_xx = buffer_pdfd[480];

    auto g_y_xz_xxx_xy = buffer_pdfd[481];

    auto g_y_xz_xxx_xz = buffer_pdfd[482];

    auto g_y_xz_xxx_yy = buffer_pdfd[483];

    auto g_y_xz_xxx_yz = buffer_pdfd[484];

    auto g_y_xz_xxx_zz = buffer_pdfd[485];

    auto g_y_xz_xxy_xx = buffer_pdfd[486];

    auto g_y_xz_xxy_xy = buffer_pdfd[487];

    auto g_y_xz_xxy_xz = buffer_pdfd[488];

    auto g_y_xz_xxy_yy = buffer_pdfd[489];

    auto g_y_xz_xxy_yz = buffer_pdfd[490];

    auto g_y_xz_xxy_zz = buffer_pdfd[491];

    auto g_y_xz_xxz_xx = buffer_pdfd[492];

    auto g_y_xz_xxz_xy = buffer_pdfd[493];

    auto g_y_xz_xxz_xz = buffer_pdfd[494];

    auto g_y_xz_xxz_yy = buffer_pdfd[495];

    auto g_y_xz_xxz_yz = buffer_pdfd[496];

    auto g_y_xz_xxz_zz = buffer_pdfd[497];

    auto g_y_xz_xyy_xx = buffer_pdfd[498];

    auto g_y_xz_xyy_xy = buffer_pdfd[499];

    auto g_y_xz_xyy_xz = buffer_pdfd[500];

    auto g_y_xz_xyy_yy = buffer_pdfd[501];

    auto g_y_xz_xyy_yz = buffer_pdfd[502];

    auto g_y_xz_xyy_zz = buffer_pdfd[503];

    auto g_y_xz_xyz_xx = buffer_pdfd[504];

    auto g_y_xz_xyz_xy = buffer_pdfd[505];

    auto g_y_xz_xyz_xz = buffer_pdfd[506];

    auto g_y_xz_xyz_yy = buffer_pdfd[507];

    auto g_y_xz_xyz_yz = buffer_pdfd[508];

    auto g_y_xz_xyz_zz = buffer_pdfd[509];

    auto g_y_xz_xzz_xx = buffer_pdfd[510];

    auto g_y_xz_xzz_xy = buffer_pdfd[511];

    auto g_y_xz_xzz_xz = buffer_pdfd[512];

    auto g_y_xz_xzz_yy = buffer_pdfd[513];

    auto g_y_xz_xzz_yz = buffer_pdfd[514];

    auto g_y_xz_xzz_zz = buffer_pdfd[515];

    auto g_y_xz_yyy_xx = buffer_pdfd[516];

    auto g_y_xz_yyy_xy = buffer_pdfd[517];

    auto g_y_xz_yyy_xz = buffer_pdfd[518];

    auto g_y_xz_yyy_yy = buffer_pdfd[519];

    auto g_y_xz_yyy_yz = buffer_pdfd[520];

    auto g_y_xz_yyy_zz = buffer_pdfd[521];

    auto g_y_xz_yyz_xx = buffer_pdfd[522];

    auto g_y_xz_yyz_xy = buffer_pdfd[523];

    auto g_y_xz_yyz_xz = buffer_pdfd[524];

    auto g_y_xz_yyz_yy = buffer_pdfd[525];

    auto g_y_xz_yyz_yz = buffer_pdfd[526];

    auto g_y_xz_yyz_zz = buffer_pdfd[527];

    auto g_y_xz_yzz_xx = buffer_pdfd[528];

    auto g_y_xz_yzz_xy = buffer_pdfd[529];

    auto g_y_xz_yzz_xz = buffer_pdfd[530];

    auto g_y_xz_yzz_yy = buffer_pdfd[531];

    auto g_y_xz_yzz_yz = buffer_pdfd[532];

    auto g_y_xz_yzz_zz = buffer_pdfd[533];

    auto g_y_xz_zzz_xx = buffer_pdfd[534];

    auto g_y_xz_zzz_xy = buffer_pdfd[535];

    auto g_y_xz_zzz_xz = buffer_pdfd[536];

    auto g_y_xz_zzz_yy = buffer_pdfd[537];

    auto g_y_xz_zzz_yz = buffer_pdfd[538];

    auto g_y_xz_zzz_zz = buffer_pdfd[539];

    auto g_y_yy_xxx_xx = buffer_pdfd[540];

    auto g_y_yy_xxx_xy = buffer_pdfd[541];

    auto g_y_yy_xxx_xz = buffer_pdfd[542];

    auto g_y_yy_xxx_yy = buffer_pdfd[543];

    auto g_y_yy_xxx_yz = buffer_pdfd[544];

    auto g_y_yy_xxx_zz = buffer_pdfd[545];

    auto g_y_yy_xxy_xx = buffer_pdfd[546];

    auto g_y_yy_xxy_xy = buffer_pdfd[547];

    auto g_y_yy_xxy_xz = buffer_pdfd[548];

    auto g_y_yy_xxy_yy = buffer_pdfd[549];

    auto g_y_yy_xxy_yz = buffer_pdfd[550];

    auto g_y_yy_xxy_zz = buffer_pdfd[551];

    auto g_y_yy_xxz_xx = buffer_pdfd[552];

    auto g_y_yy_xxz_xy = buffer_pdfd[553];

    auto g_y_yy_xxz_xz = buffer_pdfd[554];

    auto g_y_yy_xxz_yy = buffer_pdfd[555];

    auto g_y_yy_xxz_yz = buffer_pdfd[556];

    auto g_y_yy_xxz_zz = buffer_pdfd[557];

    auto g_y_yy_xyy_xx = buffer_pdfd[558];

    auto g_y_yy_xyy_xy = buffer_pdfd[559];

    auto g_y_yy_xyy_xz = buffer_pdfd[560];

    auto g_y_yy_xyy_yy = buffer_pdfd[561];

    auto g_y_yy_xyy_yz = buffer_pdfd[562];

    auto g_y_yy_xyy_zz = buffer_pdfd[563];

    auto g_y_yy_xyz_xx = buffer_pdfd[564];

    auto g_y_yy_xyz_xy = buffer_pdfd[565];

    auto g_y_yy_xyz_xz = buffer_pdfd[566];

    auto g_y_yy_xyz_yy = buffer_pdfd[567];

    auto g_y_yy_xyz_yz = buffer_pdfd[568];

    auto g_y_yy_xyz_zz = buffer_pdfd[569];

    auto g_y_yy_xzz_xx = buffer_pdfd[570];

    auto g_y_yy_xzz_xy = buffer_pdfd[571];

    auto g_y_yy_xzz_xz = buffer_pdfd[572];

    auto g_y_yy_xzz_yy = buffer_pdfd[573];

    auto g_y_yy_xzz_yz = buffer_pdfd[574];

    auto g_y_yy_xzz_zz = buffer_pdfd[575];

    auto g_y_yy_yyy_xx = buffer_pdfd[576];

    auto g_y_yy_yyy_xy = buffer_pdfd[577];

    auto g_y_yy_yyy_xz = buffer_pdfd[578];

    auto g_y_yy_yyy_yy = buffer_pdfd[579];

    auto g_y_yy_yyy_yz = buffer_pdfd[580];

    auto g_y_yy_yyy_zz = buffer_pdfd[581];

    auto g_y_yy_yyz_xx = buffer_pdfd[582];

    auto g_y_yy_yyz_xy = buffer_pdfd[583];

    auto g_y_yy_yyz_xz = buffer_pdfd[584];

    auto g_y_yy_yyz_yy = buffer_pdfd[585];

    auto g_y_yy_yyz_yz = buffer_pdfd[586];

    auto g_y_yy_yyz_zz = buffer_pdfd[587];

    auto g_y_yy_yzz_xx = buffer_pdfd[588];

    auto g_y_yy_yzz_xy = buffer_pdfd[589];

    auto g_y_yy_yzz_xz = buffer_pdfd[590];

    auto g_y_yy_yzz_yy = buffer_pdfd[591];

    auto g_y_yy_yzz_yz = buffer_pdfd[592];

    auto g_y_yy_yzz_zz = buffer_pdfd[593];

    auto g_y_yy_zzz_xx = buffer_pdfd[594];

    auto g_y_yy_zzz_xy = buffer_pdfd[595];

    auto g_y_yy_zzz_xz = buffer_pdfd[596];

    auto g_y_yy_zzz_yy = buffer_pdfd[597];

    auto g_y_yy_zzz_yz = buffer_pdfd[598];

    auto g_y_yy_zzz_zz = buffer_pdfd[599];

    auto g_y_yz_xxx_xx = buffer_pdfd[600];

    auto g_y_yz_xxx_xy = buffer_pdfd[601];

    auto g_y_yz_xxx_xz = buffer_pdfd[602];

    auto g_y_yz_xxx_yy = buffer_pdfd[603];

    auto g_y_yz_xxx_yz = buffer_pdfd[604];

    auto g_y_yz_xxx_zz = buffer_pdfd[605];

    auto g_y_yz_xxy_xx = buffer_pdfd[606];

    auto g_y_yz_xxy_xy = buffer_pdfd[607];

    auto g_y_yz_xxy_xz = buffer_pdfd[608];

    auto g_y_yz_xxy_yy = buffer_pdfd[609];

    auto g_y_yz_xxy_yz = buffer_pdfd[610];

    auto g_y_yz_xxy_zz = buffer_pdfd[611];

    auto g_y_yz_xxz_xx = buffer_pdfd[612];

    auto g_y_yz_xxz_xy = buffer_pdfd[613];

    auto g_y_yz_xxz_xz = buffer_pdfd[614];

    auto g_y_yz_xxz_yy = buffer_pdfd[615];

    auto g_y_yz_xxz_yz = buffer_pdfd[616];

    auto g_y_yz_xxz_zz = buffer_pdfd[617];

    auto g_y_yz_xyy_xx = buffer_pdfd[618];

    auto g_y_yz_xyy_xy = buffer_pdfd[619];

    auto g_y_yz_xyy_xz = buffer_pdfd[620];

    auto g_y_yz_xyy_yy = buffer_pdfd[621];

    auto g_y_yz_xyy_yz = buffer_pdfd[622];

    auto g_y_yz_xyy_zz = buffer_pdfd[623];

    auto g_y_yz_xyz_xx = buffer_pdfd[624];

    auto g_y_yz_xyz_xy = buffer_pdfd[625];

    auto g_y_yz_xyz_xz = buffer_pdfd[626];

    auto g_y_yz_xyz_yy = buffer_pdfd[627];

    auto g_y_yz_xyz_yz = buffer_pdfd[628];

    auto g_y_yz_xyz_zz = buffer_pdfd[629];

    auto g_y_yz_xzz_xx = buffer_pdfd[630];

    auto g_y_yz_xzz_xy = buffer_pdfd[631];

    auto g_y_yz_xzz_xz = buffer_pdfd[632];

    auto g_y_yz_xzz_yy = buffer_pdfd[633];

    auto g_y_yz_xzz_yz = buffer_pdfd[634];

    auto g_y_yz_xzz_zz = buffer_pdfd[635];

    auto g_y_yz_yyy_xx = buffer_pdfd[636];

    auto g_y_yz_yyy_xy = buffer_pdfd[637];

    auto g_y_yz_yyy_xz = buffer_pdfd[638];

    auto g_y_yz_yyy_yy = buffer_pdfd[639];

    auto g_y_yz_yyy_yz = buffer_pdfd[640];

    auto g_y_yz_yyy_zz = buffer_pdfd[641];

    auto g_y_yz_yyz_xx = buffer_pdfd[642];

    auto g_y_yz_yyz_xy = buffer_pdfd[643];

    auto g_y_yz_yyz_xz = buffer_pdfd[644];

    auto g_y_yz_yyz_yy = buffer_pdfd[645];

    auto g_y_yz_yyz_yz = buffer_pdfd[646];

    auto g_y_yz_yyz_zz = buffer_pdfd[647];

    auto g_y_yz_yzz_xx = buffer_pdfd[648];

    auto g_y_yz_yzz_xy = buffer_pdfd[649];

    auto g_y_yz_yzz_xz = buffer_pdfd[650];

    auto g_y_yz_yzz_yy = buffer_pdfd[651];

    auto g_y_yz_yzz_yz = buffer_pdfd[652];

    auto g_y_yz_yzz_zz = buffer_pdfd[653];

    auto g_y_yz_zzz_xx = buffer_pdfd[654];

    auto g_y_yz_zzz_xy = buffer_pdfd[655];

    auto g_y_yz_zzz_xz = buffer_pdfd[656];

    auto g_y_yz_zzz_yy = buffer_pdfd[657];

    auto g_y_yz_zzz_yz = buffer_pdfd[658];

    auto g_y_yz_zzz_zz = buffer_pdfd[659];

    auto g_y_zz_xxx_xx = buffer_pdfd[660];

    auto g_y_zz_xxx_xy = buffer_pdfd[661];

    auto g_y_zz_xxx_xz = buffer_pdfd[662];

    auto g_y_zz_xxx_yy = buffer_pdfd[663];

    auto g_y_zz_xxx_yz = buffer_pdfd[664];

    auto g_y_zz_xxx_zz = buffer_pdfd[665];

    auto g_y_zz_xxy_xx = buffer_pdfd[666];

    auto g_y_zz_xxy_xy = buffer_pdfd[667];

    auto g_y_zz_xxy_xz = buffer_pdfd[668];

    auto g_y_zz_xxy_yy = buffer_pdfd[669];

    auto g_y_zz_xxy_yz = buffer_pdfd[670];

    auto g_y_zz_xxy_zz = buffer_pdfd[671];

    auto g_y_zz_xxz_xx = buffer_pdfd[672];

    auto g_y_zz_xxz_xy = buffer_pdfd[673];

    auto g_y_zz_xxz_xz = buffer_pdfd[674];

    auto g_y_zz_xxz_yy = buffer_pdfd[675];

    auto g_y_zz_xxz_yz = buffer_pdfd[676];

    auto g_y_zz_xxz_zz = buffer_pdfd[677];

    auto g_y_zz_xyy_xx = buffer_pdfd[678];

    auto g_y_zz_xyy_xy = buffer_pdfd[679];

    auto g_y_zz_xyy_xz = buffer_pdfd[680];

    auto g_y_zz_xyy_yy = buffer_pdfd[681];

    auto g_y_zz_xyy_yz = buffer_pdfd[682];

    auto g_y_zz_xyy_zz = buffer_pdfd[683];

    auto g_y_zz_xyz_xx = buffer_pdfd[684];

    auto g_y_zz_xyz_xy = buffer_pdfd[685];

    auto g_y_zz_xyz_xz = buffer_pdfd[686];

    auto g_y_zz_xyz_yy = buffer_pdfd[687];

    auto g_y_zz_xyz_yz = buffer_pdfd[688];

    auto g_y_zz_xyz_zz = buffer_pdfd[689];

    auto g_y_zz_xzz_xx = buffer_pdfd[690];

    auto g_y_zz_xzz_xy = buffer_pdfd[691];

    auto g_y_zz_xzz_xz = buffer_pdfd[692];

    auto g_y_zz_xzz_yy = buffer_pdfd[693];

    auto g_y_zz_xzz_yz = buffer_pdfd[694];

    auto g_y_zz_xzz_zz = buffer_pdfd[695];

    auto g_y_zz_yyy_xx = buffer_pdfd[696];

    auto g_y_zz_yyy_xy = buffer_pdfd[697];

    auto g_y_zz_yyy_xz = buffer_pdfd[698];

    auto g_y_zz_yyy_yy = buffer_pdfd[699];

    auto g_y_zz_yyy_yz = buffer_pdfd[700];

    auto g_y_zz_yyy_zz = buffer_pdfd[701];

    auto g_y_zz_yyz_xx = buffer_pdfd[702];

    auto g_y_zz_yyz_xy = buffer_pdfd[703];

    auto g_y_zz_yyz_xz = buffer_pdfd[704];

    auto g_y_zz_yyz_yy = buffer_pdfd[705];

    auto g_y_zz_yyz_yz = buffer_pdfd[706];

    auto g_y_zz_yyz_zz = buffer_pdfd[707];

    auto g_y_zz_yzz_xx = buffer_pdfd[708];

    auto g_y_zz_yzz_xy = buffer_pdfd[709];

    auto g_y_zz_yzz_xz = buffer_pdfd[710];

    auto g_y_zz_yzz_yy = buffer_pdfd[711];

    auto g_y_zz_yzz_yz = buffer_pdfd[712];

    auto g_y_zz_yzz_zz = buffer_pdfd[713];

    auto g_y_zz_zzz_xx = buffer_pdfd[714];

    auto g_y_zz_zzz_xy = buffer_pdfd[715];

    auto g_y_zz_zzz_xz = buffer_pdfd[716];

    auto g_y_zz_zzz_yy = buffer_pdfd[717];

    auto g_y_zz_zzz_yz = buffer_pdfd[718];

    auto g_y_zz_zzz_zz = buffer_pdfd[719];

    auto g_z_xx_xxx_xx = buffer_pdfd[720];

    auto g_z_xx_xxx_xy = buffer_pdfd[721];

    auto g_z_xx_xxx_xz = buffer_pdfd[722];

    auto g_z_xx_xxx_yy = buffer_pdfd[723];

    auto g_z_xx_xxx_yz = buffer_pdfd[724];

    auto g_z_xx_xxx_zz = buffer_pdfd[725];

    auto g_z_xx_xxy_xx = buffer_pdfd[726];

    auto g_z_xx_xxy_xy = buffer_pdfd[727];

    auto g_z_xx_xxy_xz = buffer_pdfd[728];

    auto g_z_xx_xxy_yy = buffer_pdfd[729];

    auto g_z_xx_xxy_yz = buffer_pdfd[730];

    auto g_z_xx_xxy_zz = buffer_pdfd[731];

    auto g_z_xx_xxz_xx = buffer_pdfd[732];

    auto g_z_xx_xxz_xy = buffer_pdfd[733];

    auto g_z_xx_xxz_xz = buffer_pdfd[734];

    auto g_z_xx_xxz_yy = buffer_pdfd[735];

    auto g_z_xx_xxz_yz = buffer_pdfd[736];

    auto g_z_xx_xxz_zz = buffer_pdfd[737];

    auto g_z_xx_xyy_xx = buffer_pdfd[738];

    auto g_z_xx_xyy_xy = buffer_pdfd[739];

    auto g_z_xx_xyy_xz = buffer_pdfd[740];

    auto g_z_xx_xyy_yy = buffer_pdfd[741];

    auto g_z_xx_xyy_yz = buffer_pdfd[742];

    auto g_z_xx_xyy_zz = buffer_pdfd[743];

    auto g_z_xx_xyz_xx = buffer_pdfd[744];

    auto g_z_xx_xyz_xy = buffer_pdfd[745];

    auto g_z_xx_xyz_xz = buffer_pdfd[746];

    auto g_z_xx_xyz_yy = buffer_pdfd[747];

    auto g_z_xx_xyz_yz = buffer_pdfd[748];

    auto g_z_xx_xyz_zz = buffer_pdfd[749];

    auto g_z_xx_xzz_xx = buffer_pdfd[750];

    auto g_z_xx_xzz_xy = buffer_pdfd[751];

    auto g_z_xx_xzz_xz = buffer_pdfd[752];

    auto g_z_xx_xzz_yy = buffer_pdfd[753];

    auto g_z_xx_xzz_yz = buffer_pdfd[754];

    auto g_z_xx_xzz_zz = buffer_pdfd[755];

    auto g_z_xx_yyy_xx = buffer_pdfd[756];

    auto g_z_xx_yyy_xy = buffer_pdfd[757];

    auto g_z_xx_yyy_xz = buffer_pdfd[758];

    auto g_z_xx_yyy_yy = buffer_pdfd[759];

    auto g_z_xx_yyy_yz = buffer_pdfd[760];

    auto g_z_xx_yyy_zz = buffer_pdfd[761];

    auto g_z_xx_yyz_xx = buffer_pdfd[762];

    auto g_z_xx_yyz_xy = buffer_pdfd[763];

    auto g_z_xx_yyz_xz = buffer_pdfd[764];

    auto g_z_xx_yyz_yy = buffer_pdfd[765];

    auto g_z_xx_yyz_yz = buffer_pdfd[766];

    auto g_z_xx_yyz_zz = buffer_pdfd[767];

    auto g_z_xx_yzz_xx = buffer_pdfd[768];

    auto g_z_xx_yzz_xy = buffer_pdfd[769];

    auto g_z_xx_yzz_xz = buffer_pdfd[770];

    auto g_z_xx_yzz_yy = buffer_pdfd[771];

    auto g_z_xx_yzz_yz = buffer_pdfd[772];

    auto g_z_xx_yzz_zz = buffer_pdfd[773];

    auto g_z_xx_zzz_xx = buffer_pdfd[774];

    auto g_z_xx_zzz_xy = buffer_pdfd[775];

    auto g_z_xx_zzz_xz = buffer_pdfd[776];

    auto g_z_xx_zzz_yy = buffer_pdfd[777];

    auto g_z_xx_zzz_yz = buffer_pdfd[778];

    auto g_z_xx_zzz_zz = buffer_pdfd[779];

    auto g_z_xy_xxx_xx = buffer_pdfd[780];

    auto g_z_xy_xxx_xy = buffer_pdfd[781];

    auto g_z_xy_xxx_xz = buffer_pdfd[782];

    auto g_z_xy_xxx_yy = buffer_pdfd[783];

    auto g_z_xy_xxx_yz = buffer_pdfd[784];

    auto g_z_xy_xxx_zz = buffer_pdfd[785];

    auto g_z_xy_xxy_xx = buffer_pdfd[786];

    auto g_z_xy_xxy_xy = buffer_pdfd[787];

    auto g_z_xy_xxy_xz = buffer_pdfd[788];

    auto g_z_xy_xxy_yy = buffer_pdfd[789];

    auto g_z_xy_xxy_yz = buffer_pdfd[790];

    auto g_z_xy_xxy_zz = buffer_pdfd[791];

    auto g_z_xy_xxz_xx = buffer_pdfd[792];

    auto g_z_xy_xxz_xy = buffer_pdfd[793];

    auto g_z_xy_xxz_xz = buffer_pdfd[794];

    auto g_z_xy_xxz_yy = buffer_pdfd[795];

    auto g_z_xy_xxz_yz = buffer_pdfd[796];

    auto g_z_xy_xxz_zz = buffer_pdfd[797];

    auto g_z_xy_xyy_xx = buffer_pdfd[798];

    auto g_z_xy_xyy_xy = buffer_pdfd[799];

    auto g_z_xy_xyy_xz = buffer_pdfd[800];

    auto g_z_xy_xyy_yy = buffer_pdfd[801];

    auto g_z_xy_xyy_yz = buffer_pdfd[802];

    auto g_z_xy_xyy_zz = buffer_pdfd[803];

    auto g_z_xy_xyz_xx = buffer_pdfd[804];

    auto g_z_xy_xyz_xy = buffer_pdfd[805];

    auto g_z_xy_xyz_xz = buffer_pdfd[806];

    auto g_z_xy_xyz_yy = buffer_pdfd[807];

    auto g_z_xy_xyz_yz = buffer_pdfd[808];

    auto g_z_xy_xyz_zz = buffer_pdfd[809];

    auto g_z_xy_xzz_xx = buffer_pdfd[810];

    auto g_z_xy_xzz_xy = buffer_pdfd[811];

    auto g_z_xy_xzz_xz = buffer_pdfd[812];

    auto g_z_xy_xzz_yy = buffer_pdfd[813];

    auto g_z_xy_xzz_yz = buffer_pdfd[814];

    auto g_z_xy_xzz_zz = buffer_pdfd[815];

    auto g_z_xy_yyy_xx = buffer_pdfd[816];

    auto g_z_xy_yyy_xy = buffer_pdfd[817];

    auto g_z_xy_yyy_xz = buffer_pdfd[818];

    auto g_z_xy_yyy_yy = buffer_pdfd[819];

    auto g_z_xy_yyy_yz = buffer_pdfd[820];

    auto g_z_xy_yyy_zz = buffer_pdfd[821];

    auto g_z_xy_yyz_xx = buffer_pdfd[822];

    auto g_z_xy_yyz_xy = buffer_pdfd[823];

    auto g_z_xy_yyz_xz = buffer_pdfd[824];

    auto g_z_xy_yyz_yy = buffer_pdfd[825];

    auto g_z_xy_yyz_yz = buffer_pdfd[826];

    auto g_z_xy_yyz_zz = buffer_pdfd[827];

    auto g_z_xy_yzz_xx = buffer_pdfd[828];

    auto g_z_xy_yzz_xy = buffer_pdfd[829];

    auto g_z_xy_yzz_xz = buffer_pdfd[830];

    auto g_z_xy_yzz_yy = buffer_pdfd[831];

    auto g_z_xy_yzz_yz = buffer_pdfd[832];

    auto g_z_xy_yzz_zz = buffer_pdfd[833];

    auto g_z_xy_zzz_xx = buffer_pdfd[834];

    auto g_z_xy_zzz_xy = buffer_pdfd[835];

    auto g_z_xy_zzz_xz = buffer_pdfd[836];

    auto g_z_xy_zzz_yy = buffer_pdfd[837];

    auto g_z_xy_zzz_yz = buffer_pdfd[838];

    auto g_z_xy_zzz_zz = buffer_pdfd[839];

    auto g_z_xz_xxx_xx = buffer_pdfd[840];

    auto g_z_xz_xxx_xy = buffer_pdfd[841];

    auto g_z_xz_xxx_xz = buffer_pdfd[842];

    auto g_z_xz_xxx_yy = buffer_pdfd[843];

    auto g_z_xz_xxx_yz = buffer_pdfd[844];

    auto g_z_xz_xxx_zz = buffer_pdfd[845];

    auto g_z_xz_xxy_xx = buffer_pdfd[846];

    auto g_z_xz_xxy_xy = buffer_pdfd[847];

    auto g_z_xz_xxy_xz = buffer_pdfd[848];

    auto g_z_xz_xxy_yy = buffer_pdfd[849];

    auto g_z_xz_xxy_yz = buffer_pdfd[850];

    auto g_z_xz_xxy_zz = buffer_pdfd[851];

    auto g_z_xz_xxz_xx = buffer_pdfd[852];

    auto g_z_xz_xxz_xy = buffer_pdfd[853];

    auto g_z_xz_xxz_xz = buffer_pdfd[854];

    auto g_z_xz_xxz_yy = buffer_pdfd[855];

    auto g_z_xz_xxz_yz = buffer_pdfd[856];

    auto g_z_xz_xxz_zz = buffer_pdfd[857];

    auto g_z_xz_xyy_xx = buffer_pdfd[858];

    auto g_z_xz_xyy_xy = buffer_pdfd[859];

    auto g_z_xz_xyy_xz = buffer_pdfd[860];

    auto g_z_xz_xyy_yy = buffer_pdfd[861];

    auto g_z_xz_xyy_yz = buffer_pdfd[862];

    auto g_z_xz_xyy_zz = buffer_pdfd[863];

    auto g_z_xz_xyz_xx = buffer_pdfd[864];

    auto g_z_xz_xyz_xy = buffer_pdfd[865];

    auto g_z_xz_xyz_xz = buffer_pdfd[866];

    auto g_z_xz_xyz_yy = buffer_pdfd[867];

    auto g_z_xz_xyz_yz = buffer_pdfd[868];

    auto g_z_xz_xyz_zz = buffer_pdfd[869];

    auto g_z_xz_xzz_xx = buffer_pdfd[870];

    auto g_z_xz_xzz_xy = buffer_pdfd[871];

    auto g_z_xz_xzz_xz = buffer_pdfd[872];

    auto g_z_xz_xzz_yy = buffer_pdfd[873];

    auto g_z_xz_xzz_yz = buffer_pdfd[874];

    auto g_z_xz_xzz_zz = buffer_pdfd[875];

    auto g_z_xz_yyy_xx = buffer_pdfd[876];

    auto g_z_xz_yyy_xy = buffer_pdfd[877];

    auto g_z_xz_yyy_xz = buffer_pdfd[878];

    auto g_z_xz_yyy_yy = buffer_pdfd[879];

    auto g_z_xz_yyy_yz = buffer_pdfd[880];

    auto g_z_xz_yyy_zz = buffer_pdfd[881];

    auto g_z_xz_yyz_xx = buffer_pdfd[882];

    auto g_z_xz_yyz_xy = buffer_pdfd[883];

    auto g_z_xz_yyz_xz = buffer_pdfd[884];

    auto g_z_xz_yyz_yy = buffer_pdfd[885];

    auto g_z_xz_yyz_yz = buffer_pdfd[886];

    auto g_z_xz_yyz_zz = buffer_pdfd[887];

    auto g_z_xz_yzz_xx = buffer_pdfd[888];

    auto g_z_xz_yzz_xy = buffer_pdfd[889];

    auto g_z_xz_yzz_xz = buffer_pdfd[890];

    auto g_z_xz_yzz_yy = buffer_pdfd[891];

    auto g_z_xz_yzz_yz = buffer_pdfd[892];

    auto g_z_xz_yzz_zz = buffer_pdfd[893];

    auto g_z_xz_zzz_xx = buffer_pdfd[894];

    auto g_z_xz_zzz_xy = buffer_pdfd[895];

    auto g_z_xz_zzz_xz = buffer_pdfd[896];

    auto g_z_xz_zzz_yy = buffer_pdfd[897];

    auto g_z_xz_zzz_yz = buffer_pdfd[898];

    auto g_z_xz_zzz_zz = buffer_pdfd[899];

    auto g_z_yy_xxx_xx = buffer_pdfd[900];

    auto g_z_yy_xxx_xy = buffer_pdfd[901];

    auto g_z_yy_xxx_xz = buffer_pdfd[902];

    auto g_z_yy_xxx_yy = buffer_pdfd[903];

    auto g_z_yy_xxx_yz = buffer_pdfd[904];

    auto g_z_yy_xxx_zz = buffer_pdfd[905];

    auto g_z_yy_xxy_xx = buffer_pdfd[906];

    auto g_z_yy_xxy_xy = buffer_pdfd[907];

    auto g_z_yy_xxy_xz = buffer_pdfd[908];

    auto g_z_yy_xxy_yy = buffer_pdfd[909];

    auto g_z_yy_xxy_yz = buffer_pdfd[910];

    auto g_z_yy_xxy_zz = buffer_pdfd[911];

    auto g_z_yy_xxz_xx = buffer_pdfd[912];

    auto g_z_yy_xxz_xy = buffer_pdfd[913];

    auto g_z_yy_xxz_xz = buffer_pdfd[914];

    auto g_z_yy_xxz_yy = buffer_pdfd[915];

    auto g_z_yy_xxz_yz = buffer_pdfd[916];

    auto g_z_yy_xxz_zz = buffer_pdfd[917];

    auto g_z_yy_xyy_xx = buffer_pdfd[918];

    auto g_z_yy_xyy_xy = buffer_pdfd[919];

    auto g_z_yy_xyy_xz = buffer_pdfd[920];

    auto g_z_yy_xyy_yy = buffer_pdfd[921];

    auto g_z_yy_xyy_yz = buffer_pdfd[922];

    auto g_z_yy_xyy_zz = buffer_pdfd[923];

    auto g_z_yy_xyz_xx = buffer_pdfd[924];

    auto g_z_yy_xyz_xy = buffer_pdfd[925];

    auto g_z_yy_xyz_xz = buffer_pdfd[926];

    auto g_z_yy_xyz_yy = buffer_pdfd[927];

    auto g_z_yy_xyz_yz = buffer_pdfd[928];

    auto g_z_yy_xyz_zz = buffer_pdfd[929];

    auto g_z_yy_xzz_xx = buffer_pdfd[930];

    auto g_z_yy_xzz_xy = buffer_pdfd[931];

    auto g_z_yy_xzz_xz = buffer_pdfd[932];

    auto g_z_yy_xzz_yy = buffer_pdfd[933];

    auto g_z_yy_xzz_yz = buffer_pdfd[934];

    auto g_z_yy_xzz_zz = buffer_pdfd[935];

    auto g_z_yy_yyy_xx = buffer_pdfd[936];

    auto g_z_yy_yyy_xy = buffer_pdfd[937];

    auto g_z_yy_yyy_xz = buffer_pdfd[938];

    auto g_z_yy_yyy_yy = buffer_pdfd[939];

    auto g_z_yy_yyy_yz = buffer_pdfd[940];

    auto g_z_yy_yyy_zz = buffer_pdfd[941];

    auto g_z_yy_yyz_xx = buffer_pdfd[942];

    auto g_z_yy_yyz_xy = buffer_pdfd[943];

    auto g_z_yy_yyz_xz = buffer_pdfd[944];

    auto g_z_yy_yyz_yy = buffer_pdfd[945];

    auto g_z_yy_yyz_yz = buffer_pdfd[946];

    auto g_z_yy_yyz_zz = buffer_pdfd[947];

    auto g_z_yy_yzz_xx = buffer_pdfd[948];

    auto g_z_yy_yzz_xy = buffer_pdfd[949];

    auto g_z_yy_yzz_xz = buffer_pdfd[950];

    auto g_z_yy_yzz_yy = buffer_pdfd[951];

    auto g_z_yy_yzz_yz = buffer_pdfd[952];

    auto g_z_yy_yzz_zz = buffer_pdfd[953];

    auto g_z_yy_zzz_xx = buffer_pdfd[954];

    auto g_z_yy_zzz_xy = buffer_pdfd[955];

    auto g_z_yy_zzz_xz = buffer_pdfd[956];

    auto g_z_yy_zzz_yy = buffer_pdfd[957];

    auto g_z_yy_zzz_yz = buffer_pdfd[958];

    auto g_z_yy_zzz_zz = buffer_pdfd[959];

    auto g_z_yz_xxx_xx = buffer_pdfd[960];

    auto g_z_yz_xxx_xy = buffer_pdfd[961];

    auto g_z_yz_xxx_xz = buffer_pdfd[962];

    auto g_z_yz_xxx_yy = buffer_pdfd[963];

    auto g_z_yz_xxx_yz = buffer_pdfd[964];

    auto g_z_yz_xxx_zz = buffer_pdfd[965];

    auto g_z_yz_xxy_xx = buffer_pdfd[966];

    auto g_z_yz_xxy_xy = buffer_pdfd[967];

    auto g_z_yz_xxy_xz = buffer_pdfd[968];

    auto g_z_yz_xxy_yy = buffer_pdfd[969];

    auto g_z_yz_xxy_yz = buffer_pdfd[970];

    auto g_z_yz_xxy_zz = buffer_pdfd[971];

    auto g_z_yz_xxz_xx = buffer_pdfd[972];

    auto g_z_yz_xxz_xy = buffer_pdfd[973];

    auto g_z_yz_xxz_xz = buffer_pdfd[974];

    auto g_z_yz_xxz_yy = buffer_pdfd[975];

    auto g_z_yz_xxz_yz = buffer_pdfd[976];

    auto g_z_yz_xxz_zz = buffer_pdfd[977];

    auto g_z_yz_xyy_xx = buffer_pdfd[978];

    auto g_z_yz_xyy_xy = buffer_pdfd[979];

    auto g_z_yz_xyy_xz = buffer_pdfd[980];

    auto g_z_yz_xyy_yy = buffer_pdfd[981];

    auto g_z_yz_xyy_yz = buffer_pdfd[982];

    auto g_z_yz_xyy_zz = buffer_pdfd[983];

    auto g_z_yz_xyz_xx = buffer_pdfd[984];

    auto g_z_yz_xyz_xy = buffer_pdfd[985];

    auto g_z_yz_xyz_xz = buffer_pdfd[986];

    auto g_z_yz_xyz_yy = buffer_pdfd[987];

    auto g_z_yz_xyz_yz = buffer_pdfd[988];

    auto g_z_yz_xyz_zz = buffer_pdfd[989];

    auto g_z_yz_xzz_xx = buffer_pdfd[990];

    auto g_z_yz_xzz_xy = buffer_pdfd[991];

    auto g_z_yz_xzz_xz = buffer_pdfd[992];

    auto g_z_yz_xzz_yy = buffer_pdfd[993];

    auto g_z_yz_xzz_yz = buffer_pdfd[994];

    auto g_z_yz_xzz_zz = buffer_pdfd[995];

    auto g_z_yz_yyy_xx = buffer_pdfd[996];

    auto g_z_yz_yyy_xy = buffer_pdfd[997];

    auto g_z_yz_yyy_xz = buffer_pdfd[998];

    auto g_z_yz_yyy_yy = buffer_pdfd[999];

    auto g_z_yz_yyy_yz = buffer_pdfd[1000];

    auto g_z_yz_yyy_zz = buffer_pdfd[1001];

    auto g_z_yz_yyz_xx = buffer_pdfd[1002];

    auto g_z_yz_yyz_xy = buffer_pdfd[1003];

    auto g_z_yz_yyz_xz = buffer_pdfd[1004];

    auto g_z_yz_yyz_yy = buffer_pdfd[1005];

    auto g_z_yz_yyz_yz = buffer_pdfd[1006];

    auto g_z_yz_yyz_zz = buffer_pdfd[1007];

    auto g_z_yz_yzz_xx = buffer_pdfd[1008];

    auto g_z_yz_yzz_xy = buffer_pdfd[1009];

    auto g_z_yz_yzz_xz = buffer_pdfd[1010];

    auto g_z_yz_yzz_yy = buffer_pdfd[1011];

    auto g_z_yz_yzz_yz = buffer_pdfd[1012];

    auto g_z_yz_yzz_zz = buffer_pdfd[1013];

    auto g_z_yz_zzz_xx = buffer_pdfd[1014];

    auto g_z_yz_zzz_xy = buffer_pdfd[1015];

    auto g_z_yz_zzz_xz = buffer_pdfd[1016];

    auto g_z_yz_zzz_yy = buffer_pdfd[1017];

    auto g_z_yz_zzz_yz = buffer_pdfd[1018];

    auto g_z_yz_zzz_zz = buffer_pdfd[1019];

    auto g_z_zz_xxx_xx = buffer_pdfd[1020];

    auto g_z_zz_xxx_xy = buffer_pdfd[1021];

    auto g_z_zz_xxx_xz = buffer_pdfd[1022];

    auto g_z_zz_xxx_yy = buffer_pdfd[1023];

    auto g_z_zz_xxx_yz = buffer_pdfd[1024];

    auto g_z_zz_xxx_zz = buffer_pdfd[1025];

    auto g_z_zz_xxy_xx = buffer_pdfd[1026];

    auto g_z_zz_xxy_xy = buffer_pdfd[1027];

    auto g_z_zz_xxy_xz = buffer_pdfd[1028];

    auto g_z_zz_xxy_yy = buffer_pdfd[1029];

    auto g_z_zz_xxy_yz = buffer_pdfd[1030];

    auto g_z_zz_xxy_zz = buffer_pdfd[1031];

    auto g_z_zz_xxz_xx = buffer_pdfd[1032];

    auto g_z_zz_xxz_xy = buffer_pdfd[1033];

    auto g_z_zz_xxz_xz = buffer_pdfd[1034];

    auto g_z_zz_xxz_yy = buffer_pdfd[1035];

    auto g_z_zz_xxz_yz = buffer_pdfd[1036];

    auto g_z_zz_xxz_zz = buffer_pdfd[1037];

    auto g_z_zz_xyy_xx = buffer_pdfd[1038];

    auto g_z_zz_xyy_xy = buffer_pdfd[1039];

    auto g_z_zz_xyy_xz = buffer_pdfd[1040];

    auto g_z_zz_xyy_yy = buffer_pdfd[1041];

    auto g_z_zz_xyy_yz = buffer_pdfd[1042];

    auto g_z_zz_xyy_zz = buffer_pdfd[1043];

    auto g_z_zz_xyz_xx = buffer_pdfd[1044];

    auto g_z_zz_xyz_xy = buffer_pdfd[1045];

    auto g_z_zz_xyz_xz = buffer_pdfd[1046];

    auto g_z_zz_xyz_yy = buffer_pdfd[1047];

    auto g_z_zz_xyz_yz = buffer_pdfd[1048];

    auto g_z_zz_xyz_zz = buffer_pdfd[1049];

    auto g_z_zz_xzz_xx = buffer_pdfd[1050];

    auto g_z_zz_xzz_xy = buffer_pdfd[1051];

    auto g_z_zz_xzz_xz = buffer_pdfd[1052];

    auto g_z_zz_xzz_yy = buffer_pdfd[1053];

    auto g_z_zz_xzz_yz = buffer_pdfd[1054];

    auto g_z_zz_xzz_zz = buffer_pdfd[1055];

    auto g_z_zz_yyy_xx = buffer_pdfd[1056];

    auto g_z_zz_yyy_xy = buffer_pdfd[1057];

    auto g_z_zz_yyy_xz = buffer_pdfd[1058];

    auto g_z_zz_yyy_yy = buffer_pdfd[1059];

    auto g_z_zz_yyy_yz = buffer_pdfd[1060];

    auto g_z_zz_yyy_zz = buffer_pdfd[1061];

    auto g_z_zz_yyz_xx = buffer_pdfd[1062];

    auto g_z_zz_yyz_xy = buffer_pdfd[1063];

    auto g_z_zz_yyz_xz = buffer_pdfd[1064];

    auto g_z_zz_yyz_yy = buffer_pdfd[1065];

    auto g_z_zz_yyz_yz = buffer_pdfd[1066];

    auto g_z_zz_yyz_zz = buffer_pdfd[1067];

    auto g_z_zz_yzz_xx = buffer_pdfd[1068];

    auto g_z_zz_yzz_xy = buffer_pdfd[1069];

    auto g_z_zz_yzz_xz = buffer_pdfd[1070];

    auto g_z_zz_yzz_yy = buffer_pdfd[1071];

    auto g_z_zz_yzz_yz = buffer_pdfd[1072];

    auto g_z_zz_yzz_zz = buffer_pdfd[1073];

    auto g_z_zz_zzz_xx = buffer_pdfd[1074];

    auto g_z_zz_zzz_xy = buffer_pdfd[1075];

    auto g_z_zz_zzz_xz = buffer_pdfd[1076];

    auto g_z_zz_zzz_yy = buffer_pdfd[1077];

    auto g_z_zz_zzz_yz = buffer_pdfd[1078];

    auto g_z_zz_zzz_zz = buffer_pdfd[1079];

    /// Set up components of integrals buffer : buffer_1010_sddd

    auto g_x_0_x_0_0_xx_xx_xx = buffer_1010_sddd[0];

    auto g_x_0_x_0_0_xx_xx_xy = buffer_1010_sddd[1];

    auto g_x_0_x_0_0_xx_xx_xz = buffer_1010_sddd[2];

    auto g_x_0_x_0_0_xx_xx_yy = buffer_1010_sddd[3];

    auto g_x_0_x_0_0_xx_xx_yz = buffer_1010_sddd[4];

    auto g_x_0_x_0_0_xx_xx_zz = buffer_1010_sddd[5];

    auto g_x_0_x_0_0_xx_xy_xx = buffer_1010_sddd[6];

    auto g_x_0_x_0_0_xx_xy_xy = buffer_1010_sddd[7];

    auto g_x_0_x_0_0_xx_xy_xz = buffer_1010_sddd[8];

    auto g_x_0_x_0_0_xx_xy_yy = buffer_1010_sddd[9];

    auto g_x_0_x_0_0_xx_xy_yz = buffer_1010_sddd[10];

    auto g_x_0_x_0_0_xx_xy_zz = buffer_1010_sddd[11];

    auto g_x_0_x_0_0_xx_xz_xx = buffer_1010_sddd[12];

    auto g_x_0_x_0_0_xx_xz_xy = buffer_1010_sddd[13];

    auto g_x_0_x_0_0_xx_xz_xz = buffer_1010_sddd[14];

    auto g_x_0_x_0_0_xx_xz_yy = buffer_1010_sddd[15];

    auto g_x_0_x_0_0_xx_xz_yz = buffer_1010_sddd[16];

    auto g_x_0_x_0_0_xx_xz_zz = buffer_1010_sddd[17];

    auto g_x_0_x_0_0_xx_yy_xx = buffer_1010_sddd[18];

    auto g_x_0_x_0_0_xx_yy_xy = buffer_1010_sddd[19];

    auto g_x_0_x_0_0_xx_yy_xz = buffer_1010_sddd[20];

    auto g_x_0_x_0_0_xx_yy_yy = buffer_1010_sddd[21];

    auto g_x_0_x_0_0_xx_yy_yz = buffer_1010_sddd[22];

    auto g_x_0_x_0_0_xx_yy_zz = buffer_1010_sddd[23];

    auto g_x_0_x_0_0_xx_yz_xx = buffer_1010_sddd[24];

    auto g_x_0_x_0_0_xx_yz_xy = buffer_1010_sddd[25];

    auto g_x_0_x_0_0_xx_yz_xz = buffer_1010_sddd[26];

    auto g_x_0_x_0_0_xx_yz_yy = buffer_1010_sddd[27];

    auto g_x_0_x_0_0_xx_yz_yz = buffer_1010_sddd[28];

    auto g_x_0_x_0_0_xx_yz_zz = buffer_1010_sddd[29];

    auto g_x_0_x_0_0_xx_zz_xx = buffer_1010_sddd[30];

    auto g_x_0_x_0_0_xx_zz_xy = buffer_1010_sddd[31];

    auto g_x_0_x_0_0_xx_zz_xz = buffer_1010_sddd[32];

    auto g_x_0_x_0_0_xx_zz_yy = buffer_1010_sddd[33];

    auto g_x_0_x_0_0_xx_zz_yz = buffer_1010_sddd[34];

    auto g_x_0_x_0_0_xx_zz_zz = buffer_1010_sddd[35];

    auto g_x_0_x_0_0_xy_xx_xx = buffer_1010_sddd[36];

    auto g_x_0_x_0_0_xy_xx_xy = buffer_1010_sddd[37];

    auto g_x_0_x_0_0_xy_xx_xz = buffer_1010_sddd[38];

    auto g_x_0_x_0_0_xy_xx_yy = buffer_1010_sddd[39];

    auto g_x_0_x_0_0_xy_xx_yz = buffer_1010_sddd[40];

    auto g_x_0_x_0_0_xy_xx_zz = buffer_1010_sddd[41];

    auto g_x_0_x_0_0_xy_xy_xx = buffer_1010_sddd[42];

    auto g_x_0_x_0_0_xy_xy_xy = buffer_1010_sddd[43];

    auto g_x_0_x_0_0_xy_xy_xz = buffer_1010_sddd[44];

    auto g_x_0_x_0_0_xy_xy_yy = buffer_1010_sddd[45];

    auto g_x_0_x_0_0_xy_xy_yz = buffer_1010_sddd[46];

    auto g_x_0_x_0_0_xy_xy_zz = buffer_1010_sddd[47];

    auto g_x_0_x_0_0_xy_xz_xx = buffer_1010_sddd[48];

    auto g_x_0_x_0_0_xy_xz_xy = buffer_1010_sddd[49];

    auto g_x_0_x_0_0_xy_xz_xz = buffer_1010_sddd[50];

    auto g_x_0_x_0_0_xy_xz_yy = buffer_1010_sddd[51];

    auto g_x_0_x_0_0_xy_xz_yz = buffer_1010_sddd[52];

    auto g_x_0_x_0_0_xy_xz_zz = buffer_1010_sddd[53];

    auto g_x_0_x_0_0_xy_yy_xx = buffer_1010_sddd[54];

    auto g_x_0_x_0_0_xy_yy_xy = buffer_1010_sddd[55];

    auto g_x_0_x_0_0_xy_yy_xz = buffer_1010_sddd[56];

    auto g_x_0_x_0_0_xy_yy_yy = buffer_1010_sddd[57];

    auto g_x_0_x_0_0_xy_yy_yz = buffer_1010_sddd[58];

    auto g_x_0_x_0_0_xy_yy_zz = buffer_1010_sddd[59];

    auto g_x_0_x_0_0_xy_yz_xx = buffer_1010_sddd[60];

    auto g_x_0_x_0_0_xy_yz_xy = buffer_1010_sddd[61];

    auto g_x_0_x_0_0_xy_yz_xz = buffer_1010_sddd[62];

    auto g_x_0_x_0_0_xy_yz_yy = buffer_1010_sddd[63];

    auto g_x_0_x_0_0_xy_yz_yz = buffer_1010_sddd[64];

    auto g_x_0_x_0_0_xy_yz_zz = buffer_1010_sddd[65];

    auto g_x_0_x_0_0_xy_zz_xx = buffer_1010_sddd[66];

    auto g_x_0_x_0_0_xy_zz_xy = buffer_1010_sddd[67];

    auto g_x_0_x_0_0_xy_zz_xz = buffer_1010_sddd[68];

    auto g_x_0_x_0_0_xy_zz_yy = buffer_1010_sddd[69];

    auto g_x_0_x_0_0_xy_zz_yz = buffer_1010_sddd[70];

    auto g_x_0_x_0_0_xy_zz_zz = buffer_1010_sddd[71];

    auto g_x_0_x_0_0_xz_xx_xx = buffer_1010_sddd[72];

    auto g_x_0_x_0_0_xz_xx_xy = buffer_1010_sddd[73];

    auto g_x_0_x_0_0_xz_xx_xz = buffer_1010_sddd[74];

    auto g_x_0_x_0_0_xz_xx_yy = buffer_1010_sddd[75];

    auto g_x_0_x_0_0_xz_xx_yz = buffer_1010_sddd[76];

    auto g_x_0_x_0_0_xz_xx_zz = buffer_1010_sddd[77];

    auto g_x_0_x_0_0_xz_xy_xx = buffer_1010_sddd[78];

    auto g_x_0_x_0_0_xz_xy_xy = buffer_1010_sddd[79];

    auto g_x_0_x_0_0_xz_xy_xz = buffer_1010_sddd[80];

    auto g_x_0_x_0_0_xz_xy_yy = buffer_1010_sddd[81];

    auto g_x_0_x_0_0_xz_xy_yz = buffer_1010_sddd[82];

    auto g_x_0_x_0_0_xz_xy_zz = buffer_1010_sddd[83];

    auto g_x_0_x_0_0_xz_xz_xx = buffer_1010_sddd[84];

    auto g_x_0_x_0_0_xz_xz_xy = buffer_1010_sddd[85];

    auto g_x_0_x_0_0_xz_xz_xz = buffer_1010_sddd[86];

    auto g_x_0_x_0_0_xz_xz_yy = buffer_1010_sddd[87];

    auto g_x_0_x_0_0_xz_xz_yz = buffer_1010_sddd[88];

    auto g_x_0_x_0_0_xz_xz_zz = buffer_1010_sddd[89];

    auto g_x_0_x_0_0_xz_yy_xx = buffer_1010_sddd[90];

    auto g_x_0_x_0_0_xz_yy_xy = buffer_1010_sddd[91];

    auto g_x_0_x_0_0_xz_yy_xz = buffer_1010_sddd[92];

    auto g_x_0_x_0_0_xz_yy_yy = buffer_1010_sddd[93];

    auto g_x_0_x_0_0_xz_yy_yz = buffer_1010_sddd[94];

    auto g_x_0_x_0_0_xz_yy_zz = buffer_1010_sddd[95];

    auto g_x_0_x_0_0_xz_yz_xx = buffer_1010_sddd[96];

    auto g_x_0_x_0_0_xz_yz_xy = buffer_1010_sddd[97];

    auto g_x_0_x_0_0_xz_yz_xz = buffer_1010_sddd[98];

    auto g_x_0_x_0_0_xz_yz_yy = buffer_1010_sddd[99];

    auto g_x_0_x_0_0_xz_yz_yz = buffer_1010_sddd[100];

    auto g_x_0_x_0_0_xz_yz_zz = buffer_1010_sddd[101];

    auto g_x_0_x_0_0_xz_zz_xx = buffer_1010_sddd[102];

    auto g_x_0_x_0_0_xz_zz_xy = buffer_1010_sddd[103];

    auto g_x_0_x_0_0_xz_zz_xz = buffer_1010_sddd[104];

    auto g_x_0_x_0_0_xz_zz_yy = buffer_1010_sddd[105];

    auto g_x_0_x_0_0_xz_zz_yz = buffer_1010_sddd[106];

    auto g_x_0_x_0_0_xz_zz_zz = buffer_1010_sddd[107];

    auto g_x_0_x_0_0_yy_xx_xx = buffer_1010_sddd[108];

    auto g_x_0_x_0_0_yy_xx_xy = buffer_1010_sddd[109];

    auto g_x_0_x_0_0_yy_xx_xz = buffer_1010_sddd[110];

    auto g_x_0_x_0_0_yy_xx_yy = buffer_1010_sddd[111];

    auto g_x_0_x_0_0_yy_xx_yz = buffer_1010_sddd[112];

    auto g_x_0_x_0_0_yy_xx_zz = buffer_1010_sddd[113];

    auto g_x_0_x_0_0_yy_xy_xx = buffer_1010_sddd[114];

    auto g_x_0_x_0_0_yy_xy_xy = buffer_1010_sddd[115];

    auto g_x_0_x_0_0_yy_xy_xz = buffer_1010_sddd[116];

    auto g_x_0_x_0_0_yy_xy_yy = buffer_1010_sddd[117];

    auto g_x_0_x_0_0_yy_xy_yz = buffer_1010_sddd[118];

    auto g_x_0_x_0_0_yy_xy_zz = buffer_1010_sddd[119];

    auto g_x_0_x_0_0_yy_xz_xx = buffer_1010_sddd[120];

    auto g_x_0_x_0_0_yy_xz_xy = buffer_1010_sddd[121];

    auto g_x_0_x_0_0_yy_xz_xz = buffer_1010_sddd[122];

    auto g_x_0_x_0_0_yy_xz_yy = buffer_1010_sddd[123];

    auto g_x_0_x_0_0_yy_xz_yz = buffer_1010_sddd[124];

    auto g_x_0_x_0_0_yy_xz_zz = buffer_1010_sddd[125];

    auto g_x_0_x_0_0_yy_yy_xx = buffer_1010_sddd[126];

    auto g_x_0_x_0_0_yy_yy_xy = buffer_1010_sddd[127];

    auto g_x_0_x_0_0_yy_yy_xz = buffer_1010_sddd[128];

    auto g_x_0_x_0_0_yy_yy_yy = buffer_1010_sddd[129];

    auto g_x_0_x_0_0_yy_yy_yz = buffer_1010_sddd[130];

    auto g_x_0_x_0_0_yy_yy_zz = buffer_1010_sddd[131];

    auto g_x_0_x_0_0_yy_yz_xx = buffer_1010_sddd[132];

    auto g_x_0_x_0_0_yy_yz_xy = buffer_1010_sddd[133];

    auto g_x_0_x_0_0_yy_yz_xz = buffer_1010_sddd[134];

    auto g_x_0_x_0_0_yy_yz_yy = buffer_1010_sddd[135];

    auto g_x_0_x_0_0_yy_yz_yz = buffer_1010_sddd[136];

    auto g_x_0_x_0_0_yy_yz_zz = buffer_1010_sddd[137];

    auto g_x_0_x_0_0_yy_zz_xx = buffer_1010_sddd[138];

    auto g_x_0_x_0_0_yy_zz_xy = buffer_1010_sddd[139];

    auto g_x_0_x_0_0_yy_zz_xz = buffer_1010_sddd[140];

    auto g_x_0_x_0_0_yy_zz_yy = buffer_1010_sddd[141];

    auto g_x_0_x_0_0_yy_zz_yz = buffer_1010_sddd[142];

    auto g_x_0_x_0_0_yy_zz_zz = buffer_1010_sddd[143];

    auto g_x_0_x_0_0_yz_xx_xx = buffer_1010_sddd[144];

    auto g_x_0_x_0_0_yz_xx_xy = buffer_1010_sddd[145];

    auto g_x_0_x_0_0_yz_xx_xz = buffer_1010_sddd[146];

    auto g_x_0_x_0_0_yz_xx_yy = buffer_1010_sddd[147];

    auto g_x_0_x_0_0_yz_xx_yz = buffer_1010_sddd[148];

    auto g_x_0_x_0_0_yz_xx_zz = buffer_1010_sddd[149];

    auto g_x_0_x_0_0_yz_xy_xx = buffer_1010_sddd[150];

    auto g_x_0_x_0_0_yz_xy_xy = buffer_1010_sddd[151];

    auto g_x_0_x_0_0_yz_xy_xz = buffer_1010_sddd[152];

    auto g_x_0_x_0_0_yz_xy_yy = buffer_1010_sddd[153];

    auto g_x_0_x_0_0_yz_xy_yz = buffer_1010_sddd[154];

    auto g_x_0_x_0_0_yz_xy_zz = buffer_1010_sddd[155];

    auto g_x_0_x_0_0_yz_xz_xx = buffer_1010_sddd[156];

    auto g_x_0_x_0_0_yz_xz_xy = buffer_1010_sddd[157];

    auto g_x_0_x_0_0_yz_xz_xz = buffer_1010_sddd[158];

    auto g_x_0_x_0_0_yz_xz_yy = buffer_1010_sddd[159];

    auto g_x_0_x_0_0_yz_xz_yz = buffer_1010_sddd[160];

    auto g_x_0_x_0_0_yz_xz_zz = buffer_1010_sddd[161];

    auto g_x_0_x_0_0_yz_yy_xx = buffer_1010_sddd[162];

    auto g_x_0_x_0_0_yz_yy_xy = buffer_1010_sddd[163];

    auto g_x_0_x_0_0_yz_yy_xz = buffer_1010_sddd[164];

    auto g_x_0_x_0_0_yz_yy_yy = buffer_1010_sddd[165];

    auto g_x_0_x_0_0_yz_yy_yz = buffer_1010_sddd[166];

    auto g_x_0_x_0_0_yz_yy_zz = buffer_1010_sddd[167];

    auto g_x_0_x_0_0_yz_yz_xx = buffer_1010_sddd[168];

    auto g_x_0_x_0_0_yz_yz_xy = buffer_1010_sddd[169];

    auto g_x_0_x_0_0_yz_yz_xz = buffer_1010_sddd[170];

    auto g_x_0_x_0_0_yz_yz_yy = buffer_1010_sddd[171];

    auto g_x_0_x_0_0_yz_yz_yz = buffer_1010_sddd[172];

    auto g_x_0_x_0_0_yz_yz_zz = buffer_1010_sddd[173];

    auto g_x_0_x_0_0_yz_zz_xx = buffer_1010_sddd[174];

    auto g_x_0_x_0_0_yz_zz_xy = buffer_1010_sddd[175];

    auto g_x_0_x_0_0_yz_zz_xz = buffer_1010_sddd[176];

    auto g_x_0_x_0_0_yz_zz_yy = buffer_1010_sddd[177];

    auto g_x_0_x_0_0_yz_zz_yz = buffer_1010_sddd[178];

    auto g_x_0_x_0_0_yz_zz_zz = buffer_1010_sddd[179];

    auto g_x_0_x_0_0_zz_xx_xx = buffer_1010_sddd[180];

    auto g_x_0_x_0_0_zz_xx_xy = buffer_1010_sddd[181];

    auto g_x_0_x_0_0_zz_xx_xz = buffer_1010_sddd[182];

    auto g_x_0_x_0_0_zz_xx_yy = buffer_1010_sddd[183];

    auto g_x_0_x_0_0_zz_xx_yz = buffer_1010_sddd[184];

    auto g_x_0_x_0_0_zz_xx_zz = buffer_1010_sddd[185];

    auto g_x_0_x_0_0_zz_xy_xx = buffer_1010_sddd[186];

    auto g_x_0_x_0_0_zz_xy_xy = buffer_1010_sddd[187];

    auto g_x_0_x_0_0_zz_xy_xz = buffer_1010_sddd[188];

    auto g_x_0_x_0_0_zz_xy_yy = buffer_1010_sddd[189];

    auto g_x_0_x_0_0_zz_xy_yz = buffer_1010_sddd[190];

    auto g_x_0_x_0_0_zz_xy_zz = buffer_1010_sddd[191];

    auto g_x_0_x_0_0_zz_xz_xx = buffer_1010_sddd[192];

    auto g_x_0_x_0_0_zz_xz_xy = buffer_1010_sddd[193];

    auto g_x_0_x_0_0_zz_xz_xz = buffer_1010_sddd[194];

    auto g_x_0_x_0_0_zz_xz_yy = buffer_1010_sddd[195];

    auto g_x_0_x_0_0_zz_xz_yz = buffer_1010_sddd[196];

    auto g_x_0_x_0_0_zz_xz_zz = buffer_1010_sddd[197];

    auto g_x_0_x_0_0_zz_yy_xx = buffer_1010_sddd[198];

    auto g_x_0_x_0_0_zz_yy_xy = buffer_1010_sddd[199];

    auto g_x_0_x_0_0_zz_yy_xz = buffer_1010_sddd[200];

    auto g_x_0_x_0_0_zz_yy_yy = buffer_1010_sddd[201];

    auto g_x_0_x_0_0_zz_yy_yz = buffer_1010_sddd[202];

    auto g_x_0_x_0_0_zz_yy_zz = buffer_1010_sddd[203];

    auto g_x_0_x_0_0_zz_yz_xx = buffer_1010_sddd[204];

    auto g_x_0_x_0_0_zz_yz_xy = buffer_1010_sddd[205];

    auto g_x_0_x_0_0_zz_yz_xz = buffer_1010_sddd[206];

    auto g_x_0_x_0_0_zz_yz_yy = buffer_1010_sddd[207];

    auto g_x_0_x_0_0_zz_yz_yz = buffer_1010_sddd[208];

    auto g_x_0_x_0_0_zz_yz_zz = buffer_1010_sddd[209];

    auto g_x_0_x_0_0_zz_zz_xx = buffer_1010_sddd[210];

    auto g_x_0_x_0_0_zz_zz_xy = buffer_1010_sddd[211];

    auto g_x_0_x_0_0_zz_zz_xz = buffer_1010_sddd[212];

    auto g_x_0_x_0_0_zz_zz_yy = buffer_1010_sddd[213];

    auto g_x_0_x_0_0_zz_zz_yz = buffer_1010_sddd[214];

    auto g_x_0_x_0_0_zz_zz_zz = buffer_1010_sddd[215];

    auto g_x_0_y_0_0_xx_xx_xx = buffer_1010_sddd[216];

    auto g_x_0_y_0_0_xx_xx_xy = buffer_1010_sddd[217];

    auto g_x_0_y_0_0_xx_xx_xz = buffer_1010_sddd[218];

    auto g_x_0_y_0_0_xx_xx_yy = buffer_1010_sddd[219];

    auto g_x_0_y_0_0_xx_xx_yz = buffer_1010_sddd[220];

    auto g_x_0_y_0_0_xx_xx_zz = buffer_1010_sddd[221];

    auto g_x_0_y_0_0_xx_xy_xx = buffer_1010_sddd[222];

    auto g_x_0_y_0_0_xx_xy_xy = buffer_1010_sddd[223];

    auto g_x_0_y_0_0_xx_xy_xz = buffer_1010_sddd[224];

    auto g_x_0_y_0_0_xx_xy_yy = buffer_1010_sddd[225];

    auto g_x_0_y_0_0_xx_xy_yz = buffer_1010_sddd[226];

    auto g_x_0_y_0_0_xx_xy_zz = buffer_1010_sddd[227];

    auto g_x_0_y_0_0_xx_xz_xx = buffer_1010_sddd[228];

    auto g_x_0_y_0_0_xx_xz_xy = buffer_1010_sddd[229];

    auto g_x_0_y_0_0_xx_xz_xz = buffer_1010_sddd[230];

    auto g_x_0_y_0_0_xx_xz_yy = buffer_1010_sddd[231];

    auto g_x_0_y_0_0_xx_xz_yz = buffer_1010_sddd[232];

    auto g_x_0_y_0_0_xx_xz_zz = buffer_1010_sddd[233];

    auto g_x_0_y_0_0_xx_yy_xx = buffer_1010_sddd[234];

    auto g_x_0_y_0_0_xx_yy_xy = buffer_1010_sddd[235];

    auto g_x_0_y_0_0_xx_yy_xz = buffer_1010_sddd[236];

    auto g_x_0_y_0_0_xx_yy_yy = buffer_1010_sddd[237];

    auto g_x_0_y_0_0_xx_yy_yz = buffer_1010_sddd[238];

    auto g_x_0_y_0_0_xx_yy_zz = buffer_1010_sddd[239];

    auto g_x_0_y_0_0_xx_yz_xx = buffer_1010_sddd[240];

    auto g_x_0_y_0_0_xx_yz_xy = buffer_1010_sddd[241];

    auto g_x_0_y_0_0_xx_yz_xz = buffer_1010_sddd[242];

    auto g_x_0_y_0_0_xx_yz_yy = buffer_1010_sddd[243];

    auto g_x_0_y_0_0_xx_yz_yz = buffer_1010_sddd[244];

    auto g_x_0_y_0_0_xx_yz_zz = buffer_1010_sddd[245];

    auto g_x_0_y_0_0_xx_zz_xx = buffer_1010_sddd[246];

    auto g_x_0_y_0_0_xx_zz_xy = buffer_1010_sddd[247];

    auto g_x_0_y_0_0_xx_zz_xz = buffer_1010_sddd[248];

    auto g_x_0_y_0_0_xx_zz_yy = buffer_1010_sddd[249];

    auto g_x_0_y_0_0_xx_zz_yz = buffer_1010_sddd[250];

    auto g_x_0_y_0_0_xx_zz_zz = buffer_1010_sddd[251];

    auto g_x_0_y_0_0_xy_xx_xx = buffer_1010_sddd[252];

    auto g_x_0_y_0_0_xy_xx_xy = buffer_1010_sddd[253];

    auto g_x_0_y_0_0_xy_xx_xz = buffer_1010_sddd[254];

    auto g_x_0_y_0_0_xy_xx_yy = buffer_1010_sddd[255];

    auto g_x_0_y_0_0_xy_xx_yz = buffer_1010_sddd[256];

    auto g_x_0_y_0_0_xy_xx_zz = buffer_1010_sddd[257];

    auto g_x_0_y_0_0_xy_xy_xx = buffer_1010_sddd[258];

    auto g_x_0_y_0_0_xy_xy_xy = buffer_1010_sddd[259];

    auto g_x_0_y_0_0_xy_xy_xz = buffer_1010_sddd[260];

    auto g_x_0_y_0_0_xy_xy_yy = buffer_1010_sddd[261];

    auto g_x_0_y_0_0_xy_xy_yz = buffer_1010_sddd[262];

    auto g_x_0_y_0_0_xy_xy_zz = buffer_1010_sddd[263];

    auto g_x_0_y_0_0_xy_xz_xx = buffer_1010_sddd[264];

    auto g_x_0_y_0_0_xy_xz_xy = buffer_1010_sddd[265];

    auto g_x_0_y_0_0_xy_xz_xz = buffer_1010_sddd[266];

    auto g_x_0_y_0_0_xy_xz_yy = buffer_1010_sddd[267];

    auto g_x_0_y_0_0_xy_xz_yz = buffer_1010_sddd[268];

    auto g_x_0_y_0_0_xy_xz_zz = buffer_1010_sddd[269];

    auto g_x_0_y_0_0_xy_yy_xx = buffer_1010_sddd[270];

    auto g_x_0_y_0_0_xy_yy_xy = buffer_1010_sddd[271];

    auto g_x_0_y_0_0_xy_yy_xz = buffer_1010_sddd[272];

    auto g_x_0_y_0_0_xy_yy_yy = buffer_1010_sddd[273];

    auto g_x_0_y_0_0_xy_yy_yz = buffer_1010_sddd[274];

    auto g_x_0_y_0_0_xy_yy_zz = buffer_1010_sddd[275];

    auto g_x_0_y_0_0_xy_yz_xx = buffer_1010_sddd[276];

    auto g_x_0_y_0_0_xy_yz_xy = buffer_1010_sddd[277];

    auto g_x_0_y_0_0_xy_yz_xz = buffer_1010_sddd[278];

    auto g_x_0_y_0_0_xy_yz_yy = buffer_1010_sddd[279];

    auto g_x_0_y_0_0_xy_yz_yz = buffer_1010_sddd[280];

    auto g_x_0_y_0_0_xy_yz_zz = buffer_1010_sddd[281];

    auto g_x_0_y_0_0_xy_zz_xx = buffer_1010_sddd[282];

    auto g_x_0_y_0_0_xy_zz_xy = buffer_1010_sddd[283];

    auto g_x_0_y_0_0_xy_zz_xz = buffer_1010_sddd[284];

    auto g_x_0_y_0_0_xy_zz_yy = buffer_1010_sddd[285];

    auto g_x_0_y_0_0_xy_zz_yz = buffer_1010_sddd[286];

    auto g_x_0_y_0_0_xy_zz_zz = buffer_1010_sddd[287];

    auto g_x_0_y_0_0_xz_xx_xx = buffer_1010_sddd[288];

    auto g_x_0_y_0_0_xz_xx_xy = buffer_1010_sddd[289];

    auto g_x_0_y_0_0_xz_xx_xz = buffer_1010_sddd[290];

    auto g_x_0_y_0_0_xz_xx_yy = buffer_1010_sddd[291];

    auto g_x_0_y_0_0_xz_xx_yz = buffer_1010_sddd[292];

    auto g_x_0_y_0_0_xz_xx_zz = buffer_1010_sddd[293];

    auto g_x_0_y_0_0_xz_xy_xx = buffer_1010_sddd[294];

    auto g_x_0_y_0_0_xz_xy_xy = buffer_1010_sddd[295];

    auto g_x_0_y_0_0_xz_xy_xz = buffer_1010_sddd[296];

    auto g_x_0_y_0_0_xz_xy_yy = buffer_1010_sddd[297];

    auto g_x_0_y_0_0_xz_xy_yz = buffer_1010_sddd[298];

    auto g_x_0_y_0_0_xz_xy_zz = buffer_1010_sddd[299];

    auto g_x_0_y_0_0_xz_xz_xx = buffer_1010_sddd[300];

    auto g_x_0_y_0_0_xz_xz_xy = buffer_1010_sddd[301];

    auto g_x_0_y_0_0_xz_xz_xz = buffer_1010_sddd[302];

    auto g_x_0_y_0_0_xz_xz_yy = buffer_1010_sddd[303];

    auto g_x_0_y_0_0_xz_xz_yz = buffer_1010_sddd[304];

    auto g_x_0_y_0_0_xz_xz_zz = buffer_1010_sddd[305];

    auto g_x_0_y_0_0_xz_yy_xx = buffer_1010_sddd[306];

    auto g_x_0_y_0_0_xz_yy_xy = buffer_1010_sddd[307];

    auto g_x_0_y_0_0_xz_yy_xz = buffer_1010_sddd[308];

    auto g_x_0_y_0_0_xz_yy_yy = buffer_1010_sddd[309];

    auto g_x_0_y_0_0_xz_yy_yz = buffer_1010_sddd[310];

    auto g_x_0_y_0_0_xz_yy_zz = buffer_1010_sddd[311];

    auto g_x_0_y_0_0_xz_yz_xx = buffer_1010_sddd[312];

    auto g_x_0_y_0_0_xz_yz_xy = buffer_1010_sddd[313];

    auto g_x_0_y_0_0_xz_yz_xz = buffer_1010_sddd[314];

    auto g_x_0_y_0_0_xz_yz_yy = buffer_1010_sddd[315];

    auto g_x_0_y_0_0_xz_yz_yz = buffer_1010_sddd[316];

    auto g_x_0_y_0_0_xz_yz_zz = buffer_1010_sddd[317];

    auto g_x_0_y_0_0_xz_zz_xx = buffer_1010_sddd[318];

    auto g_x_0_y_0_0_xz_zz_xy = buffer_1010_sddd[319];

    auto g_x_0_y_0_0_xz_zz_xz = buffer_1010_sddd[320];

    auto g_x_0_y_0_0_xz_zz_yy = buffer_1010_sddd[321];

    auto g_x_0_y_0_0_xz_zz_yz = buffer_1010_sddd[322];

    auto g_x_0_y_0_0_xz_zz_zz = buffer_1010_sddd[323];

    auto g_x_0_y_0_0_yy_xx_xx = buffer_1010_sddd[324];

    auto g_x_0_y_0_0_yy_xx_xy = buffer_1010_sddd[325];

    auto g_x_0_y_0_0_yy_xx_xz = buffer_1010_sddd[326];

    auto g_x_0_y_0_0_yy_xx_yy = buffer_1010_sddd[327];

    auto g_x_0_y_0_0_yy_xx_yz = buffer_1010_sddd[328];

    auto g_x_0_y_0_0_yy_xx_zz = buffer_1010_sddd[329];

    auto g_x_0_y_0_0_yy_xy_xx = buffer_1010_sddd[330];

    auto g_x_0_y_0_0_yy_xy_xy = buffer_1010_sddd[331];

    auto g_x_0_y_0_0_yy_xy_xz = buffer_1010_sddd[332];

    auto g_x_0_y_0_0_yy_xy_yy = buffer_1010_sddd[333];

    auto g_x_0_y_0_0_yy_xy_yz = buffer_1010_sddd[334];

    auto g_x_0_y_0_0_yy_xy_zz = buffer_1010_sddd[335];

    auto g_x_0_y_0_0_yy_xz_xx = buffer_1010_sddd[336];

    auto g_x_0_y_0_0_yy_xz_xy = buffer_1010_sddd[337];

    auto g_x_0_y_0_0_yy_xz_xz = buffer_1010_sddd[338];

    auto g_x_0_y_0_0_yy_xz_yy = buffer_1010_sddd[339];

    auto g_x_0_y_0_0_yy_xz_yz = buffer_1010_sddd[340];

    auto g_x_0_y_0_0_yy_xz_zz = buffer_1010_sddd[341];

    auto g_x_0_y_0_0_yy_yy_xx = buffer_1010_sddd[342];

    auto g_x_0_y_0_0_yy_yy_xy = buffer_1010_sddd[343];

    auto g_x_0_y_0_0_yy_yy_xz = buffer_1010_sddd[344];

    auto g_x_0_y_0_0_yy_yy_yy = buffer_1010_sddd[345];

    auto g_x_0_y_0_0_yy_yy_yz = buffer_1010_sddd[346];

    auto g_x_0_y_0_0_yy_yy_zz = buffer_1010_sddd[347];

    auto g_x_0_y_0_0_yy_yz_xx = buffer_1010_sddd[348];

    auto g_x_0_y_0_0_yy_yz_xy = buffer_1010_sddd[349];

    auto g_x_0_y_0_0_yy_yz_xz = buffer_1010_sddd[350];

    auto g_x_0_y_0_0_yy_yz_yy = buffer_1010_sddd[351];

    auto g_x_0_y_0_0_yy_yz_yz = buffer_1010_sddd[352];

    auto g_x_0_y_0_0_yy_yz_zz = buffer_1010_sddd[353];

    auto g_x_0_y_0_0_yy_zz_xx = buffer_1010_sddd[354];

    auto g_x_0_y_0_0_yy_zz_xy = buffer_1010_sddd[355];

    auto g_x_0_y_0_0_yy_zz_xz = buffer_1010_sddd[356];

    auto g_x_0_y_0_0_yy_zz_yy = buffer_1010_sddd[357];

    auto g_x_0_y_0_0_yy_zz_yz = buffer_1010_sddd[358];

    auto g_x_0_y_0_0_yy_zz_zz = buffer_1010_sddd[359];

    auto g_x_0_y_0_0_yz_xx_xx = buffer_1010_sddd[360];

    auto g_x_0_y_0_0_yz_xx_xy = buffer_1010_sddd[361];

    auto g_x_0_y_0_0_yz_xx_xz = buffer_1010_sddd[362];

    auto g_x_0_y_0_0_yz_xx_yy = buffer_1010_sddd[363];

    auto g_x_0_y_0_0_yz_xx_yz = buffer_1010_sddd[364];

    auto g_x_0_y_0_0_yz_xx_zz = buffer_1010_sddd[365];

    auto g_x_0_y_0_0_yz_xy_xx = buffer_1010_sddd[366];

    auto g_x_0_y_0_0_yz_xy_xy = buffer_1010_sddd[367];

    auto g_x_0_y_0_0_yz_xy_xz = buffer_1010_sddd[368];

    auto g_x_0_y_0_0_yz_xy_yy = buffer_1010_sddd[369];

    auto g_x_0_y_0_0_yz_xy_yz = buffer_1010_sddd[370];

    auto g_x_0_y_0_0_yz_xy_zz = buffer_1010_sddd[371];

    auto g_x_0_y_0_0_yz_xz_xx = buffer_1010_sddd[372];

    auto g_x_0_y_0_0_yz_xz_xy = buffer_1010_sddd[373];

    auto g_x_0_y_0_0_yz_xz_xz = buffer_1010_sddd[374];

    auto g_x_0_y_0_0_yz_xz_yy = buffer_1010_sddd[375];

    auto g_x_0_y_0_0_yz_xz_yz = buffer_1010_sddd[376];

    auto g_x_0_y_0_0_yz_xz_zz = buffer_1010_sddd[377];

    auto g_x_0_y_0_0_yz_yy_xx = buffer_1010_sddd[378];

    auto g_x_0_y_0_0_yz_yy_xy = buffer_1010_sddd[379];

    auto g_x_0_y_0_0_yz_yy_xz = buffer_1010_sddd[380];

    auto g_x_0_y_0_0_yz_yy_yy = buffer_1010_sddd[381];

    auto g_x_0_y_0_0_yz_yy_yz = buffer_1010_sddd[382];

    auto g_x_0_y_0_0_yz_yy_zz = buffer_1010_sddd[383];

    auto g_x_0_y_0_0_yz_yz_xx = buffer_1010_sddd[384];

    auto g_x_0_y_0_0_yz_yz_xy = buffer_1010_sddd[385];

    auto g_x_0_y_0_0_yz_yz_xz = buffer_1010_sddd[386];

    auto g_x_0_y_0_0_yz_yz_yy = buffer_1010_sddd[387];

    auto g_x_0_y_0_0_yz_yz_yz = buffer_1010_sddd[388];

    auto g_x_0_y_0_0_yz_yz_zz = buffer_1010_sddd[389];

    auto g_x_0_y_0_0_yz_zz_xx = buffer_1010_sddd[390];

    auto g_x_0_y_0_0_yz_zz_xy = buffer_1010_sddd[391];

    auto g_x_0_y_0_0_yz_zz_xz = buffer_1010_sddd[392];

    auto g_x_0_y_0_0_yz_zz_yy = buffer_1010_sddd[393];

    auto g_x_0_y_0_0_yz_zz_yz = buffer_1010_sddd[394];

    auto g_x_0_y_0_0_yz_zz_zz = buffer_1010_sddd[395];

    auto g_x_0_y_0_0_zz_xx_xx = buffer_1010_sddd[396];

    auto g_x_0_y_0_0_zz_xx_xy = buffer_1010_sddd[397];

    auto g_x_0_y_0_0_zz_xx_xz = buffer_1010_sddd[398];

    auto g_x_0_y_0_0_zz_xx_yy = buffer_1010_sddd[399];

    auto g_x_0_y_0_0_zz_xx_yz = buffer_1010_sddd[400];

    auto g_x_0_y_0_0_zz_xx_zz = buffer_1010_sddd[401];

    auto g_x_0_y_0_0_zz_xy_xx = buffer_1010_sddd[402];

    auto g_x_0_y_0_0_zz_xy_xy = buffer_1010_sddd[403];

    auto g_x_0_y_0_0_zz_xy_xz = buffer_1010_sddd[404];

    auto g_x_0_y_0_0_zz_xy_yy = buffer_1010_sddd[405];

    auto g_x_0_y_0_0_zz_xy_yz = buffer_1010_sddd[406];

    auto g_x_0_y_0_0_zz_xy_zz = buffer_1010_sddd[407];

    auto g_x_0_y_0_0_zz_xz_xx = buffer_1010_sddd[408];

    auto g_x_0_y_0_0_zz_xz_xy = buffer_1010_sddd[409];

    auto g_x_0_y_0_0_zz_xz_xz = buffer_1010_sddd[410];

    auto g_x_0_y_0_0_zz_xz_yy = buffer_1010_sddd[411];

    auto g_x_0_y_0_0_zz_xz_yz = buffer_1010_sddd[412];

    auto g_x_0_y_0_0_zz_xz_zz = buffer_1010_sddd[413];

    auto g_x_0_y_0_0_zz_yy_xx = buffer_1010_sddd[414];

    auto g_x_0_y_0_0_zz_yy_xy = buffer_1010_sddd[415];

    auto g_x_0_y_0_0_zz_yy_xz = buffer_1010_sddd[416];

    auto g_x_0_y_0_0_zz_yy_yy = buffer_1010_sddd[417];

    auto g_x_0_y_0_0_zz_yy_yz = buffer_1010_sddd[418];

    auto g_x_0_y_0_0_zz_yy_zz = buffer_1010_sddd[419];

    auto g_x_0_y_0_0_zz_yz_xx = buffer_1010_sddd[420];

    auto g_x_0_y_0_0_zz_yz_xy = buffer_1010_sddd[421];

    auto g_x_0_y_0_0_zz_yz_xz = buffer_1010_sddd[422];

    auto g_x_0_y_0_0_zz_yz_yy = buffer_1010_sddd[423];

    auto g_x_0_y_0_0_zz_yz_yz = buffer_1010_sddd[424];

    auto g_x_0_y_0_0_zz_yz_zz = buffer_1010_sddd[425];

    auto g_x_0_y_0_0_zz_zz_xx = buffer_1010_sddd[426];

    auto g_x_0_y_0_0_zz_zz_xy = buffer_1010_sddd[427];

    auto g_x_0_y_0_0_zz_zz_xz = buffer_1010_sddd[428];

    auto g_x_0_y_0_0_zz_zz_yy = buffer_1010_sddd[429];

    auto g_x_0_y_0_0_zz_zz_yz = buffer_1010_sddd[430];

    auto g_x_0_y_0_0_zz_zz_zz = buffer_1010_sddd[431];

    auto g_x_0_z_0_0_xx_xx_xx = buffer_1010_sddd[432];

    auto g_x_0_z_0_0_xx_xx_xy = buffer_1010_sddd[433];

    auto g_x_0_z_0_0_xx_xx_xz = buffer_1010_sddd[434];

    auto g_x_0_z_0_0_xx_xx_yy = buffer_1010_sddd[435];

    auto g_x_0_z_0_0_xx_xx_yz = buffer_1010_sddd[436];

    auto g_x_0_z_0_0_xx_xx_zz = buffer_1010_sddd[437];

    auto g_x_0_z_0_0_xx_xy_xx = buffer_1010_sddd[438];

    auto g_x_0_z_0_0_xx_xy_xy = buffer_1010_sddd[439];

    auto g_x_0_z_0_0_xx_xy_xz = buffer_1010_sddd[440];

    auto g_x_0_z_0_0_xx_xy_yy = buffer_1010_sddd[441];

    auto g_x_0_z_0_0_xx_xy_yz = buffer_1010_sddd[442];

    auto g_x_0_z_0_0_xx_xy_zz = buffer_1010_sddd[443];

    auto g_x_0_z_0_0_xx_xz_xx = buffer_1010_sddd[444];

    auto g_x_0_z_0_0_xx_xz_xy = buffer_1010_sddd[445];

    auto g_x_0_z_0_0_xx_xz_xz = buffer_1010_sddd[446];

    auto g_x_0_z_0_0_xx_xz_yy = buffer_1010_sddd[447];

    auto g_x_0_z_0_0_xx_xz_yz = buffer_1010_sddd[448];

    auto g_x_0_z_0_0_xx_xz_zz = buffer_1010_sddd[449];

    auto g_x_0_z_0_0_xx_yy_xx = buffer_1010_sddd[450];

    auto g_x_0_z_0_0_xx_yy_xy = buffer_1010_sddd[451];

    auto g_x_0_z_0_0_xx_yy_xz = buffer_1010_sddd[452];

    auto g_x_0_z_0_0_xx_yy_yy = buffer_1010_sddd[453];

    auto g_x_0_z_0_0_xx_yy_yz = buffer_1010_sddd[454];

    auto g_x_0_z_0_0_xx_yy_zz = buffer_1010_sddd[455];

    auto g_x_0_z_0_0_xx_yz_xx = buffer_1010_sddd[456];

    auto g_x_0_z_0_0_xx_yz_xy = buffer_1010_sddd[457];

    auto g_x_0_z_0_0_xx_yz_xz = buffer_1010_sddd[458];

    auto g_x_0_z_0_0_xx_yz_yy = buffer_1010_sddd[459];

    auto g_x_0_z_0_0_xx_yz_yz = buffer_1010_sddd[460];

    auto g_x_0_z_0_0_xx_yz_zz = buffer_1010_sddd[461];

    auto g_x_0_z_0_0_xx_zz_xx = buffer_1010_sddd[462];

    auto g_x_0_z_0_0_xx_zz_xy = buffer_1010_sddd[463];

    auto g_x_0_z_0_0_xx_zz_xz = buffer_1010_sddd[464];

    auto g_x_0_z_0_0_xx_zz_yy = buffer_1010_sddd[465];

    auto g_x_0_z_0_0_xx_zz_yz = buffer_1010_sddd[466];

    auto g_x_0_z_0_0_xx_zz_zz = buffer_1010_sddd[467];

    auto g_x_0_z_0_0_xy_xx_xx = buffer_1010_sddd[468];

    auto g_x_0_z_0_0_xy_xx_xy = buffer_1010_sddd[469];

    auto g_x_0_z_0_0_xy_xx_xz = buffer_1010_sddd[470];

    auto g_x_0_z_0_0_xy_xx_yy = buffer_1010_sddd[471];

    auto g_x_0_z_0_0_xy_xx_yz = buffer_1010_sddd[472];

    auto g_x_0_z_0_0_xy_xx_zz = buffer_1010_sddd[473];

    auto g_x_0_z_0_0_xy_xy_xx = buffer_1010_sddd[474];

    auto g_x_0_z_0_0_xy_xy_xy = buffer_1010_sddd[475];

    auto g_x_0_z_0_0_xy_xy_xz = buffer_1010_sddd[476];

    auto g_x_0_z_0_0_xy_xy_yy = buffer_1010_sddd[477];

    auto g_x_0_z_0_0_xy_xy_yz = buffer_1010_sddd[478];

    auto g_x_0_z_0_0_xy_xy_zz = buffer_1010_sddd[479];

    auto g_x_0_z_0_0_xy_xz_xx = buffer_1010_sddd[480];

    auto g_x_0_z_0_0_xy_xz_xy = buffer_1010_sddd[481];

    auto g_x_0_z_0_0_xy_xz_xz = buffer_1010_sddd[482];

    auto g_x_0_z_0_0_xy_xz_yy = buffer_1010_sddd[483];

    auto g_x_0_z_0_0_xy_xz_yz = buffer_1010_sddd[484];

    auto g_x_0_z_0_0_xy_xz_zz = buffer_1010_sddd[485];

    auto g_x_0_z_0_0_xy_yy_xx = buffer_1010_sddd[486];

    auto g_x_0_z_0_0_xy_yy_xy = buffer_1010_sddd[487];

    auto g_x_0_z_0_0_xy_yy_xz = buffer_1010_sddd[488];

    auto g_x_0_z_0_0_xy_yy_yy = buffer_1010_sddd[489];

    auto g_x_0_z_0_0_xy_yy_yz = buffer_1010_sddd[490];

    auto g_x_0_z_0_0_xy_yy_zz = buffer_1010_sddd[491];

    auto g_x_0_z_0_0_xy_yz_xx = buffer_1010_sddd[492];

    auto g_x_0_z_0_0_xy_yz_xy = buffer_1010_sddd[493];

    auto g_x_0_z_0_0_xy_yz_xz = buffer_1010_sddd[494];

    auto g_x_0_z_0_0_xy_yz_yy = buffer_1010_sddd[495];

    auto g_x_0_z_0_0_xy_yz_yz = buffer_1010_sddd[496];

    auto g_x_0_z_0_0_xy_yz_zz = buffer_1010_sddd[497];

    auto g_x_0_z_0_0_xy_zz_xx = buffer_1010_sddd[498];

    auto g_x_0_z_0_0_xy_zz_xy = buffer_1010_sddd[499];

    auto g_x_0_z_0_0_xy_zz_xz = buffer_1010_sddd[500];

    auto g_x_0_z_0_0_xy_zz_yy = buffer_1010_sddd[501];

    auto g_x_0_z_0_0_xy_zz_yz = buffer_1010_sddd[502];

    auto g_x_0_z_0_0_xy_zz_zz = buffer_1010_sddd[503];

    auto g_x_0_z_0_0_xz_xx_xx = buffer_1010_sddd[504];

    auto g_x_0_z_0_0_xz_xx_xy = buffer_1010_sddd[505];

    auto g_x_0_z_0_0_xz_xx_xz = buffer_1010_sddd[506];

    auto g_x_0_z_0_0_xz_xx_yy = buffer_1010_sddd[507];

    auto g_x_0_z_0_0_xz_xx_yz = buffer_1010_sddd[508];

    auto g_x_0_z_0_0_xz_xx_zz = buffer_1010_sddd[509];

    auto g_x_0_z_0_0_xz_xy_xx = buffer_1010_sddd[510];

    auto g_x_0_z_0_0_xz_xy_xy = buffer_1010_sddd[511];

    auto g_x_0_z_0_0_xz_xy_xz = buffer_1010_sddd[512];

    auto g_x_0_z_0_0_xz_xy_yy = buffer_1010_sddd[513];

    auto g_x_0_z_0_0_xz_xy_yz = buffer_1010_sddd[514];

    auto g_x_0_z_0_0_xz_xy_zz = buffer_1010_sddd[515];

    auto g_x_0_z_0_0_xz_xz_xx = buffer_1010_sddd[516];

    auto g_x_0_z_0_0_xz_xz_xy = buffer_1010_sddd[517];

    auto g_x_0_z_0_0_xz_xz_xz = buffer_1010_sddd[518];

    auto g_x_0_z_0_0_xz_xz_yy = buffer_1010_sddd[519];

    auto g_x_0_z_0_0_xz_xz_yz = buffer_1010_sddd[520];

    auto g_x_0_z_0_0_xz_xz_zz = buffer_1010_sddd[521];

    auto g_x_0_z_0_0_xz_yy_xx = buffer_1010_sddd[522];

    auto g_x_0_z_0_0_xz_yy_xy = buffer_1010_sddd[523];

    auto g_x_0_z_0_0_xz_yy_xz = buffer_1010_sddd[524];

    auto g_x_0_z_0_0_xz_yy_yy = buffer_1010_sddd[525];

    auto g_x_0_z_0_0_xz_yy_yz = buffer_1010_sddd[526];

    auto g_x_0_z_0_0_xz_yy_zz = buffer_1010_sddd[527];

    auto g_x_0_z_0_0_xz_yz_xx = buffer_1010_sddd[528];

    auto g_x_0_z_0_0_xz_yz_xy = buffer_1010_sddd[529];

    auto g_x_0_z_0_0_xz_yz_xz = buffer_1010_sddd[530];

    auto g_x_0_z_0_0_xz_yz_yy = buffer_1010_sddd[531];

    auto g_x_0_z_0_0_xz_yz_yz = buffer_1010_sddd[532];

    auto g_x_0_z_0_0_xz_yz_zz = buffer_1010_sddd[533];

    auto g_x_0_z_0_0_xz_zz_xx = buffer_1010_sddd[534];

    auto g_x_0_z_0_0_xz_zz_xy = buffer_1010_sddd[535];

    auto g_x_0_z_0_0_xz_zz_xz = buffer_1010_sddd[536];

    auto g_x_0_z_0_0_xz_zz_yy = buffer_1010_sddd[537];

    auto g_x_0_z_0_0_xz_zz_yz = buffer_1010_sddd[538];

    auto g_x_0_z_0_0_xz_zz_zz = buffer_1010_sddd[539];

    auto g_x_0_z_0_0_yy_xx_xx = buffer_1010_sddd[540];

    auto g_x_0_z_0_0_yy_xx_xy = buffer_1010_sddd[541];

    auto g_x_0_z_0_0_yy_xx_xz = buffer_1010_sddd[542];

    auto g_x_0_z_0_0_yy_xx_yy = buffer_1010_sddd[543];

    auto g_x_0_z_0_0_yy_xx_yz = buffer_1010_sddd[544];

    auto g_x_0_z_0_0_yy_xx_zz = buffer_1010_sddd[545];

    auto g_x_0_z_0_0_yy_xy_xx = buffer_1010_sddd[546];

    auto g_x_0_z_0_0_yy_xy_xy = buffer_1010_sddd[547];

    auto g_x_0_z_0_0_yy_xy_xz = buffer_1010_sddd[548];

    auto g_x_0_z_0_0_yy_xy_yy = buffer_1010_sddd[549];

    auto g_x_0_z_0_0_yy_xy_yz = buffer_1010_sddd[550];

    auto g_x_0_z_0_0_yy_xy_zz = buffer_1010_sddd[551];

    auto g_x_0_z_0_0_yy_xz_xx = buffer_1010_sddd[552];

    auto g_x_0_z_0_0_yy_xz_xy = buffer_1010_sddd[553];

    auto g_x_0_z_0_0_yy_xz_xz = buffer_1010_sddd[554];

    auto g_x_0_z_0_0_yy_xz_yy = buffer_1010_sddd[555];

    auto g_x_0_z_0_0_yy_xz_yz = buffer_1010_sddd[556];

    auto g_x_0_z_0_0_yy_xz_zz = buffer_1010_sddd[557];

    auto g_x_0_z_0_0_yy_yy_xx = buffer_1010_sddd[558];

    auto g_x_0_z_0_0_yy_yy_xy = buffer_1010_sddd[559];

    auto g_x_0_z_0_0_yy_yy_xz = buffer_1010_sddd[560];

    auto g_x_0_z_0_0_yy_yy_yy = buffer_1010_sddd[561];

    auto g_x_0_z_0_0_yy_yy_yz = buffer_1010_sddd[562];

    auto g_x_0_z_0_0_yy_yy_zz = buffer_1010_sddd[563];

    auto g_x_0_z_0_0_yy_yz_xx = buffer_1010_sddd[564];

    auto g_x_0_z_0_0_yy_yz_xy = buffer_1010_sddd[565];

    auto g_x_0_z_0_0_yy_yz_xz = buffer_1010_sddd[566];

    auto g_x_0_z_0_0_yy_yz_yy = buffer_1010_sddd[567];

    auto g_x_0_z_0_0_yy_yz_yz = buffer_1010_sddd[568];

    auto g_x_0_z_0_0_yy_yz_zz = buffer_1010_sddd[569];

    auto g_x_0_z_0_0_yy_zz_xx = buffer_1010_sddd[570];

    auto g_x_0_z_0_0_yy_zz_xy = buffer_1010_sddd[571];

    auto g_x_0_z_0_0_yy_zz_xz = buffer_1010_sddd[572];

    auto g_x_0_z_0_0_yy_zz_yy = buffer_1010_sddd[573];

    auto g_x_0_z_0_0_yy_zz_yz = buffer_1010_sddd[574];

    auto g_x_0_z_0_0_yy_zz_zz = buffer_1010_sddd[575];

    auto g_x_0_z_0_0_yz_xx_xx = buffer_1010_sddd[576];

    auto g_x_0_z_0_0_yz_xx_xy = buffer_1010_sddd[577];

    auto g_x_0_z_0_0_yz_xx_xz = buffer_1010_sddd[578];

    auto g_x_0_z_0_0_yz_xx_yy = buffer_1010_sddd[579];

    auto g_x_0_z_0_0_yz_xx_yz = buffer_1010_sddd[580];

    auto g_x_0_z_0_0_yz_xx_zz = buffer_1010_sddd[581];

    auto g_x_0_z_0_0_yz_xy_xx = buffer_1010_sddd[582];

    auto g_x_0_z_0_0_yz_xy_xy = buffer_1010_sddd[583];

    auto g_x_0_z_0_0_yz_xy_xz = buffer_1010_sddd[584];

    auto g_x_0_z_0_0_yz_xy_yy = buffer_1010_sddd[585];

    auto g_x_0_z_0_0_yz_xy_yz = buffer_1010_sddd[586];

    auto g_x_0_z_0_0_yz_xy_zz = buffer_1010_sddd[587];

    auto g_x_0_z_0_0_yz_xz_xx = buffer_1010_sddd[588];

    auto g_x_0_z_0_0_yz_xz_xy = buffer_1010_sddd[589];

    auto g_x_0_z_0_0_yz_xz_xz = buffer_1010_sddd[590];

    auto g_x_0_z_0_0_yz_xz_yy = buffer_1010_sddd[591];

    auto g_x_0_z_0_0_yz_xz_yz = buffer_1010_sddd[592];

    auto g_x_0_z_0_0_yz_xz_zz = buffer_1010_sddd[593];

    auto g_x_0_z_0_0_yz_yy_xx = buffer_1010_sddd[594];

    auto g_x_0_z_0_0_yz_yy_xy = buffer_1010_sddd[595];

    auto g_x_0_z_0_0_yz_yy_xz = buffer_1010_sddd[596];

    auto g_x_0_z_0_0_yz_yy_yy = buffer_1010_sddd[597];

    auto g_x_0_z_0_0_yz_yy_yz = buffer_1010_sddd[598];

    auto g_x_0_z_0_0_yz_yy_zz = buffer_1010_sddd[599];

    auto g_x_0_z_0_0_yz_yz_xx = buffer_1010_sddd[600];

    auto g_x_0_z_0_0_yz_yz_xy = buffer_1010_sddd[601];

    auto g_x_0_z_0_0_yz_yz_xz = buffer_1010_sddd[602];

    auto g_x_0_z_0_0_yz_yz_yy = buffer_1010_sddd[603];

    auto g_x_0_z_0_0_yz_yz_yz = buffer_1010_sddd[604];

    auto g_x_0_z_0_0_yz_yz_zz = buffer_1010_sddd[605];

    auto g_x_0_z_0_0_yz_zz_xx = buffer_1010_sddd[606];

    auto g_x_0_z_0_0_yz_zz_xy = buffer_1010_sddd[607];

    auto g_x_0_z_0_0_yz_zz_xz = buffer_1010_sddd[608];

    auto g_x_0_z_0_0_yz_zz_yy = buffer_1010_sddd[609];

    auto g_x_0_z_0_0_yz_zz_yz = buffer_1010_sddd[610];

    auto g_x_0_z_0_0_yz_zz_zz = buffer_1010_sddd[611];

    auto g_x_0_z_0_0_zz_xx_xx = buffer_1010_sddd[612];

    auto g_x_0_z_0_0_zz_xx_xy = buffer_1010_sddd[613];

    auto g_x_0_z_0_0_zz_xx_xz = buffer_1010_sddd[614];

    auto g_x_0_z_0_0_zz_xx_yy = buffer_1010_sddd[615];

    auto g_x_0_z_0_0_zz_xx_yz = buffer_1010_sddd[616];

    auto g_x_0_z_0_0_zz_xx_zz = buffer_1010_sddd[617];

    auto g_x_0_z_0_0_zz_xy_xx = buffer_1010_sddd[618];

    auto g_x_0_z_0_0_zz_xy_xy = buffer_1010_sddd[619];

    auto g_x_0_z_0_0_zz_xy_xz = buffer_1010_sddd[620];

    auto g_x_0_z_0_0_zz_xy_yy = buffer_1010_sddd[621];

    auto g_x_0_z_0_0_zz_xy_yz = buffer_1010_sddd[622];

    auto g_x_0_z_0_0_zz_xy_zz = buffer_1010_sddd[623];

    auto g_x_0_z_0_0_zz_xz_xx = buffer_1010_sddd[624];

    auto g_x_0_z_0_0_zz_xz_xy = buffer_1010_sddd[625];

    auto g_x_0_z_0_0_zz_xz_xz = buffer_1010_sddd[626];

    auto g_x_0_z_0_0_zz_xz_yy = buffer_1010_sddd[627];

    auto g_x_0_z_0_0_zz_xz_yz = buffer_1010_sddd[628];

    auto g_x_0_z_0_0_zz_xz_zz = buffer_1010_sddd[629];

    auto g_x_0_z_0_0_zz_yy_xx = buffer_1010_sddd[630];

    auto g_x_0_z_0_0_zz_yy_xy = buffer_1010_sddd[631];

    auto g_x_0_z_0_0_zz_yy_xz = buffer_1010_sddd[632];

    auto g_x_0_z_0_0_zz_yy_yy = buffer_1010_sddd[633];

    auto g_x_0_z_0_0_zz_yy_yz = buffer_1010_sddd[634];

    auto g_x_0_z_0_0_zz_yy_zz = buffer_1010_sddd[635];

    auto g_x_0_z_0_0_zz_yz_xx = buffer_1010_sddd[636];

    auto g_x_0_z_0_0_zz_yz_xy = buffer_1010_sddd[637];

    auto g_x_0_z_0_0_zz_yz_xz = buffer_1010_sddd[638];

    auto g_x_0_z_0_0_zz_yz_yy = buffer_1010_sddd[639];

    auto g_x_0_z_0_0_zz_yz_yz = buffer_1010_sddd[640];

    auto g_x_0_z_0_0_zz_yz_zz = buffer_1010_sddd[641];

    auto g_x_0_z_0_0_zz_zz_xx = buffer_1010_sddd[642];

    auto g_x_0_z_0_0_zz_zz_xy = buffer_1010_sddd[643];

    auto g_x_0_z_0_0_zz_zz_xz = buffer_1010_sddd[644];

    auto g_x_0_z_0_0_zz_zz_yy = buffer_1010_sddd[645];

    auto g_x_0_z_0_0_zz_zz_yz = buffer_1010_sddd[646];

    auto g_x_0_z_0_0_zz_zz_zz = buffer_1010_sddd[647];

    auto g_y_0_x_0_0_xx_xx_xx = buffer_1010_sddd[648];

    auto g_y_0_x_0_0_xx_xx_xy = buffer_1010_sddd[649];

    auto g_y_0_x_0_0_xx_xx_xz = buffer_1010_sddd[650];

    auto g_y_0_x_0_0_xx_xx_yy = buffer_1010_sddd[651];

    auto g_y_0_x_0_0_xx_xx_yz = buffer_1010_sddd[652];

    auto g_y_0_x_0_0_xx_xx_zz = buffer_1010_sddd[653];

    auto g_y_0_x_0_0_xx_xy_xx = buffer_1010_sddd[654];

    auto g_y_0_x_0_0_xx_xy_xy = buffer_1010_sddd[655];

    auto g_y_0_x_0_0_xx_xy_xz = buffer_1010_sddd[656];

    auto g_y_0_x_0_0_xx_xy_yy = buffer_1010_sddd[657];

    auto g_y_0_x_0_0_xx_xy_yz = buffer_1010_sddd[658];

    auto g_y_0_x_0_0_xx_xy_zz = buffer_1010_sddd[659];

    auto g_y_0_x_0_0_xx_xz_xx = buffer_1010_sddd[660];

    auto g_y_0_x_0_0_xx_xz_xy = buffer_1010_sddd[661];

    auto g_y_0_x_0_0_xx_xz_xz = buffer_1010_sddd[662];

    auto g_y_0_x_0_0_xx_xz_yy = buffer_1010_sddd[663];

    auto g_y_0_x_0_0_xx_xz_yz = buffer_1010_sddd[664];

    auto g_y_0_x_0_0_xx_xz_zz = buffer_1010_sddd[665];

    auto g_y_0_x_0_0_xx_yy_xx = buffer_1010_sddd[666];

    auto g_y_0_x_0_0_xx_yy_xy = buffer_1010_sddd[667];

    auto g_y_0_x_0_0_xx_yy_xz = buffer_1010_sddd[668];

    auto g_y_0_x_0_0_xx_yy_yy = buffer_1010_sddd[669];

    auto g_y_0_x_0_0_xx_yy_yz = buffer_1010_sddd[670];

    auto g_y_0_x_0_0_xx_yy_zz = buffer_1010_sddd[671];

    auto g_y_0_x_0_0_xx_yz_xx = buffer_1010_sddd[672];

    auto g_y_0_x_0_0_xx_yz_xy = buffer_1010_sddd[673];

    auto g_y_0_x_0_0_xx_yz_xz = buffer_1010_sddd[674];

    auto g_y_0_x_0_0_xx_yz_yy = buffer_1010_sddd[675];

    auto g_y_0_x_0_0_xx_yz_yz = buffer_1010_sddd[676];

    auto g_y_0_x_0_0_xx_yz_zz = buffer_1010_sddd[677];

    auto g_y_0_x_0_0_xx_zz_xx = buffer_1010_sddd[678];

    auto g_y_0_x_0_0_xx_zz_xy = buffer_1010_sddd[679];

    auto g_y_0_x_0_0_xx_zz_xz = buffer_1010_sddd[680];

    auto g_y_0_x_0_0_xx_zz_yy = buffer_1010_sddd[681];

    auto g_y_0_x_0_0_xx_zz_yz = buffer_1010_sddd[682];

    auto g_y_0_x_0_0_xx_zz_zz = buffer_1010_sddd[683];

    auto g_y_0_x_0_0_xy_xx_xx = buffer_1010_sddd[684];

    auto g_y_0_x_0_0_xy_xx_xy = buffer_1010_sddd[685];

    auto g_y_0_x_0_0_xy_xx_xz = buffer_1010_sddd[686];

    auto g_y_0_x_0_0_xy_xx_yy = buffer_1010_sddd[687];

    auto g_y_0_x_0_0_xy_xx_yz = buffer_1010_sddd[688];

    auto g_y_0_x_0_0_xy_xx_zz = buffer_1010_sddd[689];

    auto g_y_0_x_0_0_xy_xy_xx = buffer_1010_sddd[690];

    auto g_y_0_x_0_0_xy_xy_xy = buffer_1010_sddd[691];

    auto g_y_0_x_0_0_xy_xy_xz = buffer_1010_sddd[692];

    auto g_y_0_x_0_0_xy_xy_yy = buffer_1010_sddd[693];

    auto g_y_0_x_0_0_xy_xy_yz = buffer_1010_sddd[694];

    auto g_y_0_x_0_0_xy_xy_zz = buffer_1010_sddd[695];

    auto g_y_0_x_0_0_xy_xz_xx = buffer_1010_sddd[696];

    auto g_y_0_x_0_0_xy_xz_xy = buffer_1010_sddd[697];

    auto g_y_0_x_0_0_xy_xz_xz = buffer_1010_sddd[698];

    auto g_y_0_x_0_0_xy_xz_yy = buffer_1010_sddd[699];

    auto g_y_0_x_0_0_xy_xz_yz = buffer_1010_sddd[700];

    auto g_y_0_x_0_0_xy_xz_zz = buffer_1010_sddd[701];

    auto g_y_0_x_0_0_xy_yy_xx = buffer_1010_sddd[702];

    auto g_y_0_x_0_0_xy_yy_xy = buffer_1010_sddd[703];

    auto g_y_0_x_0_0_xy_yy_xz = buffer_1010_sddd[704];

    auto g_y_0_x_0_0_xy_yy_yy = buffer_1010_sddd[705];

    auto g_y_0_x_0_0_xy_yy_yz = buffer_1010_sddd[706];

    auto g_y_0_x_0_0_xy_yy_zz = buffer_1010_sddd[707];

    auto g_y_0_x_0_0_xy_yz_xx = buffer_1010_sddd[708];

    auto g_y_0_x_0_0_xy_yz_xy = buffer_1010_sddd[709];

    auto g_y_0_x_0_0_xy_yz_xz = buffer_1010_sddd[710];

    auto g_y_0_x_0_0_xy_yz_yy = buffer_1010_sddd[711];

    auto g_y_0_x_0_0_xy_yz_yz = buffer_1010_sddd[712];

    auto g_y_0_x_0_0_xy_yz_zz = buffer_1010_sddd[713];

    auto g_y_0_x_0_0_xy_zz_xx = buffer_1010_sddd[714];

    auto g_y_0_x_0_0_xy_zz_xy = buffer_1010_sddd[715];

    auto g_y_0_x_0_0_xy_zz_xz = buffer_1010_sddd[716];

    auto g_y_0_x_0_0_xy_zz_yy = buffer_1010_sddd[717];

    auto g_y_0_x_0_0_xy_zz_yz = buffer_1010_sddd[718];

    auto g_y_0_x_0_0_xy_zz_zz = buffer_1010_sddd[719];

    auto g_y_0_x_0_0_xz_xx_xx = buffer_1010_sddd[720];

    auto g_y_0_x_0_0_xz_xx_xy = buffer_1010_sddd[721];

    auto g_y_0_x_0_0_xz_xx_xz = buffer_1010_sddd[722];

    auto g_y_0_x_0_0_xz_xx_yy = buffer_1010_sddd[723];

    auto g_y_0_x_0_0_xz_xx_yz = buffer_1010_sddd[724];

    auto g_y_0_x_0_0_xz_xx_zz = buffer_1010_sddd[725];

    auto g_y_0_x_0_0_xz_xy_xx = buffer_1010_sddd[726];

    auto g_y_0_x_0_0_xz_xy_xy = buffer_1010_sddd[727];

    auto g_y_0_x_0_0_xz_xy_xz = buffer_1010_sddd[728];

    auto g_y_0_x_0_0_xz_xy_yy = buffer_1010_sddd[729];

    auto g_y_0_x_0_0_xz_xy_yz = buffer_1010_sddd[730];

    auto g_y_0_x_0_0_xz_xy_zz = buffer_1010_sddd[731];

    auto g_y_0_x_0_0_xz_xz_xx = buffer_1010_sddd[732];

    auto g_y_0_x_0_0_xz_xz_xy = buffer_1010_sddd[733];

    auto g_y_0_x_0_0_xz_xz_xz = buffer_1010_sddd[734];

    auto g_y_0_x_0_0_xz_xz_yy = buffer_1010_sddd[735];

    auto g_y_0_x_0_0_xz_xz_yz = buffer_1010_sddd[736];

    auto g_y_0_x_0_0_xz_xz_zz = buffer_1010_sddd[737];

    auto g_y_0_x_0_0_xz_yy_xx = buffer_1010_sddd[738];

    auto g_y_0_x_0_0_xz_yy_xy = buffer_1010_sddd[739];

    auto g_y_0_x_0_0_xz_yy_xz = buffer_1010_sddd[740];

    auto g_y_0_x_0_0_xz_yy_yy = buffer_1010_sddd[741];

    auto g_y_0_x_0_0_xz_yy_yz = buffer_1010_sddd[742];

    auto g_y_0_x_0_0_xz_yy_zz = buffer_1010_sddd[743];

    auto g_y_0_x_0_0_xz_yz_xx = buffer_1010_sddd[744];

    auto g_y_0_x_0_0_xz_yz_xy = buffer_1010_sddd[745];

    auto g_y_0_x_0_0_xz_yz_xz = buffer_1010_sddd[746];

    auto g_y_0_x_0_0_xz_yz_yy = buffer_1010_sddd[747];

    auto g_y_0_x_0_0_xz_yz_yz = buffer_1010_sddd[748];

    auto g_y_0_x_0_0_xz_yz_zz = buffer_1010_sddd[749];

    auto g_y_0_x_0_0_xz_zz_xx = buffer_1010_sddd[750];

    auto g_y_0_x_0_0_xz_zz_xy = buffer_1010_sddd[751];

    auto g_y_0_x_0_0_xz_zz_xz = buffer_1010_sddd[752];

    auto g_y_0_x_0_0_xz_zz_yy = buffer_1010_sddd[753];

    auto g_y_0_x_0_0_xz_zz_yz = buffer_1010_sddd[754];

    auto g_y_0_x_0_0_xz_zz_zz = buffer_1010_sddd[755];

    auto g_y_0_x_0_0_yy_xx_xx = buffer_1010_sddd[756];

    auto g_y_0_x_0_0_yy_xx_xy = buffer_1010_sddd[757];

    auto g_y_0_x_0_0_yy_xx_xz = buffer_1010_sddd[758];

    auto g_y_0_x_0_0_yy_xx_yy = buffer_1010_sddd[759];

    auto g_y_0_x_0_0_yy_xx_yz = buffer_1010_sddd[760];

    auto g_y_0_x_0_0_yy_xx_zz = buffer_1010_sddd[761];

    auto g_y_0_x_0_0_yy_xy_xx = buffer_1010_sddd[762];

    auto g_y_0_x_0_0_yy_xy_xy = buffer_1010_sddd[763];

    auto g_y_0_x_0_0_yy_xy_xz = buffer_1010_sddd[764];

    auto g_y_0_x_0_0_yy_xy_yy = buffer_1010_sddd[765];

    auto g_y_0_x_0_0_yy_xy_yz = buffer_1010_sddd[766];

    auto g_y_0_x_0_0_yy_xy_zz = buffer_1010_sddd[767];

    auto g_y_0_x_0_0_yy_xz_xx = buffer_1010_sddd[768];

    auto g_y_0_x_0_0_yy_xz_xy = buffer_1010_sddd[769];

    auto g_y_0_x_0_0_yy_xz_xz = buffer_1010_sddd[770];

    auto g_y_0_x_0_0_yy_xz_yy = buffer_1010_sddd[771];

    auto g_y_0_x_0_0_yy_xz_yz = buffer_1010_sddd[772];

    auto g_y_0_x_0_0_yy_xz_zz = buffer_1010_sddd[773];

    auto g_y_0_x_0_0_yy_yy_xx = buffer_1010_sddd[774];

    auto g_y_0_x_0_0_yy_yy_xy = buffer_1010_sddd[775];

    auto g_y_0_x_0_0_yy_yy_xz = buffer_1010_sddd[776];

    auto g_y_0_x_0_0_yy_yy_yy = buffer_1010_sddd[777];

    auto g_y_0_x_0_0_yy_yy_yz = buffer_1010_sddd[778];

    auto g_y_0_x_0_0_yy_yy_zz = buffer_1010_sddd[779];

    auto g_y_0_x_0_0_yy_yz_xx = buffer_1010_sddd[780];

    auto g_y_0_x_0_0_yy_yz_xy = buffer_1010_sddd[781];

    auto g_y_0_x_0_0_yy_yz_xz = buffer_1010_sddd[782];

    auto g_y_0_x_0_0_yy_yz_yy = buffer_1010_sddd[783];

    auto g_y_0_x_0_0_yy_yz_yz = buffer_1010_sddd[784];

    auto g_y_0_x_0_0_yy_yz_zz = buffer_1010_sddd[785];

    auto g_y_0_x_0_0_yy_zz_xx = buffer_1010_sddd[786];

    auto g_y_0_x_0_0_yy_zz_xy = buffer_1010_sddd[787];

    auto g_y_0_x_0_0_yy_zz_xz = buffer_1010_sddd[788];

    auto g_y_0_x_0_0_yy_zz_yy = buffer_1010_sddd[789];

    auto g_y_0_x_0_0_yy_zz_yz = buffer_1010_sddd[790];

    auto g_y_0_x_0_0_yy_zz_zz = buffer_1010_sddd[791];

    auto g_y_0_x_0_0_yz_xx_xx = buffer_1010_sddd[792];

    auto g_y_0_x_0_0_yz_xx_xy = buffer_1010_sddd[793];

    auto g_y_0_x_0_0_yz_xx_xz = buffer_1010_sddd[794];

    auto g_y_0_x_0_0_yz_xx_yy = buffer_1010_sddd[795];

    auto g_y_0_x_0_0_yz_xx_yz = buffer_1010_sddd[796];

    auto g_y_0_x_0_0_yz_xx_zz = buffer_1010_sddd[797];

    auto g_y_0_x_0_0_yz_xy_xx = buffer_1010_sddd[798];

    auto g_y_0_x_0_0_yz_xy_xy = buffer_1010_sddd[799];

    auto g_y_0_x_0_0_yz_xy_xz = buffer_1010_sddd[800];

    auto g_y_0_x_0_0_yz_xy_yy = buffer_1010_sddd[801];

    auto g_y_0_x_0_0_yz_xy_yz = buffer_1010_sddd[802];

    auto g_y_0_x_0_0_yz_xy_zz = buffer_1010_sddd[803];

    auto g_y_0_x_0_0_yz_xz_xx = buffer_1010_sddd[804];

    auto g_y_0_x_0_0_yz_xz_xy = buffer_1010_sddd[805];

    auto g_y_0_x_0_0_yz_xz_xz = buffer_1010_sddd[806];

    auto g_y_0_x_0_0_yz_xz_yy = buffer_1010_sddd[807];

    auto g_y_0_x_0_0_yz_xz_yz = buffer_1010_sddd[808];

    auto g_y_0_x_0_0_yz_xz_zz = buffer_1010_sddd[809];

    auto g_y_0_x_0_0_yz_yy_xx = buffer_1010_sddd[810];

    auto g_y_0_x_0_0_yz_yy_xy = buffer_1010_sddd[811];

    auto g_y_0_x_0_0_yz_yy_xz = buffer_1010_sddd[812];

    auto g_y_0_x_0_0_yz_yy_yy = buffer_1010_sddd[813];

    auto g_y_0_x_0_0_yz_yy_yz = buffer_1010_sddd[814];

    auto g_y_0_x_0_0_yz_yy_zz = buffer_1010_sddd[815];

    auto g_y_0_x_0_0_yz_yz_xx = buffer_1010_sddd[816];

    auto g_y_0_x_0_0_yz_yz_xy = buffer_1010_sddd[817];

    auto g_y_0_x_0_0_yz_yz_xz = buffer_1010_sddd[818];

    auto g_y_0_x_0_0_yz_yz_yy = buffer_1010_sddd[819];

    auto g_y_0_x_0_0_yz_yz_yz = buffer_1010_sddd[820];

    auto g_y_0_x_0_0_yz_yz_zz = buffer_1010_sddd[821];

    auto g_y_0_x_0_0_yz_zz_xx = buffer_1010_sddd[822];

    auto g_y_0_x_0_0_yz_zz_xy = buffer_1010_sddd[823];

    auto g_y_0_x_0_0_yz_zz_xz = buffer_1010_sddd[824];

    auto g_y_0_x_0_0_yz_zz_yy = buffer_1010_sddd[825];

    auto g_y_0_x_0_0_yz_zz_yz = buffer_1010_sddd[826];

    auto g_y_0_x_0_0_yz_zz_zz = buffer_1010_sddd[827];

    auto g_y_0_x_0_0_zz_xx_xx = buffer_1010_sddd[828];

    auto g_y_0_x_0_0_zz_xx_xy = buffer_1010_sddd[829];

    auto g_y_0_x_0_0_zz_xx_xz = buffer_1010_sddd[830];

    auto g_y_0_x_0_0_zz_xx_yy = buffer_1010_sddd[831];

    auto g_y_0_x_0_0_zz_xx_yz = buffer_1010_sddd[832];

    auto g_y_0_x_0_0_zz_xx_zz = buffer_1010_sddd[833];

    auto g_y_0_x_0_0_zz_xy_xx = buffer_1010_sddd[834];

    auto g_y_0_x_0_0_zz_xy_xy = buffer_1010_sddd[835];

    auto g_y_0_x_0_0_zz_xy_xz = buffer_1010_sddd[836];

    auto g_y_0_x_0_0_zz_xy_yy = buffer_1010_sddd[837];

    auto g_y_0_x_0_0_zz_xy_yz = buffer_1010_sddd[838];

    auto g_y_0_x_0_0_zz_xy_zz = buffer_1010_sddd[839];

    auto g_y_0_x_0_0_zz_xz_xx = buffer_1010_sddd[840];

    auto g_y_0_x_0_0_zz_xz_xy = buffer_1010_sddd[841];

    auto g_y_0_x_0_0_zz_xz_xz = buffer_1010_sddd[842];

    auto g_y_0_x_0_0_zz_xz_yy = buffer_1010_sddd[843];

    auto g_y_0_x_0_0_zz_xz_yz = buffer_1010_sddd[844];

    auto g_y_0_x_0_0_zz_xz_zz = buffer_1010_sddd[845];

    auto g_y_0_x_0_0_zz_yy_xx = buffer_1010_sddd[846];

    auto g_y_0_x_0_0_zz_yy_xy = buffer_1010_sddd[847];

    auto g_y_0_x_0_0_zz_yy_xz = buffer_1010_sddd[848];

    auto g_y_0_x_0_0_zz_yy_yy = buffer_1010_sddd[849];

    auto g_y_0_x_0_0_zz_yy_yz = buffer_1010_sddd[850];

    auto g_y_0_x_0_0_zz_yy_zz = buffer_1010_sddd[851];

    auto g_y_0_x_0_0_zz_yz_xx = buffer_1010_sddd[852];

    auto g_y_0_x_0_0_zz_yz_xy = buffer_1010_sddd[853];

    auto g_y_0_x_0_0_zz_yz_xz = buffer_1010_sddd[854];

    auto g_y_0_x_0_0_zz_yz_yy = buffer_1010_sddd[855];

    auto g_y_0_x_0_0_zz_yz_yz = buffer_1010_sddd[856];

    auto g_y_0_x_0_0_zz_yz_zz = buffer_1010_sddd[857];

    auto g_y_0_x_0_0_zz_zz_xx = buffer_1010_sddd[858];

    auto g_y_0_x_0_0_zz_zz_xy = buffer_1010_sddd[859];

    auto g_y_0_x_0_0_zz_zz_xz = buffer_1010_sddd[860];

    auto g_y_0_x_0_0_zz_zz_yy = buffer_1010_sddd[861];

    auto g_y_0_x_0_0_zz_zz_yz = buffer_1010_sddd[862];

    auto g_y_0_x_0_0_zz_zz_zz = buffer_1010_sddd[863];

    auto g_y_0_y_0_0_xx_xx_xx = buffer_1010_sddd[864];

    auto g_y_0_y_0_0_xx_xx_xy = buffer_1010_sddd[865];

    auto g_y_0_y_0_0_xx_xx_xz = buffer_1010_sddd[866];

    auto g_y_0_y_0_0_xx_xx_yy = buffer_1010_sddd[867];

    auto g_y_0_y_0_0_xx_xx_yz = buffer_1010_sddd[868];

    auto g_y_0_y_0_0_xx_xx_zz = buffer_1010_sddd[869];

    auto g_y_0_y_0_0_xx_xy_xx = buffer_1010_sddd[870];

    auto g_y_0_y_0_0_xx_xy_xy = buffer_1010_sddd[871];

    auto g_y_0_y_0_0_xx_xy_xz = buffer_1010_sddd[872];

    auto g_y_0_y_0_0_xx_xy_yy = buffer_1010_sddd[873];

    auto g_y_0_y_0_0_xx_xy_yz = buffer_1010_sddd[874];

    auto g_y_0_y_0_0_xx_xy_zz = buffer_1010_sddd[875];

    auto g_y_0_y_0_0_xx_xz_xx = buffer_1010_sddd[876];

    auto g_y_0_y_0_0_xx_xz_xy = buffer_1010_sddd[877];

    auto g_y_0_y_0_0_xx_xz_xz = buffer_1010_sddd[878];

    auto g_y_0_y_0_0_xx_xz_yy = buffer_1010_sddd[879];

    auto g_y_0_y_0_0_xx_xz_yz = buffer_1010_sddd[880];

    auto g_y_0_y_0_0_xx_xz_zz = buffer_1010_sddd[881];

    auto g_y_0_y_0_0_xx_yy_xx = buffer_1010_sddd[882];

    auto g_y_0_y_0_0_xx_yy_xy = buffer_1010_sddd[883];

    auto g_y_0_y_0_0_xx_yy_xz = buffer_1010_sddd[884];

    auto g_y_0_y_0_0_xx_yy_yy = buffer_1010_sddd[885];

    auto g_y_0_y_0_0_xx_yy_yz = buffer_1010_sddd[886];

    auto g_y_0_y_0_0_xx_yy_zz = buffer_1010_sddd[887];

    auto g_y_0_y_0_0_xx_yz_xx = buffer_1010_sddd[888];

    auto g_y_0_y_0_0_xx_yz_xy = buffer_1010_sddd[889];

    auto g_y_0_y_0_0_xx_yz_xz = buffer_1010_sddd[890];

    auto g_y_0_y_0_0_xx_yz_yy = buffer_1010_sddd[891];

    auto g_y_0_y_0_0_xx_yz_yz = buffer_1010_sddd[892];

    auto g_y_0_y_0_0_xx_yz_zz = buffer_1010_sddd[893];

    auto g_y_0_y_0_0_xx_zz_xx = buffer_1010_sddd[894];

    auto g_y_0_y_0_0_xx_zz_xy = buffer_1010_sddd[895];

    auto g_y_0_y_0_0_xx_zz_xz = buffer_1010_sddd[896];

    auto g_y_0_y_0_0_xx_zz_yy = buffer_1010_sddd[897];

    auto g_y_0_y_0_0_xx_zz_yz = buffer_1010_sddd[898];

    auto g_y_0_y_0_0_xx_zz_zz = buffer_1010_sddd[899];

    auto g_y_0_y_0_0_xy_xx_xx = buffer_1010_sddd[900];

    auto g_y_0_y_0_0_xy_xx_xy = buffer_1010_sddd[901];

    auto g_y_0_y_0_0_xy_xx_xz = buffer_1010_sddd[902];

    auto g_y_0_y_0_0_xy_xx_yy = buffer_1010_sddd[903];

    auto g_y_0_y_0_0_xy_xx_yz = buffer_1010_sddd[904];

    auto g_y_0_y_0_0_xy_xx_zz = buffer_1010_sddd[905];

    auto g_y_0_y_0_0_xy_xy_xx = buffer_1010_sddd[906];

    auto g_y_0_y_0_0_xy_xy_xy = buffer_1010_sddd[907];

    auto g_y_0_y_0_0_xy_xy_xz = buffer_1010_sddd[908];

    auto g_y_0_y_0_0_xy_xy_yy = buffer_1010_sddd[909];

    auto g_y_0_y_0_0_xy_xy_yz = buffer_1010_sddd[910];

    auto g_y_0_y_0_0_xy_xy_zz = buffer_1010_sddd[911];

    auto g_y_0_y_0_0_xy_xz_xx = buffer_1010_sddd[912];

    auto g_y_0_y_0_0_xy_xz_xy = buffer_1010_sddd[913];

    auto g_y_0_y_0_0_xy_xz_xz = buffer_1010_sddd[914];

    auto g_y_0_y_0_0_xy_xz_yy = buffer_1010_sddd[915];

    auto g_y_0_y_0_0_xy_xz_yz = buffer_1010_sddd[916];

    auto g_y_0_y_0_0_xy_xz_zz = buffer_1010_sddd[917];

    auto g_y_0_y_0_0_xy_yy_xx = buffer_1010_sddd[918];

    auto g_y_0_y_0_0_xy_yy_xy = buffer_1010_sddd[919];

    auto g_y_0_y_0_0_xy_yy_xz = buffer_1010_sddd[920];

    auto g_y_0_y_0_0_xy_yy_yy = buffer_1010_sddd[921];

    auto g_y_0_y_0_0_xy_yy_yz = buffer_1010_sddd[922];

    auto g_y_0_y_0_0_xy_yy_zz = buffer_1010_sddd[923];

    auto g_y_0_y_0_0_xy_yz_xx = buffer_1010_sddd[924];

    auto g_y_0_y_0_0_xy_yz_xy = buffer_1010_sddd[925];

    auto g_y_0_y_0_0_xy_yz_xz = buffer_1010_sddd[926];

    auto g_y_0_y_0_0_xy_yz_yy = buffer_1010_sddd[927];

    auto g_y_0_y_0_0_xy_yz_yz = buffer_1010_sddd[928];

    auto g_y_0_y_0_0_xy_yz_zz = buffer_1010_sddd[929];

    auto g_y_0_y_0_0_xy_zz_xx = buffer_1010_sddd[930];

    auto g_y_0_y_0_0_xy_zz_xy = buffer_1010_sddd[931];

    auto g_y_0_y_0_0_xy_zz_xz = buffer_1010_sddd[932];

    auto g_y_0_y_0_0_xy_zz_yy = buffer_1010_sddd[933];

    auto g_y_0_y_0_0_xy_zz_yz = buffer_1010_sddd[934];

    auto g_y_0_y_0_0_xy_zz_zz = buffer_1010_sddd[935];

    auto g_y_0_y_0_0_xz_xx_xx = buffer_1010_sddd[936];

    auto g_y_0_y_0_0_xz_xx_xy = buffer_1010_sddd[937];

    auto g_y_0_y_0_0_xz_xx_xz = buffer_1010_sddd[938];

    auto g_y_0_y_0_0_xz_xx_yy = buffer_1010_sddd[939];

    auto g_y_0_y_0_0_xz_xx_yz = buffer_1010_sddd[940];

    auto g_y_0_y_0_0_xz_xx_zz = buffer_1010_sddd[941];

    auto g_y_0_y_0_0_xz_xy_xx = buffer_1010_sddd[942];

    auto g_y_0_y_0_0_xz_xy_xy = buffer_1010_sddd[943];

    auto g_y_0_y_0_0_xz_xy_xz = buffer_1010_sddd[944];

    auto g_y_0_y_0_0_xz_xy_yy = buffer_1010_sddd[945];

    auto g_y_0_y_0_0_xz_xy_yz = buffer_1010_sddd[946];

    auto g_y_0_y_0_0_xz_xy_zz = buffer_1010_sddd[947];

    auto g_y_0_y_0_0_xz_xz_xx = buffer_1010_sddd[948];

    auto g_y_0_y_0_0_xz_xz_xy = buffer_1010_sddd[949];

    auto g_y_0_y_0_0_xz_xz_xz = buffer_1010_sddd[950];

    auto g_y_0_y_0_0_xz_xz_yy = buffer_1010_sddd[951];

    auto g_y_0_y_0_0_xz_xz_yz = buffer_1010_sddd[952];

    auto g_y_0_y_0_0_xz_xz_zz = buffer_1010_sddd[953];

    auto g_y_0_y_0_0_xz_yy_xx = buffer_1010_sddd[954];

    auto g_y_0_y_0_0_xz_yy_xy = buffer_1010_sddd[955];

    auto g_y_0_y_0_0_xz_yy_xz = buffer_1010_sddd[956];

    auto g_y_0_y_0_0_xz_yy_yy = buffer_1010_sddd[957];

    auto g_y_0_y_0_0_xz_yy_yz = buffer_1010_sddd[958];

    auto g_y_0_y_0_0_xz_yy_zz = buffer_1010_sddd[959];

    auto g_y_0_y_0_0_xz_yz_xx = buffer_1010_sddd[960];

    auto g_y_0_y_0_0_xz_yz_xy = buffer_1010_sddd[961];

    auto g_y_0_y_0_0_xz_yz_xz = buffer_1010_sddd[962];

    auto g_y_0_y_0_0_xz_yz_yy = buffer_1010_sddd[963];

    auto g_y_0_y_0_0_xz_yz_yz = buffer_1010_sddd[964];

    auto g_y_0_y_0_0_xz_yz_zz = buffer_1010_sddd[965];

    auto g_y_0_y_0_0_xz_zz_xx = buffer_1010_sddd[966];

    auto g_y_0_y_0_0_xz_zz_xy = buffer_1010_sddd[967];

    auto g_y_0_y_0_0_xz_zz_xz = buffer_1010_sddd[968];

    auto g_y_0_y_0_0_xz_zz_yy = buffer_1010_sddd[969];

    auto g_y_0_y_0_0_xz_zz_yz = buffer_1010_sddd[970];

    auto g_y_0_y_0_0_xz_zz_zz = buffer_1010_sddd[971];

    auto g_y_0_y_0_0_yy_xx_xx = buffer_1010_sddd[972];

    auto g_y_0_y_0_0_yy_xx_xy = buffer_1010_sddd[973];

    auto g_y_0_y_0_0_yy_xx_xz = buffer_1010_sddd[974];

    auto g_y_0_y_0_0_yy_xx_yy = buffer_1010_sddd[975];

    auto g_y_0_y_0_0_yy_xx_yz = buffer_1010_sddd[976];

    auto g_y_0_y_0_0_yy_xx_zz = buffer_1010_sddd[977];

    auto g_y_0_y_0_0_yy_xy_xx = buffer_1010_sddd[978];

    auto g_y_0_y_0_0_yy_xy_xy = buffer_1010_sddd[979];

    auto g_y_0_y_0_0_yy_xy_xz = buffer_1010_sddd[980];

    auto g_y_0_y_0_0_yy_xy_yy = buffer_1010_sddd[981];

    auto g_y_0_y_0_0_yy_xy_yz = buffer_1010_sddd[982];

    auto g_y_0_y_0_0_yy_xy_zz = buffer_1010_sddd[983];

    auto g_y_0_y_0_0_yy_xz_xx = buffer_1010_sddd[984];

    auto g_y_0_y_0_0_yy_xz_xy = buffer_1010_sddd[985];

    auto g_y_0_y_0_0_yy_xz_xz = buffer_1010_sddd[986];

    auto g_y_0_y_0_0_yy_xz_yy = buffer_1010_sddd[987];

    auto g_y_0_y_0_0_yy_xz_yz = buffer_1010_sddd[988];

    auto g_y_0_y_0_0_yy_xz_zz = buffer_1010_sddd[989];

    auto g_y_0_y_0_0_yy_yy_xx = buffer_1010_sddd[990];

    auto g_y_0_y_0_0_yy_yy_xy = buffer_1010_sddd[991];

    auto g_y_0_y_0_0_yy_yy_xz = buffer_1010_sddd[992];

    auto g_y_0_y_0_0_yy_yy_yy = buffer_1010_sddd[993];

    auto g_y_0_y_0_0_yy_yy_yz = buffer_1010_sddd[994];

    auto g_y_0_y_0_0_yy_yy_zz = buffer_1010_sddd[995];

    auto g_y_0_y_0_0_yy_yz_xx = buffer_1010_sddd[996];

    auto g_y_0_y_0_0_yy_yz_xy = buffer_1010_sddd[997];

    auto g_y_0_y_0_0_yy_yz_xz = buffer_1010_sddd[998];

    auto g_y_0_y_0_0_yy_yz_yy = buffer_1010_sddd[999];

    auto g_y_0_y_0_0_yy_yz_yz = buffer_1010_sddd[1000];

    auto g_y_0_y_0_0_yy_yz_zz = buffer_1010_sddd[1001];

    auto g_y_0_y_0_0_yy_zz_xx = buffer_1010_sddd[1002];

    auto g_y_0_y_0_0_yy_zz_xy = buffer_1010_sddd[1003];

    auto g_y_0_y_0_0_yy_zz_xz = buffer_1010_sddd[1004];

    auto g_y_0_y_0_0_yy_zz_yy = buffer_1010_sddd[1005];

    auto g_y_0_y_0_0_yy_zz_yz = buffer_1010_sddd[1006];

    auto g_y_0_y_0_0_yy_zz_zz = buffer_1010_sddd[1007];

    auto g_y_0_y_0_0_yz_xx_xx = buffer_1010_sddd[1008];

    auto g_y_0_y_0_0_yz_xx_xy = buffer_1010_sddd[1009];

    auto g_y_0_y_0_0_yz_xx_xz = buffer_1010_sddd[1010];

    auto g_y_0_y_0_0_yz_xx_yy = buffer_1010_sddd[1011];

    auto g_y_0_y_0_0_yz_xx_yz = buffer_1010_sddd[1012];

    auto g_y_0_y_0_0_yz_xx_zz = buffer_1010_sddd[1013];

    auto g_y_0_y_0_0_yz_xy_xx = buffer_1010_sddd[1014];

    auto g_y_0_y_0_0_yz_xy_xy = buffer_1010_sddd[1015];

    auto g_y_0_y_0_0_yz_xy_xz = buffer_1010_sddd[1016];

    auto g_y_0_y_0_0_yz_xy_yy = buffer_1010_sddd[1017];

    auto g_y_0_y_0_0_yz_xy_yz = buffer_1010_sddd[1018];

    auto g_y_0_y_0_0_yz_xy_zz = buffer_1010_sddd[1019];

    auto g_y_0_y_0_0_yz_xz_xx = buffer_1010_sddd[1020];

    auto g_y_0_y_0_0_yz_xz_xy = buffer_1010_sddd[1021];

    auto g_y_0_y_0_0_yz_xz_xz = buffer_1010_sddd[1022];

    auto g_y_0_y_0_0_yz_xz_yy = buffer_1010_sddd[1023];

    auto g_y_0_y_0_0_yz_xz_yz = buffer_1010_sddd[1024];

    auto g_y_0_y_0_0_yz_xz_zz = buffer_1010_sddd[1025];

    auto g_y_0_y_0_0_yz_yy_xx = buffer_1010_sddd[1026];

    auto g_y_0_y_0_0_yz_yy_xy = buffer_1010_sddd[1027];

    auto g_y_0_y_0_0_yz_yy_xz = buffer_1010_sddd[1028];

    auto g_y_0_y_0_0_yz_yy_yy = buffer_1010_sddd[1029];

    auto g_y_0_y_0_0_yz_yy_yz = buffer_1010_sddd[1030];

    auto g_y_0_y_0_0_yz_yy_zz = buffer_1010_sddd[1031];

    auto g_y_0_y_0_0_yz_yz_xx = buffer_1010_sddd[1032];

    auto g_y_0_y_0_0_yz_yz_xy = buffer_1010_sddd[1033];

    auto g_y_0_y_0_0_yz_yz_xz = buffer_1010_sddd[1034];

    auto g_y_0_y_0_0_yz_yz_yy = buffer_1010_sddd[1035];

    auto g_y_0_y_0_0_yz_yz_yz = buffer_1010_sddd[1036];

    auto g_y_0_y_0_0_yz_yz_zz = buffer_1010_sddd[1037];

    auto g_y_0_y_0_0_yz_zz_xx = buffer_1010_sddd[1038];

    auto g_y_0_y_0_0_yz_zz_xy = buffer_1010_sddd[1039];

    auto g_y_0_y_0_0_yz_zz_xz = buffer_1010_sddd[1040];

    auto g_y_0_y_0_0_yz_zz_yy = buffer_1010_sddd[1041];

    auto g_y_0_y_0_0_yz_zz_yz = buffer_1010_sddd[1042];

    auto g_y_0_y_0_0_yz_zz_zz = buffer_1010_sddd[1043];

    auto g_y_0_y_0_0_zz_xx_xx = buffer_1010_sddd[1044];

    auto g_y_0_y_0_0_zz_xx_xy = buffer_1010_sddd[1045];

    auto g_y_0_y_0_0_zz_xx_xz = buffer_1010_sddd[1046];

    auto g_y_0_y_0_0_zz_xx_yy = buffer_1010_sddd[1047];

    auto g_y_0_y_0_0_zz_xx_yz = buffer_1010_sddd[1048];

    auto g_y_0_y_0_0_zz_xx_zz = buffer_1010_sddd[1049];

    auto g_y_0_y_0_0_zz_xy_xx = buffer_1010_sddd[1050];

    auto g_y_0_y_0_0_zz_xy_xy = buffer_1010_sddd[1051];

    auto g_y_0_y_0_0_zz_xy_xz = buffer_1010_sddd[1052];

    auto g_y_0_y_0_0_zz_xy_yy = buffer_1010_sddd[1053];

    auto g_y_0_y_0_0_zz_xy_yz = buffer_1010_sddd[1054];

    auto g_y_0_y_0_0_zz_xy_zz = buffer_1010_sddd[1055];

    auto g_y_0_y_0_0_zz_xz_xx = buffer_1010_sddd[1056];

    auto g_y_0_y_0_0_zz_xz_xy = buffer_1010_sddd[1057];

    auto g_y_0_y_0_0_zz_xz_xz = buffer_1010_sddd[1058];

    auto g_y_0_y_0_0_zz_xz_yy = buffer_1010_sddd[1059];

    auto g_y_0_y_0_0_zz_xz_yz = buffer_1010_sddd[1060];

    auto g_y_0_y_0_0_zz_xz_zz = buffer_1010_sddd[1061];

    auto g_y_0_y_0_0_zz_yy_xx = buffer_1010_sddd[1062];

    auto g_y_0_y_0_0_zz_yy_xy = buffer_1010_sddd[1063];

    auto g_y_0_y_0_0_zz_yy_xz = buffer_1010_sddd[1064];

    auto g_y_0_y_0_0_zz_yy_yy = buffer_1010_sddd[1065];

    auto g_y_0_y_0_0_zz_yy_yz = buffer_1010_sddd[1066];

    auto g_y_0_y_0_0_zz_yy_zz = buffer_1010_sddd[1067];

    auto g_y_0_y_0_0_zz_yz_xx = buffer_1010_sddd[1068];

    auto g_y_0_y_0_0_zz_yz_xy = buffer_1010_sddd[1069];

    auto g_y_0_y_0_0_zz_yz_xz = buffer_1010_sddd[1070];

    auto g_y_0_y_0_0_zz_yz_yy = buffer_1010_sddd[1071];

    auto g_y_0_y_0_0_zz_yz_yz = buffer_1010_sddd[1072];

    auto g_y_0_y_0_0_zz_yz_zz = buffer_1010_sddd[1073];

    auto g_y_0_y_0_0_zz_zz_xx = buffer_1010_sddd[1074];

    auto g_y_0_y_0_0_zz_zz_xy = buffer_1010_sddd[1075];

    auto g_y_0_y_0_0_zz_zz_xz = buffer_1010_sddd[1076];

    auto g_y_0_y_0_0_zz_zz_yy = buffer_1010_sddd[1077];

    auto g_y_0_y_0_0_zz_zz_yz = buffer_1010_sddd[1078];

    auto g_y_0_y_0_0_zz_zz_zz = buffer_1010_sddd[1079];

    auto g_y_0_z_0_0_xx_xx_xx = buffer_1010_sddd[1080];

    auto g_y_0_z_0_0_xx_xx_xy = buffer_1010_sddd[1081];

    auto g_y_0_z_0_0_xx_xx_xz = buffer_1010_sddd[1082];

    auto g_y_0_z_0_0_xx_xx_yy = buffer_1010_sddd[1083];

    auto g_y_0_z_0_0_xx_xx_yz = buffer_1010_sddd[1084];

    auto g_y_0_z_0_0_xx_xx_zz = buffer_1010_sddd[1085];

    auto g_y_0_z_0_0_xx_xy_xx = buffer_1010_sddd[1086];

    auto g_y_0_z_0_0_xx_xy_xy = buffer_1010_sddd[1087];

    auto g_y_0_z_0_0_xx_xy_xz = buffer_1010_sddd[1088];

    auto g_y_0_z_0_0_xx_xy_yy = buffer_1010_sddd[1089];

    auto g_y_0_z_0_0_xx_xy_yz = buffer_1010_sddd[1090];

    auto g_y_0_z_0_0_xx_xy_zz = buffer_1010_sddd[1091];

    auto g_y_0_z_0_0_xx_xz_xx = buffer_1010_sddd[1092];

    auto g_y_0_z_0_0_xx_xz_xy = buffer_1010_sddd[1093];

    auto g_y_0_z_0_0_xx_xz_xz = buffer_1010_sddd[1094];

    auto g_y_0_z_0_0_xx_xz_yy = buffer_1010_sddd[1095];

    auto g_y_0_z_0_0_xx_xz_yz = buffer_1010_sddd[1096];

    auto g_y_0_z_0_0_xx_xz_zz = buffer_1010_sddd[1097];

    auto g_y_0_z_0_0_xx_yy_xx = buffer_1010_sddd[1098];

    auto g_y_0_z_0_0_xx_yy_xy = buffer_1010_sddd[1099];

    auto g_y_0_z_0_0_xx_yy_xz = buffer_1010_sddd[1100];

    auto g_y_0_z_0_0_xx_yy_yy = buffer_1010_sddd[1101];

    auto g_y_0_z_0_0_xx_yy_yz = buffer_1010_sddd[1102];

    auto g_y_0_z_0_0_xx_yy_zz = buffer_1010_sddd[1103];

    auto g_y_0_z_0_0_xx_yz_xx = buffer_1010_sddd[1104];

    auto g_y_0_z_0_0_xx_yz_xy = buffer_1010_sddd[1105];

    auto g_y_0_z_0_0_xx_yz_xz = buffer_1010_sddd[1106];

    auto g_y_0_z_0_0_xx_yz_yy = buffer_1010_sddd[1107];

    auto g_y_0_z_0_0_xx_yz_yz = buffer_1010_sddd[1108];

    auto g_y_0_z_0_0_xx_yz_zz = buffer_1010_sddd[1109];

    auto g_y_0_z_0_0_xx_zz_xx = buffer_1010_sddd[1110];

    auto g_y_0_z_0_0_xx_zz_xy = buffer_1010_sddd[1111];

    auto g_y_0_z_0_0_xx_zz_xz = buffer_1010_sddd[1112];

    auto g_y_0_z_0_0_xx_zz_yy = buffer_1010_sddd[1113];

    auto g_y_0_z_0_0_xx_zz_yz = buffer_1010_sddd[1114];

    auto g_y_0_z_0_0_xx_zz_zz = buffer_1010_sddd[1115];

    auto g_y_0_z_0_0_xy_xx_xx = buffer_1010_sddd[1116];

    auto g_y_0_z_0_0_xy_xx_xy = buffer_1010_sddd[1117];

    auto g_y_0_z_0_0_xy_xx_xz = buffer_1010_sddd[1118];

    auto g_y_0_z_0_0_xy_xx_yy = buffer_1010_sddd[1119];

    auto g_y_0_z_0_0_xy_xx_yz = buffer_1010_sddd[1120];

    auto g_y_0_z_0_0_xy_xx_zz = buffer_1010_sddd[1121];

    auto g_y_0_z_0_0_xy_xy_xx = buffer_1010_sddd[1122];

    auto g_y_0_z_0_0_xy_xy_xy = buffer_1010_sddd[1123];

    auto g_y_0_z_0_0_xy_xy_xz = buffer_1010_sddd[1124];

    auto g_y_0_z_0_0_xy_xy_yy = buffer_1010_sddd[1125];

    auto g_y_0_z_0_0_xy_xy_yz = buffer_1010_sddd[1126];

    auto g_y_0_z_0_0_xy_xy_zz = buffer_1010_sddd[1127];

    auto g_y_0_z_0_0_xy_xz_xx = buffer_1010_sddd[1128];

    auto g_y_0_z_0_0_xy_xz_xy = buffer_1010_sddd[1129];

    auto g_y_0_z_0_0_xy_xz_xz = buffer_1010_sddd[1130];

    auto g_y_0_z_0_0_xy_xz_yy = buffer_1010_sddd[1131];

    auto g_y_0_z_0_0_xy_xz_yz = buffer_1010_sddd[1132];

    auto g_y_0_z_0_0_xy_xz_zz = buffer_1010_sddd[1133];

    auto g_y_0_z_0_0_xy_yy_xx = buffer_1010_sddd[1134];

    auto g_y_0_z_0_0_xy_yy_xy = buffer_1010_sddd[1135];

    auto g_y_0_z_0_0_xy_yy_xz = buffer_1010_sddd[1136];

    auto g_y_0_z_0_0_xy_yy_yy = buffer_1010_sddd[1137];

    auto g_y_0_z_0_0_xy_yy_yz = buffer_1010_sddd[1138];

    auto g_y_0_z_0_0_xy_yy_zz = buffer_1010_sddd[1139];

    auto g_y_0_z_0_0_xy_yz_xx = buffer_1010_sddd[1140];

    auto g_y_0_z_0_0_xy_yz_xy = buffer_1010_sddd[1141];

    auto g_y_0_z_0_0_xy_yz_xz = buffer_1010_sddd[1142];

    auto g_y_0_z_0_0_xy_yz_yy = buffer_1010_sddd[1143];

    auto g_y_0_z_0_0_xy_yz_yz = buffer_1010_sddd[1144];

    auto g_y_0_z_0_0_xy_yz_zz = buffer_1010_sddd[1145];

    auto g_y_0_z_0_0_xy_zz_xx = buffer_1010_sddd[1146];

    auto g_y_0_z_0_0_xy_zz_xy = buffer_1010_sddd[1147];

    auto g_y_0_z_0_0_xy_zz_xz = buffer_1010_sddd[1148];

    auto g_y_0_z_0_0_xy_zz_yy = buffer_1010_sddd[1149];

    auto g_y_0_z_0_0_xy_zz_yz = buffer_1010_sddd[1150];

    auto g_y_0_z_0_0_xy_zz_zz = buffer_1010_sddd[1151];

    auto g_y_0_z_0_0_xz_xx_xx = buffer_1010_sddd[1152];

    auto g_y_0_z_0_0_xz_xx_xy = buffer_1010_sddd[1153];

    auto g_y_0_z_0_0_xz_xx_xz = buffer_1010_sddd[1154];

    auto g_y_0_z_0_0_xz_xx_yy = buffer_1010_sddd[1155];

    auto g_y_0_z_0_0_xz_xx_yz = buffer_1010_sddd[1156];

    auto g_y_0_z_0_0_xz_xx_zz = buffer_1010_sddd[1157];

    auto g_y_0_z_0_0_xz_xy_xx = buffer_1010_sddd[1158];

    auto g_y_0_z_0_0_xz_xy_xy = buffer_1010_sddd[1159];

    auto g_y_0_z_0_0_xz_xy_xz = buffer_1010_sddd[1160];

    auto g_y_0_z_0_0_xz_xy_yy = buffer_1010_sddd[1161];

    auto g_y_0_z_0_0_xz_xy_yz = buffer_1010_sddd[1162];

    auto g_y_0_z_0_0_xz_xy_zz = buffer_1010_sddd[1163];

    auto g_y_0_z_0_0_xz_xz_xx = buffer_1010_sddd[1164];

    auto g_y_0_z_0_0_xz_xz_xy = buffer_1010_sddd[1165];

    auto g_y_0_z_0_0_xz_xz_xz = buffer_1010_sddd[1166];

    auto g_y_0_z_0_0_xz_xz_yy = buffer_1010_sddd[1167];

    auto g_y_0_z_0_0_xz_xz_yz = buffer_1010_sddd[1168];

    auto g_y_0_z_0_0_xz_xz_zz = buffer_1010_sddd[1169];

    auto g_y_0_z_0_0_xz_yy_xx = buffer_1010_sddd[1170];

    auto g_y_0_z_0_0_xz_yy_xy = buffer_1010_sddd[1171];

    auto g_y_0_z_0_0_xz_yy_xz = buffer_1010_sddd[1172];

    auto g_y_0_z_0_0_xz_yy_yy = buffer_1010_sddd[1173];

    auto g_y_0_z_0_0_xz_yy_yz = buffer_1010_sddd[1174];

    auto g_y_0_z_0_0_xz_yy_zz = buffer_1010_sddd[1175];

    auto g_y_0_z_0_0_xz_yz_xx = buffer_1010_sddd[1176];

    auto g_y_0_z_0_0_xz_yz_xy = buffer_1010_sddd[1177];

    auto g_y_0_z_0_0_xz_yz_xz = buffer_1010_sddd[1178];

    auto g_y_0_z_0_0_xz_yz_yy = buffer_1010_sddd[1179];

    auto g_y_0_z_0_0_xz_yz_yz = buffer_1010_sddd[1180];

    auto g_y_0_z_0_0_xz_yz_zz = buffer_1010_sddd[1181];

    auto g_y_0_z_0_0_xz_zz_xx = buffer_1010_sddd[1182];

    auto g_y_0_z_0_0_xz_zz_xy = buffer_1010_sddd[1183];

    auto g_y_0_z_0_0_xz_zz_xz = buffer_1010_sddd[1184];

    auto g_y_0_z_0_0_xz_zz_yy = buffer_1010_sddd[1185];

    auto g_y_0_z_0_0_xz_zz_yz = buffer_1010_sddd[1186];

    auto g_y_0_z_0_0_xz_zz_zz = buffer_1010_sddd[1187];

    auto g_y_0_z_0_0_yy_xx_xx = buffer_1010_sddd[1188];

    auto g_y_0_z_0_0_yy_xx_xy = buffer_1010_sddd[1189];

    auto g_y_0_z_0_0_yy_xx_xz = buffer_1010_sddd[1190];

    auto g_y_0_z_0_0_yy_xx_yy = buffer_1010_sddd[1191];

    auto g_y_0_z_0_0_yy_xx_yz = buffer_1010_sddd[1192];

    auto g_y_0_z_0_0_yy_xx_zz = buffer_1010_sddd[1193];

    auto g_y_0_z_0_0_yy_xy_xx = buffer_1010_sddd[1194];

    auto g_y_0_z_0_0_yy_xy_xy = buffer_1010_sddd[1195];

    auto g_y_0_z_0_0_yy_xy_xz = buffer_1010_sddd[1196];

    auto g_y_0_z_0_0_yy_xy_yy = buffer_1010_sddd[1197];

    auto g_y_0_z_0_0_yy_xy_yz = buffer_1010_sddd[1198];

    auto g_y_0_z_0_0_yy_xy_zz = buffer_1010_sddd[1199];

    auto g_y_0_z_0_0_yy_xz_xx = buffer_1010_sddd[1200];

    auto g_y_0_z_0_0_yy_xz_xy = buffer_1010_sddd[1201];

    auto g_y_0_z_0_0_yy_xz_xz = buffer_1010_sddd[1202];

    auto g_y_0_z_0_0_yy_xz_yy = buffer_1010_sddd[1203];

    auto g_y_0_z_0_0_yy_xz_yz = buffer_1010_sddd[1204];

    auto g_y_0_z_0_0_yy_xz_zz = buffer_1010_sddd[1205];

    auto g_y_0_z_0_0_yy_yy_xx = buffer_1010_sddd[1206];

    auto g_y_0_z_0_0_yy_yy_xy = buffer_1010_sddd[1207];

    auto g_y_0_z_0_0_yy_yy_xz = buffer_1010_sddd[1208];

    auto g_y_0_z_0_0_yy_yy_yy = buffer_1010_sddd[1209];

    auto g_y_0_z_0_0_yy_yy_yz = buffer_1010_sddd[1210];

    auto g_y_0_z_0_0_yy_yy_zz = buffer_1010_sddd[1211];

    auto g_y_0_z_0_0_yy_yz_xx = buffer_1010_sddd[1212];

    auto g_y_0_z_0_0_yy_yz_xy = buffer_1010_sddd[1213];

    auto g_y_0_z_0_0_yy_yz_xz = buffer_1010_sddd[1214];

    auto g_y_0_z_0_0_yy_yz_yy = buffer_1010_sddd[1215];

    auto g_y_0_z_0_0_yy_yz_yz = buffer_1010_sddd[1216];

    auto g_y_0_z_0_0_yy_yz_zz = buffer_1010_sddd[1217];

    auto g_y_0_z_0_0_yy_zz_xx = buffer_1010_sddd[1218];

    auto g_y_0_z_0_0_yy_zz_xy = buffer_1010_sddd[1219];

    auto g_y_0_z_0_0_yy_zz_xz = buffer_1010_sddd[1220];

    auto g_y_0_z_0_0_yy_zz_yy = buffer_1010_sddd[1221];

    auto g_y_0_z_0_0_yy_zz_yz = buffer_1010_sddd[1222];

    auto g_y_0_z_0_0_yy_zz_zz = buffer_1010_sddd[1223];

    auto g_y_0_z_0_0_yz_xx_xx = buffer_1010_sddd[1224];

    auto g_y_0_z_0_0_yz_xx_xy = buffer_1010_sddd[1225];

    auto g_y_0_z_0_0_yz_xx_xz = buffer_1010_sddd[1226];

    auto g_y_0_z_0_0_yz_xx_yy = buffer_1010_sddd[1227];

    auto g_y_0_z_0_0_yz_xx_yz = buffer_1010_sddd[1228];

    auto g_y_0_z_0_0_yz_xx_zz = buffer_1010_sddd[1229];

    auto g_y_0_z_0_0_yz_xy_xx = buffer_1010_sddd[1230];

    auto g_y_0_z_0_0_yz_xy_xy = buffer_1010_sddd[1231];

    auto g_y_0_z_0_0_yz_xy_xz = buffer_1010_sddd[1232];

    auto g_y_0_z_0_0_yz_xy_yy = buffer_1010_sddd[1233];

    auto g_y_0_z_0_0_yz_xy_yz = buffer_1010_sddd[1234];

    auto g_y_0_z_0_0_yz_xy_zz = buffer_1010_sddd[1235];

    auto g_y_0_z_0_0_yz_xz_xx = buffer_1010_sddd[1236];

    auto g_y_0_z_0_0_yz_xz_xy = buffer_1010_sddd[1237];

    auto g_y_0_z_0_0_yz_xz_xz = buffer_1010_sddd[1238];

    auto g_y_0_z_0_0_yz_xz_yy = buffer_1010_sddd[1239];

    auto g_y_0_z_0_0_yz_xz_yz = buffer_1010_sddd[1240];

    auto g_y_0_z_0_0_yz_xz_zz = buffer_1010_sddd[1241];

    auto g_y_0_z_0_0_yz_yy_xx = buffer_1010_sddd[1242];

    auto g_y_0_z_0_0_yz_yy_xy = buffer_1010_sddd[1243];

    auto g_y_0_z_0_0_yz_yy_xz = buffer_1010_sddd[1244];

    auto g_y_0_z_0_0_yz_yy_yy = buffer_1010_sddd[1245];

    auto g_y_0_z_0_0_yz_yy_yz = buffer_1010_sddd[1246];

    auto g_y_0_z_0_0_yz_yy_zz = buffer_1010_sddd[1247];

    auto g_y_0_z_0_0_yz_yz_xx = buffer_1010_sddd[1248];

    auto g_y_0_z_0_0_yz_yz_xy = buffer_1010_sddd[1249];

    auto g_y_0_z_0_0_yz_yz_xz = buffer_1010_sddd[1250];

    auto g_y_0_z_0_0_yz_yz_yy = buffer_1010_sddd[1251];

    auto g_y_0_z_0_0_yz_yz_yz = buffer_1010_sddd[1252];

    auto g_y_0_z_0_0_yz_yz_zz = buffer_1010_sddd[1253];

    auto g_y_0_z_0_0_yz_zz_xx = buffer_1010_sddd[1254];

    auto g_y_0_z_0_0_yz_zz_xy = buffer_1010_sddd[1255];

    auto g_y_0_z_0_0_yz_zz_xz = buffer_1010_sddd[1256];

    auto g_y_0_z_0_0_yz_zz_yy = buffer_1010_sddd[1257];

    auto g_y_0_z_0_0_yz_zz_yz = buffer_1010_sddd[1258];

    auto g_y_0_z_0_0_yz_zz_zz = buffer_1010_sddd[1259];

    auto g_y_0_z_0_0_zz_xx_xx = buffer_1010_sddd[1260];

    auto g_y_0_z_0_0_zz_xx_xy = buffer_1010_sddd[1261];

    auto g_y_0_z_0_0_zz_xx_xz = buffer_1010_sddd[1262];

    auto g_y_0_z_0_0_zz_xx_yy = buffer_1010_sddd[1263];

    auto g_y_0_z_0_0_zz_xx_yz = buffer_1010_sddd[1264];

    auto g_y_0_z_0_0_zz_xx_zz = buffer_1010_sddd[1265];

    auto g_y_0_z_0_0_zz_xy_xx = buffer_1010_sddd[1266];

    auto g_y_0_z_0_0_zz_xy_xy = buffer_1010_sddd[1267];

    auto g_y_0_z_0_0_zz_xy_xz = buffer_1010_sddd[1268];

    auto g_y_0_z_0_0_zz_xy_yy = buffer_1010_sddd[1269];

    auto g_y_0_z_0_0_zz_xy_yz = buffer_1010_sddd[1270];

    auto g_y_0_z_0_0_zz_xy_zz = buffer_1010_sddd[1271];

    auto g_y_0_z_0_0_zz_xz_xx = buffer_1010_sddd[1272];

    auto g_y_0_z_0_0_zz_xz_xy = buffer_1010_sddd[1273];

    auto g_y_0_z_0_0_zz_xz_xz = buffer_1010_sddd[1274];

    auto g_y_0_z_0_0_zz_xz_yy = buffer_1010_sddd[1275];

    auto g_y_0_z_0_0_zz_xz_yz = buffer_1010_sddd[1276];

    auto g_y_0_z_0_0_zz_xz_zz = buffer_1010_sddd[1277];

    auto g_y_0_z_0_0_zz_yy_xx = buffer_1010_sddd[1278];

    auto g_y_0_z_0_0_zz_yy_xy = buffer_1010_sddd[1279];

    auto g_y_0_z_0_0_zz_yy_xz = buffer_1010_sddd[1280];

    auto g_y_0_z_0_0_zz_yy_yy = buffer_1010_sddd[1281];

    auto g_y_0_z_0_0_zz_yy_yz = buffer_1010_sddd[1282];

    auto g_y_0_z_0_0_zz_yy_zz = buffer_1010_sddd[1283];

    auto g_y_0_z_0_0_zz_yz_xx = buffer_1010_sddd[1284];

    auto g_y_0_z_0_0_zz_yz_xy = buffer_1010_sddd[1285];

    auto g_y_0_z_0_0_zz_yz_xz = buffer_1010_sddd[1286];

    auto g_y_0_z_0_0_zz_yz_yy = buffer_1010_sddd[1287];

    auto g_y_0_z_0_0_zz_yz_yz = buffer_1010_sddd[1288];

    auto g_y_0_z_0_0_zz_yz_zz = buffer_1010_sddd[1289];

    auto g_y_0_z_0_0_zz_zz_xx = buffer_1010_sddd[1290];

    auto g_y_0_z_0_0_zz_zz_xy = buffer_1010_sddd[1291];

    auto g_y_0_z_0_0_zz_zz_xz = buffer_1010_sddd[1292];

    auto g_y_0_z_0_0_zz_zz_yy = buffer_1010_sddd[1293];

    auto g_y_0_z_0_0_zz_zz_yz = buffer_1010_sddd[1294];

    auto g_y_0_z_0_0_zz_zz_zz = buffer_1010_sddd[1295];

    auto g_z_0_x_0_0_xx_xx_xx = buffer_1010_sddd[1296];

    auto g_z_0_x_0_0_xx_xx_xy = buffer_1010_sddd[1297];

    auto g_z_0_x_0_0_xx_xx_xz = buffer_1010_sddd[1298];

    auto g_z_0_x_0_0_xx_xx_yy = buffer_1010_sddd[1299];

    auto g_z_0_x_0_0_xx_xx_yz = buffer_1010_sddd[1300];

    auto g_z_0_x_0_0_xx_xx_zz = buffer_1010_sddd[1301];

    auto g_z_0_x_0_0_xx_xy_xx = buffer_1010_sddd[1302];

    auto g_z_0_x_0_0_xx_xy_xy = buffer_1010_sddd[1303];

    auto g_z_0_x_0_0_xx_xy_xz = buffer_1010_sddd[1304];

    auto g_z_0_x_0_0_xx_xy_yy = buffer_1010_sddd[1305];

    auto g_z_0_x_0_0_xx_xy_yz = buffer_1010_sddd[1306];

    auto g_z_0_x_0_0_xx_xy_zz = buffer_1010_sddd[1307];

    auto g_z_0_x_0_0_xx_xz_xx = buffer_1010_sddd[1308];

    auto g_z_0_x_0_0_xx_xz_xy = buffer_1010_sddd[1309];

    auto g_z_0_x_0_0_xx_xz_xz = buffer_1010_sddd[1310];

    auto g_z_0_x_0_0_xx_xz_yy = buffer_1010_sddd[1311];

    auto g_z_0_x_0_0_xx_xz_yz = buffer_1010_sddd[1312];

    auto g_z_0_x_0_0_xx_xz_zz = buffer_1010_sddd[1313];

    auto g_z_0_x_0_0_xx_yy_xx = buffer_1010_sddd[1314];

    auto g_z_0_x_0_0_xx_yy_xy = buffer_1010_sddd[1315];

    auto g_z_0_x_0_0_xx_yy_xz = buffer_1010_sddd[1316];

    auto g_z_0_x_0_0_xx_yy_yy = buffer_1010_sddd[1317];

    auto g_z_0_x_0_0_xx_yy_yz = buffer_1010_sddd[1318];

    auto g_z_0_x_0_0_xx_yy_zz = buffer_1010_sddd[1319];

    auto g_z_0_x_0_0_xx_yz_xx = buffer_1010_sddd[1320];

    auto g_z_0_x_0_0_xx_yz_xy = buffer_1010_sddd[1321];

    auto g_z_0_x_0_0_xx_yz_xz = buffer_1010_sddd[1322];

    auto g_z_0_x_0_0_xx_yz_yy = buffer_1010_sddd[1323];

    auto g_z_0_x_0_0_xx_yz_yz = buffer_1010_sddd[1324];

    auto g_z_0_x_0_0_xx_yz_zz = buffer_1010_sddd[1325];

    auto g_z_0_x_0_0_xx_zz_xx = buffer_1010_sddd[1326];

    auto g_z_0_x_0_0_xx_zz_xy = buffer_1010_sddd[1327];

    auto g_z_0_x_0_0_xx_zz_xz = buffer_1010_sddd[1328];

    auto g_z_0_x_0_0_xx_zz_yy = buffer_1010_sddd[1329];

    auto g_z_0_x_0_0_xx_zz_yz = buffer_1010_sddd[1330];

    auto g_z_0_x_0_0_xx_zz_zz = buffer_1010_sddd[1331];

    auto g_z_0_x_0_0_xy_xx_xx = buffer_1010_sddd[1332];

    auto g_z_0_x_0_0_xy_xx_xy = buffer_1010_sddd[1333];

    auto g_z_0_x_0_0_xy_xx_xz = buffer_1010_sddd[1334];

    auto g_z_0_x_0_0_xy_xx_yy = buffer_1010_sddd[1335];

    auto g_z_0_x_0_0_xy_xx_yz = buffer_1010_sddd[1336];

    auto g_z_0_x_0_0_xy_xx_zz = buffer_1010_sddd[1337];

    auto g_z_0_x_0_0_xy_xy_xx = buffer_1010_sddd[1338];

    auto g_z_0_x_0_0_xy_xy_xy = buffer_1010_sddd[1339];

    auto g_z_0_x_0_0_xy_xy_xz = buffer_1010_sddd[1340];

    auto g_z_0_x_0_0_xy_xy_yy = buffer_1010_sddd[1341];

    auto g_z_0_x_0_0_xy_xy_yz = buffer_1010_sddd[1342];

    auto g_z_0_x_0_0_xy_xy_zz = buffer_1010_sddd[1343];

    auto g_z_0_x_0_0_xy_xz_xx = buffer_1010_sddd[1344];

    auto g_z_0_x_0_0_xy_xz_xy = buffer_1010_sddd[1345];

    auto g_z_0_x_0_0_xy_xz_xz = buffer_1010_sddd[1346];

    auto g_z_0_x_0_0_xy_xz_yy = buffer_1010_sddd[1347];

    auto g_z_0_x_0_0_xy_xz_yz = buffer_1010_sddd[1348];

    auto g_z_0_x_0_0_xy_xz_zz = buffer_1010_sddd[1349];

    auto g_z_0_x_0_0_xy_yy_xx = buffer_1010_sddd[1350];

    auto g_z_0_x_0_0_xy_yy_xy = buffer_1010_sddd[1351];

    auto g_z_0_x_0_0_xy_yy_xz = buffer_1010_sddd[1352];

    auto g_z_0_x_0_0_xy_yy_yy = buffer_1010_sddd[1353];

    auto g_z_0_x_0_0_xy_yy_yz = buffer_1010_sddd[1354];

    auto g_z_0_x_0_0_xy_yy_zz = buffer_1010_sddd[1355];

    auto g_z_0_x_0_0_xy_yz_xx = buffer_1010_sddd[1356];

    auto g_z_0_x_0_0_xy_yz_xy = buffer_1010_sddd[1357];

    auto g_z_0_x_0_0_xy_yz_xz = buffer_1010_sddd[1358];

    auto g_z_0_x_0_0_xy_yz_yy = buffer_1010_sddd[1359];

    auto g_z_0_x_0_0_xy_yz_yz = buffer_1010_sddd[1360];

    auto g_z_0_x_0_0_xy_yz_zz = buffer_1010_sddd[1361];

    auto g_z_0_x_0_0_xy_zz_xx = buffer_1010_sddd[1362];

    auto g_z_0_x_0_0_xy_zz_xy = buffer_1010_sddd[1363];

    auto g_z_0_x_0_0_xy_zz_xz = buffer_1010_sddd[1364];

    auto g_z_0_x_0_0_xy_zz_yy = buffer_1010_sddd[1365];

    auto g_z_0_x_0_0_xy_zz_yz = buffer_1010_sddd[1366];

    auto g_z_0_x_0_0_xy_zz_zz = buffer_1010_sddd[1367];

    auto g_z_0_x_0_0_xz_xx_xx = buffer_1010_sddd[1368];

    auto g_z_0_x_0_0_xz_xx_xy = buffer_1010_sddd[1369];

    auto g_z_0_x_0_0_xz_xx_xz = buffer_1010_sddd[1370];

    auto g_z_0_x_0_0_xz_xx_yy = buffer_1010_sddd[1371];

    auto g_z_0_x_0_0_xz_xx_yz = buffer_1010_sddd[1372];

    auto g_z_0_x_0_0_xz_xx_zz = buffer_1010_sddd[1373];

    auto g_z_0_x_0_0_xz_xy_xx = buffer_1010_sddd[1374];

    auto g_z_0_x_0_0_xz_xy_xy = buffer_1010_sddd[1375];

    auto g_z_0_x_0_0_xz_xy_xz = buffer_1010_sddd[1376];

    auto g_z_0_x_0_0_xz_xy_yy = buffer_1010_sddd[1377];

    auto g_z_0_x_0_0_xz_xy_yz = buffer_1010_sddd[1378];

    auto g_z_0_x_0_0_xz_xy_zz = buffer_1010_sddd[1379];

    auto g_z_0_x_0_0_xz_xz_xx = buffer_1010_sddd[1380];

    auto g_z_0_x_0_0_xz_xz_xy = buffer_1010_sddd[1381];

    auto g_z_0_x_0_0_xz_xz_xz = buffer_1010_sddd[1382];

    auto g_z_0_x_0_0_xz_xz_yy = buffer_1010_sddd[1383];

    auto g_z_0_x_0_0_xz_xz_yz = buffer_1010_sddd[1384];

    auto g_z_0_x_0_0_xz_xz_zz = buffer_1010_sddd[1385];

    auto g_z_0_x_0_0_xz_yy_xx = buffer_1010_sddd[1386];

    auto g_z_0_x_0_0_xz_yy_xy = buffer_1010_sddd[1387];

    auto g_z_0_x_0_0_xz_yy_xz = buffer_1010_sddd[1388];

    auto g_z_0_x_0_0_xz_yy_yy = buffer_1010_sddd[1389];

    auto g_z_0_x_0_0_xz_yy_yz = buffer_1010_sddd[1390];

    auto g_z_0_x_0_0_xz_yy_zz = buffer_1010_sddd[1391];

    auto g_z_0_x_0_0_xz_yz_xx = buffer_1010_sddd[1392];

    auto g_z_0_x_0_0_xz_yz_xy = buffer_1010_sddd[1393];

    auto g_z_0_x_0_0_xz_yz_xz = buffer_1010_sddd[1394];

    auto g_z_0_x_0_0_xz_yz_yy = buffer_1010_sddd[1395];

    auto g_z_0_x_0_0_xz_yz_yz = buffer_1010_sddd[1396];

    auto g_z_0_x_0_0_xz_yz_zz = buffer_1010_sddd[1397];

    auto g_z_0_x_0_0_xz_zz_xx = buffer_1010_sddd[1398];

    auto g_z_0_x_0_0_xz_zz_xy = buffer_1010_sddd[1399];

    auto g_z_0_x_0_0_xz_zz_xz = buffer_1010_sddd[1400];

    auto g_z_0_x_0_0_xz_zz_yy = buffer_1010_sddd[1401];

    auto g_z_0_x_0_0_xz_zz_yz = buffer_1010_sddd[1402];

    auto g_z_0_x_0_0_xz_zz_zz = buffer_1010_sddd[1403];

    auto g_z_0_x_0_0_yy_xx_xx = buffer_1010_sddd[1404];

    auto g_z_0_x_0_0_yy_xx_xy = buffer_1010_sddd[1405];

    auto g_z_0_x_0_0_yy_xx_xz = buffer_1010_sddd[1406];

    auto g_z_0_x_0_0_yy_xx_yy = buffer_1010_sddd[1407];

    auto g_z_0_x_0_0_yy_xx_yz = buffer_1010_sddd[1408];

    auto g_z_0_x_0_0_yy_xx_zz = buffer_1010_sddd[1409];

    auto g_z_0_x_0_0_yy_xy_xx = buffer_1010_sddd[1410];

    auto g_z_0_x_0_0_yy_xy_xy = buffer_1010_sddd[1411];

    auto g_z_0_x_0_0_yy_xy_xz = buffer_1010_sddd[1412];

    auto g_z_0_x_0_0_yy_xy_yy = buffer_1010_sddd[1413];

    auto g_z_0_x_0_0_yy_xy_yz = buffer_1010_sddd[1414];

    auto g_z_0_x_0_0_yy_xy_zz = buffer_1010_sddd[1415];

    auto g_z_0_x_0_0_yy_xz_xx = buffer_1010_sddd[1416];

    auto g_z_0_x_0_0_yy_xz_xy = buffer_1010_sddd[1417];

    auto g_z_0_x_0_0_yy_xz_xz = buffer_1010_sddd[1418];

    auto g_z_0_x_0_0_yy_xz_yy = buffer_1010_sddd[1419];

    auto g_z_0_x_0_0_yy_xz_yz = buffer_1010_sddd[1420];

    auto g_z_0_x_0_0_yy_xz_zz = buffer_1010_sddd[1421];

    auto g_z_0_x_0_0_yy_yy_xx = buffer_1010_sddd[1422];

    auto g_z_0_x_0_0_yy_yy_xy = buffer_1010_sddd[1423];

    auto g_z_0_x_0_0_yy_yy_xz = buffer_1010_sddd[1424];

    auto g_z_0_x_0_0_yy_yy_yy = buffer_1010_sddd[1425];

    auto g_z_0_x_0_0_yy_yy_yz = buffer_1010_sddd[1426];

    auto g_z_0_x_0_0_yy_yy_zz = buffer_1010_sddd[1427];

    auto g_z_0_x_0_0_yy_yz_xx = buffer_1010_sddd[1428];

    auto g_z_0_x_0_0_yy_yz_xy = buffer_1010_sddd[1429];

    auto g_z_0_x_0_0_yy_yz_xz = buffer_1010_sddd[1430];

    auto g_z_0_x_0_0_yy_yz_yy = buffer_1010_sddd[1431];

    auto g_z_0_x_0_0_yy_yz_yz = buffer_1010_sddd[1432];

    auto g_z_0_x_0_0_yy_yz_zz = buffer_1010_sddd[1433];

    auto g_z_0_x_0_0_yy_zz_xx = buffer_1010_sddd[1434];

    auto g_z_0_x_0_0_yy_zz_xy = buffer_1010_sddd[1435];

    auto g_z_0_x_0_0_yy_zz_xz = buffer_1010_sddd[1436];

    auto g_z_0_x_0_0_yy_zz_yy = buffer_1010_sddd[1437];

    auto g_z_0_x_0_0_yy_zz_yz = buffer_1010_sddd[1438];

    auto g_z_0_x_0_0_yy_zz_zz = buffer_1010_sddd[1439];

    auto g_z_0_x_0_0_yz_xx_xx = buffer_1010_sddd[1440];

    auto g_z_0_x_0_0_yz_xx_xy = buffer_1010_sddd[1441];

    auto g_z_0_x_0_0_yz_xx_xz = buffer_1010_sddd[1442];

    auto g_z_0_x_0_0_yz_xx_yy = buffer_1010_sddd[1443];

    auto g_z_0_x_0_0_yz_xx_yz = buffer_1010_sddd[1444];

    auto g_z_0_x_0_0_yz_xx_zz = buffer_1010_sddd[1445];

    auto g_z_0_x_0_0_yz_xy_xx = buffer_1010_sddd[1446];

    auto g_z_0_x_0_0_yz_xy_xy = buffer_1010_sddd[1447];

    auto g_z_0_x_0_0_yz_xy_xz = buffer_1010_sddd[1448];

    auto g_z_0_x_0_0_yz_xy_yy = buffer_1010_sddd[1449];

    auto g_z_0_x_0_0_yz_xy_yz = buffer_1010_sddd[1450];

    auto g_z_0_x_0_0_yz_xy_zz = buffer_1010_sddd[1451];

    auto g_z_0_x_0_0_yz_xz_xx = buffer_1010_sddd[1452];

    auto g_z_0_x_0_0_yz_xz_xy = buffer_1010_sddd[1453];

    auto g_z_0_x_0_0_yz_xz_xz = buffer_1010_sddd[1454];

    auto g_z_0_x_0_0_yz_xz_yy = buffer_1010_sddd[1455];

    auto g_z_0_x_0_0_yz_xz_yz = buffer_1010_sddd[1456];

    auto g_z_0_x_0_0_yz_xz_zz = buffer_1010_sddd[1457];

    auto g_z_0_x_0_0_yz_yy_xx = buffer_1010_sddd[1458];

    auto g_z_0_x_0_0_yz_yy_xy = buffer_1010_sddd[1459];

    auto g_z_0_x_0_0_yz_yy_xz = buffer_1010_sddd[1460];

    auto g_z_0_x_0_0_yz_yy_yy = buffer_1010_sddd[1461];

    auto g_z_0_x_0_0_yz_yy_yz = buffer_1010_sddd[1462];

    auto g_z_0_x_0_0_yz_yy_zz = buffer_1010_sddd[1463];

    auto g_z_0_x_0_0_yz_yz_xx = buffer_1010_sddd[1464];

    auto g_z_0_x_0_0_yz_yz_xy = buffer_1010_sddd[1465];

    auto g_z_0_x_0_0_yz_yz_xz = buffer_1010_sddd[1466];

    auto g_z_0_x_0_0_yz_yz_yy = buffer_1010_sddd[1467];

    auto g_z_0_x_0_0_yz_yz_yz = buffer_1010_sddd[1468];

    auto g_z_0_x_0_0_yz_yz_zz = buffer_1010_sddd[1469];

    auto g_z_0_x_0_0_yz_zz_xx = buffer_1010_sddd[1470];

    auto g_z_0_x_0_0_yz_zz_xy = buffer_1010_sddd[1471];

    auto g_z_0_x_0_0_yz_zz_xz = buffer_1010_sddd[1472];

    auto g_z_0_x_0_0_yz_zz_yy = buffer_1010_sddd[1473];

    auto g_z_0_x_0_0_yz_zz_yz = buffer_1010_sddd[1474];

    auto g_z_0_x_0_0_yz_zz_zz = buffer_1010_sddd[1475];

    auto g_z_0_x_0_0_zz_xx_xx = buffer_1010_sddd[1476];

    auto g_z_0_x_0_0_zz_xx_xy = buffer_1010_sddd[1477];

    auto g_z_0_x_0_0_zz_xx_xz = buffer_1010_sddd[1478];

    auto g_z_0_x_0_0_zz_xx_yy = buffer_1010_sddd[1479];

    auto g_z_0_x_0_0_zz_xx_yz = buffer_1010_sddd[1480];

    auto g_z_0_x_0_0_zz_xx_zz = buffer_1010_sddd[1481];

    auto g_z_0_x_0_0_zz_xy_xx = buffer_1010_sddd[1482];

    auto g_z_0_x_0_0_zz_xy_xy = buffer_1010_sddd[1483];

    auto g_z_0_x_0_0_zz_xy_xz = buffer_1010_sddd[1484];

    auto g_z_0_x_0_0_zz_xy_yy = buffer_1010_sddd[1485];

    auto g_z_0_x_0_0_zz_xy_yz = buffer_1010_sddd[1486];

    auto g_z_0_x_0_0_zz_xy_zz = buffer_1010_sddd[1487];

    auto g_z_0_x_0_0_zz_xz_xx = buffer_1010_sddd[1488];

    auto g_z_0_x_0_0_zz_xz_xy = buffer_1010_sddd[1489];

    auto g_z_0_x_0_0_zz_xz_xz = buffer_1010_sddd[1490];

    auto g_z_0_x_0_0_zz_xz_yy = buffer_1010_sddd[1491];

    auto g_z_0_x_0_0_zz_xz_yz = buffer_1010_sddd[1492];

    auto g_z_0_x_0_0_zz_xz_zz = buffer_1010_sddd[1493];

    auto g_z_0_x_0_0_zz_yy_xx = buffer_1010_sddd[1494];

    auto g_z_0_x_0_0_zz_yy_xy = buffer_1010_sddd[1495];

    auto g_z_0_x_0_0_zz_yy_xz = buffer_1010_sddd[1496];

    auto g_z_0_x_0_0_zz_yy_yy = buffer_1010_sddd[1497];

    auto g_z_0_x_0_0_zz_yy_yz = buffer_1010_sddd[1498];

    auto g_z_0_x_0_0_zz_yy_zz = buffer_1010_sddd[1499];

    auto g_z_0_x_0_0_zz_yz_xx = buffer_1010_sddd[1500];

    auto g_z_0_x_0_0_zz_yz_xy = buffer_1010_sddd[1501];

    auto g_z_0_x_0_0_zz_yz_xz = buffer_1010_sddd[1502];

    auto g_z_0_x_0_0_zz_yz_yy = buffer_1010_sddd[1503];

    auto g_z_0_x_0_0_zz_yz_yz = buffer_1010_sddd[1504];

    auto g_z_0_x_0_0_zz_yz_zz = buffer_1010_sddd[1505];

    auto g_z_0_x_0_0_zz_zz_xx = buffer_1010_sddd[1506];

    auto g_z_0_x_0_0_zz_zz_xy = buffer_1010_sddd[1507];

    auto g_z_0_x_0_0_zz_zz_xz = buffer_1010_sddd[1508];

    auto g_z_0_x_0_0_zz_zz_yy = buffer_1010_sddd[1509];

    auto g_z_0_x_0_0_zz_zz_yz = buffer_1010_sddd[1510];

    auto g_z_0_x_0_0_zz_zz_zz = buffer_1010_sddd[1511];

    auto g_z_0_y_0_0_xx_xx_xx = buffer_1010_sddd[1512];

    auto g_z_0_y_0_0_xx_xx_xy = buffer_1010_sddd[1513];

    auto g_z_0_y_0_0_xx_xx_xz = buffer_1010_sddd[1514];

    auto g_z_0_y_0_0_xx_xx_yy = buffer_1010_sddd[1515];

    auto g_z_0_y_0_0_xx_xx_yz = buffer_1010_sddd[1516];

    auto g_z_0_y_0_0_xx_xx_zz = buffer_1010_sddd[1517];

    auto g_z_0_y_0_0_xx_xy_xx = buffer_1010_sddd[1518];

    auto g_z_0_y_0_0_xx_xy_xy = buffer_1010_sddd[1519];

    auto g_z_0_y_0_0_xx_xy_xz = buffer_1010_sddd[1520];

    auto g_z_0_y_0_0_xx_xy_yy = buffer_1010_sddd[1521];

    auto g_z_0_y_0_0_xx_xy_yz = buffer_1010_sddd[1522];

    auto g_z_0_y_0_0_xx_xy_zz = buffer_1010_sddd[1523];

    auto g_z_0_y_0_0_xx_xz_xx = buffer_1010_sddd[1524];

    auto g_z_0_y_0_0_xx_xz_xy = buffer_1010_sddd[1525];

    auto g_z_0_y_0_0_xx_xz_xz = buffer_1010_sddd[1526];

    auto g_z_0_y_0_0_xx_xz_yy = buffer_1010_sddd[1527];

    auto g_z_0_y_0_0_xx_xz_yz = buffer_1010_sddd[1528];

    auto g_z_0_y_0_0_xx_xz_zz = buffer_1010_sddd[1529];

    auto g_z_0_y_0_0_xx_yy_xx = buffer_1010_sddd[1530];

    auto g_z_0_y_0_0_xx_yy_xy = buffer_1010_sddd[1531];

    auto g_z_0_y_0_0_xx_yy_xz = buffer_1010_sddd[1532];

    auto g_z_0_y_0_0_xx_yy_yy = buffer_1010_sddd[1533];

    auto g_z_0_y_0_0_xx_yy_yz = buffer_1010_sddd[1534];

    auto g_z_0_y_0_0_xx_yy_zz = buffer_1010_sddd[1535];

    auto g_z_0_y_0_0_xx_yz_xx = buffer_1010_sddd[1536];

    auto g_z_0_y_0_0_xx_yz_xy = buffer_1010_sddd[1537];

    auto g_z_0_y_0_0_xx_yz_xz = buffer_1010_sddd[1538];

    auto g_z_0_y_0_0_xx_yz_yy = buffer_1010_sddd[1539];

    auto g_z_0_y_0_0_xx_yz_yz = buffer_1010_sddd[1540];

    auto g_z_0_y_0_0_xx_yz_zz = buffer_1010_sddd[1541];

    auto g_z_0_y_0_0_xx_zz_xx = buffer_1010_sddd[1542];

    auto g_z_0_y_0_0_xx_zz_xy = buffer_1010_sddd[1543];

    auto g_z_0_y_0_0_xx_zz_xz = buffer_1010_sddd[1544];

    auto g_z_0_y_0_0_xx_zz_yy = buffer_1010_sddd[1545];

    auto g_z_0_y_0_0_xx_zz_yz = buffer_1010_sddd[1546];

    auto g_z_0_y_0_0_xx_zz_zz = buffer_1010_sddd[1547];

    auto g_z_0_y_0_0_xy_xx_xx = buffer_1010_sddd[1548];

    auto g_z_0_y_0_0_xy_xx_xy = buffer_1010_sddd[1549];

    auto g_z_0_y_0_0_xy_xx_xz = buffer_1010_sddd[1550];

    auto g_z_0_y_0_0_xy_xx_yy = buffer_1010_sddd[1551];

    auto g_z_0_y_0_0_xy_xx_yz = buffer_1010_sddd[1552];

    auto g_z_0_y_0_0_xy_xx_zz = buffer_1010_sddd[1553];

    auto g_z_0_y_0_0_xy_xy_xx = buffer_1010_sddd[1554];

    auto g_z_0_y_0_0_xy_xy_xy = buffer_1010_sddd[1555];

    auto g_z_0_y_0_0_xy_xy_xz = buffer_1010_sddd[1556];

    auto g_z_0_y_0_0_xy_xy_yy = buffer_1010_sddd[1557];

    auto g_z_0_y_0_0_xy_xy_yz = buffer_1010_sddd[1558];

    auto g_z_0_y_0_0_xy_xy_zz = buffer_1010_sddd[1559];

    auto g_z_0_y_0_0_xy_xz_xx = buffer_1010_sddd[1560];

    auto g_z_0_y_0_0_xy_xz_xy = buffer_1010_sddd[1561];

    auto g_z_0_y_0_0_xy_xz_xz = buffer_1010_sddd[1562];

    auto g_z_0_y_0_0_xy_xz_yy = buffer_1010_sddd[1563];

    auto g_z_0_y_0_0_xy_xz_yz = buffer_1010_sddd[1564];

    auto g_z_0_y_0_0_xy_xz_zz = buffer_1010_sddd[1565];

    auto g_z_0_y_0_0_xy_yy_xx = buffer_1010_sddd[1566];

    auto g_z_0_y_0_0_xy_yy_xy = buffer_1010_sddd[1567];

    auto g_z_0_y_0_0_xy_yy_xz = buffer_1010_sddd[1568];

    auto g_z_0_y_0_0_xy_yy_yy = buffer_1010_sddd[1569];

    auto g_z_0_y_0_0_xy_yy_yz = buffer_1010_sddd[1570];

    auto g_z_0_y_0_0_xy_yy_zz = buffer_1010_sddd[1571];

    auto g_z_0_y_0_0_xy_yz_xx = buffer_1010_sddd[1572];

    auto g_z_0_y_0_0_xy_yz_xy = buffer_1010_sddd[1573];

    auto g_z_0_y_0_0_xy_yz_xz = buffer_1010_sddd[1574];

    auto g_z_0_y_0_0_xy_yz_yy = buffer_1010_sddd[1575];

    auto g_z_0_y_0_0_xy_yz_yz = buffer_1010_sddd[1576];

    auto g_z_0_y_0_0_xy_yz_zz = buffer_1010_sddd[1577];

    auto g_z_0_y_0_0_xy_zz_xx = buffer_1010_sddd[1578];

    auto g_z_0_y_0_0_xy_zz_xy = buffer_1010_sddd[1579];

    auto g_z_0_y_0_0_xy_zz_xz = buffer_1010_sddd[1580];

    auto g_z_0_y_0_0_xy_zz_yy = buffer_1010_sddd[1581];

    auto g_z_0_y_0_0_xy_zz_yz = buffer_1010_sddd[1582];

    auto g_z_0_y_0_0_xy_zz_zz = buffer_1010_sddd[1583];

    auto g_z_0_y_0_0_xz_xx_xx = buffer_1010_sddd[1584];

    auto g_z_0_y_0_0_xz_xx_xy = buffer_1010_sddd[1585];

    auto g_z_0_y_0_0_xz_xx_xz = buffer_1010_sddd[1586];

    auto g_z_0_y_0_0_xz_xx_yy = buffer_1010_sddd[1587];

    auto g_z_0_y_0_0_xz_xx_yz = buffer_1010_sddd[1588];

    auto g_z_0_y_0_0_xz_xx_zz = buffer_1010_sddd[1589];

    auto g_z_0_y_0_0_xz_xy_xx = buffer_1010_sddd[1590];

    auto g_z_0_y_0_0_xz_xy_xy = buffer_1010_sddd[1591];

    auto g_z_0_y_0_0_xz_xy_xz = buffer_1010_sddd[1592];

    auto g_z_0_y_0_0_xz_xy_yy = buffer_1010_sddd[1593];

    auto g_z_0_y_0_0_xz_xy_yz = buffer_1010_sddd[1594];

    auto g_z_0_y_0_0_xz_xy_zz = buffer_1010_sddd[1595];

    auto g_z_0_y_0_0_xz_xz_xx = buffer_1010_sddd[1596];

    auto g_z_0_y_0_0_xz_xz_xy = buffer_1010_sddd[1597];

    auto g_z_0_y_0_0_xz_xz_xz = buffer_1010_sddd[1598];

    auto g_z_0_y_0_0_xz_xz_yy = buffer_1010_sddd[1599];

    auto g_z_0_y_0_0_xz_xz_yz = buffer_1010_sddd[1600];

    auto g_z_0_y_0_0_xz_xz_zz = buffer_1010_sddd[1601];

    auto g_z_0_y_0_0_xz_yy_xx = buffer_1010_sddd[1602];

    auto g_z_0_y_0_0_xz_yy_xy = buffer_1010_sddd[1603];

    auto g_z_0_y_0_0_xz_yy_xz = buffer_1010_sddd[1604];

    auto g_z_0_y_0_0_xz_yy_yy = buffer_1010_sddd[1605];

    auto g_z_0_y_0_0_xz_yy_yz = buffer_1010_sddd[1606];

    auto g_z_0_y_0_0_xz_yy_zz = buffer_1010_sddd[1607];

    auto g_z_0_y_0_0_xz_yz_xx = buffer_1010_sddd[1608];

    auto g_z_0_y_0_0_xz_yz_xy = buffer_1010_sddd[1609];

    auto g_z_0_y_0_0_xz_yz_xz = buffer_1010_sddd[1610];

    auto g_z_0_y_0_0_xz_yz_yy = buffer_1010_sddd[1611];

    auto g_z_0_y_0_0_xz_yz_yz = buffer_1010_sddd[1612];

    auto g_z_0_y_0_0_xz_yz_zz = buffer_1010_sddd[1613];

    auto g_z_0_y_0_0_xz_zz_xx = buffer_1010_sddd[1614];

    auto g_z_0_y_0_0_xz_zz_xy = buffer_1010_sddd[1615];

    auto g_z_0_y_0_0_xz_zz_xz = buffer_1010_sddd[1616];

    auto g_z_0_y_0_0_xz_zz_yy = buffer_1010_sddd[1617];

    auto g_z_0_y_0_0_xz_zz_yz = buffer_1010_sddd[1618];

    auto g_z_0_y_0_0_xz_zz_zz = buffer_1010_sddd[1619];

    auto g_z_0_y_0_0_yy_xx_xx = buffer_1010_sddd[1620];

    auto g_z_0_y_0_0_yy_xx_xy = buffer_1010_sddd[1621];

    auto g_z_0_y_0_0_yy_xx_xz = buffer_1010_sddd[1622];

    auto g_z_0_y_0_0_yy_xx_yy = buffer_1010_sddd[1623];

    auto g_z_0_y_0_0_yy_xx_yz = buffer_1010_sddd[1624];

    auto g_z_0_y_0_0_yy_xx_zz = buffer_1010_sddd[1625];

    auto g_z_0_y_0_0_yy_xy_xx = buffer_1010_sddd[1626];

    auto g_z_0_y_0_0_yy_xy_xy = buffer_1010_sddd[1627];

    auto g_z_0_y_0_0_yy_xy_xz = buffer_1010_sddd[1628];

    auto g_z_0_y_0_0_yy_xy_yy = buffer_1010_sddd[1629];

    auto g_z_0_y_0_0_yy_xy_yz = buffer_1010_sddd[1630];

    auto g_z_0_y_0_0_yy_xy_zz = buffer_1010_sddd[1631];

    auto g_z_0_y_0_0_yy_xz_xx = buffer_1010_sddd[1632];

    auto g_z_0_y_0_0_yy_xz_xy = buffer_1010_sddd[1633];

    auto g_z_0_y_0_0_yy_xz_xz = buffer_1010_sddd[1634];

    auto g_z_0_y_0_0_yy_xz_yy = buffer_1010_sddd[1635];

    auto g_z_0_y_0_0_yy_xz_yz = buffer_1010_sddd[1636];

    auto g_z_0_y_0_0_yy_xz_zz = buffer_1010_sddd[1637];

    auto g_z_0_y_0_0_yy_yy_xx = buffer_1010_sddd[1638];

    auto g_z_0_y_0_0_yy_yy_xy = buffer_1010_sddd[1639];

    auto g_z_0_y_0_0_yy_yy_xz = buffer_1010_sddd[1640];

    auto g_z_0_y_0_0_yy_yy_yy = buffer_1010_sddd[1641];

    auto g_z_0_y_0_0_yy_yy_yz = buffer_1010_sddd[1642];

    auto g_z_0_y_0_0_yy_yy_zz = buffer_1010_sddd[1643];

    auto g_z_0_y_0_0_yy_yz_xx = buffer_1010_sddd[1644];

    auto g_z_0_y_0_0_yy_yz_xy = buffer_1010_sddd[1645];

    auto g_z_0_y_0_0_yy_yz_xz = buffer_1010_sddd[1646];

    auto g_z_0_y_0_0_yy_yz_yy = buffer_1010_sddd[1647];

    auto g_z_0_y_0_0_yy_yz_yz = buffer_1010_sddd[1648];

    auto g_z_0_y_0_0_yy_yz_zz = buffer_1010_sddd[1649];

    auto g_z_0_y_0_0_yy_zz_xx = buffer_1010_sddd[1650];

    auto g_z_0_y_0_0_yy_zz_xy = buffer_1010_sddd[1651];

    auto g_z_0_y_0_0_yy_zz_xz = buffer_1010_sddd[1652];

    auto g_z_0_y_0_0_yy_zz_yy = buffer_1010_sddd[1653];

    auto g_z_0_y_0_0_yy_zz_yz = buffer_1010_sddd[1654];

    auto g_z_0_y_0_0_yy_zz_zz = buffer_1010_sddd[1655];

    auto g_z_0_y_0_0_yz_xx_xx = buffer_1010_sddd[1656];

    auto g_z_0_y_0_0_yz_xx_xy = buffer_1010_sddd[1657];

    auto g_z_0_y_0_0_yz_xx_xz = buffer_1010_sddd[1658];

    auto g_z_0_y_0_0_yz_xx_yy = buffer_1010_sddd[1659];

    auto g_z_0_y_0_0_yz_xx_yz = buffer_1010_sddd[1660];

    auto g_z_0_y_0_0_yz_xx_zz = buffer_1010_sddd[1661];

    auto g_z_0_y_0_0_yz_xy_xx = buffer_1010_sddd[1662];

    auto g_z_0_y_0_0_yz_xy_xy = buffer_1010_sddd[1663];

    auto g_z_0_y_0_0_yz_xy_xz = buffer_1010_sddd[1664];

    auto g_z_0_y_0_0_yz_xy_yy = buffer_1010_sddd[1665];

    auto g_z_0_y_0_0_yz_xy_yz = buffer_1010_sddd[1666];

    auto g_z_0_y_0_0_yz_xy_zz = buffer_1010_sddd[1667];

    auto g_z_0_y_0_0_yz_xz_xx = buffer_1010_sddd[1668];

    auto g_z_0_y_0_0_yz_xz_xy = buffer_1010_sddd[1669];

    auto g_z_0_y_0_0_yz_xz_xz = buffer_1010_sddd[1670];

    auto g_z_0_y_0_0_yz_xz_yy = buffer_1010_sddd[1671];

    auto g_z_0_y_0_0_yz_xz_yz = buffer_1010_sddd[1672];

    auto g_z_0_y_0_0_yz_xz_zz = buffer_1010_sddd[1673];

    auto g_z_0_y_0_0_yz_yy_xx = buffer_1010_sddd[1674];

    auto g_z_0_y_0_0_yz_yy_xy = buffer_1010_sddd[1675];

    auto g_z_0_y_0_0_yz_yy_xz = buffer_1010_sddd[1676];

    auto g_z_0_y_0_0_yz_yy_yy = buffer_1010_sddd[1677];

    auto g_z_0_y_0_0_yz_yy_yz = buffer_1010_sddd[1678];

    auto g_z_0_y_0_0_yz_yy_zz = buffer_1010_sddd[1679];

    auto g_z_0_y_0_0_yz_yz_xx = buffer_1010_sddd[1680];

    auto g_z_0_y_0_0_yz_yz_xy = buffer_1010_sddd[1681];

    auto g_z_0_y_0_0_yz_yz_xz = buffer_1010_sddd[1682];

    auto g_z_0_y_0_0_yz_yz_yy = buffer_1010_sddd[1683];

    auto g_z_0_y_0_0_yz_yz_yz = buffer_1010_sddd[1684];

    auto g_z_0_y_0_0_yz_yz_zz = buffer_1010_sddd[1685];

    auto g_z_0_y_0_0_yz_zz_xx = buffer_1010_sddd[1686];

    auto g_z_0_y_0_0_yz_zz_xy = buffer_1010_sddd[1687];

    auto g_z_0_y_0_0_yz_zz_xz = buffer_1010_sddd[1688];

    auto g_z_0_y_0_0_yz_zz_yy = buffer_1010_sddd[1689];

    auto g_z_0_y_0_0_yz_zz_yz = buffer_1010_sddd[1690];

    auto g_z_0_y_0_0_yz_zz_zz = buffer_1010_sddd[1691];

    auto g_z_0_y_0_0_zz_xx_xx = buffer_1010_sddd[1692];

    auto g_z_0_y_0_0_zz_xx_xy = buffer_1010_sddd[1693];

    auto g_z_0_y_0_0_zz_xx_xz = buffer_1010_sddd[1694];

    auto g_z_0_y_0_0_zz_xx_yy = buffer_1010_sddd[1695];

    auto g_z_0_y_0_0_zz_xx_yz = buffer_1010_sddd[1696];

    auto g_z_0_y_0_0_zz_xx_zz = buffer_1010_sddd[1697];

    auto g_z_0_y_0_0_zz_xy_xx = buffer_1010_sddd[1698];

    auto g_z_0_y_0_0_zz_xy_xy = buffer_1010_sddd[1699];

    auto g_z_0_y_0_0_zz_xy_xz = buffer_1010_sddd[1700];

    auto g_z_0_y_0_0_zz_xy_yy = buffer_1010_sddd[1701];

    auto g_z_0_y_0_0_zz_xy_yz = buffer_1010_sddd[1702];

    auto g_z_0_y_0_0_zz_xy_zz = buffer_1010_sddd[1703];

    auto g_z_0_y_0_0_zz_xz_xx = buffer_1010_sddd[1704];

    auto g_z_0_y_0_0_zz_xz_xy = buffer_1010_sddd[1705];

    auto g_z_0_y_0_0_zz_xz_xz = buffer_1010_sddd[1706];

    auto g_z_0_y_0_0_zz_xz_yy = buffer_1010_sddd[1707];

    auto g_z_0_y_0_0_zz_xz_yz = buffer_1010_sddd[1708];

    auto g_z_0_y_0_0_zz_xz_zz = buffer_1010_sddd[1709];

    auto g_z_0_y_0_0_zz_yy_xx = buffer_1010_sddd[1710];

    auto g_z_0_y_0_0_zz_yy_xy = buffer_1010_sddd[1711];

    auto g_z_0_y_0_0_zz_yy_xz = buffer_1010_sddd[1712];

    auto g_z_0_y_0_0_zz_yy_yy = buffer_1010_sddd[1713];

    auto g_z_0_y_0_0_zz_yy_yz = buffer_1010_sddd[1714];

    auto g_z_0_y_0_0_zz_yy_zz = buffer_1010_sddd[1715];

    auto g_z_0_y_0_0_zz_yz_xx = buffer_1010_sddd[1716];

    auto g_z_0_y_0_0_zz_yz_xy = buffer_1010_sddd[1717];

    auto g_z_0_y_0_0_zz_yz_xz = buffer_1010_sddd[1718];

    auto g_z_0_y_0_0_zz_yz_yy = buffer_1010_sddd[1719];

    auto g_z_0_y_0_0_zz_yz_yz = buffer_1010_sddd[1720];

    auto g_z_0_y_0_0_zz_yz_zz = buffer_1010_sddd[1721];

    auto g_z_0_y_0_0_zz_zz_xx = buffer_1010_sddd[1722];

    auto g_z_0_y_0_0_zz_zz_xy = buffer_1010_sddd[1723];

    auto g_z_0_y_0_0_zz_zz_xz = buffer_1010_sddd[1724];

    auto g_z_0_y_0_0_zz_zz_yy = buffer_1010_sddd[1725];

    auto g_z_0_y_0_0_zz_zz_yz = buffer_1010_sddd[1726];

    auto g_z_0_y_0_0_zz_zz_zz = buffer_1010_sddd[1727];

    auto g_z_0_z_0_0_xx_xx_xx = buffer_1010_sddd[1728];

    auto g_z_0_z_0_0_xx_xx_xy = buffer_1010_sddd[1729];

    auto g_z_0_z_0_0_xx_xx_xz = buffer_1010_sddd[1730];

    auto g_z_0_z_0_0_xx_xx_yy = buffer_1010_sddd[1731];

    auto g_z_0_z_0_0_xx_xx_yz = buffer_1010_sddd[1732];

    auto g_z_0_z_0_0_xx_xx_zz = buffer_1010_sddd[1733];

    auto g_z_0_z_0_0_xx_xy_xx = buffer_1010_sddd[1734];

    auto g_z_0_z_0_0_xx_xy_xy = buffer_1010_sddd[1735];

    auto g_z_0_z_0_0_xx_xy_xz = buffer_1010_sddd[1736];

    auto g_z_0_z_0_0_xx_xy_yy = buffer_1010_sddd[1737];

    auto g_z_0_z_0_0_xx_xy_yz = buffer_1010_sddd[1738];

    auto g_z_0_z_0_0_xx_xy_zz = buffer_1010_sddd[1739];

    auto g_z_0_z_0_0_xx_xz_xx = buffer_1010_sddd[1740];

    auto g_z_0_z_0_0_xx_xz_xy = buffer_1010_sddd[1741];

    auto g_z_0_z_0_0_xx_xz_xz = buffer_1010_sddd[1742];

    auto g_z_0_z_0_0_xx_xz_yy = buffer_1010_sddd[1743];

    auto g_z_0_z_0_0_xx_xz_yz = buffer_1010_sddd[1744];

    auto g_z_0_z_0_0_xx_xz_zz = buffer_1010_sddd[1745];

    auto g_z_0_z_0_0_xx_yy_xx = buffer_1010_sddd[1746];

    auto g_z_0_z_0_0_xx_yy_xy = buffer_1010_sddd[1747];

    auto g_z_0_z_0_0_xx_yy_xz = buffer_1010_sddd[1748];

    auto g_z_0_z_0_0_xx_yy_yy = buffer_1010_sddd[1749];

    auto g_z_0_z_0_0_xx_yy_yz = buffer_1010_sddd[1750];

    auto g_z_0_z_0_0_xx_yy_zz = buffer_1010_sddd[1751];

    auto g_z_0_z_0_0_xx_yz_xx = buffer_1010_sddd[1752];

    auto g_z_0_z_0_0_xx_yz_xy = buffer_1010_sddd[1753];

    auto g_z_0_z_0_0_xx_yz_xz = buffer_1010_sddd[1754];

    auto g_z_0_z_0_0_xx_yz_yy = buffer_1010_sddd[1755];

    auto g_z_0_z_0_0_xx_yz_yz = buffer_1010_sddd[1756];

    auto g_z_0_z_0_0_xx_yz_zz = buffer_1010_sddd[1757];

    auto g_z_0_z_0_0_xx_zz_xx = buffer_1010_sddd[1758];

    auto g_z_0_z_0_0_xx_zz_xy = buffer_1010_sddd[1759];

    auto g_z_0_z_0_0_xx_zz_xz = buffer_1010_sddd[1760];

    auto g_z_0_z_0_0_xx_zz_yy = buffer_1010_sddd[1761];

    auto g_z_0_z_0_0_xx_zz_yz = buffer_1010_sddd[1762];

    auto g_z_0_z_0_0_xx_zz_zz = buffer_1010_sddd[1763];

    auto g_z_0_z_0_0_xy_xx_xx = buffer_1010_sddd[1764];

    auto g_z_0_z_0_0_xy_xx_xy = buffer_1010_sddd[1765];

    auto g_z_0_z_0_0_xy_xx_xz = buffer_1010_sddd[1766];

    auto g_z_0_z_0_0_xy_xx_yy = buffer_1010_sddd[1767];

    auto g_z_0_z_0_0_xy_xx_yz = buffer_1010_sddd[1768];

    auto g_z_0_z_0_0_xy_xx_zz = buffer_1010_sddd[1769];

    auto g_z_0_z_0_0_xy_xy_xx = buffer_1010_sddd[1770];

    auto g_z_0_z_0_0_xy_xy_xy = buffer_1010_sddd[1771];

    auto g_z_0_z_0_0_xy_xy_xz = buffer_1010_sddd[1772];

    auto g_z_0_z_0_0_xy_xy_yy = buffer_1010_sddd[1773];

    auto g_z_0_z_0_0_xy_xy_yz = buffer_1010_sddd[1774];

    auto g_z_0_z_0_0_xy_xy_zz = buffer_1010_sddd[1775];

    auto g_z_0_z_0_0_xy_xz_xx = buffer_1010_sddd[1776];

    auto g_z_0_z_0_0_xy_xz_xy = buffer_1010_sddd[1777];

    auto g_z_0_z_0_0_xy_xz_xz = buffer_1010_sddd[1778];

    auto g_z_0_z_0_0_xy_xz_yy = buffer_1010_sddd[1779];

    auto g_z_0_z_0_0_xy_xz_yz = buffer_1010_sddd[1780];

    auto g_z_0_z_0_0_xy_xz_zz = buffer_1010_sddd[1781];

    auto g_z_0_z_0_0_xy_yy_xx = buffer_1010_sddd[1782];

    auto g_z_0_z_0_0_xy_yy_xy = buffer_1010_sddd[1783];

    auto g_z_0_z_0_0_xy_yy_xz = buffer_1010_sddd[1784];

    auto g_z_0_z_0_0_xy_yy_yy = buffer_1010_sddd[1785];

    auto g_z_0_z_0_0_xy_yy_yz = buffer_1010_sddd[1786];

    auto g_z_0_z_0_0_xy_yy_zz = buffer_1010_sddd[1787];

    auto g_z_0_z_0_0_xy_yz_xx = buffer_1010_sddd[1788];

    auto g_z_0_z_0_0_xy_yz_xy = buffer_1010_sddd[1789];

    auto g_z_0_z_0_0_xy_yz_xz = buffer_1010_sddd[1790];

    auto g_z_0_z_0_0_xy_yz_yy = buffer_1010_sddd[1791];

    auto g_z_0_z_0_0_xy_yz_yz = buffer_1010_sddd[1792];

    auto g_z_0_z_0_0_xy_yz_zz = buffer_1010_sddd[1793];

    auto g_z_0_z_0_0_xy_zz_xx = buffer_1010_sddd[1794];

    auto g_z_0_z_0_0_xy_zz_xy = buffer_1010_sddd[1795];

    auto g_z_0_z_0_0_xy_zz_xz = buffer_1010_sddd[1796];

    auto g_z_0_z_0_0_xy_zz_yy = buffer_1010_sddd[1797];

    auto g_z_0_z_0_0_xy_zz_yz = buffer_1010_sddd[1798];

    auto g_z_0_z_0_0_xy_zz_zz = buffer_1010_sddd[1799];

    auto g_z_0_z_0_0_xz_xx_xx = buffer_1010_sddd[1800];

    auto g_z_0_z_0_0_xz_xx_xy = buffer_1010_sddd[1801];

    auto g_z_0_z_0_0_xz_xx_xz = buffer_1010_sddd[1802];

    auto g_z_0_z_0_0_xz_xx_yy = buffer_1010_sddd[1803];

    auto g_z_0_z_0_0_xz_xx_yz = buffer_1010_sddd[1804];

    auto g_z_0_z_0_0_xz_xx_zz = buffer_1010_sddd[1805];

    auto g_z_0_z_0_0_xz_xy_xx = buffer_1010_sddd[1806];

    auto g_z_0_z_0_0_xz_xy_xy = buffer_1010_sddd[1807];

    auto g_z_0_z_0_0_xz_xy_xz = buffer_1010_sddd[1808];

    auto g_z_0_z_0_0_xz_xy_yy = buffer_1010_sddd[1809];

    auto g_z_0_z_0_0_xz_xy_yz = buffer_1010_sddd[1810];

    auto g_z_0_z_0_0_xz_xy_zz = buffer_1010_sddd[1811];

    auto g_z_0_z_0_0_xz_xz_xx = buffer_1010_sddd[1812];

    auto g_z_0_z_0_0_xz_xz_xy = buffer_1010_sddd[1813];

    auto g_z_0_z_0_0_xz_xz_xz = buffer_1010_sddd[1814];

    auto g_z_0_z_0_0_xz_xz_yy = buffer_1010_sddd[1815];

    auto g_z_0_z_0_0_xz_xz_yz = buffer_1010_sddd[1816];

    auto g_z_0_z_0_0_xz_xz_zz = buffer_1010_sddd[1817];

    auto g_z_0_z_0_0_xz_yy_xx = buffer_1010_sddd[1818];

    auto g_z_0_z_0_0_xz_yy_xy = buffer_1010_sddd[1819];

    auto g_z_0_z_0_0_xz_yy_xz = buffer_1010_sddd[1820];

    auto g_z_0_z_0_0_xz_yy_yy = buffer_1010_sddd[1821];

    auto g_z_0_z_0_0_xz_yy_yz = buffer_1010_sddd[1822];

    auto g_z_0_z_0_0_xz_yy_zz = buffer_1010_sddd[1823];

    auto g_z_0_z_0_0_xz_yz_xx = buffer_1010_sddd[1824];

    auto g_z_0_z_0_0_xz_yz_xy = buffer_1010_sddd[1825];

    auto g_z_0_z_0_0_xz_yz_xz = buffer_1010_sddd[1826];

    auto g_z_0_z_0_0_xz_yz_yy = buffer_1010_sddd[1827];

    auto g_z_0_z_0_0_xz_yz_yz = buffer_1010_sddd[1828];

    auto g_z_0_z_0_0_xz_yz_zz = buffer_1010_sddd[1829];

    auto g_z_0_z_0_0_xz_zz_xx = buffer_1010_sddd[1830];

    auto g_z_0_z_0_0_xz_zz_xy = buffer_1010_sddd[1831];

    auto g_z_0_z_0_0_xz_zz_xz = buffer_1010_sddd[1832];

    auto g_z_0_z_0_0_xz_zz_yy = buffer_1010_sddd[1833];

    auto g_z_0_z_0_0_xz_zz_yz = buffer_1010_sddd[1834];

    auto g_z_0_z_0_0_xz_zz_zz = buffer_1010_sddd[1835];

    auto g_z_0_z_0_0_yy_xx_xx = buffer_1010_sddd[1836];

    auto g_z_0_z_0_0_yy_xx_xy = buffer_1010_sddd[1837];

    auto g_z_0_z_0_0_yy_xx_xz = buffer_1010_sddd[1838];

    auto g_z_0_z_0_0_yy_xx_yy = buffer_1010_sddd[1839];

    auto g_z_0_z_0_0_yy_xx_yz = buffer_1010_sddd[1840];

    auto g_z_0_z_0_0_yy_xx_zz = buffer_1010_sddd[1841];

    auto g_z_0_z_0_0_yy_xy_xx = buffer_1010_sddd[1842];

    auto g_z_0_z_0_0_yy_xy_xy = buffer_1010_sddd[1843];

    auto g_z_0_z_0_0_yy_xy_xz = buffer_1010_sddd[1844];

    auto g_z_0_z_0_0_yy_xy_yy = buffer_1010_sddd[1845];

    auto g_z_0_z_0_0_yy_xy_yz = buffer_1010_sddd[1846];

    auto g_z_0_z_0_0_yy_xy_zz = buffer_1010_sddd[1847];

    auto g_z_0_z_0_0_yy_xz_xx = buffer_1010_sddd[1848];

    auto g_z_0_z_0_0_yy_xz_xy = buffer_1010_sddd[1849];

    auto g_z_0_z_0_0_yy_xz_xz = buffer_1010_sddd[1850];

    auto g_z_0_z_0_0_yy_xz_yy = buffer_1010_sddd[1851];

    auto g_z_0_z_0_0_yy_xz_yz = buffer_1010_sddd[1852];

    auto g_z_0_z_0_0_yy_xz_zz = buffer_1010_sddd[1853];

    auto g_z_0_z_0_0_yy_yy_xx = buffer_1010_sddd[1854];

    auto g_z_0_z_0_0_yy_yy_xy = buffer_1010_sddd[1855];

    auto g_z_0_z_0_0_yy_yy_xz = buffer_1010_sddd[1856];

    auto g_z_0_z_0_0_yy_yy_yy = buffer_1010_sddd[1857];

    auto g_z_0_z_0_0_yy_yy_yz = buffer_1010_sddd[1858];

    auto g_z_0_z_0_0_yy_yy_zz = buffer_1010_sddd[1859];

    auto g_z_0_z_0_0_yy_yz_xx = buffer_1010_sddd[1860];

    auto g_z_0_z_0_0_yy_yz_xy = buffer_1010_sddd[1861];

    auto g_z_0_z_0_0_yy_yz_xz = buffer_1010_sddd[1862];

    auto g_z_0_z_0_0_yy_yz_yy = buffer_1010_sddd[1863];

    auto g_z_0_z_0_0_yy_yz_yz = buffer_1010_sddd[1864];

    auto g_z_0_z_0_0_yy_yz_zz = buffer_1010_sddd[1865];

    auto g_z_0_z_0_0_yy_zz_xx = buffer_1010_sddd[1866];

    auto g_z_0_z_0_0_yy_zz_xy = buffer_1010_sddd[1867];

    auto g_z_0_z_0_0_yy_zz_xz = buffer_1010_sddd[1868];

    auto g_z_0_z_0_0_yy_zz_yy = buffer_1010_sddd[1869];

    auto g_z_0_z_0_0_yy_zz_yz = buffer_1010_sddd[1870];

    auto g_z_0_z_0_0_yy_zz_zz = buffer_1010_sddd[1871];

    auto g_z_0_z_0_0_yz_xx_xx = buffer_1010_sddd[1872];

    auto g_z_0_z_0_0_yz_xx_xy = buffer_1010_sddd[1873];

    auto g_z_0_z_0_0_yz_xx_xz = buffer_1010_sddd[1874];

    auto g_z_0_z_0_0_yz_xx_yy = buffer_1010_sddd[1875];

    auto g_z_0_z_0_0_yz_xx_yz = buffer_1010_sddd[1876];

    auto g_z_0_z_0_0_yz_xx_zz = buffer_1010_sddd[1877];

    auto g_z_0_z_0_0_yz_xy_xx = buffer_1010_sddd[1878];

    auto g_z_0_z_0_0_yz_xy_xy = buffer_1010_sddd[1879];

    auto g_z_0_z_0_0_yz_xy_xz = buffer_1010_sddd[1880];

    auto g_z_0_z_0_0_yz_xy_yy = buffer_1010_sddd[1881];

    auto g_z_0_z_0_0_yz_xy_yz = buffer_1010_sddd[1882];

    auto g_z_0_z_0_0_yz_xy_zz = buffer_1010_sddd[1883];

    auto g_z_0_z_0_0_yz_xz_xx = buffer_1010_sddd[1884];

    auto g_z_0_z_0_0_yz_xz_xy = buffer_1010_sddd[1885];

    auto g_z_0_z_0_0_yz_xz_xz = buffer_1010_sddd[1886];

    auto g_z_0_z_0_0_yz_xz_yy = buffer_1010_sddd[1887];

    auto g_z_0_z_0_0_yz_xz_yz = buffer_1010_sddd[1888];

    auto g_z_0_z_0_0_yz_xz_zz = buffer_1010_sddd[1889];

    auto g_z_0_z_0_0_yz_yy_xx = buffer_1010_sddd[1890];

    auto g_z_0_z_0_0_yz_yy_xy = buffer_1010_sddd[1891];

    auto g_z_0_z_0_0_yz_yy_xz = buffer_1010_sddd[1892];

    auto g_z_0_z_0_0_yz_yy_yy = buffer_1010_sddd[1893];

    auto g_z_0_z_0_0_yz_yy_yz = buffer_1010_sddd[1894];

    auto g_z_0_z_0_0_yz_yy_zz = buffer_1010_sddd[1895];

    auto g_z_0_z_0_0_yz_yz_xx = buffer_1010_sddd[1896];

    auto g_z_0_z_0_0_yz_yz_xy = buffer_1010_sddd[1897];

    auto g_z_0_z_0_0_yz_yz_xz = buffer_1010_sddd[1898];

    auto g_z_0_z_0_0_yz_yz_yy = buffer_1010_sddd[1899];

    auto g_z_0_z_0_0_yz_yz_yz = buffer_1010_sddd[1900];

    auto g_z_0_z_0_0_yz_yz_zz = buffer_1010_sddd[1901];

    auto g_z_0_z_0_0_yz_zz_xx = buffer_1010_sddd[1902];

    auto g_z_0_z_0_0_yz_zz_xy = buffer_1010_sddd[1903];

    auto g_z_0_z_0_0_yz_zz_xz = buffer_1010_sddd[1904];

    auto g_z_0_z_0_0_yz_zz_yy = buffer_1010_sddd[1905];

    auto g_z_0_z_0_0_yz_zz_yz = buffer_1010_sddd[1906];

    auto g_z_0_z_0_0_yz_zz_zz = buffer_1010_sddd[1907];

    auto g_z_0_z_0_0_zz_xx_xx = buffer_1010_sddd[1908];

    auto g_z_0_z_0_0_zz_xx_xy = buffer_1010_sddd[1909];

    auto g_z_0_z_0_0_zz_xx_xz = buffer_1010_sddd[1910];

    auto g_z_0_z_0_0_zz_xx_yy = buffer_1010_sddd[1911];

    auto g_z_0_z_0_0_zz_xx_yz = buffer_1010_sddd[1912];

    auto g_z_0_z_0_0_zz_xx_zz = buffer_1010_sddd[1913];

    auto g_z_0_z_0_0_zz_xy_xx = buffer_1010_sddd[1914];

    auto g_z_0_z_0_0_zz_xy_xy = buffer_1010_sddd[1915];

    auto g_z_0_z_0_0_zz_xy_xz = buffer_1010_sddd[1916];

    auto g_z_0_z_0_0_zz_xy_yy = buffer_1010_sddd[1917];

    auto g_z_0_z_0_0_zz_xy_yz = buffer_1010_sddd[1918];

    auto g_z_0_z_0_0_zz_xy_zz = buffer_1010_sddd[1919];

    auto g_z_0_z_0_0_zz_xz_xx = buffer_1010_sddd[1920];

    auto g_z_0_z_0_0_zz_xz_xy = buffer_1010_sddd[1921];

    auto g_z_0_z_0_0_zz_xz_xz = buffer_1010_sddd[1922];

    auto g_z_0_z_0_0_zz_xz_yy = buffer_1010_sddd[1923];

    auto g_z_0_z_0_0_zz_xz_yz = buffer_1010_sddd[1924];

    auto g_z_0_z_0_0_zz_xz_zz = buffer_1010_sddd[1925];

    auto g_z_0_z_0_0_zz_yy_xx = buffer_1010_sddd[1926];

    auto g_z_0_z_0_0_zz_yy_xy = buffer_1010_sddd[1927];

    auto g_z_0_z_0_0_zz_yy_xz = buffer_1010_sddd[1928];

    auto g_z_0_z_0_0_zz_yy_yy = buffer_1010_sddd[1929];

    auto g_z_0_z_0_0_zz_yy_yz = buffer_1010_sddd[1930];

    auto g_z_0_z_0_0_zz_yy_zz = buffer_1010_sddd[1931];

    auto g_z_0_z_0_0_zz_yz_xx = buffer_1010_sddd[1932];

    auto g_z_0_z_0_0_zz_yz_xy = buffer_1010_sddd[1933];

    auto g_z_0_z_0_0_zz_yz_xz = buffer_1010_sddd[1934];

    auto g_z_0_z_0_0_zz_yz_yy = buffer_1010_sddd[1935];

    auto g_z_0_z_0_0_zz_yz_yz = buffer_1010_sddd[1936];

    auto g_z_0_z_0_0_zz_yz_zz = buffer_1010_sddd[1937];

    auto g_z_0_z_0_0_zz_zz_xx = buffer_1010_sddd[1938];

    auto g_z_0_z_0_0_zz_zz_xy = buffer_1010_sddd[1939];

    auto g_z_0_z_0_0_zz_zz_xz = buffer_1010_sddd[1940];

    auto g_z_0_z_0_0_zz_zz_yy = buffer_1010_sddd[1941];

    auto g_z_0_z_0_0_zz_zz_yz = buffer_1010_sddd[1942];

    auto g_z_0_z_0_0_zz_zz_zz = buffer_1010_sddd[1943];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_0_x_0_0_xx_xx_xx, g_x_0_x_0_0_xx_xx_xy, g_x_0_x_0_0_xx_xx_xz, g_x_0_x_0_0_xx_xx_yy, g_x_0_x_0_0_xx_xx_yz, g_x_0_x_0_0_xx_xx_zz, g_x_xx_x_xx, g_x_xx_x_xy, g_x_xx_x_xz, g_x_xx_x_yy, g_x_xx_x_yz, g_x_xx_x_zz, g_x_xx_xxx_xx, g_x_xx_xxx_xy, g_x_xx_xxx_xz, g_x_xx_xxx_yy, g_x_xx_xxx_yz, g_x_xx_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xx_xx_xx[i] = -4.0 * g_x_xx_x_xx[i] * a_exp + 4.0 * g_x_xx_xxx_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_xx_xy[i] = -4.0 * g_x_xx_x_xy[i] * a_exp + 4.0 * g_x_xx_xxx_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_xx_xz[i] = -4.0 * g_x_xx_x_xz[i] * a_exp + 4.0 * g_x_xx_xxx_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_xx_yy[i] = -4.0 * g_x_xx_x_yy[i] * a_exp + 4.0 * g_x_xx_xxx_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_xx_yz[i] = -4.0 * g_x_xx_x_yz[i] * a_exp + 4.0 * g_x_xx_xxx_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_xx_zz[i] = -4.0 * g_x_xx_x_zz[i] * a_exp + 4.0 * g_x_xx_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_0_x_0_0_xx_xy_xx, g_x_0_x_0_0_xx_xy_xy, g_x_0_x_0_0_xx_xy_xz, g_x_0_x_0_0_xx_xy_yy, g_x_0_x_0_0_xx_xy_yz, g_x_0_x_0_0_xx_xy_zz, g_x_xx_xxy_xx, g_x_xx_xxy_xy, g_x_xx_xxy_xz, g_x_xx_xxy_yy, g_x_xx_xxy_yz, g_x_xx_xxy_zz, g_x_xx_y_xx, g_x_xx_y_xy, g_x_xx_y_xz, g_x_xx_y_yy, g_x_xx_y_yz, g_x_xx_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xx_xy_xx[i] = -2.0 * g_x_xx_y_xx[i] * a_exp + 4.0 * g_x_xx_xxy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_xy_xy[i] = -2.0 * g_x_xx_y_xy[i] * a_exp + 4.0 * g_x_xx_xxy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_xy_xz[i] = -2.0 * g_x_xx_y_xz[i] * a_exp + 4.0 * g_x_xx_xxy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_xy_yy[i] = -2.0 * g_x_xx_y_yy[i] * a_exp + 4.0 * g_x_xx_xxy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_xy_yz[i] = -2.0 * g_x_xx_y_yz[i] * a_exp + 4.0 * g_x_xx_xxy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_xy_zz[i] = -2.0 * g_x_xx_y_zz[i] * a_exp + 4.0 * g_x_xx_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_0_x_0_0_xx_xz_xx, g_x_0_x_0_0_xx_xz_xy, g_x_0_x_0_0_xx_xz_xz, g_x_0_x_0_0_xx_xz_yy, g_x_0_x_0_0_xx_xz_yz, g_x_0_x_0_0_xx_xz_zz, g_x_xx_xxz_xx, g_x_xx_xxz_xy, g_x_xx_xxz_xz, g_x_xx_xxz_yy, g_x_xx_xxz_yz, g_x_xx_xxz_zz, g_x_xx_z_xx, g_x_xx_z_xy, g_x_xx_z_xz, g_x_xx_z_yy, g_x_xx_z_yz, g_x_xx_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xx_xz_xx[i] = -2.0 * g_x_xx_z_xx[i] * a_exp + 4.0 * g_x_xx_xxz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_xz_xy[i] = -2.0 * g_x_xx_z_xy[i] * a_exp + 4.0 * g_x_xx_xxz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_xz_xz[i] = -2.0 * g_x_xx_z_xz[i] * a_exp + 4.0 * g_x_xx_xxz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_xz_yy[i] = -2.0 * g_x_xx_z_yy[i] * a_exp + 4.0 * g_x_xx_xxz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_xz_yz[i] = -2.0 * g_x_xx_z_yz[i] * a_exp + 4.0 * g_x_xx_xxz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_xz_zz[i] = -2.0 * g_x_xx_z_zz[i] * a_exp + 4.0 * g_x_xx_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_x_0_x_0_0_xx_yy_xx, g_x_0_x_0_0_xx_yy_xy, g_x_0_x_0_0_xx_yy_xz, g_x_0_x_0_0_xx_yy_yy, g_x_0_x_0_0_xx_yy_yz, g_x_0_x_0_0_xx_yy_zz, g_x_xx_xyy_xx, g_x_xx_xyy_xy, g_x_xx_xyy_xz, g_x_xx_xyy_yy, g_x_xx_xyy_yz, g_x_xx_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xx_yy_xx[i] = 4.0 * g_x_xx_xyy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_yy_xy[i] = 4.0 * g_x_xx_xyy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_yy_xz[i] = 4.0 * g_x_xx_xyy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_yy_yy[i] = 4.0 * g_x_xx_xyy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_yy_yz[i] = 4.0 * g_x_xx_xyy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_yy_zz[i] = 4.0 * g_x_xx_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_x_0_x_0_0_xx_yz_xx, g_x_0_x_0_0_xx_yz_xy, g_x_0_x_0_0_xx_yz_xz, g_x_0_x_0_0_xx_yz_yy, g_x_0_x_0_0_xx_yz_yz, g_x_0_x_0_0_xx_yz_zz, g_x_xx_xyz_xx, g_x_xx_xyz_xy, g_x_xx_xyz_xz, g_x_xx_xyz_yy, g_x_xx_xyz_yz, g_x_xx_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xx_yz_xx[i] = 4.0 * g_x_xx_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_yz_xy[i] = 4.0 * g_x_xx_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_yz_xz[i] = 4.0 * g_x_xx_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_yz_yy[i] = 4.0 * g_x_xx_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_yz_yz[i] = 4.0 * g_x_xx_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_yz_zz[i] = 4.0 * g_x_xx_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_x_0_x_0_0_xx_zz_xx, g_x_0_x_0_0_xx_zz_xy, g_x_0_x_0_0_xx_zz_xz, g_x_0_x_0_0_xx_zz_yy, g_x_0_x_0_0_xx_zz_yz, g_x_0_x_0_0_xx_zz_zz, g_x_xx_xzz_xx, g_x_xx_xzz_xy, g_x_xx_xzz_xz, g_x_xx_xzz_yy, g_x_xx_xzz_yz, g_x_xx_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xx_zz_xx[i] = 4.0 * g_x_xx_xzz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_zz_xy[i] = 4.0 * g_x_xx_xzz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_zz_xz[i] = 4.0 * g_x_xx_xzz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_zz_yy[i] = 4.0 * g_x_xx_xzz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_zz_yz[i] = 4.0 * g_x_xx_xzz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_zz_zz[i] = 4.0 * g_x_xx_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_x_0_x_0_0_xy_xx_xx, g_x_0_x_0_0_xy_xx_xy, g_x_0_x_0_0_xy_xx_xz, g_x_0_x_0_0_xy_xx_yy, g_x_0_x_0_0_xy_xx_yz, g_x_0_x_0_0_xy_xx_zz, g_x_xy_x_xx, g_x_xy_x_xy, g_x_xy_x_xz, g_x_xy_x_yy, g_x_xy_x_yz, g_x_xy_x_zz, g_x_xy_xxx_xx, g_x_xy_xxx_xy, g_x_xy_xxx_xz, g_x_xy_xxx_yy, g_x_xy_xxx_yz, g_x_xy_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xy_xx_xx[i] = -4.0 * g_x_xy_x_xx[i] * a_exp + 4.0 * g_x_xy_xxx_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_xx_xy[i] = -4.0 * g_x_xy_x_xy[i] * a_exp + 4.0 * g_x_xy_xxx_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_xx_xz[i] = -4.0 * g_x_xy_x_xz[i] * a_exp + 4.0 * g_x_xy_xxx_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_xx_yy[i] = -4.0 * g_x_xy_x_yy[i] * a_exp + 4.0 * g_x_xy_xxx_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_xx_yz[i] = -4.0 * g_x_xy_x_yz[i] * a_exp + 4.0 * g_x_xy_xxx_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_xx_zz[i] = -4.0 * g_x_xy_x_zz[i] * a_exp + 4.0 * g_x_xy_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_x_0_x_0_0_xy_xy_xx, g_x_0_x_0_0_xy_xy_xy, g_x_0_x_0_0_xy_xy_xz, g_x_0_x_0_0_xy_xy_yy, g_x_0_x_0_0_xy_xy_yz, g_x_0_x_0_0_xy_xy_zz, g_x_xy_xxy_xx, g_x_xy_xxy_xy, g_x_xy_xxy_xz, g_x_xy_xxy_yy, g_x_xy_xxy_yz, g_x_xy_xxy_zz, g_x_xy_y_xx, g_x_xy_y_xy, g_x_xy_y_xz, g_x_xy_y_yy, g_x_xy_y_yz, g_x_xy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xy_xy_xx[i] = -2.0 * g_x_xy_y_xx[i] * a_exp + 4.0 * g_x_xy_xxy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_xy_xy[i] = -2.0 * g_x_xy_y_xy[i] * a_exp + 4.0 * g_x_xy_xxy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_xy_xz[i] = -2.0 * g_x_xy_y_xz[i] * a_exp + 4.0 * g_x_xy_xxy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_xy_yy[i] = -2.0 * g_x_xy_y_yy[i] * a_exp + 4.0 * g_x_xy_xxy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_xy_yz[i] = -2.0 * g_x_xy_y_yz[i] * a_exp + 4.0 * g_x_xy_xxy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_xy_zz[i] = -2.0 * g_x_xy_y_zz[i] * a_exp + 4.0 * g_x_xy_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_x_0_x_0_0_xy_xz_xx, g_x_0_x_0_0_xy_xz_xy, g_x_0_x_0_0_xy_xz_xz, g_x_0_x_0_0_xy_xz_yy, g_x_0_x_0_0_xy_xz_yz, g_x_0_x_0_0_xy_xz_zz, g_x_xy_xxz_xx, g_x_xy_xxz_xy, g_x_xy_xxz_xz, g_x_xy_xxz_yy, g_x_xy_xxz_yz, g_x_xy_xxz_zz, g_x_xy_z_xx, g_x_xy_z_xy, g_x_xy_z_xz, g_x_xy_z_yy, g_x_xy_z_yz, g_x_xy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xy_xz_xx[i] = -2.0 * g_x_xy_z_xx[i] * a_exp + 4.0 * g_x_xy_xxz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_xz_xy[i] = -2.0 * g_x_xy_z_xy[i] * a_exp + 4.0 * g_x_xy_xxz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_xz_xz[i] = -2.0 * g_x_xy_z_xz[i] * a_exp + 4.0 * g_x_xy_xxz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_xz_yy[i] = -2.0 * g_x_xy_z_yy[i] * a_exp + 4.0 * g_x_xy_xxz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_xz_yz[i] = -2.0 * g_x_xy_z_yz[i] * a_exp + 4.0 * g_x_xy_xxz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_xz_zz[i] = -2.0 * g_x_xy_z_zz[i] * a_exp + 4.0 * g_x_xy_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_x_0_x_0_0_xy_yy_xx, g_x_0_x_0_0_xy_yy_xy, g_x_0_x_0_0_xy_yy_xz, g_x_0_x_0_0_xy_yy_yy, g_x_0_x_0_0_xy_yy_yz, g_x_0_x_0_0_xy_yy_zz, g_x_xy_xyy_xx, g_x_xy_xyy_xy, g_x_xy_xyy_xz, g_x_xy_xyy_yy, g_x_xy_xyy_yz, g_x_xy_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xy_yy_xx[i] = 4.0 * g_x_xy_xyy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_yy_xy[i] = 4.0 * g_x_xy_xyy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_yy_xz[i] = 4.0 * g_x_xy_xyy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_yy_yy[i] = 4.0 * g_x_xy_xyy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_yy_yz[i] = 4.0 * g_x_xy_xyy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_yy_zz[i] = 4.0 * g_x_xy_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_x_0_x_0_0_xy_yz_xx, g_x_0_x_0_0_xy_yz_xy, g_x_0_x_0_0_xy_yz_xz, g_x_0_x_0_0_xy_yz_yy, g_x_0_x_0_0_xy_yz_yz, g_x_0_x_0_0_xy_yz_zz, g_x_xy_xyz_xx, g_x_xy_xyz_xy, g_x_xy_xyz_xz, g_x_xy_xyz_yy, g_x_xy_xyz_yz, g_x_xy_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xy_yz_xx[i] = 4.0 * g_x_xy_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_yz_xy[i] = 4.0 * g_x_xy_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_yz_xz[i] = 4.0 * g_x_xy_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_yz_yy[i] = 4.0 * g_x_xy_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_yz_yz[i] = 4.0 * g_x_xy_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_yz_zz[i] = 4.0 * g_x_xy_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_x_0_x_0_0_xy_zz_xx, g_x_0_x_0_0_xy_zz_xy, g_x_0_x_0_0_xy_zz_xz, g_x_0_x_0_0_xy_zz_yy, g_x_0_x_0_0_xy_zz_yz, g_x_0_x_0_0_xy_zz_zz, g_x_xy_xzz_xx, g_x_xy_xzz_xy, g_x_xy_xzz_xz, g_x_xy_xzz_yy, g_x_xy_xzz_yz, g_x_xy_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xy_zz_xx[i] = 4.0 * g_x_xy_xzz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_zz_xy[i] = 4.0 * g_x_xy_xzz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_zz_xz[i] = 4.0 * g_x_xy_xzz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_zz_yy[i] = 4.0 * g_x_xy_xzz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_zz_yz[i] = 4.0 * g_x_xy_xzz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_zz_zz[i] = 4.0 * g_x_xy_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_x_0_x_0_0_xz_xx_xx, g_x_0_x_0_0_xz_xx_xy, g_x_0_x_0_0_xz_xx_xz, g_x_0_x_0_0_xz_xx_yy, g_x_0_x_0_0_xz_xx_yz, g_x_0_x_0_0_xz_xx_zz, g_x_xz_x_xx, g_x_xz_x_xy, g_x_xz_x_xz, g_x_xz_x_yy, g_x_xz_x_yz, g_x_xz_x_zz, g_x_xz_xxx_xx, g_x_xz_xxx_xy, g_x_xz_xxx_xz, g_x_xz_xxx_yy, g_x_xz_xxx_yz, g_x_xz_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xz_xx_xx[i] = -4.0 * g_x_xz_x_xx[i] * a_exp + 4.0 * g_x_xz_xxx_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_xx_xy[i] = -4.0 * g_x_xz_x_xy[i] * a_exp + 4.0 * g_x_xz_xxx_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_xx_xz[i] = -4.0 * g_x_xz_x_xz[i] * a_exp + 4.0 * g_x_xz_xxx_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_xx_yy[i] = -4.0 * g_x_xz_x_yy[i] * a_exp + 4.0 * g_x_xz_xxx_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_xx_yz[i] = -4.0 * g_x_xz_x_yz[i] * a_exp + 4.0 * g_x_xz_xxx_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_xx_zz[i] = -4.0 * g_x_xz_x_zz[i] * a_exp + 4.0 * g_x_xz_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_x_0_x_0_0_xz_xy_xx, g_x_0_x_0_0_xz_xy_xy, g_x_0_x_0_0_xz_xy_xz, g_x_0_x_0_0_xz_xy_yy, g_x_0_x_0_0_xz_xy_yz, g_x_0_x_0_0_xz_xy_zz, g_x_xz_xxy_xx, g_x_xz_xxy_xy, g_x_xz_xxy_xz, g_x_xz_xxy_yy, g_x_xz_xxy_yz, g_x_xz_xxy_zz, g_x_xz_y_xx, g_x_xz_y_xy, g_x_xz_y_xz, g_x_xz_y_yy, g_x_xz_y_yz, g_x_xz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xz_xy_xx[i] = -2.0 * g_x_xz_y_xx[i] * a_exp + 4.0 * g_x_xz_xxy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_xy_xy[i] = -2.0 * g_x_xz_y_xy[i] * a_exp + 4.0 * g_x_xz_xxy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_xy_xz[i] = -2.0 * g_x_xz_y_xz[i] * a_exp + 4.0 * g_x_xz_xxy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_xy_yy[i] = -2.0 * g_x_xz_y_yy[i] * a_exp + 4.0 * g_x_xz_xxy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_xy_yz[i] = -2.0 * g_x_xz_y_yz[i] * a_exp + 4.0 * g_x_xz_xxy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_xy_zz[i] = -2.0 * g_x_xz_y_zz[i] * a_exp + 4.0 * g_x_xz_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_x_0_x_0_0_xz_xz_xx, g_x_0_x_0_0_xz_xz_xy, g_x_0_x_0_0_xz_xz_xz, g_x_0_x_0_0_xz_xz_yy, g_x_0_x_0_0_xz_xz_yz, g_x_0_x_0_0_xz_xz_zz, g_x_xz_xxz_xx, g_x_xz_xxz_xy, g_x_xz_xxz_xz, g_x_xz_xxz_yy, g_x_xz_xxz_yz, g_x_xz_xxz_zz, g_x_xz_z_xx, g_x_xz_z_xy, g_x_xz_z_xz, g_x_xz_z_yy, g_x_xz_z_yz, g_x_xz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xz_xz_xx[i] = -2.0 * g_x_xz_z_xx[i] * a_exp + 4.0 * g_x_xz_xxz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_xz_xy[i] = -2.0 * g_x_xz_z_xy[i] * a_exp + 4.0 * g_x_xz_xxz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_xz_xz[i] = -2.0 * g_x_xz_z_xz[i] * a_exp + 4.0 * g_x_xz_xxz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_xz_yy[i] = -2.0 * g_x_xz_z_yy[i] * a_exp + 4.0 * g_x_xz_xxz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_xz_yz[i] = -2.0 * g_x_xz_z_yz[i] * a_exp + 4.0 * g_x_xz_xxz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_xz_zz[i] = -2.0 * g_x_xz_z_zz[i] * a_exp + 4.0 * g_x_xz_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_x_0_x_0_0_xz_yy_xx, g_x_0_x_0_0_xz_yy_xy, g_x_0_x_0_0_xz_yy_xz, g_x_0_x_0_0_xz_yy_yy, g_x_0_x_0_0_xz_yy_yz, g_x_0_x_0_0_xz_yy_zz, g_x_xz_xyy_xx, g_x_xz_xyy_xy, g_x_xz_xyy_xz, g_x_xz_xyy_yy, g_x_xz_xyy_yz, g_x_xz_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xz_yy_xx[i] = 4.0 * g_x_xz_xyy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_yy_xy[i] = 4.0 * g_x_xz_xyy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_yy_xz[i] = 4.0 * g_x_xz_xyy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_yy_yy[i] = 4.0 * g_x_xz_xyy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_yy_yz[i] = 4.0 * g_x_xz_xyy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_yy_zz[i] = 4.0 * g_x_xz_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_x_0_x_0_0_xz_yz_xx, g_x_0_x_0_0_xz_yz_xy, g_x_0_x_0_0_xz_yz_xz, g_x_0_x_0_0_xz_yz_yy, g_x_0_x_0_0_xz_yz_yz, g_x_0_x_0_0_xz_yz_zz, g_x_xz_xyz_xx, g_x_xz_xyz_xy, g_x_xz_xyz_xz, g_x_xz_xyz_yy, g_x_xz_xyz_yz, g_x_xz_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xz_yz_xx[i] = 4.0 * g_x_xz_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_yz_xy[i] = 4.0 * g_x_xz_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_yz_xz[i] = 4.0 * g_x_xz_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_yz_yy[i] = 4.0 * g_x_xz_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_yz_yz[i] = 4.0 * g_x_xz_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_yz_zz[i] = 4.0 * g_x_xz_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_x_0_x_0_0_xz_zz_xx, g_x_0_x_0_0_xz_zz_xy, g_x_0_x_0_0_xz_zz_xz, g_x_0_x_0_0_xz_zz_yy, g_x_0_x_0_0_xz_zz_yz, g_x_0_x_0_0_xz_zz_zz, g_x_xz_xzz_xx, g_x_xz_xzz_xy, g_x_xz_xzz_xz, g_x_xz_xzz_yy, g_x_xz_xzz_yz, g_x_xz_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xz_zz_xx[i] = 4.0 * g_x_xz_xzz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_zz_xy[i] = 4.0 * g_x_xz_xzz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_zz_xz[i] = 4.0 * g_x_xz_xzz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_zz_yy[i] = 4.0 * g_x_xz_xzz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_zz_yz[i] = 4.0 * g_x_xz_xzz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_zz_zz[i] = 4.0 * g_x_xz_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_x_0_x_0_0_yy_xx_xx, g_x_0_x_0_0_yy_xx_xy, g_x_0_x_0_0_yy_xx_xz, g_x_0_x_0_0_yy_xx_yy, g_x_0_x_0_0_yy_xx_yz, g_x_0_x_0_0_yy_xx_zz, g_x_yy_x_xx, g_x_yy_x_xy, g_x_yy_x_xz, g_x_yy_x_yy, g_x_yy_x_yz, g_x_yy_x_zz, g_x_yy_xxx_xx, g_x_yy_xxx_xy, g_x_yy_xxx_xz, g_x_yy_xxx_yy, g_x_yy_xxx_yz, g_x_yy_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_yy_xx_xx[i] = -4.0 * g_x_yy_x_xx[i] * a_exp + 4.0 * g_x_yy_xxx_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_xx_xy[i] = -4.0 * g_x_yy_x_xy[i] * a_exp + 4.0 * g_x_yy_xxx_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_xx_xz[i] = -4.0 * g_x_yy_x_xz[i] * a_exp + 4.0 * g_x_yy_xxx_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_xx_yy[i] = -4.0 * g_x_yy_x_yy[i] * a_exp + 4.0 * g_x_yy_xxx_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_xx_yz[i] = -4.0 * g_x_yy_x_yz[i] * a_exp + 4.0 * g_x_yy_xxx_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_xx_zz[i] = -4.0 * g_x_yy_x_zz[i] * a_exp + 4.0 * g_x_yy_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_x_0_x_0_0_yy_xy_xx, g_x_0_x_0_0_yy_xy_xy, g_x_0_x_0_0_yy_xy_xz, g_x_0_x_0_0_yy_xy_yy, g_x_0_x_0_0_yy_xy_yz, g_x_0_x_0_0_yy_xy_zz, g_x_yy_xxy_xx, g_x_yy_xxy_xy, g_x_yy_xxy_xz, g_x_yy_xxy_yy, g_x_yy_xxy_yz, g_x_yy_xxy_zz, g_x_yy_y_xx, g_x_yy_y_xy, g_x_yy_y_xz, g_x_yy_y_yy, g_x_yy_y_yz, g_x_yy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_yy_xy_xx[i] = -2.0 * g_x_yy_y_xx[i] * a_exp + 4.0 * g_x_yy_xxy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_xy_xy[i] = -2.0 * g_x_yy_y_xy[i] * a_exp + 4.0 * g_x_yy_xxy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_xy_xz[i] = -2.0 * g_x_yy_y_xz[i] * a_exp + 4.0 * g_x_yy_xxy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_xy_yy[i] = -2.0 * g_x_yy_y_yy[i] * a_exp + 4.0 * g_x_yy_xxy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_xy_yz[i] = -2.0 * g_x_yy_y_yz[i] * a_exp + 4.0 * g_x_yy_xxy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_xy_zz[i] = -2.0 * g_x_yy_y_zz[i] * a_exp + 4.0 * g_x_yy_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_x_0_x_0_0_yy_xz_xx, g_x_0_x_0_0_yy_xz_xy, g_x_0_x_0_0_yy_xz_xz, g_x_0_x_0_0_yy_xz_yy, g_x_0_x_0_0_yy_xz_yz, g_x_0_x_0_0_yy_xz_zz, g_x_yy_xxz_xx, g_x_yy_xxz_xy, g_x_yy_xxz_xz, g_x_yy_xxz_yy, g_x_yy_xxz_yz, g_x_yy_xxz_zz, g_x_yy_z_xx, g_x_yy_z_xy, g_x_yy_z_xz, g_x_yy_z_yy, g_x_yy_z_yz, g_x_yy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_yy_xz_xx[i] = -2.0 * g_x_yy_z_xx[i] * a_exp + 4.0 * g_x_yy_xxz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_xz_xy[i] = -2.0 * g_x_yy_z_xy[i] * a_exp + 4.0 * g_x_yy_xxz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_xz_xz[i] = -2.0 * g_x_yy_z_xz[i] * a_exp + 4.0 * g_x_yy_xxz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_xz_yy[i] = -2.0 * g_x_yy_z_yy[i] * a_exp + 4.0 * g_x_yy_xxz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_xz_yz[i] = -2.0 * g_x_yy_z_yz[i] * a_exp + 4.0 * g_x_yy_xxz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_xz_zz[i] = -2.0 * g_x_yy_z_zz[i] * a_exp + 4.0 * g_x_yy_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_x_0_x_0_0_yy_yy_xx, g_x_0_x_0_0_yy_yy_xy, g_x_0_x_0_0_yy_yy_xz, g_x_0_x_0_0_yy_yy_yy, g_x_0_x_0_0_yy_yy_yz, g_x_0_x_0_0_yy_yy_zz, g_x_yy_xyy_xx, g_x_yy_xyy_xy, g_x_yy_xyy_xz, g_x_yy_xyy_yy, g_x_yy_xyy_yz, g_x_yy_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_yy_yy_xx[i] = 4.0 * g_x_yy_xyy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_yy_xy[i] = 4.0 * g_x_yy_xyy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_yy_xz[i] = 4.0 * g_x_yy_xyy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_yy_yy[i] = 4.0 * g_x_yy_xyy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_yy_yz[i] = 4.0 * g_x_yy_xyy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_yy_zz[i] = 4.0 * g_x_yy_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_x_0_x_0_0_yy_yz_xx, g_x_0_x_0_0_yy_yz_xy, g_x_0_x_0_0_yy_yz_xz, g_x_0_x_0_0_yy_yz_yy, g_x_0_x_0_0_yy_yz_yz, g_x_0_x_0_0_yy_yz_zz, g_x_yy_xyz_xx, g_x_yy_xyz_xy, g_x_yy_xyz_xz, g_x_yy_xyz_yy, g_x_yy_xyz_yz, g_x_yy_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_yy_yz_xx[i] = 4.0 * g_x_yy_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_yz_xy[i] = 4.0 * g_x_yy_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_yz_xz[i] = 4.0 * g_x_yy_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_yz_yy[i] = 4.0 * g_x_yy_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_yz_yz[i] = 4.0 * g_x_yy_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_yz_zz[i] = 4.0 * g_x_yy_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_x_0_x_0_0_yy_zz_xx, g_x_0_x_0_0_yy_zz_xy, g_x_0_x_0_0_yy_zz_xz, g_x_0_x_0_0_yy_zz_yy, g_x_0_x_0_0_yy_zz_yz, g_x_0_x_0_0_yy_zz_zz, g_x_yy_xzz_xx, g_x_yy_xzz_xy, g_x_yy_xzz_xz, g_x_yy_xzz_yy, g_x_yy_xzz_yz, g_x_yy_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_yy_zz_xx[i] = 4.0 * g_x_yy_xzz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_zz_xy[i] = 4.0 * g_x_yy_xzz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_zz_xz[i] = 4.0 * g_x_yy_xzz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_zz_yy[i] = 4.0 * g_x_yy_xzz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_zz_yz[i] = 4.0 * g_x_yy_xzz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_zz_zz[i] = 4.0 * g_x_yy_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_x_0_x_0_0_yz_xx_xx, g_x_0_x_0_0_yz_xx_xy, g_x_0_x_0_0_yz_xx_xz, g_x_0_x_0_0_yz_xx_yy, g_x_0_x_0_0_yz_xx_yz, g_x_0_x_0_0_yz_xx_zz, g_x_yz_x_xx, g_x_yz_x_xy, g_x_yz_x_xz, g_x_yz_x_yy, g_x_yz_x_yz, g_x_yz_x_zz, g_x_yz_xxx_xx, g_x_yz_xxx_xy, g_x_yz_xxx_xz, g_x_yz_xxx_yy, g_x_yz_xxx_yz, g_x_yz_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_yz_xx_xx[i] = -4.0 * g_x_yz_x_xx[i] * a_exp + 4.0 * g_x_yz_xxx_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_xx_xy[i] = -4.0 * g_x_yz_x_xy[i] * a_exp + 4.0 * g_x_yz_xxx_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_xx_xz[i] = -4.0 * g_x_yz_x_xz[i] * a_exp + 4.0 * g_x_yz_xxx_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_xx_yy[i] = -4.0 * g_x_yz_x_yy[i] * a_exp + 4.0 * g_x_yz_xxx_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_xx_yz[i] = -4.0 * g_x_yz_x_yz[i] * a_exp + 4.0 * g_x_yz_xxx_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_xx_zz[i] = -4.0 * g_x_yz_x_zz[i] * a_exp + 4.0 * g_x_yz_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_x_0_x_0_0_yz_xy_xx, g_x_0_x_0_0_yz_xy_xy, g_x_0_x_0_0_yz_xy_xz, g_x_0_x_0_0_yz_xy_yy, g_x_0_x_0_0_yz_xy_yz, g_x_0_x_0_0_yz_xy_zz, g_x_yz_xxy_xx, g_x_yz_xxy_xy, g_x_yz_xxy_xz, g_x_yz_xxy_yy, g_x_yz_xxy_yz, g_x_yz_xxy_zz, g_x_yz_y_xx, g_x_yz_y_xy, g_x_yz_y_xz, g_x_yz_y_yy, g_x_yz_y_yz, g_x_yz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_yz_xy_xx[i] = -2.0 * g_x_yz_y_xx[i] * a_exp + 4.0 * g_x_yz_xxy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_xy_xy[i] = -2.0 * g_x_yz_y_xy[i] * a_exp + 4.0 * g_x_yz_xxy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_xy_xz[i] = -2.0 * g_x_yz_y_xz[i] * a_exp + 4.0 * g_x_yz_xxy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_xy_yy[i] = -2.0 * g_x_yz_y_yy[i] * a_exp + 4.0 * g_x_yz_xxy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_xy_yz[i] = -2.0 * g_x_yz_y_yz[i] * a_exp + 4.0 * g_x_yz_xxy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_xy_zz[i] = -2.0 * g_x_yz_y_zz[i] * a_exp + 4.0 * g_x_yz_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_x_0_x_0_0_yz_xz_xx, g_x_0_x_0_0_yz_xz_xy, g_x_0_x_0_0_yz_xz_xz, g_x_0_x_0_0_yz_xz_yy, g_x_0_x_0_0_yz_xz_yz, g_x_0_x_0_0_yz_xz_zz, g_x_yz_xxz_xx, g_x_yz_xxz_xy, g_x_yz_xxz_xz, g_x_yz_xxz_yy, g_x_yz_xxz_yz, g_x_yz_xxz_zz, g_x_yz_z_xx, g_x_yz_z_xy, g_x_yz_z_xz, g_x_yz_z_yy, g_x_yz_z_yz, g_x_yz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_yz_xz_xx[i] = -2.0 * g_x_yz_z_xx[i] * a_exp + 4.0 * g_x_yz_xxz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_xz_xy[i] = -2.0 * g_x_yz_z_xy[i] * a_exp + 4.0 * g_x_yz_xxz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_xz_xz[i] = -2.0 * g_x_yz_z_xz[i] * a_exp + 4.0 * g_x_yz_xxz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_xz_yy[i] = -2.0 * g_x_yz_z_yy[i] * a_exp + 4.0 * g_x_yz_xxz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_xz_yz[i] = -2.0 * g_x_yz_z_yz[i] * a_exp + 4.0 * g_x_yz_xxz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_xz_zz[i] = -2.0 * g_x_yz_z_zz[i] * a_exp + 4.0 * g_x_yz_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_x_0_x_0_0_yz_yy_xx, g_x_0_x_0_0_yz_yy_xy, g_x_0_x_0_0_yz_yy_xz, g_x_0_x_0_0_yz_yy_yy, g_x_0_x_0_0_yz_yy_yz, g_x_0_x_0_0_yz_yy_zz, g_x_yz_xyy_xx, g_x_yz_xyy_xy, g_x_yz_xyy_xz, g_x_yz_xyy_yy, g_x_yz_xyy_yz, g_x_yz_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_yz_yy_xx[i] = 4.0 * g_x_yz_xyy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_yy_xy[i] = 4.0 * g_x_yz_xyy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_yy_xz[i] = 4.0 * g_x_yz_xyy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_yy_yy[i] = 4.0 * g_x_yz_xyy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_yy_yz[i] = 4.0 * g_x_yz_xyy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_yy_zz[i] = 4.0 * g_x_yz_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_x_0_x_0_0_yz_yz_xx, g_x_0_x_0_0_yz_yz_xy, g_x_0_x_0_0_yz_yz_xz, g_x_0_x_0_0_yz_yz_yy, g_x_0_x_0_0_yz_yz_yz, g_x_0_x_0_0_yz_yz_zz, g_x_yz_xyz_xx, g_x_yz_xyz_xy, g_x_yz_xyz_xz, g_x_yz_xyz_yy, g_x_yz_xyz_yz, g_x_yz_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_yz_yz_xx[i] = 4.0 * g_x_yz_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_yz_xy[i] = 4.0 * g_x_yz_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_yz_xz[i] = 4.0 * g_x_yz_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_yz_yy[i] = 4.0 * g_x_yz_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_yz_yz[i] = 4.0 * g_x_yz_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_yz_zz[i] = 4.0 * g_x_yz_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_x_0_x_0_0_yz_zz_xx, g_x_0_x_0_0_yz_zz_xy, g_x_0_x_0_0_yz_zz_xz, g_x_0_x_0_0_yz_zz_yy, g_x_0_x_0_0_yz_zz_yz, g_x_0_x_0_0_yz_zz_zz, g_x_yz_xzz_xx, g_x_yz_xzz_xy, g_x_yz_xzz_xz, g_x_yz_xzz_yy, g_x_yz_xzz_yz, g_x_yz_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_yz_zz_xx[i] = 4.0 * g_x_yz_xzz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_zz_xy[i] = 4.0 * g_x_yz_xzz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_zz_xz[i] = 4.0 * g_x_yz_xzz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_zz_yy[i] = 4.0 * g_x_yz_xzz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_zz_yz[i] = 4.0 * g_x_yz_xzz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_zz_zz[i] = 4.0 * g_x_yz_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_x_0_x_0_0_zz_xx_xx, g_x_0_x_0_0_zz_xx_xy, g_x_0_x_0_0_zz_xx_xz, g_x_0_x_0_0_zz_xx_yy, g_x_0_x_0_0_zz_xx_yz, g_x_0_x_0_0_zz_xx_zz, g_x_zz_x_xx, g_x_zz_x_xy, g_x_zz_x_xz, g_x_zz_x_yy, g_x_zz_x_yz, g_x_zz_x_zz, g_x_zz_xxx_xx, g_x_zz_xxx_xy, g_x_zz_xxx_xz, g_x_zz_xxx_yy, g_x_zz_xxx_yz, g_x_zz_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_zz_xx_xx[i] = -4.0 * g_x_zz_x_xx[i] * a_exp + 4.0 * g_x_zz_xxx_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_xx_xy[i] = -4.0 * g_x_zz_x_xy[i] * a_exp + 4.0 * g_x_zz_xxx_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_xx_xz[i] = -4.0 * g_x_zz_x_xz[i] * a_exp + 4.0 * g_x_zz_xxx_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_xx_yy[i] = -4.0 * g_x_zz_x_yy[i] * a_exp + 4.0 * g_x_zz_xxx_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_xx_yz[i] = -4.0 * g_x_zz_x_yz[i] * a_exp + 4.0 * g_x_zz_xxx_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_xx_zz[i] = -4.0 * g_x_zz_x_zz[i] * a_exp + 4.0 * g_x_zz_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_x_0_x_0_0_zz_xy_xx, g_x_0_x_0_0_zz_xy_xy, g_x_0_x_0_0_zz_xy_xz, g_x_0_x_0_0_zz_xy_yy, g_x_0_x_0_0_zz_xy_yz, g_x_0_x_0_0_zz_xy_zz, g_x_zz_xxy_xx, g_x_zz_xxy_xy, g_x_zz_xxy_xz, g_x_zz_xxy_yy, g_x_zz_xxy_yz, g_x_zz_xxy_zz, g_x_zz_y_xx, g_x_zz_y_xy, g_x_zz_y_xz, g_x_zz_y_yy, g_x_zz_y_yz, g_x_zz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_zz_xy_xx[i] = -2.0 * g_x_zz_y_xx[i] * a_exp + 4.0 * g_x_zz_xxy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_xy_xy[i] = -2.0 * g_x_zz_y_xy[i] * a_exp + 4.0 * g_x_zz_xxy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_xy_xz[i] = -2.0 * g_x_zz_y_xz[i] * a_exp + 4.0 * g_x_zz_xxy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_xy_yy[i] = -2.0 * g_x_zz_y_yy[i] * a_exp + 4.0 * g_x_zz_xxy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_xy_yz[i] = -2.0 * g_x_zz_y_yz[i] * a_exp + 4.0 * g_x_zz_xxy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_xy_zz[i] = -2.0 * g_x_zz_y_zz[i] * a_exp + 4.0 * g_x_zz_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_x_0_x_0_0_zz_xz_xx, g_x_0_x_0_0_zz_xz_xy, g_x_0_x_0_0_zz_xz_xz, g_x_0_x_0_0_zz_xz_yy, g_x_0_x_0_0_zz_xz_yz, g_x_0_x_0_0_zz_xz_zz, g_x_zz_xxz_xx, g_x_zz_xxz_xy, g_x_zz_xxz_xz, g_x_zz_xxz_yy, g_x_zz_xxz_yz, g_x_zz_xxz_zz, g_x_zz_z_xx, g_x_zz_z_xy, g_x_zz_z_xz, g_x_zz_z_yy, g_x_zz_z_yz, g_x_zz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_zz_xz_xx[i] = -2.0 * g_x_zz_z_xx[i] * a_exp + 4.0 * g_x_zz_xxz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_xz_xy[i] = -2.0 * g_x_zz_z_xy[i] * a_exp + 4.0 * g_x_zz_xxz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_xz_xz[i] = -2.0 * g_x_zz_z_xz[i] * a_exp + 4.0 * g_x_zz_xxz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_xz_yy[i] = -2.0 * g_x_zz_z_yy[i] * a_exp + 4.0 * g_x_zz_xxz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_xz_yz[i] = -2.0 * g_x_zz_z_yz[i] * a_exp + 4.0 * g_x_zz_xxz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_xz_zz[i] = -2.0 * g_x_zz_z_zz[i] * a_exp + 4.0 * g_x_zz_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_x_0_x_0_0_zz_yy_xx, g_x_0_x_0_0_zz_yy_xy, g_x_0_x_0_0_zz_yy_xz, g_x_0_x_0_0_zz_yy_yy, g_x_0_x_0_0_zz_yy_yz, g_x_0_x_0_0_zz_yy_zz, g_x_zz_xyy_xx, g_x_zz_xyy_xy, g_x_zz_xyy_xz, g_x_zz_xyy_yy, g_x_zz_xyy_yz, g_x_zz_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_zz_yy_xx[i] = 4.0 * g_x_zz_xyy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_yy_xy[i] = 4.0 * g_x_zz_xyy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_yy_xz[i] = 4.0 * g_x_zz_xyy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_yy_yy[i] = 4.0 * g_x_zz_xyy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_yy_yz[i] = 4.0 * g_x_zz_xyy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_yy_zz[i] = 4.0 * g_x_zz_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_x_0_x_0_0_zz_yz_xx, g_x_0_x_0_0_zz_yz_xy, g_x_0_x_0_0_zz_yz_xz, g_x_0_x_0_0_zz_yz_yy, g_x_0_x_0_0_zz_yz_yz, g_x_0_x_0_0_zz_yz_zz, g_x_zz_xyz_xx, g_x_zz_xyz_xy, g_x_zz_xyz_xz, g_x_zz_xyz_yy, g_x_zz_xyz_yz, g_x_zz_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_zz_yz_xx[i] = 4.0 * g_x_zz_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_yz_xy[i] = 4.0 * g_x_zz_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_yz_xz[i] = 4.0 * g_x_zz_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_yz_yy[i] = 4.0 * g_x_zz_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_yz_yz[i] = 4.0 * g_x_zz_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_yz_zz[i] = 4.0 * g_x_zz_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_x_0_x_0_0_zz_zz_xx, g_x_0_x_0_0_zz_zz_xy, g_x_0_x_0_0_zz_zz_xz, g_x_0_x_0_0_zz_zz_yy, g_x_0_x_0_0_zz_zz_yz, g_x_0_x_0_0_zz_zz_zz, g_x_zz_xzz_xx, g_x_zz_xzz_xy, g_x_zz_xzz_xz, g_x_zz_xzz_yy, g_x_zz_xzz_yz, g_x_zz_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_zz_zz_xx[i] = 4.0 * g_x_zz_xzz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_zz_xy[i] = 4.0 * g_x_zz_xzz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_zz_xz[i] = 4.0 * g_x_zz_xzz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_zz_yy[i] = 4.0 * g_x_zz_xzz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_zz_yz[i] = 4.0 * g_x_zz_xzz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_zz_zz[i] = 4.0 * g_x_zz_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_x_0_y_0_0_xx_xx_xx, g_x_0_y_0_0_xx_xx_xy, g_x_0_y_0_0_xx_xx_xz, g_x_0_y_0_0_xx_xx_yy, g_x_0_y_0_0_xx_xx_yz, g_x_0_y_0_0_xx_xx_zz, g_x_xx_xxy_xx, g_x_xx_xxy_xy, g_x_xx_xxy_xz, g_x_xx_xxy_yy, g_x_xx_xxy_yz, g_x_xx_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xx_xx_xx[i] = 4.0 * g_x_xx_xxy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_xx_xy[i] = 4.0 * g_x_xx_xxy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_xx_xz[i] = 4.0 * g_x_xx_xxy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_xx_yy[i] = 4.0 * g_x_xx_xxy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_xx_yz[i] = 4.0 * g_x_xx_xxy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_xx_zz[i] = 4.0 * g_x_xx_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_x_0_y_0_0_xx_xy_xx, g_x_0_y_0_0_xx_xy_xy, g_x_0_y_0_0_xx_xy_xz, g_x_0_y_0_0_xx_xy_yy, g_x_0_y_0_0_xx_xy_yz, g_x_0_y_0_0_xx_xy_zz, g_x_xx_x_xx, g_x_xx_x_xy, g_x_xx_x_xz, g_x_xx_x_yy, g_x_xx_x_yz, g_x_xx_x_zz, g_x_xx_xyy_xx, g_x_xx_xyy_xy, g_x_xx_xyy_xz, g_x_xx_xyy_yy, g_x_xx_xyy_yz, g_x_xx_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xx_xy_xx[i] = -2.0 * g_x_xx_x_xx[i] * a_exp + 4.0 * g_x_xx_xyy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_xy_xy[i] = -2.0 * g_x_xx_x_xy[i] * a_exp + 4.0 * g_x_xx_xyy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_xy_xz[i] = -2.0 * g_x_xx_x_xz[i] * a_exp + 4.0 * g_x_xx_xyy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_xy_yy[i] = -2.0 * g_x_xx_x_yy[i] * a_exp + 4.0 * g_x_xx_xyy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_xy_yz[i] = -2.0 * g_x_xx_x_yz[i] * a_exp + 4.0 * g_x_xx_xyy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_xy_zz[i] = -2.0 * g_x_xx_x_zz[i] * a_exp + 4.0 * g_x_xx_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_x_0_y_0_0_xx_xz_xx, g_x_0_y_0_0_xx_xz_xy, g_x_0_y_0_0_xx_xz_xz, g_x_0_y_0_0_xx_xz_yy, g_x_0_y_0_0_xx_xz_yz, g_x_0_y_0_0_xx_xz_zz, g_x_xx_xyz_xx, g_x_xx_xyz_xy, g_x_xx_xyz_xz, g_x_xx_xyz_yy, g_x_xx_xyz_yz, g_x_xx_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xx_xz_xx[i] = 4.0 * g_x_xx_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_xz_xy[i] = 4.0 * g_x_xx_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_xz_xz[i] = 4.0 * g_x_xx_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_xz_yy[i] = 4.0 * g_x_xx_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_xz_yz[i] = 4.0 * g_x_xx_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_xz_zz[i] = 4.0 * g_x_xx_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_x_0_y_0_0_xx_yy_xx, g_x_0_y_0_0_xx_yy_xy, g_x_0_y_0_0_xx_yy_xz, g_x_0_y_0_0_xx_yy_yy, g_x_0_y_0_0_xx_yy_yz, g_x_0_y_0_0_xx_yy_zz, g_x_xx_y_xx, g_x_xx_y_xy, g_x_xx_y_xz, g_x_xx_y_yy, g_x_xx_y_yz, g_x_xx_y_zz, g_x_xx_yyy_xx, g_x_xx_yyy_xy, g_x_xx_yyy_xz, g_x_xx_yyy_yy, g_x_xx_yyy_yz, g_x_xx_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xx_yy_xx[i] = -4.0 * g_x_xx_y_xx[i] * a_exp + 4.0 * g_x_xx_yyy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_yy_xy[i] = -4.0 * g_x_xx_y_xy[i] * a_exp + 4.0 * g_x_xx_yyy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_yy_xz[i] = -4.0 * g_x_xx_y_xz[i] * a_exp + 4.0 * g_x_xx_yyy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_yy_yy[i] = -4.0 * g_x_xx_y_yy[i] * a_exp + 4.0 * g_x_xx_yyy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_yy_yz[i] = -4.0 * g_x_xx_y_yz[i] * a_exp + 4.0 * g_x_xx_yyy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_yy_zz[i] = -4.0 * g_x_xx_y_zz[i] * a_exp + 4.0 * g_x_xx_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_x_0_y_0_0_xx_yz_xx, g_x_0_y_0_0_xx_yz_xy, g_x_0_y_0_0_xx_yz_xz, g_x_0_y_0_0_xx_yz_yy, g_x_0_y_0_0_xx_yz_yz, g_x_0_y_0_0_xx_yz_zz, g_x_xx_yyz_xx, g_x_xx_yyz_xy, g_x_xx_yyz_xz, g_x_xx_yyz_yy, g_x_xx_yyz_yz, g_x_xx_yyz_zz, g_x_xx_z_xx, g_x_xx_z_xy, g_x_xx_z_xz, g_x_xx_z_yy, g_x_xx_z_yz, g_x_xx_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xx_yz_xx[i] = -2.0 * g_x_xx_z_xx[i] * a_exp + 4.0 * g_x_xx_yyz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_yz_xy[i] = -2.0 * g_x_xx_z_xy[i] * a_exp + 4.0 * g_x_xx_yyz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_yz_xz[i] = -2.0 * g_x_xx_z_xz[i] * a_exp + 4.0 * g_x_xx_yyz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_yz_yy[i] = -2.0 * g_x_xx_z_yy[i] * a_exp + 4.0 * g_x_xx_yyz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_yz_yz[i] = -2.0 * g_x_xx_z_yz[i] * a_exp + 4.0 * g_x_xx_yyz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_yz_zz[i] = -2.0 * g_x_xx_z_zz[i] * a_exp + 4.0 * g_x_xx_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_x_0_y_0_0_xx_zz_xx, g_x_0_y_0_0_xx_zz_xy, g_x_0_y_0_0_xx_zz_xz, g_x_0_y_0_0_xx_zz_yy, g_x_0_y_0_0_xx_zz_yz, g_x_0_y_0_0_xx_zz_zz, g_x_xx_yzz_xx, g_x_xx_yzz_xy, g_x_xx_yzz_xz, g_x_xx_yzz_yy, g_x_xx_yzz_yz, g_x_xx_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xx_zz_xx[i] = 4.0 * g_x_xx_yzz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_zz_xy[i] = 4.0 * g_x_xx_yzz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_zz_xz[i] = 4.0 * g_x_xx_yzz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_zz_yy[i] = 4.0 * g_x_xx_yzz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_zz_yz[i] = 4.0 * g_x_xx_yzz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_zz_zz[i] = 4.0 * g_x_xx_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_x_0_y_0_0_xy_xx_xx, g_x_0_y_0_0_xy_xx_xy, g_x_0_y_0_0_xy_xx_xz, g_x_0_y_0_0_xy_xx_yy, g_x_0_y_0_0_xy_xx_yz, g_x_0_y_0_0_xy_xx_zz, g_x_xy_xxy_xx, g_x_xy_xxy_xy, g_x_xy_xxy_xz, g_x_xy_xxy_yy, g_x_xy_xxy_yz, g_x_xy_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xy_xx_xx[i] = 4.0 * g_x_xy_xxy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_xx_xy[i] = 4.0 * g_x_xy_xxy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_xx_xz[i] = 4.0 * g_x_xy_xxy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_xx_yy[i] = 4.0 * g_x_xy_xxy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_xx_yz[i] = 4.0 * g_x_xy_xxy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_xx_zz[i] = 4.0 * g_x_xy_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_x_0_y_0_0_xy_xy_xx, g_x_0_y_0_0_xy_xy_xy, g_x_0_y_0_0_xy_xy_xz, g_x_0_y_0_0_xy_xy_yy, g_x_0_y_0_0_xy_xy_yz, g_x_0_y_0_0_xy_xy_zz, g_x_xy_x_xx, g_x_xy_x_xy, g_x_xy_x_xz, g_x_xy_x_yy, g_x_xy_x_yz, g_x_xy_x_zz, g_x_xy_xyy_xx, g_x_xy_xyy_xy, g_x_xy_xyy_xz, g_x_xy_xyy_yy, g_x_xy_xyy_yz, g_x_xy_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xy_xy_xx[i] = -2.0 * g_x_xy_x_xx[i] * a_exp + 4.0 * g_x_xy_xyy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_xy_xy[i] = -2.0 * g_x_xy_x_xy[i] * a_exp + 4.0 * g_x_xy_xyy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_xy_xz[i] = -2.0 * g_x_xy_x_xz[i] * a_exp + 4.0 * g_x_xy_xyy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_xy_yy[i] = -2.0 * g_x_xy_x_yy[i] * a_exp + 4.0 * g_x_xy_xyy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_xy_yz[i] = -2.0 * g_x_xy_x_yz[i] * a_exp + 4.0 * g_x_xy_xyy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_xy_zz[i] = -2.0 * g_x_xy_x_zz[i] * a_exp + 4.0 * g_x_xy_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_x_0_y_0_0_xy_xz_xx, g_x_0_y_0_0_xy_xz_xy, g_x_0_y_0_0_xy_xz_xz, g_x_0_y_0_0_xy_xz_yy, g_x_0_y_0_0_xy_xz_yz, g_x_0_y_0_0_xy_xz_zz, g_x_xy_xyz_xx, g_x_xy_xyz_xy, g_x_xy_xyz_xz, g_x_xy_xyz_yy, g_x_xy_xyz_yz, g_x_xy_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xy_xz_xx[i] = 4.0 * g_x_xy_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_xz_xy[i] = 4.0 * g_x_xy_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_xz_xz[i] = 4.0 * g_x_xy_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_xz_yy[i] = 4.0 * g_x_xy_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_xz_yz[i] = 4.0 * g_x_xy_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_xz_zz[i] = 4.0 * g_x_xy_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_x_0_y_0_0_xy_yy_xx, g_x_0_y_0_0_xy_yy_xy, g_x_0_y_0_0_xy_yy_xz, g_x_0_y_0_0_xy_yy_yy, g_x_0_y_0_0_xy_yy_yz, g_x_0_y_0_0_xy_yy_zz, g_x_xy_y_xx, g_x_xy_y_xy, g_x_xy_y_xz, g_x_xy_y_yy, g_x_xy_y_yz, g_x_xy_y_zz, g_x_xy_yyy_xx, g_x_xy_yyy_xy, g_x_xy_yyy_xz, g_x_xy_yyy_yy, g_x_xy_yyy_yz, g_x_xy_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xy_yy_xx[i] = -4.0 * g_x_xy_y_xx[i] * a_exp + 4.0 * g_x_xy_yyy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_yy_xy[i] = -4.0 * g_x_xy_y_xy[i] * a_exp + 4.0 * g_x_xy_yyy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_yy_xz[i] = -4.0 * g_x_xy_y_xz[i] * a_exp + 4.0 * g_x_xy_yyy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_yy_yy[i] = -4.0 * g_x_xy_y_yy[i] * a_exp + 4.0 * g_x_xy_yyy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_yy_yz[i] = -4.0 * g_x_xy_y_yz[i] * a_exp + 4.0 * g_x_xy_yyy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_yy_zz[i] = -4.0 * g_x_xy_y_zz[i] * a_exp + 4.0 * g_x_xy_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_x_0_y_0_0_xy_yz_xx, g_x_0_y_0_0_xy_yz_xy, g_x_0_y_0_0_xy_yz_xz, g_x_0_y_0_0_xy_yz_yy, g_x_0_y_0_0_xy_yz_yz, g_x_0_y_0_0_xy_yz_zz, g_x_xy_yyz_xx, g_x_xy_yyz_xy, g_x_xy_yyz_xz, g_x_xy_yyz_yy, g_x_xy_yyz_yz, g_x_xy_yyz_zz, g_x_xy_z_xx, g_x_xy_z_xy, g_x_xy_z_xz, g_x_xy_z_yy, g_x_xy_z_yz, g_x_xy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xy_yz_xx[i] = -2.0 * g_x_xy_z_xx[i] * a_exp + 4.0 * g_x_xy_yyz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_yz_xy[i] = -2.0 * g_x_xy_z_xy[i] * a_exp + 4.0 * g_x_xy_yyz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_yz_xz[i] = -2.0 * g_x_xy_z_xz[i] * a_exp + 4.0 * g_x_xy_yyz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_yz_yy[i] = -2.0 * g_x_xy_z_yy[i] * a_exp + 4.0 * g_x_xy_yyz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_yz_yz[i] = -2.0 * g_x_xy_z_yz[i] * a_exp + 4.0 * g_x_xy_yyz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_yz_zz[i] = -2.0 * g_x_xy_z_zz[i] * a_exp + 4.0 * g_x_xy_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_x_0_y_0_0_xy_zz_xx, g_x_0_y_0_0_xy_zz_xy, g_x_0_y_0_0_xy_zz_xz, g_x_0_y_0_0_xy_zz_yy, g_x_0_y_0_0_xy_zz_yz, g_x_0_y_0_0_xy_zz_zz, g_x_xy_yzz_xx, g_x_xy_yzz_xy, g_x_xy_yzz_xz, g_x_xy_yzz_yy, g_x_xy_yzz_yz, g_x_xy_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xy_zz_xx[i] = 4.0 * g_x_xy_yzz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_zz_xy[i] = 4.0 * g_x_xy_yzz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_zz_xz[i] = 4.0 * g_x_xy_yzz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_zz_yy[i] = 4.0 * g_x_xy_yzz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_zz_yz[i] = 4.0 * g_x_xy_yzz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_zz_zz[i] = 4.0 * g_x_xy_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_x_0_y_0_0_xz_xx_xx, g_x_0_y_0_0_xz_xx_xy, g_x_0_y_0_0_xz_xx_xz, g_x_0_y_0_0_xz_xx_yy, g_x_0_y_0_0_xz_xx_yz, g_x_0_y_0_0_xz_xx_zz, g_x_xz_xxy_xx, g_x_xz_xxy_xy, g_x_xz_xxy_xz, g_x_xz_xxy_yy, g_x_xz_xxy_yz, g_x_xz_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xz_xx_xx[i] = 4.0 * g_x_xz_xxy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_xx_xy[i] = 4.0 * g_x_xz_xxy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_xx_xz[i] = 4.0 * g_x_xz_xxy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_xx_yy[i] = 4.0 * g_x_xz_xxy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_xx_yz[i] = 4.0 * g_x_xz_xxy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_xx_zz[i] = 4.0 * g_x_xz_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_x_0_y_0_0_xz_xy_xx, g_x_0_y_0_0_xz_xy_xy, g_x_0_y_0_0_xz_xy_xz, g_x_0_y_0_0_xz_xy_yy, g_x_0_y_0_0_xz_xy_yz, g_x_0_y_0_0_xz_xy_zz, g_x_xz_x_xx, g_x_xz_x_xy, g_x_xz_x_xz, g_x_xz_x_yy, g_x_xz_x_yz, g_x_xz_x_zz, g_x_xz_xyy_xx, g_x_xz_xyy_xy, g_x_xz_xyy_xz, g_x_xz_xyy_yy, g_x_xz_xyy_yz, g_x_xz_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xz_xy_xx[i] = -2.0 * g_x_xz_x_xx[i] * a_exp + 4.0 * g_x_xz_xyy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_xy_xy[i] = -2.0 * g_x_xz_x_xy[i] * a_exp + 4.0 * g_x_xz_xyy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_xy_xz[i] = -2.0 * g_x_xz_x_xz[i] * a_exp + 4.0 * g_x_xz_xyy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_xy_yy[i] = -2.0 * g_x_xz_x_yy[i] * a_exp + 4.0 * g_x_xz_xyy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_xy_yz[i] = -2.0 * g_x_xz_x_yz[i] * a_exp + 4.0 * g_x_xz_xyy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_xy_zz[i] = -2.0 * g_x_xz_x_zz[i] * a_exp + 4.0 * g_x_xz_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_x_0_y_0_0_xz_xz_xx, g_x_0_y_0_0_xz_xz_xy, g_x_0_y_0_0_xz_xz_xz, g_x_0_y_0_0_xz_xz_yy, g_x_0_y_0_0_xz_xz_yz, g_x_0_y_0_0_xz_xz_zz, g_x_xz_xyz_xx, g_x_xz_xyz_xy, g_x_xz_xyz_xz, g_x_xz_xyz_yy, g_x_xz_xyz_yz, g_x_xz_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xz_xz_xx[i] = 4.0 * g_x_xz_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_xz_xy[i] = 4.0 * g_x_xz_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_xz_xz[i] = 4.0 * g_x_xz_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_xz_yy[i] = 4.0 * g_x_xz_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_xz_yz[i] = 4.0 * g_x_xz_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_xz_zz[i] = 4.0 * g_x_xz_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_x_0_y_0_0_xz_yy_xx, g_x_0_y_0_0_xz_yy_xy, g_x_0_y_0_0_xz_yy_xz, g_x_0_y_0_0_xz_yy_yy, g_x_0_y_0_0_xz_yy_yz, g_x_0_y_0_0_xz_yy_zz, g_x_xz_y_xx, g_x_xz_y_xy, g_x_xz_y_xz, g_x_xz_y_yy, g_x_xz_y_yz, g_x_xz_y_zz, g_x_xz_yyy_xx, g_x_xz_yyy_xy, g_x_xz_yyy_xz, g_x_xz_yyy_yy, g_x_xz_yyy_yz, g_x_xz_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xz_yy_xx[i] = -4.0 * g_x_xz_y_xx[i] * a_exp + 4.0 * g_x_xz_yyy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_yy_xy[i] = -4.0 * g_x_xz_y_xy[i] * a_exp + 4.0 * g_x_xz_yyy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_yy_xz[i] = -4.0 * g_x_xz_y_xz[i] * a_exp + 4.0 * g_x_xz_yyy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_yy_yy[i] = -4.0 * g_x_xz_y_yy[i] * a_exp + 4.0 * g_x_xz_yyy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_yy_yz[i] = -4.0 * g_x_xz_y_yz[i] * a_exp + 4.0 * g_x_xz_yyy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_yy_zz[i] = -4.0 * g_x_xz_y_zz[i] * a_exp + 4.0 * g_x_xz_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_x_0_y_0_0_xz_yz_xx, g_x_0_y_0_0_xz_yz_xy, g_x_0_y_0_0_xz_yz_xz, g_x_0_y_0_0_xz_yz_yy, g_x_0_y_0_0_xz_yz_yz, g_x_0_y_0_0_xz_yz_zz, g_x_xz_yyz_xx, g_x_xz_yyz_xy, g_x_xz_yyz_xz, g_x_xz_yyz_yy, g_x_xz_yyz_yz, g_x_xz_yyz_zz, g_x_xz_z_xx, g_x_xz_z_xy, g_x_xz_z_xz, g_x_xz_z_yy, g_x_xz_z_yz, g_x_xz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xz_yz_xx[i] = -2.0 * g_x_xz_z_xx[i] * a_exp + 4.0 * g_x_xz_yyz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_yz_xy[i] = -2.0 * g_x_xz_z_xy[i] * a_exp + 4.0 * g_x_xz_yyz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_yz_xz[i] = -2.0 * g_x_xz_z_xz[i] * a_exp + 4.0 * g_x_xz_yyz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_yz_yy[i] = -2.0 * g_x_xz_z_yy[i] * a_exp + 4.0 * g_x_xz_yyz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_yz_yz[i] = -2.0 * g_x_xz_z_yz[i] * a_exp + 4.0 * g_x_xz_yyz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_yz_zz[i] = -2.0 * g_x_xz_z_zz[i] * a_exp + 4.0 * g_x_xz_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_x_0_y_0_0_xz_zz_xx, g_x_0_y_0_0_xz_zz_xy, g_x_0_y_0_0_xz_zz_xz, g_x_0_y_0_0_xz_zz_yy, g_x_0_y_0_0_xz_zz_yz, g_x_0_y_0_0_xz_zz_zz, g_x_xz_yzz_xx, g_x_xz_yzz_xy, g_x_xz_yzz_xz, g_x_xz_yzz_yy, g_x_xz_yzz_yz, g_x_xz_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xz_zz_xx[i] = 4.0 * g_x_xz_yzz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_zz_xy[i] = 4.0 * g_x_xz_yzz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_zz_xz[i] = 4.0 * g_x_xz_yzz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_zz_yy[i] = 4.0 * g_x_xz_yzz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_zz_yz[i] = 4.0 * g_x_xz_yzz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_zz_zz[i] = 4.0 * g_x_xz_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (324-330)

    #pragma omp simd aligned(g_x_0_y_0_0_yy_xx_xx, g_x_0_y_0_0_yy_xx_xy, g_x_0_y_0_0_yy_xx_xz, g_x_0_y_0_0_yy_xx_yy, g_x_0_y_0_0_yy_xx_yz, g_x_0_y_0_0_yy_xx_zz, g_x_yy_xxy_xx, g_x_yy_xxy_xy, g_x_yy_xxy_xz, g_x_yy_xxy_yy, g_x_yy_xxy_yz, g_x_yy_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_yy_xx_xx[i] = 4.0 * g_x_yy_xxy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_xx_xy[i] = 4.0 * g_x_yy_xxy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_xx_xz[i] = 4.0 * g_x_yy_xxy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_xx_yy[i] = 4.0 * g_x_yy_xxy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_xx_yz[i] = 4.0 * g_x_yy_xxy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_xx_zz[i] = 4.0 * g_x_yy_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (330-336)

    #pragma omp simd aligned(g_x_0_y_0_0_yy_xy_xx, g_x_0_y_0_0_yy_xy_xy, g_x_0_y_0_0_yy_xy_xz, g_x_0_y_0_0_yy_xy_yy, g_x_0_y_0_0_yy_xy_yz, g_x_0_y_0_0_yy_xy_zz, g_x_yy_x_xx, g_x_yy_x_xy, g_x_yy_x_xz, g_x_yy_x_yy, g_x_yy_x_yz, g_x_yy_x_zz, g_x_yy_xyy_xx, g_x_yy_xyy_xy, g_x_yy_xyy_xz, g_x_yy_xyy_yy, g_x_yy_xyy_yz, g_x_yy_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_yy_xy_xx[i] = -2.0 * g_x_yy_x_xx[i] * a_exp + 4.0 * g_x_yy_xyy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_xy_xy[i] = -2.0 * g_x_yy_x_xy[i] * a_exp + 4.0 * g_x_yy_xyy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_xy_xz[i] = -2.0 * g_x_yy_x_xz[i] * a_exp + 4.0 * g_x_yy_xyy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_xy_yy[i] = -2.0 * g_x_yy_x_yy[i] * a_exp + 4.0 * g_x_yy_xyy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_xy_yz[i] = -2.0 * g_x_yy_x_yz[i] * a_exp + 4.0 * g_x_yy_xyy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_xy_zz[i] = -2.0 * g_x_yy_x_zz[i] * a_exp + 4.0 * g_x_yy_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (336-342)

    #pragma omp simd aligned(g_x_0_y_0_0_yy_xz_xx, g_x_0_y_0_0_yy_xz_xy, g_x_0_y_0_0_yy_xz_xz, g_x_0_y_0_0_yy_xz_yy, g_x_0_y_0_0_yy_xz_yz, g_x_0_y_0_0_yy_xz_zz, g_x_yy_xyz_xx, g_x_yy_xyz_xy, g_x_yy_xyz_xz, g_x_yy_xyz_yy, g_x_yy_xyz_yz, g_x_yy_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_yy_xz_xx[i] = 4.0 * g_x_yy_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_xz_xy[i] = 4.0 * g_x_yy_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_xz_xz[i] = 4.0 * g_x_yy_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_xz_yy[i] = 4.0 * g_x_yy_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_xz_yz[i] = 4.0 * g_x_yy_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_xz_zz[i] = 4.0 * g_x_yy_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (342-348)

    #pragma omp simd aligned(g_x_0_y_0_0_yy_yy_xx, g_x_0_y_0_0_yy_yy_xy, g_x_0_y_0_0_yy_yy_xz, g_x_0_y_0_0_yy_yy_yy, g_x_0_y_0_0_yy_yy_yz, g_x_0_y_0_0_yy_yy_zz, g_x_yy_y_xx, g_x_yy_y_xy, g_x_yy_y_xz, g_x_yy_y_yy, g_x_yy_y_yz, g_x_yy_y_zz, g_x_yy_yyy_xx, g_x_yy_yyy_xy, g_x_yy_yyy_xz, g_x_yy_yyy_yy, g_x_yy_yyy_yz, g_x_yy_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_yy_yy_xx[i] = -4.0 * g_x_yy_y_xx[i] * a_exp + 4.0 * g_x_yy_yyy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_yy_xy[i] = -4.0 * g_x_yy_y_xy[i] * a_exp + 4.0 * g_x_yy_yyy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_yy_xz[i] = -4.0 * g_x_yy_y_xz[i] * a_exp + 4.0 * g_x_yy_yyy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_yy_yy[i] = -4.0 * g_x_yy_y_yy[i] * a_exp + 4.0 * g_x_yy_yyy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_yy_yz[i] = -4.0 * g_x_yy_y_yz[i] * a_exp + 4.0 * g_x_yy_yyy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_yy_zz[i] = -4.0 * g_x_yy_y_zz[i] * a_exp + 4.0 * g_x_yy_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (348-354)

    #pragma omp simd aligned(g_x_0_y_0_0_yy_yz_xx, g_x_0_y_0_0_yy_yz_xy, g_x_0_y_0_0_yy_yz_xz, g_x_0_y_0_0_yy_yz_yy, g_x_0_y_0_0_yy_yz_yz, g_x_0_y_0_0_yy_yz_zz, g_x_yy_yyz_xx, g_x_yy_yyz_xy, g_x_yy_yyz_xz, g_x_yy_yyz_yy, g_x_yy_yyz_yz, g_x_yy_yyz_zz, g_x_yy_z_xx, g_x_yy_z_xy, g_x_yy_z_xz, g_x_yy_z_yy, g_x_yy_z_yz, g_x_yy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_yy_yz_xx[i] = -2.0 * g_x_yy_z_xx[i] * a_exp + 4.0 * g_x_yy_yyz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_yz_xy[i] = -2.0 * g_x_yy_z_xy[i] * a_exp + 4.0 * g_x_yy_yyz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_yz_xz[i] = -2.0 * g_x_yy_z_xz[i] * a_exp + 4.0 * g_x_yy_yyz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_yz_yy[i] = -2.0 * g_x_yy_z_yy[i] * a_exp + 4.0 * g_x_yy_yyz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_yz_yz[i] = -2.0 * g_x_yy_z_yz[i] * a_exp + 4.0 * g_x_yy_yyz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_yz_zz[i] = -2.0 * g_x_yy_z_zz[i] * a_exp + 4.0 * g_x_yy_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (354-360)

    #pragma omp simd aligned(g_x_0_y_0_0_yy_zz_xx, g_x_0_y_0_0_yy_zz_xy, g_x_0_y_0_0_yy_zz_xz, g_x_0_y_0_0_yy_zz_yy, g_x_0_y_0_0_yy_zz_yz, g_x_0_y_0_0_yy_zz_zz, g_x_yy_yzz_xx, g_x_yy_yzz_xy, g_x_yy_yzz_xz, g_x_yy_yzz_yy, g_x_yy_yzz_yz, g_x_yy_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_yy_zz_xx[i] = 4.0 * g_x_yy_yzz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_zz_xy[i] = 4.0 * g_x_yy_yzz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_zz_xz[i] = 4.0 * g_x_yy_yzz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_zz_yy[i] = 4.0 * g_x_yy_yzz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_zz_yz[i] = 4.0 * g_x_yy_yzz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_zz_zz[i] = 4.0 * g_x_yy_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (360-366)

    #pragma omp simd aligned(g_x_0_y_0_0_yz_xx_xx, g_x_0_y_0_0_yz_xx_xy, g_x_0_y_0_0_yz_xx_xz, g_x_0_y_0_0_yz_xx_yy, g_x_0_y_0_0_yz_xx_yz, g_x_0_y_0_0_yz_xx_zz, g_x_yz_xxy_xx, g_x_yz_xxy_xy, g_x_yz_xxy_xz, g_x_yz_xxy_yy, g_x_yz_xxy_yz, g_x_yz_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_yz_xx_xx[i] = 4.0 * g_x_yz_xxy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_xx_xy[i] = 4.0 * g_x_yz_xxy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_xx_xz[i] = 4.0 * g_x_yz_xxy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_xx_yy[i] = 4.0 * g_x_yz_xxy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_xx_yz[i] = 4.0 * g_x_yz_xxy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_xx_zz[i] = 4.0 * g_x_yz_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (366-372)

    #pragma omp simd aligned(g_x_0_y_0_0_yz_xy_xx, g_x_0_y_0_0_yz_xy_xy, g_x_0_y_0_0_yz_xy_xz, g_x_0_y_0_0_yz_xy_yy, g_x_0_y_0_0_yz_xy_yz, g_x_0_y_0_0_yz_xy_zz, g_x_yz_x_xx, g_x_yz_x_xy, g_x_yz_x_xz, g_x_yz_x_yy, g_x_yz_x_yz, g_x_yz_x_zz, g_x_yz_xyy_xx, g_x_yz_xyy_xy, g_x_yz_xyy_xz, g_x_yz_xyy_yy, g_x_yz_xyy_yz, g_x_yz_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_yz_xy_xx[i] = -2.0 * g_x_yz_x_xx[i] * a_exp + 4.0 * g_x_yz_xyy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_xy_xy[i] = -2.0 * g_x_yz_x_xy[i] * a_exp + 4.0 * g_x_yz_xyy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_xy_xz[i] = -2.0 * g_x_yz_x_xz[i] * a_exp + 4.0 * g_x_yz_xyy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_xy_yy[i] = -2.0 * g_x_yz_x_yy[i] * a_exp + 4.0 * g_x_yz_xyy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_xy_yz[i] = -2.0 * g_x_yz_x_yz[i] * a_exp + 4.0 * g_x_yz_xyy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_xy_zz[i] = -2.0 * g_x_yz_x_zz[i] * a_exp + 4.0 * g_x_yz_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (372-378)

    #pragma omp simd aligned(g_x_0_y_0_0_yz_xz_xx, g_x_0_y_0_0_yz_xz_xy, g_x_0_y_0_0_yz_xz_xz, g_x_0_y_0_0_yz_xz_yy, g_x_0_y_0_0_yz_xz_yz, g_x_0_y_0_0_yz_xz_zz, g_x_yz_xyz_xx, g_x_yz_xyz_xy, g_x_yz_xyz_xz, g_x_yz_xyz_yy, g_x_yz_xyz_yz, g_x_yz_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_yz_xz_xx[i] = 4.0 * g_x_yz_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_xz_xy[i] = 4.0 * g_x_yz_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_xz_xz[i] = 4.0 * g_x_yz_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_xz_yy[i] = 4.0 * g_x_yz_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_xz_yz[i] = 4.0 * g_x_yz_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_xz_zz[i] = 4.0 * g_x_yz_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (378-384)

    #pragma omp simd aligned(g_x_0_y_0_0_yz_yy_xx, g_x_0_y_0_0_yz_yy_xy, g_x_0_y_0_0_yz_yy_xz, g_x_0_y_0_0_yz_yy_yy, g_x_0_y_0_0_yz_yy_yz, g_x_0_y_0_0_yz_yy_zz, g_x_yz_y_xx, g_x_yz_y_xy, g_x_yz_y_xz, g_x_yz_y_yy, g_x_yz_y_yz, g_x_yz_y_zz, g_x_yz_yyy_xx, g_x_yz_yyy_xy, g_x_yz_yyy_xz, g_x_yz_yyy_yy, g_x_yz_yyy_yz, g_x_yz_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_yz_yy_xx[i] = -4.0 * g_x_yz_y_xx[i] * a_exp + 4.0 * g_x_yz_yyy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_yy_xy[i] = -4.0 * g_x_yz_y_xy[i] * a_exp + 4.0 * g_x_yz_yyy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_yy_xz[i] = -4.0 * g_x_yz_y_xz[i] * a_exp + 4.0 * g_x_yz_yyy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_yy_yy[i] = -4.0 * g_x_yz_y_yy[i] * a_exp + 4.0 * g_x_yz_yyy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_yy_yz[i] = -4.0 * g_x_yz_y_yz[i] * a_exp + 4.0 * g_x_yz_yyy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_yy_zz[i] = -4.0 * g_x_yz_y_zz[i] * a_exp + 4.0 * g_x_yz_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (384-390)

    #pragma omp simd aligned(g_x_0_y_0_0_yz_yz_xx, g_x_0_y_0_0_yz_yz_xy, g_x_0_y_0_0_yz_yz_xz, g_x_0_y_0_0_yz_yz_yy, g_x_0_y_0_0_yz_yz_yz, g_x_0_y_0_0_yz_yz_zz, g_x_yz_yyz_xx, g_x_yz_yyz_xy, g_x_yz_yyz_xz, g_x_yz_yyz_yy, g_x_yz_yyz_yz, g_x_yz_yyz_zz, g_x_yz_z_xx, g_x_yz_z_xy, g_x_yz_z_xz, g_x_yz_z_yy, g_x_yz_z_yz, g_x_yz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_yz_yz_xx[i] = -2.0 * g_x_yz_z_xx[i] * a_exp + 4.0 * g_x_yz_yyz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_yz_xy[i] = -2.0 * g_x_yz_z_xy[i] * a_exp + 4.0 * g_x_yz_yyz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_yz_xz[i] = -2.0 * g_x_yz_z_xz[i] * a_exp + 4.0 * g_x_yz_yyz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_yz_yy[i] = -2.0 * g_x_yz_z_yy[i] * a_exp + 4.0 * g_x_yz_yyz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_yz_yz[i] = -2.0 * g_x_yz_z_yz[i] * a_exp + 4.0 * g_x_yz_yyz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_yz_zz[i] = -2.0 * g_x_yz_z_zz[i] * a_exp + 4.0 * g_x_yz_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (390-396)

    #pragma omp simd aligned(g_x_0_y_0_0_yz_zz_xx, g_x_0_y_0_0_yz_zz_xy, g_x_0_y_0_0_yz_zz_xz, g_x_0_y_0_0_yz_zz_yy, g_x_0_y_0_0_yz_zz_yz, g_x_0_y_0_0_yz_zz_zz, g_x_yz_yzz_xx, g_x_yz_yzz_xy, g_x_yz_yzz_xz, g_x_yz_yzz_yy, g_x_yz_yzz_yz, g_x_yz_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_yz_zz_xx[i] = 4.0 * g_x_yz_yzz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_zz_xy[i] = 4.0 * g_x_yz_yzz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_zz_xz[i] = 4.0 * g_x_yz_yzz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_zz_yy[i] = 4.0 * g_x_yz_yzz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_zz_yz[i] = 4.0 * g_x_yz_yzz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_zz_zz[i] = 4.0 * g_x_yz_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (396-402)

    #pragma omp simd aligned(g_x_0_y_0_0_zz_xx_xx, g_x_0_y_0_0_zz_xx_xy, g_x_0_y_0_0_zz_xx_xz, g_x_0_y_0_0_zz_xx_yy, g_x_0_y_0_0_zz_xx_yz, g_x_0_y_0_0_zz_xx_zz, g_x_zz_xxy_xx, g_x_zz_xxy_xy, g_x_zz_xxy_xz, g_x_zz_xxy_yy, g_x_zz_xxy_yz, g_x_zz_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_zz_xx_xx[i] = 4.0 * g_x_zz_xxy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_xx_xy[i] = 4.0 * g_x_zz_xxy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_xx_xz[i] = 4.0 * g_x_zz_xxy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_xx_yy[i] = 4.0 * g_x_zz_xxy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_xx_yz[i] = 4.0 * g_x_zz_xxy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_xx_zz[i] = 4.0 * g_x_zz_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (402-408)

    #pragma omp simd aligned(g_x_0_y_0_0_zz_xy_xx, g_x_0_y_0_0_zz_xy_xy, g_x_0_y_0_0_zz_xy_xz, g_x_0_y_0_0_zz_xy_yy, g_x_0_y_0_0_zz_xy_yz, g_x_0_y_0_0_zz_xy_zz, g_x_zz_x_xx, g_x_zz_x_xy, g_x_zz_x_xz, g_x_zz_x_yy, g_x_zz_x_yz, g_x_zz_x_zz, g_x_zz_xyy_xx, g_x_zz_xyy_xy, g_x_zz_xyy_xz, g_x_zz_xyy_yy, g_x_zz_xyy_yz, g_x_zz_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_zz_xy_xx[i] = -2.0 * g_x_zz_x_xx[i] * a_exp + 4.0 * g_x_zz_xyy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_xy_xy[i] = -2.0 * g_x_zz_x_xy[i] * a_exp + 4.0 * g_x_zz_xyy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_xy_xz[i] = -2.0 * g_x_zz_x_xz[i] * a_exp + 4.0 * g_x_zz_xyy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_xy_yy[i] = -2.0 * g_x_zz_x_yy[i] * a_exp + 4.0 * g_x_zz_xyy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_xy_yz[i] = -2.0 * g_x_zz_x_yz[i] * a_exp + 4.0 * g_x_zz_xyy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_xy_zz[i] = -2.0 * g_x_zz_x_zz[i] * a_exp + 4.0 * g_x_zz_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (408-414)

    #pragma omp simd aligned(g_x_0_y_0_0_zz_xz_xx, g_x_0_y_0_0_zz_xz_xy, g_x_0_y_0_0_zz_xz_xz, g_x_0_y_0_0_zz_xz_yy, g_x_0_y_0_0_zz_xz_yz, g_x_0_y_0_0_zz_xz_zz, g_x_zz_xyz_xx, g_x_zz_xyz_xy, g_x_zz_xyz_xz, g_x_zz_xyz_yy, g_x_zz_xyz_yz, g_x_zz_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_zz_xz_xx[i] = 4.0 * g_x_zz_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_xz_xy[i] = 4.0 * g_x_zz_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_xz_xz[i] = 4.0 * g_x_zz_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_xz_yy[i] = 4.0 * g_x_zz_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_xz_yz[i] = 4.0 * g_x_zz_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_xz_zz[i] = 4.0 * g_x_zz_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (414-420)

    #pragma omp simd aligned(g_x_0_y_0_0_zz_yy_xx, g_x_0_y_0_0_zz_yy_xy, g_x_0_y_0_0_zz_yy_xz, g_x_0_y_0_0_zz_yy_yy, g_x_0_y_0_0_zz_yy_yz, g_x_0_y_0_0_zz_yy_zz, g_x_zz_y_xx, g_x_zz_y_xy, g_x_zz_y_xz, g_x_zz_y_yy, g_x_zz_y_yz, g_x_zz_y_zz, g_x_zz_yyy_xx, g_x_zz_yyy_xy, g_x_zz_yyy_xz, g_x_zz_yyy_yy, g_x_zz_yyy_yz, g_x_zz_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_zz_yy_xx[i] = -4.0 * g_x_zz_y_xx[i] * a_exp + 4.0 * g_x_zz_yyy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_yy_xy[i] = -4.0 * g_x_zz_y_xy[i] * a_exp + 4.0 * g_x_zz_yyy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_yy_xz[i] = -4.0 * g_x_zz_y_xz[i] * a_exp + 4.0 * g_x_zz_yyy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_yy_yy[i] = -4.0 * g_x_zz_y_yy[i] * a_exp + 4.0 * g_x_zz_yyy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_yy_yz[i] = -4.0 * g_x_zz_y_yz[i] * a_exp + 4.0 * g_x_zz_yyy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_yy_zz[i] = -4.0 * g_x_zz_y_zz[i] * a_exp + 4.0 * g_x_zz_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (420-426)

    #pragma omp simd aligned(g_x_0_y_0_0_zz_yz_xx, g_x_0_y_0_0_zz_yz_xy, g_x_0_y_0_0_zz_yz_xz, g_x_0_y_0_0_zz_yz_yy, g_x_0_y_0_0_zz_yz_yz, g_x_0_y_0_0_zz_yz_zz, g_x_zz_yyz_xx, g_x_zz_yyz_xy, g_x_zz_yyz_xz, g_x_zz_yyz_yy, g_x_zz_yyz_yz, g_x_zz_yyz_zz, g_x_zz_z_xx, g_x_zz_z_xy, g_x_zz_z_xz, g_x_zz_z_yy, g_x_zz_z_yz, g_x_zz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_zz_yz_xx[i] = -2.0 * g_x_zz_z_xx[i] * a_exp + 4.0 * g_x_zz_yyz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_yz_xy[i] = -2.0 * g_x_zz_z_xy[i] * a_exp + 4.0 * g_x_zz_yyz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_yz_xz[i] = -2.0 * g_x_zz_z_xz[i] * a_exp + 4.0 * g_x_zz_yyz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_yz_yy[i] = -2.0 * g_x_zz_z_yy[i] * a_exp + 4.0 * g_x_zz_yyz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_yz_yz[i] = -2.0 * g_x_zz_z_yz[i] * a_exp + 4.0 * g_x_zz_yyz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_yz_zz[i] = -2.0 * g_x_zz_z_zz[i] * a_exp + 4.0 * g_x_zz_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (426-432)

    #pragma omp simd aligned(g_x_0_y_0_0_zz_zz_xx, g_x_0_y_0_0_zz_zz_xy, g_x_0_y_0_0_zz_zz_xz, g_x_0_y_0_0_zz_zz_yy, g_x_0_y_0_0_zz_zz_yz, g_x_0_y_0_0_zz_zz_zz, g_x_zz_yzz_xx, g_x_zz_yzz_xy, g_x_zz_yzz_xz, g_x_zz_yzz_yy, g_x_zz_yzz_yz, g_x_zz_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_zz_zz_xx[i] = 4.0 * g_x_zz_yzz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_zz_xy[i] = 4.0 * g_x_zz_yzz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_zz_xz[i] = 4.0 * g_x_zz_yzz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_zz_yy[i] = 4.0 * g_x_zz_yzz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_zz_yz[i] = 4.0 * g_x_zz_yzz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_zz_zz[i] = 4.0 * g_x_zz_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (432-438)

    #pragma omp simd aligned(g_x_0_z_0_0_xx_xx_xx, g_x_0_z_0_0_xx_xx_xy, g_x_0_z_0_0_xx_xx_xz, g_x_0_z_0_0_xx_xx_yy, g_x_0_z_0_0_xx_xx_yz, g_x_0_z_0_0_xx_xx_zz, g_x_xx_xxz_xx, g_x_xx_xxz_xy, g_x_xx_xxz_xz, g_x_xx_xxz_yy, g_x_xx_xxz_yz, g_x_xx_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xx_xx_xx[i] = 4.0 * g_x_xx_xxz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_xx_xy[i] = 4.0 * g_x_xx_xxz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_xx_xz[i] = 4.0 * g_x_xx_xxz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_xx_yy[i] = 4.0 * g_x_xx_xxz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_xx_yz[i] = 4.0 * g_x_xx_xxz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_xx_zz[i] = 4.0 * g_x_xx_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (438-444)

    #pragma omp simd aligned(g_x_0_z_0_0_xx_xy_xx, g_x_0_z_0_0_xx_xy_xy, g_x_0_z_0_0_xx_xy_xz, g_x_0_z_0_0_xx_xy_yy, g_x_0_z_0_0_xx_xy_yz, g_x_0_z_0_0_xx_xy_zz, g_x_xx_xyz_xx, g_x_xx_xyz_xy, g_x_xx_xyz_xz, g_x_xx_xyz_yy, g_x_xx_xyz_yz, g_x_xx_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xx_xy_xx[i] = 4.0 * g_x_xx_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_xy_xy[i] = 4.0 * g_x_xx_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_xy_xz[i] = 4.0 * g_x_xx_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_xy_yy[i] = 4.0 * g_x_xx_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_xy_yz[i] = 4.0 * g_x_xx_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_xy_zz[i] = 4.0 * g_x_xx_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (444-450)

    #pragma omp simd aligned(g_x_0_z_0_0_xx_xz_xx, g_x_0_z_0_0_xx_xz_xy, g_x_0_z_0_0_xx_xz_xz, g_x_0_z_0_0_xx_xz_yy, g_x_0_z_0_0_xx_xz_yz, g_x_0_z_0_0_xx_xz_zz, g_x_xx_x_xx, g_x_xx_x_xy, g_x_xx_x_xz, g_x_xx_x_yy, g_x_xx_x_yz, g_x_xx_x_zz, g_x_xx_xzz_xx, g_x_xx_xzz_xy, g_x_xx_xzz_xz, g_x_xx_xzz_yy, g_x_xx_xzz_yz, g_x_xx_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xx_xz_xx[i] = -2.0 * g_x_xx_x_xx[i] * a_exp + 4.0 * g_x_xx_xzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_xz_xy[i] = -2.0 * g_x_xx_x_xy[i] * a_exp + 4.0 * g_x_xx_xzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_xz_xz[i] = -2.0 * g_x_xx_x_xz[i] * a_exp + 4.0 * g_x_xx_xzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_xz_yy[i] = -2.0 * g_x_xx_x_yy[i] * a_exp + 4.0 * g_x_xx_xzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_xz_yz[i] = -2.0 * g_x_xx_x_yz[i] * a_exp + 4.0 * g_x_xx_xzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_xz_zz[i] = -2.0 * g_x_xx_x_zz[i] * a_exp + 4.0 * g_x_xx_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (450-456)

    #pragma omp simd aligned(g_x_0_z_0_0_xx_yy_xx, g_x_0_z_0_0_xx_yy_xy, g_x_0_z_0_0_xx_yy_xz, g_x_0_z_0_0_xx_yy_yy, g_x_0_z_0_0_xx_yy_yz, g_x_0_z_0_0_xx_yy_zz, g_x_xx_yyz_xx, g_x_xx_yyz_xy, g_x_xx_yyz_xz, g_x_xx_yyz_yy, g_x_xx_yyz_yz, g_x_xx_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xx_yy_xx[i] = 4.0 * g_x_xx_yyz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_yy_xy[i] = 4.0 * g_x_xx_yyz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_yy_xz[i] = 4.0 * g_x_xx_yyz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_yy_yy[i] = 4.0 * g_x_xx_yyz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_yy_yz[i] = 4.0 * g_x_xx_yyz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_yy_zz[i] = 4.0 * g_x_xx_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (456-462)

    #pragma omp simd aligned(g_x_0_z_0_0_xx_yz_xx, g_x_0_z_0_0_xx_yz_xy, g_x_0_z_0_0_xx_yz_xz, g_x_0_z_0_0_xx_yz_yy, g_x_0_z_0_0_xx_yz_yz, g_x_0_z_0_0_xx_yz_zz, g_x_xx_y_xx, g_x_xx_y_xy, g_x_xx_y_xz, g_x_xx_y_yy, g_x_xx_y_yz, g_x_xx_y_zz, g_x_xx_yzz_xx, g_x_xx_yzz_xy, g_x_xx_yzz_xz, g_x_xx_yzz_yy, g_x_xx_yzz_yz, g_x_xx_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xx_yz_xx[i] = -2.0 * g_x_xx_y_xx[i] * a_exp + 4.0 * g_x_xx_yzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_yz_xy[i] = -2.0 * g_x_xx_y_xy[i] * a_exp + 4.0 * g_x_xx_yzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_yz_xz[i] = -2.0 * g_x_xx_y_xz[i] * a_exp + 4.0 * g_x_xx_yzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_yz_yy[i] = -2.0 * g_x_xx_y_yy[i] * a_exp + 4.0 * g_x_xx_yzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_yz_yz[i] = -2.0 * g_x_xx_y_yz[i] * a_exp + 4.0 * g_x_xx_yzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_yz_zz[i] = -2.0 * g_x_xx_y_zz[i] * a_exp + 4.0 * g_x_xx_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (462-468)

    #pragma omp simd aligned(g_x_0_z_0_0_xx_zz_xx, g_x_0_z_0_0_xx_zz_xy, g_x_0_z_0_0_xx_zz_xz, g_x_0_z_0_0_xx_zz_yy, g_x_0_z_0_0_xx_zz_yz, g_x_0_z_0_0_xx_zz_zz, g_x_xx_z_xx, g_x_xx_z_xy, g_x_xx_z_xz, g_x_xx_z_yy, g_x_xx_z_yz, g_x_xx_z_zz, g_x_xx_zzz_xx, g_x_xx_zzz_xy, g_x_xx_zzz_xz, g_x_xx_zzz_yy, g_x_xx_zzz_yz, g_x_xx_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xx_zz_xx[i] = -4.0 * g_x_xx_z_xx[i] * a_exp + 4.0 * g_x_xx_zzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_zz_xy[i] = -4.0 * g_x_xx_z_xy[i] * a_exp + 4.0 * g_x_xx_zzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_zz_xz[i] = -4.0 * g_x_xx_z_xz[i] * a_exp + 4.0 * g_x_xx_zzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_zz_yy[i] = -4.0 * g_x_xx_z_yy[i] * a_exp + 4.0 * g_x_xx_zzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_zz_yz[i] = -4.0 * g_x_xx_z_yz[i] * a_exp + 4.0 * g_x_xx_zzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_zz_zz[i] = -4.0 * g_x_xx_z_zz[i] * a_exp + 4.0 * g_x_xx_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (468-474)

    #pragma omp simd aligned(g_x_0_z_0_0_xy_xx_xx, g_x_0_z_0_0_xy_xx_xy, g_x_0_z_0_0_xy_xx_xz, g_x_0_z_0_0_xy_xx_yy, g_x_0_z_0_0_xy_xx_yz, g_x_0_z_0_0_xy_xx_zz, g_x_xy_xxz_xx, g_x_xy_xxz_xy, g_x_xy_xxz_xz, g_x_xy_xxz_yy, g_x_xy_xxz_yz, g_x_xy_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xy_xx_xx[i] = 4.0 * g_x_xy_xxz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_xx_xy[i] = 4.0 * g_x_xy_xxz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_xx_xz[i] = 4.0 * g_x_xy_xxz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_xx_yy[i] = 4.0 * g_x_xy_xxz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_xx_yz[i] = 4.0 * g_x_xy_xxz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_xx_zz[i] = 4.0 * g_x_xy_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (474-480)

    #pragma omp simd aligned(g_x_0_z_0_0_xy_xy_xx, g_x_0_z_0_0_xy_xy_xy, g_x_0_z_0_0_xy_xy_xz, g_x_0_z_0_0_xy_xy_yy, g_x_0_z_0_0_xy_xy_yz, g_x_0_z_0_0_xy_xy_zz, g_x_xy_xyz_xx, g_x_xy_xyz_xy, g_x_xy_xyz_xz, g_x_xy_xyz_yy, g_x_xy_xyz_yz, g_x_xy_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xy_xy_xx[i] = 4.0 * g_x_xy_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_xy_xy[i] = 4.0 * g_x_xy_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_xy_xz[i] = 4.0 * g_x_xy_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_xy_yy[i] = 4.0 * g_x_xy_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_xy_yz[i] = 4.0 * g_x_xy_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_xy_zz[i] = 4.0 * g_x_xy_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (480-486)

    #pragma omp simd aligned(g_x_0_z_0_0_xy_xz_xx, g_x_0_z_0_0_xy_xz_xy, g_x_0_z_0_0_xy_xz_xz, g_x_0_z_0_0_xy_xz_yy, g_x_0_z_0_0_xy_xz_yz, g_x_0_z_0_0_xy_xz_zz, g_x_xy_x_xx, g_x_xy_x_xy, g_x_xy_x_xz, g_x_xy_x_yy, g_x_xy_x_yz, g_x_xy_x_zz, g_x_xy_xzz_xx, g_x_xy_xzz_xy, g_x_xy_xzz_xz, g_x_xy_xzz_yy, g_x_xy_xzz_yz, g_x_xy_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xy_xz_xx[i] = -2.0 * g_x_xy_x_xx[i] * a_exp + 4.0 * g_x_xy_xzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_xz_xy[i] = -2.0 * g_x_xy_x_xy[i] * a_exp + 4.0 * g_x_xy_xzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_xz_xz[i] = -2.0 * g_x_xy_x_xz[i] * a_exp + 4.0 * g_x_xy_xzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_xz_yy[i] = -2.0 * g_x_xy_x_yy[i] * a_exp + 4.0 * g_x_xy_xzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_xz_yz[i] = -2.0 * g_x_xy_x_yz[i] * a_exp + 4.0 * g_x_xy_xzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_xz_zz[i] = -2.0 * g_x_xy_x_zz[i] * a_exp + 4.0 * g_x_xy_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (486-492)

    #pragma omp simd aligned(g_x_0_z_0_0_xy_yy_xx, g_x_0_z_0_0_xy_yy_xy, g_x_0_z_0_0_xy_yy_xz, g_x_0_z_0_0_xy_yy_yy, g_x_0_z_0_0_xy_yy_yz, g_x_0_z_0_0_xy_yy_zz, g_x_xy_yyz_xx, g_x_xy_yyz_xy, g_x_xy_yyz_xz, g_x_xy_yyz_yy, g_x_xy_yyz_yz, g_x_xy_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xy_yy_xx[i] = 4.0 * g_x_xy_yyz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_yy_xy[i] = 4.0 * g_x_xy_yyz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_yy_xz[i] = 4.0 * g_x_xy_yyz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_yy_yy[i] = 4.0 * g_x_xy_yyz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_yy_yz[i] = 4.0 * g_x_xy_yyz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_yy_zz[i] = 4.0 * g_x_xy_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (492-498)

    #pragma omp simd aligned(g_x_0_z_0_0_xy_yz_xx, g_x_0_z_0_0_xy_yz_xy, g_x_0_z_0_0_xy_yz_xz, g_x_0_z_0_0_xy_yz_yy, g_x_0_z_0_0_xy_yz_yz, g_x_0_z_0_0_xy_yz_zz, g_x_xy_y_xx, g_x_xy_y_xy, g_x_xy_y_xz, g_x_xy_y_yy, g_x_xy_y_yz, g_x_xy_y_zz, g_x_xy_yzz_xx, g_x_xy_yzz_xy, g_x_xy_yzz_xz, g_x_xy_yzz_yy, g_x_xy_yzz_yz, g_x_xy_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xy_yz_xx[i] = -2.0 * g_x_xy_y_xx[i] * a_exp + 4.0 * g_x_xy_yzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_yz_xy[i] = -2.0 * g_x_xy_y_xy[i] * a_exp + 4.0 * g_x_xy_yzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_yz_xz[i] = -2.0 * g_x_xy_y_xz[i] * a_exp + 4.0 * g_x_xy_yzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_yz_yy[i] = -2.0 * g_x_xy_y_yy[i] * a_exp + 4.0 * g_x_xy_yzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_yz_yz[i] = -2.0 * g_x_xy_y_yz[i] * a_exp + 4.0 * g_x_xy_yzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_yz_zz[i] = -2.0 * g_x_xy_y_zz[i] * a_exp + 4.0 * g_x_xy_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (498-504)

    #pragma omp simd aligned(g_x_0_z_0_0_xy_zz_xx, g_x_0_z_0_0_xy_zz_xy, g_x_0_z_0_0_xy_zz_xz, g_x_0_z_0_0_xy_zz_yy, g_x_0_z_0_0_xy_zz_yz, g_x_0_z_0_0_xy_zz_zz, g_x_xy_z_xx, g_x_xy_z_xy, g_x_xy_z_xz, g_x_xy_z_yy, g_x_xy_z_yz, g_x_xy_z_zz, g_x_xy_zzz_xx, g_x_xy_zzz_xy, g_x_xy_zzz_xz, g_x_xy_zzz_yy, g_x_xy_zzz_yz, g_x_xy_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xy_zz_xx[i] = -4.0 * g_x_xy_z_xx[i] * a_exp + 4.0 * g_x_xy_zzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_zz_xy[i] = -4.0 * g_x_xy_z_xy[i] * a_exp + 4.0 * g_x_xy_zzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_zz_xz[i] = -4.0 * g_x_xy_z_xz[i] * a_exp + 4.0 * g_x_xy_zzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_zz_yy[i] = -4.0 * g_x_xy_z_yy[i] * a_exp + 4.0 * g_x_xy_zzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_zz_yz[i] = -4.0 * g_x_xy_z_yz[i] * a_exp + 4.0 * g_x_xy_zzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_zz_zz[i] = -4.0 * g_x_xy_z_zz[i] * a_exp + 4.0 * g_x_xy_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (504-510)

    #pragma omp simd aligned(g_x_0_z_0_0_xz_xx_xx, g_x_0_z_0_0_xz_xx_xy, g_x_0_z_0_0_xz_xx_xz, g_x_0_z_0_0_xz_xx_yy, g_x_0_z_0_0_xz_xx_yz, g_x_0_z_0_0_xz_xx_zz, g_x_xz_xxz_xx, g_x_xz_xxz_xy, g_x_xz_xxz_xz, g_x_xz_xxz_yy, g_x_xz_xxz_yz, g_x_xz_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xz_xx_xx[i] = 4.0 * g_x_xz_xxz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_xx_xy[i] = 4.0 * g_x_xz_xxz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_xx_xz[i] = 4.0 * g_x_xz_xxz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_xx_yy[i] = 4.0 * g_x_xz_xxz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_xx_yz[i] = 4.0 * g_x_xz_xxz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_xx_zz[i] = 4.0 * g_x_xz_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (510-516)

    #pragma omp simd aligned(g_x_0_z_0_0_xz_xy_xx, g_x_0_z_0_0_xz_xy_xy, g_x_0_z_0_0_xz_xy_xz, g_x_0_z_0_0_xz_xy_yy, g_x_0_z_0_0_xz_xy_yz, g_x_0_z_0_0_xz_xy_zz, g_x_xz_xyz_xx, g_x_xz_xyz_xy, g_x_xz_xyz_xz, g_x_xz_xyz_yy, g_x_xz_xyz_yz, g_x_xz_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xz_xy_xx[i] = 4.0 * g_x_xz_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_xy_xy[i] = 4.0 * g_x_xz_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_xy_xz[i] = 4.0 * g_x_xz_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_xy_yy[i] = 4.0 * g_x_xz_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_xy_yz[i] = 4.0 * g_x_xz_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_xy_zz[i] = 4.0 * g_x_xz_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (516-522)

    #pragma omp simd aligned(g_x_0_z_0_0_xz_xz_xx, g_x_0_z_0_0_xz_xz_xy, g_x_0_z_0_0_xz_xz_xz, g_x_0_z_0_0_xz_xz_yy, g_x_0_z_0_0_xz_xz_yz, g_x_0_z_0_0_xz_xz_zz, g_x_xz_x_xx, g_x_xz_x_xy, g_x_xz_x_xz, g_x_xz_x_yy, g_x_xz_x_yz, g_x_xz_x_zz, g_x_xz_xzz_xx, g_x_xz_xzz_xy, g_x_xz_xzz_xz, g_x_xz_xzz_yy, g_x_xz_xzz_yz, g_x_xz_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xz_xz_xx[i] = -2.0 * g_x_xz_x_xx[i] * a_exp + 4.0 * g_x_xz_xzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_xz_xy[i] = -2.0 * g_x_xz_x_xy[i] * a_exp + 4.0 * g_x_xz_xzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_xz_xz[i] = -2.0 * g_x_xz_x_xz[i] * a_exp + 4.0 * g_x_xz_xzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_xz_yy[i] = -2.0 * g_x_xz_x_yy[i] * a_exp + 4.0 * g_x_xz_xzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_xz_yz[i] = -2.0 * g_x_xz_x_yz[i] * a_exp + 4.0 * g_x_xz_xzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_xz_zz[i] = -2.0 * g_x_xz_x_zz[i] * a_exp + 4.0 * g_x_xz_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (522-528)

    #pragma omp simd aligned(g_x_0_z_0_0_xz_yy_xx, g_x_0_z_0_0_xz_yy_xy, g_x_0_z_0_0_xz_yy_xz, g_x_0_z_0_0_xz_yy_yy, g_x_0_z_0_0_xz_yy_yz, g_x_0_z_0_0_xz_yy_zz, g_x_xz_yyz_xx, g_x_xz_yyz_xy, g_x_xz_yyz_xz, g_x_xz_yyz_yy, g_x_xz_yyz_yz, g_x_xz_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xz_yy_xx[i] = 4.0 * g_x_xz_yyz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_yy_xy[i] = 4.0 * g_x_xz_yyz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_yy_xz[i] = 4.0 * g_x_xz_yyz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_yy_yy[i] = 4.0 * g_x_xz_yyz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_yy_yz[i] = 4.0 * g_x_xz_yyz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_yy_zz[i] = 4.0 * g_x_xz_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (528-534)

    #pragma omp simd aligned(g_x_0_z_0_0_xz_yz_xx, g_x_0_z_0_0_xz_yz_xy, g_x_0_z_0_0_xz_yz_xz, g_x_0_z_0_0_xz_yz_yy, g_x_0_z_0_0_xz_yz_yz, g_x_0_z_0_0_xz_yz_zz, g_x_xz_y_xx, g_x_xz_y_xy, g_x_xz_y_xz, g_x_xz_y_yy, g_x_xz_y_yz, g_x_xz_y_zz, g_x_xz_yzz_xx, g_x_xz_yzz_xy, g_x_xz_yzz_xz, g_x_xz_yzz_yy, g_x_xz_yzz_yz, g_x_xz_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xz_yz_xx[i] = -2.0 * g_x_xz_y_xx[i] * a_exp + 4.0 * g_x_xz_yzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_yz_xy[i] = -2.0 * g_x_xz_y_xy[i] * a_exp + 4.0 * g_x_xz_yzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_yz_xz[i] = -2.0 * g_x_xz_y_xz[i] * a_exp + 4.0 * g_x_xz_yzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_yz_yy[i] = -2.0 * g_x_xz_y_yy[i] * a_exp + 4.0 * g_x_xz_yzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_yz_yz[i] = -2.0 * g_x_xz_y_yz[i] * a_exp + 4.0 * g_x_xz_yzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_yz_zz[i] = -2.0 * g_x_xz_y_zz[i] * a_exp + 4.0 * g_x_xz_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (534-540)

    #pragma omp simd aligned(g_x_0_z_0_0_xz_zz_xx, g_x_0_z_0_0_xz_zz_xy, g_x_0_z_0_0_xz_zz_xz, g_x_0_z_0_0_xz_zz_yy, g_x_0_z_0_0_xz_zz_yz, g_x_0_z_0_0_xz_zz_zz, g_x_xz_z_xx, g_x_xz_z_xy, g_x_xz_z_xz, g_x_xz_z_yy, g_x_xz_z_yz, g_x_xz_z_zz, g_x_xz_zzz_xx, g_x_xz_zzz_xy, g_x_xz_zzz_xz, g_x_xz_zzz_yy, g_x_xz_zzz_yz, g_x_xz_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xz_zz_xx[i] = -4.0 * g_x_xz_z_xx[i] * a_exp + 4.0 * g_x_xz_zzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_zz_xy[i] = -4.0 * g_x_xz_z_xy[i] * a_exp + 4.0 * g_x_xz_zzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_zz_xz[i] = -4.0 * g_x_xz_z_xz[i] * a_exp + 4.0 * g_x_xz_zzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_zz_yy[i] = -4.0 * g_x_xz_z_yy[i] * a_exp + 4.0 * g_x_xz_zzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_zz_yz[i] = -4.0 * g_x_xz_z_yz[i] * a_exp + 4.0 * g_x_xz_zzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_zz_zz[i] = -4.0 * g_x_xz_z_zz[i] * a_exp + 4.0 * g_x_xz_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (540-546)

    #pragma omp simd aligned(g_x_0_z_0_0_yy_xx_xx, g_x_0_z_0_0_yy_xx_xy, g_x_0_z_0_0_yy_xx_xz, g_x_0_z_0_0_yy_xx_yy, g_x_0_z_0_0_yy_xx_yz, g_x_0_z_0_0_yy_xx_zz, g_x_yy_xxz_xx, g_x_yy_xxz_xy, g_x_yy_xxz_xz, g_x_yy_xxz_yy, g_x_yy_xxz_yz, g_x_yy_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_yy_xx_xx[i] = 4.0 * g_x_yy_xxz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_xx_xy[i] = 4.0 * g_x_yy_xxz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_xx_xz[i] = 4.0 * g_x_yy_xxz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_xx_yy[i] = 4.0 * g_x_yy_xxz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_xx_yz[i] = 4.0 * g_x_yy_xxz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_xx_zz[i] = 4.0 * g_x_yy_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (546-552)

    #pragma omp simd aligned(g_x_0_z_0_0_yy_xy_xx, g_x_0_z_0_0_yy_xy_xy, g_x_0_z_0_0_yy_xy_xz, g_x_0_z_0_0_yy_xy_yy, g_x_0_z_0_0_yy_xy_yz, g_x_0_z_0_0_yy_xy_zz, g_x_yy_xyz_xx, g_x_yy_xyz_xy, g_x_yy_xyz_xz, g_x_yy_xyz_yy, g_x_yy_xyz_yz, g_x_yy_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_yy_xy_xx[i] = 4.0 * g_x_yy_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_xy_xy[i] = 4.0 * g_x_yy_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_xy_xz[i] = 4.0 * g_x_yy_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_xy_yy[i] = 4.0 * g_x_yy_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_xy_yz[i] = 4.0 * g_x_yy_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_xy_zz[i] = 4.0 * g_x_yy_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (552-558)

    #pragma omp simd aligned(g_x_0_z_0_0_yy_xz_xx, g_x_0_z_0_0_yy_xz_xy, g_x_0_z_0_0_yy_xz_xz, g_x_0_z_0_0_yy_xz_yy, g_x_0_z_0_0_yy_xz_yz, g_x_0_z_0_0_yy_xz_zz, g_x_yy_x_xx, g_x_yy_x_xy, g_x_yy_x_xz, g_x_yy_x_yy, g_x_yy_x_yz, g_x_yy_x_zz, g_x_yy_xzz_xx, g_x_yy_xzz_xy, g_x_yy_xzz_xz, g_x_yy_xzz_yy, g_x_yy_xzz_yz, g_x_yy_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_yy_xz_xx[i] = -2.0 * g_x_yy_x_xx[i] * a_exp + 4.0 * g_x_yy_xzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_xz_xy[i] = -2.0 * g_x_yy_x_xy[i] * a_exp + 4.0 * g_x_yy_xzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_xz_xz[i] = -2.0 * g_x_yy_x_xz[i] * a_exp + 4.0 * g_x_yy_xzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_xz_yy[i] = -2.0 * g_x_yy_x_yy[i] * a_exp + 4.0 * g_x_yy_xzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_xz_yz[i] = -2.0 * g_x_yy_x_yz[i] * a_exp + 4.0 * g_x_yy_xzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_xz_zz[i] = -2.0 * g_x_yy_x_zz[i] * a_exp + 4.0 * g_x_yy_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (558-564)

    #pragma omp simd aligned(g_x_0_z_0_0_yy_yy_xx, g_x_0_z_0_0_yy_yy_xy, g_x_0_z_0_0_yy_yy_xz, g_x_0_z_0_0_yy_yy_yy, g_x_0_z_0_0_yy_yy_yz, g_x_0_z_0_0_yy_yy_zz, g_x_yy_yyz_xx, g_x_yy_yyz_xy, g_x_yy_yyz_xz, g_x_yy_yyz_yy, g_x_yy_yyz_yz, g_x_yy_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_yy_yy_xx[i] = 4.0 * g_x_yy_yyz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_yy_xy[i] = 4.0 * g_x_yy_yyz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_yy_xz[i] = 4.0 * g_x_yy_yyz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_yy_yy[i] = 4.0 * g_x_yy_yyz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_yy_yz[i] = 4.0 * g_x_yy_yyz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_yy_zz[i] = 4.0 * g_x_yy_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (564-570)

    #pragma omp simd aligned(g_x_0_z_0_0_yy_yz_xx, g_x_0_z_0_0_yy_yz_xy, g_x_0_z_0_0_yy_yz_xz, g_x_0_z_0_0_yy_yz_yy, g_x_0_z_0_0_yy_yz_yz, g_x_0_z_0_0_yy_yz_zz, g_x_yy_y_xx, g_x_yy_y_xy, g_x_yy_y_xz, g_x_yy_y_yy, g_x_yy_y_yz, g_x_yy_y_zz, g_x_yy_yzz_xx, g_x_yy_yzz_xy, g_x_yy_yzz_xz, g_x_yy_yzz_yy, g_x_yy_yzz_yz, g_x_yy_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_yy_yz_xx[i] = -2.0 * g_x_yy_y_xx[i] * a_exp + 4.0 * g_x_yy_yzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_yz_xy[i] = -2.0 * g_x_yy_y_xy[i] * a_exp + 4.0 * g_x_yy_yzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_yz_xz[i] = -2.0 * g_x_yy_y_xz[i] * a_exp + 4.0 * g_x_yy_yzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_yz_yy[i] = -2.0 * g_x_yy_y_yy[i] * a_exp + 4.0 * g_x_yy_yzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_yz_yz[i] = -2.0 * g_x_yy_y_yz[i] * a_exp + 4.0 * g_x_yy_yzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_yz_zz[i] = -2.0 * g_x_yy_y_zz[i] * a_exp + 4.0 * g_x_yy_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (570-576)

    #pragma omp simd aligned(g_x_0_z_0_0_yy_zz_xx, g_x_0_z_0_0_yy_zz_xy, g_x_0_z_0_0_yy_zz_xz, g_x_0_z_0_0_yy_zz_yy, g_x_0_z_0_0_yy_zz_yz, g_x_0_z_0_0_yy_zz_zz, g_x_yy_z_xx, g_x_yy_z_xy, g_x_yy_z_xz, g_x_yy_z_yy, g_x_yy_z_yz, g_x_yy_z_zz, g_x_yy_zzz_xx, g_x_yy_zzz_xy, g_x_yy_zzz_xz, g_x_yy_zzz_yy, g_x_yy_zzz_yz, g_x_yy_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_yy_zz_xx[i] = -4.0 * g_x_yy_z_xx[i] * a_exp + 4.0 * g_x_yy_zzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_zz_xy[i] = -4.0 * g_x_yy_z_xy[i] * a_exp + 4.0 * g_x_yy_zzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_zz_xz[i] = -4.0 * g_x_yy_z_xz[i] * a_exp + 4.0 * g_x_yy_zzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_zz_yy[i] = -4.0 * g_x_yy_z_yy[i] * a_exp + 4.0 * g_x_yy_zzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_zz_yz[i] = -4.0 * g_x_yy_z_yz[i] * a_exp + 4.0 * g_x_yy_zzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_zz_zz[i] = -4.0 * g_x_yy_z_zz[i] * a_exp + 4.0 * g_x_yy_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (576-582)

    #pragma omp simd aligned(g_x_0_z_0_0_yz_xx_xx, g_x_0_z_0_0_yz_xx_xy, g_x_0_z_0_0_yz_xx_xz, g_x_0_z_0_0_yz_xx_yy, g_x_0_z_0_0_yz_xx_yz, g_x_0_z_0_0_yz_xx_zz, g_x_yz_xxz_xx, g_x_yz_xxz_xy, g_x_yz_xxz_xz, g_x_yz_xxz_yy, g_x_yz_xxz_yz, g_x_yz_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_yz_xx_xx[i] = 4.0 * g_x_yz_xxz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_xx_xy[i] = 4.0 * g_x_yz_xxz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_xx_xz[i] = 4.0 * g_x_yz_xxz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_xx_yy[i] = 4.0 * g_x_yz_xxz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_xx_yz[i] = 4.0 * g_x_yz_xxz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_xx_zz[i] = 4.0 * g_x_yz_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (582-588)

    #pragma omp simd aligned(g_x_0_z_0_0_yz_xy_xx, g_x_0_z_0_0_yz_xy_xy, g_x_0_z_0_0_yz_xy_xz, g_x_0_z_0_0_yz_xy_yy, g_x_0_z_0_0_yz_xy_yz, g_x_0_z_0_0_yz_xy_zz, g_x_yz_xyz_xx, g_x_yz_xyz_xy, g_x_yz_xyz_xz, g_x_yz_xyz_yy, g_x_yz_xyz_yz, g_x_yz_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_yz_xy_xx[i] = 4.0 * g_x_yz_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_xy_xy[i] = 4.0 * g_x_yz_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_xy_xz[i] = 4.0 * g_x_yz_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_xy_yy[i] = 4.0 * g_x_yz_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_xy_yz[i] = 4.0 * g_x_yz_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_xy_zz[i] = 4.0 * g_x_yz_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (588-594)

    #pragma omp simd aligned(g_x_0_z_0_0_yz_xz_xx, g_x_0_z_0_0_yz_xz_xy, g_x_0_z_0_0_yz_xz_xz, g_x_0_z_0_0_yz_xz_yy, g_x_0_z_0_0_yz_xz_yz, g_x_0_z_0_0_yz_xz_zz, g_x_yz_x_xx, g_x_yz_x_xy, g_x_yz_x_xz, g_x_yz_x_yy, g_x_yz_x_yz, g_x_yz_x_zz, g_x_yz_xzz_xx, g_x_yz_xzz_xy, g_x_yz_xzz_xz, g_x_yz_xzz_yy, g_x_yz_xzz_yz, g_x_yz_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_yz_xz_xx[i] = -2.0 * g_x_yz_x_xx[i] * a_exp + 4.0 * g_x_yz_xzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_xz_xy[i] = -2.0 * g_x_yz_x_xy[i] * a_exp + 4.0 * g_x_yz_xzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_xz_xz[i] = -2.0 * g_x_yz_x_xz[i] * a_exp + 4.0 * g_x_yz_xzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_xz_yy[i] = -2.0 * g_x_yz_x_yy[i] * a_exp + 4.0 * g_x_yz_xzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_xz_yz[i] = -2.0 * g_x_yz_x_yz[i] * a_exp + 4.0 * g_x_yz_xzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_xz_zz[i] = -2.0 * g_x_yz_x_zz[i] * a_exp + 4.0 * g_x_yz_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (594-600)

    #pragma omp simd aligned(g_x_0_z_0_0_yz_yy_xx, g_x_0_z_0_0_yz_yy_xy, g_x_0_z_0_0_yz_yy_xz, g_x_0_z_0_0_yz_yy_yy, g_x_0_z_0_0_yz_yy_yz, g_x_0_z_0_0_yz_yy_zz, g_x_yz_yyz_xx, g_x_yz_yyz_xy, g_x_yz_yyz_xz, g_x_yz_yyz_yy, g_x_yz_yyz_yz, g_x_yz_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_yz_yy_xx[i] = 4.0 * g_x_yz_yyz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_yy_xy[i] = 4.0 * g_x_yz_yyz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_yy_xz[i] = 4.0 * g_x_yz_yyz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_yy_yy[i] = 4.0 * g_x_yz_yyz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_yy_yz[i] = 4.0 * g_x_yz_yyz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_yy_zz[i] = 4.0 * g_x_yz_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (600-606)

    #pragma omp simd aligned(g_x_0_z_0_0_yz_yz_xx, g_x_0_z_0_0_yz_yz_xy, g_x_0_z_0_0_yz_yz_xz, g_x_0_z_0_0_yz_yz_yy, g_x_0_z_0_0_yz_yz_yz, g_x_0_z_0_0_yz_yz_zz, g_x_yz_y_xx, g_x_yz_y_xy, g_x_yz_y_xz, g_x_yz_y_yy, g_x_yz_y_yz, g_x_yz_y_zz, g_x_yz_yzz_xx, g_x_yz_yzz_xy, g_x_yz_yzz_xz, g_x_yz_yzz_yy, g_x_yz_yzz_yz, g_x_yz_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_yz_yz_xx[i] = -2.0 * g_x_yz_y_xx[i] * a_exp + 4.0 * g_x_yz_yzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_yz_xy[i] = -2.0 * g_x_yz_y_xy[i] * a_exp + 4.0 * g_x_yz_yzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_yz_xz[i] = -2.0 * g_x_yz_y_xz[i] * a_exp + 4.0 * g_x_yz_yzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_yz_yy[i] = -2.0 * g_x_yz_y_yy[i] * a_exp + 4.0 * g_x_yz_yzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_yz_yz[i] = -2.0 * g_x_yz_y_yz[i] * a_exp + 4.0 * g_x_yz_yzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_yz_zz[i] = -2.0 * g_x_yz_y_zz[i] * a_exp + 4.0 * g_x_yz_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (606-612)

    #pragma omp simd aligned(g_x_0_z_0_0_yz_zz_xx, g_x_0_z_0_0_yz_zz_xy, g_x_0_z_0_0_yz_zz_xz, g_x_0_z_0_0_yz_zz_yy, g_x_0_z_0_0_yz_zz_yz, g_x_0_z_0_0_yz_zz_zz, g_x_yz_z_xx, g_x_yz_z_xy, g_x_yz_z_xz, g_x_yz_z_yy, g_x_yz_z_yz, g_x_yz_z_zz, g_x_yz_zzz_xx, g_x_yz_zzz_xy, g_x_yz_zzz_xz, g_x_yz_zzz_yy, g_x_yz_zzz_yz, g_x_yz_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_yz_zz_xx[i] = -4.0 * g_x_yz_z_xx[i] * a_exp + 4.0 * g_x_yz_zzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_zz_xy[i] = -4.0 * g_x_yz_z_xy[i] * a_exp + 4.0 * g_x_yz_zzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_zz_xz[i] = -4.0 * g_x_yz_z_xz[i] * a_exp + 4.0 * g_x_yz_zzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_zz_yy[i] = -4.0 * g_x_yz_z_yy[i] * a_exp + 4.0 * g_x_yz_zzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_zz_yz[i] = -4.0 * g_x_yz_z_yz[i] * a_exp + 4.0 * g_x_yz_zzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_zz_zz[i] = -4.0 * g_x_yz_z_zz[i] * a_exp + 4.0 * g_x_yz_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (612-618)

    #pragma omp simd aligned(g_x_0_z_0_0_zz_xx_xx, g_x_0_z_0_0_zz_xx_xy, g_x_0_z_0_0_zz_xx_xz, g_x_0_z_0_0_zz_xx_yy, g_x_0_z_0_0_zz_xx_yz, g_x_0_z_0_0_zz_xx_zz, g_x_zz_xxz_xx, g_x_zz_xxz_xy, g_x_zz_xxz_xz, g_x_zz_xxz_yy, g_x_zz_xxz_yz, g_x_zz_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_zz_xx_xx[i] = 4.0 * g_x_zz_xxz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_xx_xy[i] = 4.0 * g_x_zz_xxz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_xx_xz[i] = 4.0 * g_x_zz_xxz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_xx_yy[i] = 4.0 * g_x_zz_xxz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_xx_yz[i] = 4.0 * g_x_zz_xxz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_xx_zz[i] = 4.0 * g_x_zz_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (618-624)

    #pragma omp simd aligned(g_x_0_z_0_0_zz_xy_xx, g_x_0_z_0_0_zz_xy_xy, g_x_0_z_0_0_zz_xy_xz, g_x_0_z_0_0_zz_xy_yy, g_x_0_z_0_0_zz_xy_yz, g_x_0_z_0_0_zz_xy_zz, g_x_zz_xyz_xx, g_x_zz_xyz_xy, g_x_zz_xyz_xz, g_x_zz_xyz_yy, g_x_zz_xyz_yz, g_x_zz_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_zz_xy_xx[i] = 4.0 * g_x_zz_xyz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_xy_xy[i] = 4.0 * g_x_zz_xyz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_xy_xz[i] = 4.0 * g_x_zz_xyz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_xy_yy[i] = 4.0 * g_x_zz_xyz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_xy_yz[i] = 4.0 * g_x_zz_xyz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_xy_zz[i] = 4.0 * g_x_zz_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (624-630)

    #pragma omp simd aligned(g_x_0_z_0_0_zz_xz_xx, g_x_0_z_0_0_zz_xz_xy, g_x_0_z_0_0_zz_xz_xz, g_x_0_z_0_0_zz_xz_yy, g_x_0_z_0_0_zz_xz_yz, g_x_0_z_0_0_zz_xz_zz, g_x_zz_x_xx, g_x_zz_x_xy, g_x_zz_x_xz, g_x_zz_x_yy, g_x_zz_x_yz, g_x_zz_x_zz, g_x_zz_xzz_xx, g_x_zz_xzz_xy, g_x_zz_xzz_xz, g_x_zz_xzz_yy, g_x_zz_xzz_yz, g_x_zz_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_zz_xz_xx[i] = -2.0 * g_x_zz_x_xx[i] * a_exp + 4.0 * g_x_zz_xzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_xz_xy[i] = -2.0 * g_x_zz_x_xy[i] * a_exp + 4.0 * g_x_zz_xzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_xz_xz[i] = -2.0 * g_x_zz_x_xz[i] * a_exp + 4.0 * g_x_zz_xzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_xz_yy[i] = -2.0 * g_x_zz_x_yy[i] * a_exp + 4.0 * g_x_zz_xzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_xz_yz[i] = -2.0 * g_x_zz_x_yz[i] * a_exp + 4.0 * g_x_zz_xzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_xz_zz[i] = -2.0 * g_x_zz_x_zz[i] * a_exp + 4.0 * g_x_zz_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (630-636)

    #pragma omp simd aligned(g_x_0_z_0_0_zz_yy_xx, g_x_0_z_0_0_zz_yy_xy, g_x_0_z_0_0_zz_yy_xz, g_x_0_z_0_0_zz_yy_yy, g_x_0_z_0_0_zz_yy_yz, g_x_0_z_0_0_zz_yy_zz, g_x_zz_yyz_xx, g_x_zz_yyz_xy, g_x_zz_yyz_xz, g_x_zz_yyz_yy, g_x_zz_yyz_yz, g_x_zz_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_zz_yy_xx[i] = 4.0 * g_x_zz_yyz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_yy_xy[i] = 4.0 * g_x_zz_yyz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_yy_xz[i] = 4.0 * g_x_zz_yyz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_yy_yy[i] = 4.0 * g_x_zz_yyz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_yy_yz[i] = 4.0 * g_x_zz_yyz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_yy_zz[i] = 4.0 * g_x_zz_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (636-642)

    #pragma omp simd aligned(g_x_0_z_0_0_zz_yz_xx, g_x_0_z_0_0_zz_yz_xy, g_x_0_z_0_0_zz_yz_xz, g_x_0_z_0_0_zz_yz_yy, g_x_0_z_0_0_zz_yz_yz, g_x_0_z_0_0_zz_yz_zz, g_x_zz_y_xx, g_x_zz_y_xy, g_x_zz_y_xz, g_x_zz_y_yy, g_x_zz_y_yz, g_x_zz_y_zz, g_x_zz_yzz_xx, g_x_zz_yzz_xy, g_x_zz_yzz_xz, g_x_zz_yzz_yy, g_x_zz_yzz_yz, g_x_zz_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_zz_yz_xx[i] = -2.0 * g_x_zz_y_xx[i] * a_exp + 4.0 * g_x_zz_yzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_yz_xy[i] = -2.0 * g_x_zz_y_xy[i] * a_exp + 4.0 * g_x_zz_yzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_yz_xz[i] = -2.0 * g_x_zz_y_xz[i] * a_exp + 4.0 * g_x_zz_yzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_yz_yy[i] = -2.0 * g_x_zz_y_yy[i] * a_exp + 4.0 * g_x_zz_yzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_yz_yz[i] = -2.0 * g_x_zz_y_yz[i] * a_exp + 4.0 * g_x_zz_yzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_yz_zz[i] = -2.0 * g_x_zz_y_zz[i] * a_exp + 4.0 * g_x_zz_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (642-648)

    #pragma omp simd aligned(g_x_0_z_0_0_zz_zz_xx, g_x_0_z_0_0_zz_zz_xy, g_x_0_z_0_0_zz_zz_xz, g_x_0_z_0_0_zz_zz_yy, g_x_0_z_0_0_zz_zz_yz, g_x_0_z_0_0_zz_zz_zz, g_x_zz_z_xx, g_x_zz_z_xy, g_x_zz_z_xz, g_x_zz_z_yy, g_x_zz_z_yz, g_x_zz_z_zz, g_x_zz_zzz_xx, g_x_zz_zzz_xy, g_x_zz_zzz_xz, g_x_zz_zzz_yy, g_x_zz_zzz_yz, g_x_zz_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_zz_zz_xx[i] = -4.0 * g_x_zz_z_xx[i] * a_exp + 4.0 * g_x_zz_zzz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_zz_xy[i] = -4.0 * g_x_zz_z_xy[i] * a_exp + 4.0 * g_x_zz_zzz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_zz_xz[i] = -4.0 * g_x_zz_z_xz[i] * a_exp + 4.0 * g_x_zz_zzz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_zz_yy[i] = -4.0 * g_x_zz_z_yy[i] * a_exp + 4.0 * g_x_zz_zzz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_zz_yz[i] = -4.0 * g_x_zz_z_yz[i] * a_exp + 4.0 * g_x_zz_zzz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_zz_zz[i] = -4.0 * g_x_zz_z_zz[i] * a_exp + 4.0 * g_x_zz_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (648-654)

    #pragma omp simd aligned(g_y_0_x_0_0_xx_xx_xx, g_y_0_x_0_0_xx_xx_xy, g_y_0_x_0_0_xx_xx_xz, g_y_0_x_0_0_xx_xx_yy, g_y_0_x_0_0_xx_xx_yz, g_y_0_x_0_0_xx_xx_zz, g_y_xx_x_xx, g_y_xx_x_xy, g_y_xx_x_xz, g_y_xx_x_yy, g_y_xx_x_yz, g_y_xx_x_zz, g_y_xx_xxx_xx, g_y_xx_xxx_xy, g_y_xx_xxx_xz, g_y_xx_xxx_yy, g_y_xx_xxx_yz, g_y_xx_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xx_xx_xx[i] = -4.0 * g_y_xx_x_xx[i] * a_exp + 4.0 * g_y_xx_xxx_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_xx_xy[i] = -4.0 * g_y_xx_x_xy[i] * a_exp + 4.0 * g_y_xx_xxx_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_xx_xz[i] = -4.0 * g_y_xx_x_xz[i] * a_exp + 4.0 * g_y_xx_xxx_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_xx_yy[i] = -4.0 * g_y_xx_x_yy[i] * a_exp + 4.0 * g_y_xx_xxx_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_xx_yz[i] = -4.0 * g_y_xx_x_yz[i] * a_exp + 4.0 * g_y_xx_xxx_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_xx_zz[i] = -4.0 * g_y_xx_x_zz[i] * a_exp + 4.0 * g_y_xx_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (654-660)

    #pragma omp simd aligned(g_y_0_x_0_0_xx_xy_xx, g_y_0_x_0_0_xx_xy_xy, g_y_0_x_0_0_xx_xy_xz, g_y_0_x_0_0_xx_xy_yy, g_y_0_x_0_0_xx_xy_yz, g_y_0_x_0_0_xx_xy_zz, g_y_xx_xxy_xx, g_y_xx_xxy_xy, g_y_xx_xxy_xz, g_y_xx_xxy_yy, g_y_xx_xxy_yz, g_y_xx_xxy_zz, g_y_xx_y_xx, g_y_xx_y_xy, g_y_xx_y_xz, g_y_xx_y_yy, g_y_xx_y_yz, g_y_xx_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xx_xy_xx[i] = -2.0 * g_y_xx_y_xx[i] * a_exp + 4.0 * g_y_xx_xxy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_xy_xy[i] = -2.0 * g_y_xx_y_xy[i] * a_exp + 4.0 * g_y_xx_xxy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_xy_xz[i] = -2.0 * g_y_xx_y_xz[i] * a_exp + 4.0 * g_y_xx_xxy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_xy_yy[i] = -2.0 * g_y_xx_y_yy[i] * a_exp + 4.0 * g_y_xx_xxy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_xy_yz[i] = -2.0 * g_y_xx_y_yz[i] * a_exp + 4.0 * g_y_xx_xxy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_xy_zz[i] = -2.0 * g_y_xx_y_zz[i] * a_exp + 4.0 * g_y_xx_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (660-666)

    #pragma omp simd aligned(g_y_0_x_0_0_xx_xz_xx, g_y_0_x_0_0_xx_xz_xy, g_y_0_x_0_0_xx_xz_xz, g_y_0_x_0_0_xx_xz_yy, g_y_0_x_0_0_xx_xz_yz, g_y_0_x_0_0_xx_xz_zz, g_y_xx_xxz_xx, g_y_xx_xxz_xy, g_y_xx_xxz_xz, g_y_xx_xxz_yy, g_y_xx_xxz_yz, g_y_xx_xxz_zz, g_y_xx_z_xx, g_y_xx_z_xy, g_y_xx_z_xz, g_y_xx_z_yy, g_y_xx_z_yz, g_y_xx_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xx_xz_xx[i] = -2.0 * g_y_xx_z_xx[i] * a_exp + 4.0 * g_y_xx_xxz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_xz_xy[i] = -2.0 * g_y_xx_z_xy[i] * a_exp + 4.0 * g_y_xx_xxz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_xz_xz[i] = -2.0 * g_y_xx_z_xz[i] * a_exp + 4.0 * g_y_xx_xxz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_xz_yy[i] = -2.0 * g_y_xx_z_yy[i] * a_exp + 4.0 * g_y_xx_xxz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_xz_yz[i] = -2.0 * g_y_xx_z_yz[i] * a_exp + 4.0 * g_y_xx_xxz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_xz_zz[i] = -2.0 * g_y_xx_z_zz[i] * a_exp + 4.0 * g_y_xx_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (666-672)

    #pragma omp simd aligned(g_y_0_x_0_0_xx_yy_xx, g_y_0_x_0_0_xx_yy_xy, g_y_0_x_0_0_xx_yy_xz, g_y_0_x_0_0_xx_yy_yy, g_y_0_x_0_0_xx_yy_yz, g_y_0_x_0_0_xx_yy_zz, g_y_xx_xyy_xx, g_y_xx_xyy_xy, g_y_xx_xyy_xz, g_y_xx_xyy_yy, g_y_xx_xyy_yz, g_y_xx_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xx_yy_xx[i] = 4.0 * g_y_xx_xyy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_yy_xy[i] = 4.0 * g_y_xx_xyy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_yy_xz[i] = 4.0 * g_y_xx_xyy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_yy_yy[i] = 4.0 * g_y_xx_xyy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_yy_yz[i] = 4.0 * g_y_xx_xyy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_yy_zz[i] = 4.0 * g_y_xx_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (672-678)

    #pragma omp simd aligned(g_y_0_x_0_0_xx_yz_xx, g_y_0_x_0_0_xx_yz_xy, g_y_0_x_0_0_xx_yz_xz, g_y_0_x_0_0_xx_yz_yy, g_y_0_x_0_0_xx_yz_yz, g_y_0_x_0_0_xx_yz_zz, g_y_xx_xyz_xx, g_y_xx_xyz_xy, g_y_xx_xyz_xz, g_y_xx_xyz_yy, g_y_xx_xyz_yz, g_y_xx_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xx_yz_xx[i] = 4.0 * g_y_xx_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_yz_xy[i] = 4.0 * g_y_xx_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_yz_xz[i] = 4.0 * g_y_xx_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_yz_yy[i] = 4.0 * g_y_xx_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_yz_yz[i] = 4.0 * g_y_xx_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_yz_zz[i] = 4.0 * g_y_xx_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (678-684)

    #pragma omp simd aligned(g_y_0_x_0_0_xx_zz_xx, g_y_0_x_0_0_xx_zz_xy, g_y_0_x_0_0_xx_zz_xz, g_y_0_x_0_0_xx_zz_yy, g_y_0_x_0_0_xx_zz_yz, g_y_0_x_0_0_xx_zz_zz, g_y_xx_xzz_xx, g_y_xx_xzz_xy, g_y_xx_xzz_xz, g_y_xx_xzz_yy, g_y_xx_xzz_yz, g_y_xx_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xx_zz_xx[i] = 4.0 * g_y_xx_xzz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_zz_xy[i] = 4.0 * g_y_xx_xzz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_zz_xz[i] = 4.0 * g_y_xx_xzz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_zz_yy[i] = 4.0 * g_y_xx_xzz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_zz_yz[i] = 4.0 * g_y_xx_xzz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_zz_zz[i] = 4.0 * g_y_xx_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (684-690)

    #pragma omp simd aligned(g_y_0_x_0_0_xy_xx_xx, g_y_0_x_0_0_xy_xx_xy, g_y_0_x_0_0_xy_xx_xz, g_y_0_x_0_0_xy_xx_yy, g_y_0_x_0_0_xy_xx_yz, g_y_0_x_0_0_xy_xx_zz, g_y_xy_x_xx, g_y_xy_x_xy, g_y_xy_x_xz, g_y_xy_x_yy, g_y_xy_x_yz, g_y_xy_x_zz, g_y_xy_xxx_xx, g_y_xy_xxx_xy, g_y_xy_xxx_xz, g_y_xy_xxx_yy, g_y_xy_xxx_yz, g_y_xy_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xy_xx_xx[i] = -4.0 * g_y_xy_x_xx[i] * a_exp + 4.0 * g_y_xy_xxx_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_xx_xy[i] = -4.0 * g_y_xy_x_xy[i] * a_exp + 4.0 * g_y_xy_xxx_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_xx_xz[i] = -4.0 * g_y_xy_x_xz[i] * a_exp + 4.0 * g_y_xy_xxx_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_xx_yy[i] = -4.0 * g_y_xy_x_yy[i] * a_exp + 4.0 * g_y_xy_xxx_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_xx_yz[i] = -4.0 * g_y_xy_x_yz[i] * a_exp + 4.0 * g_y_xy_xxx_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_xx_zz[i] = -4.0 * g_y_xy_x_zz[i] * a_exp + 4.0 * g_y_xy_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (690-696)

    #pragma omp simd aligned(g_y_0_x_0_0_xy_xy_xx, g_y_0_x_0_0_xy_xy_xy, g_y_0_x_0_0_xy_xy_xz, g_y_0_x_0_0_xy_xy_yy, g_y_0_x_0_0_xy_xy_yz, g_y_0_x_0_0_xy_xy_zz, g_y_xy_xxy_xx, g_y_xy_xxy_xy, g_y_xy_xxy_xz, g_y_xy_xxy_yy, g_y_xy_xxy_yz, g_y_xy_xxy_zz, g_y_xy_y_xx, g_y_xy_y_xy, g_y_xy_y_xz, g_y_xy_y_yy, g_y_xy_y_yz, g_y_xy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xy_xy_xx[i] = -2.0 * g_y_xy_y_xx[i] * a_exp + 4.0 * g_y_xy_xxy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_xy_xy[i] = -2.0 * g_y_xy_y_xy[i] * a_exp + 4.0 * g_y_xy_xxy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_xy_xz[i] = -2.0 * g_y_xy_y_xz[i] * a_exp + 4.0 * g_y_xy_xxy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_xy_yy[i] = -2.0 * g_y_xy_y_yy[i] * a_exp + 4.0 * g_y_xy_xxy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_xy_yz[i] = -2.0 * g_y_xy_y_yz[i] * a_exp + 4.0 * g_y_xy_xxy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_xy_zz[i] = -2.0 * g_y_xy_y_zz[i] * a_exp + 4.0 * g_y_xy_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (696-702)

    #pragma omp simd aligned(g_y_0_x_0_0_xy_xz_xx, g_y_0_x_0_0_xy_xz_xy, g_y_0_x_0_0_xy_xz_xz, g_y_0_x_0_0_xy_xz_yy, g_y_0_x_0_0_xy_xz_yz, g_y_0_x_0_0_xy_xz_zz, g_y_xy_xxz_xx, g_y_xy_xxz_xy, g_y_xy_xxz_xz, g_y_xy_xxz_yy, g_y_xy_xxz_yz, g_y_xy_xxz_zz, g_y_xy_z_xx, g_y_xy_z_xy, g_y_xy_z_xz, g_y_xy_z_yy, g_y_xy_z_yz, g_y_xy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xy_xz_xx[i] = -2.0 * g_y_xy_z_xx[i] * a_exp + 4.0 * g_y_xy_xxz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_xz_xy[i] = -2.0 * g_y_xy_z_xy[i] * a_exp + 4.0 * g_y_xy_xxz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_xz_xz[i] = -2.0 * g_y_xy_z_xz[i] * a_exp + 4.0 * g_y_xy_xxz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_xz_yy[i] = -2.0 * g_y_xy_z_yy[i] * a_exp + 4.0 * g_y_xy_xxz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_xz_yz[i] = -2.0 * g_y_xy_z_yz[i] * a_exp + 4.0 * g_y_xy_xxz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_xz_zz[i] = -2.0 * g_y_xy_z_zz[i] * a_exp + 4.0 * g_y_xy_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (702-708)

    #pragma omp simd aligned(g_y_0_x_0_0_xy_yy_xx, g_y_0_x_0_0_xy_yy_xy, g_y_0_x_0_0_xy_yy_xz, g_y_0_x_0_0_xy_yy_yy, g_y_0_x_0_0_xy_yy_yz, g_y_0_x_0_0_xy_yy_zz, g_y_xy_xyy_xx, g_y_xy_xyy_xy, g_y_xy_xyy_xz, g_y_xy_xyy_yy, g_y_xy_xyy_yz, g_y_xy_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xy_yy_xx[i] = 4.0 * g_y_xy_xyy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_yy_xy[i] = 4.0 * g_y_xy_xyy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_yy_xz[i] = 4.0 * g_y_xy_xyy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_yy_yy[i] = 4.0 * g_y_xy_xyy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_yy_yz[i] = 4.0 * g_y_xy_xyy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_yy_zz[i] = 4.0 * g_y_xy_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (708-714)

    #pragma omp simd aligned(g_y_0_x_0_0_xy_yz_xx, g_y_0_x_0_0_xy_yz_xy, g_y_0_x_0_0_xy_yz_xz, g_y_0_x_0_0_xy_yz_yy, g_y_0_x_0_0_xy_yz_yz, g_y_0_x_0_0_xy_yz_zz, g_y_xy_xyz_xx, g_y_xy_xyz_xy, g_y_xy_xyz_xz, g_y_xy_xyz_yy, g_y_xy_xyz_yz, g_y_xy_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xy_yz_xx[i] = 4.0 * g_y_xy_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_yz_xy[i] = 4.0 * g_y_xy_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_yz_xz[i] = 4.0 * g_y_xy_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_yz_yy[i] = 4.0 * g_y_xy_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_yz_yz[i] = 4.0 * g_y_xy_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_yz_zz[i] = 4.0 * g_y_xy_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (714-720)

    #pragma omp simd aligned(g_y_0_x_0_0_xy_zz_xx, g_y_0_x_0_0_xy_zz_xy, g_y_0_x_0_0_xy_zz_xz, g_y_0_x_0_0_xy_zz_yy, g_y_0_x_0_0_xy_zz_yz, g_y_0_x_0_0_xy_zz_zz, g_y_xy_xzz_xx, g_y_xy_xzz_xy, g_y_xy_xzz_xz, g_y_xy_xzz_yy, g_y_xy_xzz_yz, g_y_xy_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xy_zz_xx[i] = 4.0 * g_y_xy_xzz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_zz_xy[i] = 4.0 * g_y_xy_xzz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_zz_xz[i] = 4.0 * g_y_xy_xzz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_zz_yy[i] = 4.0 * g_y_xy_xzz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_zz_yz[i] = 4.0 * g_y_xy_xzz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_zz_zz[i] = 4.0 * g_y_xy_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (720-726)

    #pragma omp simd aligned(g_y_0_x_0_0_xz_xx_xx, g_y_0_x_0_0_xz_xx_xy, g_y_0_x_0_0_xz_xx_xz, g_y_0_x_0_0_xz_xx_yy, g_y_0_x_0_0_xz_xx_yz, g_y_0_x_0_0_xz_xx_zz, g_y_xz_x_xx, g_y_xz_x_xy, g_y_xz_x_xz, g_y_xz_x_yy, g_y_xz_x_yz, g_y_xz_x_zz, g_y_xz_xxx_xx, g_y_xz_xxx_xy, g_y_xz_xxx_xz, g_y_xz_xxx_yy, g_y_xz_xxx_yz, g_y_xz_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xz_xx_xx[i] = -4.0 * g_y_xz_x_xx[i] * a_exp + 4.0 * g_y_xz_xxx_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_xx_xy[i] = -4.0 * g_y_xz_x_xy[i] * a_exp + 4.0 * g_y_xz_xxx_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_xx_xz[i] = -4.0 * g_y_xz_x_xz[i] * a_exp + 4.0 * g_y_xz_xxx_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_xx_yy[i] = -4.0 * g_y_xz_x_yy[i] * a_exp + 4.0 * g_y_xz_xxx_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_xx_yz[i] = -4.0 * g_y_xz_x_yz[i] * a_exp + 4.0 * g_y_xz_xxx_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_xx_zz[i] = -4.0 * g_y_xz_x_zz[i] * a_exp + 4.0 * g_y_xz_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (726-732)

    #pragma omp simd aligned(g_y_0_x_0_0_xz_xy_xx, g_y_0_x_0_0_xz_xy_xy, g_y_0_x_0_0_xz_xy_xz, g_y_0_x_0_0_xz_xy_yy, g_y_0_x_0_0_xz_xy_yz, g_y_0_x_0_0_xz_xy_zz, g_y_xz_xxy_xx, g_y_xz_xxy_xy, g_y_xz_xxy_xz, g_y_xz_xxy_yy, g_y_xz_xxy_yz, g_y_xz_xxy_zz, g_y_xz_y_xx, g_y_xz_y_xy, g_y_xz_y_xz, g_y_xz_y_yy, g_y_xz_y_yz, g_y_xz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xz_xy_xx[i] = -2.0 * g_y_xz_y_xx[i] * a_exp + 4.0 * g_y_xz_xxy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_xy_xy[i] = -2.0 * g_y_xz_y_xy[i] * a_exp + 4.0 * g_y_xz_xxy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_xy_xz[i] = -2.0 * g_y_xz_y_xz[i] * a_exp + 4.0 * g_y_xz_xxy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_xy_yy[i] = -2.0 * g_y_xz_y_yy[i] * a_exp + 4.0 * g_y_xz_xxy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_xy_yz[i] = -2.0 * g_y_xz_y_yz[i] * a_exp + 4.0 * g_y_xz_xxy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_xy_zz[i] = -2.0 * g_y_xz_y_zz[i] * a_exp + 4.0 * g_y_xz_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (732-738)

    #pragma omp simd aligned(g_y_0_x_0_0_xz_xz_xx, g_y_0_x_0_0_xz_xz_xy, g_y_0_x_0_0_xz_xz_xz, g_y_0_x_0_0_xz_xz_yy, g_y_0_x_0_0_xz_xz_yz, g_y_0_x_0_0_xz_xz_zz, g_y_xz_xxz_xx, g_y_xz_xxz_xy, g_y_xz_xxz_xz, g_y_xz_xxz_yy, g_y_xz_xxz_yz, g_y_xz_xxz_zz, g_y_xz_z_xx, g_y_xz_z_xy, g_y_xz_z_xz, g_y_xz_z_yy, g_y_xz_z_yz, g_y_xz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xz_xz_xx[i] = -2.0 * g_y_xz_z_xx[i] * a_exp + 4.0 * g_y_xz_xxz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_xz_xy[i] = -2.0 * g_y_xz_z_xy[i] * a_exp + 4.0 * g_y_xz_xxz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_xz_xz[i] = -2.0 * g_y_xz_z_xz[i] * a_exp + 4.0 * g_y_xz_xxz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_xz_yy[i] = -2.0 * g_y_xz_z_yy[i] * a_exp + 4.0 * g_y_xz_xxz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_xz_yz[i] = -2.0 * g_y_xz_z_yz[i] * a_exp + 4.0 * g_y_xz_xxz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_xz_zz[i] = -2.0 * g_y_xz_z_zz[i] * a_exp + 4.0 * g_y_xz_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (738-744)

    #pragma omp simd aligned(g_y_0_x_0_0_xz_yy_xx, g_y_0_x_0_0_xz_yy_xy, g_y_0_x_0_0_xz_yy_xz, g_y_0_x_0_0_xz_yy_yy, g_y_0_x_0_0_xz_yy_yz, g_y_0_x_0_0_xz_yy_zz, g_y_xz_xyy_xx, g_y_xz_xyy_xy, g_y_xz_xyy_xz, g_y_xz_xyy_yy, g_y_xz_xyy_yz, g_y_xz_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xz_yy_xx[i] = 4.0 * g_y_xz_xyy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_yy_xy[i] = 4.0 * g_y_xz_xyy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_yy_xz[i] = 4.0 * g_y_xz_xyy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_yy_yy[i] = 4.0 * g_y_xz_xyy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_yy_yz[i] = 4.0 * g_y_xz_xyy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_yy_zz[i] = 4.0 * g_y_xz_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (744-750)

    #pragma omp simd aligned(g_y_0_x_0_0_xz_yz_xx, g_y_0_x_0_0_xz_yz_xy, g_y_0_x_0_0_xz_yz_xz, g_y_0_x_0_0_xz_yz_yy, g_y_0_x_0_0_xz_yz_yz, g_y_0_x_0_0_xz_yz_zz, g_y_xz_xyz_xx, g_y_xz_xyz_xy, g_y_xz_xyz_xz, g_y_xz_xyz_yy, g_y_xz_xyz_yz, g_y_xz_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xz_yz_xx[i] = 4.0 * g_y_xz_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_yz_xy[i] = 4.0 * g_y_xz_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_yz_xz[i] = 4.0 * g_y_xz_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_yz_yy[i] = 4.0 * g_y_xz_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_yz_yz[i] = 4.0 * g_y_xz_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_yz_zz[i] = 4.0 * g_y_xz_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (750-756)

    #pragma omp simd aligned(g_y_0_x_0_0_xz_zz_xx, g_y_0_x_0_0_xz_zz_xy, g_y_0_x_0_0_xz_zz_xz, g_y_0_x_0_0_xz_zz_yy, g_y_0_x_0_0_xz_zz_yz, g_y_0_x_0_0_xz_zz_zz, g_y_xz_xzz_xx, g_y_xz_xzz_xy, g_y_xz_xzz_xz, g_y_xz_xzz_yy, g_y_xz_xzz_yz, g_y_xz_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xz_zz_xx[i] = 4.0 * g_y_xz_xzz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_zz_xy[i] = 4.0 * g_y_xz_xzz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_zz_xz[i] = 4.0 * g_y_xz_xzz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_zz_yy[i] = 4.0 * g_y_xz_xzz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_zz_yz[i] = 4.0 * g_y_xz_xzz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_zz_zz[i] = 4.0 * g_y_xz_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (756-762)

    #pragma omp simd aligned(g_y_0_x_0_0_yy_xx_xx, g_y_0_x_0_0_yy_xx_xy, g_y_0_x_0_0_yy_xx_xz, g_y_0_x_0_0_yy_xx_yy, g_y_0_x_0_0_yy_xx_yz, g_y_0_x_0_0_yy_xx_zz, g_y_yy_x_xx, g_y_yy_x_xy, g_y_yy_x_xz, g_y_yy_x_yy, g_y_yy_x_yz, g_y_yy_x_zz, g_y_yy_xxx_xx, g_y_yy_xxx_xy, g_y_yy_xxx_xz, g_y_yy_xxx_yy, g_y_yy_xxx_yz, g_y_yy_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_yy_xx_xx[i] = -4.0 * g_y_yy_x_xx[i] * a_exp + 4.0 * g_y_yy_xxx_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_xx_xy[i] = -4.0 * g_y_yy_x_xy[i] * a_exp + 4.0 * g_y_yy_xxx_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_xx_xz[i] = -4.0 * g_y_yy_x_xz[i] * a_exp + 4.0 * g_y_yy_xxx_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_xx_yy[i] = -4.0 * g_y_yy_x_yy[i] * a_exp + 4.0 * g_y_yy_xxx_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_xx_yz[i] = -4.0 * g_y_yy_x_yz[i] * a_exp + 4.0 * g_y_yy_xxx_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_xx_zz[i] = -4.0 * g_y_yy_x_zz[i] * a_exp + 4.0 * g_y_yy_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (762-768)

    #pragma omp simd aligned(g_y_0_x_0_0_yy_xy_xx, g_y_0_x_0_0_yy_xy_xy, g_y_0_x_0_0_yy_xy_xz, g_y_0_x_0_0_yy_xy_yy, g_y_0_x_0_0_yy_xy_yz, g_y_0_x_0_0_yy_xy_zz, g_y_yy_xxy_xx, g_y_yy_xxy_xy, g_y_yy_xxy_xz, g_y_yy_xxy_yy, g_y_yy_xxy_yz, g_y_yy_xxy_zz, g_y_yy_y_xx, g_y_yy_y_xy, g_y_yy_y_xz, g_y_yy_y_yy, g_y_yy_y_yz, g_y_yy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_yy_xy_xx[i] = -2.0 * g_y_yy_y_xx[i] * a_exp + 4.0 * g_y_yy_xxy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_xy_xy[i] = -2.0 * g_y_yy_y_xy[i] * a_exp + 4.0 * g_y_yy_xxy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_xy_xz[i] = -2.0 * g_y_yy_y_xz[i] * a_exp + 4.0 * g_y_yy_xxy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_xy_yy[i] = -2.0 * g_y_yy_y_yy[i] * a_exp + 4.0 * g_y_yy_xxy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_xy_yz[i] = -2.0 * g_y_yy_y_yz[i] * a_exp + 4.0 * g_y_yy_xxy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_xy_zz[i] = -2.0 * g_y_yy_y_zz[i] * a_exp + 4.0 * g_y_yy_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (768-774)

    #pragma omp simd aligned(g_y_0_x_0_0_yy_xz_xx, g_y_0_x_0_0_yy_xz_xy, g_y_0_x_0_0_yy_xz_xz, g_y_0_x_0_0_yy_xz_yy, g_y_0_x_0_0_yy_xz_yz, g_y_0_x_0_0_yy_xz_zz, g_y_yy_xxz_xx, g_y_yy_xxz_xy, g_y_yy_xxz_xz, g_y_yy_xxz_yy, g_y_yy_xxz_yz, g_y_yy_xxz_zz, g_y_yy_z_xx, g_y_yy_z_xy, g_y_yy_z_xz, g_y_yy_z_yy, g_y_yy_z_yz, g_y_yy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_yy_xz_xx[i] = -2.0 * g_y_yy_z_xx[i] * a_exp + 4.0 * g_y_yy_xxz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_xz_xy[i] = -2.0 * g_y_yy_z_xy[i] * a_exp + 4.0 * g_y_yy_xxz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_xz_xz[i] = -2.0 * g_y_yy_z_xz[i] * a_exp + 4.0 * g_y_yy_xxz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_xz_yy[i] = -2.0 * g_y_yy_z_yy[i] * a_exp + 4.0 * g_y_yy_xxz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_xz_yz[i] = -2.0 * g_y_yy_z_yz[i] * a_exp + 4.0 * g_y_yy_xxz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_xz_zz[i] = -2.0 * g_y_yy_z_zz[i] * a_exp + 4.0 * g_y_yy_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (774-780)

    #pragma omp simd aligned(g_y_0_x_0_0_yy_yy_xx, g_y_0_x_0_0_yy_yy_xy, g_y_0_x_0_0_yy_yy_xz, g_y_0_x_0_0_yy_yy_yy, g_y_0_x_0_0_yy_yy_yz, g_y_0_x_0_0_yy_yy_zz, g_y_yy_xyy_xx, g_y_yy_xyy_xy, g_y_yy_xyy_xz, g_y_yy_xyy_yy, g_y_yy_xyy_yz, g_y_yy_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_yy_yy_xx[i] = 4.0 * g_y_yy_xyy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_yy_xy[i] = 4.0 * g_y_yy_xyy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_yy_xz[i] = 4.0 * g_y_yy_xyy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_yy_yy[i] = 4.0 * g_y_yy_xyy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_yy_yz[i] = 4.0 * g_y_yy_xyy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_yy_zz[i] = 4.0 * g_y_yy_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (780-786)

    #pragma omp simd aligned(g_y_0_x_0_0_yy_yz_xx, g_y_0_x_0_0_yy_yz_xy, g_y_0_x_0_0_yy_yz_xz, g_y_0_x_0_0_yy_yz_yy, g_y_0_x_0_0_yy_yz_yz, g_y_0_x_0_0_yy_yz_zz, g_y_yy_xyz_xx, g_y_yy_xyz_xy, g_y_yy_xyz_xz, g_y_yy_xyz_yy, g_y_yy_xyz_yz, g_y_yy_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_yy_yz_xx[i] = 4.0 * g_y_yy_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_yz_xy[i] = 4.0 * g_y_yy_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_yz_xz[i] = 4.0 * g_y_yy_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_yz_yy[i] = 4.0 * g_y_yy_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_yz_yz[i] = 4.0 * g_y_yy_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_yz_zz[i] = 4.0 * g_y_yy_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (786-792)

    #pragma omp simd aligned(g_y_0_x_0_0_yy_zz_xx, g_y_0_x_0_0_yy_zz_xy, g_y_0_x_0_0_yy_zz_xz, g_y_0_x_0_0_yy_zz_yy, g_y_0_x_0_0_yy_zz_yz, g_y_0_x_0_0_yy_zz_zz, g_y_yy_xzz_xx, g_y_yy_xzz_xy, g_y_yy_xzz_xz, g_y_yy_xzz_yy, g_y_yy_xzz_yz, g_y_yy_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_yy_zz_xx[i] = 4.0 * g_y_yy_xzz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_zz_xy[i] = 4.0 * g_y_yy_xzz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_zz_xz[i] = 4.0 * g_y_yy_xzz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_zz_yy[i] = 4.0 * g_y_yy_xzz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_zz_yz[i] = 4.0 * g_y_yy_xzz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_zz_zz[i] = 4.0 * g_y_yy_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (792-798)

    #pragma omp simd aligned(g_y_0_x_0_0_yz_xx_xx, g_y_0_x_0_0_yz_xx_xy, g_y_0_x_0_0_yz_xx_xz, g_y_0_x_0_0_yz_xx_yy, g_y_0_x_0_0_yz_xx_yz, g_y_0_x_0_0_yz_xx_zz, g_y_yz_x_xx, g_y_yz_x_xy, g_y_yz_x_xz, g_y_yz_x_yy, g_y_yz_x_yz, g_y_yz_x_zz, g_y_yz_xxx_xx, g_y_yz_xxx_xy, g_y_yz_xxx_xz, g_y_yz_xxx_yy, g_y_yz_xxx_yz, g_y_yz_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_yz_xx_xx[i] = -4.0 * g_y_yz_x_xx[i] * a_exp + 4.0 * g_y_yz_xxx_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_xx_xy[i] = -4.0 * g_y_yz_x_xy[i] * a_exp + 4.0 * g_y_yz_xxx_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_xx_xz[i] = -4.0 * g_y_yz_x_xz[i] * a_exp + 4.0 * g_y_yz_xxx_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_xx_yy[i] = -4.0 * g_y_yz_x_yy[i] * a_exp + 4.0 * g_y_yz_xxx_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_xx_yz[i] = -4.0 * g_y_yz_x_yz[i] * a_exp + 4.0 * g_y_yz_xxx_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_xx_zz[i] = -4.0 * g_y_yz_x_zz[i] * a_exp + 4.0 * g_y_yz_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (798-804)

    #pragma omp simd aligned(g_y_0_x_0_0_yz_xy_xx, g_y_0_x_0_0_yz_xy_xy, g_y_0_x_0_0_yz_xy_xz, g_y_0_x_0_0_yz_xy_yy, g_y_0_x_0_0_yz_xy_yz, g_y_0_x_0_0_yz_xy_zz, g_y_yz_xxy_xx, g_y_yz_xxy_xy, g_y_yz_xxy_xz, g_y_yz_xxy_yy, g_y_yz_xxy_yz, g_y_yz_xxy_zz, g_y_yz_y_xx, g_y_yz_y_xy, g_y_yz_y_xz, g_y_yz_y_yy, g_y_yz_y_yz, g_y_yz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_yz_xy_xx[i] = -2.0 * g_y_yz_y_xx[i] * a_exp + 4.0 * g_y_yz_xxy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_xy_xy[i] = -2.0 * g_y_yz_y_xy[i] * a_exp + 4.0 * g_y_yz_xxy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_xy_xz[i] = -2.0 * g_y_yz_y_xz[i] * a_exp + 4.0 * g_y_yz_xxy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_xy_yy[i] = -2.0 * g_y_yz_y_yy[i] * a_exp + 4.0 * g_y_yz_xxy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_xy_yz[i] = -2.0 * g_y_yz_y_yz[i] * a_exp + 4.0 * g_y_yz_xxy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_xy_zz[i] = -2.0 * g_y_yz_y_zz[i] * a_exp + 4.0 * g_y_yz_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (804-810)

    #pragma omp simd aligned(g_y_0_x_0_0_yz_xz_xx, g_y_0_x_0_0_yz_xz_xy, g_y_0_x_0_0_yz_xz_xz, g_y_0_x_0_0_yz_xz_yy, g_y_0_x_0_0_yz_xz_yz, g_y_0_x_0_0_yz_xz_zz, g_y_yz_xxz_xx, g_y_yz_xxz_xy, g_y_yz_xxz_xz, g_y_yz_xxz_yy, g_y_yz_xxz_yz, g_y_yz_xxz_zz, g_y_yz_z_xx, g_y_yz_z_xy, g_y_yz_z_xz, g_y_yz_z_yy, g_y_yz_z_yz, g_y_yz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_yz_xz_xx[i] = -2.0 * g_y_yz_z_xx[i] * a_exp + 4.0 * g_y_yz_xxz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_xz_xy[i] = -2.0 * g_y_yz_z_xy[i] * a_exp + 4.0 * g_y_yz_xxz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_xz_xz[i] = -2.0 * g_y_yz_z_xz[i] * a_exp + 4.0 * g_y_yz_xxz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_xz_yy[i] = -2.0 * g_y_yz_z_yy[i] * a_exp + 4.0 * g_y_yz_xxz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_xz_yz[i] = -2.0 * g_y_yz_z_yz[i] * a_exp + 4.0 * g_y_yz_xxz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_xz_zz[i] = -2.0 * g_y_yz_z_zz[i] * a_exp + 4.0 * g_y_yz_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (810-816)

    #pragma omp simd aligned(g_y_0_x_0_0_yz_yy_xx, g_y_0_x_0_0_yz_yy_xy, g_y_0_x_0_0_yz_yy_xz, g_y_0_x_0_0_yz_yy_yy, g_y_0_x_0_0_yz_yy_yz, g_y_0_x_0_0_yz_yy_zz, g_y_yz_xyy_xx, g_y_yz_xyy_xy, g_y_yz_xyy_xz, g_y_yz_xyy_yy, g_y_yz_xyy_yz, g_y_yz_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_yz_yy_xx[i] = 4.0 * g_y_yz_xyy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_yy_xy[i] = 4.0 * g_y_yz_xyy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_yy_xz[i] = 4.0 * g_y_yz_xyy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_yy_yy[i] = 4.0 * g_y_yz_xyy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_yy_yz[i] = 4.0 * g_y_yz_xyy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_yy_zz[i] = 4.0 * g_y_yz_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (816-822)

    #pragma omp simd aligned(g_y_0_x_0_0_yz_yz_xx, g_y_0_x_0_0_yz_yz_xy, g_y_0_x_0_0_yz_yz_xz, g_y_0_x_0_0_yz_yz_yy, g_y_0_x_0_0_yz_yz_yz, g_y_0_x_0_0_yz_yz_zz, g_y_yz_xyz_xx, g_y_yz_xyz_xy, g_y_yz_xyz_xz, g_y_yz_xyz_yy, g_y_yz_xyz_yz, g_y_yz_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_yz_yz_xx[i] = 4.0 * g_y_yz_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_yz_xy[i] = 4.0 * g_y_yz_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_yz_xz[i] = 4.0 * g_y_yz_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_yz_yy[i] = 4.0 * g_y_yz_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_yz_yz[i] = 4.0 * g_y_yz_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_yz_zz[i] = 4.0 * g_y_yz_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (822-828)

    #pragma omp simd aligned(g_y_0_x_0_0_yz_zz_xx, g_y_0_x_0_0_yz_zz_xy, g_y_0_x_0_0_yz_zz_xz, g_y_0_x_0_0_yz_zz_yy, g_y_0_x_0_0_yz_zz_yz, g_y_0_x_0_0_yz_zz_zz, g_y_yz_xzz_xx, g_y_yz_xzz_xy, g_y_yz_xzz_xz, g_y_yz_xzz_yy, g_y_yz_xzz_yz, g_y_yz_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_yz_zz_xx[i] = 4.0 * g_y_yz_xzz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_zz_xy[i] = 4.0 * g_y_yz_xzz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_zz_xz[i] = 4.0 * g_y_yz_xzz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_zz_yy[i] = 4.0 * g_y_yz_xzz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_zz_yz[i] = 4.0 * g_y_yz_xzz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_zz_zz[i] = 4.0 * g_y_yz_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (828-834)

    #pragma omp simd aligned(g_y_0_x_0_0_zz_xx_xx, g_y_0_x_0_0_zz_xx_xy, g_y_0_x_0_0_zz_xx_xz, g_y_0_x_0_0_zz_xx_yy, g_y_0_x_0_0_zz_xx_yz, g_y_0_x_0_0_zz_xx_zz, g_y_zz_x_xx, g_y_zz_x_xy, g_y_zz_x_xz, g_y_zz_x_yy, g_y_zz_x_yz, g_y_zz_x_zz, g_y_zz_xxx_xx, g_y_zz_xxx_xy, g_y_zz_xxx_xz, g_y_zz_xxx_yy, g_y_zz_xxx_yz, g_y_zz_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_zz_xx_xx[i] = -4.0 * g_y_zz_x_xx[i] * a_exp + 4.0 * g_y_zz_xxx_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_xx_xy[i] = -4.0 * g_y_zz_x_xy[i] * a_exp + 4.0 * g_y_zz_xxx_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_xx_xz[i] = -4.0 * g_y_zz_x_xz[i] * a_exp + 4.0 * g_y_zz_xxx_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_xx_yy[i] = -4.0 * g_y_zz_x_yy[i] * a_exp + 4.0 * g_y_zz_xxx_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_xx_yz[i] = -4.0 * g_y_zz_x_yz[i] * a_exp + 4.0 * g_y_zz_xxx_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_xx_zz[i] = -4.0 * g_y_zz_x_zz[i] * a_exp + 4.0 * g_y_zz_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (834-840)

    #pragma omp simd aligned(g_y_0_x_0_0_zz_xy_xx, g_y_0_x_0_0_zz_xy_xy, g_y_0_x_0_0_zz_xy_xz, g_y_0_x_0_0_zz_xy_yy, g_y_0_x_0_0_zz_xy_yz, g_y_0_x_0_0_zz_xy_zz, g_y_zz_xxy_xx, g_y_zz_xxy_xy, g_y_zz_xxy_xz, g_y_zz_xxy_yy, g_y_zz_xxy_yz, g_y_zz_xxy_zz, g_y_zz_y_xx, g_y_zz_y_xy, g_y_zz_y_xz, g_y_zz_y_yy, g_y_zz_y_yz, g_y_zz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_zz_xy_xx[i] = -2.0 * g_y_zz_y_xx[i] * a_exp + 4.0 * g_y_zz_xxy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_xy_xy[i] = -2.0 * g_y_zz_y_xy[i] * a_exp + 4.0 * g_y_zz_xxy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_xy_xz[i] = -2.0 * g_y_zz_y_xz[i] * a_exp + 4.0 * g_y_zz_xxy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_xy_yy[i] = -2.0 * g_y_zz_y_yy[i] * a_exp + 4.0 * g_y_zz_xxy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_xy_yz[i] = -2.0 * g_y_zz_y_yz[i] * a_exp + 4.0 * g_y_zz_xxy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_xy_zz[i] = -2.0 * g_y_zz_y_zz[i] * a_exp + 4.0 * g_y_zz_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (840-846)

    #pragma omp simd aligned(g_y_0_x_0_0_zz_xz_xx, g_y_0_x_0_0_zz_xz_xy, g_y_0_x_0_0_zz_xz_xz, g_y_0_x_0_0_zz_xz_yy, g_y_0_x_0_0_zz_xz_yz, g_y_0_x_0_0_zz_xz_zz, g_y_zz_xxz_xx, g_y_zz_xxz_xy, g_y_zz_xxz_xz, g_y_zz_xxz_yy, g_y_zz_xxz_yz, g_y_zz_xxz_zz, g_y_zz_z_xx, g_y_zz_z_xy, g_y_zz_z_xz, g_y_zz_z_yy, g_y_zz_z_yz, g_y_zz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_zz_xz_xx[i] = -2.0 * g_y_zz_z_xx[i] * a_exp + 4.0 * g_y_zz_xxz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_xz_xy[i] = -2.0 * g_y_zz_z_xy[i] * a_exp + 4.0 * g_y_zz_xxz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_xz_xz[i] = -2.0 * g_y_zz_z_xz[i] * a_exp + 4.0 * g_y_zz_xxz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_xz_yy[i] = -2.0 * g_y_zz_z_yy[i] * a_exp + 4.0 * g_y_zz_xxz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_xz_yz[i] = -2.0 * g_y_zz_z_yz[i] * a_exp + 4.0 * g_y_zz_xxz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_xz_zz[i] = -2.0 * g_y_zz_z_zz[i] * a_exp + 4.0 * g_y_zz_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (846-852)

    #pragma omp simd aligned(g_y_0_x_0_0_zz_yy_xx, g_y_0_x_0_0_zz_yy_xy, g_y_0_x_0_0_zz_yy_xz, g_y_0_x_0_0_zz_yy_yy, g_y_0_x_0_0_zz_yy_yz, g_y_0_x_0_0_zz_yy_zz, g_y_zz_xyy_xx, g_y_zz_xyy_xy, g_y_zz_xyy_xz, g_y_zz_xyy_yy, g_y_zz_xyy_yz, g_y_zz_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_zz_yy_xx[i] = 4.0 * g_y_zz_xyy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_yy_xy[i] = 4.0 * g_y_zz_xyy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_yy_xz[i] = 4.0 * g_y_zz_xyy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_yy_yy[i] = 4.0 * g_y_zz_xyy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_yy_yz[i] = 4.0 * g_y_zz_xyy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_yy_zz[i] = 4.0 * g_y_zz_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (852-858)

    #pragma omp simd aligned(g_y_0_x_0_0_zz_yz_xx, g_y_0_x_0_0_zz_yz_xy, g_y_0_x_0_0_zz_yz_xz, g_y_0_x_0_0_zz_yz_yy, g_y_0_x_0_0_zz_yz_yz, g_y_0_x_0_0_zz_yz_zz, g_y_zz_xyz_xx, g_y_zz_xyz_xy, g_y_zz_xyz_xz, g_y_zz_xyz_yy, g_y_zz_xyz_yz, g_y_zz_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_zz_yz_xx[i] = 4.0 * g_y_zz_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_yz_xy[i] = 4.0 * g_y_zz_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_yz_xz[i] = 4.0 * g_y_zz_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_yz_yy[i] = 4.0 * g_y_zz_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_yz_yz[i] = 4.0 * g_y_zz_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_yz_zz[i] = 4.0 * g_y_zz_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (858-864)

    #pragma omp simd aligned(g_y_0_x_0_0_zz_zz_xx, g_y_0_x_0_0_zz_zz_xy, g_y_0_x_0_0_zz_zz_xz, g_y_0_x_0_0_zz_zz_yy, g_y_0_x_0_0_zz_zz_yz, g_y_0_x_0_0_zz_zz_zz, g_y_zz_xzz_xx, g_y_zz_xzz_xy, g_y_zz_xzz_xz, g_y_zz_xzz_yy, g_y_zz_xzz_yz, g_y_zz_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_zz_zz_xx[i] = 4.0 * g_y_zz_xzz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_zz_xy[i] = 4.0 * g_y_zz_xzz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_zz_xz[i] = 4.0 * g_y_zz_xzz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_zz_yy[i] = 4.0 * g_y_zz_xzz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_zz_yz[i] = 4.0 * g_y_zz_xzz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_zz_zz[i] = 4.0 * g_y_zz_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (864-870)

    #pragma omp simd aligned(g_y_0_y_0_0_xx_xx_xx, g_y_0_y_0_0_xx_xx_xy, g_y_0_y_0_0_xx_xx_xz, g_y_0_y_0_0_xx_xx_yy, g_y_0_y_0_0_xx_xx_yz, g_y_0_y_0_0_xx_xx_zz, g_y_xx_xxy_xx, g_y_xx_xxy_xy, g_y_xx_xxy_xz, g_y_xx_xxy_yy, g_y_xx_xxy_yz, g_y_xx_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xx_xx_xx[i] = 4.0 * g_y_xx_xxy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_xx_xy[i] = 4.0 * g_y_xx_xxy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_xx_xz[i] = 4.0 * g_y_xx_xxy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_xx_yy[i] = 4.0 * g_y_xx_xxy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_xx_yz[i] = 4.0 * g_y_xx_xxy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_xx_zz[i] = 4.0 * g_y_xx_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (870-876)

    #pragma omp simd aligned(g_y_0_y_0_0_xx_xy_xx, g_y_0_y_0_0_xx_xy_xy, g_y_0_y_0_0_xx_xy_xz, g_y_0_y_0_0_xx_xy_yy, g_y_0_y_0_0_xx_xy_yz, g_y_0_y_0_0_xx_xy_zz, g_y_xx_x_xx, g_y_xx_x_xy, g_y_xx_x_xz, g_y_xx_x_yy, g_y_xx_x_yz, g_y_xx_x_zz, g_y_xx_xyy_xx, g_y_xx_xyy_xy, g_y_xx_xyy_xz, g_y_xx_xyy_yy, g_y_xx_xyy_yz, g_y_xx_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xx_xy_xx[i] = -2.0 * g_y_xx_x_xx[i] * a_exp + 4.0 * g_y_xx_xyy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_xy_xy[i] = -2.0 * g_y_xx_x_xy[i] * a_exp + 4.0 * g_y_xx_xyy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_xy_xz[i] = -2.0 * g_y_xx_x_xz[i] * a_exp + 4.0 * g_y_xx_xyy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_xy_yy[i] = -2.0 * g_y_xx_x_yy[i] * a_exp + 4.0 * g_y_xx_xyy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_xy_yz[i] = -2.0 * g_y_xx_x_yz[i] * a_exp + 4.0 * g_y_xx_xyy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_xy_zz[i] = -2.0 * g_y_xx_x_zz[i] * a_exp + 4.0 * g_y_xx_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (876-882)

    #pragma omp simd aligned(g_y_0_y_0_0_xx_xz_xx, g_y_0_y_0_0_xx_xz_xy, g_y_0_y_0_0_xx_xz_xz, g_y_0_y_0_0_xx_xz_yy, g_y_0_y_0_0_xx_xz_yz, g_y_0_y_0_0_xx_xz_zz, g_y_xx_xyz_xx, g_y_xx_xyz_xy, g_y_xx_xyz_xz, g_y_xx_xyz_yy, g_y_xx_xyz_yz, g_y_xx_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xx_xz_xx[i] = 4.0 * g_y_xx_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_xz_xy[i] = 4.0 * g_y_xx_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_xz_xz[i] = 4.0 * g_y_xx_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_xz_yy[i] = 4.0 * g_y_xx_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_xz_yz[i] = 4.0 * g_y_xx_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_xz_zz[i] = 4.0 * g_y_xx_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (882-888)

    #pragma omp simd aligned(g_y_0_y_0_0_xx_yy_xx, g_y_0_y_0_0_xx_yy_xy, g_y_0_y_0_0_xx_yy_xz, g_y_0_y_0_0_xx_yy_yy, g_y_0_y_0_0_xx_yy_yz, g_y_0_y_0_0_xx_yy_zz, g_y_xx_y_xx, g_y_xx_y_xy, g_y_xx_y_xz, g_y_xx_y_yy, g_y_xx_y_yz, g_y_xx_y_zz, g_y_xx_yyy_xx, g_y_xx_yyy_xy, g_y_xx_yyy_xz, g_y_xx_yyy_yy, g_y_xx_yyy_yz, g_y_xx_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xx_yy_xx[i] = -4.0 * g_y_xx_y_xx[i] * a_exp + 4.0 * g_y_xx_yyy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_yy_xy[i] = -4.0 * g_y_xx_y_xy[i] * a_exp + 4.0 * g_y_xx_yyy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_yy_xz[i] = -4.0 * g_y_xx_y_xz[i] * a_exp + 4.0 * g_y_xx_yyy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_yy_yy[i] = -4.0 * g_y_xx_y_yy[i] * a_exp + 4.0 * g_y_xx_yyy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_yy_yz[i] = -4.0 * g_y_xx_y_yz[i] * a_exp + 4.0 * g_y_xx_yyy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_yy_zz[i] = -4.0 * g_y_xx_y_zz[i] * a_exp + 4.0 * g_y_xx_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (888-894)

    #pragma omp simd aligned(g_y_0_y_0_0_xx_yz_xx, g_y_0_y_0_0_xx_yz_xy, g_y_0_y_0_0_xx_yz_xz, g_y_0_y_0_0_xx_yz_yy, g_y_0_y_0_0_xx_yz_yz, g_y_0_y_0_0_xx_yz_zz, g_y_xx_yyz_xx, g_y_xx_yyz_xy, g_y_xx_yyz_xz, g_y_xx_yyz_yy, g_y_xx_yyz_yz, g_y_xx_yyz_zz, g_y_xx_z_xx, g_y_xx_z_xy, g_y_xx_z_xz, g_y_xx_z_yy, g_y_xx_z_yz, g_y_xx_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xx_yz_xx[i] = -2.0 * g_y_xx_z_xx[i] * a_exp + 4.0 * g_y_xx_yyz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_yz_xy[i] = -2.0 * g_y_xx_z_xy[i] * a_exp + 4.0 * g_y_xx_yyz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_yz_xz[i] = -2.0 * g_y_xx_z_xz[i] * a_exp + 4.0 * g_y_xx_yyz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_yz_yy[i] = -2.0 * g_y_xx_z_yy[i] * a_exp + 4.0 * g_y_xx_yyz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_yz_yz[i] = -2.0 * g_y_xx_z_yz[i] * a_exp + 4.0 * g_y_xx_yyz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_yz_zz[i] = -2.0 * g_y_xx_z_zz[i] * a_exp + 4.0 * g_y_xx_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (894-900)

    #pragma omp simd aligned(g_y_0_y_0_0_xx_zz_xx, g_y_0_y_0_0_xx_zz_xy, g_y_0_y_0_0_xx_zz_xz, g_y_0_y_0_0_xx_zz_yy, g_y_0_y_0_0_xx_zz_yz, g_y_0_y_0_0_xx_zz_zz, g_y_xx_yzz_xx, g_y_xx_yzz_xy, g_y_xx_yzz_xz, g_y_xx_yzz_yy, g_y_xx_yzz_yz, g_y_xx_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xx_zz_xx[i] = 4.0 * g_y_xx_yzz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_zz_xy[i] = 4.0 * g_y_xx_yzz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_zz_xz[i] = 4.0 * g_y_xx_yzz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_zz_yy[i] = 4.0 * g_y_xx_yzz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_zz_yz[i] = 4.0 * g_y_xx_yzz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_zz_zz[i] = 4.0 * g_y_xx_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (900-906)

    #pragma omp simd aligned(g_y_0_y_0_0_xy_xx_xx, g_y_0_y_0_0_xy_xx_xy, g_y_0_y_0_0_xy_xx_xz, g_y_0_y_0_0_xy_xx_yy, g_y_0_y_0_0_xy_xx_yz, g_y_0_y_0_0_xy_xx_zz, g_y_xy_xxy_xx, g_y_xy_xxy_xy, g_y_xy_xxy_xz, g_y_xy_xxy_yy, g_y_xy_xxy_yz, g_y_xy_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xy_xx_xx[i] = 4.0 * g_y_xy_xxy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_xx_xy[i] = 4.0 * g_y_xy_xxy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_xx_xz[i] = 4.0 * g_y_xy_xxy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_xx_yy[i] = 4.0 * g_y_xy_xxy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_xx_yz[i] = 4.0 * g_y_xy_xxy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_xx_zz[i] = 4.0 * g_y_xy_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (906-912)

    #pragma omp simd aligned(g_y_0_y_0_0_xy_xy_xx, g_y_0_y_0_0_xy_xy_xy, g_y_0_y_0_0_xy_xy_xz, g_y_0_y_0_0_xy_xy_yy, g_y_0_y_0_0_xy_xy_yz, g_y_0_y_0_0_xy_xy_zz, g_y_xy_x_xx, g_y_xy_x_xy, g_y_xy_x_xz, g_y_xy_x_yy, g_y_xy_x_yz, g_y_xy_x_zz, g_y_xy_xyy_xx, g_y_xy_xyy_xy, g_y_xy_xyy_xz, g_y_xy_xyy_yy, g_y_xy_xyy_yz, g_y_xy_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xy_xy_xx[i] = -2.0 * g_y_xy_x_xx[i] * a_exp + 4.0 * g_y_xy_xyy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_xy_xy[i] = -2.0 * g_y_xy_x_xy[i] * a_exp + 4.0 * g_y_xy_xyy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_xy_xz[i] = -2.0 * g_y_xy_x_xz[i] * a_exp + 4.0 * g_y_xy_xyy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_xy_yy[i] = -2.0 * g_y_xy_x_yy[i] * a_exp + 4.0 * g_y_xy_xyy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_xy_yz[i] = -2.0 * g_y_xy_x_yz[i] * a_exp + 4.0 * g_y_xy_xyy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_xy_zz[i] = -2.0 * g_y_xy_x_zz[i] * a_exp + 4.0 * g_y_xy_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (912-918)

    #pragma omp simd aligned(g_y_0_y_0_0_xy_xz_xx, g_y_0_y_0_0_xy_xz_xy, g_y_0_y_0_0_xy_xz_xz, g_y_0_y_0_0_xy_xz_yy, g_y_0_y_0_0_xy_xz_yz, g_y_0_y_0_0_xy_xz_zz, g_y_xy_xyz_xx, g_y_xy_xyz_xy, g_y_xy_xyz_xz, g_y_xy_xyz_yy, g_y_xy_xyz_yz, g_y_xy_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xy_xz_xx[i] = 4.0 * g_y_xy_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_xz_xy[i] = 4.0 * g_y_xy_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_xz_xz[i] = 4.0 * g_y_xy_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_xz_yy[i] = 4.0 * g_y_xy_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_xz_yz[i] = 4.0 * g_y_xy_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_xz_zz[i] = 4.0 * g_y_xy_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (918-924)

    #pragma omp simd aligned(g_y_0_y_0_0_xy_yy_xx, g_y_0_y_0_0_xy_yy_xy, g_y_0_y_0_0_xy_yy_xz, g_y_0_y_0_0_xy_yy_yy, g_y_0_y_0_0_xy_yy_yz, g_y_0_y_0_0_xy_yy_zz, g_y_xy_y_xx, g_y_xy_y_xy, g_y_xy_y_xz, g_y_xy_y_yy, g_y_xy_y_yz, g_y_xy_y_zz, g_y_xy_yyy_xx, g_y_xy_yyy_xy, g_y_xy_yyy_xz, g_y_xy_yyy_yy, g_y_xy_yyy_yz, g_y_xy_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xy_yy_xx[i] = -4.0 * g_y_xy_y_xx[i] * a_exp + 4.0 * g_y_xy_yyy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_yy_xy[i] = -4.0 * g_y_xy_y_xy[i] * a_exp + 4.0 * g_y_xy_yyy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_yy_xz[i] = -4.0 * g_y_xy_y_xz[i] * a_exp + 4.0 * g_y_xy_yyy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_yy_yy[i] = -4.0 * g_y_xy_y_yy[i] * a_exp + 4.0 * g_y_xy_yyy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_yy_yz[i] = -4.0 * g_y_xy_y_yz[i] * a_exp + 4.0 * g_y_xy_yyy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_yy_zz[i] = -4.0 * g_y_xy_y_zz[i] * a_exp + 4.0 * g_y_xy_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (924-930)

    #pragma omp simd aligned(g_y_0_y_0_0_xy_yz_xx, g_y_0_y_0_0_xy_yz_xy, g_y_0_y_0_0_xy_yz_xz, g_y_0_y_0_0_xy_yz_yy, g_y_0_y_0_0_xy_yz_yz, g_y_0_y_0_0_xy_yz_zz, g_y_xy_yyz_xx, g_y_xy_yyz_xy, g_y_xy_yyz_xz, g_y_xy_yyz_yy, g_y_xy_yyz_yz, g_y_xy_yyz_zz, g_y_xy_z_xx, g_y_xy_z_xy, g_y_xy_z_xz, g_y_xy_z_yy, g_y_xy_z_yz, g_y_xy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xy_yz_xx[i] = -2.0 * g_y_xy_z_xx[i] * a_exp + 4.0 * g_y_xy_yyz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_yz_xy[i] = -2.0 * g_y_xy_z_xy[i] * a_exp + 4.0 * g_y_xy_yyz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_yz_xz[i] = -2.0 * g_y_xy_z_xz[i] * a_exp + 4.0 * g_y_xy_yyz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_yz_yy[i] = -2.0 * g_y_xy_z_yy[i] * a_exp + 4.0 * g_y_xy_yyz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_yz_yz[i] = -2.0 * g_y_xy_z_yz[i] * a_exp + 4.0 * g_y_xy_yyz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_yz_zz[i] = -2.0 * g_y_xy_z_zz[i] * a_exp + 4.0 * g_y_xy_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (930-936)

    #pragma omp simd aligned(g_y_0_y_0_0_xy_zz_xx, g_y_0_y_0_0_xy_zz_xy, g_y_0_y_0_0_xy_zz_xz, g_y_0_y_0_0_xy_zz_yy, g_y_0_y_0_0_xy_zz_yz, g_y_0_y_0_0_xy_zz_zz, g_y_xy_yzz_xx, g_y_xy_yzz_xy, g_y_xy_yzz_xz, g_y_xy_yzz_yy, g_y_xy_yzz_yz, g_y_xy_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xy_zz_xx[i] = 4.0 * g_y_xy_yzz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_zz_xy[i] = 4.0 * g_y_xy_yzz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_zz_xz[i] = 4.0 * g_y_xy_yzz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_zz_yy[i] = 4.0 * g_y_xy_yzz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_zz_yz[i] = 4.0 * g_y_xy_yzz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_zz_zz[i] = 4.0 * g_y_xy_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (936-942)

    #pragma omp simd aligned(g_y_0_y_0_0_xz_xx_xx, g_y_0_y_0_0_xz_xx_xy, g_y_0_y_0_0_xz_xx_xz, g_y_0_y_0_0_xz_xx_yy, g_y_0_y_0_0_xz_xx_yz, g_y_0_y_0_0_xz_xx_zz, g_y_xz_xxy_xx, g_y_xz_xxy_xy, g_y_xz_xxy_xz, g_y_xz_xxy_yy, g_y_xz_xxy_yz, g_y_xz_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xz_xx_xx[i] = 4.0 * g_y_xz_xxy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_xx_xy[i] = 4.0 * g_y_xz_xxy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_xx_xz[i] = 4.0 * g_y_xz_xxy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_xx_yy[i] = 4.0 * g_y_xz_xxy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_xx_yz[i] = 4.0 * g_y_xz_xxy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_xx_zz[i] = 4.0 * g_y_xz_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (942-948)

    #pragma omp simd aligned(g_y_0_y_0_0_xz_xy_xx, g_y_0_y_0_0_xz_xy_xy, g_y_0_y_0_0_xz_xy_xz, g_y_0_y_0_0_xz_xy_yy, g_y_0_y_0_0_xz_xy_yz, g_y_0_y_0_0_xz_xy_zz, g_y_xz_x_xx, g_y_xz_x_xy, g_y_xz_x_xz, g_y_xz_x_yy, g_y_xz_x_yz, g_y_xz_x_zz, g_y_xz_xyy_xx, g_y_xz_xyy_xy, g_y_xz_xyy_xz, g_y_xz_xyy_yy, g_y_xz_xyy_yz, g_y_xz_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xz_xy_xx[i] = -2.0 * g_y_xz_x_xx[i] * a_exp + 4.0 * g_y_xz_xyy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_xy_xy[i] = -2.0 * g_y_xz_x_xy[i] * a_exp + 4.0 * g_y_xz_xyy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_xy_xz[i] = -2.0 * g_y_xz_x_xz[i] * a_exp + 4.0 * g_y_xz_xyy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_xy_yy[i] = -2.0 * g_y_xz_x_yy[i] * a_exp + 4.0 * g_y_xz_xyy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_xy_yz[i] = -2.0 * g_y_xz_x_yz[i] * a_exp + 4.0 * g_y_xz_xyy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_xy_zz[i] = -2.0 * g_y_xz_x_zz[i] * a_exp + 4.0 * g_y_xz_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (948-954)

    #pragma omp simd aligned(g_y_0_y_0_0_xz_xz_xx, g_y_0_y_0_0_xz_xz_xy, g_y_0_y_0_0_xz_xz_xz, g_y_0_y_0_0_xz_xz_yy, g_y_0_y_0_0_xz_xz_yz, g_y_0_y_0_0_xz_xz_zz, g_y_xz_xyz_xx, g_y_xz_xyz_xy, g_y_xz_xyz_xz, g_y_xz_xyz_yy, g_y_xz_xyz_yz, g_y_xz_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xz_xz_xx[i] = 4.0 * g_y_xz_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_xz_xy[i] = 4.0 * g_y_xz_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_xz_xz[i] = 4.0 * g_y_xz_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_xz_yy[i] = 4.0 * g_y_xz_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_xz_yz[i] = 4.0 * g_y_xz_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_xz_zz[i] = 4.0 * g_y_xz_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (954-960)

    #pragma omp simd aligned(g_y_0_y_0_0_xz_yy_xx, g_y_0_y_0_0_xz_yy_xy, g_y_0_y_0_0_xz_yy_xz, g_y_0_y_0_0_xz_yy_yy, g_y_0_y_0_0_xz_yy_yz, g_y_0_y_0_0_xz_yy_zz, g_y_xz_y_xx, g_y_xz_y_xy, g_y_xz_y_xz, g_y_xz_y_yy, g_y_xz_y_yz, g_y_xz_y_zz, g_y_xz_yyy_xx, g_y_xz_yyy_xy, g_y_xz_yyy_xz, g_y_xz_yyy_yy, g_y_xz_yyy_yz, g_y_xz_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xz_yy_xx[i] = -4.0 * g_y_xz_y_xx[i] * a_exp + 4.0 * g_y_xz_yyy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_yy_xy[i] = -4.0 * g_y_xz_y_xy[i] * a_exp + 4.0 * g_y_xz_yyy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_yy_xz[i] = -4.0 * g_y_xz_y_xz[i] * a_exp + 4.0 * g_y_xz_yyy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_yy_yy[i] = -4.0 * g_y_xz_y_yy[i] * a_exp + 4.0 * g_y_xz_yyy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_yy_yz[i] = -4.0 * g_y_xz_y_yz[i] * a_exp + 4.0 * g_y_xz_yyy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_yy_zz[i] = -4.0 * g_y_xz_y_zz[i] * a_exp + 4.0 * g_y_xz_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (960-966)

    #pragma omp simd aligned(g_y_0_y_0_0_xz_yz_xx, g_y_0_y_0_0_xz_yz_xy, g_y_0_y_0_0_xz_yz_xz, g_y_0_y_0_0_xz_yz_yy, g_y_0_y_0_0_xz_yz_yz, g_y_0_y_0_0_xz_yz_zz, g_y_xz_yyz_xx, g_y_xz_yyz_xy, g_y_xz_yyz_xz, g_y_xz_yyz_yy, g_y_xz_yyz_yz, g_y_xz_yyz_zz, g_y_xz_z_xx, g_y_xz_z_xy, g_y_xz_z_xz, g_y_xz_z_yy, g_y_xz_z_yz, g_y_xz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xz_yz_xx[i] = -2.0 * g_y_xz_z_xx[i] * a_exp + 4.0 * g_y_xz_yyz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_yz_xy[i] = -2.0 * g_y_xz_z_xy[i] * a_exp + 4.0 * g_y_xz_yyz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_yz_xz[i] = -2.0 * g_y_xz_z_xz[i] * a_exp + 4.0 * g_y_xz_yyz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_yz_yy[i] = -2.0 * g_y_xz_z_yy[i] * a_exp + 4.0 * g_y_xz_yyz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_yz_yz[i] = -2.0 * g_y_xz_z_yz[i] * a_exp + 4.0 * g_y_xz_yyz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_yz_zz[i] = -2.0 * g_y_xz_z_zz[i] * a_exp + 4.0 * g_y_xz_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (966-972)

    #pragma omp simd aligned(g_y_0_y_0_0_xz_zz_xx, g_y_0_y_0_0_xz_zz_xy, g_y_0_y_0_0_xz_zz_xz, g_y_0_y_0_0_xz_zz_yy, g_y_0_y_0_0_xz_zz_yz, g_y_0_y_0_0_xz_zz_zz, g_y_xz_yzz_xx, g_y_xz_yzz_xy, g_y_xz_yzz_xz, g_y_xz_yzz_yy, g_y_xz_yzz_yz, g_y_xz_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xz_zz_xx[i] = 4.0 * g_y_xz_yzz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_zz_xy[i] = 4.0 * g_y_xz_yzz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_zz_xz[i] = 4.0 * g_y_xz_yzz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_zz_yy[i] = 4.0 * g_y_xz_yzz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_zz_yz[i] = 4.0 * g_y_xz_yzz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_zz_zz[i] = 4.0 * g_y_xz_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (972-978)

    #pragma omp simd aligned(g_y_0_y_0_0_yy_xx_xx, g_y_0_y_0_0_yy_xx_xy, g_y_0_y_0_0_yy_xx_xz, g_y_0_y_0_0_yy_xx_yy, g_y_0_y_0_0_yy_xx_yz, g_y_0_y_0_0_yy_xx_zz, g_y_yy_xxy_xx, g_y_yy_xxy_xy, g_y_yy_xxy_xz, g_y_yy_xxy_yy, g_y_yy_xxy_yz, g_y_yy_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_yy_xx_xx[i] = 4.0 * g_y_yy_xxy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_xx_xy[i] = 4.0 * g_y_yy_xxy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_xx_xz[i] = 4.0 * g_y_yy_xxy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_xx_yy[i] = 4.0 * g_y_yy_xxy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_xx_yz[i] = 4.0 * g_y_yy_xxy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_xx_zz[i] = 4.0 * g_y_yy_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (978-984)

    #pragma omp simd aligned(g_y_0_y_0_0_yy_xy_xx, g_y_0_y_0_0_yy_xy_xy, g_y_0_y_0_0_yy_xy_xz, g_y_0_y_0_0_yy_xy_yy, g_y_0_y_0_0_yy_xy_yz, g_y_0_y_0_0_yy_xy_zz, g_y_yy_x_xx, g_y_yy_x_xy, g_y_yy_x_xz, g_y_yy_x_yy, g_y_yy_x_yz, g_y_yy_x_zz, g_y_yy_xyy_xx, g_y_yy_xyy_xy, g_y_yy_xyy_xz, g_y_yy_xyy_yy, g_y_yy_xyy_yz, g_y_yy_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_yy_xy_xx[i] = -2.0 * g_y_yy_x_xx[i] * a_exp + 4.0 * g_y_yy_xyy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_xy_xy[i] = -2.0 * g_y_yy_x_xy[i] * a_exp + 4.0 * g_y_yy_xyy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_xy_xz[i] = -2.0 * g_y_yy_x_xz[i] * a_exp + 4.0 * g_y_yy_xyy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_xy_yy[i] = -2.0 * g_y_yy_x_yy[i] * a_exp + 4.0 * g_y_yy_xyy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_xy_yz[i] = -2.0 * g_y_yy_x_yz[i] * a_exp + 4.0 * g_y_yy_xyy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_xy_zz[i] = -2.0 * g_y_yy_x_zz[i] * a_exp + 4.0 * g_y_yy_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (984-990)

    #pragma omp simd aligned(g_y_0_y_0_0_yy_xz_xx, g_y_0_y_0_0_yy_xz_xy, g_y_0_y_0_0_yy_xz_xz, g_y_0_y_0_0_yy_xz_yy, g_y_0_y_0_0_yy_xz_yz, g_y_0_y_0_0_yy_xz_zz, g_y_yy_xyz_xx, g_y_yy_xyz_xy, g_y_yy_xyz_xz, g_y_yy_xyz_yy, g_y_yy_xyz_yz, g_y_yy_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_yy_xz_xx[i] = 4.0 * g_y_yy_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_xz_xy[i] = 4.0 * g_y_yy_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_xz_xz[i] = 4.0 * g_y_yy_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_xz_yy[i] = 4.0 * g_y_yy_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_xz_yz[i] = 4.0 * g_y_yy_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_xz_zz[i] = 4.0 * g_y_yy_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (990-996)

    #pragma omp simd aligned(g_y_0_y_0_0_yy_yy_xx, g_y_0_y_0_0_yy_yy_xy, g_y_0_y_0_0_yy_yy_xz, g_y_0_y_0_0_yy_yy_yy, g_y_0_y_0_0_yy_yy_yz, g_y_0_y_0_0_yy_yy_zz, g_y_yy_y_xx, g_y_yy_y_xy, g_y_yy_y_xz, g_y_yy_y_yy, g_y_yy_y_yz, g_y_yy_y_zz, g_y_yy_yyy_xx, g_y_yy_yyy_xy, g_y_yy_yyy_xz, g_y_yy_yyy_yy, g_y_yy_yyy_yz, g_y_yy_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_yy_yy_xx[i] = -4.0 * g_y_yy_y_xx[i] * a_exp + 4.0 * g_y_yy_yyy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_yy_xy[i] = -4.0 * g_y_yy_y_xy[i] * a_exp + 4.0 * g_y_yy_yyy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_yy_xz[i] = -4.0 * g_y_yy_y_xz[i] * a_exp + 4.0 * g_y_yy_yyy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_yy_yy[i] = -4.0 * g_y_yy_y_yy[i] * a_exp + 4.0 * g_y_yy_yyy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_yy_yz[i] = -4.0 * g_y_yy_y_yz[i] * a_exp + 4.0 * g_y_yy_yyy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_yy_zz[i] = -4.0 * g_y_yy_y_zz[i] * a_exp + 4.0 * g_y_yy_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (996-1002)

    #pragma omp simd aligned(g_y_0_y_0_0_yy_yz_xx, g_y_0_y_0_0_yy_yz_xy, g_y_0_y_0_0_yy_yz_xz, g_y_0_y_0_0_yy_yz_yy, g_y_0_y_0_0_yy_yz_yz, g_y_0_y_0_0_yy_yz_zz, g_y_yy_yyz_xx, g_y_yy_yyz_xy, g_y_yy_yyz_xz, g_y_yy_yyz_yy, g_y_yy_yyz_yz, g_y_yy_yyz_zz, g_y_yy_z_xx, g_y_yy_z_xy, g_y_yy_z_xz, g_y_yy_z_yy, g_y_yy_z_yz, g_y_yy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_yy_yz_xx[i] = -2.0 * g_y_yy_z_xx[i] * a_exp + 4.0 * g_y_yy_yyz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_yz_xy[i] = -2.0 * g_y_yy_z_xy[i] * a_exp + 4.0 * g_y_yy_yyz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_yz_xz[i] = -2.0 * g_y_yy_z_xz[i] * a_exp + 4.0 * g_y_yy_yyz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_yz_yy[i] = -2.0 * g_y_yy_z_yy[i] * a_exp + 4.0 * g_y_yy_yyz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_yz_yz[i] = -2.0 * g_y_yy_z_yz[i] * a_exp + 4.0 * g_y_yy_yyz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_yz_zz[i] = -2.0 * g_y_yy_z_zz[i] * a_exp + 4.0 * g_y_yy_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1002-1008)

    #pragma omp simd aligned(g_y_0_y_0_0_yy_zz_xx, g_y_0_y_0_0_yy_zz_xy, g_y_0_y_0_0_yy_zz_xz, g_y_0_y_0_0_yy_zz_yy, g_y_0_y_0_0_yy_zz_yz, g_y_0_y_0_0_yy_zz_zz, g_y_yy_yzz_xx, g_y_yy_yzz_xy, g_y_yy_yzz_xz, g_y_yy_yzz_yy, g_y_yy_yzz_yz, g_y_yy_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_yy_zz_xx[i] = 4.0 * g_y_yy_yzz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_zz_xy[i] = 4.0 * g_y_yy_yzz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_zz_xz[i] = 4.0 * g_y_yy_yzz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_zz_yy[i] = 4.0 * g_y_yy_yzz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_zz_yz[i] = 4.0 * g_y_yy_yzz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_zz_zz[i] = 4.0 * g_y_yy_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1008-1014)

    #pragma omp simd aligned(g_y_0_y_0_0_yz_xx_xx, g_y_0_y_0_0_yz_xx_xy, g_y_0_y_0_0_yz_xx_xz, g_y_0_y_0_0_yz_xx_yy, g_y_0_y_0_0_yz_xx_yz, g_y_0_y_0_0_yz_xx_zz, g_y_yz_xxy_xx, g_y_yz_xxy_xy, g_y_yz_xxy_xz, g_y_yz_xxy_yy, g_y_yz_xxy_yz, g_y_yz_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_yz_xx_xx[i] = 4.0 * g_y_yz_xxy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_xx_xy[i] = 4.0 * g_y_yz_xxy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_xx_xz[i] = 4.0 * g_y_yz_xxy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_xx_yy[i] = 4.0 * g_y_yz_xxy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_xx_yz[i] = 4.0 * g_y_yz_xxy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_xx_zz[i] = 4.0 * g_y_yz_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1014-1020)

    #pragma omp simd aligned(g_y_0_y_0_0_yz_xy_xx, g_y_0_y_0_0_yz_xy_xy, g_y_0_y_0_0_yz_xy_xz, g_y_0_y_0_0_yz_xy_yy, g_y_0_y_0_0_yz_xy_yz, g_y_0_y_0_0_yz_xy_zz, g_y_yz_x_xx, g_y_yz_x_xy, g_y_yz_x_xz, g_y_yz_x_yy, g_y_yz_x_yz, g_y_yz_x_zz, g_y_yz_xyy_xx, g_y_yz_xyy_xy, g_y_yz_xyy_xz, g_y_yz_xyy_yy, g_y_yz_xyy_yz, g_y_yz_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_yz_xy_xx[i] = -2.0 * g_y_yz_x_xx[i] * a_exp + 4.0 * g_y_yz_xyy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_xy_xy[i] = -2.0 * g_y_yz_x_xy[i] * a_exp + 4.0 * g_y_yz_xyy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_xy_xz[i] = -2.0 * g_y_yz_x_xz[i] * a_exp + 4.0 * g_y_yz_xyy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_xy_yy[i] = -2.0 * g_y_yz_x_yy[i] * a_exp + 4.0 * g_y_yz_xyy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_xy_yz[i] = -2.0 * g_y_yz_x_yz[i] * a_exp + 4.0 * g_y_yz_xyy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_xy_zz[i] = -2.0 * g_y_yz_x_zz[i] * a_exp + 4.0 * g_y_yz_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1020-1026)

    #pragma omp simd aligned(g_y_0_y_0_0_yz_xz_xx, g_y_0_y_0_0_yz_xz_xy, g_y_0_y_0_0_yz_xz_xz, g_y_0_y_0_0_yz_xz_yy, g_y_0_y_0_0_yz_xz_yz, g_y_0_y_0_0_yz_xz_zz, g_y_yz_xyz_xx, g_y_yz_xyz_xy, g_y_yz_xyz_xz, g_y_yz_xyz_yy, g_y_yz_xyz_yz, g_y_yz_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_yz_xz_xx[i] = 4.0 * g_y_yz_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_xz_xy[i] = 4.0 * g_y_yz_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_xz_xz[i] = 4.0 * g_y_yz_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_xz_yy[i] = 4.0 * g_y_yz_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_xz_yz[i] = 4.0 * g_y_yz_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_xz_zz[i] = 4.0 * g_y_yz_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1026-1032)

    #pragma omp simd aligned(g_y_0_y_0_0_yz_yy_xx, g_y_0_y_0_0_yz_yy_xy, g_y_0_y_0_0_yz_yy_xz, g_y_0_y_0_0_yz_yy_yy, g_y_0_y_0_0_yz_yy_yz, g_y_0_y_0_0_yz_yy_zz, g_y_yz_y_xx, g_y_yz_y_xy, g_y_yz_y_xz, g_y_yz_y_yy, g_y_yz_y_yz, g_y_yz_y_zz, g_y_yz_yyy_xx, g_y_yz_yyy_xy, g_y_yz_yyy_xz, g_y_yz_yyy_yy, g_y_yz_yyy_yz, g_y_yz_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_yz_yy_xx[i] = -4.0 * g_y_yz_y_xx[i] * a_exp + 4.0 * g_y_yz_yyy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_yy_xy[i] = -4.0 * g_y_yz_y_xy[i] * a_exp + 4.0 * g_y_yz_yyy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_yy_xz[i] = -4.0 * g_y_yz_y_xz[i] * a_exp + 4.0 * g_y_yz_yyy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_yy_yy[i] = -4.0 * g_y_yz_y_yy[i] * a_exp + 4.0 * g_y_yz_yyy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_yy_yz[i] = -4.0 * g_y_yz_y_yz[i] * a_exp + 4.0 * g_y_yz_yyy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_yy_zz[i] = -4.0 * g_y_yz_y_zz[i] * a_exp + 4.0 * g_y_yz_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1032-1038)

    #pragma omp simd aligned(g_y_0_y_0_0_yz_yz_xx, g_y_0_y_0_0_yz_yz_xy, g_y_0_y_0_0_yz_yz_xz, g_y_0_y_0_0_yz_yz_yy, g_y_0_y_0_0_yz_yz_yz, g_y_0_y_0_0_yz_yz_zz, g_y_yz_yyz_xx, g_y_yz_yyz_xy, g_y_yz_yyz_xz, g_y_yz_yyz_yy, g_y_yz_yyz_yz, g_y_yz_yyz_zz, g_y_yz_z_xx, g_y_yz_z_xy, g_y_yz_z_xz, g_y_yz_z_yy, g_y_yz_z_yz, g_y_yz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_yz_yz_xx[i] = -2.0 * g_y_yz_z_xx[i] * a_exp + 4.0 * g_y_yz_yyz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_yz_xy[i] = -2.0 * g_y_yz_z_xy[i] * a_exp + 4.0 * g_y_yz_yyz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_yz_xz[i] = -2.0 * g_y_yz_z_xz[i] * a_exp + 4.0 * g_y_yz_yyz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_yz_yy[i] = -2.0 * g_y_yz_z_yy[i] * a_exp + 4.0 * g_y_yz_yyz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_yz_yz[i] = -2.0 * g_y_yz_z_yz[i] * a_exp + 4.0 * g_y_yz_yyz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_yz_zz[i] = -2.0 * g_y_yz_z_zz[i] * a_exp + 4.0 * g_y_yz_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1038-1044)

    #pragma omp simd aligned(g_y_0_y_0_0_yz_zz_xx, g_y_0_y_0_0_yz_zz_xy, g_y_0_y_0_0_yz_zz_xz, g_y_0_y_0_0_yz_zz_yy, g_y_0_y_0_0_yz_zz_yz, g_y_0_y_0_0_yz_zz_zz, g_y_yz_yzz_xx, g_y_yz_yzz_xy, g_y_yz_yzz_xz, g_y_yz_yzz_yy, g_y_yz_yzz_yz, g_y_yz_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_yz_zz_xx[i] = 4.0 * g_y_yz_yzz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_zz_xy[i] = 4.0 * g_y_yz_yzz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_zz_xz[i] = 4.0 * g_y_yz_yzz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_zz_yy[i] = 4.0 * g_y_yz_yzz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_zz_yz[i] = 4.0 * g_y_yz_yzz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_zz_zz[i] = 4.0 * g_y_yz_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1044-1050)

    #pragma omp simd aligned(g_y_0_y_0_0_zz_xx_xx, g_y_0_y_0_0_zz_xx_xy, g_y_0_y_0_0_zz_xx_xz, g_y_0_y_0_0_zz_xx_yy, g_y_0_y_0_0_zz_xx_yz, g_y_0_y_0_0_zz_xx_zz, g_y_zz_xxy_xx, g_y_zz_xxy_xy, g_y_zz_xxy_xz, g_y_zz_xxy_yy, g_y_zz_xxy_yz, g_y_zz_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_zz_xx_xx[i] = 4.0 * g_y_zz_xxy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_xx_xy[i] = 4.0 * g_y_zz_xxy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_xx_xz[i] = 4.0 * g_y_zz_xxy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_xx_yy[i] = 4.0 * g_y_zz_xxy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_xx_yz[i] = 4.0 * g_y_zz_xxy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_xx_zz[i] = 4.0 * g_y_zz_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1050-1056)

    #pragma omp simd aligned(g_y_0_y_0_0_zz_xy_xx, g_y_0_y_0_0_zz_xy_xy, g_y_0_y_0_0_zz_xy_xz, g_y_0_y_0_0_zz_xy_yy, g_y_0_y_0_0_zz_xy_yz, g_y_0_y_0_0_zz_xy_zz, g_y_zz_x_xx, g_y_zz_x_xy, g_y_zz_x_xz, g_y_zz_x_yy, g_y_zz_x_yz, g_y_zz_x_zz, g_y_zz_xyy_xx, g_y_zz_xyy_xy, g_y_zz_xyy_xz, g_y_zz_xyy_yy, g_y_zz_xyy_yz, g_y_zz_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_zz_xy_xx[i] = -2.0 * g_y_zz_x_xx[i] * a_exp + 4.0 * g_y_zz_xyy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_xy_xy[i] = -2.0 * g_y_zz_x_xy[i] * a_exp + 4.0 * g_y_zz_xyy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_xy_xz[i] = -2.0 * g_y_zz_x_xz[i] * a_exp + 4.0 * g_y_zz_xyy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_xy_yy[i] = -2.0 * g_y_zz_x_yy[i] * a_exp + 4.0 * g_y_zz_xyy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_xy_yz[i] = -2.0 * g_y_zz_x_yz[i] * a_exp + 4.0 * g_y_zz_xyy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_xy_zz[i] = -2.0 * g_y_zz_x_zz[i] * a_exp + 4.0 * g_y_zz_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1056-1062)

    #pragma omp simd aligned(g_y_0_y_0_0_zz_xz_xx, g_y_0_y_0_0_zz_xz_xy, g_y_0_y_0_0_zz_xz_xz, g_y_0_y_0_0_zz_xz_yy, g_y_0_y_0_0_zz_xz_yz, g_y_0_y_0_0_zz_xz_zz, g_y_zz_xyz_xx, g_y_zz_xyz_xy, g_y_zz_xyz_xz, g_y_zz_xyz_yy, g_y_zz_xyz_yz, g_y_zz_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_zz_xz_xx[i] = 4.0 * g_y_zz_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_xz_xy[i] = 4.0 * g_y_zz_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_xz_xz[i] = 4.0 * g_y_zz_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_xz_yy[i] = 4.0 * g_y_zz_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_xz_yz[i] = 4.0 * g_y_zz_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_xz_zz[i] = 4.0 * g_y_zz_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1062-1068)

    #pragma omp simd aligned(g_y_0_y_0_0_zz_yy_xx, g_y_0_y_0_0_zz_yy_xy, g_y_0_y_0_0_zz_yy_xz, g_y_0_y_0_0_zz_yy_yy, g_y_0_y_0_0_zz_yy_yz, g_y_0_y_0_0_zz_yy_zz, g_y_zz_y_xx, g_y_zz_y_xy, g_y_zz_y_xz, g_y_zz_y_yy, g_y_zz_y_yz, g_y_zz_y_zz, g_y_zz_yyy_xx, g_y_zz_yyy_xy, g_y_zz_yyy_xz, g_y_zz_yyy_yy, g_y_zz_yyy_yz, g_y_zz_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_zz_yy_xx[i] = -4.0 * g_y_zz_y_xx[i] * a_exp + 4.0 * g_y_zz_yyy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_yy_xy[i] = -4.0 * g_y_zz_y_xy[i] * a_exp + 4.0 * g_y_zz_yyy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_yy_xz[i] = -4.0 * g_y_zz_y_xz[i] * a_exp + 4.0 * g_y_zz_yyy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_yy_yy[i] = -4.0 * g_y_zz_y_yy[i] * a_exp + 4.0 * g_y_zz_yyy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_yy_yz[i] = -4.0 * g_y_zz_y_yz[i] * a_exp + 4.0 * g_y_zz_yyy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_yy_zz[i] = -4.0 * g_y_zz_y_zz[i] * a_exp + 4.0 * g_y_zz_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1068-1074)

    #pragma omp simd aligned(g_y_0_y_0_0_zz_yz_xx, g_y_0_y_0_0_zz_yz_xy, g_y_0_y_0_0_zz_yz_xz, g_y_0_y_0_0_zz_yz_yy, g_y_0_y_0_0_zz_yz_yz, g_y_0_y_0_0_zz_yz_zz, g_y_zz_yyz_xx, g_y_zz_yyz_xy, g_y_zz_yyz_xz, g_y_zz_yyz_yy, g_y_zz_yyz_yz, g_y_zz_yyz_zz, g_y_zz_z_xx, g_y_zz_z_xy, g_y_zz_z_xz, g_y_zz_z_yy, g_y_zz_z_yz, g_y_zz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_zz_yz_xx[i] = -2.0 * g_y_zz_z_xx[i] * a_exp + 4.0 * g_y_zz_yyz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_yz_xy[i] = -2.0 * g_y_zz_z_xy[i] * a_exp + 4.0 * g_y_zz_yyz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_yz_xz[i] = -2.0 * g_y_zz_z_xz[i] * a_exp + 4.0 * g_y_zz_yyz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_yz_yy[i] = -2.0 * g_y_zz_z_yy[i] * a_exp + 4.0 * g_y_zz_yyz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_yz_yz[i] = -2.0 * g_y_zz_z_yz[i] * a_exp + 4.0 * g_y_zz_yyz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_yz_zz[i] = -2.0 * g_y_zz_z_zz[i] * a_exp + 4.0 * g_y_zz_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1074-1080)

    #pragma omp simd aligned(g_y_0_y_0_0_zz_zz_xx, g_y_0_y_0_0_zz_zz_xy, g_y_0_y_0_0_zz_zz_xz, g_y_0_y_0_0_zz_zz_yy, g_y_0_y_0_0_zz_zz_yz, g_y_0_y_0_0_zz_zz_zz, g_y_zz_yzz_xx, g_y_zz_yzz_xy, g_y_zz_yzz_xz, g_y_zz_yzz_yy, g_y_zz_yzz_yz, g_y_zz_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_zz_zz_xx[i] = 4.0 * g_y_zz_yzz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_zz_xy[i] = 4.0 * g_y_zz_yzz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_zz_xz[i] = 4.0 * g_y_zz_yzz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_zz_yy[i] = 4.0 * g_y_zz_yzz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_zz_yz[i] = 4.0 * g_y_zz_yzz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_zz_zz[i] = 4.0 * g_y_zz_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1080-1086)

    #pragma omp simd aligned(g_y_0_z_0_0_xx_xx_xx, g_y_0_z_0_0_xx_xx_xy, g_y_0_z_0_0_xx_xx_xz, g_y_0_z_0_0_xx_xx_yy, g_y_0_z_0_0_xx_xx_yz, g_y_0_z_0_0_xx_xx_zz, g_y_xx_xxz_xx, g_y_xx_xxz_xy, g_y_xx_xxz_xz, g_y_xx_xxz_yy, g_y_xx_xxz_yz, g_y_xx_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xx_xx_xx[i] = 4.0 * g_y_xx_xxz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_xx_xy[i] = 4.0 * g_y_xx_xxz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_xx_xz[i] = 4.0 * g_y_xx_xxz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_xx_yy[i] = 4.0 * g_y_xx_xxz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_xx_yz[i] = 4.0 * g_y_xx_xxz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_xx_zz[i] = 4.0 * g_y_xx_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1086-1092)

    #pragma omp simd aligned(g_y_0_z_0_0_xx_xy_xx, g_y_0_z_0_0_xx_xy_xy, g_y_0_z_0_0_xx_xy_xz, g_y_0_z_0_0_xx_xy_yy, g_y_0_z_0_0_xx_xy_yz, g_y_0_z_0_0_xx_xy_zz, g_y_xx_xyz_xx, g_y_xx_xyz_xy, g_y_xx_xyz_xz, g_y_xx_xyz_yy, g_y_xx_xyz_yz, g_y_xx_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xx_xy_xx[i] = 4.0 * g_y_xx_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_xy_xy[i] = 4.0 * g_y_xx_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_xy_xz[i] = 4.0 * g_y_xx_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_xy_yy[i] = 4.0 * g_y_xx_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_xy_yz[i] = 4.0 * g_y_xx_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_xy_zz[i] = 4.0 * g_y_xx_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1092-1098)

    #pragma omp simd aligned(g_y_0_z_0_0_xx_xz_xx, g_y_0_z_0_0_xx_xz_xy, g_y_0_z_0_0_xx_xz_xz, g_y_0_z_0_0_xx_xz_yy, g_y_0_z_0_0_xx_xz_yz, g_y_0_z_0_0_xx_xz_zz, g_y_xx_x_xx, g_y_xx_x_xy, g_y_xx_x_xz, g_y_xx_x_yy, g_y_xx_x_yz, g_y_xx_x_zz, g_y_xx_xzz_xx, g_y_xx_xzz_xy, g_y_xx_xzz_xz, g_y_xx_xzz_yy, g_y_xx_xzz_yz, g_y_xx_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xx_xz_xx[i] = -2.0 * g_y_xx_x_xx[i] * a_exp + 4.0 * g_y_xx_xzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_xz_xy[i] = -2.0 * g_y_xx_x_xy[i] * a_exp + 4.0 * g_y_xx_xzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_xz_xz[i] = -2.0 * g_y_xx_x_xz[i] * a_exp + 4.0 * g_y_xx_xzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_xz_yy[i] = -2.0 * g_y_xx_x_yy[i] * a_exp + 4.0 * g_y_xx_xzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_xz_yz[i] = -2.0 * g_y_xx_x_yz[i] * a_exp + 4.0 * g_y_xx_xzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_xz_zz[i] = -2.0 * g_y_xx_x_zz[i] * a_exp + 4.0 * g_y_xx_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1098-1104)

    #pragma omp simd aligned(g_y_0_z_0_0_xx_yy_xx, g_y_0_z_0_0_xx_yy_xy, g_y_0_z_0_0_xx_yy_xz, g_y_0_z_0_0_xx_yy_yy, g_y_0_z_0_0_xx_yy_yz, g_y_0_z_0_0_xx_yy_zz, g_y_xx_yyz_xx, g_y_xx_yyz_xy, g_y_xx_yyz_xz, g_y_xx_yyz_yy, g_y_xx_yyz_yz, g_y_xx_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xx_yy_xx[i] = 4.0 * g_y_xx_yyz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_yy_xy[i] = 4.0 * g_y_xx_yyz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_yy_xz[i] = 4.0 * g_y_xx_yyz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_yy_yy[i] = 4.0 * g_y_xx_yyz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_yy_yz[i] = 4.0 * g_y_xx_yyz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_yy_zz[i] = 4.0 * g_y_xx_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1104-1110)

    #pragma omp simd aligned(g_y_0_z_0_0_xx_yz_xx, g_y_0_z_0_0_xx_yz_xy, g_y_0_z_0_0_xx_yz_xz, g_y_0_z_0_0_xx_yz_yy, g_y_0_z_0_0_xx_yz_yz, g_y_0_z_0_0_xx_yz_zz, g_y_xx_y_xx, g_y_xx_y_xy, g_y_xx_y_xz, g_y_xx_y_yy, g_y_xx_y_yz, g_y_xx_y_zz, g_y_xx_yzz_xx, g_y_xx_yzz_xy, g_y_xx_yzz_xz, g_y_xx_yzz_yy, g_y_xx_yzz_yz, g_y_xx_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xx_yz_xx[i] = -2.0 * g_y_xx_y_xx[i] * a_exp + 4.0 * g_y_xx_yzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_yz_xy[i] = -2.0 * g_y_xx_y_xy[i] * a_exp + 4.0 * g_y_xx_yzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_yz_xz[i] = -2.0 * g_y_xx_y_xz[i] * a_exp + 4.0 * g_y_xx_yzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_yz_yy[i] = -2.0 * g_y_xx_y_yy[i] * a_exp + 4.0 * g_y_xx_yzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_yz_yz[i] = -2.0 * g_y_xx_y_yz[i] * a_exp + 4.0 * g_y_xx_yzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_yz_zz[i] = -2.0 * g_y_xx_y_zz[i] * a_exp + 4.0 * g_y_xx_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1110-1116)

    #pragma omp simd aligned(g_y_0_z_0_0_xx_zz_xx, g_y_0_z_0_0_xx_zz_xy, g_y_0_z_0_0_xx_zz_xz, g_y_0_z_0_0_xx_zz_yy, g_y_0_z_0_0_xx_zz_yz, g_y_0_z_0_0_xx_zz_zz, g_y_xx_z_xx, g_y_xx_z_xy, g_y_xx_z_xz, g_y_xx_z_yy, g_y_xx_z_yz, g_y_xx_z_zz, g_y_xx_zzz_xx, g_y_xx_zzz_xy, g_y_xx_zzz_xz, g_y_xx_zzz_yy, g_y_xx_zzz_yz, g_y_xx_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xx_zz_xx[i] = -4.0 * g_y_xx_z_xx[i] * a_exp + 4.0 * g_y_xx_zzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_zz_xy[i] = -4.0 * g_y_xx_z_xy[i] * a_exp + 4.0 * g_y_xx_zzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_zz_xz[i] = -4.0 * g_y_xx_z_xz[i] * a_exp + 4.0 * g_y_xx_zzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_zz_yy[i] = -4.0 * g_y_xx_z_yy[i] * a_exp + 4.0 * g_y_xx_zzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_zz_yz[i] = -4.0 * g_y_xx_z_yz[i] * a_exp + 4.0 * g_y_xx_zzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_zz_zz[i] = -4.0 * g_y_xx_z_zz[i] * a_exp + 4.0 * g_y_xx_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1116-1122)

    #pragma omp simd aligned(g_y_0_z_0_0_xy_xx_xx, g_y_0_z_0_0_xy_xx_xy, g_y_0_z_0_0_xy_xx_xz, g_y_0_z_0_0_xy_xx_yy, g_y_0_z_0_0_xy_xx_yz, g_y_0_z_0_0_xy_xx_zz, g_y_xy_xxz_xx, g_y_xy_xxz_xy, g_y_xy_xxz_xz, g_y_xy_xxz_yy, g_y_xy_xxz_yz, g_y_xy_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xy_xx_xx[i] = 4.0 * g_y_xy_xxz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_xx_xy[i] = 4.0 * g_y_xy_xxz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_xx_xz[i] = 4.0 * g_y_xy_xxz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_xx_yy[i] = 4.0 * g_y_xy_xxz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_xx_yz[i] = 4.0 * g_y_xy_xxz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_xx_zz[i] = 4.0 * g_y_xy_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1122-1128)

    #pragma omp simd aligned(g_y_0_z_0_0_xy_xy_xx, g_y_0_z_0_0_xy_xy_xy, g_y_0_z_0_0_xy_xy_xz, g_y_0_z_0_0_xy_xy_yy, g_y_0_z_0_0_xy_xy_yz, g_y_0_z_0_0_xy_xy_zz, g_y_xy_xyz_xx, g_y_xy_xyz_xy, g_y_xy_xyz_xz, g_y_xy_xyz_yy, g_y_xy_xyz_yz, g_y_xy_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xy_xy_xx[i] = 4.0 * g_y_xy_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_xy_xy[i] = 4.0 * g_y_xy_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_xy_xz[i] = 4.0 * g_y_xy_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_xy_yy[i] = 4.0 * g_y_xy_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_xy_yz[i] = 4.0 * g_y_xy_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_xy_zz[i] = 4.0 * g_y_xy_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1128-1134)

    #pragma omp simd aligned(g_y_0_z_0_0_xy_xz_xx, g_y_0_z_0_0_xy_xz_xy, g_y_0_z_0_0_xy_xz_xz, g_y_0_z_0_0_xy_xz_yy, g_y_0_z_0_0_xy_xz_yz, g_y_0_z_0_0_xy_xz_zz, g_y_xy_x_xx, g_y_xy_x_xy, g_y_xy_x_xz, g_y_xy_x_yy, g_y_xy_x_yz, g_y_xy_x_zz, g_y_xy_xzz_xx, g_y_xy_xzz_xy, g_y_xy_xzz_xz, g_y_xy_xzz_yy, g_y_xy_xzz_yz, g_y_xy_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xy_xz_xx[i] = -2.0 * g_y_xy_x_xx[i] * a_exp + 4.0 * g_y_xy_xzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_xz_xy[i] = -2.0 * g_y_xy_x_xy[i] * a_exp + 4.0 * g_y_xy_xzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_xz_xz[i] = -2.0 * g_y_xy_x_xz[i] * a_exp + 4.0 * g_y_xy_xzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_xz_yy[i] = -2.0 * g_y_xy_x_yy[i] * a_exp + 4.0 * g_y_xy_xzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_xz_yz[i] = -2.0 * g_y_xy_x_yz[i] * a_exp + 4.0 * g_y_xy_xzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_xz_zz[i] = -2.0 * g_y_xy_x_zz[i] * a_exp + 4.0 * g_y_xy_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1134-1140)

    #pragma omp simd aligned(g_y_0_z_0_0_xy_yy_xx, g_y_0_z_0_0_xy_yy_xy, g_y_0_z_0_0_xy_yy_xz, g_y_0_z_0_0_xy_yy_yy, g_y_0_z_0_0_xy_yy_yz, g_y_0_z_0_0_xy_yy_zz, g_y_xy_yyz_xx, g_y_xy_yyz_xy, g_y_xy_yyz_xz, g_y_xy_yyz_yy, g_y_xy_yyz_yz, g_y_xy_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xy_yy_xx[i] = 4.0 * g_y_xy_yyz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_yy_xy[i] = 4.0 * g_y_xy_yyz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_yy_xz[i] = 4.0 * g_y_xy_yyz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_yy_yy[i] = 4.0 * g_y_xy_yyz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_yy_yz[i] = 4.0 * g_y_xy_yyz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_yy_zz[i] = 4.0 * g_y_xy_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1140-1146)

    #pragma omp simd aligned(g_y_0_z_0_0_xy_yz_xx, g_y_0_z_0_0_xy_yz_xy, g_y_0_z_0_0_xy_yz_xz, g_y_0_z_0_0_xy_yz_yy, g_y_0_z_0_0_xy_yz_yz, g_y_0_z_0_0_xy_yz_zz, g_y_xy_y_xx, g_y_xy_y_xy, g_y_xy_y_xz, g_y_xy_y_yy, g_y_xy_y_yz, g_y_xy_y_zz, g_y_xy_yzz_xx, g_y_xy_yzz_xy, g_y_xy_yzz_xz, g_y_xy_yzz_yy, g_y_xy_yzz_yz, g_y_xy_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xy_yz_xx[i] = -2.0 * g_y_xy_y_xx[i] * a_exp + 4.0 * g_y_xy_yzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_yz_xy[i] = -2.0 * g_y_xy_y_xy[i] * a_exp + 4.0 * g_y_xy_yzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_yz_xz[i] = -2.0 * g_y_xy_y_xz[i] * a_exp + 4.0 * g_y_xy_yzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_yz_yy[i] = -2.0 * g_y_xy_y_yy[i] * a_exp + 4.0 * g_y_xy_yzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_yz_yz[i] = -2.0 * g_y_xy_y_yz[i] * a_exp + 4.0 * g_y_xy_yzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_yz_zz[i] = -2.0 * g_y_xy_y_zz[i] * a_exp + 4.0 * g_y_xy_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1146-1152)

    #pragma omp simd aligned(g_y_0_z_0_0_xy_zz_xx, g_y_0_z_0_0_xy_zz_xy, g_y_0_z_0_0_xy_zz_xz, g_y_0_z_0_0_xy_zz_yy, g_y_0_z_0_0_xy_zz_yz, g_y_0_z_0_0_xy_zz_zz, g_y_xy_z_xx, g_y_xy_z_xy, g_y_xy_z_xz, g_y_xy_z_yy, g_y_xy_z_yz, g_y_xy_z_zz, g_y_xy_zzz_xx, g_y_xy_zzz_xy, g_y_xy_zzz_xz, g_y_xy_zzz_yy, g_y_xy_zzz_yz, g_y_xy_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xy_zz_xx[i] = -4.0 * g_y_xy_z_xx[i] * a_exp + 4.0 * g_y_xy_zzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_zz_xy[i] = -4.0 * g_y_xy_z_xy[i] * a_exp + 4.0 * g_y_xy_zzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_zz_xz[i] = -4.0 * g_y_xy_z_xz[i] * a_exp + 4.0 * g_y_xy_zzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_zz_yy[i] = -4.0 * g_y_xy_z_yy[i] * a_exp + 4.0 * g_y_xy_zzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_zz_yz[i] = -4.0 * g_y_xy_z_yz[i] * a_exp + 4.0 * g_y_xy_zzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_zz_zz[i] = -4.0 * g_y_xy_z_zz[i] * a_exp + 4.0 * g_y_xy_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1152-1158)

    #pragma omp simd aligned(g_y_0_z_0_0_xz_xx_xx, g_y_0_z_0_0_xz_xx_xy, g_y_0_z_0_0_xz_xx_xz, g_y_0_z_0_0_xz_xx_yy, g_y_0_z_0_0_xz_xx_yz, g_y_0_z_0_0_xz_xx_zz, g_y_xz_xxz_xx, g_y_xz_xxz_xy, g_y_xz_xxz_xz, g_y_xz_xxz_yy, g_y_xz_xxz_yz, g_y_xz_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xz_xx_xx[i] = 4.0 * g_y_xz_xxz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_xx_xy[i] = 4.0 * g_y_xz_xxz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_xx_xz[i] = 4.0 * g_y_xz_xxz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_xx_yy[i] = 4.0 * g_y_xz_xxz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_xx_yz[i] = 4.0 * g_y_xz_xxz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_xx_zz[i] = 4.0 * g_y_xz_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1158-1164)

    #pragma omp simd aligned(g_y_0_z_0_0_xz_xy_xx, g_y_0_z_0_0_xz_xy_xy, g_y_0_z_0_0_xz_xy_xz, g_y_0_z_0_0_xz_xy_yy, g_y_0_z_0_0_xz_xy_yz, g_y_0_z_0_0_xz_xy_zz, g_y_xz_xyz_xx, g_y_xz_xyz_xy, g_y_xz_xyz_xz, g_y_xz_xyz_yy, g_y_xz_xyz_yz, g_y_xz_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xz_xy_xx[i] = 4.0 * g_y_xz_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_xy_xy[i] = 4.0 * g_y_xz_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_xy_xz[i] = 4.0 * g_y_xz_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_xy_yy[i] = 4.0 * g_y_xz_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_xy_yz[i] = 4.0 * g_y_xz_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_xy_zz[i] = 4.0 * g_y_xz_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1164-1170)

    #pragma omp simd aligned(g_y_0_z_0_0_xz_xz_xx, g_y_0_z_0_0_xz_xz_xy, g_y_0_z_0_0_xz_xz_xz, g_y_0_z_0_0_xz_xz_yy, g_y_0_z_0_0_xz_xz_yz, g_y_0_z_0_0_xz_xz_zz, g_y_xz_x_xx, g_y_xz_x_xy, g_y_xz_x_xz, g_y_xz_x_yy, g_y_xz_x_yz, g_y_xz_x_zz, g_y_xz_xzz_xx, g_y_xz_xzz_xy, g_y_xz_xzz_xz, g_y_xz_xzz_yy, g_y_xz_xzz_yz, g_y_xz_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xz_xz_xx[i] = -2.0 * g_y_xz_x_xx[i] * a_exp + 4.0 * g_y_xz_xzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_xz_xy[i] = -2.0 * g_y_xz_x_xy[i] * a_exp + 4.0 * g_y_xz_xzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_xz_xz[i] = -2.0 * g_y_xz_x_xz[i] * a_exp + 4.0 * g_y_xz_xzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_xz_yy[i] = -2.0 * g_y_xz_x_yy[i] * a_exp + 4.0 * g_y_xz_xzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_xz_yz[i] = -2.0 * g_y_xz_x_yz[i] * a_exp + 4.0 * g_y_xz_xzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_xz_zz[i] = -2.0 * g_y_xz_x_zz[i] * a_exp + 4.0 * g_y_xz_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1170-1176)

    #pragma omp simd aligned(g_y_0_z_0_0_xz_yy_xx, g_y_0_z_0_0_xz_yy_xy, g_y_0_z_0_0_xz_yy_xz, g_y_0_z_0_0_xz_yy_yy, g_y_0_z_0_0_xz_yy_yz, g_y_0_z_0_0_xz_yy_zz, g_y_xz_yyz_xx, g_y_xz_yyz_xy, g_y_xz_yyz_xz, g_y_xz_yyz_yy, g_y_xz_yyz_yz, g_y_xz_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xz_yy_xx[i] = 4.0 * g_y_xz_yyz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_yy_xy[i] = 4.0 * g_y_xz_yyz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_yy_xz[i] = 4.0 * g_y_xz_yyz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_yy_yy[i] = 4.0 * g_y_xz_yyz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_yy_yz[i] = 4.0 * g_y_xz_yyz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_yy_zz[i] = 4.0 * g_y_xz_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1176-1182)

    #pragma omp simd aligned(g_y_0_z_0_0_xz_yz_xx, g_y_0_z_0_0_xz_yz_xy, g_y_0_z_0_0_xz_yz_xz, g_y_0_z_0_0_xz_yz_yy, g_y_0_z_0_0_xz_yz_yz, g_y_0_z_0_0_xz_yz_zz, g_y_xz_y_xx, g_y_xz_y_xy, g_y_xz_y_xz, g_y_xz_y_yy, g_y_xz_y_yz, g_y_xz_y_zz, g_y_xz_yzz_xx, g_y_xz_yzz_xy, g_y_xz_yzz_xz, g_y_xz_yzz_yy, g_y_xz_yzz_yz, g_y_xz_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xz_yz_xx[i] = -2.0 * g_y_xz_y_xx[i] * a_exp + 4.0 * g_y_xz_yzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_yz_xy[i] = -2.0 * g_y_xz_y_xy[i] * a_exp + 4.0 * g_y_xz_yzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_yz_xz[i] = -2.0 * g_y_xz_y_xz[i] * a_exp + 4.0 * g_y_xz_yzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_yz_yy[i] = -2.0 * g_y_xz_y_yy[i] * a_exp + 4.0 * g_y_xz_yzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_yz_yz[i] = -2.0 * g_y_xz_y_yz[i] * a_exp + 4.0 * g_y_xz_yzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_yz_zz[i] = -2.0 * g_y_xz_y_zz[i] * a_exp + 4.0 * g_y_xz_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1182-1188)

    #pragma omp simd aligned(g_y_0_z_0_0_xz_zz_xx, g_y_0_z_0_0_xz_zz_xy, g_y_0_z_0_0_xz_zz_xz, g_y_0_z_0_0_xz_zz_yy, g_y_0_z_0_0_xz_zz_yz, g_y_0_z_0_0_xz_zz_zz, g_y_xz_z_xx, g_y_xz_z_xy, g_y_xz_z_xz, g_y_xz_z_yy, g_y_xz_z_yz, g_y_xz_z_zz, g_y_xz_zzz_xx, g_y_xz_zzz_xy, g_y_xz_zzz_xz, g_y_xz_zzz_yy, g_y_xz_zzz_yz, g_y_xz_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xz_zz_xx[i] = -4.0 * g_y_xz_z_xx[i] * a_exp + 4.0 * g_y_xz_zzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_zz_xy[i] = -4.0 * g_y_xz_z_xy[i] * a_exp + 4.0 * g_y_xz_zzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_zz_xz[i] = -4.0 * g_y_xz_z_xz[i] * a_exp + 4.0 * g_y_xz_zzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_zz_yy[i] = -4.0 * g_y_xz_z_yy[i] * a_exp + 4.0 * g_y_xz_zzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_zz_yz[i] = -4.0 * g_y_xz_z_yz[i] * a_exp + 4.0 * g_y_xz_zzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_zz_zz[i] = -4.0 * g_y_xz_z_zz[i] * a_exp + 4.0 * g_y_xz_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1188-1194)

    #pragma omp simd aligned(g_y_0_z_0_0_yy_xx_xx, g_y_0_z_0_0_yy_xx_xy, g_y_0_z_0_0_yy_xx_xz, g_y_0_z_0_0_yy_xx_yy, g_y_0_z_0_0_yy_xx_yz, g_y_0_z_0_0_yy_xx_zz, g_y_yy_xxz_xx, g_y_yy_xxz_xy, g_y_yy_xxz_xz, g_y_yy_xxz_yy, g_y_yy_xxz_yz, g_y_yy_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_yy_xx_xx[i] = 4.0 * g_y_yy_xxz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_xx_xy[i] = 4.0 * g_y_yy_xxz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_xx_xz[i] = 4.0 * g_y_yy_xxz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_xx_yy[i] = 4.0 * g_y_yy_xxz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_xx_yz[i] = 4.0 * g_y_yy_xxz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_xx_zz[i] = 4.0 * g_y_yy_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1194-1200)

    #pragma omp simd aligned(g_y_0_z_0_0_yy_xy_xx, g_y_0_z_0_0_yy_xy_xy, g_y_0_z_0_0_yy_xy_xz, g_y_0_z_0_0_yy_xy_yy, g_y_0_z_0_0_yy_xy_yz, g_y_0_z_0_0_yy_xy_zz, g_y_yy_xyz_xx, g_y_yy_xyz_xy, g_y_yy_xyz_xz, g_y_yy_xyz_yy, g_y_yy_xyz_yz, g_y_yy_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_yy_xy_xx[i] = 4.0 * g_y_yy_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_xy_xy[i] = 4.0 * g_y_yy_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_xy_xz[i] = 4.0 * g_y_yy_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_xy_yy[i] = 4.0 * g_y_yy_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_xy_yz[i] = 4.0 * g_y_yy_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_xy_zz[i] = 4.0 * g_y_yy_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1200-1206)

    #pragma omp simd aligned(g_y_0_z_0_0_yy_xz_xx, g_y_0_z_0_0_yy_xz_xy, g_y_0_z_0_0_yy_xz_xz, g_y_0_z_0_0_yy_xz_yy, g_y_0_z_0_0_yy_xz_yz, g_y_0_z_0_0_yy_xz_zz, g_y_yy_x_xx, g_y_yy_x_xy, g_y_yy_x_xz, g_y_yy_x_yy, g_y_yy_x_yz, g_y_yy_x_zz, g_y_yy_xzz_xx, g_y_yy_xzz_xy, g_y_yy_xzz_xz, g_y_yy_xzz_yy, g_y_yy_xzz_yz, g_y_yy_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_yy_xz_xx[i] = -2.0 * g_y_yy_x_xx[i] * a_exp + 4.0 * g_y_yy_xzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_xz_xy[i] = -2.0 * g_y_yy_x_xy[i] * a_exp + 4.0 * g_y_yy_xzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_xz_xz[i] = -2.0 * g_y_yy_x_xz[i] * a_exp + 4.0 * g_y_yy_xzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_xz_yy[i] = -2.0 * g_y_yy_x_yy[i] * a_exp + 4.0 * g_y_yy_xzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_xz_yz[i] = -2.0 * g_y_yy_x_yz[i] * a_exp + 4.0 * g_y_yy_xzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_xz_zz[i] = -2.0 * g_y_yy_x_zz[i] * a_exp + 4.0 * g_y_yy_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1206-1212)

    #pragma omp simd aligned(g_y_0_z_0_0_yy_yy_xx, g_y_0_z_0_0_yy_yy_xy, g_y_0_z_0_0_yy_yy_xz, g_y_0_z_0_0_yy_yy_yy, g_y_0_z_0_0_yy_yy_yz, g_y_0_z_0_0_yy_yy_zz, g_y_yy_yyz_xx, g_y_yy_yyz_xy, g_y_yy_yyz_xz, g_y_yy_yyz_yy, g_y_yy_yyz_yz, g_y_yy_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_yy_yy_xx[i] = 4.0 * g_y_yy_yyz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_yy_xy[i] = 4.0 * g_y_yy_yyz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_yy_xz[i] = 4.0 * g_y_yy_yyz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_yy_yy[i] = 4.0 * g_y_yy_yyz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_yy_yz[i] = 4.0 * g_y_yy_yyz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_yy_zz[i] = 4.0 * g_y_yy_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1212-1218)

    #pragma omp simd aligned(g_y_0_z_0_0_yy_yz_xx, g_y_0_z_0_0_yy_yz_xy, g_y_0_z_0_0_yy_yz_xz, g_y_0_z_0_0_yy_yz_yy, g_y_0_z_0_0_yy_yz_yz, g_y_0_z_0_0_yy_yz_zz, g_y_yy_y_xx, g_y_yy_y_xy, g_y_yy_y_xz, g_y_yy_y_yy, g_y_yy_y_yz, g_y_yy_y_zz, g_y_yy_yzz_xx, g_y_yy_yzz_xy, g_y_yy_yzz_xz, g_y_yy_yzz_yy, g_y_yy_yzz_yz, g_y_yy_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_yy_yz_xx[i] = -2.0 * g_y_yy_y_xx[i] * a_exp + 4.0 * g_y_yy_yzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_yz_xy[i] = -2.0 * g_y_yy_y_xy[i] * a_exp + 4.0 * g_y_yy_yzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_yz_xz[i] = -2.0 * g_y_yy_y_xz[i] * a_exp + 4.0 * g_y_yy_yzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_yz_yy[i] = -2.0 * g_y_yy_y_yy[i] * a_exp + 4.0 * g_y_yy_yzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_yz_yz[i] = -2.0 * g_y_yy_y_yz[i] * a_exp + 4.0 * g_y_yy_yzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_yz_zz[i] = -2.0 * g_y_yy_y_zz[i] * a_exp + 4.0 * g_y_yy_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1218-1224)

    #pragma omp simd aligned(g_y_0_z_0_0_yy_zz_xx, g_y_0_z_0_0_yy_zz_xy, g_y_0_z_0_0_yy_zz_xz, g_y_0_z_0_0_yy_zz_yy, g_y_0_z_0_0_yy_zz_yz, g_y_0_z_0_0_yy_zz_zz, g_y_yy_z_xx, g_y_yy_z_xy, g_y_yy_z_xz, g_y_yy_z_yy, g_y_yy_z_yz, g_y_yy_z_zz, g_y_yy_zzz_xx, g_y_yy_zzz_xy, g_y_yy_zzz_xz, g_y_yy_zzz_yy, g_y_yy_zzz_yz, g_y_yy_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_yy_zz_xx[i] = -4.0 * g_y_yy_z_xx[i] * a_exp + 4.0 * g_y_yy_zzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_zz_xy[i] = -4.0 * g_y_yy_z_xy[i] * a_exp + 4.0 * g_y_yy_zzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_zz_xz[i] = -4.0 * g_y_yy_z_xz[i] * a_exp + 4.0 * g_y_yy_zzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_zz_yy[i] = -4.0 * g_y_yy_z_yy[i] * a_exp + 4.0 * g_y_yy_zzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_zz_yz[i] = -4.0 * g_y_yy_z_yz[i] * a_exp + 4.0 * g_y_yy_zzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_zz_zz[i] = -4.0 * g_y_yy_z_zz[i] * a_exp + 4.0 * g_y_yy_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1224-1230)

    #pragma omp simd aligned(g_y_0_z_0_0_yz_xx_xx, g_y_0_z_0_0_yz_xx_xy, g_y_0_z_0_0_yz_xx_xz, g_y_0_z_0_0_yz_xx_yy, g_y_0_z_0_0_yz_xx_yz, g_y_0_z_0_0_yz_xx_zz, g_y_yz_xxz_xx, g_y_yz_xxz_xy, g_y_yz_xxz_xz, g_y_yz_xxz_yy, g_y_yz_xxz_yz, g_y_yz_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_yz_xx_xx[i] = 4.0 * g_y_yz_xxz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_xx_xy[i] = 4.0 * g_y_yz_xxz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_xx_xz[i] = 4.0 * g_y_yz_xxz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_xx_yy[i] = 4.0 * g_y_yz_xxz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_xx_yz[i] = 4.0 * g_y_yz_xxz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_xx_zz[i] = 4.0 * g_y_yz_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1230-1236)

    #pragma omp simd aligned(g_y_0_z_0_0_yz_xy_xx, g_y_0_z_0_0_yz_xy_xy, g_y_0_z_0_0_yz_xy_xz, g_y_0_z_0_0_yz_xy_yy, g_y_0_z_0_0_yz_xy_yz, g_y_0_z_0_0_yz_xy_zz, g_y_yz_xyz_xx, g_y_yz_xyz_xy, g_y_yz_xyz_xz, g_y_yz_xyz_yy, g_y_yz_xyz_yz, g_y_yz_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_yz_xy_xx[i] = 4.0 * g_y_yz_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_xy_xy[i] = 4.0 * g_y_yz_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_xy_xz[i] = 4.0 * g_y_yz_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_xy_yy[i] = 4.0 * g_y_yz_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_xy_yz[i] = 4.0 * g_y_yz_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_xy_zz[i] = 4.0 * g_y_yz_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1236-1242)

    #pragma omp simd aligned(g_y_0_z_0_0_yz_xz_xx, g_y_0_z_0_0_yz_xz_xy, g_y_0_z_0_0_yz_xz_xz, g_y_0_z_0_0_yz_xz_yy, g_y_0_z_0_0_yz_xz_yz, g_y_0_z_0_0_yz_xz_zz, g_y_yz_x_xx, g_y_yz_x_xy, g_y_yz_x_xz, g_y_yz_x_yy, g_y_yz_x_yz, g_y_yz_x_zz, g_y_yz_xzz_xx, g_y_yz_xzz_xy, g_y_yz_xzz_xz, g_y_yz_xzz_yy, g_y_yz_xzz_yz, g_y_yz_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_yz_xz_xx[i] = -2.0 * g_y_yz_x_xx[i] * a_exp + 4.0 * g_y_yz_xzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_xz_xy[i] = -2.0 * g_y_yz_x_xy[i] * a_exp + 4.0 * g_y_yz_xzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_xz_xz[i] = -2.0 * g_y_yz_x_xz[i] * a_exp + 4.0 * g_y_yz_xzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_xz_yy[i] = -2.0 * g_y_yz_x_yy[i] * a_exp + 4.0 * g_y_yz_xzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_xz_yz[i] = -2.0 * g_y_yz_x_yz[i] * a_exp + 4.0 * g_y_yz_xzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_xz_zz[i] = -2.0 * g_y_yz_x_zz[i] * a_exp + 4.0 * g_y_yz_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1242-1248)

    #pragma omp simd aligned(g_y_0_z_0_0_yz_yy_xx, g_y_0_z_0_0_yz_yy_xy, g_y_0_z_0_0_yz_yy_xz, g_y_0_z_0_0_yz_yy_yy, g_y_0_z_0_0_yz_yy_yz, g_y_0_z_0_0_yz_yy_zz, g_y_yz_yyz_xx, g_y_yz_yyz_xy, g_y_yz_yyz_xz, g_y_yz_yyz_yy, g_y_yz_yyz_yz, g_y_yz_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_yz_yy_xx[i] = 4.0 * g_y_yz_yyz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_yy_xy[i] = 4.0 * g_y_yz_yyz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_yy_xz[i] = 4.0 * g_y_yz_yyz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_yy_yy[i] = 4.0 * g_y_yz_yyz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_yy_yz[i] = 4.0 * g_y_yz_yyz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_yy_zz[i] = 4.0 * g_y_yz_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1248-1254)

    #pragma omp simd aligned(g_y_0_z_0_0_yz_yz_xx, g_y_0_z_0_0_yz_yz_xy, g_y_0_z_0_0_yz_yz_xz, g_y_0_z_0_0_yz_yz_yy, g_y_0_z_0_0_yz_yz_yz, g_y_0_z_0_0_yz_yz_zz, g_y_yz_y_xx, g_y_yz_y_xy, g_y_yz_y_xz, g_y_yz_y_yy, g_y_yz_y_yz, g_y_yz_y_zz, g_y_yz_yzz_xx, g_y_yz_yzz_xy, g_y_yz_yzz_xz, g_y_yz_yzz_yy, g_y_yz_yzz_yz, g_y_yz_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_yz_yz_xx[i] = -2.0 * g_y_yz_y_xx[i] * a_exp + 4.0 * g_y_yz_yzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_yz_xy[i] = -2.0 * g_y_yz_y_xy[i] * a_exp + 4.0 * g_y_yz_yzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_yz_xz[i] = -2.0 * g_y_yz_y_xz[i] * a_exp + 4.0 * g_y_yz_yzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_yz_yy[i] = -2.0 * g_y_yz_y_yy[i] * a_exp + 4.0 * g_y_yz_yzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_yz_yz[i] = -2.0 * g_y_yz_y_yz[i] * a_exp + 4.0 * g_y_yz_yzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_yz_zz[i] = -2.0 * g_y_yz_y_zz[i] * a_exp + 4.0 * g_y_yz_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1254-1260)

    #pragma omp simd aligned(g_y_0_z_0_0_yz_zz_xx, g_y_0_z_0_0_yz_zz_xy, g_y_0_z_0_0_yz_zz_xz, g_y_0_z_0_0_yz_zz_yy, g_y_0_z_0_0_yz_zz_yz, g_y_0_z_0_0_yz_zz_zz, g_y_yz_z_xx, g_y_yz_z_xy, g_y_yz_z_xz, g_y_yz_z_yy, g_y_yz_z_yz, g_y_yz_z_zz, g_y_yz_zzz_xx, g_y_yz_zzz_xy, g_y_yz_zzz_xz, g_y_yz_zzz_yy, g_y_yz_zzz_yz, g_y_yz_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_yz_zz_xx[i] = -4.0 * g_y_yz_z_xx[i] * a_exp + 4.0 * g_y_yz_zzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_zz_xy[i] = -4.0 * g_y_yz_z_xy[i] * a_exp + 4.0 * g_y_yz_zzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_zz_xz[i] = -4.0 * g_y_yz_z_xz[i] * a_exp + 4.0 * g_y_yz_zzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_zz_yy[i] = -4.0 * g_y_yz_z_yy[i] * a_exp + 4.0 * g_y_yz_zzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_zz_yz[i] = -4.0 * g_y_yz_z_yz[i] * a_exp + 4.0 * g_y_yz_zzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_zz_zz[i] = -4.0 * g_y_yz_z_zz[i] * a_exp + 4.0 * g_y_yz_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1260-1266)

    #pragma omp simd aligned(g_y_0_z_0_0_zz_xx_xx, g_y_0_z_0_0_zz_xx_xy, g_y_0_z_0_0_zz_xx_xz, g_y_0_z_0_0_zz_xx_yy, g_y_0_z_0_0_zz_xx_yz, g_y_0_z_0_0_zz_xx_zz, g_y_zz_xxz_xx, g_y_zz_xxz_xy, g_y_zz_xxz_xz, g_y_zz_xxz_yy, g_y_zz_xxz_yz, g_y_zz_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_zz_xx_xx[i] = 4.0 * g_y_zz_xxz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_xx_xy[i] = 4.0 * g_y_zz_xxz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_xx_xz[i] = 4.0 * g_y_zz_xxz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_xx_yy[i] = 4.0 * g_y_zz_xxz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_xx_yz[i] = 4.0 * g_y_zz_xxz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_xx_zz[i] = 4.0 * g_y_zz_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1266-1272)

    #pragma omp simd aligned(g_y_0_z_0_0_zz_xy_xx, g_y_0_z_0_0_zz_xy_xy, g_y_0_z_0_0_zz_xy_xz, g_y_0_z_0_0_zz_xy_yy, g_y_0_z_0_0_zz_xy_yz, g_y_0_z_0_0_zz_xy_zz, g_y_zz_xyz_xx, g_y_zz_xyz_xy, g_y_zz_xyz_xz, g_y_zz_xyz_yy, g_y_zz_xyz_yz, g_y_zz_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_zz_xy_xx[i] = 4.0 * g_y_zz_xyz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_xy_xy[i] = 4.0 * g_y_zz_xyz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_xy_xz[i] = 4.0 * g_y_zz_xyz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_xy_yy[i] = 4.0 * g_y_zz_xyz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_xy_yz[i] = 4.0 * g_y_zz_xyz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_xy_zz[i] = 4.0 * g_y_zz_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1272-1278)

    #pragma omp simd aligned(g_y_0_z_0_0_zz_xz_xx, g_y_0_z_0_0_zz_xz_xy, g_y_0_z_0_0_zz_xz_xz, g_y_0_z_0_0_zz_xz_yy, g_y_0_z_0_0_zz_xz_yz, g_y_0_z_0_0_zz_xz_zz, g_y_zz_x_xx, g_y_zz_x_xy, g_y_zz_x_xz, g_y_zz_x_yy, g_y_zz_x_yz, g_y_zz_x_zz, g_y_zz_xzz_xx, g_y_zz_xzz_xy, g_y_zz_xzz_xz, g_y_zz_xzz_yy, g_y_zz_xzz_yz, g_y_zz_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_zz_xz_xx[i] = -2.0 * g_y_zz_x_xx[i] * a_exp + 4.0 * g_y_zz_xzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_xz_xy[i] = -2.0 * g_y_zz_x_xy[i] * a_exp + 4.0 * g_y_zz_xzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_xz_xz[i] = -2.0 * g_y_zz_x_xz[i] * a_exp + 4.0 * g_y_zz_xzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_xz_yy[i] = -2.0 * g_y_zz_x_yy[i] * a_exp + 4.0 * g_y_zz_xzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_xz_yz[i] = -2.0 * g_y_zz_x_yz[i] * a_exp + 4.0 * g_y_zz_xzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_xz_zz[i] = -2.0 * g_y_zz_x_zz[i] * a_exp + 4.0 * g_y_zz_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1278-1284)

    #pragma omp simd aligned(g_y_0_z_0_0_zz_yy_xx, g_y_0_z_0_0_zz_yy_xy, g_y_0_z_0_0_zz_yy_xz, g_y_0_z_0_0_zz_yy_yy, g_y_0_z_0_0_zz_yy_yz, g_y_0_z_0_0_zz_yy_zz, g_y_zz_yyz_xx, g_y_zz_yyz_xy, g_y_zz_yyz_xz, g_y_zz_yyz_yy, g_y_zz_yyz_yz, g_y_zz_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_zz_yy_xx[i] = 4.0 * g_y_zz_yyz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_yy_xy[i] = 4.0 * g_y_zz_yyz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_yy_xz[i] = 4.0 * g_y_zz_yyz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_yy_yy[i] = 4.0 * g_y_zz_yyz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_yy_yz[i] = 4.0 * g_y_zz_yyz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_yy_zz[i] = 4.0 * g_y_zz_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1284-1290)

    #pragma omp simd aligned(g_y_0_z_0_0_zz_yz_xx, g_y_0_z_0_0_zz_yz_xy, g_y_0_z_0_0_zz_yz_xz, g_y_0_z_0_0_zz_yz_yy, g_y_0_z_0_0_zz_yz_yz, g_y_0_z_0_0_zz_yz_zz, g_y_zz_y_xx, g_y_zz_y_xy, g_y_zz_y_xz, g_y_zz_y_yy, g_y_zz_y_yz, g_y_zz_y_zz, g_y_zz_yzz_xx, g_y_zz_yzz_xy, g_y_zz_yzz_xz, g_y_zz_yzz_yy, g_y_zz_yzz_yz, g_y_zz_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_zz_yz_xx[i] = -2.0 * g_y_zz_y_xx[i] * a_exp + 4.0 * g_y_zz_yzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_yz_xy[i] = -2.0 * g_y_zz_y_xy[i] * a_exp + 4.0 * g_y_zz_yzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_yz_xz[i] = -2.0 * g_y_zz_y_xz[i] * a_exp + 4.0 * g_y_zz_yzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_yz_yy[i] = -2.0 * g_y_zz_y_yy[i] * a_exp + 4.0 * g_y_zz_yzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_yz_yz[i] = -2.0 * g_y_zz_y_yz[i] * a_exp + 4.0 * g_y_zz_yzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_yz_zz[i] = -2.0 * g_y_zz_y_zz[i] * a_exp + 4.0 * g_y_zz_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1290-1296)

    #pragma omp simd aligned(g_y_0_z_0_0_zz_zz_xx, g_y_0_z_0_0_zz_zz_xy, g_y_0_z_0_0_zz_zz_xz, g_y_0_z_0_0_zz_zz_yy, g_y_0_z_0_0_zz_zz_yz, g_y_0_z_0_0_zz_zz_zz, g_y_zz_z_xx, g_y_zz_z_xy, g_y_zz_z_xz, g_y_zz_z_yy, g_y_zz_z_yz, g_y_zz_z_zz, g_y_zz_zzz_xx, g_y_zz_zzz_xy, g_y_zz_zzz_xz, g_y_zz_zzz_yy, g_y_zz_zzz_yz, g_y_zz_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_zz_zz_xx[i] = -4.0 * g_y_zz_z_xx[i] * a_exp + 4.0 * g_y_zz_zzz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_zz_xy[i] = -4.0 * g_y_zz_z_xy[i] * a_exp + 4.0 * g_y_zz_zzz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_zz_xz[i] = -4.0 * g_y_zz_z_xz[i] * a_exp + 4.0 * g_y_zz_zzz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_zz_yy[i] = -4.0 * g_y_zz_z_yy[i] * a_exp + 4.0 * g_y_zz_zzz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_zz_yz[i] = -4.0 * g_y_zz_z_yz[i] * a_exp + 4.0 * g_y_zz_zzz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_zz_zz[i] = -4.0 * g_y_zz_z_zz[i] * a_exp + 4.0 * g_y_zz_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1296-1302)

    #pragma omp simd aligned(g_z_0_x_0_0_xx_xx_xx, g_z_0_x_0_0_xx_xx_xy, g_z_0_x_0_0_xx_xx_xz, g_z_0_x_0_0_xx_xx_yy, g_z_0_x_0_0_xx_xx_yz, g_z_0_x_0_0_xx_xx_zz, g_z_xx_x_xx, g_z_xx_x_xy, g_z_xx_x_xz, g_z_xx_x_yy, g_z_xx_x_yz, g_z_xx_x_zz, g_z_xx_xxx_xx, g_z_xx_xxx_xy, g_z_xx_xxx_xz, g_z_xx_xxx_yy, g_z_xx_xxx_yz, g_z_xx_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xx_xx_xx[i] = -4.0 * g_z_xx_x_xx[i] * a_exp + 4.0 * g_z_xx_xxx_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_xx_xy[i] = -4.0 * g_z_xx_x_xy[i] * a_exp + 4.0 * g_z_xx_xxx_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_xx_xz[i] = -4.0 * g_z_xx_x_xz[i] * a_exp + 4.0 * g_z_xx_xxx_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_xx_yy[i] = -4.0 * g_z_xx_x_yy[i] * a_exp + 4.0 * g_z_xx_xxx_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_xx_yz[i] = -4.0 * g_z_xx_x_yz[i] * a_exp + 4.0 * g_z_xx_xxx_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_xx_zz[i] = -4.0 * g_z_xx_x_zz[i] * a_exp + 4.0 * g_z_xx_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1302-1308)

    #pragma omp simd aligned(g_z_0_x_0_0_xx_xy_xx, g_z_0_x_0_0_xx_xy_xy, g_z_0_x_0_0_xx_xy_xz, g_z_0_x_0_0_xx_xy_yy, g_z_0_x_0_0_xx_xy_yz, g_z_0_x_0_0_xx_xy_zz, g_z_xx_xxy_xx, g_z_xx_xxy_xy, g_z_xx_xxy_xz, g_z_xx_xxy_yy, g_z_xx_xxy_yz, g_z_xx_xxy_zz, g_z_xx_y_xx, g_z_xx_y_xy, g_z_xx_y_xz, g_z_xx_y_yy, g_z_xx_y_yz, g_z_xx_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xx_xy_xx[i] = -2.0 * g_z_xx_y_xx[i] * a_exp + 4.0 * g_z_xx_xxy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_xy_xy[i] = -2.0 * g_z_xx_y_xy[i] * a_exp + 4.0 * g_z_xx_xxy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_xy_xz[i] = -2.0 * g_z_xx_y_xz[i] * a_exp + 4.0 * g_z_xx_xxy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_xy_yy[i] = -2.0 * g_z_xx_y_yy[i] * a_exp + 4.0 * g_z_xx_xxy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_xy_yz[i] = -2.0 * g_z_xx_y_yz[i] * a_exp + 4.0 * g_z_xx_xxy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_xy_zz[i] = -2.0 * g_z_xx_y_zz[i] * a_exp + 4.0 * g_z_xx_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1308-1314)

    #pragma omp simd aligned(g_z_0_x_0_0_xx_xz_xx, g_z_0_x_0_0_xx_xz_xy, g_z_0_x_0_0_xx_xz_xz, g_z_0_x_0_0_xx_xz_yy, g_z_0_x_0_0_xx_xz_yz, g_z_0_x_0_0_xx_xz_zz, g_z_xx_xxz_xx, g_z_xx_xxz_xy, g_z_xx_xxz_xz, g_z_xx_xxz_yy, g_z_xx_xxz_yz, g_z_xx_xxz_zz, g_z_xx_z_xx, g_z_xx_z_xy, g_z_xx_z_xz, g_z_xx_z_yy, g_z_xx_z_yz, g_z_xx_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xx_xz_xx[i] = -2.0 * g_z_xx_z_xx[i] * a_exp + 4.0 * g_z_xx_xxz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_xz_xy[i] = -2.0 * g_z_xx_z_xy[i] * a_exp + 4.0 * g_z_xx_xxz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_xz_xz[i] = -2.0 * g_z_xx_z_xz[i] * a_exp + 4.0 * g_z_xx_xxz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_xz_yy[i] = -2.0 * g_z_xx_z_yy[i] * a_exp + 4.0 * g_z_xx_xxz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_xz_yz[i] = -2.0 * g_z_xx_z_yz[i] * a_exp + 4.0 * g_z_xx_xxz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_xz_zz[i] = -2.0 * g_z_xx_z_zz[i] * a_exp + 4.0 * g_z_xx_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1314-1320)

    #pragma omp simd aligned(g_z_0_x_0_0_xx_yy_xx, g_z_0_x_0_0_xx_yy_xy, g_z_0_x_0_0_xx_yy_xz, g_z_0_x_0_0_xx_yy_yy, g_z_0_x_0_0_xx_yy_yz, g_z_0_x_0_0_xx_yy_zz, g_z_xx_xyy_xx, g_z_xx_xyy_xy, g_z_xx_xyy_xz, g_z_xx_xyy_yy, g_z_xx_xyy_yz, g_z_xx_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xx_yy_xx[i] = 4.0 * g_z_xx_xyy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_yy_xy[i] = 4.0 * g_z_xx_xyy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_yy_xz[i] = 4.0 * g_z_xx_xyy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_yy_yy[i] = 4.0 * g_z_xx_xyy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_yy_yz[i] = 4.0 * g_z_xx_xyy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_yy_zz[i] = 4.0 * g_z_xx_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1320-1326)

    #pragma omp simd aligned(g_z_0_x_0_0_xx_yz_xx, g_z_0_x_0_0_xx_yz_xy, g_z_0_x_0_0_xx_yz_xz, g_z_0_x_0_0_xx_yz_yy, g_z_0_x_0_0_xx_yz_yz, g_z_0_x_0_0_xx_yz_zz, g_z_xx_xyz_xx, g_z_xx_xyz_xy, g_z_xx_xyz_xz, g_z_xx_xyz_yy, g_z_xx_xyz_yz, g_z_xx_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xx_yz_xx[i] = 4.0 * g_z_xx_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_yz_xy[i] = 4.0 * g_z_xx_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_yz_xz[i] = 4.0 * g_z_xx_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_yz_yy[i] = 4.0 * g_z_xx_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_yz_yz[i] = 4.0 * g_z_xx_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_yz_zz[i] = 4.0 * g_z_xx_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1326-1332)

    #pragma omp simd aligned(g_z_0_x_0_0_xx_zz_xx, g_z_0_x_0_0_xx_zz_xy, g_z_0_x_0_0_xx_zz_xz, g_z_0_x_0_0_xx_zz_yy, g_z_0_x_0_0_xx_zz_yz, g_z_0_x_0_0_xx_zz_zz, g_z_xx_xzz_xx, g_z_xx_xzz_xy, g_z_xx_xzz_xz, g_z_xx_xzz_yy, g_z_xx_xzz_yz, g_z_xx_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xx_zz_xx[i] = 4.0 * g_z_xx_xzz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_zz_xy[i] = 4.0 * g_z_xx_xzz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_zz_xz[i] = 4.0 * g_z_xx_xzz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_zz_yy[i] = 4.0 * g_z_xx_xzz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_zz_yz[i] = 4.0 * g_z_xx_xzz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_zz_zz[i] = 4.0 * g_z_xx_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1332-1338)

    #pragma omp simd aligned(g_z_0_x_0_0_xy_xx_xx, g_z_0_x_0_0_xy_xx_xy, g_z_0_x_0_0_xy_xx_xz, g_z_0_x_0_0_xy_xx_yy, g_z_0_x_0_0_xy_xx_yz, g_z_0_x_0_0_xy_xx_zz, g_z_xy_x_xx, g_z_xy_x_xy, g_z_xy_x_xz, g_z_xy_x_yy, g_z_xy_x_yz, g_z_xy_x_zz, g_z_xy_xxx_xx, g_z_xy_xxx_xy, g_z_xy_xxx_xz, g_z_xy_xxx_yy, g_z_xy_xxx_yz, g_z_xy_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xy_xx_xx[i] = -4.0 * g_z_xy_x_xx[i] * a_exp + 4.0 * g_z_xy_xxx_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_xx_xy[i] = -4.0 * g_z_xy_x_xy[i] * a_exp + 4.0 * g_z_xy_xxx_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_xx_xz[i] = -4.0 * g_z_xy_x_xz[i] * a_exp + 4.0 * g_z_xy_xxx_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_xx_yy[i] = -4.0 * g_z_xy_x_yy[i] * a_exp + 4.0 * g_z_xy_xxx_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_xx_yz[i] = -4.0 * g_z_xy_x_yz[i] * a_exp + 4.0 * g_z_xy_xxx_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_xx_zz[i] = -4.0 * g_z_xy_x_zz[i] * a_exp + 4.0 * g_z_xy_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1338-1344)

    #pragma omp simd aligned(g_z_0_x_0_0_xy_xy_xx, g_z_0_x_0_0_xy_xy_xy, g_z_0_x_0_0_xy_xy_xz, g_z_0_x_0_0_xy_xy_yy, g_z_0_x_0_0_xy_xy_yz, g_z_0_x_0_0_xy_xy_zz, g_z_xy_xxy_xx, g_z_xy_xxy_xy, g_z_xy_xxy_xz, g_z_xy_xxy_yy, g_z_xy_xxy_yz, g_z_xy_xxy_zz, g_z_xy_y_xx, g_z_xy_y_xy, g_z_xy_y_xz, g_z_xy_y_yy, g_z_xy_y_yz, g_z_xy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xy_xy_xx[i] = -2.0 * g_z_xy_y_xx[i] * a_exp + 4.0 * g_z_xy_xxy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_xy_xy[i] = -2.0 * g_z_xy_y_xy[i] * a_exp + 4.0 * g_z_xy_xxy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_xy_xz[i] = -2.0 * g_z_xy_y_xz[i] * a_exp + 4.0 * g_z_xy_xxy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_xy_yy[i] = -2.0 * g_z_xy_y_yy[i] * a_exp + 4.0 * g_z_xy_xxy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_xy_yz[i] = -2.0 * g_z_xy_y_yz[i] * a_exp + 4.0 * g_z_xy_xxy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_xy_zz[i] = -2.0 * g_z_xy_y_zz[i] * a_exp + 4.0 * g_z_xy_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1344-1350)

    #pragma omp simd aligned(g_z_0_x_0_0_xy_xz_xx, g_z_0_x_0_0_xy_xz_xy, g_z_0_x_0_0_xy_xz_xz, g_z_0_x_0_0_xy_xz_yy, g_z_0_x_0_0_xy_xz_yz, g_z_0_x_0_0_xy_xz_zz, g_z_xy_xxz_xx, g_z_xy_xxz_xy, g_z_xy_xxz_xz, g_z_xy_xxz_yy, g_z_xy_xxz_yz, g_z_xy_xxz_zz, g_z_xy_z_xx, g_z_xy_z_xy, g_z_xy_z_xz, g_z_xy_z_yy, g_z_xy_z_yz, g_z_xy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xy_xz_xx[i] = -2.0 * g_z_xy_z_xx[i] * a_exp + 4.0 * g_z_xy_xxz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_xz_xy[i] = -2.0 * g_z_xy_z_xy[i] * a_exp + 4.0 * g_z_xy_xxz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_xz_xz[i] = -2.0 * g_z_xy_z_xz[i] * a_exp + 4.0 * g_z_xy_xxz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_xz_yy[i] = -2.0 * g_z_xy_z_yy[i] * a_exp + 4.0 * g_z_xy_xxz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_xz_yz[i] = -2.0 * g_z_xy_z_yz[i] * a_exp + 4.0 * g_z_xy_xxz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_xz_zz[i] = -2.0 * g_z_xy_z_zz[i] * a_exp + 4.0 * g_z_xy_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1350-1356)

    #pragma omp simd aligned(g_z_0_x_0_0_xy_yy_xx, g_z_0_x_0_0_xy_yy_xy, g_z_0_x_0_0_xy_yy_xz, g_z_0_x_0_0_xy_yy_yy, g_z_0_x_0_0_xy_yy_yz, g_z_0_x_0_0_xy_yy_zz, g_z_xy_xyy_xx, g_z_xy_xyy_xy, g_z_xy_xyy_xz, g_z_xy_xyy_yy, g_z_xy_xyy_yz, g_z_xy_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xy_yy_xx[i] = 4.0 * g_z_xy_xyy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_yy_xy[i] = 4.0 * g_z_xy_xyy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_yy_xz[i] = 4.0 * g_z_xy_xyy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_yy_yy[i] = 4.0 * g_z_xy_xyy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_yy_yz[i] = 4.0 * g_z_xy_xyy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_yy_zz[i] = 4.0 * g_z_xy_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1356-1362)

    #pragma omp simd aligned(g_z_0_x_0_0_xy_yz_xx, g_z_0_x_0_0_xy_yz_xy, g_z_0_x_0_0_xy_yz_xz, g_z_0_x_0_0_xy_yz_yy, g_z_0_x_0_0_xy_yz_yz, g_z_0_x_0_0_xy_yz_zz, g_z_xy_xyz_xx, g_z_xy_xyz_xy, g_z_xy_xyz_xz, g_z_xy_xyz_yy, g_z_xy_xyz_yz, g_z_xy_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xy_yz_xx[i] = 4.0 * g_z_xy_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_yz_xy[i] = 4.0 * g_z_xy_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_yz_xz[i] = 4.0 * g_z_xy_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_yz_yy[i] = 4.0 * g_z_xy_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_yz_yz[i] = 4.0 * g_z_xy_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_yz_zz[i] = 4.0 * g_z_xy_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1362-1368)

    #pragma omp simd aligned(g_z_0_x_0_0_xy_zz_xx, g_z_0_x_0_0_xy_zz_xy, g_z_0_x_0_0_xy_zz_xz, g_z_0_x_0_0_xy_zz_yy, g_z_0_x_0_0_xy_zz_yz, g_z_0_x_0_0_xy_zz_zz, g_z_xy_xzz_xx, g_z_xy_xzz_xy, g_z_xy_xzz_xz, g_z_xy_xzz_yy, g_z_xy_xzz_yz, g_z_xy_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xy_zz_xx[i] = 4.0 * g_z_xy_xzz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_zz_xy[i] = 4.0 * g_z_xy_xzz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_zz_xz[i] = 4.0 * g_z_xy_xzz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_zz_yy[i] = 4.0 * g_z_xy_xzz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_zz_yz[i] = 4.0 * g_z_xy_xzz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_zz_zz[i] = 4.0 * g_z_xy_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1368-1374)

    #pragma omp simd aligned(g_z_0_x_0_0_xz_xx_xx, g_z_0_x_0_0_xz_xx_xy, g_z_0_x_0_0_xz_xx_xz, g_z_0_x_0_0_xz_xx_yy, g_z_0_x_0_0_xz_xx_yz, g_z_0_x_0_0_xz_xx_zz, g_z_xz_x_xx, g_z_xz_x_xy, g_z_xz_x_xz, g_z_xz_x_yy, g_z_xz_x_yz, g_z_xz_x_zz, g_z_xz_xxx_xx, g_z_xz_xxx_xy, g_z_xz_xxx_xz, g_z_xz_xxx_yy, g_z_xz_xxx_yz, g_z_xz_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xz_xx_xx[i] = -4.0 * g_z_xz_x_xx[i] * a_exp + 4.0 * g_z_xz_xxx_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_xx_xy[i] = -4.0 * g_z_xz_x_xy[i] * a_exp + 4.0 * g_z_xz_xxx_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_xx_xz[i] = -4.0 * g_z_xz_x_xz[i] * a_exp + 4.0 * g_z_xz_xxx_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_xx_yy[i] = -4.0 * g_z_xz_x_yy[i] * a_exp + 4.0 * g_z_xz_xxx_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_xx_yz[i] = -4.0 * g_z_xz_x_yz[i] * a_exp + 4.0 * g_z_xz_xxx_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_xx_zz[i] = -4.0 * g_z_xz_x_zz[i] * a_exp + 4.0 * g_z_xz_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1374-1380)

    #pragma omp simd aligned(g_z_0_x_0_0_xz_xy_xx, g_z_0_x_0_0_xz_xy_xy, g_z_0_x_0_0_xz_xy_xz, g_z_0_x_0_0_xz_xy_yy, g_z_0_x_0_0_xz_xy_yz, g_z_0_x_0_0_xz_xy_zz, g_z_xz_xxy_xx, g_z_xz_xxy_xy, g_z_xz_xxy_xz, g_z_xz_xxy_yy, g_z_xz_xxy_yz, g_z_xz_xxy_zz, g_z_xz_y_xx, g_z_xz_y_xy, g_z_xz_y_xz, g_z_xz_y_yy, g_z_xz_y_yz, g_z_xz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xz_xy_xx[i] = -2.0 * g_z_xz_y_xx[i] * a_exp + 4.0 * g_z_xz_xxy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_xy_xy[i] = -2.0 * g_z_xz_y_xy[i] * a_exp + 4.0 * g_z_xz_xxy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_xy_xz[i] = -2.0 * g_z_xz_y_xz[i] * a_exp + 4.0 * g_z_xz_xxy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_xy_yy[i] = -2.0 * g_z_xz_y_yy[i] * a_exp + 4.0 * g_z_xz_xxy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_xy_yz[i] = -2.0 * g_z_xz_y_yz[i] * a_exp + 4.0 * g_z_xz_xxy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_xy_zz[i] = -2.0 * g_z_xz_y_zz[i] * a_exp + 4.0 * g_z_xz_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1380-1386)

    #pragma omp simd aligned(g_z_0_x_0_0_xz_xz_xx, g_z_0_x_0_0_xz_xz_xy, g_z_0_x_0_0_xz_xz_xz, g_z_0_x_0_0_xz_xz_yy, g_z_0_x_0_0_xz_xz_yz, g_z_0_x_0_0_xz_xz_zz, g_z_xz_xxz_xx, g_z_xz_xxz_xy, g_z_xz_xxz_xz, g_z_xz_xxz_yy, g_z_xz_xxz_yz, g_z_xz_xxz_zz, g_z_xz_z_xx, g_z_xz_z_xy, g_z_xz_z_xz, g_z_xz_z_yy, g_z_xz_z_yz, g_z_xz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xz_xz_xx[i] = -2.0 * g_z_xz_z_xx[i] * a_exp + 4.0 * g_z_xz_xxz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_xz_xy[i] = -2.0 * g_z_xz_z_xy[i] * a_exp + 4.0 * g_z_xz_xxz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_xz_xz[i] = -2.0 * g_z_xz_z_xz[i] * a_exp + 4.0 * g_z_xz_xxz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_xz_yy[i] = -2.0 * g_z_xz_z_yy[i] * a_exp + 4.0 * g_z_xz_xxz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_xz_yz[i] = -2.0 * g_z_xz_z_yz[i] * a_exp + 4.0 * g_z_xz_xxz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_xz_zz[i] = -2.0 * g_z_xz_z_zz[i] * a_exp + 4.0 * g_z_xz_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1386-1392)

    #pragma omp simd aligned(g_z_0_x_0_0_xz_yy_xx, g_z_0_x_0_0_xz_yy_xy, g_z_0_x_0_0_xz_yy_xz, g_z_0_x_0_0_xz_yy_yy, g_z_0_x_0_0_xz_yy_yz, g_z_0_x_0_0_xz_yy_zz, g_z_xz_xyy_xx, g_z_xz_xyy_xy, g_z_xz_xyy_xz, g_z_xz_xyy_yy, g_z_xz_xyy_yz, g_z_xz_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xz_yy_xx[i] = 4.0 * g_z_xz_xyy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_yy_xy[i] = 4.0 * g_z_xz_xyy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_yy_xz[i] = 4.0 * g_z_xz_xyy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_yy_yy[i] = 4.0 * g_z_xz_xyy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_yy_yz[i] = 4.0 * g_z_xz_xyy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_yy_zz[i] = 4.0 * g_z_xz_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1392-1398)

    #pragma omp simd aligned(g_z_0_x_0_0_xz_yz_xx, g_z_0_x_0_0_xz_yz_xy, g_z_0_x_0_0_xz_yz_xz, g_z_0_x_0_0_xz_yz_yy, g_z_0_x_0_0_xz_yz_yz, g_z_0_x_0_0_xz_yz_zz, g_z_xz_xyz_xx, g_z_xz_xyz_xy, g_z_xz_xyz_xz, g_z_xz_xyz_yy, g_z_xz_xyz_yz, g_z_xz_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xz_yz_xx[i] = 4.0 * g_z_xz_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_yz_xy[i] = 4.0 * g_z_xz_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_yz_xz[i] = 4.0 * g_z_xz_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_yz_yy[i] = 4.0 * g_z_xz_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_yz_yz[i] = 4.0 * g_z_xz_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_yz_zz[i] = 4.0 * g_z_xz_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1398-1404)

    #pragma omp simd aligned(g_z_0_x_0_0_xz_zz_xx, g_z_0_x_0_0_xz_zz_xy, g_z_0_x_0_0_xz_zz_xz, g_z_0_x_0_0_xz_zz_yy, g_z_0_x_0_0_xz_zz_yz, g_z_0_x_0_0_xz_zz_zz, g_z_xz_xzz_xx, g_z_xz_xzz_xy, g_z_xz_xzz_xz, g_z_xz_xzz_yy, g_z_xz_xzz_yz, g_z_xz_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xz_zz_xx[i] = 4.0 * g_z_xz_xzz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_zz_xy[i] = 4.0 * g_z_xz_xzz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_zz_xz[i] = 4.0 * g_z_xz_xzz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_zz_yy[i] = 4.0 * g_z_xz_xzz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_zz_yz[i] = 4.0 * g_z_xz_xzz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_zz_zz[i] = 4.0 * g_z_xz_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1404-1410)

    #pragma omp simd aligned(g_z_0_x_0_0_yy_xx_xx, g_z_0_x_0_0_yy_xx_xy, g_z_0_x_0_0_yy_xx_xz, g_z_0_x_0_0_yy_xx_yy, g_z_0_x_0_0_yy_xx_yz, g_z_0_x_0_0_yy_xx_zz, g_z_yy_x_xx, g_z_yy_x_xy, g_z_yy_x_xz, g_z_yy_x_yy, g_z_yy_x_yz, g_z_yy_x_zz, g_z_yy_xxx_xx, g_z_yy_xxx_xy, g_z_yy_xxx_xz, g_z_yy_xxx_yy, g_z_yy_xxx_yz, g_z_yy_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_yy_xx_xx[i] = -4.0 * g_z_yy_x_xx[i] * a_exp + 4.0 * g_z_yy_xxx_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_xx_xy[i] = -4.0 * g_z_yy_x_xy[i] * a_exp + 4.0 * g_z_yy_xxx_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_xx_xz[i] = -4.0 * g_z_yy_x_xz[i] * a_exp + 4.0 * g_z_yy_xxx_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_xx_yy[i] = -4.0 * g_z_yy_x_yy[i] * a_exp + 4.0 * g_z_yy_xxx_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_xx_yz[i] = -4.0 * g_z_yy_x_yz[i] * a_exp + 4.0 * g_z_yy_xxx_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_xx_zz[i] = -4.0 * g_z_yy_x_zz[i] * a_exp + 4.0 * g_z_yy_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1410-1416)

    #pragma omp simd aligned(g_z_0_x_0_0_yy_xy_xx, g_z_0_x_0_0_yy_xy_xy, g_z_0_x_0_0_yy_xy_xz, g_z_0_x_0_0_yy_xy_yy, g_z_0_x_0_0_yy_xy_yz, g_z_0_x_0_0_yy_xy_zz, g_z_yy_xxy_xx, g_z_yy_xxy_xy, g_z_yy_xxy_xz, g_z_yy_xxy_yy, g_z_yy_xxy_yz, g_z_yy_xxy_zz, g_z_yy_y_xx, g_z_yy_y_xy, g_z_yy_y_xz, g_z_yy_y_yy, g_z_yy_y_yz, g_z_yy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_yy_xy_xx[i] = -2.0 * g_z_yy_y_xx[i] * a_exp + 4.0 * g_z_yy_xxy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_xy_xy[i] = -2.0 * g_z_yy_y_xy[i] * a_exp + 4.0 * g_z_yy_xxy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_xy_xz[i] = -2.0 * g_z_yy_y_xz[i] * a_exp + 4.0 * g_z_yy_xxy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_xy_yy[i] = -2.0 * g_z_yy_y_yy[i] * a_exp + 4.0 * g_z_yy_xxy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_xy_yz[i] = -2.0 * g_z_yy_y_yz[i] * a_exp + 4.0 * g_z_yy_xxy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_xy_zz[i] = -2.0 * g_z_yy_y_zz[i] * a_exp + 4.0 * g_z_yy_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1416-1422)

    #pragma omp simd aligned(g_z_0_x_0_0_yy_xz_xx, g_z_0_x_0_0_yy_xz_xy, g_z_0_x_0_0_yy_xz_xz, g_z_0_x_0_0_yy_xz_yy, g_z_0_x_0_0_yy_xz_yz, g_z_0_x_0_0_yy_xz_zz, g_z_yy_xxz_xx, g_z_yy_xxz_xy, g_z_yy_xxz_xz, g_z_yy_xxz_yy, g_z_yy_xxz_yz, g_z_yy_xxz_zz, g_z_yy_z_xx, g_z_yy_z_xy, g_z_yy_z_xz, g_z_yy_z_yy, g_z_yy_z_yz, g_z_yy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_yy_xz_xx[i] = -2.0 * g_z_yy_z_xx[i] * a_exp + 4.0 * g_z_yy_xxz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_xz_xy[i] = -2.0 * g_z_yy_z_xy[i] * a_exp + 4.0 * g_z_yy_xxz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_xz_xz[i] = -2.0 * g_z_yy_z_xz[i] * a_exp + 4.0 * g_z_yy_xxz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_xz_yy[i] = -2.0 * g_z_yy_z_yy[i] * a_exp + 4.0 * g_z_yy_xxz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_xz_yz[i] = -2.0 * g_z_yy_z_yz[i] * a_exp + 4.0 * g_z_yy_xxz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_xz_zz[i] = -2.0 * g_z_yy_z_zz[i] * a_exp + 4.0 * g_z_yy_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1422-1428)

    #pragma omp simd aligned(g_z_0_x_0_0_yy_yy_xx, g_z_0_x_0_0_yy_yy_xy, g_z_0_x_0_0_yy_yy_xz, g_z_0_x_0_0_yy_yy_yy, g_z_0_x_0_0_yy_yy_yz, g_z_0_x_0_0_yy_yy_zz, g_z_yy_xyy_xx, g_z_yy_xyy_xy, g_z_yy_xyy_xz, g_z_yy_xyy_yy, g_z_yy_xyy_yz, g_z_yy_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_yy_yy_xx[i] = 4.0 * g_z_yy_xyy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_yy_xy[i] = 4.0 * g_z_yy_xyy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_yy_xz[i] = 4.0 * g_z_yy_xyy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_yy_yy[i] = 4.0 * g_z_yy_xyy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_yy_yz[i] = 4.0 * g_z_yy_xyy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_yy_zz[i] = 4.0 * g_z_yy_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1428-1434)

    #pragma omp simd aligned(g_z_0_x_0_0_yy_yz_xx, g_z_0_x_0_0_yy_yz_xy, g_z_0_x_0_0_yy_yz_xz, g_z_0_x_0_0_yy_yz_yy, g_z_0_x_0_0_yy_yz_yz, g_z_0_x_0_0_yy_yz_zz, g_z_yy_xyz_xx, g_z_yy_xyz_xy, g_z_yy_xyz_xz, g_z_yy_xyz_yy, g_z_yy_xyz_yz, g_z_yy_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_yy_yz_xx[i] = 4.0 * g_z_yy_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_yz_xy[i] = 4.0 * g_z_yy_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_yz_xz[i] = 4.0 * g_z_yy_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_yz_yy[i] = 4.0 * g_z_yy_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_yz_yz[i] = 4.0 * g_z_yy_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_yz_zz[i] = 4.0 * g_z_yy_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1434-1440)

    #pragma omp simd aligned(g_z_0_x_0_0_yy_zz_xx, g_z_0_x_0_0_yy_zz_xy, g_z_0_x_0_0_yy_zz_xz, g_z_0_x_0_0_yy_zz_yy, g_z_0_x_0_0_yy_zz_yz, g_z_0_x_0_0_yy_zz_zz, g_z_yy_xzz_xx, g_z_yy_xzz_xy, g_z_yy_xzz_xz, g_z_yy_xzz_yy, g_z_yy_xzz_yz, g_z_yy_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_yy_zz_xx[i] = 4.0 * g_z_yy_xzz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_zz_xy[i] = 4.0 * g_z_yy_xzz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_zz_xz[i] = 4.0 * g_z_yy_xzz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_zz_yy[i] = 4.0 * g_z_yy_xzz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_zz_yz[i] = 4.0 * g_z_yy_xzz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_zz_zz[i] = 4.0 * g_z_yy_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1440-1446)

    #pragma omp simd aligned(g_z_0_x_0_0_yz_xx_xx, g_z_0_x_0_0_yz_xx_xy, g_z_0_x_0_0_yz_xx_xz, g_z_0_x_0_0_yz_xx_yy, g_z_0_x_0_0_yz_xx_yz, g_z_0_x_0_0_yz_xx_zz, g_z_yz_x_xx, g_z_yz_x_xy, g_z_yz_x_xz, g_z_yz_x_yy, g_z_yz_x_yz, g_z_yz_x_zz, g_z_yz_xxx_xx, g_z_yz_xxx_xy, g_z_yz_xxx_xz, g_z_yz_xxx_yy, g_z_yz_xxx_yz, g_z_yz_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_yz_xx_xx[i] = -4.0 * g_z_yz_x_xx[i] * a_exp + 4.0 * g_z_yz_xxx_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_xx_xy[i] = -4.0 * g_z_yz_x_xy[i] * a_exp + 4.0 * g_z_yz_xxx_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_xx_xz[i] = -4.0 * g_z_yz_x_xz[i] * a_exp + 4.0 * g_z_yz_xxx_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_xx_yy[i] = -4.0 * g_z_yz_x_yy[i] * a_exp + 4.0 * g_z_yz_xxx_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_xx_yz[i] = -4.0 * g_z_yz_x_yz[i] * a_exp + 4.0 * g_z_yz_xxx_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_xx_zz[i] = -4.0 * g_z_yz_x_zz[i] * a_exp + 4.0 * g_z_yz_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1446-1452)

    #pragma omp simd aligned(g_z_0_x_0_0_yz_xy_xx, g_z_0_x_0_0_yz_xy_xy, g_z_0_x_0_0_yz_xy_xz, g_z_0_x_0_0_yz_xy_yy, g_z_0_x_0_0_yz_xy_yz, g_z_0_x_0_0_yz_xy_zz, g_z_yz_xxy_xx, g_z_yz_xxy_xy, g_z_yz_xxy_xz, g_z_yz_xxy_yy, g_z_yz_xxy_yz, g_z_yz_xxy_zz, g_z_yz_y_xx, g_z_yz_y_xy, g_z_yz_y_xz, g_z_yz_y_yy, g_z_yz_y_yz, g_z_yz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_yz_xy_xx[i] = -2.0 * g_z_yz_y_xx[i] * a_exp + 4.0 * g_z_yz_xxy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_xy_xy[i] = -2.0 * g_z_yz_y_xy[i] * a_exp + 4.0 * g_z_yz_xxy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_xy_xz[i] = -2.0 * g_z_yz_y_xz[i] * a_exp + 4.0 * g_z_yz_xxy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_xy_yy[i] = -2.0 * g_z_yz_y_yy[i] * a_exp + 4.0 * g_z_yz_xxy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_xy_yz[i] = -2.0 * g_z_yz_y_yz[i] * a_exp + 4.0 * g_z_yz_xxy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_xy_zz[i] = -2.0 * g_z_yz_y_zz[i] * a_exp + 4.0 * g_z_yz_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1452-1458)

    #pragma omp simd aligned(g_z_0_x_0_0_yz_xz_xx, g_z_0_x_0_0_yz_xz_xy, g_z_0_x_0_0_yz_xz_xz, g_z_0_x_0_0_yz_xz_yy, g_z_0_x_0_0_yz_xz_yz, g_z_0_x_0_0_yz_xz_zz, g_z_yz_xxz_xx, g_z_yz_xxz_xy, g_z_yz_xxz_xz, g_z_yz_xxz_yy, g_z_yz_xxz_yz, g_z_yz_xxz_zz, g_z_yz_z_xx, g_z_yz_z_xy, g_z_yz_z_xz, g_z_yz_z_yy, g_z_yz_z_yz, g_z_yz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_yz_xz_xx[i] = -2.0 * g_z_yz_z_xx[i] * a_exp + 4.0 * g_z_yz_xxz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_xz_xy[i] = -2.0 * g_z_yz_z_xy[i] * a_exp + 4.0 * g_z_yz_xxz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_xz_xz[i] = -2.0 * g_z_yz_z_xz[i] * a_exp + 4.0 * g_z_yz_xxz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_xz_yy[i] = -2.0 * g_z_yz_z_yy[i] * a_exp + 4.0 * g_z_yz_xxz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_xz_yz[i] = -2.0 * g_z_yz_z_yz[i] * a_exp + 4.0 * g_z_yz_xxz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_xz_zz[i] = -2.0 * g_z_yz_z_zz[i] * a_exp + 4.0 * g_z_yz_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1458-1464)

    #pragma omp simd aligned(g_z_0_x_0_0_yz_yy_xx, g_z_0_x_0_0_yz_yy_xy, g_z_0_x_0_0_yz_yy_xz, g_z_0_x_0_0_yz_yy_yy, g_z_0_x_0_0_yz_yy_yz, g_z_0_x_0_0_yz_yy_zz, g_z_yz_xyy_xx, g_z_yz_xyy_xy, g_z_yz_xyy_xz, g_z_yz_xyy_yy, g_z_yz_xyy_yz, g_z_yz_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_yz_yy_xx[i] = 4.0 * g_z_yz_xyy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_yy_xy[i] = 4.0 * g_z_yz_xyy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_yy_xz[i] = 4.0 * g_z_yz_xyy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_yy_yy[i] = 4.0 * g_z_yz_xyy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_yy_yz[i] = 4.0 * g_z_yz_xyy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_yy_zz[i] = 4.0 * g_z_yz_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1464-1470)

    #pragma omp simd aligned(g_z_0_x_0_0_yz_yz_xx, g_z_0_x_0_0_yz_yz_xy, g_z_0_x_0_0_yz_yz_xz, g_z_0_x_0_0_yz_yz_yy, g_z_0_x_0_0_yz_yz_yz, g_z_0_x_0_0_yz_yz_zz, g_z_yz_xyz_xx, g_z_yz_xyz_xy, g_z_yz_xyz_xz, g_z_yz_xyz_yy, g_z_yz_xyz_yz, g_z_yz_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_yz_yz_xx[i] = 4.0 * g_z_yz_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_yz_xy[i] = 4.0 * g_z_yz_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_yz_xz[i] = 4.0 * g_z_yz_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_yz_yy[i] = 4.0 * g_z_yz_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_yz_yz[i] = 4.0 * g_z_yz_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_yz_zz[i] = 4.0 * g_z_yz_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1470-1476)

    #pragma omp simd aligned(g_z_0_x_0_0_yz_zz_xx, g_z_0_x_0_0_yz_zz_xy, g_z_0_x_0_0_yz_zz_xz, g_z_0_x_0_0_yz_zz_yy, g_z_0_x_0_0_yz_zz_yz, g_z_0_x_0_0_yz_zz_zz, g_z_yz_xzz_xx, g_z_yz_xzz_xy, g_z_yz_xzz_xz, g_z_yz_xzz_yy, g_z_yz_xzz_yz, g_z_yz_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_yz_zz_xx[i] = 4.0 * g_z_yz_xzz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_zz_xy[i] = 4.0 * g_z_yz_xzz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_zz_xz[i] = 4.0 * g_z_yz_xzz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_zz_yy[i] = 4.0 * g_z_yz_xzz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_zz_yz[i] = 4.0 * g_z_yz_xzz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_zz_zz[i] = 4.0 * g_z_yz_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1476-1482)

    #pragma omp simd aligned(g_z_0_x_0_0_zz_xx_xx, g_z_0_x_0_0_zz_xx_xy, g_z_0_x_0_0_zz_xx_xz, g_z_0_x_0_0_zz_xx_yy, g_z_0_x_0_0_zz_xx_yz, g_z_0_x_0_0_zz_xx_zz, g_z_zz_x_xx, g_z_zz_x_xy, g_z_zz_x_xz, g_z_zz_x_yy, g_z_zz_x_yz, g_z_zz_x_zz, g_z_zz_xxx_xx, g_z_zz_xxx_xy, g_z_zz_xxx_xz, g_z_zz_xxx_yy, g_z_zz_xxx_yz, g_z_zz_xxx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_zz_xx_xx[i] = -4.0 * g_z_zz_x_xx[i] * a_exp + 4.0 * g_z_zz_xxx_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_xx_xy[i] = -4.0 * g_z_zz_x_xy[i] * a_exp + 4.0 * g_z_zz_xxx_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_xx_xz[i] = -4.0 * g_z_zz_x_xz[i] * a_exp + 4.0 * g_z_zz_xxx_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_xx_yy[i] = -4.0 * g_z_zz_x_yy[i] * a_exp + 4.0 * g_z_zz_xxx_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_xx_yz[i] = -4.0 * g_z_zz_x_yz[i] * a_exp + 4.0 * g_z_zz_xxx_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_xx_zz[i] = -4.0 * g_z_zz_x_zz[i] * a_exp + 4.0 * g_z_zz_xxx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1482-1488)

    #pragma omp simd aligned(g_z_0_x_0_0_zz_xy_xx, g_z_0_x_0_0_zz_xy_xy, g_z_0_x_0_0_zz_xy_xz, g_z_0_x_0_0_zz_xy_yy, g_z_0_x_0_0_zz_xy_yz, g_z_0_x_0_0_zz_xy_zz, g_z_zz_xxy_xx, g_z_zz_xxy_xy, g_z_zz_xxy_xz, g_z_zz_xxy_yy, g_z_zz_xxy_yz, g_z_zz_xxy_zz, g_z_zz_y_xx, g_z_zz_y_xy, g_z_zz_y_xz, g_z_zz_y_yy, g_z_zz_y_yz, g_z_zz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_zz_xy_xx[i] = -2.0 * g_z_zz_y_xx[i] * a_exp + 4.0 * g_z_zz_xxy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_xy_xy[i] = -2.0 * g_z_zz_y_xy[i] * a_exp + 4.0 * g_z_zz_xxy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_xy_xz[i] = -2.0 * g_z_zz_y_xz[i] * a_exp + 4.0 * g_z_zz_xxy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_xy_yy[i] = -2.0 * g_z_zz_y_yy[i] * a_exp + 4.0 * g_z_zz_xxy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_xy_yz[i] = -2.0 * g_z_zz_y_yz[i] * a_exp + 4.0 * g_z_zz_xxy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_xy_zz[i] = -2.0 * g_z_zz_y_zz[i] * a_exp + 4.0 * g_z_zz_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1488-1494)

    #pragma omp simd aligned(g_z_0_x_0_0_zz_xz_xx, g_z_0_x_0_0_zz_xz_xy, g_z_0_x_0_0_zz_xz_xz, g_z_0_x_0_0_zz_xz_yy, g_z_0_x_0_0_zz_xz_yz, g_z_0_x_0_0_zz_xz_zz, g_z_zz_xxz_xx, g_z_zz_xxz_xy, g_z_zz_xxz_xz, g_z_zz_xxz_yy, g_z_zz_xxz_yz, g_z_zz_xxz_zz, g_z_zz_z_xx, g_z_zz_z_xy, g_z_zz_z_xz, g_z_zz_z_yy, g_z_zz_z_yz, g_z_zz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_zz_xz_xx[i] = -2.0 * g_z_zz_z_xx[i] * a_exp + 4.0 * g_z_zz_xxz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_xz_xy[i] = -2.0 * g_z_zz_z_xy[i] * a_exp + 4.0 * g_z_zz_xxz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_xz_xz[i] = -2.0 * g_z_zz_z_xz[i] * a_exp + 4.0 * g_z_zz_xxz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_xz_yy[i] = -2.0 * g_z_zz_z_yy[i] * a_exp + 4.0 * g_z_zz_xxz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_xz_yz[i] = -2.0 * g_z_zz_z_yz[i] * a_exp + 4.0 * g_z_zz_xxz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_xz_zz[i] = -2.0 * g_z_zz_z_zz[i] * a_exp + 4.0 * g_z_zz_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1494-1500)

    #pragma omp simd aligned(g_z_0_x_0_0_zz_yy_xx, g_z_0_x_0_0_zz_yy_xy, g_z_0_x_0_0_zz_yy_xz, g_z_0_x_0_0_zz_yy_yy, g_z_0_x_0_0_zz_yy_yz, g_z_0_x_0_0_zz_yy_zz, g_z_zz_xyy_xx, g_z_zz_xyy_xy, g_z_zz_xyy_xz, g_z_zz_xyy_yy, g_z_zz_xyy_yz, g_z_zz_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_zz_yy_xx[i] = 4.0 * g_z_zz_xyy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_yy_xy[i] = 4.0 * g_z_zz_xyy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_yy_xz[i] = 4.0 * g_z_zz_xyy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_yy_yy[i] = 4.0 * g_z_zz_xyy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_yy_yz[i] = 4.0 * g_z_zz_xyy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_yy_zz[i] = 4.0 * g_z_zz_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1500-1506)

    #pragma omp simd aligned(g_z_0_x_0_0_zz_yz_xx, g_z_0_x_0_0_zz_yz_xy, g_z_0_x_0_0_zz_yz_xz, g_z_0_x_0_0_zz_yz_yy, g_z_0_x_0_0_zz_yz_yz, g_z_0_x_0_0_zz_yz_zz, g_z_zz_xyz_xx, g_z_zz_xyz_xy, g_z_zz_xyz_xz, g_z_zz_xyz_yy, g_z_zz_xyz_yz, g_z_zz_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_zz_yz_xx[i] = 4.0 * g_z_zz_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_yz_xy[i] = 4.0 * g_z_zz_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_yz_xz[i] = 4.0 * g_z_zz_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_yz_yy[i] = 4.0 * g_z_zz_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_yz_yz[i] = 4.0 * g_z_zz_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_yz_zz[i] = 4.0 * g_z_zz_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1506-1512)

    #pragma omp simd aligned(g_z_0_x_0_0_zz_zz_xx, g_z_0_x_0_0_zz_zz_xy, g_z_0_x_0_0_zz_zz_xz, g_z_0_x_0_0_zz_zz_yy, g_z_0_x_0_0_zz_zz_yz, g_z_0_x_0_0_zz_zz_zz, g_z_zz_xzz_xx, g_z_zz_xzz_xy, g_z_zz_xzz_xz, g_z_zz_xzz_yy, g_z_zz_xzz_yz, g_z_zz_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_zz_zz_xx[i] = 4.0 * g_z_zz_xzz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_zz_xy[i] = 4.0 * g_z_zz_xzz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_zz_xz[i] = 4.0 * g_z_zz_xzz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_zz_yy[i] = 4.0 * g_z_zz_xzz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_zz_yz[i] = 4.0 * g_z_zz_xzz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_zz_zz[i] = 4.0 * g_z_zz_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1512-1518)

    #pragma omp simd aligned(g_z_0_y_0_0_xx_xx_xx, g_z_0_y_0_0_xx_xx_xy, g_z_0_y_0_0_xx_xx_xz, g_z_0_y_0_0_xx_xx_yy, g_z_0_y_0_0_xx_xx_yz, g_z_0_y_0_0_xx_xx_zz, g_z_xx_xxy_xx, g_z_xx_xxy_xy, g_z_xx_xxy_xz, g_z_xx_xxy_yy, g_z_xx_xxy_yz, g_z_xx_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xx_xx_xx[i] = 4.0 * g_z_xx_xxy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_xx_xy[i] = 4.0 * g_z_xx_xxy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_xx_xz[i] = 4.0 * g_z_xx_xxy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_xx_yy[i] = 4.0 * g_z_xx_xxy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_xx_yz[i] = 4.0 * g_z_xx_xxy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_xx_zz[i] = 4.0 * g_z_xx_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1518-1524)

    #pragma omp simd aligned(g_z_0_y_0_0_xx_xy_xx, g_z_0_y_0_0_xx_xy_xy, g_z_0_y_0_0_xx_xy_xz, g_z_0_y_0_0_xx_xy_yy, g_z_0_y_0_0_xx_xy_yz, g_z_0_y_0_0_xx_xy_zz, g_z_xx_x_xx, g_z_xx_x_xy, g_z_xx_x_xz, g_z_xx_x_yy, g_z_xx_x_yz, g_z_xx_x_zz, g_z_xx_xyy_xx, g_z_xx_xyy_xy, g_z_xx_xyy_xz, g_z_xx_xyy_yy, g_z_xx_xyy_yz, g_z_xx_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xx_xy_xx[i] = -2.0 * g_z_xx_x_xx[i] * a_exp + 4.0 * g_z_xx_xyy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_xy_xy[i] = -2.0 * g_z_xx_x_xy[i] * a_exp + 4.0 * g_z_xx_xyy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_xy_xz[i] = -2.0 * g_z_xx_x_xz[i] * a_exp + 4.0 * g_z_xx_xyy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_xy_yy[i] = -2.0 * g_z_xx_x_yy[i] * a_exp + 4.0 * g_z_xx_xyy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_xy_yz[i] = -2.0 * g_z_xx_x_yz[i] * a_exp + 4.0 * g_z_xx_xyy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_xy_zz[i] = -2.0 * g_z_xx_x_zz[i] * a_exp + 4.0 * g_z_xx_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1524-1530)

    #pragma omp simd aligned(g_z_0_y_0_0_xx_xz_xx, g_z_0_y_0_0_xx_xz_xy, g_z_0_y_0_0_xx_xz_xz, g_z_0_y_0_0_xx_xz_yy, g_z_0_y_0_0_xx_xz_yz, g_z_0_y_0_0_xx_xz_zz, g_z_xx_xyz_xx, g_z_xx_xyz_xy, g_z_xx_xyz_xz, g_z_xx_xyz_yy, g_z_xx_xyz_yz, g_z_xx_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xx_xz_xx[i] = 4.0 * g_z_xx_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_xz_xy[i] = 4.0 * g_z_xx_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_xz_xz[i] = 4.0 * g_z_xx_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_xz_yy[i] = 4.0 * g_z_xx_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_xz_yz[i] = 4.0 * g_z_xx_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_xz_zz[i] = 4.0 * g_z_xx_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1530-1536)

    #pragma omp simd aligned(g_z_0_y_0_0_xx_yy_xx, g_z_0_y_0_0_xx_yy_xy, g_z_0_y_0_0_xx_yy_xz, g_z_0_y_0_0_xx_yy_yy, g_z_0_y_0_0_xx_yy_yz, g_z_0_y_0_0_xx_yy_zz, g_z_xx_y_xx, g_z_xx_y_xy, g_z_xx_y_xz, g_z_xx_y_yy, g_z_xx_y_yz, g_z_xx_y_zz, g_z_xx_yyy_xx, g_z_xx_yyy_xy, g_z_xx_yyy_xz, g_z_xx_yyy_yy, g_z_xx_yyy_yz, g_z_xx_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xx_yy_xx[i] = -4.0 * g_z_xx_y_xx[i] * a_exp + 4.0 * g_z_xx_yyy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_yy_xy[i] = -4.0 * g_z_xx_y_xy[i] * a_exp + 4.0 * g_z_xx_yyy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_yy_xz[i] = -4.0 * g_z_xx_y_xz[i] * a_exp + 4.0 * g_z_xx_yyy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_yy_yy[i] = -4.0 * g_z_xx_y_yy[i] * a_exp + 4.0 * g_z_xx_yyy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_yy_yz[i] = -4.0 * g_z_xx_y_yz[i] * a_exp + 4.0 * g_z_xx_yyy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_yy_zz[i] = -4.0 * g_z_xx_y_zz[i] * a_exp + 4.0 * g_z_xx_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1536-1542)

    #pragma omp simd aligned(g_z_0_y_0_0_xx_yz_xx, g_z_0_y_0_0_xx_yz_xy, g_z_0_y_0_0_xx_yz_xz, g_z_0_y_0_0_xx_yz_yy, g_z_0_y_0_0_xx_yz_yz, g_z_0_y_0_0_xx_yz_zz, g_z_xx_yyz_xx, g_z_xx_yyz_xy, g_z_xx_yyz_xz, g_z_xx_yyz_yy, g_z_xx_yyz_yz, g_z_xx_yyz_zz, g_z_xx_z_xx, g_z_xx_z_xy, g_z_xx_z_xz, g_z_xx_z_yy, g_z_xx_z_yz, g_z_xx_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xx_yz_xx[i] = -2.0 * g_z_xx_z_xx[i] * a_exp + 4.0 * g_z_xx_yyz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_yz_xy[i] = -2.0 * g_z_xx_z_xy[i] * a_exp + 4.0 * g_z_xx_yyz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_yz_xz[i] = -2.0 * g_z_xx_z_xz[i] * a_exp + 4.0 * g_z_xx_yyz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_yz_yy[i] = -2.0 * g_z_xx_z_yy[i] * a_exp + 4.0 * g_z_xx_yyz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_yz_yz[i] = -2.0 * g_z_xx_z_yz[i] * a_exp + 4.0 * g_z_xx_yyz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_yz_zz[i] = -2.0 * g_z_xx_z_zz[i] * a_exp + 4.0 * g_z_xx_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1542-1548)

    #pragma omp simd aligned(g_z_0_y_0_0_xx_zz_xx, g_z_0_y_0_0_xx_zz_xy, g_z_0_y_0_0_xx_zz_xz, g_z_0_y_0_0_xx_zz_yy, g_z_0_y_0_0_xx_zz_yz, g_z_0_y_0_0_xx_zz_zz, g_z_xx_yzz_xx, g_z_xx_yzz_xy, g_z_xx_yzz_xz, g_z_xx_yzz_yy, g_z_xx_yzz_yz, g_z_xx_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xx_zz_xx[i] = 4.0 * g_z_xx_yzz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_zz_xy[i] = 4.0 * g_z_xx_yzz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_zz_xz[i] = 4.0 * g_z_xx_yzz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_zz_yy[i] = 4.0 * g_z_xx_yzz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_zz_yz[i] = 4.0 * g_z_xx_yzz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_zz_zz[i] = 4.0 * g_z_xx_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1548-1554)

    #pragma omp simd aligned(g_z_0_y_0_0_xy_xx_xx, g_z_0_y_0_0_xy_xx_xy, g_z_0_y_0_0_xy_xx_xz, g_z_0_y_0_0_xy_xx_yy, g_z_0_y_0_0_xy_xx_yz, g_z_0_y_0_0_xy_xx_zz, g_z_xy_xxy_xx, g_z_xy_xxy_xy, g_z_xy_xxy_xz, g_z_xy_xxy_yy, g_z_xy_xxy_yz, g_z_xy_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xy_xx_xx[i] = 4.0 * g_z_xy_xxy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_xx_xy[i] = 4.0 * g_z_xy_xxy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_xx_xz[i] = 4.0 * g_z_xy_xxy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_xx_yy[i] = 4.0 * g_z_xy_xxy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_xx_yz[i] = 4.0 * g_z_xy_xxy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_xx_zz[i] = 4.0 * g_z_xy_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1554-1560)

    #pragma omp simd aligned(g_z_0_y_0_0_xy_xy_xx, g_z_0_y_0_0_xy_xy_xy, g_z_0_y_0_0_xy_xy_xz, g_z_0_y_0_0_xy_xy_yy, g_z_0_y_0_0_xy_xy_yz, g_z_0_y_0_0_xy_xy_zz, g_z_xy_x_xx, g_z_xy_x_xy, g_z_xy_x_xz, g_z_xy_x_yy, g_z_xy_x_yz, g_z_xy_x_zz, g_z_xy_xyy_xx, g_z_xy_xyy_xy, g_z_xy_xyy_xz, g_z_xy_xyy_yy, g_z_xy_xyy_yz, g_z_xy_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xy_xy_xx[i] = -2.0 * g_z_xy_x_xx[i] * a_exp + 4.0 * g_z_xy_xyy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_xy_xy[i] = -2.0 * g_z_xy_x_xy[i] * a_exp + 4.0 * g_z_xy_xyy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_xy_xz[i] = -2.0 * g_z_xy_x_xz[i] * a_exp + 4.0 * g_z_xy_xyy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_xy_yy[i] = -2.0 * g_z_xy_x_yy[i] * a_exp + 4.0 * g_z_xy_xyy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_xy_yz[i] = -2.0 * g_z_xy_x_yz[i] * a_exp + 4.0 * g_z_xy_xyy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_xy_zz[i] = -2.0 * g_z_xy_x_zz[i] * a_exp + 4.0 * g_z_xy_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1560-1566)

    #pragma omp simd aligned(g_z_0_y_0_0_xy_xz_xx, g_z_0_y_0_0_xy_xz_xy, g_z_0_y_0_0_xy_xz_xz, g_z_0_y_0_0_xy_xz_yy, g_z_0_y_0_0_xy_xz_yz, g_z_0_y_0_0_xy_xz_zz, g_z_xy_xyz_xx, g_z_xy_xyz_xy, g_z_xy_xyz_xz, g_z_xy_xyz_yy, g_z_xy_xyz_yz, g_z_xy_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xy_xz_xx[i] = 4.0 * g_z_xy_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_xz_xy[i] = 4.0 * g_z_xy_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_xz_xz[i] = 4.0 * g_z_xy_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_xz_yy[i] = 4.0 * g_z_xy_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_xz_yz[i] = 4.0 * g_z_xy_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_xz_zz[i] = 4.0 * g_z_xy_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1566-1572)

    #pragma omp simd aligned(g_z_0_y_0_0_xy_yy_xx, g_z_0_y_0_0_xy_yy_xy, g_z_0_y_0_0_xy_yy_xz, g_z_0_y_0_0_xy_yy_yy, g_z_0_y_0_0_xy_yy_yz, g_z_0_y_0_0_xy_yy_zz, g_z_xy_y_xx, g_z_xy_y_xy, g_z_xy_y_xz, g_z_xy_y_yy, g_z_xy_y_yz, g_z_xy_y_zz, g_z_xy_yyy_xx, g_z_xy_yyy_xy, g_z_xy_yyy_xz, g_z_xy_yyy_yy, g_z_xy_yyy_yz, g_z_xy_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xy_yy_xx[i] = -4.0 * g_z_xy_y_xx[i] * a_exp + 4.0 * g_z_xy_yyy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_yy_xy[i] = -4.0 * g_z_xy_y_xy[i] * a_exp + 4.0 * g_z_xy_yyy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_yy_xz[i] = -4.0 * g_z_xy_y_xz[i] * a_exp + 4.0 * g_z_xy_yyy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_yy_yy[i] = -4.0 * g_z_xy_y_yy[i] * a_exp + 4.0 * g_z_xy_yyy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_yy_yz[i] = -4.0 * g_z_xy_y_yz[i] * a_exp + 4.0 * g_z_xy_yyy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_yy_zz[i] = -4.0 * g_z_xy_y_zz[i] * a_exp + 4.0 * g_z_xy_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1572-1578)

    #pragma omp simd aligned(g_z_0_y_0_0_xy_yz_xx, g_z_0_y_0_0_xy_yz_xy, g_z_0_y_0_0_xy_yz_xz, g_z_0_y_0_0_xy_yz_yy, g_z_0_y_0_0_xy_yz_yz, g_z_0_y_0_0_xy_yz_zz, g_z_xy_yyz_xx, g_z_xy_yyz_xy, g_z_xy_yyz_xz, g_z_xy_yyz_yy, g_z_xy_yyz_yz, g_z_xy_yyz_zz, g_z_xy_z_xx, g_z_xy_z_xy, g_z_xy_z_xz, g_z_xy_z_yy, g_z_xy_z_yz, g_z_xy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xy_yz_xx[i] = -2.0 * g_z_xy_z_xx[i] * a_exp + 4.0 * g_z_xy_yyz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_yz_xy[i] = -2.0 * g_z_xy_z_xy[i] * a_exp + 4.0 * g_z_xy_yyz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_yz_xz[i] = -2.0 * g_z_xy_z_xz[i] * a_exp + 4.0 * g_z_xy_yyz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_yz_yy[i] = -2.0 * g_z_xy_z_yy[i] * a_exp + 4.0 * g_z_xy_yyz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_yz_yz[i] = -2.0 * g_z_xy_z_yz[i] * a_exp + 4.0 * g_z_xy_yyz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_yz_zz[i] = -2.0 * g_z_xy_z_zz[i] * a_exp + 4.0 * g_z_xy_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1578-1584)

    #pragma omp simd aligned(g_z_0_y_0_0_xy_zz_xx, g_z_0_y_0_0_xy_zz_xy, g_z_0_y_0_0_xy_zz_xz, g_z_0_y_0_0_xy_zz_yy, g_z_0_y_0_0_xy_zz_yz, g_z_0_y_0_0_xy_zz_zz, g_z_xy_yzz_xx, g_z_xy_yzz_xy, g_z_xy_yzz_xz, g_z_xy_yzz_yy, g_z_xy_yzz_yz, g_z_xy_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xy_zz_xx[i] = 4.0 * g_z_xy_yzz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_zz_xy[i] = 4.0 * g_z_xy_yzz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_zz_xz[i] = 4.0 * g_z_xy_yzz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_zz_yy[i] = 4.0 * g_z_xy_yzz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_zz_yz[i] = 4.0 * g_z_xy_yzz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_zz_zz[i] = 4.0 * g_z_xy_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1584-1590)

    #pragma omp simd aligned(g_z_0_y_0_0_xz_xx_xx, g_z_0_y_0_0_xz_xx_xy, g_z_0_y_0_0_xz_xx_xz, g_z_0_y_0_0_xz_xx_yy, g_z_0_y_0_0_xz_xx_yz, g_z_0_y_0_0_xz_xx_zz, g_z_xz_xxy_xx, g_z_xz_xxy_xy, g_z_xz_xxy_xz, g_z_xz_xxy_yy, g_z_xz_xxy_yz, g_z_xz_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xz_xx_xx[i] = 4.0 * g_z_xz_xxy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_xx_xy[i] = 4.0 * g_z_xz_xxy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_xx_xz[i] = 4.0 * g_z_xz_xxy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_xx_yy[i] = 4.0 * g_z_xz_xxy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_xx_yz[i] = 4.0 * g_z_xz_xxy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_xx_zz[i] = 4.0 * g_z_xz_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1590-1596)

    #pragma omp simd aligned(g_z_0_y_0_0_xz_xy_xx, g_z_0_y_0_0_xz_xy_xy, g_z_0_y_0_0_xz_xy_xz, g_z_0_y_0_0_xz_xy_yy, g_z_0_y_0_0_xz_xy_yz, g_z_0_y_0_0_xz_xy_zz, g_z_xz_x_xx, g_z_xz_x_xy, g_z_xz_x_xz, g_z_xz_x_yy, g_z_xz_x_yz, g_z_xz_x_zz, g_z_xz_xyy_xx, g_z_xz_xyy_xy, g_z_xz_xyy_xz, g_z_xz_xyy_yy, g_z_xz_xyy_yz, g_z_xz_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xz_xy_xx[i] = -2.0 * g_z_xz_x_xx[i] * a_exp + 4.0 * g_z_xz_xyy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_xy_xy[i] = -2.0 * g_z_xz_x_xy[i] * a_exp + 4.0 * g_z_xz_xyy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_xy_xz[i] = -2.0 * g_z_xz_x_xz[i] * a_exp + 4.0 * g_z_xz_xyy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_xy_yy[i] = -2.0 * g_z_xz_x_yy[i] * a_exp + 4.0 * g_z_xz_xyy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_xy_yz[i] = -2.0 * g_z_xz_x_yz[i] * a_exp + 4.0 * g_z_xz_xyy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_xy_zz[i] = -2.0 * g_z_xz_x_zz[i] * a_exp + 4.0 * g_z_xz_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1596-1602)

    #pragma omp simd aligned(g_z_0_y_0_0_xz_xz_xx, g_z_0_y_0_0_xz_xz_xy, g_z_0_y_0_0_xz_xz_xz, g_z_0_y_0_0_xz_xz_yy, g_z_0_y_0_0_xz_xz_yz, g_z_0_y_0_0_xz_xz_zz, g_z_xz_xyz_xx, g_z_xz_xyz_xy, g_z_xz_xyz_xz, g_z_xz_xyz_yy, g_z_xz_xyz_yz, g_z_xz_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xz_xz_xx[i] = 4.0 * g_z_xz_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_xz_xy[i] = 4.0 * g_z_xz_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_xz_xz[i] = 4.0 * g_z_xz_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_xz_yy[i] = 4.0 * g_z_xz_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_xz_yz[i] = 4.0 * g_z_xz_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_xz_zz[i] = 4.0 * g_z_xz_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1602-1608)

    #pragma omp simd aligned(g_z_0_y_0_0_xz_yy_xx, g_z_0_y_0_0_xz_yy_xy, g_z_0_y_0_0_xz_yy_xz, g_z_0_y_0_0_xz_yy_yy, g_z_0_y_0_0_xz_yy_yz, g_z_0_y_0_0_xz_yy_zz, g_z_xz_y_xx, g_z_xz_y_xy, g_z_xz_y_xz, g_z_xz_y_yy, g_z_xz_y_yz, g_z_xz_y_zz, g_z_xz_yyy_xx, g_z_xz_yyy_xy, g_z_xz_yyy_xz, g_z_xz_yyy_yy, g_z_xz_yyy_yz, g_z_xz_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xz_yy_xx[i] = -4.0 * g_z_xz_y_xx[i] * a_exp + 4.0 * g_z_xz_yyy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_yy_xy[i] = -4.0 * g_z_xz_y_xy[i] * a_exp + 4.0 * g_z_xz_yyy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_yy_xz[i] = -4.0 * g_z_xz_y_xz[i] * a_exp + 4.0 * g_z_xz_yyy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_yy_yy[i] = -4.0 * g_z_xz_y_yy[i] * a_exp + 4.0 * g_z_xz_yyy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_yy_yz[i] = -4.0 * g_z_xz_y_yz[i] * a_exp + 4.0 * g_z_xz_yyy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_yy_zz[i] = -4.0 * g_z_xz_y_zz[i] * a_exp + 4.0 * g_z_xz_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1608-1614)

    #pragma omp simd aligned(g_z_0_y_0_0_xz_yz_xx, g_z_0_y_0_0_xz_yz_xy, g_z_0_y_0_0_xz_yz_xz, g_z_0_y_0_0_xz_yz_yy, g_z_0_y_0_0_xz_yz_yz, g_z_0_y_0_0_xz_yz_zz, g_z_xz_yyz_xx, g_z_xz_yyz_xy, g_z_xz_yyz_xz, g_z_xz_yyz_yy, g_z_xz_yyz_yz, g_z_xz_yyz_zz, g_z_xz_z_xx, g_z_xz_z_xy, g_z_xz_z_xz, g_z_xz_z_yy, g_z_xz_z_yz, g_z_xz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xz_yz_xx[i] = -2.0 * g_z_xz_z_xx[i] * a_exp + 4.0 * g_z_xz_yyz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_yz_xy[i] = -2.0 * g_z_xz_z_xy[i] * a_exp + 4.0 * g_z_xz_yyz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_yz_xz[i] = -2.0 * g_z_xz_z_xz[i] * a_exp + 4.0 * g_z_xz_yyz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_yz_yy[i] = -2.0 * g_z_xz_z_yy[i] * a_exp + 4.0 * g_z_xz_yyz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_yz_yz[i] = -2.0 * g_z_xz_z_yz[i] * a_exp + 4.0 * g_z_xz_yyz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_yz_zz[i] = -2.0 * g_z_xz_z_zz[i] * a_exp + 4.0 * g_z_xz_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1614-1620)

    #pragma omp simd aligned(g_z_0_y_0_0_xz_zz_xx, g_z_0_y_0_0_xz_zz_xy, g_z_0_y_0_0_xz_zz_xz, g_z_0_y_0_0_xz_zz_yy, g_z_0_y_0_0_xz_zz_yz, g_z_0_y_0_0_xz_zz_zz, g_z_xz_yzz_xx, g_z_xz_yzz_xy, g_z_xz_yzz_xz, g_z_xz_yzz_yy, g_z_xz_yzz_yz, g_z_xz_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xz_zz_xx[i] = 4.0 * g_z_xz_yzz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_zz_xy[i] = 4.0 * g_z_xz_yzz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_zz_xz[i] = 4.0 * g_z_xz_yzz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_zz_yy[i] = 4.0 * g_z_xz_yzz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_zz_yz[i] = 4.0 * g_z_xz_yzz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_zz_zz[i] = 4.0 * g_z_xz_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1620-1626)

    #pragma omp simd aligned(g_z_0_y_0_0_yy_xx_xx, g_z_0_y_0_0_yy_xx_xy, g_z_0_y_0_0_yy_xx_xz, g_z_0_y_0_0_yy_xx_yy, g_z_0_y_0_0_yy_xx_yz, g_z_0_y_0_0_yy_xx_zz, g_z_yy_xxy_xx, g_z_yy_xxy_xy, g_z_yy_xxy_xz, g_z_yy_xxy_yy, g_z_yy_xxy_yz, g_z_yy_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_yy_xx_xx[i] = 4.0 * g_z_yy_xxy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_xx_xy[i] = 4.0 * g_z_yy_xxy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_xx_xz[i] = 4.0 * g_z_yy_xxy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_xx_yy[i] = 4.0 * g_z_yy_xxy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_xx_yz[i] = 4.0 * g_z_yy_xxy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_xx_zz[i] = 4.0 * g_z_yy_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1626-1632)

    #pragma omp simd aligned(g_z_0_y_0_0_yy_xy_xx, g_z_0_y_0_0_yy_xy_xy, g_z_0_y_0_0_yy_xy_xz, g_z_0_y_0_0_yy_xy_yy, g_z_0_y_0_0_yy_xy_yz, g_z_0_y_0_0_yy_xy_zz, g_z_yy_x_xx, g_z_yy_x_xy, g_z_yy_x_xz, g_z_yy_x_yy, g_z_yy_x_yz, g_z_yy_x_zz, g_z_yy_xyy_xx, g_z_yy_xyy_xy, g_z_yy_xyy_xz, g_z_yy_xyy_yy, g_z_yy_xyy_yz, g_z_yy_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_yy_xy_xx[i] = -2.0 * g_z_yy_x_xx[i] * a_exp + 4.0 * g_z_yy_xyy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_xy_xy[i] = -2.0 * g_z_yy_x_xy[i] * a_exp + 4.0 * g_z_yy_xyy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_xy_xz[i] = -2.0 * g_z_yy_x_xz[i] * a_exp + 4.0 * g_z_yy_xyy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_xy_yy[i] = -2.0 * g_z_yy_x_yy[i] * a_exp + 4.0 * g_z_yy_xyy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_xy_yz[i] = -2.0 * g_z_yy_x_yz[i] * a_exp + 4.0 * g_z_yy_xyy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_xy_zz[i] = -2.0 * g_z_yy_x_zz[i] * a_exp + 4.0 * g_z_yy_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1632-1638)

    #pragma omp simd aligned(g_z_0_y_0_0_yy_xz_xx, g_z_0_y_0_0_yy_xz_xy, g_z_0_y_0_0_yy_xz_xz, g_z_0_y_0_0_yy_xz_yy, g_z_0_y_0_0_yy_xz_yz, g_z_0_y_0_0_yy_xz_zz, g_z_yy_xyz_xx, g_z_yy_xyz_xy, g_z_yy_xyz_xz, g_z_yy_xyz_yy, g_z_yy_xyz_yz, g_z_yy_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_yy_xz_xx[i] = 4.0 * g_z_yy_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_xz_xy[i] = 4.0 * g_z_yy_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_xz_xz[i] = 4.0 * g_z_yy_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_xz_yy[i] = 4.0 * g_z_yy_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_xz_yz[i] = 4.0 * g_z_yy_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_xz_zz[i] = 4.0 * g_z_yy_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1638-1644)

    #pragma omp simd aligned(g_z_0_y_0_0_yy_yy_xx, g_z_0_y_0_0_yy_yy_xy, g_z_0_y_0_0_yy_yy_xz, g_z_0_y_0_0_yy_yy_yy, g_z_0_y_0_0_yy_yy_yz, g_z_0_y_0_0_yy_yy_zz, g_z_yy_y_xx, g_z_yy_y_xy, g_z_yy_y_xz, g_z_yy_y_yy, g_z_yy_y_yz, g_z_yy_y_zz, g_z_yy_yyy_xx, g_z_yy_yyy_xy, g_z_yy_yyy_xz, g_z_yy_yyy_yy, g_z_yy_yyy_yz, g_z_yy_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_yy_yy_xx[i] = -4.0 * g_z_yy_y_xx[i] * a_exp + 4.0 * g_z_yy_yyy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_yy_xy[i] = -4.0 * g_z_yy_y_xy[i] * a_exp + 4.0 * g_z_yy_yyy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_yy_xz[i] = -4.0 * g_z_yy_y_xz[i] * a_exp + 4.0 * g_z_yy_yyy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_yy_yy[i] = -4.0 * g_z_yy_y_yy[i] * a_exp + 4.0 * g_z_yy_yyy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_yy_yz[i] = -4.0 * g_z_yy_y_yz[i] * a_exp + 4.0 * g_z_yy_yyy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_yy_zz[i] = -4.0 * g_z_yy_y_zz[i] * a_exp + 4.0 * g_z_yy_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1644-1650)

    #pragma omp simd aligned(g_z_0_y_0_0_yy_yz_xx, g_z_0_y_0_0_yy_yz_xy, g_z_0_y_0_0_yy_yz_xz, g_z_0_y_0_0_yy_yz_yy, g_z_0_y_0_0_yy_yz_yz, g_z_0_y_0_0_yy_yz_zz, g_z_yy_yyz_xx, g_z_yy_yyz_xy, g_z_yy_yyz_xz, g_z_yy_yyz_yy, g_z_yy_yyz_yz, g_z_yy_yyz_zz, g_z_yy_z_xx, g_z_yy_z_xy, g_z_yy_z_xz, g_z_yy_z_yy, g_z_yy_z_yz, g_z_yy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_yy_yz_xx[i] = -2.0 * g_z_yy_z_xx[i] * a_exp + 4.0 * g_z_yy_yyz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_yz_xy[i] = -2.0 * g_z_yy_z_xy[i] * a_exp + 4.0 * g_z_yy_yyz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_yz_xz[i] = -2.0 * g_z_yy_z_xz[i] * a_exp + 4.0 * g_z_yy_yyz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_yz_yy[i] = -2.0 * g_z_yy_z_yy[i] * a_exp + 4.0 * g_z_yy_yyz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_yz_yz[i] = -2.0 * g_z_yy_z_yz[i] * a_exp + 4.0 * g_z_yy_yyz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_yz_zz[i] = -2.0 * g_z_yy_z_zz[i] * a_exp + 4.0 * g_z_yy_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1650-1656)

    #pragma omp simd aligned(g_z_0_y_0_0_yy_zz_xx, g_z_0_y_0_0_yy_zz_xy, g_z_0_y_0_0_yy_zz_xz, g_z_0_y_0_0_yy_zz_yy, g_z_0_y_0_0_yy_zz_yz, g_z_0_y_0_0_yy_zz_zz, g_z_yy_yzz_xx, g_z_yy_yzz_xy, g_z_yy_yzz_xz, g_z_yy_yzz_yy, g_z_yy_yzz_yz, g_z_yy_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_yy_zz_xx[i] = 4.0 * g_z_yy_yzz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_zz_xy[i] = 4.0 * g_z_yy_yzz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_zz_xz[i] = 4.0 * g_z_yy_yzz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_zz_yy[i] = 4.0 * g_z_yy_yzz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_zz_yz[i] = 4.0 * g_z_yy_yzz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_zz_zz[i] = 4.0 * g_z_yy_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1656-1662)

    #pragma omp simd aligned(g_z_0_y_0_0_yz_xx_xx, g_z_0_y_0_0_yz_xx_xy, g_z_0_y_0_0_yz_xx_xz, g_z_0_y_0_0_yz_xx_yy, g_z_0_y_0_0_yz_xx_yz, g_z_0_y_0_0_yz_xx_zz, g_z_yz_xxy_xx, g_z_yz_xxy_xy, g_z_yz_xxy_xz, g_z_yz_xxy_yy, g_z_yz_xxy_yz, g_z_yz_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_yz_xx_xx[i] = 4.0 * g_z_yz_xxy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_xx_xy[i] = 4.0 * g_z_yz_xxy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_xx_xz[i] = 4.0 * g_z_yz_xxy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_xx_yy[i] = 4.0 * g_z_yz_xxy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_xx_yz[i] = 4.0 * g_z_yz_xxy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_xx_zz[i] = 4.0 * g_z_yz_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1662-1668)

    #pragma omp simd aligned(g_z_0_y_0_0_yz_xy_xx, g_z_0_y_0_0_yz_xy_xy, g_z_0_y_0_0_yz_xy_xz, g_z_0_y_0_0_yz_xy_yy, g_z_0_y_0_0_yz_xy_yz, g_z_0_y_0_0_yz_xy_zz, g_z_yz_x_xx, g_z_yz_x_xy, g_z_yz_x_xz, g_z_yz_x_yy, g_z_yz_x_yz, g_z_yz_x_zz, g_z_yz_xyy_xx, g_z_yz_xyy_xy, g_z_yz_xyy_xz, g_z_yz_xyy_yy, g_z_yz_xyy_yz, g_z_yz_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_yz_xy_xx[i] = -2.0 * g_z_yz_x_xx[i] * a_exp + 4.0 * g_z_yz_xyy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_xy_xy[i] = -2.0 * g_z_yz_x_xy[i] * a_exp + 4.0 * g_z_yz_xyy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_xy_xz[i] = -2.0 * g_z_yz_x_xz[i] * a_exp + 4.0 * g_z_yz_xyy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_xy_yy[i] = -2.0 * g_z_yz_x_yy[i] * a_exp + 4.0 * g_z_yz_xyy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_xy_yz[i] = -2.0 * g_z_yz_x_yz[i] * a_exp + 4.0 * g_z_yz_xyy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_xy_zz[i] = -2.0 * g_z_yz_x_zz[i] * a_exp + 4.0 * g_z_yz_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1668-1674)

    #pragma omp simd aligned(g_z_0_y_0_0_yz_xz_xx, g_z_0_y_0_0_yz_xz_xy, g_z_0_y_0_0_yz_xz_xz, g_z_0_y_0_0_yz_xz_yy, g_z_0_y_0_0_yz_xz_yz, g_z_0_y_0_0_yz_xz_zz, g_z_yz_xyz_xx, g_z_yz_xyz_xy, g_z_yz_xyz_xz, g_z_yz_xyz_yy, g_z_yz_xyz_yz, g_z_yz_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_yz_xz_xx[i] = 4.0 * g_z_yz_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_xz_xy[i] = 4.0 * g_z_yz_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_xz_xz[i] = 4.0 * g_z_yz_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_xz_yy[i] = 4.0 * g_z_yz_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_xz_yz[i] = 4.0 * g_z_yz_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_xz_zz[i] = 4.0 * g_z_yz_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1674-1680)

    #pragma omp simd aligned(g_z_0_y_0_0_yz_yy_xx, g_z_0_y_0_0_yz_yy_xy, g_z_0_y_0_0_yz_yy_xz, g_z_0_y_0_0_yz_yy_yy, g_z_0_y_0_0_yz_yy_yz, g_z_0_y_0_0_yz_yy_zz, g_z_yz_y_xx, g_z_yz_y_xy, g_z_yz_y_xz, g_z_yz_y_yy, g_z_yz_y_yz, g_z_yz_y_zz, g_z_yz_yyy_xx, g_z_yz_yyy_xy, g_z_yz_yyy_xz, g_z_yz_yyy_yy, g_z_yz_yyy_yz, g_z_yz_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_yz_yy_xx[i] = -4.0 * g_z_yz_y_xx[i] * a_exp + 4.0 * g_z_yz_yyy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_yy_xy[i] = -4.0 * g_z_yz_y_xy[i] * a_exp + 4.0 * g_z_yz_yyy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_yy_xz[i] = -4.0 * g_z_yz_y_xz[i] * a_exp + 4.0 * g_z_yz_yyy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_yy_yy[i] = -4.0 * g_z_yz_y_yy[i] * a_exp + 4.0 * g_z_yz_yyy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_yy_yz[i] = -4.0 * g_z_yz_y_yz[i] * a_exp + 4.0 * g_z_yz_yyy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_yy_zz[i] = -4.0 * g_z_yz_y_zz[i] * a_exp + 4.0 * g_z_yz_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1680-1686)

    #pragma omp simd aligned(g_z_0_y_0_0_yz_yz_xx, g_z_0_y_0_0_yz_yz_xy, g_z_0_y_0_0_yz_yz_xz, g_z_0_y_0_0_yz_yz_yy, g_z_0_y_0_0_yz_yz_yz, g_z_0_y_0_0_yz_yz_zz, g_z_yz_yyz_xx, g_z_yz_yyz_xy, g_z_yz_yyz_xz, g_z_yz_yyz_yy, g_z_yz_yyz_yz, g_z_yz_yyz_zz, g_z_yz_z_xx, g_z_yz_z_xy, g_z_yz_z_xz, g_z_yz_z_yy, g_z_yz_z_yz, g_z_yz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_yz_yz_xx[i] = -2.0 * g_z_yz_z_xx[i] * a_exp + 4.0 * g_z_yz_yyz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_yz_xy[i] = -2.0 * g_z_yz_z_xy[i] * a_exp + 4.0 * g_z_yz_yyz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_yz_xz[i] = -2.0 * g_z_yz_z_xz[i] * a_exp + 4.0 * g_z_yz_yyz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_yz_yy[i] = -2.0 * g_z_yz_z_yy[i] * a_exp + 4.0 * g_z_yz_yyz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_yz_yz[i] = -2.0 * g_z_yz_z_yz[i] * a_exp + 4.0 * g_z_yz_yyz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_yz_zz[i] = -2.0 * g_z_yz_z_zz[i] * a_exp + 4.0 * g_z_yz_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1686-1692)

    #pragma omp simd aligned(g_z_0_y_0_0_yz_zz_xx, g_z_0_y_0_0_yz_zz_xy, g_z_0_y_0_0_yz_zz_xz, g_z_0_y_0_0_yz_zz_yy, g_z_0_y_0_0_yz_zz_yz, g_z_0_y_0_0_yz_zz_zz, g_z_yz_yzz_xx, g_z_yz_yzz_xy, g_z_yz_yzz_xz, g_z_yz_yzz_yy, g_z_yz_yzz_yz, g_z_yz_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_yz_zz_xx[i] = 4.0 * g_z_yz_yzz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_zz_xy[i] = 4.0 * g_z_yz_yzz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_zz_xz[i] = 4.0 * g_z_yz_yzz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_zz_yy[i] = 4.0 * g_z_yz_yzz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_zz_yz[i] = 4.0 * g_z_yz_yzz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_zz_zz[i] = 4.0 * g_z_yz_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1692-1698)

    #pragma omp simd aligned(g_z_0_y_0_0_zz_xx_xx, g_z_0_y_0_0_zz_xx_xy, g_z_0_y_0_0_zz_xx_xz, g_z_0_y_0_0_zz_xx_yy, g_z_0_y_0_0_zz_xx_yz, g_z_0_y_0_0_zz_xx_zz, g_z_zz_xxy_xx, g_z_zz_xxy_xy, g_z_zz_xxy_xz, g_z_zz_xxy_yy, g_z_zz_xxy_yz, g_z_zz_xxy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_zz_xx_xx[i] = 4.0 * g_z_zz_xxy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_xx_xy[i] = 4.0 * g_z_zz_xxy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_xx_xz[i] = 4.0 * g_z_zz_xxy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_xx_yy[i] = 4.0 * g_z_zz_xxy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_xx_yz[i] = 4.0 * g_z_zz_xxy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_xx_zz[i] = 4.0 * g_z_zz_xxy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1698-1704)

    #pragma omp simd aligned(g_z_0_y_0_0_zz_xy_xx, g_z_0_y_0_0_zz_xy_xy, g_z_0_y_0_0_zz_xy_xz, g_z_0_y_0_0_zz_xy_yy, g_z_0_y_0_0_zz_xy_yz, g_z_0_y_0_0_zz_xy_zz, g_z_zz_x_xx, g_z_zz_x_xy, g_z_zz_x_xz, g_z_zz_x_yy, g_z_zz_x_yz, g_z_zz_x_zz, g_z_zz_xyy_xx, g_z_zz_xyy_xy, g_z_zz_xyy_xz, g_z_zz_xyy_yy, g_z_zz_xyy_yz, g_z_zz_xyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_zz_xy_xx[i] = -2.0 * g_z_zz_x_xx[i] * a_exp + 4.0 * g_z_zz_xyy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_xy_xy[i] = -2.0 * g_z_zz_x_xy[i] * a_exp + 4.0 * g_z_zz_xyy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_xy_xz[i] = -2.0 * g_z_zz_x_xz[i] * a_exp + 4.0 * g_z_zz_xyy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_xy_yy[i] = -2.0 * g_z_zz_x_yy[i] * a_exp + 4.0 * g_z_zz_xyy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_xy_yz[i] = -2.0 * g_z_zz_x_yz[i] * a_exp + 4.0 * g_z_zz_xyy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_xy_zz[i] = -2.0 * g_z_zz_x_zz[i] * a_exp + 4.0 * g_z_zz_xyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1704-1710)

    #pragma omp simd aligned(g_z_0_y_0_0_zz_xz_xx, g_z_0_y_0_0_zz_xz_xy, g_z_0_y_0_0_zz_xz_xz, g_z_0_y_0_0_zz_xz_yy, g_z_0_y_0_0_zz_xz_yz, g_z_0_y_0_0_zz_xz_zz, g_z_zz_xyz_xx, g_z_zz_xyz_xy, g_z_zz_xyz_xz, g_z_zz_xyz_yy, g_z_zz_xyz_yz, g_z_zz_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_zz_xz_xx[i] = 4.0 * g_z_zz_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_xz_xy[i] = 4.0 * g_z_zz_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_xz_xz[i] = 4.0 * g_z_zz_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_xz_yy[i] = 4.0 * g_z_zz_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_xz_yz[i] = 4.0 * g_z_zz_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_xz_zz[i] = 4.0 * g_z_zz_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1710-1716)

    #pragma omp simd aligned(g_z_0_y_0_0_zz_yy_xx, g_z_0_y_0_0_zz_yy_xy, g_z_0_y_0_0_zz_yy_xz, g_z_0_y_0_0_zz_yy_yy, g_z_0_y_0_0_zz_yy_yz, g_z_0_y_0_0_zz_yy_zz, g_z_zz_y_xx, g_z_zz_y_xy, g_z_zz_y_xz, g_z_zz_y_yy, g_z_zz_y_yz, g_z_zz_y_zz, g_z_zz_yyy_xx, g_z_zz_yyy_xy, g_z_zz_yyy_xz, g_z_zz_yyy_yy, g_z_zz_yyy_yz, g_z_zz_yyy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_zz_yy_xx[i] = -4.0 * g_z_zz_y_xx[i] * a_exp + 4.0 * g_z_zz_yyy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_yy_xy[i] = -4.0 * g_z_zz_y_xy[i] * a_exp + 4.0 * g_z_zz_yyy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_yy_xz[i] = -4.0 * g_z_zz_y_xz[i] * a_exp + 4.0 * g_z_zz_yyy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_yy_yy[i] = -4.0 * g_z_zz_y_yy[i] * a_exp + 4.0 * g_z_zz_yyy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_yy_yz[i] = -4.0 * g_z_zz_y_yz[i] * a_exp + 4.0 * g_z_zz_yyy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_yy_zz[i] = -4.0 * g_z_zz_y_zz[i] * a_exp + 4.0 * g_z_zz_yyy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1716-1722)

    #pragma omp simd aligned(g_z_0_y_0_0_zz_yz_xx, g_z_0_y_0_0_zz_yz_xy, g_z_0_y_0_0_zz_yz_xz, g_z_0_y_0_0_zz_yz_yy, g_z_0_y_0_0_zz_yz_yz, g_z_0_y_0_0_zz_yz_zz, g_z_zz_yyz_xx, g_z_zz_yyz_xy, g_z_zz_yyz_xz, g_z_zz_yyz_yy, g_z_zz_yyz_yz, g_z_zz_yyz_zz, g_z_zz_z_xx, g_z_zz_z_xy, g_z_zz_z_xz, g_z_zz_z_yy, g_z_zz_z_yz, g_z_zz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_zz_yz_xx[i] = -2.0 * g_z_zz_z_xx[i] * a_exp + 4.0 * g_z_zz_yyz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_yz_xy[i] = -2.0 * g_z_zz_z_xy[i] * a_exp + 4.0 * g_z_zz_yyz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_yz_xz[i] = -2.0 * g_z_zz_z_xz[i] * a_exp + 4.0 * g_z_zz_yyz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_yz_yy[i] = -2.0 * g_z_zz_z_yy[i] * a_exp + 4.0 * g_z_zz_yyz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_yz_yz[i] = -2.0 * g_z_zz_z_yz[i] * a_exp + 4.0 * g_z_zz_yyz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_yz_zz[i] = -2.0 * g_z_zz_z_zz[i] * a_exp + 4.0 * g_z_zz_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1722-1728)

    #pragma omp simd aligned(g_z_0_y_0_0_zz_zz_xx, g_z_0_y_0_0_zz_zz_xy, g_z_0_y_0_0_zz_zz_xz, g_z_0_y_0_0_zz_zz_yy, g_z_0_y_0_0_zz_zz_yz, g_z_0_y_0_0_zz_zz_zz, g_z_zz_yzz_xx, g_z_zz_yzz_xy, g_z_zz_yzz_xz, g_z_zz_yzz_yy, g_z_zz_yzz_yz, g_z_zz_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_zz_zz_xx[i] = 4.0 * g_z_zz_yzz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_zz_xy[i] = 4.0 * g_z_zz_yzz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_zz_xz[i] = 4.0 * g_z_zz_yzz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_zz_yy[i] = 4.0 * g_z_zz_yzz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_zz_yz[i] = 4.0 * g_z_zz_yzz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_zz_zz[i] = 4.0 * g_z_zz_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1728-1734)

    #pragma omp simd aligned(g_z_0_z_0_0_xx_xx_xx, g_z_0_z_0_0_xx_xx_xy, g_z_0_z_0_0_xx_xx_xz, g_z_0_z_0_0_xx_xx_yy, g_z_0_z_0_0_xx_xx_yz, g_z_0_z_0_0_xx_xx_zz, g_z_xx_xxz_xx, g_z_xx_xxz_xy, g_z_xx_xxz_xz, g_z_xx_xxz_yy, g_z_xx_xxz_yz, g_z_xx_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xx_xx_xx[i] = 4.0 * g_z_xx_xxz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_xx_xy[i] = 4.0 * g_z_xx_xxz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_xx_xz[i] = 4.0 * g_z_xx_xxz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_xx_yy[i] = 4.0 * g_z_xx_xxz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_xx_yz[i] = 4.0 * g_z_xx_xxz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_xx_zz[i] = 4.0 * g_z_xx_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1734-1740)

    #pragma omp simd aligned(g_z_0_z_0_0_xx_xy_xx, g_z_0_z_0_0_xx_xy_xy, g_z_0_z_0_0_xx_xy_xz, g_z_0_z_0_0_xx_xy_yy, g_z_0_z_0_0_xx_xy_yz, g_z_0_z_0_0_xx_xy_zz, g_z_xx_xyz_xx, g_z_xx_xyz_xy, g_z_xx_xyz_xz, g_z_xx_xyz_yy, g_z_xx_xyz_yz, g_z_xx_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xx_xy_xx[i] = 4.0 * g_z_xx_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_xy_xy[i] = 4.0 * g_z_xx_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_xy_xz[i] = 4.0 * g_z_xx_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_xy_yy[i] = 4.0 * g_z_xx_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_xy_yz[i] = 4.0 * g_z_xx_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_xy_zz[i] = 4.0 * g_z_xx_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1740-1746)

    #pragma omp simd aligned(g_z_0_z_0_0_xx_xz_xx, g_z_0_z_0_0_xx_xz_xy, g_z_0_z_0_0_xx_xz_xz, g_z_0_z_0_0_xx_xz_yy, g_z_0_z_0_0_xx_xz_yz, g_z_0_z_0_0_xx_xz_zz, g_z_xx_x_xx, g_z_xx_x_xy, g_z_xx_x_xz, g_z_xx_x_yy, g_z_xx_x_yz, g_z_xx_x_zz, g_z_xx_xzz_xx, g_z_xx_xzz_xy, g_z_xx_xzz_xz, g_z_xx_xzz_yy, g_z_xx_xzz_yz, g_z_xx_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xx_xz_xx[i] = -2.0 * g_z_xx_x_xx[i] * a_exp + 4.0 * g_z_xx_xzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_xz_xy[i] = -2.0 * g_z_xx_x_xy[i] * a_exp + 4.0 * g_z_xx_xzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_xz_xz[i] = -2.0 * g_z_xx_x_xz[i] * a_exp + 4.0 * g_z_xx_xzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_xz_yy[i] = -2.0 * g_z_xx_x_yy[i] * a_exp + 4.0 * g_z_xx_xzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_xz_yz[i] = -2.0 * g_z_xx_x_yz[i] * a_exp + 4.0 * g_z_xx_xzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_xz_zz[i] = -2.0 * g_z_xx_x_zz[i] * a_exp + 4.0 * g_z_xx_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1746-1752)

    #pragma omp simd aligned(g_z_0_z_0_0_xx_yy_xx, g_z_0_z_0_0_xx_yy_xy, g_z_0_z_0_0_xx_yy_xz, g_z_0_z_0_0_xx_yy_yy, g_z_0_z_0_0_xx_yy_yz, g_z_0_z_0_0_xx_yy_zz, g_z_xx_yyz_xx, g_z_xx_yyz_xy, g_z_xx_yyz_xz, g_z_xx_yyz_yy, g_z_xx_yyz_yz, g_z_xx_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xx_yy_xx[i] = 4.0 * g_z_xx_yyz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_yy_xy[i] = 4.0 * g_z_xx_yyz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_yy_xz[i] = 4.0 * g_z_xx_yyz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_yy_yy[i] = 4.0 * g_z_xx_yyz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_yy_yz[i] = 4.0 * g_z_xx_yyz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_yy_zz[i] = 4.0 * g_z_xx_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1752-1758)

    #pragma omp simd aligned(g_z_0_z_0_0_xx_yz_xx, g_z_0_z_0_0_xx_yz_xy, g_z_0_z_0_0_xx_yz_xz, g_z_0_z_0_0_xx_yz_yy, g_z_0_z_0_0_xx_yz_yz, g_z_0_z_0_0_xx_yz_zz, g_z_xx_y_xx, g_z_xx_y_xy, g_z_xx_y_xz, g_z_xx_y_yy, g_z_xx_y_yz, g_z_xx_y_zz, g_z_xx_yzz_xx, g_z_xx_yzz_xy, g_z_xx_yzz_xz, g_z_xx_yzz_yy, g_z_xx_yzz_yz, g_z_xx_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xx_yz_xx[i] = -2.0 * g_z_xx_y_xx[i] * a_exp + 4.0 * g_z_xx_yzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_yz_xy[i] = -2.0 * g_z_xx_y_xy[i] * a_exp + 4.0 * g_z_xx_yzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_yz_xz[i] = -2.0 * g_z_xx_y_xz[i] * a_exp + 4.0 * g_z_xx_yzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_yz_yy[i] = -2.0 * g_z_xx_y_yy[i] * a_exp + 4.0 * g_z_xx_yzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_yz_yz[i] = -2.0 * g_z_xx_y_yz[i] * a_exp + 4.0 * g_z_xx_yzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_yz_zz[i] = -2.0 * g_z_xx_y_zz[i] * a_exp + 4.0 * g_z_xx_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1758-1764)

    #pragma omp simd aligned(g_z_0_z_0_0_xx_zz_xx, g_z_0_z_0_0_xx_zz_xy, g_z_0_z_0_0_xx_zz_xz, g_z_0_z_0_0_xx_zz_yy, g_z_0_z_0_0_xx_zz_yz, g_z_0_z_0_0_xx_zz_zz, g_z_xx_z_xx, g_z_xx_z_xy, g_z_xx_z_xz, g_z_xx_z_yy, g_z_xx_z_yz, g_z_xx_z_zz, g_z_xx_zzz_xx, g_z_xx_zzz_xy, g_z_xx_zzz_xz, g_z_xx_zzz_yy, g_z_xx_zzz_yz, g_z_xx_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xx_zz_xx[i] = -4.0 * g_z_xx_z_xx[i] * a_exp + 4.0 * g_z_xx_zzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_zz_xy[i] = -4.0 * g_z_xx_z_xy[i] * a_exp + 4.0 * g_z_xx_zzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_zz_xz[i] = -4.0 * g_z_xx_z_xz[i] * a_exp + 4.0 * g_z_xx_zzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_zz_yy[i] = -4.0 * g_z_xx_z_yy[i] * a_exp + 4.0 * g_z_xx_zzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_zz_yz[i] = -4.0 * g_z_xx_z_yz[i] * a_exp + 4.0 * g_z_xx_zzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_zz_zz[i] = -4.0 * g_z_xx_z_zz[i] * a_exp + 4.0 * g_z_xx_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1764-1770)

    #pragma omp simd aligned(g_z_0_z_0_0_xy_xx_xx, g_z_0_z_0_0_xy_xx_xy, g_z_0_z_0_0_xy_xx_xz, g_z_0_z_0_0_xy_xx_yy, g_z_0_z_0_0_xy_xx_yz, g_z_0_z_0_0_xy_xx_zz, g_z_xy_xxz_xx, g_z_xy_xxz_xy, g_z_xy_xxz_xz, g_z_xy_xxz_yy, g_z_xy_xxz_yz, g_z_xy_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xy_xx_xx[i] = 4.0 * g_z_xy_xxz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_xx_xy[i] = 4.0 * g_z_xy_xxz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_xx_xz[i] = 4.0 * g_z_xy_xxz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_xx_yy[i] = 4.0 * g_z_xy_xxz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_xx_yz[i] = 4.0 * g_z_xy_xxz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_xx_zz[i] = 4.0 * g_z_xy_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1770-1776)

    #pragma omp simd aligned(g_z_0_z_0_0_xy_xy_xx, g_z_0_z_0_0_xy_xy_xy, g_z_0_z_0_0_xy_xy_xz, g_z_0_z_0_0_xy_xy_yy, g_z_0_z_0_0_xy_xy_yz, g_z_0_z_0_0_xy_xy_zz, g_z_xy_xyz_xx, g_z_xy_xyz_xy, g_z_xy_xyz_xz, g_z_xy_xyz_yy, g_z_xy_xyz_yz, g_z_xy_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xy_xy_xx[i] = 4.0 * g_z_xy_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_xy_xy[i] = 4.0 * g_z_xy_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_xy_xz[i] = 4.0 * g_z_xy_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_xy_yy[i] = 4.0 * g_z_xy_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_xy_yz[i] = 4.0 * g_z_xy_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_xy_zz[i] = 4.0 * g_z_xy_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1776-1782)

    #pragma omp simd aligned(g_z_0_z_0_0_xy_xz_xx, g_z_0_z_0_0_xy_xz_xy, g_z_0_z_0_0_xy_xz_xz, g_z_0_z_0_0_xy_xz_yy, g_z_0_z_0_0_xy_xz_yz, g_z_0_z_0_0_xy_xz_zz, g_z_xy_x_xx, g_z_xy_x_xy, g_z_xy_x_xz, g_z_xy_x_yy, g_z_xy_x_yz, g_z_xy_x_zz, g_z_xy_xzz_xx, g_z_xy_xzz_xy, g_z_xy_xzz_xz, g_z_xy_xzz_yy, g_z_xy_xzz_yz, g_z_xy_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xy_xz_xx[i] = -2.0 * g_z_xy_x_xx[i] * a_exp + 4.0 * g_z_xy_xzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_xz_xy[i] = -2.0 * g_z_xy_x_xy[i] * a_exp + 4.0 * g_z_xy_xzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_xz_xz[i] = -2.0 * g_z_xy_x_xz[i] * a_exp + 4.0 * g_z_xy_xzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_xz_yy[i] = -2.0 * g_z_xy_x_yy[i] * a_exp + 4.0 * g_z_xy_xzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_xz_yz[i] = -2.0 * g_z_xy_x_yz[i] * a_exp + 4.0 * g_z_xy_xzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_xz_zz[i] = -2.0 * g_z_xy_x_zz[i] * a_exp + 4.0 * g_z_xy_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1782-1788)

    #pragma omp simd aligned(g_z_0_z_0_0_xy_yy_xx, g_z_0_z_0_0_xy_yy_xy, g_z_0_z_0_0_xy_yy_xz, g_z_0_z_0_0_xy_yy_yy, g_z_0_z_0_0_xy_yy_yz, g_z_0_z_0_0_xy_yy_zz, g_z_xy_yyz_xx, g_z_xy_yyz_xy, g_z_xy_yyz_xz, g_z_xy_yyz_yy, g_z_xy_yyz_yz, g_z_xy_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xy_yy_xx[i] = 4.0 * g_z_xy_yyz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_yy_xy[i] = 4.0 * g_z_xy_yyz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_yy_xz[i] = 4.0 * g_z_xy_yyz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_yy_yy[i] = 4.0 * g_z_xy_yyz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_yy_yz[i] = 4.0 * g_z_xy_yyz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_yy_zz[i] = 4.0 * g_z_xy_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1788-1794)

    #pragma omp simd aligned(g_z_0_z_0_0_xy_yz_xx, g_z_0_z_0_0_xy_yz_xy, g_z_0_z_0_0_xy_yz_xz, g_z_0_z_0_0_xy_yz_yy, g_z_0_z_0_0_xy_yz_yz, g_z_0_z_0_0_xy_yz_zz, g_z_xy_y_xx, g_z_xy_y_xy, g_z_xy_y_xz, g_z_xy_y_yy, g_z_xy_y_yz, g_z_xy_y_zz, g_z_xy_yzz_xx, g_z_xy_yzz_xy, g_z_xy_yzz_xz, g_z_xy_yzz_yy, g_z_xy_yzz_yz, g_z_xy_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xy_yz_xx[i] = -2.0 * g_z_xy_y_xx[i] * a_exp + 4.0 * g_z_xy_yzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_yz_xy[i] = -2.0 * g_z_xy_y_xy[i] * a_exp + 4.0 * g_z_xy_yzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_yz_xz[i] = -2.0 * g_z_xy_y_xz[i] * a_exp + 4.0 * g_z_xy_yzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_yz_yy[i] = -2.0 * g_z_xy_y_yy[i] * a_exp + 4.0 * g_z_xy_yzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_yz_yz[i] = -2.0 * g_z_xy_y_yz[i] * a_exp + 4.0 * g_z_xy_yzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_yz_zz[i] = -2.0 * g_z_xy_y_zz[i] * a_exp + 4.0 * g_z_xy_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1794-1800)

    #pragma omp simd aligned(g_z_0_z_0_0_xy_zz_xx, g_z_0_z_0_0_xy_zz_xy, g_z_0_z_0_0_xy_zz_xz, g_z_0_z_0_0_xy_zz_yy, g_z_0_z_0_0_xy_zz_yz, g_z_0_z_0_0_xy_zz_zz, g_z_xy_z_xx, g_z_xy_z_xy, g_z_xy_z_xz, g_z_xy_z_yy, g_z_xy_z_yz, g_z_xy_z_zz, g_z_xy_zzz_xx, g_z_xy_zzz_xy, g_z_xy_zzz_xz, g_z_xy_zzz_yy, g_z_xy_zzz_yz, g_z_xy_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xy_zz_xx[i] = -4.0 * g_z_xy_z_xx[i] * a_exp + 4.0 * g_z_xy_zzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_zz_xy[i] = -4.0 * g_z_xy_z_xy[i] * a_exp + 4.0 * g_z_xy_zzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_zz_xz[i] = -4.0 * g_z_xy_z_xz[i] * a_exp + 4.0 * g_z_xy_zzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_zz_yy[i] = -4.0 * g_z_xy_z_yy[i] * a_exp + 4.0 * g_z_xy_zzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_zz_yz[i] = -4.0 * g_z_xy_z_yz[i] * a_exp + 4.0 * g_z_xy_zzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_zz_zz[i] = -4.0 * g_z_xy_z_zz[i] * a_exp + 4.0 * g_z_xy_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1800-1806)

    #pragma omp simd aligned(g_z_0_z_0_0_xz_xx_xx, g_z_0_z_0_0_xz_xx_xy, g_z_0_z_0_0_xz_xx_xz, g_z_0_z_0_0_xz_xx_yy, g_z_0_z_0_0_xz_xx_yz, g_z_0_z_0_0_xz_xx_zz, g_z_xz_xxz_xx, g_z_xz_xxz_xy, g_z_xz_xxz_xz, g_z_xz_xxz_yy, g_z_xz_xxz_yz, g_z_xz_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xz_xx_xx[i] = 4.0 * g_z_xz_xxz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_xx_xy[i] = 4.0 * g_z_xz_xxz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_xx_xz[i] = 4.0 * g_z_xz_xxz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_xx_yy[i] = 4.0 * g_z_xz_xxz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_xx_yz[i] = 4.0 * g_z_xz_xxz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_xx_zz[i] = 4.0 * g_z_xz_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1806-1812)

    #pragma omp simd aligned(g_z_0_z_0_0_xz_xy_xx, g_z_0_z_0_0_xz_xy_xy, g_z_0_z_0_0_xz_xy_xz, g_z_0_z_0_0_xz_xy_yy, g_z_0_z_0_0_xz_xy_yz, g_z_0_z_0_0_xz_xy_zz, g_z_xz_xyz_xx, g_z_xz_xyz_xy, g_z_xz_xyz_xz, g_z_xz_xyz_yy, g_z_xz_xyz_yz, g_z_xz_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xz_xy_xx[i] = 4.0 * g_z_xz_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_xy_xy[i] = 4.0 * g_z_xz_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_xy_xz[i] = 4.0 * g_z_xz_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_xy_yy[i] = 4.0 * g_z_xz_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_xy_yz[i] = 4.0 * g_z_xz_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_xy_zz[i] = 4.0 * g_z_xz_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1812-1818)

    #pragma omp simd aligned(g_z_0_z_0_0_xz_xz_xx, g_z_0_z_0_0_xz_xz_xy, g_z_0_z_0_0_xz_xz_xz, g_z_0_z_0_0_xz_xz_yy, g_z_0_z_0_0_xz_xz_yz, g_z_0_z_0_0_xz_xz_zz, g_z_xz_x_xx, g_z_xz_x_xy, g_z_xz_x_xz, g_z_xz_x_yy, g_z_xz_x_yz, g_z_xz_x_zz, g_z_xz_xzz_xx, g_z_xz_xzz_xy, g_z_xz_xzz_xz, g_z_xz_xzz_yy, g_z_xz_xzz_yz, g_z_xz_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xz_xz_xx[i] = -2.0 * g_z_xz_x_xx[i] * a_exp + 4.0 * g_z_xz_xzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_xz_xy[i] = -2.0 * g_z_xz_x_xy[i] * a_exp + 4.0 * g_z_xz_xzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_xz_xz[i] = -2.0 * g_z_xz_x_xz[i] * a_exp + 4.0 * g_z_xz_xzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_xz_yy[i] = -2.0 * g_z_xz_x_yy[i] * a_exp + 4.0 * g_z_xz_xzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_xz_yz[i] = -2.0 * g_z_xz_x_yz[i] * a_exp + 4.0 * g_z_xz_xzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_xz_zz[i] = -2.0 * g_z_xz_x_zz[i] * a_exp + 4.0 * g_z_xz_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1818-1824)

    #pragma omp simd aligned(g_z_0_z_0_0_xz_yy_xx, g_z_0_z_0_0_xz_yy_xy, g_z_0_z_0_0_xz_yy_xz, g_z_0_z_0_0_xz_yy_yy, g_z_0_z_0_0_xz_yy_yz, g_z_0_z_0_0_xz_yy_zz, g_z_xz_yyz_xx, g_z_xz_yyz_xy, g_z_xz_yyz_xz, g_z_xz_yyz_yy, g_z_xz_yyz_yz, g_z_xz_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xz_yy_xx[i] = 4.0 * g_z_xz_yyz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_yy_xy[i] = 4.0 * g_z_xz_yyz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_yy_xz[i] = 4.0 * g_z_xz_yyz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_yy_yy[i] = 4.0 * g_z_xz_yyz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_yy_yz[i] = 4.0 * g_z_xz_yyz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_yy_zz[i] = 4.0 * g_z_xz_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1824-1830)

    #pragma omp simd aligned(g_z_0_z_0_0_xz_yz_xx, g_z_0_z_0_0_xz_yz_xy, g_z_0_z_0_0_xz_yz_xz, g_z_0_z_0_0_xz_yz_yy, g_z_0_z_0_0_xz_yz_yz, g_z_0_z_0_0_xz_yz_zz, g_z_xz_y_xx, g_z_xz_y_xy, g_z_xz_y_xz, g_z_xz_y_yy, g_z_xz_y_yz, g_z_xz_y_zz, g_z_xz_yzz_xx, g_z_xz_yzz_xy, g_z_xz_yzz_xz, g_z_xz_yzz_yy, g_z_xz_yzz_yz, g_z_xz_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xz_yz_xx[i] = -2.0 * g_z_xz_y_xx[i] * a_exp + 4.0 * g_z_xz_yzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_yz_xy[i] = -2.0 * g_z_xz_y_xy[i] * a_exp + 4.0 * g_z_xz_yzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_yz_xz[i] = -2.0 * g_z_xz_y_xz[i] * a_exp + 4.0 * g_z_xz_yzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_yz_yy[i] = -2.0 * g_z_xz_y_yy[i] * a_exp + 4.0 * g_z_xz_yzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_yz_yz[i] = -2.0 * g_z_xz_y_yz[i] * a_exp + 4.0 * g_z_xz_yzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_yz_zz[i] = -2.0 * g_z_xz_y_zz[i] * a_exp + 4.0 * g_z_xz_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1830-1836)

    #pragma omp simd aligned(g_z_0_z_0_0_xz_zz_xx, g_z_0_z_0_0_xz_zz_xy, g_z_0_z_0_0_xz_zz_xz, g_z_0_z_0_0_xz_zz_yy, g_z_0_z_0_0_xz_zz_yz, g_z_0_z_0_0_xz_zz_zz, g_z_xz_z_xx, g_z_xz_z_xy, g_z_xz_z_xz, g_z_xz_z_yy, g_z_xz_z_yz, g_z_xz_z_zz, g_z_xz_zzz_xx, g_z_xz_zzz_xy, g_z_xz_zzz_xz, g_z_xz_zzz_yy, g_z_xz_zzz_yz, g_z_xz_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xz_zz_xx[i] = -4.0 * g_z_xz_z_xx[i] * a_exp + 4.0 * g_z_xz_zzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_zz_xy[i] = -4.0 * g_z_xz_z_xy[i] * a_exp + 4.0 * g_z_xz_zzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_zz_xz[i] = -4.0 * g_z_xz_z_xz[i] * a_exp + 4.0 * g_z_xz_zzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_zz_yy[i] = -4.0 * g_z_xz_z_yy[i] * a_exp + 4.0 * g_z_xz_zzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_zz_yz[i] = -4.0 * g_z_xz_z_yz[i] * a_exp + 4.0 * g_z_xz_zzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_zz_zz[i] = -4.0 * g_z_xz_z_zz[i] * a_exp + 4.0 * g_z_xz_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1836-1842)

    #pragma omp simd aligned(g_z_0_z_0_0_yy_xx_xx, g_z_0_z_0_0_yy_xx_xy, g_z_0_z_0_0_yy_xx_xz, g_z_0_z_0_0_yy_xx_yy, g_z_0_z_0_0_yy_xx_yz, g_z_0_z_0_0_yy_xx_zz, g_z_yy_xxz_xx, g_z_yy_xxz_xy, g_z_yy_xxz_xz, g_z_yy_xxz_yy, g_z_yy_xxz_yz, g_z_yy_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_yy_xx_xx[i] = 4.0 * g_z_yy_xxz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_xx_xy[i] = 4.0 * g_z_yy_xxz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_xx_xz[i] = 4.0 * g_z_yy_xxz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_xx_yy[i] = 4.0 * g_z_yy_xxz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_xx_yz[i] = 4.0 * g_z_yy_xxz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_xx_zz[i] = 4.0 * g_z_yy_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1842-1848)

    #pragma omp simd aligned(g_z_0_z_0_0_yy_xy_xx, g_z_0_z_0_0_yy_xy_xy, g_z_0_z_0_0_yy_xy_xz, g_z_0_z_0_0_yy_xy_yy, g_z_0_z_0_0_yy_xy_yz, g_z_0_z_0_0_yy_xy_zz, g_z_yy_xyz_xx, g_z_yy_xyz_xy, g_z_yy_xyz_xz, g_z_yy_xyz_yy, g_z_yy_xyz_yz, g_z_yy_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_yy_xy_xx[i] = 4.0 * g_z_yy_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_xy_xy[i] = 4.0 * g_z_yy_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_xy_xz[i] = 4.0 * g_z_yy_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_xy_yy[i] = 4.0 * g_z_yy_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_xy_yz[i] = 4.0 * g_z_yy_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_xy_zz[i] = 4.0 * g_z_yy_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1848-1854)

    #pragma omp simd aligned(g_z_0_z_0_0_yy_xz_xx, g_z_0_z_0_0_yy_xz_xy, g_z_0_z_0_0_yy_xz_xz, g_z_0_z_0_0_yy_xz_yy, g_z_0_z_0_0_yy_xz_yz, g_z_0_z_0_0_yy_xz_zz, g_z_yy_x_xx, g_z_yy_x_xy, g_z_yy_x_xz, g_z_yy_x_yy, g_z_yy_x_yz, g_z_yy_x_zz, g_z_yy_xzz_xx, g_z_yy_xzz_xy, g_z_yy_xzz_xz, g_z_yy_xzz_yy, g_z_yy_xzz_yz, g_z_yy_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_yy_xz_xx[i] = -2.0 * g_z_yy_x_xx[i] * a_exp + 4.0 * g_z_yy_xzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_xz_xy[i] = -2.0 * g_z_yy_x_xy[i] * a_exp + 4.0 * g_z_yy_xzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_xz_xz[i] = -2.0 * g_z_yy_x_xz[i] * a_exp + 4.0 * g_z_yy_xzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_xz_yy[i] = -2.0 * g_z_yy_x_yy[i] * a_exp + 4.0 * g_z_yy_xzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_xz_yz[i] = -2.0 * g_z_yy_x_yz[i] * a_exp + 4.0 * g_z_yy_xzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_xz_zz[i] = -2.0 * g_z_yy_x_zz[i] * a_exp + 4.0 * g_z_yy_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1854-1860)

    #pragma omp simd aligned(g_z_0_z_0_0_yy_yy_xx, g_z_0_z_0_0_yy_yy_xy, g_z_0_z_0_0_yy_yy_xz, g_z_0_z_0_0_yy_yy_yy, g_z_0_z_0_0_yy_yy_yz, g_z_0_z_0_0_yy_yy_zz, g_z_yy_yyz_xx, g_z_yy_yyz_xy, g_z_yy_yyz_xz, g_z_yy_yyz_yy, g_z_yy_yyz_yz, g_z_yy_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_yy_yy_xx[i] = 4.0 * g_z_yy_yyz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_yy_xy[i] = 4.0 * g_z_yy_yyz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_yy_xz[i] = 4.0 * g_z_yy_yyz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_yy_yy[i] = 4.0 * g_z_yy_yyz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_yy_yz[i] = 4.0 * g_z_yy_yyz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_yy_zz[i] = 4.0 * g_z_yy_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1860-1866)

    #pragma omp simd aligned(g_z_0_z_0_0_yy_yz_xx, g_z_0_z_0_0_yy_yz_xy, g_z_0_z_0_0_yy_yz_xz, g_z_0_z_0_0_yy_yz_yy, g_z_0_z_0_0_yy_yz_yz, g_z_0_z_0_0_yy_yz_zz, g_z_yy_y_xx, g_z_yy_y_xy, g_z_yy_y_xz, g_z_yy_y_yy, g_z_yy_y_yz, g_z_yy_y_zz, g_z_yy_yzz_xx, g_z_yy_yzz_xy, g_z_yy_yzz_xz, g_z_yy_yzz_yy, g_z_yy_yzz_yz, g_z_yy_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_yy_yz_xx[i] = -2.0 * g_z_yy_y_xx[i] * a_exp + 4.0 * g_z_yy_yzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_yz_xy[i] = -2.0 * g_z_yy_y_xy[i] * a_exp + 4.0 * g_z_yy_yzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_yz_xz[i] = -2.0 * g_z_yy_y_xz[i] * a_exp + 4.0 * g_z_yy_yzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_yz_yy[i] = -2.0 * g_z_yy_y_yy[i] * a_exp + 4.0 * g_z_yy_yzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_yz_yz[i] = -2.0 * g_z_yy_y_yz[i] * a_exp + 4.0 * g_z_yy_yzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_yz_zz[i] = -2.0 * g_z_yy_y_zz[i] * a_exp + 4.0 * g_z_yy_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1866-1872)

    #pragma omp simd aligned(g_z_0_z_0_0_yy_zz_xx, g_z_0_z_0_0_yy_zz_xy, g_z_0_z_0_0_yy_zz_xz, g_z_0_z_0_0_yy_zz_yy, g_z_0_z_0_0_yy_zz_yz, g_z_0_z_0_0_yy_zz_zz, g_z_yy_z_xx, g_z_yy_z_xy, g_z_yy_z_xz, g_z_yy_z_yy, g_z_yy_z_yz, g_z_yy_z_zz, g_z_yy_zzz_xx, g_z_yy_zzz_xy, g_z_yy_zzz_xz, g_z_yy_zzz_yy, g_z_yy_zzz_yz, g_z_yy_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_yy_zz_xx[i] = -4.0 * g_z_yy_z_xx[i] * a_exp + 4.0 * g_z_yy_zzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_zz_xy[i] = -4.0 * g_z_yy_z_xy[i] * a_exp + 4.0 * g_z_yy_zzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_zz_xz[i] = -4.0 * g_z_yy_z_xz[i] * a_exp + 4.0 * g_z_yy_zzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_zz_yy[i] = -4.0 * g_z_yy_z_yy[i] * a_exp + 4.0 * g_z_yy_zzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_zz_yz[i] = -4.0 * g_z_yy_z_yz[i] * a_exp + 4.0 * g_z_yy_zzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_zz_zz[i] = -4.0 * g_z_yy_z_zz[i] * a_exp + 4.0 * g_z_yy_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1872-1878)

    #pragma omp simd aligned(g_z_0_z_0_0_yz_xx_xx, g_z_0_z_0_0_yz_xx_xy, g_z_0_z_0_0_yz_xx_xz, g_z_0_z_0_0_yz_xx_yy, g_z_0_z_0_0_yz_xx_yz, g_z_0_z_0_0_yz_xx_zz, g_z_yz_xxz_xx, g_z_yz_xxz_xy, g_z_yz_xxz_xz, g_z_yz_xxz_yy, g_z_yz_xxz_yz, g_z_yz_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_yz_xx_xx[i] = 4.0 * g_z_yz_xxz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_xx_xy[i] = 4.0 * g_z_yz_xxz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_xx_xz[i] = 4.0 * g_z_yz_xxz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_xx_yy[i] = 4.0 * g_z_yz_xxz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_xx_yz[i] = 4.0 * g_z_yz_xxz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_xx_zz[i] = 4.0 * g_z_yz_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1878-1884)

    #pragma omp simd aligned(g_z_0_z_0_0_yz_xy_xx, g_z_0_z_0_0_yz_xy_xy, g_z_0_z_0_0_yz_xy_xz, g_z_0_z_0_0_yz_xy_yy, g_z_0_z_0_0_yz_xy_yz, g_z_0_z_0_0_yz_xy_zz, g_z_yz_xyz_xx, g_z_yz_xyz_xy, g_z_yz_xyz_xz, g_z_yz_xyz_yy, g_z_yz_xyz_yz, g_z_yz_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_yz_xy_xx[i] = 4.0 * g_z_yz_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_xy_xy[i] = 4.0 * g_z_yz_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_xy_xz[i] = 4.0 * g_z_yz_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_xy_yy[i] = 4.0 * g_z_yz_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_xy_yz[i] = 4.0 * g_z_yz_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_xy_zz[i] = 4.0 * g_z_yz_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1884-1890)

    #pragma omp simd aligned(g_z_0_z_0_0_yz_xz_xx, g_z_0_z_0_0_yz_xz_xy, g_z_0_z_0_0_yz_xz_xz, g_z_0_z_0_0_yz_xz_yy, g_z_0_z_0_0_yz_xz_yz, g_z_0_z_0_0_yz_xz_zz, g_z_yz_x_xx, g_z_yz_x_xy, g_z_yz_x_xz, g_z_yz_x_yy, g_z_yz_x_yz, g_z_yz_x_zz, g_z_yz_xzz_xx, g_z_yz_xzz_xy, g_z_yz_xzz_xz, g_z_yz_xzz_yy, g_z_yz_xzz_yz, g_z_yz_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_yz_xz_xx[i] = -2.0 * g_z_yz_x_xx[i] * a_exp + 4.0 * g_z_yz_xzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_xz_xy[i] = -2.0 * g_z_yz_x_xy[i] * a_exp + 4.0 * g_z_yz_xzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_xz_xz[i] = -2.0 * g_z_yz_x_xz[i] * a_exp + 4.0 * g_z_yz_xzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_xz_yy[i] = -2.0 * g_z_yz_x_yy[i] * a_exp + 4.0 * g_z_yz_xzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_xz_yz[i] = -2.0 * g_z_yz_x_yz[i] * a_exp + 4.0 * g_z_yz_xzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_xz_zz[i] = -2.0 * g_z_yz_x_zz[i] * a_exp + 4.0 * g_z_yz_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1890-1896)

    #pragma omp simd aligned(g_z_0_z_0_0_yz_yy_xx, g_z_0_z_0_0_yz_yy_xy, g_z_0_z_0_0_yz_yy_xz, g_z_0_z_0_0_yz_yy_yy, g_z_0_z_0_0_yz_yy_yz, g_z_0_z_0_0_yz_yy_zz, g_z_yz_yyz_xx, g_z_yz_yyz_xy, g_z_yz_yyz_xz, g_z_yz_yyz_yy, g_z_yz_yyz_yz, g_z_yz_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_yz_yy_xx[i] = 4.0 * g_z_yz_yyz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_yy_xy[i] = 4.0 * g_z_yz_yyz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_yy_xz[i] = 4.0 * g_z_yz_yyz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_yy_yy[i] = 4.0 * g_z_yz_yyz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_yy_yz[i] = 4.0 * g_z_yz_yyz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_yy_zz[i] = 4.0 * g_z_yz_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1896-1902)

    #pragma omp simd aligned(g_z_0_z_0_0_yz_yz_xx, g_z_0_z_0_0_yz_yz_xy, g_z_0_z_0_0_yz_yz_xz, g_z_0_z_0_0_yz_yz_yy, g_z_0_z_0_0_yz_yz_yz, g_z_0_z_0_0_yz_yz_zz, g_z_yz_y_xx, g_z_yz_y_xy, g_z_yz_y_xz, g_z_yz_y_yy, g_z_yz_y_yz, g_z_yz_y_zz, g_z_yz_yzz_xx, g_z_yz_yzz_xy, g_z_yz_yzz_xz, g_z_yz_yzz_yy, g_z_yz_yzz_yz, g_z_yz_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_yz_yz_xx[i] = -2.0 * g_z_yz_y_xx[i] * a_exp + 4.0 * g_z_yz_yzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_yz_xy[i] = -2.0 * g_z_yz_y_xy[i] * a_exp + 4.0 * g_z_yz_yzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_yz_xz[i] = -2.0 * g_z_yz_y_xz[i] * a_exp + 4.0 * g_z_yz_yzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_yz_yy[i] = -2.0 * g_z_yz_y_yy[i] * a_exp + 4.0 * g_z_yz_yzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_yz_yz[i] = -2.0 * g_z_yz_y_yz[i] * a_exp + 4.0 * g_z_yz_yzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_yz_zz[i] = -2.0 * g_z_yz_y_zz[i] * a_exp + 4.0 * g_z_yz_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1902-1908)

    #pragma omp simd aligned(g_z_0_z_0_0_yz_zz_xx, g_z_0_z_0_0_yz_zz_xy, g_z_0_z_0_0_yz_zz_xz, g_z_0_z_0_0_yz_zz_yy, g_z_0_z_0_0_yz_zz_yz, g_z_0_z_0_0_yz_zz_zz, g_z_yz_z_xx, g_z_yz_z_xy, g_z_yz_z_xz, g_z_yz_z_yy, g_z_yz_z_yz, g_z_yz_z_zz, g_z_yz_zzz_xx, g_z_yz_zzz_xy, g_z_yz_zzz_xz, g_z_yz_zzz_yy, g_z_yz_zzz_yz, g_z_yz_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_yz_zz_xx[i] = -4.0 * g_z_yz_z_xx[i] * a_exp + 4.0 * g_z_yz_zzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_zz_xy[i] = -4.0 * g_z_yz_z_xy[i] * a_exp + 4.0 * g_z_yz_zzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_zz_xz[i] = -4.0 * g_z_yz_z_xz[i] * a_exp + 4.0 * g_z_yz_zzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_zz_yy[i] = -4.0 * g_z_yz_z_yy[i] * a_exp + 4.0 * g_z_yz_zzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_zz_yz[i] = -4.0 * g_z_yz_z_yz[i] * a_exp + 4.0 * g_z_yz_zzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_zz_zz[i] = -4.0 * g_z_yz_z_zz[i] * a_exp + 4.0 * g_z_yz_zzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1908-1914)

    #pragma omp simd aligned(g_z_0_z_0_0_zz_xx_xx, g_z_0_z_0_0_zz_xx_xy, g_z_0_z_0_0_zz_xx_xz, g_z_0_z_0_0_zz_xx_yy, g_z_0_z_0_0_zz_xx_yz, g_z_0_z_0_0_zz_xx_zz, g_z_zz_xxz_xx, g_z_zz_xxz_xy, g_z_zz_xxz_xz, g_z_zz_xxz_yy, g_z_zz_xxz_yz, g_z_zz_xxz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_zz_xx_xx[i] = 4.0 * g_z_zz_xxz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_xx_xy[i] = 4.0 * g_z_zz_xxz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_xx_xz[i] = 4.0 * g_z_zz_xxz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_xx_yy[i] = 4.0 * g_z_zz_xxz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_xx_yz[i] = 4.0 * g_z_zz_xxz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_xx_zz[i] = 4.0 * g_z_zz_xxz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1914-1920)

    #pragma omp simd aligned(g_z_0_z_0_0_zz_xy_xx, g_z_0_z_0_0_zz_xy_xy, g_z_0_z_0_0_zz_xy_xz, g_z_0_z_0_0_zz_xy_yy, g_z_0_z_0_0_zz_xy_yz, g_z_0_z_0_0_zz_xy_zz, g_z_zz_xyz_xx, g_z_zz_xyz_xy, g_z_zz_xyz_xz, g_z_zz_xyz_yy, g_z_zz_xyz_yz, g_z_zz_xyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_zz_xy_xx[i] = 4.0 * g_z_zz_xyz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_xy_xy[i] = 4.0 * g_z_zz_xyz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_xy_xz[i] = 4.0 * g_z_zz_xyz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_xy_yy[i] = 4.0 * g_z_zz_xyz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_xy_yz[i] = 4.0 * g_z_zz_xyz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_xy_zz[i] = 4.0 * g_z_zz_xyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1920-1926)

    #pragma omp simd aligned(g_z_0_z_0_0_zz_xz_xx, g_z_0_z_0_0_zz_xz_xy, g_z_0_z_0_0_zz_xz_xz, g_z_0_z_0_0_zz_xz_yy, g_z_0_z_0_0_zz_xz_yz, g_z_0_z_0_0_zz_xz_zz, g_z_zz_x_xx, g_z_zz_x_xy, g_z_zz_x_xz, g_z_zz_x_yy, g_z_zz_x_yz, g_z_zz_x_zz, g_z_zz_xzz_xx, g_z_zz_xzz_xy, g_z_zz_xzz_xz, g_z_zz_xzz_yy, g_z_zz_xzz_yz, g_z_zz_xzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_zz_xz_xx[i] = -2.0 * g_z_zz_x_xx[i] * a_exp + 4.0 * g_z_zz_xzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_xz_xy[i] = -2.0 * g_z_zz_x_xy[i] * a_exp + 4.0 * g_z_zz_xzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_xz_xz[i] = -2.0 * g_z_zz_x_xz[i] * a_exp + 4.0 * g_z_zz_xzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_xz_yy[i] = -2.0 * g_z_zz_x_yy[i] * a_exp + 4.0 * g_z_zz_xzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_xz_yz[i] = -2.0 * g_z_zz_x_yz[i] * a_exp + 4.0 * g_z_zz_xzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_xz_zz[i] = -2.0 * g_z_zz_x_zz[i] * a_exp + 4.0 * g_z_zz_xzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1926-1932)

    #pragma omp simd aligned(g_z_0_z_0_0_zz_yy_xx, g_z_0_z_0_0_zz_yy_xy, g_z_0_z_0_0_zz_yy_xz, g_z_0_z_0_0_zz_yy_yy, g_z_0_z_0_0_zz_yy_yz, g_z_0_z_0_0_zz_yy_zz, g_z_zz_yyz_xx, g_z_zz_yyz_xy, g_z_zz_yyz_xz, g_z_zz_yyz_yy, g_z_zz_yyz_yz, g_z_zz_yyz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_zz_yy_xx[i] = 4.0 * g_z_zz_yyz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_yy_xy[i] = 4.0 * g_z_zz_yyz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_yy_xz[i] = 4.0 * g_z_zz_yyz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_yy_yy[i] = 4.0 * g_z_zz_yyz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_yy_yz[i] = 4.0 * g_z_zz_yyz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_yy_zz[i] = 4.0 * g_z_zz_yyz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1932-1938)

    #pragma omp simd aligned(g_z_0_z_0_0_zz_yz_xx, g_z_0_z_0_0_zz_yz_xy, g_z_0_z_0_0_zz_yz_xz, g_z_0_z_0_0_zz_yz_yy, g_z_0_z_0_0_zz_yz_yz, g_z_0_z_0_0_zz_yz_zz, g_z_zz_y_xx, g_z_zz_y_xy, g_z_zz_y_xz, g_z_zz_y_yy, g_z_zz_y_yz, g_z_zz_y_zz, g_z_zz_yzz_xx, g_z_zz_yzz_xy, g_z_zz_yzz_xz, g_z_zz_yzz_yy, g_z_zz_yzz_yz, g_z_zz_yzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_zz_yz_xx[i] = -2.0 * g_z_zz_y_xx[i] * a_exp + 4.0 * g_z_zz_yzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_yz_xy[i] = -2.0 * g_z_zz_y_xy[i] * a_exp + 4.0 * g_z_zz_yzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_yz_xz[i] = -2.0 * g_z_zz_y_xz[i] * a_exp + 4.0 * g_z_zz_yzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_yz_yy[i] = -2.0 * g_z_zz_y_yy[i] * a_exp + 4.0 * g_z_zz_yzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_yz_yz[i] = -2.0 * g_z_zz_y_yz[i] * a_exp + 4.0 * g_z_zz_yzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_yz_zz[i] = -2.0 * g_z_zz_y_zz[i] * a_exp + 4.0 * g_z_zz_yzz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1938-1944)

    #pragma omp simd aligned(g_z_0_z_0_0_zz_zz_xx, g_z_0_z_0_0_zz_zz_xy, g_z_0_z_0_0_zz_zz_xz, g_z_0_z_0_0_zz_zz_yy, g_z_0_z_0_0_zz_zz_yz, g_z_0_z_0_0_zz_zz_zz, g_z_zz_z_xx, g_z_zz_z_xy, g_z_zz_z_xz, g_z_zz_z_yy, g_z_zz_z_yz, g_z_zz_z_zz, g_z_zz_zzz_xx, g_z_zz_zzz_xy, g_z_zz_zzz_xz, g_z_zz_zzz_yy, g_z_zz_zzz_yz, g_z_zz_zzz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_zz_zz_xx[i] = -4.0 * g_z_zz_z_xx[i] * a_exp + 4.0 * g_z_zz_zzz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_zz_xy[i] = -4.0 * g_z_zz_z_xy[i] * a_exp + 4.0 * g_z_zz_zzz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_zz_xz[i] = -4.0 * g_z_zz_z_xz[i] * a_exp + 4.0 * g_z_zz_zzz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_zz_yy[i] = -4.0 * g_z_zz_z_yy[i] * a_exp + 4.0 * g_z_zz_zzz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_zz_yz[i] = -4.0 * g_z_zz_z_yz[i] * a_exp + 4.0 * g_z_zz_zzz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_zz_zz[i] = -4.0 * g_z_zz_z_zz[i] * a_exp + 4.0 * g_z_zz_zzz_zz[i] * a_exp * c_exps[i];
    }
}

} // t4c_geom namespace


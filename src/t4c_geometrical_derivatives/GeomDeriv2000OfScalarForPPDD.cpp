#include "GeomDeriv2000OfScalarForPPDD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom2000_ppdd_0(CSimdArray<double>& buffer_2000_ppdd,
                     const CSimdArray<double>& buffer_ppdd,
                     const CSimdArray<double>& buffer_fpdd,
                     const double a_exp) -> void
{
    const auto ndims = buffer_2000_ppdd.number_of_columns();

    /// Set up components of auxilary buffer : buffer_ppdd

    auto g_x_x_xx_xx = buffer_ppdd[0];

    auto g_x_x_xx_xy = buffer_ppdd[1];

    auto g_x_x_xx_xz = buffer_ppdd[2];

    auto g_x_x_xx_yy = buffer_ppdd[3];

    auto g_x_x_xx_yz = buffer_ppdd[4];

    auto g_x_x_xx_zz = buffer_ppdd[5];

    auto g_x_x_xy_xx = buffer_ppdd[6];

    auto g_x_x_xy_xy = buffer_ppdd[7];

    auto g_x_x_xy_xz = buffer_ppdd[8];

    auto g_x_x_xy_yy = buffer_ppdd[9];

    auto g_x_x_xy_yz = buffer_ppdd[10];

    auto g_x_x_xy_zz = buffer_ppdd[11];

    auto g_x_x_xz_xx = buffer_ppdd[12];

    auto g_x_x_xz_xy = buffer_ppdd[13];

    auto g_x_x_xz_xz = buffer_ppdd[14];

    auto g_x_x_xz_yy = buffer_ppdd[15];

    auto g_x_x_xz_yz = buffer_ppdd[16];

    auto g_x_x_xz_zz = buffer_ppdd[17];

    auto g_x_x_yy_xx = buffer_ppdd[18];

    auto g_x_x_yy_xy = buffer_ppdd[19];

    auto g_x_x_yy_xz = buffer_ppdd[20];

    auto g_x_x_yy_yy = buffer_ppdd[21];

    auto g_x_x_yy_yz = buffer_ppdd[22];

    auto g_x_x_yy_zz = buffer_ppdd[23];

    auto g_x_x_yz_xx = buffer_ppdd[24];

    auto g_x_x_yz_xy = buffer_ppdd[25];

    auto g_x_x_yz_xz = buffer_ppdd[26];

    auto g_x_x_yz_yy = buffer_ppdd[27];

    auto g_x_x_yz_yz = buffer_ppdd[28];

    auto g_x_x_yz_zz = buffer_ppdd[29];

    auto g_x_x_zz_xx = buffer_ppdd[30];

    auto g_x_x_zz_xy = buffer_ppdd[31];

    auto g_x_x_zz_xz = buffer_ppdd[32];

    auto g_x_x_zz_yy = buffer_ppdd[33];

    auto g_x_x_zz_yz = buffer_ppdd[34];

    auto g_x_x_zz_zz = buffer_ppdd[35];

    auto g_x_y_xx_xx = buffer_ppdd[36];

    auto g_x_y_xx_xy = buffer_ppdd[37];

    auto g_x_y_xx_xz = buffer_ppdd[38];

    auto g_x_y_xx_yy = buffer_ppdd[39];

    auto g_x_y_xx_yz = buffer_ppdd[40];

    auto g_x_y_xx_zz = buffer_ppdd[41];

    auto g_x_y_xy_xx = buffer_ppdd[42];

    auto g_x_y_xy_xy = buffer_ppdd[43];

    auto g_x_y_xy_xz = buffer_ppdd[44];

    auto g_x_y_xy_yy = buffer_ppdd[45];

    auto g_x_y_xy_yz = buffer_ppdd[46];

    auto g_x_y_xy_zz = buffer_ppdd[47];

    auto g_x_y_xz_xx = buffer_ppdd[48];

    auto g_x_y_xz_xy = buffer_ppdd[49];

    auto g_x_y_xz_xz = buffer_ppdd[50];

    auto g_x_y_xz_yy = buffer_ppdd[51];

    auto g_x_y_xz_yz = buffer_ppdd[52];

    auto g_x_y_xz_zz = buffer_ppdd[53];

    auto g_x_y_yy_xx = buffer_ppdd[54];

    auto g_x_y_yy_xy = buffer_ppdd[55];

    auto g_x_y_yy_xz = buffer_ppdd[56];

    auto g_x_y_yy_yy = buffer_ppdd[57];

    auto g_x_y_yy_yz = buffer_ppdd[58];

    auto g_x_y_yy_zz = buffer_ppdd[59];

    auto g_x_y_yz_xx = buffer_ppdd[60];

    auto g_x_y_yz_xy = buffer_ppdd[61];

    auto g_x_y_yz_xz = buffer_ppdd[62];

    auto g_x_y_yz_yy = buffer_ppdd[63];

    auto g_x_y_yz_yz = buffer_ppdd[64];

    auto g_x_y_yz_zz = buffer_ppdd[65];

    auto g_x_y_zz_xx = buffer_ppdd[66];

    auto g_x_y_zz_xy = buffer_ppdd[67];

    auto g_x_y_zz_xz = buffer_ppdd[68];

    auto g_x_y_zz_yy = buffer_ppdd[69];

    auto g_x_y_zz_yz = buffer_ppdd[70];

    auto g_x_y_zz_zz = buffer_ppdd[71];

    auto g_x_z_xx_xx = buffer_ppdd[72];

    auto g_x_z_xx_xy = buffer_ppdd[73];

    auto g_x_z_xx_xz = buffer_ppdd[74];

    auto g_x_z_xx_yy = buffer_ppdd[75];

    auto g_x_z_xx_yz = buffer_ppdd[76];

    auto g_x_z_xx_zz = buffer_ppdd[77];

    auto g_x_z_xy_xx = buffer_ppdd[78];

    auto g_x_z_xy_xy = buffer_ppdd[79];

    auto g_x_z_xy_xz = buffer_ppdd[80];

    auto g_x_z_xy_yy = buffer_ppdd[81];

    auto g_x_z_xy_yz = buffer_ppdd[82];

    auto g_x_z_xy_zz = buffer_ppdd[83];

    auto g_x_z_xz_xx = buffer_ppdd[84];

    auto g_x_z_xz_xy = buffer_ppdd[85];

    auto g_x_z_xz_xz = buffer_ppdd[86];

    auto g_x_z_xz_yy = buffer_ppdd[87];

    auto g_x_z_xz_yz = buffer_ppdd[88];

    auto g_x_z_xz_zz = buffer_ppdd[89];

    auto g_x_z_yy_xx = buffer_ppdd[90];

    auto g_x_z_yy_xy = buffer_ppdd[91];

    auto g_x_z_yy_xz = buffer_ppdd[92];

    auto g_x_z_yy_yy = buffer_ppdd[93];

    auto g_x_z_yy_yz = buffer_ppdd[94];

    auto g_x_z_yy_zz = buffer_ppdd[95];

    auto g_x_z_yz_xx = buffer_ppdd[96];

    auto g_x_z_yz_xy = buffer_ppdd[97];

    auto g_x_z_yz_xz = buffer_ppdd[98];

    auto g_x_z_yz_yy = buffer_ppdd[99];

    auto g_x_z_yz_yz = buffer_ppdd[100];

    auto g_x_z_yz_zz = buffer_ppdd[101];

    auto g_x_z_zz_xx = buffer_ppdd[102];

    auto g_x_z_zz_xy = buffer_ppdd[103];

    auto g_x_z_zz_xz = buffer_ppdd[104];

    auto g_x_z_zz_yy = buffer_ppdd[105];

    auto g_x_z_zz_yz = buffer_ppdd[106];

    auto g_x_z_zz_zz = buffer_ppdd[107];

    auto g_y_x_xx_xx = buffer_ppdd[108];

    auto g_y_x_xx_xy = buffer_ppdd[109];

    auto g_y_x_xx_xz = buffer_ppdd[110];

    auto g_y_x_xx_yy = buffer_ppdd[111];

    auto g_y_x_xx_yz = buffer_ppdd[112];

    auto g_y_x_xx_zz = buffer_ppdd[113];

    auto g_y_x_xy_xx = buffer_ppdd[114];

    auto g_y_x_xy_xy = buffer_ppdd[115];

    auto g_y_x_xy_xz = buffer_ppdd[116];

    auto g_y_x_xy_yy = buffer_ppdd[117];

    auto g_y_x_xy_yz = buffer_ppdd[118];

    auto g_y_x_xy_zz = buffer_ppdd[119];

    auto g_y_x_xz_xx = buffer_ppdd[120];

    auto g_y_x_xz_xy = buffer_ppdd[121];

    auto g_y_x_xz_xz = buffer_ppdd[122];

    auto g_y_x_xz_yy = buffer_ppdd[123];

    auto g_y_x_xz_yz = buffer_ppdd[124];

    auto g_y_x_xz_zz = buffer_ppdd[125];

    auto g_y_x_yy_xx = buffer_ppdd[126];

    auto g_y_x_yy_xy = buffer_ppdd[127];

    auto g_y_x_yy_xz = buffer_ppdd[128];

    auto g_y_x_yy_yy = buffer_ppdd[129];

    auto g_y_x_yy_yz = buffer_ppdd[130];

    auto g_y_x_yy_zz = buffer_ppdd[131];

    auto g_y_x_yz_xx = buffer_ppdd[132];

    auto g_y_x_yz_xy = buffer_ppdd[133];

    auto g_y_x_yz_xz = buffer_ppdd[134];

    auto g_y_x_yz_yy = buffer_ppdd[135];

    auto g_y_x_yz_yz = buffer_ppdd[136];

    auto g_y_x_yz_zz = buffer_ppdd[137];

    auto g_y_x_zz_xx = buffer_ppdd[138];

    auto g_y_x_zz_xy = buffer_ppdd[139];

    auto g_y_x_zz_xz = buffer_ppdd[140];

    auto g_y_x_zz_yy = buffer_ppdd[141];

    auto g_y_x_zz_yz = buffer_ppdd[142];

    auto g_y_x_zz_zz = buffer_ppdd[143];

    auto g_y_y_xx_xx = buffer_ppdd[144];

    auto g_y_y_xx_xy = buffer_ppdd[145];

    auto g_y_y_xx_xz = buffer_ppdd[146];

    auto g_y_y_xx_yy = buffer_ppdd[147];

    auto g_y_y_xx_yz = buffer_ppdd[148];

    auto g_y_y_xx_zz = buffer_ppdd[149];

    auto g_y_y_xy_xx = buffer_ppdd[150];

    auto g_y_y_xy_xy = buffer_ppdd[151];

    auto g_y_y_xy_xz = buffer_ppdd[152];

    auto g_y_y_xy_yy = buffer_ppdd[153];

    auto g_y_y_xy_yz = buffer_ppdd[154];

    auto g_y_y_xy_zz = buffer_ppdd[155];

    auto g_y_y_xz_xx = buffer_ppdd[156];

    auto g_y_y_xz_xy = buffer_ppdd[157];

    auto g_y_y_xz_xz = buffer_ppdd[158];

    auto g_y_y_xz_yy = buffer_ppdd[159];

    auto g_y_y_xz_yz = buffer_ppdd[160];

    auto g_y_y_xz_zz = buffer_ppdd[161];

    auto g_y_y_yy_xx = buffer_ppdd[162];

    auto g_y_y_yy_xy = buffer_ppdd[163];

    auto g_y_y_yy_xz = buffer_ppdd[164];

    auto g_y_y_yy_yy = buffer_ppdd[165];

    auto g_y_y_yy_yz = buffer_ppdd[166];

    auto g_y_y_yy_zz = buffer_ppdd[167];

    auto g_y_y_yz_xx = buffer_ppdd[168];

    auto g_y_y_yz_xy = buffer_ppdd[169];

    auto g_y_y_yz_xz = buffer_ppdd[170];

    auto g_y_y_yz_yy = buffer_ppdd[171];

    auto g_y_y_yz_yz = buffer_ppdd[172];

    auto g_y_y_yz_zz = buffer_ppdd[173];

    auto g_y_y_zz_xx = buffer_ppdd[174];

    auto g_y_y_zz_xy = buffer_ppdd[175];

    auto g_y_y_zz_xz = buffer_ppdd[176];

    auto g_y_y_zz_yy = buffer_ppdd[177];

    auto g_y_y_zz_yz = buffer_ppdd[178];

    auto g_y_y_zz_zz = buffer_ppdd[179];

    auto g_y_z_xx_xx = buffer_ppdd[180];

    auto g_y_z_xx_xy = buffer_ppdd[181];

    auto g_y_z_xx_xz = buffer_ppdd[182];

    auto g_y_z_xx_yy = buffer_ppdd[183];

    auto g_y_z_xx_yz = buffer_ppdd[184];

    auto g_y_z_xx_zz = buffer_ppdd[185];

    auto g_y_z_xy_xx = buffer_ppdd[186];

    auto g_y_z_xy_xy = buffer_ppdd[187];

    auto g_y_z_xy_xz = buffer_ppdd[188];

    auto g_y_z_xy_yy = buffer_ppdd[189];

    auto g_y_z_xy_yz = buffer_ppdd[190];

    auto g_y_z_xy_zz = buffer_ppdd[191];

    auto g_y_z_xz_xx = buffer_ppdd[192];

    auto g_y_z_xz_xy = buffer_ppdd[193];

    auto g_y_z_xz_xz = buffer_ppdd[194];

    auto g_y_z_xz_yy = buffer_ppdd[195];

    auto g_y_z_xz_yz = buffer_ppdd[196];

    auto g_y_z_xz_zz = buffer_ppdd[197];

    auto g_y_z_yy_xx = buffer_ppdd[198];

    auto g_y_z_yy_xy = buffer_ppdd[199];

    auto g_y_z_yy_xz = buffer_ppdd[200];

    auto g_y_z_yy_yy = buffer_ppdd[201];

    auto g_y_z_yy_yz = buffer_ppdd[202];

    auto g_y_z_yy_zz = buffer_ppdd[203];

    auto g_y_z_yz_xx = buffer_ppdd[204];

    auto g_y_z_yz_xy = buffer_ppdd[205];

    auto g_y_z_yz_xz = buffer_ppdd[206];

    auto g_y_z_yz_yy = buffer_ppdd[207];

    auto g_y_z_yz_yz = buffer_ppdd[208];

    auto g_y_z_yz_zz = buffer_ppdd[209];

    auto g_y_z_zz_xx = buffer_ppdd[210];

    auto g_y_z_zz_xy = buffer_ppdd[211];

    auto g_y_z_zz_xz = buffer_ppdd[212];

    auto g_y_z_zz_yy = buffer_ppdd[213];

    auto g_y_z_zz_yz = buffer_ppdd[214];

    auto g_y_z_zz_zz = buffer_ppdd[215];

    auto g_z_x_xx_xx = buffer_ppdd[216];

    auto g_z_x_xx_xy = buffer_ppdd[217];

    auto g_z_x_xx_xz = buffer_ppdd[218];

    auto g_z_x_xx_yy = buffer_ppdd[219];

    auto g_z_x_xx_yz = buffer_ppdd[220];

    auto g_z_x_xx_zz = buffer_ppdd[221];

    auto g_z_x_xy_xx = buffer_ppdd[222];

    auto g_z_x_xy_xy = buffer_ppdd[223];

    auto g_z_x_xy_xz = buffer_ppdd[224];

    auto g_z_x_xy_yy = buffer_ppdd[225];

    auto g_z_x_xy_yz = buffer_ppdd[226];

    auto g_z_x_xy_zz = buffer_ppdd[227];

    auto g_z_x_xz_xx = buffer_ppdd[228];

    auto g_z_x_xz_xy = buffer_ppdd[229];

    auto g_z_x_xz_xz = buffer_ppdd[230];

    auto g_z_x_xz_yy = buffer_ppdd[231];

    auto g_z_x_xz_yz = buffer_ppdd[232];

    auto g_z_x_xz_zz = buffer_ppdd[233];

    auto g_z_x_yy_xx = buffer_ppdd[234];

    auto g_z_x_yy_xy = buffer_ppdd[235];

    auto g_z_x_yy_xz = buffer_ppdd[236];

    auto g_z_x_yy_yy = buffer_ppdd[237];

    auto g_z_x_yy_yz = buffer_ppdd[238];

    auto g_z_x_yy_zz = buffer_ppdd[239];

    auto g_z_x_yz_xx = buffer_ppdd[240];

    auto g_z_x_yz_xy = buffer_ppdd[241];

    auto g_z_x_yz_xz = buffer_ppdd[242];

    auto g_z_x_yz_yy = buffer_ppdd[243];

    auto g_z_x_yz_yz = buffer_ppdd[244];

    auto g_z_x_yz_zz = buffer_ppdd[245];

    auto g_z_x_zz_xx = buffer_ppdd[246];

    auto g_z_x_zz_xy = buffer_ppdd[247];

    auto g_z_x_zz_xz = buffer_ppdd[248];

    auto g_z_x_zz_yy = buffer_ppdd[249];

    auto g_z_x_zz_yz = buffer_ppdd[250];

    auto g_z_x_zz_zz = buffer_ppdd[251];

    auto g_z_y_xx_xx = buffer_ppdd[252];

    auto g_z_y_xx_xy = buffer_ppdd[253];

    auto g_z_y_xx_xz = buffer_ppdd[254];

    auto g_z_y_xx_yy = buffer_ppdd[255];

    auto g_z_y_xx_yz = buffer_ppdd[256];

    auto g_z_y_xx_zz = buffer_ppdd[257];

    auto g_z_y_xy_xx = buffer_ppdd[258];

    auto g_z_y_xy_xy = buffer_ppdd[259];

    auto g_z_y_xy_xz = buffer_ppdd[260];

    auto g_z_y_xy_yy = buffer_ppdd[261];

    auto g_z_y_xy_yz = buffer_ppdd[262];

    auto g_z_y_xy_zz = buffer_ppdd[263];

    auto g_z_y_xz_xx = buffer_ppdd[264];

    auto g_z_y_xz_xy = buffer_ppdd[265];

    auto g_z_y_xz_xz = buffer_ppdd[266];

    auto g_z_y_xz_yy = buffer_ppdd[267];

    auto g_z_y_xz_yz = buffer_ppdd[268];

    auto g_z_y_xz_zz = buffer_ppdd[269];

    auto g_z_y_yy_xx = buffer_ppdd[270];

    auto g_z_y_yy_xy = buffer_ppdd[271];

    auto g_z_y_yy_xz = buffer_ppdd[272];

    auto g_z_y_yy_yy = buffer_ppdd[273];

    auto g_z_y_yy_yz = buffer_ppdd[274];

    auto g_z_y_yy_zz = buffer_ppdd[275];

    auto g_z_y_yz_xx = buffer_ppdd[276];

    auto g_z_y_yz_xy = buffer_ppdd[277];

    auto g_z_y_yz_xz = buffer_ppdd[278];

    auto g_z_y_yz_yy = buffer_ppdd[279];

    auto g_z_y_yz_yz = buffer_ppdd[280];

    auto g_z_y_yz_zz = buffer_ppdd[281];

    auto g_z_y_zz_xx = buffer_ppdd[282];

    auto g_z_y_zz_xy = buffer_ppdd[283];

    auto g_z_y_zz_xz = buffer_ppdd[284];

    auto g_z_y_zz_yy = buffer_ppdd[285];

    auto g_z_y_zz_yz = buffer_ppdd[286];

    auto g_z_y_zz_zz = buffer_ppdd[287];

    auto g_z_z_xx_xx = buffer_ppdd[288];

    auto g_z_z_xx_xy = buffer_ppdd[289];

    auto g_z_z_xx_xz = buffer_ppdd[290];

    auto g_z_z_xx_yy = buffer_ppdd[291];

    auto g_z_z_xx_yz = buffer_ppdd[292];

    auto g_z_z_xx_zz = buffer_ppdd[293];

    auto g_z_z_xy_xx = buffer_ppdd[294];

    auto g_z_z_xy_xy = buffer_ppdd[295];

    auto g_z_z_xy_xz = buffer_ppdd[296];

    auto g_z_z_xy_yy = buffer_ppdd[297];

    auto g_z_z_xy_yz = buffer_ppdd[298];

    auto g_z_z_xy_zz = buffer_ppdd[299];

    auto g_z_z_xz_xx = buffer_ppdd[300];

    auto g_z_z_xz_xy = buffer_ppdd[301];

    auto g_z_z_xz_xz = buffer_ppdd[302];

    auto g_z_z_xz_yy = buffer_ppdd[303];

    auto g_z_z_xz_yz = buffer_ppdd[304];

    auto g_z_z_xz_zz = buffer_ppdd[305];

    auto g_z_z_yy_xx = buffer_ppdd[306];

    auto g_z_z_yy_xy = buffer_ppdd[307];

    auto g_z_z_yy_xz = buffer_ppdd[308];

    auto g_z_z_yy_yy = buffer_ppdd[309];

    auto g_z_z_yy_yz = buffer_ppdd[310];

    auto g_z_z_yy_zz = buffer_ppdd[311];

    auto g_z_z_yz_xx = buffer_ppdd[312];

    auto g_z_z_yz_xy = buffer_ppdd[313];

    auto g_z_z_yz_xz = buffer_ppdd[314];

    auto g_z_z_yz_yy = buffer_ppdd[315];

    auto g_z_z_yz_yz = buffer_ppdd[316];

    auto g_z_z_yz_zz = buffer_ppdd[317];

    auto g_z_z_zz_xx = buffer_ppdd[318];

    auto g_z_z_zz_xy = buffer_ppdd[319];

    auto g_z_z_zz_xz = buffer_ppdd[320];

    auto g_z_z_zz_yy = buffer_ppdd[321];

    auto g_z_z_zz_yz = buffer_ppdd[322];

    auto g_z_z_zz_zz = buffer_ppdd[323];

    /// Set up components of auxilary buffer : buffer_fpdd

    auto g_xxx_x_xx_xx = buffer_fpdd[0];

    auto g_xxx_x_xx_xy = buffer_fpdd[1];

    auto g_xxx_x_xx_xz = buffer_fpdd[2];

    auto g_xxx_x_xx_yy = buffer_fpdd[3];

    auto g_xxx_x_xx_yz = buffer_fpdd[4];

    auto g_xxx_x_xx_zz = buffer_fpdd[5];

    auto g_xxx_x_xy_xx = buffer_fpdd[6];

    auto g_xxx_x_xy_xy = buffer_fpdd[7];

    auto g_xxx_x_xy_xz = buffer_fpdd[8];

    auto g_xxx_x_xy_yy = buffer_fpdd[9];

    auto g_xxx_x_xy_yz = buffer_fpdd[10];

    auto g_xxx_x_xy_zz = buffer_fpdd[11];

    auto g_xxx_x_xz_xx = buffer_fpdd[12];

    auto g_xxx_x_xz_xy = buffer_fpdd[13];

    auto g_xxx_x_xz_xz = buffer_fpdd[14];

    auto g_xxx_x_xz_yy = buffer_fpdd[15];

    auto g_xxx_x_xz_yz = buffer_fpdd[16];

    auto g_xxx_x_xz_zz = buffer_fpdd[17];

    auto g_xxx_x_yy_xx = buffer_fpdd[18];

    auto g_xxx_x_yy_xy = buffer_fpdd[19];

    auto g_xxx_x_yy_xz = buffer_fpdd[20];

    auto g_xxx_x_yy_yy = buffer_fpdd[21];

    auto g_xxx_x_yy_yz = buffer_fpdd[22];

    auto g_xxx_x_yy_zz = buffer_fpdd[23];

    auto g_xxx_x_yz_xx = buffer_fpdd[24];

    auto g_xxx_x_yz_xy = buffer_fpdd[25];

    auto g_xxx_x_yz_xz = buffer_fpdd[26];

    auto g_xxx_x_yz_yy = buffer_fpdd[27];

    auto g_xxx_x_yz_yz = buffer_fpdd[28];

    auto g_xxx_x_yz_zz = buffer_fpdd[29];

    auto g_xxx_x_zz_xx = buffer_fpdd[30];

    auto g_xxx_x_zz_xy = buffer_fpdd[31];

    auto g_xxx_x_zz_xz = buffer_fpdd[32];

    auto g_xxx_x_zz_yy = buffer_fpdd[33];

    auto g_xxx_x_zz_yz = buffer_fpdd[34];

    auto g_xxx_x_zz_zz = buffer_fpdd[35];

    auto g_xxx_y_xx_xx = buffer_fpdd[36];

    auto g_xxx_y_xx_xy = buffer_fpdd[37];

    auto g_xxx_y_xx_xz = buffer_fpdd[38];

    auto g_xxx_y_xx_yy = buffer_fpdd[39];

    auto g_xxx_y_xx_yz = buffer_fpdd[40];

    auto g_xxx_y_xx_zz = buffer_fpdd[41];

    auto g_xxx_y_xy_xx = buffer_fpdd[42];

    auto g_xxx_y_xy_xy = buffer_fpdd[43];

    auto g_xxx_y_xy_xz = buffer_fpdd[44];

    auto g_xxx_y_xy_yy = buffer_fpdd[45];

    auto g_xxx_y_xy_yz = buffer_fpdd[46];

    auto g_xxx_y_xy_zz = buffer_fpdd[47];

    auto g_xxx_y_xz_xx = buffer_fpdd[48];

    auto g_xxx_y_xz_xy = buffer_fpdd[49];

    auto g_xxx_y_xz_xz = buffer_fpdd[50];

    auto g_xxx_y_xz_yy = buffer_fpdd[51];

    auto g_xxx_y_xz_yz = buffer_fpdd[52];

    auto g_xxx_y_xz_zz = buffer_fpdd[53];

    auto g_xxx_y_yy_xx = buffer_fpdd[54];

    auto g_xxx_y_yy_xy = buffer_fpdd[55];

    auto g_xxx_y_yy_xz = buffer_fpdd[56];

    auto g_xxx_y_yy_yy = buffer_fpdd[57];

    auto g_xxx_y_yy_yz = buffer_fpdd[58];

    auto g_xxx_y_yy_zz = buffer_fpdd[59];

    auto g_xxx_y_yz_xx = buffer_fpdd[60];

    auto g_xxx_y_yz_xy = buffer_fpdd[61];

    auto g_xxx_y_yz_xz = buffer_fpdd[62];

    auto g_xxx_y_yz_yy = buffer_fpdd[63];

    auto g_xxx_y_yz_yz = buffer_fpdd[64];

    auto g_xxx_y_yz_zz = buffer_fpdd[65];

    auto g_xxx_y_zz_xx = buffer_fpdd[66];

    auto g_xxx_y_zz_xy = buffer_fpdd[67];

    auto g_xxx_y_zz_xz = buffer_fpdd[68];

    auto g_xxx_y_zz_yy = buffer_fpdd[69];

    auto g_xxx_y_zz_yz = buffer_fpdd[70];

    auto g_xxx_y_zz_zz = buffer_fpdd[71];

    auto g_xxx_z_xx_xx = buffer_fpdd[72];

    auto g_xxx_z_xx_xy = buffer_fpdd[73];

    auto g_xxx_z_xx_xz = buffer_fpdd[74];

    auto g_xxx_z_xx_yy = buffer_fpdd[75];

    auto g_xxx_z_xx_yz = buffer_fpdd[76];

    auto g_xxx_z_xx_zz = buffer_fpdd[77];

    auto g_xxx_z_xy_xx = buffer_fpdd[78];

    auto g_xxx_z_xy_xy = buffer_fpdd[79];

    auto g_xxx_z_xy_xz = buffer_fpdd[80];

    auto g_xxx_z_xy_yy = buffer_fpdd[81];

    auto g_xxx_z_xy_yz = buffer_fpdd[82];

    auto g_xxx_z_xy_zz = buffer_fpdd[83];

    auto g_xxx_z_xz_xx = buffer_fpdd[84];

    auto g_xxx_z_xz_xy = buffer_fpdd[85];

    auto g_xxx_z_xz_xz = buffer_fpdd[86];

    auto g_xxx_z_xz_yy = buffer_fpdd[87];

    auto g_xxx_z_xz_yz = buffer_fpdd[88];

    auto g_xxx_z_xz_zz = buffer_fpdd[89];

    auto g_xxx_z_yy_xx = buffer_fpdd[90];

    auto g_xxx_z_yy_xy = buffer_fpdd[91];

    auto g_xxx_z_yy_xz = buffer_fpdd[92];

    auto g_xxx_z_yy_yy = buffer_fpdd[93];

    auto g_xxx_z_yy_yz = buffer_fpdd[94];

    auto g_xxx_z_yy_zz = buffer_fpdd[95];

    auto g_xxx_z_yz_xx = buffer_fpdd[96];

    auto g_xxx_z_yz_xy = buffer_fpdd[97];

    auto g_xxx_z_yz_xz = buffer_fpdd[98];

    auto g_xxx_z_yz_yy = buffer_fpdd[99];

    auto g_xxx_z_yz_yz = buffer_fpdd[100];

    auto g_xxx_z_yz_zz = buffer_fpdd[101];

    auto g_xxx_z_zz_xx = buffer_fpdd[102];

    auto g_xxx_z_zz_xy = buffer_fpdd[103];

    auto g_xxx_z_zz_xz = buffer_fpdd[104];

    auto g_xxx_z_zz_yy = buffer_fpdd[105];

    auto g_xxx_z_zz_yz = buffer_fpdd[106];

    auto g_xxx_z_zz_zz = buffer_fpdd[107];

    auto g_xxy_x_xx_xx = buffer_fpdd[108];

    auto g_xxy_x_xx_xy = buffer_fpdd[109];

    auto g_xxy_x_xx_xz = buffer_fpdd[110];

    auto g_xxy_x_xx_yy = buffer_fpdd[111];

    auto g_xxy_x_xx_yz = buffer_fpdd[112];

    auto g_xxy_x_xx_zz = buffer_fpdd[113];

    auto g_xxy_x_xy_xx = buffer_fpdd[114];

    auto g_xxy_x_xy_xy = buffer_fpdd[115];

    auto g_xxy_x_xy_xz = buffer_fpdd[116];

    auto g_xxy_x_xy_yy = buffer_fpdd[117];

    auto g_xxy_x_xy_yz = buffer_fpdd[118];

    auto g_xxy_x_xy_zz = buffer_fpdd[119];

    auto g_xxy_x_xz_xx = buffer_fpdd[120];

    auto g_xxy_x_xz_xy = buffer_fpdd[121];

    auto g_xxy_x_xz_xz = buffer_fpdd[122];

    auto g_xxy_x_xz_yy = buffer_fpdd[123];

    auto g_xxy_x_xz_yz = buffer_fpdd[124];

    auto g_xxy_x_xz_zz = buffer_fpdd[125];

    auto g_xxy_x_yy_xx = buffer_fpdd[126];

    auto g_xxy_x_yy_xy = buffer_fpdd[127];

    auto g_xxy_x_yy_xz = buffer_fpdd[128];

    auto g_xxy_x_yy_yy = buffer_fpdd[129];

    auto g_xxy_x_yy_yz = buffer_fpdd[130];

    auto g_xxy_x_yy_zz = buffer_fpdd[131];

    auto g_xxy_x_yz_xx = buffer_fpdd[132];

    auto g_xxy_x_yz_xy = buffer_fpdd[133];

    auto g_xxy_x_yz_xz = buffer_fpdd[134];

    auto g_xxy_x_yz_yy = buffer_fpdd[135];

    auto g_xxy_x_yz_yz = buffer_fpdd[136];

    auto g_xxy_x_yz_zz = buffer_fpdd[137];

    auto g_xxy_x_zz_xx = buffer_fpdd[138];

    auto g_xxy_x_zz_xy = buffer_fpdd[139];

    auto g_xxy_x_zz_xz = buffer_fpdd[140];

    auto g_xxy_x_zz_yy = buffer_fpdd[141];

    auto g_xxy_x_zz_yz = buffer_fpdd[142];

    auto g_xxy_x_zz_zz = buffer_fpdd[143];

    auto g_xxy_y_xx_xx = buffer_fpdd[144];

    auto g_xxy_y_xx_xy = buffer_fpdd[145];

    auto g_xxy_y_xx_xz = buffer_fpdd[146];

    auto g_xxy_y_xx_yy = buffer_fpdd[147];

    auto g_xxy_y_xx_yz = buffer_fpdd[148];

    auto g_xxy_y_xx_zz = buffer_fpdd[149];

    auto g_xxy_y_xy_xx = buffer_fpdd[150];

    auto g_xxy_y_xy_xy = buffer_fpdd[151];

    auto g_xxy_y_xy_xz = buffer_fpdd[152];

    auto g_xxy_y_xy_yy = buffer_fpdd[153];

    auto g_xxy_y_xy_yz = buffer_fpdd[154];

    auto g_xxy_y_xy_zz = buffer_fpdd[155];

    auto g_xxy_y_xz_xx = buffer_fpdd[156];

    auto g_xxy_y_xz_xy = buffer_fpdd[157];

    auto g_xxy_y_xz_xz = buffer_fpdd[158];

    auto g_xxy_y_xz_yy = buffer_fpdd[159];

    auto g_xxy_y_xz_yz = buffer_fpdd[160];

    auto g_xxy_y_xz_zz = buffer_fpdd[161];

    auto g_xxy_y_yy_xx = buffer_fpdd[162];

    auto g_xxy_y_yy_xy = buffer_fpdd[163];

    auto g_xxy_y_yy_xz = buffer_fpdd[164];

    auto g_xxy_y_yy_yy = buffer_fpdd[165];

    auto g_xxy_y_yy_yz = buffer_fpdd[166];

    auto g_xxy_y_yy_zz = buffer_fpdd[167];

    auto g_xxy_y_yz_xx = buffer_fpdd[168];

    auto g_xxy_y_yz_xy = buffer_fpdd[169];

    auto g_xxy_y_yz_xz = buffer_fpdd[170];

    auto g_xxy_y_yz_yy = buffer_fpdd[171];

    auto g_xxy_y_yz_yz = buffer_fpdd[172];

    auto g_xxy_y_yz_zz = buffer_fpdd[173];

    auto g_xxy_y_zz_xx = buffer_fpdd[174];

    auto g_xxy_y_zz_xy = buffer_fpdd[175];

    auto g_xxy_y_zz_xz = buffer_fpdd[176];

    auto g_xxy_y_zz_yy = buffer_fpdd[177];

    auto g_xxy_y_zz_yz = buffer_fpdd[178];

    auto g_xxy_y_zz_zz = buffer_fpdd[179];

    auto g_xxy_z_xx_xx = buffer_fpdd[180];

    auto g_xxy_z_xx_xy = buffer_fpdd[181];

    auto g_xxy_z_xx_xz = buffer_fpdd[182];

    auto g_xxy_z_xx_yy = buffer_fpdd[183];

    auto g_xxy_z_xx_yz = buffer_fpdd[184];

    auto g_xxy_z_xx_zz = buffer_fpdd[185];

    auto g_xxy_z_xy_xx = buffer_fpdd[186];

    auto g_xxy_z_xy_xy = buffer_fpdd[187];

    auto g_xxy_z_xy_xz = buffer_fpdd[188];

    auto g_xxy_z_xy_yy = buffer_fpdd[189];

    auto g_xxy_z_xy_yz = buffer_fpdd[190];

    auto g_xxy_z_xy_zz = buffer_fpdd[191];

    auto g_xxy_z_xz_xx = buffer_fpdd[192];

    auto g_xxy_z_xz_xy = buffer_fpdd[193];

    auto g_xxy_z_xz_xz = buffer_fpdd[194];

    auto g_xxy_z_xz_yy = buffer_fpdd[195];

    auto g_xxy_z_xz_yz = buffer_fpdd[196];

    auto g_xxy_z_xz_zz = buffer_fpdd[197];

    auto g_xxy_z_yy_xx = buffer_fpdd[198];

    auto g_xxy_z_yy_xy = buffer_fpdd[199];

    auto g_xxy_z_yy_xz = buffer_fpdd[200];

    auto g_xxy_z_yy_yy = buffer_fpdd[201];

    auto g_xxy_z_yy_yz = buffer_fpdd[202];

    auto g_xxy_z_yy_zz = buffer_fpdd[203];

    auto g_xxy_z_yz_xx = buffer_fpdd[204];

    auto g_xxy_z_yz_xy = buffer_fpdd[205];

    auto g_xxy_z_yz_xz = buffer_fpdd[206];

    auto g_xxy_z_yz_yy = buffer_fpdd[207];

    auto g_xxy_z_yz_yz = buffer_fpdd[208];

    auto g_xxy_z_yz_zz = buffer_fpdd[209];

    auto g_xxy_z_zz_xx = buffer_fpdd[210];

    auto g_xxy_z_zz_xy = buffer_fpdd[211];

    auto g_xxy_z_zz_xz = buffer_fpdd[212];

    auto g_xxy_z_zz_yy = buffer_fpdd[213];

    auto g_xxy_z_zz_yz = buffer_fpdd[214];

    auto g_xxy_z_zz_zz = buffer_fpdd[215];

    auto g_xxz_x_xx_xx = buffer_fpdd[216];

    auto g_xxz_x_xx_xy = buffer_fpdd[217];

    auto g_xxz_x_xx_xz = buffer_fpdd[218];

    auto g_xxz_x_xx_yy = buffer_fpdd[219];

    auto g_xxz_x_xx_yz = buffer_fpdd[220];

    auto g_xxz_x_xx_zz = buffer_fpdd[221];

    auto g_xxz_x_xy_xx = buffer_fpdd[222];

    auto g_xxz_x_xy_xy = buffer_fpdd[223];

    auto g_xxz_x_xy_xz = buffer_fpdd[224];

    auto g_xxz_x_xy_yy = buffer_fpdd[225];

    auto g_xxz_x_xy_yz = buffer_fpdd[226];

    auto g_xxz_x_xy_zz = buffer_fpdd[227];

    auto g_xxz_x_xz_xx = buffer_fpdd[228];

    auto g_xxz_x_xz_xy = buffer_fpdd[229];

    auto g_xxz_x_xz_xz = buffer_fpdd[230];

    auto g_xxz_x_xz_yy = buffer_fpdd[231];

    auto g_xxz_x_xz_yz = buffer_fpdd[232];

    auto g_xxz_x_xz_zz = buffer_fpdd[233];

    auto g_xxz_x_yy_xx = buffer_fpdd[234];

    auto g_xxz_x_yy_xy = buffer_fpdd[235];

    auto g_xxz_x_yy_xz = buffer_fpdd[236];

    auto g_xxz_x_yy_yy = buffer_fpdd[237];

    auto g_xxz_x_yy_yz = buffer_fpdd[238];

    auto g_xxz_x_yy_zz = buffer_fpdd[239];

    auto g_xxz_x_yz_xx = buffer_fpdd[240];

    auto g_xxz_x_yz_xy = buffer_fpdd[241];

    auto g_xxz_x_yz_xz = buffer_fpdd[242];

    auto g_xxz_x_yz_yy = buffer_fpdd[243];

    auto g_xxz_x_yz_yz = buffer_fpdd[244];

    auto g_xxz_x_yz_zz = buffer_fpdd[245];

    auto g_xxz_x_zz_xx = buffer_fpdd[246];

    auto g_xxz_x_zz_xy = buffer_fpdd[247];

    auto g_xxz_x_zz_xz = buffer_fpdd[248];

    auto g_xxz_x_zz_yy = buffer_fpdd[249];

    auto g_xxz_x_zz_yz = buffer_fpdd[250];

    auto g_xxz_x_zz_zz = buffer_fpdd[251];

    auto g_xxz_y_xx_xx = buffer_fpdd[252];

    auto g_xxz_y_xx_xy = buffer_fpdd[253];

    auto g_xxz_y_xx_xz = buffer_fpdd[254];

    auto g_xxz_y_xx_yy = buffer_fpdd[255];

    auto g_xxz_y_xx_yz = buffer_fpdd[256];

    auto g_xxz_y_xx_zz = buffer_fpdd[257];

    auto g_xxz_y_xy_xx = buffer_fpdd[258];

    auto g_xxz_y_xy_xy = buffer_fpdd[259];

    auto g_xxz_y_xy_xz = buffer_fpdd[260];

    auto g_xxz_y_xy_yy = buffer_fpdd[261];

    auto g_xxz_y_xy_yz = buffer_fpdd[262];

    auto g_xxz_y_xy_zz = buffer_fpdd[263];

    auto g_xxz_y_xz_xx = buffer_fpdd[264];

    auto g_xxz_y_xz_xy = buffer_fpdd[265];

    auto g_xxz_y_xz_xz = buffer_fpdd[266];

    auto g_xxz_y_xz_yy = buffer_fpdd[267];

    auto g_xxz_y_xz_yz = buffer_fpdd[268];

    auto g_xxz_y_xz_zz = buffer_fpdd[269];

    auto g_xxz_y_yy_xx = buffer_fpdd[270];

    auto g_xxz_y_yy_xy = buffer_fpdd[271];

    auto g_xxz_y_yy_xz = buffer_fpdd[272];

    auto g_xxz_y_yy_yy = buffer_fpdd[273];

    auto g_xxz_y_yy_yz = buffer_fpdd[274];

    auto g_xxz_y_yy_zz = buffer_fpdd[275];

    auto g_xxz_y_yz_xx = buffer_fpdd[276];

    auto g_xxz_y_yz_xy = buffer_fpdd[277];

    auto g_xxz_y_yz_xz = buffer_fpdd[278];

    auto g_xxz_y_yz_yy = buffer_fpdd[279];

    auto g_xxz_y_yz_yz = buffer_fpdd[280];

    auto g_xxz_y_yz_zz = buffer_fpdd[281];

    auto g_xxz_y_zz_xx = buffer_fpdd[282];

    auto g_xxz_y_zz_xy = buffer_fpdd[283];

    auto g_xxz_y_zz_xz = buffer_fpdd[284];

    auto g_xxz_y_zz_yy = buffer_fpdd[285];

    auto g_xxz_y_zz_yz = buffer_fpdd[286];

    auto g_xxz_y_zz_zz = buffer_fpdd[287];

    auto g_xxz_z_xx_xx = buffer_fpdd[288];

    auto g_xxz_z_xx_xy = buffer_fpdd[289];

    auto g_xxz_z_xx_xz = buffer_fpdd[290];

    auto g_xxz_z_xx_yy = buffer_fpdd[291];

    auto g_xxz_z_xx_yz = buffer_fpdd[292];

    auto g_xxz_z_xx_zz = buffer_fpdd[293];

    auto g_xxz_z_xy_xx = buffer_fpdd[294];

    auto g_xxz_z_xy_xy = buffer_fpdd[295];

    auto g_xxz_z_xy_xz = buffer_fpdd[296];

    auto g_xxz_z_xy_yy = buffer_fpdd[297];

    auto g_xxz_z_xy_yz = buffer_fpdd[298];

    auto g_xxz_z_xy_zz = buffer_fpdd[299];

    auto g_xxz_z_xz_xx = buffer_fpdd[300];

    auto g_xxz_z_xz_xy = buffer_fpdd[301];

    auto g_xxz_z_xz_xz = buffer_fpdd[302];

    auto g_xxz_z_xz_yy = buffer_fpdd[303];

    auto g_xxz_z_xz_yz = buffer_fpdd[304];

    auto g_xxz_z_xz_zz = buffer_fpdd[305];

    auto g_xxz_z_yy_xx = buffer_fpdd[306];

    auto g_xxz_z_yy_xy = buffer_fpdd[307];

    auto g_xxz_z_yy_xz = buffer_fpdd[308];

    auto g_xxz_z_yy_yy = buffer_fpdd[309];

    auto g_xxz_z_yy_yz = buffer_fpdd[310];

    auto g_xxz_z_yy_zz = buffer_fpdd[311];

    auto g_xxz_z_yz_xx = buffer_fpdd[312];

    auto g_xxz_z_yz_xy = buffer_fpdd[313];

    auto g_xxz_z_yz_xz = buffer_fpdd[314];

    auto g_xxz_z_yz_yy = buffer_fpdd[315];

    auto g_xxz_z_yz_yz = buffer_fpdd[316];

    auto g_xxz_z_yz_zz = buffer_fpdd[317];

    auto g_xxz_z_zz_xx = buffer_fpdd[318];

    auto g_xxz_z_zz_xy = buffer_fpdd[319];

    auto g_xxz_z_zz_xz = buffer_fpdd[320];

    auto g_xxz_z_zz_yy = buffer_fpdd[321];

    auto g_xxz_z_zz_yz = buffer_fpdd[322];

    auto g_xxz_z_zz_zz = buffer_fpdd[323];

    auto g_xyy_x_xx_xx = buffer_fpdd[324];

    auto g_xyy_x_xx_xy = buffer_fpdd[325];

    auto g_xyy_x_xx_xz = buffer_fpdd[326];

    auto g_xyy_x_xx_yy = buffer_fpdd[327];

    auto g_xyy_x_xx_yz = buffer_fpdd[328];

    auto g_xyy_x_xx_zz = buffer_fpdd[329];

    auto g_xyy_x_xy_xx = buffer_fpdd[330];

    auto g_xyy_x_xy_xy = buffer_fpdd[331];

    auto g_xyy_x_xy_xz = buffer_fpdd[332];

    auto g_xyy_x_xy_yy = buffer_fpdd[333];

    auto g_xyy_x_xy_yz = buffer_fpdd[334];

    auto g_xyy_x_xy_zz = buffer_fpdd[335];

    auto g_xyy_x_xz_xx = buffer_fpdd[336];

    auto g_xyy_x_xz_xy = buffer_fpdd[337];

    auto g_xyy_x_xz_xz = buffer_fpdd[338];

    auto g_xyy_x_xz_yy = buffer_fpdd[339];

    auto g_xyy_x_xz_yz = buffer_fpdd[340];

    auto g_xyy_x_xz_zz = buffer_fpdd[341];

    auto g_xyy_x_yy_xx = buffer_fpdd[342];

    auto g_xyy_x_yy_xy = buffer_fpdd[343];

    auto g_xyy_x_yy_xz = buffer_fpdd[344];

    auto g_xyy_x_yy_yy = buffer_fpdd[345];

    auto g_xyy_x_yy_yz = buffer_fpdd[346];

    auto g_xyy_x_yy_zz = buffer_fpdd[347];

    auto g_xyy_x_yz_xx = buffer_fpdd[348];

    auto g_xyy_x_yz_xy = buffer_fpdd[349];

    auto g_xyy_x_yz_xz = buffer_fpdd[350];

    auto g_xyy_x_yz_yy = buffer_fpdd[351];

    auto g_xyy_x_yz_yz = buffer_fpdd[352];

    auto g_xyy_x_yz_zz = buffer_fpdd[353];

    auto g_xyy_x_zz_xx = buffer_fpdd[354];

    auto g_xyy_x_zz_xy = buffer_fpdd[355];

    auto g_xyy_x_zz_xz = buffer_fpdd[356];

    auto g_xyy_x_zz_yy = buffer_fpdd[357];

    auto g_xyy_x_zz_yz = buffer_fpdd[358];

    auto g_xyy_x_zz_zz = buffer_fpdd[359];

    auto g_xyy_y_xx_xx = buffer_fpdd[360];

    auto g_xyy_y_xx_xy = buffer_fpdd[361];

    auto g_xyy_y_xx_xz = buffer_fpdd[362];

    auto g_xyy_y_xx_yy = buffer_fpdd[363];

    auto g_xyy_y_xx_yz = buffer_fpdd[364];

    auto g_xyy_y_xx_zz = buffer_fpdd[365];

    auto g_xyy_y_xy_xx = buffer_fpdd[366];

    auto g_xyy_y_xy_xy = buffer_fpdd[367];

    auto g_xyy_y_xy_xz = buffer_fpdd[368];

    auto g_xyy_y_xy_yy = buffer_fpdd[369];

    auto g_xyy_y_xy_yz = buffer_fpdd[370];

    auto g_xyy_y_xy_zz = buffer_fpdd[371];

    auto g_xyy_y_xz_xx = buffer_fpdd[372];

    auto g_xyy_y_xz_xy = buffer_fpdd[373];

    auto g_xyy_y_xz_xz = buffer_fpdd[374];

    auto g_xyy_y_xz_yy = buffer_fpdd[375];

    auto g_xyy_y_xz_yz = buffer_fpdd[376];

    auto g_xyy_y_xz_zz = buffer_fpdd[377];

    auto g_xyy_y_yy_xx = buffer_fpdd[378];

    auto g_xyy_y_yy_xy = buffer_fpdd[379];

    auto g_xyy_y_yy_xz = buffer_fpdd[380];

    auto g_xyy_y_yy_yy = buffer_fpdd[381];

    auto g_xyy_y_yy_yz = buffer_fpdd[382];

    auto g_xyy_y_yy_zz = buffer_fpdd[383];

    auto g_xyy_y_yz_xx = buffer_fpdd[384];

    auto g_xyy_y_yz_xy = buffer_fpdd[385];

    auto g_xyy_y_yz_xz = buffer_fpdd[386];

    auto g_xyy_y_yz_yy = buffer_fpdd[387];

    auto g_xyy_y_yz_yz = buffer_fpdd[388];

    auto g_xyy_y_yz_zz = buffer_fpdd[389];

    auto g_xyy_y_zz_xx = buffer_fpdd[390];

    auto g_xyy_y_zz_xy = buffer_fpdd[391];

    auto g_xyy_y_zz_xz = buffer_fpdd[392];

    auto g_xyy_y_zz_yy = buffer_fpdd[393];

    auto g_xyy_y_zz_yz = buffer_fpdd[394];

    auto g_xyy_y_zz_zz = buffer_fpdd[395];

    auto g_xyy_z_xx_xx = buffer_fpdd[396];

    auto g_xyy_z_xx_xy = buffer_fpdd[397];

    auto g_xyy_z_xx_xz = buffer_fpdd[398];

    auto g_xyy_z_xx_yy = buffer_fpdd[399];

    auto g_xyy_z_xx_yz = buffer_fpdd[400];

    auto g_xyy_z_xx_zz = buffer_fpdd[401];

    auto g_xyy_z_xy_xx = buffer_fpdd[402];

    auto g_xyy_z_xy_xy = buffer_fpdd[403];

    auto g_xyy_z_xy_xz = buffer_fpdd[404];

    auto g_xyy_z_xy_yy = buffer_fpdd[405];

    auto g_xyy_z_xy_yz = buffer_fpdd[406];

    auto g_xyy_z_xy_zz = buffer_fpdd[407];

    auto g_xyy_z_xz_xx = buffer_fpdd[408];

    auto g_xyy_z_xz_xy = buffer_fpdd[409];

    auto g_xyy_z_xz_xz = buffer_fpdd[410];

    auto g_xyy_z_xz_yy = buffer_fpdd[411];

    auto g_xyy_z_xz_yz = buffer_fpdd[412];

    auto g_xyy_z_xz_zz = buffer_fpdd[413];

    auto g_xyy_z_yy_xx = buffer_fpdd[414];

    auto g_xyy_z_yy_xy = buffer_fpdd[415];

    auto g_xyy_z_yy_xz = buffer_fpdd[416];

    auto g_xyy_z_yy_yy = buffer_fpdd[417];

    auto g_xyy_z_yy_yz = buffer_fpdd[418];

    auto g_xyy_z_yy_zz = buffer_fpdd[419];

    auto g_xyy_z_yz_xx = buffer_fpdd[420];

    auto g_xyy_z_yz_xy = buffer_fpdd[421];

    auto g_xyy_z_yz_xz = buffer_fpdd[422];

    auto g_xyy_z_yz_yy = buffer_fpdd[423];

    auto g_xyy_z_yz_yz = buffer_fpdd[424];

    auto g_xyy_z_yz_zz = buffer_fpdd[425];

    auto g_xyy_z_zz_xx = buffer_fpdd[426];

    auto g_xyy_z_zz_xy = buffer_fpdd[427];

    auto g_xyy_z_zz_xz = buffer_fpdd[428];

    auto g_xyy_z_zz_yy = buffer_fpdd[429];

    auto g_xyy_z_zz_yz = buffer_fpdd[430];

    auto g_xyy_z_zz_zz = buffer_fpdd[431];

    auto g_xyz_x_xx_xx = buffer_fpdd[432];

    auto g_xyz_x_xx_xy = buffer_fpdd[433];

    auto g_xyz_x_xx_xz = buffer_fpdd[434];

    auto g_xyz_x_xx_yy = buffer_fpdd[435];

    auto g_xyz_x_xx_yz = buffer_fpdd[436];

    auto g_xyz_x_xx_zz = buffer_fpdd[437];

    auto g_xyz_x_xy_xx = buffer_fpdd[438];

    auto g_xyz_x_xy_xy = buffer_fpdd[439];

    auto g_xyz_x_xy_xz = buffer_fpdd[440];

    auto g_xyz_x_xy_yy = buffer_fpdd[441];

    auto g_xyz_x_xy_yz = buffer_fpdd[442];

    auto g_xyz_x_xy_zz = buffer_fpdd[443];

    auto g_xyz_x_xz_xx = buffer_fpdd[444];

    auto g_xyz_x_xz_xy = buffer_fpdd[445];

    auto g_xyz_x_xz_xz = buffer_fpdd[446];

    auto g_xyz_x_xz_yy = buffer_fpdd[447];

    auto g_xyz_x_xz_yz = buffer_fpdd[448];

    auto g_xyz_x_xz_zz = buffer_fpdd[449];

    auto g_xyz_x_yy_xx = buffer_fpdd[450];

    auto g_xyz_x_yy_xy = buffer_fpdd[451];

    auto g_xyz_x_yy_xz = buffer_fpdd[452];

    auto g_xyz_x_yy_yy = buffer_fpdd[453];

    auto g_xyz_x_yy_yz = buffer_fpdd[454];

    auto g_xyz_x_yy_zz = buffer_fpdd[455];

    auto g_xyz_x_yz_xx = buffer_fpdd[456];

    auto g_xyz_x_yz_xy = buffer_fpdd[457];

    auto g_xyz_x_yz_xz = buffer_fpdd[458];

    auto g_xyz_x_yz_yy = buffer_fpdd[459];

    auto g_xyz_x_yz_yz = buffer_fpdd[460];

    auto g_xyz_x_yz_zz = buffer_fpdd[461];

    auto g_xyz_x_zz_xx = buffer_fpdd[462];

    auto g_xyz_x_zz_xy = buffer_fpdd[463];

    auto g_xyz_x_zz_xz = buffer_fpdd[464];

    auto g_xyz_x_zz_yy = buffer_fpdd[465];

    auto g_xyz_x_zz_yz = buffer_fpdd[466];

    auto g_xyz_x_zz_zz = buffer_fpdd[467];

    auto g_xyz_y_xx_xx = buffer_fpdd[468];

    auto g_xyz_y_xx_xy = buffer_fpdd[469];

    auto g_xyz_y_xx_xz = buffer_fpdd[470];

    auto g_xyz_y_xx_yy = buffer_fpdd[471];

    auto g_xyz_y_xx_yz = buffer_fpdd[472];

    auto g_xyz_y_xx_zz = buffer_fpdd[473];

    auto g_xyz_y_xy_xx = buffer_fpdd[474];

    auto g_xyz_y_xy_xy = buffer_fpdd[475];

    auto g_xyz_y_xy_xz = buffer_fpdd[476];

    auto g_xyz_y_xy_yy = buffer_fpdd[477];

    auto g_xyz_y_xy_yz = buffer_fpdd[478];

    auto g_xyz_y_xy_zz = buffer_fpdd[479];

    auto g_xyz_y_xz_xx = buffer_fpdd[480];

    auto g_xyz_y_xz_xy = buffer_fpdd[481];

    auto g_xyz_y_xz_xz = buffer_fpdd[482];

    auto g_xyz_y_xz_yy = buffer_fpdd[483];

    auto g_xyz_y_xz_yz = buffer_fpdd[484];

    auto g_xyz_y_xz_zz = buffer_fpdd[485];

    auto g_xyz_y_yy_xx = buffer_fpdd[486];

    auto g_xyz_y_yy_xy = buffer_fpdd[487];

    auto g_xyz_y_yy_xz = buffer_fpdd[488];

    auto g_xyz_y_yy_yy = buffer_fpdd[489];

    auto g_xyz_y_yy_yz = buffer_fpdd[490];

    auto g_xyz_y_yy_zz = buffer_fpdd[491];

    auto g_xyz_y_yz_xx = buffer_fpdd[492];

    auto g_xyz_y_yz_xy = buffer_fpdd[493];

    auto g_xyz_y_yz_xz = buffer_fpdd[494];

    auto g_xyz_y_yz_yy = buffer_fpdd[495];

    auto g_xyz_y_yz_yz = buffer_fpdd[496];

    auto g_xyz_y_yz_zz = buffer_fpdd[497];

    auto g_xyz_y_zz_xx = buffer_fpdd[498];

    auto g_xyz_y_zz_xy = buffer_fpdd[499];

    auto g_xyz_y_zz_xz = buffer_fpdd[500];

    auto g_xyz_y_zz_yy = buffer_fpdd[501];

    auto g_xyz_y_zz_yz = buffer_fpdd[502];

    auto g_xyz_y_zz_zz = buffer_fpdd[503];

    auto g_xyz_z_xx_xx = buffer_fpdd[504];

    auto g_xyz_z_xx_xy = buffer_fpdd[505];

    auto g_xyz_z_xx_xz = buffer_fpdd[506];

    auto g_xyz_z_xx_yy = buffer_fpdd[507];

    auto g_xyz_z_xx_yz = buffer_fpdd[508];

    auto g_xyz_z_xx_zz = buffer_fpdd[509];

    auto g_xyz_z_xy_xx = buffer_fpdd[510];

    auto g_xyz_z_xy_xy = buffer_fpdd[511];

    auto g_xyz_z_xy_xz = buffer_fpdd[512];

    auto g_xyz_z_xy_yy = buffer_fpdd[513];

    auto g_xyz_z_xy_yz = buffer_fpdd[514];

    auto g_xyz_z_xy_zz = buffer_fpdd[515];

    auto g_xyz_z_xz_xx = buffer_fpdd[516];

    auto g_xyz_z_xz_xy = buffer_fpdd[517];

    auto g_xyz_z_xz_xz = buffer_fpdd[518];

    auto g_xyz_z_xz_yy = buffer_fpdd[519];

    auto g_xyz_z_xz_yz = buffer_fpdd[520];

    auto g_xyz_z_xz_zz = buffer_fpdd[521];

    auto g_xyz_z_yy_xx = buffer_fpdd[522];

    auto g_xyz_z_yy_xy = buffer_fpdd[523];

    auto g_xyz_z_yy_xz = buffer_fpdd[524];

    auto g_xyz_z_yy_yy = buffer_fpdd[525];

    auto g_xyz_z_yy_yz = buffer_fpdd[526];

    auto g_xyz_z_yy_zz = buffer_fpdd[527];

    auto g_xyz_z_yz_xx = buffer_fpdd[528];

    auto g_xyz_z_yz_xy = buffer_fpdd[529];

    auto g_xyz_z_yz_xz = buffer_fpdd[530];

    auto g_xyz_z_yz_yy = buffer_fpdd[531];

    auto g_xyz_z_yz_yz = buffer_fpdd[532];

    auto g_xyz_z_yz_zz = buffer_fpdd[533];

    auto g_xyz_z_zz_xx = buffer_fpdd[534];

    auto g_xyz_z_zz_xy = buffer_fpdd[535];

    auto g_xyz_z_zz_xz = buffer_fpdd[536];

    auto g_xyz_z_zz_yy = buffer_fpdd[537];

    auto g_xyz_z_zz_yz = buffer_fpdd[538];

    auto g_xyz_z_zz_zz = buffer_fpdd[539];

    auto g_xzz_x_xx_xx = buffer_fpdd[540];

    auto g_xzz_x_xx_xy = buffer_fpdd[541];

    auto g_xzz_x_xx_xz = buffer_fpdd[542];

    auto g_xzz_x_xx_yy = buffer_fpdd[543];

    auto g_xzz_x_xx_yz = buffer_fpdd[544];

    auto g_xzz_x_xx_zz = buffer_fpdd[545];

    auto g_xzz_x_xy_xx = buffer_fpdd[546];

    auto g_xzz_x_xy_xy = buffer_fpdd[547];

    auto g_xzz_x_xy_xz = buffer_fpdd[548];

    auto g_xzz_x_xy_yy = buffer_fpdd[549];

    auto g_xzz_x_xy_yz = buffer_fpdd[550];

    auto g_xzz_x_xy_zz = buffer_fpdd[551];

    auto g_xzz_x_xz_xx = buffer_fpdd[552];

    auto g_xzz_x_xz_xy = buffer_fpdd[553];

    auto g_xzz_x_xz_xz = buffer_fpdd[554];

    auto g_xzz_x_xz_yy = buffer_fpdd[555];

    auto g_xzz_x_xz_yz = buffer_fpdd[556];

    auto g_xzz_x_xz_zz = buffer_fpdd[557];

    auto g_xzz_x_yy_xx = buffer_fpdd[558];

    auto g_xzz_x_yy_xy = buffer_fpdd[559];

    auto g_xzz_x_yy_xz = buffer_fpdd[560];

    auto g_xzz_x_yy_yy = buffer_fpdd[561];

    auto g_xzz_x_yy_yz = buffer_fpdd[562];

    auto g_xzz_x_yy_zz = buffer_fpdd[563];

    auto g_xzz_x_yz_xx = buffer_fpdd[564];

    auto g_xzz_x_yz_xy = buffer_fpdd[565];

    auto g_xzz_x_yz_xz = buffer_fpdd[566];

    auto g_xzz_x_yz_yy = buffer_fpdd[567];

    auto g_xzz_x_yz_yz = buffer_fpdd[568];

    auto g_xzz_x_yz_zz = buffer_fpdd[569];

    auto g_xzz_x_zz_xx = buffer_fpdd[570];

    auto g_xzz_x_zz_xy = buffer_fpdd[571];

    auto g_xzz_x_zz_xz = buffer_fpdd[572];

    auto g_xzz_x_zz_yy = buffer_fpdd[573];

    auto g_xzz_x_zz_yz = buffer_fpdd[574];

    auto g_xzz_x_zz_zz = buffer_fpdd[575];

    auto g_xzz_y_xx_xx = buffer_fpdd[576];

    auto g_xzz_y_xx_xy = buffer_fpdd[577];

    auto g_xzz_y_xx_xz = buffer_fpdd[578];

    auto g_xzz_y_xx_yy = buffer_fpdd[579];

    auto g_xzz_y_xx_yz = buffer_fpdd[580];

    auto g_xzz_y_xx_zz = buffer_fpdd[581];

    auto g_xzz_y_xy_xx = buffer_fpdd[582];

    auto g_xzz_y_xy_xy = buffer_fpdd[583];

    auto g_xzz_y_xy_xz = buffer_fpdd[584];

    auto g_xzz_y_xy_yy = buffer_fpdd[585];

    auto g_xzz_y_xy_yz = buffer_fpdd[586];

    auto g_xzz_y_xy_zz = buffer_fpdd[587];

    auto g_xzz_y_xz_xx = buffer_fpdd[588];

    auto g_xzz_y_xz_xy = buffer_fpdd[589];

    auto g_xzz_y_xz_xz = buffer_fpdd[590];

    auto g_xzz_y_xz_yy = buffer_fpdd[591];

    auto g_xzz_y_xz_yz = buffer_fpdd[592];

    auto g_xzz_y_xz_zz = buffer_fpdd[593];

    auto g_xzz_y_yy_xx = buffer_fpdd[594];

    auto g_xzz_y_yy_xy = buffer_fpdd[595];

    auto g_xzz_y_yy_xz = buffer_fpdd[596];

    auto g_xzz_y_yy_yy = buffer_fpdd[597];

    auto g_xzz_y_yy_yz = buffer_fpdd[598];

    auto g_xzz_y_yy_zz = buffer_fpdd[599];

    auto g_xzz_y_yz_xx = buffer_fpdd[600];

    auto g_xzz_y_yz_xy = buffer_fpdd[601];

    auto g_xzz_y_yz_xz = buffer_fpdd[602];

    auto g_xzz_y_yz_yy = buffer_fpdd[603];

    auto g_xzz_y_yz_yz = buffer_fpdd[604];

    auto g_xzz_y_yz_zz = buffer_fpdd[605];

    auto g_xzz_y_zz_xx = buffer_fpdd[606];

    auto g_xzz_y_zz_xy = buffer_fpdd[607];

    auto g_xzz_y_zz_xz = buffer_fpdd[608];

    auto g_xzz_y_zz_yy = buffer_fpdd[609];

    auto g_xzz_y_zz_yz = buffer_fpdd[610];

    auto g_xzz_y_zz_zz = buffer_fpdd[611];

    auto g_xzz_z_xx_xx = buffer_fpdd[612];

    auto g_xzz_z_xx_xy = buffer_fpdd[613];

    auto g_xzz_z_xx_xz = buffer_fpdd[614];

    auto g_xzz_z_xx_yy = buffer_fpdd[615];

    auto g_xzz_z_xx_yz = buffer_fpdd[616];

    auto g_xzz_z_xx_zz = buffer_fpdd[617];

    auto g_xzz_z_xy_xx = buffer_fpdd[618];

    auto g_xzz_z_xy_xy = buffer_fpdd[619];

    auto g_xzz_z_xy_xz = buffer_fpdd[620];

    auto g_xzz_z_xy_yy = buffer_fpdd[621];

    auto g_xzz_z_xy_yz = buffer_fpdd[622];

    auto g_xzz_z_xy_zz = buffer_fpdd[623];

    auto g_xzz_z_xz_xx = buffer_fpdd[624];

    auto g_xzz_z_xz_xy = buffer_fpdd[625];

    auto g_xzz_z_xz_xz = buffer_fpdd[626];

    auto g_xzz_z_xz_yy = buffer_fpdd[627];

    auto g_xzz_z_xz_yz = buffer_fpdd[628];

    auto g_xzz_z_xz_zz = buffer_fpdd[629];

    auto g_xzz_z_yy_xx = buffer_fpdd[630];

    auto g_xzz_z_yy_xy = buffer_fpdd[631];

    auto g_xzz_z_yy_xz = buffer_fpdd[632];

    auto g_xzz_z_yy_yy = buffer_fpdd[633];

    auto g_xzz_z_yy_yz = buffer_fpdd[634];

    auto g_xzz_z_yy_zz = buffer_fpdd[635];

    auto g_xzz_z_yz_xx = buffer_fpdd[636];

    auto g_xzz_z_yz_xy = buffer_fpdd[637];

    auto g_xzz_z_yz_xz = buffer_fpdd[638];

    auto g_xzz_z_yz_yy = buffer_fpdd[639];

    auto g_xzz_z_yz_yz = buffer_fpdd[640];

    auto g_xzz_z_yz_zz = buffer_fpdd[641];

    auto g_xzz_z_zz_xx = buffer_fpdd[642];

    auto g_xzz_z_zz_xy = buffer_fpdd[643];

    auto g_xzz_z_zz_xz = buffer_fpdd[644];

    auto g_xzz_z_zz_yy = buffer_fpdd[645];

    auto g_xzz_z_zz_yz = buffer_fpdd[646];

    auto g_xzz_z_zz_zz = buffer_fpdd[647];

    auto g_yyy_x_xx_xx = buffer_fpdd[648];

    auto g_yyy_x_xx_xy = buffer_fpdd[649];

    auto g_yyy_x_xx_xz = buffer_fpdd[650];

    auto g_yyy_x_xx_yy = buffer_fpdd[651];

    auto g_yyy_x_xx_yz = buffer_fpdd[652];

    auto g_yyy_x_xx_zz = buffer_fpdd[653];

    auto g_yyy_x_xy_xx = buffer_fpdd[654];

    auto g_yyy_x_xy_xy = buffer_fpdd[655];

    auto g_yyy_x_xy_xz = buffer_fpdd[656];

    auto g_yyy_x_xy_yy = buffer_fpdd[657];

    auto g_yyy_x_xy_yz = buffer_fpdd[658];

    auto g_yyy_x_xy_zz = buffer_fpdd[659];

    auto g_yyy_x_xz_xx = buffer_fpdd[660];

    auto g_yyy_x_xz_xy = buffer_fpdd[661];

    auto g_yyy_x_xz_xz = buffer_fpdd[662];

    auto g_yyy_x_xz_yy = buffer_fpdd[663];

    auto g_yyy_x_xz_yz = buffer_fpdd[664];

    auto g_yyy_x_xz_zz = buffer_fpdd[665];

    auto g_yyy_x_yy_xx = buffer_fpdd[666];

    auto g_yyy_x_yy_xy = buffer_fpdd[667];

    auto g_yyy_x_yy_xz = buffer_fpdd[668];

    auto g_yyy_x_yy_yy = buffer_fpdd[669];

    auto g_yyy_x_yy_yz = buffer_fpdd[670];

    auto g_yyy_x_yy_zz = buffer_fpdd[671];

    auto g_yyy_x_yz_xx = buffer_fpdd[672];

    auto g_yyy_x_yz_xy = buffer_fpdd[673];

    auto g_yyy_x_yz_xz = buffer_fpdd[674];

    auto g_yyy_x_yz_yy = buffer_fpdd[675];

    auto g_yyy_x_yz_yz = buffer_fpdd[676];

    auto g_yyy_x_yz_zz = buffer_fpdd[677];

    auto g_yyy_x_zz_xx = buffer_fpdd[678];

    auto g_yyy_x_zz_xy = buffer_fpdd[679];

    auto g_yyy_x_zz_xz = buffer_fpdd[680];

    auto g_yyy_x_zz_yy = buffer_fpdd[681];

    auto g_yyy_x_zz_yz = buffer_fpdd[682];

    auto g_yyy_x_zz_zz = buffer_fpdd[683];

    auto g_yyy_y_xx_xx = buffer_fpdd[684];

    auto g_yyy_y_xx_xy = buffer_fpdd[685];

    auto g_yyy_y_xx_xz = buffer_fpdd[686];

    auto g_yyy_y_xx_yy = buffer_fpdd[687];

    auto g_yyy_y_xx_yz = buffer_fpdd[688];

    auto g_yyy_y_xx_zz = buffer_fpdd[689];

    auto g_yyy_y_xy_xx = buffer_fpdd[690];

    auto g_yyy_y_xy_xy = buffer_fpdd[691];

    auto g_yyy_y_xy_xz = buffer_fpdd[692];

    auto g_yyy_y_xy_yy = buffer_fpdd[693];

    auto g_yyy_y_xy_yz = buffer_fpdd[694];

    auto g_yyy_y_xy_zz = buffer_fpdd[695];

    auto g_yyy_y_xz_xx = buffer_fpdd[696];

    auto g_yyy_y_xz_xy = buffer_fpdd[697];

    auto g_yyy_y_xz_xz = buffer_fpdd[698];

    auto g_yyy_y_xz_yy = buffer_fpdd[699];

    auto g_yyy_y_xz_yz = buffer_fpdd[700];

    auto g_yyy_y_xz_zz = buffer_fpdd[701];

    auto g_yyy_y_yy_xx = buffer_fpdd[702];

    auto g_yyy_y_yy_xy = buffer_fpdd[703];

    auto g_yyy_y_yy_xz = buffer_fpdd[704];

    auto g_yyy_y_yy_yy = buffer_fpdd[705];

    auto g_yyy_y_yy_yz = buffer_fpdd[706];

    auto g_yyy_y_yy_zz = buffer_fpdd[707];

    auto g_yyy_y_yz_xx = buffer_fpdd[708];

    auto g_yyy_y_yz_xy = buffer_fpdd[709];

    auto g_yyy_y_yz_xz = buffer_fpdd[710];

    auto g_yyy_y_yz_yy = buffer_fpdd[711];

    auto g_yyy_y_yz_yz = buffer_fpdd[712];

    auto g_yyy_y_yz_zz = buffer_fpdd[713];

    auto g_yyy_y_zz_xx = buffer_fpdd[714];

    auto g_yyy_y_zz_xy = buffer_fpdd[715];

    auto g_yyy_y_zz_xz = buffer_fpdd[716];

    auto g_yyy_y_zz_yy = buffer_fpdd[717];

    auto g_yyy_y_zz_yz = buffer_fpdd[718];

    auto g_yyy_y_zz_zz = buffer_fpdd[719];

    auto g_yyy_z_xx_xx = buffer_fpdd[720];

    auto g_yyy_z_xx_xy = buffer_fpdd[721];

    auto g_yyy_z_xx_xz = buffer_fpdd[722];

    auto g_yyy_z_xx_yy = buffer_fpdd[723];

    auto g_yyy_z_xx_yz = buffer_fpdd[724];

    auto g_yyy_z_xx_zz = buffer_fpdd[725];

    auto g_yyy_z_xy_xx = buffer_fpdd[726];

    auto g_yyy_z_xy_xy = buffer_fpdd[727];

    auto g_yyy_z_xy_xz = buffer_fpdd[728];

    auto g_yyy_z_xy_yy = buffer_fpdd[729];

    auto g_yyy_z_xy_yz = buffer_fpdd[730];

    auto g_yyy_z_xy_zz = buffer_fpdd[731];

    auto g_yyy_z_xz_xx = buffer_fpdd[732];

    auto g_yyy_z_xz_xy = buffer_fpdd[733];

    auto g_yyy_z_xz_xz = buffer_fpdd[734];

    auto g_yyy_z_xz_yy = buffer_fpdd[735];

    auto g_yyy_z_xz_yz = buffer_fpdd[736];

    auto g_yyy_z_xz_zz = buffer_fpdd[737];

    auto g_yyy_z_yy_xx = buffer_fpdd[738];

    auto g_yyy_z_yy_xy = buffer_fpdd[739];

    auto g_yyy_z_yy_xz = buffer_fpdd[740];

    auto g_yyy_z_yy_yy = buffer_fpdd[741];

    auto g_yyy_z_yy_yz = buffer_fpdd[742];

    auto g_yyy_z_yy_zz = buffer_fpdd[743];

    auto g_yyy_z_yz_xx = buffer_fpdd[744];

    auto g_yyy_z_yz_xy = buffer_fpdd[745];

    auto g_yyy_z_yz_xz = buffer_fpdd[746];

    auto g_yyy_z_yz_yy = buffer_fpdd[747];

    auto g_yyy_z_yz_yz = buffer_fpdd[748];

    auto g_yyy_z_yz_zz = buffer_fpdd[749];

    auto g_yyy_z_zz_xx = buffer_fpdd[750];

    auto g_yyy_z_zz_xy = buffer_fpdd[751];

    auto g_yyy_z_zz_xz = buffer_fpdd[752];

    auto g_yyy_z_zz_yy = buffer_fpdd[753];

    auto g_yyy_z_zz_yz = buffer_fpdd[754];

    auto g_yyy_z_zz_zz = buffer_fpdd[755];

    auto g_yyz_x_xx_xx = buffer_fpdd[756];

    auto g_yyz_x_xx_xy = buffer_fpdd[757];

    auto g_yyz_x_xx_xz = buffer_fpdd[758];

    auto g_yyz_x_xx_yy = buffer_fpdd[759];

    auto g_yyz_x_xx_yz = buffer_fpdd[760];

    auto g_yyz_x_xx_zz = buffer_fpdd[761];

    auto g_yyz_x_xy_xx = buffer_fpdd[762];

    auto g_yyz_x_xy_xy = buffer_fpdd[763];

    auto g_yyz_x_xy_xz = buffer_fpdd[764];

    auto g_yyz_x_xy_yy = buffer_fpdd[765];

    auto g_yyz_x_xy_yz = buffer_fpdd[766];

    auto g_yyz_x_xy_zz = buffer_fpdd[767];

    auto g_yyz_x_xz_xx = buffer_fpdd[768];

    auto g_yyz_x_xz_xy = buffer_fpdd[769];

    auto g_yyz_x_xz_xz = buffer_fpdd[770];

    auto g_yyz_x_xz_yy = buffer_fpdd[771];

    auto g_yyz_x_xz_yz = buffer_fpdd[772];

    auto g_yyz_x_xz_zz = buffer_fpdd[773];

    auto g_yyz_x_yy_xx = buffer_fpdd[774];

    auto g_yyz_x_yy_xy = buffer_fpdd[775];

    auto g_yyz_x_yy_xz = buffer_fpdd[776];

    auto g_yyz_x_yy_yy = buffer_fpdd[777];

    auto g_yyz_x_yy_yz = buffer_fpdd[778];

    auto g_yyz_x_yy_zz = buffer_fpdd[779];

    auto g_yyz_x_yz_xx = buffer_fpdd[780];

    auto g_yyz_x_yz_xy = buffer_fpdd[781];

    auto g_yyz_x_yz_xz = buffer_fpdd[782];

    auto g_yyz_x_yz_yy = buffer_fpdd[783];

    auto g_yyz_x_yz_yz = buffer_fpdd[784];

    auto g_yyz_x_yz_zz = buffer_fpdd[785];

    auto g_yyz_x_zz_xx = buffer_fpdd[786];

    auto g_yyz_x_zz_xy = buffer_fpdd[787];

    auto g_yyz_x_zz_xz = buffer_fpdd[788];

    auto g_yyz_x_zz_yy = buffer_fpdd[789];

    auto g_yyz_x_zz_yz = buffer_fpdd[790];

    auto g_yyz_x_zz_zz = buffer_fpdd[791];

    auto g_yyz_y_xx_xx = buffer_fpdd[792];

    auto g_yyz_y_xx_xy = buffer_fpdd[793];

    auto g_yyz_y_xx_xz = buffer_fpdd[794];

    auto g_yyz_y_xx_yy = buffer_fpdd[795];

    auto g_yyz_y_xx_yz = buffer_fpdd[796];

    auto g_yyz_y_xx_zz = buffer_fpdd[797];

    auto g_yyz_y_xy_xx = buffer_fpdd[798];

    auto g_yyz_y_xy_xy = buffer_fpdd[799];

    auto g_yyz_y_xy_xz = buffer_fpdd[800];

    auto g_yyz_y_xy_yy = buffer_fpdd[801];

    auto g_yyz_y_xy_yz = buffer_fpdd[802];

    auto g_yyz_y_xy_zz = buffer_fpdd[803];

    auto g_yyz_y_xz_xx = buffer_fpdd[804];

    auto g_yyz_y_xz_xy = buffer_fpdd[805];

    auto g_yyz_y_xz_xz = buffer_fpdd[806];

    auto g_yyz_y_xz_yy = buffer_fpdd[807];

    auto g_yyz_y_xz_yz = buffer_fpdd[808];

    auto g_yyz_y_xz_zz = buffer_fpdd[809];

    auto g_yyz_y_yy_xx = buffer_fpdd[810];

    auto g_yyz_y_yy_xy = buffer_fpdd[811];

    auto g_yyz_y_yy_xz = buffer_fpdd[812];

    auto g_yyz_y_yy_yy = buffer_fpdd[813];

    auto g_yyz_y_yy_yz = buffer_fpdd[814];

    auto g_yyz_y_yy_zz = buffer_fpdd[815];

    auto g_yyz_y_yz_xx = buffer_fpdd[816];

    auto g_yyz_y_yz_xy = buffer_fpdd[817];

    auto g_yyz_y_yz_xz = buffer_fpdd[818];

    auto g_yyz_y_yz_yy = buffer_fpdd[819];

    auto g_yyz_y_yz_yz = buffer_fpdd[820];

    auto g_yyz_y_yz_zz = buffer_fpdd[821];

    auto g_yyz_y_zz_xx = buffer_fpdd[822];

    auto g_yyz_y_zz_xy = buffer_fpdd[823];

    auto g_yyz_y_zz_xz = buffer_fpdd[824];

    auto g_yyz_y_zz_yy = buffer_fpdd[825];

    auto g_yyz_y_zz_yz = buffer_fpdd[826];

    auto g_yyz_y_zz_zz = buffer_fpdd[827];

    auto g_yyz_z_xx_xx = buffer_fpdd[828];

    auto g_yyz_z_xx_xy = buffer_fpdd[829];

    auto g_yyz_z_xx_xz = buffer_fpdd[830];

    auto g_yyz_z_xx_yy = buffer_fpdd[831];

    auto g_yyz_z_xx_yz = buffer_fpdd[832];

    auto g_yyz_z_xx_zz = buffer_fpdd[833];

    auto g_yyz_z_xy_xx = buffer_fpdd[834];

    auto g_yyz_z_xy_xy = buffer_fpdd[835];

    auto g_yyz_z_xy_xz = buffer_fpdd[836];

    auto g_yyz_z_xy_yy = buffer_fpdd[837];

    auto g_yyz_z_xy_yz = buffer_fpdd[838];

    auto g_yyz_z_xy_zz = buffer_fpdd[839];

    auto g_yyz_z_xz_xx = buffer_fpdd[840];

    auto g_yyz_z_xz_xy = buffer_fpdd[841];

    auto g_yyz_z_xz_xz = buffer_fpdd[842];

    auto g_yyz_z_xz_yy = buffer_fpdd[843];

    auto g_yyz_z_xz_yz = buffer_fpdd[844];

    auto g_yyz_z_xz_zz = buffer_fpdd[845];

    auto g_yyz_z_yy_xx = buffer_fpdd[846];

    auto g_yyz_z_yy_xy = buffer_fpdd[847];

    auto g_yyz_z_yy_xz = buffer_fpdd[848];

    auto g_yyz_z_yy_yy = buffer_fpdd[849];

    auto g_yyz_z_yy_yz = buffer_fpdd[850];

    auto g_yyz_z_yy_zz = buffer_fpdd[851];

    auto g_yyz_z_yz_xx = buffer_fpdd[852];

    auto g_yyz_z_yz_xy = buffer_fpdd[853];

    auto g_yyz_z_yz_xz = buffer_fpdd[854];

    auto g_yyz_z_yz_yy = buffer_fpdd[855];

    auto g_yyz_z_yz_yz = buffer_fpdd[856];

    auto g_yyz_z_yz_zz = buffer_fpdd[857];

    auto g_yyz_z_zz_xx = buffer_fpdd[858];

    auto g_yyz_z_zz_xy = buffer_fpdd[859];

    auto g_yyz_z_zz_xz = buffer_fpdd[860];

    auto g_yyz_z_zz_yy = buffer_fpdd[861];

    auto g_yyz_z_zz_yz = buffer_fpdd[862];

    auto g_yyz_z_zz_zz = buffer_fpdd[863];

    auto g_yzz_x_xx_xx = buffer_fpdd[864];

    auto g_yzz_x_xx_xy = buffer_fpdd[865];

    auto g_yzz_x_xx_xz = buffer_fpdd[866];

    auto g_yzz_x_xx_yy = buffer_fpdd[867];

    auto g_yzz_x_xx_yz = buffer_fpdd[868];

    auto g_yzz_x_xx_zz = buffer_fpdd[869];

    auto g_yzz_x_xy_xx = buffer_fpdd[870];

    auto g_yzz_x_xy_xy = buffer_fpdd[871];

    auto g_yzz_x_xy_xz = buffer_fpdd[872];

    auto g_yzz_x_xy_yy = buffer_fpdd[873];

    auto g_yzz_x_xy_yz = buffer_fpdd[874];

    auto g_yzz_x_xy_zz = buffer_fpdd[875];

    auto g_yzz_x_xz_xx = buffer_fpdd[876];

    auto g_yzz_x_xz_xy = buffer_fpdd[877];

    auto g_yzz_x_xz_xz = buffer_fpdd[878];

    auto g_yzz_x_xz_yy = buffer_fpdd[879];

    auto g_yzz_x_xz_yz = buffer_fpdd[880];

    auto g_yzz_x_xz_zz = buffer_fpdd[881];

    auto g_yzz_x_yy_xx = buffer_fpdd[882];

    auto g_yzz_x_yy_xy = buffer_fpdd[883];

    auto g_yzz_x_yy_xz = buffer_fpdd[884];

    auto g_yzz_x_yy_yy = buffer_fpdd[885];

    auto g_yzz_x_yy_yz = buffer_fpdd[886];

    auto g_yzz_x_yy_zz = buffer_fpdd[887];

    auto g_yzz_x_yz_xx = buffer_fpdd[888];

    auto g_yzz_x_yz_xy = buffer_fpdd[889];

    auto g_yzz_x_yz_xz = buffer_fpdd[890];

    auto g_yzz_x_yz_yy = buffer_fpdd[891];

    auto g_yzz_x_yz_yz = buffer_fpdd[892];

    auto g_yzz_x_yz_zz = buffer_fpdd[893];

    auto g_yzz_x_zz_xx = buffer_fpdd[894];

    auto g_yzz_x_zz_xy = buffer_fpdd[895];

    auto g_yzz_x_zz_xz = buffer_fpdd[896];

    auto g_yzz_x_zz_yy = buffer_fpdd[897];

    auto g_yzz_x_zz_yz = buffer_fpdd[898];

    auto g_yzz_x_zz_zz = buffer_fpdd[899];

    auto g_yzz_y_xx_xx = buffer_fpdd[900];

    auto g_yzz_y_xx_xy = buffer_fpdd[901];

    auto g_yzz_y_xx_xz = buffer_fpdd[902];

    auto g_yzz_y_xx_yy = buffer_fpdd[903];

    auto g_yzz_y_xx_yz = buffer_fpdd[904];

    auto g_yzz_y_xx_zz = buffer_fpdd[905];

    auto g_yzz_y_xy_xx = buffer_fpdd[906];

    auto g_yzz_y_xy_xy = buffer_fpdd[907];

    auto g_yzz_y_xy_xz = buffer_fpdd[908];

    auto g_yzz_y_xy_yy = buffer_fpdd[909];

    auto g_yzz_y_xy_yz = buffer_fpdd[910];

    auto g_yzz_y_xy_zz = buffer_fpdd[911];

    auto g_yzz_y_xz_xx = buffer_fpdd[912];

    auto g_yzz_y_xz_xy = buffer_fpdd[913];

    auto g_yzz_y_xz_xz = buffer_fpdd[914];

    auto g_yzz_y_xz_yy = buffer_fpdd[915];

    auto g_yzz_y_xz_yz = buffer_fpdd[916];

    auto g_yzz_y_xz_zz = buffer_fpdd[917];

    auto g_yzz_y_yy_xx = buffer_fpdd[918];

    auto g_yzz_y_yy_xy = buffer_fpdd[919];

    auto g_yzz_y_yy_xz = buffer_fpdd[920];

    auto g_yzz_y_yy_yy = buffer_fpdd[921];

    auto g_yzz_y_yy_yz = buffer_fpdd[922];

    auto g_yzz_y_yy_zz = buffer_fpdd[923];

    auto g_yzz_y_yz_xx = buffer_fpdd[924];

    auto g_yzz_y_yz_xy = buffer_fpdd[925];

    auto g_yzz_y_yz_xz = buffer_fpdd[926];

    auto g_yzz_y_yz_yy = buffer_fpdd[927];

    auto g_yzz_y_yz_yz = buffer_fpdd[928];

    auto g_yzz_y_yz_zz = buffer_fpdd[929];

    auto g_yzz_y_zz_xx = buffer_fpdd[930];

    auto g_yzz_y_zz_xy = buffer_fpdd[931];

    auto g_yzz_y_zz_xz = buffer_fpdd[932];

    auto g_yzz_y_zz_yy = buffer_fpdd[933];

    auto g_yzz_y_zz_yz = buffer_fpdd[934];

    auto g_yzz_y_zz_zz = buffer_fpdd[935];

    auto g_yzz_z_xx_xx = buffer_fpdd[936];

    auto g_yzz_z_xx_xy = buffer_fpdd[937];

    auto g_yzz_z_xx_xz = buffer_fpdd[938];

    auto g_yzz_z_xx_yy = buffer_fpdd[939];

    auto g_yzz_z_xx_yz = buffer_fpdd[940];

    auto g_yzz_z_xx_zz = buffer_fpdd[941];

    auto g_yzz_z_xy_xx = buffer_fpdd[942];

    auto g_yzz_z_xy_xy = buffer_fpdd[943];

    auto g_yzz_z_xy_xz = buffer_fpdd[944];

    auto g_yzz_z_xy_yy = buffer_fpdd[945];

    auto g_yzz_z_xy_yz = buffer_fpdd[946];

    auto g_yzz_z_xy_zz = buffer_fpdd[947];

    auto g_yzz_z_xz_xx = buffer_fpdd[948];

    auto g_yzz_z_xz_xy = buffer_fpdd[949];

    auto g_yzz_z_xz_xz = buffer_fpdd[950];

    auto g_yzz_z_xz_yy = buffer_fpdd[951];

    auto g_yzz_z_xz_yz = buffer_fpdd[952];

    auto g_yzz_z_xz_zz = buffer_fpdd[953];

    auto g_yzz_z_yy_xx = buffer_fpdd[954];

    auto g_yzz_z_yy_xy = buffer_fpdd[955];

    auto g_yzz_z_yy_xz = buffer_fpdd[956];

    auto g_yzz_z_yy_yy = buffer_fpdd[957];

    auto g_yzz_z_yy_yz = buffer_fpdd[958];

    auto g_yzz_z_yy_zz = buffer_fpdd[959];

    auto g_yzz_z_yz_xx = buffer_fpdd[960];

    auto g_yzz_z_yz_xy = buffer_fpdd[961];

    auto g_yzz_z_yz_xz = buffer_fpdd[962];

    auto g_yzz_z_yz_yy = buffer_fpdd[963];

    auto g_yzz_z_yz_yz = buffer_fpdd[964];

    auto g_yzz_z_yz_zz = buffer_fpdd[965];

    auto g_yzz_z_zz_xx = buffer_fpdd[966];

    auto g_yzz_z_zz_xy = buffer_fpdd[967];

    auto g_yzz_z_zz_xz = buffer_fpdd[968];

    auto g_yzz_z_zz_yy = buffer_fpdd[969];

    auto g_yzz_z_zz_yz = buffer_fpdd[970];

    auto g_yzz_z_zz_zz = buffer_fpdd[971];

    auto g_zzz_x_xx_xx = buffer_fpdd[972];

    auto g_zzz_x_xx_xy = buffer_fpdd[973];

    auto g_zzz_x_xx_xz = buffer_fpdd[974];

    auto g_zzz_x_xx_yy = buffer_fpdd[975];

    auto g_zzz_x_xx_yz = buffer_fpdd[976];

    auto g_zzz_x_xx_zz = buffer_fpdd[977];

    auto g_zzz_x_xy_xx = buffer_fpdd[978];

    auto g_zzz_x_xy_xy = buffer_fpdd[979];

    auto g_zzz_x_xy_xz = buffer_fpdd[980];

    auto g_zzz_x_xy_yy = buffer_fpdd[981];

    auto g_zzz_x_xy_yz = buffer_fpdd[982];

    auto g_zzz_x_xy_zz = buffer_fpdd[983];

    auto g_zzz_x_xz_xx = buffer_fpdd[984];

    auto g_zzz_x_xz_xy = buffer_fpdd[985];

    auto g_zzz_x_xz_xz = buffer_fpdd[986];

    auto g_zzz_x_xz_yy = buffer_fpdd[987];

    auto g_zzz_x_xz_yz = buffer_fpdd[988];

    auto g_zzz_x_xz_zz = buffer_fpdd[989];

    auto g_zzz_x_yy_xx = buffer_fpdd[990];

    auto g_zzz_x_yy_xy = buffer_fpdd[991];

    auto g_zzz_x_yy_xz = buffer_fpdd[992];

    auto g_zzz_x_yy_yy = buffer_fpdd[993];

    auto g_zzz_x_yy_yz = buffer_fpdd[994];

    auto g_zzz_x_yy_zz = buffer_fpdd[995];

    auto g_zzz_x_yz_xx = buffer_fpdd[996];

    auto g_zzz_x_yz_xy = buffer_fpdd[997];

    auto g_zzz_x_yz_xz = buffer_fpdd[998];

    auto g_zzz_x_yz_yy = buffer_fpdd[999];

    auto g_zzz_x_yz_yz = buffer_fpdd[1000];

    auto g_zzz_x_yz_zz = buffer_fpdd[1001];

    auto g_zzz_x_zz_xx = buffer_fpdd[1002];

    auto g_zzz_x_zz_xy = buffer_fpdd[1003];

    auto g_zzz_x_zz_xz = buffer_fpdd[1004];

    auto g_zzz_x_zz_yy = buffer_fpdd[1005];

    auto g_zzz_x_zz_yz = buffer_fpdd[1006];

    auto g_zzz_x_zz_zz = buffer_fpdd[1007];

    auto g_zzz_y_xx_xx = buffer_fpdd[1008];

    auto g_zzz_y_xx_xy = buffer_fpdd[1009];

    auto g_zzz_y_xx_xz = buffer_fpdd[1010];

    auto g_zzz_y_xx_yy = buffer_fpdd[1011];

    auto g_zzz_y_xx_yz = buffer_fpdd[1012];

    auto g_zzz_y_xx_zz = buffer_fpdd[1013];

    auto g_zzz_y_xy_xx = buffer_fpdd[1014];

    auto g_zzz_y_xy_xy = buffer_fpdd[1015];

    auto g_zzz_y_xy_xz = buffer_fpdd[1016];

    auto g_zzz_y_xy_yy = buffer_fpdd[1017];

    auto g_zzz_y_xy_yz = buffer_fpdd[1018];

    auto g_zzz_y_xy_zz = buffer_fpdd[1019];

    auto g_zzz_y_xz_xx = buffer_fpdd[1020];

    auto g_zzz_y_xz_xy = buffer_fpdd[1021];

    auto g_zzz_y_xz_xz = buffer_fpdd[1022];

    auto g_zzz_y_xz_yy = buffer_fpdd[1023];

    auto g_zzz_y_xz_yz = buffer_fpdd[1024];

    auto g_zzz_y_xz_zz = buffer_fpdd[1025];

    auto g_zzz_y_yy_xx = buffer_fpdd[1026];

    auto g_zzz_y_yy_xy = buffer_fpdd[1027];

    auto g_zzz_y_yy_xz = buffer_fpdd[1028];

    auto g_zzz_y_yy_yy = buffer_fpdd[1029];

    auto g_zzz_y_yy_yz = buffer_fpdd[1030];

    auto g_zzz_y_yy_zz = buffer_fpdd[1031];

    auto g_zzz_y_yz_xx = buffer_fpdd[1032];

    auto g_zzz_y_yz_xy = buffer_fpdd[1033];

    auto g_zzz_y_yz_xz = buffer_fpdd[1034];

    auto g_zzz_y_yz_yy = buffer_fpdd[1035];

    auto g_zzz_y_yz_yz = buffer_fpdd[1036];

    auto g_zzz_y_yz_zz = buffer_fpdd[1037];

    auto g_zzz_y_zz_xx = buffer_fpdd[1038];

    auto g_zzz_y_zz_xy = buffer_fpdd[1039];

    auto g_zzz_y_zz_xz = buffer_fpdd[1040];

    auto g_zzz_y_zz_yy = buffer_fpdd[1041];

    auto g_zzz_y_zz_yz = buffer_fpdd[1042];

    auto g_zzz_y_zz_zz = buffer_fpdd[1043];

    auto g_zzz_z_xx_xx = buffer_fpdd[1044];

    auto g_zzz_z_xx_xy = buffer_fpdd[1045];

    auto g_zzz_z_xx_xz = buffer_fpdd[1046];

    auto g_zzz_z_xx_yy = buffer_fpdd[1047];

    auto g_zzz_z_xx_yz = buffer_fpdd[1048];

    auto g_zzz_z_xx_zz = buffer_fpdd[1049];

    auto g_zzz_z_xy_xx = buffer_fpdd[1050];

    auto g_zzz_z_xy_xy = buffer_fpdd[1051];

    auto g_zzz_z_xy_xz = buffer_fpdd[1052];

    auto g_zzz_z_xy_yy = buffer_fpdd[1053];

    auto g_zzz_z_xy_yz = buffer_fpdd[1054];

    auto g_zzz_z_xy_zz = buffer_fpdd[1055];

    auto g_zzz_z_xz_xx = buffer_fpdd[1056];

    auto g_zzz_z_xz_xy = buffer_fpdd[1057];

    auto g_zzz_z_xz_xz = buffer_fpdd[1058];

    auto g_zzz_z_xz_yy = buffer_fpdd[1059];

    auto g_zzz_z_xz_yz = buffer_fpdd[1060];

    auto g_zzz_z_xz_zz = buffer_fpdd[1061];

    auto g_zzz_z_yy_xx = buffer_fpdd[1062];

    auto g_zzz_z_yy_xy = buffer_fpdd[1063];

    auto g_zzz_z_yy_xz = buffer_fpdd[1064];

    auto g_zzz_z_yy_yy = buffer_fpdd[1065];

    auto g_zzz_z_yy_yz = buffer_fpdd[1066];

    auto g_zzz_z_yy_zz = buffer_fpdd[1067];

    auto g_zzz_z_yz_xx = buffer_fpdd[1068];

    auto g_zzz_z_yz_xy = buffer_fpdd[1069];

    auto g_zzz_z_yz_xz = buffer_fpdd[1070];

    auto g_zzz_z_yz_yy = buffer_fpdd[1071];

    auto g_zzz_z_yz_yz = buffer_fpdd[1072];

    auto g_zzz_z_yz_zz = buffer_fpdd[1073];

    auto g_zzz_z_zz_xx = buffer_fpdd[1074];

    auto g_zzz_z_zz_xy = buffer_fpdd[1075];

    auto g_zzz_z_zz_xz = buffer_fpdd[1076];

    auto g_zzz_z_zz_yy = buffer_fpdd[1077];

    auto g_zzz_z_zz_yz = buffer_fpdd[1078];

    auto g_zzz_z_zz_zz = buffer_fpdd[1079];

    /// Set up components of integrals buffer : buffer_2000_ppdd

    auto g_xx_0_0_0_x_x_xx_xx = buffer_2000_ppdd[0];

    auto g_xx_0_0_0_x_x_xx_xy = buffer_2000_ppdd[1];

    auto g_xx_0_0_0_x_x_xx_xz = buffer_2000_ppdd[2];

    auto g_xx_0_0_0_x_x_xx_yy = buffer_2000_ppdd[3];

    auto g_xx_0_0_0_x_x_xx_yz = buffer_2000_ppdd[4];

    auto g_xx_0_0_0_x_x_xx_zz = buffer_2000_ppdd[5];

    auto g_xx_0_0_0_x_x_xy_xx = buffer_2000_ppdd[6];

    auto g_xx_0_0_0_x_x_xy_xy = buffer_2000_ppdd[7];

    auto g_xx_0_0_0_x_x_xy_xz = buffer_2000_ppdd[8];

    auto g_xx_0_0_0_x_x_xy_yy = buffer_2000_ppdd[9];

    auto g_xx_0_0_0_x_x_xy_yz = buffer_2000_ppdd[10];

    auto g_xx_0_0_0_x_x_xy_zz = buffer_2000_ppdd[11];

    auto g_xx_0_0_0_x_x_xz_xx = buffer_2000_ppdd[12];

    auto g_xx_0_0_0_x_x_xz_xy = buffer_2000_ppdd[13];

    auto g_xx_0_0_0_x_x_xz_xz = buffer_2000_ppdd[14];

    auto g_xx_0_0_0_x_x_xz_yy = buffer_2000_ppdd[15];

    auto g_xx_0_0_0_x_x_xz_yz = buffer_2000_ppdd[16];

    auto g_xx_0_0_0_x_x_xz_zz = buffer_2000_ppdd[17];

    auto g_xx_0_0_0_x_x_yy_xx = buffer_2000_ppdd[18];

    auto g_xx_0_0_0_x_x_yy_xy = buffer_2000_ppdd[19];

    auto g_xx_0_0_0_x_x_yy_xz = buffer_2000_ppdd[20];

    auto g_xx_0_0_0_x_x_yy_yy = buffer_2000_ppdd[21];

    auto g_xx_0_0_0_x_x_yy_yz = buffer_2000_ppdd[22];

    auto g_xx_0_0_0_x_x_yy_zz = buffer_2000_ppdd[23];

    auto g_xx_0_0_0_x_x_yz_xx = buffer_2000_ppdd[24];

    auto g_xx_0_0_0_x_x_yz_xy = buffer_2000_ppdd[25];

    auto g_xx_0_0_0_x_x_yz_xz = buffer_2000_ppdd[26];

    auto g_xx_0_0_0_x_x_yz_yy = buffer_2000_ppdd[27];

    auto g_xx_0_0_0_x_x_yz_yz = buffer_2000_ppdd[28];

    auto g_xx_0_0_0_x_x_yz_zz = buffer_2000_ppdd[29];

    auto g_xx_0_0_0_x_x_zz_xx = buffer_2000_ppdd[30];

    auto g_xx_0_0_0_x_x_zz_xy = buffer_2000_ppdd[31];

    auto g_xx_0_0_0_x_x_zz_xz = buffer_2000_ppdd[32];

    auto g_xx_0_0_0_x_x_zz_yy = buffer_2000_ppdd[33];

    auto g_xx_0_0_0_x_x_zz_yz = buffer_2000_ppdd[34];

    auto g_xx_0_0_0_x_x_zz_zz = buffer_2000_ppdd[35];

    auto g_xx_0_0_0_x_y_xx_xx = buffer_2000_ppdd[36];

    auto g_xx_0_0_0_x_y_xx_xy = buffer_2000_ppdd[37];

    auto g_xx_0_0_0_x_y_xx_xz = buffer_2000_ppdd[38];

    auto g_xx_0_0_0_x_y_xx_yy = buffer_2000_ppdd[39];

    auto g_xx_0_0_0_x_y_xx_yz = buffer_2000_ppdd[40];

    auto g_xx_0_0_0_x_y_xx_zz = buffer_2000_ppdd[41];

    auto g_xx_0_0_0_x_y_xy_xx = buffer_2000_ppdd[42];

    auto g_xx_0_0_0_x_y_xy_xy = buffer_2000_ppdd[43];

    auto g_xx_0_0_0_x_y_xy_xz = buffer_2000_ppdd[44];

    auto g_xx_0_0_0_x_y_xy_yy = buffer_2000_ppdd[45];

    auto g_xx_0_0_0_x_y_xy_yz = buffer_2000_ppdd[46];

    auto g_xx_0_0_0_x_y_xy_zz = buffer_2000_ppdd[47];

    auto g_xx_0_0_0_x_y_xz_xx = buffer_2000_ppdd[48];

    auto g_xx_0_0_0_x_y_xz_xy = buffer_2000_ppdd[49];

    auto g_xx_0_0_0_x_y_xz_xz = buffer_2000_ppdd[50];

    auto g_xx_0_0_0_x_y_xz_yy = buffer_2000_ppdd[51];

    auto g_xx_0_0_0_x_y_xz_yz = buffer_2000_ppdd[52];

    auto g_xx_0_0_0_x_y_xz_zz = buffer_2000_ppdd[53];

    auto g_xx_0_0_0_x_y_yy_xx = buffer_2000_ppdd[54];

    auto g_xx_0_0_0_x_y_yy_xy = buffer_2000_ppdd[55];

    auto g_xx_0_0_0_x_y_yy_xz = buffer_2000_ppdd[56];

    auto g_xx_0_0_0_x_y_yy_yy = buffer_2000_ppdd[57];

    auto g_xx_0_0_0_x_y_yy_yz = buffer_2000_ppdd[58];

    auto g_xx_0_0_0_x_y_yy_zz = buffer_2000_ppdd[59];

    auto g_xx_0_0_0_x_y_yz_xx = buffer_2000_ppdd[60];

    auto g_xx_0_0_0_x_y_yz_xy = buffer_2000_ppdd[61];

    auto g_xx_0_0_0_x_y_yz_xz = buffer_2000_ppdd[62];

    auto g_xx_0_0_0_x_y_yz_yy = buffer_2000_ppdd[63];

    auto g_xx_0_0_0_x_y_yz_yz = buffer_2000_ppdd[64];

    auto g_xx_0_0_0_x_y_yz_zz = buffer_2000_ppdd[65];

    auto g_xx_0_0_0_x_y_zz_xx = buffer_2000_ppdd[66];

    auto g_xx_0_0_0_x_y_zz_xy = buffer_2000_ppdd[67];

    auto g_xx_0_0_0_x_y_zz_xz = buffer_2000_ppdd[68];

    auto g_xx_0_0_0_x_y_zz_yy = buffer_2000_ppdd[69];

    auto g_xx_0_0_0_x_y_zz_yz = buffer_2000_ppdd[70];

    auto g_xx_0_0_0_x_y_zz_zz = buffer_2000_ppdd[71];

    auto g_xx_0_0_0_x_z_xx_xx = buffer_2000_ppdd[72];

    auto g_xx_0_0_0_x_z_xx_xy = buffer_2000_ppdd[73];

    auto g_xx_0_0_0_x_z_xx_xz = buffer_2000_ppdd[74];

    auto g_xx_0_0_0_x_z_xx_yy = buffer_2000_ppdd[75];

    auto g_xx_0_0_0_x_z_xx_yz = buffer_2000_ppdd[76];

    auto g_xx_0_0_0_x_z_xx_zz = buffer_2000_ppdd[77];

    auto g_xx_0_0_0_x_z_xy_xx = buffer_2000_ppdd[78];

    auto g_xx_0_0_0_x_z_xy_xy = buffer_2000_ppdd[79];

    auto g_xx_0_0_0_x_z_xy_xz = buffer_2000_ppdd[80];

    auto g_xx_0_0_0_x_z_xy_yy = buffer_2000_ppdd[81];

    auto g_xx_0_0_0_x_z_xy_yz = buffer_2000_ppdd[82];

    auto g_xx_0_0_0_x_z_xy_zz = buffer_2000_ppdd[83];

    auto g_xx_0_0_0_x_z_xz_xx = buffer_2000_ppdd[84];

    auto g_xx_0_0_0_x_z_xz_xy = buffer_2000_ppdd[85];

    auto g_xx_0_0_0_x_z_xz_xz = buffer_2000_ppdd[86];

    auto g_xx_0_0_0_x_z_xz_yy = buffer_2000_ppdd[87];

    auto g_xx_0_0_0_x_z_xz_yz = buffer_2000_ppdd[88];

    auto g_xx_0_0_0_x_z_xz_zz = buffer_2000_ppdd[89];

    auto g_xx_0_0_0_x_z_yy_xx = buffer_2000_ppdd[90];

    auto g_xx_0_0_0_x_z_yy_xy = buffer_2000_ppdd[91];

    auto g_xx_0_0_0_x_z_yy_xz = buffer_2000_ppdd[92];

    auto g_xx_0_0_0_x_z_yy_yy = buffer_2000_ppdd[93];

    auto g_xx_0_0_0_x_z_yy_yz = buffer_2000_ppdd[94];

    auto g_xx_0_0_0_x_z_yy_zz = buffer_2000_ppdd[95];

    auto g_xx_0_0_0_x_z_yz_xx = buffer_2000_ppdd[96];

    auto g_xx_0_0_0_x_z_yz_xy = buffer_2000_ppdd[97];

    auto g_xx_0_0_0_x_z_yz_xz = buffer_2000_ppdd[98];

    auto g_xx_0_0_0_x_z_yz_yy = buffer_2000_ppdd[99];

    auto g_xx_0_0_0_x_z_yz_yz = buffer_2000_ppdd[100];

    auto g_xx_0_0_0_x_z_yz_zz = buffer_2000_ppdd[101];

    auto g_xx_0_0_0_x_z_zz_xx = buffer_2000_ppdd[102];

    auto g_xx_0_0_0_x_z_zz_xy = buffer_2000_ppdd[103];

    auto g_xx_0_0_0_x_z_zz_xz = buffer_2000_ppdd[104];

    auto g_xx_0_0_0_x_z_zz_yy = buffer_2000_ppdd[105];

    auto g_xx_0_0_0_x_z_zz_yz = buffer_2000_ppdd[106];

    auto g_xx_0_0_0_x_z_zz_zz = buffer_2000_ppdd[107];

    auto g_xx_0_0_0_y_x_xx_xx = buffer_2000_ppdd[108];

    auto g_xx_0_0_0_y_x_xx_xy = buffer_2000_ppdd[109];

    auto g_xx_0_0_0_y_x_xx_xz = buffer_2000_ppdd[110];

    auto g_xx_0_0_0_y_x_xx_yy = buffer_2000_ppdd[111];

    auto g_xx_0_0_0_y_x_xx_yz = buffer_2000_ppdd[112];

    auto g_xx_0_0_0_y_x_xx_zz = buffer_2000_ppdd[113];

    auto g_xx_0_0_0_y_x_xy_xx = buffer_2000_ppdd[114];

    auto g_xx_0_0_0_y_x_xy_xy = buffer_2000_ppdd[115];

    auto g_xx_0_0_0_y_x_xy_xz = buffer_2000_ppdd[116];

    auto g_xx_0_0_0_y_x_xy_yy = buffer_2000_ppdd[117];

    auto g_xx_0_0_0_y_x_xy_yz = buffer_2000_ppdd[118];

    auto g_xx_0_0_0_y_x_xy_zz = buffer_2000_ppdd[119];

    auto g_xx_0_0_0_y_x_xz_xx = buffer_2000_ppdd[120];

    auto g_xx_0_0_0_y_x_xz_xy = buffer_2000_ppdd[121];

    auto g_xx_0_0_0_y_x_xz_xz = buffer_2000_ppdd[122];

    auto g_xx_0_0_0_y_x_xz_yy = buffer_2000_ppdd[123];

    auto g_xx_0_0_0_y_x_xz_yz = buffer_2000_ppdd[124];

    auto g_xx_0_0_0_y_x_xz_zz = buffer_2000_ppdd[125];

    auto g_xx_0_0_0_y_x_yy_xx = buffer_2000_ppdd[126];

    auto g_xx_0_0_0_y_x_yy_xy = buffer_2000_ppdd[127];

    auto g_xx_0_0_0_y_x_yy_xz = buffer_2000_ppdd[128];

    auto g_xx_0_0_0_y_x_yy_yy = buffer_2000_ppdd[129];

    auto g_xx_0_0_0_y_x_yy_yz = buffer_2000_ppdd[130];

    auto g_xx_0_0_0_y_x_yy_zz = buffer_2000_ppdd[131];

    auto g_xx_0_0_0_y_x_yz_xx = buffer_2000_ppdd[132];

    auto g_xx_0_0_0_y_x_yz_xy = buffer_2000_ppdd[133];

    auto g_xx_0_0_0_y_x_yz_xz = buffer_2000_ppdd[134];

    auto g_xx_0_0_0_y_x_yz_yy = buffer_2000_ppdd[135];

    auto g_xx_0_0_0_y_x_yz_yz = buffer_2000_ppdd[136];

    auto g_xx_0_0_0_y_x_yz_zz = buffer_2000_ppdd[137];

    auto g_xx_0_0_0_y_x_zz_xx = buffer_2000_ppdd[138];

    auto g_xx_0_0_0_y_x_zz_xy = buffer_2000_ppdd[139];

    auto g_xx_0_0_0_y_x_zz_xz = buffer_2000_ppdd[140];

    auto g_xx_0_0_0_y_x_zz_yy = buffer_2000_ppdd[141];

    auto g_xx_0_0_0_y_x_zz_yz = buffer_2000_ppdd[142];

    auto g_xx_0_0_0_y_x_zz_zz = buffer_2000_ppdd[143];

    auto g_xx_0_0_0_y_y_xx_xx = buffer_2000_ppdd[144];

    auto g_xx_0_0_0_y_y_xx_xy = buffer_2000_ppdd[145];

    auto g_xx_0_0_0_y_y_xx_xz = buffer_2000_ppdd[146];

    auto g_xx_0_0_0_y_y_xx_yy = buffer_2000_ppdd[147];

    auto g_xx_0_0_0_y_y_xx_yz = buffer_2000_ppdd[148];

    auto g_xx_0_0_0_y_y_xx_zz = buffer_2000_ppdd[149];

    auto g_xx_0_0_0_y_y_xy_xx = buffer_2000_ppdd[150];

    auto g_xx_0_0_0_y_y_xy_xy = buffer_2000_ppdd[151];

    auto g_xx_0_0_0_y_y_xy_xz = buffer_2000_ppdd[152];

    auto g_xx_0_0_0_y_y_xy_yy = buffer_2000_ppdd[153];

    auto g_xx_0_0_0_y_y_xy_yz = buffer_2000_ppdd[154];

    auto g_xx_0_0_0_y_y_xy_zz = buffer_2000_ppdd[155];

    auto g_xx_0_0_0_y_y_xz_xx = buffer_2000_ppdd[156];

    auto g_xx_0_0_0_y_y_xz_xy = buffer_2000_ppdd[157];

    auto g_xx_0_0_0_y_y_xz_xz = buffer_2000_ppdd[158];

    auto g_xx_0_0_0_y_y_xz_yy = buffer_2000_ppdd[159];

    auto g_xx_0_0_0_y_y_xz_yz = buffer_2000_ppdd[160];

    auto g_xx_0_0_0_y_y_xz_zz = buffer_2000_ppdd[161];

    auto g_xx_0_0_0_y_y_yy_xx = buffer_2000_ppdd[162];

    auto g_xx_0_0_0_y_y_yy_xy = buffer_2000_ppdd[163];

    auto g_xx_0_0_0_y_y_yy_xz = buffer_2000_ppdd[164];

    auto g_xx_0_0_0_y_y_yy_yy = buffer_2000_ppdd[165];

    auto g_xx_0_0_0_y_y_yy_yz = buffer_2000_ppdd[166];

    auto g_xx_0_0_0_y_y_yy_zz = buffer_2000_ppdd[167];

    auto g_xx_0_0_0_y_y_yz_xx = buffer_2000_ppdd[168];

    auto g_xx_0_0_0_y_y_yz_xy = buffer_2000_ppdd[169];

    auto g_xx_0_0_0_y_y_yz_xz = buffer_2000_ppdd[170];

    auto g_xx_0_0_0_y_y_yz_yy = buffer_2000_ppdd[171];

    auto g_xx_0_0_0_y_y_yz_yz = buffer_2000_ppdd[172];

    auto g_xx_0_0_0_y_y_yz_zz = buffer_2000_ppdd[173];

    auto g_xx_0_0_0_y_y_zz_xx = buffer_2000_ppdd[174];

    auto g_xx_0_0_0_y_y_zz_xy = buffer_2000_ppdd[175];

    auto g_xx_0_0_0_y_y_zz_xz = buffer_2000_ppdd[176];

    auto g_xx_0_0_0_y_y_zz_yy = buffer_2000_ppdd[177];

    auto g_xx_0_0_0_y_y_zz_yz = buffer_2000_ppdd[178];

    auto g_xx_0_0_0_y_y_zz_zz = buffer_2000_ppdd[179];

    auto g_xx_0_0_0_y_z_xx_xx = buffer_2000_ppdd[180];

    auto g_xx_0_0_0_y_z_xx_xy = buffer_2000_ppdd[181];

    auto g_xx_0_0_0_y_z_xx_xz = buffer_2000_ppdd[182];

    auto g_xx_0_0_0_y_z_xx_yy = buffer_2000_ppdd[183];

    auto g_xx_0_0_0_y_z_xx_yz = buffer_2000_ppdd[184];

    auto g_xx_0_0_0_y_z_xx_zz = buffer_2000_ppdd[185];

    auto g_xx_0_0_0_y_z_xy_xx = buffer_2000_ppdd[186];

    auto g_xx_0_0_0_y_z_xy_xy = buffer_2000_ppdd[187];

    auto g_xx_0_0_0_y_z_xy_xz = buffer_2000_ppdd[188];

    auto g_xx_0_0_0_y_z_xy_yy = buffer_2000_ppdd[189];

    auto g_xx_0_0_0_y_z_xy_yz = buffer_2000_ppdd[190];

    auto g_xx_0_0_0_y_z_xy_zz = buffer_2000_ppdd[191];

    auto g_xx_0_0_0_y_z_xz_xx = buffer_2000_ppdd[192];

    auto g_xx_0_0_0_y_z_xz_xy = buffer_2000_ppdd[193];

    auto g_xx_0_0_0_y_z_xz_xz = buffer_2000_ppdd[194];

    auto g_xx_0_0_0_y_z_xz_yy = buffer_2000_ppdd[195];

    auto g_xx_0_0_0_y_z_xz_yz = buffer_2000_ppdd[196];

    auto g_xx_0_0_0_y_z_xz_zz = buffer_2000_ppdd[197];

    auto g_xx_0_0_0_y_z_yy_xx = buffer_2000_ppdd[198];

    auto g_xx_0_0_0_y_z_yy_xy = buffer_2000_ppdd[199];

    auto g_xx_0_0_0_y_z_yy_xz = buffer_2000_ppdd[200];

    auto g_xx_0_0_0_y_z_yy_yy = buffer_2000_ppdd[201];

    auto g_xx_0_0_0_y_z_yy_yz = buffer_2000_ppdd[202];

    auto g_xx_0_0_0_y_z_yy_zz = buffer_2000_ppdd[203];

    auto g_xx_0_0_0_y_z_yz_xx = buffer_2000_ppdd[204];

    auto g_xx_0_0_0_y_z_yz_xy = buffer_2000_ppdd[205];

    auto g_xx_0_0_0_y_z_yz_xz = buffer_2000_ppdd[206];

    auto g_xx_0_0_0_y_z_yz_yy = buffer_2000_ppdd[207];

    auto g_xx_0_0_0_y_z_yz_yz = buffer_2000_ppdd[208];

    auto g_xx_0_0_0_y_z_yz_zz = buffer_2000_ppdd[209];

    auto g_xx_0_0_0_y_z_zz_xx = buffer_2000_ppdd[210];

    auto g_xx_0_0_0_y_z_zz_xy = buffer_2000_ppdd[211];

    auto g_xx_0_0_0_y_z_zz_xz = buffer_2000_ppdd[212];

    auto g_xx_0_0_0_y_z_zz_yy = buffer_2000_ppdd[213];

    auto g_xx_0_0_0_y_z_zz_yz = buffer_2000_ppdd[214];

    auto g_xx_0_0_0_y_z_zz_zz = buffer_2000_ppdd[215];

    auto g_xx_0_0_0_z_x_xx_xx = buffer_2000_ppdd[216];

    auto g_xx_0_0_0_z_x_xx_xy = buffer_2000_ppdd[217];

    auto g_xx_0_0_0_z_x_xx_xz = buffer_2000_ppdd[218];

    auto g_xx_0_0_0_z_x_xx_yy = buffer_2000_ppdd[219];

    auto g_xx_0_0_0_z_x_xx_yz = buffer_2000_ppdd[220];

    auto g_xx_0_0_0_z_x_xx_zz = buffer_2000_ppdd[221];

    auto g_xx_0_0_0_z_x_xy_xx = buffer_2000_ppdd[222];

    auto g_xx_0_0_0_z_x_xy_xy = buffer_2000_ppdd[223];

    auto g_xx_0_0_0_z_x_xy_xz = buffer_2000_ppdd[224];

    auto g_xx_0_0_0_z_x_xy_yy = buffer_2000_ppdd[225];

    auto g_xx_0_0_0_z_x_xy_yz = buffer_2000_ppdd[226];

    auto g_xx_0_0_0_z_x_xy_zz = buffer_2000_ppdd[227];

    auto g_xx_0_0_0_z_x_xz_xx = buffer_2000_ppdd[228];

    auto g_xx_0_0_0_z_x_xz_xy = buffer_2000_ppdd[229];

    auto g_xx_0_0_0_z_x_xz_xz = buffer_2000_ppdd[230];

    auto g_xx_0_0_0_z_x_xz_yy = buffer_2000_ppdd[231];

    auto g_xx_0_0_0_z_x_xz_yz = buffer_2000_ppdd[232];

    auto g_xx_0_0_0_z_x_xz_zz = buffer_2000_ppdd[233];

    auto g_xx_0_0_0_z_x_yy_xx = buffer_2000_ppdd[234];

    auto g_xx_0_0_0_z_x_yy_xy = buffer_2000_ppdd[235];

    auto g_xx_0_0_0_z_x_yy_xz = buffer_2000_ppdd[236];

    auto g_xx_0_0_0_z_x_yy_yy = buffer_2000_ppdd[237];

    auto g_xx_0_0_0_z_x_yy_yz = buffer_2000_ppdd[238];

    auto g_xx_0_0_0_z_x_yy_zz = buffer_2000_ppdd[239];

    auto g_xx_0_0_0_z_x_yz_xx = buffer_2000_ppdd[240];

    auto g_xx_0_0_0_z_x_yz_xy = buffer_2000_ppdd[241];

    auto g_xx_0_0_0_z_x_yz_xz = buffer_2000_ppdd[242];

    auto g_xx_0_0_0_z_x_yz_yy = buffer_2000_ppdd[243];

    auto g_xx_0_0_0_z_x_yz_yz = buffer_2000_ppdd[244];

    auto g_xx_0_0_0_z_x_yz_zz = buffer_2000_ppdd[245];

    auto g_xx_0_0_0_z_x_zz_xx = buffer_2000_ppdd[246];

    auto g_xx_0_0_0_z_x_zz_xy = buffer_2000_ppdd[247];

    auto g_xx_0_0_0_z_x_zz_xz = buffer_2000_ppdd[248];

    auto g_xx_0_0_0_z_x_zz_yy = buffer_2000_ppdd[249];

    auto g_xx_0_0_0_z_x_zz_yz = buffer_2000_ppdd[250];

    auto g_xx_0_0_0_z_x_zz_zz = buffer_2000_ppdd[251];

    auto g_xx_0_0_0_z_y_xx_xx = buffer_2000_ppdd[252];

    auto g_xx_0_0_0_z_y_xx_xy = buffer_2000_ppdd[253];

    auto g_xx_0_0_0_z_y_xx_xz = buffer_2000_ppdd[254];

    auto g_xx_0_0_0_z_y_xx_yy = buffer_2000_ppdd[255];

    auto g_xx_0_0_0_z_y_xx_yz = buffer_2000_ppdd[256];

    auto g_xx_0_0_0_z_y_xx_zz = buffer_2000_ppdd[257];

    auto g_xx_0_0_0_z_y_xy_xx = buffer_2000_ppdd[258];

    auto g_xx_0_0_0_z_y_xy_xy = buffer_2000_ppdd[259];

    auto g_xx_0_0_0_z_y_xy_xz = buffer_2000_ppdd[260];

    auto g_xx_0_0_0_z_y_xy_yy = buffer_2000_ppdd[261];

    auto g_xx_0_0_0_z_y_xy_yz = buffer_2000_ppdd[262];

    auto g_xx_0_0_0_z_y_xy_zz = buffer_2000_ppdd[263];

    auto g_xx_0_0_0_z_y_xz_xx = buffer_2000_ppdd[264];

    auto g_xx_0_0_0_z_y_xz_xy = buffer_2000_ppdd[265];

    auto g_xx_0_0_0_z_y_xz_xz = buffer_2000_ppdd[266];

    auto g_xx_0_0_0_z_y_xz_yy = buffer_2000_ppdd[267];

    auto g_xx_0_0_0_z_y_xz_yz = buffer_2000_ppdd[268];

    auto g_xx_0_0_0_z_y_xz_zz = buffer_2000_ppdd[269];

    auto g_xx_0_0_0_z_y_yy_xx = buffer_2000_ppdd[270];

    auto g_xx_0_0_0_z_y_yy_xy = buffer_2000_ppdd[271];

    auto g_xx_0_0_0_z_y_yy_xz = buffer_2000_ppdd[272];

    auto g_xx_0_0_0_z_y_yy_yy = buffer_2000_ppdd[273];

    auto g_xx_0_0_0_z_y_yy_yz = buffer_2000_ppdd[274];

    auto g_xx_0_0_0_z_y_yy_zz = buffer_2000_ppdd[275];

    auto g_xx_0_0_0_z_y_yz_xx = buffer_2000_ppdd[276];

    auto g_xx_0_0_0_z_y_yz_xy = buffer_2000_ppdd[277];

    auto g_xx_0_0_0_z_y_yz_xz = buffer_2000_ppdd[278];

    auto g_xx_0_0_0_z_y_yz_yy = buffer_2000_ppdd[279];

    auto g_xx_0_0_0_z_y_yz_yz = buffer_2000_ppdd[280];

    auto g_xx_0_0_0_z_y_yz_zz = buffer_2000_ppdd[281];

    auto g_xx_0_0_0_z_y_zz_xx = buffer_2000_ppdd[282];

    auto g_xx_0_0_0_z_y_zz_xy = buffer_2000_ppdd[283];

    auto g_xx_0_0_0_z_y_zz_xz = buffer_2000_ppdd[284];

    auto g_xx_0_0_0_z_y_zz_yy = buffer_2000_ppdd[285];

    auto g_xx_0_0_0_z_y_zz_yz = buffer_2000_ppdd[286];

    auto g_xx_0_0_0_z_y_zz_zz = buffer_2000_ppdd[287];

    auto g_xx_0_0_0_z_z_xx_xx = buffer_2000_ppdd[288];

    auto g_xx_0_0_0_z_z_xx_xy = buffer_2000_ppdd[289];

    auto g_xx_0_0_0_z_z_xx_xz = buffer_2000_ppdd[290];

    auto g_xx_0_0_0_z_z_xx_yy = buffer_2000_ppdd[291];

    auto g_xx_0_0_0_z_z_xx_yz = buffer_2000_ppdd[292];

    auto g_xx_0_0_0_z_z_xx_zz = buffer_2000_ppdd[293];

    auto g_xx_0_0_0_z_z_xy_xx = buffer_2000_ppdd[294];

    auto g_xx_0_0_0_z_z_xy_xy = buffer_2000_ppdd[295];

    auto g_xx_0_0_0_z_z_xy_xz = buffer_2000_ppdd[296];

    auto g_xx_0_0_0_z_z_xy_yy = buffer_2000_ppdd[297];

    auto g_xx_0_0_0_z_z_xy_yz = buffer_2000_ppdd[298];

    auto g_xx_0_0_0_z_z_xy_zz = buffer_2000_ppdd[299];

    auto g_xx_0_0_0_z_z_xz_xx = buffer_2000_ppdd[300];

    auto g_xx_0_0_0_z_z_xz_xy = buffer_2000_ppdd[301];

    auto g_xx_0_0_0_z_z_xz_xz = buffer_2000_ppdd[302];

    auto g_xx_0_0_0_z_z_xz_yy = buffer_2000_ppdd[303];

    auto g_xx_0_0_0_z_z_xz_yz = buffer_2000_ppdd[304];

    auto g_xx_0_0_0_z_z_xz_zz = buffer_2000_ppdd[305];

    auto g_xx_0_0_0_z_z_yy_xx = buffer_2000_ppdd[306];

    auto g_xx_0_0_0_z_z_yy_xy = buffer_2000_ppdd[307];

    auto g_xx_0_0_0_z_z_yy_xz = buffer_2000_ppdd[308];

    auto g_xx_0_0_0_z_z_yy_yy = buffer_2000_ppdd[309];

    auto g_xx_0_0_0_z_z_yy_yz = buffer_2000_ppdd[310];

    auto g_xx_0_0_0_z_z_yy_zz = buffer_2000_ppdd[311];

    auto g_xx_0_0_0_z_z_yz_xx = buffer_2000_ppdd[312];

    auto g_xx_0_0_0_z_z_yz_xy = buffer_2000_ppdd[313];

    auto g_xx_0_0_0_z_z_yz_xz = buffer_2000_ppdd[314];

    auto g_xx_0_0_0_z_z_yz_yy = buffer_2000_ppdd[315];

    auto g_xx_0_0_0_z_z_yz_yz = buffer_2000_ppdd[316];

    auto g_xx_0_0_0_z_z_yz_zz = buffer_2000_ppdd[317];

    auto g_xx_0_0_0_z_z_zz_xx = buffer_2000_ppdd[318];

    auto g_xx_0_0_0_z_z_zz_xy = buffer_2000_ppdd[319];

    auto g_xx_0_0_0_z_z_zz_xz = buffer_2000_ppdd[320];

    auto g_xx_0_0_0_z_z_zz_yy = buffer_2000_ppdd[321];

    auto g_xx_0_0_0_z_z_zz_yz = buffer_2000_ppdd[322];

    auto g_xx_0_0_0_z_z_zz_zz = buffer_2000_ppdd[323];

    auto g_xy_0_0_0_x_x_xx_xx = buffer_2000_ppdd[324];

    auto g_xy_0_0_0_x_x_xx_xy = buffer_2000_ppdd[325];

    auto g_xy_0_0_0_x_x_xx_xz = buffer_2000_ppdd[326];

    auto g_xy_0_0_0_x_x_xx_yy = buffer_2000_ppdd[327];

    auto g_xy_0_0_0_x_x_xx_yz = buffer_2000_ppdd[328];

    auto g_xy_0_0_0_x_x_xx_zz = buffer_2000_ppdd[329];

    auto g_xy_0_0_0_x_x_xy_xx = buffer_2000_ppdd[330];

    auto g_xy_0_0_0_x_x_xy_xy = buffer_2000_ppdd[331];

    auto g_xy_0_0_0_x_x_xy_xz = buffer_2000_ppdd[332];

    auto g_xy_0_0_0_x_x_xy_yy = buffer_2000_ppdd[333];

    auto g_xy_0_0_0_x_x_xy_yz = buffer_2000_ppdd[334];

    auto g_xy_0_0_0_x_x_xy_zz = buffer_2000_ppdd[335];

    auto g_xy_0_0_0_x_x_xz_xx = buffer_2000_ppdd[336];

    auto g_xy_0_0_0_x_x_xz_xy = buffer_2000_ppdd[337];

    auto g_xy_0_0_0_x_x_xz_xz = buffer_2000_ppdd[338];

    auto g_xy_0_0_0_x_x_xz_yy = buffer_2000_ppdd[339];

    auto g_xy_0_0_0_x_x_xz_yz = buffer_2000_ppdd[340];

    auto g_xy_0_0_0_x_x_xz_zz = buffer_2000_ppdd[341];

    auto g_xy_0_0_0_x_x_yy_xx = buffer_2000_ppdd[342];

    auto g_xy_0_0_0_x_x_yy_xy = buffer_2000_ppdd[343];

    auto g_xy_0_0_0_x_x_yy_xz = buffer_2000_ppdd[344];

    auto g_xy_0_0_0_x_x_yy_yy = buffer_2000_ppdd[345];

    auto g_xy_0_0_0_x_x_yy_yz = buffer_2000_ppdd[346];

    auto g_xy_0_0_0_x_x_yy_zz = buffer_2000_ppdd[347];

    auto g_xy_0_0_0_x_x_yz_xx = buffer_2000_ppdd[348];

    auto g_xy_0_0_0_x_x_yz_xy = buffer_2000_ppdd[349];

    auto g_xy_0_0_0_x_x_yz_xz = buffer_2000_ppdd[350];

    auto g_xy_0_0_0_x_x_yz_yy = buffer_2000_ppdd[351];

    auto g_xy_0_0_0_x_x_yz_yz = buffer_2000_ppdd[352];

    auto g_xy_0_0_0_x_x_yz_zz = buffer_2000_ppdd[353];

    auto g_xy_0_0_0_x_x_zz_xx = buffer_2000_ppdd[354];

    auto g_xy_0_0_0_x_x_zz_xy = buffer_2000_ppdd[355];

    auto g_xy_0_0_0_x_x_zz_xz = buffer_2000_ppdd[356];

    auto g_xy_0_0_0_x_x_zz_yy = buffer_2000_ppdd[357];

    auto g_xy_0_0_0_x_x_zz_yz = buffer_2000_ppdd[358];

    auto g_xy_0_0_0_x_x_zz_zz = buffer_2000_ppdd[359];

    auto g_xy_0_0_0_x_y_xx_xx = buffer_2000_ppdd[360];

    auto g_xy_0_0_0_x_y_xx_xy = buffer_2000_ppdd[361];

    auto g_xy_0_0_0_x_y_xx_xz = buffer_2000_ppdd[362];

    auto g_xy_0_0_0_x_y_xx_yy = buffer_2000_ppdd[363];

    auto g_xy_0_0_0_x_y_xx_yz = buffer_2000_ppdd[364];

    auto g_xy_0_0_0_x_y_xx_zz = buffer_2000_ppdd[365];

    auto g_xy_0_0_0_x_y_xy_xx = buffer_2000_ppdd[366];

    auto g_xy_0_0_0_x_y_xy_xy = buffer_2000_ppdd[367];

    auto g_xy_0_0_0_x_y_xy_xz = buffer_2000_ppdd[368];

    auto g_xy_0_0_0_x_y_xy_yy = buffer_2000_ppdd[369];

    auto g_xy_0_0_0_x_y_xy_yz = buffer_2000_ppdd[370];

    auto g_xy_0_0_0_x_y_xy_zz = buffer_2000_ppdd[371];

    auto g_xy_0_0_0_x_y_xz_xx = buffer_2000_ppdd[372];

    auto g_xy_0_0_0_x_y_xz_xy = buffer_2000_ppdd[373];

    auto g_xy_0_0_0_x_y_xz_xz = buffer_2000_ppdd[374];

    auto g_xy_0_0_0_x_y_xz_yy = buffer_2000_ppdd[375];

    auto g_xy_0_0_0_x_y_xz_yz = buffer_2000_ppdd[376];

    auto g_xy_0_0_0_x_y_xz_zz = buffer_2000_ppdd[377];

    auto g_xy_0_0_0_x_y_yy_xx = buffer_2000_ppdd[378];

    auto g_xy_0_0_0_x_y_yy_xy = buffer_2000_ppdd[379];

    auto g_xy_0_0_0_x_y_yy_xz = buffer_2000_ppdd[380];

    auto g_xy_0_0_0_x_y_yy_yy = buffer_2000_ppdd[381];

    auto g_xy_0_0_0_x_y_yy_yz = buffer_2000_ppdd[382];

    auto g_xy_0_0_0_x_y_yy_zz = buffer_2000_ppdd[383];

    auto g_xy_0_0_0_x_y_yz_xx = buffer_2000_ppdd[384];

    auto g_xy_0_0_0_x_y_yz_xy = buffer_2000_ppdd[385];

    auto g_xy_0_0_0_x_y_yz_xz = buffer_2000_ppdd[386];

    auto g_xy_0_0_0_x_y_yz_yy = buffer_2000_ppdd[387];

    auto g_xy_0_0_0_x_y_yz_yz = buffer_2000_ppdd[388];

    auto g_xy_0_0_0_x_y_yz_zz = buffer_2000_ppdd[389];

    auto g_xy_0_0_0_x_y_zz_xx = buffer_2000_ppdd[390];

    auto g_xy_0_0_0_x_y_zz_xy = buffer_2000_ppdd[391];

    auto g_xy_0_0_0_x_y_zz_xz = buffer_2000_ppdd[392];

    auto g_xy_0_0_0_x_y_zz_yy = buffer_2000_ppdd[393];

    auto g_xy_0_0_0_x_y_zz_yz = buffer_2000_ppdd[394];

    auto g_xy_0_0_0_x_y_zz_zz = buffer_2000_ppdd[395];

    auto g_xy_0_0_0_x_z_xx_xx = buffer_2000_ppdd[396];

    auto g_xy_0_0_0_x_z_xx_xy = buffer_2000_ppdd[397];

    auto g_xy_0_0_0_x_z_xx_xz = buffer_2000_ppdd[398];

    auto g_xy_0_0_0_x_z_xx_yy = buffer_2000_ppdd[399];

    auto g_xy_0_0_0_x_z_xx_yz = buffer_2000_ppdd[400];

    auto g_xy_0_0_0_x_z_xx_zz = buffer_2000_ppdd[401];

    auto g_xy_0_0_0_x_z_xy_xx = buffer_2000_ppdd[402];

    auto g_xy_0_0_0_x_z_xy_xy = buffer_2000_ppdd[403];

    auto g_xy_0_0_0_x_z_xy_xz = buffer_2000_ppdd[404];

    auto g_xy_0_0_0_x_z_xy_yy = buffer_2000_ppdd[405];

    auto g_xy_0_0_0_x_z_xy_yz = buffer_2000_ppdd[406];

    auto g_xy_0_0_0_x_z_xy_zz = buffer_2000_ppdd[407];

    auto g_xy_0_0_0_x_z_xz_xx = buffer_2000_ppdd[408];

    auto g_xy_0_0_0_x_z_xz_xy = buffer_2000_ppdd[409];

    auto g_xy_0_0_0_x_z_xz_xz = buffer_2000_ppdd[410];

    auto g_xy_0_0_0_x_z_xz_yy = buffer_2000_ppdd[411];

    auto g_xy_0_0_0_x_z_xz_yz = buffer_2000_ppdd[412];

    auto g_xy_0_0_0_x_z_xz_zz = buffer_2000_ppdd[413];

    auto g_xy_0_0_0_x_z_yy_xx = buffer_2000_ppdd[414];

    auto g_xy_0_0_0_x_z_yy_xy = buffer_2000_ppdd[415];

    auto g_xy_0_0_0_x_z_yy_xz = buffer_2000_ppdd[416];

    auto g_xy_0_0_0_x_z_yy_yy = buffer_2000_ppdd[417];

    auto g_xy_0_0_0_x_z_yy_yz = buffer_2000_ppdd[418];

    auto g_xy_0_0_0_x_z_yy_zz = buffer_2000_ppdd[419];

    auto g_xy_0_0_0_x_z_yz_xx = buffer_2000_ppdd[420];

    auto g_xy_0_0_0_x_z_yz_xy = buffer_2000_ppdd[421];

    auto g_xy_0_0_0_x_z_yz_xz = buffer_2000_ppdd[422];

    auto g_xy_0_0_0_x_z_yz_yy = buffer_2000_ppdd[423];

    auto g_xy_0_0_0_x_z_yz_yz = buffer_2000_ppdd[424];

    auto g_xy_0_0_0_x_z_yz_zz = buffer_2000_ppdd[425];

    auto g_xy_0_0_0_x_z_zz_xx = buffer_2000_ppdd[426];

    auto g_xy_0_0_0_x_z_zz_xy = buffer_2000_ppdd[427];

    auto g_xy_0_0_0_x_z_zz_xz = buffer_2000_ppdd[428];

    auto g_xy_0_0_0_x_z_zz_yy = buffer_2000_ppdd[429];

    auto g_xy_0_0_0_x_z_zz_yz = buffer_2000_ppdd[430];

    auto g_xy_0_0_0_x_z_zz_zz = buffer_2000_ppdd[431];

    auto g_xy_0_0_0_y_x_xx_xx = buffer_2000_ppdd[432];

    auto g_xy_0_0_0_y_x_xx_xy = buffer_2000_ppdd[433];

    auto g_xy_0_0_0_y_x_xx_xz = buffer_2000_ppdd[434];

    auto g_xy_0_0_0_y_x_xx_yy = buffer_2000_ppdd[435];

    auto g_xy_0_0_0_y_x_xx_yz = buffer_2000_ppdd[436];

    auto g_xy_0_0_0_y_x_xx_zz = buffer_2000_ppdd[437];

    auto g_xy_0_0_0_y_x_xy_xx = buffer_2000_ppdd[438];

    auto g_xy_0_0_0_y_x_xy_xy = buffer_2000_ppdd[439];

    auto g_xy_0_0_0_y_x_xy_xz = buffer_2000_ppdd[440];

    auto g_xy_0_0_0_y_x_xy_yy = buffer_2000_ppdd[441];

    auto g_xy_0_0_0_y_x_xy_yz = buffer_2000_ppdd[442];

    auto g_xy_0_0_0_y_x_xy_zz = buffer_2000_ppdd[443];

    auto g_xy_0_0_0_y_x_xz_xx = buffer_2000_ppdd[444];

    auto g_xy_0_0_0_y_x_xz_xy = buffer_2000_ppdd[445];

    auto g_xy_0_0_0_y_x_xz_xz = buffer_2000_ppdd[446];

    auto g_xy_0_0_0_y_x_xz_yy = buffer_2000_ppdd[447];

    auto g_xy_0_0_0_y_x_xz_yz = buffer_2000_ppdd[448];

    auto g_xy_0_0_0_y_x_xz_zz = buffer_2000_ppdd[449];

    auto g_xy_0_0_0_y_x_yy_xx = buffer_2000_ppdd[450];

    auto g_xy_0_0_0_y_x_yy_xy = buffer_2000_ppdd[451];

    auto g_xy_0_0_0_y_x_yy_xz = buffer_2000_ppdd[452];

    auto g_xy_0_0_0_y_x_yy_yy = buffer_2000_ppdd[453];

    auto g_xy_0_0_0_y_x_yy_yz = buffer_2000_ppdd[454];

    auto g_xy_0_0_0_y_x_yy_zz = buffer_2000_ppdd[455];

    auto g_xy_0_0_0_y_x_yz_xx = buffer_2000_ppdd[456];

    auto g_xy_0_0_0_y_x_yz_xy = buffer_2000_ppdd[457];

    auto g_xy_0_0_0_y_x_yz_xz = buffer_2000_ppdd[458];

    auto g_xy_0_0_0_y_x_yz_yy = buffer_2000_ppdd[459];

    auto g_xy_0_0_0_y_x_yz_yz = buffer_2000_ppdd[460];

    auto g_xy_0_0_0_y_x_yz_zz = buffer_2000_ppdd[461];

    auto g_xy_0_0_0_y_x_zz_xx = buffer_2000_ppdd[462];

    auto g_xy_0_0_0_y_x_zz_xy = buffer_2000_ppdd[463];

    auto g_xy_0_0_0_y_x_zz_xz = buffer_2000_ppdd[464];

    auto g_xy_0_0_0_y_x_zz_yy = buffer_2000_ppdd[465];

    auto g_xy_0_0_0_y_x_zz_yz = buffer_2000_ppdd[466];

    auto g_xy_0_0_0_y_x_zz_zz = buffer_2000_ppdd[467];

    auto g_xy_0_0_0_y_y_xx_xx = buffer_2000_ppdd[468];

    auto g_xy_0_0_0_y_y_xx_xy = buffer_2000_ppdd[469];

    auto g_xy_0_0_0_y_y_xx_xz = buffer_2000_ppdd[470];

    auto g_xy_0_0_0_y_y_xx_yy = buffer_2000_ppdd[471];

    auto g_xy_0_0_0_y_y_xx_yz = buffer_2000_ppdd[472];

    auto g_xy_0_0_0_y_y_xx_zz = buffer_2000_ppdd[473];

    auto g_xy_0_0_0_y_y_xy_xx = buffer_2000_ppdd[474];

    auto g_xy_0_0_0_y_y_xy_xy = buffer_2000_ppdd[475];

    auto g_xy_0_0_0_y_y_xy_xz = buffer_2000_ppdd[476];

    auto g_xy_0_0_0_y_y_xy_yy = buffer_2000_ppdd[477];

    auto g_xy_0_0_0_y_y_xy_yz = buffer_2000_ppdd[478];

    auto g_xy_0_0_0_y_y_xy_zz = buffer_2000_ppdd[479];

    auto g_xy_0_0_0_y_y_xz_xx = buffer_2000_ppdd[480];

    auto g_xy_0_0_0_y_y_xz_xy = buffer_2000_ppdd[481];

    auto g_xy_0_0_0_y_y_xz_xz = buffer_2000_ppdd[482];

    auto g_xy_0_0_0_y_y_xz_yy = buffer_2000_ppdd[483];

    auto g_xy_0_0_0_y_y_xz_yz = buffer_2000_ppdd[484];

    auto g_xy_0_0_0_y_y_xz_zz = buffer_2000_ppdd[485];

    auto g_xy_0_0_0_y_y_yy_xx = buffer_2000_ppdd[486];

    auto g_xy_0_0_0_y_y_yy_xy = buffer_2000_ppdd[487];

    auto g_xy_0_0_0_y_y_yy_xz = buffer_2000_ppdd[488];

    auto g_xy_0_0_0_y_y_yy_yy = buffer_2000_ppdd[489];

    auto g_xy_0_0_0_y_y_yy_yz = buffer_2000_ppdd[490];

    auto g_xy_0_0_0_y_y_yy_zz = buffer_2000_ppdd[491];

    auto g_xy_0_0_0_y_y_yz_xx = buffer_2000_ppdd[492];

    auto g_xy_0_0_0_y_y_yz_xy = buffer_2000_ppdd[493];

    auto g_xy_0_0_0_y_y_yz_xz = buffer_2000_ppdd[494];

    auto g_xy_0_0_0_y_y_yz_yy = buffer_2000_ppdd[495];

    auto g_xy_0_0_0_y_y_yz_yz = buffer_2000_ppdd[496];

    auto g_xy_0_0_0_y_y_yz_zz = buffer_2000_ppdd[497];

    auto g_xy_0_0_0_y_y_zz_xx = buffer_2000_ppdd[498];

    auto g_xy_0_0_0_y_y_zz_xy = buffer_2000_ppdd[499];

    auto g_xy_0_0_0_y_y_zz_xz = buffer_2000_ppdd[500];

    auto g_xy_0_0_0_y_y_zz_yy = buffer_2000_ppdd[501];

    auto g_xy_0_0_0_y_y_zz_yz = buffer_2000_ppdd[502];

    auto g_xy_0_0_0_y_y_zz_zz = buffer_2000_ppdd[503];

    auto g_xy_0_0_0_y_z_xx_xx = buffer_2000_ppdd[504];

    auto g_xy_0_0_0_y_z_xx_xy = buffer_2000_ppdd[505];

    auto g_xy_0_0_0_y_z_xx_xz = buffer_2000_ppdd[506];

    auto g_xy_0_0_0_y_z_xx_yy = buffer_2000_ppdd[507];

    auto g_xy_0_0_0_y_z_xx_yz = buffer_2000_ppdd[508];

    auto g_xy_0_0_0_y_z_xx_zz = buffer_2000_ppdd[509];

    auto g_xy_0_0_0_y_z_xy_xx = buffer_2000_ppdd[510];

    auto g_xy_0_0_0_y_z_xy_xy = buffer_2000_ppdd[511];

    auto g_xy_0_0_0_y_z_xy_xz = buffer_2000_ppdd[512];

    auto g_xy_0_0_0_y_z_xy_yy = buffer_2000_ppdd[513];

    auto g_xy_0_0_0_y_z_xy_yz = buffer_2000_ppdd[514];

    auto g_xy_0_0_0_y_z_xy_zz = buffer_2000_ppdd[515];

    auto g_xy_0_0_0_y_z_xz_xx = buffer_2000_ppdd[516];

    auto g_xy_0_0_0_y_z_xz_xy = buffer_2000_ppdd[517];

    auto g_xy_0_0_0_y_z_xz_xz = buffer_2000_ppdd[518];

    auto g_xy_0_0_0_y_z_xz_yy = buffer_2000_ppdd[519];

    auto g_xy_0_0_0_y_z_xz_yz = buffer_2000_ppdd[520];

    auto g_xy_0_0_0_y_z_xz_zz = buffer_2000_ppdd[521];

    auto g_xy_0_0_0_y_z_yy_xx = buffer_2000_ppdd[522];

    auto g_xy_0_0_0_y_z_yy_xy = buffer_2000_ppdd[523];

    auto g_xy_0_0_0_y_z_yy_xz = buffer_2000_ppdd[524];

    auto g_xy_0_0_0_y_z_yy_yy = buffer_2000_ppdd[525];

    auto g_xy_0_0_0_y_z_yy_yz = buffer_2000_ppdd[526];

    auto g_xy_0_0_0_y_z_yy_zz = buffer_2000_ppdd[527];

    auto g_xy_0_0_0_y_z_yz_xx = buffer_2000_ppdd[528];

    auto g_xy_0_0_0_y_z_yz_xy = buffer_2000_ppdd[529];

    auto g_xy_0_0_0_y_z_yz_xz = buffer_2000_ppdd[530];

    auto g_xy_0_0_0_y_z_yz_yy = buffer_2000_ppdd[531];

    auto g_xy_0_0_0_y_z_yz_yz = buffer_2000_ppdd[532];

    auto g_xy_0_0_0_y_z_yz_zz = buffer_2000_ppdd[533];

    auto g_xy_0_0_0_y_z_zz_xx = buffer_2000_ppdd[534];

    auto g_xy_0_0_0_y_z_zz_xy = buffer_2000_ppdd[535];

    auto g_xy_0_0_0_y_z_zz_xz = buffer_2000_ppdd[536];

    auto g_xy_0_0_0_y_z_zz_yy = buffer_2000_ppdd[537];

    auto g_xy_0_0_0_y_z_zz_yz = buffer_2000_ppdd[538];

    auto g_xy_0_0_0_y_z_zz_zz = buffer_2000_ppdd[539];

    auto g_xy_0_0_0_z_x_xx_xx = buffer_2000_ppdd[540];

    auto g_xy_0_0_0_z_x_xx_xy = buffer_2000_ppdd[541];

    auto g_xy_0_0_0_z_x_xx_xz = buffer_2000_ppdd[542];

    auto g_xy_0_0_0_z_x_xx_yy = buffer_2000_ppdd[543];

    auto g_xy_0_0_0_z_x_xx_yz = buffer_2000_ppdd[544];

    auto g_xy_0_0_0_z_x_xx_zz = buffer_2000_ppdd[545];

    auto g_xy_0_0_0_z_x_xy_xx = buffer_2000_ppdd[546];

    auto g_xy_0_0_0_z_x_xy_xy = buffer_2000_ppdd[547];

    auto g_xy_0_0_0_z_x_xy_xz = buffer_2000_ppdd[548];

    auto g_xy_0_0_0_z_x_xy_yy = buffer_2000_ppdd[549];

    auto g_xy_0_0_0_z_x_xy_yz = buffer_2000_ppdd[550];

    auto g_xy_0_0_0_z_x_xy_zz = buffer_2000_ppdd[551];

    auto g_xy_0_0_0_z_x_xz_xx = buffer_2000_ppdd[552];

    auto g_xy_0_0_0_z_x_xz_xy = buffer_2000_ppdd[553];

    auto g_xy_0_0_0_z_x_xz_xz = buffer_2000_ppdd[554];

    auto g_xy_0_0_0_z_x_xz_yy = buffer_2000_ppdd[555];

    auto g_xy_0_0_0_z_x_xz_yz = buffer_2000_ppdd[556];

    auto g_xy_0_0_0_z_x_xz_zz = buffer_2000_ppdd[557];

    auto g_xy_0_0_0_z_x_yy_xx = buffer_2000_ppdd[558];

    auto g_xy_0_0_0_z_x_yy_xy = buffer_2000_ppdd[559];

    auto g_xy_0_0_0_z_x_yy_xz = buffer_2000_ppdd[560];

    auto g_xy_0_0_0_z_x_yy_yy = buffer_2000_ppdd[561];

    auto g_xy_0_0_0_z_x_yy_yz = buffer_2000_ppdd[562];

    auto g_xy_0_0_0_z_x_yy_zz = buffer_2000_ppdd[563];

    auto g_xy_0_0_0_z_x_yz_xx = buffer_2000_ppdd[564];

    auto g_xy_0_0_0_z_x_yz_xy = buffer_2000_ppdd[565];

    auto g_xy_0_0_0_z_x_yz_xz = buffer_2000_ppdd[566];

    auto g_xy_0_0_0_z_x_yz_yy = buffer_2000_ppdd[567];

    auto g_xy_0_0_0_z_x_yz_yz = buffer_2000_ppdd[568];

    auto g_xy_0_0_0_z_x_yz_zz = buffer_2000_ppdd[569];

    auto g_xy_0_0_0_z_x_zz_xx = buffer_2000_ppdd[570];

    auto g_xy_0_0_0_z_x_zz_xy = buffer_2000_ppdd[571];

    auto g_xy_0_0_0_z_x_zz_xz = buffer_2000_ppdd[572];

    auto g_xy_0_0_0_z_x_zz_yy = buffer_2000_ppdd[573];

    auto g_xy_0_0_0_z_x_zz_yz = buffer_2000_ppdd[574];

    auto g_xy_0_0_0_z_x_zz_zz = buffer_2000_ppdd[575];

    auto g_xy_0_0_0_z_y_xx_xx = buffer_2000_ppdd[576];

    auto g_xy_0_0_0_z_y_xx_xy = buffer_2000_ppdd[577];

    auto g_xy_0_0_0_z_y_xx_xz = buffer_2000_ppdd[578];

    auto g_xy_0_0_0_z_y_xx_yy = buffer_2000_ppdd[579];

    auto g_xy_0_0_0_z_y_xx_yz = buffer_2000_ppdd[580];

    auto g_xy_0_0_0_z_y_xx_zz = buffer_2000_ppdd[581];

    auto g_xy_0_0_0_z_y_xy_xx = buffer_2000_ppdd[582];

    auto g_xy_0_0_0_z_y_xy_xy = buffer_2000_ppdd[583];

    auto g_xy_0_0_0_z_y_xy_xz = buffer_2000_ppdd[584];

    auto g_xy_0_0_0_z_y_xy_yy = buffer_2000_ppdd[585];

    auto g_xy_0_0_0_z_y_xy_yz = buffer_2000_ppdd[586];

    auto g_xy_0_0_0_z_y_xy_zz = buffer_2000_ppdd[587];

    auto g_xy_0_0_0_z_y_xz_xx = buffer_2000_ppdd[588];

    auto g_xy_0_0_0_z_y_xz_xy = buffer_2000_ppdd[589];

    auto g_xy_0_0_0_z_y_xz_xz = buffer_2000_ppdd[590];

    auto g_xy_0_0_0_z_y_xz_yy = buffer_2000_ppdd[591];

    auto g_xy_0_0_0_z_y_xz_yz = buffer_2000_ppdd[592];

    auto g_xy_0_0_0_z_y_xz_zz = buffer_2000_ppdd[593];

    auto g_xy_0_0_0_z_y_yy_xx = buffer_2000_ppdd[594];

    auto g_xy_0_0_0_z_y_yy_xy = buffer_2000_ppdd[595];

    auto g_xy_0_0_0_z_y_yy_xz = buffer_2000_ppdd[596];

    auto g_xy_0_0_0_z_y_yy_yy = buffer_2000_ppdd[597];

    auto g_xy_0_0_0_z_y_yy_yz = buffer_2000_ppdd[598];

    auto g_xy_0_0_0_z_y_yy_zz = buffer_2000_ppdd[599];

    auto g_xy_0_0_0_z_y_yz_xx = buffer_2000_ppdd[600];

    auto g_xy_0_0_0_z_y_yz_xy = buffer_2000_ppdd[601];

    auto g_xy_0_0_0_z_y_yz_xz = buffer_2000_ppdd[602];

    auto g_xy_0_0_0_z_y_yz_yy = buffer_2000_ppdd[603];

    auto g_xy_0_0_0_z_y_yz_yz = buffer_2000_ppdd[604];

    auto g_xy_0_0_0_z_y_yz_zz = buffer_2000_ppdd[605];

    auto g_xy_0_0_0_z_y_zz_xx = buffer_2000_ppdd[606];

    auto g_xy_0_0_0_z_y_zz_xy = buffer_2000_ppdd[607];

    auto g_xy_0_0_0_z_y_zz_xz = buffer_2000_ppdd[608];

    auto g_xy_0_0_0_z_y_zz_yy = buffer_2000_ppdd[609];

    auto g_xy_0_0_0_z_y_zz_yz = buffer_2000_ppdd[610];

    auto g_xy_0_0_0_z_y_zz_zz = buffer_2000_ppdd[611];

    auto g_xy_0_0_0_z_z_xx_xx = buffer_2000_ppdd[612];

    auto g_xy_0_0_0_z_z_xx_xy = buffer_2000_ppdd[613];

    auto g_xy_0_0_0_z_z_xx_xz = buffer_2000_ppdd[614];

    auto g_xy_0_0_0_z_z_xx_yy = buffer_2000_ppdd[615];

    auto g_xy_0_0_0_z_z_xx_yz = buffer_2000_ppdd[616];

    auto g_xy_0_0_0_z_z_xx_zz = buffer_2000_ppdd[617];

    auto g_xy_0_0_0_z_z_xy_xx = buffer_2000_ppdd[618];

    auto g_xy_0_0_0_z_z_xy_xy = buffer_2000_ppdd[619];

    auto g_xy_0_0_0_z_z_xy_xz = buffer_2000_ppdd[620];

    auto g_xy_0_0_0_z_z_xy_yy = buffer_2000_ppdd[621];

    auto g_xy_0_0_0_z_z_xy_yz = buffer_2000_ppdd[622];

    auto g_xy_0_0_0_z_z_xy_zz = buffer_2000_ppdd[623];

    auto g_xy_0_0_0_z_z_xz_xx = buffer_2000_ppdd[624];

    auto g_xy_0_0_0_z_z_xz_xy = buffer_2000_ppdd[625];

    auto g_xy_0_0_0_z_z_xz_xz = buffer_2000_ppdd[626];

    auto g_xy_0_0_0_z_z_xz_yy = buffer_2000_ppdd[627];

    auto g_xy_0_0_0_z_z_xz_yz = buffer_2000_ppdd[628];

    auto g_xy_0_0_0_z_z_xz_zz = buffer_2000_ppdd[629];

    auto g_xy_0_0_0_z_z_yy_xx = buffer_2000_ppdd[630];

    auto g_xy_0_0_0_z_z_yy_xy = buffer_2000_ppdd[631];

    auto g_xy_0_0_0_z_z_yy_xz = buffer_2000_ppdd[632];

    auto g_xy_0_0_0_z_z_yy_yy = buffer_2000_ppdd[633];

    auto g_xy_0_0_0_z_z_yy_yz = buffer_2000_ppdd[634];

    auto g_xy_0_0_0_z_z_yy_zz = buffer_2000_ppdd[635];

    auto g_xy_0_0_0_z_z_yz_xx = buffer_2000_ppdd[636];

    auto g_xy_0_0_0_z_z_yz_xy = buffer_2000_ppdd[637];

    auto g_xy_0_0_0_z_z_yz_xz = buffer_2000_ppdd[638];

    auto g_xy_0_0_0_z_z_yz_yy = buffer_2000_ppdd[639];

    auto g_xy_0_0_0_z_z_yz_yz = buffer_2000_ppdd[640];

    auto g_xy_0_0_0_z_z_yz_zz = buffer_2000_ppdd[641];

    auto g_xy_0_0_0_z_z_zz_xx = buffer_2000_ppdd[642];

    auto g_xy_0_0_0_z_z_zz_xy = buffer_2000_ppdd[643];

    auto g_xy_0_0_0_z_z_zz_xz = buffer_2000_ppdd[644];

    auto g_xy_0_0_0_z_z_zz_yy = buffer_2000_ppdd[645];

    auto g_xy_0_0_0_z_z_zz_yz = buffer_2000_ppdd[646];

    auto g_xy_0_0_0_z_z_zz_zz = buffer_2000_ppdd[647];

    auto g_xz_0_0_0_x_x_xx_xx = buffer_2000_ppdd[648];

    auto g_xz_0_0_0_x_x_xx_xy = buffer_2000_ppdd[649];

    auto g_xz_0_0_0_x_x_xx_xz = buffer_2000_ppdd[650];

    auto g_xz_0_0_0_x_x_xx_yy = buffer_2000_ppdd[651];

    auto g_xz_0_0_0_x_x_xx_yz = buffer_2000_ppdd[652];

    auto g_xz_0_0_0_x_x_xx_zz = buffer_2000_ppdd[653];

    auto g_xz_0_0_0_x_x_xy_xx = buffer_2000_ppdd[654];

    auto g_xz_0_0_0_x_x_xy_xy = buffer_2000_ppdd[655];

    auto g_xz_0_0_0_x_x_xy_xz = buffer_2000_ppdd[656];

    auto g_xz_0_0_0_x_x_xy_yy = buffer_2000_ppdd[657];

    auto g_xz_0_0_0_x_x_xy_yz = buffer_2000_ppdd[658];

    auto g_xz_0_0_0_x_x_xy_zz = buffer_2000_ppdd[659];

    auto g_xz_0_0_0_x_x_xz_xx = buffer_2000_ppdd[660];

    auto g_xz_0_0_0_x_x_xz_xy = buffer_2000_ppdd[661];

    auto g_xz_0_0_0_x_x_xz_xz = buffer_2000_ppdd[662];

    auto g_xz_0_0_0_x_x_xz_yy = buffer_2000_ppdd[663];

    auto g_xz_0_0_0_x_x_xz_yz = buffer_2000_ppdd[664];

    auto g_xz_0_0_0_x_x_xz_zz = buffer_2000_ppdd[665];

    auto g_xz_0_0_0_x_x_yy_xx = buffer_2000_ppdd[666];

    auto g_xz_0_0_0_x_x_yy_xy = buffer_2000_ppdd[667];

    auto g_xz_0_0_0_x_x_yy_xz = buffer_2000_ppdd[668];

    auto g_xz_0_0_0_x_x_yy_yy = buffer_2000_ppdd[669];

    auto g_xz_0_0_0_x_x_yy_yz = buffer_2000_ppdd[670];

    auto g_xz_0_0_0_x_x_yy_zz = buffer_2000_ppdd[671];

    auto g_xz_0_0_0_x_x_yz_xx = buffer_2000_ppdd[672];

    auto g_xz_0_0_0_x_x_yz_xy = buffer_2000_ppdd[673];

    auto g_xz_0_0_0_x_x_yz_xz = buffer_2000_ppdd[674];

    auto g_xz_0_0_0_x_x_yz_yy = buffer_2000_ppdd[675];

    auto g_xz_0_0_0_x_x_yz_yz = buffer_2000_ppdd[676];

    auto g_xz_0_0_0_x_x_yz_zz = buffer_2000_ppdd[677];

    auto g_xz_0_0_0_x_x_zz_xx = buffer_2000_ppdd[678];

    auto g_xz_0_0_0_x_x_zz_xy = buffer_2000_ppdd[679];

    auto g_xz_0_0_0_x_x_zz_xz = buffer_2000_ppdd[680];

    auto g_xz_0_0_0_x_x_zz_yy = buffer_2000_ppdd[681];

    auto g_xz_0_0_0_x_x_zz_yz = buffer_2000_ppdd[682];

    auto g_xz_0_0_0_x_x_zz_zz = buffer_2000_ppdd[683];

    auto g_xz_0_0_0_x_y_xx_xx = buffer_2000_ppdd[684];

    auto g_xz_0_0_0_x_y_xx_xy = buffer_2000_ppdd[685];

    auto g_xz_0_0_0_x_y_xx_xz = buffer_2000_ppdd[686];

    auto g_xz_0_0_0_x_y_xx_yy = buffer_2000_ppdd[687];

    auto g_xz_0_0_0_x_y_xx_yz = buffer_2000_ppdd[688];

    auto g_xz_0_0_0_x_y_xx_zz = buffer_2000_ppdd[689];

    auto g_xz_0_0_0_x_y_xy_xx = buffer_2000_ppdd[690];

    auto g_xz_0_0_0_x_y_xy_xy = buffer_2000_ppdd[691];

    auto g_xz_0_0_0_x_y_xy_xz = buffer_2000_ppdd[692];

    auto g_xz_0_0_0_x_y_xy_yy = buffer_2000_ppdd[693];

    auto g_xz_0_0_0_x_y_xy_yz = buffer_2000_ppdd[694];

    auto g_xz_0_0_0_x_y_xy_zz = buffer_2000_ppdd[695];

    auto g_xz_0_0_0_x_y_xz_xx = buffer_2000_ppdd[696];

    auto g_xz_0_0_0_x_y_xz_xy = buffer_2000_ppdd[697];

    auto g_xz_0_0_0_x_y_xz_xz = buffer_2000_ppdd[698];

    auto g_xz_0_0_0_x_y_xz_yy = buffer_2000_ppdd[699];

    auto g_xz_0_0_0_x_y_xz_yz = buffer_2000_ppdd[700];

    auto g_xz_0_0_0_x_y_xz_zz = buffer_2000_ppdd[701];

    auto g_xz_0_0_0_x_y_yy_xx = buffer_2000_ppdd[702];

    auto g_xz_0_0_0_x_y_yy_xy = buffer_2000_ppdd[703];

    auto g_xz_0_0_0_x_y_yy_xz = buffer_2000_ppdd[704];

    auto g_xz_0_0_0_x_y_yy_yy = buffer_2000_ppdd[705];

    auto g_xz_0_0_0_x_y_yy_yz = buffer_2000_ppdd[706];

    auto g_xz_0_0_0_x_y_yy_zz = buffer_2000_ppdd[707];

    auto g_xz_0_0_0_x_y_yz_xx = buffer_2000_ppdd[708];

    auto g_xz_0_0_0_x_y_yz_xy = buffer_2000_ppdd[709];

    auto g_xz_0_0_0_x_y_yz_xz = buffer_2000_ppdd[710];

    auto g_xz_0_0_0_x_y_yz_yy = buffer_2000_ppdd[711];

    auto g_xz_0_0_0_x_y_yz_yz = buffer_2000_ppdd[712];

    auto g_xz_0_0_0_x_y_yz_zz = buffer_2000_ppdd[713];

    auto g_xz_0_0_0_x_y_zz_xx = buffer_2000_ppdd[714];

    auto g_xz_0_0_0_x_y_zz_xy = buffer_2000_ppdd[715];

    auto g_xz_0_0_0_x_y_zz_xz = buffer_2000_ppdd[716];

    auto g_xz_0_0_0_x_y_zz_yy = buffer_2000_ppdd[717];

    auto g_xz_0_0_0_x_y_zz_yz = buffer_2000_ppdd[718];

    auto g_xz_0_0_0_x_y_zz_zz = buffer_2000_ppdd[719];

    auto g_xz_0_0_0_x_z_xx_xx = buffer_2000_ppdd[720];

    auto g_xz_0_0_0_x_z_xx_xy = buffer_2000_ppdd[721];

    auto g_xz_0_0_0_x_z_xx_xz = buffer_2000_ppdd[722];

    auto g_xz_0_0_0_x_z_xx_yy = buffer_2000_ppdd[723];

    auto g_xz_0_0_0_x_z_xx_yz = buffer_2000_ppdd[724];

    auto g_xz_0_0_0_x_z_xx_zz = buffer_2000_ppdd[725];

    auto g_xz_0_0_0_x_z_xy_xx = buffer_2000_ppdd[726];

    auto g_xz_0_0_0_x_z_xy_xy = buffer_2000_ppdd[727];

    auto g_xz_0_0_0_x_z_xy_xz = buffer_2000_ppdd[728];

    auto g_xz_0_0_0_x_z_xy_yy = buffer_2000_ppdd[729];

    auto g_xz_0_0_0_x_z_xy_yz = buffer_2000_ppdd[730];

    auto g_xz_0_0_0_x_z_xy_zz = buffer_2000_ppdd[731];

    auto g_xz_0_0_0_x_z_xz_xx = buffer_2000_ppdd[732];

    auto g_xz_0_0_0_x_z_xz_xy = buffer_2000_ppdd[733];

    auto g_xz_0_0_0_x_z_xz_xz = buffer_2000_ppdd[734];

    auto g_xz_0_0_0_x_z_xz_yy = buffer_2000_ppdd[735];

    auto g_xz_0_0_0_x_z_xz_yz = buffer_2000_ppdd[736];

    auto g_xz_0_0_0_x_z_xz_zz = buffer_2000_ppdd[737];

    auto g_xz_0_0_0_x_z_yy_xx = buffer_2000_ppdd[738];

    auto g_xz_0_0_0_x_z_yy_xy = buffer_2000_ppdd[739];

    auto g_xz_0_0_0_x_z_yy_xz = buffer_2000_ppdd[740];

    auto g_xz_0_0_0_x_z_yy_yy = buffer_2000_ppdd[741];

    auto g_xz_0_0_0_x_z_yy_yz = buffer_2000_ppdd[742];

    auto g_xz_0_0_0_x_z_yy_zz = buffer_2000_ppdd[743];

    auto g_xz_0_0_0_x_z_yz_xx = buffer_2000_ppdd[744];

    auto g_xz_0_0_0_x_z_yz_xy = buffer_2000_ppdd[745];

    auto g_xz_0_0_0_x_z_yz_xz = buffer_2000_ppdd[746];

    auto g_xz_0_0_0_x_z_yz_yy = buffer_2000_ppdd[747];

    auto g_xz_0_0_0_x_z_yz_yz = buffer_2000_ppdd[748];

    auto g_xz_0_0_0_x_z_yz_zz = buffer_2000_ppdd[749];

    auto g_xz_0_0_0_x_z_zz_xx = buffer_2000_ppdd[750];

    auto g_xz_0_0_0_x_z_zz_xy = buffer_2000_ppdd[751];

    auto g_xz_0_0_0_x_z_zz_xz = buffer_2000_ppdd[752];

    auto g_xz_0_0_0_x_z_zz_yy = buffer_2000_ppdd[753];

    auto g_xz_0_0_0_x_z_zz_yz = buffer_2000_ppdd[754];

    auto g_xz_0_0_0_x_z_zz_zz = buffer_2000_ppdd[755];

    auto g_xz_0_0_0_y_x_xx_xx = buffer_2000_ppdd[756];

    auto g_xz_0_0_0_y_x_xx_xy = buffer_2000_ppdd[757];

    auto g_xz_0_0_0_y_x_xx_xz = buffer_2000_ppdd[758];

    auto g_xz_0_0_0_y_x_xx_yy = buffer_2000_ppdd[759];

    auto g_xz_0_0_0_y_x_xx_yz = buffer_2000_ppdd[760];

    auto g_xz_0_0_0_y_x_xx_zz = buffer_2000_ppdd[761];

    auto g_xz_0_0_0_y_x_xy_xx = buffer_2000_ppdd[762];

    auto g_xz_0_0_0_y_x_xy_xy = buffer_2000_ppdd[763];

    auto g_xz_0_0_0_y_x_xy_xz = buffer_2000_ppdd[764];

    auto g_xz_0_0_0_y_x_xy_yy = buffer_2000_ppdd[765];

    auto g_xz_0_0_0_y_x_xy_yz = buffer_2000_ppdd[766];

    auto g_xz_0_0_0_y_x_xy_zz = buffer_2000_ppdd[767];

    auto g_xz_0_0_0_y_x_xz_xx = buffer_2000_ppdd[768];

    auto g_xz_0_0_0_y_x_xz_xy = buffer_2000_ppdd[769];

    auto g_xz_0_0_0_y_x_xz_xz = buffer_2000_ppdd[770];

    auto g_xz_0_0_0_y_x_xz_yy = buffer_2000_ppdd[771];

    auto g_xz_0_0_0_y_x_xz_yz = buffer_2000_ppdd[772];

    auto g_xz_0_0_0_y_x_xz_zz = buffer_2000_ppdd[773];

    auto g_xz_0_0_0_y_x_yy_xx = buffer_2000_ppdd[774];

    auto g_xz_0_0_0_y_x_yy_xy = buffer_2000_ppdd[775];

    auto g_xz_0_0_0_y_x_yy_xz = buffer_2000_ppdd[776];

    auto g_xz_0_0_0_y_x_yy_yy = buffer_2000_ppdd[777];

    auto g_xz_0_0_0_y_x_yy_yz = buffer_2000_ppdd[778];

    auto g_xz_0_0_0_y_x_yy_zz = buffer_2000_ppdd[779];

    auto g_xz_0_0_0_y_x_yz_xx = buffer_2000_ppdd[780];

    auto g_xz_0_0_0_y_x_yz_xy = buffer_2000_ppdd[781];

    auto g_xz_0_0_0_y_x_yz_xz = buffer_2000_ppdd[782];

    auto g_xz_0_0_0_y_x_yz_yy = buffer_2000_ppdd[783];

    auto g_xz_0_0_0_y_x_yz_yz = buffer_2000_ppdd[784];

    auto g_xz_0_0_0_y_x_yz_zz = buffer_2000_ppdd[785];

    auto g_xz_0_0_0_y_x_zz_xx = buffer_2000_ppdd[786];

    auto g_xz_0_0_0_y_x_zz_xy = buffer_2000_ppdd[787];

    auto g_xz_0_0_0_y_x_zz_xz = buffer_2000_ppdd[788];

    auto g_xz_0_0_0_y_x_zz_yy = buffer_2000_ppdd[789];

    auto g_xz_0_0_0_y_x_zz_yz = buffer_2000_ppdd[790];

    auto g_xz_0_0_0_y_x_zz_zz = buffer_2000_ppdd[791];

    auto g_xz_0_0_0_y_y_xx_xx = buffer_2000_ppdd[792];

    auto g_xz_0_0_0_y_y_xx_xy = buffer_2000_ppdd[793];

    auto g_xz_0_0_0_y_y_xx_xz = buffer_2000_ppdd[794];

    auto g_xz_0_0_0_y_y_xx_yy = buffer_2000_ppdd[795];

    auto g_xz_0_0_0_y_y_xx_yz = buffer_2000_ppdd[796];

    auto g_xz_0_0_0_y_y_xx_zz = buffer_2000_ppdd[797];

    auto g_xz_0_0_0_y_y_xy_xx = buffer_2000_ppdd[798];

    auto g_xz_0_0_0_y_y_xy_xy = buffer_2000_ppdd[799];

    auto g_xz_0_0_0_y_y_xy_xz = buffer_2000_ppdd[800];

    auto g_xz_0_0_0_y_y_xy_yy = buffer_2000_ppdd[801];

    auto g_xz_0_0_0_y_y_xy_yz = buffer_2000_ppdd[802];

    auto g_xz_0_0_0_y_y_xy_zz = buffer_2000_ppdd[803];

    auto g_xz_0_0_0_y_y_xz_xx = buffer_2000_ppdd[804];

    auto g_xz_0_0_0_y_y_xz_xy = buffer_2000_ppdd[805];

    auto g_xz_0_0_0_y_y_xz_xz = buffer_2000_ppdd[806];

    auto g_xz_0_0_0_y_y_xz_yy = buffer_2000_ppdd[807];

    auto g_xz_0_0_0_y_y_xz_yz = buffer_2000_ppdd[808];

    auto g_xz_0_0_0_y_y_xz_zz = buffer_2000_ppdd[809];

    auto g_xz_0_0_0_y_y_yy_xx = buffer_2000_ppdd[810];

    auto g_xz_0_0_0_y_y_yy_xy = buffer_2000_ppdd[811];

    auto g_xz_0_0_0_y_y_yy_xz = buffer_2000_ppdd[812];

    auto g_xz_0_0_0_y_y_yy_yy = buffer_2000_ppdd[813];

    auto g_xz_0_0_0_y_y_yy_yz = buffer_2000_ppdd[814];

    auto g_xz_0_0_0_y_y_yy_zz = buffer_2000_ppdd[815];

    auto g_xz_0_0_0_y_y_yz_xx = buffer_2000_ppdd[816];

    auto g_xz_0_0_0_y_y_yz_xy = buffer_2000_ppdd[817];

    auto g_xz_0_0_0_y_y_yz_xz = buffer_2000_ppdd[818];

    auto g_xz_0_0_0_y_y_yz_yy = buffer_2000_ppdd[819];

    auto g_xz_0_0_0_y_y_yz_yz = buffer_2000_ppdd[820];

    auto g_xz_0_0_0_y_y_yz_zz = buffer_2000_ppdd[821];

    auto g_xz_0_0_0_y_y_zz_xx = buffer_2000_ppdd[822];

    auto g_xz_0_0_0_y_y_zz_xy = buffer_2000_ppdd[823];

    auto g_xz_0_0_0_y_y_zz_xz = buffer_2000_ppdd[824];

    auto g_xz_0_0_0_y_y_zz_yy = buffer_2000_ppdd[825];

    auto g_xz_0_0_0_y_y_zz_yz = buffer_2000_ppdd[826];

    auto g_xz_0_0_0_y_y_zz_zz = buffer_2000_ppdd[827];

    auto g_xz_0_0_0_y_z_xx_xx = buffer_2000_ppdd[828];

    auto g_xz_0_0_0_y_z_xx_xy = buffer_2000_ppdd[829];

    auto g_xz_0_0_0_y_z_xx_xz = buffer_2000_ppdd[830];

    auto g_xz_0_0_0_y_z_xx_yy = buffer_2000_ppdd[831];

    auto g_xz_0_0_0_y_z_xx_yz = buffer_2000_ppdd[832];

    auto g_xz_0_0_0_y_z_xx_zz = buffer_2000_ppdd[833];

    auto g_xz_0_0_0_y_z_xy_xx = buffer_2000_ppdd[834];

    auto g_xz_0_0_0_y_z_xy_xy = buffer_2000_ppdd[835];

    auto g_xz_0_0_0_y_z_xy_xz = buffer_2000_ppdd[836];

    auto g_xz_0_0_0_y_z_xy_yy = buffer_2000_ppdd[837];

    auto g_xz_0_0_0_y_z_xy_yz = buffer_2000_ppdd[838];

    auto g_xz_0_0_0_y_z_xy_zz = buffer_2000_ppdd[839];

    auto g_xz_0_0_0_y_z_xz_xx = buffer_2000_ppdd[840];

    auto g_xz_0_0_0_y_z_xz_xy = buffer_2000_ppdd[841];

    auto g_xz_0_0_0_y_z_xz_xz = buffer_2000_ppdd[842];

    auto g_xz_0_0_0_y_z_xz_yy = buffer_2000_ppdd[843];

    auto g_xz_0_0_0_y_z_xz_yz = buffer_2000_ppdd[844];

    auto g_xz_0_0_0_y_z_xz_zz = buffer_2000_ppdd[845];

    auto g_xz_0_0_0_y_z_yy_xx = buffer_2000_ppdd[846];

    auto g_xz_0_0_0_y_z_yy_xy = buffer_2000_ppdd[847];

    auto g_xz_0_0_0_y_z_yy_xz = buffer_2000_ppdd[848];

    auto g_xz_0_0_0_y_z_yy_yy = buffer_2000_ppdd[849];

    auto g_xz_0_0_0_y_z_yy_yz = buffer_2000_ppdd[850];

    auto g_xz_0_0_0_y_z_yy_zz = buffer_2000_ppdd[851];

    auto g_xz_0_0_0_y_z_yz_xx = buffer_2000_ppdd[852];

    auto g_xz_0_0_0_y_z_yz_xy = buffer_2000_ppdd[853];

    auto g_xz_0_0_0_y_z_yz_xz = buffer_2000_ppdd[854];

    auto g_xz_0_0_0_y_z_yz_yy = buffer_2000_ppdd[855];

    auto g_xz_0_0_0_y_z_yz_yz = buffer_2000_ppdd[856];

    auto g_xz_0_0_0_y_z_yz_zz = buffer_2000_ppdd[857];

    auto g_xz_0_0_0_y_z_zz_xx = buffer_2000_ppdd[858];

    auto g_xz_0_0_0_y_z_zz_xy = buffer_2000_ppdd[859];

    auto g_xz_0_0_0_y_z_zz_xz = buffer_2000_ppdd[860];

    auto g_xz_0_0_0_y_z_zz_yy = buffer_2000_ppdd[861];

    auto g_xz_0_0_0_y_z_zz_yz = buffer_2000_ppdd[862];

    auto g_xz_0_0_0_y_z_zz_zz = buffer_2000_ppdd[863];

    auto g_xz_0_0_0_z_x_xx_xx = buffer_2000_ppdd[864];

    auto g_xz_0_0_0_z_x_xx_xy = buffer_2000_ppdd[865];

    auto g_xz_0_0_0_z_x_xx_xz = buffer_2000_ppdd[866];

    auto g_xz_0_0_0_z_x_xx_yy = buffer_2000_ppdd[867];

    auto g_xz_0_0_0_z_x_xx_yz = buffer_2000_ppdd[868];

    auto g_xz_0_0_0_z_x_xx_zz = buffer_2000_ppdd[869];

    auto g_xz_0_0_0_z_x_xy_xx = buffer_2000_ppdd[870];

    auto g_xz_0_0_0_z_x_xy_xy = buffer_2000_ppdd[871];

    auto g_xz_0_0_0_z_x_xy_xz = buffer_2000_ppdd[872];

    auto g_xz_0_0_0_z_x_xy_yy = buffer_2000_ppdd[873];

    auto g_xz_0_0_0_z_x_xy_yz = buffer_2000_ppdd[874];

    auto g_xz_0_0_0_z_x_xy_zz = buffer_2000_ppdd[875];

    auto g_xz_0_0_0_z_x_xz_xx = buffer_2000_ppdd[876];

    auto g_xz_0_0_0_z_x_xz_xy = buffer_2000_ppdd[877];

    auto g_xz_0_0_0_z_x_xz_xz = buffer_2000_ppdd[878];

    auto g_xz_0_0_0_z_x_xz_yy = buffer_2000_ppdd[879];

    auto g_xz_0_0_0_z_x_xz_yz = buffer_2000_ppdd[880];

    auto g_xz_0_0_0_z_x_xz_zz = buffer_2000_ppdd[881];

    auto g_xz_0_0_0_z_x_yy_xx = buffer_2000_ppdd[882];

    auto g_xz_0_0_0_z_x_yy_xy = buffer_2000_ppdd[883];

    auto g_xz_0_0_0_z_x_yy_xz = buffer_2000_ppdd[884];

    auto g_xz_0_0_0_z_x_yy_yy = buffer_2000_ppdd[885];

    auto g_xz_0_0_0_z_x_yy_yz = buffer_2000_ppdd[886];

    auto g_xz_0_0_0_z_x_yy_zz = buffer_2000_ppdd[887];

    auto g_xz_0_0_0_z_x_yz_xx = buffer_2000_ppdd[888];

    auto g_xz_0_0_0_z_x_yz_xy = buffer_2000_ppdd[889];

    auto g_xz_0_0_0_z_x_yz_xz = buffer_2000_ppdd[890];

    auto g_xz_0_0_0_z_x_yz_yy = buffer_2000_ppdd[891];

    auto g_xz_0_0_0_z_x_yz_yz = buffer_2000_ppdd[892];

    auto g_xz_0_0_0_z_x_yz_zz = buffer_2000_ppdd[893];

    auto g_xz_0_0_0_z_x_zz_xx = buffer_2000_ppdd[894];

    auto g_xz_0_0_0_z_x_zz_xy = buffer_2000_ppdd[895];

    auto g_xz_0_0_0_z_x_zz_xz = buffer_2000_ppdd[896];

    auto g_xz_0_0_0_z_x_zz_yy = buffer_2000_ppdd[897];

    auto g_xz_0_0_0_z_x_zz_yz = buffer_2000_ppdd[898];

    auto g_xz_0_0_0_z_x_zz_zz = buffer_2000_ppdd[899];

    auto g_xz_0_0_0_z_y_xx_xx = buffer_2000_ppdd[900];

    auto g_xz_0_0_0_z_y_xx_xy = buffer_2000_ppdd[901];

    auto g_xz_0_0_0_z_y_xx_xz = buffer_2000_ppdd[902];

    auto g_xz_0_0_0_z_y_xx_yy = buffer_2000_ppdd[903];

    auto g_xz_0_0_0_z_y_xx_yz = buffer_2000_ppdd[904];

    auto g_xz_0_0_0_z_y_xx_zz = buffer_2000_ppdd[905];

    auto g_xz_0_0_0_z_y_xy_xx = buffer_2000_ppdd[906];

    auto g_xz_0_0_0_z_y_xy_xy = buffer_2000_ppdd[907];

    auto g_xz_0_0_0_z_y_xy_xz = buffer_2000_ppdd[908];

    auto g_xz_0_0_0_z_y_xy_yy = buffer_2000_ppdd[909];

    auto g_xz_0_0_0_z_y_xy_yz = buffer_2000_ppdd[910];

    auto g_xz_0_0_0_z_y_xy_zz = buffer_2000_ppdd[911];

    auto g_xz_0_0_0_z_y_xz_xx = buffer_2000_ppdd[912];

    auto g_xz_0_0_0_z_y_xz_xy = buffer_2000_ppdd[913];

    auto g_xz_0_0_0_z_y_xz_xz = buffer_2000_ppdd[914];

    auto g_xz_0_0_0_z_y_xz_yy = buffer_2000_ppdd[915];

    auto g_xz_0_0_0_z_y_xz_yz = buffer_2000_ppdd[916];

    auto g_xz_0_0_0_z_y_xz_zz = buffer_2000_ppdd[917];

    auto g_xz_0_0_0_z_y_yy_xx = buffer_2000_ppdd[918];

    auto g_xz_0_0_0_z_y_yy_xy = buffer_2000_ppdd[919];

    auto g_xz_0_0_0_z_y_yy_xz = buffer_2000_ppdd[920];

    auto g_xz_0_0_0_z_y_yy_yy = buffer_2000_ppdd[921];

    auto g_xz_0_0_0_z_y_yy_yz = buffer_2000_ppdd[922];

    auto g_xz_0_0_0_z_y_yy_zz = buffer_2000_ppdd[923];

    auto g_xz_0_0_0_z_y_yz_xx = buffer_2000_ppdd[924];

    auto g_xz_0_0_0_z_y_yz_xy = buffer_2000_ppdd[925];

    auto g_xz_0_0_0_z_y_yz_xz = buffer_2000_ppdd[926];

    auto g_xz_0_0_0_z_y_yz_yy = buffer_2000_ppdd[927];

    auto g_xz_0_0_0_z_y_yz_yz = buffer_2000_ppdd[928];

    auto g_xz_0_0_0_z_y_yz_zz = buffer_2000_ppdd[929];

    auto g_xz_0_0_0_z_y_zz_xx = buffer_2000_ppdd[930];

    auto g_xz_0_0_0_z_y_zz_xy = buffer_2000_ppdd[931];

    auto g_xz_0_0_0_z_y_zz_xz = buffer_2000_ppdd[932];

    auto g_xz_0_0_0_z_y_zz_yy = buffer_2000_ppdd[933];

    auto g_xz_0_0_0_z_y_zz_yz = buffer_2000_ppdd[934];

    auto g_xz_0_0_0_z_y_zz_zz = buffer_2000_ppdd[935];

    auto g_xz_0_0_0_z_z_xx_xx = buffer_2000_ppdd[936];

    auto g_xz_0_0_0_z_z_xx_xy = buffer_2000_ppdd[937];

    auto g_xz_0_0_0_z_z_xx_xz = buffer_2000_ppdd[938];

    auto g_xz_0_0_0_z_z_xx_yy = buffer_2000_ppdd[939];

    auto g_xz_0_0_0_z_z_xx_yz = buffer_2000_ppdd[940];

    auto g_xz_0_0_0_z_z_xx_zz = buffer_2000_ppdd[941];

    auto g_xz_0_0_0_z_z_xy_xx = buffer_2000_ppdd[942];

    auto g_xz_0_0_0_z_z_xy_xy = buffer_2000_ppdd[943];

    auto g_xz_0_0_0_z_z_xy_xz = buffer_2000_ppdd[944];

    auto g_xz_0_0_0_z_z_xy_yy = buffer_2000_ppdd[945];

    auto g_xz_0_0_0_z_z_xy_yz = buffer_2000_ppdd[946];

    auto g_xz_0_0_0_z_z_xy_zz = buffer_2000_ppdd[947];

    auto g_xz_0_0_0_z_z_xz_xx = buffer_2000_ppdd[948];

    auto g_xz_0_0_0_z_z_xz_xy = buffer_2000_ppdd[949];

    auto g_xz_0_0_0_z_z_xz_xz = buffer_2000_ppdd[950];

    auto g_xz_0_0_0_z_z_xz_yy = buffer_2000_ppdd[951];

    auto g_xz_0_0_0_z_z_xz_yz = buffer_2000_ppdd[952];

    auto g_xz_0_0_0_z_z_xz_zz = buffer_2000_ppdd[953];

    auto g_xz_0_0_0_z_z_yy_xx = buffer_2000_ppdd[954];

    auto g_xz_0_0_0_z_z_yy_xy = buffer_2000_ppdd[955];

    auto g_xz_0_0_0_z_z_yy_xz = buffer_2000_ppdd[956];

    auto g_xz_0_0_0_z_z_yy_yy = buffer_2000_ppdd[957];

    auto g_xz_0_0_0_z_z_yy_yz = buffer_2000_ppdd[958];

    auto g_xz_0_0_0_z_z_yy_zz = buffer_2000_ppdd[959];

    auto g_xz_0_0_0_z_z_yz_xx = buffer_2000_ppdd[960];

    auto g_xz_0_0_0_z_z_yz_xy = buffer_2000_ppdd[961];

    auto g_xz_0_0_0_z_z_yz_xz = buffer_2000_ppdd[962];

    auto g_xz_0_0_0_z_z_yz_yy = buffer_2000_ppdd[963];

    auto g_xz_0_0_0_z_z_yz_yz = buffer_2000_ppdd[964];

    auto g_xz_0_0_0_z_z_yz_zz = buffer_2000_ppdd[965];

    auto g_xz_0_0_0_z_z_zz_xx = buffer_2000_ppdd[966];

    auto g_xz_0_0_0_z_z_zz_xy = buffer_2000_ppdd[967];

    auto g_xz_0_0_0_z_z_zz_xz = buffer_2000_ppdd[968];

    auto g_xz_0_0_0_z_z_zz_yy = buffer_2000_ppdd[969];

    auto g_xz_0_0_0_z_z_zz_yz = buffer_2000_ppdd[970];

    auto g_xz_0_0_0_z_z_zz_zz = buffer_2000_ppdd[971];

    auto g_yy_0_0_0_x_x_xx_xx = buffer_2000_ppdd[972];

    auto g_yy_0_0_0_x_x_xx_xy = buffer_2000_ppdd[973];

    auto g_yy_0_0_0_x_x_xx_xz = buffer_2000_ppdd[974];

    auto g_yy_0_0_0_x_x_xx_yy = buffer_2000_ppdd[975];

    auto g_yy_0_0_0_x_x_xx_yz = buffer_2000_ppdd[976];

    auto g_yy_0_0_0_x_x_xx_zz = buffer_2000_ppdd[977];

    auto g_yy_0_0_0_x_x_xy_xx = buffer_2000_ppdd[978];

    auto g_yy_0_0_0_x_x_xy_xy = buffer_2000_ppdd[979];

    auto g_yy_0_0_0_x_x_xy_xz = buffer_2000_ppdd[980];

    auto g_yy_0_0_0_x_x_xy_yy = buffer_2000_ppdd[981];

    auto g_yy_0_0_0_x_x_xy_yz = buffer_2000_ppdd[982];

    auto g_yy_0_0_0_x_x_xy_zz = buffer_2000_ppdd[983];

    auto g_yy_0_0_0_x_x_xz_xx = buffer_2000_ppdd[984];

    auto g_yy_0_0_0_x_x_xz_xy = buffer_2000_ppdd[985];

    auto g_yy_0_0_0_x_x_xz_xz = buffer_2000_ppdd[986];

    auto g_yy_0_0_0_x_x_xz_yy = buffer_2000_ppdd[987];

    auto g_yy_0_0_0_x_x_xz_yz = buffer_2000_ppdd[988];

    auto g_yy_0_0_0_x_x_xz_zz = buffer_2000_ppdd[989];

    auto g_yy_0_0_0_x_x_yy_xx = buffer_2000_ppdd[990];

    auto g_yy_0_0_0_x_x_yy_xy = buffer_2000_ppdd[991];

    auto g_yy_0_0_0_x_x_yy_xz = buffer_2000_ppdd[992];

    auto g_yy_0_0_0_x_x_yy_yy = buffer_2000_ppdd[993];

    auto g_yy_0_0_0_x_x_yy_yz = buffer_2000_ppdd[994];

    auto g_yy_0_0_0_x_x_yy_zz = buffer_2000_ppdd[995];

    auto g_yy_0_0_0_x_x_yz_xx = buffer_2000_ppdd[996];

    auto g_yy_0_0_0_x_x_yz_xy = buffer_2000_ppdd[997];

    auto g_yy_0_0_0_x_x_yz_xz = buffer_2000_ppdd[998];

    auto g_yy_0_0_0_x_x_yz_yy = buffer_2000_ppdd[999];

    auto g_yy_0_0_0_x_x_yz_yz = buffer_2000_ppdd[1000];

    auto g_yy_0_0_0_x_x_yz_zz = buffer_2000_ppdd[1001];

    auto g_yy_0_0_0_x_x_zz_xx = buffer_2000_ppdd[1002];

    auto g_yy_0_0_0_x_x_zz_xy = buffer_2000_ppdd[1003];

    auto g_yy_0_0_0_x_x_zz_xz = buffer_2000_ppdd[1004];

    auto g_yy_0_0_0_x_x_zz_yy = buffer_2000_ppdd[1005];

    auto g_yy_0_0_0_x_x_zz_yz = buffer_2000_ppdd[1006];

    auto g_yy_0_0_0_x_x_zz_zz = buffer_2000_ppdd[1007];

    auto g_yy_0_0_0_x_y_xx_xx = buffer_2000_ppdd[1008];

    auto g_yy_0_0_0_x_y_xx_xy = buffer_2000_ppdd[1009];

    auto g_yy_0_0_0_x_y_xx_xz = buffer_2000_ppdd[1010];

    auto g_yy_0_0_0_x_y_xx_yy = buffer_2000_ppdd[1011];

    auto g_yy_0_0_0_x_y_xx_yz = buffer_2000_ppdd[1012];

    auto g_yy_0_0_0_x_y_xx_zz = buffer_2000_ppdd[1013];

    auto g_yy_0_0_0_x_y_xy_xx = buffer_2000_ppdd[1014];

    auto g_yy_0_0_0_x_y_xy_xy = buffer_2000_ppdd[1015];

    auto g_yy_0_0_0_x_y_xy_xz = buffer_2000_ppdd[1016];

    auto g_yy_0_0_0_x_y_xy_yy = buffer_2000_ppdd[1017];

    auto g_yy_0_0_0_x_y_xy_yz = buffer_2000_ppdd[1018];

    auto g_yy_0_0_0_x_y_xy_zz = buffer_2000_ppdd[1019];

    auto g_yy_0_0_0_x_y_xz_xx = buffer_2000_ppdd[1020];

    auto g_yy_0_0_0_x_y_xz_xy = buffer_2000_ppdd[1021];

    auto g_yy_0_0_0_x_y_xz_xz = buffer_2000_ppdd[1022];

    auto g_yy_0_0_0_x_y_xz_yy = buffer_2000_ppdd[1023];

    auto g_yy_0_0_0_x_y_xz_yz = buffer_2000_ppdd[1024];

    auto g_yy_0_0_0_x_y_xz_zz = buffer_2000_ppdd[1025];

    auto g_yy_0_0_0_x_y_yy_xx = buffer_2000_ppdd[1026];

    auto g_yy_0_0_0_x_y_yy_xy = buffer_2000_ppdd[1027];

    auto g_yy_0_0_0_x_y_yy_xz = buffer_2000_ppdd[1028];

    auto g_yy_0_0_0_x_y_yy_yy = buffer_2000_ppdd[1029];

    auto g_yy_0_0_0_x_y_yy_yz = buffer_2000_ppdd[1030];

    auto g_yy_0_0_0_x_y_yy_zz = buffer_2000_ppdd[1031];

    auto g_yy_0_0_0_x_y_yz_xx = buffer_2000_ppdd[1032];

    auto g_yy_0_0_0_x_y_yz_xy = buffer_2000_ppdd[1033];

    auto g_yy_0_0_0_x_y_yz_xz = buffer_2000_ppdd[1034];

    auto g_yy_0_0_0_x_y_yz_yy = buffer_2000_ppdd[1035];

    auto g_yy_0_0_0_x_y_yz_yz = buffer_2000_ppdd[1036];

    auto g_yy_0_0_0_x_y_yz_zz = buffer_2000_ppdd[1037];

    auto g_yy_0_0_0_x_y_zz_xx = buffer_2000_ppdd[1038];

    auto g_yy_0_0_0_x_y_zz_xy = buffer_2000_ppdd[1039];

    auto g_yy_0_0_0_x_y_zz_xz = buffer_2000_ppdd[1040];

    auto g_yy_0_0_0_x_y_zz_yy = buffer_2000_ppdd[1041];

    auto g_yy_0_0_0_x_y_zz_yz = buffer_2000_ppdd[1042];

    auto g_yy_0_0_0_x_y_zz_zz = buffer_2000_ppdd[1043];

    auto g_yy_0_0_0_x_z_xx_xx = buffer_2000_ppdd[1044];

    auto g_yy_0_0_0_x_z_xx_xy = buffer_2000_ppdd[1045];

    auto g_yy_0_0_0_x_z_xx_xz = buffer_2000_ppdd[1046];

    auto g_yy_0_0_0_x_z_xx_yy = buffer_2000_ppdd[1047];

    auto g_yy_0_0_0_x_z_xx_yz = buffer_2000_ppdd[1048];

    auto g_yy_0_0_0_x_z_xx_zz = buffer_2000_ppdd[1049];

    auto g_yy_0_0_0_x_z_xy_xx = buffer_2000_ppdd[1050];

    auto g_yy_0_0_0_x_z_xy_xy = buffer_2000_ppdd[1051];

    auto g_yy_0_0_0_x_z_xy_xz = buffer_2000_ppdd[1052];

    auto g_yy_0_0_0_x_z_xy_yy = buffer_2000_ppdd[1053];

    auto g_yy_0_0_0_x_z_xy_yz = buffer_2000_ppdd[1054];

    auto g_yy_0_0_0_x_z_xy_zz = buffer_2000_ppdd[1055];

    auto g_yy_0_0_0_x_z_xz_xx = buffer_2000_ppdd[1056];

    auto g_yy_0_0_0_x_z_xz_xy = buffer_2000_ppdd[1057];

    auto g_yy_0_0_0_x_z_xz_xz = buffer_2000_ppdd[1058];

    auto g_yy_0_0_0_x_z_xz_yy = buffer_2000_ppdd[1059];

    auto g_yy_0_0_0_x_z_xz_yz = buffer_2000_ppdd[1060];

    auto g_yy_0_0_0_x_z_xz_zz = buffer_2000_ppdd[1061];

    auto g_yy_0_0_0_x_z_yy_xx = buffer_2000_ppdd[1062];

    auto g_yy_0_0_0_x_z_yy_xy = buffer_2000_ppdd[1063];

    auto g_yy_0_0_0_x_z_yy_xz = buffer_2000_ppdd[1064];

    auto g_yy_0_0_0_x_z_yy_yy = buffer_2000_ppdd[1065];

    auto g_yy_0_0_0_x_z_yy_yz = buffer_2000_ppdd[1066];

    auto g_yy_0_0_0_x_z_yy_zz = buffer_2000_ppdd[1067];

    auto g_yy_0_0_0_x_z_yz_xx = buffer_2000_ppdd[1068];

    auto g_yy_0_0_0_x_z_yz_xy = buffer_2000_ppdd[1069];

    auto g_yy_0_0_0_x_z_yz_xz = buffer_2000_ppdd[1070];

    auto g_yy_0_0_0_x_z_yz_yy = buffer_2000_ppdd[1071];

    auto g_yy_0_0_0_x_z_yz_yz = buffer_2000_ppdd[1072];

    auto g_yy_0_0_0_x_z_yz_zz = buffer_2000_ppdd[1073];

    auto g_yy_0_0_0_x_z_zz_xx = buffer_2000_ppdd[1074];

    auto g_yy_0_0_0_x_z_zz_xy = buffer_2000_ppdd[1075];

    auto g_yy_0_0_0_x_z_zz_xz = buffer_2000_ppdd[1076];

    auto g_yy_0_0_0_x_z_zz_yy = buffer_2000_ppdd[1077];

    auto g_yy_0_0_0_x_z_zz_yz = buffer_2000_ppdd[1078];

    auto g_yy_0_0_0_x_z_zz_zz = buffer_2000_ppdd[1079];

    auto g_yy_0_0_0_y_x_xx_xx = buffer_2000_ppdd[1080];

    auto g_yy_0_0_0_y_x_xx_xy = buffer_2000_ppdd[1081];

    auto g_yy_0_0_0_y_x_xx_xz = buffer_2000_ppdd[1082];

    auto g_yy_0_0_0_y_x_xx_yy = buffer_2000_ppdd[1083];

    auto g_yy_0_0_0_y_x_xx_yz = buffer_2000_ppdd[1084];

    auto g_yy_0_0_0_y_x_xx_zz = buffer_2000_ppdd[1085];

    auto g_yy_0_0_0_y_x_xy_xx = buffer_2000_ppdd[1086];

    auto g_yy_0_0_0_y_x_xy_xy = buffer_2000_ppdd[1087];

    auto g_yy_0_0_0_y_x_xy_xz = buffer_2000_ppdd[1088];

    auto g_yy_0_0_0_y_x_xy_yy = buffer_2000_ppdd[1089];

    auto g_yy_0_0_0_y_x_xy_yz = buffer_2000_ppdd[1090];

    auto g_yy_0_0_0_y_x_xy_zz = buffer_2000_ppdd[1091];

    auto g_yy_0_0_0_y_x_xz_xx = buffer_2000_ppdd[1092];

    auto g_yy_0_0_0_y_x_xz_xy = buffer_2000_ppdd[1093];

    auto g_yy_0_0_0_y_x_xz_xz = buffer_2000_ppdd[1094];

    auto g_yy_0_0_0_y_x_xz_yy = buffer_2000_ppdd[1095];

    auto g_yy_0_0_0_y_x_xz_yz = buffer_2000_ppdd[1096];

    auto g_yy_0_0_0_y_x_xz_zz = buffer_2000_ppdd[1097];

    auto g_yy_0_0_0_y_x_yy_xx = buffer_2000_ppdd[1098];

    auto g_yy_0_0_0_y_x_yy_xy = buffer_2000_ppdd[1099];

    auto g_yy_0_0_0_y_x_yy_xz = buffer_2000_ppdd[1100];

    auto g_yy_0_0_0_y_x_yy_yy = buffer_2000_ppdd[1101];

    auto g_yy_0_0_0_y_x_yy_yz = buffer_2000_ppdd[1102];

    auto g_yy_0_0_0_y_x_yy_zz = buffer_2000_ppdd[1103];

    auto g_yy_0_0_0_y_x_yz_xx = buffer_2000_ppdd[1104];

    auto g_yy_0_0_0_y_x_yz_xy = buffer_2000_ppdd[1105];

    auto g_yy_0_0_0_y_x_yz_xz = buffer_2000_ppdd[1106];

    auto g_yy_0_0_0_y_x_yz_yy = buffer_2000_ppdd[1107];

    auto g_yy_0_0_0_y_x_yz_yz = buffer_2000_ppdd[1108];

    auto g_yy_0_0_0_y_x_yz_zz = buffer_2000_ppdd[1109];

    auto g_yy_0_0_0_y_x_zz_xx = buffer_2000_ppdd[1110];

    auto g_yy_0_0_0_y_x_zz_xy = buffer_2000_ppdd[1111];

    auto g_yy_0_0_0_y_x_zz_xz = buffer_2000_ppdd[1112];

    auto g_yy_0_0_0_y_x_zz_yy = buffer_2000_ppdd[1113];

    auto g_yy_0_0_0_y_x_zz_yz = buffer_2000_ppdd[1114];

    auto g_yy_0_0_0_y_x_zz_zz = buffer_2000_ppdd[1115];

    auto g_yy_0_0_0_y_y_xx_xx = buffer_2000_ppdd[1116];

    auto g_yy_0_0_0_y_y_xx_xy = buffer_2000_ppdd[1117];

    auto g_yy_0_0_0_y_y_xx_xz = buffer_2000_ppdd[1118];

    auto g_yy_0_0_0_y_y_xx_yy = buffer_2000_ppdd[1119];

    auto g_yy_0_0_0_y_y_xx_yz = buffer_2000_ppdd[1120];

    auto g_yy_0_0_0_y_y_xx_zz = buffer_2000_ppdd[1121];

    auto g_yy_0_0_0_y_y_xy_xx = buffer_2000_ppdd[1122];

    auto g_yy_0_0_0_y_y_xy_xy = buffer_2000_ppdd[1123];

    auto g_yy_0_0_0_y_y_xy_xz = buffer_2000_ppdd[1124];

    auto g_yy_0_0_0_y_y_xy_yy = buffer_2000_ppdd[1125];

    auto g_yy_0_0_0_y_y_xy_yz = buffer_2000_ppdd[1126];

    auto g_yy_0_0_0_y_y_xy_zz = buffer_2000_ppdd[1127];

    auto g_yy_0_0_0_y_y_xz_xx = buffer_2000_ppdd[1128];

    auto g_yy_0_0_0_y_y_xz_xy = buffer_2000_ppdd[1129];

    auto g_yy_0_0_0_y_y_xz_xz = buffer_2000_ppdd[1130];

    auto g_yy_0_0_0_y_y_xz_yy = buffer_2000_ppdd[1131];

    auto g_yy_0_0_0_y_y_xz_yz = buffer_2000_ppdd[1132];

    auto g_yy_0_0_0_y_y_xz_zz = buffer_2000_ppdd[1133];

    auto g_yy_0_0_0_y_y_yy_xx = buffer_2000_ppdd[1134];

    auto g_yy_0_0_0_y_y_yy_xy = buffer_2000_ppdd[1135];

    auto g_yy_0_0_0_y_y_yy_xz = buffer_2000_ppdd[1136];

    auto g_yy_0_0_0_y_y_yy_yy = buffer_2000_ppdd[1137];

    auto g_yy_0_0_0_y_y_yy_yz = buffer_2000_ppdd[1138];

    auto g_yy_0_0_0_y_y_yy_zz = buffer_2000_ppdd[1139];

    auto g_yy_0_0_0_y_y_yz_xx = buffer_2000_ppdd[1140];

    auto g_yy_0_0_0_y_y_yz_xy = buffer_2000_ppdd[1141];

    auto g_yy_0_0_0_y_y_yz_xz = buffer_2000_ppdd[1142];

    auto g_yy_0_0_0_y_y_yz_yy = buffer_2000_ppdd[1143];

    auto g_yy_0_0_0_y_y_yz_yz = buffer_2000_ppdd[1144];

    auto g_yy_0_0_0_y_y_yz_zz = buffer_2000_ppdd[1145];

    auto g_yy_0_0_0_y_y_zz_xx = buffer_2000_ppdd[1146];

    auto g_yy_0_0_0_y_y_zz_xy = buffer_2000_ppdd[1147];

    auto g_yy_0_0_0_y_y_zz_xz = buffer_2000_ppdd[1148];

    auto g_yy_0_0_0_y_y_zz_yy = buffer_2000_ppdd[1149];

    auto g_yy_0_0_0_y_y_zz_yz = buffer_2000_ppdd[1150];

    auto g_yy_0_0_0_y_y_zz_zz = buffer_2000_ppdd[1151];

    auto g_yy_0_0_0_y_z_xx_xx = buffer_2000_ppdd[1152];

    auto g_yy_0_0_0_y_z_xx_xy = buffer_2000_ppdd[1153];

    auto g_yy_0_0_0_y_z_xx_xz = buffer_2000_ppdd[1154];

    auto g_yy_0_0_0_y_z_xx_yy = buffer_2000_ppdd[1155];

    auto g_yy_0_0_0_y_z_xx_yz = buffer_2000_ppdd[1156];

    auto g_yy_0_0_0_y_z_xx_zz = buffer_2000_ppdd[1157];

    auto g_yy_0_0_0_y_z_xy_xx = buffer_2000_ppdd[1158];

    auto g_yy_0_0_0_y_z_xy_xy = buffer_2000_ppdd[1159];

    auto g_yy_0_0_0_y_z_xy_xz = buffer_2000_ppdd[1160];

    auto g_yy_0_0_0_y_z_xy_yy = buffer_2000_ppdd[1161];

    auto g_yy_0_0_0_y_z_xy_yz = buffer_2000_ppdd[1162];

    auto g_yy_0_0_0_y_z_xy_zz = buffer_2000_ppdd[1163];

    auto g_yy_0_0_0_y_z_xz_xx = buffer_2000_ppdd[1164];

    auto g_yy_0_0_0_y_z_xz_xy = buffer_2000_ppdd[1165];

    auto g_yy_0_0_0_y_z_xz_xz = buffer_2000_ppdd[1166];

    auto g_yy_0_0_0_y_z_xz_yy = buffer_2000_ppdd[1167];

    auto g_yy_0_0_0_y_z_xz_yz = buffer_2000_ppdd[1168];

    auto g_yy_0_0_0_y_z_xz_zz = buffer_2000_ppdd[1169];

    auto g_yy_0_0_0_y_z_yy_xx = buffer_2000_ppdd[1170];

    auto g_yy_0_0_0_y_z_yy_xy = buffer_2000_ppdd[1171];

    auto g_yy_0_0_0_y_z_yy_xz = buffer_2000_ppdd[1172];

    auto g_yy_0_0_0_y_z_yy_yy = buffer_2000_ppdd[1173];

    auto g_yy_0_0_0_y_z_yy_yz = buffer_2000_ppdd[1174];

    auto g_yy_0_0_0_y_z_yy_zz = buffer_2000_ppdd[1175];

    auto g_yy_0_0_0_y_z_yz_xx = buffer_2000_ppdd[1176];

    auto g_yy_0_0_0_y_z_yz_xy = buffer_2000_ppdd[1177];

    auto g_yy_0_0_0_y_z_yz_xz = buffer_2000_ppdd[1178];

    auto g_yy_0_0_0_y_z_yz_yy = buffer_2000_ppdd[1179];

    auto g_yy_0_0_0_y_z_yz_yz = buffer_2000_ppdd[1180];

    auto g_yy_0_0_0_y_z_yz_zz = buffer_2000_ppdd[1181];

    auto g_yy_0_0_0_y_z_zz_xx = buffer_2000_ppdd[1182];

    auto g_yy_0_0_0_y_z_zz_xy = buffer_2000_ppdd[1183];

    auto g_yy_0_0_0_y_z_zz_xz = buffer_2000_ppdd[1184];

    auto g_yy_0_0_0_y_z_zz_yy = buffer_2000_ppdd[1185];

    auto g_yy_0_0_0_y_z_zz_yz = buffer_2000_ppdd[1186];

    auto g_yy_0_0_0_y_z_zz_zz = buffer_2000_ppdd[1187];

    auto g_yy_0_0_0_z_x_xx_xx = buffer_2000_ppdd[1188];

    auto g_yy_0_0_0_z_x_xx_xy = buffer_2000_ppdd[1189];

    auto g_yy_0_0_0_z_x_xx_xz = buffer_2000_ppdd[1190];

    auto g_yy_0_0_0_z_x_xx_yy = buffer_2000_ppdd[1191];

    auto g_yy_0_0_0_z_x_xx_yz = buffer_2000_ppdd[1192];

    auto g_yy_0_0_0_z_x_xx_zz = buffer_2000_ppdd[1193];

    auto g_yy_0_0_0_z_x_xy_xx = buffer_2000_ppdd[1194];

    auto g_yy_0_0_0_z_x_xy_xy = buffer_2000_ppdd[1195];

    auto g_yy_0_0_0_z_x_xy_xz = buffer_2000_ppdd[1196];

    auto g_yy_0_0_0_z_x_xy_yy = buffer_2000_ppdd[1197];

    auto g_yy_0_0_0_z_x_xy_yz = buffer_2000_ppdd[1198];

    auto g_yy_0_0_0_z_x_xy_zz = buffer_2000_ppdd[1199];

    auto g_yy_0_0_0_z_x_xz_xx = buffer_2000_ppdd[1200];

    auto g_yy_0_0_0_z_x_xz_xy = buffer_2000_ppdd[1201];

    auto g_yy_0_0_0_z_x_xz_xz = buffer_2000_ppdd[1202];

    auto g_yy_0_0_0_z_x_xz_yy = buffer_2000_ppdd[1203];

    auto g_yy_0_0_0_z_x_xz_yz = buffer_2000_ppdd[1204];

    auto g_yy_0_0_0_z_x_xz_zz = buffer_2000_ppdd[1205];

    auto g_yy_0_0_0_z_x_yy_xx = buffer_2000_ppdd[1206];

    auto g_yy_0_0_0_z_x_yy_xy = buffer_2000_ppdd[1207];

    auto g_yy_0_0_0_z_x_yy_xz = buffer_2000_ppdd[1208];

    auto g_yy_0_0_0_z_x_yy_yy = buffer_2000_ppdd[1209];

    auto g_yy_0_0_0_z_x_yy_yz = buffer_2000_ppdd[1210];

    auto g_yy_0_0_0_z_x_yy_zz = buffer_2000_ppdd[1211];

    auto g_yy_0_0_0_z_x_yz_xx = buffer_2000_ppdd[1212];

    auto g_yy_0_0_0_z_x_yz_xy = buffer_2000_ppdd[1213];

    auto g_yy_0_0_0_z_x_yz_xz = buffer_2000_ppdd[1214];

    auto g_yy_0_0_0_z_x_yz_yy = buffer_2000_ppdd[1215];

    auto g_yy_0_0_0_z_x_yz_yz = buffer_2000_ppdd[1216];

    auto g_yy_0_0_0_z_x_yz_zz = buffer_2000_ppdd[1217];

    auto g_yy_0_0_0_z_x_zz_xx = buffer_2000_ppdd[1218];

    auto g_yy_0_0_0_z_x_zz_xy = buffer_2000_ppdd[1219];

    auto g_yy_0_0_0_z_x_zz_xz = buffer_2000_ppdd[1220];

    auto g_yy_0_0_0_z_x_zz_yy = buffer_2000_ppdd[1221];

    auto g_yy_0_0_0_z_x_zz_yz = buffer_2000_ppdd[1222];

    auto g_yy_0_0_0_z_x_zz_zz = buffer_2000_ppdd[1223];

    auto g_yy_0_0_0_z_y_xx_xx = buffer_2000_ppdd[1224];

    auto g_yy_0_0_0_z_y_xx_xy = buffer_2000_ppdd[1225];

    auto g_yy_0_0_0_z_y_xx_xz = buffer_2000_ppdd[1226];

    auto g_yy_0_0_0_z_y_xx_yy = buffer_2000_ppdd[1227];

    auto g_yy_0_0_0_z_y_xx_yz = buffer_2000_ppdd[1228];

    auto g_yy_0_0_0_z_y_xx_zz = buffer_2000_ppdd[1229];

    auto g_yy_0_0_0_z_y_xy_xx = buffer_2000_ppdd[1230];

    auto g_yy_0_0_0_z_y_xy_xy = buffer_2000_ppdd[1231];

    auto g_yy_0_0_0_z_y_xy_xz = buffer_2000_ppdd[1232];

    auto g_yy_0_0_0_z_y_xy_yy = buffer_2000_ppdd[1233];

    auto g_yy_0_0_0_z_y_xy_yz = buffer_2000_ppdd[1234];

    auto g_yy_0_0_0_z_y_xy_zz = buffer_2000_ppdd[1235];

    auto g_yy_0_0_0_z_y_xz_xx = buffer_2000_ppdd[1236];

    auto g_yy_0_0_0_z_y_xz_xy = buffer_2000_ppdd[1237];

    auto g_yy_0_0_0_z_y_xz_xz = buffer_2000_ppdd[1238];

    auto g_yy_0_0_0_z_y_xz_yy = buffer_2000_ppdd[1239];

    auto g_yy_0_0_0_z_y_xz_yz = buffer_2000_ppdd[1240];

    auto g_yy_0_0_0_z_y_xz_zz = buffer_2000_ppdd[1241];

    auto g_yy_0_0_0_z_y_yy_xx = buffer_2000_ppdd[1242];

    auto g_yy_0_0_0_z_y_yy_xy = buffer_2000_ppdd[1243];

    auto g_yy_0_0_0_z_y_yy_xz = buffer_2000_ppdd[1244];

    auto g_yy_0_0_0_z_y_yy_yy = buffer_2000_ppdd[1245];

    auto g_yy_0_0_0_z_y_yy_yz = buffer_2000_ppdd[1246];

    auto g_yy_0_0_0_z_y_yy_zz = buffer_2000_ppdd[1247];

    auto g_yy_0_0_0_z_y_yz_xx = buffer_2000_ppdd[1248];

    auto g_yy_0_0_0_z_y_yz_xy = buffer_2000_ppdd[1249];

    auto g_yy_0_0_0_z_y_yz_xz = buffer_2000_ppdd[1250];

    auto g_yy_0_0_0_z_y_yz_yy = buffer_2000_ppdd[1251];

    auto g_yy_0_0_0_z_y_yz_yz = buffer_2000_ppdd[1252];

    auto g_yy_0_0_0_z_y_yz_zz = buffer_2000_ppdd[1253];

    auto g_yy_0_0_0_z_y_zz_xx = buffer_2000_ppdd[1254];

    auto g_yy_0_0_0_z_y_zz_xy = buffer_2000_ppdd[1255];

    auto g_yy_0_0_0_z_y_zz_xz = buffer_2000_ppdd[1256];

    auto g_yy_0_0_0_z_y_zz_yy = buffer_2000_ppdd[1257];

    auto g_yy_0_0_0_z_y_zz_yz = buffer_2000_ppdd[1258];

    auto g_yy_0_0_0_z_y_zz_zz = buffer_2000_ppdd[1259];

    auto g_yy_0_0_0_z_z_xx_xx = buffer_2000_ppdd[1260];

    auto g_yy_0_0_0_z_z_xx_xy = buffer_2000_ppdd[1261];

    auto g_yy_0_0_0_z_z_xx_xz = buffer_2000_ppdd[1262];

    auto g_yy_0_0_0_z_z_xx_yy = buffer_2000_ppdd[1263];

    auto g_yy_0_0_0_z_z_xx_yz = buffer_2000_ppdd[1264];

    auto g_yy_0_0_0_z_z_xx_zz = buffer_2000_ppdd[1265];

    auto g_yy_0_0_0_z_z_xy_xx = buffer_2000_ppdd[1266];

    auto g_yy_0_0_0_z_z_xy_xy = buffer_2000_ppdd[1267];

    auto g_yy_0_0_0_z_z_xy_xz = buffer_2000_ppdd[1268];

    auto g_yy_0_0_0_z_z_xy_yy = buffer_2000_ppdd[1269];

    auto g_yy_0_0_0_z_z_xy_yz = buffer_2000_ppdd[1270];

    auto g_yy_0_0_0_z_z_xy_zz = buffer_2000_ppdd[1271];

    auto g_yy_0_0_0_z_z_xz_xx = buffer_2000_ppdd[1272];

    auto g_yy_0_0_0_z_z_xz_xy = buffer_2000_ppdd[1273];

    auto g_yy_0_0_0_z_z_xz_xz = buffer_2000_ppdd[1274];

    auto g_yy_0_0_0_z_z_xz_yy = buffer_2000_ppdd[1275];

    auto g_yy_0_0_0_z_z_xz_yz = buffer_2000_ppdd[1276];

    auto g_yy_0_0_0_z_z_xz_zz = buffer_2000_ppdd[1277];

    auto g_yy_0_0_0_z_z_yy_xx = buffer_2000_ppdd[1278];

    auto g_yy_0_0_0_z_z_yy_xy = buffer_2000_ppdd[1279];

    auto g_yy_0_0_0_z_z_yy_xz = buffer_2000_ppdd[1280];

    auto g_yy_0_0_0_z_z_yy_yy = buffer_2000_ppdd[1281];

    auto g_yy_0_0_0_z_z_yy_yz = buffer_2000_ppdd[1282];

    auto g_yy_0_0_0_z_z_yy_zz = buffer_2000_ppdd[1283];

    auto g_yy_0_0_0_z_z_yz_xx = buffer_2000_ppdd[1284];

    auto g_yy_0_0_0_z_z_yz_xy = buffer_2000_ppdd[1285];

    auto g_yy_0_0_0_z_z_yz_xz = buffer_2000_ppdd[1286];

    auto g_yy_0_0_0_z_z_yz_yy = buffer_2000_ppdd[1287];

    auto g_yy_0_0_0_z_z_yz_yz = buffer_2000_ppdd[1288];

    auto g_yy_0_0_0_z_z_yz_zz = buffer_2000_ppdd[1289];

    auto g_yy_0_0_0_z_z_zz_xx = buffer_2000_ppdd[1290];

    auto g_yy_0_0_0_z_z_zz_xy = buffer_2000_ppdd[1291];

    auto g_yy_0_0_0_z_z_zz_xz = buffer_2000_ppdd[1292];

    auto g_yy_0_0_0_z_z_zz_yy = buffer_2000_ppdd[1293];

    auto g_yy_0_0_0_z_z_zz_yz = buffer_2000_ppdd[1294];

    auto g_yy_0_0_0_z_z_zz_zz = buffer_2000_ppdd[1295];

    auto g_yz_0_0_0_x_x_xx_xx = buffer_2000_ppdd[1296];

    auto g_yz_0_0_0_x_x_xx_xy = buffer_2000_ppdd[1297];

    auto g_yz_0_0_0_x_x_xx_xz = buffer_2000_ppdd[1298];

    auto g_yz_0_0_0_x_x_xx_yy = buffer_2000_ppdd[1299];

    auto g_yz_0_0_0_x_x_xx_yz = buffer_2000_ppdd[1300];

    auto g_yz_0_0_0_x_x_xx_zz = buffer_2000_ppdd[1301];

    auto g_yz_0_0_0_x_x_xy_xx = buffer_2000_ppdd[1302];

    auto g_yz_0_0_0_x_x_xy_xy = buffer_2000_ppdd[1303];

    auto g_yz_0_0_0_x_x_xy_xz = buffer_2000_ppdd[1304];

    auto g_yz_0_0_0_x_x_xy_yy = buffer_2000_ppdd[1305];

    auto g_yz_0_0_0_x_x_xy_yz = buffer_2000_ppdd[1306];

    auto g_yz_0_0_0_x_x_xy_zz = buffer_2000_ppdd[1307];

    auto g_yz_0_0_0_x_x_xz_xx = buffer_2000_ppdd[1308];

    auto g_yz_0_0_0_x_x_xz_xy = buffer_2000_ppdd[1309];

    auto g_yz_0_0_0_x_x_xz_xz = buffer_2000_ppdd[1310];

    auto g_yz_0_0_0_x_x_xz_yy = buffer_2000_ppdd[1311];

    auto g_yz_0_0_0_x_x_xz_yz = buffer_2000_ppdd[1312];

    auto g_yz_0_0_0_x_x_xz_zz = buffer_2000_ppdd[1313];

    auto g_yz_0_0_0_x_x_yy_xx = buffer_2000_ppdd[1314];

    auto g_yz_0_0_0_x_x_yy_xy = buffer_2000_ppdd[1315];

    auto g_yz_0_0_0_x_x_yy_xz = buffer_2000_ppdd[1316];

    auto g_yz_0_0_0_x_x_yy_yy = buffer_2000_ppdd[1317];

    auto g_yz_0_0_0_x_x_yy_yz = buffer_2000_ppdd[1318];

    auto g_yz_0_0_0_x_x_yy_zz = buffer_2000_ppdd[1319];

    auto g_yz_0_0_0_x_x_yz_xx = buffer_2000_ppdd[1320];

    auto g_yz_0_0_0_x_x_yz_xy = buffer_2000_ppdd[1321];

    auto g_yz_0_0_0_x_x_yz_xz = buffer_2000_ppdd[1322];

    auto g_yz_0_0_0_x_x_yz_yy = buffer_2000_ppdd[1323];

    auto g_yz_0_0_0_x_x_yz_yz = buffer_2000_ppdd[1324];

    auto g_yz_0_0_0_x_x_yz_zz = buffer_2000_ppdd[1325];

    auto g_yz_0_0_0_x_x_zz_xx = buffer_2000_ppdd[1326];

    auto g_yz_0_0_0_x_x_zz_xy = buffer_2000_ppdd[1327];

    auto g_yz_0_0_0_x_x_zz_xz = buffer_2000_ppdd[1328];

    auto g_yz_0_0_0_x_x_zz_yy = buffer_2000_ppdd[1329];

    auto g_yz_0_0_0_x_x_zz_yz = buffer_2000_ppdd[1330];

    auto g_yz_0_0_0_x_x_zz_zz = buffer_2000_ppdd[1331];

    auto g_yz_0_0_0_x_y_xx_xx = buffer_2000_ppdd[1332];

    auto g_yz_0_0_0_x_y_xx_xy = buffer_2000_ppdd[1333];

    auto g_yz_0_0_0_x_y_xx_xz = buffer_2000_ppdd[1334];

    auto g_yz_0_0_0_x_y_xx_yy = buffer_2000_ppdd[1335];

    auto g_yz_0_0_0_x_y_xx_yz = buffer_2000_ppdd[1336];

    auto g_yz_0_0_0_x_y_xx_zz = buffer_2000_ppdd[1337];

    auto g_yz_0_0_0_x_y_xy_xx = buffer_2000_ppdd[1338];

    auto g_yz_0_0_0_x_y_xy_xy = buffer_2000_ppdd[1339];

    auto g_yz_0_0_0_x_y_xy_xz = buffer_2000_ppdd[1340];

    auto g_yz_0_0_0_x_y_xy_yy = buffer_2000_ppdd[1341];

    auto g_yz_0_0_0_x_y_xy_yz = buffer_2000_ppdd[1342];

    auto g_yz_0_0_0_x_y_xy_zz = buffer_2000_ppdd[1343];

    auto g_yz_0_0_0_x_y_xz_xx = buffer_2000_ppdd[1344];

    auto g_yz_0_0_0_x_y_xz_xy = buffer_2000_ppdd[1345];

    auto g_yz_0_0_0_x_y_xz_xz = buffer_2000_ppdd[1346];

    auto g_yz_0_0_0_x_y_xz_yy = buffer_2000_ppdd[1347];

    auto g_yz_0_0_0_x_y_xz_yz = buffer_2000_ppdd[1348];

    auto g_yz_0_0_0_x_y_xz_zz = buffer_2000_ppdd[1349];

    auto g_yz_0_0_0_x_y_yy_xx = buffer_2000_ppdd[1350];

    auto g_yz_0_0_0_x_y_yy_xy = buffer_2000_ppdd[1351];

    auto g_yz_0_0_0_x_y_yy_xz = buffer_2000_ppdd[1352];

    auto g_yz_0_0_0_x_y_yy_yy = buffer_2000_ppdd[1353];

    auto g_yz_0_0_0_x_y_yy_yz = buffer_2000_ppdd[1354];

    auto g_yz_0_0_0_x_y_yy_zz = buffer_2000_ppdd[1355];

    auto g_yz_0_0_0_x_y_yz_xx = buffer_2000_ppdd[1356];

    auto g_yz_0_0_0_x_y_yz_xy = buffer_2000_ppdd[1357];

    auto g_yz_0_0_0_x_y_yz_xz = buffer_2000_ppdd[1358];

    auto g_yz_0_0_0_x_y_yz_yy = buffer_2000_ppdd[1359];

    auto g_yz_0_0_0_x_y_yz_yz = buffer_2000_ppdd[1360];

    auto g_yz_0_0_0_x_y_yz_zz = buffer_2000_ppdd[1361];

    auto g_yz_0_0_0_x_y_zz_xx = buffer_2000_ppdd[1362];

    auto g_yz_0_0_0_x_y_zz_xy = buffer_2000_ppdd[1363];

    auto g_yz_0_0_0_x_y_zz_xz = buffer_2000_ppdd[1364];

    auto g_yz_0_0_0_x_y_zz_yy = buffer_2000_ppdd[1365];

    auto g_yz_0_0_0_x_y_zz_yz = buffer_2000_ppdd[1366];

    auto g_yz_0_0_0_x_y_zz_zz = buffer_2000_ppdd[1367];

    auto g_yz_0_0_0_x_z_xx_xx = buffer_2000_ppdd[1368];

    auto g_yz_0_0_0_x_z_xx_xy = buffer_2000_ppdd[1369];

    auto g_yz_0_0_0_x_z_xx_xz = buffer_2000_ppdd[1370];

    auto g_yz_0_0_0_x_z_xx_yy = buffer_2000_ppdd[1371];

    auto g_yz_0_0_0_x_z_xx_yz = buffer_2000_ppdd[1372];

    auto g_yz_0_0_0_x_z_xx_zz = buffer_2000_ppdd[1373];

    auto g_yz_0_0_0_x_z_xy_xx = buffer_2000_ppdd[1374];

    auto g_yz_0_0_0_x_z_xy_xy = buffer_2000_ppdd[1375];

    auto g_yz_0_0_0_x_z_xy_xz = buffer_2000_ppdd[1376];

    auto g_yz_0_0_0_x_z_xy_yy = buffer_2000_ppdd[1377];

    auto g_yz_0_0_0_x_z_xy_yz = buffer_2000_ppdd[1378];

    auto g_yz_0_0_0_x_z_xy_zz = buffer_2000_ppdd[1379];

    auto g_yz_0_0_0_x_z_xz_xx = buffer_2000_ppdd[1380];

    auto g_yz_0_0_0_x_z_xz_xy = buffer_2000_ppdd[1381];

    auto g_yz_0_0_0_x_z_xz_xz = buffer_2000_ppdd[1382];

    auto g_yz_0_0_0_x_z_xz_yy = buffer_2000_ppdd[1383];

    auto g_yz_0_0_0_x_z_xz_yz = buffer_2000_ppdd[1384];

    auto g_yz_0_0_0_x_z_xz_zz = buffer_2000_ppdd[1385];

    auto g_yz_0_0_0_x_z_yy_xx = buffer_2000_ppdd[1386];

    auto g_yz_0_0_0_x_z_yy_xy = buffer_2000_ppdd[1387];

    auto g_yz_0_0_0_x_z_yy_xz = buffer_2000_ppdd[1388];

    auto g_yz_0_0_0_x_z_yy_yy = buffer_2000_ppdd[1389];

    auto g_yz_0_0_0_x_z_yy_yz = buffer_2000_ppdd[1390];

    auto g_yz_0_0_0_x_z_yy_zz = buffer_2000_ppdd[1391];

    auto g_yz_0_0_0_x_z_yz_xx = buffer_2000_ppdd[1392];

    auto g_yz_0_0_0_x_z_yz_xy = buffer_2000_ppdd[1393];

    auto g_yz_0_0_0_x_z_yz_xz = buffer_2000_ppdd[1394];

    auto g_yz_0_0_0_x_z_yz_yy = buffer_2000_ppdd[1395];

    auto g_yz_0_0_0_x_z_yz_yz = buffer_2000_ppdd[1396];

    auto g_yz_0_0_0_x_z_yz_zz = buffer_2000_ppdd[1397];

    auto g_yz_0_0_0_x_z_zz_xx = buffer_2000_ppdd[1398];

    auto g_yz_0_0_0_x_z_zz_xy = buffer_2000_ppdd[1399];

    auto g_yz_0_0_0_x_z_zz_xz = buffer_2000_ppdd[1400];

    auto g_yz_0_0_0_x_z_zz_yy = buffer_2000_ppdd[1401];

    auto g_yz_0_0_0_x_z_zz_yz = buffer_2000_ppdd[1402];

    auto g_yz_0_0_0_x_z_zz_zz = buffer_2000_ppdd[1403];

    auto g_yz_0_0_0_y_x_xx_xx = buffer_2000_ppdd[1404];

    auto g_yz_0_0_0_y_x_xx_xy = buffer_2000_ppdd[1405];

    auto g_yz_0_0_0_y_x_xx_xz = buffer_2000_ppdd[1406];

    auto g_yz_0_0_0_y_x_xx_yy = buffer_2000_ppdd[1407];

    auto g_yz_0_0_0_y_x_xx_yz = buffer_2000_ppdd[1408];

    auto g_yz_0_0_0_y_x_xx_zz = buffer_2000_ppdd[1409];

    auto g_yz_0_0_0_y_x_xy_xx = buffer_2000_ppdd[1410];

    auto g_yz_0_0_0_y_x_xy_xy = buffer_2000_ppdd[1411];

    auto g_yz_0_0_0_y_x_xy_xz = buffer_2000_ppdd[1412];

    auto g_yz_0_0_0_y_x_xy_yy = buffer_2000_ppdd[1413];

    auto g_yz_0_0_0_y_x_xy_yz = buffer_2000_ppdd[1414];

    auto g_yz_0_0_0_y_x_xy_zz = buffer_2000_ppdd[1415];

    auto g_yz_0_0_0_y_x_xz_xx = buffer_2000_ppdd[1416];

    auto g_yz_0_0_0_y_x_xz_xy = buffer_2000_ppdd[1417];

    auto g_yz_0_0_0_y_x_xz_xz = buffer_2000_ppdd[1418];

    auto g_yz_0_0_0_y_x_xz_yy = buffer_2000_ppdd[1419];

    auto g_yz_0_0_0_y_x_xz_yz = buffer_2000_ppdd[1420];

    auto g_yz_0_0_0_y_x_xz_zz = buffer_2000_ppdd[1421];

    auto g_yz_0_0_0_y_x_yy_xx = buffer_2000_ppdd[1422];

    auto g_yz_0_0_0_y_x_yy_xy = buffer_2000_ppdd[1423];

    auto g_yz_0_0_0_y_x_yy_xz = buffer_2000_ppdd[1424];

    auto g_yz_0_0_0_y_x_yy_yy = buffer_2000_ppdd[1425];

    auto g_yz_0_0_0_y_x_yy_yz = buffer_2000_ppdd[1426];

    auto g_yz_0_0_0_y_x_yy_zz = buffer_2000_ppdd[1427];

    auto g_yz_0_0_0_y_x_yz_xx = buffer_2000_ppdd[1428];

    auto g_yz_0_0_0_y_x_yz_xy = buffer_2000_ppdd[1429];

    auto g_yz_0_0_0_y_x_yz_xz = buffer_2000_ppdd[1430];

    auto g_yz_0_0_0_y_x_yz_yy = buffer_2000_ppdd[1431];

    auto g_yz_0_0_0_y_x_yz_yz = buffer_2000_ppdd[1432];

    auto g_yz_0_0_0_y_x_yz_zz = buffer_2000_ppdd[1433];

    auto g_yz_0_0_0_y_x_zz_xx = buffer_2000_ppdd[1434];

    auto g_yz_0_0_0_y_x_zz_xy = buffer_2000_ppdd[1435];

    auto g_yz_0_0_0_y_x_zz_xz = buffer_2000_ppdd[1436];

    auto g_yz_0_0_0_y_x_zz_yy = buffer_2000_ppdd[1437];

    auto g_yz_0_0_0_y_x_zz_yz = buffer_2000_ppdd[1438];

    auto g_yz_0_0_0_y_x_zz_zz = buffer_2000_ppdd[1439];

    auto g_yz_0_0_0_y_y_xx_xx = buffer_2000_ppdd[1440];

    auto g_yz_0_0_0_y_y_xx_xy = buffer_2000_ppdd[1441];

    auto g_yz_0_0_0_y_y_xx_xz = buffer_2000_ppdd[1442];

    auto g_yz_0_0_0_y_y_xx_yy = buffer_2000_ppdd[1443];

    auto g_yz_0_0_0_y_y_xx_yz = buffer_2000_ppdd[1444];

    auto g_yz_0_0_0_y_y_xx_zz = buffer_2000_ppdd[1445];

    auto g_yz_0_0_0_y_y_xy_xx = buffer_2000_ppdd[1446];

    auto g_yz_0_0_0_y_y_xy_xy = buffer_2000_ppdd[1447];

    auto g_yz_0_0_0_y_y_xy_xz = buffer_2000_ppdd[1448];

    auto g_yz_0_0_0_y_y_xy_yy = buffer_2000_ppdd[1449];

    auto g_yz_0_0_0_y_y_xy_yz = buffer_2000_ppdd[1450];

    auto g_yz_0_0_0_y_y_xy_zz = buffer_2000_ppdd[1451];

    auto g_yz_0_0_0_y_y_xz_xx = buffer_2000_ppdd[1452];

    auto g_yz_0_0_0_y_y_xz_xy = buffer_2000_ppdd[1453];

    auto g_yz_0_0_0_y_y_xz_xz = buffer_2000_ppdd[1454];

    auto g_yz_0_0_0_y_y_xz_yy = buffer_2000_ppdd[1455];

    auto g_yz_0_0_0_y_y_xz_yz = buffer_2000_ppdd[1456];

    auto g_yz_0_0_0_y_y_xz_zz = buffer_2000_ppdd[1457];

    auto g_yz_0_0_0_y_y_yy_xx = buffer_2000_ppdd[1458];

    auto g_yz_0_0_0_y_y_yy_xy = buffer_2000_ppdd[1459];

    auto g_yz_0_0_0_y_y_yy_xz = buffer_2000_ppdd[1460];

    auto g_yz_0_0_0_y_y_yy_yy = buffer_2000_ppdd[1461];

    auto g_yz_0_0_0_y_y_yy_yz = buffer_2000_ppdd[1462];

    auto g_yz_0_0_0_y_y_yy_zz = buffer_2000_ppdd[1463];

    auto g_yz_0_0_0_y_y_yz_xx = buffer_2000_ppdd[1464];

    auto g_yz_0_0_0_y_y_yz_xy = buffer_2000_ppdd[1465];

    auto g_yz_0_0_0_y_y_yz_xz = buffer_2000_ppdd[1466];

    auto g_yz_0_0_0_y_y_yz_yy = buffer_2000_ppdd[1467];

    auto g_yz_0_0_0_y_y_yz_yz = buffer_2000_ppdd[1468];

    auto g_yz_0_0_0_y_y_yz_zz = buffer_2000_ppdd[1469];

    auto g_yz_0_0_0_y_y_zz_xx = buffer_2000_ppdd[1470];

    auto g_yz_0_0_0_y_y_zz_xy = buffer_2000_ppdd[1471];

    auto g_yz_0_0_0_y_y_zz_xz = buffer_2000_ppdd[1472];

    auto g_yz_0_0_0_y_y_zz_yy = buffer_2000_ppdd[1473];

    auto g_yz_0_0_0_y_y_zz_yz = buffer_2000_ppdd[1474];

    auto g_yz_0_0_0_y_y_zz_zz = buffer_2000_ppdd[1475];

    auto g_yz_0_0_0_y_z_xx_xx = buffer_2000_ppdd[1476];

    auto g_yz_0_0_0_y_z_xx_xy = buffer_2000_ppdd[1477];

    auto g_yz_0_0_0_y_z_xx_xz = buffer_2000_ppdd[1478];

    auto g_yz_0_0_0_y_z_xx_yy = buffer_2000_ppdd[1479];

    auto g_yz_0_0_0_y_z_xx_yz = buffer_2000_ppdd[1480];

    auto g_yz_0_0_0_y_z_xx_zz = buffer_2000_ppdd[1481];

    auto g_yz_0_0_0_y_z_xy_xx = buffer_2000_ppdd[1482];

    auto g_yz_0_0_0_y_z_xy_xy = buffer_2000_ppdd[1483];

    auto g_yz_0_0_0_y_z_xy_xz = buffer_2000_ppdd[1484];

    auto g_yz_0_0_0_y_z_xy_yy = buffer_2000_ppdd[1485];

    auto g_yz_0_0_0_y_z_xy_yz = buffer_2000_ppdd[1486];

    auto g_yz_0_0_0_y_z_xy_zz = buffer_2000_ppdd[1487];

    auto g_yz_0_0_0_y_z_xz_xx = buffer_2000_ppdd[1488];

    auto g_yz_0_0_0_y_z_xz_xy = buffer_2000_ppdd[1489];

    auto g_yz_0_0_0_y_z_xz_xz = buffer_2000_ppdd[1490];

    auto g_yz_0_0_0_y_z_xz_yy = buffer_2000_ppdd[1491];

    auto g_yz_0_0_0_y_z_xz_yz = buffer_2000_ppdd[1492];

    auto g_yz_0_0_0_y_z_xz_zz = buffer_2000_ppdd[1493];

    auto g_yz_0_0_0_y_z_yy_xx = buffer_2000_ppdd[1494];

    auto g_yz_0_0_0_y_z_yy_xy = buffer_2000_ppdd[1495];

    auto g_yz_0_0_0_y_z_yy_xz = buffer_2000_ppdd[1496];

    auto g_yz_0_0_0_y_z_yy_yy = buffer_2000_ppdd[1497];

    auto g_yz_0_0_0_y_z_yy_yz = buffer_2000_ppdd[1498];

    auto g_yz_0_0_0_y_z_yy_zz = buffer_2000_ppdd[1499];

    auto g_yz_0_0_0_y_z_yz_xx = buffer_2000_ppdd[1500];

    auto g_yz_0_0_0_y_z_yz_xy = buffer_2000_ppdd[1501];

    auto g_yz_0_0_0_y_z_yz_xz = buffer_2000_ppdd[1502];

    auto g_yz_0_0_0_y_z_yz_yy = buffer_2000_ppdd[1503];

    auto g_yz_0_0_0_y_z_yz_yz = buffer_2000_ppdd[1504];

    auto g_yz_0_0_0_y_z_yz_zz = buffer_2000_ppdd[1505];

    auto g_yz_0_0_0_y_z_zz_xx = buffer_2000_ppdd[1506];

    auto g_yz_0_0_0_y_z_zz_xy = buffer_2000_ppdd[1507];

    auto g_yz_0_0_0_y_z_zz_xz = buffer_2000_ppdd[1508];

    auto g_yz_0_0_0_y_z_zz_yy = buffer_2000_ppdd[1509];

    auto g_yz_0_0_0_y_z_zz_yz = buffer_2000_ppdd[1510];

    auto g_yz_0_0_0_y_z_zz_zz = buffer_2000_ppdd[1511];

    auto g_yz_0_0_0_z_x_xx_xx = buffer_2000_ppdd[1512];

    auto g_yz_0_0_0_z_x_xx_xy = buffer_2000_ppdd[1513];

    auto g_yz_0_0_0_z_x_xx_xz = buffer_2000_ppdd[1514];

    auto g_yz_0_0_0_z_x_xx_yy = buffer_2000_ppdd[1515];

    auto g_yz_0_0_0_z_x_xx_yz = buffer_2000_ppdd[1516];

    auto g_yz_0_0_0_z_x_xx_zz = buffer_2000_ppdd[1517];

    auto g_yz_0_0_0_z_x_xy_xx = buffer_2000_ppdd[1518];

    auto g_yz_0_0_0_z_x_xy_xy = buffer_2000_ppdd[1519];

    auto g_yz_0_0_0_z_x_xy_xz = buffer_2000_ppdd[1520];

    auto g_yz_0_0_0_z_x_xy_yy = buffer_2000_ppdd[1521];

    auto g_yz_0_0_0_z_x_xy_yz = buffer_2000_ppdd[1522];

    auto g_yz_0_0_0_z_x_xy_zz = buffer_2000_ppdd[1523];

    auto g_yz_0_0_0_z_x_xz_xx = buffer_2000_ppdd[1524];

    auto g_yz_0_0_0_z_x_xz_xy = buffer_2000_ppdd[1525];

    auto g_yz_0_0_0_z_x_xz_xz = buffer_2000_ppdd[1526];

    auto g_yz_0_0_0_z_x_xz_yy = buffer_2000_ppdd[1527];

    auto g_yz_0_0_0_z_x_xz_yz = buffer_2000_ppdd[1528];

    auto g_yz_0_0_0_z_x_xz_zz = buffer_2000_ppdd[1529];

    auto g_yz_0_0_0_z_x_yy_xx = buffer_2000_ppdd[1530];

    auto g_yz_0_0_0_z_x_yy_xy = buffer_2000_ppdd[1531];

    auto g_yz_0_0_0_z_x_yy_xz = buffer_2000_ppdd[1532];

    auto g_yz_0_0_0_z_x_yy_yy = buffer_2000_ppdd[1533];

    auto g_yz_0_0_0_z_x_yy_yz = buffer_2000_ppdd[1534];

    auto g_yz_0_0_0_z_x_yy_zz = buffer_2000_ppdd[1535];

    auto g_yz_0_0_0_z_x_yz_xx = buffer_2000_ppdd[1536];

    auto g_yz_0_0_0_z_x_yz_xy = buffer_2000_ppdd[1537];

    auto g_yz_0_0_0_z_x_yz_xz = buffer_2000_ppdd[1538];

    auto g_yz_0_0_0_z_x_yz_yy = buffer_2000_ppdd[1539];

    auto g_yz_0_0_0_z_x_yz_yz = buffer_2000_ppdd[1540];

    auto g_yz_0_0_0_z_x_yz_zz = buffer_2000_ppdd[1541];

    auto g_yz_0_0_0_z_x_zz_xx = buffer_2000_ppdd[1542];

    auto g_yz_0_0_0_z_x_zz_xy = buffer_2000_ppdd[1543];

    auto g_yz_0_0_0_z_x_zz_xz = buffer_2000_ppdd[1544];

    auto g_yz_0_0_0_z_x_zz_yy = buffer_2000_ppdd[1545];

    auto g_yz_0_0_0_z_x_zz_yz = buffer_2000_ppdd[1546];

    auto g_yz_0_0_0_z_x_zz_zz = buffer_2000_ppdd[1547];

    auto g_yz_0_0_0_z_y_xx_xx = buffer_2000_ppdd[1548];

    auto g_yz_0_0_0_z_y_xx_xy = buffer_2000_ppdd[1549];

    auto g_yz_0_0_0_z_y_xx_xz = buffer_2000_ppdd[1550];

    auto g_yz_0_0_0_z_y_xx_yy = buffer_2000_ppdd[1551];

    auto g_yz_0_0_0_z_y_xx_yz = buffer_2000_ppdd[1552];

    auto g_yz_0_0_0_z_y_xx_zz = buffer_2000_ppdd[1553];

    auto g_yz_0_0_0_z_y_xy_xx = buffer_2000_ppdd[1554];

    auto g_yz_0_0_0_z_y_xy_xy = buffer_2000_ppdd[1555];

    auto g_yz_0_0_0_z_y_xy_xz = buffer_2000_ppdd[1556];

    auto g_yz_0_0_0_z_y_xy_yy = buffer_2000_ppdd[1557];

    auto g_yz_0_0_0_z_y_xy_yz = buffer_2000_ppdd[1558];

    auto g_yz_0_0_0_z_y_xy_zz = buffer_2000_ppdd[1559];

    auto g_yz_0_0_0_z_y_xz_xx = buffer_2000_ppdd[1560];

    auto g_yz_0_0_0_z_y_xz_xy = buffer_2000_ppdd[1561];

    auto g_yz_0_0_0_z_y_xz_xz = buffer_2000_ppdd[1562];

    auto g_yz_0_0_0_z_y_xz_yy = buffer_2000_ppdd[1563];

    auto g_yz_0_0_0_z_y_xz_yz = buffer_2000_ppdd[1564];

    auto g_yz_0_0_0_z_y_xz_zz = buffer_2000_ppdd[1565];

    auto g_yz_0_0_0_z_y_yy_xx = buffer_2000_ppdd[1566];

    auto g_yz_0_0_0_z_y_yy_xy = buffer_2000_ppdd[1567];

    auto g_yz_0_0_0_z_y_yy_xz = buffer_2000_ppdd[1568];

    auto g_yz_0_0_0_z_y_yy_yy = buffer_2000_ppdd[1569];

    auto g_yz_0_0_0_z_y_yy_yz = buffer_2000_ppdd[1570];

    auto g_yz_0_0_0_z_y_yy_zz = buffer_2000_ppdd[1571];

    auto g_yz_0_0_0_z_y_yz_xx = buffer_2000_ppdd[1572];

    auto g_yz_0_0_0_z_y_yz_xy = buffer_2000_ppdd[1573];

    auto g_yz_0_0_0_z_y_yz_xz = buffer_2000_ppdd[1574];

    auto g_yz_0_0_0_z_y_yz_yy = buffer_2000_ppdd[1575];

    auto g_yz_0_0_0_z_y_yz_yz = buffer_2000_ppdd[1576];

    auto g_yz_0_0_0_z_y_yz_zz = buffer_2000_ppdd[1577];

    auto g_yz_0_0_0_z_y_zz_xx = buffer_2000_ppdd[1578];

    auto g_yz_0_0_0_z_y_zz_xy = buffer_2000_ppdd[1579];

    auto g_yz_0_0_0_z_y_zz_xz = buffer_2000_ppdd[1580];

    auto g_yz_0_0_0_z_y_zz_yy = buffer_2000_ppdd[1581];

    auto g_yz_0_0_0_z_y_zz_yz = buffer_2000_ppdd[1582];

    auto g_yz_0_0_0_z_y_zz_zz = buffer_2000_ppdd[1583];

    auto g_yz_0_0_0_z_z_xx_xx = buffer_2000_ppdd[1584];

    auto g_yz_0_0_0_z_z_xx_xy = buffer_2000_ppdd[1585];

    auto g_yz_0_0_0_z_z_xx_xz = buffer_2000_ppdd[1586];

    auto g_yz_0_0_0_z_z_xx_yy = buffer_2000_ppdd[1587];

    auto g_yz_0_0_0_z_z_xx_yz = buffer_2000_ppdd[1588];

    auto g_yz_0_0_0_z_z_xx_zz = buffer_2000_ppdd[1589];

    auto g_yz_0_0_0_z_z_xy_xx = buffer_2000_ppdd[1590];

    auto g_yz_0_0_0_z_z_xy_xy = buffer_2000_ppdd[1591];

    auto g_yz_0_0_0_z_z_xy_xz = buffer_2000_ppdd[1592];

    auto g_yz_0_0_0_z_z_xy_yy = buffer_2000_ppdd[1593];

    auto g_yz_0_0_0_z_z_xy_yz = buffer_2000_ppdd[1594];

    auto g_yz_0_0_0_z_z_xy_zz = buffer_2000_ppdd[1595];

    auto g_yz_0_0_0_z_z_xz_xx = buffer_2000_ppdd[1596];

    auto g_yz_0_0_0_z_z_xz_xy = buffer_2000_ppdd[1597];

    auto g_yz_0_0_0_z_z_xz_xz = buffer_2000_ppdd[1598];

    auto g_yz_0_0_0_z_z_xz_yy = buffer_2000_ppdd[1599];

    auto g_yz_0_0_0_z_z_xz_yz = buffer_2000_ppdd[1600];

    auto g_yz_0_0_0_z_z_xz_zz = buffer_2000_ppdd[1601];

    auto g_yz_0_0_0_z_z_yy_xx = buffer_2000_ppdd[1602];

    auto g_yz_0_0_0_z_z_yy_xy = buffer_2000_ppdd[1603];

    auto g_yz_0_0_0_z_z_yy_xz = buffer_2000_ppdd[1604];

    auto g_yz_0_0_0_z_z_yy_yy = buffer_2000_ppdd[1605];

    auto g_yz_0_0_0_z_z_yy_yz = buffer_2000_ppdd[1606];

    auto g_yz_0_0_0_z_z_yy_zz = buffer_2000_ppdd[1607];

    auto g_yz_0_0_0_z_z_yz_xx = buffer_2000_ppdd[1608];

    auto g_yz_0_0_0_z_z_yz_xy = buffer_2000_ppdd[1609];

    auto g_yz_0_0_0_z_z_yz_xz = buffer_2000_ppdd[1610];

    auto g_yz_0_0_0_z_z_yz_yy = buffer_2000_ppdd[1611];

    auto g_yz_0_0_0_z_z_yz_yz = buffer_2000_ppdd[1612];

    auto g_yz_0_0_0_z_z_yz_zz = buffer_2000_ppdd[1613];

    auto g_yz_0_0_0_z_z_zz_xx = buffer_2000_ppdd[1614];

    auto g_yz_0_0_0_z_z_zz_xy = buffer_2000_ppdd[1615];

    auto g_yz_0_0_0_z_z_zz_xz = buffer_2000_ppdd[1616];

    auto g_yz_0_0_0_z_z_zz_yy = buffer_2000_ppdd[1617];

    auto g_yz_0_0_0_z_z_zz_yz = buffer_2000_ppdd[1618];

    auto g_yz_0_0_0_z_z_zz_zz = buffer_2000_ppdd[1619];

    auto g_zz_0_0_0_x_x_xx_xx = buffer_2000_ppdd[1620];

    auto g_zz_0_0_0_x_x_xx_xy = buffer_2000_ppdd[1621];

    auto g_zz_0_0_0_x_x_xx_xz = buffer_2000_ppdd[1622];

    auto g_zz_0_0_0_x_x_xx_yy = buffer_2000_ppdd[1623];

    auto g_zz_0_0_0_x_x_xx_yz = buffer_2000_ppdd[1624];

    auto g_zz_0_0_0_x_x_xx_zz = buffer_2000_ppdd[1625];

    auto g_zz_0_0_0_x_x_xy_xx = buffer_2000_ppdd[1626];

    auto g_zz_0_0_0_x_x_xy_xy = buffer_2000_ppdd[1627];

    auto g_zz_0_0_0_x_x_xy_xz = buffer_2000_ppdd[1628];

    auto g_zz_0_0_0_x_x_xy_yy = buffer_2000_ppdd[1629];

    auto g_zz_0_0_0_x_x_xy_yz = buffer_2000_ppdd[1630];

    auto g_zz_0_0_0_x_x_xy_zz = buffer_2000_ppdd[1631];

    auto g_zz_0_0_0_x_x_xz_xx = buffer_2000_ppdd[1632];

    auto g_zz_0_0_0_x_x_xz_xy = buffer_2000_ppdd[1633];

    auto g_zz_0_0_0_x_x_xz_xz = buffer_2000_ppdd[1634];

    auto g_zz_0_0_0_x_x_xz_yy = buffer_2000_ppdd[1635];

    auto g_zz_0_0_0_x_x_xz_yz = buffer_2000_ppdd[1636];

    auto g_zz_0_0_0_x_x_xz_zz = buffer_2000_ppdd[1637];

    auto g_zz_0_0_0_x_x_yy_xx = buffer_2000_ppdd[1638];

    auto g_zz_0_0_0_x_x_yy_xy = buffer_2000_ppdd[1639];

    auto g_zz_0_0_0_x_x_yy_xz = buffer_2000_ppdd[1640];

    auto g_zz_0_0_0_x_x_yy_yy = buffer_2000_ppdd[1641];

    auto g_zz_0_0_0_x_x_yy_yz = buffer_2000_ppdd[1642];

    auto g_zz_0_0_0_x_x_yy_zz = buffer_2000_ppdd[1643];

    auto g_zz_0_0_0_x_x_yz_xx = buffer_2000_ppdd[1644];

    auto g_zz_0_0_0_x_x_yz_xy = buffer_2000_ppdd[1645];

    auto g_zz_0_0_0_x_x_yz_xz = buffer_2000_ppdd[1646];

    auto g_zz_0_0_0_x_x_yz_yy = buffer_2000_ppdd[1647];

    auto g_zz_0_0_0_x_x_yz_yz = buffer_2000_ppdd[1648];

    auto g_zz_0_0_0_x_x_yz_zz = buffer_2000_ppdd[1649];

    auto g_zz_0_0_0_x_x_zz_xx = buffer_2000_ppdd[1650];

    auto g_zz_0_0_0_x_x_zz_xy = buffer_2000_ppdd[1651];

    auto g_zz_0_0_0_x_x_zz_xz = buffer_2000_ppdd[1652];

    auto g_zz_0_0_0_x_x_zz_yy = buffer_2000_ppdd[1653];

    auto g_zz_0_0_0_x_x_zz_yz = buffer_2000_ppdd[1654];

    auto g_zz_0_0_0_x_x_zz_zz = buffer_2000_ppdd[1655];

    auto g_zz_0_0_0_x_y_xx_xx = buffer_2000_ppdd[1656];

    auto g_zz_0_0_0_x_y_xx_xy = buffer_2000_ppdd[1657];

    auto g_zz_0_0_0_x_y_xx_xz = buffer_2000_ppdd[1658];

    auto g_zz_0_0_0_x_y_xx_yy = buffer_2000_ppdd[1659];

    auto g_zz_0_0_0_x_y_xx_yz = buffer_2000_ppdd[1660];

    auto g_zz_0_0_0_x_y_xx_zz = buffer_2000_ppdd[1661];

    auto g_zz_0_0_0_x_y_xy_xx = buffer_2000_ppdd[1662];

    auto g_zz_0_0_0_x_y_xy_xy = buffer_2000_ppdd[1663];

    auto g_zz_0_0_0_x_y_xy_xz = buffer_2000_ppdd[1664];

    auto g_zz_0_0_0_x_y_xy_yy = buffer_2000_ppdd[1665];

    auto g_zz_0_0_0_x_y_xy_yz = buffer_2000_ppdd[1666];

    auto g_zz_0_0_0_x_y_xy_zz = buffer_2000_ppdd[1667];

    auto g_zz_0_0_0_x_y_xz_xx = buffer_2000_ppdd[1668];

    auto g_zz_0_0_0_x_y_xz_xy = buffer_2000_ppdd[1669];

    auto g_zz_0_0_0_x_y_xz_xz = buffer_2000_ppdd[1670];

    auto g_zz_0_0_0_x_y_xz_yy = buffer_2000_ppdd[1671];

    auto g_zz_0_0_0_x_y_xz_yz = buffer_2000_ppdd[1672];

    auto g_zz_0_0_0_x_y_xz_zz = buffer_2000_ppdd[1673];

    auto g_zz_0_0_0_x_y_yy_xx = buffer_2000_ppdd[1674];

    auto g_zz_0_0_0_x_y_yy_xy = buffer_2000_ppdd[1675];

    auto g_zz_0_0_0_x_y_yy_xz = buffer_2000_ppdd[1676];

    auto g_zz_0_0_0_x_y_yy_yy = buffer_2000_ppdd[1677];

    auto g_zz_0_0_0_x_y_yy_yz = buffer_2000_ppdd[1678];

    auto g_zz_0_0_0_x_y_yy_zz = buffer_2000_ppdd[1679];

    auto g_zz_0_0_0_x_y_yz_xx = buffer_2000_ppdd[1680];

    auto g_zz_0_0_0_x_y_yz_xy = buffer_2000_ppdd[1681];

    auto g_zz_0_0_0_x_y_yz_xz = buffer_2000_ppdd[1682];

    auto g_zz_0_0_0_x_y_yz_yy = buffer_2000_ppdd[1683];

    auto g_zz_0_0_0_x_y_yz_yz = buffer_2000_ppdd[1684];

    auto g_zz_0_0_0_x_y_yz_zz = buffer_2000_ppdd[1685];

    auto g_zz_0_0_0_x_y_zz_xx = buffer_2000_ppdd[1686];

    auto g_zz_0_0_0_x_y_zz_xy = buffer_2000_ppdd[1687];

    auto g_zz_0_0_0_x_y_zz_xz = buffer_2000_ppdd[1688];

    auto g_zz_0_0_0_x_y_zz_yy = buffer_2000_ppdd[1689];

    auto g_zz_0_0_0_x_y_zz_yz = buffer_2000_ppdd[1690];

    auto g_zz_0_0_0_x_y_zz_zz = buffer_2000_ppdd[1691];

    auto g_zz_0_0_0_x_z_xx_xx = buffer_2000_ppdd[1692];

    auto g_zz_0_0_0_x_z_xx_xy = buffer_2000_ppdd[1693];

    auto g_zz_0_0_0_x_z_xx_xz = buffer_2000_ppdd[1694];

    auto g_zz_0_0_0_x_z_xx_yy = buffer_2000_ppdd[1695];

    auto g_zz_0_0_0_x_z_xx_yz = buffer_2000_ppdd[1696];

    auto g_zz_0_0_0_x_z_xx_zz = buffer_2000_ppdd[1697];

    auto g_zz_0_0_0_x_z_xy_xx = buffer_2000_ppdd[1698];

    auto g_zz_0_0_0_x_z_xy_xy = buffer_2000_ppdd[1699];

    auto g_zz_0_0_0_x_z_xy_xz = buffer_2000_ppdd[1700];

    auto g_zz_0_0_0_x_z_xy_yy = buffer_2000_ppdd[1701];

    auto g_zz_0_0_0_x_z_xy_yz = buffer_2000_ppdd[1702];

    auto g_zz_0_0_0_x_z_xy_zz = buffer_2000_ppdd[1703];

    auto g_zz_0_0_0_x_z_xz_xx = buffer_2000_ppdd[1704];

    auto g_zz_0_0_0_x_z_xz_xy = buffer_2000_ppdd[1705];

    auto g_zz_0_0_0_x_z_xz_xz = buffer_2000_ppdd[1706];

    auto g_zz_0_0_0_x_z_xz_yy = buffer_2000_ppdd[1707];

    auto g_zz_0_0_0_x_z_xz_yz = buffer_2000_ppdd[1708];

    auto g_zz_0_0_0_x_z_xz_zz = buffer_2000_ppdd[1709];

    auto g_zz_0_0_0_x_z_yy_xx = buffer_2000_ppdd[1710];

    auto g_zz_0_0_0_x_z_yy_xy = buffer_2000_ppdd[1711];

    auto g_zz_0_0_0_x_z_yy_xz = buffer_2000_ppdd[1712];

    auto g_zz_0_0_0_x_z_yy_yy = buffer_2000_ppdd[1713];

    auto g_zz_0_0_0_x_z_yy_yz = buffer_2000_ppdd[1714];

    auto g_zz_0_0_0_x_z_yy_zz = buffer_2000_ppdd[1715];

    auto g_zz_0_0_0_x_z_yz_xx = buffer_2000_ppdd[1716];

    auto g_zz_0_0_0_x_z_yz_xy = buffer_2000_ppdd[1717];

    auto g_zz_0_0_0_x_z_yz_xz = buffer_2000_ppdd[1718];

    auto g_zz_0_0_0_x_z_yz_yy = buffer_2000_ppdd[1719];

    auto g_zz_0_0_0_x_z_yz_yz = buffer_2000_ppdd[1720];

    auto g_zz_0_0_0_x_z_yz_zz = buffer_2000_ppdd[1721];

    auto g_zz_0_0_0_x_z_zz_xx = buffer_2000_ppdd[1722];

    auto g_zz_0_0_0_x_z_zz_xy = buffer_2000_ppdd[1723];

    auto g_zz_0_0_0_x_z_zz_xz = buffer_2000_ppdd[1724];

    auto g_zz_0_0_0_x_z_zz_yy = buffer_2000_ppdd[1725];

    auto g_zz_0_0_0_x_z_zz_yz = buffer_2000_ppdd[1726];

    auto g_zz_0_0_0_x_z_zz_zz = buffer_2000_ppdd[1727];

    auto g_zz_0_0_0_y_x_xx_xx = buffer_2000_ppdd[1728];

    auto g_zz_0_0_0_y_x_xx_xy = buffer_2000_ppdd[1729];

    auto g_zz_0_0_0_y_x_xx_xz = buffer_2000_ppdd[1730];

    auto g_zz_0_0_0_y_x_xx_yy = buffer_2000_ppdd[1731];

    auto g_zz_0_0_0_y_x_xx_yz = buffer_2000_ppdd[1732];

    auto g_zz_0_0_0_y_x_xx_zz = buffer_2000_ppdd[1733];

    auto g_zz_0_0_0_y_x_xy_xx = buffer_2000_ppdd[1734];

    auto g_zz_0_0_0_y_x_xy_xy = buffer_2000_ppdd[1735];

    auto g_zz_0_0_0_y_x_xy_xz = buffer_2000_ppdd[1736];

    auto g_zz_0_0_0_y_x_xy_yy = buffer_2000_ppdd[1737];

    auto g_zz_0_0_0_y_x_xy_yz = buffer_2000_ppdd[1738];

    auto g_zz_0_0_0_y_x_xy_zz = buffer_2000_ppdd[1739];

    auto g_zz_0_0_0_y_x_xz_xx = buffer_2000_ppdd[1740];

    auto g_zz_0_0_0_y_x_xz_xy = buffer_2000_ppdd[1741];

    auto g_zz_0_0_0_y_x_xz_xz = buffer_2000_ppdd[1742];

    auto g_zz_0_0_0_y_x_xz_yy = buffer_2000_ppdd[1743];

    auto g_zz_0_0_0_y_x_xz_yz = buffer_2000_ppdd[1744];

    auto g_zz_0_0_0_y_x_xz_zz = buffer_2000_ppdd[1745];

    auto g_zz_0_0_0_y_x_yy_xx = buffer_2000_ppdd[1746];

    auto g_zz_0_0_0_y_x_yy_xy = buffer_2000_ppdd[1747];

    auto g_zz_0_0_0_y_x_yy_xz = buffer_2000_ppdd[1748];

    auto g_zz_0_0_0_y_x_yy_yy = buffer_2000_ppdd[1749];

    auto g_zz_0_0_0_y_x_yy_yz = buffer_2000_ppdd[1750];

    auto g_zz_0_0_0_y_x_yy_zz = buffer_2000_ppdd[1751];

    auto g_zz_0_0_0_y_x_yz_xx = buffer_2000_ppdd[1752];

    auto g_zz_0_0_0_y_x_yz_xy = buffer_2000_ppdd[1753];

    auto g_zz_0_0_0_y_x_yz_xz = buffer_2000_ppdd[1754];

    auto g_zz_0_0_0_y_x_yz_yy = buffer_2000_ppdd[1755];

    auto g_zz_0_0_0_y_x_yz_yz = buffer_2000_ppdd[1756];

    auto g_zz_0_0_0_y_x_yz_zz = buffer_2000_ppdd[1757];

    auto g_zz_0_0_0_y_x_zz_xx = buffer_2000_ppdd[1758];

    auto g_zz_0_0_0_y_x_zz_xy = buffer_2000_ppdd[1759];

    auto g_zz_0_0_0_y_x_zz_xz = buffer_2000_ppdd[1760];

    auto g_zz_0_0_0_y_x_zz_yy = buffer_2000_ppdd[1761];

    auto g_zz_0_0_0_y_x_zz_yz = buffer_2000_ppdd[1762];

    auto g_zz_0_0_0_y_x_zz_zz = buffer_2000_ppdd[1763];

    auto g_zz_0_0_0_y_y_xx_xx = buffer_2000_ppdd[1764];

    auto g_zz_0_0_0_y_y_xx_xy = buffer_2000_ppdd[1765];

    auto g_zz_0_0_0_y_y_xx_xz = buffer_2000_ppdd[1766];

    auto g_zz_0_0_0_y_y_xx_yy = buffer_2000_ppdd[1767];

    auto g_zz_0_0_0_y_y_xx_yz = buffer_2000_ppdd[1768];

    auto g_zz_0_0_0_y_y_xx_zz = buffer_2000_ppdd[1769];

    auto g_zz_0_0_0_y_y_xy_xx = buffer_2000_ppdd[1770];

    auto g_zz_0_0_0_y_y_xy_xy = buffer_2000_ppdd[1771];

    auto g_zz_0_0_0_y_y_xy_xz = buffer_2000_ppdd[1772];

    auto g_zz_0_0_0_y_y_xy_yy = buffer_2000_ppdd[1773];

    auto g_zz_0_0_0_y_y_xy_yz = buffer_2000_ppdd[1774];

    auto g_zz_0_0_0_y_y_xy_zz = buffer_2000_ppdd[1775];

    auto g_zz_0_0_0_y_y_xz_xx = buffer_2000_ppdd[1776];

    auto g_zz_0_0_0_y_y_xz_xy = buffer_2000_ppdd[1777];

    auto g_zz_0_0_0_y_y_xz_xz = buffer_2000_ppdd[1778];

    auto g_zz_0_0_0_y_y_xz_yy = buffer_2000_ppdd[1779];

    auto g_zz_0_0_0_y_y_xz_yz = buffer_2000_ppdd[1780];

    auto g_zz_0_0_0_y_y_xz_zz = buffer_2000_ppdd[1781];

    auto g_zz_0_0_0_y_y_yy_xx = buffer_2000_ppdd[1782];

    auto g_zz_0_0_0_y_y_yy_xy = buffer_2000_ppdd[1783];

    auto g_zz_0_0_0_y_y_yy_xz = buffer_2000_ppdd[1784];

    auto g_zz_0_0_0_y_y_yy_yy = buffer_2000_ppdd[1785];

    auto g_zz_0_0_0_y_y_yy_yz = buffer_2000_ppdd[1786];

    auto g_zz_0_0_0_y_y_yy_zz = buffer_2000_ppdd[1787];

    auto g_zz_0_0_0_y_y_yz_xx = buffer_2000_ppdd[1788];

    auto g_zz_0_0_0_y_y_yz_xy = buffer_2000_ppdd[1789];

    auto g_zz_0_0_0_y_y_yz_xz = buffer_2000_ppdd[1790];

    auto g_zz_0_0_0_y_y_yz_yy = buffer_2000_ppdd[1791];

    auto g_zz_0_0_0_y_y_yz_yz = buffer_2000_ppdd[1792];

    auto g_zz_0_0_0_y_y_yz_zz = buffer_2000_ppdd[1793];

    auto g_zz_0_0_0_y_y_zz_xx = buffer_2000_ppdd[1794];

    auto g_zz_0_0_0_y_y_zz_xy = buffer_2000_ppdd[1795];

    auto g_zz_0_0_0_y_y_zz_xz = buffer_2000_ppdd[1796];

    auto g_zz_0_0_0_y_y_zz_yy = buffer_2000_ppdd[1797];

    auto g_zz_0_0_0_y_y_zz_yz = buffer_2000_ppdd[1798];

    auto g_zz_0_0_0_y_y_zz_zz = buffer_2000_ppdd[1799];

    auto g_zz_0_0_0_y_z_xx_xx = buffer_2000_ppdd[1800];

    auto g_zz_0_0_0_y_z_xx_xy = buffer_2000_ppdd[1801];

    auto g_zz_0_0_0_y_z_xx_xz = buffer_2000_ppdd[1802];

    auto g_zz_0_0_0_y_z_xx_yy = buffer_2000_ppdd[1803];

    auto g_zz_0_0_0_y_z_xx_yz = buffer_2000_ppdd[1804];

    auto g_zz_0_0_0_y_z_xx_zz = buffer_2000_ppdd[1805];

    auto g_zz_0_0_0_y_z_xy_xx = buffer_2000_ppdd[1806];

    auto g_zz_0_0_0_y_z_xy_xy = buffer_2000_ppdd[1807];

    auto g_zz_0_0_0_y_z_xy_xz = buffer_2000_ppdd[1808];

    auto g_zz_0_0_0_y_z_xy_yy = buffer_2000_ppdd[1809];

    auto g_zz_0_0_0_y_z_xy_yz = buffer_2000_ppdd[1810];

    auto g_zz_0_0_0_y_z_xy_zz = buffer_2000_ppdd[1811];

    auto g_zz_0_0_0_y_z_xz_xx = buffer_2000_ppdd[1812];

    auto g_zz_0_0_0_y_z_xz_xy = buffer_2000_ppdd[1813];

    auto g_zz_0_0_0_y_z_xz_xz = buffer_2000_ppdd[1814];

    auto g_zz_0_0_0_y_z_xz_yy = buffer_2000_ppdd[1815];

    auto g_zz_0_0_0_y_z_xz_yz = buffer_2000_ppdd[1816];

    auto g_zz_0_0_0_y_z_xz_zz = buffer_2000_ppdd[1817];

    auto g_zz_0_0_0_y_z_yy_xx = buffer_2000_ppdd[1818];

    auto g_zz_0_0_0_y_z_yy_xy = buffer_2000_ppdd[1819];

    auto g_zz_0_0_0_y_z_yy_xz = buffer_2000_ppdd[1820];

    auto g_zz_0_0_0_y_z_yy_yy = buffer_2000_ppdd[1821];

    auto g_zz_0_0_0_y_z_yy_yz = buffer_2000_ppdd[1822];

    auto g_zz_0_0_0_y_z_yy_zz = buffer_2000_ppdd[1823];

    auto g_zz_0_0_0_y_z_yz_xx = buffer_2000_ppdd[1824];

    auto g_zz_0_0_0_y_z_yz_xy = buffer_2000_ppdd[1825];

    auto g_zz_0_0_0_y_z_yz_xz = buffer_2000_ppdd[1826];

    auto g_zz_0_0_0_y_z_yz_yy = buffer_2000_ppdd[1827];

    auto g_zz_0_0_0_y_z_yz_yz = buffer_2000_ppdd[1828];

    auto g_zz_0_0_0_y_z_yz_zz = buffer_2000_ppdd[1829];

    auto g_zz_0_0_0_y_z_zz_xx = buffer_2000_ppdd[1830];

    auto g_zz_0_0_0_y_z_zz_xy = buffer_2000_ppdd[1831];

    auto g_zz_0_0_0_y_z_zz_xz = buffer_2000_ppdd[1832];

    auto g_zz_0_0_0_y_z_zz_yy = buffer_2000_ppdd[1833];

    auto g_zz_0_0_0_y_z_zz_yz = buffer_2000_ppdd[1834];

    auto g_zz_0_0_0_y_z_zz_zz = buffer_2000_ppdd[1835];

    auto g_zz_0_0_0_z_x_xx_xx = buffer_2000_ppdd[1836];

    auto g_zz_0_0_0_z_x_xx_xy = buffer_2000_ppdd[1837];

    auto g_zz_0_0_0_z_x_xx_xz = buffer_2000_ppdd[1838];

    auto g_zz_0_0_0_z_x_xx_yy = buffer_2000_ppdd[1839];

    auto g_zz_0_0_0_z_x_xx_yz = buffer_2000_ppdd[1840];

    auto g_zz_0_0_0_z_x_xx_zz = buffer_2000_ppdd[1841];

    auto g_zz_0_0_0_z_x_xy_xx = buffer_2000_ppdd[1842];

    auto g_zz_0_0_0_z_x_xy_xy = buffer_2000_ppdd[1843];

    auto g_zz_0_0_0_z_x_xy_xz = buffer_2000_ppdd[1844];

    auto g_zz_0_0_0_z_x_xy_yy = buffer_2000_ppdd[1845];

    auto g_zz_0_0_0_z_x_xy_yz = buffer_2000_ppdd[1846];

    auto g_zz_0_0_0_z_x_xy_zz = buffer_2000_ppdd[1847];

    auto g_zz_0_0_0_z_x_xz_xx = buffer_2000_ppdd[1848];

    auto g_zz_0_0_0_z_x_xz_xy = buffer_2000_ppdd[1849];

    auto g_zz_0_0_0_z_x_xz_xz = buffer_2000_ppdd[1850];

    auto g_zz_0_0_0_z_x_xz_yy = buffer_2000_ppdd[1851];

    auto g_zz_0_0_0_z_x_xz_yz = buffer_2000_ppdd[1852];

    auto g_zz_0_0_0_z_x_xz_zz = buffer_2000_ppdd[1853];

    auto g_zz_0_0_0_z_x_yy_xx = buffer_2000_ppdd[1854];

    auto g_zz_0_0_0_z_x_yy_xy = buffer_2000_ppdd[1855];

    auto g_zz_0_0_0_z_x_yy_xz = buffer_2000_ppdd[1856];

    auto g_zz_0_0_0_z_x_yy_yy = buffer_2000_ppdd[1857];

    auto g_zz_0_0_0_z_x_yy_yz = buffer_2000_ppdd[1858];

    auto g_zz_0_0_0_z_x_yy_zz = buffer_2000_ppdd[1859];

    auto g_zz_0_0_0_z_x_yz_xx = buffer_2000_ppdd[1860];

    auto g_zz_0_0_0_z_x_yz_xy = buffer_2000_ppdd[1861];

    auto g_zz_0_0_0_z_x_yz_xz = buffer_2000_ppdd[1862];

    auto g_zz_0_0_0_z_x_yz_yy = buffer_2000_ppdd[1863];

    auto g_zz_0_0_0_z_x_yz_yz = buffer_2000_ppdd[1864];

    auto g_zz_0_0_0_z_x_yz_zz = buffer_2000_ppdd[1865];

    auto g_zz_0_0_0_z_x_zz_xx = buffer_2000_ppdd[1866];

    auto g_zz_0_0_0_z_x_zz_xy = buffer_2000_ppdd[1867];

    auto g_zz_0_0_0_z_x_zz_xz = buffer_2000_ppdd[1868];

    auto g_zz_0_0_0_z_x_zz_yy = buffer_2000_ppdd[1869];

    auto g_zz_0_0_0_z_x_zz_yz = buffer_2000_ppdd[1870];

    auto g_zz_0_0_0_z_x_zz_zz = buffer_2000_ppdd[1871];

    auto g_zz_0_0_0_z_y_xx_xx = buffer_2000_ppdd[1872];

    auto g_zz_0_0_0_z_y_xx_xy = buffer_2000_ppdd[1873];

    auto g_zz_0_0_0_z_y_xx_xz = buffer_2000_ppdd[1874];

    auto g_zz_0_0_0_z_y_xx_yy = buffer_2000_ppdd[1875];

    auto g_zz_0_0_0_z_y_xx_yz = buffer_2000_ppdd[1876];

    auto g_zz_0_0_0_z_y_xx_zz = buffer_2000_ppdd[1877];

    auto g_zz_0_0_0_z_y_xy_xx = buffer_2000_ppdd[1878];

    auto g_zz_0_0_0_z_y_xy_xy = buffer_2000_ppdd[1879];

    auto g_zz_0_0_0_z_y_xy_xz = buffer_2000_ppdd[1880];

    auto g_zz_0_0_0_z_y_xy_yy = buffer_2000_ppdd[1881];

    auto g_zz_0_0_0_z_y_xy_yz = buffer_2000_ppdd[1882];

    auto g_zz_0_0_0_z_y_xy_zz = buffer_2000_ppdd[1883];

    auto g_zz_0_0_0_z_y_xz_xx = buffer_2000_ppdd[1884];

    auto g_zz_0_0_0_z_y_xz_xy = buffer_2000_ppdd[1885];

    auto g_zz_0_0_0_z_y_xz_xz = buffer_2000_ppdd[1886];

    auto g_zz_0_0_0_z_y_xz_yy = buffer_2000_ppdd[1887];

    auto g_zz_0_0_0_z_y_xz_yz = buffer_2000_ppdd[1888];

    auto g_zz_0_0_0_z_y_xz_zz = buffer_2000_ppdd[1889];

    auto g_zz_0_0_0_z_y_yy_xx = buffer_2000_ppdd[1890];

    auto g_zz_0_0_0_z_y_yy_xy = buffer_2000_ppdd[1891];

    auto g_zz_0_0_0_z_y_yy_xz = buffer_2000_ppdd[1892];

    auto g_zz_0_0_0_z_y_yy_yy = buffer_2000_ppdd[1893];

    auto g_zz_0_0_0_z_y_yy_yz = buffer_2000_ppdd[1894];

    auto g_zz_0_0_0_z_y_yy_zz = buffer_2000_ppdd[1895];

    auto g_zz_0_0_0_z_y_yz_xx = buffer_2000_ppdd[1896];

    auto g_zz_0_0_0_z_y_yz_xy = buffer_2000_ppdd[1897];

    auto g_zz_0_0_0_z_y_yz_xz = buffer_2000_ppdd[1898];

    auto g_zz_0_0_0_z_y_yz_yy = buffer_2000_ppdd[1899];

    auto g_zz_0_0_0_z_y_yz_yz = buffer_2000_ppdd[1900];

    auto g_zz_0_0_0_z_y_yz_zz = buffer_2000_ppdd[1901];

    auto g_zz_0_0_0_z_y_zz_xx = buffer_2000_ppdd[1902];

    auto g_zz_0_0_0_z_y_zz_xy = buffer_2000_ppdd[1903];

    auto g_zz_0_0_0_z_y_zz_xz = buffer_2000_ppdd[1904];

    auto g_zz_0_0_0_z_y_zz_yy = buffer_2000_ppdd[1905];

    auto g_zz_0_0_0_z_y_zz_yz = buffer_2000_ppdd[1906];

    auto g_zz_0_0_0_z_y_zz_zz = buffer_2000_ppdd[1907];

    auto g_zz_0_0_0_z_z_xx_xx = buffer_2000_ppdd[1908];

    auto g_zz_0_0_0_z_z_xx_xy = buffer_2000_ppdd[1909];

    auto g_zz_0_0_0_z_z_xx_xz = buffer_2000_ppdd[1910];

    auto g_zz_0_0_0_z_z_xx_yy = buffer_2000_ppdd[1911];

    auto g_zz_0_0_0_z_z_xx_yz = buffer_2000_ppdd[1912];

    auto g_zz_0_0_0_z_z_xx_zz = buffer_2000_ppdd[1913];

    auto g_zz_0_0_0_z_z_xy_xx = buffer_2000_ppdd[1914];

    auto g_zz_0_0_0_z_z_xy_xy = buffer_2000_ppdd[1915];

    auto g_zz_0_0_0_z_z_xy_xz = buffer_2000_ppdd[1916];

    auto g_zz_0_0_0_z_z_xy_yy = buffer_2000_ppdd[1917];

    auto g_zz_0_0_0_z_z_xy_yz = buffer_2000_ppdd[1918];

    auto g_zz_0_0_0_z_z_xy_zz = buffer_2000_ppdd[1919];

    auto g_zz_0_0_0_z_z_xz_xx = buffer_2000_ppdd[1920];

    auto g_zz_0_0_0_z_z_xz_xy = buffer_2000_ppdd[1921];

    auto g_zz_0_0_0_z_z_xz_xz = buffer_2000_ppdd[1922];

    auto g_zz_0_0_0_z_z_xz_yy = buffer_2000_ppdd[1923];

    auto g_zz_0_0_0_z_z_xz_yz = buffer_2000_ppdd[1924];

    auto g_zz_0_0_0_z_z_xz_zz = buffer_2000_ppdd[1925];

    auto g_zz_0_0_0_z_z_yy_xx = buffer_2000_ppdd[1926];

    auto g_zz_0_0_0_z_z_yy_xy = buffer_2000_ppdd[1927];

    auto g_zz_0_0_0_z_z_yy_xz = buffer_2000_ppdd[1928];

    auto g_zz_0_0_0_z_z_yy_yy = buffer_2000_ppdd[1929];

    auto g_zz_0_0_0_z_z_yy_yz = buffer_2000_ppdd[1930];

    auto g_zz_0_0_0_z_z_yy_zz = buffer_2000_ppdd[1931];

    auto g_zz_0_0_0_z_z_yz_xx = buffer_2000_ppdd[1932];

    auto g_zz_0_0_0_z_z_yz_xy = buffer_2000_ppdd[1933];

    auto g_zz_0_0_0_z_z_yz_xz = buffer_2000_ppdd[1934];

    auto g_zz_0_0_0_z_z_yz_yy = buffer_2000_ppdd[1935];

    auto g_zz_0_0_0_z_z_yz_yz = buffer_2000_ppdd[1936];

    auto g_zz_0_0_0_z_z_yz_zz = buffer_2000_ppdd[1937];

    auto g_zz_0_0_0_z_z_zz_xx = buffer_2000_ppdd[1938];

    auto g_zz_0_0_0_z_z_zz_xy = buffer_2000_ppdd[1939];

    auto g_zz_0_0_0_z_z_zz_xz = buffer_2000_ppdd[1940];

    auto g_zz_0_0_0_z_z_zz_yy = buffer_2000_ppdd[1941];

    auto g_zz_0_0_0_z_z_zz_yz = buffer_2000_ppdd[1942];

    auto g_zz_0_0_0_z_z_zz_zz = buffer_2000_ppdd[1943];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_x_xx_xx, g_x_x_xx_xy, g_x_x_xx_xz, g_x_x_xx_yy, g_x_x_xx_yz, g_x_x_xx_zz, g_xx_0_0_0_x_x_xx_xx, g_xx_0_0_0_x_x_xx_xy, g_xx_0_0_0_x_x_xx_xz, g_xx_0_0_0_x_x_xx_yy, g_xx_0_0_0_x_x_xx_yz, g_xx_0_0_0_x_x_xx_zz, g_xxx_x_xx_xx, g_xxx_x_xx_xy, g_xxx_x_xx_xz, g_xxx_x_xx_yy, g_xxx_x_xx_yz, g_xxx_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_x_xx_xx[i] = -6.0 * g_x_x_xx_xx[i] * a_exp + 4.0 * g_xxx_x_xx_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_xx_xy[i] = -6.0 * g_x_x_xx_xy[i] * a_exp + 4.0 * g_xxx_x_xx_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_xx_xz[i] = -6.0 * g_x_x_xx_xz[i] * a_exp + 4.0 * g_xxx_x_xx_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_xx_yy[i] = -6.0 * g_x_x_xx_yy[i] * a_exp + 4.0 * g_xxx_x_xx_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_xx_yz[i] = -6.0 * g_x_x_xx_yz[i] * a_exp + 4.0 * g_xxx_x_xx_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_xx_zz[i] = -6.0 * g_x_x_xx_zz[i] * a_exp + 4.0 * g_xxx_x_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_x_xy_xx, g_x_x_xy_xy, g_x_x_xy_xz, g_x_x_xy_yy, g_x_x_xy_yz, g_x_x_xy_zz, g_xx_0_0_0_x_x_xy_xx, g_xx_0_0_0_x_x_xy_xy, g_xx_0_0_0_x_x_xy_xz, g_xx_0_0_0_x_x_xy_yy, g_xx_0_0_0_x_x_xy_yz, g_xx_0_0_0_x_x_xy_zz, g_xxx_x_xy_xx, g_xxx_x_xy_xy, g_xxx_x_xy_xz, g_xxx_x_xy_yy, g_xxx_x_xy_yz, g_xxx_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_x_xy_xx[i] = -6.0 * g_x_x_xy_xx[i] * a_exp + 4.0 * g_xxx_x_xy_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_xy_xy[i] = -6.0 * g_x_x_xy_xy[i] * a_exp + 4.0 * g_xxx_x_xy_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_xy_xz[i] = -6.0 * g_x_x_xy_xz[i] * a_exp + 4.0 * g_xxx_x_xy_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_xy_yy[i] = -6.0 * g_x_x_xy_yy[i] * a_exp + 4.0 * g_xxx_x_xy_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_xy_yz[i] = -6.0 * g_x_x_xy_yz[i] * a_exp + 4.0 * g_xxx_x_xy_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_xy_zz[i] = -6.0 * g_x_x_xy_zz[i] * a_exp + 4.0 * g_xxx_x_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_x_xz_xx, g_x_x_xz_xy, g_x_x_xz_xz, g_x_x_xz_yy, g_x_x_xz_yz, g_x_x_xz_zz, g_xx_0_0_0_x_x_xz_xx, g_xx_0_0_0_x_x_xz_xy, g_xx_0_0_0_x_x_xz_xz, g_xx_0_0_0_x_x_xz_yy, g_xx_0_0_0_x_x_xz_yz, g_xx_0_0_0_x_x_xz_zz, g_xxx_x_xz_xx, g_xxx_x_xz_xy, g_xxx_x_xz_xz, g_xxx_x_xz_yy, g_xxx_x_xz_yz, g_xxx_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_x_xz_xx[i] = -6.0 * g_x_x_xz_xx[i] * a_exp + 4.0 * g_xxx_x_xz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_xz_xy[i] = -6.0 * g_x_x_xz_xy[i] * a_exp + 4.0 * g_xxx_x_xz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_xz_xz[i] = -6.0 * g_x_x_xz_xz[i] * a_exp + 4.0 * g_xxx_x_xz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_xz_yy[i] = -6.0 * g_x_x_xz_yy[i] * a_exp + 4.0 * g_xxx_x_xz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_xz_yz[i] = -6.0 * g_x_x_xz_yz[i] * a_exp + 4.0 * g_xxx_x_xz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_xz_zz[i] = -6.0 * g_x_x_xz_zz[i] * a_exp + 4.0 * g_xxx_x_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_x_x_yy_xx, g_x_x_yy_xy, g_x_x_yy_xz, g_x_x_yy_yy, g_x_x_yy_yz, g_x_x_yy_zz, g_xx_0_0_0_x_x_yy_xx, g_xx_0_0_0_x_x_yy_xy, g_xx_0_0_0_x_x_yy_xz, g_xx_0_0_0_x_x_yy_yy, g_xx_0_0_0_x_x_yy_yz, g_xx_0_0_0_x_x_yy_zz, g_xxx_x_yy_xx, g_xxx_x_yy_xy, g_xxx_x_yy_xz, g_xxx_x_yy_yy, g_xxx_x_yy_yz, g_xxx_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_x_yy_xx[i] = -6.0 * g_x_x_yy_xx[i] * a_exp + 4.0 * g_xxx_x_yy_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_yy_xy[i] = -6.0 * g_x_x_yy_xy[i] * a_exp + 4.0 * g_xxx_x_yy_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_yy_xz[i] = -6.0 * g_x_x_yy_xz[i] * a_exp + 4.0 * g_xxx_x_yy_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_yy_yy[i] = -6.0 * g_x_x_yy_yy[i] * a_exp + 4.0 * g_xxx_x_yy_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_yy_yz[i] = -6.0 * g_x_x_yy_yz[i] * a_exp + 4.0 * g_xxx_x_yy_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_yy_zz[i] = -6.0 * g_x_x_yy_zz[i] * a_exp + 4.0 * g_xxx_x_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_x_x_yz_xx, g_x_x_yz_xy, g_x_x_yz_xz, g_x_x_yz_yy, g_x_x_yz_yz, g_x_x_yz_zz, g_xx_0_0_0_x_x_yz_xx, g_xx_0_0_0_x_x_yz_xy, g_xx_0_0_0_x_x_yz_xz, g_xx_0_0_0_x_x_yz_yy, g_xx_0_0_0_x_x_yz_yz, g_xx_0_0_0_x_x_yz_zz, g_xxx_x_yz_xx, g_xxx_x_yz_xy, g_xxx_x_yz_xz, g_xxx_x_yz_yy, g_xxx_x_yz_yz, g_xxx_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_x_yz_xx[i] = -6.0 * g_x_x_yz_xx[i] * a_exp + 4.0 * g_xxx_x_yz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_yz_xy[i] = -6.0 * g_x_x_yz_xy[i] * a_exp + 4.0 * g_xxx_x_yz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_yz_xz[i] = -6.0 * g_x_x_yz_xz[i] * a_exp + 4.0 * g_xxx_x_yz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_yz_yy[i] = -6.0 * g_x_x_yz_yy[i] * a_exp + 4.0 * g_xxx_x_yz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_yz_yz[i] = -6.0 * g_x_x_yz_yz[i] * a_exp + 4.0 * g_xxx_x_yz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_yz_zz[i] = -6.0 * g_x_x_yz_zz[i] * a_exp + 4.0 * g_xxx_x_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_x_x_zz_xx, g_x_x_zz_xy, g_x_x_zz_xz, g_x_x_zz_yy, g_x_x_zz_yz, g_x_x_zz_zz, g_xx_0_0_0_x_x_zz_xx, g_xx_0_0_0_x_x_zz_xy, g_xx_0_0_0_x_x_zz_xz, g_xx_0_0_0_x_x_zz_yy, g_xx_0_0_0_x_x_zz_yz, g_xx_0_0_0_x_x_zz_zz, g_xxx_x_zz_xx, g_xxx_x_zz_xy, g_xxx_x_zz_xz, g_xxx_x_zz_yy, g_xxx_x_zz_yz, g_xxx_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_x_zz_xx[i] = -6.0 * g_x_x_zz_xx[i] * a_exp + 4.0 * g_xxx_x_zz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_zz_xy[i] = -6.0 * g_x_x_zz_xy[i] * a_exp + 4.0 * g_xxx_x_zz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_zz_xz[i] = -6.0 * g_x_x_zz_xz[i] * a_exp + 4.0 * g_xxx_x_zz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_zz_yy[i] = -6.0 * g_x_x_zz_yy[i] * a_exp + 4.0 * g_xxx_x_zz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_zz_yz[i] = -6.0 * g_x_x_zz_yz[i] * a_exp + 4.0 * g_xxx_x_zz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_zz_zz[i] = -6.0 * g_x_x_zz_zz[i] * a_exp + 4.0 * g_xxx_x_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_x_y_xx_xx, g_x_y_xx_xy, g_x_y_xx_xz, g_x_y_xx_yy, g_x_y_xx_yz, g_x_y_xx_zz, g_xx_0_0_0_x_y_xx_xx, g_xx_0_0_0_x_y_xx_xy, g_xx_0_0_0_x_y_xx_xz, g_xx_0_0_0_x_y_xx_yy, g_xx_0_0_0_x_y_xx_yz, g_xx_0_0_0_x_y_xx_zz, g_xxx_y_xx_xx, g_xxx_y_xx_xy, g_xxx_y_xx_xz, g_xxx_y_xx_yy, g_xxx_y_xx_yz, g_xxx_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_y_xx_xx[i] = -6.0 * g_x_y_xx_xx[i] * a_exp + 4.0 * g_xxx_y_xx_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_xx_xy[i] = -6.0 * g_x_y_xx_xy[i] * a_exp + 4.0 * g_xxx_y_xx_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_xx_xz[i] = -6.0 * g_x_y_xx_xz[i] * a_exp + 4.0 * g_xxx_y_xx_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_xx_yy[i] = -6.0 * g_x_y_xx_yy[i] * a_exp + 4.0 * g_xxx_y_xx_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_xx_yz[i] = -6.0 * g_x_y_xx_yz[i] * a_exp + 4.0 * g_xxx_y_xx_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_xx_zz[i] = -6.0 * g_x_y_xx_zz[i] * a_exp + 4.0 * g_xxx_y_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_x_y_xy_xx, g_x_y_xy_xy, g_x_y_xy_xz, g_x_y_xy_yy, g_x_y_xy_yz, g_x_y_xy_zz, g_xx_0_0_0_x_y_xy_xx, g_xx_0_0_0_x_y_xy_xy, g_xx_0_0_0_x_y_xy_xz, g_xx_0_0_0_x_y_xy_yy, g_xx_0_0_0_x_y_xy_yz, g_xx_0_0_0_x_y_xy_zz, g_xxx_y_xy_xx, g_xxx_y_xy_xy, g_xxx_y_xy_xz, g_xxx_y_xy_yy, g_xxx_y_xy_yz, g_xxx_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_y_xy_xx[i] = -6.0 * g_x_y_xy_xx[i] * a_exp + 4.0 * g_xxx_y_xy_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_xy_xy[i] = -6.0 * g_x_y_xy_xy[i] * a_exp + 4.0 * g_xxx_y_xy_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_xy_xz[i] = -6.0 * g_x_y_xy_xz[i] * a_exp + 4.0 * g_xxx_y_xy_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_xy_yy[i] = -6.0 * g_x_y_xy_yy[i] * a_exp + 4.0 * g_xxx_y_xy_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_xy_yz[i] = -6.0 * g_x_y_xy_yz[i] * a_exp + 4.0 * g_xxx_y_xy_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_xy_zz[i] = -6.0 * g_x_y_xy_zz[i] * a_exp + 4.0 * g_xxx_y_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_x_y_xz_xx, g_x_y_xz_xy, g_x_y_xz_xz, g_x_y_xz_yy, g_x_y_xz_yz, g_x_y_xz_zz, g_xx_0_0_0_x_y_xz_xx, g_xx_0_0_0_x_y_xz_xy, g_xx_0_0_0_x_y_xz_xz, g_xx_0_0_0_x_y_xz_yy, g_xx_0_0_0_x_y_xz_yz, g_xx_0_0_0_x_y_xz_zz, g_xxx_y_xz_xx, g_xxx_y_xz_xy, g_xxx_y_xz_xz, g_xxx_y_xz_yy, g_xxx_y_xz_yz, g_xxx_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_y_xz_xx[i] = -6.0 * g_x_y_xz_xx[i] * a_exp + 4.0 * g_xxx_y_xz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_xz_xy[i] = -6.0 * g_x_y_xz_xy[i] * a_exp + 4.0 * g_xxx_y_xz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_xz_xz[i] = -6.0 * g_x_y_xz_xz[i] * a_exp + 4.0 * g_xxx_y_xz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_xz_yy[i] = -6.0 * g_x_y_xz_yy[i] * a_exp + 4.0 * g_xxx_y_xz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_xz_yz[i] = -6.0 * g_x_y_xz_yz[i] * a_exp + 4.0 * g_xxx_y_xz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_xz_zz[i] = -6.0 * g_x_y_xz_zz[i] * a_exp + 4.0 * g_xxx_y_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_x_y_yy_xx, g_x_y_yy_xy, g_x_y_yy_xz, g_x_y_yy_yy, g_x_y_yy_yz, g_x_y_yy_zz, g_xx_0_0_0_x_y_yy_xx, g_xx_0_0_0_x_y_yy_xy, g_xx_0_0_0_x_y_yy_xz, g_xx_0_0_0_x_y_yy_yy, g_xx_0_0_0_x_y_yy_yz, g_xx_0_0_0_x_y_yy_zz, g_xxx_y_yy_xx, g_xxx_y_yy_xy, g_xxx_y_yy_xz, g_xxx_y_yy_yy, g_xxx_y_yy_yz, g_xxx_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_y_yy_xx[i] = -6.0 * g_x_y_yy_xx[i] * a_exp + 4.0 * g_xxx_y_yy_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_yy_xy[i] = -6.0 * g_x_y_yy_xy[i] * a_exp + 4.0 * g_xxx_y_yy_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_yy_xz[i] = -6.0 * g_x_y_yy_xz[i] * a_exp + 4.0 * g_xxx_y_yy_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_yy_yy[i] = -6.0 * g_x_y_yy_yy[i] * a_exp + 4.0 * g_xxx_y_yy_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_yy_yz[i] = -6.0 * g_x_y_yy_yz[i] * a_exp + 4.0 * g_xxx_y_yy_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_yy_zz[i] = -6.0 * g_x_y_yy_zz[i] * a_exp + 4.0 * g_xxx_y_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_x_y_yz_xx, g_x_y_yz_xy, g_x_y_yz_xz, g_x_y_yz_yy, g_x_y_yz_yz, g_x_y_yz_zz, g_xx_0_0_0_x_y_yz_xx, g_xx_0_0_0_x_y_yz_xy, g_xx_0_0_0_x_y_yz_xz, g_xx_0_0_0_x_y_yz_yy, g_xx_0_0_0_x_y_yz_yz, g_xx_0_0_0_x_y_yz_zz, g_xxx_y_yz_xx, g_xxx_y_yz_xy, g_xxx_y_yz_xz, g_xxx_y_yz_yy, g_xxx_y_yz_yz, g_xxx_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_y_yz_xx[i] = -6.0 * g_x_y_yz_xx[i] * a_exp + 4.0 * g_xxx_y_yz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_yz_xy[i] = -6.0 * g_x_y_yz_xy[i] * a_exp + 4.0 * g_xxx_y_yz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_yz_xz[i] = -6.0 * g_x_y_yz_xz[i] * a_exp + 4.0 * g_xxx_y_yz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_yz_yy[i] = -6.0 * g_x_y_yz_yy[i] * a_exp + 4.0 * g_xxx_y_yz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_yz_yz[i] = -6.0 * g_x_y_yz_yz[i] * a_exp + 4.0 * g_xxx_y_yz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_yz_zz[i] = -6.0 * g_x_y_yz_zz[i] * a_exp + 4.0 * g_xxx_y_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_x_y_zz_xx, g_x_y_zz_xy, g_x_y_zz_xz, g_x_y_zz_yy, g_x_y_zz_yz, g_x_y_zz_zz, g_xx_0_0_0_x_y_zz_xx, g_xx_0_0_0_x_y_zz_xy, g_xx_0_0_0_x_y_zz_xz, g_xx_0_0_0_x_y_zz_yy, g_xx_0_0_0_x_y_zz_yz, g_xx_0_0_0_x_y_zz_zz, g_xxx_y_zz_xx, g_xxx_y_zz_xy, g_xxx_y_zz_xz, g_xxx_y_zz_yy, g_xxx_y_zz_yz, g_xxx_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_y_zz_xx[i] = -6.0 * g_x_y_zz_xx[i] * a_exp + 4.0 * g_xxx_y_zz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_zz_xy[i] = -6.0 * g_x_y_zz_xy[i] * a_exp + 4.0 * g_xxx_y_zz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_zz_xz[i] = -6.0 * g_x_y_zz_xz[i] * a_exp + 4.0 * g_xxx_y_zz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_zz_yy[i] = -6.0 * g_x_y_zz_yy[i] * a_exp + 4.0 * g_xxx_y_zz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_zz_yz[i] = -6.0 * g_x_y_zz_yz[i] * a_exp + 4.0 * g_xxx_y_zz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_zz_zz[i] = -6.0 * g_x_y_zz_zz[i] * a_exp + 4.0 * g_xxx_y_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_x_z_xx_xx, g_x_z_xx_xy, g_x_z_xx_xz, g_x_z_xx_yy, g_x_z_xx_yz, g_x_z_xx_zz, g_xx_0_0_0_x_z_xx_xx, g_xx_0_0_0_x_z_xx_xy, g_xx_0_0_0_x_z_xx_xz, g_xx_0_0_0_x_z_xx_yy, g_xx_0_0_0_x_z_xx_yz, g_xx_0_0_0_x_z_xx_zz, g_xxx_z_xx_xx, g_xxx_z_xx_xy, g_xxx_z_xx_xz, g_xxx_z_xx_yy, g_xxx_z_xx_yz, g_xxx_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_z_xx_xx[i] = -6.0 * g_x_z_xx_xx[i] * a_exp + 4.0 * g_xxx_z_xx_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_xx_xy[i] = -6.0 * g_x_z_xx_xy[i] * a_exp + 4.0 * g_xxx_z_xx_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_xx_xz[i] = -6.0 * g_x_z_xx_xz[i] * a_exp + 4.0 * g_xxx_z_xx_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_xx_yy[i] = -6.0 * g_x_z_xx_yy[i] * a_exp + 4.0 * g_xxx_z_xx_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_xx_yz[i] = -6.0 * g_x_z_xx_yz[i] * a_exp + 4.0 * g_xxx_z_xx_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_xx_zz[i] = -6.0 * g_x_z_xx_zz[i] * a_exp + 4.0 * g_xxx_z_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_x_z_xy_xx, g_x_z_xy_xy, g_x_z_xy_xz, g_x_z_xy_yy, g_x_z_xy_yz, g_x_z_xy_zz, g_xx_0_0_0_x_z_xy_xx, g_xx_0_0_0_x_z_xy_xy, g_xx_0_0_0_x_z_xy_xz, g_xx_0_0_0_x_z_xy_yy, g_xx_0_0_0_x_z_xy_yz, g_xx_0_0_0_x_z_xy_zz, g_xxx_z_xy_xx, g_xxx_z_xy_xy, g_xxx_z_xy_xz, g_xxx_z_xy_yy, g_xxx_z_xy_yz, g_xxx_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_z_xy_xx[i] = -6.0 * g_x_z_xy_xx[i] * a_exp + 4.0 * g_xxx_z_xy_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_xy_xy[i] = -6.0 * g_x_z_xy_xy[i] * a_exp + 4.0 * g_xxx_z_xy_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_xy_xz[i] = -6.0 * g_x_z_xy_xz[i] * a_exp + 4.0 * g_xxx_z_xy_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_xy_yy[i] = -6.0 * g_x_z_xy_yy[i] * a_exp + 4.0 * g_xxx_z_xy_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_xy_yz[i] = -6.0 * g_x_z_xy_yz[i] * a_exp + 4.0 * g_xxx_z_xy_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_xy_zz[i] = -6.0 * g_x_z_xy_zz[i] * a_exp + 4.0 * g_xxx_z_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_x_z_xz_xx, g_x_z_xz_xy, g_x_z_xz_xz, g_x_z_xz_yy, g_x_z_xz_yz, g_x_z_xz_zz, g_xx_0_0_0_x_z_xz_xx, g_xx_0_0_0_x_z_xz_xy, g_xx_0_0_0_x_z_xz_xz, g_xx_0_0_0_x_z_xz_yy, g_xx_0_0_0_x_z_xz_yz, g_xx_0_0_0_x_z_xz_zz, g_xxx_z_xz_xx, g_xxx_z_xz_xy, g_xxx_z_xz_xz, g_xxx_z_xz_yy, g_xxx_z_xz_yz, g_xxx_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_z_xz_xx[i] = -6.0 * g_x_z_xz_xx[i] * a_exp + 4.0 * g_xxx_z_xz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_xz_xy[i] = -6.0 * g_x_z_xz_xy[i] * a_exp + 4.0 * g_xxx_z_xz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_xz_xz[i] = -6.0 * g_x_z_xz_xz[i] * a_exp + 4.0 * g_xxx_z_xz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_xz_yy[i] = -6.0 * g_x_z_xz_yy[i] * a_exp + 4.0 * g_xxx_z_xz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_xz_yz[i] = -6.0 * g_x_z_xz_yz[i] * a_exp + 4.0 * g_xxx_z_xz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_xz_zz[i] = -6.0 * g_x_z_xz_zz[i] * a_exp + 4.0 * g_xxx_z_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_x_z_yy_xx, g_x_z_yy_xy, g_x_z_yy_xz, g_x_z_yy_yy, g_x_z_yy_yz, g_x_z_yy_zz, g_xx_0_0_0_x_z_yy_xx, g_xx_0_0_0_x_z_yy_xy, g_xx_0_0_0_x_z_yy_xz, g_xx_0_0_0_x_z_yy_yy, g_xx_0_0_0_x_z_yy_yz, g_xx_0_0_0_x_z_yy_zz, g_xxx_z_yy_xx, g_xxx_z_yy_xy, g_xxx_z_yy_xz, g_xxx_z_yy_yy, g_xxx_z_yy_yz, g_xxx_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_z_yy_xx[i] = -6.0 * g_x_z_yy_xx[i] * a_exp + 4.0 * g_xxx_z_yy_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_yy_xy[i] = -6.0 * g_x_z_yy_xy[i] * a_exp + 4.0 * g_xxx_z_yy_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_yy_xz[i] = -6.0 * g_x_z_yy_xz[i] * a_exp + 4.0 * g_xxx_z_yy_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_yy_yy[i] = -6.0 * g_x_z_yy_yy[i] * a_exp + 4.0 * g_xxx_z_yy_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_yy_yz[i] = -6.0 * g_x_z_yy_yz[i] * a_exp + 4.0 * g_xxx_z_yy_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_yy_zz[i] = -6.0 * g_x_z_yy_zz[i] * a_exp + 4.0 * g_xxx_z_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_x_z_yz_xx, g_x_z_yz_xy, g_x_z_yz_xz, g_x_z_yz_yy, g_x_z_yz_yz, g_x_z_yz_zz, g_xx_0_0_0_x_z_yz_xx, g_xx_0_0_0_x_z_yz_xy, g_xx_0_0_0_x_z_yz_xz, g_xx_0_0_0_x_z_yz_yy, g_xx_0_0_0_x_z_yz_yz, g_xx_0_0_0_x_z_yz_zz, g_xxx_z_yz_xx, g_xxx_z_yz_xy, g_xxx_z_yz_xz, g_xxx_z_yz_yy, g_xxx_z_yz_yz, g_xxx_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_z_yz_xx[i] = -6.0 * g_x_z_yz_xx[i] * a_exp + 4.0 * g_xxx_z_yz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_yz_xy[i] = -6.0 * g_x_z_yz_xy[i] * a_exp + 4.0 * g_xxx_z_yz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_yz_xz[i] = -6.0 * g_x_z_yz_xz[i] * a_exp + 4.0 * g_xxx_z_yz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_yz_yy[i] = -6.0 * g_x_z_yz_yy[i] * a_exp + 4.0 * g_xxx_z_yz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_yz_yz[i] = -6.0 * g_x_z_yz_yz[i] * a_exp + 4.0 * g_xxx_z_yz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_yz_zz[i] = -6.0 * g_x_z_yz_zz[i] * a_exp + 4.0 * g_xxx_z_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_x_z_zz_xx, g_x_z_zz_xy, g_x_z_zz_xz, g_x_z_zz_yy, g_x_z_zz_yz, g_x_z_zz_zz, g_xx_0_0_0_x_z_zz_xx, g_xx_0_0_0_x_z_zz_xy, g_xx_0_0_0_x_z_zz_xz, g_xx_0_0_0_x_z_zz_yy, g_xx_0_0_0_x_z_zz_yz, g_xx_0_0_0_x_z_zz_zz, g_xxx_z_zz_xx, g_xxx_z_zz_xy, g_xxx_z_zz_xz, g_xxx_z_zz_yy, g_xxx_z_zz_yz, g_xxx_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_z_zz_xx[i] = -6.0 * g_x_z_zz_xx[i] * a_exp + 4.0 * g_xxx_z_zz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_zz_xy[i] = -6.0 * g_x_z_zz_xy[i] * a_exp + 4.0 * g_xxx_z_zz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_zz_xz[i] = -6.0 * g_x_z_zz_xz[i] * a_exp + 4.0 * g_xxx_z_zz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_zz_yy[i] = -6.0 * g_x_z_zz_yy[i] * a_exp + 4.0 * g_xxx_z_zz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_zz_yz[i] = -6.0 * g_x_z_zz_yz[i] * a_exp + 4.0 * g_xxx_z_zz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_zz_zz[i] = -6.0 * g_x_z_zz_zz[i] * a_exp + 4.0 * g_xxx_z_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_xx_0_0_0_y_x_xx_xx, g_xx_0_0_0_y_x_xx_xy, g_xx_0_0_0_y_x_xx_xz, g_xx_0_0_0_y_x_xx_yy, g_xx_0_0_0_y_x_xx_yz, g_xx_0_0_0_y_x_xx_zz, g_xxy_x_xx_xx, g_xxy_x_xx_xy, g_xxy_x_xx_xz, g_xxy_x_xx_yy, g_xxy_x_xx_yz, g_xxy_x_xx_zz, g_y_x_xx_xx, g_y_x_xx_xy, g_y_x_xx_xz, g_y_x_xx_yy, g_y_x_xx_yz, g_y_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_x_xx_xx[i] = -2.0 * g_y_x_xx_xx[i] * a_exp + 4.0 * g_xxy_x_xx_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_xx_xy[i] = -2.0 * g_y_x_xx_xy[i] * a_exp + 4.0 * g_xxy_x_xx_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_xx_xz[i] = -2.0 * g_y_x_xx_xz[i] * a_exp + 4.0 * g_xxy_x_xx_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_xx_yy[i] = -2.0 * g_y_x_xx_yy[i] * a_exp + 4.0 * g_xxy_x_xx_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_xx_yz[i] = -2.0 * g_y_x_xx_yz[i] * a_exp + 4.0 * g_xxy_x_xx_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_xx_zz[i] = -2.0 * g_y_x_xx_zz[i] * a_exp + 4.0 * g_xxy_x_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_xx_0_0_0_y_x_xy_xx, g_xx_0_0_0_y_x_xy_xy, g_xx_0_0_0_y_x_xy_xz, g_xx_0_0_0_y_x_xy_yy, g_xx_0_0_0_y_x_xy_yz, g_xx_0_0_0_y_x_xy_zz, g_xxy_x_xy_xx, g_xxy_x_xy_xy, g_xxy_x_xy_xz, g_xxy_x_xy_yy, g_xxy_x_xy_yz, g_xxy_x_xy_zz, g_y_x_xy_xx, g_y_x_xy_xy, g_y_x_xy_xz, g_y_x_xy_yy, g_y_x_xy_yz, g_y_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_x_xy_xx[i] = -2.0 * g_y_x_xy_xx[i] * a_exp + 4.0 * g_xxy_x_xy_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_xy_xy[i] = -2.0 * g_y_x_xy_xy[i] * a_exp + 4.0 * g_xxy_x_xy_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_xy_xz[i] = -2.0 * g_y_x_xy_xz[i] * a_exp + 4.0 * g_xxy_x_xy_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_xy_yy[i] = -2.0 * g_y_x_xy_yy[i] * a_exp + 4.0 * g_xxy_x_xy_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_xy_yz[i] = -2.0 * g_y_x_xy_yz[i] * a_exp + 4.0 * g_xxy_x_xy_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_xy_zz[i] = -2.0 * g_y_x_xy_zz[i] * a_exp + 4.0 * g_xxy_x_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_xx_0_0_0_y_x_xz_xx, g_xx_0_0_0_y_x_xz_xy, g_xx_0_0_0_y_x_xz_xz, g_xx_0_0_0_y_x_xz_yy, g_xx_0_0_0_y_x_xz_yz, g_xx_0_0_0_y_x_xz_zz, g_xxy_x_xz_xx, g_xxy_x_xz_xy, g_xxy_x_xz_xz, g_xxy_x_xz_yy, g_xxy_x_xz_yz, g_xxy_x_xz_zz, g_y_x_xz_xx, g_y_x_xz_xy, g_y_x_xz_xz, g_y_x_xz_yy, g_y_x_xz_yz, g_y_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_x_xz_xx[i] = -2.0 * g_y_x_xz_xx[i] * a_exp + 4.0 * g_xxy_x_xz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_xz_xy[i] = -2.0 * g_y_x_xz_xy[i] * a_exp + 4.0 * g_xxy_x_xz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_xz_xz[i] = -2.0 * g_y_x_xz_xz[i] * a_exp + 4.0 * g_xxy_x_xz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_xz_yy[i] = -2.0 * g_y_x_xz_yy[i] * a_exp + 4.0 * g_xxy_x_xz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_xz_yz[i] = -2.0 * g_y_x_xz_yz[i] * a_exp + 4.0 * g_xxy_x_xz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_xz_zz[i] = -2.0 * g_y_x_xz_zz[i] * a_exp + 4.0 * g_xxy_x_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_xx_0_0_0_y_x_yy_xx, g_xx_0_0_0_y_x_yy_xy, g_xx_0_0_0_y_x_yy_xz, g_xx_0_0_0_y_x_yy_yy, g_xx_0_0_0_y_x_yy_yz, g_xx_0_0_0_y_x_yy_zz, g_xxy_x_yy_xx, g_xxy_x_yy_xy, g_xxy_x_yy_xz, g_xxy_x_yy_yy, g_xxy_x_yy_yz, g_xxy_x_yy_zz, g_y_x_yy_xx, g_y_x_yy_xy, g_y_x_yy_xz, g_y_x_yy_yy, g_y_x_yy_yz, g_y_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_x_yy_xx[i] = -2.0 * g_y_x_yy_xx[i] * a_exp + 4.0 * g_xxy_x_yy_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_yy_xy[i] = -2.0 * g_y_x_yy_xy[i] * a_exp + 4.0 * g_xxy_x_yy_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_yy_xz[i] = -2.0 * g_y_x_yy_xz[i] * a_exp + 4.0 * g_xxy_x_yy_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_yy_yy[i] = -2.0 * g_y_x_yy_yy[i] * a_exp + 4.0 * g_xxy_x_yy_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_yy_yz[i] = -2.0 * g_y_x_yy_yz[i] * a_exp + 4.0 * g_xxy_x_yy_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_yy_zz[i] = -2.0 * g_y_x_yy_zz[i] * a_exp + 4.0 * g_xxy_x_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_xx_0_0_0_y_x_yz_xx, g_xx_0_0_0_y_x_yz_xy, g_xx_0_0_0_y_x_yz_xz, g_xx_0_0_0_y_x_yz_yy, g_xx_0_0_0_y_x_yz_yz, g_xx_0_0_0_y_x_yz_zz, g_xxy_x_yz_xx, g_xxy_x_yz_xy, g_xxy_x_yz_xz, g_xxy_x_yz_yy, g_xxy_x_yz_yz, g_xxy_x_yz_zz, g_y_x_yz_xx, g_y_x_yz_xy, g_y_x_yz_xz, g_y_x_yz_yy, g_y_x_yz_yz, g_y_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_x_yz_xx[i] = -2.0 * g_y_x_yz_xx[i] * a_exp + 4.0 * g_xxy_x_yz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_yz_xy[i] = -2.0 * g_y_x_yz_xy[i] * a_exp + 4.0 * g_xxy_x_yz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_yz_xz[i] = -2.0 * g_y_x_yz_xz[i] * a_exp + 4.0 * g_xxy_x_yz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_yz_yy[i] = -2.0 * g_y_x_yz_yy[i] * a_exp + 4.0 * g_xxy_x_yz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_yz_yz[i] = -2.0 * g_y_x_yz_yz[i] * a_exp + 4.0 * g_xxy_x_yz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_yz_zz[i] = -2.0 * g_y_x_yz_zz[i] * a_exp + 4.0 * g_xxy_x_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_xx_0_0_0_y_x_zz_xx, g_xx_0_0_0_y_x_zz_xy, g_xx_0_0_0_y_x_zz_xz, g_xx_0_0_0_y_x_zz_yy, g_xx_0_0_0_y_x_zz_yz, g_xx_0_0_0_y_x_zz_zz, g_xxy_x_zz_xx, g_xxy_x_zz_xy, g_xxy_x_zz_xz, g_xxy_x_zz_yy, g_xxy_x_zz_yz, g_xxy_x_zz_zz, g_y_x_zz_xx, g_y_x_zz_xy, g_y_x_zz_xz, g_y_x_zz_yy, g_y_x_zz_yz, g_y_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_x_zz_xx[i] = -2.0 * g_y_x_zz_xx[i] * a_exp + 4.0 * g_xxy_x_zz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_zz_xy[i] = -2.0 * g_y_x_zz_xy[i] * a_exp + 4.0 * g_xxy_x_zz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_zz_xz[i] = -2.0 * g_y_x_zz_xz[i] * a_exp + 4.0 * g_xxy_x_zz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_zz_yy[i] = -2.0 * g_y_x_zz_yy[i] * a_exp + 4.0 * g_xxy_x_zz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_zz_yz[i] = -2.0 * g_y_x_zz_yz[i] * a_exp + 4.0 * g_xxy_x_zz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_zz_zz[i] = -2.0 * g_y_x_zz_zz[i] * a_exp + 4.0 * g_xxy_x_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_xx_0_0_0_y_y_xx_xx, g_xx_0_0_0_y_y_xx_xy, g_xx_0_0_0_y_y_xx_xz, g_xx_0_0_0_y_y_xx_yy, g_xx_0_0_0_y_y_xx_yz, g_xx_0_0_0_y_y_xx_zz, g_xxy_y_xx_xx, g_xxy_y_xx_xy, g_xxy_y_xx_xz, g_xxy_y_xx_yy, g_xxy_y_xx_yz, g_xxy_y_xx_zz, g_y_y_xx_xx, g_y_y_xx_xy, g_y_y_xx_xz, g_y_y_xx_yy, g_y_y_xx_yz, g_y_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_y_xx_xx[i] = -2.0 * g_y_y_xx_xx[i] * a_exp + 4.0 * g_xxy_y_xx_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_xx_xy[i] = -2.0 * g_y_y_xx_xy[i] * a_exp + 4.0 * g_xxy_y_xx_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_xx_xz[i] = -2.0 * g_y_y_xx_xz[i] * a_exp + 4.0 * g_xxy_y_xx_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_xx_yy[i] = -2.0 * g_y_y_xx_yy[i] * a_exp + 4.0 * g_xxy_y_xx_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_xx_yz[i] = -2.0 * g_y_y_xx_yz[i] * a_exp + 4.0 * g_xxy_y_xx_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_xx_zz[i] = -2.0 * g_y_y_xx_zz[i] * a_exp + 4.0 * g_xxy_y_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_xx_0_0_0_y_y_xy_xx, g_xx_0_0_0_y_y_xy_xy, g_xx_0_0_0_y_y_xy_xz, g_xx_0_0_0_y_y_xy_yy, g_xx_0_0_0_y_y_xy_yz, g_xx_0_0_0_y_y_xy_zz, g_xxy_y_xy_xx, g_xxy_y_xy_xy, g_xxy_y_xy_xz, g_xxy_y_xy_yy, g_xxy_y_xy_yz, g_xxy_y_xy_zz, g_y_y_xy_xx, g_y_y_xy_xy, g_y_y_xy_xz, g_y_y_xy_yy, g_y_y_xy_yz, g_y_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_y_xy_xx[i] = -2.0 * g_y_y_xy_xx[i] * a_exp + 4.0 * g_xxy_y_xy_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_xy_xy[i] = -2.0 * g_y_y_xy_xy[i] * a_exp + 4.0 * g_xxy_y_xy_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_xy_xz[i] = -2.0 * g_y_y_xy_xz[i] * a_exp + 4.0 * g_xxy_y_xy_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_xy_yy[i] = -2.0 * g_y_y_xy_yy[i] * a_exp + 4.0 * g_xxy_y_xy_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_xy_yz[i] = -2.0 * g_y_y_xy_yz[i] * a_exp + 4.0 * g_xxy_y_xy_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_xy_zz[i] = -2.0 * g_y_y_xy_zz[i] * a_exp + 4.0 * g_xxy_y_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_xx_0_0_0_y_y_xz_xx, g_xx_0_0_0_y_y_xz_xy, g_xx_0_0_0_y_y_xz_xz, g_xx_0_0_0_y_y_xz_yy, g_xx_0_0_0_y_y_xz_yz, g_xx_0_0_0_y_y_xz_zz, g_xxy_y_xz_xx, g_xxy_y_xz_xy, g_xxy_y_xz_xz, g_xxy_y_xz_yy, g_xxy_y_xz_yz, g_xxy_y_xz_zz, g_y_y_xz_xx, g_y_y_xz_xy, g_y_y_xz_xz, g_y_y_xz_yy, g_y_y_xz_yz, g_y_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_y_xz_xx[i] = -2.0 * g_y_y_xz_xx[i] * a_exp + 4.0 * g_xxy_y_xz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_xz_xy[i] = -2.0 * g_y_y_xz_xy[i] * a_exp + 4.0 * g_xxy_y_xz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_xz_xz[i] = -2.0 * g_y_y_xz_xz[i] * a_exp + 4.0 * g_xxy_y_xz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_xz_yy[i] = -2.0 * g_y_y_xz_yy[i] * a_exp + 4.0 * g_xxy_y_xz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_xz_yz[i] = -2.0 * g_y_y_xz_yz[i] * a_exp + 4.0 * g_xxy_y_xz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_xz_zz[i] = -2.0 * g_y_y_xz_zz[i] * a_exp + 4.0 * g_xxy_y_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_xx_0_0_0_y_y_yy_xx, g_xx_0_0_0_y_y_yy_xy, g_xx_0_0_0_y_y_yy_xz, g_xx_0_0_0_y_y_yy_yy, g_xx_0_0_0_y_y_yy_yz, g_xx_0_0_0_y_y_yy_zz, g_xxy_y_yy_xx, g_xxy_y_yy_xy, g_xxy_y_yy_xz, g_xxy_y_yy_yy, g_xxy_y_yy_yz, g_xxy_y_yy_zz, g_y_y_yy_xx, g_y_y_yy_xy, g_y_y_yy_xz, g_y_y_yy_yy, g_y_y_yy_yz, g_y_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_y_yy_xx[i] = -2.0 * g_y_y_yy_xx[i] * a_exp + 4.0 * g_xxy_y_yy_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_yy_xy[i] = -2.0 * g_y_y_yy_xy[i] * a_exp + 4.0 * g_xxy_y_yy_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_yy_xz[i] = -2.0 * g_y_y_yy_xz[i] * a_exp + 4.0 * g_xxy_y_yy_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_yy_yy[i] = -2.0 * g_y_y_yy_yy[i] * a_exp + 4.0 * g_xxy_y_yy_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_yy_yz[i] = -2.0 * g_y_y_yy_yz[i] * a_exp + 4.0 * g_xxy_y_yy_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_yy_zz[i] = -2.0 * g_y_y_yy_zz[i] * a_exp + 4.0 * g_xxy_y_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_xx_0_0_0_y_y_yz_xx, g_xx_0_0_0_y_y_yz_xy, g_xx_0_0_0_y_y_yz_xz, g_xx_0_0_0_y_y_yz_yy, g_xx_0_0_0_y_y_yz_yz, g_xx_0_0_0_y_y_yz_zz, g_xxy_y_yz_xx, g_xxy_y_yz_xy, g_xxy_y_yz_xz, g_xxy_y_yz_yy, g_xxy_y_yz_yz, g_xxy_y_yz_zz, g_y_y_yz_xx, g_y_y_yz_xy, g_y_y_yz_xz, g_y_y_yz_yy, g_y_y_yz_yz, g_y_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_y_yz_xx[i] = -2.0 * g_y_y_yz_xx[i] * a_exp + 4.0 * g_xxy_y_yz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_yz_xy[i] = -2.0 * g_y_y_yz_xy[i] * a_exp + 4.0 * g_xxy_y_yz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_yz_xz[i] = -2.0 * g_y_y_yz_xz[i] * a_exp + 4.0 * g_xxy_y_yz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_yz_yy[i] = -2.0 * g_y_y_yz_yy[i] * a_exp + 4.0 * g_xxy_y_yz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_yz_yz[i] = -2.0 * g_y_y_yz_yz[i] * a_exp + 4.0 * g_xxy_y_yz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_yz_zz[i] = -2.0 * g_y_y_yz_zz[i] * a_exp + 4.0 * g_xxy_y_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_xx_0_0_0_y_y_zz_xx, g_xx_0_0_0_y_y_zz_xy, g_xx_0_0_0_y_y_zz_xz, g_xx_0_0_0_y_y_zz_yy, g_xx_0_0_0_y_y_zz_yz, g_xx_0_0_0_y_y_zz_zz, g_xxy_y_zz_xx, g_xxy_y_zz_xy, g_xxy_y_zz_xz, g_xxy_y_zz_yy, g_xxy_y_zz_yz, g_xxy_y_zz_zz, g_y_y_zz_xx, g_y_y_zz_xy, g_y_y_zz_xz, g_y_y_zz_yy, g_y_y_zz_yz, g_y_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_y_zz_xx[i] = -2.0 * g_y_y_zz_xx[i] * a_exp + 4.0 * g_xxy_y_zz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_zz_xy[i] = -2.0 * g_y_y_zz_xy[i] * a_exp + 4.0 * g_xxy_y_zz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_zz_xz[i] = -2.0 * g_y_y_zz_xz[i] * a_exp + 4.0 * g_xxy_y_zz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_zz_yy[i] = -2.0 * g_y_y_zz_yy[i] * a_exp + 4.0 * g_xxy_y_zz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_zz_yz[i] = -2.0 * g_y_y_zz_yz[i] * a_exp + 4.0 * g_xxy_y_zz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_zz_zz[i] = -2.0 * g_y_y_zz_zz[i] * a_exp + 4.0 * g_xxy_y_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_xx_0_0_0_y_z_xx_xx, g_xx_0_0_0_y_z_xx_xy, g_xx_0_0_0_y_z_xx_xz, g_xx_0_0_0_y_z_xx_yy, g_xx_0_0_0_y_z_xx_yz, g_xx_0_0_0_y_z_xx_zz, g_xxy_z_xx_xx, g_xxy_z_xx_xy, g_xxy_z_xx_xz, g_xxy_z_xx_yy, g_xxy_z_xx_yz, g_xxy_z_xx_zz, g_y_z_xx_xx, g_y_z_xx_xy, g_y_z_xx_xz, g_y_z_xx_yy, g_y_z_xx_yz, g_y_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_z_xx_xx[i] = -2.0 * g_y_z_xx_xx[i] * a_exp + 4.0 * g_xxy_z_xx_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_xx_xy[i] = -2.0 * g_y_z_xx_xy[i] * a_exp + 4.0 * g_xxy_z_xx_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_xx_xz[i] = -2.0 * g_y_z_xx_xz[i] * a_exp + 4.0 * g_xxy_z_xx_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_xx_yy[i] = -2.0 * g_y_z_xx_yy[i] * a_exp + 4.0 * g_xxy_z_xx_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_xx_yz[i] = -2.0 * g_y_z_xx_yz[i] * a_exp + 4.0 * g_xxy_z_xx_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_xx_zz[i] = -2.0 * g_y_z_xx_zz[i] * a_exp + 4.0 * g_xxy_z_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_xx_0_0_0_y_z_xy_xx, g_xx_0_0_0_y_z_xy_xy, g_xx_0_0_0_y_z_xy_xz, g_xx_0_0_0_y_z_xy_yy, g_xx_0_0_0_y_z_xy_yz, g_xx_0_0_0_y_z_xy_zz, g_xxy_z_xy_xx, g_xxy_z_xy_xy, g_xxy_z_xy_xz, g_xxy_z_xy_yy, g_xxy_z_xy_yz, g_xxy_z_xy_zz, g_y_z_xy_xx, g_y_z_xy_xy, g_y_z_xy_xz, g_y_z_xy_yy, g_y_z_xy_yz, g_y_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_z_xy_xx[i] = -2.0 * g_y_z_xy_xx[i] * a_exp + 4.0 * g_xxy_z_xy_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_xy_xy[i] = -2.0 * g_y_z_xy_xy[i] * a_exp + 4.0 * g_xxy_z_xy_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_xy_xz[i] = -2.0 * g_y_z_xy_xz[i] * a_exp + 4.0 * g_xxy_z_xy_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_xy_yy[i] = -2.0 * g_y_z_xy_yy[i] * a_exp + 4.0 * g_xxy_z_xy_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_xy_yz[i] = -2.0 * g_y_z_xy_yz[i] * a_exp + 4.0 * g_xxy_z_xy_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_xy_zz[i] = -2.0 * g_y_z_xy_zz[i] * a_exp + 4.0 * g_xxy_z_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_xx_0_0_0_y_z_xz_xx, g_xx_0_0_0_y_z_xz_xy, g_xx_0_0_0_y_z_xz_xz, g_xx_0_0_0_y_z_xz_yy, g_xx_0_0_0_y_z_xz_yz, g_xx_0_0_0_y_z_xz_zz, g_xxy_z_xz_xx, g_xxy_z_xz_xy, g_xxy_z_xz_xz, g_xxy_z_xz_yy, g_xxy_z_xz_yz, g_xxy_z_xz_zz, g_y_z_xz_xx, g_y_z_xz_xy, g_y_z_xz_xz, g_y_z_xz_yy, g_y_z_xz_yz, g_y_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_z_xz_xx[i] = -2.0 * g_y_z_xz_xx[i] * a_exp + 4.0 * g_xxy_z_xz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_xz_xy[i] = -2.0 * g_y_z_xz_xy[i] * a_exp + 4.0 * g_xxy_z_xz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_xz_xz[i] = -2.0 * g_y_z_xz_xz[i] * a_exp + 4.0 * g_xxy_z_xz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_xz_yy[i] = -2.0 * g_y_z_xz_yy[i] * a_exp + 4.0 * g_xxy_z_xz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_xz_yz[i] = -2.0 * g_y_z_xz_yz[i] * a_exp + 4.0 * g_xxy_z_xz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_xz_zz[i] = -2.0 * g_y_z_xz_zz[i] * a_exp + 4.0 * g_xxy_z_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_xx_0_0_0_y_z_yy_xx, g_xx_0_0_0_y_z_yy_xy, g_xx_0_0_0_y_z_yy_xz, g_xx_0_0_0_y_z_yy_yy, g_xx_0_0_0_y_z_yy_yz, g_xx_0_0_0_y_z_yy_zz, g_xxy_z_yy_xx, g_xxy_z_yy_xy, g_xxy_z_yy_xz, g_xxy_z_yy_yy, g_xxy_z_yy_yz, g_xxy_z_yy_zz, g_y_z_yy_xx, g_y_z_yy_xy, g_y_z_yy_xz, g_y_z_yy_yy, g_y_z_yy_yz, g_y_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_z_yy_xx[i] = -2.0 * g_y_z_yy_xx[i] * a_exp + 4.0 * g_xxy_z_yy_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_yy_xy[i] = -2.0 * g_y_z_yy_xy[i] * a_exp + 4.0 * g_xxy_z_yy_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_yy_xz[i] = -2.0 * g_y_z_yy_xz[i] * a_exp + 4.0 * g_xxy_z_yy_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_yy_yy[i] = -2.0 * g_y_z_yy_yy[i] * a_exp + 4.0 * g_xxy_z_yy_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_yy_yz[i] = -2.0 * g_y_z_yy_yz[i] * a_exp + 4.0 * g_xxy_z_yy_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_yy_zz[i] = -2.0 * g_y_z_yy_zz[i] * a_exp + 4.0 * g_xxy_z_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_xx_0_0_0_y_z_yz_xx, g_xx_0_0_0_y_z_yz_xy, g_xx_0_0_0_y_z_yz_xz, g_xx_0_0_0_y_z_yz_yy, g_xx_0_0_0_y_z_yz_yz, g_xx_0_0_0_y_z_yz_zz, g_xxy_z_yz_xx, g_xxy_z_yz_xy, g_xxy_z_yz_xz, g_xxy_z_yz_yy, g_xxy_z_yz_yz, g_xxy_z_yz_zz, g_y_z_yz_xx, g_y_z_yz_xy, g_y_z_yz_xz, g_y_z_yz_yy, g_y_z_yz_yz, g_y_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_z_yz_xx[i] = -2.0 * g_y_z_yz_xx[i] * a_exp + 4.0 * g_xxy_z_yz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_yz_xy[i] = -2.0 * g_y_z_yz_xy[i] * a_exp + 4.0 * g_xxy_z_yz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_yz_xz[i] = -2.0 * g_y_z_yz_xz[i] * a_exp + 4.0 * g_xxy_z_yz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_yz_yy[i] = -2.0 * g_y_z_yz_yy[i] * a_exp + 4.0 * g_xxy_z_yz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_yz_yz[i] = -2.0 * g_y_z_yz_yz[i] * a_exp + 4.0 * g_xxy_z_yz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_yz_zz[i] = -2.0 * g_y_z_yz_zz[i] * a_exp + 4.0 * g_xxy_z_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_xx_0_0_0_y_z_zz_xx, g_xx_0_0_0_y_z_zz_xy, g_xx_0_0_0_y_z_zz_xz, g_xx_0_0_0_y_z_zz_yy, g_xx_0_0_0_y_z_zz_yz, g_xx_0_0_0_y_z_zz_zz, g_xxy_z_zz_xx, g_xxy_z_zz_xy, g_xxy_z_zz_xz, g_xxy_z_zz_yy, g_xxy_z_zz_yz, g_xxy_z_zz_zz, g_y_z_zz_xx, g_y_z_zz_xy, g_y_z_zz_xz, g_y_z_zz_yy, g_y_z_zz_yz, g_y_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_z_zz_xx[i] = -2.0 * g_y_z_zz_xx[i] * a_exp + 4.0 * g_xxy_z_zz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_zz_xy[i] = -2.0 * g_y_z_zz_xy[i] * a_exp + 4.0 * g_xxy_z_zz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_zz_xz[i] = -2.0 * g_y_z_zz_xz[i] * a_exp + 4.0 * g_xxy_z_zz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_zz_yy[i] = -2.0 * g_y_z_zz_yy[i] * a_exp + 4.0 * g_xxy_z_zz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_zz_yz[i] = -2.0 * g_y_z_zz_yz[i] * a_exp + 4.0 * g_xxy_z_zz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_zz_zz[i] = -2.0 * g_y_z_zz_zz[i] * a_exp + 4.0 * g_xxy_z_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_xx_0_0_0_z_x_xx_xx, g_xx_0_0_0_z_x_xx_xy, g_xx_0_0_0_z_x_xx_xz, g_xx_0_0_0_z_x_xx_yy, g_xx_0_0_0_z_x_xx_yz, g_xx_0_0_0_z_x_xx_zz, g_xxz_x_xx_xx, g_xxz_x_xx_xy, g_xxz_x_xx_xz, g_xxz_x_xx_yy, g_xxz_x_xx_yz, g_xxz_x_xx_zz, g_z_x_xx_xx, g_z_x_xx_xy, g_z_x_xx_xz, g_z_x_xx_yy, g_z_x_xx_yz, g_z_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_x_xx_xx[i] = -2.0 * g_z_x_xx_xx[i] * a_exp + 4.0 * g_xxz_x_xx_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_xx_xy[i] = -2.0 * g_z_x_xx_xy[i] * a_exp + 4.0 * g_xxz_x_xx_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_xx_xz[i] = -2.0 * g_z_x_xx_xz[i] * a_exp + 4.0 * g_xxz_x_xx_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_xx_yy[i] = -2.0 * g_z_x_xx_yy[i] * a_exp + 4.0 * g_xxz_x_xx_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_xx_yz[i] = -2.0 * g_z_x_xx_yz[i] * a_exp + 4.0 * g_xxz_x_xx_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_xx_zz[i] = -2.0 * g_z_x_xx_zz[i] * a_exp + 4.0 * g_xxz_x_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_xx_0_0_0_z_x_xy_xx, g_xx_0_0_0_z_x_xy_xy, g_xx_0_0_0_z_x_xy_xz, g_xx_0_0_0_z_x_xy_yy, g_xx_0_0_0_z_x_xy_yz, g_xx_0_0_0_z_x_xy_zz, g_xxz_x_xy_xx, g_xxz_x_xy_xy, g_xxz_x_xy_xz, g_xxz_x_xy_yy, g_xxz_x_xy_yz, g_xxz_x_xy_zz, g_z_x_xy_xx, g_z_x_xy_xy, g_z_x_xy_xz, g_z_x_xy_yy, g_z_x_xy_yz, g_z_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_x_xy_xx[i] = -2.0 * g_z_x_xy_xx[i] * a_exp + 4.0 * g_xxz_x_xy_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_xy_xy[i] = -2.0 * g_z_x_xy_xy[i] * a_exp + 4.0 * g_xxz_x_xy_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_xy_xz[i] = -2.0 * g_z_x_xy_xz[i] * a_exp + 4.0 * g_xxz_x_xy_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_xy_yy[i] = -2.0 * g_z_x_xy_yy[i] * a_exp + 4.0 * g_xxz_x_xy_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_xy_yz[i] = -2.0 * g_z_x_xy_yz[i] * a_exp + 4.0 * g_xxz_x_xy_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_xy_zz[i] = -2.0 * g_z_x_xy_zz[i] * a_exp + 4.0 * g_xxz_x_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_xx_0_0_0_z_x_xz_xx, g_xx_0_0_0_z_x_xz_xy, g_xx_0_0_0_z_x_xz_xz, g_xx_0_0_0_z_x_xz_yy, g_xx_0_0_0_z_x_xz_yz, g_xx_0_0_0_z_x_xz_zz, g_xxz_x_xz_xx, g_xxz_x_xz_xy, g_xxz_x_xz_xz, g_xxz_x_xz_yy, g_xxz_x_xz_yz, g_xxz_x_xz_zz, g_z_x_xz_xx, g_z_x_xz_xy, g_z_x_xz_xz, g_z_x_xz_yy, g_z_x_xz_yz, g_z_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_x_xz_xx[i] = -2.0 * g_z_x_xz_xx[i] * a_exp + 4.0 * g_xxz_x_xz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_xz_xy[i] = -2.0 * g_z_x_xz_xy[i] * a_exp + 4.0 * g_xxz_x_xz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_xz_xz[i] = -2.0 * g_z_x_xz_xz[i] * a_exp + 4.0 * g_xxz_x_xz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_xz_yy[i] = -2.0 * g_z_x_xz_yy[i] * a_exp + 4.0 * g_xxz_x_xz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_xz_yz[i] = -2.0 * g_z_x_xz_yz[i] * a_exp + 4.0 * g_xxz_x_xz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_xz_zz[i] = -2.0 * g_z_x_xz_zz[i] * a_exp + 4.0 * g_xxz_x_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_xx_0_0_0_z_x_yy_xx, g_xx_0_0_0_z_x_yy_xy, g_xx_0_0_0_z_x_yy_xz, g_xx_0_0_0_z_x_yy_yy, g_xx_0_0_0_z_x_yy_yz, g_xx_0_0_0_z_x_yy_zz, g_xxz_x_yy_xx, g_xxz_x_yy_xy, g_xxz_x_yy_xz, g_xxz_x_yy_yy, g_xxz_x_yy_yz, g_xxz_x_yy_zz, g_z_x_yy_xx, g_z_x_yy_xy, g_z_x_yy_xz, g_z_x_yy_yy, g_z_x_yy_yz, g_z_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_x_yy_xx[i] = -2.0 * g_z_x_yy_xx[i] * a_exp + 4.0 * g_xxz_x_yy_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_yy_xy[i] = -2.0 * g_z_x_yy_xy[i] * a_exp + 4.0 * g_xxz_x_yy_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_yy_xz[i] = -2.0 * g_z_x_yy_xz[i] * a_exp + 4.0 * g_xxz_x_yy_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_yy_yy[i] = -2.0 * g_z_x_yy_yy[i] * a_exp + 4.0 * g_xxz_x_yy_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_yy_yz[i] = -2.0 * g_z_x_yy_yz[i] * a_exp + 4.0 * g_xxz_x_yy_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_yy_zz[i] = -2.0 * g_z_x_yy_zz[i] * a_exp + 4.0 * g_xxz_x_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_xx_0_0_0_z_x_yz_xx, g_xx_0_0_0_z_x_yz_xy, g_xx_0_0_0_z_x_yz_xz, g_xx_0_0_0_z_x_yz_yy, g_xx_0_0_0_z_x_yz_yz, g_xx_0_0_0_z_x_yz_zz, g_xxz_x_yz_xx, g_xxz_x_yz_xy, g_xxz_x_yz_xz, g_xxz_x_yz_yy, g_xxz_x_yz_yz, g_xxz_x_yz_zz, g_z_x_yz_xx, g_z_x_yz_xy, g_z_x_yz_xz, g_z_x_yz_yy, g_z_x_yz_yz, g_z_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_x_yz_xx[i] = -2.0 * g_z_x_yz_xx[i] * a_exp + 4.0 * g_xxz_x_yz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_yz_xy[i] = -2.0 * g_z_x_yz_xy[i] * a_exp + 4.0 * g_xxz_x_yz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_yz_xz[i] = -2.0 * g_z_x_yz_xz[i] * a_exp + 4.0 * g_xxz_x_yz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_yz_yy[i] = -2.0 * g_z_x_yz_yy[i] * a_exp + 4.0 * g_xxz_x_yz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_yz_yz[i] = -2.0 * g_z_x_yz_yz[i] * a_exp + 4.0 * g_xxz_x_yz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_yz_zz[i] = -2.0 * g_z_x_yz_zz[i] * a_exp + 4.0 * g_xxz_x_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_xx_0_0_0_z_x_zz_xx, g_xx_0_0_0_z_x_zz_xy, g_xx_0_0_0_z_x_zz_xz, g_xx_0_0_0_z_x_zz_yy, g_xx_0_0_0_z_x_zz_yz, g_xx_0_0_0_z_x_zz_zz, g_xxz_x_zz_xx, g_xxz_x_zz_xy, g_xxz_x_zz_xz, g_xxz_x_zz_yy, g_xxz_x_zz_yz, g_xxz_x_zz_zz, g_z_x_zz_xx, g_z_x_zz_xy, g_z_x_zz_xz, g_z_x_zz_yy, g_z_x_zz_yz, g_z_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_x_zz_xx[i] = -2.0 * g_z_x_zz_xx[i] * a_exp + 4.0 * g_xxz_x_zz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_zz_xy[i] = -2.0 * g_z_x_zz_xy[i] * a_exp + 4.0 * g_xxz_x_zz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_zz_xz[i] = -2.0 * g_z_x_zz_xz[i] * a_exp + 4.0 * g_xxz_x_zz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_zz_yy[i] = -2.0 * g_z_x_zz_yy[i] * a_exp + 4.0 * g_xxz_x_zz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_zz_yz[i] = -2.0 * g_z_x_zz_yz[i] * a_exp + 4.0 * g_xxz_x_zz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_zz_zz[i] = -2.0 * g_z_x_zz_zz[i] * a_exp + 4.0 * g_xxz_x_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_xx_0_0_0_z_y_xx_xx, g_xx_0_0_0_z_y_xx_xy, g_xx_0_0_0_z_y_xx_xz, g_xx_0_0_0_z_y_xx_yy, g_xx_0_0_0_z_y_xx_yz, g_xx_0_0_0_z_y_xx_zz, g_xxz_y_xx_xx, g_xxz_y_xx_xy, g_xxz_y_xx_xz, g_xxz_y_xx_yy, g_xxz_y_xx_yz, g_xxz_y_xx_zz, g_z_y_xx_xx, g_z_y_xx_xy, g_z_y_xx_xz, g_z_y_xx_yy, g_z_y_xx_yz, g_z_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_y_xx_xx[i] = -2.0 * g_z_y_xx_xx[i] * a_exp + 4.0 * g_xxz_y_xx_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_xx_xy[i] = -2.0 * g_z_y_xx_xy[i] * a_exp + 4.0 * g_xxz_y_xx_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_xx_xz[i] = -2.0 * g_z_y_xx_xz[i] * a_exp + 4.0 * g_xxz_y_xx_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_xx_yy[i] = -2.0 * g_z_y_xx_yy[i] * a_exp + 4.0 * g_xxz_y_xx_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_xx_yz[i] = -2.0 * g_z_y_xx_yz[i] * a_exp + 4.0 * g_xxz_y_xx_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_xx_zz[i] = -2.0 * g_z_y_xx_zz[i] * a_exp + 4.0 * g_xxz_y_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_xx_0_0_0_z_y_xy_xx, g_xx_0_0_0_z_y_xy_xy, g_xx_0_0_0_z_y_xy_xz, g_xx_0_0_0_z_y_xy_yy, g_xx_0_0_0_z_y_xy_yz, g_xx_0_0_0_z_y_xy_zz, g_xxz_y_xy_xx, g_xxz_y_xy_xy, g_xxz_y_xy_xz, g_xxz_y_xy_yy, g_xxz_y_xy_yz, g_xxz_y_xy_zz, g_z_y_xy_xx, g_z_y_xy_xy, g_z_y_xy_xz, g_z_y_xy_yy, g_z_y_xy_yz, g_z_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_y_xy_xx[i] = -2.0 * g_z_y_xy_xx[i] * a_exp + 4.0 * g_xxz_y_xy_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_xy_xy[i] = -2.0 * g_z_y_xy_xy[i] * a_exp + 4.0 * g_xxz_y_xy_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_xy_xz[i] = -2.0 * g_z_y_xy_xz[i] * a_exp + 4.0 * g_xxz_y_xy_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_xy_yy[i] = -2.0 * g_z_y_xy_yy[i] * a_exp + 4.0 * g_xxz_y_xy_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_xy_yz[i] = -2.0 * g_z_y_xy_yz[i] * a_exp + 4.0 * g_xxz_y_xy_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_xy_zz[i] = -2.0 * g_z_y_xy_zz[i] * a_exp + 4.0 * g_xxz_y_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_xx_0_0_0_z_y_xz_xx, g_xx_0_0_0_z_y_xz_xy, g_xx_0_0_0_z_y_xz_xz, g_xx_0_0_0_z_y_xz_yy, g_xx_0_0_0_z_y_xz_yz, g_xx_0_0_0_z_y_xz_zz, g_xxz_y_xz_xx, g_xxz_y_xz_xy, g_xxz_y_xz_xz, g_xxz_y_xz_yy, g_xxz_y_xz_yz, g_xxz_y_xz_zz, g_z_y_xz_xx, g_z_y_xz_xy, g_z_y_xz_xz, g_z_y_xz_yy, g_z_y_xz_yz, g_z_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_y_xz_xx[i] = -2.0 * g_z_y_xz_xx[i] * a_exp + 4.0 * g_xxz_y_xz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_xz_xy[i] = -2.0 * g_z_y_xz_xy[i] * a_exp + 4.0 * g_xxz_y_xz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_xz_xz[i] = -2.0 * g_z_y_xz_xz[i] * a_exp + 4.0 * g_xxz_y_xz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_xz_yy[i] = -2.0 * g_z_y_xz_yy[i] * a_exp + 4.0 * g_xxz_y_xz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_xz_yz[i] = -2.0 * g_z_y_xz_yz[i] * a_exp + 4.0 * g_xxz_y_xz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_xz_zz[i] = -2.0 * g_z_y_xz_zz[i] * a_exp + 4.0 * g_xxz_y_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_xx_0_0_0_z_y_yy_xx, g_xx_0_0_0_z_y_yy_xy, g_xx_0_0_0_z_y_yy_xz, g_xx_0_0_0_z_y_yy_yy, g_xx_0_0_0_z_y_yy_yz, g_xx_0_0_0_z_y_yy_zz, g_xxz_y_yy_xx, g_xxz_y_yy_xy, g_xxz_y_yy_xz, g_xxz_y_yy_yy, g_xxz_y_yy_yz, g_xxz_y_yy_zz, g_z_y_yy_xx, g_z_y_yy_xy, g_z_y_yy_xz, g_z_y_yy_yy, g_z_y_yy_yz, g_z_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_y_yy_xx[i] = -2.0 * g_z_y_yy_xx[i] * a_exp + 4.0 * g_xxz_y_yy_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_yy_xy[i] = -2.0 * g_z_y_yy_xy[i] * a_exp + 4.0 * g_xxz_y_yy_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_yy_xz[i] = -2.0 * g_z_y_yy_xz[i] * a_exp + 4.0 * g_xxz_y_yy_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_yy_yy[i] = -2.0 * g_z_y_yy_yy[i] * a_exp + 4.0 * g_xxz_y_yy_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_yy_yz[i] = -2.0 * g_z_y_yy_yz[i] * a_exp + 4.0 * g_xxz_y_yy_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_yy_zz[i] = -2.0 * g_z_y_yy_zz[i] * a_exp + 4.0 * g_xxz_y_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_xx_0_0_0_z_y_yz_xx, g_xx_0_0_0_z_y_yz_xy, g_xx_0_0_0_z_y_yz_xz, g_xx_0_0_0_z_y_yz_yy, g_xx_0_0_0_z_y_yz_yz, g_xx_0_0_0_z_y_yz_zz, g_xxz_y_yz_xx, g_xxz_y_yz_xy, g_xxz_y_yz_xz, g_xxz_y_yz_yy, g_xxz_y_yz_yz, g_xxz_y_yz_zz, g_z_y_yz_xx, g_z_y_yz_xy, g_z_y_yz_xz, g_z_y_yz_yy, g_z_y_yz_yz, g_z_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_y_yz_xx[i] = -2.0 * g_z_y_yz_xx[i] * a_exp + 4.0 * g_xxz_y_yz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_yz_xy[i] = -2.0 * g_z_y_yz_xy[i] * a_exp + 4.0 * g_xxz_y_yz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_yz_xz[i] = -2.0 * g_z_y_yz_xz[i] * a_exp + 4.0 * g_xxz_y_yz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_yz_yy[i] = -2.0 * g_z_y_yz_yy[i] * a_exp + 4.0 * g_xxz_y_yz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_yz_yz[i] = -2.0 * g_z_y_yz_yz[i] * a_exp + 4.0 * g_xxz_y_yz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_yz_zz[i] = -2.0 * g_z_y_yz_zz[i] * a_exp + 4.0 * g_xxz_y_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_xx_0_0_0_z_y_zz_xx, g_xx_0_0_0_z_y_zz_xy, g_xx_0_0_0_z_y_zz_xz, g_xx_0_0_0_z_y_zz_yy, g_xx_0_0_0_z_y_zz_yz, g_xx_0_0_0_z_y_zz_zz, g_xxz_y_zz_xx, g_xxz_y_zz_xy, g_xxz_y_zz_xz, g_xxz_y_zz_yy, g_xxz_y_zz_yz, g_xxz_y_zz_zz, g_z_y_zz_xx, g_z_y_zz_xy, g_z_y_zz_xz, g_z_y_zz_yy, g_z_y_zz_yz, g_z_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_y_zz_xx[i] = -2.0 * g_z_y_zz_xx[i] * a_exp + 4.0 * g_xxz_y_zz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_zz_xy[i] = -2.0 * g_z_y_zz_xy[i] * a_exp + 4.0 * g_xxz_y_zz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_zz_xz[i] = -2.0 * g_z_y_zz_xz[i] * a_exp + 4.0 * g_xxz_y_zz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_zz_yy[i] = -2.0 * g_z_y_zz_yy[i] * a_exp + 4.0 * g_xxz_y_zz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_zz_yz[i] = -2.0 * g_z_y_zz_yz[i] * a_exp + 4.0 * g_xxz_y_zz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_zz_zz[i] = -2.0 * g_z_y_zz_zz[i] * a_exp + 4.0 * g_xxz_y_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_xx_0_0_0_z_z_xx_xx, g_xx_0_0_0_z_z_xx_xy, g_xx_0_0_0_z_z_xx_xz, g_xx_0_0_0_z_z_xx_yy, g_xx_0_0_0_z_z_xx_yz, g_xx_0_0_0_z_z_xx_zz, g_xxz_z_xx_xx, g_xxz_z_xx_xy, g_xxz_z_xx_xz, g_xxz_z_xx_yy, g_xxz_z_xx_yz, g_xxz_z_xx_zz, g_z_z_xx_xx, g_z_z_xx_xy, g_z_z_xx_xz, g_z_z_xx_yy, g_z_z_xx_yz, g_z_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_z_xx_xx[i] = -2.0 * g_z_z_xx_xx[i] * a_exp + 4.0 * g_xxz_z_xx_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_xx_xy[i] = -2.0 * g_z_z_xx_xy[i] * a_exp + 4.0 * g_xxz_z_xx_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_xx_xz[i] = -2.0 * g_z_z_xx_xz[i] * a_exp + 4.0 * g_xxz_z_xx_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_xx_yy[i] = -2.0 * g_z_z_xx_yy[i] * a_exp + 4.0 * g_xxz_z_xx_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_xx_yz[i] = -2.0 * g_z_z_xx_yz[i] * a_exp + 4.0 * g_xxz_z_xx_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_xx_zz[i] = -2.0 * g_z_z_xx_zz[i] * a_exp + 4.0 * g_xxz_z_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_xx_0_0_0_z_z_xy_xx, g_xx_0_0_0_z_z_xy_xy, g_xx_0_0_0_z_z_xy_xz, g_xx_0_0_0_z_z_xy_yy, g_xx_0_0_0_z_z_xy_yz, g_xx_0_0_0_z_z_xy_zz, g_xxz_z_xy_xx, g_xxz_z_xy_xy, g_xxz_z_xy_xz, g_xxz_z_xy_yy, g_xxz_z_xy_yz, g_xxz_z_xy_zz, g_z_z_xy_xx, g_z_z_xy_xy, g_z_z_xy_xz, g_z_z_xy_yy, g_z_z_xy_yz, g_z_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_z_xy_xx[i] = -2.0 * g_z_z_xy_xx[i] * a_exp + 4.0 * g_xxz_z_xy_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_xy_xy[i] = -2.0 * g_z_z_xy_xy[i] * a_exp + 4.0 * g_xxz_z_xy_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_xy_xz[i] = -2.0 * g_z_z_xy_xz[i] * a_exp + 4.0 * g_xxz_z_xy_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_xy_yy[i] = -2.0 * g_z_z_xy_yy[i] * a_exp + 4.0 * g_xxz_z_xy_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_xy_yz[i] = -2.0 * g_z_z_xy_yz[i] * a_exp + 4.0 * g_xxz_z_xy_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_xy_zz[i] = -2.0 * g_z_z_xy_zz[i] * a_exp + 4.0 * g_xxz_z_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_xx_0_0_0_z_z_xz_xx, g_xx_0_0_0_z_z_xz_xy, g_xx_0_0_0_z_z_xz_xz, g_xx_0_0_0_z_z_xz_yy, g_xx_0_0_0_z_z_xz_yz, g_xx_0_0_0_z_z_xz_zz, g_xxz_z_xz_xx, g_xxz_z_xz_xy, g_xxz_z_xz_xz, g_xxz_z_xz_yy, g_xxz_z_xz_yz, g_xxz_z_xz_zz, g_z_z_xz_xx, g_z_z_xz_xy, g_z_z_xz_xz, g_z_z_xz_yy, g_z_z_xz_yz, g_z_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_z_xz_xx[i] = -2.0 * g_z_z_xz_xx[i] * a_exp + 4.0 * g_xxz_z_xz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_xz_xy[i] = -2.0 * g_z_z_xz_xy[i] * a_exp + 4.0 * g_xxz_z_xz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_xz_xz[i] = -2.0 * g_z_z_xz_xz[i] * a_exp + 4.0 * g_xxz_z_xz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_xz_yy[i] = -2.0 * g_z_z_xz_yy[i] * a_exp + 4.0 * g_xxz_z_xz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_xz_yz[i] = -2.0 * g_z_z_xz_yz[i] * a_exp + 4.0 * g_xxz_z_xz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_xz_zz[i] = -2.0 * g_z_z_xz_zz[i] * a_exp + 4.0 * g_xxz_z_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_xx_0_0_0_z_z_yy_xx, g_xx_0_0_0_z_z_yy_xy, g_xx_0_0_0_z_z_yy_xz, g_xx_0_0_0_z_z_yy_yy, g_xx_0_0_0_z_z_yy_yz, g_xx_0_0_0_z_z_yy_zz, g_xxz_z_yy_xx, g_xxz_z_yy_xy, g_xxz_z_yy_xz, g_xxz_z_yy_yy, g_xxz_z_yy_yz, g_xxz_z_yy_zz, g_z_z_yy_xx, g_z_z_yy_xy, g_z_z_yy_xz, g_z_z_yy_yy, g_z_z_yy_yz, g_z_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_z_yy_xx[i] = -2.0 * g_z_z_yy_xx[i] * a_exp + 4.0 * g_xxz_z_yy_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_yy_xy[i] = -2.0 * g_z_z_yy_xy[i] * a_exp + 4.0 * g_xxz_z_yy_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_yy_xz[i] = -2.0 * g_z_z_yy_xz[i] * a_exp + 4.0 * g_xxz_z_yy_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_yy_yy[i] = -2.0 * g_z_z_yy_yy[i] * a_exp + 4.0 * g_xxz_z_yy_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_yy_yz[i] = -2.0 * g_z_z_yy_yz[i] * a_exp + 4.0 * g_xxz_z_yy_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_yy_zz[i] = -2.0 * g_z_z_yy_zz[i] * a_exp + 4.0 * g_xxz_z_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_xx_0_0_0_z_z_yz_xx, g_xx_0_0_0_z_z_yz_xy, g_xx_0_0_0_z_z_yz_xz, g_xx_0_0_0_z_z_yz_yy, g_xx_0_0_0_z_z_yz_yz, g_xx_0_0_0_z_z_yz_zz, g_xxz_z_yz_xx, g_xxz_z_yz_xy, g_xxz_z_yz_xz, g_xxz_z_yz_yy, g_xxz_z_yz_yz, g_xxz_z_yz_zz, g_z_z_yz_xx, g_z_z_yz_xy, g_z_z_yz_xz, g_z_z_yz_yy, g_z_z_yz_yz, g_z_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_z_yz_xx[i] = -2.0 * g_z_z_yz_xx[i] * a_exp + 4.0 * g_xxz_z_yz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_yz_xy[i] = -2.0 * g_z_z_yz_xy[i] * a_exp + 4.0 * g_xxz_z_yz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_yz_xz[i] = -2.0 * g_z_z_yz_xz[i] * a_exp + 4.0 * g_xxz_z_yz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_yz_yy[i] = -2.0 * g_z_z_yz_yy[i] * a_exp + 4.0 * g_xxz_z_yz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_yz_yz[i] = -2.0 * g_z_z_yz_yz[i] * a_exp + 4.0 * g_xxz_z_yz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_yz_zz[i] = -2.0 * g_z_z_yz_zz[i] * a_exp + 4.0 * g_xxz_z_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_xx_0_0_0_z_z_zz_xx, g_xx_0_0_0_z_z_zz_xy, g_xx_0_0_0_z_z_zz_xz, g_xx_0_0_0_z_z_zz_yy, g_xx_0_0_0_z_z_zz_yz, g_xx_0_0_0_z_z_zz_zz, g_xxz_z_zz_xx, g_xxz_z_zz_xy, g_xxz_z_zz_xz, g_xxz_z_zz_yy, g_xxz_z_zz_yz, g_xxz_z_zz_zz, g_z_z_zz_xx, g_z_z_zz_xy, g_z_z_zz_xz, g_z_z_zz_yy, g_z_z_zz_yz, g_z_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_z_zz_xx[i] = -2.0 * g_z_z_zz_xx[i] * a_exp + 4.0 * g_xxz_z_zz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_zz_xy[i] = -2.0 * g_z_z_zz_xy[i] * a_exp + 4.0 * g_xxz_z_zz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_zz_xz[i] = -2.0 * g_z_z_zz_xz[i] * a_exp + 4.0 * g_xxz_z_zz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_zz_yy[i] = -2.0 * g_z_z_zz_yy[i] * a_exp + 4.0 * g_xxz_z_zz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_zz_yz[i] = -2.0 * g_z_z_zz_yz[i] * a_exp + 4.0 * g_xxz_z_zz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_zz_zz[i] = -2.0 * g_z_z_zz_zz[i] * a_exp + 4.0 * g_xxz_z_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (324-330)

    #pragma omp simd aligned(g_xxy_x_xx_xx, g_xxy_x_xx_xy, g_xxy_x_xx_xz, g_xxy_x_xx_yy, g_xxy_x_xx_yz, g_xxy_x_xx_zz, g_xy_0_0_0_x_x_xx_xx, g_xy_0_0_0_x_x_xx_xy, g_xy_0_0_0_x_x_xx_xz, g_xy_0_0_0_x_x_xx_yy, g_xy_0_0_0_x_x_xx_yz, g_xy_0_0_0_x_x_xx_zz, g_y_x_xx_xx, g_y_x_xx_xy, g_y_x_xx_xz, g_y_x_xx_yy, g_y_x_xx_yz, g_y_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_x_xx_xx[i] = -2.0 * g_y_x_xx_xx[i] * a_exp + 4.0 * g_xxy_x_xx_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_xx_xy[i] = -2.0 * g_y_x_xx_xy[i] * a_exp + 4.0 * g_xxy_x_xx_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_xx_xz[i] = -2.0 * g_y_x_xx_xz[i] * a_exp + 4.0 * g_xxy_x_xx_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_xx_yy[i] = -2.0 * g_y_x_xx_yy[i] * a_exp + 4.0 * g_xxy_x_xx_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_xx_yz[i] = -2.0 * g_y_x_xx_yz[i] * a_exp + 4.0 * g_xxy_x_xx_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_xx_zz[i] = -2.0 * g_y_x_xx_zz[i] * a_exp + 4.0 * g_xxy_x_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (330-336)

    #pragma omp simd aligned(g_xxy_x_xy_xx, g_xxy_x_xy_xy, g_xxy_x_xy_xz, g_xxy_x_xy_yy, g_xxy_x_xy_yz, g_xxy_x_xy_zz, g_xy_0_0_0_x_x_xy_xx, g_xy_0_0_0_x_x_xy_xy, g_xy_0_0_0_x_x_xy_xz, g_xy_0_0_0_x_x_xy_yy, g_xy_0_0_0_x_x_xy_yz, g_xy_0_0_0_x_x_xy_zz, g_y_x_xy_xx, g_y_x_xy_xy, g_y_x_xy_xz, g_y_x_xy_yy, g_y_x_xy_yz, g_y_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_x_xy_xx[i] = -2.0 * g_y_x_xy_xx[i] * a_exp + 4.0 * g_xxy_x_xy_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_xy_xy[i] = -2.0 * g_y_x_xy_xy[i] * a_exp + 4.0 * g_xxy_x_xy_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_xy_xz[i] = -2.0 * g_y_x_xy_xz[i] * a_exp + 4.0 * g_xxy_x_xy_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_xy_yy[i] = -2.0 * g_y_x_xy_yy[i] * a_exp + 4.0 * g_xxy_x_xy_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_xy_yz[i] = -2.0 * g_y_x_xy_yz[i] * a_exp + 4.0 * g_xxy_x_xy_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_xy_zz[i] = -2.0 * g_y_x_xy_zz[i] * a_exp + 4.0 * g_xxy_x_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (336-342)

    #pragma omp simd aligned(g_xxy_x_xz_xx, g_xxy_x_xz_xy, g_xxy_x_xz_xz, g_xxy_x_xz_yy, g_xxy_x_xz_yz, g_xxy_x_xz_zz, g_xy_0_0_0_x_x_xz_xx, g_xy_0_0_0_x_x_xz_xy, g_xy_0_0_0_x_x_xz_xz, g_xy_0_0_0_x_x_xz_yy, g_xy_0_0_0_x_x_xz_yz, g_xy_0_0_0_x_x_xz_zz, g_y_x_xz_xx, g_y_x_xz_xy, g_y_x_xz_xz, g_y_x_xz_yy, g_y_x_xz_yz, g_y_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_x_xz_xx[i] = -2.0 * g_y_x_xz_xx[i] * a_exp + 4.0 * g_xxy_x_xz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_xz_xy[i] = -2.0 * g_y_x_xz_xy[i] * a_exp + 4.0 * g_xxy_x_xz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_xz_xz[i] = -2.0 * g_y_x_xz_xz[i] * a_exp + 4.0 * g_xxy_x_xz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_xz_yy[i] = -2.0 * g_y_x_xz_yy[i] * a_exp + 4.0 * g_xxy_x_xz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_xz_yz[i] = -2.0 * g_y_x_xz_yz[i] * a_exp + 4.0 * g_xxy_x_xz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_xz_zz[i] = -2.0 * g_y_x_xz_zz[i] * a_exp + 4.0 * g_xxy_x_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (342-348)

    #pragma omp simd aligned(g_xxy_x_yy_xx, g_xxy_x_yy_xy, g_xxy_x_yy_xz, g_xxy_x_yy_yy, g_xxy_x_yy_yz, g_xxy_x_yy_zz, g_xy_0_0_0_x_x_yy_xx, g_xy_0_0_0_x_x_yy_xy, g_xy_0_0_0_x_x_yy_xz, g_xy_0_0_0_x_x_yy_yy, g_xy_0_0_0_x_x_yy_yz, g_xy_0_0_0_x_x_yy_zz, g_y_x_yy_xx, g_y_x_yy_xy, g_y_x_yy_xz, g_y_x_yy_yy, g_y_x_yy_yz, g_y_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_x_yy_xx[i] = -2.0 * g_y_x_yy_xx[i] * a_exp + 4.0 * g_xxy_x_yy_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_yy_xy[i] = -2.0 * g_y_x_yy_xy[i] * a_exp + 4.0 * g_xxy_x_yy_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_yy_xz[i] = -2.0 * g_y_x_yy_xz[i] * a_exp + 4.0 * g_xxy_x_yy_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_yy_yy[i] = -2.0 * g_y_x_yy_yy[i] * a_exp + 4.0 * g_xxy_x_yy_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_yy_yz[i] = -2.0 * g_y_x_yy_yz[i] * a_exp + 4.0 * g_xxy_x_yy_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_yy_zz[i] = -2.0 * g_y_x_yy_zz[i] * a_exp + 4.0 * g_xxy_x_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (348-354)

    #pragma omp simd aligned(g_xxy_x_yz_xx, g_xxy_x_yz_xy, g_xxy_x_yz_xz, g_xxy_x_yz_yy, g_xxy_x_yz_yz, g_xxy_x_yz_zz, g_xy_0_0_0_x_x_yz_xx, g_xy_0_0_0_x_x_yz_xy, g_xy_0_0_0_x_x_yz_xz, g_xy_0_0_0_x_x_yz_yy, g_xy_0_0_0_x_x_yz_yz, g_xy_0_0_0_x_x_yz_zz, g_y_x_yz_xx, g_y_x_yz_xy, g_y_x_yz_xz, g_y_x_yz_yy, g_y_x_yz_yz, g_y_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_x_yz_xx[i] = -2.0 * g_y_x_yz_xx[i] * a_exp + 4.0 * g_xxy_x_yz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_yz_xy[i] = -2.0 * g_y_x_yz_xy[i] * a_exp + 4.0 * g_xxy_x_yz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_yz_xz[i] = -2.0 * g_y_x_yz_xz[i] * a_exp + 4.0 * g_xxy_x_yz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_yz_yy[i] = -2.0 * g_y_x_yz_yy[i] * a_exp + 4.0 * g_xxy_x_yz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_yz_yz[i] = -2.0 * g_y_x_yz_yz[i] * a_exp + 4.0 * g_xxy_x_yz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_yz_zz[i] = -2.0 * g_y_x_yz_zz[i] * a_exp + 4.0 * g_xxy_x_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (354-360)

    #pragma omp simd aligned(g_xxy_x_zz_xx, g_xxy_x_zz_xy, g_xxy_x_zz_xz, g_xxy_x_zz_yy, g_xxy_x_zz_yz, g_xxy_x_zz_zz, g_xy_0_0_0_x_x_zz_xx, g_xy_0_0_0_x_x_zz_xy, g_xy_0_0_0_x_x_zz_xz, g_xy_0_0_0_x_x_zz_yy, g_xy_0_0_0_x_x_zz_yz, g_xy_0_0_0_x_x_zz_zz, g_y_x_zz_xx, g_y_x_zz_xy, g_y_x_zz_xz, g_y_x_zz_yy, g_y_x_zz_yz, g_y_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_x_zz_xx[i] = -2.0 * g_y_x_zz_xx[i] * a_exp + 4.0 * g_xxy_x_zz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_zz_xy[i] = -2.0 * g_y_x_zz_xy[i] * a_exp + 4.0 * g_xxy_x_zz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_zz_xz[i] = -2.0 * g_y_x_zz_xz[i] * a_exp + 4.0 * g_xxy_x_zz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_zz_yy[i] = -2.0 * g_y_x_zz_yy[i] * a_exp + 4.0 * g_xxy_x_zz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_zz_yz[i] = -2.0 * g_y_x_zz_yz[i] * a_exp + 4.0 * g_xxy_x_zz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_zz_zz[i] = -2.0 * g_y_x_zz_zz[i] * a_exp + 4.0 * g_xxy_x_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (360-366)

    #pragma omp simd aligned(g_xxy_y_xx_xx, g_xxy_y_xx_xy, g_xxy_y_xx_xz, g_xxy_y_xx_yy, g_xxy_y_xx_yz, g_xxy_y_xx_zz, g_xy_0_0_0_x_y_xx_xx, g_xy_0_0_0_x_y_xx_xy, g_xy_0_0_0_x_y_xx_xz, g_xy_0_0_0_x_y_xx_yy, g_xy_0_0_0_x_y_xx_yz, g_xy_0_0_0_x_y_xx_zz, g_y_y_xx_xx, g_y_y_xx_xy, g_y_y_xx_xz, g_y_y_xx_yy, g_y_y_xx_yz, g_y_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_y_xx_xx[i] = -2.0 * g_y_y_xx_xx[i] * a_exp + 4.0 * g_xxy_y_xx_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_xx_xy[i] = -2.0 * g_y_y_xx_xy[i] * a_exp + 4.0 * g_xxy_y_xx_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_xx_xz[i] = -2.0 * g_y_y_xx_xz[i] * a_exp + 4.0 * g_xxy_y_xx_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_xx_yy[i] = -2.0 * g_y_y_xx_yy[i] * a_exp + 4.0 * g_xxy_y_xx_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_xx_yz[i] = -2.0 * g_y_y_xx_yz[i] * a_exp + 4.0 * g_xxy_y_xx_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_xx_zz[i] = -2.0 * g_y_y_xx_zz[i] * a_exp + 4.0 * g_xxy_y_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (366-372)

    #pragma omp simd aligned(g_xxy_y_xy_xx, g_xxy_y_xy_xy, g_xxy_y_xy_xz, g_xxy_y_xy_yy, g_xxy_y_xy_yz, g_xxy_y_xy_zz, g_xy_0_0_0_x_y_xy_xx, g_xy_0_0_0_x_y_xy_xy, g_xy_0_0_0_x_y_xy_xz, g_xy_0_0_0_x_y_xy_yy, g_xy_0_0_0_x_y_xy_yz, g_xy_0_0_0_x_y_xy_zz, g_y_y_xy_xx, g_y_y_xy_xy, g_y_y_xy_xz, g_y_y_xy_yy, g_y_y_xy_yz, g_y_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_y_xy_xx[i] = -2.0 * g_y_y_xy_xx[i] * a_exp + 4.0 * g_xxy_y_xy_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_xy_xy[i] = -2.0 * g_y_y_xy_xy[i] * a_exp + 4.0 * g_xxy_y_xy_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_xy_xz[i] = -2.0 * g_y_y_xy_xz[i] * a_exp + 4.0 * g_xxy_y_xy_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_xy_yy[i] = -2.0 * g_y_y_xy_yy[i] * a_exp + 4.0 * g_xxy_y_xy_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_xy_yz[i] = -2.0 * g_y_y_xy_yz[i] * a_exp + 4.0 * g_xxy_y_xy_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_xy_zz[i] = -2.0 * g_y_y_xy_zz[i] * a_exp + 4.0 * g_xxy_y_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (372-378)

    #pragma omp simd aligned(g_xxy_y_xz_xx, g_xxy_y_xz_xy, g_xxy_y_xz_xz, g_xxy_y_xz_yy, g_xxy_y_xz_yz, g_xxy_y_xz_zz, g_xy_0_0_0_x_y_xz_xx, g_xy_0_0_0_x_y_xz_xy, g_xy_0_0_0_x_y_xz_xz, g_xy_0_0_0_x_y_xz_yy, g_xy_0_0_0_x_y_xz_yz, g_xy_0_0_0_x_y_xz_zz, g_y_y_xz_xx, g_y_y_xz_xy, g_y_y_xz_xz, g_y_y_xz_yy, g_y_y_xz_yz, g_y_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_y_xz_xx[i] = -2.0 * g_y_y_xz_xx[i] * a_exp + 4.0 * g_xxy_y_xz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_xz_xy[i] = -2.0 * g_y_y_xz_xy[i] * a_exp + 4.0 * g_xxy_y_xz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_xz_xz[i] = -2.0 * g_y_y_xz_xz[i] * a_exp + 4.0 * g_xxy_y_xz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_xz_yy[i] = -2.0 * g_y_y_xz_yy[i] * a_exp + 4.0 * g_xxy_y_xz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_xz_yz[i] = -2.0 * g_y_y_xz_yz[i] * a_exp + 4.0 * g_xxy_y_xz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_xz_zz[i] = -2.0 * g_y_y_xz_zz[i] * a_exp + 4.0 * g_xxy_y_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (378-384)

    #pragma omp simd aligned(g_xxy_y_yy_xx, g_xxy_y_yy_xy, g_xxy_y_yy_xz, g_xxy_y_yy_yy, g_xxy_y_yy_yz, g_xxy_y_yy_zz, g_xy_0_0_0_x_y_yy_xx, g_xy_0_0_0_x_y_yy_xy, g_xy_0_0_0_x_y_yy_xz, g_xy_0_0_0_x_y_yy_yy, g_xy_0_0_0_x_y_yy_yz, g_xy_0_0_0_x_y_yy_zz, g_y_y_yy_xx, g_y_y_yy_xy, g_y_y_yy_xz, g_y_y_yy_yy, g_y_y_yy_yz, g_y_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_y_yy_xx[i] = -2.0 * g_y_y_yy_xx[i] * a_exp + 4.0 * g_xxy_y_yy_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_yy_xy[i] = -2.0 * g_y_y_yy_xy[i] * a_exp + 4.0 * g_xxy_y_yy_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_yy_xz[i] = -2.0 * g_y_y_yy_xz[i] * a_exp + 4.0 * g_xxy_y_yy_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_yy_yy[i] = -2.0 * g_y_y_yy_yy[i] * a_exp + 4.0 * g_xxy_y_yy_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_yy_yz[i] = -2.0 * g_y_y_yy_yz[i] * a_exp + 4.0 * g_xxy_y_yy_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_yy_zz[i] = -2.0 * g_y_y_yy_zz[i] * a_exp + 4.0 * g_xxy_y_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (384-390)

    #pragma omp simd aligned(g_xxy_y_yz_xx, g_xxy_y_yz_xy, g_xxy_y_yz_xz, g_xxy_y_yz_yy, g_xxy_y_yz_yz, g_xxy_y_yz_zz, g_xy_0_0_0_x_y_yz_xx, g_xy_0_0_0_x_y_yz_xy, g_xy_0_0_0_x_y_yz_xz, g_xy_0_0_0_x_y_yz_yy, g_xy_0_0_0_x_y_yz_yz, g_xy_0_0_0_x_y_yz_zz, g_y_y_yz_xx, g_y_y_yz_xy, g_y_y_yz_xz, g_y_y_yz_yy, g_y_y_yz_yz, g_y_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_y_yz_xx[i] = -2.0 * g_y_y_yz_xx[i] * a_exp + 4.0 * g_xxy_y_yz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_yz_xy[i] = -2.0 * g_y_y_yz_xy[i] * a_exp + 4.0 * g_xxy_y_yz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_yz_xz[i] = -2.0 * g_y_y_yz_xz[i] * a_exp + 4.0 * g_xxy_y_yz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_yz_yy[i] = -2.0 * g_y_y_yz_yy[i] * a_exp + 4.0 * g_xxy_y_yz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_yz_yz[i] = -2.0 * g_y_y_yz_yz[i] * a_exp + 4.0 * g_xxy_y_yz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_yz_zz[i] = -2.0 * g_y_y_yz_zz[i] * a_exp + 4.0 * g_xxy_y_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (390-396)

    #pragma omp simd aligned(g_xxy_y_zz_xx, g_xxy_y_zz_xy, g_xxy_y_zz_xz, g_xxy_y_zz_yy, g_xxy_y_zz_yz, g_xxy_y_zz_zz, g_xy_0_0_0_x_y_zz_xx, g_xy_0_0_0_x_y_zz_xy, g_xy_0_0_0_x_y_zz_xz, g_xy_0_0_0_x_y_zz_yy, g_xy_0_0_0_x_y_zz_yz, g_xy_0_0_0_x_y_zz_zz, g_y_y_zz_xx, g_y_y_zz_xy, g_y_y_zz_xz, g_y_y_zz_yy, g_y_y_zz_yz, g_y_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_y_zz_xx[i] = -2.0 * g_y_y_zz_xx[i] * a_exp + 4.0 * g_xxy_y_zz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_zz_xy[i] = -2.0 * g_y_y_zz_xy[i] * a_exp + 4.0 * g_xxy_y_zz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_zz_xz[i] = -2.0 * g_y_y_zz_xz[i] * a_exp + 4.0 * g_xxy_y_zz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_zz_yy[i] = -2.0 * g_y_y_zz_yy[i] * a_exp + 4.0 * g_xxy_y_zz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_zz_yz[i] = -2.0 * g_y_y_zz_yz[i] * a_exp + 4.0 * g_xxy_y_zz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_zz_zz[i] = -2.0 * g_y_y_zz_zz[i] * a_exp + 4.0 * g_xxy_y_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (396-402)

    #pragma omp simd aligned(g_xxy_z_xx_xx, g_xxy_z_xx_xy, g_xxy_z_xx_xz, g_xxy_z_xx_yy, g_xxy_z_xx_yz, g_xxy_z_xx_zz, g_xy_0_0_0_x_z_xx_xx, g_xy_0_0_0_x_z_xx_xy, g_xy_0_0_0_x_z_xx_xz, g_xy_0_0_0_x_z_xx_yy, g_xy_0_0_0_x_z_xx_yz, g_xy_0_0_0_x_z_xx_zz, g_y_z_xx_xx, g_y_z_xx_xy, g_y_z_xx_xz, g_y_z_xx_yy, g_y_z_xx_yz, g_y_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_z_xx_xx[i] = -2.0 * g_y_z_xx_xx[i] * a_exp + 4.0 * g_xxy_z_xx_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_xx_xy[i] = -2.0 * g_y_z_xx_xy[i] * a_exp + 4.0 * g_xxy_z_xx_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_xx_xz[i] = -2.0 * g_y_z_xx_xz[i] * a_exp + 4.0 * g_xxy_z_xx_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_xx_yy[i] = -2.0 * g_y_z_xx_yy[i] * a_exp + 4.0 * g_xxy_z_xx_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_xx_yz[i] = -2.0 * g_y_z_xx_yz[i] * a_exp + 4.0 * g_xxy_z_xx_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_xx_zz[i] = -2.0 * g_y_z_xx_zz[i] * a_exp + 4.0 * g_xxy_z_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (402-408)

    #pragma omp simd aligned(g_xxy_z_xy_xx, g_xxy_z_xy_xy, g_xxy_z_xy_xz, g_xxy_z_xy_yy, g_xxy_z_xy_yz, g_xxy_z_xy_zz, g_xy_0_0_0_x_z_xy_xx, g_xy_0_0_0_x_z_xy_xy, g_xy_0_0_0_x_z_xy_xz, g_xy_0_0_0_x_z_xy_yy, g_xy_0_0_0_x_z_xy_yz, g_xy_0_0_0_x_z_xy_zz, g_y_z_xy_xx, g_y_z_xy_xy, g_y_z_xy_xz, g_y_z_xy_yy, g_y_z_xy_yz, g_y_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_z_xy_xx[i] = -2.0 * g_y_z_xy_xx[i] * a_exp + 4.0 * g_xxy_z_xy_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_xy_xy[i] = -2.0 * g_y_z_xy_xy[i] * a_exp + 4.0 * g_xxy_z_xy_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_xy_xz[i] = -2.0 * g_y_z_xy_xz[i] * a_exp + 4.0 * g_xxy_z_xy_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_xy_yy[i] = -2.0 * g_y_z_xy_yy[i] * a_exp + 4.0 * g_xxy_z_xy_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_xy_yz[i] = -2.0 * g_y_z_xy_yz[i] * a_exp + 4.0 * g_xxy_z_xy_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_xy_zz[i] = -2.0 * g_y_z_xy_zz[i] * a_exp + 4.0 * g_xxy_z_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (408-414)

    #pragma omp simd aligned(g_xxy_z_xz_xx, g_xxy_z_xz_xy, g_xxy_z_xz_xz, g_xxy_z_xz_yy, g_xxy_z_xz_yz, g_xxy_z_xz_zz, g_xy_0_0_0_x_z_xz_xx, g_xy_0_0_0_x_z_xz_xy, g_xy_0_0_0_x_z_xz_xz, g_xy_0_0_0_x_z_xz_yy, g_xy_0_0_0_x_z_xz_yz, g_xy_0_0_0_x_z_xz_zz, g_y_z_xz_xx, g_y_z_xz_xy, g_y_z_xz_xz, g_y_z_xz_yy, g_y_z_xz_yz, g_y_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_z_xz_xx[i] = -2.0 * g_y_z_xz_xx[i] * a_exp + 4.0 * g_xxy_z_xz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_xz_xy[i] = -2.0 * g_y_z_xz_xy[i] * a_exp + 4.0 * g_xxy_z_xz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_xz_xz[i] = -2.0 * g_y_z_xz_xz[i] * a_exp + 4.0 * g_xxy_z_xz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_xz_yy[i] = -2.0 * g_y_z_xz_yy[i] * a_exp + 4.0 * g_xxy_z_xz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_xz_yz[i] = -2.0 * g_y_z_xz_yz[i] * a_exp + 4.0 * g_xxy_z_xz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_xz_zz[i] = -2.0 * g_y_z_xz_zz[i] * a_exp + 4.0 * g_xxy_z_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (414-420)

    #pragma omp simd aligned(g_xxy_z_yy_xx, g_xxy_z_yy_xy, g_xxy_z_yy_xz, g_xxy_z_yy_yy, g_xxy_z_yy_yz, g_xxy_z_yy_zz, g_xy_0_0_0_x_z_yy_xx, g_xy_0_0_0_x_z_yy_xy, g_xy_0_0_0_x_z_yy_xz, g_xy_0_0_0_x_z_yy_yy, g_xy_0_0_0_x_z_yy_yz, g_xy_0_0_0_x_z_yy_zz, g_y_z_yy_xx, g_y_z_yy_xy, g_y_z_yy_xz, g_y_z_yy_yy, g_y_z_yy_yz, g_y_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_z_yy_xx[i] = -2.0 * g_y_z_yy_xx[i] * a_exp + 4.0 * g_xxy_z_yy_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_yy_xy[i] = -2.0 * g_y_z_yy_xy[i] * a_exp + 4.0 * g_xxy_z_yy_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_yy_xz[i] = -2.0 * g_y_z_yy_xz[i] * a_exp + 4.0 * g_xxy_z_yy_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_yy_yy[i] = -2.0 * g_y_z_yy_yy[i] * a_exp + 4.0 * g_xxy_z_yy_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_yy_yz[i] = -2.0 * g_y_z_yy_yz[i] * a_exp + 4.0 * g_xxy_z_yy_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_yy_zz[i] = -2.0 * g_y_z_yy_zz[i] * a_exp + 4.0 * g_xxy_z_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (420-426)

    #pragma omp simd aligned(g_xxy_z_yz_xx, g_xxy_z_yz_xy, g_xxy_z_yz_xz, g_xxy_z_yz_yy, g_xxy_z_yz_yz, g_xxy_z_yz_zz, g_xy_0_0_0_x_z_yz_xx, g_xy_0_0_0_x_z_yz_xy, g_xy_0_0_0_x_z_yz_xz, g_xy_0_0_0_x_z_yz_yy, g_xy_0_0_0_x_z_yz_yz, g_xy_0_0_0_x_z_yz_zz, g_y_z_yz_xx, g_y_z_yz_xy, g_y_z_yz_xz, g_y_z_yz_yy, g_y_z_yz_yz, g_y_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_z_yz_xx[i] = -2.0 * g_y_z_yz_xx[i] * a_exp + 4.0 * g_xxy_z_yz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_yz_xy[i] = -2.0 * g_y_z_yz_xy[i] * a_exp + 4.0 * g_xxy_z_yz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_yz_xz[i] = -2.0 * g_y_z_yz_xz[i] * a_exp + 4.0 * g_xxy_z_yz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_yz_yy[i] = -2.0 * g_y_z_yz_yy[i] * a_exp + 4.0 * g_xxy_z_yz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_yz_yz[i] = -2.0 * g_y_z_yz_yz[i] * a_exp + 4.0 * g_xxy_z_yz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_yz_zz[i] = -2.0 * g_y_z_yz_zz[i] * a_exp + 4.0 * g_xxy_z_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (426-432)

    #pragma omp simd aligned(g_xxy_z_zz_xx, g_xxy_z_zz_xy, g_xxy_z_zz_xz, g_xxy_z_zz_yy, g_xxy_z_zz_yz, g_xxy_z_zz_zz, g_xy_0_0_0_x_z_zz_xx, g_xy_0_0_0_x_z_zz_xy, g_xy_0_0_0_x_z_zz_xz, g_xy_0_0_0_x_z_zz_yy, g_xy_0_0_0_x_z_zz_yz, g_xy_0_0_0_x_z_zz_zz, g_y_z_zz_xx, g_y_z_zz_xy, g_y_z_zz_xz, g_y_z_zz_yy, g_y_z_zz_yz, g_y_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_z_zz_xx[i] = -2.0 * g_y_z_zz_xx[i] * a_exp + 4.0 * g_xxy_z_zz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_zz_xy[i] = -2.0 * g_y_z_zz_xy[i] * a_exp + 4.0 * g_xxy_z_zz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_zz_xz[i] = -2.0 * g_y_z_zz_xz[i] * a_exp + 4.0 * g_xxy_z_zz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_zz_yy[i] = -2.0 * g_y_z_zz_yy[i] * a_exp + 4.0 * g_xxy_z_zz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_zz_yz[i] = -2.0 * g_y_z_zz_yz[i] * a_exp + 4.0 * g_xxy_z_zz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_zz_zz[i] = -2.0 * g_y_z_zz_zz[i] * a_exp + 4.0 * g_xxy_z_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (432-438)

    #pragma omp simd aligned(g_x_x_xx_xx, g_x_x_xx_xy, g_x_x_xx_xz, g_x_x_xx_yy, g_x_x_xx_yz, g_x_x_xx_zz, g_xy_0_0_0_y_x_xx_xx, g_xy_0_0_0_y_x_xx_xy, g_xy_0_0_0_y_x_xx_xz, g_xy_0_0_0_y_x_xx_yy, g_xy_0_0_0_y_x_xx_yz, g_xy_0_0_0_y_x_xx_zz, g_xyy_x_xx_xx, g_xyy_x_xx_xy, g_xyy_x_xx_xz, g_xyy_x_xx_yy, g_xyy_x_xx_yz, g_xyy_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_x_xx_xx[i] = -2.0 * g_x_x_xx_xx[i] * a_exp + 4.0 * g_xyy_x_xx_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_xx_xy[i] = -2.0 * g_x_x_xx_xy[i] * a_exp + 4.0 * g_xyy_x_xx_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_xx_xz[i] = -2.0 * g_x_x_xx_xz[i] * a_exp + 4.0 * g_xyy_x_xx_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_xx_yy[i] = -2.0 * g_x_x_xx_yy[i] * a_exp + 4.0 * g_xyy_x_xx_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_xx_yz[i] = -2.0 * g_x_x_xx_yz[i] * a_exp + 4.0 * g_xyy_x_xx_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_xx_zz[i] = -2.0 * g_x_x_xx_zz[i] * a_exp + 4.0 * g_xyy_x_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (438-444)

    #pragma omp simd aligned(g_x_x_xy_xx, g_x_x_xy_xy, g_x_x_xy_xz, g_x_x_xy_yy, g_x_x_xy_yz, g_x_x_xy_zz, g_xy_0_0_0_y_x_xy_xx, g_xy_0_0_0_y_x_xy_xy, g_xy_0_0_0_y_x_xy_xz, g_xy_0_0_0_y_x_xy_yy, g_xy_0_0_0_y_x_xy_yz, g_xy_0_0_0_y_x_xy_zz, g_xyy_x_xy_xx, g_xyy_x_xy_xy, g_xyy_x_xy_xz, g_xyy_x_xy_yy, g_xyy_x_xy_yz, g_xyy_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_x_xy_xx[i] = -2.0 * g_x_x_xy_xx[i] * a_exp + 4.0 * g_xyy_x_xy_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_xy_xy[i] = -2.0 * g_x_x_xy_xy[i] * a_exp + 4.0 * g_xyy_x_xy_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_xy_xz[i] = -2.0 * g_x_x_xy_xz[i] * a_exp + 4.0 * g_xyy_x_xy_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_xy_yy[i] = -2.0 * g_x_x_xy_yy[i] * a_exp + 4.0 * g_xyy_x_xy_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_xy_yz[i] = -2.0 * g_x_x_xy_yz[i] * a_exp + 4.0 * g_xyy_x_xy_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_xy_zz[i] = -2.0 * g_x_x_xy_zz[i] * a_exp + 4.0 * g_xyy_x_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (444-450)

    #pragma omp simd aligned(g_x_x_xz_xx, g_x_x_xz_xy, g_x_x_xz_xz, g_x_x_xz_yy, g_x_x_xz_yz, g_x_x_xz_zz, g_xy_0_0_0_y_x_xz_xx, g_xy_0_0_0_y_x_xz_xy, g_xy_0_0_0_y_x_xz_xz, g_xy_0_0_0_y_x_xz_yy, g_xy_0_0_0_y_x_xz_yz, g_xy_0_0_0_y_x_xz_zz, g_xyy_x_xz_xx, g_xyy_x_xz_xy, g_xyy_x_xz_xz, g_xyy_x_xz_yy, g_xyy_x_xz_yz, g_xyy_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_x_xz_xx[i] = -2.0 * g_x_x_xz_xx[i] * a_exp + 4.0 * g_xyy_x_xz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_xz_xy[i] = -2.0 * g_x_x_xz_xy[i] * a_exp + 4.0 * g_xyy_x_xz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_xz_xz[i] = -2.0 * g_x_x_xz_xz[i] * a_exp + 4.0 * g_xyy_x_xz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_xz_yy[i] = -2.0 * g_x_x_xz_yy[i] * a_exp + 4.0 * g_xyy_x_xz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_xz_yz[i] = -2.0 * g_x_x_xz_yz[i] * a_exp + 4.0 * g_xyy_x_xz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_xz_zz[i] = -2.0 * g_x_x_xz_zz[i] * a_exp + 4.0 * g_xyy_x_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (450-456)

    #pragma omp simd aligned(g_x_x_yy_xx, g_x_x_yy_xy, g_x_x_yy_xz, g_x_x_yy_yy, g_x_x_yy_yz, g_x_x_yy_zz, g_xy_0_0_0_y_x_yy_xx, g_xy_0_0_0_y_x_yy_xy, g_xy_0_0_0_y_x_yy_xz, g_xy_0_0_0_y_x_yy_yy, g_xy_0_0_0_y_x_yy_yz, g_xy_0_0_0_y_x_yy_zz, g_xyy_x_yy_xx, g_xyy_x_yy_xy, g_xyy_x_yy_xz, g_xyy_x_yy_yy, g_xyy_x_yy_yz, g_xyy_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_x_yy_xx[i] = -2.0 * g_x_x_yy_xx[i] * a_exp + 4.0 * g_xyy_x_yy_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_yy_xy[i] = -2.0 * g_x_x_yy_xy[i] * a_exp + 4.0 * g_xyy_x_yy_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_yy_xz[i] = -2.0 * g_x_x_yy_xz[i] * a_exp + 4.0 * g_xyy_x_yy_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_yy_yy[i] = -2.0 * g_x_x_yy_yy[i] * a_exp + 4.0 * g_xyy_x_yy_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_yy_yz[i] = -2.0 * g_x_x_yy_yz[i] * a_exp + 4.0 * g_xyy_x_yy_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_yy_zz[i] = -2.0 * g_x_x_yy_zz[i] * a_exp + 4.0 * g_xyy_x_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (456-462)

    #pragma omp simd aligned(g_x_x_yz_xx, g_x_x_yz_xy, g_x_x_yz_xz, g_x_x_yz_yy, g_x_x_yz_yz, g_x_x_yz_zz, g_xy_0_0_0_y_x_yz_xx, g_xy_0_0_0_y_x_yz_xy, g_xy_0_0_0_y_x_yz_xz, g_xy_0_0_0_y_x_yz_yy, g_xy_0_0_0_y_x_yz_yz, g_xy_0_0_0_y_x_yz_zz, g_xyy_x_yz_xx, g_xyy_x_yz_xy, g_xyy_x_yz_xz, g_xyy_x_yz_yy, g_xyy_x_yz_yz, g_xyy_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_x_yz_xx[i] = -2.0 * g_x_x_yz_xx[i] * a_exp + 4.0 * g_xyy_x_yz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_yz_xy[i] = -2.0 * g_x_x_yz_xy[i] * a_exp + 4.0 * g_xyy_x_yz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_yz_xz[i] = -2.0 * g_x_x_yz_xz[i] * a_exp + 4.0 * g_xyy_x_yz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_yz_yy[i] = -2.0 * g_x_x_yz_yy[i] * a_exp + 4.0 * g_xyy_x_yz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_yz_yz[i] = -2.0 * g_x_x_yz_yz[i] * a_exp + 4.0 * g_xyy_x_yz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_yz_zz[i] = -2.0 * g_x_x_yz_zz[i] * a_exp + 4.0 * g_xyy_x_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (462-468)

    #pragma omp simd aligned(g_x_x_zz_xx, g_x_x_zz_xy, g_x_x_zz_xz, g_x_x_zz_yy, g_x_x_zz_yz, g_x_x_zz_zz, g_xy_0_0_0_y_x_zz_xx, g_xy_0_0_0_y_x_zz_xy, g_xy_0_0_0_y_x_zz_xz, g_xy_0_0_0_y_x_zz_yy, g_xy_0_0_0_y_x_zz_yz, g_xy_0_0_0_y_x_zz_zz, g_xyy_x_zz_xx, g_xyy_x_zz_xy, g_xyy_x_zz_xz, g_xyy_x_zz_yy, g_xyy_x_zz_yz, g_xyy_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_x_zz_xx[i] = -2.0 * g_x_x_zz_xx[i] * a_exp + 4.0 * g_xyy_x_zz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_zz_xy[i] = -2.0 * g_x_x_zz_xy[i] * a_exp + 4.0 * g_xyy_x_zz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_zz_xz[i] = -2.0 * g_x_x_zz_xz[i] * a_exp + 4.0 * g_xyy_x_zz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_zz_yy[i] = -2.0 * g_x_x_zz_yy[i] * a_exp + 4.0 * g_xyy_x_zz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_zz_yz[i] = -2.0 * g_x_x_zz_yz[i] * a_exp + 4.0 * g_xyy_x_zz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_zz_zz[i] = -2.0 * g_x_x_zz_zz[i] * a_exp + 4.0 * g_xyy_x_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (468-474)

    #pragma omp simd aligned(g_x_y_xx_xx, g_x_y_xx_xy, g_x_y_xx_xz, g_x_y_xx_yy, g_x_y_xx_yz, g_x_y_xx_zz, g_xy_0_0_0_y_y_xx_xx, g_xy_0_0_0_y_y_xx_xy, g_xy_0_0_0_y_y_xx_xz, g_xy_0_0_0_y_y_xx_yy, g_xy_0_0_0_y_y_xx_yz, g_xy_0_0_0_y_y_xx_zz, g_xyy_y_xx_xx, g_xyy_y_xx_xy, g_xyy_y_xx_xz, g_xyy_y_xx_yy, g_xyy_y_xx_yz, g_xyy_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_y_xx_xx[i] = -2.0 * g_x_y_xx_xx[i] * a_exp + 4.0 * g_xyy_y_xx_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_xx_xy[i] = -2.0 * g_x_y_xx_xy[i] * a_exp + 4.0 * g_xyy_y_xx_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_xx_xz[i] = -2.0 * g_x_y_xx_xz[i] * a_exp + 4.0 * g_xyy_y_xx_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_xx_yy[i] = -2.0 * g_x_y_xx_yy[i] * a_exp + 4.0 * g_xyy_y_xx_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_xx_yz[i] = -2.0 * g_x_y_xx_yz[i] * a_exp + 4.0 * g_xyy_y_xx_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_xx_zz[i] = -2.0 * g_x_y_xx_zz[i] * a_exp + 4.0 * g_xyy_y_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (474-480)

    #pragma omp simd aligned(g_x_y_xy_xx, g_x_y_xy_xy, g_x_y_xy_xz, g_x_y_xy_yy, g_x_y_xy_yz, g_x_y_xy_zz, g_xy_0_0_0_y_y_xy_xx, g_xy_0_0_0_y_y_xy_xy, g_xy_0_0_0_y_y_xy_xz, g_xy_0_0_0_y_y_xy_yy, g_xy_0_0_0_y_y_xy_yz, g_xy_0_0_0_y_y_xy_zz, g_xyy_y_xy_xx, g_xyy_y_xy_xy, g_xyy_y_xy_xz, g_xyy_y_xy_yy, g_xyy_y_xy_yz, g_xyy_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_y_xy_xx[i] = -2.0 * g_x_y_xy_xx[i] * a_exp + 4.0 * g_xyy_y_xy_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_xy_xy[i] = -2.0 * g_x_y_xy_xy[i] * a_exp + 4.0 * g_xyy_y_xy_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_xy_xz[i] = -2.0 * g_x_y_xy_xz[i] * a_exp + 4.0 * g_xyy_y_xy_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_xy_yy[i] = -2.0 * g_x_y_xy_yy[i] * a_exp + 4.0 * g_xyy_y_xy_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_xy_yz[i] = -2.0 * g_x_y_xy_yz[i] * a_exp + 4.0 * g_xyy_y_xy_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_xy_zz[i] = -2.0 * g_x_y_xy_zz[i] * a_exp + 4.0 * g_xyy_y_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (480-486)

    #pragma omp simd aligned(g_x_y_xz_xx, g_x_y_xz_xy, g_x_y_xz_xz, g_x_y_xz_yy, g_x_y_xz_yz, g_x_y_xz_zz, g_xy_0_0_0_y_y_xz_xx, g_xy_0_0_0_y_y_xz_xy, g_xy_0_0_0_y_y_xz_xz, g_xy_0_0_0_y_y_xz_yy, g_xy_0_0_0_y_y_xz_yz, g_xy_0_0_0_y_y_xz_zz, g_xyy_y_xz_xx, g_xyy_y_xz_xy, g_xyy_y_xz_xz, g_xyy_y_xz_yy, g_xyy_y_xz_yz, g_xyy_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_y_xz_xx[i] = -2.0 * g_x_y_xz_xx[i] * a_exp + 4.0 * g_xyy_y_xz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_xz_xy[i] = -2.0 * g_x_y_xz_xy[i] * a_exp + 4.0 * g_xyy_y_xz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_xz_xz[i] = -2.0 * g_x_y_xz_xz[i] * a_exp + 4.0 * g_xyy_y_xz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_xz_yy[i] = -2.0 * g_x_y_xz_yy[i] * a_exp + 4.0 * g_xyy_y_xz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_xz_yz[i] = -2.0 * g_x_y_xz_yz[i] * a_exp + 4.0 * g_xyy_y_xz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_xz_zz[i] = -2.0 * g_x_y_xz_zz[i] * a_exp + 4.0 * g_xyy_y_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (486-492)

    #pragma omp simd aligned(g_x_y_yy_xx, g_x_y_yy_xy, g_x_y_yy_xz, g_x_y_yy_yy, g_x_y_yy_yz, g_x_y_yy_zz, g_xy_0_0_0_y_y_yy_xx, g_xy_0_0_0_y_y_yy_xy, g_xy_0_0_0_y_y_yy_xz, g_xy_0_0_0_y_y_yy_yy, g_xy_0_0_0_y_y_yy_yz, g_xy_0_0_0_y_y_yy_zz, g_xyy_y_yy_xx, g_xyy_y_yy_xy, g_xyy_y_yy_xz, g_xyy_y_yy_yy, g_xyy_y_yy_yz, g_xyy_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_y_yy_xx[i] = -2.0 * g_x_y_yy_xx[i] * a_exp + 4.0 * g_xyy_y_yy_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_yy_xy[i] = -2.0 * g_x_y_yy_xy[i] * a_exp + 4.0 * g_xyy_y_yy_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_yy_xz[i] = -2.0 * g_x_y_yy_xz[i] * a_exp + 4.0 * g_xyy_y_yy_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_yy_yy[i] = -2.0 * g_x_y_yy_yy[i] * a_exp + 4.0 * g_xyy_y_yy_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_yy_yz[i] = -2.0 * g_x_y_yy_yz[i] * a_exp + 4.0 * g_xyy_y_yy_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_yy_zz[i] = -2.0 * g_x_y_yy_zz[i] * a_exp + 4.0 * g_xyy_y_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (492-498)

    #pragma omp simd aligned(g_x_y_yz_xx, g_x_y_yz_xy, g_x_y_yz_xz, g_x_y_yz_yy, g_x_y_yz_yz, g_x_y_yz_zz, g_xy_0_0_0_y_y_yz_xx, g_xy_0_0_0_y_y_yz_xy, g_xy_0_0_0_y_y_yz_xz, g_xy_0_0_0_y_y_yz_yy, g_xy_0_0_0_y_y_yz_yz, g_xy_0_0_0_y_y_yz_zz, g_xyy_y_yz_xx, g_xyy_y_yz_xy, g_xyy_y_yz_xz, g_xyy_y_yz_yy, g_xyy_y_yz_yz, g_xyy_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_y_yz_xx[i] = -2.0 * g_x_y_yz_xx[i] * a_exp + 4.0 * g_xyy_y_yz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_yz_xy[i] = -2.0 * g_x_y_yz_xy[i] * a_exp + 4.0 * g_xyy_y_yz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_yz_xz[i] = -2.0 * g_x_y_yz_xz[i] * a_exp + 4.0 * g_xyy_y_yz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_yz_yy[i] = -2.0 * g_x_y_yz_yy[i] * a_exp + 4.0 * g_xyy_y_yz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_yz_yz[i] = -2.0 * g_x_y_yz_yz[i] * a_exp + 4.0 * g_xyy_y_yz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_yz_zz[i] = -2.0 * g_x_y_yz_zz[i] * a_exp + 4.0 * g_xyy_y_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (498-504)

    #pragma omp simd aligned(g_x_y_zz_xx, g_x_y_zz_xy, g_x_y_zz_xz, g_x_y_zz_yy, g_x_y_zz_yz, g_x_y_zz_zz, g_xy_0_0_0_y_y_zz_xx, g_xy_0_0_0_y_y_zz_xy, g_xy_0_0_0_y_y_zz_xz, g_xy_0_0_0_y_y_zz_yy, g_xy_0_0_0_y_y_zz_yz, g_xy_0_0_0_y_y_zz_zz, g_xyy_y_zz_xx, g_xyy_y_zz_xy, g_xyy_y_zz_xz, g_xyy_y_zz_yy, g_xyy_y_zz_yz, g_xyy_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_y_zz_xx[i] = -2.0 * g_x_y_zz_xx[i] * a_exp + 4.0 * g_xyy_y_zz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_zz_xy[i] = -2.0 * g_x_y_zz_xy[i] * a_exp + 4.0 * g_xyy_y_zz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_zz_xz[i] = -2.0 * g_x_y_zz_xz[i] * a_exp + 4.0 * g_xyy_y_zz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_zz_yy[i] = -2.0 * g_x_y_zz_yy[i] * a_exp + 4.0 * g_xyy_y_zz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_zz_yz[i] = -2.0 * g_x_y_zz_yz[i] * a_exp + 4.0 * g_xyy_y_zz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_zz_zz[i] = -2.0 * g_x_y_zz_zz[i] * a_exp + 4.0 * g_xyy_y_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (504-510)

    #pragma omp simd aligned(g_x_z_xx_xx, g_x_z_xx_xy, g_x_z_xx_xz, g_x_z_xx_yy, g_x_z_xx_yz, g_x_z_xx_zz, g_xy_0_0_0_y_z_xx_xx, g_xy_0_0_0_y_z_xx_xy, g_xy_0_0_0_y_z_xx_xz, g_xy_0_0_0_y_z_xx_yy, g_xy_0_0_0_y_z_xx_yz, g_xy_0_0_0_y_z_xx_zz, g_xyy_z_xx_xx, g_xyy_z_xx_xy, g_xyy_z_xx_xz, g_xyy_z_xx_yy, g_xyy_z_xx_yz, g_xyy_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_z_xx_xx[i] = -2.0 * g_x_z_xx_xx[i] * a_exp + 4.0 * g_xyy_z_xx_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_xx_xy[i] = -2.0 * g_x_z_xx_xy[i] * a_exp + 4.0 * g_xyy_z_xx_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_xx_xz[i] = -2.0 * g_x_z_xx_xz[i] * a_exp + 4.0 * g_xyy_z_xx_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_xx_yy[i] = -2.0 * g_x_z_xx_yy[i] * a_exp + 4.0 * g_xyy_z_xx_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_xx_yz[i] = -2.0 * g_x_z_xx_yz[i] * a_exp + 4.0 * g_xyy_z_xx_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_xx_zz[i] = -2.0 * g_x_z_xx_zz[i] * a_exp + 4.0 * g_xyy_z_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (510-516)

    #pragma omp simd aligned(g_x_z_xy_xx, g_x_z_xy_xy, g_x_z_xy_xz, g_x_z_xy_yy, g_x_z_xy_yz, g_x_z_xy_zz, g_xy_0_0_0_y_z_xy_xx, g_xy_0_0_0_y_z_xy_xy, g_xy_0_0_0_y_z_xy_xz, g_xy_0_0_0_y_z_xy_yy, g_xy_0_0_0_y_z_xy_yz, g_xy_0_0_0_y_z_xy_zz, g_xyy_z_xy_xx, g_xyy_z_xy_xy, g_xyy_z_xy_xz, g_xyy_z_xy_yy, g_xyy_z_xy_yz, g_xyy_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_z_xy_xx[i] = -2.0 * g_x_z_xy_xx[i] * a_exp + 4.0 * g_xyy_z_xy_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_xy_xy[i] = -2.0 * g_x_z_xy_xy[i] * a_exp + 4.0 * g_xyy_z_xy_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_xy_xz[i] = -2.0 * g_x_z_xy_xz[i] * a_exp + 4.0 * g_xyy_z_xy_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_xy_yy[i] = -2.0 * g_x_z_xy_yy[i] * a_exp + 4.0 * g_xyy_z_xy_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_xy_yz[i] = -2.0 * g_x_z_xy_yz[i] * a_exp + 4.0 * g_xyy_z_xy_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_xy_zz[i] = -2.0 * g_x_z_xy_zz[i] * a_exp + 4.0 * g_xyy_z_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (516-522)

    #pragma omp simd aligned(g_x_z_xz_xx, g_x_z_xz_xy, g_x_z_xz_xz, g_x_z_xz_yy, g_x_z_xz_yz, g_x_z_xz_zz, g_xy_0_0_0_y_z_xz_xx, g_xy_0_0_0_y_z_xz_xy, g_xy_0_0_0_y_z_xz_xz, g_xy_0_0_0_y_z_xz_yy, g_xy_0_0_0_y_z_xz_yz, g_xy_0_0_0_y_z_xz_zz, g_xyy_z_xz_xx, g_xyy_z_xz_xy, g_xyy_z_xz_xz, g_xyy_z_xz_yy, g_xyy_z_xz_yz, g_xyy_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_z_xz_xx[i] = -2.0 * g_x_z_xz_xx[i] * a_exp + 4.0 * g_xyy_z_xz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_xz_xy[i] = -2.0 * g_x_z_xz_xy[i] * a_exp + 4.0 * g_xyy_z_xz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_xz_xz[i] = -2.0 * g_x_z_xz_xz[i] * a_exp + 4.0 * g_xyy_z_xz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_xz_yy[i] = -2.0 * g_x_z_xz_yy[i] * a_exp + 4.0 * g_xyy_z_xz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_xz_yz[i] = -2.0 * g_x_z_xz_yz[i] * a_exp + 4.0 * g_xyy_z_xz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_xz_zz[i] = -2.0 * g_x_z_xz_zz[i] * a_exp + 4.0 * g_xyy_z_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (522-528)

    #pragma omp simd aligned(g_x_z_yy_xx, g_x_z_yy_xy, g_x_z_yy_xz, g_x_z_yy_yy, g_x_z_yy_yz, g_x_z_yy_zz, g_xy_0_0_0_y_z_yy_xx, g_xy_0_0_0_y_z_yy_xy, g_xy_0_0_0_y_z_yy_xz, g_xy_0_0_0_y_z_yy_yy, g_xy_0_0_0_y_z_yy_yz, g_xy_0_0_0_y_z_yy_zz, g_xyy_z_yy_xx, g_xyy_z_yy_xy, g_xyy_z_yy_xz, g_xyy_z_yy_yy, g_xyy_z_yy_yz, g_xyy_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_z_yy_xx[i] = -2.0 * g_x_z_yy_xx[i] * a_exp + 4.0 * g_xyy_z_yy_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_yy_xy[i] = -2.0 * g_x_z_yy_xy[i] * a_exp + 4.0 * g_xyy_z_yy_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_yy_xz[i] = -2.0 * g_x_z_yy_xz[i] * a_exp + 4.0 * g_xyy_z_yy_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_yy_yy[i] = -2.0 * g_x_z_yy_yy[i] * a_exp + 4.0 * g_xyy_z_yy_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_yy_yz[i] = -2.0 * g_x_z_yy_yz[i] * a_exp + 4.0 * g_xyy_z_yy_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_yy_zz[i] = -2.0 * g_x_z_yy_zz[i] * a_exp + 4.0 * g_xyy_z_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (528-534)

    #pragma omp simd aligned(g_x_z_yz_xx, g_x_z_yz_xy, g_x_z_yz_xz, g_x_z_yz_yy, g_x_z_yz_yz, g_x_z_yz_zz, g_xy_0_0_0_y_z_yz_xx, g_xy_0_0_0_y_z_yz_xy, g_xy_0_0_0_y_z_yz_xz, g_xy_0_0_0_y_z_yz_yy, g_xy_0_0_0_y_z_yz_yz, g_xy_0_0_0_y_z_yz_zz, g_xyy_z_yz_xx, g_xyy_z_yz_xy, g_xyy_z_yz_xz, g_xyy_z_yz_yy, g_xyy_z_yz_yz, g_xyy_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_z_yz_xx[i] = -2.0 * g_x_z_yz_xx[i] * a_exp + 4.0 * g_xyy_z_yz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_yz_xy[i] = -2.0 * g_x_z_yz_xy[i] * a_exp + 4.0 * g_xyy_z_yz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_yz_xz[i] = -2.0 * g_x_z_yz_xz[i] * a_exp + 4.0 * g_xyy_z_yz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_yz_yy[i] = -2.0 * g_x_z_yz_yy[i] * a_exp + 4.0 * g_xyy_z_yz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_yz_yz[i] = -2.0 * g_x_z_yz_yz[i] * a_exp + 4.0 * g_xyy_z_yz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_yz_zz[i] = -2.0 * g_x_z_yz_zz[i] * a_exp + 4.0 * g_xyy_z_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (534-540)

    #pragma omp simd aligned(g_x_z_zz_xx, g_x_z_zz_xy, g_x_z_zz_xz, g_x_z_zz_yy, g_x_z_zz_yz, g_x_z_zz_zz, g_xy_0_0_0_y_z_zz_xx, g_xy_0_0_0_y_z_zz_xy, g_xy_0_0_0_y_z_zz_xz, g_xy_0_0_0_y_z_zz_yy, g_xy_0_0_0_y_z_zz_yz, g_xy_0_0_0_y_z_zz_zz, g_xyy_z_zz_xx, g_xyy_z_zz_xy, g_xyy_z_zz_xz, g_xyy_z_zz_yy, g_xyy_z_zz_yz, g_xyy_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_z_zz_xx[i] = -2.0 * g_x_z_zz_xx[i] * a_exp + 4.0 * g_xyy_z_zz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_zz_xy[i] = -2.0 * g_x_z_zz_xy[i] * a_exp + 4.0 * g_xyy_z_zz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_zz_xz[i] = -2.0 * g_x_z_zz_xz[i] * a_exp + 4.0 * g_xyy_z_zz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_zz_yy[i] = -2.0 * g_x_z_zz_yy[i] * a_exp + 4.0 * g_xyy_z_zz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_zz_yz[i] = -2.0 * g_x_z_zz_yz[i] * a_exp + 4.0 * g_xyy_z_zz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_zz_zz[i] = -2.0 * g_x_z_zz_zz[i] * a_exp + 4.0 * g_xyy_z_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (540-546)

    #pragma omp simd aligned(g_xy_0_0_0_z_x_xx_xx, g_xy_0_0_0_z_x_xx_xy, g_xy_0_0_0_z_x_xx_xz, g_xy_0_0_0_z_x_xx_yy, g_xy_0_0_0_z_x_xx_yz, g_xy_0_0_0_z_x_xx_zz, g_xyz_x_xx_xx, g_xyz_x_xx_xy, g_xyz_x_xx_xz, g_xyz_x_xx_yy, g_xyz_x_xx_yz, g_xyz_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_x_xx_xx[i] = 4.0 * g_xyz_x_xx_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_xx_xy[i] = 4.0 * g_xyz_x_xx_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_xx_xz[i] = 4.0 * g_xyz_x_xx_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_xx_yy[i] = 4.0 * g_xyz_x_xx_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_xx_yz[i] = 4.0 * g_xyz_x_xx_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_xx_zz[i] = 4.0 * g_xyz_x_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (546-552)

    #pragma omp simd aligned(g_xy_0_0_0_z_x_xy_xx, g_xy_0_0_0_z_x_xy_xy, g_xy_0_0_0_z_x_xy_xz, g_xy_0_0_0_z_x_xy_yy, g_xy_0_0_0_z_x_xy_yz, g_xy_0_0_0_z_x_xy_zz, g_xyz_x_xy_xx, g_xyz_x_xy_xy, g_xyz_x_xy_xz, g_xyz_x_xy_yy, g_xyz_x_xy_yz, g_xyz_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_x_xy_xx[i] = 4.0 * g_xyz_x_xy_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_xy_xy[i] = 4.0 * g_xyz_x_xy_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_xy_xz[i] = 4.0 * g_xyz_x_xy_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_xy_yy[i] = 4.0 * g_xyz_x_xy_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_xy_yz[i] = 4.0 * g_xyz_x_xy_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_xy_zz[i] = 4.0 * g_xyz_x_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (552-558)

    #pragma omp simd aligned(g_xy_0_0_0_z_x_xz_xx, g_xy_0_0_0_z_x_xz_xy, g_xy_0_0_0_z_x_xz_xz, g_xy_0_0_0_z_x_xz_yy, g_xy_0_0_0_z_x_xz_yz, g_xy_0_0_0_z_x_xz_zz, g_xyz_x_xz_xx, g_xyz_x_xz_xy, g_xyz_x_xz_xz, g_xyz_x_xz_yy, g_xyz_x_xz_yz, g_xyz_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_x_xz_xx[i] = 4.0 * g_xyz_x_xz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_xz_xy[i] = 4.0 * g_xyz_x_xz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_xz_xz[i] = 4.0 * g_xyz_x_xz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_xz_yy[i] = 4.0 * g_xyz_x_xz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_xz_yz[i] = 4.0 * g_xyz_x_xz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_xz_zz[i] = 4.0 * g_xyz_x_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (558-564)

    #pragma omp simd aligned(g_xy_0_0_0_z_x_yy_xx, g_xy_0_0_0_z_x_yy_xy, g_xy_0_0_0_z_x_yy_xz, g_xy_0_0_0_z_x_yy_yy, g_xy_0_0_0_z_x_yy_yz, g_xy_0_0_0_z_x_yy_zz, g_xyz_x_yy_xx, g_xyz_x_yy_xy, g_xyz_x_yy_xz, g_xyz_x_yy_yy, g_xyz_x_yy_yz, g_xyz_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_x_yy_xx[i] = 4.0 * g_xyz_x_yy_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_yy_xy[i] = 4.0 * g_xyz_x_yy_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_yy_xz[i] = 4.0 * g_xyz_x_yy_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_yy_yy[i] = 4.0 * g_xyz_x_yy_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_yy_yz[i] = 4.0 * g_xyz_x_yy_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_yy_zz[i] = 4.0 * g_xyz_x_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (564-570)

    #pragma omp simd aligned(g_xy_0_0_0_z_x_yz_xx, g_xy_0_0_0_z_x_yz_xy, g_xy_0_0_0_z_x_yz_xz, g_xy_0_0_0_z_x_yz_yy, g_xy_0_0_0_z_x_yz_yz, g_xy_0_0_0_z_x_yz_zz, g_xyz_x_yz_xx, g_xyz_x_yz_xy, g_xyz_x_yz_xz, g_xyz_x_yz_yy, g_xyz_x_yz_yz, g_xyz_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_x_yz_xx[i] = 4.0 * g_xyz_x_yz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_yz_xy[i] = 4.0 * g_xyz_x_yz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_yz_xz[i] = 4.0 * g_xyz_x_yz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_yz_yy[i] = 4.0 * g_xyz_x_yz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_yz_yz[i] = 4.0 * g_xyz_x_yz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_yz_zz[i] = 4.0 * g_xyz_x_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (570-576)

    #pragma omp simd aligned(g_xy_0_0_0_z_x_zz_xx, g_xy_0_0_0_z_x_zz_xy, g_xy_0_0_0_z_x_zz_xz, g_xy_0_0_0_z_x_zz_yy, g_xy_0_0_0_z_x_zz_yz, g_xy_0_0_0_z_x_zz_zz, g_xyz_x_zz_xx, g_xyz_x_zz_xy, g_xyz_x_zz_xz, g_xyz_x_zz_yy, g_xyz_x_zz_yz, g_xyz_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_x_zz_xx[i] = 4.0 * g_xyz_x_zz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_zz_xy[i] = 4.0 * g_xyz_x_zz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_zz_xz[i] = 4.0 * g_xyz_x_zz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_zz_yy[i] = 4.0 * g_xyz_x_zz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_zz_yz[i] = 4.0 * g_xyz_x_zz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_zz_zz[i] = 4.0 * g_xyz_x_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (576-582)

    #pragma omp simd aligned(g_xy_0_0_0_z_y_xx_xx, g_xy_0_0_0_z_y_xx_xy, g_xy_0_0_0_z_y_xx_xz, g_xy_0_0_0_z_y_xx_yy, g_xy_0_0_0_z_y_xx_yz, g_xy_0_0_0_z_y_xx_zz, g_xyz_y_xx_xx, g_xyz_y_xx_xy, g_xyz_y_xx_xz, g_xyz_y_xx_yy, g_xyz_y_xx_yz, g_xyz_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_y_xx_xx[i] = 4.0 * g_xyz_y_xx_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_xx_xy[i] = 4.0 * g_xyz_y_xx_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_xx_xz[i] = 4.0 * g_xyz_y_xx_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_xx_yy[i] = 4.0 * g_xyz_y_xx_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_xx_yz[i] = 4.0 * g_xyz_y_xx_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_xx_zz[i] = 4.0 * g_xyz_y_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (582-588)

    #pragma omp simd aligned(g_xy_0_0_0_z_y_xy_xx, g_xy_0_0_0_z_y_xy_xy, g_xy_0_0_0_z_y_xy_xz, g_xy_0_0_0_z_y_xy_yy, g_xy_0_0_0_z_y_xy_yz, g_xy_0_0_0_z_y_xy_zz, g_xyz_y_xy_xx, g_xyz_y_xy_xy, g_xyz_y_xy_xz, g_xyz_y_xy_yy, g_xyz_y_xy_yz, g_xyz_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_y_xy_xx[i] = 4.0 * g_xyz_y_xy_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_xy_xy[i] = 4.0 * g_xyz_y_xy_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_xy_xz[i] = 4.0 * g_xyz_y_xy_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_xy_yy[i] = 4.0 * g_xyz_y_xy_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_xy_yz[i] = 4.0 * g_xyz_y_xy_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_xy_zz[i] = 4.0 * g_xyz_y_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (588-594)

    #pragma omp simd aligned(g_xy_0_0_0_z_y_xz_xx, g_xy_0_0_0_z_y_xz_xy, g_xy_0_0_0_z_y_xz_xz, g_xy_0_0_0_z_y_xz_yy, g_xy_0_0_0_z_y_xz_yz, g_xy_0_0_0_z_y_xz_zz, g_xyz_y_xz_xx, g_xyz_y_xz_xy, g_xyz_y_xz_xz, g_xyz_y_xz_yy, g_xyz_y_xz_yz, g_xyz_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_y_xz_xx[i] = 4.0 * g_xyz_y_xz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_xz_xy[i] = 4.0 * g_xyz_y_xz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_xz_xz[i] = 4.0 * g_xyz_y_xz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_xz_yy[i] = 4.0 * g_xyz_y_xz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_xz_yz[i] = 4.0 * g_xyz_y_xz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_xz_zz[i] = 4.0 * g_xyz_y_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (594-600)

    #pragma omp simd aligned(g_xy_0_0_0_z_y_yy_xx, g_xy_0_0_0_z_y_yy_xy, g_xy_0_0_0_z_y_yy_xz, g_xy_0_0_0_z_y_yy_yy, g_xy_0_0_0_z_y_yy_yz, g_xy_0_0_0_z_y_yy_zz, g_xyz_y_yy_xx, g_xyz_y_yy_xy, g_xyz_y_yy_xz, g_xyz_y_yy_yy, g_xyz_y_yy_yz, g_xyz_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_y_yy_xx[i] = 4.0 * g_xyz_y_yy_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_yy_xy[i] = 4.0 * g_xyz_y_yy_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_yy_xz[i] = 4.0 * g_xyz_y_yy_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_yy_yy[i] = 4.0 * g_xyz_y_yy_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_yy_yz[i] = 4.0 * g_xyz_y_yy_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_yy_zz[i] = 4.0 * g_xyz_y_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (600-606)

    #pragma omp simd aligned(g_xy_0_0_0_z_y_yz_xx, g_xy_0_0_0_z_y_yz_xy, g_xy_0_0_0_z_y_yz_xz, g_xy_0_0_0_z_y_yz_yy, g_xy_0_0_0_z_y_yz_yz, g_xy_0_0_0_z_y_yz_zz, g_xyz_y_yz_xx, g_xyz_y_yz_xy, g_xyz_y_yz_xz, g_xyz_y_yz_yy, g_xyz_y_yz_yz, g_xyz_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_y_yz_xx[i] = 4.0 * g_xyz_y_yz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_yz_xy[i] = 4.0 * g_xyz_y_yz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_yz_xz[i] = 4.0 * g_xyz_y_yz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_yz_yy[i] = 4.0 * g_xyz_y_yz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_yz_yz[i] = 4.0 * g_xyz_y_yz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_yz_zz[i] = 4.0 * g_xyz_y_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (606-612)

    #pragma omp simd aligned(g_xy_0_0_0_z_y_zz_xx, g_xy_0_0_0_z_y_zz_xy, g_xy_0_0_0_z_y_zz_xz, g_xy_0_0_0_z_y_zz_yy, g_xy_0_0_0_z_y_zz_yz, g_xy_0_0_0_z_y_zz_zz, g_xyz_y_zz_xx, g_xyz_y_zz_xy, g_xyz_y_zz_xz, g_xyz_y_zz_yy, g_xyz_y_zz_yz, g_xyz_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_y_zz_xx[i] = 4.0 * g_xyz_y_zz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_zz_xy[i] = 4.0 * g_xyz_y_zz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_zz_xz[i] = 4.0 * g_xyz_y_zz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_zz_yy[i] = 4.0 * g_xyz_y_zz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_zz_yz[i] = 4.0 * g_xyz_y_zz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_zz_zz[i] = 4.0 * g_xyz_y_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (612-618)

    #pragma omp simd aligned(g_xy_0_0_0_z_z_xx_xx, g_xy_0_0_0_z_z_xx_xy, g_xy_0_0_0_z_z_xx_xz, g_xy_0_0_0_z_z_xx_yy, g_xy_0_0_0_z_z_xx_yz, g_xy_0_0_0_z_z_xx_zz, g_xyz_z_xx_xx, g_xyz_z_xx_xy, g_xyz_z_xx_xz, g_xyz_z_xx_yy, g_xyz_z_xx_yz, g_xyz_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_z_xx_xx[i] = 4.0 * g_xyz_z_xx_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_xx_xy[i] = 4.0 * g_xyz_z_xx_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_xx_xz[i] = 4.0 * g_xyz_z_xx_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_xx_yy[i] = 4.0 * g_xyz_z_xx_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_xx_yz[i] = 4.0 * g_xyz_z_xx_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_xx_zz[i] = 4.0 * g_xyz_z_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (618-624)

    #pragma omp simd aligned(g_xy_0_0_0_z_z_xy_xx, g_xy_0_0_0_z_z_xy_xy, g_xy_0_0_0_z_z_xy_xz, g_xy_0_0_0_z_z_xy_yy, g_xy_0_0_0_z_z_xy_yz, g_xy_0_0_0_z_z_xy_zz, g_xyz_z_xy_xx, g_xyz_z_xy_xy, g_xyz_z_xy_xz, g_xyz_z_xy_yy, g_xyz_z_xy_yz, g_xyz_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_z_xy_xx[i] = 4.0 * g_xyz_z_xy_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_xy_xy[i] = 4.0 * g_xyz_z_xy_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_xy_xz[i] = 4.0 * g_xyz_z_xy_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_xy_yy[i] = 4.0 * g_xyz_z_xy_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_xy_yz[i] = 4.0 * g_xyz_z_xy_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_xy_zz[i] = 4.0 * g_xyz_z_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (624-630)

    #pragma omp simd aligned(g_xy_0_0_0_z_z_xz_xx, g_xy_0_0_0_z_z_xz_xy, g_xy_0_0_0_z_z_xz_xz, g_xy_0_0_0_z_z_xz_yy, g_xy_0_0_0_z_z_xz_yz, g_xy_0_0_0_z_z_xz_zz, g_xyz_z_xz_xx, g_xyz_z_xz_xy, g_xyz_z_xz_xz, g_xyz_z_xz_yy, g_xyz_z_xz_yz, g_xyz_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_z_xz_xx[i] = 4.0 * g_xyz_z_xz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_xz_xy[i] = 4.0 * g_xyz_z_xz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_xz_xz[i] = 4.0 * g_xyz_z_xz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_xz_yy[i] = 4.0 * g_xyz_z_xz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_xz_yz[i] = 4.0 * g_xyz_z_xz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_xz_zz[i] = 4.0 * g_xyz_z_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (630-636)

    #pragma omp simd aligned(g_xy_0_0_0_z_z_yy_xx, g_xy_0_0_0_z_z_yy_xy, g_xy_0_0_0_z_z_yy_xz, g_xy_0_0_0_z_z_yy_yy, g_xy_0_0_0_z_z_yy_yz, g_xy_0_0_0_z_z_yy_zz, g_xyz_z_yy_xx, g_xyz_z_yy_xy, g_xyz_z_yy_xz, g_xyz_z_yy_yy, g_xyz_z_yy_yz, g_xyz_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_z_yy_xx[i] = 4.0 * g_xyz_z_yy_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_yy_xy[i] = 4.0 * g_xyz_z_yy_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_yy_xz[i] = 4.0 * g_xyz_z_yy_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_yy_yy[i] = 4.0 * g_xyz_z_yy_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_yy_yz[i] = 4.0 * g_xyz_z_yy_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_yy_zz[i] = 4.0 * g_xyz_z_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (636-642)

    #pragma omp simd aligned(g_xy_0_0_0_z_z_yz_xx, g_xy_0_0_0_z_z_yz_xy, g_xy_0_0_0_z_z_yz_xz, g_xy_0_0_0_z_z_yz_yy, g_xy_0_0_0_z_z_yz_yz, g_xy_0_0_0_z_z_yz_zz, g_xyz_z_yz_xx, g_xyz_z_yz_xy, g_xyz_z_yz_xz, g_xyz_z_yz_yy, g_xyz_z_yz_yz, g_xyz_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_z_yz_xx[i] = 4.0 * g_xyz_z_yz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_yz_xy[i] = 4.0 * g_xyz_z_yz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_yz_xz[i] = 4.0 * g_xyz_z_yz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_yz_yy[i] = 4.0 * g_xyz_z_yz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_yz_yz[i] = 4.0 * g_xyz_z_yz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_yz_zz[i] = 4.0 * g_xyz_z_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (642-648)

    #pragma omp simd aligned(g_xy_0_0_0_z_z_zz_xx, g_xy_0_0_0_z_z_zz_xy, g_xy_0_0_0_z_z_zz_xz, g_xy_0_0_0_z_z_zz_yy, g_xy_0_0_0_z_z_zz_yz, g_xy_0_0_0_z_z_zz_zz, g_xyz_z_zz_xx, g_xyz_z_zz_xy, g_xyz_z_zz_xz, g_xyz_z_zz_yy, g_xyz_z_zz_yz, g_xyz_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_z_zz_xx[i] = 4.0 * g_xyz_z_zz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_zz_xy[i] = 4.0 * g_xyz_z_zz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_zz_xz[i] = 4.0 * g_xyz_z_zz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_zz_yy[i] = 4.0 * g_xyz_z_zz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_zz_yz[i] = 4.0 * g_xyz_z_zz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_zz_zz[i] = 4.0 * g_xyz_z_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (648-654)

    #pragma omp simd aligned(g_xxz_x_xx_xx, g_xxz_x_xx_xy, g_xxz_x_xx_xz, g_xxz_x_xx_yy, g_xxz_x_xx_yz, g_xxz_x_xx_zz, g_xz_0_0_0_x_x_xx_xx, g_xz_0_0_0_x_x_xx_xy, g_xz_0_0_0_x_x_xx_xz, g_xz_0_0_0_x_x_xx_yy, g_xz_0_0_0_x_x_xx_yz, g_xz_0_0_0_x_x_xx_zz, g_z_x_xx_xx, g_z_x_xx_xy, g_z_x_xx_xz, g_z_x_xx_yy, g_z_x_xx_yz, g_z_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_x_xx_xx[i] = -2.0 * g_z_x_xx_xx[i] * a_exp + 4.0 * g_xxz_x_xx_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_xx_xy[i] = -2.0 * g_z_x_xx_xy[i] * a_exp + 4.0 * g_xxz_x_xx_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_xx_xz[i] = -2.0 * g_z_x_xx_xz[i] * a_exp + 4.0 * g_xxz_x_xx_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_xx_yy[i] = -2.0 * g_z_x_xx_yy[i] * a_exp + 4.0 * g_xxz_x_xx_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_xx_yz[i] = -2.0 * g_z_x_xx_yz[i] * a_exp + 4.0 * g_xxz_x_xx_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_xx_zz[i] = -2.0 * g_z_x_xx_zz[i] * a_exp + 4.0 * g_xxz_x_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (654-660)

    #pragma omp simd aligned(g_xxz_x_xy_xx, g_xxz_x_xy_xy, g_xxz_x_xy_xz, g_xxz_x_xy_yy, g_xxz_x_xy_yz, g_xxz_x_xy_zz, g_xz_0_0_0_x_x_xy_xx, g_xz_0_0_0_x_x_xy_xy, g_xz_0_0_0_x_x_xy_xz, g_xz_0_0_0_x_x_xy_yy, g_xz_0_0_0_x_x_xy_yz, g_xz_0_0_0_x_x_xy_zz, g_z_x_xy_xx, g_z_x_xy_xy, g_z_x_xy_xz, g_z_x_xy_yy, g_z_x_xy_yz, g_z_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_x_xy_xx[i] = -2.0 * g_z_x_xy_xx[i] * a_exp + 4.0 * g_xxz_x_xy_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_xy_xy[i] = -2.0 * g_z_x_xy_xy[i] * a_exp + 4.0 * g_xxz_x_xy_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_xy_xz[i] = -2.0 * g_z_x_xy_xz[i] * a_exp + 4.0 * g_xxz_x_xy_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_xy_yy[i] = -2.0 * g_z_x_xy_yy[i] * a_exp + 4.0 * g_xxz_x_xy_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_xy_yz[i] = -2.0 * g_z_x_xy_yz[i] * a_exp + 4.0 * g_xxz_x_xy_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_xy_zz[i] = -2.0 * g_z_x_xy_zz[i] * a_exp + 4.0 * g_xxz_x_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (660-666)

    #pragma omp simd aligned(g_xxz_x_xz_xx, g_xxz_x_xz_xy, g_xxz_x_xz_xz, g_xxz_x_xz_yy, g_xxz_x_xz_yz, g_xxz_x_xz_zz, g_xz_0_0_0_x_x_xz_xx, g_xz_0_0_0_x_x_xz_xy, g_xz_0_0_0_x_x_xz_xz, g_xz_0_0_0_x_x_xz_yy, g_xz_0_0_0_x_x_xz_yz, g_xz_0_0_0_x_x_xz_zz, g_z_x_xz_xx, g_z_x_xz_xy, g_z_x_xz_xz, g_z_x_xz_yy, g_z_x_xz_yz, g_z_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_x_xz_xx[i] = -2.0 * g_z_x_xz_xx[i] * a_exp + 4.0 * g_xxz_x_xz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_xz_xy[i] = -2.0 * g_z_x_xz_xy[i] * a_exp + 4.0 * g_xxz_x_xz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_xz_xz[i] = -2.0 * g_z_x_xz_xz[i] * a_exp + 4.0 * g_xxz_x_xz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_xz_yy[i] = -2.0 * g_z_x_xz_yy[i] * a_exp + 4.0 * g_xxz_x_xz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_xz_yz[i] = -2.0 * g_z_x_xz_yz[i] * a_exp + 4.0 * g_xxz_x_xz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_xz_zz[i] = -2.0 * g_z_x_xz_zz[i] * a_exp + 4.0 * g_xxz_x_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (666-672)

    #pragma omp simd aligned(g_xxz_x_yy_xx, g_xxz_x_yy_xy, g_xxz_x_yy_xz, g_xxz_x_yy_yy, g_xxz_x_yy_yz, g_xxz_x_yy_zz, g_xz_0_0_0_x_x_yy_xx, g_xz_0_0_0_x_x_yy_xy, g_xz_0_0_0_x_x_yy_xz, g_xz_0_0_0_x_x_yy_yy, g_xz_0_0_0_x_x_yy_yz, g_xz_0_0_0_x_x_yy_zz, g_z_x_yy_xx, g_z_x_yy_xy, g_z_x_yy_xz, g_z_x_yy_yy, g_z_x_yy_yz, g_z_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_x_yy_xx[i] = -2.0 * g_z_x_yy_xx[i] * a_exp + 4.0 * g_xxz_x_yy_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_yy_xy[i] = -2.0 * g_z_x_yy_xy[i] * a_exp + 4.0 * g_xxz_x_yy_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_yy_xz[i] = -2.0 * g_z_x_yy_xz[i] * a_exp + 4.0 * g_xxz_x_yy_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_yy_yy[i] = -2.0 * g_z_x_yy_yy[i] * a_exp + 4.0 * g_xxz_x_yy_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_yy_yz[i] = -2.0 * g_z_x_yy_yz[i] * a_exp + 4.0 * g_xxz_x_yy_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_yy_zz[i] = -2.0 * g_z_x_yy_zz[i] * a_exp + 4.0 * g_xxz_x_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (672-678)

    #pragma omp simd aligned(g_xxz_x_yz_xx, g_xxz_x_yz_xy, g_xxz_x_yz_xz, g_xxz_x_yz_yy, g_xxz_x_yz_yz, g_xxz_x_yz_zz, g_xz_0_0_0_x_x_yz_xx, g_xz_0_0_0_x_x_yz_xy, g_xz_0_0_0_x_x_yz_xz, g_xz_0_0_0_x_x_yz_yy, g_xz_0_0_0_x_x_yz_yz, g_xz_0_0_0_x_x_yz_zz, g_z_x_yz_xx, g_z_x_yz_xy, g_z_x_yz_xz, g_z_x_yz_yy, g_z_x_yz_yz, g_z_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_x_yz_xx[i] = -2.0 * g_z_x_yz_xx[i] * a_exp + 4.0 * g_xxz_x_yz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_yz_xy[i] = -2.0 * g_z_x_yz_xy[i] * a_exp + 4.0 * g_xxz_x_yz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_yz_xz[i] = -2.0 * g_z_x_yz_xz[i] * a_exp + 4.0 * g_xxz_x_yz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_yz_yy[i] = -2.0 * g_z_x_yz_yy[i] * a_exp + 4.0 * g_xxz_x_yz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_yz_yz[i] = -2.0 * g_z_x_yz_yz[i] * a_exp + 4.0 * g_xxz_x_yz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_yz_zz[i] = -2.0 * g_z_x_yz_zz[i] * a_exp + 4.0 * g_xxz_x_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (678-684)

    #pragma omp simd aligned(g_xxz_x_zz_xx, g_xxz_x_zz_xy, g_xxz_x_zz_xz, g_xxz_x_zz_yy, g_xxz_x_zz_yz, g_xxz_x_zz_zz, g_xz_0_0_0_x_x_zz_xx, g_xz_0_0_0_x_x_zz_xy, g_xz_0_0_0_x_x_zz_xz, g_xz_0_0_0_x_x_zz_yy, g_xz_0_0_0_x_x_zz_yz, g_xz_0_0_0_x_x_zz_zz, g_z_x_zz_xx, g_z_x_zz_xy, g_z_x_zz_xz, g_z_x_zz_yy, g_z_x_zz_yz, g_z_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_x_zz_xx[i] = -2.0 * g_z_x_zz_xx[i] * a_exp + 4.0 * g_xxz_x_zz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_zz_xy[i] = -2.0 * g_z_x_zz_xy[i] * a_exp + 4.0 * g_xxz_x_zz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_zz_xz[i] = -2.0 * g_z_x_zz_xz[i] * a_exp + 4.0 * g_xxz_x_zz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_zz_yy[i] = -2.0 * g_z_x_zz_yy[i] * a_exp + 4.0 * g_xxz_x_zz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_zz_yz[i] = -2.0 * g_z_x_zz_yz[i] * a_exp + 4.0 * g_xxz_x_zz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_zz_zz[i] = -2.0 * g_z_x_zz_zz[i] * a_exp + 4.0 * g_xxz_x_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (684-690)

    #pragma omp simd aligned(g_xxz_y_xx_xx, g_xxz_y_xx_xy, g_xxz_y_xx_xz, g_xxz_y_xx_yy, g_xxz_y_xx_yz, g_xxz_y_xx_zz, g_xz_0_0_0_x_y_xx_xx, g_xz_0_0_0_x_y_xx_xy, g_xz_0_0_0_x_y_xx_xz, g_xz_0_0_0_x_y_xx_yy, g_xz_0_0_0_x_y_xx_yz, g_xz_0_0_0_x_y_xx_zz, g_z_y_xx_xx, g_z_y_xx_xy, g_z_y_xx_xz, g_z_y_xx_yy, g_z_y_xx_yz, g_z_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_y_xx_xx[i] = -2.0 * g_z_y_xx_xx[i] * a_exp + 4.0 * g_xxz_y_xx_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_xx_xy[i] = -2.0 * g_z_y_xx_xy[i] * a_exp + 4.0 * g_xxz_y_xx_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_xx_xz[i] = -2.0 * g_z_y_xx_xz[i] * a_exp + 4.0 * g_xxz_y_xx_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_xx_yy[i] = -2.0 * g_z_y_xx_yy[i] * a_exp + 4.0 * g_xxz_y_xx_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_xx_yz[i] = -2.0 * g_z_y_xx_yz[i] * a_exp + 4.0 * g_xxz_y_xx_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_xx_zz[i] = -2.0 * g_z_y_xx_zz[i] * a_exp + 4.0 * g_xxz_y_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (690-696)

    #pragma omp simd aligned(g_xxz_y_xy_xx, g_xxz_y_xy_xy, g_xxz_y_xy_xz, g_xxz_y_xy_yy, g_xxz_y_xy_yz, g_xxz_y_xy_zz, g_xz_0_0_0_x_y_xy_xx, g_xz_0_0_0_x_y_xy_xy, g_xz_0_0_0_x_y_xy_xz, g_xz_0_0_0_x_y_xy_yy, g_xz_0_0_0_x_y_xy_yz, g_xz_0_0_0_x_y_xy_zz, g_z_y_xy_xx, g_z_y_xy_xy, g_z_y_xy_xz, g_z_y_xy_yy, g_z_y_xy_yz, g_z_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_y_xy_xx[i] = -2.0 * g_z_y_xy_xx[i] * a_exp + 4.0 * g_xxz_y_xy_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_xy_xy[i] = -2.0 * g_z_y_xy_xy[i] * a_exp + 4.0 * g_xxz_y_xy_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_xy_xz[i] = -2.0 * g_z_y_xy_xz[i] * a_exp + 4.0 * g_xxz_y_xy_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_xy_yy[i] = -2.0 * g_z_y_xy_yy[i] * a_exp + 4.0 * g_xxz_y_xy_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_xy_yz[i] = -2.0 * g_z_y_xy_yz[i] * a_exp + 4.0 * g_xxz_y_xy_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_xy_zz[i] = -2.0 * g_z_y_xy_zz[i] * a_exp + 4.0 * g_xxz_y_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (696-702)

    #pragma omp simd aligned(g_xxz_y_xz_xx, g_xxz_y_xz_xy, g_xxz_y_xz_xz, g_xxz_y_xz_yy, g_xxz_y_xz_yz, g_xxz_y_xz_zz, g_xz_0_0_0_x_y_xz_xx, g_xz_0_0_0_x_y_xz_xy, g_xz_0_0_0_x_y_xz_xz, g_xz_0_0_0_x_y_xz_yy, g_xz_0_0_0_x_y_xz_yz, g_xz_0_0_0_x_y_xz_zz, g_z_y_xz_xx, g_z_y_xz_xy, g_z_y_xz_xz, g_z_y_xz_yy, g_z_y_xz_yz, g_z_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_y_xz_xx[i] = -2.0 * g_z_y_xz_xx[i] * a_exp + 4.0 * g_xxz_y_xz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_xz_xy[i] = -2.0 * g_z_y_xz_xy[i] * a_exp + 4.0 * g_xxz_y_xz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_xz_xz[i] = -2.0 * g_z_y_xz_xz[i] * a_exp + 4.0 * g_xxz_y_xz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_xz_yy[i] = -2.0 * g_z_y_xz_yy[i] * a_exp + 4.0 * g_xxz_y_xz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_xz_yz[i] = -2.0 * g_z_y_xz_yz[i] * a_exp + 4.0 * g_xxz_y_xz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_xz_zz[i] = -2.0 * g_z_y_xz_zz[i] * a_exp + 4.0 * g_xxz_y_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (702-708)

    #pragma omp simd aligned(g_xxz_y_yy_xx, g_xxz_y_yy_xy, g_xxz_y_yy_xz, g_xxz_y_yy_yy, g_xxz_y_yy_yz, g_xxz_y_yy_zz, g_xz_0_0_0_x_y_yy_xx, g_xz_0_0_0_x_y_yy_xy, g_xz_0_0_0_x_y_yy_xz, g_xz_0_0_0_x_y_yy_yy, g_xz_0_0_0_x_y_yy_yz, g_xz_0_0_0_x_y_yy_zz, g_z_y_yy_xx, g_z_y_yy_xy, g_z_y_yy_xz, g_z_y_yy_yy, g_z_y_yy_yz, g_z_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_y_yy_xx[i] = -2.0 * g_z_y_yy_xx[i] * a_exp + 4.0 * g_xxz_y_yy_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_yy_xy[i] = -2.0 * g_z_y_yy_xy[i] * a_exp + 4.0 * g_xxz_y_yy_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_yy_xz[i] = -2.0 * g_z_y_yy_xz[i] * a_exp + 4.0 * g_xxz_y_yy_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_yy_yy[i] = -2.0 * g_z_y_yy_yy[i] * a_exp + 4.0 * g_xxz_y_yy_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_yy_yz[i] = -2.0 * g_z_y_yy_yz[i] * a_exp + 4.0 * g_xxz_y_yy_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_yy_zz[i] = -2.0 * g_z_y_yy_zz[i] * a_exp + 4.0 * g_xxz_y_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (708-714)

    #pragma omp simd aligned(g_xxz_y_yz_xx, g_xxz_y_yz_xy, g_xxz_y_yz_xz, g_xxz_y_yz_yy, g_xxz_y_yz_yz, g_xxz_y_yz_zz, g_xz_0_0_0_x_y_yz_xx, g_xz_0_0_0_x_y_yz_xy, g_xz_0_0_0_x_y_yz_xz, g_xz_0_0_0_x_y_yz_yy, g_xz_0_0_0_x_y_yz_yz, g_xz_0_0_0_x_y_yz_zz, g_z_y_yz_xx, g_z_y_yz_xy, g_z_y_yz_xz, g_z_y_yz_yy, g_z_y_yz_yz, g_z_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_y_yz_xx[i] = -2.0 * g_z_y_yz_xx[i] * a_exp + 4.0 * g_xxz_y_yz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_yz_xy[i] = -2.0 * g_z_y_yz_xy[i] * a_exp + 4.0 * g_xxz_y_yz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_yz_xz[i] = -2.0 * g_z_y_yz_xz[i] * a_exp + 4.0 * g_xxz_y_yz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_yz_yy[i] = -2.0 * g_z_y_yz_yy[i] * a_exp + 4.0 * g_xxz_y_yz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_yz_yz[i] = -2.0 * g_z_y_yz_yz[i] * a_exp + 4.0 * g_xxz_y_yz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_yz_zz[i] = -2.0 * g_z_y_yz_zz[i] * a_exp + 4.0 * g_xxz_y_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (714-720)

    #pragma omp simd aligned(g_xxz_y_zz_xx, g_xxz_y_zz_xy, g_xxz_y_zz_xz, g_xxz_y_zz_yy, g_xxz_y_zz_yz, g_xxz_y_zz_zz, g_xz_0_0_0_x_y_zz_xx, g_xz_0_0_0_x_y_zz_xy, g_xz_0_0_0_x_y_zz_xz, g_xz_0_0_0_x_y_zz_yy, g_xz_0_0_0_x_y_zz_yz, g_xz_0_0_0_x_y_zz_zz, g_z_y_zz_xx, g_z_y_zz_xy, g_z_y_zz_xz, g_z_y_zz_yy, g_z_y_zz_yz, g_z_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_y_zz_xx[i] = -2.0 * g_z_y_zz_xx[i] * a_exp + 4.0 * g_xxz_y_zz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_zz_xy[i] = -2.0 * g_z_y_zz_xy[i] * a_exp + 4.0 * g_xxz_y_zz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_zz_xz[i] = -2.0 * g_z_y_zz_xz[i] * a_exp + 4.0 * g_xxz_y_zz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_zz_yy[i] = -2.0 * g_z_y_zz_yy[i] * a_exp + 4.0 * g_xxz_y_zz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_zz_yz[i] = -2.0 * g_z_y_zz_yz[i] * a_exp + 4.0 * g_xxz_y_zz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_zz_zz[i] = -2.0 * g_z_y_zz_zz[i] * a_exp + 4.0 * g_xxz_y_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (720-726)

    #pragma omp simd aligned(g_xxz_z_xx_xx, g_xxz_z_xx_xy, g_xxz_z_xx_xz, g_xxz_z_xx_yy, g_xxz_z_xx_yz, g_xxz_z_xx_zz, g_xz_0_0_0_x_z_xx_xx, g_xz_0_0_0_x_z_xx_xy, g_xz_0_0_0_x_z_xx_xz, g_xz_0_0_0_x_z_xx_yy, g_xz_0_0_0_x_z_xx_yz, g_xz_0_0_0_x_z_xx_zz, g_z_z_xx_xx, g_z_z_xx_xy, g_z_z_xx_xz, g_z_z_xx_yy, g_z_z_xx_yz, g_z_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_z_xx_xx[i] = -2.0 * g_z_z_xx_xx[i] * a_exp + 4.0 * g_xxz_z_xx_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_xx_xy[i] = -2.0 * g_z_z_xx_xy[i] * a_exp + 4.0 * g_xxz_z_xx_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_xx_xz[i] = -2.0 * g_z_z_xx_xz[i] * a_exp + 4.0 * g_xxz_z_xx_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_xx_yy[i] = -2.0 * g_z_z_xx_yy[i] * a_exp + 4.0 * g_xxz_z_xx_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_xx_yz[i] = -2.0 * g_z_z_xx_yz[i] * a_exp + 4.0 * g_xxz_z_xx_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_xx_zz[i] = -2.0 * g_z_z_xx_zz[i] * a_exp + 4.0 * g_xxz_z_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (726-732)

    #pragma omp simd aligned(g_xxz_z_xy_xx, g_xxz_z_xy_xy, g_xxz_z_xy_xz, g_xxz_z_xy_yy, g_xxz_z_xy_yz, g_xxz_z_xy_zz, g_xz_0_0_0_x_z_xy_xx, g_xz_0_0_0_x_z_xy_xy, g_xz_0_0_0_x_z_xy_xz, g_xz_0_0_0_x_z_xy_yy, g_xz_0_0_0_x_z_xy_yz, g_xz_0_0_0_x_z_xy_zz, g_z_z_xy_xx, g_z_z_xy_xy, g_z_z_xy_xz, g_z_z_xy_yy, g_z_z_xy_yz, g_z_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_z_xy_xx[i] = -2.0 * g_z_z_xy_xx[i] * a_exp + 4.0 * g_xxz_z_xy_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_xy_xy[i] = -2.0 * g_z_z_xy_xy[i] * a_exp + 4.0 * g_xxz_z_xy_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_xy_xz[i] = -2.0 * g_z_z_xy_xz[i] * a_exp + 4.0 * g_xxz_z_xy_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_xy_yy[i] = -2.0 * g_z_z_xy_yy[i] * a_exp + 4.0 * g_xxz_z_xy_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_xy_yz[i] = -2.0 * g_z_z_xy_yz[i] * a_exp + 4.0 * g_xxz_z_xy_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_xy_zz[i] = -2.0 * g_z_z_xy_zz[i] * a_exp + 4.0 * g_xxz_z_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (732-738)

    #pragma omp simd aligned(g_xxz_z_xz_xx, g_xxz_z_xz_xy, g_xxz_z_xz_xz, g_xxz_z_xz_yy, g_xxz_z_xz_yz, g_xxz_z_xz_zz, g_xz_0_0_0_x_z_xz_xx, g_xz_0_0_0_x_z_xz_xy, g_xz_0_0_0_x_z_xz_xz, g_xz_0_0_0_x_z_xz_yy, g_xz_0_0_0_x_z_xz_yz, g_xz_0_0_0_x_z_xz_zz, g_z_z_xz_xx, g_z_z_xz_xy, g_z_z_xz_xz, g_z_z_xz_yy, g_z_z_xz_yz, g_z_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_z_xz_xx[i] = -2.0 * g_z_z_xz_xx[i] * a_exp + 4.0 * g_xxz_z_xz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_xz_xy[i] = -2.0 * g_z_z_xz_xy[i] * a_exp + 4.0 * g_xxz_z_xz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_xz_xz[i] = -2.0 * g_z_z_xz_xz[i] * a_exp + 4.0 * g_xxz_z_xz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_xz_yy[i] = -2.0 * g_z_z_xz_yy[i] * a_exp + 4.0 * g_xxz_z_xz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_xz_yz[i] = -2.0 * g_z_z_xz_yz[i] * a_exp + 4.0 * g_xxz_z_xz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_xz_zz[i] = -2.0 * g_z_z_xz_zz[i] * a_exp + 4.0 * g_xxz_z_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (738-744)

    #pragma omp simd aligned(g_xxz_z_yy_xx, g_xxz_z_yy_xy, g_xxz_z_yy_xz, g_xxz_z_yy_yy, g_xxz_z_yy_yz, g_xxz_z_yy_zz, g_xz_0_0_0_x_z_yy_xx, g_xz_0_0_0_x_z_yy_xy, g_xz_0_0_0_x_z_yy_xz, g_xz_0_0_0_x_z_yy_yy, g_xz_0_0_0_x_z_yy_yz, g_xz_0_0_0_x_z_yy_zz, g_z_z_yy_xx, g_z_z_yy_xy, g_z_z_yy_xz, g_z_z_yy_yy, g_z_z_yy_yz, g_z_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_z_yy_xx[i] = -2.0 * g_z_z_yy_xx[i] * a_exp + 4.0 * g_xxz_z_yy_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_yy_xy[i] = -2.0 * g_z_z_yy_xy[i] * a_exp + 4.0 * g_xxz_z_yy_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_yy_xz[i] = -2.0 * g_z_z_yy_xz[i] * a_exp + 4.0 * g_xxz_z_yy_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_yy_yy[i] = -2.0 * g_z_z_yy_yy[i] * a_exp + 4.0 * g_xxz_z_yy_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_yy_yz[i] = -2.0 * g_z_z_yy_yz[i] * a_exp + 4.0 * g_xxz_z_yy_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_yy_zz[i] = -2.0 * g_z_z_yy_zz[i] * a_exp + 4.0 * g_xxz_z_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (744-750)

    #pragma omp simd aligned(g_xxz_z_yz_xx, g_xxz_z_yz_xy, g_xxz_z_yz_xz, g_xxz_z_yz_yy, g_xxz_z_yz_yz, g_xxz_z_yz_zz, g_xz_0_0_0_x_z_yz_xx, g_xz_0_0_0_x_z_yz_xy, g_xz_0_0_0_x_z_yz_xz, g_xz_0_0_0_x_z_yz_yy, g_xz_0_0_0_x_z_yz_yz, g_xz_0_0_0_x_z_yz_zz, g_z_z_yz_xx, g_z_z_yz_xy, g_z_z_yz_xz, g_z_z_yz_yy, g_z_z_yz_yz, g_z_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_z_yz_xx[i] = -2.0 * g_z_z_yz_xx[i] * a_exp + 4.0 * g_xxz_z_yz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_yz_xy[i] = -2.0 * g_z_z_yz_xy[i] * a_exp + 4.0 * g_xxz_z_yz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_yz_xz[i] = -2.0 * g_z_z_yz_xz[i] * a_exp + 4.0 * g_xxz_z_yz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_yz_yy[i] = -2.0 * g_z_z_yz_yy[i] * a_exp + 4.0 * g_xxz_z_yz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_yz_yz[i] = -2.0 * g_z_z_yz_yz[i] * a_exp + 4.0 * g_xxz_z_yz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_yz_zz[i] = -2.0 * g_z_z_yz_zz[i] * a_exp + 4.0 * g_xxz_z_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (750-756)

    #pragma omp simd aligned(g_xxz_z_zz_xx, g_xxz_z_zz_xy, g_xxz_z_zz_xz, g_xxz_z_zz_yy, g_xxz_z_zz_yz, g_xxz_z_zz_zz, g_xz_0_0_0_x_z_zz_xx, g_xz_0_0_0_x_z_zz_xy, g_xz_0_0_0_x_z_zz_xz, g_xz_0_0_0_x_z_zz_yy, g_xz_0_0_0_x_z_zz_yz, g_xz_0_0_0_x_z_zz_zz, g_z_z_zz_xx, g_z_z_zz_xy, g_z_z_zz_xz, g_z_z_zz_yy, g_z_z_zz_yz, g_z_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_z_zz_xx[i] = -2.0 * g_z_z_zz_xx[i] * a_exp + 4.0 * g_xxz_z_zz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_zz_xy[i] = -2.0 * g_z_z_zz_xy[i] * a_exp + 4.0 * g_xxz_z_zz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_zz_xz[i] = -2.0 * g_z_z_zz_xz[i] * a_exp + 4.0 * g_xxz_z_zz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_zz_yy[i] = -2.0 * g_z_z_zz_yy[i] * a_exp + 4.0 * g_xxz_z_zz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_zz_yz[i] = -2.0 * g_z_z_zz_yz[i] * a_exp + 4.0 * g_xxz_z_zz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_zz_zz[i] = -2.0 * g_z_z_zz_zz[i] * a_exp + 4.0 * g_xxz_z_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (756-762)

    #pragma omp simd aligned(g_xyz_x_xx_xx, g_xyz_x_xx_xy, g_xyz_x_xx_xz, g_xyz_x_xx_yy, g_xyz_x_xx_yz, g_xyz_x_xx_zz, g_xz_0_0_0_y_x_xx_xx, g_xz_0_0_0_y_x_xx_xy, g_xz_0_0_0_y_x_xx_xz, g_xz_0_0_0_y_x_xx_yy, g_xz_0_0_0_y_x_xx_yz, g_xz_0_0_0_y_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_x_xx_xx[i] = 4.0 * g_xyz_x_xx_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_xx_xy[i] = 4.0 * g_xyz_x_xx_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_xx_xz[i] = 4.0 * g_xyz_x_xx_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_xx_yy[i] = 4.0 * g_xyz_x_xx_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_xx_yz[i] = 4.0 * g_xyz_x_xx_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_xx_zz[i] = 4.0 * g_xyz_x_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (762-768)

    #pragma omp simd aligned(g_xyz_x_xy_xx, g_xyz_x_xy_xy, g_xyz_x_xy_xz, g_xyz_x_xy_yy, g_xyz_x_xy_yz, g_xyz_x_xy_zz, g_xz_0_0_0_y_x_xy_xx, g_xz_0_0_0_y_x_xy_xy, g_xz_0_0_0_y_x_xy_xz, g_xz_0_0_0_y_x_xy_yy, g_xz_0_0_0_y_x_xy_yz, g_xz_0_0_0_y_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_x_xy_xx[i] = 4.0 * g_xyz_x_xy_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_xy_xy[i] = 4.0 * g_xyz_x_xy_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_xy_xz[i] = 4.0 * g_xyz_x_xy_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_xy_yy[i] = 4.0 * g_xyz_x_xy_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_xy_yz[i] = 4.0 * g_xyz_x_xy_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_xy_zz[i] = 4.0 * g_xyz_x_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (768-774)

    #pragma omp simd aligned(g_xyz_x_xz_xx, g_xyz_x_xz_xy, g_xyz_x_xz_xz, g_xyz_x_xz_yy, g_xyz_x_xz_yz, g_xyz_x_xz_zz, g_xz_0_0_0_y_x_xz_xx, g_xz_0_0_0_y_x_xz_xy, g_xz_0_0_0_y_x_xz_xz, g_xz_0_0_0_y_x_xz_yy, g_xz_0_0_0_y_x_xz_yz, g_xz_0_0_0_y_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_x_xz_xx[i] = 4.0 * g_xyz_x_xz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_xz_xy[i] = 4.0 * g_xyz_x_xz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_xz_xz[i] = 4.0 * g_xyz_x_xz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_xz_yy[i] = 4.0 * g_xyz_x_xz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_xz_yz[i] = 4.0 * g_xyz_x_xz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_xz_zz[i] = 4.0 * g_xyz_x_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (774-780)

    #pragma omp simd aligned(g_xyz_x_yy_xx, g_xyz_x_yy_xy, g_xyz_x_yy_xz, g_xyz_x_yy_yy, g_xyz_x_yy_yz, g_xyz_x_yy_zz, g_xz_0_0_0_y_x_yy_xx, g_xz_0_0_0_y_x_yy_xy, g_xz_0_0_0_y_x_yy_xz, g_xz_0_0_0_y_x_yy_yy, g_xz_0_0_0_y_x_yy_yz, g_xz_0_0_0_y_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_x_yy_xx[i] = 4.0 * g_xyz_x_yy_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_yy_xy[i] = 4.0 * g_xyz_x_yy_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_yy_xz[i] = 4.0 * g_xyz_x_yy_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_yy_yy[i] = 4.0 * g_xyz_x_yy_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_yy_yz[i] = 4.0 * g_xyz_x_yy_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_yy_zz[i] = 4.0 * g_xyz_x_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (780-786)

    #pragma omp simd aligned(g_xyz_x_yz_xx, g_xyz_x_yz_xy, g_xyz_x_yz_xz, g_xyz_x_yz_yy, g_xyz_x_yz_yz, g_xyz_x_yz_zz, g_xz_0_0_0_y_x_yz_xx, g_xz_0_0_0_y_x_yz_xy, g_xz_0_0_0_y_x_yz_xz, g_xz_0_0_0_y_x_yz_yy, g_xz_0_0_0_y_x_yz_yz, g_xz_0_0_0_y_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_x_yz_xx[i] = 4.0 * g_xyz_x_yz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_yz_xy[i] = 4.0 * g_xyz_x_yz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_yz_xz[i] = 4.0 * g_xyz_x_yz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_yz_yy[i] = 4.0 * g_xyz_x_yz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_yz_yz[i] = 4.0 * g_xyz_x_yz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_yz_zz[i] = 4.0 * g_xyz_x_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (786-792)

    #pragma omp simd aligned(g_xyz_x_zz_xx, g_xyz_x_zz_xy, g_xyz_x_zz_xz, g_xyz_x_zz_yy, g_xyz_x_zz_yz, g_xyz_x_zz_zz, g_xz_0_0_0_y_x_zz_xx, g_xz_0_0_0_y_x_zz_xy, g_xz_0_0_0_y_x_zz_xz, g_xz_0_0_0_y_x_zz_yy, g_xz_0_0_0_y_x_zz_yz, g_xz_0_0_0_y_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_x_zz_xx[i] = 4.0 * g_xyz_x_zz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_zz_xy[i] = 4.0 * g_xyz_x_zz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_zz_xz[i] = 4.0 * g_xyz_x_zz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_zz_yy[i] = 4.0 * g_xyz_x_zz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_zz_yz[i] = 4.0 * g_xyz_x_zz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_zz_zz[i] = 4.0 * g_xyz_x_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (792-798)

    #pragma omp simd aligned(g_xyz_y_xx_xx, g_xyz_y_xx_xy, g_xyz_y_xx_xz, g_xyz_y_xx_yy, g_xyz_y_xx_yz, g_xyz_y_xx_zz, g_xz_0_0_0_y_y_xx_xx, g_xz_0_0_0_y_y_xx_xy, g_xz_0_0_0_y_y_xx_xz, g_xz_0_0_0_y_y_xx_yy, g_xz_0_0_0_y_y_xx_yz, g_xz_0_0_0_y_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_y_xx_xx[i] = 4.0 * g_xyz_y_xx_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_xx_xy[i] = 4.0 * g_xyz_y_xx_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_xx_xz[i] = 4.0 * g_xyz_y_xx_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_xx_yy[i] = 4.0 * g_xyz_y_xx_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_xx_yz[i] = 4.0 * g_xyz_y_xx_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_xx_zz[i] = 4.0 * g_xyz_y_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (798-804)

    #pragma omp simd aligned(g_xyz_y_xy_xx, g_xyz_y_xy_xy, g_xyz_y_xy_xz, g_xyz_y_xy_yy, g_xyz_y_xy_yz, g_xyz_y_xy_zz, g_xz_0_0_0_y_y_xy_xx, g_xz_0_0_0_y_y_xy_xy, g_xz_0_0_0_y_y_xy_xz, g_xz_0_0_0_y_y_xy_yy, g_xz_0_0_0_y_y_xy_yz, g_xz_0_0_0_y_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_y_xy_xx[i] = 4.0 * g_xyz_y_xy_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_xy_xy[i] = 4.0 * g_xyz_y_xy_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_xy_xz[i] = 4.0 * g_xyz_y_xy_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_xy_yy[i] = 4.0 * g_xyz_y_xy_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_xy_yz[i] = 4.0 * g_xyz_y_xy_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_xy_zz[i] = 4.0 * g_xyz_y_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (804-810)

    #pragma omp simd aligned(g_xyz_y_xz_xx, g_xyz_y_xz_xy, g_xyz_y_xz_xz, g_xyz_y_xz_yy, g_xyz_y_xz_yz, g_xyz_y_xz_zz, g_xz_0_0_0_y_y_xz_xx, g_xz_0_0_0_y_y_xz_xy, g_xz_0_0_0_y_y_xz_xz, g_xz_0_0_0_y_y_xz_yy, g_xz_0_0_0_y_y_xz_yz, g_xz_0_0_0_y_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_y_xz_xx[i] = 4.0 * g_xyz_y_xz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_xz_xy[i] = 4.0 * g_xyz_y_xz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_xz_xz[i] = 4.0 * g_xyz_y_xz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_xz_yy[i] = 4.0 * g_xyz_y_xz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_xz_yz[i] = 4.0 * g_xyz_y_xz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_xz_zz[i] = 4.0 * g_xyz_y_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (810-816)

    #pragma omp simd aligned(g_xyz_y_yy_xx, g_xyz_y_yy_xy, g_xyz_y_yy_xz, g_xyz_y_yy_yy, g_xyz_y_yy_yz, g_xyz_y_yy_zz, g_xz_0_0_0_y_y_yy_xx, g_xz_0_0_0_y_y_yy_xy, g_xz_0_0_0_y_y_yy_xz, g_xz_0_0_0_y_y_yy_yy, g_xz_0_0_0_y_y_yy_yz, g_xz_0_0_0_y_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_y_yy_xx[i] = 4.0 * g_xyz_y_yy_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_yy_xy[i] = 4.0 * g_xyz_y_yy_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_yy_xz[i] = 4.0 * g_xyz_y_yy_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_yy_yy[i] = 4.0 * g_xyz_y_yy_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_yy_yz[i] = 4.0 * g_xyz_y_yy_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_yy_zz[i] = 4.0 * g_xyz_y_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (816-822)

    #pragma omp simd aligned(g_xyz_y_yz_xx, g_xyz_y_yz_xy, g_xyz_y_yz_xz, g_xyz_y_yz_yy, g_xyz_y_yz_yz, g_xyz_y_yz_zz, g_xz_0_0_0_y_y_yz_xx, g_xz_0_0_0_y_y_yz_xy, g_xz_0_0_0_y_y_yz_xz, g_xz_0_0_0_y_y_yz_yy, g_xz_0_0_0_y_y_yz_yz, g_xz_0_0_0_y_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_y_yz_xx[i] = 4.0 * g_xyz_y_yz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_yz_xy[i] = 4.0 * g_xyz_y_yz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_yz_xz[i] = 4.0 * g_xyz_y_yz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_yz_yy[i] = 4.0 * g_xyz_y_yz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_yz_yz[i] = 4.0 * g_xyz_y_yz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_yz_zz[i] = 4.0 * g_xyz_y_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (822-828)

    #pragma omp simd aligned(g_xyz_y_zz_xx, g_xyz_y_zz_xy, g_xyz_y_zz_xz, g_xyz_y_zz_yy, g_xyz_y_zz_yz, g_xyz_y_zz_zz, g_xz_0_0_0_y_y_zz_xx, g_xz_0_0_0_y_y_zz_xy, g_xz_0_0_0_y_y_zz_xz, g_xz_0_0_0_y_y_zz_yy, g_xz_0_0_0_y_y_zz_yz, g_xz_0_0_0_y_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_y_zz_xx[i] = 4.0 * g_xyz_y_zz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_zz_xy[i] = 4.0 * g_xyz_y_zz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_zz_xz[i] = 4.0 * g_xyz_y_zz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_zz_yy[i] = 4.0 * g_xyz_y_zz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_zz_yz[i] = 4.0 * g_xyz_y_zz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_zz_zz[i] = 4.0 * g_xyz_y_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (828-834)

    #pragma omp simd aligned(g_xyz_z_xx_xx, g_xyz_z_xx_xy, g_xyz_z_xx_xz, g_xyz_z_xx_yy, g_xyz_z_xx_yz, g_xyz_z_xx_zz, g_xz_0_0_0_y_z_xx_xx, g_xz_0_0_0_y_z_xx_xy, g_xz_0_0_0_y_z_xx_xz, g_xz_0_0_0_y_z_xx_yy, g_xz_0_0_0_y_z_xx_yz, g_xz_0_0_0_y_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_z_xx_xx[i] = 4.0 * g_xyz_z_xx_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_xx_xy[i] = 4.0 * g_xyz_z_xx_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_xx_xz[i] = 4.0 * g_xyz_z_xx_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_xx_yy[i] = 4.0 * g_xyz_z_xx_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_xx_yz[i] = 4.0 * g_xyz_z_xx_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_xx_zz[i] = 4.0 * g_xyz_z_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (834-840)

    #pragma omp simd aligned(g_xyz_z_xy_xx, g_xyz_z_xy_xy, g_xyz_z_xy_xz, g_xyz_z_xy_yy, g_xyz_z_xy_yz, g_xyz_z_xy_zz, g_xz_0_0_0_y_z_xy_xx, g_xz_0_0_0_y_z_xy_xy, g_xz_0_0_0_y_z_xy_xz, g_xz_0_0_0_y_z_xy_yy, g_xz_0_0_0_y_z_xy_yz, g_xz_0_0_0_y_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_z_xy_xx[i] = 4.0 * g_xyz_z_xy_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_xy_xy[i] = 4.0 * g_xyz_z_xy_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_xy_xz[i] = 4.0 * g_xyz_z_xy_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_xy_yy[i] = 4.0 * g_xyz_z_xy_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_xy_yz[i] = 4.0 * g_xyz_z_xy_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_xy_zz[i] = 4.0 * g_xyz_z_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (840-846)

    #pragma omp simd aligned(g_xyz_z_xz_xx, g_xyz_z_xz_xy, g_xyz_z_xz_xz, g_xyz_z_xz_yy, g_xyz_z_xz_yz, g_xyz_z_xz_zz, g_xz_0_0_0_y_z_xz_xx, g_xz_0_0_0_y_z_xz_xy, g_xz_0_0_0_y_z_xz_xz, g_xz_0_0_0_y_z_xz_yy, g_xz_0_0_0_y_z_xz_yz, g_xz_0_0_0_y_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_z_xz_xx[i] = 4.0 * g_xyz_z_xz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_xz_xy[i] = 4.0 * g_xyz_z_xz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_xz_xz[i] = 4.0 * g_xyz_z_xz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_xz_yy[i] = 4.0 * g_xyz_z_xz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_xz_yz[i] = 4.0 * g_xyz_z_xz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_xz_zz[i] = 4.0 * g_xyz_z_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (846-852)

    #pragma omp simd aligned(g_xyz_z_yy_xx, g_xyz_z_yy_xy, g_xyz_z_yy_xz, g_xyz_z_yy_yy, g_xyz_z_yy_yz, g_xyz_z_yy_zz, g_xz_0_0_0_y_z_yy_xx, g_xz_0_0_0_y_z_yy_xy, g_xz_0_0_0_y_z_yy_xz, g_xz_0_0_0_y_z_yy_yy, g_xz_0_0_0_y_z_yy_yz, g_xz_0_0_0_y_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_z_yy_xx[i] = 4.0 * g_xyz_z_yy_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_yy_xy[i] = 4.0 * g_xyz_z_yy_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_yy_xz[i] = 4.0 * g_xyz_z_yy_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_yy_yy[i] = 4.0 * g_xyz_z_yy_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_yy_yz[i] = 4.0 * g_xyz_z_yy_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_yy_zz[i] = 4.0 * g_xyz_z_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (852-858)

    #pragma omp simd aligned(g_xyz_z_yz_xx, g_xyz_z_yz_xy, g_xyz_z_yz_xz, g_xyz_z_yz_yy, g_xyz_z_yz_yz, g_xyz_z_yz_zz, g_xz_0_0_0_y_z_yz_xx, g_xz_0_0_0_y_z_yz_xy, g_xz_0_0_0_y_z_yz_xz, g_xz_0_0_0_y_z_yz_yy, g_xz_0_0_0_y_z_yz_yz, g_xz_0_0_0_y_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_z_yz_xx[i] = 4.0 * g_xyz_z_yz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_yz_xy[i] = 4.0 * g_xyz_z_yz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_yz_xz[i] = 4.0 * g_xyz_z_yz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_yz_yy[i] = 4.0 * g_xyz_z_yz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_yz_yz[i] = 4.0 * g_xyz_z_yz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_yz_zz[i] = 4.0 * g_xyz_z_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (858-864)

    #pragma omp simd aligned(g_xyz_z_zz_xx, g_xyz_z_zz_xy, g_xyz_z_zz_xz, g_xyz_z_zz_yy, g_xyz_z_zz_yz, g_xyz_z_zz_zz, g_xz_0_0_0_y_z_zz_xx, g_xz_0_0_0_y_z_zz_xy, g_xz_0_0_0_y_z_zz_xz, g_xz_0_0_0_y_z_zz_yy, g_xz_0_0_0_y_z_zz_yz, g_xz_0_0_0_y_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_z_zz_xx[i] = 4.0 * g_xyz_z_zz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_zz_xy[i] = 4.0 * g_xyz_z_zz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_zz_xz[i] = 4.0 * g_xyz_z_zz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_zz_yy[i] = 4.0 * g_xyz_z_zz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_zz_yz[i] = 4.0 * g_xyz_z_zz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_zz_zz[i] = 4.0 * g_xyz_z_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (864-870)

    #pragma omp simd aligned(g_x_x_xx_xx, g_x_x_xx_xy, g_x_x_xx_xz, g_x_x_xx_yy, g_x_x_xx_yz, g_x_x_xx_zz, g_xz_0_0_0_z_x_xx_xx, g_xz_0_0_0_z_x_xx_xy, g_xz_0_0_0_z_x_xx_xz, g_xz_0_0_0_z_x_xx_yy, g_xz_0_0_0_z_x_xx_yz, g_xz_0_0_0_z_x_xx_zz, g_xzz_x_xx_xx, g_xzz_x_xx_xy, g_xzz_x_xx_xz, g_xzz_x_xx_yy, g_xzz_x_xx_yz, g_xzz_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_x_xx_xx[i] = -2.0 * g_x_x_xx_xx[i] * a_exp + 4.0 * g_xzz_x_xx_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_xx_xy[i] = -2.0 * g_x_x_xx_xy[i] * a_exp + 4.0 * g_xzz_x_xx_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_xx_xz[i] = -2.0 * g_x_x_xx_xz[i] * a_exp + 4.0 * g_xzz_x_xx_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_xx_yy[i] = -2.0 * g_x_x_xx_yy[i] * a_exp + 4.0 * g_xzz_x_xx_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_xx_yz[i] = -2.0 * g_x_x_xx_yz[i] * a_exp + 4.0 * g_xzz_x_xx_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_xx_zz[i] = -2.0 * g_x_x_xx_zz[i] * a_exp + 4.0 * g_xzz_x_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (870-876)

    #pragma omp simd aligned(g_x_x_xy_xx, g_x_x_xy_xy, g_x_x_xy_xz, g_x_x_xy_yy, g_x_x_xy_yz, g_x_x_xy_zz, g_xz_0_0_0_z_x_xy_xx, g_xz_0_0_0_z_x_xy_xy, g_xz_0_0_0_z_x_xy_xz, g_xz_0_0_0_z_x_xy_yy, g_xz_0_0_0_z_x_xy_yz, g_xz_0_0_0_z_x_xy_zz, g_xzz_x_xy_xx, g_xzz_x_xy_xy, g_xzz_x_xy_xz, g_xzz_x_xy_yy, g_xzz_x_xy_yz, g_xzz_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_x_xy_xx[i] = -2.0 * g_x_x_xy_xx[i] * a_exp + 4.0 * g_xzz_x_xy_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_xy_xy[i] = -2.0 * g_x_x_xy_xy[i] * a_exp + 4.0 * g_xzz_x_xy_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_xy_xz[i] = -2.0 * g_x_x_xy_xz[i] * a_exp + 4.0 * g_xzz_x_xy_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_xy_yy[i] = -2.0 * g_x_x_xy_yy[i] * a_exp + 4.0 * g_xzz_x_xy_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_xy_yz[i] = -2.0 * g_x_x_xy_yz[i] * a_exp + 4.0 * g_xzz_x_xy_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_xy_zz[i] = -2.0 * g_x_x_xy_zz[i] * a_exp + 4.0 * g_xzz_x_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (876-882)

    #pragma omp simd aligned(g_x_x_xz_xx, g_x_x_xz_xy, g_x_x_xz_xz, g_x_x_xz_yy, g_x_x_xz_yz, g_x_x_xz_zz, g_xz_0_0_0_z_x_xz_xx, g_xz_0_0_0_z_x_xz_xy, g_xz_0_0_0_z_x_xz_xz, g_xz_0_0_0_z_x_xz_yy, g_xz_0_0_0_z_x_xz_yz, g_xz_0_0_0_z_x_xz_zz, g_xzz_x_xz_xx, g_xzz_x_xz_xy, g_xzz_x_xz_xz, g_xzz_x_xz_yy, g_xzz_x_xz_yz, g_xzz_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_x_xz_xx[i] = -2.0 * g_x_x_xz_xx[i] * a_exp + 4.0 * g_xzz_x_xz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_xz_xy[i] = -2.0 * g_x_x_xz_xy[i] * a_exp + 4.0 * g_xzz_x_xz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_xz_xz[i] = -2.0 * g_x_x_xz_xz[i] * a_exp + 4.0 * g_xzz_x_xz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_xz_yy[i] = -2.0 * g_x_x_xz_yy[i] * a_exp + 4.0 * g_xzz_x_xz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_xz_yz[i] = -2.0 * g_x_x_xz_yz[i] * a_exp + 4.0 * g_xzz_x_xz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_xz_zz[i] = -2.0 * g_x_x_xz_zz[i] * a_exp + 4.0 * g_xzz_x_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (882-888)

    #pragma omp simd aligned(g_x_x_yy_xx, g_x_x_yy_xy, g_x_x_yy_xz, g_x_x_yy_yy, g_x_x_yy_yz, g_x_x_yy_zz, g_xz_0_0_0_z_x_yy_xx, g_xz_0_0_0_z_x_yy_xy, g_xz_0_0_0_z_x_yy_xz, g_xz_0_0_0_z_x_yy_yy, g_xz_0_0_0_z_x_yy_yz, g_xz_0_0_0_z_x_yy_zz, g_xzz_x_yy_xx, g_xzz_x_yy_xy, g_xzz_x_yy_xz, g_xzz_x_yy_yy, g_xzz_x_yy_yz, g_xzz_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_x_yy_xx[i] = -2.0 * g_x_x_yy_xx[i] * a_exp + 4.0 * g_xzz_x_yy_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_yy_xy[i] = -2.0 * g_x_x_yy_xy[i] * a_exp + 4.0 * g_xzz_x_yy_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_yy_xz[i] = -2.0 * g_x_x_yy_xz[i] * a_exp + 4.0 * g_xzz_x_yy_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_yy_yy[i] = -2.0 * g_x_x_yy_yy[i] * a_exp + 4.0 * g_xzz_x_yy_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_yy_yz[i] = -2.0 * g_x_x_yy_yz[i] * a_exp + 4.0 * g_xzz_x_yy_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_yy_zz[i] = -2.0 * g_x_x_yy_zz[i] * a_exp + 4.0 * g_xzz_x_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (888-894)

    #pragma omp simd aligned(g_x_x_yz_xx, g_x_x_yz_xy, g_x_x_yz_xz, g_x_x_yz_yy, g_x_x_yz_yz, g_x_x_yz_zz, g_xz_0_0_0_z_x_yz_xx, g_xz_0_0_0_z_x_yz_xy, g_xz_0_0_0_z_x_yz_xz, g_xz_0_0_0_z_x_yz_yy, g_xz_0_0_0_z_x_yz_yz, g_xz_0_0_0_z_x_yz_zz, g_xzz_x_yz_xx, g_xzz_x_yz_xy, g_xzz_x_yz_xz, g_xzz_x_yz_yy, g_xzz_x_yz_yz, g_xzz_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_x_yz_xx[i] = -2.0 * g_x_x_yz_xx[i] * a_exp + 4.0 * g_xzz_x_yz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_yz_xy[i] = -2.0 * g_x_x_yz_xy[i] * a_exp + 4.0 * g_xzz_x_yz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_yz_xz[i] = -2.0 * g_x_x_yz_xz[i] * a_exp + 4.0 * g_xzz_x_yz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_yz_yy[i] = -2.0 * g_x_x_yz_yy[i] * a_exp + 4.0 * g_xzz_x_yz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_yz_yz[i] = -2.0 * g_x_x_yz_yz[i] * a_exp + 4.0 * g_xzz_x_yz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_yz_zz[i] = -2.0 * g_x_x_yz_zz[i] * a_exp + 4.0 * g_xzz_x_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (894-900)

    #pragma omp simd aligned(g_x_x_zz_xx, g_x_x_zz_xy, g_x_x_zz_xz, g_x_x_zz_yy, g_x_x_zz_yz, g_x_x_zz_zz, g_xz_0_0_0_z_x_zz_xx, g_xz_0_0_0_z_x_zz_xy, g_xz_0_0_0_z_x_zz_xz, g_xz_0_0_0_z_x_zz_yy, g_xz_0_0_0_z_x_zz_yz, g_xz_0_0_0_z_x_zz_zz, g_xzz_x_zz_xx, g_xzz_x_zz_xy, g_xzz_x_zz_xz, g_xzz_x_zz_yy, g_xzz_x_zz_yz, g_xzz_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_x_zz_xx[i] = -2.0 * g_x_x_zz_xx[i] * a_exp + 4.0 * g_xzz_x_zz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_zz_xy[i] = -2.0 * g_x_x_zz_xy[i] * a_exp + 4.0 * g_xzz_x_zz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_zz_xz[i] = -2.0 * g_x_x_zz_xz[i] * a_exp + 4.0 * g_xzz_x_zz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_zz_yy[i] = -2.0 * g_x_x_zz_yy[i] * a_exp + 4.0 * g_xzz_x_zz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_zz_yz[i] = -2.0 * g_x_x_zz_yz[i] * a_exp + 4.0 * g_xzz_x_zz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_zz_zz[i] = -2.0 * g_x_x_zz_zz[i] * a_exp + 4.0 * g_xzz_x_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (900-906)

    #pragma omp simd aligned(g_x_y_xx_xx, g_x_y_xx_xy, g_x_y_xx_xz, g_x_y_xx_yy, g_x_y_xx_yz, g_x_y_xx_zz, g_xz_0_0_0_z_y_xx_xx, g_xz_0_0_0_z_y_xx_xy, g_xz_0_0_0_z_y_xx_xz, g_xz_0_0_0_z_y_xx_yy, g_xz_0_0_0_z_y_xx_yz, g_xz_0_0_0_z_y_xx_zz, g_xzz_y_xx_xx, g_xzz_y_xx_xy, g_xzz_y_xx_xz, g_xzz_y_xx_yy, g_xzz_y_xx_yz, g_xzz_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_y_xx_xx[i] = -2.0 * g_x_y_xx_xx[i] * a_exp + 4.0 * g_xzz_y_xx_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_xx_xy[i] = -2.0 * g_x_y_xx_xy[i] * a_exp + 4.0 * g_xzz_y_xx_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_xx_xz[i] = -2.0 * g_x_y_xx_xz[i] * a_exp + 4.0 * g_xzz_y_xx_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_xx_yy[i] = -2.0 * g_x_y_xx_yy[i] * a_exp + 4.0 * g_xzz_y_xx_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_xx_yz[i] = -2.0 * g_x_y_xx_yz[i] * a_exp + 4.0 * g_xzz_y_xx_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_xx_zz[i] = -2.0 * g_x_y_xx_zz[i] * a_exp + 4.0 * g_xzz_y_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (906-912)

    #pragma omp simd aligned(g_x_y_xy_xx, g_x_y_xy_xy, g_x_y_xy_xz, g_x_y_xy_yy, g_x_y_xy_yz, g_x_y_xy_zz, g_xz_0_0_0_z_y_xy_xx, g_xz_0_0_0_z_y_xy_xy, g_xz_0_0_0_z_y_xy_xz, g_xz_0_0_0_z_y_xy_yy, g_xz_0_0_0_z_y_xy_yz, g_xz_0_0_0_z_y_xy_zz, g_xzz_y_xy_xx, g_xzz_y_xy_xy, g_xzz_y_xy_xz, g_xzz_y_xy_yy, g_xzz_y_xy_yz, g_xzz_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_y_xy_xx[i] = -2.0 * g_x_y_xy_xx[i] * a_exp + 4.0 * g_xzz_y_xy_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_xy_xy[i] = -2.0 * g_x_y_xy_xy[i] * a_exp + 4.0 * g_xzz_y_xy_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_xy_xz[i] = -2.0 * g_x_y_xy_xz[i] * a_exp + 4.0 * g_xzz_y_xy_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_xy_yy[i] = -2.0 * g_x_y_xy_yy[i] * a_exp + 4.0 * g_xzz_y_xy_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_xy_yz[i] = -2.0 * g_x_y_xy_yz[i] * a_exp + 4.0 * g_xzz_y_xy_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_xy_zz[i] = -2.0 * g_x_y_xy_zz[i] * a_exp + 4.0 * g_xzz_y_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (912-918)

    #pragma omp simd aligned(g_x_y_xz_xx, g_x_y_xz_xy, g_x_y_xz_xz, g_x_y_xz_yy, g_x_y_xz_yz, g_x_y_xz_zz, g_xz_0_0_0_z_y_xz_xx, g_xz_0_0_0_z_y_xz_xy, g_xz_0_0_0_z_y_xz_xz, g_xz_0_0_0_z_y_xz_yy, g_xz_0_0_0_z_y_xz_yz, g_xz_0_0_0_z_y_xz_zz, g_xzz_y_xz_xx, g_xzz_y_xz_xy, g_xzz_y_xz_xz, g_xzz_y_xz_yy, g_xzz_y_xz_yz, g_xzz_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_y_xz_xx[i] = -2.0 * g_x_y_xz_xx[i] * a_exp + 4.0 * g_xzz_y_xz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_xz_xy[i] = -2.0 * g_x_y_xz_xy[i] * a_exp + 4.0 * g_xzz_y_xz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_xz_xz[i] = -2.0 * g_x_y_xz_xz[i] * a_exp + 4.0 * g_xzz_y_xz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_xz_yy[i] = -2.0 * g_x_y_xz_yy[i] * a_exp + 4.0 * g_xzz_y_xz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_xz_yz[i] = -2.0 * g_x_y_xz_yz[i] * a_exp + 4.0 * g_xzz_y_xz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_xz_zz[i] = -2.0 * g_x_y_xz_zz[i] * a_exp + 4.0 * g_xzz_y_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (918-924)

    #pragma omp simd aligned(g_x_y_yy_xx, g_x_y_yy_xy, g_x_y_yy_xz, g_x_y_yy_yy, g_x_y_yy_yz, g_x_y_yy_zz, g_xz_0_0_0_z_y_yy_xx, g_xz_0_0_0_z_y_yy_xy, g_xz_0_0_0_z_y_yy_xz, g_xz_0_0_0_z_y_yy_yy, g_xz_0_0_0_z_y_yy_yz, g_xz_0_0_0_z_y_yy_zz, g_xzz_y_yy_xx, g_xzz_y_yy_xy, g_xzz_y_yy_xz, g_xzz_y_yy_yy, g_xzz_y_yy_yz, g_xzz_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_y_yy_xx[i] = -2.0 * g_x_y_yy_xx[i] * a_exp + 4.0 * g_xzz_y_yy_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_yy_xy[i] = -2.0 * g_x_y_yy_xy[i] * a_exp + 4.0 * g_xzz_y_yy_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_yy_xz[i] = -2.0 * g_x_y_yy_xz[i] * a_exp + 4.0 * g_xzz_y_yy_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_yy_yy[i] = -2.0 * g_x_y_yy_yy[i] * a_exp + 4.0 * g_xzz_y_yy_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_yy_yz[i] = -2.0 * g_x_y_yy_yz[i] * a_exp + 4.0 * g_xzz_y_yy_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_yy_zz[i] = -2.0 * g_x_y_yy_zz[i] * a_exp + 4.0 * g_xzz_y_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (924-930)

    #pragma omp simd aligned(g_x_y_yz_xx, g_x_y_yz_xy, g_x_y_yz_xz, g_x_y_yz_yy, g_x_y_yz_yz, g_x_y_yz_zz, g_xz_0_0_0_z_y_yz_xx, g_xz_0_0_0_z_y_yz_xy, g_xz_0_0_0_z_y_yz_xz, g_xz_0_0_0_z_y_yz_yy, g_xz_0_0_0_z_y_yz_yz, g_xz_0_0_0_z_y_yz_zz, g_xzz_y_yz_xx, g_xzz_y_yz_xy, g_xzz_y_yz_xz, g_xzz_y_yz_yy, g_xzz_y_yz_yz, g_xzz_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_y_yz_xx[i] = -2.0 * g_x_y_yz_xx[i] * a_exp + 4.0 * g_xzz_y_yz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_yz_xy[i] = -2.0 * g_x_y_yz_xy[i] * a_exp + 4.0 * g_xzz_y_yz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_yz_xz[i] = -2.0 * g_x_y_yz_xz[i] * a_exp + 4.0 * g_xzz_y_yz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_yz_yy[i] = -2.0 * g_x_y_yz_yy[i] * a_exp + 4.0 * g_xzz_y_yz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_yz_yz[i] = -2.0 * g_x_y_yz_yz[i] * a_exp + 4.0 * g_xzz_y_yz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_yz_zz[i] = -2.0 * g_x_y_yz_zz[i] * a_exp + 4.0 * g_xzz_y_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (930-936)

    #pragma omp simd aligned(g_x_y_zz_xx, g_x_y_zz_xy, g_x_y_zz_xz, g_x_y_zz_yy, g_x_y_zz_yz, g_x_y_zz_zz, g_xz_0_0_0_z_y_zz_xx, g_xz_0_0_0_z_y_zz_xy, g_xz_0_0_0_z_y_zz_xz, g_xz_0_0_0_z_y_zz_yy, g_xz_0_0_0_z_y_zz_yz, g_xz_0_0_0_z_y_zz_zz, g_xzz_y_zz_xx, g_xzz_y_zz_xy, g_xzz_y_zz_xz, g_xzz_y_zz_yy, g_xzz_y_zz_yz, g_xzz_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_y_zz_xx[i] = -2.0 * g_x_y_zz_xx[i] * a_exp + 4.0 * g_xzz_y_zz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_zz_xy[i] = -2.0 * g_x_y_zz_xy[i] * a_exp + 4.0 * g_xzz_y_zz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_zz_xz[i] = -2.0 * g_x_y_zz_xz[i] * a_exp + 4.0 * g_xzz_y_zz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_zz_yy[i] = -2.0 * g_x_y_zz_yy[i] * a_exp + 4.0 * g_xzz_y_zz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_zz_yz[i] = -2.0 * g_x_y_zz_yz[i] * a_exp + 4.0 * g_xzz_y_zz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_zz_zz[i] = -2.0 * g_x_y_zz_zz[i] * a_exp + 4.0 * g_xzz_y_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (936-942)

    #pragma omp simd aligned(g_x_z_xx_xx, g_x_z_xx_xy, g_x_z_xx_xz, g_x_z_xx_yy, g_x_z_xx_yz, g_x_z_xx_zz, g_xz_0_0_0_z_z_xx_xx, g_xz_0_0_0_z_z_xx_xy, g_xz_0_0_0_z_z_xx_xz, g_xz_0_0_0_z_z_xx_yy, g_xz_0_0_0_z_z_xx_yz, g_xz_0_0_0_z_z_xx_zz, g_xzz_z_xx_xx, g_xzz_z_xx_xy, g_xzz_z_xx_xz, g_xzz_z_xx_yy, g_xzz_z_xx_yz, g_xzz_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_z_xx_xx[i] = -2.0 * g_x_z_xx_xx[i] * a_exp + 4.0 * g_xzz_z_xx_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_xx_xy[i] = -2.0 * g_x_z_xx_xy[i] * a_exp + 4.0 * g_xzz_z_xx_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_xx_xz[i] = -2.0 * g_x_z_xx_xz[i] * a_exp + 4.0 * g_xzz_z_xx_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_xx_yy[i] = -2.0 * g_x_z_xx_yy[i] * a_exp + 4.0 * g_xzz_z_xx_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_xx_yz[i] = -2.0 * g_x_z_xx_yz[i] * a_exp + 4.0 * g_xzz_z_xx_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_xx_zz[i] = -2.0 * g_x_z_xx_zz[i] * a_exp + 4.0 * g_xzz_z_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (942-948)

    #pragma omp simd aligned(g_x_z_xy_xx, g_x_z_xy_xy, g_x_z_xy_xz, g_x_z_xy_yy, g_x_z_xy_yz, g_x_z_xy_zz, g_xz_0_0_0_z_z_xy_xx, g_xz_0_0_0_z_z_xy_xy, g_xz_0_0_0_z_z_xy_xz, g_xz_0_0_0_z_z_xy_yy, g_xz_0_0_0_z_z_xy_yz, g_xz_0_0_0_z_z_xy_zz, g_xzz_z_xy_xx, g_xzz_z_xy_xy, g_xzz_z_xy_xz, g_xzz_z_xy_yy, g_xzz_z_xy_yz, g_xzz_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_z_xy_xx[i] = -2.0 * g_x_z_xy_xx[i] * a_exp + 4.0 * g_xzz_z_xy_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_xy_xy[i] = -2.0 * g_x_z_xy_xy[i] * a_exp + 4.0 * g_xzz_z_xy_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_xy_xz[i] = -2.0 * g_x_z_xy_xz[i] * a_exp + 4.0 * g_xzz_z_xy_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_xy_yy[i] = -2.0 * g_x_z_xy_yy[i] * a_exp + 4.0 * g_xzz_z_xy_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_xy_yz[i] = -2.0 * g_x_z_xy_yz[i] * a_exp + 4.0 * g_xzz_z_xy_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_xy_zz[i] = -2.0 * g_x_z_xy_zz[i] * a_exp + 4.0 * g_xzz_z_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (948-954)

    #pragma omp simd aligned(g_x_z_xz_xx, g_x_z_xz_xy, g_x_z_xz_xz, g_x_z_xz_yy, g_x_z_xz_yz, g_x_z_xz_zz, g_xz_0_0_0_z_z_xz_xx, g_xz_0_0_0_z_z_xz_xy, g_xz_0_0_0_z_z_xz_xz, g_xz_0_0_0_z_z_xz_yy, g_xz_0_0_0_z_z_xz_yz, g_xz_0_0_0_z_z_xz_zz, g_xzz_z_xz_xx, g_xzz_z_xz_xy, g_xzz_z_xz_xz, g_xzz_z_xz_yy, g_xzz_z_xz_yz, g_xzz_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_z_xz_xx[i] = -2.0 * g_x_z_xz_xx[i] * a_exp + 4.0 * g_xzz_z_xz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_xz_xy[i] = -2.0 * g_x_z_xz_xy[i] * a_exp + 4.0 * g_xzz_z_xz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_xz_xz[i] = -2.0 * g_x_z_xz_xz[i] * a_exp + 4.0 * g_xzz_z_xz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_xz_yy[i] = -2.0 * g_x_z_xz_yy[i] * a_exp + 4.0 * g_xzz_z_xz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_xz_yz[i] = -2.0 * g_x_z_xz_yz[i] * a_exp + 4.0 * g_xzz_z_xz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_xz_zz[i] = -2.0 * g_x_z_xz_zz[i] * a_exp + 4.0 * g_xzz_z_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (954-960)

    #pragma omp simd aligned(g_x_z_yy_xx, g_x_z_yy_xy, g_x_z_yy_xz, g_x_z_yy_yy, g_x_z_yy_yz, g_x_z_yy_zz, g_xz_0_0_0_z_z_yy_xx, g_xz_0_0_0_z_z_yy_xy, g_xz_0_0_0_z_z_yy_xz, g_xz_0_0_0_z_z_yy_yy, g_xz_0_0_0_z_z_yy_yz, g_xz_0_0_0_z_z_yy_zz, g_xzz_z_yy_xx, g_xzz_z_yy_xy, g_xzz_z_yy_xz, g_xzz_z_yy_yy, g_xzz_z_yy_yz, g_xzz_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_z_yy_xx[i] = -2.0 * g_x_z_yy_xx[i] * a_exp + 4.0 * g_xzz_z_yy_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_yy_xy[i] = -2.0 * g_x_z_yy_xy[i] * a_exp + 4.0 * g_xzz_z_yy_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_yy_xz[i] = -2.0 * g_x_z_yy_xz[i] * a_exp + 4.0 * g_xzz_z_yy_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_yy_yy[i] = -2.0 * g_x_z_yy_yy[i] * a_exp + 4.0 * g_xzz_z_yy_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_yy_yz[i] = -2.0 * g_x_z_yy_yz[i] * a_exp + 4.0 * g_xzz_z_yy_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_yy_zz[i] = -2.0 * g_x_z_yy_zz[i] * a_exp + 4.0 * g_xzz_z_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (960-966)

    #pragma omp simd aligned(g_x_z_yz_xx, g_x_z_yz_xy, g_x_z_yz_xz, g_x_z_yz_yy, g_x_z_yz_yz, g_x_z_yz_zz, g_xz_0_0_0_z_z_yz_xx, g_xz_0_0_0_z_z_yz_xy, g_xz_0_0_0_z_z_yz_xz, g_xz_0_0_0_z_z_yz_yy, g_xz_0_0_0_z_z_yz_yz, g_xz_0_0_0_z_z_yz_zz, g_xzz_z_yz_xx, g_xzz_z_yz_xy, g_xzz_z_yz_xz, g_xzz_z_yz_yy, g_xzz_z_yz_yz, g_xzz_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_z_yz_xx[i] = -2.0 * g_x_z_yz_xx[i] * a_exp + 4.0 * g_xzz_z_yz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_yz_xy[i] = -2.0 * g_x_z_yz_xy[i] * a_exp + 4.0 * g_xzz_z_yz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_yz_xz[i] = -2.0 * g_x_z_yz_xz[i] * a_exp + 4.0 * g_xzz_z_yz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_yz_yy[i] = -2.0 * g_x_z_yz_yy[i] * a_exp + 4.0 * g_xzz_z_yz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_yz_yz[i] = -2.0 * g_x_z_yz_yz[i] * a_exp + 4.0 * g_xzz_z_yz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_yz_zz[i] = -2.0 * g_x_z_yz_zz[i] * a_exp + 4.0 * g_xzz_z_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (966-972)

    #pragma omp simd aligned(g_x_z_zz_xx, g_x_z_zz_xy, g_x_z_zz_xz, g_x_z_zz_yy, g_x_z_zz_yz, g_x_z_zz_zz, g_xz_0_0_0_z_z_zz_xx, g_xz_0_0_0_z_z_zz_xy, g_xz_0_0_0_z_z_zz_xz, g_xz_0_0_0_z_z_zz_yy, g_xz_0_0_0_z_z_zz_yz, g_xz_0_0_0_z_z_zz_zz, g_xzz_z_zz_xx, g_xzz_z_zz_xy, g_xzz_z_zz_xz, g_xzz_z_zz_yy, g_xzz_z_zz_yz, g_xzz_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_z_zz_xx[i] = -2.0 * g_x_z_zz_xx[i] * a_exp + 4.0 * g_xzz_z_zz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_zz_xy[i] = -2.0 * g_x_z_zz_xy[i] * a_exp + 4.0 * g_xzz_z_zz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_zz_xz[i] = -2.0 * g_x_z_zz_xz[i] * a_exp + 4.0 * g_xzz_z_zz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_zz_yy[i] = -2.0 * g_x_z_zz_yy[i] * a_exp + 4.0 * g_xzz_z_zz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_zz_yz[i] = -2.0 * g_x_z_zz_yz[i] * a_exp + 4.0 * g_xzz_z_zz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_zz_zz[i] = -2.0 * g_x_z_zz_zz[i] * a_exp + 4.0 * g_xzz_z_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (972-978)

    #pragma omp simd aligned(g_x_x_xx_xx, g_x_x_xx_xy, g_x_x_xx_xz, g_x_x_xx_yy, g_x_x_xx_yz, g_x_x_xx_zz, g_xyy_x_xx_xx, g_xyy_x_xx_xy, g_xyy_x_xx_xz, g_xyy_x_xx_yy, g_xyy_x_xx_yz, g_xyy_x_xx_zz, g_yy_0_0_0_x_x_xx_xx, g_yy_0_0_0_x_x_xx_xy, g_yy_0_0_0_x_x_xx_xz, g_yy_0_0_0_x_x_xx_yy, g_yy_0_0_0_x_x_xx_yz, g_yy_0_0_0_x_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_x_xx_xx[i] = -2.0 * g_x_x_xx_xx[i] * a_exp + 4.0 * g_xyy_x_xx_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_xx_xy[i] = -2.0 * g_x_x_xx_xy[i] * a_exp + 4.0 * g_xyy_x_xx_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_xx_xz[i] = -2.0 * g_x_x_xx_xz[i] * a_exp + 4.0 * g_xyy_x_xx_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_xx_yy[i] = -2.0 * g_x_x_xx_yy[i] * a_exp + 4.0 * g_xyy_x_xx_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_xx_yz[i] = -2.0 * g_x_x_xx_yz[i] * a_exp + 4.0 * g_xyy_x_xx_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_xx_zz[i] = -2.0 * g_x_x_xx_zz[i] * a_exp + 4.0 * g_xyy_x_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (978-984)

    #pragma omp simd aligned(g_x_x_xy_xx, g_x_x_xy_xy, g_x_x_xy_xz, g_x_x_xy_yy, g_x_x_xy_yz, g_x_x_xy_zz, g_xyy_x_xy_xx, g_xyy_x_xy_xy, g_xyy_x_xy_xz, g_xyy_x_xy_yy, g_xyy_x_xy_yz, g_xyy_x_xy_zz, g_yy_0_0_0_x_x_xy_xx, g_yy_0_0_0_x_x_xy_xy, g_yy_0_0_0_x_x_xy_xz, g_yy_0_0_0_x_x_xy_yy, g_yy_0_0_0_x_x_xy_yz, g_yy_0_0_0_x_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_x_xy_xx[i] = -2.0 * g_x_x_xy_xx[i] * a_exp + 4.0 * g_xyy_x_xy_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_xy_xy[i] = -2.0 * g_x_x_xy_xy[i] * a_exp + 4.0 * g_xyy_x_xy_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_xy_xz[i] = -2.0 * g_x_x_xy_xz[i] * a_exp + 4.0 * g_xyy_x_xy_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_xy_yy[i] = -2.0 * g_x_x_xy_yy[i] * a_exp + 4.0 * g_xyy_x_xy_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_xy_yz[i] = -2.0 * g_x_x_xy_yz[i] * a_exp + 4.0 * g_xyy_x_xy_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_xy_zz[i] = -2.0 * g_x_x_xy_zz[i] * a_exp + 4.0 * g_xyy_x_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (984-990)

    #pragma omp simd aligned(g_x_x_xz_xx, g_x_x_xz_xy, g_x_x_xz_xz, g_x_x_xz_yy, g_x_x_xz_yz, g_x_x_xz_zz, g_xyy_x_xz_xx, g_xyy_x_xz_xy, g_xyy_x_xz_xz, g_xyy_x_xz_yy, g_xyy_x_xz_yz, g_xyy_x_xz_zz, g_yy_0_0_0_x_x_xz_xx, g_yy_0_0_0_x_x_xz_xy, g_yy_0_0_0_x_x_xz_xz, g_yy_0_0_0_x_x_xz_yy, g_yy_0_0_0_x_x_xz_yz, g_yy_0_0_0_x_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_x_xz_xx[i] = -2.0 * g_x_x_xz_xx[i] * a_exp + 4.0 * g_xyy_x_xz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_xz_xy[i] = -2.0 * g_x_x_xz_xy[i] * a_exp + 4.0 * g_xyy_x_xz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_xz_xz[i] = -2.0 * g_x_x_xz_xz[i] * a_exp + 4.0 * g_xyy_x_xz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_xz_yy[i] = -2.0 * g_x_x_xz_yy[i] * a_exp + 4.0 * g_xyy_x_xz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_xz_yz[i] = -2.0 * g_x_x_xz_yz[i] * a_exp + 4.0 * g_xyy_x_xz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_xz_zz[i] = -2.0 * g_x_x_xz_zz[i] * a_exp + 4.0 * g_xyy_x_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (990-996)

    #pragma omp simd aligned(g_x_x_yy_xx, g_x_x_yy_xy, g_x_x_yy_xz, g_x_x_yy_yy, g_x_x_yy_yz, g_x_x_yy_zz, g_xyy_x_yy_xx, g_xyy_x_yy_xy, g_xyy_x_yy_xz, g_xyy_x_yy_yy, g_xyy_x_yy_yz, g_xyy_x_yy_zz, g_yy_0_0_0_x_x_yy_xx, g_yy_0_0_0_x_x_yy_xy, g_yy_0_0_0_x_x_yy_xz, g_yy_0_0_0_x_x_yy_yy, g_yy_0_0_0_x_x_yy_yz, g_yy_0_0_0_x_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_x_yy_xx[i] = -2.0 * g_x_x_yy_xx[i] * a_exp + 4.0 * g_xyy_x_yy_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_yy_xy[i] = -2.0 * g_x_x_yy_xy[i] * a_exp + 4.0 * g_xyy_x_yy_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_yy_xz[i] = -2.0 * g_x_x_yy_xz[i] * a_exp + 4.0 * g_xyy_x_yy_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_yy_yy[i] = -2.0 * g_x_x_yy_yy[i] * a_exp + 4.0 * g_xyy_x_yy_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_yy_yz[i] = -2.0 * g_x_x_yy_yz[i] * a_exp + 4.0 * g_xyy_x_yy_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_yy_zz[i] = -2.0 * g_x_x_yy_zz[i] * a_exp + 4.0 * g_xyy_x_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (996-1002)

    #pragma omp simd aligned(g_x_x_yz_xx, g_x_x_yz_xy, g_x_x_yz_xz, g_x_x_yz_yy, g_x_x_yz_yz, g_x_x_yz_zz, g_xyy_x_yz_xx, g_xyy_x_yz_xy, g_xyy_x_yz_xz, g_xyy_x_yz_yy, g_xyy_x_yz_yz, g_xyy_x_yz_zz, g_yy_0_0_0_x_x_yz_xx, g_yy_0_0_0_x_x_yz_xy, g_yy_0_0_0_x_x_yz_xz, g_yy_0_0_0_x_x_yz_yy, g_yy_0_0_0_x_x_yz_yz, g_yy_0_0_0_x_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_x_yz_xx[i] = -2.0 * g_x_x_yz_xx[i] * a_exp + 4.0 * g_xyy_x_yz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_yz_xy[i] = -2.0 * g_x_x_yz_xy[i] * a_exp + 4.0 * g_xyy_x_yz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_yz_xz[i] = -2.0 * g_x_x_yz_xz[i] * a_exp + 4.0 * g_xyy_x_yz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_yz_yy[i] = -2.0 * g_x_x_yz_yy[i] * a_exp + 4.0 * g_xyy_x_yz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_yz_yz[i] = -2.0 * g_x_x_yz_yz[i] * a_exp + 4.0 * g_xyy_x_yz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_yz_zz[i] = -2.0 * g_x_x_yz_zz[i] * a_exp + 4.0 * g_xyy_x_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1002-1008)

    #pragma omp simd aligned(g_x_x_zz_xx, g_x_x_zz_xy, g_x_x_zz_xz, g_x_x_zz_yy, g_x_x_zz_yz, g_x_x_zz_zz, g_xyy_x_zz_xx, g_xyy_x_zz_xy, g_xyy_x_zz_xz, g_xyy_x_zz_yy, g_xyy_x_zz_yz, g_xyy_x_zz_zz, g_yy_0_0_0_x_x_zz_xx, g_yy_0_0_0_x_x_zz_xy, g_yy_0_0_0_x_x_zz_xz, g_yy_0_0_0_x_x_zz_yy, g_yy_0_0_0_x_x_zz_yz, g_yy_0_0_0_x_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_x_zz_xx[i] = -2.0 * g_x_x_zz_xx[i] * a_exp + 4.0 * g_xyy_x_zz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_zz_xy[i] = -2.0 * g_x_x_zz_xy[i] * a_exp + 4.0 * g_xyy_x_zz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_zz_xz[i] = -2.0 * g_x_x_zz_xz[i] * a_exp + 4.0 * g_xyy_x_zz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_zz_yy[i] = -2.0 * g_x_x_zz_yy[i] * a_exp + 4.0 * g_xyy_x_zz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_zz_yz[i] = -2.0 * g_x_x_zz_yz[i] * a_exp + 4.0 * g_xyy_x_zz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_zz_zz[i] = -2.0 * g_x_x_zz_zz[i] * a_exp + 4.0 * g_xyy_x_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1008-1014)

    #pragma omp simd aligned(g_x_y_xx_xx, g_x_y_xx_xy, g_x_y_xx_xz, g_x_y_xx_yy, g_x_y_xx_yz, g_x_y_xx_zz, g_xyy_y_xx_xx, g_xyy_y_xx_xy, g_xyy_y_xx_xz, g_xyy_y_xx_yy, g_xyy_y_xx_yz, g_xyy_y_xx_zz, g_yy_0_0_0_x_y_xx_xx, g_yy_0_0_0_x_y_xx_xy, g_yy_0_0_0_x_y_xx_xz, g_yy_0_0_0_x_y_xx_yy, g_yy_0_0_0_x_y_xx_yz, g_yy_0_0_0_x_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_y_xx_xx[i] = -2.0 * g_x_y_xx_xx[i] * a_exp + 4.0 * g_xyy_y_xx_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_xx_xy[i] = -2.0 * g_x_y_xx_xy[i] * a_exp + 4.0 * g_xyy_y_xx_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_xx_xz[i] = -2.0 * g_x_y_xx_xz[i] * a_exp + 4.0 * g_xyy_y_xx_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_xx_yy[i] = -2.0 * g_x_y_xx_yy[i] * a_exp + 4.0 * g_xyy_y_xx_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_xx_yz[i] = -2.0 * g_x_y_xx_yz[i] * a_exp + 4.0 * g_xyy_y_xx_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_xx_zz[i] = -2.0 * g_x_y_xx_zz[i] * a_exp + 4.0 * g_xyy_y_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (1014-1020)

    #pragma omp simd aligned(g_x_y_xy_xx, g_x_y_xy_xy, g_x_y_xy_xz, g_x_y_xy_yy, g_x_y_xy_yz, g_x_y_xy_zz, g_xyy_y_xy_xx, g_xyy_y_xy_xy, g_xyy_y_xy_xz, g_xyy_y_xy_yy, g_xyy_y_xy_yz, g_xyy_y_xy_zz, g_yy_0_0_0_x_y_xy_xx, g_yy_0_0_0_x_y_xy_xy, g_yy_0_0_0_x_y_xy_xz, g_yy_0_0_0_x_y_xy_yy, g_yy_0_0_0_x_y_xy_yz, g_yy_0_0_0_x_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_y_xy_xx[i] = -2.0 * g_x_y_xy_xx[i] * a_exp + 4.0 * g_xyy_y_xy_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_xy_xy[i] = -2.0 * g_x_y_xy_xy[i] * a_exp + 4.0 * g_xyy_y_xy_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_xy_xz[i] = -2.0 * g_x_y_xy_xz[i] * a_exp + 4.0 * g_xyy_y_xy_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_xy_yy[i] = -2.0 * g_x_y_xy_yy[i] * a_exp + 4.0 * g_xyy_y_xy_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_xy_yz[i] = -2.0 * g_x_y_xy_yz[i] * a_exp + 4.0 * g_xyy_y_xy_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_xy_zz[i] = -2.0 * g_x_y_xy_zz[i] * a_exp + 4.0 * g_xyy_y_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1020-1026)

    #pragma omp simd aligned(g_x_y_xz_xx, g_x_y_xz_xy, g_x_y_xz_xz, g_x_y_xz_yy, g_x_y_xz_yz, g_x_y_xz_zz, g_xyy_y_xz_xx, g_xyy_y_xz_xy, g_xyy_y_xz_xz, g_xyy_y_xz_yy, g_xyy_y_xz_yz, g_xyy_y_xz_zz, g_yy_0_0_0_x_y_xz_xx, g_yy_0_0_0_x_y_xz_xy, g_yy_0_0_0_x_y_xz_xz, g_yy_0_0_0_x_y_xz_yy, g_yy_0_0_0_x_y_xz_yz, g_yy_0_0_0_x_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_y_xz_xx[i] = -2.0 * g_x_y_xz_xx[i] * a_exp + 4.0 * g_xyy_y_xz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_xz_xy[i] = -2.0 * g_x_y_xz_xy[i] * a_exp + 4.0 * g_xyy_y_xz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_xz_xz[i] = -2.0 * g_x_y_xz_xz[i] * a_exp + 4.0 * g_xyy_y_xz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_xz_yy[i] = -2.0 * g_x_y_xz_yy[i] * a_exp + 4.0 * g_xyy_y_xz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_xz_yz[i] = -2.0 * g_x_y_xz_yz[i] * a_exp + 4.0 * g_xyy_y_xz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_xz_zz[i] = -2.0 * g_x_y_xz_zz[i] * a_exp + 4.0 * g_xyy_y_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1026-1032)

    #pragma omp simd aligned(g_x_y_yy_xx, g_x_y_yy_xy, g_x_y_yy_xz, g_x_y_yy_yy, g_x_y_yy_yz, g_x_y_yy_zz, g_xyy_y_yy_xx, g_xyy_y_yy_xy, g_xyy_y_yy_xz, g_xyy_y_yy_yy, g_xyy_y_yy_yz, g_xyy_y_yy_zz, g_yy_0_0_0_x_y_yy_xx, g_yy_0_0_0_x_y_yy_xy, g_yy_0_0_0_x_y_yy_xz, g_yy_0_0_0_x_y_yy_yy, g_yy_0_0_0_x_y_yy_yz, g_yy_0_0_0_x_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_y_yy_xx[i] = -2.0 * g_x_y_yy_xx[i] * a_exp + 4.0 * g_xyy_y_yy_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_yy_xy[i] = -2.0 * g_x_y_yy_xy[i] * a_exp + 4.0 * g_xyy_y_yy_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_yy_xz[i] = -2.0 * g_x_y_yy_xz[i] * a_exp + 4.0 * g_xyy_y_yy_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_yy_yy[i] = -2.0 * g_x_y_yy_yy[i] * a_exp + 4.0 * g_xyy_y_yy_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_yy_yz[i] = -2.0 * g_x_y_yy_yz[i] * a_exp + 4.0 * g_xyy_y_yy_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_yy_zz[i] = -2.0 * g_x_y_yy_zz[i] * a_exp + 4.0 * g_xyy_y_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1032-1038)

    #pragma omp simd aligned(g_x_y_yz_xx, g_x_y_yz_xy, g_x_y_yz_xz, g_x_y_yz_yy, g_x_y_yz_yz, g_x_y_yz_zz, g_xyy_y_yz_xx, g_xyy_y_yz_xy, g_xyy_y_yz_xz, g_xyy_y_yz_yy, g_xyy_y_yz_yz, g_xyy_y_yz_zz, g_yy_0_0_0_x_y_yz_xx, g_yy_0_0_0_x_y_yz_xy, g_yy_0_0_0_x_y_yz_xz, g_yy_0_0_0_x_y_yz_yy, g_yy_0_0_0_x_y_yz_yz, g_yy_0_0_0_x_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_y_yz_xx[i] = -2.0 * g_x_y_yz_xx[i] * a_exp + 4.0 * g_xyy_y_yz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_yz_xy[i] = -2.0 * g_x_y_yz_xy[i] * a_exp + 4.0 * g_xyy_y_yz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_yz_xz[i] = -2.0 * g_x_y_yz_xz[i] * a_exp + 4.0 * g_xyy_y_yz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_yz_yy[i] = -2.0 * g_x_y_yz_yy[i] * a_exp + 4.0 * g_xyy_y_yz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_yz_yz[i] = -2.0 * g_x_y_yz_yz[i] * a_exp + 4.0 * g_xyy_y_yz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_yz_zz[i] = -2.0 * g_x_y_yz_zz[i] * a_exp + 4.0 * g_xyy_y_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1038-1044)

    #pragma omp simd aligned(g_x_y_zz_xx, g_x_y_zz_xy, g_x_y_zz_xz, g_x_y_zz_yy, g_x_y_zz_yz, g_x_y_zz_zz, g_xyy_y_zz_xx, g_xyy_y_zz_xy, g_xyy_y_zz_xz, g_xyy_y_zz_yy, g_xyy_y_zz_yz, g_xyy_y_zz_zz, g_yy_0_0_0_x_y_zz_xx, g_yy_0_0_0_x_y_zz_xy, g_yy_0_0_0_x_y_zz_xz, g_yy_0_0_0_x_y_zz_yy, g_yy_0_0_0_x_y_zz_yz, g_yy_0_0_0_x_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_y_zz_xx[i] = -2.0 * g_x_y_zz_xx[i] * a_exp + 4.0 * g_xyy_y_zz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_zz_xy[i] = -2.0 * g_x_y_zz_xy[i] * a_exp + 4.0 * g_xyy_y_zz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_zz_xz[i] = -2.0 * g_x_y_zz_xz[i] * a_exp + 4.0 * g_xyy_y_zz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_zz_yy[i] = -2.0 * g_x_y_zz_yy[i] * a_exp + 4.0 * g_xyy_y_zz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_zz_yz[i] = -2.0 * g_x_y_zz_yz[i] * a_exp + 4.0 * g_xyy_y_zz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_zz_zz[i] = -2.0 * g_x_y_zz_zz[i] * a_exp + 4.0 * g_xyy_y_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1044-1050)

    #pragma omp simd aligned(g_x_z_xx_xx, g_x_z_xx_xy, g_x_z_xx_xz, g_x_z_xx_yy, g_x_z_xx_yz, g_x_z_xx_zz, g_xyy_z_xx_xx, g_xyy_z_xx_xy, g_xyy_z_xx_xz, g_xyy_z_xx_yy, g_xyy_z_xx_yz, g_xyy_z_xx_zz, g_yy_0_0_0_x_z_xx_xx, g_yy_0_0_0_x_z_xx_xy, g_yy_0_0_0_x_z_xx_xz, g_yy_0_0_0_x_z_xx_yy, g_yy_0_0_0_x_z_xx_yz, g_yy_0_0_0_x_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_z_xx_xx[i] = -2.0 * g_x_z_xx_xx[i] * a_exp + 4.0 * g_xyy_z_xx_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_xx_xy[i] = -2.0 * g_x_z_xx_xy[i] * a_exp + 4.0 * g_xyy_z_xx_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_xx_xz[i] = -2.0 * g_x_z_xx_xz[i] * a_exp + 4.0 * g_xyy_z_xx_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_xx_yy[i] = -2.0 * g_x_z_xx_yy[i] * a_exp + 4.0 * g_xyy_z_xx_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_xx_yz[i] = -2.0 * g_x_z_xx_yz[i] * a_exp + 4.0 * g_xyy_z_xx_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_xx_zz[i] = -2.0 * g_x_z_xx_zz[i] * a_exp + 4.0 * g_xyy_z_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (1050-1056)

    #pragma omp simd aligned(g_x_z_xy_xx, g_x_z_xy_xy, g_x_z_xy_xz, g_x_z_xy_yy, g_x_z_xy_yz, g_x_z_xy_zz, g_xyy_z_xy_xx, g_xyy_z_xy_xy, g_xyy_z_xy_xz, g_xyy_z_xy_yy, g_xyy_z_xy_yz, g_xyy_z_xy_zz, g_yy_0_0_0_x_z_xy_xx, g_yy_0_0_0_x_z_xy_xy, g_yy_0_0_0_x_z_xy_xz, g_yy_0_0_0_x_z_xy_yy, g_yy_0_0_0_x_z_xy_yz, g_yy_0_0_0_x_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_z_xy_xx[i] = -2.0 * g_x_z_xy_xx[i] * a_exp + 4.0 * g_xyy_z_xy_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_xy_xy[i] = -2.0 * g_x_z_xy_xy[i] * a_exp + 4.0 * g_xyy_z_xy_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_xy_xz[i] = -2.0 * g_x_z_xy_xz[i] * a_exp + 4.0 * g_xyy_z_xy_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_xy_yy[i] = -2.0 * g_x_z_xy_yy[i] * a_exp + 4.0 * g_xyy_z_xy_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_xy_yz[i] = -2.0 * g_x_z_xy_yz[i] * a_exp + 4.0 * g_xyy_z_xy_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_xy_zz[i] = -2.0 * g_x_z_xy_zz[i] * a_exp + 4.0 * g_xyy_z_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1056-1062)

    #pragma omp simd aligned(g_x_z_xz_xx, g_x_z_xz_xy, g_x_z_xz_xz, g_x_z_xz_yy, g_x_z_xz_yz, g_x_z_xz_zz, g_xyy_z_xz_xx, g_xyy_z_xz_xy, g_xyy_z_xz_xz, g_xyy_z_xz_yy, g_xyy_z_xz_yz, g_xyy_z_xz_zz, g_yy_0_0_0_x_z_xz_xx, g_yy_0_0_0_x_z_xz_xy, g_yy_0_0_0_x_z_xz_xz, g_yy_0_0_0_x_z_xz_yy, g_yy_0_0_0_x_z_xz_yz, g_yy_0_0_0_x_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_z_xz_xx[i] = -2.0 * g_x_z_xz_xx[i] * a_exp + 4.0 * g_xyy_z_xz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_xz_xy[i] = -2.0 * g_x_z_xz_xy[i] * a_exp + 4.0 * g_xyy_z_xz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_xz_xz[i] = -2.0 * g_x_z_xz_xz[i] * a_exp + 4.0 * g_xyy_z_xz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_xz_yy[i] = -2.0 * g_x_z_xz_yy[i] * a_exp + 4.0 * g_xyy_z_xz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_xz_yz[i] = -2.0 * g_x_z_xz_yz[i] * a_exp + 4.0 * g_xyy_z_xz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_xz_zz[i] = -2.0 * g_x_z_xz_zz[i] * a_exp + 4.0 * g_xyy_z_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1062-1068)

    #pragma omp simd aligned(g_x_z_yy_xx, g_x_z_yy_xy, g_x_z_yy_xz, g_x_z_yy_yy, g_x_z_yy_yz, g_x_z_yy_zz, g_xyy_z_yy_xx, g_xyy_z_yy_xy, g_xyy_z_yy_xz, g_xyy_z_yy_yy, g_xyy_z_yy_yz, g_xyy_z_yy_zz, g_yy_0_0_0_x_z_yy_xx, g_yy_0_0_0_x_z_yy_xy, g_yy_0_0_0_x_z_yy_xz, g_yy_0_0_0_x_z_yy_yy, g_yy_0_0_0_x_z_yy_yz, g_yy_0_0_0_x_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_z_yy_xx[i] = -2.0 * g_x_z_yy_xx[i] * a_exp + 4.0 * g_xyy_z_yy_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_yy_xy[i] = -2.0 * g_x_z_yy_xy[i] * a_exp + 4.0 * g_xyy_z_yy_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_yy_xz[i] = -2.0 * g_x_z_yy_xz[i] * a_exp + 4.0 * g_xyy_z_yy_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_yy_yy[i] = -2.0 * g_x_z_yy_yy[i] * a_exp + 4.0 * g_xyy_z_yy_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_yy_yz[i] = -2.0 * g_x_z_yy_yz[i] * a_exp + 4.0 * g_xyy_z_yy_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_yy_zz[i] = -2.0 * g_x_z_yy_zz[i] * a_exp + 4.0 * g_xyy_z_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1068-1074)

    #pragma omp simd aligned(g_x_z_yz_xx, g_x_z_yz_xy, g_x_z_yz_xz, g_x_z_yz_yy, g_x_z_yz_yz, g_x_z_yz_zz, g_xyy_z_yz_xx, g_xyy_z_yz_xy, g_xyy_z_yz_xz, g_xyy_z_yz_yy, g_xyy_z_yz_yz, g_xyy_z_yz_zz, g_yy_0_0_0_x_z_yz_xx, g_yy_0_0_0_x_z_yz_xy, g_yy_0_0_0_x_z_yz_xz, g_yy_0_0_0_x_z_yz_yy, g_yy_0_0_0_x_z_yz_yz, g_yy_0_0_0_x_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_z_yz_xx[i] = -2.0 * g_x_z_yz_xx[i] * a_exp + 4.0 * g_xyy_z_yz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_yz_xy[i] = -2.0 * g_x_z_yz_xy[i] * a_exp + 4.0 * g_xyy_z_yz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_yz_xz[i] = -2.0 * g_x_z_yz_xz[i] * a_exp + 4.0 * g_xyy_z_yz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_yz_yy[i] = -2.0 * g_x_z_yz_yy[i] * a_exp + 4.0 * g_xyy_z_yz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_yz_yz[i] = -2.0 * g_x_z_yz_yz[i] * a_exp + 4.0 * g_xyy_z_yz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_yz_zz[i] = -2.0 * g_x_z_yz_zz[i] * a_exp + 4.0 * g_xyy_z_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1074-1080)

    #pragma omp simd aligned(g_x_z_zz_xx, g_x_z_zz_xy, g_x_z_zz_xz, g_x_z_zz_yy, g_x_z_zz_yz, g_x_z_zz_zz, g_xyy_z_zz_xx, g_xyy_z_zz_xy, g_xyy_z_zz_xz, g_xyy_z_zz_yy, g_xyy_z_zz_yz, g_xyy_z_zz_zz, g_yy_0_0_0_x_z_zz_xx, g_yy_0_0_0_x_z_zz_xy, g_yy_0_0_0_x_z_zz_xz, g_yy_0_0_0_x_z_zz_yy, g_yy_0_0_0_x_z_zz_yz, g_yy_0_0_0_x_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_z_zz_xx[i] = -2.0 * g_x_z_zz_xx[i] * a_exp + 4.0 * g_xyy_z_zz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_zz_xy[i] = -2.0 * g_x_z_zz_xy[i] * a_exp + 4.0 * g_xyy_z_zz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_zz_xz[i] = -2.0 * g_x_z_zz_xz[i] * a_exp + 4.0 * g_xyy_z_zz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_zz_yy[i] = -2.0 * g_x_z_zz_yy[i] * a_exp + 4.0 * g_xyy_z_zz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_zz_yz[i] = -2.0 * g_x_z_zz_yz[i] * a_exp + 4.0 * g_xyy_z_zz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_zz_zz[i] = -2.0 * g_x_z_zz_zz[i] * a_exp + 4.0 * g_xyy_z_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1080-1086)

    #pragma omp simd aligned(g_y_x_xx_xx, g_y_x_xx_xy, g_y_x_xx_xz, g_y_x_xx_yy, g_y_x_xx_yz, g_y_x_xx_zz, g_yy_0_0_0_y_x_xx_xx, g_yy_0_0_0_y_x_xx_xy, g_yy_0_0_0_y_x_xx_xz, g_yy_0_0_0_y_x_xx_yy, g_yy_0_0_0_y_x_xx_yz, g_yy_0_0_0_y_x_xx_zz, g_yyy_x_xx_xx, g_yyy_x_xx_xy, g_yyy_x_xx_xz, g_yyy_x_xx_yy, g_yyy_x_xx_yz, g_yyy_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_x_xx_xx[i] = -6.0 * g_y_x_xx_xx[i] * a_exp + 4.0 * g_yyy_x_xx_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_xx_xy[i] = -6.0 * g_y_x_xx_xy[i] * a_exp + 4.0 * g_yyy_x_xx_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_xx_xz[i] = -6.0 * g_y_x_xx_xz[i] * a_exp + 4.0 * g_yyy_x_xx_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_xx_yy[i] = -6.0 * g_y_x_xx_yy[i] * a_exp + 4.0 * g_yyy_x_xx_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_xx_yz[i] = -6.0 * g_y_x_xx_yz[i] * a_exp + 4.0 * g_yyy_x_xx_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_xx_zz[i] = -6.0 * g_y_x_xx_zz[i] * a_exp + 4.0 * g_yyy_x_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (1086-1092)

    #pragma omp simd aligned(g_y_x_xy_xx, g_y_x_xy_xy, g_y_x_xy_xz, g_y_x_xy_yy, g_y_x_xy_yz, g_y_x_xy_zz, g_yy_0_0_0_y_x_xy_xx, g_yy_0_0_0_y_x_xy_xy, g_yy_0_0_0_y_x_xy_xz, g_yy_0_0_0_y_x_xy_yy, g_yy_0_0_0_y_x_xy_yz, g_yy_0_0_0_y_x_xy_zz, g_yyy_x_xy_xx, g_yyy_x_xy_xy, g_yyy_x_xy_xz, g_yyy_x_xy_yy, g_yyy_x_xy_yz, g_yyy_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_x_xy_xx[i] = -6.0 * g_y_x_xy_xx[i] * a_exp + 4.0 * g_yyy_x_xy_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_xy_xy[i] = -6.0 * g_y_x_xy_xy[i] * a_exp + 4.0 * g_yyy_x_xy_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_xy_xz[i] = -6.0 * g_y_x_xy_xz[i] * a_exp + 4.0 * g_yyy_x_xy_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_xy_yy[i] = -6.0 * g_y_x_xy_yy[i] * a_exp + 4.0 * g_yyy_x_xy_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_xy_yz[i] = -6.0 * g_y_x_xy_yz[i] * a_exp + 4.0 * g_yyy_x_xy_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_xy_zz[i] = -6.0 * g_y_x_xy_zz[i] * a_exp + 4.0 * g_yyy_x_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1092-1098)

    #pragma omp simd aligned(g_y_x_xz_xx, g_y_x_xz_xy, g_y_x_xz_xz, g_y_x_xz_yy, g_y_x_xz_yz, g_y_x_xz_zz, g_yy_0_0_0_y_x_xz_xx, g_yy_0_0_0_y_x_xz_xy, g_yy_0_0_0_y_x_xz_xz, g_yy_0_0_0_y_x_xz_yy, g_yy_0_0_0_y_x_xz_yz, g_yy_0_0_0_y_x_xz_zz, g_yyy_x_xz_xx, g_yyy_x_xz_xy, g_yyy_x_xz_xz, g_yyy_x_xz_yy, g_yyy_x_xz_yz, g_yyy_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_x_xz_xx[i] = -6.0 * g_y_x_xz_xx[i] * a_exp + 4.0 * g_yyy_x_xz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_xz_xy[i] = -6.0 * g_y_x_xz_xy[i] * a_exp + 4.0 * g_yyy_x_xz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_xz_xz[i] = -6.0 * g_y_x_xz_xz[i] * a_exp + 4.0 * g_yyy_x_xz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_xz_yy[i] = -6.0 * g_y_x_xz_yy[i] * a_exp + 4.0 * g_yyy_x_xz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_xz_yz[i] = -6.0 * g_y_x_xz_yz[i] * a_exp + 4.0 * g_yyy_x_xz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_xz_zz[i] = -6.0 * g_y_x_xz_zz[i] * a_exp + 4.0 * g_yyy_x_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1098-1104)

    #pragma omp simd aligned(g_y_x_yy_xx, g_y_x_yy_xy, g_y_x_yy_xz, g_y_x_yy_yy, g_y_x_yy_yz, g_y_x_yy_zz, g_yy_0_0_0_y_x_yy_xx, g_yy_0_0_0_y_x_yy_xy, g_yy_0_0_0_y_x_yy_xz, g_yy_0_0_0_y_x_yy_yy, g_yy_0_0_0_y_x_yy_yz, g_yy_0_0_0_y_x_yy_zz, g_yyy_x_yy_xx, g_yyy_x_yy_xy, g_yyy_x_yy_xz, g_yyy_x_yy_yy, g_yyy_x_yy_yz, g_yyy_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_x_yy_xx[i] = -6.0 * g_y_x_yy_xx[i] * a_exp + 4.0 * g_yyy_x_yy_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_yy_xy[i] = -6.0 * g_y_x_yy_xy[i] * a_exp + 4.0 * g_yyy_x_yy_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_yy_xz[i] = -6.0 * g_y_x_yy_xz[i] * a_exp + 4.0 * g_yyy_x_yy_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_yy_yy[i] = -6.0 * g_y_x_yy_yy[i] * a_exp + 4.0 * g_yyy_x_yy_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_yy_yz[i] = -6.0 * g_y_x_yy_yz[i] * a_exp + 4.0 * g_yyy_x_yy_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_yy_zz[i] = -6.0 * g_y_x_yy_zz[i] * a_exp + 4.0 * g_yyy_x_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1104-1110)

    #pragma omp simd aligned(g_y_x_yz_xx, g_y_x_yz_xy, g_y_x_yz_xz, g_y_x_yz_yy, g_y_x_yz_yz, g_y_x_yz_zz, g_yy_0_0_0_y_x_yz_xx, g_yy_0_0_0_y_x_yz_xy, g_yy_0_0_0_y_x_yz_xz, g_yy_0_0_0_y_x_yz_yy, g_yy_0_0_0_y_x_yz_yz, g_yy_0_0_0_y_x_yz_zz, g_yyy_x_yz_xx, g_yyy_x_yz_xy, g_yyy_x_yz_xz, g_yyy_x_yz_yy, g_yyy_x_yz_yz, g_yyy_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_x_yz_xx[i] = -6.0 * g_y_x_yz_xx[i] * a_exp + 4.0 * g_yyy_x_yz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_yz_xy[i] = -6.0 * g_y_x_yz_xy[i] * a_exp + 4.0 * g_yyy_x_yz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_yz_xz[i] = -6.0 * g_y_x_yz_xz[i] * a_exp + 4.0 * g_yyy_x_yz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_yz_yy[i] = -6.0 * g_y_x_yz_yy[i] * a_exp + 4.0 * g_yyy_x_yz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_yz_yz[i] = -6.0 * g_y_x_yz_yz[i] * a_exp + 4.0 * g_yyy_x_yz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_yz_zz[i] = -6.0 * g_y_x_yz_zz[i] * a_exp + 4.0 * g_yyy_x_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1110-1116)

    #pragma omp simd aligned(g_y_x_zz_xx, g_y_x_zz_xy, g_y_x_zz_xz, g_y_x_zz_yy, g_y_x_zz_yz, g_y_x_zz_zz, g_yy_0_0_0_y_x_zz_xx, g_yy_0_0_0_y_x_zz_xy, g_yy_0_0_0_y_x_zz_xz, g_yy_0_0_0_y_x_zz_yy, g_yy_0_0_0_y_x_zz_yz, g_yy_0_0_0_y_x_zz_zz, g_yyy_x_zz_xx, g_yyy_x_zz_xy, g_yyy_x_zz_xz, g_yyy_x_zz_yy, g_yyy_x_zz_yz, g_yyy_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_x_zz_xx[i] = -6.0 * g_y_x_zz_xx[i] * a_exp + 4.0 * g_yyy_x_zz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_zz_xy[i] = -6.0 * g_y_x_zz_xy[i] * a_exp + 4.0 * g_yyy_x_zz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_zz_xz[i] = -6.0 * g_y_x_zz_xz[i] * a_exp + 4.0 * g_yyy_x_zz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_zz_yy[i] = -6.0 * g_y_x_zz_yy[i] * a_exp + 4.0 * g_yyy_x_zz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_zz_yz[i] = -6.0 * g_y_x_zz_yz[i] * a_exp + 4.0 * g_yyy_x_zz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_zz_zz[i] = -6.0 * g_y_x_zz_zz[i] * a_exp + 4.0 * g_yyy_x_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1116-1122)

    #pragma omp simd aligned(g_y_y_xx_xx, g_y_y_xx_xy, g_y_y_xx_xz, g_y_y_xx_yy, g_y_y_xx_yz, g_y_y_xx_zz, g_yy_0_0_0_y_y_xx_xx, g_yy_0_0_0_y_y_xx_xy, g_yy_0_0_0_y_y_xx_xz, g_yy_0_0_0_y_y_xx_yy, g_yy_0_0_0_y_y_xx_yz, g_yy_0_0_0_y_y_xx_zz, g_yyy_y_xx_xx, g_yyy_y_xx_xy, g_yyy_y_xx_xz, g_yyy_y_xx_yy, g_yyy_y_xx_yz, g_yyy_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_y_xx_xx[i] = -6.0 * g_y_y_xx_xx[i] * a_exp + 4.0 * g_yyy_y_xx_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_xx_xy[i] = -6.0 * g_y_y_xx_xy[i] * a_exp + 4.0 * g_yyy_y_xx_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_xx_xz[i] = -6.0 * g_y_y_xx_xz[i] * a_exp + 4.0 * g_yyy_y_xx_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_xx_yy[i] = -6.0 * g_y_y_xx_yy[i] * a_exp + 4.0 * g_yyy_y_xx_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_xx_yz[i] = -6.0 * g_y_y_xx_yz[i] * a_exp + 4.0 * g_yyy_y_xx_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_xx_zz[i] = -6.0 * g_y_y_xx_zz[i] * a_exp + 4.0 * g_yyy_y_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (1122-1128)

    #pragma omp simd aligned(g_y_y_xy_xx, g_y_y_xy_xy, g_y_y_xy_xz, g_y_y_xy_yy, g_y_y_xy_yz, g_y_y_xy_zz, g_yy_0_0_0_y_y_xy_xx, g_yy_0_0_0_y_y_xy_xy, g_yy_0_0_0_y_y_xy_xz, g_yy_0_0_0_y_y_xy_yy, g_yy_0_0_0_y_y_xy_yz, g_yy_0_0_0_y_y_xy_zz, g_yyy_y_xy_xx, g_yyy_y_xy_xy, g_yyy_y_xy_xz, g_yyy_y_xy_yy, g_yyy_y_xy_yz, g_yyy_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_y_xy_xx[i] = -6.0 * g_y_y_xy_xx[i] * a_exp + 4.0 * g_yyy_y_xy_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_xy_xy[i] = -6.0 * g_y_y_xy_xy[i] * a_exp + 4.0 * g_yyy_y_xy_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_xy_xz[i] = -6.0 * g_y_y_xy_xz[i] * a_exp + 4.0 * g_yyy_y_xy_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_xy_yy[i] = -6.0 * g_y_y_xy_yy[i] * a_exp + 4.0 * g_yyy_y_xy_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_xy_yz[i] = -6.0 * g_y_y_xy_yz[i] * a_exp + 4.0 * g_yyy_y_xy_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_xy_zz[i] = -6.0 * g_y_y_xy_zz[i] * a_exp + 4.0 * g_yyy_y_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1128-1134)

    #pragma omp simd aligned(g_y_y_xz_xx, g_y_y_xz_xy, g_y_y_xz_xz, g_y_y_xz_yy, g_y_y_xz_yz, g_y_y_xz_zz, g_yy_0_0_0_y_y_xz_xx, g_yy_0_0_0_y_y_xz_xy, g_yy_0_0_0_y_y_xz_xz, g_yy_0_0_0_y_y_xz_yy, g_yy_0_0_0_y_y_xz_yz, g_yy_0_0_0_y_y_xz_zz, g_yyy_y_xz_xx, g_yyy_y_xz_xy, g_yyy_y_xz_xz, g_yyy_y_xz_yy, g_yyy_y_xz_yz, g_yyy_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_y_xz_xx[i] = -6.0 * g_y_y_xz_xx[i] * a_exp + 4.0 * g_yyy_y_xz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_xz_xy[i] = -6.0 * g_y_y_xz_xy[i] * a_exp + 4.0 * g_yyy_y_xz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_xz_xz[i] = -6.0 * g_y_y_xz_xz[i] * a_exp + 4.0 * g_yyy_y_xz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_xz_yy[i] = -6.0 * g_y_y_xz_yy[i] * a_exp + 4.0 * g_yyy_y_xz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_xz_yz[i] = -6.0 * g_y_y_xz_yz[i] * a_exp + 4.0 * g_yyy_y_xz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_xz_zz[i] = -6.0 * g_y_y_xz_zz[i] * a_exp + 4.0 * g_yyy_y_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1134-1140)

    #pragma omp simd aligned(g_y_y_yy_xx, g_y_y_yy_xy, g_y_y_yy_xz, g_y_y_yy_yy, g_y_y_yy_yz, g_y_y_yy_zz, g_yy_0_0_0_y_y_yy_xx, g_yy_0_0_0_y_y_yy_xy, g_yy_0_0_0_y_y_yy_xz, g_yy_0_0_0_y_y_yy_yy, g_yy_0_0_0_y_y_yy_yz, g_yy_0_0_0_y_y_yy_zz, g_yyy_y_yy_xx, g_yyy_y_yy_xy, g_yyy_y_yy_xz, g_yyy_y_yy_yy, g_yyy_y_yy_yz, g_yyy_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_y_yy_xx[i] = -6.0 * g_y_y_yy_xx[i] * a_exp + 4.0 * g_yyy_y_yy_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_yy_xy[i] = -6.0 * g_y_y_yy_xy[i] * a_exp + 4.0 * g_yyy_y_yy_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_yy_xz[i] = -6.0 * g_y_y_yy_xz[i] * a_exp + 4.0 * g_yyy_y_yy_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_yy_yy[i] = -6.0 * g_y_y_yy_yy[i] * a_exp + 4.0 * g_yyy_y_yy_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_yy_yz[i] = -6.0 * g_y_y_yy_yz[i] * a_exp + 4.0 * g_yyy_y_yy_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_yy_zz[i] = -6.0 * g_y_y_yy_zz[i] * a_exp + 4.0 * g_yyy_y_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1140-1146)

    #pragma omp simd aligned(g_y_y_yz_xx, g_y_y_yz_xy, g_y_y_yz_xz, g_y_y_yz_yy, g_y_y_yz_yz, g_y_y_yz_zz, g_yy_0_0_0_y_y_yz_xx, g_yy_0_0_0_y_y_yz_xy, g_yy_0_0_0_y_y_yz_xz, g_yy_0_0_0_y_y_yz_yy, g_yy_0_0_0_y_y_yz_yz, g_yy_0_0_0_y_y_yz_zz, g_yyy_y_yz_xx, g_yyy_y_yz_xy, g_yyy_y_yz_xz, g_yyy_y_yz_yy, g_yyy_y_yz_yz, g_yyy_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_y_yz_xx[i] = -6.0 * g_y_y_yz_xx[i] * a_exp + 4.0 * g_yyy_y_yz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_yz_xy[i] = -6.0 * g_y_y_yz_xy[i] * a_exp + 4.0 * g_yyy_y_yz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_yz_xz[i] = -6.0 * g_y_y_yz_xz[i] * a_exp + 4.0 * g_yyy_y_yz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_yz_yy[i] = -6.0 * g_y_y_yz_yy[i] * a_exp + 4.0 * g_yyy_y_yz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_yz_yz[i] = -6.0 * g_y_y_yz_yz[i] * a_exp + 4.0 * g_yyy_y_yz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_yz_zz[i] = -6.0 * g_y_y_yz_zz[i] * a_exp + 4.0 * g_yyy_y_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1146-1152)

    #pragma omp simd aligned(g_y_y_zz_xx, g_y_y_zz_xy, g_y_y_zz_xz, g_y_y_zz_yy, g_y_y_zz_yz, g_y_y_zz_zz, g_yy_0_0_0_y_y_zz_xx, g_yy_0_0_0_y_y_zz_xy, g_yy_0_0_0_y_y_zz_xz, g_yy_0_0_0_y_y_zz_yy, g_yy_0_0_0_y_y_zz_yz, g_yy_0_0_0_y_y_zz_zz, g_yyy_y_zz_xx, g_yyy_y_zz_xy, g_yyy_y_zz_xz, g_yyy_y_zz_yy, g_yyy_y_zz_yz, g_yyy_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_y_zz_xx[i] = -6.0 * g_y_y_zz_xx[i] * a_exp + 4.0 * g_yyy_y_zz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_zz_xy[i] = -6.0 * g_y_y_zz_xy[i] * a_exp + 4.0 * g_yyy_y_zz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_zz_xz[i] = -6.0 * g_y_y_zz_xz[i] * a_exp + 4.0 * g_yyy_y_zz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_zz_yy[i] = -6.0 * g_y_y_zz_yy[i] * a_exp + 4.0 * g_yyy_y_zz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_zz_yz[i] = -6.0 * g_y_y_zz_yz[i] * a_exp + 4.0 * g_yyy_y_zz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_zz_zz[i] = -6.0 * g_y_y_zz_zz[i] * a_exp + 4.0 * g_yyy_y_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1152-1158)

    #pragma omp simd aligned(g_y_z_xx_xx, g_y_z_xx_xy, g_y_z_xx_xz, g_y_z_xx_yy, g_y_z_xx_yz, g_y_z_xx_zz, g_yy_0_0_0_y_z_xx_xx, g_yy_0_0_0_y_z_xx_xy, g_yy_0_0_0_y_z_xx_xz, g_yy_0_0_0_y_z_xx_yy, g_yy_0_0_0_y_z_xx_yz, g_yy_0_0_0_y_z_xx_zz, g_yyy_z_xx_xx, g_yyy_z_xx_xy, g_yyy_z_xx_xz, g_yyy_z_xx_yy, g_yyy_z_xx_yz, g_yyy_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_z_xx_xx[i] = -6.0 * g_y_z_xx_xx[i] * a_exp + 4.0 * g_yyy_z_xx_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_xx_xy[i] = -6.0 * g_y_z_xx_xy[i] * a_exp + 4.0 * g_yyy_z_xx_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_xx_xz[i] = -6.0 * g_y_z_xx_xz[i] * a_exp + 4.0 * g_yyy_z_xx_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_xx_yy[i] = -6.0 * g_y_z_xx_yy[i] * a_exp + 4.0 * g_yyy_z_xx_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_xx_yz[i] = -6.0 * g_y_z_xx_yz[i] * a_exp + 4.0 * g_yyy_z_xx_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_xx_zz[i] = -6.0 * g_y_z_xx_zz[i] * a_exp + 4.0 * g_yyy_z_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (1158-1164)

    #pragma omp simd aligned(g_y_z_xy_xx, g_y_z_xy_xy, g_y_z_xy_xz, g_y_z_xy_yy, g_y_z_xy_yz, g_y_z_xy_zz, g_yy_0_0_0_y_z_xy_xx, g_yy_0_0_0_y_z_xy_xy, g_yy_0_0_0_y_z_xy_xz, g_yy_0_0_0_y_z_xy_yy, g_yy_0_0_0_y_z_xy_yz, g_yy_0_0_0_y_z_xy_zz, g_yyy_z_xy_xx, g_yyy_z_xy_xy, g_yyy_z_xy_xz, g_yyy_z_xy_yy, g_yyy_z_xy_yz, g_yyy_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_z_xy_xx[i] = -6.0 * g_y_z_xy_xx[i] * a_exp + 4.0 * g_yyy_z_xy_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_xy_xy[i] = -6.0 * g_y_z_xy_xy[i] * a_exp + 4.0 * g_yyy_z_xy_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_xy_xz[i] = -6.0 * g_y_z_xy_xz[i] * a_exp + 4.0 * g_yyy_z_xy_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_xy_yy[i] = -6.0 * g_y_z_xy_yy[i] * a_exp + 4.0 * g_yyy_z_xy_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_xy_yz[i] = -6.0 * g_y_z_xy_yz[i] * a_exp + 4.0 * g_yyy_z_xy_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_xy_zz[i] = -6.0 * g_y_z_xy_zz[i] * a_exp + 4.0 * g_yyy_z_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1164-1170)

    #pragma omp simd aligned(g_y_z_xz_xx, g_y_z_xz_xy, g_y_z_xz_xz, g_y_z_xz_yy, g_y_z_xz_yz, g_y_z_xz_zz, g_yy_0_0_0_y_z_xz_xx, g_yy_0_0_0_y_z_xz_xy, g_yy_0_0_0_y_z_xz_xz, g_yy_0_0_0_y_z_xz_yy, g_yy_0_0_0_y_z_xz_yz, g_yy_0_0_0_y_z_xz_zz, g_yyy_z_xz_xx, g_yyy_z_xz_xy, g_yyy_z_xz_xz, g_yyy_z_xz_yy, g_yyy_z_xz_yz, g_yyy_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_z_xz_xx[i] = -6.0 * g_y_z_xz_xx[i] * a_exp + 4.0 * g_yyy_z_xz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_xz_xy[i] = -6.0 * g_y_z_xz_xy[i] * a_exp + 4.0 * g_yyy_z_xz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_xz_xz[i] = -6.0 * g_y_z_xz_xz[i] * a_exp + 4.0 * g_yyy_z_xz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_xz_yy[i] = -6.0 * g_y_z_xz_yy[i] * a_exp + 4.0 * g_yyy_z_xz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_xz_yz[i] = -6.0 * g_y_z_xz_yz[i] * a_exp + 4.0 * g_yyy_z_xz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_xz_zz[i] = -6.0 * g_y_z_xz_zz[i] * a_exp + 4.0 * g_yyy_z_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1170-1176)

    #pragma omp simd aligned(g_y_z_yy_xx, g_y_z_yy_xy, g_y_z_yy_xz, g_y_z_yy_yy, g_y_z_yy_yz, g_y_z_yy_zz, g_yy_0_0_0_y_z_yy_xx, g_yy_0_0_0_y_z_yy_xy, g_yy_0_0_0_y_z_yy_xz, g_yy_0_0_0_y_z_yy_yy, g_yy_0_0_0_y_z_yy_yz, g_yy_0_0_0_y_z_yy_zz, g_yyy_z_yy_xx, g_yyy_z_yy_xy, g_yyy_z_yy_xz, g_yyy_z_yy_yy, g_yyy_z_yy_yz, g_yyy_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_z_yy_xx[i] = -6.0 * g_y_z_yy_xx[i] * a_exp + 4.0 * g_yyy_z_yy_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_yy_xy[i] = -6.0 * g_y_z_yy_xy[i] * a_exp + 4.0 * g_yyy_z_yy_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_yy_xz[i] = -6.0 * g_y_z_yy_xz[i] * a_exp + 4.0 * g_yyy_z_yy_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_yy_yy[i] = -6.0 * g_y_z_yy_yy[i] * a_exp + 4.0 * g_yyy_z_yy_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_yy_yz[i] = -6.0 * g_y_z_yy_yz[i] * a_exp + 4.0 * g_yyy_z_yy_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_yy_zz[i] = -6.0 * g_y_z_yy_zz[i] * a_exp + 4.0 * g_yyy_z_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1176-1182)

    #pragma omp simd aligned(g_y_z_yz_xx, g_y_z_yz_xy, g_y_z_yz_xz, g_y_z_yz_yy, g_y_z_yz_yz, g_y_z_yz_zz, g_yy_0_0_0_y_z_yz_xx, g_yy_0_0_0_y_z_yz_xy, g_yy_0_0_0_y_z_yz_xz, g_yy_0_0_0_y_z_yz_yy, g_yy_0_0_0_y_z_yz_yz, g_yy_0_0_0_y_z_yz_zz, g_yyy_z_yz_xx, g_yyy_z_yz_xy, g_yyy_z_yz_xz, g_yyy_z_yz_yy, g_yyy_z_yz_yz, g_yyy_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_z_yz_xx[i] = -6.0 * g_y_z_yz_xx[i] * a_exp + 4.0 * g_yyy_z_yz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_yz_xy[i] = -6.0 * g_y_z_yz_xy[i] * a_exp + 4.0 * g_yyy_z_yz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_yz_xz[i] = -6.0 * g_y_z_yz_xz[i] * a_exp + 4.0 * g_yyy_z_yz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_yz_yy[i] = -6.0 * g_y_z_yz_yy[i] * a_exp + 4.0 * g_yyy_z_yz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_yz_yz[i] = -6.0 * g_y_z_yz_yz[i] * a_exp + 4.0 * g_yyy_z_yz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_yz_zz[i] = -6.0 * g_y_z_yz_zz[i] * a_exp + 4.0 * g_yyy_z_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1182-1188)

    #pragma omp simd aligned(g_y_z_zz_xx, g_y_z_zz_xy, g_y_z_zz_xz, g_y_z_zz_yy, g_y_z_zz_yz, g_y_z_zz_zz, g_yy_0_0_0_y_z_zz_xx, g_yy_0_0_0_y_z_zz_xy, g_yy_0_0_0_y_z_zz_xz, g_yy_0_0_0_y_z_zz_yy, g_yy_0_0_0_y_z_zz_yz, g_yy_0_0_0_y_z_zz_zz, g_yyy_z_zz_xx, g_yyy_z_zz_xy, g_yyy_z_zz_xz, g_yyy_z_zz_yy, g_yyy_z_zz_yz, g_yyy_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_z_zz_xx[i] = -6.0 * g_y_z_zz_xx[i] * a_exp + 4.0 * g_yyy_z_zz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_zz_xy[i] = -6.0 * g_y_z_zz_xy[i] * a_exp + 4.0 * g_yyy_z_zz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_zz_xz[i] = -6.0 * g_y_z_zz_xz[i] * a_exp + 4.0 * g_yyy_z_zz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_zz_yy[i] = -6.0 * g_y_z_zz_yy[i] * a_exp + 4.0 * g_yyy_z_zz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_zz_yz[i] = -6.0 * g_y_z_zz_yz[i] * a_exp + 4.0 * g_yyy_z_zz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_zz_zz[i] = -6.0 * g_y_z_zz_zz[i] * a_exp + 4.0 * g_yyy_z_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1188-1194)

    #pragma omp simd aligned(g_yy_0_0_0_z_x_xx_xx, g_yy_0_0_0_z_x_xx_xy, g_yy_0_0_0_z_x_xx_xz, g_yy_0_0_0_z_x_xx_yy, g_yy_0_0_0_z_x_xx_yz, g_yy_0_0_0_z_x_xx_zz, g_yyz_x_xx_xx, g_yyz_x_xx_xy, g_yyz_x_xx_xz, g_yyz_x_xx_yy, g_yyz_x_xx_yz, g_yyz_x_xx_zz, g_z_x_xx_xx, g_z_x_xx_xy, g_z_x_xx_xz, g_z_x_xx_yy, g_z_x_xx_yz, g_z_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_x_xx_xx[i] = -2.0 * g_z_x_xx_xx[i] * a_exp + 4.0 * g_yyz_x_xx_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_xx_xy[i] = -2.0 * g_z_x_xx_xy[i] * a_exp + 4.0 * g_yyz_x_xx_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_xx_xz[i] = -2.0 * g_z_x_xx_xz[i] * a_exp + 4.0 * g_yyz_x_xx_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_xx_yy[i] = -2.0 * g_z_x_xx_yy[i] * a_exp + 4.0 * g_yyz_x_xx_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_xx_yz[i] = -2.0 * g_z_x_xx_yz[i] * a_exp + 4.0 * g_yyz_x_xx_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_xx_zz[i] = -2.0 * g_z_x_xx_zz[i] * a_exp + 4.0 * g_yyz_x_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (1194-1200)

    #pragma omp simd aligned(g_yy_0_0_0_z_x_xy_xx, g_yy_0_0_0_z_x_xy_xy, g_yy_0_0_0_z_x_xy_xz, g_yy_0_0_0_z_x_xy_yy, g_yy_0_0_0_z_x_xy_yz, g_yy_0_0_0_z_x_xy_zz, g_yyz_x_xy_xx, g_yyz_x_xy_xy, g_yyz_x_xy_xz, g_yyz_x_xy_yy, g_yyz_x_xy_yz, g_yyz_x_xy_zz, g_z_x_xy_xx, g_z_x_xy_xy, g_z_x_xy_xz, g_z_x_xy_yy, g_z_x_xy_yz, g_z_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_x_xy_xx[i] = -2.0 * g_z_x_xy_xx[i] * a_exp + 4.0 * g_yyz_x_xy_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_xy_xy[i] = -2.0 * g_z_x_xy_xy[i] * a_exp + 4.0 * g_yyz_x_xy_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_xy_xz[i] = -2.0 * g_z_x_xy_xz[i] * a_exp + 4.0 * g_yyz_x_xy_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_xy_yy[i] = -2.0 * g_z_x_xy_yy[i] * a_exp + 4.0 * g_yyz_x_xy_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_xy_yz[i] = -2.0 * g_z_x_xy_yz[i] * a_exp + 4.0 * g_yyz_x_xy_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_xy_zz[i] = -2.0 * g_z_x_xy_zz[i] * a_exp + 4.0 * g_yyz_x_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1200-1206)

    #pragma omp simd aligned(g_yy_0_0_0_z_x_xz_xx, g_yy_0_0_0_z_x_xz_xy, g_yy_0_0_0_z_x_xz_xz, g_yy_0_0_0_z_x_xz_yy, g_yy_0_0_0_z_x_xz_yz, g_yy_0_0_0_z_x_xz_zz, g_yyz_x_xz_xx, g_yyz_x_xz_xy, g_yyz_x_xz_xz, g_yyz_x_xz_yy, g_yyz_x_xz_yz, g_yyz_x_xz_zz, g_z_x_xz_xx, g_z_x_xz_xy, g_z_x_xz_xz, g_z_x_xz_yy, g_z_x_xz_yz, g_z_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_x_xz_xx[i] = -2.0 * g_z_x_xz_xx[i] * a_exp + 4.0 * g_yyz_x_xz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_xz_xy[i] = -2.0 * g_z_x_xz_xy[i] * a_exp + 4.0 * g_yyz_x_xz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_xz_xz[i] = -2.0 * g_z_x_xz_xz[i] * a_exp + 4.0 * g_yyz_x_xz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_xz_yy[i] = -2.0 * g_z_x_xz_yy[i] * a_exp + 4.0 * g_yyz_x_xz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_xz_yz[i] = -2.0 * g_z_x_xz_yz[i] * a_exp + 4.0 * g_yyz_x_xz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_xz_zz[i] = -2.0 * g_z_x_xz_zz[i] * a_exp + 4.0 * g_yyz_x_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1206-1212)

    #pragma omp simd aligned(g_yy_0_0_0_z_x_yy_xx, g_yy_0_0_0_z_x_yy_xy, g_yy_0_0_0_z_x_yy_xz, g_yy_0_0_0_z_x_yy_yy, g_yy_0_0_0_z_x_yy_yz, g_yy_0_0_0_z_x_yy_zz, g_yyz_x_yy_xx, g_yyz_x_yy_xy, g_yyz_x_yy_xz, g_yyz_x_yy_yy, g_yyz_x_yy_yz, g_yyz_x_yy_zz, g_z_x_yy_xx, g_z_x_yy_xy, g_z_x_yy_xz, g_z_x_yy_yy, g_z_x_yy_yz, g_z_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_x_yy_xx[i] = -2.0 * g_z_x_yy_xx[i] * a_exp + 4.0 * g_yyz_x_yy_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_yy_xy[i] = -2.0 * g_z_x_yy_xy[i] * a_exp + 4.0 * g_yyz_x_yy_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_yy_xz[i] = -2.0 * g_z_x_yy_xz[i] * a_exp + 4.0 * g_yyz_x_yy_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_yy_yy[i] = -2.0 * g_z_x_yy_yy[i] * a_exp + 4.0 * g_yyz_x_yy_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_yy_yz[i] = -2.0 * g_z_x_yy_yz[i] * a_exp + 4.0 * g_yyz_x_yy_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_yy_zz[i] = -2.0 * g_z_x_yy_zz[i] * a_exp + 4.0 * g_yyz_x_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1212-1218)

    #pragma omp simd aligned(g_yy_0_0_0_z_x_yz_xx, g_yy_0_0_0_z_x_yz_xy, g_yy_0_0_0_z_x_yz_xz, g_yy_0_0_0_z_x_yz_yy, g_yy_0_0_0_z_x_yz_yz, g_yy_0_0_0_z_x_yz_zz, g_yyz_x_yz_xx, g_yyz_x_yz_xy, g_yyz_x_yz_xz, g_yyz_x_yz_yy, g_yyz_x_yz_yz, g_yyz_x_yz_zz, g_z_x_yz_xx, g_z_x_yz_xy, g_z_x_yz_xz, g_z_x_yz_yy, g_z_x_yz_yz, g_z_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_x_yz_xx[i] = -2.0 * g_z_x_yz_xx[i] * a_exp + 4.0 * g_yyz_x_yz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_yz_xy[i] = -2.0 * g_z_x_yz_xy[i] * a_exp + 4.0 * g_yyz_x_yz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_yz_xz[i] = -2.0 * g_z_x_yz_xz[i] * a_exp + 4.0 * g_yyz_x_yz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_yz_yy[i] = -2.0 * g_z_x_yz_yy[i] * a_exp + 4.0 * g_yyz_x_yz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_yz_yz[i] = -2.0 * g_z_x_yz_yz[i] * a_exp + 4.0 * g_yyz_x_yz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_yz_zz[i] = -2.0 * g_z_x_yz_zz[i] * a_exp + 4.0 * g_yyz_x_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1218-1224)

    #pragma omp simd aligned(g_yy_0_0_0_z_x_zz_xx, g_yy_0_0_0_z_x_zz_xy, g_yy_0_0_0_z_x_zz_xz, g_yy_0_0_0_z_x_zz_yy, g_yy_0_0_0_z_x_zz_yz, g_yy_0_0_0_z_x_zz_zz, g_yyz_x_zz_xx, g_yyz_x_zz_xy, g_yyz_x_zz_xz, g_yyz_x_zz_yy, g_yyz_x_zz_yz, g_yyz_x_zz_zz, g_z_x_zz_xx, g_z_x_zz_xy, g_z_x_zz_xz, g_z_x_zz_yy, g_z_x_zz_yz, g_z_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_x_zz_xx[i] = -2.0 * g_z_x_zz_xx[i] * a_exp + 4.0 * g_yyz_x_zz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_zz_xy[i] = -2.0 * g_z_x_zz_xy[i] * a_exp + 4.0 * g_yyz_x_zz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_zz_xz[i] = -2.0 * g_z_x_zz_xz[i] * a_exp + 4.0 * g_yyz_x_zz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_zz_yy[i] = -2.0 * g_z_x_zz_yy[i] * a_exp + 4.0 * g_yyz_x_zz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_zz_yz[i] = -2.0 * g_z_x_zz_yz[i] * a_exp + 4.0 * g_yyz_x_zz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_zz_zz[i] = -2.0 * g_z_x_zz_zz[i] * a_exp + 4.0 * g_yyz_x_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1224-1230)

    #pragma omp simd aligned(g_yy_0_0_0_z_y_xx_xx, g_yy_0_0_0_z_y_xx_xy, g_yy_0_0_0_z_y_xx_xz, g_yy_0_0_0_z_y_xx_yy, g_yy_0_0_0_z_y_xx_yz, g_yy_0_0_0_z_y_xx_zz, g_yyz_y_xx_xx, g_yyz_y_xx_xy, g_yyz_y_xx_xz, g_yyz_y_xx_yy, g_yyz_y_xx_yz, g_yyz_y_xx_zz, g_z_y_xx_xx, g_z_y_xx_xy, g_z_y_xx_xz, g_z_y_xx_yy, g_z_y_xx_yz, g_z_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_y_xx_xx[i] = -2.0 * g_z_y_xx_xx[i] * a_exp + 4.0 * g_yyz_y_xx_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_xx_xy[i] = -2.0 * g_z_y_xx_xy[i] * a_exp + 4.0 * g_yyz_y_xx_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_xx_xz[i] = -2.0 * g_z_y_xx_xz[i] * a_exp + 4.0 * g_yyz_y_xx_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_xx_yy[i] = -2.0 * g_z_y_xx_yy[i] * a_exp + 4.0 * g_yyz_y_xx_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_xx_yz[i] = -2.0 * g_z_y_xx_yz[i] * a_exp + 4.0 * g_yyz_y_xx_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_xx_zz[i] = -2.0 * g_z_y_xx_zz[i] * a_exp + 4.0 * g_yyz_y_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (1230-1236)

    #pragma omp simd aligned(g_yy_0_0_0_z_y_xy_xx, g_yy_0_0_0_z_y_xy_xy, g_yy_0_0_0_z_y_xy_xz, g_yy_0_0_0_z_y_xy_yy, g_yy_0_0_0_z_y_xy_yz, g_yy_0_0_0_z_y_xy_zz, g_yyz_y_xy_xx, g_yyz_y_xy_xy, g_yyz_y_xy_xz, g_yyz_y_xy_yy, g_yyz_y_xy_yz, g_yyz_y_xy_zz, g_z_y_xy_xx, g_z_y_xy_xy, g_z_y_xy_xz, g_z_y_xy_yy, g_z_y_xy_yz, g_z_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_y_xy_xx[i] = -2.0 * g_z_y_xy_xx[i] * a_exp + 4.0 * g_yyz_y_xy_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_xy_xy[i] = -2.0 * g_z_y_xy_xy[i] * a_exp + 4.0 * g_yyz_y_xy_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_xy_xz[i] = -2.0 * g_z_y_xy_xz[i] * a_exp + 4.0 * g_yyz_y_xy_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_xy_yy[i] = -2.0 * g_z_y_xy_yy[i] * a_exp + 4.0 * g_yyz_y_xy_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_xy_yz[i] = -2.0 * g_z_y_xy_yz[i] * a_exp + 4.0 * g_yyz_y_xy_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_xy_zz[i] = -2.0 * g_z_y_xy_zz[i] * a_exp + 4.0 * g_yyz_y_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1236-1242)

    #pragma omp simd aligned(g_yy_0_0_0_z_y_xz_xx, g_yy_0_0_0_z_y_xz_xy, g_yy_0_0_0_z_y_xz_xz, g_yy_0_0_0_z_y_xz_yy, g_yy_0_0_0_z_y_xz_yz, g_yy_0_0_0_z_y_xz_zz, g_yyz_y_xz_xx, g_yyz_y_xz_xy, g_yyz_y_xz_xz, g_yyz_y_xz_yy, g_yyz_y_xz_yz, g_yyz_y_xz_zz, g_z_y_xz_xx, g_z_y_xz_xy, g_z_y_xz_xz, g_z_y_xz_yy, g_z_y_xz_yz, g_z_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_y_xz_xx[i] = -2.0 * g_z_y_xz_xx[i] * a_exp + 4.0 * g_yyz_y_xz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_xz_xy[i] = -2.0 * g_z_y_xz_xy[i] * a_exp + 4.0 * g_yyz_y_xz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_xz_xz[i] = -2.0 * g_z_y_xz_xz[i] * a_exp + 4.0 * g_yyz_y_xz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_xz_yy[i] = -2.0 * g_z_y_xz_yy[i] * a_exp + 4.0 * g_yyz_y_xz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_xz_yz[i] = -2.0 * g_z_y_xz_yz[i] * a_exp + 4.0 * g_yyz_y_xz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_xz_zz[i] = -2.0 * g_z_y_xz_zz[i] * a_exp + 4.0 * g_yyz_y_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1242-1248)

    #pragma omp simd aligned(g_yy_0_0_0_z_y_yy_xx, g_yy_0_0_0_z_y_yy_xy, g_yy_0_0_0_z_y_yy_xz, g_yy_0_0_0_z_y_yy_yy, g_yy_0_0_0_z_y_yy_yz, g_yy_0_0_0_z_y_yy_zz, g_yyz_y_yy_xx, g_yyz_y_yy_xy, g_yyz_y_yy_xz, g_yyz_y_yy_yy, g_yyz_y_yy_yz, g_yyz_y_yy_zz, g_z_y_yy_xx, g_z_y_yy_xy, g_z_y_yy_xz, g_z_y_yy_yy, g_z_y_yy_yz, g_z_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_y_yy_xx[i] = -2.0 * g_z_y_yy_xx[i] * a_exp + 4.0 * g_yyz_y_yy_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_yy_xy[i] = -2.0 * g_z_y_yy_xy[i] * a_exp + 4.0 * g_yyz_y_yy_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_yy_xz[i] = -2.0 * g_z_y_yy_xz[i] * a_exp + 4.0 * g_yyz_y_yy_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_yy_yy[i] = -2.0 * g_z_y_yy_yy[i] * a_exp + 4.0 * g_yyz_y_yy_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_yy_yz[i] = -2.0 * g_z_y_yy_yz[i] * a_exp + 4.0 * g_yyz_y_yy_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_yy_zz[i] = -2.0 * g_z_y_yy_zz[i] * a_exp + 4.0 * g_yyz_y_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1248-1254)

    #pragma omp simd aligned(g_yy_0_0_0_z_y_yz_xx, g_yy_0_0_0_z_y_yz_xy, g_yy_0_0_0_z_y_yz_xz, g_yy_0_0_0_z_y_yz_yy, g_yy_0_0_0_z_y_yz_yz, g_yy_0_0_0_z_y_yz_zz, g_yyz_y_yz_xx, g_yyz_y_yz_xy, g_yyz_y_yz_xz, g_yyz_y_yz_yy, g_yyz_y_yz_yz, g_yyz_y_yz_zz, g_z_y_yz_xx, g_z_y_yz_xy, g_z_y_yz_xz, g_z_y_yz_yy, g_z_y_yz_yz, g_z_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_y_yz_xx[i] = -2.0 * g_z_y_yz_xx[i] * a_exp + 4.0 * g_yyz_y_yz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_yz_xy[i] = -2.0 * g_z_y_yz_xy[i] * a_exp + 4.0 * g_yyz_y_yz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_yz_xz[i] = -2.0 * g_z_y_yz_xz[i] * a_exp + 4.0 * g_yyz_y_yz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_yz_yy[i] = -2.0 * g_z_y_yz_yy[i] * a_exp + 4.0 * g_yyz_y_yz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_yz_yz[i] = -2.0 * g_z_y_yz_yz[i] * a_exp + 4.0 * g_yyz_y_yz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_yz_zz[i] = -2.0 * g_z_y_yz_zz[i] * a_exp + 4.0 * g_yyz_y_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1254-1260)

    #pragma omp simd aligned(g_yy_0_0_0_z_y_zz_xx, g_yy_0_0_0_z_y_zz_xy, g_yy_0_0_0_z_y_zz_xz, g_yy_0_0_0_z_y_zz_yy, g_yy_0_0_0_z_y_zz_yz, g_yy_0_0_0_z_y_zz_zz, g_yyz_y_zz_xx, g_yyz_y_zz_xy, g_yyz_y_zz_xz, g_yyz_y_zz_yy, g_yyz_y_zz_yz, g_yyz_y_zz_zz, g_z_y_zz_xx, g_z_y_zz_xy, g_z_y_zz_xz, g_z_y_zz_yy, g_z_y_zz_yz, g_z_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_y_zz_xx[i] = -2.0 * g_z_y_zz_xx[i] * a_exp + 4.0 * g_yyz_y_zz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_zz_xy[i] = -2.0 * g_z_y_zz_xy[i] * a_exp + 4.0 * g_yyz_y_zz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_zz_xz[i] = -2.0 * g_z_y_zz_xz[i] * a_exp + 4.0 * g_yyz_y_zz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_zz_yy[i] = -2.0 * g_z_y_zz_yy[i] * a_exp + 4.0 * g_yyz_y_zz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_zz_yz[i] = -2.0 * g_z_y_zz_yz[i] * a_exp + 4.0 * g_yyz_y_zz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_zz_zz[i] = -2.0 * g_z_y_zz_zz[i] * a_exp + 4.0 * g_yyz_y_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1260-1266)

    #pragma omp simd aligned(g_yy_0_0_0_z_z_xx_xx, g_yy_0_0_0_z_z_xx_xy, g_yy_0_0_0_z_z_xx_xz, g_yy_0_0_0_z_z_xx_yy, g_yy_0_0_0_z_z_xx_yz, g_yy_0_0_0_z_z_xx_zz, g_yyz_z_xx_xx, g_yyz_z_xx_xy, g_yyz_z_xx_xz, g_yyz_z_xx_yy, g_yyz_z_xx_yz, g_yyz_z_xx_zz, g_z_z_xx_xx, g_z_z_xx_xy, g_z_z_xx_xz, g_z_z_xx_yy, g_z_z_xx_yz, g_z_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_z_xx_xx[i] = -2.0 * g_z_z_xx_xx[i] * a_exp + 4.0 * g_yyz_z_xx_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_xx_xy[i] = -2.0 * g_z_z_xx_xy[i] * a_exp + 4.0 * g_yyz_z_xx_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_xx_xz[i] = -2.0 * g_z_z_xx_xz[i] * a_exp + 4.0 * g_yyz_z_xx_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_xx_yy[i] = -2.0 * g_z_z_xx_yy[i] * a_exp + 4.0 * g_yyz_z_xx_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_xx_yz[i] = -2.0 * g_z_z_xx_yz[i] * a_exp + 4.0 * g_yyz_z_xx_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_xx_zz[i] = -2.0 * g_z_z_xx_zz[i] * a_exp + 4.0 * g_yyz_z_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (1266-1272)

    #pragma omp simd aligned(g_yy_0_0_0_z_z_xy_xx, g_yy_0_0_0_z_z_xy_xy, g_yy_0_0_0_z_z_xy_xz, g_yy_0_0_0_z_z_xy_yy, g_yy_0_0_0_z_z_xy_yz, g_yy_0_0_0_z_z_xy_zz, g_yyz_z_xy_xx, g_yyz_z_xy_xy, g_yyz_z_xy_xz, g_yyz_z_xy_yy, g_yyz_z_xy_yz, g_yyz_z_xy_zz, g_z_z_xy_xx, g_z_z_xy_xy, g_z_z_xy_xz, g_z_z_xy_yy, g_z_z_xy_yz, g_z_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_z_xy_xx[i] = -2.0 * g_z_z_xy_xx[i] * a_exp + 4.0 * g_yyz_z_xy_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_xy_xy[i] = -2.0 * g_z_z_xy_xy[i] * a_exp + 4.0 * g_yyz_z_xy_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_xy_xz[i] = -2.0 * g_z_z_xy_xz[i] * a_exp + 4.0 * g_yyz_z_xy_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_xy_yy[i] = -2.0 * g_z_z_xy_yy[i] * a_exp + 4.0 * g_yyz_z_xy_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_xy_yz[i] = -2.0 * g_z_z_xy_yz[i] * a_exp + 4.0 * g_yyz_z_xy_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_xy_zz[i] = -2.0 * g_z_z_xy_zz[i] * a_exp + 4.0 * g_yyz_z_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1272-1278)

    #pragma omp simd aligned(g_yy_0_0_0_z_z_xz_xx, g_yy_0_0_0_z_z_xz_xy, g_yy_0_0_0_z_z_xz_xz, g_yy_0_0_0_z_z_xz_yy, g_yy_0_0_0_z_z_xz_yz, g_yy_0_0_0_z_z_xz_zz, g_yyz_z_xz_xx, g_yyz_z_xz_xy, g_yyz_z_xz_xz, g_yyz_z_xz_yy, g_yyz_z_xz_yz, g_yyz_z_xz_zz, g_z_z_xz_xx, g_z_z_xz_xy, g_z_z_xz_xz, g_z_z_xz_yy, g_z_z_xz_yz, g_z_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_z_xz_xx[i] = -2.0 * g_z_z_xz_xx[i] * a_exp + 4.0 * g_yyz_z_xz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_xz_xy[i] = -2.0 * g_z_z_xz_xy[i] * a_exp + 4.0 * g_yyz_z_xz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_xz_xz[i] = -2.0 * g_z_z_xz_xz[i] * a_exp + 4.0 * g_yyz_z_xz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_xz_yy[i] = -2.0 * g_z_z_xz_yy[i] * a_exp + 4.0 * g_yyz_z_xz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_xz_yz[i] = -2.0 * g_z_z_xz_yz[i] * a_exp + 4.0 * g_yyz_z_xz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_xz_zz[i] = -2.0 * g_z_z_xz_zz[i] * a_exp + 4.0 * g_yyz_z_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1278-1284)

    #pragma omp simd aligned(g_yy_0_0_0_z_z_yy_xx, g_yy_0_0_0_z_z_yy_xy, g_yy_0_0_0_z_z_yy_xz, g_yy_0_0_0_z_z_yy_yy, g_yy_0_0_0_z_z_yy_yz, g_yy_0_0_0_z_z_yy_zz, g_yyz_z_yy_xx, g_yyz_z_yy_xy, g_yyz_z_yy_xz, g_yyz_z_yy_yy, g_yyz_z_yy_yz, g_yyz_z_yy_zz, g_z_z_yy_xx, g_z_z_yy_xy, g_z_z_yy_xz, g_z_z_yy_yy, g_z_z_yy_yz, g_z_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_z_yy_xx[i] = -2.0 * g_z_z_yy_xx[i] * a_exp + 4.0 * g_yyz_z_yy_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_yy_xy[i] = -2.0 * g_z_z_yy_xy[i] * a_exp + 4.0 * g_yyz_z_yy_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_yy_xz[i] = -2.0 * g_z_z_yy_xz[i] * a_exp + 4.0 * g_yyz_z_yy_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_yy_yy[i] = -2.0 * g_z_z_yy_yy[i] * a_exp + 4.0 * g_yyz_z_yy_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_yy_yz[i] = -2.0 * g_z_z_yy_yz[i] * a_exp + 4.0 * g_yyz_z_yy_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_yy_zz[i] = -2.0 * g_z_z_yy_zz[i] * a_exp + 4.0 * g_yyz_z_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1284-1290)

    #pragma omp simd aligned(g_yy_0_0_0_z_z_yz_xx, g_yy_0_0_0_z_z_yz_xy, g_yy_0_0_0_z_z_yz_xz, g_yy_0_0_0_z_z_yz_yy, g_yy_0_0_0_z_z_yz_yz, g_yy_0_0_0_z_z_yz_zz, g_yyz_z_yz_xx, g_yyz_z_yz_xy, g_yyz_z_yz_xz, g_yyz_z_yz_yy, g_yyz_z_yz_yz, g_yyz_z_yz_zz, g_z_z_yz_xx, g_z_z_yz_xy, g_z_z_yz_xz, g_z_z_yz_yy, g_z_z_yz_yz, g_z_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_z_yz_xx[i] = -2.0 * g_z_z_yz_xx[i] * a_exp + 4.0 * g_yyz_z_yz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_yz_xy[i] = -2.0 * g_z_z_yz_xy[i] * a_exp + 4.0 * g_yyz_z_yz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_yz_xz[i] = -2.0 * g_z_z_yz_xz[i] * a_exp + 4.0 * g_yyz_z_yz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_yz_yy[i] = -2.0 * g_z_z_yz_yy[i] * a_exp + 4.0 * g_yyz_z_yz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_yz_yz[i] = -2.0 * g_z_z_yz_yz[i] * a_exp + 4.0 * g_yyz_z_yz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_yz_zz[i] = -2.0 * g_z_z_yz_zz[i] * a_exp + 4.0 * g_yyz_z_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1290-1296)

    #pragma omp simd aligned(g_yy_0_0_0_z_z_zz_xx, g_yy_0_0_0_z_z_zz_xy, g_yy_0_0_0_z_z_zz_xz, g_yy_0_0_0_z_z_zz_yy, g_yy_0_0_0_z_z_zz_yz, g_yy_0_0_0_z_z_zz_zz, g_yyz_z_zz_xx, g_yyz_z_zz_xy, g_yyz_z_zz_xz, g_yyz_z_zz_yy, g_yyz_z_zz_yz, g_yyz_z_zz_zz, g_z_z_zz_xx, g_z_z_zz_xy, g_z_z_zz_xz, g_z_z_zz_yy, g_z_z_zz_yz, g_z_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_z_zz_xx[i] = -2.0 * g_z_z_zz_xx[i] * a_exp + 4.0 * g_yyz_z_zz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_zz_xy[i] = -2.0 * g_z_z_zz_xy[i] * a_exp + 4.0 * g_yyz_z_zz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_zz_xz[i] = -2.0 * g_z_z_zz_xz[i] * a_exp + 4.0 * g_yyz_z_zz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_zz_yy[i] = -2.0 * g_z_z_zz_yy[i] * a_exp + 4.0 * g_yyz_z_zz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_zz_yz[i] = -2.0 * g_z_z_zz_yz[i] * a_exp + 4.0 * g_yyz_z_zz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_zz_zz[i] = -2.0 * g_z_z_zz_zz[i] * a_exp + 4.0 * g_yyz_z_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1296-1302)

    #pragma omp simd aligned(g_xyz_x_xx_xx, g_xyz_x_xx_xy, g_xyz_x_xx_xz, g_xyz_x_xx_yy, g_xyz_x_xx_yz, g_xyz_x_xx_zz, g_yz_0_0_0_x_x_xx_xx, g_yz_0_0_0_x_x_xx_xy, g_yz_0_0_0_x_x_xx_xz, g_yz_0_0_0_x_x_xx_yy, g_yz_0_0_0_x_x_xx_yz, g_yz_0_0_0_x_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_x_xx_xx[i] = 4.0 * g_xyz_x_xx_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_xx_xy[i] = 4.0 * g_xyz_x_xx_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_xx_xz[i] = 4.0 * g_xyz_x_xx_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_xx_yy[i] = 4.0 * g_xyz_x_xx_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_xx_yz[i] = 4.0 * g_xyz_x_xx_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_xx_zz[i] = 4.0 * g_xyz_x_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (1302-1308)

    #pragma omp simd aligned(g_xyz_x_xy_xx, g_xyz_x_xy_xy, g_xyz_x_xy_xz, g_xyz_x_xy_yy, g_xyz_x_xy_yz, g_xyz_x_xy_zz, g_yz_0_0_0_x_x_xy_xx, g_yz_0_0_0_x_x_xy_xy, g_yz_0_0_0_x_x_xy_xz, g_yz_0_0_0_x_x_xy_yy, g_yz_0_0_0_x_x_xy_yz, g_yz_0_0_0_x_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_x_xy_xx[i] = 4.0 * g_xyz_x_xy_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_xy_xy[i] = 4.0 * g_xyz_x_xy_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_xy_xz[i] = 4.0 * g_xyz_x_xy_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_xy_yy[i] = 4.0 * g_xyz_x_xy_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_xy_yz[i] = 4.0 * g_xyz_x_xy_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_xy_zz[i] = 4.0 * g_xyz_x_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1308-1314)

    #pragma omp simd aligned(g_xyz_x_xz_xx, g_xyz_x_xz_xy, g_xyz_x_xz_xz, g_xyz_x_xz_yy, g_xyz_x_xz_yz, g_xyz_x_xz_zz, g_yz_0_0_0_x_x_xz_xx, g_yz_0_0_0_x_x_xz_xy, g_yz_0_0_0_x_x_xz_xz, g_yz_0_0_0_x_x_xz_yy, g_yz_0_0_0_x_x_xz_yz, g_yz_0_0_0_x_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_x_xz_xx[i] = 4.0 * g_xyz_x_xz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_xz_xy[i] = 4.0 * g_xyz_x_xz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_xz_xz[i] = 4.0 * g_xyz_x_xz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_xz_yy[i] = 4.0 * g_xyz_x_xz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_xz_yz[i] = 4.0 * g_xyz_x_xz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_xz_zz[i] = 4.0 * g_xyz_x_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1314-1320)

    #pragma omp simd aligned(g_xyz_x_yy_xx, g_xyz_x_yy_xy, g_xyz_x_yy_xz, g_xyz_x_yy_yy, g_xyz_x_yy_yz, g_xyz_x_yy_zz, g_yz_0_0_0_x_x_yy_xx, g_yz_0_0_0_x_x_yy_xy, g_yz_0_0_0_x_x_yy_xz, g_yz_0_0_0_x_x_yy_yy, g_yz_0_0_0_x_x_yy_yz, g_yz_0_0_0_x_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_x_yy_xx[i] = 4.0 * g_xyz_x_yy_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_yy_xy[i] = 4.0 * g_xyz_x_yy_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_yy_xz[i] = 4.0 * g_xyz_x_yy_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_yy_yy[i] = 4.0 * g_xyz_x_yy_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_yy_yz[i] = 4.0 * g_xyz_x_yy_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_yy_zz[i] = 4.0 * g_xyz_x_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1320-1326)

    #pragma omp simd aligned(g_xyz_x_yz_xx, g_xyz_x_yz_xy, g_xyz_x_yz_xz, g_xyz_x_yz_yy, g_xyz_x_yz_yz, g_xyz_x_yz_zz, g_yz_0_0_0_x_x_yz_xx, g_yz_0_0_0_x_x_yz_xy, g_yz_0_0_0_x_x_yz_xz, g_yz_0_0_0_x_x_yz_yy, g_yz_0_0_0_x_x_yz_yz, g_yz_0_0_0_x_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_x_yz_xx[i] = 4.0 * g_xyz_x_yz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_yz_xy[i] = 4.0 * g_xyz_x_yz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_yz_xz[i] = 4.0 * g_xyz_x_yz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_yz_yy[i] = 4.0 * g_xyz_x_yz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_yz_yz[i] = 4.0 * g_xyz_x_yz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_yz_zz[i] = 4.0 * g_xyz_x_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1326-1332)

    #pragma omp simd aligned(g_xyz_x_zz_xx, g_xyz_x_zz_xy, g_xyz_x_zz_xz, g_xyz_x_zz_yy, g_xyz_x_zz_yz, g_xyz_x_zz_zz, g_yz_0_0_0_x_x_zz_xx, g_yz_0_0_0_x_x_zz_xy, g_yz_0_0_0_x_x_zz_xz, g_yz_0_0_0_x_x_zz_yy, g_yz_0_0_0_x_x_zz_yz, g_yz_0_0_0_x_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_x_zz_xx[i] = 4.0 * g_xyz_x_zz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_zz_xy[i] = 4.0 * g_xyz_x_zz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_zz_xz[i] = 4.0 * g_xyz_x_zz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_zz_yy[i] = 4.0 * g_xyz_x_zz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_zz_yz[i] = 4.0 * g_xyz_x_zz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_zz_zz[i] = 4.0 * g_xyz_x_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1332-1338)

    #pragma omp simd aligned(g_xyz_y_xx_xx, g_xyz_y_xx_xy, g_xyz_y_xx_xz, g_xyz_y_xx_yy, g_xyz_y_xx_yz, g_xyz_y_xx_zz, g_yz_0_0_0_x_y_xx_xx, g_yz_0_0_0_x_y_xx_xy, g_yz_0_0_0_x_y_xx_xz, g_yz_0_0_0_x_y_xx_yy, g_yz_0_0_0_x_y_xx_yz, g_yz_0_0_0_x_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_y_xx_xx[i] = 4.0 * g_xyz_y_xx_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_xx_xy[i] = 4.0 * g_xyz_y_xx_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_xx_xz[i] = 4.0 * g_xyz_y_xx_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_xx_yy[i] = 4.0 * g_xyz_y_xx_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_xx_yz[i] = 4.0 * g_xyz_y_xx_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_xx_zz[i] = 4.0 * g_xyz_y_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (1338-1344)

    #pragma omp simd aligned(g_xyz_y_xy_xx, g_xyz_y_xy_xy, g_xyz_y_xy_xz, g_xyz_y_xy_yy, g_xyz_y_xy_yz, g_xyz_y_xy_zz, g_yz_0_0_0_x_y_xy_xx, g_yz_0_0_0_x_y_xy_xy, g_yz_0_0_0_x_y_xy_xz, g_yz_0_0_0_x_y_xy_yy, g_yz_0_0_0_x_y_xy_yz, g_yz_0_0_0_x_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_y_xy_xx[i] = 4.0 * g_xyz_y_xy_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_xy_xy[i] = 4.0 * g_xyz_y_xy_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_xy_xz[i] = 4.0 * g_xyz_y_xy_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_xy_yy[i] = 4.0 * g_xyz_y_xy_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_xy_yz[i] = 4.0 * g_xyz_y_xy_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_xy_zz[i] = 4.0 * g_xyz_y_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1344-1350)

    #pragma omp simd aligned(g_xyz_y_xz_xx, g_xyz_y_xz_xy, g_xyz_y_xz_xz, g_xyz_y_xz_yy, g_xyz_y_xz_yz, g_xyz_y_xz_zz, g_yz_0_0_0_x_y_xz_xx, g_yz_0_0_0_x_y_xz_xy, g_yz_0_0_0_x_y_xz_xz, g_yz_0_0_0_x_y_xz_yy, g_yz_0_0_0_x_y_xz_yz, g_yz_0_0_0_x_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_y_xz_xx[i] = 4.0 * g_xyz_y_xz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_xz_xy[i] = 4.0 * g_xyz_y_xz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_xz_xz[i] = 4.0 * g_xyz_y_xz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_xz_yy[i] = 4.0 * g_xyz_y_xz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_xz_yz[i] = 4.0 * g_xyz_y_xz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_xz_zz[i] = 4.0 * g_xyz_y_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1350-1356)

    #pragma omp simd aligned(g_xyz_y_yy_xx, g_xyz_y_yy_xy, g_xyz_y_yy_xz, g_xyz_y_yy_yy, g_xyz_y_yy_yz, g_xyz_y_yy_zz, g_yz_0_0_0_x_y_yy_xx, g_yz_0_0_0_x_y_yy_xy, g_yz_0_0_0_x_y_yy_xz, g_yz_0_0_0_x_y_yy_yy, g_yz_0_0_0_x_y_yy_yz, g_yz_0_0_0_x_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_y_yy_xx[i] = 4.0 * g_xyz_y_yy_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_yy_xy[i] = 4.0 * g_xyz_y_yy_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_yy_xz[i] = 4.0 * g_xyz_y_yy_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_yy_yy[i] = 4.0 * g_xyz_y_yy_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_yy_yz[i] = 4.0 * g_xyz_y_yy_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_yy_zz[i] = 4.0 * g_xyz_y_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1356-1362)

    #pragma omp simd aligned(g_xyz_y_yz_xx, g_xyz_y_yz_xy, g_xyz_y_yz_xz, g_xyz_y_yz_yy, g_xyz_y_yz_yz, g_xyz_y_yz_zz, g_yz_0_0_0_x_y_yz_xx, g_yz_0_0_0_x_y_yz_xy, g_yz_0_0_0_x_y_yz_xz, g_yz_0_0_0_x_y_yz_yy, g_yz_0_0_0_x_y_yz_yz, g_yz_0_0_0_x_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_y_yz_xx[i] = 4.0 * g_xyz_y_yz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_yz_xy[i] = 4.0 * g_xyz_y_yz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_yz_xz[i] = 4.0 * g_xyz_y_yz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_yz_yy[i] = 4.0 * g_xyz_y_yz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_yz_yz[i] = 4.0 * g_xyz_y_yz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_yz_zz[i] = 4.0 * g_xyz_y_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1362-1368)

    #pragma omp simd aligned(g_xyz_y_zz_xx, g_xyz_y_zz_xy, g_xyz_y_zz_xz, g_xyz_y_zz_yy, g_xyz_y_zz_yz, g_xyz_y_zz_zz, g_yz_0_0_0_x_y_zz_xx, g_yz_0_0_0_x_y_zz_xy, g_yz_0_0_0_x_y_zz_xz, g_yz_0_0_0_x_y_zz_yy, g_yz_0_0_0_x_y_zz_yz, g_yz_0_0_0_x_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_y_zz_xx[i] = 4.0 * g_xyz_y_zz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_zz_xy[i] = 4.0 * g_xyz_y_zz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_zz_xz[i] = 4.0 * g_xyz_y_zz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_zz_yy[i] = 4.0 * g_xyz_y_zz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_zz_yz[i] = 4.0 * g_xyz_y_zz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_zz_zz[i] = 4.0 * g_xyz_y_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1368-1374)

    #pragma omp simd aligned(g_xyz_z_xx_xx, g_xyz_z_xx_xy, g_xyz_z_xx_xz, g_xyz_z_xx_yy, g_xyz_z_xx_yz, g_xyz_z_xx_zz, g_yz_0_0_0_x_z_xx_xx, g_yz_0_0_0_x_z_xx_xy, g_yz_0_0_0_x_z_xx_xz, g_yz_0_0_0_x_z_xx_yy, g_yz_0_0_0_x_z_xx_yz, g_yz_0_0_0_x_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_z_xx_xx[i] = 4.0 * g_xyz_z_xx_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_xx_xy[i] = 4.0 * g_xyz_z_xx_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_xx_xz[i] = 4.0 * g_xyz_z_xx_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_xx_yy[i] = 4.0 * g_xyz_z_xx_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_xx_yz[i] = 4.0 * g_xyz_z_xx_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_xx_zz[i] = 4.0 * g_xyz_z_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (1374-1380)

    #pragma omp simd aligned(g_xyz_z_xy_xx, g_xyz_z_xy_xy, g_xyz_z_xy_xz, g_xyz_z_xy_yy, g_xyz_z_xy_yz, g_xyz_z_xy_zz, g_yz_0_0_0_x_z_xy_xx, g_yz_0_0_0_x_z_xy_xy, g_yz_0_0_0_x_z_xy_xz, g_yz_0_0_0_x_z_xy_yy, g_yz_0_0_0_x_z_xy_yz, g_yz_0_0_0_x_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_z_xy_xx[i] = 4.0 * g_xyz_z_xy_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_xy_xy[i] = 4.0 * g_xyz_z_xy_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_xy_xz[i] = 4.0 * g_xyz_z_xy_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_xy_yy[i] = 4.0 * g_xyz_z_xy_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_xy_yz[i] = 4.0 * g_xyz_z_xy_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_xy_zz[i] = 4.0 * g_xyz_z_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1380-1386)

    #pragma omp simd aligned(g_xyz_z_xz_xx, g_xyz_z_xz_xy, g_xyz_z_xz_xz, g_xyz_z_xz_yy, g_xyz_z_xz_yz, g_xyz_z_xz_zz, g_yz_0_0_0_x_z_xz_xx, g_yz_0_0_0_x_z_xz_xy, g_yz_0_0_0_x_z_xz_xz, g_yz_0_0_0_x_z_xz_yy, g_yz_0_0_0_x_z_xz_yz, g_yz_0_0_0_x_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_z_xz_xx[i] = 4.0 * g_xyz_z_xz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_xz_xy[i] = 4.0 * g_xyz_z_xz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_xz_xz[i] = 4.0 * g_xyz_z_xz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_xz_yy[i] = 4.0 * g_xyz_z_xz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_xz_yz[i] = 4.0 * g_xyz_z_xz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_xz_zz[i] = 4.0 * g_xyz_z_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1386-1392)

    #pragma omp simd aligned(g_xyz_z_yy_xx, g_xyz_z_yy_xy, g_xyz_z_yy_xz, g_xyz_z_yy_yy, g_xyz_z_yy_yz, g_xyz_z_yy_zz, g_yz_0_0_0_x_z_yy_xx, g_yz_0_0_0_x_z_yy_xy, g_yz_0_0_0_x_z_yy_xz, g_yz_0_0_0_x_z_yy_yy, g_yz_0_0_0_x_z_yy_yz, g_yz_0_0_0_x_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_z_yy_xx[i] = 4.0 * g_xyz_z_yy_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_yy_xy[i] = 4.0 * g_xyz_z_yy_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_yy_xz[i] = 4.0 * g_xyz_z_yy_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_yy_yy[i] = 4.0 * g_xyz_z_yy_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_yy_yz[i] = 4.0 * g_xyz_z_yy_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_yy_zz[i] = 4.0 * g_xyz_z_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1392-1398)

    #pragma omp simd aligned(g_xyz_z_yz_xx, g_xyz_z_yz_xy, g_xyz_z_yz_xz, g_xyz_z_yz_yy, g_xyz_z_yz_yz, g_xyz_z_yz_zz, g_yz_0_0_0_x_z_yz_xx, g_yz_0_0_0_x_z_yz_xy, g_yz_0_0_0_x_z_yz_xz, g_yz_0_0_0_x_z_yz_yy, g_yz_0_0_0_x_z_yz_yz, g_yz_0_0_0_x_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_z_yz_xx[i] = 4.0 * g_xyz_z_yz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_yz_xy[i] = 4.0 * g_xyz_z_yz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_yz_xz[i] = 4.0 * g_xyz_z_yz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_yz_yy[i] = 4.0 * g_xyz_z_yz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_yz_yz[i] = 4.0 * g_xyz_z_yz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_yz_zz[i] = 4.0 * g_xyz_z_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1398-1404)

    #pragma omp simd aligned(g_xyz_z_zz_xx, g_xyz_z_zz_xy, g_xyz_z_zz_xz, g_xyz_z_zz_yy, g_xyz_z_zz_yz, g_xyz_z_zz_zz, g_yz_0_0_0_x_z_zz_xx, g_yz_0_0_0_x_z_zz_xy, g_yz_0_0_0_x_z_zz_xz, g_yz_0_0_0_x_z_zz_yy, g_yz_0_0_0_x_z_zz_yz, g_yz_0_0_0_x_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_z_zz_xx[i] = 4.0 * g_xyz_z_zz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_zz_xy[i] = 4.0 * g_xyz_z_zz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_zz_xz[i] = 4.0 * g_xyz_z_zz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_zz_yy[i] = 4.0 * g_xyz_z_zz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_zz_yz[i] = 4.0 * g_xyz_z_zz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_zz_zz[i] = 4.0 * g_xyz_z_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1404-1410)

    #pragma omp simd aligned(g_yyz_x_xx_xx, g_yyz_x_xx_xy, g_yyz_x_xx_xz, g_yyz_x_xx_yy, g_yyz_x_xx_yz, g_yyz_x_xx_zz, g_yz_0_0_0_y_x_xx_xx, g_yz_0_0_0_y_x_xx_xy, g_yz_0_0_0_y_x_xx_xz, g_yz_0_0_0_y_x_xx_yy, g_yz_0_0_0_y_x_xx_yz, g_yz_0_0_0_y_x_xx_zz, g_z_x_xx_xx, g_z_x_xx_xy, g_z_x_xx_xz, g_z_x_xx_yy, g_z_x_xx_yz, g_z_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_x_xx_xx[i] = -2.0 * g_z_x_xx_xx[i] * a_exp + 4.0 * g_yyz_x_xx_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_xx_xy[i] = -2.0 * g_z_x_xx_xy[i] * a_exp + 4.0 * g_yyz_x_xx_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_xx_xz[i] = -2.0 * g_z_x_xx_xz[i] * a_exp + 4.0 * g_yyz_x_xx_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_xx_yy[i] = -2.0 * g_z_x_xx_yy[i] * a_exp + 4.0 * g_yyz_x_xx_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_xx_yz[i] = -2.0 * g_z_x_xx_yz[i] * a_exp + 4.0 * g_yyz_x_xx_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_xx_zz[i] = -2.0 * g_z_x_xx_zz[i] * a_exp + 4.0 * g_yyz_x_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (1410-1416)

    #pragma omp simd aligned(g_yyz_x_xy_xx, g_yyz_x_xy_xy, g_yyz_x_xy_xz, g_yyz_x_xy_yy, g_yyz_x_xy_yz, g_yyz_x_xy_zz, g_yz_0_0_0_y_x_xy_xx, g_yz_0_0_0_y_x_xy_xy, g_yz_0_0_0_y_x_xy_xz, g_yz_0_0_0_y_x_xy_yy, g_yz_0_0_0_y_x_xy_yz, g_yz_0_0_0_y_x_xy_zz, g_z_x_xy_xx, g_z_x_xy_xy, g_z_x_xy_xz, g_z_x_xy_yy, g_z_x_xy_yz, g_z_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_x_xy_xx[i] = -2.0 * g_z_x_xy_xx[i] * a_exp + 4.0 * g_yyz_x_xy_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_xy_xy[i] = -2.0 * g_z_x_xy_xy[i] * a_exp + 4.0 * g_yyz_x_xy_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_xy_xz[i] = -2.0 * g_z_x_xy_xz[i] * a_exp + 4.0 * g_yyz_x_xy_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_xy_yy[i] = -2.0 * g_z_x_xy_yy[i] * a_exp + 4.0 * g_yyz_x_xy_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_xy_yz[i] = -2.0 * g_z_x_xy_yz[i] * a_exp + 4.0 * g_yyz_x_xy_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_xy_zz[i] = -2.0 * g_z_x_xy_zz[i] * a_exp + 4.0 * g_yyz_x_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1416-1422)

    #pragma omp simd aligned(g_yyz_x_xz_xx, g_yyz_x_xz_xy, g_yyz_x_xz_xz, g_yyz_x_xz_yy, g_yyz_x_xz_yz, g_yyz_x_xz_zz, g_yz_0_0_0_y_x_xz_xx, g_yz_0_0_0_y_x_xz_xy, g_yz_0_0_0_y_x_xz_xz, g_yz_0_0_0_y_x_xz_yy, g_yz_0_0_0_y_x_xz_yz, g_yz_0_0_0_y_x_xz_zz, g_z_x_xz_xx, g_z_x_xz_xy, g_z_x_xz_xz, g_z_x_xz_yy, g_z_x_xz_yz, g_z_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_x_xz_xx[i] = -2.0 * g_z_x_xz_xx[i] * a_exp + 4.0 * g_yyz_x_xz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_xz_xy[i] = -2.0 * g_z_x_xz_xy[i] * a_exp + 4.0 * g_yyz_x_xz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_xz_xz[i] = -2.0 * g_z_x_xz_xz[i] * a_exp + 4.0 * g_yyz_x_xz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_xz_yy[i] = -2.0 * g_z_x_xz_yy[i] * a_exp + 4.0 * g_yyz_x_xz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_xz_yz[i] = -2.0 * g_z_x_xz_yz[i] * a_exp + 4.0 * g_yyz_x_xz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_xz_zz[i] = -2.0 * g_z_x_xz_zz[i] * a_exp + 4.0 * g_yyz_x_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1422-1428)

    #pragma omp simd aligned(g_yyz_x_yy_xx, g_yyz_x_yy_xy, g_yyz_x_yy_xz, g_yyz_x_yy_yy, g_yyz_x_yy_yz, g_yyz_x_yy_zz, g_yz_0_0_0_y_x_yy_xx, g_yz_0_0_0_y_x_yy_xy, g_yz_0_0_0_y_x_yy_xz, g_yz_0_0_0_y_x_yy_yy, g_yz_0_0_0_y_x_yy_yz, g_yz_0_0_0_y_x_yy_zz, g_z_x_yy_xx, g_z_x_yy_xy, g_z_x_yy_xz, g_z_x_yy_yy, g_z_x_yy_yz, g_z_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_x_yy_xx[i] = -2.0 * g_z_x_yy_xx[i] * a_exp + 4.0 * g_yyz_x_yy_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_yy_xy[i] = -2.0 * g_z_x_yy_xy[i] * a_exp + 4.0 * g_yyz_x_yy_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_yy_xz[i] = -2.0 * g_z_x_yy_xz[i] * a_exp + 4.0 * g_yyz_x_yy_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_yy_yy[i] = -2.0 * g_z_x_yy_yy[i] * a_exp + 4.0 * g_yyz_x_yy_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_yy_yz[i] = -2.0 * g_z_x_yy_yz[i] * a_exp + 4.0 * g_yyz_x_yy_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_yy_zz[i] = -2.0 * g_z_x_yy_zz[i] * a_exp + 4.0 * g_yyz_x_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1428-1434)

    #pragma omp simd aligned(g_yyz_x_yz_xx, g_yyz_x_yz_xy, g_yyz_x_yz_xz, g_yyz_x_yz_yy, g_yyz_x_yz_yz, g_yyz_x_yz_zz, g_yz_0_0_0_y_x_yz_xx, g_yz_0_0_0_y_x_yz_xy, g_yz_0_0_0_y_x_yz_xz, g_yz_0_0_0_y_x_yz_yy, g_yz_0_0_0_y_x_yz_yz, g_yz_0_0_0_y_x_yz_zz, g_z_x_yz_xx, g_z_x_yz_xy, g_z_x_yz_xz, g_z_x_yz_yy, g_z_x_yz_yz, g_z_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_x_yz_xx[i] = -2.0 * g_z_x_yz_xx[i] * a_exp + 4.0 * g_yyz_x_yz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_yz_xy[i] = -2.0 * g_z_x_yz_xy[i] * a_exp + 4.0 * g_yyz_x_yz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_yz_xz[i] = -2.0 * g_z_x_yz_xz[i] * a_exp + 4.0 * g_yyz_x_yz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_yz_yy[i] = -2.0 * g_z_x_yz_yy[i] * a_exp + 4.0 * g_yyz_x_yz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_yz_yz[i] = -2.0 * g_z_x_yz_yz[i] * a_exp + 4.0 * g_yyz_x_yz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_yz_zz[i] = -2.0 * g_z_x_yz_zz[i] * a_exp + 4.0 * g_yyz_x_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1434-1440)

    #pragma omp simd aligned(g_yyz_x_zz_xx, g_yyz_x_zz_xy, g_yyz_x_zz_xz, g_yyz_x_zz_yy, g_yyz_x_zz_yz, g_yyz_x_zz_zz, g_yz_0_0_0_y_x_zz_xx, g_yz_0_0_0_y_x_zz_xy, g_yz_0_0_0_y_x_zz_xz, g_yz_0_0_0_y_x_zz_yy, g_yz_0_0_0_y_x_zz_yz, g_yz_0_0_0_y_x_zz_zz, g_z_x_zz_xx, g_z_x_zz_xy, g_z_x_zz_xz, g_z_x_zz_yy, g_z_x_zz_yz, g_z_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_x_zz_xx[i] = -2.0 * g_z_x_zz_xx[i] * a_exp + 4.0 * g_yyz_x_zz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_zz_xy[i] = -2.0 * g_z_x_zz_xy[i] * a_exp + 4.0 * g_yyz_x_zz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_zz_xz[i] = -2.0 * g_z_x_zz_xz[i] * a_exp + 4.0 * g_yyz_x_zz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_zz_yy[i] = -2.0 * g_z_x_zz_yy[i] * a_exp + 4.0 * g_yyz_x_zz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_zz_yz[i] = -2.0 * g_z_x_zz_yz[i] * a_exp + 4.0 * g_yyz_x_zz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_zz_zz[i] = -2.0 * g_z_x_zz_zz[i] * a_exp + 4.0 * g_yyz_x_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1440-1446)

    #pragma omp simd aligned(g_yyz_y_xx_xx, g_yyz_y_xx_xy, g_yyz_y_xx_xz, g_yyz_y_xx_yy, g_yyz_y_xx_yz, g_yyz_y_xx_zz, g_yz_0_0_0_y_y_xx_xx, g_yz_0_0_0_y_y_xx_xy, g_yz_0_0_0_y_y_xx_xz, g_yz_0_0_0_y_y_xx_yy, g_yz_0_0_0_y_y_xx_yz, g_yz_0_0_0_y_y_xx_zz, g_z_y_xx_xx, g_z_y_xx_xy, g_z_y_xx_xz, g_z_y_xx_yy, g_z_y_xx_yz, g_z_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_y_xx_xx[i] = -2.0 * g_z_y_xx_xx[i] * a_exp + 4.0 * g_yyz_y_xx_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_xx_xy[i] = -2.0 * g_z_y_xx_xy[i] * a_exp + 4.0 * g_yyz_y_xx_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_xx_xz[i] = -2.0 * g_z_y_xx_xz[i] * a_exp + 4.0 * g_yyz_y_xx_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_xx_yy[i] = -2.0 * g_z_y_xx_yy[i] * a_exp + 4.0 * g_yyz_y_xx_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_xx_yz[i] = -2.0 * g_z_y_xx_yz[i] * a_exp + 4.0 * g_yyz_y_xx_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_xx_zz[i] = -2.0 * g_z_y_xx_zz[i] * a_exp + 4.0 * g_yyz_y_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (1446-1452)

    #pragma omp simd aligned(g_yyz_y_xy_xx, g_yyz_y_xy_xy, g_yyz_y_xy_xz, g_yyz_y_xy_yy, g_yyz_y_xy_yz, g_yyz_y_xy_zz, g_yz_0_0_0_y_y_xy_xx, g_yz_0_0_0_y_y_xy_xy, g_yz_0_0_0_y_y_xy_xz, g_yz_0_0_0_y_y_xy_yy, g_yz_0_0_0_y_y_xy_yz, g_yz_0_0_0_y_y_xy_zz, g_z_y_xy_xx, g_z_y_xy_xy, g_z_y_xy_xz, g_z_y_xy_yy, g_z_y_xy_yz, g_z_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_y_xy_xx[i] = -2.0 * g_z_y_xy_xx[i] * a_exp + 4.0 * g_yyz_y_xy_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_xy_xy[i] = -2.0 * g_z_y_xy_xy[i] * a_exp + 4.0 * g_yyz_y_xy_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_xy_xz[i] = -2.0 * g_z_y_xy_xz[i] * a_exp + 4.0 * g_yyz_y_xy_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_xy_yy[i] = -2.0 * g_z_y_xy_yy[i] * a_exp + 4.0 * g_yyz_y_xy_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_xy_yz[i] = -2.0 * g_z_y_xy_yz[i] * a_exp + 4.0 * g_yyz_y_xy_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_xy_zz[i] = -2.0 * g_z_y_xy_zz[i] * a_exp + 4.0 * g_yyz_y_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1452-1458)

    #pragma omp simd aligned(g_yyz_y_xz_xx, g_yyz_y_xz_xy, g_yyz_y_xz_xz, g_yyz_y_xz_yy, g_yyz_y_xz_yz, g_yyz_y_xz_zz, g_yz_0_0_0_y_y_xz_xx, g_yz_0_0_0_y_y_xz_xy, g_yz_0_0_0_y_y_xz_xz, g_yz_0_0_0_y_y_xz_yy, g_yz_0_0_0_y_y_xz_yz, g_yz_0_0_0_y_y_xz_zz, g_z_y_xz_xx, g_z_y_xz_xy, g_z_y_xz_xz, g_z_y_xz_yy, g_z_y_xz_yz, g_z_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_y_xz_xx[i] = -2.0 * g_z_y_xz_xx[i] * a_exp + 4.0 * g_yyz_y_xz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_xz_xy[i] = -2.0 * g_z_y_xz_xy[i] * a_exp + 4.0 * g_yyz_y_xz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_xz_xz[i] = -2.0 * g_z_y_xz_xz[i] * a_exp + 4.0 * g_yyz_y_xz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_xz_yy[i] = -2.0 * g_z_y_xz_yy[i] * a_exp + 4.0 * g_yyz_y_xz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_xz_yz[i] = -2.0 * g_z_y_xz_yz[i] * a_exp + 4.0 * g_yyz_y_xz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_xz_zz[i] = -2.0 * g_z_y_xz_zz[i] * a_exp + 4.0 * g_yyz_y_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1458-1464)

    #pragma omp simd aligned(g_yyz_y_yy_xx, g_yyz_y_yy_xy, g_yyz_y_yy_xz, g_yyz_y_yy_yy, g_yyz_y_yy_yz, g_yyz_y_yy_zz, g_yz_0_0_0_y_y_yy_xx, g_yz_0_0_0_y_y_yy_xy, g_yz_0_0_0_y_y_yy_xz, g_yz_0_0_0_y_y_yy_yy, g_yz_0_0_0_y_y_yy_yz, g_yz_0_0_0_y_y_yy_zz, g_z_y_yy_xx, g_z_y_yy_xy, g_z_y_yy_xz, g_z_y_yy_yy, g_z_y_yy_yz, g_z_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_y_yy_xx[i] = -2.0 * g_z_y_yy_xx[i] * a_exp + 4.0 * g_yyz_y_yy_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_yy_xy[i] = -2.0 * g_z_y_yy_xy[i] * a_exp + 4.0 * g_yyz_y_yy_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_yy_xz[i] = -2.0 * g_z_y_yy_xz[i] * a_exp + 4.0 * g_yyz_y_yy_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_yy_yy[i] = -2.0 * g_z_y_yy_yy[i] * a_exp + 4.0 * g_yyz_y_yy_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_yy_yz[i] = -2.0 * g_z_y_yy_yz[i] * a_exp + 4.0 * g_yyz_y_yy_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_yy_zz[i] = -2.0 * g_z_y_yy_zz[i] * a_exp + 4.0 * g_yyz_y_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1464-1470)

    #pragma omp simd aligned(g_yyz_y_yz_xx, g_yyz_y_yz_xy, g_yyz_y_yz_xz, g_yyz_y_yz_yy, g_yyz_y_yz_yz, g_yyz_y_yz_zz, g_yz_0_0_0_y_y_yz_xx, g_yz_0_0_0_y_y_yz_xy, g_yz_0_0_0_y_y_yz_xz, g_yz_0_0_0_y_y_yz_yy, g_yz_0_0_0_y_y_yz_yz, g_yz_0_0_0_y_y_yz_zz, g_z_y_yz_xx, g_z_y_yz_xy, g_z_y_yz_xz, g_z_y_yz_yy, g_z_y_yz_yz, g_z_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_y_yz_xx[i] = -2.0 * g_z_y_yz_xx[i] * a_exp + 4.0 * g_yyz_y_yz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_yz_xy[i] = -2.0 * g_z_y_yz_xy[i] * a_exp + 4.0 * g_yyz_y_yz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_yz_xz[i] = -2.0 * g_z_y_yz_xz[i] * a_exp + 4.0 * g_yyz_y_yz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_yz_yy[i] = -2.0 * g_z_y_yz_yy[i] * a_exp + 4.0 * g_yyz_y_yz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_yz_yz[i] = -2.0 * g_z_y_yz_yz[i] * a_exp + 4.0 * g_yyz_y_yz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_yz_zz[i] = -2.0 * g_z_y_yz_zz[i] * a_exp + 4.0 * g_yyz_y_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1470-1476)

    #pragma omp simd aligned(g_yyz_y_zz_xx, g_yyz_y_zz_xy, g_yyz_y_zz_xz, g_yyz_y_zz_yy, g_yyz_y_zz_yz, g_yyz_y_zz_zz, g_yz_0_0_0_y_y_zz_xx, g_yz_0_0_0_y_y_zz_xy, g_yz_0_0_0_y_y_zz_xz, g_yz_0_0_0_y_y_zz_yy, g_yz_0_0_0_y_y_zz_yz, g_yz_0_0_0_y_y_zz_zz, g_z_y_zz_xx, g_z_y_zz_xy, g_z_y_zz_xz, g_z_y_zz_yy, g_z_y_zz_yz, g_z_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_y_zz_xx[i] = -2.0 * g_z_y_zz_xx[i] * a_exp + 4.0 * g_yyz_y_zz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_zz_xy[i] = -2.0 * g_z_y_zz_xy[i] * a_exp + 4.0 * g_yyz_y_zz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_zz_xz[i] = -2.0 * g_z_y_zz_xz[i] * a_exp + 4.0 * g_yyz_y_zz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_zz_yy[i] = -2.0 * g_z_y_zz_yy[i] * a_exp + 4.0 * g_yyz_y_zz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_zz_yz[i] = -2.0 * g_z_y_zz_yz[i] * a_exp + 4.0 * g_yyz_y_zz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_zz_zz[i] = -2.0 * g_z_y_zz_zz[i] * a_exp + 4.0 * g_yyz_y_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1476-1482)

    #pragma omp simd aligned(g_yyz_z_xx_xx, g_yyz_z_xx_xy, g_yyz_z_xx_xz, g_yyz_z_xx_yy, g_yyz_z_xx_yz, g_yyz_z_xx_zz, g_yz_0_0_0_y_z_xx_xx, g_yz_0_0_0_y_z_xx_xy, g_yz_0_0_0_y_z_xx_xz, g_yz_0_0_0_y_z_xx_yy, g_yz_0_0_0_y_z_xx_yz, g_yz_0_0_0_y_z_xx_zz, g_z_z_xx_xx, g_z_z_xx_xy, g_z_z_xx_xz, g_z_z_xx_yy, g_z_z_xx_yz, g_z_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_z_xx_xx[i] = -2.0 * g_z_z_xx_xx[i] * a_exp + 4.0 * g_yyz_z_xx_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_xx_xy[i] = -2.0 * g_z_z_xx_xy[i] * a_exp + 4.0 * g_yyz_z_xx_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_xx_xz[i] = -2.0 * g_z_z_xx_xz[i] * a_exp + 4.0 * g_yyz_z_xx_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_xx_yy[i] = -2.0 * g_z_z_xx_yy[i] * a_exp + 4.0 * g_yyz_z_xx_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_xx_yz[i] = -2.0 * g_z_z_xx_yz[i] * a_exp + 4.0 * g_yyz_z_xx_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_xx_zz[i] = -2.0 * g_z_z_xx_zz[i] * a_exp + 4.0 * g_yyz_z_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (1482-1488)

    #pragma omp simd aligned(g_yyz_z_xy_xx, g_yyz_z_xy_xy, g_yyz_z_xy_xz, g_yyz_z_xy_yy, g_yyz_z_xy_yz, g_yyz_z_xy_zz, g_yz_0_0_0_y_z_xy_xx, g_yz_0_0_0_y_z_xy_xy, g_yz_0_0_0_y_z_xy_xz, g_yz_0_0_0_y_z_xy_yy, g_yz_0_0_0_y_z_xy_yz, g_yz_0_0_0_y_z_xy_zz, g_z_z_xy_xx, g_z_z_xy_xy, g_z_z_xy_xz, g_z_z_xy_yy, g_z_z_xy_yz, g_z_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_z_xy_xx[i] = -2.0 * g_z_z_xy_xx[i] * a_exp + 4.0 * g_yyz_z_xy_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_xy_xy[i] = -2.0 * g_z_z_xy_xy[i] * a_exp + 4.0 * g_yyz_z_xy_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_xy_xz[i] = -2.0 * g_z_z_xy_xz[i] * a_exp + 4.0 * g_yyz_z_xy_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_xy_yy[i] = -2.0 * g_z_z_xy_yy[i] * a_exp + 4.0 * g_yyz_z_xy_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_xy_yz[i] = -2.0 * g_z_z_xy_yz[i] * a_exp + 4.0 * g_yyz_z_xy_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_xy_zz[i] = -2.0 * g_z_z_xy_zz[i] * a_exp + 4.0 * g_yyz_z_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1488-1494)

    #pragma omp simd aligned(g_yyz_z_xz_xx, g_yyz_z_xz_xy, g_yyz_z_xz_xz, g_yyz_z_xz_yy, g_yyz_z_xz_yz, g_yyz_z_xz_zz, g_yz_0_0_0_y_z_xz_xx, g_yz_0_0_0_y_z_xz_xy, g_yz_0_0_0_y_z_xz_xz, g_yz_0_0_0_y_z_xz_yy, g_yz_0_0_0_y_z_xz_yz, g_yz_0_0_0_y_z_xz_zz, g_z_z_xz_xx, g_z_z_xz_xy, g_z_z_xz_xz, g_z_z_xz_yy, g_z_z_xz_yz, g_z_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_z_xz_xx[i] = -2.0 * g_z_z_xz_xx[i] * a_exp + 4.0 * g_yyz_z_xz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_xz_xy[i] = -2.0 * g_z_z_xz_xy[i] * a_exp + 4.0 * g_yyz_z_xz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_xz_xz[i] = -2.0 * g_z_z_xz_xz[i] * a_exp + 4.0 * g_yyz_z_xz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_xz_yy[i] = -2.0 * g_z_z_xz_yy[i] * a_exp + 4.0 * g_yyz_z_xz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_xz_yz[i] = -2.0 * g_z_z_xz_yz[i] * a_exp + 4.0 * g_yyz_z_xz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_xz_zz[i] = -2.0 * g_z_z_xz_zz[i] * a_exp + 4.0 * g_yyz_z_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1494-1500)

    #pragma omp simd aligned(g_yyz_z_yy_xx, g_yyz_z_yy_xy, g_yyz_z_yy_xz, g_yyz_z_yy_yy, g_yyz_z_yy_yz, g_yyz_z_yy_zz, g_yz_0_0_0_y_z_yy_xx, g_yz_0_0_0_y_z_yy_xy, g_yz_0_0_0_y_z_yy_xz, g_yz_0_0_0_y_z_yy_yy, g_yz_0_0_0_y_z_yy_yz, g_yz_0_0_0_y_z_yy_zz, g_z_z_yy_xx, g_z_z_yy_xy, g_z_z_yy_xz, g_z_z_yy_yy, g_z_z_yy_yz, g_z_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_z_yy_xx[i] = -2.0 * g_z_z_yy_xx[i] * a_exp + 4.0 * g_yyz_z_yy_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_yy_xy[i] = -2.0 * g_z_z_yy_xy[i] * a_exp + 4.0 * g_yyz_z_yy_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_yy_xz[i] = -2.0 * g_z_z_yy_xz[i] * a_exp + 4.0 * g_yyz_z_yy_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_yy_yy[i] = -2.0 * g_z_z_yy_yy[i] * a_exp + 4.0 * g_yyz_z_yy_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_yy_yz[i] = -2.0 * g_z_z_yy_yz[i] * a_exp + 4.0 * g_yyz_z_yy_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_yy_zz[i] = -2.0 * g_z_z_yy_zz[i] * a_exp + 4.0 * g_yyz_z_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1500-1506)

    #pragma omp simd aligned(g_yyz_z_yz_xx, g_yyz_z_yz_xy, g_yyz_z_yz_xz, g_yyz_z_yz_yy, g_yyz_z_yz_yz, g_yyz_z_yz_zz, g_yz_0_0_0_y_z_yz_xx, g_yz_0_0_0_y_z_yz_xy, g_yz_0_0_0_y_z_yz_xz, g_yz_0_0_0_y_z_yz_yy, g_yz_0_0_0_y_z_yz_yz, g_yz_0_0_0_y_z_yz_zz, g_z_z_yz_xx, g_z_z_yz_xy, g_z_z_yz_xz, g_z_z_yz_yy, g_z_z_yz_yz, g_z_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_z_yz_xx[i] = -2.0 * g_z_z_yz_xx[i] * a_exp + 4.0 * g_yyz_z_yz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_yz_xy[i] = -2.0 * g_z_z_yz_xy[i] * a_exp + 4.0 * g_yyz_z_yz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_yz_xz[i] = -2.0 * g_z_z_yz_xz[i] * a_exp + 4.0 * g_yyz_z_yz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_yz_yy[i] = -2.0 * g_z_z_yz_yy[i] * a_exp + 4.0 * g_yyz_z_yz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_yz_yz[i] = -2.0 * g_z_z_yz_yz[i] * a_exp + 4.0 * g_yyz_z_yz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_yz_zz[i] = -2.0 * g_z_z_yz_zz[i] * a_exp + 4.0 * g_yyz_z_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1506-1512)

    #pragma omp simd aligned(g_yyz_z_zz_xx, g_yyz_z_zz_xy, g_yyz_z_zz_xz, g_yyz_z_zz_yy, g_yyz_z_zz_yz, g_yyz_z_zz_zz, g_yz_0_0_0_y_z_zz_xx, g_yz_0_0_0_y_z_zz_xy, g_yz_0_0_0_y_z_zz_xz, g_yz_0_0_0_y_z_zz_yy, g_yz_0_0_0_y_z_zz_yz, g_yz_0_0_0_y_z_zz_zz, g_z_z_zz_xx, g_z_z_zz_xy, g_z_z_zz_xz, g_z_z_zz_yy, g_z_z_zz_yz, g_z_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_z_zz_xx[i] = -2.0 * g_z_z_zz_xx[i] * a_exp + 4.0 * g_yyz_z_zz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_zz_xy[i] = -2.0 * g_z_z_zz_xy[i] * a_exp + 4.0 * g_yyz_z_zz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_zz_xz[i] = -2.0 * g_z_z_zz_xz[i] * a_exp + 4.0 * g_yyz_z_zz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_zz_yy[i] = -2.0 * g_z_z_zz_yy[i] * a_exp + 4.0 * g_yyz_z_zz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_zz_yz[i] = -2.0 * g_z_z_zz_yz[i] * a_exp + 4.0 * g_yyz_z_zz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_zz_zz[i] = -2.0 * g_z_z_zz_zz[i] * a_exp + 4.0 * g_yyz_z_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1512-1518)

    #pragma omp simd aligned(g_y_x_xx_xx, g_y_x_xx_xy, g_y_x_xx_xz, g_y_x_xx_yy, g_y_x_xx_yz, g_y_x_xx_zz, g_yz_0_0_0_z_x_xx_xx, g_yz_0_0_0_z_x_xx_xy, g_yz_0_0_0_z_x_xx_xz, g_yz_0_0_0_z_x_xx_yy, g_yz_0_0_0_z_x_xx_yz, g_yz_0_0_0_z_x_xx_zz, g_yzz_x_xx_xx, g_yzz_x_xx_xy, g_yzz_x_xx_xz, g_yzz_x_xx_yy, g_yzz_x_xx_yz, g_yzz_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_x_xx_xx[i] = -2.0 * g_y_x_xx_xx[i] * a_exp + 4.0 * g_yzz_x_xx_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_xx_xy[i] = -2.0 * g_y_x_xx_xy[i] * a_exp + 4.0 * g_yzz_x_xx_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_xx_xz[i] = -2.0 * g_y_x_xx_xz[i] * a_exp + 4.0 * g_yzz_x_xx_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_xx_yy[i] = -2.0 * g_y_x_xx_yy[i] * a_exp + 4.0 * g_yzz_x_xx_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_xx_yz[i] = -2.0 * g_y_x_xx_yz[i] * a_exp + 4.0 * g_yzz_x_xx_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_xx_zz[i] = -2.0 * g_y_x_xx_zz[i] * a_exp + 4.0 * g_yzz_x_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (1518-1524)

    #pragma omp simd aligned(g_y_x_xy_xx, g_y_x_xy_xy, g_y_x_xy_xz, g_y_x_xy_yy, g_y_x_xy_yz, g_y_x_xy_zz, g_yz_0_0_0_z_x_xy_xx, g_yz_0_0_0_z_x_xy_xy, g_yz_0_0_0_z_x_xy_xz, g_yz_0_0_0_z_x_xy_yy, g_yz_0_0_0_z_x_xy_yz, g_yz_0_0_0_z_x_xy_zz, g_yzz_x_xy_xx, g_yzz_x_xy_xy, g_yzz_x_xy_xz, g_yzz_x_xy_yy, g_yzz_x_xy_yz, g_yzz_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_x_xy_xx[i] = -2.0 * g_y_x_xy_xx[i] * a_exp + 4.0 * g_yzz_x_xy_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_xy_xy[i] = -2.0 * g_y_x_xy_xy[i] * a_exp + 4.0 * g_yzz_x_xy_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_xy_xz[i] = -2.0 * g_y_x_xy_xz[i] * a_exp + 4.0 * g_yzz_x_xy_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_xy_yy[i] = -2.0 * g_y_x_xy_yy[i] * a_exp + 4.0 * g_yzz_x_xy_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_xy_yz[i] = -2.0 * g_y_x_xy_yz[i] * a_exp + 4.0 * g_yzz_x_xy_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_xy_zz[i] = -2.0 * g_y_x_xy_zz[i] * a_exp + 4.0 * g_yzz_x_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1524-1530)

    #pragma omp simd aligned(g_y_x_xz_xx, g_y_x_xz_xy, g_y_x_xz_xz, g_y_x_xz_yy, g_y_x_xz_yz, g_y_x_xz_zz, g_yz_0_0_0_z_x_xz_xx, g_yz_0_0_0_z_x_xz_xy, g_yz_0_0_0_z_x_xz_xz, g_yz_0_0_0_z_x_xz_yy, g_yz_0_0_0_z_x_xz_yz, g_yz_0_0_0_z_x_xz_zz, g_yzz_x_xz_xx, g_yzz_x_xz_xy, g_yzz_x_xz_xz, g_yzz_x_xz_yy, g_yzz_x_xz_yz, g_yzz_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_x_xz_xx[i] = -2.0 * g_y_x_xz_xx[i] * a_exp + 4.0 * g_yzz_x_xz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_xz_xy[i] = -2.0 * g_y_x_xz_xy[i] * a_exp + 4.0 * g_yzz_x_xz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_xz_xz[i] = -2.0 * g_y_x_xz_xz[i] * a_exp + 4.0 * g_yzz_x_xz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_xz_yy[i] = -2.0 * g_y_x_xz_yy[i] * a_exp + 4.0 * g_yzz_x_xz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_xz_yz[i] = -2.0 * g_y_x_xz_yz[i] * a_exp + 4.0 * g_yzz_x_xz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_xz_zz[i] = -2.0 * g_y_x_xz_zz[i] * a_exp + 4.0 * g_yzz_x_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1530-1536)

    #pragma omp simd aligned(g_y_x_yy_xx, g_y_x_yy_xy, g_y_x_yy_xz, g_y_x_yy_yy, g_y_x_yy_yz, g_y_x_yy_zz, g_yz_0_0_0_z_x_yy_xx, g_yz_0_0_0_z_x_yy_xy, g_yz_0_0_0_z_x_yy_xz, g_yz_0_0_0_z_x_yy_yy, g_yz_0_0_0_z_x_yy_yz, g_yz_0_0_0_z_x_yy_zz, g_yzz_x_yy_xx, g_yzz_x_yy_xy, g_yzz_x_yy_xz, g_yzz_x_yy_yy, g_yzz_x_yy_yz, g_yzz_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_x_yy_xx[i] = -2.0 * g_y_x_yy_xx[i] * a_exp + 4.0 * g_yzz_x_yy_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_yy_xy[i] = -2.0 * g_y_x_yy_xy[i] * a_exp + 4.0 * g_yzz_x_yy_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_yy_xz[i] = -2.0 * g_y_x_yy_xz[i] * a_exp + 4.0 * g_yzz_x_yy_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_yy_yy[i] = -2.0 * g_y_x_yy_yy[i] * a_exp + 4.0 * g_yzz_x_yy_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_yy_yz[i] = -2.0 * g_y_x_yy_yz[i] * a_exp + 4.0 * g_yzz_x_yy_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_yy_zz[i] = -2.0 * g_y_x_yy_zz[i] * a_exp + 4.0 * g_yzz_x_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1536-1542)

    #pragma omp simd aligned(g_y_x_yz_xx, g_y_x_yz_xy, g_y_x_yz_xz, g_y_x_yz_yy, g_y_x_yz_yz, g_y_x_yz_zz, g_yz_0_0_0_z_x_yz_xx, g_yz_0_0_0_z_x_yz_xy, g_yz_0_0_0_z_x_yz_xz, g_yz_0_0_0_z_x_yz_yy, g_yz_0_0_0_z_x_yz_yz, g_yz_0_0_0_z_x_yz_zz, g_yzz_x_yz_xx, g_yzz_x_yz_xy, g_yzz_x_yz_xz, g_yzz_x_yz_yy, g_yzz_x_yz_yz, g_yzz_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_x_yz_xx[i] = -2.0 * g_y_x_yz_xx[i] * a_exp + 4.0 * g_yzz_x_yz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_yz_xy[i] = -2.0 * g_y_x_yz_xy[i] * a_exp + 4.0 * g_yzz_x_yz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_yz_xz[i] = -2.0 * g_y_x_yz_xz[i] * a_exp + 4.0 * g_yzz_x_yz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_yz_yy[i] = -2.0 * g_y_x_yz_yy[i] * a_exp + 4.0 * g_yzz_x_yz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_yz_yz[i] = -2.0 * g_y_x_yz_yz[i] * a_exp + 4.0 * g_yzz_x_yz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_yz_zz[i] = -2.0 * g_y_x_yz_zz[i] * a_exp + 4.0 * g_yzz_x_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1542-1548)

    #pragma omp simd aligned(g_y_x_zz_xx, g_y_x_zz_xy, g_y_x_zz_xz, g_y_x_zz_yy, g_y_x_zz_yz, g_y_x_zz_zz, g_yz_0_0_0_z_x_zz_xx, g_yz_0_0_0_z_x_zz_xy, g_yz_0_0_0_z_x_zz_xz, g_yz_0_0_0_z_x_zz_yy, g_yz_0_0_0_z_x_zz_yz, g_yz_0_0_0_z_x_zz_zz, g_yzz_x_zz_xx, g_yzz_x_zz_xy, g_yzz_x_zz_xz, g_yzz_x_zz_yy, g_yzz_x_zz_yz, g_yzz_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_x_zz_xx[i] = -2.0 * g_y_x_zz_xx[i] * a_exp + 4.0 * g_yzz_x_zz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_zz_xy[i] = -2.0 * g_y_x_zz_xy[i] * a_exp + 4.0 * g_yzz_x_zz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_zz_xz[i] = -2.0 * g_y_x_zz_xz[i] * a_exp + 4.0 * g_yzz_x_zz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_zz_yy[i] = -2.0 * g_y_x_zz_yy[i] * a_exp + 4.0 * g_yzz_x_zz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_zz_yz[i] = -2.0 * g_y_x_zz_yz[i] * a_exp + 4.0 * g_yzz_x_zz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_zz_zz[i] = -2.0 * g_y_x_zz_zz[i] * a_exp + 4.0 * g_yzz_x_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1548-1554)

    #pragma omp simd aligned(g_y_y_xx_xx, g_y_y_xx_xy, g_y_y_xx_xz, g_y_y_xx_yy, g_y_y_xx_yz, g_y_y_xx_zz, g_yz_0_0_0_z_y_xx_xx, g_yz_0_0_0_z_y_xx_xy, g_yz_0_0_0_z_y_xx_xz, g_yz_0_0_0_z_y_xx_yy, g_yz_0_0_0_z_y_xx_yz, g_yz_0_0_0_z_y_xx_zz, g_yzz_y_xx_xx, g_yzz_y_xx_xy, g_yzz_y_xx_xz, g_yzz_y_xx_yy, g_yzz_y_xx_yz, g_yzz_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_y_xx_xx[i] = -2.0 * g_y_y_xx_xx[i] * a_exp + 4.0 * g_yzz_y_xx_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_xx_xy[i] = -2.0 * g_y_y_xx_xy[i] * a_exp + 4.0 * g_yzz_y_xx_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_xx_xz[i] = -2.0 * g_y_y_xx_xz[i] * a_exp + 4.0 * g_yzz_y_xx_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_xx_yy[i] = -2.0 * g_y_y_xx_yy[i] * a_exp + 4.0 * g_yzz_y_xx_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_xx_yz[i] = -2.0 * g_y_y_xx_yz[i] * a_exp + 4.0 * g_yzz_y_xx_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_xx_zz[i] = -2.0 * g_y_y_xx_zz[i] * a_exp + 4.0 * g_yzz_y_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (1554-1560)

    #pragma omp simd aligned(g_y_y_xy_xx, g_y_y_xy_xy, g_y_y_xy_xz, g_y_y_xy_yy, g_y_y_xy_yz, g_y_y_xy_zz, g_yz_0_0_0_z_y_xy_xx, g_yz_0_0_0_z_y_xy_xy, g_yz_0_0_0_z_y_xy_xz, g_yz_0_0_0_z_y_xy_yy, g_yz_0_0_0_z_y_xy_yz, g_yz_0_0_0_z_y_xy_zz, g_yzz_y_xy_xx, g_yzz_y_xy_xy, g_yzz_y_xy_xz, g_yzz_y_xy_yy, g_yzz_y_xy_yz, g_yzz_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_y_xy_xx[i] = -2.0 * g_y_y_xy_xx[i] * a_exp + 4.0 * g_yzz_y_xy_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_xy_xy[i] = -2.0 * g_y_y_xy_xy[i] * a_exp + 4.0 * g_yzz_y_xy_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_xy_xz[i] = -2.0 * g_y_y_xy_xz[i] * a_exp + 4.0 * g_yzz_y_xy_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_xy_yy[i] = -2.0 * g_y_y_xy_yy[i] * a_exp + 4.0 * g_yzz_y_xy_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_xy_yz[i] = -2.0 * g_y_y_xy_yz[i] * a_exp + 4.0 * g_yzz_y_xy_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_xy_zz[i] = -2.0 * g_y_y_xy_zz[i] * a_exp + 4.0 * g_yzz_y_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1560-1566)

    #pragma omp simd aligned(g_y_y_xz_xx, g_y_y_xz_xy, g_y_y_xz_xz, g_y_y_xz_yy, g_y_y_xz_yz, g_y_y_xz_zz, g_yz_0_0_0_z_y_xz_xx, g_yz_0_0_0_z_y_xz_xy, g_yz_0_0_0_z_y_xz_xz, g_yz_0_0_0_z_y_xz_yy, g_yz_0_0_0_z_y_xz_yz, g_yz_0_0_0_z_y_xz_zz, g_yzz_y_xz_xx, g_yzz_y_xz_xy, g_yzz_y_xz_xz, g_yzz_y_xz_yy, g_yzz_y_xz_yz, g_yzz_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_y_xz_xx[i] = -2.0 * g_y_y_xz_xx[i] * a_exp + 4.0 * g_yzz_y_xz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_xz_xy[i] = -2.0 * g_y_y_xz_xy[i] * a_exp + 4.0 * g_yzz_y_xz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_xz_xz[i] = -2.0 * g_y_y_xz_xz[i] * a_exp + 4.0 * g_yzz_y_xz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_xz_yy[i] = -2.0 * g_y_y_xz_yy[i] * a_exp + 4.0 * g_yzz_y_xz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_xz_yz[i] = -2.0 * g_y_y_xz_yz[i] * a_exp + 4.0 * g_yzz_y_xz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_xz_zz[i] = -2.0 * g_y_y_xz_zz[i] * a_exp + 4.0 * g_yzz_y_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1566-1572)

    #pragma omp simd aligned(g_y_y_yy_xx, g_y_y_yy_xy, g_y_y_yy_xz, g_y_y_yy_yy, g_y_y_yy_yz, g_y_y_yy_zz, g_yz_0_0_0_z_y_yy_xx, g_yz_0_0_0_z_y_yy_xy, g_yz_0_0_0_z_y_yy_xz, g_yz_0_0_0_z_y_yy_yy, g_yz_0_0_0_z_y_yy_yz, g_yz_0_0_0_z_y_yy_zz, g_yzz_y_yy_xx, g_yzz_y_yy_xy, g_yzz_y_yy_xz, g_yzz_y_yy_yy, g_yzz_y_yy_yz, g_yzz_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_y_yy_xx[i] = -2.0 * g_y_y_yy_xx[i] * a_exp + 4.0 * g_yzz_y_yy_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_yy_xy[i] = -2.0 * g_y_y_yy_xy[i] * a_exp + 4.0 * g_yzz_y_yy_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_yy_xz[i] = -2.0 * g_y_y_yy_xz[i] * a_exp + 4.0 * g_yzz_y_yy_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_yy_yy[i] = -2.0 * g_y_y_yy_yy[i] * a_exp + 4.0 * g_yzz_y_yy_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_yy_yz[i] = -2.0 * g_y_y_yy_yz[i] * a_exp + 4.0 * g_yzz_y_yy_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_yy_zz[i] = -2.0 * g_y_y_yy_zz[i] * a_exp + 4.0 * g_yzz_y_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1572-1578)

    #pragma omp simd aligned(g_y_y_yz_xx, g_y_y_yz_xy, g_y_y_yz_xz, g_y_y_yz_yy, g_y_y_yz_yz, g_y_y_yz_zz, g_yz_0_0_0_z_y_yz_xx, g_yz_0_0_0_z_y_yz_xy, g_yz_0_0_0_z_y_yz_xz, g_yz_0_0_0_z_y_yz_yy, g_yz_0_0_0_z_y_yz_yz, g_yz_0_0_0_z_y_yz_zz, g_yzz_y_yz_xx, g_yzz_y_yz_xy, g_yzz_y_yz_xz, g_yzz_y_yz_yy, g_yzz_y_yz_yz, g_yzz_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_y_yz_xx[i] = -2.0 * g_y_y_yz_xx[i] * a_exp + 4.0 * g_yzz_y_yz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_yz_xy[i] = -2.0 * g_y_y_yz_xy[i] * a_exp + 4.0 * g_yzz_y_yz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_yz_xz[i] = -2.0 * g_y_y_yz_xz[i] * a_exp + 4.0 * g_yzz_y_yz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_yz_yy[i] = -2.0 * g_y_y_yz_yy[i] * a_exp + 4.0 * g_yzz_y_yz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_yz_yz[i] = -2.0 * g_y_y_yz_yz[i] * a_exp + 4.0 * g_yzz_y_yz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_yz_zz[i] = -2.0 * g_y_y_yz_zz[i] * a_exp + 4.0 * g_yzz_y_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1578-1584)

    #pragma omp simd aligned(g_y_y_zz_xx, g_y_y_zz_xy, g_y_y_zz_xz, g_y_y_zz_yy, g_y_y_zz_yz, g_y_y_zz_zz, g_yz_0_0_0_z_y_zz_xx, g_yz_0_0_0_z_y_zz_xy, g_yz_0_0_0_z_y_zz_xz, g_yz_0_0_0_z_y_zz_yy, g_yz_0_0_0_z_y_zz_yz, g_yz_0_0_0_z_y_zz_zz, g_yzz_y_zz_xx, g_yzz_y_zz_xy, g_yzz_y_zz_xz, g_yzz_y_zz_yy, g_yzz_y_zz_yz, g_yzz_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_y_zz_xx[i] = -2.0 * g_y_y_zz_xx[i] * a_exp + 4.0 * g_yzz_y_zz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_zz_xy[i] = -2.0 * g_y_y_zz_xy[i] * a_exp + 4.0 * g_yzz_y_zz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_zz_xz[i] = -2.0 * g_y_y_zz_xz[i] * a_exp + 4.0 * g_yzz_y_zz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_zz_yy[i] = -2.0 * g_y_y_zz_yy[i] * a_exp + 4.0 * g_yzz_y_zz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_zz_yz[i] = -2.0 * g_y_y_zz_yz[i] * a_exp + 4.0 * g_yzz_y_zz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_zz_zz[i] = -2.0 * g_y_y_zz_zz[i] * a_exp + 4.0 * g_yzz_y_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1584-1590)

    #pragma omp simd aligned(g_y_z_xx_xx, g_y_z_xx_xy, g_y_z_xx_xz, g_y_z_xx_yy, g_y_z_xx_yz, g_y_z_xx_zz, g_yz_0_0_0_z_z_xx_xx, g_yz_0_0_0_z_z_xx_xy, g_yz_0_0_0_z_z_xx_xz, g_yz_0_0_0_z_z_xx_yy, g_yz_0_0_0_z_z_xx_yz, g_yz_0_0_0_z_z_xx_zz, g_yzz_z_xx_xx, g_yzz_z_xx_xy, g_yzz_z_xx_xz, g_yzz_z_xx_yy, g_yzz_z_xx_yz, g_yzz_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_z_xx_xx[i] = -2.0 * g_y_z_xx_xx[i] * a_exp + 4.0 * g_yzz_z_xx_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_xx_xy[i] = -2.0 * g_y_z_xx_xy[i] * a_exp + 4.0 * g_yzz_z_xx_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_xx_xz[i] = -2.0 * g_y_z_xx_xz[i] * a_exp + 4.0 * g_yzz_z_xx_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_xx_yy[i] = -2.0 * g_y_z_xx_yy[i] * a_exp + 4.0 * g_yzz_z_xx_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_xx_yz[i] = -2.0 * g_y_z_xx_yz[i] * a_exp + 4.0 * g_yzz_z_xx_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_xx_zz[i] = -2.0 * g_y_z_xx_zz[i] * a_exp + 4.0 * g_yzz_z_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (1590-1596)

    #pragma omp simd aligned(g_y_z_xy_xx, g_y_z_xy_xy, g_y_z_xy_xz, g_y_z_xy_yy, g_y_z_xy_yz, g_y_z_xy_zz, g_yz_0_0_0_z_z_xy_xx, g_yz_0_0_0_z_z_xy_xy, g_yz_0_0_0_z_z_xy_xz, g_yz_0_0_0_z_z_xy_yy, g_yz_0_0_0_z_z_xy_yz, g_yz_0_0_0_z_z_xy_zz, g_yzz_z_xy_xx, g_yzz_z_xy_xy, g_yzz_z_xy_xz, g_yzz_z_xy_yy, g_yzz_z_xy_yz, g_yzz_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_z_xy_xx[i] = -2.0 * g_y_z_xy_xx[i] * a_exp + 4.0 * g_yzz_z_xy_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_xy_xy[i] = -2.0 * g_y_z_xy_xy[i] * a_exp + 4.0 * g_yzz_z_xy_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_xy_xz[i] = -2.0 * g_y_z_xy_xz[i] * a_exp + 4.0 * g_yzz_z_xy_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_xy_yy[i] = -2.0 * g_y_z_xy_yy[i] * a_exp + 4.0 * g_yzz_z_xy_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_xy_yz[i] = -2.0 * g_y_z_xy_yz[i] * a_exp + 4.0 * g_yzz_z_xy_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_xy_zz[i] = -2.0 * g_y_z_xy_zz[i] * a_exp + 4.0 * g_yzz_z_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1596-1602)

    #pragma omp simd aligned(g_y_z_xz_xx, g_y_z_xz_xy, g_y_z_xz_xz, g_y_z_xz_yy, g_y_z_xz_yz, g_y_z_xz_zz, g_yz_0_0_0_z_z_xz_xx, g_yz_0_0_0_z_z_xz_xy, g_yz_0_0_0_z_z_xz_xz, g_yz_0_0_0_z_z_xz_yy, g_yz_0_0_0_z_z_xz_yz, g_yz_0_0_0_z_z_xz_zz, g_yzz_z_xz_xx, g_yzz_z_xz_xy, g_yzz_z_xz_xz, g_yzz_z_xz_yy, g_yzz_z_xz_yz, g_yzz_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_z_xz_xx[i] = -2.0 * g_y_z_xz_xx[i] * a_exp + 4.0 * g_yzz_z_xz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_xz_xy[i] = -2.0 * g_y_z_xz_xy[i] * a_exp + 4.0 * g_yzz_z_xz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_xz_xz[i] = -2.0 * g_y_z_xz_xz[i] * a_exp + 4.0 * g_yzz_z_xz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_xz_yy[i] = -2.0 * g_y_z_xz_yy[i] * a_exp + 4.0 * g_yzz_z_xz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_xz_yz[i] = -2.0 * g_y_z_xz_yz[i] * a_exp + 4.0 * g_yzz_z_xz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_xz_zz[i] = -2.0 * g_y_z_xz_zz[i] * a_exp + 4.0 * g_yzz_z_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1602-1608)

    #pragma omp simd aligned(g_y_z_yy_xx, g_y_z_yy_xy, g_y_z_yy_xz, g_y_z_yy_yy, g_y_z_yy_yz, g_y_z_yy_zz, g_yz_0_0_0_z_z_yy_xx, g_yz_0_0_0_z_z_yy_xy, g_yz_0_0_0_z_z_yy_xz, g_yz_0_0_0_z_z_yy_yy, g_yz_0_0_0_z_z_yy_yz, g_yz_0_0_0_z_z_yy_zz, g_yzz_z_yy_xx, g_yzz_z_yy_xy, g_yzz_z_yy_xz, g_yzz_z_yy_yy, g_yzz_z_yy_yz, g_yzz_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_z_yy_xx[i] = -2.0 * g_y_z_yy_xx[i] * a_exp + 4.0 * g_yzz_z_yy_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_yy_xy[i] = -2.0 * g_y_z_yy_xy[i] * a_exp + 4.0 * g_yzz_z_yy_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_yy_xz[i] = -2.0 * g_y_z_yy_xz[i] * a_exp + 4.0 * g_yzz_z_yy_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_yy_yy[i] = -2.0 * g_y_z_yy_yy[i] * a_exp + 4.0 * g_yzz_z_yy_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_yy_yz[i] = -2.0 * g_y_z_yy_yz[i] * a_exp + 4.0 * g_yzz_z_yy_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_yy_zz[i] = -2.0 * g_y_z_yy_zz[i] * a_exp + 4.0 * g_yzz_z_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1608-1614)

    #pragma omp simd aligned(g_y_z_yz_xx, g_y_z_yz_xy, g_y_z_yz_xz, g_y_z_yz_yy, g_y_z_yz_yz, g_y_z_yz_zz, g_yz_0_0_0_z_z_yz_xx, g_yz_0_0_0_z_z_yz_xy, g_yz_0_0_0_z_z_yz_xz, g_yz_0_0_0_z_z_yz_yy, g_yz_0_0_0_z_z_yz_yz, g_yz_0_0_0_z_z_yz_zz, g_yzz_z_yz_xx, g_yzz_z_yz_xy, g_yzz_z_yz_xz, g_yzz_z_yz_yy, g_yzz_z_yz_yz, g_yzz_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_z_yz_xx[i] = -2.0 * g_y_z_yz_xx[i] * a_exp + 4.0 * g_yzz_z_yz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_yz_xy[i] = -2.0 * g_y_z_yz_xy[i] * a_exp + 4.0 * g_yzz_z_yz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_yz_xz[i] = -2.0 * g_y_z_yz_xz[i] * a_exp + 4.0 * g_yzz_z_yz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_yz_yy[i] = -2.0 * g_y_z_yz_yy[i] * a_exp + 4.0 * g_yzz_z_yz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_yz_yz[i] = -2.0 * g_y_z_yz_yz[i] * a_exp + 4.0 * g_yzz_z_yz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_yz_zz[i] = -2.0 * g_y_z_yz_zz[i] * a_exp + 4.0 * g_yzz_z_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1614-1620)

    #pragma omp simd aligned(g_y_z_zz_xx, g_y_z_zz_xy, g_y_z_zz_xz, g_y_z_zz_yy, g_y_z_zz_yz, g_y_z_zz_zz, g_yz_0_0_0_z_z_zz_xx, g_yz_0_0_0_z_z_zz_xy, g_yz_0_0_0_z_z_zz_xz, g_yz_0_0_0_z_z_zz_yy, g_yz_0_0_0_z_z_zz_yz, g_yz_0_0_0_z_z_zz_zz, g_yzz_z_zz_xx, g_yzz_z_zz_xy, g_yzz_z_zz_xz, g_yzz_z_zz_yy, g_yzz_z_zz_yz, g_yzz_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_z_zz_xx[i] = -2.0 * g_y_z_zz_xx[i] * a_exp + 4.0 * g_yzz_z_zz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_zz_xy[i] = -2.0 * g_y_z_zz_xy[i] * a_exp + 4.0 * g_yzz_z_zz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_zz_xz[i] = -2.0 * g_y_z_zz_xz[i] * a_exp + 4.0 * g_yzz_z_zz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_zz_yy[i] = -2.0 * g_y_z_zz_yy[i] * a_exp + 4.0 * g_yzz_z_zz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_zz_yz[i] = -2.0 * g_y_z_zz_yz[i] * a_exp + 4.0 * g_yzz_z_zz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_zz_zz[i] = -2.0 * g_y_z_zz_zz[i] * a_exp + 4.0 * g_yzz_z_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1620-1626)

    #pragma omp simd aligned(g_x_x_xx_xx, g_x_x_xx_xy, g_x_x_xx_xz, g_x_x_xx_yy, g_x_x_xx_yz, g_x_x_xx_zz, g_xzz_x_xx_xx, g_xzz_x_xx_xy, g_xzz_x_xx_xz, g_xzz_x_xx_yy, g_xzz_x_xx_yz, g_xzz_x_xx_zz, g_zz_0_0_0_x_x_xx_xx, g_zz_0_0_0_x_x_xx_xy, g_zz_0_0_0_x_x_xx_xz, g_zz_0_0_0_x_x_xx_yy, g_zz_0_0_0_x_x_xx_yz, g_zz_0_0_0_x_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_x_xx_xx[i] = -2.0 * g_x_x_xx_xx[i] * a_exp + 4.0 * g_xzz_x_xx_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_xx_xy[i] = -2.0 * g_x_x_xx_xy[i] * a_exp + 4.0 * g_xzz_x_xx_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_xx_xz[i] = -2.0 * g_x_x_xx_xz[i] * a_exp + 4.0 * g_xzz_x_xx_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_xx_yy[i] = -2.0 * g_x_x_xx_yy[i] * a_exp + 4.0 * g_xzz_x_xx_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_xx_yz[i] = -2.0 * g_x_x_xx_yz[i] * a_exp + 4.0 * g_xzz_x_xx_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_xx_zz[i] = -2.0 * g_x_x_xx_zz[i] * a_exp + 4.0 * g_xzz_x_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (1626-1632)

    #pragma omp simd aligned(g_x_x_xy_xx, g_x_x_xy_xy, g_x_x_xy_xz, g_x_x_xy_yy, g_x_x_xy_yz, g_x_x_xy_zz, g_xzz_x_xy_xx, g_xzz_x_xy_xy, g_xzz_x_xy_xz, g_xzz_x_xy_yy, g_xzz_x_xy_yz, g_xzz_x_xy_zz, g_zz_0_0_0_x_x_xy_xx, g_zz_0_0_0_x_x_xy_xy, g_zz_0_0_0_x_x_xy_xz, g_zz_0_0_0_x_x_xy_yy, g_zz_0_0_0_x_x_xy_yz, g_zz_0_0_0_x_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_x_xy_xx[i] = -2.0 * g_x_x_xy_xx[i] * a_exp + 4.0 * g_xzz_x_xy_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_xy_xy[i] = -2.0 * g_x_x_xy_xy[i] * a_exp + 4.0 * g_xzz_x_xy_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_xy_xz[i] = -2.0 * g_x_x_xy_xz[i] * a_exp + 4.0 * g_xzz_x_xy_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_xy_yy[i] = -2.0 * g_x_x_xy_yy[i] * a_exp + 4.0 * g_xzz_x_xy_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_xy_yz[i] = -2.0 * g_x_x_xy_yz[i] * a_exp + 4.0 * g_xzz_x_xy_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_xy_zz[i] = -2.0 * g_x_x_xy_zz[i] * a_exp + 4.0 * g_xzz_x_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1632-1638)

    #pragma omp simd aligned(g_x_x_xz_xx, g_x_x_xz_xy, g_x_x_xz_xz, g_x_x_xz_yy, g_x_x_xz_yz, g_x_x_xz_zz, g_xzz_x_xz_xx, g_xzz_x_xz_xy, g_xzz_x_xz_xz, g_xzz_x_xz_yy, g_xzz_x_xz_yz, g_xzz_x_xz_zz, g_zz_0_0_0_x_x_xz_xx, g_zz_0_0_0_x_x_xz_xy, g_zz_0_0_0_x_x_xz_xz, g_zz_0_0_0_x_x_xz_yy, g_zz_0_0_0_x_x_xz_yz, g_zz_0_0_0_x_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_x_xz_xx[i] = -2.0 * g_x_x_xz_xx[i] * a_exp + 4.0 * g_xzz_x_xz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_xz_xy[i] = -2.0 * g_x_x_xz_xy[i] * a_exp + 4.0 * g_xzz_x_xz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_xz_xz[i] = -2.0 * g_x_x_xz_xz[i] * a_exp + 4.0 * g_xzz_x_xz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_xz_yy[i] = -2.0 * g_x_x_xz_yy[i] * a_exp + 4.0 * g_xzz_x_xz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_xz_yz[i] = -2.0 * g_x_x_xz_yz[i] * a_exp + 4.0 * g_xzz_x_xz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_xz_zz[i] = -2.0 * g_x_x_xz_zz[i] * a_exp + 4.0 * g_xzz_x_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1638-1644)

    #pragma omp simd aligned(g_x_x_yy_xx, g_x_x_yy_xy, g_x_x_yy_xz, g_x_x_yy_yy, g_x_x_yy_yz, g_x_x_yy_zz, g_xzz_x_yy_xx, g_xzz_x_yy_xy, g_xzz_x_yy_xz, g_xzz_x_yy_yy, g_xzz_x_yy_yz, g_xzz_x_yy_zz, g_zz_0_0_0_x_x_yy_xx, g_zz_0_0_0_x_x_yy_xy, g_zz_0_0_0_x_x_yy_xz, g_zz_0_0_0_x_x_yy_yy, g_zz_0_0_0_x_x_yy_yz, g_zz_0_0_0_x_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_x_yy_xx[i] = -2.0 * g_x_x_yy_xx[i] * a_exp + 4.0 * g_xzz_x_yy_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_yy_xy[i] = -2.0 * g_x_x_yy_xy[i] * a_exp + 4.0 * g_xzz_x_yy_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_yy_xz[i] = -2.0 * g_x_x_yy_xz[i] * a_exp + 4.0 * g_xzz_x_yy_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_yy_yy[i] = -2.0 * g_x_x_yy_yy[i] * a_exp + 4.0 * g_xzz_x_yy_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_yy_yz[i] = -2.0 * g_x_x_yy_yz[i] * a_exp + 4.0 * g_xzz_x_yy_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_yy_zz[i] = -2.0 * g_x_x_yy_zz[i] * a_exp + 4.0 * g_xzz_x_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1644-1650)

    #pragma omp simd aligned(g_x_x_yz_xx, g_x_x_yz_xy, g_x_x_yz_xz, g_x_x_yz_yy, g_x_x_yz_yz, g_x_x_yz_zz, g_xzz_x_yz_xx, g_xzz_x_yz_xy, g_xzz_x_yz_xz, g_xzz_x_yz_yy, g_xzz_x_yz_yz, g_xzz_x_yz_zz, g_zz_0_0_0_x_x_yz_xx, g_zz_0_0_0_x_x_yz_xy, g_zz_0_0_0_x_x_yz_xz, g_zz_0_0_0_x_x_yz_yy, g_zz_0_0_0_x_x_yz_yz, g_zz_0_0_0_x_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_x_yz_xx[i] = -2.0 * g_x_x_yz_xx[i] * a_exp + 4.0 * g_xzz_x_yz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_yz_xy[i] = -2.0 * g_x_x_yz_xy[i] * a_exp + 4.0 * g_xzz_x_yz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_yz_xz[i] = -2.0 * g_x_x_yz_xz[i] * a_exp + 4.0 * g_xzz_x_yz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_yz_yy[i] = -2.0 * g_x_x_yz_yy[i] * a_exp + 4.0 * g_xzz_x_yz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_yz_yz[i] = -2.0 * g_x_x_yz_yz[i] * a_exp + 4.0 * g_xzz_x_yz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_yz_zz[i] = -2.0 * g_x_x_yz_zz[i] * a_exp + 4.0 * g_xzz_x_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1650-1656)

    #pragma omp simd aligned(g_x_x_zz_xx, g_x_x_zz_xy, g_x_x_zz_xz, g_x_x_zz_yy, g_x_x_zz_yz, g_x_x_zz_zz, g_xzz_x_zz_xx, g_xzz_x_zz_xy, g_xzz_x_zz_xz, g_xzz_x_zz_yy, g_xzz_x_zz_yz, g_xzz_x_zz_zz, g_zz_0_0_0_x_x_zz_xx, g_zz_0_0_0_x_x_zz_xy, g_zz_0_0_0_x_x_zz_xz, g_zz_0_0_0_x_x_zz_yy, g_zz_0_0_0_x_x_zz_yz, g_zz_0_0_0_x_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_x_zz_xx[i] = -2.0 * g_x_x_zz_xx[i] * a_exp + 4.0 * g_xzz_x_zz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_zz_xy[i] = -2.0 * g_x_x_zz_xy[i] * a_exp + 4.0 * g_xzz_x_zz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_zz_xz[i] = -2.0 * g_x_x_zz_xz[i] * a_exp + 4.0 * g_xzz_x_zz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_zz_yy[i] = -2.0 * g_x_x_zz_yy[i] * a_exp + 4.0 * g_xzz_x_zz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_zz_yz[i] = -2.0 * g_x_x_zz_yz[i] * a_exp + 4.0 * g_xzz_x_zz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_zz_zz[i] = -2.0 * g_x_x_zz_zz[i] * a_exp + 4.0 * g_xzz_x_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1656-1662)

    #pragma omp simd aligned(g_x_y_xx_xx, g_x_y_xx_xy, g_x_y_xx_xz, g_x_y_xx_yy, g_x_y_xx_yz, g_x_y_xx_zz, g_xzz_y_xx_xx, g_xzz_y_xx_xy, g_xzz_y_xx_xz, g_xzz_y_xx_yy, g_xzz_y_xx_yz, g_xzz_y_xx_zz, g_zz_0_0_0_x_y_xx_xx, g_zz_0_0_0_x_y_xx_xy, g_zz_0_0_0_x_y_xx_xz, g_zz_0_0_0_x_y_xx_yy, g_zz_0_0_0_x_y_xx_yz, g_zz_0_0_0_x_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_y_xx_xx[i] = -2.0 * g_x_y_xx_xx[i] * a_exp + 4.0 * g_xzz_y_xx_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_xx_xy[i] = -2.0 * g_x_y_xx_xy[i] * a_exp + 4.0 * g_xzz_y_xx_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_xx_xz[i] = -2.0 * g_x_y_xx_xz[i] * a_exp + 4.0 * g_xzz_y_xx_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_xx_yy[i] = -2.0 * g_x_y_xx_yy[i] * a_exp + 4.0 * g_xzz_y_xx_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_xx_yz[i] = -2.0 * g_x_y_xx_yz[i] * a_exp + 4.0 * g_xzz_y_xx_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_xx_zz[i] = -2.0 * g_x_y_xx_zz[i] * a_exp + 4.0 * g_xzz_y_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (1662-1668)

    #pragma omp simd aligned(g_x_y_xy_xx, g_x_y_xy_xy, g_x_y_xy_xz, g_x_y_xy_yy, g_x_y_xy_yz, g_x_y_xy_zz, g_xzz_y_xy_xx, g_xzz_y_xy_xy, g_xzz_y_xy_xz, g_xzz_y_xy_yy, g_xzz_y_xy_yz, g_xzz_y_xy_zz, g_zz_0_0_0_x_y_xy_xx, g_zz_0_0_0_x_y_xy_xy, g_zz_0_0_0_x_y_xy_xz, g_zz_0_0_0_x_y_xy_yy, g_zz_0_0_0_x_y_xy_yz, g_zz_0_0_0_x_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_y_xy_xx[i] = -2.0 * g_x_y_xy_xx[i] * a_exp + 4.0 * g_xzz_y_xy_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_xy_xy[i] = -2.0 * g_x_y_xy_xy[i] * a_exp + 4.0 * g_xzz_y_xy_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_xy_xz[i] = -2.0 * g_x_y_xy_xz[i] * a_exp + 4.0 * g_xzz_y_xy_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_xy_yy[i] = -2.0 * g_x_y_xy_yy[i] * a_exp + 4.0 * g_xzz_y_xy_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_xy_yz[i] = -2.0 * g_x_y_xy_yz[i] * a_exp + 4.0 * g_xzz_y_xy_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_xy_zz[i] = -2.0 * g_x_y_xy_zz[i] * a_exp + 4.0 * g_xzz_y_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1668-1674)

    #pragma omp simd aligned(g_x_y_xz_xx, g_x_y_xz_xy, g_x_y_xz_xz, g_x_y_xz_yy, g_x_y_xz_yz, g_x_y_xz_zz, g_xzz_y_xz_xx, g_xzz_y_xz_xy, g_xzz_y_xz_xz, g_xzz_y_xz_yy, g_xzz_y_xz_yz, g_xzz_y_xz_zz, g_zz_0_0_0_x_y_xz_xx, g_zz_0_0_0_x_y_xz_xy, g_zz_0_0_0_x_y_xz_xz, g_zz_0_0_0_x_y_xz_yy, g_zz_0_0_0_x_y_xz_yz, g_zz_0_0_0_x_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_y_xz_xx[i] = -2.0 * g_x_y_xz_xx[i] * a_exp + 4.0 * g_xzz_y_xz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_xz_xy[i] = -2.0 * g_x_y_xz_xy[i] * a_exp + 4.0 * g_xzz_y_xz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_xz_xz[i] = -2.0 * g_x_y_xz_xz[i] * a_exp + 4.0 * g_xzz_y_xz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_xz_yy[i] = -2.0 * g_x_y_xz_yy[i] * a_exp + 4.0 * g_xzz_y_xz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_xz_yz[i] = -2.0 * g_x_y_xz_yz[i] * a_exp + 4.0 * g_xzz_y_xz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_xz_zz[i] = -2.0 * g_x_y_xz_zz[i] * a_exp + 4.0 * g_xzz_y_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1674-1680)

    #pragma omp simd aligned(g_x_y_yy_xx, g_x_y_yy_xy, g_x_y_yy_xz, g_x_y_yy_yy, g_x_y_yy_yz, g_x_y_yy_zz, g_xzz_y_yy_xx, g_xzz_y_yy_xy, g_xzz_y_yy_xz, g_xzz_y_yy_yy, g_xzz_y_yy_yz, g_xzz_y_yy_zz, g_zz_0_0_0_x_y_yy_xx, g_zz_0_0_0_x_y_yy_xy, g_zz_0_0_0_x_y_yy_xz, g_zz_0_0_0_x_y_yy_yy, g_zz_0_0_0_x_y_yy_yz, g_zz_0_0_0_x_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_y_yy_xx[i] = -2.0 * g_x_y_yy_xx[i] * a_exp + 4.0 * g_xzz_y_yy_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_yy_xy[i] = -2.0 * g_x_y_yy_xy[i] * a_exp + 4.0 * g_xzz_y_yy_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_yy_xz[i] = -2.0 * g_x_y_yy_xz[i] * a_exp + 4.0 * g_xzz_y_yy_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_yy_yy[i] = -2.0 * g_x_y_yy_yy[i] * a_exp + 4.0 * g_xzz_y_yy_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_yy_yz[i] = -2.0 * g_x_y_yy_yz[i] * a_exp + 4.0 * g_xzz_y_yy_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_yy_zz[i] = -2.0 * g_x_y_yy_zz[i] * a_exp + 4.0 * g_xzz_y_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1680-1686)

    #pragma omp simd aligned(g_x_y_yz_xx, g_x_y_yz_xy, g_x_y_yz_xz, g_x_y_yz_yy, g_x_y_yz_yz, g_x_y_yz_zz, g_xzz_y_yz_xx, g_xzz_y_yz_xy, g_xzz_y_yz_xz, g_xzz_y_yz_yy, g_xzz_y_yz_yz, g_xzz_y_yz_zz, g_zz_0_0_0_x_y_yz_xx, g_zz_0_0_0_x_y_yz_xy, g_zz_0_0_0_x_y_yz_xz, g_zz_0_0_0_x_y_yz_yy, g_zz_0_0_0_x_y_yz_yz, g_zz_0_0_0_x_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_y_yz_xx[i] = -2.0 * g_x_y_yz_xx[i] * a_exp + 4.0 * g_xzz_y_yz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_yz_xy[i] = -2.0 * g_x_y_yz_xy[i] * a_exp + 4.0 * g_xzz_y_yz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_yz_xz[i] = -2.0 * g_x_y_yz_xz[i] * a_exp + 4.0 * g_xzz_y_yz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_yz_yy[i] = -2.0 * g_x_y_yz_yy[i] * a_exp + 4.0 * g_xzz_y_yz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_yz_yz[i] = -2.0 * g_x_y_yz_yz[i] * a_exp + 4.0 * g_xzz_y_yz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_yz_zz[i] = -2.0 * g_x_y_yz_zz[i] * a_exp + 4.0 * g_xzz_y_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1686-1692)

    #pragma omp simd aligned(g_x_y_zz_xx, g_x_y_zz_xy, g_x_y_zz_xz, g_x_y_zz_yy, g_x_y_zz_yz, g_x_y_zz_zz, g_xzz_y_zz_xx, g_xzz_y_zz_xy, g_xzz_y_zz_xz, g_xzz_y_zz_yy, g_xzz_y_zz_yz, g_xzz_y_zz_zz, g_zz_0_0_0_x_y_zz_xx, g_zz_0_0_0_x_y_zz_xy, g_zz_0_0_0_x_y_zz_xz, g_zz_0_0_0_x_y_zz_yy, g_zz_0_0_0_x_y_zz_yz, g_zz_0_0_0_x_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_y_zz_xx[i] = -2.0 * g_x_y_zz_xx[i] * a_exp + 4.0 * g_xzz_y_zz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_zz_xy[i] = -2.0 * g_x_y_zz_xy[i] * a_exp + 4.0 * g_xzz_y_zz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_zz_xz[i] = -2.0 * g_x_y_zz_xz[i] * a_exp + 4.0 * g_xzz_y_zz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_zz_yy[i] = -2.0 * g_x_y_zz_yy[i] * a_exp + 4.0 * g_xzz_y_zz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_zz_yz[i] = -2.0 * g_x_y_zz_yz[i] * a_exp + 4.0 * g_xzz_y_zz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_zz_zz[i] = -2.0 * g_x_y_zz_zz[i] * a_exp + 4.0 * g_xzz_y_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1692-1698)

    #pragma omp simd aligned(g_x_z_xx_xx, g_x_z_xx_xy, g_x_z_xx_xz, g_x_z_xx_yy, g_x_z_xx_yz, g_x_z_xx_zz, g_xzz_z_xx_xx, g_xzz_z_xx_xy, g_xzz_z_xx_xz, g_xzz_z_xx_yy, g_xzz_z_xx_yz, g_xzz_z_xx_zz, g_zz_0_0_0_x_z_xx_xx, g_zz_0_0_0_x_z_xx_xy, g_zz_0_0_0_x_z_xx_xz, g_zz_0_0_0_x_z_xx_yy, g_zz_0_0_0_x_z_xx_yz, g_zz_0_0_0_x_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_z_xx_xx[i] = -2.0 * g_x_z_xx_xx[i] * a_exp + 4.0 * g_xzz_z_xx_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_xx_xy[i] = -2.0 * g_x_z_xx_xy[i] * a_exp + 4.0 * g_xzz_z_xx_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_xx_xz[i] = -2.0 * g_x_z_xx_xz[i] * a_exp + 4.0 * g_xzz_z_xx_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_xx_yy[i] = -2.0 * g_x_z_xx_yy[i] * a_exp + 4.0 * g_xzz_z_xx_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_xx_yz[i] = -2.0 * g_x_z_xx_yz[i] * a_exp + 4.0 * g_xzz_z_xx_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_xx_zz[i] = -2.0 * g_x_z_xx_zz[i] * a_exp + 4.0 * g_xzz_z_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (1698-1704)

    #pragma omp simd aligned(g_x_z_xy_xx, g_x_z_xy_xy, g_x_z_xy_xz, g_x_z_xy_yy, g_x_z_xy_yz, g_x_z_xy_zz, g_xzz_z_xy_xx, g_xzz_z_xy_xy, g_xzz_z_xy_xz, g_xzz_z_xy_yy, g_xzz_z_xy_yz, g_xzz_z_xy_zz, g_zz_0_0_0_x_z_xy_xx, g_zz_0_0_0_x_z_xy_xy, g_zz_0_0_0_x_z_xy_xz, g_zz_0_0_0_x_z_xy_yy, g_zz_0_0_0_x_z_xy_yz, g_zz_0_0_0_x_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_z_xy_xx[i] = -2.0 * g_x_z_xy_xx[i] * a_exp + 4.0 * g_xzz_z_xy_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_xy_xy[i] = -2.0 * g_x_z_xy_xy[i] * a_exp + 4.0 * g_xzz_z_xy_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_xy_xz[i] = -2.0 * g_x_z_xy_xz[i] * a_exp + 4.0 * g_xzz_z_xy_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_xy_yy[i] = -2.0 * g_x_z_xy_yy[i] * a_exp + 4.0 * g_xzz_z_xy_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_xy_yz[i] = -2.0 * g_x_z_xy_yz[i] * a_exp + 4.0 * g_xzz_z_xy_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_xy_zz[i] = -2.0 * g_x_z_xy_zz[i] * a_exp + 4.0 * g_xzz_z_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1704-1710)

    #pragma omp simd aligned(g_x_z_xz_xx, g_x_z_xz_xy, g_x_z_xz_xz, g_x_z_xz_yy, g_x_z_xz_yz, g_x_z_xz_zz, g_xzz_z_xz_xx, g_xzz_z_xz_xy, g_xzz_z_xz_xz, g_xzz_z_xz_yy, g_xzz_z_xz_yz, g_xzz_z_xz_zz, g_zz_0_0_0_x_z_xz_xx, g_zz_0_0_0_x_z_xz_xy, g_zz_0_0_0_x_z_xz_xz, g_zz_0_0_0_x_z_xz_yy, g_zz_0_0_0_x_z_xz_yz, g_zz_0_0_0_x_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_z_xz_xx[i] = -2.0 * g_x_z_xz_xx[i] * a_exp + 4.0 * g_xzz_z_xz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_xz_xy[i] = -2.0 * g_x_z_xz_xy[i] * a_exp + 4.0 * g_xzz_z_xz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_xz_xz[i] = -2.0 * g_x_z_xz_xz[i] * a_exp + 4.0 * g_xzz_z_xz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_xz_yy[i] = -2.0 * g_x_z_xz_yy[i] * a_exp + 4.0 * g_xzz_z_xz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_xz_yz[i] = -2.0 * g_x_z_xz_yz[i] * a_exp + 4.0 * g_xzz_z_xz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_xz_zz[i] = -2.0 * g_x_z_xz_zz[i] * a_exp + 4.0 * g_xzz_z_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1710-1716)

    #pragma omp simd aligned(g_x_z_yy_xx, g_x_z_yy_xy, g_x_z_yy_xz, g_x_z_yy_yy, g_x_z_yy_yz, g_x_z_yy_zz, g_xzz_z_yy_xx, g_xzz_z_yy_xy, g_xzz_z_yy_xz, g_xzz_z_yy_yy, g_xzz_z_yy_yz, g_xzz_z_yy_zz, g_zz_0_0_0_x_z_yy_xx, g_zz_0_0_0_x_z_yy_xy, g_zz_0_0_0_x_z_yy_xz, g_zz_0_0_0_x_z_yy_yy, g_zz_0_0_0_x_z_yy_yz, g_zz_0_0_0_x_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_z_yy_xx[i] = -2.0 * g_x_z_yy_xx[i] * a_exp + 4.0 * g_xzz_z_yy_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_yy_xy[i] = -2.0 * g_x_z_yy_xy[i] * a_exp + 4.0 * g_xzz_z_yy_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_yy_xz[i] = -2.0 * g_x_z_yy_xz[i] * a_exp + 4.0 * g_xzz_z_yy_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_yy_yy[i] = -2.0 * g_x_z_yy_yy[i] * a_exp + 4.0 * g_xzz_z_yy_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_yy_yz[i] = -2.0 * g_x_z_yy_yz[i] * a_exp + 4.0 * g_xzz_z_yy_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_yy_zz[i] = -2.0 * g_x_z_yy_zz[i] * a_exp + 4.0 * g_xzz_z_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1716-1722)

    #pragma omp simd aligned(g_x_z_yz_xx, g_x_z_yz_xy, g_x_z_yz_xz, g_x_z_yz_yy, g_x_z_yz_yz, g_x_z_yz_zz, g_xzz_z_yz_xx, g_xzz_z_yz_xy, g_xzz_z_yz_xz, g_xzz_z_yz_yy, g_xzz_z_yz_yz, g_xzz_z_yz_zz, g_zz_0_0_0_x_z_yz_xx, g_zz_0_0_0_x_z_yz_xy, g_zz_0_0_0_x_z_yz_xz, g_zz_0_0_0_x_z_yz_yy, g_zz_0_0_0_x_z_yz_yz, g_zz_0_0_0_x_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_z_yz_xx[i] = -2.0 * g_x_z_yz_xx[i] * a_exp + 4.0 * g_xzz_z_yz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_yz_xy[i] = -2.0 * g_x_z_yz_xy[i] * a_exp + 4.0 * g_xzz_z_yz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_yz_xz[i] = -2.0 * g_x_z_yz_xz[i] * a_exp + 4.0 * g_xzz_z_yz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_yz_yy[i] = -2.0 * g_x_z_yz_yy[i] * a_exp + 4.0 * g_xzz_z_yz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_yz_yz[i] = -2.0 * g_x_z_yz_yz[i] * a_exp + 4.0 * g_xzz_z_yz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_yz_zz[i] = -2.0 * g_x_z_yz_zz[i] * a_exp + 4.0 * g_xzz_z_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1722-1728)

    #pragma omp simd aligned(g_x_z_zz_xx, g_x_z_zz_xy, g_x_z_zz_xz, g_x_z_zz_yy, g_x_z_zz_yz, g_x_z_zz_zz, g_xzz_z_zz_xx, g_xzz_z_zz_xy, g_xzz_z_zz_xz, g_xzz_z_zz_yy, g_xzz_z_zz_yz, g_xzz_z_zz_zz, g_zz_0_0_0_x_z_zz_xx, g_zz_0_0_0_x_z_zz_xy, g_zz_0_0_0_x_z_zz_xz, g_zz_0_0_0_x_z_zz_yy, g_zz_0_0_0_x_z_zz_yz, g_zz_0_0_0_x_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_z_zz_xx[i] = -2.0 * g_x_z_zz_xx[i] * a_exp + 4.0 * g_xzz_z_zz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_zz_xy[i] = -2.0 * g_x_z_zz_xy[i] * a_exp + 4.0 * g_xzz_z_zz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_zz_xz[i] = -2.0 * g_x_z_zz_xz[i] * a_exp + 4.0 * g_xzz_z_zz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_zz_yy[i] = -2.0 * g_x_z_zz_yy[i] * a_exp + 4.0 * g_xzz_z_zz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_zz_yz[i] = -2.0 * g_x_z_zz_yz[i] * a_exp + 4.0 * g_xzz_z_zz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_zz_zz[i] = -2.0 * g_x_z_zz_zz[i] * a_exp + 4.0 * g_xzz_z_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1728-1734)

    #pragma omp simd aligned(g_y_x_xx_xx, g_y_x_xx_xy, g_y_x_xx_xz, g_y_x_xx_yy, g_y_x_xx_yz, g_y_x_xx_zz, g_yzz_x_xx_xx, g_yzz_x_xx_xy, g_yzz_x_xx_xz, g_yzz_x_xx_yy, g_yzz_x_xx_yz, g_yzz_x_xx_zz, g_zz_0_0_0_y_x_xx_xx, g_zz_0_0_0_y_x_xx_xy, g_zz_0_0_0_y_x_xx_xz, g_zz_0_0_0_y_x_xx_yy, g_zz_0_0_0_y_x_xx_yz, g_zz_0_0_0_y_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_x_xx_xx[i] = -2.0 * g_y_x_xx_xx[i] * a_exp + 4.0 * g_yzz_x_xx_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_xx_xy[i] = -2.0 * g_y_x_xx_xy[i] * a_exp + 4.0 * g_yzz_x_xx_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_xx_xz[i] = -2.0 * g_y_x_xx_xz[i] * a_exp + 4.0 * g_yzz_x_xx_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_xx_yy[i] = -2.0 * g_y_x_xx_yy[i] * a_exp + 4.0 * g_yzz_x_xx_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_xx_yz[i] = -2.0 * g_y_x_xx_yz[i] * a_exp + 4.0 * g_yzz_x_xx_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_xx_zz[i] = -2.0 * g_y_x_xx_zz[i] * a_exp + 4.0 * g_yzz_x_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (1734-1740)

    #pragma omp simd aligned(g_y_x_xy_xx, g_y_x_xy_xy, g_y_x_xy_xz, g_y_x_xy_yy, g_y_x_xy_yz, g_y_x_xy_zz, g_yzz_x_xy_xx, g_yzz_x_xy_xy, g_yzz_x_xy_xz, g_yzz_x_xy_yy, g_yzz_x_xy_yz, g_yzz_x_xy_zz, g_zz_0_0_0_y_x_xy_xx, g_zz_0_0_0_y_x_xy_xy, g_zz_0_0_0_y_x_xy_xz, g_zz_0_0_0_y_x_xy_yy, g_zz_0_0_0_y_x_xy_yz, g_zz_0_0_0_y_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_x_xy_xx[i] = -2.0 * g_y_x_xy_xx[i] * a_exp + 4.0 * g_yzz_x_xy_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_xy_xy[i] = -2.0 * g_y_x_xy_xy[i] * a_exp + 4.0 * g_yzz_x_xy_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_xy_xz[i] = -2.0 * g_y_x_xy_xz[i] * a_exp + 4.0 * g_yzz_x_xy_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_xy_yy[i] = -2.0 * g_y_x_xy_yy[i] * a_exp + 4.0 * g_yzz_x_xy_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_xy_yz[i] = -2.0 * g_y_x_xy_yz[i] * a_exp + 4.0 * g_yzz_x_xy_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_xy_zz[i] = -2.0 * g_y_x_xy_zz[i] * a_exp + 4.0 * g_yzz_x_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1740-1746)

    #pragma omp simd aligned(g_y_x_xz_xx, g_y_x_xz_xy, g_y_x_xz_xz, g_y_x_xz_yy, g_y_x_xz_yz, g_y_x_xz_zz, g_yzz_x_xz_xx, g_yzz_x_xz_xy, g_yzz_x_xz_xz, g_yzz_x_xz_yy, g_yzz_x_xz_yz, g_yzz_x_xz_zz, g_zz_0_0_0_y_x_xz_xx, g_zz_0_0_0_y_x_xz_xy, g_zz_0_0_0_y_x_xz_xz, g_zz_0_0_0_y_x_xz_yy, g_zz_0_0_0_y_x_xz_yz, g_zz_0_0_0_y_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_x_xz_xx[i] = -2.0 * g_y_x_xz_xx[i] * a_exp + 4.0 * g_yzz_x_xz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_xz_xy[i] = -2.0 * g_y_x_xz_xy[i] * a_exp + 4.0 * g_yzz_x_xz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_xz_xz[i] = -2.0 * g_y_x_xz_xz[i] * a_exp + 4.0 * g_yzz_x_xz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_xz_yy[i] = -2.0 * g_y_x_xz_yy[i] * a_exp + 4.0 * g_yzz_x_xz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_xz_yz[i] = -2.0 * g_y_x_xz_yz[i] * a_exp + 4.0 * g_yzz_x_xz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_xz_zz[i] = -2.0 * g_y_x_xz_zz[i] * a_exp + 4.0 * g_yzz_x_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1746-1752)

    #pragma omp simd aligned(g_y_x_yy_xx, g_y_x_yy_xy, g_y_x_yy_xz, g_y_x_yy_yy, g_y_x_yy_yz, g_y_x_yy_zz, g_yzz_x_yy_xx, g_yzz_x_yy_xy, g_yzz_x_yy_xz, g_yzz_x_yy_yy, g_yzz_x_yy_yz, g_yzz_x_yy_zz, g_zz_0_0_0_y_x_yy_xx, g_zz_0_0_0_y_x_yy_xy, g_zz_0_0_0_y_x_yy_xz, g_zz_0_0_0_y_x_yy_yy, g_zz_0_0_0_y_x_yy_yz, g_zz_0_0_0_y_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_x_yy_xx[i] = -2.0 * g_y_x_yy_xx[i] * a_exp + 4.0 * g_yzz_x_yy_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_yy_xy[i] = -2.0 * g_y_x_yy_xy[i] * a_exp + 4.0 * g_yzz_x_yy_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_yy_xz[i] = -2.0 * g_y_x_yy_xz[i] * a_exp + 4.0 * g_yzz_x_yy_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_yy_yy[i] = -2.0 * g_y_x_yy_yy[i] * a_exp + 4.0 * g_yzz_x_yy_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_yy_yz[i] = -2.0 * g_y_x_yy_yz[i] * a_exp + 4.0 * g_yzz_x_yy_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_yy_zz[i] = -2.0 * g_y_x_yy_zz[i] * a_exp + 4.0 * g_yzz_x_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1752-1758)

    #pragma omp simd aligned(g_y_x_yz_xx, g_y_x_yz_xy, g_y_x_yz_xz, g_y_x_yz_yy, g_y_x_yz_yz, g_y_x_yz_zz, g_yzz_x_yz_xx, g_yzz_x_yz_xy, g_yzz_x_yz_xz, g_yzz_x_yz_yy, g_yzz_x_yz_yz, g_yzz_x_yz_zz, g_zz_0_0_0_y_x_yz_xx, g_zz_0_0_0_y_x_yz_xy, g_zz_0_0_0_y_x_yz_xz, g_zz_0_0_0_y_x_yz_yy, g_zz_0_0_0_y_x_yz_yz, g_zz_0_0_0_y_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_x_yz_xx[i] = -2.0 * g_y_x_yz_xx[i] * a_exp + 4.0 * g_yzz_x_yz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_yz_xy[i] = -2.0 * g_y_x_yz_xy[i] * a_exp + 4.0 * g_yzz_x_yz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_yz_xz[i] = -2.0 * g_y_x_yz_xz[i] * a_exp + 4.0 * g_yzz_x_yz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_yz_yy[i] = -2.0 * g_y_x_yz_yy[i] * a_exp + 4.0 * g_yzz_x_yz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_yz_yz[i] = -2.0 * g_y_x_yz_yz[i] * a_exp + 4.0 * g_yzz_x_yz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_yz_zz[i] = -2.0 * g_y_x_yz_zz[i] * a_exp + 4.0 * g_yzz_x_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1758-1764)

    #pragma omp simd aligned(g_y_x_zz_xx, g_y_x_zz_xy, g_y_x_zz_xz, g_y_x_zz_yy, g_y_x_zz_yz, g_y_x_zz_zz, g_yzz_x_zz_xx, g_yzz_x_zz_xy, g_yzz_x_zz_xz, g_yzz_x_zz_yy, g_yzz_x_zz_yz, g_yzz_x_zz_zz, g_zz_0_0_0_y_x_zz_xx, g_zz_0_0_0_y_x_zz_xy, g_zz_0_0_0_y_x_zz_xz, g_zz_0_0_0_y_x_zz_yy, g_zz_0_0_0_y_x_zz_yz, g_zz_0_0_0_y_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_x_zz_xx[i] = -2.0 * g_y_x_zz_xx[i] * a_exp + 4.0 * g_yzz_x_zz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_zz_xy[i] = -2.0 * g_y_x_zz_xy[i] * a_exp + 4.0 * g_yzz_x_zz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_zz_xz[i] = -2.0 * g_y_x_zz_xz[i] * a_exp + 4.0 * g_yzz_x_zz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_zz_yy[i] = -2.0 * g_y_x_zz_yy[i] * a_exp + 4.0 * g_yzz_x_zz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_zz_yz[i] = -2.0 * g_y_x_zz_yz[i] * a_exp + 4.0 * g_yzz_x_zz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_zz_zz[i] = -2.0 * g_y_x_zz_zz[i] * a_exp + 4.0 * g_yzz_x_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1764-1770)

    #pragma omp simd aligned(g_y_y_xx_xx, g_y_y_xx_xy, g_y_y_xx_xz, g_y_y_xx_yy, g_y_y_xx_yz, g_y_y_xx_zz, g_yzz_y_xx_xx, g_yzz_y_xx_xy, g_yzz_y_xx_xz, g_yzz_y_xx_yy, g_yzz_y_xx_yz, g_yzz_y_xx_zz, g_zz_0_0_0_y_y_xx_xx, g_zz_0_0_0_y_y_xx_xy, g_zz_0_0_0_y_y_xx_xz, g_zz_0_0_0_y_y_xx_yy, g_zz_0_0_0_y_y_xx_yz, g_zz_0_0_0_y_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_y_xx_xx[i] = -2.0 * g_y_y_xx_xx[i] * a_exp + 4.0 * g_yzz_y_xx_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_xx_xy[i] = -2.0 * g_y_y_xx_xy[i] * a_exp + 4.0 * g_yzz_y_xx_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_xx_xz[i] = -2.0 * g_y_y_xx_xz[i] * a_exp + 4.0 * g_yzz_y_xx_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_xx_yy[i] = -2.0 * g_y_y_xx_yy[i] * a_exp + 4.0 * g_yzz_y_xx_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_xx_yz[i] = -2.0 * g_y_y_xx_yz[i] * a_exp + 4.0 * g_yzz_y_xx_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_xx_zz[i] = -2.0 * g_y_y_xx_zz[i] * a_exp + 4.0 * g_yzz_y_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (1770-1776)

    #pragma omp simd aligned(g_y_y_xy_xx, g_y_y_xy_xy, g_y_y_xy_xz, g_y_y_xy_yy, g_y_y_xy_yz, g_y_y_xy_zz, g_yzz_y_xy_xx, g_yzz_y_xy_xy, g_yzz_y_xy_xz, g_yzz_y_xy_yy, g_yzz_y_xy_yz, g_yzz_y_xy_zz, g_zz_0_0_0_y_y_xy_xx, g_zz_0_0_0_y_y_xy_xy, g_zz_0_0_0_y_y_xy_xz, g_zz_0_0_0_y_y_xy_yy, g_zz_0_0_0_y_y_xy_yz, g_zz_0_0_0_y_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_y_xy_xx[i] = -2.0 * g_y_y_xy_xx[i] * a_exp + 4.0 * g_yzz_y_xy_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_xy_xy[i] = -2.0 * g_y_y_xy_xy[i] * a_exp + 4.0 * g_yzz_y_xy_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_xy_xz[i] = -2.0 * g_y_y_xy_xz[i] * a_exp + 4.0 * g_yzz_y_xy_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_xy_yy[i] = -2.0 * g_y_y_xy_yy[i] * a_exp + 4.0 * g_yzz_y_xy_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_xy_yz[i] = -2.0 * g_y_y_xy_yz[i] * a_exp + 4.0 * g_yzz_y_xy_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_xy_zz[i] = -2.0 * g_y_y_xy_zz[i] * a_exp + 4.0 * g_yzz_y_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1776-1782)

    #pragma omp simd aligned(g_y_y_xz_xx, g_y_y_xz_xy, g_y_y_xz_xz, g_y_y_xz_yy, g_y_y_xz_yz, g_y_y_xz_zz, g_yzz_y_xz_xx, g_yzz_y_xz_xy, g_yzz_y_xz_xz, g_yzz_y_xz_yy, g_yzz_y_xz_yz, g_yzz_y_xz_zz, g_zz_0_0_0_y_y_xz_xx, g_zz_0_0_0_y_y_xz_xy, g_zz_0_0_0_y_y_xz_xz, g_zz_0_0_0_y_y_xz_yy, g_zz_0_0_0_y_y_xz_yz, g_zz_0_0_0_y_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_y_xz_xx[i] = -2.0 * g_y_y_xz_xx[i] * a_exp + 4.0 * g_yzz_y_xz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_xz_xy[i] = -2.0 * g_y_y_xz_xy[i] * a_exp + 4.0 * g_yzz_y_xz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_xz_xz[i] = -2.0 * g_y_y_xz_xz[i] * a_exp + 4.0 * g_yzz_y_xz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_xz_yy[i] = -2.0 * g_y_y_xz_yy[i] * a_exp + 4.0 * g_yzz_y_xz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_xz_yz[i] = -2.0 * g_y_y_xz_yz[i] * a_exp + 4.0 * g_yzz_y_xz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_xz_zz[i] = -2.0 * g_y_y_xz_zz[i] * a_exp + 4.0 * g_yzz_y_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1782-1788)

    #pragma omp simd aligned(g_y_y_yy_xx, g_y_y_yy_xy, g_y_y_yy_xz, g_y_y_yy_yy, g_y_y_yy_yz, g_y_y_yy_zz, g_yzz_y_yy_xx, g_yzz_y_yy_xy, g_yzz_y_yy_xz, g_yzz_y_yy_yy, g_yzz_y_yy_yz, g_yzz_y_yy_zz, g_zz_0_0_0_y_y_yy_xx, g_zz_0_0_0_y_y_yy_xy, g_zz_0_0_0_y_y_yy_xz, g_zz_0_0_0_y_y_yy_yy, g_zz_0_0_0_y_y_yy_yz, g_zz_0_0_0_y_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_y_yy_xx[i] = -2.0 * g_y_y_yy_xx[i] * a_exp + 4.0 * g_yzz_y_yy_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_yy_xy[i] = -2.0 * g_y_y_yy_xy[i] * a_exp + 4.0 * g_yzz_y_yy_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_yy_xz[i] = -2.0 * g_y_y_yy_xz[i] * a_exp + 4.0 * g_yzz_y_yy_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_yy_yy[i] = -2.0 * g_y_y_yy_yy[i] * a_exp + 4.0 * g_yzz_y_yy_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_yy_yz[i] = -2.0 * g_y_y_yy_yz[i] * a_exp + 4.0 * g_yzz_y_yy_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_yy_zz[i] = -2.0 * g_y_y_yy_zz[i] * a_exp + 4.0 * g_yzz_y_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1788-1794)

    #pragma omp simd aligned(g_y_y_yz_xx, g_y_y_yz_xy, g_y_y_yz_xz, g_y_y_yz_yy, g_y_y_yz_yz, g_y_y_yz_zz, g_yzz_y_yz_xx, g_yzz_y_yz_xy, g_yzz_y_yz_xz, g_yzz_y_yz_yy, g_yzz_y_yz_yz, g_yzz_y_yz_zz, g_zz_0_0_0_y_y_yz_xx, g_zz_0_0_0_y_y_yz_xy, g_zz_0_0_0_y_y_yz_xz, g_zz_0_0_0_y_y_yz_yy, g_zz_0_0_0_y_y_yz_yz, g_zz_0_0_0_y_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_y_yz_xx[i] = -2.0 * g_y_y_yz_xx[i] * a_exp + 4.0 * g_yzz_y_yz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_yz_xy[i] = -2.0 * g_y_y_yz_xy[i] * a_exp + 4.0 * g_yzz_y_yz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_yz_xz[i] = -2.0 * g_y_y_yz_xz[i] * a_exp + 4.0 * g_yzz_y_yz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_yz_yy[i] = -2.0 * g_y_y_yz_yy[i] * a_exp + 4.0 * g_yzz_y_yz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_yz_yz[i] = -2.0 * g_y_y_yz_yz[i] * a_exp + 4.0 * g_yzz_y_yz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_yz_zz[i] = -2.0 * g_y_y_yz_zz[i] * a_exp + 4.0 * g_yzz_y_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1794-1800)

    #pragma omp simd aligned(g_y_y_zz_xx, g_y_y_zz_xy, g_y_y_zz_xz, g_y_y_zz_yy, g_y_y_zz_yz, g_y_y_zz_zz, g_yzz_y_zz_xx, g_yzz_y_zz_xy, g_yzz_y_zz_xz, g_yzz_y_zz_yy, g_yzz_y_zz_yz, g_yzz_y_zz_zz, g_zz_0_0_0_y_y_zz_xx, g_zz_0_0_0_y_y_zz_xy, g_zz_0_0_0_y_y_zz_xz, g_zz_0_0_0_y_y_zz_yy, g_zz_0_0_0_y_y_zz_yz, g_zz_0_0_0_y_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_y_zz_xx[i] = -2.0 * g_y_y_zz_xx[i] * a_exp + 4.0 * g_yzz_y_zz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_zz_xy[i] = -2.0 * g_y_y_zz_xy[i] * a_exp + 4.0 * g_yzz_y_zz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_zz_xz[i] = -2.0 * g_y_y_zz_xz[i] * a_exp + 4.0 * g_yzz_y_zz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_zz_yy[i] = -2.0 * g_y_y_zz_yy[i] * a_exp + 4.0 * g_yzz_y_zz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_zz_yz[i] = -2.0 * g_y_y_zz_yz[i] * a_exp + 4.0 * g_yzz_y_zz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_zz_zz[i] = -2.0 * g_y_y_zz_zz[i] * a_exp + 4.0 * g_yzz_y_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1800-1806)

    #pragma omp simd aligned(g_y_z_xx_xx, g_y_z_xx_xy, g_y_z_xx_xz, g_y_z_xx_yy, g_y_z_xx_yz, g_y_z_xx_zz, g_yzz_z_xx_xx, g_yzz_z_xx_xy, g_yzz_z_xx_xz, g_yzz_z_xx_yy, g_yzz_z_xx_yz, g_yzz_z_xx_zz, g_zz_0_0_0_y_z_xx_xx, g_zz_0_0_0_y_z_xx_xy, g_zz_0_0_0_y_z_xx_xz, g_zz_0_0_0_y_z_xx_yy, g_zz_0_0_0_y_z_xx_yz, g_zz_0_0_0_y_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_z_xx_xx[i] = -2.0 * g_y_z_xx_xx[i] * a_exp + 4.0 * g_yzz_z_xx_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_xx_xy[i] = -2.0 * g_y_z_xx_xy[i] * a_exp + 4.0 * g_yzz_z_xx_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_xx_xz[i] = -2.0 * g_y_z_xx_xz[i] * a_exp + 4.0 * g_yzz_z_xx_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_xx_yy[i] = -2.0 * g_y_z_xx_yy[i] * a_exp + 4.0 * g_yzz_z_xx_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_xx_yz[i] = -2.0 * g_y_z_xx_yz[i] * a_exp + 4.0 * g_yzz_z_xx_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_xx_zz[i] = -2.0 * g_y_z_xx_zz[i] * a_exp + 4.0 * g_yzz_z_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (1806-1812)

    #pragma omp simd aligned(g_y_z_xy_xx, g_y_z_xy_xy, g_y_z_xy_xz, g_y_z_xy_yy, g_y_z_xy_yz, g_y_z_xy_zz, g_yzz_z_xy_xx, g_yzz_z_xy_xy, g_yzz_z_xy_xz, g_yzz_z_xy_yy, g_yzz_z_xy_yz, g_yzz_z_xy_zz, g_zz_0_0_0_y_z_xy_xx, g_zz_0_0_0_y_z_xy_xy, g_zz_0_0_0_y_z_xy_xz, g_zz_0_0_0_y_z_xy_yy, g_zz_0_0_0_y_z_xy_yz, g_zz_0_0_0_y_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_z_xy_xx[i] = -2.0 * g_y_z_xy_xx[i] * a_exp + 4.0 * g_yzz_z_xy_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_xy_xy[i] = -2.0 * g_y_z_xy_xy[i] * a_exp + 4.0 * g_yzz_z_xy_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_xy_xz[i] = -2.0 * g_y_z_xy_xz[i] * a_exp + 4.0 * g_yzz_z_xy_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_xy_yy[i] = -2.0 * g_y_z_xy_yy[i] * a_exp + 4.0 * g_yzz_z_xy_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_xy_yz[i] = -2.0 * g_y_z_xy_yz[i] * a_exp + 4.0 * g_yzz_z_xy_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_xy_zz[i] = -2.0 * g_y_z_xy_zz[i] * a_exp + 4.0 * g_yzz_z_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1812-1818)

    #pragma omp simd aligned(g_y_z_xz_xx, g_y_z_xz_xy, g_y_z_xz_xz, g_y_z_xz_yy, g_y_z_xz_yz, g_y_z_xz_zz, g_yzz_z_xz_xx, g_yzz_z_xz_xy, g_yzz_z_xz_xz, g_yzz_z_xz_yy, g_yzz_z_xz_yz, g_yzz_z_xz_zz, g_zz_0_0_0_y_z_xz_xx, g_zz_0_0_0_y_z_xz_xy, g_zz_0_0_0_y_z_xz_xz, g_zz_0_0_0_y_z_xz_yy, g_zz_0_0_0_y_z_xz_yz, g_zz_0_0_0_y_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_z_xz_xx[i] = -2.0 * g_y_z_xz_xx[i] * a_exp + 4.0 * g_yzz_z_xz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_xz_xy[i] = -2.0 * g_y_z_xz_xy[i] * a_exp + 4.0 * g_yzz_z_xz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_xz_xz[i] = -2.0 * g_y_z_xz_xz[i] * a_exp + 4.0 * g_yzz_z_xz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_xz_yy[i] = -2.0 * g_y_z_xz_yy[i] * a_exp + 4.0 * g_yzz_z_xz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_xz_yz[i] = -2.0 * g_y_z_xz_yz[i] * a_exp + 4.0 * g_yzz_z_xz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_xz_zz[i] = -2.0 * g_y_z_xz_zz[i] * a_exp + 4.0 * g_yzz_z_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1818-1824)

    #pragma omp simd aligned(g_y_z_yy_xx, g_y_z_yy_xy, g_y_z_yy_xz, g_y_z_yy_yy, g_y_z_yy_yz, g_y_z_yy_zz, g_yzz_z_yy_xx, g_yzz_z_yy_xy, g_yzz_z_yy_xz, g_yzz_z_yy_yy, g_yzz_z_yy_yz, g_yzz_z_yy_zz, g_zz_0_0_0_y_z_yy_xx, g_zz_0_0_0_y_z_yy_xy, g_zz_0_0_0_y_z_yy_xz, g_zz_0_0_0_y_z_yy_yy, g_zz_0_0_0_y_z_yy_yz, g_zz_0_0_0_y_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_z_yy_xx[i] = -2.0 * g_y_z_yy_xx[i] * a_exp + 4.0 * g_yzz_z_yy_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_yy_xy[i] = -2.0 * g_y_z_yy_xy[i] * a_exp + 4.0 * g_yzz_z_yy_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_yy_xz[i] = -2.0 * g_y_z_yy_xz[i] * a_exp + 4.0 * g_yzz_z_yy_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_yy_yy[i] = -2.0 * g_y_z_yy_yy[i] * a_exp + 4.0 * g_yzz_z_yy_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_yy_yz[i] = -2.0 * g_y_z_yy_yz[i] * a_exp + 4.0 * g_yzz_z_yy_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_yy_zz[i] = -2.0 * g_y_z_yy_zz[i] * a_exp + 4.0 * g_yzz_z_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1824-1830)

    #pragma omp simd aligned(g_y_z_yz_xx, g_y_z_yz_xy, g_y_z_yz_xz, g_y_z_yz_yy, g_y_z_yz_yz, g_y_z_yz_zz, g_yzz_z_yz_xx, g_yzz_z_yz_xy, g_yzz_z_yz_xz, g_yzz_z_yz_yy, g_yzz_z_yz_yz, g_yzz_z_yz_zz, g_zz_0_0_0_y_z_yz_xx, g_zz_0_0_0_y_z_yz_xy, g_zz_0_0_0_y_z_yz_xz, g_zz_0_0_0_y_z_yz_yy, g_zz_0_0_0_y_z_yz_yz, g_zz_0_0_0_y_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_z_yz_xx[i] = -2.0 * g_y_z_yz_xx[i] * a_exp + 4.0 * g_yzz_z_yz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_yz_xy[i] = -2.0 * g_y_z_yz_xy[i] * a_exp + 4.0 * g_yzz_z_yz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_yz_xz[i] = -2.0 * g_y_z_yz_xz[i] * a_exp + 4.0 * g_yzz_z_yz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_yz_yy[i] = -2.0 * g_y_z_yz_yy[i] * a_exp + 4.0 * g_yzz_z_yz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_yz_yz[i] = -2.0 * g_y_z_yz_yz[i] * a_exp + 4.0 * g_yzz_z_yz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_yz_zz[i] = -2.0 * g_y_z_yz_zz[i] * a_exp + 4.0 * g_yzz_z_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1830-1836)

    #pragma omp simd aligned(g_y_z_zz_xx, g_y_z_zz_xy, g_y_z_zz_xz, g_y_z_zz_yy, g_y_z_zz_yz, g_y_z_zz_zz, g_yzz_z_zz_xx, g_yzz_z_zz_xy, g_yzz_z_zz_xz, g_yzz_z_zz_yy, g_yzz_z_zz_yz, g_yzz_z_zz_zz, g_zz_0_0_0_y_z_zz_xx, g_zz_0_0_0_y_z_zz_xy, g_zz_0_0_0_y_z_zz_xz, g_zz_0_0_0_y_z_zz_yy, g_zz_0_0_0_y_z_zz_yz, g_zz_0_0_0_y_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_z_zz_xx[i] = -2.0 * g_y_z_zz_xx[i] * a_exp + 4.0 * g_yzz_z_zz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_zz_xy[i] = -2.0 * g_y_z_zz_xy[i] * a_exp + 4.0 * g_yzz_z_zz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_zz_xz[i] = -2.0 * g_y_z_zz_xz[i] * a_exp + 4.0 * g_yzz_z_zz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_zz_yy[i] = -2.0 * g_y_z_zz_yy[i] * a_exp + 4.0 * g_yzz_z_zz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_zz_yz[i] = -2.0 * g_y_z_zz_yz[i] * a_exp + 4.0 * g_yzz_z_zz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_zz_zz[i] = -2.0 * g_y_z_zz_zz[i] * a_exp + 4.0 * g_yzz_z_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1836-1842)

    #pragma omp simd aligned(g_z_x_xx_xx, g_z_x_xx_xy, g_z_x_xx_xz, g_z_x_xx_yy, g_z_x_xx_yz, g_z_x_xx_zz, g_zz_0_0_0_z_x_xx_xx, g_zz_0_0_0_z_x_xx_xy, g_zz_0_0_0_z_x_xx_xz, g_zz_0_0_0_z_x_xx_yy, g_zz_0_0_0_z_x_xx_yz, g_zz_0_0_0_z_x_xx_zz, g_zzz_x_xx_xx, g_zzz_x_xx_xy, g_zzz_x_xx_xz, g_zzz_x_xx_yy, g_zzz_x_xx_yz, g_zzz_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_x_xx_xx[i] = -6.0 * g_z_x_xx_xx[i] * a_exp + 4.0 * g_zzz_x_xx_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_xx_xy[i] = -6.0 * g_z_x_xx_xy[i] * a_exp + 4.0 * g_zzz_x_xx_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_xx_xz[i] = -6.0 * g_z_x_xx_xz[i] * a_exp + 4.0 * g_zzz_x_xx_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_xx_yy[i] = -6.0 * g_z_x_xx_yy[i] * a_exp + 4.0 * g_zzz_x_xx_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_xx_yz[i] = -6.0 * g_z_x_xx_yz[i] * a_exp + 4.0 * g_zzz_x_xx_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_xx_zz[i] = -6.0 * g_z_x_xx_zz[i] * a_exp + 4.0 * g_zzz_x_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (1842-1848)

    #pragma omp simd aligned(g_z_x_xy_xx, g_z_x_xy_xy, g_z_x_xy_xz, g_z_x_xy_yy, g_z_x_xy_yz, g_z_x_xy_zz, g_zz_0_0_0_z_x_xy_xx, g_zz_0_0_0_z_x_xy_xy, g_zz_0_0_0_z_x_xy_xz, g_zz_0_0_0_z_x_xy_yy, g_zz_0_0_0_z_x_xy_yz, g_zz_0_0_0_z_x_xy_zz, g_zzz_x_xy_xx, g_zzz_x_xy_xy, g_zzz_x_xy_xz, g_zzz_x_xy_yy, g_zzz_x_xy_yz, g_zzz_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_x_xy_xx[i] = -6.0 * g_z_x_xy_xx[i] * a_exp + 4.0 * g_zzz_x_xy_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_xy_xy[i] = -6.0 * g_z_x_xy_xy[i] * a_exp + 4.0 * g_zzz_x_xy_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_xy_xz[i] = -6.0 * g_z_x_xy_xz[i] * a_exp + 4.0 * g_zzz_x_xy_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_xy_yy[i] = -6.0 * g_z_x_xy_yy[i] * a_exp + 4.0 * g_zzz_x_xy_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_xy_yz[i] = -6.0 * g_z_x_xy_yz[i] * a_exp + 4.0 * g_zzz_x_xy_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_xy_zz[i] = -6.0 * g_z_x_xy_zz[i] * a_exp + 4.0 * g_zzz_x_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1848-1854)

    #pragma omp simd aligned(g_z_x_xz_xx, g_z_x_xz_xy, g_z_x_xz_xz, g_z_x_xz_yy, g_z_x_xz_yz, g_z_x_xz_zz, g_zz_0_0_0_z_x_xz_xx, g_zz_0_0_0_z_x_xz_xy, g_zz_0_0_0_z_x_xz_xz, g_zz_0_0_0_z_x_xz_yy, g_zz_0_0_0_z_x_xz_yz, g_zz_0_0_0_z_x_xz_zz, g_zzz_x_xz_xx, g_zzz_x_xz_xy, g_zzz_x_xz_xz, g_zzz_x_xz_yy, g_zzz_x_xz_yz, g_zzz_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_x_xz_xx[i] = -6.0 * g_z_x_xz_xx[i] * a_exp + 4.0 * g_zzz_x_xz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_xz_xy[i] = -6.0 * g_z_x_xz_xy[i] * a_exp + 4.0 * g_zzz_x_xz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_xz_xz[i] = -6.0 * g_z_x_xz_xz[i] * a_exp + 4.0 * g_zzz_x_xz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_xz_yy[i] = -6.0 * g_z_x_xz_yy[i] * a_exp + 4.0 * g_zzz_x_xz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_xz_yz[i] = -6.0 * g_z_x_xz_yz[i] * a_exp + 4.0 * g_zzz_x_xz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_xz_zz[i] = -6.0 * g_z_x_xz_zz[i] * a_exp + 4.0 * g_zzz_x_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1854-1860)

    #pragma omp simd aligned(g_z_x_yy_xx, g_z_x_yy_xy, g_z_x_yy_xz, g_z_x_yy_yy, g_z_x_yy_yz, g_z_x_yy_zz, g_zz_0_0_0_z_x_yy_xx, g_zz_0_0_0_z_x_yy_xy, g_zz_0_0_0_z_x_yy_xz, g_zz_0_0_0_z_x_yy_yy, g_zz_0_0_0_z_x_yy_yz, g_zz_0_0_0_z_x_yy_zz, g_zzz_x_yy_xx, g_zzz_x_yy_xy, g_zzz_x_yy_xz, g_zzz_x_yy_yy, g_zzz_x_yy_yz, g_zzz_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_x_yy_xx[i] = -6.0 * g_z_x_yy_xx[i] * a_exp + 4.0 * g_zzz_x_yy_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_yy_xy[i] = -6.0 * g_z_x_yy_xy[i] * a_exp + 4.0 * g_zzz_x_yy_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_yy_xz[i] = -6.0 * g_z_x_yy_xz[i] * a_exp + 4.0 * g_zzz_x_yy_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_yy_yy[i] = -6.0 * g_z_x_yy_yy[i] * a_exp + 4.0 * g_zzz_x_yy_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_yy_yz[i] = -6.0 * g_z_x_yy_yz[i] * a_exp + 4.0 * g_zzz_x_yy_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_yy_zz[i] = -6.0 * g_z_x_yy_zz[i] * a_exp + 4.0 * g_zzz_x_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1860-1866)

    #pragma omp simd aligned(g_z_x_yz_xx, g_z_x_yz_xy, g_z_x_yz_xz, g_z_x_yz_yy, g_z_x_yz_yz, g_z_x_yz_zz, g_zz_0_0_0_z_x_yz_xx, g_zz_0_0_0_z_x_yz_xy, g_zz_0_0_0_z_x_yz_xz, g_zz_0_0_0_z_x_yz_yy, g_zz_0_0_0_z_x_yz_yz, g_zz_0_0_0_z_x_yz_zz, g_zzz_x_yz_xx, g_zzz_x_yz_xy, g_zzz_x_yz_xz, g_zzz_x_yz_yy, g_zzz_x_yz_yz, g_zzz_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_x_yz_xx[i] = -6.0 * g_z_x_yz_xx[i] * a_exp + 4.0 * g_zzz_x_yz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_yz_xy[i] = -6.0 * g_z_x_yz_xy[i] * a_exp + 4.0 * g_zzz_x_yz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_yz_xz[i] = -6.0 * g_z_x_yz_xz[i] * a_exp + 4.0 * g_zzz_x_yz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_yz_yy[i] = -6.0 * g_z_x_yz_yy[i] * a_exp + 4.0 * g_zzz_x_yz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_yz_yz[i] = -6.0 * g_z_x_yz_yz[i] * a_exp + 4.0 * g_zzz_x_yz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_yz_zz[i] = -6.0 * g_z_x_yz_zz[i] * a_exp + 4.0 * g_zzz_x_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1866-1872)

    #pragma omp simd aligned(g_z_x_zz_xx, g_z_x_zz_xy, g_z_x_zz_xz, g_z_x_zz_yy, g_z_x_zz_yz, g_z_x_zz_zz, g_zz_0_0_0_z_x_zz_xx, g_zz_0_0_0_z_x_zz_xy, g_zz_0_0_0_z_x_zz_xz, g_zz_0_0_0_z_x_zz_yy, g_zz_0_0_0_z_x_zz_yz, g_zz_0_0_0_z_x_zz_zz, g_zzz_x_zz_xx, g_zzz_x_zz_xy, g_zzz_x_zz_xz, g_zzz_x_zz_yy, g_zzz_x_zz_yz, g_zzz_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_x_zz_xx[i] = -6.0 * g_z_x_zz_xx[i] * a_exp + 4.0 * g_zzz_x_zz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_zz_xy[i] = -6.0 * g_z_x_zz_xy[i] * a_exp + 4.0 * g_zzz_x_zz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_zz_xz[i] = -6.0 * g_z_x_zz_xz[i] * a_exp + 4.0 * g_zzz_x_zz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_zz_yy[i] = -6.0 * g_z_x_zz_yy[i] * a_exp + 4.0 * g_zzz_x_zz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_zz_yz[i] = -6.0 * g_z_x_zz_yz[i] * a_exp + 4.0 * g_zzz_x_zz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_zz_zz[i] = -6.0 * g_z_x_zz_zz[i] * a_exp + 4.0 * g_zzz_x_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1872-1878)

    #pragma omp simd aligned(g_z_y_xx_xx, g_z_y_xx_xy, g_z_y_xx_xz, g_z_y_xx_yy, g_z_y_xx_yz, g_z_y_xx_zz, g_zz_0_0_0_z_y_xx_xx, g_zz_0_0_0_z_y_xx_xy, g_zz_0_0_0_z_y_xx_xz, g_zz_0_0_0_z_y_xx_yy, g_zz_0_0_0_z_y_xx_yz, g_zz_0_0_0_z_y_xx_zz, g_zzz_y_xx_xx, g_zzz_y_xx_xy, g_zzz_y_xx_xz, g_zzz_y_xx_yy, g_zzz_y_xx_yz, g_zzz_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_y_xx_xx[i] = -6.0 * g_z_y_xx_xx[i] * a_exp + 4.0 * g_zzz_y_xx_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_xx_xy[i] = -6.0 * g_z_y_xx_xy[i] * a_exp + 4.0 * g_zzz_y_xx_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_xx_xz[i] = -6.0 * g_z_y_xx_xz[i] * a_exp + 4.0 * g_zzz_y_xx_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_xx_yy[i] = -6.0 * g_z_y_xx_yy[i] * a_exp + 4.0 * g_zzz_y_xx_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_xx_yz[i] = -6.0 * g_z_y_xx_yz[i] * a_exp + 4.0 * g_zzz_y_xx_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_xx_zz[i] = -6.0 * g_z_y_xx_zz[i] * a_exp + 4.0 * g_zzz_y_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (1878-1884)

    #pragma omp simd aligned(g_z_y_xy_xx, g_z_y_xy_xy, g_z_y_xy_xz, g_z_y_xy_yy, g_z_y_xy_yz, g_z_y_xy_zz, g_zz_0_0_0_z_y_xy_xx, g_zz_0_0_0_z_y_xy_xy, g_zz_0_0_0_z_y_xy_xz, g_zz_0_0_0_z_y_xy_yy, g_zz_0_0_0_z_y_xy_yz, g_zz_0_0_0_z_y_xy_zz, g_zzz_y_xy_xx, g_zzz_y_xy_xy, g_zzz_y_xy_xz, g_zzz_y_xy_yy, g_zzz_y_xy_yz, g_zzz_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_y_xy_xx[i] = -6.0 * g_z_y_xy_xx[i] * a_exp + 4.0 * g_zzz_y_xy_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_xy_xy[i] = -6.0 * g_z_y_xy_xy[i] * a_exp + 4.0 * g_zzz_y_xy_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_xy_xz[i] = -6.0 * g_z_y_xy_xz[i] * a_exp + 4.0 * g_zzz_y_xy_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_xy_yy[i] = -6.0 * g_z_y_xy_yy[i] * a_exp + 4.0 * g_zzz_y_xy_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_xy_yz[i] = -6.0 * g_z_y_xy_yz[i] * a_exp + 4.0 * g_zzz_y_xy_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_xy_zz[i] = -6.0 * g_z_y_xy_zz[i] * a_exp + 4.0 * g_zzz_y_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1884-1890)

    #pragma omp simd aligned(g_z_y_xz_xx, g_z_y_xz_xy, g_z_y_xz_xz, g_z_y_xz_yy, g_z_y_xz_yz, g_z_y_xz_zz, g_zz_0_0_0_z_y_xz_xx, g_zz_0_0_0_z_y_xz_xy, g_zz_0_0_0_z_y_xz_xz, g_zz_0_0_0_z_y_xz_yy, g_zz_0_0_0_z_y_xz_yz, g_zz_0_0_0_z_y_xz_zz, g_zzz_y_xz_xx, g_zzz_y_xz_xy, g_zzz_y_xz_xz, g_zzz_y_xz_yy, g_zzz_y_xz_yz, g_zzz_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_y_xz_xx[i] = -6.0 * g_z_y_xz_xx[i] * a_exp + 4.0 * g_zzz_y_xz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_xz_xy[i] = -6.0 * g_z_y_xz_xy[i] * a_exp + 4.0 * g_zzz_y_xz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_xz_xz[i] = -6.0 * g_z_y_xz_xz[i] * a_exp + 4.0 * g_zzz_y_xz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_xz_yy[i] = -6.0 * g_z_y_xz_yy[i] * a_exp + 4.0 * g_zzz_y_xz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_xz_yz[i] = -6.0 * g_z_y_xz_yz[i] * a_exp + 4.0 * g_zzz_y_xz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_xz_zz[i] = -6.0 * g_z_y_xz_zz[i] * a_exp + 4.0 * g_zzz_y_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1890-1896)

    #pragma omp simd aligned(g_z_y_yy_xx, g_z_y_yy_xy, g_z_y_yy_xz, g_z_y_yy_yy, g_z_y_yy_yz, g_z_y_yy_zz, g_zz_0_0_0_z_y_yy_xx, g_zz_0_0_0_z_y_yy_xy, g_zz_0_0_0_z_y_yy_xz, g_zz_0_0_0_z_y_yy_yy, g_zz_0_0_0_z_y_yy_yz, g_zz_0_0_0_z_y_yy_zz, g_zzz_y_yy_xx, g_zzz_y_yy_xy, g_zzz_y_yy_xz, g_zzz_y_yy_yy, g_zzz_y_yy_yz, g_zzz_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_y_yy_xx[i] = -6.0 * g_z_y_yy_xx[i] * a_exp + 4.0 * g_zzz_y_yy_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_yy_xy[i] = -6.0 * g_z_y_yy_xy[i] * a_exp + 4.0 * g_zzz_y_yy_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_yy_xz[i] = -6.0 * g_z_y_yy_xz[i] * a_exp + 4.0 * g_zzz_y_yy_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_yy_yy[i] = -6.0 * g_z_y_yy_yy[i] * a_exp + 4.0 * g_zzz_y_yy_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_yy_yz[i] = -6.0 * g_z_y_yy_yz[i] * a_exp + 4.0 * g_zzz_y_yy_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_yy_zz[i] = -6.0 * g_z_y_yy_zz[i] * a_exp + 4.0 * g_zzz_y_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1896-1902)

    #pragma omp simd aligned(g_z_y_yz_xx, g_z_y_yz_xy, g_z_y_yz_xz, g_z_y_yz_yy, g_z_y_yz_yz, g_z_y_yz_zz, g_zz_0_0_0_z_y_yz_xx, g_zz_0_0_0_z_y_yz_xy, g_zz_0_0_0_z_y_yz_xz, g_zz_0_0_0_z_y_yz_yy, g_zz_0_0_0_z_y_yz_yz, g_zz_0_0_0_z_y_yz_zz, g_zzz_y_yz_xx, g_zzz_y_yz_xy, g_zzz_y_yz_xz, g_zzz_y_yz_yy, g_zzz_y_yz_yz, g_zzz_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_y_yz_xx[i] = -6.0 * g_z_y_yz_xx[i] * a_exp + 4.0 * g_zzz_y_yz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_yz_xy[i] = -6.0 * g_z_y_yz_xy[i] * a_exp + 4.0 * g_zzz_y_yz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_yz_xz[i] = -6.0 * g_z_y_yz_xz[i] * a_exp + 4.0 * g_zzz_y_yz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_yz_yy[i] = -6.0 * g_z_y_yz_yy[i] * a_exp + 4.0 * g_zzz_y_yz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_yz_yz[i] = -6.0 * g_z_y_yz_yz[i] * a_exp + 4.0 * g_zzz_y_yz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_yz_zz[i] = -6.0 * g_z_y_yz_zz[i] * a_exp + 4.0 * g_zzz_y_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1902-1908)

    #pragma omp simd aligned(g_z_y_zz_xx, g_z_y_zz_xy, g_z_y_zz_xz, g_z_y_zz_yy, g_z_y_zz_yz, g_z_y_zz_zz, g_zz_0_0_0_z_y_zz_xx, g_zz_0_0_0_z_y_zz_xy, g_zz_0_0_0_z_y_zz_xz, g_zz_0_0_0_z_y_zz_yy, g_zz_0_0_0_z_y_zz_yz, g_zz_0_0_0_z_y_zz_zz, g_zzz_y_zz_xx, g_zzz_y_zz_xy, g_zzz_y_zz_xz, g_zzz_y_zz_yy, g_zzz_y_zz_yz, g_zzz_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_y_zz_xx[i] = -6.0 * g_z_y_zz_xx[i] * a_exp + 4.0 * g_zzz_y_zz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_zz_xy[i] = -6.0 * g_z_y_zz_xy[i] * a_exp + 4.0 * g_zzz_y_zz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_zz_xz[i] = -6.0 * g_z_y_zz_xz[i] * a_exp + 4.0 * g_zzz_y_zz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_zz_yy[i] = -6.0 * g_z_y_zz_yy[i] * a_exp + 4.0 * g_zzz_y_zz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_zz_yz[i] = -6.0 * g_z_y_zz_yz[i] * a_exp + 4.0 * g_zzz_y_zz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_zz_zz[i] = -6.0 * g_z_y_zz_zz[i] * a_exp + 4.0 * g_zzz_y_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1908-1914)

    #pragma omp simd aligned(g_z_z_xx_xx, g_z_z_xx_xy, g_z_z_xx_xz, g_z_z_xx_yy, g_z_z_xx_yz, g_z_z_xx_zz, g_zz_0_0_0_z_z_xx_xx, g_zz_0_0_0_z_z_xx_xy, g_zz_0_0_0_z_z_xx_xz, g_zz_0_0_0_z_z_xx_yy, g_zz_0_0_0_z_z_xx_yz, g_zz_0_0_0_z_z_xx_zz, g_zzz_z_xx_xx, g_zzz_z_xx_xy, g_zzz_z_xx_xz, g_zzz_z_xx_yy, g_zzz_z_xx_yz, g_zzz_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_z_xx_xx[i] = -6.0 * g_z_z_xx_xx[i] * a_exp + 4.0 * g_zzz_z_xx_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_xx_xy[i] = -6.0 * g_z_z_xx_xy[i] * a_exp + 4.0 * g_zzz_z_xx_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_xx_xz[i] = -6.0 * g_z_z_xx_xz[i] * a_exp + 4.0 * g_zzz_z_xx_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_xx_yy[i] = -6.0 * g_z_z_xx_yy[i] * a_exp + 4.0 * g_zzz_z_xx_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_xx_yz[i] = -6.0 * g_z_z_xx_yz[i] * a_exp + 4.0 * g_zzz_z_xx_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_xx_zz[i] = -6.0 * g_z_z_xx_zz[i] * a_exp + 4.0 * g_zzz_z_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (1914-1920)

    #pragma omp simd aligned(g_z_z_xy_xx, g_z_z_xy_xy, g_z_z_xy_xz, g_z_z_xy_yy, g_z_z_xy_yz, g_z_z_xy_zz, g_zz_0_0_0_z_z_xy_xx, g_zz_0_0_0_z_z_xy_xy, g_zz_0_0_0_z_z_xy_xz, g_zz_0_0_0_z_z_xy_yy, g_zz_0_0_0_z_z_xy_yz, g_zz_0_0_0_z_z_xy_zz, g_zzz_z_xy_xx, g_zzz_z_xy_xy, g_zzz_z_xy_xz, g_zzz_z_xy_yy, g_zzz_z_xy_yz, g_zzz_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_z_xy_xx[i] = -6.0 * g_z_z_xy_xx[i] * a_exp + 4.0 * g_zzz_z_xy_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_xy_xy[i] = -6.0 * g_z_z_xy_xy[i] * a_exp + 4.0 * g_zzz_z_xy_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_xy_xz[i] = -6.0 * g_z_z_xy_xz[i] * a_exp + 4.0 * g_zzz_z_xy_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_xy_yy[i] = -6.0 * g_z_z_xy_yy[i] * a_exp + 4.0 * g_zzz_z_xy_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_xy_yz[i] = -6.0 * g_z_z_xy_yz[i] * a_exp + 4.0 * g_zzz_z_xy_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_xy_zz[i] = -6.0 * g_z_z_xy_zz[i] * a_exp + 4.0 * g_zzz_z_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1920-1926)

    #pragma omp simd aligned(g_z_z_xz_xx, g_z_z_xz_xy, g_z_z_xz_xz, g_z_z_xz_yy, g_z_z_xz_yz, g_z_z_xz_zz, g_zz_0_0_0_z_z_xz_xx, g_zz_0_0_0_z_z_xz_xy, g_zz_0_0_0_z_z_xz_xz, g_zz_0_0_0_z_z_xz_yy, g_zz_0_0_0_z_z_xz_yz, g_zz_0_0_0_z_z_xz_zz, g_zzz_z_xz_xx, g_zzz_z_xz_xy, g_zzz_z_xz_xz, g_zzz_z_xz_yy, g_zzz_z_xz_yz, g_zzz_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_z_xz_xx[i] = -6.0 * g_z_z_xz_xx[i] * a_exp + 4.0 * g_zzz_z_xz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_xz_xy[i] = -6.0 * g_z_z_xz_xy[i] * a_exp + 4.0 * g_zzz_z_xz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_xz_xz[i] = -6.0 * g_z_z_xz_xz[i] * a_exp + 4.0 * g_zzz_z_xz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_xz_yy[i] = -6.0 * g_z_z_xz_yy[i] * a_exp + 4.0 * g_zzz_z_xz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_xz_yz[i] = -6.0 * g_z_z_xz_yz[i] * a_exp + 4.0 * g_zzz_z_xz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_xz_zz[i] = -6.0 * g_z_z_xz_zz[i] * a_exp + 4.0 * g_zzz_z_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1926-1932)

    #pragma omp simd aligned(g_z_z_yy_xx, g_z_z_yy_xy, g_z_z_yy_xz, g_z_z_yy_yy, g_z_z_yy_yz, g_z_z_yy_zz, g_zz_0_0_0_z_z_yy_xx, g_zz_0_0_0_z_z_yy_xy, g_zz_0_0_0_z_z_yy_xz, g_zz_0_0_0_z_z_yy_yy, g_zz_0_0_0_z_z_yy_yz, g_zz_0_0_0_z_z_yy_zz, g_zzz_z_yy_xx, g_zzz_z_yy_xy, g_zzz_z_yy_xz, g_zzz_z_yy_yy, g_zzz_z_yy_yz, g_zzz_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_z_yy_xx[i] = -6.0 * g_z_z_yy_xx[i] * a_exp + 4.0 * g_zzz_z_yy_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_yy_xy[i] = -6.0 * g_z_z_yy_xy[i] * a_exp + 4.0 * g_zzz_z_yy_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_yy_xz[i] = -6.0 * g_z_z_yy_xz[i] * a_exp + 4.0 * g_zzz_z_yy_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_yy_yy[i] = -6.0 * g_z_z_yy_yy[i] * a_exp + 4.0 * g_zzz_z_yy_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_yy_yz[i] = -6.0 * g_z_z_yy_yz[i] * a_exp + 4.0 * g_zzz_z_yy_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_yy_zz[i] = -6.0 * g_z_z_yy_zz[i] * a_exp + 4.0 * g_zzz_z_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1932-1938)

    #pragma omp simd aligned(g_z_z_yz_xx, g_z_z_yz_xy, g_z_z_yz_xz, g_z_z_yz_yy, g_z_z_yz_yz, g_z_z_yz_zz, g_zz_0_0_0_z_z_yz_xx, g_zz_0_0_0_z_z_yz_xy, g_zz_0_0_0_z_z_yz_xz, g_zz_0_0_0_z_z_yz_yy, g_zz_0_0_0_z_z_yz_yz, g_zz_0_0_0_z_z_yz_zz, g_zzz_z_yz_xx, g_zzz_z_yz_xy, g_zzz_z_yz_xz, g_zzz_z_yz_yy, g_zzz_z_yz_yz, g_zzz_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_z_yz_xx[i] = -6.0 * g_z_z_yz_xx[i] * a_exp + 4.0 * g_zzz_z_yz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_yz_xy[i] = -6.0 * g_z_z_yz_xy[i] * a_exp + 4.0 * g_zzz_z_yz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_yz_xz[i] = -6.0 * g_z_z_yz_xz[i] * a_exp + 4.0 * g_zzz_z_yz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_yz_yy[i] = -6.0 * g_z_z_yz_yy[i] * a_exp + 4.0 * g_zzz_z_yz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_yz_yz[i] = -6.0 * g_z_z_yz_yz[i] * a_exp + 4.0 * g_zzz_z_yz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_yz_zz[i] = -6.0 * g_z_z_yz_zz[i] * a_exp + 4.0 * g_zzz_z_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1938-1944)

    #pragma omp simd aligned(g_z_z_zz_xx, g_z_z_zz_xy, g_z_z_zz_xz, g_z_z_zz_yy, g_z_z_zz_yz, g_z_z_zz_zz, g_zz_0_0_0_z_z_zz_xx, g_zz_0_0_0_z_z_zz_xy, g_zz_0_0_0_z_z_zz_xz, g_zz_0_0_0_z_z_zz_yy, g_zz_0_0_0_z_z_zz_yz, g_zz_0_0_0_z_z_zz_zz, g_zzz_z_zz_xx, g_zzz_z_zz_xy, g_zzz_z_zz_xz, g_zzz_z_zz_yy, g_zzz_z_zz_yz, g_zzz_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_z_zz_xx[i] = -6.0 * g_z_z_zz_xx[i] * a_exp + 4.0 * g_zzz_z_zz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_zz_xy[i] = -6.0 * g_z_z_zz_xy[i] * a_exp + 4.0 * g_zzz_z_zz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_zz_xz[i] = -6.0 * g_z_z_zz_xz[i] * a_exp + 4.0 * g_zzz_z_zz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_zz_yy[i] = -6.0 * g_z_z_zz_yy[i] * a_exp + 4.0 * g_zzz_z_zz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_zz_yz[i] = -6.0 * g_z_z_zz_yz[i] * a_exp + 4.0 * g_zzz_z_zz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_zz_zz[i] = -6.0 * g_z_z_zz_zz[i] * a_exp + 4.0 * g_zzz_z_zz_zz[i] * a_exp * a_exp;
    }
}

} // t4c_geom namespace


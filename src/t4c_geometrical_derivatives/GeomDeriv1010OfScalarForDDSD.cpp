#include "GeomDeriv1010OfScalarForDDSD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1010_ddsd_0(CSimdArray<double>& buffer_1010_ddsd,
                     const CSimdArray<double>& buffer_pdpd,
                     const CSimdArray<double>& buffer_fdpd,
                     const double a_exp,
                     const double* c_exps) -> void
{
    const auto ndims = buffer_1010_ddsd.number_of_columns();

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

    /// Set up components of auxilary buffer : buffer_fdpd

    auto g_xxx_xx_x_xx = buffer_fdpd[0];

    auto g_xxx_xx_x_xy = buffer_fdpd[1];

    auto g_xxx_xx_x_xz = buffer_fdpd[2];

    auto g_xxx_xx_x_yy = buffer_fdpd[3];

    auto g_xxx_xx_x_yz = buffer_fdpd[4];

    auto g_xxx_xx_x_zz = buffer_fdpd[5];

    auto g_xxx_xx_y_xx = buffer_fdpd[6];

    auto g_xxx_xx_y_xy = buffer_fdpd[7];

    auto g_xxx_xx_y_xz = buffer_fdpd[8];

    auto g_xxx_xx_y_yy = buffer_fdpd[9];

    auto g_xxx_xx_y_yz = buffer_fdpd[10];

    auto g_xxx_xx_y_zz = buffer_fdpd[11];

    auto g_xxx_xx_z_xx = buffer_fdpd[12];

    auto g_xxx_xx_z_xy = buffer_fdpd[13];

    auto g_xxx_xx_z_xz = buffer_fdpd[14];

    auto g_xxx_xx_z_yy = buffer_fdpd[15];

    auto g_xxx_xx_z_yz = buffer_fdpd[16];

    auto g_xxx_xx_z_zz = buffer_fdpd[17];

    auto g_xxx_xy_x_xx = buffer_fdpd[18];

    auto g_xxx_xy_x_xy = buffer_fdpd[19];

    auto g_xxx_xy_x_xz = buffer_fdpd[20];

    auto g_xxx_xy_x_yy = buffer_fdpd[21];

    auto g_xxx_xy_x_yz = buffer_fdpd[22];

    auto g_xxx_xy_x_zz = buffer_fdpd[23];

    auto g_xxx_xy_y_xx = buffer_fdpd[24];

    auto g_xxx_xy_y_xy = buffer_fdpd[25];

    auto g_xxx_xy_y_xz = buffer_fdpd[26];

    auto g_xxx_xy_y_yy = buffer_fdpd[27];

    auto g_xxx_xy_y_yz = buffer_fdpd[28];

    auto g_xxx_xy_y_zz = buffer_fdpd[29];

    auto g_xxx_xy_z_xx = buffer_fdpd[30];

    auto g_xxx_xy_z_xy = buffer_fdpd[31];

    auto g_xxx_xy_z_xz = buffer_fdpd[32];

    auto g_xxx_xy_z_yy = buffer_fdpd[33];

    auto g_xxx_xy_z_yz = buffer_fdpd[34];

    auto g_xxx_xy_z_zz = buffer_fdpd[35];

    auto g_xxx_xz_x_xx = buffer_fdpd[36];

    auto g_xxx_xz_x_xy = buffer_fdpd[37];

    auto g_xxx_xz_x_xz = buffer_fdpd[38];

    auto g_xxx_xz_x_yy = buffer_fdpd[39];

    auto g_xxx_xz_x_yz = buffer_fdpd[40];

    auto g_xxx_xz_x_zz = buffer_fdpd[41];

    auto g_xxx_xz_y_xx = buffer_fdpd[42];

    auto g_xxx_xz_y_xy = buffer_fdpd[43];

    auto g_xxx_xz_y_xz = buffer_fdpd[44];

    auto g_xxx_xz_y_yy = buffer_fdpd[45];

    auto g_xxx_xz_y_yz = buffer_fdpd[46];

    auto g_xxx_xz_y_zz = buffer_fdpd[47];

    auto g_xxx_xz_z_xx = buffer_fdpd[48];

    auto g_xxx_xz_z_xy = buffer_fdpd[49];

    auto g_xxx_xz_z_xz = buffer_fdpd[50];

    auto g_xxx_xz_z_yy = buffer_fdpd[51];

    auto g_xxx_xz_z_yz = buffer_fdpd[52];

    auto g_xxx_xz_z_zz = buffer_fdpd[53];

    auto g_xxx_yy_x_xx = buffer_fdpd[54];

    auto g_xxx_yy_x_xy = buffer_fdpd[55];

    auto g_xxx_yy_x_xz = buffer_fdpd[56];

    auto g_xxx_yy_x_yy = buffer_fdpd[57];

    auto g_xxx_yy_x_yz = buffer_fdpd[58];

    auto g_xxx_yy_x_zz = buffer_fdpd[59];

    auto g_xxx_yy_y_xx = buffer_fdpd[60];

    auto g_xxx_yy_y_xy = buffer_fdpd[61];

    auto g_xxx_yy_y_xz = buffer_fdpd[62];

    auto g_xxx_yy_y_yy = buffer_fdpd[63];

    auto g_xxx_yy_y_yz = buffer_fdpd[64];

    auto g_xxx_yy_y_zz = buffer_fdpd[65];

    auto g_xxx_yy_z_xx = buffer_fdpd[66];

    auto g_xxx_yy_z_xy = buffer_fdpd[67];

    auto g_xxx_yy_z_xz = buffer_fdpd[68];

    auto g_xxx_yy_z_yy = buffer_fdpd[69];

    auto g_xxx_yy_z_yz = buffer_fdpd[70];

    auto g_xxx_yy_z_zz = buffer_fdpd[71];

    auto g_xxx_yz_x_xx = buffer_fdpd[72];

    auto g_xxx_yz_x_xy = buffer_fdpd[73];

    auto g_xxx_yz_x_xz = buffer_fdpd[74];

    auto g_xxx_yz_x_yy = buffer_fdpd[75];

    auto g_xxx_yz_x_yz = buffer_fdpd[76];

    auto g_xxx_yz_x_zz = buffer_fdpd[77];

    auto g_xxx_yz_y_xx = buffer_fdpd[78];

    auto g_xxx_yz_y_xy = buffer_fdpd[79];

    auto g_xxx_yz_y_xz = buffer_fdpd[80];

    auto g_xxx_yz_y_yy = buffer_fdpd[81];

    auto g_xxx_yz_y_yz = buffer_fdpd[82];

    auto g_xxx_yz_y_zz = buffer_fdpd[83];

    auto g_xxx_yz_z_xx = buffer_fdpd[84];

    auto g_xxx_yz_z_xy = buffer_fdpd[85];

    auto g_xxx_yz_z_xz = buffer_fdpd[86];

    auto g_xxx_yz_z_yy = buffer_fdpd[87];

    auto g_xxx_yz_z_yz = buffer_fdpd[88];

    auto g_xxx_yz_z_zz = buffer_fdpd[89];

    auto g_xxx_zz_x_xx = buffer_fdpd[90];

    auto g_xxx_zz_x_xy = buffer_fdpd[91];

    auto g_xxx_zz_x_xz = buffer_fdpd[92];

    auto g_xxx_zz_x_yy = buffer_fdpd[93];

    auto g_xxx_zz_x_yz = buffer_fdpd[94];

    auto g_xxx_zz_x_zz = buffer_fdpd[95];

    auto g_xxx_zz_y_xx = buffer_fdpd[96];

    auto g_xxx_zz_y_xy = buffer_fdpd[97];

    auto g_xxx_zz_y_xz = buffer_fdpd[98];

    auto g_xxx_zz_y_yy = buffer_fdpd[99];

    auto g_xxx_zz_y_yz = buffer_fdpd[100];

    auto g_xxx_zz_y_zz = buffer_fdpd[101];

    auto g_xxx_zz_z_xx = buffer_fdpd[102];

    auto g_xxx_zz_z_xy = buffer_fdpd[103];

    auto g_xxx_zz_z_xz = buffer_fdpd[104];

    auto g_xxx_zz_z_yy = buffer_fdpd[105];

    auto g_xxx_zz_z_yz = buffer_fdpd[106];

    auto g_xxx_zz_z_zz = buffer_fdpd[107];

    auto g_xxy_xx_x_xx = buffer_fdpd[108];

    auto g_xxy_xx_x_xy = buffer_fdpd[109];

    auto g_xxy_xx_x_xz = buffer_fdpd[110];

    auto g_xxy_xx_x_yy = buffer_fdpd[111];

    auto g_xxy_xx_x_yz = buffer_fdpd[112];

    auto g_xxy_xx_x_zz = buffer_fdpd[113];

    auto g_xxy_xx_y_xx = buffer_fdpd[114];

    auto g_xxy_xx_y_xy = buffer_fdpd[115];

    auto g_xxy_xx_y_xz = buffer_fdpd[116];

    auto g_xxy_xx_y_yy = buffer_fdpd[117];

    auto g_xxy_xx_y_yz = buffer_fdpd[118];

    auto g_xxy_xx_y_zz = buffer_fdpd[119];

    auto g_xxy_xx_z_xx = buffer_fdpd[120];

    auto g_xxy_xx_z_xy = buffer_fdpd[121];

    auto g_xxy_xx_z_xz = buffer_fdpd[122];

    auto g_xxy_xx_z_yy = buffer_fdpd[123];

    auto g_xxy_xx_z_yz = buffer_fdpd[124];

    auto g_xxy_xx_z_zz = buffer_fdpd[125];

    auto g_xxy_xy_x_xx = buffer_fdpd[126];

    auto g_xxy_xy_x_xy = buffer_fdpd[127];

    auto g_xxy_xy_x_xz = buffer_fdpd[128];

    auto g_xxy_xy_x_yy = buffer_fdpd[129];

    auto g_xxy_xy_x_yz = buffer_fdpd[130];

    auto g_xxy_xy_x_zz = buffer_fdpd[131];

    auto g_xxy_xy_y_xx = buffer_fdpd[132];

    auto g_xxy_xy_y_xy = buffer_fdpd[133];

    auto g_xxy_xy_y_xz = buffer_fdpd[134];

    auto g_xxy_xy_y_yy = buffer_fdpd[135];

    auto g_xxy_xy_y_yz = buffer_fdpd[136];

    auto g_xxy_xy_y_zz = buffer_fdpd[137];

    auto g_xxy_xy_z_xx = buffer_fdpd[138];

    auto g_xxy_xy_z_xy = buffer_fdpd[139];

    auto g_xxy_xy_z_xz = buffer_fdpd[140];

    auto g_xxy_xy_z_yy = buffer_fdpd[141];

    auto g_xxy_xy_z_yz = buffer_fdpd[142];

    auto g_xxy_xy_z_zz = buffer_fdpd[143];

    auto g_xxy_xz_x_xx = buffer_fdpd[144];

    auto g_xxy_xz_x_xy = buffer_fdpd[145];

    auto g_xxy_xz_x_xz = buffer_fdpd[146];

    auto g_xxy_xz_x_yy = buffer_fdpd[147];

    auto g_xxy_xz_x_yz = buffer_fdpd[148];

    auto g_xxy_xz_x_zz = buffer_fdpd[149];

    auto g_xxy_xz_y_xx = buffer_fdpd[150];

    auto g_xxy_xz_y_xy = buffer_fdpd[151];

    auto g_xxy_xz_y_xz = buffer_fdpd[152];

    auto g_xxy_xz_y_yy = buffer_fdpd[153];

    auto g_xxy_xz_y_yz = buffer_fdpd[154];

    auto g_xxy_xz_y_zz = buffer_fdpd[155];

    auto g_xxy_xz_z_xx = buffer_fdpd[156];

    auto g_xxy_xz_z_xy = buffer_fdpd[157];

    auto g_xxy_xz_z_xz = buffer_fdpd[158];

    auto g_xxy_xz_z_yy = buffer_fdpd[159];

    auto g_xxy_xz_z_yz = buffer_fdpd[160];

    auto g_xxy_xz_z_zz = buffer_fdpd[161];

    auto g_xxy_yy_x_xx = buffer_fdpd[162];

    auto g_xxy_yy_x_xy = buffer_fdpd[163];

    auto g_xxy_yy_x_xz = buffer_fdpd[164];

    auto g_xxy_yy_x_yy = buffer_fdpd[165];

    auto g_xxy_yy_x_yz = buffer_fdpd[166];

    auto g_xxy_yy_x_zz = buffer_fdpd[167];

    auto g_xxy_yy_y_xx = buffer_fdpd[168];

    auto g_xxy_yy_y_xy = buffer_fdpd[169];

    auto g_xxy_yy_y_xz = buffer_fdpd[170];

    auto g_xxy_yy_y_yy = buffer_fdpd[171];

    auto g_xxy_yy_y_yz = buffer_fdpd[172];

    auto g_xxy_yy_y_zz = buffer_fdpd[173];

    auto g_xxy_yy_z_xx = buffer_fdpd[174];

    auto g_xxy_yy_z_xy = buffer_fdpd[175];

    auto g_xxy_yy_z_xz = buffer_fdpd[176];

    auto g_xxy_yy_z_yy = buffer_fdpd[177];

    auto g_xxy_yy_z_yz = buffer_fdpd[178];

    auto g_xxy_yy_z_zz = buffer_fdpd[179];

    auto g_xxy_yz_x_xx = buffer_fdpd[180];

    auto g_xxy_yz_x_xy = buffer_fdpd[181];

    auto g_xxy_yz_x_xz = buffer_fdpd[182];

    auto g_xxy_yz_x_yy = buffer_fdpd[183];

    auto g_xxy_yz_x_yz = buffer_fdpd[184];

    auto g_xxy_yz_x_zz = buffer_fdpd[185];

    auto g_xxy_yz_y_xx = buffer_fdpd[186];

    auto g_xxy_yz_y_xy = buffer_fdpd[187];

    auto g_xxy_yz_y_xz = buffer_fdpd[188];

    auto g_xxy_yz_y_yy = buffer_fdpd[189];

    auto g_xxy_yz_y_yz = buffer_fdpd[190];

    auto g_xxy_yz_y_zz = buffer_fdpd[191];

    auto g_xxy_yz_z_xx = buffer_fdpd[192];

    auto g_xxy_yz_z_xy = buffer_fdpd[193];

    auto g_xxy_yz_z_xz = buffer_fdpd[194];

    auto g_xxy_yz_z_yy = buffer_fdpd[195];

    auto g_xxy_yz_z_yz = buffer_fdpd[196];

    auto g_xxy_yz_z_zz = buffer_fdpd[197];

    auto g_xxy_zz_x_xx = buffer_fdpd[198];

    auto g_xxy_zz_x_xy = buffer_fdpd[199];

    auto g_xxy_zz_x_xz = buffer_fdpd[200];

    auto g_xxy_zz_x_yy = buffer_fdpd[201];

    auto g_xxy_zz_x_yz = buffer_fdpd[202];

    auto g_xxy_zz_x_zz = buffer_fdpd[203];

    auto g_xxy_zz_y_xx = buffer_fdpd[204];

    auto g_xxy_zz_y_xy = buffer_fdpd[205];

    auto g_xxy_zz_y_xz = buffer_fdpd[206];

    auto g_xxy_zz_y_yy = buffer_fdpd[207];

    auto g_xxy_zz_y_yz = buffer_fdpd[208];

    auto g_xxy_zz_y_zz = buffer_fdpd[209];

    auto g_xxy_zz_z_xx = buffer_fdpd[210];

    auto g_xxy_zz_z_xy = buffer_fdpd[211];

    auto g_xxy_zz_z_xz = buffer_fdpd[212];

    auto g_xxy_zz_z_yy = buffer_fdpd[213];

    auto g_xxy_zz_z_yz = buffer_fdpd[214];

    auto g_xxy_zz_z_zz = buffer_fdpd[215];

    auto g_xxz_xx_x_xx = buffer_fdpd[216];

    auto g_xxz_xx_x_xy = buffer_fdpd[217];

    auto g_xxz_xx_x_xz = buffer_fdpd[218];

    auto g_xxz_xx_x_yy = buffer_fdpd[219];

    auto g_xxz_xx_x_yz = buffer_fdpd[220];

    auto g_xxz_xx_x_zz = buffer_fdpd[221];

    auto g_xxz_xx_y_xx = buffer_fdpd[222];

    auto g_xxz_xx_y_xy = buffer_fdpd[223];

    auto g_xxz_xx_y_xz = buffer_fdpd[224];

    auto g_xxz_xx_y_yy = buffer_fdpd[225];

    auto g_xxz_xx_y_yz = buffer_fdpd[226];

    auto g_xxz_xx_y_zz = buffer_fdpd[227];

    auto g_xxz_xx_z_xx = buffer_fdpd[228];

    auto g_xxz_xx_z_xy = buffer_fdpd[229];

    auto g_xxz_xx_z_xz = buffer_fdpd[230];

    auto g_xxz_xx_z_yy = buffer_fdpd[231];

    auto g_xxz_xx_z_yz = buffer_fdpd[232];

    auto g_xxz_xx_z_zz = buffer_fdpd[233];

    auto g_xxz_xy_x_xx = buffer_fdpd[234];

    auto g_xxz_xy_x_xy = buffer_fdpd[235];

    auto g_xxz_xy_x_xz = buffer_fdpd[236];

    auto g_xxz_xy_x_yy = buffer_fdpd[237];

    auto g_xxz_xy_x_yz = buffer_fdpd[238];

    auto g_xxz_xy_x_zz = buffer_fdpd[239];

    auto g_xxz_xy_y_xx = buffer_fdpd[240];

    auto g_xxz_xy_y_xy = buffer_fdpd[241];

    auto g_xxz_xy_y_xz = buffer_fdpd[242];

    auto g_xxz_xy_y_yy = buffer_fdpd[243];

    auto g_xxz_xy_y_yz = buffer_fdpd[244];

    auto g_xxz_xy_y_zz = buffer_fdpd[245];

    auto g_xxz_xy_z_xx = buffer_fdpd[246];

    auto g_xxz_xy_z_xy = buffer_fdpd[247];

    auto g_xxz_xy_z_xz = buffer_fdpd[248];

    auto g_xxz_xy_z_yy = buffer_fdpd[249];

    auto g_xxz_xy_z_yz = buffer_fdpd[250];

    auto g_xxz_xy_z_zz = buffer_fdpd[251];

    auto g_xxz_xz_x_xx = buffer_fdpd[252];

    auto g_xxz_xz_x_xy = buffer_fdpd[253];

    auto g_xxz_xz_x_xz = buffer_fdpd[254];

    auto g_xxz_xz_x_yy = buffer_fdpd[255];

    auto g_xxz_xz_x_yz = buffer_fdpd[256];

    auto g_xxz_xz_x_zz = buffer_fdpd[257];

    auto g_xxz_xz_y_xx = buffer_fdpd[258];

    auto g_xxz_xz_y_xy = buffer_fdpd[259];

    auto g_xxz_xz_y_xz = buffer_fdpd[260];

    auto g_xxz_xz_y_yy = buffer_fdpd[261];

    auto g_xxz_xz_y_yz = buffer_fdpd[262];

    auto g_xxz_xz_y_zz = buffer_fdpd[263];

    auto g_xxz_xz_z_xx = buffer_fdpd[264];

    auto g_xxz_xz_z_xy = buffer_fdpd[265];

    auto g_xxz_xz_z_xz = buffer_fdpd[266];

    auto g_xxz_xz_z_yy = buffer_fdpd[267];

    auto g_xxz_xz_z_yz = buffer_fdpd[268];

    auto g_xxz_xz_z_zz = buffer_fdpd[269];

    auto g_xxz_yy_x_xx = buffer_fdpd[270];

    auto g_xxz_yy_x_xy = buffer_fdpd[271];

    auto g_xxz_yy_x_xz = buffer_fdpd[272];

    auto g_xxz_yy_x_yy = buffer_fdpd[273];

    auto g_xxz_yy_x_yz = buffer_fdpd[274];

    auto g_xxz_yy_x_zz = buffer_fdpd[275];

    auto g_xxz_yy_y_xx = buffer_fdpd[276];

    auto g_xxz_yy_y_xy = buffer_fdpd[277];

    auto g_xxz_yy_y_xz = buffer_fdpd[278];

    auto g_xxz_yy_y_yy = buffer_fdpd[279];

    auto g_xxz_yy_y_yz = buffer_fdpd[280];

    auto g_xxz_yy_y_zz = buffer_fdpd[281];

    auto g_xxz_yy_z_xx = buffer_fdpd[282];

    auto g_xxz_yy_z_xy = buffer_fdpd[283];

    auto g_xxz_yy_z_xz = buffer_fdpd[284];

    auto g_xxz_yy_z_yy = buffer_fdpd[285];

    auto g_xxz_yy_z_yz = buffer_fdpd[286];

    auto g_xxz_yy_z_zz = buffer_fdpd[287];

    auto g_xxz_yz_x_xx = buffer_fdpd[288];

    auto g_xxz_yz_x_xy = buffer_fdpd[289];

    auto g_xxz_yz_x_xz = buffer_fdpd[290];

    auto g_xxz_yz_x_yy = buffer_fdpd[291];

    auto g_xxz_yz_x_yz = buffer_fdpd[292];

    auto g_xxz_yz_x_zz = buffer_fdpd[293];

    auto g_xxz_yz_y_xx = buffer_fdpd[294];

    auto g_xxz_yz_y_xy = buffer_fdpd[295];

    auto g_xxz_yz_y_xz = buffer_fdpd[296];

    auto g_xxz_yz_y_yy = buffer_fdpd[297];

    auto g_xxz_yz_y_yz = buffer_fdpd[298];

    auto g_xxz_yz_y_zz = buffer_fdpd[299];

    auto g_xxz_yz_z_xx = buffer_fdpd[300];

    auto g_xxz_yz_z_xy = buffer_fdpd[301];

    auto g_xxz_yz_z_xz = buffer_fdpd[302];

    auto g_xxz_yz_z_yy = buffer_fdpd[303];

    auto g_xxz_yz_z_yz = buffer_fdpd[304];

    auto g_xxz_yz_z_zz = buffer_fdpd[305];

    auto g_xxz_zz_x_xx = buffer_fdpd[306];

    auto g_xxz_zz_x_xy = buffer_fdpd[307];

    auto g_xxz_zz_x_xz = buffer_fdpd[308];

    auto g_xxz_zz_x_yy = buffer_fdpd[309];

    auto g_xxz_zz_x_yz = buffer_fdpd[310];

    auto g_xxz_zz_x_zz = buffer_fdpd[311];

    auto g_xxz_zz_y_xx = buffer_fdpd[312];

    auto g_xxz_zz_y_xy = buffer_fdpd[313];

    auto g_xxz_zz_y_xz = buffer_fdpd[314];

    auto g_xxz_zz_y_yy = buffer_fdpd[315];

    auto g_xxz_zz_y_yz = buffer_fdpd[316];

    auto g_xxz_zz_y_zz = buffer_fdpd[317];

    auto g_xxz_zz_z_xx = buffer_fdpd[318];

    auto g_xxz_zz_z_xy = buffer_fdpd[319];

    auto g_xxz_zz_z_xz = buffer_fdpd[320];

    auto g_xxz_zz_z_yy = buffer_fdpd[321];

    auto g_xxz_zz_z_yz = buffer_fdpd[322];

    auto g_xxz_zz_z_zz = buffer_fdpd[323];

    auto g_xyy_xx_x_xx = buffer_fdpd[324];

    auto g_xyy_xx_x_xy = buffer_fdpd[325];

    auto g_xyy_xx_x_xz = buffer_fdpd[326];

    auto g_xyy_xx_x_yy = buffer_fdpd[327];

    auto g_xyy_xx_x_yz = buffer_fdpd[328];

    auto g_xyy_xx_x_zz = buffer_fdpd[329];

    auto g_xyy_xx_y_xx = buffer_fdpd[330];

    auto g_xyy_xx_y_xy = buffer_fdpd[331];

    auto g_xyy_xx_y_xz = buffer_fdpd[332];

    auto g_xyy_xx_y_yy = buffer_fdpd[333];

    auto g_xyy_xx_y_yz = buffer_fdpd[334];

    auto g_xyy_xx_y_zz = buffer_fdpd[335];

    auto g_xyy_xx_z_xx = buffer_fdpd[336];

    auto g_xyy_xx_z_xy = buffer_fdpd[337];

    auto g_xyy_xx_z_xz = buffer_fdpd[338];

    auto g_xyy_xx_z_yy = buffer_fdpd[339];

    auto g_xyy_xx_z_yz = buffer_fdpd[340];

    auto g_xyy_xx_z_zz = buffer_fdpd[341];

    auto g_xyy_xy_x_xx = buffer_fdpd[342];

    auto g_xyy_xy_x_xy = buffer_fdpd[343];

    auto g_xyy_xy_x_xz = buffer_fdpd[344];

    auto g_xyy_xy_x_yy = buffer_fdpd[345];

    auto g_xyy_xy_x_yz = buffer_fdpd[346];

    auto g_xyy_xy_x_zz = buffer_fdpd[347];

    auto g_xyy_xy_y_xx = buffer_fdpd[348];

    auto g_xyy_xy_y_xy = buffer_fdpd[349];

    auto g_xyy_xy_y_xz = buffer_fdpd[350];

    auto g_xyy_xy_y_yy = buffer_fdpd[351];

    auto g_xyy_xy_y_yz = buffer_fdpd[352];

    auto g_xyy_xy_y_zz = buffer_fdpd[353];

    auto g_xyy_xy_z_xx = buffer_fdpd[354];

    auto g_xyy_xy_z_xy = buffer_fdpd[355];

    auto g_xyy_xy_z_xz = buffer_fdpd[356];

    auto g_xyy_xy_z_yy = buffer_fdpd[357];

    auto g_xyy_xy_z_yz = buffer_fdpd[358];

    auto g_xyy_xy_z_zz = buffer_fdpd[359];

    auto g_xyy_xz_x_xx = buffer_fdpd[360];

    auto g_xyy_xz_x_xy = buffer_fdpd[361];

    auto g_xyy_xz_x_xz = buffer_fdpd[362];

    auto g_xyy_xz_x_yy = buffer_fdpd[363];

    auto g_xyy_xz_x_yz = buffer_fdpd[364];

    auto g_xyy_xz_x_zz = buffer_fdpd[365];

    auto g_xyy_xz_y_xx = buffer_fdpd[366];

    auto g_xyy_xz_y_xy = buffer_fdpd[367];

    auto g_xyy_xz_y_xz = buffer_fdpd[368];

    auto g_xyy_xz_y_yy = buffer_fdpd[369];

    auto g_xyy_xz_y_yz = buffer_fdpd[370];

    auto g_xyy_xz_y_zz = buffer_fdpd[371];

    auto g_xyy_xz_z_xx = buffer_fdpd[372];

    auto g_xyy_xz_z_xy = buffer_fdpd[373];

    auto g_xyy_xz_z_xz = buffer_fdpd[374];

    auto g_xyy_xz_z_yy = buffer_fdpd[375];

    auto g_xyy_xz_z_yz = buffer_fdpd[376];

    auto g_xyy_xz_z_zz = buffer_fdpd[377];

    auto g_xyy_yy_x_xx = buffer_fdpd[378];

    auto g_xyy_yy_x_xy = buffer_fdpd[379];

    auto g_xyy_yy_x_xz = buffer_fdpd[380];

    auto g_xyy_yy_x_yy = buffer_fdpd[381];

    auto g_xyy_yy_x_yz = buffer_fdpd[382];

    auto g_xyy_yy_x_zz = buffer_fdpd[383];

    auto g_xyy_yy_y_xx = buffer_fdpd[384];

    auto g_xyy_yy_y_xy = buffer_fdpd[385];

    auto g_xyy_yy_y_xz = buffer_fdpd[386];

    auto g_xyy_yy_y_yy = buffer_fdpd[387];

    auto g_xyy_yy_y_yz = buffer_fdpd[388];

    auto g_xyy_yy_y_zz = buffer_fdpd[389];

    auto g_xyy_yy_z_xx = buffer_fdpd[390];

    auto g_xyy_yy_z_xy = buffer_fdpd[391];

    auto g_xyy_yy_z_xz = buffer_fdpd[392];

    auto g_xyy_yy_z_yy = buffer_fdpd[393];

    auto g_xyy_yy_z_yz = buffer_fdpd[394];

    auto g_xyy_yy_z_zz = buffer_fdpd[395];

    auto g_xyy_yz_x_xx = buffer_fdpd[396];

    auto g_xyy_yz_x_xy = buffer_fdpd[397];

    auto g_xyy_yz_x_xz = buffer_fdpd[398];

    auto g_xyy_yz_x_yy = buffer_fdpd[399];

    auto g_xyy_yz_x_yz = buffer_fdpd[400];

    auto g_xyy_yz_x_zz = buffer_fdpd[401];

    auto g_xyy_yz_y_xx = buffer_fdpd[402];

    auto g_xyy_yz_y_xy = buffer_fdpd[403];

    auto g_xyy_yz_y_xz = buffer_fdpd[404];

    auto g_xyy_yz_y_yy = buffer_fdpd[405];

    auto g_xyy_yz_y_yz = buffer_fdpd[406];

    auto g_xyy_yz_y_zz = buffer_fdpd[407];

    auto g_xyy_yz_z_xx = buffer_fdpd[408];

    auto g_xyy_yz_z_xy = buffer_fdpd[409];

    auto g_xyy_yz_z_xz = buffer_fdpd[410];

    auto g_xyy_yz_z_yy = buffer_fdpd[411];

    auto g_xyy_yz_z_yz = buffer_fdpd[412];

    auto g_xyy_yz_z_zz = buffer_fdpd[413];

    auto g_xyy_zz_x_xx = buffer_fdpd[414];

    auto g_xyy_zz_x_xy = buffer_fdpd[415];

    auto g_xyy_zz_x_xz = buffer_fdpd[416];

    auto g_xyy_zz_x_yy = buffer_fdpd[417];

    auto g_xyy_zz_x_yz = buffer_fdpd[418];

    auto g_xyy_zz_x_zz = buffer_fdpd[419];

    auto g_xyy_zz_y_xx = buffer_fdpd[420];

    auto g_xyy_zz_y_xy = buffer_fdpd[421];

    auto g_xyy_zz_y_xz = buffer_fdpd[422];

    auto g_xyy_zz_y_yy = buffer_fdpd[423];

    auto g_xyy_zz_y_yz = buffer_fdpd[424];

    auto g_xyy_zz_y_zz = buffer_fdpd[425];

    auto g_xyy_zz_z_xx = buffer_fdpd[426];

    auto g_xyy_zz_z_xy = buffer_fdpd[427];

    auto g_xyy_zz_z_xz = buffer_fdpd[428];

    auto g_xyy_zz_z_yy = buffer_fdpd[429];

    auto g_xyy_zz_z_yz = buffer_fdpd[430];

    auto g_xyy_zz_z_zz = buffer_fdpd[431];

    auto g_xyz_xx_x_xx = buffer_fdpd[432];

    auto g_xyz_xx_x_xy = buffer_fdpd[433];

    auto g_xyz_xx_x_xz = buffer_fdpd[434];

    auto g_xyz_xx_x_yy = buffer_fdpd[435];

    auto g_xyz_xx_x_yz = buffer_fdpd[436];

    auto g_xyz_xx_x_zz = buffer_fdpd[437];

    auto g_xyz_xx_y_xx = buffer_fdpd[438];

    auto g_xyz_xx_y_xy = buffer_fdpd[439];

    auto g_xyz_xx_y_xz = buffer_fdpd[440];

    auto g_xyz_xx_y_yy = buffer_fdpd[441];

    auto g_xyz_xx_y_yz = buffer_fdpd[442];

    auto g_xyz_xx_y_zz = buffer_fdpd[443];

    auto g_xyz_xx_z_xx = buffer_fdpd[444];

    auto g_xyz_xx_z_xy = buffer_fdpd[445];

    auto g_xyz_xx_z_xz = buffer_fdpd[446];

    auto g_xyz_xx_z_yy = buffer_fdpd[447];

    auto g_xyz_xx_z_yz = buffer_fdpd[448];

    auto g_xyz_xx_z_zz = buffer_fdpd[449];

    auto g_xyz_xy_x_xx = buffer_fdpd[450];

    auto g_xyz_xy_x_xy = buffer_fdpd[451];

    auto g_xyz_xy_x_xz = buffer_fdpd[452];

    auto g_xyz_xy_x_yy = buffer_fdpd[453];

    auto g_xyz_xy_x_yz = buffer_fdpd[454];

    auto g_xyz_xy_x_zz = buffer_fdpd[455];

    auto g_xyz_xy_y_xx = buffer_fdpd[456];

    auto g_xyz_xy_y_xy = buffer_fdpd[457];

    auto g_xyz_xy_y_xz = buffer_fdpd[458];

    auto g_xyz_xy_y_yy = buffer_fdpd[459];

    auto g_xyz_xy_y_yz = buffer_fdpd[460];

    auto g_xyz_xy_y_zz = buffer_fdpd[461];

    auto g_xyz_xy_z_xx = buffer_fdpd[462];

    auto g_xyz_xy_z_xy = buffer_fdpd[463];

    auto g_xyz_xy_z_xz = buffer_fdpd[464];

    auto g_xyz_xy_z_yy = buffer_fdpd[465];

    auto g_xyz_xy_z_yz = buffer_fdpd[466];

    auto g_xyz_xy_z_zz = buffer_fdpd[467];

    auto g_xyz_xz_x_xx = buffer_fdpd[468];

    auto g_xyz_xz_x_xy = buffer_fdpd[469];

    auto g_xyz_xz_x_xz = buffer_fdpd[470];

    auto g_xyz_xz_x_yy = buffer_fdpd[471];

    auto g_xyz_xz_x_yz = buffer_fdpd[472];

    auto g_xyz_xz_x_zz = buffer_fdpd[473];

    auto g_xyz_xz_y_xx = buffer_fdpd[474];

    auto g_xyz_xz_y_xy = buffer_fdpd[475];

    auto g_xyz_xz_y_xz = buffer_fdpd[476];

    auto g_xyz_xz_y_yy = buffer_fdpd[477];

    auto g_xyz_xz_y_yz = buffer_fdpd[478];

    auto g_xyz_xz_y_zz = buffer_fdpd[479];

    auto g_xyz_xz_z_xx = buffer_fdpd[480];

    auto g_xyz_xz_z_xy = buffer_fdpd[481];

    auto g_xyz_xz_z_xz = buffer_fdpd[482];

    auto g_xyz_xz_z_yy = buffer_fdpd[483];

    auto g_xyz_xz_z_yz = buffer_fdpd[484];

    auto g_xyz_xz_z_zz = buffer_fdpd[485];

    auto g_xyz_yy_x_xx = buffer_fdpd[486];

    auto g_xyz_yy_x_xy = buffer_fdpd[487];

    auto g_xyz_yy_x_xz = buffer_fdpd[488];

    auto g_xyz_yy_x_yy = buffer_fdpd[489];

    auto g_xyz_yy_x_yz = buffer_fdpd[490];

    auto g_xyz_yy_x_zz = buffer_fdpd[491];

    auto g_xyz_yy_y_xx = buffer_fdpd[492];

    auto g_xyz_yy_y_xy = buffer_fdpd[493];

    auto g_xyz_yy_y_xz = buffer_fdpd[494];

    auto g_xyz_yy_y_yy = buffer_fdpd[495];

    auto g_xyz_yy_y_yz = buffer_fdpd[496];

    auto g_xyz_yy_y_zz = buffer_fdpd[497];

    auto g_xyz_yy_z_xx = buffer_fdpd[498];

    auto g_xyz_yy_z_xy = buffer_fdpd[499];

    auto g_xyz_yy_z_xz = buffer_fdpd[500];

    auto g_xyz_yy_z_yy = buffer_fdpd[501];

    auto g_xyz_yy_z_yz = buffer_fdpd[502];

    auto g_xyz_yy_z_zz = buffer_fdpd[503];

    auto g_xyz_yz_x_xx = buffer_fdpd[504];

    auto g_xyz_yz_x_xy = buffer_fdpd[505];

    auto g_xyz_yz_x_xz = buffer_fdpd[506];

    auto g_xyz_yz_x_yy = buffer_fdpd[507];

    auto g_xyz_yz_x_yz = buffer_fdpd[508];

    auto g_xyz_yz_x_zz = buffer_fdpd[509];

    auto g_xyz_yz_y_xx = buffer_fdpd[510];

    auto g_xyz_yz_y_xy = buffer_fdpd[511];

    auto g_xyz_yz_y_xz = buffer_fdpd[512];

    auto g_xyz_yz_y_yy = buffer_fdpd[513];

    auto g_xyz_yz_y_yz = buffer_fdpd[514];

    auto g_xyz_yz_y_zz = buffer_fdpd[515];

    auto g_xyz_yz_z_xx = buffer_fdpd[516];

    auto g_xyz_yz_z_xy = buffer_fdpd[517];

    auto g_xyz_yz_z_xz = buffer_fdpd[518];

    auto g_xyz_yz_z_yy = buffer_fdpd[519];

    auto g_xyz_yz_z_yz = buffer_fdpd[520];

    auto g_xyz_yz_z_zz = buffer_fdpd[521];

    auto g_xyz_zz_x_xx = buffer_fdpd[522];

    auto g_xyz_zz_x_xy = buffer_fdpd[523];

    auto g_xyz_zz_x_xz = buffer_fdpd[524];

    auto g_xyz_zz_x_yy = buffer_fdpd[525];

    auto g_xyz_zz_x_yz = buffer_fdpd[526];

    auto g_xyz_zz_x_zz = buffer_fdpd[527];

    auto g_xyz_zz_y_xx = buffer_fdpd[528];

    auto g_xyz_zz_y_xy = buffer_fdpd[529];

    auto g_xyz_zz_y_xz = buffer_fdpd[530];

    auto g_xyz_zz_y_yy = buffer_fdpd[531];

    auto g_xyz_zz_y_yz = buffer_fdpd[532];

    auto g_xyz_zz_y_zz = buffer_fdpd[533];

    auto g_xyz_zz_z_xx = buffer_fdpd[534];

    auto g_xyz_zz_z_xy = buffer_fdpd[535];

    auto g_xyz_zz_z_xz = buffer_fdpd[536];

    auto g_xyz_zz_z_yy = buffer_fdpd[537];

    auto g_xyz_zz_z_yz = buffer_fdpd[538];

    auto g_xyz_zz_z_zz = buffer_fdpd[539];

    auto g_xzz_xx_x_xx = buffer_fdpd[540];

    auto g_xzz_xx_x_xy = buffer_fdpd[541];

    auto g_xzz_xx_x_xz = buffer_fdpd[542];

    auto g_xzz_xx_x_yy = buffer_fdpd[543];

    auto g_xzz_xx_x_yz = buffer_fdpd[544];

    auto g_xzz_xx_x_zz = buffer_fdpd[545];

    auto g_xzz_xx_y_xx = buffer_fdpd[546];

    auto g_xzz_xx_y_xy = buffer_fdpd[547];

    auto g_xzz_xx_y_xz = buffer_fdpd[548];

    auto g_xzz_xx_y_yy = buffer_fdpd[549];

    auto g_xzz_xx_y_yz = buffer_fdpd[550];

    auto g_xzz_xx_y_zz = buffer_fdpd[551];

    auto g_xzz_xx_z_xx = buffer_fdpd[552];

    auto g_xzz_xx_z_xy = buffer_fdpd[553];

    auto g_xzz_xx_z_xz = buffer_fdpd[554];

    auto g_xzz_xx_z_yy = buffer_fdpd[555];

    auto g_xzz_xx_z_yz = buffer_fdpd[556];

    auto g_xzz_xx_z_zz = buffer_fdpd[557];

    auto g_xzz_xy_x_xx = buffer_fdpd[558];

    auto g_xzz_xy_x_xy = buffer_fdpd[559];

    auto g_xzz_xy_x_xz = buffer_fdpd[560];

    auto g_xzz_xy_x_yy = buffer_fdpd[561];

    auto g_xzz_xy_x_yz = buffer_fdpd[562];

    auto g_xzz_xy_x_zz = buffer_fdpd[563];

    auto g_xzz_xy_y_xx = buffer_fdpd[564];

    auto g_xzz_xy_y_xy = buffer_fdpd[565];

    auto g_xzz_xy_y_xz = buffer_fdpd[566];

    auto g_xzz_xy_y_yy = buffer_fdpd[567];

    auto g_xzz_xy_y_yz = buffer_fdpd[568];

    auto g_xzz_xy_y_zz = buffer_fdpd[569];

    auto g_xzz_xy_z_xx = buffer_fdpd[570];

    auto g_xzz_xy_z_xy = buffer_fdpd[571];

    auto g_xzz_xy_z_xz = buffer_fdpd[572];

    auto g_xzz_xy_z_yy = buffer_fdpd[573];

    auto g_xzz_xy_z_yz = buffer_fdpd[574];

    auto g_xzz_xy_z_zz = buffer_fdpd[575];

    auto g_xzz_xz_x_xx = buffer_fdpd[576];

    auto g_xzz_xz_x_xy = buffer_fdpd[577];

    auto g_xzz_xz_x_xz = buffer_fdpd[578];

    auto g_xzz_xz_x_yy = buffer_fdpd[579];

    auto g_xzz_xz_x_yz = buffer_fdpd[580];

    auto g_xzz_xz_x_zz = buffer_fdpd[581];

    auto g_xzz_xz_y_xx = buffer_fdpd[582];

    auto g_xzz_xz_y_xy = buffer_fdpd[583];

    auto g_xzz_xz_y_xz = buffer_fdpd[584];

    auto g_xzz_xz_y_yy = buffer_fdpd[585];

    auto g_xzz_xz_y_yz = buffer_fdpd[586];

    auto g_xzz_xz_y_zz = buffer_fdpd[587];

    auto g_xzz_xz_z_xx = buffer_fdpd[588];

    auto g_xzz_xz_z_xy = buffer_fdpd[589];

    auto g_xzz_xz_z_xz = buffer_fdpd[590];

    auto g_xzz_xz_z_yy = buffer_fdpd[591];

    auto g_xzz_xz_z_yz = buffer_fdpd[592];

    auto g_xzz_xz_z_zz = buffer_fdpd[593];

    auto g_xzz_yy_x_xx = buffer_fdpd[594];

    auto g_xzz_yy_x_xy = buffer_fdpd[595];

    auto g_xzz_yy_x_xz = buffer_fdpd[596];

    auto g_xzz_yy_x_yy = buffer_fdpd[597];

    auto g_xzz_yy_x_yz = buffer_fdpd[598];

    auto g_xzz_yy_x_zz = buffer_fdpd[599];

    auto g_xzz_yy_y_xx = buffer_fdpd[600];

    auto g_xzz_yy_y_xy = buffer_fdpd[601];

    auto g_xzz_yy_y_xz = buffer_fdpd[602];

    auto g_xzz_yy_y_yy = buffer_fdpd[603];

    auto g_xzz_yy_y_yz = buffer_fdpd[604];

    auto g_xzz_yy_y_zz = buffer_fdpd[605];

    auto g_xzz_yy_z_xx = buffer_fdpd[606];

    auto g_xzz_yy_z_xy = buffer_fdpd[607];

    auto g_xzz_yy_z_xz = buffer_fdpd[608];

    auto g_xzz_yy_z_yy = buffer_fdpd[609];

    auto g_xzz_yy_z_yz = buffer_fdpd[610];

    auto g_xzz_yy_z_zz = buffer_fdpd[611];

    auto g_xzz_yz_x_xx = buffer_fdpd[612];

    auto g_xzz_yz_x_xy = buffer_fdpd[613];

    auto g_xzz_yz_x_xz = buffer_fdpd[614];

    auto g_xzz_yz_x_yy = buffer_fdpd[615];

    auto g_xzz_yz_x_yz = buffer_fdpd[616];

    auto g_xzz_yz_x_zz = buffer_fdpd[617];

    auto g_xzz_yz_y_xx = buffer_fdpd[618];

    auto g_xzz_yz_y_xy = buffer_fdpd[619];

    auto g_xzz_yz_y_xz = buffer_fdpd[620];

    auto g_xzz_yz_y_yy = buffer_fdpd[621];

    auto g_xzz_yz_y_yz = buffer_fdpd[622];

    auto g_xzz_yz_y_zz = buffer_fdpd[623];

    auto g_xzz_yz_z_xx = buffer_fdpd[624];

    auto g_xzz_yz_z_xy = buffer_fdpd[625];

    auto g_xzz_yz_z_xz = buffer_fdpd[626];

    auto g_xzz_yz_z_yy = buffer_fdpd[627];

    auto g_xzz_yz_z_yz = buffer_fdpd[628];

    auto g_xzz_yz_z_zz = buffer_fdpd[629];

    auto g_xzz_zz_x_xx = buffer_fdpd[630];

    auto g_xzz_zz_x_xy = buffer_fdpd[631];

    auto g_xzz_zz_x_xz = buffer_fdpd[632];

    auto g_xzz_zz_x_yy = buffer_fdpd[633];

    auto g_xzz_zz_x_yz = buffer_fdpd[634];

    auto g_xzz_zz_x_zz = buffer_fdpd[635];

    auto g_xzz_zz_y_xx = buffer_fdpd[636];

    auto g_xzz_zz_y_xy = buffer_fdpd[637];

    auto g_xzz_zz_y_xz = buffer_fdpd[638];

    auto g_xzz_zz_y_yy = buffer_fdpd[639];

    auto g_xzz_zz_y_yz = buffer_fdpd[640];

    auto g_xzz_zz_y_zz = buffer_fdpd[641];

    auto g_xzz_zz_z_xx = buffer_fdpd[642];

    auto g_xzz_zz_z_xy = buffer_fdpd[643];

    auto g_xzz_zz_z_xz = buffer_fdpd[644];

    auto g_xzz_zz_z_yy = buffer_fdpd[645];

    auto g_xzz_zz_z_yz = buffer_fdpd[646];

    auto g_xzz_zz_z_zz = buffer_fdpd[647];

    auto g_yyy_xx_x_xx = buffer_fdpd[648];

    auto g_yyy_xx_x_xy = buffer_fdpd[649];

    auto g_yyy_xx_x_xz = buffer_fdpd[650];

    auto g_yyy_xx_x_yy = buffer_fdpd[651];

    auto g_yyy_xx_x_yz = buffer_fdpd[652];

    auto g_yyy_xx_x_zz = buffer_fdpd[653];

    auto g_yyy_xx_y_xx = buffer_fdpd[654];

    auto g_yyy_xx_y_xy = buffer_fdpd[655];

    auto g_yyy_xx_y_xz = buffer_fdpd[656];

    auto g_yyy_xx_y_yy = buffer_fdpd[657];

    auto g_yyy_xx_y_yz = buffer_fdpd[658];

    auto g_yyy_xx_y_zz = buffer_fdpd[659];

    auto g_yyy_xx_z_xx = buffer_fdpd[660];

    auto g_yyy_xx_z_xy = buffer_fdpd[661];

    auto g_yyy_xx_z_xz = buffer_fdpd[662];

    auto g_yyy_xx_z_yy = buffer_fdpd[663];

    auto g_yyy_xx_z_yz = buffer_fdpd[664];

    auto g_yyy_xx_z_zz = buffer_fdpd[665];

    auto g_yyy_xy_x_xx = buffer_fdpd[666];

    auto g_yyy_xy_x_xy = buffer_fdpd[667];

    auto g_yyy_xy_x_xz = buffer_fdpd[668];

    auto g_yyy_xy_x_yy = buffer_fdpd[669];

    auto g_yyy_xy_x_yz = buffer_fdpd[670];

    auto g_yyy_xy_x_zz = buffer_fdpd[671];

    auto g_yyy_xy_y_xx = buffer_fdpd[672];

    auto g_yyy_xy_y_xy = buffer_fdpd[673];

    auto g_yyy_xy_y_xz = buffer_fdpd[674];

    auto g_yyy_xy_y_yy = buffer_fdpd[675];

    auto g_yyy_xy_y_yz = buffer_fdpd[676];

    auto g_yyy_xy_y_zz = buffer_fdpd[677];

    auto g_yyy_xy_z_xx = buffer_fdpd[678];

    auto g_yyy_xy_z_xy = buffer_fdpd[679];

    auto g_yyy_xy_z_xz = buffer_fdpd[680];

    auto g_yyy_xy_z_yy = buffer_fdpd[681];

    auto g_yyy_xy_z_yz = buffer_fdpd[682];

    auto g_yyy_xy_z_zz = buffer_fdpd[683];

    auto g_yyy_xz_x_xx = buffer_fdpd[684];

    auto g_yyy_xz_x_xy = buffer_fdpd[685];

    auto g_yyy_xz_x_xz = buffer_fdpd[686];

    auto g_yyy_xz_x_yy = buffer_fdpd[687];

    auto g_yyy_xz_x_yz = buffer_fdpd[688];

    auto g_yyy_xz_x_zz = buffer_fdpd[689];

    auto g_yyy_xz_y_xx = buffer_fdpd[690];

    auto g_yyy_xz_y_xy = buffer_fdpd[691];

    auto g_yyy_xz_y_xz = buffer_fdpd[692];

    auto g_yyy_xz_y_yy = buffer_fdpd[693];

    auto g_yyy_xz_y_yz = buffer_fdpd[694];

    auto g_yyy_xz_y_zz = buffer_fdpd[695];

    auto g_yyy_xz_z_xx = buffer_fdpd[696];

    auto g_yyy_xz_z_xy = buffer_fdpd[697];

    auto g_yyy_xz_z_xz = buffer_fdpd[698];

    auto g_yyy_xz_z_yy = buffer_fdpd[699];

    auto g_yyy_xz_z_yz = buffer_fdpd[700];

    auto g_yyy_xz_z_zz = buffer_fdpd[701];

    auto g_yyy_yy_x_xx = buffer_fdpd[702];

    auto g_yyy_yy_x_xy = buffer_fdpd[703];

    auto g_yyy_yy_x_xz = buffer_fdpd[704];

    auto g_yyy_yy_x_yy = buffer_fdpd[705];

    auto g_yyy_yy_x_yz = buffer_fdpd[706];

    auto g_yyy_yy_x_zz = buffer_fdpd[707];

    auto g_yyy_yy_y_xx = buffer_fdpd[708];

    auto g_yyy_yy_y_xy = buffer_fdpd[709];

    auto g_yyy_yy_y_xz = buffer_fdpd[710];

    auto g_yyy_yy_y_yy = buffer_fdpd[711];

    auto g_yyy_yy_y_yz = buffer_fdpd[712];

    auto g_yyy_yy_y_zz = buffer_fdpd[713];

    auto g_yyy_yy_z_xx = buffer_fdpd[714];

    auto g_yyy_yy_z_xy = buffer_fdpd[715];

    auto g_yyy_yy_z_xz = buffer_fdpd[716];

    auto g_yyy_yy_z_yy = buffer_fdpd[717];

    auto g_yyy_yy_z_yz = buffer_fdpd[718];

    auto g_yyy_yy_z_zz = buffer_fdpd[719];

    auto g_yyy_yz_x_xx = buffer_fdpd[720];

    auto g_yyy_yz_x_xy = buffer_fdpd[721];

    auto g_yyy_yz_x_xz = buffer_fdpd[722];

    auto g_yyy_yz_x_yy = buffer_fdpd[723];

    auto g_yyy_yz_x_yz = buffer_fdpd[724];

    auto g_yyy_yz_x_zz = buffer_fdpd[725];

    auto g_yyy_yz_y_xx = buffer_fdpd[726];

    auto g_yyy_yz_y_xy = buffer_fdpd[727];

    auto g_yyy_yz_y_xz = buffer_fdpd[728];

    auto g_yyy_yz_y_yy = buffer_fdpd[729];

    auto g_yyy_yz_y_yz = buffer_fdpd[730];

    auto g_yyy_yz_y_zz = buffer_fdpd[731];

    auto g_yyy_yz_z_xx = buffer_fdpd[732];

    auto g_yyy_yz_z_xy = buffer_fdpd[733];

    auto g_yyy_yz_z_xz = buffer_fdpd[734];

    auto g_yyy_yz_z_yy = buffer_fdpd[735];

    auto g_yyy_yz_z_yz = buffer_fdpd[736];

    auto g_yyy_yz_z_zz = buffer_fdpd[737];

    auto g_yyy_zz_x_xx = buffer_fdpd[738];

    auto g_yyy_zz_x_xy = buffer_fdpd[739];

    auto g_yyy_zz_x_xz = buffer_fdpd[740];

    auto g_yyy_zz_x_yy = buffer_fdpd[741];

    auto g_yyy_zz_x_yz = buffer_fdpd[742];

    auto g_yyy_zz_x_zz = buffer_fdpd[743];

    auto g_yyy_zz_y_xx = buffer_fdpd[744];

    auto g_yyy_zz_y_xy = buffer_fdpd[745];

    auto g_yyy_zz_y_xz = buffer_fdpd[746];

    auto g_yyy_zz_y_yy = buffer_fdpd[747];

    auto g_yyy_zz_y_yz = buffer_fdpd[748];

    auto g_yyy_zz_y_zz = buffer_fdpd[749];

    auto g_yyy_zz_z_xx = buffer_fdpd[750];

    auto g_yyy_zz_z_xy = buffer_fdpd[751];

    auto g_yyy_zz_z_xz = buffer_fdpd[752];

    auto g_yyy_zz_z_yy = buffer_fdpd[753];

    auto g_yyy_zz_z_yz = buffer_fdpd[754];

    auto g_yyy_zz_z_zz = buffer_fdpd[755];

    auto g_yyz_xx_x_xx = buffer_fdpd[756];

    auto g_yyz_xx_x_xy = buffer_fdpd[757];

    auto g_yyz_xx_x_xz = buffer_fdpd[758];

    auto g_yyz_xx_x_yy = buffer_fdpd[759];

    auto g_yyz_xx_x_yz = buffer_fdpd[760];

    auto g_yyz_xx_x_zz = buffer_fdpd[761];

    auto g_yyz_xx_y_xx = buffer_fdpd[762];

    auto g_yyz_xx_y_xy = buffer_fdpd[763];

    auto g_yyz_xx_y_xz = buffer_fdpd[764];

    auto g_yyz_xx_y_yy = buffer_fdpd[765];

    auto g_yyz_xx_y_yz = buffer_fdpd[766];

    auto g_yyz_xx_y_zz = buffer_fdpd[767];

    auto g_yyz_xx_z_xx = buffer_fdpd[768];

    auto g_yyz_xx_z_xy = buffer_fdpd[769];

    auto g_yyz_xx_z_xz = buffer_fdpd[770];

    auto g_yyz_xx_z_yy = buffer_fdpd[771];

    auto g_yyz_xx_z_yz = buffer_fdpd[772];

    auto g_yyz_xx_z_zz = buffer_fdpd[773];

    auto g_yyz_xy_x_xx = buffer_fdpd[774];

    auto g_yyz_xy_x_xy = buffer_fdpd[775];

    auto g_yyz_xy_x_xz = buffer_fdpd[776];

    auto g_yyz_xy_x_yy = buffer_fdpd[777];

    auto g_yyz_xy_x_yz = buffer_fdpd[778];

    auto g_yyz_xy_x_zz = buffer_fdpd[779];

    auto g_yyz_xy_y_xx = buffer_fdpd[780];

    auto g_yyz_xy_y_xy = buffer_fdpd[781];

    auto g_yyz_xy_y_xz = buffer_fdpd[782];

    auto g_yyz_xy_y_yy = buffer_fdpd[783];

    auto g_yyz_xy_y_yz = buffer_fdpd[784];

    auto g_yyz_xy_y_zz = buffer_fdpd[785];

    auto g_yyz_xy_z_xx = buffer_fdpd[786];

    auto g_yyz_xy_z_xy = buffer_fdpd[787];

    auto g_yyz_xy_z_xz = buffer_fdpd[788];

    auto g_yyz_xy_z_yy = buffer_fdpd[789];

    auto g_yyz_xy_z_yz = buffer_fdpd[790];

    auto g_yyz_xy_z_zz = buffer_fdpd[791];

    auto g_yyz_xz_x_xx = buffer_fdpd[792];

    auto g_yyz_xz_x_xy = buffer_fdpd[793];

    auto g_yyz_xz_x_xz = buffer_fdpd[794];

    auto g_yyz_xz_x_yy = buffer_fdpd[795];

    auto g_yyz_xz_x_yz = buffer_fdpd[796];

    auto g_yyz_xz_x_zz = buffer_fdpd[797];

    auto g_yyz_xz_y_xx = buffer_fdpd[798];

    auto g_yyz_xz_y_xy = buffer_fdpd[799];

    auto g_yyz_xz_y_xz = buffer_fdpd[800];

    auto g_yyz_xz_y_yy = buffer_fdpd[801];

    auto g_yyz_xz_y_yz = buffer_fdpd[802];

    auto g_yyz_xz_y_zz = buffer_fdpd[803];

    auto g_yyz_xz_z_xx = buffer_fdpd[804];

    auto g_yyz_xz_z_xy = buffer_fdpd[805];

    auto g_yyz_xz_z_xz = buffer_fdpd[806];

    auto g_yyz_xz_z_yy = buffer_fdpd[807];

    auto g_yyz_xz_z_yz = buffer_fdpd[808];

    auto g_yyz_xz_z_zz = buffer_fdpd[809];

    auto g_yyz_yy_x_xx = buffer_fdpd[810];

    auto g_yyz_yy_x_xy = buffer_fdpd[811];

    auto g_yyz_yy_x_xz = buffer_fdpd[812];

    auto g_yyz_yy_x_yy = buffer_fdpd[813];

    auto g_yyz_yy_x_yz = buffer_fdpd[814];

    auto g_yyz_yy_x_zz = buffer_fdpd[815];

    auto g_yyz_yy_y_xx = buffer_fdpd[816];

    auto g_yyz_yy_y_xy = buffer_fdpd[817];

    auto g_yyz_yy_y_xz = buffer_fdpd[818];

    auto g_yyz_yy_y_yy = buffer_fdpd[819];

    auto g_yyz_yy_y_yz = buffer_fdpd[820];

    auto g_yyz_yy_y_zz = buffer_fdpd[821];

    auto g_yyz_yy_z_xx = buffer_fdpd[822];

    auto g_yyz_yy_z_xy = buffer_fdpd[823];

    auto g_yyz_yy_z_xz = buffer_fdpd[824];

    auto g_yyz_yy_z_yy = buffer_fdpd[825];

    auto g_yyz_yy_z_yz = buffer_fdpd[826];

    auto g_yyz_yy_z_zz = buffer_fdpd[827];

    auto g_yyz_yz_x_xx = buffer_fdpd[828];

    auto g_yyz_yz_x_xy = buffer_fdpd[829];

    auto g_yyz_yz_x_xz = buffer_fdpd[830];

    auto g_yyz_yz_x_yy = buffer_fdpd[831];

    auto g_yyz_yz_x_yz = buffer_fdpd[832];

    auto g_yyz_yz_x_zz = buffer_fdpd[833];

    auto g_yyz_yz_y_xx = buffer_fdpd[834];

    auto g_yyz_yz_y_xy = buffer_fdpd[835];

    auto g_yyz_yz_y_xz = buffer_fdpd[836];

    auto g_yyz_yz_y_yy = buffer_fdpd[837];

    auto g_yyz_yz_y_yz = buffer_fdpd[838];

    auto g_yyz_yz_y_zz = buffer_fdpd[839];

    auto g_yyz_yz_z_xx = buffer_fdpd[840];

    auto g_yyz_yz_z_xy = buffer_fdpd[841];

    auto g_yyz_yz_z_xz = buffer_fdpd[842];

    auto g_yyz_yz_z_yy = buffer_fdpd[843];

    auto g_yyz_yz_z_yz = buffer_fdpd[844];

    auto g_yyz_yz_z_zz = buffer_fdpd[845];

    auto g_yyz_zz_x_xx = buffer_fdpd[846];

    auto g_yyz_zz_x_xy = buffer_fdpd[847];

    auto g_yyz_zz_x_xz = buffer_fdpd[848];

    auto g_yyz_zz_x_yy = buffer_fdpd[849];

    auto g_yyz_zz_x_yz = buffer_fdpd[850];

    auto g_yyz_zz_x_zz = buffer_fdpd[851];

    auto g_yyz_zz_y_xx = buffer_fdpd[852];

    auto g_yyz_zz_y_xy = buffer_fdpd[853];

    auto g_yyz_zz_y_xz = buffer_fdpd[854];

    auto g_yyz_zz_y_yy = buffer_fdpd[855];

    auto g_yyz_zz_y_yz = buffer_fdpd[856];

    auto g_yyz_zz_y_zz = buffer_fdpd[857];

    auto g_yyz_zz_z_xx = buffer_fdpd[858];

    auto g_yyz_zz_z_xy = buffer_fdpd[859];

    auto g_yyz_zz_z_xz = buffer_fdpd[860];

    auto g_yyz_zz_z_yy = buffer_fdpd[861];

    auto g_yyz_zz_z_yz = buffer_fdpd[862];

    auto g_yyz_zz_z_zz = buffer_fdpd[863];

    auto g_yzz_xx_x_xx = buffer_fdpd[864];

    auto g_yzz_xx_x_xy = buffer_fdpd[865];

    auto g_yzz_xx_x_xz = buffer_fdpd[866];

    auto g_yzz_xx_x_yy = buffer_fdpd[867];

    auto g_yzz_xx_x_yz = buffer_fdpd[868];

    auto g_yzz_xx_x_zz = buffer_fdpd[869];

    auto g_yzz_xx_y_xx = buffer_fdpd[870];

    auto g_yzz_xx_y_xy = buffer_fdpd[871];

    auto g_yzz_xx_y_xz = buffer_fdpd[872];

    auto g_yzz_xx_y_yy = buffer_fdpd[873];

    auto g_yzz_xx_y_yz = buffer_fdpd[874];

    auto g_yzz_xx_y_zz = buffer_fdpd[875];

    auto g_yzz_xx_z_xx = buffer_fdpd[876];

    auto g_yzz_xx_z_xy = buffer_fdpd[877];

    auto g_yzz_xx_z_xz = buffer_fdpd[878];

    auto g_yzz_xx_z_yy = buffer_fdpd[879];

    auto g_yzz_xx_z_yz = buffer_fdpd[880];

    auto g_yzz_xx_z_zz = buffer_fdpd[881];

    auto g_yzz_xy_x_xx = buffer_fdpd[882];

    auto g_yzz_xy_x_xy = buffer_fdpd[883];

    auto g_yzz_xy_x_xz = buffer_fdpd[884];

    auto g_yzz_xy_x_yy = buffer_fdpd[885];

    auto g_yzz_xy_x_yz = buffer_fdpd[886];

    auto g_yzz_xy_x_zz = buffer_fdpd[887];

    auto g_yzz_xy_y_xx = buffer_fdpd[888];

    auto g_yzz_xy_y_xy = buffer_fdpd[889];

    auto g_yzz_xy_y_xz = buffer_fdpd[890];

    auto g_yzz_xy_y_yy = buffer_fdpd[891];

    auto g_yzz_xy_y_yz = buffer_fdpd[892];

    auto g_yzz_xy_y_zz = buffer_fdpd[893];

    auto g_yzz_xy_z_xx = buffer_fdpd[894];

    auto g_yzz_xy_z_xy = buffer_fdpd[895];

    auto g_yzz_xy_z_xz = buffer_fdpd[896];

    auto g_yzz_xy_z_yy = buffer_fdpd[897];

    auto g_yzz_xy_z_yz = buffer_fdpd[898];

    auto g_yzz_xy_z_zz = buffer_fdpd[899];

    auto g_yzz_xz_x_xx = buffer_fdpd[900];

    auto g_yzz_xz_x_xy = buffer_fdpd[901];

    auto g_yzz_xz_x_xz = buffer_fdpd[902];

    auto g_yzz_xz_x_yy = buffer_fdpd[903];

    auto g_yzz_xz_x_yz = buffer_fdpd[904];

    auto g_yzz_xz_x_zz = buffer_fdpd[905];

    auto g_yzz_xz_y_xx = buffer_fdpd[906];

    auto g_yzz_xz_y_xy = buffer_fdpd[907];

    auto g_yzz_xz_y_xz = buffer_fdpd[908];

    auto g_yzz_xz_y_yy = buffer_fdpd[909];

    auto g_yzz_xz_y_yz = buffer_fdpd[910];

    auto g_yzz_xz_y_zz = buffer_fdpd[911];

    auto g_yzz_xz_z_xx = buffer_fdpd[912];

    auto g_yzz_xz_z_xy = buffer_fdpd[913];

    auto g_yzz_xz_z_xz = buffer_fdpd[914];

    auto g_yzz_xz_z_yy = buffer_fdpd[915];

    auto g_yzz_xz_z_yz = buffer_fdpd[916];

    auto g_yzz_xz_z_zz = buffer_fdpd[917];

    auto g_yzz_yy_x_xx = buffer_fdpd[918];

    auto g_yzz_yy_x_xy = buffer_fdpd[919];

    auto g_yzz_yy_x_xz = buffer_fdpd[920];

    auto g_yzz_yy_x_yy = buffer_fdpd[921];

    auto g_yzz_yy_x_yz = buffer_fdpd[922];

    auto g_yzz_yy_x_zz = buffer_fdpd[923];

    auto g_yzz_yy_y_xx = buffer_fdpd[924];

    auto g_yzz_yy_y_xy = buffer_fdpd[925];

    auto g_yzz_yy_y_xz = buffer_fdpd[926];

    auto g_yzz_yy_y_yy = buffer_fdpd[927];

    auto g_yzz_yy_y_yz = buffer_fdpd[928];

    auto g_yzz_yy_y_zz = buffer_fdpd[929];

    auto g_yzz_yy_z_xx = buffer_fdpd[930];

    auto g_yzz_yy_z_xy = buffer_fdpd[931];

    auto g_yzz_yy_z_xz = buffer_fdpd[932];

    auto g_yzz_yy_z_yy = buffer_fdpd[933];

    auto g_yzz_yy_z_yz = buffer_fdpd[934];

    auto g_yzz_yy_z_zz = buffer_fdpd[935];

    auto g_yzz_yz_x_xx = buffer_fdpd[936];

    auto g_yzz_yz_x_xy = buffer_fdpd[937];

    auto g_yzz_yz_x_xz = buffer_fdpd[938];

    auto g_yzz_yz_x_yy = buffer_fdpd[939];

    auto g_yzz_yz_x_yz = buffer_fdpd[940];

    auto g_yzz_yz_x_zz = buffer_fdpd[941];

    auto g_yzz_yz_y_xx = buffer_fdpd[942];

    auto g_yzz_yz_y_xy = buffer_fdpd[943];

    auto g_yzz_yz_y_xz = buffer_fdpd[944];

    auto g_yzz_yz_y_yy = buffer_fdpd[945];

    auto g_yzz_yz_y_yz = buffer_fdpd[946];

    auto g_yzz_yz_y_zz = buffer_fdpd[947];

    auto g_yzz_yz_z_xx = buffer_fdpd[948];

    auto g_yzz_yz_z_xy = buffer_fdpd[949];

    auto g_yzz_yz_z_xz = buffer_fdpd[950];

    auto g_yzz_yz_z_yy = buffer_fdpd[951];

    auto g_yzz_yz_z_yz = buffer_fdpd[952];

    auto g_yzz_yz_z_zz = buffer_fdpd[953];

    auto g_yzz_zz_x_xx = buffer_fdpd[954];

    auto g_yzz_zz_x_xy = buffer_fdpd[955];

    auto g_yzz_zz_x_xz = buffer_fdpd[956];

    auto g_yzz_zz_x_yy = buffer_fdpd[957];

    auto g_yzz_zz_x_yz = buffer_fdpd[958];

    auto g_yzz_zz_x_zz = buffer_fdpd[959];

    auto g_yzz_zz_y_xx = buffer_fdpd[960];

    auto g_yzz_zz_y_xy = buffer_fdpd[961];

    auto g_yzz_zz_y_xz = buffer_fdpd[962];

    auto g_yzz_zz_y_yy = buffer_fdpd[963];

    auto g_yzz_zz_y_yz = buffer_fdpd[964];

    auto g_yzz_zz_y_zz = buffer_fdpd[965];

    auto g_yzz_zz_z_xx = buffer_fdpd[966];

    auto g_yzz_zz_z_xy = buffer_fdpd[967];

    auto g_yzz_zz_z_xz = buffer_fdpd[968];

    auto g_yzz_zz_z_yy = buffer_fdpd[969];

    auto g_yzz_zz_z_yz = buffer_fdpd[970];

    auto g_yzz_zz_z_zz = buffer_fdpd[971];

    auto g_zzz_xx_x_xx = buffer_fdpd[972];

    auto g_zzz_xx_x_xy = buffer_fdpd[973];

    auto g_zzz_xx_x_xz = buffer_fdpd[974];

    auto g_zzz_xx_x_yy = buffer_fdpd[975];

    auto g_zzz_xx_x_yz = buffer_fdpd[976];

    auto g_zzz_xx_x_zz = buffer_fdpd[977];

    auto g_zzz_xx_y_xx = buffer_fdpd[978];

    auto g_zzz_xx_y_xy = buffer_fdpd[979];

    auto g_zzz_xx_y_xz = buffer_fdpd[980];

    auto g_zzz_xx_y_yy = buffer_fdpd[981];

    auto g_zzz_xx_y_yz = buffer_fdpd[982];

    auto g_zzz_xx_y_zz = buffer_fdpd[983];

    auto g_zzz_xx_z_xx = buffer_fdpd[984];

    auto g_zzz_xx_z_xy = buffer_fdpd[985];

    auto g_zzz_xx_z_xz = buffer_fdpd[986];

    auto g_zzz_xx_z_yy = buffer_fdpd[987];

    auto g_zzz_xx_z_yz = buffer_fdpd[988];

    auto g_zzz_xx_z_zz = buffer_fdpd[989];

    auto g_zzz_xy_x_xx = buffer_fdpd[990];

    auto g_zzz_xy_x_xy = buffer_fdpd[991];

    auto g_zzz_xy_x_xz = buffer_fdpd[992];

    auto g_zzz_xy_x_yy = buffer_fdpd[993];

    auto g_zzz_xy_x_yz = buffer_fdpd[994];

    auto g_zzz_xy_x_zz = buffer_fdpd[995];

    auto g_zzz_xy_y_xx = buffer_fdpd[996];

    auto g_zzz_xy_y_xy = buffer_fdpd[997];

    auto g_zzz_xy_y_xz = buffer_fdpd[998];

    auto g_zzz_xy_y_yy = buffer_fdpd[999];

    auto g_zzz_xy_y_yz = buffer_fdpd[1000];

    auto g_zzz_xy_y_zz = buffer_fdpd[1001];

    auto g_zzz_xy_z_xx = buffer_fdpd[1002];

    auto g_zzz_xy_z_xy = buffer_fdpd[1003];

    auto g_zzz_xy_z_xz = buffer_fdpd[1004];

    auto g_zzz_xy_z_yy = buffer_fdpd[1005];

    auto g_zzz_xy_z_yz = buffer_fdpd[1006];

    auto g_zzz_xy_z_zz = buffer_fdpd[1007];

    auto g_zzz_xz_x_xx = buffer_fdpd[1008];

    auto g_zzz_xz_x_xy = buffer_fdpd[1009];

    auto g_zzz_xz_x_xz = buffer_fdpd[1010];

    auto g_zzz_xz_x_yy = buffer_fdpd[1011];

    auto g_zzz_xz_x_yz = buffer_fdpd[1012];

    auto g_zzz_xz_x_zz = buffer_fdpd[1013];

    auto g_zzz_xz_y_xx = buffer_fdpd[1014];

    auto g_zzz_xz_y_xy = buffer_fdpd[1015];

    auto g_zzz_xz_y_xz = buffer_fdpd[1016];

    auto g_zzz_xz_y_yy = buffer_fdpd[1017];

    auto g_zzz_xz_y_yz = buffer_fdpd[1018];

    auto g_zzz_xz_y_zz = buffer_fdpd[1019];

    auto g_zzz_xz_z_xx = buffer_fdpd[1020];

    auto g_zzz_xz_z_xy = buffer_fdpd[1021];

    auto g_zzz_xz_z_xz = buffer_fdpd[1022];

    auto g_zzz_xz_z_yy = buffer_fdpd[1023];

    auto g_zzz_xz_z_yz = buffer_fdpd[1024];

    auto g_zzz_xz_z_zz = buffer_fdpd[1025];

    auto g_zzz_yy_x_xx = buffer_fdpd[1026];

    auto g_zzz_yy_x_xy = buffer_fdpd[1027];

    auto g_zzz_yy_x_xz = buffer_fdpd[1028];

    auto g_zzz_yy_x_yy = buffer_fdpd[1029];

    auto g_zzz_yy_x_yz = buffer_fdpd[1030];

    auto g_zzz_yy_x_zz = buffer_fdpd[1031];

    auto g_zzz_yy_y_xx = buffer_fdpd[1032];

    auto g_zzz_yy_y_xy = buffer_fdpd[1033];

    auto g_zzz_yy_y_xz = buffer_fdpd[1034];

    auto g_zzz_yy_y_yy = buffer_fdpd[1035];

    auto g_zzz_yy_y_yz = buffer_fdpd[1036];

    auto g_zzz_yy_y_zz = buffer_fdpd[1037];

    auto g_zzz_yy_z_xx = buffer_fdpd[1038];

    auto g_zzz_yy_z_xy = buffer_fdpd[1039];

    auto g_zzz_yy_z_xz = buffer_fdpd[1040];

    auto g_zzz_yy_z_yy = buffer_fdpd[1041];

    auto g_zzz_yy_z_yz = buffer_fdpd[1042];

    auto g_zzz_yy_z_zz = buffer_fdpd[1043];

    auto g_zzz_yz_x_xx = buffer_fdpd[1044];

    auto g_zzz_yz_x_xy = buffer_fdpd[1045];

    auto g_zzz_yz_x_xz = buffer_fdpd[1046];

    auto g_zzz_yz_x_yy = buffer_fdpd[1047];

    auto g_zzz_yz_x_yz = buffer_fdpd[1048];

    auto g_zzz_yz_x_zz = buffer_fdpd[1049];

    auto g_zzz_yz_y_xx = buffer_fdpd[1050];

    auto g_zzz_yz_y_xy = buffer_fdpd[1051];

    auto g_zzz_yz_y_xz = buffer_fdpd[1052];

    auto g_zzz_yz_y_yy = buffer_fdpd[1053];

    auto g_zzz_yz_y_yz = buffer_fdpd[1054];

    auto g_zzz_yz_y_zz = buffer_fdpd[1055];

    auto g_zzz_yz_z_xx = buffer_fdpd[1056];

    auto g_zzz_yz_z_xy = buffer_fdpd[1057];

    auto g_zzz_yz_z_xz = buffer_fdpd[1058];

    auto g_zzz_yz_z_yy = buffer_fdpd[1059];

    auto g_zzz_yz_z_yz = buffer_fdpd[1060];

    auto g_zzz_yz_z_zz = buffer_fdpd[1061];

    auto g_zzz_zz_x_xx = buffer_fdpd[1062];

    auto g_zzz_zz_x_xy = buffer_fdpd[1063];

    auto g_zzz_zz_x_xz = buffer_fdpd[1064];

    auto g_zzz_zz_x_yy = buffer_fdpd[1065];

    auto g_zzz_zz_x_yz = buffer_fdpd[1066];

    auto g_zzz_zz_x_zz = buffer_fdpd[1067];

    auto g_zzz_zz_y_xx = buffer_fdpd[1068];

    auto g_zzz_zz_y_xy = buffer_fdpd[1069];

    auto g_zzz_zz_y_xz = buffer_fdpd[1070];

    auto g_zzz_zz_y_yy = buffer_fdpd[1071];

    auto g_zzz_zz_y_yz = buffer_fdpd[1072];

    auto g_zzz_zz_y_zz = buffer_fdpd[1073];

    auto g_zzz_zz_z_xx = buffer_fdpd[1074];

    auto g_zzz_zz_z_xy = buffer_fdpd[1075];

    auto g_zzz_zz_z_xz = buffer_fdpd[1076];

    auto g_zzz_zz_z_yy = buffer_fdpd[1077];

    auto g_zzz_zz_z_yz = buffer_fdpd[1078];

    auto g_zzz_zz_z_zz = buffer_fdpd[1079];

    /// Set up components of integrals buffer : buffer_1010_ddsd

    auto g_x_0_x_0_xx_xx_0_xx = buffer_1010_ddsd[0];

    auto g_x_0_x_0_xx_xx_0_xy = buffer_1010_ddsd[1];

    auto g_x_0_x_0_xx_xx_0_xz = buffer_1010_ddsd[2];

    auto g_x_0_x_0_xx_xx_0_yy = buffer_1010_ddsd[3];

    auto g_x_0_x_0_xx_xx_0_yz = buffer_1010_ddsd[4];

    auto g_x_0_x_0_xx_xx_0_zz = buffer_1010_ddsd[5];

    auto g_x_0_x_0_xx_xy_0_xx = buffer_1010_ddsd[6];

    auto g_x_0_x_0_xx_xy_0_xy = buffer_1010_ddsd[7];

    auto g_x_0_x_0_xx_xy_0_xz = buffer_1010_ddsd[8];

    auto g_x_0_x_0_xx_xy_0_yy = buffer_1010_ddsd[9];

    auto g_x_0_x_0_xx_xy_0_yz = buffer_1010_ddsd[10];

    auto g_x_0_x_0_xx_xy_0_zz = buffer_1010_ddsd[11];

    auto g_x_0_x_0_xx_xz_0_xx = buffer_1010_ddsd[12];

    auto g_x_0_x_0_xx_xz_0_xy = buffer_1010_ddsd[13];

    auto g_x_0_x_0_xx_xz_0_xz = buffer_1010_ddsd[14];

    auto g_x_0_x_0_xx_xz_0_yy = buffer_1010_ddsd[15];

    auto g_x_0_x_0_xx_xz_0_yz = buffer_1010_ddsd[16];

    auto g_x_0_x_0_xx_xz_0_zz = buffer_1010_ddsd[17];

    auto g_x_0_x_0_xx_yy_0_xx = buffer_1010_ddsd[18];

    auto g_x_0_x_0_xx_yy_0_xy = buffer_1010_ddsd[19];

    auto g_x_0_x_0_xx_yy_0_xz = buffer_1010_ddsd[20];

    auto g_x_0_x_0_xx_yy_0_yy = buffer_1010_ddsd[21];

    auto g_x_0_x_0_xx_yy_0_yz = buffer_1010_ddsd[22];

    auto g_x_0_x_0_xx_yy_0_zz = buffer_1010_ddsd[23];

    auto g_x_0_x_0_xx_yz_0_xx = buffer_1010_ddsd[24];

    auto g_x_0_x_0_xx_yz_0_xy = buffer_1010_ddsd[25];

    auto g_x_0_x_0_xx_yz_0_xz = buffer_1010_ddsd[26];

    auto g_x_0_x_0_xx_yz_0_yy = buffer_1010_ddsd[27];

    auto g_x_0_x_0_xx_yz_0_yz = buffer_1010_ddsd[28];

    auto g_x_0_x_0_xx_yz_0_zz = buffer_1010_ddsd[29];

    auto g_x_0_x_0_xx_zz_0_xx = buffer_1010_ddsd[30];

    auto g_x_0_x_0_xx_zz_0_xy = buffer_1010_ddsd[31];

    auto g_x_0_x_0_xx_zz_0_xz = buffer_1010_ddsd[32];

    auto g_x_0_x_0_xx_zz_0_yy = buffer_1010_ddsd[33];

    auto g_x_0_x_0_xx_zz_0_yz = buffer_1010_ddsd[34];

    auto g_x_0_x_0_xx_zz_0_zz = buffer_1010_ddsd[35];

    auto g_x_0_x_0_xy_xx_0_xx = buffer_1010_ddsd[36];

    auto g_x_0_x_0_xy_xx_0_xy = buffer_1010_ddsd[37];

    auto g_x_0_x_0_xy_xx_0_xz = buffer_1010_ddsd[38];

    auto g_x_0_x_0_xy_xx_0_yy = buffer_1010_ddsd[39];

    auto g_x_0_x_0_xy_xx_0_yz = buffer_1010_ddsd[40];

    auto g_x_0_x_0_xy_xx_0_zz = buffer_1010_ddsd[41];

    auto g_x_0_x_0_xy_xy_0_xx = buffer_1010_ddsd[42];

    auto g_x_0_x_0_xy_xy_0_xy = buffer_1010_ddsd[43];

    auto g_x_0_x_0_xy_xy_0_xz = buffer_1010_ddsd[44];

    auto g_x_0_x_0_xy_xy_0_yy = buffer_1010_ddsd[45];

    auto g_x_0_x_0_xy_xy_0_yz = buffer_1010_ddsd[46];

    auto g_x_0_x_0_xy_xy_0_zz = buffer_1010_ddsd[47];

    auto g_x_0_x_0_xy_xz_0_xx = buffer_1010_ddsd[48];

    auto g_x_0_x_0_xy_xz_0_xy = buffer_1010_ddsd[49];

    auto g_x_0_x_0_xy_xz_0_xz = buffer_1010_ddsd[50];

    auto g_x_0_x_0_xy_xz_0_yy = buffer_1010_ddsd[51];

    auto g_x_0_x_0_xy_xz_0_yz = buffer_1010_ddsd[52];

    auto g_x_0_x_0_xy_xz_0_zz = buffer_1010_ddsd[53];

    auto g_x_0_x_0_xy_yy_0_xx = buffer_1010_ddsd[54];

    auto g_x_0_x_0_xy_yy_0_xy = buffer_1010_ddsd[55];

    auto g_x_0_x_0_xy_yy_0_xz = buffer_1010_ddsd[56];

    auto g_x_0_x_0_xy_yy_0_yy = buffer_1010_ddsd[57];

    auto g_x_0_x_0_xy_yy_0_yz = buffer_1010_ddsd[58];

    auto g_x_0_x_0_xy_yy_0_zz = buffer_1010_ddsd[59];

    auto g_x_0_x_0_xy_yz_0_xx = buffer_1010_ddsd[60];

    auto g_x_0_x_0_xy_yz_0_xy = buffer_1010_ddsd[61];

    auto g_x_0_x_0_xy_yz_0_xz = buffer_1010_ddsd[62];

    auto g_x_0_x_0_xy_yz_0_yy = buffer_1010_ddsd[63];

    auto g_x_0_x_0_xy_yz_0_yz = buffer_1010_ddsd[64];

    auto g_x_0_x_0_xy_yz_0_zz = buffer_1010_ddsd[65];

    auto g_x_0_x_0_xy_zz_0_xx = buffer_1010_ddsd[66];

    auto g_x_0_x_0_xy_zz_0_xy = buffer_1010_ddsd[67];

    auto g_x_0_x_0_xy_zz_0_xz = buffer_1010_ddsd[68];

    auto g_x_0_x_0_xy_zz_0_yy = buffer_1010_ddsd[69];

    auto g_x_0_x_0_xy_zz_0_yz = buffer_1010_ddsd[70];

    auto g_x_0_x_0_xy_zz_0_zz = buffer_1010_ddsd[71];

    auto g_x_0_x_0_xz_xx_0_xx = buffer_1010_ddsd[72];

    auto g_x_0_x_0_xz_xx_0_xy = buffer_1010_ddsd[73];

    auto g_x_0_x_0_xz_xx_0_xz = buffer_1010_ddsd[74];

    auto g_x_0_x_0_xz_xx_0_yy = buffer_1010_ddsd[75];

    auto g_x_0_x_0_xz_xx_0_yz = buffer_1010_ddsd[76];

    auto g_x_0_x_0_xz_xx_0_zz = buffer_1010_ddsd[77];

    auto g_x_0_x_0_xz_xy_0_xx = buffer_1010_ddsd[78];

    auto g_x_0_x_0_xz_xy_0_xy = buffer_1010_ddsd[79];

    auto g_x_0_x_0_xz_xy_0_xz = buffer_1010_ddsd[80];

    auto g_x_0_x_0_xz_xy_0_yy = buffer_1010_ddsd[81];

    auto g_x_0_x_0_xz_xy_0_yz = buffer_1010_ddsd[82];

    auto g_x_0_x_0_xz_xy_0_zz = buffer_1010_ddsd[83];

    auto g_x_0_x_0_xz_xz_0_xx = buffer_1010_ddsd[84];

    auto g_x_0_x_0_xz_xz_0_xy = buffer_1010_ddsd[85];

    auto g_x_0_x_0_xz_xz_0_xz = buffer_1010_ddsd[86];

    auto g_x_0_x_0_xz_xz_0_yy = buffer_1010_ddsd[87];

    auto g_x_0_x_0_xz_xz_0_yz = buffer_1010_ddsd[88];

    auto g_x_0_x_0_xz_xz_0_zz = buffer_1010_ddsd[89];

    auto g_x_0_x_0_xz_yy_0_xx = buffer_1010_ddsd[90];

    auto g_x_0_x_0_xz_yy_0_xy = buffer_1010_ddsd[91];

    auto g_x_0_x_0_xz_yy_0_xz = buffer_1010_ddsd[92];

    auto g_x_0_x_0_xz_yy_0_yy = buffer_1010_ddsd[93];

    auto g_x_0_x_0_xz_yy_0_yz = buffer_1010_ddsd[94];

    auto g_x_0_x_0_xz_yy_0_zz = buffer_1010_ddsd[95];

    auto g_x_0_x_0_xz_yz_0_xx = buffer_1010_ddsd[96];

    auto g_x_0_x_0_xz_yz_0_xy = buffer_1010_ddsd[97];

    auto g_x_0_x_0_xz_yz_0_xz = buffer_1010_ddsd[98];

    auto g_x_0_x_0_xz_yz_0_yy = buffer_1010_ddsd[99];

    auto g_x_0_x_0_xz_yz_0_yz = buffer_1010_ddsd[100];

    auto g_x_0_x_0_xz_yz_0_zz = buffer_1010_ddsd[101];

    auto g_x_0_x_0_xz_zz_0_xx = buffer_1010_ddsd[102];

    auto g_x_0_x_0_xz_zz_0_xy = buffer_1010_ddsd[103];

    auto g_x_0_x_0_xz_zz_0_xz = buffer_1010_ddsd[104];

    auto g_x_0_x_0_xz_zz_0_yy = buffer_1010_ddsd[105];

    auto g_x_0_x_0_xz_zz_0_yz = buffer_1010_ddsd[106];

    auto g_x_0_x_0_xz_zz_0_zz = buffer_1010_ddsd[107];

    auto g_x_0_x_0_yy_xx_0_xx = buffer_1010_ddsd[108];

    auto g_x_0_x_0_yy_xx_0_xy = buffer_1010_ddsd[109];

    auto g_x_0_x_0_yy_xx_0_xz = buffer_1010_ddsd[110];

    auto g_x_0_x_0_yy_xx_0_yy = buffer_1010_ddsd[111];

    auto g_x_0_x_0_yy_xx_0_yz = buffer_1010_ddsd[112];

    auto g_x_0_x_0_yy_xx_0_zz = buffer_1010_ddsd[113];

    auto g_x_0_x_0_yy_xy_0_xx = buffer_1010_ddsd[114];

    auto g_x_0_x_0_yy_xy_0_xy = buffer_1010_ddsd[115];

    auto g_x_0_x_0_yy_xy_0_xz = buffer_1010_ddsd[116];

    auto g_x_0_x_0_yy_xy_0_yy = buffer_1010_ddsd[117];

    auto g_x_0_x_0_yy_xy_0_yz = buffer_1010_ddsd[118];

    auto g_x_0_x_0_yy_xy_0_zz = buffer_1010_ddsd[119];

    auto g_x_0_x_0_yy_xz_0_xx = buffer_1010_ddsd[120];

    auto g_x_0_x_0_yy_xz_0_xy = buffer_1010_ddsd[121];

    auto g_x_0_x_0_yy_xz_0_xz = buffer_1010_ddsd[122];

    auto g_x_0_x_0_yy_xz_0_yy = buffer_1010_ddsd[123];

    auto g_x_0_x_0_yy_xz_0_yz = buffer_1010_ddsd[124];

    auto g_x_0_x_0_yy_xz_0_zz = buffer_1010_ddsd[125];

    auto g_x_0_x_0_yy_yy_0_xx = buffer_1010_ddsd[126];

    auto g_x_0_x_0_yy_yy_0_xy = buffer_1010_ddsd[127];

    auto g_x_0_x_0_yy_yy_0_xz = buffer_1010_ddsd[128];

    auto g_x_0_x_0_yy_yy_0_yy = buffer_1010_ddsd[129];

    auto g_x_0_x_0_yy_yy_0_yz = buffer_1010_ddsd[130];

    auto g_x_0_x_0_yy_yy_0_zz = buffer_1010_ddsd[131];

    auto g_x_0_x_0_yy_yz_0_xx = buffer_1010_ddsd[132];

    auto g_x_0_x_0_yy_yz_0_xy = buffer_1010_ddsd[133];

    auto g_x_0_x_0_yy_yz_0_xz = buffer_1010_ddsd[134];

    auto g_x_0_x_0_yy_yz_0_yy = buffer_1010_ddsd[135];

    auto g_x_0_x_0_yy_yz_0_yz = buffer_1010_ddsd[136];

    auto g_x_0_x_0_yy_yz_0_zz = buffer_1010_ddsd[137];

    auto g_x_0_x_0_yy_zz_0_xx = buffer_1010_ddsd[138];

    auto g_x_0_x_0_yy_zz_0_xy = buffer_1010_ddsd[139];

    auto g_x_0_x_0_yy_zz_0_xz = buffer_1010_ddsd[140];

    auto g_x_0_x_0_yy_zz_0_yy = buffer_1010_ddsd[141];

    auto g_x_0_x_0_yy_zz_0_yz = buffer_1010_ddsd[142];

    auto g_x_0_x_0_yy_zz_0_zz = buffer_1010_ddsd[143];

    auto g_x_0_x_0_yz_xx_0_xx = buffer_1010_ddsd[144];

    auto g_x_0_x_0_yz_xx_0_xy = buffer_1010_ddsd[145];

    auto g_x_0_x_0_yz_xx_0_xz = buffer_1010_ddsd[146];

    auto g_x_0_x_0_yz_xx_0_yy = buffer_1010_ddsd[147];

    auto g_x_0_x_0_yz_xx_0_yz = buffer_1010_ddsd[148];

    auto g_x_0_x_0_yz_xx_0_zz = buffer_1010_ddsd[149];

    auto g_x_0_x_0_yz_xy_0_xx = buffer_1010_ddsd[150];

    auto g_x_0_x_0_yz_xy_0_xy = buffer_1010_ddsd[151];

    auto g_x_0_x_0_yz_xy_0_xz = buffer_1010_ddsd[152];

    auto g_x_0_x_0_yz_xy_0_yy = buffer_1010_ddsd[153];

    auto g_x_0_x_0_yz_xy_0_yz = buffer_1010_ddsd[154];

    auto g_x_0_x_0_yz_xy_0_zz = buffer_1010_ddsd[155];

    auto g_x_0_x_0_yz_xz_0_xx = buffer_1010_ddsd[156];

    auto g_x_0_x_0_yz_xz_0_xy = buffer_1010_ddsd[157];

    auto g_x_0_x_0_yz_xz_0_xz = buffer_1010_ddsd[158];

    auto g_x_0_x_0_yz_xz_0_yy = buffer_1010_ddsd[159];

    auto g_x_0_x_0_yz_xz_0_yz = buffer_1010_ddsd[160];

    auto g_x_0_x_0_yz_xz_0_zz = buffer_1010_ddsd[161];

    auto g_x_0_x_0_yz_yy_0_xx = buffer_1010_ddsd[162];

    auto g_x_0_x_0_yz_yy_0_xy = buffer_1010_ddsd[163];

    auto g_x_0_x_0_yz_yy_0_xz = buffer_1010_ddsd[164];

    auto g_x_0_x_0_yz_yy_0_yy = buffer_1010_ddsd[165];

    auto g_x_0_x_0_yz_yy_0_yz = buffer_1010_ddsd[166];

    auto g_x_0_x_0_yz_yy_0_zz = buffer_1010_ddsd[167];

    auto g_x_0_x_0_yz_yz_0_xx = buffer_1010_ddsd[168];

    auto g_x_0_x_0_yz_yz_0_xy = buffer_1010_ddsd[169];

    auto g_x_0_x_0_yz_yz_0_xz = buffer_1010_ddsd[170];

    auto g_x_0_x_0_yz_yz_0_yy = buffer_1010_ddsd[171];

    auto g_x_0_x_0_yz_yz_0_yz = buffer_1010_ddsd[172];

    auto g_x_0_x_0_yz_yz_0_zz = buffer_1010_ddsd[173];

    auto g_x_0_x_0_yz_zz_0_xx = buffer_1010_ddsd[174];

    auto g_x_0_x_0_yz_zz_0_xy = buffer_1010_ddsd[175];

    auto g_x_0_x_0_yz_zz_0_xz = buffer_1010_ddsd[176];

    auto g_x_0_x_0_yz_zz_0_yy = buffer_1010_ddsd[177];

    auto g_x_0_x_0_yz_zz_0_yz = buffer_1010_ddsd[178];

    auto g_x_0_x_0_yz_zz_0_zz = buffer_1010_ddsd[179];

    auto g_x_0_x_0_zz_xx_0_xx = buffer_1010_ddsd[180];

    auto g_x_0_x_0_zz_xx_0_xy = buffer_1010_ddsd[181];

    auto g_x_0_x_0_zz_xx_0_xz = buffer_1010_ddsd[182];

    auto g_x_0_x_0_zz_xx_0_yy = buffer_1010_ddsd[183];

    auto g_x_0_x_0_zz_xx_0_yz = buffer_1010_ddsd[184];

    auto g_x_0_x_0_zz_xx_0_zz = buffer_1010_ddsd[185];

    auto g_x_0_x_0_zz_xy_0_xx = buffer_1010_ddsd[186];

    auto g_x_0_x_0_zz_xy_0_xy = buffer_1010_ddsd[187];

    auto g_x_0_x_0_zz_xy_0_xz = buffer_1010_ddsd[188];

    auto g_x_0_x_0_zz_xy_0_yy = buffer_1010_ddsd[189];

    auto g_x_0_x_0_zz_xy_0_yz = buffer_1010_ddsd[190];

    auto g_x_0_x_0_zz_xy_0_zz = buffer_1010_ddsd[191];

    auto g_x_0_x_0_zz_xz_0_xx = buffer_1010_ddsd[192];

    auto g_x_0_x_0_zz_xz_0_xy = buffer_1010_ddsd[193];

    auto g_x_0_x_0_zz_xz_0_xz = buffer_1010_ddsd[194];

    auto g_x_0_x_0_zz_xz_0_yy = buffer_1010_ddsd[195];

    auto g_x_0_x_0_zz_xz_0_yz = buffer_1010_ddsd[196];

    auto g_x_0_x_0_zz_xz_0_zz = buffer_1010_ddsd[197];

    auto g_x_0_x_0_zz_yy_0_xx = buffer_1010_ddsd[198];

    auto g_x_0_x_0_zz_yy_0_xy = buffer_1010_ddsd[199];

    auto g_x_0_x_0_zz_yy_0_xz = buffer_1010_ddsd[200];

    auto g_x_0_x_0_zz_yy_0_yy = buffer_1010_ddsd[201];

    auto g_x_0_x_0_zz_yy_0_yz = buffer_1010_ddsd[202];

    auto g_x_0_x_0_zz_yy_0_zz = buffer_1010_ddsd[203];

    auto g_x_0_x_0_zz_yz_0_xx = buffer_1010_ddsd[204];

    auto g_x_0_x_0_zz_yz_0_xy = buffer_1010_ddsd[205];

    auto g_x_0_x_0_zz_yz_0_xz = buffer_1010_ddsd[206];

    auto g_x_0_x_0_zz_yz_0_yy = buffer_1010_ddsd[207];

    auto g_x_0_x_0_zz_yz_0_yz = buffer_1010_ddsd[208];

    auto g_x_0_x_0_zz_yz_0_zz = buffer_1010_ddsd[209];

    auto g_x_0_x_0_zz_zz_0_xx = buffer_1010_ddsd[210];

    auto g_x_0_x_0_zz_zz_0_xy = buffer_1010_ddsd[211];

    auto g_x_0_x_0_zz_zz_0_xz = buffer_1010_ddsd[212];

    auto g_x_0_x_0_zz_zz_0_yy = buffer_1010_ddsd[213];

    auto g_x_0_x_0_zz_zz_0_yz = buffer_1010_ddsd[214];

    auto g_x_0_x_0_zz_zz_0_zz = buffer_1010_ddsd[215];

    auto g_x_0_y_0_xx_xx_0_xx = buffer_1010_ddsd[216];

    auto g_x_0_y_0_xx_xx_0_xy = buffer_1010_ddsd[217];

    auto g_x_0_y_0_xx_xx_0_xz = buffer_1010_ddsd[218];

    auto g_x_0_y_0_xx_xx_0_yy = buffer_1010_ddsd[219];

    auto g_x_0_y_0_xx_xx_0_yz = buffer_1010_ddsd[220];

    auto g_x_0_y_0_xx_xx_0_zz = buffer_1010_ddsd[221];

    auto g_x_0_y_0_xx_xy_0_xx = buffer_1010_ddsd[222];

    auto g_x_0_y_0_xx_xy_0_xy = buffer_1010_ddsd[223];

    auto g_x_0_y_0_xx_xy_0_xz = buffer_1010_ddsd[224];

    auto g_x_0_y_0_xx_xy_0_yy = buffer_1010_ddsd[225];

    auto g_x_0_y_0_xx_xy_0_yz = buffer_1010_ddsd[226];

    auto g_x_0_y_0_xx_xy_0_zz = buffer_1010_ddsd[227];

    auto g_x_0_y_0_xx_xz_0_xx = buffer_1010_ddsd[228];

    auto g_x_0_y_0_xx_xz_0_xy = buffer_1010_ddsd[229];

    auto g_x_0_y_0_xx_xz_0_xz = buffer_1010_ddsd[230];

    auto g_x_0_y_0_xx_xz_0_yy = buffer_1010_ddsd[231];

    auto g_x_0_y_0_xx_xz_0_yz = buffer_1010_ddsd[232];

    auto g_x_0_y_0_xx_xz_0_zz = buffer_1010_ddsd[233];

    auto g_x_0_y_0_xx_yy_0_xx = buffer_1010_ddsd[234];

    auto g_x_0_y_0_xx_yy_0_xy = buffer_1010_ddsd[235];

    auto g_x_0_y_0_xx_yy_0_xz = buffer_1010_ddsd[236];

    auto g_x_0_y_0_xx_yy_0_yy = buffer_1010_ddsd[237];

    auto g_x_0_y_0_xx_yy_0_yz = buffer_1010_ddsd[238];

    auto g_x_0_y_0_xx_yy_0_zz = buffer_1010_ddsd[239];

    auto g_x_0_y_0_xx_yz_0_xx = buffer_1010_ddsd[240];

    auto g_x_0_y_0_xx_yz_0_xy = buffer_1010_ddsd[241];

    auto g_x_0_y_0_xx_yz_0_xz = buffer_1010_ddsd[242];

    auto g_x_0_y_0_xx_yz_0_yy = buffer_1010_ddsd[243];

    auto g_x_0_y_0_xx_yz_0_yz = buffer_1010_ddsd[244];

    auto g_x_0_y_0_xx_yz_0_zz = buffer_1010_ddsd[245];

    auto g_x_0_y_0_xx_zz_0_xx = buffer_1010_ddsd[246];

    auto g_x_0_y_0_xx_zz_0_xy = buffer_1010_ddsd[247];

    auto g_x_0_y_0_xx_zz_0_xz = buffer_1010_ddsd[248];

    auto g_x_0_y_0_xx_zz_0_yy = buffer_1010_ddsd[249];

    auto g_x_0_y_0_xx_zz_0_yz = buffer_1010_ddsd[250];

    auto g_x_0_y_0_xx_zz_0_zz = buffer_1010_ddsd[251];

    auto g_x_0_y_0_xy_xx_0_xx = buffer_1010_ddsd[252];

    auto g_x_0_y_0_xy_xx_0_xy = buffer_1010_ddsd[253];

    auto g_x_0_y_0_xy_xx_0_xz = buffer_1010_ddsd[254];

    auto g_x_0_y_0_xy_xx_0_yy = buffer_1010_ddsd[255];

    auto g_x_0_y_0_xy_xx_0_yz = buffer_1010_ddsd[256];

    auto g_x_0_y_0_xy_xx_0_zz = buffer_1010_ddsd[257];

    auto g_x_0_y_0_xy_xy_0_xx = buffer_1010_ddsd[258];

    auto g_x_0_y_0_xy_xy_0_xy = buffer_1010_ddsd[259];

    auto g_x_0_y_0_xy_xy_0_xz = buffer_1010_ddsd[260];

    auto g_x_0_y_0_xy_xy_0_yy = buffer_1010_ddsd[261];

    auto g_x_0_y_0_xy_xy_0_yz = buffer_1010_ddsd[262];

    auto g_x_0_y_0_xy_xy_0_zz = buffer_1010_ddsd[263];

    auto g_x_0_y_0_xy_xz_0_xx = buffer_1010_ddsd[264];

    auto g_x_0_y_0_xy_xz_0_xy = buffer_1010_ddsd[265];

    auto g_x_0_y_0_xy_xz_0_xz = buffer_1010_ddsd[266];

    auto g_x_0_y_0_xy_xz_0_yy = buffer_1010_ddsd[267];

    auto g_x_0_y_0_xy_xz_0_yz = buffer_1010_ddsd[268];

    auto g_x_0_y_0_xy_xz_0_zz = buffer_1010_ddsd[269];

    auto g_x_0_y_0_xy_yy_0_xx = buffer_1010_ddsd[270];

    auto g_x_0_y_0_xy_yy_0_xy = buffer_1010_ddsd[271];

    auto g_x_0_y_0_xy_yy_0_xz = buffer_1010_ddsd[272];

    auto g_x_0_y_0_xy_yy_0_yy = buffer_1010_ddsd[273];

    auto g_x_0_y_0_xy_yy_0_yz = buffer_1010_ddsd[274];

    auto g_x_0_y_0_xy_yy_0_zz = buffer_1010_ddsd[275];

    auto g_x_0_y_0_xy_yz_0_xx = buffer_1010_ddsd[276];

    auto g_x_0_y_0_xy_yz_0_xy = buffer_1010_ddsd[277];

    auto g_x_0_y_0_xy_yz_0_xz = buffer_1010_ddsd[278];

    auto g_x_0_y_0_xy_yz_0_yy = buffer_1010_ddsd[279];

    auto g_x_0_y_0_xy_yz_0_yz = buffer_1010_ddsd[280];

    auto g_x_0_y_0_xy_yz_0_zz = buffer_1010_ddsd[281];

    auto g_x_0_y_0_xy_zz_0_xx = buffer_1010_ddsd[282];

    auto g_x_0_y_0_xy_zz_0_xy = buffer_1010_ddsd[283];

    auto g_x_0_y_0_xy_zz_0_xz = buffer_1010_ddsd[284];

    auto g_x_0_y_0_xy_zz_0_yy = buffer_1010_ddsd[285];

    auto g_x_0_y_0_xy_zz_0_yz = buffer_1010_ddsd[286];

    auto g_x_0_y_0_xy_zz_0_zz = buffer_1010_ddsd[287];

    auto g_x_0_y_0_xz_xx_0_xx = buffer_1010_ddsd[288];

    auto g_x_0_y_0_xz_xx_0_xy = buffer_1010_ddsd[289];

    auto g_x_0_y_0_xz_xx_0_xz = buffer_1010_ddsd[290];

    auto g_x_0_y_0_xz_xx_0_yy = buffer_1010_ddsd[291];

    auto g_x_0_y_0_xz_xx_0_yz = buffer_1010_ddsd[292];

    auto g_x_0_y_0_xz_xx_0_zz = buffer_1010_ddsd[293];

    auto g_x_0_y_0_xz_xy_0_xx = buffer_1010_ddsd[294];

    auto g_x_0_y_0_xz_xy_0_xy = buffer_1010_ddsd[295];

    auto g_x_0_y_0_xz_xy_0_xz = buffer_1010_ddsd[296];

    auto g_x_0_y_0_xz_xy_0_yy = buffer_1010_ddsd[297];

    auto g_x_0_y_0_xz_xy_0_yz = buffer_1010_ddsd[298];

    auto g_x_0_y_0_xz_xy_0_zz = buffer_1010_ddsd[299];

    auto g_x_0_y_0_xz_xz_0_xx = buffer_1010_ddsd[300];

    auto g_x_0_y_0_xz_xz_0_xy = buffer_1010_ddsd[301];

    auto g_x_0_y_0_xz_xz_0_xz = buffer_1010_ddsd[302];

    auto g_x_0_y_0_xz_xz_0_yy = buffer_1010_ddsd[303];

    auto g_x_0_y_0_xz_xz_0_yz = buffer_1010_ddsd[304];

    auto g_x_0_y_0_xz_xz_0_zz = buffer_1010_ddsd[305];

    auto g_x_0_y_0_xz_yy_0_xx = buffer_1010_ddsd[306];

    auto g_x_0_y_0_xz_yy_0_xy = buffer_1010_ddsd[307];

    auto g_x_0_y_0_xz_yy_0_xz = buffer_1010_ddsd[308];

    auto g_x_0_y_0_xz_yy_0_yy = buffer_1010_ddsd[309];

    auto g_x_0_y_0_xz_yy_0_yz = buffer_1010_ddsd[310];

    auto g_x_0_y_0_xz_yy_0_zz = buffer_1010_ddsd[311];

    auto g_x_0_y_0_xz_yz_0_xx = buffer_1010_ddsd[312];

    auto g_x_0_y_0_xz_yz_0_xy = buffer_1010_ddsd[313];

    auto g_x_0_y_0_xz_yz_0_xz = buffer_1010_ddsd[314];

    auto g_x_0_y_0_xz_yz_0_yy = buffer_1010_ddsd[315];

    auto g_x_0_y_0_xz_yz_0_yz = buffer_1010_ddsd[316];

    auto g_x_0_y_0_xz_yz_0_zz = buffer_1010_ddsd[317];

    auto g_x_0_y_0_xz_zz_0_xx = buffer_1010_ddsd[318];

    auto g_x_0_y_0_xz_zz_0_xy = buffer_1010_ddsd[319];

    auto g_x_0_y_0_xz_zz_0_xz = buffer_1010_ddsd[320];

    auto g_x_0_y_0_xz_zz_0_yy = buffer_1010_ddsd[321];

    auto g_x_0_y_0_xz_zz_0_yz = buffer_1010_ddsd[322];

    auto g_x_0_y_0_xz_zz_0_zz = buffer_1010_ddsd[323];

    auto g_x_0_y_0_yy_xx_0_xx = buffer_1010_ddsd[324];

    auto g_x_0_y_0_yy_xx_0_xy = buffer_1010_ddsd[325];

    auto g_x_0_y_0_yy_xx_0_xz = buffer_1010_ddsd[326];

    auto g_x_0_y_0_yy_xx_0_yy = buffer_1010_ddsd[327];

    auto g_x_0_y_0_yy_xx_0_yz = buffer_1010_ddsd[328];

    auto g_x_0_y_0_yy_xx_0_zz = buffer_1010_ddsd[329];

    auto g_x_0_y_0_yy_xy_0_xx = buffer_1010_ddsd[330];

    auto g_x_0_y_0_yy_xy_0_xy = buffer_1010_ddsd[331];

    auto g_x_0_y_0_yy_xy_0_xz = buffer_1010_ddsd[332];

    auto g_x_0_y_0_yy_xy_0_yy = buffer_1010_ddsd[333];

    auto g_x_0_y_0_yy_xy_0_yz = buffer_1010_ddsd[334];

    auto g_x_0_y_0_yy_xy_0_zz = buffer_1010_ddsd[335];

    auto g_x_0_y_0_yy_xz_0_xx = buffer_1010_ddsd[336];

    auto g_x_0_y_0_yy_xz_0_xy = buffer_1010_ddsd[337];

    auto g_x_0_y_0_yy_xz_0_xz = buffer_1010_ddsd[338];

    auto g_x_0_y_0_yy_xz_0_yy = buffer_1010_ddsd[339];

    auto g_x_0_y_0_yy_xz_0_yz = buffer_1010_ddsd[340];

    auto g_x_0_y_0_yy_xz_0_zz = buffer_1010_ddsd[341];

    auto g_x_0_y_0_yy_yy_0_xx = buffer_1010_ddsd[342];

    auto g_x_0_y_0_yy_yy_0_xy = buffer_1010_ddsd[343];

    auto g_x_0_y_0_yy_yy_0_xz = buffer_1010_ddsd[344];

    auto g_x_0_y_0_yy_yy_0_yy = buffer_1010_ddsd[345];

    auto g_x_0_y_0_yy_yy_0_yz = buffer_1010_ddsd[346];

    auto g_x_0_y_0_yy_yy_0_zz = buffer_1010_ddsd[347];

    auto g_x_0_y_0_yy_yz_0_xx = buffer_1010_ddsd[348];

    auto g_x_0_y_0_yy_yz_0_xy = buffer_1010_ddsd[349];

    auto g_x_0_y_0_yy_yz_0_xz = buffer_1010_ddsd[350];

    auto g_x_0_y_0_yy_yz_0_yy = buffer_1010_ddsd[351];

    auto g_x_0_y_0_yy_yz_0_yz = buffer_1010_ddsd[352];

    auto g_x_0_y_0_yy_yz_0_zz = buffer_1010_ddsd[353];

    auto g_x_0_y_0_yy_zz_0_xx = buffer_1010_ddsd[354];

    auto g_x_0_y_0_yy_zz_0_xy = buffer_1010_ddsd[355];

    auto g_x_0_y_0_yy_zz_0_xz = buffer_1010_ddsd[356];

    auto g_x_0_y_0_yy_zz_0_yy = buffer_1010_ddsd[357];

    auto g_x_0_y_0_yy_zz_0_yz = buffer_1010_ddsd[358];

    auto g_x_0_y_0_yy_zz_0_zz = buffer_1010_ddsd[359];

    auto g_x_0_y_0_yz_xx_0_xx = buffer_1010_ddsd[360];

    auto g_x_0_y_0_yz_xx_0_xy = buffer_1010_ddsd[361];

    auto g_x_0_y_0_yz_xx_0_xz = buffer_1010_ddsd[362];

    auto g_x_0_y_0_yz_xx_0_yy = buffer_1010_ddsd[363];

    auto g_x_0_y_0_yz_xx_0_yz = buffer_1010_ddsd[364];

    auto g_x_0_y_0_yz_xx_0_zz = buffer_1010_ddsd[365];

    auto g_x_0_y_0_yz_xy_0_xx = buffer_1010_ddsd[366];

    auto g_x_0_y_0_yz_xy_0_xy = buffer_1010_ddsd[367];

    auto g_x_0_y_0_yz_xy_0_xz = buffer_1010_ddsd[368];

    auto g_x_0_y_0_yz_xy_0_yy = buffer_1010_ddsd[369];

    auto g_x_0_y_0_yz_xy_0_yz = buffer_1010_ddsd[370];

    auto g_x_0_y_0_yz_xy_0_zz = buffer_1010_ddsd[371];

    auto g_x_0_y_0_yz_xz_0_xx = buffer_1010_ddsd[372];

    auto g_x_0_y_0_yz_xz_0_xy = buffer_1010_ddsd[373];

    auto g_x_0_y_0_yz_xz_0_xz = buffer_1010_ddsd[374];

    auto g_x_0_y_0_yz_xz_0_yy = buffer_1010_ddsd[375];

    auto g_x_0_y_0_yz_xz_0_yz = buffer_1010_ddsd[376];

    auto g_x_0_y_0_yz_xz_0_zz = buffer_1010_ddsd[377];

    auto g_x_0_y_0_yz_yy_0_xx = buffer_1010_ddsd[378];

    auto g_x_0_y_0_yz_yy_0_xy = buffer_1010_ddsd[379];

    auto g_x_0_y_0_yz_yy_0_xz = buffer_1010_ddsd[380];

    auto g_x_0_y_0_yz_yy_0_yy = buffer_1010_ddsd[381];

    auto g_x_0_y_0_yz_yy_0_yz = buffer_1010_ddsd[382];

    auto g_x_0_y_0_yz_yy_0_zz = buffer_1010_ddsd[383];

    auto g_x_0_y_0_yz_yz_0_xx = buffer_1010_ddsd[384];

    auto g_x_0_y_0_yz_yz_0_xy = buffer_1010_ddsd[385];

    auto g_x_0_y_0_yz_yz_0_xz = buffer_1010_ddsd[386];

    auto g_x_0_y_0_yz_yz_0_yy = buffer_1010_ddsd[387];

    auto g_x_0_y_0_yz_yz_0_yz = buffer_1010_ddsd[388];

    auto g_x_0_y_0_yz_yz_0_zz = buffer_1010_ddsd[389];

    auto g_x_0_y_0_yz_zz_0_xx = buffer_1010_ddsd[390];

    auto g_x_0_y_0_yz_zz_0_xy = buffer_1010_ddsd[391];

    auto g_x_0_y_0_yz_zz_0_xz = buffer_1010_ddsd[392];

    auto g_x_0_y_0_yz_zz_0_yy = buffer_1010_ddsd[393];

    auto g_x_0_y_0_yz_zz_0_yz = buffer_1010_ddsd[394];

    auto g_x_0_y_0_yz_zz_0_zz = buffer_1010_ddsd[395];

    auto g_x_0_y_0_zz_xx_0_xx = buffer_1010_ddsd[396];

    auto g_x_0_y_0_zz_xx_0_xy = buffer_1010_ddsd[397];

    auto g_x_0_y_0_zz_xx_0_xz = buffer_1010_ddsd[398];

    auto g_x_0_y_0_zz_xx_0_yy = buffer_1010_ddsd[399];

    auto g_x_0_y_0_zz_xx_0_yz = buffer_1010_ddsd[400];

    auto g_x_0_y_0_zz_xx_0_zz = buffer_1010_ddsd[401];

    auto g_x_0_y_0_zz_xy_0_xx = buffer_1010_ddsd[402];

    auto g_x_0_y_0_zz_xy_0_xy = buffer_1010_ddsd[403];

    auto g_x_0_y_0_zz_xy_0_xz = buffer_1010_ddsd[404];

    auto g_x_0_y_0_zz_xy_0_yy = buffer_1010_ddsd[405];

    auto g_x_0_y_0_zz_xy_0_yz = buffer_1010_ddsd[406];

    auto g_x_0_y_0_zz_xy_0_zz = buffer_1010_ddsd[407];

    auto g_x_0_y_0_zz_xz_0_xx = buffer_1010_ddsd[408];

    auto g_x_0_y_0_zz_xz_0_xy = buffer_1010_ddsd[409];

    auto g_x_0_y_0_zz_xz_0_xz = buffer_1010_ddsd[410];

    auto g_x_0_y_0_zz_xz_0_yy = buffer_1010_ddsd[411];

    auto g_x_0_y_0_zz_xz_0_yz = buffer_1010_ddsd[412];

    auto g_x_0_y_0_zz_xz_0_zz = buffer_1010_ddsd[413];

    auto g_x_0_y_0_zz_yy_0_xx = buffer_1010_ddsd[414];

    auto g_x_0_y_0_zz_yy_0_xy = buffer_1010_ddsd[415];

    auto g_x_0_y_0_zz_yy_0_xz = buffer_1010_ddsd[416];

    auto g_x_0_y_0_zz_yy_0_yy = buffer_1010_ddsd[417];

    auto g_x_0_y_0_zz_yy_0_yz = buffer_1010_ddsd[418];

    auto g_x_0_y_0_zz_yy_0_zz = buffer_1010_ddsd[419];

    auto g_x_0_y_0_zz_yz_0_xx = buffer_1010_ddsd[420];

    auto g_x_0_y_0_zz_yz_0_xy = buffer_1010_ddsd[421];

    auto g_x_0_y_0_zz_yz_0_xz = buffer_1010_ddsd[422];

    auto g_x_0_y_0_zz_yz_0_yy = buffer_1010_ddsd[423];

    auto g_x_0_y_0_zz_yz_0_yz = buffer_1010_ddsd[424];

    auto g_x_0_y_0_zz_yz_0_zz = buffer_1010_ddsd[425];

    auto g_x_0_y_0_zz_zz_0_xx = buffer_1010_ddsd[426];

    auto g_x_0_y_0_zz_zz_0_xy = buffer_1010_ddsd[427];

    auto g_x_0_y_0_zz_zz_0_xz = buffer_1010_ddsd[428];

    auto g_x_0_y_0_zz_zz_0_yy = buffer_1010_ddsd[429];

    auto g_x_0_y_0_zz_zz_0_yz = buffer_1010_ddsd[430];

    auto g_x_0_y_0_zz_zz_0_zz = buffer_1010_ddsd[431];

    auto g_x_0_z_0_xx_xx_0_xx = buffer_1010_ddsd[432];

    auto g_x_0_z_0_xx_xx_0_xy = buffer_1010_ddsd[433];

    auto g_x_0_z_0_xx_xx_0_xz = buffer_1010_ddsd[434];

    auto g_x_0_z_0_xx_xx_0_yy = buffer_1010_ddsd[435];

    auto g_x_0_z_0_xx_xx_0_yz = buffer_1010_ddsd[436];

    auto g_x_0_z_0_xx_xx_0_zz = buffer_1010_ddsd[437];

    auto g_x_0_z_0_xx_xy_0_xx = buffer_1010_ddsd[438];

    auto g_x_0_z_0_xx_xy_0_xy = buffer_1010_ddsd[439];

    auto g_x_0_z_0_xx_xy_0_xz = buffer_1010_ddsd[440];

    auto g_x_0_z_0_xx_xy_0_yy = buffer_1010_ddsd[441];

    auto g_x_0_z_0_xx_xy_0_yz = buffer_1010_ddsd[442];

    auto g_x_0_z_0_xx_xy_0_zz = buffer_1010_ddsd[443];

    auto g_x_0_z_0_xx_xz_0_xx = buffer_1010_ddsd[444];

    auto g_x_0_z_0_xx_xz_0_xy = buffer_1010_ddsd[445];

    auto g_x_0_z_0_xx_xz_0_xz = buffer_1010_ddsd[446];

    auto g_x_0_z_0_xx_xz_0_yy = buffer_1010_ddsd[447];

    auto g_x_0_z_0_xx_xz_0_yz = buffer_1010_ddsd[448];

    auto g_x_0_z_0_xx_xz_0_zz = buffer_1010_ddsd[449];

    auto g_x_0_z_0_xx_yy_0_xx = buffer_1010_ddsd[450];

    auto g_x_0_z_0_xx_yy_0_xy = buffer_1010_ddsd[451];

    auto g_x_0_z_0_xx_yy_0_xz = buffer_1010_ddsd[452];

    auto g_x_0_z_0_xx_yy_0_yy = buffer_1010_ddsd[453];

    auto g_x_0_z_0_xx_yy_0_yz = buffer_1010_ddsd[454];

    auto g_x_0_z_0_xx_yy_0_zz = buffer_1010_ddsd[455];

    auto g_x_0_z_0_xx_yz_0_xx = buffer_1010_ddsd[456];

    auto g_x_0_z_0_xx_yz_0_xy = buffer_1010_ddsd[457];

    auto g_x_0_z_0_xx_yz_0_xz = buffer_1010_ddsd[458];

    auto g_x_0_z_0_xx_yz_0_yy = buffer_1010_ddsd[459];

    auto g_x_0_z_0_xx_yz_0_yz = buffer_1010_ddsd[460];

    auto g_x_0_z_0_xx_yz_0_zz = buffer_1010_ddsd[461];

    auto g_x_0_z_0_xx_zz_0_xx = buffer_1010_ddsd[462];

    auto g_x_0_z_0_xx_zz_0_xy = buffer_1010_ddsd[463];

    auto g_x_0_z_0_xx_zz_0_xz = buffer_1010_ddsd[464];

    auto g_x_0_z_0_xx_zz_0_yy = buffer_1010_ddsd[465];

    auto g_x_0_z_0_xx_zz_0_yz = buffer_1010_ddsd[466];

    auto g_x_0_z_0_xx_zz_0_zz = buffer_1010_ddsd[467];

    auto g_x_0_z_0_xy_xx_0_xx = buffer_1010_ddsd[468];

    auto g_x_0_z_0_xy_xx_0_xy = buffer_1010_ddsd[469];

    auto g_x_0_z_0_xy_xx_0_xz = buffer_1010_ddsd[470];

    auto g_x_0_z_0_xy_xx_0_yy = buffer_1010_ddsd[471];

    auto g_x_0_z_0_xy_xx_0_yz = buffer_1010_ddsd[472];

    auto g_x_0_z_0_xy_xx_0_zz = buffer_1010_ddsd[473];

    auto g_x_0_z_0_xy_xy_0_xx = buffer_1010_ddsd[474];

    auto g_x_0_z_0_xy_xy_0_xy = buffer_1010_ddsd[475];

    auto g_x_0_z_0_xy_xy_0_xz = buffer_1010_ddsd[476];

    auto g_x_0_z_0_xy_xy_0_yy = buffer_1010_ddsd[477];

    auto g_x_0_z_0_xy_xy_0_yz = buffer_1010_ddsd[478];

    auto g_x_0_z_0_xy_xy_0_zz = buffer_1010_ddsd[479];

    auto g_x_0_z_0_xy_xz_0_xx = buffer_1010_ddsd[480];

    auto g_x_0_z_0_xy_xz_0_xy = buffer_1010_ddsd[481];

    auto g_x_0_z_0_xy_xz_0_xz = buffer_1010_ddsd[482];

    auto g_x_0_z_0_xy_xz_0_yy = buffer_1010_ddsd[483];

    auto g_x_0_z_0_xy_xz_0_yz = buffer_1010_ddsd[484];

    auto g_x_0_z_0_xy_xz_0_zz = buffer_1010_ddsd[485];

    auto g_x_0_z_0_xy_yy_0_xx = buffer_1010_ddsd[486];

    auto g_x_0_z_0_xy_yy_0_xy = buffer_1010_ddsd[487];

    auto g_x_0_z_0_xy_yy_0_xz = buffer_1010_ddsd[488];

    auto g_x_0_z_0_xy_yy_0_yy = buffer_1010_ddsd[489];

    auto g_x_0_z_0_xy_yy_0_yz = buffer_1010_ddsd[490];

    auto g_x_0_z_0_xy_yy_0_zz = buffer_1010_ddsd[491];

    auto g_x_0_z_0_xy_yz_0_xx = buffer_1010_ddsd[492];

    auto g_x_0_z_0_xy_yz_0_xy = buffer_1010_ddsd[493];

    auto g_x_0_z_0_xy_yz_0_xz = buffer_1010_ddsd[494];

    auto g_x_0_z_0_xy_yz_0_yy = buffer_1010_ddsd[495];

    auto g_x_0_z_0_xy_yz_0_yz = buffer_1010_ddsd[496];

    auto g_x_0_z_0_xy_yz_0_zz = buffer_1010_ddsd[497];

    auto g_x_0_z_0_xy_zz_0_xx = buffer_1010_ddsd[498];

    auto g_x_0_z_0_xy_zz_0_xy = buffer_1010_ddsd[499];

    auto g_x_0_z_0_xy_zz_0_xz = buffer_1010_ddsd[500];

    auto g_x_0_z_0_xy_zz_0_yy = buffer_1010_ddsd[501];

    auto g_x_0_z_0_xy_zz_0_yz = buffer_1010_ddsd[502];

    auto g_x_0_z_0_xy_zz_0_zz = buffer_1010_ddsd[503];

    auto g_x_0_z_0_xz_xx_0_xx = buffer_1010_ddsd[504];

    auto g_x_0_z_0_xz_xx_0_xy = buffer_1010_ddsd[505];

    auto g_x_0_z_0_xz_xx_0_xz = buffer_1010_ddsd[506];

    auto g_x_0_z_0_xz_xx_0_yy = buffer_1010_ddsd[507];

    auto g_x_0_z_0_xz_xx_0_yz = buffer_1010_ddsd[508];

    auto g_x_0_z_0_xz_xx_0_zz = buffer_1010_ddsd[509];

    auto g_x_0_z_0_xz_xy_0_xx = buffer_1010_ddsd[510];

    auto g_x_0_z_0_xz_xy_0_xy = buffer_1010_ddsd[511];

    auto g_x_0_z_0_xz_xy_0_xz = buffer_1010_ddsd[512];

    auto g_x_0_z_0_xz_xy_0_yy = buffer_1010_ddsd[513];

    auto g_x_0_z_0_xz_xy_0_yz = buffer_1010_ddsd[514];

    auto g_x_0_z_0_xz_xy_0_zz = buffer_1010_ddsd[515];

    auto g_x_0_z_0_xz_xz_0_xx = buffer_1010_ddsd[516];

    auto g_x_0_z_0_xz_xz_0_xy = buffer_1010_ddsd[517];

    auto g_x_0_z_0_xz_xz_0_xz = buffer_1010_ddsd[518];

    auto g_x_0_z_0_xz_xz_0_yy = buffer_1010_ddsd[519];

    auto g_x_0_z_0_xz_xz_0_yz = buffer_1010_ddsd[520];

    auto g_x_0_z_0_xz_xz_0_zz = buffer_1010_ddsd[521];

    auto g_x_0_z_0_xz_yy_0_xx = buffer_1010_ddsd[522];

    auto g_x_0_z_0_xz_yy_0_xy = buffer_1010_ddsd[523];

    auto g_x_0_z_0_xz_yy_0_xz = buffer_1010_ddsd[524];

    auto g_x_0_z_0_xz_yy_0_yy = buffer_1010_ddsd[525];

    auto g_x_0_z_0_xz_yy_0_yz = buffer_1010_ddsd[526];

    auto g_x_0_z_0_xz_yy_0_zz = buffer_1010_ddsd[527];

    auto g_x_0_z_0_xz_yz_0_xx = buffer_1010_ddsd[528];

    auto g_x_0_z_0_xz_yz_0_xy = buffer_1010_ddsd[529];

    auto g_x_0_z_0_xz_yz_0_xz = buffer_1010_ddsd[530];

    auto g_x_0_z_0_xz_yz_0_yy = buffer_1010_ddsd[531];

    auto g_x_0_z_0_xz_yz_0_yz = buffer_1010_ddsd[532];

    auto g_x_0_z_0_xz_yz_0_zz = buffer_1010_ddsd[533];

    auto g_x_0_z_0_xz_zz_0_xx = buffer_1010_ddsd[534];

    auto g_x_0_z_0_xz_zz_0_xy = buffer_1010_ddsd[535];

    auto g_x_0_z_0_xz_zz_0_xz = buffer_1010_ddsd[536];

    auto g_x_0_z_0_xz_zz_0_yy = buffer_1010_ddsd[537];

    auto g_x_0_z_0_xz_zz_0_yz = buffer_1010_ddsd[538];

    auto g_x_0_z_0_xz_zz_0_zz = buffer_1010_ddsd[539];

    auto g_x_0_z_0_yy_xx_0_xx = buffer_1010_ddsd[540];

    auto g_x_0_z_0_yy_xx_0_xy = buffer_1010_ddsd[541];

    auto g_x_0_z_0_yy_xx_0_xz = buffer_1010_ddsd[542];

    auto g_x_0_z_0_yy_xx_0_yy = buffer_1010_ddsd[543];

    auto g_x_0_z_0_yy_xx_0_yz = buffer_1010_ddsd[544];

    auto g_x_0_z_0_yy_xx_0_zz = buffer_1010_ddsd[545];

    auto g_x_0_z_0_yy_xy_0_xx = buffer_1010_ddsd[546];

    auto g_x_0_z_0_yy_xy_0_xy = buffer_1010_ddsd[547];

    auto g_x_0_z_0_yy_xy_0_xz = buffer_1010_ddsd[548];

    auto g_x_0_z_0_yy_xy_0_yy = buffer_1010_ddsd[549];

    auto g_x_0_z_0_yy_xy_0_yz = buffer_1010_ddsd[550];

    auto g_x_0_z_0_yy_xy_0_zz = buffer_1010_ddsd[551];

    auto g_x_0_z_0_yy_xz_0_xx = buffer_1010_ddsd[552];

    auto g_x_0_z_0_yy_xz_0_xy = buffer_1010_ddsd[553];

    auto g_x_0_z_0_yy_xz_0_xz = buffer_1010_ddsd[554];

    auto g_x_0_z_0_yy_xz_0_yy = buffer_1010_ddsd[555];

    auto g_x_0_z_0_yy_xz_0_yz = buffer_1010_ddsd[556];

    auto g_x_0_z_0_yy_xz_0_zz = buffer_1010_ddsd[557];

    auto g_x_0_z_0_yy_yy_0_xx = buffer_1010_ddsd[558];

    auto g_x_0_z_0_yy_yy_0_xy = buffer_1010_ddsd[559];

    auto g_x_0_z_0_yy_yy_0_xz = buffer_1010_ddsd[560];

    auto g_x_0_z_0_yy_yy_0_yy = buffer_1010_ddsd[561];

    auto g_x_0_z_0_yy_yy_0_yz = buffer_1010_ddsd[562];

    auto g_x_0_z_0_yy_yy_0_zz = buffer_1010_ddsd[563];

    auto g_x_0_z_0_yy_yz_0_xx = buffer_1010_ddsd[564];

    auto g_x_0_z_0_yy_yz_0_xy = buffer_1010_ddsd[565];

    auto g_x_0_z_0_yy_yz_0_xz = buffer_1010_ddsd[566];

    auto g_x_0_z_0_yy_yz_0_yy = buffer_1010_ddsd[567];

    auto g_x_0_z_0_yy_yz_0_yz = buffer_1010_ddsd[568];

    auto g_x_0_z_0_yy_yz_0_zz = buffer_1010_ddsd[569];

    auto g_x_0_z_0_yy_zz_0_xx = buffer_1010_ddsd[570];

    auto g_x_0_z_0_yy_zz_0_xy = buffer_1010_ddsd[571];

    auto g_x_0_z_0_yy_zz_0_xz = buffer_1010_ddsd[572];

    auto g_x_0_z_0_yy_zz_0_yy = buffer_1010_ddsd[573];

    auto g_x_0_z_0_yy_zz_0_yz = buffer_1010_ddsd[574];

    auto g_x_0_z_0_yy_zz_0_zz = buffer_1010_ddsd[575];

    auto g_x_0_z_0_yz_xx_0_xx = buffer_1010_ddsd[576];

    auto g_x_0_z_0_yz_xx_0_xy = buffer_1010_ddsd[577];

    auto g_x_0_z_0_yz_xx_0_xz = buffer_1010_ddsd[578];

    auto g_x_0_z_0_yz_xx_0_yy = buffer_1010_ddsd[579];

    auto g_x_0_z_0_yz_xx_0_yz = buffer_1010_ddsd[580];

    auto g_x_0_z_0_yz_xx_0_zz = buffer_1010_ddsd[581];

    auto g_x_0_z_0_yz_xy_0_xx = buffer_1010_ddsd[582];

    auto g_x_0_z_0_yz_xy_0_xy = buffer_1010_ddsd[583];

    auto g_x_0_z_0_yz_xy_0_xz = buffer_1010_ddsd[584];

    auto g_x_0_z_0_yz_xy_0_yy = buffer_1010_ddsd[585];

    auto g_x_0_z_0_yz_xy_0_yz = buffer_1010_ddsd[586];

    auto g_x_0_z_0_yz_xy_0_zz = buffer_1010_ddsd[587];

    auto g_x_0_z_0_yz_xz_0_xx = buffer_1010_ddsd[588];

    auto g_x_0_z_0_yz_xz_0_xy = buffer_1010_ddsd[589];

    auto g_x_0_z_0_yz_xz_0_xz = buffer_1010_ddsd[590];

    auto g_x_0_z_0_yz_xz_0_yy = buffer_1010_ddsd[591];

    auto g_x_0_z_0_yz_xz_0_yz = buffer_1010_ddsd[592];

    auto g_x_0_z_0_yz_xz_0_zz = buffer_1010_ddsd[593];

    auto g_x_0_z_0_yz_yy_0_xx = buffer_1010_ddsd[594];

    auto g_x_0_z_0_yz_yy_0_xy = buffer_1010_ddsd[595];

    auto g_x_0_z_0_yz_yy_0_xz = buffer_1010_ddsd[596];

    auto g_x_0_z_0_yz_yy_0_yy = buffer_1010_ddsd[597];

    auto g_x_0_z_0_yz_yy_0_yz = buffer_1010_ddsd[598];

    auto g_x_0_z_0_yz_yy_0_zz = buffer_1010_ddsd[599];

    auto g_x_0_z_0_yz_yz_0_xx = buffer_1010_ddsd[600];

    auto g_x_0_z_0_yz_yz_0_xy = buffer_1010_ddsd[601];

    auto g_x_0_z_0_yz_yz_0_xz = buffer_1010_ddsd[602];

    auto g_x_0_z_0_yz_yz_0_yy = buffer_1010_ddsd[603];

    auto g_x_0_z_0_yz_yz_0_yz = buffer_1010_ddsd[604];

    auto g_x_0_z_0_yz_yz_0_zz = buffer_1010_ddsd[605];

    auto g_x_0_z_0_yz_zz_0_xx = buffer_1010_ddsd[606];

    auto g_x_0_z_0_yz_zz_0_xy = buffer_1010_ddsd[607];

    auto g_x_0_z_0_yz_zz_0_xz = buffer_1010_ddsd[608];

    auto g_x_0_z_0_yz_zz_0_yy = buffer_1010_ddsd[609];

    auto g_x_0_z_0_yz_zz_0_yz = buffer_1010_ddsd[610];

    auto g_x_0_z_0_yz_zz_0_zz = buffer_1010_ddsd[611];

    auto g_x_0_z_0_zz_xx_0_xx = buffer_1010_ddsd[612];

    auto g_x_0_z_0_zz_xx_0_xy = buffer_1010_ddsd[613];

    auto g_x_0_z_0_zz_xx_0_xz = buffer_1010_ddsd[614];

    auto g_x_0_z_0_zz_xx_0_yy = buffer_1010_ddsd[615];

    auto g_x_0_z_0_zz_xx_0_yz = buffer_1010_ddsd[616];

    auto g_x_0_z_0_zz_xx_0_zz = buffer_1010_ddsd[617];

    auto g_x_0_z_0_zz_xy_0_xx = buffer_1010_ddsd[618];

    auto g_x_0_z_0_zz_xy_0_xy = buffer_1010_ddsd[619];

    auto g_x_0_z_0_zz_xy_0_xz = buffer_1010_ddsd[620];

    auto g_x_0_z_0_zz_xy_0_yy = buffer_1010_ddsd[621];

    auto g_x_0_z_0_zz_xy_0_yz = buffer_1010_ddsd[622];

    auto g_x_0_z_0_zz_xy_0_zz = buffer_1010_ddsd[623];

    auto g_x_0_z_0_zz_xz_0_xx = buffer_1010_ddsd[624];

    auto g_x_0_z_0_zz_xz_0_xy = buffer_1010_ddsd[625];

    auto g_x_0_z_0_zz_xz_0_xz = buffer_1010_ddsd[626];

    auto g_x_0_z_0_zz_xz_0_yy = buffer_1010_ddsd[627];

    auto g_x_0_z_0_zz_xz_0_yz = buffer_1010_ddsd[628];

    auto g_x_0_z_0_zz_xz_0_zz = buffer_1010_ddsd[629];

    auto g_x_0_z_0_zz_yy_0_xx = buffer_1010_ddsd[630];

    auto g_x_0_z_0_zz_yy_0_xy = buffer_1010_ddsd[631];

    auto g_x_0_z_0_zz_yy_0_xz = buffer_1010_ddsd[632];

    auto g_x_0_z_0_zz_yy_0_yy = buffer_1010_ddsd[633];

    auto g_x_0_z_0_zz_yy_0_yz = buffer_1010_ddsd[634];

    auto g_x_0_z_0_zz_yy_0_zz = buffer_1010_ddsd[635];

    auto g_x_0_z_0_zz_yz_0_xx = buffer_1010_ddsd[636];

    auto g_x_0_z_0_zz_yz_0_xy = buffer_1010_ddsd[637];

    auto g_x_0_z_0_zz_yz_0_xz = buffer_1010_ddsd[638];

    auto g_x_0_z_0_zz_yz_0_yy = buffer_1010_ddsd[639];

    auto g_x_0_z_0_zz_yz_0_yz = buffer_1010_ddsd[640];

    auto g_x_0_z_0_zz_yz_0_zz = buffer_1010_ddsd[641];

    auto g_x_0_z_0_zz_zz_0_xx = buffer_1010_ddsd[642];

    auto g_x_0_z_0_zz_zz_0_xy = buffer_1010_ddsd[643];

    auto g_x_0_z_0_zz_zz_0_xz = buffer_1010_ddsd[644];

    auto g_x_0_z_0_zz_zz_0_yy = buffer_1010_ddsd[645];

    auto g_x_0_z_0_zz_zz_0_yz = buffer_1010_ddsd[646];

    auto g_x_0_z_0_zz_zz_0_zz = buffer_1010_ddsd[647];

    auto g_y_0_x_0_xx_xx_0_xx = buffer_1010_ddsd[648];

    auto g_y_0_x_0_xx_xx_0_xy = buffer_1010_ddsd[649];

    auto g_y_0_x_0_xx_xx_0_xz = buffer_1010_ddsd[650];

    auto g_y_0_x_0_xx_xx_0_yy = buffer_1010_ddsd[651];

    auto g_y_0_x_0_xx_xx_0_yz = buffer_1010_ddsd[652];

    auto g_y_0_x_0_xx_xx_0_zz = buffer_1010_ddsd[653];

    auto g_y_0_x_0_xx_xy_0_xx = buffer_1010_ddsd[654];

    auto g_y_0_x_0_xx_xy_0_xy = buffer_1010_ddsd[655];

    auto g_y_0_x_0_xx_xy_0_xz = buffer_1010_ddsd[656];

    auto g_y_0_x_0_xx_xy_0_yy = buffer_1010_ddsd[657];

    auto g_y_0_x_0_xx_xy_0_yz = buffer_1010_ddsd[658];

    auto g_y_0_x_0_xx_xy_0_zz = buffer_1010_ddsd[659];

    auto g_y_0_x_0_xx_xz_0_xx = buffer_1010_ddsd[660];

    auto g_y_0_x_0_xx_xz_0_xy = buffer_1010_ddsd[661];

    auto g_y_0_x_0_xx_xz_0_xz = buffer_1010_ddsd[662];

    auto g_y_0_x_0_xx_xz_0_yy = buffer_1010_ddsd[663];

    auto g_y_0_x_0_xx_xz_0_yz = buffer_1010_ddsd[664];

    auto g_y_0_x_0_xx_xz_0_zz = buffer_1010_ddsd[665];

    auto g_y_0_x_0_xx_yy_0_xx = buffer_1010_ddsd[666];

    auto g_y_0_x_0_xx_yy_0_xy = buffer_1010_ddsd[667];

    auto g_y_0_x_0_xx_yy_0_xz = buffer_1010_ddsd[668];

    auto g_y_0_x_0_xx_yy_0_yy = buffer_1010_ddsd[669];

    auto g_y_0_x_0_xx_yy_0_yz = buffer_1010_ddsd[670];

    auto g_y_0_x_0_xx_yy_0_zz = buffer_1010_ddsd[671];

    auto g_y_0_x_0_xx_yz_0_xx = buffer_1010_ddsd[672];

    auto g_y_0_x_0_xx_yz_0_xy = buffer_1010_ddsd[673];

    auto g_y_0_x_0_xx_yz_0_xz = buffer_1010_ddsd[674];

    auto g_y_0_x_0_xx_yz_0_yy = buffer_1010_ddsd[675];

    auto g_y_0_x_0_xx_yz_0_yz = buffer_1010_ddsd[676];

    auto g_y_0_x_0_xx_yz_0_zz = buffer_1010_ddsd[677];

    auto g_y_0_x_0_xx_zz_0_xx = buffer_1010_ddsd[678];

    auto g_y_0_x_0_xx_zz_0_xy = buffer_1010_ddsd[679];

    auto g_y_0_x_0_xx_zz_0_xz = buffer_1010_ddsd[680];

    auto g_y_0_x_0_xx_zz_0_yy = buffer_1010_ddsd[681];

    auto g_y_0_x_0_xx_zz_0_yz = buffer_1010_ddsd[682];

    auto g_y_0_x_0_xx_zz_0_zz = buffer_1010_ddsd[683];

    auto g_y_0_x_0_xy_xx_0_xx = buffer_1010_ddsd[684];

    auto g_y_0_x_0_xy_xx_0_xy = buffer_1010_ddsd[685];

    auto g_y_0_x_0_xy_xx_0_xz = buffer_1010_ddsd[686];

    auto g_y_0_x_0_xy_xx_0_yy = buffer_1010_ddsd[687];

    auto g_y_0_x_0_xy_xx_0_yz = buffer_1010_ddsd[688];

    auto g_y_0_x_0_xy_xx_0_zz = buffer_1010_ddsd[689];

    auto g_y_0_x_0_xy_xy_0_xx = buffer_1010_ddsd[690];

    auto g_y_0_x_0_xy_xy_0_xy = buffer_1010_ddsd[691];

    auto g_y_0_x_0_xy_xy_0_xz = buffer_1010_ddsd[692];

    auto g_y_0_x_0_xy_xy_0_yy = buffer_1010_ddsd[693];

    auto g_y_0_x_0_xy_xy_0_yz = buffer_1010_ddsd[694];

    auto g_y_0_x_0_xy_xy_0_zz = buffer_1010_ddsd[695];

    auto g_y_0_x_0_xy_xz_0_xx = buffer_1010_ddsd[696];

    auto g_y_0_x_0_xy_xz_0_xy = buffer_1010_ddsd[697];

    auto g_y_0_x_0_xy_xz_0_xz = buffer_1010_ddsd[698];

    auto g_y_0_x_0_xy_xz_0_yy = buffer_1010_ddsd[699];

    auto g_y_0_x_0_xy_xz_0_yz = buffer_1010_ddsd[700];

    auto g_y_0_x_0_xy_xz_0_zz = buffer_1010_ddsd[701];

    auto g_y_0_x_0_xy_yy_0_xx = buffer_1010_ddsd[702];

    auto g_y_0_x_0_xy_yy_0_xy = buffer_1010_ddsd[703];

    auto g_y_0_x_0_xy_yy_0_xz = buffer_1010_ddsd[704];

    auto g_y_0_x_0_xy_yy_0_yy = buffer_1010_ddsd[705];

    auto g_y_0_x_0_xy_yy_0_yz = buffer_1010_ddsd[706];

    auto g_y_0_x_0_xy_yy_0_zz = buffer_1010_ddsd[707];

    auto g_y_0_x_0_xy_yz_0_xx = buffer_1010_ddsd[708];

    auto g_y_0_x_0_xy_yz_0_xy = buffer_1010_ddsd[709];

    auto g_y_0_x_0_xy_yz_0_xz = buffer_1010_ddsd[710];

    auto g_y_0_x_0_xy_yz_0_yy = buffer_1010_ddsd[711];

    auto g_y_0_x_0_xy_yz_0_yz = buffer_1010_ddsd[712];

    auto g_y_0_x_0_xy_yz_0_zz = buffer_1010_ddsd[713];

    auto g_y_0_x_0_xy_zz_0_xx = buffer_1010_ddsd[714];

    auto g_y_0_x_0_xy_zz_0_xy = buffer_1010_ddsd[715];

    auto g_y_0_x_0_xy_zz_0_xz = buffer_1010_ddsd[716];

    auto g_y_0_x_0_xy_zz_0_yy = buffer_1010_ddsd[717];

    auto g_y_0_x_0_xy_zz_0_yz = buffer_1010_ddsd[718];

    auto g_y_0_x_0_xy_zz_0_zz = buffer_1010_ddsd[719];

    auto g_y_0_x_0_xz_xx_0_xx = buffer_1010_ddsd[720];

    auto g_y_0_x_0_xz_xx_0_xy = buffer_1010_ddsd[721];

    auto g_y_0_x_0_xz_xx_0_xz = buffer_1010_ddsd[722];

    auto g_y_0_x_0_xz_xx_0_yy = buffer_1010_ddsd[723];

    auto g_y_0_x_0_xz_xx_0_yz = buffer_1010_ddsd[724];

    auto g_y_0_x_0_xz_xx_0_zz = buffer_1010_ddsd[725];

    auto g_y_0_x_0_xz_xy_0_xx = buffer_1010_ddsd[726];

    auto g_y_0_x_0_xz_xy_0_xy = buffer_1010_ddsd[727];

    auto g_y_0_x_0_xz_xy_0_xz = buffer_1010_ddsd[728];

    auto g_y_0_x_0_xz_xy_0_yy = buffer_1010_ddsd[729];

    auto g_y_0_x_0_xz_xy_0_yz = buffer_1010_ddsd[730];

    auto g_y_0_x_0_xz_xy_0_zz = buffer_1010_ddsd[731];

    auto g_y_0_x_0_xz_xz_0_xx = buffer_1010_ddsd[732];

    auto g_y_0_x_0_xz_xz_0_xy = buffer_1010_ddsd[733];

    auto g_y_0_x_0_xz_xz_0_xz = buffer_1010_ddsd[734];

    auto g_y_0_x_0_xz_xz_0_yy = buffer_1010_ddsd[735];

    auto g_y_0_x_0_xz_xz_0_yz = buffer_1010_ddsd[736];

    auto g_y_0_x_0_xz_xz_0_zz = buffer_1010_ddsd[737];

    auto g_y_0_x_0_xz_yy_0_xx = buffer_1010_ddsd[738];

    auto g_y_0_x_0_xz_yy_0_xy = buffer_1010_ddsd[739];

    auto g_y_0_x_0_xz_yy_0_xz = buffer_1010_ddsd[740];

    auto g_y_0_x_0_xz_yy_0_yy = buffer_1010_ddsd[741];

    auto g_y_0_x_0_xz_yy_0_yz = buffer_1010_ddsd[742];

    auto g_y_0_x_0_xz_yy_0_zz = buffer_1010_ddsd[743];

    auto g_y_0_x_0_xz_yz_0_xx = buffer_1010_ddsd[744];

    auto g_y_0_x_0_xz_yz_0_xy = buffer_1010_ddsd[745];

    auto g_y_0_x_0_xz_yz_0_xz = buffer_1010_ddsd[746];

    auto g_y_0_x_0_xz_yz_0_yy = buffer_1010_ddsd[747];

    auto g_y_0_x_0_xz_yz_0_yz = buffer_1010_ddsd[748];

    auto g_y_0_x_0_xz_yz_0_zz = buffer_1010_ddsd[749];

    auto g_y_0_x_0_xz_zz_0_xx = buffer_1010_ddsd[750];

    auto g_y_0_x_0_xz_zz_0_xy = buffer_1010_ddsd[751];

    auto g_y_0_x_0_xz_zz_0_xz = buffer_1010_ddsd[752];

    auto g_y_0_x_0_xz_zz_0_yy = buffer_1010_ddsd[753];

    auto g_y_0_x_0_xz_zz_0_yz = buffer_1010_ddsd[754];

    auto g_y_0_x_0_xz_zz_0_zz = buffer_1010_ddsd[755];

    auto g_y_0_x_0_yy_xx_0_xx = buffer_1010_ddsd[756];

    auto g_y_0_x_0_yy_xx_0_xy = buffer_1010_ddsd[757];

    auto g_y_0_x_0_yy_xx_0_xz = buffer_1010_ddsd[758];

    auto g_y_0_x_0_yy_xx_0_yy = buffer_1010_ddsd[759];

    auto g_y_0_x_0_yy_xx_0_yz = buffer_1010_ddsd[760];

    auto g_y_0_x_0_yy_xx_0_zz = buffer_1010_ddsd[761];

    auto g_y_0_x_0_yy_xy_0_xx = buffer_1010_ddsd[762];

    auto g_y_0_x_0_yy_xy_0_xy = buffer_1010_ddsd[763];

    auto g_y_0_x_0_yy_xy_0_xz = buffer_1010_ddsd[764];

    auto g_y_0_x_0_yy_xy_0_yy = buffer_1010_ddsd[765];

    auto g_y_0_x_0_yy_xy_0_yz = buffer_1010_ddsd[766];

    auto g_y_0_x_0_yy_xy_0_zz = buffer_1010_ddsd[767];

    auto g_y_0_x_0_yy_xz_0_xx = buffer_1010_ddsd[768];

    auto g_y_0_x_0_yy_xz_0_xy = buffer_1010_ddsd[769];

    auto g_y_0_x_0_yy_xz_0_xz = buffer_1010_ddsd[770];

    auto g_y_0_x_0_yy_xz_0_yy = buffer_1010_ddsd[771];

    auto g_y_0_x_0_yy_xz_0_yz = buffer_1010_ddsd[772];

    auto g_y_0_x_0_yy_xz_0_zz = buffer_1010_ddsd[773];

    auto g_y_0_x_0_yy_yy_0_xx = buffer_1010_ddsd[774];

    auto g_y_0_x_0_yy_yy_0_xy = buffer_1010_ddsd[775];

    auto g_y_0_x_0_yy_yy_0_xz = buffer_1010_ddsd[776];

    auto g_y_0_x_0_yy_yy_0_yy = buffer_1010_ddsd[777];

    auto g_y_0_x_0_yy_yy_0_yz = buffer_1010_ddsd[778];

    auto g_y_0_x_0_yy_yy_0_zz = buffer_1010_ddsd[779];

    auto g_y_0_x_0_yy_yz_0_xx = buffer_1010_ddsd[780];

    auto g_y_0_x_0_yy_yz_0_xy = buffer_1010_ddsd[781];

    auto g_y_0_x_0_yy_yz_0_xz = buffer_1010_ddsd[782];

    auto g_y_0_x_0_yy_yz_0_yy = buffer_1010_ddsd[783];

    auto g_y_0_x_0_yy_yz_0_yz = buffer_1010_ddsd[784];

    auto g_y_0_x_0_yy_yz_0_zz = buffer_1010_ddsd[785];

    auto g_y_0_x_0_yy_zz_0_xx = buffer_1010_ddsd[786];

    auto g_y_0_x_0_yy_zz_0_xy = buffer_1010_ddsd[787];

    auto g_y_0_x_0_yy_zz_0_xz = buffer_1010_ddsd[788];

    auto g_y_0_x_0_yy_zz_0_yy = buffer_1010_ddsd[789];

    auto g_y_0_x_0_yy_zz_0_yz = buffer_1010_ddsd[790];

    auto g_y_0_x_0_yy_zz_0_zz = buffer_1010_ddsd[791];

    auto g_y_0_x_0_yz_xx_0_xx = buffer_1010_ddsd[792];

    auto g_y_0_x_0_yz_xx_0_xy = buffer_1010_ddsd[793];

    auto g_y_0_x_0_yz_xx_0_xz = buffer_1010_ddsd[794];

    auto g_y_0_x_0_yz_xx_0_yy = buffer_1010_ddsd[795];

    auto g_y_0_x_0_yz_xx_0_yz = buffer_1010_ddsd[796];

    auto g_y_0_x_0_yz_xx_0_zz = buffer_1010_ddsd[797];

    auto g_y_0_x_0_yz_xy_0_xx = buffer_1010_ddsd[798];

    auto g_y_0_x_0_yz_xy_0_xy = buffer_1010_ddsd[799];

    auto g_y_0_x_0_yz_xy_0_xz = buffer_1010_ddsd[800];

    auto g_y_0_x_0_yz_xy_0_yy = buffer_1010_ddsd[801];

    auto g_y_0_x_0_yz_xy_0_yz = buffer_1010_ddsd[802];

    auto g_y_0_x_0_yz_xy_0_zz = buffer_1010_ddsd[803];

    auto g_y_0_x_0_yz_xz_0_xx = buffer_1010_ddsd[804];

    auto g_y_0_x_0_yz_xz_0_xy = buffer_1010_ddsd[805];

    auto g_y_0_x_0_yz_xz_0_xz = buffer_1010_ddsd[806];

    auto g_y_0_x_0_yz_xz_0_yy = buffer_1010_ddsd[807];

    auto g_y_0_x_0_yz_xz_0_yz = buffer_1010_ddsd[808];

    auto g_y_0_x_0_yz_xz_0_zz = buffer_1010_ddsd[809];

    auto g_y_0_x_0_yz_yy_0_xx = buffer_1010_ddsd[810];

    auto g_y_0_x_0_yz_yy_0_xy = buffer_1010_ddsd[811];

    auto g_y_0_x_0_yz_yy_0_xz = buffer_1010_ddsd[812];

    auto g_y_0_x_0_yz_yy_0_yy = buffer_1010_ddsd[813];

    auto g_y_0_x_0_yz_yy_0_yz = buffer_1010_ddsd[814];

    auto g_y_0_x_0_yz_yy_0_zz = buffer_1010_ddsd[815];

    auto g_y_0_x_0_yz_yz_0_xx = buffer_1010_ddsd[816];

    auto g_y_0_x_0_yz_yz_0_xy = buffer_1010_ddsd[817];

    auto g_y_0_x_0_yz_yz_0_xz = buffer_1010_ddsd[818];

    auto g_y_0_x_0_yz_yz_0_yy = buffer_1010_ddsd[819];

    auto g_y_0_x_0_yz_yz_0_yz = buffer_1010_ddsd[820];

    auto g_y_0_x_0_yz_yz_0_zz = buffer_1010_ddsd[821];

    auto g_y_0_x_0_yz_zz_0_xx = buffer_1010_ddsd[822];

    auto g_y_0_x_0_yz_zz_0_xy = buffer_1010_ddsd[823];

    auto g_y_0_x_0_yz_zz_0_xz = buffer_1010_ddsd[824];

    auto g_y_0_x_0_yz_zz_0_yy = buffer_1010_ddsd[825];

    auto g_y_0_x_0_yz_zz_0_yz = buffer_1010_ddsd[826];

    auto g_y_0_x_0_yz_zz_0_zz = buffer_1010_ddsd[827];

    auto g_y_0_x_0_zz_xx_0_xx = buffer_1010_ddsd[828];

    auto g_y_0_x_0_zz_xx_0_xy = buffer_1010_ddsd[829];

    auto g_y_0_x_0_zz_xx_0_xz = buffer_1010_ddsd[830];

    auto g_y_0_x_0_zz_xx_0_yy = buffer_1010_ddsd[831];

    auto g_y_0_x_0_zz_xx_0_yz = buffer_1010_ddsd[832];

    auto g_y_0_x_0_zz_xx_0_zz = buffer_1010_ddsd[833];

    auto g_y_0_x_0_zz_xy_0_xx = buffer_1010_ddsd[834];

    auto g_y_0_x_0_zz_xy_0_xy = buffer_1010_ddsd[835];

    auto g_y_0_x_0_zz_xy_0_xz = buffer_1010_ddsd[836];

    auto g_y_0_x_0_zz_xy_0_yy = buffer_1010_ddsd[837];

    auto g_y_0_x_0_zz_xy_0_yz = buffer_1010_ddsd[838];

    auto g_y_0_x_0_zz_xy_0_zz = buffer_1010_ddsd[839];

    auto g_y_0_x_0_zz_xz_0_xx = buffer_1010_ddsd[840];

    auto g_y_0_x_0_zz_xz_0_xy = buffer_1010_ddsd[841];

    auto g_y_0_x_0_zz_xz_0_xz = buffer_1010_ddsd[842];

    auto g_y_0_x_0_zz_xz_0_yy = buffer_1010_ddsd[843];

    auto g_y_0_x_0_zz_xz_0_yz = buffer_1010_ddsd[844];

    auto g_y_0_x_0_zz_xz_0_zz = buffer_1010_ddsd[845];

    auto g_y_0_x_0_zz_yy_0_xx = buffer_1010_ddsd[846];

    auto g_y_0_x_0_zz_yy_0_xy = buffer_1010_ddsd[847];

    auto g_y_0_x_0_zz_yy_0_xz = buffer_1010_ddsd[848];

    auto g_y_0_x_0_zz_yy_0_yy = buffer_1010_ddsd[849];

    auto g_y_0_x_0_zz_yy_0_yz = buffer_1010_ddsd[850];

    auto g_y_0_x_0_zz_yy_0_zz = buffer_1010_ddsd[851];

    auto g_y_0_x_0_zz_yz_0_xx = buffer_1010_ddsd[852];

    auto g_y_0_x_0_zz_yz_0_xy = buffer_1010_ddsd[853];

    auto g_y_0_x_0_zz_yz_0_xz = buffer_1010_ddsd[854];

    auto g_y_0_x_0_zz_yz_0_yy = buffer_1010_ddsd[855];

    auto g_y_0_x_0_zz_yz_0_yz = buffer_1010_ddsd[856];

    auto g_y_0_x_0_zz_yz_0_zz = buffer_1010_ddsd[857];

    auto g_y_0_x_0_zz_zz_0_xx = buffer_1010_ddsd[858];

    auto g_y_0_x_0_zz_zz_0_xy = buffer_1010_ddsd[859];

    auto g_y_0_x_0_zz_zz_0_xz = buffer_1010_ddsd[860];

    auto g_y_0_x_0_zz_zz_0_yy = buffer_1010_ddsd[861];

    auto g_y_0_x_0_zz_zz_0_yz = buffer_1010_ddsd[862];

    auto g_y_0_x_0_zz_zz_0_zz = buffer_1010_ddsd[863];

    auto g_y_0_y_0_xx_xx_0_xx = buffer_1010_ddsd[864];

    auto g_y_0_y_0_xx_xx_0_xy = buffer_1010_ddsd[865];

    auto g_y_0_y_0_xx_xx_0_xz = buffer_1010_ddsd[866];

    auto g_y_0_y_0_xx_xx_0_yy = buffer_1010_ddsd[867];

    auto g_y_0_y_0_xx_xx_0_yz = buffer_1010_ddsd[868];

    auto g_y_0_y_0_xx_xx_0_zz = buffer_1010_ddsd[869];

    auto g_y_0_y_0_xx_xy_0_xx = buffer_1010_ddsd[870];

    auto g_y_0_y_0_xx_xy_0_xy = buffer_1010_ddsd[871];

    auto g_y_0_y_0_xx_xy_0_xz = buffer_1010_ddsd[872];

    auto g_y_0_y_0_xx_xy_0_yy = buffer_1010_ddsd[873];

    auto g_y_0_y_0_xx_xy_0_yz = buffer_1010_ddsd[874];

    auto g_y_0_y_0_xx_xy_0_zz = buffer_1010_ddsd[875];

    auto g_y_0_y_0_xx_xz_0_xx = buffer_1010_ddsd[876];

    auto g_y_0_y_0_xx_xz_0_xy = buffer_1010_ddsd[877];

    auto g_y_0_y_0_xx_xz_0_xz = buffer_1010_ddsd[878];

    auto g_y_0_y_0_xx_xz_0_yy = buffer_1010_ddsd[879];

    auto g_y_0_y_0_xx_xz_0_yz = buffer_1010_ddsd[880];

    auto g_y_0_y_0_xx_xz_0_zz = buffer_1010_ddsd[881];

    auto g_y_0_y_0_xx_yy_0_xx = buffer_1010_ddsd[882];

    auto g_y_0_y_0_xx_yy_0_xy = buffer_1010_ddsd[883];

    auto g_y_0_y_0_xx_yy_0_xz = buffer_1010_ddsd[884];

    auto g_y_0_y_0_xx_yy_0_yy = buffer_1010_ddsd[885];

    auto g_y_0_y_0_xx_yy_0_yz = buffer_1010_ddsd[886];

    auto g_y_0_y_0_xx_yy_0_zz = buffer_1010_ddsd[887];

    auto g_y_0_y_0_xx_yz_0_xx = buffer_1010_ddsd[888];

    auto g_y_0_y_0_xx_yz_0_xy = buffer_1010_ddsd[889];

    auto g_y_0_y_0_xx_yz_0_xz = buffer_1010_ddsd[890];

    auto g_y_0_y_0_xx_yz_0_yy = buffer_1010_ddsd[891];

    auto g_y_0_y_0_xx_yz_0_yz = buffer_1010_ddsd[892];

    auto g_y_0_y_0_xx_yz_0_zz = buffer_1010_ddsd[893];

    auto g_y_0_y_0_xx_zz_0_xx = buffer_1010_ddsd[894];

    auto g_y_0_y_0_xx_zz_0_xy = buffer_1010_ddsd[895];

    auto g_y_0_y_0_xx_zz_0_xz = buffer_1010_ddsd[896];

    auto g_y_0_y_0_xx_zz_0_yy = buffer_1010_ddsd[897];

    auto g_y_0_y_0_xx_zz_0_yz = buffer_1010_ddsd[898];

    auto g_y_0_y_0_xx_zz_0_zz = buffer_1010_ddsd[899];

    auto g_y_0_y_0_xy_xx_0_xx = buffer_1010_ddsd[900];

    auto g_y_0_y_0_xy_xx_0_xy = buffer_1010_ddsd[901];

    auto g_y_0_y_0_xy_xx_0_xz = buffer_1010_ddsd[902];

    auto g_y_0_y_0_xy_xx_0_yy = buffer_1010_ddsd[903];

    auto g_y_0_y_0_xy_xx_0_yz = buffer_1010_ddsd[904];

    auto g_y_0_y_0_xy_xx_0_zz = buffer_1010_ddsd[905];

    auto g_y_0_y_0_xy_xy_0_xx = buffer_1010_ddsd[906];

    auto g_y_0_y_0_xy_xy_0_xy = buffer_1010_ddsd[907];

    auto g_y_0_y_0_xy_xy_0_xz = buffer_1010_ddsd[908];

    auto g_y_0_y_0_xy_xy_0_yy = buffer_1010_ddsd[909];

    auto g_y_0_y_0_xy_xy_0_yz = buffer_1010_ddsd[910];

    auto g_y_0_y_0_xy_xy_0_zz = buffer_1010_ddsd[911];

    auto g_y_0_y_0_xy_xz_0_xx = buffer_1010_ddsd[912];

    auto g_y_0_y_0_xy_xz_0_xy = buffer_1010_ddsd[913];

    auto g_y_0_y_0_xy_xz_0_xz = buffer_1010_ddsd[914];

    auto g_y_0_y_0_xy_xz_0_yy = buffer_1010_ddsd[915];

    auto g_y_0_y_0_xy_xz_0_yz = buffer_1010_ddsd[916];

    auto g_y_0_y_0_xy_xz_0_zz = buffer_1010_ddsd[917];

    auto g_y_0_y_0_xy_yy_0_xx = buffer_1010_ddsd[918];

    auto g_y_0_y_0_xy_yy_0_xy = buffer_1010_ddsd[919];

    auto g_y_0_y_0_xy_yy_0_xz = buffer_1010_ddsd[920];

    auto g_y_0_y_0_xy_yy_0_yy = buffer_1010_ddsd[921];

    auto g_y_0_y_0_xy_yy_0_yz = buffer_1010_ddsd[922];

    auto g_y_0_y_0_xy_yy_0_zz = buffer_1010_ddsd[923];

    auto g_y_0_y_0_xy_yz_0_xx = buffer_1010_ddsd[924];

    auto g_y_0_y_0_xy_yz_0_xy = buffer_1010_ddsd[925];

    auto g_y_0_y_0_xy_yz_0_xz = buffer_1010_ddsd[926];

    auto g_y_0_y_0_xy_yz_0_yy = buffer_1010_ddsd[927];

    auto g_y_0_y_0_xy_yz_0_yz = buffer_1010_ddsd[928];

    auto g_y_0_y_0_xy_yz_0_zz = buffer_1010_ddsd[929];

    auto g_y_0_y_0_xy_zz_0_xx = buffer_1010_ddsd[930];

    auto g_y_0_y_0_xy_zz_0_xy = buffer_1010_ddsd[931];

    auto g_y_0_y_0_xy_zz_0_xz = buffer_1010_ddsd[932];

    auto g_y_0_y_0_xy_zz_0_yy = buffer_1010_ddsd[933];

    auto g_y_0_y_0_xy_zz_0_yz = buffer_1010_ddsd[934];

    auto g_y_0_y_0_xy_zz_0_zz = buffer_1010_ddsd[935];

    auto g_y_0_y_0_xz_xx_0_xx = buffer_1010_ddsd[936];

    auto g_y_0_y_0_xz_xx_0_xy = buffer_1010_ddsd[937];

    auto g_y_0_y_0_xz_xx_0_xz = buffer_1010_ddsd[938];

    auto g_y_0_y_0_xz_xx_0_yy = buffer_1010_ddsd[939];

    auto g_y_0_y_0_xz_xx_0_yz = buffer_1010_ddsd[940];

    auto g_y_0_y_0_xz_xx_0_zz = buffer_1010_ddsd[941];

    auto g_y_0_y_0_xz_xy_0_xx = buffer_1010_ddsd[942];

    auto g_y_0_y_0_xz_xy_0_xy = buffer_1010_ddsd[943];

    auto g_y_0_y_0_xz_xy_0_xz = buffer_1010_ddsd[944];

    auto g_y_0_y_0_xz_xy_0_yy = buffer_1010_ddsd[945];

    auto g_y_0_y_0_xz_xy_0_yz = buffer_1010_ddsd[946];

    auto g_y_0_y_0_xz_xy_0_zz = buffer_1010_ddsd[947];

    auto g_y_0_y_0_xz_xz_0_xx = buffer_1010_ddsd[948];

    auto g_y_0_y_0_xz_xz_0_xy = buffer_1010_ddsd[949];

    auto g_y_0_y_0_xz_xz_0_xz = buffer_1010_ddsd[950];

    auto g_y_0_y_0_xz_xz_0_yy = buffer_1010_ddsd[951];

    auto g_y_0_y_0_xz_xz_0_yz = buffer_1010_ddsd[952];

    auto g_y_0_y_0_xz_xz_0_zz = buffer_1010_ddsd[953];

    auto g_y_0_y_0_xz_yy_0_xx = buffer_1010_ddsd[954];

    auto g_y_0_y_0_xz_yy_0_xy = buffer_1010_ddsd[955];

    auto g_y_0_y_0_xz_yy_0_xz = buffer_1010_ddsd[956];

    auto g_y_0_y_0_xz_yy_0_yy = buffer_1010_ddsd[957];

    auto g_y_0_y_0_xz_yy_0_yz = buffer_1010_ddsd[958];

    auto g_y_0_y_0_xz_yy_0_zz = buffer_1010_ddsd[959];

    auto g_y_0_y_0_xz_yz_0_xx = buffer_1010_ddsd[960];

    auto g_y_0_y_0_xz_yz_0_xy = buffer_1010_ddsd[961];

    auto g_y_0_y_0_xz_yz_0_xz = buffer_1010_ddsd[962];

    auto g_y_0_y_0_xz_yz_0_yy = buffer_1010_ddsd[963];

    auto g_y_0_y_0_xz_yz_0_yz = buffer_1010_ddsd[964];

    auto g_y_0_y_0_xz_yz_0_zz = buffer_1010_ddsd[965];

    auto g_y_0_y_0_xz_zz_0_xx = buffer_1010_ddsd[966];

    auto g_y_0_y_0_xz_zz_0_xy = buffer_1010_ddsd[967];

    auto g_y_0_y_0_xz_zz_0_xz = buffer_1010_ddsd[968];

    auto g_y_0_y_0_xz_zz_0_yy = buffer_1010_ddsd[969];

    auto g_y_0_y_0_xz_zz_0_yz = buffer_1010_ddsd[970];

    auto g_y_0_y_0_xz_zz_0_zz = buffer_1010_ddsd[971];

    auto g_y_0_y_0_yy_xx_0_xx = buffer_1010_ddsd[972];

    auto g_y_0_y_0_yy_xx_0_xy = buffer_1010_ddsd[973];

    auto g_y_0_y_0_yy_xx_0_xz = buffer_1010_ddsd[974];

    auto g_y_0_y_0_yy_xx_0_yy = buffer_1010_ddsd[975];

    auto g_y_0_y_0_yy_xx_0_yz = buffer_1010_ddsd[976];

    auto g_y_0_y_0_yy_xx_0_zz = buffer_1010_ddsd[977];

    auto g_y_0_y_0_yy_xy_0_xx = buffer_1010_ddsd[978];

    auto g_y_0_y_0_yy_xy_0_xy = buffer_1010_ddsd[979];

    auto g_y_0_y_0_yy_xy_0_xz = buffer_1010_ddsd[980];

    auto g_y_0_y_0_yy_xy_0_yy = buffer_1010_ddsd[981];

    auto g_y_0_y_0_yy_xy_0_yz = buffer_1010_ddsd[982];

    auto g_y_0_y_0_yy_xy_0_zz = buffer_1010_ddsd[983];

    auto g_y_0_y_0_yy_xz_0_xx = buffer_1010_ddsd[984];

    auto g_y_0_y_0_yy_xz_0_xy = buffer_1010_ddsd[985];

    auto g_y_0_y_0_yy_xz_0_xz = buffer_1010_ddsd[986];

    auto g_y_0_y_0_yy_xz_0_yy = buffer_1010_ddsd[987];

    auto g_y_0_y_0_yy_xz_0_yz = buffer_1010_ddsd[988];

    auto g_y_0_y_0_yy_xz_0_zz = buffer_1010_ddsd[989];

    auto g_y_0_y_0_yy_yy_0_xx = buffer_1010_ddsd[990];

    auto g_y_0_y_0_yy_yy_0_xy = buffer_1010_ddsd[991];

    auto g_y_0_y_0_yy_yy_0_xz = buffer_1010_ddsd[992];

    auto g_y_0_y_0_yy_yy_0_yy = buffer_1010_ddsd[993];

    auto g_y_0_y_0_yy_yy_0_yz = buffer_1010_ddsd[994];

    auto g_y_0_y_0_yy_yy_0_zz = buffer_1010_ddsd[995];

    auto g_y_0_y_0_yy_yz_0_xx = buffer_1010_ddsd[996];

    auto g_y_0_y_0_yy_yz_0_xy = buffer_1010_ddsd[997];

    auto g_y_0_y_0_yy_yz_0_xz = buffer_1010_ddsd[998];

    auto g_y_0_y_0_yy_yz_0_yy = buffer_1010_ddsd[999];

    auto g_y_0_y_0_yy_yz_0_yz = buffer_1010_ddsd[1000];

    auto g_y_0_y_0_yy_yz_0_zz = buffer_1010_ddsd[1001];

    auto g_y_0_y_0_yy_zz_0_xx = buffer_1010_ddsd[1002];

    auto g_y_0_y_0_yy_zz_0_xy = buffer_1010_ddsd[1003];

    auto g_y_0_y_0_yy_zz_0_xz = buffer_1010_ddsd[1004];

    auto g_y_0_y_0_yy_zz_0_yy = buffer_1010_ddsd[1005];

    auto g_y_0_y_0_yy_zz_0_yz = buffer_1010_ddsd[1006];

    auto g_y_0_y_0_yy_zz_0_zz = buffer_1010_ddsd[1007];

    auto g_y_0_y_0_yz_xx_0_xx = buffer_1010_ddsd[1008];

    auto g_y_0_y_0_yz_xx_0_xy = buffer_1010_ddsd[1009];

    auto g_y_0_y_0_yz_xx_0_xz = buffer_1010_ddsd[1010];

    auto g_y_0_y_0_yz_xx_0_yy = buffer_1010_ddsd[1011];

    auto g_y_0_y_0_yz_xx_0_yz = buffer_1010_ddsd[1012];

    auto g_y_0_y_0_yz_xx_0_zz = buffer_1010_ddsd[1013];

    auto g_y_0_y_0_yz_xy_0_xx = buffer_1010_ddsd[1014];

    auto g_y_0_y_0_yz_xy_0_xy = buffer_1010_ddsd[1015];

    auto g_y_0_y_0_yz_xy_0_xz = buffer_1010_ddsd[1016];

    auto g_y_0_y_0_yz_xy_0_yy = buffer_1010_ddsd[1017];

    auto g_y_0_y_0_yz_xy_0_yz = buffer_1010_ddsd[1018];

    auto g_y_0_y_0_yz_xy_0_zz = buffer_1010_ddsd[1019];

    auto g_y_0_y_0_yz_xz_0_xx = buffer_1010_ddsd[1020];

    auto g_y_0_y_0_yz_xz_0_xy = buffer_1010_ddsd[1021];

    auto g_y_0_y_0_yz_xz_0_xz = buffer_1010_ddsd[1022];

    auto g_y_0_y_0_yz_xz_0_yy = buffer_1010_ddsd[1023];

    auto g_y_0_y_0_yz_xz_0_yz = buffer_1010_ddsd[1024];

    auto g_y_0_y_0_yz_xz_0_zz = buffer_1010_ddsd[1025];

    auto g_y_0_y_0_yz_yy_0_xx = buffer_1010_ddsd[1026];

    auto g_y_0_y_0_yz_yy_0_xy = buffer_1010_ddsd[1027];

    auto g_y_0_y_0_yz_yy_0_xz = buffer_1010_ddsd[1028];

    auto g_y_0_y_0_yz_yy_0_yy = buffer_1010_ddsd[1029];

    auto g_y_0_y_0_yz_yy_0_yz = buffer_1010_ddsd[1030];

    auto g_y_0_y_0_yz_yy_0_zz = buffer_1010_ddsd[1031];

    auto g_y_0_y_0_yz_yz_0_xx = buffer_1010_ddsd[1032];

    auto g_y_0_y_0_yz_yz_0_xy = buffer_1010_ddsd[1033];

    auto g_y_0_y_0_yz_yz_0_xz = buffer_1010_ddsd[1034];

    auto g_y_0_y_0_yz_yz_0_yy = buffer_1010_ddsd[1035];

    auto g_y_0_y_0_yz_yz_0_yz = buffer_1010_ddsd[1036];

    auto g_y_0_y_0_yz_yz_0_zz = buffer_1010_ddsd[1037];

    auto g_y_0_y_0_yz_zz_0_xx = buffer_1010_ddsd[1038];

    auto g_y_0_y_0_yz_zz_0_xy = buffer_1010_ddsd[1039];

    auto g_y_0_y_0_yz_zz_0_xz = buffer_1010_ddsd[1040];

    auto g_y_0_y_0_yz_zz_0_yy = buffer_1010_ddsd[1041];

    auto g_y_0_y_0_yz_zz_0_yz = buffer_1010_ddsd[1042];

    auto g_y_0_y_0_yz_zz_0_zz = buffer_1010_ddsd[1043];

    auto g_y_0_y_0_zz_xx_0_xx = buffer_1010_ddsd[1044];

    auto g_y_0_y_0_zz_xx_0_xy = buffer_1010_ddsd[1045];

    auto g_y_0_y_0_zz_xx_0_xz = buffer_1010_ddsd[1046];

    auto g_y_0_y_0_zz_xx_0_yy = buffer_1010_ddsd[1047];

    auto g_y_0_y_0_zz_xx_0_yz = buffer_1010_ddsd[1048];

    auto g_y_0_y_0_zz_xx_0_zz = buffer_1010_ddsd[1049];

    auto g_y_0_y_0_zz_xy_0_xx = buffer_1010_ddsd[1050];

    auto g_y_0_y_0_zz_xy_0_xy = buffer_1010_ddsd[1051];

    auto g_y_0_y_0_zz_xy_0_xz = buffer_1010_ddsd[1052];

    auto g_y_0_y_0_zz_xy_0_yy = buffer_1010_ddsd[1053];

    auto g_y_0_y_0_zz_xy_0_yz = buffer_1010_ddsd[1054];

    auto g_y_0_y_0_zz_xy_0_zz = buffer_1010_ddsd[1055];

    auto g_y_0_y_0_zz_xz_0_xx = buffer_1010_ddsd[1056];

    auto g_y_0_y_0_zz_xz_0_xy = buffer_1010_ddsd[1057];

    auto g_y_0_y_0_zz_xz_0_xz = buffer_1010_ddsd[1058];

    auto g_y_0_y_0_zz_xz_0_yy = buffer_1010_ddsd[1059];

    auto g_y_0_y_0_zz_xz_0_yz = buffer_1010_ddsd[1060];

    auto g_y_0_y_0_zz_xz_0_zz = buffer_1010_ddsd[1061];

    auto g_y_0_y_0_zz_yy_0_xx = buffer_1010_ddsd[1062];

    auto g_y_0_y_0_zz_yy_0_xy = buffer_1010_ddsd[1063];

    auto g_y_0_y_0_zz_yy_0_xz = buffer_1010_ddsd[1064];

    auto g_y_0_y_0_zz_yy_0_yy = buffer_1010_ddsd[1065];

    auto g_y_0_y_0_zz_yy_0_yz = buffer_1010_ddsd[1066];

    auto g_y_0_y_0_zz_yy_0_zz = buffer_1010_ddsd[1067];

    auto g_y_0_y_0_zz_yz_0_xx = buffer_1010_ddsd[1068];

    auto g_y_0_y_0_zz_yz_0_xy = buffer_1010_ddsd[1069];

    auto g_y_0_y_0_zz_yz_0_xz = buffer_1010_ddsd[1070];

    auto g_y_0_y_0_zz_yz_0_yy = buffer_1010_ddsd[1071];

    auto g_y_0_y_0_zz_yz_0_yz = buffer_1010_ddsd[1072];

    auto g_y_0_y_0_zz_yz_0_zz = buffer_1010_ddsd[1073];

    auto g_y_0_y_0_zz_zz_0_xx = buffer_1010_ddsd[1074];

    auto g_y_0_y_0_zz_zz_0_xy = buffer_1010_ddsd[1075];

    auto g_y_0_y_0_zz_zz_0_xz = buffer_1010_ddsd[1076];

    auto g_y_0_y_0_zz_zz_0_yy = buffer_1010_ddsd[1077];

    auto g_y_0_y_0_zz_zz_0_yz = buffer_1010_ddsd[1078];

    auto g_y_0_y_0_zz_zz_0_zz = buffer_1010_ddsd[1079];

    auto g_y_0_z_0_xx_xx_0_xx = buffer_1010_ddsd[1080];

    auto g_y_0_z_0_xx_xx_0_xy = buffer_1010_ddsd[1081];

    auto g_y_0_z_0_xx_xx_0_xz = buffer_1010_ddsd[1082];

    auto g_y_0_z_0_xx_xx_0_yy = buffer_1010_ddsd[1083];

    auto g_y_0_z_0_xx_xx_0_yz = buffer_1010_ddsd[1084];

    auto g_y_0_z_0_xx_xx_0_zz = buffer_1010_ddsd[1085];

    auto g_y_0_z_0_xx_xy_0_xx = buffer_1010_ddsd[1086];

    auto g_y_0_z_0_xx_xy_0_xy = buffer_1010_ddsd[1087];

    auto g_y_0_z_0_xx_xy_0_xz = buffer_1010_ddsd[1088];

    auto g_y_0_z_0_xx_xy_0_yy = buffer_1010_ddsd[1089];

    auto g_y_0_z_0_xx_xy_0_yz = buffer_1010_ddsd[1090];

    auto g_y_0_z_0_xx_xy_0_zz = buffer_1010_ddsd[1091];

    auto g_y_0_z_0_xx_xz_0_xx = buffer_1010_ddsd[1092];

    auto g_y_0_z_0_xx_xz_0_xy = buffer_1010_ddsd[1093];

    auto g_y_0_z_0_xx_xz_0_xz = buffer_1010_ddsd[1094];

    auto g_y_0_z_0_xx_xz_0_yy = buffer_1010_ddsd[1095];

    auto g_y_0_z_0_xx_xz_0_yz = buffer_1010_ddsd[1096];

    auto g_y_0_z_0_xx_xz_0_zz = buffer_1010_ddsd[1097];

    auto g_y_0_z_0_xx_yy_0_xx = buffer_1010_ddsd[1098];

    auto g_y_0_z_0_xx_yy_0_xy = buffer_1010_ddsd[1099];

    auto g_y_0_z_0_xx_yy_0_xz = buffer_1010_ddsd[1100];

    auto g_y_0_z_0_xx_yy_0_yy = buffer_1010_ddsd[1101];

    auto g_y_0_z_0_xx_yy_0_yz = buffer_1010_ddsd[1102];

    auto g_y_0_z_0_xx_yy_0_zz = buffer_1010_ddsd[1103];

    auto g_y_0_z_0_xx_yz_0_xx = buffer_1010_ddsd[1104];

    auto g_y_0_z_0_xx_yz_0_xy = buffer_1010_ddsd[1105];

    auto g_y_0_z_0_xx_yz_0_xz = buffer_1010_ddsd[1106];

    auto g_y_0_z_0_xx_yz_0_yy = buffer_1010_ddsd[1107];

    auto g_y_0_z_0_xx_yz_0_yz = buffer_1010_ddsd[1108];

    auto g_y_0_z_0_xx_yz_0_zz = buffer_1010_ddsd[1109];

    auto g_y_0_z_0_xx_zz_0_xx = buffer_1010_ddsd[1110];

    auto g_y_0_z_0_xx_zz_0_xy = buffer_1010_ddsd[1111];

    auto g_y_0_z_0_xx_zz_0_xz = buffer_1010_ddsd[1112];

    auto g_y_0_z_0_xx_zz_0_yy = buffer_1010_ddsd[1113];

    auto g_y_0_z_0_xx_zz_0_yz = buffer_1010_ddsd[1114];

    auto g_y_0_z_0_xx_zz_0_zz = buffer_1010_ddsd[1115];

    auto g_y_0_z_0_xy_xx_0_xx = buffer_1010_ddsd[1116];

    auto g_y_0_z_0_xy_xx_0_xy = buffer_1010_ddsd[1117];

    auto g_y_0_z_0_xy_xx_0_xz = buffer_1010_ddsd[1118];

    auto g_y_0_z_0_xy_xx_0_yy = buffer_1010_ddsd[1119];

    auto g_y_0_z_0_xy_xx_0_yz = buffer_1010_ddsd[1120];

    auto g_y_0_z_0_xy_xx_0_zz = buffer_1010_ddsd[1121];

    auto g_y_0_z_0_xy_xy_0_xx = buffer_1010_ddsd[1122];

    auto g_y_0_z_0_xy_xy_0_xy = buffer_1010_ddsd[1123];

    auto g_y_0_z_0_xy_xy_0_xz = buffer_1010_ddsd[1124];

    auto g_y_0_z_0_xy_xy_0_yy = buffer_1010_ddsd[1125];

    auto g_y_0_z_0_xy_xy_0_yz = buffer_1010_ddsd[1126];

    auto g_y_0_z_0_xy_xy_0_zz = buffer_1010_ddsd[1127];

    auto g_y_0_z_0_xy_xz_0_xx = buffer_1010_ddsd[1128];

    auto g_y_0_z_0_xy_xz_0_xy = buffer_1010_ddsd[1129];

    auto g_y_0_z_0_xy_xz_0_xz = buffer_1010_ddsd[1130];

    auto g_y_0_z_0_xy_xz_0_yy = buffer_1010_ddsd[1131];

    auto g_y_0_z_0_xy_xz_0_yz = buffer_1010_ddsd[1132];

    auto g_y_0_z_0_xy_xz_0_zz = buffer_1010_ddsd[1133];

    auto g_y_0_z_0_xy_yy_0_xx = buffer_1010_ddsd[1134];

    auto g_y_0_z_0_xy_yy_0_xy = buffer_1010_ddsd[1135];

    auto g_y_0_z_0_xy_yy_0_xz = buffer_1010_ddsd[1136];

    auto g_y_0_z_0_xy_yy_0_yy = buffer_1010_ddsd[1137];

    auto g_y_0_z_0_xy_yy_0_yz = buffer_1010_ddsd[1138];

    auto g_y_0_z_0_xy_yy_0_zz = buffer_1010_ddsd[1139];

    auto g_y_0_z_0_xy_yz_0_xx = buffer_1010_ddsd[1140];

    auto g_y_0_z_0_xy_yz_0_xy = buffer_1010_ddsd[1141];

    auto g_y_0_z_0_xy_yz_0_xz = buffer_1010_ddsd[1142];

    auto g_y_0_z_0_xy_yz_0_yy = buffer_1010_ddsd[1143];

    auto g_y_0_z_0_xy_yz_0_yz = buffer_1010_ddsd[1144];

    auto g_y_0_z_0_xy_yz_0_zz = buffer_1010_ddsd[1145];

    auto g_y_0_z_0_xy_zz_0_xx = buffer_1010_ddsd[1146];

    auto g_y_0_z_0_xy_zz_0_xy = buffer_1010_ddsd[1147];

    auto g_y_0_z_0_xy_zz_0_xz = buffer_1010_ddsd[1148];

    auto g_y_0_z_0_xy_zz_0_yy = buffer_1010_ddsd[1149];

    auto g_y_0_z_0_xy_zz_0_yz = buffer_1010_ddsd[1150];

    auto g_y_0_z_0_xy_zz_0_zz = buffer_1010_ddsd[1151];

    auto g_y_0_z_0_xz_xx_0_xx = buffer_1010_ddsd[1152];

    auto g_y_0_z_0_xz_xx_0_xy = buffer_1010_ddsd[1153];

    auto g_y_0_z_0_xz_xx_0_xz = buffer_1010_ddsd[1154];

    auto g_y_0_z_0_xz_xx_0_yy = buffer_1010_ddsd[1155];

    auto g_y_0_z_0_xz_xx_0_yz = buffer_1010_ddsd[1156];

    auto g_y_0_z_0_xz_xx_0_zz = buffer_1010_ddsd[1157];

    auto g_y_0_z_0_xz_xy_0_xx = buffer_1010_ddsd[1158];

    auto g_y_0_z_0_xz_xy_0_xy = buffer_1010_ddsd[1159];

    auto g_y_0_z_0_xz_xy_0_xz = buffer_1010_ddsd[1160];

    auto g_y_0_z_0_xz_xy_0_yy = buffer_1010_ddsd[1161];

    auto g_y_0_z_0_xz_xy_0_yz = buffer_1010_ddsd[1162];

    auto g_y_0_z_0_xz_xy_0_zz = buffer_1010_ddsd[1163];

    auto g_y_0_z_0_xz_xz_0_xx = buffer_1010_ddsd[1164];

    auto g_y_0_z_0_xz_xz_0_xy = buffer_1010_ddsd[1165];

    auto g_y_0_z_0_xz_xz_0_xz = buffer_1010_ddsd[1166];

    auto g_y_0_z_0_xz_xz_0_yy = buffer_1010_ddsd[1167];

    auto g_y_0_z_0_xz_xz_0_yz = buffer_1010_ddsd[1168];

    auto g_y_0_z_0_xz_xz_0_zz = buffer_1010_ddsd[1169];

    auto g_y_0_z_0_xz_yy_0_xx = buffer_1010_ddsd[1170];

    auto g_y_0_z_0_xz_yy_0_xy = buffer_1010_ddsd[1171];

    auto g_y_0_z_0_xz_yy_0_xz = buffer_1010_ddsd[1172];

    auto g_y_0_z_0_xz_yy_0_yy = buffer_1010_ddsd[1173];

    auto g_y_0_z_0_xz_yy_0_yz = buffer_1010_ddsd[1174];

    auto g_y_0_z_0_xz_yy_0_zz = buffer_1010_ddsd[1175];

    auto g_y_0_z_0_xz_yz_0_xx = buffer_1010_ddsd[1176];

    auto g_y_0_z_0_xz_yz_0_xy = buffer_1010_ddsd[1177];

    auto g_y_0_z_0_xz_yz_0_xz = buffer_1010_ddsd[1178];

    auto g_y_0_z_0_xz_yz_0_yy = buffer_1010_ddsd[1179];

    auto g_y_0_z_0_xz_yz_0_yz = buffer_1010_ddsd[1180];

    auto g_y_0_z_0_xz_yz_0_zz = buffer_1010_ddsd[1181];

    auto g_y_0_z_0_xz_zz_0_xx = buffer_1010_ddsd[1182];

    auto g_y_0_z_0_xz_zz_0_xy = buffer_1010_ddsd[1183];

    auto g_y_0_z_0_xz_zz_0_xz = buffer_1010_ddsd[1184];

    auto g_y_0_z_0_xz_zz_0_yy = buffer_1010_ddsd[1185];

    auto g_y_0_z_0_xz_zz_0_yz = buffer_1010_ddsd[1186];

    auto g_y_0_z_0_xz_zz_0_zz = buffer_1010_ddsd[1187];

    auto g_y_0_z_0_yy_xx_0_xx = buffer_1010_ddsd[1188];

    auto g_y_0_z_0_yy_xx_0_xy = buffer_1010_ddsd[1189];

    auto g_y_0_z_0_yy_xx_0_xz = buffer_1010_ddsd[1190];

    auto g_y_0_z_0_yy_xx_0_yy = buffer_1010_ddsd[1191];

    auto g_y_0_z_0_yy_xx_0_yz = buffer_1010_ddsd[1192];

    auto g_y_0_z_0_yy_xx_0_zz = buffer_1010_ddsd[1193];

    auto g_y_0_z_0_yy_xy_0_xx = buffer_1010_ddsd[1194];

    auto g_y_0_z_0_yy_xy_0_xy = buffer_1010_ddsd[1195];

    auto g_y_0_z_0_yy_xy_0_xz = buffer_1010_ddsd[1196];

    auto g_y_0_z_0_yy_xy_0_yy = buffer_1010_ddsd[1197];

    auto g_y_0_z_0_yy_xy_0_yz = buffer_1010_ddsd[1198];

    auto g_y_0_z_0_yy_xy_0_zz = buffer_1010_ddsd[1199];

    auto g_y_0_z_0_yy_xz_0_xx = buffer_1010_ddsd[1200];

    auto g_y_0_z_0_yy_xz_0_xy = buffer_1010_ddsd[1201];

    auto g_y_0_z_0_yy_xz_0_xz = buffer_1010_ddsd[1202];

    auto g_y_0_z_0_yy_xz_0_yy = buffer_1010_ddsd[1203];

    auto g_y_0_z_0_yy_xz_0_yz = buffer_1010_ddsd[1204];

    auto g_y_0_z_0_yy_xz_0_zz = buffer_1010_ddsd[1205];

    auto g_y_0_z_0_yy_yy_0_xx = buffer_1010_ddsd[1206];

    auto g_y_0_z_0_yy_yy_0_xy = buffer_1010_ddsd[1207];

    auto g_y_0_z_0_yy_yy_0_xz = buffer_1010_ddsd[1208];

    auto g_y_0_z_0_yy_yy_0_yy = buffer_1010_ddsd[1209];

    auto g_y_0_z_0_yy_yy_0_yz = buffer_1010_ddsd[1210];

    auto g_y_0_z_0_yy_yy_0_zz = buffer_1010_ddsd[1211];

    auto g_y_0_z_0_yy_yz_0_xx = buffer_1010_ddsd[1212];

    auto g_y_0_z_0_yy_yz_0_xy = buffer_1010_ddsd[1213];

    auto g_y_0_z_0_yy_yz_0_xz = buffer_1010_ddsd[1214];

    auto g_y_0_z_0_yy_yz_0_yy = buffer_1010_ddsd[1215];

    auto g_y_0_z_0_yy_yz_0_yz = buffer_1010_ddsd[1216];

    auto g_y_0_z_0_yy_yz_0_zz = buffer_1010_ddsd[1217];

    auto g_y_0_z_0_yy_zz_0_xx = buffer_1010_ddsd[1218];

    auto g_y_0_z_0_yy_zz_0_xy = buffer_1010_ddsd[1219];

    auto g_y_0_z_0_yy_zz_0_xz = buffer_1010_ddsd[1220];

    auto g_y_0_z_0_yy_zz_0_yy = buffer_1010_ddsd[1221];

    auto g_y_0_z_0_yy_zz_0_yz = buffer_1010_ddsd[1222];

    auto g_y_0_z_0_yy_zz_0_zz = buffer_1010_ddsd[1223];

    auto g_y_0_z_0_yz_xx_0_xx = buffer_1010_ddsd[1224];

    auto g_y_0_z_0_yz_xx_0_xy = buffer_1010_ddsd[1225];

    auto g_y_0_z_0_yz_xx_0_xz = buffer_1010_ddsd[1226];

    auto g_y_0_z_0_yz_xx_0_yy = buffer_1010_ddsd[1227];

    auto g_y_0_z_0_yz_xx_0_yz = buffer_1010_ddsd[1228];

    auto g_y_0_z_0_yz_xx_0_zz = buffer_1010_ddsd[1229];

    auto g_y_0_z_0_yz_xy_0_xx = buffer_1010_ddsd[1230];

    auto g_y_0_z_0_yz_xy_0_xy = buffer_1010_ddsd[1231];

    auto g_y_0_z_0_yz_xy_0_xz = buffer_1010_ddsd[1232];

    auto g_y_0_z_0_yz_xy_0_yy = buffer_1010_ddsd[1233];

    auto g_y_0_z_0_yz_xy_0_yz = buffer_1010_ddsd[1234];

    auto g_y_0_z_0_yz_xy_0_zz = buffer_1010_ddsd[1235];

    auto g_y_0_z_0_yz_xz_0_xx = buffer_1010_ddsd[1236];

    auto g_y_0_z_0_yz_xz_0_xy = buffer_1010_ddsd[1237];

    auto g_y_0_z_0_yz_xz_0_xz = buffer_1010_ddsd[1238];

    auto g_y_0_z_0_yz_xz_0_yy = buffer_1010_ddsd[1239];

    auto g_y_0_z_0_yz_xz_0_yz = buffer_1010_ddsd[1240];

    auto g_y_0_z_0_yz_xz_0_zz = buffer_1010_ddsd[1241];

    auto g_y_0_z_0_yz_yy_0_xx = buffer_1010_ddsd[1242];

    auto g_y_0_z_0_yz_yy_0_xy = buffer_1010_ddsd[1243];

    auto g_y_0_z_0_yz_yy_0_xz = buffer_1010_ddsd[1244];

    auto g_y_0_z_0_yz_yy_0_yy = buffer_1010_ddsd[1245];

    auto g_y_0_z_0_yz_yy_0_yz = buffer_1010_ddsd[1246];

    auto g_y_0_z_0_yz_yy_0_zz = buffer_1010_ddsd[1247];

    auto g_y_0_z_0_yz_yz_0_xx = buffer_1010_ddsd[1248];

    auto g_y_0_z_0_yz_yz_0_xy = buffer_1010_ddsd[1249];

    auto g_y_0_z_0_yz_yz_0_xz = buffer_1010_ddsd[1250];

    auto g_y_0_z_0_yz_yz_0_yy = buffer_1010_ddsd[1251];

    auto g_y_0_z_0_yz_yz_0_yz = buffer_1010_ddsd[1252];

    auto g_y_0_z_0_yz_yz_0_zz = buffer_1010_ddsd[1253];

    auto g_y_0_z_0_yz_zz_0_xx = buffer_1010_ddsd[1254];

    auto g_y_0_z_0_yz_zz_0_xy = buffer_1010_ddsd[1255];

    auto g_y_0_z_0_yz_zz_0_xz = buffer_1010_ddsd[1256];

    auto g_y_0_z_0_yz_zz_0_yy = buffer_1010_ddsd[1257];

    auto g_y_0_z_0_yz_zz_0_yz = buffer_1010_ddsd[1258];

    auto g_y_0_z_0_yz_zz_0_zz = buffer_1010_ddsd[1259];

    auto g_y_0_z_0_zz_xx_0_xx = buffer_1010_ddsd[1260];

    auto g_y_0_z_0_zz_xx_0_xy = buffer_1010_ddsd[1261];

    auto g_y_0_z_0_zz_xx_0_xz = buffer_1010_ddsd[1262];

    auto g_y_0_z_0_zz_xx_0_yy = buffer_1010_ddsd[1263];

    auto g_y_0_z_0_zz_xx_0_yz = buffer_1010_ddsd[1264];

    auto g_y_0_z_0_zz_xx_0_zz = buffer_1010_ddsd[1265];

    auto g_y_0_z_0_zz_xy_0_xx = buffer_1010_ddsd[1266];

    auto g_y_0_z_0_zz_xy_0_xy = buffer_1010_ddsd[1267];

    auto g_y_0_z_0_zz_xy_0_xz = buffer_1010_ddsd[1268];

    auto g_y_0_z_0_zz_xy_0_yy = buffer_1010_ddsd[1269];

    auto g_y_0_z_0_zz_xy_0_yz = buffer_1010_ddsd[1270];

    auto g_y_0_z_0_zz_xy_0_zz = buffer_1010_ddsd[1271];

    auto g_y_0_z_0_zz_xz_0_xx = buffer_1010_ddsd[1272];

    auto g_y_0_z_0_zz_xz_0_xy = buffer_1010_ddsd[1273];

    auto g_y_0_z_0_zz_xz_0_xz = buffer_1010_ddsd[1274];

    auto g_y_0_z_0_zz_xz_0_yy = buffer_1010_ddsd[1275];

    auto g_y_0_z_0_zz_xz_0_yz = buffer_1010_ddsd[1276];

    auto g_y_0_z_0_zz_xz_0_zz = buffer_1010_ddsd[1277];

    auto g_y_0_z_0_zz_yy_0_xx = buffer_1010_ddsd[1278];

    auto g_y_0_z_0_zz_yy_0_xy = buffer_1010_ddsd[1279];

    auto g_y_0_z_0_zz_yy_0_xz = buffer_1010_ddsd[1280];

    auto g_y_0_z_0_zz_yy_0_yy = buffer_1010_ddsd[1281];

    auto g_y_0_z_0_zz_yy_0_yz = buffer_1010_ddsd[1282];

    auto g_y_0_z_0_zz_yy_0_zz = buffer_1010_ddsd[1283];

    auto g_y_0_z_0_zz_yz_0_xx = buffer_1010_ddsd[1284];

    auto g_y_0_z_0_zz_yz_0_xy = buffer_1010_ddsd[1285];

    auto g_y_0_z_0_zz_yz_0_xz = buffer_1010_ddsd[1286];

    auto g_y_0_z_0_zz_yz_0_yy = buffer_1010_ddsd[1287];

    auto g_y_0_z_0_zz_yz_0_yz = buffer_1010_ddsd[1288];

    auto g_y_0_z_0_zz_yz_0_zz = buffer_1010_ddsd[1289];

    auto g_y_0_z_0_zz_zz_0_xx = buffer_1010_ddsd[1290];

    auto g_y_0_z_0_zz_zz_0_xy = buffer_1010_ddsd[1291];

    auto g_y_0_z_0_zz_zz_0_xz = buffer_1010_ddsd[1292];

    auto g_y_0_z_0_zz_zz_0_yy = buffer_1010_ddsd[1293];

    auto g_y_0_z_0_zz_zz_0_yz = buffer_1010_ddsd[1294];

    auto g_y_0_z_0_zz_zz_0_zz = buffer_1010_ddsd[1295];

    auto g_z_0_x_0_xx_xx_0_xx = buffer_1010_ddsd[1296];

    auto g_z_0_x_0_xx_xx_0_xy = buffer_1010_ddsd[1297];

    auto g_z_0_x_0_xx_xx_0_xz = buffer_1010_ddsd[1298];

    auto g_z_0_x_0_xx_xx_0_yy = buffer_1010_ddsd[1299];

    auto g_z_0_x_0_xx_xx_0_yz = buffer_1010_ddsd[1300];

    auto g_z_0_x_0_xx_xx_0_zz = buffer_1010_ddsd[1301];

    auto g_z_0_x_0_xx_xy_0_xx = buffer_1010_ddsd[1302];

    auto g_z_0_x_0_xx_xy_0_xy = buffer_1010_ddsd[1303];

    auto g_z_0_x_0_xx_xy_0_xz = buffer_1010_ddsd[1304];

    auto g_z_0_x_0_xx_xy_0_yy = buffer_1010_ddsd[1305];

    auto g_z_0_x_0_xx_xy_0_yz = buffer_1010_ddsd[1306];

    auto g_z_0_x_0_xx_xy_0_zz = buffer_1010_ddsd[1307];

    auto g_z_0_x_0_xx_xz_0_xx = buffer_1010_ddsd[1308];

    auto g_z_0_x_0_xx_xz_0_xy = buffer_1010_ddsd[1309];

    auto g_z_0_x_0_xx_xz_0_xz = buffer_1010_ddsd[1310];

    auto g_z_0_x_0_xx_xz_0_yy = buffer_1010_ddsd[1311];

    auto g_z_0_x_0_xx_xz_0_yz = buffer_1010_ddsd[1312];

    auto g_z_0_x_0_xx_xz_0_zz = buffer_1010_ddsd[1313];

    auto g_z_0_x_0_xx_yy_0_xx = buffer_1010_ddsd[1314];

    auto g_z_0_x_0_xx_yy_0_xy = buffer_1010_ddsd[1315];

    auto g_z_0_x_0_xx_yy_0_xz = buffer_1010_ddsd[1316];

    auto g_z_0_x_0_xx_yy_0_yy = buffer_1010_ddsd[1317];

    auto g_z_0_x_0_xx_yy_0_yz = buffer_1010_ddsd[1318];

    auto g_z_0_x_0_xx_yy_0_zz = buffer_1010_ddsd[1319];

    auto g_z_0_x_0_xx_yz_0_xx = buffer_1010_ddsd[1320];

    auto g_z_0_x_0_xx_yz_0_xy = buffer_1010_ddsd[1321];

    auto g_z_0_x_0_xx_yz_0_xz = buffer_1010_ddsd[1322];

    auto g_z_0_x_0_xx_yz_0_yy = buffer_1010_ddsd[1323];

    auto g_z_0_x_0_xx_yz_0_yz = buffer_1010_ddsd[1324];

    auto g_z_0_x_0_xx_yz_0_zz = buffer_1010_ddsd[1325];

    auto g_z_0_x_0_xx_zz_0_xx = buffer_1010_ddsd[1326];

    auto g_z_0_x_0_xx_zz_0_xy = buffer_1010_ddsd[1327];

    auto g_z_0_x_0_xx_zz_0_xz = buffer_1010_ddsd[1328];

    auto g_z_0_x_0_xx_zz_0_yy = buffer_1010_ddsd[1329];

    auto g_z_0_x_0_xx_zz_0_yz = buffer_1010_ddsd[1330];

    auto g_z_0_x_0_xx_zz_0_zz = buffer_1010_ddsd[1331];

    auto g_z_0_x_0_xy_xx_0_xx = buffer_1010_ddsd[1332];

    auto g_z_0_x_0_xy_xx_0_xy = buffer_1010_ddsd[1333];

    auto g_z_0_x_0_xy_xx_0_xz = buffer_1010_ddsd[1334];

    auto g_z_0_x_0_xy_xx_0_yy = buffer_1010_ddsd[1335];

    auto g_z_0_x_0_xy_xx_0_yz = buffer_1010_ddsd[1336];

    auto g_z_0_x_0_xy_xx_0_zz = buffer_1010_ddsd[1337];

    auto g_z_0_x_0_xy_xy_0_xx = buffer_1010_ddsd[1338];

    auto g_z_0_x_0_xy_xy_0_xy = buffer_1010_ddsd[1339];

    auto g_z_0_x_0_xy_xy_0_xz = buffer_1010_ddsd[1340];

    auto g_z_0_x_0_xy_xy_0_yy = buffer_1010_ddsd[1341];

    auto g_z_0_x_0_xy_xy_0_yz = buffer_1010_ddsd[1342];

    auto g_z_0_x_0_xy_xy_0_zz = buffer_1010_ddsd[1343];

    auto g_z_0_x_0_xy_xz_0_xx = buffer_1010_ddsd[1344];

    auto g_z_0_x_0_xy_xz_0_xy = buffer_1010_ddsd[1345];

    auto g_z_0_x_0_xy_xz_0_xz = buffer_1010_ddsd[1346];

    auto g_z_0_x_0_xy_xz_0_yy = buffer_1010_ddsd[1347];

    auto g_z_0_x_0_xy_xz_0_yz = buffer_1010_ddsd[1348];

    auto g_z_0_x_0_xy_xz_0_zz = buffer_1010_ddsd[1349];

    auto g_z_0_x_0_xy_yy_0_xx = buffer_1010_ddsd[1350];

    auto g_z_0_x_0_xy_yy_0_xy = buffer_1010_ddsd[1351];

    auto g_z_0_x_0_xy_yy_0_xz = buffer_1010_ddsd[1352];

    auto g_z_0_x_0_xy_yy_0_yy = buffer_1010_ddsd[1353];

    auto g_z_0_x_0_xy_yy_0_yz = buffer_1010_ddsd[1354];

    auto g_z_0_x_0_xy_yy_0_zz = buffer_1010_ddsd[1355];

    auto g_z_0_x_0_xy_yz_0_xx = buffer_1010_ddsd[1356];

    auto g_z_0_x_0_xy_yz_0_xy = buffer_1010_ddsd[1357];

    auto g_z_0_x_0_xy_yz_0_xz = buffer_1010_ddsd[1358];

    auto g_z_0_x_0_xy_yz_0_yy = buffer_1010_ddsd[1359];

    auto g_z_0_x_0_xy_yz_0_yz = buffer_1010_ddsd[1360];

    auto g_z_0_x_0_xy_yz_0_zz = buffer_1010_ddsd[1361];

    auto g_z_0_x_0_xy_zz_0_xx = buffer_1010_ddsd[1362];

    auto g_z_0_x_0_xy_zz_0_xy = buffer_1010_ddsd[1363];

    auto g_z_0_x_0_xy_zz_0_xz = buffer_1010_ddsd[1364];

    auto g_z_0_x_0_xy_zz_0_yy = buffer_1010_ddsd[1365];

    auto g_z_0_x_0_xy_zz_0_yz = buffer_1010_ddsd[1366];

    auto g_z_0_x_0_xy_zz_0_zz = buffer_1010_ddsd[1367];

    auto g_z_0_x_0_xz_xx_0_xx = buffer_1010_ddsd[1368];

    auto g_z_0_x_0_xz_xx_0_xy = buffer_1010_ddsd[1369];

    auto g_z_0_x_0_xz_xx_0_xz = buffer_1010_ddsd[1370];

    auto g_z_0_x_0_xz_xx_0_yy = buffer_1010_ddsd[1371];

    auto g_z_0_x_0_xz_xx_0_yz = buffer_1010_ddsd[1372];

    auto g_z_0_x_0_xz_xx_0_zz = buffer_1010_ddsd[1373];

    auto g_z_0_x_0_xz_xy_0_xx = buffer_1010_ddsd[1374];

    auto g_z_0_x_0_xz_xy_0_xy = buffer_1010_ddsd[1375];

    auto g_z_0_x_0_xz_xy_0_xz = buffer_1010_ddsd[1376];

    auto g_z_0_x_0_xz_xy_0_yy = buffer_1010_ddsd[1377];

    auto g_z_0_x_0_xz_xy_0_yz = buffer_1010_ddsd[1378];

    auto g_z_0_x_0_xz_xy_0_zz = buffer_1010_ddsd[1379];

    auto g_z_0_x_0_xz_xz_0_xx = buffer_1010_ddsd[1380];

    auto g_z_0_x_0_xz_xz_0_xy = buffer_1010_ddsd[1381];

    auto g_z_0_x_0_xz_xz_0_xz = buffer_1010_ddsd[1382];

    auto g_z_0_x_0_xz_xz_0_yy = buffer_1010_ddsd[1383];

    auto g_z_0_x_0_xz_xz_0_yz = buffer_1010_ddsd[1384];

    auto g_z_0_x_0_xz_xz_0_zz = buffer_1010_ddsd[1385];

    auto g_z_0_x_0_xz_yy_0_xx = buffer_1010_ddsd[1386];

    auto g_z_0_x_0_xz_yy_0_xy = buffer_1010_ddsd[1387];

    auto g_z_0_x_0_xz_yy_0_xz = buffer_1010_ddsd[1388];

    auto g_z_0_x_0_xz_yy_0_yy = buffer_1010_ddsd[1389];

    auto g_z_0_x_0_xz_yy_0_yz = buffer_1010_ddsd[1390];

    auto g_z_0_x_0_xz_yy_0_zz = buffer_1010_ddsd[1391];

    auto g_z_0_x_0_xz_yz_0_xx = buffer_1010_ddsd[1392];

    auto g_z_0_x_0_xz_yz_0_xy = buffer_1010_ddsd[1393];

    auto g_z_0_x_0_xz_yz_0_xz = buffer_1010_ddsd[1394];

    auto g_z_0_x_0_xz_yz_0_yy = buffer_1010_ddsd[1395];

    auto g_z_0_x_0_xz_yz_0_yz = buffer_1010_ddsd[1396];

    auto g_z_0_x_0_xz_yz_0_zz = buffer_1010_ddsd[1397];

    auto g_z_0_x_0_xz_zz_0_xx = buffer_1010_ddsd[1398];

    auto g_z_0_x_0_xz_zz_0_xy = buffer_1010_ddsd[1399];

    auto g_z_0_x_0_xz_zz_0_xz = buffer_1010_ddsd[1400];

    auto g_z_0_x_0_xz_zz_0_yy = buffer_1010_ddsd[1401];

    auto g_z_0_x_0_xz_zz_0_yz = buffer_1010_ddsd[1402];

    auto g_z_0_x_0_xz_zz_0_zz = buffer_1010_ddsd[1403];

    auto g_z_0_x_0_yy_xx_0_xx = buffer_1010_ddsd[1404];

    auto g_z_0_x_0_yy_xx_0_xy = buffer_1010_ddsd[1405];

    auto g_z_0_x_0_yy_xx_0_xz = buffer_1010_ddsd[1406];

    auto g_z_0_x_0_yy_xx_0_yy = buffer_1010_ddsd[1407];

    auto g_z_0_x_0_yy_xx_0_yz = buffer_1010_ddsd[1408];

    auto g_z_0_x_0_yy_xx_0_zz = buffer_1010_ddsd[1409];

    auto g_z_0_x_0_yy_xy_0_xx = buffer_1010_ddsd[1410];

    auto g_z_0_x_0_yy_xy_0_xy = buffer_1010_ddsd[1411];

    auto g_z_0_x_0_yy_xy_0_xz = buffer_1010_ddsd[1412];

    auto g_z_0_x_0_yy_xy_0_yy = buffer_1010_ddsd[1413];

    auto g_z_0_x_0_yy_xy_0_yz = buffer_1010_ddsd[1414];

    auto g_z_0_x_0_yy_xy_0_zz = buffer_1010_ddsd[1415];

    auto g_z_0_x_0_yy_xz_0_xx = buffer_1010_ddsd[1416];

    auto g_z_0_x_0_yy_xz_0_xy = buffer_1010_ddsd[1417];

    auto g_z_0_x_0_yy_xz_0_xz = buffer_1010_ddsd[1418];

    auto g_z_0_x_0_yy_xz_0_yy = buffer_1010_ddsd[1419];

    auto g_z_0_x_0_yy_xz_0_yz = buffer_1010_ddsd[1420];

    auto g_z_0_x_0_yy_xz_0_zz = buffer_1010_ddsd[1421];

    auto g_z_0_x_0_yy_yy_0_xx = buffer_1010_ddsd[1422];

    auto g_z_0_x_0_yy_yy_0_xy = buffer_1010_ddsd[1423];

    auto g_z_0_x_0_yy_yy_0_xz = buffer_1010_ddsd[1424];

    auto g_z_0_x_0_yy_yy_0_yy = buffer_1010_ddsd[1425];

    auto g_z_0_x_0_yy_yy_0_yz = buffer_1010_ddsd[1426];

    auto g_z_0_x_0_yy_yy_0_zz = buffer_1010_ddsd[1427];

    auto g_z_0_x_0_yy_yz_0_xx = buffer_1010_ddsd[1428];

    auto g_z_0_x_0_yy_yz_0_xy = buffer_1010_ddsd[1429];

    auto g_z_0_x_0_yy_yz_0_xz = buffer_1010_ddsd[1430];

    auto g_z_0_x_0_yy_yz_0_yy = buffer_1010_ddsd[1431];

    auto g_z_0_x_0_yy_yz_0_yz = buffer_1010_ddsd[1432];

    auto g_z_0_x_0_yy_yz_0_zz = buffer_1010_ddsd[1433];

    auto g_z_0_x_0_yy_zz_0_xx = buffer_1010_ddsd[1434];

    auto g_z_0_x_0_yy_zz_0_xy = buffer_1010_ddsd[1435];

    auto g_z_0_x_0_yy_zz_0_xz = buffer_1010_ddsd[1436];

    auto g_z_0_x_0_yy_zz_0_yy = buffer_1010_ddsd[1437];

    auto g_z_0_x_0_yy_zz_0_yz = buffer_1010_ddsd[1438];

    auto g_z_0_x_0_yy_zz_0_zz = buffer_1010_ddsd[1439];

    auto g_z_0_x_0_yz_xx_0_xx = buffer_1010_ddsd[1440];

    auto g_z_0_x_0_yz_xx_0_xy = buffer_1010_ddsd[1441];

    auto g_z_0_x_0_yz_xx_0_xz = buffer_1010_ddsd[1442];

    auto g_z_0_x_0_yz_xx_0_yy = buffer_1010_ddsd[1443];

    auto g_z_0_x_0_yz_xx_0_yz = buffer_1010_ddsd[1444];

    auto g_z_0_x_0_yz_xx_0_zz = buffer_1010_ddsd[1445];

    auto g_z_0_x_0_yz_xy_0_xx = buffer_1010_ddsd[1446];

    auto g_z_0_x_0_yz_xy_0_xy = buffer_1010_ddsd[1447];

    auto g_z_0_x_0_yz_xy_0_xz = buffer_1010_ddsd[1448];

    auto g_z_0_x_0_yz_xy_0_yy = buffer_1010_ddsd[1449];

    auto g_z_0_x_0_yz_xy_0_yz = buffer_1010_ddsd[1450];

    auto g_z_0_x_0_yz_xy_0_zz = buffer_1010_ddsd[1451];

    auto g_z_0_x_0_yz_xz_0_xx = buffer_1010_ddsd[1452];

    auto g_z_0_x_0_yz_xz_0_xy = buffer_1010_ddsd[1453];

    auto g_z_0_x_0_yz_xz_0_xz = buffer_1010_ddsd[1454];

    auto g_z_0_x_0_yz_xz_0_yy = buffer_1010_ddsd[1455];

    auto g_z_0_x_0_yz_xz_0_yz = buffer_1010_ddsd[1456];

    auto g_z_0_x_0_yz_xz_0_zz = buffer_1010_ddsd[1457];

    auto g_z_0_x_0_yz_yy_0_xx = buffer_1010_ddsd[1458];

    auto g_z_0_x_0_yz_yy_0_xy = buffer_1010_ddsd[1459];

    auto g_z_0_x_0_yz_yy_0_xz = buffer_1010_ddsd[1460];

    auto g_z_0_x_0_yz_yy_0_yy = buffer_1010_ddsd[1461];

    auto g_z_0_x_0_yz_yy_0_yz = buffer_1010_ddsd[1462];

    auto g_z_0_x_0_yz_yy_0_zz = buffer_1010_ddsd[1463];

    auto g_z_0_x_0_yz_yz_0_xx = buffer_1010_ddsd[1464];

    auto g_z_0_x_0_yz_yz_0_xy = buffer_1010_ddsd[1465];

    auto g_z_0_x_0_yz_yz_0_xz = buffer_1010_ddsd[1466];

    auto g_z_0_x_0_yz_yz_0_yy = buffer_1010_ddsd[1467];

    auto g_z_0_x_0_yz_yz_0_yz = buffer_1010_ddsd[1468];

    auto g_z_0_x_0_yz_yz_0_zz = buffer_1010_ddsd[1469];

    auto g_z_0_x_0_yz_zz_0_xx = buffer_1010_ddsd[1470];

    auto g_z_0_x_0_yz_zz_0_xy = buffer_1010_ddsd[1471];

    auto g_z_0_x_0_yz_zz_0_xz = buffer_1010_ddsd[1472];

    auto g_z_0_x_0_yz_zz_0_yy = buffer_1010_ddsd[1473];

    auto g_z_0_x_0_yz_zz_0_yz = buffer_1010_ddsd[1474];

    auto g_z_0_x_0_yz_zz_0_zz = buffer_1010_ddsd[1475];

    auto g_z_0_x_0_zz_xx_0_xx = buffer_1010_ddsd[1476];

    auto g_z_0_x_0_zz_xx_0_xy = buffer_1010_ddsd[1477];

    auto g_z_0_x_0_zz_xx_0_xz = buffer_1010_ddsd[1478];

    auto g_z_0_x_0_zz_xx_0_yy = buffer_1010_ddsd[1479];

    auto g_z_0_x_0_zz_xx_0_yz = buffer_1010_ddsd[1480];

    auto g_z_0_x_0_zz_xx_0_zz = buffer_1010_ddsd[1481];

    auto g_z_0_x_0_zz_xy_0_xx = buffer_1010_ddsd[1482];

    auto g_z_0_x_0_zz_xy_0_xy = buffer_1010_ddsd[1483];

    auto g_z_0_x_0_zz_xy_0_xz = buffer_1010_ddsd[1484];

    auto g_z_0_x_0_zz_xy_0_yy = buffer_1010_ddsd[1485];

    auto g_z_0_x_0_zz_xy_0_yz = buffer_1010_ddsd[1486];

    auto g_z_0_x_0_zz_xy_0_zz = buffer_1010_ddsd[1487];

    auto g_z_0_x_0_zz_xz_0_xx = buffer_1010_ddsd[1488];

    auto g_z_0_x_0_zz_xz_0_xy = buffer_1010_ddsd[1489];

    auto g_z_0_x_0_zz_xz_0_xz = buffer_1010_ddsd[1490];

    auto g_z_0_x_0_zz_xz_0_yy = buffer_1010_ddsd[1491];

    auto g_z_0_x_0_zz_xz_0_yz = buffer_1010_ddsd[1492];

    auto g_z_0_x_0_zz_xz_0_zz = buffer_1010_ddsd[1493];

    auto g_z_0_x_0_zz_yy_0_xx = buffer_1010_ddsd[1494];

    auto g_z_0_x_0_zz_yy_0_xy = buffer_1010_ddsd[1495];

    auto g_z_0_x_0_zz_yy_0_xz = buffer_1010_ddsd[1496];

    auto g_z_0_x_0_zz_yy_0_yy = buffer_1010_ddsd[1497];

    auto g_z_0_x_0_zz_yy_0_yz = buffer_1010_ddsd[1498];

    auto g_z_0_x_0_zz_yy_0_zz = buffer_1010_ddsd[1499];

    auto g_z_0_x_0_zz_yz_0_xx = buffer_1010_ddsd[1500];

    auto g_z_0_x_0_zz_yz_0_xy = buffer_1010_ddsd[1501];

    auto g_z_0_x_0_zz_yz_0_xz = buffer_1010_ddsd[1502];

    auto g_z_0_x_0_zz_yz_0_yy = buffer_1010_ddsd[1503];

    auto g_z_0_x_0_zz_yz_0_yz = buffer_1010_ddsd[1504];

    auto g_z_0_x_0_zz_yz_0_zz = buffer_1010_ddsd[1505];

    auto g_z_0_x_0_zz_zz_0_xx = buffer_1010_ddsd[1506];

    auto g_z_0_x_0_zz_zz_0_xy = buffer_1010_ddsd[1507];

    auto g_z_0_x_0_zz_zz_0_xz = buffer_1010_ddsd[1508];

    auto g_z_0_x_0_zz_zz_0_yy = buffer_1010_ddsd[1509];

    auto g_z_0_x_0_zz_zz_0_yz = buffer_1010_ddsd[1510];

    auto g_z_0_x_0_zz_zz_0_zz = buffer_1010_ddsd[1511];

    auto g_z_0_y_0_xx_xx_0_xx = buffer_1010_ddsd[1512];

    auto g_z_0_y_0_xx_xx_0_xy = buffer_1010_ddsd[1513];

    auto g_z_0_y_0_xx_xx_0_xz = buffer_1010_ddsd[1514];

    auto g_z_0_y_0_xx_xx_0_yy = buffer_1010_ddsd[1515];

    auto g_z_0_y_0_xx_xx_0_yz = buffer_1010_ddsd[1516];

    auto g_z_0_y_0_xx_xx_0_zz = buffer_1010_ddsd[1517];

    auto g_z_0_y_0_xx_xy_0_xx = buffer_1010_ddsd[1518];

    auto g_z_0_y_0_xx_xy_0_xy = buffer_1010_ddsd[1519];

    auto g_z_0_y_0_xx_xy_0_xz = buffer_1010_ddsd[1520];

    auto g_z_0_y_0_xx_xy_0_yy = buffer_1010_ddsd[1521];

    auto g_z_0_y_0_xx_xy_0_yz = buffer_1010_ddsd[1522];

    auto g_z_0_y_0_xx_xy_0_zz = buffer_1010_ddsd[1523];

    auto g_z_0_y_0_xx_xz_0_xx = buffer_1010_ddsd[1524];

    auto g_z_0_y_0_xx_xz_0_xy = buffer_1010_ddsd[1525];

    auto g_z_0_y_0_xx_xz_0_xz = buffer_1010_ddsd[1526];

    auto g_z_0_y_0_xx_xz_0_yy = buffer_1010_ddsd[1527];

    auto g_z_0_y_0_xx_xz_0_yz = buffer_1010_ddsd[1528];

    auto g_z_0_y_0_xx_xz_0_zz = buffer_1010_ddsd[1529];

    auto g_z_0_y_0_xx_yy_0_xx = buffer_1010_ddsd[1530];

    auto g_z_0_y_0_xx_yy_0_xy = buffer_1010_ddsd[1531];

    auto g_z_0_y_0_xx_yy_0_xz = buffer_1010_ddsd[1532];

    auto g_z_0_y_0_xx_yy_0_yy = buffer_1010_ddsd[1533];

    auto g_z_0_y_0_xx_yy_0_yz = buffer_1010_ddsd[1534];

    auto g_z_0_y_0_xx_yy_0_zz = buffer_1010_ddsd[1535];

    auto g_z_0_y_0_xx_yz_0_xx = buffer_1010_ddsd[1536];

    auto g_z_0_y_0_xx_yz_0_xy = buffer_1010_ddsd[1537];

    auto g_z_0_y_0_xx_yz_0_xz = buffer_1010_ddsd[1538];

    auto g_z_0_y_0_xx_yz_0_yy = buffer_1010_ddsd[1539];

    auto g_z_0_y_0_xx_yz_0_yz = buffer_1010_ddsd[1540];

    auto g_z_0_y_0_xx_yz_0_zz = buffer_1010_ddsd[1541];

    auto g_z_0_y_0_xx_zz_0_xx = buffer_1010_ddsd[1542];

    auto g_z_0_y_0_xx_zz_0_xy = buffer_1010_ddsd[1543];

    auto g_z_0_y_0_xx_zz_0_xz = buffer_1010_ddsd[1544];

    auto g_z_0_y_0_xx_zz_0_yy = buffer_1010_ddsd[1545];

    auto g_z_0_y_0_xx_zz_0_yz = buffer_1010_ddsd[1546];

    auto g_z_0_y_0_xx_zz_0_zz = buffer_1010_ddsd[1547];

    auto g_z_0_y_0_xy_xx_0_xx = buffer_1010_ddsd[1548];

    auto g_z_0_y_0_xy_xx_0_xy = buffer_1010_ddsd[1549];

    auto g_z_0_y_0_xy_xx_0_xz = buffer_1010_ddsd[1550];

    auto g_z_0_y_0_xy_xx_0_yy = buffer_1010_ddsd[1551];

    auto g_z_0_y_0_xy_xx_0_yz = buffer_1010_ddsd[1552];

    auto g_z_0_y_0_xy_xx_0_zz = buffer_1010_ddsd[1553];

    auto g_z_0_y_0_xy_xy_0_xx = buffer_1010_ddsd[1554];

    auto g_z_0_y_0_xy_xy_0_xy = buffer_1010_ddsd[1555];

    auto g_z_0_y_0_xy_xy_0_xz = buffer_1010_ddsd[1556];

    auto g_z_0_y_0_xy_xy_0_yy = buffer_1010_ddsd[1557];

    auto g_z_0_y_0_xy_xy_0_yz = buffer_1010_ddsd[1558];

    auto g_z_0_y_0_xy_xy_0_zz = buffer_1010_ddsd[1559];

    auto g_z_0_y_0_xy_xz_0_xx = buffer_1010_ddsd[1560];

    auto g_z_0_y_0_xy_xz_0_xy = buffer_1010_ddsd[1561];

    auto g_z_0_y_0_xy_xz_0_xz = buffer_1010_ddsd[1562];

    auto g_z_0_y_0_xy_xz_0_yy = buffer_1010_ddsd[1563];

    auto g_z_0_y_0_xy_xz_0_yz = buffer_1010_ddsd[1564];

    auto g_z_0_y_0_xy_xz_0_zz = buffer_1010_ddsd[1565];

    auto g_z_0_y_0_xy_yy_0_xx = buffer_1010_ddsd[1566];

    auto g_z_0_y_0_xy_yy_0_xy = buffer_1010_ddsd[1567];

    auto g_z_0_y_0_xy_yy_0_xz = buffer_1010_ddsd[1568];

    auto g_z_0_y_0_xy_yy_0_yy = buffer_1010_ddsd[1569];

    auto g_z_0_y_0_xy_yy_0_yz = buffer_1010_ddsd[1570];

    auto g_z_0_y_0_xy_yy_0_zz = buffer_1010_ddsd[1571];

    auto g_z_0_y_0_xy_yz_0_xx = buffer_1010_ddsd[1572];

    auto g_z_0_y_0_xy_yz_0_xy = buffer_1010_ddsd[1573];

    auto g_z_0_y_0_xy_yz_0_xz = buffer_1010_ddsd[1574];

    auto g_z_0_y_0_xy_yz_0_yy = buffer_1010_ddsd[1575];

    auto g_z_0_y_0_xy_yz_0_yz = buffer_1010_ddsd[1576];

    auto g_z_0_y_0_xy_yz_0_zz = buffer_1010_ddsd[1577];

    auto g_z_0_y_0_xy_zz_0_xx = buffer_1010_ddsd[1578];

    auto g_z_0_y_0_xy_zz_0_xy = buffer_1010_ddsd[1579];

    auto g_z_0_y_0_xy_zz_0_xz = buffer_1010_ddsd[1580];

    auto g_z_0_y_0_xy_zz_0_yy = buffer_1010_ddsd[1581];

    auto g_z_0_y_0_xy_zz_0_yz = buffer_1010_ddsd[1582];

    auto g_z_0_y_0_xy_zz_0_zz = buffer_1010_ddsd[1583];

    auto g_z_0_y_0_xz_xx_0_xx = buffer_1010_ddsd[1584];

    auto g_z_0_y_0_xz_xx_0_xy = buffer_1010_ddsd[1585];

    auto g_z_0_y_0_xz_xx_0_xz = buffer_1010_ddsd[1586];

    auto g_z_0_y_0_xz_xx_0_yy = buffer_1010_ddsd[1587];

    auto g_z_0_y_0_xz_xx_0_yz = buffer_1010_ddsd[1588];

    auto g_z_0_y_0_xz_xx_0_zz = buffer_1010_ddsd[1589];

    auto g_z_0_y_0_xz_xy_0_xx = buffer_1010_ddsd[1590];

    auto g_z_0_y_0_xz_xy_0_xy = buffer_1010_ddsd[1591];

    auto g_z_0_y_0_xz_xy_0_xz = buffer_1010_ddsd[1592];

    auto g_z_0_y_0_xz_xy_0_yy = buffer_1010_ddsd[1593];

    auto g_z_0_y_0_xz_xy_0_yz = buffer_1010_ddsd[1594];

    auto g_z_0_y_0_xz_xy_0_zz = buffer_1010_ddsd[1595];

    auto g_z_0_y_0_xz_xz_0_xx = buffer_1010_ddsd[1596];

    auto g_z_0_y_0_xz_xz_0_xy = buffer_1010_ddsd[1597];

    auto g_z_0_y_0_xz_xz_0_xz = buffer_1010_ddsd[1598];

    auto g_z_0_y_0_xz_xz_0_yy = buffer_1010_ddsd[1599];

    auto g_z_0_y_0_xz_xz_0_yz = buffer_1010_ddsd[1600];

    auto g_z_0_y_0_xz_xz_0_zz = buffer_1010_ddsd[1601];

    auto g_z_0_y_0_xz_yy_0_xx = buffer_1010_ddsd[1602];

    auto g_z_0_y_0_xz_yy_0_xy = buffer_1010_ddsd[1603];

    auto g_z_0_y_0_xz_yy_0_xz = buffer_1010_ddsd[1604];

    auto g_z_0_y_0_xz_yy_0_yy = buffer_1010_ddsd[1605];

    auto g_z_0_y_0_xz_yy_0_yz = buffer_1010_ddsd[1606];

    auto g_z_0_y_0_xz_yy_0_zz = buffer_1010_ddsd[1607];

    auto g_z_0_y_0_xz_yz_0_xx = buffer_1010_ddsd[1608];

    auto g_z_0_y_0_xz_yz_0_xy = buffer_1010_ddsd[1609];

    auto g_z_0_y_0_xz_yz_0_xz = buffer_1010_ddsd[1610];

    auto g_z_0_y_0_xz_yz_0_yy = buffer_1010_ddsd[1611];

    auto g_z_0_y_0_xz_yz_0_yz = buffer_1010_ddsd[1612];

    auto g_z_0_y_0_xz_yz_0_zz = buffer_1010_ddsd[1613];

    auto g_z_0_y_0_xz_zz_0_xx = buffer_1010_ddsd[1614];

    auto g_z_0_y_0_xz_zz_0_xy = buffer_1010_ddsd[1615];

    auto g_z_0_y_0_xz_zz_0_xz = buffer_1010_ddsd[1616];

    auto g_z_0_y_0_xz_zz_0_yy = buffer_1010_ddsd[1617];

    auto g_z_0_y_0_xz_zz_0_yz = buffer_1010_ddsd[1618];

    auto g_z_0_y_0_xz_zz_0_zz = buffer_1010_ddsd[1619];

    auto g_z_0_y_0_yy_xx_0_xx = buffer_1010_ddsd[1620];

    auto g_z_0_y_0_yy_xx_0_xy = buffer_1010_ddsd[1621];

    auto g_z_0_y_0_yy_xx_0_xz = buffer_1010_ddsd[1622];

    auto g_z_0_y_0_yy_xx_0_yy = buffer_1010_ddsd[1623];

    auto g_z_0_y_0_yy_xx_0_yz = buffer_1010_ddsd[1624];

    auto g_z_0_y_0_yy_xx_0_zz = buffer_1010_ddsd[1625];

    auto g_z_0_y_0_yy_xy_0_xx = buffer_1010_ddsd[1626];

    auto g_z_0_y_0_yy_xy_0_xy = buffer_1010_ddsd[1627];

    auto g_z_0_y_0_yy_xy_0_xz = buffer_1010_ddsd[1628];

    auto g_z_0_y_0_yy_xy_0_yy = buffer_1010_ddsd[1629];

    auto g_z_0_y_0_yy_xy_0_yz = buffer_1010_ddsd[1630];

    auto g_z_0_y_0_yy_xy_0_zz = buffer_1010_ddsd[1631];

    auto g_z_0_y_0_yy_xz_0_xx = buffer_1010_ddsd[1632];

    auto g_z_0_y_0_yy_xz_0_xy = buffer_1010_ddsd[1633];

    auto g_z_0_y_0_yy_xz_0_xz = buffer_1010_ddsd[1634];

    auto g_z_0_y_0_yy_xz_0_yy = buffer_1010_ddsd[1635];

    auto g_z_0_y_0_yy_xz_0_yz = buffer_1010_ddsd[1636];

    auto g_z_0_y_0_yy_xz_0_zz = buffer_1010_ddsd[1637];

    auto g_z_0_y_0_yy_yy_0_xx = buffer_1010_ddsd[1638];

    auto g_z_0_y_0_yy_yy_0_xy = buffer_1010_ddsd[1639];

    auto g_z_0_y_0_yy_yy_0_xz = buffer_1010_ddsd[1640];

    auto g_z_0_y_0_yy_yy_0_yy = buffer_1010_ddsd[1641];

    auto g_z_0_y_0_yy_yy_0_yz = buffer_1010_ddsd[1642];

    auto g_z_0_y_0_yy_yy_0_zz = buffer_1010_ddsd[1643];

    auto g_z_0_y_0_yy_yz_0_xx = buffer_1010_ddsd[1644];

    auto g_z_0_y_0_yy_yz_0_xy = buffer_1010_ddsd[1645];

    auto g_z_0_y_0_yy_yz_0_xz = buffer_1010_ddsd[1646];

    auto g_z_0_y_0_yy_yz_0_yy = buffer_1010_ddsd[1647];

    auto g_z_0_y_0_yy_yz_0_yz = buffer_1010_ddsd[1648];

    auto g_z_0_y_0_yy_yz_0_zz = buffer_1010_ddsd[1649];

    auto g_z_0_y_0_yy_zz_0_xx = buffer_1010_ddsd[1650];

    auto g_z_0_y_0_yy_zz_0_xy = buffer_1010_ddsd[1651];

    auto g_z_0_y_0_yy_zz_0_xz = buffer_1010_ddsd[1652];

    auto g_z_0_y_0_yy_zz_0_yy = buffer_1010_ddsd[1653];

    auto g_z_0_y_0_yy_zz_0_yz = buffer_1010_ddsd[1654];

    auto g_z_0_y_0_yy_zz_0_zz = buffer_1010_ddsd[1655];

    auto g_z_0_y_0_yz_xx_0_xx = buffer_1010_ddsd[1656];

    auto g_z_0_y_0_yz_xx_0_xy = buffer_1010_ddsd[1657];

    auto g_z_0_y_0_yz_xx_0_xz = buffer_1010_ddsd[1658];

    auto g_z_0_y_0_yz_xx_0_yy = buffer_1010_ddsd[1659];

    auto g_z_0_y_0_yz_xx_0_yz = buffer_1010_ddsd[1660];

    auto g_z_0_y_0_yz_xx_0_zz = buffer_1010_ddsd[1661];

    auto g_z_0_y_0_yz_xy_0_xx = buffer_1010_ddsd[1662];

    auto g_z_0_y_0_yz_xy_0_xy = buffer_1010_ddsd[1663];

    auto g_z_0_y_0_yz_xy_0_xz = buffer_1010_ddsd[1664];

    auto g_z_0_y_0_yz_xy_0_yy = buffer_1010_ddsd[1665];

    auto g_z_0_y_0_yz_xy_0_yz = buffer_1010_ddsd[1666];

    auto g_z_0_y_0_yz_xy_0_zz = buffer_1010_ddsd[1667];

    auto g_z_0_y_0_yz_xz_0_xx = buffer_1010_ddsd[1668];

    auto g_z_0_y_0_yz_xz_0_xy = buffer_1010_ddsd[1669];

    auto g_z_0_y_0_yz_xz_0_xz = buffer_1010_ddsd[1670];

    auto g_z_0_y_0_yz_xz_0_yy = buffer_1010_ddsd[1671];

    auto g_z_0_y_0_yz_xz_0_yz = buffer_1010_ddsd[1672];

    auto g_z_0_y_0_yz_xz_0_zz = buffer_1010_ddsd[1673];

    auto g_z_0_y_0_yz_yy_0_xx = buffer_1010_ddsd[1674];

    auto g_z_0_y_0_yz_yy_0_xy = buffer_1010_ddsd[1675];

    auto g_z_0_y_0_yz_yy_0_xz = buffer_1010_ddsd[1676];

    auto g_z_0_y_0_yz_yy_0_yy = buffer_1010_ddsd[1677];

    auto g_z_0_y_0_yz_yy_0_yz = buffer_1010_ddsd[1678];

    auto g_z_0_y_0_yz_yy_0_zz = buffer_1010_ddsd[1679];

    auto g_z_0_y_0_yz_yz_0_xx = buffer_1010_ddsd[1680];

    auto g_z_0_y_0_yz_yz_0_xy = buffer_1010_ddsd[1681];

    auto g_z_0_y_0_yz_yz_0_xz = buffer_1010_ddsd[1682];

    auto g_z_0_y_0_yz_yz_0_yy = buffer_1010_ddsd[1683];

    auto g_z_0_y_0_yz_yz_0_yz = buffer_1010_ddsd[1684];

    auto g_z_0_y_0_yz_yz_0_zz = buffer_1010_ddsd[1685];

    auto g_z_0_y_0_yz_zz_0_xx = buffer_1010_ddsd[1686];

    auto g_z_0_y_0_yz_zz_0_xy = buffer_1010_ddsd[1687];

    auto g_z_0_y_0_yz_zz_0_xz = buffer_1010_ddsd[1688];

    auto g_z_0_y_0_yz_zz_0_yy = buffer_1010_ddsd[1689];

    auto g_z_0_y_0_yz_zz_0_yz = buffer_1010_ddsd[1690];

    auto g_z_0_y_0_yz_zz_0_zz = buffer_1010_ddsd[1691];

    auto g_z_0_y_0_zz_xx_0_xx = buffer_1010_ddsd[1692];

    auto g_z_0_y_0_zz_xx_0_xy = buffer_1010_ddsd[1693];

    auto g_z_0_y_0_zz_xx_0_xz = buffer_1010_ddsd[1694];

    auto g_z_0_y_0_zz_xx_0_yy = buffer_1010_ddsd[1695];

    auto g_z_0_y_0_zz_xx_0_yz = buffer_1010_ddsd[1696];

    auto g_z_0_y_0_zz_xx_0_zz = buffer_1010_ddsd[1697];

    auto g_z_0_y_0_zz_xy_0_xx = buffer_1010_ddsd[1698];

    auto g_z_0_y_0_zz_xy_0_xy = buffer_1010_ddsd[1699];

    auto g_z_0_y_0_zz_xy_0_xz = buffer_1010_ddsd[1700];

    auto g_z_0_y_0_zz_xy_0_yy = buffer_1010_ddsd[1701];

    auto g_z_0_y_0_zz_xy_0_yz = buffer_1010_ddsd[1702];

    auto g_z_0_y_0_zz_xy_0_zz = buffer_1010_ddsd[1703];

    auto g_z_0_y_0_zz_xz_0_xx = buffer_1010_ddsd[1704];

    auto g_z_0_y_0_zz_xz_0_xy = buffer_1010_ddsd[1705];

    auto g_z_0_y_0_zz_xz_0_xz = buffer_1010_ddsd[1706];

    auto g_z_0_y_0_zz_xz_0_yy = buffer_1010_ddsd[1707];

    auto g_z_0_y_0_zz_xz_0_yz = buffer_1010_ddsd[1708];

    auto g_z_0_y_0_zz_xz_0_zz = buffer_1010_ddsd[1709];

    auto g_z_0_y_0_zz_yy_0_xx = buffer_1010_ddsd[1710];

    auto g_z_0_y_0_zz_yy_0_xy = buffer_1010_ddsd[1711];

    auto g_z_0_y_0_zz_yy_0_xz = buffer_1010_ddsd[1712];

    auto g_z_0_y_0_zz_yy_0_yy = buffer_1010_ddsd[1713];

    auto g_z_0_y_0_zz_yy_0_yz = buffer_1010_ddsd[1714];

    auto g_z_0_y_0_zz_yy_0_zz = buffer_1010_ddsd[1715];

    auto g_z_0_y_0_zz_yz_0_xx = buffer_1010_ddsd[1716];

    auto g_z_0_y_0_zz_yz_0_xy = buffer_1010_ddsd[1717];

    auto g_z_0_y_0_zz_yz_0_xz = buffer_1010_ddsd[1718];

    auto g_z_0_y_0_zz_yz_0_yy = buffer_1010_ddsd[1719];

    auto g_z_0_y_0_zz_yz_0_yz = buffer_1010_ddsd[1720];

    auto g_z_0_y_0_zz_yz_0_zz = buffer_1010_ddsd[1721];

    auto g_z_0_y_0_zz_zz_0_xx = buffer_1010_ddsd[1722];

    auto g_z_0_y_0_zz_zz_0_xy = buffer_1010_ddsd[1723];

    auto g_z_0_y_0_zz_zz_0_xz = buffer_1010_ddsd[1724];

    auto g_z_0_y_0_zz_zz_0_yy = buffer_1010_ddsd[1725];

    auto g_z_0_y_0_zz_zz_0_yz = buffer_1010_ddsd[1726];

    auto g_z_0_y_0_zz_zz_0_zz = buffer_1010_ddsd[1727];

    auto g_z_0_z_0_xx_xx_0_xx = buffer_1010_ddsd[1728];

    auto g_z_0_z_0_xx_xx_0_xy = buffer_1010_ddsd[1729];

    auto g_z_0_z_0_xx_xx_0_xz = buffer_1010_ddsd[1730];

    auto g_z_0_z_0_xx_xx_0_yy = buffer_1010_ddsd[1731];

    auto g_z_0_z_0_xx_xx_0_yz = buffer_1010_ddsd[1732];

    auto g_z_0_z_0_xx_xx_0_zz = buffer_1010_ddsd[1733];

    auto g_z_0_z_0_xx_xy_0_xx = buffer_1010_ddsd[1734];

    auto g_z_0_z_0_xx_xy_0_xy = buffer_1010_ddsd[1735];

    auto g_z_0_z_0_xx_xy_0_xz = buffer_1010_ddsd[1736];

    auto g_z_0_z_0_xx_xy_0_yy = buffer_1010_ddsd[1737];

    auto g_z_0_z_0_xx_xy_0_yz = buffer_1010_ddsd[1738];

    auto g_z_0_z_0_xx_xy_0_zz = buffer_1010_ddsd[1739];

    auto g_z_0_z_0_xx_xz_0_xx = buffer_1010_ddsd[1740];

    auto g_z_0_z_0_xx_xz_0_xy = buffer_1010_ddsd[1741];

    auto g_z_0_z_0_xx_xz_0_xz = buffer_1010_ddsd[1742];

    auto g_z_0_z_0_xx_xz_0_yy = buffer_1010_ddsd[1743];

    auto g_z_0_z_0_xx_xz_0_yz = buffer_1010_ddsd[1744];

    auto g_z_0_z_0_xx_xz_0_zz = buffer_1010_ddsd[1745];

    auto g_z_0_z_0_xx_yy_0_xx = buffer_1010_ddsd[1746];

    auto g_z_0_z_0_xx_yy_0_xy = buffer_1010_ddsd[1747];

    auto g_z_0_z_0_xx_yy_0_xz = buffer_1010_ddsd[1748];

    auto g_z_0_z_0_xx_yy_0_yy = buffer_1010_ddsd[1749];

    auto g_z_0_z_0_xx_yy_0_yz = buffer_1010_ddsd[1750];

    auto g_z_0_z_0_xx_yy_0_zz = buffer_1010_ddsd[1751];

    auto g_z_0_z_0_xx_yz_0_xx = buffer_1010_ddsd[1752];

    auto g_z_0_z_0_xx_yz_0_xy = buffer_1010_ddsd[1753];

    auto g_z_0_z_0_xx_yz_0_xz = buffer_1010_ddsd[1754];

    auto g_z_0_z_0_xx_yz_0_yy = buffer_1010_ddsd[1755];

    auto g_z_0_z_0_xx_yz_0_yz = buffer_1010_ddsd[1756];

    auto g_z_0_z_0_xx_yz_0_zz = buffer_1010_ddsd[1757];

    auto g_z_0_z_0_xx_zz_0_xx = buffer_1010_ddsd[1758];

    auto g_z_0_z_0_xx_zz_0_xy = buffer_1010_ddsd[1759];

    auto g_z_0_z_0_xx_zz_0_xz = buffer_1010_ddsd[1760];

    auto g_z_0_z_0_xx_zz_0_yy = buffer_1010_ddsd[1761];

    auto g_z_0_z_0_xx_zz_0_yz = buffer_1010_ddsd[1762];

    auto g_z_0_z_0_xx_zz_0_zz = buffer_1010_ddsd[1763];

    auto g_z_0_z_0_xy_xx_0_xx = buffer_1010_ddsd[1764];

    auto g_z_0_z_0_xy_xx_0_xy = buffer_1010_ddsd[1765];

    auto g_z_0_z_0_xy_xx_0_xz = buffer_1010_ddsd[1766];

    auto g_z_0_z_0_xy_xx_0_yy = buffer_1010_ddsd[1767];

    auto g_z_0_z_0_xy_xx_0_yz = buffer_1010_ddsd[1768];

    auto g_z_0_z_0_xy_xx_0_zz = buffer_1010_ddsd[1769];

    auto g_z_0_z_0_xy_xy_0_xx = buffer_1010_ddsd[1770];

    auto g_z_0_z_0_xy_xy_0_xy = buffer_1010_ddsd[1771];

    auto g_z_0_z_0_xy_xy_0_xz = buffer_1010_ddsd[1772];

    auto g_z_0_z_0_xy_xy_0_yy = buffer_1010_ddsd[1773];

    auto g_z_0_z_0_xy_xy_0_yz = buffer_1010_ddsd[1774];

    auto g_z_0_z_0_xy_xy_0_zz = buffer_1010_ddsd[1775];

    auto g_z_0_z_0_xy_xz_0_xx = buffer_1010_ddsd[1776];

    auto g_z_0_z_0_xy_xz_0_xy = buffer_1010_ddsd[1777];

    auto g_z_0_z_0_xy_xz_0_xz = buffer_1010_ddsd[1778];

    auto g_z_0_z_0_xy_xz_0_yy = buffer_1010_ddsd[1779];

    auto g_z_0_z_0_xy_xz_0_yz = buffer_1010_ddsd[1780];

    auto g_z_0_z_0_xy_xz_0_zz = buffer_1010_ddsd[1781];

    auto g_z_0_z_0_xy_yy_0_xx = buffer_1010_ddsd[1782];

    auto g_z_0_z_0_xy_yy_0_xy = buffer_1010_ddsd[1783];

    auto g_z_0_z_0_xy_yy_0_xz = buffer_1010_ddsd[1784];

    auto g_z_0_z_0_xy_yy_0_yy = buffer_1010_ddsd[1785];

    auto g_z_0_z_0_xy_yy_0_yz = buffer_1010_ddsd[1786];

    auto g_z_0_z_0_xy_yy_0_zz = buffer_1010_ddsd[1787];

    auto g_z_0_z_0_xy_yz_0_xx = buffer_1010_ddsd[1788];

    auto g_z_0_z_0_xy_yz_0_xy = buffer_1010_ddsd[1789];

    auto g_z_0_z_0_xy_yz_0_xz = buffer_1010_ddsd[1790];

    auto g_z_0_z_0_xy_yz_0_yy = buffer_1010_ddsd[1791];

    auto g_z_0_z_0_xy_yz_0_yz = buffer_1010_ddsd[1792];

    auto g_z_0_z_0_xy_yz_0_zz = buffer_1010_ddsd[1793];

    auto g_z_0_z_0_xy_zz_0_xx = buffer_1010_ddsd[1794];

    auto g_z_0_z_0_xy_zz_0_xy = buffer_1010_ddsd[1795];

    auto g_z_0_z_0_xy_zz_0_xz = buffer_1010_ddsd[1796];

    auto g_z_0_z_0_xy_zz_0_yy = buffer_1010_ddsd[1797];

    auto g_z_0_z_0_xy_zz_0_yz = buffer_1010_ddsd[1798];

    auto g_z_0_z_0_xy_zz_0_zz = buffer_1010_ddsd[1799];

    auto g_z_0_z_0_xz_xx_0_xx = buffer_1010_ddsd[1800];

    auto g_z_0_z_0_xz_xx_0_xy = buffer_1010_ddsd[1801];

    auto g_z_0_z_0_xz_xx_0_xz = buffer_1010_ddsd[1802];

    auto g_z_0_z_0_xz_xx_0_yy = buffer_1010_ddsd[1803];

    auto g_z_0_z_0_xz_xx_0_yz = buffer_1010_ddsd[1804];

    auto g_z_0_z_0_xz_xx_0_zz = buffer_1010_ddsd[1805];

    auto g_z_0_z_0_xz_xy_0_xx = buffer_1010_ddsd[1806];

    auto g_z_0_z_0_xz_xy_0_xy = buffer_1010_ddsd[1807];

    auto g_z_0_z_0_xz_xy_0_xz = buffer_1010_ddsd[1808];

    auto g_z_0_z_0_xz_xy_0_yy = buffer_1010_ddsd[1809];

    auto g_z_0_z_0_xz_xy_0_yz = buffer_1010_ddsd[1810];

    auto g_z_0_z_0_xz_xy_0_zz = buffer_1010_ddsd[1811];

    auto g_z_0_z_0_xz_xz_0_xx = buffer_1010_ddsd[1812];

    auto g_z_0_z_0_xz_xz_0_xy = buffer_1010_ddsd[1813];

    auto g_z_0_z_0_xz_xz_0_xz = buffer_1010_ddsd[1814];

    auto g_z_0_z_0_xz_xz_0_yy = buffer_1010_ddsd[1815];

    auto g_z_0_z_0_xz_xz_0_yz = buffer_1010_ddsd[1816];

    auto g_z_0_z_0_xz_xz_0_zz = buffer_1010_ddsd[1817];

    auto g_z_0_z_0_xz_yy_0_xx = buffer_1010_ddsd[1818];

    auto g_z_0_z_0_xz_yy_0_xy = buffer_1010_ddsd[1819];

    auto g_z_0_z_0_xz_yy_0_xz = buffer_1010_ddsd[1820];

    auto g_z_0_z_0_xz_yy_0_yy = buffer_1010_ddsd[1821];

    auto g_z_0_z_0_xz_yy_0_yz = buffer_1010_ddsd[1822];

    auto g_z_0_z_0_xz_yy_0_zz = buffer_1010_ddsd[1823];

    auto g_z_0_z_0_xz_yz_0_xx = buffer_1010_ddsd[1824];

    auto g_z_0_z_0_xz_yz_0_xy = buffer_1010_ddsd[1825];

    auto g_z_0_z_0_xz_yz_0_xz = buffer_1010_ddsd[1826];

    auto g_z_0_z_0_xz_yz_0_yy = buffer_1010_ddsd[1827];

    auto g_z_0_z_0_xz_yz_0_yz = buffer_1010_ddsd[1828];

    auto g_z_0_z_0_xz_yz_0_zz = buffer_1010_ddsd[1829];

    auto g_z_0_z_0_xz_zz_0_xx = buffer_1010_ddsd[1830];

    auto g_z_0_z_0_xz_zz_0_xy = buffer_1010_ddsd[1831];

    auto g_z_0_z_0_xz_zz_0_xz = buffer_1010_ddsd[1832];

    auto g_z_0_z_0_xz_zz_0_yy = buffer_1010_ddsd[1833];

    auto g_z_0_z_0_xz_zz_0_yz = buffer_1010_ddsd[1834];

    auto g_z_0_z_0_xz_zz_0_zz = buffer_1010_ddsd[1835];

    auto g_z_0_z_0_yy_xx_0_xx = buffer_1010_ddsd[1836];

    auto g_z_0_z_0_yy_xx_0_xy = buffer_1010_ddsd[1837];

    auto g_z_0_z_0_yy_xx_0_xz = buffer_1010_ddsd[1838];

    auto g_z_0_z_0_yy_xx_0_yy = buffer_1010_ddsd[1839];

    auto g_z_0_z_0_yy_xx_0_yz = buffer_1010_ddsd[1840];

    auto g_z_0_z_0_yy_xx_0_zz = buffer_1010_ddsd[1841];

    auto g_z_0_z_0_yy_xy_0_xx = buffer_1010_ddsd[1842];

    auto g_z_0_z_0_yy_xy_0_xy = buffer_1010_ddsd[1843];

    auto g_z_0_z_0_yy_xy_0_xz = buffer_1010_ddsd[1844];

    auto g_z_0_z_0_yy_xy_0_yy = buffer_1010_ddsd[1845];

    auto g_z_0_z_0_yy_xy_0_yz = buffer_1010_ddsd[1846];

    auto g_z_0_z_0_yy_xy_0_zz = buffer_1010_ddsd[1847];

    auto g_z_0_z_0_yy_xz_0_xx = buffer_1010_ddsd[1848];

    auto g_z_0_z_0_yy_xz_0_xy = buffer_1010_ddsd[1849];

    auto g_z_0_z_0_yy_xz_0_xz = buffer_1010_ddsd[1850];

    auto g_z_0_z_0_yy_xz_0_yy = buffer_1010_ddsd[1851];

    auto g_z_0_z_0_yy_xz_0_yz = buffer_1010_ddsd[1852];

    auto g_z_0_z_0_yy_xz_0_zz = buffer_1010_ddsd[1853];

    auto g_z_0_z_0_yy_yy_0_xx = buffer_1010_ddsd[1854];

    auto g_z_0_z_0_yy_yy_0_xy = buffer_1010_ddsd[1855];

    auto g_z_0_z_0_yy_yy_0_xz = buffer_1010_ddsd[1856];

    auto g_z_0_z_0_yy_yy_0_yy = buffer_1010_ddsd[1857];

    auto g_z_0_z_0_yy_yy_0_yz = buffer_1010_ddsd[1858];

    auto g_z_0_z_0_yy_yy_0_zz = buffer_1010_ddsd[1859];

    auto g_z_0_z_0_yy_yz_0_xx = buffer_1010_ddsd[1860];

    auto g_z_0_z_0_yy_yz_0_xy = buffer_1010_ddsd[1861];

    auto g_z_0_z_0_yy_yz_0_xz = buffer_1010_ddsd[1862];

    auto g_z_0_z_0_yy_yz_0_yy = buffer_1010_ddsd[1863];

    auto g_z_0_z_0_yy_yz_0_yz = buffer_1010_ddsd[1864];

    auto g_z_0_z_0_yy_yz_0_zz = buffer_1010_ddsd[1865];

    auto g_z_0_z_0_yy_zz_0_xx = buffer_1010_ddsd[1866];

    auto g_z_0_z_0_yy_zz_0_xy = buffer_1010_ddsd[1867];

    auto g_z_0_z_0_yy_zz_0_xz = buffer_1010_ddsd[1868];

    auto g_z_0_z_0_yy_zz_0_yy = buffer_1010_ddsd[1869];

    auto g_z_0_z_0_yy_zz_0_yz = buffer_1010_ddsd[1870];

    auto g_z_0_z_0_yy_zz_0_zz = buffer_1010_ddsd[1871];

    auto g_z_0_z_0_yz_xx_0_xx = buffer_1010_ddsd[1872];

    auto g_z_0_z_0_yz_xx_0_xy = buffer_1010_ddsd[1873];

    auto g_z_0_z_0_yz_xx_0_xz = buffer_1010_ddsd[1874];

    auto g_z_0_z_0_yz_xx_0_yy = buffer_1010_ddsd[1875];

    auto g_z_0_z_0_yz_xx_0_yz = buffer_1010_ddsd[1876];

    auto g_z_0_z_0_yz_xx_0_zz = buffer_1010_ddsd[1877];

    auto g_z_0_z_0_yz_xy_0_xx = buffer_1010_ddsd[1878];

    auto g_z_0_z_0_yz_xy_0_xy = buffer_1010_ddsd[1879];

    auto g_z_0_z_0_yz_xy_0_xz = buffer_1010_ddsd[1880];

    auto g_z_0_z_0_yz_xy_0_yy = buffer_1010_ddsd[1881];

    auto g_z_0_z_0_yz_xy_0_yz = buffer_1010_ddsd[1882];

    auto g_z_0_z_0_yz_xy_0_zz = buffer_1010_ddsd[1883];

    auto g_z_0_z_0_yz_xz_0_xx = buffer_1010_ddsd[1884];

    auto g_z_0_z_0_yz_xz_0_xy = buffer_1010_ddsd[1885];

    auto g_z_0_z_0_yz_xz_0_xz = buffer_1010_ddsd[1886];

    auto g_z_0_z_0_yz_xz_0_yy = buffer_1010_ddsd[1887];

    auto g_z_0_z_0_yz_xz_0_yz = buffer_1010_ddsd[1888];

    auto g_z_0_z_0_yz_xz_0_zz = buffer_1010_ddsd[1889];

    auto g_z_0_z_0_yz_yy_0_xx = buffer_1010_ddsd[1890];

    auto g_z_0_z_0_yz_yy_0_xy = buffer_1010_ddsd[1891];

    auto g_z_0_z_0_yz_yy_0_xz = buffer_1010_ddsd[1892];

    auto g_z_0_z_0_yz_yy_0_yy = buffer_1010_ddsd[1893];

    auto g_z_0_z_0_yz_yy_0_yz = buffer_1010_ddsd[1894];

    auto g_z_0_z_0_yz_yy_0_zz = buffer_1010_ddsd[1895];

    auto g_z_0_z_0_yz_yz_0_xx = buffer_1010_ddsd[1896];

    auto g_z_0_z_0_yz_yz_0_xy = buffer_1010_ddsd[1897];

    auto g_z_0_z_0_yz_yz_0_xz = buffer_1010_ddsd[1898];

    auto g_z_0_z_0_yz_yz_0_yy = buffer_1010_ddsd[1899];

    auto g_z_0_z_0_yz_yz_0_yz = buffer_1010_ddsd[1900];

    auto g_z_0_z_0_yz_yz_0_zz = buffer_1010_ddsd[1901];

    auto g_z_0_z_0_yz_zz_0_xx = buffer_1010_ddsd[1902];

    auto g_z_0_z_0_yz_zz_0_xy = buffer_1010_ddsd[1903];

    auto g_z_0_z_0_yz_zz_0_xz = buffer_1010_ddsd[1904];

    auto g_z_0_z_0_yz_zz_0_yy = buffer_1010_ddsd[1905];

    auto g_z_0_z_0_yz_zz_0_yz = buffer_1010_ddsd[1906];

    auto g_z_0_z_0_yz_zz_0_zz = buffer_1010_ddsd[1907];

    auto g_z_0_z_0_zz_xx_0_xx = buffer_1010_ddsd[1908];

    auto g_z_0_z_0_zz_xx_0_xy = buffer_1010_ddsd[1909];

    auto g_z_0_z_0_zz_xx_0_xz = buffer_1010_ddsd[1910];

    auto g_z_0_z_0_zz_xx_0_yy = buffer_1010_ddsd[1911];

    auto g_z_0_z_0_zz_xx_0_yz = buffer_1010_ddsd[1912];

    auto g_z_0_z_0_zz_xx_0_zz = buffer_1010_ddsd[1913];

    auto g_z_0_z_0_zz_xy_0_xx = buffer_1010_ddsd[1914];

    auto g_z_0_z_0_zz_xy_0_xy = buffer_1010_ddsd[1915];

    auto g_z_0_z_0_zz_xy_0_xz = buffer_1010_ddsd[1916];

    auto g_z_0_z_0_zz_xy_0_yy = buffer_1010_ddsd[1917];

    auto g_z_0_z_0_zz_xy_0_yz = buffer_1010_ddsd[1918];

    auto g_z_0_z_0_zz_xy_0_zz = buffer_1010_ddsd[1919];

    auto g_z_0_z_0_zz_xz_0_xx = buffer_1010_ddsd[1920];

    auto g_z_0_z_0_zz_xz_0_xy = buffer_1010_ddsd[1921];

    auto g_z_0_z_0_zz_xz_0_xz = buffer_1010_ddsd[1922];

    auto g_z_0_z_0_zz_xz_0_yy = buffer_1010_ddsd[1923];

    auto g_z_0_z_0_zz_xz_0_yz = buffer_1010_ddsd[1924];

    auto g_z_0_z_0_zz_xz_0_zz = buffer_1010_ddsd[1925];

    auto g_z_0_z_0_zz_yy_0_xx = buffer_1010_ddsd[1926];

    auto g_z_0_z_0_zz_yy_0_xy = buffer_1010_ddsd[1927];

    auto g_z_0_z_0_zz_yy_0_xz = buffer_1010_ddsd[1928];

    auto g_z_0_z_0_zz_yy_0_yy = buffer_1010_ddsd[1929];

    auto g_z_0_z_0_zz_yy_0_yz = buffer_1010_ddsd[1930];

    auto g_z_0_z_0_zz_yy_0_zz = buffer_1010_ddsd[1931];

    auto g_z_0_z_0_zz_yz_0_xx = buffer_1010_ddsd[1932];

    auto g_z_0_z_0_zz_yz_0_xy = buffer_1010_ddsd[1933];

    auto g_z_0_z_0_zz_yz_0_xz = buffer_1010_ddsd[1934];

    auto g_z_0_z_0_zz_yz_0_yy = buffer_1010_ddsd[1935];

    auto g_z_0_z_0_zz_yz_0_yz = buffer_1010_ddsd[1936];

    auto g_z_0_z_0_zz_yz_0_zz = buffer_1010_ddsd[1937];

    auto g_z_0_z_0_zz_zz_0_xx = buffer_1010_ddsd[1938];

    auto g_z_0_z_0_zz_zz_0_xy = buffer_1010_ddsd[1939];

    auto g_z_0_z_0_zz_zz_0_xz = buffer_1010_ddsd[1940];

    auto g_z_0_z_0_zz_zz_0_yy = buffer_1010_ddsd[1941];

    auto g_z_0_z_0_zz_zz_0_yz = buffer_1010_ddsd[1942];

    auto g_z_0_z_0_zz_zz_0_zz = buffer_1010_ddsd[1943];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_0_x_0_xx_xx_0_xx, g_x_0_x_0_xx_xx_0_xy, g_x_0_x_0_xx_xx_0_xz, g_x_0_x_0_xx_xx_0_yy, g_x_0_x_0_xx_xx_0_yz, g_x_0_x_0_xx_xx_0_zz, g_x_xx_x_xx, g_x_xx_x_xy, g_x_xx_x_xz, g_x_xx_x_yy, g_x_xx_x_yz, g_x_xx_x_zz, g_xxx_xx_x_xx, g_xxx_xx_x_xy, g_xxx_xx_x_xz, g_xxx_xx_x_yy, g_xxx_xx_x_yz, g_xxx_xx_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_xx_xx_0_xx[i] = -4.0 * g_x_xx_x_xx[i] * c_exps[i] + 4.0 * g_xxx_xx_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_xx_0_xy[i] = -4.0 * g_x_xx_x_xy[i] * c_exps[i] + 4.0 * g_xxx_xx_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_xx_0_xz[i] = -4.0 * g_x_xx_x_xz[i] * c_exps[i] + 4.0 * g_xxx_xx_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_xx_0_yy[i] = -4.0 * g_x_xx_x_yy[i] * c_exps[i] + 4.0 * g_xxx_xx_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_xx_0_yz[i] = -4.0 * g_x_xx_x_yz[i] * c_exps[i] + 4.0 * g_xxx_xx_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_xx_0_zz[i] = -4.0 * g_x_xx_x_zz[i] * c_exps[i] + 4.0 * g_xxx_xx_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_0_x_0_xx_xy_0_xx, g_x_0_x_0_xx_xy_0_xy, g_x_0_x_0_xx_xy_0_xz, g_x_0_x_0_xx_xy_0_yy, g_x_0_x_0_xx_xy_0_yz, g_x_0_x_0_xx_xy_0_zz, g_x_xy_x_xx, g_x_xy_x_xy, g_x_xy_x_xz, g_x_xy_x_yy, g_x_xy_x_yz, g_x_xy_x_zz, g_xxx_xy_x_xx, g_xxx_xy_x_xy, g_xxx_xy_x_xz, g_xxx_xy_x_yy, g_xxx_xy_x_yz, g_xxx_xy_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_xx_xy_0_xx[i] = -4.0 * g_x_xy_x_xx[i] * c_exps[i] + 4.0 * g_xxx_xy_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_xy_0_xy[i] = -4.0 * g_x_xy_x_xy[i] * c_exps[i] + 4.0 * g_xxx_xy_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_xy_0_xz[i] = -4.0 * g_x_xy_x_xz[i] * c_exps[i] + 4.0 * g_xxx_xy_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_xy_0_yy[i] = -4.0 * g_x_xy_x_yy[i] * c_exps[i] + 4.0 * g_xxx_xy_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_xy_0_yz[i] = -4.0 * g_x_xy_x_yz[i] * c_exps[i] + 4.0 * g_xxx_xy_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_xy_0_zz[i] = -4.0 * g_x_xy_x_zz[i] * c_exps[i] + 4.0 * g_xxx_xy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_0_x_0_xx_xz_0_xx, g_x_0_x_0_xx_xz_0_xy, g_x_0_x_0_xx_xz_0_xz, g_x_0_x_0_xx_xz_0_yy, g_x_0_x_0_xx_xz_0_yz, g_x_0_x_0_xx_xz_0_zz, g_x_xz_x_xx, g_x_xz_x_xy, g_x_xz_x_xz, g_x_xz_x_yy, g_x_xz_x_yz, g_x_xz_x_zz, g_xxx_xz_x_xx, g_xxx_xz_x_xy, g_xxx_xz_x_xz, g_xxx_xz_x_yy, g_xxx_xz_x_yz, g_xxx_xz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_xx_xz_0_xx[i] = -4.0 * g_x_xz_x_xx[i] * c_exps[i] + 4.0 * g_xxx_xz_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_xz_0_xy[i] = -4.0 * g_x_xz_x_xy[i] * c_exps[i] + 4.0 * g_xxx_xz_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_xz_0_xz[i] = -4.0 * g_x_xz_x_xz[i] * c_exps[i] + 4.0 * g_xxx_xz_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_xz_0_yy[i] = -4.0 * g_x_xz_x_yy[i] * c_exps[i] + 4.0 * g_xxx_xz_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_xz_0_yz[i] = -4.0 * g_x_xz_x_yz[i] * c_exps[i] + 4.0 * g_xxx_xz_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_xz_0_zz[i] = -4.0 * g_x_xz_x_zz[i] * c_exps[i] + 4.0 * g_xxx_xz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_x_0_x_0_xx_yy_0_xx, g_x_0_x_0_xx_yy_0_xy, g_x_0_x_0_xx_yy_0_xz, g_x_0_x_0_xx_yy_0_yy, g_x_0_x_0_xx_yy_0_yz, g_x_0_x_0_xx_yy_0_zz, g_x_yy_x_xx, g_x_yy_x_xy, g_x_yy_x_xz, g_x_yy_x_yy, g_x_yy_x_yz, g_x_yy_x_zz, g_xxx_yy_x_xx, g_xxx_yy_x_xy, g_xxx_yy_x_xz, g_xxx_yy_x_yy, g_xxx_yy_x_yz, g_xxx_yy_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_xx_yy_0_xx[i] = -4.0 * g_x_yy_x_xx[i] * c_exps[i] + 4.0 * g_xxx_yy_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_yy_0_xy[i] = -4.0 * g_x_yy_x_xy[i] * c_exps[i] + 4.0 * g_xxx_yy_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_yy_0_xz[i] = -4.0 * g_x_yy_x_xz[i] * c_exps[i] + 4.0 * g_xxx_yy_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_yy_0_yy[i] = -4.0 * g_x_yy_x_yy[i] * c_exps[i] + 4.0 * g_xxx_yy_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_yy_0_yz[i] = -4.0 * g_x_yy_x_yz[i] * c_exps[i] + 4.0 * g_xxx_yy_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_yy_0_zz[i] = -4.0 * g_x_yy_x_zz[i] * c_exps[i] + 4.0 * g_xxx_yy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_x_0_x_0_xx_yz_0_xx, g_x_0_x_0_xx_yz_0_xy, g_x_0_x_0_xx_yz_0_xz, g_x_0_x_0_xx_yz_0_yy, g_x_0_x_0_xx_yz_0_yz, g_x_0_x_0_xx_yz_0_zz, g_x_yz_x_xx, g_x_yz_x_xy, g_x_yz_x_xz, g_x_yz_x_yy, g_x_yz_x_yz, g_x_yz_x_zz, g_xxx_yz_x_xx, g_xxx_yz_x_xy, g_xxx_yz_x_xz, g_xxx_yz_x_yy, g_xxx_yz_x_yz, g_xxx_yz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_xx_yz_0_xx[i] = -4.0 * g_x_yz_x_xx[i] * c_exps[i] + 4.0 * g_xxx_yz_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_yz_0_xy[i] = -4.0 * g_x_yz_x_xy[i] * c_exps[i] + 4.0 * g_xxx_yz_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_yz_0_xz[i] = -4.0 * g_x_yz_x_xz[i] * c_exps[i] + 4.0 * g_xxx_yz_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_yz_0_yy[i] = -4.0 * g_x_yz_x_yy[i] * c_exps[i] + 4.0 * g_xxx_yz_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_yz_0_yz[i] = -4.0 * g_x_yz_x_yz[i] * c_exps[i] + 4.0 * g_xxx_yz_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_yz_0_zz[i] = -4.0 * g_x_yz_x_zz[i] * c_exps[i] + 4.0 * g_xxx_yz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_x_0_x_0_xx_zz_0_xx, g_x_0_x_0_xx_zz_0_xy, g_x_0_x_0_xx_zz_0_xz, g_x_0_x_0_xx_zz_0_yy, g_x_0_x_0_xx_zz_0_yz, g_x_0_x_0_xx_zz_0_zz, g_x_zz_x_xx, g_x_zz_x_xy, g_x_zz_x_xz, g_x_zz_x_yy, g_x_zz_x_yz, g_x_zz_x_zz, g_xxx_zz_x_xx, g_xxx_zz_x_xy, g_xxx_zz_x_xz, g_xxx_zz_x_yy, g_xxx_zz_x_yz, g_xxx_zz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_xx_zz_0_xx[i] = -4.0 * g_x_zz_x_xx[i] * c_exps[i] + 4.0 * g_xxx_zz_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_zz_0_xy[i] = -4.0 * g_x_zz_x_xy[i] * c_exps[i] + 4.0 * g_xxx_zz_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_zz_0_xz[i] = -4.0 * g_x_zz_x_xz[i] * c_exps[i] + 4.0 * g_xxx_zz_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_zz_0_yy[i] = -4.0 * g_x_zz_x_yy[i] * c_exps[i] + 4.0 * g_xxx_zz_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_zz_0_yz[i] = -4.0 * g_x_zz_x_yz[i] * c_exps[i] + 4.0 * g_xxx_zz_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_zz_0_zz[i] = -4.0 * g_x_zz_x_zz[i] * c_exps[i] + 4.0 * g_xxx_zz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_x_0_x_0_xy_xx_0_xx, g_x_0_x_0_xy_xx_0_xy, g_x_0_x_0_xy_xx_0_xz, g_x_0_x_0_xy_xx_0_yy, g_x_0_x_0_xy_xx_0_yz, g_x_0_x_0_xy_xx_0_zz, g_xxy_xx_x_xx, g_xxy_xx_x_xy, g_xxy_xx_x_xz, g_xxy_xx_x_yy, g_xxy_xx_x_yz, g_xxy_xx_x_zz, g_y_xx_x_xx, g_y_xx_x_xy, g_y_xx_x_xz, g_y_xx_x_yy, g_y_xx_x_yz, g_y_xx_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_xy_xx_0_xx[i] = -2.0 * g_y_xx_x_xx[i] * c_exps[i] + 4.0 * g_xxy_xx_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_xx_0_xy[i] = -2.0 * g_y_xx_x_xy[i] * c_exps[i] + 4.0 * g_xxy_xx_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_xx_0_xz[i] = -2.0 * g_y_xx_x_xz[i] * c_exps[i] + 4.0 * g_xxy_xx_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_xx_0_yy[i] = -2.0 * g_y_xx_x_yy[i] * c_exps[i] + 4.0 * g_xxy_xx_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_xx_0_yz[i] = -2.0 * g_y_xx_x_yz[i] * c_exps[i] + 4.0 * g_xxy_xx_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_xx_0_zz[i] = -2.0 * g_y_xx_x_zz[i] * c_exps[i] + 4.0 * g_xxy_xx_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_x_0_x_0_xy_xy_0_xx, g_x_0_x_0_xy_xy_0_xy, g_x_0_x_0_xy_xy_0_xz, g_x_0_x_0_xy_xy_0_yy, g_x_0_x_0_xy_xy_0_yz, g_x_0_x_0_xy_xy_0_zz, g_xxy_xy_x_xx, g_xxy_xy_x_xy, g_xxy_xy_x_xz, g_xxy_xy_x_yy, g_xxy_xy_x_yz, g_xxy_xy_x_zz, g_y_xy_x_xx, g_y_xy_x_xy, g_y_xy_x_xz, g_y_xy_x_yy, g_y_xy_x_yz, g_y_xy_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_xy_xy_0_xx[i] = -2.0 * g_y_xy_x_xx[i] * c_exps[i] + 4.0 * g_xxy_xy_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_xy_0_xy[i] = -2.0 * g_y_xy_x_xy[i] * c_exps[i] + 4.0 * g_xxy_xy_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_xy_0_xz[i] = -2.0 * g_y_xy_x_xz[i] * c_exps[i] + 4.0 * g_xxy_xy_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_xy_0_yy[i] = -2.0 * g_y_xy_x_yy[i] * c_exps[i] + 4.0 * g_xxy_xy_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_xy_0_yz[i] = -2.0 * g_y_xy_x_yz[i] * c_exps[i] + 4.0 * g_xxy_xy_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_xy_0_zz[i] = -2.0 * g_y_xy_x_zz[i] * c_exps[i] + 4.0 * g_xxy_xy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_x_0_x_0_xy_xz_0_xx, g_x_0_x_0_xy_xz_0_xy, g_x_0_x_0_xy_xz_0_xz, g_x_0_x_0_xy_xz_0_yy, g_x_0_x_0_xy_xz_0_yz, g_x_0_x_0_xy_xz_0_zz, g_xxy_xz_x_xx, g_xxy_xz_x_xy, g_xxy_xz_x_xz, g_xxy_xz_x_yy, g_xxy_xz_x_yz, g_xxy_xz_x_zz, g_y_xz_x_xx, g_y_xz_x_xy, g_y_xz_x_xz, g_y_xz_x_yy, g_y_xz_x_yz, g_y_xz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_xy_xz_0_xx[i] = -2.0 * g_y_xz_x_xx[i] * c_exps[i] + 4.0 * g_xxy_xz_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_xz_0_xy[i] = -2.0 * g_y_xz_x_xy[i] * c_exps[i] + 4.0 * g_xxy_xz_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_xz_0_xz[i] = -2.0 * g_y_xz_x_xz[i] * c_exps[i] + 4.0 * g_xxy_xz_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_xz_0_yy[i] = -2.0 * g_y_xz_x_yy[i] * c_exps[i] + 4.0 * g_xxy_xz_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_xz_0_yz[i] = -2.0 * g_y_xz_x_yz[i] * c_exps[i] + 4.0 * g_xxy_xz_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_xz_0_zz[i] = -2.0 * g_y_xz_x_zz[i] * c_exps[i] + 4.0 * g_xxy_xz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_x_0_x_0_xy_yy_0_xx, g_x_0_x_0_xy_yy_0_xy, g_x_0_x_0_xy_yy_0_xz, g_x_0_x_0_xy_yy_0_yy, g_x_0_x_0_xy_yy_0_yz, g_x_0_x_0_xy_yy_0_zz, g_xxy_yy_x_xx, g_xxy_yy_x_xy, g_xxy_yy_x_xz, g_xxy_yy_x_yy, g_xxy_yy_x_yz, g_xxy_yy_x_zz, g_y_yy_x_xx, g_y_yy_x_xy, g_y_yy_x_xz, g_y_yy_x_yy, g_y_yy_x_yz, g_y_yy_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_xy_yy_0_xx[i] = -2.0 * g_y_yy_x_xx[i] * c_exps[i] + 4.0 * g_xxy_yy_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_yy_0_xy[i] = -2.0 * g_y_yy_x_xy[i] * c_exps[i] + 4.0 * g_xxy_yy_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_yy_0_xz[i] = -2.0 * g_y_yy_x_xz[i] * c_exps[i] + 4.0 * g_xxy_yy_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_yy_0_yy[i] = -2.0 * g_y_yy_x_yy[i] * c_exps[i] + 4.0 * g_xxy_yy_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_yy_0_yz[i] = -2.0 * g_y_yy_x_yz[i] * c_exps[i] + 4.0 * g_xxy_yy_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_yy_0_zz[i] = -2.0 * g_y_yy_x_zz[i] * c_exps[i] + 4.0 * g_xxy_yy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_x_0_x_0_xy_yz_0_xx, g_x_0_x_0_xy_yz_0_xy, g_x_0_x_0_xy_yz_0_xz, g_x_0_x_0_xy_yz_0_yy, g_x_0_x_0_xy_yz_0_yz, g_x_0_x_0_xy_yz_0_zz, g_xxy_yz_x_xx, g_xxy_yz_x_xy, g_xxy_yz_x_xz, g_xxy_yz_x_yy, g_xxy_yz_x_yz, g_xxy_yz_x_zz, g_y_yz_x_xx, g_y_yz_x_xy, g_y_yz_x_xz, g_y_yz_x_yy, g_y_yz_x_yz, g_y_yz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_xy_yz_0_xx[i] = -2.0 * g_y_yz_x_xx[i] * c_exps[i] + 4.0 * g_xxy_yz_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_yz_0_xy[i] = -2.0 * g_y_yz_x_xy[i] * c_exps[i] + 4.0 * g_xxy_yz_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_yz_0_xz[i] = -2.0 * g_y_yz_x_xz[i] * c_exps[i] + 4.0 * g_xxy_yz_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_yz_0_yy[i] = -2.0 * g_y_yz_x_yy[i] * c_exps[i] + 4.0 * g_xxy_yz_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_yz_0_yz[i] = -2.0 * g_y_yz_x_yz[i] * c_exps[i] + 4.0 * g_xxy_yz_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_yz_0_zz[i] = -2.0 * g_y_yz_x_zz[i] * c_exps[i] + 4.0 * g_xxy_yz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_x_0_x_0_xy_zz_0_xx, g_x_0_x_0_xy_zz_0_xy, g_x_0_x_0_xy_zz_0_xz, g_x_0_x_0_xy_zz_0_yy, g_x_0_x_0_xy_zz_0_yz, g_x_0_x_0_xy_zz_0_zz, g_xxy_zz_x_xx, g_xxy_zz_x_xy, g_xxy_zz_x_xz, g_xxy_zz_x_yy, g_xxy_zz_x_yz, g_xxy_zz_x_zz, g_y_zz_x_xx, g_y_zz_x_xy, g_y_zz_x_xz, g_y_zz_x_yy, g_y_zz_x_yz, g_y_zz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_xy_zz_0_xx[i] = -2.0 * g_y_zz_x_xx[i] * c_exps[i] + 4.0 * g_xxy_zz_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_zz_0_xy[i] = -2.0 * g_y_zz_x_xy[i] * c_exps[i] + 4.0 * g_xxy_zz_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_zz_0_xz[i] = -2.0 * g_y_zz_x_xz[i] * c_exps[i] + 4.0 * g_xxy_zz_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_zz_0_yy[i] = -2.0 * g_y_zz_x_yy[i] * c_exps[i] + 4.0 * g_xxy_zz_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_zz_0_yz[i] = -2.0 * g_y_zz_x_yz[i] * c_exps[i] + 4.0 * g_xxy_zz_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_zz_0_zz[i] = -2.0 * g_y_zz_x_zz[i] * c_exps[i] + 4.0 * g_xxy_zz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_x_0_x_0_xz_xx_0_xx, g_x_0_x_0_xz_xx_0_xy, g_x_0_x_0_xz_xx_0_xz, g_x_0_x_0_xz_xx_0_yy, g_x_0_x_0_xz_xx_0_yz, g_x_0_x_0_xz_xx_0_zz, g_xxz_xx_x_xx, g_xxz_xx_x_xy, g_xxz_xx_x_xz, g_xxz_xx_x_yy, g_xxz_xx_x_yz, g_xxz_xx_x_zz, g_z_xx_x_xx, g_z_xx_x_xy, g_z_xx_x_xz, g_z_xx_x_yy, g_z_xx_x_yz, g_z_xx_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_xz_xx_0_xx[i] = -2.0 * g_z_xx_x_xx[i] * c_exps[i] + 4.0 * g_xxz_xx_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_xx_0_xy[i] = -2.0 * g_z_xx_x_xy[i] * c_exps[i] + 4.0 * g_xxz_xx_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_xx_0_xz[i] = -2.0 * g_z_xx_x_xz[i] * c_exps[i] + 4.0 * g_xxz_xx_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_xx_0_yy[i] = -2.0 * g_z_xx_x_yy[i] * c_exps[i] + 4.0 * g_xxz_xx_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_xx_0_yz[i] = -2.0 * g_z_xx_x_yz[i] * c_exps[i] + 4.0 * g_xxz_xx_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_xx_0_zz[i] = -2.0 * g_z_xx_x_zz[i] * c_exps[i] + 4.0 * g_xxz_xx_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_x_0_x_0_xz_xy_0_xx, g_x_0_x_0_xz_xy_0_xy, g_x_0_x_0_xz_xy_0_xz, g_x_0_x_0_xz_xy_0_yy, g_x_0_x_0_xz_xy_0_yz, g_x_0_x_0_xz_xy_0_zz, g_xxz_xy_x_xx, g_xxz_xy_x_xy, g_xxz_xy_x_xz, g_xxz_xy_x_yy, g_xxz_xy_x_yz, g_xxz_xy_x_zz, g_z_xy_x_xx, g_z_xy_x_xy, g_z_xy_x_xz, g_z_xy_x_yy, g_z_xy_x_yz, g_z_xy_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_xz_xy_0_xx[i] = -2.0 * g_z_xy_x_xx[i] * c_exps[i] + 4.0 * g_xxz_xy_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_xy_0_xy[i] = -2.0 * g_z_xy_x_xy[i] * c_exps[i] + 4.0 * g_xxz_xy_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_xy_0_xz[i] = -2.0 * g_z_xy_x_xz[i] * c_exps[i] + 4.0 * g_xxz_xy_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_xy_0_yy[i] = -2.0 * g_z_xy_x_yy[i] * c_exps[i] + 4.0 * g_xxz_xy_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_xy_0_yz[i] = -2.0 * g_z_xy_x_yz[i] * c_exps[i] + 4.0 * g_xxz_xy_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_xy_0_zz[i] = -2.0 * g_z_xy_x_zz[i] * c_exps[i] + 4.0 * g_xxz_xy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_x_0_x_0_xz_xz_0_xx, g_x_0_x_0_xz_xz_0_xy, g_x_0_x_0_xz_xz_0_xz, g_x_0_x_0_xz_xz_0_yy, g_x_0_x_0_xz_xz_0_yz, g_x_0_x_0_xz_xz_0_zz, g_xxz_xz_x_xx, g_xxz_xz_x_xy, g_xxz_xz_x_xz, g_xxz_xz_x_yy, g_xxz_xz_x_yz, g_xxz_xz_x_zz, g_z_xz_x_xx, g_z_xz_x_xy, g_z_xz_x_xz, g_z_xz_x_yy, g_z_xz_x_yz, g_z_xz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_xz_xz_0_xx[i] = -2.0 * g_z_xz_x_xx[i] * c_exps[i] + 4.0 * g_xxz_xz_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_xz_0_xy[i] = -2.0 * g_z_xz_x_xy[i] * c_exps[i] + 4.0 * g_xxz_xz_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_xz_0_xz[i] = -2.0 * g_z_xz_x_xz[i] * c_exps[i] + 4.0 * g_xxz_xz_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_xz_0_yy[i] = -2.0 * g_z_xz_x_yy[i] * c_exps[i] + 4.0 * g_xxz_xz_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_xz_0_yz[i] = -2.0 * g_z_xz_x_yz[i] * c_exps[i] + 4.0 * g_xxz_xz_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_xz_0_zz[i] = -2.0 * g_z_xz_x_zz[i] * c_exps[i] + 4.0 * g_xxz_xz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_x_0_x_0_xz_yy_0_xx, g_x_0_x_0_xz_yy_0_xy, g_x_0_x_0_xz_yy_0_xz, g_x_0_x_0_xz_yy_0_yy, g_x_0_x_0_xz_yy_0_yz, g_x_0_x_0_xz_yy_0_zz, g_xxz_yy_x_xx, g_xxz_yy_x_xy, g_xxz_yy_x_xz, g_xxz_yy_x_yy, g_xxz_yy_x_yz, g_xxz_yy_x_zz, g_z_yy_x_xx, g_z_yy_x_xy, g_z_yy_x_xz, g_z_yy_x_yy, g_z_yy_x_yz, g_z_yy_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_xz_yy_0_xx[i] = -2.0 * g_z_yy_x_xx[i] * c_exps[i] + 4.0 * g_xxz_yy_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_yy_0_xy[i] = -2.0 * g_z_yy_x_xy[i] * c_exps[i] + 4.0 * g_xxz_yy_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_yy_0_xz[i] = -2.0 * g_z_yy_x_xz[i] * c_exps[i] + 4.0 * g_xxz_yy_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_yy_0_yy[i] = -2.0 * g_z_yy_x_yy[i] * c_exps[i] + 4.0 * g_xxz_yy_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_yy_0_yz[i] = -2.0 * g_z_yy_x_yz[i] * c_exps[i] + 4.0 * g_xxz_yy_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_yy_0_zz[i] = -2.0 * g_z_yy_x_zz[i] * c_exps[i] + 4.0 * g_xxz_yy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_x_0_x_0_xz_yz_0_xx, g_x_0_x_0_xz_yz_0_xy, g_x_0_x_0_xz_yz_0_xz, g_x_0_x_0_xz_yz_0_yy, g_x_0_x_0_xz_yz_0_yz, g_x_0_x_0_xz_yz_0_zz, g_xxz_yz_x_xx, g_xxz_yz_x_xy, g_xxz_yz_x_xz, g_xxz_yz_x_yy, g_xxz_yz_x_yz, g_xxz_yz_x_zz, g_z_yz_x_xx, g_z_yz_x_xy, g_z_yz_x_xz, g_z_yz_x_yy, g_z_yz_x_yz, g_z_yz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_xz_yz_0_xx[i] = -2.0 * g_z_yz_x_xx[i] * c_exps[i] + 4.0 * g_xxz_yz_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_yz_0_xy[i] = -2.0 * g_z_yz_x_xy[i] * c_exps[i] + 4.0 * g_xxz_yz_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_yz_0_xz[i] = -2.0 * g_z_yz_x_xz[i] * c_exps[i] + 4.0 * g_xxz_yz_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_yz_0_yy[i] = -2.0 * g_z_yz_x_yy[i] * c_exps[i] + 4.0 * g_xxz_yz_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_yz_0_yz[i] = -2.0 * g_z_yz_x_yz[i] * c_exps[i] + 4.0 * g_xxz_yz_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_yz_0_zz[i] = -2.0 * g_z_yz_x_zz[i] * c_exps[i] + 4.0 * g_xxz_yz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_x_0_x_0_xz_zz_0_xx, g_x_0_x_0_xz_zz_0_xy, g_x_0_x_0_xz_zz_0_xz, g_x_0_x_0_xz_zz_0_yy, g_x_0_x_0_xz_zz_0_yz, g_x_0_x_0_xz_zz_0_zz, g_xxz_zz_x_xx, g_xxz_zz_x_xy, g_xxz_zz_x_xz, g_xxz_zz_x_yy, g_xxz_zz_x_yz, g_xxz_zz_x_zz, g_z_zz_x_xx, g_z_zz_x_xy, g_z_zz_x_xz, g_z_zz_x_yy, g_z_zz_x_yz, g_z_zz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_xz_zz_0_xx[i] = -2.0 * g_z_zz_x_xx[i] * c_exps[i] + 4.0 * g_xxz_zz_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_zz_0_xy[i] = -2.0 * g_z_zz_x_xy[i] * c_exps[i] + 4.0 * g_xxz_zz_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_zz_0_xz[i] = -2.0 * g_z_zz_x_xz[i] * c_exps[i] + 4.0 * g_xxz_zz_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_zz_0_yy[i] = -2.0 * g_z_zz_x_yy[i] * c_exps[i] + 4.0 * g_xxz_zz_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_zz_0_yz[i] = -2.0 * g_z_zz_x_yz[i] * c_exps[i] + 4.0 * g_xxz_zz_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_zz_0_zz[i] = -2.0 * g_z_zz_x_zz[i] * c_exps[i] + 4.0 * g_xxz_zz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_x_0_x_0_yy_xx_0_xx, g_x_0_x_0_yy_xx_0_xy, g_x_0_x_0_yy_xx_0_xz, g_x_0_x_0_yy_xx_0_yy, g_x_0_x_0_yy_xx_0_yz, g_x_0_x_0_yy_xx_0_zz, g_xyy_xx_x_xx, g_xyy_xx_x_xy, g_xyy_xx_x_xz, g_xyy_xx_x_yy, g_xyy_xx_x_yz, g_xyy_xx_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_yy_xx_0_xx[i] = 4.0 * g_xyy_xx_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_xx_0_xy[i] = 4.0 * g_xyy_xx_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_xx_0_xz[i] = 4.0 * g_xyy_xx_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_xx_0_yy[i] = 4.0 * g_xyy_xx_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_xx_0_yz[i] = 4.0 * g_xyy_xx_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_xx_0_zz[i] = 4.0 * g_xyy_xx_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_x_0_x_0_yy_xy_0_xx, g_x_0_x_0_yy_xy_0_xy, g_x_0_x_0_yy_xy_0_xz, g_x_0_x_0_yy_xy_0_yy, g_x_0_x_0_yy_xy_0_yz, g_x_0_x_0_yy_xy_0_zz, g_xyy_xy_x_xx, g_xyy_xy_x_xy, g_xyy_xy_x_xz, g_xyy_xy_x_yy, g_xyy_xy_x_yz, g_xyy_xy_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_yy_xy_0_xx[i] = 4.0 * g_xyy_xy_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_xy_0_xy[i] = 4.0 * g_xyy_xy_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_xy_0_xz[i] = 4.0 * g_xyy_xy_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_xy_0_yy[i] = 4.0 * g_xyy_xy_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_xy_0_yz[i] = 4.0 * g_xyy_xy_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_xy_0_zz[i] = 4.0 * g_xyy_xy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_x_0_x_0_yy_xz_0_xx, g_x_0_x_0_yy_xz_0_xy, g_x_0_x_0_yy_xz_0_xz, g_x_0_x_0_yy_xz_0_yy, g_x_0_x_0_yy_xz_0_yz, g_x_0_x_0_yy_xz_0_zz, g_xyy_xz_x_xx, g_xyy_xz_x_xy, g_xyy_xz_x_xz, g_xyy_xz_x_yy, g_xyy_xz_x_yz, g_xyy_xz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_yy_xz_0_xx[i] = 4.0 * g_xyy_xz_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_xz_0_xy[i] = 4.0 * g_xyy_xz_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_xz_0_xz[i] = 4.0 * g_xyy_xz_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_xz_0_yy[i] = 4.0 * g_xyy_xz_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_xz_0_yz[i] = 4.0 * g_xyy_xz_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_xz_0_zz[i] = 4.0 * g_xyy_xz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_x_0_x_0_yy_yy_0_xx, g_x_0_x_0_yy_yy_0_xy, g_x_0_x_0_yy_yy_0_xz, g_x_0_x_0_yy_yy_0_yy, g_x_0_x_0_yy_yy_0_yz, g_x_0_x_0_yy_yy_0_zz, g_xyy_yy_x_xx, g_xyy_yy_x_xy, g_xyy_yy_x_xz, g_xyy_yy_x_yy, g_xyy_yy_x_yz, g_xyy_yy_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_yy_yy_0_xx[i] = 4.0 * g_xyy_yy_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_yy_0_xy[i] = 4.0 * g_xyy_yy_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_yy_0_xz[i] = 4.0 * g_xyy_yy_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_yy_0_yy[i] = 4.0 * g_xyy_yy_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_yy_0_yz[i] = 4.0 * g_xyy_yy_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_yy_0_zz[i] = 4.0 * g_xyy_yy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_x_0_x_0_yy_yz_0_xx, g_x_0_x_0_yy_yz_0_xy, g_x_0_x_0_yy_yz_0_xz, g_x_0_x_0_yy_yz_0_yy, g_x_0_x_0_yy_yz_0_yz, g_x_0_x_0_yy_yz_0_zz, g_xyy_yz_x_xx, g_xyy_yz_x_xy, g_xyy_yz_x_xz, g_xyy_yz_x_yy, g_xyy_yz_x_yz, g_xyy_yz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_yy_yz_0_xx[i] = 4.0 * g_xyy_yz_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_yz_0_xy[i] = 4.0 * g_xyy_yz_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_yz_0_xz[i] = 4.0 * g_xyy_yz_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_yz_0_yy[i] = 4.0 * g_xyy_yz_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_yz_0_yz[i] = 4.0 * g_xyy_yz_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_yz_0_zz[i] = 4.0 * g_xyy_yz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_x_0_x_0_yy_zz_0_xx, g_x_0_x_0_yy_zz_0_xy, g_x_0_x_0_yy_zz_0_xz, g_x_0_x_0_yy_zz_0_yy, g_x_0_x_0_yy_zz_0_yz, g_x_0_x_0_yy_zz_0_zz, g_xyy_zz_x_xx, g_xyy_zz_x_xy, g_xyy_zz_x_xz, g_xyy_zz_x_yy, g_xyy_zz_x_yz, g_xyy_zz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_yy_zz_0_xx[i] = 4.0 * g_xyy_zz_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_zz_0_xy[i] = 4.0 * g_xyy_zz_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_zz_0_xz[i] = 4.0 * g_xyy_zz_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_zz_0_yy[i] = 4.0 * g_xyy_zz_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_zz_0_yz[i] = 4.0 * g_xyy_zz_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_zz_0_zz[i] = 4.0 * g_xyy_zz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_x_0_x_0_yz_xx_0_xx, g_x_0_x_0_yz_xx_0_xy, g_x_0_x_0_yz_xx_0_xz, g_x_0_x_0_yz_xx_0_yy, g_x_0_x_0_yz_xx_0_yz, g_x_0_x_0_yz_xx_0_zz, g_xyz_xx_x_xx, g_xyz_xx_x_xy, g_xyz_xx_x_xz, g_xyz_xx_x_yy, g_xyz_xx_x_yz, g_xyz_xx_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_yz_xx_0_xx[i] = 4.0 * g_xyz_xx_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_xx_0_xy[i] = 4.0 * g_xyz_xx_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_xx_0_xz[i] = 4.0 * g_xyz_xx_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_xx_0_yy[i] = 4.0 * g_xyz_xx_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_xx_0_yz[i] = 4.0 * g_xyz_xx_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_xx_0_zz[i] = 4.0 * g_xyz_xx_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_x_0_x_0_yz_xy_0_xx, g_x_0_x_0_yz_xy_0_xy, g_x_0_x_0_yz_xy_0_xz, g_x_0_x_0_yz_xy_0_yy, g_x_0_x_0_yz_xy_0_yz, g_x_0_x_0_yz_xy_0_zz, g_xyz_xy_x_xx, g_xyz_xy_x_xy, g_xyz_xy_x_xz, g_xyz_xy_x_yy, g_xyz_xy_x_yz, g_xyz_xy_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_yz_xy_0_xx[i] = 4.0 * g_xyz_xy_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_xy_0_xy[i] = 4.0 * g_xyz_xy_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_xy_0_xz[i] = 4.0 * g_xyz_xy_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_xy_0_yy[i] = 4.0 * g_xyz_xy_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_xy_0_yz[i] = 4.0 * g_xyz_xy_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_xy_0_zz[i] = 4.0 * g_xyz_xy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_x_0_x_0_yz_xz_0_xx, g_x_0_x_0_yz_xz_0_xy, g_x_0_x_0_yz_xz_0_xz, g_x_0_x_0_yz_xz_0_yy, g_x_0_x_0_yz_xz_0_yz, g_x_0_x_0_yz_xz_0_zz, g_xyz_xz_x_xx, g_xyz_xz_x_xy, g_xyz_xz_x_xz, g_xyz_xz_x_yy, g_xyz_xz_x_yz, g_xyz_xz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_yz_xz_0_xx[i] = 4.0 * g_xyz_xz_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_xz_0_xy[i] = 4.0 * g_xyz_xz_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_xz_0_xz[i] = 4.0 * g_xyz_xz_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_xz_0_yy[i] = 4.0 * g_xyz_xz_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_xz_0_yz[i] = 4.0 * g_xyz_xz_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_xz_0_zz[i] = 4.0 * g_xyz_xz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_x_0_x_0_yz_yy_0_xx, g_x_0_x_0_yz_yy_0_xy, g_x_0_x_0_yz_yy_0_xz, g_x_0_x_0_yz_yy_0_yy, g_x_0_x_0_yz_yy_0_yz, g_x_0_x_0_yz_yy_0_zz, g_xyz_yy_x_xx, g_xyz_yy_x_xy, g_xyz_yy_x_xz, g_xyz_yy_x_yy, g_xyz_yy_x_yz, g_xyz_yy_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_yz_yy_0_xx[i] = 4.0 * g_xyz_yy_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_yy_0_xy[i] = 4.0 * g_xyz_yy_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_yy_0_xz[i] = 4.0 * g_xyz_yy_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_yy_0_yy[i] = 4.0 * g_xyz_yy_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_yy_0_yz[i] = 4.0 * g_xyz_yy_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_yy_0_zz[i] = 4.0 * g_xyz_yy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_x_0_x_0_yz_yz_0_xx, g_x_0_x_0_yz_yz_0_xy, g_x_0_x_0_yz_yz_0_xz, g_x_0_x_0_yz_yz_0_yy, g_x_0_x_0_yz_yz_0_yz, g_x_0_x_0_yz_yz_0_zz, g_xyz_yz_x_xx, g_xyz_yz_x_xy, g_xyz_yz_x_xz, g_xyz_yz_x_yy, g_xyz_yz_x_yz, g_xyz_yz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_yz_yz_0_xx[i] = 4.0 * g_xyz_yz_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_yz_0_xy[i] = 4.0 * g_xyz_yz_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_yz_0_xz[i] = 4.0 * g_xyz_yz_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_yz_0_yy[i] = 4.0 * g_xyz_yz_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_yz_0_yz[i] = 4.0 * g_xyz_yz_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_yz_0_zz[i] = 4.0 * g_xyz_yz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_x_0_x_0_yz_zz_0_xx, g_x_0_x_0_yz_zz_0_xy, g_x_0_x_0_yz_zz_0_xz, g_x_0_x_0_yz_zz_0_yy, g_x_0_x_0_yz_zz_0_yz, g_x_0_x_0_yz_zz_0_zz, g_xyz_zz_x_xx, g_xyz_zz_x_xy, g_xyz_zz_x_xz, g_xyz_zz_x_yy, g_xyz_zz_x_yz, g_xyz_zz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_yz_zz_0_xx[i] = 4.0 * g_xyz_zz_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_zz_0_xy[i] = 4.0 * g_xyz_zz_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_zz_0_xz[i] = 4.0 * g_xyz_zz_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_zz_0_yy[i] = 4.0 * g_xyz_zz_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_zz_0_yz[i] = 4.0 * g_xyz_zz_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_zz_0_zz[i] = 4.0 * g_xyz_zz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_x_0_x_0_zz_xx_0_xx, g_x_0_x_0_zz_xx_0_xy, g_x_0_x_0_zz_xx_0_xz, g_x_0_x_0_zz_xx_0_yy, g_x_0_x_0_zz_xx_0_yz, g_x_0_x_0_zz_xx_0_zz, g_xzz_xx_x_xx, g_xzz_xx_x_xy, g_xzz_xx_x_xz, g_xzz_xx_x_yy, g_xzz_xx_x_yz, g_xzz_xx_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_zz_xx_0_xx[i] = 4.0 * g_xzz_xx_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_xx_0_xy[i] = 4.0 * g_xzz_xx_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_xx_0_xz[i] = 4.0 * g_xzz_xx_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_xx_0_yy[i] = 4.0 * g_xzz_xx_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_xx_0_yz[i] = 4.0 * g_xzz_xx_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_xx_0_zz[i] = 4.0 * g_xzz_xx_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_x_0_x_0_zz_xy_0_xx, g_x_0_x_0_zz_xy_0_xy, g_x_0_x_0_zz_xy_0_xz, g_x_0_x_0_zz_xy_0_yy, g_x_0_x_0_zz_xy_0_yz, g_x_0_x_0_zz_xy_0_zz, g_xzz_xy_x_xx, g_xzz_xy_x_xy, g_xzz_xy_x_xz, g_xzz_xy_x_yy, g_xzz_xy_x_yz, g_xzz_xy_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_zz_xy_0_xx[i] = 4.0 * g_xzz_xy_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_xy_0_xy[i] = 4.0 * g_xzz_xy_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_xy_0_xz[i] = 4.0 * g_xzz_xy_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_xy_0_yy[i] = 4.0 * g_xzz_xy_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_xy_0_yz[i] = 4.0 * g_xzz_xy_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_xy_0_zz[i] = 4.0 * g_xzz_xy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_x_0_x_0_zz_xz_0_xx, g_x_0_x_0_zz_xz_0_xy, g_x_0_x_0_zz_xz_0_xz, g_x_0_x_0_zz_xz_0_yy, g_x_0_x_0_zz_xz_0_yz, g_x_0_x_0_zz_xz_0_zz, g_xzz_xz_x_xx, g_xzz_xz_x_xy, g_xzz_xz_x_xz, g_xzz_xz_x_yy, g_xzz_xz_x_yz, g_xzz_xz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_zz_xz_0_xx[i] = 4.0 * g_xzz_xz_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_xz_0_xy[i] = 4.0 * g_xzz_xz_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_xz_0_xz[i] = 4.0 * g_xzz_xz_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_xz_0_yy[i] = 4.0 * g_xzz_xz_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_xz_0_yz[i] = 4.0 * g_xzz_xz_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_xz_0_zz[i] = 4.0 * g_xzz_xz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_x_0_x_0_zz_yy_0_xx, g_x_0_x_0_zz_yy_0_xy, g_x_0_x_0_zz_yy_0_xz, g_x_0_x_0_zz_yy_0_yy, g_x_0_x_0_zz_yy_0_yz, g_x_0_x_0_zz_yy_0_zz, g_xzz_yy_x_xx, g_xzz_yy_x_xy, g_xzz_yy_x_xz, g_xzz_yy_x_yy, g_xzz_yy_x_yz, g_xzz_yy_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_zz_yy_0_xx[i] = 4.0 * g_xzz_yy_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_yy_0_xy[i] = 4.0 * g_xzz_yy_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_yy_0_xz[i] = 4.0 * g_xzz_yy_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_yy_0_yy[i] = 4.0 * g_xzz_yy_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_yy_0_yz[i] = 4.0 * g_xzz_yy_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_yy_0_zz[i] = 4.0 * g_xzz_yy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_x_0_x_0_zz_yz_0_xx, g_x_0_x_0_zz_yz_0_xy, g_x_0_x_0_zz_yz_0_xz, g_x_0_x_0_zz_yz_0_yy, g_x_0_x_0_zz_yz_0_yz, g_x_0_x_0_zz_yz_0_zz, g_xzz_yz_x_xx, g_xzz_yz_x_xy, g_xzz_yz_x_xz, g_xzz_yz_x_yy, g_xzz_yz_x_yz, g_xzz_yz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_zz_yz_0_xx[i] = 4.0 * g_xzz_yz_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_yz_0_xy[i] = 4.0 * g_xzz_yz_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_yz_0_xz[i] = 4.0 * g_xzz_yz_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_yz_0_yy[i] = 4.0 * g_xzz_yz_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_yz_0_yz[i] = 4.0 * g_xzz_yz_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_yz_0_zz[i] = 4.0 * g_xzz_yz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_x_0_x_0_zz_zz_0_xx, g_x_0_x_0_zz_zz_0_xy, g_x_0_x_0_zz_zz_0_xz, g_x_0_x_0_zz_zz_0_yy, g_x_0_x_0_zz_zz_0_yz, g_x_0_x_0_zz_zz_0_zz, g_xzz_zz_x_xx, g_xzz_zz_x_xy, g_xzz_zz_x_xz, g_xzz_zz_x_yy, g_xzz_zz_x_yz, g_xzz_zz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_zz_zz_0_xx[i] = 4.0 * g_xzz_zz_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_zz_0_xy[i] = 4.0 * g_xzz_zz_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_zz_0_xz[i] = 4.0 * g_xzz_zz_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_zz_0_yy[i] = 4.0 * g_xzz_zz_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_zz_0_yz[i] = 4.0 * g_xzz_zz_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_zz_0_zz[i] = 4.0 * g_xzz_zz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_x_0_y_0_xx_xx_0_xx, g_x_0_y_0_xx_xx_0_xy, g_x_0_y_0_xx_xx_0_xz, g_x_0_y_0_xx_xx_0_yy, g_x_0_y_0_xx_xx_0_yz, g_x_0_y_0_xx_xx_0_zz, g_x_xx_y_xx, g_x_xx_y_xy, g_x_xx_y_xz, g_x_xx_y_yy, g_x_xx_y_yz, g_x_xx_y_zz, g_xxx_xx_y_xx, g_xxx_xx_y_xy, g_xxx_xx_y_xz, g_xxx_xx_y_yy, g_xxx_xx_y_yz, g_xxx_xx_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_xx_xx_0_xx[i] = -4.0 * g_x_xx_y_xx[i] * c_exps[i] + 4.0 * g_xxx_xx_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_xx_0_xy[i] = -4.0 * g_x_xx_y_xy[i] * c_exps[i] + 4.0 * g_xxx_xx_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_xx_0_xz[i] = -4.0 * g_x_xx_y_xz[i] * c_exps[i] + 4.0 * g_xxx_xx_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_xx_0_yy[i] = -4.0 * g_x_xx_y_yy[i] * c_exps[i] + 4.0 * g_xxx_xx_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_xx_0_yz[i] = -4.0 * g_x_xx_y_yz[i] * c_exps[i] + 4.0 * g_xxx_xx_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_xx_0_zz[i] = -4.0 * g_x_xx_y_zz[i] * c_exps[i] + 4.0 * g_xxx_xx_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_x_0_y_0_xx_xy_0_xx, g_x_0_y_0_xx_xy_0_xy, g_x_0_y_0_xx_xy_0_xz, g_x_0_y_0_xx_xy_0_yy, g_x_0_y_0_xx_xy_0_yz, g_x_0_y_0_xx_xy_0_zz, g_x_xy_y_xx, g_x_xy_y_xy, g_x_xy_y_xz, g_x_xy_y_yy, g_x_xy_y_yz, g_x_xy_y_zz, g_xxx_xy_y_xx, g_xxx_xy_y_xy, g_xxx_xy_y_xz, g_xxx_xy_y_yy, g_xxx_xy_y_yz, g_xxx_xy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_xx_xy_0_xx[i] = -4.0 * g_x_xy_y_xx[i] * c_exps[i] + 4.0 * g_xxx_xy_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_xy_0_xy[i] = -4.0 * g_x_xy_y_xy[i] * c_exps[i] + 4.0 * g_xxx_xy_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_xy_0_xz[i] = -4.0 * g_x_xy_y_xz[i] * c_exps[i] + 4.0 * g_xxx_xy_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_xy_0_yy[i] = -4.0 * g_x_xy_y_yy[i] * c_exps[i] + 4.0 * g_xxx_xy_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_xy_0_yz[i] = -4.0 * g_x_xy_y_yz[i] * c_exps[i] + 4.0 * g_xxx_xy_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_xy_0_zz[i] = -4.0 * g_x_xy_y_zz[i] * c_exps[i] + 4.0 * g_xxx_xy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_x_0_y_0_xx_xz_0_xx, g_x_0_y_0_xx_xz_0_xy, g_x_0_y_0_xx_xz_0_xz, g_x_0_y_0_xx_xz_0_yy, g_x_0_y_0_xx_xz_0_yz, g_x_0_y_0_xx_xz_0_zz, g_x_xz_y_xx, g_x_xz_y_xy, g_x_xz_y_xz, g_x_xz_y_yy, g_x_xz_y_yz, g_x_xz_y_zz, g_xxx_xz_y_xx, g_xxx_xz_y_xy, g_xxx_xz_y_xz, g_xxx_xz_y_yy, g_xxx_xz_y_yz, g_xxx_xz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_xx_xz_0_xx[i] = -4.0 * g_x_xz_y_xx[i] * c_exps[i] + 4.0 * g_xxx_xz_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_xz_0_xy[i] = -4.0 * g_x_xz_y_xy[i] * c_exps[i] + 4.0 * g_xxx_xz_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_xz_0_xz[i] = -4.0 * g_x_xz_y_xz[i] * c_exps[i] + 4.0 * g_xxx_xz_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_xz_0_yy[i] = -4.0 * g_x_xz_y_yy[i] * c_exps[i] + 4.0 * g_xxx_xz_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_xz_0_yz[i] = -4.0 * g_x_xz_y_yz[i] * c_exps[i] + 4.0 * g_xxx_xz_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_xz_0_zz[i] = -4.0 * g_x_xz_y_zz[i] * c_exps[i] + 4.0 * g_xxx_xz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_x_0_y_0_xx_yy_0_xx, g_x_0_y_0_xx_yy_0_xy, g_x_0_y_0_xx_yy_0_xz, g_x_0_y_0_xx_yy_0_yy, g_x_0_y_0_xx_yy_0_yz, g_x_0_y_0_xx_yy_0_zz, g_x_yy_y_xx, g_x_yy_y_xy, g_x_yy_y_xz, g_x_yy_y_yy, g_x_yy_y_yz, g_x_yy_y_zz, g_xxx_yy_y_xx, g_xxx_yy_y_xy, g_xxx_yy_y_xz, g_xxx_yy_y_yy, g_xxx_yy_y_yz, g_xxx_yy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_xx_yy_0_xx[i] = -4.0 * g_x_yy_y_xx[i] * c_exps[i] + 4.0 * g_xxx_yy_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_yy_0_xy[i] = -4.0 * g_x_yy_y_xy[i] * c_exps[i] + 4.0 * g_xxx_yy_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_yy_0_xz[i] = -4.0 * g_x_yy_y_xz[i] * c_exps[i] + 4.0 * g_xxx_yy_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_yy_0_yy[i] = -4.0 * g_x_yy_y_yy[i] * c_exps[i] + 4.0 * g_xxx_yy_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_yy_0_yz[i] = -4.0 * g_x_yy_y_yz[i] * c_exps[i] + 4.0 * g_xxx_yy_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_yy_0_zz[i] = -4.0 * g_x_yy_y_zz[i] * c_exps[i] + 4.0 * g_xxx_yy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_x_0_y_0_xx_yz_0_xx, g_x_0_y_0_xx_yz_0_xy, g_x_0_y_0_xx_yz_0_xz, g_x_0_y_0_xx_yz_0_yy, g_x_0_y_0_xx_yz_0_yz, g_x_0_y_0_xx_yz_0_zz, g_x_yz_y_xx, g_x_yz_y_xy, g_x_yz_y_xz, g_x_yz_y_yy, g_x_yz_y_yz, g_x_yz_y_zz, g_xxx_yz_y_xx, g_xxx_yz_y_xy, g_xxx_yz_y_xz, g_xxx_yz_y_yy, g_xxx_yz_y_yz, g_xxx_yz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_xx_yz_0_xx[i] = -4.0 * g_x_yz_y_xx[i] * c_exps[i] + 4.0 * g_xxx_yz_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_yz_0_xy[i] = -4.0 * g_x_yz_y_xy[i] * c_exps[i] + 4.0 * g_xxx_yz_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_yz_0_xz[i] = -4.0 * g_x_yz_y_xz[i] * c_exps[i] + 4.0 * g_xxx_yz_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_yz_0_yy[i] = -4.0 * g_x_yz_y_yy[i] * c_exps[i] + 4.0 * g_xxx_yz_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_yz_0_yz[i] = -4.0 * g_x_yz_y_yz[i] * c_exps[i] + 4.0 * g_xxx_yz_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_yz_0_zz[i] = -4.0 * g_x_yz_y_zz[i] * c_exps[i] + 4.0 * g_xxx_yz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_x_0_y_0_xx_zz_0_xx, g_x_0_y_0_xx_zz_0_xy, g_x_0_y_0_xx_zz_0_xz, g_x_0_y_0_xx_zz_0_yy, g_x_0_y_0_xx_zz_0_yz, g_x_0_y_0_xx_zz_0_zz, g_x_zz_y_xx, g_x_zz_y_xy, g_x_zz_y_xz, g_x_zz_y_yy, g_x_zz_y_yz, g_x_zz_y_zz, g_xxx_zz_y_xx, g_xxx_zz_y_xy, g_xxx_zz_y_xz, g_xxx_zz_y_yy, g_xxx_zz_y_yz, g_xxx_zz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_xx_zz_0_xx[i] = -4.0 * g_x_zz_y_xx[i] * c_exps[i] + 4.0 * g_xxx_zz_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_zz_0_xy[i] = -4.0 * g_x_zz_y_xy[i] * c_exps[i] + 4.0 * g_xxx_zz_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_zz_0_xz[i] = -4.0 * g_x_zz_y_xz[i] * c_exps[i] + 4.0 * g_xxx_zz_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_zz_0_yy[i] = -4.0 * g_x_zz_y_yy[i] * c_exps[i] + 4.0 * g_xxx_zz_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_zz_0_yz[i] = -4.0 * g_x_zz_y_yz[i] * c_exps[i] + 4.0 * g_xxx_zz_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_zz_0_zz[i] = -4.0 * g_x_zz_y_zz[i] * c_exps[i] + 4.0 * g_xxx_zz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_x_0_y_0_xy_xx_0_xx, g_x_0_y_0_xy_xx_0_xy, g_x_0_y_0_xy_xx_0_xz, g_x_0_y_0_xy_xx_0_yy, g_x_0_y_0_xy_xx_0_yz, g_x_0_y_0_xy_xx_0_zz, g_xxy_xx_y_xx, g_xxy_xx_y_xy, g_xxy_xx_y_xz, g_xxy_xx_y_yy, g_xxy_xx_y_yz, g_xxy_xx_y_zz, g_y_xx_y_xx, g_y_xx_y_xy, g_y_xx_y_xz, g_y_xx_y_yy, g_y_xx_y_yz, g_y_xx_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_xy_xx_0_xx[i] = -2.0 * g_y_xx_y_xx[i] * c_exps[i] + 4.0 * g_xxy_xx_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_xx_0_xy[i] = -2.0 * g_y_xx_y_xy[i] * c_exps[i] + 4.0 * g_xxy_xx_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_xx_0_xz[i] = -2.0 * g_y_xx_y_xz[i] * c_exps[i] + 4.0 * g_xxy_xx_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_xx_0_yy[i] = -2.0 * g_y_xx_y_yy[i] * c_exps[i] + 4.0 * g_xxy_xx_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_xx_0_yz[i] = -2.0 * g_y_xx_y_yz[i] * c_exps[i] + 4.0 * g_xxy_xx_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_xx_0_zz[i] = -2.0 * g_y_xx_y_zz[i] * c_exps[i] + 4.0 * g_xxy_xx_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_x_0_y_0_xy_xy_0_xx, g_x_0_y_0_xy_xy_0_xy, g_x_0_y_0_xy_xy_0_xz, g_x_0_y_0_xy_xy_0_yy, g_x_0_y_0_xy_xy_0_yz, g_x_0_y_0_xy_xy_0_zz, g_xxy_xy_y_xx, g_xxy_xy_y_xy, g_xxy_xy_y_xz, g_xxy_xy_y_yy, g_xxy_xy_y_yz, g_xxy_xy_y_zz, g_y_xy_y_xx, g_y_xy_y_xy, g_y_xy_y_xz, g_y_xy_y_yy, g_y_xy_y_yz, g_y_xy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_xy_xy_0_xx[i] = -2.0 * g_y_xy_y_xx[i] * c_exps[i] + 4.0 * g_xxy_xy_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_xy_0_xy[i] = -2.0 * g_y_xy_y_xy[i] * c_exps[i] + 4.0 * g_xxy_xy_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_xy_0_xz[i] = -2.0 * g_y_xy_y_xz[i] * c_exps[i] + 4.0 * g_xxy_xy_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_xy_0_yy[i] = -2.0 * g_y_xy_y_yy[i] * c_exps[i] + 4.0 * g_xxy_xy_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_xy_0_yz[i] = -2.0 * g_y_xy_y_yz[i] * c_exps[i] + 4.0 * g_xxy_xy_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_xy_0_zz[i] = -2.0 * g_y_xy_y_zz[i] * c_exps[i] + 4.0 * g_xxy_xy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_x_0_y_0_xy_xz_0_xx, g_x_0_y_0_xy_xz_0_xy, g_x_0_y_0_xy_xz_0_xz, g_x_0_y_0_xy_xz_0_yy, g_x_0_y_0_xy_xz_0_yz, g_x_0_y_0_xy_xz_0_zz, g_xxy_xz_y_xx, g_xxy_xz_y_xy, g_xxy_xz_y_xz, g_xxy_xz_y_yy, g_xxy_xz_y_yz, g_xxy_xz_y_zz, g_y_xz_y_xx, g_y_xz_y_xy, g_y_xz_y_xz, g_y_xz_y_yy, g_y_xz_y_yz, g_y_xz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_xy_xz_0_xx[i] = -2.0 * g_y_xz_y_xx[i] * c_exps[i] + 4.0 * g_xxy_xz_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_xz_0_xy[i] = -2.0 * g_y_xz_y_xy[i] * c_exps[i] + 4.0 * g_xxy_xz_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_xz_0_xz[i] = -2.0 * g_y_xz_y_xz[i] * c_exps[i] + 4.0 * g_xxy_xz_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_xz_0_yy[i] = -2.0 * g_y_xz_y_yy[i] * c_exps[i] + 4.0 * g_xxy_xz_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_xz_0_yz[i] = -2.0 * g_y_xz_y_yz[i] * c_exps[i] + 4.0 * g_xxy_xz_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_xz_0_zz[i] = -2.0 * g_y_xz_y_zz[i] * c_exps[i] + 4.0 * g_xxy_xz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_x_0_y_0_xy_yy_0_xx, g_x_0_y_0_xy_yy_0_xy, g_x_0_y_0_xy_yy_0_xz, g_x_0_y_0_xy_yy_0_yy, g_x_0_y_0_xy_yy_0_yz, g_x_0_y_0_xy_yy_0_zz, g_xxy_yy_y_xx, g_xxy_yy_y_xy, g_xxy_yy_y_xz, g_xxy_yy_y_yy, g_xxy_yy_y_yz, g_xxy_yy_y_zz, g_y_yy_y_xx, g_y_yy_y_xy, g_y_yy_y_xz, g_y_yy_y_yy, g_y_yy_y_yz, g_y_yy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_xy_yy_0_xx[i] = -2.0 * g_y_yy_y_xx[i] * c_exps[i] + 4.0 * g_xxy_yy_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_yy_0_xy[i] = -2.0 * g_y_yy_y_xy[i] * c_exps[i] + 4.0 * g_xxy_yy_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_yy_0_xz[i] = -2.0 * g_y_yy_y_xz[i] * c_exps[i] + 4.0 * g_xxy_yy_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_yy_0_yy[i] = -2.0 * g_y_yy_y_yy[i] * c_exps[i] + 4.0 * g_xxy_yy_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_yy_0_yz[i] = -2.0 * g_y_yy_y_yz[i] * c_exps[i] + 4.0 * g_xxy_yy_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_yy_0_zz[i] = -2.0 * g_y_yy_y_zz[i] * c_exps[i] + 4.0 * g_xxy_yy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_x_0_y_0_xy_yz_0_xx, g_x_0_y_0_xy_yz_0_xy, g_x_0_y_0_xy_yz_0_xz, g_x_0_y_0_xy_yz_0_yy, g_x_0_y_0_xy_yz_0_yz, g_x_0_y_0_xy_yz_0_zz, g_xxy_yz_y_xx, g_xxy_yz_y_xy, g_xxy_yz_y_xz, g_xxy_yz_y_yy, g_xxy_yz_y_yz, g_xxy_yz_y_zz, g_y_yz_y_xx, g_y_yz_y_xy, g_y_yz_y_xz, g_y_yz_y_yy, g_y_yz_y_yz, g_y_yz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_xy_yz_0_xx[i] = -2.0 * g_y_yz_y_xx[i] * c_exps[i] + 4.0 * g_xxy_yz_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_yz_0_xy[i] = -2.0 * g_y_yz_y_xy[i] * c_exps[i] + 4.0 * g_xxy_yz_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_yz_0_xz[i] = -2.0 * g_y_yz_y_xz[i] * c_exps[i] + 4.0 * g_xxy_yz_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_yz_0_yy[i] = -2.0 * g_y_yz_y_yy[i] * c_exps[i] + 4.0 * g_xxy_yz_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_yz_0_yz[i] = -2.0 * g_y_yz_y_yz[i] * c_exps[i] + 4.0 * g_xxy_yz_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_yz_0_zz[i] = -2.0 * g_y_yz_y_zz[i] * c_exps[i] + 4.0 * g_xxy_yz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_x_0_y_0_xy_zz_0_xx, g_x_0_y_0_xy_zz_0_xy, g_x_0_y_0_xy_zz_0_xz, g_x_0_y_0_xy_zz_0_yy, g_x_0_y_0_xy_zz_0_yz, g_x_0_y_0_xy_zz_0_zz, g_xxy_zz_y_xx, g_xxy_zz_y_xy, g_xxy_zz_y_xz, g_xxy_zz_y_yy, g_xxy_zz_y_yz, g_xxy_zz_y_zz, g_y_zz_y_xx, g_y_zz_y_xy, g_y_zz_y_xz, g_y_zz_y_yy, g_y_zz_y_yz, g_y_zz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_xy_zz_0_xx[i] = -2.0 * g_y_zz_y_xx[i] * c_exps[i] + 4.0 * g_xxy_zz_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_zz_0_xy[i] = -2.0 * g_y_zz_y_xy[i] * c_exps[i] + 4.0 * g_xxy_zz_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_zz_0_xz[i] = -2.0 * g_y_zz_y_xz[i] * c_exps[i] + 4.0 * g_xxy_zz_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_zz_0_yy[i] = -2.0 * g_y_zz_y_yy[i] * c_exps[i] + 4.0 * g_xxy_zz_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_zz_0_yz[i] = -2.0 * g_y_zz_y_yz[i] * c_exps[i] + 4.0 * g_xxy_zz_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_zz_0_zz[i] = -2.0 * g_y_zz_y_zz[i] * c_exps[i] + 4.0 * g_xxy_zz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_x_0_y_0_xz_xx_0_xx, g_x_0_y_0_xz_xx_0_xy, g_x_0_y_0_xz_xx_0_xz, g_x_0_y_0_xz_xx_0_yy, g_x_0_y_0_xz_xx_0_yz, g_x_0_y_0_xz_xx_0_zz, g_xxz_xx_y_xx, g_xxz_xx_y_xy, g_xxz_xx_y_xz, g_xxz_xx_y_yy, g_xxz_xx_y_yz, g_xxz_xx_y_zz, g_z_xx_y_xx, g_z_xx_y_xy, g_z_xx_y_xz, g_z_xx_y_yy, g_z_xx_y_yz, g_z_xx_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_xz_xx_0_xx[i] = -2.0 * g_z_xx_y_xx[i] * c_exps[i] + 4.0 * g_xxz_xx_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_xx_0_xy[i] = -2.0 * g_z_xx_y_xy[i] * c_exps[i] + 4.0 * g_xxz_xx_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_xx_0_xz[i] = -2.0 * g_z_xx_y_xz[i] * c_exps[i] + 4.0 * g_xxz_xx_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_xx_0_yy[i] = -2.0 * g_z_xx_y_yy[i] * c_exps[i] + 4.0 * g_xxz_xx_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_xx_0_yz[i] = -2.0 * g_z_xx_y_yz[i] * c_exps[i] + 4.0 * g_xxz_xx_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_xx_0_zz[i] = -2.0 * g_z_xx_y_zz[i] * c_exps[i] + 4.0 * g_xxz_xx_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_x_0_y_0_xz_xy_0_xx, g_x_0_y_0_xz_xy_0_xy, g_x_0_y_0_xz_xy_0_xz, g_x_0_y_0_xz_xy_0_yy, g_x_0_y_0_xz_xy_0_yz, g_x_0_y_0_xz_xy_0_zz, g_xxz_xy_y_xx, g_xxz_xy_y_xy, g_xxz_xy_y_xz, g_xxz_xy_y_yy, g_xxz_xy_y_yz, g_xxz_xy_y_zz, g_z_xy_y_xx, g_z_xy_y_xy, g_z_xy_y_xz, g_z_xy_y_yy, g_z_xy_y_yz, g_z_xy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_xz_xy_0_xx[i] = -2.0 * g_z_xy_y_xx[i] * c_exps[i] + 4.0 * g_xxz_xy_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_xy_0_xy[i] = -2.0 * g_z_xy_y_xy[i] * c_exps[i] + 4.0 * g_xxz_xy_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_xy_0_xz[i] = -2.0 * g_z_xy_y_xz[i] * c_exps[i] + 4.0 * g_xxz_xy_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_xy_0_yy[i] = -2.0 * g_z_xy_y_yy[i] * c_exps[i] + 4.0 * g_xxz_xy_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_xy_0_yz[i] = -2.0 * g_z_xy_y_yz[i] * c_exps[i] + 4.0 * g_xxz_xy_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_xy_0_zz[i] = -2.0 * g_z_xy_y_zz[i] * c_exps[i] + 4.0 * g_xxz_xy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_x_0_y_0_xz_xz_0_xx, g_x_0_y_0_xz_xz_0_xy, g_x_0_y_0_xz_xz_0_xz, g_x_0_y_0_xz_xz_0_yy, g_x_0_y_0_xz_xz_0_yz, g_x_0_y_0_xz_xz_0_zz, g_xxz_xz_y_xx, g_xxz_xz_y_xy, g_xxz_xz_y_xz, g_xxz_xz_y_yy, g_xxz_xz_y_yz, g_xxz_xz_y_zz, g_z_xz_y_xx, g_z_xz_y_xy, g_z_xz_y_xz, g_z_xz_y_yy, g_z_xz_y_yz, g_z_xz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_xz_xz_0_xx[i] = -2.0 * g_z_xz_y_xx[i] * c_exps[i] + 4.0 * g_xxz_xz_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_xz_0_xy[i] = -2.0 * g_z_xz_y_xy[i] * c_exps[i] + 4.0 * g_xxz_xz_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_xz_0_xz[i] = -2.0 * g_z_xz_y_xz[i] * c_exps[i] + 4.0 * g_xxz_xz_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_xz_0_yy[i] = -2.0 * g_z_xz_y_yy[i] * c_exps[i] + 4.0 * g_xxz_xz_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_xz_0_yz[i] = -2.0 * g_z_xz_y_yz[i] * c_exps[i] + 4.0 * g_xxz_xz_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_xz_0_zz[i] = -2.0 * g_z_xz_y_zz[i] * c_exps[i] + 4.0 * g_xxz_xz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_x_0_y_0_xz_yy_0_xx, g_x_0_y_0_xz_yy_0_xy, g_x_0_y_0_xz_yy_0_xz, g_x_0_y_0_xz_yy_0_yy, g_x_0_y_0_xz_yy_0_yz, g_x_0_y_0_xz_yy_0_zz, g_xxz_yy_y_xx, g_xxz_yy_y_xy, g_xxz_yy_y_xz, g_xxz_yy_y_yy, g_xxz_yy_y_yz, g_xxz_yy_y_zz, g_z_yy_y_xx, g_z_yy_y_xy, g_z_yy_y_xz, g_z_yy_y_yy, g_z_yy_y_yz, g_z_yy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_xz_yy_0_xx[i] = -2.0 * g_z_yy_y_xx[i] * c_exps[i] + 4.0 * g_xxz_yy_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_yy_0_xy[i] = -2.0 * g_z_yy_y_xy[i] * c_exps[i] + 4.0 * g_xxz_yy_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_yy_0_xz[i] = -2.0 * g_z_yy_y_xz[i] * c_exps[i] + 4.0 * g_xxz_yy_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_yy_0_yy[i] = -2.0 * g_z_yy_y_yy[i] * c_exps[i] + 4.0 * g_xxz_yy_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_yy_0_yz[i] = -2.0 * g_z_yy_y_yz[i] * c_exps[i] + 4.0 * g_xxz_yy_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_yy_0_zz[i] = -2.0 * g_z_yy_y_zz[i] * c_exps[i] + 4.0 * g_xxz_yy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_x_0_y_0_xz_yz_0_xx, g_x_0_y_0_xz_yz_0_xy, g_x_0_y_0_xz_yz_0_xz, g_x_0_y_0_xz_yz_0_yy, g_x_0_y_0_xz_yz_0_yz, g_x_0_y_0_xz_yz_0_zz, g_xxz_yz_y_xx, g_xxz_yz_y_xy, g_xxz_yz_y_xz, g_xxz_yz_y_yy, g_xxz_yz_y_yz, g_xxz_yz_y_zz, g_z_yz_y_xx, g_z_yz_y_xy, g_z_yz_y_xz, g_z_yz_y_yy, g_z_yz_y_yz, g_z_yz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_xz_yz_0_xx[i] = -2.0 * g_z_yz_y_xx[i] * c_exps[i] + 4.0 * g_xxz_yz_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_yz_0_xy[i] = -2.0 * g_z_yz_y_xy[i] * c_exps[i] + 4.0 * g_xxz_yz_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_yz_0_xz[i] = -2.0 * g_z_yz_y_xz[i] * c_exps[i] + 4.0 * g_xxz_yz_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_yz_0_yy[i] = -2.0 * g_z_yz_y_yy[i] * c_exps[i] + 4.0 * g_xxz_yz_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_yz_0_yz[i] = -2.0 * g_z_yz_y_yz[i] * c_exps[i] + 4.0 * g_xxz_yz_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_yz_0_zz[i] = -2.0 * g_z_yz_y_zz[i] * c_exps[i] + 4.0 * g_xxz_yz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_x_0_y_0_xz_zz_0_xx, g_x_0_y_0_xz_zz_0_xy, g_x_0_y_0_xz_zz_0_xz, g_x_0_y_0_xz_zz_0_yy, g_x_0_y_0_xz_zz_0_yz, g_x_0_y_0_xz_zz_0_zz, g_xxz_zz_y_xx, g_xxz_zz_y_xy, g_xxz_zz_y_xz, g_xxz_zz_y_yy, g_xxz_zz_y_yz, g_xxz_zz_y_zz, g_z_zz_y_xx, g_z_zz_y_xy, g_z_zz_y_xz, g_z_zz_y_yy, g_z_zz_y_yz, g_z_zz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_xz_zz_0_xx[i] = -2.0 * g_z_zz_y_xx[i] * c_exps[i] + 4.0 * g_xxz_zz_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_zz_0_xy[i] = -2.0 * g_z_zz_y_xy[i] * c_exps[i] + 4.0 * g_xxz_zz_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_zz_0_xz[i] = -2.0 * g_z_zz_y_xz[i] * c_exps[i] + 4.0 * g_xxz_zz_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_zz_0_yy[i] = -2.0 * g_z_zz_y_yy[i] * c_exps[i] + 4.0 * g_xxz_zz_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_zz_0_yz[i] = -2.0 * g_z_zz_y_yz[i] * c_exps[i] + 4.0 * g_xxz_zz_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_zz_0_zz[i] = -2.0 * g_z_zz_y_zz[i] * c_exps[i] + 4.0 * g_xxz_zz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (324-330)

    #pragma omp simd aligned(g_x_0_y_0_yy_xx_0_xx, g_x_0_y_0_yy_xx_0_xy, g_x_0_y_0_yy_xx_0_xz, g_x_0_y_0_yy_xx_0_yy, g_x_0_y_0_yy_xx_0_yz, g_x_0_y_0_yy_xx_0_zz, g_xyy_xx_y_xx, g_xyy_xx_y_xy, g_xyy_xx_y_xz, g_xyy_xx_y_yy, g_xyy_xx_y_yz, g_xyy_xx_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_yy_xx_0_xx[i] = 4.0 * g_xyy_xx_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_xx_0_xy[i] = 4.0 * g_xyy_xx_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_xx_0_xz[i] = 4.0 * g_xyy_xx_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_xx_0_yy[i] = 4.0 * g_xyy_xx_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_xx_0_yz[i] = 4.0 * g_xyy_xx_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_xx_0_zz[i] = 4.0 * g_xyy_xx_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (330-336)

    #pragma omp simd aligned(g_x_0_y_0_yy_xy_0_xx, g_x_0_y_0_yy_xy_0_xy, g_x_0_y_0_yy_xy_0_xz, g_x_0_y_0_yy_xy_0_yy, g_x_0_y_0_yy_xy_0_yz, g_x_0_y_0_yy_xy_0_zz, g_xyy_xy_y_xx, g_xyy_xy_y_xy, g_xyy_xy_y_xz, g_xyy_xy_y_yy, g_xyy_xy_y_yz, g_xyy_xy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_yy_xy_0_xx[i] = 4.0 * g_xyy_xy_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_xy_0_xy[i] = 4.0 * g_xyy_xy_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_xy_0_xz[i] = 4.0 * g_xyy_xy_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_xy_0_yy[i] = 4.0 * g_xyy_xy_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_xy_0_yz[i] = 4.0 * g_xyy_xy_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_xy_0_zz[i] = 4.0 * g_xyy_xy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (336-342)

    #pragma omp simd aligned(g_x_0_y_0_yy_xz_0_xx, g_x_0_y_0_yy_xz_0_xy, g_x_0_y_0_yy_xz_0_xz, g_x_0_y_0_yy_xz_0_yy, g_x_0_y_0_yy_xz_0_yz, g_x_0_y_0_yy_xz_0_zz, g_xyy_xz_y_xx, g_xyy_xz_y_xy, g_xyy_xz_y_xz, g_xyy_xz_y_yy, g_xyy_xz_y_yz, g_xyy_xz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_yy_xz_0_xx[i] = 4.0 * g_xyy_xz_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_xz_0_xy[i] = 4.0 * g_xyy_xz_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_xz_0_xz[i] = 4.0 * g_xyy_xz_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_xz_0_yy[i] = 4.0 * g_xyy_xz_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_xz_0_yz[i] = 4.0 * g_xyy_xz_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_xz_0_zz[i] = 4.0 * g_xyy_xz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (342-348)

    #pragma omp simd aligned(g_x_0_y_0_yy_yy_0_xx, g_x_0_y_0_yy_yy_0_xy, g_x_0_y_0_yy_yy_0_xz, g_x_0_y_0_yy_yy_0_yy, g_x_0_y_0_yy_yy_0_yz, g_x_0_y_0_yy_yy_0_zz, g_xyy_yy_y_xx, g_xyy_yy_y_xy, g_xyy_yy_y_xz, g_xyy_yy_y_yy, g_xyy_yy_y_yz, g_xyy_yy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_yy_yy_0_xx[i] = 4.0 * g_xyy_yy_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_yy_0_xy[i] = 4.0 * g_xyy_yy_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_yy_0_xz[i] = 4.0 * g_xyy_yy_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_yy_0_yy[i] = 4.0 * g_xyy_yy_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_yy_0_yz[i] = 4.0 * g_xyy_yy_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_yy_0_zz[i] = 4.0 * g_xyy_yy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (348-354)

    #pragma omp simd aligned(g_x_0_y_0_yy_yz_0_xx, g_x_0_y_0_yy_yz_0_xy, g_x_0_y_0_yy_yz_0_xz, g_x_0_y_0_yy_yz_0_yy, g_x_0_y_0_yy_yz_0_yz, g_x_0_y_0_yy_yz_0_zz, g_xyy_yz_y_xx, g_xyy_yz_y_xy, g_xyy_yz_y_xz, g_xyy_yz_y_yy, g_xyy_yz_y_yz, g_xyy_yz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_yy_yz_0_xx[i] = 4.0 * g_xyy_yz_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_yz_0_xy[i] = 4.0 * g_xyy_yz_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_yz_0_xz[i] = 4.0 * g_xyy_yz_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_yz_0_yy[i] = 4.0 * g_xyy_yz_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_yz_0_yz[i] = 4.0 * g_xyy_yz_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_yz_0_zz[i] = 4.0 * g_xyy_yz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (354-360)

    #pragma omp simd aligned(g_x_0_y_0_yy_zz_0_xx, g_x_0_y_0_yy_zz_0_xy, g_x_0_y_0_yy_zz_0_xz, g_x_0_y_0_yy_zz_0_yy, g_x_0_y_0_yy_zz_0_yz, g_x_0_y_0_yy_zz_0_zz, g_xyy_zz_y_xx, g_xyy_zz_y_xy, g_xyy_zz_y_xz, g_xyy_zz_y_yy, g_xyy_zz_y_yz, g_xyy_zz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_yy_zz_0_xx[i] = 4.0 * g_xyy_zz_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_zz_0_xy[i] = 4.0 * g_xyy_zz_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_zz_0_xz[i] = 4.0 * g_xyy_zz_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_zz_0_yy[i] = 4.0 * g_xyy_zz_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_zz_0_yz[i] = 4.0 * g_xyy_zz_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_zz_0_zz[i] = 4.0 * g_xyy_zz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (360-366)

    #pragma omp simd aligned(g_x_0_y_0_yz_xx_0_xx, g_x_0_y_0_yz_xx_0_xy, g_x_0_y_0_yz_xx_0_xz, g_x_0_y_0_yz_xx_0_yy, g_x_0_y_0_yz_xx_0_yz, g_x_0_y_0_yz_xx_0_zz, g_xyz_xx_y_xx, g_xyz_xx_y_xy, g_xyz_xx_y_xz, g_xyz_xx_y_yy, g_xyz_xx_y_yz, g_xyz_xx_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_yz_xx_0_xx[i] = 4.0 * g_xyz_xx_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_xx_0_xy[i] = 4.0 * g_xyz_xx_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_xx_0_xz[i] = 4.0 * g_xyz_xx_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_xx_0_yy[i] = 4.0 * g_xyz_xx_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_xx_0_yz[i] = 4.0 * g_xyz_xx_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_xx_0_zz[i] = 4.0 * g_xyz_xx_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (366-372)

    #pragma omp simd aligned(g_x_0_y_0_yz_xy_0_xx, g_x_0_y_0_yz_xy_0_xy, g_x_0_y_0_yz_xy_0_xz, g_x_0_y_0_yz_xy_0_yy, g_x_0_y_0_yz_xy_0_yz, g_x_0_y_0_yz_xy_0_zz, g_xyz_xy_y_xx, g_xyz_xy_y_xy, g_xyz_xy_y_xz, g_xyz_xy_y_yy, g_xyz_xy_y_yz, g_xyz_xy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_yz_xy_0_xx[i] = 4.0 * g_xyz_xy_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_xy_0_xy[i] = 4.0 * g_xyz_xy_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_xy_0_xz[i] = 4.0 * g_xyz_xy_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_xy_0_yy[i] = 4.0 * g_xyz_xy_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_xy_0_yz[i] = 4.0 * g_xyz_xy_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_xy_0_zz[i] = 4.0 * g_xyz_xy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (372-378)

    #pragma omp simd aligned(g_x_0_y_0_yz_xz_0_xx, g_x_0_y_0_yz_xz_0_xy, g_x_0_y_0_yz_xz_0_xz, g_x_0_y_0_yz_xz_0_yy, g_x_0_y_0_yz_xz_0_yz, g_x_0_y_0_yz_xz_0_zz, g_xyz_xz_y_xx, g_xyz_xz_y_xy, g_xyz_xz_y_xz, g_xyz_xz_y_yy, g_xyz_xz_y_yz, g_xyz_xz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_yz_xz_0_xx[i] = 4.0 * g_xyz_xz_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_xz_0_xy[i] = 4.0 * g_xyz_xz_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_xz_0_xz[i] = 4.0 * g_xyz_xz_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_xz_0_yy[i] = 4.0 * g_xyz_xz_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_xz_0_yz[i] = 4.0 * g_xyz_xz_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_xz_0_zz[i] = 4.0 * g_xyz_xz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (378-384)

    #pragma omp simd aligned(g_x_0_y_0_yz_yy_0_xx, g_x_0_y_0_yz_yy_0_xy, g_x_0_y_0_yz_yy_0_xz, g_x_0_y_0_yz_yy_0_yy, g_x_0_y_0_yz_yy_0_yz, g_x_0_y_0_yz_yy_0_zz, g_xyz_yy_y_xx, g_xyz_yy_y_xy, g_xyz_yy_y_xz, g_xyz_yy_y_yy, g_xyz_yy_y_yz, g_xyz_yy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_yz_yy_0_xx[i] = 4.0 * g_xyz_yy_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_yy_0_xy[i] = 4.0 * g_xyz_yy_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_yy_0_xz[i] = 4.0 * g_xyz_yy_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_yy_0_yy[i] = 4.0 * g_xyz_yy_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_yy_0_yz[i] = 4.0 * g_xyz_yy_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_yy_0_zz[i] = 4.0 * g_xyz_yy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (384-390)

    #pragma omp simd aligned(g_x_0_y_0_yz_yz_0_xx, g_x_0_y_0_yz_yz_0_xy, g_x_0_y_0_yz_yz_0_xz, g_x_0_y_0_yz_yz_0_yy, g_x_0_y_0_yz_yz_0_yz, g_x_0_y_0_yz_yz_0_zz, g_xyz_yz_y_xx, g_xyz_yz_y_xy, g_xyz_yz_y_xz, g_xyz_yz_y_yy, g_xyz_yz_y_yz, g_xyz_yz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_yz_yz_0_xx[i] = 4.0 * g_xyz_yz_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_yz_0_xy[i] = 4.0 * g_xyz_yz_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_yz_0_xz[i] = 4.0 * g_xyz_yz_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_yz_0_yy[i] = 4.0 * g_xyz_yz_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_yz_0_yz[i] = 4.0 * g_xyz_yz_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_yz_0_zz[i] = 4.0 * g_xyz_yz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (390-396)

    #pragma omp simd aligned(g_x_0_y_0_yz_zz_0_xx, g_x_0_y_0_yz_zz_0_xy, g_x_0_y_0_yz_zz_0_xz, g_x_0_y_0_yz_zz_0_yy, g_x_0_y_0_yz_zz_0_yz, g_x_0_y_0_yz_zz_0_zz, g_xyz_zz_y_xx, g_xyz_zz_y_xy, g_xyz_zz_y_xz, g_xyz_zz_y_yy, g_xyz_zz_y_yz, g_xyz_zz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_yz_zz_0_xx[i] = 4.0 * g_xyz_zz_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_zz_0_xy[i] = 4.0 * g_xyz_zz_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_zz_0_xz[i] = 4.0 * g_xyz_zz_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_zz_0_yy[i] = 4.0 * g_xyz_zz_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_zz_0_yz[i] = 4.0 * g_xyz_zz_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_zz_0_zz[i] = 4.0 * g_xyz_zz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (396-402)

    #pragma omp simd aligned(g_x_0_y_0_zz_xx_0_xx, g_x_0_y_0_zz_xx_0_xy, g_x_0_y_0_zz_xx_0_xz, g_x_0_y_0_zz_xx_0_yy, g_x_0_y_0_zz_xx_0_yz, g_x_0_y_0_zz_xx_0_zz, g_xzz_xx_y_xx, g_xzz_xx_y_xy, g_xzz_xx_y_xz, g_xzz_xx_y_yy, g_xzz_xx_y_yz, g_xzz_xx_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_zz_xx_0_xx[i] = 4.0 * g_xzz_xx_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_xx_0_xy[i] = 4.0 * g_xzz_xx_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_xx_0_xz[i] = 4.0 * g_xzz_xx_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_xx_0_yy[i] = 4.0 * g_xzz_xx_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_xx_0_yz[i] = 4.0 * g_xzz_xx_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_xx_0_zz[i] = 4.0 * g_xzz_xx_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (402-408)

    #pragma omp simd aligned(g_x_0_y_0_zz_xy_0_xx, g_x_0_y_0_zz_xy_0_xy, g_x_0_y_0_zz_xy_0_xz, g_x_0_y_0_zz_xy_0_yy, g_x_0_y_0_zz_xy_0_yz, g_x_0_y_0_zz_xy_0_zz, g_xzz_xy_y_xx, g_xzz_xy_y_xy, g_xzz_xy_y_xz, g_xzz_xy_y_yy, g_xzz_xy_y_yz, g_xzz_xy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_zz_xy_0_xx[i] = 4.0 * g_xzz_xy_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_xy_0_xy[i] = 4.0 * g_xzz_xy_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_xy_0_xz[i] = 4.0 * g_xzz_xy_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_xy_0_yy[i] = 4.0 * g_xzz_xy_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_xy_0_yz[i] = 4.0 * g_xzz_xy_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_xy_0_zz[i] = 4.0 * g_xzz_xy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (408-414)

    #pragma omp simd aligned(g_x_0_y_0_zz_xz_0_xx, g_x_0_y_0_zz_xz_0_xy, g_x_0_y_0_zz_xz_0_xz, g_x_0_y_0_zz_xz_0_yy, g_x_0_y_0_zz_xz_0_yz, g_x_0_y_0_zz_xz_0_zz, g_xzz_xz_y_xx, g_xzz_xz_y_xy, g_xzz_xz_y_xz, g_xzz_xz_y_yy, g_xzz_xz_y_yz, g_xzz_xz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_zz_xz_0_xx[i] = 4.0 * g_xzz_xz_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_xz_0_xy[i] = 4.0 * g_xzz_xz_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_xz_0_xz[i] = 4.0 * g_xzz_xz_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_xz_0_yy[i] = 4.0 * g_xzz_xz_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_xz_0_yz[i] = 4.0 * g_xzz_xz_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_xz_0_zz[i] = 4.0 * g_xzz_xz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (414-420)

    #pragma omp simd aligned(g_x_0_y_0_zz_yy_0_xx, g_x_0_y_0_zz_yy_0_xy, g_x_0_y_0_zz_yy_0_xz, g_x_0_y_0_zz_yy_0_yy, g_x_0_y_0_zz_yy_0_yz, g_x_0_y_0_zz_yy_0_zz, g_xzz_yy_y_xx, g_xzz_yy_y_xy, g_xzz_yy_y_xz, g_xzz_yy_y_yy, g_xzz_yy_y_yz, g_xzz_yy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_zz_yy_0_xx[i] = 4.0 * g_xzz_yy_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_yy_0_xy[i] = 4.0 * g_xzz_yy_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_yy_0_xz[i] = 4.0 * g_xzz_yy_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_yy_0_yy[i] = 4.0 * g_xzz_yy_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_yy_0_yz[i] = 4.0 * g_xzz_yy_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_yy_0_zz[i] = 4.0 * g_xzz_yy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (420-426)

    #pragma omp simd aligned(g_x_0_y_0_zz_yz_0_xx, g_x_0_y_0_zz_yz_0_xy, g_x_0_y_0_zz_yz_0_xz, g_x_0_y_0_zz_yz_0_yy, g_x_0_y_0_zz_yz_0_yz, g_x_0_y_0_zz_yz_0_zz, g_xzz_yz_y_xx, g_xzz_yz_y_xy, g_xzz_yz_y_xz, g_xzz_yz_y_yy, g_xzz_yz_y_yz, g_xzz_yz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_zz_yz_0_xx[i] = 4.0 * g_xzz_yz_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_yz_0_xy[i] = 4.0 * g_xzz_yz_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_yz_0_xz[i] = 4.0 * g_xzz_yz_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_yz_0_yy[i] = 4.0 * g_xzz_yz_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_yz_0_yz[i] = 4.0 * g_xzz_yz_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_yz_0_zz[i] = 4.0 * g_xzz_yz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (426-432)

    #pragma omp simd aligned(g_x_0_y_0_zz_zz_0_xx, g_x_0_y_0_zz_zz_0_xy, g_x_0_y_0_zz_zz_0_xz, g_x_0_y_0_zz_zz_0_yy, g_x_0_y_0_zz_zz_0_yz, g_x_0_y_0_zz_zz_0_zz, g_xzz_zz_y_xx, g_xzz_zz_y_xy, g_xzz_zz_y_xz, g_xzz_zz_y_yy, g_xzz_zz_y_yz, g_xzz_zz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_zz_zz_0_xx[i] = 4.0 * g_xzz_zz_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_zz_0_xy[i] = 4.0 * g_xzz_zz_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_zz_0_xz[i] = 4.0 * g_xzz_zz_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_zz_0_yy[i] = 4.0 * g_xzz_zz_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_zz_0_yz[i] = 4.0 * g_xzz_zz_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_zz_0_zz[i] = 4.0 * g_xzz_zz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (432-438)

    #pragma omp simd aligned(g_x_0_z_0_xx_xx_0_xx, g_x_0_z_0_xx_xx_0_xy, g_x_0_z_0_xx_xx_0_xz, g_x_0_z_0_xx_xx_0_yy, g_x_0_z_0_xx_xx_0_yz, g_x_0_z_0_xx_xx_0_zz, g_x_xx_z_xx, g_x_xx_z_xy, g_x_xx_z_xz, g_x_xx_z_yy, g_x_xx_z_yz, g_x_xx_z_zz, g_xxx_xx_z_xx, g_xxx_xx_z_xy, g_xxx_xx_z_xz, g_xxx_xx_z_yy, g_xxx_xx_z_yz, g_xxx_xx_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_xx_xx_0_xx[i] = -4.0 * g_x_xx_z_xx[i] * c_exps[i] + 4.0 * g_xxx_xx_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_xx_0_xy[i] = -4.0 * g_x_xx_z_xy[i] * c_exps[i] + 4.0 * g_xxx_xx_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_xx_0_xz[i] = -4.0 * g_x_xx_z_xz[i] * c_exps[i] + 4.0 * g_xxx_xx_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_xx_0_yy[i] = -4.0 * g_x_xx_z_yy[i] * c_exps[i] + 4.0 * g_xxx_xx_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_xx_0_yz[i] = -4.0 * g_x_xx_z_yz[i] * c_exps[i] + 4.0 * g_xxx_xx_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_xx_0_zz[i] = -4.0 * g_x_xx_z_zz[i] * c_exps[i] + 4.0 * g_xxx_xx_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (438-444)

    #pragma omp simd aligned(g_x_0_z_0_xx_xy_0_xx, g_x_0_z_0_xx_xy_0_xy, g_x_0_z_0_xx_xy_0_xz, g_x_0_z_0_xx_xy_0_yy, g_x_0_z_0_xx_xy_0_yz, g_x_0_z_0_xx_xy_0_zz, g_x_xy_z_xx, g_x_xy_z_xy, g_x_xy_z_xz, g_x_xy_z_yy, g_x_xy_z_yz, g_x_xy_z_zz, g_xxx_xy_z_xx, g_xxx_xy_z_xy, g_xxx_xy_z_xz, g_xxx_xy_z_yy, g_xxx_xy_z_yz, g_xxx_xy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_xx_xy_0_xx[i] = -4.0 * g_x_xy_z_xx[i] * c_exps[i] + 4.0 * g_xxx_xy_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_xy_0_xy[i] = -4.0 * g_x_xy_z_xy[i] * c_exps[i] + 4.0 * g_xxx_xy_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_xy_0_xz[i] = -4.0 * g_x_xy_z_xz[i] * c_exps[i] + 4.0 * g_xxx_xy_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_xy_0_yy[i] = -4.0 * g_x_xy_z_yy[i] * c_exps[i] + 4.0 * g_xxx_xy_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_xy_0_yz[i] = -4.0 * g_x_xy_z_yz[i] * c_exps[i] + 4.0 * g_xxx_xy_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_xy_0_zz[i] = -4.0 * g_x_xy_z_zz[i] * c_exps[i] + 4.0 * g_xxx_xy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (444-450)

    #pragma omp simd aligned(g_x_0_z_0_xx_xz_0_xx, g_x_0_z_0_xx_xz_0_xy, g_x_0_z_0_xx_xz_0_xz, g_x_0_z_0_xx_xz_0_yy, g_x_0_z_0_xx_xz_0_yz, g_x_0_z_0_xx_xz_0_zz, g_x_xz_z_xx, g_x_xz_z_xy, g_x_xz_z_xz, g_x_xz_z_yy, g_x_xz_z_yz, g_x_xz_z_zz, g_xxx_xz_z_xx, g_xxx_xz_z_xy, g_xxx_xz_z_xz, g_xxx_xz_z_yy, g_xxx_xz_z_yz, g_xxx_xz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_xx_xz_0_xx[i] = -4.0 * g_x_xz_z_xx[i] * c_exps[i] + 4.0 * g_xxx_xz_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_xz_0_xy[i] = -4.0 * g_x_xz_z_xy[i] * c_exps[i] + 4.0 * g_xxx_xz_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_xz_0_xz[i] = -4.0 * g_x_xz_z_xz[i] * c_exps[i] + 4.0 * g_xxx_xz_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_xz_0_yy[i] = -4.0 * g_x_xz_z_yy[i] * c_exps[i] + 4.0 * g_xxx_xz_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_xz_0_yz[i] = -4.0 * g_x_xz_z_yz[i] * c_exps[i] + 4.0 * g_xxx_xz_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_xz_0_zz[i] = -4.0 * g_x_xz_z_zz[i] * c_exps[i] + 4.0 * g_xxx_xz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (450-456)

    #pragma omp simd aligned(g_x_0_z_0_xx_yy_0_xx, g_x_0_z_0_xx_yy_0_xy, g_x_0_z_0_xx_yy_0_xz, g_x_0_z_0_xx_yy_0_yy, g_x_0_z_0_xx_yy_0_yz, g_x_0_z_0_xx_yy_0_zz, g_x_yy_z_xx, g_x_yy_z_xy, g_x_yy_z_xz, g_x_yy_z_yy, g_x_yy_z_yz, g_x_yy_z_zz, g_xxx_yy_z_xx, g_xxx_yy_z_xy, g_xxx_yy_z_xz, g_xxx_yy_z_yy, g_xxx_yy_z_yz, g_xxx_yy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_xx_yy_0_xx[i] = -4.0 * g_x_yy_z_xx[i] * c_exps[i] + 4.0 * g_xxx_yy_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_yy_0_xy[i] = -4.0 * g_x_yy_z_xy[i] * c_exps[i] + 4.0 * g_xxx_yy_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_yy_0_xz[i] = -4.0 * g_x_yy_z_xz[i] * c_exps[i] + 4.0 * g_xxx_yy_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_yy_0_yy[i] = -4.0 * g_x_yy_z_yy[i] * c_exps[i] + 4.0 * g_xxx_yy_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_yy_0_yz[i] = -4.0 * g_x_yy_z_yz[i] * c_exps[i] + 4.0 * g_xxx_yy_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_yy_0_zz[i] = -4.0 * g_x_yy_z_zz[i] * c_exps[i] + 4.0 * g_xxx_yy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (456-462)

    #pragma omp simd aligned(g_x_0_z_0_xx_yz_0_xx, g_x_0_z_0_xx_yz_0_xy, g_x_0_z_0_xx_yz_0_xz, g_x_0_z_0_xx_yz_0_yy, g_x_0_z_0_xx_yz_0_yz, g_x_0_z_0_xx_yz_0_zz, g_x_yz_z_xx, g_x_yz_z_xy, g_x_yz_z_xz, g_x_yz_z_yy, g_x_yz_z_yz, g_x_yz_z_zz, g_xxx_yz_z_xx, g_xxx_yz_z_xy, g_xxx_yz_z_xz, g_xxx_yz_z_yy, g_xxx_yz_z_yz, g_xxx_yz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_xx_yz_0_xx[i] = -4.0 * g_x_yz_z_xx[i] * c_exps[i] + 4.0 * g_xxx_yz_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_yz_0_xy[i] = -4.0 * g_x_yz_z_xy[i] * c_exps[i] + 4.0 * g_xxx_yz_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_yz_0_xz[i] = -4.0 * g_x_yz_z_xz[i] * c_exps[i] + 4.0 * g_xxx_yz_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_yz_0_yy[i] = -4.0 * g_x_yz_z_yy[i] * c_exps[i] + 4.0 * g_xxx_yz_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_yz_0_yz[i] = -4.0 * g_x_yz_z_yz[i] * c_exps[i] + 4.0 * g_xxx_yz_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_yz_0_zz[i] = -4.0 * g_x_yz_z_zz[i] * c_exps[i] + 4.0 * g_xxx_yz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (462-468)

    #pragma omp simd aligned(g_x_0_z_0_xx_zz_0_xx, g_x_0_z_0_xx_zz_0_xy, g_x_0_z_0_xx_zz_0_xz, g_x_0_z_0_xx_zz_0_yy, g_x_0_z_0_xx_zz_0_yz, g_x_0_z_0_xx_zz_0_zz, g_x_zz_z_xx, g_x_zz_z_xy, g_x_zz_z_xz, g_x_zz_z_yy, g_x_zz_z_yz, g_x_zz_z_zz, g_xxx_zz_z_xx, g_xxx_zz_z_xy, g_xxx_zz_z_xz, g_xxx_zz_z_yy, g_xxx_zz_z_yz, g_xxx_zz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_xx_zz_0_xx[i] = -4.0 * g_x_zz_z_xx[i] * c_exps[i] + 4.0 * g_xxx_zz_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_zz_0_xy[i] = -4.0 * g_x_zz_z_xy[i] * c_exps[i] + 4.0 * g_xxx_zz_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_zz_0_xz[i] = -4.0 * g_x_zz_z_xz[i] * c_exps[i] + 4.0 * g_xxx_zz_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_zz_0_yy[i] = -4.0 * g_x_zz_z_yy[i] * c_exps[i] + 4.0 * g_xxx_zz_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_zz_0_yz[i] = -4.0 * g_x_zz_z_yz[i] * c_exps[i] + 4.0 * g_xxx_zz_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_zz_0_zz[i] = -4.0 * g_x_zz_z_zz[i] * c_exps[i] + 4.0 * g_xxx_zz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (468-474)

    #pragma omp simd aligned(g_x_0_z_0_xy_xx_0_xx, g_x_0_z_0_xy_xx_0_xy, g_x_0_z_0_xy_xx_0_xz, g_x_0_z_0_xy_xx_0_yy, g_x_0_z_0_xy_xx_0_yz, g_x_0_z_0_xy_xx_0_zz, g_xxy_xx_z_xx, g_xxy_xx_z_xy, g_xxy_xx_z_xz, g_xxy_xx_z_yy, g_xxy_xx_z_yz, g_xxy_xx_z_zz, g_y_xx_z_xx, g_y_xx_z_xy, g_y_xx_z_xz, g_y_xx_z_yy, g_y_xx_z_yz, g_y_xx_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_xy_xx_0_xx[i] = -2.0 * g_y_xx_z_xx[i] * c_exps[i] + 4.0 * g_xxy_xx_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_xx_0_xy[i] = -2.0 * g_y_xx_z_xy[i] * c_exps[i] + 4.0 * g_xxy_xx_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_xx_0_xz[i] = -2.0 * g_y_xx_z_xz[i] * c_exps[i] + 4.0 * g_xxy_xx_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_xx_0_yy[i] = -2.0 * g_y_xx_z_yy[i] * c_exps[i] + 4.0 * g_xxy_xx_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_xx_0_yz[i] = -2.0 * g_y_xx_z_yz[i] * c_exps[i] + 4.0 * g_xxy_xx_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_xx_0_zz[i] = -2.0 * g_y_xx_z_zz[i] * c_exps[i] + 4.0 * g_xxy_xx_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (474-480)

    #pragma omp simd aligned(g_x_0_z_0_xy_xy_0_xx, g_x_0_z_0_xy_xy_0_xy, g_x_0_z_0_xy_xy_0_xz, g_x_0_z_0_xy_xy_0_yy, g_x_0_z_0_xy_xy_0_yz, g_x_0_z_0_xy_xy_0_zz, g_xxy_xy_z_xx, g_xxy_xy_z_xy, g_xxy_xy_z_xz, g_xxy_xy_z_yy, g_xxy_xy_z_yz, g_xxy_xy_z_zz, g_y_xy_z_xx, g_y_xy_z_xy, g_y_xy_z_xz, g_y_xy_z_yy, g_y_xy_z_yz, g_y_xy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_xy_xy_0_xx[i] = -2.0 * g_y_xy_z_xx[i] * c_exps[i] + 4.0 * g_xxy_xy_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_xy_0_xy[i] = -2.0 * g_y_xy_z_xy[i] * c_exps[i] + 4.0 * g_xxy_xy_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_xy_0_xz[i] = -2.0 * g_y_xy_z_xz[i] * c_exps[i] + 4.0 * g_xxy_xy_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_xy_0_yy[i] = -2.0 * g_y_xy_z_yy[i] * c_exps[i] + 4.0 * g_xxy_xy_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_xy_0_yz[i] = -2.0 * g_y_xy_z_yz[i] * c_exps[i] + 4.0 * g_xxy_xy_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_xy_0_zz[i] = -2.0 * g_y_xy_z_zz[i] * c_exps[i] + 4.0 * g_xxy_xy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (480-486)

    #pragma omp simd aligned(g_x_0_z_0_xy_xz_0_xx, g_x_0_z_0_xy_xz_0_xy, g_x_0_z_0_xy_xz_0_xz, g_x_0_z_0_xy_xz_0_yy, g_x_0_z_0_xy_xz_0_yz, g_x_0_z_0_xy_xz_0_zz, g_xxy_xz_z_xx, g_xxy_xz_z_xy, g_xxy_xz_z_xz, g_xxy_xz_z_yy, g_xxy_xz_z_yz, g_xxy_xz_z_zz, g_y_xz_z_xx, g_y_xz_z_xy, g_y_xz_z_xz, g_y_xz_z_yy, g_y_xz_z_yz, g_y_xz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_xy_xz_0_xx[i] = -2.0 * g_y_xz_z_xx[i] * c_exps[i] + 4.0 * g_xxy_xz_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_xz_0_xy[i] = -2.0 * g_y_xz_z_xy[i] * c_exps[i] + 4.0 * g_xxy_xz_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_xz_0_xz[i] = -2.0 * g_y_xz_z_xz[i] * c_exps[i] + 4.0 * g_xxy_xz_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_xz_0_yy[i] = -2.0 * g_y_xz_z_yy[i] * c_exps[i] + 4.0 * g_xxy_xz_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_xz_0_yz[i] = -2.0 * g_y_xz_z_yz[i] * c_exps[i] + 4.0 * g_xxy_xz_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_xz_0_zz[i] = -2.0 * g_y_xz_z_zz[i] * c_exps[i] + 4.0 * g_xxy_xz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (486-492)

    #pragma omp simd aligned(g_x_0_z_0_xy_yy_0_xx, g_x_0_z_0_xy_yy_0_xy, g_x_0_z_0_xy_yy_0_xz, g_x_0_z_0_xy_yy_0_yy, g_x_0_z_0_xy_yy_0_yz, g_x_0_z_0_xy_yy_0_zz, g_xxy_yy_z_xx, g_xxy_yy_z_xy, g_xxy_yy_z_xz, g_xxy_yy_z_yy, g_xxy_yy_z_yz, g_xxy_yy_z_zz, g_y_yy_z_xx, g_y_yy_z_xy, g_y_yy_z_xz, g_y_yy_z_yy, g_y_yy_z_yz, g_y_yy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_xy_yy_0_xx[i] = -2.0 * g_y_yy_z_xx[i] * c_exps[i] + 4.0 * g_xxy_yy_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_yy_0_xy[i] = -2.0 * g_y_yy_z_xy[i] * c_exps[i] + 4.0 * g_xxy_yy_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_yy_0_xz[i] = -2.0 * g_y_yy_z_xz[i] * c_exps[i] + 4.0 * g_xxy_yy_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_yy_0_yy[i] = -2.0 * g_y_yy_z_yy[i] * c_exps[i] + 4.0 * g_xxy_yy_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_yy_0_yz[i] = -2.0 * g_y_yy_z_yz[i] * c_exps[i] + 4.0 * g_xxy_yy_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_yy_0_zz[i] = -2.0 * g_y_yy_z_zz[i] * c_exps[i] + 4.0 * g_xxy_yy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (492-498)

    #pragma omp simd aligned(g_x_0_z_0_xy_yz_0_xx, g_x_0_z_0_xy_yz_0_xy, g_x_0_z_0_xy_yz_0_xz, g_x_0_z_0_xy_yz_0_yy, g_x_0_z_0_xy_yz_0_yz, g_x_0_z_0_xy_yz_0_zz, g_xxy_yz_z_xx, g_xxy_yz_z_xy, g_xxy_yz_z_xz, g_xxy_yz_z_yy, g_xxy_yz_z_yz, g_xxy_yz_z_zz, g_y_yz_z_xx, g_y_yz_z_xy, g_y_yz_z_xz, g_y_yz_z_yy, g_y_yz_z_yz, g_y_yz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_xy_yz_0_xx[i] = -2.0 * g_y_yz_z_xx[i] * c_exps[i] + 4.0 * g_xxy_yz_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_yz_0_xy[i] = -2.0 * g_y_yz_z_xy[i] * c_exps[i] + 4.0 * g_xxy_yz_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_yz_0_xz[i] = -2.0 * g_y_yz_z_xz[i] * c_exps[i] + 4.0 * g_xxy_yz_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_yz_0_yy[i] = -2.0 * g_y_yz_z_yy[i] * c_exps[i] + 4.0 * g_xxy_yz_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_yz_0_yz[i] = -2.0 * g_y_yz_z_yz[i] * c_exps[i] + 4.0 * g_xxy_yz_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_yz_0_zz[i] = -2.0 * g_y_yz_z_zz[i] * c_exps[i] + 4.0 * g_xxy_yz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (498-504)

    #pragma omp simd aligned(g_x_0_z_0_xy_zz_0_xx, g_x_0_z_0_xy_zz_0_xy, g_x_0_z_0_xy_zz_0_xz, g_x_0_z_0_xy_zz_0_yy, g_x_0_z_0_xy_zz_0_yz, g_x_0_z_0_xy_zz_0_zz, g_xxy_zz_z_xx, g_xxy_zz_z_xy, g_xxy_zz_z_xz, g_xxy_zz_z_yy, g_xxy_zz_z_yz, g_xxy_zz_z_zz, g_y_zz_z_xx, g_y_zz_z_xy, g_y_zz_z_xz, g_y_zz_z_yy, g_y_zz_z_yz, g_y_zz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_xy_zz_0_xx[i] = -2.0 * g_y_zz_z_xx[i] * c_exps[i] + 4.0 * g_xxy_zz_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_zz_0_xy[i] = -2.0 * g_y_zz_z_xy[i] * c_exps[i] + 4.0 * g_xxy_zz_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_zz_0_xz[i] = -2.0 * g_y_zz_z_xz[i] * c_exps[i] + 4.0 * g_xxy_zz_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_zz_0_yy[i] = -2.0 * g_y_zz_z_yy[i] * c_exps[i] + 4.0 * g_xxy_zz_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_zz_0_yz[i] = -2.0 * g_y_zz_z_yz[i] * c_exps[i] + 4.0 * g_xxy_zz_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_zz_0_zz[i] = -2.0 * g_y_zz_z_zz[i] * c_exps[i] + 4.0 * g_xxy_zz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (504-510)

    #pragma omp simd aligned(g_x_0_z_0_xz_xx_0_xx, g_x_0_z_0_xz_xx_0_xy, g_x_0_z_0_xz_xx_0_xz, g_x_0_z_0_xz_xx_0_yy, g_x_0_z_0_xz_xx_0_yz, g_x_0_z_0_xz_xx_0_zz, g_xxz_xx_z_xx, g_xxz_xx_z_xy, g_xxz_xx_z_xz, g_xxz_xx_z_yy, g_xxz_xx_z_yz, g_xxz_xx_z_zz, g_z_xx_z_xx, g_z_xx_z_xy, g_z_xx_z_xz, g_z_xx_z_yy, g_z_xx_z_yz, g_z_xx_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_xz_xx_0_xx[i] = -2.0 * g_z_xx_z_xx[i] * c_exps[i] + 4.0 * g_xxz_xx_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_xx_0_xy[i] = -2.0 * g_z_xx_z_xy[i] * c_exps[i] + 4.0 * g_xxz_xx_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_xx_0_xz[i] = -2.0 * g_z_xx_z_xz[i] * c_exps[i] + 4.0 * g_xxz_xx_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_xx_0_yy[i] = -2.0 * g_z_xx_z_yy[i] * c_exps[i] + 4.0 * g_xxz_xx_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_xx_0_yz[i] = -2.0 * g_z_xx_z_yz[i] * c_exps[i] + 4.0 * g_xxz_xx_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_xx_0_zz[i] = -2.0 * g_z_xx_z_zz[i] * c_exps[i] + 4.0 * g_xxz_xx_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (510-516)

    #pragma omp simd aligned(g_x_0_z_0_xz_xy_0_xx, g_x_0_z_0_xz_xy_0_xy, g_x_0_z_0_xz_xy_0_xz, g_x_0_z_0_xz_xy_0_yy, g_x_0_z_0_xz_xy_0_yz, g_x_0_z_0_xz_xy_0_zz, g_xxz_xy_z_xx, g_xxz_xy_z_xy, g_xxz_xy_z_xz, g_xxz_xy_z_yy, g_xxz_xy_z_yz, g_xxz_xy_z_zz, g_z_xy_z_xx, g_z_xy_z_xy, g_z_xy_z_xz, g_z_xy_z_yy, g_z_xy_z_yz, g_z_xy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_xz_xy_0_xx[i] = -2.0 * g_z_xy_z_xx[i] * c_exps[i] + 4.0 * g_xxz_xy_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_xy_0_xy[i] = -2.0 * g_z_xy_z_xy[i] * c_exps[i] + 4.0 * g_xxz_xy_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_xy_0_xz[i] = -2.0 * g_z_xy_z_xz[i] * c_exps[i] + 4.0 * g_xxz_xy_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_xy_0_yy[i] = -2.0 * g_z_xy_z_yy[i] * c_exps[i] + 4.0 * g_xxz_xy_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_xy_0_yz[i] = -2.0 * g_z_xy_z_yz[i] * c_exps[i] + 4.0 * g_xxz_xy_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_xy_0_zz[i] = -2.0 * g_z_xy_z_zz[i] * c_exps[i] + 4.0 * g_xxz_xy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (516-522)

    #pragma omp simd aligned(g_x_0_z_0_xz_xz_0_xx, g_x_0_z_0_xz_xz_0_xy, g_x_0_z_0_xz_xz_0_xz, g_x_0_z_0_xz_xz_0_yy, g_x_0_z_0_xz_xz_0_yz, g_x_0_z_0_xz_xz_0_zz, g_xxz_xz_z_xx, g_xxz_xz_z_xy, g_xxz_xz_z_xz, g_xxz_xz_z_yy, g_xxz_xz_z_yz, g_xxz_xz_z_zz, g_z_xz_z_xx, g_z_xz_z_xy, g_z_xz_z_xz, g_z_xz_z_yy, g_z_xz_z_yz, g_z_xz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_xz_xz_0_xx[i] = -2.0 * g_z_xz_z_xx[i] * c_exps[i] + 4.0 * g_xxz_xz_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_xz_0_xy[i] = -2.0 * g_z_xz_z_xy[i] * c_exps[i] + 4.0 * g_xxz_xz_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_xz_0_xz[i] = -2.0 * g_z_xz_z_xz[i] * c_exps[i] + 4.0 * g_xxz_xz_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_xz_0_yy[i] = -2.0 * g_z_xz_z_yy[i] * c_exps[i] + 4.0 * g_xxz_xz_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_xz_0_yz[i] = -2.0 * g_z_xz_z_yz[i] * c_exps[i] + 4.0 * g_xxz_xz_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_xz_0_zz[i] = -2.0 * g_z_xz_z_zz[i] * c_exps[i] + 4.0 * g_xxz_xz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (522-528)

    #pragma omp simd aligned(g_x_0_z_0_xz_yy_0_xx, g_x_0_z_0_xz_yy_0_xy, g_x_0_z_0_xz_yy_0_xz, g_x_0_z_0_xz_yy_0_yy, g_x_0_z_0_xz_yy_0_yz, g_x_0_z_0_xz_yy_0_zz, g_xxz_yy_z_xx, g_xxz_yy_z_xy, g_xxz_yy_z_xz, g_xxz_yy_z_yy, g_xxz_yy_z_yz, g_xxz_yy_z_zz, g_z_yy_z_xx, g_z_yy_z_xy, g_z_yy_z_xz, g_z_yy_z_yy, g_z_yy_z_yz, g_z_yy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_xz_yy_0_xx[i] = -2.0 * g_z_yy_z_xx[i] * c_exps[i] + 4.0 * g_xxz_yy_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_yy_0_xy[i] = -2.0 * g_z_yy_z_xy[i] * c_exps[i] + 4.0 * g_xxz_yy_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_yy_0_xz[i] = -2.0 * g_z_yy_z_xz[i] * c_exps[i] + 4.0 * g_xxz_yy_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_yy_0_yy[i] = -2.0 * g_z_yy_z_yy[i] * c_exps[i] + 4.0 * g_xxz_yy_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_yy_0_yz[i] = -2.0 * g_z_yy_z_yz[i] * c_exps[i] + 4.0 * g_xxz_yy_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_yy_0_zz[i] = -2.0 * g_z_yy_z_zz[i] * c_exps[i] + 4.0 * g_xxz_yy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (528-534)

    #pragma omp simd aligned(g_x_0_z_0_xz_yz_0_xx, g_x_0_z_0_xz_yz_0_xy, g_x_0_z_0_xz_yz_0_xz, g_x_0_z_0_xz_yz_0_yy, g_x_0_z_0_xz_yz_0_yz, g_x_0_z_0_xz_yz_0_zz, g_xxz_yz_z_xx, g_xxz_yz_z_xy, g_xxz_yz_z_xz, g_xxz_yz_z_yy, g_xxz_yz_z_yz, g_xxz_yz_z_zz, g_z_yz_z_xx, g_z_yz_z_xy, g_z_yz_z_xz, g_z_yz_z_yy, g_z_yz_z_yz, g_z_yz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_xz_yz_0_xx[i] = -2.0 * g_z_yz_z_xx[i] * c_exps[i] + 4.0 * g_xxz_yz_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_yz_0_xy[i] = -2.0 * g_z_yz_z_xy[i] * c_exps[i] + 4.0 * g_xxz_yz_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_yz_0_xz[i] = -2.0 * g_z_yz_z_xz[i] * c_exps[i] + 4.0 * g_xxz_yz_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_yz_0_yy[i] = -2.0 * g_z_yz_z_yy[i] * c_exps[i] + 4.0 * g_xxz_yz_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_yz_0_yz[i] = -2.0 * g_z_yz_z_yz[i] * c_exps[i] + 4.0 * g_xxz_yz_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_yz_0_zz[i] = -2.0 * g_z_yz_z_zz[i] * c_exps[i] + 4.0 * g_xxz_yz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (534-540)

    #pragma omp simd aligned(g_x_0_z_0_xz_zz_0_xx, g_x_0_z_0_xz_zz_0_xy, g_x_0_z_0_xz_zz_0_xz, g_x_0_z_0_xz_zz_0_yy, g_x_0_z_0_xz_zz_0_yz, g_x_0_z_0_xz_zz_0_zz, g_xxz_zz_z_xx, g_xxz_zz_z_xy, g_xxz_zz_z_xz, g_xxz_zz_z_yy, g_xxz_zz_z_yz, g_xxz_zz_z_zz, g_z_zz_z_xx, g_z_zz_z_xy, g_z_zz_z_xz, g_z_zz_z_yy, g_z_zz_z_yz, g_z_zz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_xz_zz_0_xx[i] = -2.0 * g_z_zz_z_xx[i] * c_exps[i] + 4.0 * g_xxz_zz_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_zz_0_xy[i] = -2.0 * g_z_zz_z_xy[i] * c_exps[i] + 4.0 * g_xxz_zz_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_zz_0_xz[i] = -2.0 * g_z_zz_z_xz[i] * c_exps[i] + 4.0 * g_xxz_zz_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_zz_0_yy[i] = -2.0 * g_z_zz_z_yy[i] * c_exps[i] + 4.0 * g_xxz_zz_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_zz_0_yz[i] = -2.0 * g_z_zz_z_yz[i] * c_exps[i] + 4.0 * g_xxz_zz_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_zz_0_zz[i] = -2.0 * g_z_zz_z_zz[i] * c_exps[i] + 4.0 * g_xxz_zz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (540-546)

    #pragma omp simd aligned(g_x_0_z_0_yy_xx_0_xx, g_x_0_z_0_yy_xx_0_xy, g_x_0_z_0_yy_xx_0_xz, g_x_0_z_0_yy_xx_0_yy, g_x_0_z_0_yy_xx_0_yz, g_x_0_z_0_yy_xx_0_zz, g_xyy_xx_z_xx, g_xyy_xx_z_xy, g_xyy_xx_z_xz, g_xyy_xx_z_yy, g_xyy_xx_z_yz, g_xyy_xx_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_yy_xx_0_xx[i] = 4.0 * g_xyy_xx_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_xx_0_xy[i] = 4.0 * g_xyy_xx_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_xx_0_xz[i] = 4.0 * g_xyy_xx_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_xx_0_yy[i] = 4.0 * g_xyy_xx_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_xx_0_yz[i] = 4.0 * g_xyy_xx_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_xx_0_zz[i] = 4.0 * g_xyy_xx_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (546-552)

    #pragma omp simd aligned(g_x_0_z_0_yy_xy_0_xx, g_x_0_z_0_yy_xy_0_xy, g_x_0_z_0_yy_xy_0_xz, g_x_0_z_0_yy_xy_0_yy, g_x_0_z_0_yy_xy_0_yz, g_x_0_z_0_yy_xy_0_zz, g_xyy_xy_z_xx, g_xyy_xy_z_xy, g_xyy_xy_z_xz, g_xyy_xy_z_yy, g_xyy_xy_z_yz, g_xyy_xy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_yy_xy_0_xx[i] = 4.0 * g_xyy_xy_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_xy_0_xy[i] = 4.0 * g_xyy_xy_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_xy_0_xz[i] = 4.0 * g_xyy_xy_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_xy_0_yy[i] = 4.0 * g_xyy_xy_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_xy_0_yz[i] = 4.0 * g_xyy_xy_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_xy_0_zz[i] = 4.0 * g_xyy_xy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (552-558)

    #pragma omp simd aligned(g_x_0_z_0_yy_xz_0_xx, g_x_0_z_0_yy_xz_0_xy, g_x_0_z_0_yy_xz_0_xz, g_x_0_z_0_yy_xz_0_yy, g_x_0_z_0_yy_xz_0_yz, g_x_0_z_0_yy_xz_0_zz, g_xyy_xz_z_xx, g_xyy_xz_z_xy, g_xyy_xz_z_xz, g_xyy_xz_z_yy, g_xyy_xz_z_yz, g_xyy_xz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_yy_xz_0_xx[i] = 4.0 * g_xyy_xz_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_xz_0_xy[i] = 4.0 * g_xyy_xz_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_xz_0_xz[i] = 4.0 * g_xyy_xz_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_xz_0_yy[i] = 4.0 * g_xyy_xz_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_xz_0_yz[i] = 4.0 * g_xyy_xz_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_xz_0_zz[i] = 4.0 * g_xyy_xz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (558-564)

    #pragma omp simd aligned(g_x_0_z_0_yy_yy_0_xx, g_x_0_z_0_yy_yy_0_xy, g_x_0_z_0_yy_yy_0_xz, g_x_0_z_0_yy_yy_0_yy, g_x_0_z_0_yy_yy_0_yz, g_x_0_z_0_yy_yy_0_zz, g_xyy_yy_z_xx, g_xyy_yy_z_xy, g_xyy_yy_z_xz, g_xyy_yy_z_yy, g_xyy_yy_z_yz, g_xyy_yy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_yy_yy_0_xx[i] = 4.0 * g_xyy_yy_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_yy_0_xy[i] = 4.0 * g_xyy_yy_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_yy_0_xz[i] = 4.0 * g_xyy_yy_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_yy_0_yy[i] = 4.0 * g_xyy_yy_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_yy_0_yz[i] = 4.0 * g_xyy_yy_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_yy_0_zz[i] = 4.0 * g_xyy_yy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (564-570)

    #pragma omp simd aligned(g_x_0_z_0_yy_yz_0_xx, g_x_0_z_0_yy_yz_0_xy, g_x_0_z_0_yy_yz_0_xz, g_x_0_z_0_yy_yz_0_yy, g_x_0_z_0_yy_yz_0_yz, g_x_0_z_0_yy_yz_0_zz, g_xyy_yz_z_xx, g_xyy_yz_z_xy, g_xyy_yz_z_xz, g_xyy_yz_z_yy, g_xyy_yz_z_yz, g_xyy_yz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_yy_yz_0_xx[i] = 4.0 * g_xyy_yz_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_yz_0_xy[i] = 4.0 * g_xyy_yz_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_yz_0_xz[i] = 4.0 * g_xyy_yz_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_yz_0_yy[i] = 4.0 * g_xyy_yz_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_yz_0_yz[i] = 4.0 * g_xyy_yz_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_yz_0_zz[i] = 4.0 * g_xyy_yz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (570-576)

    #pragma omp simd aligned(g_x_0_z_0_yy_zz_0_xx, g_x_0_z_0_yy_zz_0_xy, g_x_0_z_0_yy_zz_0_xz, g_x_0_z_0_yy_zz_0_yy, g_x_0_z_0_yy_zz_0_yz, g_x_0_z_0_yy_zz_0_zz, g_xyy_zz_z_xx, g_xyy_zz_z_xy, g_xyy_zz_z_xz, g_xyy_zz_z_yy, g_xyy_zz_z_yz, g_xyy_zz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_yy_zz_0_xx[i] = 4.0 * g_xyy_zz_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_zz_0_xy[i] = 4.0 * g_xyy_zz_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_zz_0_xz[i] = 4.0 * g_xyy_zz_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_zz_0_yy[i] = 4.0 * g_xyy_zz_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_zz_0_yz[i] = 4.0 * g_xyy_zz_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_zz_0_zz[i] = 4.0 * g_xyy_zz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (576-582)

    #pragma omp simd aligned(g_x_0_z_0_yz_xx_0_xx, g_x_0_z_0_yz_xx_0_xy, g_x_0_z_0_yz_xx_0_xz, g_x_0_z_0_yz_xx_0_yy, g_x_0_z_0_yz_xx_0_yz, g_x_0_z_0_yz_xx_0_zz, g_xyz_xx_z_xx, g_xyz_xx_z_xy, g_xyz_xx_z_xz, g_xyz_xx_z_yy, g_xyz_xx_z_yz, g_xyz_xx_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_yz_xx_0_xx[i] = 4.0 * g_xyz_xx_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_xx_0_xy[i] = 4.0 * g_xyz_xx_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_xx_0_xz[i] = 4.0 * g_xyz_xx_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_xx_0_yy[i] = 4.0 * g_xyz_xx_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_xx_0_yz[i] = 4.0 * g_xyz_xx_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_xx_0_zz[i] = 4.0 * g_xyz_xx_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (582-588)

    #pragma omp simd aligned(g_x_0_z_0_yz_xy_0_xx, g_x_0_z_0_yz_xy_0_xy, g_x_0_z_0_yz_xy_0_xz, g_x_0_z_0_yz_xy_0_yy, g_x_0_z_0_yz_xy_0_yz, g_x_0_z_0_yz_xy_0_zz, g_xyz_xy_z_xx, g_xyz_xy_z_xy, g_xyz_xy_z_xz, g_xyz_xy_z_yy, g_xyz_xy_z_yz, g_xyz_xy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_yz_xy_0_xx[i] = 4.0 * g_xyz_xy_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_xy_0_xy[i] = 4.0 * g_xyz_xy_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_xy_0_xz[i] = 4.0 * g_xyz_xy_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_xy_0_yy[i] = 4.0 * g_xyz_xy_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_xy_0_yz[i] = 4.0 * g_xyz_xy_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_xy_0_zz[i] = 4.0 * g_xyz_xy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (588-594)

    #pragma omp simd aligned(g_x_0_z_0_yz_xz_0_xx, g_x_0_z_0_yz_xz_0_xy, g_x_0_z_0_yz_xz_0_xz, g_x_0_z_0_yz_xz_0_yy, g_x_0_z_0_yz_xz_0_yz, g_x_0_z_0_yz_xz_0_zz, g_xyz_xz_z_xx, g_xyz_xz_z_xy, g_xyz_xz_z_xz, g_xyz_xz_z_yy, g_xyz_xz_z_yz, g_xyz_xz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_yz_xz_0_xx[i] = 4.0 * g_xyz_xz_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_xz_0_xy[i] = 4.0 * g_xyz_xz_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_xz_0_xz[i] = 4.0 * g_xyz_xz_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_xz_0_yy[i] = 4.0 * g_xyz_xz_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_xz_0_yz[i] = 4.0 * g_xyz_xz_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_xz_0_zz[i] = 4.0 * g_xyz_xz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (594-600)

    #pragma omp simd aligned(g_x_0_z_0_yz_yy_0_xx, g_x_0_z_0_yz_yy_0_xy, g_x_0_z_0_yz_yy_0_xz, g_x_0_z_0_yz_yy_0_yy, g_x_0_z_0_yz_yy_0_yz, g_x_0_z_0_yz_yy_0_zz, g_xyz_yy_z_xx, g_xyz_yy_z_xy, g_xyz_yy_z_xz, g_xyz_yy_z_yy, g_xyz_yy_z_yz, g_xyz_yy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_yz_yy_0_xx[i] = 4.0 * g_xyz_yy_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_yy_0_xy[i] = 4.0 * g_xyz_yy_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_yy_0_xz[i] = 4.0 * g_xyz_yy_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_yy_0_yy[i] = 4.0 * g_xyz_yy_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_yy_0_yz[i] = 4.0 * g_xyz_yy_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_yy_0_zz[i] = 4.0 * g_xyz_yy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (600-606)

    #pragma omp simd aligned(g_x_0_z_0_yz_yz_0_xx, g_x_0_z_0_yz_yz_0_xy, g_x_0_z_0_yz_yz_0_xz, g_x_0_z_0_yz_yz_0_yy, g_x_0_z_0_yz_yz_0_yz, g_x_0_z_0_yz_yz_0_zz, g_xyz_yz_z_xx, g_xyz_yz_z_xy, g_xyz_yz_z_xz, g_xyz_yz_z_yy, g_xyz_yz_z_yz, g_xyz_yz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_yz_yz_0_xx[i] = 4.0 * g_xyz_yz_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_yz_0_xy[i] = 4.0 * g_xyz_yz_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_yz_0_xz[i] = 4.0 * g_xyz_yz_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_yz_0_yy[i] = 4.0 * g_xyz_yz_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_yz_0_yz[i] = 4.0 * g_xyz_yz_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_yz_0_zz[i] = 4.0 * g_xyz_yz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (606-612)

    #pragma omp simd aligned(g_x_0_z_0_yz_zz_0_xx, g_x_0_z_0_yz_zz_0_xy, g_x_0_z_0_yz_zz_0_xz, g_x_0_z_0_yz_zz_0_yy, g_x_0_z_0_yz_zz_0_yz, g_x_0_z_0_yz_zz_0_zz, g_xyz_zz_z_xx, g_xyz_zz_z_xy, g_xyz_zz_z_xz, g_xyz_zz_z_yy, g_xyz_zz_z_yz, g_xyz_zz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_yz_zz_0_xx[i] = 4.0 * g_xyz_zz_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_zz_0_xy[i] = 4.0 * g_xyz_zz_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_zz_0_xz[i] = 4.0 * g_xyz_zz_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_zz_0_yy[i] = 4.0 * g_xyz_zz_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_zz_0_yz[i] = 4.0 * g_xyz_zz_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_zz_0_zz[i] = 4.0 * g_xyz_zz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (612-618)

    #pragma omp simd aligned(g_x_0_z_0_zz_xx_0_xx, g_x_0_z_0_zz_xx_0_xy, g_x_0_z_0_zz_xx_0_xz, g_x_0_z_0_zz_xx_0_yy, g_x_0_z_0_zz_xx_0_yz, g_x_0_z_0_zz_xx_0_zz, g_xzz_xx_z_xx, g_xzz_xx_z_xy, g_xzz_xx_z_xz, g_xzz_xx_z_yy, g_xzz_xx_z_yz, g_xzz_xx_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_zz_xx_0_xx[i] = 4.0 * g_xzz_xx_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_xx_0_xy[i] = 4.0 * g_xzz_xx_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_xx_0_xz[i] = 4.0 * g_xzz_xx_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_xx_0_yy[i] = 4.0 * g_xzz_xx_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_xx_0_yz[i] = 4.0 * g_xzz_xx_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_xx_0_zz[i] = 4.0 * g_xzz_xx_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (618-624)

    #pragma omp simd aligned(g_x_0_z_0_zz_xy_0_xx, g_x_0_z_0_zz_xy_0_xy, g_x_0_z_0_zz_xy_0_xz, g_x_0_z_0_zz_xy_0_yy, g_x_0_z_0_zz_xy_0_yz, g_x_0_z_0_zz_xy_0_zz, g_xzz_xy_z_xx, g_xzz_xy_z_xy, g_xzz_xy_z_xz, g_xzz_xy_z_yy, g_xzz_xy_z_yz, g_xzz_xy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_zz_xy_0_xx[i] = 4.0 * g_xzz_xy_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_xy_0_xy[i] = 4.0 * g_xzz_xy_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_xy_0_xz[i] = 4.0 * g_xzz_xy_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_xy_0_yy[i] = 4.0 * g_xzz_xy_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_xy_0_yz[i] = 4.0 * g_xzz_xy_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_xy_0_zz[i] = 4.0 * g_xzz_xy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (624-630)

    #pragma omp simd aligned(g_x_0_z_0_zz_xz_0_xx, g_x_0_z_0_zz_xz_0_xy, g_x_0_z_0_zz_xz_0_xz, g_x_0_z_0_zz_xz_0_yy, g_x_0_z_0_zz_xz_0_yz, g_x_0_z_0_zz_xz_0_zz, g_xzz_xz_z_xx, g_xzz_xz_z_xy, g_xzz_xz_z_xz, g_xzz_xz_z_yy, g_xzz_xz_z_yz, g_xzz_xz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_zz_xz_0_xx[i] = 4.0 * g_xzz_xz_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_xz_0_xy[i] = 4.0 * g_xzz_xz_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_xz_0_xz[i] = 4.0 * g_xzz_xz_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_xz_0_yy[i] = 4.0 * g_xzz_xz_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_xz_0_yz[i] = 4.0 * g_xzz_xz_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_xz_0_zz[i] = 4.0 * g_xzz_xz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (630-636)

    #pragma omp simd aligned(g_x_0_z_0_zz_yy_0_xx, g_x_0_z_0_zz_yy_0_xy, g_x_0_z_0_zz_yy_0_xz, g_x_0_z_0_zz_yy_0_yy, g_x_0_z_0_zz_yy_0_yz, g_x_0_z_0_zz_yy_0_zz, g_xzz_yy_z_xx, g_xzz_yy_z_xy, g_xzz_yy_z_xz, g_xzz_yy_z_yy, g_xzz_yy_z_yz, g_xzz_yy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_zz_yy_0_xx[i] = 4.0 * g_xzz_yy_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_yy_0_xy[i] = 4.0 * g_xzz_yy_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_yy_0_xz[i] = 4.0 * g_xzz_yy_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_yy_0_yy[i] = 4.0 * g_xzz_yy_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_yy_0_yz[i] = 4.0 * g_xzz_yy_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_yy_0_zz[i] = 4.0 * g_xzz_yy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (636-642)

    #pragma omp simd aligned(g_x_0_z_0_zz_yz_0_xx, g_x_0_z_0_zz_yz_0_xy, g_x_0_z_0_zz_yz_0_xz, g_x_0_z_0_zz_yz_0_yy, g_x_0_z_0_zz_yz_0_yz, g_x_0_z_0_zz_yz_0_zz, g_xzz_yz_z_xx, g_xzz_yz_z_xy, g_xzz_yz_z_xz, g_xzz_yz_z_yy, g_xzz_yz_z_yz, g_xzz_yz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_zz_yz_0_xx[i] = 4.0 * g_xzz_yz_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_yz_0_xy[i] = 4.0 * g_xzz_yz_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_yz_0_xz[i] = 4.0 * g_xzz_yz_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_yz_0_yy[i] = 4.0 * g_xzz_yz_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_yz_0_yz[i] = 4.0 * g_xzz_yz_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_yz_0_zz[i] = 4.0 * g_xzz_yz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (642-648)

    #pragma omp simd aligned(g_x_0_z_0_zz_zz_0_xx, g_x_0_z_0_zz_zz_0_xy, g_x_0_z_0_zz_zz_0_xz, g_x_0_z_0_zz_zz_0_yy, g_x_0_z_0_zz_zz_0_yz, g_x_0_z_0_zz_zz_0_zz, g_xzz_zz_z_xx, g_xzz_zz_z_xy, g_xzz_zz_z_xz, g_xzz_zz_z_yy, g_xzz_zz_z_yz, g_xzz_zz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_zz_zz_0_xx[i] = 4.0 * g_xzz_zz_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_zz_0_xy[i] = 4.0 * g_xzz_zz_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_zz_0_xz[i] = 4.0 * g_xzz_zz_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_zz_0_yy[i] = 4.0 * g_xzz_zz_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_zz_0_yz[i] = 4.0 * g_xzz_zz_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_zz_0_zz[i] = 4.0 * g_xzz_zz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (648-654)

    #pragma omp simd aligned(g_xxy_xx_x_xx, g_xxy_xx_x_xy, g_xxy_xx_x_xz, g_xxy_xx_x_yy, g_xxy_xx_x_yz, g_xxy_xx_x_zz, g_y_0_x_0_xx_xx_0_xx, g_y_0_x_0_xx_xx_0_xy, g_y_0_x_0_xx_xx_0_xz, g_y_0_x_0_xx_xx_0_yy, g_y_0_x_0_xx_xx_0_yz, g_y_0_x_0_xx_xx_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_xx_xx_0_xx[i] = 4.0 * g_xxy_xx_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_xx_0_xy[i] = 4.0 * g_xxy_xx_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_xx_0_xz[i] = 4.0 * g_xxy_xx_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_xx_0_yy[i] = 4.0 * g_xxy_xx_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_xx_0_yz[i] = 4.0 * g_xxy_xx_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_xx_0_zz[i] = 4.0 * g_xxy_xx_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (654-660)

    #pragma omp simd aligned(g_xxy_xy_x_xx, g_xxy_xy_x_xy, g_xxy_xy_x_xz, g_xxy_xy_x_yy, g_xxy_xy_x_yz, g_xxy_xy_x_zz, g_y_0_x_0_xx_xy_0_xx, g_y_0_x_0_xx_xy_0_xy, g_y_0_x_0_xx_xy_0_xz, g_y_0_x_0_xx_xy_0_yy, g_y_0_x_0_xx_xy_0_yz, g_y_0_x_0_xx_xy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_xx_xy_0_xx[i] = 4.0 * g_xxy_xy_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_xy_0_xy[i] = 4.0 * g_xxy_xy_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_xy_0_xz[i] = 4.0 * g_xxy_xy_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_xy_0_yy[i] = 4.0 * g_xxy_xy_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_xy_0_yz[i] = 4.0 * g_xxy_xy_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_xy_0_zz[i] = 4.0 * g_xxy_xy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (660-666)

    #pragma omp simd aligned(g_xxy_xz_x_xx, g_xxy_xz_x_xy, g_xxy_xz_x_xz, g_xxy_xz_x_yy, g_xxy_xz_x_yz, g_xxy_xz_x_zz, g_y_0_x_0_xx_xz_0_xx, g_y_0_x_0_xx_xz_0_xy, g_y_0_x_0_xx_xz_0_xz, g_y_0_x_0_xx_xz_0_yy, g_y_0_x_0_xx_xz_0_yz, g_y_0_x_0_xx_xz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_xx_xz_0_xx[i] = 4.0 * g_xxy_xz_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_xz_0_xy[i] = 4.0 * g_xxy_xz_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_xz_0_xz[i] = 4.0 * g_xxy_xz_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_xz_0_yy[i] = 4.0 * g_xxy_xz_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_xz_0_yz[i] = 4.0 * g_xxy_xz_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_xz_0_zz[i] = 4.0 * g_xxy_xz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (666-672)

    #pragma omp simd aligned(g_xxy_yy_x_xx, g_xxy_yy_x_xy, g_xxy_yy_x_xz, g_xxy_yy_x_yy, g_xxy_yy_x_yz, g_xxy_yy_x_zz, g_y_0_x_0_xx_yy_0_xx, g_y_0_x_0_xx_yy_0_xy, g_y_0_x_0_xx_yy_0_xz, g_y_0_x_0_xx_yy_0_yy, g_y_0_x_0_xx_yy_0_yz, g_y_0_x_0_xx_yy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_xx_yy_0_xx[i] = 4.0 * g_xxy_yy_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_yy_0_xy[i] = 4.0 * g_xxy_yy_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_yy_0_xz[i] = 4.0 * g_xxy_yy_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_yy_0_yy[i] = 4.0 * g_xxy_yy_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_yy_0_yz[i] = 4.0 * g_xxy_yy_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_yy_0_zz[i] = 4.0 * g_xxy_yy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (672-678)

    #pragma omp simd aligned(g_xxy_yz_x_xx, g_xxy_yz_x_xy, g_xxy_yz_x_xz, g_xxy_yz_x_yy, g_xxy_yz_x_yz, g_xxy_yz_x_zz, g_y_0_x_0_xx_yz_0_xx, g_y_0_x_0_xx_yz_0_xy, g_y_0_x_0_xx_yz_0_xz, g_y_0_x_0_xx_yz_0_yy, g_y_0_x_0_xx_yz_0_yz, g_y_0_x_0_xx_yz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_xx_yz_0_xx[i] = 4.0 * g_xxy_yz_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_yz_0_xy[i] = 4.0 * g_xxy_yz_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_yz_0_xz[i] = 4.0 * g_xxy_yz_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_yz_0_yy[i] = 4.0 * g_xxy_yz_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_yz_0_yz[i] = 4.0 * g_xxy_yz_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_yz_0_zz[i] = 4.0 * g_xxy_yz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (678-684)

    #pragma omp simd aligned(g_xxy_zz_x_xx, g_xxy_zz_x_xy, g_xxy_zz_x_xz, g_xxy_zz_x_yy, g_xxy_zz_x_yz, g_xxy_zz_x_zz, g_y_0_x_0_xx_zz_0_xx, g_y_0_x_0_xx_zz_0_xy, g_y_0_x_0_xx_zz_0_xz, g_y_0_x_0_xx_zz_0_yy, g_y_0_x_0_xx_zz_0_yz, g_y_0_x_0_xx_zz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_xx_zz_0_xx[i] = 4.0 * g_xxy_zz_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_zz_0_xy[i] = 4.0 * g_xxy_zz_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_zz_0_xz[i] = 4.0 * g_xxy_zz_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_zz_0_yy[i] = 4.0 * g_xxy_zz_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_zz_0_yz[i] = 4.0 * g_xxy_zz_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_zz_0_zz[i] = 4.0 * g_xxy_zz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (684-690)

    #pragma omp simd aligned(g_x_xx_x_xx, g_x_xx_x_xy, g_x_xx_x_xz, g_x_xx_x_yy, g_x_xx_x_yz, g_x_xx_x_zz, g_xyy_xx_x_xx, g_xyy_xx_x_xy, g_xyy_xx_x_xz, g_xyy_xx_x_yy, g_xyy_xx_x_yz, g_xyy_xx_x_zz, g_y_0_x_0_xy_xx_0_xx, g_y_0_x_0_xy_xx_0_xy, g_y_0_x_0_xy_xx_0_xz, g_y_0_x_0_xy_xx_0_yy, g_y_0_x_0_xy_xx_0_yz, g_y_0_x_0_xy_xx_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_xy_xx_0_xx[i] = -2.0 * g_x_xx_x_xx[i] * c_exps[i] + 4.0 * g_xyy_xx_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_xx_0_xy[i] = -2.0 * g_x_xx_x_xy[i] * c_exps[i] + 4.0 * g_xyy_xx_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_xx_0_xz[i] = -2.0 * g_x_xx_x_xz[i] * c_exps[i] + 4.0 * g_xyy_xx_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_xx_0_yy[i] = -2.0 * g_x_xx_x_yy[i] * c_exps[i] + 4.0 * g_xyy_xx_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_xx_0_yz[i] = -2.0 * g_x_xx_x_yz[i] * c_exps[i] + 4.0 * g_xyy_xx_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_xx_0_zz[i] = -2.0 * g_x_xx_x_zz[i] * c_exps[i] + 4.0 * g_xyy_xx_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (690-696)

    #pragma omp simd aligned(g_x_xy_x_xx, g_x_xy_x_xy, g_x_xy_x_xz, g_x_xy_x_yy, g_x_xy_x_yz, g_x_xy_x_zz, g_xyy_xy_x_xx, g_xyy_xy_x_xy, g_xyy_xy_x_xz, g_xyy_xy_x_yy, g_xyy_xy_x_yz, g_xyy_xy_x_zz, g_y_0_x_0_xy_xy_0_xx, g_y_0_x_0_xy_xy_0_xy, g_y_0_x_0_xy_xy_0_xz, g_y_0_x_0_xy_xy_0_yy, g_y_0_x_0_xy_xy_0_yz, g_y_0_x_0_xy_xy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_xy_xy_0_xx[i] = -2.0 * g_x_xy_x_xx[i] * c_exps[i] + 4.0 * g_xyy_xy_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_xy_0_xy[i] = -2.0 * g_x_xy_x_xy[i] * c_exps[i] + 4.0 * g_xyy_xy_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_xy_0_xz[i] = -2.0 * g_x_xy_x_xz[i] * c_exps[i] + 4.0 * g_xyy_xy_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_xy_0_yy[i] = -2.0 * g_x_xy_x_yy[i] * c_exps[i] + 4.0 * g_xyy_xy_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_xy_0_yz[i] = -2.0 * g_x_xy_x_yz[i] * c_exps[i] + 4.0 * g_xyy_xy_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_xy_0_zz[i] = -2.0 * g_x_xy_x_zz[i] * c_exps[i] + 4.0 * g_xyy_xy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (696-702)

    #pragma omp simd aligned(g_x_xz_x_xx, g_x_xz_x_xy, g_x_xz_x_xz, g_x_xz_x_yy, g_x_xz_x_yz, g_x_xz_x_zz, g_xyy_xz_x_xx, g_xyy_xz_x_xy, g_xyy_xz_x_xz, g_xyy_xz_x_yy, g_xyy_xz_x_yz, g_xyy_xz_x_zz, g_y_0_x_0_xy_xz_0_xx, g_y_0_x_0_xy_xz_0_xy, g_y_0_x_0_xy_xz_0_xz, g_y_0_x_0_xy_xz_0_yy, g_y_0_x_0_xy_xz_0_yz, g_y_0_x_0_xy_xz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_xy_xz_0_xx[i] = -2.0 * g_x_xz_x_xx[i] * c_exps[i] + 4.0 * g_xyy_xz_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_xz_0_xy[i] = -2.0 * g_x_xz_x_xy[i] * c_exps[i] + 4.0 * g_xyy_xz_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_xz_0_xz[i] = -2.0 * g_x_xz_x_xz[i] * c_exps[i] + 4.0 * g_xyy_xz_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_xz_0_yy[i] = -2.0 * g_x_xz_x_yy[i] * c_exps[i] + 4.0 * g_xyy_xz_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_xz_0_yz[i] = -2.0 * g_x_xz_x_yz[i] * c_exps[i] + 4.0 * g_xyy_xz_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_xz_0_zz[i] = -2.0 * g_x_xz_x_zz[i] * c_exps[i] + 4.0 * g_xyy_xz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (702-708)

    #pragma omp simd aligned(g_x_yy_x_xx, g_x_yy_x_xy, g_x_yy_x_xz, g_x_yy_x_yy, g_x_yy_x_yz, g_x_yy_x_zz, g_xyy_yy_x_xx, g_xyy_yy_x_xy, g_xyy_yy_x_xz, g_xyy_yy_x_yy, g_xyy_yy_x_yz, g_xyy_yy_x_zz, g_y_0_x_0_xy_yy_0_xx, g_y_0_x_0_xy_yy_0_xy, g_y_0_x_0_xy_yy_0_xz, g_y_0_x_0_xy_yy_0_yy, g_y_0_x_0_xy_yy_0_yz, g_y_0_x_0_xy_yy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_xy_yy_0_xx[i] = -2.0 * g_x_yy_x_xx[i] * c_exps[i] + 4.0 * g_xyy_yy_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_yy_0_xy[i] = -2.0 * g_x_yy_x_xy[i] * c_exps[i] + 4.0 * g_xyy_yy_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_yy_0_xz[i] = -2.0 * g_x_yy_x_xz[i] * c_exps[i] + 4.0 * g_xyy_yy_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_yy_0_yy[i] = -2.0 * g_x_yy_x_yy[i] * c_exps[i] + 4.0 * g_xyy_yy_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_yy_0_yz[i] = -2.0 * g_x_yy_x_yz[i] * c_exps[i] + 4.0 * g_xyy_yy_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_yy_0_zz[i] = -2.0 * g_x_yy_x_zz[i] * c_exps[i] + 4.0 * g_xyy_yy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (708-714)

    #pragma omp simd aligned(g_x_yz_x_xx, g_x_yz_x_xy, g_x_yz_x_xz, g_x_yz_x_yy, g_x_yz_x_yz, g_x_yz_x_zz, g_xyy_yz_x_xx, g_xyy_yz_x_xy, g_xyy_yz_x_xz, g_xyy_yz_x_yy, g_xyy_yz_x_yz, g_xyy_yz_x_zz, g_y_0_x_0_xy_yz_0_xx, g_y_0_x_0_xy_yz_0_xy, g_y_0_x_0_xy_yz_0_xz, g_y_0_x_0_xy_yz_0_yy, g_y_0_x_0_xy_yz_0_yz, g_y_0_x_0_xy_yz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_xy_yz_0_xx[i] = -2.0 * g_x_yz_x_xx[i] * c_exps[i] + 4.0 * g_xyy_yz_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_yz_0_xy[i] = -2.0 * g_x_yz_x_xy[i] * c_exps[i] + 4.0 * g_xyy_yz_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_yz_0_xz[i] = -2.0 * g_x_yz_x_xz[i] * c_exps[i] + 4.0 * g_xyy_yz_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_yz_0_yy[i] = -2.0 * g_x_yz_x_yy[i] * c_exps[i] + 4.0 * g_xyy_yz_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_yz_0_yz[i] = -2.0 * g_x_yz_x_yz[i] * c_exps[i] + 4.0 * g_xyy_yz_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_yz_0_zz[i] = -2.0 * g_x_yz_x_zz[i] * c_exps[i] + 4.0 * g_xyy_yz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (714-720)

    #pragma omp simd aligned(g_x_zz_x_xx, g_x_zz_x_xy, g_x_zz_x_xz, g_x_zz_x_yy, g_x_zz_x_yz, g_x_zz_x_zz, g_xyy_zz_x_xx, g_xyy_zz_x_xy, g_xyy_zz_x_xz, g_xyy_zz_x_yy, g_xyy_zz_x_yz, g_xyy_zz_x_zz, g_y_0_x_0_xy_zz_0_xx, g_y_0_x_0_xy_zz_0_xy, g_y_0_x_0_xy_zz_0_xz, g_y_0_x_0_xy_zz_0_yy, g_y_0_x_0_xy_zz_0_yz, g_y_0_x_0_xy_zz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_xy_zz_0_xx[i] = -2.0 * g_x_zz_x_xx[i] * c_exps[i] + 4.0 * g_xyy_zz_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_zz_0_xy[i] = -2.0 * g_x_zz_x_xy[i] * c_exps[i] + 4.0 * g_xyy_zz_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_zz_0_xz[i] = -2.0 * g_x_zz_x_xz[i] * c_exps[i] + 4.0 * g_xyy_zz_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_zz_0_yy[i] = -2.0 * g_x_zz_x_yy[i] * c_exps[i] + 4.0 * g_xyy_zz_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_zz_0_yz[i] = -2.0 * g_x_zz_x_yz[i] * c_exps[i] + 4.0 * g_xyy_zz_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_zz_0_zz[i] = -2.0 * g_x_zz_x_zz[i] * c_exps[i] + 4.0 * g_xyy_zz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (720-726)

    #pragma omp simd aligned(g_xyz_xx_x_xx, g_xyz_xx_x_xy, g_xyz_xx_x_xz, g_xyz_xx_x_yy, g_xyz_xx_x_yz, g_xyz_xx_x_zz, g_y_0_x_0_xz_xx_0_xx, g_y_0_x_0_xz_xx_0_xy, g_y_0_x_0_xz_xx_0_xz, g_y_0_x_0_xz_xx_0_yy, g_y_0_x_0_xz_xx_0_yz, g_y_0_x_0_xz_xx_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_xz_xx_0_xx[i] = 4.0 * g_xyz_xx_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_xx_0_xy[i] = 4.0 * g_xyz_xx_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_xx_0_xz[i] = 4.0 * g_xyz_xx_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_xx_0_yy[i] = 4.0 * g_xyz_xx_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_xx_0_yz[i] = 4.0 * g_xyz_xx_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_xx_0_zz[i] = 4.0 * g_xyz_xx_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (726-732)

    #pragma omp simd aligned(g_xyz_xy_x_xx, g_xyz_xy_x_xy, g_xyz_xy_x_xz, g_xyz_xy_x_yy, g_xyz_xy_x_yz, g_xyz_xy_x_zz, g_y_0_x_0_xz_xy_0_xx, g_y_0_x_0_xz_xy_0_xy, g_y_0_x_0_xz_xy_0_xz, g_y_0_x_0_xz_xy_0_yy, g_y_0_x_0_xz_xy_0_yz, g_y_0_x_0_xz_xy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_xz_xy_0_xx[i] = 4.0 * g_xyz_xy_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_xy_0_xy[i] = 4.0 * g_xyz_xy_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_xy_0_xz[i] = 4.0 * g_xyz_xy_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_xy_0_yy[i] = 4.0 * g_xyz_xy_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_xy_0_yz[i] = 4.0 * g_xyz_xy_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_xy_0_zz[i] = 4.0 * g_xyz_xy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (732-738)

    #pragma omp simd aligned(g_xyz_xz_x_xx, g_xyz_xz_x_xy, g_xyz_xz_x_xz, g_xyz_xz_x_yy, g_xyz_xz_x_yz, g_xyz_xz_x_zz, g_y_0_x_0_xz_xz_0_xx, g_y_0_x_0_xz_xz_0_xy, g_y_0_x_0_xz_xz_0_xz, g_y_0_x_0_xz_xz_0_yy, g_y_0_x_0_xz_xz_0_yz, g_y_0_x_0_xz_xz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_xz_xz_0_xx[i] = 4.0 * g_xyz_xz_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_xz_0_xy[i] = 4.0 * g_xyz_xz_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_xz_0_xz[i] = 4.0 * g_xyz_xz_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_xz_0_yy[i] = 4.0 * g_xyz_xz_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_xz_0_yz[i] = 4.0 * g_xyz_xz_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_xz_0_zz[i] = 4.0 * g_xyz_xz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (738-744)

    #pragma omp simd aligned(g_xyz_yy_x_xx, g_xyz_yy_x_xy, g_xyz_yy_x_xz, g_xyz_yy_x_yy, g_xyz_yy_x_yz, g_xyz_yy_x_zz, g_y_0_x_0_xz_yy_0_xx, g_y_0_x_0_xz_yy_0_xy, g_y_0_x_0_xz_yy_0_xz, g_y_0_x_0_xz_yy_0_yy, g_y_0_x_0_xz_yy_0_yz, g_y_0_x_0_xz_yy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_xz_yy_0_xx[i] = 4.0 * g_xyz_yy_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_yy_0_xy[i] = 4.0 * g_xyz_yy_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_yy_0_xz[i] = 4.0 * g_xyz_yy_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_yy_0_yy[i] = 4.0 * g_xyz_yy_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_yy_0_yz[i] = 4.0 * g_xyz_yy_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_yy_0_zz[i] = 4.0 * g_xyz_yy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (744-750)

    #pragma omp simd aligned(g_xyz_yz_x_xx, g_xyz_yz_x_xy, g_xyz_yz_x_xz, g_xyz_yz_x_yy, g_xyz_yz_x_yz, g_xyz_yz_x_zz, g_y_0_x_0_xz_yz_0_xx, g_y_0_x_0_xz_yz_0_xy, g_y_0_x_0_xz_yz_0_xz, g_y_0_x_0_xz_yz_0_yy, g_y_0_x_0_xz_yz_0_yz, g_y_0_x_0_xz_yz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_xz_yz_0_xx[i] = 4.0 * g_xyz_yz_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_yz_0_xy[i] = 4.0 * g_xyz_yz_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_yz_0_xz[i] = 4.0 * g_xyz_yz_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_yz_0_yy[i] = 4.0 * g_xyz_yz_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_yz_0_yz[i] = 4.0 * g_xyz_yz_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_yz_0_zz[i] = 4.0 * g_xyz_yz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (750-756)

    #pragma omp simd aligned(g_xyz_zz_x_xx, g_xyz_zz_x_xy, g_xyz_zz_x_xz, g_xyz_zz_x_yy, g_xyz_zz_x_yz, g_xyz_zz_x_zz, g_y_0_x_0_xz_zz_0_xx, g_y_0_x_0_xz_zz_0_xy, g_y_0_x_0_xz_zz_0_xz, g_y_0_x_0_xz_zz_0_yy, g_y_0_x_0_xz_zz_0_yz, g_y_0_x_0_xz_zz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_xz_zz_0_xx[i] = 4.0 * g_xyz_zz_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_zz_0_xy[i] = 4.0 * g_xyz_zz_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_zz_0_xz[i] = 4.0 * g_xyz_zz_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_zz_0_yy[i] = 4.0 * g_xyz_zz_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_zz_0_yz[i] = 4.0 * g_xyz_zz_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_zz_0_zz[i] = 4.0 * g_xyz_zz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (756-762)

    #pragma omp simd aligned(g_y_0_x_0_yy_xx_0_xx, g_y_0_x_0_yy_xx_0_xy, g_y_0_x_0_yy_xx_0_xz, g_y_0_x_0_yy_xx_0_yy, g_y_0_x_0_yy_xx_0_yz, g_y_0_x_0_yy_xx_0_zz, g_y_xx_x_xx, g_y_xx_x_xy, g_y_xx_x_xz, g_y_xx_x_yy, g_y_xx_x_yz, g_y_xx_x_zz, g_yyy_xx_x_xx, g_yyy_xx_x_xy, g_yyy_xx_x_xz, g_yyy_xx_x_yy, g_yyy_xx_x_yz, g_yyy_xx_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_yy_xx_0_xx[i] = -4.0 * g_y_xx_x_xx[i] * c_exps[i] + 4.0 * g_yyy_xx_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_xx_0_xy[i] = -4.0 * g_y_xx_x_xy[i] * c_exps[i] + 4.0 * g_yyy_xx_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_xx_0_xz[i] = -4.0 * g_y_xx_x_xz[i] * c_exps[i] + 4.0 * g_yyy_xx_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_xx_0_yy[i] = -4.0 * g_y_xx_x_yy[i] * c_exps[i] + 4.0 * g_yyy_xx_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_xx_0_yz[i] = -4.0 * g_y_xx_x_yz[i] * c_exps[i] + 4.0 * g_yyy_xx_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_xx_0_zz[i] = -4.0 * g_y_xx_x_zz[i] * c_exps[i] + 4.0 * g_yyy_xx_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (762-768)

    #pragma omp simd aligned(g_y_0_x_0_yy_xy_0_xx, g_y_0_x_0_yy_xy_0_xy, g_y_0_x_0_yy_xy_0_xz, g_y_0_x_0_yy_xy_0_yy, g_y_0_x_0_yy_xy_0_yz, g_y_0_x_0_yy_xy_0_zz, g_y_xy_x_xx, g_y_xy_x_xy, g_y_xy_x_xz, g_y_xy_x_yy, g_y_xy_x_yz, g_y_xy_x_zz, g_yyy_xy_x_xx, g_yyy_xy_x_xy, g_yyy_xy_x_xz, g_yyy_xy_x_yy, g_yyy_xy_x_yz, g_yyy_xy_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_yy_xy_0_xx[i] = -4.0 * g_y_xy_x_xx[i] * c_exps[i] + 4.0 * g_yyy_xy_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_xy_0_xy[i] = -4.0 * g_y_xy_x_xy[i] * c_exps[i] + 4.0 * g_yyy_xy_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_xy_0_xz[i] = -4.0 * g_y_xy_x_xz[i] * c_exps[i] + 4.0 * g_yyy_xy_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_xy_0_yy[i] = -4.0 * g_y_xy_x_yy[i] * c_exps[i] + 4.0 * g_yyy_xy_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_xy_0_yz[i] = -4.0 * g_y_xy_x_yz[i] * c_exps[i] + 4.0 * g_yyy_xy_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_xy_0_zz[i] = -4.0 * g_y_xy_x_zz[i] * c_exps[i] + 4.0 * g_yyy_xy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (768-774)

    #pragma omp simd aligned(g_y_0_x_0_yy_xz_0_xx, g_y_0_x_0_yy_xz_0_xy, g_y_0_x_0_yy_xz_0_xz, g_y_0_x_0_yy_xz_0_yy, g_y_0_x_0_yy_xz_0_yz, g_y_0_x_0_yy_xz_0_zz, g_y_xz_x_xx, g_y_xz_x_xy, g_y_xz_x_xz, g_y_xz_x_yy, g_y_xz_x_yz, g_y_xz_x_zz, g_yyy_xz_x_xx, g_yyy_xz_x_xy, g_yyy_xz_x_xz, g_yyy_xz_x_yy, g_yyy_xz_x_yz, g_yyy_xz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_yy_xz_0_xx[i] = -4.0 * g_y_xz_x_xx[i] * c_exps[i] + 4.0 * g_yyy_xz_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_xz_0_xy[i] = -4.0 * g_y_xz_x_xy[i] * c_exps[i] + 4.0 * g_yyy_xz_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_xz_0_xz[i] = -4.0 * g_y_xz_x_xz[i] * c_exps[i] + 4.0 * g_yyy_xz_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_xz_0_yy[i] = -4.0 * g_y_xz_x_yy[i] * c_exps[i] + 4.0 * g_yyy_xz_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_xz_0_yz[i] = -4.0 * g_y_xz_x_yz[i] * c_exps[i] + 4.0 * g_yyy_xz_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_xz_0_zz[i] = -4.0 * g_y_xz_x_zz[i] * c_exps[i] + 4.0 * g_yyy_xz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (774-780)

    #pragma omp simd aligned(g_y_0_x_0_yy_yy_0_xx, g_y_0_x_0_yy_yy_0_xy, g_y_0_x_0_yy_yy_0_xz, g_y_0_x_0_yy_yy_0_yy, g_y_0_x_0_yy_yy_0_yz, g_y_0_x_0_yy_yy_0_zz, g_y_yy_x_xx, g_y_yy_x_xy, g_y_yy_x_xz, g_y_yy_x_yy, g_y_yy_x_yz, g_y_yy_x_zz, g_yyy_yy_x_xx, g_yyy_yy_x_xy, g_yyy_yy_x_xz, g_yyy_yy_x_yy, g_yyy_yy_x_yz, g_yyy_yy_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_yy_yy_0_xx[i] = -4.0 * g_y_yy_x_xx[i] * c_exps[i] + 4.0 * g_yyy_yy_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_yy_0_xy[i] = -4.0 * g_y_yy_x_xy[i] * c_exps[i] + 4.0 * g_yyy_yy_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_yy_0_xz[i] = -4.0 * g_y_yy_x_xz[i] * c_exps[i] + 4.0 * g_yyy_yy_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_yy_0_yy[i] = -4.0 * g_y_yy_x_yy[i] * c_exps[i] + 4.0 * g_yyy_yy_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_yy_0_yz[i] = -4.0 * g_y_yy_x_yz[i] * c_exps[i] + 4.0 * g_yyy_yy_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_yy_0_zz[i] = -4.0 * g_y_yy_x_zz[i] * c_exps[i] + 4.0 * g_yyy_yy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (780-786)

    #pragma omp simd aligned(g_y_0_x_0_yy_yz_0_xx, g_y_0_x_0_yy_yz_0_xy, g_y_0_x_0_yy_yz_0_xz, g_y_0_x_0_yy_yz_0_yy, g_y_0_x_0_yy_yz_0_yz, g_y_0_x_0_yy_yz_0_zz, g_y_yz_x_xx, g_y_yz_x_xy, g_y_yz_x_xz, g_y_yz_x_yy, g_y_yz_x_yz, g_y_yz_x_zz, g_yyy_yz_x_xx, g_yyy_yz_x_xy, g_yyy_yz_x_xz, g_yyy_yz_x_yy, g_yyy_yz_x_yz, g_yyy_yz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_yy_yz_0_xx[i] = -4.0 * g_y_yz_x_xx[i] * c_exps[i] + 4.0 * g_yyy_yz_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_yz_0_xy[i] = -4.0 * g_y_yz_x_xy[i] * c_exps[i] + 4.0 * g_yyy_yz_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_yz_0_xz[i] = -4.0 * g_y_yz_x_xz[i] * c_exps[i] + 4.0 * g_yyy_yz_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_yz_0_yy[i] = -4.0 * g_y_yz_x_yy[i] * c_exps[i] + 4.0 * g_yyy_yz_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_yz_0_yz[i] = -4.0 * g_y_yz_x_yz[i] * c_exps[i] + 4.0 * g_yyy_yz_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_yz_0_zz[i] = -4.0 * g_y_yz_x_zz[i] * c_exps[i] + 4.0 * g_yyy_yz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (786-792)

    #pragma omp simd aligned(g_y_0_x_0_yy_zz_0_xx, g_y_0_x_0_yy_zz_0_xy, g_y_0_x_0_yy_zz_0_xz, g_y_0_x_0_yy_zz_0_yy, g_y_0_x_0_yy_zz_0_yz, g_y_0_x_0_yy_zz_0_zz, g_y_zz_x_xx, g_y_zz_x_xy, g_y_zz_x_xz, g_y_zz_x_yy, g_y_zz_x_yz, g_y_zz_x_zz, g_yyy_zz_x_xx, g_yyy_zz_x_xy, g_yyy_zz_x_xz, g_yyy_zz_x_yy, g_yyy_zz_x_yz, g_yyy_zz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_yy_zz_0_xx[i] = -4.0 * g_y_zz_x_xx[i] * c_exps[i] + 4.0 * g_yyy_zz_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_zz_0_xy[i] = -4.0 * g_y_zz_x_xy[i] * c_exps[i] + 4.0 * g_yyy_zz_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_zz_0_xz[i] = -4.0 * g_y_zz_x_xz[i] * c_exps[i] + 4.0 * g_yyy_zz_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_zz_0_yy[i] = -4.0 * g_y_zz_x_yy[i] * c_exps[i] + 4.0 * g_yyy_zz_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_zz_0_yz[i] = -4.0 * g_y_zz_x_yz[i] * c_exps[i] + 4.0 * g_yyy_zz_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_zz_0_zz[i] = -4.0 * g_y_zz_x_zz[i] * c_exps[i] + 4.0 * g_yyy_zz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (792-798)

    #pragma omp simd aligned(g_y_0_x_0_yz_xx_0_xx, g_y_0_x_0_yz_xx_0_xy, g_y_0_x_0_yz_xx_0_xz, g_y_0_x_0_yz_xx_0_yy, g_y_0_x_0_yz_xx_0_yz, g_y_0_x_0_yz_xx_0_zz, g_yyz_xx_x_xx, g_yyz_xx_x_xy, g_yyz_xx_x_xz, g_yyz_xx_x_yy, g_yyz_xx_x_yz, g_yyz_xx_x_zz, g_z_xx_x_xx, g_z_xx_x_xy, g_z_xx_x_xz, g_z_xx_x_yy, g_z_xx_x_yz, g_z_xx_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_yz_xx_0_xx[i] = -2.0 * g_z_xx_x_xx[i] * c_exps[i] + 4.0 * g_yyz_xx_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_xx_0_xy[i] = -2.0 * g_z_xx_x_xy[i] * c_exps[i] + 4.0 * g_yyz_xx_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_xx_0_xz[i] = -2.0 * g_z_xx_x_xz[i] * c_exps[i] + 4.0 * g_yyz_xx_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_xx_0_yy[i] = -2.0 * g_z_xx_x_yy[i] * c_exps[i] + 4.0 * g_yyz_xx_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_xx_0_yz[i] = -2.0 * g_z_xx_x_yz[i] * c_exps[i] + 4.0 * g_yyz_xx_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_xx_0_zz[i] = -2.0 * g_z_xx_x_zz[i] * c_exps[i] + 4.0 * g_yyz_xx_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (798-804)

    #pragma omp simd aligned(g_y_0_x_0_yz_xy_0_xx, g_y_0_x_0_yz_xy_0_xy, g_y_0_x_0_yz_xy_0_xz, g_y_0_x_0_yz_xy_0_yy, g_y_0_x_0_yz_xy_0_yz, g_y_0_x_0_yz_xy_0_zz, g_yyz_xy_x_xx, g_yyz_xy_x_xy, g_yyz_xy_x_xz, g_yyz_xy_x_yy, g_yyz_xy_x_yz, g_yyz_xy_x_zz, g_z_xy_x_xx, g_z_xy_x_xy, g_z_xy_x_xz, g_z_xy_x_yy, g_z_xy_x_yz, g_z_xy_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_yz_xy_0_xx[i] = -2.0 * g_z_xy_x_xx[i] * c_exps[i] + 4.0 * g_yyz_xy_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_xy_0_xy[i] = -2.0 * g_z_xy_x_xy[i] * c_exps[i] + 4.0 * g_yyz_xy_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_xy_0_xz[i] = -2.0 * g_z_xy_x_xz[i] * c_exps[i] + 4.0 * g_yyz_xy_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_xy_0_yy[i] = -2.0 * g_z_xy_x_yy[i] * c_exps[i] + 4.0 * g_yyz_xy_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_xy_0_yz[i] = -2.0 * g_z_xy_x_yz[i] * c_exps[i] + 4.0 * g_yyz_xy_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_xy_0_zz[i] = -2.0 * g_z_xy_x_zz[i] * c_exps[i] + 4.0 * g_yyz_xy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (804-810)

    #pragma omp simd aligned(g_y_0_x_0_yz_xz_0_xx, g_y_0_x_0_yz_xz_0_xy, g_y_0_x_0_yz_xz_0_xz, g_y_0_x_0_yz_xz_0_yy, g_y_0_x_0_yz_xz_0_yz, g_y_0_x_0_yz_xz_0_zz, g_yyz_xz_x_xx, g_yyz_xz_x_xy, g_yyz_xz_x_xz, g_yyz_xz_x_yy, g_yyz_xz_x_yz, g_yyz_xz_x_zz, g_z_xz_x_xx, g_z_xz_x_xy, g_z_xz_x_xz, g_z_xz_x_yy, g_z_xz_x_yz, g_z_xz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_yz_xz_0_xx[i] = -2.0 * g_z_xz_x_xx[i] * c_exps[i] + 4.0 * g_yyz_xz_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_xz_0_xy[i] = -2.0 * g_z_xz_x_xy[i] * c_exps[i] + 4.0 * g_yyz_xz_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_xz_0_xz[i] = -2.0 * g_z_xz_x_xz[i] * c_exps[i] + 4.0 * g_yyz_xz_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_xz_0_yy[i] = -2.0 * g_z_xz_x_yy[i] * c_exps[i] + 4.0 * g_yyz_xz_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_xz_0_yz[i] = -2.0 * g_z_xz_x_yz[i] * c_exps[i] + 4.0 * g_yyz_xz_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_xz_0_zz[i] = -2.0 * g_z_xz_x_zz[i] * c_exps[i] + 4.0 * g_yyz_xz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (810-816)

    #pragma omp simd aligned(g_y_0_x_0_yz_yy_0_xx, g_y_0_x_0_yz_yy_0_xy, g_y_0_x_0_yz_yy_0_xz, g_y_0_x_0_yz_yy_0_yy, g_y_0_x_0_yz_yy_0_yz, g_y_0_x_0_yz_yy_0_zz, g_yyz_yy_x_xx, g_yyz_yy_x_xy, g_yyz_yy_x_xz, g_yyz_yy_x_yy, g_yyz_yy_x_yz, g_yyz_yy_x_zz, g_z_yy_x_xx, g_z_yy_x_xy, g_z_yy_x_xz, g_z_yy_x_yy, g_z_yy_x_yz, g_z_yy_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_yz_yy_0_xx[i] = -2.0 * g_z_yy_x_xx[i] * c_exps[i] + 4.0 * g_yyz_yy_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_yy_0_xy[i] = -2.0 * g_z_yy_x_xy[i] * c_exps[i] + 4.0 * g_yyz_yy_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_yy_0_xz[i] = -2.0 * g_z_yy_x_xz[i] * c_exps[i] + 4.0 * g_yyz_yy_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_yy_0_yy[i] = -2.0 * g_z_yy_x_yy[i] * c_exps[i] + 4.0 * g_yyz_yy_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_yy_0_yz[i] = -2.0 * g_z_yy_x_yz[i] * c_exps[i] + 4.0 * g_yyz_yy_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_yy_0_zz[i] = -2.0 * g_z_yy_x_zz[i] * c_exps[i] + 4.0 * g_yyz_yy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (816-822)

    #pragma omp simd aligned(g_y_0_x_0_yz_yz_0_xx, g_y_0_x_0_yz_yz_0_xy, g_y_0_x_0_yz_yz_0_xz, g_y_0_x_0_yz_yz_0_yy, g_y_0_x_0_yz_yz_0_yz, g_y_0_x_0_yz_yz_0_zz, g_yyz_yz_x_xx, g_yyz_yz_x_xy, g_yyz_yz_x_xz, g_yyz_yz_x_yy, g_yyz_yz_x_yz, g_yyz_yz_x_zz, g_z_yz_x_xx, g_z_yz_x_xy, g_z_yz_x_xz, g_z_yz_x_yy, g_z_yz_x_yz, g_z_yz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_yz_yz_0_xx[i] = -2.0 * g_z_yz_x_xx[i] * c_exps[i] + 4.0 * g_yyz_yz_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_yz_0_xy[i] = -2.0 * g_z_yz_x_xy[i] * c_exps[i] + 4.0 * g_yyz_yz_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_yz_0_xz[i] = -2.0 * g_z_yz_x_xz[i] * c_exps[i] + 4.0 * g_yyz_yz_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_yz_0_yy[i] = -2.0 * g_z_yz_x_yy[i] * c_exps[i] + 4.0 * g_yyz_yz_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_yz_0_yz[i] = -2.0 * g_z_yz_x_yz[i] * c_exps[i] + 4.0 * g_yyz_yz_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_yz_0_zz[i] = -2.0 * g_z_yz_x_zz[i] * c_exps[i] + 4.0 * g_yyz_yz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (822-828)

    #pragma omp simd aligned(g_y_0_x_0_yz_zz_0_xx, g_y_0_x_0_yz_zz_0_xy, g_y_0_x_0_yz_zz_0_xz, g_y_0_x_0_yz_zz_0_yy, g_y_0_x_0_yz_zz_0_yz, g_y_0_x_0_yz_zz_0_zz, g_yyz_zz_x_xx, g_yyz_zz_x_xy, g_yyz_zz_x_xz, g_yyz_zz_x_yy, g_yyz_zz_x_yz, g_yyz_zz_x_zz, g_z_zz_x_xx, g_z_zz_x_xy, g_z_zz_x_xz, g_z_zz_x_yy, g_z_zz_x_yz, g_z_zz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_yz_zz_0_xx[i] = -2.0 * g_z_zz_x_xx[i] * c_exps[i] + 4.0 * g_yyz_zz_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_zz_0_xy[i] = -2.0 * g_z_zz_x_xy[i] * c_exps[i] + 4.0 * g_yyz_zz_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_zz_0_xz[i] = -2.0 * g_z_zz_x_xz[i] * c_exps[i] + 4.0 * g_yyz_zz_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_zz_0_yy[i] = -2.0 * g_z_zz_x_yy[i] * c_exps[i] + 4.0 * g_yyz_zz_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_zz_0_yz[i] = -2.0 * g_z_zz_x_yz[i] * c_exps[i] + 4.0 * g_yyz_zz_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_zz_0_zz[i] = -2.0 * g_z_zz_x_zz[i] * c_exps[i] + 4.0 * g_yyz_zz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (828-834)

    #pragma omp simd aligned(g_y_0_x_0_zz_xx_0_xx, g_y_0_x_0_zz_xx_0_xy, g_y_0_x_0_zz_xx_0_xz, g_y_0_x_0_zz_xx_0_yy, g_y_0_x_0_zz_xx_0_yz, g_y_0_x_0_zz_xx_0_zz, g_yzz_xx_x_xx, g_yzz_xx_x_xy, g_yzz_xx_x_xz, g_yzz_xx_x_yy, g_yzz_xx_x_yz, g_yzz_xx_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_zz_xx_0_xx[i] = 4.0 * g_yzz_xx_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_xx_0_xy[i] = 4.0 * g_yzz_xx_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_xx_0_xz[i] = 4.0 * g_yzz_xx_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_xx_0_yy[i] = 4.0 * g_yzz_xx_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_xx_0_yz[i] = 4.0 * g_yzz_xx_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_xx_0_zz[i] = 4.0 * g_yzz_xx_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (834-840)

    #pragma omp simd aligned(g_y_0_x_0_zz_xy_0_xx, g_y_0_x_0_zz_xy_0_xy, g_y_0_x_0_zz_xy_0_xz, g_y_0_x_0_zz_xy_0_yy, g_y_0_x_0_zz_xy_0_yz, g_y_0_x_0_zz_xy_0_zz, g_yzz_xy_x_xx, g_yzz_xy_x_xy, g_yzz_xy_x_xz, g_yzz_xy_x_yy, g_yzz_xy_x_yz, g_yzz_xy_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_zz_xy_0_xx[i] = 4.0 * g_yzz_xy_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_xy_0_xy[i] = 4.0 * g_yzz_xy_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_xy_0_xz[i] = 4.0 * g_yzz_xy_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_xy_0_yy[i] = 4.0 * g_yzz_xy_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_xy_0_yz[i] = 4.0 * g_yzz_xy_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_xy_0_zz[i] = 4.0 * g_yzz_xy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (840-846)

    #pragma omp simd aligned(g_y_0_x_0_zz_xz_0_xx, g_y_0_x_0_zz_xz_0_xy, g_y_0_x_0_zz_xz_0_xz, g_y_0_x_0_zz_xz_0_yy, g_y_0_x_0_zz_xz_0_yz, g_y_0_x_0_zz_xz_0_zz, g_yzz_xz_x_xx, g_yzz_xz_x_xy, g_yzz_xz_x_xz, g_yzz_xz_x_yy, g_yzz_xz_x_yz, g_yzz_xz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_zz_xz_0_xx[i] = 4.0 * g_yzz_xz_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_xz_0_xy[i] = 4.0 * g_yzz_xz_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_xz_0_xz[i] = 4.0 * g_yzz_xz_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_xz_0_yy[i] = 4.0 * g_yzz_xz_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_xz_0_yz[i] = 4.0 * g_yzz_xz_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_xz_0_zz[i] = 4.0 * g_yzz_xz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (846-852)

    #pragma omp simd aligned(g_y_0_x_0_zz_yy_0_xx, g_y_0_x_0_zz_yy_0_xy, g_y_0_x_0_zz_yy_0_xz, g_y_0_x_0_zz_yy_0_yy, g_y_0_x_0_zz_yy_0_yz, g_y_0_x_0_zz_yy_0_zz, g_yzz_yy_x_xx, g_yzz_yy_x_xy, g_yzz_yy_x_xz, g_yzz_yy_x_yy, g_yzz_yy_x_yz, g_yzz_yy_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_zz_yy_0_xx[i] = 4.0 * g_yzz_yy_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_yy_0_xy[i] = 4.0 * g_yzz_yy_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_yy_0_xz[i] = 4.0 * g_yzz_yy_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_yy_0_yy[i] = 4.0 * g_yzz_yy_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_yy_0_yz[i] = 4.0 * g_yzz_yy_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_yy_0_zz[i] = 4.0 * g_yzz_yy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (852-858)

    #pragma omp simd aligned(g_y_0_x_0_zz_yz_0_xx, g_y_0_x_0_zz_yz_0_xy, g_y_0_x_0_zz_yz_0_xz, g_y_0_x_0_zz_yz_0_yy, g_y_0_x_0_zz_yz_0_yz, g_y_0_x_0_zz_yz_0_zz, g_yzz_yz_x_xx, g_yzz_yz_x_xy, g_yzz_yz_x_xz, g_yzz_yz_x_yy, g_yzz_yz_x_yz, g_yzz_yz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_zz_yz_0_xx[i] = 4.0 * g_yzz_yz_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_yz_0_xy[i] = 4.0 * g_yzz_yz_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_yz_0_xz[i] = 4.0 * g_yzz_yz_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_yz_0_yy[i] = 4.0 * g_yzz_yz_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_yz_0_yz[i] = 4.0 * g_yzz_yz_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_yz_0_zz[i] = 4.0 * g_yzz_yz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (858-864)

    #pragma omp simd aligned(g_y_0_x_0_zz_zz_0_xx, g_y_0_x_0_zz_zz_0_xy, g_y_0_x_0_zz_zz_0_xz, g_y_0_x_0_zz_zz_0_yy, g_y_0_x_0_zz_zz_0_yz, g_y_0_x_0_zz_zz_0_zz, g_yzz_zz_x_xx, g_yzz_zz_x_xy, g_yzz_zz_x_xz, g_yzz_zz_x_yy, g_yzz_zz_x_yz, g_yzz_zz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_zz_zz_0_xx[i] = 4.0 * g_yzz_zz_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_zz_0_xy[i] = 4.0 * g_yzz_zz_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_zz_0_xz[i] = 4.0 * g_yzz_zz_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_zz_0_yy[i] = 4.0 * g_yzz_zz_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_zz_0_yz[i] = 4.0 * g_yzz_zz_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_zz_0_zz[i] = 4.0 * g_yzz_zz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (864-870)

    #pragma omp simd aligned(g_xxy_xx_y_xx, g_xxy_xx_y_xy, g_xxy_xx_y_xz, g_xxy_xx_y_yy, g_xxy_xx_y_yz, g_xxy_xx_y_zz, g_y_0_y_0_xx_xx_0_xx, g_y_0_y_0_xx_xx_0_xy, g_y_0_y_0_xx_xx_0_xz, g_y_0_y_0_xx_xx_0_yy, g_y_0_y_0_xx_xx_0_yz, g_y_0_y_0_xx_xx_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_xx_xx_0_xx[i] = 4.0 * g_xxy_xx_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_xx_0_xy[i] = 4.0 * g_xxy_xx_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_xx_0_xz[i] = 4.0 * g_xxy_xx_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_xx_0_yy[i] = 4.0 * g_xxy_xx_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_xx_0_yz[i] = 4.0 * g_xxy_xx_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_xx_0_zz[i] = 4.0 * g_xxy_xx_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (870-876)

    #pragma omp simd aligned(g_xxy_xy_y_xx, g_xxy_xy_y_xy, g_xxy_xy_y_xz, g_xxy_xy_y_yy, g_xxy_xy_y_yz, g_xxy_xy_y_zz, g_y_0_y_0_xx_xy_0_xx, g_y_0_y_0_xx_xy_0_xy, g_y_0_y_0_xx_xy_0_xz, g_y_0_y_0_xx_xy_0_yy, g_y_0_y_0_xx_xy_0_yz, g_y_0_y_0_xx_xy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_xx_xy_0_xx[i] = 4.0 * g_xxy_xy_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_xy_0_xy[i] = 4.0 * g_xxy_xy_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_xy_0_xz[i] = 4.0 * g_xxy_xy_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_xy_0_yy[i] = 4.0 * g_xxy_xy_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_xy_0_yz[i] = 4.0 * g_xxy_xy_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_xy_0_zz[i] = 4.0 * g_xxy_xy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (876-882)

    #pragma omp simd aligned(g_xxy_xz_y_xx, g_xxy_xz_y_xy, g_xxy_xz_y_xz, g_xxy_xz_y_yy, g_xxy_xz_y_yz, g_xxy_xz_y_zz, g_y_0_y_0_xx_xz_0_xx, g_y_0_y_0_xx_xz_0_xy, g_y_0_y_0_xx_xz_0_xz, g_y_0_y_0_xx_xz_0_yy, g_y_0_y_0_xx_xz_0_yz, g_y_0_y_0_xx_xz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_xx_xz_0_xx[i] = 4.0 * g_xxy_xz_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_xz_0_xy[i] = 4.0 * g_xxy_xz_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_xz_0_xz[i] = 4.0 * g_xxy_xz_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_xz_0_yy[i] = 4.0 * g_xxy_xz_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_xz_0_yz[i] = 4.0 * g_xxy_xz_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_xz_0_zz[i] = 4.0 * g_xxy_xz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (882-888)

    #pragma omp simd aligned(g_xxy_yy_y_xx, g_xxy_yy_y_xy, g_xxy_yy_y_xz, g_xxy_yy_y_yy, g_xxy_yy_y_yz, g_xxy_yy_y_zz, g_y_0_y_0_xx_yy_0_xx, g_y_0_y_0_xx_yy_0_xy, g_y_0_y_0_xx_yy_0_xz, g_y_0_y_0_xx_yy_0_yy, g_y_0_y_0_xx_yy_0_yz, g_y_0_y_0_xx_yy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_xx_yy_0_xx[i] = 4.0 * g_xxy_yy_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_yy_0_xy[i] = 4.0 * g_xxy_yy_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_yy_0_xz[i] = 4.0 * g_xxy_yy_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_yy_0_yy[i] = 4.0 * g_xxy_yy_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_yy_0_yz[i] = 4.0 * g_xxy_yy_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_yy_0_zz[i] = 4.0 * g_xxy_yy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (888-894)

    #pragma omp simd aligned(g_xxy_yz_y_xx, g_xxy_yz_y_xy, g_xxy_yz_y_xz, g_xxy_yz_y_yy, g_xxy_yz_y_yz, g_xxy_yz_y_zz, g_y_0_y_0_xx_yz_0_xx, g_y_0_y_0_xx_yz_0_xy, g_y_0_y_0_xx_yz_0_xz, g_y_0_y_0_xx_yz_0_yy, g_y_0_y_0_xx_yz_0_yz, g_y_0_y_0_xx_yz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_xx_yz_0_xx[i] = 4.0 * g_xxy_yz_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_yz_0_xy[i] = 4.0 * g_xxy_yz_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_yz_0_xz[i] = 4.0 * g_xxy_yz_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_yz_0_yy[i] = 4.0 * g_xxy_yz_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_yz_0_yz[i] = 4.0 * g_xxy_yz_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_yz_0_zz[i] = 4.0 * g_xxy_yz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (894-900)

    #pragma omp simd aligned(g_xxy_zz_y_xx, g_xxy_zz_y_xy, g_xxy_zz_y_xz, g_xxy_zz_y_yy, g_xxy_zz_y_yz, g_xxy_zz_y_zz, g_y_0_y_0_xx_zz_0_xx, g_y_0_y_0_xx_zz_0_xy, g_y_0_y_0_xx_zz_0_xz, g_y_0_y_0_xx_zz_0_yy, g_y_0_y_0_xx_zz_0_yz, g_y_0_y_0_xx_zz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_xx_zz_0_xx[i] = 4.0 * g_xxy_zz_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_zz_0_xy[i] = 4.0 * g_xxy_zz_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_zz_0_xz[i] = 4.0 * g_xxy_zz_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_zz_0_yy[i] = 4.0 * g_xxy_zz_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_zz_0_yz[i] = 4.0 * g_xxy_zz_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_zz_0_zz[i] = 4.0 * g_xxy_zz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (900-906)

    #pragma omp simd aligned(g_x_xx_y_xx, g_x_xx_y_xy, g_x_xx_y_xz, g_x_xx_y_yy, g_x_xx_y_yz, g_x_xx_y_zz, g_xyy_xx_y_xx, g_xyy_xx_y_xy, g_xyy_xx_y_xz, g_xyy_xx_y_yy, g_xyy_xx_y_yz, g_xyy_xx_y_zz, g_y_0_y_0_xy_xx_0_xx, g_y_0_y_0_xy_xx_0_xy, g_y_0_y_0_xy_xx_0_xz, g_y_0_y_0_xy_xx_0_yy, g_y_0_y_0_xy_xx_0_yz, g_y_0_y_0_xy_xx_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_xy_xx_0_xx[i] = -2.0 * g_x_xx_y_xx[i] * c_exps[i] + 4.0 * g_xyy_xx_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_xx_0_xy[i] = -2.0 * g_x_xx_y_xy[i] * c_exps[i] + 4.0 * g_xyy_xx_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_xx_0_xz[i] = -2.0 * g_x_xx_y_xz[i] * c_exps[i] + 4.0 * g_xyy_xx_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_xx_0_yy[i] = -2.0 * g_x_xx_y_yy[i] * c_exps[i] + 4.0 * g_xyy_xx_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_xx_0_yz[i] = -2.0 * g_x_xx_y_yz[i] * c_exps[i] + 4.0 * g_xyy_xx_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_xx_0_zz[i] = -2.0 * g_x_xx_y_zz[i] * c_exps[i] + 4.0 * g_xyy_xx_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (906-912)

    #pragma omp simd aligned(g_x_xy_y_xx, g_x_xy_y_xy, g_x_xy_y_xz, g_x_xy_y_yy, g_x_xy_y_yz, g_x_xy_y_zz, g_xyy_xy_y_xx, g_xyy_xy_y_xy, g_xyy_xy_y_xz, g_xyy_xy_y_yy, g_xyy_xy_y_yz, g_xyy_xy_y_zz, g_y_0_y_0_xy_xy_0_xx, g_y_0_y_0_xy_xy_0_xy, g_y_0_y_0_xy_xy_0_xz, g_y_0_y_0_xy_xy_0_yy, g_y_0_y_0_xy_xy_0_yz, g_y_0_y_0_xy_xy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_xy_xy_0_xx[i] = -2.0 * g_x_xy_y_xx[i] * c_exps[i] + 4.0 * g_xyy_xy_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_xy_0_xy[i] = -2.0 * g_x_xy_y_xy[i] * c_exps[i] + 4.0 * g_xyy_xy_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_xy_0_xz[i] = -2.0 * g_x_xy_y_xz[i] * c_exps[i] + 4.0 * g_xyy_xy_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_xy_0_yy[i] = -2.0 * g_x_xy_y_yy[i] * c_exps[i] + 4.0 * g_xyy_xy_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_xy_0_yz[i] = -2.0 * g_x_xy_y_yz[i] * c_exps[i] + 4.0 * g_xyy_xy_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_xy_0_zz[i] = -2.0 * g_x_xy_y_zz[i] * c_exps[i] + 4.0 * g_xyy_xy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (912-918)

    #pragma omp simd aligned(g_x_xz_y_xx, g_x_xz_y_xy, g_x_xz_y_xz, g_x_xz_y_yy, g_x_xz_y_yz, g_x_xz_y_zz, g_xyy_xz_y_xx, g_xyy_xz_y_xy, g_xyy_xz_y_xz, g_xyy_xz_y_yy, g_xyy_xz_y_yz, g_xyy_xz_y_zz, g_y_0_y_0_xy_xz_0_xx, g_y_0_y_0_xy_xz_0_xy, g_y_0_y_0_xy_xz_0_xz, g_y_0_y_0_xy_xz_0_yy, g_y_0_y_0_xy_xz_0_yz, g_y_0_y_0_xy_xz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_xy_xz_0_xx[i] = -2.0 * g_x_xz_y_xx[i] * c_exps[i] + 4.0 * g_xyy_xz_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_xz_0_xy[i] = -2.0 * g_x_xz_y_xy[i] * c_exps[i] + 4.0 * g_xyy_xz_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_xz_0_xz[i] = -2.0 * g_x_xz_y_xz[i] * c_exps[i] + 4.0 * g_xyy_xz_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_xz_0_yy[i] = -2.0 * g_x_xz_y_yy[i] * c_exps[i] + 4.0 * g_xyy_xz_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_xz_0_yz[i] = -2.0 * g_x_xz_y_yz[i] * c_exps[i] + 4.0 * g_xyy_xz_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_xz_0_zz[i] = -2.0 * g_x_xz_y_zz[i] * c_exps[i] + 4.0 * g_xyy_xz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (918-924)

    #pragma omp simd aligned(g_x_yy_y_xx, g_x_yy_y_xy, g_x_yy_y_xz, g_x_yy_y_yy, g_x_yy_y_yz, g_x_yy_y_zz, g_xyy_yy_y_xx, g_xyy_yy_y_xy, g_xyy_yy_y_xz, g_xyy_yy_y_yy, g_xyy_yy_y_yz, g_xyy_yy_y_zz, g_y_0_y_0_xy_yy_0_xx, g_y_0_y_0_xy_yy_0_xy, g_y_0_y_0_xy_yy_0_xz, g_y_0_y_0_xy_yy_0_yy, g_y_0_y_0_xy_yy_0_yz, g_y_0_y_0_xy_yy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_xy_yy_0_xx[i] = -2.0 * g_x_yy_y_xx[i] * c_exps[i] + 4.0 * g_xyy_yy_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_yy_0_xy[i] = -2.0 * g_x_yy_y_xy[i] * c_exps[i] + 4.0 * g_xyy_yy_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_yy_0_xz[i] = -2.0 * g_x_yy_y_xz[i] * c_exps[i] + 4.0 * g_xyy_yy_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_yy_0_yy[i] = -2.0 * g_x_yy_y_yy[i] * c_exps[i] + 4.0 * g_xyy_yy_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_yy_0_yz[i] = -2.0 * g_x_yy_y_yz[i] * c_exps[i] + 4.0 * g_xyy_yy_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_yy_0_zz[i] = -2.0 * g_x_yy_y_zz[i] * c_exps[i] + 4.0 * g_xyy_yy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (924-930)

    #pragma omp simd aligned(g_x_yz_y_xx, g_x_yz_y_xy, g_x_yz_y_xz, g_x_yz_y_yy, g_x_yz_y_yz, g_x_yz_y_zz, g_xyy_yz_y_xx, g_xyy_yz_y_xy, g_xyy_yz_y_xz, g_xyy_yz_y_yy, g_xyy_yz_y_yz, g_xyy_yz_y_zz, g_y_0_y_0_xy_yz_0_xx, g_y_0_y_0_xy_yz_0_xy, g_y_0_y_0_xy_yz_0_xz, g_y_0_y_0_xy_yz_0_yy, g_y_0_y_0_xy_yz_0_yz, g_y_0_y_0_xy_yz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_xy_yz_0_xx[i] = -2.0 * g_x_yz_y_xx[i] * c_exps[i] + 4.0 * g_xyy_yz_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_yz_0_xy[i] = -2.0 * g_x_yz_y_xy[i] * c_exps[i] + 4.0 * g_xyy_yz_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_yz_0_xz[i] = -2.0 * g_x_yz_y_xz[i] * c_exps[i] + 4.0 * g_xyy_yz_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_yz_0_yy[i] = -2.0 * g_x_yz_y_yy[i] * c_exps[i] + 4.0 * g_xyy_yz_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_yz_0_yz[i] = -2.0 * g_x_yz_y_yz[i] * c_exps[i] + 4.0 * g_xyy_yz_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_yz_0_zz[i] = -2.0 * g_x_yz_y_zz[i] * c_exps[i] + 4.0 * g_xyy_yz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (930-936)

    #pragma omp simd aligned(g_x_zz_y_xx, g_x_zz_y_xy, g_x_zz_y_xz, g_x_zz_y_yy, g_x_zz_y_yz, g_x_zz_y_zz, g_xyy_zz_y_xx, g_xyy_zz_y_xy, g_xyy_zz_y_xz, g_xyy_zz_y_yy, g_xyy_zz_y_yz, g_xyy_zz_y_zz, g_y_0_y_0_xy_zz_0_xx, g_y_0_y_0_xy_zz_0_xy, g_y_0_y_0_xy_zz_0_xz, g_y_0_y_0_xy_zz_0_yy, g_y_0_y_0_xy_zz_0_yz, g_y_0_y_0_xy_zz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_xy_zz_0_xx[i] = -2.0 * g_x_zz_y_xx[i] * c_exps[i] + 4.0 * g_xyy_zz_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_zz_0_xy[i] = -2.0 * g_x_zz_y_xy[i] * c_exps[i] + 4.0 * g_xyy_zz_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_zz_0_xz[i] = -2.0 * g_x_zz_y_xz[i] * c_exps[i] + 4.0 * g_xyy_zz_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_zz_0_yy[i] = -2.0 * g_x_zz_y_yy[i] * c_exps[i] + 4.0 * g_xyy_zz_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_zz_0_yz[i] = -2.0 * g_x_zz_y_yz[i] * c_exps[i] + 4.0 * g_xyy_zz_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_zz_0_zz[i] = -2.0 * g_x_zz_y_zz[i] * c_exps[i] + 4.0 * g_xyy_zz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (936-942)

    #pragma omp simd aligned(g_xyz_xx_y_xx, g_xyz_xx_y_xy, g_xyz_xx_y_xz, g_xyz_xx_y_yy, g_xyz_xx_y_yz, g_xyz_xx_y_zz, g_y_0_y_0_xz_xx_0_xx, g_y_0_y_0_xz_xx_0_xy, g_y_0_y_0_xz_xx_0_xz, g_y_0_y_0_xz_xx_0_yy, g_y_0_y_0_xz_xx_0_yz, g_y_0_y_0_xz_xx_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_xz_xx_0_xx[i] = 4.0 * g_xyz_xx_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_xx_0_xy[i] = 4.0 * g_xyz_xx_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_xx_0_xz[i] = 4.0 * g_xyz_xx_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_xx_0_yy[i] = 4.0 * g_xyz_xx_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_xx_0_yz[i] = 4.0 * g_xyz_xx_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_xx_0_zz[i] = 4.0 * g_xyz_xx_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (942-948)

    #pragma omp simd aligned(g_xyz_xy_y_xx, g_xyz_xy_y_xy, g_xyz_xy_y_xz, g_xyz_xy_y_yy, g_xyz_xy_y_yz, g_xyz_xy_y_zz, g_y_0_y_0_xz_xy_0_xx, g_y_0_y_0_xz_xy_0_xy, g_y_0_y_0_xz_xy_0_xz, g_y_0_y_0_xz_xy_0_yy, g_y_0_y_0_xz_xy_0_yz, g_y_0_y_0_xz_xy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_xz_xy_0_xx[i] = 4.0 * g_xyz_xy_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_xy_0_xy[i] = 4.0 * g_xyz_xy_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_xy_0_xz[i] = 4.0 * g_xyz_xy_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_xy_0_yy[i] = 4.0 * g_xyz_xy_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_xy_0_yz[i] = 4.0 * g_xyz_xy_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_xy_0_zz[i] = 4.0 * g_xyz_xy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (948-954)

    #pragma omp simd aligned(g_xyz_xz_y_xx, g_xyz_xz_y_xy, g_xyz_xz_y_xz, g_xyz_xz_y_yy, g_xyz_xz_y_yz, g_xyz_xz_y_zz, g_y_0_y_0_xz_xz_0_xx, g_y_0_y_0_xz_xz_0_xy, g_y_0_y_0_xz_xz_0_xz, g_y_0_y_0_xz_xz_0_yy, g_y_0_y_0_xz_xz_0_yz, g_y_0_y_0_xz_xz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_xz_xz_0_xx[i] = 4.0 * g_xyz_xz_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_xz_0_xy[i] = 4.0 * g_xyz_xz_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_xz_0_xz[i] = 4.0 * g_xyz_xz_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_xz_0_yy[i] = 4.0 * g_xyz_xz_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_xz_0_yz[i] = 4.0 * g_xyz_xz_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_xz_0_zz[i] = 4.0 * g_xyz_xz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (954-960)

    #pragma omp simd aligned(g_xyz_yy_y_xx, g_xyz_yy_y_xy, g_xyz_yy_y_xz, g_xyz_yy_y_yy, g_xyz_yy_y_yz, g_xyz_yy_y_zz, g_y_0_y_0_xz_yy_0_xx, g_y_0_y_0_xz_yy_0_xy, g_y_0_y_0_xz_yy_0_xz, g_y_0_y_0_xz_yy_0_yy, g_y_0_y_0_xz_yy_0_yz, g_y_0_y_0_xz_yy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_xz_yy_0_xx[i] = 4.0 * g_xyz_yy_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_yy_0_xy[i] = 4.0 * g_xyz_yy_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_yy_0_xz[i] = 4.0 * g_xyz_yy_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_yy_0_yy[i] = 4.0 * g_xyz_yy_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_yy_0_yz[i] = 4.0 * g_xyz_yy_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_yy_0_zz[i] = 4.0 * g_xyz_yy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (960-966)

    #pragma omp simd aligned(g_xyz_yz_y_xx, g_xyz_yz_y_xy, g_xyz_yz_y_xz, g_xyz_yz_y_yy, g_xyz_yz_y_yz, g_xyz_yz_y_zz, g_y_0_y_0_xz_yz_0_xx, g_y_0_y_0_xz_yz_0_xy, g_y_0_y_0_xz_yz_0_xz, g_y_0_y_0_xz_yz_0_yy, g_y_0_y_0_xz_yz_0_yz, g_y_0_y_0_xz_yz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_xz_yz_0_xx[i] = 4.0 * g_xyz_yz_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_yz_0_xy[i] = 4.0 * g_xyz_yz_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_yz_0_xz[i] = 4.0 * g_xyz_yz_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_yz_0_yy[i] = 4.0 * g_xyz_yz_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_yz_0_yz[i] = 4.0 * g_xyz_yz_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_yz_0_zz[i] = 4.0 * g_xyz_yz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (966-972)

    #pragma omp simd aligned(g_xyz_zz_y_xx, g_xyz_zz_y_xy, g_xyz_zz_y_xz, g_xyz_zz_y_yy, g_xyz_zz_y_yz, g_xyz_zz_y_zz, g_y_0_y_0_xz_zz_0_xx, g_y_0_y_0_xz_zz_0_xy, g_y_0_y_0_xz_zz_0_xz, g_y_0_y_0_xz_zz_0_yy, g_y_0_y_0_xz_zz_0_yz, g_y_0_y_0_xz_zz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_xz_zz_0_xx[i] = 4.0 * g_xyz_zz_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_zz_0_xy[i] = 4.0 * g_xyz_zz_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_zz_0_xz[i] = 4.0 * g_xyz_zz_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_zz_0_yy[i] = 4.0 * g_xyz_zz_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_zz_0_yz[i] = 4.0 * g_xyz_zz_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_zz_0_zz[i] = 4.0 * g_xyz_zz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (972-978)

    #pragma omp simd aligned(g_y_0_y_0_yy_xx_0_xx, g_y_0_y_0_yy_xx_0_xy, g_y_0_y_0_yy_xx_0_xz, g_y_0_y_0_yy_xx_0_yy, g_y_0_y_0_yy_xx_0_yz, g_y_0_y_0_yy_xx_0_zz, g_y_xx_y_xx, g_y_xx_y_xy, g_y_xx_y_xz, g_y_xx_y_yy, g_y_xx_y_yz, g_y_xx_y_zz, g_yyy_xx_y_xx, g_yyy_xx_y_xy, g_yyy_xx_y_xz, g_yyy_xx_y_yy, g_yyy_xx_y_yz, g_yyy_xx_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_yy_xx_0_xx[i] = -4.0 * g_y_xx_y_xx[i] * c_exps[i] + 4.0 * g_yyy_xx_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_xx_0_xy[i] = -4.0 * g_y_xx_y_xy[i] * c_exps[i] + 4.0 * g_yyy_xx_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_xx_0_xz[i] = -4.0 * g_y_xx_y_xz[i] * c_exps[i] + 4.0 * g_yyy_xx_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_xx_0_yy[i] = -4.0 * g_y_xx_y_yy[i] * c_exps[i] + 4.0 * g_yyy_xx_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_xx_0_yz[i] = -4.0 * g_y_xx_y_yz[i] * c_exps[i] + 4.0 * g_yyy_xx_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_xx_0_zz[i] = -4.0 * g_y_xx_y_zz[i] * c_exps[i] + 4.0 * g_yyy_xx_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (978-984)

    #pragma omp simd aligned(g_y_0_y_0_yy_xy_0_xx, g_y_0_y_0_yy_xy_0_xy, g_y_0_y_0_yy_xy_0_xz, g_y_0_y_0_yy_xy_0_yy, g_y_0_y_0_yy_xy_0_yz, g_y_0_y_0_yy_xy_0_zz, g_y_xy_y_xx, g_y_xy_y_xy, g_y_xy_y_xz, g_y_xy_y_yy, g_y_xy_y_yz, g_y_xy_y_zz, g_yyy_xy_y_xx, g_yyy_xy_y_xy, g_yyy_xy_y_xz, g_yyy_xy_y_yy, g_yyy_xy_y_yz, g_yyy_xy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_yy_xy_0_xx[i] = -4.0 * g_y_xy_y_xx[i] * c_exps[i] + 4.0 * g_yyy_xy_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_xy_0_xy[i] = -4.0 * g_y_xy_y_xy[i] * c_exps[i] + 4.0 * g_yyy_xy_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_xy_0_xz[i] = -4.0 * g_y_xy_y_xz[i] * c_exps[i] + 4.0 * g_yyy_xy_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_xy_0_yy[i] = -4.0 * g_y_xy_y_yy[i] * c_exps[i] + 4.0 * g_yyy_xy_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_xy_0_yz[i] = -4.0 * g_y_xy_y_yz[i] * c_exps[i] + 4.0 * g_yyy_xy_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_xy_0_zz[i] = -4.0 * g_y_xy_y_zz[i] * c_exps[i] + 4.0 * g_yyy_xy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (984-990)

    #pragma omp simd aligned(g_y_0_y_0_yy_xz_0_xx, g_y_0_y_0_yy_xz_0_xy, g_y_0_y_0_yy_xz_0_xz, g_y_0_y_0_yy_xz_0_yy, g_y_0_y_0_yy_xz_0_yz, g_y_0_y_0_yy_xz_0_zz, g_y_xz_y_xx, g_y_xz_y_xy, g_y_xz_y_xz, g_y_xz_y_yy, g_y_xz_y_yz, g_y_xz_y_zz, g_yyy_xz_y_xx, g_yyy_xz_y_xy, g_yyy_xz_y_xz, g_yyy_xz_y_yy, g_yyy_xz_y_yz, g_yyy_xz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_yy_xz_0_xx[i] = -4.0 * g_y_xz_y_xx[i] * c_exps[i] + 4.0 * g_yyy_xz_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_xz_0_xy[i] = -4.0 * g_y_xz_y_xy[i] * c_exps[i] + 4.0 * g_yyy_xz_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_xz_0_xz[i] = -4.0 * g_y_xz_y_xz[i] * c_exps[i] + 4.0 * g_yyy_xz_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_xz_0_yy[i] = -4.0 * g_y_xz_y_yy[i] * c_exps[i] + 4.0 * g_yyy_xz_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_xz_0_yz[i] = -4.0 * g_y_xz_y_yz[i] * c_exps[i] + 4.0 * g_yyy_xz_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_xz_0_zz[i] = -4.0 * g_y_xz_y_zz[i] * c_exps[i] + 4.0 * g_yyy_xz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (990-996)

    #pragma omp simd aligned(g_y_0_y_0_yy_yy_0_xx, g_y_0_y_0_yy_yy_0_xy, g_y_0_y_0_yy_yy_0_xz, g_y_0_y_0_yy_yy_0_yy, g_y_0_y_0_yy_yy_0_yz, g_y_0_y_0_yy_yy_0_zz, g_y_yy_y_xx, g_y_yy_y_xy, g_y_yy_y_xz, g_y_yy_y_yy, g_y_yy_y_yz, g_y_yy_y_zz, g_yyy_yy_y_xx, g_yyy_yy_y_xy, g_yyy_yy_y_xz, g_yyy_yy_y_yy, g_yyy_yy_y_yz, g_yyy_yy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_yy_yy_0_xx[i] = -4.0 * g_y_yy_y_xx[i] * c_exps[i] + 4.0 * g_yyy_yy_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_yy_0_xy[i] = -4.0 * g_y_yy_y_xy[i] * c_exps[i] + 4.0 * g_yyy_yy_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_yy_0_xz[i] = -4.0 * g_y_yy_y_xz[i] * c_exps[i] + 4.0 * g_yyy_yy_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_yy_0_yy[i] = -4.0 * g_y_yy_y_yy[i] * c_exps[i] + 4.0 * g_yyy_yy_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_yy_0_yz[i] = -4.0 * g_y_yy_y_yz[i] * c_exps[i] + 4.0 * g_yyy_yy_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_yy_0_zz[i] = -4.0 * g_y_yy_y_zz[i] * c_exps[i] + 4.0 * g_yyy_yy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (996-1002)

    #pragma omp simd aligned(g_y_0_y_0_yy_yz_0_xx, g_y_0_y_0_yy_yz_0_xy, g_y_0_y_0_yy_yz_0_xz, g_y_0_y_0_yy_yz_0_yy, g_y_0_y_0_yy_yz_0_yz, g_y_0_y_0_yy_yz_0_zz, g_y_yz_y_xx, g_y_yz_y_xy, g_y_yz_y_xz, g_y_yz_y_yy, g_y_yz_y_yz, g_y_yz_y_zz, g_yyy_yz_y_xx, g_yyy_yz_y_xy, g_yyy_yz_y_xz, g_yyy_yz_y_yy, g_yyy_yz_y_yz, g_yyy_yz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_yy_yz_0_xx[i] = -4.0 * g_y_yz_y_xx[i] * c_exps[i] + 4.0 * g_yyy_yz_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_yz_0_xy[i] = -4.0 * g_y_yz_y_xy[i] * c_exps[i] + 4.0 * g_yyy_yz_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_yz_0_xz[i] = -4.0 * g_y_yz_y_xz[i] * c_exps[i] + 4.0 * g_yyy_yz_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_yz_0_yy[i] = -4.0 * g_y_yz_y_yy[i] * c_exps[i] + 4.0 * g_yyy_yz_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_yz_0_yz[i] = -4.0 * g_y_yz_y_yz[i] * c_exps[i] + 4.0 * g_yyy_yz_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_yz_0_zz[i] = -4.0 * g_y_yz_y_zz[i] * c_exps[i] + 4.0 * g_yyy_yz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1002-1008)

    #pragma omp simd aligned(g_y_0_y_0_yy_zz_0_xx, g_y_0_y_0_yy_zz_0_xy, g_y_0_y_0_yy_zz_0_xz, g_y_0_y_0_yy_zz_0_yy, g_y_0_y_0_yy_zz_0_yz, g_y_0_y_0_yy_zz_0_zz, g_y_zz_y_xx, g_y_zz_y_xy, g_y_zz_y_xz, g_y_zz_y_yy, g_y_zz_y_yz, g_y_zz_y_zz, g_yyy_zz_y_xx, g_yyy_zz_y_xy, g_yyy_zz_y_xz, g_yyy_zz_y_yy, g_yyy_zz_y_yz, g_yyy_zz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_yy_zz_0_xx[i] = -4.0 * g_y_zz_y_xx[i] * c_exps[i] + 4.0 * g_yyy_zz_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_zz_0_xy[i] = -4.0 * g_y_zz_y_xy[i] * c_exps[i] + 4.0 * g_yyy_zz_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_zz_0_xz[i] = -4.0 * g_y_zz_y_xz[i] * c_exps[i] + 4.0 * g_yyy_zz_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_zz_0_yy[i] = -4.0 * g_y_zz_y_yy[i] * c_exps[i] + 4.0 * g_yyy_zz_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_zz_0_yz[i] = -4.0 * g_y_zz_y_yz[i] * c_exps[i] + 4.0 * g_yyy_zz_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_zz_0_zz[i] = -4.0 * g_y_zz_y_zz[i] * c_exps[i] + 4.0 * g_yyy_zz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1008-1014)

    #pragma omp simd aligned(g_y_0_y_0_yz_xx_0_xx, g_y_0_y_0_yz_xx_0_xy, g_y_0_y_0_yz_xx_0_xz, g_y_0_y_0_yz_xx_0_yy, g_y_0_y_0_yz_xx_0_yz, g_y_0_y_0_yz_xx_0_zz, g_yyz_xx_y_xx, g_yyz_xx_y_xy, g_yyz_xx_y_xz, g_yyz_xx_y_yy, g_yyz_xx_y_yz, g_yyz_xx_y_zz, g_z_xx_y_xx, g_z_xx_y_xy, g_z_xx_y_xz, g_z_xx_y_yy, g_z_xx_y_yz, g_z_xx_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_yz_xx_0_xx[i] = -2.0 * g_z_xx_y_xx[i] * c_exps[i] + 4.0 * g_yyz_xx_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_xx_0_xy[i] = -2.0 * g_z_xx_y_xy[i] * c_exps[i] + 4.0 * g_yyz_xx_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_xx_0_xz[i] = -2.0 * g_z_xx_y_xz[i] * c_exps[i] + 4.0 * g_yyz_xx_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_xx_0_yy[i] = -2.0 * g_z_xx_y_yy[i] * c_exps[i] + 4.0 * g_yyz_xx_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_xx_0_yz[i] = -2.0 * g_z_xx_y_yz[i] * c_exps[i] + 4.0 * g_yyz_xx_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_xx_0_zz[i] = -2.0 * g_z_xx_y_zz[i] * c_exps[i] + 4.0 * g_yyz_xx_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1014-1020)

    #pragma omp simd aligned(g_y_0_y_0_yz_xy_0_xx, g_y_0_y_0_yz_xy_0_xy, g_y_0_y_0_yz_xy_0_xz, g_y_0_y_0_yz_xy_0_yy, g_y_0_y_0_yz_xy_0_yz, g_y_0_y_0_yz_xy_0_zz, g_yyz_xy_y_xx, g_yyz_xy_y_xy, g_yyz_xy_y_xz, g_yyz_xy_y_yy, g_yyz_xy_y_yz, g_yyz_xy_y_zz, g_z_xy_y_xx, g_z_xy_y_xy, g_z_xy_y_xz, g_z_xy_y_yy, g_z_xy_y_yz, g_z_xy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_yz_xy_0_xx[i] = -2.0 * g_z_xy_y_xx[i] * c_exps[i] + 4.0 * g_yyz_xy_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_xy_0_xy[i] = -2.0 * g_z_xy_y_xy[i] * c_exps[i] + 4.0 * g_yyz_xy_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_xy_0_xz[i] = -2.0 * g_z_xy_y_xz[i] * c_exps[i] + 4.0 * g_yyz_xy_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_xy_0_yy[i] = -2.0 * g_z_xy_y_yy[i] * c_exps[i] + 4.0 * g_yyz_xy_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_xy_0_yz[i] = -2.0 * g_z_xy_y_yz[i] * c_exps[i] + 4.0 * g_yyz_xy_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_xy_0_zz[i] = -2.0 * g_z_xy_y_zz[i] * c_exps[i] + 4.0 * g_yyz_xy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1020-1026)

    #pragma omp simd aligned(g_y_0_y_0_yz_xz_0_xx, g_y_0_y_0_yz_xz_0_xy, g_y_0_y_0_yz_xz_0_xz, g_y_0_y_0_yz_xz_0_yy, g_y_0_y_0_yz_xz_0_yz, g_y_0_y_0_yz_xz_0_zz, g_yyz_xz_y_xx, g_yyz_xz_y_xy, g_yyz_xz_y_xz, g_yyz_xz_y_yy, g_yyz_xz_y_yz, g_yyz_xz_y_zz, g_z_xz_y_xx, g_z_xz_y_xy, g_z_xz_y_xz, g_z_xz_y_yy, g_z_xz_y_yz, g_z_xz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_yz_xz_0_xx[i] = -2.0 * g_z_xz_y_xx[i] * c_exps[i] + 4.0 * g_yyz_xz_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_xz_0_xy[i] = -2.0 * g_z_xz_y_xy[i] * c_exps[i] + 4.0 * g_yyz_xz_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_xz_0_xz[i] = -2.0 * g_z_xz_y_xz[i] * c_exps[i] + 4.0 * g_yyz_xz_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_xz_0_yy[i] = -2.0 * g_z_xz_y_yy[i] * c_exps[i] + 4.0 * g_yyz_xz_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_xz_0_yz[i] = -2.0 * g_z_xz_y_yz[i] * c_exps[i] + 4.0 * g_yyz_xz_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_xz_0_zz[i] = -2.0 * g_z_xz_y_zz[i] * c_exps[i] + 4.0 * g_yyz_xz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1026-1032)

    #pragma omp simd aligned(g_y_0_y_0_yz_yy_0_xx, g_y_0_y_0_yz_yy_0_xy, g_y_0_y_0_yz_yy_0_xz, g_y_0_y_0_yz_yy_0_yy, g_y_0_y_0_yz_yy_0_yz, g_y_0_y_0_yz_yy_0_zz, g_yyz_yy_y_xx, g_yyz_yy_y_xy, g_yyz_yy_y_xz, g_yyz_yy_y_yy, g_yyz_yy_y_yz, g_yyz_yy_y_zz, g_z_yy_y_xx, g_z_yy_y_xy, g_z_yy_y_xz, g_z_yy_y_yy, g_z_yy_y_yz, g_z_yy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_yz_yy_0_xx[i] = -2.0 * g_z_yy_y_xx[i] * c_exps[i] + 4.0 * g_yyz_yy_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_yy_0_xy[i] = -2.0 * g_z_yy_y_xy[i] * c_exps[i] + 4.0 * g_yyz_yy_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_yy_0_xz[i] = -2.0 * g_z_yy_y_xz[i] * c_exps[i] + 4.0 * g_yyz_yy_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_yy_0_yy[i] = -2.0 * g_z_yy_y_yy[i] * c_exps[i] + 4.0 * g_yyz_yy_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_yy_0_yz[i] = -2.0 * g_z_yy_y_yz[i] * c_exps[i] + 4.0 * g_yyz_yy_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_yy_0_zz[i] = -2.0 * g_z_yy_y_zz[i] * c_exps[i] + 4.0 * g_yyz_yy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1032-1038)

    #pragma omp simd aligned(g_y_0_y_0_yz_yz_0_xx, g_y_0_y_0_yz_yz_0_xy, g_y_0_y_0_yz_yz_0_xz, g_y_0_y_0_yz_yz_0_yy, g_y_0_y_0_yz_yz_0_yz, g_y_0_y_0_yz_yz_0_zz, g_yyz_yz_y_xx, g_yyz_yz_y_xy, g_yyz_yz_y_xz, g_yyz_yz_y_yy, g_yyz_yz_y_yz, g_yyz_yz_y_zz, g_z_yz_y_xx, g_z_yz_y_xy, g_z_yz_y_xz, g_z_yz_y_yy, g_z_yz_y_yz, g_z_yz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_yz_yz_0_xx[i] = -2.0 * g_z_yz_y_xx[i] * c_exps[i] + 4.0 * g_yyz_yz_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_yz_0_xy[i] = -2.0 * g_z_yz_y_xy[i] * c_exps[i] + 4.0 * g_yyz_yz_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_yz_0_xz[i] = -2.0 * g_z_yz_y_xz[i] * c_exps[i] + 4.0 * g_yyz_yz_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_yz_0_yy[i] = -2.0 * g_z_yz_y_yy[i] * c_exps[i] + 4.0 * g_yyz_yz_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_yz_0_yz[i] = -2.0 * g_z_yz_y_yz[i] * c_exps[i] + 4.0 * g_yyz_yz_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_yz_0_zz[i] = -2.0 * g_z_yz_y_zz[i] * c_exps[i] + 4.0 * g_yyz_yz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1038-1044)

    #pragma omp simd aligned(g_y_0_y_0_yz_zz_0_xx, g_y_0_y_0_yz_zz_0_xy, g_y_0_y_0_yz_zz_0_xz, g_y_0_y_0_yz_zz_0_yy, g_y_0_y_0_yz_zz_0_yz, g_y_0_y_0_yz_zz_0_zz, g_yyz_zz_y_xx, g_yyz_zz_y_xy, g_yyz_zz_y_xz, g_yyz_zz_y_yy, g_yyz_zz_y_yz, g_yyz_zz_y_zz, g_z_zz_y_xx, g_z_zz_y_xy, g_z_zz_y_xz, g_z_zz_y_yy, g_z_zz_y_yz, g_z_zz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_yz_zz_0_xx[i] = -2.0 * g_z_zz_y_xx[i] * c_exps[i] + 4.0 * g_yyz_zz_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_zz_0_xy[i] = -2.0 * g_z_zz_y_xy[i] * c_exps[i] + 4.0 * g_yyz_zz_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_zz_0_xz[i] = -2.0 * g_z_zz_y_xz[i] * c_exps[i] + 4.0 * g_yyz_zz_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_zz_0_yy[i] = -2.0 * g_z_zz_y_yy[i] * c_exps[i] + 4.0 * g_yyz_zz_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_zz_0_yz[i] = -2.0 * g_z_zz_y_yz[i] * c_exps[i] + 4.0 * g_yyz_zz_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_zz_0_zz[i] = -2.0 * g_z_zz_y_zz[i] * c_exps[i] + 4.0 * g_yyz_zz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1044-1050)

    #pragma omp simd aligned(g_y_0_y_0_zz_xx_0_xx, g_y_0_y_0_zz_xx_0_xy, g_y_0_y_0_zz_xx_0_xz, g_y_0_y_0_zz_xx_0_yy, g_y_0_y_0_zz_xx_0_yz, g_y_0_y_0_zz_xx_0_zz, g_yzz_xx_y_xx, g_yzz_xx_y_xy, g_yzz_xx_y_xz, g_yzz_xx_y_yy, g_yzz_xx_y_yz, g_yzz_xx_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_zz_xx_0_xx[i] = 4.0 * g_yzz_xx_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_xx_0_xy[i] = 4.0 * g_yzz_xx_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_xx_0_xz[i] = 4.0 * g_yzz_xx_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_xx_0_yy[i] = 4.0 * g_yzz_xx_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_xx_0_yz[i] = 4.0 * g_yzz_xx_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_xx_0_zz[i] = 4.0 * g_yzz_xx_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1050-1056)

    #pragma omp simd aligned(g_y_0_y_0_zz_xy_0_xx, g_y_0_y_0_zz_xy_0_xy, g_y_0_y_0_zz_xy_0_xz, g_y_0_y_0_zz_xy_0_yy, g_y_0_y_0_zz_xy_0_yz, g_y_0_y_0_zz_xy_0_zz, g_yzz_xy_y_xx, g_yzz_xy_y_xy, g_yzz_xy_y_xz, g_yzz_xy_y_yy, g_yzz_xy_y_yz, g_yzz_xy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_zz_xy_0_xx[i] = 4.0 * g_yzz_xy_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_xy_0_xy[i] = 4.0 * g_yzz_xy_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_xy_0_xz[i] = 4.0 * g_yzz_xy_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_xy_0_yy[i] = 4.0 * g_yzz_xy_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_xy_0_yz[i] = 4.0 * g_yzz_xy_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_xy_0_zz[i] = 4.0 * g_yzz_xy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1056-1062)

    #pragma omp simd aligned(g_y_0_y_0_zz_xz_0_xx, g_y_0_y_0_zz_xz_0_xy, g_y_0_y_0_zz_xz_0_xz, g_y_0_y_0_zz_xz_0_yy, g_y_0_y_0_zz_xz_0_yz, g_y_0_y_0_zz_xz_0_zz, g_yzz_xz_y_xx, g_yzz_xz_y_xy, g_yzz_xz_y_xz, g_yzz_xz_y_yy, g_yzz_xz_y_yz, g_yzz_xz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_zz_xz_0_xx[i] = 4.0 * g_yzz_xz_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_xz_0_xy[i] = 4.0 * g_yzz_xz_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_xz_0_xz[i] = 4.0 * g_yzz_xz_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_xz_0_yy[i] = 4.0 * g_yzz_xz_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_xz_0_yz[i] = 4.0 * g_yzz_xz_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_xz_0_zz[i] = 4.0 * g_yzz_xz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1062-1068)

    #pragma omp simd aligned(g_y_0_y_0_zz_yy_0_xx, g_y_0_y_0_zz_yy_0_xy, g_y_0_y_0_zz_yy_0_xz, g_y_0_y_0_zz_yy_0_yy, g_y_0_y_0_zz_yy_0_yz, g_y_0_y_0_zz_yy_0_zz, g_yzz_yy_y_xx, g_yzz_yy_y_xy, g_yzz_yy_y_xz, g_yzz_yy_y_yy, g_yzz_yy_y_yz, g_yzz_yy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_zz_yy_0_xx[i] = 4.0 * g_yzz_yy_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_yy_0_xy[i] = 4.0 * g_yzz_yy_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_yy_0_xz[i] = 4.0 * g_yzz_yy_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_yy_0_yy[i] = 4.0 * g_yzz_yy_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_yy_0_yz[i] = 4.0 * g_yzz_yy_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_yy_0_zz[i] = 4.0 * g_yzz_yy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1068-1074)

    #pragma omp simd aligned(g_y_0_y_0_zz_yz_0_xx, g_y_0_y_0_zz_yz_0_xy, g_y_0_y_0_zz_yz_0_xz, g_y_0_y_0_zz_yz_0_yy, g_y_0_y_0_zz_yz_0_yz, g_y_0_y_0_zz_yz_0_zz, g_yzz_yz_y_xx, g_yzz_yz_y_xy, g_yzz_yz_y_xz, g_yzz_yz_y_yy, g_yzz_yz_y_yz, g_yzz_yz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_zz_yz_0_xx[i] = 4.0 * g_yzz_yz_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_yz_0_xy[i] = 4.0 * g_yzz_yz_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_yz_0_xz[i] = 4.0 * g_yzz_yz_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_yz_0_yy[i] = 4.0 * g_yzz_yz_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_yz_0_yz[i] = 4.0 * g_yzz_yz_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_yz_0_zz[i] = 4.0 * g_yzz_yz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1074-1080)

    #pragma omp simd aligned(g_y_0_y_0_zz_zz_0_xx, g_y_0_y_0_zz_zz_0_xy, g_y_0_y_0_zz_zz_0_xz, g_y_0_y_0_zz_zz_0_yy, g_y_0_y_0_zz_zz_0_yz, g_y_0_y_0_zz_zz_0_zz, g_yzz_zz_y_xx, g_yzz_zz_y_xy, g_yzz_zz_y_xz, g_yzz_zz_y_yy, g_yzz_zz_y_yz, g_yzz_zz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_zz_zz_0_xx[i] = 4.0 * g_yzz_zz_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_zz_0_xy[i] = 4.0 * g_yzz_zz_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_zz_0_xz[i] = 4.0 * g_yzz_zz_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_zz_0_yy[i] = 4.0 * g_yzz_zz_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_zz_0_yz[i] = 4.0 * g_yzz_zz_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_zz_0_zz[i] = 4.0 * g_yzz_zz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1080-1086)

    #pragma omp simd aligned(g_xxy_xx_z_xx, g_xxy_xx_z_xy, g_xxy_xx_z_xz, g_xxy_xx_z_yy, g_xxy_xx_z_yz, g_xxy_xx_z_zz, g_y_0_z_0_xx_xx_0_xx, g_y_0_z_0_xx_xx_0_xy, g_y_0_z_0_xx_xx_0_xz, g_y_0_z_0_xx_xx_0_yy, g_y_0_z_0_xx_xx_0_yz, g_y_0_z_0_xx_xx_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_xx_xx_0_xx[i] = 4.0 * g_xxy_xx_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_xx_0_xy[i] = 4.0 * g_xxy_xx_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_xx_0_xz[i] = 4.0 * g_xxy_xx_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_xx_0_yy[i] = 4.0 * g_xxy_xx_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_xx_0_yz[i] = 4.0 * g_xxy_xx_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_xx_0_zz[i] = 4.0 * g_xxy_xx_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1086-1092)

    #pragma omp simd aligned(g_xxy_xy_z_xx, g_xxy_xy_z_xy, g_xxy_xy_z_xz, g_xxy_xy_z_yy, g_xxy_xy_z_yz, g_xxy_xy_z_zz, g_y_0_z_0_xx_xy_0_xx, g_y_0_z_0_xx_xy_0_xy, g_y_0_z_0_xx_xy_0_xz, g_y_0_z_0_xx_xy_0_yy, g_y_0_z_0_xx_xy_0_yz, g_y_0_z_0_xx_xy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_xx_xy_0_xx[i] = 4.0 * g_xxy_xy_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_xy_0_xy[i] = 4.0 * g_xxy_xy_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_xy_0_xz[i] = 4.0 * g_xxy_xy_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_xy_0_yy[i] = 4.0 * g_xxy_xy_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_xy_0_yz[i] = 4.0 * g_xxy_xy_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_xy_0_zz[i] = 4.0 * g_xxy_xy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1092-1098)

    #pragma omp simd aligned(g_xxy_xz_z_xx, g_xxy_xz_z_xy, g_xxy_xz_z_xz, g_xxy_xz_z_yy, g_xxy_xz_z_yz, g_xxy_xz_z_zz, g_y_0_z_0_xx_xz_0_xx, g_y_0_z_0_xx_xz_0_xy, g_y_0_z_0_xx_xz_0_xz, g_y_0_z_0_xx_xz_0_yy, g_y_0_z_0_xx_xz_0_yz, g_y_0_z_0_xx_xz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_xx_xz_0_xx[i] = 4.0 * g_xxy_xz_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_xz_0_xy[i] = 4.0 * g_xxy_xz_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_xz_0_xz[i] = 4.0 * g_xxy_xz_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_xz_0_yy[i] = 4.0 * g_xxy_xz_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_xz_0_yz[i] = 4.0 * g_xxy_xz_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_xz_0_zz[i] = 4.0 * g_xxy_xz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1098-1104)

    #pragma omp simd aligned(g_xxy_yy_z_xx, g_xxy_yy_z_xy, g_xxy_yy_z_xz, g_xxy_yy_z_yy, g_xxy_yy_z_yz, g_xxy_yy_z_zz, g_y_0_z_0_xx_yy_0_xx, g_y_0_z_0_xx_yy_0_xy, g_y_0_z_0_xx_yy_0_xz, g_y_0_z_0_xx_yy_0_yy, g_y_0_z_0_xx_yy_0_yz, g_y_0_z_0_xx_yy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_xx_yy_0_xx[i] = 4.0 * g_xxy_yy_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_yy_0_xy[i] = 4.0 * g_xxy_yy_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_yy_0_xz[i] = 4.0 * g_xxy_yy_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_yy_0_yy[i] = 4.0 * g_xxy_yy_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_yy_0_yz[i] = 4.0 * g_xxy_yy_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_yy_0_zz[i] = 4.0 * g_xxy_yy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1104-1110)

    #pragma omp simd aligned(g_xxy_yz_z_xx, g_xxy_yz_z_xy, g_xxy_yz_z_xz, g_xxy_yz_z_yy, g_xxy_yz_z_yz, g_xxy_yz_z_zz, g_y_0_z_0_xx_yz_0_xx, g_y_0_z_0_xx_yz_0_xy, g_y_0_z_0_xx_yz_0_xz, g_y_0_z_0_xx_yz_0_yy, g_y_0_z_0_xx_yz_0_yz, g_y_0_z_0_xx_yz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_xx_yz_0_xx[i] = 4.0 * g_xxy_yz_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_yz_0_xy[i] = 4.0 * g_xxy_yz_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_yz_0_xz[i] = 4.0 * g_xxy_yz_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_yz_0_yy[i] = 4.0 * g_xxy_yz_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_yz_0_yz[i] = 4.0 * g_xxy_yz_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_yz_0_zz[i] = 4.0 * g_xxy_yz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1110-1116)

    #pragma omp simd aligned(g_xxy_zz_z_xx, g_xxy_zz_z_xy, g_xxy_zz_z_xz, g_xxy_zz_z_yy, g_xxy_zz_z_yz, g_xxy_zz_z_zz, g_y_0_z_0_xx_zz_0_xx, g_y_0_z_0_xx_zz_0_xy, g_y_0_z_0_xx_zz_0_xz, g_y_0_z_0_xx_zz_0_yy, g_y_0_z_0_xx_zz_0_yz, g_y_0_z_0_xx_zz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_xx_zz_0_xx[i] = 4.0 * g_xxy_zz_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_zz_0_xy[i] = 4.0 * g_xxy_zz_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_zz_0_xz[i] = 4.0 * g_xxy_zz_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_zz_0_yy[i] = 4.0 * g_xxy_zz_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_zz_0_yz[i] = 4.0 * g_xxy_zz_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_zz_0_zz[i] = 4.0 * g_xxy_zz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1116-1122)

    #pragma omp simd aligned(g_x_xx_z_xx, g_x_xx_z_xy, g_x_xx_z_xz, g_x_xx_z_yy, g_x_xx_z_yz, g_x_xx_z_zz, g_xyy_xx_z_xx, g_xyy_xx_z_xy, g_xyy_xx_z_xz, g_xyy_xx_z_yy, g_xyy_xx_z_yz, g_xyy_xx_z_zz, g_y_0_z_0_xy_xx_0_xx, g_y_0_z_0_xy_xx_0_xy, g_y_0_z_0_xy_xx_0_xz, g_y_0_z_0_xy_xx_0_yy, g_y_0_z_0_xy_xx_0_yz, g_y_0_z_0_xy_xx_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_xy_xx_0_xx[i] = -2.0 * g_x_xx_z_xx[i] * c_exps[i] + 4.0 * g_xyy_xx_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_xx_0_xy[i] = -2.0 * g_x_xx_z_xy[i] * c_exps[i] + 4.0 * g_xyy_xx_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_xx_0_xz[i] = -2.0 * g_x_xx_z_xz[i] * c_exps[i] + 4.0 * g_xyy_xx_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_xx_0_yy[i] = -2.0 * g_x_xx_z_yy[i] * c_exps[i] + 4.0 * g_xyy_xx_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_xx_0_yz[i] = -2.0 * g_x_xx_z_yz[i] * c_exps[i] + 4.0 * g_xyy_xx_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_xx_0_zz[i] = -2.0 * g_x_xx_z_zz[i] * c_exps[i] + 4.0 * g_xyy_xx_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1122-1128)

    #pragma omp simd aligned(g_x_xy_z_xx, g_x_xy_z_xy, g_x_xy_z_xz, g_x_xy_z_yy, g_x_xy_z_yz, g_x_xy_z_zz, g_xyy_xy_z_xx, g_xyy_xy_z_xy, g_xyy_xy_z_xz, g_xyy_xy_z_yy, g_xyy_xy_z_yz, g_xyy_xy_z_zz, g_y_0_z_0_xy_xy_0_xx, g_y_0_z_0_xy_xy_0_xy, g_y_0_z_0_xy_xy_0_xz, g_y_0_z_0_xy_xy_0_yy, g_y_0_z_0_xy_xy_0_yz, g_y_0_z_0_xy_xy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_xy_xy_0_xx[i] = -2.0 * g_x_xy_z_xx[i] * c_exps[i] + 4.0 * g_xyy_xy_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_xy_0_xy[i] = -2.0 * g_x_xy_z_xy[i] * c_exps[i] + 4.0 * g_xyy_xy_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_xy_0_xz[i] = -2.0 * g_x_xy_z_xz[i] * c_exps[i] + 4.0 * g_xyy_xy_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_xy_0_yy[i] = -2.0 * g_x_xy_z_yy[i] * c_exps[i] + 4.0 * g_xyy_xy_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_xy_0_yz[i] = -2.0 * g_x_xy_z_yz[i] * c_exps[i] + 4.0 * g_xyy_xy_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_xy_0_zz[i] = -2.0 * g_x_xy_z_zz[i] * c_exps[i] + 4.0 * g_xyy_xy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1128-1134)

    #pragma omp simd aligned(g_x_xz_z_xx, g_x_xz_z_xy, g_x_xz_z_xz, g_x_xz_z_yy, g_x_xz_z_yz, g_x_xz_z_zz, g_xyy_xz_z_xx, g_xyy_xz_z_xy, g_xyy_xz_z_xz, g_xyy_xz_z_yy, g_xyy_xz_z_yz, g_xyy_xz_z_zz, g_y_0_z_0_xy_xz_0_xx, g_y_0_z_0_xy_xz_0_xy, g_y_0_z_0_xy_xz_0_xz, g_y_0_z_0_xy_xz_0_yy, g_y_0_z_0_xy_xz_0_yz, g_y_0_z_0_xy_xz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_xy_xz_0_xx[i] = -2.0 * g_x_xz_z_xx[i] * c_exps[i] + 4.0 * g_xyy_xz_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_xz_0_xy[i] = -2.0 * g_x_xz_z_xy[i] * c_exps[i] + 4.0 * g_xyy_xz_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_xz_0_xz[i] = -2.0 * g_x_xz_z_xz[i] * c_exps[i] + 4.0 * g_xyy_xz_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_xz_0_yy[i] = -2.0 * g_x_xz_z_yy[i] * c_exps[i] + 4.0 * g_xyy_xz_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_xz_0_yz[i] = -2.0 * g_x_xz_z_yz[i] * c_exps[i] + 4.0 * g_xyy_xz_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_xz_0_zz[i] = -2.0 * g_x_xz_z_zz[i] * c_exps[i] + 4.0 * g_xyy_xz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1134-1140)

    #pragma omp simd aligned(g_x_yy_z_xx, g_x_yy_z_xy, g_x_yy_z_xz, g_x_yy_z_yy, g_x_yy_z_yz, g_x_yy_z_zz, g_xyy_yy_z_xx, g_xyy_yy_z_xy, g_xyy_yy_z_xz, g_xyy_yy_z_yy, g_xyy_yy_z_yz, g_xyy_yy_z_zz, g_y_0_z_0_xy_yy_0_xx, g_y_0_z_0_xy_yy_0_xy, g_y_0_z_0_xy_yy_0_xz, g_y_0_z_0_xy_yy_0_yy, g_y_0_z_0_xy_yy_0_yz, g_y_0_z_0_xy_yy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_xy_yy_0_xx[i] = -2.0 * g_x_yy_z_xx[i] * c_exps[i] + 4.0 * g_xyy_yy_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_yy_0_xy[i] = -2.0 * g_x_yy_z_xy[i] * c_exps[i] + 4.0 * g_xyy_yy_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_yy_0_xz[i] = -2.0 * g_x_yy_z_xz[i] * c_exps[i] + 4.0 * g_xyy_yy_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_yy_0_yy[i] = -2.0 * g_x_yy_z_yy[i] * c_exps[i] + 4.0 * g_xyy_yy_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_yy_0_yz[i] = -2.0 * g_x_yy_z_yz[i] * c_exps[i] + 4.0 * g_xyy_yy_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_yy_0_zz[i] = -2.0 * g_x_yy_z_zz[i] * c_exps[i] + 4.0 * g_xyy_yy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1140-1146)

    #pragma omp simd aligned(g_x_yz_z_xx, g_x_yz_z_xy, g_x_yz_z_xz, g_x_yz_z_yy, g_x_yz_z_yz, g_x_yz_z_zz, g_xyy_yz_z_xx, g_xyy_yz_z_xy, g_xyy_yz_z_xz, g_xyy_yz_z_yy, g_xyy_yz_z_yz, g_xyy_yz_z_zz, g_y_0_z_0_xy_yz_0_xx, g_y_0_z_0_xy_yz_0_xy, g_y_0_z_0_xy_yz_0_xz, g_y_0_z_0_xy_yz_0_yy, g_y_0_z_0_xy_yz_0_yz, g_y_0_z_0_xy_yz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_xy_yz_0_xx[i] = -2.0 * g_x_yz_z_xx[i] * c_exps[i] + 4.0 * g_xyy_yz_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_yz_0_xy[i] = -2.0 * g_x_yz_z_xy[i] * c_exps[i] + 4.0 * g_xyy_yz_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_yz_0_xz[i] = -2.0 * g_x_yz_z_xz[i] * c_exps[i] + 4.0 * g_xyy_yz_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_yz_0_yy[i] = -2.0 * g_x_yz_z_yy[i] * c_exps[i] + 4.0 * g_xyy_yz_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_yz_0_yz[i] = -2.0 * g_x_yz_z_yz[i] * c_exps[i] + 4.0 * g_xyy_yz_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_yz_0_zz[i] = -2.0 * g_x_yz_z_zz[i] * c_exps[i] + 4.0 * g_xyy_yz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1146-1152)

    #pragma omp simd aligned(g_x_zz_z_xx, g_x_zz_z_xy, g_x_zz_z_xz, g_x_zz_z_yy, g_x_zz_z_yz, g_x_zz_z_zz, g_xyy_zz_z_xx, g_xyy_zz_z_xy, g_xyy_zz_z_xz, g_xyy_zz_z_yy, g_xyy_zz_z_yz, g_xyy_zz_z_zz, g_y_0_z_0_xy_zz_0_xx, g_y_0_z_0_xy_zz_0_xy, g_y_0_z_0_xy_zz_0_xz, g_y_0_z_0_xy_zz_0_yy, g_y_0_z_0_xy_zz_0_yz, g_y_0_z_0_xy_zz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_xy_zz_0_xx[i] = -2.0 * g_x_zz_z_xx[i] * c_exps[i] + 4.0 * g_xyy_zz_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_zz_0_xy[i] = -2.0 * g_x_zz_z_xy[i] * c_exps[i] + 4.0 * g_xyy_zz_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_zz_0_xz[i] = -2.0 * g_x_zz_z_xz[i] * c_exps[i] + 4.0 * g_xyy_zz_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_zz_0_yy[i] = -2.0 * g_x_zz_z_yy[i] * c_exps[i] + 4.0 * g_xyy_zz_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_zz_0_yz[i] = -2.0 * g_x_zz_z_yz[i] * c_exps[i] + 4.0 * g_xyy_zz_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_zz_0_zz[i] = -2.0 * g_x_zz_z_zz[i] * c_exps[i] + 4.0 * g_xyy_zz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1152-1158)

    #pragma omp simd aligned(g_xyz_xx_z_xx, g_xyz_xx_z_xy, g_xyz_xx_z_xz, g_xyz_xx_z_yy, g_xyz_xx_z_yz, g_xyz_xx_z_zz, g_y_0_z_0_xz_xx_0_xx, g_y_0_z_0_xz_xx_0_xy, g_y_0_z_0_xz_xx_0_xz, g_y_0_z_0_xz_xx_0_yy, g_y_0_z_0_xz_xx_0_yz, g_y_0_z_0_xz_xx_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_xz_xx_0_xx[i] = 4.0 * g_xyz_xx_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_xx_0_xy[i] = 4.0 * g_xyz_xx_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_xx_0_xz[i] = 4.0 * g_xyz_xx_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_xx_0_yy[i] = 4.0 * g_xyz_xx_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_xx_0_yz[i] = 4.0 * g_xyz_xx_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_xx_0_zz[i] = 4.0 * g_xyz_xx_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1158-1164)

    #pragma omp simd aligned(g_xyz_xy_z_xx, g_xyz_xy_z_xy, g_xyz_xy_z_xz, g_xyz_xy_z_yy, g_xyz_xy_z_yz, g_xyz_xy_z_zz, g_y_0_z_0_xz_xy_0_xx, g_y_0_z_0_xz_xy_0_xy, g_y_0_z_0_xz_xy_0_xz, g_y_0_z_0_xz_xy_0_yy, g_y_0_z_0_xz_xy_0_yz, g_y_0_z_0_xz_xy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_xz_xy_0_xx[i] = 4.0 * g_xyz_xy_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_xy_0_xy[i] = 4.0 * g_xyz_xy_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_xy_0_xz[i] = 4.0 * g_xyz_xy_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_xy_0_yy[i] = 4.0 * g_xyz_xy_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_xy_0_yz[i] = 4.0 * g_xyz_xy_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_xy_0_zz[i] = 4.0 * g_xyz_xy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1164-1170)

    #pragma omp simd aligned(g_xyz_xz_z_xx, g_xyz_xz_z_xy, g_xyz_xz_z_xz, g_xyz_xz_z_yy, g_xyz_xz_z_yz, g_xyz_xz_z_zz, g_y_0_z_0_xz_xz_0_xx, g_y_0_z_0_xz_xz_0_xy, g_y_0_z_0_xz_xz_0_xz, g_y_0_z_0_xz_xz_0_yy, g_y_0_z_0_xz_xz_0_yz, g_y_0_z_0_xz_xz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_xz_xz_0_xx[i] = 4.0 * g_xyz_xz_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_xz_0_xy[i] = 4.0 * g_xyz_xz_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_xz_0_xz[i] = 4.0 * g_xyz_xz_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_xz_0_yy[i] = 4.0 * g_xyz_xz_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_xz_0_yz[i] = 4.0 * g_xyz_xz_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_xz_0_zz[i] = 4.0 * g_xyz_xz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1170-1176)

    #pragma omp simd aligned(g_xyz_yy_z_xx, g_xyz_yy_z_xy, g_xyz_yy_z_xz, g_xyz_yy_z_yy, g_xyz_yy_z_yz, g_xyz_yy_z_zz, g_y_0_z_0_xz_yy_0_xx, g_y_0_z_0_xz_yy_0_xy, g_y_0_z_0_xz_yy_0_xz, g_y_0_z_0_xz_yy_0_yy, g_y_0_z_0_xz_yy_0_yz, g_y_0_z_0_xz_yy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_xz_yy_0_xx[i] = 4.0 * g_xyz_yy_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_yy_0_xy[i] = 4.0 * g_xyz_yy_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_yy_0_xz[i] = 4.0 * g_xyz_yy_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_yy_0_yy[i] = 4.0 * g_xyz_yy_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_yy_0_yz[i] = 4.0 * g_xyz_yy_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_yy_0_zz[i] = 4.0 * g_xyz_yy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1176-1182)

    #pragma omp simd aligned(g_xyz_yz_z_xx, g_xyz_yz_z_xy, g_xyz_yz_z_xz, g_xyz_yz_z_yy, g_xyz_yz_z_yz, g_xyz_yz_z_zz, g_y_0_z_0_xz_yz_0_xx, g_y_0_z_0_xz_yz_0_xy, g_y_0_z_0_xz_yz_0_xz, g_y_0_z_0_xz_yz_0_yy, g_y_0_z_0_xz_yz_0_yz, g_y_0_z_0_xz_yz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_xz_yz_0_xx[i] = 4.0 * g_xyz_yz_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_yz_0_xy[i] = 4.0 * g_xyz_yz_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_yz_0_xz[i] = 4.0 * g_xyz_yz_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_yz_0_yy[i] = 4.0 * g_xyz_yz_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_yz_0_yz[i] = 4.0 * g_xyz_yz_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_yz_0_zz[i] = 4.0 * g_xyz_yz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1182-1188)

    #pragma omp simd aligned(g_xyz_zz_z_xx, g_xyz_zz_z_xy, g_xyz_zz_z_xz, g_xyz_zz_z_yy, g_xyz_zz_z_yz, g_xyz_zz_z_zz, g_y_0_z_0_xz_zz_0_xx, g_y_0_z_0_xz_zz_0_xy, g_y_0_z_0_xz_zz_0_xz, g_y_0_z_0_xz_zz_0_yy, g_y_0_z_0_xz_zz_0_yz, g_y_0_z_0_xz_zz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_xz_zz_0_xx[i] = 4.0 * g_xyz_zz_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_zz_0_xy[i] = 4.0 * g_xyz_zz_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_zz_0_xz[i] = 4.0 * g_xyz_zz_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_zz_0_yy[i] = 4.0 * g_xyz_zz_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_zz_0_yz[i] = 4.0 * g_xyz_zz_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_zz_0_zz[i] = 4.0 * g_xyz_zz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1188-1194)

    #pragma omp simd aligned(g_y_0_z_0_yy_xx_0_xx, g_y_0_z_0_yy_xx_0_xy, g_y_0_z_0_yy_xx_0_xz, g_y_0_z_0_yy_xx_0_yy, g_y_0_z_0_yy_xx_0_yz, g_y_0_z_0_yy_xx_0_zz, g_y_xx_z_xx, g_y_xx_z_xy, g_y_xx_z_xz, g_y_xx_z_yy, g_y_xx_z_yz, g_y_xx_z_zz, g_yyy_xx_z_xx, g_yyy_xx_z_xy, g_yyy_xx_z_xz, g_yyy_xx_z_yy, g_yyy_xx_z_yz, g_yyy_xx_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_yy_xx_0_xx[i] = -4.0 * g_y_xx_z_xx[i] * c_exps[i] + 4.0 * g_yyy_xx_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_xx_0_xy[i] = -4.0 * g_y_xx_z_xy[i] * c_exps[i] + 4.0 * g_yyy_xx_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_xx_0_xz[i] = -4.0 * g_y_xx_z_xz[i] * c_exps[i] + 4.0 * g_yyy_xx_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_xx_0_yy[i] = -4.0 * g_y_xx_z_yy[i] * c_exps[i] + 4.0 * g_yyy_xx_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_xx_0_yz[i] = -4.0 * g_y_xx_z_yz[i] * c_exps[i] + 4.0 * g_yyy_xx_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_xx_0_zz[i] = -4.0 * g_y_xx_z_zz[i] * c_exps[i] + 4.0 * g_yyy_xx_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1194-1200)

    #pragma omp simd aligned(g_y_0_z_0_yy_xy_0_xx, g_y_0_z_0_yy_xy_0_xy, g_y_0_z_0_yy_xy_0_xz, g_y_0_z_0_yy_xy_0_yy, g_y_0_z_0_yy_xy_0_yz, g_y_0_z_0_yy_xy_0_zz, g_y_xy_z_xx, g_y_xy_z_xy, g_y_xy_z_xz, g_y_xy_z_yy, g_y_xy_z_yz, g_y_xy_z_zz, g_yyy_xy_z_xx, g_yyy_xy_z_xy, g_yyy_xy_z_xz, g_yyy_xy_z_yy, g_yyy_xy_z_yz, g_yyy_xy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_yy_xy_0_xx[i] = -4.0 * g_y_xy_z_xx[i] * c_exps[i] + 4.0 * g_yyy_xy_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_xy_0_xy[i] = -4.0 * g_y_xy_z_xy[i] * c_exps[i] + 4.0 * g_yyy_xy_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_xy_0_xz[i] = -4.0 * g_y_xy_z_xz[i] * c_exps[i] + 4.0 * g_yyy_xy_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_xy_0_yy[i] = -4.0 * g_y_xy_z_yy[i] * c_exps[i] + 4.0 * g_yyy_xy_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_xy_0_yz[i] = -4.0 * g_y_xy_z_yz[i] * c_exps[i] + 4.0 * g_yyy_xy_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_xy_0_zz[i] = -4.0 * g_y_xy_z_zz[i] * c_exps[i] + 4.0 * g_yyy_xy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1200-1206)

    #pragma omp simd aligned(g_y_0_z_0_yy_xz_0_xx, g_y_0_z_0_yy_xz_0_xy, g_y_0_z_0_yy_xz_0_xz, g_y_0_z_0_yy_xz_0_yy, g_y_0_z_0_yy_xz_0_yz, g_y_0_z_0_yy_xz_0_zz, g_y_xz_z_xx, g_y_xz_z_xy, g_y_xz_z_xz, g_y_xz_z_yy, g_y_xz_z_yz, g_y_xz_z_zz, g_yyy_xz_z_xx, g_yyy_xz_z_xy, g_yyy_xz_z_xz, g_yyy_xz_z_yy, g_yyy_xz_z_yz, g_yyy_xz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_yy_xz_0_xx[i] = -4.0 * g_y_xz_z_xx[i] * c_exps[i] + 4.0 * g_yyy_xz_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_xz_0_xy[i] = -4.0 * g_y_xz_z_xy[i] * c_exps[i] + 4.0 * g_yyy_xz_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_xz_0_xz[i] = -4.0 * g_y_xz_z_xz[i] * c_exps[i] + 4.0 * g_yyy_xz_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_xz_0_yy[i] = -4.0 * g_y_xz_z_yy[i] * c_exps[i] + 4.0 * g_yyy_xz_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_xz_0_yz[i] = -4.0 * g_y_xz_z_yz[i] * c_exps[i] + 4.0 * g_yyy_xz_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_xz_0_zz[i] = -4.0 * g_y_xz_z_zz[i] * c_exps[i] + 4.0 * g_yyy_xz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1206-1212)

    #pragma omp simd aligned(g_y_0_z_0_yy_yy_0_xx, g_y_0_z_0_yy_yy_0_xy, g_y_0_z_0_yy_yy_0_xz, g_y_0_z_0_yy_yy_0_yy, g_y_0_z_0_yy_yy_0_yz, g_y_0_z_0_yy_yy_0_zz, g_y_yy_z_xx, g_y_yy_z_xy, g_y_yy_z_xz, g_y_yy_z_yy, g_y_yy_z_yz, g_y_yy_z_zz, g_yyy_yy_z_xx, g_yyy_yy_z_xy, g_yyy_yy_z_xz, g_yyy_yy_z_yy, g_yyy_yy_z_yz, g_yyy_yy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_yy_yy_0_xx[i] = -4.0 * g_y_yy_z_xx[i] * c_exps[i] + 4.0 * g_yyy_yy_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_yy_0_xy[i] = -4.0 * g_y_yy_z_xy[i] * c_exps[i] + 4.0 * g_yyy_yy_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_yy_0_xz[i] = -4.0 * g_y_yy_z_xz[i] * c_exps[i] + 4.0 * g_yyy_yy_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_yy_0_yy[i] = -4.0 * g_y_yy_z_yy[i] * c_exps[i] + 4.0 * g_yyy_yy_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_yy_0_yz[i] = -4.0 * g_y_yy_z_yz[i] * c_exps[i] + 4.0 * g_yyy_yy_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_yy_0_zz[i] = -4.0 * g_y_yy_z_zz[i] * c_exps[i] + 4.0 * g_yyy_yy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1212-1218)

    #pragma omp simd aligned(g_y_0_z_0_yy_yz_0_xx, g_y_0_z_0_yy_yz_0_xy, g_y_0_z_0_yy_yz_0_xz, g_y_0_z_0_yy_yz_0_yy, g_y_0_z_0_yy_yz_0_yz, g_y_0_z_0_yy_yz_0_zz, g_y_yz_z_xx, g_y_yz_z_xy, g_y_yz_z_xz, g_y_yz_z_yy, g_y_yz_z_yz, g_y_yz_z_zz, g_yyy_yz_z_xx, g_yyy_yz_z_xy, g_yyy_yz_z_xz, g_yyy_yz_z_yy, g_yyy_yz_z_yz, g_yyy_yz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_yy_yz_0_xx[i] = -4.0 * g_y_yz_z_xx[i] * c_exps[i] + 4.0 * g_yyy_yz_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_yz_0_xy[i] = -4.0 * g_y_yz_z_xy[i] * c_exps[i] + 4.0 * g_yyy_yz_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_yz_0_xz[i] = -4.0 * g_y_yz_z_xz[i] * c_exps[i] + 4.0 * g_yyy_yz_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_yz_0_yy[i] = -4.0 * g_y_yz_z_yy[i] * c_exps[i] + 4.0 * g_yyy_yz_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_yz_0_yz[i] = -4.0 * g_y_yz_z_yz[i] * c_exps[i] + 4.0 * g_yyy_yz_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_yz_0_zz[i] = -4.0 * g_y_yz_z_zz[i] * c_exps[i] + 4.0 * g_yyy_yz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1218-1224)

    #pragma omp simd aligned(g_y_0_z_0_yy_zz_0_xx, g_y_0_z_0_yy_zz_0_xy, g_y_0_z_0_yy_zz_0_xz, g_y_0_z_0_yy_zz_0_yy, g_y_0_z_0_yy_zz_0_yz, g_y_0_z_0_yy_zz_0_zz, g_y_zz_z_xx, g_y_zz_z_xy, g_y_zz_z_xz, g_y_zz_z_yy, g_y_zz_z_yz, g_y_zz_z_zz, g_yyy_zz_z_xx, g_yyy_zz_z_xy, g_yyy_zz_z_xz, g_yyy_zz_z_yy, g_yyy_zz_z_yz, g_yyy_zz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_yy_zz_0_xx[i] = -4.0 * g_y_zz_z_xx[i] * c_exps[i] + 4.0 * g_yyy_zz_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_zz_0_xy[i] = -4.0 * g_y_zz_z_xy[i] * c_exps[i] + 4.0 * g_yyy_zz_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_zz_0_xz[i] = -4.0 * g_y_zz_z_xz[i] * c_exps[i] + 4.0 * g_yyy_zz_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_zz_0_yy[i] = -4.0 * g_y_zz_z_yy[i] * c_exps[i] + 4.0 * g_yyy_zz_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_zz_0_yz[i] = -4.0 * g_y_zz_z_yz[i] * c_exps[i] + 4.0 * g_yyy_zz_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_zz_0_zz[i] = -4.0 * g_y_zz_z_zz[i] * c_exps[i] + 4.0 * g_yyy_zz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1224-1230)

    #pragma omp simd aligned(g_y_0_z_0_yz_xx_0_xx, g_y_0_z_0_yz_xx_0_xy, g_y_0_z_0_yz_xx_0_xz, g_y_0_z_0_yz_xx_0_yy, g_y_0_z_0_yz_xx_0_yz, g_y_0_z_0_yz_xx_0_zz, g_yyz_xx_z_xx, g_yyz_xx_z_xy, g_yyz_xx_z_xz, g_yyz_xx_z_yy, g_yyz_xx_z_yz, g_yyz_xx_z_zz, g_z_xx_z_xx, g_z_xx_z_xy, g_z_xx_z_xz, g_z_xx_z_yy, g_z_xx_z_yz, g_z_xx_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_yz_xx_0_xx[i] = -2.0 * g_z_xx_z_xx[i] * c_exps[i] + 4.0 * g_yyz_xx_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_xx_0_xy[i] = -2.0 * g_z_xx_z_xy[i] * c_exps[i] + 4.0 * g_yyz_xx_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_xx_0_xz[i] = -2.0 * g_z_xx_z_xz[i] * c_exps[i] + 4.0 * g_yyz_xx_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_xx_0_yy[i] = -2.0 * g_z_xx_z_yy[i] * c_exps[i] + 4.0 * g_yyz_xx_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_xx_0_yz[i] = -2.0 * g_z_xx_z_yz[i] * c_exps[i] + 4.0 * g_yyz_xx_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_xx_0_zz[i] = -2.0 * g_z_xx_z_zz[i] * c_exps[i] + 4.0 * g_yyz_xx_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1230-1236)

    #pragma omp simd aligned(g_y_0_z_0_yz_xy_0_xx, g_y_0_z_0_yz_xy_0_xy, g_y_0_z_0_yz_xy_0_xz, g_y_0_z_0_yz_xy_0_yy, g_y_0_z_0_yz_xy_0_yz, g_y_0_z_0_yz_xy_0_zz, g_yyz_xy_z_xx, g_yyz_xy_z_xy, g_yyz_xy_z_xz, g_yyz_xy_z_yy, g_yyz_xy_z_yz, g_yyz_xy_z_zz, g_z_xy_z_xx, g_z_xy_z_xy, g_z_xy_z_xz, g_z_xy_z_yy, g_z_xy_z_yz, g_z_xy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_yz_xy_0_xx[i] = -2.0 * g_z_xy_z_xx[i] * c_exps[i] + 4.0 * g_yyz_xy_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_xy_0_xy[i] = -2.0 * g_z_xy_z_xy[i] * c_exps[i] + 4.0 * g_yyz_xy_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_xy_0_xz[i] = -2.0 * g_z_xy_z_xz[i] * c_exps[i] + 4.0 * g_yyz_xy_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_xy_0_yy[i] = -2.0 * g_z_xy_z_yy[i] * c_exps[i] + 4.0 * g_yyz_xy_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_xy_0_yz[i] = -2.0 * g_z_xy_z_yz[i] * c_exps[i] + 4.0 * g_yyz_xy_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_xy_0_zz[i] = -2.0 * g_z_xy_z_zz[i] * c_exps[i] + 4.0 * g_yyz_xy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1236-1242)

    #pragma omp simd aligned(g_y_0_z_0_yz_xz_0_xx, g_y_0_z_0_yz_xz_0_xy, g_y_0_z_0_yz_xz_0_xz, g_y_0_z_0_yz_xz_0_yy, g_y_0_z_0_yz_xz_0_yz, g_y_0_z_0_yz_xz_0_zz, g_yyz_xz_z_xx, g_yyz_xz_z_xy, g_yyz_xz_z_xz, g_yyz_xz_z_yy, g_yyz_xz_z_yz, g_yyz_xz_z_zz, g_z_xz_z_xx, g_z_xz_z_xy, g_z_xz_z_xz, g_z_xz_z_yy, g_z_xz_z_yz, g_z_xz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_yz_xz_0_xx[i] = -2.0 * g_z_xz_z_xx[i] * c_exps[i] + 4.0 * g_yyz_xz_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_xz_0_xy[i] = -2.0 * g_z_xz_z_xy[i] * c_exps[i] + 4.0 * g_yyz_xz_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_xz_0_xz[i] = -2.0 * g_z_xz_z_xz[i] * c_exps[i] + 4.0 * g_yyz_xz_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_xz_0_yy[i] = -2.0 * g_z_xz_z_yy[i] * c_exps[i] + 4.0 * g_yyz_xz_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_xz_0_yz[i] = -2.0 * g_z_xz_z_yz[i] * c_exps[i] + 4.0 * g_yyz_xz_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_xz_0_zz[i] = -2.0 * g_z_xz_z_zz[i] * c_exps[i] + 4.0 * g_yyz_xz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1242-1248)

    #pragma omp simd aligned(g_y_0_z_0_yz_yy_0_xx, g_y_0_z_0_yz_yy_0_xy, g_y_0_z_0_yz_yy_0_xz, g_y_0_z_0_yz_yy_0_yy, g_y_0_z_0_yz_yy_0_yz, g_y_0_z_0_yz_yy_0_zz, g_yyz_yy_z_xx, g_yyz_yy_z_xy, g_yyz_yy_z_xz, g_yyz_yy_z_yy, g_yyz_yy_z_yz, g_yyz_yy_z_zz, g_z_yy_z_xx, g_z_yy_z_xy, g_z_yy_z_xz, g_z_yy_z_yy, g_z_yy_z_yz, g_z_yy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_yz_yy_0_xx[i] = -2.0 * g_z_yy_z_xx[i] * c_exps[i] + 4.0 * g_yyz_yy_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_yy_0_xy[i] = -2.0 * g_z_yy_z_xy[i] * c_exps[i] + 4.0 * g_yyz_yy_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_yy_0_xz[i] = -2.0 * g_z_yy_z_xz[i] * c_exps[i] + 4.0 * g_yyz_yy_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_yy_0_yy[i] = -2.0 * g_z_yy_z_yy[i] * c_exps[i] + 4.0 * g_yyz_yy_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_yy_0_yz[i] = -2.0 * g_z_yy_z_yz[i] * c_exps[i] + 4.0 * g_yyz_yy_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_yy_0_zz[i] = -2.0 * g_z_yy_z_zz[i] * c_exps[i] + 4.0 * g_yyz_yy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1248-1254)

    #pragma omp simd aligned(g_y_0_z_0_yz_yz_0_xx, g_y_0_z_0_yz_yz_0_xy, g_y_0_z_0_yz_yz_0_xz, g_y_0_z_0_yz_yz_0_yy, g_y_0_z_0_yz_yz_0_yz, g_y_0_z_0_yz_yz_0_zz, g_yyz_yz_z_xx, g_yyz_yz_z_xy, g_yyz_yz_z_xz, g_yyz_yz_z_yy, g_yyz_yz_z_yz, g_yyz_yz_z_zz, g_z_yz_z_xx, g_z_yz_z_xy, g_z_yz_z_xz, g_z_yz_z_yy, g_z_yz_z_yz, g_z_yz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_yz_yz_0_xx[i] = -2.0 * g_z_yz_z_xx[i] * c_exps[i] + 4.0 * g_yyz_yz_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_yz_0_xy[i] = -2.0 * g_z_yz_z_xy[i] * c_exps[i] + 4.0 * g_yyz_yz_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_yz_0_xz[i] = -2.0 * g_z_yz_z_xz[i] * c_exps[i] + 4.0 * g_yyz_yz_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_yz_0_yy[i] = -2.0 * g_z_yz_z_yy[i] * c_exps[i] + 4.0 * g_yyz_yz_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_yz_0_yz[i] = -2.0 * g_z_yz_z_yz[i] * c_exps[i] + 4.0 * g_yyz_yz_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_yz_0_zz[i] = -2.0 * g_z_yz_z_zz[i] * c_exps[i] + 4.0 * g_yyz_yz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1254-1260)

    #pragma omp simd aligned(g_y_0_z_0_yz_zz_0_xx, g_y_0_z_0_yz_zz_0_xy, g_y_0_z_0_yz_zz_0_xz, g_y_0_z_0_yz_zz_0_yy, g_y_0_z_0_yz_zz_0_yz, g_y_0_z_0_yz_zz_0_zz, g_yyz_zz_z_xx, g_yyz_zz_z_xy, g_yyz_zz_z_xz, g_yyz_zz_z_yy, g_yyz_zz_z_yz, g_yyz_zz_z_zz, g_z_zz_z_xx, g_z_zz_z_xy, g_z_zz_z_xz, g_z_zz_z_yy, g_z_zz_z_yz, g_z_zz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_yz_zz_0_xx[i] = -2.0 * g_z_zz_z_xx[i] * c_exps[i] + 4.0 * g_yyz_zz_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_zz_0_xy[i] = -2.0 * g_z_zz_z_xy[i] * c_exps[i] + 4.0 * g_yyz_zz_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_zz_0_xz[i] = -2.0 * g_z_zz_z_xz[i] * c_exps[i] + 4.0 * g_yyz_zz_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_zz_0_yy[i] = -2.0 * g_z_zz_z_yy[i] * c_exps[i] + 4.0 * g_yyz_zz_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_zz_0_yz[i] = -2.0 * g_z_zz_z_yz[i] * c_exps[i] + 4.0 * g_yyz_zz_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_zz_0_zz[i] = -2.0 * g_z_zz_z_zz[i] * c_exps[i] + 4.0 * g_yyz_zz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1260-1266)

    #pragma omp simd aligned(g_y_0_z_0_zz_xx_0_xx, g_y_0_z_0_zz_xx_0_xy, g_y_0_z_0_zz_xx_0_xz, g_y_0_z_0_zz_xx_0_yy, g_y_0_z_0_zz_xx_0_yz, g_y_0_z_0_zz_xx_0_zz, g_yzz_xx_z_xx, g_yzz_xx_z_xy, g_yzz_xx_z_xz, g_yzz_xx_z_yy, g_yzz_xx_z_yz, g_yzz_xx_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_zz_xx_0_xx[i] = 4.0 * g_yzz_xx_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_xx_0_xy[i] = 4.0 * g_yzz_xx_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_xx_0_xz[i] = 4.0 * g_yzz_xx_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_xx_0_yy[i] = 4.0 * g_yzz_xx_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_xx_0_yz[i] = 4.0 * g_yzz_xx_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_xx_0_zz[i] = 4.0 * g_yzz_xx_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1266-1272)

    #pragma omp simd aligned(g_y_0_z_0_zz_xy_0_xx, g_y_0_z_0_zz_xy_0_xy, g_y_0_z_0_zz_xy_0_xz, g_y_0_z_0_zz_xy_0_yy, g_y_0_z_0_zz_xy_0_yz, g_y_0_z_0_zz_xy_0_zz, g_yzz_xy_z_xx, g_yzz_xy_z_xy, g_yzz_xy_z_xz, g_yzz_xy_z_yy, g_yzz_xy_z_yz, g_yzz_xy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_zz_xy_0_xx[i] = 4.0 * g_yzz_xy_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_xy_0_xy[i] = 4.0 * g_yzz_xy_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_xy_0_xz[i] = 4.0 * g_yzz_xy_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_xy_0_yy[i] = 4.0 * g_yzz_xy_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_xy_0_yz[i] = 4.0 * g_yzz_xy_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_xy_0_zz[i] = 4.0 * g_yzz_xy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1272-1278)

    #pragma omp simd aligned(g_y_0_z_0_zz_xz_0_xx, g_y_0_z_0_zz_xz_0_xy, g_y_0_z_0_zz_xz_0_xz, g_y_0_z_0_zz_xz_0_yy, g_y_0_z_0_zz_xz_0_yz, g_y_0_z_0_zz_xz_0_zz, g_yzz_xz_z_xx, g_yzz_xz_z_xy, g_yzz_xz_z_xz, g_yzz_xz_z_yy, g_yzz_xz_z_yz, g_yzz_xz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_zz_xz_0_xx[i] = 4.0 * g_yzz_xz_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_xz_0_xy[i] = 4.0 * g_yzz_xz_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_xz_0_xz[i] = 4.0 * g_yzz_xz_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_xz_0_yy[i] = 4.0 * g_yzz_xz_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_xz_0_yz[i] = 4.0 * g_yzz_xz_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_xz_0_zz[i] = 4.0 * g_yzz_xz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1278-1284)

    #pragma omp simd aligned(g_y_0_z_0_zz_yy_0_xx, g_y_0_z_0_zz_yy_0_xy, g_y_0_z_0_zz_yy_0_xz, g_y_0_z_0_zz_yy_0_yy, g_y_0_z_0_zz_yy_0_yz, g_y_0_z_0_zz_yy_0_zz, g_yzz_yy_z_xx, g_yzz_yy_z_xy, g_yzz_yy_z_xz, g_yzz_yy_z_yy, g_yzz_yy_z_yz, g_yzz_yy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_zz_yy_0_xx[i] = 4.0 * g_yzz_yy_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_yy_0_xy[i] = 4.0 * g_yzz_yy_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_yy_0_xz[i] = 4.0 * g_yzz_yy_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_yy_0_yy[i] = 4.0 * g_yzz_yy_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_yy_0_yz[i] = 4.0 * g_yzz_yy_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_yy_0_zz[i] = 4.0 * g_yzz_yy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1284-1290)

    #pragma omp simd aligned(g_y_0_z_0_zz_yz_0_xx, g_y_0_z_0_zz_yz_0_xy, g_y_0_z_0_zz_yz_0_xz, g_y_0_z_0_zz_yz_0_yy, g_y_0_z_0_zz_yz_0_yz, g_y_0_z_0_zz_yz_0_zz, g_yzz_yz_z_xx, g_yzz_yz_z_xy, g_yzz_yz_z_xz, g_yzz_yz_z_yy, g_yzz_yz_z_yz, g_yzz_yz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_zz_yz_0_xx[i] = 4.0 * g_yzz_yz_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_yz_0_xy[i] = 4.0 * g_yzz_yz_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_yz_0_xz[i] = 4.0 * g_yzz_yz_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_yz_0_yy[i] = 4.0 * g_yzz_yz_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_yz_0_yz[i] = 4.0 * g_yzz_yz_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_yz_0_zz[i] = 4.0 * g_yzz_yz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1290-1296)

    #pragma omp simd aligned(g_y_0_z_0_zz_zz_0_xx, g_y_0_z_0_zz_zz_0_xy, g_y_0_z_0_zz_zz_0_xz, g_y_0_z_0_zz_zz_0_yy, g_y_0_z_0_zz_zz_0_yz, g_y_0_z_0_zz_zz_0_zz, g_yzz_zz_z_xx, g_yzz_zz_z_xy, g_yzz_zz_z_xz, g_yzz_zz_z_yy, g_yzz_zz_z_yz, g_yzz_zz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_zz_zz_0_xx[i] = 4.0 * g_yzz_zz_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_zz_0_xy[i] = 4.0 * g_yzz_zz_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_zz_0_xz[i] = 4.0 * g_yzz_zz_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_zz_0_yy[i] = 4.0 * g_yzz_zz_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_zz_0_yz[i] = 4.0 * g_yzz_zz_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_zz_0_zz[i] = 4.0 * g_yzz_zz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1296-1302)

    #pragma omp simd aligned(g_xxz_xx_x_xx, g_xxz_xx_x_xy, g_xxz_xx_x_xz, g_xxz_xx_x_yy, g_xxz_xx_x_yz, g_xxz_xx_x_zz, g_z_0_x_0_xx_xx_0_xx, g_z_0_x_0_xx_xx_0_xy, g_z_0_x_0_xx_xx_0_xz, g_z_0_x_0_xx_xx_0_yy, g_z_0_x_0_xx_xx_0_yz, g_z_0_x_0_xx_xx_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_xx_xx_0_xx[i] = 4.0 * g_xxz_xx_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_xx_0_xy[i] = 4.0 * g_xxz_xx_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_xx_0_xz[i] = 4.0 * g_xxz_xx_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_xx_0_yy[i] = 4.0 * g_xxz_xx_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_xx_0_yz[i] = 4.0 * g_xxz_xx_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_xx_0_zz[i] = 4.0 * g_xxz_xx_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1302-1308)

    #pragma omp simd aligned(g_xxz_xy_x_xx, g_xxz_xy_x_xy, g_xxz_xy_x_xz, g_xxz_xy_x_yy, g_xxz_xy_x_yz, g_xxz_xy_x_zz, g_z_0_x_0_xx_xy_0_xx, g_z_0_x_0_xx_xy_0_xy, g_z_0_x_0_xx_xy_0_xz, g_z_0_x_0_xx_xy_0_yy, g_z_0_x_0_xx_xy_0_yz, g_z_0_x_0_xx_xy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_xx_xy_0_xx[i] = 4.0 * g_xxz_xy_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_xy_0_xy[i] = 4.0 * g_xxz_xy_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_xy_0_xz[i] = 4.0 * g_xxz_xy_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_xy_0_yy[i] = 4.0 * g_xxz_xy_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_xy_0_yz[i] = 4.0 * g_xxz_xy_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_xy_0_zz[i] = 4.0 * g_xxz_xy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1308-1314)

    #pragma omp simd aligned(g_xxz_xz_x_xx, g_xxz_xz_x_xy, g_xxz_xz_x_xz, g_xxz_xz_x_yy, g_xxz_xz_x_yz, g_xxz_xz_x_zz, g_z_0_x_0_xx_xz_0_xx, g_z_0_x_0_xx_xz_0_xy, g_z_0_x_0_xx_xz_0_xz, g_z_0_x_0_xx_xz_0_yy, g_z_0_x_0_xx_xz_0_yz, g_z_0_x_0_xx_xz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_xx_xz_0_xx[i] = 4.0 * g_xxz_xz_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_xz_0_xy[i] = 4.0 * g_xxz_xz_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_xz_0_xz[i] = 4.0 * g_xxz_xz_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_xz_0_yy[i] = 4.0 * g_xxz_xz_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_xz_0_yz[i] = 4.0 * g_xxz_xz_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_xz_0_zz[i] = 4.0 * g_xxz_xz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1314-1320)

    #pragma omp simd aligned(g_xxz_yy_x_xx, g_xxz_yy_x_xy, g_xxz_yy_x_xz, g_xxz_yy_x_yy, g_xxz_yy_x_yz, g_xxz_yy_x_zz, g_z_0_x_0_xx_yy_0_xx, g_z_0_x_0_xx_yy_0_xy, g_z_0_x_0_xx_yy_0_xz, g_z_0_x_0_xx_yy_0_yy, g_z_0_x_0_xx_yy_0_yz, g_z_0_x_0_xx_yy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_xx_yy_0_xx[i] = 4.0 * g_xxz_yy_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_yy_0_xy[i] = 4.0 * g_xxz_yy_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_yy_0_xz[i] = 4.0 * g_xxz_yy_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_yy_0_yy[i] = 4.0 * g_xxz_yy_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_yy_0_yz[i] = 4.0 * g_xxz_yy_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_yy_0_zz[i] = 4.0 * g_xxz_yy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1320-1326)

    #pragma omp simd aligned(g_xxz_yz_x_xx, g_xxz_yz_x_xy, g_xxz_yz_x_xz, g_xxz_yz_x_yy, g_xxz_yz_x_yz, g_xxz_yz_x_zz, g_z_0_x_0_xx_yz_0_xx, g_z_0_x_0_xx_yz_0_xy, g_z_0_x_0_xx_yz_0_xz, g_z_0_x_0_xx_yz_0_yy, g_z_0_x_0_xx_yz_0_yz, g_z_0_x_0_xx_yz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_xx_yz_0_xx[i] = 4.0 * g_xxz_yz_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_yz_0_xy[i] = 4.0 * g_xxz_yz_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_yz_0_xz[i] = 4.0 * g_xxz_yz_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_yz_0_yy[i] = 4.0 * g_xxz_yz_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_yz_0_yz[i] = 4.0 * g_xxz_yz_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_yz_0_zz[i] = 4.0 * g_xxz_yz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1326-1332)

    #pragma omp simd aligned(g_xxz_zz_x_xx, g_xxz_zz_x_xy, g_xxz_zz_x_xz, g_xxz_zz_x_yy, g_xxz_zz_x_yz, g_xxz_zz_x_zz, g_z_0_x_0_xx_zz_0_xx, g_z_0_x_0_xx_zz_0_xy, g_z_0_x_0_xx_zz_0_xz, g_z_0_x_0_xx_zz_0_yy, g_z_0_x_0_xx_zz_0_yz, g_z_0_x_0_xx_zz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_xx_zz_0_xx[i] = 4.0 * g_xxz_zz_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_zz_0_xy[i] = 4.0 * g_xxz_zz_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_zz_0_xz[i] = 4.0 * g_xxz_zz_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_zz_0_yy[i] = 4.0 * g_xxz_zz_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_zz_0_yz[i] = 4.0 * g_xxz_zz_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_zz_0_zz[i] = 4.0 * g_xxz_zz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1332-1338)

    #pragma omp simd aligned(g_xyz_xx_x_xx, g_xyz_xx_x_xy, g_xyz_xx_x_xz, g_xyz_xx_x_yy, g_xyz_xx_x_yz, g_xyz_xx_x_zz, g_z_0_x_0_xy_xx_0_xx, g_z_0_x_0_xy_xx_0_xy, g_z_0_x_0_xy_xx_0_xz, g_z_0_x_0_xy_xx_0_yy, g_z_0_x_0_xy_xx_0_yz, g_z_0_x_0_xy_xx_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_xy_xx_0_xx[i] = 4.0 * g_xyz_xx_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_xx_0_xy[i] = 4.0 * g_xyz_xx_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_xx_0_xz[i] = 4.0 * g_xyz_xx_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_xx_0_yy[i] = 4.0 * g_xyz_xx_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_xx_0_yz[i] = 4.0 * g_xyz_xx_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_xx_0_zz[i] = 4.0 * g_xyz_xx_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1338-1344)

    #pragma omp simd aligned(g_xyz_xy_x_xx, g_xyz_xy_x_xy, g_xyz_xy_x_xz, g_xyz_xy_x_yy, g_xyz_xy_x_yz, g_xyz_xy_x_zz, g_z_0_x_0_xy_xy_0_xx, g_z_0_x_0_xy_xy_0_xy, g_z_0_x_0_xy_xy_0_xz, g_z_0_x_0_xy_xy_0_yy, g_z_0_x_0_xy_xy_0_yz, g_z_0_x_0_xy_xy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_xy_xy_0_xx[i] = 4.0 * g_xyz_xy_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_xy_0_xy[i] = 4.0 * g_xyz_xy_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_xy_0_xz[i] = 4.0 * g_xyz_xy_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_xy_0_yy[i] = 4.0 * g_xyz_xy_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_xy_0_yz[i] = 4.0 * g_xyz_xy_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_xy_0_zz[i] = 4.0 * g_xyz_xy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1344-1350)

    #pragma omp simd aligned(g_xyz_xz_x_xx, g_xyz_xz_x_xy, g_xyz_xz_x_xz, g_xyz_xz_x_yy, g_xyz_xz_x_yz, g_xyz_xz_x_zz, g_z_0_x_0_xy_xz_0_xx, g_z_0_x_0_xy_xz_0_xy, g_z_0_x_0_xy_xz_0_xz, g_z_0_x_0_xy_xz_0_yy, g_z_0_x_0_xy_xz_0_yz, g_z_0_x_0_xy_xz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_xy_xz_0_xx[i] = 4.0 * g_xyz_xz_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_xz_0_xy[i] = 4.0 * g_xyz_xz_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_xz_0_xz[i] = 4.0 * g_xyz_xz_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_xz_0_yy[i] = 4.0 * g_xyz_xz_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_xz_0_yz[i] = 4.0 * g_xyz_xz_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_xz_0_zz[i] = 4.0 * g_xyz_xz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1350-1356)

    #pragma omp simd aligned(g_xyz_yy_x_xx, g_xyz_yy_x_xy, g_xyz_yy_x_xz, g_xyz_yy_x_yy, g_xyz_yy_x_yz, g_xyz_yy_x_zz, g_z_0_x_0_xy_yy_0_xx, g_z_0_x_0_xy_yy_0_xy, g_z_0_x_0_xy_yy_0_xz, g_z_0_x_0_xy_yy_0_yy, g_z_0_x_0_xy_yy_0_yz, g_z_0_x_0_xy_yy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_xy_yy_0_xx[i] = 4.0 * g_xyz_yy_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_yy_0_xy[i] = 4.0 * g_xyz_yy_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_yy_0_xz[i] = 4.0 * g_xyz_yy_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_yy_0_yy[i] = 4.0 * g_xyz_yy_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_yy_0_yz[i] = 4.0 * g_xyz_yy_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_yy_0_zz[i] = 4.0 * g_xyz_yy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1356-1362)

    #pragma omp simd aligned(g_xyz_yz_x_xx, g_xyz_yz_x_xy, g_xyz_yz_x_xz, g_xyz_yz_x_yy, g_xyz_yz_x_yz, g_xyz_yz_x_zz, g_z_0_x_0_xy_yz_0_xx, g_z_0_x_0_xy_yz_0_xy, g_z_0_x_0_xy_yz_0_xz, g_z_0_x_0_xy_yz_0_yy, g_z_0_x_0_xy_yz_0_yz, g_z_0_x_0_xy_yz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_xy_yz_0_xx[i] = 4.0 * g_xyz_yz_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_yz_0_xy[i] = 4.0 * g_xyz_yz_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_yz_0_xz[i] = 4.0 * g_xyz_yz_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_yz_0_yy[i] = 4.0 * g_xyz_yz_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_yz_0_yz[i] = 4.0 * g_xyz_yz_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_yz_0_zz[i] = 4.0 * g_xyz_yz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1362-1368)

    #pragma omp simd aligned(g_xyz_zz_x_xx, g_xyz_zz_x_xy, g_xyz_zz_x_xz, g_xyz_zz_x_yy, g_xyz_zz_x_yz, g_xyz_zz_x_zz, g_z_0_x_0_xy_zz_0_xx, g_z_0_x_0_xy_zz_0_xy, g_z_0_x_0_xy_zz_0_xz, g_z_0_x_0_xy_zz_0_yy, g_z_0_x_0_xy_zz_0_yz, g_z_0_x_0_xy_zz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_xy_zz_0_xx[i] = 4.0 * g_xyz_zz_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_zz_0_xy[i] = 4.0 * g_xyz_zz_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_zz_0_xz[i] = 4.0 * g_xyz_zz_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_zz_0_yy[i] = 4.0 * g_xyz_zz_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_zz_0_yz[i] = 4.0 * g_xyz_zz_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_zz_0_zz[i] = 4.0 * g_xyz_zz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1368-1374)

    #pragma omp simd aligned(g_x_xx_x_xx, g_x_xx_x_xy, g_x_xx_x_xz, g_x_xx_x_yy, g_x_xx_x_yz, g_x_xx_x_zz, g_xzz_xx_x_xx, g_xzz_xx_x_xy, g_xzz_xx_x_xz, g_xzz_xx_x_yy, g_xzz_xx_x_yz, g_xzz_xx_x_zz, g_z_0_x_0_xz_xx_0_xx, g_z_0_x_0_xz_xx_0_xy, g_z_0_x_0_xz_xx_0_xz, g_z_0_x_0_xz_xx_0_yy, g_z_0_x_0_xz_xx_0_yz, g_z_0_x_0_xz_xx_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_xz_xx_0_xx[i] = -2.0 * g_x_xx_x_xx[i] * c_exps[i] + 4.0 * g_xzz_xx_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_xx_0_xy[i] = -2.0 * g_x_xx_x_xy[i] * c_exps[i] + 4.0 * g_xzz_xx_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_xx_0_xz[i] = -2.0 * g_x_xx_x_xz[i] * c_exps[i] + 4.0 * g_xzz_xx_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_xx_0_yy[i] = -2.0 * g_x_xx_x_yy[i] * c_exps[i] + 4.0 * g_xzz_xx_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_xx_0_yz[i] = -2.0 * g_x_xx_x_yz[i] * c_exps[i] + 4.0 * g_xzz_xx_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_xx_0_zz[i] = -2.0 * g_x_xx_x_zz[i] * c_exps[i] + 4.0 * g_xzz_xx_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1374-1380)

    #pragma omp simd aligned(g_x_xy_x_xx, g_x_xy_x_xy, g_x_xy_x_xz, g_x_xy_x_yy, g_x_xy_x_yz, g_x_xy_x_zz, g_xzz_xy_x_xx, g_xzz_xy_x_xy, g_xzz_xy_x_xz, g_xzz_xy_x_yy, g_xzz_xy_x_yz, g_xzz_xy_x_zz, g_z_0_x_0_xz_xy_0_xx, g_z_0_x_0_xz_xy_0_xy, g_z_0_x_0_xz_xy_0_xz, g_z_0_x_0_xz_xy_0_yy, g_z_0_x_0_xz_xy_0_yz, g_z_0_x_0_xz_xy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_xz_xy_0_xx[i] = -2.0 * g_x_xy_x_xx[i] * c_exps[i] + 4.0 * g_xzz_xy_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_xy_0_xy[i] = -2.0 * g_x_xy_x_xy[i] * c_exps[i] + 4.0 * g_xzz_xy_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_xy_0_xz[i] = -2.0 * g_x_xy_x_xz[i] * c_exps[i] + 4.0 * g_xzz_xy_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_xy_0_yy[i] = -2.0 * g_x_xy_x_yy[i] * c_exps[i] + 4.0 * g_xzz_xy_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_xy_0_yz[i] = -2.0 * g_x_xy_x_yz[i] * c_exps[i] + 4.0 * g_xzz_xy_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_xy_0_zz[i] = -2.0 * g_x_xy_x_zz[i] * c_exps[i] + 4.0 * g_xzz_xy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1380-1386)

    #pragma omp simd aligned(g_x_xz_x_xx, g_x_xz_x_xy, g_x_xz_x_xz, g_x_xz_x_yy, g_x_xz_x_yz, g_x_xz_x_zz, g_xzz_xz_x_xx, g_xzz_xz_x_xy, g_xzz_xz_x_xz, g_xzz_xz_x_yy, g_xzz_xz_x_yz, g_xzz_xz_x_zz, g_z_0_x_0_xz_xz_0_xx, g_z_0_x_0_xz_xz_0_xy, g_z_0_x_0_xz_xz_0_xz, g_z_0_x_0_xz_xz_0_yy, g_z_0_x_0_xz_xz_0_yz, g_z_0_x_0_xz_xz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_xz_xz_0_xx[i] = -2.0 * g_x_xz_x_xx[i] * c_exps[i] + 4.0 * g_xzz_xz_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_xz_0_xy[i] = -2.0 * g_x_xz_x_xy[i] * c_exps[i] + 4.0 * g_xzz_xz_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_xz_0_xz[i] = -2.0 * g_x_xz_x_xz[i] * c_exps[i] + 4.0 * g_xzz_xz_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_xz_0_yy[i] = -2.0 * g_x_xz_x_yy[i] * c_exps[i] + 4.0 * g_xzz_xz_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_xz_0_yz[i] = -2.0 * g_x_xz_x_yz[i] * c_exps[i] + 4.0 * g_xzz_xz_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_xz_0_zz[i] = -2.0 * g_x_xz_x_zz[i] * c_exps[i] + 4.0 * g_xzz_xz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1386-1392)

    #pragma omp simd aligned(g_x_yy_x_xx, g_x_yy_x_xy, g_x_yy_x_xz, g_x_yy_x_yy, g_x_yy_x_yz, g_x_yy_x_zz, g_xzz_yy_x_xx, g_xzz_yy_x_xy, g_xzz_yy_x_xz, g_xzz_yy_x_yy, g_xzz_yy_x_yz, g_xzz_yy_x_zz, g_z_0_x_0_xz_yy_0_xx, g_z_0_x_0_xz_yy_0_xy, g_z_0_x_0_xz_yy_0_xz, g_z_0_x_0_xz_yy_0_yy, g_z_0_x_0_xz_yy_0_yz, g_z_0_x_0_xz_yy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_xz_yy_0_xx[i] = -2.0 * g_x_yy_x_xx[i] * c_exps[i] + 4.0 * g_xzz_yy_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_yy_0_xy[i] = -2.0 * g_x_yy_x_xy[i] * c_exps[i] + 4.0 * g_xzz_yy_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_yy_0_xz[i] = -2.0 * g_x_yy_x_xz[i] * c_exps[i] + 4.0 * g_xzz_yy_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_yy_0_yy[i] = -2.0 * g_x_yy_x_yy[i] * c_exps[i] + 4.0 * g_xzz_yy_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_yy_0_yz[i] = -2.0 * g_x_yy_x_yz[i] * c_exps[i] + 4.0 * g_xzz_yy_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_yy_0_zz[i] = -2.0 * g_x_yy_x_zz[i] * c_exps[i] + 4.0 * g_xzz_yy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1392-1398)

    #pragma omp simd aligned(g_x_yz_x_xx, g_x_yz_x_xy, g_x_yz_x_xz, g_x_yz_x_yy, g_x_yz_x_yz, g_x_yz_x_zz, g_xzz_yz_x_xx, g_xzz_yz_x_xy, g_xzz_yz_x_xz, g_xzz_yz_x_yy, g_xzz_yz_x_yz, g_xzz_yz_x_zz, g_z_0_x_0_xz_yz_0_xx, g_z_0_x_0_xz_yz_0_xy, g_z_0_x_0_xz_yz_0_xz, g_z_0_x_0_xz_yz_0_yy, g_z_0_x_0_xz_yz_0_yz, g_z_0_x_0_xz_yz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_xz_yz_0_xx[i] = -2.0 * g_x_yz_x_xx[i] * c_exps[i] + 4.0 * g_xzz_yz_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_yz_0_xy[i] = -2.0 * g_x_yz_x_xy[i] * c_exps[i] + 4.0 * g_xzz_yz_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_yz_0_xz[i] = -2.0 * g_x_yz_x_xz[i] * c_exps[i] + 4.0 * g_xzz_yz_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_yz_0_yy[i] = -2.0 * g_x_yz_x_yy[i] * c_exps[i] + 4.0 * g_xzz_yz_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_yz_0_yz[i] = -2.0 * g_x_yz_x_yz[i] * c_exps[i] + 4.0 * g_xzz_yz_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_yz_0_zz[i] = -2.0 * g_x_yz_x_zz[i] * c_exps[i] + 4.0 * g_xzz_yz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1398-1404)

    #pragma omp simd aligned(g_x_zz_x_xx, g_x_zz_x_xy, g_x_zz_x_xz, g_x_zz_x_yy, g_x_zz_x_yz, g_x_zz_x_zz, g_xzz_zz_x_xx, g_xzz_zz_x_xy, g_xzz_zz_x_xz, g_xzz_zz_x_yy, g_xzz_zz_x_yz, g_xzz_zz_x_zz, g_z_0_x_0_xz_zz_0_xx, g_z_0_x_0_xz_zz_0_xy, g_z_0_x_0_xz_zz_0_xz, g_z_0_x_0_xz_zz_0_yy, g_z_0_x_0_xz_zz_0_yz, g_z_0_x_0_xz_zz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_xz_zz_0_xx[i] = -2.0 * g_x_zz_x_xx[i] * c_exps[i] + 4.0 * g_xzz_zz_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_zz_0_xy[i] = -2.0 * g_x_zz_x_xy[i] * c_exps[i] + 4.0 * g_xzz_zz_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_zz_0_xz[i] = -2.0 * g_x_zz_x_xz[i] * c_exps[i] + 4.0 * g_xzz_zz_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_zz_0_yy[i] = -2.0 * g_x_zz_x_yy[i] * c_exps[i] + 4.0 * g_xzz_zz_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_zz_0_yz[i] = -2.0 * g_x_zz_x_yz[i] * c_exps[i] + 4.0 * g_xzz_zz_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_zz_0_zz[i] = -2.0 * g_x_zz_x_zz[i] * c_exps[i] + 4.0 * g_xzz_zz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1404-1410)

    #pragma omp simd aligned(g_yyz_xx_x_xx, g_yyz_xx_x_xy, g_yyz_xx_x_xz, g_yyz_xx_x_yy, g_yyz_xx_x_yz, g_yyz_xx_x_zz, g_z_0_x_0_yy_xx_0_xx, g_z_0_x_0_yy_xx_0_xy, g_z_0_x_0_yy_xx_0_xz, g_z_0_x_0_yy_xx_0_yy, g_z_0_x_0_yy_xx_0_yz, g_z_0_x_0_yy_xx_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_yy_xx_0_xx[i] = 4.0 * g_yyz_xx_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_xx_0_xy[i] = 4.0 * g_yyz_xx_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_xx_0_xz[i] = 4.0 * g_yyz_xx_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_xx_0_yy[i] = 4.0 * g_yyz_xx_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_xx_0_yz[i] = 4.0 * g_yyz_xx_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_xx_0_zz[i] = 4.0 * g_yyz_xx_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1410-1416)

    #pragma omp simd aligned(g_yyz_xy_x_xx, g_yyz_xy_x_xy, g_yyz_xy_x_xz, g_yyz_xy_x_yy, g_yyz_xy_x_yz, g_yyz_xy_x_zz, g_z_0_x_0_yy_xy_0_xx, g_z_0_x_0_yy_xy_0_xy, g_z_0_x_0_yy_xy_0_xz, g_z_0_x_0_yy_xy_0_yy, g_z_0_x_0_yy_xy_0_yz, g_z_0_x_0_yy_xy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_yy_xy_0_xx[i] = 4.0 * g_yyz_xy_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_xy_0_xy[i] = 4.0 * g_yyz_xy_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_xy_0_xz[i] = 4.0 * g_yyz_xy_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_xy_0_yy[i] = 4.0 * g_yyz_xy_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_xy_0_yz[i] = 4.0 * g_yyz_xy_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_xy_0_zz[i] = 4.0 * g_yyz_xy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1416-1422)

    #pragma omp simd aligned(g_yyz_xz_x_xx, g_yyz_xz_x_xy, g_yyz_xz_x_xz, g_yyz_xz_x_yy, g_yyz_xz_x_yz, g_yyz_xz_x_zz, g_z_0_x_0_yy_xz_0_xx, g_z_0_x_0_yy_xz_0_xy, g_z_0_x_0_yy_xz_0_xz, g_z_0_x_0_yy_xz_0_yy, g_z_0_x_0_yy_xz_0_yz, g_z_0_x_0_yy_xz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_yy_xz_0_xx[i] = 4.0 * g_yyz_xz_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_xz_0_xy[i] = 4.0 * g_yyz_xz_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_xz_0_xz[i] = 4.0 * g_yyz_xz_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_xz_0_yy[i] = 4.0 * g_yyz_xz_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_xz_0_yz[i] = 4.0 * g_yyz_xz_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_xz_0_zz[i] = 4.0 * g_yyz_xz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1422-1428)

    #pragma omp simd aligned(g_yyz_yy_x_xx, g_yyz_yy_x_xy, g_yyz_yy_x_xz, g_yyz_yy_x_yy, g_yyz_yy_x_yz, g_yyz_yy_x_zz, g_z_0_x_0_yy_yy_0_xx, g_z_0_x_0_yy_yy_0_xy, g_z_0_x_0_yy_yy_0_xz, g_z_0_x_0_yy_yy_0_yy, g_z_0_x_0_yy_yy_0_yz, g_z_0_x_0_yy_yy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_yy_yy_0_xx[i] = 4.0 * g_yyz_yy_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_yy_0_xy[i] = 4.0 * g_yyz_yy_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_yy_0_xz[i] = 4.0 * g_yyz_yy_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_yy_0_yy[i] = 4.0 * g_yyz_yy_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_yy_0_yz[i] = 4.0 * g_yyz_yy_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_yy_0_zz[i] = 4.0 * g_yyz_yy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1428-1434)

    #pragma omp simd aligned(g_yyz_yz_x_xx, g_yyz_yz_x_xy, g_yyz_yz_x_xz, g_yyz_yz_x_yy, g_yyz_yz_x_yz, g_yyz_yz_x_zz, g_z_0_x_0_yy_yz_0_xx, g_z_0_x_0_yy_yz_0_xy, g_z_0_x_0_yy_yz_0_xz, g_z_0_x_0_yy_yz_0_yy, g_z_0_x_0_yy_yz_0_yz, g_z_0_x_0_yy_yz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_yy_yz_0_xx[i] = 4.0 * g_yyz_yz_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_yz_0_xy[i] = 4.0 * g_yyz_yz_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_yz_0_xz[i] = 4.0 * g_yyz_yz_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_yz_0_yy[i] = 4.0 * g_yyz_yz_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_yz_0_yz[i] = 4.0 * g_yyz_yz_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_yz_0_zz[i] = 4.0 * g_yyz_yz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1434-1440)

    #pragma omp simd aligned(g_yyz_zz_x_xx, g_yyz_zz_x_xy, g_yyz_zz_x_xz, g_yyz_zz_x_yy, g_yyz_zz_x_yz, g_yyz_zz_x_zz, g_z_0_x_0_yy_zz_0_xx, g_z_0_x_0_yy_zz_0_xy, g_z_0_x_0_yy_zz_0_xz, g_z_0_x_0_yy_zz_0_yy, g_z_0_x_0_yy_zz_0_yz, g_z_0_x_0_yy_zz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_yy_zz_0_xx[i] = 4.0 * g_yyz_zz_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_zz_0_xy[i] = 4.0 * g_yyz_zz_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_zz_0_xz[i] = 4.0 * g_yyz_zz_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_zz_0_yy[i] = 4.0 * g_yyz_zz_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_zz_0_yz[i] = 4.0 * g_yyz_zz_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_zz_0_zz[i] = 4.0 * g_yyz_zz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1440-1446)

    #pragma omp simd aligned(g_y_xx_x_xx, g_y_xx_x_xy, g_y_xx_x_xz, g_y_xx_x_yy, g_y_xx_x_yz, g_y_xx_x_zz, g_yzz_xx_x_xx, g_yzz_xx_x_xy, g_yzz_xx_x_xz, g_yzz_xx_x_yy, g_yzz_xx_x_yz, g_yzz_xx_x_zz, g_z_0_x_0_yz_xx_0_xx, g_z_0_x_0_yz_xx_0_xy, g_z_0_x_0_yz_xx_0_xz, g_z_0_x_0_yz_xx_0_yy, g_z_0_x_0_yz_xx_0_yz, g_z_0_x_0_yz_xx_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_yz_xx_0_xx[i] = -2.0 * g_y_xx_x_xx[i] * c_exps[i] + 4.0 * g_yzz_xx_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_xx_0_xy[i] = -2.0 * g_y_xx_x_xy[i] * c_exps[i] + 4.0 * g_yzz_xx_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_xx_0_xz[i] = -2.0 * g_y_xx_x_xz[i] * c_exps[i] + 4.0 * g_yzz_xx_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_xx_0_yy[i] = -2.0 * g_y_xx_x_yy[i] * c_exps[i] + 4.0 * g_yzz_xx_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_xx_0_yz[i] = -2.0 * g_y_xx_x_yz[i] * c_exps[i] + 4.0 * g_yzz_xx_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_xx_0_zz[i] = -2.0 * g_y_xx_x_zz[i] * c_exps[i] + 4.0 * g_yzz_xx_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1446-1452)

    #pragma omp simd aligned(g_y_xy_x_xx, g_y_xy_x_xy, g_y_xy_x_xz, g_y_xy_x_yy, g_y_xy_x_yz, g_y_xy_x_zz, g_yzz_xy_x_xx, g_yzz_xy_x_xy, g_yzz_xy_x_xz, g_yzz_xy_x_yy, g_yzz_xy_x_yz, g_yzz_xy_x_zz, g_z_0_x_0_yz_xy_0_xx, g_z_0_x_0_yz_xy_0_xy, g_z_0_x_0_yz_xy_0_xz, g_z_0_x_0_yz_xy_0_yy, g_z_0_x_0_yz_xy_0_yz, g_z_0_x_0_yz_xy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_yz_xy_0_xx[i] = -2.0 * g_y_xy_x_xx[i] * c_exps[i] + 4.0 * g_yzz_xy_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_xy_0_xy[i] = -2.0 * g_y_xy_x_xy[i] * c_exps[i] + 4.0 * g_yzz_xy_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_xy_0_xz[i] = -2.0 * g_y_xy_x_xz[i] * c_exps[i] + 4.0 * g_yzz_xy_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_xy_0_yy[i] = -2.0 * g_y_xy_x_yy[i] * c_exps[i] + 4.0 * g_yzz_xy_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_xy_0_yz[i] = -2.0 * g_y_xy_x_yz[i] * c_exps[i] + 4.0 * g_yzz_xy_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_xy_0_zz[i] = -2.0 * g_y_xy_x_zz[i] * c_exps[i] + 4.0 * g_yzz_xy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1452-1458)

    #pragma omp simd aligned(g_y_xz_x_xx, g_y_xz_x_xy, g_y_xz_x_xz, g_y_xz_x_yy, g_y_xz_x_yz, g_y_xz_x_zz, g_yzz_xz_x_xx, g_yzz_xz_x_xy, g_yzz_xz_x_xz, g_yzz_xz_x_yy, g_yzz_xz_x_yz, g_yzz_xz_x_zz, g_z_0_x_0_yz_xz_0_xx, g_z_0_x_0_yz_xz_0_xy, g_z_0_x_0_yz_xz_0_xz, g_z_0_x_0_yz_xz_0_yy, g_z_0_x_0_yz_xz_0_yz, g_z_0_x_0_yz_xz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_yz_xz_0_xx[i] = -2.0 * g_y_xz_x_xx[i] * c_exps[i] + 4.0 * g_yzz_xz_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_xz_0_xy[i] = -2.0 * g_y_xz_x_xy[i] * c_exps[i] + 4.0 * g_yzz_xz_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_xz_0_xz[i] = -2.0 * g_y_xz_x_xz[i] * c_exps[i] + 4.0 * g_yzz_xz_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_xz_0_yy[i] = -2.0 * g_y_xz_x_yy[i] * c_exps[i] + 4.0 * g_yzz_xz_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_xz_0_yz[i] = -2.0 * g_y_xz_x_yz[i] * c_exps[i] + 4.0 * g_yzz_xz_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_xz_0_zz[i] = -2.0 * g_y_xz_x_zz[i] * c_exps[i] + 4.0 * g_yzz_xz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1458-1464)

    #pragma omp simd aligned(g_y_yy_x_xx, g_y_yy_x_xy, g_y_yy_x_xz, g_y_yy_x_yy, g_y_yy_x_yz, g_y_yy_x_zz, g_yzz_yy_x_xx, g_yzz_yy_x_xy, g_yzz_yy_x_xz, g_yzz_yy_x_yy, g_yzz_yy_x_yz, g_yzz_yy_x_zz, g_z_0_x_0_yz_yy_0_xx, g_z_0_x_0_yz_yy_0_xy, g_z_0_x_0_yz_yy_0_xz, g_z_0_x_0_yz_yy_0_yy, g_z_0_x_0_yz_yy_0_yz, g_z_0_x_0_yz_yy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_yz_yy_0_xx[i] = -2.0 * g_y_yy_x_xx[i] * c_exps[i] + 4.0 * g_yzz_yy_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_yy_0_xy[i] = -2.0 * g_y_yy_x_xy[i] * c_exps[i] + 4.0 * g_yzz_yy_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_yy_0_xz[i] = -2.0 * g_y_yy_x_xz[i] * c_exps[i] + 4.0 * g_yzz_yy_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_yy_0_yy[i] = -2.0 * g_y_yy_x_yy[i] * c_exps[i] + 4.0 * g_yzz_yy_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_yy_0_yz[i] = -2.0 * g_y_yy_x_yz[i] * c_exps[i] + 4.0 * g_yzz_yy_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_yy_0_zz[i] = -2.0 * g_y_yy_x_zz[i] * c_exps[i] + 4.0 * g_yzz_yy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1464-1470)

    #pragma omp simd aligned(g_y_yz_x_xx, g_y_yz_x_xy, g_y_yz_x_xz, g_y_yz_x_yy, g_y_yz_x_yz, g_y_yz_x_zz, g_yzz_yz_x_xx, g_yzz_yz_x_xy, g_yzz_yz_x_xz, g_yzz_yz_x_yy, g_yzz_yz_x_yz, g_yzz_yz_x_zz, g_z_0_x_0_yz_yz_0_xx, g_z_0_x_0_yz_yz_0_xy, g_z_0_x_0_yz_yz_0_xz, g_z_0_x_0_yz_yz_0_yy, g_z_0_x_0_yz_yz_0_yz, g_z_0_x_0_yz_yz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_yz_yz_0_xx[i] = -2.0 * g_y_yz_x_xx[i] * c_exps[i] + 4.0 * g_yzz_yz_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_yz_0_xy[i] = -2.0 * g_y_yz_x_xy[i] * c_exps[i] + 4.0 * g_yzz_yz_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_yz_0_xz[i] = -2.0 * g_y_yz_x_xz[i] * c_exps[i] + 4.0 * g_yzz_yz_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_yz_0_yy[i] = -2.0 * g_y_yz_x_yy[i] * c_exps[i] + 4.0 * g_yzz_yz_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_yz_0_yz[i] = -2.0 * g_y_yz_x_yz[i] * c_exps[i] + 4.0 * g_yzz_yz_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_yz_0_zz[i] = -2.0 * g_y_yz_x_zz[i] * c_exps[i] + 4.0 * g_yzz_yz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1470-1476)

    #pragma omp simd aligned(g_y_zz_x_xx, g_y_zz_x_xy, g_y_zz_x_xz, g_y_zz_x_yy, g_y_zz_x_yz, g_y_zz_x_zz, g_yzz_zz_x_xx, g_yzz_zz_x_xy, g_yzz_zz_x_xz, g_yzz_zz_x_yy, g_yzz_zz_x_yz, g_yzz_zz_x_zz, g_z_0_x_0_yz_zz_0_xx, g_z_0_x_0_yz_zz_0_xy, g_z_0_x_0_yz_zz_0_xz, g_z_0_x_0_yz_zz_0_yy, g_z_0_x_0_yz_zz_0_yz, g_z_0_x_0_yz_zz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_yz_zz_0_xx[i] = -2.0 * g_y_zz_x_xx[i] * c_exps[i] + 4.0 * g_yzz_zz_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_zz_0_xy[i] = -2.0 * g_y_zz_x_xy[i] * c_exps[i] + 4.0 * g_yzz_zz_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_zz_0_xz[i] = -2.0 * g_y_zz_x_xz[i] * c_exps[i] + 4.0 * g_yzz_zz_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_zz_0_yy[i] = -2.0 * g_y_zz_x_yy[i] * c_exps[i] + 4.0 * g_yzz_zz_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_zz_0_yz[i] = -2.0 * g_y_zz_x_yz[i] * c_exps[i] + 4.0 * g_yzz_zz_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_zz_0_zz[i] = -2.0 * g_y_zz_x_zz[i] * c_exps[i] + 4.0 * g_yzz_zz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1476-1482)

    #pragma omp simd aligned(g_z_0_x_0_zz_xx_0_xx, g_z_0_x_0_zz_xx_0_xy, g_z_0_x_0_zz_xx_0_xz, g_z_0_x_0_zz_xx_0_yy, g_z_0_x_0_zz_xx_0_yz, g_z_0_x_0_zz_xx_0_zz, g_z_xx_x_xx, g_z_xx_x_xy, g_z_xx_x_xz, g_z_xx_x_yy, g_z_xx_x_yz, g_z_xx_x_zz, g_zzz_xx_x_xx, g_zzz_xx_x_xy, g_zzz_xx_x_xz, g_zzz_xx_x_yy, g_zzz_xx_x_yz, g_zzz_xx_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_zz_xx_0_xx[i] = -4.0 * g_z_xx_x_xx[i] * c_exps[i] + 4.0 * g_zzz_xx_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_xx_0_xy[i] = -4.0 * g_z_xx_x_xy[i] * c_exps[i] + 4.0 * g_zzz_xx_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_xx_0_xz[i] = -4.0 * g_z_xx_x_xz[i] * c_exps[i] + 4.0 * g_zzz_xx_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_xx_0_yy[i] = -4.0 * g_z_xx_x_yy[i] * c_exps[i] + 4.0 * g_zzz_xx_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_xx_0_yz[i] = -4.0 * g_z_xx_x_yz[i] * c_exps[i] + 4.0 * g_zzz_xx_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_xx_0_zz[i] = -4.0 * g_z_xx_x_zz[i] * c_exps[i] + 4.0 * g_zzz_xx_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1482-1488)

    #pragma omp simd aligned(g_z_0_x_0_zz_xy_0_xx, g_z_0_x_0_zz_xy_0_xy, g_z_0_x_0_zz_xy_0_xz, g_z_0_x_0_zz_xy_0_yy, g_z_0_x_0_zz_xy_0_yz, g_z_0_x_0_zz_xy_0_zz, g_z_xy_x_xx, g_z_xy_x_xy, g_z_xy_x_xz, g_z_xy_x_yy, g_z_xy_x_yz, g_z_xy_x_zz, g_zzz_xy_x_xx, g_zzz_xy_x_xy, g_zzz_xy_x_xz, g_zzz_xy_x_yy, g_zzz_xy_x_yz, g_zzz_xy_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_zz_xy_0_xx[i] = -4.0 * g_z_xy_x_xx[i] * c_exps[i] + 4.0 * g_zzz_xy_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_xy_0_xy[i] = -4.0 * g_z_xy_x_xy[i] * c_exps[i] + 4.0 * g_zzz_xy_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_xy_0_xz[i] = -4.0 * g_z_xy_x_xz[i] * c_exps[i] + 4.0 * g_zzz_xy_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_xy_0_yy[i] = -4.0 * g_z_xy_x_yy[i] * c_exps[i] + 4.0 * g_zzz_xy_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_xy_0_yz[i] = -4.0 * g_z_xy_x_yz[i] * c_exps[i] + 4.0 * g_zzz_xy_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_xy_0_zz[i] = -4.0 * g_z_xy_x_zz[i] * c_exps[i] + 4.0 * g_zzz_xy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1488-1494)

    #pragma omp simd aligned(g_z_0_x_0_zz_xz_0_xx, g_z_0_x_0_zz_xz_0_xy, g_z_0_x_0_zz_xz_0_xz, g_z_0_x_0_zz_xz_0_yy, g_z_0_x_0_zz_xz_0_yz, g_z_0_x_0_zz_xz_0_zz, g_z_xz_x_xx, g_z_xz_x_xy, g_z_xz_x_xz, g_z_xz_x_yy, g_z_xz_x_yz, g_z_xz_x_zz, g_zzz_xz_x_xx, g_zzz_xz_x_xy, g_zzz_xz_x_xz, g_zzz_xz_x_yy, g_zzz_xz_x_yz, g_zzz_xz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_zz_xz_0_xx[i] = -4.0 * g_z_xz_x_xx[i] * c_exps[i] + 4.0 * g_zzz_xz_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_xz_0_xy[i] = -4.0 * g_z_xz_x_xy[i] * c_exps[i] + 4.0 * g_zzz_xz_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_xz_0_xz[i] = -4.0 * g_z_xz_x_xz[i] * c_exps[i] + 4.0 * g_zzz_xz_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_xz_0_yy[i] = -4.0 * g_z_xz_x_yy[i] * c_exps[i] + 4.0 * g_zzz_xz_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_xz_0_yz[i] = -4.0 * g_z_xz_x_yz[i] * c_exps[i] + 4.0 * g_zzz_xz_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_xz_0_zz[i] = -4.0 * g_z_xz_x_zz[i] * c_exps[i] + 4.0 * g_zzz_xz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1494-1500)

    #pragma omp simd aligned(g_z_0_x_0_zz_yy_0_xx, g_z_0_x_0_zz_yy_0_xy, g_z_0_x_0_zz_yy_0_xz, g_z_0_x_0_zz_yy_0_yy, g_z_0_x_0_zz_yy_0_yz, g_z_0_x_0_zz_yy_0_zz, g_z_yy_x_xx, g_z_yy_x_xy, g_z_yy_x_xz, g_z_yy_x_yy, g_z_yy_x_yz, g_z_yy_x_zz, g_zzz_yy_x_xx, g_zzz_yy_x_xy, g_zzz_yy_x_xz, g_zzz_yy_x_yy, g_zzz_yy_x_yz, g_zzz_yy_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_zz_yy_0_xx[i] = -4.0 * g_z_yy_x_xx[i] * c_exps[i] + 4.0 * g_zzz_yy_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_yy_0_xy[i] = -4.0 * g_z_yy_x_xy[i] * c_exps[i] + 4.0 * g_zzz_yy_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_yy_0_xz[i] = -4.0 * g_z_yy_x_xz[i] * c_exps[i] + 4.0 * g_zzz_yy_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_yy_0_yy[i] = -4.0 * g_z_yy_x_yy[i] * c_exps[i] + 4.0 * g_zzz_yy_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_yy_0_yz[i] = -4.0 * g_z_yy_x_yz[i] * c_exps[i] + 4.0 * g_zzz_yy_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_yy_0_zz[i] = -4.0 * g_z_yy_x_zz[i] * c_exps[i] + 4.0 * g_zzz_yy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1500-1506)

    #pragma omp simd aligned(g_z_0_x_0_zz_yz_0_xx, g_z_0_x_0_zz_yz_0_xy, g_z_0_x_0_zz_yz_0_xz, g_z_0_x_0_zz_yz_0_yy, g_z_0_x_0_zz_yz_0_yz, g_z_0_x_0_zz_yz_0_zz, g_z_yz_x_xx, g_z_yz_x_xy, g_z_yz_x_xz, g_z_yz_x_yy, g_z_yz_x_yz, g_z_yz_x_zz, g_zzz_yz_x_xx, g_zzz_yz_x_xy, g_zzz_yz_x_xz, g_zzz_yz_x_yy, g_zzz_yz_x_yz, g_zzz_yz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_zz_yz_0_xx[i] = -4.0 * g_z_yz_x_xx[i] * c_exps[i] + 4.0 * g_zzz_yz_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_yz_0_xy[i] = -4.0 * g_z_yz_x_xy[i] * c_exps[i] + 4.0 * g_zzz_yz_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_yz_0_xz[i] = -4.0 * g_z_yz_x_xz[i] * c_exps[i] + 4.0 * g_zzz_yz_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_yz_0_yy[i] = -4.0 * g_z_yz_x_yy[i] * c_exps[i] + 4.0 * g_zzz_yz_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_yz_0_yz[i] = -4.0 * g_z_yz_x_yz[i] * c_exps[i] + 4.0 * g_zzz_yz_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_yz_0_zz[i] = -4.0 * g_z_yz_x_zz[i] * c_exps[i] + 4.0 * g_zzz_yz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1506-1512)

    #pragma omp simd aligned(g_z_0_x_0_zz_zz_0_xx, g_z_0_x_0_zz_zz_0_xy, g_z_0_x_0_zz_zz_0_xz, g_z_0_x_0_zz_zz_0_yy, g_z_0_x_0_zz_zz_0_yz, g_z_0_x_0_zz_zz_0_zz, g_z_zz_x_xx, g_z_zz_x_xy, g_z_zz_x_xz, g_z_zz_x_yy, g_z_zz_x_yz, g_z_zz_x_zz, g_zzz_zz_x_xx, g_zzz_zz_x_xy, g_zzz_zz_x_xz, g_zzz_zz_x_yy, g_zzz_zz_x_yz, g_zzz_zz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_zz_zz_0_xx[i] = -4.0 * g_z_zz_x_xx[i] * c_exps[i] + 4.0 * g_zzz_zz_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_zz_0_xy[i] = -4.0 * g_z_zz_x_xy[i] * c_exps[i] + 4.0 * g_zzz_zz_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_zz_0_xz[i] = -4.0 * g_z_zz_x_xz[i] * c_exps[i] + 4.0 * g_zzz_zz_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_zz_0_yy[i] = -4.0 * g_z_zz_x_yy[i] * c_exps[i] + 4.0 * g_zzz_zz_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_zz_0_yz[i] = -4.0 * g_z_zz_x_yz[i] * c_exps[i] + 4.0 * g_zzz_zz_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_zz_0_zz[i] = -4.0 * g_z_zz_x_zz[i] * c_exps[i] + 4.0 * g_zzz_zz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1512-1518)

    #pragma omp simd aligned(g_xxz_xx_y_xx, g_xxz_xx_y_xy, g_xxz_xx_y_xz, g_xxz_xx_y_yy, g_xxz_xx_y_yz, g_xxz_xx_y_zz, g_z_0_y_0_xx_xx_0_xx, g_z_0_y_0_xx_xx_0_xy, g_z_0_y_0_xx_xx_0_xz, g_z_0_y_0_xx_xx_0_yy, g_z_0_y_0_xx_xx_0_yz, g_z_0_y_0_xx_xx_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_xx_xx_0_xx[i] = 4.0 * g_xxz_xx_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_xx_0_xy[i] = 4.0 * g_xxz_xx_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_xx_0_xz[i] = 4.0 * g_xxz_xx_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_xx_0_yy[i] = 4.0 * g_xxz_xx_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_xx_0_yz[i] = 4.0 * g_xxz_xx_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_xx_0_zz[i] = 4.0 * g_xxz_xx_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1518-1524)

    #pragma omp simd aligned(g_xxz_xy_y_xx, g_xxz_xy_y_xy, g_xxz_xy_y_xz, g_xxz_xy_y_yy, g_xxz_xy_y_yz, g_xxz_xy_y_zz, g_z_0_y_0_xx_xy_0_xx, g_z_0_y_0_xx_xy_0_xy, g_z_0_y_0_xx_xy_0_xz, g_z_0_y_0_xx_xy_0_yy, g_z_0_y_0_xx_xy_0_yz, g_z_0_y_0_xx_xy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_xx_xy_0_xx[i] = 4.0 * g_xxz_xy_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_xy_0_xy[i] = 4.0 * g_xxz_xy_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_xy_0_xz[i] = 4.0 * g_xxz_xy_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_xy_0_yy[i] = 4.0 * g_xxz_xy_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_xy_0_yz[i] = 4.0 * g_xxz_xy_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_xy_0_zz[i] = 4.0 * g_xxz_xy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1524-1530)

    #pragma omp simd aligned(g_xxz_xz_y_xx, g_xxz_xz_y_xy, g_xxz_xz_y_xz, g_xxz_xz_y_yy, g_xxz_xz_y_yz, g_xxz_xz_y_zz, g_z_0_y_0_xx_xz_0_xx, g_z_0_y_0_xx_xz_0_xy, g_z_0_y_0_xx_xz_0_xz, g_z_0_y_0_xx_xz_0_yy, g_z_0_y_0_xx_xz_0_yz, g_z_0_y_0_xx_xz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_xx_xz_0_xx[i] = 4.0 * g_xxz_xz_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_xz_0_xy[i] = 4.0 * g_xxz_xz_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_xz_0_xz[i] = 4.0 * g_xxz_xz_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_xz_0_yy[i] = 4.0 * g_xxz_xz_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_xz_0_yz[i] = 4.0 * g_xxz_xz_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_xz_0_zz[i] = 4.0 * g_xxz_xz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1530-1536)

    #pragma omp simd aligned(g_xxz_yy_y_xx, g_xxz_yy_y_xy, g_xxz_yy_y_xz, g_xxz_yy_y_yy, g_xxz_yy_y_yz, g_xxz_yy_y_zz, g_z_0_y_0_xx_yy_0_xx, g_z_0_y_0_xx_yy_0_xy, g_z_0_y_0_xx_yy_0_xz, g_z_0_y_0_xx_yy_0_yy, g_z_0_y_0_xx_yy_0_yz, g_z_0_y_0_xx_yy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_xx_yy_0_xx[i] = 4.0 * g_xxz_yy_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_yy_0_xy[i] = 4.0 * g_xxz_yy_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_yy_0_xz[i] = 4.0 * g_xxz_yy_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_yy_0_yy[i] = 4.0 * g_xxz_yy_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_yy_0_yz[i] = 4.0 * g_xxz_yy_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_yy_0_zz[i] = 4.0 * g_xxz_yy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1536-1542)

    #pragma omp simd aligned(g_xxz_yz_y_xx, g_xxz_yz_y_xy, g_xxz_yz_y_xz, g_xxz_yz_y_yy, g_xxz_yz_y_yz, g_xxz_yz_y_zz, g_z_0_y_0_xx_yz_0_xx, g_z_0_y_0_xx_yz_0_xy, g_z_0_y_0_xx_yz_0_xz, g_z_0_y_0_xx_yz_0_yy, g_z_0_y_0_xx_yz_0_yz, g_z_0_y_0_xx_yz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_xx_yz_0_xx[i] = 4.0 * g_xxz_yz_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_yz_0_xy[i] = 4.0 * g_xxz_yz_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_yz_0_xz[i] = 4.0 * g_xxz_yz_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_yz_0_yy[i] = 4.0 * g_xxz_yz_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_yz_0_yz[i] = 4.0 * g_xxz_yz_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_yz_0_zz[i] = 4.0 * g_xxz_yz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1542-1548)

    #pragma omp simd aligned(g_xxz_zz_y_xx, g_xxz_zz_y_xy, g_xxz_zz_y_xz, g_xxz_zz_y_yy, g_xxz_zz_y_yz, g_xxz_zz_y_zz, g_z_0_y_0_xx_zz_0_xx, g_z_0_y_0_xx_zz_0_xy, g_z_0_y_0_xx_zz_0_xz, g_z_0_y_0_xx_zz_0_yy, g_z_0_y_0_xx_zz_0_yz, g_z_0_y_0_xx_zz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_xx_zz_0_xx[i] = 4.0 * g_xxz_zz_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_zz_0_xy[i] = 4.0 * g_xxz_zz_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_zz_0_xz[i] = 4.0 * g_xxz_zz_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_zz_0_yy[i] = 4.0 * g_xxz_zz_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_zz_0_yz[i] = 4.0 * g_xxz_zz_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_zz_0_zz[i] = 4.0 * g_xxz_zz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1548-1554)

    #pragma omp simd aligned(g_xyz_xx_y_xx, g_xyz_xx_y_xy, g_xyz_xx_y_xz, g_xyz_xx_y_yy, g_xyz_xx_y_yz, g_xyz_xx_y_zz, g_z_0_y_0_xy_xx_0_xx, g_z_0_y_0_xy_xx_0_xy, g_z_0_y_0_xy_xx_0_xz, g_z_0_y_0_xy_xx_0_yy, g_z_0_y_0_xy_xx_0_yz, g_z_0_y_0_xy_xx_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_xy_xx_0_xx[i] = 4.0 * g_xyz_xx_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_xx_0_xy[i] = 4.0 * g_xyz_xx_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_xx_0_xz[i] = 4.0 * g_xyz_xx_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_xx_0_yy[i] = 4.0 * g_xyz_xx_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_xx_0_yz[i] = 4.0 * g_xyz_xx_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_xx_0_zz[i] = 4.0 * g_xyz_xx_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1554-1560)

    #pragma omp simd aligned(g_xyz_xy_y_xx, g_xyz_xy_y_xy, g_xyz_xy_y_xz, g_xyz_xy_y_yy, g_xyz_xy_y_yz, g_xyz_xy_y_zz, g_z_0_y_0_xy_xy_0_xx, g_z_0_y_0_xy_xy_0_xy, g_z_0_y_0_xy_xy_0_xz, g_z_0_y_0_xy_xy_0_yy, g_z_0_y_0_xy_xy_0_yz, g_z_0_y_0_xy_xy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_xy_xy_0_xx[i] = 4.0 * g_xyz_xy_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_xy_0_xy[i] = 4.0 * g_xyz_xy_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_xy_0_xz[i] = 4.0 * g_xyz_xy_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_xy_0_yy[i] = 4.0 * g_xyz_xy_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_xy_0_yz[i] = 4.0 * g_xyz_xy_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_xy_0_zz[i] = 4.0 * g_xyz_xy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1560-1566)

    #pragma omp simd aligned(g_xyz_xz_y_xx, g_xyz_xz_y_xy, g_xyz_xz_y_xz, g_xyz_xz_y_yy, g_xyz_xz_y_yz, g_xyz_xz_y_zz, g_z_0_y_0_xy_xz_0_xx, g_z_0_y_0_xy_xz_0_xy, g_z_0_y_0_xy_xz_0_xz, g_z_0_y_0_xy_xz_0_yy, g_z_0_y_0_xy_xz_0_yz, g_z_0_y_0_xy_xz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_xy_xz_0_xx[i] = 4.0 * g_xyz_xz_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_xz_0_xy[i] = 4.0 * g_xyz_xz_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_xz_0_xz[i] = 4.0 * g_xyz_xz_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_xz_0_yy[i] = 4.0 * g_xyz_xz_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_xz_0_yz[i] = 4.0 * g_xyz_xz_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_xz_0_zz[i] = 4.0 * g_xyz_xz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1566-1572)

    #pragma omp simd aligned(g_xyz_yy_y_xx, g_xyz_yy_y_xy, g_xyz_yy_y_xz, g_xyz_yy_y_yy, g_xyz_yy_y_yz, g_xyz_yy_y_zz, g_z_0_y_0_xy_yy_0_xx, g_z_0_y_0_xy_yy_0_xy, g_z_0_y_0_xy_yy_0_xz, g_z_0_y_0_xy_yy_0_yy, g_z_0_y_0_xy_yy_0_yz, g_z_0_y_0_xy_yy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_xy_yy_0_xx[i] = 4.0 * g_xyz_yy_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_yy_0_xy[i] = 4.0 * g_xyz_yy_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_yy_0_xz[i] = 4.0 * g_xyz_yy_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_yy_0_yy[i] = 4.0 * g_xyz_yy_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_yy_0_yz[i] = 4.0 * g_xyz_yy_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_yy_0_zz[i] = 4.0 * g_xyz_yy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1572-1578)

    #pragma omp simd aligned(g_xyz_yz_y_xx, g_xyz_yz_y_xy, g_xyz_yz_y_xz, g_xyz_yz_y_yy, g_xyz_yz_y_yz, g_xyz_yz_y_zz, g_z_0_y_0_xy_yz_0_xx, g_z_0_y_0_xy_yz_0_xy, g_z_0_y_0_xy_yz_0_xz, g_z_0_y_0_xy_yz_0_yy, g_z_0_y_0_xy_yz_0_yz, g_z_0_y_0_xy_yz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_xy_yz_0_xx[i] = 4.0 * g_xyz_yz_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_yz_0_xy[i] = 4.0 * g_xyz_yz_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_yz_0_xz[i] = 4.0 * g_xyz_yz_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_yz_0_yy[i] = 4.0 * g_xyz_yz_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_yz_0_yz[i] = 4.0 * g_xyz_yz_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_yz_0_zz[i] = 4.0 * g_xyz_yz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1578-1584)

    #pragma omp simd aligned(g_xyz_zz_y_xx, g_xyz_zz_y_xy, g_xyz_zz_y_xz, g_xyz_zz_y_yy, g_xyz_zz_y_yz, g_xyz_zz_y_zz, g_z_0_y_0_xy_zz_0_xx, g_z_0_y_0_xy_zz_0_xy, g_z_0_y_0_xy_zz_0_xz, g_z_0_y_0_xy_zz_0_yy, g_z_0_y_0_xy_zz_0_yz, g_z_0_y_0_xy_zz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_xy_zz_0_xx[i] = 4.0 * g_xyz_zz_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_zz_0_xy[i] = 4.0 * g_xyz_zz_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_zz_0_xz[i] = 4.0 * g_xyz_zz_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_zz_0_yy[i] = 4.0 * g_xyz_zz_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_zz_0_yz[i] = 4.0 * g_xyz_zz_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_zz_0_zz[i] = 4.0 * g_xyz_zz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1584-1590)

    #pragma omp simd aligned(g_x_xx_y_xx, g_x_xx_y_xy, g_x_xx_y_xz, g_x_xx_y_yy, g_x_xx_y_yz, g_x_xx_y_zz, g_xzz_xx_y_xx, g_xzz_xx_y_xy, g_xzz_xx_y_xz, g_xzz_xx_y_yy, g_xzz_xx_y_yz, g_xzz_xx_y_zz, g_z_0_y_0_xz_xx_0_xx, g_z_0_y_0_xz_xx_0_xy, g_z_0_y_0_xz_xx_0_xz, g_z_0_y_0_xz_xx_0_yy, g_z_0_y_0_xz_xx_0_yz, g_z_0_y_0_xz_xx_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_xz_xx_0_xx[i] = -2.0 * g_x_xx_y_xx[i] * c_exps[i] + 4.0 * g_xzz_xx_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_xx_0_xy[i] = -2.0 * g_x_xx_y_xy[i] * c_exps[i] + 4.0 * g_xzz_xx_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_xx_0_xz[i] = -2.0 * g_x_xx_y_xz[i] * c_exps[i] + 4.0 * g_xzz_xx_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_xx_0_yy[i] = -2.0 * g_x_xx_y_yy[i] * c_exps[i] + 4.0 * g_xzz_xx_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_xx_0_yz[i] = -2.0 * g_x_xx_y_yz[i] * c_exps[i] + 4.0 * g_xzz_xx_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_xx_0_zz[i] = -2.0 * g_x_xx_y_zz[i] * c_exps[i] + 4.0 * g_xzz_xx_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1590-1596)

    #pragma omp simd aligned(g_x_xy_y_xx, g_x_xy_y_xy, g_x_xy_y_xz, g_x_xy_y_yy, g_x_xy_y_yz, g_x_xy_y_zz, g_xzz_xy_y_xx, g_xzz_xy_y_xy, g_xzz_xy_y_xz, g_xzz_xy_y_yy, g_xzz_xy_y_yz, g_xzz_xy_y_zz, g_z_0_y_0_xz_xy_0_xx, g_z_0_y_0_xz_xy_0_xy, g_z_0_y_0_xz_xy_0_xz, g_z_0_y_0_xz_xy_0_yy, g_z_0_y_0_xz_xy_0_yz, g_z_0_y_0_xz_xy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_xz_xy_0_xx[i] = -2.0 * g_x_xy_y_xx[i] * c_exps[i] + 4.0 * g_xzz_xy_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_xy_0_xy[i] = -2.0 * g_x_xy_y_xy[i] * c_exps[i] + 4.0 * g_xzz_xy_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_xy_0_xz[i] = -2.0 * g_x_xy_y_xz[i] * c_exps[i] + 4.0 * g_xzz_xy_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_xy_0_yy[i] = -2.0 * g_x_xy_y_yy[i] * c_exps[i] + 4.0 * g_xzz_xy_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_xy_0_yz[i] = -2.0 * g_x_xy_y_yz[i] * c_exps[i] + 4.0 * g_xzz_xy_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_xy_0_zz[i] = -2.0 * g_x_xy_y_zz[i] * c_exps[i] + 4.0 * g_xzz_xy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1596-1602)

    #pragma omp simd aligned(g_x_xz_y_xx, g_x_xz_y_xy, g_x_xz_y_xz, g_x_xz_y_yy, g_x_xz_y_yz, g_x_xz_y_zz, g_xzz_xz_y_xx, g_xzz_xz_y_xy, g_xzz_xz_y_xz, g_xzz_xz_y_yy, g_xzz_xz_y_yz, g_xzz_xz_y_zz, g_z_0_y_0_xz_xz_0_xx, g_z_0_y_0_xz_xz_0_xy, g_z_0_y_0_xz_xz_0_xz, g_z_0_y_0_xz_xz_0_yy, g_z_0_y_0_xz_xz_0_yz, g_z_0_y_0_xz_xz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_xz_xz_0_xx[i] = -2.0 * g_x_xz_y_xx[i] * c_exps[i] + 4.0 * g_xzz_xz_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_xz_0_xy[i] = -2.0 * g_x_xz_y_xy[i] * c_exps[i] + 4.0 * g_xzz_xz_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_xz_0_xz[i] = -2.0 * g_x_xz_y_xz[i] * c_exps[i] + 4.0 * g_xzz_xz_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_xz_0_yy[i] = -2.0 * g_x_xz_y_yy[i] * c_exps[i] + 4.0 * g_xzz_xz_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_xz_0_yz[i] = -2.0 * g_x_xz_y_yz[i] * c_exps[i] + 4.0 * g_xzz_xz_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_xz_0_zz[i] = -2.0 * g_x_xz_y_zz[i] * c_exps[i] + 4.0 * g_xzz_xz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1602-1608)

    #pragma omp simd aligned(g_x_yy_y_xx, g_x_yy_y_xy, g_x_yy_y_xz, g_x_yy_y_yy, g_x_yy_y_yz, g_x_yy_y_zz, g_xzz_yy_y_xx, g_xzz_yy_y_xy, g_xzz_yy_y_xz, g_xzz_yy_y_yy, g_xzz_yy_y_yz, g_xzz_yy_y_zz, g_z_0_y_0_xz_yy_0_xx, g_z_0_y_0_xz_yy_0_xy, g_z_0_y_0_xz_yy_0_xz, g_z_0_y_0_xz_yy_0_yy, g_z_0_y_0_xz_yy_0_yz, g_z_0_y_0_xz_yy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_xz_yy_0_xx[i] = -2.0 * g_x_yy_y_xx[i] * c_exps[i] + 4.0 * g_xzz_yy_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_yy_0_xy[i] = -2.0 * g_x_yy_y_xy[i] * c_exps[i] + 4.0 * g_xzz_yy_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_yy_0_xz[i] = -2.0 * g_x_yy_y_xz[i] * c_exps[i] + 4.0 * g_xzz_yy_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_yy_0_yy[i] = -2.0 * g_x_yy_y_yy[i] * c_exps[i] + 4.0 * g_xzz_yy_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_yy_0_yz[i] = -2.0 * g_x_yy_y_yz[i] * c_exps[i] + 4.0 * g_xzz_yy_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_yy_0_zz[i] = -2.0 * g_x_yy_y_zz[i] * c_exps[i] + 4.0 * g_xzz_yy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1608-1614)

    #pragma omp simd aligned(g_x_yz_y_xx, g_x_yz_y_xy, g_x_yz_y_xz, g_x_yz_y_yy, g_x_yz_y_yz, g_x_yz_y_zz, g_xzz_yz_y_xx, g_xzz_yz_y_xy, g_xzz_yz_y_xz, g_xzz_yz_y_yy, g_xzz_yz_y_yz, g_xzz_yz_y_zz, g_z_0_y_0_xz_yz_0_xx, g_z_0_y_0_xz_yz_0_xy, g_z_0_y_0_xz_yz_0_xz, g_z_0_y_0_xz_yz_0_yy, g_z_0_y_0_xz_yz_0_yz, g_z_0_y_0_xz_yz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_xz_yz_0_xx[i] = -2.0 * g_x_yz_y_xx[i] * c_exps[i] + 4.0 * g_xzz_yz_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_yz_0_xy[i] = -2.0 * g_x_yz_y_xy[i] * c_exps[i] + 4.0 * g_xzz_yz_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_yz_0_xz[i] = -2.0 * g_x_yz_y_xz[i] * c_exps[i] + 4.0 * g_xzz_yz_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_yz_0_yy[i] = -2.0 * g_x_yz_y_yy[i] * c_exps[i] + 4.0 * g_xzz_yz_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_yz_0_yz[i] = -2.0 * g_x_yz_y_yz[i] * c_exps[i] + 4.0 * g_xzz_yz_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_yz_0_zz[i] = -2.0 * g_x_yz_y_zz[i] * c_exps[i] + 4.0 * g_xzz_yz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1614-1620)

    #pragma omp simd aligned(g_x_zz_y_xx, g_x_zz_y_xy, g_x_zz_y_xz, g_x_zz_y_yy, g_x_zz_y_yz, g_x_zz_y_zz, g_xzz_zz_y_xx, g_xzz_zz_y_xy, g_xzz_zz_y_xz, g_xzz_zz_y_yy, g_xzz_zz_y_yz, g_xzz_zz_y_zz, g_z_0_y_0_xz_zz_0_xx, g_z_0_y_0_xz_zz_0_xy, g_z_0_y_0_xz_zz_0_xz, g_z_0_y_0_xz_zz_0_yy, g_z_0_y_0_xz_zz_0_yz, g_z_0_y_0_xz_zz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_xz_zz_0_xx[i] = -2.0 * g_x_zz_y_xx[i] * c_exps[i] + 4.0 * g_xzz_zz_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_zz_0_xy[i] = -2.0 * g_x_zz_y_xy[i] * c_exps[i] + 4.0 * g_xzz_zz_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_zz_0_xz[i] = -2.0 * g_x_zz_y_xz[i] * c_exps[i] + 4.0 * g_xzz_zz_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_zz_0_yy[i] = -2.0 * g_x_zz_y_yy[i] * c_exps[i] + 4.0 * g_xzz_zz_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_zz_0_yz[i] = -2.0 * g_x_zz_y_yz[i] * c_exps[i] + 4.0 * g_xzz_zz_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_zz_0_zz[i] = -2.0 * g_x_zz_y_zz[i] * c_exps[i] + 4.0 * g_xzz_zz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1620-1626)

    #pragma omp simd aligned(g_yyz_xx_y_xx, g_yyz_xx_y_xy, g_yyz_xx_y_xz, g_yyz_xx_y_yy, g_yyz_xx_y_yz, g_yyz_xx_y_zz, g_z_0_y_0_yy_xx_0_xx, g_z_0_y_0_yy_xx_0_xy, g_z_0_y_0_yy_xx_0_xz, g_z_0_y_0_yy_xx_0_yy, g_z_0_y_0_yy_xx_0_yz, g_z_0_y_0_yy_xx_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_yy_xx_0_xx[i] = 4.0 * g_yyz_xx_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_xx_0_xy[i] = 4.0 * g_yyz_xx_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_xx_0_xz[i] = 4.0 * g_yyz_xx_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_xx_0_yy[i] = 4.0 * g_yyz_xx_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_xx_0_yz[i] = 4.0 * g_yyz_xx_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_xx_0_zz[i] = 4.0 * g_yyz_xx_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1626-1632)

    #pragma omp simd aligned(g_yyz_xy_y_xx, g_yyz_xy_y_xy, g_yyz_xy_y_xz, g_yyz_xy_y_yy, g_yyz_xy_y_yz, g_yyz_xy_y_zz, g_z_0_y_0_yy_xy_0_xx, g_z_0_y_0_yy_xy_0_xy, g_z_0_y_0_yy_xy_0_xz, g_z_0_y_0_yy_xy_0_yy, g_z_0_y_0_yy_xy_0_yz, g_z_0_y_0_yy_xy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_yy_xy_0_xx[i] = 4.0 * g_yyz_xy_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_xy_0_xy[i] = 4.0 * g_yyz_xy_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_xy_0_xz[i] = 4.0 * g_yyz_xy_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_xy_0_yy[i] = 4.0 * g_yyz_xy_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_xy_0_yz[i] = 4.0 * g_yyz_xy_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_xy_0_zz[i] = 4.0 * g_yyz_xy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1632-1638)

    #pragma omp simd aligned(g_yyz_xz_y_xx, g_yyz_xz_y_xy, g_yyz_xz_y_xz, g_yyz_xz_y_yy, g_yyz_xz_y_yz, g_yyz_xz_y_zz, g_z_0_y_0_yy_xz_0_xx, g_z_0_y_0_yy_xz_0_xy, g_z_0_y_0_yy_xz_0_xz, g_z_0_y_0_yy_xz_0_yy, g_z_0_y_0_yy_xz_0_yz, g_z_0_y_0_yy_xz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_yy_xz_0_xx[i] = 4.0 * g_yyz_xz_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_xz_0_xy[i] = 4.0 * g_yyz_xz_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_xz_0_xz[i] = 4.0 * g_yyz_xz_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_xz_0_yy[i] = 4.0 * g_yyz_xz_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_xz_0_yz[i] = 4.0 * g_yyz_xz_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_xz_0_zz[i] = 4.0 * g_yyz_xz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1638-1644)

    #pragma omp simd aligned(g_yyz_yy_y_xx, g_yyz_yy_y_xy, g_yyz_yy_y_xz, g_yyz_yy_y_yy, g_yyz_yy_y_yz, g_yyz_yy_y_zz, g_z_0_y_0_yy_yy_0_xx, g_z_0_y_0_yy_yy_0_xy, g_z_0_y_0_yy_yy_0_xz, g_z_0_y_0_yy_yy_0_yy, g_z_0_y_0_yy_yy_0_yz, g_z_0_y_0_yy_yy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_yy_yy_0_xx[i] = 4.0 * g_yyz_yy_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_yy_0_xy[i] = 4.0 * g_yyz_yy_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_yy_0_xz[i] = 4.0 * g_yyz_yy_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_yy_0_yy[i] = 4.0 * g_yyz_yy_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_yy_0_yz[i] = 4.0 * g_yyz_yy_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_yy_0_zz[i] = 4.0 * g_yyz_yy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1644-1650)

    #pragma omp simd aligned(g_yyz_yz_y_xx, g_yyz_yz_y_xy, g_yyz_yz_y_xz, g_yyz_yz_y_yy, g_yyz_yz_y_yz, g_yyz_yz_y_zz, g_z_0_y_0_yy_yz_0_xx, g_z_0_y_0_yy_yz_0_xy, g_z_0_y_0_yy_yz_0_xz, g_z_0_y_0_yy_yz_0_yy, g_z_0_y_0_yy_yz_0_yz, g_z_0_y_0_yy_yz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_yy_yz_0_xx[i] = 4.0 * g_yyz_yz_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_yz_0_xy[i] = 4.0 * g_yyz_yz_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_yz_0_xz[i] = 4.0 * g_yyz_yz_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_yz_0_yy[i] = 4.0 * g_yyz_yz_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_yz_0_yz[i] = 4.0 * g_yyz_yz_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_yz_0_zz[i] = 4.0 * g_yyz_yz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1650-1656)

    #pragma omp simd aligned(g_yyz_zz_y_xx, g_yyz_zz_y_xy, g_yyz_zz_y_xz, g_yyz_zz_y_yy, g_yyz_zz_y_yz, g_yyz_zz_y_zz, g_z_0_y_0_yy_zz_0_xx, g_z_0_y_0_yy_zz_0_xy, g_z_0_y_0_yy_zz_0_xz, g_z_0_y_0_yy_zz_0_yy, g_z_0_y_0_yy_zz_0_yz, g_z_0_y_0_yy_zz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_yy_zz_0_xx[i] = 4.0 * g_yyz_zz_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_zz_0_xy[i] = 4.0 * g_yyz_zz_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_zz_0_xz[i] = 4.0 * g_yyz_zz_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_zz_0_yy[i] = 4.0 * g_yyz_zz_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_zz_0_yz[i] = 4.0 * g_yyz_zz_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_zz_0_zz[i] = 4.0 * g_yyz_zz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1656-1662)

    #pragma omp simd aligned(g_y_xx_y_xx, g_y_xx_y_xy, g_y_xx_y_xz, g_y_xx_y_yy, g_y_xx_y_yz, g_y_xx_y_zz, g_yzz_xx_y_xx, g_yzz_xx_y_xy, g_yzz_xx_y_xz, g_yzz_xx_y_yy, g_yzz_xx_y_yz, g_yzz_xx_y_zz, g_z_0_y_0_yz_xx_0_xx, g_z_0_y_0_yz_xx_0_xy, g_z_0_y_0_yz_xx_0_xz, g_z_0_y_0_yz_xx_0_yy, g_z_0_y_0_yz_xx_0_yz, g_z_0_y_0_yz_xx_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_yz_xx_0_xx[i] = -2.0 * g_y_xx_y_xx[i] * c_exps[i] + 4.0 * g_yzz_xx_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_xx_0_xy[i] = -2.0 * g_y_xx_y_xy[i] * c_exps[i] + 4.0 * g_yzz_xx_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_xx_0_xz[i] = -2.0 * g_y_xx_y_xz[i] * c_exps[i] + 4.0 * g_yzz_xx_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_xx_0_yy[i] = -2.0 * g_y_xx_y_yy[i] * c_exps[i] + 4.0 * g_yzz_xx_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_xx_0_yz[i] = -2.0 * g_y_xx_y_yz[i] * c_exps[i] + 4.0 * g_yzz_xx_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_xx_0_zz[i] = -2.0 * g_y_xx_y_zz[i] * c_exps[i] + 4.0 * g_yzz_xx_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1662-1668)

    #pragma omp simd aligned(g_y_xy_y_xx, g_y_xy_y_xy, g_y_xy_y_xz, g_y_xy_y_yy, g_y_xy_y_yz, g_y_xy_y_zz, g_yzz_xy_y_xx, g_yzz_xy_y_xy, g_yzz_xy_y_xz, g_yzz_xy_y_yy, g_yzz_xy_y_yz, g_yzz_xy_y_zz, g_z_0_y_0_yz_xy_0_xx, g_z_0_y_0_yz_xy_0_xy, g_z_0_y_0_yz_xy_0_xz, g_z_0_y_0_yz_xy_0_yy, g_z_0_y_0_yz_xy_0_yz, g_z_0_y_0_yz_xy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_yz_xy_0_xx[i] = -2.0 * g_y_xy_y_xx[i] * c_exps[i] + 4.0 * g_yzz_xy_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_xy_0_xy[i] = -2.0 * g_y_xy_y_xy[i] * c_exps[i] + 4.0 * g_yzz_xy_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_xy_0_xz[i] = -2.0 * g_y_xy_y_xz[i] * c_exps[i] + 4.0 * g_yzz_xy_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_xy_0_yy[i] = -2.0 * g_y_xy_y_yy[i] * c_exps[i] + 4.0 * g_yzz_xy_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_xy_0_yz[i] = -2.0 * g_y_xy_y_yz[i] * c_exps[i] + 4.0 * g_yzz_xy_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_xy_0_zz[i] = -2.0 * g_y_xy_y_zz[i] * c_exps[i] + 4.0 * g_yzz_xy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1668-1674)

    #pragma omp simd aligned(g_y_xz_y_xx, g_y_xz_y_xy, g_y_xz_y_xz, g_y_xz_y_yy, g_y_xz_y_yz, g_y_xz_y_zz, g_yzz_xz_y_xx, g_yzz_xz_y_xy, g_yzz_xz_y_xz, g_yzz_xz_y_yy, g_yzz_xz_y_yz, g_yzz_xz_y_zz, g_z_0_y_0_yz_xz_0_xx, g_z_0_y_0_yz_xz_0_xy, g_z_0_y_0_yz_xz_0_xz, g_z_0_y_0_yz_xz_0_yy, g_z_0_y_0_yz_xz_0_yz, g_z_0_y_0_yz_xz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_yz_xz_0_xx[i] = -2.0 * g_y_xz_y_xx[i] * c_exps[i] + 4.0 * g_yzz_xz_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_xz_0_xy[i] = -2.0 * g_y_xz_y_xy[i] * c_exps[i] + 4.0 * g_yzz_xz_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_xz_0_xz[i] = -2.0 * g_y_xz_y_xz[i] * c_exps[i] + 4.0 * g_yzz_xz_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_xz_0_yy[i] = -2.0 * g_y_xz_y_yy[i] * c_exps[i] + 4.0 * g_yzz_xz_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_xz_0_yz[i] = -2.0 * g_y_xz_y_yz[i] * c_exps[i] + 4.0 * g_yzz_xz_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_xz_0_zz[i] = -2.0 * g_y_xz_y_zz[i] * c_exps[i] + 4.0 * g_yzz_xz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1674-1680)

    #pragma omp simd aligned(g_y_yy_y_xx, g_y_yy_y_xy, g_y_yy_y_xz, g_y_yy_y_yy, g_y_yy_y_yz, g_y_yy_y_zz, g_yzz_yy_y_xx, g_yzz_yy_y_xy, g_yzz_yy_y_xz, g_yzz_yy_y_yy, g_yzz_yy_y_yz, g_yzz_yy_y_zz, g_z_0_y_0_yz_yy_0_xx, g_z_0_y_0_yz_yy_0_xy, g_z_0_y_0_yz_yy_0_xz, g_z_0_y_0_yz_yy_0_yy, g_z_0_y_0_yz_yy_0_yz, g_z_0_y_0_yz_yy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_yz_yy_0_xx[i] = -2.0 * g_y_yy_y_xx[i] * c_exps[i] + 4.0 * g_yzz_yy_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_yy_0_xy[i] = -2.0 * g_y_yy_y_xy[i] * c_exps[i] + 4.0 * g_yzz_yy_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_yy_0_xz[i] = -2.0 * g_y_yy_y_xz[i] * c_exps[i] + 4.0 * g_yzz_yy_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_yy_0_yy[i] = -2.0 * g_y_yy_y_yy[i] * c_exps[i] + 4.0 * g_yzz_yy_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_yy_0_yz[i] = -2.0 * g_y_yy_y_yz[i] * c_exps[i] + 4.0 * g_yzz_yy_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_yy_0_zz[i] = -2.0 * g_y_yy_y_zz[i] * c_exps[i] + 4.0 * g_yzz_yy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1680-1686)

    #pragma omp simd aligned(g_y_yz_y_xx, g_y_yz_y_xy, g_y_yz_y_xz, g_y_yz_y_yy, g_y_yz_y_yz, g_y_yz_y_zz, g_yzz_yz_y_xx, g_yzz_yz_y_xy, g_yzz_yz_y_xz, g_yzz_yz_y_yy, g_yzz_yz_y_yz, g_yzz_yz_y_zz, g_z_0_y_0_yz_yz_0_xx, g_z_0_y_0_yz_yz_0_xy, g_z_0_y_0_yz_yz_0_xz, g_z_0_y_0_yz_yz_0_yy, g_z_0_y_0_yz_yz_0_yz, g_z_0_y_0_yz_yz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_yz_yz_0_xx[i] = -2.0 * g_y_yz_y_xx[i] * c_exps[i] + 4.0 * g_yzz_yz_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_yz_0_xy[i] = -2.0 * g_y_yz_y_xy[i] * c_exps[i] + 4.0 * g_yzz_yz_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_yz_0_xz[i] = -2.0 * g_y_yz_y_xz[i] * c_exps[i] + 4.0 * g_yzz_yz_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_yz_0_yy[i] = -2.0 * g_y_yz_y_yy[i] * c_exps[i] + 4.0 * g_yzz_yz_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_yz_0_yz[i] = -2.0 * g_y_yz_y_yz[i] * c_exps[i] + 4.0 * g_yzz_yz_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_yz_0_zz[i] = -2.0 * g_y_yz_y_zz[i] * c_exps[i] + 4.0 * g_yzz_yz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1686-1692)

    #pragma omp simd aligned(g_y_zz_y_xx, g_y_zz_y_xy, g_y_zz_y_xz, g_y_zz_y_yy, g_y_zz_y_yz, g_y_zz_y_zz, g_yzz_zz_y_xx, g_yzz_zz_y_xy, g_yzz_zz_y_xz, g_yzz_zz_y_yy, g_yzz_zz_y_yz, g_yzz_zz_y_zz, g_z_0_y_0_yz_zz_0_xx, g_z_0_y_0_yz_zz_0_xy, g_z_0_y_0_yz_zz_0_xz, g_z_0_y_0_yz_zz_0_yy, g_z_0_y_0_yz_zz_0_yz, g_z_0_y_0_yz_zz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_yz_zz_0_xx[i] = -2.0 * g_y_zz_y_xx[i] * c_exps[i] + 4.0 * g_yzz_zz_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_zz_0_xy[i] = -2.0 * g_y_zz_y_xy[i] * c_exps[i] + 4.0 * g_yzz_zz_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_zz_0_xz[i] = -2.0 * g_y_zz_y_xz[i] * c_exps[i] + 4.0 * g_yzz_zz_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_zz_0_yy[i] = -2.0 * g_y_zz_y_yy[i] * c_exps[i] + 4.0 * g_yzz_zz_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_zz_0_yz[i] = -2.0 * g_y_zz_y_yz[i] * c_exps[i] + 4.0 * g_yzz_zz_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_zz_0_zz[i] = -2.0 * g_y_zz_y_zz[i] * c_exps[i] + 4.0 * g_yzz_zz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1692-1698)

    #pragma omp simd aligned(g_z_0_y_0_zz_xx_0_xx, g_z_0_y_0_zz_xx_0_xy, g_z_0_y_0_zz_xx_0_xz, g_z_0_y_0_zz_xx_0_yy, g_z_0_y_0_zz_xx_0_yz, g_z_0_y_0_zz_xx_0_zz, g_z_xx_y_xx, g_z_xx_y_xy, g_z_xx_y_xz, g_z_xx_y_yy, g_z_xx_y_yz, g_z_xx_y_zz, g_zzz_xx_y_xx, g_zzz_xx_y_xy, g_zzz_xx_y_xz, g_zzz_xx_y_yy, g_zzz_xx_y_yz, g_zzz_xx_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_zz_xx_0_xx[i] = -4.0 * g_z_xx_y_xx[i] * c_exps[i] + 4.0 * g_zzz_xx_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_xx_0_xy[i] = -4.0 * g_z_xx_y_xy[i] * c_exps[i] + 4.0 * g_zzz_xx_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_xx_0_xz[i] = -4.0 * g_z_xx_y_xz[i] * c_exps[i] + 4.0 * g_zzz_xx_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_xx_0_yy[i] = -4.0 * g_z_xx_y_yy[i] * c_exps[i] + 4.0 * g_zzz_xx_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_xx_0_yz[i] = -4.0 * g_z_xx_y_yz[i] * c_exps[i] + 4.0 * g_zzz_xx_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_xx_0_zz[i] = -4.0 * g_z_xx_y_zz[i] * c_exps[i] + 4.0 * g_zzz_xx_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1698-1704)

    #pragma omp simd aligned(g_z_0_y_0_zz_xy_0_xx, g_z_0_y_0_zz_xy_0_xy, g_z_0_y_0_zz_xy_0_xz, g_z_0_y_0_zz_xy_0_yy, g_z_0_y_0_zz_xy_0_yz, g_z_0_y_0_zz_xy_0_zz, g_z_xy_y_xx, g_z_xy_y_xy, g_z_xy_y_xz, g_z_xy_y_yy, g_z_xy_y_yz, g_z_xy_y_zz, g_zzz_xy_y_xx, g_zzz_xy_y_xy, g_zzz_xy_y_xz, g_zzz_xy_y_yy, g_zzz_xy_y_yz, g_zzz_xy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_zz_xy_0_xx[i] = -4.0 * g_z_xy_y_xx[i] * c_exps[i] + 4.0 * g_zzz_xy_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_xy_0_xy[i] = -4.0 * g_z_xy_y_xy[i] * c_exps[i] + 4.0 * g_zzz_xy_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_xy_0_xz[i] = -4.0 * g_z_xy_y_xz[i] * c_exps[i] + 4.0 * g_zzz_xy_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_xy_0_yy[i] = -4.0 * g_z_xy_y_yy[i] * c_exps[i] + 4.0 * g_zzz_xy_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_xy_0_yz[i] = -4.0 * g_z_xy_y_yz[i] * c_exps[i] + 4.0 * g_zzz_xy_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_xy_0_zz[i] = -4.0 * g_z_xy_y_zz[i] * c_exps[i] + 4.0 * g_zzz_xy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1704-1710)

    #pragma omp simd aligned(g_z_0_y_0_zz_xz_0_xx, g_z_0_y_0_zz_xz_0_xy, g_z_0_y_0_zz_xz_0_xz, g_z_0_y_0_zz_xz_0_yy, g_z_0_y_0_zz_xz_0_yz, g_z_0_y_0_zz_xz_0_zz, g_z_xz_y_xx, g_z_xz_y_xy, g_z_xz_y_xz, g_z_xz_y_yy, g_z_xz_y_yz, g_z_xz_y_zz, g_zzz_xz_y_xx, g_zzz_xz_y_xy, g_zzz_xz_y_xz, g_zzz_xz_y_yy, g_zzz_xz_y_yz, g_zzz_xz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_zz_xz_0_xx[i] = -4.0 * g_z_xz_y_xx[i] * c_exps[i] + 4.0 * g_zzz_xz_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_xz_0_xy[i] = -4.0 * g_z_xz_y_xy[i] * c_exps[i] + 4.0 * g_zzz_xz_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_xz_0_xz[i] = -4.0 * g_z_xz_y_xz[i] * c_exps[i] + 4.0 * g_zzz_xz_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_xz_0_yy[i] = -4.0 * g_z_xz_y_yy[i] * c_exps[i] + 4.0 * g_zzz_xz_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_xz_0_yz[i] = -4.0 * g_z_xz_y_yz[i] * c_exps[i] + 4.0 * g_zzz_xz_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_xz_0_zz[i] = -4.0 * g_z_xz_y_zz[i] * c_exps[i] + 4.0 * g_zzz_xz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1710-1716)

    #pragma omp simd aligned(g_z_0_y_0_zz_yy_0_xx, g_z_0_y_0_zz_yy_0_xy, g_z_0_y_0_zz_yy_0_xz, g_z_0_y_0_zz_yy_0_yy, g_z_0_y_0_zz_yy_0_yz, g_z_0_y_0_zz_yy_0_zz, g_z_yy_y_xx, g_z_yy_y_xy, g_z_yy_y_xz, g_z_yy_y_yy, g_z_yy_y_yz, g_z_yy_y_zz, g_zzz_yy_y_xx, g_zzz_yy_y_xy, g_zzz_yy_y_xz, g_zzz_yy_y_yy, g_zzz_yy_y_yz, g_zzz_yy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_zz_yy_0_xx[i] = -4.0 * g_z_yy_y_xx[i] * c_exps[i] + 4.0 * g_zzz_yy_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_yy_0_xy[i] = -4.0 * g_z_yy_y_xy[i] * c_exps[i] + 4.0 * g_zzz_yy_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_yy_0_xz[i] = -4.0 * g_z_yy_y_xz[i] * c_exps[i] + 4.0 * g_zzz_yy_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_yy_0_yy[i] = -4.0 * g_z_yy_y_yy[i] * c_exps[i] + 4.0 * g_zzz_yy_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_yy_0_yz[i] = -4.0 * g_z_yy_y_yz[i] * c_exps[i] + 4.0 * g_zzz_yy_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_yy_0_zz[i] = -4.0 * g_z_yy_y_zz[i] * c_exps[i] + 4.0 * g_zzz_yy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1716-1722)

    #pragma omp simd aligned(g_z_0_y_0_zz_yz_0_xx, g_z_0_y_0_zz_yz_0_xy, g_z_0_y_0_zz_yz_0_xz, g_z_0_y_0_zz_yz_0_yy, g_z_0_y_0_zz_yz_0_yz, g_z_0_y_0_zz_yz_0_zz, g_z_yz_y_xx, g_z_yz_y_xy, g_z_yz_y_xz, g_z_yz_y_yy, g_z_yz_y_yz, g_z_yz_y_zz, g_zzz_yz_y_xx, g_zzz_yz_y_xy, g_zzz_yz_y_xz, g_zzz_yz_y_yy, g_zzz_yz_y_yz, g_zzz_yz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_zz_yz_0_xx[i] = -4.0 * g_z_yz_y_xx[i] * c_exps[i] + 4.0 * g_zzz_yz_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_yz_0_xy[i] = -4.0 * g_z_yz_y_xy[i] * c_exps[i] + 4.0 * g_zzz_yz_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_yz_0_xz[i] = -4.0 * g_z_yz_y_xz[i] * c_exps[i] + 4.0 * g_zzz_yz_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_yz_0_yy[i] = -4.0 * g_z_yz_y_yy[i] * c_exps[i] + 4.0 * g_zzz_yz_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_yz_0_yz[i] = -4.0 * g_z_yz_y_yz[i] * c_exps[i] + 4.0 * g_zzz_yz_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_yz_0_zz[i] = -4.0 * g_z_yz_y_zz[i] * c_exps[i] + 4.0 * g_zzz_yz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1722-1728)

    #pragma omp simd aligned(g_z_0_y_0_zz_zz_0_xx, g_z_0_y_0_zz_zz_0_xy, g_z_0_y_0_zz_zz_0_xz, g_z_0_y_0_zz_zz_0_yy, g_z_0_y_0_zz_zz_0_yz, g_z_0_y_0_zz_zz_0_zz, g_z_zz_y_xx, g_z_zz_y_xy, g_z_zz_y_xz, g_z_zz_y_yy, g_z_zz_y_yz, g_z_zz_y_zz, g_zzz_zz_y_xx, g_zzz_zz_y_xy, g_zzz_zz_y_xz, g_zzz_zz_y_yy, g_zzz_zz_y_yz, g_zzz_zz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_zz_zz_0_xx[i] = -4.0 * g_z_zz_y_xx[i] * c_exps[i] + 4.0 * g_zzz_zz_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_zz_0_xy[i] = -4.0 * g_z_zz_y_xy[i] * c_exps[i] + 4.0 * g_zzz_zz_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_zz_0_xz[i] = -4.0 * g_z_zz_y_xz[i] * c_exps[i] + 4.0 * g_zzz_zz_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_zz_0_yy[i] = -4.0 * g_z_zz_y_yy[i] * c_exps[i] + 4.0 * g_zzz_zz_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_zz_0_yz[i] = -4.0 * g_z_zz_y_yz[i] * c_exps[i] + 4.0 * g_zzz_zz_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_zz_0_zz[i] = -4.0 * g_z_zz_y_zz[i] * c_exps[i] + 4.0 * g_zzz_zz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1728-1734)

    #pragma omp simd aligned(g_xxz_xx_z_xx, g_xxz_xx_z_xy, g_xxz_xx_z_xz, g_xxz_xx_z_yy, g_xxz_xx_z_yz, g_xxz_xx_z_zz, g_z_0_z_0_xx_xx_0_xx, g_z_0_z_0_xx_xx_0_xy, g_z_0_z_0_xx_xx_0_xz, g_z_0_z_0_xx_xx_0_yy, g_z_0_z_0_xx_xx_0_yz, g_z_0_z_0_xx_xx_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_xx_xx_0_xx[i] = 4.0 * g_xxz_xx_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_xx_0_xy[i] = 4.0 * g_xxz_xx_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_xx_0_xz[i] = 4.0 * g_xxz_xx_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_xx_0_yy[i] = 4.0 * g_xxz_xx_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_xx_0_yz[i] = 4.0 * g_xxz_xx_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_xx_0_zz[i] = 4.0 * g_xxz_xx_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1734-1740)

    #pragma omp simd aligned(g_xxz_xy_z_xx, g_xxz_xy_z_xy, g_xxz_xy_z_xz, g_xxz_xy_z_yy, g_xxz_xy_z_yz, g_xxz_xy_z_zz, g_z_0_z_0_xx_xy_0_xx, g_z_0_z_0_xx_xy_0_xy, g_z_0_z_0_xx_xy_0_xz, g_z_0_z_0_xx_xy_0_yy, g_z_0_z_0_xx_xy_0_yz, g_z_0_z_0_xx_xy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_xx_xy_0_xx[i] = 4.0 * g_xxz_xy_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_xy_0_xy[i] = 4.0 * g_xxz_xy_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_xy_0_xz[i] = 4.0 * g_xxz_xy_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_xy_0_yy[i] = 4.0 * g_xxz_xy_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_xy_0_yz[i] = 4.0 * g_xxz_xy_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_xy_0_zz[i] = 4.0 * g_xxz_xy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1740-1746)

    #pragma omp simd aligned(g_xxz_xz_z_xx, g_xxz_xz_z_xy, g_xxz_xz_z_xz, g_xxz_xz_z_yy, g_xxz_xz_z_yz, g_xxz_xz_z_zz, g_z_0_z_0_xx_xz_0_xx, g_z_0_z_0_xx_xz_0_xy, g_z_0_z_0_xx_xz_0_xz, g_z_0_z_0_xx_xz_0_yy, g_z_0_z_0_xx_xz_0_yz, g_z_0_z_0_xx_xz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_xx_xz_0_xx[i] = 4.0 * g_xxz_xz_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_xz_0_xy[i] = 4.0 * g_xxz_xz_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_xz_0_xz[i] = 4.0 * g_xxz_xz_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_xz_0_yy[i] = 4.0 * g_xxz_xz_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_xz_0_yz[i] = 4.0 * g_xxz_xz_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_xz_0_zz[i] = 4.0 * g_xxz_xz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1746-1752)

    #pragma omp simd aligned(g_xxz_yy_z_xx, g_xxz_yy_z_xy, g_xxz_yy_z_xz, g_xxz_yy_z_yy, g_xxz_yy_z_yz, g_xxz_yy_z_zz, g_z_0_z_0_xx_yy_0_xx, g_z_0_z_0_xx_yy_0_xy, g_z_0_z_0_xx_yy_0_xz, g_z_0_z_0_xx_yy_0_yy, g_z_0_z_0_xx_yy_0_yz, g_z_0_z_0_xx_yy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_xx_yy_0_xx[i] = 4.0 * g_xxz_yy_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_yy_0_xy[i] = 4.0 * g_xxz_yy_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_yy_0_xz[i] = 4.0 * g_xxz_yy_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_yy_0_yy[i] = 4.0 * g_xxz_yy_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_yy_0_yz[i] = 4.0 * g_xxz_yy_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_yy_0_zz[i] = 4.0 * g_xxz_yy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1752-1758)

    #pragma omp simd aligned(g_xxz_yz_z_xx, g_xxz_yz_z_xy, g_xxz_yz_z_xz, g_xxz_yz_z_yy, g_xxz_yz_z_yz, g_xxz_yz_z_zz, g_z_0_z_0_xx_yz_0_xx, g_z_0_z_0_xx_yz_0_xy, g_z_0_z_0_xx_yz_0_xz, g_z_0_z_0_xx_yz_0_yy, g_z_0_z_0_xx_yz_0_yz, g_z_0_z_0_xx_yz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_xx_yz_0_xx[i] = 4.0 * g_xxz_yz_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_yz_0_xy[i] = 4.0 * g_xxz_yz_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_yz_0_xz[i] = 4.0 * g_xxz_yz_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_yz_0_yy[i] = 4.0 * g_xxz_yz_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_yz_0_yz[i] = 4.0 * g_xxz_yz_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_yz_0_zz[i] = 4.0 * g_xxz_yz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1758-1764)

    #pragma omp simd aligned(g_xxz_zz_z_xx, g_xxz_zz_z_xy, g_xxz_zz_z_xz, g_xxz_zz_z_yy, g_xxz_zz_z_yz, g_xxz_zz_z_zz, g_z_0_z_0_xx_zz_0_xx, g_z_0_z_0_xx_zz_0_xy, g_z_0_z_0_xx_zz_0_xz, g_z_0_z_0_xx_zz_0_yy, g_z_0_z_0_xx_zz_0_yz, g_z_0_z_0_xx_zz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_xx_zz_0_xx[i] = 4.0 * g_xxz_zz_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_zz_0_xy[i] = 4.0 * g_xxz_zz_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_zz_0_xz[i] = 4.0 * g_xxz_zz_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_zz_0_yy[i] = 4.0 * g_xxz_zz_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_zz_0_yz[i] = 4.0 * g_xxz_zz_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_zz_0_zz[i] = 4.0 * g_xxz_zz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1764-1770)

    #pragma omp simd aligned(g_xyz_xx_z_xx, g_xyz_xx_z_xy, g_xyz_xx_z_xz, g_xyz_xx_z_yy, g_xyz_xx_z_yz, g_xyz_xx_z_zz, g_z_0_z_0_xy_xx_0_xx, g_z_0_z_0_xy_xx_0_xy, g_z_0_z_0_xy_xx_0_xz, g_z_0_z_0_xy_xx_0_yy, g_z_0_z_0_xy_xx_0_yz, g_z_0_z_0_xy_xx_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_xy_xx_0_xx[i] = 4.0 * g_xyz_xx_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_xx_0_xy[i] = 4.0 * g_xyz_xx_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_xx_0_xz[i] = 4.0 * g_xyz_xx_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_xx_0_yy[i] = 4.0 * g_xyz_xx_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_xx_0_yz[i] = 4.0 * g_xyz_xx_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_xx_0_zz[i] = 4.0 * g_xyz_xx_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1770-1776)

    #pragma omp simd aligned(g_xyz_xy_z_xx, g_xyz_xy_z_xy, g_xyz_xy_z_xz, g_xyz_xy_z_yy, g_xyz_xy_z_yz, g_xyz_xy_z_zz, g_z_0_z_0_xy_xy_0_xx, g_z_0_z_0_xy_xy_0_xy, g_z_0_z_0_xy_xy_0_xz, g_z_0_z_0_xy_xy_0_yy, g_z_0_z_0_xy_xy_0_yz, g_z_0_z_0_xy_xy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_xy_xy_0_xx[i] = 4.0 * g_xyz_xy_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_xy_0_xy[i] = 4.0 * g_xyz_xy_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_xy_0_xz[i] = 4.0 * g_xyz_xy_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_xy_0_yy[i] = 4.0 * g_xyz_xy_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_xy_0_yz[i] = 4.0 * g_xyz_xy_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_xy_0_zz[i] = 4.0 * g_xyz_xy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1776-1782)

    #pragma omp simd aligned(g_xyz_xz_z_xx, g_xyz_xz_z_xy, g_xyz_xz_z_xz, g_xyz_xz_z_yy, g_xyz_xz_z_yz, g_xyz_xz_z_zz, g_z_0_z_0_xy_xz_0_xx, g_z_0_z_0_xy_xz_0_xy, g_z_0_z_0_xy_xz_0_xz, g_z_0_z_0_xy_xz_0_yy, g_z_0_z_0_xy_xz_0_yz, g_z_0_z_0_xy_xz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_xy_xz_0_xx[i] = 4.0 * g_xyz_xz_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_xz_0_xy[i] = 4.0 * g_xyz_xz_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_xz_0_xz[i] = 4.0 * g_xyz_xz_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_xz_0_yy[i] = 4.0 * g_xyz_xz_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_xz_0_yz[i] = 4.0 * g_xyz_xz_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_xz_0_zz[i] = 4.0 * g_xyz_xz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1782-1788)

    #pragma omp simd aligned(g_xyz_yy_z_xx, g_xyz_yy_z_xy, g_xyz_yy_z_xz, g_xyz_yy_z_yy, g_xyz_yy_z_yz, g_xyz_yy_z_zz, g_z_0_z_0_xy_yy_0_xx, g_z_0_z_0_xy_yy_0_xy, g_z_0_z_0_xy_yy_0_xz, g_z_0_z_0_xy_yy_0_yy, g_z_0_z_0_xy_yy_0_yz, g_z_0_z_0_xy_yy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_xy_yy_0_xx[i] = 4.0 * g_xyz_yy_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_yy_0_xy[i] = 4.0 * g_xyz_yy_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_yy_0_xz[i] = 4.0 * g_xyz_yy_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_yy_0_yy[i] = 4.0 * g_xyz_yy_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_yy_0_yz[i] = 4.0 * g_xyz_yy_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_yy_0_zz[i] = 4.0 * g_xyz_yy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1788-1794)

    #pragma omp simd aligned(g_xyz_yz_z_xx, g_xyz_yz_z_xy, g_xyz_yz_z_xz, g_xyz_yz_z_yy, g_xyz_yz_z_yz, g_xyz_yz_z_zz, g_z_0_z_0_xy_yz_0_xx, g_z_0_z_0_xy_yz_0_xy, g_z_0_z_0_xy_yz_0_xz, g_z_0_z_0_xy_yz_0_yy, g_z_0_z_0_xy_yz_0_yz, g_z_0_z_0_xy_yz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_xy_yz_0_xx[i] = 4.0 * g_xyz_yz_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_yz_0_xy[i] = 4.0 * g_xyz_yz_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_yz_0_xz[i] = 4.0 * g_xyz_yz_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_yz_0_yy[i] = 4.0 * g_xyz_yz_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_yz_0_yz[i] = 4.0 * g_xyz_yz_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_yz_0_zz[i] = 4.0 * g_xyz_yz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1794-1800)

    #pragma omp simd aligned(g_xyz_zz_z_xx, g_xyz_zz_z_xy, g_xyz_zz_z_xz, g_xyz_zz_z_yy, g_xyz_zz_z_yz, g_xyz_zz_z_zz, g_z_0_z_0_xy_zz_0_xx, g_z_0_z_0_xy_zz_0_xy, g_z_0_z_0_xy_zz_0_xz, g_z_0_z_0_xy_zz_0_yy, g_z_0_z_0_xy_zz_0_yz, g_z_0_z_0_xy_zz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_xy_zz_0_xx[i] = 4.0 * g_xyz_zz_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_zz_0_xy[i] = 4.0 * g_xyz_zz_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_zz_0_xz[i] = 4.0 * g_xyz_zz_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_zz_0_yy[i] = 4.0 * g_xyz_zz_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_zz_0_yz[i] = 4.0 * g_xyz_zz_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_zz_0_zz[i] = 4.0 * g_xyz_zz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1800-1806)

    #pragma omp simd aligned(g_x_xx_z_xx, g_x_xx_z_xy, g_x_xx_z_xz, g_x_xx_z_yy, g_x_xx_z_yz, g_x_xx_z_zz, g_xzz_xx_z_xx, g_xzz_xx_z_xy, g_xzz_xx_z_xz, g_xzz_xx_z_yy, g_xzz_xx_z_yz, g_xzz_xx_z_zz, g_z_0_z_0_xz_xx_0_xx, g_z_0_z_0_xz_xx_0_xy, g_z_0_z_0_xz_xx_0_xz, g_z_0_z_0_xz_xx_0_yy, g_z_0_z_0_xz_xx_0_yz, g_z_0_z_0_xz_xx_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_xz_xx_0_xx[i] = -2.0 * g_x_xx_z_xx[i] * c_exps[i] + 4.0 * g_xzz_xx_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_xx_0_xy[i] = -2.0 * g_x_xx_z_xy[i] * c_exps[i] + 4.0 * g_xzz_xx_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_xx_0_xz[i] = -2.0 * g_x_xx_z_xz[i] * c_exps[i] + 4.0 * g_xzz_xx_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_xx_0_yy[i] = -2.0 * g_x_xx_z_yy[i] * c_exps[i] + 4.0 * g_xzz_xx_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_xx_0_yz[i] = -2.0 * g_x_xx_z_yz[i] * c_exps[i] + 4.0 * g_xzz_xx_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_xx_0_zz[i] = -2.0 * g_x_xx_z_zz[i] * c_exps[i] + 4.0 * g_xzz_xx_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1806-1812)

    #pragma omp simd aligned(g_x_xy_z_xx, g_x_xy_z_xy, g_x_xy_z_xz, g_x_xy_z_yy, g_x_xy_z_yz, g_x_xy_z_zz, g_xzz_xy_z_xx, g_xzz_xy_z_xy, g_xzz_xy_z_xz, g_xzz_xy_z_yy, g_xzz_xy_z_yz, g_xzz_xy_z_zz, g_z_0_z_0_xz_xy_0_xx, g_z_0_z_0_xz_xy_0_xy, g_z_0_z_0_xz_xy_0_xz, g_z_0_z_0_xz_xy_0_yy, g_z_0_z_0_xz_xy_0_yz, g_z_0_z_0_xz_xy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_xz_xy_0_xx[i] = -2.0 * g_x_xy_z_xx[i] * c_exps[i] + 4.0 * g_xzz_xy_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_xy_0_xy[i] = -2.0 * g_x_xy_z_xy[i] * c_exps[i] + 4.0 * g_xzz_xy_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_xy_0_xz[i] = -2.0 * g_x_xy_z_xz[i] * c_exps[i] + 4.0 * g_xzz_xy_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_xy_0_yy[i] = -2.0 * g_x_xy_z_yy[i] * c_exps[i] + 4.0 * g_xzz_xy_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_xy_0_yz[i] = -2.0 * g_x_xy_z_yz[i] * c_exps[i] + 4.0 * g_xzz_xy_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_xy_0_zz[i] = -2.0 * g_x_xy_z_zz[i] * c_exps[i] + 4.0 * g_xzz_xy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1812-1818)

    #pragma omp simd aligned(g_x_xz_z_xx, g_x_xz_z_xy, g_x_xz_z_xz, g_x_xz_z_yy, g_x_xz_z_yz, g_x_xz_z_zz, g_xzz_xz_z_xx, g_xzz_xz_z_xy, g_xzz_xz_z_xz, g_xzz_xz_z_yy, g_xzz_xz_z_yz, g_xzz_xz_z_zz, g_z_0_z_0_xz_xz_0_xx, g_z_0_z_0_xz_xz_0_xy, g_z_0_z_0_xz_xz_0_xz, g_z_0_z_0_xz_xz_0_yy, g_z_0_z_0_xz_xz_0_yz, g_z_0_z_0_xz_xz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_xz_xz_0_xx[i] = -2.0 * g_x_xz_z_xx[i] * c_exps[i] + 4.0 * g_xzz_xz_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_xz_0_xy[i] = -2.0 * g_x_xz_z_xy[i] * c_exps[i] + 4.0 * g_xzz_xz_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_xz_0_xz[i] = -2.0 * g_x_xz_z_xz[i] * c_exps[i] + 4.0 * g_xzz_xz_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_xz_0_yy[i] = -2.0 * g_x_xz_z_yy[i] * c_exps[i] + 4.0 * g_xzz_xz_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_xz_0_yz[i] = -2.0 * g_x_xz_z_yz[i] * c_exps[i] + 4.0 * g_xzz_xz_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_xz_0_zz[i] = -2.0 * g_x_xz_z_zz[i] * c_exps[i] + 4.0 * g_xzz_xz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1818-1824)

    #pragma omp simd aligned(g_x_yy_z_xx, g_x_yy_z_xy, g_x_yy_z_xz, g_x_yy_z_yy, g_x_yy_z_yz, g_x_yy_z_zz, g_xzz_yy_z_xx, g_xzz_yy_z_xy, g_xzz_yy_z_xz, g_xzz_yy_z_yy, g_xzz_yy_z_yz, g_xzz_yy_z_zz, g_z_0_z_0_xz_yy_0_xx, g_z_0_z_0_xz_yy_0_xy, g_z_0_z_0_xz_yy_0_xz, g_z_0_z_0_xz_yy_0_yy, g_z_0_z_0_xz_yy_0_yz, g_z_0_z_0_xz_yy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_xz_yy_0_xx[i] = -2.0 * g_x_yy_z_xx[i] * c_exps[i] + 4.0 * g_xzz_yy_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_yy_0_xy[i] = -2.0 * g_x_yy_z_xy[i] * c_exps[i] + 4.0 * g_xzz_yy_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_yy_0_xz[i] = -2.0 * g_x_yy_z_xz[i] * c_exps[i] + 4.0 * g_xzz_yy_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_yy_0_yy[i] = -2.0 * g_x_yy_z_yy[i] * c_exps[i] + 4.0 * g_xzz_yy_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_yy_0_yz[i] = -2.0 * g_x_yy_z_yz[i] * c_exps[i] + 4.0 * g_xzz_yy_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_yy_0_zz[i] = -2.0 * g_x_yy_z_zz[i] * c_exps[i] + 4.0 * g_xzz_yy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1824-1830)

    #pragma omp simd aligned(g_x_yz_z_xx, g_x_yz_z_xy, g_x_yz_z_xz, g_x_yz_z_yy, g_x_yz_z_yz, g_x_yz_z_zz, g_xzz_yz_z_xx, g_xzz_yz_z_xy, g_xzz_yz_z_xz, g_xzz_yz_z_yy, g_xzz_yz_z_yz, g_xzz_yz_z_zz, g_z_0_z_0_xz_yz_0_xx, g_z_0_z_0_xz_yz_0_xy, g_z_0_z_0_xz_yz_0_xz, g_z_0_z_0_xz_yz_0_yy, g_z_0_z_0_xz_yz_0_yz, g_z_0_z_0_xz_yz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_xz_yz_0_xx[i] = -2.0 * g_x_yz_z_xx[i] * c_exps[i] + 4.0 * g_xzz_yz_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_yz_0_xy[i] = -2.0 * g_x_yz_z_xy[i] * c_exps[i] + 4.0 * g_xzz_yz_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_yz_0_xz[i] = -2.0 * g_x_yz_z_xz[i] * c_exps[i] + 4.0 * g_xzz_yz_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_yz_0_yy[i] = -2.0 * g_x_yz_z_yy[i] * c_exps[i] + 4.0 * g_xzz_yz_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_yz_0_yz[i] = -2.0 * g_x_yz_z_yz[i] * c_exps[i] + 4.0 * g_xzz_yz_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_yz_0_zz[i] = -2.0 * g_x_yz_z_zz[i] * c_exps[i] + 4.0 * g_xzz_yz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1830-1836)

    #pragma omp simd aligned(g_x_zz_z_xx, g_x_zz_z_xy, g_x_zz_z_xz, g_x_zz_z_yy, g_x_zz_z_yz, g_x_zz_z_zz, g_xzz_zz_z_xx, g_xzz_zz_z_xy, g_xzz_zz_z_xz, g_xzz_zz_z_yy, g_xzz_zz_z_yz, g_xzz_zz_z_zz, g_z_0_z_0_xz_zz_0_xx, g_z_0_z_0_xz_zz_0_xy, g_z_0_z_0_xz_zz_0_xz, g_z_0_z_0_xz_zz_0_yy, g_z_0_z_0_xz_zz_0_yz, g_z_0_z_0_xz_zz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_xz_zz_0_xx[i] = -2.0 * g_x_zz_z_xx[i] * c_exps[i] + 4.0 * g_xzz_zz_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_zz_0_xy[i] = -2.0 * g_x_zz_z_xy[i] * c_exps[i] + 4.0 * g_xzz_zz_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_zz_0_xz[i] = -2.0 * g_x_zz_z_xz[i] * c_exps[i] + 4.0 * g_xzz_zz_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_zz_0_yy[i] = -2.0 * g_x_zz_z_yy[i] * c_exps[i] + 4.0 * g_xzz_zz_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_zz_0_yz[i] = -2.0 * g_x_zz_z_yz[i] * c_exps[i] + 4.0 * g_xzz_zz_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_zz_0_zz[i] = -2.0 * g_x_zz_z_zz[i] * c_exps[i] + 4.0 * g_xzz_zz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1836-1842)

    #pragma omp simd aligned(g_yyz_xx_z_xx, g_yyz_xx_z_xy, g_yyz_xx_z_xz, g_yyz_xx_z_yy, g_yyz_xx_z_yz, g_yyz_xx_z_zz, g_z_0_z_0_yy_xx_0_xx, g_z_0_z_0_yy_xx_0_xy, g_z_0_z_0_yy_xx_0_xz, g_z_0_z_0_yy_xx_0_yy, g_z_0_z_0_yy_xx_0_yz, g_z_0_z_0_yy_xx_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_yy_xx_0_xx[i] = 4.0 * g_yyz_xx_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_xx_0_xy[i] = 4.0 * g_yyz_xx_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_xx_0_xz[i] = 4.0 * g_yyz_xx_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_xx_0_yy[i] = 4.0 * g_yyz_xx_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_xx_0_yz[i] = 4.0 * g_yyz_xx_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_xx_0_zz[i] = 4.0 * g_yyz_xx_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1842-1848)

    #pragma omp simd aligned(g_yyz_xy_z_xx, g_yyz_xy_z_xy, g_yyz_xy_z_xz, g_yyz_xy_z_yy, g_yyz_xy_z_yz, g_yyz_xy_z_zz, g_z_0_z_0_yy_xy_0_xx, g_z_0_z_0_yy_xy_0_xy, g_z_0_z_0_yy_xy_0_xz, g_z_0_z_0_yy_xy_0_yy, g_z_0_z_0_yy_xy_0_yz, g_z_0_z_0_yy_xy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_yy_xy_0_xx[i] = 4.0 * g_yyz_xy_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_xy_0_xy[i] = 4.0 * g_yyz_xy_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_xy_0_xz[i] = 4.0 * g_yyz_xy_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_xy_0_yy[i] = 4.0 * g_yyz_xy_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_xy_0_yz[i] = 4.0 * g_yyz_xy_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_xy_0_zz[i] = 4.0 * g_yyz_xy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1848-1854)

    #pragma omp simd aligned(g_yyz_xz_z_xx, g_yyz_xz_z_xy, g_yyz_xz_z_xz, g_yyz_xz_z_yy, g_yyz_xz_z_yz, g_yyz_xz_z_zz, g_z_0_z_0_yy_xz_0_xx, g_z_0_z_0_yy_xz_0_xy, g_z_0_z_0_yy_xz_0_xz, g_z_0_z_0_yy_xz_0_yy, g_z_0_z_0_yy_xz_0_yz, g_z_0_z_0_yy_xz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_yy_xz_0_xx[i] = 4.0 * g_yyz_xz_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_xz_0_xy[i] = 4.0 * g_yyz_xz_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_xz_0_xz[i] = 4.0 * g_yyz_xz_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_xz_0_yy[i] = 4.0 * g_yyz_xz_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_xz_0_yz[i] = 4.0 * g_yyz_xz_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_xz_0_zz[i] = 4.0 * g_yyz_xz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1854-1860)

    #pragma omp simd aligned(g_yyz_yy_z_xx, g_yyz_yy_z_xy, g_yyz_yy_z_xz, g_yyz_yy_z_yy, g_yyz_yy_z_yz, g_yyz_yy_z_zz, g_z_0_z_0_yy_yy_0_xx, g_z_0_z_0_yy_yy_0_xy, g_z_0_z_0_yy_yy_0_xz, g_z_0_z_0_yy_yy_0_yy, g_z_0_z_0_yy_yy_0_yz, g_z_0_z_0_yy_yy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_yy_yy_0_xx[i] = 4.0 * g_yyz_yy_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_yy_0_xy[i] = 4.0 * g_yyz_yy_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_yy_0_xz[i] = 4.0 * g_yyz_yy_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_yy_0_yy[i] = 4.0 * g_yyz_yy_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_yy_0_yz[i] = 4.0 * g_yyz_yy_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_yy_0_zz[i] = 4.0 * g_yyz_yy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1860-1866)

    #pragma omp simd aligned(g_yyz_yz_z_xx, g_yyz_yz_z_xy, g_yyz_yz_z_xz, g_yyz_yz_z_yy, g_yyz_yz_z_yz, g_yyz_yz_z_zz, g_z_0_z_0_yy_yz_0_xx, g_z_0_z_0_yy_yz_0_xy, g_z_0_z_0_yy_yz_0_xz, g_z_0_z_0_yy_yz_0_yy, g_z_0_z_0_yy_yz_0_yz, g_z_0_z_0_yy_yz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_yy_yz_0_xx[i] = 4.0 * g_yyz_yz_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_yz_0_xy[i] = 4.0 * g_yyz_yz_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_yz_0_xz[i] = 4.0 * g_yyz_yz_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_yz_0_yy[i] = 4.0 * g_yyz_yz_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_yz_0_yz[i] = 4.0 * g_yyz_yz_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_yz_0_zz[i] = 4.0 * g_yyz_yz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1866-1872)

    #pragma omp simd aligned(g_yyz_zz_z_xx, g_yyz_zz_z_xy, g_yyz_zz_z_xz, g_yyz_zz_z_yy, g_yyz_zz_z_yz, g_yyz_zz_z_zz, g_z_0_z_0_yy_zz_0_xx, g_z_0_z_0_yy_zz_0_xy, g_z_0_z_0_yy_zz_0_xz, g_z_0_z_0_yy_zz_0_yy, g_z_0_z_0_yy_zz_0_yz, g_z_0_z_0_yy_zz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_yy_zz_0_xx[i] = 4.0 * g_yyz_zz_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_zz_0_xy[i] = 4.0 * g_yyz_zz_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_zz_0_xz[i] = 4.0 * g_yyz_zz_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_zz_0_yy[i] = 4.0 * g_yyz_zz_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_zz_0_yz[i] = 4.0 * g_yyz_zz_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_zz_0_zz[i] = 4.0 * g_yyz_zz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1872-1878)

    #pragma omp simd aligned(g_y_xx_z_xx, g_y_xx_z_xy, g_y_xx_z_xz, g_y_xx_z_yy, g_y_xx_z_yz, g_y_xx_z_zz, g_yzz_xx_z_xx, g_yzz_xx_z_xy, g_yzz_xx_z_xz, g_yzz_xx_z_yy, g_yzz_xx_z_yz, g_yzz_xx_z_zz, g_z_0_z_0_yz_xx_0_xx, g_z_0_z_0_yz_xx_0_xy, g_z_0_z_0_yz_xx_0_xz, g_z_0_z_0_yz_xx_0_yy, g_z_0_z_0_yz_xx_0_yz, g_z_0_z_0_yz_xx_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_yz_xx_0_xx[i] = -2.0 * g_y_xx_z_xx[i] * c_exps[i] + 4.0 * g_yzz_xx_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_xx_0_xy[i] = -2.0 * g_y_xx_z_xy[i] * c_exps[i] + 4.0 * g_yzz_xx_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_xx_0_xz[i] = -2.0 * g_y_xx_z_xz[i] * c_exps[i] + 4.0 * g_yzz_xx_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_xx_0_yy[i] = -2.0 * g_y_xx_z_yy[i] * c_exps[i] + 4.0 * g_yzz_xx_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_xx_0_yz[i] = -2.0 * g_y_xx_z_yz[i] * c_exps[i] + 4.0 * g_yzz_xx_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_xx_0_zz[i] = -2.0 * g_y_xx_z_zz[i] * c_exps[i] + 4.0 * g_yzz_xx_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1878-1884)

    #pragma omp simd aligned(g_y_xy_z_xx, g_y_xy_z_xy, g_y_xy_z_xz, g_y_xy_z_yy, g_y_xy_z_yz, g_y_xy_z_zz, g_yzz_xy_z_xx, g_yzz_xy_z_xy, g_yzz_xy_z_xz, g_yzz_xy_z_yy, g_yzz_xy_z_yz, g_yzz_xy_z_zz, g_z_0_z_0_yz_xy_0_xx, g_z_0_z_0_yz_xy_0_xy, g_z_0_z_0_yz_xy_0_xz, g_z_0_z_0_yz_xy_0_yy, g_z_0_z_0_yz_xy_0_yz, g_z_0_z_0_yz_xy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_yz_xy_0_xx[i] = -2.0 * g_y_xy_z_xx[i] * c_exps[i] + 4.0 * g_yzz_xy_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_xy_0_xy[i] = -2.0 * g_y_xy_z_xy[i] * c_exps[i] + 4.0 * g_yzz_xy_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_xy_0_xz[i] = -2.0 * g_y_xy_z_xz[i] * c_exps[i] + 4.0 * g_yzz_xy_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_xy_0_yy[i] = -2.0 * g_y_xy_z_yy[i] * c_exps[i] + 4.0 * g_yzz_xy_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_xy_0_yz[i] = -2.0 * g_y_xy_z_yz[i] * c_exps[i] + 4.0 * g_yzz_xy_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_xy_0_zz[i] = -2.0 * g_y_xy_z_zz[i] * c_exps[i] + 4.0 * g_yzz_xy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1884-1890)

    #pragma omp simd aligned(g_y_xz_z_xx, g_y_xz_z_xy, g_y_xz_z_xz, g_y_xz_z_yy, g_y_xz_z_yz, g_y_xz_z_zz, g_yzz_xz_z_xx, g_yzz_xz_z_xy, g_yzz_xz_z_xz, g_yzz_xz_z_yy, g_yzz_xz_z_yz, g_yzz_xz_z_zz, g_z_0_z_0_yz_xz_0_xx, g_z_0_z_0_yz_xz_0_xy, g_z_0_z_0_yz_xz_0_xz, g_z_0_z_0_yz_xz_0_yy, g_z_0_z_0_yz_xz_0_yz, g_z_0_z_0_yz_xz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_yz_xz_0_xx[i] = -2.0 * g_y_xz_z_xx[i] * c_exps[i] + 4.0 * g_yzz_xz_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_xz_0_xy[i] = -2.0 * g_y_xz_z_xy[i] * c_exps[i] + 4.0 * g_yzz_xz_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_xz_0_xz[i] = -2.0 * g_y_xz_z_xz[i] * c_exps[i] + 4.0 * g_yzz_xz_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_xz_0_yy[i] = -2.0 * g_y_xz_z_yy[i] * c_exps[i] + 4.0 * g_yzz_xz_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_xz_0_yz[i] = -2.0 * g_y_xz_z_yz[i] * c_exps[i] + 4.0 * g_yzz_xz_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_xz_0_zz[i] = -2.0 * g_y_xz_z_zz[i] * c_exps[i] + 4.0 * g_yzz_xz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1890-1896)

    #pragma omp simd aligned(g_y_yy_z_xx, g_y_yy_z_xy, g_y_yy_z_xz, g_y_yy_z_yy, g_y_yy_z_yz, g_y_yy_z_zz, g_yzz_yy_z_xx, g_yzz_yy_z_xy, g_yzz_yy_z_xz, g_yzz_yy_z_yy, g_yzz_yy_z_yz, g_yzz_yy_z_zz, g_z_0_z_0_yz_yy_0_xx, g_z_0_z_0_yz_yy_0_xy, g_z_0_z_0_yz_yy_0_xz, g_z_0_z_0_yz_yy_0_yy, g_z_0_z_0_yz_yy_0_yz, g_z_0_z_0_yz_yy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_yz_yy_0_xx[i] = -2.0 * g_y_yy_z_xx[i] * c_exps[i] + 4.0 * g_yzz_yy_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_yy_0_xy[i] = -2.0 * g_y_yy_z_xy[i] * c_exps[i] + 4.0 * g_yzz_yy_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_yy_0_xz[i] = -2.0 * g_y_yy_z_xz[i] * c_exps[i] + 4.0 * g_yzz_yy_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_yy_0_yy[i] = -2.0 * g_y_yy_z_yy[i] * c_exps[i] + 4.0 * g_yzz_yy_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_yy_0_yz[i] = -2.0 * g_y_yy_z_yz[i] * c_exps[i] + 4.0 * g_yzz_yy_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_yy_0_zz[i] = -2.0 * g_y_yy_z_zz[i] * c_exps[i] + 4.0 * g_yzz_yy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1896-1902)

    #pragma omp simd aligned(g_y_yz_z_xx, g_y_yz_z_xy, g_y_yz_z_xz, g_y_yz_z_yy, g_y_yz_z_yz, g_y_yz_z_zz, g_yzz_yz_z_xx, g_yzz_yz_z_xy, g_yzz_yz_z_xz, g_yzz_yz_z_yy, g_yzz_yz_z_yz, g_yzz_yz_z_zz, g_z_0_z_0_yz_yz_0_xx, g_z_0_z_0_yz_yz_0_xy, g_z_0_z_0_yz_yz_0_xz, g_z_0_z_0_yz_yz_0_yy, g_z_0_z_0_yz_yz_0_yz, g_z_0_z_0_yz_yz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_yz_yz_0_xx[i] = -2.0 * g_y_yz_z_xx[i] * c_exps[i] + 4.0 * g_yzz_yz_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_yz_0_xy[i] = -2.0 * g_y_yz_z_xy[i] * c_exps[i] + 4.0 * g_yzz_yz_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_yz_0_xz[i] = -2.0 * g_y_yz_z_xz[i] * c_exps[i] + 4.0 * g_yzz_yz_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_yz_0_yy[i] = -2.0 * g_y_yz_z_yy[i] * c_exps[i] + 4.0 * g_yzz_yz_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_yz_0_yz[i] = -2.0 * g_y_yz_z_yz[i] * c_exps[i] + 4.0 * g_yzz_yz_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_yz_0_zz[i] = -2.0 * g_y_yz_z_zz[i] * c_exps[i] + 4.0 * g_yzz_yz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1902-1908)

    #pragma omp simd aligned(g_y_zz_z_xx, g_y_zz_z_xy, g_y_zz_z_xz, g_y_zz_z_yy, g_y_zz_z_yz, g_y_zz_z_zz, g_yzz_zz_z_xx, g_yzz_zz_z_xy, g_yzz_zz_z_xz, g_yzz_zz_z_yy, g_yzz_zz_z_yz, g_yzz_zz_z_zz, g_z_0_z_0_yz_zz_0_xx, g_z_0_z_0_yz_zz_0_xy, g_z_0_z_0_yz_zz_0_xz, g_z_0_z_0_yz_zz_0_yy, g_z_0_z_0_yz_zz_0_yz, g_z_0_z_0_yz_zz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_yz_zz_0_xx[i] = -2.0 * g_y_zz_z_xx[i] * c_exps[i] + 4.0 * g_yzz_zz_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_zz_0_xy[i] = -2.0 * g_y_zz_z_xy[i] * c_exps[i] + 4.0 * g_yzz_zz_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_zz_0_xz[i] = -2.0 * g_y_zz_z_xz[i] * c_exps[i] + 4.0 * g_yzz_zz_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_zz_0_yy[i] = -2.0 * g_y_zz_z_yy[i] * c_exps[i] + 4.0 * g_yzz_zz_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_zz_0_yz[i] = -2.0 * g_y_zz_z_yz[i] * c_exps[i] + 4.0 * g_yzz_zz_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_zz_0_zz[i] = -2.0 * g_y_zz_z_zz[i] * c_exps[i] + 4.0 * g_yzz_zz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1908-1914)

    #pragma omp simd aligned(g_z_0_z_0_zz_xx_0_xx, g_z_0_z_0_zz_xx_0_xy, g_z_0_z_0_zz_xx_0_xz, g_z_0_z_0_zz_xx_0_yy, g_z_0_z_0_zz_xx_0_yz, g_z_0_z_0_zz_xx_0_zz, g_z_xx_z_xx, g_z_xx_z_xy, g_z_xx_z_xz, g_z_xx_z_yy, g_z_xx_z_yz, g_z_xx_z_zz, g_zzz_xx_z_xx, g_zzz_xx_z_xy, g_zzz_xx_z_xz, g_zzz_xx_z_yy, g_zzz_xx_z_yz, g_zzz_xx_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_zz_xx_0_xx[i] = -4.0 * g_z_xx_z_xx[i] * c_exps[i] + 4.0 * g_zzz_xx_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_xx_0_xy[i] = -4.0 * g_z_xx_z_xy[i] * c_exps[i] + 4.0 * g_zzz_xx_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_xx_0_xz[i] = -4.0 * g_z_xx_z_xz[i] * c_exps[i] + 4.0 * g_zzz_xx_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_xx_0_yy[i] = -4.0 * g_z_xx_z_yy[i] * c_exps[i] + 4.0 * g_zzz_xx_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_xx_0_yz[i] = -4.0 * g_z_xx_z_yz[i] * c_exps[i] + 4.0 * g_zzz_xx_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_xx_0_zz[i] = -4.0 * g_z_xx_z_zz[i] * c_exps[i] + 4.0 * g_zzz_xx_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1914-1920)

    #pragma omp simd aligned(g_z_0_z_0_zz_xy_0_xx, g_z_0_z_0_zz_xy_0_xy, g_z_0_z_0_zz_xy_0_xz, g_z_0_z_0_zz_xy_0_yy, g_z_0_z_0_zz_xy_0_yz, g_z_0_z_0_zz_xy_0_zz, g_z_xy_z_xx, g_z_xy_z_xy, g_z_xy_z_xz, g_z_xy_z_yy, g_z_xy_z_yz, g_z_xy_z_zz, g_zzz_xy_z_xx, g_zzz_xy_z_xy, g_zzz_xy_z_xz, g_zzz_xy_z_yy, g_zzz_xy_z_yz, g_zzz_xy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_zz_xy_0_xx[i] = -4.0 * g_z_xy_z_xx[i] * c_exps[i] + 4.0 * g_zzz_xy_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_xy_0_xy[i] = -4.0 * g_z_xy_z_xy[i] * c_exps[i] + 4.0 * g_zzz_xy_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_xy_0_xz[i] = -4.0 * g_z_xy_z_xz[i] * c_exps[i] + 4.0 * g_zzz_xy_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_xy_0_yy[i] = -4.0 * g_z_xy_z_yy[i] * c_exps[i] + 4.0 * g_zzz_xy_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_xy_0_yz[i] = -4.0 * g_z_xy_z_yz[i] * c_exps[i] + 4.0 * g_zzz_xy_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_xy_0_zz[i] = -4.0 * g_z_xy_z_zz[i] * c_exps[i] + 4.0 * g_zzz_xy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1920-1926)

    #pragma omp simd aligned(g_z_0_z_0_zz_xz_0_xx, g_z_0_z_0_zz_xz_0_xy, g_z_0_z_0_zz_xz_0_xz, g_z_0_z_0_zz_xz_0_yy, g_z_0_z_0_zz_xz_0_yz, g_z_0_z_0_zz_xz_0_zz, g_z_xz_z_xx, g_z_xz_z_xy, g_z_xz_z_xz, g_z_xz_z_yy, g_z_xz_z_yz, g_z_xz_z_zz, g_zzz_xz_z_xx, g_zzz_xz_z_xy, g_zzz_xz_z_xz, g_zzz_xz_z_yy, g_zzz_xz_z_yz, g_zzz_xz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_zz_xz_0_xx[i] = -4.0 * g_z_xz_z_xx[i] * c_exps[i] + 4.0 * g_zzz_xz_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_xz_0_xy[i] = -4.0 * g_z_xz_z_xy[i] * c_exps[i] + 4.0 * g_zzz_xz_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_xz_0_xz[i] = -4.0 * g_z_xz_z_xz[i] * c_exps[i] + 4.0 * g_zzz_xz_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_xz_0_yy[i] = -4.0 * g_z_xz_z_yy[i] * c_exps[i] + 4.0 * g_zzz_xz_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_xz_0_yz[i] = -4.0 * g_z_xz_z_yz[i] * c_exps[i] + 4.0 * g_zzz_xz_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_xz_0_zz[i] = -4.0 * g_z_xz_z_zz[i] * c_exps[i] + 4.0 * g_zzz_xz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1926-1932)

    #pragma omp simd aligned(g_z_0_z_0_zz_yy_0_xx, g_z_0_z_0_zz_yy_0_xy, g_z_0_z_0_zz_yy_0_xz, g_z_0_z_0_zz_yy_0_yy, g_z_0_z_0_zz_yy_0_yz, g_z_0_z_0_zz_yy_0_zz, g_z_yy_z_xx, g_z_yy_z_xy, g_z_yy_z_xz, g_z_yy_z_yy, g_z_yy_z_yz, g_z_yy_z_zz, g_zzz_yy_z_xx, g_zzz_yy_z_xy, g_zzz_yy_z_xz, g_zzz_yy_z_yy, g_zzz_yy_z_yz, g_zzz_yy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_zz_yy_0_xx[i] = -4.0 * g_z_yy_z_xx[i] * c_exps[i] + 4.0 * g_zzz_yy_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_yy_0_xy[i] = -4.0 * g_z_yy_z_xy[i] * c_exps[i] + 4.0 * g_zzz_yy_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_yy_0_xz[i] = -4.0 * g_z_yy_z_xz[i] * c_exps[i] + 4.0 * g_zzz_yy_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_yy_0_yy[i] = -4.0 * g_z_yy_z_yy[i] * c_exps[i] + 4.0 * g_zzz_yy_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_yy_0_yz[i] = -4.0 * g_z_yy_z_yz[i] * c_exps[i] + 4.0 * g_zzz_yy_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_yy_0_zz[i] = -4.0 * g_z_yy_z_zz[i] * c_exps[i] + 4.0 * g_zzz_yy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1932-1938)

    #pragma omp simd aligned(g_z_0_z_0_zz_yz_0_xx, g_z_0_z_0_zz_yz_0_xy, g_z_0_z_0_zz_yz_0_xz, g_z_0_z_0_zz_yz_0_yy, g_z_0_z_0_zz_yz_0_yz, g_z_0_z_0_zz_yz_0_zz, g_z_yz_z_xx, g_z_yz_z_xy, g_z_yz_z_xz, g_z_yz_z_yy, g_z_yz_z_yz, g_z_yz_z_zz, g_zzz_yz_z_xx, g_zzz_yz_z_xy, g_zzz_yz_z_xz, g_zzz_yz_z_yy, g_zzz_yz_z_yz, g_zzz_yz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_zz_yz_0_xx[i] = -4.0 * g_z_yz_z_xx[i] * c_exps[i] + 4.0 * g_zzz_yz_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_yz_0_xy[i] = -4.0 * g_z_yz_z_xy[i] * c_exps[i] + 4.0 * g_zzz_yz_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_yz_0_xz[i] = -4.0 * g_z_yz_z_xz[i] * c_exps[i] + 4.0 * g_zzz_yz_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_yz_0_yy[i] = -4.0 * g_z_yz_z_yy[i] * c_exps[i] + 4.0 * g_zzz_yz_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_yz_0_yz[i] = -4.0 * g_z_yz_z_yz[i] * c_exps[i] + 4.0 * g_zzz_yz_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_yz_0_zz[i] = -4.0 * g_z_yz_z_zz[i] * c_exps[i] + 4.0 * g_zzz_yz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1938-1944)

    #pragma omp simd aligned(g_z_0_z_0_zz_zz_0_xx, g_z_0_z_0_zz_zz_0_xy, g_z_0_z_0_zz_zz_0_xz, g_z_0_z_0_zz_zz_0_yy, g_z_0_z_0_zz_zz_0_yz, g_z_0_z_0_zz_zz_0_zz, g_z_zz_z_xx, g_z_zz_z_xy, g_z_zz_z_xz, g_z_zz_z_yy, g_z_zz_z_yz, g_z_zz_z_zz, g_zzz_zz_z_xx, g_zzz_zz_z_xy, g_zzz_zz_z_xz, g_zzz_zz_z_yy, g_zzz_zz_z_yz, g_zzz_zz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_zz_zz_0_xx[i] = -4.0 * g_z_zz_z_xx[i] * c_exps[i] + 4.0 * g_zzz_zz_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_zz_0_xy[i] = -4.0 * g_z_zz_z_xy[i] * c_exps[i] + 4.0 * g_zzz_zz_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_zz_0_xz[i] = -4.0 * g_z_zz_z_xz[i] * c_exps[i] + 4.0 * g_zzz_zz_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_zz_0_yy[i] = -4.0 * g_z_zz_z_yy[i] * c_exps[i] + 4.0 * g_zzz_zz_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_zz_0_yz[i] = -4.0 * g_z_zz_z_yz[i] * c_exps[i] + 4.0 * g_zzz_zz_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_zz_0_zz[i] = -4.0 * g_z_zz_z_zz[i] * c_exps[i] + 4.0 * g_zzz_zz_z_zz[i] * a_exp * c_exps[i];
    }
}

} // t4c_geom namespace


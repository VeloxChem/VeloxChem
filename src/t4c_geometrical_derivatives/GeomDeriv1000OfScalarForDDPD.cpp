#include "GeomDeriv1000OfScalarForDDPD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1000_ddpd_0(CSimdArray<double>& buffer_1000_ddpd,
                     const CSimdArray<double>& buffer_pdpd,
                     const CSimdArray<double>& buffer_fdpd,
                     const double a_exp) -> void
{
    const auto ndims = buffer_1000_ddpd.number_of_columns();

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

    /// Set up components of integrals buffer : buffer_1000_ddpd

    auto g_x_0_0_0_xx_xx_x_xx = buffer_1000_ddpd[0];

    auto g_x_0_0_0_xx_xx_x_xy = buffer_1000_ddpd[1];

    auto g_x_0_0_0_xx_xx_x_xz = buffer_1000_ddpd[2];

    auto g_x_0_0_0_xx_xx_x_yy = buffer_1000_ddpd[3];

    auto g_x_0_0_0_xx_xx_x_yz = buffer_1000_ddpd[4];

    auto g_x_0_0_0_xx_xx_x_zz = buffer_1000_ddpd[5];

    auto g_x_0_0_0_xx_xx_y_xx = buffer_1000_ddpd[6];

    auto g_x_0_0_0_xx_xx_y_xy = buffer_1000_ddpd[7];

    auto g_x_0_0_0_xx_xx_y_xz = buffer_1000_ddpd[8];

    auto g_x_0_0_0_xx_xx_y_yy = buffer_1000_ddpd[9];

    auto g_x_0_0_0_xx_xx_y_yz = buffer_1000_ddpd[10];

    auto g_x_0_0_0_xx_xx_y_zz = buffer_1000_ddpd[11];

    auto g_x_0_0_0_xx_xx_z_xx = buffer_1000_ddpd[12];

    auto g_x_0_0_0_xx_xx_z_xy = buffer_1000_ddpd[13];

    auto g_x_0_0_0_xx_xx_z_xz = buffer_1000_ddpd[14];

    auto g_x_0_0_0_xx_xx_z_yy = buffer_1000_ddpd[15];

    auto g_x_0_0_0_xx_xx_z_yz = buffer_1000_ddpd[16];

    auto g_x_0_0_0_xx_xx_z_zz = buffer_1000_ddpd[17];

    auto g_x_0_0_0_xx_xy_x_xx = buffer_1000_ddpd[18];

    auto g_x_0_0_0_xx_xy_x_xy = buffer_1000_ddpd[19];

    auto g_x_0_0_0_xx_xy_x_xz = buffer_1000_ddpd[20];

    auto g_x_0_0_0_xx_xy_x_yy = buffer_1000_ddpd[21];

    auto g_x_0_0_0_xx_xy_x_yz = buffer_1000_ddpd[22];

    auto g_x_0_0_0_xx_xy_x_zz = buffer_1000_ddpd[23];

    auto g_x_0_0_0_xx_xy_y_xx = buffer_1000_ddpd[24];

    auto g_x_0_0_0_xx_xy_y_xy = buffer_1000_ddpd[25];

    auto g_x_0_0_0_xx_xy_y_xz = buffer_1000_ddpd[26];

    auto g_x_0_0_0_xx_xy_y_yy = buffer_1000_ddpd[27];

    auto g_x_0_0_0_xx_xy_y_yz = buffer_1000_ddpd[28];

    auto g_x_0_0_0_xx_xy_y_zz = buffer_1000_ddpd[29];

    auto g_x_0_0_0_xx_xy_z_xx = buffer_1000_ddpd[30];

    auto g_x_0_0_0_xx_xy_z_xy = buffer_1000_ddpd[31];

    auto g_x_0_0_0_xx_xy_z_xz = buffer_1000_ddpd[32];

    auto g_x_0_0_0_xx_xy_z_yy = buffer_1000_ddpd[33];

    auto g_x_0_0_0_xx_xy_z_yz = buffer_1000_ddpd[34];

    auto g_x_0_0_0_xx_xy_z_zz = buffer_1000_ddpd[35];

    auto g_x_0_0_0_xx_xz_x_xx = buffer_1000_ddpd[36];

    auto g_x_0_0_0_xx_xz_x_xy = buffer_1000_ddpd[37];

    auto g_x_0_0_0_xx_xz_x_xz = buffer_1000_ddpd[38];

    auto g_x_0_0_0_xx_xz_x_yy = buffer_1000_ddpd[39];

    auto g_x_0_0_0_xx_xz_x_yz = buffer_1000_ddpd[40];

    auto g_x_0_0_0_xx_xz_x_zz = buffer_1000_ddpd[41];

    auto g_x_0_0_0_xx_xz_y_xx = buffer_1000_ddpd[42];

    auto g_x_0_0_0_xx_xz_y_xy = buffer_1000_ddpd[43];

    auto g_x_0_0_0_xx_xz_y_xz = buffer_1000_ddpd[44];

    auto g_x_0_0_0_xx_xz_y_yy = buffer_1000_ddpd[45];

    auto g_x_0_0_0_xx_xz_y_yz = buffer_1000_ddpd[46];

    auto g_x_0_0_0_xx_xz_y_zz = buffer_1000_ddpd[47];

    auto g_x_0_0_0_xx_xz_z_xx = buffer_1000_ddpd[48];

    auto g_x_0_0_0_xx_xz_z_xy = buffer_1000_ddpd[49];

    auto g_x_0_0_0_xx_xz_z_xz = buffer_1000_ddpd[50];

    auto g_x_0_0_0_xx_xz_z_yy = buffer_1000_ddpd[51];

    auto g_x_0_0_0_xx_xz_z_yz = buffer_1000_ddpd[52];

    auto g_x_0_0_0_xx_xz_z_zz = buffer_1000_ddpd[53];

    auto g_x_0_0_0_xx_yy_x_xx = buffer_1000_ddpd[54];

    auto g_x_0_0_0_xx_yy_x_xy = buffer_1000_ddpd[55];

    auto g_x_0_0_0_xx_yy_x_xz = buffer_1000_ddpd[56];

    auto g_x_0_0_0_xx_yy_x_yy = buffer_1000_ddpd[57];

    auto g_x_0_0_0_xx_yy_x_yz = buffer_1000_ddpd[58];

    auto g_x_0_0_0_xx_yy_x_zz = buffer_1000_ddpd[59];

    auto g_x_0_0_0_xx_yy_y_xx = buffer_1000_ddpd[60];

    auto g_x_0_0_0_xx_yy_y_xy = buffer_1000_ddpd[61];

    auto g_x_0_0_0_xx_yy_y_xz = buffer_1000_ddpd[62];

    auto g_x_0_0_0_xx_yy_y_yy = buffer_1000_ddpd[63];

    auto g_x_0_0_0_xx_yy_y_yz = buffer_1000_ddpd[64];

    auto g_x_0_0_0_xx_yy_y_zz = buffer_1000_ddpd[65];

    auto g_x_0_0_0_xx_yy_z_xx = buffer_1000_ddpd[66];

    auto g_x_0_0_0_xx_yy_z_xy = buffer_1000_ddpd[67];

    auto g_x_0_0_0_xx_yy_z_xz = buffer_1000_ddpd[68];

    auto g_x_0_0_0_xx_yy_z_yy = buffer_1000_ddpd[69];

    auto g_x_0_0_0_xx_yy_z_yz = buffer_1000_ddpd[70];

    auto g_x_0_0_0_xx_yy_z_zz = buffer_1000_ddpd[71];

    auto g_x_0_0_0_xx_yz_x_xx = buffer_1000_ddpd[72];

    auto g_x_0_0_0_xx_yz_x_xy = buffer_1000_ddpd[73];

    auto g_x_0_0_0_xx_yz_x_xz = buffer_1000_ddpd[74];

    auto g_x_0_0_0_xx_yz_x_yy = buffer_1000_ddpd[75];

    auto g_x_0_0_0_xx_yz_x_yz = buffer_1000_ddpd[76];

    auto g_x_0_0_0_xx_yz_x_zz = buffer_1000_ddpd[77];

    auto g_x_0_0_0_xx_yz_y_xx = buffer_1000_ddpd[78];

    auto g_x_0_0_0_xx_yz_y_xy = buffer_1000_ddpd[79];

    auto g_x_0_0_0_xx_yz_y_xz = buffer_1000_ddpd[80];

    auto g_x_0_0_0_xx_yz_y_yy = buffer_1000_ddpd[81];

    auto g_x_0_0_0_xx_yz_y_yz = buffer_1000_ddpd[82];

    auto g_x_0_0_0_xx_yz_y_zz = buffer_1000_ddpd[83];

    auto g_x_0_0_0_xx_yz_z_xx = buffer_1000_ddpd[84];

    auto g_x_0_0_0_xx_yz_z_xy = buffer_1000_ddpd[85];

    auto g_x_0_0_0_xx_yz_z_xz = buffer_1000_ddpd[86];

    auto g_x_0_0_0_xx_yz_z_yy = buffer_1000_ddpd[87];

    auto g_x_0_0_0_xx_yz_z_yz = buffer_1000_ddpd[88];

    auto g_x_0_0_0_xx_yz_z_zz = buffer_1000_ddpd[89];

    auto g_x_0_0_0_xx_zz_x_xx = buffer_1000_ddpd[90];

    auto g_x_0_0_0_xx_zz_x_xy = buffer_1000_ddpd[91];

    auto g_x_0_0_0_xx_zz_x_xz = buffer_1000_ddpd[92];

    auto g_x_0_0_0_xx_zz_x_yy = buffer_1000_ddpd[93];

    auto g_x_0_0_0_xx_zz_x_yz = buffer_1000_ddpd[94];

    auto g_x_0_0_0_xx_zz_x_zz = buffer_1000_ddpd[95];

    auto g_x_0_0_0_xx_zz_y_xx = buffer_1000_ddpd[96];

    auto g_x_0_0_0_xx_zz_y_xy = buffer_1000_ddpd[97];

    auto g_x_0_0_0_xx_zz_y_xz = buffer_1000_ddpd[98];

    auto g_x_0_0_0_xx_zz_y_yy = buffer_1000_ddpd[99];

    auto g_x_0_0_0_xx_zz_y_yz = buffer_1000_ddpd[100];

    auto g_x_0_0_0_xx_zz_y_zz = buffer_1000_ddpd[101];

    auto g_x_0_0_0_xx_zz_z_xx = buffer_1000_ddpd[102];

    auto g_x_0_0_0_xx_zz_z_xy = buffer_1000_ddpd[103];

    auto g_x_0_0_0_xx_zz_z_xz = buffer_1000_ddpd[104];

    auto g_x_0_0_0_xx_zz_z_yy = buffer_1000_ddpd[105];

    auto g_x_0_0_0_xx_zz_z_yz = buffer_1000_ddpd[106];

    auto g_x_0_0_0_xx_zz_z_zz = buffer_1000_ddpd[107];

    auto g_x_0_0_0_xy_xx_x_xx = buffer_1000_ddpd[108];

    auto g_x_0_0_0_xy_xx_x_xy = buffer_1000_ddpd[109];

    auto g_x_0_0_0_xy_xx_x_xz = buffer_1000_ddpd[110];

    auto g_x_0_0_0_xy_xx_x_yy = buffer_1000_ddpd[111];

    auto g_x_0_0_0_xy_xx_x_yz = buffer_1000_ddpd[112];

    auto g_x_0_0_0_xy_xx_x_zz = buffer_1000_ddpd[113];

    auto g_x_0_0_0_xy_xx_y_xx = buffer_1000_ddpd[114];

    auto g_x_0_0_0_xy_xx_y_xy = buffer_1000_ddpd[115];

    auto g_x_0_0_0_xy_xx_y_xz = buffer_1000_ddpd[116];

    auto g_x_0_0_0_xy_xx_y_yy = buffer_1000_ddpd[117];

    auto g_x_0_0_0_xy_xx_y_yz = buffer_1000_ddpd[118];

    auto g_x_0_0_0_xy_xx_y_zz = buffer_1000_ddpd[119];

    auto g_x_0_0_0_xy_xx_z_xx = buffer_1000_ddpd[120];

    auto g_x_0_0_0_xy_xx_z_xy = buffer_1000_ddpd[121];

    auto g_x_0_0_0_xy_xx_z_xz = buffer_1000_ddpd[122];

    auto g_x_0_0_0_xy_xx_z_yy = buffer_1000_ddpd[123];

    auto g_x_0_0_0_xy_xx_z_yz = buffer_1000_ddpd[124];

    auto g_x_0_0_0_xy_xx_z_zz = buffer_1000_ddpd[125];

    auto g_x_0_0_0_xy_xy_x_xx = buffer_1000_ddpd[126];

    auto g_x_0_0_0_xy_xy_x_xy = buffer_1000_ddpd[127];

    auto g_x_0_0_0_xy_xy_x_xz = buffer_1000_ddpd[128];

    auto g_x_0_0_0_xy_xy_x_yy = buffer_1000_ddpd[129];

    auto g_x_0_0_0_xy_xy_x_yz = buffer_1000_ddpd[130];

    auto g_x_0_0_0_xy_xy_x_zz = buffer_1000_ddpd[131];

    auto g_x_0_0_0_xy_xy_y_xx = buffer_1000_ddpd[132];

    auto g_x_0_0_0_xy_xy_y_xy = buffer_1000_ddpd[133];

    auto g_x_0_0_0_xy_xy_y_xz = buffer_1000_ddpd[134];

    auto g_x_0_0_0_xy_xy_y_yy = buffer_1000_ddpd[135];

    auto g_x_0_0_0_xy_xy_y_yz = buffer_1000_ddpd[136];

    auto g_x_0_0_0_xy_xy_y_zz = buffer_1000_ddpd[137];

    auto g_x_0_0_0_xy_xy_z_xx = buffer_1000_ddpd[138];

    auto g_x_0_0_0_xy_xy_z_xy = buffer_1000_ddpd[139];

    auto g_x_0_0_0_xy_xy_z_xz = buffer_1000_ddpd[140];

    auto g_x_0_0_0_xy_xy_z_yy = buffer_1000_ddpd[141];

    auto g_x_0_0_0_xy_xy_z_yz = buffer_1000_ddpd[142];

    auto g_x_0_0_0_xy_xy_z_zz = buffer_1000_ddpd[143];

    auto g_x_0_0_0_xy_xz_x_xx = buffer_1000_ddpd[144];

    auto g_x_0_0_0_xy_xz_x_xy = buffer_1000_ddpd[145];

    auto g_x_0_0_0_xy_xz_x_xz = buffer_1000_ddpd[146];

    auto g_x_0_0_0_xy_xz_x_yy = buffer_1000_ddpd[147];

    auto g_x_0_0_0_xy_xz_x_yz = buffer_1000_ddpd[148];

    auto g_x_0_0_0_xy_xz_x_zz = buffer_1000_ddpd[149];

    auto g_x_0_0_0_xy_xz_y_xx = buffer_1000_ddpd[150];

    auto g_x_0_0_0_xy_xz_y_xy = buffer_1000_ddpd[151];

    auto g_x_0_0_0_xy_xz_y_xz = buffer_1000_ddpd[152];

    auto g_x_0_0_0_xy_xz_y_yy = buffer_1000_ddpd[153];

    auto g_x_0_0_0_xy_xz_y_yz = buffer_1000_ddpd[154];

    auto g_x_0_0_0_xy_xz_y_zz = buffer_1000_ddpd[155];

    auto g_x_0_0_0_xy_xz_z_xx = buffer_1000_ddpd[156];

    auto g_x_0_0_0_xy_xz_z_xy = buffer_1000_ddpd[157];

    auto g_x_0_0_0_xy_xz_z_xz = buffer_1000_ddpd[158];

    auto g_x_0_0_0_xy_xz_z_yy = buffer_1000_ddpd[159];

    auto g_x_0_0_0_xy_xz_z_yz = buffer_1000_ddpd[160];

    auto g_x_0_0_0_xy_xz_z_zz = buffer_1000_ddpd[161];

    auto g_x_0_0_0_xy_yy_x_xx = buffer_1000_ddpd[162];

    auto g_x_0_0_0_xy_yy_x_xy = buffer_1000_ddpd[163];

    auto g_x_0_0_0_xy_yy_x_xz = buffer_1000_ddpd[164];

    auto g_x_0_0_0_xy_yy_x_yy = buffer_1000_ddpd[165];

    auto g_x_0_0_0_xy_yy_x_yz = buffer_1000_ddpd[166];

    auto g_x_0_0_0_xy_yy_x_zz = buffer_1000_ddpd[167];

    auto g_x_0_0_0_xy_yy_y_xx = buffer_1000_ddpd[168];

    auto g_x_0_0_0_xy_yy_y_xy = buffer_1000_ddpd[169];

    auto g_x_0_0_0_xy_yy_y_xz = buffer_1000_ddpd[170];

    auto g_x_0_0_0_xy_yy_y_yy = buffer_1000_ddpd[171];

    auto g_x_0_0_0_xy_yy_y_yz = buffer_1000_ddpd[172];

    auto g_x_0_0_0_xy_yy_y_zz = buffer_1000_ddpd[173];

    auto g_x_0_0_0_xy_yy_z_xx = buffer_1000_ddpd[174];

    auto g_x_0_0_0_xy_yy_z_xy = buffer_1000_ddpd[175];

    auto g_x_0_0_0_xy_yy_z_xz = buffer_1000_ddpd[176];

    auto g_x_0_0_0_xy_yy_z_yy = buffer_1000_ddpd[177];

    auto g_x_0_0_0_xy_yy_z_yz = buffer_1000_ddpd[178];

    auto g_x_0_0_0_xy_yy_z_zz = buffer_1000_ddpd[179];

    auto g_x_0_0_0_xy_yz_x_xx = buffer_1000_ddpd[180];

    auto g_x_0_0_0_xy_yz_x_xy = buffer_1000_ddpd[181];

    auto g_x_0_0_0_xy_yz_x_xz = buffer_1000_ddpd[182];

    auto g_x_0_0_0_xy_yz_x_yy = buffer_1000_ddpd[183];

    auto g_x_0_0_0_xy_yz_x_yz = buffer_1000_ddpd[184];

    auto g_x_0_0_0_xy_yz_x_zz = buffer_1000_ddpd[185];

    auto g_x_0_0_0_xy_yz_y_xx = buffer_1000_ddpd[186];

    auto g_x_0_0_0_xy_yz_y_xy = buffer_1000_ddpd[187];

    auto g_x_0_0_0_xy_yz_y_xz = buffer_1000_ddpd[188];

    auto g_x_0_0_0_xy_yz_y_yy = buffer_1000_ddpd[189];

    auto g_x_0_0_0_xy_yz_y_yz = buffer_1000_ddpd[190];

    auto g_x_0_0_0_xy_yz_y_zz = buffer_1000_ddpd[191];

    auto g_x_0_0_0_xy_yz_z_xx = buffer_1000_ddpd[192];

    auto g_x_0_0_0_xy_yz_z_xy = buffer_1000_ddpd[193];

    auto g_x_0_0_0_xy_yz_z_xz = buffer_1000_ddpd[194];

    auto g_x_0_0_0_xy_yz_z_yy = buffer_1000_ddpd[195];

    auto g_x_0_0_0_xy_yz_z_yz = buffer_1000_ddpd[196];

    auto g_x_0_0_0_xy_yz_z_zz = buffer_1000_ddpd[197];

    auto g_x_0_0_0_xy_zz_x_xx = buffer_1000_ddpd[198];

    auto g_x_0_0_0_xy_zz_x_xy = buffer_1000_ddpd[199];

    auto g_x_0_0_0_xy_zz_x_xz = buffer_1000_ddpd[200];

    auto g_x_0_0_0_xy_zz_x_yy = buffer_1000_ddpd[201];

    auto g_x_0_0_0_xy_zz_x_yz = buffer_1000_ddpd[202];

    auto g_x_0_0_0_xy_zz_x_zz = buffer_1000_ddpd[203];

    auto g_x_0_0_0_xy_zz_y_xx = buffer_1000_ddpd[204];

    auto g_x_0_0_0_xy_zz_y_xy = buffer_1000_ddpd[205];

    auto g_x_0_0_0_xy_zz_y_xz = buffer_1000_ddpd[206];

    auto g_x_0_0_0_xy_zz_y_yy = buffer_1000_ddpd[207];

    auto g_x_0_0_0_xy_zz_y_yz = buffer_1000_ddpd[208];

    auto g_x_0_0_0_xy_zz_y_zz = buffer_1000_ddpd[209];

    auto g_x_0_0_0_xy_zz_z_xx = buffer_1000_ddpd[210];

    auto g_x_0_0_0_xy_zz_z_xy = buffer_1000_ddpd[211];

    auto g_x_0_0_0_xy_zz_z_xz = buffer_1000_ddpd[212];

    auto g_x_0_0_0_xy_zz_z_yy = buffer_1000_ddpd[213];

    auto g_x_0_0_0_xy_zz_z_yz = buffer_1000_ddpd[214];

    auto g_x_0_0_0_xy_zz_z_zz = buffer_1000_ddpd[215];

    auto g_x_0_0_0_xz_xx_x_xx = buffer_1000_ddpd[216];

    auto g_x_0_0_0_xz_xx_x_xy = buffer_1000_ddpd[217];

    auto g_x_0_0_0_xz_xx_x_xz = buffer_1000_ddpd[218];

    auto g_x_0_0_0_xz_xx_x_yy = buffer_1000_ddpd[219];

    auto g_x_0_0_0_xz_xx_x_yz = buffer_1000_ddpd[220];

    auto g_x_0_0_0_xz_xx_x_zz = buffer_1000_ddpd[221];

    auto g_x_0_0_0_xz_xx_y_xx = buffer_1000_ddpd[222];

    auto g_x_0_0_0_xz_xx_y_xy = buffer_1000_ddpd[223];

    auto g_x_0_0_0_xz_xx_y_xz = buffer_1000_ddpd[224];

    auto g_x_0_0_0_xz_xx_y_yy = buffer_1000_ddpd[225];

    auto g_x_0_0_0_xz_xx_y_yz = buffer_1000_ddpd[226];

    auto g_x_0_0_0_xz_xx_y_zz = buffer_1000_ddpd[227];

    auto g_x_0_0_0_xz_xx_z_xx = buffer_1000_ddpd[228];

    auto g_x_0_0_0_xz_xx_z_xy = buffer_1000_ddpd[229];

    auto g_x_0_0_0_xz_xx_z_xz = buffer_1000_ddpd[230];

    auto g_x_0_0_0_xz_xx_z_yy = buffer_1000_ddpd[231];

    auto g_x_0_0_0_xz_xx_z_yz = buffer_1000_ddpd[232];

    auto g_x_0_0_0_xz_xx_z_zz = buffer_1000_ddpd[233];

    auto g_x_0_0_0_xz_xy_x_xx = buffer_1000_ddpd[234];

    auto g_x_0_0_0_xz_xy_x_xy = buffer_1000_ddpd[235];

    auto g_x_0_0_0_xz_xy_x_xz = buffer_1000_ddpd[236];

    auto g_x_0_0_0_xz_xy_x_yy = buffer_1000_ddpd[237];

    auto g_x_0_0_0_xz_xy_x_yz = buffer_1000_ddpd[238];

    auto g_x_0_0_0_xz_xy_x_zz = buffer_1000_ddpd[239];

    auto g_x_0_0_0_xz_xy_y_xx = buffer_1000_ddpd[240];

    auto g_x_0_0_0_xz_xy_y_xy = buffer_1000_ddpd[241];

    auto g_x_0_0_0_xz_xy_y_xz = buffer_1000_ddpd[242];

    auto g_x_0_0_0_xz_xy_y_yy = buffer_1000_ddpd[243];

    auto g_x_0_0_0_xz_xy_y_yz = buffer_1000_ddpd[244];

    auto g_x_0_0_0_xz_xy_y_zz = buffer_1000_ddpd[245];

    auto g_x_0_0_0_xz_xy_z_xx = buffer_1000_ddpd[246];

    auto g_x_0_0_0_xz_xy_z_xy = buffer_1000_ddpd[247];

    auto g_x_0_0_0_xz_xy_z_xz = buffer_1000_ddpd[248];

    auto g_x_0_0_0_xz_xy_z_yy = buffer_1000_ddpd[249];

    auto g_x_0_0_0_xz_xy_z_yz = buffer_1000_ddpd[250];

    auto g_x_0_0_0_xz_xy_z_zz = buffer_1000_ddpd[251];

    auto g_x_0_0_0_xz_xz_x_xx = buffer_1000_ddpd[252];

    auto g_x_0_0_0_xz_xz_x_xy = buffer_1000_ddpd[253];

    auto g_x_0_0_0_xz_xz_x_xz = buffer_1000_ddpd[254];

    auto g_x_0_0_0_xz_xz_x_yy = buffer_1000_ddpd[255];

    auto g_x_0_0_0_xz_xz_x_yz = buffer_1000_ddpd[256];

    auto g_x_0_0_0_xz_xz_x_zz = buffer_1000_ddpd[257];

    auto g_x_0_0_0_xz_xz_y_xx = buffer_1000_ddpd[258];

    auto g_x_0_0_0_xz_xz_y_xy = buffer_1000_ddpd[259];

    auto g_x_0_0_0_xz_xz_y_xz = buffer_1000_ddpd[260];

    auto g_x_0_0_0_xz_xz_y_yy = buffer_1000_ddpd[261];

    auto g_x_0_0_0_xz_xz_y_yz = buffer_1000_ddpd[262];

    auto g_x_0_0_0_xz_xz_y_zz = buffer_1000_ddpd[263];

    auto g_x_0_0_0_xz_xz_z_xx = buffer_1000_ddpd[264];

    auto g_x_0_0_0_xz_xz_z_xy = buffer_1000_ddpd[265];

    auto g_x_0_0_0_xz_xz_z_xz = buffer_1000_ddpd[266];

    auto g_x_0_0_0_xz_xz_z_yy = buffer_1000_ddpd[267];

    auto g_x_0_0_0_xz_xz_z_yz = buffer_1000_ddpd[268];

    auto g_x_0_0_0_xz_xz_z_zz = buffer_1000_ddpd[269];

    auto g_x_0_0_0_xz_yy_x_xx = buffer_1000_ddpd[270];

    auto g_x_0_0_0_xz_yy_x_xy = buffer_1000_ddpd[271];

    auto g_x_0_0_0_xz_yy_x_xz = buffer_1000_ddpd[272];

    auto g_x_0_0_0_xz_yy_x_yy = buffer_1000_ddpd[273];

    auto g_x_0_0_0_xz_yy_x_yz = buffer_1000_ddpd[274];

    auto g_x_0_0_0_xz_yy_x_zz = buffer_1000_ddpd[275];

    auto g_x_0_0_0_xz_yy_y_xx = buffer_1000_ddpd[276];

    auto g_x_0_0_0_xz_yy_y_xy = buffer_1000_ddpd[277];

    auto g_x_0_0_0_xz_yy_y_xz = buffer_1000_ddpd[278];

    auto g_x_0_0_0_xz_yy_y_yy = buffer_1000_ddpd[279];

    auto g_x_0_0_0_xz_yy_y_yz = buffer_1000_ddpd[280];

    auto g_x_0_0_0_xz_yy_y_zz = buffer_1000_ddpd[281];

    auto g_x_0_0_0_xz_yy_z_xx = buffer_1000_ddpd[282];

    auto g_x_0_0_0_xz_yy_z_xy = buffer_1000_ddpd[283];

    auto g_x_0_0_0_xz_yy_z_xz = buffer_1000_ddpd[284];

    auto g_x_0_0_0_xz_yy_z_yy = buffer_1000_ddpd[285];

    auto g_x_0_0_0_xz_yy_z_yz = buffer_1000_ddpd[286];

    auto g_x_0_0_0_xz_yy_z_zz = buffer_1000_ddpd[287];

    auto g_x_0_0_0_xz_yz_x_xx = buffer_1000_ddpd[288];

    auto g_x_0_0_0_xz_yz_x_xy = buffer_1000_ddpd[289];

    auto g_x_0_0_0_xz_yz_x_xz = buffer_1000_ddpd[290];

    auto g_x_0_0_0_xz_yz_x_yy = buffer_1000_ddpd[291];

    auto g_x_0_0_0_xz_yz_x_yz = buffer_1000_ddpd[292];

    auto g_x_0_0_0_xz_yz_x_zz = buffer_1000_ddpd[293];

    auto g_x_0_0_0_xz_yz_y_xx = buffer_1000_ddpd[294];

    auto g_x_0_0_0_xz_yz_y_xy = buffer_1000_ddpd[295];

    auto g_x_0_0_0_xz_yz_y_xz = buffer_1000_ddpd[296];

    auto g_x_0_0_0_xz_yz_y_yy = buffer_1000_ddpd[297];

    auto g_x_0_0_0_xz_yz_y_yz = buffer_1000_ddpd[298];

    auto g_x_0_0_0_xz_yz_y_zz = buffer_1000_ddpd[299];

    auto g_x_0_0_0_xz_yz_z_xx = buffer_1000_ddpd[300];

    auto g_x_0_0_0_xz_yz_z_xy = buffer_1000_ddpd[301];

    auto g_x_0_0_0_xz_yz_z_xz = buffer_1000_ddpd[302];

    auto g_x_0_0_0_xz_yz_z_yy = buffer_1000_ddpd[303];

    auto g_x_0_0_0_xz_yz_z_yz = buffer_1000_ddpd[304];

    auto g_x_0_0_0_xz_yz_z_zz = buffer_1000_ddpd[305];

    auto g_x_0_0_0_xz_zz_x_xx = buffer_1000_ddpd[306];

    auto g_x_0_0_0_xz_zz_x_xy = buffer_1000_ddpd[307];

    auto g_x_0_0_0_xz_zz_x_xz = buffer_1000_ddpd[308];

    auto g_x_0_0_0_xz_zz_x_yy = buffer_1000_ddpd[309];

    auto g_x_0_0_0_xz_zz_x_yz = buffer_1000_ddpd[310];

    auto g_x_0_0_0_xz_zz_x_zz = buffer_1000_ddpd[311];

    auto g_x_0_0_0_xz_zz_y_xx = buffer_1000_ddpd[312];

    auto g_x_0_0_0_xz_zz_y_xy = buffer_1000_ddpd[313];

    auto g_x_0_0_0_xz_zz_y_xz = buffer_1000_ddpd[314];

    auto g_x_0_0_0_xz_zz_y_yy = buffer_1000_ddpd[315];

    auto g_x_0_0_0_xz_zz_y_yz = buffer_1000_ddpd[316];

    auto g_x_0_0_0_xz_zz_y_zz = buffer_1000_ddpd[317];

    auto g_x_0_0_0_xz_zz_z_xx = buffer_1000_ddpd[318];

    auto g_x_0_0_0_xz_zz_z_xy = buffer_1000_ddpd[319];

    auto g_x_0_0_0_xz_zz_z_xz = buffer_1000_ddpd[320];

    auto g_x_0_0_0_xz_zz_z_yy = buffer_1000_ddpd[321];

    auto g_x_0_0_0_xz_zz_z_yz = buffer_1000_ddpd[322];

    auto g_x_0_0_0_xz_zz_z_zz = buffer_1000_ddpd[323];

    auto g_x_0_0_0_yy_xx_x_xx = buffer_1000_ddpd[324];

    auto g_x_0_0_0_yy_xx_x_xy = buffer_1000_ddpd[325];

    auto g_x_0_0_0_yy_xx_x_xz = buffer_1000_ddpd[326];

    auto g_x_0_0_0_yy_xx_x_yy = buffer_1000_ddpd[327];

    auto g_x_0_0_0_yy_xx_x_yz = buffer_1000_ddpd[328];

    auto g_x_0_0_0_yy_xx_x_zz = buffer_1000_ddpd[329];

    auto g_x_0_0_0_yy_xx_y_xx = buffer_1000_ddpd[330];

    auto g_x_0_0_0_yy_xx_y_xy = buffer_1000_ddpd[331];

    auto g_x_0_0_0_yy_xx_y_xz = buffer_1000_ddpd[332];

    auto g_x_0_0_0_yy_xx_y_yy = buffer_1000_ddpd[333];

    auto g_x_0_0_0_yy_xx_y_yz = buffer_1000_ddpd[334];

    auto g_x_0_0_0_yy_xx_y_zz = buffer_1000_ddpd[335];

    auto g_x_0_0_0_yy_xx_z_xx = buffer_1000_ddpd[336];

    auto g_x_0_0_0_yy_xx_z_xy = buffer_1000_ddpd[337];

    auto g_x_0_0_0_yy_xx_z_xz = buffer_1000_ddpd[338];

    auto g_x_0_0_0_yy_xx_z_yy = buffer_1000_ddpd[339];

    auto g_x_0_0_0_yy_xx_z_yz = buffer_1000_ddpd[340];

    auto g_x_0_0_0_yy_xx_z_zz = buffer_1000_ddpd[341];

    auto g_x_0_0_0_yy_xy_x_xx = buffer_1000_ddpd[342];

    auto g_x_0_0_0_yy_xy_x_xy = buffer_1000_ddpd[343];

    auto g_x_0_0_0_yy_xy_x_xz = buffer_1000_ddpd[344];

    auto g_x_0_0_0_yy_xy_x_yy = buffer_1000_ddpd[345];

    auto g_x_0_0_0_yy_xy_x_yz = buffer_1000_ddpd[346];

    auto g_x_0_0_0_yy_xy_x_zz = buffer_1000_ddpd[347];

    auto g_x_0_0_0_yy_xy_y_xx = buffer_1000_ddpd[348];

    auto g_x_0_0_0_yy_xy_y_xy = buffer_1000_ddpd[349];

    auto g_x_0_0_0_yy_xy_y_xz = buffer_1000_ddpd[350];

    auto g_x_0_0_0_yy_xy_y_yy = buffer_1000_ddpd[351];

    auto g_x_0_0_0_yy_xy_y_yz = buffer_1000_ddpd[352];

    auto g_x_0_0_0_yy_xy_y_zz = buffer_1000_ddpd[353];

    auto g_x_0_0_0_yy_xy_z_xx = buffer_1000_ddpd[354];

    auto g_x_0_0_0_yy_xy_z_xy = buffer_1000_ddpd[355];

    auto g_x_0_0_0_yy_xy_z_xz = buffer_1000_ddpd[356];

    auto g_x_0_0_0_yy_xy_z_yy = buffer_1000_ddpd[357];

    auto g_x_0_0_0_yy_xy_z_yz = buffer_1000_ddpd[358];

    auto g_x_0_0_0_yy_xy_z_zz = buffer_1000_ddpd[359];

    auto g_x_0_0_0_yy_xz_x_xx = buffer_1000_ddpd[360];

    auto g_x_0_0_0_yy_xz_x_xy = buffer_1000_ddpd[361];

    auto g_x_0_0_0_yy_xz_x_xz = buffer_1000_ddpd[362];

    auto g_x_0_0_0_yy_xz_x_yy = buffer_1000_ddpd[363];

    auto g_x_0_0_0_yy_xz_x_yz = buffer_1000_ddpd[364];

    auto g_x_0_0_0_yy_xz_x_zz = buffer_1000_ddpd[365];

    auto g_x_0_0_0_yy_xz_y_xx = buffer_1000_ddpd[366];

    auto g_x_0_0_0_yy_xz_y_xy = buffer_1000_ddpd[367];

    auto g_x_0_0_0_yy_xz_y_xz = buffer_1000_ddpd[368];

    auto g_x_0_0_0_yy_xz_y_yy = buffer_1000_ddpd[369];

    auto g_x_0_0_0_yy_xz_y_yz = buffer_1000_ddpd[370];

    auto g_x_0_0_0_yy_xz_y_zz = buffer_1000_ddpd[371];

    auto g_x_0_0_0_yy_xz_z_xx = buffer_1000_ddpd[372];

    auto g_x_0_0_0_yy_xz_z_xy = buffer_1000_ddpd[373];

    auto g_x_0_0_0_yy_xz_z_xz = buffer_1000_ddpd[374];

    auto g_x_0_0_0_yy_xz_z_yy = buffer_1000_ddpd[375];

    auto g_x_0_0_0_yy_xz_z_yz = buffer_1000_ddpd[376];

    auto g_x_0_0_0_yy_xz_z_zz = buffer_1000_ddpd[377];

    auto g_x_0_0_0_yy_yy_x_xx = buffer_1000_ddpd[378];

    auto g_x_0_0_0_yy_yy_x_xy = buffer_1000_ddpd[379];

    auto g_x_0_0_0_yy_yy_x_xz = buffer_1000_ddpd[380];

    auto g_x_0_0_0_yy_yy_x_yy = buffer_1000_ddpd[381];

    auto g_x_0_0_0_yy_yy_x_yz = buffer_1000_ddpd[382];

    auto g_x_0_0_0_yy_yy_x_zz = buffer_1000_ddpd[383];

    auto g_x_0_0_0_yy_yy_y_xx = buffer_1000_ddpd[384];

    auto g_x_0_0_0_yy_yy_y_xy = buffer_1000_ddpd[385];

    auto g_x_0_0_0_yy_yy_y_xz = buffer_1000_ddpd[386];

    auto g_x_0_0_0_yy_yy_y_yy = buffer_1000_ddpd[387];

    auto g_x_0_0_0_yy_yy_y_yz = buffer_1000_ddpd[388];

    auto g_x_0_0_0_yy_yy_y_zz = buffer_1000_ddpd[389];

    auto g_x_0_0_0_yy_yy_z_xx = buffer_1000_ddpd[390];

    auto g_x_0_0_0_yy_yy_z_xy = buffer_1000_ddpd[391];

    auto g_x_0_0_0_yy_yy_z_xz = buffer_1000_ddpd[392];

    auto g_x_0_0_0_yy_yy_z_yy = buffer_1000_ddpd[393];

    auto g_x_0_0_0_yy_yy_z_yz = buffer_1000_ddpd[394];

    auto g_x_0_0_0_yy_yy_z_zz = buffer_1000_ddpd[395];

    auto g_x_0_0_0_yy_yz_x_xx = buffer_1000_ddpd[396];

    auto g_x_0_0_0_yy_yz_x_xy = buffer_1000_ddpd[397];

    auto g_x_0_0_0_yy_yz_x_xz = buffer_1000_ddpd[398];

    auto g_x_0_0_0_yy_yz_x_yy = buffer_1000_ddpd[399];

    auto g_x_0_0_0_yy_yz_x_yz = buffer_1000_ddpd[400];

    auto g_x_0_0_0_yy_yz_x_zz = buffer_1000_ddpd[401];

    auto g_x_0_0_0_yy_yz_y_xx = buffer_1000_ddpd[402];

    auto g_x_0_0_0_yy_yz_y_xy = buffer_1000_ddpd[403];

    auto g_x_0_0_0_yy_yz_y_xz = buffer_1000_ddpd[404];

    auto g_x_0_0_0_yy_yz_y_yy = buffer_1000_ddpd[405];

    auto g_x_0_0_0_yy_yz_y_yz = buffer_1000_ddpd[406];

    auto g_x_0_0_0_yy_yz_y_zz = buffer_1000_ddpd[407];

    auto g_x_0_0_0_yy_yz_z_xx = buffer_1000_ddpd[408];

    auto g_x_0_0_0_yy_yz_z_xy = buffer_1000_ddpd[409];

    auto g_x_0_0_0_yy_yz_z_xz = buffer_1000_ddpd[410];

    auto g_x_0_0_0_yy_yz_z_yy = buffer_1000_ddpd[411];

    auto g_x_0_0_0_yy_yz_z_yz = buffer_1000_ddpd[412];

    auto g_x_0_0_0_yy_yz_z_zz = buffer_1000_ddpd[413];

    auto g_x_0_0_0_yy_zz_x_xx = buffer_1000_ddpd[414];

    auto g_x_0_0_0_yy_zz_x_xy = buffer_1000_ddpd[415];

    auto g_x_0_0_0_yy_zz_x_xz = buffer_1000_ddpd[416];

    auto g_x_0_0_0_yy_zz_x_yy = buffer_1000_ddpd[417];

    auto g_x_0_0_0_yy_zz_x_yz = buffer_1000_ddpd[418];

    auto g_x_0_0_0_yy_zz_x_zz = buffer_1000_ddpd[419];

    auto g_x_0_0_0_yy_zz_y_xx = buffer_1000_ddpd[420];

    auto g_x_0_0_0_yy_zz_y_xy = buffer_1000_ddpd[421];

    auto g_x_0_0_0_yy_zz_y_xz = buffer_1000_ddpd[422];

    auto g_x_0_0_0_yy_zz_y_yy = buffer_1000_ddpd[423];

    auto g_x_0_0_0_yy_zz_y_yz = buffer_1000_ddpd[424];

    auto g_x_0_0_0_yy_zz_y_zz = buffer_1000_ddpd[425];

    auto g_x_0_0_0_yy_zz_z_xx = buffer_1000_ddpd[426];

    auto g_x_0_0_0_yy_zz_z_xy = buffer_1000_ddpd[427];

    auto g_x_0_0_0_yy_zz_z_xz = buffer_1000_ddpd[428];

    auto g_x_0_0_0_yy_zz_z_yy = buffer_1000_ddpd[429];

    auto g_x_0_0_0_yy_zz_z_yz = buffer_1000_ddpd[430];

    auto g_x_0_0_0_yy_zz_z_zz = buffer_1000_ddpd[431];

    auto g_x_0_0_0_yz_xx_x_xx = buffer_1000_ddpd[432];

    auto g_x_0_0_0_yz_xx_x_xy = buffer_1000_ddpd[433];

    auto g_x_0_0_0_yz_xx_x_xz = buffer_1000_ddpd[434];

    auto g_x_0_0_0_yz_xx_x_yy = buffer_1000_ddpd[435];

    auto g_x_0_0_0_yz_xx_x_yz = buffer_1000_ddpd[436];

    auto g_x_0_0_0_yz_xx_x_zz = buffer_1000_ddpd[437];

    auto g_x_0_0_0_yz_xx_y_xx = buffer_1000_ddpd[438];

    auto g_x_0_0_0_yz_xx_y_xy = buffer_1000_ddpd[439];

    auto g_x_0_0_0_yz_xx_y_xz = buffer_1000_ddpd[440];

    auto g_x_0_0_0_yz_xx_y_yy = buffer_1000_ddpd[441];

    auto g_x_0_0_0_yz_xx_y_yz = buffer_1000_ddpd[442];

    auto g_x_0_0_0_yz_xx_y_zz = buffer_1000_ddpd[443];

    auto g_x_0_0_0_yz_xx_z_xx = buffer_1000_ddpd[444];

    auto g_x_0_0_0_yz_xx_z_xy = buffer_1000_ddpd[445];

    auto g_x_0_0_0_yz_xx_z_xz = buffer_1000_ddpd[446];

    auto g_x_0_0_0_yz_xx_z_yy = buffer_1000_ddpd[447];

    auto g_x_0_0_0_yz_xx_z_yz = buffer_1000_ddpd[448];

    auto g_x_0_0_0_yz_xx_z_zz = buffer_1000_ddpd[449];

    auto g_x_0_0_0_yz_xy_x_xx = buffer_1000_ddpd[450];

    auto g_x_0_0_0_yz_xy_x_xy = buffer_1000_ddpd[451];

    auto g_x_0_0_0_yz_xy_x_xz = buffer_1000_ddpd[452];

    auto g_x_0_0_0_yz_xy_x_yy = buffer_1000_ddpd[453];

    auto g_x_0_0_0_yz_xy_x_yz = buffer_1000_ddpd[454];

    auto g_x_0_0_0_yz_xy_x_zz = buffer_1000_ddpd[455];

    auto g_x_0_0_0_yz_xy_y_xx = buffer_1000_ddpd[456];

    auto g_x_0_0_0_yz_xy_y_xy = buffer_1000_ddpd[457];

    auto g_x_0_0_0_yz_xy_y_xz = buffer_1000_ddpd[458];

    auto g_x_0_0_0_yz_xy_y_yy = buffer_1000_ddpd[459];

    auto g_x_0_0_0_yz_xy_y_yz = buffer_1000_ddpd[460];

    auto g_x_0_0_0_yz_xy_y_zz = buffer_1000_ddpd[461];

    auto g_x_0_0_0_yz_xy_z_xx = buffer_1000_ddpd[462];

    auto g_x_0_0_0_yz_xy_z_xy = buffer_1000_ddpd[463];

    auto g_x_0_0_0_yz_xy_z_xz = buffer_1000_ddpd[464];

    auto g_x_0_0_0_yz_xy_z_yy = buffer_1000_ddpd[465];

    auto g_x_0_0_0_yz_xy_z_yz = buffer_1000_ddpd[466];

    auto g_x_0_0_0_yz_xy_z_zz = buffer_1000_ddpd[467];

    auto g_x_0_0_0_yz_xz_x_xx = buffer_1000_ddpd[468];

    auto g_x_0_0_0_yz_xz_x_xy = buffer_1000_ddpd[469];

    auto g_x_0_0_0_yz_xz_x_xz = buffer_1000_ddpd[470];

    auto g_x_0_0_0_yz_xz_x_yy = buffer_1000_ddpd[471];

    auto g_x_0_0_0_yz_xz_x_yz = buffer_1000_ddpd[472];

    auto g_x_0_0_0_yz_xz_x_zz = buffer_1000_ddpd[473];

    auto g_x_0_0_0_yz_xz_y_xx = buffer_1000_ddpd[474];

    auto g_x_0_0_0_yz_xz_y_xy = buffer_1000_ddpd[475];

    auto g_x_0_0_0_yz_xz_y_xz = buffer_1000_ddpd[476];

    auto g_x_0_0_0_yz_xz_y_yy = buffer_1000_ddpd[477];

    auto g_x_0_0_0_yz_xz_y_yz = buffer_1000_ddpd[478];

    auto g_x_0_0_0_yz_xz_y_zz = buffer_1000_ddpd[479];

    auto g_x_0_0_0_yz_xz_z_xx = buffer_1000_ddpd[480];

    auto g_x_0_0_0_yz_xz_z_xy = buffer_1000_ddpd[481];

    auto g_x_0_0_0_yz_xz_z_xz = buffer_1000_ddpd[482];

    auto g_x_0_0_0_yz_xz_z_yy = buffer_1000_ddpd[483];

    auto g_x_0_0_0_yz_xz_z_yz = buffer_1000_ddpd[484];

    auto g_x_0_0_0_yz_xz_z_zz = buffer_1000_ddpd[485];

    auto g_x_0_0_0_yz_yy_x_xx = buffer_1000_ddpd[486];

    auto g_x_0_0_0_yz_yy_x_xy = buffer_1000_ddpd[487];

    auto g_x_0_0_0_yz_yy_x_xz = buffer_1000_ddpd[488];

    auto g_x_0_0_0_yz_yy_x_yy = buffer_1000_ddpd[489];

    auto g_x_0_0_0_yz_yy_x_yz = buffer_1000_ddpd[490];

    auto g_x_0_0_0_yz_yy_x_zz = buffer_1000_ddpd[491];

    auto g_x_0_0_0_yz_yy_y_xx = buffer_1000_ddpd[492];

    auto g_x_0_0_0_yz_yy_y_xy = buffer_1000_ddpd[493];

    auto g_x_0_0_0_yz_yy_y_xz = buffer_1000_ddpd[494];

    auto g_x_0_0_0_yz_yy_y_yy = buffer_1000_ddpd[495];

    auto g_x_0_0_0_yz_yy_y_yz = buffer_1000_ddpd[496];

    auto g_x_0_0_0_yz_yy_y_zz = buffer_1000_ddpd[497];

    auto g_x_0_0_0_yz_yy_z_xx = buffer_1000_ddpd[498];

    auto g_x_0_0_0_yz_yy_z_xy = buffer_1000_ddpd[499];

    auto g_x_0_0_0_yz_yy_z_xz = buffer_1000_ddpd[500];

    auto g_x_0_0_0_yz_yy_z_yy = buffer_1000_ddpd[501];

    auto g_x_0_0_0_yz_yy_z_yz = buffer_1000_ddpd[502];

    auto g_x_0_0_0_yz_yy_z_zz = buffer_1000_ddpd[503];

    auto g_x_0_0_0_yz_yz_x_xx = buffer_1000_ddpd[504];

    auto g_x_0_0_0_yz_yz_x_xy = buffer_1000_ddpd[505];

    auto g_x_0_0_0_yz_yz_x_xz = buffer_1000_ddpd[506];

    auto g_x_0_0_0_yz_yz_x_yy = buffer_1000_ddpd[507];

    auto g_x_0_0_0_yz_yz_x_yz = buffer_1000_ddpd[508];

    auto g_x_0_0_0_yz_yz_x_zz = buffer_1000_ddpd[509];

    auto g_x_0_0_0_yz_yz_y_xx = buffer_1000_ddpd[510];

    auto g_x_0_0_0_yz_yz_y_xy = buffer_1000_ddpd[511];

    auto g_x_0_0_0_yz_yz_y_xz = buffer_1000_ddpd[512];

    auto g_x_0_0_0_yz_yz_y_yy = buffer_1000_ddpd[513];

    auto g_x_0_0_0_yz_yz_y_yz = buffer_1000_ddpd[514];

    auto g_x_0_0_0_yz_yz_y_zz = buffer_1000_ddpd[515];

    auto g_x_0_0_0_yz_yz_z_xx = buffer_1000_ddpd[516];

    auto g_x_0_0_0_yz_yz_z_xy = buffer_1000_ddpd[517];

    auto g_x_0_0_0_yz_yz_z_xz = buffer_1000_ddpd[518];

    auto g_x_0_0_0_yz_yz_z_yy = buffer_1000_ddpd[519];

    auto g_x_0_0_0_yz_yz_z_yz = buffer_1000_ddpd[520];

    auto g_x_0_0_0_yz_yz_z_zz = buffer_1000_ddpd[521];

    auto g_x_0_0_0_yz_zz_x_xx = buffer_1000_ddpd[522];

    auto g_x_0_0_0_yz_zz_x_xy = buffer_1000_ddpd[523];

    auto g_x_0_0_0_yz_zz_x_xz = buffer_1000_ddpd[524];

    auto g_x_0_0_0_yz_zz_x_yy = buffer_1000_ddpd[525];

    auto g_x_0_0_0_yz_zz_x_yz = buffer_1000_ddpd[526];

    auto g_x_0_0_0_yz_zz_x_zz = buffer_1000_ddpd[527];

    auto g_x_0_0_0_yz_zz_y_xx = buffer_1000_ddpd[528];

    auto g_x_0_0_0_yz_zz_y_xy = buffer_1000_ddpd[529];

    auto g_x_0_0_0_yz_zz_y_xz = buffer_1000_ddpd[530];

    auto g_x_0_0_0_yz_zz_y_yy = buffer_1000_ddpd[531];

    auto g_x_0_0_0_yz_zz_y_yz = buffer_1000_ddpd[532];

    auto g_x_0_0_0_yz_zz_y_zz = buffer_1000_ddpd[533];

    auto g_x_0_0_0_yz_zz_z_xx = buffer_1000_ddpd[534];

    auto g_x_0_0_0_yz_zz_z_xy = buffer_1000_ddpd[535];

    auto g_x_0_0_0_yz_zz_z_xz = buffer_1000_ddpd[536];

    auto g_x_0_0_0_yz_zz_z_yy = buffer_1000_ddpd[537];

    auto g_x_0_0_0_yz_zz_z_yz = buffer_1000_ddpd[538];

    auto g_x_0_0_0_yz_zz_z_zz = buffer_1000_ddpd[539];

    auto g_x_0_0_0_zz_xx_x_xx = buffer_1000_ddpd[540];

    auto g_x_0_0_0_zz_xx_x_xy = buffer_1000_ddpd[541];

    auto g_x_0_0_0_zz_xx_x_xz = buffer_1000_ddpd[542];

    auto g_x_0_0_0_zz_xx_x_yy = buffer_1000_ddpd[543];

    auto g_x_0_0_0_zz_xx_x_yz = buffer_1000_ddpd[544];

    auto g_x_0_0_0_zz_xx_x_zz = buffer_1000_ddpd[545];

    auto g_x_0_0_0_zz_xx_y_xx = buffer_1000_ddpd[546];

    auto g_x_0_0_0_zz_xx_y_xy = buffer_1000_ddpd[547];

    auto g_x_0_0_0_zz_xx_y_xz = buffer_1000_ddpd[548];

    auto g_x_0_0_0_zz_xx_y_yy = buffer_1000_ddpd[549];

    auto g_x_0_0_0_zz_xx_y_yz = buffer_1000_ddpd[550];

    auto g_x_0_0_0_zz_xx_y_zz = buffer_1000_ddpd[551];

    auto g_x_0_0_0_zz_xx_z_xx = buffer_1000_ddpd[552];

    auto g_x_0_0_0_zz_xx_z_xy = buffer_1000_ddpd[553];

    auto g_x_0_0_0_zz_xx_z_xz = buffer_1000_ddpd[554];

    auto g_x_0_0_0_zz_xx_z_yy = buffer_1000_ddpd[555];

    auto g_x_0_0_0_zz_xx_z_yz = buffer_1000_ddpd[556];

    auto g_x_0_0_0_zz_xx_z_zz = buffer_1000_ddpd[557];

    auto g_x_0_0_0_zz_xy_x_xx = buffer_1000_ddpd[558];

    auto g_x_0_0_0_zz_xy_x_xy = buffer_1000_ddpd[559];

    auto g_x_0_0_0_zz_xy_x_xz = buffer_1000_ddpd[560];

    auto g_x_0_0_0_zz_xy_x_yy = buffer_1000_ddpd[561];

    auto g_x_0_0_0_zz_xy_x_yz = buffer_1000_ddpd[562];

    auto g_x_0_0_0_zz_xy_x_zz = buffer_1000_ddpd[563];

    auto g_x_0_0_0_zz_xy_y_xx = buffer_1000_ddpd[564];

    auto g_x_0_0_0_zz_xy_y_xy = buffer_1000_ddpd[565];

    auto g_x_0_0_0_zz_xy_y_xz = buffer_1000_ddpd[566];

    auto g_x_0_0_0_zz_xy_y_yy = buffer_1000_ddpd[567];

    auto g_x_0_0_0_zz_xy_y_yz = buffer_1000_ddpd[568];

    auto g_x_0_0_0_zz_xy_y_zz = buffer_1000_ddpd[569];

    auto g_x_0_0_0_zz_xy_z_xx = buffer_1000_ddpd[570];

    auto g_x_0_0_0_zz_xy_z_xy = buffer_1000_ddpd[571];

    auto g_x_0_0_0_zz_xy_z_xz = buffer_1000_ddpd[572];

    auto g_x_0_0_0_zz_xy_z_yy = buffer_1000_ddpd[573];

    auto g_x_0_0_0_zz_xy_z_yz = buffer_1000_ddpd[574];

    auto g_x_0_0_0_zz_xy_z_zz = buffer_1000_ddpd[575];

    auto g_x_0_0_0_zz_xz_x_xx = buffer_1000_ddpd[576];

    auto g_x_0_0_0_zz_xz_x_xy = buffer_1000_ddpd[577];

    auto g_x_0_0_0_zz_xz_x_xz = buffer_1000_ddpd[578];

    auto g_x_0_0_0_zz_xz_x_yy = buffer_1000_ddpd[579];

    auto g_x_0_0_0_zz_xz_x_yz = buffer_1000_ddpd[580];

    auto g_x_0_0_0_zz_xz_x_zz = buffer_1000_ddpd[581];

    auto g_x_0_0_0_zz_xz_y_xx = buffer_1000_ddpd[582];

    auto g_x_0_0_0_zz_xz_y_xy = buffer_1000_ddpd[583];

    auto g_x_0_0_0_zz_xz_y_xz = buffer_1000_ddpd[584];

    auto g_x_0_0_0_zz_xz_y_yy = buffer_1000_ddpd[585];

    auto g_x_0_0_0_zz_xz_y_yz = buffer_1000_ddpd[586];

    auto g_x_0_0_0_zz_xz_y_zz = buffer_1000_ddpd[587];

    auto g_x_0_0_0_zz_xz_z_xx = buffer_1000_ddpd[588];

    auto g_x_0_0_0_zz_xz_z_xy = buffer_1000_ddpd[589];

    auto g_x_0_0_0_zz_xz_z_xz = buffer_1000_ddpd[590];

    auto g_x_0_0_0_zz_xz_z_yy = buffer_1000_ddpd[591];

    auto g_x_0_0_0_zz_xz_z_yz = buffer_1000_ddpd[592];

    auto g_x_0_0_0_zz_xz_z_zz = buffer_1000_ddpd[593];

    auto g_x_0_0_0_zz_yy_x_xx = buffer_1000_ddpd[594];

    auto g_x_0_0_0_zz_yy_x_xy = buffer_1000_ddpd[595];

    auto g_x_0_0_0_zz_yy_x_xz = buffer_1000_ddpd[596];

    auto g_x_0_0_0_zz_yy_x_yy = buffer_1000_ddpd[597];

    auto g_x_0_0_0_zz_yy_x_yz = buffer_1000_ddpd[598];

    auto g_x_0_0_0_zz_yy_x_zz = buffer_1000_ddpd[599];

    auto g_x_0_0_0_zz_yy_y_xx = buffer_1000_ddpd[600];

    auto g_x_0_0_0_zz_yy_y_xy = buffer_1000_ddpd[601];

    auto g_x_0_0_0_zz_yy_y_xz = buffer_1000_ddpd[602];

    auto g_x_0_0_0_zz_yy_y_yy = buffer_1000_ddpd[603];

    auto g_x_0_0_0_zz_yy_y_yz = buffer_1000_ddpd[604];

    auto g_x_0_0_0_zz_yy_y_zz = buffer_1000_ddpd[605];

    auto g_x_0_0_0_zz_yy_z_xx = buffer_1000_ddpd[606];

    auto g_x_0_0_0_zz_yy_z_xy = buffer_1000_ddpd[607];

    auto g_x_0_0_0_zz_yy_z_xz = buffer_1000_ddpd[608];

    auto g_x_0_0_0_zz_yy_z_yy = buffer_1000_ddpd[609];

    auto g_x_0_0_0_zz_yy_z_yz = buffer_1000_ddpd[610];

    auto g_x_0_0_0_zz_yy_z_zz = buffer_1000_ddpd[611];

    auto g_x_0_0_0_zz_yz_x_xx = buffer_1000_ddpd[612];

    auto g_x_0_0_0_zz_yz_x_xy = buffer_1000_ddpd[613];

    auto g_x_0_0_0_zz_yz_x_xz = buffer_1000_ddpd[614];

    auto g_x_0_0_0_zz_yz_x_yy = buffer_1000_ddpd[615];

    auto g_x_0_0_0_zz_yz_x_yz = buffer_1000_ddpd[616];

    auto g_x_0_0_0_zz_yz_x_zz = buffer_1000_ddpd[617];

    auto g_x_0_0_0_zz_yz_y_xx = buffer_1000_ddpd[618];

    auto g_x_0_0_0_zz_yz_y_xy = buffer_1000_ddpd[619];

    auto g_x_0_0_0_zz_yz_y_xz = buffer_1000_ddpd[620];

    auto g_x_0_0_0_zz_yz_y_yy = buffer_1000_ddpd[621];

    auto g_x_0_0_0_zz_yz_y_yz = buffer_1000_ddpd[622];

    auto g_x_0_0_0_zz_yz_y_zz = buffer_1000_ddpd[623];

    auto g_x_0_0_0_zz_yz_z_xx = buffer_1000_ddpd[624];

    auto g_x_0_0_0_zz_yz_z_xy = buffer_1000_ddpd[625];

    auto g_x_0_0_0_zz_yz_z_xz = buffer_1000_ddpd[626];

    auto g_x_0_0_0_zz_yz_z_yy = buffer_1000_ddpd[627];

    auto g_x_0_0_0_zz_yz_z_yz = buffer_1000_ddpd[628];

    auto g_x_0_0_0_zz_yz_z_zz = buffer_1000_ddpd[629];

    auto g_x_0_0_0_zz_zz_x_xx = buffer_1000_ddpd[630];

    auto g_x_0_0_0_zz_zz_x_xy = buffer_1000_ddpd[631];

    auto g_x_0_0_0_zz_zz_x_xz = buffer_1000_ddpd[632];

    auto g_x_0_0_0_zz_zz_x_yy = buffer_1000_ddpd[633];

    auto g_x_0_0_0_zz_zz_x_yz = buffer_1000_ddpd[634];

    auto g_x_0_0_0_zz_zz_x_zz = buffer_1000_ddpd[635];

    auto g_x_0_0_0_zz_zz_y_xx = buffer_1000_ddpd[636];

    auto g_x_0_0_0_zz_zz_y_xy = buffer_1000_ddpd[637];

    auto g_x_0_0_0_zz_zz_y_xz = buffer_1000_ddpd[638];

    auto g_x_0_0_0_zz_zz_y_yy = buffer_1000_ddpd[639];

    auto g_x_0_0_0_zz_zz_y_yz = buffer_1000_ddpd[640];

    auto g_x_0_0_0_zz_zz_y_zz = buffer_1000_ddpd[641];

    auto g_x_0_0_0_zz_zz_z_xx = buffer_1000_ddpd[642];

    auto g_x_0_0_0_zz_zz_z_xy = buffer_1000_ddpd[643];

    auto g_x_0_0_0_zz_zz_z_xz = buffer_1000_ddpd[644];

    auto g_x_0_0_0_zz_zz_z_yy = buffer_1000_ddpd[645];

    auto g_x_0_0_0_zz_zz_z_yz = buffer_1000_ddpd[646];

    auto g_x_0_0_0_zz_zz_z_zz = buffer_1000_ddpd[647];

    auto g_y_0_0_0_xx_xx_x_xx = buffer_1000_ddpd[648];

    auto g_y_0_0_0_xx_xx_x_xy = buffer_1000_ddpd[649];

    auto g_y_0_0_0_xx_xx_x_xz = buffer_1000_ddpd[650];

    auto g_y_0_0_0_xx_xx_x_yy = buffer_1000_ddpd[651];

    auto g_y_0_0_0_xx_xx_x_yz = buffer_1000_ddpd[652];

    auto g_y_0_0_0_xx_xx_x_zz = buffer_1000_ddpd[653];

    auto g_y_0_0_0_xx_xx_y_xx = buffer_1000_ddpd[654];

    auto g_y_0_0_0_xx_xx_y_xy = buffer_1000_ddpd[655];

    auto g_y_0_0_0_xx_xx_y_xz = buffer_1000_ddpd[656];

    auto g_y_0_0_0_xx_xx_y_yy = buffer_1000_ddpd[657];

    auto g_y_0_0_0_xx_xx_y_yz = buffer_1000_ddpd[658];

    auto g_y_0_0_0_xx_xx_y_zz = buffer_1000_ddpd[659];

    auto g_y_0_0_0_xx_xx_z_xx = buffer_1000_ddpd[660];

    auto g_y_0_0_0_xx_xx_z_xy = buffer_1000_ddpd[661];

    auto g_y_0_0_0_xx_xx_z_xz = buffer_1000_ddpd[662];

    auto g_y_0_0_0_xx_xx_z_yy = buffer_1000_ddpd[663];

    auto g_y_0_0_0_xx_xx_z_yz = buffer_1000_ddpd[664];

    auto g_y_0_0_0_xx_xx_z_zz = buffer_1000_ddpd[665];

    auto g_y_0_0_0_xx_xy_x_xx = buffer_1000_ddpd[666];

    auto g_y_0_0_0_xx_xy_x_xy = buffer_1000_ddpd[667];

    auto g_y_0_0_0_xx_xy_x_xz = buffer_1000_ddpd[668];

    auto g_y_0_0_0_xx_xy_x_yy = buffer_1000_ddpd[669];

    auto g_y_0_0_0_xx_xy_x_yz = buffer_1000_ddpd[670];

    auto g_y_0_0_0_xx_xy_x_zz = buffer_1000_ddpd[671];

    auto g_y_0_0_0_xx_xy_y_xx = buffer_1000_ddpd[672];

    auto g_y_0_0_0_xx_xy_y_xy = buffer_1000_ddpd[673];

    auto g_y_0_0_0_xx_xy_y_xz = buffer_1000_ddpd[674];

    auto g_y_0_0_0_xx_xy_y_yy = buffer_1000_ddpd[675];

    auto g_y_0_0_0_xx_xy_y_yz = buffer_1000_ddpd[676];

    auto g_y_0_0_0_xx_xy_y_zz = buffer_1000_ddpd[677];

    auto g_y_0_0_0_xx_xy_z_xx = buffer_1000_ddpd[678];

    auto g_y_0_0_0_xx_xy_z_xy = buffer_1000_ddpd[679];

    auto g_y_0_0_0_xx_xy_z_xz = buffer_1000_ddpd[680];

    auto g_y_0_0_0_xx_xy_z_yy = buffer_1000_ddpd[681];

    auto g_y_0_0_0_xx_xy_z_yz = buffer_1000_ddpd[682];

    auto g_y_0_0_0_xx_xy_z_zz = buffer_1000_ddpd[683];

    auto g_y_0_0_0_xx_xz_x_xx = buffer_1000_ddpd[684];

    auto g_y_0_0_0_xx_xz_x_xy = buffer_1000_ddpd[685];

    auto g_y_0_0_0_xx_xz_x_xz = buffer_1000_ddpd[686];

    auto g_y_0_0_0_xx_xz_x_yy = buffer_1000_ddpd[687];

    auto g_y_0_0_0_xx_xz_x_yz = buffer_1000_ddpd[688];

    auto g_y_0_0_0_xx_xz_x_zz = buffer_1000_ddpd[689];

    auto g_y_0_0_0_xx_xz_y_xx = buffer_1000_ddpd[690];

    auto g_y_0_0_0_xx_xz_y_xy = buffer_1000_ddpd[691];

    auto g_y_0_0_0_xx_xz_y_xz = buffer_1000_ddpd[692];

    auto g_y_0_0_0_xx_xz_y_yy = buffer_1000_ddpd[693];

    auto g_y_0_0_0_xx_xz_y_yz = buffer_1000_ddpd[694];

    auto g_y_0_0_0_xx_xz_y_zz = buffer_1000_ddpd[695];

    auto g_y_0_0_0_xx_xz_z_xx = buffer_1000_ddpd[696];

    auto g_y_0_0_0_xx_xz_z_xy = buffer_1000_ddpd[697];

    auto g_y_0_0_0_xx_xz_z_xz = buffer_1000_ddpd[698];

    auto g_y_0_0_0_xx_xz_z_yy = buffer_1000_ddpd[699];

    auto g_y_0_0_0_xx_xz_z_yz = buffer_1000_ddpd[700];

    auto g_y_0_0_0_xx_xz_z_zz = buffer_1000_ddpd[701];

    auto g_y_0_0_0_xx_yy_x_xx = buffer_1000_ddpd[702];

    auto g_y_0_0_0_xx_yy_x_xy = buffer_1000_ddpd[703];

    auto g_y_0_0_0_xx_yy_x_xz = buffer_1000_ddpd[704];

    auto g_y_0_0_0_xx_yy_x_yy = buffer_1000_ddpd[705];

    auto g_y_0_0_0_xx_yy_x_yz = buffer_1000_ddpd[706];

    auto g_y_0_0_0_xx_yy_x_zz = buffer_1000_ddpd[707];

    auto g_y_0_0_0_xx_yy_y_xx = buffer_1000_ddpd[708];

    auto g_y_0_0_0_xx_yy_y_xy = buffer_1000_ddpd[709];

    auto g_y_0_0_0_xx_yy_y_xz = buffer_1000_ddpd[710];

    auto g_y_0_0_0_xx_yy_y_yy = buffer_1000_ddpd[711];

    auto g_y_0_0_0_xx_yy_y_yz = buffer_1000_ddpd[712];

    auto g_y_0_0_0_xx_yy_y_zz = buffer_1000_ddpd[713];

    auto g_y_0_0_0_xx_yy_z_xx = buffer_1000_ddpd[714];

    auto g_y_0_0_0_xx_yy_z_xy = buffer_1000_ddpd[715];

    auto g_y_0_0_0_xx_yy_z_xz = buffer_1000_ddpd[716];

    auto g_y_0_0_0_xx_yy_z_yy = buffer_1000_ddpd[717];

    auto g_y_0_0_0_xx_yy_z_yz = buffer_1000_ddpd[718];

    auto g_y_0_0_0_xx_yy_z_zz = buffer_1000_ddpd[719];

    auto g_y_0_0_0_xx_yz_x_xx = buffer_1000_ddpd[720];

    auto g_y_0_0_0_xx_yz_x_xy = buffer_1000_ddpd[721];

    auto g_y_0_0_0_xx_yz_x_xz = buffer_1000_ddpd[722];

    auto g_y_0_0_0_xx_yz_x_yy = buffer_1000_ddpd[723];

    auto g_y_0_0_0_xx_yz_x_yz = buffer_1000_ddpd[724];

    auto g_y_0_0_0_xx_yz_x_zz = buffer_1000_ddpd[725];

    auto g_y_0_0_0_xx_yz_y_xx = buffer_1000_ddpd[726];

    auto g_y_0_0_0_xx_yz_y_xy = buffer_1000_ddpd[727];

    auto g_y_0_0_0_xx_yz_y_xz = buffer_1000_ddpd[728];

    auto g_y_0_0_0_xx_yz_y_yy = buffer_1000_ddpd[729];

    auto g_y_0_0_0_xx_yz_y_yz = buffer_1000_ddpd[730];

    auto g_y_0_0_0_xx_yz_y_zz = buffer_1000_ddpd[731];

    auto g_y_0_0_0_xx_yz_z_xx = buffer_1000_ddpd[732];

    auto g_y_0_0_0_xx_yz_z_xy = buffer_1000_ddpd[733];

    auto g_y_0_0_0_xx_yz_z_xz = buffer_1000_ddpd[734];

    auto g_y_0_0_0_xx_yz_z_yy = buffer_1000_ddpd[735];

    auto g_y_0_0_0_xx_yz_z_yz = buffer_1000_ddpd[736];

    auto g_y_0_0_0_xx_yz_z_zz = buffer_1000_ddpd[737];

    auto g_y_0_0_0_xx_zz_x_xx = buffer_1000_ddpd[738];

    auto g_y_0_0_0_xx_zz_x_xy = buffer_1000_ddpd[739];

    auto g_y_0_0_0_xx_zz_x_xz = buffer_1000_ddpd[740];

    auto g_y_0_0_0_xx_zz_x_yy = buffer_1000_ddpd[741];

    auto g_y_0_0_0_xx_zz_x_yz = buffer_1000_ddpd[742];

    auto g_y_0_0_0_xx_zz_x_zz = buffer_1000_ddpd[743];

    auto g_y_0_0_0_xx_zz_y_xx = buffer_1000_ddpd[744];

    auto g_y_0_0_0_xx_zz_y_xy = buffer_1000_ddpd[745];

    auto g_y_0_0_0_xx_zz_y_xz = buffer_1000_ddpd[746];

    auto g_y_0_0_0_xx_zz_y_yy = buffer_1000_ddpd[747];

    auto g_y_0_0_0_xx_zz_y_yz = buffer_1000_ddpd[748];

    auto g_y_0_0_0_xx_zz_y_zz = buffer_1000_ddpd[749];

    auto g_y_0_0_0_xx_zz_z_xx = buffer_1000_ddpd[750];

    auto g_y_0_0_0_xx_zz_z_xy = buffer_1000_ddpd[751];

    auto g_y_0_0_0_xx_zz_z_xz = buffer_1000_ddpd[752];

    auto g_y_0_0_0_xx_zz_z_yy = buffer_1000_ddpd[753];

    auto g_y_0_0_0_xx_zz_z_yz = buffer_1000_ddpd[754];

    auto g_y_0_0_0_xx_zz_z_zz = buffer_1000_ddpd[755];

    auto g_y_0_0_0_xy_xx_x_xx = buffer_1000_ddpd[756];

    auto g_y_0_0_0_xy_xx_x_xy = buffer_1000_ddpd[757];

    auto g_y_0_0_0_xy_xx_x_xz = buffer_1000_ddpd[758];

    auto g_y_0_0_0_xy_xx_x_yy = buffer_1000_ddpd[759];

    auto g_y_0_0_0_xy_xx_x_yz = buffer_1000_ddpd[760];

    auto g_y_0_0_0_xy_xx_x_zz = buffer_1000_ddpd[761];

    auto g_y_0_0_0_xy_xx_y_xx = buffer_1000_ddpd[762];

    auto g_y_0_0_0_xy_xx_y_xy = buffer_1000_ddpd[763];

    auto g_y_0_0_0_xy_xx_y_xz = buffer_1000_ddpd[764];

    auto g_y_0_0_0_xy_xx_y_yy = buffer_1000_ddpd[765];

    auto g_y_0_0_0_xy_xx_y_yz = buffer_1000_ddpd[766];

    auto g_y_0_0_0_xy_xx_y_zz = buffer_1000_ddpd[767];

    auto g_y_0_0_0_xy_xx_z_xx = buffer_1000_ddpd[768];

    auto g_y_0_0_0_xy_xx_z_xy = buffer_1000_ddpd[769];

    auto g_y_0_0_0_xy_xx_z_xz = buffer_1000_ddpd[770];

    auto g_y_0_0_0_xy_xx_z_yy = buffer_1000_ddpd[771];

    auto g_y_0_0_0_xy_xx_z_yz = buffer_1000_ddpd[772];

    auto g_y_0_0_0_xy_xx_z_zz = buffer_1000_ddpd[773];

    auto g_y_0_0_0_xy_xy_x_xx = buffer_1000_ddpd[774];

    auto g_y_0_0_0_xy_xy_x_xy = buffer_1000_ddpd[775];

    auto g_y_0_0_0_xy_xy_x_xz = buffer_1000_ddpd[776];

    auto g_y_0_0_0_xy_xy_x_yy = buffer_1000_ddpd[777];

    auto g_y_0_0_0_xy_xy_x_yz = buffer_1000_ddpd[778];

    auto g_y_0_0_0_xy_xy_x_zz = buffer_1000_ddpd[779];

    auto g_y_0_0_0_xy_xy_y_xx = buffer_1000_ddpd[780];

    auto g_y_0_0_0_xy_xy_y_xy = buffer_1000_ddpd[781];

    auto g_y_0_0_0_xy_xy_y_xz = buffer_1000_ddpd[782];

    auto g_y_0_0_0_xy_xy_y_yy = buffer_1000_ddpd[783];

    auto g_y_0_0_0_xy_xy_y_yz = buffer_1000_ddpd[784];

    auto g_y_0_0_0_xy_xy_y_zz = buffer_1000_ddpd[785];

    auto g_y_0_0_0_xy_xy_z_xx = buffer_1000_ddpd[786];

    auto g_y_0_0_0_xy_xy_z_xy = buffer_1000_ddpd[787];

    auto g_y_0_0_0_xy_xy_z_xz = buffer_1000_ddpd[788];

    auto g_y_0_0_0_xy_xy_z_yy = buffer_1000_ddpd[789];

    auto g_y_0_0_0_xy_xy_z_yz = buffer_1000_ddpd[790];

    auto g_y_0_0_0_xy_xy_z_zz = buffer_1000_ddpd[791];

    auto g_y_0_0_0_xy_xz_x_xx = buffer_1000_ddpd[792];

    auto g_y_0_0_0_xy_xz_x_xy = buffer_1000_ddpd[793];

    auto g_y_0_0_0_xy_xz_x_xz = buffer_1000_ddpd[794];

    auto g_y_0_0_0_xy_xz_x_yy = buffer_1000_ddpd[795];

    auto g_y_0_0_0_xy_xz_x_yz = buffer_1000_ddpd[796];

    auto g_y_0_0_0_xy_xz_x_zz = buffer_1000_ddpd[797];

    auto g_y_0_0_0_xy_xz_y_xx = buffer_1000_ddpd[798];

    auto g_y_0_0_0_xy_xz_y_xy = buffer_1000_ddpd[799];

    auto g_y_0_0_0_xy_xz_y_xz = buffer_1000_ddpd[800];

    auto g_y_0_0_0_xy_xz_y_yy = buffer_1000_ddpd[801];

    auto g_y_0_0_0_xy_xz_y_yz = buffer_1000_ddpd[802];

    auto g_y_0_0_0_xy_xz_y_zz = buffer_1000_ddpd[803];

    auto g_y_0_0_0_xy_xz_z_xx = buffer_1000_ddpd[804];

    auto g_y_0_0_0_xy_xz_z_xy = buffer_1000_ddpd[805];

    auto g_y_0_0_0_xy_xz_z_xz = buffer_1000_ddpd[806];

    auto g_y_0_0_0_xy_xz_z_yy = buffer_1000_ddpd[807];

    auto g_y_0_0_0_xy_xz_z_yz = buffer_1000_ddpd[808];

    auto g_y_0_0_0_xy_xz_z_zz = buffer_1000_ddpd[809];

    auto g_y_0_0_0_xy_yy_x_xx = buffer_1000_ddpd[810];

    auto g_y_0_0_0_xy_yy_x_xy = buffer_1000_ddpd[811];

    auto g_y_0_0_0_xy_yy_x_xz = buffer_1000_ddpd[812];

    auto g_y_0_0_0_xy_yy_x_yy = buffer_1000_ddpd[813];

    auto g_y_0_0_0_xy_yy_x_yz = buffer_1000_ddpd[814];

    auto g_y_0_0_0_xy_yy_x_zz = buffer_1000_ddpd[815];

    auto g_y_0_0_0_xy_yy_y_xx = buffer_1000_ddpd[816];

    auto g_y_0_0_0_xy_yy_y_xy = buffer_1000_ddpd[817];

    auto g_y_0_0_0_xy_yy_y_xz = buffer_1000_ddpd[818];

    auto g_y_0_0_0_xy_yy_y_yy = buffer_1000_ddpd[819];

    auto g_y_0_0_0_xy_yy_y_yz = buffer_1000_ddpd[820];

    auto g_y_0_0_0_xy_yy_y_zz = buffer_1000_ddpd[821];

    auto g_y_0_0_0_xy_yy_z_xx = buffer_1000_ddpd[822];

    auto g_y_0_0_0_xy_yy_z_xy = buffer_1000_ddpd[823];

    auto g_y_0_0_0_xy_yy_z_xz = buffer_1000_ddpd[824];

    auto g_y_0_0_0_xy_yy_z_yy = buffer_1000_ddpd[825];

    auto g_y_0_0_0_xy_yy_z_yz = buffer_1000_ddpd[826];

    auto g_y_0_0_0_xy_yy_z_zz = buffer_1000_ddpd[827];

    auto g_y_0_0_0_xy_yz_x_xx = buffer_1000_ddpd[828];

    auto g_y_0_0_0_xy_yz_x_xy = buffer_1000_ddpd[829];

    auto g_y_0_0_0_xy_yz_x_xz = buffer_1000_ddpd[830];

    auto g_y_0_0_0_xy_yz_x_yy = buffer_1000_ddpd[831];

    auto g_y_0_0_0_xy_yz_x_yz = buffer_1000_ddpd[832];

    auto g_y_0_0_0_xy_yz_x_zz = buffer_1000_ddpd[833];

    auto g_y_0_0_0_xy_yz_y_xx = buffer_1000_ddpd[834];

    auto g_y_0_0_0_xy_yz_y_xy = buffer_1000_ddpd[835];

    auto g_y_0_0_0_xy_yz_y_xz = buffer_1000_ddpd[836];

    auto g_y_0_0_0_xy_yz_y_yy = buffer_1000_ddpd[837];

    auto g_y_0_0_0_xy_yz_y_yz = buffer_1000_ddpd[838];

    auto g_y_0_0_0_xy_yz_y_zz = buffer_1000_ddpd[839];

    auto g_y_0_0_0_xy_yz_z_xx = buffer_1000_ddpd[840];

    auto g_y_0_0_0_xy_yz_z_xy = buffer_1000_ddpd[841];

    auto g_y_0_0_0_xy_yz_z_xz = buffer_1000_ddpd[842];

    auto g_y_0_0_0_xy_yz_z_yy = buffer_1000_ddpd[843];

    auto g_y_0_0_0_xy_yz_z_yz = buffer_1000_ddpd[844];

    auto g_y_0_0_0_xy_yz_z_zz = buffer_1000_ddpd[845];

    auto g_y_0_0_0_xy_zz_x_xx = buffer_1000_ddpd[846];

    auto g_y_0_0_0_xy_zz_x_xy = buffer_1000_ddpd[847];

    auto g_y_0_0_0_xy_zz_x_xz = buffer_1000_ddpd[848];

    auto g_y_0_0_0_xy_zz_x_yy = buffer_1000_ddpd[849];

    auto g_y_0_0_0_xy_zz_x_yz = buffer_1000_ddpd[850];

    auto g_y_0_0_0_xy_zz_x_zz = buffer_1000_ddpd[851];

    auto g_y_0_0_0_xy_zz_y_xx = buffer_1000_ddpd[852];

    auto g_y_0_0_0_xy_zz_y_xy = buffer_1000_ddpd[853];

    auto g_y_0_0_0_xy_zz_y_xz = buffer_1000_ddpd[854];

    auto g_y_0_0_0_xy_zz_y_yy = buffer_1000_ddpd[855];

    auto g_y_0_0_0_xy_zz_y_yz = buffer_1000_ddpd[856];

    auto g_y_0_0_0_xy_zz_y_zz = buffer_1000_ddpd[857];

    auto g_y_0_0_0_xy_zz_z_xx = buffer_1000_ddpd[858];

    auto g_y_0_0_0_xy_zz_z_xy = buffer_1000_ddpd[859];

    auto g_y_0_0_0_xy_zz_z_xz = buffer_1000_ddpd[860];

    auto g_y_0_0_0_xy_zz_z_yy = buffer_1000_ddpd[861];

    auto g_y_0_0_0_xy_zz_z_yz = buffer_1000_ddpd[862];

    auto g_y_0_0_0_xy_zz_z_zz = buffer_1000_ddpd[863];

    auto g_y_0_0_0_xz_xx_x_xx = buffer_1000_ddpd[864];

    auto g_y_0_0_0_xz_xx_x_xy = buffer_1000_ddpd[865];

    auto g_y_0_0_0_xz_xx_x_xz = buffer_1000_ddpd[866];

    auto g_y_0_0_0_xz_xx_x_yy = buffer_1000_ddpd[867];

    auto g_y_0_0_0_xz_xx_x_yz = buffer_1000_ddpd[868];

    auto g_y_0_0_0_xz_xx_x_zz = buffer_1000_ddpd[869];

    auto g_y_0_0_0_xz_xx_y_xx = buffer_1000_ddpd[870];

    auto g_y_0_0_0_xz_xx_y_xy = buffer_1000_ddpd[871];

    auto g_y_0_0_0_xz_xx_y_xz = buffer_1000_ddpd[872];

    auto g_y_0_0_0_xz_xx_y_yy = buffer_1000_ddpd[873];

    auto g_y_0_0_0_xz_xx_y_yz = buffer_1000_ddpd[874];

    auto g_y_0_0_0_xz_xx_y_zz = buffer_1000_ddpd[875];

    auto g_y_0_0_0_xz_xx_z_xx = buffer_1000_ddpd[876];

    auto g_y_0_0_0_xz_xx_z_xy = buffer_1000_ddpd[877];

    auto g_y_0_0_0_xz_xx_z_xz = buffer_1000_ddpd[878];

    auto g_y_0_0_0_xz_xx_z_yy = buffer_1000_ddpd[879];

    auto g_y_0_0_0_xz_xx_z_yz = buffer_1000_ddpd[880];

    auto g_y_0_0_0_xz_xx_z_zz = buffer_1000_ddpd[881];

    auto g_y_0_0_0_xz_xy_x_xx = buffer_1000_ddpd[882];

    auto g_y_0_0_0_xz_xy_x_xy = buffer_1000_ddpd[883];

    auto g_y_0_0_0_xz_xy_x_xz = buffer_1000_ddpd[884];

    auto g_y_0_0_0_xz_xy_x_yy = buffer_1000_ddpd[885];

    auto g_y_0_0_0_xz_xy_x_yz = buffer_1000_ddpd[886];

    auto g_y_0_0_0_xz_xy_x_zz = buffer_1000_ddpd[887];

    auto g_y_0_0_0_xz_xy_y_xx = buffer_1000_ddpd[888];

    auto g_y_0_0_0_xz_xy_y_xy = buffer_1000_ddpd[889];

    auto g_y_0_0_0_xz_xy_y_xz = buffer_1000_ddpd[890];

    auto g_y_0_0_0_xz_xy_y_yy = buffer_1000_ddpd[891];

    auto g_y_0_0_0_xz_xy_y_yz = buffer_1000_ddpd[892];

    auto g_y_0_0_0_xz_xy_y_zz = buffer_1000_ddpd[893];

    auto g_y_0_0_0_xz_xy_z_xx = buffer_1000_ddpd[894];

    auto g_y_0_0_0_xz_xy_z_xy = buffer_1000_ddpd[895];

    auto g_y_0_0_0_xz_xy_z_xz = buffer_1000_ddpd[896];

    auto g_y_0_0_0_xz_xy_z_yy = buffer_1000_ddpd[897];

    auto g_y_0_0_0_xz_xy_z_yz = buffer_1000_ddpd[898];

    auto g_y_0_0_0_xz_xy_z_zz = buffer_1000_ddpd[899];

    auto g_y_0_0_0_xz_xz_x_xx = buffer_1000_ddpd[900];

    auto g_y_0_0_0_xz_xz_x_xy = buffer_1000_ddpd[901];

    auto g_y_0_0_0_xz_xz_x_xz = buffer_1000_ddpd[902];

    auto g_y_0_0_0_xz_xz_x_yy = buffer_1000_ddpd[903];

    auto g_y_0_0_0_xz_xz_x_yz = buffer_1000_ddpd[904];

    auto g_y_0_0_0_xz_xz_x_zz = buffer_1000_ddpd[905];

    auto g_y_0_0_0_xz_xz_y_xx = buffer_1000_ddpd[906];

    auto g_y_0_0_0_xz_xz_y_xy = buffer_1000_ddpd[907];

    auto g_y_0_0_0_xz_xz_y_xz = buffer_1000_ddpd[908];

    auto g_y_0_0_0_xz_xz_y_yy = buffer_1000_ddpd[909];

    auto g_y_0_0_0_xz_xz_y_yz = buffer_1000_ddpd[910];

    auto g_y_0_0_0_xz_xz_y_zz = buffer_1000_ddpd[911];

    auto g_y_0_0_0_xz_xz_z_xx = buffer_1000_ddpd[912];

    auto g_y_0_0_0_xz_xz_z_xy = buffer_1000_ddpd[913];

    auto g_y_0_0_0_xz_xz_z_xz = buffer_1000_ddpd[914];

    auto g_y_0_0_0_xz_xz_z_yy = buffer_1000_ddpd[915];

    auto g_y_0_0_0_xz_xz_z_yz = buffer_1000_ddpd[916];

    auto g_y_0_0_0_xz_xz_z_zz = buffer_1000_ddpd[917];

    auto g_y_0_0_0_xz_yy_x_xx = buffer_1000_ddpd[918];

    auto g_y_0_0_0_xz_yy_x_xy = buffer_1000_ddpd[919];

    auto g_y_0_0_0_xz_yy_x_xz = buffer_1000_ddpd[920];

    auto g_y_0_0_0_xz_yy_x_yy = buffer_1000_ddpd[921];

    auto g_y_0_0_0_xz_yy_x_yz = buffer_1000_ddpd[922];

    auto g_y_0_0_0_xz_yy_x_zz = buffer_1000_ddpd[923];

    auto g_y_0_0_0_xz_yy_y_xx = buffer_1000_ddpd[924];

    auto g_y_0_0_0_xz_yy_y_xy = buffer_1000_ddpd[925];

    auto g_y_0_0_0_xz_yy_y_xz = buffer_1000_ddpd[926];

    auto g_y_0_0_0_xz_yy_y_yy = buffer_1000_ddpd[927];

    auto g_y_0_0_0_xz_yy_y_yz = buffer_1000_ddpd[928];

    auto g_y_0_0_0_xz_yy_y_zz = buffer_1000_ddpd[929];

    auto g_y_0_0_0_xz_yy_z_xx = buffer_1000_ddpd[930];

    auto g_y_0_0_0_xz_yy_z_xy = buffer_1000_ddpd[931];

    auto g_y_0_0_0_xz_yy_z_xz = buffer_1000_ddpd[932];

    auto g_y_0_0_0_xz_yy_z_yy = buffer_1000_ddpd[933];

    auto g_y_0_0_0_xz_yy_z_yz = buffer_1000_ddpd[934];

    auto g_y_0_0_0_xz_yy_z_zz = buffer_1000_ddpd[935];

    auto g_y_0_0_0_xz_yz_x_xx = buffer_1000_ddpd[936];

    auto g_y_0_0_0_xz_yz_x_xy = buffer_1000_ddpd[937];

    auto g_y_0_0_0_xz_yz_x_xz = buffer_1000_ddpd[938];

    auto g_y_0_0_0_xz_yz_x_yy = buffer_1000_ddpd[939];

    auto g_y_0_0_0_xz_yz_x_yz = buffer_1000_ddpd[940];

    auto g_y_0_0_0_xz_yz_x_zz = buffer_1000_ddpd[941];

    auto g_y_0_0_0_xz_yz_y_xx = buffer_1000_ddpd[942];

    auto g_y_0_0_0_xz_yz_y_xy = buffer_1000_ddpd[943];

    auto g_y_0_0_0_xz_yz_y_xz = buffer_1000_ddpd[944];

    auto g_y_0_0_0_xz_yz_y_yy = buffer_1000_ddpd[945];

    auto g_y_0_0_0_xz_yz_y_yz = buffer_1000_ddpd[946];

    auto g_y_0_0_0_xz_yz_y_zz = buffer_1000_ddpd[947];

    auto g_y_0_0_0_xz_yz_z_xx = buffer_1000_ddpd[948];

    auto g_y_0_0_0_xz_yz_z_xy = buffer_1000_ddpd[949];

    auto g_y_0_0_0_xz_yz_z_xz = buffer_1000_ddpd[950];

    auto g_y_0_0_0_xz_yz_z_yy = buffer_1000_ddpd[951];

    auto g_y_0_0_0_xz_yz_z_yz = buffer_1000_ddpd[952];

    auto g_y_0_0_0_xz_yz_z_zz = buffer_1000_ddpd[953];

    auto g_y_0_0_0_xz_zz_x_xx = buffer_1000_ddpd[954];

    auto g_y_0_0_0_xz_zz_x_xy = buffer_1000_ddpd[955];

    auto g_y_0_0_0_xz_zz_x_xz = buffer_1000_ddpd[956];

    auto g_y_0_0_0_xz_zz_x_yy = buffer_1000_ddpd[957];

    auto g_y_0_0_0_xz_zz_x_yz = buffer_1000_ddpd[958];

    auto g_y_0_0_0_xz_zz_x_zz = buffer_1000_ddpd[959];

    auto g_y_0_0_0_xz_zz_y_xx = buffer_1000_ddpd[960];

    auto g_y_0_0_0_xz_zz_y_xy = buffer_1000_ddpd[961];

    auto g_y_0_0_0_xz_zz_y_xz = buffer_1000_ddpd[962];

    auto g_y_0_0_0_xz_zz_y_yy = buffer_1000_ddpd[963];

    auto g_y_0_0_0_xz_zz_y_yz = buffer_1000_ddpd[964];

    auto g_y_0_0_0_xz_zz_y_zz = buffer_1000_ddpd[965];

    auto g_y_0_0_0_xz_zz_z_xx = buffer_1000_ddpd[966];

    auto g_y_0_0_0_xz_zz_z_xy = buffer_1000_ddpd[967];

    auto g_y_0_0_0_xz_zz_z_xz = buffer_1000_ddpd[968];

    auto g_y_0_0_0_xz_zz_z_yy = buffer_1000_ddpd[969];

    auto g_y_0_0_0_xz_zz_z_yz = buffer_1000_ddpd[970];

    auto g_y_0_0_0_xz_zz_z_zz = buffer_1000_ddpd[971];

    auto g_y_0_0_0_yy_xx_x_xx = buffer_1000_ddpd[972];

    auto g_y_0_0_0_yy_xx_x_xy = buffer_1000_ddpd[973];

    auto g_y_0_0_0_yy_xx_x_xz = buffer_1000_ddpd[974];

    auto g_y_0_0_0_yy_xx_x_yy = buffer_1000_ddpd[975];

    auto g_y_0_0_0_yy_xx_x_yz = buffer_1000_ddpd[976];

    auto g_y_0_0_0_yy_xx_x_zz = buffer_1000_ddpd[977];

    auto g_y_0_0_0_yy_xx_y_xx = buffer_1000_ddpd[978];

    auto g_y_0_0_0_yy_xx_y_xy = buffer_1000_ddpd[979];

    auto g_y_0_0_0_yy_xx_y_xz = buffer_1000_ddpd[980];

    auto g_y_0_0_0_yy_xx_y_yy = buffer_1000_ddpd[981];

    auto g_y_0_0_0_yy_xx_y_yz = buffer_1000_ddpd[982];

    auto g_y_0_0_0_yy_xx_y_zz = buffer_1000_ddpd[983];

    auto g_y_0_0_0_yy_xx_z_xx = buffer_1000_ddpd[984];

    auto g_y_0_0_0_yy_xx_z_xy = buffer_1000_ddpd[985];

    auto g_y_0_0_0_yy_xx_z_xz = buffer_1000_ddpd[986];

    auto g_y_0_0_0_yy_xx_z_yy = buffer_1000_ddpd[987];

    auto g_y_0_0_0_yy_xx_z_yz = buffer_1000_ddpd[988];

    auto g_y_0_0_0_yy_xx_z_zz = buffer_1000_ddpd[989];

    auto g_y_0_0_0_yy_xy_x_xx = buffer_1000_ddpd[990];

    auto g_y_0_0_0_yy_xy_x_xy = buffer_1000_ddpd[991];

    auto g_y_0_0_0_yy_xy_x_xz = buffer_1000_ddpd[992];

    auto g_y_0_0_0_yy_xy_x_yy = buffer_1000_ddpd[993];

    auto g_y_0_0_0_yy_xy_x_yz = buffer_1000_ddpd[994];

    auto g_y_0_0_0_yy_xy_x_zz = buffer_1000_ddpd[995];

    auto g_y_0_0_0_yy_xy_y_xx = buffer_1000_ddpd[996];

    auto g_y_0_0_0_yy_xy_y_xy = buffer_1000_ddpd[997];

    auto g_y_0_0_0_yy_xy_y_xz = buffer_1000_ddpd[998];

    auto g_y_0_0_0_yy_xy_y_yy = buffer_1000_ddpd[999];

    auto g_y_0_0_0_yy_xy_y_yz = buffer_1000_ddpd[1000];

    auto g_y_0_0_0_yy_xy_y_zz = buffer_1000_ddpd[1001];

    auto g_y_0_0_0_yy_xy_z_xx = buffer_1000_ddpd[1002];

    auto g_y_0_0_0_yy_xy_z_xy = buffer_1000_ddpd[1003];

    auto g_y_0_0_0_yy_xy_z_xz = buffer_1000_ddpd[1004];

    auto g_y_0_0_0_yy_xy_z_yy = buffer_1000_ddpd[1005];

    auto g_y_0_0_0_yy_xy_z_yz = buffer_1000_ddpd[1006];

    auto g_y_0_0_0_yy_xy_z_zz = buffer_1000_ddpd[1007];

    auto g_y_0_0_0_yy_xz_x_xx = buffer_1000_ddpd[1008];

    auto g_y_0_0_0_yy_xz_x_xy = buffer_1000_ddpd[1009];

    auto g_y_0_0_0_yy_xz_x_xz = buffer_1000_ddpd[1010];

    auto g_y_0_0_0_yy_xz_x_yy = buffer_1000_ddpd[1011];

    auto g_y_0_0_0_yy_xz_x_yz = buffer_1000_ddpd[1012];

    auto g_y_0_0_0_yy_xz_x_zz = buffer_1000_ddpd[1013];

    auto g_y_0_0_0_yy_xz_y_xx = buffer_1000_ddpd[1014];

    auto g_y_0_0_0_yy_xz_y_xy = buffer_1000_ddpd[1015];

    auto g_y_0_0_0_yy_xz_y_xz = buffer_1000_ddpd[1016];

    auto g_y_0_0_0_yy_xz_y_yy = buffer_1000_ddpd[1017];

    auto g_y_0_0_0_yy_xz_y_yz = buffer_1000_ddpd[1018];

    auto g_y_0_0_0_yy_xz_y_zz = buffer_1000_ddpd[1019];

    auto g_y_0_0_0_yy_xz_z_xx = buffer_1000_ddpd[1020];

    auto g_y_0_0_0_yy_xz_z_xy = buffer_1000_ddpd[1021];

    auto g_y_0_0_0_yy_xz_z_xz = buffer_1000_ddpd[1022];

    auto g_y_0_0_0_yy_xz_z_yy = buffer_1000_ddpd[1023];

    auto g_y_0_0_0_yy_xz_z_yz = buffer_1000_ddpd[1024];

    auto g_y_0_0_0_yy_xz_z_zz = buffer_1000_ddpd[1025];

    auto g_y_0_0_0_yy_yy_x_xx = buffer_1000_ddpd[1026];

    auto g_y_0_0_0_yy_yy_x_xy = buffer_1000_ddpd[1027];

    auto g_y_0_0_0_yy_yy_x_xz = buffer_1000_ddpd[1028];

    auto g_y_0_0_0_yy_yy_x_yy = buffer_1000_ddpd[1029];

    auto g_y_0_0_0_yy_yy_x_yz = buffer_1000_ddpd[1030];

    auto g_y_0_0_0_yy_yy_x_zz = buffer_1000_ddpd[1031];

    auto g_y_0_0_0_yy_yy_y_xx = buffer_1000_ddpd[1032];

    auto g_y_0_0_0_yy_yy_y_xy = buffer_1000_ddpd[1033];

    auto g_y_0_0_0_yy_yy_y_xz = buffer_1000_ddpd[1034];

    auto g_y_0_0_0_yy_yy_y_yy = buffer_1000_ddpd[1035];

    auto g_y_0_0_0_yy_yy_y_yz = buffer_1000_ddpd[1036];

    auto g_y_0_0_0_yy_yy_y_zz = buffer_1000_ddpd[1037];

    auto g_y_0_0_0_yy_yy_z_xx = buffer_1000_ddpd[1038];

    auto g_y_0_0_0_yy_yy_z_xy = buffer_1000_ddpd[1039];

    auto g_y_0_0_0_yy_yy_z_xz = buffer_1000_ddpd[1040];

    auto g_y_0_0_0_yy_yy_z_yy = buffer_1000_ddpd[1041];

    auto g_y_0_0_0_yy_yy_z_yz = buffer_1000_ddpd[1042];

    auto g_y_0_0_0_yy_yy_z_zz = buffer_1000_ddpd[1043];

    auto g_y_0_0_0_yy_yz_x_xx = buffer_1000_ddpd[1044];

    auto g_y_0_0_0_yy_yz_x_xy = buffer_1000_ddpd[1045];

    auto g_y_0_0_0_yy_yz_x_xz = buffer_1000_ddpd[1046];

    auto g_y_0_0_0_yy_yz_x_yy = buffer_1000_ddpd[1047];

    auto g_y_0_0_0_yy_yz_x_yz = buffer_1000_ddpd[1048];

    auto g_y_0_0_0_yy_yz_x_zz = buffer_1000_ddpd[1049];

    auto g_y_0_0_0_yy_yz_y_xx = buffer_1000_ddpd[1050];

    auto g_y_0_0_0_yy_yz_y_xy = buffer_1000_ddpd[1051];

    auto g_y_0_0_0_yy_yz_y_xz = buffer_1000_ddpd[1052];

    auto g_y_0_0_0_yy_yz_y_yy = buffer_1000_ddpd[1053];

    auto g_y_0_0_0_yy_yz_y_yz = buffer_1000_ddpd[1054];

    auto g_y_0_0_0_yy_yz_y_zz = buffer_1000_ddpd[1055];

    auto g_y_0_0_0_yy_yz_z_xx = buffer_1000_ddpd[1056];

    auto g_y_0_0_0_yy_yz_z_xy = buffer_1000_ddpd[1057];

    auto g_y_0_0_0_yy_yz_z_xz = buffer_1000_ddpd[1058];

    auto g_y_0_0_0_yy_yz_z_yy = buffer_1000_ddpd[1059];

    auto g_y_0_0_0_yy_yz_z_yz = buffer_1000_ddpd[1060];

    auto g_y_0_0_0_yy_yz_z_zz = buffer_1000_ddpd[1061];

    auto g_y_0_0_0_yy_zz_x_xx = buffer_1000_ddpd[1062];

    auto g_y_0_0_0_yy_zz_x_xy = buffer_1000_ddpd[1063];

    auto g_y_0_0_0_yy_zz_x_xz = buffer_1000_ddpd[1064];

    auto g_y_0_0_0_yy_zz_x_yy = buffer_1000_ddpd[1065];

    auto g_y_0_0_0_yy_zz_x_yz = buffer_1000_ddpd[1066];

    auto g_y_0_0_0_yy_zz_x_zz = buffer_1000_ddpd[1067];

    auto g_y_0_0_0_yy_zz_y_xx = buffer_1000_ddpd[1068];

    auto g_y_0_0_0_yy_zz_y_xy = buffer_1000_ddpd[1069];

    auto g_y_0_0_0_yy_zz_y_xz = buffer_1000_ddpd[1070];

    auto g_y_0_0_0_yy_zz_y_yy = buffer_1000_ddpd[1071];

    auto g_y_0_0_0_yy_zz_y_yz = buffer_1000_ddpd[1072];

    auto g_y_0_0_0_yy_zz_y_zz = buffer_1000_ddpd[1073];

    auto g_y_0_0_0_yy_zz_z_xx = buffer_1000_ddpd[1074];

    auto g_y_0_0_0_yy_zz_z_xy = buffer_1000_ddpd[1075];

    auto g_y_0_0_0_yy_zz_z_xz = buffer_1000_ddpd[1076];

    auto g_y_0_0_0_yy_zz_z_yy = buffer_1000_ddpd[1077];

    auto g_y_0_0_0_yy_zz_z_yz = buffer_1000_ddpd[1078];

    auto g_y_0_0_0_yy_zz_z_zz = buffer_1000_ddpd[1079];

    auto g_y_0_0_0_yz_xx_x_xx = buffer_1000_ddpd[1080];

    auto g_y_0_0_0_yz_xx_x_xy = buffer_1000_ddpd[1081];

    auto g_y_0_0_0_yz_xx_x_xz = buffer_1000_ddpd[1082];

    auto g_y_0_0_0_yz_xx_x_yy = buffer_1000_ddpd[1083];

    auto g_y_0_0_0_yz_xx_x_yz = buffer_1000_ddpd[1084];

    auto g_y_0_0_0_yz_xx_x_zz = buffer_1000_ddpd[1085];

    auto g_y_0_0_0_yz_xx_y_xx = buffer_1000_ddpd[1086];

    auto g_y_0_0_0_yz_xx_y_xy = buffer_1000_ddpd[1087];

    auto g_y_0_0_0_yz_xx_y_xz = buffer_1000_ddpd[1088];

    auto g_y_0_0_0_yz_xx_y_yy = buffer_1000_ddpd[1089];

    auto g_y_0_0_0_yz_xx_y_yz = buffer_1000_ddpd[1090];

    auto g_y_0_0_0_yz_xx_y_zz = buffer_1000_ddpd[1091];

    auto g_y_0_0_0_yz_xx_z_xx = buffer_1000_ddpd[1092];

    auto g_y_0_0_0_yz_xx_z_xy = buffer_1000_ddpd[1093];

    auto g_y_0_0_0_yz_xx_z_xz = buffer_1000_ddpd[1094];

    auto g_y_0_0_0_yz_xx_z_yy = buffer_1000_ddpd[1095];

    auto g_y_0_0_0_yz_xx_z_yz = buffer_1000_ddpd[1096];

    auto g_y_0_0_0_yz_xx_z_zz = buffer_1000_ddpd[1097];

    auto g_y_0_0_0_yz_xy_x_xx = buffer_1000_ddpd[1098];

    auto g_y_0_0_0_yz_xy_x_xy = buffer_1000_ddpd[1099];

    auto g_y_0_0_0_yz_xy_x_xz = buffer_1000_ddpd[1100];

    auto g_y_0_0_0_yz_xy_x_yy = buffer_1000_ddpd[1101];

    auto g_y_0_0_0_yz_xy_x_yz = buffer_1000_ddpd[1102];

    auto g_y_0_0_0_yz_xy_x_zz = buffer_1000_ddpd[1103];

    auto g_y_0_0_0_yz_xy_y_xx = buffer_1000_ddpd[1104];

    auto g_y_0_0_0_yz_xy_y_xy = buffer_1000_ddpd[1105];

    auto g_y_0_0_0_yz_xy_y_xz = buffer_1000_ddpd[1106];

    auto g_y_0_0_0_yz_xy_y_yy = buffer_1000_ddpd[1107];

    auto g_y_0_0_0_yz_xy_y_yz = buffer_1000_ddpd[1108];

    auto g_y_0_0_0_yz_xy_y_zz = buffer_1000_ddpd[1109];

    auto g_y_0_0_0_yz_xy_z_xx = buffer_1000_ddpd[1110];

    auto g_y_0_0_0_yz_xy_z_xy = buffer_1000_ddpd[1111];

    auto g_y_0_0_0_yz_xy_z_xz = buffer_1000_ddpd[1112];

    auto g_y_0_0_0_yz_xy_z_yy = buffer_1000_ddpd[1113];

    auto g_y_0_0_0_yz_xy_z_yz = buffer_1000_ddpd[1114];

    auto g_y_0_0_0_yz_xy_z_zz = buffer_1000_ddpd[1115];

    auto g_y_0_0_0_yz_xz_x_xx = buffer_1000_ddpd[1116];

    auto g_y_0_0_0_yz_xz_x_xy = buffer_1000_ddpd[1117];

    auto g_y_0_0_0_yz_xz_x_xz = buffer_1000_ddpd[1118];

    auto g_y_0_0_0_yz_xz_x_yy = buffer_1000_ddpd[1119];

    auto g_y_0_0_0_yz_xz_x_yz = buffer_1000_ddpd[1120];

    auto g_y_0_0_0_yz_xz_x_zz = buffer_1000_ddpd[1121];

    auto g_y_0_0_0_yz_xz_y_xx = buffer_1000_ddpd[1122];

    auto g_y_0_0_0_yz_xz_y_xy = buffer_1000_ddpd[1123];

    auto g_y_0_0_0_yz_xz_y_xz = buffer_1000_ddpd[1124];

    auto g_y_0_0_0_yz_xz_y_yy = buffer_1000_ddpd[1125];

    auto g_y_0_0_0_yz_xz_y_yz = buffer_1000_ddpd[1126];

    auto g_y_0_0_0_yz_xz_y_zz = buffer_1000_ddpd[1127];

    auto g_y_0_0_0_yz_xz_z_xx = buffer_1000_ddpd[1128];

    auto g_y_0_0_0_yz_xz_z_xy = buffer_1000_ddpd[1129];

    auto g_y_0_0_0_yz_xz_z_xz = buffer_1000_ddpd[1130];

    auto g_y_0_0_0_yz_xz_z_yy = buffer_1000_ddpd[1131];

    auto g_y_0_0_0_yz_xz_z_yz = buffer_1000_ddpd[1132];

    auto g_y_0_0_0_yz_xz_z_zz = buffer_1000_ddpd[1133];

    auto g_y_0_0_0_yz_yy_x_xx = buffer_1000_ddpd[1134];

    auto g_y_0_0_0_yz_yy_x_xy = buffer_1000_ddpd[1135];

    auto g_y_0_0_0_yz_yy_x_xz = buffer_1000_ddpd[1136];

    auto g_y_0_0_0_yz_yy_x_yy = buffer_1000_ddpd[1137];

    auto g_y_0_0_0_yz_yy_x_yz = buffer_1000_ddpd[1138];

    auto g_y_0_0_0_yz_yy_x_zz = buffer_1000_ddpd[1139];

    auto g_y_0_0_0_yz_yy_y_xx = buffer_1000_ddpd[1140];

    auto g_y_0_0_0_yz_yy_y_xy = buffer_1000_ddpd[1141];

    auto g_y_0_0_0_yz_yy_y_xz = buffer_1000_ddpd[1142];

    auto g_y_0_0_0_yz_yy_y_yy = buffer_1000_ddpd[1143];

    auto g_y_0_0_0_yz_yy_y_yz = buffer_1000_ddpd[1144];

    auto g_y_0_0_0_yz_yy_y_zz = buffer_1000_ddpd[1145];

    auto g_y_0_0_0_yz_yy_z_xx = buffer_1000_ddpd[1146];

    auto g_y_0_0_0_yz_yy_z_xy = buffer_1000_ddpd[1147];

    auto g_y_0_0_0_yz_yy_z_xz = buffer_1000_ddpd[1148];

    auto g_y_0_0_0_yz_yy_z_yy = buffer_1000_ddpd[1149];

    auto g_y_0_0_0_yz_yy_z_yz = buffer_1000_ddpd[1150];

    auto g_y_0_0_0_yz_yy_z_zz = buffer_1000_ddpd[1151];

    auto g_y_0_0_0_yz_yz_x_xx = buffer_1000_ddpd[1152];

    auto g_y_0_0_0_yz_yz_x_xy = buffer_1000_ddpd[1153];

    auto g_y_0_0_0_yz_yz_x_xz = buffer_1000_ddpd[1154];

    auto g_y_0_0_0_yz_yz_x_yy = buffer_1000_ddpd[1155];

    auto g_y_0_0_0_yz_yz_x_yz = buffer_1000_ddpd[1156];

    auto g_y_0_0_0_yz_yz_x_zz = buffer_1000_ddpd[1157];

    auto g_y_0_0_0_yz_yz_y_xx = buffer_1000_ddpd[1158];

    auto g_y_0_0_0_yz_yz_y_xy = buffer_1000_ddpd[1159];

    auto g_y_0_0_0_yz_yz_y_xz = buffer_1000_ddpd[1160];

    auto g_y_0_0_0_yz_yz_y_yy = buffer_1000_ddpd[1161];

    auto g_y_0_0_0_yz_yz_y_yz = buffer_1000_ddpd[1162];

    auto g_y_0_0_0_yz_yz_y_zz = buffer_1000_ddpd[1163];

    auto g_y_0_0_0_yz_yz_z_xx = buffer_1000_ddpd[1164];

    auto g_y_0_0_0_yz_yz_z_xy = buffer_1000_ddpd[1165];

    auto g_y_0_0_0_yz_yz_z_xz = buffer_1000_ddpd[1166];

    auto g_y_0_0_0_yz_yz_z_yy = buffer_1000_ddpd[1167];

    auto g_y_0_0_0_yz_yz_z_yz = buffer_1000_ddpd[1168];

    auto g_y_0_0_0_yz_yz_z_zz = buffer_1000_ddpd[1169];

    auto g_y_0_0_0_yz_zz_x_xx = buffer_1000_ddpd[1170];

    auto g_y_0_0_0_yz_zz_x_xy = buffer_1000_ddpd[1171];

    auto g_y_0_0_0_yz_zz_x_xz = buffer_1000_ddpd[1172];

    auto g_y_0_0_0_yz_zz_x_yy = buffer_1000_ddpd[1173];

    auto g_y_0_0_0_yz_zz_x_yz = buffer_1000_ddpd[1174];

    auto g_y_0_0_0_yz_zz_x_zz = buffer_1000_ddpd[1175];

    auto g_y_0_0_0_yz_zz_y_xx = buffer_1000_ddpd[1176];

    auto g_y_0_0_0_yz_zz_y_xy = buffer_1000_ddpd[1177];

    auto g_y_0_0_0_yz_zz_y_xz = buffer_1000_ddpd[1178];

    auto g_y_0_0_0_yz_zz_y_yy = buffer_1000_ddpd[1179];

    auto g_y_0_0_0_yz_zz_y_yz = buffer_1000_ddpd[1180];

    auto g_y_0_0_0_yz_zz_y_zz = buffer_1000_ddpd[1181];

    auto g_y_0_0_0_yz_zz_z_xx = buffer_1000_ddpd[1182];

    auto g_y_0_0_0_yz_zz_z_xy = buffer_1000_ddpd[1183];

    auto g_y_0_0_0_yz_zz_z_xz = buffer_1000_ddpd[1184];

    auto g_y_0_0_0_yz_zz_z_yy = buffer_1000_ddpd[1185];

    auto g_y_0_0_0_yz_zz_z_yz = buffer_1000_ddpd[1186];

    auto g_y_0_0_0_yz_zz_z_zz = buffer_1000_ddpd[1187];

    auto g_y_0_0_0_zz_xx_x_xx = buffer_1000_ddpd[1188];

    auto g_y_0_0_0_zz_xx_x_xy = buffer_1000_ddpd[1189];

    auto g_y_0_0_0_zz_xx_x_xz = buffer_1000_ddpd[1190];

    auto g_y_0_0_0_zz_xx_x_yy = buffer_1000_ddpd[1191];

    auto g_y_0_0_0_zz_xx_x_yz = buffer_1000_ddpd[1192];

    auto g_y_0_0_0_zz_xx_x_zz = buffer_1000_ddpd[1193];

    auto g_y_0_0_0_zz_xx_y_xx = buffer_1000_ddpd[1194];

    auto g_y_0_0_0_zz_xx_y_xy = buffer_1000_ddpd[1195];

    auto g_y_0_0_0_zz_xx_y_xz = buffer_1000_ddpd[1196];

    auto g_y_0_0_0_zz_xx_y_yy = buffer_1000_ddpd[1197];

    auto g_y_0_0_0_zz_xx_y_yz = buffer_1000_ddpd[1198];

    auto g_y_0_0_0_zz_xx_y_zz = buffer_1000_ddpd[1199];

    auto g_y_0_0_0_zz_xx_z_xx = buffer_1000_ddpd[1200];

    auto g_y_0_0_0_zz_xx_z_xy = buffer_1000_ddpd[1201];

    auto g_y_0_0_0_zz_xx_z_xz = buffer_1000_ddpd[1202];

    auto g_y_0_0_0_zz_xx_z_yy = buffer_1000_ddpd[1203];

    auto g_y_0_0_0_zz_xx_z_yz = buffer_1000_ddpd[1204];

    auto g_y_0_0_0_zz_xx_z_zz = buffer_1000_ddpd[1205];

    auto g_y_0_0_0_zz_xy_x_xx = buffer_1000_ddpd[1206];

    auto g_y_0_0_0_zz_xy_x_xy = buffer_1000_ddpd[1207];

    auto g_y_0_0_0_zz_xy_x_xz = buffer_1000_ddpd[1208];

    auto g_y_0_0_0_zz_xy_x_yy = buffer_1000_ddpd[1209];

    auto g_y_0_0_0_zz_xy_x_yz = buffer_1000_ddpd[1210];

    auto g_y_0_0_0_zz_xy_x_zz = buffer_1000_ddpd[1211];

    auto g_y_0_0_0_zz_xy_y_xx = buffer_1000_ddpd[1212];

    auto g_y_0_0_0_zz_xy_y_xy = buffer_1000_ddpd[1213];

    auto g_y_0_0_0_zz_xy_y_xz = buffer_1000_ddpd[1214];

    auto g_y_0_0_0_zz_xy_y_yy = buffer_1000_ddpd[1215];

    auto g_y_0_0_0_zz_xy_y_yz = buffer_1000_ddpd[1216];

    auto g_y_0_0_0_zz_xy_y_zz = buffer_1000_ddpd[1217];

    auto g_y_0_0_0_zz_xy_z_xx = buffer_1000_ddpd[1218];

    auto g_y_0_0_0_zz_xy_z_xy = buffer_1000_ddpd[1219];

    auto g_y_0_0_0_zz_xy_z_xz = buffer_1000_ddpd[1220];

    auto g_y_0_0_0_zz_xy_z_yy = buffer_1000_ddpd[1221];

    auto g_y_0_0_0_zz_xy_z_yz = buffer_1000_ddpd[1222];

    auto g_y_0_0_0_zz_xy_z_zz = buffer_1000_ddpd[1223];

    auto g_y_0_0_0_zz_xz_x_xx = buffer_1000_ddpd[1224];

    auto g_y_0_0_0_zz_xz_x_xy = buffer_1000_ddpd[1225];

    auto g_y_0_0_0_zz_xz_x_xz = buffer_1000_ddpd[1226];

    auto g_y_0_0_0_zz_xz_x_yy = buffer_1000_ddpd[1227];

    auto g_y_0_0_0_zz_xz_x_yz = buffer_1000_ddpd[1228];

    auto g_y_0_0_0_zz_xz_x_zz = buffer_1000_ddpd[1229];

    auto g_y_0_0_0_zz_xz_y_xx = buffer_1000_ddpd[1230];

    auto g_y_0_0_0_zz_xz_y_xy = buffer_1000_ddpd[1231];

    auto g_y_0_0_0_zz_xz_y_xz = buffer_1000_ddpd[1232];

    auto g_y_0_0_0_zz_xz_y_yy = buffer_1000_ddpd[1233];

    auto g_y_0_0_0_zz_xz_y_yz = buffer_1000_ddpd[1234];

    auto g_y_0_0_0_zz_xz_y_zz = buffer_1000_ddpd[1235];

    auto g_y_0_0_0_zz_xz_z_xx = buffer_1000_ddpd[1236];

    auto g_y_0_0_0_zz_xz_z_xy = buffer_1000_ddpd[1237];

    auto g_y_0_0_0_zz_xz_z_xz = buffer_1000_ddpd[1238];

    auto g_y_0_0_0_zz_xz_z_yy = buffer_1000_ddpd[1239];

    auto g_y_0_0_0_zz_xz_z_yz = buffer_1000_ddpd[1240];

    auto g_y_0_0_0_zz_xz_z_zz = buffer_1000_ddpd[1241];

    auto g_y_0_0_0_zz_yy_x_xx = buffer_1000_ddpd[1242];

    auto g_y_0_0_0_zz_yy_x_xy = buffer_1000_ddpd[1243];

    auto g_y_0_0_0_zz_yy_x_xz = buffer_1000_ddpd[1244];

    auto g_y_0_0_0_zz_yy_x_yy = buffer_1000_ddpd[1245];

    auto g_y_0_0_0_zz_yy_x_yz = buffer_1000_ddpd[1246];

    auto g_y_0_0_0_zz_yy_x_zz = buffer_1000_ddpd[1247];

    auto g_y_0_0_0_zz_yy_y_xx = buffer_1000_ddpd[1248];

    auto g_y_0_0_0_zz_yy_y_xy = buffer_1000_ddpd[1249];

    auto g_y_0_0_0_zz_yy_y_xz = buffer_1000_ddpd[1250];

    auto g_y_0_0_0_zz_yy_y_yy = buffer_1000_ddpd[1251];

    auto g_y_0_0_0_zz_yy_y_yz = buffer_1000_ddpd[1252];

    auto g_y_0_0_0_zz_yy_y_zz = buffer_1000_ddpd[1253];

    auto g_y_0_0_0_zz_yy_z_xx = buffer_1000_ddpd[1254];

    auto g_y_0_0_0_zz_yy_z_xy = buffer_1000_ddpd[1255];

    auto g_y_0_0_0_zz_yy_z_xz = buffer_1000_ddpd[1256];

    auto g_y_0_0_0_zz_yy_z_yy = buffer_1000_ddpd[1257];

    auto g_y_0_0_0_zz_yy_z_yz = buffer_1000_ddpd[1258];

    auto g_y_0_0_0_zz_yy_z_zz = buffer_1000_ddpd[1259];

    auto g_y_0_0_0_zz_yz_x_xx = buffer_1000_ddpd[1260];

    auto g_y_0_0_0_zz_yz_x_xy = buffer_1000_ddpd[1261];

    auto g_y_0_0_0_zz_yz_x_xz = buffer_1000_ddpd[1262];

    auto g_y_0_0_0_zz_yz_x_yy = buffer_1000_ddpd[1263];

    auto g_y_0_0_0_zz_yz_x_yz = buffer_1000_ddpd[1264];

    auto g_y_0_0_0_zz_yz_x_zz = buffer_1000_ddpd[1265];

    auto g_y_0_0_0_zz_yz_y_xx = buffer_1000_ddpd[1266];

    auto g_y_0_0_0_zz_yz_y_xy = buffer_1000_ddpd[1267];

    auto g_y_0_0_0_zz_yz_y_xz = buffer_1000_ddpd[1268];

    auto g_y_0_0_0_zz_yz_y_yy = buffer_1000_ddpd[1269];

    auto g_y_0_0_0_zz_yz_y_yz = buffer_1000_ddpd[1270];

    auto g_y_0_0_0_zz_yz_y_zz = buffer_1000_ddpd[1271];

    auto g_y_0_0_0_zz_yz_z_xx = buffer_1000_ddpd[1272];

    auto g_y_0_0_0_zz_yz_z_xy = buffer_1000_ddpd[1273];

    auto g_y_0_0_0_zz_yz_z_xz = buffer_1000_ddpd[1274];

    auto g_y_0_0_0_zz_yz_z_yy = buffer_1000_ddpd[1275];

    auto g_y_0_0_0_zz_yz_z_yz = buffer_1000_ddpd[1276];

    auto g_y_0_0_0_zz_yz_z_zz = buffer_1000_ddpd[1277];

    auto g_y_0_0_0_zz_zz_x_xx = buffer_1000_ddpd[1278];

    auto g_y_0_0_0_zz_zz_x_xy = buffer_1000_ddpd[1279];

    auto g_y_0_0_0_zz_zz_x_xz = buffer_1000_ddpd[1280];

    auto g_y_0_0_0_zz_zz_x_yy = buffer_1000_ddpd[1281];

    auto g_y_0_0_0_zz_zz_x_yz = buffer_1000_ddpd[1282];

    auto g_y_0_0_0_zz_zz_x_zz = buffer_1000_ddpd[1283];

    auto g_y_0_0_0_zz_zz_y_xx = buffer_1000_ddpd[1284];

    auto g_y_0_0_0_zz_zz_y_xy = buffer_1000_ddpd[1285];

    auto g_y_0_0_0_zz_zz_y_xz = buffer_1000_ddpd[1286];

    auto g_y_0_0_0_zz_zz_y_yy = buffer_1000_ddpd[1287];

    auto g_y_0_0_0_zz_zz_y_yz = buffer_1000_ddpd[1288];

    auto g_y_0_0_0_zz_zz_y_zz = buffer_1000_ddpd[1289];

    auto g_y_0_0_0_zz_zz_z_xx = buffer_1000_ddpd[1290];

    auto g_y_0_0_0_zz_zz_z_xy = buffer_1000_ddpd[1291];

    auto g_y_0_0_0_zz_zz_z_xz = buffer_1000_ddpd[1292];

    auto g_y_0_0_0_zz_zz_z_yy = buffer_1000_ddpd[1293];

    auto g_y_0_0_0_zz_zz_z_yz = buffer_1000_ddpd[1294];

    auto g_y_0_0_0_zz_zz_z_zz = buffer_1000_ddpd[1295];

    auto g_z_0_0_0_xx_xx_x_xx = buffer_1000_ddpd[1296];

    auto g_z_0_0_0_xx_xx_x_xy = buffer_1000_ddpd[1297];

    auto g_z_0_0_0_xx_xx_x_xz = buffer_1000_ddpd[1298];

    auto g_z_0_0_0_xx_xx_x_yy = buffer_1000_ddpd[1299];

    auto g_z_0_0_0_xx_xx_x_yz = buffer_1000_ddpd[1300];

    auto g_z_0_0_0_xx_xx_x_zz = buffer_1000_ddpd[1301];

    auto g_z_0_0_0_xx_xx_y_xx = buffer_1000_ddpd[1302];

    auto g_z_0_0_0_xx_xx_y_xy = buffer_1000_ddpd[1303];

    auto g_z_0_0_0_xx_xx_y_xz = buffer_1000_ddpd[1304];

    auto g_z_0_0_0_xx_xx_y_yy = buffer_1000_ddpd[1305];

    auto g_z_0_0_0_xx_xx_y_yz = buffer_1000_ddpd[1306];

    auto g_z_0_0_0_xx_xx_y_zz = buffer_1000_ddpd[1307];

    auto g_z_0_0_0_xx_xx_z_xx = buffer_1000_ddpd[1308];

    auto g_z_0_0_0_xx_xx_z_xy = buffer_1000_ddpd[1309];

    auto g_z_0_0_0_xx_xx_z_xz = buffer_1000_ddpd[1310];

    auto g_z_0_0_0_xx_xx_z_yy = buffer_1000_ddpd[1311];

    auto g_z_0_0_0_xx_xx_z_yz = buffer_1000_ddpd[1312];

    auto g_z_0_0_0_xx_xx_z_zz = buffer_1000_ddpd[1313];

    auto g_z_0_0_0_xx_xy_x_xx = buffer_1000_ddpd[1314];

    auto g_z_0_0_0_xx_xy_x_xy = buffer_1000_ddpd[1315];

    auto g_z_0_0_0_xx_xy_x_xz = buffer_1000_ddpd[1316];

    auto g_z_0_0_0_xx_xy_x_yy = buffer_1000_ddpd[1317];

    auto g_z_0_0_0_xx_xy_x_yz = buffer_1000_ddpd[1318];

    auto g_z_0_0_0_xx_xy_x_zz = buffer_1000_ddpd[1319];

    auto g_z_0_0_0_xx_xy_y_xx = buffer_1000_ddpd[1320];

    auto g_z_0_0_0_xx_xy_y_xy = buffer_1000_ddpd[1321];

    auto g_z_0_0_0_xx_xy_y_xz = buffer_1000_ddpd[1322];

    auto g_z_0_0_0_xx_xy_y_yy = buffer_1000_ddpd[1323];

    auto g_z_0_0_0_xx_xy_y_yz = buffer_1000_ddpd[1324];

    auto g_z_0_0_0_xx_xy_y_zz = buffer_1000_ddpd[1325];

    auto g_z_0_0_0_xx_xy_z_xx = buffer_1000_ddpd[1326];

    auto g_z_0_0_0_xx_xy_z_xy = buffer_1000_ddpd[1327];

    auto g_z_0_0_0_xx_xy_z_xz = buffer_1000_ddpd[1328];

    auto g_z_0_0_0_xx_xy_z_yy = buffer_1000_ddpd[1329];

    auto g_z_0_0_0_xx_xy_z_yz = buffer_1000_ddpd[1330];

    auto g_z_0_0_0_xx_xy_z_zz = buffer_1000_ddpd[1331];

    auto g_z_0_0_0_xx_xz_x_xx = buffer_1000_ddpd[1332];

    auto g_z_0_0_0_xx_xz_x_xy = buffer_1000_ddpd[1333];

    auto g_z_0_0_0_xx_xz_x_xz = buffer_1000_ddpd[1334];

    auto g_z_0_0_0_xx_xz_x_yy = buffer_1000_ddpd[1335];

    auto g_z_0_0_0_xx_xz_x_yz = buffer_1000_ddpd[1336];

    auto g_z_0_0_0_xx_xz_x_zz = buffer_1000_ddpd[1337];

    auto g_z_0_0_0_xx_xz_y_xx = buffer_1000_ddpd[1338];

    auto g_z_0_0_0_xx_xz_y_xy = buffer_1000_ddpd[1339];

    auto g_z_0_0_0_xx_xz_y_xz = buffer_1000_ddpd[1340];

    auto g_z_0_0_0_xx_xz_y_yy = buffer_1000_ddpd[1341];

    auto g_z_0_0_0_xx_xz_y_yz = buffer_1000_ddpd[1342];

    auto g_z_0_0_0_xx_xz_y_zz = buffer_1000_ddpd[1343];

    auto g_z_0_0_0_xx_xz_z_xx = buffer_1000_ddpd[1344];

    auto g_z_0_0_0_xx_xz_z_xy = buffer_1000_ddpd[1345];

    auto g_z_0_0_0_xx_xz_z_xz = buffer_1000_ddpd[1346];

    auto g_z_0_0_0_xx_xz_z_yy = buffer_1000_ddpd[1347];

    auto g_z_0_0_0_xx_xz_z_yz = buffer_1000_ddpd[1348];

    auto g_z_0_0_0_xx_xz_z_zz = buffer_1000_ddpd[1349];

    auto g_z_0_0_0_xx_yy_x_xx = buffer_1000_ddpd[1350];

    auto g_z_0_0_0_xx_yy_x_xy = buffer_1000_ddpd[1351];

    auto g_z_0_0_0_xx_yy_x_xz = buffer_1000_ddpd[1352];

    auto g_z_0_0_0_xx_yy_x_yy = buffer_1000_ddpd[1353];

    auto g_z_0_0_0_xx_yy_x_yz = buffer_1000_ddpd[1354];

    auto g_z_0_0_0_xx_yy_x_zz = buffer_1000_ddpd[1355];

    auto g_z_0_0_0_xx_yy_y_xx = buffer_1000_ddpd[1356];

    auto g_z_0_0_0_xx_yy_y_xy = buffer_1000_ddpd[1357];

    auto g_z_0_0_0_xx_yy_y_xz = buffer_1000_ddpd[1358];

    auto g_z_0_0_0_xx_yy_y_yy = buffer_1000_ddpd[1359];

    auto g_z_0_0_0_xx_yy_y_yz = buffer_1000_ddpd[1360];

    auto g_z_0_0_0_xx_yy_y_zz = buffer_1000_ddpd[1361];

    auto g_z_0_0_0_xx_yy_z_xx = buffer_1000_ddpd[1362];

    auto g_z_0_0_0_xx_yy_z_xy = buffer_1000_ddpd[1363];

    auto g_z_0_0_0_xx_yy_z_xz = buffer_1000_ddpd[1364];

    auto g_z_0_0_0_xx_yy_z_yy = buffer_1000_ddpd[1365];

    auto g_z_0_0_0_xx_yy_z_yz = buffer_1000_ddpd[1366];

    auto g_z_0_0_0_xx_yy_z_zz = buffer_1000_ddpd[1367];

    auto g_z_0_0_0_xx_yz_x_xx = buffer_1000_ddpd[1368];

    auto g_z_0_0_0_xx_yz_x_xy = buffer_1000_ddpd[1369];

    auto g_z_0_0_0_xx_yz_x_xz = buffer_1000_ddpd[1370];

    auto g_z_0_0_0_xx_yz_x_yy = buffer_1000_ddpd[1371];

    auto g_z_0_0_0_xx_yz_x_yz = buffer_1000_ddpd[1372];

    auto g_z_0_0_0_xx_yz_x_zz = buffer_1000_ddpd[1373];

    auto g_z_0_0_0_xx_yz_y_xx = buffer_1000_ddpd[1374];

    auto g_z_0_0_0_xx_yz_y_xy = buffer_1000_ddpd[1375];

    auto g_z_0_0_0_xx_yz_y_xz = buffer_1000_ddpd[1376];

    auto g_z_0_0_0_xx_yz_y_yy = buffer_1000_ddpd[1377];

    auto g_z_0_0_0_xx_yz_y_yz = buffer_1000_ddpd[1378];

    auto g_z_0_0_0_xx_yz_y_zz = buffer_1000_ddpd[1379];

    auto g_z_0_0_0_xx_yz_z_xx = buffer_1000_ddpd[1380];

    auto g_z_0_0_0_xx_yz_z_xy = buffer_1000_ddpd[1381];

    auto g_z_0_0_0_xx_yz_z_xz = buffer_1000_ddpd[1382];

    auto g_z_0_0_0_xx_yz_z_yy = buffer_1000_ddpd[1383];

    auto g_z_0_0_0_xx_yz_z_yz = buffer_1000_ddpd[1384];

    auto g_z_0_0_0_xx_yz_z_zz = buffer_1000_ddpd[1385];

    auto g_z_0_0_0_xx_zz_x_xx = buffer_1000_ddpd[1386];

    auto g_z_0_0_0_xx_zz_x_xy = buffer_1000_ddpd[1387];

    auto g_z_0_0_0_xx_zz_x_xz = buffer_1000_ddpd[1388];

    auto g_z_0_0_0_xx_zz_x_yy = buffer_1000_ddpd[1389];

    auto g_z_0_0_0_xx_zz_x_yz = buffer_1000_ddpd[1390];

    auto g_z_0_0_0_xx_zz_x_zz = buffer_1000_ddpd[1391];

    auto g_z_0_0_0_xx_zz_y_xx = buffer_1000_ddpd[1392];

    auto g_z_0_0_0_xx_zz_y_xy = buffer_1000_ddpd[1393];

    auto g_z_0_0_0_xx_zz_y_xz = buffer_1000_ddpd[1394];

    auto g_z_0_0_0_xx_zz_y_yy = buffer_1000_ddpd[1395];

    auto g_z_0_0_0_xx_zz_y_yz = buffer_1000_ddpd[1396];

    auto g_z_0_0_0_xx_zz_y_zz = buffer_1000_ddpd[1397];

    auto g_z_0_0_0_xx_zz_z_xx = buffer_1000_ddpd[1398];

    auto g_z_0_0_0_xx_zz_z_xy = buffer_1000_ddpd[1399];

    auto g_z_0_0_0_xx_zz_z_xz = buffer_1000_ddpd[1400];

    auto g_z_0_0_0_xx_zz_z_yy = buffer_1000_ddpd[1401];

    auto g_z_0_0_0_xx_zz_z_yz = buffer_1000_ddpd[1402];

    auto g_z_0_0_0_xx_zz_z_zz = buffer_1000_ddpd[1403];

    auto g_z_0_0_0_xy_xx_x_xx = buffer_1000_ddpd[1404];

    auto g_z_0_0_0_xy_xx_x_xy = buffer_1000_ddpd[1405];

    auto g_z_0_0_0_xy_xx_x_xz = buffer_1000_ddpd[1406];

    auto g_z_0_0_0_xy_xx_x_yy = buffer_1000_ddpd[1407];

    auto g_z_0_0_0_xy_xx_x_yz = buffer_1000_ddpd[1408];

    auto g_z_0_0_0_xy_xx_x_zz = buffer_1000_ddpd[1409];

    auto g_z_0_0_0_xy_xx_y_xx = buffer_1000_ddpd[1410];

    auto g_z_0_0_0_xy_xx_y_xy = buffer_1000_ddpd[1411];

    auto g_z_0_0_0_xy_xx_y_xz = buffer_1000_ddpd[1412];

    auto g_z_0_0_0_xy_xx_y_yy = buffer_1000_ddpd[1413];

    auto g_z_0_0_0_xy_xx_y_yz = buffer_1000_ddpd[1414];

    auto g_z_0_0_0_xy_xx_y_zz = buffer_1000_ddpd[1415];

    auto g_z_0_0_0_xy_xx_z_xx = buffer_1000_ddpd[1416];

    auto g_z_0_0_0_xy_xx_z_xy = buffer_1000_ddpd[1417];

    auto g_z_0_0_0_xy_xx_z_xz = buffer_1000_ddpd[1418];

    auto g_z_0_0_0_xy_xx_z_yy = buffer_1000_ddpd[1419];

    auto g_z_0_0_0_xy_xx_z_yz = buffer_1000_ddpd[1420];

    auto g_z_0_0_0_xy_xx_z_zz = buffer_1000_ddpd[1421];

    auto g_z_0_0_0_xy_xy_x_xx = buffer_1000_ddpd[1422];

    auto g_z_0_0_0_xy_xy_x_xy = buffer_1000_ddpd[1423];

    auto g_z_0_0_0_xy_xy_x_xz = buffer_1000_ddpd[1424];

    auto g_z_0_0_0_xy_xy_x_yy = buffer_1000_ddpd[1425];

    auto g_z_0_0_0_xy_xy_x_yz = buffer_1000_ddpd[1426];

    auto g_z_0_0_0_xy_xy_x_zz = buffer_1000_ddpd[1427];

    auto g_z_0_0_0_xy_xy_y_xx = buffer_1000_ddpd[1428];

    auto g_z_0_0_0_xy_xy_y_xy = buffer_1000_ddpd[1429];

    auto g_z_0_0_0_xy_xy_y_xz = buffer_1000_ddpd[1430];

    auto g_z_0_0_0_xy_xy_y_yy = buffer_1000_ddpd[1431];

    auto g_z_0_0_0_xy_xy_y_yz = buffer_1000_ddpd[1432];

    auto g_z_0_0_0_xy_xy_y_zz = buffer_1000_ddpd[1433];

    auto g_z_0_0_0_xy_xy_z_xx = buffer_1000_ddpd[1434];

    auto g_z_0_0_0_xy_xy_z_xy = buffer_1000_ddpd[1435];

    auto g_z_0_0_0_xy_xy_z_xz = buffer_1000_ddpd[1436];

    auto g_z_0_0_0_xy_xy_z_yy = buffer_1000_ddpd[1437];

    auto g_z_0_0_0_xy_xy_z_yz = buffer_1000_ddpd[1438];

    auto g_z_0_0_0_xy_xy_z_zz = buffer_1000_ddpd[1439];

    auto g_z_0_0_0_xy_xz_x_xx = buffer_1000_ddpd[1440];

    auto g_z_0_0_0_xy_xz_x_xy = buffer_1000_ddpd[1441];

    auto g_z_0_0_0_xy_xz_x_xz = buffer_1000_ddpd[1442];

    auto g_z_0_0_0_xy_xz_x_yy = buffer_1000_ddpd[1443];

    auto g_z_0_0_0_xy_xz_x_yz = buffer_1000_ddpd[1444];

    auto g_z_0_0_0_xy_xz_x_zz = buffer_1000_ddpd[1445];

    auto g_z_0_0_0_xy_xz_y_xx = buffer_1000_ddpd[1446];

    auto g_z_0_0_0_xy_xz_y_xy = buffer_1000_ddpd[1447];

    auto g_z_0_0_0_xy_xz_y_xz = buffer_1000_ddpd[1448];

    auto g_z_0_0_0_xy_xz_y_yy = buffer_1000_ddpd[1449];

    auto g_z_0_0_0_xy_xz_y_yz = buffer_1000_ddpd[1450];

    auto g_z_0_0_0_xy_xz_y_zz = buffer_1000_ddpd[1451];

    auto g_z_0_0_0_xy_xz_z_xx = buffer_1000_ddpd[1452];

    auto g_z_0_0_0_xy_xz_z_xy = buffer_1000_ddpd[1453];

    auto g_z_0_0_0_xy_xz_z_xz = buffer_1000_ddpd[1454];

    auto g_z_0_0_0_xy_xz_z_yy = buffer_1000_ddpd[1455];

    auto g_z_0_0_0_xy_xz_z_yz = buffer_1000_ddpd[1456];

    auto g_z_0_0_0_xy_xz_z_zz = buffer_1000_ddpd[1457];

    auto g_z_0_0_0_xy_yy_x_xx = buffer_1000_ddpd[1458];

    auto g_z_0_0_0_xy_yy_x_xy = buffer_1000_ddpd[1459];

    auto g_z_0_0_0_xy_yy_x_xz = buffer_1000_ddpd[1460];

    auto g_z_0_0_0_xy_yy_x_yy = buffer_1000_ddpd[1461];

    auto g_z_0_0_0_xy_yy_x_yz = buffer_1000_ddpd[1462];

    auto g_z_0_0_0_xy_yy_x_zz = buffer_1000_ddpd[1463];

    auto g_z_0_0_0_xy_yy_y_xx = buffer_1000_ddpd[1464];

    auto g_z_0_0_0_xy_yy_y_xy = buffer_1000_ddpd[1465];

    auto g_z_0_0_0_xy_yy_y_xz = buffer_1000_ddpd[1466];

    auto g_z_0_0_0_xy_yy_y_yy = buffer_1000_ddpd[1467];

    auto g_z_0_0_0_xy_yy_y_yz = buffer_1000_ddpd[1468];

    auto g_z_0_0_0_xy_yy_y_zz = buffer_1000_ddpd[1469];

    auto g_z_0_0_0_xy_yy_z_xx = buffer_1000_ddpd[1470];

    auto g_z_0_0_0_xy_yy_z_xy = buffer_1000_ddpd[1471];

    auto g_z_0_0_0_xy_yy_z_xz = buffer_1000_ddpd[1472];

    auto g_z_0_0_0_xy_yy_z_yy = buffer_1000_ddpd[1473];

    auto g_z_0_0_0_xy_yy_z_yz = buffer_1000_ddpd[1474];

    auto g_z_0_0_0_xy_yy_z_zz = buffer_1000_ddpd[1475];

    auto g_z_0_0_0_xy_yz_x_xx = buffer_1000_ddpd[1476];

    auto g_z_0_0_0_xy_yz_x_xy = buffer_1000_ddpd[1477];

    auto g_z_0_0_0_xy_yz_x_xz = buffer_1000_ddpd[1478];

    auto g_z_0_0_0_xy_yz_x_yy = buffer_1000_ddpd[1479];

    auto g_z_0_0_0_xy_yz_x_yz = buffer_1000_ddpd[1480];

    auto g_z_0_0_0_xy_yz_x_zz = buffer_1000_ddpd[1481];

    auto g_z_0_0_0_xy_yz_y_xx = buffer_1000_ddpd[1482];

    auto g_z_0_0_0_xy_yz_y_xy = buffer_1000_ddpd[1483];

    auto g_z_0_0_0_xy_yz_y_xz = buffer_1000_ddpd[1484];

    auto g_z_0_0_0_xy_yz_y_yy = buffer_1000_ddpd[1485];

    auto g_z_0_0_0_xy_yz_y_yz = buffer_1000_ddpd[1486];

    auto g_z_0_0_0_xy_yz_y_zz = buffer_1000_ddpd[1487];

    auto g_z_0_0_0_xy_yz_z_xx = buffer_1000_ddpd[1488];

    auto g_z_0_0_0_xy_yz_z_xy = buffer_1000_ddpd[1489];

    auto g_z_0_0_0_xy_yz_z_xz = buffer_1000_ddpd[1490];

    auto g_z_0_0_0_xy_yz_z_yy = buffer_1000_ddpd[1491];

    auto g_z_0_0_0_xy_yz_z_yz = buffer_1000_ddpd[1492];

    auto g_z_0_0_0_xy_yz_z_zz = buffer_1000_ddpd[1493];

    auto g_z_0_0_0_xy_zz_x_xx = buffer_1000_ddpd[1494];

    auto g_z_0_0_0_xy_zz_x_xy = buffer_1000_ddpd[1495];

    auto g_z_0_0_0_xy_zz_x_xz = buffer_1000_ddpd[1496];

    auto g_z_0_0_0_xy_zz_x_yy = buffer_1000_ddpd[1497];

    auto g_z_0_0_0_xy_zz_x_yz = buffer_1000_ddpd[1498];

    auto g_z_0_0_0_xy_zz_x_zz = buffer_1000_ddpd[1499];

    auto g_z_0_0_0_xy_zz_y_xx = buffer_1000_ddpd[1500];

    auto g_z_0_0_0_xy_zz_y_xy = buffer_1000_ddpd[1501];

    auto g_z_0_0_0_xy_zz_y_xz = buffer_1000_ddpd[1502];

    auto g_z_0_0_0_xy_zz_y_yy = buffer_1000_ddpd[1503];

    auto g_z_0_0_0_xy_zz_y_yz = buffer_1000_ddpd[1504];

    auto g_z_0_0_0_xy_zz_y_zz = buffer_1000_ddpd[1505];

    auto g_z_0_0_0_xy_zz_z_xx = buffer_1000_ddpd[1506];

    auto g_z_0_0_0_xy_zz_z_xy = buffer_1000_ddpd[1507];

    auto g_z_0_0_0_xy_zz_z_xz = buffer_1000_ddpd[1508];

    auto g_z_0_0_0_xy_zz_z_yy = buffer_1000_ddpd[1509];

    auto g_z_0_0_0_xy_zz_z_yz = buffer_1000_ddpd[1510];

    auto g_z_0_0_0_xy_zz_z_zz = buffer_1000_ddpd[1511];

    auto g_z_0_0_0_xz_xx_x_xx = buffer_1000_ddpd[1512];

    auto g_z_0_0_0_xz_xx_x_xy = buffer_1000_ddpd[1513];

    auto g_z_0_0_0_xz_xx_x_xz = buffer_1000_ddpd[1514];

    auto g_z_0_0_0_xz_xx_x_yy = buffer_1000_ddpd[1515];

    auto g_z_0_0_0_xz_xx_x_yz = buffer_1000_ddpd[1516];

    auto g_z_0_0_0_xz_xx_x_zz = buffer_1000_ddpd[1517];

    auto g_z_0_0_0_xz_xx_y_xx = buffer_1000_ddpd[1518];

    auto g_z_0_0_0_xz_xx_y_xy = buffer_1000_ddpd[1519];

    auto g_z_0_0_0_xz_xx_y_xz = buffer_1000_ddpd[1520];

    auto g_z_0_0_0_xz_xx_y_yy = buffer_1000_ddpd[1521];

    auto g_z_0_0_0_xz_xx_y_yz = buffer_1000_ddpd[1522];

    auto g_z_0_0_0_xz_xx_y_zz = buffer_1000_ddpd[1523];

    auto g_z_0_0_0_xz_xx_z_xx = buffer_1000_ddpd[1524];

    auto g_z_0_0_0_xz_xx_z_xy = buffer_1000_ddpd[1525];

    auto g_z_0_0_0_xz_xx_z_xz = buffer_1000_ddpd[1526];

    auto g_z_0_0_0_xz_xx_z_yy = buffer_1000_ddpd[1527];

    auto g_z_0_0_0_xz_xx_z_yz = buffer_1000_ddpd[1528];

    auto g_z_0_0_0_xz_xx_z_zz = buffer_1000_ddpd[1529];

    auto g_z_0_0_0_xz_xy_x_xx = buffer_1000_ddpd[1530];

    auto g_z_0_0_0_xz_xy_x_xy = buffer_1000_ddpd[1531];

    auto g_z_0_0_0_xz_xy_x_xz = buffer_1000_ddpd[1532];

    auto g_z_0_0_0_xz_xy_x_yy = buffer_1000_ddpd[1533];

    auto g_z_0_0_0_xz_xy_x_yz = buffer_1000_ddpd[1534];

    auto g_z_0_0_0_xz_xy_x_zz = buffer_1000_ddpd[1535];

    auto g_z_0_0_0_xz_xy_y_xx = buffer_1000_ddpd[1536];

    auto g_z_0_0_0_xz_xy_y_xy = buffer_1000_ddpd[1537];

    auto g_z_0_0_0_xz_xy_y_xz = buffer_1000_ddpd[1538];

    auto g_z_0_0_0_xz_xy_y_yy = buffer_1000_ddpd[1539];

    auto g_z_0_0_0_xz_xy_y_yz = buffer_1000_ddpd[1540];

    auto g_z_0_0_0_xz_xy_y_zz = buffer_1000_ddpd[1541];

    auto g_z_0_0_0_xz_xy_z_xx = buffer_1000_ddpd[1542];

    auto g_z_0_0_0_xz_xy_z_xy = buffer_1000_ddpd[1543];

    auto g_z_0_0_0_xz_xy_z_xz = buffer_1000_ddpd[1544];

    auto g_z_0_0_0_xz_xy_z_yy = buffer_1000_ddpd[1545];

    auto g_z_0_0_0_xz_xy_z_yz = buffer_1000_ddpd[1546];

    auto g_z_0_0_0_xz_xy_z_zz = buffer_1000_ddpd[1547];

    auto g_z_0_0_0_xz_xz_x_xx = buffer_1000_ddpd[1548];

    auto g_z_0_0_0_xz_xz_x_xy = buffer_1000_ddpd[1549];

    auto g_z_0_0_0_xz_xz_x_xz = buffer_1000_ddpd[1550];

    auto g_z_0_0_0_xz_xz_x_yy = buffer_1000_ddpd[1551];

    auto g_z_0_0_0_xz_xz_x_yz = buffer_1000_ddpd[1552];

    auto g_z_0_0_0_xz_xz_x_zz = buffer_1000_ddpd[1553];

    auto g_z_0_0_0_xz_xz_y_xx = buffer_1000_ddpd[1554];

    auto g_z_0_0_0_xz_xz_y_xy = buffer_1000_ddpd[1555];

    auto g_z_0_0_0_xz_xz_y_xz = buffer_1000_ddpd[1556];

    auto g_z_0_0_0_xz_xz_y_yy = buffer_1000_ddpd[1557];

    auto g_z_0_0_0_xz_xz_y_yz = buffer_1000_ddpd[1558];

    auto g_z_0_0_0_xz_xz_y_zz = buffer_1000_ddpd[1559];

    auto g_z_0_0_0_xz_xz_z_xx = buffer_1000_ddpd[1560];

    auto g_z_0_0_0_xz_xz_z_xy = buffer_1000_ddpd[1561];

    auto g_z_0_0_0_xz_xz_z_xz = buffer_1000_ddpd[1562];

    auto g_z_0_0_0_xz_xz_z_yy = buffer_1000_ddpd[1563];

    auto g_z_0_0_0_xz_xz_z_yz = buffer_1000_ddpd[1564];

    auto g_z_0_0_0_xz_xz_z_zz = buffer_1000_ddpd[1565];

    auto g_z_0_0_0_xz_yy_x_xx = buffer_1000_ddpd[1566];

    auto g_z_0_0_0_xz_yy_x_xy = buffer_1000_ddpd[1567];

    auto g_z_0_0_0_xz_yy_x_xz = buffer_1000_ddpd[1568];

    auto g_z_0_0_0_xz_yy_x_yy = buffer_1000_ddpd[1569];

    auto g_z_0_0_0_xz_yy_x_yz = buffer_1000_ddpd[1570];

    auto g_z_0_0_0_xz_yy_x_zz = buffer_1000_ddpd[1571];

    auto g_z_0_0_0_xz_yy_y_xx = buffer_1000_ddpd[1572];

    auto g_z_0_0_0_xz_yy_y_xy = buffer_1000_ddpd[1573];

    auto g_z_0_0_0_xz_yy_y_xz = buffer_1000_ddpd[1574];

    auto g_z_0_0_0_xz_yy_y_yy = buffer_1000_ddpd[1575];

    auto g_z_0_0_0_xz_yy_y_yz = buffer_1000_ddpd[1576];

    auto g_z_0_0_0_xz_yy_y_zz = buffer_1000_ddpd[1577];

    auto g_z_0_0_0_xz_yy_z_xx = buffer_1000_ddpd[1578];

    auto g_z_0_0_0_xz_yy_z_xy = buffer_1000_ddpd[1579];

    auto g_z_0_0_0_xz_yy_z_xz = buffer_1000_ddpd[1580];

    auto g_z_0_0_0_xz_yy_z_yy = buffer_1000_ddpd[1581];

    auto g_z_0_0_0_xz_yy_z_yz = buffer_1000_ddpd[1582];

    auto g_z_0_0_0_xz_yy_z_zz = buffer_1000_ddpd[1583];

    auto g_z_0_0_0_xz_yz_x_xx = buffer_1000_ddpd[1584];

    auto g_z_0_0_0_xz_yz_x_xy = buffer_1000_ddpd[1585];

    auto g_z_0_0_0_xz_yz_x_xz = buffer_1000_ddpd[1586];

    auto g_z_0_0_0_xz_yz_x_yy = buffer_1000_ddpd[1587];

    auto g_z_0_0_0_xz_yz_x_yz = buffer_1000_ddpd[1588];

    auto g_z_0_0_0_xz_yz_x_zz = buffer_1000_ddpd[1589];

    auto g_z_0_0_0_xz_yz_y_xx = buffer_1000_ddpd[1590];

    auto g_z_0_0_0_xz_yz_y_xy = buffer_1000_ddpd[1591];

    auto g_z_0_0_0_xz_yz_y_xz = buffer_1000_ddpd[1592];

    auto g_z_0_0_0_xz_yz_y_yy = buffer_1000_ddpd[1593];

    auto g_z_0_0_0_xz_yz_y_yz = buffer_1000_ddpd[1594];

    auto g_z_0_0_0_xz_yz_y_zz = buffer_1000_ddpd[1595];

    auto g_z_0_0_0_xz_yz_z_xx = buffer_1000_ddpd[1596];

    auto g_z_0_0_0_xz_yz_z_xy = buffer_1000_ddpd[1597];

    auto g_z_0_0_0_xz_yz_z_xz = buffer_1000_ddpd[1598];

    auto g_z_0_0_0_xz_yz_z_yy = buffer_1000_ddpd[1599];

    auto g_z_0_0_0_xz_yz_z_yz = buffer_1000_ddpd[1600];

    auto g_z_0_0_0_xz_yz_z_zz = buffer_1000_ddpd[1601];

    auto g_z_0_0_0_xz_zz_x_xx = buffer_1000_ddpd[1602];

    auto g_z_0_0_0_xz_zz_x_xy = buffer_1000_ddpd[1603];

    auto g_z_0_0_0_xz_zz_x_xz = buffer_1000_ddpd[1604];

    auto g_z_0_0_0_xz_zz_x_yy = buffer_1000_ddpd[1605];

    auto g_z_0_0_0_xz_zz_x_yz = buffer_1000_ddpd[1606];

    auto g_z_0_0_0_xz_zz_x_zz = buffer_1000_ddpd[1607];

    auto g_z_0_0_0_xz_zz_y_xx = buffer_1000_ddpd[1608];

    auto g_z_0_0_0_xz_zz_y_xy = buffer_1000_ddpd[1609];

    auto g_z_0_0_0_xz_zz_y_xz = buffer_1000_ddpd[1610];

    auto g_z_0_0_0_xz_zz_y_yy = buffer_1000_ddpd[1611];

    auto g_z_0_0_0_xz_zz_y_yz = buffer_1000_ddpd[1612];

    auto g_z_0_0_0_xz_zz_y_zz = buffer_1000_ddpd[1613];

    auto g_z_0_0_0_xz_zz_z_xx = buffer_1000_ddpd[1614];

    auto g_z_0_0_0_xz_zz_z_xy = buffer_1000_ddpd[1615];

    auto g_z_0_0_0_xz_zz_z_xz = buffer_1000_ddpd[1616];

    auto g_z_0_0_0_xz_zz_z_yy = buffer_1000_ddpd[1617];

    auto g_z_0_0_0_xz_zz_z_yz = buffer_1000_ddpd[1618];

    auto g_z_0_0_0_xz_zz_z_zz = buffer_1000_ddpd[1619];

    auto g_z_0_0_0_yy_xx_x_xx = buffer_1000_ddpd[1620];

    auto g_z_0_0_0_yy_xx_x_xy = buffer_1000_ddpd[1621];

    auto g_z_0_0_0_yy_xx_x_xz = buffer_1000_ddpd[1622];

    auto g_z_0_0_0_yy_xx_x_yy = buffer_1000_ddpd[1623];

    auto g_z_0_0_0_yy_xx_x_yz = buffer_1000_ddpd[1624];

    auto g_z_0_0_0_yy_xx_x_zz = buffer_1000_ddpd[1625];

    auto g_z_0_0_0_yy_xx_y_xx = buffer_1000_ddpd[1626];

    auto g_z_0_0_0_yy_xx_y_xy = buffer_1000_ddpd[1627];

    auto g_z_0_0_0_yy_xx_y_xz = buffer_1000_ddpd[1628];

    auto g_z_0_0_0_yy_xx_y_yy = buffer_1000_ddpd[1629];

    auto g_z_0_0_0_yy_xx_y_yz = buffer_1000_ddpd[1630];

    auto g_z_0_0_0_yy_xx_y_zz = buffer_1000_ddpd[1631];

    auto g_z_0_0_0_yy_xx_z_xx = buffer_1000_ddpd[1632];

    auto g_z_0_0_0_yy_xx_z_xy = buffer_1000_ddpd[1633];

    auto g_z_0_0_0_yy_xx_z_xz = buffer_1000_ddpd[1634];

    auto g_z_0_0_0_yy_xx_z_yy = buffer_1000_ddpd[1635];

    auto g_z_0_0_0_yy_xx_z_yz = buffer_1000_ddpd[1636];

    auto g_z_0_0_0_yy_xx_z_zz = buffer_1000_ddpd[1637];

    auto g_z_0_0_0_yy_xy_x_xx = buffer_1000_ddpd[1638];

    auto g_z_0_0_0_yy_xy_x_xy = buffer_1000_ddpd[1639];

    auto g_z_0_0_0_yy_xy_x_xz = buffer_1000_ddpd[1640];

    auto g_z_0_0_0_yy_xy_x_yy = buffer_1000_ddpd[1641];

    auto g_z_0_0_0_yy_xy_x_yz = buffer_1000_ddpd[1642];

    auto g_z_0_0_0_yy_xy_x_zz = buffer_1000_ddpd[1643];

    auto g_z_0_0_0_yy_xy_y_xx = buffer_1000_ddpd[1644];

    auto g_z_0_0_0_yy_xy_y_xy = buffer_1000_ddpd[1645];

    auto g_z_0_0_0_yy_xy_y_xz = buffer_1000_ddpd[1646];

    auto g_z_0_0_0_yy_xy_y_yy = buffer_1000_ddpd[1647];

    auto g_z_0_0_0_yy_xy_y_yz = buffer_1000_ddpd[1648];

    auto g_z_0_0_0_yy_xy_y_zz = buffer_1000_ddpd[1649];

    auto g_z_0_0_0_yy_xy_z_xx = buffer_1000_ddpd[1650];

    auto g_z_0_0_0_yy_xy_z_xy = buffer_1000_ddpd[1651];

    auto g_z_0_0_0_yy_xy_z_xz = buffer_1000_ddpd[1652];

    auto g_z_0_0_0_yy_xy_z_yy = buffer_1000_ddpd[1653];

    auto g_z_0_0_0_yy_xy_z_yz = buffer_1000_ddpd[1654];

    auto g_z_0_0_0_yy_xy_z_zz = buffer_1000_ddpd[1655];

    auto g_z_0_0_0_yy_xz_x_xx = buffer_1000_ddpd[1656];

    auto g_z_0_0_0_yy_xz_x_xy = buffer_1000_ddpd[1657];

    auto g_z_0_0_0_yy_xz_x_xz = buffer_1000_ddpd[1658];

    auto g_z_0_0_0_yy_xz_x_yy = buffer_1000_ddpd[1659];

    auto g_z_0_0_0_yy_xz_x_yz = buffer_1000_ddpd[1660];

    auto g_z_0_0_0_yy_xz_x_zz = buffer_1000_ddpd[1661];

    auto g_z_0_0_0_yy_xz_y_xx = buffer_1000_ddpd[1662];

    auto g_z_0_0_0_yy_xz_y_xy = buffer_1000_ddpd[1663];

    auto g_z_0_0_0_yy_xz_y_xz = buffer_1000_ddpd[1664];

    auto g_z_0_0_0_yy_xz_y_yy = buffer_1000_ddpd[1665];

    auto g_z_0_0_0_yy_xz_y_yz = buffer_1000_ddpd[1666];

    auto g_z_0_0_0_yy_xz_y_zz = buffer_1000_ddpd[1667];

    auto g_z_0_0_0_yy_xz_z_xx = buffer_1000_ddpd[1668];

    auto g_z_0_0_0_yy_xz_z_xy = buffer_1000_ddpd[1669];

    auto g_z_0_0_0_yy_xz_z_xz = buffer_1000_ddpd[1670];

    auto g_z_0_0_0_yy_xz_z_yy = buffer_1000_ddpd[1671];

    auto g_z_0_0_0_yy_xz_z_yz = buffer_1000_ddpd[1672];

    auto g_z_0_0_0_yy_xz_z_zz = buffer_1000_ddpd[1673];

    auto g_z_0_0_0_yy_yy_x_xx = buffer_1000_ddpd[1674];

    auto g_z_0_0_0_yy_yy_x_xy = buffer_1000_ddpd[1675];

    auto g_z_0_0_0_yy_yy_x_xz = buffer_1000_ddpd[1676];

    auto g_z_0_0_0_yy_yy_x_yy = buffer_1000_ddpd[1677];

    auto g_z_0_0_0_yy_yy_x_yz = buffer_1000_ddpd[1678];

    auto g_z_0_0_0_yy_yy_x_zz = buffer_1000_ddpd[1679];

    auto g_z_0_0_0_yy_yy_y_xx = buffer_1000_ddpd[1680];

    auto g_z_0_0_0_yy_yy_y_xy = buffer_1000_ddpd[1681];

    auto g_z_0_0_0_yy_yy_y_xz = buffer_1000_ddpd[1682];

    auto g_z_0_0_0_yy_yy_y_yy = buffer_1000_ddpd[1683];

    auto g_z_0_0_0_yy_yy_y_yz = buffer_1000_ddpd[1684];

    auto g_z_0_0_0_yy_yy_y_zz = buffer_1000_ddpd[1685];

    auto g_z_0_0_0_yy_yy_z_xx = buffer_1000_ddpd[1686];

    auto g_z_0_0_0_yy_yy_z_xy = buffer_1000_ddpd[1687];

    auto g_z_0_0_0_yy_yy_z_xz = buffer_1000_ddpd[1688];

    auto g_z_0_0_0_yy_yy_z_yy = buffer_1000_ddpd[1689];

    auto g_z_0_0_0_yy_yy_z_yz = buffer_1000_ddpd[1690];

    auto g_z_0_0_0_yy_yy_z_zz = buffer_1000_ddpd[1691];

    auto g_z_0_0_0_yy_yz_x_xx = buffer_1000_ddpd[1692];

    auto g_z_0_0_0_yy_yz_x_xy = buffer_1000_ddpd[1693];

    auto g_z_0_0_0_yy_yz_x_xz = buffer_1000_ddpd[1694];

    auto g_z_0_0_0_yy_yz_x_yy = buffer_1000_ddpd[1695];

    auto g_z_0_0_0_yy_yz_x_yz = buffer_1000_ddpd[1696];

    auto g_z_0_0_0_yy_yz_x_zz = buffer_1000_ddpd[1697];

    auto g_z_0_0_0_yy_yz_y_xx = buffer_1000_ddpd[1698];

    auto g_z_0_0_0_yy_yz_y_xy = buffer_1000_ddpd[1699];

    auto g_z_0_0_0_yy_yz_y_xz = buffer_1000_ddpd[1700];

    auto g_z_0_0_0_yy_yz_y_yy = buffer_1000_ddpd[1701];

    auto g_z_0_0_0_yy_yz_y_yz = buffer_1000_ddpd[1702];

    auto g_z_0_0_0_yy_yz_y_zz = buffer_1000_ddpd[1703];

    auto g_z_0_0_0_yy_yz_z_xx = buffer_1000_ddpd[1704];

    auto g_z_0_0_0_yy_yz_z_xy = buffer_1000_ddpd[1705];

    auto g_z_0_0_0_yy_yz_z_xz = buffer_1000_ddpd[1706];

    auto g_z_0_0_0_yy_yz_z_yy = buffer_1000_ddpd[1707];

    auto g_z_0_0_0_yy_yz_z_yz = buffer_1000_ddpd[1708];

    auto g_z_0_0_0_yy_yz_z_zz = buffer_1000_ddpd[1709];

    auto g_z_0_0_0_yy_zz_x_xx = buffer_1000_ddpd[1710];

    auto g_z_0_0_0_yy_zz_x_xy = buffer_1000_ddpd[1711];

    auto g_z_0_0_0_yy_zz_x_xz = buffer_1000_ddpd[1712];

    auto g_z_0_0_0_yy_zz_x_yy = buffer_1000_ddpd[1713];

    auto g_z_0_0_0_yy_zz_x_yz = buffer_1000_ddpd[1714];

    auto g_z_0_0_0_yy_zz_x_zz = buffer_1000_ddpd[1715];

    auto g_z_0_0_0_yy_zz_y_xx = buffer_1000_ddpd[1716];

    auto g_z_0_0_0_yy_zz_y_xy = buffer_1000_ddpd[1717];

    auto g_z_0_0_0_yy_zz_y_xz = buffer_1000_ddpd[1718];

    auto g_z_0_0_0_yy_zz_y_yy = buffer_1000_ddpd[1719];

    auto g_z_0_0_0_yy_zz_y_yz = buffer_1000_ddpd[1720];

    auto g_z_0_0_0_yy_zz_y_zz = buffer_1000_ddpd[1721];

    auto g_z_0_0_0_yy_zz_z_xx = buffer_1000_ddpd[1722];

    auto g_z_0_0_0_yy_zz_z_xy = buffer_1000_ddpd[1723];

    auto g_z_0_0_0_yy_zz_z_xz = buffer_1000_ddpd[1724];

    auto g_z_0_0_0_yy_zz_z_yy = buffer_1000_ddpd[1725];

    auto g_z_0_0_0_yy_zz_z_yz = buffer_1000_ddpd[1726];

    auto g_z_0_0_0_yy_zz_z_zz = buffer_1000_ddpd[1727];

    auto g_z_0_0_0_yz_xx_x_xx = buffer_1000_ddpd[1728];

    auto g_z_0_0_0_yz_xx_x_xy = buffer_1000_ddpd[1729];

    auto g_z_0_0_0_yz_xx_x_xz = buffer_1000_ddpd[1730];

    auto g_z_0_0_0_yz_xx_x_yy = buffer_1000_ddpd[1731];

    auto g_z_0_0_0_yz_xx_x_yz = buffer_1000_ddpd[1732];

    auto g_z_0_0_0_yz_xx_x_zz = buffer_1000_ddpd[1733];

    auto g_z_0_0_0_yz_xx_y_xx = buffer_1000_ddpd[1734];

    auto g_z_0_0_0_yz_xx_y_xy = buffer_1000_ddpd[1735];

    auto g_z_0_0_0_yz_xx_y_xz = buffer_1000_ddpd[1736];

    auto g_z_0_0_0_yz_xx_y_yy = buffer_1000_ddpd[1737];

    auto g_z_0_0_0_yz_xx_y_yz = buffer_1000_ddpd[1738];

    auto g_z_0_0_0_yz_xx_y_zz = buffer_1000_ddpd[1739];

    auto g_z_0_0_0_yz_xx_z_xx = buffer_1000_ddpd[1740];

    auto g_z_0_0_0_yz_xx_z_xy = buffer_1000_ddpd[1741];

    auto g_z_0_0_0_yz_xx_z_xz = buffer_1000_ddpd[1742];

    auto g_z_0_0_0_yz_xx_z_yy = buffer_1000_ddpd[1743];

    auto g_z_0_0_0_yz_xx_z_yz = buffer_1000_ddpd[1744];

    auto g_z_0_0_0_yz_xx_z_zz = buffer_1000_ddpd[1745];

    auto g_z_0_0_0_yz_xy_x_xx = buffer_1000_ddpd[1746];

    auto g_z_0_0_0_yz_xy_x_xy = buffer_1000_ddpd[1747];

    auto g_z_0_0_0_yz_xy_x_xz = buffer_1000_ddpd[1748];

    auto g_z_0_0_0_yz_xy_x_yy = buffer_1000_ddpd[1749];

    auto g_z_0_0_0_yz_xy_x_yz = buffer_1000_ddpd[1750];

    auto g_z_0_0_0_yz_xy_x_zz = buffer_1000_ddpd[1751];

    auto g_z_0_0_0_yz_xy_y_xx = buffer_1000_ddpd[1752];

    auto g_z_0_0_0_yz_xy_y_xy = buffer_1000_ddpd[1753];

    auto g_z_0_0_0_yz_xy_y_xz = buffer_1000_ddpd[1754];

    auto g_z_0_0_0_yz_xy_y_yy = buffer_1000_ddpd[1755];

    auto g_z_0_0_0_yz_xy_y_yz = buffer_1000_ddpd[1756];

    auto g_z_0_0_0_yz_xy_y_zz = buffer_1000_ddpd[1757];

    auto g_z_0_0_0_yz_xy_z_xx = buffer_1000_ddpd[1758];

    auto g_z_0_0_0_yz_xy_z_xy = buffer_1000_ddpd[1759];

    auto g_z_0_0_0_yz_xy_z_xz = buffer_1000_ddpd[1760];

    auto g_z_0_0_0_yz_xy_z_yy = buffer_1000_ddpd[1761];

    auto g_z_0_0_0_yz_xy_z_yz = buffer_1000_ddpd[1762];

    auto g_z_0_0_0_yz_xy_z_zz = buffer_1000_ddpd[1763];

    auto g_z_0_0_0_yz_xz_x_xx = buffer_1000_ddpd[1764];

    auto g_z_0_0_0_yz_xz_x_xy = buffer_1000_ddpd[1765];

    auto g_z_0_0_0_yz_xz_x_xz = buffer_1000_ddpd[1766];

    auto g_z_0_0_0_yz_xz_x_yy = buffer_1000_ddpd[1767];

    auto g_z_0_0_0_yz_xz_x_yz = buffer_1000_ddpd[1768];

    auto g_z_0_0_0_yz_xz_x_zz = buffer_1000_ddpd[1769];

    auto g_z_0_0_0_yz_xz_y_xx = buffer_1000_ddpd[1770];

    auto g_z_0_0_0_yz_xz_y_xy = buffer_1000_ddpd[1771];

    auto g_z_0_0_0_yz_xz_y_xz = buffer_1000_ddpd[1772];

    auto g_z_0_0_0_yz_xz_y_yy = buffer_1000_ddpd[1773];

    auto g_z_0_0_0_yz_xz_y_yz = buffer_1000_ddpd[1774];

    auto g_z_0_0_0_yz_xz_y_zz = buffer_1000_ddpd[1775];

    auto g_z_0_0_0_yz_xz_z_xx = buffer_1000_ddpd[1776];

    auto g_z_0_0_0_yz_xz_z_xy = buffer_1000_ddpd[1777];

    auto g_z_0_0_0_yz_xz_z_xz = buffer_1000_ddpd[1778];

    auto g_z_0_0_0_yz_xz_z_yy = buffer_1000_ddpd[1779];

    auto g_z_0_0_0_yz_xz_z_yz = buffer_1000_ddpd[1780];

    auto g_z_0_0_0_yz_xz_z_zz = buffer_1000_ddpd[1781];

    auto g_z_0_0_0_yz_yy_x_xx = buffer_1000_ddpd[1782];

    auto g_z_0_0_0_yz_yy_x_xy = buffer_1000_ddpd[1783];

    auto g_z_0_0_0_yz_yy_x_xz = buffer_1000_ddpd[1784];

    auto g_z_0_0_0_yz_yy_x_yy = buffer_1000_ddpd[1785];

    auto g_z_0_0_0_yz_yy_x_yz = buffer_1000_ddpd[1786];

    auto g_z_0_0_0_yz_yy_x_zz = buffer_1000_ddpd[1787];

    auto g_z_0_0_0_yz_yy_y_xx = buffer_1000_ddpd[1788];

    auto g_z_0_0_0_yz_yy_y_xy = buffer_1000_ddpd[1789];

    auto g_z_0_0_0_yz_yy_y_xz = buffer_1000_ddpd[1790];

    auto g_z_0_0_0_yz_yy_y_yy = buffer_1000_ddpd[1791];

    auto g_z_0_0_0_yz_yy_y_yz = buffer_1000_ddpd[1792];

    auto g_z_0_0_0_yz_yy_y_zz = buffer_1000_ddpd[1793];

    auto g_z_0_0_0_yz_yy_z_xx = buffer_1000_ddpd[1794];

    auto g_z_0_0_0_yz_yy_z_xy = buffer_1000_ddpd[1795];

    auto g_z_0_0_0_yz_yy_z_xz = buffer_1000_ddpd[1796];

    auto g_z_0_0_0_yz_yy_z_yy = buffer_1000_ddpd[1797];

    auto g_z_0_0_0_yz_yy_z_yz = buffer_1000_ddpd[1798];

    auto g_z_0_0_0_yz_yy_z_zz = buffer_1000_ddpd[1799];

    auto g_z_0_0_0_yz_yz_x_xx = buffer_1000_ddpd[1800];

    auto g_z_0_0_0_yz_yz_x_xy = buffer_1000_ddpd[1801];

    auto g_z_0_0_0_yz_yz_x_xz = buffer_1000_ddpd[1802];

    auto g_z_0_0_0_yz_yz_x_yy = buffer_1000_ddpd[1803];

    auto g_z_0_0_0_yz_yz_x_yz = buffer_1000_ddpd[1804];

    auto g_z_0_0_0_yz_yz_x_zz = buffer_1000_ddpd[1805];

    auto g_z_0_0_0_yz_yz_y_xx = buffer_1000_ddpd[1806];

    auto g_z_0_0_0_yz_yz_y_xy = buffer_1000_ddpd[1807];

    auto g_z_0_0_0_yz_yz_y_xz = buffer_1000_ddpd[1808];

    auto g_z_0_0_0_yz_yz_y_yy = buffer_1000_ddpd[1809];

    auto g_z_0_0_0_yz_yz_y_yz = buffer_1000_ddpd[1810];

    auto g_z_0_0_0_yz_yz_y_zz = buffer_1000_ddpd[1811];

    auto g_z_0_0_0_yz_yz_z_xx = buffer_1000_ddpd[1812];

    auto g_z_0_0_0_yz_yz_z_xy = buffer_1000_ddpd[1813];

    auto g_z_0_0_0_yz_yz_z_xz = buffer_1000_ddpd[1814];

    auto g_z_0_0_0_yz_yz_z_yy = buffer_1000_ddpd[1815];

    auto g_z_0_0_0_yz_yz_z_yz = buffer_1000_ddpd[1816];

    auto g_z_0_0_0_yz_yz_z_zz = buffer_1000_ddpd[1817];

    auto g_z_0_0_0_yz_zz_x_xx = buffer_1000_ddpd[1818];

    auto g_z_0_0_0_yz_zz_x_xy = buffer_1000_ddpd[1819];

    auto g_z_0_0_0_yz_zz_x_xz = buffer_1000_ddpd[1820];

    auto g_z_0_0_0_yz_zz_x_yy = buffer_1000_ddpd[1821];

    auto g_z_0_0_0_yz_zz_x_yz = buffer_1000_ddpd[1822];

    auto g_z_0_0_0_yz_zz_x_zz = buffer_1000_ddpd[1823];

    auto g_z_0_0_0_yz_zz_y_xx = buffer_1000_ddpd[1824];

    auto g_z_0_0_0_yz_zz_y_xy = buffer_1000_ddpd[1825];

    auto g_z_0_0_0_yz_zz_y_xz = buffer_1000_ddpd[1826];

    auto g_z_0_0_0_yz_zz_y_yy = buffer_1000_ddpd[1827];

    auto g_z_0_0_0_yz_zz_y_yz = buffer_1000_ddpd[1828];

    auto g_z_0_0_0_yz_zz_y_zz = buffer_1000_ddpd[1829];

    auto g_z_0_0_0_yz_zz_z_xx = buffer_1000_ddpd[1830];

    auto g_z_0_0_0_yz_zz_z_xy = buffer_1000_ddpd[1831];

    auto g_z_0_0_0_yz_zz_z_xz = buffer_1000_ddpd[1832];

    auto g_z_0_0_0_yz_zz_z_yy = buffer_1000_ddpd[1833];

    auto g_z_0_0_0_yz_zz_z_yz = buffer_1000_ddpd[1834];

    auto g_z_0_0_0_yz_zz_z_zz = buffer_1000_ddpd[1835];

    auto g_z_0_0_0_zz_xx_x_xx = buffer_1000_ddpd[1836];

    auto g_z_0_0_0_zz_xx_x_xy = buffer_1000_ddpd[1837];

    auto g_z_0_0_0_zz_xx_x_xz = buffer_1000_ddpd[1838];

    auto g_z_0_0_0_zz_xx_x_yy = buffer_1000_ddpd[1839];

    auto g_z_0_0_0_zz_xx_x_yz = buffer_1000_ddpd[1840];

    auto g_z_0_0_0_zz_xx_x_zz = buffer_1000_ddpd[1841];

    auto g_z_0_0_0_zz_xx_y_xx = buffer_1000_ddpd[1842];

    auto g_z_0_0_0_zz_xx_y_xy = buffer_1000_ddpd[1843];

    auto g_z_0_0_0_zz_xx_y_xz = buffer_1000_ddpd[1844];

    auto g_z_0_0_0_zz_xx_y_yy = buffer_1000_ddpd[1845];

    auto g_z_0_0_0_zz_xx_y_yz = buffer_1000_ddpd[1846];

    auto g_z_0_0_0_zz_xx_y_zz = buffer_1000_ddpd[1847];

    auto g_z_0_0_0_zz_xx_z_xx = buffer_1000_ddpd[1848];

    auto g_z_0_0_0_zz_xx_z_xy = buffer_1000_ddpd[1849];

    auto g_z_0_0_0_zz_xx_z_xz = buffer_1000_ddpd[1850];

    auto g_z_0_0_0_zz_xx_z_yy = buffer_1000_ddpd[1851];

    auto g_z_0_0_0_zz_xx_z_yz = buffer_1000_ddpd[1852];

    auto g_z_0_0_0_zz_xx_z_zz = buffer_1000_ddpd[1853];

    auto g_z_0_0_0_zz_xy_x_xx = buffer_1000_ddpd[1854];

    auto g_z_0_0_0_zz_xy_x_xy = buffer_1000_ddpd[1855];

    auto g_z_0_0_0_zz_xy_x_xz = buffer_1000_ddpd[1856];

    auto g_z_0_0_0_zz_xy_x_yy = buffer_1000_ddpd[1857];

    auto g_z_0_0_0_zz_xy_x_yz = buffer_1000_ddpd[1858];

    auto g_z_0_0_0_zz_xy_x_zz = buffer_1000_ddpd[1859];

    auto g_z_0_0_0_zz_xy_y_xx = buffer_1000_ddpd[1860];

    auto g_z_0_0_0_zz_xy_y_xy = buffer_1000_ddpd[1861];

    auto g_z_0_0_0_zz_xy_y_xz = buffer_1000_ddpd[1862];

    auto g_z_0_0_0_zz_xy_y_yy = buffer_1000_ddpd[1863];

    auto g_z_0_0_0_zz_xy_y_yz = buffer_1000_ddpd[1864];

    auto g_z_0_0_0_zz_xy_y_zz = buffer_1000_ddpd[1865];

    auto g_z_0_0_0_zz_xy_z_xx = buffer_1000_ddpd[1866];

    auto g_z_0_0_0_zz_xy_z_xy = buffer_1000_ddpd[1867];

    auto g_z_0_0_0_zz_xy_z_xz = buffer_1000_ddpd[1868];

    auto g_z_0_0_0_zz_xy_z_yy = buffer_1000_ddpd[1869];

    auto g_z_0_0_0_zz_xy_z_yz = buffer_1000_ddpd[1870];

    auto g_z_0_0_0_zz_xy_z_zz = buffer_1000_ddpd[1871];

    auto g_z_0_0_0_zz_xz_x_xx = buffer_1000_ddpd[1872];

    auto g_z_0_0_0_zz_xz_x_xy = buffer_1000_ddpd[1873];

    auto g_z_0_0_0_zz_xz_x_xz = buffer_1000_ddpd[1874];

    auto g_z_0_0_0_zz_xz_x_yy = buffer_1000_ddpd[1875];

    auto g_z_0_0_0_zz_xz_x_yz = buffer_1000_ddpd[1876];

    auto g_z_0_0_0_zz_xz_x_zz = buffer_1000_ddpd[1877];

    auto g_z_0_0_0_zz_xz_y_xx = buffer_1000_ddpd[1878];

    auto g_z_0_0_0_zz_xz_y_xy = buffer_1000_ddpd[1879];

    auto g_z_0_0_0_zz_xz_y_xz = buffer_1000_ddpd[1880];

    auto g_z_0_0_0_zz_xz_y_yy = buffer_1000_ddpd[1881];

    auto g_z_0_0_0_zz_xz_y_yz = buffer_1000_ddpd[1882];

    auto g_z_0_0_0_zz_xz_y_zz = buffer_1000_ddpd[1883];

    auto g_z_0_0_0_zz_xz_z_xx = buffer_1000_ddpd[1884];

    auto g_z_0_0_0_zz_xz_z_xy = buffer_1000_ddpd[1885];

    auto g_z_0_0_0_zz_xz_z_xz = buffer_1000_ddpd[1886];

    auto g_z_0_0_0_zz_xz_z_yy = buffer_1000_ddpd[1887];

    auto g_z_0_0_0_zz_xz_z_yz = buffer_1000_ddpd[1888];

    auto g_z_0_0_0_zz_xz_z_zz = buffer_1000_ddpd[1889];

    auto g_z_0_0_0_zz_yy_x_xx = buffer_1000_ddpd[1890];

    auto g_z_0_0_0_zz_yy_x_xy = buffer_1000_ddpd[1891];

    auto g_z_0_0_0_zz_yy_x_xz = buffer_1000_ddpd[1892];

    auto g_z_0_0_0_zz_yy_x_yy = buffer_1000_ddpd[1893];

    auto g_z_0_0_0_zz_yy_x_yz = buffer_1000_ddpd[1894];

    auto g_z_0_0_0_zz_yy_x_zz = buffer_1000_ddpd[1895];

    auto g_z_0_0_0_zz_yy_y_xx = buffer_1000_ddpd[1896];

    auto g_z_0_0_0_zz_yy_y_xy = buffer_1000_ddpd[1897];

    auto g_z_0_0_0_zz_yy_y_xz = buffer_1000_ddpd[1898];

    auto g_z_0_0_0_zz_yy_y_yy = buffer_1000_ddpd[1899];

    auto g_z_0_0_0_zz_yy_y_yz = buffer_1000_ddpd[1900];

    auto g_z_0_0_0_zz_yy_y_zz = buffer_1000_ddpd[1901];

    auto g_z_0_0_0_zz_yy_z_xx = buffer_1000_ddpd[1902];

    auto g_z_0_0_0_zz_yy_z_xy = buffer_1000_ddpd[1903];

    auto g_z_0_0_0_zz_yy_z_xz = buffer_1000_ddpd[1904];

    auto g_z_0_0_0_zz_yy_z_yy = buffer_1000_ddpd[1905];

    auto g_z_0_0_0_zz_yy_z_yz = buffer_1000_ddpd[1906];

    auto g_z_0_0_0_zz_yy_z_zz = buffer_1000_ddpd[1907];

    auto g_z_0_0_0_zz_yz_x_xx = buffer_1000_ddpd[1908];

    auto g_z_0_0_0_zz_yz_x_xy = buffer_1000_ddpd[1909];

    auto g_z_0_0_0_zz_yz_x_xz = buffer_1000_ddpd[1910];

    auto g_z_0_0_0_zz_yz_x_yy = buffer_1000_ddpd[1911];

    auto g_z_0_0_0_zz_yz_x_yz = buffer_1000_ddpd[1912];

    auto g_z_0_0_0_zz_yz_x_zz = buffer_1000_ddpd[1913];

    auto g_z_0_0_0_zz_yz_y_xx = buffer_1000_ddpd[1914];

    auto g_z_0_0_0_zz_yz_y_xy = buffer_1000_ddpd[1915];

    auto g_z_0_0_0_zz_yz_y_xz = buffer_1000_ddpd[1916];

    auto g_z_0_0_0_zz_yz_y_yy = buffer_1000_ddpd[1917];

    auto g_z_0_0_0_zz_yz_y_yz = buffer_1000_ddpd[1918];

    auto g_z_0_0_0_zz_yz_y_zz = buffer_1000_ddpd[1919];

    auto g_z_0_0_0_zz_yz_z_xx = buffer_1000_ddpd[1920];

    auto g_z_0_0_0_zz_yz_z_xy = buffer_1000_ddpd[1921];

    auto g_z_0_0_0_zz_yz_z_xz = buffer_1000_ddpd[1922];

    auto g_z_0_0_0_zz_yz_z_yy = buffer_1000_ddpd[1923];

    auto g_z_0_0_0_zz_yz_z_yz = buffer_1000_ddpd[1924];

    auto g_z_0_0_0_zz_yz_z_zz = buffer_1000_ddpd[1925];

    auto g_z_0_0_0_zz_zz_x_xx = buffer_1000_ddpd[1926];

    auto g_z_0_0_0_zz_zz_x_xy = buffer_1000_ddpd[1927];

    auto g_z_0_0_0_zz_zz_x_xz = buffer_1000_ddpd[1928];

    auto g_z_0_0_0_zz_zz_x_yy = buffer_1000_ddpd[1929];

    auto g_z_0_0_0_zz_zz_x_yz = buffer_1000_ddpd[1930];

    auto g_z_0_0_0_zz_zz_x_zz = buffer_1000_ddpd[1931];

    auto g_z_0_0_0_zz_zz_y_xx = buffer_1000_ddpd[1932];

    auto g_z_0_0_0_zz_zz_y_xy = buffer_1000_ddpd[1933];

    auto g_z_0_0_0_zz_zz_y_xz = buffer_1000_ddpd[1934];

    auto g_z_0_0_0_zz_zz_y_yy = buffer_1000_ddpd[1935];

    auto g_z_0_0_0_zz_zz_y_yz = buffer_1000_ddpd[1936];

    auto g_z_0_0_0_zz_zz_y_zz = buffer_1000_ddpd[1937];

    auto g_z_0_0_0_zz_zz_z_xx = buffer_1000_ddpd[1938];

    auto g_z_0_0_0_zz_zz_z_xy = buffer_1000_ddpd[1939];

    auto g_z_0_0_0_zz_zz_z_xz = buffer_1000_ddpd[1940];

    auto g_z_0_0_0_zz_zz_z_yy = buffer_1000_ddpd[1941];

    auto g_z_0_0_0_zz_zz_z_yz = buffer_1000_ddpd[1942];

    auto g_z_0_0_0_zz_zz_z_zz = buffer_1000_ddpd[1943];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_0_0_0_xx_xx_x_xx, g_x_0_0_0_xx_xx_x_xy, g_x_0_0_0_xx_xx_x_xz, g_x_0_0_0_xx_xx_x_yy, g_x_0_0_0_xx_xx_x_yz, g_x_0_0_0_xx_xx_x_zz, g_x_xx_x_xx, g_x_xx_x_xy, g_x_xx_x_xz, g_x_xx_x_yy, g_x_xx_x_yz, g_x_xx_x_zz, g_xxx_xx_x_xx, g_xxx_xx_x_xy, g_xxx_xx_x_xz, g_xxx_xx_x_yy, g_xxx_xx_x_yz, g_xxx_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_xx_x_xx[i] = -2.0 * g_x_xx_x_xx[i] + 2.0 * g_xxx_xx_x_xx[i] * a_exp;

        g_x_0_0_0_xx_xx_x_xy[i] = -2.0 * g_x_xx_x_xy[i] + 2.0 * g_xxx_xx_x_xy[i] * a_exp;

        g_x_0_0_0_xx_xx_x_xz[i] = -2.0 * g_x_xx_x_xz[i] + 2.0 * g_xxx_xx_x_xz[i] * a_exp;

        g_x_0_0_0_xx_xx_x_yy[i] = -2.0 * g_x_xx_x_yy[i] + 2.0 * g_xxx_xx_x_yy[i] * a_exp;

        g_x_0_0_0_xx_xx_x_yz[i] = -2.0 * g_x_xx_x_yz[i] + 2.0 * g_xxx_xx_x_yz[i] * a_exp;

        g_x_0_0_0_xx_xx_x_zz[i] = -2.0 * g_x_xx_x_zz[i] + 2.0 * g_xxx_xx_x_zz[i] * a_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_0_0_0_xx_xx_y_xx, g_x_0_0_0_xx_xx_y_xy, g_x_0_0_0_xx_xx_y_xz, g_x_0_0_0_xx_xx_y_yy, g_x_0_0_0_xx_xx_y_yz, g_x_0_0_0_xx_xx_y_zz, g_x_xx_y_xx, g_x_xx_y_xy, g_x_xx_y_xz, g_x_xx_y_yy, g_x_xx_y_yz, g_x_xx_y_zz, g_xxx_xx_y_xx, g_xxx_xx_y_xy, g_xxx_xx_y_xz, g_xxx_xx_y_yy, g_xxx_xx_y_yz, g_xxx_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_xx_y_xx[i] = -2.0 * g_x_xx_y_xx[i] + 2.0 * g_xxx_xx_y_xx[i] * a_exp;

        g_x_0_0_0_xx_xx_y_xy[i] = -2.0 * g_x_xx_y_xy[i] + 2.0 * g_xxx_xx_y_xy[i] * a_exp;

        g_x_0_0_0_xx_xx_y_xz[i] = -2.0 * g_x_xx_y_xz[i] + 2.0 * g_xxx_xx_y_xz[i] * a_exp;

        g_x_0_0_0_xx_xx_y_yy[i] = -2.0 * g_x_xx_y_yy[i] + 2.0 * g_xxx_xx_y_yy[i] * a_exp;

        g_x_0_0_0_xx_xx_y_yz[i] = -2.0 * g_x_xx_y_yz[i] + 2.0 * g_xxx_xx_y_yz[i] * a_exp;

        g_x_0_0_0_xx_xx_y_zz[i] = -2.0 * g_x_xx_y_zz[i] + 2.0 * g_xxx_xx_y_zz[i] * a_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_0_0_0_xx_xx_z_xx, g_x_0_0_0_xx_xx_z_xy, g_x_0_0_0_xx_xx_z_xz, g_x_0_0_0_xx_xx_z_yy, g_x_0_0_0_xx_xx_z_yz, g_x_0_0_0_xx_xx_z_zz, g_x_xx_z_xx, g_x_xx_z_xy, g_x_xx_z_xz, g_x_xx_z_yy, g_x_xx_z_yz, g_x_xx_z_zz, g_xxx_xx_z_xx, g_xxx_xx_z_xy, g_xxx_xx_z_xz, g_xxx_xx_z_yy, g_xxx_xx_z_yz, g_xxx_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_xx_z_xx[i] = -2.0 * g_x_xx_z_xx[i] + 2.0 * g_xxx_xx_z_xx[i] * a_exp;

        g_x_0_0_0_xx_xx_z_xy[i] = -2.0 * g_x_xx_z_xy[i] + 2.0 * g_xxx_xx_z_xy[i] * a_exp;

        g_x_0_0_0_xx_xx_z_xz[i] = -2.0 * g_x_xx_z_xz[i] + 2.0 * g_xxx_xx_z_xz[i] * a_exp;

        g_x_0_0_0_xx_xx_z_yy[i] = -2.0 * g_x_xx_z_yy[i] + 2.0 * g_xxx_xx_z_yy[i] * a_exp;

        g_x_0_0_0_xx_xx_z_yz[i] = -2.0 * g_x_xx_z_yz[i] + 2.0 * g_xxx_xx_z_yz[i] * a_exp;

        g_x_0_0_0_xx_xx_z_zz[i] = -2.0 * g_x_xx_z_zz[i] + 2.0 * g_xxx_xx_z_zz[i] * a_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_x_0_0_0_xx_xy_x_xx, g_x_0_0_0_xx_xy_x_xy, g_x_0_0_0_xx_xy_x_xz, g_x_0_0_0_xx_xy_x_yy, g_x_0_0_0_xx_xy_x_yz, g_x_0_0_0_xx_xy_x_zz, g_x_xy_x_xx, g_x_xy_x_xy, g_x_xy_x_xz, g_x_xy_x_yy, g_x_xy_x_yz, g_x_xy_x_zz, g_xxx_xy_x_xx, g_xxx_xy_x_xy, g_xxx_xy_x_xz, g_xxx_xy_x_yy, g_xxx_xy_x_yz, g_xxx_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_xy_x_xx[i] = -2.0 * g_x_xy_x_xx[i] + 2.0 * g_xxx_xy_x_xx[i] * a_exp;

        g_x_0_0_0_xx_xy_x_xy[i] = -2.0 * g_x_xy_x_xy[i] + 2.0 * g_xxx_xy_x_xy[i] * a_exp;

        g_x_0_0_0_xx_xy_x_xz[i] = -2.0 * g_x_xy_x_xz[i] + 2.0 * g_xxx_xy_x_xz[i] * a_exp;

        g_x_0_0_0_xx_xy_x_yy[i] = -2.0 * g_x_xy_x_yy[i] + 2.0 * g_xxx_xy_x_yy[i] * a_exp;

        g_x_0_0_0_xx_xy_x_yz[i] = -2.0 * g_x_xy_x_yz[i] + 2.0 * g_xxx_xy_x_yz[i] * a_exp;

        g_x_0_0_0_xx_xy_x_zz[i] = -2.0 * g_x_xy_x_zz[i] + 2.0 * g_xxx_xy_x_zz[i] * a_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_x_0_0_0_xx_xy_y_xx, g_x_0_0_0_xx_xy_y_xy, g_x_0_0_0_xx_xy_y_xz, g_x_0_0_0_xx_xy_y_yy, g_x_0_0_0_xx_xy_y_yz, g_x_0_0_0_xx_xy_y_zz, g_x_xy_y_xx, g_x_xy_y_xy, g_x_xy_y_xz, g_x_xy_y_yy, g_x_xy_y_yz, g_x_xy_y_zz, g_xxx_xy_y_xx, g_xxx_xy_y_xy, g_xxx_xy_y_xz, g_xxx_xy_y_yy, g_xxx_xy_y_yz, g_xxx_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_xy_y_xx[i] = -2.0 * g_x_xy_y_xx[i] + 2.0 * g_xxx_xy_y_xx[i] * a_exp;

        g_x_0_0_0_xx_xy_y_xy[i] = -2.0 * g_x_xy_y_xy[i] + 2.0 * g_xxx_xy_y_xy[i] * a_exp;

        g_x_0_0_0_xx_xy_y_xz[i] = -2.0 * g_x_xy_y_xz[i] + 2.0 * g_xxx_xy_y_xz[i] * a_exp;

        g_x_0_0_0_xx_xy_y_yy[i] = -2.0 * g_x_xy_y_yy[i] + 2.0 * g_xxx_xy_y_yy[i] * a_exp;

        g_x_0_0_0_xx_xy_y_yz[i] = -2.0 * g_x_xy_y_yz[i] + 2.0 * g_xxx_xy_y_yz[i] * a_exp;

        g_x_0_0_0_xx_xy_y_zz[i] = -2.0 * g_x_xy_y_zz[i] + 2.0 * g_xxx_xy_y_zz[i] * a_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_x_0_0_0_xx_xy_z_xx, g_x_0_0_0_xx_xy_z_xy, g_x_0_0_0_xx_xy_z_xz, g_x_0_0_0_xx_xy_z_yy, g_x_0_0_0_xx_xy_z_yz, g_x_0_0_0_xx_xy_z_zz, g_x_xy_z_xx, g_x_xy_z_xy, g_x_xy_z_xz, g_x_xy_z_yy, g_x_xy_z_yz, g_x_xy_z_zz, g_xxx_xy_z_xx, g_xxx_xy_z_xy, g_xxx_xy_z_xz, g_xxx_xy_z_yy, g_xxx_xy_z_yz, g_xxx_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_xy_z_xx[i] = -2.0 * g_x_xy_z_xx[i] + 2.0 * g_xxx_xy_z_xx[i] * a_exp;

        g_x_0_0_0_xx_xy_z_xy[i] = -2.0 * g_x_xy_z_xy[i] + 2.0 * g_xxx_xy_z_xy[i] * a_exp;

        g_x_0_0_0_xx_xy_z_xz[i] = -2.0 * g_x_xy_z_xz[i] + 2.0 * g_xxx_xy_z_xz[i] * a_exp;

        g_x_0_0_0_xx_xy_z_yy[i] = -2.0 * g_x_xy_z_yy[i] + 2.0 * g_xxx_xy_z_yy[i] * a_exp;

        g_x_0_0_0_xx_xy_z_yz[i] = -2.0 * g_x_xy_z_yz[i] + 2.0 * g_xxx_xy_z_yz[i] * a_exp;

        g_x_0_0_0_xx_xy_z_zz[i] = -2.0 * g_x_xy_z_zz[i] + 2.0 * g_xxx_xy_z_zz[i] * a_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_x_0_0_0_xx_xz_x_xx, g_x_0_0_0_xx_xz_x_xy, g_x_0_0_0_xx_xz_x_xz, g_x_0_0_0_xx_xz_x_yy, g_x_0_0_0_xx_xz_x_yz, g_x_0_0_0_xx_xz_x_zz, g_x_xz_x_xx, g_x_xz_x_xy, g_x_xz_x_xz, g_x_xz_x_yy, g_x_xz_x_yz, g_x_xz_x_zz, g_xxx_xz_x_xx, g_xxx_xz_x_xy, g_xxx_xz_x_xz, g_xxx_xz_x_yy, g_xxx_xz_x_yz, g_xxx_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_xz_x_xx[i] = -2.0 * g_x_xz_x_xx[i] + 2.0 * g_xxx_xz_x_xx[i] * a_exp;

        g_x_0_0_0_xx_xz_x_xy[i] = -2.0 * g_x_xz_x_xy[i] + 2.0 * g_xxx_xz_x_xy[i] * a_exp;

        g_x_0_0_0_xx_xz_x_xz[i] = -2.0 * g_x_xz_x_xz[i] + 2.0 * g_xxx_xz_x_xz[i] * a_exp;

        g_x_0_0_0_xx_xz_x_yy[i] = -2.0 * g_x_xz_x_yy[i] + 2.0 * g_xxx_xz_x_yy[i] * a_exp;

        g_x_0_0_0_xx_xz_x_yz[i] = -2.0 * g_x_xz_x_yz[i] + 2.0 * g_xxx_xz_x_yz[i] * a_exp;

        g_x_0_0_0_xx_xz_x_zz[i] = -2.0 * g_x_xz_x_zz[i] + 2.0 * g_xxx_xz_x_zz[i] * a_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_x_0_0_0_xx_xz_y_xx, g_x_0_0_0_xx_xz_y_xy, g_x_0_0_0_xx_xz_y_xz, g_x_0_0_0_xx_xz_y_yy, g_x_0_0_0_xx_xz_y_yz, g_x_0_0_0_xx_xz_y_zz, g_x_xz_y_xx, g_x_xz_y_xy, g_x_xz_y_xz, g_x_xz_y_yy, g_x_xz_y_yz, g_x_xz_y_zz, g_xxx_xz_y_xx, g_xxx_xz_y_xy, g_xxx_xz_y_xz, g_xxx_xz_y_yy, g_xxx_xz_y_yz, g_xxx_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_xz_y_xx[i] = -2.0 * g_x_xz_y_xx[i] + 2.0 * g_xxx_xz_y_xx[i] * a_exp;

        g_x_0_0_0_xx_xz_y_xy[i] = -2.0 * g_x_xz_y_xy[i] + 2.0 * g_xxx_xz_y_xy[i] * a_exp;

        g_x_0_0_0_xx_xz_y_xz[i] = -2.0 * g_x_xz_y_xz[i] + 2.0 * g_xxx_xz_y_xz[i] * a_exp;

        g_x_0_0_0_xx_xz_y_yy[i] = -2.0 * g_x_xz_y_yy[i] + 2.0 * g_xxx_xz_y_yy[i] * a_exp;

        g_x_0_0_0_xx_xz_y_yz[i] = -2.0 * g_x_xz_y_yz[i] + 2.0 * g_xxx_xz_y_yz[i] * a_exp;

        g_x_0_0_0_xx_xz_y_zz[i] = -2.0 * g_x_xz_y_zz[i] + 2.0 * g_xxx_xz_y_zz[i] * a_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_x_0_0_0_xx_xz_z_xx, g_x_0_0_0_xx_xz_z_xy, g_x_0_0_0_xx_xz_z_xz, g_x_0_0_0_xx_xz_z_yy, g_x_0_0_0_xx_xz_z_yz, g_x_0_0_0_xx_xz_z_zz, g_x_xz_z_xx, g_x_xz_z_xy, g_x_xz_z_xz, g_x_xz_z_yy, g_x_xz_z_yz, g_x_xz_z_zz, g_xxx_xz_z_xx, g_xxx_xz_z_xy, g_xxx_xz_z_xz, g_xxx_xz_z_yy, g_xxx_xz_z_yz, g_xxx_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_xz_z_xx[i] = -2.0 * g_x_xz_z_xx[i] + 2.0 * g_xxx_xz_z_xx[i] * a_exp;

        g_x_0_0_0_xx_xz_z_xy[i] = -2.0 * g_x_xz_z_xy[i] + 2.0 * g_xxx_xz_z_xy[i] * a_exp;

        g_x_0_0_0_xx_xz_z_xz[i] = -2.0 * g_x_xz_z_xz[i] + 2.0 * g_xxx_xz_z_xz[i] * a_exp;

        g_x_0_0_0_xx_xz_z_yy[i] = -2.0 * g_x_xz_z_yy[i] + 2.0 * g_xxx_xz_z_yy[i] * a_exp;

        g_x_0_0_0_xx_xz_z_yz[i] = -2.0 * g_x_xz_z_yz[i] + 2.0 * g_xxx_xz_z_yz[i] * a_exp;

        g_x_0_0_0_xx_xz_z_zz[i] = -2.0 * g_x_xz_z_zz[i] + 2.0 * g_xxx_xz_z_zz[i] * a_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_x_0_0_0_xx_yy_x_xx, g_x_0_0_0_xx_yy_x_xy, g_x_0_0_0_xx_yy_x_xz, g_x_0_0_0_xx_yy_x_yy, g_x_0_0_0_xx_yy_x_yz, g_x_0_0_0_xx_yy_x_zz, g_x_yy_x_xx, g_x_yy_x_xy, g_x_yy_x_xz, g_x_yy_x_yy, g_x_yy_x_yz, g_x_yy_x_zz, g_xxx_yy_x_xx, g_xxx_yy_x_xy, g_xxx_yy_x_xz, g_xxx_yy_x_yy, g_xxx_yy_x_yz, g_xxx_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_yy_x_xx[i] = -2.0 * g_x_yy_x_xx[i] + 2.0 * g_xxx_yy_x_xx[i] * a_exp;

        g_x_0_0_0_xx_yy_x_xy[i] = -2.0 * g_x_yy_x_xy[i] + 2.0 * g_xxx_yy_x_xy[i] * a_exp;

        g_x_0_0_0_xx_yy_x_xz[i] = -2.0 * g_x_yy_x_xz[i] + 2.0 * g_xxx_yy_x_xz[i] * a_exp;

        g_x_0_0_0_xx_yy_x_yy[i] = -2.0 * g_x_yy_x_yy[i] + 2.0 * g_xxx_yy_x_yy[i] * a_exp;

        g_x_0_0_0_xx_yy_x_yz[i] = -2.0 * g_x_yy_x_yz[i] + 2.0 * g_xxx_yy_x_yz[i] * a_exp;

        g_x_0_0_0_xx_yy_x_zz[i] = -2.0 * g_x_yy_x_zz[i] + 2.0 * g_xxx_yy_x_zz[i] * a_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_x_0_0_0_xx_yy_y_xx, g_x_0_0_0_xx_yy_y_xy, g_x_0_0_0_xx_yy_y_xz, g_x_0_0_0_xx_yy_y_yy, g_x_0_0_0_xx_yy_y_yz, g_x_0_0_0_xx_yy_y_zz, g_x_yy_y_xx, g_x_yy_y_xy, g_x_yy_y_xz, g_x_yy_y_yy, g_x_yy_y_yz, g_x_yy_y_zz, g_xxx_yy_y_xx, g_xxx_yy_y_xy, g_xxx_yy_y_xz, g_xxx_yy_y_yy, g_xxx_yy_y_yz, g_xxx_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_yy_y_xx[i] = -2.0 * g_x_yy_y_xx[i] + 2.0 * g_xxx_yy_y_xx[i] * a_exp;

        g_x_0_0_0_xx_yy_y_xy[i] = -2.0 * g_x_yy_y_xy[i] + 2.0 * g_xxx_yy_y_xy[i] * a_exp;

        g_x_0_0_0_xx_yy_y_xz[i] = -2.0 * g_x_yy_y_xz[i] + 2.0 * g_xxx_yy_y_xz[i] * a_exp;

        g_x_0_0_0_xx_yy_y_yy[i] = -2.0 * g_x_yy_y_yy[i] + 2.0 * g_xxx_yy_y_yy[i] * a_exp;

        g_x_0_0_0_xx_yy_y_yz[i] = -2.0 * g_x_yy_y_yz[i] + 2.0 * g_xxx_yy_y_yz[i] * a_exp;

        g_x_0_0_0_xx_yy_y_zz[i] = -2.0 * g_x_yy_y_zz[i] + 2.0 * g_xxx_yy_y_zz[i] * a_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_x_0_0_0_xx_yy_z_xx, g_x_0_0_0_xx_yy_z_xy, g_x_0_0_0_xx_yy_z_xz, g_x_0_0_0_xx_yy_z_yy, g_x_0_0_0_xx_yy_z_yz, g_x_0_0_0_xx_yy_z_zz, g_x_yy_z_xx, g_x_yy_z_xy, g_x_yy_z_xz, g_x_yy_z_yy, g_x_yy_z_yz, g_x_yy_z_zz, g_xxx_yy_z_xx, g_xxx_yy_z_xy, g_xxx_yy_z_xz, g_xxx_yy_z_yy, g_xxx_yy_z_yz, g_xxx_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_yy_z_xx[i] = -2.0 * g_x_yy_z_xx[i] + 2.0 * g_xxx_yy_z_xx[i] * a_exp;

        g_x_0_0_0_xx_yy_z_xy[i] = -2.0 * g_x_yy_z_xy[i] + 2.0 * g_xxx_yy_z_xy[i] * a_exp;

        g_x_0_0_0_xx_yy_z_xz[i] = -2.0 * g_x_yy_z_xz[i] + 2.0 * g_xxx_yy_z_xz[i] * a_exp;

        g_x_0_0_0_xx_yy_z_yy[i] = -2.0 * g_x_yy_z_yy[i] + 2.0 * g_xxx_yy_z_yy[i] * a_exp;

        g_x_0_0_0_xx_yy_z_yz[i] = -2.0 * g_x_yy_z_yz[i] + 2.0 * g_xxx_yy_z_yz[i] * a_exp;

        g_x_0_0_0_xx_yy_z_zz[i] = -2.0 * g_x_yy_z_zz[i] + 2.0 * g_xxx_yy_z_zz[i] * a_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_x_0_0_0_xx_yz_x_xx, g_x_0_0_0_xx_yz_x_xy, g_x_0_0_0_xx_yz_x_xz, g_x_0_0_0_xx_yz_x_yy, g_x_0_0_0_xx_yz_x_yz, g_x_0_0_0_xx_yz_x_zz, g_x_yz_x_xx, g_x_yz_x_xy, g_x_yz_x_xz, g_x_yz_x_yy, g_x_yz_x_yz, g_x_yz_x_zz, g_xxx_yz_x_xx, g_xxx_yz_x_xy, g_xxx_yz_x_xz, g_xxx_yz_x_yy, g_xxx_yz_x_yz, g_xxx_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_yz_x_xx[i] = -2.0 * g_x_yz_x_xx[i] + 2.0 * g_xxx_yz_x_xx[i] * a_exp;

        g_x_0_0_0_xx_yz_x_xy[i] = -2.0 * g_x_yz_x_xy[i] + 2.0 * g_xxx_yz_x_xy[i] * a_exp;

        g_x_0_0_0_xx_yz_x_xz[i] = -2.0 * g_x_yz_x_xz[i] + 2.0 * g_xxx_yz_x_xz[i] * a_exp;

        g_x_0_0_0_xx_yz_x_yy[i] = -2.0 * g_x_yz_x_yy[i] + 2.0 * g_xxx_yz_x_yy[i] * a_exp;

        g_x_0_0_0_xx_yz_x_yz[i] = -2.0 * g_x_yz_x_yz[i] + 2.0 * g_xxx_yz_x_yz[i] * a_exp;

        g_x_0_0_0_xx_yz_x_zz[i] = -2.0 * g_x_yz_x_zz[i] + 2.0 * g_xxx_yz_x_zz[i] * a_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_x_0_0_0_xx_yz_y_xx, g_x_0_0_0_xx_yz_y_xy, g_x_0_0_0_xx_yz_y_xz, g_x_0_0_0_xx_yz_y_yy, g_x_0_0_0_xx_yz_y_yz, g_x_0_0_0_xx_yz_y_zz, g_x_yz_y_xx, g_x_yz_y_xy, g_x_yz_y_xz, g_x_yz_y_yy, g_x_yz_y_yz, g_x_yz_y_zz, g_xxx_yz_y_xx, g_xxx_yz_y_xy, g_xxx_yz_y_xz, g_xxx_yz_y_yy, g_xxx_yz_y_yz, g_xxx_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_yz_y_xx[i] = -2.0 * g_x_yz_y_xx[i] + 2.0 * g_xxx_yz_y_xx[i] * a_exp;

        g_x_0_0_0_xx_yz_y_xy[i] = -2.0 * g_x_yz_y_xy[i] + 2.0 * g_xxx_yz_y_xy[i] * a_exp;

        g_x_0_0_0_xx_yz_y_xz[i] = -2.0 * g_x_yz_y_xz[i] + 2.0 * g_xxx_yz_y_xz[i] * a_exp;

        g_x_0_0_0_xx_yz_y_yy[i] = -2.0 * g_x_yz_y_yy[i] + 2.0 * g_xxx_yz_y_yy[i] * a_exp;

        g_x_0_0_0_xx_yz_y_yz[i] = -2.0 * g_x_yz_y_yz[i] + 2.0 * g_xxx_yz_y_yz[i] * a_exp;

        g_x_0_0_0_xx_yz_y_zz[i] = -2.0 * g_x_yz_y_zz[i] + 2.0 * g_xxx_yz_y_zz[i] * a_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_x_0_0_0_xx_yz_z_xx, g_x_0_0_0_xx_yz_z_xy, g_x_0_0_0_xx_yz_z_xz, g_x_0_0_0_xx_yz_z_yy, g_x_0_0_0_xx_yz_z_yz, g_x_0_0_0_xx_yz_z_zz, g_x_yz_z_xx, g_x_yz_z_xy, g_x_yz_z_xz, g_x_yz_z_yy, g_x_yz_z_yz, g_x_yz_z_zz, g_xxx_yz_z_xx, g_xxx_yz_z_xy, g_xxx_yz_z_xz, g_xxx_yz_z_yy, g_xxx_yz_z_yz, g_xxx_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_yz_z_xx[i] = -2.0 * g_x_yz_z_xx[i] + 2.0 * g_xxx_yz_z_xx[i] * a_exp;

        g_x_0_0_0_xx_yz_z_xy[i] = -2.0 * g_x_yz_z_xy[i] + 2.0 * g_xxx_yz_z_xy[i] * a_exp;

        g_x_0_0_0_xx_yz_z_xz[i] = -2.0 * g_x_yz_z_xz[i] + 2.0 * g_xxx_yz_z_xz[i] * a_exp;

        g_x_0_0_0_xx_yz_z_yy[i] = -2.0 * g_x_yz_z_yy[i] + 2.0 * g_xxx_yz_z_yy[i] * a_exp;

        g_x_0_0_0_xx_yz_z_yz[i] = -2.0 * g_x_yz_z_yz[i] + 2.0 * g_xxx_yz_z_yz[i] * a_exp;

        g_x_0_0_0_xx_yz_z_zz[i] = -2.0 * g_x_yz_z_zz[i] + 2.0 * g_xxx_yz_z_zz[i] * a_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_x_0_0_0_xx_zz_x_xx, g_x_0_0_0_xx_zz_x_xy, g_x_0_0_0_xx_zz_x_xz, g_x_0_0_0_xx_zz_x_yy, g_x_0_0_0_xx_zz_x_yz, g_x_0_0_0_xx_zz_x_zz, g_x_zz_x_xx, g_x_zz_x_xy, g_x_zz_x_xz, g_x_zz_x_yy, g_x_zz_x_yz, g_x_zz_x_zz, g_xxx_zz_x_xx, g_xxx_zz_x_xy, g_xxx_zz_x_xz, g_xxx_zz_x_yy, g_xxx_zz_x_yz, g_xxx_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_zz_x_xx[i] = -2.0 * g_x_zz_x_xx[i] + 2.0 * g_xxx_zz_x_xx[i] * a_exp;

        g_x_0_0_0_xx_zz_x_xy[i] = -2.0 * g_x_zz_x_xy[i] + 2.0 * g_xxx_zz_x_xy[i] * a_exp;

        g_x_0_0_0_xx_zz_x_xz[i] = -2.0 * g_x_zz_x_xz[i] + 2.0 * g_xxx_zz_x_xz[i] * a_exp;

        g_x_0_0_0_xx_zz_x_yy[i] = -2.0 * g_x_zz_x_yy[i] + 2.0 * g_xxx_zz_x_yy[i] * a_exp;

        g_x_0_0_0_xx_zz_x_yz[i] = -2.0 * g_x_zz_x_yz[i] + 2.0 * g_xxx_zz_x_yz[i] * a_exp;

        g_x_0_0_0_xx_zz_x_zz[i] = -2.0 * g_x_zz_x_zz[i] + 2.0 * g_xxx_zz_x_zz[i] * a_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_x_0_0_0_xx_zz_y_xx, g_x_0_0_0_xx_zz_y_xy, g_x_0_0_0_xx_zz_y_xz, g_x_0_0_0_xx_zz_y_yy, g_x_0_0_0_xx_zz_y_yz, g_x_0_0_0_xx_zz_y_zz, g_x_zz_y_xx, g_x_zz_y_xy, g_x_zz_y_xz, g_x_zz_y_yy, g_x_zz_y_yz, g_x_zz_y_zz, g_xxx_zz_y_xx, g_xxx_zz_y_xy, g_xxx_zz_y_xz, g_xxx_zz_y_yy, g_xxx_zz_y_yz, g_xxx_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_zz_y_xx[i] = -2.0 * g_x_zz_y_xx[i] + 2.0 * g_xxx_zz_y_xx[i] * a_exp;

        g_x_0_0_0_xx_zz_y_xy[i] = -2.0 * g_x_zz_y_xy[i] + 2.0 * g_xxx_zz_y_xy[i] * a_exp;

        g_x_0_0_0_xx_zz_y_xz[i] = -2.0 * g_x_zz_y_xz[i] + 2.0 * g_xxx_zz_y_xz[i] * a_exp;

        g_x_0_0_0_xx_zz_y_yy[i] = -2.0 * g_x_zz_y_yy[i] + 2.0 * g_xxx_zz_y_yy[i] * a_exp;

        g_x_0_0_0_xx_zz_y_yz[i] = -2.0 * g_x_zz_y_yz[i] + 2.0 * g_xxx_zz_y_yz[i] * a_exp;

        g_x_0_0_0_xx_zz_y_zz[i] = -2.0 * g_x_zz_y_zz[i] + 2.0 * g_xxx_zz_y_zz[i] * a_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_x_0_0_0_xx_zz_z_xx, g_x_0_0_0_xx_zz_z_xy, g_x_0_0_0_xx_zz_z_xz, g_x_0_0_0_xx_zz_z_yy, g_x_0_0_0_xx_zz_z_yz, g_x_0_0_0_xx_zz_z_zz, g_x_zz_z_xx, g_x_zz_z_xy, g_x_zz_z_xz, g_x_zz_z_yy, g_x_zz_z_yz, g_x_zz_z_zz, g_xxx_zz_z_xx, g_xxx_zz_z_xy, g_xxx_zz_z_xz, g_xxx_zz_z_yy, g_xxx_zz_z_yz, g_xxx_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_zz_z_xx[i] = -2.0 * g_x_zz_z_xx[i] + 2.0 * g_xxx_zz_z_xx[i] * a_exp;

        g_x_0_0_0_xx_zz_z_xy[i] = -2.0 * g_x_zz_z_xy[i] + 2.0 * g_xxx_zz_z_xy[i] * a_exp;

        g_x_0_0_0_xx_zz_z_xz[i] = -2.0 * g_x_zz_z_xz[i] + 2.0 * g_xxx_zz_z_xz[i] * a_exp;

        g_x_0_0_0_xx_zz_z_yy[i] = -2.0 * g_x_zz_z_yy[i] + 2.0 * g_xxx_zz_z_yy[i] * a_exp;

        g_x_0_0_0_xx_zz_z_yz[i] = -2.0 * g_x_zz_z_yz[i] + 2.0 * g_xxx_zz_z_yz[i] * a_exp;

        g_x_0_0_0_xx_zz_z_zz[i] = -2.0 * g_x_zz_z_zz[i] + 2.0 * g_xxx_zz_z_zz[i] * a_exp;
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_x_0_0_0_xy_xx_x_xx, g_x_0_0_0_xy_xx_x_xy, g_x_0_0_0_xy_xx_x_xz, g_x_0_0_0_xy_xx_x_yy, g_x_0_0_0_xy_xx_x_yz, g_x_0_0_0_xy_xx_x_zz, g_xxy_xx_x_xx, g_xxy_xx_x_xy, g_xxy_xx_x_xz, g_xxy_xx_x_yy, g_xxy_xx_x_yz, g_xxy_xx_x_zz, g_y_xx_x_xx, g_y_xx_x_xy, g_y_xx_x_xz, g_y_xx_x_yy, g_y_xx_x_yz, g_y_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_xx_x_xx[i] = -g_y_xx_x_xx[i] + 2.0 * g_xxy_xx_x_xx[i] * a_exp;

        g_x_0_0_0_xy_xx_x_xy[i] = -g_y_xx_x_xy[i] + 2.0 * g_xxy_xx_x_xy[i] * a_exp;

        g_x_0_0_0_xy_xx_x_xz[i] = -g_y_xx_x_xz[i] + 2.0 * g_xxy_xx_x_xz[i] * a_exp;

        g_x_0_0_0_xy_xx_x_yy[i] = -g_y_xx_x_yy[i] + 2.0 * g_xxy_xx_x_yy[i] * a_exp;

        g_x_0_0_0_xy_xx_x_yz[i] = -g_y_xx_x_yz[i] + 2.0 * g_xxy_xx_x_yz[i] * a_exp;

        g_x_0_0_0_xy_xx_x_zz[i] = -g_y_xx_x_zz[i] + 2.0 * g_xxy_xx_x_zz[i] * a_exp;
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_x_0_0_0_xy_xx_y_xx, g_x_0_0_0_xy_xx_y_xy, g_x_0_0_0_xy_xx_y_xz, g_x_0_0_0_xy_xx_y_yy, g_x_0_0_0_xy_xx_y_yz, g_x_0_0_0_xy_xx_y_zz, g_xxy_xx_y_xx, g_xxy_xx_y_xy, g_xxy_xx_y_xz, g_xxy_xx_y_yy, g_xxy_xx_y_yz, g_xxy_xx_y_zz, g_y_xx_y_xx, g_y_xx_y_xy, g_y_xx_y_xz, g_y_xx_y_yy, g_y_xx_y_yz, g_y_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_xx_y_xx[i] = -g_y_xx_y_xx[i] + 2.0 * g_xxy_xx_y_xx[i] * a_exp;

        g_x_0_0_0_xy_xx_y_xy[i] = -g_y_xx_y_xy[i] + 2.0 * g_xxy_xx_y_xy[i] * a_exp;

        g_x_0_0_0_xy_xx_y_xz[i] = -g_y_xx_y_xz[i] + 2.0 * g_xxy_xx_y_xz[i] * a_exp;

        g_x_0_0_0_xy_xx_y_yy[i] = -g_y_xx_y_yy[i] + 2.0 * g_xxy_xx_y_yy[i] * a_exp;

        g_x_0_0_0_xy_xx_y_yz[i] = -g_y_xx_y_yz[i] + 2.0 * g_xxy_xx_y_yz[i] * a_exp;

        g_x_0_0_0_xy_xx_y_zz[i] = -g_y_xx_y_zz[i] + 2.0 * g_xxy_xx_y_zz[i] * a_exp;
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_x_0_0_0_xy_xx_z_xx, g_x_0_0_0_xy_xx_z_xy, g_x_0_0_0_xy_xx_z_xz, g_x_0_0_0_xy_xx_z_yy, g_x_0_0_0_xy_xx_z_yz, g_x_0_0_0_xy_xx_z_zz, g_xxy_xx_z_xx, g_xxy_xx_z_xy, g_xxy_xx_z_xz, g_xxy_xx_z_yy, g_xxy_xx_z_yz, g_xxy_xx_z_zz, g_y_xx_z_xx, g_y_xx_z_xy, g_y_xx_z_xz, g_y_xx_z_yy, g_y_xx_z_yz, g_y_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_xx_z_xx[i] = -g_y_xx_z_xx[i] + 2.0 * g_xxy_xx_z_xx[i] * a_exp;

        g_x_0_0_0_xy_xx_z_xy[i] = -g_y_xx_z_xy[i] + 2.0 * g_xxy_xx_z_xy[i] * a_exp;

        g_x_0_0_0_xy_xx_z_xz[i] = -g_y_xx_z_xz[i] + 2.0 * g_xxy_xx_z_xz[i] * a_exp;

        g_x_0_0_0_xy_xx_z_yy[i] = -g_y_xx_z_yy[i] + 2.0 * g_xxy_xx_z_yy[i] * a_exp;

        g_x_0_0_0_xy_xx_z_yz[i] = -g_y_xx_z_yz[i] + 2.0 * g_xxy_xx_z_yz[i] * a_exp;

        g_x_0_0_0_xy_xx_z_zz[i] = -g_y_xx_z_zz[i] + 2.0 * g_xxy_xx_z_zz[i] * a_exp;
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_x_0_0_0_xy_xy_x_xx, g_x_0_0_0_xy_xy_x_xy, g_x_0_0_0_xy_xy_x_xz, g_x_0_0_0_xy_xy_x_yy, g_x_0_0_0_xy_xy_x_yz, g_x_0_0_0_xy_xy_x_zz, g_xxy_xy_x_xx, g_xxy_xy_x_xy, g_xxy_xy_x_xz, g_xxy_xy_x_yy, g_xxy_xy_x_yz, g_xxy_xy_x_zz, g_y_xy_x_xx, g_y_xy_x_xy, g_y_xy_x_xz, g_y_xy_x_yy, g_y_xy_x_yz, g_y_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_xy_x_xx[i] = -g_y_xy_x_xx[i] + 2.0 * g_xxy_xy_x_xx[i] * a_exp;

        g_x_0_0_0_xy_xy_x_xy[i] = -g_y_xy_x_xy[i] + 2.0 * g_xxy_xy_x_xy[i] * a_exp;

        g_x_0_0_0_xy_xy_x_xz[i] = -g_y_xy_x_xz[i] + 2.0 * g_xxy_xy_x_xz[i] * a_exp;

        g_x_0_0_0_xy_xy_x_yy[i] = -g_y_xy_x_yy[i] + 2.0 * g_xxy_xy_x_yy[i] * a_exp;

        g_x_0_0_0_xy_xy_x_yz[i] = -g_y_xy_x_yz[i] + 2.0 * g_xxy_xy_x_yz[i] * a_exp;

        g_x_0_0_0_xy_xy_x_zz[i] = -g_y_xy_x_zz[i] + 2.0 * g_xxy_xy_x_zz[i] * a_exp;
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_x_0_0_0_xy_xy_y_xx, g_x_0_0_0_xy_xy_y_xy, g_x_0_0_0_xy_xy_y_xz, g_x_0_0_0_xy_xy_y_yy, g_x_0_0_0_xy_xy_y_yz, g_x_0_0_0_xy_xy_y_zz, g_xxy_xy_y_xx, g_xxy_xy_y_xy, g_xxy_xy_y_xz, g_xxy_xy_y_yy, g_xxy_xy_y_yz, g_xxy_xy_y_zz, g_y_xy_y_xx, g_y_xy_y_xy, g_y_xy_y_xz, g_y_xy_y_yy, g_y_xy_y_yz, g_y_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_xy_y_xx[i] = -g_y_xy_y_xx[i] + 2.0 * g_xxy_xy_y_xx[i] * a_exp;

        g_x_0_0_0_xy_xy_y_xy[i] = -g_y_xy_y_xy[i] + 2.0 * g_xxy_xy_y_xy[i] * a_exp;

        g_x_0_0_0_xy_xy_y_xz[i] = -g_y_xy_y_xz[i] + 2.0 * g_xxy_xy_y_xz[i] * a_exp;

        g_x_0_0_0_xy_xy_y_yy[i] = -g_y_xy_y_yy[i] + 2.0 * g_xxy_xy_y_yy[i] * a_exp;

        g_x_0_0_0_xy_xy_y_yz[i] = -g_y_xy_y_yz[i] + 2.0 * g_xxy_xy_y_yz[i] * a_exp;

        g_x_0_0_0_xy_xy_y_zz[i] = -g_y_xy_y_zz[i] + 2.0 * g_xxy_xy_y_zz[i] * a_exp;
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_x_0_0_0_xy_xy_z_xx, g_x_0_0_0_xy_xy_z_xy, g_x_0_0_0_xy_xy_z_xz, g_x_0_0_0_xy_xy_z_yy, g_x_0_0_0_xy_xy_z_yz, g_x_0_0_0_xy_xy_z_zz, g_xxy_xy_z_xx, g_xxy_xy_z_xy, g_xxy_xy_z_xz, g_xxy_xy_z_yy, g_xxy_xy_z_yz, g_xxy_xy_z_zz, g_y_xy_z_xx, g_y_xy_z_xy, g_y_xy_z_xz, g_y_xy_z_yy, g_y_xy_z_yz, g_y_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_xy_z_xx[i] = -g_y_xy_z_xx[i] + 2.0 * g_xxy_xy_z_xx[i] * a_exp;

        g_x_0_0_0_xy_xy_z_xy[i] = -g_y_xy_z_xy[i] + 2.0 * g_xxy_xy_z_xy[i] * a_exp;

        g_x_0_0_0_xy_xy_z_xz[i] = -g_y_xy_z_xz[i] + 2.0 * g_xxy_xy_z_xz[i] * a_exp;

        g_x_0_0_0_xy_xy_z_yy[i] = -g_y_xy_z_yy[i] + 2.0 * g_xxy_xy_z_yy[i] * a_exp;

        g_x_0_0_0_xy_xy_z_yz[i] = -g_y_xy_z_yz[i] + 2.0 * g_xxy_xy_z_yz[i] * a_exp;

        g_x_0_0_0_xy_xy_z_zz[i] = -g_y_xy_z_zz[i] + 2.0 * g_xxy_xy_z_zz[i] * a_exp;
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_x_0_0_0_xy_xz_x_xx, g_x_0_0_0_xy_xz_x_xy, g_x_0_0_0_xy_xz_x_xz, g_x_0_0_0_xy_xz_x_yy, g_x_0_0_0_xy_xz_x_yz, g_x_0_0_0_xy_xz_x_zz, g_xxy_xz_x_xx, g_xxy_xz_x_xy, g_xxy_xz_x_xz, g_xxy_xz_x_yy, g_xxy_xz_x_yz, g_xxy_xz_x_zz, g_y_xz_x_xx, g_y_xz_x_xy, g_y_xz_x_xz, g_y_xz_x_yy, g_y_xz_x_yz, g_y_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_xz_x_xx[i] = -g_y_xz_x_xx[i] + 2.0 * g_xxy_xz_x_xx[i] * a_exp;

        g_x_0_0_0_xy_xz_x_xy[i] = -g_y_xz_x_xy[i] + 2.0 * g_xxy_xz_x_xy[i] * a_exp;

        g_x_0_0_0_xy_xz_x_xz[i] = -g_y_xz_x_xz[i] + 2.0 * g_xxy_xz_x_xz[i] * a_exp;

        g_x_0_0_0_xy_xz_x_yy[i] = -g_y_xz_x_yy[i] + 2.0 * g_xxy_xz_x_yy[i] * a_exp;

        g_x_0_0_0_xy_xz_x_yz[i] = -g_y_xz_x_yz[i] + 2.0 * g_xxy_xz_x_yz[i] * a_exp;

        g_x_0_0_0_xy_xz_x_zz[i] = -g_y_xz_x_zz[i] + 2.0 * g_xxy_xz_x_zz[i] * a_exp;
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_x_0_0_0_xy_xz_y_xx, g_x_0_0_0_xy_xz_y_xy, g_x_0_0_0_xy_xz_y_xz, g_x_0_0_0_xy_xz_y_yy, g_x_0_0_0_xy_xz_y_yz, g_x_0_0_0_xy_xz_y_zz, g_xxy_xz_y_xx, g_xxy_xz_y_xy, g_xxy_xz_y_xz, g_xxy_xz_y_yy, g_xxy_xz_y_yz, g_xxy_xz_y_zz, g_y_xz_y_xx, g_y_xz_y_xy, g_y_xz_y_xz, g_y_xz_y_yy, g_y_xz_y_yz, g_y_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_xz_y_xx[i] = -g_y_xz_y_xx[i] + 2.0 * g_xxy_xz_y_xx[i] * a_exp;

        g_x_0_0_0_xy_xz_y_xy[i] = -g_y_xz_y_xy[i] + 2.0 * g_xxy_xz_y_xy[i] * a_exp;

        g_x_0_0_0_xy_xz_y_xz[i] = -g_y_xz_y_xz[i] + 2.0 * g_xxy_xz_y_xz[i] * a_exp;

        g_x_0_0_0_xy_xz_y_yy[i] = -g_y_xz_y_yy[i] + 2.0 * g_xxy_xz_y_yy[i] * a_exp;

        g_x_0_0_0_xy_xz_y_yz[i] = -g_y_xz_y_yz[i] + 2.0 * g_xxy_xz_y_yz[i] * a_exp;

        g_x_0_0_0_xy_xz_y_zz[i] = -g_y_xz_y_zz[i] + 2.0 * g_xxy_xz_y_zz[i] * a_exp;
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_x_0_0_0_xy_xz_z_xx, g_x_0_0_0_xy_xz_z_xy, g_x_0_0_0_xy_xz_z_xz, g_x_0_0_0_xy_xz_z_yy, g_x_0_0_0_xy_xz_z_yz, g_x_0_0_0_xy_xz_z_zz, g_xxy_xz_z_xx, g_xxy_xz_z_xy, g_xxy_xz_z_xz, g_xxy_xz_z_yy, g_xxy_xz_z_yz, g_xxy_xz_z_zz, g_y_xz_z_xx, g_y_xz_z_xy, g_y_xz_z_xz, g_y_xz_z_yy, g_y_xz_z_yz, g_y_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_xz_z_xx[i] = -g_y_xz_z_xx[i] + 2.0 * g_xxy_xz_z_xx[i] * a_exp;

        g_x_0_0_0_xy_xz_z_xy[i] = -g_y_xz_z_xy[i] + 2.0 * g_xxy_xz_z_xy[i] * a_exp;

        g_x_0_0_0_xy_xz_z_xz[i] = -g_y_xz_z_xz[i] + 2.0 * g_xxy_xz_z_xz[i] * a_exp;

        g_x_0_0_0_xy_xz_z_yy[i] = -g_y_xz_z_yy[i] + 2.0 * g_xxy_xz_z_yy[i] * a_exp;

        g_x_0_0_0_xy_xz_z_yz[i] = -g_y_xz_z_yz[i] + 2.0 * g_xxy_xz_z_yz[i] * a_exp;

        g_x_0_0_0_xy_xz_z_zz[i] = -g_y_xz_z_zz[i] + 2.0 * g_xxy_xz_z_zz[i] * a_exp;
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_x_0_0_0_xy_yy_x_xx, g_x_0_0_0_xy_yy_x_xy, g_x_0_0_0_xy_yy_x_xz, g_x_0_0_0_xy_yy_x_yy, g_x_0_0_0_xy_yy_x_yz, g_x_0_0_0_xy_yy_x_zz, g_xxy_yy_x_xx, g_xxy_yy_x_xy, g_xxy_yy_x_xz, g_xxy_yy_x_yy, g_xxy_yy_x_yz, g_xxy_yy_x_zz, g_y_yy_x_xx, g_y_yy_x_xy, g_y_yy_x_xz, g_y_yy_x_yy, g_y_yy_x_yz, g_y_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_yy_x_xx[i] = -g_y_yy_x_xx[i] + 2.0 * g_xxy_yy_x_xx[i] * a_exp;

        g_x_0_0_0_xy_yy_x_xy[i] = -g_y_yy_x_xy[i] + 2.0 * g_xxy_yy_x_xy[i] * a_exp;

        g_x_0_0_0_xy_yy_x_xz[i] = -g_y_yy_x_xz[i] + 2.0 * g_xxy_yy_x_xz[i] * a_exp;

        g_x_0_0_0_xy_yy_x_yy[i] = -g_y_yy_x_yy[i] + 2.0 * g_xxy_yy_x_yy[i] * a_exp;

        g_x_0_0_0_xy_yy_x_yz[i] = -g_y_yy_x_yz[i] + 2.0 * g_xxy_yy_x_yz[i] * a_exp;

        g_x_0_0_0_xy_yy_x_zz[i] = -g_y_yy_x_zz[i] + 2.0 * g_xxy_yy_x_zz[i] * a_exp;
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_x_0_0_0_xy_yy_y_xx, g_x_0_0_0_xy_yy_y_xy, g_x_0_0_0_xy_yy_y_xz, g_x_0_0_0_xy_yy_y_yy, g_x_0_0_0_xy_yy_y_yz, g_x_0_0_0_xy_yy_y_zz, g_xxy_yy_y_xx, g_xxy_yy_y_xy, g_xxy_yy_y_xz, g_xxy_yy_y_yy, g_xxy_yy_y_yz, g_xxy_yy_y_zz, g_y_yy_y_xx, g_y_yy_y_xy, g_y_yy_y_xz, g_y_yy_y_yy, g_y_yy_y_yz, g_y_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_yy_y_xx[i] = -g_y_yy_y_xx[i] + 2.0 * g_xxy_yy_y_xx[i] * a_exp;

        g_x_0_0_0_xy_yy_y_xy[i] = -g_y_yy_y_xy[i] + 2.0 * g_xxy_yy_y_xy[i] * a_exp;

        g_x_0_0_0_xy_yy_y_xz[i] = -g_y_yy_y_xz[i] + 2.0 * g_xxy_yy_y_xz[i] * a_exp;

        g_x_0_0_0_xy_yy_y_yy[i] = -g_y_yy_y_yy[i] + 2.0 * g_xxy_yy_y_yy[i] * a_exp;

        g_x_0_0_0_xy_yy_y_yz[i] = -g_y_yy_y_yz[i] + 2.0 * g_xxy_yy_y_yz[i] * a_exp;

        g_x_0_0_0_xy_yy_y_zz[i] = -g_y_yy_y_zz[i] + 2.0 * g_xxy_yy_y_zz[i] * a_exp;
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_x_0_0_0_xy_yy_z_xx, g_x_0_0_0_xy_yy_z_xy, g_x_0_0_0_xy_yy_z_xz, g_x_0_0_0_xy_yy_z_yy, g_x_0_0_0_xy_yy_z_yz, g_x_0_0_0_xy_yy_z_zz, g_xxy_yy_z_xx, g_xxy_yy_z_xy, g_xxy_yy_z_xz, g_xxy_yy_z_yy, g_xxy_yy_z_yz, g_xxy_yy_z_zz, g_y_yy_z_xx, g_y_yy_z_xy, g_y_yy_z_xz, g_y_yy_z_yy, g_y_yy_z_yz, g_y_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_yy_z_xx[i] = -g_y_yy_z_xx[i] + 2.0 * g_xxy_yy_z_xx[i] * a_exp;

        g_x_0_0_0_xy_yy_z_xy[i] = -g_y_yy_z_xy[i] + 2.0 * g_xxy_yy_z_xy[i] * a_exp;

        g_x_0_0_0_xy_yy_z_xz[i] = -g_y_yy_z_xz[i] + 2.0 * g_xxy_yy_z_xz[i] * a_exp;

        g_x_0_0_0_xy_yy_z_yy[i] = -g_y_yy_z_yy[i] + 2.0 * g_xxy_yy_z_yy[i] * a_exp;

        g_x_0_0_0_xy_yy_z_yz[i] = -g_y_yy_z_yz[i] + 2.0 * g_xxy_yy_z_yz[i] * a_exp;

        g_x_0_0_0_xy_yy_z_zz[i] = -g_y_yy_z_zz[i] + 2.0 * g_xxy_yy_z_zz[i] * a_exp;
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_x_0_0_0_xy_yz_x_xx, g_x_0_0_0_xy_yz_x_xy, g_x_0_0_0_xy_yz_x_xz, g_x_0_0_0_xy_yz_x_yy, g_x_0_0_0_xy_yz_x_yz, g_x_0_0_0_xy_yz_x_zz, g_xxy_yz_x_xx, g_xxy_yz_x_xy, g_xxy_yz_x_xz, g_xxy_yz_x_yy, g_xxy_yz_x_yz, g_xxy_yz_x_zz, g_y_yz_x_xx, g_y_yz_x_xy, g_y_yz_x_xz, g_y_yz_x_yy, g_y_yz_x_yz, g_y_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_yz_x_xx[i] = -g_y_yz_x_xx[i] + 2.0 * g_xxy_yz_x_xx[i] * a_exp;

        g_x_0_0_0_xy_yz_x_xy[i] = -g_y_yz_x_xy[i] + 2.0 * g_xxy_yz_x_xy[i] * a_exp;

        g_x_0_0_0_xy_yz_x_xz[i] = -g_y_yz_x_xz[i] + 2.0 * g_xxy_yz_x_xz[i] * a_exp;

        g_x_0_0_0_xy_yz_x_yy[i] = -g_y_yz_x_yy[i] + 2.0 * g_xxy_yz_x_yy[i] * a_exp;

        g_x_0_0_0_xy_yz_x_yz[i] = -g_y_yz_x_yz[i] + 2.0 * g_xxy_yz_x_yz[i] * a_exp;

        g_x_0_0_0_xy_yz_x_zz[i] = -g_y_yz_x_zz[i] + 2.0 * g_xxy_yz_x_zz[i] * a_exp;
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_x_0_0_0_xy_yz_y_xx, g_x_0_0_0_xy_yz_y_xy, g_x_0_0_0_xy_yz_y_xz, g_x_0_0_0_xy_yz_y_yy, g_x_0_0_0_xy_yz_y_yz, g_x_0_0_0_xy_yz_y_zz, g_xxy_yz_y_xx, g_xxy_yz_y_xy, g_xxy_yz_y_xz, g_xxy_yz_y_yy, g_xxy_yz_y_yz, g_xxy_yz_y_zz, g_y_yz_y_xx, g_y_yz_y_xy, g_y_yz_y_xz, g_y_yz_y_yy, g_y_yz_y_yz, g_y_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_yz_y_xx[i] = -g_y_yz_y_xx[i] + 2.0 * g_xxy_yz_y_xx[i] * a_exp;

        g_x_0_0_0_xy_yz_y_xy[i] = -g_y_yz_y_xy[i] + 2.0 * g_xxy_yz_y_xy[i] * a_exp;

        g_x_0_0_0_xy_yz_y_xz[i] = -g_y_yz_y_xz[i] + 2.0 * g_xxy_yz_y_xz[i] * a_exp;

        g_x_0_0_0_xy_yz_y_yy[i] = -g_y_yz_y_yy[i] + 2.0 * g_xxy_yz_y_yy[i] * a_exp;

        g_x_0_0_0_xy_yz_y_yz[i] = -g_y_yz_y_yz[i] + 2.0 * g_xxy_yz_y_yz[i] * a_exp;

        g_x_0_0_0_xy_yz_y_zz[i] = -g_y_yz_y_zz[i] + 2.0 * g_xxy_yz_y_zz[i] * a_exp;
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_x_0_0_0_xy_yz_z_xx, g_x_0_0_0_xy_yz_z_xy, g_x_0_0_0_xy_yz_z_xz, g_x_0_0_0_xy_yz_z_yy, g_x_0_0_0_xy_yz_z_yz, g_x_0_0_0_xy_yz_z_zz, g_xxy_yz_z_xx, g_xxy_yz_z_xy, g_xxy_yz_z_xz, g_xxy_yz_z_yy, g_xxy_yz_z_yz, g_xxy_yz_z_zz, g_y_yz_z_xx, g_y_yz_z_xy, g_y_yz_z_xz, g_y_yz_z_yy, g_y_yz_z_yz, g_y_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_yz_z_xx[i] = -g_y_yz_z_xx[i] + 2.0 * g_xxy_yz_z_xx[i] * a_exp;

        g_x_0_0_0_xy_yz_z_xy[i] = -g_y_yz_z_xy[i] + 2.0 * g_xxy_yz_z_xy[i] * a_exp;

        g_x_0_0_0_xy_yz_z_xz[i] = -g_y_yz_z_xz[i] + 2.0 * g_xxy_yz_z_xz[i] * a_exp;

        g_x_0_0_0_xy_yz_z_yy[i] = -g_y_yz_z_yy[i] + 2.0 * g_xxy_yz_z_yy[i] * a_exp;

        g_x_0_0_0_xy_yz_z_yz[i] = -g_y_yz_z_yz[i] + 2.0 * g_xxy_yz_z_yz[i] * a_exp;

        g_x_0_0_0_xy_yz_z_zz[i] = -g_y_yz_z_zz[i] + 2.0 * g_xxy_yz_z_zz[i] * a_exp;
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_x_0_0_0_xy_zz_x_xx, g_x_0_0_0_xy_zz_x_xy, g_x_0_0_0_xy_zz_x_xz, g_x_0_0_0_xy_zz_x_yy, g_x_0_0_0_xy_zz_x_yz, g_x_0_0_0_xy_zz_x_zz, g_xxy_zz_x_xx, g_xxy_zz_x_xy, g_xxy_zz_x_xz, g_xxy_zz_x_yy, g_xxy_zz_x_yz, g_xxy_zz_x_zz, g_y_zz_x_xx, g_y_zz_x_xy, g_y_zz_x_xz, g_y_zz_x_yy, g_y_zz_x_yz, g_y_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_zz_x_xx[i] = -g_y_zz_x_xx[i] + 2.0 * g_xxy_zz_x_xx[i] * a_exp;

        g_x_0_0_0_xy_zz_x_xy[i] = -g_y_zz_x_xy[i] + 2.0 * g_xxy_zz_x_xy[i] * a_exp;

        g_x_0_0_0_xy_zz_x_xz[i] = -g_y_zz_x_xz[i] + 2.0 * g_xxy_zz_x_xz[i] * a_exp;

        g_x_0_0_0_xy_zz_x_yy[i] = -g_y_zz_x_yy[i] + 2.0 * g_xxy_zz_x_yy[i] * a_exp;

        g_x_0_0_0_xy_zz_x_yz[i] = -g_y_zz_x_yz[i] + 2.0 * g_xxy_zz_x_yz[i] * a_exp;

        g_x_0_0_0_xy_zz_x_zz[i] = -g_y_zz_x_zz[i] + 2.0 * g_xxy_zz_x_zz[i] * a_exp;
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_x_0_0_0_xy_zz_y_xx, g_x_0_0_0_xy_zz_y_xy, g_x_0_0_0_xy_zz_y_xz, g_x_0_0_0_xy_zz_y_yy, g_x_0_0_0_xy_zz_y_yz, g_x_0_0_0_xy_zz_y_zz, g_xxy_zz_y_xx, g_xxy_zz_y_xy, g_xxy_zz_y_xz, g_xxy_zz_y_yy, g_xxy_zz_y_yz, g_xxy_zz_y_zz, g_y_zz_y_xx, g_y_zz_y_xy, g_y_zz_y_xz, g_y_zz_y_yy, g_y_zz_y_yz, g_y_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_zz_y_xx[i] = -g_y_zz_y_xx[i] + 2.0 * g_xxy_zz_y_xx[i] * a_exp;

        g_x_0_0_0_xy_zz_y_xy[i] = -g_y_zz_y_xy[i] + 2.0 * g_xxy_zz_y_xy[i] * a_exp;

        g_x_0_0_0_xy_zz_y_xz[i] = -g_y_zz_y_xz[i] + 2.0 * g_xxy_zz_y_xz[i] * a_exp;

        g_x_0_0_0_xy_zz_y_yy[i] = -g_y_zz_y_yy[i] + 2.0 * g_xxy_zz_y_yy[i] * a_exp;

        g_x_0_0_0_xy_zz_y_yz[i] = -g_y_zz_y_yz[i] + 2.0 * g_xxy_zz_y_yz[i] * a_exp;

        g_x_0_0_0_xy_zz_y_zz[i] = -g_y_zz_y_zz[i] + 2.0 * g_xxy_zz_y_zz[i] * a_exp;
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_x_0_0_0_xy_zz_z_xx, g_x_0_0_0_xy_zz_z_xy, g_x_0_0_0_xy_zz_z_xz, g_x_0_0_0_xy_zz_z_yy, g_x_0_0_0_xy_zz_z_yz, g_x_0_0_0_xy_zz_z_zz, g_xxy_zz_z_xx, g_xxy_zz_z_xy, g_xxy_zz_z_xz, g_xxy_zz_z_yy, g_xxy_zz_z_yz, g_xxy_zz_z_zz, g_y_zz_z_xx, g_y_zz_z_xy, g_y_zz_z_xz, g_y_zz_z_yy, g_y_zz_z_yz, g_y_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_zz_z_xx[i] = -g_y_zz_z_xx[i] + 2.0 * g_xxy_zz_z_xx[i] * a_exp;

        g_x_0_0_0_xy_zz_z_xy[i] = -g_y_zz_z_xy[i] + 2.0 * g_xxy_zz_z_xy[i] * a_exp;

        g_x_0_0_0_xy_zz_z_xz[i] = -g_y_zz_z_xz[i] + 2.0 * g_xxy_zz_z_xz[i] * a_exp;

        g_x_0_0_0_xy_zz_z_yy[i] = -g_y_zz_z_yy[i] + 2.0 * g_xxy_zz_z_yy[i] * a_exp;

        g_x_0_0_0_xy_zz_z_yz[i] = -g_y_zz_z_yz[i] + 2.0 * g_xxy_zz_z_yz[i] * a_exp;

        g_x_0_0_0_xy_zz_z_zz[i] = -g_y_zz_z_zz[i] + 2.0 * g_xxy_zz_z_zz[i] * a_exp;
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_x_0_0_0_xz_xx_x_xx, g_x_0_0_0_xz_xx_x_xy, g_x_0_0_0_xz_xx_x_xz, g_x_0_0_0_xz_xx_x_yy, g_x_0_0_0_xz_xx_x_yz, g_x_0_0_0_xz_xx_x_zz, g_xxz_xx_x_xx, g_xxz_xx_x_xy, g_xxz_xx_x_xz, g_xxz_xx_x_yy, g_xxz_xx_x_yz, g_xxz_xx_x_zz, g_z_xx_x_xx, g_z_xx_x_xy, g_z_xx_x_xz, g_z_xx_x_yy, g_z_xx_x_yz, g_z_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_xx_x_xx[i] = -g_z_xx_x_xx[i] + 2.0 * g_xxz_xx_x_xx[i] * a_exp;

        g_x_0_0_0_xz_xx_x_xy[i] = -g_z_xx_x_xy[i] + 2.0 * g_xxz_xx_x_xy[i] * a_exp;

        g_x_0_0_0_xz_xx_x_xz[i] = -g_z_xx_x_xz[i] + 2.0 * g_xxz_xx_x_xz[i] * a_exp;

        g_x_0_0_0_xz_xx_x_yy[i] = -g_z_xx_x_yy[i] + 2.0 * g_xxz_xx_x_yy[i] * a_exp;

        g_x_0_0_0_xz_xx_x_yz[i] = -g_z_xx_x_yz[i] + 2.0 * g_xxz_xx_x_yz[i] * a_exp;

        g_x_0_0_0_xz_xx_x_zz[i] = -g_z_xx_x_zz[i] + 2.0 * g_xxz_xx_x_zz[i] * a_exp;
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_x_0_0_0_xz_xx_y_xx, g_x_0_0_0_xz_xx_y_xy, g_x_0_0_0_xz_xx_y_xz, g_x_0_0_0_xz_xx_y_yy, g_x_0_0_0_xz_xx_y_yz, g_x_0_0_0_xz_xx_y_zz, g_xxz_xx_y_xx, g_xxz_xx_y_xy, g_xxz_xx_y_xz, g_xxz_xx_y_yy, g_xxz_xx_y_yz, g_xxz_xx_y_zz, g_z_xx_y_xx, g_z_xx_y_xy, g_z_xx_y_xz, g_z_xx_y_yy, g_z_xx_y_yz, g_z_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_xx_y_xx[i] = -g_z_xx_y_xx[i] + 2.0 * g_xxz_xx_y_xx[i] * a_exp;

        g_x_0_0_0_xz_xx_y_xy[i] = -g_z_xx_y_xy[i] + 2.0 * g_xxz_xx_y_xy[i] * a_exp;

        g_x_0_0_0_xz_xx_y_xz[i] = -g_z_xx_y_xz[i] + 2.0 * g_xxz_xx_y_xz[i] * a_exp;

        g_x_0_0_0_xz_xx_y_yy[i] = -g_z_xx_y_yy[i] + 2.0 * g_xxz_xx_y_yy[i] * a_exp;

        g_x_0_0_0_xz_xx_y_yz[i] = -g_z_xx_y_yz[i] + 2.0 * g_xxz_xx_y_yz[i] * a_exp;

        g_x_0_0_0_xz_xx_y_zz[i] = -g_z_xx_y_zz[i] + 2.0 * g_xxz_xx_y_zz[i] * a_exp;
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_x_0_0_0_xz_xx_z_xx, g_x_0_0_0_xz_xx_z_xy, g_x_0_0_0_xz_xx_z_xz, g_x_0_0_0_xz_xx_z_yy, g_x_0_0_0_xz_xx_z_yz, g_x_0_0_0_xz_xx_z_zz, g_xxz_xx_z_xx, g_xxz_xx_z_xy, g_xxz_xx_z_xz, g_xxz_xx_z_yy, g_xxz_xx_z_yz, g_xxz_xx_z_zz, g_z_xx_z_xx, g_z_xx_z_xy, g_z_xx_z_xz, g_z_xx_z_yy, g_z_xx_z_yz, g_z_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_xx_z_xx[i] = -g_z_xx_z_xx[i] + 2.0 * g_xxz_xx_z_xx[i] * a_exp;

        g_x_0_0_0_xz_xx_z_xy[i] = -g_z_xx_z_xy[i] + 2.0 * g_xxz_xx_z_xy[i] * a_exp;

        g_x_0_0_0_xz_xx_z_xz[i] = -g_z_xx_z_xz[i] + 2.0 * g_xxz_xx_z_xz[i] * a_exp;

        g_x_0_0_0_xz_xx_z_yy[i] = -g_z_xx_z_yy[i] + 2.0 * g_xxz_xx_z_yy[i] * a_exp;

        g_x_0_0_0_xz_xx_z_yz[i] = -g_z_xx_z_yz[i] + 2.0 * g_xxz_xx_z_yz[i] * a_exp;

        g_x_0_0_0_xz_xx_z_zz[i] = -g_z_xx_z_zz[i] + 2.0 * g_xxz_xx_z_zz[i] * a_exp;
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_x_0_0_0_xz_xy_x_xx, g_x_0_0_0_xz_xy_x_xy, g_x_0_0_0_xz_xy_x_xz, g_x_0_0_0_xz_xy_x_yy, g_x_0_0_0_xz_xy_x_yz, g_x_0_0_0_xz_xy_x_zz, g_xxz_xy_x_xx, g_xxz_xy_x_xy, g_xxz_xy_x_xz, g_xxz_xy_x_yy, g_xxz_xy_x_yz, g_xxz_xy_x_zz, g_z_xy_x_xx, g_z_xy_x_xy, g_z_xy_x_xz, g_z_xy_x_yy, g_z_xy_x_yz, g_z_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_xy_x_xx[i] = -g_z_xy_x_xx[i] + 2.0 * g_xxz_xy_x_xx[i] * a_exp;

        g_x_0_0_0_xz_xy_x_xy[i] = -g_z_xy_x_xy[i] + 2.0 * g_xxz_xy_x_xy[i] * a_exp;

        g_x_0_0_0_xz_xy_x_xz[i] = -g_z_xy_x_xz[i] + 2.0 * g_xxz_xy_x_xz[i] * a_exp;

        g_x_0_0_0_xz_xy_x_yy[i] = -g_z_xy_x_yy[i] + 2.0 * g_xxz_xy_x_yy[i] * a_exp;

        g_x_0_0_0_xz_xy_x_yz[i] = -g_z_xy_x_yz[i] + 2.0 * g_xxz_xy_x_yz[i] * a_exp;

        g_x_0_0_0_xz_xy_x_zz[i] = -g_z_xy_x_zz[i] + 2.0 * g_xxz_xy_x_zz[i] * a_exp;
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_x_0_0_0_xz_xy_y_xx, g_x_0_0_0_xz_xy_y_xy, g_x_0_0_0_xz_xy_y_xz, g_x_0_0_0_xz_xy_y_yy, g_x_0_0_0_xz_xy_y_yz, g_x_0_0_0_xz_xy_y_zz, g_xxz_xy_y_xx, g_xxz_xy_y_xy, g_xxz_xy_y_xz, g_xxz_xy_y_yy, g_xxz_xy_y_yz, g_xxz_xy_y_zz, g_z_xy_y_xx, g_z_xy_y_xy, g_z_xy_y_xz, g_z_xy_y_yy, g_z_xy_y_yz, g_z_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_xy_y_xx[i] = -g_z_xy_y_xx[i] + 2.0 * g_xxz_xy_y_xx[i] * a_exp;

        g_x_0_0_0_xz_xy_y_xy[i] = -g_z_xy_y_xy[i] + 2.0 * g_xxz_xy_y_xy[i] * a_exp;

        g_x_0_0_0_xz_xy_y_xz[i] = -g_z_xy_y_xz[i] + 2.0 * g_xxz_xy_y_xz[i] * a_exp;

        g_x_0_0_0_xz_xy_y_yy[i] = -g_z_xy_y_yy[i] + 2.0 * g_xxz_xy_y_yy[i] * a_exp;

        g_x_0_0_0_xz_xy_y_yz[i] = -g_z_xy_y_yz[i] + 2.0 * g_xxz_xy_y_yz[i] * a_exp;

        g_x_0_0_0_xz_xy_y_zz[i] = -g_z_xy_y_zz[i] + 2.0 * g_xxz_xy_y_zz[i] * a_exp;
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_x_0_0_0_xz_xy_z_xx, g_x_0_0_0_xz_xy_z_xy, g_x_0_0_0_xz_xy_z_xz, g_x_0_0_0_xz_xy_z_yy, g_x_0_0_0_xz_xy_z_yz, g_x_0_0_0_xz_xy_z_zz, g_xxz_xy_z_xx, g_xxz_xy_z_xy, g_xxz_xy_z_xz, g_xxz_xy_z_yy, g_xxz_xy_z_yz, g_xxz_xy_z_zz, g_z_xy_z_xx, g_z_xy_z_xy, g_z_xy_z_xz, g_z_xy_z_yy, g_z_xy_z_yz, g_z_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_xy_z_xx[i] = -g_z_xy_z_xx[i] + 2.0 * g_xxz_xy_z_xx[i] * a_exp;

        g_x_0_0_0_xz_xy_z_xy[i] = -g_z_xy_z_xy[i] + 2.0 * g_xxz_xy_z_xy[i] * a_exp;

        g_x_0_0_0_xz_xy_z_xz[i] = -g_z_xy_z_xz[i] + 2.0 * g_xxz_xy_z_xz[i] * a_exp;

        g_x_0_0_0_xz_xy_z_yy[i] = -g_z_xy_z_yy[i] + 2.0 * g_xxz_xy_z_yy[i] * a_exp;

        g_x_0_0_0_xz_xy_z_yz[i] = -g_z_xy_z_yz[i] + 2.0 * g_xxz_xy_z_yz[i] * a_exp;

        g_x_0_0_0_xz_xy_z_zz[i] = -g_z_xy_z_zz[i] + 2.0 * g_xxz_xy_z_zz[i] * a_exp;
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_x_0_0_0_xz_xz_x_xx, g_x_0_0_0_xz_xz_x_xy, g_x_0_0_0_xz_xz_x_xz, g_x_0_0_0_xz_xz_x_yy, g_x_0_0_0_xz_xz_x_yz, g_x_0_0_0_xz_xz_x_zz, g_xxz_xz_x_xx, g_xxz_xz_x_xy, g_xxz_xz_x_xz, g_xxz_xz_x_yy, g_xxz_xz_x_yz, g_xxz_xz_x_zz, g_z_xz_x_xx, g_z_xz_x_xy, g_z_xz_x_xz, g_z_xz_x_yy, g_z_xz_x_yz, g_z_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_xz_x_xx[i] = -g_z_xz_x_xx[i] + 2.0 * g_xxz_xz_x_xx[i] * a_exp;

        g_x_0_0_0_xz_xz_x_xy[i] = -g_z_xz_x_xy[i] + 2.0 * g_xxz_xz_x_xy[i] * a_exp;

        g_x_0_0_0_xz_xz_x_xz[i] = -g_z_xz_x_xz[i] + 2.0 * g_xxz_xz_x_xz[i] * a_exp;

        g_x_0_0_0_xz_xz_x_yy[i] = -g_z_xz_x_yy[i] + 2.0 * g_xxz_xz_x_yy[i] * a_exp;

        g_x_0_0_0_xz_xz_x_yz[i] = -g_z_xz_x_yz[i] + 2.0 * g_xxz_xz_x_yz[i] * a_exp;

        g_x_0_0_0_xz_xz_x_zz[i] = -g_z_xz_x_zz[i] + 2.0 * g_xxz_xz_x_zz[i] * a_exp;
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_x_0_0_0_xz_xz_y_xx, g_x_0_0_0_xz_xz_y_xy, g_x_0_0_0_xz_xz_y_xz, g_x_0_0_0_xz_xz_y_yy, g_x_0_0_0_xz_xz_y_yz, g_x_0_0_0_xz_xz_y_zz, g_xxz_xz_y_xx, g_xxz_xz_y_xy, g_xxz_xz_y_xz, g_xxz_xz_y_yy, g_xxz_xz_y_yz, g_xxz_xz_y_zz, g_z_xz_y_xx, g_z_xz_y_xy, g_z_xz_y_xz, g_z_xz_y_yy, g_z_xz_y_yz, g_z_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_xz_y_xx[i] = -g_z_xz_y_xx[i] + 2.0 * g_xxz_xz_y_xx[i] * a_exp;

        g_x_0_0_0_xz_xz_y_xy[i] = -g_z_xz_y_xy[i] + 2.0 * g_xxz_xz_y_xy[i] * a_exp;

        g_x_0_0_0_xz_xz_y_xz[i] = -g_z_xz_y_xz[i] + 2.0 * g_xxz_xz_y_xz[i] * a_exp;

        g_x_0_0_0_xz_xz_y_yy[i] = -g_z_xz_y_yy[i] + 2.0 * g_xxz_xz_y_yy[i] * a_exp;

        g_x_0_0_0_xz_xz_y_yz[i] = -g_z_xz_y_yz[i] + 2.0 * g_xxz_xz_y_yz[i] * a_exp;

        g_x_0_0_0_xz_xz_y_zz[i] = -g_z_xz_y_zz[i] + 2.0 * g_xxz_xz_y_zz[i] * a_exp;
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_x_0_0_0_xz_xz_z_xx, g_x_0_0_0_xz_xz_z_xy, g_x_0_0_0_xz_xz_z_xz, g_x_0_0_0_xz_xz_z_yy, g_x_0_0_0_xz_xz_z_yz, g_x_0_0_0_xz_xz_z_zz, g_xxz_xz_z_xx, g_xxz_xz_z_xy, g_xxz_xz_z_xz, g_xxz_xz_z_yy, g_xxz_xz_z_yz, g_xxz_xz_z_zz, g_z_xz_z_xx, g_z_xz_z_xy, g_z_xz_z_xz, g_z_xz_z_yy, g_z_xz_z_yz, g_z_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_xz_z_xx[i] = -g_z_xz_z_xx[i] + 2.0 * g_xxz_xz_z_xx[i] * a_exp;

        g_x_0_0_0_xz_xz_z_xy[i] = -g_z_xz_z_xy[i] + 2.0 * g_xxz_xz_z_xy[i] * a_exp;

        g_x_0_0_0_xz_xz_z_xz[i] = -g_z_xz_z_xz[i] + 2.0 * g_xxz_xz_z_xz[i] * a_exp;

        g_x_0_0_0_xz_xz_z_yy[i] = -g_z_xz_z_yy[i] + 2.0 * g_xxz_xz_z_yy[i] * a_exp;

        g_x_0_0_0_xz_xz_z_yz[i] = -g_z_xz_z_yz[i] + 2.0 * g_xxz_xz_z_yz[i] * a_exp;

        g_x_0_0_0_xz_xz_z_zz[i] = -g_z_xz_z_zz[i] + 2.0 * g_xxz_xz_z_zz[i] * a_exp;
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_x_0_0_0_xz_yy_x_xx, g_x_0_0_0_xz_yy_x_xy, g_x_0_0_0_xz_yy_x_xz, g_x_0_0_0_xz_yy_x_yy, g_x_0_0_0_xz_yy_x_yz, g_x_0_0_0_xz_yy_x_zz, g_xxz_yy_x_xx, g_xxz_yy_x_xy, g_xxz_yy_x_xz, g_xxz_yy_x_yy, g_xxz_yy_x_yz, g_xxz_yy_x_zz, g_z_yy_x_xx, g_z_yy_x_xy, g_z_yy_x_xz, g_z_yy_x_yy, g_z_yy_x_yz, g_z_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_yy_x_xx[i] = -g_z_yy_x_xx[i] + 2.0 * g_xxz_yy_x_xx[i] * a_exp;

        g_x_0_0_0_xz_yy_x_xy[i] = -g_z_yy_x_xy[i] + 2.0 * g_xxz_yy_x_xy[i] * a_exp;

        g_x_0_0_0_xz_yy_x_xz[i] = -g_z_yy_x_xz[i] + 2.0 * g_xxz_yy_x_xz[i] * a_exp;

        g_x_0_0_0_xz_yy_x_yy[i] = -g_z_yy_x_yy[i] + 2.0 * g_xxz_yy_x_yy[i] * a_exp;

        g_x_0_0_0_xz_yy_x_yz[i] = -g_z_yy_x_yz[i] + 2.0 * g_xxz_yy_x_yz[i] * a_exp;

        g_x_0_0_0_xz_yy_x_zz[i] = -g_z_yy_x_zz[i] + 2.0 * g_xxz_yy_x_zz[i] * a_exp;
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_x_0_0_0_xz_yy_y_xx, g_x_0_0_0_xz_yy_y_xy, g_x_0_0_0_xz_yy_y_xz, g_x_0_0_0_xz_yy_y_yy, g_x_0_0_0_xz_yy_y_yz, g_x_0_0_0_xz_yy_y_zz, g_xxz_yy_y_xx, g_xxz_yy_y_xy, g_xxz_yy_y_xz, g_xxz_yy_y_yy, g_xxz_yy_y_yz, g_xxz_yy_y_zz, g_z_yy_y_xx, g_z_yy_y_xy, g_z_yy_y_xz, g_z_yy_y_yy, g_z_yy_y_yz, g_z_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_yy_y_xx[i] = -g_z_yy_y_xx[i] + 2.0 * g_xxz_yy_y_xx[i] * a_exp;

        g_x_0_0_0_xz_yy_y_xy[i] = -g_z_yy_y_xy[i] + 2.0 * g_xxz_yy_y_xy[i] * a_exp;

        g_x_0_0_0_xz_yy_y_xz[i] = -g_z_yy_y_xz[i] + 2.0 * g_xxz_yy_y_xz[i] * a_exp;

        g_x_0_0_0_xz_yy_y_yy[i] = -g_z_yy_y_yy[i] + 2.0 * g_xxz_yy_y_yy[i] * a_exp;

        g_x_0_0_0_xz_yy_y_yz[i] = -g_z_yy_y_yz[i] + 2.0 * g_xxz_yy_y_yz[i] * a_exp;

        g_x_0_0_0_xz_yy_y_zz[i] = -g_z_yy_y_zz[i] + 2.0 * g_xxz_yy_y_zz[i] * a_exp;
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_x_0_0_0_xz_yy_z_xx, g_x_0_0_0_xz_yy_z_xy, g_x_0_0_0_xz_yy_z_xz, g_x_0_0_0_xz_yy_z_yy, g_x_0_0_0_xz_yy_z_yz, g_x_0_0_0_xz_yy_z_zz, g_xxz_yy_z_xx, g_xxz_yy_z_xy, g_xxz_yy_z_xz, g_xxz_yy_z_yy, g_xxz_yy_z_yz, g_xxz_yy_z_zz, g_z_yy_z_xx, g_z_yy_z_xy, g_z_yy_z_xz, g_z_yy_z_yy, g_z_yy_z_yz, g_z_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_yy_z_xx[i] = -g_z_yy_z_xx[i] + 2.0 * g_xxz_yy_z_xx[i] * a_exp;

        g_x_0_0_0_xz_yy_z_xy[i] = -g_z_yy_z_xy[i] + 2.0 * g_xxz_yy_z_xy[i] * a_exp;

        g_x_0_0_0_xz_yy_z_xz[i] = -g_z_yy_z_xz[i] + 2.0 * g_xxz_yy_z_xz[i] * a_exp;

        g_x_0_0_0_xz_yy_z_yy[i] = -g_z_yy_z_yy[i] + 2.0 * g_xxz_yy_z_yy[i] * a_exp;

        g_x_0_0_0_xz_yy_z_yz[i] = -g_z_yy_z_yz[i] + 2.0 * g_xxz_yy_z_yz[i] * a_exp;

        g_x_0_0_0_xz_yy_z_zz[i] = -g_z_yy_z_zz[i] + 2.0 * g_xxz_yy_z_zz[i] * a_exp;
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_x_0_0_0_xz_yz_x_xx, g_x_0_0_0_xz_yz_x_xy, g_x_0_0_0_xz_yz_x_xz, g_x_0_0_0_xz_yz_x_yy, g_x_0_0_0_xz_yz_x_yz, g_x_0_0_0_xz_yz_x_zz, g_xxz_yz_x_xx, g_xxz_yz_x_xy, g_xxz_yz_x_xz, g_xxz_yz_x_yy, g_xxz_yz_x_yz, g_xxz_yz_x_zz, g_z_yz_x_xx, g_z_yz_x_xy, g_z_yz_x_xz, g_z_yz_x_yy, g_z_yz_x_yz, g_z_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_yz_x_xx[i] = -g_z_yz_x_xx[i] + 2.0 * g_xxz_yz_x_xx[i] * a_exp;

        g_x_0_0_0_xz_yz_x_xy[i] = -g_z_yz_x_xy[i] + 2.0 * g_xxz_yz_x_xy[i] * a_exp;

        g_x_0_0_0_xz_yz_x_xz[i] = -g_z_yz_x_xz[i] + 2.0 * g_xxz_yz_x_xz[i] * a_exp;

        g_x_0_0_0_xz_yz_x_yy[i] = -g_z_yz_x_yy[i] + 2.0 * g_xxz_yz_x_yy[i] * a_exp;

        g_x_0_0_0_xz_yz_x_yz[i] = -g_z_yz_x_yz[i] + 2.0 * g_xxz_yz_x_yz[i] * a_exp;

        g_x_0_0_0_xz_yz_x_zz[i] = -g_z_yz_x_zz[i] + 2.0 * g_xxz_yz_x_zz[i] * a_exp;
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_x_0_0_0_xz_yz_y_xx, g_x_0_0_0_xz_yz_y_xy, g_x_0_0_0_xz_yz_y_xz, g_x_0_0_0_xz_yz_y_yy, g_x_0_0_0_xz_yz_y_yz, g_x_0_0_0_xz_yz_y_zz, g_xxz_yz_y_xx, g_xxz_yz_y_xy, g_xxz_yz_y_xz, g_xxz_yz_y_yy, g_xxz_yz_y_yz, g_xxz_yz_y_zz, g_z_yz_y_xx, g_z_yz_y_xy, g_z_yz_y_xz, g_z_yz_y_yy, g_z_yz_y_yz, g_z_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_yz_y_xx[i] = -g_z_yz_y_xx[i] + 2.0 * g_xxz_yz_y_xx[i] * a_exp;

        g_x_0_0_0_xz_yz_y_xy[i] = -g_z_yz_y_xy[i] + 2.0 * g_xxz_yz_y_xy[i] * a_exp;

        g_x_0_0_0_xz_yz_y_xz[i] = -g_z_yz_y_xz[i] + 2.0 * g_xxz_yz_y_xz[i] * a_exp;

        g_x_0_0_0_xz_yz_y_yy[i] = -g_z_yz_y_yy[i] + 2.0 * g_xxz_yz_y_yy[i] * a_exp;

        g_x_0_0_0_xz_yz_y_yz[i] = -g_z_yz_y_yz[i] + 2.0 * g_xxz_yz_y_yz[i] * a_exp;

        g_x_0_0_0_xz_yz_y_zz[i] = -g_z_yz_y_zz[i] + 2.0 * g_xxz_yz_y_zz[i] * a_exp;
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_x_0_0_0_xz_yz_z_xx, g_x_0_0_0_xz_yz_z_xy, g_x_0_0_0_xz_yz_z_xz, g_x_0_0_0_xz_yz_z_yy, g_x_0_0_0_xz_yz_z_yz, g_x_0_0_0_xz_yz_z_zz, g_xxz_yz_z_xx, g_xxz_yz_z_xy, g_xxz_yz_z_xz, g_xxz_yz_z_yy, g_xxz_yz_z_yz, g_xxz_yz_z_zz, g_z_yz_z_xx, g_z_yz_z_xy, g_z_yz_z_xz, g_z_yz_z_yy, g_z_yz_z_yz, g_z_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_yz_z_xx[i] = -g_z_yz_z_xx[i] + 2.0 * g_xxz_yz_z_xx[i] * a_exp;

        g_x_0_0_0_xz_yz_z_xy[i] = -g_z_yz_z_xy[i] + 2.0 * g_xxz_yz_z_xy[i] * a_exp;

        g_x_0_0_0_xz_yz_z_xz[i] = -g_z_yz_z_xz[i] + 2.0 * g_xxz_yz_z_xz[i] * a_exp;

        g_x_0_0_0_xz_yz_z_yy[i] = -g_z_yz_z_yy[i] + 2.0 * g_xxz_yz_z_yy[i] * a_exp;

        g_x_0_0_0_xz_yz_z_yz[i] = -g_z_yz_z_yz[i] + 2.0 * g_xxz_yz_z_yz[i] * a_exp;

        g_x_0_0_0_xz_yz_z_zz[i] = -g_z_yz_z_zz[i] + 2.0 * g_xxz_yz_z_zz[i] * a_exp;
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_x_0_0_0_xz_zz_x_xx, g_x_0_0_0_xz_zz_x_xy, g_x_0_0_0_xz_zz_x_xz, g_x_0_0_0_xz_zz_x_yy, g_x_0_0_0_xz_zz_x_yz, g_x_0_0_0_xz_zz_x_zz, g_xxz_zz_x_xx, g_xxz_zz_x_xy, g_xxz_zz_x_xz, g_xxz_zz_x_yy, g_xxz_zz_x_yz, g_xxz_zz_x_zz, g_z_zz_x_xx, g_z_zz_x_xy, g_z_zz_x_xz, g_z_zz_x_yy, g_z_zz_x_yz, g_z_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_zz_x_xx[i] = -g_z_zz_x_xx[i] + 2.0 * g_xxz_zz_x_xx[i] * a_exp;

        g_x_0_0_0_xz_zz_x_xy[i] = -g_z_zz_x_xy[i] + 2.0 * g_xxz_zz_x_xy[i] * a_exp;

        g_x_0_0_0_xz_zz_x_xz[i] = -g_z_zz_x_xz[i] + 2.0 * g_xxz_zz_x_xz[i] * a_exp;

        g_x_0_0_0_xz_zz_x_yy[i] = -g_z_zz_x_yy[i] + 2.0 * g_xxz_zz_x_yy[i] * a_exp;

        g_x_0_0_0_xz_zz_x_yz[i] = -g_z_zz_x_yz[i] + 2.0 * g_xxz_zz_x_yz[i] * a_exp;

        g_x_0_0_0_xz_zz_x_zz[i] = -g_z_zz_x_zz[i] + 2.0 * g_xxz_zz_x_zz[i] * a_exp;
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_x_0_0_0_xz_zz_y_xx, g_x_0_0_0_xz_zz_y_xy, g_x_0_0_0_xz_zz_y_xz, g_x_0_0_0_xz_zz_y_yy, g_x_0_0_0_xz_zz_y_yz, g_x_0_0_0_xz_zz_y_zz, g_xxz_zz_y_xx, g_xxz_zz_y_xy, g_xxz_zz_y_xz, g_xxz_zz_y_yy, g_xxz_zz_y_yz, g_xxz_zz_y_zz, g_z_zz_y_xx, g_z_zz_y_xy, g_z_zz_y_xz, g_z_zz_y_yy, g_z_zz_y_yz, g_z_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_zz_y_xx[i] = -g_z_zz_y_xx[i] + 2.0 * g_xxz_zz_y_xx[i] * a_exp;

        g_x_0_0_0_xz_zz_y_xy[i] = -g_z_zz_y_xy[i] + 2.0 * g_xxz_zz_y_xy[i] * a_exp;

        g_x_0_0_0_xz_zz_y_xz[i] = -g_z_zz_y_xz[i] + 2.0 * g_xxz_zz_y_xz[i] * a_exp;

        g_x_0_0_0_xz_zz_y_yy[i] = -g_z_zz_y_yy[i] + 2.0 * g_xxz_zz_y_yy[i] * a_exp;

        g_x_0_0_0_xz_zz_y_yz[i] = -g_z_zz_y_yz[i] + 2.0 * g_xxz_zz_y_yz[i] * a_exp;

        g_x_0_0_0_xz_zz_y_zz[i] = -g_z_zz_y_zz[i] + 2.0 * g_xxz_zz_y_zz[i] * a_exp;
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_x_0_0_0_xz_zz_z_xx, g_x_0_0_0_xz_zz_z_xy, g_x_0_0_0_xz_zz_z_xz, g_x_0_0_0_xz_zz_z_yy, g_x_0_0_0_xz_zz_z_yz, g_x_0_0_0_xz_zz_z_zz, g_xxz_zz_z_xx, g_xxz_zz_z_xy, g_xxz_zz_z_xz, g_xxz_zz_z_yy, g_xxz_zz_z_yz, g_xxz_zz_z_zz, g_z_zz_z_xx, g_z_zz_z_xy, g_z_zz_z_xz, g_z_zz_z_yy, g_z_zz_z_yz, g_z_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_zz_z_xx[i] = -g_z_zz_z_xx[i] + 2.0 * g_xxz_zz_z_xx[i] * a_exp;

        g_x_0_0_0_xz_zz_z_xy[i] = -g_z_zz_z_xy[i] + 2.0 * g_xxz_zz_z_xy[i] * a_exp;

        g_x_0_0_0_xz_zz_z_xz[i] = -g_z_zz_z_xz[i] + 2.0 * g_xxz_zz_z_xz[i] * a_exp;

        g_x_0_0_0_xz_zz_z_yy[i] = -g_z_zz_z_yy[i] + 2.0 * g_xxz_zz_z_yy[i] * a_exp;

        g_x_0_0_0_xz_zz_z_yz[i] = -g_z_zz_z_yz[i] + 2.0 * g_xxz_zz_z_yz[i] * a_exp;

        g_x_0_0_0_xz_zz_z_zz[i] = -g_z_zz_z_zz[i] + 2.0 * g_xxz_zz_z_zz[i] * a_exp;
    }
    // integrals block (324-330)

    #pragma omp simd aligned(g_x_0_0_0_yy_xx_x_xx, g_x_0_0_0_yy_xx_x_xy, g_x_0_0_0_yy_xx_x_xz, g_x_0_0_0_yy_xx_x_yy, g_x_0_0_0_yy_xx_x_yz, g_x_0_0_0_yy_xx_x_zz, g_xyy_xx_x_xx, g_xyy_xx_x_xy, g_xyy_xx_x_xz, g_xyy_xx_x_yy, g_xyy_xx_x_yz, g_xyy_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_xx_x_xx[i] = 2.0 * g_xyy_xx_x_xx[i] * a_exp;

        g_x_0_0_0_yy_xx_x_xy[i] = 2.0 * g_xyy_xx_x_xy[i] * a_exp;

        g_x_0_0_0_yy_xx_x_xz[i] = 2.0 * g_xyy_xx_x_xz[i] * a_exp;

        g_x_0_0_0_yy_xx_x_yy[i] = 2.0 * g_xyy_xx_x_yy[i] * a_exp;

        g_x_0_0_0_yy_xx_x_yz[i] = 2.0 * g_xyy_xx_x_yz[i] * a_exp;

        g_x_0_0_0_yy_xx_x_zz[i] = 2.0 * g_xyy_xx_x_zz[i] * a_exp;
    }
    // integrals block (330-336)

    #pragma omp simd aligned(g_x_0_0_0_yy_xx_y_xx, g_x_0_0_0_yy_xx_y_xy, g_x_0_0_0_yy_xx_y_xz, g_x_0_0_0_yy_xx_y_yy, g_x_0_0_0_yy_xx_y_yz, g_x_0_0_0_yy_xx_y_zz, g_xyy_xx_y_xx, g_xyy_xx_y_xy, g_xyy_xx_y_xz, g_xyy_xx_y_yy, g_xyy_xx_y_yz, g_xyy_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_xx_y_xx[i] = 2.0 * g_xyy_xx_y_xx[i] * a_exp;

        g_x_0_0_0_yy_xx_y_xy[i] = 2.0 * g_xyy_xx_y_xy[i] * a_exp;

        g_x_0_0_0_yy_xx_y_xz[i] = 2.0 * g_xyy_xx_y_xz[i] * a_exp;

        g_x_0_0_0_yy_xx_y_yy[i] = 2.0 * g_xyy_xx_y_yy[i] * a_exp;

        g_x_0_0_0_yy_xx_y_yz[i] = 2.0 * g_xyy_xx_y_yz[i] * a_exp;

        g_x_0_0_0_yy_xx_y_zz[i] = 2.0 * g_xyy_xx_y_zz[i] * a_exp;
    }
    // integrals block (336-342)

    #pragma omp simd aligned(g_x_0_0_0_yy_xx_z_xx, g_x_0_0_0_yy_xx_z_xy, g_x_0_0_0_yy_xx_z_xz, g_x_0_0_0_yy_xx_z_yy, g_x_0_0_0_yy_xx_z_yz, g_x_0_0_0_yy_xx_z_zz, g_xyy_xx_z_xx, g_xyy_xx_z_xy, g_xyy_xx_z_xz, g_xyy_xx_z_yy, g_xyy_xx_z_yz, g_xyy_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_xx_z_xx[i] = 2.0 * g_xyy_xx_z_xx[i] * a_exp;

        g_x_0_0_0_yy_xx_z_xy[i] = 2.0 * g_xyy_xx_z_xy[i] * a_exp;

        g_x_0_0_0_yy_xx_z_xz[i] = 2.0 * g_xyy_xx_z_xz[i] * a_exp;

        g_x_0_0_0_yy_xx_z_yy[i] = 2.0 * g_xyy_xx_z_yy[i] * a_exp;

        g_x_0_0_0_yy_xx_z_yz[i] = 2.0 * g_xyy_xx_z_yz[i] * a_exp;

        g_x_0_0_0_yy_xx_z_zz[i] = 2.0 * g_xyy_xx_z_zz[i] * a_exp;
    }
    // integrals block (342-348)

    #pragma omp simd aligned(g_x_0_0_0_yy_xy_x_xx, g_x_0_0_0_yy_xy_x_xy, g_x_0_0_0_yy_xy_x_xz, g_x_0_0_0_yy_xy_x_yy, g_x_0_0_0_yy_xy_x_yz, g_x_0_0_0_yy_xy_x_zz, g_xyy_xy_x_xx, g_xyy_xy_x_xy, g_xyy_xy_x_xz, g_xyy_xy_x_yy, g_xyy_xy_x_yz, g_xyy_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_xy_x_xx[i] = 2.0 * g_xyy_xy_x_xx[i] * a_exp;

        g_x_0_0_0_yy_xy_x_xy[i] = 2.0 * g_xyy_xy_x_xy[i] * a_exp;

        g_x_0_0_0_yy_xy_x_xz[i] = 2.0 * g_xyy_xy_x_xz[i] * a_exp;

        g_x_0_0_0_yy_xy_x_yy[i] = 2.0 * g_xyy_xy_x_yy[i] * a_exp;

        g_x_0_0_0_yy_xy_x_yz[i] = 2.0 * g_xyy_xy_x_yz[i] * a_exp;

        g_x_0_0_0_yy_xy_x_zz[i] = 2.0 * g_xyy_xy_x_zz[i] * a_exp;
    }
    // integrals block (348-354)

    #pragma omp simd aligned(g_x_0_0_0_yy_xy_y_xx, g_x_0_0_0_yy_xy_y_xy, g_x_0_0_0_yy_xy_y_xz, g_x_0_0_0_yy_xy_y_yy, g_x_0_0_0_yy_xy_y_yz, g_x_0_0_0_yy_xy_y_zz, g_xyy_xy_y_xx, g_xyy_xy_y_xy, g_xyy_xy_y_xz, g_xyy_xy_y_yy, g_xyy_xy_y_yz, g_xyy_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_xy_y_xx[i] = 2.0 * g_xyy_xy_y_xx[i] * a_exp;

        g_x_0_0_0_yy_xy_y_xy[i] = 2.0 * g_xyy_xy_y_xy[i] * a_exp;

        g_x_0_0_0_yy_xy_y_xz[i] = 2.0 * g_xyy_xy_y_xz[i] * a_exp;

        g_x_0_0_0_yy_xy_y_yy[i] = 2.0 * g_xyy_xy_y_yy[i] * a_exp;

        g_x_0_0_0_yy_xy_y_yz[i] = 2.0 * g_xyy_xy_y_yz[i] * a_exp;

        g_x_0_0_0_yy_xy_y_zz[i] = 2.0 * g_xyy_xy_y_zz[i] * a_exp;
    }
    // integrals block (354-360)

    #pragma omp simd aligned(g_x_0_0_0_yy_xy_z_xx, g_x_0_0_0_yy_xy_z_xy, g_x_0_0_0_yy_xy_z_xz, g_x_0_0_0_yy_xy_z_yy, g_x_0_0_0_yy_xy_z_yz, g_x_0_0_0_yy_xy_z_zz, g_xyy_xy_z_xx, g_xyy_xy_z_xy, g_xyy_xy_z_xz, g_xyy_xy_z_yy, g_xyy_xy_z_yz, g_xyy_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_xy_z_xx[i] = 2.0 * g_xyy_xy_z_xx[i] * a_exp;

        g_x_0_0_0_yy_xy_z_xy[i] = 2.0 * g_xyy_xy_z_xy[i] * a_exp;

        g_x_0_0_0_yy_xy_z_xz[i] = 2.0 * g_xyy_xy_z_xz[i] * a_exp;

        g_x_0_0_0_yy_xy_z_yy[i] = 2.0 * g_xyy_xy_z_yy[i] * a_exp;

        g_x_0_0_0_yy_xy_z_yz[i] = 2.0 * g_xyy_xy_z_yz[i] * a_exp;

        g_x_0_0_0_yy_xy_z_zz[i] = 2.0 * g_xyy_xy_z_zz[i] * a_exp;
    }
    // integrals block (360-366)

    #pragma omp simd aligned(g_x_0_0_0_yy_xz_x_xx, g_x_0_0_0_yy_xz_x_xy, g_x_0_0_0_yy_xz_x_xz, g_x_0_0_0_yy_xz_x_yy, g_x_0_0_0_yy_xz_x_yz, g_x_0_0_0_yy_xz_x_zz, g_xyy_xz_x_xx, g_xyy_xz_x_xy, g_xyy_xz_x_xz, g_xyy_xz_x_yy, g_xyy_xz_x_yz, g_xyy_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_xz_x_xx[i] = 2.0 * g_xyy_xz_x_xx[i] * a_exp;

        g_x_0_0_0_yy_xz_x_xy[i] = 2.0 * g_xyy_xz_x_xy[i] * a_exp;

        g_x_0_0_0_yy_xz_x_xz[i] = 2.0 * g_xyy_xz_x_xz[i] * a_exp;

        g_x_0_0_0_yy_xz_x_yy[i] = 2.0 * g_xyy_xz_x_yy[i] * a_exp;

        g_x_0_0_0_yy_xz_x_yz[i] = 2.0 * g_xyy_xz_x_yz[i] * a_exp;

        g_x_0_0_0_yy_xz_x_zz[i] = 2.0 * g_xyy_xz_x_zz[i] * a_exp;
    }
    // integrals block (366-372)

    #pragma omp simd aligned(g_x_0_0_0_yy_xz_y_xx, g_x_0_0_0_yy_xz_y_xy, g_x_0_0_0_yy_xz_y_xz, g_x_0_0_0_yy_xz_y_yy, g_x_0_0_0_yy_xz_y_yz, g_x_0_0_0_yy_xz_y_zz, g_xyy_xz_y_xx, g_xyy_xz_y_xy, g_xyy_xz_y_xz, g_xyy_xz_y_yy, g_xyy_xz_y_yz, g_xyy_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_xz_y_xx[i] = 2.0 * g_xyy_xz_y_xx[i] * a_exp;

        g_x_0_0_0_yy_xz_y_xy[i] = 2.0 * g_xyy_xz_y_xy[i] * a_exp;

        g_x_0_0_0_yy_xz_y_xz[i] = 2.0 * g_xyy_xz_y_xz[i] * a_exp;

        g_x_0_0_0_yy_xz_y_yy[i] = 2.0 * g_xyy_xz_y_yy[i] * a_exp;

        g_x_0_0_0_yy_xz_y_yz[i] = 2.0 * g_xyy_xz_y_yz[i] * a_exp;

        g_x_0_0_0_yy_xz_y_zz[i] = 2.0 * g_xyy_xz_y_zz[i] * a_exp;
    }
    // integrals block (372-378)

    #pragma omp simd aligned(g_x_0_0_0_yy_xz_z_xx, g_x_0_0_0_yy_xz_z_xy, g_x_0_0_0_yy_xz_z_xz, g_x_0_0_0_yy_xz_z_yy, g_x_0_0_0_yy_xz_z_yz, g_x_0_0_0_yy_xz_z_zz, g_xyy_xz_z_xx, g_xyy_xz_z_xy, g_xyy_xz_z_xz, g_xyy_xz_z_yy, g_xyy_xz_z_yz, g_xyy_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_xz_z_xx[i] = 2.0 * g_xyy_xz_z_xx[i] * a_exp;

        g_x_0_0_0_yy_xz_z_xy[i] = 2.0 * g_xyy_xz_z_xy[i] * a_exp;

        g_x_0_0_0_yy_xz_z_xz[i] = 2.0 * g_xyy_xz_z_xz[i] * a_exp;

        g_x_0_0_0_yy_xz_z_yy[i] = 2.0 * g_xyy_xz_z_yy[i] * a_exp;

        g_x_0_0_0_yy_xz_z_yz[i] = 2.0 * g_xyy_xz_z_yz[i] * a_exp;

        g_x_0_0_0_yy_xz_z_zz[i] = 2.0 * g_xyy_xz_z_zz[i] * a_exp;
    }
    // integrals block (378-384)

    #pragma omp simd aligned(g_x_0_0_0_yy_yy_x_xx, g_x_0_0_0_yy_yy_x_xy, g_x_0_0_0_yy_yy_x_xz, g_x_0_0_0_yy_yy_x_yy, g_x_0_0_0_yy_yy_x_yz, g_x_0_0_0_yy_yy_x_zz, g_xyy_yy_x_xx, g_xyy_yy_x_xy, g_xyy_yy_x_xz, g_xyy_yy_x_yy, g_xyy_yy_x_yz, g_xyy_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_yy_x_xx[i] = 2.0 * g_xyy_yy_x_xx[i] * a_exp;

        g_x_0_0_0_yy_yy_x_xy[i] = 2.0 * g_xyy_yy_x_xy[i] * a_exp;

        g_x_0_0_0_yy_yy_x_xz[i] = 2.0 * g_xyy_yy_x_xz[i] * a_exp;

        g_x_0_0_0_yy_yy_x_yy[i] = 2.0 * g_xyy_yy_x_yy[i] * a_exp;

        g_x_0_0_0_yy_yy_x_yz[i] = 2.0 * g_xyy_yy_x_yz[i] * a_exp;

        g_x_0_0_0_yy_yy_x_zz[i] = 2.0 * g_xyy_yy_x_zz[i] * a_exp;
    }
    // integrals block (384-390)

    #pragma omp simd aligned(g_x_0_0_0_yy_yy_y_xx, g_x_0_0_0_yy_yy_y_xy, g_x_0_0_0_yy_yy_y_xz, g_x_0_0_0_yy_yy_y_yy, g_x_0_0_0_yy_yy_y_yz, g_x_0_0_0_yy_yy_y_zz, g_xyy_yy_y_xx, g_xyy_yy_y_xy, g_xyy_yy_y_xz, g_xyy_yy_y_yy, g_xyy_yy_y_yz, g_xyy_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_yy_y_xx[i] = 2.0 * g_xyy_yy_y_xx[i] * a_exp;

        g_x_0_0_0_yy_yy_y_xy[i] = 2.0 * g_xyy_yy_y_xy[i] * a_exp;

        g_x_0_0_0_yy_yy_y_xz[i] = 2.0 * g_xyy_yy_y_xz[i] * a_exp;

        g_x_0_0_0_yy_yy_y_yy[i] = 2.0 * g_xyy_yy_y_yy[i] * a_exp;

        g_x_0_0_0_yy_yy_y_yz[i] = 2.0 * g_xyy_yy_y_yz[i] * a_exp;

        g_x_0_0_0_yy_yy_y_zz[i] = 2.0 * g_xyy_yy_y_zz[i] * a_exp;
    }
    // integrals block (390-396)

    #pragma omp simd aligned(g_x_0_0_0_yy_yy_z_xx, g_x_0_0_0_yy_yy_z_xy, g_x_0_0_0_yy_yy_z_xz, g_x_0_0_0_yy_yy_z_yy, g_x_0_0_0_yy_yy_z_yz, g_x_0_0_0_yy_yy_z_zz, g_xyy_yy_z_xx, g_xyy_yy_z_xy, g_xyy_yy_z_xz, g_xyy_yy_z_yy, g_xyy_yy_z_yz, g_xyy_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_yy_z_xx[i] = 2.0 * g_xyy_yy_z_xx[i] * a_exp;

        g_x_0_0_0_yy_yy_z_xy[i] = 2.0 * g_xyy_yy_z_xy[i] * a_exp;

        g_x_0_0_0_yy_yy_z_xz[i] = 2.0 * g_xyy_yy_z_xz[i] * a_exp;

        g_x_0_0_0_yy_yy_z_yy[i] = 2.0 * g_xyy_yy_z_yy[i] * a_exp;

        g_x_0_0_0_yy_yy_z_yz[i] = 2.0 * g_xyy_yy_z_yz[i] * a_exp;

        g_x_0_0_0_yy_yy_z_zz[i] = 2.0 * g_xyy_yy_z_zz[i] * a_exp;
    }
    // integrals block (396-402)

    #pragma omp simd aligned(g_x_0_0_0_yy_yz_x_xx, g_x_0_0_0_yy_yz_x_xy, g_x_0_0_0_yy_yz_x_xz, g_x_0_0_0_yy_yz_x_yy, g_x_0_0_0_yy_yz_x_yz, g_x_0_0_0_yy_yz_x_zz, g_xyy_yz_x_xx, g_xyy_yz_x_xy, g_xyy_yz_x_xz, g_xyy_yz_x_yy, g_xyy_yz_x_yz, g_xyy_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_yz_x_xx[i] = 2.0 * g_xyy_yz_x_xx[i] * a_exp;

        g_x_0_0_0_yy_yz_x_xy[i] = 2.0 * g_xyy_yz_x_xy[i] * a_exp;

        g_x_0_0_0_yy_yz_x_xz[i] = 2.0 * g_xyy_yz_x_xz[i] * a_exp;

        g_x_0_0_0_yy_yz_x_yy[i] = 2.0 * g_xyy_yz_x_yy[i] * a_exp;

        g_x_0_0_0_yy_yz_x_yz[i] = 2.0 * g_xyy_yz_x_yz[i] * a_exp;

        g_x_0_0_0_yy_yz_x_zz[i] = 2.0 * g_xyy_yz_x_zz[i] * a_exp;
    }
    // integrals block (402-408)

    #pragma omp simd aligned(g_x_0_0_0_yy_yz_y_xx, g_x_0_0_0_yy_yz_y_xy, g_x_0_0_0_yy_yz_y_xz, g_x_0_0_0_yy_yz_y_yy, g_x_0_0_0_yy_yz_y_yz, g_x_0_0_0_yy_yz_y_zz, g_xyy_yz_y_xx, g_xyy_yz_y_xy, g_xyy_yz_y_xz, g_xyy_yz_y_yy, g_xyy_yz_y_yz, g_xyy_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_yz_y_xx[i] = 2.0 * g_xyy_yz_y_xx[i] * a_exp;

        g_x_0_0_0_yy_yz_y_xy[i] = 2.0 * g_xyy_yz_y_xy[i] * a_exp;

        g_x_0_0_0_yy_yz_y_xz[i] = 2.0 * g_xyy_yz_y_xz[i] * a_exp;

        g_x_0_0_0_yy_yz_y_yy[i] = 2.0 * g_xyy_yz_y_yy[i] * a_exp;

        g_x_0_0_0_yy_yz_y_yz[i] = 2.0 * g_xyy_yz_y_yz[i] * a_exp;

        g_x_0_0_0_yy_yz_y_zz[i] = 2.0 * g_xyy_yz_y_zz[i] * a_exp;
    }
    // integrals block (408-414)

    #pragma omp simd aligned(g_x_0_0_0_yy_yz_z_xx, g_x_0_0_0_yy_yz_z_xy, g_x_0_0_0_yy_yz_z_xz, g_x_0_0_0_yy_yz_z_yy, g_x_0_0_0_yy_yz_z_yz, g_x_0_0_0_yy_yz_z_zz, g_xyy_yz_z_xx, g_xyy_yz_z_xy, g_xyy_yz_z_xz, g_xyy_yz_z_yy, g_xyy_yz_z_yz, g_xyy_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_yz_z_xx[i] = 2.0 * g_xyy_yz_z_xx[i] * a_exp;

        g_x_0_0_0_yy_yz_z_xy[i] = 2.0 * g_xyy_yz_z_xy[i] * a_exp;

        g_x_0_0_0_yy_yz_z_xz[i] = 2.0 * g_xyy_yz_z_xz[i] * a_exp;

        g_x_0_0_0_yy_yz_z_yy[i] = 2.0 * g_xyy_yz_z_yy[i] * a_exp;

        g_x_0_0_0_yy_yz_z_yz[i] = 2.0 * g_xyy_yz_z_yz[i] * a_exp;

        g_x_0_0_0_yy_yz_z_zz[i] = 2.0 * g_xyy_yz_z_zz[i] * a_exp;
    }
    // integrals block (414-420)

    #pragma omp simd aligned(g_x_0_0_0_yy_zz_x_xx, g_x_0_0_0_yy_zz_x_xy, g_x_0_0_0_yy_zz_x_xz, g_x_0_0_0_yy_zz_x_yy, g_x_0_0_0_yy_zz_x_yz, g_x_0_0_0_yy_zz_x_zz, g_xyy_zz_x_xx, g_xyy_zz_x_xy, g_xyy_zz_x_xz, g_xyy_zz_x_yy, g_xyy_zz_x_yz, g_xyy_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_zz_x_xx[i] = 2.0 * g_xyy_zz_x_xx[i] * a_exp;

        g_x_0_0_0_yy_zz_x_xy[i] = 2.0 * g_xyy_zz_x_xy[i] * a_exp;

        g_x_0_0_0_yy_zz_x_xz[i] = 2.0 * g_xyy_zz_x_xz[i] * a_exp;

        g_x_0_0_0_yy_zz_x_yy[i] = 2.0 * g_xyy_zz_x_yy[i] * a_exp;

        g_x_0_0_0_yy_zz_x_yz[i] = 2.0 * g_xyy_zz_x_yz[i] * a_exp;

        g_x_0_0_0_yy_zz_x_zz[i] = 2.0 * g_xyy_zz_x_zz[i] * a_exp;
    }
    // integrals block (420-426)

    #pragma omp simd aligned(g_x_0_0_0_yy_zz_y_xx, g_x_0_0_0_yy_zz_y_xy, g_x_0_0_0_yy_zz_y_xz, g_x_0_0_0_yy_zz_y_yy, g_x_0_0_0_yy_zz_y_yz, g_x_0_0_0_yy_zz_y_zz, g_xyy_zz_y_xx, g_xyy_zz_y_xy, g_xyy_zz_y_xz, g_xyy_zz_y_yy, g_xyy_zz_y_yz, g_xyy_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_zz_y_xx[i] = 2.0 * g_xyy_zz_y_xx[i] * a_exp;

        g_x_0_0_0_yy_zz_y_xy[i] = 2.0 * g_xyy_zz_y_xy[i] * a_exp;

        g_x_0_0_0_yy_zz_y_xz[i] = 2.0 * g_xyy_zz_y_xz[i] * a_exp;

        g_x_0_0_0_yy_zz_y_yy[i] = 2.0 * g_xyy_zz_y_yy[i] * a_exp;

        g_x_0_0_0_yy_zz_y_yz[i] = 2.0 * g_xyy_zz_y_yz[i] * a_exp;

        g_x_0_0_0_yy_zz_y_zz[i] = 2.0 * g_xyy_zz_y_zz[i] * a_exp;
    }
    // integrals block (426-432)

    #pragma omp simd aligned(g_x_0_0_0_yy_zz_z_xx, g_x_0_0_0_yy_zz_z_xy, g_x_0_0_0_yy_zz_z_xz, g_x_0_0_0_yy_zz_z_yy, g_x_0_0_0_yy_zz_z_yz, g_x_0_0_0_yy_zz_z_zz, g_xyy_zz_z_xx, g_xyy_zz_z_xy, g_xyy_zz_z_xz, g_xyy_zz_z_yy, g_xyy_zz_z_yz, g_xyy_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_zz_z_xx[i] = 2.0 * g_xyy_zz_z_xx[i] * a_exp;

        g_x_0_0_0_yy_zz_z_xy[i] = 2.0 * g_xyy_zz_z_xy[i] * a_exp;

        g_x_0_0_0_yy_zz_z_xz[i] = 2.0 * g_xyy_zz_z_xz[i] * a_exp;

        g_x_0_0_0_yy_zz_z_yy[i] = 2.0 * g_xyy_zz_z_yy[i] * a_exp;

        g_x_0_0_0_yy_zz_z_yz[i] = 2.0 * g_xyy_zz_z_yz[i] * a_exp;

        g_x_0_0_0_yy_zz_z_zz[i] = 2.0 * g_xyy_zz_z_zz[i] * a_exp;
    }
    // integrals block (432-438)

    #pragma omp simd aligned(g_x_0_0_0_yz_xx_x_xx, g_x_0_0_0_yz_xx_x_xy, g_x_0_0_0_yz_xx_x_xz, g_x_0_0_0_yz_xx_x_yy, g_x_0_0_0_yz_xx_x_yz, g_x_0_0_0_yz_xx_x_zz, g_xyz_xx_x_xx, g_xyz_xx_x_xy, g_xyz_xx_x_xz, g_xyz_xx_x_yy, g_xyz_xx_x_yz, g_xyz_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_xx_x_xx[i] = 2.0 * g_xyz_xx_x_xx[i] * a_exp;

        g_x_0_0_0_yz_xx_x_xy[i] = 2.0 * g_xyz_xx_x_xy[i] * a_exp;

        g_x_0_0_0_yz_xx_x_xz[i] = 2.0 * g_xyz_xx_x_xz[i] * a_exp;

        g_x_0_0_0_yz_xx_x_yy[i] = 2.0 * g_xyz_xx_x_yy[i] * a_exp;

        g_x_0_0_0_yz_xx_x_yz[i] = 2.0 * g_xyz_xx_x_yz[i] * a_exp;

        g_x_0_0_0_yz_xx_x_zz[i] = 2.0 * g_xyz_xx_x_zz[i] * a_exp;
    }
    // integrals block (438-444)

    #pragma omp simd aligned(g_x_0_0_0_yz_xx_y_xx, g_x_0_0_0_yz_xx_y_xy, g_x_0_0_0_yz_xx_y_xz, g_x_0_0_0_yz_xx_y_yy, g_x_0_0_0_yz_xx_y_yz, g_x_0_0_0_yz_xx_y_zz, g_xyz_xx_y_xx, g_xyz_xx_y_xy, g_xyz_xx_y_xz, g_xyz_xx_y_yy, g_xyz_xx_y_yz, g_xyz_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_xx_y_xx[i] = 2.0 * g_xyz_xx_y_xx[i] * a_exp;

        g_x_0_0_0_yz_xx_y_xy[i] = 2.0 * g_xyz_xx_y_xy[i] * a_exp;

        g_x_0_0_0_yz_xx_y_xz[i] = 2.0 * g_xyz_xx_y_xz[i] * a_exp;

        g_x_0_0_0_yz_xx_y_yy[i] = 2.0 * g_xyz_xx_y_yy[i] * a_exp;

        g_x_0_0_0_yz_xx_y_yz[i] = 2.0 * g_xyz_xx_y_yz[i] * a_exp;

        g_x_0_0_0_yz_xx_y_zz[i] = 2.0 * g_xyz_xx_y_zz[i] * a_exp;
    }
    // integrals block (444-450)

    #pragma omp simd aligned(g_x_0_0_0_yz_xx_z_xx, g_x_0_0_0_yz_xx_z_xy, g_x_0_0_0_yz_xx_z_xz, g_x_0_0_0_yz_xx_z_yy, g_x_0_0_0_yz_xx_z_yz, g_x_0_0_0_yz_xx_z_zz, g_xyz_xx_z_xx, g_xyz_xx_z_xy, g_xyz_xx_z_xz, g_xyz_xx_z_yy, g_xyz_xx_z_yz, g_xyz_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_xx_z_xx[i] = 2.0 * g_xyz_xx_z_xx[i] * a_exp;

        g_x_0_0_0_yz_xx_z_xy[i] = 2.0 * g_xyz_xx_z_xy[i] * a_exp;

        g_x_0_0_0_yz_xx_z_xz[i] = 2.0 * g_xyz_xx_z_xz[i] * a_exp;

        g_x_0_0_0_yz_xx_z_yy[i] = 2.0 * g_xyz_xx_z_yy[i] * a_exp;

        g_x_0_0_0_yz_xx_z_yz[i] = 2.0 * g_xyz_xx_z_yz[i] * a_exp;

        g_x_0_0_0_yz_xx_z_zz[i] = 2.0 * g_xyz_xx_z_zz[i] * a_exp;
    }
    // integrals block (450-456)

    #pragma omp simd aligned(g_x_0_0_0_yz_xy_x_xx, g_x_0_0_0_yz_xy_x_xy, g_x_0_0_0_yz_xy_x_xz, g_x_0_0_0_yz_xy_x_yy, g_x_0_0_0_yz_xy_x_yz, g_x_0_0_0_yz_xy_x_zz, g_xyz_xy_x_xx, g_xyz_xy_x_xy, g_xyz_xy_x_xz, g_xyz_xy_x_yy, g_xyz_xy_x_yz, g_xyz_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_xy_x_xx[i] = 2.0 * g_xyz_xy_x_xx[i] * a_exp;

        g_x_0_0_0_yz_xy_x_xy[i] = 2.0 * g_xyz_xy_x_xy[i] * a_exp;

        g_x_0_0_0_yz_xy_x_xz[i] = 2.0 * g_xyz_xy_x_xz[i] * a_exp;

        g_x_0_0_0_yz_xy_x_yy[i] = 2.0 * g_xyz_xy_x_yy[i] * a_exp;

        g_x_0_0_0_yz_xy_x_yz[i] = 2.0 * g_xyz_xy_x_yz[i] * a_exp;

        g_x_0_0_0_yz_xy_x_zz[i] = 2.0 * g_xyz_xy_x_zz[i] * a_exp;
    }
    // integrals block (456-462)

    #pragma omp simd aligned(g_x_0_0_0_yz_xy_y_xx, g_x_0_0_0_yz_xy_y_xy, g_x_0_0_0_yz_xy_y_xz, g_x_0_0_0_yz_xy_y_yy, g_x_0_0_0_yz_xy_y_yz, g_x_0_0_0_yz_xy_y_zz, g_xyz_xy_y_xx, g_xyz_xy_y_xy, g_xyz_xy_y_xz, g_xyz_xy_y_yy, g_xyz_xy_y_yz, g_xyz_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_xy_y_xx[i] = 2.0 * g_xyz_xy_y_xx[i] * a_exp;

        g_x_0_0_0_yz_xy_y_xy[i] = 2.0 * g_xyz_xy_y_xy[i] * a_exp;

        g_x_0_0_0_yz_xy_y_xz[i] = 2.0 * g_xyz_xy_y_xz[i] * a_exp;

        g_x_0_0_0_yz_xy_y_yy[i] = 2.0 * g_xyz_xy_y_yy[i] * a_exp;

        g_x_0_0_0_yz_xy_y_yz[i] = 2.0 * g_xyz_xy_y_yz[i] * a_exp;

        g_x_0_0_0_yz_xy_y_zz[i] = 2.0 * g_xyz_xy_y_zz[i] * a_exp;
    }
    // integrals block (462-468)

    #pragma omp simd aligned(g_x_0_0_0_yz_xy_z_xx, g_x_0_0_0_yz_xy_z_xy, g_x_0_0_0_yz_xy_z_xz, g_x_0_0_0_yz_xy_z_yy, g_x_0_0_0_yz_xy_z_yz, g_x_0_0_0_yz_xy_z_zz, g_xyz_xy_z_xx, g_xyz_xy_z_xy, g_xyz_xy_z_xz, g_xyz_xy_z_yy, g_xyz_xy_z_yz, g_xyz_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_xy_z_xx[i] = 2.0 * g_xyz_xy_z_xx[i] * a_exp;

        g_x_0_0_0_yz_xy_z_xy[i] = 2.0 * g_xyz_xy_z_xy[i] * a_exp;

        g_x_0_0_0_yz_xy_z_xz[i] = 2.0 * g_xyz_xy_z_xz[i] * a_exp;

        g_x_0_0_0_yz_xy_z_yy[i] = 2.0 * g_xyz_xy_z_yy[i] * a_exp;

        g_x_0_0_0_yz_xy_z_yz[i] = 2.0 * g_xyz_xy_z_yz[i] * a_exp;

        g_x_0_0_0_yz_xy_z_zz[i] = 2.0 * g_xyz_xy_z_zz[i] * a_exp;
    }
    // integrals block (468-474)

    #pragma omp simd aligned(g_x_0_0_0_yz_xz_x_xx, g_x_0_0_0_yz_xz_x_xy, g_x_0_0_0_yz_xz_x_xz, g_x_0_0_0_yz_xz_x_yy, g_x_0_0_0_yz_xz_x_yz, g_x_0_0_0_yz_xz_x_zz, g_xyz_xz_x_xx, g_xyz_xz_x_xy, g_xyz_xz_x_xz, g_xyz_xz_x_yy, g_xyz_xz_x_yz, g_xyz_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_xz_x_xx[i] = 2.0 * g_xyz_xz_x_xx[i] * a_exp;

        g_x_0_0_0_yz_xz_x_xy[i] = 2.0 * g_xyz_xz_x_xy[i] * a_exp;

        g_x_0_0_0_yz_xz_x_xz[i] = 2.0 * g_xyz_xz_x_xz[i] * a_exp;

        g_x_0_0_0_yz_xz_x_yy[i] = 2.0 * g_xyz_xz_x_yy[i] * a_exp;

        g_x_0_0_0_yz_xz_x_yz[i] = 2.0 * g_xyz_xz_x_yz[i] * a_exp;

        g_x_0_0_0_yz_xz_x_zz[i] = 2.0 * g_xyz_xz_x_zz[i] * a_exp;
    }
    // integrals block (474-480)

    #pragma omp simd aligned(g_x_0_0_0_yz_xz_y_xx, g_x_0_0_0_yz_xz_y_xy, g_x_0_0_0_yz_xz_y_xz, g_x_0_0_0_yz_xz_y_yy, g_x_0_0_0_yz_xz_y_yz, g_x_0_0_0_yz_xz_y_zz, g_xyz_xz_y_xx, g_xyz_xz_y_xy, g_xyz_xz_y_xz, g_xyz_xz_y_yy, g_xyz_xz_y_yz, g_xyz_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_xz_y_xx[i] = 2.0 * g_xyz_xz_y_xx[i] * a_exp;

        g_x_0_0_0_yz_xz_y_xy[i] = 2.0 * g_xyz_xz_y_xy[i] * a_exp;

        g_x_0_0_0_yz_xz_y_xz[i] = 2.0 * g_xyz_xz_y_xz[i] * a_exp;

        g_x_0_0_0_yz_xz_y_yy[i] = 2.0 * g_xyz_xz_y_yy[i] * a_exp;

        g_x_0_0_0_yz_xz_y_yz[i] = 2.0 * g_xyz_xz_y_yz[i] * a_exp;

        g_x_0_0_0_yz_xz_y_zz[i] = 2.0 * g_xyz_xz_y_zz[i] * a_exp;
    }
    // integrals block (480-486)

    #pragma omp simd aligned(g_x_0_0_0_yz_xz_z_xx, g_x_0_0_0_yz_xz_z_xy, g_x_0_0_0_yz_xz_z_xz, g_x_0_0_0_yz_xz_z_yy, g_x_0_0_0_yz_xz_z_yz, g_x_0_0_0_yz_xz_z_zz, g_xyz_xz_z_xx, g_xyz_xz_z_xy, g_xyz_xz_z_xz, g_xyz_xz_z_yy, g_xyz_xz_z_yz, g_xyz_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_xz_z_xx[i] = 2.0 * g_xyz_xz_z_xx[i] * a_exp;

        g_x_0_0_0_yz_xz_z_xy[i] = 2.0 * g_xyz_xz_z_xy[i] * a_exp;

        g_x_0_0_0_yz_xz_z_xz[i] = 2.0 * g_xyz_xz_z_xz[i] * a_exp;

        g_x_0_0_0_yz_xz_z_yy[i] = 2.0 * g_xyz_xz_z_yy[i] * a_exp;

        g_x_0_0_0_yz_xz_z_yz[i] = 2.0 * g_xyz_xz_z_yz[i] * a_exp;

        g_x_0_0_0_yz_xz_z_zz[i] = 2.0 * g_xyz_xz_z_zz[i] * a_exp;
    }
    // integrals block (486-492)

    #pragma omp simd aligned(g_x_0_0_0_yz_yy_x_xx, g_x_0_0_0_yz_yy_x_xy, g_x_0_0_0_yz_yy_x_xz, g_x_0_0_0_yz_yy_x_yy, g_x_0_0_0_yz_yy_x_yz, g_x_0_0_0_yz_yy_x_zz, g_xyz_yy_x_xx, g_xyz_yy_x_xy, g_xyz_yy_x_xz, g_xyz_yy_x_yy, g_xyz_yy_x_yz, g_xyz_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_yy_x_xx[i] = 2.0 * g_xyz_yy_x_xx[i] * a_exp;

        g_x_0_0_0_yz_yy_x_xy[i] = 2.0 * g_xyz_yy_x_xy[i] * a_exp;

        g_x_0_0_0_yz_yy_x_xz[i] = 2.0 * g_xyz_yy_x_xz[i] * a_exp;

        g_x_0_0_0_yz_yy_x_yy[i] = 2.0 * g_xyz_yy_x_yy[i] * a_exp;

        g_x_0_0_0_yz_yy_x_yz[i] = 2.0 * g_xyz_yy_x_yz[i] * a_exp;

        g_x_0_0_0_yz_yy_x_zz[i] = 2.0 * g_xyz_yy_x_zz[i] * a_exp;
    }
    // integrals block (492-498)

    #pragma omp simd aligned(g_x_0_0_0_yz_yy_y_xx, g_x_0_0_0_yz_yy_y_xy, g_x_0_0_0_yz_yy_y_xz, g_x_0_0_0_yz_yy_y_yy, g_x_0_0_0_yz_yy_y_yz, g_x_0_0_0_yz_yy_y_zz, g_xyz_yy_y_xx, g_xyz_yy_y_xy, g_xyz_yy_y_xz, g_xyz_yy_y_yy, g_xyz_yy_y_yz, g_xyz_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_yy_y_xx[i] = 2.0 * g_xyz_yy_y_xx[i] * a_exp;

        g_x_0_0_0_yz_yy_y_xy[i] = 2.0 * g_xyz_yy_y_xy[i] * a_exp;

        g_x_0_0_0_yz_yy_y_xz[i] = 2.0 * g_xyz_yy_y_xz[i] * a_exp;

        g_x_0_0_0_yz_yy_y_yy[i] = 2.0 * g_xyz_yy_y_yy[i] * a_exp;

        g_x_0_0_0_yz_yy_y_yz[i] = 2.0 * g_xyz_yy_y_yz[i] * a_exp;

        g_x_0_0_0_yz_yy_y_zz[i] = 2.0 * g_xyz_yy_y_zz[i] * a_exp;
    }
    // integrals block (498-504)

    #pragma omp simd aligned(g_x_0_0_0_yz_yy_z_xx, g_x_0_0_0_yz_yy_z_xy, g_x_0_0_0_yz_yy_z_xz, g_x_0_0_0_yz_yy_z_yy, g_x_0_0_0_yz_yy_z_yz, g_x_0_0_0_yz_yy_z_zz, g_xyz_yy_z_xx, g_xyz_yy_z_xy, g_xyz_yy_z_xz, g_xyz_yy_z_yy, g_xyz_yy_z_yz, g_xyz_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_yy_z_xx[i] = 2.0 * g_xyz_yy_z_xx[i] * a_exp;

        g_x_0_0_0_yz_yy_z_xy[i] = 2.0 * g_xyz_yy_z_xy[i] * a_exp;

        g_x_0_0_0_yz_yy_z_xz[i] = 2.0 * g_xyz_yy_z_xz[i] * a_exp;

        g_x_0_0_0_yz_yy_z_yy[i] = 2.0 * g_xyz_yy_z_yy[i] * a_exp;

        g_x_0_0_0_yz_yy_z_yz[i] = 2.0 * g_xyz_yy_z_yz[i] * a_exp;

        g_x_0_0_0_yz_yy_z_zz[i] = 2.0 * g_xyz_yy_z_zz[i] * a_exp;
    }
    // integrals block (504-510)

    #pragma omp simd aligned(g_x_0_0_0_yz_yz_x_xx, g_x_0_0_0_yz_yz_x_xy, g_x_0_0_0_yz_yz_x_xz, g_x_0_0_0_yz_yz_x_yy, g_x_0_0_0_yz_yz_x_yz, g_x_0_0_0_yz_yz_x_zz, g_xyz_yz_x_xx, g_xyz_yz_x_xy, g_xyz_yz_x_xz, g_xyz_yz_x_yy, g_xyz_yz_x_yz, g_xyz_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_yz_x_xx[i] = 2.0 * g_xyz_yz_x_xx[i] * a_exp;

        g_x_0_0_0_yz_yz_x_xy[i] = 2.0 * g_xyz_yz_x_xy[i] * a_exp;

        g_x_0_0_0_yz_yz_x_xz[i] = 2.0 * g_xyz_yz_x_xz[i] * a_exp;

        g_x_0_0_0_yz_yz_x_yy[i] = 2.0 * g_xyz_yz_x_yy[i] * a_exp;

        g_x_0_0_0_yz_yz_x_yz[i] = 2.0 * g_xyz_yz_x_yz[i] * a_exp;

        g_x_0_0_0_yz_yz_x_zz[i] = 2.0 * g_xyz_yz_x_zz[i] * a_exp;
    }
    // integrals block (510-516)

    #pragma omp simd aligned(g_x_0_0_0_yz_yz_y_xx, g_x_0_0_0_yz_yz_y_xy, g_x_0_0_0_yz_yz_y_xz, g_x_0_0_0_yz_yz_y_yy, g_x_0_0_0_yz_yz_y_yz, g_x_0_0_0_yz_yz_y_zz, g_xyz_yz_y_xx, g_xyz_yz_y_xy, g_xyz_yz_y_xz, g_xyz_yz_y_yy, g_xyz_yz_y_yz, g_xyz_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_yz_y_xx[i] = 2.0 * g_xyz_yz_y_xx[i] * a_exp;

        g_x_0_0_0_yz_yz_y_xy[i] = 2.0 * g_xyz_yz_y_xy[i] * a_exp;

        g_x_0_0_0_yz_yz_y_xz[i] = 2.0 * g_xyz_yz_y_xz[i] * a_exp;

        g_x_0_0_0_yz_yz_y_yy[i] = 2.0 * g_xyz_yz_y_yy[i] * a_exp;

        g_x_0_0_0_yz_yz_y_yz[i] = 2.0 * g_xyz_yz_y_yz[i] * a_exp;

        g_x_0_0_0_yz_yz_y_zz[i] = 2.0 * g_xyz_yz_y_zz[i] * a_exp;
    }
    // integrals block (516-522)

    #pragma omp simd aligned(g_x_0_0_0_yz_yz_z_xx, g_x_0_0_0_yz_yz_z_xy, g_x_0_0_0_yz_yz_z_xz, g_x_0_0_0_yz_yz_z_yy, g_x_0_0_0_yz_yz_z_yz, g_x_0_0_0_yz_yz_z_zz, g_xyz_yz_z_xx, g_xyz_yz_z_xy, g_xyz_yz_z_xz, g_xyz_yz_z_yy, g_xyz_yz_z_yz, g_xyz_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_yz_z_xx[i] = 2.0 * g_xyz_yz_z_xx[i] * a_exp;

        g_x_0_0_0_yz_yz_z_xy[i] = 2.0 * g_xyz_yz_z_xy[i] * a_exp;

        g_x_0_0_0_yz_yz_z_xz[i] = 2.0 * g_xyz_yz_z_xz[i] * a_exp;

        g_x_0_0_0_yz_yz_z_yy[i] = 2.0 * g_xyz_yz_z_yy[i] * a_exp;

        g_x_0_0_0_yz_yz_z_yz[i] = 2.0 * g_xyz_yz_z_yz[i] * a_exp;

        g_x_0_0_0_yz_yz_z_zz[i] = 2.0 * g_xyz_yz_z_zz[i] * a_exp;
    }
    // integrals block (522-528)

    #pragma omp simd aligned(g_x_0_0_0_yz_zz_x_xx, g_x_0_0_0_yz_zz_x_xy, g_x_0_0_0_yz_zz_x_xz, g_x_0_0_0_yz_zz_x_yy, g_x_0_0_0_yz_zz_x_yz, g_x_0_0_0_yz_zz_x_zz, g_xyz_zz_x_xx, g_xyz_zz_x_xy, g_xyz_zz_x_xz, g_xyz_zz_x_yy, g_xyz_zz_x_yz, g_xyz_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_zz_x_xx[i] = 2.0 * g_xyz_zz_x_xx[i] * a_exp;

        g_x_0_0_0_yz_zz_x_xy[i] = 2.0 * g_xyz_zz_x_xy[i] * a_exp;

        g_x_0_0_0_yz_zz_x_xz[i] = 2.0 * g_xyz_zz_x_xz[i] * a_exp;

        g_x_0_0_0_yz_zz_x_yy[i] = 2.0 * g_xyz_zz_x_yy[i] * a_exp;

        g_x_0_0_0_yz_zz_x_yz[i] = 2.0 * g_xyz_zz_x_yz[i] * a_exp;

        g_x_0_0_0_yz_zz_x_zz[i] = 2.0 * g_xyz_zz_x_zz[i] * a_exp;
    }
    // integrals block (528-534)

    #pragma omp simd aligned(g_x_0_0_0_yz_zz_y_xx, g_x_0_0_0_yz_zz_y_xy, g_x_0_0_0_yz_zz_y_xz, g_x_0_0_0_yz_zz_y_yy, g_x_0_0_0_yz_zz_y_yz, g_x_0_0_0_yz_zz_y_zz, g_xyz_zz_y_xx, g_xyz_zz_y_xy, g_xyz_zz_y_xz, g_xyz_zz_y_yy, g_xyz_zz_y_yz, g_xyz_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_zz_y_xx[i] = 2.0 * g_xyz_zz_y_xx[i] * a_exp;

        g_x_0_0_0_yz_zz_y_xy[i] = 2.0 * g_xyz_zz_y_xy[i] * a_exp;

        g_x_0_0_0_yz_zz_y_xz[i] = 2.0 * g_xyz_zz_y_xz[i] * a_exp;

        g_x_0_0_0_yz_zz_y_yy[i] = 2.0 * g_xyz_zz_y_yy[i] * a_exp;

        g_x_0_0_0_yz_zz_y_yz[i] = 2.0 * g_xyz_zz_y_yz[i] * a_exp;

        g_x_0_0_0_yz_zz_y_zz[i] = 2.0 * g_xyz_zz_y_zz[i] * a_exp;
    }
    // integrals block (534-540)

    #pragma omp simd aligned(g_x_0_0_0_yz_zz_z_xx, g_x_0_0_0_yz_zz_z_xy, g_x_0_0_0_yz_zz_z_xz, g_x_0_0_0_yz_zz_z_yy, g_x_0_0_0_yz_zz_z_yz, g_x_0_0_0_yz_zz_z_zz, g_xyz_zz_z_xx, g_xyz_zz_z_xy, g_xyz_zz_z_xz, g_xyz_zz_z_yy, g_xyz_zz_z_yz, g_xyz_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_zz_z_xx[i] = 2.0 * g_xyz_zz_z_xx[i] * a_exp;

        g_x_0_0_0_yz_zz_z_xy[i] = 2.0 * g_xyz_zz_z_xy[i] * a_exp;

        g_x_0_0_0_yz_zz_z_xz[i] = 2.0 * g_xyz_zz_z_xz[i] * a_exp;

        g_x_0_0_0_yz_zz_z_yy[i] = 2.0 * g_xyz_zz_z_yy[i] * a_exp;

        g_x_0_0_0_yz_zz_z_yz[i] = 2.0 * g_xyz_zz_z_yz[i] * a_exp;

        g_x_0_0_0_yz_zz_z_zz[i] = 2.0 * g_xyz_zz_z_zz[i] * a_exp;
    }
    // integrals block (540-546)

    #pragma omp simd aligned(g_x_0_0_0_zz_xx_x_xx, g_x_0_0_0_zz_xx_x_xy, g_x_0_0_0_zz_xx_x_xz, g_x_0_0_0_zz_xx_x_yy, g_x_0_0_0_zz_xx_x_yz, g_x_0_0_0_zz_xx_x_zz, g_xzz_xx_x_xx, g_xzz_xx_x_xy, g_xzz_xx_x_xz, g_xzz_xx_x_yy, g_xzz_xx_x_yz, g_xzz_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_xx_x_xx[i] = 2.0 * g_xzz_xx_x_xx[i] * a_exp;

        g_x_0_0_0_zz_xx_x_xy[i] = 2.0 * g_xzz_xx_x_xy[i] * a_exp;

        g_x_0_0_0_zz_xx_x_xz[i] = 2.0 * g_xzz_xx_x_xz[i] * a_exp;

        g_x_0_0_0_zz_xx_x_yy[i] = 2.0 * g_xzz_xx_x_yy[i] * a_exp;

        g_x_0_0_0_zz_xx_x_yz[i] = 2.0 * g_xzz_xx_x_yz[i] * a_exp;

        g_x_0_0_0_zz_xx_x_zz[i] = 2.0 * g_xzz_xx_x_zz[i] * a_exp;
    }
    // integrals block (546-552)

    #pragma omp simd aligned(g_x_0_0_0_zz_xx_y_xx, g_x_0_0_0_zz_xx_y_xy, g_x_0_0_0_zz_xx_y_xz, g_x_0_0_0_zz_xx_y_yy, g_x_0_0_0_zz_xx_y_yz, g_x_0_0_0_zz_xx_y_zz, g_xzz_xx_y_xx, g_xzz_xx_y_xy, g_xzz_xx_y_xz, g_xzz_xx_y_yy, g_xzz_xx_y_yz, g_xzz_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_xx_y_xx[i] = 2.0 * g_xzz_xx_y_xx[i] * a_exp;

        g_x_0_0_0_zz_xx_y_xy[i] = 2.0 * g_xzz_xx_y_xy[i] * a_exp;

        g_x_0_0_0_zz_xx_y_xz[i] = 2.0 * g_xzz_xx_y_xz[i] * a_exp;

        g_x_0_0_0_zz_xx_y_yy[i] = 2.0 * g_xzz_xx_y_yy[i] * a_exp;

        g_x_0_0_0_zz_xx_y_yz[i] = 2.0 * g_xzz_xx_y_yz[i] * a_exp;

        g_x_0_0_0_zz_xx_y_zz[i] = 2.0 * g_xzz_xx_y_zz[i] * a_exp;
    }
    // integrals block (552-558)

    #pragma omp simd aligned(g_x_0_0_0_zz_xx_z_xx, g_x_0_0_0_zz_xx_z_xy, g_x_0_0_0_zz_xx_z_xz, g_x_0_0_0_zz_xx_z_yy, g_x_0_0_0_zz_xx_z_yz, g_x_0_0_0_zz_xx_z_zz, g_xzz_xx_z_xx, g_xzz_xx_z_xy, g_xzz_xx_z_xz, g_xzz_xx_z_yy, g_xzz_xx_z_yz, g_xzz_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_xx_z_xx[i] = 2.0 * g_xzz_xx_z_xx[i] * a_exp;

        g_x_0_0_0_zz_xx_z_xy[i] = 2.0 * g_xzz_xx_z_xy[i] * a_exp;

        g_x_0_0_0_zz_xx_z_xz[i] = 2.0 * g_xzz_xx_z_xz[i] * a_exp;

        g_x_0_0_0_zz_xx_z_yy[i] = 2.0 * g_xzz_xx_z_yy[i] * a_exp;

        g_x_0_0_0_zz_xx_z_yz[i] = 2.0 * g_xzz_xx_z_yz[i] * a_exp;

        g_x_0_0_0_zz_xx_z_zz[i] = 2.0 * g_xzz_xx_z_zz[i] * a_exp;
    }
    // integrals block (558-564)

    #pragma omp simd aligned(g_x_0_0_0_zz_xy_x_xx, g_x_0_0_0_zz_xy_x_xy, g_x_0_0_0_zz_xy_x_xz, g_x_0_0_0_zz_xy_x_yy, g_x_0_0_0_zz_xy_x_yz, g_x_0_0_0_zz_xy_x_zz, g_xzz_xy_x_xx, g_xzz_xy_x_xy, g_xzz_xy_x_xz, g_xzz_xy_x_yy, g_xzz_xy_x_yz, g_xzz_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_xy_x_xx[i] = 2.0 * g_xzz_xy_x_xx[i] * a_exp;

        g_x_0_0_0_zz_xy_x_xy[i] = 2.0 * g_xzz_xy_x_xy[i] * a_exp;

        g_x_0_0_0_zz_xy_x_xz[i] = 2.0 * g_xzz_xy_x_xz[i] * a_exp;

        g_x_0_0_0_zz_xy_x_yy[i] = 2.0 * g_xzz_xy_x_yy[i] * a_exp;

        g_x_0_0_0_zz_xy_x_yz[i] = 2.0 * g_xzz_xy_x_yz[i] * a_exp;

        g_x_0_0_0_zz_xy_x_zz[i] = 2.0 * g_xzz_xy_x_zz[i] * a_exp;
    }
    // integrals block (564-570)

    #pragma omp simd aligned(g_x_0_0_0_zz_xy_y_xx, g_x_0_0_0_zz_xy_y_xy, g_x_0_0_0_zz_xy_y_xz, g_x_0_0_0_zz_xy_y_yy, g_x_0_0_0_zz_xy_y_yz, g_x_0_0_0_zz_xy_y_zz, g_xzz_xy_y_xx, g_xzz_xy_y_xy, g_xzz_xy_y_xz, g_xzz_xy_y_yy, g_xzz_xy_y_yz, g_xzz_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_xy_y_xx[i] = 2.0 * g_xzz_xy_y_xx[i] * a_exp;

        g_x_0_0_0_zz_xy_y_xy[i] = 2.0 * g_xzz_xy_y_xy[i] * a_exp;

        g_x_0_0_0_zz_xy_y_xz[i] = 2.0 * g_xzz_xy_y_xz[i] * a_exp;

        g_x_0_0_0_zz_xy_y_yy[i] = 2.0 * g_xzz_xy_y_yy[i] * a_exp;

        g_x_0_0_0_zz_xy_y_yz[i] = 2.0 * g_xzz_xy_y_yz[i] * a_exp;

        g_x_0_0_0_zz_xy_y_zz[i] = 2.0 * g_xzz_xy_y_zz[i] * a_exp;
    }
    // integrals block (570-576)

    #pragma omp simd aligned(g_x_0_0_0_zz_xy_z_xx, g_x_0_0_0_zz_xy_z_xy, g_x_0_0_0_zz_xy_z_xz, g_x_0_0_0_zz_xy_z_yy, g_x_0_0_0_zz_xy_z_yz, g_x_0_0_0_zz_xy_z_zz, g_xzz_xy_z_xx, g_xzz_xy_z_xy, g_xzz_xy_z_xz, g_xzz_xy_z_yy, g_xzz_xy_z_yz, g_xzz_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_xy_z_xx[i] = 2.0 * g_xzz_xy_z_xx[i] * a_exp;

        g_x_0_0_0_zz_xy_z_xy[i] = 2.0 * g_xzz_xy_z_xy[i] * a_exp;

        g_x_0_0_0_zz_xy_z_xz[i] = 2.0 * g_xzz_xy_z_xz[i] * a_exp;

        g_x_0_0_0_zz_xy_z_yy[i] = 2.0 * g_xzz_xy_z_yy[i] * a_exp;

        g_x_0_0_0_zz_xy_z_yz[i] = 2.0 * g_xzz_xy_z_yz[i] * a_exp;

        g_x_0_0_0_zz_xy_z_zz[i] = 2.0 * g_xzz_xy_z_zz[i] * a_exp;
    }
    // integrals block (576-582)

    #pragma omp simd aligned(g_x_0_0_0_zz_xz_x_xx, g_x_0_0_0_zz_xz_x_xy, g_x_0_0_0_zz_xz_x_xz, g_x_0_0_0_zz_xz_x_yy, g_x_0_0_0_zz_xz_x_yz, g_x_0_0_0_zz_xz_x_zz, g_xzz_xz_x_xx, g_xzz_xz_x_xy, g_xzz_xz_x_xz, g_xzz_xz_x_yy, g_xzz_xz_x_yz, g_xzz_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_xz_x_xx[i] = 2.0 * g_xzz_xz_x_xx[i] * a_exp;

        g_x_0_0_0_zz_xz_x_xy[i] = 2.0 * g_xzz_xz_x_xy[i] * a_exp;

        g_x_0_0_0_zz_xz_x_xz[i] = 2.0 * g_xzz_xz_x_xz[i] * a_exp;

        g_x_0_0_0_zz_xz_x_yy[i] = 2.0 * g_xzz_xz_x_yy[i] * a_exp;

        g_x_0_0_0_zz_xz_x_yz[i] = 2.0 * g_xzz_xz_x_yz[i] * a_exp;

        g_x_0_0_0_zz_xz_x_zz[i] = 2.0 * g_xzz_xz_x_zz[i] * a_exp;
    }
    // integrals block (582-588)

    #pragma omp simd aligned(g_x_0_0_0_zz_xz_y_xx, g_x_0_0_0_zz_xz_y_xy, g_x_0_0_0_zz_xz_y_xz, g_x_0_0_0_zz_xz_y_yy, g_x_0_0_0_zz_xz_y_yz, g_x_0_0_0_zz_xz_y_zz, g_xzz_xz_y_xx, g_xzz_xz_y_xy, g_xzz_xz_y_xz, g_xzz_xz_y_yy, g_xzz_xz_y_yz, g_xzz_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_xz_y_xx[i] = 2.0 * g_xzz_xz_y_xx[i] * a_exp;

        g_x_0_0_0_zz_xz_y_xy[i] = 2.0 * g_xzz_xz_y_xy[i] * a_exp;

        g_x_0_0_0_zz_xz_y_xz[i] = 2.0 * g_xzz_xz_y_xz[i] * a_exp;

        g_x_0_0_0_zz_xz_y_yy[i] = 2.0 * g_xzz_xz_y_yy[i] * a_exp;

        g_x_0_0_0_zz_xz_y_yz[i] = 2.0 * g_xzz_xz_y_yz[i] * a_exp;

        g_x_0_0_0_zz_xz_y_zz[i] = 2.0 * g_xzz_xz_y_zz[i] * a_exp;
    }
    // integrals block (588-594)

    #pragma omp simd aligned(g_x_0_0_0_zz_xz_z_xx, g_x_0_0_0_zz_xz_z_xy, g_x_0_0_0_zz_xz_z_xz, g_x_0_0_0_zz_xz_z_yy, g_x_0_0_0_zz_xz_z_yz, g_x_0_0_0_zz_xz_z_zz, g_xzz_xz_z_xx, g_xzz_xz_z_xy, g_xzz_xz_z_xz, g_xzz_xz_z_yy, g_xzz_xz_z_yz, g_xzz_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_xz_z_xx[i] = 2.0 * g_xzz_xz_z_xx[i] * a_exp;

        g_x_0_0_0_zz_xz_z_xy[i] = 2.0 * g_xzz_xz_z_xy[i] * a_exp;

        g_x_0_0_0_zz_xz_z_xz[i] = 2.0 * g_xzz_xz_z_xz[i] * a_exp;

        g_x_0_0_0_zz_xz_z_yy[i] = 2.0 * g_xzz_xz_z_yy[i] * a_exp;

        g_x_0_0_0_zz_xz_z_yz[i] = 2.0 * g_xzz_xz_z_yz[i] * a_exp;

        g_x_0_0_0_zz_xz_z_zz[i] = 2.0 * g_xzz_xz_z_zz[i] * a_exp;
    }
    // integrals block (594-600)

    #pragma omp simd aligned(g_x_0_0_0_zz_yy_x_xx, g_x_0_0_0_zz_yy_x_xy, g_x_0_0_0_zz_yy_x_xz, g_x_0_0_0_zz_yy_x_yy, g_x_0_0_0_zz_yy_x_yz, g_x_0_0_0_zz_yy_x_zz, g_xzz_yy_x_xx, g_xzz_yy_x_xy, g_xzz_yy_x_xz, g_xzz_yy_x_yy, g_xzz_yy_x_yz, g_xzz_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_yy_x_xx[i] = 2.0 * g_xzz_yy_x_xx[i] * a_exp;

        g_x_0_0_0_zz_yy_x_xy[i] = 2.0 * g_xzz_yy_x_xy[i] * a_exp;

        g_x_0_0_0_zz_yy_x_xz[i] = 2.0 * g_xzz_yy_x_xz[i] * a_exp;

        g_x_0_0_0_zz_yy_x_yy[i] = 2.0 * g_xzz_yy_x_yy[i] * a_exp;

        g_x_0_0_0_zz_yy_x_yz[i] = 2.0 * g_xzz_yy_x_yz[i] * a_exp;

        g_x_0_0_0_zz_yy_x_zz[i] = 2.0 * g_xzz_yy_x_zz[i] * a_exp;
    }
    // integrals block (600-606)

    #pragma omp simd aligned(g_x_0_0_0_zz_yy_y_xx, g_x_0_0_0_zz_yy_y_xy, g_x_0_0_0_zz_yy_y_xz, g_x_0_0_0_zz_yy_y_yy, g_x_0_0_0_zz_yy_y_yz, g_x_0_0_0_zz_yy_y_zz, g_xzz_yy_y_xx, g_xzz_yy_y_xy, g_xzz_yy_y_xz, g_xzz_yy_y_yy, g_xzz_yy_y_yz, g_xzz_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_yy_y_xx[i] = 2.0 * g_xzz_yy_y_xx[i] * a_exp;

        g_x_0_0_0_zz_yy_y_xy[i] = 2.0 * g_xzz_yy_y_xy[i] * a_exp;

        g_x_0_0_0_zz_yy_y_xz[i] = 2.0 * g_xzz_yy_y_xz[i] * a_exp;

        g_x_0_0_0_zz_yy_y_yy[i] = 2.0 * g_xzz_yy_y_yy[i] * a_exp;

        g_x_0_0_0_zz_yy_y_yz[i] = 2.0 * g_xzz_yy_y_yz[i] * a_exp;

        g_x_0_0_0_zz_yy_y_zz[i] = 2.0 * g_xzz_yy_y_zz[i] * a_exp;
    }
    // integrals block (606-612)

    #pragma omp simd aligned(g_x_0_0_0_zz_yy_z_xx, g_x_0_0_0_zz_yy_z_xy, g_x_0_0_0_zz_yy_z_xz, g_x_0_0_0_zz_yy_z_yy, g_x_0_0_0_zz_yy_z_yz, g_x_0_0_0_zz_yy_z_zz, g_xzz_yy_z_xx, g_xzz_yy_z_xy, g_xzz_yy_z_xz, g_xzz_yy_z_yy, g_xzz_yy_z_yz, g_xzz_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_yy_z_xx[i] = 2.0 * g_xzz_yy_z_xx[i] * a_exp;

        g_x_0_0_0_zz_yy_z_xy[i] = 2.0 * g_xzz_yy_z_xy[i] * a_exp;

        g_x_0_0_0_zz_yy_z_xz[i] = 2.0 * g_xzz_yy_z_xz[i] * a_exp;

        g_x_0_0_0_zz_yy_z_yy[i] = 2.0 * g_xzz_yy_z_yy[i] * a_exp;

        g_x_0_0_0_zz_yy_z_yz[i] = 2.0 * g_xzz_yy_z_yz[i] * a_exp;

        g_x_0_0_0_zz_yy_z_zz[i] = 2.0 * g_xzz_yy_z_zz[i] * a_exp;
    }
    // integrals block (612-618)

    #pragma omp simd aligned(g_x_0_0_0_zz_yz_x_xx, g_x_0_0_0_zz_yz_x_xy, g_x_0_0_0_zz_yz_x_xz, g_x_0_0_0_zz_yz_x_yy, g_x_0_0_0_zz_yz_x_yz, g_x_0_0_0_zz_yz_x_zz, g_xzz_yz_x_xx, g_xzz_yz_x_xy, g_xzz_yz_x_xz, g_xzz_yz_x_yy, g_xzz_yz_x_yz, g_xzz_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_yz_x_xx[i] = 2.0 * g_xzz_yz_x_xx[i] * a_exp;

        g_x_0_0_0_zz_yz_x_xy[i] = 2.0 * g_xzz_yz_x_xy[i] * a_exp;

        g_x_0_0_0_zz_yz_x_xz[i] = 2.0 * g_xzz_yz_x_xz[i] * a_exp;

        g_x_0_0_0_zz_yz_x_yy[i] = 2.0 * g_xzz_yz_x_yy[i] * a_exp;

        g_x_0_0_0_zz_yz_x_yz[i] = 2.0 * g_xzz_yz_x_yz[i] * a_exp;

        g_x_0_0_0_zz_yz_x_zz[i] = 2.0 * g_xzz_yz_x_zz[i] * a_exp;
    }
    // integrals block (618-624)

    #pragma omp simd aligned(g_x_0_0_0_zz_yz_y_xx, g_x_0_0_0_zz_yz_y_xy, g_x_0_0_0_zz_yz_y_xz, g_x_0_0_0_zz_yz_y_yy, g_x_0_0_0_zz_yz_y_yz, g_x_0_0_0_zz_yz_y_zz, g_xzz_yz_y_xx, g_xzz_yz_y_xy, g_xzz_yz_y_xz, g_xzz_yz_y_yy, g_xzz_yz_y_yz, g_xzz_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_yz_y_xx[i] = 2.0 * g_xzz_yz_y_xx[i] * a_exp;

        g_x_0_0_0_zz_yz_y_xy[i] = 2.0 * g_xzz_yz_y_xy[i] * a_exp;

        g_x_0_0_0_zz_yz_y_xz[i] = 2.0 * g_xzz_yz_y_xz[i] * a_exp;

        g_x_0_0_0_zz_yz_y_yy[i] = 2.0 * g_xzz_yz_y_yy[i] * a_exp;

        g_x_0_0_0_zz_yz_y_yz[i] = 2.0 * g_xzz_yz_y_yz[i] * a_exp;

        g_x_0_0_0_zz_yz_y_zz[i] = 2.0 * g_xzz_yz_y_zz[i] * a_exp;
    }
    // integrals block (624-630)

    #pragma omp simd aligned(g_x_0_0_0_zz_yz_z_xx, g_x_0_0_0_zz_yz_z_xy, g_x_0_0_0_zz_yz_z_xz, g_x_0_0_0_zz_yz_z_yy, g_x_0_0_0_zz_yz_z_yz, g_x_0_0_0_zz_yz_z_zz, g_xzz_yz_z_xx, g_xzz_yz_z_xy, g_xzz_yz_z_xz, g_xzz_yz_z_yy, g_xzz_yz_z_yz, g_xzz_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_yz_z_xx[i] = 2.0 * g_xzz_yz_z_xx[i] * a_exp;

        g_x_0_0_0_zz_yz_z_xy[i] = 2.0 * g_xzz_yz_z_xy[i] * a_exp;

        g_x_0_0_0_zz_yz_z_xz[i] = 2.0 * g_xzz_yz_z_xz[i] * a_exp;

        g_x_0_0_0_zz_yz_z_yy[i] = 2.0 * g_xzz_yz_z_yy[i] * a_exp;

        g_x_0_0_0_zz_yz_z_yz[i] = 2.0 * g_xzz_yz_z_yz[i] * a_exp;

        g_x_0_0_0_zz_yz_z_zz[i] = 2.0 * g_xzz_yz_z_zz[i] * a_exp;
    }
    // integrals block (630-636)

    #pragma omp simd aligned(g_x_0_0_0_zz_zz_x_xx, g_x_0_0_0_zz_zz_x_xy, g_x_0_0_0_zz_zz_x_xz, g_x_0_0_0_zz_zz_x_yy, g_x_0_0_0_zz_zz_x_yz, g_x_0_0_0_zz_zz_x_zz, g_xzz_zz_x_xx, g_xzz_zz_x_xy, g_xzz_zz_x_xz, g_xzz_zz_x_yy, g_xzz_zz_x_yz, g_xzz_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_zz_x_xx[i] = 2.0 * g_xzz_zz_x_xx[i] * a_exp;

        g_x_0_0_0_zz_zz_x_xy[i] = 2.0 * g_xzz_zz_x_xy[i] * a_exp;

        g_x_0_0_0_zz_zz_x_xz[i] = 2.0 * g_xzz_zz_x_xz[i] * a_exp;

        g_x_0_0_0_zz_zz_x_yy[i] = 2.0 * g_xzz_zz_x_yy[i] * a_exp;

        g_x_0_0_0_zz_zz_x_yz[i] = 2.0 * g_xzz_zz_x_yz[i] * a_exp;

        g_x_0_0_0_zz_zz_x_zz[i] = 2.0 * g_xzz_zz_x_zz[i] * a_exp;
    }
    // integrals block (636-642)

    #pragma omp simd aligned(g_x_0_0_0_zz_zz_y_xx, g_x_0_0_0_zz_zz_y_xy, g_x_0_0_0_zz_zz_y_xz, g_x_0_0_0_zz_zz_y_yy, g_x_0_0_0_zz_zz_y_yz, g_x_0_0_0_zz_zz_y_zz, g_xzz_zz_y_xx, g_xzz_zz_y_xy, g_xzz_zz_y_xz, g_xzz_zz_y_yy, g_xzz_zz_y_yz, g_xzz_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_zz_y_xx[i] = 2.0 * g_xzz_zz_y_xx[i] * a_exp;

        g_x_0_0_0_zz_zz_y_xy[i] = 2.0 * g_xzz_zz_y_xy[i] * a_exp;

        g_x_0_0_0_zz_zz_y_xz[i] = 2.0 * g_xzz_zz_y_xz[i] * a_exp;

        g_x_0_0_0_zz_zz_y_yy[i] = 2.0 * g_xzz_zz_y_yy[i] * a_exp;

        g_x_0_0_0_zz_zz_y_yz[i] = 2.0 * g_xzz_zz_y_yz[i] * a_exp;

        g_x_0_0_0_zz_zz_y_zz[i] = 2.0 * g_xzz_zz_y_zz[i] * a_exp;
    }
    // integrals block (642-648)

    #pragma omp simd aligned(g_x_0_0_0_zz_zz_z_xx, g_x_0_0_0_zz_zz_z_xy, g_x_0_0_0_zz_zz_z_xz, g_x_0_0_0_zz_zz_z_yy, g_x_0_0_0_zz_zz_z_yz, g_x_0_0_0_zz_zz_z_zz, g_xzz_zz_z_xx, g_xzz_zz_z_xy, g_xzz_zz_z_xz, g_xzz_zz_z_yy, g_xzz_zz_z_yz, g_xzz_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_zz_z_xx[i] = 2.0 * g_xzz_zz_z_xx[i] * a_exp;

        g_x_0_0_0_zz_zz_z_xy[i] = 2.0 * g_xzz_zz_z_xy[i] * a_exp;

        g_x_0_0_0_zz_zz_z_xz[i] = 2.0 * g_xzz_zz_z_xz[i] * a_exp;

        g_x_0_0_0_zz_zz_z_yy[i] = 2.0 * g_xzz_zz_z_yy[i] * a_exp;

        g_x_0_0_0_zz_zz_z_yz[i] = 2.0 * g_xzz_zz_z_yz[i] * a_exp;

        g_x_0_0_0_zz_zz_z_zz[i] = 2.0 * g_xzz_zz_z_zz[i] * a_exp;
    }
    // integrals block (648-654)

    #pragma omp simd aligned(g_xxy_xx_x_xx, g_xxy_xx_x_xy, g_xxy_xx_x_xz, g_xxy_xx_x_yy, g_xxy_xx_x_yz, g_xxy_xx_x_zz, g_y_0_0_0_xx_xx_x_xx, g_y_0_0_0_xx_xx_x_xy, g_y_0_0_0_xx_xx_x_xz, g_y_0_0_0_xx_xx_x_yy, g_y_0_0_0_xx_xx_x_yz, g_y_0_0_0_xx_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_xx_x_xx[i] = 2.0 * g_xxy_xx_x_xx[i] * a_exp;

        g_y_0_0_0_xx_xx_x_xy[i] = 2.0 * g_xxy_xx_x_xy[i] * a_exp;

        g_y_0_0_0_xx_xx_x_xz[i] = 2.0 * g_xxy_xx_x_xz[i] * a_exp;

        g_y_0_0_0_xx_xx_x_yy[i] = 2.0 * g_xxy_xx_x_yy[i] * a_exp;

        g_y_0_0_0_xx_xx_x_yz[i] = 2.0 * g_xxy_xx_x_yz[i] * a_exp;

        g_y_0_0_0_xx_xx_x_zz[i] = 2.0 * g_xxy_xx_x_zz[i] * a_exp;
    }
    // integrals block (654-660)

    #pragma omp simd aligned(g_xxy_xx_y_xx, g_xxy_xx_y_xy, g_xxy_xx_y_xz, g_xxy_xx_y_yy, g_xxy_xx_y_yz, g_xxy_xx_y_zz, g_y_0_0_0_xx_xx_y_xx, g_y_0_0_0_xx_xx_y_xy, g_y_0_0_0_xx_xx_y_xz, g_y_0_0_0_xx_xx_y_yy, g_y_0_0_0_xx_xx_y_yz, g_y_0_0_0_xx_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_xx_y_xx[i] = 2.0 * g_xxy_xx_y_xx[i] * a_exp;

        g_y_0_0_0_xx_xx_y_xy[i] = 2.0 * g_xxy_xx_y_xy[i] * a_exp;

        g_y_0_0_0_xx_xx_y_xz[i] = 2.0 * g_xxy_xx_y_xz[i] * a_exp;

        g_y_0_0_0_xx_xx_y_yy[i] = 2.0 * g_xxy_xx_y_yy[i] * a_exp;

        g_y_0_0_0_xx_xx_y_yz[i] = 2.0 * g_xxy_xx_y_yz[i] * a_exp;

        g_y_0_0_0_xx_xx_y_zz[i] = 2.0 * g_xxy_xx_y_zz[i] * a_exp;
    }
    // integrals block (660-666)

    #pragma omp simd aligned(g_xxy_xx_z_xx, g_xxy_xx_z_xy, g_xxy_xx_z_xz, g_xxy_xx_z_yy, g_xxy_xx_z_yz, g_xxy_xx_z_zz, g_y_0_0_0_xx_xx_z_xx, g_y_0_0_0_xx_xx_z_xy, g_y_0_0_0_xx_xx_z_xz, g_y_0_0_0_xx_xx_z_yy, g_y_0_0_0_xx_xx_z_yz, g_y_0_0_0_xx_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_xx_z_xx[i] = 2.0 * g_xxy_xx_z_xx[i] * a_exp;

        g_y_0_0_0_xx_xx_z_xy[i] = 2.0 * g_xxy_xx_z_xy[i] * a_exp;

        g_y_0_0_0_xx_xx_z_xz[i] = 2.0 * g_xxy_xx_z_xz[i] * a_exp;

        g_y_0_0_0_xx_xx_z_yy[i] = 2.0 * g_xxy_xx_z_yy[i] * a_exp;

        g_y_0_0_0_xx_xx_z_yz[i] = 2.0 * g_xxy_xx_z_yz[i] * a_exp;

        g_y_0_0_0_xx_xx_z_zz[i] = 2.0 * g_xxy_xx_z_zz[i] * a_exp;
    }
    // integrals block (666-672)

    #pragma omp simd aligned(g_xxy_xy_x_xx, g_xxy_xy_x_xy, g_xxy_xy_x_xz, g_xxy_xy_x_yy, g_xxy_xy_x_yz, g_xxy_xy_x_zz, g_y_0_0_0_xx_xy_x_xx, g_y_0_0_0_xx_xy_x_xy, g_y_0_0_0_xx_xy_x_xz, g_y_0_0_0_xx_xy_x_yy, g_y_0_0_0_xx_xy_x_yz, g_y_0_0_0_xx_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_xy_x_xx[i] = 2.0 * g_xxy_xy_x_xx[i] * a_exp;

        g_y_0_0_0_xx_xy_x_xy[i] = 2.0 * g_xxy_xy_x_xy[i] * a_exp;

        g_y_0_0_0_xx_xy_x_xz[i] = 2.0 * g_xxy_xy_x_xz[i] * a_exp;

        g_y_0_0_0_xx_xy_x_yy[i] = 2.0 * g_xxy_xy_x_yy[i] * a_exp;

        g_y_0_0_0_xx_xy_x_yz[i] = 2.0 * g_xxy_xy_x_yz[i] * a_exp;

        g_y_0_0_0_xx_xy_x_zz[i] = 2.0 * g_xxy_xy_x_zz[i] * a_exp;
    }
    // integrals block (672-678)

    #pragma omp simd aligned(g_xxy_xy_y_xx, g_xxy_xy_y_xy, g_xxy_xy_y_xz, g_xxy_xy_y_yy, g_xxy_xy_y_yz, g_xxy_xy_y_zz, g_y_0_0_0_xx_xy_y_xx, g_y_0_0_0_xx_xy_y_xy, g_y_0_0_0_xx_xy_y_xz, g_y_0_0_0_xx_xy_y_yy, g_y_0_0_0_xx_xy_y_yz, g_y_0_0_0_xx_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_xy_y_xx[i] = 2.0 * g_xxy_xy_y_xx[i] * a_exp;

        g_y_0_0_0_xx_xy_y_xy[i] = 2.0 * g_xxy_xy_y_xy[i] * a_exp;

        g_y_0_0_0_xx_xy_y_xz[i] = 2.0 * g_xxy_xy_y_xz[i] * a_exp;

        g_y_0_0_0_xx_xy_y_yy[i] = 2.0 * g_xxy_xy_y_yy[i] * a_exp;

        g_y_0_0_0_xx_xy_y_yz[i] = 2.0 * g_xxy_xy_y_yz[i] * a_exp;

        g_y_0_0_0_xx_xy_y_zz[i] = 2.0 * g_xxy_xy_y_zz[i] * a_exp;
    }
    // integrals block (678-684)

    #pragma omp simd aligned(g_xxy_xy_z_xx, g_xxy_xy_z_xy, g_xxy_xy_z_xz, g_xxy_xy_z_yy, g_xxy_xy_z_yz, g_xxy_xy_z_zz, g_y_0_0_0_xx_xy_z_xx, g_y_0_0_0_xx_xy_z_xy, g_y_0_0_0_xx_xy_z_xz, g_y_0_0_0_xx_xy_z_yy, g_y_0_0_0_xx_xy_z_yz, g_y_0_0_0_xx_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_xy_z_xx[i] = 2.0 * g_xxy_xy_z_xx[i] * a_exp;

        g_y_0_0_0_xx_xy_z_xy[i] = 2.0 * g_xxy_xy_z_xy[i] * a_exp;

        g_y_0_0_0_xx_xy_z_xz[i] = 2.0 * g_xxy_xy_z_xz[i] * a_exp;

        g_y_0_0_0_xx_xy_z_yy[i] = 2.0 * g_xxy_xy_z_yy[i] * a_exp;

        g_y_0_0_0_xx_xy_z_yz[i] = 2.0 * g_xxy_xy_z_yz[i] * a_exp;

        g_y_0_0_0_xx_xy_z_zz[i] = 2.0 * g_xxy_xy_z_zz[i] * a_exp;
    }
    // integrals block (684-690)

    #pragma omp simd aligned(g_xxy_xz_x_xx, g_xxy_xz_x_xy, g_xxy_xz_x_xz, g_xxy_xz_x_yy, g_xxy_xz_x_yz, g_xxy_xz_x_zz, g_y_0_0_0_xx_xz_x_xx, g_y_0_0_0_xx_xz_x_xy, g_y_0_0_0_xx_xz_x_xz, g_y_0_0_0_xx_xz_x_yy, g_y_0_0_0_xx_xz_x_yz, g_y_0_0_0_xx_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_xz_x_xx[i] = 2.0 * g_xxy_xz_x_xx[i] * a_exp;

        g_y_0_0_0_xx_xz_x_xy[i] = 2.0 * g_xxy_xz_x_xy[i] * a_exp;

        g_y_0_0_0_xx_xz_x_xz[i] = 2.0 * g_xxy_xz_x_xz[i] * a_exp;

        g_y_0_0_0_xx_xz_x_yy[i] = 2.0 * g_xxy_xz_x_yy[i] * a_exp;

        g_y_0_0_0_xx_xz_x_yz[i] = 2.0 * g_xxy_xz_x_yz[i] * a_exp;

        g_y_0_0_0_xx_xz_x_zz[i] = 2.0 * g_xxy_xz_x_zz[i] * a_exp;
    }
    // integrals block (690-696)

    #pragma omp simd aligned(g_xxy_xz_y_xx, g_xxy_xz_y_xy, g_xxy_xz_y_xz, g_xxy_xz_y_yy, g_xxy_xz_y_yz, g_xxy_xz_y_zz, g_y_0_0_0_xx_xz_y_xx, g_y_0_0_0_xx_xz_y_xy, g_y_0_0_0_xx_xz_y_xz, g_y_0_0_0_xx_xz_y_yy, g_y_0_0_0_xx_xz_y_yz, g_y_0_0_0_xx_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_xz_y_xx[i] = 2.0 * g_xxy_xz_y_xx[i] * a_exp;

        g_y_0_0_0_xx_xz_y_xy[i] = 2.0 * g_xxy_xz_y_xy[i] * a_exp;

        g_y_0_0_0_xx_xz_y_xz[i] = 2.0 * g_xxy_xz_y_xz[i] * a_exp;

        g_y_0_0_0_xx_xz_y_yy[i] = 2.0 * g_xxy_xz_y_yy[i] * a_exp;

        g_y_0_0_0_xx_xz_y_yz[i] = 2.0 * g_xxy_xz_y_yz[i] * a_exp;

        g_y_0_0_0_xx_xz_y_zz[i] = 2.0 * g_xxy_xz_y_zz[i] * a_exp;
    }
    // integrals block (696-702)

    #pragma omp simd aligned(g_xxy_xz_z_xx, g_xxy_xz_z_xy, g_xxy_xz_z_xz, g_xxy_xz_z_yy, g_xxy_xz_z_yz, g_xxy_xz_z_zz, g_y_0_0_0_xx_xz_z_xx, g_y_0_0_0_xx_xz_z_xy, g_y_0_0_0_xx_xz_z_xz, g_y_0_0_0_xx_xz_z_yy, g_y_0_0_0_xx_xz_z_yz, g_y_0_0_0_xx_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_xz_z_xx[i] = 2.0 * g_xxy_xz_z_xx[i] * a_exp;

        g_y_0_0_0_xx_xz_z_xy[i] = 2.0 * g_xxy_xz_z_xy[i] * a_exp;

        g_y_0_0_0_xx_xz_z_xz[i] = 2.0 * g_xxy_xz_z_xz[i] * a_exp;

        g_y_0_0_0_xx_xz_z_yy[i] = 2.0 * g_xxy_xz_z_yy[i] * a_exp;

        g_y_0_0_0_xx_xz_z_yz[i] = 2.0 * g_xxy_xz_z_yz[i] * a_exp;

        g_y_0_0_0_xx_xz_z_zz[i] = 2.0 * g_xxy_xz_z_zz[i] * a_exp;
    }
    // integrals block (702-708)

    #pragma omp simd aligned(g_xxy_yy_x_xx, g_xxy_yy_x_xy, g_xxy_yy_x_xz, g_xxy_yy_x_yy, g_xxy_yy_x_yz, g_xxy_yy_x_zz, g_y_0_0_0_xx_yy_x_xx, g_y_0_0_0_xx_yy_x_xy, g_y_0_0_0_xx_yy_x_xz, g_y_0_0_0_xx_yy_x_yy, g_y_0_0_0_xx_yy_x_yz, g_y_0_0_0_xx_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_yy_x_xx[i] = 2.0 * g_xxy_yy_x_xx[i] * a_exp;

        g_y_0_0_0_xx_yy_x_xy[i] = 2.0 * g_xxy_yy_x_xy[i] * a_exp;

        g_y_0_0_0_xx_yy_x_xz[i] = 2.0 * g_xxy_yy_x_xz[i] * a_exp;

        g_y_0_0_0_xx_yy_x_yy[i] = 2.0 * g_xxy_yy_x_yy[i] * a_exp;

        g_y_0_0_0_xx_yy_x_yz[i] = 2.0 * g_xxy_yy_x_yz[i] * a_exp;

        g_y_0_0_0_xx_yy_x_zz[i] = 2.0 * g_xxy_yy_x_zz[i] * a_exp;
    }
    // integrals block (708-714)

    #pragma omp simd aligned(g_xxy_yy_y_xx, g_xxy_yy_y_xy, g_xxy_yy_y_xz, g_xxy_yy_y_yy, g_xxy_yy_y_yz, g_xxy_yy_y_zz, g_y_0_0_0_xx_yy_y_xx, g_y_0_0_0_xx_yy_y_xy, g_y_0_0_0_xx_yy_y_xz, g_y_0_0_0_xx_yy_y_yy, g_y_0_0_0_xx_yy_y_yz, g_y_0_0_0_xx_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_yy_y_xx[i] = 2.0 * g_xxy_yy_y_xx[i] * a_exp;

        g_y_0_0_0_xx_yy_y_xy[i] = 2.0 * g_xxy_yy_y_xy[i] * a_exp;

        g_y_0_0_0_xx_yy_y_xz[i] = 2.0 * g_xxy_yy_y_xz[i] * a_exp;

        g_y_0_0_0_xx_yy_y_yy[i] = 2.0 * g_xxy_yy_y_yy[i] * a_exp;

        g_y_0_0_0_xx_yy_y_yz[i] = 2.0 * g_xxy_yy_y_yz[i] * a_exp;

        g_y_0_0_0_xx_yy_y_zz[i] = 2.0 * g_xxy_yy_y_zz[i] * a_exp;
    }
    // integrals block (714-720)

    #pragma omp simd aligned(g_xxy_yy_z_xx, g_xxy_yy_z_xy, g_xxy_yy_z_xz, g_xxy_yy_z_yy, g_xxy_yy_z_yz, g_xxy_yy_z_zz, g_y_0_0_0_xx_yy_z_xx, g_y_0_0_0_xx_yy_z_xy, g_y_0_0_0_xx_yy_z_xz, g_y_0_0_0_xx_yy_z_yy, g_y_0_0_0_xx_yy_z_yz, g_y_0_0_0_xx_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_yy_z_xx[i] = 2.0 * g_xxy_yy_z_xx[i] * a_exp;

        g_y_0_0_0_xx_yy_z_xy[i] = 2.0 * g_xxy_yy_z_xy[i] * a_exp;

        g_y_0_0_0_xx_yy_z_xz[i] = 2.0 * g_xxy_yy_z_xz[i] * a_exp;

        g_y_0_0_0_xx_yy_z_yy[i] = 2.0 * g_xxy_yy_z_yy[i] * a_exp;

        g_y_0_0_0_xx_yy_z_yz[i] = 2.0 * g_xxy_yy_z_yz[i] * a_exp;

        g_y_0_0_0_xx_yy_z_zz[i] = 2.0 * g_xxy_yy_z_zz[i] * a_exp;
    }
    // integrals block (720-726)

    #pragma omp simd aligned(g_xxy_yz_x_xx, g_xxy_yz_x_xy, g_xxy_yz_x_xz, g_xxy_yz_x_yy, g_xxy_yz_x_yz, g_xxy_yz_x_zz, g_y_0_0_0_xx_yz_x_xx, g_y_0_0_0_xx_yz_x_xy, g_y_0_0_0_xx_yz_x_xz, g_y_0_0_0_xx_yz_x_yy, g_y_0_0_0_xx_yz_x_yz, g_y_0_0_0_xx_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_yz_x_xx[i] = 2.0 * g_xxy_yz_x_xx[i] * a_exp;

        g_y_0_0_0_xx_yz_x_xy[i] = 2.0 * g_xxy_yz_x_xy[i] * a_exp;

        g_y_0_0_0_xx_yz_x_xz[i] = 2.0 * g_xxy_yz_x_xz[i] * a_exp;

        g_y_0_0_0_xx_yz_x_yy[i] = 2.0 * g_xxy_yz_x_yy[i] * a_exp;

        g_y_0_0_0_xx_yz_x_yz[i] = 2.0 * g_xxy_yz_x_yz[i] * a_exp;

        g_y_0_0_0_xx_yz_x_zz[i] = 2.0 * g_xxy_yz_x_zz[i] * a_exp;
    }
    // integrals block (726-732)

    #pragma omp simd aligned(g_xxy_yz_y_xx, g_xxy_yz_y_xy, g_xxy_yz_y_xz, g_xxy_yz_y_yy, g_xxy_yz_y_yz, g_xxy_yz_y_zz, g_y_0_0_0_xx_yz_y_xx, g_y_0_0_0_xx_yz_y_xy, g_y_0_0_0_xx_yz_y_xz, g_y_0_0_0_xx_yz_y_yy, g_y_0_0_0_xx_yz_y_yz, g_y_0_0_0_xx_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_yz_y_xx[i] = 2.0 * g_xxy_yz_y_xx[i] * a_exp;

        g_y_0_0_0_xx_yz_y_xy[i] = 2.0 * g_xxy_yz_y_xy[i] * a_exp;

        g_y_0_0_0_xx_yz_y_xz[i] = 2.0 * g_xxy_yz_y_xz[i] * a_exp;

        g_y_0_0_0_xx_yz_y_yy[i] = 2.0 * g_xxy_yz_y_yy[i] * a_exp;

        g_y_0_0_0_xx_yz_y_yz[i] = 2.0 * g_xxy_yz_y_yz[i] * a_exp;

        g_y_0_0_0_xx_yz_y_zz[i] = 2.0 * g_xxy_yz_y_zz[i] * a_exp;
    }
    // integrals block (732-738)

    #pragma omp simd aligned(g_xxy_yz_z_xx, g_xxy_yz_z_xy, g_xxy_yz_z_xz, g_xxy_yz_z_yy, g_xxy_yz_z_yz, g_xxy_yz_z_zz, g_y_0_0_0_xx_yz_z_xx, g_y_0_0_0_xx_yz_z_xy, g_y_0_0_0_xx_yz_z_xz, g_y_0_0_0_xx_yz_z_yy, g_y_0_0_0_xx_yz_z_yz, g_y_0_0_0_xx_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_yz_z_xx[i] = 2.0 * g_xxy_yz_z_xx[i] * a_exp;

        g_y_0_0_0_xx_yz_z_xy[i] = 2.0 * g_xxy_yz_z_xy[i] * a_exp;

        g_y_0_0_0_xx_yz_z_xz[i] = 2.0 * g_xxy_yz_z_xz[i] * a_exp;

        g_y_0_0_0_xx_yz_z_yy[i] = 2.0 * g_xxy_yz_z_yy[i] * a_exp;

        g_y_0_0_0_xx_yz_z_yz[i] = 2.0 * g_xxy_yz_z_yz[i] * a_exp;

        g_y_0_0_0_xx_yz_z_zz[i] = 2.0 * g_xxy_yz_z_zz[i] * a_exp;
    }
    // integrals block (738-744)

    #pragma omp simd aligned(g_xxy_zz_x_xx, g_xxy_zz_x_xy, g_xxy_zz_x_xz, g_xxy_zz_x_yy, g_xxy_zz_x_yz, g_xxy_zz_x_zz, g_y_0_0_0_xx_zz_x_xx, g_y_0_0_0_xx_zz_x_xy, g_y_0_0_0_xx_zz_x_xz, g_y_0_0_0_xx_zz_x_yy, g_y_0_0_0_xx_zz_x_yz, g_y_0_0_0_xx_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_zz_x_xx[i] = 2.0 * g_xxy_zz_x_xx[i] * a_exp;

        g_y_0_0_0_xx_zz_x_xy[i] = 2.0 * g_xxy_zz_x_xy[i] * a_exp;

        g_y_0_0_0_xx_zz_x_xz[i] = 2.0 * g_xxy_zz_x_xz[i] * a_exp;

        g_y_0_0_0_xx_zz_x_yy[i] = 2.0 * g_xxy_zz_x_yy[i] * a_exp;

        g_y_0_0_0_xx_zz_x_yz[i] = 2.0 * g_xxy_zz_x_yz[i] * a_exp;

        g_y_0_0_0_xx_zz_x_zz[i] = 2.0 * g_xxy_zz_x_zz[i] * a_exp;
    }
    // integrals block (744-750)

    #pragma omp simd aligned(g_xxy_zz_y_xx, g_xxy_zz_y_xy, g_xxy_zz_y_xz, g_xxy_zz_y_yy, g_xxy_zz_y_yz, g_xxy_zz_y_zz, g_y_0_0_0_xx_zz_y_xx, g_y_0_0_0_xx_zz_y_xy, g_y_0_0_0_xx_zz_y_xz, g_y_0_0_0_xx_zz_y_yy, g_y_0_0_0_xx_zz_y_yz, g_y_0_0_0_xx_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_zz_y_xx[i] = 2.0 * g_xxy_zz_y_xx[i] * a_exp;

        g_y_0_0_0_xx_zz_y_xy[i] = 2.0 * g_xxy_zz_y_xy[i] * a_exp;

        g_y_0_0_0_xx_zz_y_xz[i] = 2.0 * g_xxy_zz_y_xz[i] * a_exp;

        g_y_0_0_0_xx_zz_y_yy[i] = 2.0 * g_xxy_zz_y_yy[i] * a_exp;

        g_y_0_0_0_xx_zz_y_yz[i] = 2.0 * g_xxy_zz_y_yz[i] * a_exp;

        g_y_0_0_0_xx_zz_y_zz[i] = 2.0 * g_xxy_zz_y_zz[i] * a_exp;
    }
    // integrals block (750-756)

    #pragma omp simd aligned(g_xxy_zz_z_xx, g_xxy_zz_z_xy, g_xxy_zz_z_xz, g_xxy_zz_z_yy, g_xxy_zz_z_yz, g_xxy_zz_z_zz, g_y_0_0_0_xx_zz_z_xx, g_y_0_0_0_xx_zz_z_xy, g_y_0_0_0_xx_zz_z_xz, g_y_0_0_0_xx_zz_z_yy, g_y_0_0_0_xx_zz_z_yz, g_y_0_0_0_xx_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_zz_z_xx[i] = 2.0 * g_xxy_zz_z_xx[i] * a_exp;

        g_y_0_0_0_xx_zz_z_xy[i] = 2.0 * g_xxy_zz_z_xy[i] * a_exp;

        g_y_0_0_0_xx_zz_z_xz[i] = 2.0 * g_xxy_zz_z_xz[i] * a_exp;

        g_y_0_0_0_xx_zz_z_yy[i] = 2.0 * g_xxy_zz_z_yy[i] * a_exp;

        g_y_0_0_0_xx_zz_z_yz[i] = 2.0 * g_xxy_zz_z_yz[i] * a_exp;

        g_y_0_0_0_xx_zz_z_zz[i] = 2.0 * g_xxy_zz_z_zz[i] * a_exp;
    }
    // integrals block (756-762)

    #pragma omp simd aligned(g_x_xx_x_xx, g_x_xx_x_xy, g_x_xx_x_xz, g_x_xx_x_yy, g_x_xx_x_yz, g_x_xx_x_zz, g_xyy_xx_x_xx, g_xyy_xx_x_xy, g_xyy_xx_x_xz, g_xyy_xx_x_yy, g_xyy_xx_x_yz, g_xyy_xx_x_zz, g_y_0_0_0_xy_xx_x_xx, g_y_0_0_0_xy_xx_x_xy, g_y_0_0_0_xy_xx_x_xz, g_y_0_0_0_xy_xx_x_yy, g_y_0_0_0_xy_xx_x_yz, g_y_0_0_0_xy_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_xx_x_xx[i] = -g_x_xx_x_xx[i] + 2.0 * g_xyy_xx_x_xx[i] * a_exp;

        g_y_0_0_0_xy_xx_x_xy[i] = -g_x_xx_x_xy[i] + 2.0 * g_xyy_xx_x_xy[i] * a_exp;

        g_y_0_0_0_xy_xx_x_xz[i] = -g_x_xx_x_xz[i] + 2.0 * g_xyy_xx_x_xz[i] * a_exp;

        g_y_0_0_0_xy_xx_x_yy[i] = -g_x_xx_x_yy[i] + 2.0 * g_xyy_xx_x_yy[i] * a_exp;

        g_y_0_0_0_xy_xx_x_yz[i] = -g_x_xx_x_yz[i] + 2.0 * g_xyy_xx_x_yz[i] * a_exp;

        g_y_0_0_0_xy_xx_x_zz[i] = -g_x_xx_x_zz[i] + 2.0 * g_xyy_xx_x_zz[i] * a_exp;
    }
    // integrals block (762-768)

    #pragma omp simd aligned(g_x_xx_y_xx, g_x_xx_y_xy, g_x_xx_y_xz, g_x_xx_y_yy, g_x_xx_y_yz, g_x_xx_y_zz, g_xyy_xx_y_xx, g_xyy_xx_y_xy, g_xyy_xx_y_xz, g_xyy_xx_y_yy, g_xyy_xx_y_yz, g_xyy_xx_y_zz, g_y_0_0_0_xy_xx_y_xx, g_y_0_0_0_xy_xx_y_xy, g_y_0_0_0_xy_xx_y_xz, g_y_0_0_0_xy_xx_y_yy, g_y_0_0_0_xy_xx_y_yz, g_y_0_0_0_xy_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_xx_y_xx[i] = -g_x_xx_y_xx[i] + 2.0 * g_xyy_xx_y_xx[i] * a_exp;

        g_y_0_0_0_xy_xx_y_xy[i] = -g_x_xx_y_xy[i] + 2.0 * g_xyy_xx_y_xy[i] * a_exp;

        g_y_0_0_0_xy_xx_y_xz[i] = -g_x_xx_y_xz[i] + 2.0 * g_xyy_xx_y_xz[i] * a_exp;

        g_y_0_0_0_xy_xx_y_yy[i] = -g_x_xx_y_yy[i] + 2.0 * g_xyy_xx_y_yy[i] * a_exp;

        g_y_0_0_0_xy_xx_y_yz[i] = -g_x_xx_y_yz[i] + 2.0 * g_xyy_xx_y_yz[i] * a_exp;

        g_y_0_0_0_xy_xx_y_zz[i] = -g_x_xx_y_zz[i] + 2.0 * g_xyy_xx_y_zz[i] * a_exp;
    }
    // integrals block (768-774)

    #pragma omp simd aligned(g_x_xx_z_xx, g_x_xx_z_xy, g_x_xx_z_xz, g_x_xx_z_yy, g_x_xx_z_yz, g_x_xx_z_zz, g_xyy_xx_z_xx, g_xyy_xx_z_xy, g_xyy_xx_z_xz, g_xyy_xx_z_yy, g_xyy_xx_z_yz, g_xyy_xx_z_zz, g_y_0_0_0_xy_xx_z_xx, g_y_0_0_0_xy_xx_z_xy, g_y_0_0_0_xy_xx_z_xz, g_y_0_0_0_xy_xx_z_yy, g_y_0_0_0_xy_xx_z_yz, g_y_0_0_0_xy_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_xx_z_xx[i] = -g_x_xx_z_xx[i] + 2.0 * g_xyy_xx_z_xx[i] * a_exp;

        g_y_0_0_0_xy_xx_z_xy[i] = -g_x_xx_z_xy[i] + 2.0 * g_xyy_xx_z_xy[i] * a_exp;

        g_y_0_0_0_xy_xx_z_xz[i] = -g_x_xx_z_xz[i] + 2.0 * g_xyy_xx_z_xz[i] * a_exp;

        g_y_0_0_0_xy_xx_z_yy[i] = -g_x_xx_z_yy[i] + 2.0 * g_xyy_xx_z_yy[i] * a_exp;

        g_y_0_0_0_xy_xx_z_yz[i] = -g_x_xx_z_yz[i] + 2.0 * g_xyy_xx_z_yz[i] * a_exp;

        g_y_0_0_0_xy_xx_z_zz[i] = -g_x_xx_z_zz[i] + 2.0 * g_xyy_xx_z_zz[i] * a_exp;
    }
    // integrals block (774-780)

    #pragma omp simd aligned(g_x_xy_x_xx, g_x_xy_x_xy, g_x_xy_x_xz, g_x_xy_x_yy, g_x_xy_x_yz, g_x_xy_x_zz, g_xyy_xy_x_xx, g_xyy_xy_x_xy, g_xyy_xy_x_xz, g_xyy_xy_x_yy, g_xyy_xy_x_yz, g_xyy_xy_x_zz, g_y_0_0_0_xy_xy_x_xx, g_y_0_0_0_xy_xy_x_xy, g_y_0_0_0_xy_xy_x_xz, g_y_0_0_0_xy_xy_x_yy, g_y_0_0_0_xy_xy_x_yz, g_y_0_0_0_xy_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_xy_x_xx[i] = -g_x_xy_x_xx[i] + 2.0 * g_xyy_xy_x_xx[i] * a_exp;

        g_y_0_0_0_xy_xy_x_xy[i] = -g_x_xy_x_xy[i] + 2.0 * g_xyy_xy_x_xy[i] * a_exp;

        g_y_0_0_0_xy_xy_x_xz[i] = -g_x_xy_x_xz[i] + 2.0 * g_xyy_xy_x_xz[i] * a_exp;

        g_y_0_0_0_xy_xy_x_yy[i] = -g_x_xy_x_yy[i] + 2.0 * g_xyy_xy_x_yy[i] * a_exp;

        g_y_0_0_0_xy_xy_x_yz[i] = -g_x_xy_x_yz[i] + 2.0 * g_xyy_xy_x_yz[i] * a_exp;

        g_y_0_0_0_xy_xy_x_zz[i] = -g_x_xy_x_zz[i] + 2.0 * g_xyy_xy_x_zz[i] * a_exp;
    }
    // integrals block (780-786)

    #pragma omp simd aligned(g_x_xy_y_xx, g_x_xy_y_xy, g_x_xy_y_xz, g_x_xy_y_yy, g_x_xy_y_yz, g_x_xy_y_zz, g_xyy_xy_y_xx, g_xyy_xy_y_xy, g_xyy_xy_y_xz, g_xyy_xy_y_yy, g_xyy_xy_y_yz, g_xyy_xy_y_zz, g_y_0_0_0_xy_xy_y_xx, g_y_0_0_0_xy_xy_y_xy, g_y_0_0_0_xy_xy_y_xz, g_y_0_0_0_xy_xy_y_yy, g_y_0_0_0_xy_xy_y_yz, g_y_0_0_0_xy_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_xy_y_xx[i] = -g_x_xy_y_xx[i] + 2.0 * g_xyy_xy_y_xx[i] * a_exp;

        g_y_0_0_0_xy_xy_y_xy[i] = -g_x_xy_y_xy[i] + 2.0 * g_xyy_xy_y_xy[i] * a_exp;

        g_y_0_0_0_xy_xy_y_xz[i] = -g_x_xy_y_xz[i] + 2.0 * g_xyy_xy_y_xz[i] * a_exp;

        g_y_0_0_0_xy_xy_y_yy[i] = -g_x_xy_y_yy[i] + 2.0 * g_xyy_xy_y_yy[i] * a_exp;

        g_y_0_0_0_xy_xy_y_yz[i] = -g_x_xy_y_yz[i] + 2.0 * g_xyy_xy_y_yz[i] * a_exp;

        g_y_0_0_0_xy_xy_y_zz[i] = -g_x_xy_y_zz[i] + 2.0 * g_xyy_xy_y_zz[i] * a_exp;
    }
    // integrals block (786-792)

    #pragma omp simd aligned(g_x_xy_z_xx, g_x_xy_z_xy, g_x_xy_z_xz, g_x_xy_z_yy, g_x_xy_z_yz, g_x_xy_z_zz, g_xyy_xy_z_xx, g_xyy_xy_z_xy, g_xyy_xy_z_xz, g_xyy_xy_z_yy, g_xyy_xy_z_yz, g_xyy_xy_z_zz, g_y_0_0_0_xy_xy_z_xx, g_y_0_0_0_xy_xy_z_xy, g_y_0_0_0_xy_xy_z_xz, g_y_0_0_0_xy_xy_z_yy, g_y_0_0_0_xy_xy_z_yz, g_y_0_0_0_xy_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_xy_z_xx[i] = -g_x_xy_z_xx[i] + 2.0 * g_xyy_xy_z_xx[i] * a_exp;

        g_y_0_0_0_xy_xy_z_xy[i] = -g_x_xy_z_xy[i] + 2.0 * g_xyy_xy_z_xy[i] * a_exp;

        g_y_0_0_0_xy_xy_z_xz[i] = -g_x_xy_z_xz[i] + 2.0 * g_xyy_xy_z_xz[i] * a_exp;

        g_y_0_0_0_xy_xy_z_yy[i] = -g_x_xy_z_yy[i] + 2.0 * g_xyy_xy_z_yy[i] * a_exp;

        g_y_0_0_0_xy_xy_z_yz[i] = -g_x_xy_z_yz[i] + 2.0 * g_xyy_xy_z_yz[i] * a_exp;

        g_y_0_0_0_xy_xy_z_zz[i] = -g_x_xy_z_zz[i] + 2.0 * g_xyy_xy_z_zz[i] * a_exp;
    }
    // integrals block (792-798)

    #pragma omp simd aligned(g_x_xz_x_xx, g_x_xz_x_xy, g_x_xz_x_xz, g_x_xz_x_yy, g_x_xz_x_yz, g_x_xz_x_zz, g_xyy_xz_x_xx, g_xyy_xz_x_xy, g_xyy_xz_x_xz, g_xyy_xz_x_yy, g_xyy_xz_x_yz, g_xyy_xz_x_zz, g_y_0_0_0_xy_xz_x_xx, g_y_0_0_0_xy_xz_x_xy, g_y_0_0_0_xy_xz_x_xz, g_y_0_0_0_xy_xz_x_yy, g_y_0_0_0_xy_xz_x_yz, g_y_0_0_0_xy_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_xz_x_xx[i] = -g_x_xz_x_xx[i] + 2.0 * g_xyy_xz_x_xx[i] * a_exp;

        g_y_0_0_0_xy_xz_x_xy[i] = -g_x_xz_x_xy[i] + 2.0 * g_xyy_xz_x_xy[i] * a_exp;

        g_y_0_0_0_xy_xz_x_xz[i] = -g_x_xz_x_xz[i] + 2.0 * g_xyy_xz_x_xz[i] * a_exp;

        g_y_0_0_0_xy_xz_x_yy[i] = -g_x_xz_x_yy[i] + 2.0 * g_xyy_xz_x_yy[i] * a_exp;

        g_y_0_0_0_xy_xz_x_yz[i] = -g_x_xz_x_yz[i] + 2.0 * g_xyy_xz_x_yz[i] * a_exp;

        g_y_0_0_0_xy_xz_x_zz[i] = -g_x_xz_x_zz[i] + 2.0 * g_xyy_xz_x_zz[i] * a_exp;
    }
    // integrals block (798-804)

    #pragma omp simd aligned(g_x_xz_y_xx, g_x_xz_y_xy, g_x_xz_y_xz, g_x_xz_y_yy, g_x_xz_y_yz, g_x_xz_y_zz, g_xyy_xz_y_xx, g_xyy_xz_y_xy, g_xyy_xz_y_xz, g_xyy_xz_y_yy, g_xyy_xz_y_yz, g_xyy_xz_y_zz, g_y_0_0_0_xy_xz_y_xx, g_y_0_0_0_xy_xz_y_xy, g_y_0_0_0_xy_xz_y_xz, g_y_0_0_0_xy_xz_y_yy, g_y_0_0_0_xy_xz_y_yz, g_y_0_0_0_xy_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_xz_y_xx[i] = -g_x_xz_y_xx[i] + 2.0 * g_xyy_xz_y_xx[i] * a_exp;

        g_y_0_0_0_xy_xz_y_xy[i] = -g_x_xz_y_xy[i] + 2.0 * g_xyy_xz_y_xy[i] * a_exp;

        g_y_0_0_0_xy_xz_y_xz[i] = -g_x_xz_y_xz[i] + 2.0 * g_xyy_xz_y_xz[i] * a_exp;

        g_y_0_0_0_xy_xz_y_yy[i] = -g_x_xz_y_yy[i] + 2.0 * g_xyy_xz_y_yy[i] * a_exp;

        g_y_0_0_0_xy_xz_y_yz[i] = -g_x_xz_y_yz[i] + 2.0 * g_xyy_xz_y_yz[i] * a_exp;

        g_y_0_0_0_xy_xz_y_zz[i] = -g_x_xz_y_zz[i] + 2.0 * g_xyy_xz_y_zz[i] * a_exp;
    }
    // integrals block (804-810)

    #pragma omp simd aligned(g_x_xz_z_xx, g_x_xz_z_xy, g_x_xz_z_xz, g_x_xz_z_yy, g_x_xz_z_yz, g_x_xz_z_zz, g_xyy_xz_z_xx, g_xyy_xz_z_xy, g_xyy_xz_z_xz, g_xyy_xz_z_yy, g_xyy_xz_z_yz, g_xyy_xz_z_zz, g_y_0_0_0_xy_xz_z_xx, g_y_0_0_0_xy_xz_z_xy, g_y_0_0_0_xy_xz_z_xz, g_y_0_0_0_xy_xz_z_yy, g_y_0_0_0_xy_xz_z_yz, g_y_0_0_0_xy_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_xz_z_xx[i] = -g_x_xz_z_xx[i] + 2.0 * g_xyy_xz_z_xx[i] * a_exp;

        g_y_0_0_0_xy_xz_z_xy[i] = -g_x_xz_z_xy[i] + 2.0 * g_xyy_xz_z_xy[i] * a_exp;

        g_y_0_0_0_xy_xz_z_xz[i] = -g_x_xz_z_xz[i] + 2.0 * g_xyy_xz_z_xz[i] * a_exp;

        g_y_0_0_0_xy_xz_z_yy[i] = -g_x_xz_z_yy[i] + 2.0 * g_xyy_xz_z_yy[i] * a_exp;

        g_y_0_0_0_xy_xz_z_yz[i] = -g_x_xz_z_yz[i] + 2.0 * g_xyy_xz_z_yz[i] * a_exp;

        g_y_0_0_0_xy_xz_z_zz[i] = -g_x_xz_z_zz[i] + 2.0 * g_xyy_xz_z_zz[i] * a_exp;
    }
    // integrals block (810-816)

    #pragma omp simd aligned(g_x_yy_x_xx, g_x_yy_x_xy, g_x_yy_x_xz, g_x_yy_x_yy, g_x_yy_x_yz, g_x_yy_x_zz, g_xyy_yy_x_xx, g_xyy_yy_x_xy, g_xyy_yy_x_xz, g_xyy_yy_x_yy, g_xyy_yy_x_yz, g_xyy_yy_x_zz, g_y_0_0_0_xy_yy_x_xx, g_y_0_0_0_xy_yy_x_xy, g_y_0_0_0_xy_yy_x_xz, g_y_0_0_0_xy_yy_x_yy, g_y_0_0_0_xy_yy_x_yz, g_y_0_0_0_xy_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_yy_x_xx[i] = -g_x_yy_x_xx[i] + 2.0 * g_xyy_yy_x_xx[i] * a_exp;

        g_y_0_0_0_xy_yy_x_xy[i] = -g_x_yy_x_xy[i] + 2.0 * g_xyy_yy_x_xy[i] * a_exp;

        g_y_0_0_0_xy_yy_x_xz[i] = -g_x_yy_x_xz[i] + 2.0 * g_xyy_yy_x_xz[i] * a_exp;

        g_y_0_0_0_xy_yy_x_yy[i] = -g_x_yy_x_yy[i] + 2.0 * g_xyy_yy_x_yy[i] * a_exp;

        g_y_0_0_0_xy_yy_x_yz[i] = -g_x_yy_x_yz[i] + 2.0 * g_xyy_yy_x_yz[i] * a_exp;

        g_y_0_0_0_xy_yy_x_zz[i] = -g_x_yy_x_zz[i] + 2.0 * g_xyy_yy_x_zz[i] * a_exp;
    }
    // integrals block (816-822)

    #pragma omp simd aligned(g_x_yy_y_xx, g_x_yy_y_xy, g_x_yy_y_xz, g_x_yy_y_yy, g_x_yy_y_yz, g_x_yy_y_zz, g_xyy_yy_y_xx, g_xyy_yy_y_xy, g_xyy_yy_y_xz, g_xyy_yy_y_yy, g_xyy_yy_y_yz, g_xyy_yy_y_zz, g_y_0_0_0_xy_yy_y_xx, g_y_0_0_0_xy_yy_y_xy, g_y_0_0_0_xy_yy_y_xz, g_y_0_0_0_xy_yy_y_yy, g_y_0_0_0_xy_yy_y_yz, g_y_0_0_0_xy_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_yy_y_xx[i] = -g_x_yy_y_xx[i] + 2.0 * g_xyy_yy_y_xx[i] * a_exp;

        g_y_0_0_0_xy_yy_y_xy[i] = -g_x_yy_y_xy[i] + 2.0 * g_xyy_yy_y_xy[i] * a_exp;

        g_y_0_0_0_xy_yy_y_xz[i] = -g_x_yy_y_xz[i] + 2.0 * g_xyy_yy_y_xz[i] * a_exp;

        g_y_0_0_0_xy_yy_y_yy[i] = -g_x_yy_y_yy[i] + 2.0 * g_xyy_yy_y_yy[i] * a_exp;

        g_y_0_0_0_xy_yy_y_yz[i] = -g_x_yy_y_yz[i] + 2.0 * g_xyy_yy_y_yz[i] * a_exp;

        g_y_0_0_0_xy_yy_y_zz[i] = -g_x_yy_y_zz[i] + 2.0 * g_xyy_yy_y_zz[i] * a_exp;
    }
    // integrals block (822-828)

    #pragma omp simd aligned(g_x_yy_z_xx, g_x_yy_z_xy, g_x_yy_z_xz, g_x_yy_z_yy, g_x_yy_z_yz, g_x_yy_z_zz, g_xyy_yy_z_xx, g_xyy_yy_z_xy, g_xyy_yy_z_xz, g_xyy_yy_z_yy, g_xyy_yy_z_yz, g_xyy_yy_z_zz, g_y_0_0_0_xy_yy_z_xx, g_y_0_0_0_xy_yy_z_xy, g_y_0_0_0_xy_yy_z_xz, g_y_0_0_0_xy_yy_z_yy, g_y_0_0_0_xy_yy_z_yz, g_y_0_0_0_xy_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_yy_z_xx[i] = -g_x_yy_z_xx[i] + 2.0 * g_xyy_yy_z_xx[i] * a_exp;

        g_y_0_0_0_xy_yy_z_xy[i] = -g_x_yy_z_xy[i] + 2.0 * g_xyy_yy_z_xy[i] * a_exp;

        g_y_0_0_0_xy_yy_z_xz[i] = -g_x_yy_z_xz[i] + 2.0 * g_xyy_yy_z_xz[i] * a_exp;

        g_y_0_0_0_xy_yy_z_yy[i] = -g_x_yy_z_yy[i] + 2.0 * g_xyy_yy_z_yy[i] * a_exp;

        g_y_0_0_0_xy_yy_z_yz[i] = -g_x_yy_z_yz[i] + 2.0 * g_xyy_yy_z_yz[i] * a_exp;

        g_y_0_0_0_xy_yy_z_zz[i] = -g_x_yy_z_zz[i] + 2.0 * g_xyy_yy_z_zz[i] * a_exp;
    }
    // integrals block (828-834)

    #pragma omp simd aligned(g_x_yz_x_xx, g_x_yz_x_xy, g_x_yz_x_xz, g_x_yz_x_yy, g_x_yz_x_yz, g_x_yz_x_zz, g_xyy_yz_x_xx, g_xyy_yz_x_xy, g_xyy_yz_x_xz, g_xyy_yz_x_yy, g_xyy_yz_x_yz, g_xyy_yz_x_zz, g_y_0_0_0_xy_yz_x_xx, g_y_0_0_0_xy_yz_x_xy, g_y_0_0_0_xy_yz_x_xz, g_y_0_0_0_xy_yz_x_yy, g_y_0_0_0_xy_yz_x_yz, g_y_0_0_0_xy_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_yz_x_xx[i] = -g_x_yz_x_xx[i] + 2.0 * g_xyy_yz_x_xx[i] * a_exp;

        g_y_0_0_0_xy_yz_x_xy[i] = -g_x_yz_x_xy[i] + 2.0 * g_xyy_yz_x_xy[i] * a_exp;

        g_y_0_0_0_xy_yz_x_xz[i] = -g_x_yz_x_xz[i] + 2.0 * g_xyy_yz_x_xz[i] * a_exp;

        g_y_0_0_0_xy_yz_x_yy[i] = -g_x_yz_x_yy[i] + 2.0 * g_xyy_yz_x_yy[i] * a_exp;

        g_y_0_0_0_xy_yz_x_yz[i] = -g_x_yz_x_yz[i] + 2.0 * g_xyy_yz_x_yz[i] * a_exp;

        g_y_0_0_0_xy_yz_x_zz[i] = -g_x_yz_x_zz[i] + 2.0 * g_xyy_yz_x_zz[i] * a_exp;
    }
    // integrals block (834-840)

    #pragma omp simd aligned(g_x_yz_y_xx, g_x_yz_y_xy, g_x_yz_y_xz, g_x_yz_y_yy, g_x_yz_y_yz, g_x_yz_y_zz, g_xyy_yz_y_xx, g_xyy_yz_y_xy, g_xyy_yz_y_xz, g_xyy_yz_y_yy, g_xyy_yz_y_yz, g_xyy_yz_y_zz, g_y_0_0_0_xy_yz_y_xx, g_y_0_0_0_xy_yz_y_xy, g_y_0_0_0_xy_yz_y_xz, g_y_0_0_0_xy_yz_y_yy, g_y_0_0_0_xy_yz_y_yz, g_y_0_0_0_xy_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_yz_y_xx[i] = -g_x_yz_y_xx[i] + 2.0 * g_xyy_yz_y_xx[i] * a_exp;

        g_y_0_0_0_xy_yz_y_xy[i] = -g_x_yz_y_xy[i] + 2.0 * g_xyy_yz_y_xy[i] * a_exp;

        g_y_0_0_0_xy_yz_y_xz[i] = -g_x_yz_y_xz[i] + 2.0 * g_xyy_yz_y_xz[i] * a_exp;

        g_y_0_0_0_xy_yz_y_yy[i] = -g_x_yz_y_yy[i] + 2.0 * g_xyy_yz_y_yy[i] * a_exp;

        g_y_0_0_0_xy_yz_y_yz[i] = -g_x_yz_y_yz[i] + 2.0 * g_xyy_yz_y_yz[i] * a_exp;

        g_y_0_0_0_xy_yz_y_zz[i] = -g_x_yz_y_zz[i] + 2.0 * g_xyy_yz_y_zz[i] * a_exp;
    }
    // integrals block (840-846)

    #pragma omp simd aligned(g_x_yz_z_xx, g_x_yz_z_xy, g_x_yz_z_xz, g_x_yz_z_yy, g_x_yz_z_yz, g_x_yz_z_zz, g_xyy_yz_z_xx, g_xyy_yz_z_xy, g_xyy_yz_z_xz, g_xyy_yz_z_yy, g_xyy_yz_z_yz, g_xyy_yz_z_zz, g_y_0_0_0_xy_yz_z_xx, g_y_0_0_0_xy_yz_z_xy, g_y_0_0_0_xy_yz_z_xz, g_y_0_0_0_xy_yz_z_yy, g_y_0_0_0_xy_yz_z_yz, g_y_0_0_0_xy_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_yz_z_xx[i] = -g_x_yz_z_xx[i] + 2.0 * g_xyy_yz_z_xx[i] * a_exp;

        g_y_0_0_0_xy_yz_z_xy[i] = -g_x_yz_z_xy[i] + 2.0 * g_xyy_yz_z_xy[i] * a_exp;

        g_y_0_0_0_xy_yz_z_xz[i] = -g_x_yz_z_xz[i] + 2.0 * g_xyy_yz_z_xz[i] * a_exp;

        g_y_0_0_0_xy_yz_z_yy[i] = -g_x_yz_z_yy[i] + 2.0 * g_xyy_yz_z_yy[i] * a_exp;

        g_y_0_0_0_xy_yz_z_yz[i] = -g_x_yz_z_yz[i] + 2.0 * g_xyy_yz_z_yz[i] * a_exp;

        g_y_0_0_0_xy_yz_z_zz[i] = -g_x_yz_z_zz[i] + 2.0 * g_xyy_yz_z_zz[i] * a_exp;
    }
    // integrals block (846-852)

    #pragma omp simd aligned(g_x_zz_x_xx, g_x_zz_x_xy, g_x_zz_x_xz, g_x_zz_x_yy, g_x_zz_x_yz, g_x_zz_x_zz, g_xyy_zz_x_xx, g_xyy_zz_x_xy, g_xyy_zz_x_xz, g_xyy_zz_x_yy, g_xyy_zz_x_yz, g_xyy_zz_x_zz, g_y_0_0_0_xy_zz_x_xx, g_y_0_0_0_xy_zz_x_xy, g_y_0_0_0_xy_zz_x_xz, g_y_0_0_0_xy_zz_x_yy, g_y_0_0_0_xy_zz_x_yz, g_y_0_0_0_xy_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_zz_x_xx[i] = -g_x_zz_x_xx[i] + 2.0 * g_xyy_zz_x_xx[i] * a_exp;

        g_y_0_0_0_xy_zz_x_xy[i] = -g_x_zz_x_xy[i] + 2.0 * g_xyy_zz_x_xy[i] * a_exp;

        g_y_0_0_0_xy_zz_x_xz[i] = -g_x_zz_x_xz[i] + 2.0 * g_xyy_zz_x_xz[i] * a_exp;

        g_y_0_0_0_xy_zz_x_yy[i] = -g_x_zz_x_yy[i] + 2.0 * g_xyy_zz_x_yy[i] * a_exp;

        g_y_0_0_0_xy_zz_x_yz[i] = -g_x_zz_x_yz[i] + 2.0 * g_xyy_zz_x_yz[i] * a_exp;

        g_y_0_0_0_xy_zz_x_zz[i] = -g_x_zz_x_zz[i] + 2.0 * g_xyy_zz_x_zz[i] * a_exp;
    }
    // integrals block (852-858)

    #pragma omp simd aligned(g_x_zz_y_xx, g_x_zz_y_xy, g_x_zz_y_xz, g_x_zz_y_yy, g_x_zz_y_yz, g_x_zz_y_zz, g_xyy_zz_y_xx, g_xyy_zz_y_xy, g_xyy_zz_y_xz, g_xyy_zz_y_yy, g_xyy_zz_y_yz, g_xyy_zz_y_zz, g_y_0_0_0_xy_zz_y_xx, g_y_0_0_0_xy_zz_y_xy, g_y_0_0_0_xy_zz_y_xz, g_y_0_0_0_xy_zz_y_yy, g_y_0_0_0_xy_zz_y_yz, g_y_0_0_0_xy_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_zz_y_xx[i] = -g_x_zz_y_xx[i] + 2.0 * g_xyy_zz_y_xx[i] * a_exp;

        g_y_0_0_0_xy_zz_y_xy[i] = -g_x_zz_y_xy[i] + 2.0 * g_xyy_zz_y_xy[i] * a_exp;

        g_y_0_0_0_xy_zz_y_xz[i] = -g_x_zz_y_xz[i] + 2.0 * g_xyy_zz_y_xz[i] * a_exp;

        g_y_0_0_0_xy_zz_y_yy[i] = -g_x_zz_y_yy[i] + 2.0 * g_xyy_zz_y_yy[i] * a_exp;

        g_y_0_0_0_xy_zz_y_yz[i] = -g_x_zz_y_yz[i] + 2.0 * g_xyy_zz_y_yz[i] * a_exp;

        g_y_0_0_0_xy_zz_y_zz[i] = -g_x_zz_y_zz[i] + 2.0 * g_xyy_zz_y_zz[i] * a_exp;
    }
    // integrals block (858-864)

    #pragma omp simd aligned(g_x_zz_z_xx, g_x_zz_z_xy, g_x_zz_z_xz, g_x_zz_z_yy, g_x_zz_z_yz, g_x_zz_z_zz, g_xyy_zz_z_xx, g_xyy_zz_z_xy, g_xyy_zz_z_xz, g_xyy_zz_z_yy, g_xyy_zz_z_yz, g_xyy_zz_z_zz, g_y_0_0_0_xy_zz_z_xx, g_y_0_0_0_xy_zz_z_xy, g_y_0_0_0_xy_zz_z_xz, g_y_0_0_0_xy_zz_z_yy, g_y_0_0_0_xy_zz_z_yz, g_y_0_0_0_xy_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_zz_z_xx[i] = -g_x_zz_z_xx[i] + 2.0 * g_xyy_zz_z_xx[i] * a_exp;

        g_y_0_0_0_xy_zz_z_xy[i] = -g_x_zz_z_xy[i] + 2.0 * g_xyy_zz_z_xy[i] * a_exp;

        g_y_0_0_0_xy_zz_z_xz[i] = -g_x_zz_z_xz[i] + 2.0 * g_xyy_zz_z_xz[i] * a_exp;

        g_y_0_0_0_xy_zz_z_yy[i] = -g_x_zz_z_yy[i] + 2.0 * g_xyy_zz_z_yy[i] * a_exp;

        g_y_0_0_0_xy_zz_z_yz[i] = -g_x_zz_z_yz[i] + 2.0 * g_xyy_zz_z_yz[i] * a_exp;

        g_y_0_0_0_xy_zz_z_zz[i] = -g_x_zz_z_zz[i] + 2.0 * g_xyy_zz_z_zz[i] * a_exp;
    }
    // integrals block (864-870)

    #pragma omp simd aligned(g_xyz_xx_x_xx, g_xyz_xx_x_xy, g_xyz_xx_x_xz, g_xyz_xx_x_yy, g_xyz_xx_x_yz, g_xyz_xx_x_zz, g_y_0_0_0_xz_xx_x_xx, g_y_0_0_0_xz_xx_x_xy, g_y_0_0_0_xz_xx_x_xz, g_y_0_0_0_xz_xx_x_yy, g_y_0_0_0_xz_xx_x_yz, g_y_0_0_0_xz_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_xx_x_xx[i] = 2.0 * g_xyz_xx_x_xx[i] * a_exp;

        g_y_0_0_0_xz_xx_x_xy[i] = 2.0 * g_xyz_xx_x_xy[i] * a_exp;

        g_y_0_0_0_xz_xx_x_xz[i] = 2.0 * g_xyz_xx_x_xz[i] * a_exp;

        g_y_0_0_0_xz_xx_x_yy[i] = 2.0 * g_xyz_xx_x_yy[i] * a_exp;

        g_y_0_0_0_xz_xx_x_yz[i] = 2.0 * g_xyz_xx_x_yz[i] * a_exp;

        g_y_0_0_0_xz_xx_x_zz[i] = 2.0 * g_xyz_xx_x_zz[i] * a_exp;
    }
    // integrals block (870-876)

    #pragma omp simd aligned(g_xyz_xx_y_xx, g_xyz_xx_y_xy, g_xyz_xx_y_xz, g_xyz_xx_y_yy, g_xyz_xx_y_yz, g_xyz_xx_y_zz, g_y_0_0_0_xz_xx_y_xx, g_y_0_0_0_xz_xx_y_xy, g_y_0_0_0_xz_xx_y_xz, g_y_0_0_0_xz_xx_y_yy, g_y_0_0_0_xz_xx_y_yz, g_y_0_0_0_xz_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_xx_y_xx[i] = 2.0 * g_xyz_xx_y_xx[i] * a_exp;

        g_y_0_0_0_xz_xx_y_xy[i] = 2.0 * g_xyz_xx_y_xy[i] * a_exp;

        g_y_0_0_0_xz_xx_y_xz[i] = 2.0 * g_xyz_xx_y_xz[i] * a_exp;

        g_y_0_0_0_xz_xx_y_yy[i] = 2.0 * g_xyz_xx_y_yy[i] * a_exp;

        g_y_0_0_0_xz_xx_y_yz[i] = 2.0 * g_xyz_xx_y_yz[i] * a_exp;

        g_y_0_0_0_xz_xx_y_zz[i] = 2.0 * g_xyz_xx_y_zz[i] * a_exp;
    }
    // integrals block (876-882)

    #pragma omp simd aligned(g_xyz_xx_z_xx, g_xyz_xx_z_xy, g_xyz_xx_z_xz, g_xyz_xx_z_yy, g_xyz_xx_z_yz, g_xyz_xx_z_zz, g_y_0_0_0_xz_xx_z_xx, g_y_0_0_0_xz_xx_z_xy, g_y_0_0_0_xz_xx_z_xz, g_y_0_0_0_xz_xx_z_yy, g_y_0_0_0_xz_xx_z_yz, g_y_0_0_0_xz_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_xx_z_xx[i] = 2.0 * g_xyz_xx_z_xx[i] * a_exp;

        g_y_0_0_0_xz_xx_z_xy[i] = 2.0 * g_xyz_xx_z_xy[i] * a_exp;

        g_y_0_0_0_xz_xx_z_xz[i] = 2.0 * g_xyz_xx_z_xz[i] * a_exp;

        g_y_0_0_0_xz_xx_z_yy[i] = 2.0 * g_xyz_xx_z_yy[i] * a_exp;

        g_y_0_0_0_xz_xx_z_yz[i] = 2.0 * g_xyz_xx_z_yz[i] * a_exp;

        g_y_0_0_0_xz_xx_z_zz[i] = 2.0 * g_xyz_xx_z_zz[i] * a_exp;
    }
    // integrals block (882-888)

    #pragma omp simd aligned(g_xyz_xy_x_xx, g_xyz_xy_x_xy, g_xyz_xy_x_xz, g_xyz_xy_x_yy, g_xyz_xy_x_yz, g_xyz_xy_x_zz, g_y_0_0_0_xz_xy_x_xx, g_y_0_0_0_xz_xy_x_xy, g_y_0_0_0_xz_xy_x_xz, g_y_0_0_0_xz_xy_x_yy, g_y_0_0_0_xz_xy_x_yz, g_y_0_0_0_xz_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_xy_x_xx[i] = 2.0 * g_xyz_xy_x_xx[i] * a_exp;

        g_y_0_0_0_xz_xy_x_xy[i] = 2.0 * g_xyz_xy_x_xy[i] * a_exp;

        g_y_0_0_0_xz_xy_x_xz[i] = 2.0 * g_xyz_xy_x_xz[i] * a_exp;

        g_y_0_0_0_xz_xy_x_yy[i] = 2.0 * g_xyz_xy_x_yy[i] * a_exp;

        g_y_0_0_0_xz_xy_x_yz[i] = 2.0 * g_xyz_xy_x_yz[i] * a_exp;

        g_y_0_0_0_xz_xy_x_zz[i] = 2.0 * g_xyz_xy_x_zz[i] * a_exp;
    }
    // integrals block (888-894)

    #pragma omp simd aligned(g_xyz_xy_y_xx, g_xyz_xy_y_xy, g_xyz_xy_y_xz, g_xyz_xy_y_yy, g_xyz_xy_y_yz, g_xyz_xy_y_zz, g_y_0_0_0_xz_xy_y_xx, g_y_0_0_0_xz_xy_y_xy, g_y_0_0_0_xz_xy_y_xz, g_y_0_0_0_xz_xy_y_yy, g_y_0_0_0_xz_xy_y_yz, g_y_0_0_0_xz_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_xy_y_xx[i] = 2.0 * g_xyz_xy_y_xx[i] * a_exp;

        g_y_0_0_0_xz_xy_y_xy[i] = 2.0 * g_xyz_xy_y_xy[i] * a_exp;

        g_y_0_0_0_xz_xy_y_xz[i] = 2.0 * g_xyz_xy_y_xz[i] * a_exp;

        g_y_0_0_0_xz_xy_y_yy[i] = 2.0 * g_xyz_xy_y_yy[i] * a_exp;

        g_y_0_0_0_xz_xy_y_yz[i] = 2.0 * g_xyz_xy_y_yz[i] * a_exp;

        g_y_0_0_0_xz_xy_y_zz[i] = 2.0 * g_xyz_xy_y_zz[i] * a_exp;
    }
    // integrals block (894-900)

    #pragma omp simd aligned(g_xyz_xy_z_xx, g_xyz_xy_z_xy, g_xyz_xy_z_xz, g_xyz_xy_z_yy, g_xyz_xy_z_yz, g_xyz_xy_z_zz, g_y_0_0_0_xz_xy_z_xx, g_y_0_0_0_xz_xy_z_xy, g_y_0_0_0_xz_xy_z_xz, g_y_0_0_0_xz_xy_z_yy, g_y_0_0_0_xz_xy_z_yz, g_y_0_0_0_xz_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_xy_z_xx[i] = 2.0 * g_xyz_xy_z_xx[i] * a_exp;

        g_y_0_0_0_xz_xy_z_xy[i] = 2.0 * g_xyz_xy_z_xy[i] * a_exp;

        g_y_0_0_0_xz_xy_z_xz[i] = 2.0 * g_xyz_xy_z_xz[i] * a_exp;

        g_y_0_0_0_xz_xy_z_yy[i] = 2.0 * g_xyz_xy_z_yy[i] * a_exp;

        g_y_0_0_0_xz_xy_z_yz[i] = 2.0 * g_xyz_xy_z_yz[i] * a_exp;

        g_y_0_0_0_xz_xy_z_zz[i] = 2.0 * g_xyz_xy_z_zz[i] * a_exp;
    }
    // integrals block (900-906)

    #pragma omp simd aligned(g_xyz_xz_x_xx, g_xyz_xz_x_xy, g_xyz_xz_x_xz, g_xyz_xz_x_yy, g_xyz_xz_x_yz, g_xyz_xz_x_zz, g_y_0_0_0_xz_xz_x_xx, g_y_0_0_0_xz_xz_x_xy, g_y_0_0_0_xz_xz_x_xz, g_y_0_0_0_xz_xz_x_yy, g_y_0_0_0_xz_xz_x_yz, g_y_0_0_0_xz_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_xz_x_xx[i] = 2.0 * g_xyz_xz_x_xx[i] * a_exp;

        g_y_0_0_0_xz_xz_x_xy[i] = 2.0 * g_xyz_xz_x_xy[i] * a_exp;

        g_y_0_0_0_xz_xz_x_xz[i] = 2.0 * g_xyz_xz_x_xz[i] * a_exp;

        g_y_0_0_0_xz_xz_x_yy[i] = 2.0 * g_xyz_xz_x_yy[i] * a_exp;

        g_y_0_0_0_xz_xz_x_yz[i] = 2.0 * g_xyz_xz_x_yz[i] * a_exp;

        g_y_0_0_0_xz_xz_x_zz[i] = 2.0 * g_xyz_xz_x_zz[i] * a_exp;
    }
    // integrals block (906-912)

    #pragma omp simd aligned(g_xyz_xz_y_xx, g_xyz_xz_y_xy, g_xyz_xz_y_xz, g_xyz_xz_y_yy, g_xyz_xz_y_yz, g_xyz_xz_y_zz, g_y_0_0_0_xz_xz_y_xx, g_y_0_0_0_xz_xz_y_xy, g_y_0_0_0_xz_xz_y_xz, g_y_0_0_0_xz_xz_y_yy, g_y_0_0_0_xz_xz_y_yz, g_y_0_0_0_xz_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_xz_y_xx[i] = 2.0 * g_xyz_xz_y_xx[i] * a_exp;

        g_y_0_0_0_xz_xz_y_xy[i] = 2.0 * g_xyz_xz_y_xy[i] * a_exp;

        g_y_0_0_0_xz_xz_y_xz[i] = 2.0 * g_xyz_xz_y_xz[i] * a_exp;

        g_y_0_0_0_xz_xz_y_yy[i] = 2.0 * g_xyz_xz_y_yy[i] * a_exp;

        g_y_0_0_0_xz_xz_y_yz[i] = 2.0 * g_xyz_xz_y_yz[i] * a_exp;

        g_y_0_0_0_xz_xz_y_zz[i] = 2.0 * g_xyz_xz_y_zz[i] * a_exp;
    }
    // integrals block (912-918)

    #pragma omp simd aligned(g_xyz_xz_z_xx, g_xyz_xz_z_xy, g_xyz_xz_z_xz, g_xyz_xz_z_yy, g_xyz_xz_z_yz, g_xyz_xz_z_zz, g_y_0_0_0_xz_xz_z_xx, g_y_0_0_0_xz_xz_z_xy, g_y_0_0_0_xz_xz_z_xz, g_y_0_0_0_xz_xz_z_yy, g_y_0_0_0_xz_xz_z_yz, g_y_0_0_0_xz_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_xz_z_xx[i] = 2.0 * g_xyz_xz_z_xx[i] * a_exp;

        g_y_0_0_0_xz_xz_z_xy[i] = 2.0 * g_xyz_xz_z_xy[i] * a_exp;

        g_y_0_0_0_xz_xz_z_xz[i] = 2.0 * g_xyz_xz_z_xz[i] * a_exp;

        g_y_0_0_0_xz_xz_z_yy[i] = 2.0 * g_xyz_xz_z_yy[i] * a_exp;

        g_y_0_0_0_xz_xz_z_yz[i] = 2.0 * g_xyz_xz_z_yz[i] * a_exp;

        g_y_0_0_0_xz_xz_z_zz[i] = 2.0 * g_xyz_xz_z_zz[i] * a_exp;
    }
    // integrals block (918-924)

    #pragma omp simd aligned(g_xyz_yy_x_xx, g_xyz_yy_x_xy, g_xyz_yy_x_xz, g_xyz_yy_x_yy, g_xyz_yy_x_yz, g_xyz_yy_x_zz, g_y_0_0_0_xz_yy_x_xx, g_y_0_0_0_xz_yy_x_xy, g_y_0_0_0_xz_yy_x_xz, g_y_0_0_0_xz_yy_x_yy, g_y_0_0_0_xz_yy_x_yz, g_y_0_0_0_xz_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_yy_x_xx[i] = 2.0 * g_xyz_yy_x_xx[i] * a_exp;

        g_y_0_0_0_xz_yy_x_xy[i] = 2.0 * g_xyz_yy_x_xy[i] * a_exp;

        g_y_0_0_0_xz_yy_x_xz[i] = 2.0 * g_xyz_yy_x_xz[i] * a_exp;

        g_y_0_0_0_xz_yy_x_yy[i] = 2.0 * g_xyz_yy_x_yy[i] * a_exp;

        g_y_0_0_0_xz_yy_x_yz[i] = 2.0 * g_xyz_yy_x_yz[i] * a_exp;

        g_y_0_0_0_xz_yy_x_zz[i] = 2.0 * g_xyz_yy_x_zz[i] * a_exp;
    }
    // integrals block (924-930)

    #pragma omp simd aligned(g_xyz_yy_y_xx, g_xyz_yy_y_xy, g_xyz_yy_y_xz, g_xyz_yy_y_yy, g_xyz_yy_y_yz, g_xyz_yy_y_zz, g_y_0_0_0_xz_yy_y_xx, g_y_0_0_0_xz_yy_y_xy, g_y_0_0_0_xz_yy_y_xz, g_y_0_0_0_xz_yy_y_yy, g_y_0_0_0_xz_yy_y_yz, g_y_0_0_0_xz_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_yy_y_xx[i] = 2.0 * g_xyz_yy_y_xx[i] * a_exp;

        g_y_0_0_0_xz_yy_y_xy[i] = 2.0 * g_xyz_yy_y_xy[i] * a_exp;

        g_y_0_0_0_xz_yy_y_xz[i] = 2.0 * g_xyz_yy_y_xz[i] * a_exp;

        g_y_0_0_0_xz_yy_y_yy[i] = 2.0 * g_xyz_yy_y_yy[i] * a_exp;

        g_y_0_0_0_xz_yy_y_yz[i] = 2.0 * g_xyz_yy_y_yz[i] * a_exp;

        g_y_0_0_0_xz_yy_y_zz[i] = 2.0 * g_xyz_yy_y_zz[i] * a_exp;
    }
    // integrals block (930-936)

    #pragma omp simd aligned(g_xyz_yy_z_xx, g_xyz_yy_z_xy, g_xyz_yy_z_xz, g_xyz_yy_z_yy, g_xyz_yy_z_yz, g_xyz_yy_z_zz, g_y_0_0_0_xz_yy_z_xx, g_y_0_0_0_xz_yy_z_xy, g_y_0_0_0_xz_yy_z_xz, g_y_0_0_0_xz_yy_z_yy, g_y_0_0_0_xz_yy_z_yz, g_y_0_0_0_xz_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_yy_z_xx[i] = 2.0 * g_xyz_yy_z_xx[i] * a_exp;

        g_y_0_0_0_xz_yy_z_xy[i] = 2.0 * g_xyz_yy_z_xy[i] * a_exp;

        g_y_0_0_0_xz_yy_z_xz[i] = 2.0 * g_xyz_yy_z_xz[i] * a_exp;

        g_y_0_0_0_xz_yy_z_yy[i] = 2.0 * g_xyz_yy_z_yy[i] * a_exp;

        g_y_0_0_0_xz_yy_z_yz[i] = 2.0 * g_xyz_yy_z_yz[i] * a_exp;

        g_y_0_0_0_xz_yy_z_zz[i] = 2.0 * g_xyz_yy_z_zz[i] * a_exp;
    }
    // integrals block (936-942)

    #pragma omp simd aligned(g_xyz_yz_x_xx, g_xyz_yz_x_xy, g_xyz_yz_x_xz, g_xyz_yz_x_yy, g_xyz_yz_x_yz, g_xyz_yz_x_zz, g_y_0_0_0_xz_yz_x_xx, g_y_0_0_0_xz_yz_x_xy, g_y_0_0_0_xz_yz_x_xz, g_y_0_0_0_xz_yz_x_yy, g_y_0_0_0_xz_yz_x_yz, g_y_0_0_0_xz_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_yz_x_xx[i] = 2.0 * g_xyz_yz_x_xx[i] * a_exp;

        g_y_0_0_0_xz_yz_x_xy[i] = 2.0 * g_xyz_yz_x_xy[i] * a_exp;

        g_y_0_0_0_xz_yz_x_xz[i] = 2.0 * g_xyz_yz_x_xz[i] * a_exp;

        g_y_0_0_0_xz_yz_x_yy[i] = 2.0 * g_xyz_yz_x_yy[i] * a_exp;

        g_y_0_0_0_xz_yz_x_yz[i] = 2.0 * g_xyz_yz_x_yz[i] * a_exp;

        g_y_0_0_0_xz_yz_x_zz[i] = 2.0 * g_xyz_yz_x_zz[i] * a_exp;
    }
    // integrals block (942-948)

    #pragma omp simd aligned(g_xyz_yz_y_xx, g_xyz_yz_y_xy, g_xyz_yz_y_xz, g_xyz_yz_y_yy, g_xyz_yz_y_yz, g_xyz_yz_y_zz, g_y_0_0_0_xz_yz_y_xx, g_y_0_0_0_xz_yz_y_xy, g_y_0_0_0_xz_yz_y_xz, g_y_0_0_0_xz_yz_y_yy, g_y_0_0_0_xz_yz_y_yz, g_y_0_0_0_xz_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_yz_y_xx[i] = 2.0 * g_xyz_yz_y_xx[i] * a_exp;

        g_y_0_0_0_xz_yz_y_xy[i] = 2.0 * g_xyz_yz_y_xy[i] * a_exp;

        g_y_0_0_0_xz_yz_y_xz[i] = 2.0 * g_xyz_yz_y_xz[i] * a_exp;

        g_y_0_0_0_xz_yz_y_yy[i] = 2.0 * g_xyz_yz_y_yy[i] * a_exp;

        g_y_0_0_0_xz_yz_y_yz[i] = 2.0 * g_xyz_yz_y_yz[i] * a_exp;

        g_y_0_0_0_xz_yz_y_zz[i] = 2.0 * g_xyz_yz_y_zz[i] * a_exp;
    }
    // integrals block (948-954)

    #pragma omp simd aligned(g_xyz_yz_z_xx, g_xyz_yz_z_xy, g_xyz_yz_z_xz, g_xyz_yz_z_yy, g_xyz_yz_z_yz, g_xyz_yz_z_zz, g_y_0_0_0_xz_yz_z_xx, g_y_0_0_0_xz_yz_z_xy, g_y_0_0_0_xz_yz_z_xz, g_y_0_0_0_xz_yz_z_yy, g_y_0_0_0_xz_yz_z_yz, g_y_0_0_0_xz_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_yz_z_xx[i] = 2.0 * g_xyz_yz_z_xx[i] * a_exp;

        g_y_0_0_0_xz_yz_z_xy[i] = 2.0 * g_xyz_yz_z_xy[i] * a_exp;

        g_y_0_0_0_xz_yz_z_xz[i] = 2.0 * g_xyz_yz_z_xz[i] * a_exp;

        g_y_0_0_0_xz_yz_z_yy[i] = 2.0 * g_xyz_yz_z_yy[i] * a_exp;

        g_y_0_0_0_xz_yz_z_yz[i] = 2.0 * g_xyz_yz_z_yz[i] * a_exp;

        g_y_0_0_0_xz_yz_z_zz[i] = 2.0 * g_xyz_yz_z_zz[i] * a_exp;
    }
    // integrals block (954-960)

    #pragma omp simd aligned(g_xyz_zz_x_xx, g_xyz_zz_x_xy, g_xyz_zz_x_xz, g_xyz_zz_x_yy, g_xyz_zz_x_yz, g_xyz_zz_x_zz, g_y_0_0_0_xz_zz_x_xx, g_y_0_0_0_xz_zz_x_xy, g_y_0_0_0_xz_zz_x_xz, g_y_0_0_0_xz_zz_x_yy, g_y_0_0_0_xz_zz_x_yz, g_y_0_0_0_xz_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_zz_x_xx[i] = 2.0 * g_xyz_zz_x_xx[i] * a_exp;

        g_y_0_0_0_xz_zz_x_xy[i] = 2.0 * g_xyz_zz_x_xy[i] * a_exp;

        g_y_0_0_0_xz_zz_x_xz[i] = 2.0 * g_xyz_zz_x_xz[i] * a_exp;

        g_y_0_0_0_xz_zz_x_yy[i] = 2.0 * g_xyz_zz_x_yy[i] * a_exp;

        g_y_0_0_0_xz_zz_x_yz[i] = 2.0 * g_xyz_zz_x_yz[i] * a_exp;

        g_y_0_0_0_xz_zz_x_zz[i] = 2.0 * g_xyz_zz_x_zz[i] * a_exp;
    }
    // integrals block (960-966)

    #pragma omp simd aligned(g_xyz_zz_y_xx, g_xyz_zz_y_xy, g_xyz_zz_y_xz, g_xyz_zz_y_yy, g_xyz_zz_y_yz, g_xyz_zz_y_zz, g_y_0_0_0_xz_zz_y_xx, g_y_0_0_0_xz_zz_y_xy, g_y_0_0_0_xz_zz_y_xz, g_y_0_0_0_xz_zz_y_yy, g_y_0_0_0_xz_zz_y_yz, g_y_0_0_0_xz_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_zz_y_xx[i] = 2.0 * g_xyz_zz_y_xx[i] * a_exp;

        g_y_0_0_0_xz_zz_y_xy[i] = 2.0 * g_xyz_zz_y_xy[i] * a_exp;

        g_y_0_0_0_xz_zz_y_xz[i] = 2.0 * g_xyz_zz_y_xz[i] * a_exp;

        g_y_0_0_0_xz_zz_y_yy[i] = 2.0 * g_xyz_zz_y_yy[i] * a_exp;

        g_y_0_0_0_xz_zz_y_yz[i] = 2.0 * g_xyz_zz_y_yz[i] * a_exp;

        g_y_0_0_0_xz_zz_y_zz[i] = 2.0 * g_xyz_zz_y_zz[i] * a_exp;
    }
    // integrals block (966-972)

    #pragma omp simd aligned(g_xyz_zz_z_xx, g_xyz_zz_z_xy, g_xyz_zz_z_xz, g_xyz_zz_z_yy, g_xyz_zz_z_yz, g_xyz_zz_z_zz, g_y_0_0_0_xz_zz_z_xx, g_y_0_0_0_xz_zz_z_xy, g_y_0_0_0_xz_zz_z_xz, g_y_0_0_0_xz_zz_z_yy, g_y_0_0_0_xz_zz_z_yz, g_y_0_0_0_xz_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_zz_z_xx[i] = 2.0 * g_xyz_zz_z_xx[i] * a_exp;

        g_y_0_0_0_xz_zz_z_xy[i] = 2.0 * g_xyz_zz_z_xy[i] * a_exp;

        g_y_0_0_0_xz_zz_z_xz[i] = 2.0 * g_xyz_zz_z_xz[i] * a_exp;

        g_y_0_0_0_xz_zz_z_yy[i] = 2.0 * g_xyz_zz_z_yy[i] * a_exp;

        g_y_0_0_0_xz_zz_z_yz[i] = 2.0 * g_xyz_zz_z_yz[i] * a_exp;

        g_y_0_0_0_xz_zz_z_zz[i] = 2.0 * g_xyz_zz_z_zz[i] * a_exp;
    }
    // integrals block (972-978)

    #pragma omp simd aligned(g_y_0_0_0_yy_xx_x_xx, g_y_0_0_0_yy_xx_x_xy, g_y_0_0_0_yy_xx_x_xz, g_y_0_0_0_yy_xx_x_yy, g_y_0_0_0_yy_xx_x_yz, g_y_0_0_0_yy_xx_x_zz, g_y_xx_x_xx, g_y_xx_x_xy, g_y_xx_x_xz, g_y_xx_x_yy, g_y_xx_x_yz, g_y_xx_x_zz, g_yyy_xx_x_xx, g_yyy_xx_x_xy, g_yyy_xx_x_xz, g_yyy_xx_x_yy, g_yyy_xx_x_yz, g_yyy_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_xx_x_xx[i] = -2.0 * g_y_xx_x_xx[i] + 2.0 * g_yyy_xx_x_xx[i] * a_exp;

        g_y_0_0_0_yy_xx_x_xy[i] = -2.0 * g_y_xx_x_xy[i] + 2.0 * g_yyy_xx_x_xy[i] * a_exp;

        g_y_0_0_0_yy_xx_x_xz[i] = -2.0 * g_y_xx_x_xz[i] + 2.0 * g_yyy_xx_x_xz[i] * a_exp;

        g_y_0_0_0_yy_xx_x_yy[i] = -2.0 * g_y_xx_x_yy[i] + 2.0 * g_yyy_xx_x_yy[i] * a_exp;

        g_y_0_0_0_yy_xx_x_yz[i] = -2.0 * g_y_xx_x_yz[i] + 2.0 * g_yyy_xx_x_yz[i] * a_exp;

        g_y_0_0_0_yy_xx_x_zz[i] = -2.0 * g_y_xx_x_zz[i] + 2.0 * g_yyy_xx_x_zz[i] * a_exp;
    }
    // integrals block (978-984)

    #pragma omp simd aligned(g_y_0_0_0_yy_xx_y_xx, g_y_0_0_0_yy_xx_y_xy, g_y_0_0_0_yy_xx_y_xz, g_y_0_0_0_yy_xx_y_yy, g_y_0_0_0_yy_xx_y_yz, g_y_0_0_0_yy_xx_y_zz, g_y_xx_y_xx, g_y_xx_y_xy, g_y_xx_y_xz, g_y_xx_y_yy, g_y_xx_y_yz, g_y_xx_y_zz, g_yyy_xx_y_xx, g_yyy_xx_y_xy, g_yyy_xx_y_xz, g_yyy_xx_y_yy, g_yyy_xx_y_yz, g_yyy_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_xx_y_xx[i] = -2.0 * g_y_xx_y_xx[i] + 2.0 * g_yyy_xx_y_xx[i] * a_exp;

        g_y_0_0_0_yy_xx_y_xy[i] = -2.0 * g_y_xx_y_xy[i] + 2.0 * g_yyy_xx_y_xy[i] * a_exp;

        g_y_0_0_0_yy_xx_y_xz[i] = -2.0 * g_y_xx_y_xz[i] + 2.0 * g_yyy_xx_y_xz[i] * a_exp;

        g_y_0_0_0_yy_xx_y_yy[i] = -2.0 * g_y_xx_y_yy[i] + 2.0 * g_yyy_xx_y_yy[i] * a_exp;

        g_y_0_0_0_yy_xx_y_yz[i] = -2.0 * g_y_xx_y_yz[i] + 2.0 * g_yyy_xx_y_yz[i] * a_exp;

        g_y_0_0_0_yy_xx_y_zz[i] = -2.0 * g_y_xx_y_zz[i] + 2.0 * g_yyy_xx_y_zz[i] * a_exp;
    }
    // integrals block (984-990)

    #pragma omp simd aligned(g_y_0_0_0_yy_xx_z_xx, g_y_0_0_0_yy_xx_z_xy, g_y_0_0_0_yy_xx_z_xz, g_y_0_0_0_yy_xx_z_yy, g_y_0_0_0_yy_xx_z_yz, g_y_0_0_0_yy_xx_z_zz, g_y_xx_z_xx, g_y_xx_z_xy, g_y_xx_z_xz, g_y_xx_z_yy, g_y_xx_z_yz, g_y_xx_z_zz, g_yyy_xx_z_xx, g_yyy_xx_z_xy, g_yyy_xx_z_xz, g_yyy_xx_z_yy, g_yyy_xx_z_yz, g_yyy_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_xx_z_xx[i] = -2.0 * g_y_xx_z_xx[i] + 2.0 * g_yyy_xx_z_xx[i] * a_exp;

        g_y_0_0_0_yy_xx_z_xy[i] = -2.0 * g_y_xx_z_xy[i] + 2.0 * g_yyy_xx_z_xy[i] * a_exp;

        g_y_0_0_0_yy_xx_z_xz[i] = -2.0 * g_y_xx_z_xz[i] + 2.0 * g_yyy_xx_z_xz[i] * a_exp;

        g_y_0_0_0_yy_xx_z_yy[i] = -2.0 * g_y_xx_z_yy[i] + 2.0 * g_yyy_xx_z_yy[i] * a_exp;

        g_y_0_0_0_yy_xx_z_yz[i] = -2.0 * g_y_xx_z_yz[i] + 2.0 * g_yyy_xx_z_yz[i] * a_exp;

        g_y_0_0_0_yy_xx_z_zz[i] = -2.0 * g_y_xx_z_zz[i] + 2.0 * g_yyy_xx_z_zz[i] * a_exp;
    }
    // integrals block (990-996)

    #pragma omp simd aligned(g_y_0_0_0_yy_xy_x_xx, g_y_0_0_0_yy_xy_x_xy, g_y_0_0_0_yy_xy_x_xz, g_y_0_0_0_yy_xy_x_yy, g_y_0_0_0_yy_xy_x_yz, g_y_0_0_0_yy_xy_x_zz, g_y_xy_x_xx, g_y_xy_x_xy, g_y_xy_x_xz, g_y_xy_x_yy, g_y_xy_x_yz, g_y_xy_x_zz, g_yyy_xy_x_xx, g_yyy_xy_x_xy, g_yyy_xy_x_xz, g_yyy_xy_x_yy, g_yyy_xy_x_yz, g_yyy_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_xy_x_xx[i] = -2.0 * g_y_xy_x_xx[i] + 2.0 * g_yyy_xy_x_xx[i] * a_exp;

        g_y_0_0_0_yy_xy_x_xy[i] = -2.0 * g_y_xy_x_xy[i] + 2.0 * g_yyy_xy_x_xy[i] * a_exp;

        g_y_0_0_0_yy_xy_x_xz[i] = -2.0 * g_y_xy_x_xz[i] + 2.0 * g_yyy_xy_x_xz[i] * a_exp;

        g_y_0_0_0_yy_xy_x_yy[i] = -2.0 * g_y_xy_x_yy[i] + 2.0 * g_yyy_xy_x_yy[i] * a_exp;

        g_y_0_0_0_yy_xy_x_yz[i] = -2.0 * g_y_xy_x_yz[i] + 2.0 * g_yyy_xy_x_yz[i] * a_exp;

        g_y_0_0_0_yy_xy_x_zz[i] = -2.0 * g_y_xy_x_zz[i] + 2.0 * g_yyy_xy_x_zz[i] * a_exp;
    }
    // integrals block (996-1002)

    #pragma omp simd aligned(g_y_0_0_0_yy_xy_y_xx, g_y_0_0_0_yy_xy_y_xy, g_y_0_0_0_yy_xy_y_xz, g_y_0_0_0_yy_xy_y_yy, g_y_0_0_0_yy_xy_y_yz, g_y_0_0_0_yy_xy_y_zz, g_y_xy_y_xx, g_y_xy_y_xy, g_y_xy_y_xz, g_y_xy_y_yy, g_y_xy_y_yz, g_y_xy_y_zz, g_yyy_xy_y_xx, g_yyy_xy_y_xy, g_yyy_xy_y_xz, g_yyy_xy_y_yy, g_yyy_xy_y_yz, g_yyy_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_xy_y_xx[i] = -2.0 * g_y_xy_y_xx[i] + 2.0 * g_yyy_xy_y_xx[i] * a_exp;

        g_y_0_0_0_yy_xy_y_xy[i] = -2.0 * g_y_xy_y_xy[i] + 2.0 * g_yyy_xy_y_xy[i] * a_exp;

        g_y_0_0_0_yy_xy_y_xz[i] = -2.0 * g_y_xy_y_xz[i] + 2.0 * g_yyy_xy_y_xz[i] * a_exp;

        g_y_0_0_0_yy_xy_y_yy[i] = -2.0 * g_y_xy_y_yy[i] + 2.0 * g_yyy_xy_y_yy[i] * a_exp;

        g_y_0_0_0_yy_xy_y_yz[i] = -2.0 * g_y_xy_y_yz[i] + 2.0 * g_yyy_xy_y_yz[i] * a_exp;

        g_y_0_0_0_yy_xy_y_zz[i] = -2.0 * g_y_xy_y_zz[i] + 2.0 * g_yyy_xy_y_zz[i] * a_exp;
    }
    // integrals block (1002-1008)

    #pragma omp simd aligned(g_y_0_0_0_yy_xy_z_xx, g_y_0_0_0_yy_xy_z_xy, g_y_0_0_0_yy_xy_z_xz, g_y_0_0_0_yy_xy_z_yy, g_y_0_0_0_yy_xy_z_yz, g_y_0_0_0_yy_xy_z_zz, g_y_xy_z_xx, g_y_xy_z_xy, g_y_xy_z_xz, g_y_xy_z_yy, g_y_xy_z_yz, g_y_xy_z_zz, g_yyy_xy_z_xx, g_yyy_xy_z_xy, g_yyy_xy_z_xz, g_yyy_xy_z_yy, g_yyy_xy_z_yz, g_yyy_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_xy_z_xx[i] = -2.0 * g_y_xy_z_xx[i] + 2.0 * g_yyy_xy_z_xx[i] * a_exp;

        g_y_0_0_0_yy_xy_z_xy[i] = -2.0 * g_y_xy_z_xy[i] + 2.0 * g_yyy_xy_z_xy[i] * a_exp;

        g_y_0_0_0_yy_xy_z_xz[i] = -2.0 * g_y_xy_z_xz[i] + 2.0 * g_yyy_xy_z_xz[i] * a_exp;

        g_y_0_0_0_yy_xy_z_yy[i] = -2.0 * g_y_xy_z_yy[i] + 2.0 * g_yyy_xy_z_yy[i] * a_exp;

        g_y_0_0_0_yy_xy_z_yz[i] = -2.0 * g_y_xy_z_yz[i] + 2.0 * g_yyy_xy_z_yz[i] * a_exp;

        g_y_0_0_0_yy_xy_z_zz[i] = -2.0 * g_y_xy_z_zz[i] + 2.0 * g_yyy_xy_z_zz[i] * a_exp;
    }
    // integrals block (1008-1014)

    #pragma omp simd aligned(g_y_0_0_0_yy_xz_x_xx, g_y_0_0_0_yy_xz_x_xy, g_y_0_0_0_yy_xz_x_xz, g_y_0_0_0_yy_xz_x_yy, g_y_0_0_0_yy_xz_x_yz, g_y_0_0_0_yy_xz_x_zz, g_y_xz_x_xx, g_y_xz_x_xy, g_y_xz_x_xz, g_y_xz_x_yy, g_y_xz_x_yz, g_y_xz_x_zz, g_yyy_xz_x_xx, g_yyy_xz_x_xy, g_yyy_xz_x_xz, g_yyy_xz_x_yy, g_yyy_xz_x_yz, g_yyy_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_xz_x_xx[i] = -2.0 * g_y_xz_x_xx[i] + 2.0 * g_yyy_xz_x_xx[i] * a_exp;

        g_y_0_0_0_yy_xz_x_xy[i] = -2.0 * g_y_xz_x_xy[i] + 2.0 * g_yyy_xz_x_xy[i] * a_exp;

        g_y_0_0_0_yy_xz_x_xz[i] = -2.0 * g_y_xz_x_xz[i] + 2.0 * g_yyy_xz_x_xz[i] * a_exp;

        g_y_0_0_0_yy_xz_x_yy[i] = -2.0 * g_y_xz_x_yy[i] + 2.0 * g_yyy_xz_x_yy[i] * a_exp;

        g_y_0_0_0_yy_xz_x_yz[i] = -2.0 * g_y_xz_x_yz[i] + 2.0 * g_yyy_xz_x_yz[i] * a_exp;

        g_y_0_0_0_yy_xz_x_zz[i] = -2.0 * g_y_xz_x_zz[i] + 2.0 * g_yyy_xz_x_zz[i] * a_exp;
    }
    // integrals block (1014-1020)

    #pragma omp simd aligned(g_y_0_0_0_yy_xz_y_xx, g_y_0_0_0_yy_xz_y_xy, g_y_0_0_0_yy_xz_y_xz, g_y_0_0_0_yy_xz_y_yy, g_y_0_0_0_yy_xz_y_yz, g_y_0_0_0_yy_xz_y_zz, g_y_xz_y_xx, g_y_xz_y_xy, g_y_xz_y_xz, g_y_xz_y_yy, g_y_xz_y_yz, g_y_xz_y_zz, g_yyy_xz_y_xx, g_yyy_xz_y_xy, g_yyy_xz_y_xz, g_yyy_xz_y_yy, g_yyy_xz_y_yz, g_yyy_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_xz_y_xx[i] = -2.0 * g_y_xz_y_xx[i] + 2.0 * g_yyy_xz_y_xx[i] * a_exp;

        g_y_0_0_0_yy_xz_y_xy[i] = -2.0 * g_y_xz_y_xy[i] + 2.0 * g_yyy_xz_y_xy[i] * a_exp;

        g_y_0_0_0_yy_xz_y_xz[i] = -2.0 * g_y_xz_y_xz[i] + 2.0 * g_yyy_xz_y_xz[i] * a_exp;

        g_y_0_0_0_yy_xz_y_yy[i] = -2.0 * g_y_xz_y_yy[i] + 2.0 * g_yyy_xz_y_yy[i] * a_exp;

        g_y_0_0_0_yy_xz_y_yz[i] = -2.0 * g_y_xz_y_yz[i] + 2.0 * g_yyy_xz_y_yz[i] * a_exp;

        g_y_0_0_0_yy_xz_y_zz[i] = -2.0 * g_y_xz_y_zz[i] + 2.0 * g_yyy_xz_y_zz[i] * a_exp;
    }
    // integrals block (1020-1026)

    #pragma omp simd aligned(g_y_0_0_0_yy_xz_z_xx, g_y_0_0_0_yy_xz_z_xy, g_y_0_0_0_yy_xz_z_xz, g_y_0_0_0_yy_xz_z_yy, g_y_0_0_0_yy_xz_z_yz, g_y_0_0_0_yy_xz_z_zz, g_y_xz_z_xx, g_y_xz_z_xy, g_y_xz_z_xz, g_y_xz_z_yy, g_y_xz_z_yz, g_y_xz_z_zz, g_yyy_xz_z_xx, g_yyy_xz_z_xy, g_yyy_xz_z_xz, g_yyy_xz_z_yy, g_yyy_xz_z_yz, g_yyy_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_xz_z_xx[i] = -2.0 * g_y_xz_z_xx[i] + 2.0 * g_yyy_xz_z_xx[i] * a_exp;

        g_y_0_0_0_yy_xz_z_xy[i] = -2.0 * g_y_xz_z_xy[i] + 2.0 * g_yyy_xz_z_xy[i] * a_exp;

        g_y_0_0_0_yy_xz_z_xz[i] = -2.0 * g_y_xz_z_xz[i] + 2.0 * g_yyy_xz_z_xz[i] * a_exp;

        g_y_0_0_0_yy_xz_z_yy[i] = -2.0 * g_y_xz_z_yy[i] + 2.0 * g_yyy_xz_z_yy[i] * a_exp;

        g_y_0_0_0_yy_xz_z_yz[i] = -2.0 * g_y_xz_z_yz[i] + 2.0 * g_yyy_xz_z_yz[i] * a_exp;

        g_y_0_0_0_yy_xz_z_zz[i] = -2.0 * g_y_xz_z_zz[i] + 2.0 * g_yyy_xz_z_zz[i] * a_exp;
    }
    // integrals block (1026-1032)

    #pragma omp simd aligned(g_y_0_0_0_yy_yy_x_xx, g_y_0_0_0_yy_yy_x_xy, g_y_0_0_0_yy_yy_x_xz, g_y_0_0_0_yy_yy_x_yy, g_y_0_0_0_yy_yy_x_yz, g_y_0_0_0_yy_yy_x_zz, g_y_yy_x_xx, g_y_yy_x_xy, g_y_yy_x_xz, g_y_yy_x_yy, g_y_yy_x_yz, g_y_yy_x_zz, g_yyy_yy_x_xx, g_yyy_yy_x_xy, g_yyy_yy_x_xz, g_yyy_yy_x_yy, g_yyy_yy_x_yz, g_yyy_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_yy_x_xx[i] = -2.0 * g_y_yy_x_xx[i] + 2.0 * g_yyy_yy_x_xx[i] * a_exp;

        g_y_0_0_0_yy_yy_x_xy[i] = -2.0 * g_y_yy_x_xy[i] + 2.0 * g_yyy_yy_x_xy[i] * a_exp;

        g_y_0_0_0_yy_yy_x_xz[i] = -2.0 * g_y_yy_x_xz[i] + 2.0 * g_yyy_yy_x_xz[i] * a_exp;

        g_y_0_0_0_yy_yy_x_yy[i] = -2.0 * g_y_yy_x_yy[i] + 2.0 * g_yyy_yy_x_yy[i] * a_exp;

        g_y_0_0_0_yy_yy_x_yz[i] = -2.0 * g_y_yy_x_yz[i] + 2.0 * g_yyy_yy_x_yz[i] * a_exp;

        g_y_0_0_0_yy_yy_x_zz[i] = -2.0 * g_y_yy_x_zz[i] + 2.0 * g_yyy_yy_x_zz[i] * a_exp;
    }
    // integrals block (1032-1038)

    #pragma omp simd aligned(g_y_0_0_0_yy_yy_y_xx, g_y_0_0_0_yy_yy_y_xy, g_y_0_0_0_yy_yy_y_xz, g_y_0_0_0_yy_yy_y_yy, g_y_0_0_0_yy_yy_y_yz, g_y_0_0_0_yy_yy_y_zz, g_y_yy_y_xx, g_y_yy_y_xy, g_y_yy_y_xz, g_y_yy_y_yy, g_y_yy_y_yz, g_y_yy_y_zz, g_yyy_yy_y_xx, g_yyy_yy_y_xy, g_yyy_yy_y_xz, g_yyy_yy_y_yy, g_yyy_yy_y_yz, g_yyy_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_yy_y_xx[i] = -2.0 * g_y_yy_y_xx[i] + 2.0 * g_yyy_yy_y_xx[i] * a_exp;

        g_y_0_0_0_yy_yy_y_xy[i] = -2.0 * g_y_yy_y_xy[i] + 2.0 * g_yyy_yy_y_xy[i] * a_exp;

        g_y_0_0_0_yy_yy_y_xz[i] = -2.0 * g_y_yy_y_xz[i] + 2.0 * g_yyy_yy_y_xz[i] * a_exp;

        g_y_0_0_0_yy_yy_y_yy[i] = -2.0 * g_y_yy_y_yy[i] + 2.0 * g_yyy_yy_y_yy[i] * a_exp;

        g_y_0_0_0_yy_yy_y_yz[i] = -2.0 * g_y_yy_y_yz[i] + 2.0 * g_yyy_yy_y_yz[i] * a_exp;

        g_y_0_0_0_yy_yy_y_zz[i] = -2.0 * g_y_yy_y_zz[i] + 2.0 * g_yyy_yy_y_zz[i] * a_exp;
    }
    // integrals block (1038-1044)

    #pragma omp simd aligned(g_y_0_0_0_yy_yy_z_xx, g_y_0_0_0_yy_yy_z_xy, g_y_0_0_0_yy_yy_z_xz, g_y_0_0_0_yy_yy_z_yy, g_y_0_0_0_yy_yy_z_yz, g_y_0_0_0_yy_yy_z_zz, g_y_yy_z_xx, g_y_yy_z_xy, g_y_yy_z_xz, g_y_yy_z_yy, g_y_yy_z_yz, g_y_yy_z_zz, g_yyy_yy_z_xx, g_yyy_yy_z_xy, g_yyy_yy_z_xz, g_yyy_yy_z_yy, g_yyy_yy_z_yz, g_yyy_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_yy_z_xx[i] = -2.0 * g_y_yy_z_xx[i] + 2.0 * g_yyy_yy_z_xx[i] * a_exp;

        g_y_0_0_0_yy_yy_z_xy[i] = -2.0 * g_y_yy_z_xy[i] + 2.0 * g_yyy_yy_z_xy[i] * a_exp;

        g_y_0_0_0_yy_yy_z_xz[i] = -2.0 * g_y_yy_z_xz[i] + 2.0 * g_yyy_yy_z_xz[i] * a_exp;

        g_y_0_0_0_yy_yy_z_yy[i] = -2.0 * g_y_yy_z_yy[i] + 2.0 * g_yyy_yy_z_yy[i] * a_exp;

        g_y_0_0_0_yy_yy_z_yz[i] = -2.0 * g_y_yy_z_yz[i] + 2.0 * g_yyy_yy_z_yz[i] * a_exp;

        g_y_0_0_0_yy_yy_z_zz[i] = -2.0 * g_y_yy_z_zz[i] + 2.0 * g_yyy_yy_z_zz[i] * a_exp;
    }
    // integrals block (1044-1050)

    #pragma omp simd aligned(g_y_0_0_0_yy_yz_x_xx, g_y_0_0_0_yy_yz_x_xy, g_y_0_0_0_yy_yz_x_xz, g_y_0_0_0_yy_yz_x_yy, g_y_0_0_0_yy_yz_x_yz, g_y_0_0_0_yy_yz_x_zz, g_y_yz_x_xx, g_y_yz_x_xy, g_y_yz_x_xz, g_y_yz_x_yy, g_y_yz_x_yz, g_y_yz_x_zz, g_yyy_yz_x_xx, g_yyy_yz_x_xy, g_yyy_yz_x_xz, g_yyy_yz_x_yy, g_yyy_yz_x_yz, g_yyy_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_yz_x_xx[i] = -2.0 * g_y_yz_x_xx[i] + 2.0 * g_yyy_yz_x_xx[i] * a_exp;

        g_y_0_0_0_yy_yz_x_xy[i] = -2.0 * g_y_yz_x_xy[i] + 2.0 * g_yyy_yz_x_xy[i] * a_exp;

        g_y_0_0_0_yy_yz_x_xz[i] = -2.0 * g_y_yz_x_xz[i] + 2.0 * g_yyy_yz_x_xz[i] * a_exp;

        g_y_0_0_0_yy_yz_x_yy[i] = -2.0 * g_y_yz_x_yy[i] + 2.0 * g_yyy_yz_x_yy[i] * a_exp;

        g_y_0_0_0_yy_yz_x_yz[i] = -2.0 * g_y_yz_x_yz[i] + 2.0 * g_yyy_yz_x_yz[i] * a_exp;

        g_y_0_0_0_yy_yz_x_zz[i] = -2.0 * g_y_yz_x_zz[i] + 2.0 * g_yyy_yz_x_zz[i] * a_exp;
    }
    // integrals block (1050-1056)

    #pragma omp simd aligned(g_y_0_0_0_yy_yz_y_xx, g_y_0_0_0_yy_yz_y_xy, g_y_0_0_0_yy_yz_y_xz, g_y_0_0_0_yy_yz_y_yy, g_y_0_0_0_yy_yz_y_yz, g_y_0_0_0_yy_yz_y_zz, g_y_yz_y_xx, g_y_yz_y_xy, g_y_yz_y_xz, g_y_yz_y_yy, g_y_yz_y_yz, g_y_yz_y_zz, g_yyy_yz_y_xx, g_yyy_yz_y_xy, g_yyy_yz_y_xz, g_yyy_yz_y_yy, g_yyy_yz_y_yz, g_yyy_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_yz_y_xx[i] = -2.0 * g_y_yz_y_xx[i] + 2.0 * g_yyy_yz_y_xx[i] * a_exp;

        g_y_0_0_0_yy_yz_y_xy[i] = -2.0 * g_y_yz_y_xy[i] + 2.0 * g_yyy_yz_y_xy[i] * a_exp;

        g_y_0_0_0_yy_yz_y_xz[i] = -2.0 * g_y_yz_y_xz[i] + 2.0 * g_yyy_yz_y_xz[i] * a_exp;

        g_y_0_0_0_yy_yz_y_yy[i] = -2.0 * g_y_yz_y_yy[i] + 2.0 * g_yyy_yz_y_yy[i] * a_exp;

        g_y_0_0_0_yy_yz_y_yz[i] = -2.0 * g_y_yz_y_yz[i] + 2.0 * g_yyy_yz_y_yz[i] * a_exp;

        g_y_0_0_0_yy_yz_y_zz[i] = -2.0 * g_y_yz_y_zz[i] + 2.0 * g_yyy_yz_y_zz[i] * a_exp;
    }
    // integrals block (1056-1062)

    #pragma omp simd aligned(g_y_0_0_0_yy_yz_z_xx, g_y_0_0_0_yy_yz_z_xy, g_y_0_0_0_yy_yz_z_xz, g_y_0_0_0_yy_yz_z_yy, g_y_0_0_0_yy_yz_z_yz, g_y_0_0_0_yy_yz_z_zz, g_y_yz_z_xx, g_y_yz_z_xy, g_y_yz_z_xz, g_y_yz_z_yy, g_y_yz_z_yz, g_y_yz_z_zz, g_yyy_yz_z_xx, g_yyy_yz_z_xy, g_yyy_yz_z_xz, g_yyy_yz_z_yy, g_yyy_yz_z_yz, g_yyy_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_yz_z_xx[i] = -2.0 * g_y_yz_z_xx[i] + 2.0 * g_yyy_yz_z_xx[i] * a_exp;

        g_y_0_0_0_yy_yz_z_xy[i] = -2.0 * g_y_yz_z_xy[i] + 2.0 * g_yyy_yz_z_xy[i] * a_exp;

        g_y_0_0_0_yy_yz_z_xz[i] = -2.0 * g_y_yz_z_xz[i] + 2.0 * g_yyy_yz_z_xz[i] * a_exp;

        g_y_0_0_0_yy_yz_z_yy[i] = -2.0 * g_y_yz_z_yy[i] + 2.0 * g_yyy_yz_z_yy[i] * a_exp;

        g_y_0_0_0_yy_yz_z_yz[i] = -2.0 * g_y_yz_z_yz[i] + 2.0 * g_yyy_yz_z_yz[i] * a_exp;

        g_y_0_0_0_yy_yz_z_zz[i] = -2.0 * g_y_yz_z_zz[i] + 2.0 * g_yyy_yz_z_zz[i] * a_exp;
    }
    // integrals block (1062-1068)

    #pragma omp simd aligned(g_y_0_0_0_yy_zz_x_xx, g_y_0_0_0_yy_zz_x_xy, g_y_0_0_0_yy_zz_x_xz, g_y_0_0_0_yy_zz_x_yy, g_y_0_0_0_yy_zz_x_yz, g_y_0_0_0_yy_zz_x_zz, g_y_zz_x_xx, g_y_zz_x_xy, g_y_zz_x_xz, g_y_zz_x_yy, g_y_zz_x_yz, g_y_zz_x_zz, g_yyy_zz_x_xx, g_yyy_zz_x_xy, g_yyy_zz_x_xz, g_yyy_zz_x_yy, g_yyy_zz_x_yz, g_yyy_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_zz_x_xx[i] = -2.0 * g_y_zz_x_xx[i] + 2.0 * g_yyy_zz_x_xx[i] * a_exp;

        g_y_0_0_0_yy_zz_x_xy[i] = -2.0 * g_y_zz_x_xy[i] + 2.0 * g_yyy_zz_x_xy[i] * a_exp;

        g_y_0_0_0_yy_zz_x_xz[i] = -2.0 * g_y_zz_x_xz[i] + 2.0 * g_yyy_zz_x_xz[i] * a_exp;

        g_y_0_0_0_yy_zz_x_yy[i] = -2.0 * g_y_zz_x_yy[i] + 2.0 * g_yyy_zz_x_yy[i] * a_exp;

        g_y_0_0_0_yy_zz_x_yz[i] = -2.0 * g_y_zz_x_yz[i] + 2.0 * g_yyy_zz_x_yz[i] * a_exp;

        g_y_0_0_0_yy_zz_x_zz[i] = -2.0 * g_y_zz_x_zz[i] + 2.0 * g_yyy_zz_x_zz[i] * a_exp;
    }
    // integrals block (1068-1074)

    #pragma omp simd aligned(g_y_0_0_0_yy_zz_y_xx, g_y_0_0_0_yy_zz_y_xy, g_y_0_0_0_yy_zz_y_xz, g_y_0_0_0_yy_zz_y_yy, g_y_0_0_0_yy_zz_y_yz, g_y_0_0_0_yy_zz_y_zz, g_y_zz_y_xx, g_y_zz_y_xy, g_y_zz_y_xz, g_y_zz_y_yy, g_y_zz_y_yz, g_y_zz_y_zz, g_yyy_zz_y_xx, g_yyy_zz_y_xy, g_yyy_zz_y_xz, g_yyy_zz_y_yy, g_yyy_zz_y_yz, g_yyy_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_zz_y_xx[i] = -2.0 * g_y_zz_y_xx[i] + 2.0 * g_yyy_zz_y_xx[i] * a_exp;

        g_y_0_0_0_yy_zz_y_xy[i] = -2.0 * g_y_zz_y_xy[i] + 2.0 * g_yyy_zz_y_xy[i] * a_exp;

        g_y_0_0_0_yy_zz_y_xz[i] = -2.0 * g_y_zz_y_xz[i] + 2.0 * g_yyy_zz_y_xz[i] * a_exp;

        g_y_0_0_0_yy_zz_y_yy[i] = -2.0 * g_y_zz_y_yy[i] + 2.0 * g_yyy_zz_y_yy[i] * a_exp;

        g_y_0_0_0_yy_zz_y_yz[i] = -2.0 * g_y_zz_y_yz[i] + 2.0 * g_yyy_zz_y_yz[i] * a_exp;

        g_y_0_0_0_yy_zz_y_zz[i] = -2.0 * g_y_zz_y_zz[i] + 2.0 * g_yyy_zz_y_zz[i] * a_exp;
    }
    // integrals block (1074-1080)

    #pragma omp simd aligned(g_y_0_0_0_yy_zz_z_xx, g_y_0_0_0_yy_zz_z_xy, g_y_0_0_0_yy_zz_z_xz, g_y_0_0_0_yy_zz_z_yy, g_y_0_0_0_yy_zz_z_yz, g_y_0_0_0_yy_zz_z_zz, g_y_zz_z_xx, g_y_zz_z_xy, g_y_zz_z_xz, g_y_zz_z_yy, g_y_zz_z_yz, g_y_zz_z_zz, g_yyy_zz_z_xx, g_yyy_zz_z_xy, g_yyy_zz_z_xz, g_yyy_zz_z_yy, g_yyy_zz_z_yz, g_yyy_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_zz_z_xx[i] = -2.0 * g_y_zz_z_xx[i] + 2.0 * g_yyy_zz_z_xx[i] * a_exp;

        g_y_0_0_0_yy_zz_z_xy[i] = -2.0 * g_y_zz_z_xy[i] + 2.0 * g_yyy_zz_z_xy[i] * a_exp;

        g_y_0_0_0_yy_zz_z_xz[i] = -2.0 * g_y_zz_z_xz[i] + 2.0 * g_yyy_zz_z_xz[i] * a_exp;

        g_y_0_0_0_yy_zz_z_yy[i] = -2.0 * g_y_zz_z_yy[i] + 2.0 * g_yyy_zz_z_yy[i] * a_exp;

        g_y_0_0_0_yy_zz_z_yz[i] = -2.0 * g_y_zz_z_yz[i] + 2.0 * g_yyy_zz_z_yz[i] * a_exp;

        g_y_0_0_0_yy_zz_z_zz[i] = -2.0 * g_y_zz_z_zz[i] + 2.0 * g_yyy_zz_z_zz[i] * a_exp;
    }
    // integrals block (1080-1086)

    #pragma omp simd aligned(g_y_0_0_0_yz_xx_x_xx, g_y_0_0_0_yz_xx_x_xy, g_y_0_0_0_yz_xx_x_xz, g_y_0_0_0_yz_xx_x_yy, g_y_0_0_0_yz_xx_x_yz, g_y_0_0_0_yz_xx_x_zz, g_yyz_xx_x_xx, g_yyz_xx_x_xy, g_yyz_xx_x_xz, g_yyz_xx_x_yy, g_yyz_xx_x_yz, g_yyz_xx_x_zz, g_z_xx_x_xx, g_z_xx_x_xy, g_z_xx_x_xz, g_z_xx_x_yy, g_z_xx_x_yz, g_z_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_xx_x_xx[i] = -g_z_xx_x_xx[i] + 2.0 * g_yyz_xx_x_xx[i] * a_exp;

        g_y_0_0_0_yz_xx_x_xy[i] = -g_z_xx_x_xy[i] + 2.0 * g_yyz_xx_x_xy[i] * a_exp;

        g_y_0_0_0_yz_xx_x_xz[i] = -g_z_xx_x_xz[i] + 2.0 * g_yyz_xx_x_xz[i] * a_exp;

        g_y_0_0_0_yz_xx_x_yy[i] = -g_z_xx_x_yy[i] + 2.0 * g_yyz_xx_x_yy[i] * a_exp;

        g_y_0_0_0_yz_xx_x_yz[i] = -g_z_xx_x_yz[i] + 2.0 * g_yyz_xx_x_yz[i] * a_exp;

        g_y_0_0_0_yz_xx_x_zz[i] = -g_z_xx_x_zz[i] + 2.0 * g_yyz_xx_x_zz[i] * a_exp;
    }
    // integrals block (1086-1092)

    #pragma omp simd aligned(g_y_0_0_0_yz_xx_y_xx, g_y_0_0_0_yz_xx_y_xy, g_y_0_0_0_yz_xx_y_xz, g_y_0_0_0_yz_xx_y_yy, g_y_0_0_0_yz_xx_y_yz, g_y_0_0_0_yz_xx_y_zz, g_yyz_xx_y_xx, g_yyz_xx_y_xy, g_yyz_xx_y_xz, g_yyz_xx_y_yy, g_yyz_xx_y_yz, g_yyz_xx_y_zz, g_z_xx_y_xx, g_z_xx_y_xy, g_z_xx_y_xz, g_z_xx_y_yy, g_z_xx_y_yz, g_z_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_xx_y_xx[i] = -g_z_xx_y_xx[i] + 2.0 * g_yyz_xx_y_xx[i] * a_exp;

        g_y_0_0_0_yz_xx_y_xy[i] = -g_z_xx_y_xy[i] + 2.0 * g_yyz_xx_y_xy[i] * a_exp;

        g_y_0_0_0_yz_xx_y_xz[i] = -g_z_xx_y_xz[i] + 2.0 * g_yyz_xx_y_xz[i] * a_exp;

        g_y_0_0_0_yz_xx_y_yy[i] = -g_z_xx_y_yy[i] + 2.0 * g_yyz_xx_y_yy[i] * a_exp;

        g_y_0_0_0_yz_xx_y_yz[i] = -g_z_xx_y_yz[i] + 2.0 * g_yyz_xx_y_yz[i] * a_exp;

        g_y_0_0_0_yz_xx_y_zz[i] = -g_z_xx_y_zz[i] + 2.0 * g_yyz_xx_y_zz[i] * a_exp;
    }
    // integrals block (1092-1098)

    #pragma omp simd aligned(g_y_0_0_0_yz_xx_z_xx, g_y_0_0_0_yz_xx_z_xy, g_y_0_0_0_yz_xx_z_xz, g_y_0_0_0_yz_xx_z_yy, g_y_0_0_0_yz_xx_z_yz, g_y_0_0_0_yz_xx_z_zz, g_yyz_xx_z_xx, g_yyz_xx_z_xy, g_yyz_xx_z_xz, g_yyz_xx_z_yy, g_yyz_xx_z_yz, g_yyz_xx_z_zz, g_z_xx_z_xx, g_z_xx_z_xy, g_z_xx_z_xz, g_z_xx_z_yy, g_z_xx_z_yz, g_z_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_xx_z_xx[i] = -g_z_xx_z_xx[i] + 2.0 * g_yyz_xx_z_xx[i] * a_exp;

        g_y_0_0_0_yz_xx_z_xy[i] = -g_z_xx_z_xy[i] + 2.0 * g_yyz_xx_z_xy[i] * a_exp;

        g_y_0_0_0_yz_xx_z_xz[i] = -g_z_xx_z_xz[i] + 2.0 * g_yyz_xx_z_xz[i] * a_exp;

        g_y_0_0_0_yz_xx_z_yy[i] = -g_z_xx_z_yy[i] + 2.0 * g_yyz_xx_z_yy[i] * a_exp;

        g_y_0_0_0_yz_xx_z_yz[i] = -g_z_xx_z_yz[i] + 2.0 * g_yyz_xx_z_yz[i] * a_exp;

        g_y_0_0_0_yz_xx_z_zz[i] = -g_z_xx_z_zz[i] + 2.0 * g_yyz_xx_z_zz[i] * a_exp;
    }
    // integrals block (1098-1104)

    #pragma omp simd aligned(g_y_0_0_0_yz_xy_x_xx, g_y_0_0_0_yz_xy_x_xy, g_y_0_0_0_yz_xy_x_xz, g_y_0_0_0_yz_xy_x_yy, g_y_0_0_0_yz_xy_x_yz, g_y_0_0_0_yz_xy_x_zz, g_yyz_xy_x_xx, g_yyz_xy_x_xy, g_yyz_xy_x_xz, g_yyz_xy_x_yy, g_yyz_xy_x_yz, g_yyz_xy_x_zz, g_z_xy_x_xx, g_z_xy_x_xy, g_z_xy_x_xz, g_z_xy_x_yy, g_z_xy_x_yz, g_z_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_xy_x_xx[i] = -g_z_xy_x_xx[i] + 2.0 * g_yyz_xy_x_xx[i] * a_exp;

        g_y_0_0_0_yz_xy_x_xy[i] = -g_z_xy_x_xy[i] + 2.0 * g_yyz_xy_x_xy[i] * a_exp;

        g_y_0_0_0_yz_xy_x_xz[i] = -g_z_xy_x_xz[i] + 2.0 * g_yyz_xy_x_xz[i] * a_exp;

        g_y_0_0_0_yz_xy_x_yy[i] = -g_z_xy_x_yy[i] + 2.0 * g_yyz_xy_x_yy[i] * a_exp;

        g_y_0_0_0_yz_xy_x_yz[i] = -g_z_xy_x_yz[i] + 2.0 * g_yyz_xy_x_yz[i] * a_exp;

        g_y_0_0_0_yz_xy_x_zz[i] = -g_z_xy_x_zz[i] + 2.0 * g_yyz_xy_x_zz[i] * a_exp;
    }
    // integrals block (1104-1110)

    #pragma omp simd aligned(g_y_0_0_0_yz_xy_y_xx, g_y_0_0_0_yz_xy_y_xy, g_y_0_0_0_yz_xy_y_xz, g_y_0_0_0_yz_xy_y_yy, g_y_0_0_0_yz_xy_y_yz, g_y_0_0_0_yz_xy_y_zz, g_yyz_xy_y_xx, g_yyz_xy_y_xy, g_yyz_xy_y_xz, g_yyz_xy_y_yy, g_yyz_xy_y_yz, g_yyz_xy_y_zz, g_z_xy_y_xx, g_z_xy_y_xy, g_z_xy_y_xz, g_z_xy_y_yy, g_z_xy_y_yz, g_z_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_xy_y_xx[i] = -g_z_xy_y_xx[i] + 2.0 * g_yyz_xy_y_xx[i] * a_exp;

        g_y_0_0_0_yz_xy_y_xy[i] = -g_z_xy_y_xy[i] + 2.0 * g_yyz_xy_y_xy[i] * a_exp;

        g_y_0_0_0_yz_xy_y_xz[i] = -g_z_xy_y_xz[i] + 2.0 * g_yyz_xy_y_xz[i] * a_exp;

        g_y_0_0_0_yz_xy_y_yy[i] = -g_z_xy_y_yy[i] + 2.0 * g_yyz_xy_y_yy[i] * a_exp;

        g_y_0_0_0_yz_xy_y_yz[i] = -g_z_xy_y_yz[i] + 2.0 * g_yyz_xy_y_yz[i] * a_exp;

        g_y_0_0_0_yz_xy_y_zz[i] = -g_z_xy_y_zz[i] + 2.0 * g_yyz_xy_y_zz[i] * a_exp;
    }
    // integrals block (1110-1116)

    #pragma omp simd aligned(g_y_0_0_0_yz_xy_z_xx, g_y_0_0_0_yz_xy_z_xy, g_y_0_0_0_yz_xy_z_xz, g_y_0_0_0_yz_xy_z_yy, g_y_0_0_0_yz_xy_z_yz, g_y_0_0_0_yz_xy_z_zz, g_yyz_xy_z_xx, g_yyz_xy_z_xy, g_yyz_xy_z_xz, g_yyz_xy_z_yy, g_yyz_xy_z_yz, g_yyz_xy_z_zz, g_z_xy_z_xx, g_z_xy_z_xy, g_z_xy_z_xz, g_z_xy_z_yy, g_z_xy_z_yz, g_z_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_xy_z_xx[i] = -g_z_xy_z_xx[i] + 2.0 * g_yyz_xy_z_xx[i] * a_exp;

        g_y_0_0_0_yz_xy_z_xy[i] = -g_z_xy_z_xy[i] + 2.0 * g_yyz_xy_z_xy[i] * a_exp;

        g_y_0_0_0_yz_xy_z_xz[i] = -g_z_xy_z_xz[i] + 2.0 * g_yyz_xy_z_xz[i] * a_exp;

        g_y_0_0_0_yz_xy_z_yy[i] = -g_z_xy_z_yy[i] + 2.0 * g_yyz_xy_z_yy[i] * a_exp;

        g_y_0_0_0_yz_xy_z_yz[i] = -g_z_xy_z_yz[i] + 2.0 * g_yyz_xy_z_yz[i] * a_exp;

        g_y_0_0_0_yz_xy_z_zz[i] = -g_z_xy_z_zz[i] + 2.0 * g_yyz_xy_z_zz[i] * a_exp;
    }
    // integrals block (1116-1122)

    #pragma omp simd aligned(g_y_0_0_0_yz_xz_x_xx, g_y_0_0_0_yz_xz_x_xy, g_y_0_0_0_yz_xz_x_xz, g_y_0_0_0_yz_xz_x_yy, g_y_0_0_0_yz_xz_x_yz, g_y_0_0_0_yz_xz_x_zz, g_yyz_xz_x_xx, g_yyz_xz_x_xy, g_yyz_xz_x_xz, g_yyz_xz_x_yy, g_yyz_xz_x_yz, g_yyz_xz_x_zz, g_z_xz_x_xx, g_z_xz_x_xy, g_z_xz_x_xz, g_z_xz_x_yy, g_z_xz_x_yz, g_z_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_xz_x_xx[i] = -g_z_xz_x_xx[i] + 2.0 * g_yyz_xz_x_xx[i] * a_exp;

        g_y_0_0_0_yz_xz_x_xy[i] = -g_z_xz_x_xy[i] + 2.0 * g_yyz_xz_x_xy[i] * a_exp;

        g_y_0_0_0_yz_xz_x_xz[i] = -g_z_xz_x_xz[i] + 2.0 * g_yyz_xz_x_xz[i] * a_exp;

        g_y_0_0_0_yz_xz_x_yy[i] = -g_z_xz_x_yy[i] + 2.0 * g_yyz_xz_x_yy[i] * a_exp;

        g_y_0_0_0_yz_xz_x_yz[i] = -g_z_xz_x_yz[i] + 2.0 * g_yyz_xz_x_yz[i] * a_exp;

        g_y_0_0_0_yz_xz_x_zz[i] = -g_z_xz_x_zz[i] + 2.0 * g_yyz_xz_x_zz[i] * a_exp;
    }
    // integrals block (1122-1128)

    #pragma omp simd aligned(g_y_0_0_0_yz_xz_y_xx, g_y_0_0_0_yz_xz_y_xy, g_y_0_0_0_yz_xz_y_xz, g_y_0_0_0_yz_xz_y_yy, g_y_0_0_0_yz_xz_y_yz, g_y_0_0_0_yz_xz_y_zz, g_yyz_xz_y_xx, g_yyz_xz_y_xy, g_yyz_xz_y_xz, g_yyz_xz_y_yy, g_yyz_xz_y_yz, g_yyz_xz_y_zz, g_z_xz_y_xx, g_z_xz_y_xy, g_z_xz_y_xz, g_z_xz_y_yy, g_z_xz_y_yz, g_z_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_xz_y_xx[i] = -g_z_xz_y_xx[i] + 2.0 * g_yyz_xz_y_xx[i] * a_exp;

        g_y_0_0_0_yz_xz_y_xy[i] = -g_z_xz_y_xy[i] + 2.0 * g_yyz_xz_y_xy[i] * a_exp;

        g_y_0_0_0_yz_xz_y_xz[i] = -g_z_xz_y_xz[i] + 2.0 * g_yyz_xz_y_xz[i] * a_exp;

        g_y_0_0_0_yz_xz_y_yy[i] = -g_z_xz_y_yy[i] + 2.0 * g_yyz_xz_y_yy[i] * a_exp;

        g_y_0_0_0_yz_xz_y_yz[i] = -g_z_xz_y_yz[i] + 2.0 * g_yyz_xz_y_yz[i] * a_exp;

        g_y_0_0_0_yz_xz_y_zz[i] = -g_z_xz_y_zz[i] + 2.0 * g_yyz_xz_y_zz[i] * a_exp;
    }
    // integrals block (1128-1134)

    #pragma omp simd aligned(g_y_0_0_0_yz_xz_z_xx, g_y_0_0_0_yz_xz_z_xy, g_y_0_0_0_yz_xz_z_xz, g_y_0_0_0_yz_xz_z_yy, g_y_0_0_0_yz_xz_z_yz, g_y_0_0_0_yz_xz_z_zz, g_yyz_xz_z_xx, g_yyz_xz_z_xy, g_yyz_xz_z_xz, g_yyz_xz_z_yy, g_yyz_xz_z_yz, g_yyz_xz_z_zz, g_z_xz_z_xx, g_z_xz_z_xy, g_z_xz_z_xz, g_z_xz_z_yy, g_z_xz_z_yz, g_z_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_xz_z_xx[i] = -g_z_xz_z_xx[i] + 2.0 * g_yyz_xz_z_xx[i] * a_exp;

        g_y_0_0_0_yz_xz_z_xy[i] = -g_z_xz_z_xy[i] + 2.0 * g_yyz_xz_z_xy[i] * a_exp;

        g_y_0_0_0_yz_xz_z_xz[i] = -g_z_xz_z_xz[i] + 2.0 * g_yyz_xz_z_xz[i] * a_exp;

        g_y_0_0_0_yz_xz_z_yy[i] = -g_z_xz_z_yy[i] + 2.0 * g_yyz_xz_z_yy[i] * a_exp;

        g_y_0_0_0_yz_xz_z_yz[i] = -g_z_xz_z_yz[i] + 2.0 * g_yyz_xz_z_yz[i] * a_exp;

        g_y_0_0_0_yz_xz_z_zz[i] = -g_z_xz_z_zz[i] + 2.0 * g_yyz_xz_z_zz[i] * a_exp;
    }
    // integrals block (1134-1140)

    #pragma omp simd aligned(g_y_0_0_0_yz_yy_x_xx, g_y_0_0_0_yz_yy_x_xy, g_y_0_0_0_yz_yy_x_xz, g_y_0_0_0_yz_yy_x_yy, g_y_0_0_0_yz_yy_x_yz, g_y_0_0_0_yz_yy_x_zz, g_yyz_yy_x_xx, g_yyz_yy_x_xy, g_yyz_yy_x_xz, g_yyz_yy_x_yy, g_yyz_yy_x_yz, g_yyz_yy_x_zz, g_z_yy_x_xx, g_z_yy_x_xy, g_z_yy_x_xz, g_z_yy_x_yy, g_z_yy_x_yz, g_z_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_yy_x_xx[i] = -g_z_yy_x_xx[i] + 2.0 * g_yyz_yy_x_xx[i] * a_exp;

        g_y_0_0_0_yz_yy_x_xy[i] = -g_z_yy_x_xy[i] + 2.0 * g_yyz_yy_x_xy[i] * a_exp;

        g_y_0_0_0_yz_yy_x_xz[i] = -g_z_yy_x_xz[i] + 2.0 * g_yyz_yy_x_xz[i] * a_exp;

        g_y_0_0_0_yz_yy_x_yy[i] = -g_z_yy_x_yy[i] + 2.0 * g_yyz_yy_x_yy[i] * a_exp;

        g_y_0_0_0_yz_yy_x_yz[i] = -g_z_yy_x_yz[i] + 2.0 * g_yyz_yy_x_yz[i] * a_exp;

        g_y_0_0_0_yz_yy_x_zz[i] = -g_z_yy_x_zz[i] + 2.0 * g_yyz_yy_x_zz[i] * a_exp;
    }
    // integrals block (1140-1146)

    #pragma omp simd aligned(g_y_0_0_0_yz_yy_y_xx, g_y_0_0_0_yz_yy_y_xy, g_y_0_0_0_yz_yy_y_xz, g_y_0_0_0_yz_yy_y_yy, g_y_0_0_0_yz_yy_y_yz, g_y_0_0_0_yz_yy_y_zz, g_yyz_yy_y_xx, g_yyz_yy_y_xy, g_yyz_yy_y_xz, g_yyz_yy_y_yy, g_yyz_yy_y_yz, g_yyz_yy_y_zz, g_z_yy_y_xx, g_z_yy_y_xy, g_z_yy_y_xz, g_z_yy_y_yy, g_z_yy_y_yz, g_z_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_yy_y_xx[i] = -g_z_yy_y_xx[i] + 2.0 * g_yyz_yy_y_xx[i] * a_exp;

        g_y_0_0_0_yz_yy_y_xy[i] = -g_z_yy_y_xy[i] + 2.0 * g_yyz_yy_y_xy[i] * a_exp;

        g_y_0_0_0_yz_yy_y_xz[i] = -g_z_yy_y_xz[i] + 2.0 * g_yyz_yy_y_xz[i] * a_exp;

        g_y_0_0_0_yz_yy_y_yy[i] = -g_z_yy_y_yy[i] + 2.0 * g_yyz_yy_y_yy[i] * a_exp;

        g_y_0_0_0_yz_yy_y_yz[i] = -g_z_yy_y_yz[i] + 2.0 * g_yyz_yy_y_yz[i] * a_exp;

        g_y_0_0_0_yz_yy_y_zz[i] = -g_z_yy_y_zz[i] + 2.0 * g_yyz_yy_y_zz[i] * a_exp;
    }
    // integrals block (1146-1152)

    #pragma omp simd aligned(g_y_0_0_0_yz_yy_z_xx, g_y_0_0_0_yz_yy_z_xy, g_y_0_0_0_yz_yy_z_xz, g_y_0_0_0_yz_yy_z_yy, g_y_0_0_0_yz_yy_z_yz, g_y_0_0_0_yz_yy_z_zz, g_yyz_yy_z_xx, g_yyz_yy_z_xy, g_yyz_yy_z_xz, g_yyz_yy_z_yy, g_yyz_yy_z_yz, g_yyz_yy_z_zz, g_z_yy_z_xx, g_z_yy_z_xy, g_z_yy_z_xz, g_z_yy_z_yy, g_z_yy_z_yz, g_z_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_yy_z_xx[i] = -g_z_yy_z_xx[i] + 2.0 * g_yyz_yy_z_xx[i] * a_exp;

        g_y_0_0_0_yz_yy_z_xy[i] = -g_z_yy_z_xy[i] + 2.0 * g_yyz_yy_z_xy[i] * a_exp;

        g_y_0_0_0_yz_yy_z_xz[i] = -g_z_yy_z_xz[i] + 2.0 * g_yyz_yy_z_xz[i] * a_exp;

        g_y_0_0_0_yz_yy_z_yy[i] = -g_z_yy_z_yy[i] + 2.0 * g_yyz_yy_z_yy[i] * a_exp;

        g_y_0_0_0_yz_yy_z_yz[i] = -g_z_yy_z_yz[i] + 2.0 * g_yyz_yy_z_yz[i] * a_exp;

        g_y_0_0_0_yz_yy_z_zz[i] = -g_z_yy_z_zz[i] + 2.0 * g_yyz_yy_z_zz[i] * a_exp;
    }
    // integrals block (1152-1158)

    #pragma omp simd aligned(g_y_0_0_0_yz_yz_x_xx, g_y_0_0_0_yz_yz_x_xy, g_y_0_0_0_yz_yz_x_xz, g_y_0_0_0_yz_yz_x_yy, g_y_0_0_0_yz_yz_x_yz, g_y_0_0_0_yz_yz_x_zz, g_yyz_yz_x_xx, g_yyz_yz_x_xy, g_yyz_yz_x_xz, g_yyz_yz_x_yy, g_yyz_yz_x_yz, g_yyz_yz_x_zz, g_z_yz_x_xx, g_z_yz_x_xy, g_z_yz_x_xz, g_z_yz_x_yy, g_z_yz_x_yz, g_z_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_yz_x_xx[i] = -g_z_yz_x_xx[i] + 2.0 * g_yyz_yz_x_xx[i] * a_exp;

        g_y_0_0_0_yz_yz_x_xy[i] = -g_z_yz_x_xy[i] + 2.0 * g_yyz_yz_x_xy[i] * a_exp;

        g_y_0_0_0_yz_yz_x_xz[i] = -g_z_yz_x_xz[i] + 2.0 * g_yyz_yz_x_xz[i] * a_exp;

        g_y_0_0_0_yz_yz_x_yy[i] = -g_z_yz_x_yy[i] + 2.0 * g_yyz_yz_x_yy[i] * a_exp;

        g_y_0_0_0_yz_yz_x_yz[i] = -g_z_yz_x_yz[i] + 2.0 * g_yyz_yz_x_yz[i] * a_exp;

        g_y_0_0_0_yz_yz_x_zz[i] = -g_z_yz_x_zz[i] + 2.0 * g_yyz_yz_x_zz[i] * a_exp;
    }
    // integrals block (1158-1164)

    #pragma omp simd aligned(g_y_0_0_0_yz_yz_y_xx, g_y_0_0_0_yz_yz_y_xy, g_y_0_0_0_yz_yz_y_xz, g_y_0_0_0_yz_yz_y_yy, g_y_0_0_0_yz_yz_y_yz, g_y_0_0_0_yz_yz_y_zz, g_yyz_yz_y_xx, g_yyz_yz_y_xy, g_yyz_yz_y_xz, g_yyz_yz_y_yy, g_yyz_yz_y_yz, g_yyz_yz_y_zz, g_z_yz_y_xx, g_z_yz_y_xy, g_z_yz_y_xz, g_z_yz_y_yy, g_z_yz_y_yz, g_z_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_yz_y_xx[i] = -g_z_yz_y_xx[i] + 2.0 * g_yyz_yz_y_xx[i] * a_exp;

        g_y_0_0_0_yz_yz_y_xy[i] = -g_z_yz_y_xy[i] + 2.0 * g_yyz_yz_y_xy[i] * a_exp;

        g_y_0_0_0_yz_yz_y_xz[i] = -g_z_yz_y_xz[i] + 2.0 * g_yyz_yz_y_xz[i] * a_exp;

        g_y_0_0_0_yz_yz_y_yy[i] = -g_z_yz_y_yy[i] + 2.0 * g_yyz_yz_y_yy[i] * a_exp;

        g_y_0_0_0_yz_yz_y_yz[i] = -g_z_yz_y_yz[i] + 2.0 * g_yyz_yz_y_yz[i] * a_exp;

        g_y_0_0_0_yz_yz_y_zz[i] = -g_z_yz_y_zz[i] + 2.0 * g_yyz_yz_y_zz[i] * a_exp;
    }
    // integrals block (1164-1170)

    #pragma omp simd aligned(g_y_0_0_0_yz_yz_z_xx, g_y_0_0_0_yz_yz_z_xy, g_y_0_0_0_yz_yz_z_xz, g_y_0_0_0_yz_yz_z_yy, g_y_0_0_0_yz_yz_z_yz, g_y_0_0_0_yz_yz_z_zz, g_yyz_yz_z_xx, g_yyz_yz_z_xy, g_yyz_yz_z_xz, g_yyz_yz_z_yy, g_yyz_yz_z_yz, g_yyz_yz_z_zz, g_z_yz_z_xx, g_z_yz_z_xy, g_z_yz_z_xz, g_z_yz_z_yy, g_z_yz_z_yz, g_z_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_yz_z_xx[i] = -g_z_yz_z_xx[i] + 2.0 * g_yyz_yz_z_xx[i] * a_exp;

        g_y_0_0_0_yz_yz_z_xy[i] = -g_z_yz_z_xy[i] + 2.0 * g_yyz_yz_z_xy[i] * a_exp;

        g_y_0_0_0_yz_yz_z_xz[i] = -g_z_yz_z_xz[i] + 2.0 * g_yyz_yz_z_xz[i] * a_exp;

        g_y_0_0_0_yz_yz_z_yy[i] = -g_z_yz_z_yy[i] + 2.0 * g_yyz_yz_z_yy[i] * a_exp;

        g_y_0_0_0_yz_yz_z_yz[i] = -g_z_yz_z_yz[i] + 2.0 * g_yyz_yz_z_yz[i] * a_exp;

        g_y_0_0_0_yz_yz_z_zz[i] = -g_z_yz_z_zz[i] + 2.0 * g_yyz_yz_z_zz[i] * a_exp;
    }
    // integrals block (1170-1176)

    #pragma omp simd aligned(g_y_0_0_0_yz_zz_x_xx, g_y_0_0_0_yz_zz_x_xy, g_y_0_0_0_yz_zz_x_xz, g_y_0_0_0_yz_zz_x_yy, g_y_0_0_0_yz_zz_x_yz, g_y_0_0_0_yz_zz_x_zz, g_yyz_zz_x_xx, g_yyz_zz_x_xy, g_yyz_zz_x_xz, g_yyz_zz_x_yy, g_yyz_zz_x_yz, g_yyz_zz_x_zz, g_z_zz_x_xx, g_z_zz_x_xy, g_z_zz_x_xz, g_z_zz_x_yy, g_z_zz_x_yz, g_z_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_zz_x_xx[i] = -g_z_zz_x_xx[i] + 2.0 * g_yyz_zz_x_xx[i] * a_exp;

        g_y_0_0_0_yz_zz_x_xy[i] = -g_z_zz_x_xy[i] + 2.0 * g_yyz_zz_x_xy[i] * a_exp;

        g_y_0_0_0_yz_zz_x_xz[i] = -g_z_zz_x_xz[i] + 2.0 * g_yyz_zz_x_xz[i] * a_exp;

        g_y_0_0_0_yz_zz_x_yy[i] = -g_z_zz_x_yy[i] + 2.0 * g_yyz_zz_x_yy[i] * a_exp;

        g_y_0_0_0_yz_zz_x_yz[i] = -g_z_zz_x_yz[i] + 2.0 * g_yyz_zz_x_yz[i] * a_exp;

        g_y_0_0_0_yz_zz_x_zz[i] = -g_z_zz_x_zz[i] + 2.0 * g_yyz_zz_x_zz[i] * a_exp;
    }
    // integrals block (1176-1182)

    #pragma omp simd aligned(g_y_0_0_0_yz_zz_y_xx, g_y_0_0_0_yz_zz_y_xy, g_y_0_0_0_yz_zz_y_xz, g_y_0_0_0_yz_zz_y_yy, g_y_0_0_0_yz_zz_y_yz, g_y_0_0_0_yz_zz_y_zz, g_yyz_zz_y_xx, g_yyz_zz_y_xy, g_yyz_zz_y_xz, g_yyz_zz_y_yy, g_yyz_zz_y_yz, g_yyz_zz_y_zz, g_z_zz_y_xx, g_z_zz_y_xy, g_z_zz_y_xz, g_z_zz_y_yy, g_z_zz_y_yz, g_z_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_zz_y_xx[i] = -g_z_zz_y_xx[i] + 2.0 * g_yyz_zz_y_xx[i] * a_exp;

        g_y_0_0_0_yz_zz_y_xy[i] = -g_z_zz_y_xy[i] + 2.0 * g_yyz_zz_y_xy[i] * a_exp;

        g_y_0_0_0_yz_zz_y_xz[i] = -g_z_zz_y_xz[i] + 2.0 * g_yyz_zz_y_xz[i] * a_exp;

        g_y_0_0_0_yz_zz_y_yy[i] = -g_z_zz_y_yy[i] + 2.0 * g_yyz_zz_y_yy[i] * a_exp;

        g_y_0_0_0_yz_zz_y_yz[i] = -g_z_zz_y_yz[i] + 2.0 * g_yyz_zz_y_yz[i] * a_exp;

        g_y_0_0_0_yz_zz_y_zz[i] = -g_z_zz_y_zz[i] + 2.0 * g_yyz_zz_y_zz[i] * a_exp;
    }
    // integrals block (1182-1188)

    #pragma omp simd aligned(g_y_0_0_0_yz_zz_z_xx, g_y_0_0_0_yz_zz_z_xy, g_y_0_0_0_yz_zz_z_xz, g_y_0_0_0_yz_zz_z_yy, g_y_0_0_0_yz_zz_z_yz, g_y_0_0_0_yz_zz_z_zz, g_yyz_zz_z_xx, g_yyz_zz_z_xy, g_yyz_zz_z_xz, g_yyz_zz_z_yy, g_yyz_zz_z_yz, g_yyz_zz_z_zz, g_z_zz_z_xx, g_z_zz_z_xy, g_z_zz_z_xz, g_z_zz_z_yy, g_z_zz_z_yz, g_z_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_zz_z_xx[i] = -g_z_zz_z_xx[i] + 2.0 * g_yyz_zz_z_xx[i] * a_exp;

        g_y_0_0_0_yz_zz_z_xy[i] = -g_z_zz_z_xy[i] + 2.0 * g_yyz_zz_z_xy[i] * a_exp;

        g_y_0_0_0_yz_zz_z_xz[i] = -g_z_zz_z_xz[i] + 2.0 * g_yyz_zz_z_xz[i] * a_exp;

        g_y_0_0_0_yz_zz_z_yy[i] = -g_z_zz_z_yy[i] + 2.0 * g_yyz_zz_z_yy[i] * a_exp;

        g_y_0_0_0_yz_zz_z_yz[i] = -g_z_zz_z_yz[i] + 2.0 * g_yyz_zz_z_yz[i] * a_exp;

        g_y_0_0_0_yz_zz_z_zz[i] = -g_z_zz_z_zz[i] + 2.0 * g_yyz_zz_z_zz[i] * a_exp;
    }
    // integrals block (1188-1194)

    #pragma omp simd aligned(g_y_0_0_0_zz_xx_x_xx, g_y_0_0_0_zz_xx_x_xy, g_y_0_0_0_zz_xx_x_xz, g_y_0_0_0_zz_xx_x_yy, g_y_0_0_0_zz_xx_x_yz, g_y_0_0_0_zz_xx_x_zz, g_yzz_xx_x_xx, g_yzz_xx_x_xy, g_yzz_xx_x_xz, g_yzz_xx_x_yy, g_yzz_xx_x_yz, g_yzz_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_xx_x_xx[i] = 2.0 * g_yzz_xx_x_xx[i] * a_exp;

        g_y_0_0_0_zz_xx_x_xy[i] = 2.0 * g_yzz_xx_x_xy[i] * a_exp;

        g_y_0_0_0_zz_xx_x_xz[i] = 2.0 * g_yzz_xx_x_xz[i] * a_exp;

        g_y_0_0_0_zz_xx_x_yy[i] = 2.0 * g_yzz_xx_x_yy[i] * a_exp;

        g_y_0_0_0_zz_xx_x_yz[i] = 2.0 * g_yzz_xx_x_yz[i] * a_exp;

        g_y_0_0_0_zz_xx_x_zz[i] = 2.0 * g_yzz_xx_x_zz[i] * a_exp;
    }
    // integrals block (1194-1200)

    #pragma omp simd aligned(g_y_0_0_0_zz_xx_y_xx, g_y_0_0_0_zz_xx_y_xy, g_y_0_0_0_zz_xx_y_xz, g_y_0_0_0_zz_xx_y_yy, g_y_0_0_0_zz_xx_y_yz, g_y_0_0_0_zz_xx_y_zz, g_yzz_xx_y_xx, g_yzz_xx_y_xy, g_yzz_xx_y_xz, g_yzz_xx_y_yy, g_yzz_xx_y_yz, g_yzz_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_xx_y_xx[i] = 2.0 * g_yzz_xx_y_xx[i] * a_exp;

        g_y_0_0_0_zz_xx_y_xy[i] = 2.0 * g_yzz_xx_y_xy[i] * a_exp;

        g_y_0_0_0_zz_xx_y_xz[i] = 2.0 * g_yzz_xx_y_xz[i] * a_exp;

        g_y_0_0_0_zz_xx_y_yy[i] = 2.0 * g_yzz_xx_y_yy[i] * a_exp;

        g_y_0_0_0_zz_xx_y_yz[i] = 2.0 * g_yzz_xx_y_yz[i] * a_exp;

        g_y_0_0_0_zz_xx_y_zz[i] = 2.0 * g_yzz_xx_y_zz[i] * a_exp;
    }
    // integrals block (1200-1206)

    #pragma omp simd aligned(g_y_0_0_0_zz_xx_z_xx, g_y_0_0_0_zz_xx_z_xy, g_y_0_0_0_zz_xx_z_xz, g_y_0_0_0_zz_xx_z_yy, g_y_0_0_0_zz_xx_z_yz, g_y_0_0_0_zz_xx_z_zz, g_yzz_xx_z_xx, g_yzz_xx_z_xy, g_yzz_xx_z_xz, g_yzz_xx_z_yy, g_yzz_xx_z_yz, g_yzz_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_xx_z_xx[i] = 2.0 * g_yzz_xx_z_xx[i] * a_exp;

        g_y_0_0_0_zz_xx_z_xy[i] = 2.0 * g_yzz_xx_z_xy[i] * a_exp;

        g_y_0_0_0_zz_xx_z_xz[i] = 2.0 * g_yzz_xx_z_xz[i] * a_exp;

        g_y_0_0_0_zz_xx_z_yy[i] = 2.0 * g_yzz_xx_z_yy[i] * a_exp;

        g_y_0_0_0_zz_xx_z_yz[i] = 2.0 * g_yzz_xx_z_yz[i] * a_exp;

        g_y_0_0_0_zz_xx_z_zz[i] = 2.0 * g_yzz_xx_z_zz[i] * a_exp;
    }
    // integrals block (1206-1212)

    #pragma omp simd aligned(g_y_0_0_0_zz_xy_x_xx, g_y_0_0_0_zz_xy_x_xy, g_y_0_0_0_zz_xy_x_xz, g_y_0_0_0_zz_xy_x_yy, g_y_0_0_0_zz_xy_x_yz, g_y_0_0_0_zz_xy_x_zz, g_yzz_xy_x_xx, g_yzz_xy_x_xy, g_yzz_xy_x_xz, g_yzz_xy_x_yy, g_yzz_xy_x_yz, g_yzz_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_xy_x_xx[i] = 2.0 * g_yzz_xy_x_xx[i] * a_exp;

        g_y_0_0_0_zz_xy_x_xy[i] = 2.0 * g_yzz_xy_x_xy[i] * a_exp;

        g_y_0_0_0_zz_xy_x_xz[i] = 2.0 * g_yzz_xy_x_xz[i] * a_exp;

        g_y_0_0_0_zz_xy_x_yy[i] = 2.0 * g_yzz_xy_x_yy[i] * a_exp;

        g_y_0_0_0_zz_xy_x_yz[i] = 2.0 * g_yzz_xy_x_yz[i] * a_exp;

        g_y_0_0_0_zz_xy_x_zz[i] = 2.0 * g_yzz_xy_x_zz[i] * a_exp;
    }
    // integrals block (1212-1218)

    #pragma omp simd aligned(g_y_0_0_0_zz_xy_y_xx, g_y_0_0_0_zz_xy_y_xy, g_y_0_0_0_zz_xy_y_xz, g_y_0_0_0_zz_xy_y_yy, g_y_0_0_0_zz_xy_y_yz, g_y_0_0_0_zz_xy_y_zz, g_yzz_xy_y_xx, g_yzz_xy_y_xy, g_yzz_xy_y_xz, g_yzz_xy_y_yy, g_yzz_xy_y_yz, g_yzz_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_xy_y_xx[i] = 2.0 * g_yzz_xy_y_xx[i] * a_exp;

        g_y_0_0_0_zz_xy_y_xy[i] = 2.0 * g_yzz_xy_y_xy[i] * a_exp;

        g_y_0_0_0_zz_xy_y_xz[i] = 2.0 * g_yzz_xy_y_xz[i] * a_exp;

        g_y_0_0_0_zz_xy_y_yy[i] = 2.0 * g_yzz_xy_y_yy[i] * a_exp;

        g_y_0_0_0_zz_xy_y_yz[i] = 2.0 * g_yzz_xy_y_yz[i] * a_exp;

        g_y_0_0_0_zz_xy_y_zz[i] = 2.0 * g_yzz_xy_y_zz[i] * a_exp;
    }
    // integrals block (1218-1224)

    #pragma omp simd aligned(g_y_0_0_0_zz_xy_z_xx, g_y_0_0_0_zz_xy_z_xy, g_y_0_0_0_zz_xy_z_xz, g_y_0_0_0_zz_xy_z_yy, g_y_0_0_0_zz_xy_z_yz, g_y_0_0_0_zz_xy_z_zz, g_yzz_xy_z_xx, g_yzz_xy_z_xy, g_yzz_xy_z_xz, g_yzz_xy_z_yy, g_yzz_xy_z_yz, g_yzz_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_xy_z_xx[i] = 2.0 * g_yzz_xy_z_xx[i] * a_exp;

        g_y_0_0_0_zz_xy_z_xy[i] = 2.0 * g_yzz_xy_z_xy[i] * a_exp;

        g_y_0_0_0_zz_xy_z_xz[i] = 2.0 * g_yzz_xy_z_xz[i] * a_exp;

        g_y_0_0_0_zz_xy_z_yy[i] = 2.0 * g_yzz_xy_z_yy[i] * a_exp;

        g_y_0_0_0_zz_xy_z_yz[i] = 2.0 * g_yzz_xy_z_yz[i] * a_exp;

        g_y_0_0_0_zz_xy_z_zz[i] = 2.0 * g_yzz_xy_z_zz[i] * a_exp;
    }
    // integrals block (1224-1230)

    #pragma omp simd aligned(g_y_0_0_0_zz_xz_x_xx, g_y_0_0_0_zz_xz_x_xy, g_y_0_0_0_zz_xz_x_xz, g_y_0_0_0_zz_xz_x_yy, g_y_0_0_0_zz_xz_x_yz, g_y_0_0_0_zz_xz_x_zz, g_yzz_xz_x_xx, g_yzz_xz_x_xy, g_yzz_xz_x_xz, g_yzz_xz_x_yy, g_yzz_xz_x_yz, g_yzz_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_xz_x_xx[i] = 2.0 * g_yzz_xz_x_xx[i] * a_exp;

        g_y_0_0_0_zz_xz_x_xy[i] = 2.0 * g_yzz_xz_x_xy[i] * a_exp;

        g_y_0_0_0_zz_xz_x_xz[i] = 2.0 * g_yzz_xz_x_xz[i] * a_exp;

        g_y_0_0_0_zz_xz_x_yy[i] = 2.0 * g_yzz_xz_x_yy[i] * a_exp;

        g_y_0_0_0_zz_xz_x_yz[i] = 2.0 * g_yzz_xz_x_yz[i] * a_exp;

        g_y_0_0_0_zz_xz_x_zz[i] = 2.0 * g_yzz_xz_x_zz[i] * a_exp;
    }
    // integrals block (1230-1236)

    #pragma omp simd aligned(g_y_0_0_0_zz_xz_y_xx, g_y_0_0_0_zz_xz_y_xy, g_y_0_0_0_zz_xz_y_xz, g_y_0_0_0_zz_xz_y_yy, g_y_0_0_0_zz_xz_y_yz, g_y_0_0_0_zz_xz_y_zz, g_yzz_xz_y_xx, g_yzz_xz_y_xy, g_yzz_xz_y_xz, g_yzz_xz_y_yy, g_yzz_xz_y_yz, g_yzz_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_xz_y_xx[i] = 2.0 * g_yzz_xz_y_xx[i] * a_exp;

        g_y_0_0_0_zz_xz_y_xy[i] = 2.0 * g_yzz_xz_y_xy[i] * a_exp;

        g_y_0_0_0_zz_xz_y_xz[i] = 2.0 * g_yzz_xz_y_xz[i] * a_exp;

        g_y_0_0_0_zz_xz_y_yy[i] = 2.0 * g_yzz_xz_y_yy[i] * a_exp;

        g_y_0_0_0_zz_xz_y_yz[i] = 2.0 * g_yzz_xz_y_yz[i] * a_exp;

        g_y_0_0_0_zz_xz_y_zz[i] = 2.0 * g_yzz_xz_y_zz[i] * a_exp;
    }
    // integrals block (1236-1242)

    #pragma omp simd aligned(g_y_0_0_0_zz_xz_z_xx, g_y_0_0_0_zz_xz_z_xy, g_y_0_0_0_zz_xz_z_xz, g_y_0_0_0_zz_xz_z_yy, g_y_0_0_0_zz_xz_z_yz, g_y_0_0_0_zz_xz_z_zz, g_yzz_xz_z_xx, g_yzz_xz_z_xy, g_yzz_xz_z_xz, g_yzz_xz_z_yy, g_yzz_xz_z_yz, g_yzz_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_xz_z_xx[i] = 2.0 * g_yzz_xz_z_xx[i] * a_exp;

        g_y_0_0_0_zz_xz_z_xy[i] = 2.0 * g_yzz_xz_z_xy[i] * a_exp;

        g_y_0_0_0_zz_xz_z_xz[i] = 2.0 * g_yzz_xz_z_xz[i] * a_exp;

        g_y_0_0_0_zz_xz_z_yy[i] = 2.0 * g_yzz_xz_z_yy[i] * a_exp;

        g_y_0_0_0_zz_xz_z_yz[i] = 2.0 * g_yzz_xz_z_yz[i] * a_exp;

        g_y_0_0_0_zz_xz_z_zz[i] = 2.0 * g_yzz_xz_z_zz[i] * a_exp;
    }
    // integrals block (1242-1248)

    #pragma omp simd aligned(g_y_0_0_0_zz_yy_x_xx, g_y_0_0_0_zz_yy_x_xy, g_y_0_0_0_zz_yy_x_xz, g_y_0_0_0_zz_yy_x_yy, g_y_0_0_0_zz_yy_x_yz, g_y_0_0_0_zz_yy_x_zz, g_yzz_yy_x_xx, g_yzz_yy_x_xy, g_yzz_yy_x_xz, g_yzz_yy_x_yy, g_yzz_yy_x_yz, g_yzz_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_yy_x_xx[i] = 2.0 * g_yzz_yy_x_xx[i] * a_exp;

        g_y_0_0_0_zz_yy_x_xy[i] = 2.0 * g_yzz_yy_x_xy[i] * a_exp;

        g_y_0_0_0_zz_yy_x_xz[i] = 2.0 * g_yzz_yy_x_xz[i] * a_exp;

        g_y_0_0_0_zz_yy_x_yy[i] = 2.0 * g_yzz_yy_x_yy[i] * a_exp;

        g_y_0_0_0_zz_yy_x_yz[i] = 2.0 * g_yzz_yy_x_yz[i] * a_exp;

        g_y_0_0_0_zz_yy_x_zz[i] = 2.0 * g_yzz_yy_x_zz[i] * a_exp;
    }
    // integrals block (1248-1254)

    #pragma omp simd aligned(g_y_0_0_0_zz_yy_y_xx, g_y_0_0_0_zz_yy_y_xy, g_y_0_0_0_zz_yy_y_xz, g_y_0_0_0_zz_yy_y_yy, g_y_0_0_0_zz_yy_y_yz, g_y_0_0_0_zz_yy_y_zz, g_yzz_yy_y_xx, g_yzz_yy_y_xy, g_yzz_yy_y_xz, g_yzz_yy_y_yy, g_yzz_yy_y_yz, g_yzz_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_yy_y_xx[i] = 2.0 * g_yzz_yy_y_xx[i] * a_exp;

        g_y_0_0_0_zz_yy_y_xy[i] = 2.0 * g_yzz_yy_y_xy[i] * a_exp;

        g_y_0_0_0_zz_yy_y_xz[i] = 2.0 * g_yzz_yy_y_xz[i] * a_exp;

        g_y_0_0_0_zz_yy_y_yy[i] = 2.0 * g_yzz_yy_y_yy[i] * a_exp;

        g_y_0_0_0_zz_yy_y_yz[i] = 2.0 * g_yzz_yy_y_yz[i] * a_exp;

        g_y_0_0_0_zz_yy_y_zz[i] = 2.0 * g_yzz_yy_y_zz[i] * a_exp;
    }
    // integrals block (1254-1260)

    #pragma omp simd aligned(g_y_0_0_0_zz_yy_z_xx, g_y_0_0_0_zz_yy_z_xy, g_y_0_0_0_zz_yy_z_xz, g_y_0_0_0_zz_yy_z_yy, g_y_0_0_0_zz_yy_z_yz, g_y_0_0_0_zz_yy_z_zz, g_yzz_yy_z_xx, g_yzz_yy_z_xy, g_yzz_yy_z_xz, g_yzz_yy_z_yy, g_yzz_yy_z_yz, g_yzz_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_yy_z_xx[i] = 2.0 * g_yzz_yy_z_xx[i] * a_exp;

        g_y_0_0_0_zz_yy_z_xy[i] = 2.0 * g_yzz_yy_z_xy[i] * a_exp;

        g_y_0_0_0_zz_yy_z_xz[i] = 2.0 * g_yzz_yy_z_xz[i] * a_exp;

        g_y_0_0_0_zz_yy_z_yy[i] = 2.0 * g_yzz_yy_z_yy[i] * a_exp;

        g_y_0_0_0_zz_yy_z_yz[i] = 2.0 * g_yzz_yy_z_yz[i] * a_exp;

        g_y_0_0_0_zz_yy_z_zz[i] = 2.0 * g_yzz_yy_z_zz[i] * a_exp;
    }
    // integrals block (1260-1266)

    #pragma omp simd aligned(g_y_0_0_0_zz_yz_x_xx, g_y_0_0_0_zz_yz_x_xy, g_y_0_0_0_zz_yz_x_xz, g_y_0_0_0_zz_yz_x_yy, g_y_0_0_0_zz_yz_x_yz, g_y_0_0_0_zz_yz_x_zz, g_yzz_yz_x_xx, g_yzz_yz_x_xy, g_yzz_yz_x_xz, g_yzz_yz_x_yy, g_yzz_yz_x_yz, g_yzz_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_yz_x_xx[i] = 2.0 * g_yzz_yz_x_xx[i] * a_exp;

        g_y_0_0_0_zz_yz_x_xy[i] = 2.0 * g_yzz_yz_x_xy[i] * a_exp;

        g_y_0_0_0_zz_yz_x_xz[i] = 2.0 * g_yzz_yz_x_xz[i] * a_exp;

        g_y_0_0_0_zz_yz_x_yy[i] = 2.0 * g_yzz_yz_x_yy[i] * a_exp;

        g_y_0_0_0_zz_yz_x_yz[i] = 2.0 * g_yzz_yz_x_yz[i] * a_exp;

        g_y_0_0_0_zz_yz_x_zz[i] = 2.0 * g_yzz_yz_x_zz[i] * a_exp;
    }
    // integrals block (1266-1272)

    #pragma omp simd aligned(g_y_0_0_0_zz_yz_y_xx, g_y_0_0_0_zz_yz_y_xy, g_y_0_0_0_zz_yz_y_xz, g_y_0_0_0_zz_yz_y_yy, g_y_0_0_0_zz_yz_y_yz, g_y_0_0_0_zz_yz_y_zz, g_yzz_yz_y_xx, g_yzz_yz_y_xy, g_yzz_yz_y_xz, g_yzz_yz_y_yy, g_yzz_yz_y_yz, g_yzz_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_yz_y_xx[i] = 2.0 * g_yzz_yz_y_xx[i] * a_exp;

        g_y_0_0_0_zz_yz_y_xy[i] = 2.0 * g_yzz_yz_y_xy[i] * a_exp;

        g_y_0_0_0_zz_yz_y_xz[i] = 2.0 * g_yzz_yz_y_xz[i] * a_exp;

        g_y_0_0_0_zz_yz_y_yy[i] = 2.0 * g_yzz_yz_y_yy[i] * a_exp;

        g_y_0_0_0_zz_yz_y_yz[i] = 2.0 * g_yzz_yz_y_yz[i] * a_exp;

        g_y_0_0_0_zz_yz_y_zz[i] = 2.0 * g_yzz_yz_y_zz[i] * a_exp;
    }
    // integrals block (1272-1278)

    #pragma omp simd aligned(g_y_0_0_0_zz_yz_z_xx, g_y_0_0_0_zz_yz_z_xy, g_y_0_0_0_zz_yz_z_xz, g_y_0_0_0_zz_yz_z_yy, g_y_0_0_0_zz_yz_z_yz, g_y_0_0_0_zz_yz_z_zz, g_yzz_yz_z_xx, g_yzz_yz_z_xy, g_yzz_yz_z_xz, g_yzz_yz_z_yy, g_yzz_yz_z_yz, g_yzz_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_yz_z_xx[i] = 2.0 * g_yzz_yz_z_xx[i] * a_exp;

        g_y_0_0_0_zz_yz_z_xy[i] = 2.0 * g_yzz_yz_z_xy[i] * a_exp;

        g_y_0_0_0_zz_yz_z_xz[i] = 2.0 * g_yzz_yz_z_xz[i] * a_exp;

        g_y_0_0_0_zz_yz_z_yy[i] = 2.0 * g_yzz_yz_z_yy[i] * a_exp;

        g_y_0_0_0_zz_yz_z_yz[i] = 2.0 * g_yzz_yz_z_yz[i] * a_exp;

        g_y_0_0_0_zz_yz_z_zz[i] = 2.0 * g_yzz_yz_z_zz[i] * a_exp;
    }
    // integrals block (1278-1284)

    #pragma omp simd aligned(g_y_0_0_0_zz_zz_x_xx, g_y_0_0_0_zz_zz_x_xy, g_y_0_0_0_zz_zz_x_xz, g_y_0_0_0_zz_zz_x_yy, g_y_0_0_0_zz_zz_x_yz, g_y_0_0_0_zz_zz_x_zz, g_yzz_zz_x_xx, g_yzz_zz_x_xy, g_yzz_zz_x_xz, g_yzz_zz_x_yy, g_yzz_zz_x_yz, g_yzz_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_zz_x_xx[i] = 2.0 * g_yzz_zz_x_xx[i] * a_exp;

        g_y_0_0_0_zz_zz_x_xy[i] = 2.0 * g_yzz_zz_x_xy[i] * a_exp;

        g_y_0_0_0_zz_zz_x_xz[i] = 2.0 * g_yzz_zz_x_xz[i] * a_exp;

        g_y_0_0_0_zz_zz_x_yy[i] = 2.0 * g_yzz_zz_x_yy[i] * a_exp;

        g_y_0_0_0_zz_zz_x_yz[i] = 2.0 * g_yzz_zz_x_yz[i] * a_exp;

        g_y_0_0_0_zz_zz_x_zz[i] = 2.0 * g_yzz_zz_x_zz[i] * a_exp;
    }
    // integrals block (1284-1290)

    #pragma omp simd aligned(g_y_0_0_0_zz_zz_y_xx, g_y_0_0_0_zz_zz_y_xy, g_y_0_0_0_zz_zz_y_xz, g_y_0_0_0_zz_zz_y_yy, g_y_0_0_0_zz_zz_y_yz, g_y_0_0_0_zz_zz_y_zz, g_yzz_zz_y_xx, g_yzz_zz_y_xy, g_yzz_zz_y_xz, g_yzz_zz_y_yy, g_yzz_zz_y_yz, g_yzz_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_zz_y_xx[i] = 2.0 * g_yzz_zz_y_xx[i] * a_exp;

        g_y_0_0_0_zz_zz_y_xy[i] = 2.0 * g_yzz_zz_y_xy[i] * a_exp;

        g_y_0_0_0_zz_zz_y_xz[i] = 2.0 * g_yzz_zz_y_xz[i] * a_exp;

        g_y_0_0_0_zz_zz_y_yy[i] = 2.0 * g_yzz_zz_y_yy[i] * a_exp;

        g_y_0_0_0_zz_zz_y_yz[i] = 2.0 * g_yzz_zz_y_yz[i] * a_exp;

        g_y_0_0_0_zz_zz_y_zz[i] = 2.0 * g_yzz_zz_y_zz[i] * a_exp;
    }
    // integrals block (1290-1296)

    #pragma omp simd aligned(g_y_0_0_0_zz_zz_z_xx, g_y_0_0_0_zz_zz_z_xy, g_y_0_0_0_zz_zz_z_xz, g_y_0_0_0_zz_zz_z_yy, g_y_0_0_0_zz_zz_z_yz, g_y_0_0_0_zz_zz_z_zz, g_yzz_zz_z_xx, g_yzz_zz_z_xy, g_yzz_zz_z_xz, g_yzz_zz_z_yy, g_yzz_zz_z_yz, g_yzz_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_zz_z_xx[i] = 2.0 * g_yzz_zz_z_xx[i] * a_exp;

        g_y_0_0_0_zz_zz_z_xy[i] = 2.0 * g_yzz_zz_z_xy[i] * a_exp;

        g_y_0_0_0_zz_zz_z_xz[i] = 2.0 * g_yzz_zz_z_xz[i] * a_exp;

        g_y_0_0_0_zz_zz_z_yy[i] = 2.0 * g_yzz_zz_z_yy[i] * a_exp;

        g_y_0_0_0_zz_zz_z_yz[i] = 2.0 * g_yzz_zz_z_yz[i] * a_exp;

        g_y_0_0_0_zz_zz_z_zz[i] = 2.0 * g_yzz_zz_z_zz[i] * a_exp;
    }
    // integrals block (1296-1302)

    #pragma omp simd aligned(g_xxz_xx_x_xx, g_xxz_xx_x_xy, g_xxz_xx_x_xz, g_xxz_xx_x_yy, g_xxz_xx_x_yz, g_xxz_xx_x_zz, g_z_0_0_0_xx_xx_x_xx, g_z_0_0_0_xx_xx_x_xy, g_z_0_0_0_xx_xx_x_xz, g_z_0_0_0_xx_xx_x_yy, g_z_0_0_0_xx_xx_x_yz, g_z_0_0_0_xx_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_xx_x_xx[i] = 2.0 * g_xxz_xx_x_xx[i] * a_exp;

        g_z_0_0_0_xx_xx_x_xy[i] = 2.0 * g_xxz_xx_x_xy[i] * a_exp;

        g_z_0_0_0_xx_xx_x_xz[i] = 2.0 * g_xxz_xx_x_xz[i] * a_exp;

        g_z_0_0_0_xx_xx_x_yy[i] = 2.0 * g_xxz_xx_x_yy[i] * a_exp;

        g_z_0_0_0_xx_xx_x_yz[i] = 2.0 * g_xxz_xx_x_yz[i] * a_exp;

        g_z_0_0_0_xx_xx_x_zz[i] = 2.0 * g_xxz_xx_x_zz[i] * a_exp;
    }
    // integrals block (1302-1308)

    #pragma omp simd aligned(g_xxz_xx_y_xx, g_xxz_xx_y_xy, g_xxz_xx_y_xz, g_xxz_xx_y_yy, g_xxz_xx_y_yz, g_xxz_xx_y_zz, g_z_0_0_0_xx_xx_y_xx, g_z_0_0_0_xx_xx_y_xy, g_z_0_0_0_xx_xx_y_xz, g_z_0_0_0_xx_xx_y_yy, g_z_0_0_0_xx_xx_y_yz, g_z_0_0_0_xx_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_xx_y_xx[i] = 2.0 * g_xxz_xx_y_xx[i] * a_exp;

        g_z_0_0_0_xx_xx_y_xy[i] = 2.0 * g_xxz_xx_y_xy[i] * a_exp;

        g_z_0_0_0_xx_xx_y_xz[i] = 2.0 * g_xxz_xx_y_xz[i] * a_exp;

        g_z_0_0_0_xx_xx_y_yy[i] = 2.0 * g_xxz_xx_y_yy[i] * a_exp;

        g_z_0_0_0_xx_xx_y_yz[i] = 2.0 * g_xxz_xx_y_yz[i] * a_exp;

        g_z_0_0_0_xx_xx_y_zz[i] = 2.0 * g_xxz_xx_y_zz[i] * a_exp;
    }
    // integrals block (1308-1314)

    #pragma omp simd aligned(g_xxz_xx_z_xx, g_xxz_xx_z_xy, g_xxz_xx_z_xz, g_xxz_xx_z_yy, g_xxz_xx_z_yz, g_xxz_xx_z_zz, g_z_0_0_0_xx_xx_z_xx, g_z_0_0_0_xx_xx_z_xy, g_z_0_0_0_xx_xx_z_xz, g_z_0_0_0_xx_xx_z_yy, g_z_0_0_0_xx_xx_z_yz, g_z_0_0_0_xx_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_xx_z_xx[i] = 2.0 * g_xxz_xx_z_xx[i] * a_exp;

        g_z_0_0_0_xx_xx_z_xy[i] = 2.0 * g_xxz_xx_z_xy[i] * a_exp;

        g_z_0_0_0_xx_xx_z_xz[i] = 2.0 * g_xxz_xx_z_xz[i] * a_exp;

        g_z_0_0_0_xx_xx_z_yy[i] = 2.0 * g_xxz_xx_z_yy[i] * a_exp;

        g_z_0_0_0_xx_xx_z_yz[i] = 2.0 * g_xxz_xx_z_yz[i] * a_exp;

        g_z_0_0_0_xx_xx_z_zz[i] = 2.0 * g_xxz_xx_z_zz[i] * a_exp;
    }
    // integrals block (1314-1320)

    #pragma omp simd aligned(g_xxz_xy_x_xx, g_xxz_xy_x_xy, g_xxz_xy_x_xz, g_xxz_xy_x_yy, g_xxz_xy_x_yz, g_xxz_xy_x_zz, g_z_0_0_0_xx_xy_x_xx, g_z_0_0_0_xx_xy_x_xy, g_z_0_0_0_xx_xy_x_xz, g_z_0_0_0_xx_xy_x_yy, g_z_0_0_0_xx_xy_x_yz, g_z_0_0_0_xx_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_xy_x_xx[i] = 2.0 * g_xxz_xy_x_xx[i] * a_exp;

        g_z_0_0_0_xx_xy_x_xy[i] = 2.0 * g_xxz_xy_x_xy[i] * a_exp;

        g_z_0_0_0_xx_xy_x_xz[i] = 2.0 * g_xxz_xy_x_xz[i] * a_exp;

        g_z_0_0_0_xx_xy_x_yy[i] = 2.0 * g_xxz_xy_x_yy[i] * a_exp;

        g_z_0_0_0_xx_xy_x_yz[i] = 2.0 * g_xxz_xy_x_yz[i] * a_exp;

        g_z_0_0_0_xx_xy_x_zz[i] = 2.0 * g_xxz_xy_x_zz[i] * a_exp;
    }
    // integrals block (1320-1326)

    #pragma omp simd aligned(g_xxz_xy_y_xx, g_xxz_xy_y_xy, g_xxz_xy_y_xz, g_xxz_xy_y_yy, g_xxz_xy_y_yz, g_xxz_xy_y_zz, g_z_0_0_0_xx_xy_y_xx, g_z_0_0_0_xx_xy_y_xy, g_z_0_0_0_xx_xy_y_xz, g_z_0_0_0_xx_xy_y_yy, g_z_0_0_0_xx_xy_y_yz, g_z_0_0_0_xx_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_xy_y_xx[i] = 2.0 * g_xxz_xy_y_xx[i] * a_exp;

        g_z_0_0_0_xx_xy_y_xy[i] = 2.0 * g_xxz_xy_y_xy[i] * a_exp;

        g_z_0_0_0_xx_xy_y_xz[i] = 2.0 * g_xxz_xy_y_xz[i] * a_exp;

        g_z_0_0_0_xx_xy_y_yy[i] = 2.0 * g_xxz_xy_y_yy[i] * a_exp;

        g_z_0_0_0_xx_xy_y_yz[i] = 2.0 * g_xxz_xy_y_yz[i] * a_exp;

        g_z_0_0_0_xx_xy_y_zz[i] = 2.0 * g_xxz_xy_y_zz[i] * a_exp;
    }
    // integrals block (1326-1332)

    #pragma omp simd aligned(g_xxz_xy_z_xx, g_xxz_xy_z_xy, g_xxz_xy_z_xz, g_xxz_xy_z_yy, g_xxz_xy_z_yz, g_xxz_xy_z_zz, g_z_0_0_0_xx_xy_z_xx, g_z_0_0_0_xx_xy_z_xy, g_z_0_0_0_xx_xy_z_xz, g_z_0_0_0_xx_xy_z_yy, g_z_0_0_0_xx_xy_z_yz, g_z_0_0_0_xx_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_xy_z_xx[i] = 2.0 * g_xxz_xy_z_xx[i] * a_exp;

        g_z_0_0_0_xx_xy_z_xy[i] = 2.0 * g_xxz_xy_z_xy[i] * a_exp;

        g_z_0_0_0_xx_xy_z_xz[i] = 2.0 * g_xxz_xy_z_xz[i] * a_exp;

        g_z_0_0_0_xx_xy_z_yy[i] = 2.0 * g_xxz_xy_z_yy[i] * a_exp;

        g_z_0_0_0_xx_xy_z_yz[i] = 2.0 * g_xxz_xy_z_yz[i] * a_exp;

        g_z_0_0_0_xx_xy_z_zz[i] = 2.0 * g_xxz_xy_z_zz[i] * a_exp;
    }
    // integrals block (1332-1338)

    #pragma omp simd aligned(g_xxz_xz_x_xx, g_xxz_xz_x_xy, g_xxz_xz_x_xz, g_xxz_xz_x_yy, g_xxz_xz_x_yz, g_xxz_xz_x_zz, g_z_0_0_0_xx_xz_x_xx, g_z_0_0_0_xx_xz_x_xy, g_z_0_0_0_xx_xz_x_xz, g_z_0_0_0_xx_xz_x_yy, g_z_0_0_0_xx_xz_x_yz, g_z_0_0_0_xx_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_xz_x_xx[i] = 2.0 * g_xxz_xz_x_xx[i] * a_exp;

        g_z_0_0_0_xx_xz_x_xy[i] = 2.0 * g_xxz_xz_x_xy[i] * a_exp;

        g_z_0_0_0_xx_xz_x_xz[i] = 2.0 * g_xxz_xz_x_xz[i] * a_exp;

        g_z_0_0_0_xx_xz_x_yy[i] = 2.0 * g_xxz_xz_x_yy[i] * a_exp;

        g_z_0_0_0_xx_xz_x_yz[i] = 2.0 * g_xxz_xz_x_yz[i] * a_exp;

        g_z_0_0_0_xx_xz_x_zz[i] = 2.0 * g_xxz_xz_x_zz[i] * a_exp;
    }
    // integrals block (1338-1344)

    #pragma omp simd aligned(g_xxz_xz_y_xx, g_xxz_xz_y_xy, g_xxz_xz_y_xz, g_xxz_xz_y_yy, g_xxz_xz_y_yz, g_xxz_xz_y_zz, g_z_0_0_0_xx_xz_y_xx, g_z_0_0_0_xx_xz_y_xy, g_z_0_0_0_xx_xz_y_xz, g_z_0_0_0_xx_xz_y_yy, g_z_0_0_0_xx_xz_y_yz, g_z_0_0_0_xx_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_xz_y_xx[i] = 2.0 * g_xxz_xz_y_xx[i] * a_exp;

        g_z_0_0_0_xx_xz_y_xy[i] = 2.0 * g_xxz_xz_y_xy[i] * a_exp;

        g_z_0_0_0_xx_xz_y_xz[i] = 2.0 * g_xxz_xz_y_xz[i] * a_exp;

        g_z_0_0_0_xx_xz_y_yy[i] = 2.0 * g_xxz_xz_y_yy[i] * a_exp;

        g_z_0_0_0_xx_xz_y_yz[i] = 2.0 * g_xxz_xz_y_yz[i] * a_exp;

        g_z_0_0_0_xx_xz_y_zz[i] = 2.0 * g_xxz_xz_y_zz[i] * a_exp;
    }
    // integrals block (1344-1350)

    #pragma omp simd aligned(g_xxz_xz_z_xx, g_xxz_xz_z_xy, g_xxz_xz_z_xz, g_xxz_xz_z_yy, g_xxz_xz_z_yz, g_xxz_xz_z_zz, g_z_0_0_0_xx_xz_z_xx, g_z_0_0_0_xx_xz_z_xy, g_z_0_0_0_xx_xz_z_xz, g_z_0_0_0_xx_xz_z_yy, g_z_0_0_0_xx_xz_z_yz, g_z_0_0_0_xx_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_xz_z_xx[i] = 2.0 * g_xxz_xz_z_xx[i] * a_exp;

        g_z_0_0_0_xx_xz_z_xy[i] = 2.0 * g_xxz_xz_z_xy[i] * a_exp;

        g_z_0_0_0_xx_xz_z_xz[i] = 2.0 * g_xxz_xz_z_xz[i] * a_exp;

        g_z_0_0_0_xx_xz_z_yy[i] = 2.0 * g_xxz_xz_z_yy[i] * a_exp;

        g_z_0_0_0_xx_xz_z_yz[i] = 2.0 * g_xxz_xz_z_yz[i] * a_exp;

        g_z_0_0_0_xx_xz_z_zz[i] = 2.0 * g_xxz_xz_z_zz[i] * a_exp;
    }
    // integrals block (1350-1356)

    #pragma omp simd aligned(g_xxz_yy_x_xx, g_xxz_yy_x_xy, g_xxz_yy_x_xz, g_xxz_yy_x_yy, g_xxz_yy_x_yz, g_xxz_yy_x_zz, g_z_0_0_0_xx_yy_x_xx, g_z_0_0_0_xx_yy_x_xy, g_z_0_0_0_xx_yy_x_xz, g_z_0_0_0_xx_yy_x_yy, g_z_0_0_0_xx_yy_x_yz, g_z_0_0_0_xx_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_yy_x_xx[i] = 2.0 * g_xxz_yy_x_xx[i] * a_exp;

        g_z_0_0_0_xx_yy_x_xy[i] = 2.0 * g_xxz_yy_x_xy[i] * a_exp;

        g_z_0_0_0_xx_yy_x_xz[i] = 2.0 * g_xxz_yy_x_xz[i] * a_exp;

        g_z_0_0_0_xx_yy_x_yy[i] = 2.0 * g_xxz_yy_x_yy[i] * a_exp;

        g_z_0_0_0_xx_yy_x_yz[i] = 2.0 * g_xxz_yy_x_yz[i] * a_exp;

        g_z_0_0_0_xx_yy_x_zz[i] = 2.0 * g_xxz_yy_x_zz[i] * a_exp;
    }
    // integrals block (1356-1362)

    #pragma omp simd aligned(g_xxz_yy_y_xx, g_xxz_yy_y_xy, g_xxz_yy_y_xz, g_xxz_yy_y_yy, g_xxz_yy_y_yz, g_xxz_yy_y_zz, g_z_0_0_0_xx_yy_y_xx, g_z_0_0_0_xx_yy_y_xy, g_z_0_0_0_xx_yy_y_xz, g_z_0_0_0_xx_yy_y_yy, g_z_0_0_0_xx_yy_y_yz, g_z_0_0_0_xx_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_yy_y_xx[i] = 2.0 * g_xxz_yy_y_xx[i] * a_exp;

        g_z_0_0_0_xx_yy_y_xy[i] = 2.0 * g_xxz_yy_y_xy[i] * a_exp;

        g_z_0_0_0_xx_yy_y_xz[i] = 2.0 * g_xxz_yy_y_xz[i] * a_exp;

        g_z_0_0_0_xx_yy_y_yy[i] = 2.0 * g_xxz_yy_y_yy[i] * a_exp;

        g_z_0_0_0_xx_yy_y_yz[i] = 2.0 * g_xxz_yy_y_yz[i] * a_exp;

        g_z_0_0_0_xx_yy_y_zz[i] = 2.0 * g_xxz_yy_y_zz[i] * a_exp;
    }
    // integrals block (1362-1368)

    #pragma omp simd aligned(g_xxz_yy_z_xx, g_xxz_yy_z_xy, g_xxz_yy_z_xz, g_xxz_yy_z_yy, g_xxz_yy_z_yz, g_xxz_yy_z_zz, g_z_0_0_0_xx_yy_z_xx, g_z_0_0_0_xx_yy_z_xy, g_z_0_0_0_xx_yy_z_xz, g_z_0_0_0_xx_yy_z_yy, g_z_0_0_0_xx_yy_z_yz, g_z_0_0_0_xx_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_yy_z_xx[i] = 2.0 * g_xxz_yy_z_xx[i] * a_exp;

        g_z_0_0_0_xx_yy_z_xy[i] = 2.0 * g_xxz_yy_z_xy[i] * a_exp;

        g_z_0_0_0_xx_yy_z_xz[i] = 2.0 * g_xxz_yy_z_xz[i] * a_exp;

        g_z_0_0_0_xx_yy_z_yy[i] = 2.0 * g_xxz_yy_z_yy[i] * a_exp;

        g_z_0_0_0_xx_yy_z_yz[i] = 2.0 * g_xxz_yy_z_yz[i] * a_exp;

        g_z_0_0_0_xx_yy_z_zz[i] = 2.0 * g_xxz_yy_z_zz[i] * a_exp;
    }
    // integrals block (1368-1374)

    #pragma omp simd aligned(g_xxz_yz_x_xx, g_xxz_yz_x_xy, g_xxz_yz_x_xz, g_xxz_yz_x_yy, g_xxz_yz_x_yz, g_xxz_yz_x_zz, g_z_0_0_0_xx_yz_x_xx, g_z_0_0_0_xx_yz_x_xy, g_z_0_0_0_xx_yz_x_xz, g_z_0_0_0_xx_yz_x_yy, g_z_0_0_0_xx_yz_x_yz, g_z_0_0_0_xx_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_yz_x_xx[i] = 2.0 * g_xxz_yz_x_xx[i] * a_exp;

        g_z_0_0_0_xx_yz_x_xy[i] = 2.0 * g_xxz_yz_x_xy[i] * a_exp;

        g_z_0_0_0_xx_yz_x_xz[i] = 2.0 * g_xxz_yz_x_xz[i] * a_exp;

        g_z_0_0_0_xx_yz_x_yy[i] = 2.0 * g_xxz_yz_x_yy[i] * a_exp;

        g_z_0_0_0_xx_yz_x_yz[i] = 2.0 * g_xxz_yz_x_yz[i] * a_exp;

        g_z_0_0_0_xx_yz_x_zz[i] = 2.0 * g_xxz_yz_x_zz[i] * a_exp;
    }
    // integrals block (1374-1380)

    #pragma omp simd aligned(g_xxz_yz_y_xx, g_xxz_yz_y_xy, g_xxz_yz_y_xz, g_xxz_yz_y_yy, g_xxz_yz_y_yz, g_xxz_yz_y_zz, g_z_0_0_0_xx_yz_y_xx, g_z_0_0_0_xx_yz_y_xy, g_z_0_0_0_xx_yz_y_xz, g_z_0_0_0_xx_yz_y_yy, g_z_0_0_0_xx_yz_y_yz, g_z_0_0_0_xx_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_yz_y_xx[i] = 2.0 * g_xxz_yz_y_xx[i] * a_exp;

        g_z_0_0_0_xx_yz_y_xy[i] = 2.0 * g_xxz_yz_y_xy[i] * a_exp;

        g_z_0_0_0_xx_yz_y_xz[i] = 2.0 * g_xxz_yz_y_xz[i] * a_exp;

        g_z_0_0_0_xx_yz_y_yy[i] = 2.0 * g_xxz_yz_y_yy[i] * a_exp;

        g_z_0_0_0_xx_yz_y_yz[i] = 2.0 * g_xxz_yz_y_yz[i] * a_exp;

        g_z_0_0_0_xx_yz_y_zz[i] = 2.0 * g_xxz_yz_y_zz[i] * a_exp;
    }
    // integrals block (1380-1386)

    #pragma omp simd aligned(g_xxz_yz_z_xx, g_xxz_yz_z_xy, g_xxz_yz_z_xz, g_xxz_yz_z_yy, g_xxz_yz_z_yz, g_xxz_yz_z_zz, g_z_0_0_0_xx_yz_z_xx, g_z_0_0_0_xx_yz_z_xy, g_z_0_0_0_xx_yz_z_xz, g_z_0_0_0_xx_yz_z_yy, g_z_0_0_0_xx_yz_z_yz, g_z_0_0_0_xx_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_yz_z_xx[i] = 2.0 * g_xxz_yz_z_xx[i] * a_exp;

        g_z_0_0_0_xx_yz_z_xy[i] = 2.0 * g_xxz_yz_z_xy[i] * a_exp;

        g_z_0_0_0_xx_yz_z_xz[i] = 2.0 * g_xxz_yz_z_xz[i] * a_exp;

        g_z_0_0_0_xx_yz_z_yy[i] = 2.0 * g_xxz_yz_z_yy[i] * a_exp;

        g_z_0_0_0_xx_yz_z_yz[i] = 2.0 * g_xxz_yz_z_yz[i] * a_exp;

        g_z_0_0_0_xx_yz_z_zz[i] = 2.0 * g_xxz_yz_z_zz[i] * a_exp;
    }
    // integrals block (1386-1392)

    #pragma omp simd aligned(g_xxz_zz_x_xx, g_xxz_zz_x_xy, g_xxz_zz_x_xz, g_xxz_zz_x_yy, g_xxz_zz_x_yz, g_xxz_zz_x_zz, g_z_0_0_0_xx_zz_x_xx, g_z_0_0_0_xx_zz_x_xy, g_z_0_0_0_xx_zz_x_xz, g_z_0_0_0_xx_zz_x_yy, g_z_0_0_0_xx_zz_x_yz, g_z_0_0_0_xx_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_zz_x_xx[i] = 2.0 * g_xxz_zz_x_xx[i] * a_exp;

        g_z_0_0_0_xx_zz_x_xy[i] = 2.0 * g_xxz_zz_x_xy[i] * a_exp;

        g_z_0_0_0_xx_zz_x_xz[i] = 2.0 * g_xxz_zz_x_xz[i] * a_exp;

        g_z_0_0_0_xx_zz_x_yy[i] = 2.0 * g_xxz_zz_x_yy[i] * a_exp;

        g_z_0_0_0_xx_zz_x_yz[i] = 2.0 * g_xxz_zz_x_yz[i] * a_exp;

        g_z_0_0_0_xx_zz_x_zz[i] = 2.0 * g_xxz_zz_x_zz[i] * a_exp;
    }
    // integrals block (1392-1398)

    #pragma omp simd aligned(g_xxz_zz_y_xx, g_xxz_zz_y_xy, g_xxz_zz_y_xz, g_xxz_zz_y_yy, g_xxz_zz_y_yz, g_xxz_zz_y_zz, g_z_0_0_0_xx_zz_y_xx, g_z_0_0_0_xx_zz_y_xy, g_z_0_0_0_xx_zz_y_xz, g_z_0_0_0_xx_zz_y_yy, g_z_0_0_0_xx_zz_y_yz, g_z_0_0_0_xx_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_zz_y_xx[i] = 2.0 * g_xxz_zz_y_xx[i] * a_exp;

        g_z_0_0_0_xx_zz_y_xy[i] = 2.0 * g_xxz_zz_y_xy[i] * a_exp;

        g_z_0_0_0_xx_zz_y_xz[i] = 2.0 * g_xxz_zz_y_xz[i] * a_exp;

        g_z_0_0_0_xx_zz_y_yy[i] = 2.0 * g_xxz_zz_y_yy[i] * a_exp;

        g_z_0_0_0_xx_zz_y_yz[i] = 2.0 * g_xxz_zz_y_yz[i] * a_exp;

        g_z_0_0_0_xx_zz_y_zz[i] = 2.0 * g_xxz_zz_y_zz[i] * a_exp;
    }
    // integrals block (1398-1404)

    #pragma omp simd aligned(g_xxz_zz_z_xx, g_xxz_zz_z_xy, g_xxz_zz_z_xz, g_xxz_zz_z_yy, g_xxz_zz_z_yz, g_xxz_zz_z_zz, g_z_0_0_0_xx_zz_z_xx, g_z_0_0_0_xx_zz_z_xy, g_z_0_0_0_xx_zz_z_xz, g_z_0_0_0_xx_zz_z_yy, g_z_0_0_0_xx_zz_z_yz, g_z_0_0_0_xx_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_zz_z_xx[i] = 2.0 * g_xxz_zz_z_xx[i] * a_exp;

        g_z_0_0_0_xx_zz_z_xy[i] = 2.0 * g_xxz_zz_z_xy[i] * a_exp;

        g_z_0_0_0_xx_zz_z_xz[i] = 2.0 * g_xxz_zz_z_xz[i] * a_exp;

        g_z_0_0_0_xx_zz_z_yy[i] = 2.0 * g_xxz_zz_z_yy[i] * a_exp;

        g_z_0_0_0_xx_zz_z_yz[i] = 2.0 * g_xxz_zz_z_yz[i] * a_exp;

        g_z_0_0_0_xx_zz_z_zz[i] = 2.0 * g_xxz_zz_z_zz[i] * a_exp;
    }
    // integrals block (1404-1410)

    #pragma omp simd aligned(g_xyz_xx_x_xx, g_xyz_xx_x_xy, g_xyz_xx_x_xz, g_xyz_xx_x_yy, g_xyz_xx_x_yz, g_xyz_xx_x_zz, g_z_0_0_0_xy_xx_x_xx, g_z_0_0_0_xy_xx_x_xy, g_z_0_0_0_xy_xx_x_xz, g_z_0_0_0_xy_xx_x_yy, g_z_0_0_0_xy_xx_x_yz, g_z_0_0_0_xy_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_xx_x_xx[i] = 2.0 * g_xyz_xx_x_xx[i] * a_exp;

        g_z_0_0_0_xy_xx_x_xy[i] = 2.0 * g_xyz_xx_x_xy[i] * a_exp;

        g_z_0_0_0_xy_xx_x_xz[i] = 2.0 * g_xyz_xx_x_xz[i] * a_exp;

        g_z_0_0_0_xy_xx_x_yy[i] = 2.0 * g_xyz_xx_x_yy[i] * a_exp;

        g_z_0_0_0_xy_xx_x_yz[i] = 2.0 * g_xyz_xx_x_yz[i] * a_exp;

        g_z_0_0_0_xy_xx_x_zz[i] = 2.0 * g_xyz_xx_x_zz[i] * a_exp;
    }
    // integrals block (1410-1416)

    #pragma omp simd aligned(g_xyz_xx_y_xx, g_xyz_xx_y_xy, g_xyz_xx_y_xz, g_xyz_xx_y_yy, g_xyz_xx_y_yz, g_xyz_xx_y_zz, g_z_0_0_0_xy_xx_y_xx, g_z_0_0_0_xy_xx_y_xy, g_z_0_0_0_xy_xx_y_xz, g_z_0_0_0_xy_xx_y_yy, g_z_0_0_0_xy_xx_y_yz, g_z_0_0_0_xy_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_xx_y_xx[i] = 2.0 * g_xyz_xx_y_xx[i] * a_exp;

        g_z_0_0_0_xy_xx_y_xy[i] = 2.0 * g_xyz_xx_y_xy[i] * a_exp;

        g_z_0_0_0_xy_xx_y_xz[i] = 2.0 * g_xyz_xx_y_xz[i] * a_exp;

        g_z_0_0_0_xy_xx_y_yy[i] = 2.0 * g_xyz_xx_y_yy[i] * a_exp;

        g_z_0_0_0_xy_xx_y_yz[i] = 2.0 * g_xyz_xx_y_yz[i] * a_exp;

        g_z_0_0_0_xy_xx_y_zz[i] = 2.0 * g_xyz_xx_y_zz[i] * a_exp;
    }
    // integrals block (1416-1422)

    #pragma omp simd aligned(g_xyz_xx_z_xx, g_xyz_xx_z_xy, g_xyz_xx_z_xz, g_xyz_xx_z_yy, g_xyz_xx_z_yz, g_xyz_xx_z_zz, g_z_0_0_0_xy_xx_z_xx, g_z_0_0_0_xy_xx_z_xy, g_z_0_0_0_xy_xx_z_xz, g_z_0_0_0_xy_xx_z_yy, g_z_0_0_0_xy_xx_z_yz, g_z_0_0_0_xy_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_xx_z_xx[i] = 2.0 * g_xyz_xx_z_xx[i] * a_exp;

        g_z_0_0_0_xy_xx_z_xy[i] = 2.0 * g_xyz_xx_z_xy[i] * a_exp;

        g_z_0_0_0_xy_xx_z_xz[i] = 2.0 * g_xyz_xx_z_xz[i] * a_exp;

        g_z_0_0_0_xy_xx_z_yy[i] = 2.0 * g_xyz_xx_z_yy[i] * a_exp;

        g_z_0_0_0_xy_xx_z_yz[i] = 2.0 * g_xyz_xx_z_yz[i] * a_exp;

        g_z_0_0_0_xy_xx_z_zz[i] = 2.0 * g_xyz_xx_z_zz[i] * a_exp;
    }
    // integrals block (1422-1428)

    #pragma omp simd aligned(g_xyz_xy_x_xx, g_xyz_xy_x_xy, g_xyz_xy_x_xz, g_xyz_xy_x_yy, g_xyz_xy_x_yz, g_xyz_xy_x_zz, g_z_0_0_0_xy_xy_x_xx, g_z_0_0_0_xy_xy_x_xy, g_z_0_0_0_xy_xy_x_xz, g_z_0_0_0_xy_xy_x_yy, g_z_0_0_0_xy_xy_x_yz, g_z_0_0_0_xy_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_xy_x_xx[i] = 2.0 * g_xyz_xy_x_xx[i] * a_exp;

        g_z_0_0_0_xy_xy_x_xy[i] = 2.0 * g_xyz_xy_x_xy[i] * a_exp;

        g_z_0_0_0_xy_xy_x_xz[i] = 2.0 * g_xyz_xy_x_xz[i] * a_exp;

        g_z_0_0_0_xy_xy_x_yy[i] = 2.0 * g_xyz_xy_x_yy[i] * a_exp;

        g_z_0_0_0_xy_xy_x_yz[i] = 2.0 * g_xyz_xy_x_yz[i] * a_exp;

        g_z_0_0_0_xy_xy_x_zz[i] = 2.0 * g_xyz_xy_x_zz[i] * a_exp;
    }
    // integrals block (1428-1434)

    #pragma omp simd aligned(g_xyz_xy_y_xx, g_xyz_xy_y_xy, g_xyz_xy_y_xz, g_xyz_xy_y_yy, g_xyz_xy_y_yz, g_xyz_xy_y_zz, g_z_0_0_0_xy_xy_y_xx, g_z_0_0_0_xy_xy_y_xy, g_z_0_0_0_xy_xy_y_xz, g_z_0_0_0_xy_xy_y_yy, g_z_0_0_0_xy_xy_y_yz, g_z_0_0_0_xy_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_xy_y_xx[i] = 2.0 * g_xyz_xy_y_xx[i] * a_exp;

        g_z_0_0_0_xy_xy_y_xy[i] = 2.0 * g_xyz_xy_y_xy[i] * a_exp;

        g_z_0_0_0_xy_xy_y_xz[i] = 2.0 * g_xyz_xy_y_xz[i] * a_exp;

        g_z_0_0_0_xy_xy_y_yy[i] = 2.0 * g_xyz_xy_y_yy[i] * a_exp;

        g_z_0_0_0_xy_xy_y_yz[i] = 2.0 * g_xyz_xy_y_yz[i] * a_exp;

        g_z_0_0_0_xy_xy_y_zz[i] = 2.0 * g_xyz_xy_y_zz[i] * a_exp;
    }
    // integrals block (1434-1440)

    #pragma omp simd aligned(g_xyz_xy_z_xx, g_xyz_xy_z_xy, g_xyz_xy_z_xz, g_xyz_xy_z_yy, g_xyz_xy_z_yz, g_xyz_xy_z_zz, g_z_0_0_0_xy_xy_z_xx, g_z_0_0_0_xy_xy_z_xy, g_z_0_0_0_xy_xy_z_xz, g_z_0_0_0_xy_xy_z_yy, g_z_0_0_0_xy_xy_z_yz, g_z_0_0_0_xy_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_xy_z_xx[i] = 2.0 * g_xyz_xy_z_xx[i] * a_exp;

        g_z_0_0_0_xy_xy_z_xy[i] = 2.0 * g_xyz_xy_z_xy[i] * a_exp;

        g_z_0_0_0_xy_xy_z_xz[i] = 2.0 * g_xyz_xy_z_xz[i] * a_exp;

        g_z_0_0_0_xy_xy_z_yy[i] = 2.0 * g_xyz_xy_z_yy[i] * a_exp;

        g_z_0_0_0_xy_xy_z_yz[i] = 2.0 * g_xyz_xy_z_yz[i] * a_exp;

        g_z_0_0_0_xy_xy_z_zz[i] = 2.0 * g_xyz_xy_z_zz[i] * a_exp;
    }
    // integrals block (1440-1446)

    #pragma omp simd aligned(g_xyz_xz_x_xx, g_xyz_xz_x_xy, g_xyz_xz_x_xz, g_xyz_xz_x_yy, g_xyz_xz_x_yz, g_xyz_xz_x_zz, g_z_0_0_0_xy_xz_x_xx, g_z_0_0_0_xy_xz_x_xy, g_z_0_0_0_xy_xz_x_xz, g_z_0_0_0_xy_xz_x_yy, g_z_0_0_0_xy_xz_x_yz, g_z_0_0_0_xy_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_xz_x_xx[i] = 2.0 * g_xyz_xz_x_xx[i] * a_exp;

        g_z_0_0_0_xy_xz_x_xy[i] = 2.0 * g_xyz_xz_x_xy[i] * a_exp;

        g_z_0_0_0_xy_xz_x_xz[i] = 2.0 * g_xyz_xz_x_xz[i] * a_exp;

        g_z_0_0_0_xy_xz_x_yy[i] = 2.0 * g_xyz_xz_x_yy[i] * a_exp;

        g_z_0_0_0_xy_xz_x_yz[i] = 2.0 * g_xyz_xz_x_yz[i] * a_exp;

        g_z_0_0_0_xy_xz_x_zz[i] = 2.0 * g_xyz_xz_x_zz[i] * a_exp;
    }
    // integrals block (1446-1452)

    #pragma omp simd aligned(g_xyz_xz_y_xx, g_xyz_xz_y_xy, g_xyz_xz_y_xz, g_xyz_xz_y_yy, g_xyz_xz_y_yz, g_xyz_xz_y_zz, g_z_0_0_0_xy_xz_y_xx, g_z_0_0_0_xy_xz_y_xy, g_z_0_0_0_xy_xz_y_xz, g_z_0_0_0_xy_xz_y_yy, g_z_0_0_0_xy_xz_y_yz, g_z_0_0_0_xy_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_xz_y_xx[i] = 2.0 * g_xyz_xz_y_xx[i] * a_exp;

        g_z_0_0_0_xy_xz_y_xy[i] = 2.0 * g_xyz_xz_y_xy[i] * a_exp;

        g_z_0_0_0_xy_xz_y_xz[i] = 2.0 * g_xyz_xz_y_xz[i] * a_exp;

        g_z_0_0_0_xy_xz_y_yy[i] = 2.0 * g_xyz_xz_y_yy[i] * a_exp;

        g_z_0_0_0_xy_xz_y_yz[i] = 2.0 * g_xyz_xz_y_yz[i] * a_exp;

        g_z_0_0_0_xy_xz_y_zz[i] = 2.0 * g_xyz_xz_y_zz[i] * a_exp;
    }
    // integrals block (1452-1458)

    #pragma omp simd aligned(g_xyz_xz_z_xx, g_xyz_xz_z_xy, g_xyz_xz_z_xz, g_xyz_xz_z_yy, g_xyz_xz_z_yz, g_xyz_xz_z_zz, g_z_0_0_0_xy_xz_z_xx, g_z_0_0_0_xy_xz_z_xy, g_z_0_0_0_xy_xz_z_xz, g_z_0_0_0_xy_xz_z_yy, g_z_0_0_0_xy_xz_z_yz, g_z_0_0_0_xy_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_xz_z_xx[i] = 2.0 * g_xyz_xz_z_xx[i] * a_exp;

        g_z_0_0_0_xy_xz_z_xy[i] = 2.0 * g_xyz_xz_z_xy[i] * a_exp;

        g_z_0_0_0_xy_xz_z_xz[i] = 2.0 * g_xyz_xz_z_xz[i] * a_exp;

        g_z_0_0_0_xy_xz_z_yy[i] = 2.0 * g_xyz_xz_z_yy[i] * a_exp;

        g_z_0_0_0_xy_xz_z_yz[i] = 2.0 * g_xyz_xz_z_yz[i] * a_exp;

        g_z_0_0_0_xy_xz_z_zz[i] = 2.0 * g_xyz_xz_z_zz[i] * a_exp;
    }
    // integrals block (1458-1464)

    #pragma omp simd aligned(g_xyz_yy_x_xx, g_xyz_yy_x_xy, g_xyz_yy_x_xz, g_xyz_yy_x_yy, g_xyz_yy_x_yz, g_xyz_yy_x_zz, g_z_0_0_0_xy_yy_x_xx, g_z_0_0_0_xy_yy_x_xy, g_z_0_0_0_xy_yy_x_xz, g_z_0_0_0_xy_yy_x_yy, g_z_0_0_0_xy_yy_x_yz, g_z_0_0_0_xy_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_yy_x_xx[i] = 2.0 * g_xyz_yy_x_xx[i] * a_exp;

        g_z_0_0_0_xy_yy_x_xy[i] = 2.0 * g_xyz_yy_x_xy[i] * a_exp;

        g_z_0_0_0_xy_yy_x_xz[i] = 2.0 * g_xyz_yy_x_xz[i] * a_exp;

        g_z_0_0_0_xy_yy_x_yy[i] = 2.0 * g_xyz_yy_x_yy[i] * a_exp;

        g_z_0_0_0_xy_yy_x_yz[i] = 2.0 * g_xyz_yy_x_yz[i] * a_exp;

        g_z_0_0_0_xy_yy_x_zz[i] = 2.0 * g_xyz_yy_x_zz[i] * a_exp;
    }
    // integrals block (1464-1470)

    #pragma omp simd aligned(g_xyz_yy_y_xx, g_xyz_yy_y_xy, g_xyz_yy_y_xz, g_xyz_yy_y_yy, g_xyz_yy_y_yz, g_xyz_yy_y_zz, g_z_0_0_0_xy_yy_y_xx, g_z_0_0_0_xy_yy_y_xy, g_z_0_0_0_xy_yy_y_xz, g_z_0_0_0_xy_yy_y_yy, g_z_0_0_0_xy_yy_y_yz, g_z_0_0_0_xy_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_yy_y_xx[i] = 2.0 * g_xyz_yy_y_xx[i] * a_exp;

        g_z_0_0_0_xy_yy_y_xy[i] = 2.0 * g_xyz_yy_y_xy[i] * a_exp;

        g_z_0_0_0_xy_yy_y_xz[i] = 2.0 * g_xyz_yy_y_xz[i] * a_exp;

        g_z_0_0_0_xy_yy_y_yy[i] = 2.0 * g_xyz_yy_y_yy[i] * a_exp;

        g_z_0_0_0_xy_yy_y_yz[i] = 2.0 * g_xyz_yy_y_yz[i] * a_exp;

        g_z_0_0_0_xy_yy_y_zz[i] = 2.0 * g_xyz_yy_y_zz[i] * a_exp;
    }
    // integrals block (1470-1476)

    #pragma omp simd aligned(g_xyz_yy_z_xx, g_xyz_yy_z_xy, g_xyz_yy_z_xz, g_xyz_yy_z_yy, g_xyz_yy_z_yz, g_xyz_yy_z_zz, g_z_0_0_0_xy_yy_z_xx, g_z_0_0_0_xy_yy_z_xy, g_z_0_0_0_xy_yy_z_xz, g_z_0_0_0_xy_yy_z_yy, g_z_0_0_0_xy_yy_z_yz, g_z_0_0_0_xy_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_yy_z_xx[i] = 2.0 * g_xyz_yy_z_xx[i] * a_exp;

        g_z_0_0_0_xy_yy_z_xy[i] = 2.0 * g_xyz_yy_z_xy[i] * a_exp;

        g_z_0_0_0_xy_yy_z_xz[i] = 2.0 * g_xyz_yy_z_xz[i] * a_exp;

        g_z_0_0_0_xy_yy_z_yy[i] = 2.0 * g_xyz_yy_z_yy[i] * a_exp;

        g_z_0_0_0_xy_yy_z_yz[i] = 2.0 * g_xyz_yy_z_yz[i] * a_exp;

        g_z_0_0_0_xy_yy_z_zz[i] = 2.0 * g_xyz_yy_z_zz[i] * a_exp;
    }
    // integrals block (1476-1482)

    #pragma omp simd aligned(g_xyz_yz_x_xx, g_xyz_yz_x_xy, g_xyz_yz_x_xz, g_xyz_yz_x_yy, g_xyz_yz_x_yz, g_xyz_yz_x_zz, g_z_0_0_0_xy_yz_x_xx, g_z_0_0_0_xy_yz_x_xy, g_z_0_0_0_xy_yz_x_xz, g_z_0_0_0_xy_yz_x_yy, g_z_0_0_0_xy_yz_x_yz, g_z_0_0_0_xy_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_yz_x_xx[i] = 2.0 * g_xyz_yz_x_xx[i] * a_exp;

        g_z_0_0_0_xy_yz_x_xy[i] = 2.0 * g_xyz_yz_x_xy[i] * a_exp;

        g_z_0_0_0_xy_yz_x_xz[i] = 2.0 * g_xyz_yz_x_xz[i] * a_exp;

        g_z_0_0_0_xy_yz_x_yy[i] = 2.0 * g_xyz_yz_x_yy[i] * a_exp;

        g_z_0_0_0_xy_yz_x_yz[i] = 2.0 * g_xyz_yz_x_yz[i] * a_exp;

        g_z_0_0_0_xy_yz_x_zz[i] = 2.0 * g_xyz_yz_x_zz[i] * a_exp;
    }
    // integrals block (1482-1488)

    #pragma omp simd aligned(g_xyz_yz_y_xx, g_xyz_yz_y_xy, g_xyz_yz_y_xz, g_xyz_yz_y_yy, g_xyz_yz_y_yz, g_xyz_yz_y_zz, g_z_0_0_0_xy_yz_y_xx, g_z_0_0_0_xy_yz_y_xy, g_z_0_0_0_xy_yz_y_xz, g_z_0_0_0_xy_yz_y_yy, g_z_0_0_0_xy_yz_y_yz, g_z_0_0_0_xy_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_yz_y_xx[i] = 2.0 * g_xyz_yz_y_xx[i] * a_exp;

        g_z_0_0_0_xy_yz_y_xy[i] = 2.0 * g_xyz_yz_y_xy[i] * a_exp;

        g_z_0_0_0_xy_yz_y_xz[i] = 2.0 * g_xyz_yz_y_xz[i] * a_exp;

        g_z_0_0_0_xy_yz_y_yy[i] = 2.0 * g_xyz_yz_y_yy[i] * a_exp;

        g_z_0_0_0_xy_yz_y_yz[i] = 2.0 * g_xyz_yz_y_yz[i] * a_exp;

        g_z_0_0_0_xy_yz_y_zz[i] = 2.0 * g_xyz_yz_y_zz[i] * a_exp;
    }
    // integrals block (1488-1494)

    #pragma omp simd aligned(g_xyz_yz_z_xx, g_xyz_yz_z_xy, g_xyz_yz_z_xz, g_xyz_yz_z_yy, g_xyz_yz_z_yz, g_xyz_yz_z_zz, g_z_0_0_0_xy_yz_z_xx, g_z_0_0_0_xy_yz_z_xy, g_z_0_0_0_xy_yz_z_xz, g_z_0_0_0_xy_yz_z_yy, g_z_0_0_0_xy_yz_z_yz, g_z_0_0_0_xy_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_yz_z_xx[i] = 2.0 * g_xyz_yz_z_xx[i] * a_exp;

        g_z_0_0_0_xy_yz_z_xy[i] = 2.0 * g_xyz_yz_z_xy[i] * a_exp;

        g_z_0_0_0_xy_yz_z_xz[i] = 2.0 * g_xyz_yz_z_xz[i] * a_exp;

        g_z_0_0_0_xy_yz_z_yy[i] = 2.0 * g_xyz_yz_z_yy[i] * a_exp;

        g_z_0_0_0_xy_yz_z_yz[i] = 2.0 * g_xyz_yz_z_yz[i] * a_exp;

        g_z_0_0_0_xy_yz_z_zz[i] = 2.0 * g_xyz_yz_z_zz[i] * a_exp;
    }
    // integrals block (1494-1500)

    #pragma omp simd aligned(g_xyz_zz_x_xx, g_xyz_zz_x_xy, g_xyz_zz_x_xz, g_xyz_zz_x_yy, g_xyz_zz_x_yz, g_xyz_zz_x_zz, g_z_0_0_0_xy_zz_x_xx, g_z_0_0_0_xy_zz_x_xy, g_z_0_0_0_xy_zz_x_xz, g_z_0_0_0_xy_zz_x_yy, g_z_0_0_0_xy_zz_x_yz, g_z_0_0_0_xy_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_zz_x_xx[i] = 2.0 * g_xyz_zz_x_xx[i] * a_exp;

        g_z_0_0_0_xy_zz_x_xy[i] = 2.0 * g_xyz_zz_x_xy[i] * a_exp;

        g_z_0_0_0_xy_zz_x_xz[i] = 2.0 * g_xyz_zz_x_xz[i] * a_exp;

        g_z_0_0_0_xy_zz_x_yy[i] = 2.0 * g_xyz_zz_x_yy[i] * a_exp;

        g_z_0_0_0_xy_zz_x_yz[i] = 2.0 * g_xyz_zz_x_yz[i] * a_exp;

        g_z_0_0_0_xy_zz_x_zz[i] = 2.0 * g_xyz_zz_x_zz[i] * a_exp;
    }
    // integrals block (1500-1506)

    #pragma omp simd aligned(g_xyz_zz_y_xx, g_xyz_zz_y_xy, g_xyz_zz_y_xz, g_xyz_zz_y_yy, g_xyz_zz_y_yz, g_xyz_zz_y_zz, g_z_0_0_0_xy_zz_y_xx, g_z_0_0_0_xy_zz_y_xy, g_z_0_0_0_xy_zz_y_xz, g_z_0_0_0_xy_zz_y_yy, g_z_0_0_0_xy_zz_y_yz, g_z_0_0_0_xy_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_zz_y_xx[i] = 2.0 * g_xyz_zz_y_xx[i] * a_exp;

        g_z_0_0_0_xy_zz_y_xy[i] = 2.0 * g_xyz_zz_y_xy[i] * a_exp;

        g_z_0_0_0_xy_zz_y_xz[i] = 2.0 * g_xyz_zz_y_xz[i] * a_exp;

        g_z_0_0_0_xy_zz_y_yy[i] = 2.0 * g_xyz_zz_y_yy[i] * a_exp;

        g_z_0_0_0_xy_zz_y_yz[i] = 2.0 * g_xyz_zz_y_yz[i] * a_exp;

        g_z_0_0_0_xy_zz_y_zz[i] = 2.0 * g_xyz_zz_y_zz[i] * a_exp;
    }
    // integrals block (1506-1512)

    #pragma omp simd aligned(g_xyz_zz_z_xx, g_xyz_zz_z_xy, g_xyz_zz_z_xz, g_xyz_zz_z_yy, g_xyz_zz_z_yz, g_xyz_zz_z_zz, g_z_0_0_0_xy_zz_z_xx, g_z_0_0_0_xy_zz_z_xy, g_z_0_0_0_xy_zz_z_xz, g_z_0_0_0_xy_zz_z_yy, g_z_0_0_0_xy_zz_z_yz, g_z_0_0_0_xy_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_zz_z_xx[i] = 2.0 * g_xyz_zz_z_xx[i] * a_exp;

        g_z_0_0_0_xy_zz_z_xy[i] = 2.0 * g_xyz_zz_z_xy[i] * a_exp;

        g_z_0_0_0_xy_zz_z_xz[i] = 2.0 * g_xyz_zz_z_xz[i] * a_exp;

        g_z_0_0_0_xy_zz_z_yy[i] = 2.0 * g_xyz_zz_z_yy[i] * a_exp;

        g_z_0_0_0_xy_zz_z_yz[i] = 2.0 * g_xyz_zz_z_yz[i] * a_exp;

        g_z_0_0_0_xy_zz_z_zz[i] = 2.0 * g_xyz_zz_z_zz[i] * a_exp;
    }
    // integrals block (1512-1518)

    #pragma omp simd aligned(g_x_xx_x_xx, g_x_xx_x_xy, g_x_xx_x_xz, g_x_xx_x_yy, g_x_xx_x_yz, g_x_xx_x_zz, g_xzz_xx_x_xx, g_xzz_xx_x_xy, g_xzz_xx_x_xz, g_xzz_xx_x_yy, g_xzz_xx_x_yz, g_xzz_xx_x_zz, g_z_0_0_0_xz_xx_x_xx, g_z_0_0_0_xz_xx_x_xy, g_z_0_0_0_xz_xx_x_xz, g_z_0_0_0_xz_xx_x_yy, g_z_0_0_0_xz_xx_x_yz, g_z_0_0_0_xz_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_xx_x_xx[i] = -g_x_xx_x_xx[i] + 2.0 * g_xzz_xx_x_xx[i] * a_exp;

        g_z_0_0_0_xz_xx_x_xy[i] = -g_x_xx_x_xy[i] + 2.0 * g_xzz_xx_x_xy[i] * a_exp;

        g_z_0_0_0_xz_xx_x_xz[i] = -g_x_xx_x_xz[i] + 2.0 * g_xzz_xx_x_xz[i] * a_exp;

        g_z_0_0_0_xz_xx_x_yy[i] = -g_x_xx_x_yy[i] + 2.0 * g_xzz_xx_x_yy[i] * a_exp;

        g_z_0_0_0_xz_xx_x_yz[i] = -g_x_xx_x_yz[i] + 2.0 * g_xzz_xx_x_yz[i] * a_exp;

        g_z_0_0_0_xz_xx_x_zz[i] = -g_x_xx_x_zz[i] + 2.0 * g_xzz_xx_x_zz[i] * a_exp;
    }
    // integrals block (1518-1524)

    #pragma omp simd aligned(g_x_xx_y_xx, g_x_xx_y_xy, g_x_xx_y_xz, g_x_xx_y_yy, g_x_xx_y_yz, g_x_xx_y_zz, g_xzz_xx_y_xx, g_xzz_xx_y_xy, g_xzz_xx_y_xz, g_xzz_xx_y_yy, g_xzz_xx_y_yz, g_xzz_xx_y_zz, g_z_0_0_0_xz_xx_y_xx, g_z_0_0_0_xz_xx_y_xy, g_z_0_0_0_xz_xx_y_xz, g_z_0_0_0_xz_xx_y_yy, g_z_0_0_0_xz_xx_y_yz, g_z_0_0_0_xz_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_xx_y_xx[i] = -g_x_xx_y_xx[i] + 2.0 * g_xzz_xx_y_xx[i] * a_exp;

        g_z_0_0_0_xz_xx_y_xy[i] = -g_x_xx_y_xy[i] + 2.0 * g_xzz_xx_y_xy[i] * a_exp;

        g_z_0_0_0_xz_xx_y_xz[i] = -g_x_xx_y_xz[i] + 2.0 * g_xzz_xx_y_xz[i] * a_exp;

        g_z_0_0_0_xz_xx_y_yy[i] = -g_x_xx_y_yy[i] + 2.0 * g_xzz_xx_y_yy[i] * a_exp;

        g_z_0_0_0_xz_xx_y_yz[i] = -g_x_xx_y_yz[i] + 2.0 * g_xzz_xx_y_yz[i] * a_exp;

        g_z_0_0_0_xz_xx_y_zz[i] = -g_x_xx_y_zz[i] + 2.0 * g_xzz_xx_y_zz[i] * a_exp;
    }
    // integrals block (1524-1530)

    #pragma omp simd aligned(g_x_xx_z_xx, g_x_xx_z_xy, g_x_xx_z_xz, g_x_xx_z_yy, g_x_xx_z_yz, g_x_xx_z_zz, g_xzz_xx_z_xx, g_xzz_xx_z_xy, g_xzz_xx_z_xz, g_xzz_xx_z_yy, g_xzz_xx_z_yz, g_xzz_xx_z_zz, g_z_0_0_0_xz_xx_z_xx, g_z_0_0_0_xz_xx_z_xy, g_z_0_0_0_xz_xx_z_xz, g_z_0_0_0_xz_xx_z_yy, g_z_0_0_0_xz_xx_z_yz, g_z_0_0_0_xz_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_xx_z_xx[i] = -g_x_xx_z_xx[i] + 2.0 * g_xzz_xx_z_xx[i] * a_exp;

        g_z_0_0_0_xz_xx_z_xy[i] = -g_x_xx_z_xy[i] + 2.0 * g_xzz_xx_z_xy[i] * a_exp;

        g_z_0_0_0_xz_xx_z_xz[i] = -g_x_xx_z_xz[i] + 2.0 * g_xzz_xx_z_xz[i] * a_exp;

        g_z_0_0_0_xz_xx_z_yy[i] = -g_x_xx_z_yy[i] + 2.0 * g_xzz_xx_z_yy[i] * a_exp;

        g_z_0_0_0_xz_xx_z_yz[i] = -g_x_xx_z_yz[i] + 2.0 * g_xzz_xx_z_yz[i] * a_exp;

        g_z_0_0_0_xz_xx_z_zz[i] = -g_x_xx_z_zz[i] + 2.0 * g_xzz_xx_z_zz[i] * a_exp;
    }
    // integrals block (1530-1536)

    #pragma omp simd aligned(g_x_xy_x_xx, g_x_xy_x_xy, g_x_xy_x_xz, g_x_xy_x_yy, g_x_xy_x_yz, g_x_xy_x_zz, g_xzz_xy_x_xx, g_xzz_xy_x_xy, g_xzz_xy_x_xz, g_xzz_xy_x_yy, g_xzz_xy_x_yz, g_xzz_xy_x_zz, g_z_0_0_0_xz_xy_x_xx, g_z_0_0_0_xz_xy_x_xy, g_z_0_0_0_xz_xy_x_xz, g_z_0_0_0_xz_xy_x_yy, g_z_0_0_0_xz_xy_x_yz, g_z_0_0_0_xz_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_xy_x_xx[i] = -g_x_xy_x_xx[i] + 2.0 * g_xzz_xy_x_xx[i] * a_exp;

        g_z_0_0_0_xz_xy_x_xy[i] = -g_x_xy_x_xy[i] + 2.0 * g_xzz_xy_x_xy[i] * a_exp;

        g_z_0_0_0_xz_xy_x_xz[i] = -g_x_xy_x_xz[i] + 2.0 * g_xzz_xy_x_xz[i] * a_exp;

        g_z_0_0_0_xz_xy_x_yy[i] = -g_x_xy_x_yy[i] + 2.0 * g_xzz_xy_x_yy[i] * a_exp;

        g_z_0_0_0_xz_xy_x_yz[i] = -g_x_xy_x_yz[i] + 2.0 * g_xzz_xy_x_yz[i] * a_exp;

        g_z_0_0_0_xz_xy_x_zz[i] = -g_x_xy_x_zz[i] + 2.0 * g_xzz_xy_x_zz[i] * a_exp;
    }
    // integrals block (1536-1542)

    #pragma omp simd aligned(g_x_xy_y_xx, g_x_xy_y_xy, g_x_xy_y_xz, g_x_xy_y_yy, g_x_xy_y_yz, g_x_xy_y_zz, g_xzz_xy_y_xx, g_xzz_xy_y_xy, g_xzz_xy_y_xz, g_xzz_xy_y_yy, g_xzz_xy_y_yz, g_xzz_xy_y_zz, g_z_0_0_0_xz_xy_y_xx, g_z_0_0_0_xz_xy_y_xy, g_z_0_0_0_xz_xy_y_xz, g_z_0_0_0_xz_xy_y_yy, g_z_0_0_0_xz_xy_y_yz, g_z_0_0_0_xz_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_xy_y_xx[i] = -g_x_xy_y_xx[i] + 2.0 * g_xzz_xy_y_xx[i] * a_exp;

        g_z_0_0_0_xz_xy_y_xy[i] = -g_x_xy_y_xy[i] + 2.0 * g_xzz_xy_y_xy[i] * a_exp;

        g_z_0_0_0_xz_xy_y_xz[i] = -g_x_xy_y_xz[i] + 2.0 * g_xzz_xy_y_xz[i] * a_exp;

        g_z_0_0_0_xz_xy_y_yy[i] = -g_x_xy_y_yy[i] + 2.0 * g_xzz_xy_y_yy[i] * a_exp;

        g_z_0_0_0_xz_xy_y_yz[i] = -g_x_xy_y_yz[i] + 2.0 * g_xzz_xy_y_yz[i] * a_exp;

        g_z_0_0_0_xz_xy_y_zz[i] = -g_x_xy_y_zz[i] + 2.0 * g_xzz_xy_y_zz[i] * a_exp;
    }
    // integrals block (1542-1548)

    #pragma omp simd aligned(g_x_xy_z_xx, g_x_xy_z_xy, g_x_xy_z_xz, g_x_xy_z_yy, g_x_xy_z_yz, g_x_xy_z_zz, g_xzz_xy_z_xx, g_xzz_xy_z_xy, g_xzz_xy_z_xz, g_xzz_xy_z_yy, g_xzz_xy_z_yz, g_xzz_xy_z_zz, g_z_0_0_0_xz_xy_z_xx, g_z_0_0_0_xz_xy_z_xy, g_z_0_0_0_xz_xy_z_xz, g_z_0_0_0_xz_xy_z_yy, g_z_0_0_0_xz_xy_z_yz, g_z_0_0_0_xz_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_xy_z_xx[i] = -g_x_xy_z_xx[i] + 2.0 * g_xzz_xy_z_xx[i] * a_exp;

        g_z_0_0_0_xz_xy_z_xy[i] = -g_x_xy_z_xy[i] + 2.0 * g_xzz_xy_z_xy[i] * a_exp;

        g_z_0_0_0_xz_xy_z_xz[i] = -g_x_xy_z_xz[i] + 2.0 * g_xzz_xy_z_xz[i] * a_exp;

        g_z_0_0_0_xz_xy_z_yy[i] = -g_x_xy_z_yy[i] + 2.0 * g_xzz_xy_z_yy[i] * a_exp;

        g_z_0_0_0_xz_xy_z_yz[i] = -g_x_xy_z_yz[i] + 2.0 * g_xzz_xy_z_yz[i] * a_exp;

        g_z_0_0_0_xz_xy_z_zz[i] = -g_x_xy_z_zz[i] + 2.0 * g_xzz_xy_z_zz[i] * a_exp;
    }
    // integrals block (1548-1554)

    #pragma omp simd aligned(g_x_xz_x_xx, g_x_xz_x_xy, g_x_xz_x_xz, g_x_xz_x_yy, g_x_xz_x_yz, g_x_xz_x_zz, g_xzz_xz_x_xx, g_xzz_xz_x_xy, g_xzz_xz_x_xz, g_xzz_xz_x_yy, g_xzz_xz_x_yz, g_xzz_xz_x_zz, g_z_0_0_0_xz_xz_x_xx, g_z_0_0_0_xz_xz_x_xy, g_z_0_0_0_xz_xz_x_xz, g_z_0_0_0_xz_xz_x_yy, g_z_0_0_0_xz_xz_x_yz, g_z_0_0_0_xz_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_xz_x_xx[i] = -g_x_xz_x_xx[i] + 2.0 * g_xzz_xz_x_xx[i] * a_exp;

        g_z_0_0_0_xz_xz_x_xy[i] = -g_x_xz_x_xy[i] + 2.0 * g_xzz_xz_x_xy[i] * a_exp;

        g_z_0_0_0_xz_xz_x_xz[i] = -g_x_xz_x_xz[i] + 2.0 * g_xzz_xz_x_xz[i] * a_exp;

        g_z_0_0_0_xz_xz_x_yy[i] = -g_x_xz_x_yy[i] + 2.0 * g_xzz_xz_x_yy[i] * a_exp;

        g_z_0_0_0_xz_xz_x_yz[i] = -g_x_xz_x_yz[i] + 2.0 * g_xzz_xz_x_yz[i] * a_exp;

        g_z_0_0_0_xz_xz_x_zz[i] = -g_x_xz_x_zz[i] + 2.0 * g_xzz_xz_x_zz[i] * a_exp;
    }
    // integrals block (1554-1560)

    #pragma omp simd aligned(g_x_xz_y_xx, g_x_xz_y_xy, g_x_xz_y_xz, g_x_xz_y_yy, g_x_xz_y_yz, g_x_xz_y_zz, g_xzz_xz_y_xx, g_xzz_xz_y_xy, g_xzz_xz_y_xz, g_xzz_xz_y_yy, g_xzz_xz_y_yz, g_xzz_xz_y_zz, g_z_0_0_0_xz_xz_y_xx, g_z_0_0_0_xz_xz_y_xy, g_z_0_0_0_xz_xz_y_xz, g_z_0_0_0_xz_xz_y_yy, g_z_0_0_0_xz_xz_y_yz, g_z_0_0_0_xz_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_xz_y_xx[i] = -g_x_xz_y_xx[i] + 2.0 * g_xzz_xz_y_xx[i] * a_exp;

        g_z_0_0_0_xz_xz_y_xy[i] = -g_x_xz_y_xy[i] + 2.0 * g_xzz_xz_y_xy[i] * a_exp;

        g_z_0_0_0_xz_xz_y_xz[i] = -g_x_xz_y_xz[i] + 2.0 * g_xzz_xz_y_xz[i] * a_exp;

        g_z_0_0_0_xz_xz_y_yy[i] = -g_x_xz_y_yy[i] + 2.0 * g_xzz_xz_y_yy[i] * a_exp;

        g_z_0_0_0_xz_xz_y_yz[i] = -g_x_xz_y_yz[i] + 2.0 * g_xzz_xz_y_yz[i] * a_exp;

        g_z_0_0_0_xz_xz_y_zz[i] = -g_x_xz_y_zz[i] + 2.0 * g_xzz_xz_y_zz[i] * a_exp;
    }
    // integrals block (1560-1566)

    #pragma omp simd aligned(g_x_xz_z_xx, g_x_xz_z_xy, g_x_xz_z_xz, g_x_xz_z_yy, g_x_xz_z_yz, g_x_xz_z_zz, g_xzz_xz_z_xx, g_xzz_xz_z_xy, g_xzz_xz_z_xz, g_xzz_xz_z_yy, g_xzz_xz_z_yz, g_xzz_xz_z_zz, g_z_0_0_0_xz_xz_z_xx, g_z_0_0_0_xz_xz_z_xy, g_z_0_0_0_xz_xz_z_xz, g_z_0_0_0_xz_xz_z_yy, g_z_0_0_0_xz_xz_z_yz, g_z_0_0_0_xz_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_xz_z_xx[i] = -g_x_xz_z_xx[i] + 2.0 * g_xzz_xz_z_xx[i] * a_exp;

        g_z_0_0_0_xz_xz_z_xy[i] = -g_x_xz_z_xy[i] + 2.0 * g_xzz_xz_z_xy[i] * a_exp;

        g_z_0_0_0_xz_xz_z_xz[i] = -g_x_xz_z_xz[i] + 2.0 * g_xzz_xz_z_xz[i] * a_exp;

        g_z_0_0_0_xz_xz_z_yy[i] = -g_x_xz_z_yy[i] + 2.0 * g_xzz_xz_z_yy[i] * a_exp;

        g_z_0_0_0_xz_xz_z_yz[i] = -g_x_xz_z_yz[i] + 2.0 * g_xzz_xz_z_yz[i] * a_exp;

        g_z_0_0_0_xz_xz_z_zz[i] = -g_x_xz_z_zz[i] + 2.0 * g_xzz_xz_z_zz[i] * a_exp;
    }
    // integrals block (1566-1572)

    #pragma omp simd aligned(g_x_yy_x_xx, g_x_yy_x_xy, g_x_yy_x_xz, g_x_yy_x_yy, g_x_yy_x_yz, g_x_yy_x_zz, g_xzz_yy_x_xx, g_xzz_yy_x_xy, g_xzz_yy_x_xz, g_xzz_yy_x_yy, g_xzz_yy_x_yz, g_xzz_yy_x_zz, g_z_0_0_0_xz_yy_x_xx, g_z_0_0_0_xz_yy_x_xy, g_z_0_0_0_xz_yy_x_xz, g_z_0_0_0_xz_yy_x_yy, g_z_0_0_0_xz_yy_x_yz, g_z_0_0_0_xz_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_yy_x_xx[i] = -g_x_yy_x_xx[i] + 2.0 * g_xzz_yy_x_xx[i] * a_exp;

        g_z_0_0_0_xz_yy_x_xy[i] = -g_x_yy_x_xy[i] + 2.0 * g_xzz_yy_x_xy[i] * a_exp;

        g_z_0_0_0_xz_yy_x_xz[i] = -g_x_yy_x_xz[i] + 2.0 * g_xzz_yy_x_xz[i] * a_exp;

        g_z_0_0_0_xz_yy_x_yy[i] = -g_x_yy_x_yy[i] + 2.0 * g_xzz_yy_x_yy[i] * a_exp;

        g_z_0_0_0_xz_yy_x_yz[i] = -g_x_yy_x_yz[i] + 2.0 * g_xzz_yy_x_yz[i] * a_exp;

        g_z_0_0_0_xz_yy_x_zz[i] = -g_x_yy_x_zz[i] + 2.0 * g_xzz_yy_x_zz[i] * a_exp;
    }
    // integrals block (1572-1578)

    #pragma omp simd aligned(g_x_yy_y_xx, g_x_yy_y_xy, g_x_yy_y_xz, g_x_yy_y_yy, g_x_yy_y_yz, g_x_yy_y_zz, g_xzz_yy_y_xx, g_xzz_yy_y_xy, g_xzz_yy_y_xz, g_xzz_yy_y_yy, g_xzz_yy_y_yz, g_xzz_yy_y_zz, g_z_0_0_0_xz_yy_y_xx, g_z_0_0_0_xz_yy_y_xy, g_z_0_0_0_xz_yy_y_xz, g_z_0_0_0_xz_yy_y_yy, g_z_0_0_0_xz_yy_y_yz, g_z_0_0_0_xz_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_yy_y_xx[i] = -g_x_yy_y_xx[i] + 2.0 * g_xzz_yy_y_xx[i] * a_exp;

        g_z_0_0_0_xz_yy_y_xy[i] = -g_x_yy_y_xy[i] + 2.0 * g_xzz_yy_y_xy[i] * a_exp;

        g_z_0_0_0_xz_yy_y_xz[i] = -g_x_yy_y_xz[i] + 2.0 * g_xzz_yy_y_xz[i] * a_exp;

        g_z_0_0_0_xz_yy_y_yy[i] = -g_x_yy_y_yy[i] + 2.0 * g_xzz_yy_y_yy[i] * a_exp;

        g_z_0_0_0_xz_yy_y_yz[i] = -g_x_yy_y_yz[i] + 2.0 * g_xzz_yy_y_yz[i] * a_exp;

        g_z_0_0_0_xz_yy_y_zz[i] = -g_x_yy_y_zz[i] + 2.0 * g_xzz_yy_y_zz[i] * a_exp;
    }
    // integrals block (1578-1584)

    #pragma omp simd aligned(g_x_yy_z_xx, g_x_yy_z_xy, g_x_yy_z_xz, g_x_yy_z_yy, g_x_yy_z_yz, g_x_yy_z_zz, g_xzz_yy_z_xx, g_xzz_yy_z_xy, g_xzz_yy_z_xz, g_xzz_yy_z_yy, g_xzz_yy_z_yz, g_xzz_yy_z_zz, g_z_0_0_0_xz_yy_z_xx, g_z_0_0_0_xz_yy_z_xy, g_z_0_0_0_xz_yy_z_xz, g_z_0_0_0_xz_yy_z_yy, g_z_0_0_0_xz_yy_z_yz, g_z_0_0_0_xz_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_yy_z_xx[i] = -g_x_yy_z_xx[i] + 2.0 * g_xzz_yy_z_xx[i] * a_exp;

        g_z_0_0_0_xz_yy_z_xy[i] = -g_x_yy_z_xy[i] + 2.0 * g_xzz_yy_z_xy[i] * a_exp;

        g_z_0_0_0_xz_yy_z_xz[i] = -g_x_yy_z_xz[i] + 2.0 * g_xzz_yy_z_xz[i] * a_exp;

        g_z_0_0_0_xz_yy_z_yy[i] = -g_x_yy_z_yy[i] + 2.0 * g_xzz_yy_z_yy[i] * a_exp;

        g_z_0_0_0_xz_yy_z_yz[i] = -g_x_yy_z_yz[i] + 2.0 * g_xzz_yy_z_yz[i] * a_exp;

        g_z_0_0_0_xz_yy_z_zz[i] = -g_x_yy_z_zz[i] + 2.0 * g_xzz_yy_z_zz[i] * a_exp;
    }
    // integrals block (1584-1590)

    #pragma omp simd aligned(g_x_yz_x_xx, g_x_yz_x_xy, g_x_yz_x_xz, g_x_yz_x_yy, g_x_yz_x_yz, g_x_yz_x_zz, g_xzz_yz_x_xx, g_xzz_yz_x_xy, g_xzz_yz_x_xz, g_xzz_yz_x_yy, g_xzz_yz_x_yz, g_xzz_yz_x_zz, g_z_0_0_0_xz_yz_x_xx, g_z_0_0_0_xz_yz_x_xy, g_z_0_0_0_xz_yz_x_xz, g_z_0_0_0_xz_yz_x_yy, g_z_0_0_0_xz_yz_x_yz, g_z_0_0_0_xz_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_yz_x_xx[i] = -g_x_yz_x_xx[i] + 2.0 * g_xzz_yz_x_xx[i] * a_exp;

        g_z_0_0_0_xz_yz_x_xy[i] = -g_x_yz_x_xy[i] + 2.0 * g_xzz_yz_x_xy[i] * a_exp;

        g_z_0_0_0_xz_yz_x_xz[i] = -g_x_yz_x_xz[i] + 2.0 * g_xzz_yz_x_xz[i] * a_exp;

        g_z_0_0_0_xz_yz_x_yy[i] = -g_x_yz_x_yy[i] + 2.0 * g_xzz_yz_x_yy[i] * a_exp;

        g_z_0_0_0_xz_yz_x_yz[i] = -g_x_yz_x_yz[i] + 2.0 * g_xzz_yz_x_yz[i] * a_exp;

        g_z_0_0_0_xz_yz_x_zz[i] = -g_x_yz_x_zz[i] + 2.0 * g_xzz_yz_x_zz[i] * a_exp;
    }
    // integrals block (1590-1596)

    #pragma omp simd aligned(g_x_yz_y_xx, g_x_yz_y_xy, g_x_yz_y_xz, g_x_yz_y_yy, g_x_yz_y_yz, g_x_yz_y_zz, g_xzz_yz_y_xx, g_xzz_yz_y_xy, g_xzz_yz_y_xz, g_xzz_yz_y_yy, g_xzz_yz_y_yz, g_xzz_yz_y_zz, g_z_0_0_0_xz_yz_y_xx, g_z_0_0_0_xz_yz_y_xy, g_z_0_0_0_xz_yz_y_xz, g_z_0_0_0_xz_yz_y_yy, g_z_0_0_0_xz_yz_y_yz, g_z_0_0_0_xz_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_yz_y_xx[i] = -g_x_yz_y_xx[i] + 2.0 * g_xzz_yz_y_xx[i] * a_exp;

        g_z_0_0_0_xz_yz_y_xy[i] = -g_x_yz_y_xy[i] + 2.0 * g_xzz_yz_y_xy[i] * a_exp;

        g_z_0_0_0_xz_yz_y_xz[i] = -g_x_yz_y_xz[i] + 2.0 * g_xzz_yz_y_xz[i] * a_exp;

        g_z_0_0_0_xz_yz_y_yy[i] = -g_x_yz_y_yy[i] + 2.0 * g_xzz_yz_y_yy[i] * a_exp;

        g_z_0_0_0_xz_yz_y_yz[i] = -g_x_yz_y_yz[i] + 2.0 * g_xzz_yz_y_yz[i] * a_exp;

        g_z_0_0_0_xz_yz_y_zz[i] = -g_x_yz_y_zz[i] + 2.0 * g_xzz_yz_y_zz[i] * a_exp;
    }
    // integrals block (1596-1602)

    #pragma omp simd aligned(g_x_yz_z_xx, g_x_yz_z_xy, g_x_yz_z_xz, g_x_yz_z_yy, g_x_yz_z_yz, g_x_yz_z_zz, g_xzz_yz_z_xx, g_xzz_yz_z_xy, g_xzz_yz_z_xz, g_xzz_yz_z_yy, g_xzz_yz_z_yz, g_xzz_yz_z_zz, g_z_0_0_0_xz_yz_z_xx, g_z_0_0_0_xz_yz_z_xy, g_z_0_0_0_xz_yz_z_xz, g_z_0_0_0_xz_yz_z_yy, g_z_0_0_0_xz_yz_z_yz, g_z_0_0_0_xz_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_yz_z_xx[i] = -g_x_yz_z_xx[i] + 2.0 * g_xzz_yz_z_xx[i] * a_exp;

        g_z_0_0_0_xz_yz_z_xy[i] = -g_x_yz_z_xy[i] + 2.0 * g_xzz_yz_z_xy[i] * a_exp;

        g_z_0_0_0_xz_yz_z_xz[i] = -g_x_yz_z_xz[i] + 2.0 * g_xzz_yz_z_xz[i] * a_exp;

        g_z_0_0_0_xz_yz_z_yy[i] = -g_x_yz_z_yy[i] + 2.0 * g_xzz_yz_z_yy[i] * a_exp;

        g_z_0_0_0_xz_yz_z_yz[i] = -g_x_yz_z_yz[i] + 2.0 * g_xzz_yz_z_yz[i] * a_exp;

        g_z_0_0_0_xz_yz_z_zz[i] = -g_x_yz_z_zz[i] + 2.0 * g_xzz_yz_z_zz[i] * a_exp;
    }
    // integrals block (1602-1608)

    #pragma omp simd aligned(g_x_zz_x_xx, g_x_zz_x_xy, g_x_zz_x_xz, g_x_zz_x_yy, g_x_zz_x_yz, g_x_zz_x_zz, g_xzz_zz_x_xx, g_xzz_zz_x_xy, g_xzz_zz_x_xz, g_xzz_zz_x_yy, g_xzz_zz_x_yz, g_xzz_zz_x_zz, g_z_0_0_0_xz_zz_x_xx, g_z_0_0_0_xz_zz_x_xy, g_z_0_0_0_xz_zz_x_xz, g_z_0_0_0_xz_zz_x_yy, g_z_0_0_0_xz_zz_x_yz, g_z_0_0_0_xz_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_zz_x_xx[i] = -g_x_zz_x_xx[i] + 2.0 * g_xzz_zz_x_xx[i] * a_exp;

        g_z_0_0_0_xz_zz_x_xy[i] = -g_x_zz_x_xy[i] + 2.0 * g_xzz_zz_x_xy[i] * a_exp;

        g_z_0_0_0_xz_zz_x_xz[i] = -g_x_zz_x_xz[i] + 2.0 * g_xzz_zz_x_xz[i] * a_exp;

        g_z_0_0_0_xz_zz_x_yy[i] = -g_x_zz_x_yy[i] + 2.0 * g_xzz_zz_x_yy[i] * a_exp;

        g_z_0_0_0_xz_zz_x_yz[i] = -g_x_zz_x_yz[i] + 2.0 * g_xzz_zz_x_yz[i] * a_exp;

        g_z_0_0_0_xz_zz_x_zz[i] = -g_x_zz_x_zz[i] + 2.0 * g_xzz_zz_x_zz[i] * a_exp;
    }
    // integrals block (1608-1614)

    #pragma omp simd aligned(g_x_zz_y_xx, g_x_zz_y_xy, g_x_zz_y_xz, g_x_zz_y_yy, g_x_zz_y_yz, g_x_zz_y_zz, g_xzz_zz_y_xx, g_xzz_zz_y_xy, g_xzz_zz_y_xz, g_xzz_zz_y_yy, g_xzz_zz_y_yz, g_xzz_zz_y_zz, g_z_0_0_0_xz_zz_y_xx, g_z_0_0_0_xz_zz_y_xy, g_z_0_0_0_xz_zz_y_xz, g_z_0_0_0_xz_zz_y_yy, g_z_0_0_0_xz_zz_y_yz, g_z_0_0_0_xz_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_zz_y_xx[i] = -g_x_zz_y_xx[i] + 2.0 * g_xzz_zz_y_xx[i] * a_exp;

        g_z_0_0_0_xz_zz_y_xy[i] = -g_x_zz_y_xy[i] + 2.0 * g_xzz_zz_y_xy[i] * a_exp;

        g_z_0_0_0_xz_zz_y_xz[i] = -g_x_zz_y_xz[i] + 2.0 * g_xzz_zz_y_xz[i] * a_exp;

        g_z_0_0_0_xz_zz_y_yy[i] = -g_x_zz_y_yy[i] + 2.0 * g_xzz_zz_y_yy[i] * a_exp;

        g_z_0_0_0_xz_zz_y_yz[i] = -g_x_zz_y_yz[i] + 2.0 * g_xzz_zz_y_yz[i] * a_exp;

        g_z_0_0_0_xz_zz_y_zz[i] = -g_x_zz_y_zz[i] + 2.0 * g_xzz_zz_y_zz[i] * a_exp;
    }
    // integrals block (1614-1620)

    #pragma omp simd aligned(g_x_zz_z_xx, g_x_zz_z_xy, g_x_zz_z_xz, g_x_zz_z_yy, g_x_zz_z_yz, g_x_zz_z_zz, g_xzz_zz_z_xx, g_xzz_zz_z_xy, g_xzz_zz_z_xz, g_xzz_zz_z_yy, g_xzz_zz_z_yz, g_xzz_zz_z_zz, g_z_0_0_0_xz_zz_z_xx, g_z_0_0_0_xz_zz_z_xy, g_z_0_0_0_xz_zz_z_xz, g_z_0_0_0_xz_zz_z_yy, g_z_0_0_0_xz_zz_z_yz, g_z_0_0_0_xz_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_zz_z_xx[i] = -g_x_zz_z_xx[i] + 2.0 * g_xzz_zz_z_xx[i] * a_exp;

        g_z_0_0_0_xz_zz_z_xy[i] = -g_x_zz_z_xy[i] + 2.0 * g_xzz_zz_z_xy[i] * a_exp;

        g_z_0_0_0_xz_zz_z_xz[i] = -g_x_zz_z_xz[i] + 2.0 * g_xzz_zz_z_xz[i] * a_exp;

        g_z_0_0_0_xz_zz_z_yy[i] = -g_x_zz_z_yy[i] + 2.0 * g_xzz_zz_z_yy[i] * a_exp;

        g_z_0_0_0_xz_zz_z_yz[i] = -g_x_zz_z_yz[i] + 2.0 * g_xzz_zz_z_yz[i] * a_exp;

        g_z_0_0_0_xz_zz_z_zz[i] = -g_x_zz_z_zz[i] + 2.0 * g_xzz_zz_z_zz[i] * a_exp;
    }
    // integrals block (1620-1626)

    #pragma omp simd aligned(g_yyz_xx_x_xx, g_yyz_xx_x_xy, g_yyz_xx_x_xz, g_yyz_xx_x_yy, g_yyz_xx_x_yz, g_yyz_xx_x_zz, g_z_0_0_0_yy_xx_x_xx, g_z_0_0_0_yy_xx_x_xy, g_z_0_0_0_yy_xx_x_xz, g_z_0_0_0_yy_xx_x_yy, g_z_0_0_0_yy_xx_x_yz, g_z_0_0_0_yy_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_xx_x_xx[i] = 2.0 * g_yyz_xx_x_xx[i] * a_exp;

        g_z_0_0_0_yy_xx_x_xy[i] = 2.0 * g_yyz_xx_x_xy[i] * a_exp;

        g_z_0_0_0_yy_xx_x_xz[i] = 2.0 * g_yyz_xx_x_xz[i] * a_exp;

        g_z_0_0_0_yy_xx_x_yy[i] = 2.0 * g_yyz_xx_x_yy[i] * a_exp;

        g_z_0_0_0_yy_xx_x_yz[i] = 2.0 * g_yyz_xx_x_yz[i] * a_exp;

        g_z_0_0_0_yy_xx_x_zz[i] = 2.0 * g_yyz_xx_x_zz[i] * a_exp;
    }
    // integrals block (1626-1632)

    #pragma omp simd aligned(g_yyz_xx_y_xx, g_yyz_xx_y_xy, g_yyz_xx_y_xz, g_yyz_xx_y_yy, g_yyz_xx_y_yz, g_yyz_xx_y_zz, g_z_0_0_0_yy_xx_y_xx, g_z_0_0_0_yy_xx_y_xy, g_z_0_0_0_yy_xx_y_xz, g_z_0_0_0_yy_xx_y_yy, g_z_0_0_0_yy_xx_y_yz, g_z_0_0_0_yy_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_xx_y_xx[i] = 2.0 * g_yyz_xx_y_xx[i] * a_exp;

        g_z_0_0_0_yy_xx_y_xy[i] = 2.0 * g_yyz_xx_y_xy[i] * a_exp;

        g_z_0_0_0_yy_xx_y_xz[i] = 2.0 * g_yyz_xx_y_xz[i] * a_exp;

        g_z_0_0_0_yy_xx_y_yy[i] = 2.0 * g_yyz_xx_y_yy[i] * a_exp;

        g_z_0_0_0_yy_xx_y_yz[i] = 2.0 * g_yyz_xx_y_yz[i] * a_exp;

        g_z_0_0_0_yy_xx_y_zz[i] = 2.0 * g_yyz_xx_y_zz[i] * a_exp;
    }
    // integrals block (1632-1638)

    #pragma omp simd aligned(g_yyz_xx_z_xx, g_yyz_xx_z_xy, g_yyz_xx_z_xz, g_yyz_xx_z_yy, g_yyz_xx_z_yz, g_yyz_xx_z_zz, g_z_0_0_0_yy_xx_z_xx, g_z_0_0_0_yy_xx_z_xy, g_z_0_0_0_yy_xx_z_xz, g_z_0_0_0_yy_xx_z_yy, g_z_0_0_0_yy_xx_z_yz, g_z_0_0_0_yy_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_xx_z_xx[i] = 2.0 * g_yyz_xx_z_xx[i] * a_exp;

        g_z_0_0_0_yy_xx_z_xy[i] = 2.0 * g_yyz_xx_z_xy[i] * a_exp;

        g_z_0_0_0_yy_xx_z_xz[i] = 2.0 * g_yyz_xx_z_xz[i] * a_exp;

        g_z_0_0_0_yy_xx_z_yy[i] = 2.0 * g_yyz_xx_z_yy[i] * a_exp;

        g_z_0_0_0_yy_xx_z_yz[i] = 2.0 * g_yyz_xx_z_yz[i] * a_exp;

        g_z_0_0_0_yy_xx_z_zz[i] = 2.0 * g_yyz_xx_z_zz[i] * a_exp;
    }
    // integrals block (1638-1644)

    #pragma omp simd aligned(g_yyz_xy_x_xx, g_yyz_xy_x_xy, g_yyz_xy_x_xz, g_yyz_xy_x_yy, g_yyz_xy_x_yz, g_yyz_xy_x_zz, g_z_0_0_0_yy_xy_x_xx, g_z_0_0_0_yy_xy_x_xy, g_z_0_0_0_yy_xy_x_xz, g_z_0_0_0_yy_xy_x_yy, g_z_0_0_0_yy_xy_x_yz, g_z_0_0_0_yy_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_xy_x_xx[i] = 2.0 * g_yyz_xy_x_xx[i] * a_exp;

        g_z_0_0_0_yy_xy_x_xy[i] = 2.0 * g_yyz_xy_x_xy[i] * a_exp;

        g_z_0_0_0_yy_xy_x_xz[i] = 2.0 * g_yyz_xy_x_xz[i] * a_exp;

        g_z_0_0_0_yy_xy_x_yy[i] = 2.0 * g_yyz_xy_x_yy[i] * a_exp;

        g_z_0_0_0_yy_xy_x_yz[i] = 2.0 * g_yyz_xy_x_yz[i] * a_exp;

        g_z_0_0_0_yy_xy_x_zz[i] = 2.0 * g_yyz_xy_x_zz[i] * a_exp;
    }
    // integrals block (1644-1650)

    #pragma omp simd aligned(g_yyz_xy_y_xx, g_yyz_xy_y_xy, g_yyz_xy_y_xz, g_yyz_xy_y_yy, g_yyz_xy_y_yz, g_yyz_xy_y_zz, g_z_0_0_0_yy_xy_y_xx, g_z_0_0_0_yy_xy_y_xy, g_z_0_0_0_yy_xy_y_xz, g_z_0_0_0_yy_xy_y_yy, g_z_0_0_0_yy_xy_y_yz, g_z_0_0_0_yy_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_xy_y_xx[i] = 2.0 * g_yyz_xy_y_xx[i] * a_exp;

        g_z_0_0_0_yy_xy_y_xy[i] = 2.0 * g_yyz_xy_y_xy[i] * a_exp;

        g_z_0_0_0_yy_xy_y_xz[i] = 2.0 * g_yyz_xy_y_xz[i] * a_exp;

        g_z_0_0_0_yy_xy_y_yy[i] = 2.0 * g_yyz_xy_y_yy[i] * a_exp;

        g_z_0_0_0_yy_xy_y_yz[i] = 2.0 * g_yyz_xy_y_yz[i] * a_exp;

        g_z_0_0_0_yy_xy_y_zz[i] = 2.0 * g_yyz_xy_y_zz[i] * a_exp;
    }
    // integrals block (1650-1656)

    #pragma omp simd aligned(g_yyz_xy_z_xx, g_yyz_xy_z_xy, g_yyz_xy_z_xz, g_yyz_xy_z_yy, g_yyz_xy_z_yz, g_yyz_xy_z_zz, g_z_0_0_0_yy_xy_z_xx, g_z_0_0_0_yy_xy_z_xy, g_z_0_0_0_yy_xy_z_xz, g_z_0_0_0_yy_xy_z_yy, g_z_0_0_0_yy_xy_z_yz, g_z_0_0_0_yy_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_xy_z_xx[i] = 2.0 * g_yyz_xy_z_xx[i] * a_exp;

        g_z_0_0_0_yy_xy_z_xy[i] = 2.0 * g_yyz_xy_z_xy[i] * a_exp;

        g_z_0_0_0_yy_xy_z_xz[i] = 2.0 * g_yyz_xy_z_xz[i] * a_exp;

        g_z_0_0_0_yy_xy_z_yy[i] = 2.0 * g_yyz_xy_z_yy[i] * a_exp;

        g_z_0_0_0_yy_xy_z_yz[i] = 2.0 * g_yyz_xy_z_yz[i] * a_exp;

        g_z_0_0_0_yy_xy_z_zz[i] = 2.0 * g_yyz_xy_z_zz[i] * a_exp;
    }
    // integrals block (1656-1662)

    #pragma omp simd aligned(g_yyz_xz_x_xx, g_yyz_xz_x_xy, g_yyz_xz_x_xz, g_yyz_xz_x_yy, g_yyz_xz_x_yz, g_yyz_xz_x_zz, g_z_0_0_0_yy_xz_x_xx, g_z_0_0_0_yy_xz_x_xy, g_z_0_0_0_yy_xz_x_xz, g_z_0_0_0_yy_xz_x_yy, g_z_0_0_0_yy_xz_x_yz, g_z_0_0_0_yy_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_xz_x_xx[i] = 2.0 * g_yyz_xz_x_xx[i] * a_exp;

        g_z_0_0_0_yy_xz_x_xy[i] = 2.0 * g_yyz_xz_x_xy[i] * a_exp;

        g_z_0_0_0_yy_xz_x_xz[i] = 2.0 * g_yyz_xz_x_xz[i] * a_exp;

        g_z_0_0_0_yy_xz_x_yy[i] = 2.0 * g_yyz_xz_x_yy[i] * a_exp;

        g_z_0_0_0_yy_xz_x_yz[i] = 2.0 * g_yyz_xz_x_yz[i] * a_exp;

        g_z_0_0_0_yy_xz_x_zz[i] = 2.0 * g_yyz_xz_x_zz[i] * a_exp;
    }
    // integrals block (1662-1668)

    #pragma omp simd aligned(g_yyz_xz_y_xx, g_yyz_xz_y_xy, g_yyz_xz_y_xz, g_yyz_xz_y_yy, g_yyz_xz_y_yz, g_yyz_xz_y_zz, g_z_0_0_0_yy_xz_y_xx, g_z_0_0_0_yy_xz_y_xy, g_z_0_0_0_yy_xz_y_xz, g_z_0_0_0_yy_xz_y_yy, g_z_0_0_0_yy_xz_y_yz, g_z_0_0_0_yy_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_xz_y_xx[i] = 2.0 * g_yyz_xz_y_xx[i] * a_exp;

        g_z_0_0_0_yy_xz_y_xy[i] = 2.0 * g_yyz_xz_y_xy[i] * a_exp;

        g_z_0_0_0_yy_xz_y_xz[i] = 2.0 * g_yyz_xz_y_xz[i] * a_exp;

        g_z_0_0_0_yy_xz_y_yy[i] = 2.0 * g_yyz_xz_y_yy[i] * a_exp;

        g_z_0_0_0_yy_xz_y_yz[i] = 2.0 * g_yyz_xz_y_yz[i] * a_exp;

        g_z_0_0_0_yy_xz_y_zz[i] = 2.0 * g_yyz_xz_y_zz[i] * a_exp;
    }
    // integrals block (1668-1674)

    #pragma omp simd aligned(g_yyz_xz_z_xx, g_yyz_xz_z_xy, g_yyz_xz_z_xz, g_yyz_xz_z_yy, g_yyz_xz_z_yz, g_yyz_xz_z_zz, g_z_0_0_0_yy_xz_z_xx, g_z_0_0_0_yy_xz_z_xy, g_z_0_0_0_yy_xz_z_xz, g_z_0_0_0_yy_xz_z_yy, g_z_0_0_0_yy_xz_z_yz, g_z_0_0_0_yy_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_xz_z_xx[i] = 2.0 * g_yyz_xz_z_xx[i] * a_exp;

        g_z_0_0_0_yy_xz_z_xy[i] = 2.0 * g_yyz_xz_z_xy[i] * a_exp;

        g_z_0_0_0_yy_xz_z_xz[i] = 2.0 * g_yyz_xz_z_xz[i] * a_exp;

        g_z_0_0_0_yy_xz_z_yy[i] = 2.0 * g_yyz_xz_z_yy[i] * a_exp;

        g_z_0_0_0_yy_xz_z_yz[i] = 2.0 * g_yyz_xz_z_yz[i] * a_exp;

        g_z_0_0_0_yy_xz_z_zz[i] = 2.0 * g_yyz_xz_z_zz[i] * a_exp;
    }
    // integrals block (1674-1680)

    #pragma omp simd aligned(g_yyz_yy_x_xx, g_yyz_yy_x_xy, g_yyz_yy_x_xz, g_yyz_yy_x_yy, g_yyz_yy_x_yz, g_yyz_yy_x_zz, g_z_0_0_0_yy_yy_x_xx, g_z_0_0_0_yy_yy_x_xy, g_z_0_0_0_yy_yy_x_xz, g_z_0_0_0_yy_yy_x_yy, g_z_0_0_0_yy_yy_x_yz, g_z_0_0_0_yy_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_yy_x_xx[i] = 2.0 * g_yyz_yy_x_xx[i] * a_exp;

        g_z_0_0_0_yy_yy_x_xy[i] = 2.0 * g_yyz_yy_x_xy[i] * a_exp;

        g_z_0_0_0_yy_yy_x_xz[i] = 2.0 * g_yyz_yy_x_xz[i] * a_exp;

        g_z_0_0_0_yy_yy_x_yy[i] = 2.0 * g_yyz_yy_x_yy[i] * a_exp;

        g_z_0_0_0_yy_yy_x_yz[i] = 2.0 * g_yyz_yy_x_yz[i] * a_exp;

        g_z_0_0_0_yy_yy_x_zz[i] = 2.0 * g_yyz_yy_x_zz[i] * a_exp;
    }
    // integrals block (1680-1686)

    #pragma omp simd aligned(g_yyz_yy_y_xx, g_yyz_yy_y_xy, g_yyz_yy_y_xz, g_yyz_yy_y_yy, g_yyz_yy_y_yz, g_yyz_yy_y_zz, g_z_0_0_0_yy_yy_y_xx, g_z_0_0_0_yy_yy_y_xy, g_z_0_0_0_yy_yy_y_xz, g_z_0_0_0_yy_yy_y_yy, g_z_0_0_0_yy_yy_y_yz, g_z_0_0_0_yy_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_yy_y_xx[i] = 2.0 * g_yyz_yy_y_xx[i] * a_exp;

        g_z_0_0_0_yy_yy_y_xy[i] = 2.0 * g_yyz_yy_y_xy[i] * a_exp;

        g_z_0_0_0_yy_yy_y_xz[i] = 2.0 * g_yyz_yy_y_xz[i] * a_exp;

        g_z_0_0_0_yy_yy_y_yy[i] = 2.0 * g_yyz_yy_y_yy[i] * a_exp;

        g_z_0_0_0_yy_yy_y_yz[i] = 2.0 * g_yyz_yy_y_yz[i] * a_exp;

        g_z_0_0_0_yy_yy_y_zz[i] = 2.0 * g_yyz_yy_y_zz[i] * a_exp;
    }
    // integrals block (1686-1692)

    #pragma omp simd aligned(g_yyz_yy_z_xx, g_yyz_yy_z_xy, g_yyz_yy_z_xz, g_yyz_yy_z_yy, g_yyz_yy_z_yz, g_yyz_yy_z_zz, g_z_0_0_0_yy_yy_z_xx, g_z_0_0_0_yy_yy_z_xy, g_z_0_0_0_yy_yy_z_xz, g_z_0_0_0_yy_yy_z_yy, g_z_0_0_0_yy_yy_z_yz, g_z_0_0_0_yy_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_yy_z_xx[i] = 2.0 * g_yyz_yy_z_xx[i] * a_exp;

        g_z_0_0_0_yy_yy_z_xy[i] = 2.0 * g_yyz_yy_z_xy[i] * a_exp;

        g_z_0_0_0_yy_yy_z_xz[i] = 2.0 * g_yyz_yy_z_xz[i] * a_exp;

        g_z_0_0_0_yy_yy_z_yy[i] = 2.0 * g_yyz_yy_z_yy[i] * a_exp;

        g_z_0_0_0_yy_yy_z_yz[i] = 2.0 * g_yyz_yy_z_yz[i] * a_exp;

        g_z_0_0_0_yy_yy_z_zz[i] = 2.0 * g_yyz_yy_z_zz[i] * a_exp;
    }
    // integrals block (1692-1698)

    #pragma omp simd aligned(g_yyz_yz_x_xx, g_yyz_yz_x_xy, g_yyz_yz_x_xz, g_yyz_yz_x_yy, g_yyz_yz_x_yz, g_yyz_yz_x_zz, g_z_0_0_0_yy_yz_x_xx, g_z_0_0_0_yy_yz_x_xy, g_z_0_0_0_yy_yz_x_xz, g_z_0_0_0_yy_yz_x_yy, g_z_0_0_0_yy_yz_x_yz, g_z_0_0_0_yy_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_yz_x_xx[i] = 2.0 * g_yyz_yz_x_xx[i] * a_exp;

        g_z_0_0_0_yy_yz_x_xy[i] = 2.0 * g_yyz_yz_x_xy[i] * a_exp;

        g_z_0_0_0_yy_yz_x_xz[i] = 2.0 * g_yyz_yz_x_xz[i] * a_exp;

        g_z_0_0_0_yy_yz_x_yy[i] = 2.0 * g_yyz_yz_x_yy[i] * a_exp;

        g_z_0_0_0_yy_yz_x_yz[i] = 2.0 * g_yyz_yz_x_yz[i] * a_exp;

        g_z_0_0_0_yy_yz_x_zz[i] = 2.0 * g_yyz_yz_x_zz[i] * a_exp;
    }
    // integrals block (1698-1704)

    #pragma omp simd aligned(g_yyz_yz_y_xx, g_yyz_yz_y_xy, g_yyz_yz_y_xz, g_yyz_yz_y_yy, g_yyz_yz_y_yz, g_yyz_yz_y_zz, g_z_0_0_0_yy_yz_y_xx, g_z_0_0_0_yy_yz_y_xy, g_z_0_0_0_yy_yz_y_xz, g_z_0_0_0_yy_yz_y_yy, g_z_0_0_0_yy_yz_y_yz, g_z_0_0_0_yy_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_yz_y_xx[i] = 2.0 * g_yyz_yz_y_xx[i] * a_exp;

        g_z_0_0_0_yy_yz_y_xy[i] = 2.0 * g_yyz_yz_y_xy[i] * a_exp;

        g_z_0_0_0_yy_yz_y_xz[i] = 2.0 * g_yyz_yz_y_xz[i] * a_exp;

        g_z_0_0_0_yy_yz_y_yy[i] = 2.0 * g_yyz_yz_y_yy[i] * a_exp;

        g_z_0_0_0_yy_yz_y_yz[i] = 2.0 * g_yyz_yz_y_yz[i] * a_exp;

        g_z_0_0_0_yy_yz_y_zz[i] = 2.0 * g_yyz_yz_y_zz[i] * a_exp;
    }
    // integrals block (1704-1710)

    #pragma omp simd aligned(g_yyz_yz_z_xx, g_yyz_yz_z_xy, g_yyz_yz_z_xz, g_yyz_yz_z_yy, g_yyz_yz_z_yz, g_yyz_yz_z_zz, g_z_0_0_0_yy_yz_z_xx, g_z_0_0_0_yy_yz_z_xy, g_z_0_0_0_yy_yz_z_xz, g_z_0_0_0_yy_yz_z_yy, g_z_0_0_0_yy_yz_z_yz, g_z_0_0_0_yy_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_yz_z_xx[i] = 2.0 * g_yyz_yz_z_xx[i] * a_exp;

        g_z_0_0_0_yy_yz_z_xy[i] = 2.0 * g_yyz_yz_z_xy[i] * a_exp;

        g_z_0_0_0_yy_yz_z_xz[i] = 2.0 * g_yyz_yz_z_xz[i] * a_exp;

        g_z_0_0_0_yy_yz_z_yy[i] = 2.0 * g_yyz_yz_z_yy[i] * a_exp;

        g_z_0_0_0_yy_yz_z_yz[i] = 2.0 * g_yyz_yz_z_yz[i] * a_exp;

        g_z_0_0_0_yy_yz_z_zz[i] = 2.0 * g_yyz_yz_z_zz[i] * a_exp;
    }
    // integrals block (1710-1716)

    #pragma omp simd aligned(g_yyz_zz_x_xx, g_yyz_zz_x_xy, g_yyz_zz_x_xz, g_yyz_zz_x_yy, g_yyz_zz_x_yz, g_yyz_zz_x_zz, g_z_0_0_0_yy_zz_x_xx, g_z_0_0_0_yy_zz_x_xy, g_z_0_0_0_yy_zz_x_xz, g_z_0_0_0_yy_zz_x_yy, g_z_0_0_0_yy_zz_x_yz, g_z_0_0_0_yy_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_zz_x_xx[i] = 2.0 * g_yyz_zz_x_xx[i] * a_exp;

        g_z_0_0_0_yy_zz_x_xy[i] = 2.0 * g_yyz_zz_x_xy[i] * a_exp;

        g_z_0_0_0_yy_zz_x_xz[i] = 2.0 * g_yyz_zz_x_xz[i] * a_exp;

        g_z_0_0_0_yy_zz_x_yy[i] = 2.0 * g_yyz_zz_x_yy[i] * a_exp;

        g_z_0_0_0_yy_zz_x_yz[i] = 2.0 * g_yyz_zz_x_yz[i] * a_exp;

        g_z_0_0_0_yy_zz_x_zz[i] = 2.0 * g_yyz_zz_x_zz[i] * a_exp;
    }
    // integrals block (1716-1722)

    #pragma omp simd aligned(g_yyz_zz_y_xx, g_yyz_zz_y_xy, g_yyz_zz_y_xz, g_yyz_zz_y_yy, g_yyz_zz_y_yz, g_yyz_zz_y_zz, g_z_0_0_0_yy_zz_y_xx, g_z_0_0_0_yy_zz_y_xy, g_z_0_0_0_yy_zz_y_xz, g_z_0_0_0_yy_zz_y_yy, g_z_0_0_0_yy_zz_y_yz, g_z_0_0_0_yy_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_zz_y_xx[i] = 2.0 * g_yyz_zz_y_xx[i] * a_exp;

        g_z_0_0_0_yy_zz_y_xy[i] = 2.0 * g_yyz_zz_y_xy[i] * a_exp;

        g_z_0_0_0_yy_zz_y_xz[i] = 2.0 * g_yyz_zz_y_xz[i] * a_exp;

        g_z_0_0_0_yy_zz_y_yy[i] = 2.0 * g_yyz_zz_y_yy[i] * a_exp;

        g_z_0_0_0_yy_zz_y_yz[i] = 2.0 * g_yyz_zz_y_yz[i] * a_exp;

        g_z_0_0_0_yy_zz_y_zz[i] = 2.0 * g_yyz_zz_y_zz[i] * a_exp;
    }
    // integrals block (1722-1728)

    #pragma omp simd aligned(g_yyz_zz_z_xx, g_yyz_zz_z_xy, g_yyz_zz_z_xz, g_yyz_zz_z_yy, g_yyz_zz_z_yz, g_yyz_zz_z_zz, g_z_0_0_0_yy_zz_z_xx, g_z_0_0_0_yy_zz_z_xy, g_z_0_0_0_yy_zz_z_xz, g_z_0_0_0_yy_zz_z_yy, g_z_0_0_0_yy_zz_z_yz, g_z_0_0_0_yy_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_zz_z_xx[i] = 2.0 * g_yyz_zz_z_xx[i] * a_exp;

        g_z_0_0_0_yy_zz_z_xy[i] = 2.0 * g_yyz_zz_z_xy[i] * a_exp;

        g_z_0_0_0_yy_zz_z_xz[i] = 2.0 * g_yyz_zz_z_xz[i] * a_exp;

        g_z_0_0_0_yy_zz_z_yy[i] = 2.0 * g_yyz_zz_z_yy[i] * a_exp;

        g_z_0_0_0_yy_zz_z_yz[i] = 2.0 * g_yyz_zz_z_yz[i] * a_exp;

        g_z_0_0_0_yy_zz_z_zz[i] = 2.0 * g_yyz_zz_z_zz[i] * a_exp;
    }
    // integrals block (1728-1734)

    #pragma omp simd aligned(g_y_xx_x_xx, g_y_xx_x_xy, g_y_xx_x_xz, g_y_xx_x_yy, g_y_xx_x_yz, g_y_xx_x_zz, g_yzz_xx_x_xx, g_yzz_xx_x_xy, g_yzz_xx_x_xz, g_yzz_xx_x_yy, g_yzz_xx_x_yz, g_yzz_xx_x_zz, g_z_0_0_0_yz_xx_x_xx, g_z_0_0_0_yz_xx_x_xy, g_z_0_0_0_yz_xx_x_xz, g_z_0_0_0_yz_xx_x_yy, g_z_0_0_0_yz_xx_x_yz, g_z_0_0_0_yz_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_xx_x_xx[i] = -g_y_xx_x_xx[i] + 2.0 * g_yzz_xx_x_xx[i] * a_exp;

        g_z_0_0_0_yz_xx_x_xy[i] = -g_y_xx_x_xy[i] + 2.0 * g_yzz_xx_x_xy[i] * a_exp;

        g_z_0_0_0_yz_xx_x_xz[i] = -g_y_xx_x_xz[i] + 2.0 * g_yzz_xx_x_xz[i] * a_exp;

        g_z_0_0_0_yz_xx_x_yy[i] = -g_y_xx_x_yy[i] + 2.0 * g_yzz_xx_x_yy[i] * a_exp;

        g_z_0_0_0_yz_xx_x_yz[i] = -g_y_xx_x_yz[i] + 2.0 * g_yzz_xx_x_yz[i] * a_exp;

        g_z_0_0_0_yz_xx_x_zz[i] = -g_y_xx_x_zz[i] + 2.0 * g_yzz_xx_x_zz[i] * a_exp;
    }
    // integrals block (1734-1740)

    #pragma omp simd aligned(g_y_xx_y_xx, g_y_xx_y_xy, g_y_xx_y_xz, g_y_xx_y_yy, g_y_xx_y_yz, g_y_xx_y_zz, g_yzz_xx_y_xx, g_yzz_xx_y_xy, g_yzz_xx_y_xz, g_yzz_xx_y_yy, g_yzz_xx_y_yz, g_yzz_xx_y_zz, g_z_0_0_0_yz_xx_y_xx, g_z_0_0_0_yz_xx_y_xy, g_z_0_0_0_yz_xx_y_xz, g_z_0_0_0_yz_xx_y_yy, g_z_0_0_0_yz_xx_y_yz, g_z_0_0_0_yz_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_xx_y_xx[i] = -g_y_xx_y_xx[i] + 2.0 * g_yzz_xx_y_xx[i] * a_exp;

        g_z_0_0_0_yz_xx_y_xy[i] = -g_y_xx_y_xy[i] + 2.0 * g_yzz_xx_y_xy[i] * a_exp;

        g_z_0_0_0_yz_xx_y_xz[i] = -g_y_xx_y_xz[i] + 2.0 * g_yzz_xx_y_xz[i] * a_exp;

        g_z_0_0_0_yz_xx_y_yy[i] = -g_y_xx_y_yy[i] + 2.0 * g_yzz_xx_y_yy[i] * a_exp;

        g_z_0_0_0_yz_xx_y_yz[i] = -g_y_xx_y_yz[i] + 2.0 * g_yzz_xx_y_yz[i] * a_exp;

        g_z_0_0_0_yz_xx_y_zz[i] = -g_y_xx_y_zz[i] + 2.0 * g_yzz_xx_y_zz[i] * a_exp;
    }
    // integrals block (1740-1746)

    #pragma omp simd aligned(g_y_xx_z_xx, g_y_xx_z_xy, g_y_xx_z_xz, g_y_xx_z_yy, g_y_xx_z_yz, g_y_xx_z_zz, g_yzz_xx_z_xx, g_yzz_xx_z_xy, g_yzz_xx_z_xz, g_yzz_xx_z_yy, g_yzz_xx_z_yz, g_yzz_xx_z_zz, g_z_0_0_0_yz_xx_z_xx, g_z_0_0_0_yz_xx_z_xy, g_z_0_0_0_yz_xx_z_xz, g_z_0_0_0_yz_xx_z_yy, g_z_0_0_0_yz_xx_z_yz, g_z_0_0_0_yz_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_xx_z_xx[i] = -g_y_xx_z_xx[i] + 2.0 * g_yzz_xx_z_xx[i] * a_exp;

        g_z_0_0_0_yz_xx_z_xy[i] = -g_y_xx_z_xy[i] + 2.0 * g_yzz_xx_z_xy[i] * a_exp;

        g_z_0_0_0_yz_xx_z_xz[i] = -g_y_xx_z_xz[i] + 2.0 * g_yzz_xx_z_xz[i] * a_exp;

        g_z_0_0_0_yz_xx_z_yy[i] = -g_y_xx_z_yy[i] + 2.0 * g_yzz_xx_z_yy[i] * a_exp;

        g_z_0_0_0_yz_xx_z_yz[i] = -g_y_xx_z_yz[i] + 2.0 * g_yzz_xx_z_yz[i] * a_exp;

        g_z_0_0_0_yz_xx_z_zz[i] = -g_y_xx_z_zz[i] + 2.0 * g_yzz_xx_z_zz[i] * a_exp;
    }
    // integrals block (1746-1752)

    #pragma omp simd aligned(g_y_xy_x_xx, g_y_xy_x_xy, g_y_xy_x_xz, g_y_xy_x_yy, g_y_xy_x_yz, g_y_xy_x_zz, g_yzz_xy_x_xx, g_yzz_xy_x_xy, g_yzz_xy_x_xz, g_yzz_xy_x_yy, g_yzz_xy_x_yz, g_yzz_xy_x_zz, g_z_0_0_0_yz_xy_x_xx, g_z_0_0_0_yz_xy_x_xy, g_z_0_0_0_yz_xy_x_xz, g_z_0_0_0_yz_xy_x_yy, g_z_0_0_0_yz_xy_x_yz, g_z_0_0_0_yz_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_xy_x_xx[i] = -g_y_xy_x_xx[i] + 2.0 * g_yzz_xy_x_xx[i] * a_exp;

        g_z_0_0_0_yz_xy_x_xy[i] = -g_y_xy_x_xy[i] + 2.0 * g_yzz_xy_x_xy[i] * a_exp;

        g_z_0_0_0_yz_xy_x_xz[i] = -g_y_xy_x_xz[i] + 2.0 * g_yzz_xy_x_xz[i] * a_exp;

        g_z_0_0_0_yz_xy_x_yy[i] = -g_y_xy_x_yy[i] + 2.0 * g_yzz_xy_x_yy[i] * a_exp;

        g_z_0_0_0_yz_xy_x_yz[i] = -g_y_xy_x_yz[i] + 2.0 * g_yzz_xy_x_yz[i] * a_exp;

        g_z_0_0_0_yz_xy_x_zz[i] = -g_y_xy_x_zz[i] + 2.0 * g_yzz_xy_x_zz[i] * a_exp;
    }
    // integrals block (1752-1758)

    #pragma omp simd aligned(g_y_xy_y_xx, g_y_xy_y_xy, g_y_xy_y_xz, g_y_xy_y_yy, g_y_xy_y_yz, g_y_xy_y_zz, g_yzz_xy_y_xx, g_yzz_xy_y_xy, g_yzz_xy_y_xz, g_yzz_xy_y_yy, g_yzz_xy_y_yz, g_yzz_xy_y_zz, g_z_0_0_0_yz_xy_y_xx, g_z_0_0_0_yz_xy_y_xy, g_z_0_0_0_yz_xy_y_xz, g_z_0_0_0_yz_xy_y_yy, g_z_0_0_0_yz_xy_y_yz, g_z_0_0_0_yz_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_xy_y_xx[i] = -g_y_xy_y_xx[i] + 2.0 * g_yzz_xy_y_xx[i] * a_exp;

        g_z_0_0_0_yz_xy_y_xy[i] = -g_y_xy_y_xy[i] + 2.0 * g_yzz_xy_y_xy[i] * a_exp;

        g_z_0_0_0_yz_xy_y_xz[i] = -g_y_xy_y_xz[i] + 2.0 * g_yzz_xy_y_xz[i] * a_exp;

        g_z_0_0_0_yz_xy_y_yy[i] = -g_y_xy_y_yy[i] + 2.0 * g_yzz_xy_y_yy[i] * a_exp;

        g_z_0_0_0_yz_xy_y_yz[i] = -g_y_xy_y_yz[i] + 2.0 * g_yzz_xy_y_yz[i] * a_exp;

        g_z_0_0_0_yz_xy_y_zz[i] = -g_y_xy_y_zz[i] + 2.0 * g_yzz_xy_y_zz[i] * a_exp;
    }
    // integrals block (1758-1764)

    #pragma omp simd aligned(g_y_xy_z_xx, g_y_xy_z_xy, g_y_xy_z_xz, g_y_xy_z_yy, g_y_xy_z_yz, g_y_xy_z_zz, g_yzz_xy_z_xx, g_yzz_xy_z_xy, g_yzz_xy_z_xz, g_yzz_xy_z_yy, g_yzz_xy_z_yz, g_yzz_xy_z_zz, g_z_0_0_0_yz_xy_z_xx, g_z_0_0_0_yz_xy_z_xy, g_z_0_0_0_yz_xy_z_xz, g_z_0_0_0_yz_xy_z_yy, g_z_0_0_0_yz_xy_z_yz, g_z_0_0_0_yz_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_xy_z_xx[i] = -g_y_xy_z_xx[i] + 2.0 * g_yzz_xy_z_xx[i] * a_exp;

        g_z_0_0_0_yz_xy_z_xy[i] = -g_y_xy_z_xy[i] + 2.0 * g_yzz_xy_z_xy[i] * a_exp;

        g_z_0_0_0_yz_xy_z_xz[i] = -g_y_xy_z_xz[i] + 2.0 * g_yzz_xy_z_xz[i] * a_exp;

        g_z_0_0_0_yz_xy_z_yy[i] = -g_y_xy_z_yy[i] + 2.0 * g_yzz_xy_z_yy[i] * a_exp;

        g_z_0_0_0_yz_xy_z_yz[i] = -g_y_xy_z_yz[i] + 2.0 * g_yzz_xy_z_yz[i] * a_exp;

        g_z_0_0_0_yz_xy_z_zz[i] = -g_y_xy_z_zz[i] + 2.0 * g_yzz_xy_z_zz[i] * a_exp;
    }
    // integrals block (1764-1770)

    #pragma omp simd aligned(g_y_xz_x_xx, g_y_xz_x_xy, g_y_xz_x_xz, g_y_xz_x_yy, g_y_xz_x_yz, g_y_xz_x_zz, g_yzz_xz_x_xx, g_yzz_xz_x_xy, g_yzz_xz_x_xz, g_yzz_xz_x_yy, g_yzz_xz_x_yz, g_yzz_xz_x_zz, g_z_0_0_0_yz_xz_x_xx, g_z_0_0_0_yz_xz_x_xy, g_z_0_0_0_yz_xz_x_xz, g_z_0_0_0_yz_xz_x_yy, g_z_0_0_0_yz_xz_x_yz, g_z_0_0_0_yz_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_xz_x_xx[i] = -g_y_xz_x_xx[i] + 2.0 * g_yzz_xz_x_xx[i] * a_exp;

        g_z_0_0_0_yz_xz_x_xy[i] = -g_y_xz_x_xy[i] + 2.0 * g_yzz_xz_x_xy[i] * a_exp;

        g_z_0_0_0_yz_xz_x_xz[i] = -g_y_xz_x_xz[i] + 2.0 * g_yzz_xz_x_xz[i] * a_exp;

        g_z_0_0_0_yz_xz_x_yy[i] = -g_y_xz_x_yy[i] + 2.0 * g_yzz_xz_x_yy[i] * a_exp;

        g_z_0_0_0_yz_xz_x_yz[i] = -g_y_xz_x_yz[i] + 2.0 * g_yzz_xz_x_yz[i] * a_exp;

        g_z_0_0_0_yz_xz_x_zz[i] = -g_y_xz_x_zz[i] + 2.0 * g_yzz_xz_x_zz[i] * a_exp;
    }
    // integrals block (1770-1776)

    #pragma omp simd aligned(g_y_xz_y_xx, g_y_xz_y_xy, g_y_xz_y_xz, g_y_xz_y_yy, g_y_xz_y_yz, g_y_xz_y_zz, g_yzz_xz_y_xx, g_yzz_xz_y_xy, g_yzz_xz_y_xz, g_yzz_xz_y_yy, g_yzz_xz_y_yz, g_yzz_xz_y_zz, g_z_0_0_0_yz_xz_y_xx, g_z_0_0_0_yz_xz_y_xy, g_z_0_0_0_yz_xz_y_xz, g_z_0_0_0_yz_xz_y_yy, g_z_0_0_0_yz_xz_y_yz, g_z_0_0_0_yz_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_xz_y_xx[i] = -g_y_xz_y_xx[i] + 2.0 * g_yzz_xz_y_xx[i] * a_exp;

        g_z_0_0_0_yz_xz_y_xy[i] = -g_y_xz_y_xy[i] + 2.0 * g_yzz_xz_y_xy[i] * a_exp;

        g_z_0_0_0_yz_xz_y_xz[i] = -g_y_xz_y_xz[i] + 2.0 * g_yzz_xz_y_xz[i] * a_exp;

        g_z_0_0_0_yz_xz_y_yy[i] = -g_y_xz_y_yy[i] + 2.0 * g_yzz_xz_y_yy[i] * a_exp;

        g_z_0_0_0_yz_xz_y_yz[i] = -g_y_xz_y_yz[i] + 2.0 * g_yzz_xz_y_yz[i] * a_exp;

        g_z_0_0_0_yz_xz_y_zz[i] = -g_y_xz_y_zz[i] + 2.0 * g_yzz_xz_y_zz[i] * a_exp;
    }
    // integrals block (1776-1782)

    #pragma omp simd aligned(g_y_xz_z_xx, g_y_xz_z_xy, g_y_xz_z_xz, g_y_xz_z_yy, g_y_xz_z_yz, g_y_xz_z_zz, g_yzz_xz_z_xx, g_yzz_xz_z_xy, g_yzz_xz_z_xz, g_yzz_xz_z_yy, g_yzz_xz_z_yz, g_yzz_xz_z_zz, g_z_0_0_0_yz_xz_z_xx, g_z_0_0_0_yz_xz_z_xy, g_z_0_0_0_yz_xz_z_xz, g_z_0_0_0_yz_xz_z_yy, g_z_0_0_0_yz_xz_z_yz, g_z_0_0_0_yz_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_xz_z_xx[i] = -g_y_xz_z_xx[i] + 2.0 * g_yzz_xz_z_xx[i] * a_exp;

        g_z_0_0_0_yz_xz_z_xy[i] = -g_y_xz_z_xy[i] + 2.0 * g_yzz_xz_z_xy[i] * a_exp;

        g_z_0_0_0_yz_xz_z_xz[i] = -g_y_xz_z_xz[i] + 2.0 * g_yzz_xz_z_xz[i] * a_exp;

        g_z_0_0_0_yz_xz_z_yy[i] = -g_y_xz_z_yy[i] + 2.0 * g_yzz_xz_z_yy[i] * a_exp;

        g_z_0_0_0_yz_xz_z_yz[i] = -g_y_xz_z_yz[i] + 2.0 * g_yzz_xz_z_yz[i] * a_exp;

        g_z_0_0_0_yz_xz_z_zz[i] = -g_y_xz_z_zz[i] + 2.0 * g_yzz_xz_z_zz[i] * a_exp;
    }
    // integrals block (1782-1788)

    #pragma omp simd aligned(g_y_yy_x_xx, g_y_yy_x_xy, g_y_yy_x_xz, g_y_yy_x_yy, g_y_yy_x_yz, g_y_yy_x_zz, g_yzz_yy_x_xx, g_yzz_yy_x_xy, g_yzz_yy_x_xz, g_yzz_yy_x_yy, g_yzz_yy_x_yz, g_yzz_yy_x_zz, g_z_0_0_0_yz_yy_x_xx, g_z_0_0_0_yz_yy_x_xy, g_z_0_0_0_yz_yy_x_xz, g_z_0_0_0_yz_yy_x_yy, g_z_0_0_0_yz_yy_x_yz, g_z_0_0_0_yz_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_yy_x_xx[i] = -g_y_yy_x_xx[i] + 2.0 * g_yzz_yy_x_xx[i] * a_exp;

        g_z_0_0_0_yz_yy_x_xy[i] = -g_y_yy_x_xy[i] + 2.0 * g_yzz_yy_x_xy[i] * a_exp;

        g_z_0_0_0_yz_yy_x_xz[i] = -g_y_yy_x_xz[i] + 2.0 * g_yzz_yy_x_xz[i] * a_exp;

        g_z_0_0_0_yz_yy_x_yy[i] = -g_y_yy_x_yy[i] + 2.0 * g_yzz_yy_x_yy[i] * a_exp;

        g_z_0_0_0_yz_yy_x_yz[i] = -g_y_yy_x_yz[i] + 2.0 * g_yzz_yy_x_yz[i] * a_exp;

        g_z_0_0_0_yz_yy_x_zz[i] = -g_y_yy_x_zz[i] + 2.0 * g_yzz_yy_x_zz[i] * a_exp;
    }
    // integrals block (1788-1794)

    #pragma omp simd aligned(g_y_yy_y_xx, g_y_yy_y_xy, g_y_yy_y_xz, g_y_yy_y_yy, g_y_yy_y_yz, g_y_yy_y_zz, g_yzz_yy_y_xx, g_yzz_yy_y_xy, g_yzz_yy_y_xz, g_yzz_yy_y_yy, g_yzz_yy_y_yz, g_yzz_yy_y_zz, g_z_0_0_0_yz_yy_y_xx, g_z_0_0_0_yz_yy_y_xy, g_z_0_0_0_yz_yy_y_xz, g_z_0_0_0_yz_yy_y_yy, g_z_0_0_0_yz_yy_y_yz, g_z_0_0_0_yz_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_yy_y_xx[i] = -g_y_yy_y_xx[i] + 2.0 * g_yzz_yy_y_xx[i] * a_exp;

        g_z_0_0_0_yz_yy_y_xy[i] = -g_y_yy_y_xy[i] + 2.0 * g_yzz_yy_y_xy[i] * a_exp;

        g_z_0_0_0_yz_yy_y_xz[i] = -g_y_yy_y_xz[i] + 2.0 * g_yzz_yy_y_xz[i] * a_exp;

        g_z_0_0_0_yz_yy_y_yy[i] = -g_y_yy_y_yy[i] + 2.0 * g_yzz_yy_y_yy[i] * a_exp;

        g_z_0_0_0_yz_yy_y_yz[i] = -g_y_yy_y_yz[i] + 2.0 * g_yzz_yy_y_yz[i] * a_exp;

        g_z_0_0_0_yz_yy_y_zz[i] = -g_y_yy_y_zz[i] + 2.0 * g_yzz_yy_y_zz[i] * a_exp;
    }
    // integrals block (1794-1800)

    #pragma omp simd aligned(g_y_yy_z_xx, g_y_yy_z_xy, g_y_yy_z_xz, g_y_yy_z_yy, g_y_yy_z_yz, g_y_yy_z_zz, g_yzz_yy_z_xx, g_yzz_yy_z_xy, g_yzz_yy_z_xz, g_yzz_yy_z_yy, g_yzz_yy_z_yz, g_yzz_yy_z_zz, g_z_0_0_0_yz_yy_z_xx, g_z_0_0_0_yz_yy_z_xy, g_z_0_0_0_yz_yy_z_xz, g_z_0_0_0_yz_yy_z_yy, g_z_0_0_0_yz_yy_z_yz, g_z_0_0_0_yz_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_yy_z_xx[i] = -g_y_yy_z_xx[i] + 2.0 * g_yzz_yy_z_xx[i] * a_exp;

        g_z_0_0_0_yz_yy_z_xy[i] = -g_y_yy_z_xy[i] + 2.0 * g_yzz_yy_z_xy[i] * a_exp;

        g_z_0_0_0_yz_yy_z_xz[i] = -g_y_yy_z_xz[i] + 2.0 * g_yzz_yy_z_xz[i] * a_exp;

        g_z_0_0_0_yz_yy_z_yy[i] = -g_y_yy_z_yy[i] + 2.0 * g_yzz_yy_z_yy[i] * a_exp;

        g_z_0_0_0_yz_yy_z_yz[i] = -g_y_yy_z_yz[i] + 2.0 * g_yzz_yy_z_yz[i] * a_exp;

        g_z_0_0_0_yz_yy_z_zz[i] = -g_y_yy_z_zz[i] + 2.0 * g_yzz_yy_z_zz[i] * a_exp;
    }
    // integrals block (1800-1806)

    #pragma omp simd aligned(g_y_yz_x_xx, g_y_yz_x_xy, g_y_yz_x_xz, g_y_yz_x_yy, g_y_yz_x_yz, g_y_yz_x_zz, g_yzz_yz_x_xx, g_yzz_yz_x_xy, g_yzz_yz_x_xz, g_yzz_yz_x_yy, g_yzz_yz_x_yz, g_yzz_yz_x_zz, g_z_0_0_0_yz_yz_x_xx, g_z_0_0_0_yz_yz_x_xy, g_z_0_0_0_yz_yz_x_xz, g_z_0_0_0_yz_yz_x_yy, g_z_0_0_0_yz_yz_x_yz, g_z_0_0_0_yz_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_yz_x_xx[i] = -g_y_yz_x_xx[i] + 2.0 * g_yzz_yz_x_xx[i] * a_exp;

        g_z_0_0_0_yz_yz_x_xy[i] = -g_y_yz_x_xy[i] + 2.0 * g_yzz_yz_x_xy[i] * a_exp;

        g_z_0_0_0_yz_yz_x_xz[i] = -g_y_yz_x_xz[i] + 2.0 * g_yzz_yz_x_xz[i] * a_exp;

        g_z_0_0_0_yz_yz_x_yy[i] = -g_y_yz_x_yy[i] + 2.0 * g_yzz_yz_x_yy[i] * a_exp;

        g_z_0_0_0_yz_yz_x_yz[i] = -g_y_yz_x_yz[i] + 2.0 * g_yzz_yz_x_yz[i] * a_exp;

        g_z_0_0_0_yz_yz_x_zz[i] = -g_y_yz_x_zz[i] + 2.0 * g_yzz_yz_x_zz[i] * a_exp;
    }
    // integrals block (1806-1812)

    #pragma omp simd aligned(g_y_yz_y_xx, g_y_yz_y_xy, g_y_yz_y_xz, g_y_yz_y_yy, g_y_yz_y_yz, g_y_yz_y_zz, g_yzz_yz_y_xx, g_yzz_yz_y_xy, g_yzz_yz_y_xz, g_yzz_yz_y_yy, g_yzz_yz_y_yz, g_yzz_yz_y_zz, g_z_0_0_0_yz_yz_y_xx, g_z_0_0_0_yz_yz_y_xy, g_z_0_0_0_yz_yz_y_xz, g_z_0_0_0_yz_yz_y_yy, g_z_0_0_0_yz_yz_y_yz, g_z_0_0_0_yz_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_yz_y_xx[i] = -g_y_yz_y_xx[i] + 2.0 * g_yzz_yz_y_xx[i] * a_exp;

        g_z_0_0_0_yz_yz_y_xy[i] = -g_y_yz_y_xy[i] + 2.0 * g_yzz_yz_y_xy[i] * a_exp;

        g_z_0_0_0_yz_yz_y_xz[i] = -g_y_yz_y_xz[i] + 2.0 * g_yzz_yz_y_xz[i] * a_exp;

        g_z_0_0_0_yz_yz_y_yy[i] = -g_y_yz_y_yy[i] + 2.0 * g_yzz_yz_y_yy[i] * a_exp;

        g_z_0_0_0_yz_yz_y_yz[i] = -g_y_yz_y_yz[i] + 2.0 * g_yzz_yz_y_yz[i] * a_exp;

        g_z_0_0_0_yz_yz_y_zz[i] = -g_y_yz_y_zz[i] + 2.0 * g_yzz_yz_y_zz[i] * a_exp;
    }
    // integrals block (1812-1818)

    #pragma omp simd aligned(g_y_yz_z_xx, g_y_yz_z_xy, g_y_yz_z_xz, g_y_yz_z_yy, g_y_yz_z_yz, g_y_yz_z_zz, g_yzz_yz_z_xx, g_yzz_yz_z_xy, g_yzz_yz_z_xz, g_yzz_yz_z_yy, g_yzz_yz_z_yz, g_yzz_yz_z_zz, g_z_0_0_0_yz_yz_z_xx, g_z_0_0_0_yz_yz_z_xy, g_z_0_0_0_yz_yz_z_xz, g_z_0_0_0_yz_yz_z_yy, g_z_0_0_0_yz_yz_z_yz, g_z_0_0_0_yz_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_yz_z_xx[i] = -g_y_yz_z_xx[i] + 2.0 * g_yzz_yz_z_xx[i] * a_exp;

        g_z_0_0_0_yz_yz_z_xy[i] = -g_y_yz_z_xy[i] + 2.0 * g_yzz_yz_z_xy[i] * a_exp;

        g_z_0_0_0_yz_yz_z_xz[i] = -g_y_yz_z_xz[i] + 2.0 * g_yzz_yz_z_xz[i] * a_exp;

        g_z_0_0_0_yz_yz_z_yy[i] = -g_y_yz_z_yy[i] + 2.0 * g_yzz_yz_z_yy[i] * a_exp;

        g_z_0_0_0_yz_yz_z_yz[i] = -g_y_yz_z_yz[i] + 2.0 * g_yzz_yz_z_yz[i] * a_exp;

        g_z_0_0_0_yz_yz_z_zz[i] = -g_y_yz_z_zz[i] + 2.0 * g_yzz_yz_z_zz[i] * a_exp;
    }
    // integrals block (1818-1824)

    #pragma omp simd aligned(g_y_zz_x_xx, g_y_zz_x_xy, g_y_zz_x_xz, g_y_zz_x_yy, g_y_zz_x_yz, g_y_zz_x_zz, g_yzz_zz_x_xx, g_yzz_zz_x_xy, g_yzz_zz_x_xz, g_yzz_zz_x_yy, g_yzz_zz_x_yz, g_yzz_zz_x_zz, g_z_0_0_0_yz_zz_x_xx, g_z_0_0_0_yz_zz_x_xy, g_z_0_0_0_yz_zz_x_xz, g_z_0_0_0_yz_zz_x_yy, g_z_0_0_0_yz_zz_x_yz, g_z_0_0_0_yz_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_zz_x_xx[i] = -g_y_zz_x_xx[i] + 2.0 * g_yzz_zz_x_xx[i] * a_exp;

        g_z_0_0_0_yz_zz_x_xy[i] = -g_y_zz_x_xy[i] + 2.0 * g_yzz_zz_x_xy[i] * a_exp;

        g_z_0_0_0_yz_zz_x_xz[i] = -g_y_zz_x_xz[i] + 2.0 * g_yzz_zz_x_xz[i] * a_exp;

        g_z_0_0_0_yz_zz_x_yy[i] = -g_y_zz_x_yy[i] + 2.0 * g_yzz_zz_x_yy[i] * a_exp;

        g_z_0_0_0_yz_zz_x_yz[i] = -g_y_zz_x_yz[i] + 2.0 * g_yzz_zz_x_yz[i] * a_exp;

        g_z_0_0_0_yz_zz_x_zz[i] = -g_y_zz_x_zz[i] + 2.0 * g_yzz_zz_x_zz[i] * a_exp;
    }
    // integrals block (1824-1830)

    #pragma omp simd aligned(g_y_zz_y_xx, g_y_zz_y_xy, g_y_zz_y_xz, g_y_zz_y_yy, g_y_zz_y_yz, g_y_zz_y_zz, g_yzz_zz_y_xx, g_yzz_zz_y_xy, g_yzz_zz_y_xz, g_yzz_zz_y_yy, g_yzz_zz_y_yz, g_yzz_zz_y_zz, g_z_0_0_0_yz_zz_y_xx, g_z_0_0_0_yz_zz_y_xy, g_z_0_0_0_yz_zz_y_xz, g_z_0_0_0_yz_zz_y_yy, g_z_0_0_0_yz_zz_y_yz, g_z_0_0_0_yz_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_zz_y_xx[i] = -g_y_zz_y_xx[i] + 2.0 * g_yzz_zz_y_xx[i] * a_exp;

        g_z_0_0_0_yz_zz_y_xy[i] = -g_y_zz_y_xy[i] + 2.0 * g_yzz_zz_y_xy[i] * a_exp;

        g_z_0_0_0_yz_zz_y_xz[i] = -g_y_zz_y_xz[i] + 2.0 * g_yzz_zz_y_xz[i] * a_exp;

        g_z_0_0_0_yz_zz_y_yy[i] = -g_y_zz_y_yy[i] + 2.0 * g_yzz_zz_y_yy[i] * a_exp;

        g_z_0_0_0_yz_zz_y_yz[i] = -g_y_zz_y_yz[i] + 2.0 * g_yzz_zz_y_yz[i] * a_exp;

        g_z_0_0_0_yz_zz_y_zz[i] = -g_y_zz_y_zz[i] + 2.0 * g_yzz_zz_y_zz[i] * a_exp;
    }
    // integrals block (1830-1836)

    #pragma omp simd aligned(g_y_zz_z_xx, g_y_zz_z_xy, g_y_zz_z_xz, g_y_zz_z_yy, g_y_zz_z_yz, g_y_zz_z_zz, g_yzz_zz_z_xx, g_yzz_zz_z_xy, g_yzz_zz_z_xz, g_yzz_zz_z_yy, g_yzz_zz_z_yz, g_yzz_zz_z_zz, g_z_0_0_0_yz_zz_z_xx, g_z_0_0_0_yz_zz_z_xy, g_z_0_0_0_yz_zz_z_xz, g_z_0_0_0_yz_zz_z_yy, g_z_0_0_0_yz_zz_z_yz, g_z_0_0_0_yz_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_zz_z_xx[i] = -g_y_zz_z_xx[i] + 2.0 * g_yzz_zz_z_xx[i] * a_exp;

        g_z_0_0_0_yz_zz_z_xy[i] = -g_y_zz_z_xy[i] + 2.0 * g_yzz_zz_z_xy[i] * a_exp;

        g_z_0_0_0_yz_zz_z_xz[i] = -g_y_zz_z_xz[i] + 2.0 * g_yzz_zz_z_xz[i] * a_exp;

        g_z_0_0_0_yz_zz_z_yy[i] = -g_y_zz_z_yy[i] + 2.0 * g_yzz_zz_z_yy[i] * a_exp;

        g_z_0_0_0_yz_zz_z_yz[i] = -g_y_zz_z_yz[i] + 2.0 * g_yzz_zz_z_yz[i] * a_exp;

        g_z_0_0_0_yz_zz_z_zz[i] = -g_y_zz_z_zz[i] + 2.0 * g_yzz_zz_z_zz[i] * a_exp;
    }
    // integrals block (1836-1842)

    #pragma omp simd aligned(g_z_0_0_0_zz_xx_x_xx, g_z_0_0_0_zz_xx_x_xy, g_z_0_0_0_zz_xx_x_xz, g_z_0_0_0_zz_xx_x_yy, g_z_0_0_0_zz_xx_x_yz, g_z_0_0_0_zz_xx_x_zz, g_z_xx_x_xx, g_z_xx_x_xy, g_z_xx_x_xz, g_z_xx_x_yy, g_z_xx_x_yz, g_z_xx_x_zz, g_zzz_xx_x_xx, g_zzz_xx_x_xy, g_zzz_xx_x_xz, g_zzz_xx_x_yy, g_zzz_xx_x_yz, g_zzz_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_xx_x_xx[i] = -2.0 * g_z_xx_x_xx[i] + 2.0 * g_zzz_xx_x_xx[i] * a_exp;

        g_z_0_0_0_zz_xx_x_xy[i] = -2.0 * g_z_xx_x_xy[i] + 2.0 * g_zzz_xx_x_xy[i] * a_exp;

        g_z_0_0_0_zz_xx_x_xz[i] = -2.0 * g_z_xx_x_xz[i] + 2.0 * g_zzz_xx_x_xz[i] * a_exp;

        g_z_0_0_0_zz_xx_x_yy[i] = -2.0 * g_z_xx_x_yy[i] + 2.0 * g_zzz_xx_x_yy[i] * a_exp;

        g_z_0_0_0_zz_xx_x_yz[i] = -2.0 * g_z_xx_x_yz[i] + 2.0 * g_zzz_xx_x_yz[i] * a_exp;

        g_z_0_0_0_zz_xx_x_zz[i] = -2.0 * g_z_xx_x_zz[i] + 2.0 * g_zzz_xx_x_zz[i] * a_exp;
    }
    // integrals block (1842-1848)

    #pragma omp simd aligned(g_z_0_0_0_zz_xx_y_xx, g_z_0_0_0_zz_xx_y_xy, g_z_0_0_0_zz_xx_y_xz, g_z_0_0_0_zz_xx_y_yy, g_z_0_0_0_zz_xx_y_yz, g_z_0_0_0_zz_xx_y_zz, g_z_xx_y_xx, g_z_xx_y_xy, g_z_xx_y_xz, g_z_xx_y_yy, g_z_xx_y_yz, g_z_xx_y_zz, g_zzz_xx_y_xx, g_zzz_xx_y_xy, g_zzz_xx_y_xz, g_zzz_xx_y_yy, g_zzz_xx_y_yz, g_zzz_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_xx_y_xx[i] = -2.0 * g_z_xx_y_xx[i] + 2.0 * g_zzz_xx_y_xx[i] * a_exp;

        g_z_0_0_0_zz_xx_y_xy[i] = -2.0 * g_z_xx_y_xy[i] + 2.0 * g_zzz_xx_y_xy[i] * a_exp;

        g_z_0_0_0_zz_xx_y_xz[i] = -2.0 * g_z_xx_y_xz[i] + 2.0 * g_zzz_xx_y_xz[i] * a_exp;

        g_z_0_0_0_zz_xx_y_yy[i] = -2.0 * g_z_xx_y_yy[i] + 2.0 * g_zzz_xx_y_yy[i] * a_exp;

        g_z_0_0_0_zz_xx_y_yz[i] = -2.0 * g_z_xx_y_yz[i] + 2.0 * g_zzz_xx_y_yz[i] * a_exp;

        g_z_0_0_0_zz_xx_y_zz[i] = -2.0 * g_z_xx_y_zz[i] + 2.0 * g_zzz_xx_y_zz[i] * a_exp;
    }
    // integrals block (1848-1854)

    #pragma omp simd aligned(g_z_0_0_0_zz_xx_z_xx, g_z_0_0_0_zz_xx_z_xy, g_z_0_0_0_zz_xx_z_xz, g_z_0_0_0_zz_xx_z_yy, g_z_0_0_0_zz_xx_z_yz, g_z_0_0_0_zz_xx_z_zz, g_z_xx_z_xx, g_z_xx_z_xy, g_z_xx_z_xz, g_z_xx_z_yy, g_z_xx_z_yz, g_z_xx_z_zz, g_zzz_xx_z_xx, g_zzz_xx_z_xy, g_zzz_xx_z_xz, g_zzz_xx_z_yy, g_zzz_xx_z_yz, g_zzz_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_xx_z_xx[i] = -2.0 * g_z_xx_z_xx[i] + 2.0 * g_zzz_xx_z_xx[i] * a_exp;

        g_z_0_0_0_zz_xx_z_xy[i] = -2.0 * g_z_xx_z_xy[i] + 2.0 * g_zzz_xx_z_xy[i] * a_exp;

        g_z_0_0_0_zz_xx_z_xz[i] = -2.0 * g_z_xx_z_xz[i] + 2.0 * g_zzz_xx_z_xz[i] * a_exp;

        g_z_0_0_0_zz_xx_z_yy[i] = -2.0 * g_z_xx_z_yy[i] + 2.0 * g_zzz_xx_z_yy[i] * a_exp;

        g_z_0_0_0_zz_xx_z_yz[i] = -2.0 * g_z_xx_z_yz[i] + 2.0 * g_zzz_xx_z_yz[i] * a_exp;

        g_z_0_0_0_zz_xx_z_zz[i] = -2.0 * g_z_xx_z_zz[i] + 2.0 * g_zzz_xx_z_zz[i] * a_exp;
    }
    // integrals block (1854-1860)

    #pragma omp simd aligned(g_z_0_0_0_zz_xy_x_xx, g_z_0_0_0_zz_xy_x_xy, g_z_0_0_0_zz_xy_x_xz, g_z_0_0_0_zz_xy_x_yy, g_z_0_0_0_zz_xy_x_yz, g_z_0_0_0_zz_xy_x_zz, g_z_xy_x_xx, g_z_xy_x_xy, g_z_xy_x_xz, g_z_xy_x_yy, g_z_xy_x_yz, g_z_xy_x_zz, g_zzz_xy_x_xx, g_zzz_xy_x_xy, g_zzz_xy_x_xz, g_zzz_xy_x_yy, g_zzz_xy_x_yz, g_zzz_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_xy_x_xx[i] = -2.0 * g_z_xy_x_xx[i] + 2.0 * g_zzz_xy_x_xx[i] * a_exp;

        g_z_0_0_0_zz_xy_x_xy[i] = -2.0 * g_z_xy_x_xy[i] + 2.0 * g_zzz_xy_x_xy[i] * a_exp;

        g_z_0_0_0_zz_xy_x_xz[i] = -2.0 * g_z_xy_x_xz[i] + 2.0 * g_zzz_xy_x_xz[i] * a_exp;

        g_z_0_0_0_zz_xy_x_yy[i] = -2.0 * g_z_xy_x_yy[i] + 2.0 * g_zzz_xy_x_yy[i] * a_exp;

        g_z_0_0_0_zz_xy_x_yz[i] = -2.0 * g_z_xy_x_yz[i] + 2.0 * g_zzz_xy_x_yz[i] * a_exp;

        g_z_0_0_0_zz_xy_x_zz[i] = -2.0 * g_z_xy_x_zz[i] + 2.0 * g_zzz_xy_x_zz[i] * a_exp;
    }
    // integrals block (1860-1866)

    #pragma omp simd aligned(g_z_0_0_0_zz_xy_y_xx, g_z_0_0_0_zz_xy_y_xy, g_z_0_0_0_zz_xy_y_xz, g_z_0_0_0_zz_xy_y_yy, g_z_0_0_0_zz_xy_y_yz, g_z_0_0_0_zz_xy_y_zz, g_z_xy_y_xx, g_z_xy_y_xy, g_z_xy_y_xz, g_z_xy_y_yy, g_z_xy_y_yz, g_z_xy_y_zz, g_zzz_xy_y_xx, g_zzz_xy_y_xy, g_zzz_xy_y_xz, g_zzz_xy_y_yy, g_zzz_xy_y_yz, g_zzz_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_xy_y_xx[i] = -2.0 * g_z_xy_y_xx[i] + 2.0 * g_zzz_xy_y_xx[i] * a_exp;

        g_z_0_0_0_zz_xy_y_xy[i] = -2.0 * g_z_xy_y_xy[i] + 2.0 * g_zzz_xy_y_xy[i] * a_exp;

        g_z_0_0_0_zz_xy_y_xz[i] = -2.0 * g_z_xy_y_xz[i] + 2.0 * g_zzz_xy_y_xz[i] * a_exp;

        g_z_0_0_0_zz_xy_y_yy[i] = -2.0 * g_z_xy_y_yy[i] + 2.0 * g_zzz_xy_y_yy[i] * a_exp;

        g_z_0_0_0_zz_xy_y_yz[i] = -2.0 * g_z_xy_y_yz[i] + 2.0 * g_zzz_xy_y_yz[i] * a_exp;

        g_z_0_0_0_zz_xy_y_zz[i] = -2.0 * g_z_xy_y_zz[i] + 2.0 * g_zzz_xy_y_zz[i] * a_exp;
    }
    // integrals block (1866-1872)

    #pragma omp simd aligned(g_z_0_0_0_zz_xy_z_xx, g_z_0_0_0_zz_xy_z_xy, g_z_0_0_0_zz_xy_z_xz, g_z_0_0_0_zz_xy_z_yy, g_z_0_0_0_zz_xy_z_yz, g_z_0_0_0_zz_xy_z_zz, g_z_xy_z_xx, g_z_xy_z_xy, g_z_xy_z_xz, g_z_xy_z_yy, g_z_xy_z_yz, g_z_xy_z_zz, g_zzz_xy_z_xx, g_zzz_xy_z_xy, g_zzz_xy_z_xz, g_zzz_xy_z_yy, g_zzz_xy_z_yz, g_zzz_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_xy_z_xx[i] = -2.0 * g_z_xy_z_xx[i] + 2.0 * g_zzz_xy_z_xx[i] * a_exp;

        g_z_0_0_0_zz_xy_z_xy[i] = -2.0 * g_z_xy_z_xy[i] + 2.0 * g_zzz_xy_z_xy[i] * a_exp;

        g_z_0_0_0_zz_xy_z_xz[i] = -2.0 * g_z_xy_z_xz[i] + 2.0 * g_zzz_xy_z_xz[i] * a_exp;

        g_z_0_0_0_zz_xy_z_yy[i] = -2.0 * g_z_xy_z_yy[i] + 2.0 * g_zzz_xy_z_yy[i] * a_exp;

        g_z_0_0_0_zz_xy_z_yz[i] = -2.0 * g_z_xy_z_yz[i] + 2.0 * g_zzz_xy_z_yz[i] * a_exp;

        g_z_0_0_0_zz_xy_z_zz[i] = -2.0 * g_z_xy_z_zz[i] + 2.0 * g_zzz_xy_z_zz[i] * a_exp;
    }
    // integrals block (1872-1878)

    #pragma omp simd aligned(g_z_0_0_0_zz_xz_x_xx, g_z_0_0_0_zz_xz_x_xy, g_z_0_0_0_zz_xz_x_xz, g_z_0_0_0_zz_xz_x_yy, g_z_0_0_0_zz_xz_x_yz, g_z_0_0_0_zz_xz_x_zz, g_z_xz_x_xx, g_z_xz_x_xy, g_z_xz_x_xz, g_z_xz_x_yy, g_z_xz_x_yz, g_z_xz_x_zz, g_zzz_xz_x_xx, g_zzz_xz_x_xy, g_zzz_xz_x_xz, g_zzz_xz_x_yy, g_zzz_xz_x_yz, g_zzz_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_xz_x_xx[i] = -2.0 * g_z_xz_x_xx[i] + 2.0 * g_zzz_xz_x_xx[i] * a_exp;

        g_z_0_0_0_zz_xz_x_xy[i] = -2.0 * g_z_xz_x_xy[i] + 2.0 * g_zzz_xz_x_xy[i] * a_exp;

        g_z_0_0_0_zz_xz_x_xz[i] = -2.0 * g_z_xz_x_xz[i] + 2.0 * g_zzz_xz_x_xz[i] * a_exp;

        g_z_0_0_0_zz_xz_x_yy[i] = -2.0 * g_z_xz_x_yy[i] + 2.0 * g_zzz_xz_x_yy[i] * a_exp;

        g_z_0_0_0_zz_xz_x_yz[i] = -2.0 * g_z_xz_x_yz[i] + 2.0 * g_zzz_xz_x_yz[i] * a_exp;

        g_z_0_0_0_zz_xz_x_zz[i] = -2.0 * g_z_xz_x_zz[i] + 2.0 * g_zzz_xz_x_zz[i] * a_exp;
    }
    // integrals block (1878-1884)

    #pragma omp simd aligned(g_z_0_0_0_zz_xz_y_xx, g_z_0_0_0_zz_xz_y_xy, g_z_0_0_0_zz_xz_y_xz, g_z_0_0_0_zz_xz_y_yy, g_z_0_0_0_zz_xz_y_yz, g_z_0_0_0_zz_xz_y_zz, g_z_xz_y_xx, g_z_xz_y_xy, g_z_xz_y_xz, g_z_xz_y_yy, g_z_xz_y_yz, g_z_xz_y_zz, g_zzz_xz_y_xx, g_zzz_xz_y_xy, g_zzz_xz_y_xz, g_zzz_xz_y_yy, g_zzz_xz_y_yz, g_zzz_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_xz_y_xx[i] = -2.0 * g_z_xz_y_xx[i] + 2.0 * g_zzz_xz_y_xx[i] * a_exp;

        g_z_0_0_0_zz_xz_y_xy[i] = -2.0 * g_z_xz_y_xy[i] + 2.0 * g_zzz_xz_y_xy[i] * a_exp;

        g_z_0_0_0_zz_xz_y_xz[i] = -2.0 * g_z_xz_y_xz[i] + 2.0 * g_zzz_xz_y_xz[i] * a_exp;

        g_z_0_0_0_zz_xz_y_yy[i] = -2.0 * g_z_xz_y_yy[i] + 2.0 * g_zzz_xz_y_yy[i] * a_exp;

        g_z_0_0_0_zz_xz_y_yz[i] = -2.0 * g_z_xz_y_yz[i] + 2.0 * g_zzz_xz_y_yz[i] * a_exp;

        g_z_0_0_0_zz_xz_y_zz[i] = -2.0 * g_z_xz_y_zz[i] + 2.0 * g_zzz_xz_y_zz[i] * a_exp;
    }
    // integrals block (1884-1890)

    #pragma omp simd aligned(g_z_0_0_0_zz_xz_z_xx, g_z_0_0_0_zz_xz_z_xy, g_z_0_0_0_zz_xz_z_xz, g_z_0_0_0_zz_xz_z_yy, g_z_0_0_0_zz_xz_z_yz, g_z_0_0_0_zz_xz_z_zz, g_z_xz_z_xx, g_z_xz_z_xy, g_z_xz_z_xz, g_z_xz_z_yy, g_z_xz_z_yz, g_z_xz_z_zz, g_zzz_xz_z_xx, g_zzz_xz_z_xy, g_zzz_xz_z_xz, g_zzz_xz_z_yy, g_zzz_xz_z_yz, g_zzz_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_xz_z_xx[i] = -2.0 * g_z_xz_z_xx[i] + 2.0 * g_zzz_xz_z_xx[i] * a_exp;

        g_z_0_0_0_zz_xz_z_xy[i] = -2.0 * g_z_xz_z_xy[i] + 2.0 * g_zzz_xz_z_xy[i] * a_exp;

        g_z_0_0_0_zz_xz_z_xz[i] = -2.0 * g_z_xz_z_xz[i] + 2.0 * g_zzz_xz_z_xz[i] * a_exp;

        g_z_0_0_0_zz_xz_z_yy[i] = -2.0 * g_z_xz_z_yy[i] + 2.0 * g_zzz_xz_z_yy[i] * a_exp;

        g_z_0_0_0_zz_xz_z_yz[i] = -2.0 * g_z_xz_z_yz[i] + 2.0 * g_zzz_xz_z_yz[i] * a_exp;

        g_z_0_0_0_zz_xz_z_zz[i] = -2.0 * g_z_xz_z_zz[i] + 2.0 * g_zzz_xz_z_zz[i] * a_exp;
    }
    // integrals block (1890-1896)

    #pragma omp simd aligned(g_z_0_0_0_zz_yy_x_xx, g_z_0_0_0_zz_yy_x_xy, g_z_0_0_0_zz_yy_x_xz, g_z_0_0_0_zz_yy_x_yy, g_z_0_0_0_zz_yy_x_yz, g_z_0_0_0_zz_yy_x_zz, g_z_yy_x_xx, g_z_yy_x_xy, g_z_yy_x_xz, g_z_yy_x_yy, g_z_yy_x_yz, g_z_yy_x_zz, g_zzz_yy_x_xx, g_zzz_yy_x_xy, g_zzz_yy_x_xz, g_zzz_yy_x_yy, g_zzz_yy_x_yz, g_zzz_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_yy_x_xx[i] = -2.0 * g_z_yy_x_xx[i] + 2.0 * g_zzz_yy_x_xx[i] * a_exp;

        g_z_0_0_0_zz_yy_x_xy[i] = -2.0 * g_z_yy_x_xy[i] + 2.0 * g_zzz_yy_x_xy[i] * a_exp;

        g_z_0_0_0_zz_yy_x_xz[i] = -2.0 * g_z_yy_x_xz[i] + 2.0 * g_zzz_yy_x_xz[i] * a_exp;

        g_z_0_0_0_zz_yy_x_yy[i] = -2.0 * g_z_yy_x_yy[i] + 2.0 * g_zzz_yy_x_yy[i] * a_exp;

        g_z_0_0_0_zz_yy_x_yz[i] = -2.0 * g_z_yy_x_yz[i] + 2.0 * g_zzz_yy_x_yz[i] * a_exp;

        g_z_0_0_0_zz_yy_x_zz[i] = -2.0 * g_z_yy_x_zz[i] + 2.0 * g_zzz_yy_x_zz[i] * a_exp;
    }
    // integrals block (1896-1902)

    #pragma omp simd aligned(g_z_0_0_0_zz_yy_y_xx, g_z_0_0_0_zz_yy_y_xy, g_z_0_0_0_zz_yy_y_xz, g_z_0_0_0_zz_yy_y_yy, g_z_0_0_0_zz_yy_y_yz, g_z_0_0_0_zz_yy_y_zz, g_z_yy_y_xx, g_z_yy_y_xy, g_z_yy_y_xz, g_z_yy_y_yy, g_z_yy_y_yz, g_z_yy_y_zz, g_zzz_yy_y_xx, g_zzz_yy_y_xy, g_zzz_yy_y_xz, g_zzz_yy_y_yy, g_zzz_yy_y_yz, g_zzz_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_yy_y_xx[i] = -2.0 * g_z_yy_y_xx[i] + 2.0 * g_zzz_yy_y_xx[i] * a_exp;

        g_z_0_0_0_zz_yy_y_xy[i] = -2.0 * g_z_yy_y_xy[i] + 2.0 * g_zzz_yy_y_xy[i] * a_exp;

        g_z_0_0_0_zz_yy_y_xz[i] = -2.0 * g_z_yy_y_xz[i] + 2.0 * g_zzz_yy_y_xz[i] * a_exp;

        g_z_0_0_0_zz_yy_y_yy[i] = -2.0 * g_z_yy_y_yy[i] + 2.0 * g_zzz_yy_y_yy[i] * a_exp;

        g_z_0_0_0_zz_yy_y_yz[i] = -2.0 * g_z_yy_y_yz[i] + 2.0 * g_zzz_yy_y_yz[i] * a_exp;

        g_z_0_0_0_zz_yy_y_zz[i] = -2.0 * g_z_yy_y_zz[i] + 2.0 * g_zzz_yy_y_zz[i] * a_exp;
    }
    // integrals block (1902-1908)

    #pragma omp simd aligned(g_z_0_0_0_zz_yy_z_xx, g_z_0_0_0_zz_yy_z_xy, g_z_0_0_0_zz_yy_z_xz, g_z_0_0_0_zz_yy_z_yy, g_z_0_0_0_zz_yy_z_yz, g_z_0_0_0_zz_yy_z_zz, g_z_yy_z_xx, g_z_yy_z_xy, g_z_yy_z_xz, g_z_yy_z_yy, g_z_yy_z_yz, g_z_yy_z_zz, g_zzz_yy_z_xx, g_zzz_yy_z_xy, g_zzz_yy_z_xz, g_zzz_yy_z_yy, g_zzz_yy_z_yz, g_zzz_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_yy_z_xx[i] = -2.0 * g_z_yy_z_xx[i] + 2.0 * g_zzz_yy_z_xx[i] * a_exp;

        g_z_0_0_0_zz_yy_z_xy[i] = -2.0 * g_z_yy_z_xy[i] + 2.0 * g_zzz_yy_z_xy[i] * a_exp;

        g_z_0_0_0_zz_yy_z_xz[i] = -2.0 * g_z_yy_z_xz[i] + 2.0 * g_zzz_yy_z_xz[i] * a_exp;

        g_z_0_0_0_zz_yy_z_yy[i] = -2.0 * g_z_yy_z_yy[i] + 2.0 * g_zzz_yy_z_yy[i] * a_exp;

        g_z_0_0_0_zz_yy_z_yz[i] = -2.0 * g_z_yy_z_yz[i] + 2.0 * g_zzz_yy_z_yz[i] * a_exp;

        g_z_0_0_0_zz_yy_z_zz[i] = -2.0 * g_z_yy_z_zz[i] + 2.0 * g_zzz_yy_z_zz[i] * a_exp;
    }
    // integrals block (1908-1914)

    #pragma omp simd aligned(g_z_0_0_0_zz_yz_x_xx, g_z_0_0_0_zz_yz_x_xy, g_z_0_0_0_zz_yz_x_xz, g_z_0_0_0_zz_yz_x_yy, g_z_0_0_0_zz_yz_x_yz, g_z_0_0_0_zz_yz_x_zz, g_z_yz_x_xx, g_z_yz_x_xy, g_z_yz_x_xz, g_z_yz_x_yy, g_z_yz_x_yz, g_z_yz_x_zz, g_zzz_yz_x_xx, g_zzz_yz_x_xy, g_zzz_yz_x_xz, g_zzz_yz_x_yy, g_zzz_yz_x_yz, g_zzz_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_yz_x_xx[i] = -2.0 * g_z_yz_x_xx[i] + 2.0 * g_zzz_yz_x_xx[i] * a_exp;

        g_z_0_0_0_zz_yz_x_xy[i] = -2.0 * g_z_yz_x_xy[i] + 2.0 * g_zzz_yz_x_xy[i] * a_exp;

        g_z_0_0_0_zz_yz_x_xz[i] = -2.0 * g_z_yz_x_xz[i] + 2.0 * g_zzz_yz_x_xz[i] * a_exp;

        g_z_0_0_0_zz_yz_x_yy[i] = -2.0 * g_z_yz_x_yy[i] + 2.0 * g_zzz_yz_x_yy[i] * a_exp;

        g_z_0_0_0_zz_yz_x_yz[i] = -2.0 * g_z_yz_x_yz[i] + 2.0 * g_zzz_yz_x_yz[i] * a_exp;

        g_z_0_0_0_zz_yz_x_zz[i] = -2.0 * g_z_yz_x_zz[i] + 2.0 * g_zzz_yz_x_zz[i] * a_exp;
    }
    // integrals block (1914-1920)

    #pragma omp simd aligned(g_z_0_0_0_zz_yz_y_xx, g_z_0_0_0_zz_yz_y_xy, g_z_0_0_0_zz_yz_y_xz, g_z_0_0_0_zz_yz_y_yy, g_z_0_0_0_zz_yz_y_yz, g_z_0_0_0_zz_yz_y_zz, g_z_yz_y_xx, g_z_yz_y_xy, g_z_yz_y_xz, g_z_yz_y_yy, g_z_yz_y_yz, g_z_yz_y_zz, g_zzz_yz_y_xx, g_zzz_yz_y_xy, g_zzz_yz_y_xz, g_zzz_yz_y_yy, g_zzz_yz_y_yz, g_zzz_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_yz_y_xx[i] = -2.0 * g_z_yz_y_xx[i] + 2.0 * g_zzz_yz_y_xx[i] * a_exp;

        g_z_0_0_0_zz_yz_y_xy[i] = -2.0 * g_z_yz_y_xy[i] + 2.0 * g_zzz_yz_y_xy[i] * a_exp;

        g_z_0_0_0_zz_yz_y_xz[i] = -2.0 * g_z_yz_y_xz[i] + 2.0 * g_zzz_yz_y_xz[i] * a_exp;

        g_z_0_0_0_zz_yz_y_yy[i] = -2.0 * g_z_yz_y_yy[i] + 2.0 * g_zzz_yz_y_yy[i] * a_exp;

        g_z_0_0_0_zz_yz_y_yz[i] = -2.0 * g_z_yz_y_yz[i] + 2.0 * g_zzz_yz_y_yz[i] * a_exp;

        g_z_0_0_0_zz_yz_y_zz[i] = -2.0 * g_z_yz_y_zz[i] + 2.0 * g_zzz_yz_y_zz[i] * a_exp;
    }
    // integrals block (1920-1926)

    #pragma omp simd aligned(g_z_0_0_0_zz_yz_z_xx, g_z_0_0_0_zz_yz_z_xy, g_z_0_0_0_zz_yz_z_xz, g_z_0_0_0_zz_yz_z_yy, g_z_0_0_0_zz_yz_z_yz, g_z_0_0_0_zz_yz_z_zz, g_z_yz_z_xx, g_z_yz_z_xy, g_z_yz_z_xz, g_z_yz_z_yy, g_z_yz_z_yz, g_z_yz_z_zz, g_zzz_yz_z_xx, g_zzz_yz_z_xy, g_zzz_yz_z_xz, g_zzz_yz_z_yy, g_zzz_yz_z_yz, g_zzz_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_yz_z_xx[i] = -2.0 * g_z_yz_z_xx[i] + 2.0 * g_zzz_yz_z_xx[i] * a_exp;

        g_z_0_0_0_zz_yz_z_xy[i] = -2.0 * g_z_yz_z_xy[i] + 2.0 * g_zzz_yz_z_xy[i] * a_exp;

        g_z_0_0_0_zz_yz_z_xz[i] = -2.0 * g_z_yz_z_xz[i] + 2.0 * g_zzz_yz_z_xz[i] * a_exp;

        g_z_0_0_0_zz_yz_z_yy[i] = -2.0 * g_z_yz_z_yy[i] + 2.0 * g_zzz_yz_z_yy[i] * a_exp;

        g_z_0_0_0_zz_yz_z_yz[i] = -2.0 * g_z_yz_z_yz[i] + 2.0 * g_zzz_yz_z_yz[i] * a_exp;

        g_z_0_0_0_zz_yz_z_zz[i] = -2.0 * g_z_yz_z_zz[i] + 2.0 * g_zzz_yz_z_zz[i] * a_exp;
    }
    // integrals block (1926-1932)

    #pragma omp simd aligned(g_z_0_0_0_zz_zz_x_xx, g_z_0_0_0_zz_zz_x_xy, g_z_0_0_0_zz_zz_x_xz, g_z_0_0_0_zz_zz_x_yy, g_z_0_0_0_zz_zz_x_yz, g_z_0_0_0_zz_zz_x_zz, g_z_zz_x_xx, g_z_zz_x_xy, g_z_zz_x_xz, g_z_zz_x_yy, g_z_zz_x_yz, g_z_zz_x_zz, g_zzz_zz_x_xx, g_zzz_zz_x_xy, g_zzz_zz_x_xz, g_zzz_zz_x_yy, g_zzz_zz_x_yz, g_zzz_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_zz_x_xx[i] = -2.0 * g_z_zz_x_xx[i] + 2.0 * g_zzz_zz_x_xx[i] * a_exp;

        g_z_0_0_0_zz_zz_x_xy[i] = -2.0 * g_z_zz_x_xy[i] + 2.0 * g_zzz_zz_x_xy[i] * a_exp;

        g_z_0_0_0_zz_zz_x_xz[i] = -2.0 * g_z_zz_x_xz[i] + 2.0 * g_zzz_zz_x_xz[i] * a_exp;

        g_z_0_0_0_zz_zz_x_yy[i] = -2.0 * g_z_zz_x_yy[i] + 2.0 * g_zzz_zz_x_yy[i] * a_exp;

        g_z_0_0_0_zz_zz_x_yz[i] = -2.0 * g_z_zz_x_yz[i] + 2.0 * g_zzz_zz_x_yz[i] * a_exp;

        g_z_0_0_0_zz_zz_x_zz[i] = -2.0 * g_z_zz_x_zz[i] + 2.0 * g_zzz_zz_x_zz[i] * a_exp;
    }
    // integrals block (1932-1938)

    #pragma omp simd aligned(g_z_0_0_0_zz_zz_y_xx, g_z_0_0_0_zz_zz_y_xy, g_z_0_0_0_zz_zz_y_xz, g_z_0_0_0_zz_zz_y_yy, g_z_0_0_0_zz_zz_y_yz, g_z_0_0_0_zz_zz_y_zz, g_z_zz_y_xx, g_z_zz_y_xy, g_z_zz_y_xz, g_z_zz_y_yy, g_z_zz_y_yz, g_z_zz_y_zz, g_zzz_zz_y_xx, g_zzz_zz_y_xy, g_zzz_zz_y_xz, g_zzz_zz_y_yy, g_zzz_zz_y_yz, g_zzz_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_zz_y_xx[i] = -2.0 * g_z_zz_y_xx[i] + 2.0 * g_zzz_zz_y_xx[i] * a_exp;

        g_z_0_0_0_zz_zz_y_xy[i] = -2.0 * g_z_zz_y_xy[i] + 2.0 * g_zzz_zz_y_xy[i] * a_exp;

        g_z_0_0_0_zz_zz_y_xz[i] = -2.0 * g_z_zz_y_xz[i] + 2.0 * g_zzz_zz_y_xz[i] * a_exp;

        g_z_0_0_0_zz_zz_y_yy[i] = -2.0 * g_z_zz_y_yy[i] + 2.0 * g_zzz_zz_y_yy[i] * a_exp;

        g_z_0_0_0_zz_zz_y_yz[i] = -2.0 * g_z_zz_y_yz[i] + 2.0 * g_zzz_zz_y_yz[i] * a_exp;

        g_z_0_0_0_zz_zz_y_zz[i] = -2.0 * g_z_zz_y_zz[i] + 2.0 * g_zzz_zz_y_zz[i] * a_exp;
    }
    // integrals block (1938-1944)

    #pragma omp simd aligned(g_z_0_0_0_zz_zz_z_xx, g_z_0_0_0_zz_zz_z_xy, g_z_0_0_0_zz_zz_z_xz, g_z_0_0_0_zz_zz_z_yy, g_z_0_0_0_zz_zz_z_yz, g_z_0_0_0_zz_zz_z_zz, g_z_zz_z_xx, g_z_zz_z_xy, g_z_zz_z_xz, g_z_zz_z_yy, g_z_zz_z_yz, g_z_zz_z_zz, g_zzz_zz_z_xx, g_zzz_zz_z_xy, g_zzz_zz_z_xz, g_zzz_zz_z_yy, g_zzz_zz_z_yz, g_zzz_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_zz_z_xx[i] = -2.0 * g_z_zz_z_xx[i] + 2.0 * g_zzz_zz_z_xx[i] * a_exp;

        g_z_0_0_0_zz_zz_z_xy[i] = -2.0 * g_z_zz_z_xy[i] + 2.0 * g_zzz_zz_z_xy[i] * a_exp;

        g_z_0_0_0_zz_zz_z_xz[i] = -2.0 * g_z_zz_z_xz[i] + 2.0 * g_zzz_zz_z_xz[i] * a_exp;

        g_z_0_0_0_zz_zz_z_yy[i] = -2.0 * g_z_zz_z_yy[i] + 2.0 * g_zzz_zz_z_yy[i] * a_exp;

        g_z_0_0_0_zz_zz_z_yz[i] = -2.0 * g_z_zz_z_yz[i] + 2.0 * g_zzz_zz_z_yz[i] * a_exp;

        g_z_0_0_0_zz_zz_z_zz[i] = -2.0 * g_z_zz_z_zz[i] + 2.0 * g_zzz_zz_z_zz[i] * a_exp;
    }
}

} // t4c_geom namespace


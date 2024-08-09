#include "GeomDeriv2000OfScalarForPDPD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom2000_pdpd_0(CSimdArray<double>& buffer_2000_pdpd,
                     const CSimdArray<double>& buffer_pdpd,
                     const CSimdArray<double>& buffer_fdpd,
                     const double a_exp) -> void
{
    const auto ndims = buffer_2000_pdpd.number_of_columns();

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

    /// Set up components of integrals buffer : buffer_2000_pdpd

    auto g_xx_0_0_0_x_xx_x_xx = buffer_2000_pdpd[0];

    auto g_xx_0_0_0_x_xx_x_xy = buffer_2000_pdpd[1];

    auto g_xx_0_0_0_x_xx_x_xz = buffer_2000_pdpd[2];

    auto g_xx_0_0_0_x_xx_x_yy = buffer_2000_pdpd[3];

    auto g_xx_0_0_0_x_xx_x_yz = buffer_2000_pdpd[4];

    auto g_xx_0_0_0_x_xx_x_zz = buffer_2000_pdpd[5];

    auto g_xx_0_0_0_x_xx_y_xx = buffer_2000_pdpd[6];

    auto g_xx_0_0_0_x_xx_y_xy = buffer_2000_pdpd[7];

    auto g_xx_0_0_0_x_xx_y_xz = buffer_2000_pdpd[8];

    auto g_xx_0_0_0_x_xx_y_yy = buffer_2000_pdpd[9];

    auto g_xx_0_0_0_x_xx_y_yz = buffer_2000_pdpd[10];

    auto g_xx_0_0_0_x_xx_y_zz = buffer_2000_pdpd[11];

    auto g_xx_0_0_0_x_xx_z_xx = buffer_2000_pdpd[12];

    auto g_xx_0_0_0_x_xx_z_xy = buffer_2000_pdpd[13];

    auto g_xx_0_0_0_x_xx_z_xz = buffer_2000_pdpd[14];

    auto g_xx_0_0_0_x_xx_z_yy = buffer_2000_pdpd[15];

    auto g_xx_0_0_0_x_xx_z_yz = buffer_2000_pdpd[16];

    auto g_xx_0_0_0_x_xx_z_zz = buffer_2000_pdpd[17];

    auto g_xx_0_0_0_x_xy_x_xx = buffer_2000_pdpd[18];

    auto g_xx_0_0_0_x_xy_x_xy = buffer_2000_pdpd[19];

    auto g_xx_0_0_0_x_xy_x_xz = buffer_2000_pdpd[20];

    auto g_xx_0_0_0_x_xy_x_yy = buffer_2000_pdpd[21];

    auto g_xx_0_0_0_x_xy_x_yz = buffer_2000_pdpd[22];

    auto g_xx_0_0_0_x_xy_x_zz = buffer_2000_pdpd[23];

    auto g_xx_0_0_0_x_xy_y_xx = buffer_2000_pdpd[24];

    auto g_xx_0_0_0_x_xy_y_xy = buffer_2000_pdpd[25];

    auto g_xx_0_0_0_x_xy_y_xz = buffer_2000_pdpd[26];

    auto g_xx_0_0_0_x_xy_y_yy = buffer_2000_pdpd[27];

    auto g_xx_0_0_0_x_xy_y_yz = buffer_2000_pdpd[28];

    auto g_xx_0_0_0_x_xy_y_zz = buffer_2000_pdpd[29];

    auto g_xx_0_0_0_x_xy_z_xx = buffer_2000_pdpd[30];

    auto g_xx_0_0_0_x_xy_z_xy = buffer_2000_pdpd[31];

    auto g_xx_0_0_0_x_xy_z_xz = buffer_2000_pdpd[32];

    auto g_xx_0_0_0_x_xy_z_yy = buffer_2000_pdpd[33];

    auto g_xx_0_0_0_x_xy_z_yz = buffer_2000_pdpd[34];

    auto g_xx_0_0_0_x_xy_z_zz = buffer_2000_pdpd[35];

    auto g_xx_0_0_0_x_xz_x_xx = buffer_2000_pdpd[36];

    auto g_xx_0_0_0_x_xz_x_xy = buffer_2000_pdpd[37];

    auto g_xx_0_0_0_x_xz_x_xz = buffer_2000_pdpd[38];

    auto g_xx_0_0_0_x_xz_x_yy = buffer_2000_pdpd[39];

    auto g_xx_0_0_0_x_xz_x_yz = buffer_2000_pdpd[40];

    auto g_xx_0_0_0_x_xz_x_zz = buffer_2000_pdpd[41];

    auto g_xx_0_0_0_x_xz_y_xx = buffer_2000_pdpd[42];

    auto g_xx_0_0_0_x_xz_y_xy = buffer_2000_pdpd[43];

    auto g_xx_0_0_0_x_xz_y_xz = buffer_2000_pdpd[44];

    auto g_xx_0_0_0_x_xz_y_yy = buffer_2000_pdpd[45];

    auto g_xx_0_0_0_x_xz_y_yz = buffer_2000_pdpd[46];

    auto g_xx_0_0_0_x_xz_y_zz = buffer_2000_pdpd[47];

    auto g_xx_0_0_0_x_xz_z_xx = buffer_2000_pdpd[48];

    auto g_xx_0_0_0_x_xz_z_xy = buffer_2000_pdpd[49];

    auto g_xx_0_0_0_x_xz_z_xz = buffer_2000_pdpd[50];

    auto g_xx_0_0_0_x_xz_z_yy = buffer_2000_pdpd[51];

    auto g_xx_0_0_0_x_xz_z_yz = buffer_2000_pdpd[52];

    auto g_xx_0_0_0_x_xz_z_zz = buffer_2000_pdpd[53];

    auto g_xx_0_0_0_x_yy_x_xx = buffer_2000_pdpd[54];

    auto g_xx_0_0_0_x_yy_x_xy = buffer_2000_pdpd[55];

    auto g_xx_0_0_0_x_yy_x_xz = buffer_2000_pdpd[56];

    auto g_xx_0_0_0_x_yy_x_yy = buffer_2000_pdpd[57];

    auto g_xx_0_0_0_x_yy_x_yz = buffer_2000_pdpd[58];

    auto g_xx_0_0_0_x_yy_x_zz = buffer_2000_pdpd[59];

    auto g_xx_0_0_0_x_yy_y_xx = buffer_2000_pdpd[60];

    auto g_xx_0_0_0_x_yy_y_xy = buffer_2000_pdpd[61];

    auto g_xx_0_0_0_x_yy_y_xz = buffer_2000_pdpd[62];

    auto g_xx_0_0_0_x_yy_y_yy = buffer_2000_pdpd[63];

    auto g_xx_0_0_0_x_yy_y_yz = buffer_2000_pdpd[64];

    auto g_xx_0_0_0_x_yy_y_zz = buffer_2000_pdpd[65];

    auto g_xx_0_0_0_x_yy_z_xx = buffer_2000_pdpd[66];

    auto g_xx_0_0_0_x_yy_z_xy = buffer_2000_pdpd[67];

    auto g_xx_0_0_0_x_yy_z_xz = buffer_2000_pdpd[68];

    auto g_xx_0_0_0_x_yy_z_yy = buffer_2000_pdpd[69];

    auto g_xx_0_0_0_x_yy_z_yz = buffer_2000_pdpd[70];

    auto g_xx_0_0_0_x_yy_z_zz = buffer_2000_pdpd[71];

    auto g_xx_0_0_0_x_yz_x_xx = buffer_2000_pdpd[72];

    auto g_xx_0_0_0_x_yz_x_xy = buffer_2000_pdpd[73];

    auto g_xx_0_0_0_x_yz_x_xz = buffer_2000_pdpd[74];

    auto g_xx_0_0_0_x_yz_x_yy = buffer_2000_pdpd[75];

    auto g_xx_0_0_0_x_yz_x_yz = buffer_2000_pdpd[76];

    auto g_xx_0_0_0_x_yz_x_zz = buffer_2000_pdpd[77];

    auto g_xx_0_0_0_x_yz_y_xx = buffer_2000_pdpd[78];

    auto g_xx_0_0_0_x_yz_y_xy = buffer_2000_pdpd[79];

    auto g_xx_0_0_0_x_yz_y_xz = buffer_2000_pdpd[80];

    auto g_xx_0_0_0_x_yz_y_yy = buffer_2000_pdpd[81];

    auto g_xx_0_0_0_x_yz_y_yz = buffer_2000_pdpd[82];

    auto g_xx_0_0_0_x_yz_y_zz = buffer_2000_pdpd[83];

    auto g_xx_0_0_0_x_yz_z_xx = buffer_2000_pdpd[84];

    auto g_xx_0_0_0_x_yz_z_xy = buffer_2000_pdpd[85];

    auto g_xx_0_0_0_x_yz_z_xz = buffer_2000_pdpd[86];

    auto g_xx_0_0_0_x_yz_z_yy = buffer_2000_pdpd[87];

    auto g_xx_0_0_0_x_yz_z_yz = buffer_2000_pdpd[88];

    auto g_xx_0_0_0_x_yz_z_zz = buffer_2000_pdpd[89];

    auto g_xx_0_0_0_x_zz_x_xx = buffer_2000_pdpd[90];

    auto g_xx_0_0_0_x_zz_x_xy = buffer_2000_pdpd[91];

    auto g_xx_0_0_0_x_zz_x_xz = buffer_2000_pdpd[92];

    auto g_xx_0_0_0_x_zz_x_yy = buffer_2000_pdpd[93];

    auto g_xx_0_0_0_x_zz_x_yz = buffer_2000_pdpd[94];

    auto g_xx_0_0_0_x_zz_x_zz = buffer_2000_pdpd[95];

    auto g_xx_0_0_0_x_zz_y_xx = buffer_2000_pdpd[96];

    auto g_xx_0_0_0_x_zz_y_xy = buffer_2000_pdpd[97];

    auto g_xx_0_0_0_x_zz_y_xz = buffer_2000_pdpd[98];

    auto g_xx_0_0_0_x_zz_y_yy = buffer_2000_pdpd[99];

    auto g_xx_0_0_0_x_zz_y_yz = buffer_2000_pdpd[100];

    auto g_xx_0_0_0_x_zz_y_zz = buffer_2000_pdpd[101];

    auto g_xx_0_0_0_x_zz_z_xx = buffer_2000_pdpd[102];

    auto g_xx_0_0_0_x_zz_z_xy = buffer_2000_pdpd[103];

    auto g_xx_0_0_0_x_zz_z_xz = buffer_2000_pdpd[104];

    auto g_xx_0_0_0_x_zz_z_yy = buffer_2000_pdpd[105];

    auto g_xx_0_0_0_x_zz_z_yz = buffer_2000_pdpd[106];

    auto g_xx_0_0_0_x_zz_z_zz = buffer_2000_pdpd[107];

    auto g_xx_0_0_0_y_xx_x_xx = buffer_2000_pdpd[108];

    auto g_xx_0_0_0_y_xx_x_xy = buffer_2000_pdpd[109];

    auto g_xx_0_0_0_y_xx_x_xz = buffer_2000_pdpd[110];

    auto g_xx_0_0_0_y_xx_x_yy = buffer_2000_pdpd[111];

    auto g_xx_0_0_0_y_xx_x_yz = buffer_2000_pdpd[112];

    auto g_xx_0_0_0_y_xx_x_zz = buffer_2000_pdpd[113];

    auto g_xx_0_0_0_y_xx_y_xx = buffer_2000_pdpd[114];

    auto g_xx_0_0_0_y_xx_y_xy = buffer_2000_pdpd[115];

    auto g_xx_0_0_0_y_xx_y_xz = buffer_2000_pdpd[116];

    auto g_xx_0_0_0_y_xx_y_yy = buffer_2000_pdpd[117];

    auto g_xx_0_0_0_y_xx_y_yz = buffer_2000_pdpd[118];

    auto g_xx_0_0_0_y_xx_y_zz = buffer_2000_pdpd[119];

    auto g_xx_0_0_0_y_xx_z_xx = buffer_2000_pdpd[120];

    auto g_xx_0_0_0_y_xx_z_xy = buffer_2000_pdpd[121];

    auto g_xx_0_0_0_y_xx_z_xz = buffer_2000_pdpd[122];

    auto g_xx_0_0_0_y_xx_z_yy = buffer_2000_pdpd[123];

    auto g_xx_0_0_0_y_xx_z_yz = buffer_2000_pdpd[124];

    auto g_xx_0_0_0_y_xx_z_zz = buffer_2000_pdpd[125];

    auto g_xx_0_0_0_y_xy_x_xx = buffer_2000_pdpd[126];

    auto g_xx_0_0_0_y_xy_x_xy = buffer_2000_pdpd[127];

    auto g_xx_0_0_0_y_xy_x_xz = buffer_2000_pdpd[128];

    auto g_xx_0_0_0_y_xy_x_yy = buffer_2000_pdpd[129];

    auto g_xx_0_0_0_y_xy_x_yz = buffer_2000_pdpd[130];

    auto g_xx_0_0_0_y_xy_x_zz = buffer_2000_pdpd[131];

    auto g_xx_0_0_0_y_xy_y_xx = buffer_2000_pdpd[132];

    auto g_xx_0_0_0_y_xy_y_xy = buffer_2000_pdpd[133];

    auto g_xx_0_0_0_y_xy_y_xz = buffer_2000_pdpd[134];

    auto g_xx_0_0_0_y_xy_y_yy = buffer_2000_pdpd[135];

    auto g_xx_0_0_0_y_xy_y_yz = buffer_2000_pdpd[136];

    auto g_xx_0_0_0_y_xy_y_zz = buffer_2000_pdpd[137];

    auto g_xx_0_0_0_y_xy_z_xx = buffer_2000_pdpd[138];

    auto g_xx_0_0_0_y_xy_z_xy = buffer_2000_pdpd[139];

    auto g_xx_0_0_0_y_xy_z_xz = buffer_2000_pdpd[140];

    auto g_xx_0_0_0_y_xy_z_yy = buffer_2000_pdpd[141];

    auto g_xx_0_0_0_y_xy_z_yz = buffer_2000_pdpd[142];

    auto g_xx_0_0_0_y_xy_z_zz = buffer_2000_pdpd[143];

    auto g_xx_0_0_0_y_xz_x_xx = buffer_2000_pdpd[144];

    auto g_xx_0_0_0_y_xz_x_xy = buffer_2000_pdpd[145];

    auto g_xx_0_0_0_y_xz_x_xz = buffer_2000_pdpd[146];

    auto g_xx_0_0_0_y_xz_x_yy = buffer_2000_pdpd[147];

    auto g_xx_0_0_0_y_xz_x_yz = buffer_2000_pdpd[148];

    auto g_xx_0_0_0_y_xz_x_zz = buffer_2000_pdpd[149];

    auto g_xx_0_0_0_y_xz_y_xx = buffer_2000_pdpd[150];

    auto g_xx_0_0_0_y_xz_y_xy = buffer_2000_pdpd[151];

    auto g_xx_0_0_0_y_xz_y_xz = buffer_2000_pdpd[152];

    auto g_xx_0_0_0_y_xz_y_yy = buffer_2000_pdpd[153];

    auto g_xx_0_0_0_y_xz_y_yz = buffer_2000_pdpd[154];

    auto g_xx_0_0_0_y_xz_y_zz = buffer_2000_pdpd[155];

    auto g_xx_0_0_0_y_xz_z_xx = buffer_2000_pdpd[156];

    auto g_xx_0_0_0_y_xz_z_xy = buffer_2000_pdpd[157];

    auto g_xx_0_0_0_y_xz_z_xz = buffer_2000_pdpd[158];

    auto g_xx_0_0_0_y_xz_z_yy = buffer_2000_pdpd[159];

    auto g_xx_0_0_0_y_xz_z_yz = buffer_2000_pdpd[160];

    auto g_xx_0_0_0_y_xz_z_zz = buffer_2000_pdpd[161];

    auto g_xx_0_0_0_y_yy_x_xx = buffer_2000_pdpd[162];

    auto g_xx_0_0_0_y_yy_x_xy = buffer_2000_pdpd[163];

    auto g_xx_0_0_0_y_yy_x_xz = buffer_2000_pdpd[164];

    auto g_xx_0_0_0_y_yy_x_yy = buffer_2000_pdpd[165];

    auto g_xx_0_0_0_y_yy_x_yz = buffer_2000_pdpd[166];

    auto g_xx_0_0_0_y_yy_x_zz = buffer_2000_pdpd[167];

    auto g_xx_0_0_0_y_yy_y_xx = buffer_2000_pdpd[168];

    auto g_xx_0_0_0_y_yy_y_xy = buffer_2000_pdpd[169];

    auto g_xx_0_0_0_y_yy_y_xz = buffer_2000_pdpd[170];

    auto g_xx_0_0_0_y_yy_y_yy = buffer_2000_pdpd[171];

    auto g_xx_0_0_0_y_yy_y_yz = buffer_2000_pdpd[172];

    auto g_xx_0_0_0_y_yy_y_zz = buffer_2000_pdpd[173];

    auto g_xx_0_0_0_y_yy_z_xx = buffer_2000_pdpd[174];

    auto g_xx_0_0_0_y_yy_z_xy = buffer_2000_pdpd[175];

    auto g_xx_0_0_0_y_yy_z_xz = buffer_2000_pdpd[176];

    auto g_xx_0_0_0_y_yy_z_yy = buffer_2000_pdpd[177];

    auto g_xx_0_0_0_y_yy_z_yz = buffer_2000_pdpd[178];

    auto g_xx_0_0_0_y_yy_z_zz = buffer_2000_pdpd[179];

    auto g_xx_0_0_0_y_yz_x_xx = buffer_2000_pdpd[180];

    auto g_xx_0_0_0_y_yz_x_xy = buffer_2000_pdpd[181];

    auto g_xx_0_0_0_y_yz_x_xz = buffer_2000_pdpd[182];

    auto g_xx_0_0_0_y_yz_x_yy = buffer_2000_pdpd[183];

    auto g_xx_0_0_0_y_yz_x_yz = buffer_2000_pdpd[184];

    auto g_xx_0_0_0_y_yz_x_zz = buffer_2000_pdpd[185];

    auto g_xx_0_0_0_y_yz_y_xx = buffer_2000_pdpd[186];

    auto g_xx_0_0_0_y_yz_y_xy = buffer_2000_pdpd[187];

    auto g_xx_0_0_0_y_yz_y_xz = buffer_2000_pdpd[188];

    auto g_xx_0_0_0_y_yz_y_yy = buffer_2000_pdpd[189];

    auto g_xx_0_0_0_y_yz_y_yz = buffer_2000_pdpd[190];

    auto g_xx_0_0_0_y_yz_y_zz = buffer_2000_pdpd[191];

    auto g_xx_0_0_0_y_yz_z_xx = buffer_2000_pdpd[192];

    auto g_xx_0_0_0_y_yz_z_xy = buffer_2000_pdpd[193];

    auto g_xx_0_0_0_y_yz_z_xz = buffer_2000_pdpd[194];

    auto g_xx_0_0_0_y_yz_z_yy = buffer_2000_pdpd[195];

    auto g_xx_0_0_0_y_yz_z_yz = buffer_2000_pdpd[196];

    auto g_xx_0_0_0_y_yz_z_zz = buffer_2000_pdpd[197];

    auto g_xx_0_0_0_y_zz_x_xx = buffer_2000_pdpd[198];

    auto g_xx_0_0_0_y_zz_x_xy = buffer_2000_pdpd[199];

    auto g_xx_0_0_0_y_zz_x_xz = buffer_2000_pdpd[200];

    auto g_xx_0_0_0_y_zz_x_yy = buffer_2000_pdpd[201];

    auto g_xx_0_0_0_y_zz_x_yz = buffer_2000_pdpd[202];

    auto g_xx_0_0_0_y_zz_x_zz = buffer_2000_pdpd[203];

    auto g_xx_0_0_0_y_zz_y_xx = buffer_2000_pdpd[204];

    auto g_xx_0_0_0_y_zz_y_xy = buffer_2000_pdpd[205];

    auto g_xx_0_0_0_y_zz_y_xz = buffer_2000_pdpd[206];

    auto g_xx_0_0_0_y_zz_y_yy = buffer_2000_pdpd[207];

    auto g_xx_0_0_0_y_zz_y_yz = buffer_2000_pdpd[208];

    auto g_xx_0_0_0_y_zz_y_zz = buffer_2000_pdpd[209];

    auto g_xx_0_0_0_y_zz_z_xx = buffer_2000_pdpd[210];

    auto g_xx_0_0_0_y_zz_z_xy = buffer_2000_pdpd[211];

    auto g_xx_0_0_0_y_zz_z_xz = buffer_2000_pdpd[212];

    auto g_xx_0_0_0_y_zz_z_yy = buffer_2000_pdpd[213];

    auto g_xx_0_0_0_y_zz_z_yz = buffer_2000_pdpd[214];

    auto g_xx_0_0_0_y_zz_z_zz = buffer_2000_pdpd[215];

    auto g_xx_0_0_0_z_xx_x_xx = buffer_2000_pdpd[216];

    auto g_xx_0_0_0_z_xx_x_xy = buffer_2000_pdpd[217];

    auto g_xx_0_0_0_z_xx_x_xz = buffer_2000_pdpd[218];

    auto g_xx_0_0_0_z_xx_x_yy = buffer_2000_pdpd[219];

    auto g_xx_0_0_0_z_xx_x_yz = buffer_2000_pdpd[220];

    auto g_xx_0_0_0_z_xx_x_zz = buffer_2000_pdpd[221];

    auto g_xx_0_0_0_z_xx_y_xx = buffer_2000_pdpd[222];

    auto g_xx_0_0_0_z_xx_y_xy = buffer_2000_pdpd[223];

    auto g_xx_0_0_0_z_xx_y_xz = buffer_2000_pdpd[224];

    auto g_xx_0_0_0_z_xx_y_yy = buffer_2000_pdpd[225];

    auto g_xx_0_0_0_z_xx_y_yz = buffer_2000_pdpd[226];

    auto g_xx_0_0_0_z_xx_y_zz = buffer_2000_pdpd[227];

    auto g_xx_0_0_0_z_xx_z_xx = buffer_2000_pdpd[228];

    auto g_xx_0_0_0_z_xx_z_xy = buffer_2000_pdpd[229];

    auto g_xx_0_0_0_z_xx_z_xz = buffer_2000_pdpd[230];

    auto g_xx_0_0_0_z_xx_z_yy = buffer_2000_pdpd[231];

    auto g_xx_0_0_0_z_xx_z_yz = buffer_2000_pdpd[232];

    auto g_xx_0_0_0_z_xx_z_zz = buffer_2000_pdpd[233];

    auto g_xx_0_0_0_z_xy_x_xx = buffer_2000_pdpd[234];

    auto g_xx_0_0_0_z_xy_x_xy = buffer_2000_pdpd[235];

    auto g_xx_0_0_0_z_xy_x_xz = buffer_2000_pdpd[236];

    auto g_xx_0_0_0_z_xy_x_yy = buffer_2000_pdpd[237];

    auto g_xx_0_0_0_z_xy_x_yz = buffer_2000_pdpd[238];

    auto g_xx_0_0_0_z_xy_x_zz = buffer_2000_pdpd[239];

    auto g_xx_0_0_0_z_xy_y_xx = buffer_2000_pdpd[240];

    auto g_xx_0_0_0_z_xy_y_xy = buffer_2000_pdpd[241];

    auto g_xx_0_0_0_z_xy_y_xz = buffer_2000_pdpd[242];

    auto g_xx_0_0_0_z_xy_y_yy = buffer_2000_pdpd[243];

    auto g_xx_0_0_0_z_xy_y_yz = buffer_2000_pdpd[244];

    auto g_xx_0_0_0_z_xy_y_zz = buffer_2000_pdpd[245];

    auto g_xx_0_0_0_z_xy_z_xx = buffer_2000_pdpd[246];

    auto g_xx_0_0_0_z_xy_z_xy = buffer_2000_pdpd[247];

    auto g_xx_0_0_0_z_xy_z_xz = buffer_2000_pdpd[248];

    auto g_xx_0_0_0_z_xy_z_yy = buffer_2000_pdpd[249];

    auto g_xx_0_0_0_z_xy_z_yz = buffer_2000_pdpd[250];

    auto g_xx_0_0_0_z_xy_z_zz = buffer_2000_pdpd[251];

    auto g_xx_0_0_0_z_xz_x_xx = buffer_2000_pdpd[252];

    auto g_xx_0_0_0_z_xz_x_xy = buffer_2000_pdpd[253];

    auto g_xx_0_0_0_z_xz_x_xz = buffer_2000_pdpd[254];

    auto g_xx_0_0_0_z_xz_x_yy = buffer_2000_pdpd[255];

    auto g_xx_0_0_0_z_xz_x_yz = buffer_2000_pdpd[256];

    auto g_xx_0_0_0_z_xz_x_zz = buffer_2000_pdpd[257];

    auto g_xx_0_0_0_z_xz_y_xx = buffer_2000_pdpd[258];

    auto g_xx_0_0_0_z_xz_y_xy = buffer_2000_pdpd[259];

    auto g_xx_0_0_0_z_xz_y_xz = buffer_2000_pdpd[260];

    auto g_xx_0_0_0_z_xz_y_yy = buffer_2000_pdpd[261];

    auto g_xx_0_0_0_z_xz_y_yz = buffer_2000_pdpd[262];

    auto g_xx_0_0_0_z_xz_y_zz = buffer_2000_pdpd[263];

    auto g_xx_0_0_0_z_xz_z_xx = buffer_2000_pdpd[264];

    auto g_xx_0_0_0_z_xz_z_xy = buffer_2000_pdpd[265];

    auto g_xx_0_0_0_z_xz_z_xz = buffer_2000_pdpd[266];

    auto g_xx_0_0_0_z_xz_z_yy = buffer_2000_pdpd[267];

    auto g_xx_0_0_0_z_xz_z_yz = buffer_2000_pdpd[268];

    auto g_xx_0_0_0_z_xz_z_zz = buffer_2000_pdpd[269];

    auto g_xx_0_0_0_z_yy_x_xx = buffer_2000_pdpd[270];

    auto g_xx_0_0_0_z_yy_x_xy = buffer_2000_pdpd[271];

    auto g_xx_0_0_0_z_yy_x_xz = buffer_2000_pdpd[272];

    auto g_xx_0_0_0_z_yy_x_yy = buffer_2000_pdpd[273];

    auto g_xx_0_0_0_z_yy_x_yz = buffer_2000_pdpd[274];

    auto g_xx_0_0_0_z_yy_x_zz = buffer_2000_pdpd[275];

    auto g_xx_0_0_0_z_yy_y_xx = buffer_2000_pdpd[276];

    auto g_xx_0_0_0_z_yy_y_xy = buffer_2000_pdpd[277];

    auto g_xx_0_0_0_z_yy_y_xz = buffer_2000_pdpd[278];

    auto g_xx_0_0_0_z_yy_y_yy = buffer_2000_pdpd[279];

    auto g_xx_0_0_0_z_yy_y_yz = buffer_2000_pdpd[280];

    auto g_xx_0_0_0_z_yy_y_zz = buffer_2000_pdpd[281];

    auto g_xx_0_0_0_z_yy_z_xx = buffer_2000_pdpd[282];

    auto g_xx_0_0_0_z_yy_z_xy = buffer_2000_pdpd[283];

    auto g_xx_0_0_0_z_yy_z_xz = buffer_2000_pdpd[284];

    auto g_xx_0_0_0_z_yy_z_yy = buffer_2000_pdpd[285];

    auto g_xx_0_0_0_z_yy_z_yz = buffer_2000_pdpd[286];

    auto g_xx_0_0_0_z_yy_z_zz = buffer_2000_pdpd[287];

    auto g_xx_0_0_0_z_yz_x_xx = buffer_2000_pdpd[288];

    auto g_xx_0_0_0_z_yz_x_xy = buffer_2000_pdpd[289];

    auto g_xx_0_0_0_z_yz_x_xz = buffer_2000_pdpd[290];

    auto g_xx_0_0_0_z_yz_x_yy = buffer_2000_pdpd[291];

    auto g_xx_0_0_0_z_yz_x_yz = buffer_2000_pdpd[292];

    auto g_xx_0_0_0_z_yz_x_zz = buffer_2000_pdpd[293];

    auto g_xx_0_0_0_z_yz_y_xx = buffer_2000_pdpd[294];

    auto g_xx_0_0_0_z_yz_y_xy = buffer_2000_pdpd[295];

    auto g_xx_0_0_0_z_yz_y_xz = buffer_2000_pdpd[296];

    auto g_xx_0_0_0_z_yz_y_yy = buffer_2000_pdpd[297];

    auto g_xx_0_0_0_z_yz_y_yz = buffer_2000_pdpd[298];

    auto g_xx_0_0_0_z_yz_y_zz = buffer_2000_pdpd[299];

    auto g_xx_0_0_0_z_yz_z_xx = buffer_2000_pdpd[300];

    auto g_xx_0_0_0_z_yz_z_xy = buffer_2000_pdpd[301];

    auto g_xx_0_0_0_z_yz_z_xz = buffer_2000_pdpd[302];

    auto g_xx_0_0_0_z_yz_z_yy = buffer_2000_pdpd[303];

    auto g_xx_0_0_0_z_yz_z_yz = buffer_2000_pdpd[304];

    auto g_xx_0_0_0_z_yz_z_zz = buffer_2000_pdpd[305];

    auto g_xx_0_0_0_z_zz_x_xx = buffer_2000_pdpd[306];

    auto g_xx_0_0_0_z_zz_x_xy = buffer_2000_pdpd[307];

    auto g_xx_0_0_0_z_zz_x_xz = buffer_2000_pdpd[308];

    auto g_xx_0_0_0_z_zz_x_yy = buffer_2000_pdpd[309];

    auto g_xx_0_0_0_z_zz_x_yz = buffer_2000_pdpd[310];

    auto g_xx_0_0_0_z_zz_x_zz = buffer_2000_pdpd[311];

    auto g_xx_0_0_0_z_zz_y_xx = buffer_2000_pdpd[312];

    auto g_xx_0_0_0_z_zz_y_xy = buffer_2000_pdpd[313];

    auto g_xx_0_0_0_z_zz_y_xz = buffer_2000_pdpd[314];

    auto g_xx_0_0_0_z_zz_y_yy = buffer_2000_pdpd[315];

    auto g_xx_0_0_0_z_zz_y_yz = buffer_2000_pdpd[316];

    auto g_xx_0_0_0_z_zz_y_zz = buffer_2000_pdpd[317];

    auto g_xx_0_0_0_z_zz_z_xx = buffer_2000_pdpd[318];

    auto g_xx_0_0_0_z_zz_z_xy = buffer_2000_pdpd[319];

    auto g_xx_0_0_0_z_zz_z_xz = buffer_2000_pdpd[320];

    auto g_xx_0_0_0_z_zz_z_yy = buffer_2000_pdpd[321];

    auto g_xx_0_0_0_z_zz_z_yz = buffer_2000_pdpd[322];

    auto g_xx_0_0_0_z_zz_z_zz = buffer_2000_pdpd[323];

    auto g_xy_0_0_0_x_xx_x_xx = buffer_2000_pdpd[324];

    auto g_xy_0_0_0_x_xx_x_xy = buffer_2000_pdpd[325];

    auto g_xy_0_0_0_x_xx_x_xz = buffer_2000_pdpd[326];

    auto g_xy_0_0_0_x_xx_x_yy = buffer_2000_pdpd[327];

    auto g_xy_0_0_0_x_xx_x_yz = buffer_2000_pdpd[328];

    auto g_xy_0_0_0_x_xx_x_zz = buffer_2000_pdpd[329];

    auto g_xy_0_0_0_x_xx_y_xx = buffer_2000_pdpd[330];

    auto g_xy_0_0_0_x_xx_y_xy = buffer_2000_pdpd[331];

    auto g_xy_0_0_0_x_xx_y_xz = buffer_2000_pdpd[332];

    auto g_xy_0_0_0_x_xx_y_yy = buffer_2000_pdpd[333];

    auto g_xy_0_0_0_x_xx_y_yz = buffer_2000_pdpd[334];

    auto g_xy_0_0_0_x_xx_y_zz = buffer_2000_pdpd[335];

    auto g_xy_0_0_0_x_xx_z_xx = buffer_2000_pdpd[336];

    auto g_xy_0_0_0_x_xx_z_xy = buffer_2000_pdpd[337];

    auto g_xy_0_0_0_x_xx_z_xz = buffer_2000_pdpd[338];

    auto g_xy_0_0_0_x_xx_z_yy = buffer_2000_pdpd[339];

    auto g_xy_0_0_0_x_xx_z_yz = buffer_2000_pdpd[340];

    auto g_xy_0_0_0_x_xx_z_zz = buffer_2000_pdpd[341];

    auto g_xy_0_0_0_x_xy_x_xx = buffer_2000_pdpd[342];

    auto g_xy_0_0_0_x_xy_x_xy = buffer_2000_pdpd[343];

    auto g_xy_0_0_0_x_xy_x_xz = buffer_2000_pdpd[344];

    auto g_xy_0_0_0_x_xy_x_yy = buffer_2000_pdpd[345];

    auto g_xy_0_0_0_x_xy_x_yz = buffer_2000_pdpd[346];

    auto g_xy_0_0_0_x_xy_x_zz = buffer_2000_pdpd[347];

    auto g_xy_0_0_0_x_xy_y_xx = buffer_2000_pdpd[348];

    auto g_xy_0_0_0_x_xy_y_xy = buffer_2000_pdpd[349];

    auto g_xy_0_0_0_x_xy_y_xz = buffer_2000_pdpd[350];

    auto g_xy_0_0_0_x_xy_y_yy = buffer_2000_pdpd[351];

    auto g_xy_0_0_0_x_xy_y_yz = buffer_2000_pdpd[352];

    auto g_xy_0_0_0_x_xy_y_zz = buffer_2000_pdpd[353];

    auto g_xy_0_0_0_x_xy_z_xx = buffer_2000_pdpd[354];

    auto g_xy_0_0_0_x_xy_z_xy = buffer_2000_pdpd[355];

    auto g_xy_0_0_0_x_xy_z_xz = buffer_2000_pdpd[356];

    auto g_xy_0_0_0_x_xy_z_yy = buffer_2000_pdpd[357];

    auto g_xy_0_0_0_x_xy_z_yz = buffer_2000_pdpd[358];

    auto g_xy_0_0_0_x_xy_z_zz = buffer_2000_pdpd[359];

    auto g_xy_0_0_0_x_xz_x_xx = buffer_2000_pdpd[360];

    auto g_xy_0_0_0_x_xz_x_xy = buffer_2000_pdpd[361];

    auto g_xy_0_0_0_x_xz_x_xz = buffer_2000_pdpd[362];

    auto g_xy_0_0_0_x_xz_x_yy = buffer_2000_pdpd[363];

    auto g_xy_0_0_0_x_xz_x_yz = buffer_2000_pdpd[364];

    auto g_xy_0_0_0_x_xz_x_zz = buffer_2000_pdpd[365];

    auto g_xy_0_0_0_x_xz_y_xx = buffer_2000_pdpd[366];

    auto g_xy_0_0_0_x_xz_y_xy = buffer_2000_pdpd[367];

    auto g_xy_0_0_0_x_xz_y_xz = buffer_2000_pdpd[368];

    auto g_xy_0_0_0_x_xz_y_yy = buffer_2000_pdpd[369];

    auto g_xy_0_0_0_x_xz_y_yz = buffer_2000_pdpd[370];

    auto g_xy_0_0_0_x_xz_y_zz = buffer_2000_pdpd[371];

    auto g_xy_0_0_0_x_xz_z_xx = buffer_2000_pdpd[372];

    auto g_xy_0_0_0_x_xz_z_xy = buffer_2000_pdpd[373];

    auto g_xy_0_0_0_x_xz_z_xz = buffer_2000_pdpd[374];

    auto g_xy_0_0_0_x_xz_z_yy = buffer_2000_pdpd[375];

    auto g_xy_0_0_0_x_xz_z_yz = buffer_2000_pdpd[376];

    auto g_xy_0_0_0_x_xz_z_zz = buffer_2000_pdpd[377];

    auto g_xy_0_0_0_x_yy_x_xx = buffer_2000_pdpd[378];

    auto g_xy_0_0_0_x_yy_x_xy = buffer_2000_pdpd[379];

    auto g_xy_0_0_0_x_yy_x_xz = buffer_2000_pdpd[380];

    auto g_xy_0_0_0_x_yy_x_yy = buffer_2000_pdpd[381];

    auto g_xy_0_0_0_x_yy_x_yz = buffer_2000_pdpd[382];

    auto g_xy_0_0_0_x_yy_x_zz = buffer_2000_pdpd[383];

    auto g_xy_0_0_0_x_yy_y_xx = buffer_2000_pdpd[384];

    auto g_xy_0_0_0_x_yy_y_xy = buffer_2000_pdpd[385];

    auto g_xy_0_0_0_x_yy_y_xz = buffer_2000_pdpd[386];

    auto g_xy_0_0_0_x_yy_y_yy = buffer_2000_pdpd[387];

    auto g_xy_0_0_0_x_yy_y_yz = buffer_2000_pdpd[388];

    auto g_xy_0_0_0_x_yy_y_zz = buffer_2000_pdpd[389];

    auto g_xy_0_0_0_x_yy_z_xx = buffer_2000_pdpd[390];

    auto g_xy_0_0_0_x_yy_z_xy = buffer_2000_pdpd[391];

    auto g_xy_0_0_0_x_yy_z_xz = buffer_2000_pdpd[392];

    auto g_xy_0_0_0_x_yy_z_yy = buffer_2000_pdpd[393];

    auto g_xy_0_0_0_x_yy_z_yz = buffer_2000_pdpd[394];

    auto g_xy_0_0_0_x_yy_z_zz = buffer_2000_pdpd[395];

    auto g_xy_0_0_0_x_yz_x_xx = buffer_2000_pdpd[396];

    auto g_xy_0_0_0_x_yz_x_xy = buffer_2000_pdpd[397];

    auto g_xy_0_0_0_x_yz_x_xz = buffer_2000_pdpd[398];

    auto g_xy_0_0_0_x_yz_x_yy = buffer_2000_pdpd[399];

    auto g_xy_0_0_0_x_yz_x_yz = buffer_2000_pdpd[400];

    auto g_xy_0_0_0_x_yz_x_zz = buffer_2000_pdpd[401];

    auto g_xy_0_0_0_x_yz_y_xx = buffer_2000_pdpd[402];

    auto g_xy_0_0_0_x_yz_y_xy = buffer_2000_pdpd[403];

    auto g_xy_0_0_0_x_yz_y_xz = buffer_2000_pdpd[404];

    auto g_xy_0_0_0_x_yz_y_yy = buffer_2000_pdpd[405];

    auto g_xy_0_0_0_x_yz_y_yz = buffer_2000_pdpd[406];

    auto g_xy_0_0_0_x_yz_y_zz = buffer_2000_pdpd[407];

    auto g_xy_0_0_0_x_yz_z_xx = buffer_2000_pdpd[408];

    auto g_xy_0_0_0_x_yz_z_xy = buffer_2000_pdpd[409];

    auto g_xy_0_0_0_x_yz_z_xz = buffer_2000_pdpd[410];

    auto g_xy_0_0_0_x_yz_z_yy = buffer_2000_pdpd[411];

    auto g_xy_0_0_0_x_yz_z_yz = buffer_2000_pdpd[412];

    auto g_xy_0_0_0_x_yz_z_zz = buffer_2000_pdpd[413];

    auto g_xy_0_0_0_x_zz_x_xx = buffer_2000_pdpd[414];

    auto g_xy_0_0_0_x_zz_x_xy = buffer_2000_pdpd[415];

    auto g_xy_0_0_0_x_zz_x_xz = buffer_2000_pdpd[416];

    auto g_xy_0_0_0_x_zz_x_yy = buffer_2000_pdpd[417];

    auto g_xy_0_0_0_x_zz_x_yz = buffer_2000_pdpd[418];

    auto g_xy_0_0_0_x_zz_x_zz = buffer_2000_pdpd[419];

    auto g_xy_0_0_0_x_zz_y_xx = buffer_2000_pdpd[420];

    auto g_xy_0_0_0_x_zz_y_xy = buffer_2000_pdpd[421];

    auto g_xy_0_0_0_x_zz_y_xz = buffer_2000_pdpd[422];

    auto g_xy_0_0_0_x_zz_y_yy = buffer_2000_pdpd[423];

    auto g_xy_0_0_0_x_zz_y_yz = buffer_2000_pdpd[424];

    auto g_xy_0_0_0_x_zz_y_zz = buffer_2000_pdpd[425];

    auto g_xy_0_0_0_x_zz_z_xx = buffer_2000_pdpd[426];

    auto g_xy_0_0_0_x_zz_z_xy = buffer_2000_pdpd[427];

    auto g_xy_0_0_0_x_zz_z_xz = buffer_2000_pdpd[428];

    auto g_xy_0_0_0_x_zz_z_yy = buffer_2000_pdpd[429];

    auto g_xy_0_0_0_x_zz_z_yz = buffer_2000_pdpd[430];

    auto g_xy_0_0_0_x_zz_z_zz = buffer_2000_pdpd[431];

    auto g_xy_0_0_0_y_xx_x_xx = buffer_2000_pdpd[432];

    auto g_xy_0_0_0_y_xx_x_xy = buffer_2000_pdpd[433];

    auto g_xy_0_0_0_y_xx_x_xz = buffer_2000_pdpd[434];

    auto g_xy_0_0_0_y_xx_x_yy = buffer_2000_pdpd[435];

    auto g_xy_0_0_0_y_xx_x_yz = buffer_2000_pdpd[436];

    auto g_xy_0_0_0_y_xx_x_zz = buffer_2000_pdpd[437];

    auto g_xy_0_0_0_y_xx_y_xx = buffer_2000_pdpd[438];

    auto g_xy_0_0_0_y_xx_y_xy = buffer_2000_pdpd[439];

    auto g_xy_0_0_0_y_xx_y_xz = buffer_2000_pdpd[440];

    auto g_xy_0_0_0_y_xx_y_yy = buffer_2000_pdpd[441];

    auto g_xy_0_0_0_y_xx_y_yz = buffer_2000_pdpd[442];

    auto g_xy_0_0_0_y_xx_y_zz = buffer_2000_pdpd[443];

    auto g_xy_0_0_0_y_xx_z_xx = buffer_2000_pdpd[444];

    auto g_xy_0_0_0_y_xx_z_xy = buffer_2000_pdpd[445];

    auto g_xy_0_0_0_y_xx_z_xz = buffer_2000_pdpd[446];

    auto g_xy_0_0_0_y_xx_z_yy = buffer_2000_pdpd[447];

    auto g_xy_0_0_0_y_xx_z_yz = buffer_2000_pdpd[448];

    auto g_xy_0_0_0_y_xx_z_zz = buffer_2000_pdpd[449];

    auto g_xy_0_0_0_y_xy_x_xx = buffer_2000_pdpd[450];

    auto g_xy_0_0_0_y_xy_x_xy = buffer_2000_pdpd[451];

    auto g_xy_0_0_0_y_xy_x_xz = buffer_2000_pdpd[452];

    auto g_xy_0_0_0_y_xy_x_yy = buffer_2000_pdpd[453];

    auto g_xy_0_0_0_y_xy_x_yz = buffer_2000_pdpd[454];

    auto g_xy_0_0_0_y_xy_x_zz = buffer_2000_pdpd[455];

    auto g_xy_0_0_0_y_xy_y_xx = buffer_2000_pdpd[456];

    auto g_xy_0_0_0_y_xy_y_xy = buffer_2000_pdpd[457];

    auto g_xy_0_0_0_y_xy_y_xz = buffer_2000_pdpd[458];

    auto g_xy_0_0_0_y_xy_y_yy = buffer_2000_pdpd[459];

    auto g_xy_0_0_0_y_xy_y_yz = buffer_2000_pdpd[460];

    auto g_xy_0_0_0_y_xy_y_zz = buffer_2000_pdpd[461];

    auto g_xy_0_0_0_y_xy_z_xx = buffer_2000_pdpd[462];

    auto g_xy_0_0_0_y_xy_z_xy = buffer_2000_pdpd[463];

    auto g_xy_0_0_0_y_xy_z_xz = buffer_2000_pdpd[464];

    auto g_xy_0_0_0_y_xy_z_yy = buffer_2000_pdpd[465];

    auto g_xy_0_0_0_y_xy_z_yz = buffer_2000_pdpd[466];

    auto g_xy_0_0_0_y_xy_z_zz = buffer_2000_pdpd[467];

    auto g_xy_0_0_0_y_xz_x_xx = buffer_2000_pdpd[468];

    auto g_xy_0_0_0_y_xz_x_xy = buffer_2000_pdpd[469];

    auto g_xy_0_0_0_y_xz_x_xz = buffer_2000_pdpd[470];

    auto g_xy_0_0_0_y_xz_x_yy = buffer_2000_pdpd[471];

    auto g_xy_0_0_0_y_xz_x_yz = buffer_2000_pdpd[472];

    auto g_xy_0_0_0_y_xz_x_zz = buffer_2000_pdpd[473];

    auto g_xy_0_0_0_y_xz_y_xx = buffer_2000_pdpd[474];

    auto g_xy_0_0_0_y_xz_y_xy = buffer_2000_pdpd[475];

    auto g_xy_0_0_0_y_xz_y_xz = buffer_2000_pdpd[476];

    auto g_xy_0_0_0_y_xz_y_yy = buffer_2000_pdpd[477];

    auto g_xy_0_0_0_y_xz_y_yz = buffer_2000_pdpd[478];

    auto g_xy_0_0_0_y_xz_y_zz = buffer_2000_pdpd[479];

    auto g_xy_0_0_0_y_xz_z_xx = buffer_2000_pdpd[480];

    auto g_xy_0_0_0_y_xz_z_xy = buffer_2000_pdpd[481];

    auto g_xy_0_0_0_y_xz_z_xz = buffer_2000_pdpd[482];

    auto g_xy_0_0_0_y_xz_z_yy = buffer_2000_pdpd[483];

    auto g_xy_0_0_0_y_xz_z_yz = buffer_2000_pdpd[484];

    auto g_xy_0_0_0_y_xz_z_zz = buffer_2000_pdpd[485];

    auto g_xy_0_0_0_y_yy_x_xx = buffer_2000_pdpd[486];

    auto g_xy_0_0_0_y_yy_x_xy = buffer_2000_pdpd[487];

    auto g_xy_0_0_0_y_yy_x_xz = buffer_2000_pdpd[488];

    auto g_xy_0_0_0_y_yy_x_yy = buffer_2000_pdpd[489];

    auto g_xy_0_0_0_y_yy_x_yz = buffer_2000_pdpd[490];

    auto g_xy_0_0_0_y_yy_x_zz = buffer_2000_pdpd[491];

    auto g_xy_0_0_0_y_yy_y_xx = buffer_2000_pdpd[492];

    auto g_xy_0_0_0_y_yy_y_xy = buffer_2000_pdpd[493];

    auto g_xy_0_0_0_y_yy_y_xz = buffer_2000_pdpd[494];

    auto g_xy_0_0_0_y_yy_y_yy = buffer_2000_pdpd[495];

    auto g_xy_0_0_0_y_yy_y_yz = buffer_2000_pdpd[496];

    auto g_xy_0_0_0_y_yy_y_zz = buffer_2000_pdpd[497];

    auto g_xy_0_0_0_y_yy_z_xx = buffer_2000_pdpd[498];

    auto g_xy_0_0_0_y_yy_z_xy = buffer_2000_pdpd[499];

    auto g_xy_0_0_0_y_yy_z_xz = buffer_2000_pdpd[500];

    auto g_xy_0_0_0_y_yy_z_yy = buffer_2000_pdpd[501];

    auto g_xy_0_0_0_y_yy_z_yz = buffer_2000_pdpd[502];

    auto g_xy_0_0_0_y_yy_z_zz = buffer_2000_pdpd[503];

    auto g_xy_0_0_0_y_yz_x_xx = buffer_2000_pdpd[504];

    auto g_xy_0_0_0_y_yz_x_xy = buffer_2000_pdpd[505];

    auto g_xy_0_0_0_y_yz_x_xz = buffer_2000_pdpd[506];

    auto g_xy_0_0_0_y_yz_x_yy = buffer_2000_pdpd[507];

    auto g_xy_0_0_0_y_yz_x_yz = buffer_2000_pdpd[508];

    auto g_xy_0_0_0_y_yz_x_zz = buffer_2000_pdpd[509];

    auto g_xy_0_0_0_y_yz_y_xx = buffer_2000_pdpd[510];

    auto g_xy_0_0_0_y_yz_y_xy = buffer_2000_pdpd[511];

    auto g_xy_0_0_0_y_yz_y_xz = buffer_2000_pdpd[512];

    auto g_xy_0_0_0_y_yz_y_yy = buffer_2000_pdpd[513];

    auto g_xy_0_0_0_y_yz_y_yz = buffer_2000_pdpd[514];

    auto g_xy_0_0_0_y_yz_y_zz = buffer_2000_pdpd[515];

    auto g_xy_0_0_0_y_yz_z_xx = buffer_2000_pdpd[516];

    auto g_xy_0_0_0_y_yz_z_xy = buffer_2000_pdpd[517];

    auto g_xy_0_0_0_y_yz_z_xz = buffer_2000_pdpd[518];

    auto g_xy_0_0_0_y_yz_z_yy = buffer_2000_pdpd[519];

    auto g_xy_0_0_0_y_yz_z_yz = buffer_2000_pdpd[520];

    auto g_xy_0_0_0_y_yz_z_zz = buffer_2000_pdpd[521];

    auto g_xy_0_0_0_y_zz_x_xx = buffer_2000_pdpd[522];

    auto g_xy_0_0_0_y_zz_x_xy = buffer_2000_pdpd[523];

    auto g_xy_0_0_0_y_zz_x_xz = buffer_2000_pdpd[524];

    auto g_xy_0_0_0_y_zz_x_yy = buffer_2000_pdpd[525];

    auto g_xy_0_0_0_y_zz_x_yz = buffer_2000_pdpd[526];

    auto g_xy_0_0_0_y_zz_x_zz = buffer_2000_pdpd[527];

    auto g_xy_0_0_0_y_zz_y_xx = buffer_2000_pdpd[528];

    auto g_xy_0_0_0_y_zz_y_xy = buffer_2000_pdpd[529];

    auto g_xy_0_0_0_y_zz_y_xz = buffer_2000_pdpd[530];

    auto g_xy_0_0_0_y_zz_y_yy = buffer_2000_pdpd[531];

    auto g_xy_0_0_0_y_zz_y_yz = buffer_2000_pdpd[532];

    auto g_xy_0_0_0_y_zz_y_zz = buffer_2000_pdpd[533];

    auto g_xy_0_0_0_y_zz_z_xx = buffer_2000_pdpd[534];

    auto g_xy_0_0_0_y_zz_z_xy = buffer_2000_pdpd[535];

    auto g_xy_0_0_0_y_zz_z_xz = buffer_2000_pdpd[536];

    auto g_xy_0_0_0_y_zz_z_yy = buffer_2000_pdpd[537];

    auto g_xy_0_0_0_y_zz_z_yz = buffer_2000_pdpd[538];

    auto g_xy_0_0_0_y_zz_z_zz = buffer_2000_pdpd[539];

    auto g_xy_0_0_0_z_xx_x_xx = buffer_2000_pdpd[540];

    auto g_xy_0_0_0_z_xx_x_xy = buffer_2000_pdpd[541];

    auto g_xy_0_0_0_z_xx_x_xz = buffer_2000_pdpd[542];

    auto g_xy_0_0_0_z_xx_x_yy = buffer_2000_pdpd[543];

    auto g_xy_0_0_0_z_xx_x_yz = buffer_2000_pdpd[544];

    auto g_xy_0_0_0_z_xx_x_zz = buffer_2000_pdpd[545];

    auto g_xy_0_0_0_z_xx_y_xx = buffer_2000_pdpd[546];

    auto g_xy_0_0_0_z_xx_y_xy = buffer_2000_pdpd[547];

    auto g_xy_0_0_0_z_xx_y_xz = buffer_2000_pdpd[548];

    auto g_xy_0_0_0_z_xx_y_yy = buffer_2000_pdpd[549];

    auto g_xy_0_0_0_z_xx_y_yz = buffer_2000_pdpd[550];

    auto g_xy_0_0_0_z_xx_y_zz = buffer_2000_pdpd[551];

    auto g_xy_0_0_0_z_xx_z_xx = buffer_2000_pdpd[552];

    auto g_xy_0_0_0_z_xx_z_xy = buffer_2000_pdpd[553];

    auto g_xy_0_0_0_z_xx_z_xz = buffer_2000_pdpd[554];

    auto g_xy_0_0_0_z_xx_z_yy = buffer_2000_pdpd[555];

    auto g_xy_0_0_0_z_xx_z_yz = buffer_2000_pdpd[556];

    auto g_xy_0_0_0_z_xx_z_zz = buffer_2000_pdpd[557];

    auto g_xy_0_0_0_z_xy_x_xx = buffer_2000_pdpd[558];

    auto g_xy_0_0_0_z_xy_x_xy = buffer_2000_pdpd[559];

    auto g_xy_0_0_0_z_xy_x_xz = buffer_2000_pdpd[560];

    auto g_xy_0_0_0_z_xy_x_yy = buffer_2000_pdpd[561];

    auto g_xy_0_0_0_z_xy_x_yz = buffer_2000_pdpd[562];

    auto g_xy_0_0_0_z_xy_x_zz = buffer_2000_pdpd[563];

    auto g_xy_0_0_0_z_xy_y_xx = buffer_2000_pdpd[564];

    auto g_xy_0_0_0_z_xy_y_xy = buffer_2000_pdpd[565];

    auto g_xy_0_0_0_z_xy_y_xz = buffer_2000_pdpd[566];

    auto g_xy_0_0_0_z_xy_y_yy = buffer_2000_pdpd[567];

    auto g_xy_0_0_0_z_xy_y_yz = buffer_2000_pdpd[568];

    auto g_xy_0_0_0_z_xy_y_zz = buffer_2000_pdpd[569];

    auto g_xy_0_0_0_z_xy_z_xx = buffer_2000_pdpd[570];

    auto g_xy_0_0_0_z_xy_z_xy = buffer_2000_pdpd[571];

    auto g_xy_0_0_0_z_xy_z_xz = buffer_2000_pdpd[572];

    auto g_xy_0_0_0_z_xy_z_yy = buffer_2000_pdpd[573];

    auto g_xy_0_0_0_z_xy_z_yz = buffer_2000_pdpd[574];

    auto g_xy_0_0_0_z_xy_z_zz = buffer_2000_pdpd[575];

    auto g_xy_0_0_0_z_xz_x_xx = buffer_2000_pdpd[576];

    auto g_xy_0_0_0_z_xz_x_xy = buffer_2000_pdpd[577];

    auto g_xy_0_0_0_z_xz_x_xz = buffer_2000_pdpd[578];

    auto g_xy_0_0_0_z_xz_x_yy = buffer_2000_pdpd[579];

    auto g_xy_0_0_0_z_xz_x_yz = buffer_2000_pdpd[580];

    auto g_xy_0_0_0_z_xz_x_zz = buffer_2000_pdpd[581];

    auto g_xy_0_0_0_z_xz_y_xx = buffer_2000_pdpd[582];

    auto g_xy_0_0_0_z_xz_y_xy = buffer_2000_pdpd[583];

    auto g_xy_0_0_0_z_xz_y_xz = buffer_2000_pdpd[584];

    auto g_xy_0_0_0_z_xz_y_yy = buffer_2000_pdpd[585];

    auto g_xy_0_0_0_z_xz_y_yz = buffer_2000_pdpd[586];

    auto g_xy_0_0_0_z_xz_y_zz = buffer_2000_pdpd[587];

    auto g_xy_0_0_0_z_xz_z_xx = buffer_2000_pdpd[588];

    auto g_xy_0_0_0_z_xz_z_xy = buffer_2000_pdpd[589];

    auto g_xy_0_0_0_z_xz_z_xz = buffer_2000_pdpd[590];

    auto g_xy_0_0_0_z_xz_z_yy = buffer_2000_pdpd[591];

    auto g_xy_0_0_0_z_xz_z_yz = buffer_2000_pdpd[592];

    auto g_xy_0_0_0_z_xz_z_zz = buffer_2000_pdpd[593];

    auto g_xy_0_0_0_z_yy_x_xx = buffer_2000_pdpd[594];

    auto g_xy_0_0_0_z_yy_x_xy = buffer_2000_pdpd[595];

    auto g_xy_0_0_0_z_yy_x_xz = buffer_2000_pdpd[596];

    auto g_xy_0_0_0_z_yy_x_yy = buffer_2000_pdpd[597];

    auto g_xy_0_0_0_z_yy_x_yz = buffer_2000_pdpd[598];

    auto g_xy_0_0_0_z_yy_x_zz = buffer_2000_pdpd[599];

    auto g_xy_0_0_0_z_yy_y_xx = buffer_2000_pdpd[600];

    auto g_xy_0_0_0_z_yy_y_xy = buffer_2000_pdpd[601];

    auto g_xy_0_0_0_z_yy_y_xz = buffer_2000_pdpd[602];

    auto g_xy_0_0_0_z_yy_y_yy = buffer_2000_pdpd[603];

    auto g_xy_0_0_0_z_yy_y_yz = buffer_2000_pdpd[604];

    auto g_xy_0_0_0_z_yy_y_zz = buffer_2000_pdpd[605];

    auto g_xy_0_0_0_z_yy_z_xx = buffer_2000_pdpd[606];

    auto g_xy_0_0_0_z_yy_z_xy = buffer_2000_pdpd[607];

    auto g_xy_0_0_0_z_yy_z_xz = buffer_2000_pdpd[608];

    auto g_xy_0_0_0_z_yy_z_yy = buffer_2000_pdpd[609];

    auto g_xy_0_0_0_z_yy_z_yz = buffer_2000_pdpd[610];

    auto g_xy_0_0_0_z_yy_z_zz = buffer_2000_pdpd[611];

    auto g_xy_0_0_0_z_yz_x_xx = buffer_2000_pdpd[612];

    auto g_xy_0_0_0_z_yz_x_xy = buffer_2000_pdpd[613];

    auto g_xy_0_0_0_z_yz_x_xz = buffer_2000_pdpd[614];

    auto g_xy_0_0_0_z_yz_x_yy = buffer_2000_pdpd[615];

    auto g_xy_0_0_0_z_yz_x_yz = buffer_2000_pdpd[616];

    auto g_xy_0_0_0_z_yz_x_zz = buffer_2000_pdpd[617];

    auto g_xy_0_0_0_z_yz_y_xx = buffer_2000_pdpd[618];

    auto g_xy_0_0_0_z_yz_y_xy = buffer_2000_pdpd[619];

    auto g_xy_0_0_0_z_yz_y_xz = buffer_2000_pdpd[620];

    auto g_xy_0_0_0_z_yz_y_yy = buffer_2000_pdpd[621];

    auto g_xy_0_0_0_z_yz_y_yz = buffer_2000_pdpd[622];

    auto g_xy_0_0_0_z_yz_y_zz = buffer_2000_pdpd[623];

    auto g_xy_0_0_0_z_yz_z_xx = buffer_2000_pdpd[624];

    auto g_xy_0_0_0_z_yz_z_xy = buffer_2000_pdpd[625];

    auto g_xy_0_0_0_z_yz_z_xz = buffer_2000_pdpd[626];

    auto g_xy_0_0_0_z_yz_z_yy = buffer_2000_pdpd[627];

    auto g_xy_0_0_0_z_yz_z_yz = buffer_2000_pdpd[628];

    auto g_xy_0_0_0_z_yz_z_zz = buffer_2000_pdpd[629];

    auto g_xy_0_0_0_z_zz_x_xx = buffer_2000_pdpd[630];

    auto g_xy_0_0_0_z_zz_x_xy = buffer_2000_pdpd[631];

    auto g_xy_0_0_0_z_zz_x_xz = buffer_2000_pdpd[632];

    auto g_xy_0_0_0_z_zz_x_yy = buffer_2000_pdpd[633];

    auto g_xy_0_0_0_z_zz_x_yz = buffer_2000_pdpd[634];

    auto g_xy_0_0_0_z_zz_x_zz = buffer_2000_pdpd[635];

    auto g_xy_0_0_0_z_zz_y_xx = buffer_2000_pdpd[636];

    auto g_xy_0_0_0_z_zz_y_xy = buffer_2000_pdpd[637];

    auto g_xy_0_0_0_z_zz_y_xz = buffer_2000_pdpd[638];

    auto g_xy_0_0_0_z_zz_y_yy = buffer_2000_pdpd[639];

    auto g_xy_0_0_0_z_zz_y_yz = buffer_2000_pdpd[640];

    auto g_xy_0_0_0_z_zz_y_zz = buffer_2000_pdpd[641];

    auto g_xy_0_0_0_z_zz_z_xx = buffer_2000_pdpd[642];

    auto g_xy_0_0_0_z_zz_z_xy = buffer_2000_pdpd[643];

    auto g_xy_0_0_0_z_zz_z_xz = buffer_2000_pdpd[644];

    auto g_xy_0_0_0_z_zz_z_yy = buffer_2000_pdpd[645];

    auto g_xy_0_0_0_z_zz_z_yz = buffer_2000_pdpd[646];

    auto g_xy_0_0_0_z_zz_z_zz = buffer_2000_pdpd[647];

    auto g_xz_0_0_0_x_xx_x_xx = buffer_2000_pdpd[648];

    auto g_xz_0_0_0_x_xx_x_xy = buffer_2000_pdpd[649];

    auto g_xz_0_0_0_x_xx_x_xz = buffer_2000_pdpd[650];

    auto g_xz_0_0_0_x_xx_x_yy = buffer_2000_pdpd[651];

    auto g_xz_0_0_0_x_xx_x_yz = buffer_2000_pdpd[652];

    auto g_xz_0_0_0_x_xx_x_zz = buffer_2000_pdpd[653];

    auto g_xz_0_0_0_x_xx_y_xx = buffer_2000_pdpd[654];

    auto g_xz_0_0_0_x_xx_y_xy = buffer_2000_pdpd[655];

    auto g_xz_0_0_0_x_xx_y_xz = buffer_2000_pdpd[656];

    auto g_xz_0_0_0_x_xx_y_yy = buffer_2000_pdpd[657];

    auto g_xz_0_0_0_x_xx_y_yz = buffer_2000_pdpd[658];

    auto g_xz_0_0_0_x_xx_y_zz = buffer_2000_pdpd[659];

    auto g_xz_0_0_0_x_xx_z_xx = buffer_2000_pdpd[660];

    auto g_xz_0_0_0_x_xx_z_xy = buffer_2000_pdpd[661];

    auto g_xz_0_0_0_x_xx_z_xz = buffer_2000_pdpd[662];

    auto g_xz_0_0_0_x_xx_z_yy = buffer_2000_pdpd[663];

    auto g_xz_0_0_0_x_xx_z_yz = buffer_2000_pdpd[664];

    auto g_xz_0_0_0_x_xx_z_zz = buffer_2000_pdpd[665];

    auto g_xz_0_0_0_x_xy_x_xx = buffer_2000_pdpd[666];

    auto g_xz_0_0_0_x_xy_x_xy = buffer_2000_pdpd[667];

    auto g_xz_0_0_0_x_xy_x_xz = buffer_2000_pdpd[668];

    auto g_xz_0_0_0_x_xy_x_yy = buffer_2000_pdpd[669];

    auto g_xz_0_0_0_x_xy_x_yz = buffer_2000_pdpd[670];

    auto g_xz_0_0_0_x_xy_x_zz = buffer_2000_pdpd[671];

    auto g_xz_0_0_0_x_xy_y_xx = buffer_2000_pdpd[672];

    auto g_xz_0_0_0_x_xy_y_xy = buffer_2000_pdpd[673];

    auto g_xz_0_0_0_x_xy_y_xz = buffer_2000_pdpd[674];

    auto g_xz_0_0_0_x_xy_y_yy = buffer_2000_pdpd[675];

    auto g_xz_0_0_0_x_xy_y_yz = buffer_2000_pdpd[676];

    auto g_xz_0_0_0_x_xy_y_zz = buffer_2000_pdpd[677];

    auto g_xz_0_0_0_x_xy_z_xx = buffer_2000_pdpd[678];

    auto g_xz_0_0_0_x_xy_z_xy = buffer_2000_pdpd[679];

    auto g_xz_0_0_0_x_xy_z_xz = buffer_2000_pdpd[680];

    auto g_xz_0_0_0_x_xy_z_yy = buffer_2000_pdpd[681];

    auto g_xz_0_0_0_x_xy_z_yz = buffer_2000_pdpd[682];

    auto g_xz_0_0_0_x_xy_z_zz = buffer_2000_pdpd[683];

    auto g_xz_0_0_0_x_xz_x_xx = buffer_2000_pdpd[684];

    auto g_xz_0_0_0_x_xz_x_xy = buffer_2000_pdpd[685];

    auto g_xz_0_0_0_x_xz_x_xz = buffer_2000_pdpd[686];

    auto g_xz_0_0_0_x_xz_x_yy = buffer_2000_pdpd[687];

    auto g_xz_0_0_0_x_xz_x_yz = buffer_2000_pdpd[688];

    auto g_xz_0_0_0_x_xz_x_zz = buffer_2000_pdpd[689];

    auto g_xz_0_0_0_x_xz_y_xx = buffer_2000_pdpd[690];

    auto g_xz_0_0_0_x_xz_y_xy = buffer_2000_pdpd[691];

    auto g_xz_0_0_0_x_xz_y_xz = buffer_2000_pdpd[692];

    auto g_xz_0_0_0_x_xz_y_yy = buffer_2000_pdpd[693];

    auto g_xz_0_0_0_x_xz_y_yz = buffer_2000_pdpd[694];

    auto g_xz_0_0_0_x_xz_y_zz = buffer_2000_pdpd[695];

    auto g_xz_0_0_0_x_xz_z_xx = buffer_2000_pdpd[696];

    auto g_xz_0_0_0_x_xz_z_xy = buffer_2000_pdpd[697];

    auto g_xz_0_0_0_x_xz_z_xz = buffer_2000_pdpd[698];

    auto g_xz_0_0_0_x_xz_z_yy = buffer_2000_pdpd[699];

    auto g_xz_0_0_0_x_xz_z_yz = buffer_2000_pdpd[700];

    auto g_xz_0_0_0_x_xz_z_zz = buffer_2000_pdpd[701];

    auto g_xz_0_0_0_x_yy_x_xx = buffer_2000_pdpd[702];

    auto g_xz_0_0_0_x_yy_x_xy = buffer_2000_pdpd[703];

    auto g_xz_0_0_0_x_yy_x_xz = buffer_2000_pdpd[704];

    auto g_xz_0_0_0_x_yy_x_yy = buffer_2000_pdpd[705];

    auto g_xz_0_0_0_x_yy_x_yz = buffer_2000_pdpd[706];

    auto g_xz_0_0_0_x_yy_x_zz = buffer_2000_pdpd[707];

    auto g_xz_0_0_0_x_yy_y_xx = buffer_2000_pdpd[708];

    auto g_xz_0_0_0_x_yy_y_xy = buffer_2000_pdpd[709];

    auto g_xz_0_0_0_x_yy_y_xz = buffer_2000_pdpd[710];

    auto g_xz_0_0_0_x_yy_y_yy = buffer_2000_pdpd[711];

    auto g_xz_0_0_0_x_yy_y_yz = buffer_2000_pdpd[712];

    auto g_xz_0_0_0_x_yy_y_zz = buffer_2000_pdpd[713];

    auto g_xz_0_0_0_x_yy_z_xx = buffer_2000_pdpd[714];

    auto g_xz_0_0_0_x_yy_z_xy = buffer_2000_pdpd[715];

    auto g_xz_0_0_0_x_yy_z_xz = buffer_2000_pdpd[716];

    auto g_xz_0_0_0_x_yy_z_yy = buffer_2000_pdpd[717];

    auto g_xz_0_0_0_x_yy_z_yz = buffer_2000_pdpd[718];

    auto g_xz_0_0_0_x_yy_z_zz = buffer_2000_pdpd[719];

    auto g_xz_0_0_0_x_yz_x_xx = buffer_2000_pdpd[720];

    auto g_xz_0_0_0_x_yz_x_xy = buffer_2000_pdpd[721];

    auto g_xz_0_0_0_x_yz_x_xz = buffer_2000_pdpd[722];

    auto g_xz_0_0_0_x_yz_x_yy = buffer_2000_pdpd[723];

    auto g_xz_0_0_0_x_yz_x_yz = buffer_2000_pdpd[724];

    auto g_xz_0_0_0_x_yz_x_zz = buffer_2000_pdpd[725];

    auto g_xz_0_0_0_x_yz_y_xx = buffer_2000_pdpd[726];

    auto g_xz_0_0_0_x_yz_y_xy = buffer_2000_pdpd[727];

    auto g_xz_0_0_0_x_yz_y_xz = buffer_2000_pdpd[728];

    auto g_xz_0_0_0_x_yz_y_yy = buffer_2000_pdpd[729];

    auto g_xz_0_0_0_x_yz_y_yz = buffer_2000_pdpd[730];

    auto g_xz_0_0_0_x_yz_y_zz = buffer_2000_pdpd[731];

    auto g_xz_0_0_0_x_yz_z_xx = buffer_2000_pdpd[732];

    auto g_xz_0_0_0_x_yz_z_xy = buffer_2000_pdpd[733];

    auto g_xz_0_0_0_x_yz_z_xz = buffer_2000_pdpd[734];

    auto g_xz_0_0_0_x_yz_z_yy = buffer_2000_pdpd[735];

    auto g_xz_0_0_0_x_yz_z_yz = buffer_2000_pdpd[736];

    auto g_xz_0_0_0_x_yz_z_zz = buffer_2000_pdpd[737];

    auto g_xz_0_0_0_x_zz_x_xx = buffer_2000_pdpd[738];

    auto g_xz_0_0_0_x_zz_x_xy = buffer_2000_pdpd[739];

    auto g_xz_0_0_0_x_zz_x_xz = buffer_2000_pdpd[740];

    auto g_xz_0_0_0_x_zz_x_yy = buffer_2000_pdpd[741];

    auto g_xz_0_0_0_x_zz_x_yz = buffer_2000_pdpd[742];

    auto g_xz_0_0_0_x_zz_x_zz = buffer_2000_pdpd[743];

    auto g_xz_0_0_0_x_zz_y_xx = buffer_2000_pdpd[744];

    auto g_xz_0_0_0_x_zz_y_xy = buffer_2000_pdpd[745];

    auto g_xz_0_0_0_x_zz_y_xz = buffer_2000_pdpd[746];

    auto g_xz_0_0_0_x_zz_y_yy = buffer_2000_pdpd[747];

    auto g_xz_0_0_0_x_zz_y_yz = buffer_2000_pdpd[748];

    auto g_xz_0_0_0_x_zz_y_zz = buffer_2000_pdpd[749];

    auto g_xz_0_0_0_x_zz_z_xx = buffer_2000_pdpd[750];

    auto g_xz_0_0_0_x_zz_z_xy = buffer_2000_pdpd[751];

    auto g_xz_0_0_0_x_zz_z_xz = buffer_2000_pdpd[752];

    auto g_xz_0_0_0_x_zz_z_yy = buffer_2000_pdpd[753];

    auto g_xz_0_0_0_x_zz_z_yz = buffer_2000_pdpd[754];

    auto g_xz_0_0_0_x_zz_z_zz = buffer_2000_pdpd[755];

    auto g_xz_0_0_0_y_xx_x_xx = buffer_2000_pdpd[756];

    auto g_xz_0_0_0_y_xx_x_xy = buffer_2000_pdpd[757];

    auto g_xz_0_0_0_y_xx_x_xz = buffer_2000_pdpd[758];

    auto g_xz_0_0_0_y_xx_x_yy = buffer_2000_pdpd[759];

    auto g_xz_0_0_0_y_xx_x_yz = buffer_2000_pdpd[760];

    auto g_xz_0_0_0_y_xx_x_zz = buffer_2000_pdpd[761];

    auto g_xz_0_0_0_y_xx_y_xx = buffer_2000_pdpd[762];

    auto g_xz_0_0_0_y_xx_y_xy = buffer_2000_pdpd[763];

    auto g_xz_0_0_0_y_xx_y_xz = buffer_2000_pdpd[764];

    auto g_xz_0_0_0_y_xx_y_yy = buffer_2000_pdpd[765];

    auto g_xz_0_0_0_y_xx_y_yz = buffer_2000_pdpd[766];

    auto g_xz_0_0_0_y_xx_y_zz = buffer_2000_pdpd[767];

    auto g_xz_0_0_0_y_xx_z_xx = buffer_2000_pdpd[768];

    auto g_xz_0_0_0_y_xx_z_xy = buffer_2000_pdpd[769];

    auto g_xz_0_0_0_y_xx_z_xz = buffer_2000_pdpd[770];

    auto g_xz_0_0_0_y_xx_z_yy = buffer_2000_pdpd[771];

    auto g_xz_0_0_0_y_xx_z_yz = buffer_2000_pdpd[772];

    auto g_xz_0_0_0_y_xx_z_zz = buffer_2000_pdpd[773];

    auto g_xz_0_0_0_y_xy_x_xx = buffer_2000_pdpd[774];

    auto g_xz_0_0_0_y_xy_x_xy = buffer_2000_pdpd[775];

    auto g_xz_0_0_0_y_xy_x_xz = buffer_2000_pdpd[776];

    auto g_xz_0_0_0_y_xy_x_yy = buffer_2000_pdpd[777];

    auto g_xz_0_0_0_y_xy_x_yz = buffer_2000_pdpd[778];

    auto g_xz_0_0_0_y_xy_x_zz = buffer_2000_pdpd[779];

    auto g_xz_0_0_0_y_xy_y_xx = buffer_2000_pdpd[780];

    auto g_xz_0_0_0_y_xy_y_xy = buffer_2000_pdpd[781];

    auto g_xz_0_0_0_y_xy_y_xz = buffer_2000_pdpd[782];

    auto g_xz_0_0_0_y_xy_y_yy = buffer_2000_pdpd[783];

    auto g_xz_0_0_0_y_xy_y_yz = buffer_2000_pdpd[784];

    auto g_xz_0_0_0_y_xy_y_zz = buffer_2000_pdpd[785];

    auto g_xz_0_0_0_y_xy_z_xx = buffer_2000_pdpd[786];

    auto g_xz_0_0_0_y_xy_z_xy = buffer_2000_pdpd[787];

    auto g_xz_0_0_0_y_xy_z_xz = buffer_2000_pdpd[788];

    auto g_xz_0_0_0_y_xy_z_yy = buffer_2000_pdpd[789];

    auto g_xz_0_0_0_y_xy_z_yz = buffer_2000_pdpd[790];

    auto g_xz_0_0_0_y_xy_z_zz = buffer_2000_pdpd[791];

    auto g_xz_0_0_0_y_xz_x_xx = buffer_2000_pdpd[792];

    auto g_xz_0_0_0_y_xz_x_xy = buffer_2000_pdpd[793];

    auto g_xz_0_0_0_y_xz_x_xz = buffer_2000_pdpd[794];

    auto g_xz_0_0_0_y_xz_x_yy = buffer_2000_pdpd[795];

    auto g_xz_0_0_0_y_xz_x_yz = buffer_2000_pdpd[796];

    auto g_xz_0_0_0_y_xz_x_zz = buffer_2000_pdpd[797];

    auto g_xz_0_0_0_y_xz_y_xx = buffer_2000_pdpd[798];

    auto g_xz_0_0_0_y_xz_y_xy = buffer_2000_pdpd[799];

    auto g_xz_0_0_0_y_xz_y_xz = buffer_2000_pdpd[800];

    auto g_xz_0_0_0_y_xz_y_yy = buffer_2000_pdpd[801];

    auto g_xz_0_0_0_y_xz_y_yz = buffer_2000_pdpd[802];

    auto g_xz_0_0_0_y_xz_y_zz = buffer_2000_pdpd[803];

    auto g_xz_0_0_0_y_xz_z_xx = buffer_2000_pdpd[804];

    auto g_xz_0_0_0_y_xz_z_xy = buffer_2000_pdpd[805];

    auto g_xz_0_0_0_y_xz_z_xz = buffer_2000_pdpd[806];

    auto g_xz_0_0_0_y_xz_z_yy = buffer_2000_pdpd[807];

    auto g_xz_0_0_0_y_xz_z_yz = buffer_2000_pdpd[808];

    auto g_xz_0_0_0_y_xz_z_zz = buffer_2000_pdpd[809];

    auto g_xz_0_0_0_y_yy_x_xx = buffer_2000_pdpd[810];

    auto g_xz_0_0_0_y_yy_x_xy = buffer_2000_pdpd[811];

    auto g_xz_0_0_0_y_yy_x_xz = buffer_2000_pdpd[812];

    auto g_xz_0_0_0_y_yy_x_yy = buffer_2000_pdpd[813];

    auto g_xz_0_0_0_y_yy_x_yz = buffer_2000_pdpd[814];

    auto g_xz_0_0_0_y_yy_x_zz = buffer_2000_pdpd[815];

    auto g_xz_0_0_0_y_yy_y_xx = buffer_2000_pdpd[816];

    auto g_xz_0_0_0_y_yy_y_xy = buffer_2000_pdpd[817];

    auto g_xz_0_0_0_y_yy_y_xz = buffer_2000_pdpd[818];

    auto g_xz_0_0_0_y_yy_y_yy = buffer_2000_pdpd[819];

    auto g_xz_0_0_0_y_yy_y_yz = buffer_2000_pdpd[820];

    auto g_xz_0_0_0_y_yy_y_zz = buffer_2000_pdpd[821];

    auto g_xz_0_0_0_y_yy_z_xx = buffer_2000_pdpd[822];

    auto g_xz_0_0_0_y_yy_z_xy = buffer_2000_pdpd[823];

    auto g_xz_0_0_0_y_yy_z_xz = buffer_2000_pdpd[824];

    auto g_xz_0_0_0_y_yy_z_yy = buffer_2000_pdpd[825];

    auto g_xz_0_0_0_y_yy_z_yz = buffer_2000_pdpd[826];

    auto g_xz_0_0_0_y_yy_z_zz = buffer_2000_pdpd[827];

    auto g_xz_0_0_0_y_yz_x_xx = buffer_2000_pdpd[828];

    auto g_xz_0_0_0_y_yz_x_xy = buffer_2000_pdpd[829];

    auto g_xz_0_0_0_y_yz_x_xz = buffer_2000_pdpd[830];

    auto g_xz_0_0_0_y_yz_x_yy = buffer_2000_pdpd[831];

    auto g_xz_0_0_0_y_yz_x_yz = buffer_2000_pdpd[832];

    auto g_xz_0_0_0_y_yz_x_zz = buffer_2000_pdpd[833];

    auto g_xz_0_0_0_y_yz_y_xx = buffer_2000_pdpd[834];

    auto g_xz_0_0_0_y_yz_y_xy = buffer_2000_pdpd[835];

    auto g_xz_0_0_0_y_yz_y_xz = buffer_2000_pdpd[836];

    auto g_xz_0_0_0_y_yz_y_yy = buffer_2000_pdpd[837];

    auto g_xz_0_0_0_y_yz_y_yz = buffer_2000_pdpd[838];

    auto g_xz_0_0_0_y_yz_y_zz = buffer_2000_pdpd[839];

    auto g_xz_0_0_0_y_yz_z_xx = buffer_2000_pdpd[840];

    auto g_xz_0_0_0_y_yz_z_xy = buffer_2000_pdpd[841];

    auto g_xz_0_0_0_y_yz_z_xz = buffer_2000_pdpd[842];

    auto g_xz_0_0_0_y_yz_z_yy = buffer_2000_pdpd[843];

    auto g_xz_0_0_0_y_yz_z_yz = buffer_2000_pdpd[844];

    auto g_xz_0_0_0_y_yz_z_zz = buffer_2000_pdpd[845];

    auto g_xz_0_0_0_y_zz_x_xx = buffer_2000_pdpd[846];

    auto g_xz_0_0_0_y_zz_x_xy = buffer_2000_pdpd[847];

    auto g_xz_0_0_0_y_zz_x_xz = buffer_2000_pdpd[848];

    auto g_xz_0_0_0_y_zz_x_yy = buffer_2000_pdpd[849];

    auto g_xz_0_0_0_y_zz_x_yz = buffer_2000_pdpd[850];

    auto g_xz_0_0_0_y_zz_x_zz = buffer_2000_pdpd[851];

    auto g_xz_0_0_0_y_zz_y_xx = buffer_2000_pdpd[852];

    auto g_xz_0_0_0_y_zz_y_xy = buffer_2000_pdpd[853];

    auto g_xz_0_0_0_y_zz_y_xz = buffer_2000_pdpd[854];

    auto g_xz_0_0_0_y_zz_y_yy = buffer_2000_pdpd[855];

    auto g_xz_0_0_0_y_zz_y_yz = buffer_2000_pdpd[856];

    auto g_xz_0_0_0_y_zz_y_zz = buffer_2000_pdpd[857];

    auto g_xz_0_0_0_y_zz_z_xx = buffer_2000_pdpd[858];

    auto g_xz_0_0_0_y_zz_z_xy = buffer_2000_pdpd[859];

    auto g_xz_0_0_0_y_zz_z_xz = buffer_2000_pdpd[860];

    auto g_xz_0_0_0_y_zz_z_yy = buffer_2000_pdpd[861];

    auto g_xz_0_0_0_y_zz_z_yz = buffer_2000_pdpd[862];

    auto g_xz_0_0_0_y_zz_z_zz = buffer_2000_pdpd[863];

    auto g_xz_0_0_0_z_xx_x_xx = buffer_2000_pdpd[864];

    auto g_xz_0_0_0_z_xx_x_xy = buffer_2000_pdpd[865];

    auto g_xz_0_0_0_z_xx_x_xz = buffer_2000_pdpd[866];

    auto g_xz_0_0_0_z_xx_x_yy = buffer_2000_pdpd[867];

    auto g_xz_0_0_0_z_xx_x_yz = buffer_2000_pdpd[868];

    auto g_xz_0_0_0_z_xx_x_zz = buffer_2000_pdpd[869];

    auto g_xz_0_0_0_z_xx_y_xx = buffer_2000_pdpd[870];

    auto g_xz_0_0_0_z_xx_y_xy = buffer_2000_pdpd[871];

    auto g_xz_0_0_0_z_xx_y_xz = buffer_2000_pdpd[872];

    auto g_xz_0_0_0_z_xx_y_yy = buffer_2000_pdpd[873];

    auto g_xz_0_0_0_z_xx_y_yz = buffer_2000_pdpd[874];

    auto g_xz_0_0_0_z_xx_y_zz = buffer_2000_pdpd[875];

    auto g_xz_0_0_0_z_xx_z_xx = buffer_2000_pdpd[876];

    auto g_xz_0_0_0_z_xx_z_xy = buffer_2000_pdpd[877];

    auto g_xz_0_0_0_z_xx_z_xz = buffer_2000_pdpd[878];

    auto g_xz_0_0_0_z_xx_z_yy = buffer_2000_pdpd[879];

    auto g_xz_0_0_0_z_xx_z_yz = buffer_2000_pdpd[880];

    auto g_xz_0_0_0_z_xx_z_zz = buffer_2000_pdpd[881];

    auto g_xz_0_0_0_z_xy_x_xx = buffer_2000_pdpd[882];

    auto g_xz_0_0_0_z_xy_x_xy = buffer_2000_pdpd[883];

    auto g_xz_0_0_0_z_xy_x_xz = buffer_2000_pdpd[884];

    auto g_xz_0_0_0_z_xy_x_yy = buffer_2000_pdpd[885];

    auto g_xz_0_0_0_z_xy_x_yz = buffer_2000_pdpd[886];

    auto g_xz_0_0_0_z_xy_x_zz = buffer_2000_pdpd[887];

    auto g_xz_0_0_0_z_xy_y_xx = buffer_2000_pdpd[888];

    auto g_xz_0_0_0_z_xy_y_xy = buffer_2000_pdpd[889];

    auto g_xz_0_0_0_z_xy_y_xz = buffer_2000_pdpd[890];

    auto g_xz_0_0_0_z_xy_y_yy = buffer_2000_pdpd[891];

    auto g_xz_0_0_0_z_xy_y_yz = buffer_2000_pdpd[892];

    auto g_xz_0_0_0_z_xy_y_zz = buffer_2000_pdpd[893];

    auto g_xz_0_0_0_z_xy_z_xx = buffer_2000_pdpd[894];

    auto g_xz_0_0_0_z_xy_z_xy = buffer_2000_pdpd[895];

    auto g_xz_0_0_0_z_xy_z_xz = buffer_2000_pdpd[896];

    auto g_xz_0_0_0_z_xy_z_yy = buffer_2000_pdpd[897];

    auto g_xz_0_0_0_z_xy_z_yz = buffer_2000_pdpd[898];

    auto g_xz_0_0_0_z_xy_z_zz = buffer_2000_pdpd[899];

    auto g_xz_0_0_0_z_xz_x_xx = buffer_2000_pdpd[900];

    auto g_xz_0_0_0_z_xz_x_xy = buffer_2000_pdpd[901];

    auto g_xz_0_0_0_z_xz_x_xz = buffer_2000_pdpd[902];

    auto g_xz_0_0_0_z_xz_x_yy = buffer_2000_pdpd[903];

    auto g_xz_0_0_0_z_xz_x_yz = buffer_2000_pdpd[904];

    auto g_xz_0_0_0_z_xz_x_zz = buffer_2000_pdpd[905];

    auto g_xz_0_0_0_z_xz_y_xx = buffer_2000_pdpd[906];

    auto g_xz_0_0_0_z_xz_y_xy = buffer_2000_pdpd[907];

    auto g_xz_0_0_0_z_xz_y_xz = buffer_2000_pdpd[908];

    auto g_xz_0_0_0_z_xz_y_yy = buffer_2000_pdpd[909];

    auto g_xz_0_0_0_z_xz_y_yz = buffer_2000_pdpd[910];

    auto g_xz_0_0_0_z_xz_y_zz = buffer_2000_pdpd[911];

    auto g_xz_0_0_0_z_xz_z_xx = buffer_2000_pdpd[912];

    auto g_xz_0_0_0_z_xz_z_xy = buffer_2000_pdpd[913];

    auto g_xz_0_0_0_z_xz_z_xz = buffer_2000_pdpd[914];

    auto g_xz_0_0_0_z_xz_z_yy = buffer_2000_pdpd[915];

    auto g_xz_0_0_0_z_xz_z_yz = buffer_2000_pdpd[916];

    auto g_xz_0_0_0_z_xz_z_zz = buffer_2000_pdpd[917];

    auto g_xz_0_0_0_z_yy_x_xx = buffer_2000_pdpd[918];

    auto g_xz_0_0_0_z_yy_x_xy = buffer_2000_pdpd[919];

    auto g_xz_0_0_0_z_yy_x_xz = buffer_2000_pdpd[920];

    auto g_xz_0_0_0_z_yy_x_yy = buffer_2000_pdpd[921];

    auto g_xz_0_0_0_z_yy_x_yz = buffer_2000_pdpd[922];

    auto g_xz_0_0_0_z_yy_x_zz = buffer_2000_pdpd[923];

    auto g_xz_0_0_0_z_yy_y_xx = buffer_2000_pdpd[924];

    auto g_xz_0_0_0_z_yy_y_xy = buffer_2000_pdpd[925];

    auto g_xz_0_0_0_z_yy_y_xz = buffer_2000_pdpd[926];

    auto g_xz_0_0_0_z_yy_y_yy = buffer_2000_pdpd[927];

    auto g_xz_0_0_0_z_yy_y_yz = buffer_2000_pdpd[928];

    auto g_xz_0_0_0_z_yy_y_zz = buffer_2000_pdpd[929];

    auto g_xz_0_0_0_z_yy_z_xx = buffer_2000_pdpd[930];

    auto g_xz_0_0_0_z_yy_z_xy = buffer_2000_pdpd[931];

    auto g_xz_0_0_0_z_yy_z_xz = buffer_2000_pdpd[932];

    auto g_xz_0_0_0_z_yy_z_yy = buffer_2000_pdpd[933];

    auto g_xz_0_0_0_z_yy_z_yz = buffer_2000_pdpd[934];

    auto g_xz_0_0_0_z_yy_z_zz = buffer_2000_pdpd[935];

    auto g_xz_0_0_0_z_yz_x_xx = buffer_2000_pdpd[936];

    auto g_xz_0_0_0_z_yz_x_xy = buffer_2000_pdpd[937];

    auto g_xz_0_0_0_z_yz_x_xz = buffer_2000_pdpd[938];

    auto g_xz_0_0_0_z_yz_x_yy = buffer_2000_pdpd[939];

    auto g_xz_0_0_0_z_yz_x_yz = buffer_2000_pdpd[940];

    auto g_xz_0_0_0_z_yz_x_zz = buffer_2000_pdpd[941];

    auto g_xz_0_0_0_z_yz_y_xx = buffer_2000_pdpd[942];

    auto g_xz_0_0_0_z_yz_y_xy = buffer_2000_pdpd[943];

    auto g_xz_0_0_0_z_yz_y_xz = buffer_2000_pdpd[944];

    auto g_xz_0_0_0_z_yz_y_yy = buffer_2000_pdpd[945];

    auto g_xz_0_0_0_z_yz_y_yz = buffer_2000_pdpd[946];

    auto g_xz_0_0_0_z_yz_y_zz = buffer_2000_pdpd[947];

    auto g_xz_0_0_0_z_yz_z_xx = buffer_2000_pdpd[948];

    auto g_xz_0_0_0_z_yz_z_xy = buffer_2000_pdpd[949];

    auto g_xz_0_0_0_z_yz_z_xz = buffer_2000_pdpd[950];

    auto g_xz_0_0_0_z_yz_z_yy = buffer_2000_pdpd[951];

    auto g_xz_0_0_0_z_yz_z_yz = buffer_2000_pdpd[952];

    auto g_xz_0_0_0_z_yz_z_zz = buffer_2000_pdpd[953];

    auto g_xz_0_0_0_z_zz_x_xx = buffer_2000_pdpd[954];

    auto g_xz_0_0_0_z_zz_x_xy = buffer_2000_pdpd[955];

    auto g_xz_0_0_0_z_zz_x_xz = buffer_2000_pdpd[956];

    auto g_xz_0_0_0_z_zz_x_yy = buffer_2000_pdpd[957];

    auto g_xz_0_0_0_z_zz_x_yz = buffer_2000_pdpd[958];

    auto g_xz_0_0_0_z_zz_x_zz = buffer_2000_pdpd[959];

    auto g_xz_0_0_0_z_zz_y_xx = buffer_2000_pdpd[960];

    auto g_xz_0_0_0_z_zz_y_xy = buffer_2000_pdpd[961];

    auto g_xz_0_0_0_z_zz_y_xz = buffer_2000_pdpd[962];

    auto g_xz_0_0_0_z_zz_y_yy = buffer_2000_pdpd[963];

    auto g_xz_0_0_0_z_zz_y_yz = buffer_2000_pdpd[964];

    auto g_xz_0_0_0_z_zz_y_zz = buffer_2000_pdpd[965];

    auto g_xz_0_0_0_z_zz_z_xx = buffer_2000_pdpd[966];

    auto g_xz_0_0_0_z_zz_z_xy = buffer_2000_pdpd[967];

    auto g_xz_0_0_0_z_zz_z_xz = buffer_2000_pdpd[968];

    auto g_xz_0_0_0_z_zz_z_yy = buffer_2000_pdpd[969];

    auto g_xz_0_0_0_z_zz_z_yz = buffer_2000_pdpd[970];

    auto g_xz_0_0_0_z_zz_z_zz = buffer_2000_pdpd[971];

    auto g_yy_0_0_0_x_xx_x_xx = buffer_2000_pdpd[972];

    auto g_yy_0_0_0_x_xx_x_xy = buffer_2000_pdpd[973];

    auto g_yy_0_0_0_x_xx_x_xz = buffer_2000_pdpd[974];

    auto g_yy_0_0_0_x_xx_x_yy = buffer_2000_pdpd[975];

    auto g_yy_0_0_0_x_xx_x_yz = buffer_2000_pdpd[976];

    auto g_yy_0_0_0_x_xx_x_zz = buffer_2000_pdpd[977];

    auto g_yy_0_0_0_x_xx_y_xx = buffer_2000_pdpd[978];

    auto g_yy_0_0_0_x_xx_y_xy = buffer_2000_pdpd[979];

    auto g_yy_0_0_0_x_xx_y_xz = buffer_2000_pdpd[980];

    auto g_yy_0_0_0_x_xx_y_yy = buffer_2000_pdpd[981];

    auto g_yy_0_0_0_x_xx_y_yz = buffer_2000_pdpd[982];

    auto g_yy_0_0_0_x_xx_y_zz = buffer_2000_pdpd[983];

    auto g_yy_0_0_0_x_xx_z_xx = buffer_2000_pdpd[984];

    auto g_yy_0_0_0_x_xx_z_xy = buffer_2000_pdpd[985];

    auto g_yy_0_0_0_x_xx_z_xz = buffer_2000_pdpd[986];

    auto g_yy_0_0_0_x_xx_z_yy = buffer_2000_pdpd[987];

    auto g_yy_0_0_0_x_xx_z_yz = buffer_2000_pdpd[988];

    auto g_yy_0_0_0_x_xx_z_zz = buffer_2000_pdpd[989];

    auto g_yy_0_0_0_x_xy_x_xx = buffer_2000_pdpd[990];

    auto g_yy_0_0_0_x_xy_x_xy = buffer_2000_pdpd[991];

    auto g_yy_0_0_0_x_xy_x_xz = buffer_2000_pdpd[992];

    auto g_yy_0_0_0_x_xy_x_yy = buffer_2000_pdpd[993];

    auto g_yy_0_0_0_x_xy_x_yz = buffer_2000_pdpd[994];

    auto g_yy_0_0_0_x_xy_x_zz = buffer_2000_pdpd[995];

    auto g_yy_0_0_0_x_xy_y_xx = buffer_2000_pdpd[996];

    auto g_yy_0_0_0_x_xy_y_xy = buffer_2000_pdpd[997];

    auto g_yy_0_0_0_x_xy_y_xz = buffer_2000_pdpd[998];

    auto g_yy_0_0_0_x_xy_y_yy = buffer_2000_pdpd[999];

    auto g_yy_0_0_0_x_xy_y_yz = buffer_2000_pdpd[1000];

    auto g_yy_0_0_0_x_xy_y_zz = buffer_2000_pdpd[1001];

    auto g_yy_0_0_0_x_xy_z_xx = buffer_2000_pdpd[1002];

    auto g_yy_0_0_0_x_xy_z_xy = buffer_2000_pdpd[1003];

    auto g_yy_0_0_0_x_xy_z_xz = buffer_2000_pdpd[1004];

    auto g_yy_0_0_0_x_xy_z_yy = buffer_2000_pdpd[1005];

    auto g_yy_0_0_0_x_xy_z_yz = buffer_2000_pdpd[1006];

    auto g_yy_0_0_0_x_xy_z_zz = buffer_2000_pdpd[1007];

    auto g_yy_0_0_0_x_xz_x_xx = buffer_2000_pdpd[1008];

    auto g_yy_0_0_0_x_xz_x_xy = buffer_2000_pdpd[1009];

    auto g_yy_0_0_0_x_xz_x_xz = buffer_2000_pdpd[1010];

    auto g_yy_0_0_0_x_xz_x_yy = buffer_2000_pdpd[1011];

    auto g_yy_0_0_0_x_xz_x_yz = buffer_2000_pdpd[1012];

    auto g_yy_0_0_0_x_xz_x_zz = buffer_2000_pdpd[1013];

    auto g_yy_0_0_0_x_xz_y_xx = buffer_2000_pdpd[1014];

    auto g_yy_0_0_0_x_xz_y_xy = buffer_2000_pdpd[1015];

    auto g_yy_0_0_0_x_xz_y_xz = buffer_2000_pdpd[1016];

    auto g_yy_0_0_0_x_xz_y_yy = buffer_2000_pdpd[1017];

    auto g_yy_0_0_0_x_xz_y_yz = buffer_2000_pdpd[1018];

    auto g_yy_0_0_0_x_xz_y_zz = buffer_2000_pdpd[1019];

    auto g_yy_0_0_0_x_xz_z_xx = buffer_2000_pdpd[1020];

    auto g_yy_0_0_0_x_xz_z_xy = buffer_2000_pdpd[1021];

    auto g_yy_0_0_0_x_xz_z_xz = buffer_2000_pdpd[1022];

    auto g_yy_0_0_0_x_xz_z_yy = buffer_2000_pdpd[1023];

    auto g_yy_0_0_0_x_xz_z_yz = buffer_2000_pdpd[1024];

    auto g_yy_0_0_0_x_xz_z_zz = buffer_2000_pdpd[1025];

    auto g_yy_0_0_0_x_yy_x_xx = buffer_2000_pdpd[1026];

    auto g_yy_0_0_0_x_yy_x_xy = buffer_2000_pdpd[1027];

    auto g_yy_0_0_0_x_yy_x_xz = buffer_2000_pdpd[1028];

    auto g_yy_0_0_0_x_yy_x_yy = buffer_2000_pdpd[1029];

    auto g_yy_0_0_0_x_yy_x_yz = buffer_2000_pdpd[1030];

    auto g_yy_0_0_0_x_yy_x_zz = buffer_2000_pdpd[1031];

    auto g_yy_0_0_0_x_yy_y_xx = buffer_2000_pdpd[1032];

    auto g_yy_0_0_0_x_yy_y_xy = buffer_2000_pdpd[1033];

    auto g_yy_0_0_0_x_yy_y_xz = buffer_2000_pdpd[1034];

    auto g_yy_0_0_0_x_yy_y_yy = buffer_2000_pdpd[1035];

    auto g_yy_0_0_0_x_yy_y_yz = buffer_2000_pdpd[1036];

    auto g_yy_0_0_0_x_yy_y_zz = buffer_2000_pdpd[1037];

    auto g_yy_0_0_0_x_yy_z_xx = buffer_2000_pdpd[1038];

    auto g_yy_0_0_0_x_yy_z_xy = buffer_2000_pdpd[1039];

    auto g_yy_0_0_0_x_yy_z_xz = buffer_2000_pdpd[1040];

    auto g_yy_0_0_0_x_yy_z_yy = buffer_2000_pdpd[1041];

    auto g_yy_0_0_0_x_yy_z_yz = buffer_2000_pdpd[1042];

    auto g_yy_0_0_0_x_yy_z_zz = buffer_2000_pdpd[1043];

    auto g_yy_0_0_0_x_yz_x_xx = buffer_2000_pdpd[1044];

    auto g_yy_0_0_0_x_yz_x_xy = buffer_2000_pdpd[1045];

    auto g_yy_0_0_0_x_yz_x_xz = buffer_2000_pdpd[1046];

    auto g_yy_0_0_0_x_yz_x_yy = buffer_2000_pdpd[1047];

    auto g_yy_0_0_0_x_yz_x_yz = buffer_2000_pdpd[1048];

    auto g_yy_0_0_0_x_yz_x_zz = buffer_2000_pdpd[1049];

    auto g_yy_0_0_0_x_yz_y_xx = buffer_2000_pdpd[1050];

    auto g_yy_0_0_0_x_yz_y_xy = buffer_2000_pdpd[1051];

    auto g_yy_0_0_0_x_yz_y_xz = buffer_2000_pdpd[1052];

    auto g_yy_0_0_0_x_yz_y_yy = buffer_2000_pdpd[1053];

    auto g_yy_0_0_0_x_yz_y_yz = buffer_2000_pdpd[1054];

    auto g_yy_0_0_0_x_yz_y_zz = buffer_2000_pdpd[1055];

    auto g_yy_0_0_0_x_yz_z_xx = buffer_2000_pdpd[1056];

    auto g_yy_0_0_0_x_yz_z_xy = buffer_2000_pdpd[1057];

    auto g_yy_0_0_0_x_yz_z_xz = buffer_2000_pdpd[1058];

    auto g_yy_0_0_0_x_yz_z_yy = buffer_2000_pdpd[1059];

    auto g_yy_0_0_0_x_yz_z_yz = buffer_2000_pdpd[1060];

    auto g_yy_0_0_0_x_yz_z_zz = buffer_2000_pdpd[1061];

    auto g_yy_0_0_0_x_zz_x_xx = buffer_2000_pdpd[1062];

    auto g_yy_0_0_0_x_zz_x_xy = buffer_2000_pdpd[1063];

    auto g_yy_0_0_0_x_zz_x_xz = buffer_2000_pdpd[1064];

    auto g_yy_0_0_0_x_zz_x_yy = buffer_2000_pdpd[1065];

    auto g_yy_0_0_0_x_zz_x_yz = buffer_2000_pdpd[1066];

    auto g_yy_0_0_0_x_zz_x_zz = buffer_2000_pdpd[1067];

    auto g_yy_0_0_0_x_zz_y_xx = buffer_2000_pdpd[1068];

    auto g_yy_0_0_0_x_zz_y_xy = buffer_2000_pdpd[1069];

    auto g_yy_0_0_0_x_zz_y_xz = buffer_2000_pdpd[1070];

    auto g_yy_0_0_0_x_zz_y_yy = buffer_2000_pdpd[1071];

    auto g_yy_0_0_0_x_zz_y_yz = buffer_2000_pdpd[1072];

    auto g_yy_0_0_0_x_zz_y_zz = buffer_2000_pdpd[1073];

    auto g_yy_0_0_0_x_zz_z_xx = buffer_2000_pdpd[1074];

    auto g_yy_0_0_0_x_zz_z_xy = buffer_2000_pdpd[1075];

    auto g_yy_0_0_0_x_zz_z_xz = buffer_2000_pdpd[1076];

    auto g_yy_0_0_0_x_zz_z_yy = buffer_2000_pdpd[1077];

    auto g_yy_0_0_0_x_zz_z_yz = buffer_2000_pdpd[1078];

    auto g_yy_0_0_0_x_zz_z_zz = buffer_2000_pdpd[1079];

    auto g_yy_0_0_0_y_xx_x_xx = buffer_2000_pdpd[1080];

    auto g_yy_0_0_0_y_xx_x_xy = buffer_2000_pdpd[1081];

    auto g_yy_0_0_0_y_xx_x_xz = buffer_2000_pdpd[1082];

    auto g_yy_0_0_0_y_xx_x_yy = buffer_2000_pdpd[1083];

    auto g_yy_0_0_0_y_xx_x_yz = buffer_2000_pdpd[1084];

    auto g_yy_0_0_0_y_xx_x_zz = buffer_2000_pdpd[1085];

    auto g_yy_0_0_0_y_xx_y_xx = buffer_2000_pdpd[1086];

    auto g_yy_0_0_0_y_xx_y_xy = buffer_2000_pdpd[1087];

    auto g_yy_0_0_0_y_xx_y_xz = buffer_2000_pdpd[1088];

    auto g_yy_0_0_0_y_xx_y_yy = buffer_2000_pdpd[1089];

    auto g_yy_0_0_0_y_xx_y_yz = buffer_2000_pdpd[1090];

    auto g_yy_0_0_0_y_xx_y_zz = buffer_2000_pdpd[1091];

    auto g_yy_0_0_0_y_xx_z_xx = buffer_2000_pdpd[1092];

    auto g_yy_0_0_0_y_xx_z_xy = buffer_2000_pdpd[1093];

    auto g_yy_0_0_0_y_xx_z_xz = buffer_2000_pdpd[1094];

    auto g_yy_0_0_0_y_xx_z_yy = buffer_2000_pdpd[1095];

    auto g_yy_0_0_0_y_xx_z_yz = buffer_2000_pdpd[1096];

    auto g_yy_0_0_0_y_xx_z_zz = buffer_2000_pdpd[1097];

    auto g_yy_0_0_0_y_xy_x_xx = buffer_2000_pdpd[1098];

    auto g_yy_0_0_0_y_xy_x_xy = buffer_2000_pdpd[1099];

    auto g_yy_0_0_0_y_xy_x_xz = buffer_2000_pdpd[1100];

    auto g_yy_0_0_0_y_xy_x_yy = buffer_2000_pdpd[1101];

    auto g_yy_0_0_0_y_xy_x_yz = buffer_2000_pdpd[1102];

    auto g_yy_0_0_0_y_xy_x_zz = buffer_2000_pdpd[1103];

    auto g_yy_0_0_0_y_xy_y_xx = buffer_2000_pdpd[1104];

    auto g_yy_0_0_0_y_xy_y_xy = buffer_2000_pdpd[1105];

    auto g_yy_0_0_0_y_xy_y_xz = buffer_2000_pdpd[1106];

    auto g_yy_0_0_0_y_xy_y_yy = buffer_2000_pdpd[1107];

    auto g_yy_0_0_0_y_xy_y_yz = buffer_2000_pdpd[1108];

    auto g_yy_0_0_0_y_xy_y_zz = buffer_2000_pdpd[1109];

    auto g_yy_0_0_0_y_xy_z_xx = buffer_2000_pdpd[1110];

    auto g_yy_0_0_0_y_xy_z_xy = buffer_2000_pdpd[1111];

    auto g_yy_0_0_0_y_xy_z_xz = buffer_2000_pdpd[1112];

    auto g_yy_0_0_0_y_xy_z_yy = buffer_2000_pdpd[1113];

    auto g_yy_0_0_0_y_xy_z_yz = buffer_2000_pdpd[1114];

    auto g_yy_0_0_0_y_xy_z_zz = buffer_2000_pdpd[1115];

    auto g_yy_0_0_0_y_xz_x_xx = buffer_2000_pdpd[1116];

    auto g_yy_0_0_0_y_xz_x_xy = buffer_2000_pdpd[1117];

    auto g_yy_0_0_0_y_xz_x_xz = buffer_2000_pdpd[1118];

    auto g_yy_0_0_0_y_xz_x_yy = buffer_2000_pdpd[1119];

    auto g_yy_0_0_0_y_xz_x_yz = buffer_2000_pdpd[1120];

    auto g_yy_0_0_0_y_xz_x_zz = buffer_2000_pdpd[1121];

    auto g_yy_0_0_0_y_xz_y_xx = buffer_2000_pdpd[1122];

    auto g_yy_0_0_0_y_xz_y_xy = buffer_2000_pdpd[1123];

    auto g_yy_0_0_0_y_xz_y_xz = buffer_2000_pdpd[1124];

    auto g_yy_0_0_0_y_xz_y_yy = buffer_2000_pdpd[1125];

    auto g_yy_0_0_0_y_xz_y_yz = buffer_2000_pdpd[1126];

    auto g_yy_0_0_0_y_xz_y_zz = buffer_2000_pdpd[1127];

    auto g_yy_0_0_0_y_xz_z_xx = buffer_2000_pdpd[1128];

    auto g_yy_0_0_0_y_xz_z_xy = buffer_2000_pdpd[1129];

    auto g_yy_0_0_0_y_xz_z_xz = buffer_2000_pdpd[1130];

    auto g_yy_0_0_0_y_xz_z_yy = buffer_2000_pdpd[1131];

    auto g_yy_0_0_0_y_xz_z_yz = buffer_2000_pdpd[1132];

    auto g_yy_0_0_0_y_xz_z_zz = buffer_2000_pdpd[1133];

    auto g_yy_0_0_0_y_yy_x_xx = buffer_2000_pdpd[1134];

    auto g_yy_0_0_0_y_yy_x_xy = buffer_2000_pdpd[1135];

    auto g_yy_0_0_0_y_yy_x_xz = buffer_2000_pdpd[1136];

    auto g_yy_0_0_0_y_yy_x_yy = buffer_2000_pdpd[1137];

    auto g_yy_0_0_0_y_yy_x_yz = buffer_2000_pdpd[1138];

    auto g_yy_0_0_0_y_yy_x_zz = buffer_2000_pdpd[1139];

    auto g_yy_0_0_0_y_yy_y_xx = buffer_2000_pdpd[1140];

    auto g_yy_0_0_0_y_yy_y_xy = buffer_2000_pdpd[1141];

    auto g_yy_0_0_0_y_yy_y_xz = buffer_2000_pdpd[1142];

    auto g_yy_0_0_0_y_yy_y_yy = buffer_2000_pdpd[1143];

    auto g_yy_0_0_0_y_yy_y_yz = buffer_2000_pdpd[1144];

    auto g_yy_0_0_0_y_yy_y_zz = buffer_2000_pdpd[1145];

    auto g_yy_0_0_0_y_yy_z_xx = buffer_2000_pdpd[1146];

    auto g_yy_0_0_0_y_yy_z_xy = buffer_2000_pdpd[1147];

    auto g_yy_0_0_0_y_yy_z_xz = buffer_2000_pdpd[1148];

    auto g_yy_0_0_0_y_yy_z_yy = buffer_2000_pdpd[1149];

    auto g_yy_0_0_0_y_yy_z_yz = buffer_2000_pdpd[1150];

    auto g_yy_0_0_0_y_yy_z_zz = buffer_2000_pdpd[1151];

    auto g_yy_0_0_0_y_yz_x_xx = buffer_2000_pdpd[1152];

    auto g_yy_0_0_0_y_yz_x_xy = buffer_2000_pdpd[1153];

    auto g_yy_0_0_0_y_yz_x_xz = buffer_2000_pdpd[1154];

    auto g_yy_0_0_0_y_yz_x_yy = buffer_2000_pdpd[1155];

    auto g_yy_0_0_0_y_yz_x_yz = buffer_2000_pdpd[1156];

    auto g_yy_0_0_0_y_yz_x_zz = buffer_2000_pdpd[1157];

    auto g_yy_0_0_0_y_yz_y_xx = buffer_2000_pdpd[1158];

    auto g_yy_0_0_0_y_yz_y_xy = buffer_2000_pdpd[1159];

    auto g_yy_0_0_0_y_yz_y_xz = buffer_2000_pdpd[1160];

    auto g_yy_0_0_0_y_yz_y_yy = buffer_2000_pdpd[1161];

    auto g_yy_0_0_0_y_yz_y_yz = buffer_2000_pdpd[1162];

    auto g_yy_0_0_0_y_yz_y_zz = buffer_2000_pdpd[1163];

    auto g_yy_0_0_0_y_yz_z_xx = buffer_2000_pdpd[1164];

    auto g_yy_0_0_0_y_yz_z_xy = buffer_2000_pdpd[1165];

    auto g_yy_0_0_0_y_yz_z_xz = buffer_2000_pdpd[1166];

    auto g_yy_0_0_0_y_yz_z_yy = buffer_2000_pdpd[1167];

    auto g_yy_0_0_0_y_yz_z_yz = buffer_2000_pdpd[1168];

    auto g_yy_0_0_0_y_yz_z_zz = buffer_2000_pdpd[1169];

    auto g_yy_0_0_0_y_zz_x_xx = buffer_2000_pdpd[1170];

    auto g_yy_0_0_0_y_zz_x_xy = buffer_2000_pdpd[1171];

    auto g_yy_0_0_0_y_zz_x_xz = buffer_2000_pdpd[1172];

    auto g_yy_0_0_0_y_zz_x_yy = buffer_2000_pdpd[1173];

    auto g_yy_0_0_0_y_zz_x_yz = buffer_2000_pdpd[1174];

    auto g_yy_0_0_0_y_zz_x_zz = buffer_2000_pdpd[1175];

    auto g_yy_0_0_0_y_zz_y_xx = buffer_2000_pdpd[1176];

    auto g_yy_0_0_0_y_zz_y_xy = buffer_2000_pdpd[1177];

    auto g_yy_0_0_0_y_zz_y_xz = buffer_2000_pdpd[1178];

    auto g_yy_0_0_0_y_zz_y_yy = buffer_2000_pdpd[1179];

    auto g_yy_0_0_0_y_zz_y_yz = buffer_2000_pdpd[1180];

    auto g_yy_0_0_0_y_zz_y_zz = buffer_2000_pdpd[1181];

    auto g_yy_0_0_0_y_zz_z_xx = buffer_2000_pdpd[1182];

    auto g_yy_0_0_0_y_zz_z_xy = buffer_2000_pdpd[1183];

    auto g_yy_0_0_0_y_zz_z_xz = buffer_2000_pdpd[1184];

    auto g_yy_0_0_0_y_zz_z_yy = buffer_2000_pdpd[1185];

    auto g_yy_0_0_0_y_zz_z_yz = buffer_2000_pdpd[1186];

    auto g_yy_0_0_0_y_zz_z_zz = buffer_2000_pdpd[1187];

    auto g_yy_0_0_0_z_xx_x_xx = buffer_2000_pdpd[1188];

    auto g_yy_0_0_0_z_xx_x_xy = buffer_2000_pdpd[1189];

    auto g_yy_0_0_0_z_xx_x_xz = buffer_2000_pdpd[1190];

    auto g_yy_0_0_0_z_xx_x_yy = buffer_2000_pdpd[1191];

    auto g_yy_0_0_0_z_xx_x_yz = buffer_2000_pdpd[1192];

    auto g_yy_0_0_0_z_xx_x_zz = buffer_2000_pdpd[1193];

    auto g_yy_0_0_0_z_xx_y_xx = buffer_2000_pdpd[1194];

    auto g_yy_0_0_0_z_xx_y_xy = buffer_2000_pdpd[1195];

    auto g_yy_0_0_0_z_xx_y_xz = buffer_2000_pdpd[1196];

    auto g_yy_0_0_0_z_xx_y_yy = buffer_2000_pdpd[1197];

    auto g_yy_0_0_0_z_xx_y_yz = buffer_2000_pdpd[1198];

    auto g_yy_0_0_0_z_xx_y_zz = buffer_2000_pdpd[1199];

    auto g_yy_0_0_0_z_xx_z_xx = buffer_2000_pdpd[1200];

    auto g_yy_0_0_0_z_xx_z_xy = buffer_2000_pdpd[1201];

    auto g_yy_0_0_0_z_xx_z_xz = buffer_2000_pdpd[1202];

    auto g_yy_0_0_0_z_xx_z_yy = buffer_2000_pdpd[1203];

    auto g_yy_0_0_0_z_xx_z_yz = buffer_2000_pdpd[1204];

    auto g_yy_0_0_0_z_xx_z_zz = buffer_2000_pdpd[1205];

    auto g_yy_0_0_0_z_xy_x_xx = buffer_2000_pdpd[1206];

    auto g_yy_0_0_0_z_xy_x_xy = buffer_2000_pdpd[1207];

    auto g_yy_0_0_0_z_xy_x_xz = buffer_2000_pdpd[1208];

    auto g_yy_0_0_0_z_xy_x_yy = buffer_2000_pdpd[1209];

    auto g_yy_0_0_0_z_xy_x_yz = buffer_2000_pdpd[1210];

    auto g_yy_0_0_0_z_xy_x_zz = buffer_2000_pdpd[1211];

    auto g_yy_0_0_0_z_xy_y_xx = buffer_2000_pdpd[1212];

    auto g_yy_0_0_0_z_xy_y_xy = buffer_2000_pdpd[1213];

    auto g_yy_0_0_0_z_xy_y_xz = buffer_2000_pdpd[1214];

    auto g_yy_0_0_0_z_xy_y_yy = buffer_2000_pdpd[1215];

    auto g_yy_0_0_0_z_xy_y_yz = buffer_2000_pdpd[1216];

    auto g_yy_0_0_0_z_xy_y_zz = buffer_2000_pdpd[1217];

    auto g_yy_0_0_0_z_xy_z_xx = buffer_2000_pdpd[1218];

    auto g_yy_0_0_0_z_xy_z_xy = buffer_2000_pdpd[1219];

    auto g_yy_0_0_0_z_xy_z_xz = buffer_2000_pdpd[1220];

    auto g_yy_0_0_0_z_xy_z_yy = buffer_2000_pdpd[1221];

    auto g_yy_0_0_0_z_xy_z_yz = buffer_2000_pdpd[1222];

    auto g_yy_0_0_0_z_xy_z_zz = buffer_2000_pdpd[1223];

    auto g_yy_0_0_0_z_xz_x_xx = buffer_2000_pdpd[1224];

    auto g_yy_0_0_0_z_xz_x_xy = buffer_2000_pdpd[1225];

    auto g_yy_0_0_0_z_xz_x_xz = buffer_2000_pdpd[1226];

    auto g_yy_0_0_0_z_xz_x_yy = buffer_2000_pdpd[1227];

    auto g_yy_0_0_0_z_xz_x_yz = buffer_2000_pdpd[1228];

    auto g_yy_0_0_0_z_xz_x_zz = buffer_2000_pdpd[1229];

    auto g_yy_0_0_0_z_xz_y_xx = buffer_2000_pdpd[1230];

    auto g_yy_0_0_0_z_xz_y_xy = buffer_2000_pdpd[1231];

    auto g_yy_0_0_0_z_xz_y_xz = buffer_2000_pdpd[1232];

    auto g_yy_0_0_0_z_xz_y_yy = buffer_2000_pdpd[1233];

    auto g_yy_0_0_0_z_xz_y_yz = buffer_2000_pdpd[1234];

    auto g_yy_0_0_0_z_xz_y_zz = buffer_2000_pdpd[1235];

    auto g_yy_0_0_0_z_xz_z_xx = buffer_2000_pdpd[1236];

    auto g_yy_0_0_0_z_xz_z_xy = buffer_2000_pdpd[1237];

    auto g_yy_0_0_0_z_xz_z_xz = buffer_2000_pdpd[1238];

    auto g_yy_0_0_0_z_xz_z_yy = buffer_2000_pdpd[1239];

    auto g_yy_0_0_0_z_xz_z_yz = buffer_2000_pdpd[1240];

    auto g_yy_0_0_0_z_xz_z_zz = buffer_2000_pdpd[1241];

    auto g_yy_0_0_0_z_yy_x_xx = buffer_2000_pdpd[1242];

    auto g_yy_0_0_0_z_yy_x_xy = buffer_2000_pdpd[1243];

    auto g_yy_0_0_0_z_yy_x_xz = buffer_2000_pdpd[1244];

    auto g_yy_0_0_0_z_yy_x_yy = buffer_2000_pdpd[1245];

    auto g_yy_0_0_0_z_yy_x_yz = buffer_2000_pdpd[1246];

    auto g_yy_0_0_0_z_yy_x_zz = buffer_2000_pdpd[1247];

    auto g_yy_0_0_0_z_yy_y_xx = buffer_2000_pdpd[1248];

    auto g_yy_0_0_0_z_yy_y_xy = buffer_2000_pdpd[1249];

    auto g_yy_0_0_0_z_yy_y_xz = buffer_2000_pdpd[1250];

    auto g_yy_0_0_0_z_yy_y_yy = buffer_2000_pdpd[1251];

    auto g_yy_0_0_0_z_yy_y_yz = buffer_2000_pdpd[1252];

    auto g_yy_0_0_0_z_yy_y_zz = buffer_2000_pdpd[1253];

    auto g_yy_0_0_0_z_yy_z_xx = buffer_2000_pdpd[1254];

    auto g_yy_0_0_0_z_yy_z_xy = buffer_2000_pdpd[1255];

    auto g_yy_0_0_0_z_yy_z_xz = buffer_2000_pdpd[1256];

    auto g_yy_0_0_0_z_yy_z_yy = buffer_2000_pdpd[1257];

    auto g_yy_0_0_0_z_yy_z_yz = buffer_2000_pdpd[1258];

    auto g_yy_0_0_0_z_yy_z_zz = buffer_2000_pdpd[1259];

    auto g_yy_0_0_0_z_yz_x_xx = buffer_2000_pdpd[1260];

    auto g_yy_0_0_0_z_yz_x_xy = buffer_2000_pdpd[1261];

    auto g_yy_0_0_0_z_yz_x_xz = buffer_2000_pdpd[1262];

    auto g_yy_0_0_0_z_yz_x_yy = buffer_2000_pdpd[1263];

    auto g_yy_0_0_0_z_yz_x_yz = buffer_2000_pdpd[1264];

    auto g_yy_0_0_0_z_yz_x_zz = buffer_2000_pdpd[1265];

    auto g_yy_0_0_0_z_yz_y_xx = buffer_2000_pdpd[1266];

    auto g_yy_0_0_0_z_yz_y_xy = buffer_2000_pdpd[1267];

    auto g_yy_0_0_0_z_yz_y_xz = buffer_2000_pdpd[1268];

    auto g_yy_0_0_0_z_yz_y_yy = buffer_2000_pdpd[1269];

    auto g_yy_0_0_0_z_yz_y_yz = buffer_2000_pdpd[1270];

    auto g_yy_0_0_0_z_yz_y_zz = buffer_2000_pdpd[1271];

    auto g_yy_0_0_0_z_yz_z_xx = buffer_2000_pdpd[1272];

    auto g_yy_0_0_0_z_yz_z_xy = buffer_2000_pdpd[1273];

    auto g_yy_0_0_0_z_yz_z_xz = buffer_2000_pdpd[1274];

    auto g_yy_0_0_0_z_yz_z_yy = buffer_2000_pdpd[1275];

    auto g_yy_0_0_0_z_yz_z_yz = buffer_2000_pdpd[1276];

    auto g_yy_0_0_0_z_yz_z_zz = buffer_2000_pdpd[1277];

    auto g_yy_0_0_0_z_zz_x_xx = buffer_2000_pdpd[1278];

    auto g_yy_0_0_0_z_zz_x_xy = buffer_2000_pdpd[1279];

    auto g_yy_0_0_0_z_zz_x_xz = buffer_2000_pdpd[1280];

    auto g_yy_0_0_0_z_zz_x_yy = buffer_2000_pdpd[1281];

    auto g_yy_0_0_0_z_zz_x_yz = buffer_2000_pdpd[1282];

    auto g_yy_0_0_0_z_zz_x_zz = buffer_2000_pdpd[1283];

    auto g_yy_0_0_0_z_zz_y_xx = buffer_2000_pdpd[1284];

    auto g_yy_0_0_0_z_zz_y_xy = buffer_2000_pdpd[1285];

    auto g_yy_0_0_0_z_zz_y_xz = buffer_2000_pdpd[1286];

    auto g_yy_0_0_0_z_zz_y_yy = buffer_2000_pdpd[1287];

    auto g_yy_0_0_0_z_zz_y_yz = buffer_2000_pdpd[1288];

    auto g_yy_0_0_0_z_zz_y_zz = buffer_2000_pdpd[1289];

    auto g_yy_0_0_0_z_zz_z_xx = buffer_2000_pdpd[1290];

    auto g_yy_0_0_0_z_zz_z_xy = buffer_2000_pdpd[1291];

    auto g_yy_0_0_0_z_zz_z_xz = buffer_2000_pdpd[1292];

    auto g_yy_0_0_0_z_zz_z_yy = buffer_2000_pdpd[1293];

    auto g_yy_0_0_0_z_zz_z_yz = buffer_2000_pdpd[1294];

    auto g_yy_0_0_0_z_zz_z_zz = buffer_2000_pdpd[1295];

    auto g_yz_0_0_0_x_xx_x_xx = buffer_2000_pdpd[1296];

    auto g_yz_0_0_0_x_xx_x_xy = buffer_2000_pdpd[1297];

    auto g_yz_0_0_0_x_xx_x_xz = buffer_2000_pdpd[1298];

    auto g_yz_0_0_0_x_xx_x_yy = buffer_2000_pdpd[1299];

    auto g_yz_0_0_0_x_xx_x_yz = buffer_2000_pdpd[1300];

    auto g_yz_0_0_0_x_xx_x_zz = buffer_2000_pdpd[1301];

    auto g_yz_0_0_0_x_xx_y_xx = buffer_2000_pdpd[1302];

    auto g_yz_0_0_0_x_xx_y_xy = buffer_2000_pdpd[1303];

    auto g_yz_0_0_0_x_xx_y_xz = buffer_2000_pdpd[1304];

    auto g_yz_0_0_0_x_xx_y_yy = buffer_2000_pdpd[1305];

    auto g_yz_0_0_0_x_xx_y_yz = buffer_2000_pdpd[1306];

    auto g_yz_0_0_0_x_xx_y_zz = buffer_2000_pdpd[1307];

    auto g_yz_0_0_0_x_xx_z_xx = buffer_2000_pdpd[1308];

    auto g_yz_0_0_0_x_xx_z_xy = buffer_2000_pdpd[1309];

    auto g_yz_0_0_0_x_xx_z_xz = buffer_2000_pdpd[1310];

    auto g_yz_0_0_0_x_xx_z_yy = buffer_2000_pdpd[1311];

    auto g_yz_0_0_0_x_xx_z_yz = buffer_2000_pdpd[1312];

    auto g_yz_0_0_0_x_xx_z_zz = buffer_2000_pdpd[1313];

    auto g_yz_0_0_0_x_xy_x_xx = buffer_2000_pdpd[1314];

    auto g_yz_0_0_0_x_xy_x_xy = buffer_2000_pdpd[1315];

    auto g_yz_0_0_0_x_xy_x_xz = buffer_2000_pdpd[1316];

    auto g_yz_0_0_0_x_xy_x_yy = buffer_2000_pdpd[1317];

    auto g_yz_0_0_0_x_xy_x_yz = buffer_2000_pdpd[1318];

    auto g_yz_0_0_0_x_xy_x_zz = buffer_2000_pdpd[1319];

    auto g_yz_0_0_0_x_xy_y_xx = buffer_2000_pdpd[1320];

    auto g_yz_0_0_0_x_xy_y_xy = buffer_2000_pdpd[1321];

    auto g_yz_0_0_0_x_xy_y_xz = buffer_2000_pdpd[1322];

    auto g_yz_0_0_0_x_xy_y_yy = buffer_2000_pdpd[1323];

    auto g_yz_0_0_0_x_xy_y_yz = buffer_2000_pdpd[1324];

    auto g_yz_0_0_0_x_xy_y_zz = buffer_2000_pdpd[1325];

    auto g_yz_0_0_0_x_xy_z_xx = buffer_2000_pdpd[1326];

    auto g_yz_0_0_0_x_xy_z_xy = buffer_2000_pdpd[1327];

    auto g_yz_0_0_0_x_xy_z_xz = buffer_2000_pdpd[1328];

    auto g_yz_0_0_0_x_xy_z_yy = buffer_2000_pdpd[1329];

    auto g_yz_0_0_0_x_xy_z_yz = buffer_2000_pdpd[1330];

    auto g_yz_0_0_0_x_xy_z_zz = buffer_2000_pdpd[1331];

    auto g_yz_0_0_0_x_xz_x_xx = buffer_2000_pdpd[1332];

    auto g_yz_0_0_0_x_xz_x_xy = buffer_2000_pdpd[1333];

    auto g_yz_0_0_0_x_xz_x_xz = buffer_2000_pdpd[1334];

    auto g_yz_0_0_0_x_xz_x_yy = buffer_2000_pdpd[1335];

    auto g_yz_0_0_0_x_xz_x_yz = buffer_2000_pdpd[1336];

    auto g_yz_0_0_0_x_xz_x_zz = buffer_2000_pdpd[1337];

    auto g_yz_0_0_0_x_xz_y_xx = buffer_2000_pdpd[1338];

    auto g_yz_0_0_0_x_xz_y_xy = buffer_2000_pdpd[1339];

    auto g_yz_0_0_0_x_xz_y_xz = buffer_2000_pdpd[1340];

    auto g_yz_0_0_0_x_xz_y_yy = buffer_2000_pdpd[1341];

    auto g_yz_0_0_0_x_xz_y_yz = buffer_2000_pdpd[1342];

    auto g_yz_0_0_0_x_xz_y_zz = buffer_2000_pdpd[1343];

    auto g_yz_0_0_0_x_xz_z_xx = buffer_2000_pdpd[1344];

    auto g_yz_0_0_0_x_xz_z_xy = buffer_2000_pdpd[1345];

    auto g_yz_0_0_0_x_xz_z_xz = buffer_2000_pdpd[1346];

    auto g_yz_0_0_0_x_xz_z_yy = buffer_2000_pdpd[1347];

    auto g_yz_0_0_0_x_xz_z_yz = buffer_2000_pdpd[1348];

    auto g_yz_0_0_0_x_xz_z_zz = buffer_2000_pdpd[1349];

    auto g_yz_0_0_0_x_yy_x_xx = buffer_2000_pdpd[1350];

    auto g_yz_0_0_0_x_yy_x_xy = buffer_2000_pdpd[1351];

    auto g_yz_0_0_0_x_yy_x_xz = buffer_2000_pdpd[1352];

    auto g_yz_0_0_0_x_yy_x_yy = buffer_2000_pdpd[1353];

    auto g_yz_0_0_0_x_yy_x_yz = buffer_2000_pdpd[1354];

    auto g_yz_0_0_0_x_yy_x_zz = buffer_2000_pdpd[1355];

    auto g_yz_0_0_0_x_yy_y_xx = buffer_2000_pdpd[1356];

    auto g_yz_0_0_0_x_yy_y_xy = buffer_2000_pdpd[1357];

    auto g_yz_0_0_0_x_yy_y_xz = buffer_2000_pdpd[1358];

    auto g_yz_0_0_0_x_yy_y_yy = buffer_2000_pdpd[1359];

    auto g_yz_0_0_0_x_yy_y_yz = buffer_2000_pdpd[1360];

    auto g_yz_0_0_0_x_yy_y_zz = buffer_2000_pdpd[1361];

    auto g_yz_0_0_0_x_yy_z_xx = buffer_2000_pdpd[1362];

    auto g_yz_0_0_0_x_yy_z_xy = buffer_2000_pdpd[1363];

    auto g_yz_0_0_0_x_yy_z_xz = buffer_2000_pdpd[1364];

    auto g_yz_0_0_0_x_yy_z_yy = buffer_2000_pdpd[1365];

    auto g_yz_0_0_0_x_yy_z_yz = buffer_2000_pdpd[1366];

    auto g_yz_0_0_0_x_yy_z_zz = buffer_2000_pdpd[1367];

    auto g_yz_0_0_0_x_yz_x_xx = buffer_2000_pdpd[1368];

    auto g_yz_0_0_0_x_yz_x_xy = buffer_2000_pdpd[1369];

    auto g_yz_0_0_0_x_yz_x_xz = buffer_2000_pdpd[1370];

    auto g_yz_0_0_0_x_yz_x_yy = buffer_2000_pdpd[1371];

    auto g_yz_0_0_0_x_yz_x_yz = buffer_2000_pdpd[1372];

    auto g_yz_0_0_0_x_yz_x_zz = buffer_2000_pdpd[1373];

    auto g_yz_0_0_0_x_yz_y_xx = buffer_2000_pdpd[1374];

    auto g_yz_0_0_0_x_yz_y_xy = buffer_2000_pdpd[1375];

    auto g_yz_0_0_0_x_yz_y_xz = buffer_2000_pdpd[1376];

    auto g_yz_0_0_0_x_yz_y_yy = buffer_2000_pdpd[1377];

    auto g_yz_0_0_0_x_yz_y_yz = buffer_2000_pdpd[1378];

    auto g_yz_0_0_0_x_yz_y_zz = buffer_2000_pdpd[1379];

    auto g_yz_0_0_0_x_yz_z_xx = buffer_2000_pdpd[1380];

    auto g_yz_0_0_0_x_yz_z_xy = buffer_2000_pdpd[1381];

    auto g_yz_0_0_0_x_yz_z_xz = buffer_2000_pdpd[1382];

    auto g_yz_0_0_0_x_yz_z_yy = buffer_2000_pdpd[1383];

    auto g_yz_0_0_0_x_yz_z_yz = buffer_2000_pdpd[1384];

    auto g_yz_0_0_0_x_yz_z_zz = buffer_2000_pdpd[1385];

    auto g_yz_0_0_0_x_zz_x_xx = buffer_2000_pdpd[1386];

    auto g_yz_0_0_0_x_zz_x_xy = buffer_2000_pdpd[1387];

    auto g_yz_0_0_0_x_zz_x_xz = buffer_2000_pdpd[1388];

    auto g_yz_0_0_0_x_zz_x_yy = buffer_2000_pdpd[1389];

    auto g_yz_0_0_0_x_zz_x_yz = buffer_2000_pdpd[1390];

    auto g_yz_0_0_0_x_zz_x_zz = buffer_2000_pdpd[1391];

    auto g_yz_0_0_0_x_zz_y_xx = buffer_2000_pdpd[1392];

    auto g_yz_0_0_0_x_zz_y_xy = buffer_2000_pdpd[1393];

    auto g_yz_0_0_0_x_zz_y_xz = buffer_2000_pdpd[1394];

    auto g_yz_0_0_0_x_zz_y_yy = buffer_2000_pdpd[1395];

    auto g_yz_0_0_0_x_zz_y_yz = buffer_2000_pdpd[1396];

    auto g_yz_0_0_0_x_zz_y_zz = buffer_2000_pdpd[1397];

    auto g_yz_0_0_0_x_zz_z_xx = buffer_2000_pdpd[1398];

    auto g_yz_0_0_0_x_zz_z_xy = buffer_2000_pdpd[1399];

    auto g_yz_0_0_0_x_zz_z_xz = buffer_2000_pdpd[1400];

    auto g_yz_0_0_0_x_zz_z_yy = buffer_2000_pdpd[1401];

    auto g_yz_0_0_0_x_zz_z_yz = buffer_2000_pdpd[1402];

    auto g_yz_0_0_0_x_zz_z_zz = buffer_2000_pdpd[1403];

    auto g_yz_0_0_0_y_xx_x_xx = buffer_2000_pdpd[1404];

    auto g_yz_0_0_0_y_xx_x_xy = buffer_2000_pdpd[1405];

    auto g_yz_0_0_0_y_xx_x_xz = buffer_2000_pdpd[1406];

    auto g_yz_0_0_0_y_xx_x_yy = buffer_2000_pdpd[1407];

    auto g_yz_0_0_0_y_xx_x_yz = buffer_2000_pdpd[1408];

    auto g_yz_0_0_0_y_xx_x_zz = buffer_2000_pdpd[1409];

    auto g_yz_0_0_0_y_xx_y_xx = buffer_2000_pdpd[1410];

    auto g_yz_0_0_0_y_xx_y_xy = buffer_2000_pdpd[1411];

    auto g_yz_0_0_0_y_xx_y_xz = buffer_2000_pdpd[1412];

    auto g_yz_0_0_0_y_xx_y_yy = buffer_2000_pdpd[1413];

    auto g_yz_0_0_0_y_xx_y_yz = buffer_2000_pdpd[1414];

    auto g_yz_0_0_0_y_xx_y_zz = buffer_2000_pdpd[1415];

    auto g_yz_0_0_0_y_xx_z_xx = buffer_2000_pdpd[1416];

    auto g_yz_0_0_0_y_xx_z_xy = buffer_2000_pdpd[1417];

    auto g_yz_0_0_0_y_xx_z_xz = buffer_2000_pdpd[1418];

    auto g_yz_0_0_0_y_xx_z_yy = buffer_2000_pdpd[1419];

    auto g_yz_0_0_0_y_xx_z_yz = buffer_2000_pdpd[1420];

    auto g_yz_0_0_0_y_xx_z_zz = buffer_2000_pdpd[1421];

    auto g_yz_0_0_0_y_xy_x_xx = buffer_2000_pdpd[1422];

    auto g_yz_0_0_0_y_xy_x_xy = buffer_2000_pdpd[1423];

    auto g_yz_0_0_0_y_xy_x_xz = buffer_2000_pdpd[1424];

    auto g_yz_0_0_0_y_xy_x_yy = buffer_2000_pdpd[1425];

    auto g_yz_0_0_0_y_xy_x_yz = buffer_2000_pdpd[1426];

    auto g_yz_0_0_0_y_xy_x_zz = buffer_2000_pdpd[1427];

    auto g_yz_0_0_0_y_xy_y_xx = buffer_2000_pdpd[1428];

    auto g_yz_0_0_0_y_xy_y_xy = buffer_2000_pdpd[1429];

    auto g_yz_0_0_0_y_xy_y_xz = buffer_2000_pdpd[1430];

    auto g_yz_0_0_0_y_xy_y_yy = buffer_2000_pdpd[1431];

    auto g_yz_0_0_0_y_xy_y_yz = buffer_2000_pdpd[1432];

    auto g_yz_0_0_0_y_xy_y_zz = buffer_2000_pdpd[1433];

    auto g_yz_0_0_0_y_xy_z_xx = buffer_2000_pdpd[1434];

    auto g_yz_0_0_0_y_xy_z_xy = buffer_2000_pdpd[1435];

    auto g_yz_0_0_0_y_xy_z_xz = buffer_2000_pdpd[1436];

    auto g_yz_0_0_0_y_xy_z_yy = buffer_2000_pdpd[1437];

    auto g_yz_0_0_0_y_xy_z_yz = buffer_2000_pdpd[1438];

    auto g_yz_0_0_0_y_xy_z_zz = buffer_2000_pdpd[1439];

    auto g_yz_0_0_0_y_xz_x_xx = buffer_2000_pdpd[1440];

    auto g_yz_0_0_0_y_xz_x_xy = buffer_2000_pdpd[1441];

    auto g_yz_0_0_0_y_xz_x_xz = buffer_2000_pdpd[1442];

    auto g_yz_0_0_0_y_xz_x_yy = buffer_2000_pdpd[1443];

    auto g_yz_0_0_0_y_xz_x_yz = buffer_2000_pdpd[1444];

    auto g_yz_0_0_0_y_xz_x_zz = buffer_2000_pdpd[1445];

    auto g_yz_0_0_0_y_xz_y_xx = buffer_2000_pdpd[1446];

    auto g_yz_0_0_0_y_xz_y_xy = buffer_2000_pdpd[1447];

    auto g_yz_0_0_0_y_xz_y_xz = buffer_2000_pdpd[1448];

    auto g_yz_0_0_0_y_xz_y_yy = buffer_2000_pdpd[1449];

    auto g_yz_0_0_0_y_xz_y_yz = buffer_2000_pdpd[1450];

    auto g_yz_0_0_0_y_xz_y_zz = buffer_2000_pdpd[1451];

    auto g_yz_0_0_0_y_xz_z_xx = buffer_2000_pdpd[1452];

    auto g_yz_0_0_0_y_xz_z_xy = buffer_2000_pdpd[1453];

    auto g_yz_0_0_0_y_xz_z_xz = buffer_2000_pdpd[1454];

    auto g_yz_0_0_0_y_xz_z_yy = buffer_2000_pdpd[1455];

    auto g_yz_0_0_0_y_xz_z_yz = buffer_2000_pdpd[1456];

    auto g_yz_0_0_0_y_xz_z_zz = buffer_2000_pdpd[1457];

    auto g_yz_0_0_0_y_yy_x_xx = buffer_2000_pdpd[1458];

    auto g_yz_0_0_0_y_yy_x_xy = buffer_2000_pdpd[1459];

    auto g_yz_0_0_0_y_yy_x_xz = buffer_2000_pdpd[1460];

    auto g_yz_0_0_0_y_yy_x_yy = buffer_2000_pdpd[1461];

    auto g_yz_0_0_0_y_yy_x_yz = buffer_2000_pdpd[1462];

    auto g_yz_0_0_0_y_yy_x_zz = buffer_2000_pdpd[1463];

    auto g_yz_0_0_0_y_yy_y_xx = buffer_2000_pdpd[1464];

    auto g_yz_0_0_0_y_yy_y_xy = buffer_2000_pdpd[1465];

    auto g_yz_0_0_0_y_yy_y_xz = buffer_2000_pdpd[1466];

    auto g_yz_0_0_0_y_yy_y_yy = buffer_2000_pdpd[1467];

    auto g_yz_0_0_0_y_yy_y_yz = buffer_2000_pdpd[1468];

    auto g_yz_0_0_0_y_yy_y_zz = buffer_2000_pdpd[1469];

    auto g_yz_0_0_0_y_yy_z_xx = buffer_2000_pdpd[1470];

    auto g_yz_0_0_0_y_yy_z_xy = buffer_2000_pdpd[1471];

    auto g_yz_0_0_0_y_yy_z_xz = buffer_2000_pdpd[1472];

    auto g_yz_0_0_0_y_yy_z_yy = buffer_2000_pdpd[1473];

    auto g_yz_0_0_0_y_yy_z_yz = buffer_2000_pdpd[1474];

    auto g_yz_0_0_0_y_yy_z_zz = buffer_2000_pdpd[1475];

    auto g_yz_0_0_0_y_yz_x_xx = buffer_2000_pdpd[1476];

    auto g_yz_0_0_0_y_yz_x_xy = buffer_2000_pdpd[1477];

    auto g_yz_0_0_0_y_yz_x_xz = buffer_2000_pdpd[1478];

    auto g_yz_0_0_0_y_yz_x_yy = buffer_2000_pdpd[1479];

    auto g_yz_0_0_0_y_yz_x_yz = buffer_2000_pdpd[1480];

    auto g_yz_0_0_0_y_yz_x_zz = buffer_2000_pdpd[1481];

    auto g_yz_0_0_0_y_yz_y_xx = buffer_2000_pdpd[1482];

    auto g_yz_0_0_0_y_yz_y_xy = buffer_2000_pdpd[1483];

    auto g_yz_0_0_0_y_yz_y_xz = buffer_2000_pdpd[1484];

    auto g_yz_0_0_0_y_yz_y_yy = buffer_2000_pdpd[1485];

    auto g_yz_0_0_0_y_yz_y_yz = buffer_2000_pdpd[1486];

    auto g_yz_0_0_0_y_yz_y_zz = buffer_2000_pdpd[1487];

    auto g_yz_0_0_0_y_yz_z_xx = buffer_2000_pdpd[1488];

    auto g_yz_0_0_0_y_yz_z_xy = buffer_2000_pdpd[1489];

    auto g_yz_0_0_0_y_yz_z_xz = buffer_2000_pdpd[1490];

    auto g_yz_0_0_0_y_yz_z_yy = buffer_2000_pdpd[1491];

    auto g_yz_0_0_0_y_yz_z_yz = buffer_2000_pdpd[1492];

    auto g_yz_0_0_0_y_yz_z_zz = buffer_2000_pdpd[1493];

    auto g_yz_0_0_0_y_zz_x_xx = buffer_2000_pdpd[1494];

    auto g_yz_0_0_0_y_zz_x_xy = buffer_2000_pdpd[1495];

    auto g_yz_0_0_0_y_zz_x_xz = buffer_2000_pdpd[1496];

    auto g_yz_0_0_0_y_zz_x_yy = buffer_2000_pdpd[1497];

    auto g_yz_0_0_0_y_zz_x_yz = buffer_2000_pdpd[1498];

    auto g_yz_0_0_0_y_zz_x_zz = buffer_2000_pdpd[1499];

    auto g_yz_0_0_0_y_zz_y_xx = buffer_2000_pdpd[1500];

    auto g_yz_0_0_0_y_zz_y_xy = buffer_2000_pdpd[1501];

    auto g_yz_0_0_0_y_zz_y_xz = buffer_2000_pdpd[1502];

    auto g_yz_0_0_0_y_zz_y_yy = buffer_2000_pdpd[1503];

    auto g_yz_0_0_0_y_zz_y_yz = buffer_2000_pdpd[1504];

    auto g_yz_0_0_0_y_zz_y_zz = buffer_2000_pdpd[1505];

    auto g_yz_0_0_0_y_zz_z_xx = buffer_2000_pdpd[1506];

    auto g_yz_0_0_0_y_zz_z_xy = buffer_2000_pdpd[1507];

    auto g_yz_0_0_0_y_zz_z_xz = buffer_2000_pdpd[1508];

    auto g_yz_0_0_0_y_zz_z_yy = buffer_2000_pdpd[1509];

    auto g_yz_0_0_0_y_zz_z_yz = buffer_2000_pdpd[1510];

    auto g_yz_0_0_0_y_zz_z_zz = buffer_2000_pdpd[1511];

    auto g_yz_0_0_0_z_xx_x_xx = buffer_2000_pdpd[1512];

    auto g_yz_0_0_0_z_xx_x_xy = buffer_2000_pdpd[1513];

    auto g_yz_0_0_0_z_xx_x_xz = buffer_2000_pdpd[1514];

    auto g_yz_0_0_0_z_xx_x_yy = buffer_2000_pdpd[1515];

    auto g_yz_0_0_0_z_xx_x_yz = buffer_2000_pdpd[1516];

    auto g_yz_0_0_0_z_xx_x_zz = buffer_2000_pdpd[1517];

    auto g_yz_0_0_0_z_xx_y_xx = buffer_2000_pdpd[1518];

    auto g_yz_0_0_0_z_xx_y_xy = buffer_2000_pdpd[1519];

    auto g_yz_0_0_0_z_xx_y_xz = buffer_2000_pdpd[1520];

    auto g_yz_0_0_0_z_xx_y_yy = buffer_2000_pdpd[1521];

    auto g_yz_0_0_0_z_xx_y_yz = buffer_2000_pdpd[1522];

    auto g_yz_0_0_0_z_xx_y_zz = buffer_2000_pdpd[1523];

    auto g_yz_0_0_0_z_xx_z_xx = buffer_2000_pdpd[1524];

    auto g_yz_0_0_0_z_xx_z_xy = buffer_2000_pdpd[1525];

    auto g_yz_0_0_0_z_xx_z_xz = buffer_2000_pdpd[1526];

    auto g_yz_0_0_0_z_xx_z_yy = buffer_2000_pdpd[1527];

    auto g_yz_0_0_0_z_xx_z_yz = buffer_2000_pdpd[1528];

    auto g_yz_0_0_0_z_xx_z_zz = buffer_2000_pdpd[1529];

    auto g_yz_0_0_0_z_xy_x_xx = buffer_2000_pdpd[1530];

    auto g_yz_0_0_0_z_xy_x_xy = buffer_2000_pdpd[1531];

    auto g_yz_0_0_0_z_xy_x_xz = buffer_2000_pdpd[1532];

    auto g_yz_0_0_0_z_xy_x_yy = buffer_2000_pdpd[1533];

    auto g_yz_0_0_0_z_xy_x_yz = buffer_2000_pdpd[1534];

    auto g_yz_0_0_0_z_xy_x_zz = buffer_2000_pdpd[1535];

    auto g_yz_0_0_0_z_xy_y_xx = buffer_2000_pdpd[1536];

    auto g_yz_0_0_0_z_xy_y_xy = buffer_2000_pdpd[1537];

    auto g_yz_0_0_0_z_xy_y_xz = buffer_2000_pdpd[1538];

    auto g_yz_0_0_0_z_xy_y_yy = buffer_2000_pdpd[1539];

    auto g_yz_0_0_0_z_xy_y_yz = buffer_2000_pdpd[1540];

    auto g_yz_0_0_0_z_xy_y_zz = buffer_2000_pdpd[1541];

    auto g_yz_0_0_0_z_xy_z_xx = buffer_2000_pdpd[1542];

    auto g_yz_0_0_0_z_xy_z_xy = buffer_2000_pdpd[1543];

    auto g_yz_0_0_0_z_xy_z_xz = buffer_2000_pdpd[1544];

    auto g_yz_0_0_0_z_xy_z_yy = buffer_2000_pdpd[1545];

    auto g_yz_0_0_0_z_xy_z_yz = buffer_2000_pdpd[1546];

    auto g_yz_0_0_0_z_xy_z_zz = buffer_2000_pdpd[1547];

    auto g_yz_0_0_0_z_xz_x_xx = buffer_2000_pdpd[1548];

    auto g_yz_0_0_0_z_xz_x_xy = buffer_2000_pdpd[1549];

    auto g_yz_0_0_0_z_xz_x_xz = buffer_2000_pdpd[1550];

    auto g_yz_0_0_0_z_xz_x_yy = buffer_2000_pdpd[1551];

    auto g_yz_0_0_0_z_xz_x_yz = buffer_2000_pdpd[1552];

    auto g_yz_0_0_0_z_xz_x_zz = buffer_2000_pdpd[1553];

    auto g_yz_0_0_0_z_xz_y_xx = buffer_2000_pdpd[1554];

    auto g_yz_0_0_0_z_xz_y_xy = buffer_2000_pdpd[1555];

    auto g_yz_0_0_0_z_xz_y_xz = buffer_2000_pdpd[1556];

    auto g_yz_0_0_0_z_xz_y_yy = buffer_2000_pdpd[1557];

    auto g_yz_0_0_0_z_xz_y_yz = buffer_2000_pdpd[1558];

    auto g_yz_0_0_0_z_xz_y_zz = buffer_2000_pdpd[1559];

    auto g_yz_0_0_0_z_xz_z_xx = buffer_2000_pdpd[1560];

    auto g_yz_0_0_0_z_xz_z_xy = buffer_2000_pdpd[1561];

    auto g_yz_0_0_0_z_xz_z_xz = buffer_2000_pdpd[1562];

    auto g_yz_0_0_0_z_xz_z_yy = buffer_2000_pdpd[1563];

    auto g_yz_0_0_0_z_xz_z_yz = buffer_2000_pdpd[1564];

    auto g_yz_0_0_0_z_xz_z_zz = buffer_2000_pdpd[1565];

    auto g_yz_0_0_0_z_yy_x_xx = buffer_2000_pdpd[1566];

    auto g_yz_0_0_0_z_yy_x_xy = buffer_2000_pdpd[1567];

    auto g_yz_0_0_0_z_yy_x_xz = buffer_2000_pdpd[1568];

    auto g_yz_0_0_0_z_yy_x_yy = buffer_2000_pdpd[1569];

    auto g_yz_0_0_0_z_yy_x_yz = buffer_2000_pdpd[1570];

    auto g_yz_0_0_0_z_yy_x_zz = buffer_2000_pdpd[1571];

    auto g_yz_0_0_0_z_yy_y_xx = buffer_2000_pdpd[1572];

    auto g_yz_0_0_0_z_yy_y_xy = buffer_2000_pdpd[1573];

    auto g_yz_0_0_0_z_yy_y_xz = buffer_2000_pdpd[1574];

    auto g_yz_0_0_0_z_yy_y_yy = buffer_2000_pdpd[1575];

    auto g_yz_0_0_0_z_yy_y_yz = buffer_2000_pdpd[1576];

    auto g_yz_0_0_0_z_yy_y_zz = buffer_2000_pdpd[1577];

    auto g_yz_0_0_0_z_yy_z_xx = buffer_2000_pdpd[1578];

    auto g_yz_0_0_0_z_yy_z_xy = buffer_2000_pdpd[1579];

    auto g_yz_0_0_0_z_yy_z_xz = buffer_2000_pdpd[1580];

    auto g_yz_0_0_0_z_yy_z_yy = buffer_2000_pdpd[1581];

    auto g_yz_0_0_0_z_yy_z_yz = buffer_2000_pdpd[1582];

    auto g_yz_0_0_0_z_yy_z_zz = buffer_2000_pdpd[1583];

    auto g_yz_0_0_0_z_yz_x_xx = buffer_2000_pdpd[1584];

    auto g_yz_0_0_0_z_yz_x_xy = buffer_2000_pdpd[1585];

    auto g_yz_0_0_0_z_yz_x_xz = buffer_2000_pdpd[1586];

    auto g_yz_0_0_0_z_yz_x_yy = buffer_2000_pdpd[1587];

    auto g_yz_0_0_0_z_yz_x_yz = buffer_2000_pdpd[1588];

    auto g_yz_0_0_0_z_yz_x_zz = buffer_2000_pdpd[1589];

    auto g_yz_0_0_0_z_yz_y_xx = buffer_2000_pdpd[1590];

    auto g_yz_0_0_0_z_yz_y_xy = buffer_2000_pdpd[1591];

    auto g_yz_0_0_0_z_yz_y_xz = buffer_2000_pdpd[1592];

    auto g_yz_0_0_0_z_yz_y_yy = buffer_2000_pdpd[1593];

    auto g_yz_0_0_0_z_yz_y_yz = buffer_2000_pdpd[1594];

    auto g_yz_0_0_0_z_yz_y_zz = buffer_2000_pdpd[1595];

    auto g_yz_0_0_0_z_yz_z_xx = buffer_2000_pdpd[1596];

    auto g_yz_0_0_0_z_yz_z_xy = buffer_2000_pdpd[1597];

    auto g_yz_0_0_0_z_yz_z_xz = buffer_2000_pdpd[1598];

    auto g_yz_0_0_0_z_yz_z_yy = buffer_2000_pdpd[1599];

    auto g_yz_0_0_0_z_yz_z_yz = buffer_2000_pdpd[1600];

    auto g_yz_0_0_0_z_yz_z_zz = buffer_2000_pdpd[1601];

    auto g_yz_0_0_0_z_zz_x_xx = buffer_2000_pdpd[1602];

    auto g_yz_0_0_0_z_zz_x_xy = buffer_2000_pdpd[1603];

    auto g_yz_0_0_0_z_zz_x_xz = buffer_2000_pdpd[1604];

    auto g_yz_0_0_0_z_zz_x_yy = buffer_2000_pdpd[1605];

    auto g_yz_0_0_0_z_zz_x_yz = buffer_2000_pdpd[1606];

    auto g_yz_0_0_0_z_zz_x_zz = buffer_2000_pdpd[1607];

    auto g_yz_0_0_0_z_zz_y_xx = buffer_2000_pdpd[1608];

    auto g_yz_0_0_0_z_zz_y_xy = buffer_2000_pdpd[1609];

    auto g_yz_0_0_0_z_zz_y_xz = buffer_2000_pdpd[1610];

    auto g_yz_0_0_0_z_zz_y_yy = buffer_2000_pdpd[1611];

    auto g_yz_0_0_0_z_zz_y_yz = buffer_2000_pdpd[1612];

    auto g_yz_0_0_0_z_zz_y_zz = buffer_2000_pdpd[1613];

    auto g_yz_0_0_0_z_zz_z_xx = buffer_2000_pdpd[1614];

    auto g_yz_0_0_0_z_zz_z_xy = buffer_2000_pdpd[1615];

    auto g_yz_0_0_0_z_zz_z_xz = buffer_2000_pdpd[1616];

    auto g_yz_0_0_0_z_zz_z_yy = buffer_2000_pdpd[1617];

    auto g_yz_0_0_0_z_zz_z_yz = buffer_2000_pdpd[1618];

    auto g_yz_0_0_0_z_zz_z_zz = buffer_2000_pdpd[1619];

    auto g_zz_0_0_0_x_xx_x_xx = buffer_2000_pdpd[1620];

    auto g_zz_0_0_0_x_xx_x_xy = buffer_2000_pdpd[1621];

    auto g_zz_0_0_0_x_xx_x_xz = buffer_2000_pdpd[1622];

    auto g_zz_0_0_0_x_xx_x_yy = buffer_2000_pdpd[1623];

    auto g_zz_0_0_0_x_xx_x_yz = buffer_2000_pdpd[1624];

    auto g_zz_0_0_0_x_xx_x_zz = buffer_2000_pdpd[1625];

    auto g_zz_0_0_0_x_xx_y_xx = buffer_2000_pdpd[1626];

    auto g_zz_0_0_0_x_xx_y_xy = buffer_2000_pdpd[1627];

    auto g_zz_0_0_0_x_xx_y_xz = buffer_2000_pdpd[1628];

    auto g_zz_0_0_0_x_xx_y_yy = buffer_2000_pdpd[1629];

    auto g_zz_0_0_0_x_xx_y_yz = buffer_2000_pdpd[1630];

    auto g_zz_0_0_0_x_xx_y_zz = buffer_2000_pdpd[1631];

    auto g_zz_0_0_0_x_xx_z_xx = buffer_2000_pdpd[1632];

    auto g_zz_0_0_0_x_xx_z_xy = buffer_2000_pdpd[1633];

    auto g_zz_0_0_0_x_xx_z_xz = buffer_2000_pdpd[1634];

    auto g_zz_0_0_0_x_xx_z_yy = buffer_2000_pdpd[1635];

    auto g_zz_0_0_0_x_xx_z_yz = buffer_2000_pdpd[1636];

    auto g_zz_0_0_0_x_xx_z_zz = buffer_2000_pdpd[1637];

    auto g_zz_0_0_0_x_xy_x_xx = buffer_2000_pdpd[1638];

    auto g_zz_0_0_0_x_xy_x_xy = buffer_2000_pdpd[1639];

    auto g_zz_0_0_0_x_xy_x_xz = buffer_2000_pdpd[1640];

    auto g_zz_0_0_0_x_xy_x_yy = buffer_2000_pdpd[1641];

    auto g_zz_0_0_0_x_xy_x_yz = buffer_2000_pdpd[1642];

    auto g_zz_0_0_0_x_xy_x_zz = buffer_2000_pdpd[1643];

    auto g_zz_0_0_0_x_xy_y_xx = buffer_2000_pdpd[1644];

    auto g_zz_0_0_0_x_xy_y_xy = buffer_2000_pdpd[1645];

    auto g_zz_0_0_0_x_xy_y_xz = buffer_2000_pdpd[1646];

    auto g_zz_0_0_0_x_xy_y_yy = buffer_2000_pdpd[1647];

    auto g_zz_0_0_0_x_xy_y_yz = buffer_2000_pdpd[1648];

    auto g_zz_0_0_0_x_xy_y_zz = buffer_2000_pdpd[1649];

    auto g_zz_0_0_0_x_xy_z_xx = buffer_2000_pdpd[1650];

    auto g_zz_0_0_0_x_xy_z_xy = buffer_2000_pdpd[1651];

    auto g_zz_0_0_0_x_xy_z_xz = buffer_2000_pdpd[1652];

    auto g_zz_0_0_0_x_xy_z_yy = buffer_2000_pdpd[1653];

    auto g_zz_0_0_0_x_xy_z_yz = buffer_2000_pdpd[1654];

    auto g_zz_0_0_0_x_xy_z_zz = buffer_2000_pdpd[1655];

    auto g_zz_0_0_0_x_xz_x_xx = buffer_2000_pdpd[1656];

    auto g_zz_0_0_0_x_xz_x_xy = buffer_2000_pdpd[1657];

    auto g_zz_0_0_0_x_xz_x_xz = buffer_2000_pdpd[1658];

    auto g_zz_0_0_0_x_xz_x_yy = buffer_2000_pdpd[1659];

    auto g_zz_0_0_0_x_xz_x_yz = buffer_2000_pdpd[1660];

    auto g_zz_0_0_0_x_xz_x_zz = buffer_2000_pdpd[1661];

    auto g_zz_0_0_0_x_xz_y_xx = buffer_2000_pdpd[1662];

    auto g_zz_0_0_0_x_xz_y_xy = buffer_2000_pdpd[1663];

    auto g_zz_0_0_0_x_xz_y_xz = buffer_2000_pdpd[1664];

    auto g_zz_0_0_0_x_xz_y_yy = buffer_2000_pdpd[1665];

    auto g_zz_0_0_0_x_xz_y_yz = buffer_2000_pdpd[1666];

    auto g_zz_0_0_0_x_xz_y_zz = buffer_2000_pdpd[1667];

    auto g_zz_0_0_0_x_xz_z_xx = buffer_2000_pdpd[1668];

    auto g_zz_0_0_0_x_xz_z_xy = buffer_2000_pdpd[1669];

    auto g_zz_0_0_0_x_xz_z_xz = buffer_2000_pdpd[1670];

    auto g_zz_0_0_0_x_xz_z_yy = buffer_2000_pdpd[1671];

    auto g_zz_0_0_0_x_xz_z_yz = buffer_2000_pdpd[1672];

    auto g_zz_0_0_0_x_xz_z_zz = buffer_2000_pdpd[1673];

    auto g_zz_0_0_0_x_yy_x_xx = buffer_2000_pdpd[1674];

    auto g_zz_0_0_0_x_yy_x_xy = buffer_2000_pdpd[1675];

    auto g_zz_0_0_0_x_yy_x_xz = buffer_2000_pdpd[1676];

    auto g_zz_0_0_0_x_yy_x_yy = buffer_2000_pdpd[1677];

    auto g_zz_0_0_0_x_yy_x_yz = buffer_2000_pdpd[1678];

    auto g_zz_0_0_0_x_yy_x_zz = buffer_2000_pdpd[1679];

    auto g_zz_0_0_0_x_yy_y_xx = buffer_2000_pdpd[1680];

    auto g_zz_0_0_0_x_yy_y_xy = buffer_2000_pdpd[1681];

    auto g_zz_0_0_0_x_yy_y_xz = buffer_2000_pdpd[1682];

    auto g_zz_0_0_0_x_yy_y_yy = buffer_2000_pdpd[1683];

    auto g_zz_0_0_0_x_yy_y_yz = buffer_2000_pdpd[1684];

    auto g_zz_0_0_0_x_yy_y_zz = buffer_2000_pdpd[1685];

    auto g_zz_0_0_0_x_yy_z_xx = buffer_2000_pdpd[1686];

    auto g_zz_0_0_0_x_yy_z_xy = buffer_2000_pdpd[1687];

    auto g_zz_0_0_0_x_yy_z_xz = buffer_2000_pdpd[1688];

    auto g_zz_0_0_0_x_yy_z_yy = buffer_2000_pdpd[1689];

    auto g_zz_0_0_0_x_yy_z_yz = buffer_2000_pdpd[1690];

    auto g_zz_0_0_0_x_yy_z_zz = buffer_2000_pdpd[1691];

    auto g_zz_0_0_0_x_yz_x_xx = buffer_2000_pdpd[1692];

    auto g_zz_0_0_0_x_yz_x_xy = buffer_2000_pdpd[1693];

    auto g_zz_0_0_0_x_yz_x_xz = buffer_2000_pdpd[1694];

    auto g_zz_0_0_0_x_yz_x_yy = buffer_2000_pdpd[1695];

    auto g_zz_0_0_0_x_yz_x_yz = buffer_2000_pdpd[1696];

    auto g_zz_0_0_0_x_yz_x_zz = buffer_2000_pdpd[1697];

    auto g_zz_0_0_0_x_yz_y_xx = buffer_2000_pdpd[1698];

    auto g_zz_0_0_0_x_yz_y_xy = buffer_2000_pdpd[1699];

    auto g_zz_0_0_0_x_yz_y_xz = buffer_2000_pdpd[1700];

    auto g_zz_0_0_0_x_yz_y_yy = buffer_2000_pdpd[1701];

    auto g_zz_0_0_0_x_yz_y_yz = buffer_2000_pdpd[1702];

    auto g_zz_0_0_0_x_yz_y_zz = buffer_2000_pdpd[1703];

    auto g_zz_0_0_0_x_yz_z_xx = buffer_2000_pdpd[1704];

    auto g_zz_0_0_0_x_yz_z_xy = buffer_2000_pdpd[1705];

    auto g_zz_0_0_0_x_yz_z_xz = buffer_2000_pdpd[1706];

    auto g_zz_0_0_0_x_yz_z_yy = buffer_2000_pdpd[1707];

    auto g_zz_0_0_0_x_yz_z_yz = buffer_2000_pdpd[1708];

    auto g_zz_0_0_0_x_yz_z_zz = buffer_2000_pdpd[1709];

    auto g_zz_0_0_0_x_zz_x_xx = buffer_2000_pdpd[1710];

    auto g_zz_0_0_0_x_zz_x_xy = buffer_2000_pdpd[1711];

    auto g_zz_0_0_0_x_zz_x_xz = buffer_2000_pdpd[1712];

    auto g_zz_0_0_0_x_zz_x_yy = buffer_2000_pdpd[1713];

    auto g_zz_0_0_0_x_zz_x_yz = buffer_2000_pdpd[1714];

    auto g_zz_0_0_0_x_zz_x_zz = buffer_2000_pdpd[1715];

    auto g_zz_0_0_0_x_zz_y_xx = buffer_2000_pdpd[1716];

    auto g_zz_0_0_0_x_zz_y_xy = buffer_2000_pdpd[1717];

    auto g_zz_0_0_0_x_zz_y_xz = buffer_2000_pdpd[1718];

    auto g_zz_0_0_0_x_zz_y_yy = buffer_2000_pdpd[1719];

    auto g_zz_0_0_0_x_zz_y_yz = buffer_2000_pdpd[1720];

    auto g_zz_0_0_0_x_zz_y_zz = buffer_2000_pdpd[1721];

    auto g_zz_0_0_0_x_zz_z_xx = buffer_2000_pdpd[1722];

    auto g_zz_0_0_0_x_zz_z_xy = buffer_2000_pdpd[1723];

    auto g_zz_0_0_0_x_zz_z_xz = buffer_2000_pdpd[1724];

    auto g_zz_0_0_0_x_zz_z_yy = buffer_2000_pdpd[1725];

    auto g_zz_0_0_0_x_zz_z_yz = buffer_2000_pdpd[1726];

    auto g_zz_0_0_0_x_zz_z_zz = buffer_2000_pdpd[1727];

    auto g_zz_0_0_0_y_xx_x_xx = buffer_2000_pdpd[1728];

    auto g_zz_0_0_0_y_xx_x_xy = buffer_2000_pdpd[1729];

    auto g_zz_0_0_0_y_xx_x_xz = buffer_2000_pdpd[1730];

    auto g_zz_0_0_0_y_xx_x_yy = buffer_2000_pdpd[1731];

    auto g_zz_0_0_0_y_xx_x_yz = buffer_2000_pdpd[1732];

    auto g_zz_0_0_0_y_xx_x_zz = buffer_2000_pdpd[1733];

    auto g_zz_0_0_0_y_xx_y_xx = buffer_2000_pdpd[1734];

    auto g_zz_0_0_0_y_xx_y_xy = buffer_2000_pdpd[1735];

    auto g_zz_0_0_0_y_xx_y_xz = buffer_2000_pdpd[1736];

    auto g_zz_0_0_0_y_xx_y_yy = buffer_2000_pdpd[1737];

    auto g_zz_0_0_0_y_xx_y_yz = buffer_2000_pdpd[1738];

    auto g_zz_0_0_0_y_xx_y_zz = buffer_2000_pdpd[1739];

    auto g_zz_0_0_0_y_xx_z_xx = buffer_2000_pdpd[1740];

    auto g_zz_0_0_0_y_xx_z_xy = buffer_2000_pdpd[1741];

    auto g_zz_0_0_0_y_xx_z_xz = buffer_2000_pdpd[1742];

    auto g_zz_0_0_0_y_xx_z_yy = buffer_2000_pdpd[1743];

    auto g_zz_0_0_0_y_xx_z_yz = buffer_2000_pdpd[1744];

    auto g_zz_0_0_0_y_xx_z_zz = buffer_2000_pdpd[1745];

    auto g_zz_0_0_0_y_xy_x_xx = buffer_2000_pdpd[1746];

    auto g_zz_0_0_0_y_xy_x_xy = buffer_2000_pdpd[1747];

    auto g_zz_0_0_0_y_xy_x_xz = buffer_2000_pdpd[1748];

    auto g_zz_0_0_0_y_xy_x_yy = buffer_2000_pdpd[1749];

    auto g_zz_0_0_0_y_xy_x_yz = buffer_2000_pdpd[1750];

    auto g_zz_0_0_0_y_xy_x_zz = buffer_2000_pdpd[1751];

    auto g_zz_0_0_0_y_xy_y_xx = buffer_2000_pdpd[1752];

    auto g_zz_0_0_0_y_xy_y_xy = buffer_2000_pdpd[1753];

    auto g_zz_0_0_0_y_xy_y_xz = buffer_2000_pdpd[1754];

    auto g_zz_0_0_0_y_xy_y_yy = buffer_2000_pdpd[1755];

    auto g_zz_0_0_0_y_xy_y_yz = buffer_2000_pdpd[1756];

    auto g_zz_0_0_0_y_xy_y_zz = buffer_2000_pdpd[1757];

    auto g_zz_0_0_0_y_xy_z_xx = buffer_2000_pdpd[1758];

    auto g_zz_0_0_0_y_xy_z_xy = buffer_2000_pdpd[1759];

    auto g_zz_0_0_0_y_xy_z_xz = buffer_2000_pdpd[1760];

    auto g_zz_0_0_0_y_xy_z_yy = buffer_2000_pdpd[1761];

    auto g_zz_0_0_0_y_xy_z_yz = buffer_2000_pdpd[1762];

    auto g_zz_0_0_0_y_xy_z_zz = buffer_2000_pdpd[1763];

    auto g_zz_0_0_0_y_xz_x_xx = buffer_2000_pdpd[1764];

    auto g_zz_0_0_0_y_xz_x_xy = buffer_2000_pdpd[1765];

    auto g_zz_0_0_0_y_xz_x_xz = buffer_2000_pdpd[1766];

    auto g_zz_0_0_0_y_xz_x_yy = buffer_2000_pdpd[1767];

    auto g_zz_0_0_0_y_xz_x_yz = buffer_2000_pdpd[1768];

    auto g_zz_0_0_0_y_xz_x_zz = buffer_2000_pdpd[1769];

    auto g_zz_0_0_0_y_xz_y_xx = buffer_2000_pdpd[1770];

    auto g_zz_0_0_0_y_xz_y_xy = buffer_2000_pdpd[1771];

    auto g_zz_0_0_0_y_xz_y_xz = buffer_2000_pdpd[1772];

    auto g_zz_0_0_0_y_xz_y_yy = buffer_2000_pdpd[1773];

    auto g_zz_0_0_0_y_xz_y_yz = buffer_2000_pdpd[1774];

    auto g_zz_0_0_0_y_xz_y_zz = buffer_2000_pdpd[1775];

    auto g_zz_0_0_0_y_xz_z_xx = buffer_2000_pdpd[1776];

    auto g_zz_0_0_0_y_xz_z_xy = buffer_2000_pdpd[1777];

    auto g_zz_0_0_0_y_xz_z_xz = buffer_2000_pdpd[1778];

    auto g_zz_0_0_0_y_xz_z_yy = buffer_2000_pdpd[1779];

    auto g_zz_0_0_0_y_xz_z_yz = buffer_2000_pdpd[1780];

    auto g_zz_0_0_0_y_xz_z_zz = buffer_2000_pdpd[1781];

    auto g_zz_0_0_0_y_yy_x_xx = buffer_2000_pdpd[1782];

    auto g_zz_0_0_0_y_yy_x_xy = buffer_2000_pdpd[1783];

    auto g_zz_0_0_0_y_yy_x_xz = buffer_2000_pdpd[1784];

    auto g_zz_0_0_0_y_yy_x_yy = buffer_2000_pdpd[1785];

    auto g_zz_0_0_0_y_yy_x_yz = buffer_2000_pdpd[1786];

    auto g_zz_0_0_0_y_yy_x_zz = buffer_2000_pdpd[1787];

    auto g_zz_0_0_0_y_yy_y_xx = buffer_2000_pdpd[1788];

    auto g_zz_0_0_0_y_yy_y_xy = buffer_2000_pdpd[1789];

    auto g_zz_0_0_0_y_yy_y_xz = buffer_2000_pdpd[1790];

    auto g_zz_0_0_0_y_yy_y_yy = buffer_2000_pdpd[1791];

    auto g_zz_0_0_0_y_yy_y_yz = buffer_2000_pdpd[1792];

    auto g_zz_0_0_0_y_yy_y_zz = buffer_2000_pdpd[1793];

    auto g_zz_0_0_0_y_yy_z_xx = buffer_2000_pdpd[1794];

    auto g_zz_0_0_0_y_yy_z_xy = buffer_2000_pdpd[1795];

    auto g_zz_0_0_0_y_yy_z_xz = buffer_2000_pdpd[1796];

    auto g_zz_0_0_0_y_yy_z_yy = buffer_2000_pdpd[1797];

    auto g_zz_0_0_0_y_yy_z_yz = buffer_2000_pdpd[1798];

    auto g_zz_0_0_0_y_yy_z_zz = buffer_2000_pdpd[1799];

    auto g_zz_0_0_0_y_yz_x_xx = buffer_2000_pdpd[1800];

    auto g_zz_0_0_0_y_yz_x_xy = buffer_2000_pdpd[1801];

    auto g_zz_0_0_0_y_yz_x_xz = buffer_2000_pdpd[1802];

    auto g_zz_0_0_0_y_yz_x_yy = buffer_2000_pdpd[1803];

    auto g_zz_0_0_0_y_yz_x_yz = buffer_2000_pdpd[1804];

    auto g_zz_0_0_0_y_yz_x_zz = buffer_2000_pdpd[1805];

    auto g_zz_0_0_0_y_yz_y_xx = buffer_2000_pdpd[1806];

    auto g_zz_0_0_0_y_yz_y_xy = buffer_2000_pdpd[1807];

    auto g_zz_0_0_0_y_yz_y_xz = buffer_2000_pdpd[1808];

    auto g_zz_0_0_0_y_yz_y_yy = buffer_2000_pdpd[1809];

    auto g_zz_0_0_0_y_yz_y_yz = buffer_2000_pdpd[1810];

    auto g_zz_0_0_0_y_yz_y_zz = buffer_2000_pdpd[1811];

    auto g_zz_0_0_0_y_yz_z_xx = buffer_2000_pdpd[1812];

    auto g_zz_0_0_0_y_yz_z_xy = buffer_2000_pdpd[1813];

    auto g_zz_0_0_0_y_yz_z_xz = buffer_2000_pdpd[1814];

    auto g_zz_0_0_0_y_yz_z_yy = buffer_2000_pdpd[1815];

    auto g_zz_0_0_0_y_yz_z_yz = buffer_2000_pdpd[1816];

    auto g_zz_0_0_0_y_yz_z_zz = buffer_2000_pdpd[1817];

    auto g_zz_0_0_0_y_zz_x_xx = buffer_2000_pdpd[1818];

    auto g_zz_0_0_0_y_zz_x_xy = buffer_2000_pdpd[1819];

    auto g_zz_0_0_0_y_zz_x_xz = buffer_2000_pdpd[1820];

    auto g_zz_0_0_0_y_zz_x_yy = buffer_2000_pdpd[1821];

    auto g_zz_0_0_0_y_zz_x_yz = buffer_2000_pdpd[1822];

    auto g_zz_0_0_0_y_zz_x_zz = buffer_2000_pdpd[1823];

    auto g_zz_0_0_0_y_zz_y_xx = buffer_2000_pdpd[1824];

    auto g_zz_0_0_0_y_zz_y_xy = buffer_2000_pdpd[1825];

    auto g_zz_0_0_0_y_zz_y_xz = buffer_2000_pdpd[1826];

    auto g_zz_0_0_0_y_zz_y_yy = buffer_2000_pdpd[1827];

    auto g_zz_0_0_0_y_zz_y_yz = buffer_2000_pdpd[1828];

    auto g_zz_0_0_0_y_zz_y_zz = buffer_2000_pdpd[1829];

    auto g_zz_0_0_0_y_zz_z_xx = buffer_2000_pdpd[1830];

    auto g_zz_0_0_0_y_zz_z_xy = buffer_2000_pdpd[1831];

    auto g_zz_0_0_0_y_zz_z_xz = buffer_2000_pdpd[1832];

    auto g_zz_0_0_0_y_zz_z_yy = buffer_2000_pdpd[1833];

    auto g_zz_0_0_0_y_zz_z_yz = buffer_2000_pdpd[1834];

    auto g_zz_0_0_0_y_zz_z_zz = buffer_2000_pdpd[1835];

    auto g_zz_0_0_0_z_xx_x_xx = buffer_2000_pdpd[1836];

    auto g_zz_0_0_0_z_xx_x_xy = buffer_2000_pdpd[1837];

    auto g_zz_0_0_0_z_xx_x_xz = buffer_2000_pdpd[1838];

    auto g_zz_0_0_0_z_xx_x_yy = buffer_2000_pdpd[1839];

    auto g_zz_0_0_0_z_xx_x_yz = buffer_2000_pdpd[1840];

    auto g_zz_0_0_0_z_xx_x_zz = buffer_2000_pdpd[1841];

    auto g_zz_0_0_0_z_xx_y_xx = buffer_2000_pdpd[1842];

    auto g_zz_0_0_0_z_xx_y_xy = buffer_2000_pdpd[1843];

    auto g_zz_0_0_0_z_xx_y_xz = buffer_2000_pdpd[1844];

    auto g_zz_0_0_0_z_xx_y_yy = buffer_2000_pdpd[1845];

    auto g_zz_0_0_0_z_xx_y_yz = buffer_2000_pdpd[1846];

    auto g_zz_0_0_0_z_xx_y_zz = buffer_2000_pdpd[1847];

    auto g_zz_0_0_0_z_xx_z_xx = buffer_2000_pdpd[1848];

    auto g_zz_0_0_0_z_xx_z_xy = buffer_2000_pdpd[1849];

    auto g_zz_0_0_0_z_xx_z_xz = buffer_2000_pdpd[1850];

    auto g_zz_0_0_0_z_xx_z_yy = buffer_2000_pdpd[1851];

    auto g_zz_0_0_0_z_xx_z_yz = buffer_2000_pdpd[1852];

    auto g_zz_0_0_0_z_xx_z_zz = buffer_2000_pdpd[1853];

    auto g_zz_0_0_0_z_xy_x_xx = buffer_2000_pdpd[1854];

    auto g_zz_0_0_0_z_xy_x_xy = buffer_2000_pdpd[1855];

    auto g_zz_0_0_0_z_xy_x_xz = buffer_2000_pdpd[1856];

    auto g_zz_0_0_0_z_xy_x_yy = buffer_2000_pdpd[1857];

    auto g_zz_0_0_0_z_xy_x_yz = buffer_2000_pdpd[1858];

    auto g_zz_0_0_0_z_xy_x_zz = buffer_2000_pdpd[1859];

    auto g_zz_0_0_0_z_xy_y_xx = buffer_2000_pdpd[1860];

    auto g_zz_0_0_0_z_xy_y_xy = buffer_2000_pdpd[1861];

    auto g_zz_0_0_0_z_xy_y_xz = buffer_2000_pdpd[1862];

    auto g_zz_0_0_0_z_xy_y_yy = buffer_2000_pdpd[1863];

    auto g_zz_0_0_0_z_xy_y_yz = buffer_2000_pdpd[1864];

    auto g_zz_0_0_0_z_xy_y_zz = buffer_2000_pdpd[1865];

    auto g_zz_0_0_0_z_xy_z_xx = buffer_2000_pdpd[1866];

    auto g_zz_0_0_0_z_xy_z_xy = buffer_2000_pdpd[1867];

    auto g_zz_0_0_0_z_xy_z_xz = buffer_2000_pdpd[1868];

    auto g_zz_0_0_0_z_xy_z_yy = buffer_2000_pdpd[1869];

    auto g_zz_0_0_0_z_xy_z_yz = buffer_2000_pdpd[1870];

    auto g_zz_0_0_0_z_xy_z_zz = buffer_2000_pdpd[1871];

    auto g_zz_0_0_0_z_xz_x_xx = buffer_2000_pdpd[1872];

    auto g_zz_0_0_0_z_xz_x_xy = buffer_2000_pdpd[1873];

    auto g_zz_0_0_0_z_xz_x_xz = buffer_2000_pdpd[1874];

    auto g_zz_0_0_0_z_xz_x_yy = buffer_2000_pdpd[1875];

    auto g_zz_0_0_0_z_xz_x_yz = buffer_2000_pdpd[1876];

    auto g_zz_0_0_0_z_xz_x_zz = buffer_2000_pdpd[1877];

    auto g_zz_0_0_0_z_xz_y_xx = buffer_2000_pdpd[1878];

    auto g_zz_0_0_0_z_xz_y_xy = buffer_2000_pdpd[1879];

    auto g_zz_0_0_0_z_xz_y_xz = buffer_2000_pdpd[1880];

    auto g_zz_0_0_0_z_xz_y_yy = buffer_2000_pdpd[1881];

    auto g_zz_0_0_0_z_xz_y_yz = buffer_2000_pdpd[1882];

    auto g_zz_0_0_0_z_xz_y_zz = buffer_2000_pdpd[1883];

    auto g_zz_0_0_0_z_xz_z_xx = buffer_2000_pdpd[1884];

    auto g_zz_0_0_0_z_xz_z_xy = buffer_2000_pdpd[1885];

    auto g_zz_0_0_0_z_xz_z_xz = buffer_2000_pdpd[1886];

    auto g_zz_0_0_0_z_xz_z_yy = buffer_2000_pdpd[1887];

    auto g_zz_0_0_0_z_xz_z_yz = buffer_2000_pdpd[1888];

    auto g_zz_0_0_0_z_xz_z_zz = buffer_2000_pdpd[1889];

    auto g_zz_0_0_0_z_yy_x_xx = buffer_2000_pdpd[1890];

    auto g_zz_0_0_0_z_yy_x_xy = buffer_2000_pdpd[1891];

    auto g_zz_0_0_0_z_yy_x_xz = buffer_2000_pdpd[1892];

    auto g_zz_0_0_0_z_yy_x_yy = buffer_2000_pdpd[1893];

    auto g_zz_0_0_0_z_yy_x_yz = buffer_2000_pdpd[1894];

    auto g_zz_0_0_0_z_yy_x_zz = buffer_2000_pdpd[1895];

    auto g_zz_0_0_0_z_yy_y_xx = buffer_2000_pdpd[1896];

    auto g_zz_0_0_0_z_yy_y_xy = buffer_2000_pdpd[1897];

    auto g_zz_0_0_0_z_yy_y_xz = buffer_2000_pdpd[1898];

    auto g_zz_0_0_0_z_yy_y_yy = buffer_2000_pdpd[1899];

    auto g_zz_0_0_0_z_yy_y_yz = buffer_2000_pdpd[1900];

    auto g_zz_0_0_0_z_yy_y_zz = buffer_2000_pdpd[1901];

    auto g_zz_0_0_0_z_yy_z_xx = buffer_2000_pdpd[1902];

    auto g_zz_0_0_0_z_yy_z_xy = buffer_2000_pdpd[1903];

    auto g_zz_0_0_0_z_yy_z_xz = buffer_2000_pdpd[1904];

    auto g_zz_0_0_0_z_yy_z_yy = buffer_2000_pdpd[1905];

    auto g_zz_0_0_0_z_yy_z_yz = buffer_2000_pdpd[1906];

    auto g_zz_0_0_0_z_yy_z_zz = buffer_2000_pdpd[1907];

    auto g_zz_0_0_0_z_yz_x_xx = buffer_2000_pdpd[1908];

    auto g_zz_0_0_0_z_yz_x_xy = buffer_2000_pdpd[1909];

    auto g_zz_0_0_0_z_yz_x_xz = buffer_2000_pdpd[1910];

    auto g_zz_0_0_0_z_yz_x_yy = buffer_2000_pdpd[1911];

    auto g_zz_0_0_0_z_yz_x_yz = buffer_2000_pdpd[1912];

    auto g_zz_0_0_0_z_yz_x_zz = buffer_2000_pdpd[1913];

    auto g_zz_0_0_0_z_yz_y_xx = buffer_2000_pdpd[1914];

    auto g_zz_0_0_0_z_yz_y_xy = buffer_2000_pdpd[1915];

    auto g_zz_0_0_0_z_yz_y_xz = buffer_2000_pdpd[1916];

    auto g_zz_0_0_0_z_yz_y_yy = buffer_2000_pdpd[1917];

    auto g_zz_0_0_0_z_yz_y_yz = buffer_2000_pdpd[1918];

    auto g_zz_0_0_0_z_yz_y_zz = buffer_2000_pdpd[1919];

    auto g_zz_0_0_0_z_yz_z_xx = buffer_2000_pdpd[1920];

    auto g_zz_0_0_0_z_yz_z_xy = buffer_2000_pdpd[1921];

    auto g_zz_0_0_0_z_yz_z_xz = buffer_2000_pdpd[1922];

    auto g_zz_0_0_0_z_yz_z_yy = buffer_2000_pdpd[1923];

    auto g_zz_0_0_0_z_yz_z_yz = buffer_2000_pdpd[1924];

    auto g_zz_0_0_0_z_yz_z_zz = buffer_2000_pdpd[1925];

    auto g_zz_0_0_0_z_zz_x_xx = buffer_2000_pdpd[1926];

    auto g_zz_0_0_0_z_zz_x_xy = buffer_2000_pdpd[1927];

    auto g_zz_0_0_0_z_zz_x_xz = buffer_2000_pdpd[1928];

    auto g_zz_0_0_0_z_zz_x_yy = buffer_2000_pdpd[1929];

    auto g_zz_0_0_0_z_zz_x_yz = buffer_2000_pdpd[1930];

    auto g_zz_0_0_0_z_zz_x_zz = buffer_2000_pdpd[1931];

    auto g_zz_0_0_0_z_zz_y_xx = buffer_2000_pdpd[1932];

    auto g_zz_0_0_0_z_zz_y_xy = buffer_2000_pdpd[1933];

    auto g_zz_0_0_0_z_zz_y_xz = buffer_2000_pdpd[1934];

    auto g_zz_0_0_0_z_zz_y_yy = buffer_2000_pdpd[1935];

    auto g_zz_0_0_0_z_zz_y_yz = buffer_2000_pdpd[1936];

    auto g_zz_0_0_0_z_zz_y_zz = buffer_2000_pdpd[1937];

    auto g_zz_0_0_0_z_zz_z_xx = buffer_2000_pdpd[1938];

    auto g_zz_0_0_0_z_zz_z_xy = buffer_2000_pdpd[1939];

    auto g_zz_0_0_0_z_zz_z_xz = buffer_2000_pdpd[1940];

    auto g_zz_0_0_0_z_zz_z_yy = buffer_2000_pdpd[1941];

    auto g_zz_0_0_0_z_zz_z_yz = buffer_2000_pdpd[1942];

    auto g_zz_0_0_0_z_zz_z_zz = buffer_2000_pdpd[1943];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_xx_x_xx, g_x_xx_x_xy, g_x_xx_x_xz, g_x_xx_x_yy, g_x_xx_x_yz, g_x_xx_x_zz, g_xx_0_0_0_x_xx_x_xx, g_xx_0_0_0_x_xx_x_xy, g_xx_0_0_0_x_xx_x_xz, g_xx_0_0_0_x_xx_x_yy, g_xx_0_0_0_x_xx_x_yz, g_xx_0_0_0_x_xx_x_zz, g_xxx_xx_x_xx, g_xxx_xx_x_xy, g_xxx_xx_x_xz, g_xxx_xx_x_yy, g_xxx_xx_x_yz, g_xxx_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_xx_x_xx[i] = -6.0 * g_x_xx_x_xx[i] * a_exp + 4.0 * g_xxx_xx_x_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xx_x_xy[i] = -6.0 * g_x_xx_x_xy[i] * a_exp + 4.0 * g_xxx_xx_x_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xx_x_xz[i] = -6.0 * g_x_xx_x_xz[i] * a_exp + 4.0 * g_xxx_xx_x_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xx_x_yy[i] = -6.0 * g_x_xx_x_yy[i] * a_exp + 4.0 * g_xxx_xx_x_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xx_x_yz[i] = -6.0 * g_x_xx_x_yz[i] * a_exp + 4.0 * g_xxx_xx_x_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xx_x_zz[i] = -6.0 * g_x_xx_x_zz[i] * a_exp + 4.0 * g_xxx_xx_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_xx_y_xx, g_x_xx_y_xy, g_x_xx_y_xz, g_x_xx_y_yy, g_x_xx_y_yz, g_x_xx_y_zz, g_xx_0_0_0_x_xx_y_xx, g_xx_0_0_0_x_xx_y_xy, g_xx_0_0_0_x_xx_y_xz, g_xx_0_0_0_x_xx_y_yy, g_xx_0_0_0_x_xx_y_yz, g_xx_0_0_0_x_xx_y_zz, g_xxx_xx_y_xx, g_xxx_xx_y_xy, g_xxx_xx_y_xz, g_xxx_xx_y_yy, g_xxx_xx_y_yz, g_xxx_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_xx_y_xx[i] = -6.0 * g_x_xx_y_xx[i] * a_exp + 4.0 * g_xxx_xx_y_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xx_y_xy[i] = -6.0 * g_x_xx_y_xy[i] * a_exp + 4.0 * g_xxx_xx_y_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xx_y_xz[i] = -6.0 * g_x_xx_y_xz[i] * a_exp + 4.0 * g_xxx_xx_y_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xx_y_yy[i] = -6.0 * g_x_xx_y_yy[i] * a_exp + 4.0 * g_xxx_xx_y_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xx_y_yz[i] = -6.0 * g_x_xx_y_yz[i] * a_exp + 4.0 * g_xxx_xx_y_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xx_y_zz[i] = -6.0 * g_x_xx_y_zz[i] * a_exp + 4.0 * g_xxx_xx_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_xx_z_xx, g_x_xx_z_xy, g_x_xx_z_xz, g_x_xx_z_yy, g_x_xx_z_yz, g_x_xx_z_zz, g_xx_0_0_0_x_xx_z_xx, g_xx_0_0_0_x_xx_z_xy, g_xx_0_0_0_x_xx_z_xz, g_xx_0_0_0_x_xx_z_yy, g_xx_0_0_0_x_xx_z_yz, g_xx_0_0_0_x_xx_z_zz, g_xxx_xx_z_xx, g_xxx_xx_z_xy, g_xxx_xx_z_xz, g_xxx_xx_z_yy, g_xxx_xx_z_yz, g_xxx_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_xx_z_xx[i] = -6.0 * g_x_xx_z_xx[i] * a_exp + 4.0 * g_xxx_xx_z_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xx_z_xy[i] = -6.0 * g_x_xx_z_xy[i] * a_exp + 4.0 * g_xxx_xx_z_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xx_z_xz[i] = -6.0 * g_x_xx_z_xz[i] * a_exp + 4.0 * g_xxx_xx_z_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xx_z_yy[i] = -6.0 * g_x_xx_z_yy[i] * a_exp + 4.0 * g_xxx_xx_z_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xx_z_yz[i] = -6.0 * g_x_xx_z_yz[i] * a_exp + 4.0 * g_xxx_xx_z_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xx_z_zz[i] = -6.0 * g_x_xx_z_zz[i] * a_exp + 4.0 * g_xxx_xx_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_x_xy_x_xx, g_x_xy_x_xy, g_x_xy_x_xz, g_x_xy_x_yy, g_x_xy_x_yz, g_x_xy_x_zz, g_xx_0_0_0_x_xy_x_xx, g_xx_0_0_0_x_xy_x_xy, g_xx_0_0_0_x_xy_x_xz, g_xx_0_0_0_x_xy_x_yy, g_xx_0_0_0_x_xy_x_yz, g_xx_0_0_0_x_xy_x_zz, g_xxx_xy_x_xx, g_xxx_xy_x_xy, g_xxx_xy_x_xz, g_xxx_xy_x_yy, g_xxx_xy_x_yz, g_xxx_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_xy_x_xx[i] = -6.0 * g_x_xy_x_xx[i] * a_exp + 4.0 * g_xxx_xy_x_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xy_x_xy[i] = -6.0 * g_x_xy_x_xy[i] * a_exp + 4.0 * g_xxx_xy_x_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xy_x_xz[i] = -6.0 * g_x_xy_x_xz[i] * a_exp + 4.0 * g_xxx_xy_x_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xy_x_yy[i] = -6.0 * g_x_xy_x_yy[i] * a_exp + 4.0 * g_xxx_xy_x_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xy_x_yz[i] = -6.0 * g_x_xy_x_yz[i] * a_exp + 4.0 * g_xxx_xy_x_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xy_x_zz[i] = -6.0 * g_x_xy_x_zz[i] * a_exp + 4.0 * g_xxx_xy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_x_xy_y_xx, g_x_xy_y_xy, g_x_xy_y_xz, g_x_xy_y_yy, g_x_xy_y_yz, g_x_xy_y_zz, g_xx_0_0_0_x_xy_y_xx, g_xx_0_0_0_x_xy_y_xy, g_xx_0_0_0_x_xy_y_xz, g_xx_0_0_0_x_xy_y_yy, g_xx_0_0_0_x_xy_y_yz, g_xx_0_0_0_x_xy_y_zz, g_xxx_xy_y_xx, g_xxx_xy_y_xy, g_xxx_xy_y_xz, g_xxx_xy_y_yy, g_xxx_xy_y_yz, g_xxx_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_xy_y_xx[i] = -6.0 * g_x_xy_y_xx[i] * a_exp + 4.0 * g_xxx_xy_y_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xy_y_xy[i] = -6.0 * g_x_xy_y_xy[i] * a_exp + 4.0 * g_xxx_xy_y_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xy_y_xz[i] = -6.0 * g_x_xy_y_xz[i] * a_exp + 4.0 * g_xxx_xy_y_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xy_y_yy[i] = -6.0 * g_x_xy_y_yy[i] * a_exp + 4.0 * g_xxx_xy_y_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xy_y_yz[i] = -6.0 * g_x_xy_y_yz[i] * a_exp + 4.0 * g_xxx_xy_y_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xy_y_zz[i] = -6.0 * g_x_xy_y_zz[i] * a_exp + 4.0 * g_xxx_xy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_x_xy_z_xx, g_x_xy_z_xy, g_x_xy_z_xz, g_x_xy_z_yy, g_x_xy_z_yz, g_x_xy_z_zz, g_xx_0_0_0_x_xy_z_xx, g_xx_0_0_0_x_xy_z_xy, g_xx_0_0_0_x_xy_z_xz, g_xx_0_0_0_x_xy_z_yy, g_xx_0_0_0_x_xy_z_yz, g_xx_0_0_0_x_xy_z_zz, g_xxx_xy_z_xx, g_xxx_xy_z_xy, g_xxx_xy_z_xz, g_xxx_xy_z_yy, g_xxx_xy_z_yz, g_xxx_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_xy_z_xx[i] = -6.0 * g_x_xy_z_xx[i] * a_exp + 4.0 * g_xxx_xy_z_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xy_z_xy[i] = -6.0 * g_x_xy_z_xy[i] * a_exp + 4.0 * g_xxx_xy_z_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xy_z_xz[i] = -6.0 * g_x_xy_z_xz[i] * a_exp + 4.0 * g_xxx_xy_z_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xy_z_yy[i] = -6.0 * g_x_xy_z_yy[i] * a_exp + 4.0 * g_xxx_xy_z_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xy_z_yz[i] = -6.0 * g_x_xy_z_yz[i] * a_exp + 4.0 * g_xxx_xy_z_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xy_z_zz[i] = -6.0 * g_x_xy_z_zz[i] * a_exp + 4.0 * g_xxx_xy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_x_xz_x_xx, g_x_xz_x_xy, g_x_xz_x_xz, g_x_xz_x_yy, g_x_xz_x_yz, g_x_xz_x_zz, g_xx_0_0_0_x_xz_x_xx, g_xx_0_0_0_x_xz_x_xy, g_xx_0_0_0_x_xz_x_xz, g_xx_0_0_0_x_xz_x_yy, g_xx_0_0_0_x_xz_x_yz, g_xx_0_0_0_x_xz_x_zz, g_xxx_xz_x_xx, g_xxx_xz_x_xy, g_xxx_xz_x_xz, g_xxx_xz_x_yy, g_xxx_xz_x_yz, g_xxx_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_xz_x_xx[i] = -6.0 * g_x_xz_x_xx[i] * a_exp + 4.0 * g_xxx_xz_x_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xz_x_xy[i] = -6.0 * g_x_xz_x_xy[i] * a_exp + 4.0 * g_xxx_xz_x_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xz_x_xz[i] = -6.0 * g_x_xz_x_xz[i] * a_exp + 4.0 * g_xxx_xz_x_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xz_x_yy[i] = -6.0 * g_x_xz_x_yy[i] * a_exp + 4.0 * g_xxx_xz_x_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xz_x_yz[i] = -6.0 * g_x_xz_x_yz[i] * a_exp + 4.0 * g_xxx_xz_x_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xz_x_zz[i] = -6.0 * g_x_xz_x_zz[i] * a_exp + 4.0 * g_xxx_xz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_x_xz_y_xx, g_x_xz_y_xy, g_x_xz_y_xz, g_x_xz_y_yy, g_x_xz_y_yz, g_x_xz_y_zz, g_xx_0_0_0_x_xz_y_xx, g_xx_0_0_0_x_xz_y_xy, g_xx_0_0_0_x_xz_y_xz, g_xx_0_0_0_x_xz_y_yy, g_xx_0_0_0_x_xz_y_yz, g_xx_0_0_0_x_xz_y_zz, g_xxx_xz_y_xx, g_xxx_xz_y_xy, g_xxx_xz_y_xz, g_xxx_xz_y_yy, g_xxx_xz_y_yz, g_xxx_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_xz_y_xx[i] = -6.0 * g_x_xz_y_xx[i] * a_exp + 4.0 * g_xxx_xz_y_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xz_y_xy[i] = -6.0 * g_x_xz_y_xy[i] * a_exp + 4.0 * g_xxx_xz_y_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xz_y_xz[i] = -6.0 * g_x_xz_y_xz[i] * a_exp + 4.0 * g_xxx_xz_y_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xz_y_yy[i] = -6.0 * g_x_xz_y_yy[i] * a_exp + 4.0 * g_xxx_xz_y_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xz_y_yz[i] = -6.0 * g_x_xz_y_yz[i] * a_exp + 4.0 * g_xxx_xz_y_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xz_y_zz[i] = -6.0 * g_x_xz_y_zz[i] * a_exp + 4.0 * g_xxx_xz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_x_xz_z_xx, g_x_xz_z_xy, g_x_xz_z_xz, g_x_xz_z_yy, g_x_xz_z_yz, g_x_xz_z_zz, g_xx_0_0_0_x_xz_z_xx, g_xx_0_0_0_x_xz_z_xy, g_xx_0_0_0_x_xz_z_xz, g_xx_0_0_0_x_xz_z_yy, g_xx_0_0_0_x_xz_z_yz, g_xx_0_0_0_x_xz_z_zz, g_xxx_xz_z_xx, g_xxx_xz_z_xy, g_xxx_xz_z_xz, g_xxx_xz_z_yy, g_xxx_xz_z_yz, g_xxx_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_xz_z_xx[i] = -6.0 * g_x_xz_z_xx[i] * a_exp + 4.0 * g_xxx_xz_z_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xz_z_xy[i] = -6.0 * g_x_xz_z_xy[i] * a_exp + 4.0 * g_xxx_xz_z_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xz_z_xz[i] = -6.0 * g_x_xz_z_xz[i] * a_exp + 4.0 * g_xxx_xz_z_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xz_z_yy[i] = -6.0 * g_x_xz_z_yy[i] * a_exp + 4.0 * g_xxx_xz_z_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xz_z_yz[i] = -6.0 * g_x_xz_z_yz[i] * a_exp + 4.0 * g_xxx_xz_z_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xz_z_zz[i] = -6.0 * g_x_xz_z_zz[i] * a_exp + 4.0 * g_xxx_xz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_x_yy_x_xx, g_x_yy_x_xy, g_x_yy_x_xz, g_x_yy_x_yy, g_x_yy_x_yz, g_x_yy_x_zz, g_xx_0_0_0_x_yy_x_xx, g_xx_0_0_0_x_yy_x_xy, g_xx_0_0_0_x_yy_x_xz, g_xx_0_0_0_x_yy_x_yy, g_xx_0_0_0_x_yy_x_yz, g_xx_0_0_0_x_yy_x_zz, g_xxx_yy_x_xx, g_xxx_yy_x_xy, g_xxx_yy_x_xz, g_xxx_yy_x_yy, g_xxx_yy_x_yz, g_xxx_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_yy_x_xx[i] = -6.0 * g_x_yy_x_xx[i] * a_exp + 4.0 * g_xxx_yy_x_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yy_x_xy[i] = -6.0 * g_x_yy_x_xy[i] * a_exp + 4.0 * g_xxx_yy_x_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yy_x_xz[i] = -6.0 * g_x_yy_x_xz[i] * a_exp + 4.0 * g_xxx_yy_x_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yy_x_yy[i] = -6.0 * g_x_yy_x_yy[i] * a_exp + 4.0 * g_xxx_yy_x_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yy_x_yz[i] = -6.0 * g_x_yy_x_yz[i] * a_exp + 4.0 * g_xxx_yy_x_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yy_x_zz[i] = -6.0 * g_x_yy_x_zz[i] * a_exp + 4.0 * g_xxx_yy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_x_yy_y_xx, g_x_yy_y_xy, g_x_yy_y_xz, g_x_yy_y_yy, g_x_yy_y_yz, g_x_yy_y_zz, g_xx_0_0_0_x_yy_y_xx, g_xx_0_0_0_x_yy_y_xy, g_xx_0_0_0_x_yy_y_xz, g_xx_0_0_0_x_yy_y_yy, g_xx_0_0_0_x_yy_y_yz, g_xx_0_0_0_x_yy_y_zz, g_xxx_yy_y_xx, g_xxx_yy_y_xy, g_xxx_yy_y_xz, g_xxx_yy_y_yy, g_xxx_yy_y_yz, g_xxx_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_yy_y_xx[i] = -6.0 * g_x_yy_y_xx[i] * a_exp + 4.0 * g_xxx_yy_y_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yy_y_xy[i] = -6.0 * g_x_yy_y_xy[i] * a_exp + 4.0 * g_xxx_yy_y_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yy_y_xz[i] = -6.0 * g_x_yy_y_xz[i] * a_exp + 4.0 * g_xxx_yy_y_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yy_y_yy[i] = -6.0 * g_x_yy_y_yy[i] * a_exp + 4.0 * g_xxx_yy_y_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yy_y_yz[i] = -6.0 * g_x_yy_y_yz[i] * a_exp + 4.0 * g_xxx_yy_y_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yy_y_zz[i] = -6.0 * g_x_yy_y_zz[i] * a_exp + 4.0 * g_xxx_yy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_x_yy_z_xx, g_x_yy_z_xy, g_x_yy_z_xz, g_x_yy_z_yy, g_x_yy_z_yz, g_x_yy_z_zz, g_xx_0_0_0_x_yy_z_xx, g_xx_0_0_0_x_yy_z_xy, g_xx_0_0_0_x_yy_z_xz, g_xx_0_0_0_x_yy_z_yy, g_xx_0_0_0_x_yy_z_yz, g_xx_0_0_0_x_yy_z_zz, g_xxx_yy_z_xx, g_xxx_yy_z_xy, g_xxx_yy_z_xz, g_xxx_yy_z_yy, g_xxx_yy_z_yz, g_xxx_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_yy_z_xx[i] = -6.0 * g_x_yy_z_xx[i] * a_exp + 4.0 * g_xxx_yy_z_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yy_z_xy[i] = -6.0 * g_x_yy_z_xy[i] * a_exp + 4.0 * g_xxx_yy_z_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yy_z_xz[i] = -6.0 * g_x_yy_z_xz[i] * a_exp + 4.0 * g_xxx_yy_z_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yy_z_yy[i] = -6.0 * g_x_yy_z_yy[i] * a_exp + 4.0 * g_xxx_yy_z_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yy_z_yz[i] = -6.0 * g_x_yy_z_yz[i] * a_exp + 4.0 * g_xxx_yy_z_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yy_z_zz[i] = -6.0 * g_x_yy_z_zz[i] * a_exp + 4.0 * g_xxx_yy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_x_yz_x_xx, g_x_yz_x_xy, g_x_yz_x_xz, g_x_yz_x_yy, g_x_yz_x_yz, g_x_yz_x_zz, g_xx_0_0_0_x_yz_x_xx, g_xx_0_0_0_x_yz_x_xy, g_xx_0_0_0_x_yz_x_xz, g_xx_0_0_0_x_yz_x_yy, g_xx_0_0_0_x_yz_x_yz, g_xx_0_0_0_x_yz_x_zz, g_xxx_yz_x_xx, g_xxx_yz_x_xy, g_xxx_yz_x_xz, g_xxx_yz_x_yy, g_xxx_yz_x_yz, g_xxx_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_yz_x_xx[i] = -6.0 * g_x_yz_x_xx[i] * a_exp + 4.0 * g_xxx_yz_x_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yz_x_xy[i] = -6.0 * g_x_yz_x_xy[i] * a_exp + 4.0 * g_xxx_yz_x_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yz_x_xz[i] = -6.0 * g_x_yz_x_xz[i] * a_exp + 4.0 * g_xxx_yz_x_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yz_x_yy[i] = -6.0 * g_x_yz_x_yy[i] * a_exp + 4.0 * g_xxx_yz_x_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yz_x_yz[i] = -6.0 * g_x_yz_x_yz[i] * a_exp + 4.0 * g_xxx_yz_x_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yz_x_zz[i] = -6.0 * g_x_yz_x_zz[i] * a_exp + 4.0 * g_xxx_yz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_x_yz_y_xx, g_x_yz_y_xy, g_x_yz_y_xz, g_x_yz_y_yy, g_x_yz_y_yz, g_x_yz_y_zz, g_xx_0_0_0_x_yz_y_xx, g_xx_0_0_0_x_yz_y_xy, g_xx_0_0_0_x_yz_y_xz, g_xx_0_0_0_x_yz_y_yy, g_xx_0_0_0_x_yz_y_yz, g_xx_0_0_0_x_yz_y_zz, g_xxx_yz_y_xx, g_xxx_yz_y_xy, g_xxx_yz_y_xz, g_xxx_yz_y_yy, g_xxx_yz_y_yz, g_xxx_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_yz_y_xx[i] = -6.0 * g_x_yz_y_xx[i] * a_exp + 4.0 * g_xxx_yz_y_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yz_y_xy[i] = -6.0 * g_x_yz_y_xy[i] * a_exp + 4.0 * g_xxx_yz_y_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yz_y_xz[i] = -6.0 * g_x_yz_y_xz[i] * a_exp + 4.0 * g_xxx_yz_y_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yz_y_yy[i] = -6.0 * g_x_yz_y_yy[i] * a_exp + 4.0 * g_xxx_yz_y_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yz_y_yz[i] = -6.0 * g_x_yz_y_yz[i] * a_exp + 4.0 * g_xxx_yz_y_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yz_y_zz[i] = -6.0 * g_x_yz_y_zz[i] * a_exp + 4.0 * g_xxx_yz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_x_yz_z_xx, g_x_yz_z_xy, g_x_yz_z_xz, g_x_yz_z_yy, g_x_yz_z_yz, g_x_yz_z_zz, g_xx_0_0_0_x_yz_z_xx, g_xx_0_0_0_x_yz_z_xy, g_xx_0_0_0_x_yz_z_xz, g_xx_0_0_0_x_yz_z_yy, g_xx_0_0_0_x_yz_z_yz, g_xx_0_0_0_x_yz_z_zz, g_xxx_yz_z_xx, g_xxx_yz_z_xy, g_xxx_yz_z_xz, g_xxx_yz_z_yy, g_xxx_yz_z_yz, g_xxx_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_yz_z_xx[i] = -6.0 * g_x_yz_z_xx[i] * a_exp + 4.0 * g_xxx_yz_z_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yz_z_xy[i] = -6.0 * g_x_yz_z_xy[i] * a_exp + 4.0 * g_xxx_yz_z_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yz_z_xz[i] = -6.0 * g_x_yz_z_xz[i] * a_exp + 4.0 * g_xxx_yz_z_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yz_z_yy[i] = -6.0 * g_x_yz_z_yy[i] * a_exp + 4.0 * g_xxx_yz_z_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yz_z_yz[i] = -6.0 * g_x_yz_z_yz[i] * a_exp + 4.0 * g_xxx_yz_z_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yz_z_zz[i] = -6.0 * g_x_yz_z_zz[i] * a_exp + 4.0 * g_xxx_yz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_x_zz_x_xx, g_x_zz_x_xy, g_x_zz_x_xz, g_x_zz_x_yy, g_x_zz_x_yz, g_x_zz_x_zz, g_xx_0_0_0_x_zz_x_xx, g_xx_0_0_0_x_zz_x_xy, g_xx_0_0_0_x_zz_x_xz, g_xx_0_0_0_x_zz_x_yy, g_xx_0_0_0_x_zz_x_yz, g_xx_0_0_0_x_zz_x_zz, g_xxx_zz_x_xx, g_xxx_zz_x_xy, g_xxx_zz_x_xz, g_xxx_zz_x_yy, g_xxx_zz_x_yz, g_xxx_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_zz_x_xx[i] = -6.0 * g_x_zz_x_xx[i] * a_exp + 4.0 * g_xxx_zz_x_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_zz_x_xy[i] = -6.0 * g_x_zz_x_xy[i] * a_exp + 4.0 * g_xxx_zz_x_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_zz_x_xz[i] = -6.0 * g_x_zz_x_xz[i] * a_exp + 4.0 * g_xxx_zz_x_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_zz_x_yy[i] = -6.0 * g_x_zz_x_yy[i] * a_exp + 4.0 * g_xxx_zz_x_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_zz_x_yz[i] = -6.0 * g_x_zz_x_yz[i] * a_exp + 4.0 * g_xxx_zz_x_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_zz_x_zz[i] = -6.0 * g_x_zz_x_zz[i] * a_exp + 4.0 * g_xxx_zz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_x_zz_y_xx, g_x_zz_y_xy, g_x_zz_y_xz, g_x_zz_y_yy, g_x_zz_y_yz, g_x_zz_y_zz, g_xx_0_0_0_x_zz_y_xx, g_xx_0_0_0_x_zz_y_xy, g_xx_0_0_0_x_zz_y_xz, g_xx_0_0_0_x_zz_y_yy, g_xx_0_0_0_x_zz_y_yz, g_xx_0_0_0_x_zz_y_zz, g_xxx_zz_y_xx, g_xxx_zz_y_xy, g_xxx_zz_y_xz, g_xxx_zz_y_yy, g_xxx_zz_y_yz, g_xxx_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_zz_y_xx[i] = -6.0 * g_x_zz_y_xx[i] * a_exp + 4.0 * g_xxx_zz_y_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_zz_y_xy[i] = -6.0 * g_x_zz_y_xy[i] * a_exp + 4.0 * g_xxx_zz_y_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_zz_y_xz[i] = -6.0 * g_x_zz_y_xz[i] * a_exp + 4.0 * g_xxx_zz_y_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_zz_y_yy[i] = -6.0 * g_x_zz_y_yy[i] * a_exp + 4.0 * g_xxx_zz_y_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_zz_y_yz[i] = -6.0 * g_x_zz_y_yz[i] * a_exp + 4.0 * g_xxx_zz_y_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_zz_y_zz[i] = -6.0 * g_x_zz_y_zz[i] * a_exp + 4.0 * g_xxx_zz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_x_zz_z_xx, g_x_zz_z_xy, g_x_zz_z_xz, g_x_zz_z_yy, g_x_zz_z_yz, g_x_zz_z_zz, g_xx_0_0_0_x_zz_z_xx, g_xx_0_0_0_x_zz_z_xy, g_xx_0_0_0_x_zz_z_xz, g_xx_0_0_0_x_zz_z_yy, g_xx_0_0_0_x_zz_z_yz, g_xx_0_0_0_x_zz_z_zz, g_xxx_zz_z_xx, g_xxx_zz_z_xy, g_xxx_zz_z_xz, g_xxx_zz_z_yy, g_xxx_zz_z_yz, g_xxx_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_zz_z_xx[i] = -6.0 * g_x_zz_z_xx[i] * a_exp + 4.0 * g_xxx_zz_z_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_x_zz_z_xy[i] = -6.0 * g_x_zz_z_xy[i] * a_exp + 4.0 * g_xxx_zz_z_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_zz_z_xz[i] = -6.0 * g_x_zz_z_xz[i] * a_exp + 4.0 * g_xxx_zz_z_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_zz_z_yy[i] = -6.0 * g_x_zz_z_yy[i] * a_exp + 4.0 * g_xxx_zz_z_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_x_zz_z_yz[i] = -6.0 * g_x_zz_z_yz[i] * a_exp + 4.0 * g_xxx_zz_z_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_x_zz_z_zz[i] = -6.0 * g_x_zz_z_zz[i] * a_exp + 4.0 * g_xxx_zz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_xx_0_0_0_y_xx_x_xx, g_xx_0_0_0_y_xx_x_xy, g_xx_0_0_0_y_xx_x_xz, g_xx_0_0_0_y_xx_x_yy, g_xx_0_0_0_y_xx_x_yz, g_xx_0_0_0_y_xx_x_zz, g_xxy_xx_x_xx, g_xxy_xx_x_xy, g_xxy_xx_x_xz, g_xxy_xx_x_yy, g_xxy_xx_x_yz, g_xxy_xx_x_zz, g_y_xx_x_xx, g_y_xx_x_xy, g_y_xx_x_xz, g_y_xx_x_yy, g_y_xx_x_yz, g_y_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_xx_x_xx[i] = -2.0 * g_y_xx_x_xx[i] * a_exp + 4.0 * g_xxy_xx_x_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xx_x_xy[i] = -2.0 * g_y_xx_x_xy[i] * a_exp + 4.0 * g_xxy_xx_x_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xx_x_xz[i] = -2.0 * g_y_xx_x_xz[i] * a_exp + 4.0 * g_xxy_xx_x_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xx_x_yy[i] = -2.0 * g_y_xx_x_yy[i] * a_exp + 4.0 * g_xxy_xx_x_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xx_x_yz[i] = -2.0 * g_y_xx_x_yz[i] * a_exp + 4.0 * g_xxy_xx_x_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xx_x_zz[i] = -2.0 * g_y_xx_x_zz[i] * a_exp + 4.0 * g_xxy_xx_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_xx_0_0_0_y_xx_y_xx, g_xx_0_0_0_y_xx_y_xy, g_xx_0_0_0_y_xx_y_xz, g_xx_0_0_0_y_xx_y_yy, g_xx_0_0_0_y_xx_y_yz, g_xx_0_0_0_y_xx_y_zz, g_xxy_xx_y_xx, g_xxy_xx_y_xy, g_xxy_xx_y_xz, g_xxy_xx_y_yy, g_xxy_xx_y_yz, g_xxy_xx_y_zz, g_y_xx_y_xx, g_y_xx_y_xy, g_y_xx_y_xz, g_y_xx_y_yy, g_y_xx_y_yz, g_y_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_xx_y_xx[i] = -2.0 * g_y_xx_y_xx[i] * a_exp + 4.0 * g_xxy_xx_y_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xx_y_xy[i] = -2.0 * g_y_xx_y_xy[i] * a_exp + 4.0 * g_xxy_xx_y_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xx_y_xz[i] = -2.0 * g_y_xx_y_xz[i] * a_exp + 4.0 * g_xxy_xx_y_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xx_y_yy[i] = -2.0 * g_y_xx_y_yy[i] * a_exp + 4.0 * g_xxy_xx_y_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xx_y_yz[i] = -2.0 * g_y_xx_y_yz[i] * a_exp + 4.0 * g_xxy_xx_y_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xx_y_zz[i] = -2.0 * g_y_xx_y_zz[i] * a_exp + 4.0 * g_xxy_xx_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_xx_0_0_0_y_xx_z_xx, g_xx_0_0_0_y_xx_z_xy, g_xx_0_0_0_y_xx_z_xz, g_xx_0_0_0_y_xx_z_yy, g_xx_0_0_0_y_xx_z_yz, g_xx_0_0_0_y_xx_z_zz, g_xxy_xx_z_xx, g_xxy_xx_z_xy, g_xxy_xx_z_xz, g_xxy_xx_z_yy, g_xxy_xx_z_yz, g_xxy_xx_z_zz, g_y_xx_z_xx, g_y_xx_z_xy, g_y_xx_z_xz, g_y_xx_z_yy, g_y_xx_z_yz, g_y_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_xx_z_xx[i] = -2.0 * g_y_xx_z_xx[i] * a_exp + 4.0 * g_xxy_xx_z_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xx_z_xy[i] = -2.0 * g_y_xx_z_xy[i] * a_exp + 4.0 * g_xxy_xx_z_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xx_z_xz[i] = -2.0 * g_y_xx_z_xz[i] * a_exp + 4.0 * g_xxy_xx_z_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xx_z_yy[i] = -2.0 * g_y_xx_z_yy[i] * a_exp + 4.0 * g_xxy_xx_z_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xx_z_yz[i] = -2.0 * g_y_xx_z_yz[i] * a_exp + 4.0 * g_xxy_xx_z_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xx_z_zz[i] = -2.0 * g_y_xx_z_zz[i] * a_exp + 4.0 * g_xxy_xx_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_xx_0_0_0_y_xy_x_xx, g_xx_0_0_0_y_xy_x_xy, g_xx_0_0_0_y_xy_x_xz, g_xx_0_0_0_y_xy_x_yy, g_xx_0_0_0_y_xy_x_yz, g_xx_0_0_0_y_xy_x_zz, g_xxy_xy_x_xx, g_xxy_xy_x_xy, g_xxy_xy_x_xz, g_xxy_xy_x_yy, g_xxy_xy_x_yz, g_xxy_xy_x_zz, g_y_xy_x_xx, g_y_xy_x_xy, g_y_xy_x_xz, g_y_xy_x_yy, g_y_xy_x_yz, g_y_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_xy_x_xx[i] = -2.0 * g_y_xy_x_xx[i] * a_exp + 4.0 * g_xxy_xy_x_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xy_x_xy[i] = -2.0 * g_y_xy_x_xy[i] * a_exp + 4.0 * g_xxy_xy_x_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xy_x_xz[i] = -2.0 * g_y_xy_x_xz[i] * a_exp + 4.0 * g_xxy_xy_x_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xy_x_yy[i] = -2.0 * g_y_xy_x_yy[i] * a_exp + 4.0 * g_xxy_xy_x_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xy_x_yz[i] = -2.0 * g_y_xy_x_yz[i] * a_exp + 4.0 * g_xxy_xy_x_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xy_x_zz[i] = -2.0 * g_y_xy_x_zz[i] * a_exp + 4.0 * g_xxy_xy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_xx_0_0_0_y_xy_y_xx, g_xx_0_0_0_y_xy_y_xy, g_xx_0_0_0_y_xy_y_xz, g_xx_0_0_0_y_xy_y_yy, g_xx_0_0_0_y_xy_y_yz, g_xx_0_0_0_y_xy_y_zz, g_xxy_xy_y_xx, g_xxy_xy_y_xy, g_xxy_xy_y_xz, g_xxy_xy_y_yy, g_xxy_xy_y_yz, g_xxy_xy_y_zz, g_y_xy_y_xx, g_y_xy_y_xy, g_y_xy_y_xz, g_y_xy_y_yy, g_y_xy_y_yz, g_y_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_xy_y_xx[i] = -2.0 * g_y_xy_y_xx[i] * a_exp + 4.0 * g_xxy_xy_y_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xy_y_xy[i] = -2.0 * g_y_xy_y_xy[i] * a_exp + 4.0 * g_xxy_xy_y_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xy_y_xz[i] = -2.0 * g_y_xy_y_xz[i] * a_exp + 4.0 * g_xxy_xy_y_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xy_y_yy[i] = -2.0 * g_y_xy_y_yy[i] * a_exp + 4.0 * g_xxy_xy_y_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xy_y_yz[i] = -2.0 * g_y_xy_y_yz[i] * a_exp + 4.0 * g_xxy_xy_y_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xy_y_zz[i] = -2.0 * g_y_xy_y_zz[i] * a_exp + 4.0 * g_xxy_xy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_xx_0_0_0_y_xy_z_xx, g_xx_0_0_0_y_xy_z_xy, g_xx_0_0_0_y_xy_z_xz, g_xx_0_0_0_y_xy_z_yy, g_xx_0_0_0_y_xy_z_yz, g_xx_0_0_0_y_xy_z_zz, g_xxy_xy_z_xx, g_xxy_xy_z_xy, g_xxy_xy_z_xz, g_xxy_xy_z_yy, g_xxy_xy_z_yz, g_xxy_xy_z_zz, g_y_xy_z_xx, g_y_xy_z_xy, g_y_xy_z_xz, g_y_xy_z_yy, g_y_xy_z_yz, g_y_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_xy_z_xx[i] = -2.0 * g_y_xy_z_xx[i] * a_exp + 4.0 * g_xxy_xy_z_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xy_z_xy[i] = -2.0 * g_y_xy_z_xy[i] * a_exp + 4.0 * g_xxy_xy_z_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xy_z_xz[i] = -2.0 * g_y_xy_z_xz[i] * a_exp + 4.0 * g_xxy_xy_z_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xy_z_yy[i] = -2.0 * g_y_xy_z_yy[i] * a_exp + 4.0 * g_xxy_xy_z_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xy_z_yz[i] = -2.0 * g_y_xy_z_yz[i] * a_exp + 4.0 * g_xxy_xy_z_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xy_z_zz[i] = -2.0 * g_y_xy_z_zz[i] * a_exp + 4.0 * g_xxy_xy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_xx_0_0_0_y_xz_x_xx, g_xx_0_0_0_y_xz_x_xy, g_xx_0_0_0_y_xz_x_xz, g_xx_0_0_0_y_xz_x_yy, g_xx_0_0_0_y_xz_x_yz, g_xx_0_0_0_y_xz_x_zz, g_xxy_xz_x_xx, g_xxy_xz_x_xy, g_xxy_xz_x_xz, g_xxy_xz_x_yy, g_xxy_xz_x_yz, g_xxy_xz_x_zz, g_y_xz_x_xx, g_y_xz_x_xy, g_y_xz_x_xz, g_y_xz_x_yy, g_y_xz_x_yz, g_y_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_xz_x_xx[i] = -2.0 * g_y_xz_x_xx[i] * a_exp + 4.0 * g_xxy_xz_x_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xz_x_xy[i] = -2.0 * g_y_xz_x_xy[i] * a_exp + 4.0 * g_xxy_xz_x_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xz_x_xz[i] = -2.0 * g_y_xz_x_xz[i] * a_exp + 4.0 * g_xxy_xz_x_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xz_x_yy[i] = -2.0 * g_y_xz_x_yy[i] * a_exp + 4.0 * g_xxy_xz_x_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xz_x_yz[i] = -2.0 * g_y_xz_x_yz[i] * a_exp + 4.0 * g_xxy_xz_x_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xz_x_zz[i] = -2.0 * g_y_xz_x_zz[i] * a_exp + 4.0 * g_xxy_xz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_xx_0_0_0_y_xz_y_xx, g_xx_0_0_0_y_xz_y_xy, g_xx_0_0_0_y_xz_y_xz, g_xx_0_0_0_y_xz_y_yy, g_xx_0_0_0_y_xz_y_yz, g_xx_0_0_0_y_xz_y_zz, g_xxy_xz_y_xx, g_xxy_xz_y_xy, g_xxy_xz_y_xz, g_xxy_xz_y_yy, g_xxy_xz_y_yz, g_xxy_xz_y_zz, g_y_xz_y_xx, g_y_xz_y_xy, g_y_xz_y_xz, g_y_xz_y_yy, g_y_xz_y_yz, g_y_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_xz_y_xx[i] = -2.0 * g_y_xz_y_xx[i] * a_exp + 4.0 * g_xxy_xz_y_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xz_y_xy[i] = -2.0 * g_y_xz_y_xy[i] * a_exp + 4.0 * g_xxy_xz_y_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xz_y_xz[i] = -2.0 * g_y_xz_y_xz[i] * a_exp + 4.0 * g_xxy_xz_y_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xz_y_yy[i] = -2.0 * g_y_xz_y_yy[i] * a_exp + 4.0 * g_xxy_xz_y_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xz_y_yz[i] = -2.0 * g_y_xz_y_yz[i] * a_exp + 4.0 * g_xxy_xz_y_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xz_y_zz[i] = -2.0 * g_y_xz_y_zz[i] * a_exp + 4.0 * g_xxy_xz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_xx_0_0_0_y_xz_z_xx, g_xx_0_0_0_y_xz_z_xy, g_xx_0_0_0_y_xz_z_xz, g_xx_0_0_0_y_xz_z_yy, g_xx_0_0_0_y_xz_z_yz, g_xx_0_0_0_y_xz_z_zz, g_xxy_xz_z_xx, g_xxy_xz_z_xy, g_xxy_xz_z_xz, g_xxy_xz_z_yy, g_xxy_xz_z_yz, g_xxy_xz_z_zz, g_y_xz_z_xx, g_y_xz_z_xy, g_y_xz_z_xz, g_y_xz_z_yy, g_y_xz_z_yz, g_y_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_xz_z_xx[i] = -2.0 * g_y_xz_z_xx[i] * a_exp + 4.0 * g_xxy_xz_z_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xz_z_xy[i] = -2.0 * g_y_xz_z_xy[i] * a_exp + 4.0 * g_xxy_xz_z_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xz_z_xz[i] = -2.0 * g_y_xz_z_xz[i] * a_exp + 4.0 * g_xxy_xz_z_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xz_z_yy[i] = -2.0 * g_y_xz_z_yy[i] * a_exp + 4.0 * g_xxy_xz_z_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xz_z_yz[i] = -2.0 * g_y_xz_z_yz[i] * a_exp + 4.0 * g_xxy_xz_z_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xz_z_zz[i] = -2.0 * g_y_xz_z_zz[i] * a_exp + 4.0 * g_xxy_xz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_xx_0_0_0_y_yy_x_xx, g_xx_0_0_0_y_yy_x_xy, g_xx_0_0_0_y_yy_x_xz, g_xx_0_0_0_y_yy_x_yy, g_xx_0_0_0_y_yy_x_yz, g_xx_0_0_0_y_yy_x_zz, g_xxy_yy_x_xx, g_xxy_yy_x_xy, g_xxy_yy_x_xz, g_xxy_yy_x_yy, g_xxy_yy_x_yz, g_xxy_yy_x_zz, g_y_yy_x_xx, g_y_yy_x_xy, g_y_yy_x_xz, g_y_yy_x_yy, g_y_yy_x_yz, g_y_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_yy_x_xx[i] = -2.0 * g_y_yy_x_xx[i] * a_exp + 4.0 * g_xxy_yy_x_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yy_x_xy[i] = -2.0 * g_y_yy_x_xy[i] * a_exp + 4.0 * g_xxy_yy_x_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yy_x_xz[i] = -2.0 * g_y_yy_x_xz[i] * a_exp + 4.0 * g_xxy_yy_x_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yy_x_yy[i] = -2.0 * g_y_yy_x_yy[i] * a_exp + 4.0 * g_xxy_yy_x_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yy_x_yz[i] = -2.0 * g_y_yy_x_yz[i] * a_exp + 4.0 * g_xxy_yy_x_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yy_x_zz[i] = -2.0 * g_y_yy_x_zz[i] * a_exp + 4.0 * g_xxy_yy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_xx_0_0_0_y_yy_y_xx, g_xx_0_0_0_y_yy_y_xy, g_xx_0_0_0_y_yy_y_xz, g_xx_0_0_0_y_yy_y_yy, g_xx_0_0_0_y_yy_y_yz, g_xx_0_0_0_y_yy_y_zz, g_xxy_yy_y_xx, g_xxy_yy_y_xy, g_xxy_yy_y_xz, g_xxy_yy_y_yy, g_xxy_yy_y_yz, g_xxy_yy_y_zz, g_y_yy_y_xx, g_y_yy_y_xy, g_y_yy_y_xz, g_y_yy_y_yy, g_y_yy_y_yz, g_y_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_yy_y_xx[i] = -2.0 * g_y_yy_y_xx[i] * a_exp + 4.0 * g_xxy_yy_y_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yy_y_xy[i] = -2.0 * g_y_yy_y_xy[i] * a_exp + 4.0 * g_xxy_yy_y_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yy_y_xz[i] = -2.0 * g_y_yy_y_xz[i] * a_exp + 4.0 * g_xxy_yy_y_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yy_y_yy[i] = -2.0 * g_y_yy_y_yy[i] * a_exp + 4.0 * g_xxy_yy_y_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yy_y_yz[i] = -2.0 * g_y_yy_y_yz[i] * a_exp + 4.0 * g_xxy_yy_y_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yy_y_zz[i] = -2.0 * g_y_yy_y_zz[i] * a_exp + 4.0 * g_xxy_yy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_xx_0_0_0_y_yy_z_xx, g_xx_0_0_0_y_yy_z_xy, g_xx_0_0_0_y_yy_z_xz, g_xx_0_0_0_y_yy_z_yy, g_xx_0_0_0_y_yy_z_yz, g_xx_0_0_0_y_yy_z_zz, g_xxy_yy_z_xx, g_xxy_yy_z_xy, g_xxy_yy_z_xz, g_xxy_yy_z_yy, g_xxy_yy_z_yz, g_xxy_yy_z_zz, g_y_yy_z_xx, g_y_yy_z_xy, g_y_yy_z_xz, g_y_yy_z_yy, g_y_yy_z_yz, g_y_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_yy_z_xx[i] = -2.0 * g_y_yy_z_xx[i] * a_exp + 4.0 * g_xxy_yy_z_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yy_z_xy[i] = -2.0 * g_y_yy_z_xy[i] * a_exp + 4.0 * g_xxy_yy_z_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yy_z_xz[i] = -2.0 * g_y_yy_z_xz[i] * a_exp + 4.0 * g_xxy_yy_z_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yy_z_yy[i] = -2.0 * g_y_yy_z_yy[i] * a_exp + 4.0 * g_xxy_yy_z_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yy_z_yz[i] = -2.0 * g_y_yy_z_yz[i] * a_exp + 4.0 * g_xxy_yy_z_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yy_z_zz[i] = -2.0 * g_y_yy_z_zz[i] * a_exp + 4.0 * g_xxy_yy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_xx_0_0_0_y_yz_x_xx, g_xx_0_0_0_y_yz_x_xy, g_xx_0_0_0_y_yz_x_xz, g_xx_0_0_0_y_yz_x_yy, g_xx_0_0_0_y_yz_x_yz, g_xx_0_0_0_y_yz_x_zz, g_xxy_yz_x_xx, g_xxy_yz_x_xy, g_xxy_yz_x_xz, g_xxy_yz_x_yy, g_xxy_yz_x_yz, g_xxy_yz_x_zz, g_y_yz_x_xx, g_y_yz_x_xy, g_y_yz_x_xz, g_y_yz_x_yy, g_y_yz_x_yz, g_y_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_yz_x_xx[i] = -2.0 * g_y_yz_x_xx[i] * a_exp + 4.0 * g_xxy_yz_x_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yz_x_xy[i] = -2.0 * g_y_yz_x_xy[i] * a_exp + 4.0 * g_xxy_yz_x_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yz_x_xz[i] = -2.0 * g_y_yz_x_xz[i] * a_exp + 4.0 * g_xxy_yz_x_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yz_x_yy[i] = -2.0 * g_y_yz_x_yy[i] * a_exp + 4.0 * g_xxy_yz_x_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yz_x_yz[i] = -2.0 * g_y_yz_x_yz[i] * a_exp + 4.0 * g_xxy_yz_x_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yz_x_zz[i] = -2.0 * g_y_yz_x_zz[i] * a_exp + 4.0 * g_xxy_yz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_xx_0_0_0_y_yz_y_xx, g_xx_0_0_0_y_yz_y_xy, g_xx_0_0_0_y_yz_y_xz, g_xx_0_0_0_y_yz_y_yy, g_xx_0_0_0_y_yz_y_yz, g_xx_0_0_0_y_yz_y_zz, g_xxy_yz_y_xx, g_xxy_yz_y_xy, g_xxy_yz_y_xz, g_xxy_yz_y_yy, g_xxy_yz_y_yz, g_xxy_yz_y_zz, g_y_yz_y_xx, g_y_yz_y_xy, g_y_yz_y_xz, g_y_yz_y_yy, g_y_yz_y_yz, g_y_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_yz_y_xx[i] = -2.0 * g_y_yz_y_xx[i] * a_exp + 4.0 * g_xxy_yz_y_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yz_y_xy[i] = -2.0 * g_y_yz_y_xy[i] * a_exp + 4.0 * g_xxy_yz_y_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yz_y_xz[i] = -2.0 * g_y_yz_y_xz[i] * a_exp + 4.0 * g_xxy_yz_y_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yz_y_yy[i] = -2.0 * g_y_yz_y_yy[i] * a_exp + 4.0 * g_xxy_yz_y_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yz_y_yz[i] = -2.0 * g_y_yz_y_yz[i] * a_exp + 4.0 * g_xxy_yz_y_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yz_y_zz[i] = -2.0 * g_y_yz_y_zz[i] * a_exp + 4.0 * g_xxy_yz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_xx_0_0_0_y_yz_z_xx, g_xx_0_0_0_y_yz_z_xy, g_xx_0_0_0_y_yz_z_xz, g_xx_0_0_0_y_yz_z_yy, g_xx_0_0_0_y_yz_z_yz, g_xx_0_0_0_y_yz_z_zz, g_xxy_yz_z_xx, g_xxy_yz_z_xy, g_xxy_yz_z_xz, g_xxy_yz_z_yy, g_xxy_yz_z_yz, g_xxy_yz_z_zz, g_y_yz_z_xx, g_y_yz_z_xy, g_y_yz_z_xz, g_y_yz_z_yy, g_y_yz_z_yz, g_y_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_yz_z_xx[i] = -2.0 * g_y_yz_z_xx[i] * a_exp + 4.0 * g_xxy_yz_z_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yz_z_xy[i] = -2.0 * g_y_yz_z_xy[i] * a_exp + 4.0 * g_xxy_yz_z_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yz_z_xz[i] = -2.0 * g_y_yz_z_xz[i] * a_exp + 4.0 * g_xxy_yz_z_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yz_z_yy[i] = -2.0 * g_y_yz_z_yy[i] * a_exp + 4.0 * g_xxy_yz_z_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yz_z_yz[i] = -2.0 * g_y_yz_z_yz[i] * a_exp + 4.0 * g_xxy_yz_z_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yz_z_zz[i] = -2.0 * g_y_yz_z_zz[i] * a_exp + 4.0 * g_xxy_yz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_xx_0_0_0_y_zz_x_xx, g_xx_0_0_0_y_zz_x_xy, g_xx_0_0_0_y_zz_x_xz, g_xx_0_0_0_y_zz_x_yy, g_xx_0_0_0_y_zz_x_yz, g_xx_0_0_0_y_zz_x_zz, g_xxy_zz_x_xx, g_xxy_zz_x_xy, g_xxy_zz_x_xz, g_xxy_zz_x_yy, g_xxy_zz_x_yz, g_xxy_zz_x_zz, g_y_zz_x_xx, g_y_zz_x_xy, g_y_zz_x_xz, g_y_zz_x_yy, g_y_zz_x_yz, g_y_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_zz_x_xx[i] = -2.0 * g_y_zz_x_xx[i] * a_exp + 4.0 * g_xxy_zz_x_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_zz_x_xy[i] = -2.0 * g_y_zz_x_xy[i] * a_exp + 4.0 * g_xxy_zz_x_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_zz_x_xz[i] = -2.0 * g_y_zz_x_xz[i] * a_exp + 4.0 * g_xxy_zz_x_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_zz_x_yy[i] = -2.0 * g_y_zz_x_yy[i] * a_exp + 4.0 * g_xxy_zz_x_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_zz_x_yz[i] = -2.0 * g_y_zz_x_yz[i] * a_exp + 4.0 * g_xxy_zz_x_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_zz_x_zz[i] = -2.0 * g_y_zz_x_zz[i] * a_exp + 4.0 * g_xxy_zz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_xx_0_0_0_y_zz_y_xx, g_xx_0_0_0_y_zz_y_xy, g_xx_0_0_0_y_zz_y_xz, g_xx_0_0_0_y_zz_y_yy, g_xx_0_0_0_y_zz_y_yz, g_xx_0_0_0_y_zz_y_zz, g_xxy_zz_y_xx, g_xxy_zz_y_xy, g_xxy_zz_y_xz, g_xxy_zz_y_yy, g_xxy_zz_y_yz, g_xxy_zz_y_zz, g_y_zz_y_xx, g_y_zz_y_xy, g_y_zz_y_xz, g_y_zz_y_yy, g_y_zz_y_yz, g_y_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_zz_y_xx[i] = -2.0 * g_y_zz_y_xx[i] * a_exp + 4.0 * g_xxy_zz_y_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_zz_y_xy[i] = -2.0 * g_y_zz_y_xy[i] * a_exp + 4.0 * g_xxy_zz_y_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_zz_y_xz[i] = -2.0 * g_y_zz_y_xz[i] * a_exp + 4.0 * g_xxy_zz_y_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_zz_y_yy[i] = -2.0 * g_y_zz_y_yy[i] * a_exp + 4.0 * g_xxy_zz_y_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_zz_y_yz[i] = -2.0 * g_y_zz_y_yz[i] * a_exp + 4.0 * g_xxy_zz_y_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_zz_y_zz[i] = -2.0 * g_y_zz_y_zz[i] * a_exp + 4.0 * g_xxy_zz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_xx_0_0_0_y_zz_z_xx, g_xx_0_0_0_y_zz_z_xy, g_xx_0_0_0_y_zz_z_xz, g_xx_0_0_0_y_zz_z_yy, g_xx_0_0_0_y_zz_z_yz, g_xx_0_0_0_y_zz_z_zz, g_xxy_zz_z_xx, g_xxy_zz_z_xy, g_xxy_zz_z_xz, g_xxy_zz_z_yy, g_xxy_zz_z_yz, g_xxy_zz_z_zz, g_y_zz_z_xx, g_y_zz_z_xy, g_y_zz_z_xz, g_y_zz_z_yy, g_y_zz_z_yz, g_y_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_zz_z_xx[i] = -2.0 * g_y_zz_z_xx[i] * a_exp + 4.0 * g_xxy_zz_z_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_y_zz_z_xy[i] = -2.0 * g_y_zz_z_xy[i] * a_exp + 4.0 * g_xxy_zz_z_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_zz_z_xz[i] = -2.0 * g_y_zz_z_xz[i] * a_exp + 4.0 * g_xxy_zz_z_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_zz_z_yy[i] = -2.0 * g_y_zz_z_yy[i] * a_exp + 4.0 * g_xxy_zz_z_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_y_zz_z_yz[i] = -2.0 * g_y_zz_z_yz[i] * a_exp + 4.0 * g_xxy_zz_z_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_y_zz_z_zz[i] = -2.0 * g_y_zz_z_zz[i] * a_exp + 4.0 * g_xxy_zz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_xx_0_0_0_z_xx_x_xx, g_xx_0_0_0_z_xx_x_xy, g_xx_0_0_0_z_xx_x_xz, g_xx_0_0_0_z_xx_x_yy, g_xx_0_0_0_z_xx_x_yz, g_xx_0_0_0_z_xx_x_zz, g_xxz_xx_x_xx, g_xxz_xx_x_xy, g_xxz_xx_x_xz, g_xxz_xx_x_yy, g_xxz_xx_x_yz, g_xxz_xx_x_zz, g_z_xx_x_xx, g_z_xx_x_xy, g_z_xx_x_xz, g_z_xx_x_yy, g_z_xx_x_yz, g_z_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_xx_x_xx[i] = -2.0 * g_z_xx_x_xx[i] * a_exp + 4.0 * g_xxz_xx_x_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xx_x_xy[i] = -2.0 * g_z_xx_x_xy[i] * a_exp + 4.0 * g_xxz_xx_x_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xx_x_xz[i] = -2.0 * g_z_xx_x_xz[i] * a_exp + 4.0 * g_xxz_xx_x_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xx_x_yy[i] = -2.0 * g_z_xx_x_yy[i] * a_exp + 4.0 * g_xxz_xx_x_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xx_x_yz[i] = -2.0 * g_z_xx_x_yz[i] * a_exp + 4.0 * g_xxz_xx_x_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xx_x_zz[i] = -2.0 * g_z_xx_x_zz[i] * a_exp + 4.0 * g_xxz_xx_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_xx_0_0_0_z_xx_y_xx, g_xx_0_0_0_z_xx_y_xy, g_xx_0_0_0_z_xx_y_xz, g_xx_0_0_0_z_xx_y_yy, g_xx_0_0_0_z_xx_y_yz, g_xx_0_0_0_z_xx_y_zz, g_xxz_xx_y_xx, g_xxz_xx_y_xy, g_xxz_xx_y_xz, g_xxz_xx_y_yy, g_xxz_xx_y_yz, g_xxz_xx_y_zz, g_z_xx_y_xx, g_z_xx_y_xy, g_z_xx_y_xz, g_z_xx_y_yy, g_z_xx_y_yz, g_z_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_xx_y_xx[i] = -2.0 * g_z_xx_y_xx[i] * a_exp + 4.0 * g_xxz_xx_y_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xx_y_xy[i] = -2.0 * g_z_xx_y_xy[i] * a_exp + 4.0 * g_xxz_xx_y_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xx_y_xz[i] = -2.0 * g_z_xx_y_xz[i] * a_exp + 4.0 * g_xxz_xx_y_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xx_y_yy[i] = -2.0 * g_z_xx_y_yy[i] * a_exp + 4.0 * g_xxz_xx_y_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xx_y_yz[i] = -2.0 * g_z_xx_y_yz[i] * a_exp + 4.0 * g_xxz_xx_y_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xx_y_zz[i] = -2.0 * g_z_xx_y_zz[i] * a_exp + 4.0 * g_xxz_xx_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_xx_0_0_0_z_xx_z_xx, g_xx_0_0_0_z_xx_z_xy, g_xx_0_0_0_z_xx_z_xz, g_xx_0_0_0_z_xx_z_yy, g_xx_0_0_0_z_xx_z_yz, g_xx_0_0_0_z_xx_z_zz, g_xxz_xx_z_xx, g_xxz_xx_z_xy, g_xxz_xx_z_xz, g_xxz_xx_z_yy, g_xxz_xx_z_yz, g_xxz_xx_z_zz, g_z_xx_z_xx, g_z_xx_z_xy, g_z_xx_z_xz, g_z_xx_z_yy, g_z_xx_z_yz, g_z_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_xx_z_xx[i] = -2.0 * g_z_xx_z_xx[i] * a_exp + 4.0 * g_xxz_xx_z_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xx_z_xy[i] = -2.0 * g_z_xx_z_xy[i] * a_exp + 4.0 * g_xxz_xx_z_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xx_z_xz[i] = -2.0 * g_z_xx_z_xz[i] * a_exp + 4.0 * g_xxz_xx_z_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xx_z_yy[i] = -2.0 * g_z_xx_z_yy[i] * a_exp + 4.0 * g_xxz_xx_z_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xx_z_yz[i] = -2.0 * g_z_xx_z_yz[i] * a_exp + 4.0 * g_xxz_xx_z_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xx_z_zz[i] = -2.0 * g_z_xx_z_zz[i] * a_exp + 4.0 * g_xxz_xx_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_xx_0_0_0_z_xy_x_xx, g_xx_0_0_0_z_xy_x_xy, g_xx_0_0_0_z_xy_x_xz, g_xx_0_0_0_z_xy_x_yy, g_xx_0_0_0_z_xy_x_yz, g_xx_0_0_0_z_xy_x_zz, g_xxz_xy_x_xx, g_xxz_xy_x_xy, g_xxz_xy_x_xz, g_xxz_xy_x_yy, g_xxz_xy_x_yz, g_xxz_xy_x_zz, g_z_xy_x_xx, g_z_xy_x_xy, g_z_xy_x_xz, g_z_xy_x_yy, g_z_xy_x_yz, g_z_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_xy_x_xx[i] = -2.0 * g_z_xy_x_xx[i] * a_exp + 4.0 * g_xxz_xy_x_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xy_x_xy[i] = -2.0 * g_z_xy_x_xy[i] * a_exp + 4.0 * g_xxz_xy_x_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xy_x_xz[i] = -2.0 * g_z_xy_x_xz[i] * a_exp + 4.0 * g_xxz_xy_x_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xy_x_yy[i] = -2.0 * g_z_xy_x_yy[i] * a_exp + 4.0 * g_xxz_xy_x_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xy_x_yz[i] = -2.0 * g_z_xy_x_yz[i] * a_exp + 4.0 * g_xxz_xy_x_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xy_x_zz[i] = -2.0 * g_z_xy_x_zz[i] * a_exp + 4.0 * g_xxz_xy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_xx_0_0_0_z_xy_y_xx, g_xx_0_0_0_z_xy_y_xy, g_xx_0_0_0_z_xy_y_xz, g_xx_0_0_0_z_xy_y_yy, g_xx_0_0_0_z_xy_y_yz, g_xx_0_0_0_z_xy_y_zz, g_xxz_xy_y_xx, g_xxz_xy_y_xy, g_xxz_xy_y_xz, g_xxz_xy_y_yy, g_xxz_xy_y_yz, g_xxz_xy_y_zz, g_z_xy_y_xx, g_z_xy_y_xy, g_z_xy_y_xz, g_z_xy_y_yy, g_z_xy_y_yz, g_z_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_xy_y_xx[i] = -2.0 * g_z_xy_y_xx[i] * a_exp + 4.0 * g_xxz_xy_y_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xy_y_xy[i] = -2.0 * g_z_xy_y_xy[i] * a_exp + 4.0 * g_xxz_xy_y_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xy_y_xz[i] = -2.0 * g_z_xy_y_xz[i] * a_exp + 4.0 * g_xxz_xy_y_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xy_y_yy[i] = -2.0 * g_z_xy_y_yy[i] * a_exp + 4.0 * g_xxz_xy_y_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xy_y_yz[i] = -2.0 * g_z_xy_y_yz[i] * a_exp + 4.0 * g_xxz_xy_y_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xy_y_zz[i] = -2.0 * g_z_xy_y_zz[i] * a_exp + 4.0 * g_xxz_xy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_xx_0_0_0_z_xy_z_xx, g_xx_0_0_0_z_xy_z_xy, g_xx_0_0_0_z_xy_z_xz, g_xx_0_0_0_z_xy_z_yy, g_xx_0_0_0_z_xy_z_yz, g_xx_0_0_0_z_xy_z_zz, g_xxz_xy_z_xx, g_xxz_xy_z_xy, g_xxz_xy_z_xz, g_xxz_xy_z_yy, g_xxz_xy_z_yz, g_xxz_xy_z_zz, g_z_xy_z_xx, g_z_xy_z_xy, g_z_xy_z_xz, g_z_xy_z_yy, g_z_xy_z_yz, g_z_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_xy_z_xx[i] = -2.0 * g_z_xy_z_xx[i] * a_exp + 4.0 * g_xxz_xy_z_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xy_z_xy[i] = -2.0 * g_z_xy_z_xy[i] * a_exp + 4.0 * g_xxz_xy_z_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xy_z_xz[i] = -2.0 * g_z_xy_z_xz[i] * a_exp + 4.0 * g_xxz_xy_z_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xy_z_yy[i] = -2.0 * g_z_xy_z_yy[i] * a_exp + 4.0 * g_xxz_xy_z_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xy_z_yz[i] = -2.0 * g_z_xy_z_yz[i] * a_exp + 4.0 * g_xxz_xy_z_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xy_z_zz[i] = -2.0 * g_z_xy_z_zz[i] * a_exp + 4.0 * g_xxz_xy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_xx_0_0_0_z_xz_x_xx, g_xx_0_0_0_z_xz_x_xy, g_xx_0_0_0_z_xz_x_xz, g_xx_0_0_0_z_xz_x_yy, g_xx_0_0_0_z_xz_x_yz, g_xx_0_0_0_z_xz_x_zz, g_xxz_xz_x_xx, g_xxz_xz_x_xy, g_xxz_xz_x_xz, g_xxz_xz_x_yy, g_xxz_xz_x_yz, g_xxz_xz_x_zz, g_z_xz_x_xx, g_z_xz_x_xy, g_z_xz_x_xz, g_z_xz_x_yy, g_z_xz_x_yz, g_z_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_xz_x_xx[i] = -2.0 * g_z_xz_x_xx[i] * a_exp + 4.0 * g_xxz_xz_x_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xz_x_xy[i] = -2.0 * g_z_xz_x_xy[i] * a_exp + 4.0 * g_xxz_xz_x_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xz_x_xz[i] = -2.0 * g_z_xz_x_xz[i] * a_exp + 4.0 * g_xxz_xz_x_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xz_x_yy[i] = -2.0 * g_z_xz_x_yy[i] * a_exp + 4.0 * g_xxz_xz_x_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xz_x_yz[i] = -2.0 * g_z_xz_x_yz[i] * a_exp + 4.0 * g_xxz_xz_x_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xz_x_zz[i] = -2.0 * g_z_xz_x_zz[i] * a_exp + 4.0 * g_xxz_xz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_xx_0_0_0_z_xz_y_xx, g_xx_0_0_0_z_xz_y_xy, g_xx_0_0_0_z_xz_y_xz, g_xx_0_0_0_z_xz_y_yy, g_xx_0_0_0_z_xz_y_yz, g_xx_0_0_0_z_xz_y_zz, g_xxz_xz_y_xx, g_xxz_xz_y_xy, g_xxz_xz_y_xz, g_xxz_xz_y_yy, g_xxz_xz_y_yz, g_xxz_xz_y_zz, g_z_xz_y_xx, g_z_xz_y_xy, g_z_xz_y_xz, g_z_xz_y_yy, g_z_xz_y_yz, g_z_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_xz_y_xx[i] = -2.0 * g_z_xz_y_xx[i] * a_exp + 4.0 * g_xxz_xz_y_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xz_y_xy[i] = -2.0 * g_z_xz_y_xy[i] * a_exp + 4.0 * g_xxz_xz_y_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xz_y_xz[i] = -2.0 * g_z_xz_y_xz[i] * a_exp + 4.0 * g_xxz_xz_y_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xz_y_yy[i] = -2.0 * g_z_xz_y_yy[i] * a_exp + 4.0 * g_xxz_xz_y_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xz_y_yz[i] = -2.0 * g_z_xz_y_yz[i] * a_exp + 4.0 * g_xxz_xz_y_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xz_y_zz[i] = -2.0 * g_z_xz_y_zz[i] * a_exp + 4.0 * g_xxz_xz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_xx_0_0_0_z_xz_z_xx, g_xx_0_0_0_z_xz_z_xy, g_xx_0_0_0_z_xz_z_xz, g_xx_0_0_0_z_xz_z_yy, g_xx_0_0_0_z_xz_z_yz, g_xx_0_0_0_z_xz_z_zz, g_xxz_xz_z_xx, g_xxz_xz_z_xy, g_xxz_xz_z_xz, g_xxz_xz_z_yy, g_xxz_xz_z_yz, g_xxz_xz_z_zz, g_z_xz_z_xx, g_z_xz_z_xy, g_z_xz_z_xz, g_z_xz_z_yy, g_z_xz_z_yz, g_z_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_xz_z_xx[i] = -2.0 * g_z_xz_z_xx[i] * a_exp + 4.0 * g_xxz_xz_z_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xz_z_xy[i] = -2.0 * g_z_xz_z_xy[i] * a_exp + 4.0 * g_xxz_xz_z_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xz_z_xz[i] = -2.0 * g_z_xz_z_xz[i] * a_exp + 4.0 * g_xxz_xz_z_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xz_z_yy[i] = -2.0 * g_z_xz_z_yy[i] * a_exp + 4.0 * g_xxz_xz_z_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xz_z_yz[i] = -2.0 * g_z_xz_z_yz[i] * a_exp + 4.0 * g_xxz_xz_z_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xz_z_zz[i] = -2.0 * g_z_xz_z_zz[i] * a_exp + 4.0 * g_xxz_xz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_xx_0_0_0_z_yy_x_xx, g_xx_0_0_0_z_yy_x_xy, g_xx_0_0_0_z_yy_x_xz, g_xx_0_0_0_z_yy_x_yy, g_xx_0_0_0_z_yy_x_yz, g_xx_0_0_0_z_yy_x_zz, g_xxz_yy_x_xx, g_xxz_yy_x_xy, g_xxz_yy_x_xz, g_xxz_yy_x_yy, g_xxz_yy_x_yz, g_xxz_yy_x_zz, g_z_yy_x_xx, g_z_yy_x_xy, g_z_yy_x_xz, g_z_yy_x_yy, g_z_yy_x_yz, g_z_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_yy_x_xx[i] = -2.0 * g_z_yy_x_xx[i] * a_exp + 4.0 * g_xxz_yy_x_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yy_x_xy[i] = -2.0 * g_z_yy_x_xy[i] * a_exp + 4.0 * g_xxz_yy_x_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yy_x_xz[i] = -2.0 * g_z_yy_x_xz[i] * a_exp + 4.0 * g_xxz_yy_x_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yy_x_yy[i] = -2.0 * g_z_yy_x_yy[i] * a_exp + 4.0 * g_xxz_yy_x_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yy_x_yz[i] = -2.0 * g_z_yy_x_yz[i] * a_exp + 4.0 * g_xxz_yy_x_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yy_x_zz[i] = -2.0 * g_z_yy_x_zz[i] * a_exp + 4.0 * g_xxz_yy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_xx_0_0_0_z_yy_y_xx, g_xx_0_0_0_z_yy_y_xy, g_xx_0_0_0_z_yy_y_xz, g_xx_0_0_0_z_yy_y_yy, g_xx_0_0_0_z_yy_y_yz, g_xx_0_0_0_z_yy_y_zz, g_xxz_yy_y_xx, g_xxz_yy_y_xy, g_xxz_yy_y_xz, g_xxz_yy_y_yy, g_xxz_yy_y_yz, g_xxz_yy_y_zz, g_z_yy_y_xx, g_z_yy_y_xy, g_z_yy_y_xz, g_z_yy_y_yy, g_z_yy_y_yz, g_z_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_yy_y_xx[i] = -2.0 * g_z_yy_y_xx[i] * a_exp + 4.0 * g_xxz_yy_y_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yy_y_xy[i] = -2.0 * g_z_yy_y_xy[i] * a_exp + 4.0 * g_xxz_yy_y_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yy_y_xz[i] = -2.0 * g_z_yy_y_xz[i] * a_exp + 4.0 * g_xxz_yy_y_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yy_y_yy[i] = -2.0 * g_z_yy_y_yy[i] * a_exp + 4.0 * g_xxz_yy_y_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yy_y_yz[i] = -2.0 * g_z_yy_y_yz[i] * a_exp + 4.0 * g_xxz_yy_y_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yy_y_zz[i] = -2.0 * g_z_yy_y_zz[i] * a_exp + 4.0 * g_xxz_yy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_xx_0_0_0_z_yy_z_xx, g_xx_0_0_0_z_yy_z_xy, g_xx_0_0_0_z_yy_z_xz, g_xx_0_0_0_z_yy_z_yy, g_xx_0_0_0_z_yy_z_yz, g_xx_0_0_0_z_yy_z_zz, g_xxz_yy_z_xx, g_xxz_yy_z_xy, g_xxz_yy_z_xz, g_xxz_yy_z_yy, g_xxz_yy_z_yz, g_xxz_yy_z_zz, g_z_yy_z_xx, g_z_yy_z_xy, g_z_yy_z_xz, g_z_yy_z_yy, g_z_yy_z_yz, g_z_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_yy_z_xx[i] = -2.0 * g_z_yy_z_xx[i] * a_exp + 4.0 * g_xxz_yy_z_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yy_z_xy[i] = -2.0 * g_z_yy_z_xy[i] * a_exp + 4.0 * g_xxz_yy_z_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yy_z_xz[i] = -2.0 * g_z_yy_z_xz[i] * a_exp + 4.0 * g_xxz_yy_z_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yy_z_yy[i] = -2.0 * g_z_yy_z_yy[i] * a_exp + 4.0 * g_xxz_yy_z_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yy_z_yz[i] = -2.0 * g_z_yy_z_yz[i] * a_exp + 4.0 * g_xxz_yy_z_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yy_z_zz[i] = -2.0 * g_z_yy_z_zz[i] * a_exp + 4.0 * g_xxz_yy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_xx_0_0_0_z_yz_x_xx, g_xx_0_0_0_z_yz_x_xy, g_xx_0_0_0_z_yz_x_xz, g_xx_0_0_0_z_yz_x_yy, g_xx_0_0_0_z_yz_x_yz, g_xx_0_0_0_z_yz_x_zz, g_xxz_yz_x_xx, g_xxz_yz_x_xy, g_xxz_yz_x_xz, g_xxz_yz_x_yy, g_xxz_yz_x_yz, g_xxz_yz_x_zz, g_z_yz_x_xx, g_z_yz_x_xy, g_z_yz_x_xz, g_z_yz_x_yy, g_z_yz_x_yz, g_z_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_yz_x_xx[i] = -2.0 * g_z_yz_x_xx[i] * a_exp + 4.0 * g_xxz_yz_x_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yz_x_xy[i] = -2.0 * g_z_yz_x_xy[i] * a_exp + 4.0 * g_xxz_yz_x_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yz_x_xz[i] = -2.0 * g_z_yz_x_xz[i] * a_exp + 4.0 * g_xxz_yz_x_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yz_x_yy[i] = -2.0 * g_z_yz_x_yy[i] * a_exp + 4.0 * g_xxz_yz_x_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yz_x_yz[i] = -2.0 * g_z_yz_x_yz[i] * a_exp + 4.0 * g_xxz_yz_x_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yz_x_zz[i] = -2.0 * g_z_yz_x_zz[i] * a_exp + 4.0 * g_xxz_yz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_xx_0_0_0_z_yz_y_xx, g_xx_0_0_0_z_yz_y_xy, g_xx_0_0_0_z_yz_y_xz, g_xx_0_0_0_z_yz_y_yy, g_xx_0_0_0_z_yz_y_yz, g_xx_0_0_0_z_yz_y_zz, g_xxz_yz_y_xx, g_xxz_yz_y_xy, g_xxz_yz_y_xz, g_xxz_yz_y_yy, g_xxz_yz_y_yz, g_xxz_yz_y_zz, g_z_yz_y_xx, g_z_yz_y_xy, g_z_yz_y_xz, g_z_yz_y_yy, g_z_yz_y_yz, g_z_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_yz_y_xx[i] = -2.0 * g_z_yz_y_xx[i] * a_exp + 4.0 * g_xxz_yz_y_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yz_y_xy[i] = -2.0 * g_z_yz_y_xy[i] * a_exp + 4.0 * g_xxz_yz_y_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yz_y_xz[i] = -2.0 * g_z_yz_y_xz[i] * a_exp + 4.0 * g_xxz_yz_y_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yz_y_yy[i] = -2.0 * g_z_yz_y_yy[i] * a_exp + 4.0 * g_xxz_yz_y_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yz_y_yz[i] = -2.0 * g_z_yz_y_yz[i] * a_exp + 4.0 * g_xxz_yz_y_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yz_y_zz[i] = -2.0 * g_z_yz_y_zz[i] * a_exp + 4.0 * g_xxz_yz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_xx_0_0_0_z_yz_z_xx, g_xx_0_0_0_z_yz_z_xy, g_xx_0_0_0_z_yz_z_xz, g_xx_0_0_0_z_yz_z_yy, g_xx_0_0_0_z_yz_z_yz, g_xx_0_0_0_z_yz_z_zz, g_xxz_yz_z_xx, g_xxz_yz_z_xy, g_xxz_yz_z_xz, g_xxz_yz_z_yy, g_xxz_yz_z_yz, g_xxz_yz_z_zz, g_z_yz_z_xx, g_z_yz_z_xy, g_z_yz_z_xz, g_z_yz_z_yy, g_z_yz_z_yz, g_z_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_yz_z_xx[i] = -2.0 * g_z_yz_z_xx[i] * a_exp + 4.0 * g_xxz_yz_z_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yz_z_xy[i] = -2.0 * g_z_yz_z_xy[i] * a_exp + 4.0 * g_xxz_yz_z_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yz_z_xz[i] = -2.0 * g_z_yz_z_xz[i] * a_exp + 4.0 * g_xxz_yz_z_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yz_z_yy[i] = -2.0 * g_z_yz_z_yy[i] * a_exp + 4.0 * g_xxz_yz_z_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yz_z_yz[i] = -2.0 * g_z_yz_z_yz[i] * a_exp + 4.0 * g_xxz_yz_z_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yz_z_zz[i] = -2.0 * g_z_yz_z_zz[i] * a_exp + 4.0 * g_xxz_yz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_xx_0_0_0_z_zz_x_xx, g_xx_0_0_0_z_zz_x_xy, g_xx_0_0_0_z_zz_x_xz, g_xx_0_0_0_z_zz_x_yy, g_xx_0_0_0_z_zz_x_yz, g_xx_0_0_0_z_zz_x_zz, g_xxz_zz_x_xx, g_xxz_zz_x_xy, g_xxz_zz_x_xz, g_xxz_zz_x_yy, g_xxz_zz_x_yz, g_xxz_zz_x_zz, g_z_zz_x_xx, g_z_zz_x_xy, g_z_zz_x_xz, g_z_zz_x_yy, g_z_zz_x_yz, g_z_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_zz_x_xx[i] = -2.0 * g_z_zz_x_xx[i] * a_exp + 4.0 * g_xxz_zz_x_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_zz_x_xy[i] = -2.0 * g_z_zz_x_xy[i] * a_exp + 4.0 * g_xxz_zz_x_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_zz_x_xz[i] = -2.0 * g_z_zz_x_xz[i] * a_exp + 4.0 * g_xxz_zz_x_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_zz_x_yy[i] = -2.0 * g_z_zz_x_yy[i] * a_exp + 4.0 * g_xxz_zz_x_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_zz_x_yz[i] = -2.0 * g_z_zz_x_yz[i] * a_exp + 4.0 * g_xxz_zz_x_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_zz_x_zz[i] = -2.0 * g_z_zz_x_zz[i] * a_exp + 4.0 * g_xxz_zz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_xx_0_0_0_z_zz_y_xx, g_xx_0_0_0_z_zz_y_xy, g_xx_0_0_0_z_zz_y_xz, g_xx_0_0_0_z_zz_y_yy, g_xx_0_0_0_z_zz_y_yz, g_xx_0_0_0_z_zz_y_zz, g_xxz_zz_y_xx, g_xxz_zz_y_xy, g_xxz_zz_y_xz, g_xxz_zz_y_yy, g_xxz_zz_y_yz, g_xxz_zz_y_zz, g_z_zz_y_xx, g_z_zz_y_xy, g_z_zz_y_xz, g_z_zz_y_yy, g_z_zz_y_yz, g_z_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_zz_y_xx[i] = -2.0 * g_z_zz_y_xx[i] * a_exp + 4.0 * g_xxz_zz_y_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_zz_y_xy[i] = -2.0 * g_z_zz_y_xy[i] * a_exp + 4.0 * g_xxz_zz_y_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_zz_y_xz[i] = -2.0 * g_z_zz_y_xz[i] * a_exp + 4.0 * g_xxz_zz_y_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_zz_y_yy[i] = -2.0 * g_z_zz_y_yy[i] * a_exp + 4.0 * g_xxz_zz_y_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_zz_y_yz[i] = -2.0 * g_z_zz_y_yz[i] * a_exp + 4.0 * g_xxz_zz_y_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_zz_y_zz[i] = -2.0 * g_z_zz_y_zz[i] * a_exp + 4.0 * g_xxz_zz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_xx_0_0_0_z_zz_z_xx, g_xx_0_0_0_z_zz_z_xy, g_xx_0_0_0_z_zz_z_xz, g_xx_0_0_0_z_zz_z_yy, g_xx_0_0_0_z_zz_z_yz, g_xx_0_0_0_z_zz_z_zz, g_xxz_zz_z_xx, g_xxz_zz_z_xy, g_xxz_zz_z_xz, g_xxz_zz_z_yy, g_xxz_zz_z_yz, g_xxz_zz_z_zz, g_z_zz_z_xx, g_z_zz_z_xy, g_z_zz_z_xz, g_z_zz_z_yy, g_z_zz_z_yz, g_z_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_zz_z_xx[i] = -2.0 * g_z_zz_z_xx[i] * a_exp + 4.0 * g_xxz_zz_z_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_z_zz_z_xy[i] = -2.0 * g_z_zz_z_xy[i] * a_exp + 4.0 * g_xxz_zz_z_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_zz_z_xz[i] = -2.0 * g_z_zz_z_xz[i] * a_exp + 4.0 * g_xxz_zz_z_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_zz_z_yy[i] = -2.0 * g_z_zz_z_yy[i] * a_exp + 4.0 * g_xxz_zz_z_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_z_zz_z_yz[i] = -2.0 * g_z_zz_z_yz[i] * a_exp + 4.0 * g_xxz_zz_z_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_z_zz_z_zz[i] = -2.0 * g_z_zz_z_zz[i] * a_exp + 4.0 * g_xxz_zz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (324-330)

    #pragma omp simd aligned(g_xxy_xx_x_xx, g_xxy_xx_x_xy, g_xxy_xx_x_xz, g_xxy_xx_x_yy, g_xxy_xx_x_yz, g_xxy_xx_x_zz, g_xy_0_0_0_x_xx_x_xx, g_xy_0_0_0_x_xx_x_xy, g_xy_0_0_0_x_xx_x_xz, g_xy_0_0_0_x_xx_x_yy, g_xy_0_0_0_x_xx_x_yz, g_xy_0_0_0_x_xx_x_zz, g_y_xx_x_xx, g_y_xx_x_xy, g_y_xx_x_xz, g_y_xx_x_yy, g_y_xx_x_yz, g_y_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_xx_x_xx[i] = -2.0 * g_y_xx_x_xx[i] * a_exp + 4.0 * g_xxy_xx_x_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xx_x_xy[i] = -2.0 * g_y_xx_x_xy[i] * a_exp + 4.0 * g_xxy_xx_x_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xx_x_xz[i] = -2.0 * g_y_xx_x_xz[i] * a_exp + 4.0 * g_xxy_xx_x_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xx_x_yy[i] = -2.0 * g_y_xx_x_yy[i] * a_exp + 4.0 * g_xxy_xx_x_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xx_x_yz[i] = -2.0 * g_y_xx_x_yz[i] * a_exp + 4.0 * g_xxy_xx_x_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xx_x_zz[i] = -2.0 * g_y_xx_x_zz[i] * a_exp + 4.0 * g_xxy_xx_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (330-336)

    #pragma omp simd aligned(g_xxy_xx_y_xx, g_xxy_xx_y_xy, g_xxy_xx_y_xz, g_xxy_xx_y_yy, g_xxy_xx_y_yz, g_xxy_xx_y_zz, g_xy_0_0_0_x_xx_y_xx, g_xy_0_0_0_x_xx_y_xy, g_xy_0_0_0_x_xx_y_xz, g_xy_0_0_0_x_xx_y_yy, g_xy_0_0_0_x_xx_y_yz, g_xy_0_0_0_x_xx_y_zz, g_y_xx_y_xx, g_y_xx_y_xy, g_y_xx_y_xz, g_y_xx_y_yy, g_y_xx_y_yz, g_y_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_xx_y_xx[i] = -2.0 * g_y_xx_y_xx[i] * a_exp + 4.0 * g_xxy_xx_y_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xx_y_xy[i] = -2.0 * g_y_xx_y_xy[i] * a_exp + 4.0 * g_xxy_xx_y_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xx_y_xz[i] = -2.0 * g_y_xx_y_xz[i] * a_exp + 4.0 * g_xxy_xx_y_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xx_y_yy[i] = -2.0 * g_y_xx_y_yy[i] * a_exp + 4.0 * g_xxy_xx_y_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xx_y_yz[i] = -2.0 * g_y_xx_y_yz[i] * a_exp + 4.0 * g_xxy_xx_y_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xx_y_zz[i] = -2.0 * g_y_xx_y_zz[i] * a_exp + 4.0 * g_xxy_xx_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (336-342)

    #pragma omp simd aligned(g_xxy_xx_z_xx, g_xxy_xx_z_xy, g_xxy_xx_z_xz, g_xxy_xx_z_yy, g_xxy_xx_z_yz, g_xxy_xx_z_zz, g_xy_0_0_0_x_xx_z_xx, g_xy_0_0_0_x_xx_z_xy, g_xy_0_0_0_x_xx_z_xz, g_xy_0_0_0_x_xx_z_yy, g_xy_0_0_0_x_xx_z_yz, g_xy_0_0_0_x_xx_z_zz, g_y_xx_z_xx, g_y_xx_z_xy, g_y_xx_z_xz, g_y_xx_z_yy, g_y_xx_z_yz, g_y_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_xx_z_xx[i] = -2.0 * g_y_xx_z_xx[i] * a_exp + 4.0 * g_xxy_xx_z_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xx_z_xy[i] = -2.0 * g_y_xx_z_xy[i] * a_exp + 4.0 * g_xxy_xx_z_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xx_z_xz[i] = -2.0 * g_y_xx_z_xz[i] * a_exp + 4.0 * g_xxy_xx_z_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xx_z_yy[i] = -2.0 * g_y_xx_z_yy[i] * a_exp + 4.0 * g_xxy_xx_z_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xx_z_yz[i] = -2.0 * g_y_xx_z_yz[i] * a_exp + 4.0 * g_xxy_xx_z_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xx_z_zz[i] = -2.0 * g_y_xx_z_zz[i] * a_exp + 4.0 * g_xxy_xx_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (342-348)

    #pragma omp simd aligned(g_xxy_xy_x_xx, g_xxy_xy_x_xy, g_xxy_xy_x_xz, g_xxy_xy_x_yy, g_xxy_xy_x_yz, g_xxy_xy_x_zz, g_xy_0_0_0_x_xy_x_xx, g_xy_0_0_0_x_xy_x_xy, g_xy_0_0_0_x_xy_x_xz, g_xy_0_0_0_x_xy_x_yy, g_xy_0_0_0_x_xy_x_yz, g_xy_0_0_0_x_xy_x_zz, g_y_xy_x_xx, g_y_xy_x_xy, g_y_xy_x_xz, g_y_xy_x_yy, g_y_xy_x_yz, g_y_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_xy_x_xx[i] = -2.0 * g_y_xy_x_xx[i] * a_exp + 4.0 * g_xxy_xy_x_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xy_x_xy[i] = -2.0 * g_y_xy_x_xy[i] * a_exp + 4.0 * g_xxy_xy_x_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xy_x_xz[i] = -2.0 * g_y_xy_x_xz[i] * a_exp + 4.0 * g_xxy_xy_x_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xy_x_yy[i] = -2.0 * g_y_xy_x_yy[i] * a_exp + 4.0 * g_xxy_xy_x_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xy_x_yz[i] = -2.0 * g_y_xy_x_yz[i] * a_exp + 4.0 * g_xxy_xy_x_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xy_x_zz[i] = -2.0 * g_y_xy_x_zz[i] * a_exp + 4.0 * g_xxy_xy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (348-354)

    #pragma omp simd aligned(g_xxy_xy_y_xx, g_xxy_xy_y_xy, g_xxy_xy_y_xz, g_xxy_xy_y_yy, g_xxy_xy_y_yz, g_xxy_xy_y_zz, g_xy_0_0_0_x_xy_y_xx, g_xy_0_0_0_x_xy_y_xy, g_xy_0_0_0_x_xy_y_xz, g_xy_0_0_0_x_xy_y_yy, g_xy_0_0_0_x_xy_y_yz, g_xy_0_0_0_x_xy_y_zz, g_y_xy_y_xx, g_y_xy_y_xy, g_y_xy_y_xz, g_y_xy_y_yy, g_y_xy_y_yz, g_y_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_xy_y_xx[i] = -2.0 * g_y_xy_y_xx[i] * a_exp + 4.0 * g_xxy_xy_y_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xy_y_xy[i] = -2.0 * g_y_xy_y_xy[i] * a_exp + 4.0 * g_xxy_xy_y_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xy_y_xz[i] = -2.0 * g_y_xy_y_xz[i] * a_exp + 4.0 * g_xxy_xy_y_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xy_y_yy[i] = -2.0 * g_y_xy_y_yy[i] * a_exp + 4.0 * g_xxy_xy_y_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xy_y_yz[i] = -2.0 * g_y_xy_y_yz[i] * a_exp + 4.0 * g_xxy_xy_y_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xy_y_zz[i] = -2.0 * g_y_xy_y_zz[i] * a_exp + 4.0 * g_xxy_xy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (354-360)

    #pragma omp simd aligned(g_xxy_xy_z_xx, g_xxy_xy_z_xy, g_xxy_xy_z_xz, g_xxy_xy_z_yy, g_xxy_xy_z_yz, g_xxy_xy_z_zz, g_xy_0_0_0_x_xy_z_xx, g_xy_0_0_0_x_xy_z_xy, g_xy_0_0_0_x_xy_z_xz, g_xy_0_0_0_x_xy_z_yy, g_xy_0_0_0_x_xy_z_yz, g_xy_0_0_0_x_xy_z_zz, g_y_xy_z_xx, g_y_xy_z_xy, g_y_xy_z_xz, g_y_xy_z_yy, g_y_xy_z_yz, g_y_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_xy_z_xx[i] = -2.0 * g_y_xy_z_xx[i] * a_exp + 4.0 * g_xxy_xy_z_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xy_z_xy[i] = -2.0 * g_y_xy_z_xy[i] * a_exp + 4.0 * g_xxy_xy_z_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xy_z_xz[i] = -2.0 * g_y_xy_z_xz[i] * a_exp + 4.0 * g_xxy_xy_z_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xy_z_yy[i] = -2.0 * g_y_xy_z_yy[i] * a_exp + 4.0 * g_xxy_xy_z_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xy_z_yz[i] = -2.0 * g_y_xy_z_yz[i] * a_exp + 4.0 * g_xxy_xy_z_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xy_z_zz[i] = -2.0 * g_y_xy_z_zz[i] * a_exp + 4.0 * g_xxy_xy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (360-366)

    #pragma omp simd aligned(g_xxy_xz_x_xx, g_xxy_xz_x_xy, g_xxy_xz_x_xz, g_xxy_xz_x_yy, g_xxy_xz_x_yz, g_xxy_xz_x_zz, g_xy_0_0_0_x_xz_x_xx, g_xy_0_0_0_x_xz_x_xy, g_xy_0_0_0_x_xz_x_xz, g_xy_0_0_0_x_xz_x_yy, g_xy_0_0_0_x_xz_x_yz, g_xy_0_0_0_x_xz_x_zz, g_y_xz_x_xx, g_y_xz_x_xy, g_y_xz_x_xz, g_y_xz_x_yy, g_y_xz_x_yz, g_y_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_xz_x_xx[i] = -2.0 * g_y_xz_x_xx[i] * a_exp + 4.0 * g_xxy_xz_x_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xz_x_xy[i] = -2.0 * g_y_xz_x_xy[i] * a_exp + 4.0 * g_xxy_xz_x_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xz_x_xz[i] = -2.0 * g_y_xz_x_xz[i] * a_exp + 4.0 * g_xxy_xz_x_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xz_x_yy[i] = -2.0 * g_y_xz_x_yy[i] * a_exp + 4.0 * g_xxy_xz_x_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xz_x_yz[i] = -2.0 * g_y_xz_x_yz[i] * a_exp + 4.0 * g_xxy_xz_x_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xz_x_zz[i] = -2.0 * g_y_xz_x_zz[i] * a_exp + 4.0 * g_xxy_xz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (366-372)

    #pragma omp simd aligned(g_xxy_xz_y_xx, g_xxy_xz_y_xy, g_xxy_xz_y_xz, g_xxy_xz_y_yy, g_xxy_xz_y_yz, g_xxy_xz_y_zz, g_xy_0_0_0_x_xz_y_xx, g_xy_0_0_0_x_xz_y_xy, g_xy_0_0_0_x_xz_y_xz, g_xy_0_0_0_x_xz_y_yy, g_xy_0_0_0_x_xz_y_yz, g_xy_0_0_0_x_xz_y_zz, g_y_xz_y_xx, g_y_xz_y_xy, g_y_xz_y_xz, g_y_xz_y_yy, g_y_xz_y_yz, g_y_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_xz_y_xx[i] = -2.0 * g_y_xz_y_xx[i] * a_exp + 4.0 * g_xxy_xz_y_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xz_y_xy[i] = -2.0 * g_y_xz_y_xy[i] * a_exp + 4.0 * g_xxy_xz_y_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xz_y_xz[i] = -2.0 * g_y_xz_y_xz[i] * a_exp + 4.0 * g_xxy_xz_y_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xz_y_yy[i] = -2.0 * g_y_xz_y_yy[i] * a_exp + 4.0 * g_xxy_xz_y_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xz_y_yz[i] = -2.0 * g_y_xz_y_yz[i] * a_exp + 4.0 * g_xxy_xz_y_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xz_y_zz[i] = -2.0 * g_y_xz_y_zz[i] * a_exp + 4.0 * g_xxy_xz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (372-378)

    #pragma omp simd aligned(g_xxy_xz_z_xx, g_xxy_xz_z_xy, g_xxy_xz_z_xz, g_xxy_xz_z_yy, g_xxy_xz_z_yz, g_xxy_xz_z_zz, g_xy_0_0_0_x_xz_z_xx, g_xy_0_0_0_x_xz_z_xy, g_xy_0_0_0_x_xz_z_xz, g_xy_0_0_0_x_xz_z_yy, g_xy_0_0_0_x_xz_z_yz, g_xy_0_0_0_x_xz_z_zz, g_y_xz_z_xx, g_y_xz_z_xy, g_y_xz_z_xz, g_y_xz_z_yy, g_y_xz_z_yz, g_y_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_xz_z_xx[i] = -2.0 * g_y_xz_z_xx[i] * a_exp + 4.0 * g_xxy_xz_z_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xz_z_xy[i] = -2.0 * g_y_xz_z_xy[i] * a_exp + 4.0 * g_xxy_xz_z_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xz_z_xz[i] = -2.0 * g_y_xz_z_xz[i] * a_exp + 4.0 * g_xxy_xz_z_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xz_z_yy[i] = -2.0 * g_y_xz_z_yy[i] * a_exp + 4.0 * g_xxy_xz_z_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xz_z_yz[i] = -2.0 * g_y_xz_z_yz[i] * a_exp + 4.0 * g_xxy_xz_z_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xz_z_zz[i] = -2.0 * g_y_xz_z_zz[i] * a_exp + 4.0 * g_xxy_xz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (378-384)

    #pragma omp simd aligned(g_xxy_yy_x_xx, g_xxy_yy_x_xy, g_xxy_yy_x_xz, g_xxy_yy_x_yy, g_xxy_yy_x_yz, g_xxy_yy_x_zz, g_xy_0_0_0_x_yy_x_xx, g_xy_0_0_0_x_yy_x_xy, g_xy_0_0_0_x_yy_x_xz, g_xy_0_0_0_x_yy_x_yy, g_xy_0_0_0_x_yy_x_yz, g_xy_0_0_0_x_yy_x_zz, g_y_yy_x_xx, g_y_yy_x_xy, g_y_yy_x_xz, g_y_yy_x_yy, g_y_yy_x_yz, g_y_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_yy_x_xx[i] = -2.0 * g_y_yy_x_xx[i] * a_exp + 4.0 * g_xxy_yy_x_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yy_x_xy[i] = -2.0 * g_y_yy_x_xy[i] * a_exp + 4.0 * g_xxy_yy_x_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yy_x_xz[i] = -2.0 * g_y_yy_x_xz[i] * a_exp + 4.0 * g_xxy_yy_x_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yy_x_yy[i] = -2.0 * g_y_yy_x_yy[i] * a_exp + 4.0 * g_xxy_yy_x_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yy_x_yz[i] = -2.0 * g_y_yy_x_yz[i] * a_exp + 4.0 * g_xxy_yy_x_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yy_x_zz[i] = -2.0 * g_y_yy_x_zz[i] * a_exp + 4.0 * g_xxy_yy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (384-390)

    #pragma omp simd aligned(g_xxy_yy_y_xx, g_xxy_yy_y_xy, g_xxy_yy_y_xz, g_xxy_yy_y_yy, g_xxy_yy_y_yz, g_xxy_yy_y_zz, g_xy_0_0_0_x_yy_y_xx, g_xy_0_0_0_x_yy_y_xy, g_xy_0_0_0_x_yy_y_xz, g_xy_0_0_0_x_yy_y_yy, g_xy_0_0_0_x_yy_y_yz, g_xy_0_0_0_x_yy_y_zz, g_y_yy_y_xx, g_y_yy_y_xy, g_y_yy_y_xz, g_y_yy_y_yy, g_y_yy_y_yz, g_y_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_yy_y_xx[i] = -2.0 * g_y_yy_y_xx[i] * a_exp + 4.0 * g_xxy_yy_y_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yy_y_xy[i] = -2.0 * g_y_yy_y_xy[i] * a_exp + 4.0 * g_xxy_yy_y_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yy_y_xz[i] = -2.0 * g_y_yy_y_xz[i] * a_exp + 4.0 * g_xxy_yy_y_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yy_y_yy[i] = -2.0 * g_y_yy_y_yy[i] * a_exp + 4.0 * g_xxy_yy_y_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yy_y_yz[i] = -2.0 * g_y_yy_y_yz[i] * a_exp + 4.0 * g_xxy_yy_y_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yy_y_zz[i] = -2.0 * g_y_yy_y_zz[i] * a_exp + 4.0 * g_xxy_yy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (390-396)

    #pragma omp simd aligned(g_xxy_yy_z_xx, g_xxy_yy_z_xy, g_xxy_yy_z_xz, g_xxy_yy_z_yy, g_xxy_yy_z_yz, g_xxy_yy_z_zz, g_xy_0_0_0_x_yy_z_xx, g_xy_0_0_0_x_yy_z_xy, g_xy_0_0_0_x_yy_z_xz, g_xy_0_0_0_x_yy_z_yy, g_xy_0_0_0_x_yy_z_yz, g_xy_0_0_0_x_yy_z_zz, g_y_yy_z_xx, g_y_yy_z_xy, g_y_yy_z_xz, g_y_yy_z_yy, g_y_yy_z_yz, g_y_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_yy_z_xx[i] = -2.0 * g_y_yy_z_xx[i] * a_exp + 4.0 * g_xxy_yy_z_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yy_z_xy[i] = -2.0 * g_y_yy_z_xy[i] * a_exp + 4.0 * g_xxy_yy_z_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yy_z_xz[i] = -2.0 * g_y_yy_z_xz[i] * a_exp + 4.0 * g_xxy_yy_z_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yy_z_yy[i] = -2.0 * g_y_yy_z_yy[i] * a_exp + 4.0 * g_xxy_yy_z_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yy_z_yz[i] = -2.0 * g_y_yy_z_yz[i] * a_exp + 4.0 * g_xxy_yy_z_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yy_z_zz[i] = -2.0 * g_y_yy_z_zz[i] * a_exp + 4.0 * g_xxy_yy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (396-402)

    #pragma omp simd aligned(g_xxy_yz_x_xx, g_xxy_yz_x_xy, g_xxy_yz_x_xz, g_xxy_yz_x_yy, g_xxy_yz_x_yz, g_xxy_yz_x_zz, g_xy_0_0_0_x_yz_x_xx, g_xy_0_0_0_x_yz_x_xy, g_xy_0_0_0_x_yz_x_xz, g_xy_0_0_0_x_yz_x_yy, g_xy_0_0_0_x_yz_x_yz, g_xy_0_0_0_x_yz_x_zz, g_y_yz_x_xx, g_y_yz_x_xy, g_y_yz_x_xz, g_y_yz_x_yy, g_y_yz_x_yz, g_y_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_yz_x_xx[i] = -2.0 * g_y_yz_x_xx[i] * a_exp + 4.0 * g_xxy_yz_x_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yz_x_xy[i] = -2.0 * g_y_yz_x_xy[i] * a_exp + 4.0 * g_xxy_yz_x_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yz_x_xz[i] = -2.0 * g_y_yz_x_xz[i] * a_exp + 4.0 * g_xxy_yz_x_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yz_x_yy[i] = -2.0 * g_y_yz_x_yy[i] * a_exp + 4.0 * g_xxy_yz_x_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yz_x_yz[i] = -2.0 * g_y_yz_x_yz[i] * a_exp + 4.0 * g_xxy_yz_x_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yz_x_zz[i] = -2.0 * g_y_yz_x_zz[i] * a_exp + 4.0 * g_xxy_yz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (402-408)

    #pragma omp simd aligned(g_xxy_yz_y_xx, g_xxy_yz_y_xy, g_xxy_yz_y_xz, g_xxy_yz_y_yy, g_xxy_yz_y_yz, g_xxy_yz_y_zz, g_xy_0_0_0_x_yz_y_xx, g_xy_0_0_0_x_yz_y_xy, g_xy_0_0_0_x_yz_y_xz, g_xy_0_0_0_x_yz_y_yy, g_xy_0_0_0_x_yz_y_yz, g_xy_0_0_0_x_yz_y_zz, g_y_yz_y_xx, g_y_yz_y_xy, g_y_yz_y_xz, g_y_yz_y_yy, g_y_yz_y_yz, g_y_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_yz_y_xx[i] = -2.0 * g_y_yz_y_xx[i] * a_exp + 4.0 * g_xxy_yz_y_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yz_y_xy[i] = -2.0 * g_y_yz_y_xy[i] * a_exp + 4.0 * g_xxy_yz_y_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yz_y_xz[i] = -2.0 * g_y_yz_y_xz[i] * a_exp + 4.0 * g_xxy_yz_y_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yz_y_yy[i] = -2.0 * g_y_yz_y_yy[i] * a_exp + 4.0 * g_xxy_yz_y_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yz_y_yz[i] = -2.0 * g_y_yz_y_yz[i] * a_exp + 4.0 * g_xxy_yz_y_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yz_y_zz[i] = -2.0 * g_y_yz_y_zz[i] * a_exp + 4.0 * g_xxy_yz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (408-414)

    #pragma omp simd aligned(g_xxy_yz_z_xx, g_xxy_yz_z_xy, g_xxy_yz_z_xz, g_xxy_yz_z_yy, g_xxy_yz_z_yz, g_xxy_yz_z_zz, g_xy_0_0_0_x_yz_z_xx, g_xy_0_0_0_x_yz_z_xy, g_xy_0_0_0_x_yz_z_xz, g_xy_0_0_0_x_yz_z_yy, g_xy_0_0_0_x_yz_z_yz, g_xy_0_0_0_x_yz_z_zz, g_y_yz_z_xx, g_y_yz_z_xy, g_y_yz_z_xz, g_y_yz_z_yy, g_y_yz_z_yz, g_y_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_yz_z_xx[i] = -2.0 * g_y_yz_z_xx[i] * a_exp + 4.0 * g_xxy_yz_z_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yz_z_xy[i] = -2.0 * g_y_yz_z_xy[i] * a_exp + 4.0 * g_xxy_yz_z_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yz_z_xz[i] = -2.0 * g_y_yz_z_xz[i] * a_exp + 4.0 * g_xxy_yz_z_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yz_z_yy[i] = -2.0 * g_y_yz_z_yy[i] * a_exp + 4.0 * g_xxy_yz_z_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yz_z_yz[i] = -2.0 * g_y_yz_z_yz[i] * a_exp + 4.0 * g_xxy_yz_z_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yz_z_zz[i] = -2.0 * g_y_yz_z_zz[i] * a_exp + 4.0 * g_xxy_yz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (414-420)

    #pragma omp simd aligned(g_xxy_zz_x_xx, g_xxy_zz_x_xy, g_xxy_zz_x_xz, g_xxy_zz_x_yy, g_xxy_zz_x_yz, g_xxy_zz_x_zz, g_xy_0_0_0_x_zz_x_xx, g_xy_0_0_0_x_zz_x_xy, g_xy_0_0_0_x_zz_x_xz, g_xy_0_0_0_x_zz_x_yy, g_xy_0_0_0_x_zz_x_yz, g_xy_0_0_0_x_zz_x_zz, g_y_zz_x_xx, g_y_zz_x_xy, g_y_zz_x_xz, g_y_zz_x_yy, g_y_zz_x_yz, g_y_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_zz_x_xx[i] = -2.0 * g_y_zz_x_xx[i] * a_exp + 4.0 * g_xxy_zz_x_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_zz_x_xy[i] = -2.0 * g_y_zz_x_xy[i] * a_exp + 4.0 * g_xxy_zz_x_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_zz_x_xz[i] = -2.0 * g_y_zz_x_xz[i] * a_exp + 4.0 * g_xxy_zz_x_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_zz_x_yy[i] = -2.0 * g_y_zz_x_yy[i] * a_exp + 4.0 * g_xxy_zz_x_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_zz_x_yz[i] = -2.0 * g_y_zz_x_yz[i] * a_exp + 4.0 * g_xxy_zz_x_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_zz_x_zz[i] = -2.0 * g_y_zz_x_zz[i] * a_exp + 4.0 * g_xxy_zz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (420-426)

    #pragma omp simd aligned(g_xxy_zz_y_xx, g_xxy_zz_y_xy, g_xxy_zz_y_xz, g_xxy_zz_y_yy, g_xxy_zz_y_yz, g_xxy_zz_y_zz, g_xy_0_0_0_x_zz_y_xx, g_xy_0_0_0_x_zz_y_xy, g_xy_0_0_0_x_zz_y_xz, g_xy_0_0_0_x_zz_y_yy, g_xy_0_0_0_x_zz_y_yz, g_xy_0_0_0_x_zz_y_zz, g_y_zz_y_xx, g_y_zz_y_xy, g_y_zz_y_xz, g_y_zz_y_yy, g_y_zz_y_yz, g_y_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_zz_y_xx[i] = -2.0 * g_y_zz_y_xx[i] * a_exp + 4.0 * g_xxy_zz_y_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_zz_y_xy[i] = -2.0 * g_y_zz_y_xy[i] * a_exp + 4.0 * g_xxy_zz_y_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_zz_y_xz[i] = -2.0 * g_y_zz_y_xz[i] * a_exp + 4.0 * g_xxy_zz_y_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_zz_y_yy[i] = -2.0 * g_y_zz_y_yy[i] * a_exp + 4.0 * g_xxy_zz_y_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_zz_y_yz[i] = -2.0 * g_y_zz_y_yz[i] * a_exp + 4.0 * g_xxy_zz_y_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_zz_y_zz[i] = -2.0 * g_y_zz_y_zz[i] * a_exp + 4.0 * g_xxy_zz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (426-432)

    #pragma omp simd aligned(g_xxy_zz_z_xx, g_xxy_zz_z_xy, g_xxy_zz_z_xz, g_xxy_zz_z_yy, g_xxy_zz_z_yz, g_xxy_zz_z_zz, g_xy_0_0_0_x_zz_z_xx, g_xy_0_0_0_x_zz_z_xy, g_xy_0_0_0_x_zz_z_xz, g_xy_0_0_0_x_zz_z_yy, g_xy_0_0_0_x_zz_z_yz, g_xy_0_0_0_x_zz_z_zz, g_y_zz_z_xx, g_y_zz_z_xy, g_y_zz_z_xz, g_y_zz_z_yy, g_y_zz_z_yz, g_y_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_zz_z_xx[i] = -2.0 * g_y_zz_z_xx[i] * a_exp + 4.0 * g_xxy_zz_z_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_x_zz_z_xy[i] = -2.0 * g_y_zz_z_xy[i] * a_exp + 4.0 * g_xxy_zz_z_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_zz_z_xz[i] = -2.0 * g_y_zz_z_xz[i] * a_exp + 4.0 * g_xxy_zz_z_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_zz_z_yy[i] = -2.0 * g_y_zz_z_yy[i] * a_exp + 4.0 * g_xxy_zz_z_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_x_zz_z_yz[i] = -2.0 * g_y_zz_z_yz[i] * a_exp + 4.0 * g_xxy_zz_z_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_x_zz_z_zz[i] = -2.0 * g_y_zz_z_zz[i] * a_exp + 4.0 * g_xxy_zz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (432-438)

    #pragma omp simd aligned(g_x_xx_x_xx, g_x_xx_x_xy, g_x_xx_x_xz, g_x_xx_x_yy, g_x_xx_x_yz, g_x_xx_x_zz, g_xy_0_0_0_y_xx_x_xx, g_xy_0_0_0_y_xx_x_xy, g_xy_0_0_0_y_xx_x_xz, g_xy_0_0_0_y_xx_x_yy, g_xy_0_0_0_y_xx_x_yz, g_xy_0_0_0_y_xx_x_zz, g_xyy_xx_x_xx, g_xyy_xx_x_xy, g_xyy_xx_x_xz, g_xyy_xx_x_yy, g_xyy_xx_x_yz, g_xyy_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_xx_x_xx[i] = -2.0 * g_x_xx_x_xx[i] * a_exp + 4.0 * g_xyy_xx_x_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xx_x_xy[i] = -2.0 * g_x_xx_x_xy[i] * a_exp + 4.0 * g_xyy_xx_x_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xx_x_xz[i] = -2.0 * g_x_xx_x_xz[i] * a_exp + 4.0 * g_xyy_xx_x_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xx_x_yy[i] = -2.0 * g_x_xx_x_yy[i] * a_exp + 4.0 * g_xyy_xx_x_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xx_x_yz[i] = -2.0 * g_x_xx_x_yz[i] * a_exp + 4.0 * g_xyy_xx_x_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xx_x_zz[i] = -2.0 * g_x_xx_x_zz[i] * a_exp + 4.0 * g_xyy_xx_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (438-444)

    #pragma omp simd aligned(g_x_xx_y_xx, g_x_xx_y_xy, g_x_xx_y_xz, g_x_xx_y_yy, g_x_xx_y_yz, g_x_xx_y_zz, g_xy_0_0_0_y_xx_y_xx, g_xy_0_0_0_y_xx_y_xy, g_xy_0_0_0_y_xx_y_xz, g_xy_0_0_0_y_xx_y_yy, g_xy_0_0_0_y_xx_y_yz, g_xy_0_0_0_y_xx_y_zz, g_xyy_xx_y_xx, g_xyy_xx_y_xy, g_xyy_xx_y_xz, g_xyy_xx_y_yy, g_xyy_xx_y_yz, g_xyy_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_xx_y_xx[i] = -2.0 * g_x_xx_y_xx[i] * a_exp + 4.0 * g_xyy_xx_y_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xx_y_xy[i] = -2.0 * g_x_xx_y_xy[i] * a_exp + 4.0 * g_xyy_xx_y_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xx_y_xz[i] = -2.0 * g_x_xx_y_xz[i] * a_exp + 4.0 * g_xyy_xx_y_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xx_y_yy[i] = -2.0 * g_x_xx_y_yy[i] * a_exp + 4.0 * g_xyy_xx_y_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xx_y_yz[i] = -2.0 * g_x_xx_y_yz[i] * a_exp + 4.0 * g_xyy_xx_y_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xx_y_zz[i] = -2.0 * g_x_xx_y_zz[i] * a_exp + 4.0 * g_xyy_xx_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (444-450)

    #pragma omp simd aligned(g_x_xx_z_xx, g_x_xx_z_xy, g_x_xx_z_xz, g_x_xx_z_yy, g_x_xx_z_yz, g_x_xx_z_zz, g_xy_0_0_0_y_xx_z_xx, g_xy_0_0_0_y_xx_z_xy, g_xy_0_0_0_y_xx_z_xz, g_xy_0_0_0_y_xx_z_yy, g_xy_0_0_0_y_xx_z_yz, g_xy_0_0_0_y_xx_z_zz, g_xyy_xx_z_xx, g_xyy_xx_z_xy, g_xyy_xx_z_xz, g_xyy_xx_z_yy, g_xyy_xx_z_yz, g_xyy_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_xx_z_xx[i] = -2.0 * g_x_xx_z_xx[i] * a_exp + 4.0 * g_xyy_xx_z_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xx_z_xy[i] = -2.0 * g_x_xx_z_xy[i] * a_exp + 4.0 * g_xyy_xx_z_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xx_z_xz[i] = -2.0 * g_x_xx_z_xz[i] * a_exp + 4.0 * g_xyy_xx_z_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xx_z_yy[i] = -2.0 * g_x_xx_z_yy[i] * a_exp + 4.0 * g_xyy_xx_z_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xx_z_yz[i] = -2.0 * g_x_xx_z_yz[i] * a_exp + 4.0 * g_xyy_xx_z_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xx_z_zz[i] = -2.0 * g_x_xx_z_zz[i] * a_exp + 4.0 * g_xyy_xx_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (450-456)

    #pragma omp simd aligned(g_x_xy_x_xx, g_x_xy_x_xy, g_x_xy_x_xz, g_x_xy_x_yy, g_x_xy_x_yz, g_x_xy_x_zz, g_xy_0_0_0_y_xy_x_xx, g_xy_0_0_0_y_xy_x_xy, g_xy_0_0_0_y_xy_x_xz, g_xy_0_0_0_y_xy_x_yy, g_xy_0_0_0_y_xy_x_yz, g_xy_0_0_0_y_xy_x_zz, g_xyy_xy_x_xx, g_xyy_xy_x_xy, g_xyy_xy_x_xz, g_xyy_xy_x_yy, g_xyy_xy_x_yz, g_xyy_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_xy_x_xx[i] = -2.0 * g_x_xy_x_xx[i] * a_exp + 4.0 * g_xyy_xy_x_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xy_x_xy[i] = -2.0 * g_x_xy_x_xy[i] * a_exp + 4.0 * g_xyy_xy_x_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xy_x_xz[i] = -2.0 * g_x_xy_x_xz[i] * a_exp + 4.0 * g_xyy_xy_x_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xy_x_yy[i] = -2.0 * g_x_xy_x_yy[i] * a_exp + 4.0 * g_xyy_xy_x_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xy_x_yz[i] = -2.0 * g_x_xy_x_yz[i] * a_exp + 4.0 * g_xyy_xy_x_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xy_x_zz[i] = -2.0 * g_x_xy_x_zz[i] * a_exp + 4.0 * g_xyy_xy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (456-462)

    #pragma omp simd aligned(g_x_xy_y_xx, g_x_xy_y_xy, g_x_xy_y_xz, g_x_xy_y_yy, g_x_xy_y_yz, g_x_xy_y_zz, g_xy_0_0_0_y_xy_y_xx, g_xy_0_0_0_y_xy_y_xy, g_xy_0_0_0_y_xy_y_xz, g_xy_0_0_0_y_xy_y_yy, g_xy_0_0_0_y_xy_y_yz, g_xy_0_0_0_y_xy_y_zz, g_xyy_xy_y_xx, g_xyy_xy_y_xy, g_xyy_xy_y_xz, g_xyy_xy_y_yy, g_xyy_xy_y_yz, g_xyy_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_xy_y_xx[i] = -2.0 * g_x_xy_y_xx[i] * a_exp + 4.0 * g_xyy_xy_y_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xy_y_xy[i] = -2.0 * g_x_xy_y_xy[i] * a_exp + 4.0 * g_xyy_xy_y_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xy_y_xz[i] = -2.0 * g_x_xy_y_xz[i] * a_exp + 4.0 * g_xyy_xy_y_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xy_y_yy[i] = -2.0 * g_x_xy_y_yy[i] * a_exp + 4.0 * g_xyy_xy_y_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xy_y_yz[i] = -2.0 * g_x_xy_y_yz[i] * a_exp + 4.0 * g_xyy_xy_y_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xy_y_zz[i] = -2.0 * g_x_xy_y_zz[i] * a_exp + 4.0 * g_xyy_xy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (462-468)

    #pragma omp simd aligned(g_x_xy_z_xx, g_x_xy_z_xy, g_x_xy_z_xz, g_x_xy_z_yy, g_x_xy_z_yz, g_x_xy_z_zz, g_xy_0_0_0_y_xy_z_xx, g_xy_0_0_0_y_xy_z_xy, g_xy_0_0_0_y_xy_z_xz, g_xy_0_0_0_y_xy_z_yy, g_xy_0_0_0_y_xy_z_yz, g_xy_0_0_0_y_xy_z_zz, g_xyy_xy_z_xx, g_xyy_xy_z_xy, g_xyy_xy_z_xz, g_xyy_xy_z_yy, g_xyy_xy_z_yz, g_xyy_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_xy_z_xx[i] = -2.0 * g_x_xy_z_xx[i] * a_exp + 4.0 * g_xyy_xy_z_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xy_z_xy[i] = -2.0 * g_x_xy_z_xy[i] * a_exp + 4.0 * g_xyy_xy_z_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xy_z_xz[i] = -2.0 * g_x_xy_z_xz[i] * a_exp + 4.0 * g_xyy_xy_z_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xy_z_yy[i] = -2.0 * g_x_xy_z_yy[i] * a_exp + 4.0 * g_xyy_xy_z_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xy_z_yz[i] = -2.0 * g_x_xy_z_yz[i] * a_exp + 4.0 * g_xyy_xy_z_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xy_z_zz[i] = -2.0 * g_x_xy_z_zz[i] * a_exp + 4.0 * g_xyy_xy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (468-474)

    #pragma omp simd aligned(g_x_xz_x_xx, g_x_xz_x_xy, g_x_xz_x_xz, g_x_xz_x_yy, g_x_xz_x_yz, g_x_xz_x_zz, g_xy_0_0_0_y_xz_x_xx, g_xy_0_0_0_y_xz_x_xy, g_xy_0_0_0_y_xz_x_xz, g_xy_0_0_0_y_xz_x_yy, g_xy_0_0_0_y_xz_x_yz, g_xy_0_0_0_y_xz_x_zz, g_xyy_xz_x_xx, g_xyy_xz_x_xy, g_xyy_xz_x_xz, g_xyy_xz_x_yy, g_xyy_xz_x_yz, g_xyy_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_xz_x_xx[i] = -2.0 * g_x_xz_x_xx[i] * a_exp + 4.0 * g_xyy_xz_x_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xz_x_xy[i] = -2.0 * g_x_xz_x_xy[i] * a_exp + 4.0 * g_xyy_xz_x_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xz_x_xz[i] = -2.0 * g_x_xz_x_xz[i] * a_exp + 4.0 * g_xyy_xz_x_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xz_x_yy[i] = -2.0 * g_x_xz_x_yy[i] * a_exp + 4.0 * g_xyy_xz_x_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xz_x_yz[i] = -2.0 * g_x_xz_x_yz[i] * a_exp + 4.0 * g_xyy_xz_x_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xz_x_zz[i] = -2.0 * g_x_xz_x_zz[i] * a_exp + 4.0 * g_xyy_xz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (474-480)

    #pragma omp simd aligned(g_x_xz_y_xx, g_x_xz_y_xy, g_x_xz_y_xz, g_x_xz_y_yy, g_x_xz_y_yz, g_x_xz_y_zz, g_xy_0_0_0_y_xz_y_xx, g_xy_0_0_0_y_xz_y_xy, g_xy_0_0_0_y_xz_y_xz, g_xy_0_0_0_y_xz_y_yy, g_xy_0_0_0_y_xz_y_yz, g_xy_0_0_0_y_xz_y_zz, g_xyy_xz_y_xx, g_xyy_xz_y_xy, g_xyy_xz_y_xz, g_xyy_xz_y_yy, g_xyy_xz_y_yz, g_xyy_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_xz_y_xx[i] = -2.0 * g_x_xz_y_xx[i] * a_exp + 4.0 * g_xyy_xz_y_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xz_y_xy[i] = -2.0 * g_x_xz_y_xy[i] * a_exp + 4.0 * g_xyy_xz_y_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xz_y_xz[i] = -2.0 * g_x_xz_y_xz[i] * a_exp + 4.0 * g_xyy_xz_y_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xz_y_yy[i] = -2.0 * g_x_xz_y_yy[i] * a_exp + 4.0 * g_xyy_xz_y_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xz_y_yz[i] = -2.0 * g_x_xz_y_yz[i] * a_exp + 4.0 * g_xyy_xz_y_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xz_y_zz[i] = -2.0 * g_x_xz_y_zz[i] * a_exp + 4.0 * g_xyy_xz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (480-486)

    #pragma omp simd aligned(g_x_xz_z_xx, g_x_xz_z_xy, g_x_xz_z_xz, g_x_xz_z_yy, g_x_xz_z_yz, g_x_xz_z_zz, g_xy_0_0_0_y_xz_z_xx, g_xy_0_0_0_y_xz_z_xy, g_xy_0_0_0_y_xz_z_xz, g_xy_0_0_0_y_xz_z_yy, g_xy_0_0_0_y_xz_z_yz, g_xy_0_0_0_y_xz_z_zz, g_xyy_xz_z_xx, g_xyy_xz_z_xy, g_xyy_xz_z_xz, g_xyy_xz_z_yy, g_xyy_xz_z_yz, g_xyy_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_xz_z_xx[i] = -2.0 * g_x_xz_z_xx[i] * a_exp + 4.0 * g_xyy_xz_z_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xz_z_xy[i] = -2.0 * g_x_xz_z_xy[i] * a_exp + 4.0 * g_xyy_xz_z_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xz_z_xz[i] = -2.0 * g_x_xz_z_xz[i] * a_exp + 4.0 * g_xyy_xz_z_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xz_z_yy[i] = -2.0 * g_x_xz_z_yy[i] * a_exp + 4.0 * g_xyy_xz_z_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xz_z_yz[i] = -2.0 * g_x_xz_z_yz[i] * a_exp + 4.0 * g_xyy_xz_z_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xz_z_zz[i] = -2.0 * g_x_xz_z_zz[i] * a_exp + 4.0 * g_xyy_xz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (486-492)

    #pragma omp simd aligned(g_x_yy_x_xx, g_x_yy_x_xy, g_x_yy_x_xz, g_x_yy_x_yy, g_x_yy_x_yz, g_x_yy_x_zz, g_xy_0_0_0_y_yy_x_xx, g_xy_0_0_0_y_yy_x_xy, g_xy_0_0_0_y_yy_x_xz, g_xy_0_0_0_y_yy_x_yy, g_xy_0_0_0_y_yy_x_yz, g_xy_0_0_0_y_yy_x_zz, g_xyy_yy_x_xx, g_xyy_yy_x_xy, g_xyy_yy_x_xz, g_xyy_yy_x_yy, g_xyy_yy_x_yz, g_xyy_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_yy_x_xx[i] = -2.0 * g_x_yy_x_xx[i] * a_exp + 4.0 * g_xyy_yy_x_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yy_x_xy[i] = -2.0 * g_x_yy_x_xy[i] * a_exp + 4.0 * g_xyy_yy_x_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yy_x_xz[i] = -2.0 * g_x_yy_x_xz[i] * a_exp + 4.0 * g_xyy_yy_x_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yy_x_yy[i] = -2.0 * g_x_yy_x_yy[i] * a_exp + 4.0 * g_xyy_yy_x_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yy_x_yz[i] = -2.0 * g_x_yy_x_yz[i] * a_exp + 4.0 * g_xyy_yy_x_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yy_x_zz[i] = -2.0 * g_x_yy_x_zz[i] * a_exp + 4.0 * g_xyy_yy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (492-498)

    #pragma omp simd aligned(g_x_yy_y_xx, g_x_yy_y_xy, g_x_yy_y_xz, g_x_yy_y_yy, g_x_yy_y_yz, g_x_yy_y_zz, g_xy_0_0_0_y_yy_y_xx, g_xy_0_0_0_y_yy_y_xy, g_xy_0_0_0_y_yy_y_xz, g_xy_0_0_0_y_yy_y_yy, g_xy_0_0_0_y_yy_y_yz, g_xy_0_0_0_y_yy_y_zz, g_xyy_yy_y_xx, g_xyy_yy_y_xy, g_xyy_yy_y_xz, g_xyy_yy_y_yy, g_xyy_yy_y_yz, g_xyy_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_yy_y_xx[i] = -2.0 * g_x_yy_y_xx[i] * a_exp + 4.0 * g_xyy_yy_y_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yy_y_xy[i] = -2.0 * g_x_yy_y_xy[i] * a_exp + 4.0 * g_xyy_yy_y_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yy_y_xz[i] = -2.0 * g_x_yy_y_xz[i] * a_exp + 4.0 * g_xyy_yy_y_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yy_y_yy[i] = -2.0 * g_x_yy_y_yy[i] * a_exp + 4.0 * g_xyy_yy_y_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yy_y_yz[i] = -2.0 * g_x_yy_y_yz[i] * a_exp + 4.0 * g_xyy_yy_y_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yy_y_zz[i] = -2.0 * g_x_yy_y_zz[i] * a_exp + 4.0 * g_xyy_yy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (498-504)

    #pragma omp simd aligned(g_x_yy_z_xx, g_x_yy_z_xy, g_x_yy_z_xz, g_x_yy_z_yy, g_x_yy_z_yz, g_x_yy_z_zz, g_xy_0_0_0_y_yy_z_xx, g_xy_0_0_0_y_yy_z_xy, g_xy_0_0_0_y_yy_z_xz, g_xy_0_0_0_y_yy_z_yy, g_xy_0_0_0_y_yy_z_yz, g_xy_0_0_0_y_yy_z_zz, g_xyy_yy_z_xx, g_xyy_yy_z_xy, g_xyy_yy_z_xz, g_xyy_yy_z_yy, g_xyy_yy_z_yz, g_xyy_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_yy_z_xx[i] = -2.0 * g_x_yy_z_xx[i] * a_exp + 4.0 * g_xyy_yy_z_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yy_z_xy[i] = -2.0 * g_x_yy_z_xy[i] * a_exp + 4.0 * g_xyy_yy_z_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yy_z_xz[i] = -2.0 * g_x_yy_z_xz[i] * a_exp + 4.0 * g_xyy_yy_z_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yy_z_yy[i] = -2.0 * g_x_yy_z_yy[i] * a_exp + 4.0 * g_xyy_yy_z_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yy_z_yz[i] = -2.0 * g_x_yy_z_yz[i] * a_exp + 4.0 * g_xyy_yy_z_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yy_z_zz[i] = -2.0 * g_x_yy_z_zz[i] * a_exp + 4.0 * g_xyy_yy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (504-510)

    #pragma omp simd aligned(g_x_yz_x_xx, g_x_yz_x_xy, g_x_yz_x_xz, g_x_yz_x_yy, g_x_yz_x_yz, g_x_yz_x_zz, g_xy_0_0_0_y_yz_x_xx, g_xy_0_0_0_y_yz_x_xy, g_xy_0_0_0_y_yz_x_xz, g_xy_0_0_0_y_yz_x_yy, g_xy_0_0_0_y_yz_x_yz, g_xy_0_0_0_y_yz_x_zz, g_xyy_yz_x_xx, g_xyy_yz_x_xy, g_xyy_yz_x_xz, g_xyy_yz_x_yy, g_xyy_yz_x_yz, g_xyy_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_yz_x_xx[i] = -2.0 * g_x_yz_x_xx[i] * a_exp + 4.0 * g_xyy_yz_x_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yz_x_xy[i] = -2.0 * g_x_yz_x_xy[i] * a_exp + 4.0 * g_xyy_yz_x_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yz_x_xz[i] = -2.0 * g_x_yz_x_xz[i] * a_exp + 4.0 * g_xyy_yz_x_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yz_x_yy[i] = -2.0 * g_x_yz_x_yy[i] * a_exp + 4.0 * g_xyy_yz_x_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yz_x_yz[i] = -2.0 * g_x_yz_x_yz[i] * a_exp + 4.0 * g_xyy_yz_x_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yz_x_zz[i] = -2.0 * g_x_yz_x_zz[i] * a_exp + 4.0 * g_xyy_yz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (510-516)

    #pragma omp simd aligned(g_x_yz_y_xx, g_x_yz_y_xy, g_x_yz_y_xz, g_x_yz_y_yy, g_x_yz_y_yz, g_x_yz_y_zz, g_xy_0_0_0_y_yz_y_xx, g_xy_0_0_0_y_yz_y_xy, g_xy_0_0_0_y_yz_y_xz, g_xy_0_0_0_y_yz_y_yy, g_xy_0_0_0_y_yz_y_yz, g_xy_0_0_0_y_yz_y_zz, g_xyy_yz_y_xx, g_xyy_yz_y_xy, g_xyy_yz_y_xz, g_xyy_yz_y_yy, g_xyy_yz_y_yz, g_xyy_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_yz_y_xx[i] = -2.0 * g_x_yz_y_xx[i] * a_exp + 4.0 * g_xyy_yz_y_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yz_y_xy[i] = -2.0 * g_x_yz_y_xy[i] * a_exp + 4.0 * g_xyy_yz_y_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yz_y_xz[i] = -2.0 * g_x_yz_y_xz[i] * a_exp + 4.0 * g_xyy_yz_y_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yz_y_yy[i] = -2.0 * g_x_yz_y_yy[i] * a_exp + 4.0 * g_xyy_yz_y_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yz_y_yz[i] = -2.0 * g_x_yz_y_yz[i] * a_exp + 4.0 * g_xyy_yz_y_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yz_y_zz[i] = -2.0 * g_x_yz_y_zz[i] * a_exp + 4.0 * g_xyy_yz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (516-522)

    #pragma omp simd aligned(g_x_yz_z_xx, g_x_yz_z_xy, g_x_yz_z_xz, g_x_yz_z_yy, g_x_yz_z_yz, g_x_yz_z_zz, g_xy_0_0_0_y_yz_z_xx, g_xy_0_0_0_y_yz_z_xy, g_xy_0_0_0_y_yz_z_xz, g_xy_0_0_0_y_yz_z_yy, g_xy_0_0_0_y_yz_z_yz, g_xy_0_0_0_y_yz_z_zz, g_xyy_yz_z_xx, g_xyy_yz_z_xy, g_xyy_yz_z_xz, g_xyy_yz_z_yy, g_xyy_yz_z_yz, g_xyy_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_yz_z_xx[i] = -2.0 * g_x_yz_z_xx[i] * a_exp + 4.0 * g_xyy_yz_z_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yz_z_xy[i] = -2.0 * g_x_yz_z_xy[i] * a_exp + 4.0 * g_xyy_yz_z_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yz_z_xz[i] = -2.0 * g_x_yz_z_xz[i] * a_exp + 4.0 * g_xyy_yz_z_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yz_z_yy[i] = -2.0 * g_x_yz_z_yy[i] * a_exp + 4.0 * g_xyy_yz_z_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yz_z_yz[i] = -2.0 * g_x_yz_z_yz[i] * a_exp + 4.0 * g_xyy_yz_z_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yz_z_zz[i] = -2.0 * g_x_yz_z_zz[i] * a_exp + 4.0 * g_xyy_yz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (522-528)

    #pragma omp simd aligned(g_x_zz_x_xx, g_x_zz_x_xy, g_x_zz_x_xz, g_x_zz_x_yy, g_x_zz_x_yz, g_x_zz_x_zz, g_xy_0_0_0_y_zz_x_xx, g_xy_0_0_0_y_zz_x_xy, g_xy_0_0_0_y_zz_x_xz, g_xy_0_0_0_y_zz_x_yy, g_xy_0_0_0_y_zz_x_yz, g_xy_0_0_0_y_zz_x_zz, g_xyy_zz_x_xx, g_xyy_zz_x_xy, g_xyy_zz_x_xz, g_xyy_zz_x_yy, g_xyy_zz_x_yz, g_xyy_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_zz_x_xx[i] = -2.0 * g_x_zz_x_xx[i] * a_exp + 4.0 * g_xyy_zz_x_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_zz_x_xy[i] = -2.0 * g_x_zz_x_xy[i] * a_exp + 4.0 * g_xyy_zz_x_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_zz_x_xz[i] = -2.0 * g_x_zz_x_xz[i] * a_exp + 4.0 * g_xyy_zz_x_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_zz_x_yy[i] = -2.0 * g_x_zz_x_yy[i] * a_exp + 4.0 * g_xyy_zz_x_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_zz_x_yz[i] = -2.0 * g_x_zz_x_yz[i] * a_exp + 4.0 * g_xyy_zz_x_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_zz_x_zz[i] = -2.0 * g_x_zz_x_zz[i] * a_exp + 4.0 * g_xyy_zz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (528-534)

    #pragma omp simd aligned(g_x_zz_y_xx, g_x_zz_y_xy, g_x_zz_y_xz, g_x_zz_y_yy, g_x_zz_y_yz, g_x_zz_y_zz, g_xy_0_0_0_y_zz_y_xx, g_xy_0_0_0_y_zz_y_xy, g_xy_0_0_0_y_zz_y_xz, g_xy_0_0_0_y_zz_y_yy, g_xy_0_0_0_y_zz_y_yz, g_xy_0_0_0_y_zz_y_zz, g_xyy_zz_y_xx, g_xyy_zz_y_xy, g_xyy_zz_y_xz, g_xyy_zz_y_yy, g_xyy_zz_y_yz, g_xyy_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_zz_y_xx[i] = -2.0 * g_x_zz_y_xx[i] * a_exp + 4.0 * g_xyy_zz_y_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_zz_y_xy[i] = -2.0 * g_x_zz_y_xy[i] * a_exp + 4.0 * g_xyy_zz_y_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_zz_y_xz[i] = -2.0 * g_x_zz_y_xz[i] * a_exp + 4.0 * g_xyy_zz_y_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_zz_y_yy[i] = -2.0 * g_x_zz_y_yy[i] * a_exp + 4.0 * g_xyy_zz_y_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_zz_y_yz[i] = -2.0 * g_x_zz_y_yz[i] * a_exp + 4.0 * g_xyy_zz_y_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_zz_y_zz[i] = -2.0 * g_x_zz_y_zz[i] * a_exp + 4.0 * g_xyy_zz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (534-540)

    #pragma omp simd aligned(g_x_zz_z_xx, g_x_zz_z_xy, g_x_zz_z_xz, g_x_zz_z_yy, g_x_zz_z_yz, g_x_zz_z_zz, g_xy_0_0_0_y_zz_z_xx, g_xy_0_0_0_y_zz_z_xy, g_xy_0_0_0_y_zz_z_xz, g_xy_0_0_0_y_zz_z_yy, g_xy_0_0_0_y_zz_z_yz, g_xy_0_0_0_y_zz_z_zz, g_xyy_zz_z_xx, g_xyy_zz_z_xy, g_xyy_zz_z_xz, g_xyy_zz_z_yy, g_xyy_zz_z_yz, g_xyy_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_zz_z_xx[i] = -2.0 * g_x_zz_z_xx[i] * a_exp + 4.0 * g_xyy_zz_z_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_y_zz_z_xy[i] = -2.0 * g_x_zz_z_xy[i] * a_exp + 4.0 * g_xyy_zz_z_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_zz_z_xz[i] = -2.0 * g_x_zz_z_xz[i] * a_exp + 4.0 * g_xyy_zz_z_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_zz_z_yy[i] = -2.0 * g_x_zz_z_yy[i] * a_exp + 4.0 * g_xyy_zz_z_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_y_zz_z_yz[i] = -2.0 * g_x_zz_z_yz[i] * a_exp + 4.0 * g_xyy_zz_z_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_y_zz_z_zz[i] = -2.0 * g_x_zz_z_zz[i] * a_exp + 4.0 * g_xyy_zz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (540-546)

    #pragma omp simd aligned(g_xy_0_0_0_z_xx_x_xx, g_xy_0_0_0_z_xx_x_xy, g_xy_0_0_0_z_xx_x_xz, g_xy_0_0_0_z_xx_x_yy, g_xy_0_0_0_z_xx_x_yz, g_xy_0_0_0_z_xx_x_zz, g_xyz_xx_x_xx, g_xyz_xx_x_xy, g_xyz_xx_x_xz, g_xyz_xx_x_yy, g_xyz_xx_x_yz, g_xyz_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_xx_x_xx[i] = 4.0 * g_xyz_xx_x_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xx_x_xy[i] = 4.0 * g_xyz_xx_x_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xx_x_xz[i] = 4.0 * g_xyz_xx_x_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xx_x_yy[i] = 4.0 * g_xyz_xx_x_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xx_x_yz[i] = 4.0 * g_xyz_xx_x_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xx_x_zz[i] = 4.0 * g_xyz_xx_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (546-552)

    #pragma omp simd aligned(g_xy_0_0_0_z_xx_y_xx, g_xy_0_0_0_z_xx_y_xy, g_xy_0_0_0_z_xx_y_xz, g_xy_0_0_0_z_xx_y_yy, g_xy_0_0_0_z_xx_y_yz, g_xy_0_0_0_z_xx_y_zz, g_xyz_xx_y_xx, g_xyz_xx_y_xy, g_xyz_xx_y_xz, g_xyz_xx_y_yy, g_xyz_xx_y_yz, g_xyz_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_xx_y_xx[i] = 4.0 * g_xyz_xx_y_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xx_y_xy[i] = 4.0 * g_xyz_xx_y_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xx_y_xz[i] = 4.0 * g_xyz_xx_y_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xx_y_yy[i] = 4.0 * g_xyz_xx_y_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xx_y_yz[i] = 4.0 * g_xyz_xx_y_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xx_y_zz[i] = 4.0 * g_xyz_xx_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (552-558)

    #pragma omp simd aligned(g_xy_0_0_0_z_xx_z_xx, g_xy_0_0_0_z_xx_z_xy, g_xy_0_0_0_z_xx_z_xz, g_xy_0_0_0_z_xx_z_yy, g_xy_0_0_0_z_xx_z_yz, g_xy_0_0_0_z_xx_z_zz, g_xyz_xx_z_xx, g_xyz_xx_z_xy, g_xyz_xx_z_xz, g_xyz_xx_z_yy, g_xyz_xx_z_yz, g_xyz_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_xx_z_xx[i] = 4.0 * g_xyz_xx_z_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xx_z_xy[i] = 4.0 * g_xyz_xx_z_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xx_z_xz[i] = 4.0 * g_xyz_xx_z_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xx_z_yy[i] = 4.0 * g_xyz_xx_z_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xx_z_yz[i] = 4.0 * g_xyz_xx_z_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xx_z_zz[i] = 4.0 * g_xyz_xx_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (558-564)

    #pragma omp simd aligned(g_xy_0_0_0_z_xy_x_xx, g_xy_0_0_0_z_xy_x_xy, g_xy_0_0_0_z_xy_x_xz, g_xy_0_0_0_z_xy_x_yy, g_xy_0_0_0_z_xy_x_yz, g_xy_0_0_0_z_xy_x_zz, g_xyz_xy_x_xx, g_xyz_xy_x_xy, g_xyz_xy_x_xz, g_xyz_xy_x_yy, g_xyz_xy_x_yz, g_xyz_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_xy_x_xx[i] = 4.0 * g_xyz_xy_x_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xy_x_xy[i] = 4.0 * g_xyz_xy_x_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xy_x_xz[i] = 4.0 * g_xyz_xy_x_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xy_x_yy[i] = 4.0 * g_xyz_xy_x_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xy_x_yz[i] = 4.0 * g_xyz_xy_x_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xy_x_zz[i] = 4.0 * g_xyz_xy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (564-570)

    #pragma omp simd aligned(g_xy_0_0_0_z_xy_y_xx, g_xy_0_0_0_z_xy_y_xy, g_xy_0_0_0_z_xy_y_xz, g_xy_0_0_0_z_xy_y_yy, g_xy_0_0_0_z_xy_y_yz, g_xy_0_0_0_z_xy_y_zz, g_xyz_xy_y_xx, g_xyz_xy_y_xy, g_xyz_xy_y_xz, g_xyz_xy_y_yy, g_xyz_xy_y_yz, g_xyz_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_xy_y_xx[i] = 4.0 * g_xyz_xy_y_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xy_y_xy[i] = 4.0 * g_xyz_xy_y_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xy_y_xz[i] = 4.0 * g_xyz_xy_y_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xy_y_yy[i] = 4.0 * g_xyz_xy_y_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xy_y_yz[i] = 4.0 * g_xyz_xy_y_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xy_y_zz[i] = 4.0 * g_xyz_xy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (570-576)

    #pragma omp simd aligned(g_xy_0_0_0_z_xy_z_xx, g_xy_0_0_0_z_xy_z_xy, g_xy_0_0_0_z_xy_z_xz, g_xy_0_0_0_z_xy_z_yy, g_xy_0_0_0_z_xy_z_yz, g_xy_0_0_0_z_xy_z_zz, g_xyz_xy_z_xx, g_xyz_xy_z_xy, g_xyz_xy_z_xz, g_xyz_xy_z_yy, g_xyz_xy_z_yz, g_xyz_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_xy_z_xx[i] = 4.0 * g_xyz_xy_z_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xy_z_xy[i] = 4.0 * g_xyz_xy_z_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xy_z_xz[i] = 4.0 * g_xyz_xy_z_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xy_z_yy[i] = 4.0 * g_xyz_xy_z_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xy_z_yz[i] = 4.0 * g_xyz_xy_z_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xy_z_zz[i] = 4.0 * g_xyz_xy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (576-582)

    #pragma omp simd aligned(g_xy_0_0_0_z_xz_x_xx, g_xy_0_0_0_z_xz_x_xy, g_xy_0_0_0_z_xz_x_xz, g_xy_0_0_0_z_xz_x_yy, g_xy_0_0_0_z_xz_x_yz, g_xy_0_0_0_z_xz_x_zz, g_xyz_xz_x_xx, g_xyz_xz_x_xy, g_xyz_xz_x_xz, g_xyz_xz_x_yy, g_xyz_xz_x_yz, g_xyz_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_xz_x_xx[i] = 4.0 * g_xyz_xz_x_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xz_x_xy[i] = 4.0 * g_xyz_xz_x_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xz_x_xz[i] = 4.0 * g_xyz_xz_x_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xz_x_yy[i] = 4.0 * g_xyz_xz_x_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xz_x_yz[i] = 4.0 * g_xyz_xz_x_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xz_x_zz[i] = 4.0 * g_xyz_xz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (582-588)

    #pragma omp simd aligned(g_xy_0_0_0_z_xz_y_xx, g_xy_0_0_0_z_xz_y_xy, g_xy_0_0_0_z_xz_y_xz, g_xy_0_0_0_z_xz_y_yy, g_xy_0_0_0_z_xz_y_yz, g_xy_0_0_0_z_xz_y_zz, g_xyz_xz_y_xx, g_xyz_xz_y_xy, g_xyz_xz_y_xz, g_xyz_xz_y_yy, g_xyz_xz_y_yz, g_xyz_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_xz_y_xx[i] = 4.0 * g_xyz_xz_y_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xz_y_xy[i] = 4.0 * g_xyz_xz_y_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xz_y_xz[i] = 4.0 * g_xyz_xz_y_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xz_y_yy[i] = 4.0 * g_xyz_xz_y_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xz_y_yz[i] = 4.0 * g_xyz_xz_y_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xz_y_zz[i] = 4.0 * g_xyz_xz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (588-594)

    #pragma omp simd aligned(g_xy_0_0_0_z_xz_z_xx, g_xy_0_0_0_z_xz_z_xy, g_xy_0_0_0_z_xz_z_xz, g_xy_0_0_0_z_xz_z_yy, g_xy_0_0_0_z_xz_z_yz, g_xy_0_0_0_z_xz_z_zz, g_xyz_xz_z_xx, g_xyz_xz_z_xy, g_xyz_xz_z_xz, g_xyz_xz_z_yy, g_xyz_xz_z_yz, g_xyz_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_xz_z_xx[i] = 4.0 * g_xyz_xz_z_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xz_z_xy[i] = 4.0 * g_xyz_xz_z_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xz_z_xz[i] = 4.0 * g_xyz_xz_z_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xz_z_yy[i] = 4.0 * g_xyz_xz_z_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xz_z_yz[i] = 4.0 * g_xyz_xz_z_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xz_z_zz[i] = 4.0 * g_xyz_xz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (594-600)

    #pragma omp simd aligned(g_xy_0_0_0_z_yy_x_xx, g_xy_0_0_0_z_yy_x_xy, g_xy_0_0_0_z_yy_x_xz, g_xy_0_0_0_z_yy_x_yy, g_xy_0_0_0_z_yy_x_yz, g_xy_0_0_0_z_yy_x_zz, g_xyz_yy_x_xx, g_xyz_yy_x_xy, g_xyz_yy_x_xz, g_xyz_yy_x_yy, g_xyz_yy_x_yz, g_xyz_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_yy_x_xx[i] = 4.0 * g_xyz_yy_x_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yy_x_xy[i] = 4.0 * g_xyz_yy_x_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yy_x_xz[i] = 4.0 * g_xyz_yy_x_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yy_x_yy[i] = 4.0 * g_xyz_yy_x_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yy_x_yz[i] = 4.0 * g_xyz_yy_x_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yy_x_zz[i] = 4.0 * g_xyz_yy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (600-606)

    #pragma omp simd aligned(g_xy_0_0_0_z_yy_y_xx, g_xy_0_0_0_z_yy_y_xy, g_xy_0_0_0_z_yy_y_xz, g_xy_0_0_0_z_yy_y_yy, g_xy_0_0_0_z_yy_y_yz, g_xy_0_0_0_z_yy_y_zz, g_xyz_yy_y_xx, g_xyz_yy_y_xy, g_xyz_yy_y_xz, g_xyz_yy_y_yy, g_xyz_yy_y_yz, g_xyz_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_yy_y_xx[i] = 4.0 * g_xyz_yy_y_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yy_y_xy[i] = 4.0 * g_xyz_yy_y_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yy_y_xz[i] = 4.0 * g_xyz_yy_y_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yy_y_yy[i] = 4.0 * g_xyz_yy_y_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yy_y_yz[i] = 4.0 * g_xyz_yy_y_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yy_y_zz[i] = 4.0 * g_xyz_yy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (606-612)

    #pragma omp simd aligned(g_xy_0_0_0_z_yy_z_xx, g_xy_0_0_0_z_yy_z_xy, g_xy_0_0_0_z_yy_z_xz, g_xy_0_0_0_z_yy_z_yy, g_xy_0_0_0_z_yy_z_yz, g_xy_0_0_0_z_yy_z_zz, g_xyz_yy_z_xx, g_xyz_yy_z_xy, g_xyz_yy_z_xz, g_xyz_yy_z_yy, g_xyz_yy_z_yz, g_xyz_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_yy_z_xx[i] = 4.0 * g_xyz_yy_z_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yy_z_xy[i] = 4.0 * g_xyz_yy_z_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yy_z_xz[i] = 4.0 * g_xyz_yy_z_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yy_z_yy[i] = 4.0 * g_xyz_yy_z_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yy_z_yz[i] = 4.0 * g_xyz_yy_z_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yy_z_zz[i] = 4.0 * g_xyz_yy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (612-618)

    #pragma omp simd aligned(g_xy_0_0_0_z_yz_x_xx, g_xy_0_0_0_z_yz_x_xy, g_xy_0_0_0_z_yz_x_xz, g_xy_0_0_0_z_yz_x_yy, g_xy_0_0_0_z_yz_x_yz, g_xy_0_0_0_z_yz_x_zz, g_xyz_yz_x_xx, g_xyz_yz_x_xy, g_xyz_yz_x_xz, g_xyz_yz_x_yy, g_xyz_yz_x_yz, g_xyz_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_yz_x_xx[i] = 4.0 * g_xyz_yz_x_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yz_x_xy[i] = 4.0 * g_xyz_yz_x_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yz_x_xz[i] = 4.0 * g_xyz_yz_x_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yz_x_yy[i] = 4.0 * g_xyz_yz_x_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yz_x_yz[i] = 4.0 * g_xyz_yz_x_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yz_x_zz[i] = 4.0 * g_xyz_yz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (618-624)

    #pragma omp simd aligned(g_xy_0_0_0_z_yz_y_xx, g_xy_0_0_0_z_yz_y_xy, g_xy_0_0_0_z_yz_y_xz, g_xy_0_0_0_z_yz_y_yy, g_xy_0_0_0_z_yz_y_yz, g_xy_0_0_0_z_yz_y_zz, g_xyz_yz_y_xx, g_xyz_yz_y_xy, g_xyz_yz_y_xz, g_xyz_yz_y_yy, g_xyz_yz_y_yz, g_xyz_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_yz_y_xx[i] = 4.0 * g_xyz_yz_y_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yz_y_xy[i] = 4.0 * g_xyz_yz_y_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yz_y_xz[i] = 4.0 * g_xyz_yz_y_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yz_y_yy[i] = 4.0 * g_xyz_yz_y_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yz_y_yz[i] = 4.0 * g_xyz_yz_y_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yz_y_zz[i] = 4.0 * g_xyz_yz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (624-630)

    #pragma omp simd aligned(g_xy_0_0_0_z_yz_z_xx, g_xy_0_0_0_z_yz_z_xy, g_xy_0_0_0_z_yz_z_xz, g_xy_0_0_0_z_yz_z_yy, g_xy_0_0_0_z_yz_z_yz, g_xy_0_0_0_z_yz_z_zz, g_xyz_yz_z_xx, g_xyz_yz_z_xy, g_xyz_yz_z_xz, g_xyz_yz_z_yy, g_xyz_yz_z_yz, g_xyz_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_yz_z_xx[i] = 4.0 * g_xyz_yz_z_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yz_z_xy[i] = 4.0 * g_xyz_yz_z_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yz_z_xz[i] = 4.0 * g_xyz_yz_z_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yz_z_yy[i] = 4.0 * g_xyz_yz_z_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yz_z_yz[i] = 4.0 * g_xyz_yz_z_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yz_z_zz[i] = 4.0 * g_xyz_yz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (630-636)

    #pragma omp simd aligned(g_xy_0_0_0_z_zz_x_xx, g_xy_0_0_0_z_zz_x_xy, g_xy_0_0_0_z_zz_x_xz, g_xy_0_0_0_z_zz_x_yy, g_xy_0_0_0_z_zz_x_yz, g_xy_0_0_0_z_zz_x_zz, g_xyz_zz_x_xx, g_xyz_zz_x_xy, g_xyz_zz_x_xz, g_xyz_zz_x_yy, g_xyz_zz_x_yz, g_xyz_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_zz_x_xx[i] = 4.0 * g_xyz_zz_x_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_zz_x_xy[i] = 4.0 * g_xyz_zz_x_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_zz_x_xz[i] = 4.0 * g_xyz_zz_x_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_zz_x_yy[i] = 4.0 * g_xyz_zz_x_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_zz_x_yz[i] = 4.0 * g_xyz_zz_x_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_zz_x_zz[i] = 4.0 * g_xyz_zz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (636-642)

    #pragma omp simd aligned(g_xy_0_0_0_z_zz_y_xx, g_xy_0_0_0_z_zz_y_xy, g_xy_0_0_0_z_zz_y_xz, g_xy_0_0_0_z_zz_y_yy, g_xy_0_0_0_z_zz_y_yz, g_xy_0_0_0_z_zz_y_zz, g_xyz_zz_y_xx, g_xyz_zz_y_xy, g_xyz_zz_y_xz, g_xyz_zz_y_yy, g_xyz_zz_y_yz, g_xyz_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_zz_y_xx[i] = 4.0 * g_xyz_zz_y_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_zz_y_xy[i] = 4.0 * g_xyz_zz_y_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_zz_y_xz[i] = 4.0 * g_xyz_zz_y_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_zz_y_yy[i] = 4.0 * g_xyz_zz_y_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_zz_y_yz[i] = 4.0 * g_xyz_zz_y_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_zz_y_zz[i] = 4.0 * g_xyz_zz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (642-648)

    #pragma omp simd aligned(g_xy_0_0_0_z_zz_z_xx, g_xy_0_0_0_z_zz_z_xy, g_xy_0_0_0_z_zz_z_xz, g_xy_0_0_0_z_zz_z_yy, g_xy_0_0_0_z_zz_z_yz, g_xy_0_0_0_z_zz_z_zz, g_xyz_zz_z_xx, g_xyz_zz_z_xy, g_xyz_zz_z_xz, g_xyz_zz_z_yy, g_xyz_zz_z_yz, g_xyz_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_zz_z_xx[i] = 4.0 * g_xyz_zz_z_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_z_zz_z_xy[i] = 4.0 * g_xyz_zz_z_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_zz_z_xz[i] = 4.0 * g_xyz_zz_z_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_zz_z_yy[i] = 4.0 * g_xyz_zz_z_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_z_zz_z_yz[i] = 4.0 * g_xyz_zz_z_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_z_zz_z_zz[i] = 4.0 * g_xyz_zz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (648-654)

    #pragma omp simd aligned(g_xxz_xx_x_xx, g_xxz_xx_x_xy, g_xxz_xx_x_xz, g_xxz_xx_x_yy, g_xxz_xx_x_yz, g_xxz_xx_x_zz, g_xz_0_0_0_x_xx_x_xx, g_xz_0_0_0_x_xx_x_xy, g_xz_0_0_0_x_xx_x_xz, g_xz_0_0_0_x_xx_x_yy, g_xz_0_0_0_x_xx_x_yz, g_xz_0_0_0_x_xx_x_zz, g_z_xx_x_xx, g_z_xx_x_xy, g_z_xx_x_xz, g_z_xx_x_yy, g_z_xx_x_yz, g_z_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_xx_x_xx[i] = -2.0 * g_z_xx_x_xx[i] * a_exp + 4.0 * g_xxz_xx_x_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xx_x_xy[i] = -2.0 * g_z_xx_x_xy[i] * a_exp + 4.0 * g_xxz_xx_x_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xx_x_xz[i] = -2.0 * g_z_xx_x_xz[i] * a_exp + 4.0 * g_xxz_xx_x_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xx_x_yy[i] = -2.0 * g_z_xx_x_yy[i] * a_exp + 4.0 * g_xxz_xx_x_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xx_x_yz[i] = -2.0 * g_z_xx_x_yz[i] * a_exp + 4.0 * g_xxz_xx_x_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xx_x_zz[i] = -2.0 * g_z_xx_x_zz[i] * a_exp + 4.0 * g_xxz_xx_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (654-660)

    #pragma omp simd aligned(g_xxz_xx_y_xx, g_xxz_xx_y_xy, g_xxz_xx_y_xz, g_xxz_xx_y_yy, g_xxz_xx_y_yz, g_xxz_xx_y_zz, g_xz_0_0_0_x_xx_y_xx, g_xz_0_0_0_x_xx_y_xy, g_xz_0_0_0_x_xx_y_xz, g_xz_0_0_0_x_xx_y_yy, g_xz_0_0_0_x_xx_y_yz, g_xz_0_0_0_x_xx_y_zz, g_z_xx_y_xx, g_z_xx_y_xy, g_z_xx_y_xz, g_z_xx_y_yy, g_z_xx_y_yz, g_z_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_xx_y_xx[i] = -2.0 * g_z_xx_y_xx[i] * a_exp + 4.0 * g_xxz_xx_y_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xx_y_xy[i] = -2.0 * g_z_xx_y_xy[i] * a_exp + 4.0 * g_xxz_xx_y_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xx_y_xz[i] = -2.0 * g_z_xx_y_xz[i] * a_exp + 4.0 * g_xxz_xx_y_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xx_y_yy[i] = -2.0 * g_z_xx_y_yy[i] * a_exp + 4.0 * g_xxz_xx_y_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xx_y_yz[i] = -2.0 * g_z_xx_y_yz[i] * a_exp + 4.0 * g_xxz_xx_y_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xx_y_zz[i] = -2.0 * g_z_xx_y_zz[i] * a_exp + 4.0 * g_xxz_xx_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (660-666)

    #pragma omp simd aligned(g_xxz_xx_z_xx, g_xxz_xx_z_xy, g_xxz_xx_z_xz, g_xxz_xx_z_yy, g_xxz_xx_z_yz, g_xxz_xx_z_zz, g_xz_0_0_0_x_xx_z_xx, g_xz_0_0_0_x_xx_z_xy, g_xz_0_0_0_x_xx_z_xz, g_xz_0_0_0_x_xx_z_yy, g_xz_0_0_0_x_xx_z_yz, g_xz_0_0_0_x_xx_z_zz, g_z_xx_z_xx, g_z_xx_z_xy, g_z_xx_z_xz, g_z_xx_z_yy, g_z_xx_z_yz, g_z_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_xx_z_xx[i] = -2.0 * g_z_xx_z_xx[i] * a_exp + 4.0 * g_xxz_xx_z_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xx_z_xy[i] = -2.0 * g_z_xx_z_xy[i] * a_exp + 4.0 * g_xxz_xx_z_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xx_z_xz[i] = -2.0 * g_z_xx_z_xz[i] * a_exp + 4.0 * g_xxz_xx_z_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xx_z_yy[i] = -2.0 * g_z_xx_z_yy[i] * a_exp + 4.0 * g_xxz_xx_z_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xx_z_yz[i] = -2.0 * g_z_xx_z_yz[i] * a_exp + 4.0 * g_xxz_xx_z_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xx_z_zz[i] = -2.0 * g_z_xx_z_zz[i] * a_exp + 4.0 * g_xxz_xx_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (666-672)

    #pragma omp simd aligned(g_xxz_xy_x_xx, g_xxz_xy_x_xy, g_xxz_xy_x_xz, g_xxz_xy_x_yy, g_xxz_xy_x_yz, g_xxz_xy_x_zz, g_xz_0_0_0_x_xy_x_xx, g_xz_0_0_0_x_xy_x_xy, g_xz_0_0_0_x_xy_x_xz, g_xz_0_0_0_x_xy_x_yy, g_xz_0_0_0_x_xy_x_yz, g_xz_0_0_0_x_xy_x_zz, g_z_xy_x_xx, g_z_xy_x_xy, g_z_xy_x_xz, g_z_xy_x_yy, g_z_xy_x_yz, g_z_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_xy_x_xx[i] = -2.0 * g_z_xy_x_xx[i] * a_exp + 4.0 * g_xxz_xy_x_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xy_x_xy[i] = -2.0 * g_z_xy_x_xy[i] * a_exp + 4.0 * g_xxz_xy_x_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xy_x_xz[i] = -2.0 * g_z_xy_x_xz[i] * a_exp + 4.0 * g_xxz_xy_x_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xy_x_yy[i] = -2.0 * g_z_xy_x_yy[i] * a_exp + 4.0 * g_xxz_xy_x_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xy_x_yz[i] = -2.0 * g_z_xy_x_yz[i] * a_exp + 4.0 * g_xxz_xy_x_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xy_x_zz[i] = -2.0 * g_z_xy_x_zz[i] * a_exp + 4.0 * g_xxz_xy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (672-678)

    #pragma omp simd aligned(g_xxz_xy_y_xx, g_xxz_xy_y_xy, g_xxz_xy_y_xz, g_xxz_xy_y_yy, g_xxz_xy_y_yz, g_xxz_xy_y_zz, g_xz_0_0_0_x_xy_y_xx, g_xz_0_0_0_x_xy_y_xy, g_xz_0_0_0_x_xy_y_xz, g_xz_0_0_0_x_xy_y_yy, g_xz_0_0_0_x_xy_y_yz, g_xz_0_0_0_x_xy_y_zz, g_z_xy_y_xx, g_z_xy_y_xy, g_z_xy_y_xz, g_z_xy_y_yy, g_z_xy_y_yz, g_z_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_xy_y_xx[i] = -2.0 * g_z_xy_y_xx[i] * a_exp + 4.0 * g_xxz_xy_y_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xy_y_xy[i] = -2.0 * g_z_xy_y_xy[i] * a_exp + 4.0 * g_xxz_xy_y_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xy_y_xz[i] = -2.0 * g_z_xy_y_xz[i] * a_exp + 4.0 * g_xxz_xy_y_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xy_y_yy[i] = -2.0 * g_z_xy_y_yy[i] * a_exp + 4.0 * g_xxz_xy_y_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xy_y_yz[i] = -2.0 * g_z_xy_y_yz[i] * a_exp + 4.0 * g_xxz_xy_y_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xy_y_zz[i] = -2.0 * g_z_xy_y_zz[i] * a_exp + 4.0 * g_xxz_xy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (678-684)

    #pragma omp simd aligned(g_xxz_xy_z_xx, g_xxz_xy_z_xy, g_xxz_xy_z_xz, g_xxz_xy_z_yy, g_xxz_xy_z_yz, g_xxz_xy_z_zz, g_xz_0_0_0_x_xy_z_xx, g_xz_0_0_0_x_xy_z_xy, g_xz_0_0_0_x_xy_z_xz, g_xz_0_0_0_x_xy_z_yy, g_xz_0_0_0_x_xy_z_yz, g_xz_0_0_0_x_xy_z_zz, g_z_xy_z_xx, g_z_xy_z_xy, g_z_xy_z_xz, g_z_xy_z_yy, g_z_xy_z_yz, g_z_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_xy_z_xx[i] = -2.0 * g_z_xy_z_xx[i] * a_exp + 4.0 * g_xxz_xy_z_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xy_z_xy[i] = -2.0 * g_z_xy_z_xy[i] * a_exp + 4.0 * g_xxz_xy_z_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xy_z_xz[i] = -2.0 * g_z_xy_z_xz[i] * a_exp + 4.0 * g_xxz_xy_z_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xy_z_yy[i] = -2.0 * g_z_xy_z_yy[i] * a_exp + 4.0 * g_xxz_xy_z_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xy_z_yz[i] = -2.0 * g_z_xy_z_yz[i] * a_exp + 4.0 * g_xxz_xy_z_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xy_z_zz[i] = -2.0 * g_z_xy_z_zz[i] * a_exp + 4.0 * g_xxz_xy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (684-690)

    #pragma omp simd aligned(g_xxz_xz_x_xx, g_xxz_xz_x_xy, g_xxz_xz_x_xz, g_xxz_xz_x_yy, g_xxz_xz_x_yz, g_xxz_xz_x_zz, g_xz_0_0_0_x_xz_x_xx, g_xz_0_0_0_x_xz_x_xy, g_xz_0_0_0_x_xz_x_xz, g_xz_0_0_0_x_xz_x_yy, g_xz_0_0_0_x_xz_x_yz, g_xz_0_0_0_x_xz_x_zz, g_z_xz_x_xx, g_z_xz_x_xy, g_z_xz_x_xz, g_z_xz_x_yy, g_z_xz_x_yz, g_z_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_xz_x_xx[i] = -2.0 * g_z_xz_x_xx[i] * a_exp + 4.0 * g_xxz_xz_x_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xz_x_xy[i] = -2.0 * g_z_xz_x_xy[i] * a_exp + 4.0 * g_xxz_xz_x_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xz_x_xz[i] = -2.0 * g_z_xz_x_xz[i] * a_exp + 4.0 * g_xxz_xz_x_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xz_x_yy[i] = -2.0 * g_z_xz_x_yy[i] * a_exp + 4.0 * g_xxz_xz_x_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xz_x_yz[i] = -2.0 * g_z_xz_x_yz[i] * a_exp + 4.0 * g_xxz_xz_x_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xz_x_zz[i] = -2.0 * g_z_xz_x_zz[i] * a_exp + 4.0 * g_xxz_xz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (690-696)

    #pragma omp simd aligned(g_xxz_xz_y_xx, g_xxz_xz_y_xy, g_xxz_xz_y_xz, g_xxz_xz_y_yy, g_xxz_xz_y_yz, g_xxz_xz_y_zz, g_xz_0_0_0_x_xz_y_xx, g_xz_0_0_0_x_xz_y_xy, g_xz_0_0_0_x_xz_y_xz, g_xz_0_0_0_x_xz_y_yy, g_xz_0_0_0_x_xz_y_yz, g_xz_0_0_0_x_xz_y_zz, g_z_xz_y_xx, g_z_xz_y_xy, g_z_xz_y_xz, g_z_xz_y_yy, g_z_xz_y_yz, g_z_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_xz_y_xx[i] = -2.0 * g_z_xz_y_xx[i] * a_exp + 4.0 * g_xxz_xz_y_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xz_y_xy[i] = -2.0 * g_z_xz_y_xy[i] * a_exp + 4.0 * g_xxz_xz_y_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xz_y_xz[i] = -2.0 * g_z_xz_y_xz[i] * a_exp + 4.0 * g_xxz_xz_y_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xz_y_yy[i] = -2.0 * g_z_xz_y_yy[i] * a_exp + 4.0 * g_xxz_xz_y_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xz_y_yz[i] = -2.0 * g_z_xz_y_yz[i] * a_exp + 4.0 * g_xxz_xz_y_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xz_y_zz[i] = -2.0 * g_z_xz_y_zz[i] * a_exp + 4.0 * g_xxz_xz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (696-702)

    #pragma omp simd aligned(g_xxz_xz_z_xx, g_xxz_xz_z_xy, g_xxz_xz_z_xz, g_xxz_xz_z_yy, g_xxz_xz_z_yz, g_xxz_xz_z_zz, g_xz_0_0_0_x_xz_z_xx, g_xz_0_0_0_x_xz_z_xy, g_xz_0_0_0_x_xz_z_xz, g_xz_0_0_0_x_xz_z_yy, g_xz_0_0_0_x_xz_z_yz, g_xz_0_0_0_x_xz_z_zz, g_z_xz_z_xx, g_z_xz_z_xy, g_z_xz_z_xz, g_z_xz_z_yy, g_z_xz_z_yz, g_z_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_xz_z_xx[i] = -2.0 * g_z_xz_z_xx[i] * a_exp + 4.0 * g_xxz_xz_z_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xz_z_xy[i] = -2.0 * g_z_xz_z_xy[i] * a_exp + 4.0 * g_xxz_xz_z_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xz_z_xz[i] = -2.0 * g_z_xz_z_xz[i] * a_exp + 4.0 * g_xxz_xz_z_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xz_z_yy[i] = -2.0 * g_z_xz_z_yy[i] * a_exp + 4.0 * g_xxz_xz_z_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xz_z_yz[i] = -2.0 * g_z_xz_z_yz[i] * a_exp + 4.0 * g_xxz_xz_z_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xz_z_zz[i] = -2.0 * g_z_xz_z_zz[i] * a_exp + 4.0 * g_xxz_xz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (702-708)

    #pragma omp simd aligned(g_xxz_yy_x_xx, g_xxz_yy_x_xy, g_xxz_yy_x_xz, g_xxz_yy_x_yy, g_xxz_yy_x_yz, g_xxz_yy_x_zz, g_xz_0_0_0_x_yy_x_xx, g_xz_0_0_0_x_yy_x_xy, g_xz_0_0_0_x_yy_x_xz, g_xz_0_0_0_x_yy_x_yy, g_xz_0_0_0_x_yy_x_yz, g_xz_0_0_0_x_yy_x_zz, g_z_yy_x_xx, g_z_yy_x_xy, g_z_yy_x_xz, g_z_yy_x_yy, g_z_yy_x_yz, g_z_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_yy_x_xx[i] = -2.0 * g_z_yy_x_xx[i] * a_exp + 4.0 * g_xxz_yy_x_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yy_x_xy[i] = -2.0 * g_z_yy_x_xy[i] * a_exp + 4.0 * g_xxz_yy_x_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yy_x_xz[i] = -2.0 * g_z_yy_x_xz[i] * a_exp + 4.0 * g_xxz_yy_x_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yy_x_yy[i] = -2.0 * g_z_yy_x_yy[i] * a_exp + 4.0 * g_xxz_yy_x_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yy_x_yz[i] = -2.0 * g_z_yy_x_yz[i] * a_exp + 4.0 * g_xxz_yy_x_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yy_x_zz[i] = -2.0 * g_z_yy_x_zz[i] * a_exp + 4.0 * g_xxz_yy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (708-714)

    #pragma omp simd aligned(g_xxz_yy_y_xx, g_xxz_yy_y_xy, g_xxz_yy_y_xz, g_xxz_yy_y_yy, g_xxz_yy_y_yz, g_xxz_yy_y_zz, g_xz_0_0_0_x_yy_y_xx, g_xz_0_0_0_x_yy_y_xy, g_xz_0_0_0_x_yy_y_xz, g_xz_0_0_0_x_yy_y_yy, g_xz_0_0_0_x_yy_y_yz, g_xz_0_0_0_x_yy_y_zz, g_z_yy_y_xx, g_z_yy_y_xy, g_z_yy_y_xz, g_z_yy_y_yy, g_z_yy_y_yz, g_z_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_yy_y_xx[i] = -2.0 * g_z_yy_y_xx[i] * a_exp + 4.0 * g_xxz_yy_y_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yy_y_xy[i] = -2.0 * g_z_yy_y_xy[i] * a_exp + 4.0 * g_xxz_yy_y_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yy_y_xz[i] = -2.0 * g_z_yy_y_xz[i] * a_exp + 4.0 * g_xxz_yy_y_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yy_y_yy[i] = -2.0 * g_z_yy_y_yy[i] * a_exp + 4.0 * g_xxz_yy_y_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yy_y_yz[i] = -2.0 * g_z_yy_y_yz[i] * a_exp + 4.0 * g_xxz_yy_y_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yy_y_zz[i] = -2.0 * g_z_yy_y_zz[i] * a_exp + 4.0 * g_xxz_yy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (714-720)

    #pragma omp simd aligned(g_xxz_yy_z_xx, g_xxz_yy_z_xy, g_xxz_yy_z_xz, g_xxz_yy_z_yy, g_xxz_yy_z_yz, g_xxz_yy_z_zz, g_xz_0_0_0_x_yy_z_xx, g_xz_0_0_0_x_yy_z_xy, g_xz_0_0_0_x_yy_z_xz, g_xz_0_0_0_x_yy_z_yy, g_xz_0_0_0_x_yy_z_yz, g_xz_0_0_0_x_yy_z_zz, g_z_yy_z_xx, g_z_yy_z_xy, g_z_yy_z_xz, g_z_yy_z_yy, g_z_yy_z_yz, g_z_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_yy_z_xx[i] = -2.0 * g_z_yy_z_xx[i] * a_exp + 4.0 * g_xxz_yy_z_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yy_z_xy[i] = -2.0 * g_z_yy_z_xy[i] * a_exp + 4.0 * g_xxz_yy_z_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yy_z_xz[i] = -2.0 * g_z_yy_z_xz[i] * a_exp + 4.0 * g_xxz_yy_z_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yy_z_yy[i] = -2.0 * g_z_yy_z_yy[i] * a_exp + 4.0 * g_xxz_yy_z_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yy_z_yz[i] = -2.0 * g_z_yy_z_yz[i] * a_exp + 4.0 * g_xxz_yy_z_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yy_z_zz[i] = -2.0 * g_z_yy_z_zz[i] * a_exp + 4.0 * g_xxz_yy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (720-726)

    #pragma omp simd aligned(g_xxz_yz_x_xx, g_xxz_yz_x_xy, g_xxz_yz_x_xz, g_xxz_yz_x_yy, g_xxz_yz_x_yz, g_xxz_yz_x_zz, g_xz_0_0_0_x_yz_x_xx, g_xz_0_0_0_x_yz_x_xy, g_xz_0_0_0_x_yz_x_xz, g_xz_0_0_0_x_yz_x_yy, g_xz_0_0_0_x_yz_x_yz, g_xz_0_0_0_x_yz_x_zz, g_z_yz_x_xx, g_z_yz_x_xy, g_z_yz_x_xz, g_z_yz_x_yy, g_z_yz_x_yz, g_z_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_yz_x_xx[i] = -2.0 * g_z_yz_x_xx[i] * a_exp + 4.0 * g_xxz_yz_x_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yz_x_xy[i] = -2.0 * g_z_yz_x_xy[i] * a_exp + 4.0 * g_xxz_yz_x_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yz_x_xz[i] = -2.0 * g_z_yz_x_xz[i] * a_exp + 4.0 * g_xxz_yz_x_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yz_x_yy[i] = -2.0 * g_z_yz_x_yy[i] * a_exp + 4.0 * g_xxz_yz_x_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yz_x_yz[i] = -2.0 * g_z_yz_x_yz[i] * a_exp + 4.0 * g_xxz_yz_x_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yz_x_zz[i] = -2.0 * g_z_yz_x_zz[i] * a_exp + 4.0 * g_xxz_yz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (726-732)

    #pragma omp simd aligned(g_xxz_yz_y_xx, g_xxz_yz_y_xy, g_xxz_yz_y_xz, g_xxz_yz_y_yy, g_xxz_yz_y_yz, g_xxz_yz_y_zz, g_xz_0_0_0_x_yz_y_xx, g_xz_0_0_0_x_yz_y_xy, g_xz_0_0_0_x_yz_y_xz, g_xz_0_0_0_x_yz_y_yy, g_xz_0_0_0_x_yz_y_yz, g_xz_0_0_0_x_yz_y_zz, g_z_yz_y_xx, g_z_yz_y_xy, g_z_yz_y_xz, g_z_yz_y_yy, g_z_yz_y_yz, g_z_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_yz_y_xx[i] = -2.0 * g_z_yz_y_xx[i] * a_exp + 4.0 * g_xxz_yz_y_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yz_y_xy[i] = -2.0 * g_z_yz_y_xy[i] * a_exp + 4.0 * g_xxz_yz_y_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yz_y_xz[i] = -2.0 * g_z_yz_y_xz[i] * a_exp + 4.0 * g_xxz_yz_y_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yz_y_yy[i] = -2.0 * g_z_yz_y_yy[i] * a_exp + 4.0 * g_xxz_yz_y_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yz_y_yz[i] = -2.0 * g_z_yz_y_yz[i] * a_exp + 4.0 * g_xxz_yz_y_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yz_y_zz[i] = -2.0 * g_z_yz_y_zz[i] * a_exp + 4.0 * g_xxz_yz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (732-738)

    #pragma omp simd aligned(g_xxz_yz_z_xx, g_xxz_yz_z_xy, g_xxz_yz_z_xz, g_xxz_yz_z_yy, g_xxz_yz_z_yz, g_xxz_yz_z_zz, g_xz_0_0_0_x_yz_z_xx, g_xz_0_0_0_x_yz_z_xy, g_xz_0_0_0_x_yz_z_xz, g_xz_0_0_0_x_yz_z_yy, g_xz_0_0_0_x_yz_z_yz, g_xz_0_0_0_x_yz_z_zz, g_z_yz_z_xx, g_z_yz_z_xy, g_z_yz_z_xz, g_z_yz_z_yy, g_z_yz_z_yz, g_z_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_yz_z_xx[i] = -2.0 * g_z_yz_z_xx[i] * a_exp + 4.0 * g_xxz_yz_z_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yz_z_xy[i] = -2.0 * g_z_yz_z_xy[i] * a_exp + 4.0 * g_xxz_yz_z_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yz_z_xz[i] = -2.0 * g_z_yz_z_xz[i] * a_exp + 4.0 * g_xxz_yz_z_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yz_z_yy[i] = -2.0 * g_z_yz_z_yy[i] * a_exp + 4.0 * g_xxz_yz_z_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yz_z_yz[i] = -2.0 * g_z_yz_z_yz[i] * a_exp + 4.0 * g_xxz_yz_z_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yz_z_zz[i] = -2.0 * g_z_yz_z_zz[i] * a_exp + 4.0 * g_xxz_yz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (738-744)

    #pragma omp simd aligned(g_xxz_zz_x_xx, g_xxz_zz_x_xy, g_xxz_zz_x_xz, g_xxz_zz_x_yy, g_xxz_zz_x_yz, g_xxz_zz_x_zz, g_xz_0_0_0_x_zz_x_xx, g_xz_0_0_0_x_zz_x_xy, g_xz_0_0_0_x_zz_x_xz, g_xz_0_0_0_x_zz_x_yy, g_xz_0_0_0_x_zz_x_yz, g_xz_0_0_0_x_zz_x_zz, g_z_zz_x_xx, g_z_zz_x_xy, g_z_zz_x_xz, g_z_zz_x_yy, g_z_zz_x_yz, g_z_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_zz_x_xx[i] = -2.0 * g_z_zz_x_xx[i] * a_exp + 4.0 * g_xxz_zz_x_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_zz_x_xy[i] = -2.0 * g_z_zz_x_xy[i] * a_exp + 4.0 * g_xxz_zz_x_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_zz_x_xz[i] = -2.0 * g_z_zz_x_xz[i] * a_exp + 4.0 * g_xxz_zz_x_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_zz_x_yy[i] = -2.0 * g_z_zz_x_yy[i] * a_exp + 4.0 * g_xxz_zz_x_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_zz_x_yz[i] = -2.0 * g_z_zz_x_yz[i] * a_exp + 4.0 * g_xxz_zz_x_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_zz_x_zz[i] = -2.0 * g_z_zz_x_zz[i] * a_exp + 4.0 * g_xxz_zz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (744-750)

    #pragma omp simd aligned(g_xxz_zz_y_xx, g_xxz_zz_y_xy, g_xxz_zz_y_xz, g_xxz_zz_y_yy, g_xxz_zz_y_yz, g_xxz_zz_y_zz, g_xz_0_0_0_x_zz_y_xx, g_xz_0_0_0_x_zz_y_xy, g_xz_0_0_0_x_zz_y_xz, g_xz_0_0_0_x_zz_y_yy, g_xz_0_0_0_x_zz_y_yz, g_xz_0_0_0_x_zz_y_zz, g_z_zz_y_xx, g_z_zz_y_xy, g_z_zz_y_xz, g_z_zz_y_yy, g_z_zz_y_yz, g_z_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_zz_y_xx[i] = -2.0 * g_z_zz_y_xx[i] * a_exp + 4.0 * g_xxz_zz_y_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_zz_y_xy[i] = -2.0 * g_z_zz_y_xy[i] * a_exp + 4.0 * g_xxz_zz_y_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_zz_y_xz[i] = -2.0 * g_z_zz_y_xz[i] * a_exp + 4.0 * g_xxz_zz_y_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_zz_y_yy[i] = -2.0 * g_z_zz_y_yy[i] * a_exp + 4.0 * g_xxz_zz_y_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_zz_y_yz[i] = -2.0 * g_z_zz_y_yz[i] * a_exp + 4.0 * g_xxz_zz_y_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_zz_y_zz[i] = -2.0 * g_z_zz_y_zz[i] * a_exp + 4.0 * g_xxz_zz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (750-756)

    #pragma omp simd aligned(g_xxz_zz_z_xx, g_xxz_zz_z_xy, g_xxz_zz_z_xz, g_xxz_zz_z_yy, g_xxz_zz_z_yz, g_xxz_zz_z_zz, g_xz_0_0_0_x_zz_z_xx, g_xz_0_0_0_x_zz_z_xy, g_xz_0_0_0_x_zz_z_xz, g_xz_0_0_0_x_zz_z_yy, g_xz_0_0_0_x_zz_z_yz, g_xz_0_0_0_x_zz_z_zz, g_z_zz_z_xx, g_z_zz_z_xy, g_z_zz_z_xz, g_z_zz_z_yy, g_z_zz_z_yz, g_z_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_zz_z_xx[i] = -2.0 * g_z_zz_z_xx[i] * a_exp + 4.0 * g_xxz_zz_z_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_x_zz_z_xy[i] = -2.0 * g_z_zz_z_xy[i] * a_exp + 4.0 * g_xxz_zz_z_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_zz_z_xz[i] = -2.0 * g_z_zz_z_xz[i] * a_exp + 4.0 * g_xxz_zz_z_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_zz_z_yy[i] = -2.0 * g_z_zz_z_yy[i] * a_exp + 4.0 * g_xxz_zz_z_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_x_zz_z_yz[i] = -2.0 * g_z_zz_z_yz[i] * a_exp + 4.0 * g_xxz_zz_z_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_x_zz_z_zz[i] = -2.0 * g_z_zz_z_zz[i] * a_exp + 4.0 * g_xxz_zz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (756-762)

    #pragma omp simd aligned(g_xyz_xx_x_xx, g_xyz_xx_x_xy, g_xyz_xx_x_xz, g_xyz_xx_x_yy, g_xyz_xx_x_yz, g_xyz_xx_x_zz, g_xz_0_0_0_y_xx_x_xx, g_xz_0_0_0_y_xx_x_xy, g_xz_0_0_0_y_xx_x_xz, g_xz_0_0_0_y_xx_x_yy, g_xz_0_0_0_y_xx_x_yz, g_xz_0_0_0_y_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_xx_x_xx[i] = 4.0 * g_xyz_xx_x_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xx_x_xy[i] = 4.0 * g_xyz_xx_x_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xx_x_xz[i] = 4.0 * g_xyz_xx_x_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xx_x_yy[i] = 4.0 * g_xyz_xx_x_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xx_x_yz[i] = 4.0 * g_xyz_xx_x_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xx_x_zz[i] = 4.0 * g_xyz_xx_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (762-768)

    #pragma omp simd aligned(g_xyz_xx_y_xx, g_xyz_xx_y_xy, g_xyz_xx_y_xz, g_xyz_xx_y_yy, g_xyz_xx_y_yz, g_xyz_xx_y_zz, g_xz_0_0_0_y_xx_y_xx, g_xz_0_0_0_y_xx_y_xy, g_xz_0_0_0_y_xx_y_xz, g_xz_0_0_0_y_xx_y_yy, g_xz_0_0_0_y_xx_y_yz, g_xz_0_0_0_y_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_xx_y_xx[i] = 4.0 * g_xyz_xx_y_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xx_y_xy[i] = 4.0 * g_xyz_xx_y_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xx_y_xz[i] = 4.0 * g_xyz_xx_y_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xx_y_yy[i] = 4.0 * g_xyz_xx_y_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xx_y_yz[i] = 4.0 * g_xyz_xx_y_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xx_y_zz[i] = 4.0 * g_xyz_xx_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (768-774)

    #pragma omp simd aligned(g_xyz_xx_z_xx, g_xyz_xx_z_xy, g_xyz_xx_z_xz, g_xyz_xx_z_yy, g_xyz_xx_z_yz, g_xyz_xx_z_zz, g_xz_0_0_0_y_xx_z_xx, g_xz_0_0_0_y_xx_z_xy, g_xz_0_0_0_y_xx_z_xz, g_xz_0_0_0_y_xx_z_yy, g_xz_0_0_0_y_xx_z_yz, g_xz_0_0_0_y_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_xx_z_xx[i] = 4.0 * g_xyz_xx_z_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xx_z_xy[i] = 4.0 * g_xyz_xx_z_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xx_z_xz[i] = 4.0 * g_xyz_xx_z_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xx_z_yy[i] = 4.0 * g_xyz_xx_z_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xx_z_yz[i] = 4.0 * g_xyz_xx_z_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xx_z_zz[i] = 4.0 * g_xyz_xx_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (774-780)

    #pragma omp simd aligned(g_xyz_xy_x_xx, g_xyz_xy_x_xy, g_xyz_xy_x_xz, g_xyz_xy_x_yy, g_xyz_xy_x_yz, g_xyz_xy_x_zz, g_xz_0_0_0_y_xy_x_xx, g_xz_0_0_0_y_xy_x_xy, g_xz_0_0_0_y_xy_x_xz, g_xz_0_0_0_y_xy_x_yy, g_xz_0_0_0_y_xy_x_yz, g_xz_0_0_0_y_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_xy_x_xx[i] = 4.0 * g_xyz_xy_x_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xy_x_xy[i] = 4.0 * g_xyz_xy_x_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xy_x_xz[i] = 4.0 * g_xyz_xy_x_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xy_x_yy[i] = 4.0 * g_xyz_xy_x_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xy_x_yz[i] = 4.0 * g_xyz_xy_x_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xy_x_zz[i] = 4.0 * g_xyz_xy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (780-786)

    #pragma omp simd aligned(g_xyz_xy_y_xx, g_xyz_xy_y_xy, g_xyz_xy_y_xz, g_xyz_xy_y_yy, g_xyz_xy_y_yz, g_xyz_xy_y_zz, g_xz_0_0_0_y_xy_y_xx, g_xz_0_0_0_y_xy_y_xy, g_xz_0_0_0_y_xy_y_xz, g_xz_0_0_0_y_xy_y_yy, g_xz_0_0_0_y_xy_y_yz, g_xz_0_0_0_y_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_xy_y_xx[i] = 4.0 * g_xyz_xy_y_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xy_y_xy[i] = 4.0 * g_xyz_xy_y_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xy_y_xz[i] = 4.0 * g_xyz_xy_y_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xy_y_yy[i] = 4.0 * g_xyz_xy_y_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xy_y_yz[i] = 4.0 * g_xyz_xy_y_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xy_y_zz[i] = 4.0 * g_xyz_xy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (786-792)

    #pragma omp simd aligned(g_xyz_xy_z_xx, g_xyz_xy_z_xy, g_xyz_xy_z_xz, g_xyz_xy_z_yy, g_xyz_xy_z_yz, g_xyz_xy_z_zz, g_xz_0_0_0_y_xy_z_xx, g_xz_0_0_0_y_xy_z_xy, g_xz_0_0_0_y_xy_z_xz, g_xz_0_0_0_y_xy_z_yy, g_xz_0_0_0_y_xy_z_yz, g_xz_0_0_0_y_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_xy_z_xx[i] = 4.0 * g_xyz_xy_z_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xy_z_xy[i] = 4.0 * g_xyz_xy_z_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xy_z_xz[i] = 4.0 * g_xyz_xy_z_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xy_z_yy[i] = 4.0 * g_xyz_xy_z_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xy_z_yz[i] = 4.0 * g_xyz_xy_z_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xy_z_zz[i] = 4.0 * g_xyz_xy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (792-798)

    #pragma omp simd aligned(g_xyz_xz_x_xx, g_xyz_xz_x_xy, g_xyz_xz_x_xz, g_xyz_xz_x_yy, g_xyz_xz_x_yz, g_xyz_xz_x_zz, g_xz_0_0_0_y_xz_x_xx, g_xz_0_0_0_y_xz_x_xy, g_xz_0_0_0_y_xz_x_xz, g_xz_0_0_0_y_xz_x_yy, g_xz_0_0_0_y_xz_x_yz, g_xz_0_0_0_y_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_xz_x_xx[i] = 4.0 * g_xyz_xz_x_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xz_x_xy[i] = 4.0 * g_xyz_xz_x_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xz_x_xz[i] = 4.0 * g_xyz_xz_x_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xz_x_yy[i] = 4.0 * g_xyz_xz_x_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xz_x_yz[i] = 4.0 * g_xyz_xz_x_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xz_x_zz[i] = 4.0 * g_xyz_xz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (798-804)

    #pragma omp simd aligned(g_xyz_xz_y_xx, g_xyz_xz_y_xy, g_xyz_xz_y_xz, g_xyz_xz_y_yy, g_xyz_xz_y_yz, g_xyz_xz_y_zz, g_xz_0_0_0_y_xz_y_xx, g_xz_0_0_0_y_xz_y_xy, g_xz_0_0_0_y_xz_y_xz, g_xz_0_0_0_y_xz_y_yy, g_xz_0_0_0_y_xz_y_yz, g_xz_0_0_0_y_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_xz_y_xx[i] = 4.0 * g_xyz_xz_y_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xz_y_xy[i] = 4.0 * g_xyz_xz_y_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xz_y_xz[i] = 4.0 * g_xyz_xz_y_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xz_y_yy[i] = 4.0 * g_xyz_xz_y_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xz_y_yz[i] = 4.0 * g_xyz_xz_y_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xz_y_zz[i] = 4.0 * g_xyz_xz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (804-810)

    #pragma omp simd aligned(g_xyz_xz_z_xx, g_xyz_xz_z_xy, g_xyz_xz_z_xz, g_xyz_xz_z_yy, g_xyz_xz_z_yz, g_xyz_xz_z_zz, g_xz_0_0_0_y_xz_z_xx, g_xz_0_0_0_y_xz_z_xy, g_xz_0_0_0_y_xz_z_xz, g_xz_0_0_0_y_xz_z_yy, g_xz_0_0_0_y_xz_z_yz, g_xz_0_0_0_y_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_xz_z_xx[i] = 4.0 * g_xyz_xz_z_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xz_z_xy[i] = 4.0 * g_xyz_xz_z_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xz_z_xz[i] = 4.0 * g_xyz_xz_z_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xz_z_yy[i] = 4.0 * g_xyz_xz_z_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xz_z_yz[i] = 4.0 * g_xyz_xz_z_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xz_z_zz[i] = 4.0 * g_xyz_xz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (810-816)

    #pragma omp simd aligned(g_xyz_yy_x_xx, g_xyz_yy_x_xy, g_xyz_yy_x_xz, g_xyz_yy_x_yy, g_xyz_yy_x_yz, g_xyz_yy_x_zz, g_xz_0_0_0_y_yy_x_xx, g_xz_0_0_0_y_yy_x_xy, g_xz_0_0_0_y_yy_x_xz, g_xz_0_0_0_y_yy_x_yy, g_xz_0_0_0_y_yy_x_yz, g_xz_0_0_0_y_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_yy_x_xx[i] = 4.0 * g_xyz_yy_x_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yy_x_xy[i] = 4.0 * g_xyz_yy_x_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yy_x_xz[i] = 4.0 * g_xyz_yy_x_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yy_x_yy[i] = 4.0 * g_xyz_yy_x_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yy_x_yz[i] = 4.0 * g_xyz_yy_x_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yy_x_zz[i] = 4.0 * g_xyz_yy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (816-822)

    #pragma omp simd aligned(g_xyz_yy_y_xx, g_xyz_yy_y_xy, g_xyz_yy_y_xz, g_xyz_yy_y_yy, g_xyz_yy_y_yz, g_xyz_yy_y_zz, g_xz_0_0_0_y_yy_y_xx, g_xz_0_0_0_y_yy_y_xy, g_xz_0_0_0_y_yy_y_xz, g_xz_0_0_0_y_yy_y_yy, g_xz_0_0_0_y_yy_y_yz, g_xz_0_0_0_y_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_yy_y_xx[i] = 4.0 * g_xyz_yy_y_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yy_y_xy[i] = 4.0 * g_xyz_yy_y_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yy_y_xz[i] = 4.0 * g_xyz_yy_y_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yy_y_yy[i] = 4.0 * g_xyz_yy_y_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yy_y_yz[i] = 4.0 * g_xyz_yy_y_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yy_y_zz[i] = 4.0 * g_xyz_yy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (822-828)

    #pragma omp simd aligned(g_xyz_yy_z_xx, g_xyz_yy_z_xy, g_xyz_yy_z_xz, g_xyz_yy_z_yy, g_xyz_yy_z_yz, g_xyz_yy_z_zz, g_xz_0_0_0_y_yy_z_xx, g_xz_0_0_0_y_yy_z_xy, g_xz_0_0_0_y_yy_z_xz, g_xz_0_0_0_y_yy_z_yy, g_xz_0_0_0_y_yy_z_yz, g_xz_0_0_0_y_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_yy_z_xx[i] = 4.0 * g_xyz_yy_z_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yy_z_xy[i] = 4.0 * g_xyz_yy_z_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yy_z_xz[i] = 4.0 * g_xyz_yy_z_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yy_z_yy[i] = 4.0 * g_xyz_yy_z_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yy_z_yz[i] = 4.0 * g_xyz_yy_z_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yy_z_zz[i] = 4.0 * g_xyz_yy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (828-834)

    #pragma omp simd aligned(g_xyz_yz_x_xx, g_xyz_yz_x_xy, g_xyz_yz_x_xz, g_xyz_yz_x_yy, g_xyz_yz_x_yz, g_xyz_yz_x_zz, g_xz_0_0_0_y_yz_x_xx, g_xz_0_0_0_y_yz_x_xy, g_xz_0_0_0_y_yz_x_xz, g_xz_0_0_0_y_yz_x_yy, g_xz_0_0_0_y_yz_x_yz, g_xz_0_0_0_y_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_yz_x_xx[i] = 4.0 * g_xyz_yz_x_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yz_x_xy[i] = 4.0 * g_xyz_yz_x_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yz_x_xz[i] = 4.0 * g_xyz_yz_x_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yz_x_yy[i] = 4.0 * g_xyz_yz_x_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yz_x_yz[i] = 4.0 * g_xyz_yz_x_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yz_x_zz[i] = 4.0 * g_xyz_yz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (834-840)

    #pragma omp simd aligned(g_xyz_yz_y_xx, g_xyz_yz_y_xy, g_xyz_yz_y_xz, g_xyz_yz_y_yy, g_xyz_yz_y_yz, g_xyz_yz_y_zz, g_xz_0_0_0_y_yz_y_xx, g_xz_0_0_0_y_yz_y_xy, g_xz_0_0_0_y_yz_y_xz, g_xz_0_0_0_y_yz_y_yy, g_xz_0_0_0_y_yz_y_yz, g_xz_0_0_0_y_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_yz_y_xx[i] = 4.0 * g_xyz_yz_y_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yz_y_xy[i] = 4.0 * g_xyz_yz_y_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yz_y_xz[i] = 4.0 * g_xyz_yz_y_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yz_y_yy[i] = 4.0 * g_xyz_yz_y_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yz_y_yz[i] = 4.0 * g_xyz_yz_y_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yz_y_zz[i] = 4.0 * g_xyz_yz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (840-846)

    #pragma omp simd aligned(g_xyz_yz_z_xx, g_xyz_yz_z_xy, g_xyz_yz_z_xz, g_xyz_yz_z_yy, g_xyz_yz_z_yz, g_xyz_yz_z_zz, g_xz_0_0_0_y_yz_z_xx, g_xz_0_0_0_y_yz_z_xy, g_xz_0_0_0_y_yz_z_xz, g_xz_0_0_0_y_yz_z_yy, g_xz_0_0_0_y_yz_z_yz, g_xz_0_0_0_y_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_yz_z_xx[i] = 4.0 * g_xyz_yz_z_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yz_z_xy[i] = 4.0 * g_xyz_yz_z_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yz_z_xz[i] = 4.0 * g_xyz_yz_z_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yz_z_yy[i] = 4.0 * g_xyz_yz_z_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yz_z_yz[i] = 4.0 * g_xyz_yz_z_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yz_z_zz[i] = 4.0 * g_xyz_yz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (846-852)

    #pragma omp simd aligned(g_xyz_zz_x_xx, g_xyz_zz_x_xy, g_xyz_zz_x_xz, g_xyz_zz_x_yy, g_xyz_zz_x_yz, g_xyz_zz_x_zz, g_xz_0_0_0_y_zz_x_xx, g_xz_0_0_0_y_zz_x_xy, g_xz_0_0_0_y_zz_x_xz, g_xz_0_0_0_y_zz_x_yy, g_xz_0_0_0_y_zz_x_yz, g_xz_0_0_0_y_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_zz_x_xx[i] = 4.0 * g_xyz_zz_x_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_zz_x_xy[i] = 4.0 * g_xyz_zz_x_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_zz_x_xz[i] = 4.0 * g_xyz_zz_x_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_zz_x_yy[i] = 4.0 * g_xyz_zz_x_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_zz_x_yz[i] = 4.0 * g_xyz_zz_x_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_zz_x_zz[i] = 4.0 * g_xyz_zz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (852-858)

    #pragma omp simd aligned(g_xyz_zz_y_xx, g_xyz_zz_y_xy, g_xyz_zz_y_xz, g_xyz_zz_y_yy, g_xyz_zz_y_yz, g_xyz_zz_y_zz, g_xz_0_0_0_y_zz_y_xx, g_xz_0_0_0_y_zz_y_xy, g_xz_0_0_0_y_zz_y_xz, g_xz_0_0_0_y_zz_y_yy, g_xz_0_0_0_y_zz_y_yz, g_xz_0_0_0_y_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_zz_y_xx[i] = 4.0 * g_xyz_zz_y_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_zz_y_xy[i] = 4.0 * g_xyz_zz_y_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_zz_y_xz[i] = 4.0 * g_xyz_zz_y_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_zz_y_yy[i] = 4.0 * g_xyz_zz_y_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_zz_y_yz[i] = 4.0 * g_xyz_zz_y_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_zz_y_zz[i] = 4.0 * g_xyz_zz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (858-864)

    #pragma omp simd aligned(g_xyz_zz_z_xx, g_xyz_zz_z_xy, g_xyz_zz_z_xz, g_xyz_zz_z_yy, g_xyz_zz_z_yz, g_xyz_zz_z_zz, g_xz_0_0_0_y_zz_z_xx, g_xz_0_0_0_y_zz_z_xy, g_xz_0_0_0_y_zz_z_xz, g_xz_0_0_0_y_zz_z_yy, g_xz_0_0_0_y_zz_z_yz, g_xz_0_0_0_y_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_zz_z_xx[i] = 4.0 * g_xyz_zz_z_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_y_zz_z_xy[i] = 4.0 * g_xyz_zz_z_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_zz_z_xz[i] = 4.0 * g_xyz_zz_z_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_zz_z_yy[i] = 4.0 * g_xyz_zz_z_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_y_zz_z_yz[i] = 4.0 * g_xyz_zz_z_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_y_zz_z_zz[i] = 4.0 * g_xyz_zz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (864-870)

    #pragma omp simd aligned(g_x_xx_x_xx, g_x_xx_x_xy, g_x_xx_x_xz, g_x_xx_x_yy, g_x_xx_x_yz, g_x_xx_x_zz, g_xz_0_0_0_z_xx_x_xx, g_xz_0_0_0_z_xx_x_xy, g_xz_0_0_0_z_xx_x_xz, g_xz_0_0_0_z_xx_x_yy, g_xz_0_0_0_z_xx_x_yz, g_xz_0_0_0_z_xx_x_zz, g_xzz_xx_x_xx, g_xzz_xx_x_xy, g_xzz_xx_x_xz, g_xzz_xx_x_yy, g_xzz_xx_x_yz, g_xzz_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_xx_x_xx[i] = -2.0 * g_x_xx_x_xx[i] * a_exp + 4.0 * g_xzz_xx_x_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xx_x_xy[i] = -2.0 * g_x_xx_x_xy[i] * a_exp + 4.0 * g_xzz_xx_x_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xx_x_xz[i] = -2.0 * g_x_xx_x_xz[i] * a_exp + 4.0 * g_xzz_xx_x_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xx_x_yy[i] = -2.0 * g_x_xx_x_yy[i] * a_exp + 4.0 * g_xzz_xx_x_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xx_x_yz[i] = -2.0 * g_x_xx_x_yz[i] * a_exp + 4.0 * g_xzz_xx_x_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xx_x_zz[i] = -2.0 * g_x_xx_x_zz[i] * a_exp + 4.0 * g_xzz_xx_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (870-876)

    #pragma omp simd aligned(g_x_xx_y_xx, g_x_xx_y_xy, g_x_xx_y_xz, g_x_xx_y_yy, g_x_xx_y_yz, g_x_xx_y_zz, g_xz_0_0_0_z_xx_y_xx, g_xz_0_0_0_z_xx_y_xy, g_xz_0_0_0_z_xx_y_xz, g_xz_0_0_0_z_xx_y_yy, g_xz_0_0_0_z_xx_y_yz, g_xz_0_0_0_z_xx_y_zz, g_xzz_xx_y_xx, g_xzz_xx_y_xy, g_xzz_xx_y_xz, g_xzz_xx_y_yy, g_xzz_xx_y_yz, g_xzz_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_xx_y_xx[i] = -2.0 * g_x_xx_y_xx[i] * a_exp + 4.0 * g_xzz_xx_y_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xx_y_xy[i] = -2.0 * g_x_xx_y_xy[i] * a_exp + 4.0 * g_xzz_xx_y_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xx_y_xz[i] = -2.0 * g_x_xx_y_xz[i] * a_exp + 4.0 * g_xzz_xx_y_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xx_y_yy[i] = -2.0 * g_x_xx_y_yy[i] * a_exp + 4.0 * g_xzz_xx_y_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xx_y_yz[i] = -2.0 * g_x_xx_y_yz[i] * a_exp + 4.0 * g_xzz_xx_y_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xx_y_zz[i] = -2.0 * g_x_xx_y_zz[i] * a_exp + 4.0 * g_xzz_xx_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (876-882)

    #pragma omp simd aligned(g_x_xx_z_xx, g_x_xx_z_xy, g_x_xx_z_xz, g_x_xx_z_yy, g_x_xx_z_yz, g_x_xx_z_zz, g_xz_0_0_0_z_xx_z_xx, g_xz_0_0_0_z_xx_z_xy, g_xz_0_0_0_z_xx_z_xz, g_xz_0_0_0_z_xx_z_yy, g_xz_0_0_0_z_xx_z_yz, g_xz_0_0_0_z_xx_z_zz, g_xzz_xx_z_xx, g_xzz_xx_z_xy, g_xzz_xx_z_xz, g_xzz_xx_z_yy, g_xzz_xx_z_yz, g_xzz_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_xx_z_xx[i] = -2.0 * g_x_xx_z_xx[i] * a_exp + 4.0 * g_xzz_xx_z_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xx_z_xy[i] = -2.0 * g_x_xx_z_xy[i] * a_exp + 4.0 * g_xzz_xx_z_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xx_z_xz[i] = -2.0 * g_x_xx_z_xz[i] * a_exp + 4.0 * g_xzz_xx_z_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xx_z_yy[i] = -2.0 * g_x_xx_z_yy[i] * a_exp + 4.0 * g_xzz_xx_z_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xx_z_yz[i] = -2.0 * g_x_xx_z_yz[i] * a_exp + 4.0 * g_xzz_xx_z_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xx_z_zz[i] = -2.0 * g_x_xx_z_zz[i] * a_exp + 4.0 * g_xzz_xx_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (882-888)

    #pragma omp simd aligned(g_x_xy_x_xx, g_x_xy_x_xy, g_x_xy_x_xz, g_x_xy_x_yy, g_x_xy_x_yz, g_x_xy_x_zz, g_xz_0_0_0_z_xy_x_xx, g_xz_0_0_0_z_xy_x_xy, g_xz_0_0_0_z_xy_x_xz, g_xz_0_0_0_z_xy_x_yy, g_xz_0_0_0_z_xy_x_yz, g_xz_0_0_0_z_xy_x_zz, g_xzz_xy_x_xx, g_xzz_xy_x_xy, g_xzz_xy_x_xz, g_xzz_xy_x_yy, g_xzz_xy_x_yz, g_xzz_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_xy_x_xx[i] = -2.0 * g_x_xy_x_xx[i] * a_exp + 4.0 * g_xzz_xy_x_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xy_x_xy[i] = -2.0 * g_x_xy_x_xy[i] * a_exp + 4.0 * g_xzz_xy_x_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xy_x_xz[i] = -2.0 * g_x_xy_x_xz[i] * a_exp + 4.0 * g_xzz_xy_x_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xy_x_yy[i] = -2.0 * g_x_xy_x_yy[i] * a_exp + 4.0 * g_xzz_xy_x_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xy_x_yz[i] = -2.0 * g_x_xy_x_yz[i] * a_exp + 4.0 * g_xzz_xy_x_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xy_x_zz[i] = -2.0 * g_x_xy_x_zz[i] * a_exp + 4.0 * g_xzz_xy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (888-894)

    #pragma omp simd aligned(g_x_xy_y_xx, g_x_xy_y_xy, g_x_xy_y_xz, g_x_xy_y_yy, g_x_xy_y_yz, g_x_xy_y_zz, g_xz_0_0_0_z_xy_y_xx, g_xz_0_0_0_z_xy_y_xy, g_xz_0_0_0_z_xy_y_xz, g_xz_0_0_0_z_xy_y_yy, g_xz_0_0_0_z_xy_y_yz, g_xz_0_0_0_z_xy_y_zz, g_xzz_xy_y_xx, g_xzz_xy_y_xy, g_xzz_xy_y_xz, g_xzz_xy_y_yy, g_xzz_xy_y_yz, g_xzz_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_xy_y_xx[i] = -2.0 * g_x_xy_y_xx[i] * a_exp + 4.0 * g_xzz_xy_y_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xy_y_xy[i] = -2.0 * g_x_xy_y_xy[i] * a_exp + 4.0 * g_xzz_xy_y_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xy_y_xz[i] = -2.0 * g_x_xy_y_xz[i] * a_exp + 4.0 * g_xzz_xy_y_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xy_y_yy[i] = -2.0 * g_x_xy_y_yy[i] * a_exp + 4.0 * g_xzz_xy_y_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xy_y_yz[i] = -2.0 * g_x_xy_y_yz[i] * a_exp + 4.0 * g_xzz_xy_y_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xy_y_zz[i] = -2.0 * g_x_xy_y_zz[i] * a_exp + 4.0 * g_xzz_xy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (894-900)

    #pragma omp simd aligned(g_x_xy_z_xx, g_x_xy_z_xy, g_x_xy_z_xz, g_x_xy_z_yy, g_x_xy_z_yz, g_x_xy_z_zz, g_xz_0_0_0_z_xy_z_xx, g_xz_0_0_0_z_xy_z_xy, g_xz_0_0_0_z_xy_z_xz, g_xz_0_0_0_z_xy_z_yy, g_xz_0_0_0_z_xy_z_yz, g_xz_0_0_0_z_xy_z_zz, g_xzz_xy_z_xx, g_xzz_xy_z_xy, g_xzz_xy_z_xz, g_xzz_xy_z_yy, g_xzz_xy_z_yz, g_xzz_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_xy_z_xx[i] = -2.0 * g_x_xy_z_xx[i] * a_exp + 4.0 * g_xzz_xy_z_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xy_z_xy[i] = -2.0 * g_x_xy_z_xy[i] * a_exp + 4.0 * g_xzz_xy_z_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xy_z_xz[i] = -2.0 * g_x_xy_z_xz[i] * a_exp + 4.0 * g_xzz_xy_z_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xy_z_yy[i] = -2.0 * g_x_xy_z_yy[i] * a_exp + 4.0 * g_xzz_xy_z_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xy_z_yz[i] = -2.0 * g_x_xy_z_yz[i] * a_exp + 4.0 * g_xzz_xy_z_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xy_z_zz[i] = -2.0 * g_x_xy_z_zz[i] * a_exp + 4.0 * g_xzz_xy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (900-906)

    #pragma omp simd aligned(g_x_xz_x_xx, g_x_xz_x_xy, g_x_xz_x_xz, g_x_xz_x_yy, g_x_xz_x_yz, g_x_xz_x_zz, g_xz_0_0_0_z_xz_x_xx, g_xz_0_0_0_z_xz_x_xy, g_xz_0_0_0_z_xz_x_xz, g_xz_0_0_0_z_xz_x_yy, g_xz_0_0_0_z_xz_x_yz, g_xz_0_0_0_z_xz_x_zz, g_xzz_xz_x_xx, g_xzz_xz_x_xy, g_xzz_xz_x_xz, g_xzz_xz_x_yy, g_xzz_xz_x_yz, g_xzz_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_xz_x_xx[i] = -2.0 * g_x_xz_x_xx[i] * a_exp + 4.0 * g_xzz_xz_x_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xz_x_xy[i] = -2.0 * g_x_xz_x_xy[i] * a_exp + 4.0 * g_xzz_xz_x_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xz_x_xz[i] = -2.0 * g_x_xz_x_xz[i] * a_exp + 4.0 * g_xzz_xz_x_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xz_x_yy[i] = -2.0 * g_x_xz_x_yy[i] * a_exp + 4.0 * g_xzz_xz_x_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xz_x_yz[i] = -2.0 * g_x_xz_x_yz[i] * a_exp + 4.0 * g_xzz_xz_x_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xz_x_zz[i] = -2.0 * g_x_xz_x_zz[i] * a_exp + 4.0 * g_xzz_xz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (906-912)

    #pragma omp simd aligned(g_x_xz_y_xx, g_x_xz_y_xy, g_x_xz_y_xz, g_x_xz_y_yy, g_x_xz_y_yz, g_x_xz_y_zz, g_xz_0_0_0_z_xz_y_xx, g_xz_0_0_0_z_xz_y_xy, g_xz_0_0_0_z_xz_y_xz, g_xz_0_0_0_z_xz_y_yy, g_xz_0_0_0_z_xz_y_yz, g_xz_0_0_0_z_xz_y_zz, g_xzz_xz_y_xx, g_xzz_xz_y_xy, g_xzz_xz_y_xz, g_xzz_xz_y_yy, g_xzz_xz_y_yz, g_xzz_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_xz_y_xx[i] = -2.0 * g_x_xz_y_xx[i] * a_exp + 4.0 * g_xzz_xz_y_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xz_y_xy[i] = -2.0 * g_x_xz_y_xy[i] * a_exp + 4.0 * g_xzz_xz_y_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xz_y_xz[i] = -2.0 * g_x_xz_y_xz[i] * a_exp + 4.0 * g_xzz_xz_y_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xz_y_yy[i] = -2.0 * g_x_xz_y_yy[i] * a_exp + 4.0 * g_xzz_xz_y_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xz_y_yz[i] = -2.0 * g_x_xz_y_yz[i] * a_exp + 4.0 * g_xzz_xz_y_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xz_y_zz[i] = -2.0 * g_x_xz_y_zz[i] * a_exp + 4.0 * g_xzz_xz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (912-918)

    #pragma omp simd aligned(g_x_xz_z_xx, g_x_xz_z_xy, g_x_xz_z_xz, g_x_xz_z_yy, g_x_xz_z_yz, g_x_xz_z_zz, g_xz_0_0_0_z_xz_z_xx, g_xz_0_0_0_z_xz_z_xy, g_xz_0_0_0_z_xz_z_xz, g_xz_0_0_0_z_xz_z_yy, g_xz_0_0_0_z_xz_z_yz, g_xz_0_0_0_z_xz_z_zz, g_xzz_xz_z_xx, g_xzz_xz_z_xy, g_xzz_xz_z_xz, g_xzz_xz_z_yy, g_xzz_xz_z_yz, g_xzz_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_xz_z_xx[i] = -2.0 * g_x_xz_z_xx[i] * a_exp + 4.0 * g_xzz_xz_z_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xz_z_xy[i] = -2.0 * g_x_xz_z_xy[i] * a_exp + 4.0 * g_xzz_xz_z_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xz_z_xz[i] = -2.0 * g_x_xz_z_xz[i] * a_exp + 4.0 * g_xzz_xz_z_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xz_z_yy[i] = -2.0 * g_x_xz_z_yy[i] * a_exp + 4.0 * g_xzz_xz_z_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xz_z_yz[i] = -2.0 * g_x_xz_z_yz[i] * a_exp + 4.0 * g_xzz_xz_z_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xz_z_zz[i] = -2.0 * g_x_xz_z_zz[i] * a_exp + 4.0 * g_xzz_xz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (918-924)

    #pragma omp simd aligned(g_x_yy_x_xx, g_x_yy_x_xy, g_x_yy_x_xz, g_x_yy_x_yy, g_x_yy_x_yz, g_x_yy_x_zz, g_xz_0_0_0_z_yy_x_xx, g_xz_0_0_0_z_yy_x_xy, g_xz_0_0_0_z_yy_x_xz, g_xz_0_0_0_z_yy_x_yy, g_xz_0_0_0_z_yy_x_yz, g_xz_0_0_0_z_yy_x_zz, g_xzz_yy_x_xx, g_xzz_yy_x_xy, g_xzz_yy_x_xz, g_xzz_yy_x_yy, g_xzz_yy_x_yz, g_xzz_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_yy_x_xx[i] = -2.0 * g_x_yy_x_xx[i] * a_exp + 4.0 * g_xzz_yy_x_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yy_x_xy[i] = -2.0 * g_x_yy_x_xy[i] * a_exp + 4.0 * g_xzz_yy_x_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yy_x_xz[i] = -2.0 * g_x_yy_x_xz[i] * a_exp + 4.0 * g_xzz_yy_x_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yy_x_yy[i] = -2.0 * g_x_yy_x_yy[i] * a_exp + 4.0 * g_xzz_yy_x_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yy_x_yz[i] = -2.0 * g_x_yy_x_yz[i] * a_exp + 4.0 * g_xzz_yy_x_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yy_x_zz[i] = -2.0 * g_x_yy_x_zz[i] * a_exp + 4.0 * g_xzz_yy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (924-930)

    #pragma omp simd aligned(g_x_yy_y_xx, g_x_yy_y_xy, g_x_yy_y_xz, g_x_yy_y_yy, g_x_yy_y_yz, g_x_yy_y_zz, g_xz_0_0_0_z_yy_y_xx, g_xz_0_0_0_z_yy_y_xy, g_xz_0_0_0_z_yy_y_xz, g_xz_0_0_0_z_yy_y_yy, g_xz_0_0_0_z_yy_y_yz, g_xz_0_0_0_z_yy_y_zz, g_xzz_yy_y_xx, g_xzz_yy_y_xy, g_xzz_yy_y_xz, g_xzz_yy_y_yy, g_xzz_yy_y_yz, g_xzz_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_yy_y_xx[i] = -2.0 * g_x_yy_y_xx[i] * a_exp + 4.0 * g_xzz_yy_y_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yy_y_xy[i] = -2.0 * g_x_yy_y_xy[i] * a_exp + 4.0 * g_xzz_yy_y_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yy_y_xz[i] = -2.0 * g_x_yy_y_xz[i] * a_exp + 4.0 * g_xzz_yy_y_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yy_y_yy[i] = -2.0 * g_x_yy_y_yy[i] * a_exp + 4.0 * g_xzz_yy_y_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yy_y_yz[i] = -2.0 * g_x_yy_y_yz[i] * a_exp + 4.0 * g_xzz_yy_y_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yy_y_zz[i] = -2.0 * g_x_yy_y_zz[i] * a_exp + 4.0 * g_xzz_yy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (930-936)

    #pragma omp simd aligned(g_x_yy_z_xx, g_x_yy_z_xy, g_x_yy_z_xz, g_x_yy_z_yy, g_x_yy_z_yz, g_x_yy_z_zz, g_xz_0_0_0_z_yy_z_xx, g_xz_0_0_0_z_yy_z_xy, g_xz_0_0_0_z_yy_z_xz, g_xz_0_0_0_z_yy_z_yy, g_xz_0_0_0_z_yy_z_yz, g_xz_0_0_0_z_yy_z_zz, g_xzz_yy_z_xx, g_xzz_yy_z_xy, g_xzz_yy_z_xz, g_xzz_yy_z_yy, g_xzz_yy_z_yz, g_xzz_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_yy_z_xx[i] = -2.0 * g_x_yy_z_xx[i] * a_exp + 4.0 * g_xzz_yy_z_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yy_z_xy[i] = -2.0 * g_x_yy_z_xy[i] * a_exp + 4.0 * g_xzz_yy_z_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yy_z_xz[i] = -2.0 * g_x_yy_z_xz[i] * a_exp + 4.0 * g_xzz_yy_z_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yy_z_yy[i] = -2.0 * g_x_yy_z_yy[i] * a_exp + 4.0 * g_xzz_yy_z_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yy_z_yz[i] = -2.0 * g_x_yy_z_yz[i] * a_exp + 4.0 * g_xzz_yy_z_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yy_z_zz[i] = -2.0 * g_x_yy_z_zz[i] * a_exp + 4.0 * g_xzz_yy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (936-942)

    #pragma omp simd aligned(g_x_yz_x_xx, g_x_yz_x_xy, g_x_yz_x_xz, g_x_yz_x_yy, g_x_yz_x_yz, g_x_yz_x_zz, g_xz_0_0_0_z_yz_x_xx, g_xz_0_0_0_z_yz_x_xy, g_xz_0_0_0_z_yz_x_xz, g_xz_0_0_0_z_yz_x_yy, g_xz_0_0_0_z_yz_x_yz, g_xz_0_0_0_z_yz_x_zz, g_xzz_yz_x_xx, g_xzz_yz_x_xy, g_xzz_yz_x_xz, g_xzz_yz_x_yy, g_xzz_yz_x_yz, g_xzz_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_yz_x_xx[i] = -2.0 * g_x_yz_x_xx[i] * a_exp + 4.0 * g_xzz_yz_x_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yz_x_xy[i] = -2.0 * g_x_yz_x_xy[i] * a_exp + 4.0 * g_xzz_yz_x_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yz_x_xz[i] = -2.0 * g_x_yz_x_xz[i] * a_exp + 4.0 * g_xzz_yz_x_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yz_x_yy[i] = -2.0 * g_x_yz_x_yy[i] * a_exp + 4.0 * g_xzz_yz_x_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yz_x_yz[i] = -2.0 * g_x_yz_x_yz[i] * a_exp + 4.0 * g_xzz_yz_x_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yz_x_zz[i] = -2.0 * g_x_yz_x_zz[i] * a_exp + 4.0 * g_xzz_yz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (942-948)

    #pragma omp simd aligned(g_x_yz_y_xx, g_x_yz_y_xy, g_x_yz_y_xz, g_x_yz_y_yy, g_x_yz_y_yz, g_x_yz_y_zz, g_xz_0_0_0_z_yz_y_xx, g_xz_0_0_0_z_yz_y_xy, g_xz_0_0_0_z_yz_y_xz, g_xz_0_0_0_z_yz_y_yy, g_xz_0_0_0_z_yz_y_yz, g_xz_0_0_0_z_yz_y_zz, g_xzz_yz_y_xx, g_xzz_yz_y_xy, g_xzz_yz_y_xz, g_xzz_yz_y_yy, g_xzz_yz_y_yz, g_xzz_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_yz_y_xx[i] = -2.0 * g_x_yz_y_xx[i] * a_exp + 4.0 * g_xzz_yz_y_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yz_y_xy[i] = -2.0 * g_x_yz_y_xy[i] * a_exp + 4.0 * g_xzz_yz_y_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yz_y_xz[i] = -2.0 * g_x_yz_y_xz[i] * a_exp + 4.0 * g_xzz_yz_y_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yz_y_yy[i] = -2.0 * g_x_yz_y_yy[i] * a_exp + 4.0 * g_xzz_yz_y_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yz_y_yz[i] = -2.0 * g_x_yz_y_yz[i] * a_exp + 4.0 * g_xzz_yz_y_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yz_y_zz[i] = -2.0 * g_x_yz_y_zz[i] * a_exp + 4.0 * g_xzz_yz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (948-954)

    #pragma omp simd aligned(g_x_yz_z_xx, g_x_yz_z_xy, g_x_yz_z_xz, g_x_yz_z_yy, g_x_yz_z_yz, g_x_yz_z_zz, g_xz_0_0_0_z_yz_z_xx, g_xz_0_0_0_z_yz_z_xy, g_xz_0_0_0_z_yz_z_xz, g_xz_0_0_0_z_yz_z_yy, g_xz_0_0_0_z_yz_z_yz, g_xz_0_0_0_z_yz_z_zz, g_xzz_yz_z_xx, g_xzz_yz_z_xy, g_xzz_yz_z_xz, g_xzz_yz_z_yy, g_xzz_yz_z_yz, g_xzz_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_yz_z_xx[i] = -2.0 * g_x_yz_z_xx[i] * a_exp + 4.0 * g_xzz_yz_z_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yz_z_xy[i] = -2.0 * g_x_yz_z_xy[i] * a_exp + 4.0 * g_xzz_yz_z_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yz_z_xz[i] = -2.0 * g_x_yz_z_xz[i] * a_exp + 4.0 * g_xzz_yz_z_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yz_z_yy[i] = -2.0 * g_x_yz_z_yy[i] * a_exp + 4.0 * g_xzz_yz_z_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yz_z_yz[i] = -2.0 * g_x_yz_z_yz[i] * a_exp + 4.0 * g_xzz_yz_z_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yz_z_zz[i] = -2.0 * g_x_yz_z_zz[i] * a_exp + 4.0 * g_xzz_yz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (954-960)

    #pragma omp simd aligned(g_x_zz_x_xx, g_x_zz_x_xy, g_x_zz_x_xz, g_x_zz_x_yy, g_x_zz_x_yz, g_x_zz_x_zz, g_xz_0_0_0_z_zz_x_xx, g_xz_0_0_0_z_zz_x_xy, g_xz_0_0_0_z_zz_x_xz, g_xz_0_0_0_z_zz_x_yy, g_xz_0_0_0_z_zz_x_yz, g_xz_0_0_0_z_zz_x_zz, g_xzz_zz_x_xx, g_xzz_zz_x_xy, g_xzz_zz_x_xz, g_xzz_zz_x_yy, g_xzz_zz_x_yz, g_xzz_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_zz_x_xx[i] = -2.0 * g_x_zz_x_xx[i] * a_exp + 4.0 * g_xzz_zz_x_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_zz_x_xy[i] = -2.0 * g_x_zz_x_xy[i] * a_exp + 4.0 * g_xzz_zz_x_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_zz_x_xz[i] = -2.0 * g_x_zz_x_xz[i] * a_exp + 4.0 * g_xzz_zz_x_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_zz_x_yy[i] = -2.0 * g_x_zz_x_yy[i] * a_exp + 4.0 * g_xzz_zz_x_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_zz_x_yz[i] = -2.0 * g_x_zz_x_yz[i] * a_exp + 4.0 * g_xzz_zz_x_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_zz_x_zz[i] = -2.0 * g_x_zz_x_zz[i] * a_exp + 4.0 * g_xzz_zz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (960-966)

    #pragma omp simd aligned(g_x_zz_y_xx, g_x_zz_y_xy, g_x_zz_y_xz, g_x_zz_y_yy, g_x_zz_y_yz, g_x_zz_y_zz, g_xz_0_0_0_z_zz_y_xx, g_xz_0_0_0_z_zz_y_xy, g_xz_0_0_0_z_zz_y_xz, g_xz_0_0_0_z_zz_y_yy, g_xz_0_0_0_z_zz_y_yz, g_xz_0_0_0_z_zz_y_zz, g_xzz_zz_y_xx, g_xzz_zz_y_xy, g_xzz_zz_y_xz, g_xzz_zz_y_yy, g_xzz_zz_y_yz, g_xzz_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_zz_y_xx[i] = -2.0 * g_x_zz_y_xx[i] * a_exp + 4.0 * g_xzz_zz_y_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_zz_y_xy[i] = -2.0 * g_x_zz_y_xy[i] * a_exp + 4.0 * g_xzz_zz_y_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_zz_y_xz[i] = -2.0 * g_x_zz_y_xz[i] * a_exp + 4.0 * g_xzz_zz_y_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_zz_y_yy[i] = -2.0 * g_x_zz_y_yy[i] * a_exp + 4.0 * g_xzz_zz_y_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_zz_y_yz[i] = -2.0 * g_x_zz_y_yz[i] * a_exp + 4.0 * g_xzz_zz_y_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_zz_y_zz[i] = -2.0 * g_x_zz_y_zz[i] * a_exp + 4.0 * g_xzz_zz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (966-972)

    #pragma omp simd aligned(g_x_zz_z_xx, g_x_zz_z_xy, g_x_zz_z_xz, g_x_zz_z_yy, g_x_zz_z_yz, g_x_zz_z_zz, g_xz_0_0_0_z_zz_z_xx, g_xz_0_0_0_z_zz_z_xy, g_xz_0_0_0_z_zz_z_xz, g_xz_0_0_0_z_zz_z_yy, g_xz_0_0_0_z_zz_z_yz, g_xz_0_0_0_z_zz_z_zz, g_xzz_zz_z_xx, g_xzz_zz_z_xy, g_xzz_zz_z_xz, g_xzz_zz_z_yy, g_xzz_zz_z_yz, g_xzz_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_zz_z_xx[i] = -2.0 * g_x_zz_z_xx[i] * a_exp + 4.0 * g_xzz_zz_z_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_z_zz_z_xy[i] = -2.0 * g_x_zz_z_xy[i] * a_exp + 4.0 * g_xzz_zz_z_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_zz_z_xz[i] = -2.0 * g_x_zz_z_xz[i] * a_exp + 4.0 * g_xzz_zz_z_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_zz_z_yy[i] = -2.0 * g_x_zz_z_yy[i] * a_exp + 4.0 * g_xzz_zz_z_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_z_zz_z_yz[i] = -2.0 * g_x_zz_z_yz[i] * a_exp + 4.0 * g_xzz_zz_z_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_z_zz_z_zz[i] = -2.0 * g_x_zz_z_zz[i] * a_exp + 4.0 * g_xzz_zz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (972-978)

    #pragma omp simd aligned(g_x_xx_x_xx, g_x_xx_x_xy, g_x_xx_x_xz, g_x_xx_x_yy, g_x_xx_x_yz, g_x_xx_x_zz, g_xyy_xx_x_xx, g_xyy_xx_x_xy, g_xyy_xx_x_xz, g_xyy_xx_x_yy, g_xyy_xx_x_yz, g_xyy_xx_x_zz, g_yy_0_0_0_x_xx_x_xx, g_yy_0_0_0_x_xx_x_xy, g_yy_0_0_0_x_xx_x_xz, g_yy_0_0_0_x_xx_x_yy, g_yy_0_0_0_x_xx_x_yz, g_yy_0_0_0_x_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_xx_x_xx[i] = -2.0 * g_x_xx_x_xx[i] * a_exp + 4.0 * g_xyy_xx_x_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xx_x_xy[i] = -2.0 * g_x_xx_x_xy[i] * a_exp + 4.0 * g_xyy_xx_x_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xx_x_xz[i] = -2.0 * g_x_xx_x_xz[i] * a_exp + 4.0 * g_xyy_xx_x_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xx_x_yy[i] = -2.0 * g_x_xx_x_yy[i] * a_exp + 4.0 * g_xyy_xx_x_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xx_x_yz[i] = -2.0 * g_x_xx_x_yz[i] * a_exp + 4.0 * g_xyy_xx_x_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xx_x_zz[i] = -2.0 * g_x_xx_x_zz[i] * a_exp + 4.0 * g_xyy_xx_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (978-984)

    #pragma omp simd aligned(g_x_xx_y_xx, g_x_xx_y_xy, g_x_xx_y_xz, g_x_xx_y_yy, g_x_xx_y_yz, g_x_xx_y_zz, g_xyy_xx_y_xx, g_xyy_xx_y_xy, g_xyy_xx_y_xz, g_xyy_xx_y_yy, g_xyy_xx_y_yz, g_xyy_xx_y_zz, g_yy_0_0_0_x_xx_y_xx, g_yy_0_0_0_x_xx_y_xy, g_yy_0_0_0_x_xx_y_xz, g_yy_0_0_0_x_xx_y_yy, g_yy_0_0_0_x_xx_y_yz, g_yy_0_0_0_x_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_xx_y_xx[i] = -2.0 * g_x_xx_y_xx[i] * a_exp + 4.0 * g_xyy_xx_y_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xx_y_xy[i] = -2.0 * g_x_xx_y_xy[i] * a_exp + 4.0 * g_xyy_xx_y_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xx_y_xz[i] = -2.0 * g_x_xx_y_xz[i] * a_exp + 4.0 * g_xyy_xx_y_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xx_y_yy[i] = -2.0 * g_x_xx_y_yy[i] * a_exp + 4.0 * g_xyy_xx_y_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xx_y_yz[i] = -2.0 * g_x_xx_y_yz[i] * a_exp + 4.0 * g_xyy_xx_y_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xx_y_zz[i] = -2.0 * g_x_xx_y_zz[i] * a_exp + 4.0 * g_xyy_xx_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (984-990)

    #pragma omp simd aligned(g_x_xx_z_xx, g_x_xx_z_xy, g_x_xx_z_xz, g_x_xx_z_yy, g_x_xx_z_yz, g_x_xx_z_zz, g_xyy_xx_z_xx, g_xyy_xx_z_xy, g_xyy_xx_z_xz, g_xyy_xx_z_yy, g_xyy_xx_z_yz, g_xyy_xx_z_zz, g_yy_0_0_0_x_xx_z_xx, g_yy_0_0_0_x_xx_z_xy, g_yy_0_0_0_x_xx_z_xz, g_yy_0_0_0_x_xx_z_yy, g_yy_0_0_0_x_xx_z_yz, g_yy_0_0_0_x_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_xx_z_xx[i] = -2.0 * g_x_xx_z_xx[i] * a_exp + 4.0 * g_xyy_xx_z_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xx_z_xy[i] = -2.0 * g_x_xx_z_xy[i] * a_exp + 4.0 * g_xyy_xx_z_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xx_z_xz[i] = -2.0 * g_x_xx_z_xz[i] * a_exp + 4.0 * g_xyy_xx_z_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xx_z_yy[i] = -2.0 * g_x_xx_z_yy[i] * a_exp + 4.0 * g_xyy_xx_z_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xx_z_yz[i] = -2.0 * g_x_xx_z_yz[i] * a_exp + 4.0 * g_xyy_xx_z_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xx_z_zz[i] = -2.0 * g_x_xx_z_zz[i] * a_exp + 4.0 * g_xyy_xx_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (990-996)

    #pragma omp simd aligned(g_x_xy_x_xx, g_x_xy_x_xy, g_x_xy_x_xz, g_x_xy_x_yy, g_x_xy_x_yz, g_x_xy_x_zz, g_xyy_xy_x_xx, g_xyy_xy_x_xy, g_xyy_xy_x_xz, g_xyy_xy_x_yy, g_xyy_xy_x_yz, g_xyy_xy_x_zz, g_yy_0_0_0_x_xy_x_xx, g_yy_0_0_0_x_xy_x_xy, g_yy_0_0_0_x_xy_x_xz, g_yy_0_0_0_x_xy_x_yy, g_yy_0_0_0_x_xy_x_yz, g_yy_0_0_0_x_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_xy_x_xx[i] = -2.0 * g_x_xy_x_xx[i] * a_exp + 4.0 * g_xyy_xy_x_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xy_x_xy[i] = -2.0 * g_x_xy_x_xy[i] * a_exp + 4.0 * g_xyy_xy_x_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xy_x_xz[i] = -2.0 * g_x_xy_x_xz[i] * a_exp + 4.0 * g_xyy_xy_x_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xy_x_yy[i] = -2.0 * g_x_xy_x_yy[i] * a_exp + 4.0 * g_xyy_xy_x_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xy_x_yz[i] = -2.0 * g_x_xy_x_yz[i] * a_exp + 4.0 * g_xyy_xy_x_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xy_x_zz[i] = -2.0 * g_x_xy_x_zz[i] * a_exp + 4.0 * g_xyy_xy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (996-1002)

    #pragma omp simd aligned(g_x_xy_y_xx, g_x_xy_y_xy, g_x_xy_y_xz, g_x_xy_y_yy, g_x_xy_y_yz, g_x_xy_y_zz, g_xyy_xy_y_xx, g_xyy_xy_y_xy, g_xyy_xy_y_xz, g_xyy_xy_y_yy, g_xyy_xy_y_yz, g_xyy_xy_y_zz, g_yy_0_0_0_x_xy_y_xx, g_yy_0_0_0_x_xy_y_xy, g_yy_0_0_0_x_xy_y_xz, g_yy_0_0_0_x_xy_y_yy, g_yy_0_0_0_x_xy_y_yz, g_yy_0_0_0_x_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_xy_y_xx[i] = -2.0 * g_x_xy_y_xx[i] * a_exp + 4.0 * g_xyy_xy_y_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xy_y_xy[i] = -2.0 * g_x_xy_y_xy[i] * a_exp + 4.0 * g_xyy_xy_y_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xy_y_xz[i] = -2.0 * g_x_xy_y_xz[i] * a_exp + 4.0 * g_xyy_xy_y_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xy_y_yy[i] = -2.0 * g_x_xy_y_yy[i] * a_exp + 4.0 * g_xyy_xy_y_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xy_y_yz[i] = -2.0 * g_x_xy_y_yz[i] * a_exp + 4.0 * g_xyy_xy_y_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xy_y_zz[i] = -2.0 * g_x_xy_y_zz[i] * a_exp + 4.0 * g_xyy_xy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1002-1008)

    #pragma omp simd aligned(g_x_xy_z_xx, g_x_xy_z_xy, g_x_xy_z_xz, g_x_xy_z_yy, g_x_xy_z_yz, g_x_xy_z_zz, g_xyy_xy_z_xx, g_xyy_xy_z_xy, g_xyy_xy_z_xz, g_xyy_xy_z_yy, g_xyy_xy_z_yz, g_xyy_xy_z_zz, g_yy_0_0_0_x_xy_z_xx, g_yy_0_0_0_x_xy_z_xy, g_yy_0_0_0_x_xy_z_xz, g_yy_0_0_0_x_xy_z_yy, g_yy_0_0_0_x_xy_z_yz, g_yy_0_0_0_x_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_xy_z_xx[i] = -2.0 * g_x_xy_z_xx[i] * a_exp + 4.0 * g_xyy_xy_z_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xy_z_xy[i] = -2.0 * g_x_xy_z_xy[i] * a_exp + 4.0 * g_xyy_xy_z_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xy_z_xz[i] = -2.0 * g_x_xy_z_xz[i] * a_exp + 4.0 * g_xyy_xy_z_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xy_z_yy[i] = -2.0 * g_x_xy_z_yy[i] * a_exp + 4.0 * g_xyy_xy_z_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xy_z_yz[i] = -2.0 * g_x_xy_z_yz[i] * a_exp + 4.0 * g_xyy_xy_z_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xy_z_zz[i] = -2.0 * g_x_xy_z_zz[i] * a_exp + 4.0 * g_xyy_xy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1008-1014)

    #pragma omp simd aligned(g_x_xz_x_xx, g_x_xz_x_xy, g_x_xz_x_xz, g_x_xz_x_yy, g_x_xz_x_yz, g_x_xz_x_zz, g_xyy_xz_x_xx, g_xyy_xz_x_xy, g_xyy_xz_x_xz, g_xyy_xz_x_yy, g_xyy_xz_x_yz, g_xyy_xz_x_zz, g_yy_0_0_0_x_xz_x_xx, g_yy_0_0_0_x_xz_x_xy, g_yy_0_0_0_x_xz_x_xz, g_yy_0_0_0_x_xz_x_yy, g_yy_0_0_0_x_xz_x_yz, g_yy_0_0_0_x_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_xz_x_xx[i] = -2.0 * g_x_xz_x_xx[i] * a_exp + 4.0 * g_xyy_xz_x_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xz_x_xy[i] = -2.0 * g_x_xz_x_xy[i] * a_exp + 4.0 * g_xyy_xz_x_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xz_x_xz[i] = -2.0 * g_x_xz_x_xz[i] * a_exp + 4.0 * g_xyy_xz_x_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xz_x_yy[i] = -2.0 * g_x_xz_x_yy[i] * a_exp + 4.0 * g_xyy_xz_x_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xz_x_yz[i] = -2.0 * g_x_xz_x_yz[i] * a_exp + 4.0 * g_xyy_xz_x_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xz_x_zz[i] = -2.0 * g_x_xz_x_zz[i] * a_exp + 4.0 * g_xyy_xz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1014-1020)

    #pragma omp simd aligned(g_x_xz_y_xx, g_x_xz_y_xy, g_x_xz_y_xz, g_x_xz_y_yy, g_x_xz_y_yz, g_x_xz_y_zz, g_xyy_xz_y_xx, g_xyy_xz_y_xy, g_xyy_xz_y_xz, g_xyy_xz_y_yy, g_xyy_xz_y_yz, g_xyy_xz_y_zz, g_yy_0_0_0_x_xz_y_xx, g_yy_0_0_0_x_xz_y_xy, g_yy_0_0_0_x_xz_y_xz, g_yy_0_0_0_x_xz_y_yy, g_yy_0_0_0_x_xz_y_yz, g_yy_0_0_0_x_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_xz_y_xx[i] = -2.0 * g_x_xz_y_xx[i] * a_exp + 4.0 * g_xyy_xz_y_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xz_y_xy[i] = -2.0 * g_x_xz_y_xy[i] * a_exp + 4.0 * g_xyy_xz_y_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xz_y_xz[i] = -2.0 * g_x_xz_y_xz[i] * a_exp + 4.0 * g_xyy_xz_y_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xz_y_yy[i] = -2.0 * g_x_xz_y_yy[i] * a_exp + 4.0 * g_xyy_xz_y_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xz_y_yz[i] = -2.0 * g_x_xz_y_yz[i] * a_exp + 4.0 * g_xyy_xz_y_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xz_y_zz[i] = -2.0 * g_x_xz_y_zz[i] * a_exp + 4.0 * g_xyy_xz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1020-1026)

    #pragma omp simd aligned(g_x_xz_z_xx, g_x_xz_z_xy, g_x_xz_z_xz, g_x_xz_z_yy, g_x_xz_z_yz, g_x_xz_z_zz, g_xyy_xz_z_xx, g_xyy_xz_z_xy, g_xyy_xz_z_xz, g_xyy_xz_z_yy, g_xyy_xz_z_yz, g_xyy_xz_z_zz, g_yy_0_0_0_x_xz_z_xx, g_yy_0_0_0_x_xz_z_xy, g_yy_0_0_0_x_xz_z_xz, g_yy_0_0_0_x_xz_z_yy, g_yy_0_0_0_x_xz_z_yz, g_yy_0_0_0_x_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_xz_z_xx[i] = -2.0 * g_x_xz_z_xx[i] * a_exp + 4.0 * g_xyy_xz_z_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xz_z_xy[i] = -2.0 * g_x_xz_z_xy[i] * a_exp + 4.0 * g_xyy_xz_z_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xz_z_xz[i] = -2.0 * g_x_xz_z_xz[i] * a_exp + 4.0 * g_xyy_xz_z_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xz_z_yy[i] = -2.0 * g_x_xz_z_yy[i] * a_exp + 4.0 * g_xyy_xz_z_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xz_z_yz[i] = -2.0 * g_x_xz_z_yz[i] * a_exp + 4.0 * g_xyy_xz_z_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xz_z_zz[i] = -2.0 * g_x_xz_z_zz[i] * a_exp + 4.0 * g_xyy_xz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1026-1032)

    #pragma omp simd aligned(g_x_yy_x_xx, g_x_yy_x_xy, g_x_yy_x_xz, g_x_yy_x_yy, g_x_yy_x_yz, g_x_yy_x_zz, g_xyy_yy_x_xx, g_xyy_yy_x_xy, g_xyy_yy_x_xz, g_xyy_yy_x_yy, g_xyy_yy_x_yz, g_xyy_yy_x_zz, g_yy_0_0_0_x_yy_x_xx, g_yy_0_0_0_x_yy_x_xy, g_yy_0_0_0_x_yy_x_xz, g_yy_0_0_0_x_yy_x_yy, g_yy_0_0_0_x_yy_x_yz, g_yy_0_0_0_x_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_yy_x_xx[i] = -2.0 * g_x_yy_x_xx[i] * a_exp + 4.0 * g_xyy_yy_x_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yy_x_xy[i] = -2.0 * g_x_yy_x_xy[i] * a_exp + 4.0 * g_xyy_yy_x_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yy_x_xz[i] = -2.0 * g_x_yy_x_xz[i] * a_exp + 4.0 * g_xyy_yy_x_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yy_x_yy[i] = -2.0 * g_x_yy_x_yy[i] * a_exp + 4.0 * g_xyy_yy_x_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yy_x_yz[i] = -2.0 * g_x_yy_x_yz[i] * a_exp + 4.0 * g_xyy_yy_x_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yy_x_zz[i] = -2.0 * g_x_yy_x_zz[i] * a_exp + 4.0 * g_xyy_yy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1032-1038)

    #pragma omp simd aligned(g_x_yy_y_xx, g_x_yy_y_xy, g_x_yy_y_xz, g_x_yy_y_yy, g_x_yy_y_yz, g_x_yy_y_zz, g_xyy_yy_y_xx, g_xyy_yy_y_xy, g_xyy_yy_y_xz, g_xyy_yy_y_yy, g_xyy_yy_y_yz, g_xyy_yy_y_zz, g_yy_0_0_0_x_yy_y_xx, g_yy_0_0_0_x_yy_y_xy, g_yy_0_0_0_x_yy_y_xz, g_yy_0_0_0_x_yy_y_yy, g_yy_0_0_0_x_yy_y_yz, g_yy_0_0_0_x_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_yy_y_xx[i] = -2.0 * g_x_yy_y_xx[i] * a_exp + 4.0 * g_xyy_yy_y_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yy_y_xy[i] = -2.0 * g_x_yy_y_xy[i] * a_exp + 4.0 * g_xyy_yy_y_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yy_y_xz[i] = -2.0 * g_x_yy_y_xz[i] * a_exp + 4.0 * g_xyy_yy_y_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yy_y_yy[i] = -2.0 * g_x_yy_y_yy[i] * a_exp + 4.0 * g_xyy_yy_y_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yy_y_yz[i] = -2.0 * g_x_yy_y_yz[i] * a_exp + 4.0 * g_xyy_yy_y_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yy_y_zz[i] = -2.0 * g_x_yy_y_zz[i] * a_exp + 4.0 * g_xyy_yy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1038-1044)

    #pragma omp simd aligned(g_x_yy_z_xx, g_x_yy_z_xy, g_x_yy_z_xz, g_x_yy_z_yy, g_x_yy_z_yz, g_x_yy_z_zz, g_xyy_yy_z_xx, g_xyy_yy_z_xy, g_xyy_yy_z_xz, g_xyy_yy_z_yy, g_xyy_yy_z_yz, g_xyy_yy_z_zz, g_yy_0_0_0_x_yy_z_xx, g_yy_0_0_0_x_yy_z_xy, g_yy_0_0_0_x_yy_z_xz, g_yy_0_0_0_x_yy_z_yy, g_yy_0_0_0_x_yy_z_yz, g_yy_0_0_0_x_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_yy_z_xx[i] = -2.0 * g_x_yy_z_xx[i] * a_exp + 4.0 * g_xyy_yy_z_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yy_z_xy[i] = -2.0 * g_x_yy_z_xy[i] * a_exp + 4.0 * g_xyy_yy_z_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yy_z_xz[i] = -2.0 * g_x_yy_z_xz[i] * a_exp + 4.0 * g_xyy_yy_z_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yy_z_yy[i] = -2.0 * g_x_yy_z_yy[i] * a_exp + 4.0 * g_xyy_yy_z_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yy_z_yz[i] = -2.0 * g_x_yy_z_yz[i] * a_exp + 4.0 * g_xyy_yy_z_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yy_z_zz[i] = -2.0 * g_x_yy_z_zz[i] * a_exp + 4.0 * g_xyy_yy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1044-1050)

    #pragma omp simd aligned(g_x_yz_x_xx, g_x_yz_x_xy, g_x_yz_x_xz, g_x_yz_x_yy, g_x_yz_x_yz, g_x_yz_x_zz, g_xyy_yz_x_xx, g_xyy_yz_x_xy, g_xyy_yz_x_xz, g_xyy_yz_x_yy, g_xyy_yz_x_yz, g_xyy_yz_x_zz, g_yy_0_0_0_x_yz_x_xx, g_yy_0_0_0_x_yz_x_xy, g_yy_0_0_0_x_yz_x_xz, g_yy_0_0_0_x_yz_x_yy, g_yy_0_0_0_x_yz_x_yz, g_yy_0_0_0_x_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_yz_x_xx[i] = -2.0 * g_x_yz_x_xx[i] * a_exp + 4.0 * g_xyy_yz_x_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yz_x_xy[i] = -2.0 * g_x_yz_x_xy[i] * a_exp + 4.0 * g_xyy_yz_x_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yz_x_xz[i] = -2.0 * g_x_yz_x_xz[i] * a_exp + 4.0 * g_xyy_yz_x_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yz_x_yy[i] = -2.0 * g_x_yz_x_yy[i] * a_exp + 4.0 * g_xyy_yz_x_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yz_x_yz[i] = -2.0 * g_x_yz_x_yz[i] * a_exp + 4.0 * g_xyy_yz_x_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yz_x_zz[i] = -2.0 * g_x_yz_x_zz[i] * a_exp + 4.0 * g_xyy_yz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1050-1056)

    #pragma omp simd aligned(g_x_yz_y_xx, g_x_yz_y_xy, g_x_yz_y_xz, g_x_yz_y_yy, g_x_yz_y_yz, g_x_yz_y_zz, g_xyy_yz_y_xx, g_xyy_yz_y_xy, g_xyy_yz_y_xz, g_xyy_yz_y_yy, g_xyy_yz_y_yz, g_xyy_yz_y_zz, g_yy_0_0_0_x_yz_y_xx, g_yy_0_0_0_x_yz_y_xy, g_yy_0_0_0_x_yz_y_xz, g_yy_0_0_0_x_yz_y_yy, g_yy_0_0_0_x_yz_y_yz, g_yy_0_0_0_x_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_yz_y_xx[i] = -2.0 * g_x_yz_y_xx[i] * a_exp + 4.0 * g_xyy_yz_y_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yz_y_xy[i] = -2.0 * g_x_yz_y_xy[i] * a_exp + 4.0 * g_xyy_yz_y_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yz_y_xz[i] = -2.0 * g_x_yz_y_xz[i] * a_exp + 4.0 * g_xyy_yz_y_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yz_y_yy[i] = -2.0 * g_x_yz_y_yy[i] * a_exp + 4.0 * g_xyy_yz_y_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yz_y_yz[i] = -2.0 * g_x_yz_y_yz[i] * a_exp + 4.0 * g_xyy_yz_y_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yz_y_zz[i] = -2.0 * g_x_yz_y_zz[i] * a_exp + 4.0 * g_xyy_yz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1056-1062)

    #pragma omp simd aligned(g_x_yz_z_xx, g_x_yz_z_xy, g_x_yz_z_xz, g_x_yz_z_yy, g_x_yz_z_yz, g_x_yz_z_zz, g_xyy_yz_z_xx, g_xyy_yz_z_xy, g_xyy_yz_z_xz, g_xyy_yz_z_yy, g_xyy_yz_z_yz, g_xyy_yz_z_zz, g_yy_0_0_0_x_yz_z_xx, g_yy_0_0_0_x_yz_z_xy, g_yy_0_0_0_x_yz_z_xz, g_yy_0_0_0_x_yz_z_yy, g_yy_0_0_0_x_yz_z_yz, g_yy_0_0_0_x_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_yz_z_xx[i] = -2.0 * g_x_yz_z_xx[i] * a_exp + 4.0 * g_xyy_yz_z_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yz_z_xy[i] = -2.0 * g_x_yz_z_xy[i] * a_exp + 4.0 * g_xyy_yz_z_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yz_z_xz[i] = -2.0 * g_x_yz_z_xz[i] * a_exp + 4.0 * g_xyy_yz_z_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yz_z_yy[i] = -2.0 * g_x_yz_z_yy[i] * a_exp + 4.0 * g_xyy_yz_z_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yz_z_yz[i] = -2.0 * g_x_yz_z_yz[i] * a_exp + 4.0 * g_xyy_yz_z_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yz_z_zz[i] = -2.0 * g_x_yz_z_zz[i] * a_exp + 4.0 * g_xyy_yz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1062-1068)

    #pragma omp simd aligned(g_x_zz_x_xx, g_x_zz_x_xy, g_x_zz_x_xz, g_x_zz_x_yy, g_x_zz_x_yz, g_x_zz_x_zz, g_xyy_zz_x_xx, g_xyy_zz_x_xy, g_xyy_zz_x_xz, g_xyy_zz_x_yy, g_xyy_zz_x_yz, g_xyy_zz_x_zz, g_yy_0_0_0_x_zz_x_xx, g_yy_0_0_0_x_zz_x_xy, g_yy_0_0_0_x_zz_x_xz, g_yy_0_0_0_x_zz_x_yy, g_yy_0_0_0_x_zz_x_yz, g_yy_0_0_0_x_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_zz_x_xx[i] = -2.0 * g_x_zz_x_xx[i] * a_exp + 4.0 * g_xyy_zz_x_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_zz_x_xy[i] = -2.0 * g_x_zz_x_xy[i] * a_exp + 4.0 * g_xyy_zz_x_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_zz_x_xz[i] = -2.0 * g_x_zz_x_xz[i] * a_exp + 4.0 * g_xyy_zz_x_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_zz_x_yy[i] = -2.0 * g_x_zz_x_yy[i] * a_exp + 4.0 * g_xyy_zz_x_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_zz_x_yz[i] = -2.0 * g_x_zz_x_yz[i] * a_exp + 4.0 * g_xyy_zz_x_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_zz_x_zz[i] = -2.0 * g_x_zz_x_zz[i] * a_exp + 4.0 * g_xyy_zz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1068-1074)

    #pragma omp simd aligned(g_x_zz_y_xx, g_x_zz_y_xy, g_x_zz_y_xz, g_x_zz_y_yy, g_x_zz_y_yz, g_x_zz_y_zz, g_xyy_zz_y_xx, g_xyy_zz_y_xy, g_xyy_zz_y_xz, g_xyy_zz_y_yy, g_xyy_zz_y_yz, g_xyy_zz_y_zz, g_yy_0_0_0_x_zz_y_xx, g_yy_0_0_0_x_zz_y_xy, g_yy_0_0_0_x_zz_y_xz, g_yy_0_0_0_x_zz_y_yy, g_yy_0_0_0_x_zz_y_yz, g_yy_0_0_0_x_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_zz_y_xx[i] = -2.0 * g_x_zz_y_xx[i] * a_exp + 4.0 * g_xyy_zz_y_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_zz_y_xy[i] = -2.0 * g_x_zz_y_xy[i] * a_exp + 4.0 * g_xyy_zz_y_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_zz_y_xz[i] = -2.0 * g_x_zz_y_xz[i] * a_exp + 4.0 * g_xyy_zz_y_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_zz_y_yy[i] = -2.0 * g_x_zz_y_yy[i] * a_exp + 4.0 * g_xyy_zz_y_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_zz_y_yz[i] = -2.0 * g_x_zz_y_yz[i] * a_exp + 4.0 * g_xyy_zz_y_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_zz_y_zz[i] = -2.0 * g_x_zz_y_zz[i] * a_exp + 4.0 * g_xyy_zz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1074-1080)

    #pragma omp simd aligned(g_x_zz_z_xx, g_x_zz_z_xy, g_x_zz_z_xz, g_x_zz_z_yy, g_x_zz_z_yz, g_x_zz_z_zz, g_xyy_zz_z_xx, g_xyy_zz_z_xy, g_xyy_zz_z_xz, g_xyy_zz_z_yy, g_xyy_zz_z_yz, g_xyy_zz_z_zz, g_yy_0_0_0_x_zz_z_xx, g_yy_0_0_0_x_zz_z_xy, g_yy_0_0_0_x_zz_z_xz, g_yy_0_0_0_x_zz_z_yy, g_yy_0_0_0_x_zz_z_yz, g_yy_0_0_0_x_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_zz_z_xx[i] = -2.0 * g_x_zz_z_xx[i] * a_exp + 4.0 * g_xyy_zz_z_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_x_zz_z_xy[i] = -2.0 * g_x_zz_z_xy[i] * a_exp + 4.0 * g_xyy_zz_z_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_zz_z_xz[i] = -2.0 * g_x_zz_z_xz[i] * a_exp + 4.0 * g_xyy_zz_z_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_zz_z_yy[i] = -2.0 * g_x_zz_z_yy[i] * a_exp + 4.0 * g_xyy_zz_z_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_x_zz_z_yz[i] = -2.0 * g_x_zz_z_yz[i] * a_exp + 4.0 * g_xyy_zz_z_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_x_zz_z_zz[i] = -2.0 * g_x_zz_z_zz[i] * a_exp + 4.0 * g_xyy_zz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1080-1086)

    #pragma omp simd aligned(g_y_xx_x_xx, g_y_xx_x_xy, g_y_xx_x_xz, g_y_xx_x_yy, g_y_xx_x_yz, g_y_xx_x_zz, g_yy_0_0_0_y_xx_x_xx, g_yy_0_0_0_y_xx_x_xy, g_yy_0_0_0_y_xx_x_xz, g_yy_0_0_0_y_xx_x_yy, g_yy_0_0_0_y_xx_x_yz, g_yy_0_0_0_y_xx_x_zz, g_yyy_xx_x_xx, g_yyy_xx_x_xy, g_yyy_xx_x_xz, g_yyy_xx_x_yy, g_yyy_xx_x_yz, g_yyy_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_xx_x_xx[i] = -6.0 * g_y_xx_x_xx[i] * a_exp + 4.0 * g_yyy_xx_x_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xx_x_xy[i] = -6.0 * g_y_xx_x_xy[i] * a_exp + 4.0 * g_yyy_xx_x_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xx_x_xz[i] = -6.0 * g_y_xx_x_xz[i] * a_exp + 4.0 * g_yyy_xx_x_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xx_x_yy[i] = -6.0 * g_y_xx_x_yy[i] * a_exp + 4.0 * g_yyy_xx_x_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xx_x_yz[i] = -6.0 * g_y_xx_x_yz[i] * a_exp + 4.0 * g_yyy_xx_x_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xx_x_zz[i] = -6.0 * g_y_xx_x_zz[i] * a_exp + 4.0 * g_yyy_xx_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1086-1092)

    #pragma omp simd aligned(g_y_xx_y_xx, g_y_xx_y_xy, g_y_xx_y_xz, g_y_xx_y_yy, g_y_xx_y_yz, g_y_xx_y_zz, g_yy_0_0_0_y_xx_y_xx, g_yy_0_0_0_y_xx_y_xy, g_yy_0_0_0_y_xx_y_xz, g_yy_0_0_0_y_xx_y_yy, g_yy_0_0_0_y_xx_y_yz, g_yy_0_0_0_y_xx_y_zz, g_yyy_xx_y_xx, g_yyy_xx_y_xy, g_yyy_xx_y_xz, g_yyy_xx_y_yy, g_yyy_xx_y_yz, g_yyy_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_xx_y_xx[i] = -6.0 * g_y_xx_y_xx[i] * a_exp + 4.0 * g_yyy_xx_y_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xx_y_xy[i] = -6.0 * g_y_xx_y_xy[i] * a_exp + 4.0 * g_yyy_xx_y_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xx_y_xz[i] = -6.0 * g_y_xx_y_xz[i] * a_exp + 4.0 * g_yyy_xx_y_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xx_y_yy[i] = -6.0 * g_y_xx_y_yy[i] * a_exp + 4.0 * g_yyy_xx_y_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xx_y_yz[i] = -6.0 * g_y_xx_y_yz[i] * a_exp + 4.0 * g_yyy_xx_y_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xx_y_zz[i] = -6.0 * g_y_xx_y_zz[i] * a_exp + 4.0 * g_yyy_xx_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1092-1098)

    #pragma omp simd aligned(g_y_xx_z_xx, g_y_xx_z_xy, g_y_xx_z_xz, g_y_xx_z_yy, g_y_xx_z_yz, g_y_xx_z_zz, g_yy_0_0_0_y_xx_z_xx, g_yy_0_0_0_y_xx_z_xy, g_yy_0_0_0_y_xx_z_xz, g_yy_0_0_0_y_xx_z_yy, g_yy_0_0_0_y_xx_z_yz, g_yy_0_0_0_y_xx_z_zz, g_yyy_xx_z_xx, g_yyy_xx_z_xy, g_yyy_xx_z_xz, g_yyy_xx_z_yy, g_yyy_xx_z_yz, g_yyy_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_xx_z_xx[i] = -6.0 * g_y_xx_z_xx[i] * a_exp + 4.0 * g_yyy_xx_z_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xx_z_xy[i] = -6.0 * g_y_xx_z_xy[i] * a_exp + 4.0 * g_yyy_xx_z_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xx_z_xz[i] = -6.0 * g_y_xx_z_xz[i] * a_exp + 4.0 * g_yyy_xx_z_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xx_z_yy[i] = -6.0 * g_y_xx_z_yy[i] * a_exp + 4.0 * g_yyy_xx_z_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xx_z_yz[i] = -6.0 * g_y_xx_z_yz[i] * a_exp + 4.0 * g_yyy_xx_z_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xx_z_zz[i] = -6.0 * g_y_xx_z_zz[i] * a_exp + 4.0 * g_yyy_xx_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1098-1104)

    #pragma omp simd aligned(g_y_xy_x_xx, g_y_xy_x_xy, g_y_xy_x_xz, g_y_xy_x_yy, g_y_xy_x_yz, g_y_xy_x_zz, g_yy_0_0_0_y_xy_x_xx, g_yy_0_0_0_y_xy_x_xy, g_yy_0_0_0_y_xy_x_xz, g_yy_0_0_0_y_xy_x_yy, g_yy_0_0_0_y_xy_x_yz, g_yy_0_0_0_y_xy_x_zz, g_yyy_xy_x_xx, g_yyy_xy_x_xy, g_yyy_xy_x_xz, g_yyy_xy_x_yy, g_yyy_xy_x_yz, g_yyy_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_xy_x_xx[i] = -6.0 * g_y_xy_x_xx[i] * a_exp + 4.0 * g_yyy_xy_x_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xy_x_xy[i] = -6.0 * g_y_xy_x_xy[i] * a_exp + 4.0 * g_yyy_xy_x_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xy_x_xz[i] = -6.0 * g_y_xy_x_xz[i] * a_exp + 4.0 * g_yyy_xy_x_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xy_x_yy[i] = -6.0 * g_y_xy_x_yy[i] * a_exp + 4.0 * g_yyy_xy_x_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xy_x_yz[i] = -6.0 * g_y_xy_x_yz[i] * a_exp + 4.0 * g_yyy_xy_x_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xy_x_zz[i] = -6.0 * g_y_xy_x_zz[i] * a_exp + 4.0 * g_yyy_xy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1104-1110)

    #pragma omp simd aligned(g_y_xy_y_xx, g_y_xy_y_xy, g_y_xy_y_xz, g_y_xy_y_yy, g_y_xy_y_yz, g_y_xy_y_zz, g_yy_0_0_0_y_xy_y_xx, g_yy_0_0_0_y_xy_y_xy, g_yy_0_0_0_y_xy_y_xz, g_yy_0_0_0_y_xy_y_yy, g_yy_0_0_0_y_xy_y_yz, g_yy_0_0_0_y_xy_y_zz, g_yyy_xy_y_xx, g_yyy_xy_y_xy, g_yyy_xy_y_xz, g_yyy_xy_y_yy, g_yyy_xy_y_yz, g_yyy_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_xy_y_xx[i] = -6.0 * g_y_xy_y_xx[i] * a_exp + 4.0 * g_yyy_xy_y_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xy_y_xy[i] = -6.0 * g_y_xy_y_xy[i] * a_exp + 4.0 * g_yyy_xy_y_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xy_y_xz[i] = -6.0 * g_y_xy_y_xz[i] * a_exp + 4.0 * g_yyy_xy_y_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xy_y_yy[i] = -6.0 * g_y_xy_y_yy[i] * a_exp + 4.0 * g_yyy_xy_y_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xy_y_yz[i] = -6.0 * g_y_xy_y_yz[i] * a_exp + 4.0 * g_yyy_xy_y_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xy_y_zz[i] = -6.0 * g_y_xy_y_zz[i] * a_exp + 4.0 * g_yyy_xy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1110-1116)

    #pragma omp simd aligned(g_y_xy_z_xx, g_y_xy_z_xy, g_y_xy_z_xz, g_y_xy_z_yy, g_y_xy_z_yz, g_y_xy_z_zz, g_yy_0_0_0_y_xy_z_xx, g_yy_0_0_0_y_xy_z_xy, g_yy_0_0_0_y_xy_z_xz, g_yy_0_0_0_y_xy_z_yy, g_yy_0_0_0_y_xy_z_yz, g_yy_0_0_0_y_xy_z_zz, g_yyy_xy_z_xx, g_yyy_xy_z_xy, g_yyy_xy_z_xz, g_yyy_xy_z_yy, g_yyy_xy_z_yz, g_yyy_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_xy_z_xx[i] = -6.0 * g_y_xy_z_xx[i] * a_exp + 4.0 * g_yyy_xy_z_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xy_z_xy[i] = -6.0 * g_y_xy_z_xy[i] * a_exp + 4.0 * g_yyy_xy_z_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xy_z_xz[i] = -6.0 * g_y_xy_z_xz[i] * a_exp + 4.0 * g_yyy_xy_z_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xy_z_yy[i] = -6.0 * g_y_xy_z_yy[i] * a_exp + 4.0 * g_yyy_xy_z_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xy_z_yz[i] = -6.0 * g_y_xy_z_yz[i] * a_exp + 4.0 * g_yyy_xy_z_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xy_z_zz[i] = -6.0 * g_y_xy_z_zz[i] * a_exp + 4.0 * g_yyy_xy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1116-1122)

    #pragma omp simd aligned(g_y_xz_x_xx, g_y_xz_x_xy, g_y_xz_x_xz, g_y_xz_x_yy, g_y_xz_x_yz, g_y_xz_x_zz, g_yy_0_0_0_y_xz_x_xx, g_yy_0_0_0_y_xz_x_xy, g_yy_0_0_0_y_xz_x_xz, g_yy_0_0_0_y_xz_x_yy, g_yy_0_0_0_y_xz_x_yz, g_yy_0_0_0_y_xz_x_zz, g_yyy_xz_x_xx, g_yyy_xz_x_xy, g_yyy_xz_x_xz, g_yyy_xz_x_yy, g_yyy_xz_x_yz, g_yyy_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_xz_x_xx[i] = -6.0 * g_y_xz_x_xx[i] * a_exp + 4.0 * g_yyy_xz_x_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xz_x_xy[i] = -6.0 * g_y_xz_x_xy[i] * a_exp + 4.0 * g_yyy_xz_x_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xz_x_xz[i] = -6.0 * g_y_xz_x_xz[i] * a_exp + 4.0 * g_yyy_xz_x_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xz_x_yy[i] = -6.0 * g_y_xz_x_yy[i] * a_exp + 4.0 * g_yyy_xz_x_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xz_x_yz[i] = -6.0 * g_y_xz_x_yz[i] * a_exp + 4.0 * g_yyy_xz_x_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xz_x_zz[i] = -6.0 * g_y_xz_x_zz[i] * a_exp + 4.0 * g_yyy_xz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1122-1128)

    #pragma omp simd aligned(g_y_xz_y_xx, g_y_xz_y_xy, g_y_xz_y_xz, g_y_xz_y_yy, g_y_xz_y_yz, g_y_xz_y_zz, g_yy_0_0_0_y_xz_y_xx, g_yy_0_0_0_y_xz_y_xy, g_yy_0_0_0_y_xz_y_xz, g_yy_0_0_0_y_xz_y_yy, g_yy_0_0_0_y_xz_y_yz, g_yy_0_0_0_y_xz_y_zz, g_yyy_xz_y_xx, g_yyy_xz_y_xy, g_yyy_xz_y_xz, g_yyy_xz_y_yy, g_yyy_xz_y_yz, g_yyy_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_xz_y_xx[i] = -6.0 * g_y_xz_y_xx[i] * a_exp + 4.0 * g_yyy_xz_y_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xz_y_xy[i] = -6.0 * g_y_xz_y_xy[i] * a_exp + 4.0 * g_yyy_xz_y_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xz_y_xz[i] = -6.0 * g_y_xz_y_xz[i] * a_exp + 4.0 * g_yyy_xz_y_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xz_y_yy[i] = -6.0 * g_y_xz_y_yy[i] * a_exp + 4.0 * g_yyy_xz_y_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xz_y_yz[i] = -6.0 * g_y_xz_y_yz[i] * a_exp + 4.0 * g_yyy_xz_y_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xz_y_zz[i] = -6.0 * g_y_xz_y_zz[i] * a_exp + 4.0 * g_yyy_xz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1128-1134)

    #pragma omp simd aligned(g_y_xz_z_xx, g_y_xz_z_xy, g_y_xz_z_xz, g_y_xz_z_yy, g_y_xz_z_yz, g_y_xz_z_zz, g_yy_0_0_0_y_xz_z_xx, g_yy_0_0_0_y_xz_z_xy, g_yy_0_0_0_y_xz_z_xz, g_yy_0_0_0_y_xz_z_yy, g_yy_0_0_0_y_xz_z_yz, g_yy_0_0_0_y_xz_z_zz, g_yyy_xz_z_xx, g_yyy_xz_z_xy, g_yyy_xz_z_xz, g_yyy_xz_z_yy, g_yyy_xz_z_yz, g_yyy_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_xz_z_xx[i] = -6.0 * g_y_xz_z_xx[i] * a_exp + 4.0 * g_yyy_xz_z_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xz_z_xy[i] = -6.0 * g_y_xz_z_xy[i] * a_exp + 4.0 * g_yyy_xz_z_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xz_z_xz[i] = -6.0 * g_y_xz_z_xz[i] * a_exp + 4.0 * g_yyy_xz_z_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xz_z_yy[i] = -6.0 * g_y_xz_z_yy[i] * a_exp + 4.0 * g_yyy_xz_z_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xz_z_yz[i] = -6.0 * g_y_xz_z_yz[i] * a_exp + 4.0 * g_yyy_xz_z_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xz_z_zz[i] = -6.0 * g_y_xz_z_zz[i] * a_exp + 4.0 * g_yyy_xz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1134-1140)

    #pragma omp simd aligned(g_y_yy_x_xx, g_y_yy_x_xy, g_y_yy_x_xz, g_y_yy_x_yy, g_y_yy_x_yz, g_y_yy_x_zz, g_yy_0_0_0_y_yy_x_xx, g_yy_0_0_0_y_yy_x_xy, g_yy_0_0_0_y_yy_x_xz, g_yy_0_0_0_y_yy_x_yy, g_yy_0_0_0_y_yy_x_yz, g_yy_0_0_0_y_yy_x_zz, g_yyy_yy_x_xx, g_yyy_yy_x_xy, g_yyy_yy_x_xz, g_yyy_yy_x_yy, g_yyy_yy_x_yz, g_yyy_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_yy_x_xx[i] = -6.0 * g_y_yy_x_xx[i] * a_exp + 4.0 * g_yyy_yy_x_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yy_x_xy[i] = -6.0 * g_y_yy_x_xy[i] * a_exp + 4.0 * g_yyy_yy_x_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yy_x_xz[i] = -6.0 * g_y_yy_x_xz[i] * a_exp + 4.0 * g_yyy_yy_x_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yy_x_yy[i] = -6.0 * g_y_yy_x_yy[i] * a_exp + 4.0 * g_yyy_yy_x_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yy_x_yz[i] = -6.0 * g_y_yy_x_yz[i] * a_exp + 4.0 * g_yyy_yy_x_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yy_x_zz[i] = -6.0 * g_y_yy_x_zz[i] * a_exp + 4.0 * g_yyy_yy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1140-1146)

    #pragma omp simd aligned(g_y_yy_y_xx, g_y_yy_y_xy, g_y_yy_y_xz, g_y_yy_y_yy, g_y_yy_y_yz, g_y_yy_y_zz, g_yy_0_0_0_y_yy_y_xx, g_yy_0_0_0_y_yy_y_xy, g_yy_0_0_0_y_yy_y_xz, g_yy_0_0_0_y_yy_y_yy, g_yy_0_0_0_y_yy_y_yz, g_yy_0_0_0_y_yy_y_zz, g_yyy_yy_y_xx, g_yyy_yy_y_xy, g_yyy_yy_y_xz, g_yyy_yy_y_yy, g_yyy_yy_y_yz, g_yyy_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_yy_y_xx[i] = -6.0 * g_y_yy_y_xx[i] * a_exp + 4.0 * g_yyy_yy_y_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yy_y_xy[i] = -6.0 * g_y_yy_y_xy[i] * a_exp + 4.0 * g_yyy_yy_y_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yy_y_xz[i] = -6.0 * g_y_yy_y_xz[i] * a_exp + 4.0 * g_yyy_yy_y_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yy_y_yy[i] = -6.0 * g_y_yy_y_yy[i] * a_exp + 4.0 * g_yyy_yy_y_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yy_y_yz[i] = -6.0 * g_y_yy_y_yz[i] * a_exp + 4.0 * g_yyy_yy_y_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yy_y_zz[i] = -6.0 * g_y_yy_y_zz[i] * a_exp + 4.0 * g_yyy_yy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1146-1152)

    #pragma omp simd aligned(g_y_yy_z_xx, g_y_yy_z_xy, g_y_yy_z_xz, g_y_yy_z_yy, g_y_yy_z_yz, g_y_yy_z_zz, g_yy_0_0_0_y_yy_z_xx, g_yy_0_0_0_y_yy_z_xy, g_yy_0_0_0_y_yy_z_xz, g_yy_0_0_0_y_yy_z_yy, g_yy_0_0_0_y_yy_z_yz, g_yy_0_0_0_y_yy_z_zz, g_yyy_yy_z_xx, g_yyy_yy_z_xy, g_yyy_yy_z_xz, g_yyy_yy_z_yy, g_yyy_yy_z_yz, g_yyy_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_yy_z_xx[i] = -6.0 * g_y_yy_z_xx[i] * a_exp + 4.0 * g_yyy_yy_z_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yy_z_xy[i] = -6.0 * g_y_yy_z_xy[i] * a_exp + 4.0 * g_yyy_yy_z_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yy_z_xz[i] = -6.0 * g_y_yy_z_xz[i] * a_exp + 4.0 * g_yyy_yy_z_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yy_z_yy[i] = -6.0 * g_y_yy_z_yy[i] * a_exp + 4.0 * g_yyy_yy_z_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yy_z_yz[i] = -6.0 * g_y_yy_z_yz[i] * a_exp + 4.0 * g_yyy_yy_z_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yy_z_zz[i] = -6.0 * g_y_yy_z_zz[i] * a_exp + 4.0 * g_yyy_yy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1152-1158)

    #pragma omp simd aligned(g_y_yz_x_xx, g_y_yz_x_xy, g_y_yz_x_xz, g_y_yz_x_yy, g_y_yz_x_yz, g_y_yz_x_zz, g_yy_0_0_0_y_yz_x_xx, g_yy_0_0_0_y_yz_x_xy, g_yy_0_0_0_y_yz_x_xz, g_yy_0_0_0_y_yz_x_yy, g_yy_0_0_0_y_yz_x_yz, g_yy_0_0_0_y_yz_x_zz, g_yyy_yz_x_xx, g_yyy_yz_x_xy, g_yyy_yz_x_xz, g_yyy_yz_x_yy, g_yyy_yz_x_yz, g_yyy_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_yz_x_xx[i] = -6.0 * g_y_yz_x_xx[i] * a_exp + 4.0 * g_yyy_yz_x_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yz_x_xy[i] = -6.0 * g_y_yz_x_xy[i] * a_exp + 4.0 * g_yyy_yz_x_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yz_x_xz[i] = -6.0 * g_y_yz_x_xz[i] * a_exp + 4.0 * g_yyy_yz_x_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yz_x_yy[i] = -6.0 * g_y_yz_x_yy[i] * a_exp + 4.0 * g_yyy_yz_x_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yz_x_yz[i] = -6.0 * g_y_yz_x_yz[i] * a_exp + 4.0 * g_yyy_yz_x_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yz_x_zz[i] = -6.0 * g_y_yz_x_zz[i] * a_exp + 4.0 * g_yyy_yz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1158-1164)

    #pragma omp simd aligned(g_y_yz_y_xx, g_y_yz_y_xy, g_y_yz_y_xz, g_y_yz_y_yy, g_y_yz_y_yz, g_y_yz_y_zz, g_yy_0_0_0_y_yz_y_xx, g_yy_0_0_0_y_yz_y_xy, g_yy_0_0_0_y_yz_y_xz, g_yy_0_0_0_y_yz_y_yy, g_yy_0_0_0_y_yz_y_yz, g_yy_0_0_0_y_yz_y_zz, g_yyy_yz_y_xx, g_yyy_yz_y_xy, g_yyy_yz_y_xz, g_yyy_yz_y_yy, g_yyy_yz_y_yz, g_yyy_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_yz_y_xx[i] = -6.0 * g_y_yz_y_xx[i] * a_exp + 4.0 * g_yyy_yz_y_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yz_y_xy[i] = -6.0 * g_y_yz_y_xy[i] * a_exp + 4.0 * g_yyy_yz_y_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yz_y_xz[i] = -6.0 * g_y_yz_y_xz[i] * a_exp + 4.0 * g_yyy_yz_y_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yz_y_yy[i] = -6.0 * g_y_yz_y_yy[i] * a_exp + 4.0 * g_yyy_yz_y_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yz_y_yz[i] = -6.0 * g_y_yz_y_yz[i] * a_exp + 4.0 * g_yyy_yz_y_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yz_y_zz[i] = -6.0 * g_y_yz_y_zz[i] * a_exp + 4.0 * g_yyy_yz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1164-1170)

    #pragma omp simd aligned(g_y_yz_z_xx, g_y_yz_z_xy, g_y_yz_z_xz, g_y_yz_z_yy, g_y_yz_z_yz, g_y_yz_z_zz, g_yy_0_0_0_y_yz_z_xx, g_yy_0_0_0_y_yz_z_xy, g_yy_0_0_0_y_yz_z_xz, g_yy_0_0_0_y_yz_z_yy, g_yy_0_0_0_y_yz_z_yz, g_yy_0_0_0_y_yz_z_zz, g_yyy_yz_z_xx, g_yyy_yz_z_xy, g_yyy_yz_z_xz, g_yyy_yz_z_yy, g_yyy_yz_z_yz, g_yyy_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_yz_z_xx[i] = -6.0 * g_y_yz_z_xx[i] * a_exp + 4.0 * g_yyy_yz_z_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yz_z_xy[i] = -6.0 * g_y_yz_z_xy[i] * a_exp + 4.0 * g_yyy_yz_z_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yz_z_xz[i] = -6.0 * g_y_yz_z_xz[i] * a_exp + 4.0 * g_yyy_yz_z_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yz_z_yy[i] = -6.0 * g_y_yz_z_yy[i] * a_exp + 4.0 * g_yyy_yz_z_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yz_z_yz[i] = -6.0 * g_y_yz_z_yz[i] * a_exp + 4.0 * g_yyy_yz_z_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yz_z_zz[i] = -6.0 * g_y_yz_z_zz[i] * a_exp + 4.0 * g_yyy_yz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1170-1176)

    #pragma omp simd aligned(g_y_zz_x_xx, g_y_zz_x_xy, g_y_zz_x_xz, g_y_zz_x_yy, g_y_zz_x_yz, g_y_zz_x_zz, g_yy_0_0_0_y_zz_x_xx, g_yy_0_0_0_y_zz_x_xy, g_yy_0_0_0_y_zz_x_xz, g_yy_0_0_0_y_zz_x_yy, g_yy_0_0_0_y_zz_x_yz, g_yy_0_0_0_y_zz_x_zz, g_yyy_zz_x_xx, g_yyy_zz_x_xy, g_yyy_zz_x_xz, g_yyy_zz_x_yy, g_yyy_zz_x_yz, g_yyy_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_zz_x_xx[i] = -6.0 * g_y_zz_x_xx[i] * a_exp + 4.0 * g_yyy_zz_x_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_zz_x_xy[i] = -6.0 * g_y_zz_x_xy[i] * a_exp + 4.0 * g_yyy_zz_x_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_zz_x_xz[i] = -6.0 * g_y_zz_x_xz[i] * a_exp + 4.0 * g_yyy_zz_x_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_zz_x_yy[i] = -6.0 * g_y_zz_x_yy[i] * a_exp + 4.0 * g_yyy_zz_x_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_zz_x_yz[i] = -6.0 * g_y_zz_x_yz[i] * a_exp + 4.0 * g_yyy_zz_x_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_zz_x_zz[i] = -6.0 * g_y_zz_x_zz[i] * a_exp + 4.0 * g_yyy_zz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1176-1182)

    #pragma omp simd aligned(g_y_zz_y_xx, g_y_zz_y_xy, g_y_zz_y_xz, g_y_zz_y_yy, g_y_zz_y_yz, g_y_zz_y_zz, g_yy_0_0_0_y_zz_y_xx, g_yy_0_0_0_y_zz_y_xy, g_yy_0_0_0_y_zz_y_xz, g_yy_0_0_0_y_zz_y_yy, g_yy_0_0_0_y_zz_y_yz, g_yy_0_0_0_y_zz_y_zz, g_yyy_zz_y_xx, g_yyy_zz_y_xy, g_yyy_zz_y_xz, g_yyy_zz_y_yy, g_yyy_zz_y_yz, g_yyy_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_zz_y_xx[i] = -6.0 * g_y_zz_y_xx[i] * a_exp + 4.0 * g_yyy_zz_y_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_zz_y_xy[i] = -6.0 * g_y_zz_y_xy[i] * a_exp + 4.0 * g_yyy_zz_y_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_zz_y_xz[i] = -6.0 * g_y_zz_y_xz[i] * a_exp + 4.0 * g_yyy_zz_y_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_zz_y_yy[i] = -6.0 * g_y_zz_y_yy[i] * a_exp + 4.0 * g_yyy_zz_y_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_zz_y_yz[i] = -6.0 * g_y_zz_y_yz[i] * a_exp + 4.0 * g_yyy_zz_y_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_zz_y_zz[i] = -6.0 * g_y_zz_y_zz[i] * a_exp + 4.0 * g_yyy_zz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1182-1188)

    #pragma omp simd aligned(g_y_zz_z_xx, g_y_zz_z_xy, g_y_zz_z_xz, g_y_zz_z_yy, g_y_zz_z_yz, g_y_zz_z_zz, g_yy_0_0_0_y_zz_z_xx, g_yy_0_0_0_y_zz_z_xy, g_yy_0_0_0_y_zz_z_xz, g_yy_0_0_0_y_zz_z_yy, g_yy_0_0_0_y_zz_z_yz, g_yy_0_0_0_y_zz_z_zz, g_yyy_zz_z_xx, g_yyy_zz_z_xy, g_yyy_zz_z_xz, g_yyy_zz_z_yy, g_yyy_zz_z_yz, g_yyy_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_zz_z_xx[i] = -6.0 * g_y_zz_z_xx[i] * a_exp + 4.0 * g_yyy_zz_z_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_y_zz_z_xy[i] = -6.0 * g_y_zz_z_xy[i] * a_exp + 4.0 * g_yyy_zz_z_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_zz_z_xz[i] = -6.0 * g_y_zz_z_xz[i] * a_exp + 4.0 * g_yyy_zz_z_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_zz_z_yy[i] = -6.0 * g_y_zz_z_yy[i] * a_exp + 4.0 * g_yyy_zz_z_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_y_zz_z_yz[i] = -6.0 * g_y_zz_z_yz[i] * a_exp + 4.0 * g_yyy_zz_z_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_y_zz_z_zz[i] = -6.0 * g_y_zz_z_zz[i] * a_exp + 4.0 * g_yyy_zz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1188-1194)

    #pragma omp simd aligned(g_yy_0_0_0_z_xx_x_xx, g_yy_0_0_0_z_xx_x_xy, g_yy_0_0_0_z_xx_x_xz, g_yy_0_0_0_z_xx_x_yy, g_yy_0_0_0_z_xx_x_yz, g_yy_0_0_0_z_xx_x_zz, g_yyz_xx_x_xx, g_yyz_xx_x_xy, g_yyz_xx_x_xz, g_yyz_xx_x_yy, g_yyz_xx_x_yz, g_yyz_xx_x_zz, g_z_xx_x_xx, g_z_xx_x_xy, g_z_xx_x_xz, g_z_xx_x_yy, g_z_xx_x_yz, g_z_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_xx_x_xx[i] = -2.0 * g_z_xx_x_xx[i] * a_exp + 4.0 * g_yyz_xx_x_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xx_x_xy[i] = -2.0 * g_z_xx_x_xy[i] * a_exp + 4.0 * g_yyz_xx_x_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xx_x_xz[i] = -2.0 * g_z_xx_x_xz[i] * a_exp + 4.0 * g_yyz_xx_x_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xx_x_yy[i] = -2.0 * g_z_xx_x_yy[i] * a_exp + 4.0 * g_yyz_xx_x_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xx_x_yz[i] = -2.0 * g_z_xx_x_yz[i] * a_exp + 4.0 * g_yyz_xx_x_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xx_x_zz[i] = -2.0 * g_z_xx_x_zz[i] * a_exp + 4.0 * g_yyz_xx_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1194-1200)

    #pragma omp simd aligned(g_yy_0_0_0_z_xx_y_xx, g_yy_0_0_0_z_xx_y_xy, g_yy_0_0_0_z_xx_y_xz, g_yy_0_0_0_z_xx_y_yy, g_yy_0_0_0_z_xx_y_yz, g_yy_0_0_0_z_xx_y_zz, g_yyz_xx_y_xx, g_yyz_xx_y_xy, g_yyz_xx_y_xz, g_yyz_xx_y_yy, g_yyz_xx_y_yz, g_yyz_xx_y_zz, g_z_xx_y_xx, g_z_xx_y_xy, g_z_xx_y_xz, g_z_xx_y_yy, g_z_xx_y_yz, g_z_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_xx_y_xx[i] = -2.0 * g_z_xx_y_xx[i] * a_exp + 4.0 * g_yyz_xx_y_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xx_y_xy[i] = -2.0 * g_z_xx_y_xy[i] * a_exp + 4.0 * g_yyz_xx_y_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xx_y_xz[i] = -2.0 * g_z_xx_y_xz[i] * a_exp + 4.0 * g_yyz_xx_y_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xx_y_yy[i] = -2.0 * g_z_xx_y_yy[i] * a_exp + 4.0 * g_yyz_xx_y_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xx_y_yz[i] = -2.0 * g_z_xx_y_yz[i] * a_exp + 4.0 * g_yyz_xx_y_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xx_y_zz[i] = -2.0 * g_z_xx_y_zz[i] * a_exp + 4.0 * g_yyz_xx_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1200-1206)

    #pragma omp simd aligned(g_yy_0_0_0_z_xx_z_xx, g_yy_0_0_0_z_xx_z_xy, g_yy_0_0_0_z_xx_z_xz, g_yy_0_0_0_z_xx_z_yy, g_yy_0_0_0_z_xx_z_yz, g_yy_0_0_0_z_xx_z_zz, g_yyz_xx_z_xx, g_yyz_xx_z_xy, g_yyz_xx_z_xz, g_yyz_xx_z_yy, g_yyz_xx_z_yz, g_yyz_xx_z_zz, g_z_xx_z_xx, g_z_xx_z_xy, g_z_xx_z_xz, g_z_xx_z_yy, g_z_xx_z_yz, g_z_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_xx_z_xx[i] = -2.0 * g_z_xx_z_xx[i] * a_exp + 4.0 * g_yyz_xx_z_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xx_z_xy[i] = -2.0 * g_z_xx_z_xy[i] * a_exp + 4.0 * g_yyz_xx_z_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xx_z_xz[i] = -2.0 * g_z_xx_z_xz[i] * a_exp + 4.0 * g_yyz_xx_z_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xx_z_yy[i] = -2.0 * g_z_xx_z_yy[i] * a_exp + 4.0 * g_yyz_xx_z_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xx_z_yz[i] = -2.0 * g_z_xx_z_yz[i] * a_exp + 4.0 * g_yyz_xx_z_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xx_z_zz[i] = -2.0 * g_z_xx_z_zz[i] * a_exp + 4.0 * g_yyz_xx_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1206-1212)

    #pragma omp simd aligned(g_yy_0_0_0_z_xy_x_xx, g_yy_0_0_0_z_xy_x_xy, g_yy_0_0_0_z_xy_x_xz, g_yy_0_0_0_z_xy_x_yy, g_yy_0_0_0_z_xy_x_yz, g_yy_0_0_0_z_xy_x_zz, g_yyz_xy_x_xx, g_yyz_xy_x_xy, g_yyz_xy_x_xz, g_yyz_xy_x_yy, g_yyz_xy_x_yz, g_yyz_xy_x_zz, g_z_xy_x_xx, g_z_xy_x_xy, g_z_xy_x_xz, g_z_xy_x_yy, g_z_xy_x_yz, g_z_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_xy_x_xx[i] = -2.0 * g_z_xy_x_xx[i] * a_exp + 4.0 * g_yyz_xy_x_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xy_x_xy[i] = -2.0 * g_z_xy_x_xy[i] * a_exp + 4.0 * g_yyz_xy_x_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xy_x_xz[i] = -2.0 * g_z_xy_x_xz[i] * a_exp + 4.0 * g_yyz_xy_x_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xy_x_yy[i] = -2.0 * g_z_xy_x_yy[i] * a_exp + 4.0 * g_yyz_xy_x_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xy_x_yz[i] = -2.0 * g_z_xy_x_yz[i] * a_exp + 4.0 * g_yyz_xy_x_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xy_x_zz[i] = -2.0 * g_z_xy_x_zz[i] * a_exp + 4.0 * g_yyz_xy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1212-1218)

    #pragma omp simd aligned(g_yy_0_0_0_z_xy_y_xx, g_yy_0_0_0_z_xy_y_xy, g_yy_0_0_0_z_xy_y_xz, g_yy_0_0_0_z_xy_y_yy, g_yy_0_0_0_z_xy_y_yz, g_yy_0_0_0_z_xy_y_zz, g_yyz_xy_y_xx, g_yyz_xy_y_xy, g_yyz_xy_y_xz, g_yyz_xy_y_yy, g_yyz_xy_y_yz, g_yyz_xy_y_zz, g_z_xy_y_xx, g_z_xy_y_xy, g_z_xy_y_xz, g_z_xy_y_yy, g_z_xy_y_yz, g_z_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_xy_y_xx[i] = -2.0 * g_z_xy_y_xx[i] * a_exp + 4.0 * g_yyz_xy_y_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xy_y_xy[i] = -2.0 * g_z_xy_y_xy[i] * a_exp + 4.0 * g_yyz_xy_y_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xy_y_xz[i] = -2.0 * g_z_xy_y_xz[i] * a_exp + 4.0 * g_yyz_xy_y_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xy_y_yy[i] = -2.0 * g_z_xy_y_yy[i] * a_exp + 4.0 * g_yyz_xy_y_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xy_y_yz[i] = -2.0 * g_z_xy_y_yz[i] * a_exp + 4.0 * g_yyz_xy_y_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xy_y_zz[i] = -2.0 * g_z_xy_y_zz[i] * a_exp + 4.0 * g_yyz_xy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1218-1224)

    #pragma omp simd aligned(g_yy_0_0_0_z_xy_z_xx, g_yy_0_0_0_z_xy_z_xy, g_yy_0_0_0_z_xy_z_xz, g_yy_0_0_0_z_xy_z_yy, g_yy_0_0_0_z_xy_z_yz, g_yy_0_0_0_z_xy_z_zz, g_yyz_xy_z_xx, g_yyz_xy_z_xy, g_yyz_xy_z_xz, g_yyz_xy_z_yy, g_yyz_xy_z_yz, g_yyz_xy_z_zz, g_z_xy_z_xx, g_z_xy_z_xy, g_z_xy_z_xz, g_z_xy_z_yy, g_z_xy_z_yz, g_z_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_xy_z_xx[i] = -2.0 * g_z_xy_z_xx[i] * a_exp + 4.0 * g_yyz_xy_z_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xy_z_xy[i] = -2.0 * g_z_xy_z_xy[i] * a_exp + 4.0 * g_yyz_xy_z_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xy_z_xz[i] = -2.0 * g_z_xy_z_xz[i] * a_exp + 4.0 * g_yyz_xy_z_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xy_z_yy[i] = -2.0 * g_z_xy_z_yy[i] * a_exp + 4.0 * g_yyz_xy_z_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xy_z_yz[i] = -2.0 * g_z_xy_z_yz[i] * a_exp + 4.0 * g_yyz_xy_z_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xy_z_zz[i] = -2.0 * g_z_xy_z_zz[i] * a_exp + 4.0 * g_yyz_xy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1224-1230)

    #pragma omp simd aligned(g_yy_0_0_0_z_xz_x_xx, g_yy_0_0_0_z_xz_x_xy, g_yy_0_0_0_z_xz_x_xz, g_yy_0_0_0_z_xz_x_yy, g_yy_0_0_0_z_xz_x_yz, g_yy_0_0_0_z_xz_x_zz, g_yyz_xz_x_xx, g_yyz_xz_x_xy, g_yyz_xz_x_xz, g_yyz_xz_x_yy, g_yyz_xz_x_yz, g_yyz_xz_x_zz, g_z_xz_x_xx, g_z_xz_x_xy, g_z_xz_x_xz, g_z_xz_x_yy, g_z_xz_x_yz, g_z_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_xz_x_xx[i] = -2.0 * g_z_xz_x_xx[i] * a_exp + 4.0 * g_yyz_xz_x_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xz_x_xy[i] = -2.0 * g_z_xz_x_xy[i] * a_exp + 4.0 * g_yyz_xz_x_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xz_x_xz[i] = -2.0 * g_z_xz_x_xz[i] * a_exp + 4.0 * g_yyz_xz_x_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xz_x_yy[i] = -2.0 * g_z_xz_x_yy[i] * a_exp + 4.0 * g_yyz_xz_x_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xz_x_yz[i] = -2.0 * g_z_xz_x_yz[i] * a_exp + 4.0 * g_yyz_xz_x_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xz_x_zz[i] = -2.0 * g_z_xz_x_zz[i] * a_exp + 4.0 * g_yyz_xz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1230-1236)

    #pragma omp simd aligned(g_yy_0_0_0_z_xz_y_xx, g_yy_0_0_0_z_xz_y_xy, g_yy_0_0_0_z_xz_y_xz, g_yy_0_0_0_z_xz_y_yy, g_yy_0_0_0_z_xz_y_yz, g_yy_0_0_0_z_xz_y_zz, g_yyz_xz_y_xx, g_yyz_xz_y_xy, g_yyz_xz_y_xz, g_yyz_xz_y_yy, g_yyz_xz_y_yz, g_yyz_xz_y_zz, g_z_xz_y_xx, g_z_xz_y_xy, g_z_xz_y_xz, g_z_xz_y_yy, g_z_xz_y_yz, g_z_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_xz_y_xx[i] = -2.0 * g_z_xz_y_xx[i] * a_exp + 4.0 * g_yyz_xz_y_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xz_y_xy[i] = -2.0 * g_z_xz_y_xy[i] * a_exp + 4.0 * g_yyz_xz_y_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xz_y_xz[i] = -2.0 * g_z_xz_y_xz[i] * a_exp + 4.0 * g_yyz_xz_y_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xz_y_yy[i] = -2.0 * g_z_xz_y_yy[i] * a_exp + 4.0 * g_yyz_xz_y_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xz_y_yz[i] = -2.0 * g_z_xz_y_yz[i] * a_exp + 4.0 * g_yyz_xz_y_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xz_y_zz[i] = -2.0 * g_z_xz_y_zz[i] * a_exp + 4.0 * g_yyz_xz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1236-1242)

    #pragma omp simd aligned(g_yy_0_0_0_z_xz_z_xx, g_yy_0_0_0_z_xz_z_xy, g_yy_0_0_0_z_xz_z_xz, g_yy_0_0_0_z_xz_z_yy, g_yy_0_0_0_z_xz_z_yz, g_yy_0_0_0_z_xz_z_zz, g_yyz_xz_z_xx, g_yyz_xz_z_xy, g_yyz_xz_z_xz, g_yyz_xz_z_yy, g_yyz_xz_z_yz, g_yyz_xz_z_zz, g_z_xz_z_xx, g_z_xz_z_xy, g_z_xz_z_xz, g_z_xz_z_yy, g_z_xz_z_yz, g_z_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_xz_z_xx[i] = -2.0 * g_z_xz_z_xx[i] * a_exp + 4.0 * g_yyz_xz_z_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xz_z_xy[i] = -2.0 * g_z_xz_z_xy[i] * a_exp + 4.0 * g_yyz_xz_z_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xz_z_xz[i] = -2.0 * g_z_xz_z_xz[i] * a_exp + 4.0 * g_yyz_xz_z_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xz_z_yy[i] = -2.0 * g_z_xz_z_yy[i] * a_exp + 4.0 * g_yyz_xz_z_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xz_z_yz[i] = -2.0 * g_z_xz_z_yz[i] * a_exp + 4.0 * g_yyz_xz_z_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xz_z_zz[i] = -2.0 * g_z_xz_z_zz[i] * a_exp + 4.0 * g_yyz_xz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1242-1248)

    #pragma omp simd aligned(g_yy_0_0_0_z_yy_x_xx, g_yy_0_0_0_z_yy_x_xy, g_yy_0_0_0_z_yy_x_xz, g_yy_0_0_0_z_yy_x_yy, g_yy_0_0_0_z_yy_x_yz, g_yy_0_0_0_z_yy_x_zz, g_yyz_yy_x_xx, g_yyz_yy_x_xy, g_yyz_yy_x_xz, g_yyz_yy_x_yy, g_yyz_yy_x_yz, g_yyz_yy_x_zz, g_z_yy_x_xx, g_z_yy_x_xy, g_z_yy_x_xz, g_z_yy_x_yy, g_z_yy_x_yz, g_z_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_yy_x_xx[i] = -2.0 * g_z_yy_x_xx[i] * a_exp + 4.0 * g_yyz_yy_x_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yy_x_xy[i] = -2.0 * g_z_yy_x_xy[i] * a_exp + 4.0 * g_yyz_yy_x_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yy_x_xz[i] = -2.0 * g_z_yy_x_xz[i] * a_exp + 4.0 * g_yyz_yy_x_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yy_x_yy[i] = -2.0 * g_z_yy_x_yy[i] * a_exp + 4.0 * g_yyz_yy_x_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yy_x_yz[i] = -2.0 * g_z_yy_x_yz[i] * a_exp + 4.0 * g_yyz_yy_x_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yy_x_zz[i] = -2.0 * g_z_yy_x_zz[i] * a_exp + 4.0 * g_yyz_yy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1248-1254)

    #pragma omp simd aligned(g_yy_0_0_0_z_yy_y_xx, g_yy_0_0_0_z_yy_y_xy, g_yy_0_0_0_z_yy_y_xz, g_yy_0_0_0_z_yy_y_yy, g_yy_0_0_0_z_yy_y_yz, g_yy_0_0_0_z_yy_y_zz, g_yyz_yy_y_xx, g_yyz_yy_y_xy, g_yyz_yy_y_xz, g_yyz_yy_y_yy, g_yyz_yy_y_yz, g_yyz_yy_y_zz, g_z_yy_y_xx, g_z_yy_y_xy, g_z_yy_y_xz, g_z_yy_y_yy, g_z_yy_y_yz, g_z_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_yy_y_xx[i] = -2.0 * g_z_yy_y_xx[i] * a_exp + 4.0 * g_yyz_yy_y_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yy_y_xy[i] = -2.0 * g_z_yy_y_xy[i] * a_exp + 4.0 * g_yyz_yy_y_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yy_y_xz[i] = -2.0 * g_z_yy_y_xz[i] * a_exp + 4.0 * g_yyz_yy_y_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yy_y_yy[i] = -2.0 * g_z_yy_y_yy[i] * a_exp + 4.0 * g_yyz_yy_y_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yy_y_yz[i] = -2.0 * g_z_yy_y_yz[i] * a_exp + 4.0 * g_yyz_yy_y_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yy_y_zz[i] = -2.0 * g_z_yy_y_zz[i] * a_exp + 4.0 * g_yyz_yy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1254-1260)

    #pragma omp simd aligned(g_yy_0_0_0_z_yy_z_xx, g_yy_0_0_0_z_yy_z_xy, g_yy_0_0_0_z_yy_z_xz, g_yy_0_0_0_z_yy_z_yy, g_yy_0_0_0_z_yy_z_yz, g_yy_0_0_0_z_yy_z_zz, g_yyz_yy_z_xx, g_yyz_yy_z_xy, g_yyz_yy_z_xz, g_yyz_yy_z_yy, g_yyz_yy_z_yz, g_yyz_yy_z_zz, g_z_yy_z_xx, g_z_yy_z_xy, g_z_yy_z_xz, g_z_yy_z_yy, g_z_yy_z_yz, g_z_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_yy_z_xx[i] = -2.0 * g_z_yy_z_xx[i] * a_exp + 4.0 * g_yyz_yy_z_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yy_z_xy[i] = -2.0 * g_z_yy_z_xy[i] * a_exp + 4.0 * g_yyz_yy_z_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yy_z_xz[i] = -2.0 * g_z_yy_z_xz[i] * a_exp + 4.0 * g_yyz_yy_z_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yy_z_yy[i] = -2.0 * g_z_yy_z_yy[i] * a_exp + 4.0 * g_yyz_yy_z_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yy_z_yz[i] = -2.0 * g_z_yy_z_yz[i] * a_exp + 4.0 * g_yyz_yy_z_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yy_z_zz[i] = -2.0 * g_z_yy_z_zz[i] * a_exp + 4.0 * g_yyz_yy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1260-1266)

    #pragma omp simd aligned(g_yy_0_0_0_z_yz_x_xx, g_yy_0_0_0_z_yz_x_xy, g_yy_0_0_0_z_yz_x_xz, g_yy_0_0_0_z_yz_x_yy, g_yy_0_0_0_z_yz_x_yz, g_yy_0_0_0_z_yz_x_zz, g_yyz_yz_x_xx, g_yyz_yz_x_xy, g_yyz_yz_x_xz, g_yyz_yz_x_yy, g_yyz_yz_x_yz, g_yyz_yz_x_zz, g_z_yz_x_xx, g_z_yz_x_xy, g_z_yz_x_xz, g_z_yz_x_yy, g_z_yz_x_yz, g_z_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_yz_x_xx[i] = -2.0 * g_z_yz_x_xx[i] * a_exp + 4.0 * g_yyz_yz_x_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yz_x_xy[i] = -2.0 * g_z_yz_x_xy[i] * a_exp + 4.0 * g_yyz_yz_x_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yz_x_xz[i] = -2.0 * g_z_yz_x_xz[i] * a_exp + 4.0 * g_yyz_yz_x_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yz_x_yy[i] = -2.0 * g_z_yz_x_yy[i] * a_exp + 4.0 * g_yyz_yz_x_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yz_x_yz[i] = -2.0 * g_z_yz_x_yz[i] * a_exp + 4.0 * g_yyz_yz_x_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yz_x_zz[i] = -2.0 * g_z_yz_x_zz[i] * a_exp + 4.0 * g_yyz_yz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1266-1272)

    #pragma omp simd aligned(g_yy_0_0_0_z_yz_y_xx, g_yy_0_0_0_z_yz_y_xy, g_yy_0_0_0_z_yz_y_xz, g_yy_0_0_0_z_yz_y_yy, g_yy_0_0_0_z_yz_y_yz, g_yy_0_0_0_z_yz_y_zz, g_yyz_yz_y_xx, g_yyz_yz_y_xy, g_yyz_yz_y_xz, g_yyz_yz_y_yy, g_yyz_yz_y_yz, g_yyz_yz_y_zz, g_z_yz_y_xx, g_z_yz_y_xy, g_z_yz_y_xz, g_z_yz_y_yy, g_z_yz_y_yz, g_z_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_yz_y_xx[i] = -2.0 * g_z_yz_y_xx[i] * a_exp + 4.0 * g_yyz_yz_y_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yz_y_xy[i] = -2.0 * g_z_yz_y_xy[i] * a_exp + 4.0 * g_yyz_yz_y_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yz_y_xz[i] = -2.0 * g_z_yz_y_xz[i] * a_exp + 4.0 * g_yyz_yz_y_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yz_y_yy[i] = -2.0 * g_z_yz_y_yy[i] * a_exp + 4.0 * g_yyz_yz_y_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yz_y_yz[i] = -2.0 * g_z_yz_y_yz[i] * a_exp + 4.0 * g_yyz_yz_y_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yz_y_zz[i] = -2.0 * g_z_yz_y_zz[i] * a_exp + 4.0 * g_yyz_yz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1272-1278)

    #pragma omp simd aligned(g_yy_0_0_0_z_yz_z_xx, g_yy_0_0_0_z_yz_z_xy, g_yy_0_0_0_z_yz_z_xz, g_yy_0_0_0_z_yz_z_yy, g_yy_0_0_0_z_yz_z_yz, g_yy_0_0_0_z_yz_z_zz, g_yyz_yz_z_xx, g_yyz_yz_z_xy, g_yyz_yz_z_xz, g_yyz_yz_z_yy, g_yyz_yz_z_yz, g_yyz_yz_z_zz, g_z_yz_z_xx, g_z_yz_z_xy, g_z_yz_z_xz, g_z_yz_z_yy, g_z_yz_z_yz, g_z_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_yz_z_xx[i] = -2.0 * g_z_yz_z_xx[i] * a_exp + 4.0 * g_yyz_yz_z_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yz_z_xy[i] = -2.0 * g_z_yz_z_xy[i] * a_exp + 4.0 * g_yyz_yz_z_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yz_z_xz[i] = -2.0 * g_z_yz_z_xz[i] * a_exp + 4.0 * g_yyz_yz_z_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yz_z_yy[i] = -2.0 * g_z_yz_z_yy[i] * a_exp + 4.0 * g_yyz_yz_z_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yz_z_yz[i] = -2.0 * g_z_yz_z_yz[i] * a_exp + 4.0 * g_yyz_yz_z_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yz_z_zz[i] = -2.0 * g_z_yz_z_zz[i] * a_exp + 4.0 * g_yyz_yz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1278-1284)

    #pragma omp simd aligned(g_yy_0_0_0_z_zz_x_xx, g_yy_0_0_0_z_zz_x_xy, g_yy_0_0_0_z_zz_x_xz, g_yy_0_0_0_z_zz_x_yy, g_yy_0_0_0_z_zz_x_yz, g_yy_0_0_0_z_zz_x_zz, g_yyz_zz_x_xx, g_yyz_zz_x_xy, g_yyz_zz_x_xz, g_yyz_zz_x_yy, g_yyz_zz_x_yz, g_yyz_zz_x_zz, g_z_zz_x_xx, g_z_zz_x_xy, g_z_zz_x_xz, g_z_zz_x_yy, g_z_zz_x_yz, g_z_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_zz_x_xx[i] = -2.0 * g_z_zz_x_xx[i] * a_exp + 4.0 * g_yyz_zz_x_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_zz_x_xy[i] = -2.0 * g_z_zz_x_xy[i] * a_exp + 4.0 * g_yyz_zz_x_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_zz_x_xz[i] = -2.0 * g_z_zz_x_xz[i] * a_exp + 4.0 * g_yyz_zz_x_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_zz_x_yy[i] = -2.0 * g_z_zz_x_yy[i] * a_exp + 4.0 * g_yyz_zz_x_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_zz_x_yz[i] = -2.0 * g_z_zz_x_yz[i] * a_exp + 4.0 * g_yyz_zz_x_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_zz_x_zz[i] = -2.0 * g_z_zz_x_zz[i] * a_exp + 4.0 * g_yyz_zz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1284-1290)

    #pragma omp simd aligned(g_yy_0_0_0_z_zz_y_xx, g_yy_0_0_0_z_zz_y_xy, g_yy_0_0_0_z_zz_y_xz, g_yy_0_0_0_z_zz_y_yy, g_yy_0_0_0_z_zz_y_yz, g_yy_0_0_0_z_zz_y_zz, g_yyz_zz_y_xx, g_yyz_zz_y_xy, g_yyz_zz_y_xz, g_yyz_zz_y_yy, g_yyz_zz_y_yz, g_yyz_zz_y_zz, g_z_zz_y_xx, g_z_zz_y_xy, g_z_zz_y_xz, g_z_zz_y_yy, g_z_zz_y_yz, g_z_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_zz_y_xx[i] = -2.0 * g_z_zz_y_xx[i] * a_exp + 4.0 * g_yyz_zz_y_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_zz_y_xy[i] = -2.0 * g_z_zz_y_xy[i] * a_exp + 4.0 * g_yyz_zz_y_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_zz_y_xz[i] = -2.0 * g_z_zz_y_xz[i] * a_exp + 4.0 * g_yyz_zz_y_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_zz_y_yy[i] = -2.0 * g_z_zz_y_yy[i] * a_exp + 4.0 * g_yyz_zz_y_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_zz_y_yz[i] = -2.0 * g_z_zz_y_yz[i] * a_exp + 4.0 * g_yyz_zz_y_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_zz_y_zz[i] = -2.0 * g_z_zz_y_zz[i] * a_exp + 4.0 * g_yyz_zz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1290-1296)

    #pragma omp simd aligned(g_yy_0_0_0_z_zz_z_xx, g_yy_0_0_0_z_zz_z_xy, g_yy_0_0_0_z_zz_z_xz, g_yy_0_0_0_z_zz_z_yy, g_yy_0_0_0_z_zz_z_yz, g_yy_0_0_0_z_zz_z_zz, g_yyz_zz_z_xx, g_yyz_zz_z_xy, g_yyz_zz_z_xz, g_yyz_zz_z_yy, g_yyz_zz_z_yz, g_yyz_zz_z_zz, g_z_zz_z_xx, g_z_zz_z_xy, g_z_zz_z_xz, g_z_zz_z_yy, g_z_zz_z_yz, g_z_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_zz_z_xx[i] = -2.0 * g_z_zz_z_xx[i] * a_exp + 4.0 * g_yyz_zz_z_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_z_zz_z_xy[i] = -2.0 * g_z_zz_z_xy[i] * a_exp + 4.0 * g_yyz_zz_z_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_zz_z_xz[i] = -2.0 * g_z_zz_z_xz[i] * a_exp + 4.0 * g_yyz_zz_z_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_zz_z_yy[i] = -2.0 * g_z_zz_z_yy[i] * a_exp + 4.0 * g_yyz_zz_z_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_z_zz_z_yz[i] = -2.0 * g_z_zz_z_yz[i] * a_exp + 4.0 * g_yyz_zz_z_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_z_zz_z_zz[i] = -2.0 * g_z_zz_z_zz[i] * a_exp + 4.0 * g_yyz_zz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1296-1302)

    #pragma omp simd aligned(g_xyz_xx_x_xx, g_xyz_xx_x_xy, g_xyz_xx_x_xz, g_xyz_xx_x_yy, g_xyz_xx_x_yz, g_xyz_xx_x_zz, g_yz_0_0_0_x_xx_x_xx, g_yz_0_0_0_x_xx_x_xy, g_yz_0_0_0_x_xx_x_xz, g_yz_0_0_0_x_xx_x_yy, g_yz_0_0_0_x_xx_x_yz, g_yz_0_0_0_x_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_xx_x_xx[i] = 4.0 * g_xyz_xx_x_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xx_x_xy[i] = 4.0 * g_xyz_xx_x_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xx_x_xz[i] = 4.0 * g_xyz_xx_x_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xx_x_yy[i] = 4.0 * g_xyz_xx_x_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xx_x_yz[i] = 4.0 * g_xyz_xx_x_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xx_x_zz[i] = 4.0 * g_xyz_xx_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1302-1308)

    #pragma omp simd aligned(g_xyz_xx_y_xx, g_xyz_xx_y_xy, g_xyz_xx_y_xz, g_xyz_xx_y_yy, g_xyz_xx_y_yz, g_xyz_xx_y_zz, g_yz_0_0_0_x_xx_y_xx, g_yz_0_0_0_x_xx_y_xy, g_yz_0_0_0_x_xx_y_xz, g_yz_0_0_0_x_xx_y_yy, g_yz_0_0_0_x_xx_y_yz, g_yz_0_0_0_x_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_xx_y_xx[i] = 4.0 * g_xyz_xx_y_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xx_y_xy[i] = 4.0 * g_xyz_xx_y_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xx_y_xz[i] = 4.0 * g_xyz_xx_y_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xx_y_yy[i] = 4.0 * g_xyz_xx_y_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xx_y_yz[i] = 4.0 * g_xyz_xx_y_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xx_y_zz[i] = 4.0 * g_xyz_xx_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1308-1314)

    #pragma omp simd aligned(g_xyz_xx_z_xx, g_xyz_xx_z_xy, g_xyz_xx_z_xz, g_xyz_xx_z_yy, g_xyz_xx_z_yz, g_xyz_xx_z_zz, g_yz_0_0_0_x_xx_z_xx, g_yz_0_0_0_x_xx_z_xy, g_yz_0_0_0_x_xx_z_xz, g_yz_0_0_0_x_xx_z_yy, g_yz_0_0_0_x_xx_z_yz, g_yz_0_0_0_x_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_xx_z_xx[i] = 4.0 * g_xyz_xx_z_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xx_z_xy[i] = 4.0 * g_xyz_xx_z_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xx_z_xz[i] = 4.0 * g_xyz_xx_z_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xx_z_yy[i] = 4.0 * g_xyz_xx_z_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xx_z_yz[i] = 4.0 * g_xyz_xx_z_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xx_z_zz[i] = 4.0 * g_xyz_xx_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1314-1320)

    #pragma omp simd aligned(g_xyz_xy_x_xx, g_xyz_xy_x_xy, g_xyz_xy_x_xz, g_xyz_xy_x_yy, g_xyz_xy_x_yz, g_xyz_xy_x_zz, g_yz_0_0_0_x_xy_x_xx, g_yz_0_0_0_x_xy_x_xy, g_yz_0_0_0_x_xy_x_xz, g_yz_0_0_0_x_xy_x_yy, g_yz_0_0_0_x_xy_x_yz, g_yz_0_0_0_x_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_xy_x_xx[i] = 4.0 * g_xyz_xy_x_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xy_x_xy[i] = 4.0 * g_xyz_xy_x_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xy_x_xz[i] = 4.0 * g_xyz_xy_x_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xy_x_yy[i] = 4.0 * g_xyz_xy_x_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xy_x_yz[i] = 4.0 * g_xyz_xy_x_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xy_x_zz[i] = 4.0 * g_xyz_xy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1320-1326)

    #pragma omp simd aligned(g_xyz_xy_y_xx, g_xyz_xy_y_xy, g_xyz_xy_y_xz, g_xyz_xy_y_yy, g_xyz_xy_y_yz, g_xyz_xy_y_zz, g_yz_0_0_0_x_xy_y_xx, g_yz_0_0_0_x_xy_y_xy, g_yz_0_0_0_x_xy_y_xz, g_yz_0_0_0_x_xy_y_yy, g_yz_0_0_0_x_xy_y_yz, g_yz_0_0_0_x_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_xy_y_xx[i] = 4.0 * g_xyz_xy_y_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xy_y_xy[i] = 4.0 * g_xyz_xy_y_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xy_y_xz[i] = 4.0 * g_xyz_xy_y_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xy_y_yy[i] = 4.0 * g_xyz_xy_y_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xy_y_yz[i] = 4.0 * g_xyz_xy_y_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xy_y_zz[i] = 4.0 * g_xyz_xy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1326-1332)

    #pragma omp simd aligned(g_xyz_xy_z_xx, g_xyz_xy_z_xy, g_xyz_xy_z_xz, g_xyz_xy_z_yy, g_xyz_xy_z_yz, g_xyz_xy_z_zz, g_yz_0_0_0_x_xy_z_xx, g_yz_0_0_0_x_xy_z_xy, g_yz_0_0_0_x_xy_z_xz, g_yz_0_0_0_x_xy_z_yy, g_yz_0_0_0_x_xy_z_yz, g_yz_0_0_0_x_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_xy_z_xx[i] = 4.0 * g_xyz_xy_z_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xy_z_xy[i] = 4.0 * g_xyz_xy_z_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xy_z_xz[i] = 4.0 * g_xyz_xy_z_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xy_z_yy[i] = 4.0 * g_xyz_xy_z_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xy_z_yz[i] = 4.0 * g_xyz_xy_z_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xy_z_zz[i] = 4.0 * g_xyz_xy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1332-1338)

    #pragma omp simd aligned(g_xyz_xz_x_xx, g_xyz_xz_x_xy, g_xyz_xz_x_xz, g_xyz_xz_x_yy, g_xyz_xz_x_yz, g_xyz_xz_x_zz, g_yz_0_0_0_x_xz_x_xx, g_yz_0_0_0_x_xz_x_xy, g_yz_0_0_0_x_xz_x_xz, g_yz_0_0_0_x_xz_x_yy, g_yz_0_0_0_x_xz_x_yz, g_yz_0_0_0_x_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_xz_x_xx[i] = 4.0 * g_xyz_xz_x_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xz_x_xy[i] = 4.0 * g_xyz_xz_x_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xz_x_xz[i] = 4.0 * g_xyz_xz_x_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xz_x_yy[i] = 4.0 * g_xyz_xz_x_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xz_x_yz[i] = 4.0 * g_xyz_xz_x_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xz_x_zz[i] = 4.0 * g_xyz_xz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1338-1344)

    #pragma omp simd aligned(g_xyz_xz_y_xx, g_xyz_xz_y_xy, g_xyz_xz_y_xz, g_xyz_xz_y_yy, g_xyz_xz_y_yz, g_xyz_xz_y_zz, g_yz_0_0_0_x_xz_y_xx, g_yz_0_0_0_x_xz_y_xy, g_yz_0_0_0_x_xz_y_xz, g_yz_0_0_0_x_xz_y_yy, g_yz_0_0_0_x_xz_y_yz, g_yz_0_0_0_x_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_xz_y_xx[i] = 4.0 * g_xyz_xz_y_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xz_y_xy[i] = 4.0 * g_xyz_xz_y_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xz_y_xz[i] = 4.0 * g_xyz_xz_y_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xz_y_yy[i] = 4.0 * g_xyz_xz_y_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xz_y_yz[i] = 4.0 * g_xyz_xz_y_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xz_y_zz[i] = 4.0 * g_xyz_xz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1344-1350)

    #pragma omp simd aligned(g_xyz_xz_z_xx, g_xyz_xz_z_xy, g_xyz_xz_z_xz, g_xyz_xz_z_yy, g_xyz_xz_z_yz, g_xyz_xz_z_zz, g_yz_0_0_0_x_xz_z_xx, g_yz_0_0_0_x_xz_z_xy, g_yz_0_0_0_x_xz_z_xz, g_yz_0_0_0_x_xz_z_yy, g_yz_0_0_0_x_xz_z_yz, g_yz_0_0_0_x_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_xz_z_xx[i] = 4.0 * g_xyz_xz_z_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xz_z_xy[i] = 4.0 * g_xyz_xz_z_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xz_z_xz[i] = 4.0 * g_xyz_xz_z_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xz_z_yy[i] = 4.0 * g_xyz_xz_z_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xz_z_yz[i] = 4.0 * g_xyz_xz_z_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xz_z_zz[i] = 4.0 * g_xyz_xz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1350-1356)

    #pragma omp simd aligned(g_xyz_yy_x_xx, g_xyz_yy_x_xy, g_xyz_yy_x_xz, g_xyz_yy_x_yy, g_xyz_yy_x_yz, g_xyz_yy_x_zz, g_yz_0_0_0_x_yy_x_xx, g_yz_0_0_0_x_yy_x_xy, g_yz_0_0_0_x_yy_x_xz, g_yz_0_0_0_x_yy_x_yy, g_yz_0_0_0_x_yy_x_yz, g_yz_0_0_0_x_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_yy_x_xx[i] = 4.0 * g_xyz_yy_x_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yy_x_xy[i] = 4.0 * g_xyz_yy_x_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yy_x_xz[i] = 4.0 * g_xyz_yy_x_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yy_x_yy[i] = 4.0 * g_xyz_yy_x_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yy_x_yz[i] = 4.0 * g_xyz_yy_x_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yy_x_zz[i] = 4.0 * g_xyz_yy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1356-1362)

    #pragma omp simd aligned(g_xyz_yy_y_xx, g_xyz_yy_y_xy, g_xyz_yy_y_xz, g_xyz_yy_y_yy, g_xyz_yy_y_yz, g_xyz_yy_y_zz, g_yz_0_0_0_x_yy_y_xx, g_yz_0_0_0_x_yy_y_xy, g_yz_0_0_0_x_yy_y_xz, g_yz_0_0_0_x_yy_y_yy, g_yz_0_0_0_x_yy_y_yz, g_yz_0_0_0_x_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_yy_y_xx[i] = 4.0 * g_xyz_yy_y_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yy_y_xy[i] = 4.0 * g_xyz_yy_y_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yy_y_xz[i] = 4.0 * g_xyz_yy_y_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yy_y_yy[i] = 4.0 * g_xyz_yy_y_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yy_y_yz[i] = 4.0 * g_xyz_yy_y_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yy_y_zz[i] = 4.0 * g_xyz_yy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1362-1368)

    #pragma omp simd aligned(g_xyz_yy_z_xx, g_xyz_yy_z_xy, g_xyz_yy_z_xz, g_xyz_yy_z_yy, g_xyz_yy_z_yz, g_xyz_yy_z_zz, g_yz_0_0_0_x_yy_z_xx, g_yz_0_0_0_x_yy_z_xy, g_yz_0_0_0_x_yy_z_xz, g_yz_0_0_0_x_yy_z_yy, g_yz_0_0_0_x_yy_z_yz, g_yz_0_0_0_x_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_yy_z_xx[i] = 4.0 * g_xyz_yy_z_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yy_z_xy[i] = 4.0 * g_xyz_yy_z_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yy_z_xz[i] = 4.0 * g_xyz_yy_z_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yy_z_yy[i] = 4.0 * g_xyz_yy_z_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yy_z_yz[i] = 4.0 * g_xyz_yy_z_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yy_z_zz[i] = 4.0 * g_xyz_yy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1368-1374)

    #pragma omp simd aligned(g_xyz_yz_x_xx, g_xyz_yz_x_xy, g_xyz_yz_x_xz, g_xyz_yz_x_yy, g_xyz_yz_x_yz, g_xyz_yz_x_zz, g_yz_0_0_0_x_yz_x_xx, g_yz_0_0_0_x_yz_x_xy, g_yz_0_0_0_x_yz_x_xz, g_yz_0_0_0_x_yz_x_yy, g_yz_0_0_0_x_yz_x_yz, g_yz_0_0_0_x_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_yz_x_xx[i] = 4.0 * g_xyz_yz_x_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yz_x_xy[i] = 4.0 * g_xyz_yz_x_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yz_x_xz[i] = 4.0 * g_xyz_yz_x_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yz_x_yy[i] = 4.0 * g_xyz_yz_x_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yz_x_yz[i] = 4.0 * g_xyz_yz_x_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yz_x_zz[i] = 4.0 * g_xyz_yz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1374-1380)

    #pragma omp simd aligned(g_xyz_yz_y_xx, g_xyz_yz_y_xy, g_xyz_yz_y_xz, g_xyz_yz_y_yy, g_xyz_yz_y_yz, g_xyz_yz_y_zz, g_yz_0_0_0_x_yz_y_xx, g_yz_0_0_0_x_yz_y_xy, g_yz_0_0_0_x_yz_y_xz, g_yz_0_0_0_x_yz_y_yy, g_yz_0_0_0_x_yz_y_yz, g_yz_0_0_0_x_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_yz_y_xx[i] = 4.0 * g_xyz_yz_y_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yz_y_xy[i] = 4.0 * g_xyz_yz_y_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yz_y_xz[i] = 4.0 * g_xyz_yz_y_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yz_y_yy[i] = 4.0 * g_xyz_yz_y_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yz_y_yz[i] = 4.0 * g_xyz_yz_y_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yz_y_zz[i] = 4.0 * g_xyz_yz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1380-1386)

    #pragma omp simd aligned(g_xyz_yz_z_xx, g_xyz_yz_z_xy, g_xyz_yz_z_xz, g_xyz_yz_z_yy, g_xyz_yz_z_yz, g_xyz_yz_z_zz, g_yz_0_0_0_x_yz_z_xx, g_yz_0_0_0_x_yz_z_xy, g_yz_0_0_0_x_yz_z_xz, g_yz_0_0_0_x_yz_z_yy, g_yz_0_0_0_x_yz_z_yz, g_yz_0_0_0_x_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_yz_z_xx[i] = 4.0 * g_xyz_yz_z_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yz_z_xy[i] = 4.0 * g_xyz_yz_z_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yz_z_xz[i] = 4.0 * g_xyz_yz_z_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yz_z_yy[i] = 4.0 * g_xyz_yz_z_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yz_z_yz[i] = 4.0 * g_xyz_yz_z_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yz_z_zz[i] = 4.0 * g_xyz_yz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1386-1392)

    #pragma omp simd aligned(g_xyz_zz_x_xx, g_xyz_zz_x_xy, g_xyz_zz_x_xz, g_xyz_zz_x_yy, g_xyz_zz_x_yz, g_xyz_zz_x_zz, g_yz_0_0_0_x_zz_x_xx, g_yz_0_0_0_x_zz_x_xy, g_yz_0_0_0_x_zz_x_xz, g_yz_0_0_0_x_zz_x_yy, g_yz_0_0_0_x_zz_x_yz, g_yz_0_0_0_x_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_zz_x_xx[i] = 4.0 * g_xyz_zz_x_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_zz_x_xy[i] = 4.0 * g_xyz_zz_x_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_zz_x_xz[i] = 4.0 * g_xyz_zz_x_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_zz_x_yy[i] = 4.0 * g_xyz_zz_x_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_zz_x_yz[i] = 4.0 * g_xyz_zz_x_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_zz_x_zz[i] = 4.0 * g_xyz_zz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1392-1398)

    #pragma omp simd aligned(g_xyz_zz_y_xx, g_xyz_zz_y_xy, g_xyz_zz_y_xz, g_xyz_zz_y_yy, g_xyz_zz_y_yz, g_xyz_zz_y_zz, g_yz_0_0_0_x_zz_y_xx, g_yz_0_0_0_x_zz_y_xy, g_yz_0_0_0_x_zz_y_xz, g_yz_0_0_0_x_zz_y_yy, g_yz_0_0_0_x_zz_y_yz, g_yz_0_0_0_x_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_zz_y_xx[i] = 4.0 * g_xyz_zz_y_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_zz_y_xy[i] = 4.0 * g_xyz_zz_y_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_zz_y_xz[i] = 4.0 * g_xyz_zz_y_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_zz_y_yy[i] = 4.0 * g_xyz_zz_y_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_zz_y_yz[i] = 4.0 * g_xyz_zz_y_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_zz_y_zz[i] = 4.0 * g_xyz_zz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1398-1404)

    #pragma omp simd aligned(g_xyz_zz_z_xx, g_xyz_zz_z_xy, g_xyz_zz_z_xz, g_xyz_zz_z_yy, g_xyz_zz_z_yz, g_xyz_zz_z_zz, g_yz_0_0_0_x_zz_z_xx, g_yz_0_0_0_x_zz_z_xy, g_yz_0_0_0_x_zz_z_xz, g_yz_0_0_0_x_zz_z_yy, g_yz_0_0_0_x_zz_z_yz, g_yz_0_0_0_x_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_zz_z_xx[i] = 4.0 * g_xyz_zz_z_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_x_zz_z_xy[i] = 4.0 * g_xyz_zz_z_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_zz_z_xz[i] = 4.0 * g_xyz_zz_z_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_zz_z_yy[i] = 4.0 * g_xyz_zz_z_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_x_zz_z_yz[i] = 4.0 * g_xyz_zz_z_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_x_zz_z_zz[i] = 4.0 * g_xyz_zz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1404-1410)

    #pragma omp simd aligned(g_yyz_xx_x_xx, g_yyz_xx_x_xy, g_yyz_xx_x_xz, g_yyz_xx_x_yy, g_yyz_xx_x_yz, g_yyz_xx_x_zz, g_yz_0_0_0_y_xx_x_xx, g_yz_0_0_0_y_xx_x_xy, g_yz_0_0_0_y_xx_x_xz, g_yz_0_0_0_y_xx_x_yy, g_yz_0_0_0_y_xx_x_yz, g_yz_0_0_0_y_xx_x_zz, g_z_xx_x_xx, g_z_xx_x_xy, g_z_xx_x_xz, g_z_xx_x_yy, g_z_xx_x_yz, g_z_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_xx_x_xx[i] = -2.0 * g_z_xx_x_xx[i] * a_exp + 4.0 * g_yyz_xx_x_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xx_x_xy[i] = -2.0 * g_z_xx_x_xy[i] * a_exp + 4.0 * g_yyz_xx_x_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xx_x_xz[i] = -2.0 * g_z_xx_x_xz[i] * a_exp + 4.0 * g_yyz_xx_x_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xx_x_yy[i] = -2.0 * g_z_xx_x_yy[i] * a_exp + 4.0 * g_yyz_xx_x_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xx_x_yz[i] = -2.0 * g_z_xx_x_yz[i] * a_exp + 4.0 * g_yyz_xx_x_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xx_x_zz[i] = -2.0 * g_z_xx_x_zz[i] * a_exp + 4.0 * g_yyz_xx_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1410-1416)

    #pragma omp simd aligned(g_yyz_xx_y_xx, g_yyz_xx_y_xy, g_yyz_xx_y_xz, g_yyz_xx_y_yy, g_yyz_xx_y_yz, g_yyz_xx_y_zz, g_yz_0_0_0_y_xx_y_xx, g_yz_0_0_0_y_xx_y_xy, g_yz_0_0_0_y_xx_y_xz, g_yz_0_0_0_y_xx_y_yy, g_yz_0_0_0_y_xx_y_yz, g_yz_0_0_0_y_xx_y_zz, g_z_xx_y_xx, g_z_xx_y_xy, g_z_xx_y_xz, g_z_xx_y_yy, g_z_xx_y_yz, g_z_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_xx_y_xx[i] = -2.0 * g_z_xx_y_xx[i] * a_exp + 4.0 * g_yyz_xx_y_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xx_y_xy[i] = -2.0 * g_z_xx_y_xy[i] * a_exp + 4.0 * g_yyz_xx_y_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xx_y_xz[i] = -2.0 * g_z_xx_y_xz[i] * a_exp + 4.0 * g_yyz_xx_y_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xx_y_yy[i] = -2.0 * g_z_xx_y_yy[i] * a_exp + 4.0 * g_yyz_xx_y_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xx_y_yz[i] = -2.0 * g_z_xx_y_yz[i] * a_exp + 4.0 * g_yyz_xx_y_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xx_y_zz[i] = -2.0 * g_z_xx_y_zz[i] * a_exp + 4.0 * g_yyz_xx_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1416-1422)

    #pragma omp simd aligned(g_yyz_xx_z_xx, g_yyz_xx_z_xy, g_yyz_xx_z_xz, g_yyz_xx_z_yy, g_yyz_xx_z_yz, g_yyz_xx_z_zz, g_yz_0_0_0_y_xx_z_xx, g_yz_0_0_0_y_xx_z_xy, g_yz_0_0_0_y_xx_z_xz, g_yz_0_0_0_y_xx_z_yy, g_yz_0_0_0_y_xx_z_yz, g_yz_0_0_0_y_xx_z_zz, g_z_xx_z_xx, g_z_xx_z_xy, g_z_xx_z_xz, g_z_xx_z_yy, g_z_xx_z_yz, g_z_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_xx_z_xx[i] = -2.0 * g_z_xx_z_xx[i] * a_exp + 4.0 * g_yyz_xx_z_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xx_z_xy[i] = -2.0 * g_z_xx_z_xy[i] * a_exp + 4.0 * g_yyz_xx_z_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xx_z_xz[i] = -2.0 * g_z_xx_z_xz[i] * a_exp + 4.0 * g_yyz_xx_z_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xx_z_yy[i] = -2.0 * g_z_xx_z_yy[i] * a_exp + 4.0 * g_yyz_xx_z_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xx_z_yz[i] = -2.0 * g_z_xx_z_yz[i] * a_exp + 4.0 * g_yyz_xx_z_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xx_z_zz[i] = -2.0 * g_z_xx_z_zz[i] * a_exp + 4.0 * g_yyz_xx_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1422-1428)

    #pragma omp simd aligned(g_yyz_xy_x_xx, g_yyz_xy_x_xy, g_yyz_xy_x_xz, g_yyz_xy_x_yy, g_yyz_xy_x_yz, g_yyz_xy_x_zz, g_yz_0_0_0_y_xy_x_xx, g_yz_0_0_0_y_xy_x_xy, g_yz_0_0_0_y_xy_x_xz, g_yz_0_0_0_y_xy_x_yy, g_yz_0_0_0_y_xy_x_yz, g_yz_0_0_0_y_xy_x_zz, g_z_xy_x_xx, g_z_xy_x_xy, g_z_xy_x_xz, g_z_xy_x_yy, g_z_xy_x_yz, g_z_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_xy_x_xx[i] = -2.0 * g_z_xy_x_xx[i] * a_exp + 4.0 * g_yyz_xy_x_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xy_x_xy[i] = -2.0 * g_z_xy_x_xy[i] * a_exp + 4.0 * g_yyz_xy_x_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xy_x_xz[i] = -2.0 * g_z_xy_x_xz[i] * a_exp + 4.0 * g_yyz_xy_x_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xy_x_yy[i] = -2.0 * g_z_xy_x_yy[i] * a_exp + 4.0 * g_yyz_xy_x_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xy_x_yz[i] = -2.0 * g_z_xy_x_yz[i] * a_exp + 4.0 * g_yyz_xy_x_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xy_x_zz[i] = -2.0 * g_z_xy_x_zz[i] * a_exp + 4.0 * g_yyz_xy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1428-1434)

    #pragma omp simd aligned(g_yyz_xy_y_xx, g_yyz_xy_y_xy, g_yyz_xy_y_xz, g_yyz_xy_y_yy, g_yyz_xy_y_yz, g_yyz_xy_y_zz, g_yz_0_0_0_y_xy_y_xx, g_yz_0_0_0_y_xy_y_xy, g_yz_0_0_0_y_xy_y_xz, g_yz_0_0_0_y_xy_y_yy, g_yz_0_0_0_y_xy_y_yz, g_yz_0_0_0_y_xy_y_zz, g_z_xy_y_xx, g_z_xy_y_xy, g_z_xy_y_xz, g_z_xy_y_yy, g_z_xy_y_yz, g_z_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_xy_y_xx[i] = -2.0 * g_z_xy_y_xx[i] * a_exp + 4.0 * g_yyz_xy_y_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xy_y_xy[i] = -2.0 * g_z_xy_y_xy[i] * a_exp + 4.0 * g_yyz_xy_y_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xy_y_xz[i] = -2.0 * g_z_xy_y_xz[i] * a_exp + 4.0 * g_yyz_xy_y_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xy_y_yy[i] = -2.0 * g_z_xy_y_yy[i] * a_exp + 4.0 * g_yyz_xy_y_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xy_y_yz[i] = -2.0 * g_z_xy_y_yz[i] * a_exp + 4.0 * g_yyz_xy_y_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xy_y_zz[i] = -2.0 * g_z_xy_y_zz[i] * a_exp + 4.0 * g_yyz_xy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1434-1440)

    #pragma omp simd aligned(g_yyz_xy_z_xx, g_yyz_xy_z_xy, g_yyz_xy_z_xz, g_yyz_xy_z_yy, g_yyz_xy_z_yz, g_yyz_xy_z_zz, g_yz_0_0_0_y_xy_z_xx, g_yz_0_0_0_y_xy_z_xy, g_yz_0_0_0_y_xy_z_xz, g_yz_0_0_0_y_xy_z_yy, g_yz_0_0_0_y_xy_z_yz, g_yz_0_0_0_y_xy_z_zz, g_z_xy_z_xx, g_z_xy_z_xy, g_z_xy_z_xz, g_z_xy_z_yy, g_z_xy_z_yz, g_z_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_xy_z_xx[i] = -2.0 * g_z_xy_z_xx[i] * a_exp + 4.0 * g_yyz_xy_z_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xy_z_xy[i] = -2.0 * g_z_xy_z_xy[i] * a_exp + 4.0 * g_yyz_xy_z_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xy_z_xz[i] = -2.0 * g_z_xy_z_xz[i] * a_exp + 4.0 * g_yyz_xy_z_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xy_z_yy[i] = -2.0 * g_z_xy_z_yy[i] * a_exp + 4.0 * g_yyz_xy_z_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xy_z_yz[i] = -2.0 * g_z_xy_z_yz[i] * a_exp + 4.0 * g_yyz_xy_z_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xy_z_zz[i] = -2.0 * g_z_xy_z_zz[i] * a_exp + 4.0 * g_yyz_xy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1440-1446)

    #pragma omp simd aligned(g_yyz_xz_x_xx, g_yyz_xz_x_xy, g_yyz_xz_x_xz, g_yyz_xz_x_yy, g_yyz_xz_x_yz, g_yyz_xz_x_zz, g_yz_0_0_0_y_xz_x_xx, g_yz_0_0_0_y_xz_x_xy, g_yz_0_0_0_y_xz_x_xz, g_yz_0_0_0_y_xz_x_yy, g_yz_0_0_0_y_xz_x_yz, g_yz_0_0_0_y_xz_x_zz, g_z_xz_x_xx, g_z_xz_x_xy, g_z_xz_x_xz, g_z_xz_x_yy, g_z_xz_x_yz, g_z_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_xz_x_xx[i] = -2.0 * g_z_xz_x_xx[i] * a_exp + 4.0 * g_yyz_xz_x_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xz_x_xy[i] = -2.0 * g_z_xz_x_xy[i] * a_exp + 4.0 * g_yyz_xz_x_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xz_x_xz[i] = -2.0 * g_z_xz_x_xz[i] * a_exp + 4.0 * g_yyz_xz_x_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xz_x_yy[i] = -2.0 * g_z_xz_x_yy[i] * a_exp + 4.0 * g_yyz_xz_x_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xz_x_yz[i] = -2.0 * g_z_xz_x_yz[i] * a_exp + 4.0 * g_yyz_xz_x_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xz_x_zz[i] = -2.0 * g_z_xz_x_zz[i] * a_exp + 4.0 * g_yyz_xz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1446-1452)

    #pragma omp simd aligned(g_yyz_xz_y_xx, g_yyz_xz_y_xy, g_yyz_xz_y_xz, g_yyz_xz_y_yy, g_yyz_xz_y_yz, g_yyz_xz_y_zz, g_yz_0_0_0_y_xz_y_xx, g_yz_0_0_0_y_xz_y_xy, g_yz_0_0_0_y_xz_y_xz, g_yz_0_0_0_y_xz_y_yy, g_yz_0_0_0_y_xz_y_yz, g_yz_0_0_0_y_xz_y_zz, g_z_xz_y_xx, g_z_xz_y_xy, g_z_xz_y_xz, g_z_xz_y_yy, g_z_xz_y_yz, g_z_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_xz_y_xx[i] = -2.0 * g_z_xz_y_xx[i] * a_exp + 4.0 * g_yyz_xz_y_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xz_y_xy[i] = -2.0 * g_z_xz_y_xy[i] * a_exp + 4.0 * g_yyz_xz_y_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xz_y_xz[i] = -2.0 * g_z_xz_y_xz[i] * a_exp + 4.0 * g_yyz_xz_y_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xz_y_yy[i] = -2.0 * g_z_xz_y_yy[i] * a_exp + 4.0 * g_yyz_xz_y_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xz_y_yz[i] = -2.0 * g_z_xz_y_yz[i] * a_exp + 4.0 * g_yyz_xz_y_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xz_y_zz[i] = -2.0 * g_z_xz_y_zz[i] * a_exp + 4.0 * g_yyz_xz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1452-1458)

    #pragma omp simd aligned(g_yyz_xz_z_xx, g_yyz_xz_z_xy, g_yyz_xz_z_xz, g_yyz_xz_z_yy, g_yyz_xz_z_yz, g_yyz_xz_z_zz, g_yz_0_0_0_y_xz_z_xx, g_yz_0_0_0_y_xz_z_xy, g_yz_0_0_0_y_xz_z_xz, g_yz_0_0_0_y_xz_z_yy, g_yz_0_0_0_y_xz_z_yz, g_yz_0_0_0_y_xz_z_zz, g_z_xz_z_xx, g_z_xz_z_xy, g_z_xz_z_xz, g_z_xz_z_yy, g_z_xz_z_yz, g_z_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_xz_z_xx[i] = -2.0 * g_z_xz_z_xx[i] * a_exp + 4.0 * g_yyz_xz_z_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xz_z_xy[i] = -2.0 * g_z_xz_z_xy[i] * a_exp + 4.0 * g_yyz_xz_z_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xz_z_xz[i] = -2.0 * g_z_xz_z_xz[i] * a_exp + 4.0 * g_yyz_xz_z_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xz_z_yy[i] = -2.0 * g_z_xz_z_yy[i] * a_exp + 4.0 * g_yyz_xz_z_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xz_z_yz[i] = -2.0 * g_z_xz_z_yz[i] * a_exp + 4.0 * g_yyz_xz_z_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xz_z_zz[i] = -2.0 * g_z_xz_z_zz[i] * a_exp + 4.0 * g_yyz_xz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1458-1464)

    #pragma omp simd aligned(g_yyz_yy_x_xx, g_yyz_yy_x_xy, g_yyz_yy_x_xz, g_yyz_yy_x_yy, g_yyz_yy_x_yz, g_yyz_yy_x_zz, g_yz_0_0_0_y_yy_x_xx, g_yz_0_0_0_y_yy_x_xy, g_yz_0_0_0_y_yy_x_xz, g_yz_0_0_0_y_yy_x_yy, g_yz_0_0_0_y_yy_x_yz, g_yz_0_0_0_y_yy_x_zz, g_z_yy_x_xx, g_z_yy_x_xy, g_z_yy_x_xz, g_z_yy_x_yy, g_z_yy_x_yz, g_z_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_yy_x_xx[i] = -2.0 * g_z_yy_x_xx[i] * a_exp + 4.0 * g_yyz_yy_x_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yy_x_xy[i] = -2.0 * g_z_yy_x_xy[i] * a_exp + 4.0 * g_yyz_yy_x_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yy_x_xz[i] = -2.0 * g_z_yy_x_xz[i] * a_exp + 4.0 * g_yyz_yy_x_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yy_x_yy[i] = -2.0 * g_z_yy_x_yy[i] * a_exp + 4.0 * g_yyz_yy_x_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yy_x_yz[i] = -2.0 * g_z_yy_x_yz[i] * a_exp + 4.0 * g_yyz_yy_x_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yy_x_zz[i] = -2.0 * g_z_yy_x_zz[i] * a_exp + 4.0 * g_yyz_yy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1464-1470)

    #pragma omp simd aligned(g_yyz_yy_y_xx, g_yyz_yy_y_xy, g_yyz_yy_y_xz, g_yyz_yy_y_yy, g_yyz_yy_y_yz, g_yyz_yy_y_zz, g_yz_0_0_0_y_yy_y_xx, g_yz_0_0_0_y_yy_y_xy, g_yz_0_0_0_y_yy_y_xz, g_yz_0_0_0_y_yy_y_yy, g_yz_0_0_0_y_yy_y_yz, g_yz_0_0_0_y_yy_y_zz, g_z_yy_y_xx, g_z_yy_y_xy, g_z_yy_y_xz, g_z_yy_y_yy, g_z_yy_y_yz, g_z_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_yy_y_xx[i] = -2.0 * g_z_yy_y_xx[i] * a_exp + 4.0 * g_yyz_yy_y_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yy_y_xy[i] = -2.0 * g_z_yy_y_xy[i] * a_exp + 4.0 * g_yyz_yy_y_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yy_y_xz[i] = -2.0 * g_z_yy_y_xz[i] * a_exp + 4.0 * g_yyz_yy_y_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yy_y_yy[i] = -2.0 * g_z_yy_y_yy[i] * a_exp + 4.0 * g_yyz_yy_y_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yy_y_yz[i] = -2.0 * g_z_yy_y_yz[i] * a_exp + 4.0 * g_yyz_yy_y_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yy_y_zz[i] = -2.0 * g_z_yy_y_zz[i] * a_exp + 4.0 * g_yyz_yy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1470-1476)

    #pragma omp simd aligned(g_yyz_yy_z_xx, g_yyz_yy_z_xy, g_yyz_yy_z_xz, g_yyz_yy_z_yy, g_yyz_yy_z_yz, g_yyz_yy_z_zz, g_yz_0_0_0_y_yy_z_xx, g_yz_0_0_0_y_yy_z_xy, g_yz_0_0_0_y_yy_z_xz, g_yz_0_0_0_y_yy_z_yy, g_yz_0_0_0_y_yy_z_yz, g_yz_0_0_0_y_yy_z_zz, g_z_yy_z_xx, g_z_yy_z_xy, g_z_yy_z_xz, g_z_yy_z_yy, g_z_yy_z_yz, g_z_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_yy_z_xx[i] = -2.0 * g_z_yy_z_xx[i] * a_exp + 4.0 * g_yyz_yy_z_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yy_z_xy[i] = -2.0 * g_z_yy_z_xy[i] * a_exp + 4.0 * g_yyz_yy_z_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yy_z_xz[i] = -2.0 * g_z_yy_z_xz[i] * a_exp + 4.0 * g_yyz_yy_z_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yy_z_yy[i] = -2.0 * g_z_yy_z_yy[i] * a_exp + 4.0 * g_yyz_yy_z_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yy_z_yz[i] = -2.0 * g_z_yy_z_yz[i] * a_exp + 4.0 * g_yyz_yy_z_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yy_z_zz[i] = -2.0 * g_z_yy_z_zz[i] * a_exp + 4.0 * g_yyz_yy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1476-1482)

    #pragma omp simd aligned(g_yyz_yz_x_xx, g_yyz_yz_x_xy, g_yyz_yz_x_xz, g_yyz_yz_x_yy, g_yyz_yz_x_yz, g_yyz_yz_x_zz, g_yz_0_0_0_y_yz_x_xx, g_yz_0_0_0_y_yz_x_xy, g_yz_0_0_0_y_yz_x_xz, g_yz_0_0_0_y_yz_x_yy, g_yz_0_0_0_y_yz_x_yz, g_yz_0_0_0_y_yz_x_zz, g_z_yz_x_xx, g_z_yz_x_xy, g_z_yz_x_xz, g_z_yz_x_yy, g_z_yz_x_yz, g_z_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_yz_x_xx[i] = -2.0 * g_z_yz_x_xx[i] * a_exp + 4.0 * g_yyz_yz_x_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yz_x_xy[i] = -2.0 * g_z_yz_x_xy[i] * a_exp + 4.0 * g_yyz_yz_x_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yz_x_xz[i] = -2.0 * g_z_yz_x_xz[i] * a_exp + 4.0 * g_yyz_yz_x_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yz_x_yy[i] = -2.0 * g_z_yz_x_yy[i] * a_exp + 4.0 * g_yyz_yz_x_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yz_x_yz[i] = -2.0 * g_z_yz_x_yz[i] * a_exp + 4.0 * g_yyz_yz_x_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yz_x_zz[i] = -2.0 * g_z_yz_x_zz[i] * a_exp + 4.0 * g_yyz_yz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1482-1488)

    #pragma omp simd aligned(g_yyz_yz_y_xx, g_yyz_yz_y_xy, g_yyz_yz_y_xz, g_yyz_yz_y_yy, g_yyz_yz_y_yz, g_yyz_yz_y_zz, g_yz_0_0_0_y_yz_y_xx, g_yz_0_0_0_y_yz_y_xy, g_yz_0_0_0_y_yz_y_xz, g_yz_0_0_0_y_yz_y_yy, g_yz_0_0_0_y_yz_y_yz, g_yz_0_0_0_y_yz_y_zz, g_z_yz_y_xx, g_z_yz_y_xy, g_z_yz_y_xz, g_z_yz_y_yy, g_z_yz_y_yz, g_z_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_yz_y_xx[i] = -2.0 * g_z_yz_y_xx[i] * a_exp + 4.0 * g_yyz_yz_y_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yz_y_xy[i] = -2.0 * g_z_yz_y_xy[i] * a_exp + 4.0 * g_yyz_yz_y_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yz_y_xz[i] = -2.0 * g_z_yz_y_xz[i] * a_exp + 4.0 * g_yyz_yz_y_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yz_y_yy[i] = -2.0 * g_z_yz_y_yy[i] * a_exp + 4.0 * g_yyz_yz_y_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yz_y_yz[i] = -2.0 * g_z_yz_y_yz[i] * a_exp + 4.0 * g_yyz_yz_y_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yz_y_zz[i] = -2.0 * g_z_yz_y_zz[i] * a_exp + 4.0 * g_yyz_yz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1488-1494)

    #pragma omp simd aligned(g_yyz_yz_z_xx, g_yyz_yz_z_xy, g_yyz_yz_z_xz, g_yyz_yz_z_yy, g_yyz_yz_z_yz, g_yyz_yz_z_zz, g_yz_0_0_0_y_yz_z_xx, g_yz_0_0_0_y_yz_z_xy, g_yz_0_0_0_y_yz_z_xz, g_yz_0_0_0_y_yz_z_yy, g_yz_0_0_0_y_yz_z_yz, g_yz_0_0_0_y_yz_z_zz, g_z_yz_z_xx, g_z_yz_z_xy, g_z_yz_z_xz, g_z_yz_z_yy, g_z_yz_z_yz, g_z_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_yz_z_xx[i] = -2.0 * g_z_yz_z_xx[i] * a_exp + 4.0 * g_yyz_yz_z_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yz_z_xy[i] = -2.0 * g_z_yz_z_xy[i] * a_exp + 4.0 * g_yyz_yz_z_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yz_z_xz[i] = -2.0 * g_z_yz_z_xz[i] * a_exp + 4.0 * g_yyz_yz_z_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yz_z_yy[i] = -2.0 * g_z_yz_z_yy[i] * a_exp + 4.0 * g_yyz_yz_z_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yz_z_yz[i] = -2.0 * g_z_yz_z_yz[i] * a_exp + 4.0 * g_yyz_yz_z_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yz_z_zz[i] = -2.0 * g_z_yz_z_zz[i] * a_exp + 4.0 * g_yyz_yz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1494-1500)

    #pragma omp simd aligned(g_yyz_zz_x_xx, g_yyz_zz_x_xy, g_yyz_zz_x_xz, g_yyz_zz_x_yy, g_yyz_zz_x_yz, g_yyz_zz_x_zz, g_yz_0_0_0_y_zz_x_xx, g_yz_0_0_0_y_zz_x_xy, g_yz_0_0_0_y_zz_x_xz, g_yz_0_0_0_y_zz_x_yy, g_yz_0_0_0_y_zz_x_yz, g_yz_0_0_0_y_zz_x_zz, g_z_zz_x_xx, g_z_zz_x_xy, g_z_zz_x_xz, g_z_zz_x_yy, g_z_zz_x_yz, g_z_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_zz_x_xx[i] = -2.0 * g_z_zz_x_xx[i] * a_exp + 4.0 * g_yyz_zz_x_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_zz_x_xy[i] = -2.0 * g_z_zz_x_xy[i] * a_exp + 4.0 * g_yyz_zz_x_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_zz_x_xz[i] = -2.0 * g_z_zz_x_xz[i] * a_exp + 4.0 * g_yyz_zz_x_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_zz_x_yy[i] = -2.0 * g_z_zz_x_yy[i] * a_exp + 4.0 * g_yyz_zz_x_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_zz_x_yz[i] = -2.0 * g_z_zz_x_yz[i] * a_exp + 4.0 * g_yyz_zz_x_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_zz_x_zz[i] = -2.0 * g_z_zz_x_zz[i] * a_exp + 4.0 * g_yyz_zz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1500-1506)

    #pragma omp simd aligned(g_yyz_zz_y_xx, g_yyz_zz_y_xy, g_yyz_zz_y_xz, g_yyz_zz_y_yy, g_yyz_zz_y_yz, g_yyz_zz_y_zz, g_yz_0_0_0_y_zz_y_xx, g_yz_0_0_0_y_zz_y_xy, g_yz_0_0_0_y_zz_y_xz, g_yz_0_0_0_y_zz_y_yy, g_yz_0_0_0_y_zz_y_yz, g_yz_0_0_0_y_zz_y_zz, g_z_zz_y_xx, g_z_zz_y_xy, g_z_zz_y_xz, g_z_zz_y_yy, g_z_zz_y_yz, g_z_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_zz_y_xx[i] = -2.0 * g_z_zz_y_xx[i] * a_exp + 4.0 * g_yyz_zz_y_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_zz_y_xy[i] = -2.0 * g_z_zz_y_xy[i] * a_exp + 4.0 * g_yyz_zz_y_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_zz_y_xz[i] = -2.0 * g_z_zz_y_xz[i] * a_exp + 4.0 * g_yyz_zz_y_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_zz_y_yy[i] = -2.0 * g_z_zz_y_yy[i] * a_exp + 4.0 * g_yyz_zz_y_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_zz_y_yz[i] = -2.0 * g_z_zz_y_yz[i] * a_exp + 4.0 * g_yyz_zz_y_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_zz_y_zz[i] = -2.0 * g_z_zz_y_zz[i] * a_exp + 4.0 * g_yyz_zz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1506-1512)

    #pragma omp simd aligned(g_yyz_zz_z_xx, g_yyz_zz_z_xy, g_yyz_zz_z_xz, g_yyz_zz_z_yy, g_yyz_zz_z_yz, g_yyz_zz_z_zz, g_yz_0_0_0_y_zz_z_xx, g_yz_0_0_0_y_zz_z_xy, g_yz_0_0_0_y_zz_z_xz, g_yz_0_0_0_y_zz_z_yy, g_yz_0_0_0_y_zz_z_yz, g_yz_0_0_0_y_zz_z_zz, g_z_zz_z_xx, g_z_zz_z_xy, g_z_zz_z_xz, g_z_zz_z_yy, g_z_zz_z_yz, g_z_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_zz_z_xx[i] = -2.0 * g_z_zz_z_xx[i] * a_exp + 4.0 * g_yyz_zz_z_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_y_zz_z_xy[i] = -2.0 * g_z_zz_z_xy[i] * a_exp + 4.0 * g_yyz_zz_z_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_zz_z_xz[i] = -2.0 * g_z_zz_z_xz[i] * a_exp + 4.0 * g_yyz_zz_z_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_zz_z_yy[i] = -2.0 * g_z_zz_z_yy[i] * a_exp + 4.0 * g_yyz_zz_z_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_y_zz_z_yz[i] = -2.0 * g_z_zz_z_yz[i] * a_exp + 4.0 * g_yyz_zz_z_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_y_zz_z_zz[i] = -2.0 * g_z_zz_z_zz[i] * a_exp + 4.0 * g_yyz_zz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1512-1518)

    #pragma omp simd aligned(g_y_xx_x_xx, g_y_xx_x_xy, g_y_xx_x_xz, g_y_xx_x_yy, g_y_xx_x_yz, g_y_xx_x_zz, g_yz_0_0_0_z_xx_x_xx, g_yz_0_0_0_z_xx_x_xy, g_yz_0_0_0_z_xx_x_xz, g_yz_0_0_0_z_xx_x_yy, g_yz_0_0_0_z_xx_x_yz, g_yz_0_0_0_z_xx_x_zz, g_yzz_xx_x_xx, g_yzz_xx_x_xy, g_yzz_xx_x_xz, g_yzz_xx_x_yy, g_yzz_xx_x_yz, g_yzz_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_xx_x_xx[i] = -2.0 * g_y_xx_x_xx[i] * a_exp + 4.0 * g_yzz_xx_x_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xx_x_xy[i] = -2.0 * g_y_xx_x_xy[i] * a_exp + 4.0 * g_yzz_xx_x_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xx_x_xz[i] = -2.0 * g_y_xx_x_xz[i] * a_exp + 4.0 * g_yzz_xx_x_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xx_x_yy[i] = -2.0 * g_y_xx_x_yy[i] * a_exp + 4.0 * g_yzz_xx_x_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xx_x_yz[i] = -2.0 * g_y_xx_x_yz[i] * a_exp + 4.0 * g_yzz_xx_x_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xx_x_zz[i] = -2.0 * g_y_xx_x_zz[i] * a_exp + 4.0 * g_yzz_xx_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1518-1524)

    #pragma omp simd aligned(g_y_xx_y_xx, g_y_xx_y_xy, g_y_xx_y_xz, g_y_xx_y_yy, g_y_xx_y_yz, g_y_xx_y_zz, g_yz_0_0_0_z_xx_y_xx, g_yz_0_0_0_z_xx_y_xy, g_yz_0_0_0_z_xx_y_xz, g_yz_0_0_0_z_xx_y_yy, g_yz_0_0_0_z_xx_y_yz, g_yz_0_0_0_z_xx_y_zz, g_yzz_xx_y_xx, g_yzz_xx_y_xy, g_yzz_xx_y_xz, g_yzz_xx_y_yy, g_yzz_xx_y_yz, g_yzz_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_xx_y_xx[i] = -2.0 * g_y_xx_y_xx[i] * a_exp + 4.0 * g_yzz_xx_y_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xx_y_xy[i] = -2.0 * g_y_xx_y_xy[i] * a_exp + 4.0 * g_yzz_xx_y_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xx_y_xz[i] = -2.0 * g_y_xx_y_xz[i] * a_exp + 4.0 * g_yzz_xx_y_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xx_y_yy[i] = -2.0 * g_y_xx_y_yy[i] * a_exp + 4.0 * g_yzz_xx_y_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xx_y_yz[i] = -2.0 * g_y_xx_y_yz[i] * a_exp + 4.0 * g_yzz_xx_y_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xx_y_zz[i] = -2.0 * g_y_xx_y_zz[i] * a_exp + 4.0 * g_yzz_xx_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1524-1530)

    #pragma omp simd aligned(g_y_xx_z_xx, g_y_xx_z_xy, g_y_xx_z_xz, g_y_xx_z_yy, g_y_xx_z_yz, g_y_xx_z_zz, g_yz_0_0_0_z_xx_z_xx, g_yz_0_0_0_z_xx_z_xy, g_yz_0_0_0_z_xx_z_xz, g_yz_0_0_0_z_xx_z_yy, g_yz_0_0_0_z_xx_z_yz, g_yz_0_0_0_z_xx_z_zz, g_yzz_xx_z_xx, g_yzz_xx_z_xy, g_yzz_xx_z_xz, g_yzz_xx_z_yy, g_yzz_xx_z_yz, g_yzz_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_xx_z_xx[i] = -2.0 * g_y_xx_z_xx[i] * a_exp + 4.0 * g_yzz_xx_z_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xx_z_xy[i] = -2.0 * g_y_xx_z_xy[i] * a_exp + 4.0 * g_yzz_xx_z_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xx_z_xz[i] = -2.0 * g_y_xx_z_xz[i] * a_exp + 4.0 * g_yzz_xx_z_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xx_z_yy[i] = -2.0 * g_y_xx_z_yy[i] * a_exp + 4.0 * g_yzz_xx_z_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xx_z_yz[i] = -2.0 * g_y_xx_z_yz[i] * a_exp + 4.0 * g_yzz_xx_z_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xx_z_zz[i] = -2.0 * g_y_xx_z_zz[i] * a_exp + 4.0 * g_yzz_xx_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1530-1536)

    #pragma omp simd aligned(g_y_xy_x_xx, g_y_xy_x_xy, g_y_xy_x_xz, g_y_xy_x_yy, g_y_xy_x_yz, g_y_xy_x_zz, g_yz_0_0_0_z_xy_x_xx, g_yz_0_0_0_z_xy_x_xy, g_yz_0_0_0_z_xy_x_xz, g_yz_0_0_0_z_xy_x_yy, g_yz_0_0_0_z_xy_x_yz, g_yz_0_0_0_z_xy_x_zz, g_yzz_xy_x_xx, g_yzz_xy_x_xy, g_yzz_xy_x_xz, g_yzz_xy_x_yy, g_yzz_xy_x_yz, g_yzz_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_xy_x_xx[i] = -2.0 * g_y_xy_x_xx[i] * a_exp + 4.0 * g_yzz_xy_x_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xy_x_xy[i] = -2.0 * g_y_xy_x_xy[i] * a_exp + 4.0 * g_yzz_xy_x_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xy_x_xz[i] = -2.0 * g_y_xy_x_xz[i] * a_exp + 4.0 * g_yzz_xy_x_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xy_x_yy[i] = -2.0 * g_y_xy_x_yy[i] * a_exp + 4.0 * g_yzz_xy_x_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xy_x_yz[i] = -2.0 * g_y_xy_x_yz[i] * a_exp + 4.0 * g_yzz_xy_x_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xy_x_zz[i] = -2.0 * g_y_xy_x_zz[i] * a_exp + 4.0 * g_yzz_xy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1536-1542)

    #pragma omp simd aligned(g_y_xy_y_xx, g_y_xy_y_xy, g_y_xy_y_xz, g_y_xy_y_yy, g_y_xy_y_yz, g_y_xy_y_zz, g_yz_0_0_0_z_xy_y_xx, g_yz_0_0_0_z_xy_y_xy, g_yz_0_0_0_z_xy_y_xz, g_yz_0_0_0_z_xy_y_yy, g_yz_0_0_0_z_xy_y_yz, g_yz_0_0_0_z_xy_y_zz, g_yzz_xy_y_xx, g_yzz_xy_y_xy, g_yzz_xy_y_xz, g_yzz_xy_y_yy, g_yzz_xy_y_yz, g_yzz_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_xy_y_xx[i] = -2.0 * g_y_xy_y_xx[i] * a_exp + 4.0 * g_yzz_xy_y_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xy_y_xy[i] = -2.0 * g_y_xy_y_xy[i] * a_exp + 4.0 * g_yzz_xy_y_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xy_y_xz[i] = -2.0 * g_y_xy_y_xz[i] * a_exp + 4.0 * g_yzz_xy_y_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xy_y_yy[i] = -2.0 * g_y_xy_y_yy[i] * a_exp + 4.0 * g_yzz_xy_y_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xy_y_yz[i] = -2.0 * g_y_xy_y_yz[i] * a_exp + 4.0 * g_yzz_xy_y_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xy_y_zz[i] = -2.0 * g_y_xy_y_zz[i] * a_exp + 4.0 * g_yzz_xy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1542-1548)

    #pragma omp simd aligned(g_y_xy_z_xx, g_y_xy_z_xy, g_y_xy_z_xz, g_y_xy_z_yy, g_y_xy_z_yz, g_y_xy_z_zz, g_yz_0_0_0_z_xy_z_xx, g_yz_0_0_0_z_xy_z_xy, g_yz_0_0_0_z_xy_z_xz, g_yz_0_0_0_z_xy_z_yy, g_yz_0_0_0_z_xy_z_yz, g_yz_0_0_0_z_xy_z_zz, g_yzz_xy_z_xx, g_yzz_xy_z_xy, g_yzz_xy_z_xz, g_yzz_xy_z_yy, g_yzz_xy_z_yz, g_yzz_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_xy_z_xx[i] = -2.0 * g_y_xy_z_xx[i] * a_exp + 4.0 * g_yzz_xy_z_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xy_z_xy[i] = -2.0 * g_y_xy_z_xy[i] * a_exp + 4.0 * g_yzz_xy_z_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xy_z_xz[i] = -2.0 * g_y_xy_z_xz[i] * a_exp + 4.0 * g_yzz_xy_z_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xy_z_yy[i] = -2.0 * g_y_xy_z_yy[i] * a_exp + 4.0 * g_yzz_xy_z_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xy_z_yz[i] = -2.0 * g_y_xy_z_yz[i] * a_exp + 4.0 * g_yzz_xy_z_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xy_z_zz[i] = -2.0 * g_y_xy_z_zz[i] * a_exp + 4.0 * g_yzz_xy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1548-1554)

    #pragma omp simd aligned(g_y_xz_x_xx, g_y_xz_x_xy, g_y_xz_x_xz, g_y_xz_x_yy, g_y_xz_x_yz, g_y_xz_x_zz, g_yz_0_0_0_z_xz_x_xx, g_yz_0_0_0_z_xz_x_xy, g_yz_0_0_0_z_xz_x_xz, g_yz_0_0_0_z_xz_x_yy, g_yz_0_0_0_z_xz_x_yz, g_yz_0_0_0_z_xz_x_zz, g_yzz_xz_x_xx, g_yzz_xz_x_xy, g_yzz_xz_x_xz, g_yzz_xz_x_yy, g_yzz_xz_x_yz, g_yzz_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_xz_x_xx[i] = -2.0 * g_y_xz_x_xx[i] * a_exp + 4.0 * g_yzz_xz_x_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xz_x_xy[i] = -2.0 * g_y_xz_x_xy[i] * a_exp + 4.0 * g_yzz_xz_x_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xz_x_xz[i] = -2.0 * g_y_xz_x_xz[i] * a_exp + 4.0 * g_yzz_xz_x_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xz_x_yy[i] = -2.0 * g_y_xz_x_yy[i] * a_exp + 4.0 * g_yzz_xz_x_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xz_x_yz[i] = -2.0 * g_y_xz_x_yz[i] * a_exp + 4.0 * g_yzz_xz_x_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xz_x_zz[i] = -2.0 * g_y_xz_x_zz[i] * a_exp + 4.0 * g_yzz_xz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1554-1560)

    #pragma omp simd aligned(g_y_xz_y_xx, g_y_xz_y_xy, g_y_xz_y_xz, g_y_xz_y_yy, g_y_xz_y_yz, g_y_xz_y_zz, g_yz_0_0_0_z_xz_y_xx, g_yz_0_0_0_z_xz_y_xy, g_yz_0_0_0_z_xz_y_xz, g_yz_0_0_0_z_xz_y_yy, g_yz_0_0_0_z_xz_y_yz, g_yz_0_0_0_z_xz_y_zz, g_yzz_xz_y_xx, g_yzz_xz_y_xy, g_yzz_xz_y_xz, g_yzz_xz_y_yy, g_yzz_xz_y_yz, g_yzz_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_xz_y_xx[i] = -2.0 * g_y_xz_y_xx[i] * a_exp + 4.0 * g_yzz_xz_y_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xz_y_xy[i] = -2.0 * g_y_xz_y_xy[i] * a_exp + 4.0 * g_yzz_xz_y_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xz_y_xz[i] = -2.0 * g_y_xz_y_xz[i] * a_exp + 4.0 * g_yzz_xz_y_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xz_y_yy[i] = -2.0 * g_y_xz_y_yy[i] * a_exp + 4.0 * g_yzz_xz_y_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xz_y_yz[i] = -2.0 * g_y_xz_y_yz[i] * a_exp + 4.0 * g_yzz_xz_y_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xz_y_zz[i] = -2.0 * g_y_xz_y_zz[i] * a_exp + 4.0 * g_yzz_xz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1560-1566)

    #pragma omp simd aligned(g_y_xz_z_xx, g_y_xz_z_xy, g_y_xz_z_xz, g_y_xz_z_yy, g_y_xz_z_yz, g_y_xz_z_zz, g_yz_0_0_0_z_xz_z_xx, g_yz_0_0_0_z_xz_z_xy, g_yz_0_0_0_z_xz_z_xz, g_yz_0_0_0_z_xz_z_yy, g_yz_0_0_0_z_xz_z_yz, g_yz_0_0_0_z_xz_z_zz, g_yzz_xz_z_xx, g_yzz_xz_z_xy, g_yzz_xz_z_xz, g_yzz_xz_z_yy, g_yzz_xz_z_yz, g_yzz_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_xz_z_xx[i] = -2.0 * g_y_xz_z_xx[i] * a_exp + 4.0 * g_yzz_xz_z_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xz_z_xy[i] = -2.0 * g_y_xz_z_xy[i] * a_exp + 4.0 * g_yzz_xz_z_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xz_z_xz[i] = -2.0 * g_y_xz_z_xz[i] * a_exp + 4.0 * g_yzz_xz_z_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xz_z_yy[i] = -2.0 * g_y_xz_z_yy[i] * a_exp + 4.0 * g_yzz_xz_z_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xz_z_yz[i] = -2.0 * g_y_xz_z_yz[i] * a_exp + 4.0 * g_yzz_xz_z_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xz_z_zz[i] = -2.0 * g_y_xz_z_zz[i] * a_exp + 4.0 * g_yzz_xz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1566-1572)

    #pragma omp simd aligned(g_y_yy_x_xx, g_y_yy_x_xy, g_y_yy_x_xz, g_y_yy_x_yy, g_y_yy_x_yz, g_y_yy_x_zz, g_yz_0_0_0_z_yy_x_xx, g_yz_0_0_0_z_yy_x_xy, g_yz_0_0_0_z_yy_x_xz, g_yz_0_0_0_z_yy_x_yy, g_yz_0_0_0_z_yy_x_yz, g_yz_0_0_0_z_yy_x_zz, g_yzz_yy_x_xx, g_yzz_yy_x_xy, g_yzz_yy_x_xz, g_yzz_yy_x_yy, g_yzz_yy_x_yz, g_yzz_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_yy_x_xx[i] = -2.0 * g_y_yy_x_xx[i] * a_exp + 4.0 * g_yzz_yy_x_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yy_x_xy[i] = -2.0 * g_y_yy_x_xy[i] * a_exp + 4.0 * g_yzz_yy_x_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yy_x_xz[i] = -2.0 * g_y_yy_x_xz[i] * a_exp + 4.0 * g_yzz_yy_x_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yy_x_yy[i] = -2.0 * g_y_yy_x_yy[i] * a_exp + 4.0 * g_yzz_yy_x_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yy_x_yz[i] = -2.0 * g_y_yy_x_yz[i] * a_exp + 4.0 * g_yzz_yy_x_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yy_x_zz[i] = -2.0 * g_y_yy_x_zz[i] * a_exp + 4.0 * g_yzz_yy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1572-1578)

    #pragma omp simd aligned(g_y_yy_y_xx, g_y_yy_y_xy, g_y_yy_y_xz, g_y_yy_y_yy, g_y_yy_y_yz, g_y_yy_y_zz, g_yz_0_0_0_z_yy_y_xx, g_yz_0_0_0_z_yy_y_xy, g_yz_0_0_0_z_yy_y_xz, g_yz_0_0_0_z_yy_y_yy, g_yz_0_0_0_z_yy_y_yz, g_yz_0_0_0_z_yy_y_zz, g_yzz_yy_y_xx, g_yzz_yy_y_xy, g_yzz_yy_y_xz, g_yzz_yy_y_yy, g_yzz_yy_y_yz, g_yzz_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_yy_y_xx[i] = -2.0 * g_y_yy_y_xx[i] * a_exp + 4.0 * g_yzz_yy_y_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yy_y_xy[i] = -2.0 * g_y_yy_y_xy[i] * a_exp + 4.0 * g_yzz_yy_y_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yy_y_xz[i] = -2.0 * g_y_yy_y_xz[i] * a_exp + 4.0 * g_yzz_yy_y_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yy_y_yy[i] = -2.0 * g_y_yy_y_yy[i] * a_exp + 4.0 * g_yzz_yy_y_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yy_y_yz[i] = -2.0 * g_y_yy_y_yz[i] * a_exp + 4.0 * g_yzz_yy_y_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yy_y_zz[i] = -2.0 * g_y_yy_y_zz[i] * a_exp + 4.0 * g_yzz_yy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1578-1584)

    #pragma omp simd aligned(g_y_yy_z_xx, g_y_yy_z_xy, g_y_yy_z_xz, g_y_yy_z_yy, g_y_yy_z_yz, g_y_yy_z_zz, g_yz_0_0_0_z_yy_z_xx, g_yz_0_0_0_z_yy_z_xy, g_yz_0_0_0_z_yy_z_xz, g_yz_0_0_0_z_yy_z_yy, g_yz_0_0_0_z_yy_z_yz, g_yz_0_0_0_z_yy_z_zz, g_yzz_yy_z_xx, g_yzz_yy_z_xy, g_yzz_yy_z_xz, g_yzz_yy_z_yy, g_yzz_yy_z_yz, g_yzz_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_yy_z_xx[i] = -2.0 * g_y_yy_z_xx[i] * a_exp + 4.0 * g_yzz_yy_z_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yy_z_xy[i] = -2.0 * g_y_yy_z_xy[i] * a_exp + 4.0 * g_yzz_yy_z_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yy_z_xz[i] = -2.0 * g_y_yy_z_xz[i] * a_exp + 4.0 * g_yzz_yy_z_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yy_z_yy[i] = -2.0 * g_y_yy_z_yy[i] * a_exp + 4.0 * g_yzz_yy_z_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yy_z_yz[i] = -2.0 * g_y_yy_z_yz[i] * a_exp + 4.0 * g_yzz_yy_z_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yy_z_zz[i] = -2.0 * g_y_yy_z_zz[i] * a_exp + 4.0 * g_yzz_yy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1584-1590)

    #pragma omp simd aligned(g_y_yz_x_xx, g_y_yz_x_xy, g_y_yz_x_xz, g_y_yz_x_yy, g_y_yz_x_yz, g_y_yz_x_zz, g_yz_0_0_0_z_yz_x_xx, g_yz_0_0_0_z_yz_x_xy, g_yz_0_0_0_z_yz_x_xz, g_yz_0_0_0_z_yz_x_yy, g_yz_0_0_0_z_yz_x_yz, g_yz_0_0_0_z_yz_x_zz, g_yzz_yz_x_xx, g_yzz_yz_x_xy, g_yzz_yz_x_xz, g_yzz_yz_x_yy, g_yzz_yz_x_yz, g_yzz_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_yz_x_xx[i] = -2.0 * g_y_yz_x_xx[i] * a_exp + 4.0 * g_yzz_yz_x_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yz_x_xy[i] = -2.0 * g_y_yz_x_xy[i] * a_exp + 4.0 * g_yzz_yz_x_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yz_x_xz[i] = -2.0 * g_y_yz_x_xz[i] * a_exp + 4.0 * g_yzz_yz_x_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yz_x_yy[i] = -2.0 * g_y_yz_x_yy[i] * a_exp + 4.0 * g_yzz_yz_x_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yz_x_yz[i] = -2.0 * g_y_yz_x_yz[i] * a_exp + 4.0 * g_yzz_yz_x_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yz_x_zz[i] = -2.0 * g_y_yz_x_zz[i] * a_exp + 4.0 * g_yzz_yz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1590-1596)

    #pragma omp simd aligned(g_y_yz_y_xx, g_y_yz_y_xy, g_y_yz_y_xz, g_y_yz_y_yy, g_y_yz_y_yz, g_y_yz_y_zz, g_yz_0_0_0_z_yz_y_xx, g_yz_0_0_0_z_yz_y_xy, g_yz_0_0_0_z_yz_y_xz, g_yz_0_0_0_z_yz_y_yy, g_yz_0_0_0_z_yz_y_yz, g_yz_0_0_0_z_yz_y_zz, g_yzz_yz_y_xx, g_yzz_yz_y_xy, g_yzz_yz_y_xz, g_yzz_yz_y_yy, g_yzz_yz_y_yz, g_yzz_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_yz_y_xx[i] = -2.0 * g_y_yz_y_xx[i] * a_exp + 4.0 * g_yzz_yz_y_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yz_y_xy[i] = -2.0 * g_y_yz_y_xy[i] * a_exp + 4.0 * g_yzz_yz_y_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yz_y_xz[i] = -2.0 * g_y_yz_y_xz[i] * a_exp + 4.0 * g_yzz_yz_y_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yz_y_yy[i] = -2.0 * g_y_yz_y_yy[i] * a_exp + 4.0 * g_yzz_yz_y_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yz_y_yz[i] = -2.0 * g_y_yz_y_yz[i] * a_exp + 4.0 * g_yzz_yz_y_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yz_y_zz[i] = -2.0 * g_y_yz_y_zz[i] * a_exp + 4.0 * g_yzz_yz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1596-1602)

    #pragma omp simd aligned(g_y_yz_z_xx, g_y_yz_z_xy, g_y_yz_z_xz, g_y_yz_z_yy, g_y_yz_z_yz, g_y_yz_z_zz, g_yz_0_0_0_z_yz_z_xx, g_yz_0_0_0_z_yz_z_xy, g_yz_0_0_0_z_yz_z_xz, g_yz_0_0_0_z_yz_z_yy, g_yz_0_0_0_z_yz_z_yz, g_yz_0_0_0_z_yz_z_zz, g_yzz_yz_z_xx, g_yzz_yz_z_xy, g_yzz_yz_z_xz, g_yzz_yz_z_yy, g_yzz_yz_z_yz, g_yzz_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_yz_z_xx[i] = -2.0 * g_y_yz_z_xx[i] * a_exp + 4.0 * g_yzz_yz_z_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yz_z_xy[i] = -2.0 * g_y_yz_z_xy[i] * a_exp + 4.0 * g_yzz_yz_z_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yz_z_xz[i] = -2.0 * g_y_yz_z_xz[i] * a_exp + 4.0 * g_yzz_yz_z_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yz_z_yy[i] = -2.0 * g_y_yz_z_yy[i] * a_exp + 4.0 * g_yzz_yz_z_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yz_z_yz[i] = -2.0 * g_y_yz_z_yz[i] * a_exp + 4.0 * g_yzz_yz_z_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yz_z_zz[i] = -2.0 * g_y_yz_z_zz[i] * a_exp + 4.0 * g_yzz_yz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1602-1608)

    #pragma omp simd aligned(g_y_zz_x_xx, g_y_zz_x_xy, g_y_zz_x_xz, g_y_zz_x_yy, g_y_zz_x_yz, g_y_zz_x_zz, g_yz_0_0_0_z_zz_x_xx, g_yz_0_0_0_z_zz_x_xy, g_yz_0_0_0_z_zz_x_xz, g_yz_0_0_0_z_zz_x_yy, g_yz_0_0_0_z_zz_x_yz, g_yz_0_0_0_z_zz_x_zz, g_yzz_zz_x_xx, g_yzz_zz_x_xy, g_yzz_zz_x_xz, g_yzz_zz_x_yy, g_yzz_zz_x_yz, g_yzz_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_zz_x_xx[i] = -2.0 * g_y_zz_x_xx[i] * a_exp + 4.0 * g_yzz_zz_x_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_zz_x_xy[i] = -2.0 * g_y_zz_x_xy[i] * a_exp + 4.0 * g_yzz_zz_x_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_zz_x_xz[i] = -2.0 * g_y_zz_x_xz[i] * a_exp + 4.0 * g_yzz_zz_x_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_zz_x_yy[i] = -2.0 * g_y_zz_x_yy[i] * a_exp + 4.0 * g_yzz_zz_x_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_zz_x_yz[i] = -2.0 * g_y_zz_x_yz[i] * a_exp + 4.0 * g_yzz_zz_x_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_zz_x_zz[i] = -2.0 * g_y_zz_x_zz[i] * a_exp + 4.0 * g_yzz_zz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1608-1614)

    #pragma omp simd aligned(g_y_zz_y_xx, g_y_zz_y_xy, g_y_zz_y_xz, g_y_zz_y_yy, g_y_zz_y_yz, g_y_zz_y_zz, g_yz_0_0_0_z_zz_y_xx, g_yz_0_0_0_z_zz_y_xy, g_yz_0_0_0_z_zz_y_xz, g_yz_0_0_0_z_zz_y_yy, g_yz_0_0_0_z_zz_y_yz, g_yz_0_0_0_z_zz_y_zz, g_yzz_zz_y_xx, g_yzz_zz_y_xy, g_yzz_zz_y_xz, g_yzz_zz_y_yy, g_yzz_zz_y_yz, g_yzz_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_zz_y_xx[i] = -2.0 * g_y_zz_y_xx[i] * a_exp + 4.0 * g_yzz_zz_y_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_zz_y_xy[i] = -2.0 * g_y_zz_y_xy[i] * a_exp + 4.0 * g_yzz_zz_y_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_zz_y_xz[i] = -2.0 * g_y_zz_y_xz[i] * a_exp + 4.0 * g_yzz_zz_y_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_zz_y_yy[i] = -2.0 * g_y_zz_y_yy[i] * a_exp + 4.0 * g_yzz_zz_y_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_zz_y_yz[i] = -2.0 * g_y_zz_y_yz[i] * a_exp + 4.0 * g_yzz_zz_y_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_zz_y_zz[i] = -2.0 * g_y_zz_y_zz[i] * a_exp + 4.0 * g_yzz_zz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1614-1620)

    #pragma omp simd aligned(g_y_zz_z_xx, g_y_zz_z_xy, g_y_zz_z_xz, g_y_zz_z_yy, g_y_zz_z_yz, g_y_zz_z_zz, g_yz_0_0_0_z_zz_z_xx, g_yz_0_0_0_z_zz_z_xy, g_yz_0_0_0_z_zz_z_xz, g_yz_0_0_0_z_zz_z_yy, g_yz_0_0_0_z_zz_z_yz, g_yz_0_0_0_z_zz_z_zz, g_yzz_zz_z_xx, g_yzz_zz_z_xy, g_yzz_zz_z_xz, g_yzz_zz_z_yy, g_yzz_zz_z_yz, g_yzz_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_zz_z_xx[i] = -2.0 * g_y_zz_z_xx[i] * a_exp + 4.0 * g_yzz_zz_z_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_z_zz_z_xy[i] = -2.0 * g_y_zz_z_xy[i] * a_exp + 4.0 * g_yzz_zz_z_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_zz_z_xz[i] = -2.0 * g_y_zz_z_xz[i] * a_exp + 4.0 * g_yzz_zz_z_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_zz_z_yy[i] = -2.0 * g_y_zz_z_yy[i] * a_exp + 4.0 * g_yzz_zz_z_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_z_zz_z_yz[i] = -2.0 * g_y_zz_z_yz[i] * a_exp + 4.0 * g_yzz_zz_z_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_z_zz_z_zz[i] = -2.0 * g_y_zz_z_zz[i] * a_exp + 4.0 * g_yzz_zz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1620-1626)

    #pragma omp simd aligned(g_x_xx_x_xx, g_x_xx_x_xy, g_x_xx_x_xz, g_x_xx_x_yy, g_x_xx_x_yz, g_x_xx_x_zz, g_xzz_xx_x_xx, g_xzz_xx_x_xy, g_xzz_xx_x_xz, g_xzz_xx_x_yy, g_xzz_xx_x_yz, g_xzz_xx_x_zz, g_zz_0_0_0_x_xx_x_xx, g_zz_0_0_0_x_xx_x_xy, g_zz_0_0_0_x_xx_x_xz, g_zz_0_0_0_x_xx_x_yy, g_zz_0_0_0_x_xx_x_yz, g_zz_0_0_0_x_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_xx_x_xx[i] = -2.0 * g_x_xx_x_xx[i] * a_exp + 4.0 * g_xzz_xx_x_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xx_x_xy[i] = -2.0 * g_x_xx_x_xy[i] * a_exp + 4.0 * g_xzz_xx_x_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xx_x_xz[i] = -2.0 * g_x_xx_x_xz[i] * a_exp + 4.0 * g_xzz_xx_x_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xx_x_yy[i] = -2.0 * g_x_xx_x_yy[i] * a_exp + 4.0 * g_xzz_xx_x_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xx_x_yz[i] = -2.0 * g_x_xx_x_yz[i] * a_exp + 4.0 * g_xzz_xx_x_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xx_x_zz[i] = -2.0 * g_x_xx_x_zz[i] * a_exp + 4.0 * g_xzz_xx_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1626-1632)

    #pragma omp simd aligned(g_x_xx_y_xx, g_x_xx_y_xy, g_x_xx_y_xz, g_x_xx_y_yy, g_x_xx_y_yz, g_x_xx_y_zz, g_xzz_xx_y_xx, g_xzz_xx_y_xy, g_xzz_xx_y_xz, g_xzz_xx_y_yy, g_xzz_xx_y_yz, g_xzz_xx_y_zz, g_zz_0_0_0_x_xx_y_xx, g_zz_0_0_0_x_xx_y_xy, g_zz_0_0_0_x_xx_y_xz, g_zz_0_0_0_x_xx_y_yy, g_zz_0_0_0_x_xx_y_yz, g_zz_0_0_0_x_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_xx_y_xx[i] = -2.0 * g_x_xx_y_xx[i] * a_exp + 4.0 * g_xzz_xx_y_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xx_y_xy[i] = -2.0 * g_x_xx_y_xy[i] * a_exp + 4.0 * g_xzz_xx_y_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xx_y_xz[i] = -2.0 * g_x_xx_y_xz[i] * a_exp + 4.0 * g_xzz_xx_y_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xx_y_yy[i] = -2.0 * g_x_xx_y_yy[i] * a_exp + 4.0 * g_xzz_xx_y_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xx_y_yz[i] = -2.0 * g_x_xx_y_yz[i] * a_exp + 4.0 * g_xzz_xx_y_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xx_y_zz[i] = -2.0 * g_x_xx_y_zz[i] * a_exp + 4.0 * g_xzz_xx_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1632-1638)

    #pragma omp simd aligned(g_x_xx_z_xx, g_x_xx_z_xy, g_x_xx_z_xz, g_x_xx_z_yy, g_x_xx_z_yz, g_x_xx_z_zz, g_xzz_xx_z_xx, g_xzz_xx_z_xy, g_xzz_xx_z_xz, g_xzz_xx_z_yy, g_xzz_xx_z_yz, g_xzz_xx_z_zz, g_zz_0_0_0_x_xx_z_xx, g_zz_0_0_0_x_xx_z_xy, g_zz_0_0_0_x_xx_z_xz, g_zz_0_0_0_x_xx_z_yy, g_zz_0_0_0_x_xx_z_yz, g_zz_0_0_0_x_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_xx_z_xx[i] = -2.0 * g_x_xx_z_xx[i] * a_exp + 4.0 * g_xzz_xx_z_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xx_z_xy[i] = -2.0 * g_x_xx_z_xy[i] * a_exp + 4.0 * g_xzz_xx_z_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xx_z_xz[i] = -2.0 * g_x_xx_z_xz[i] * a_exp + 4.0 * g_xzz_xx_z_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xx_z_yy[i] = -2.0 * g_x_xx_z_yy[i] * a_exp + 4.0 * g_xzz_xx_z_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xx_z_yz[i] = -2.0 * g_x_xx_z_yz[i] * a_exp + 4.0 * g_xzz_xx_z_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xx_z_zz[i] = -2.0 * g_x_xx_z_zz[i] * a_exp + 4.0 * g_xzz_xx_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1638-1644)

    #pragma omp simd aligned(g_x_xy_x_xx, g_x_xy_x_xy, g_x_xy_x_xz, g_x_xy_x_yy, g_x_xy_x_yz, g_x_xy_x_zz, g_xzz_xy_x_xx, g_xzz_xy_x_xy, g_xzz_xy_x_xz, g_xzz_xy_x_yy, g_xzz_xy_x_yz, g_xzz_xy_x_zz, g_zz_0_0_0_x_xy_x_xx, g_zz_0_0_0_x_xy_x_xy, g_zz_0_0_0_x_xy_x_xz, g_zz_0_0_0_x_xy_x_yy, g_zz_0_0_0_x_xy_x_yz, g_zz_0_0_0_x_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_xy_x_xx[i] = -2.0 * g_x_xy_x_xx[i] * a_exp + 4.0 * g_xzz_xy_x_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xy_x_xy[i] = -2.0 * g_x_xy_x_xy[i] * a_exp + 4.0 * g_xzz_xy_x_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xy_x_xz[i] = -2.0 * g_x_xy_x_xz[i] * a_exp + 4.0 * g_xzz_xy_x_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xy_x_yy[i] = -2.0 * g_x_xy_x_yy[i] * a_exp + 4.0 * g_xzz_xy_x_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xy_x_yz[i] = -2.0 * g_x_xy_x_yz[i] * a_exp + 4.0 * g_xzz_xy_x_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xy_x_zz[i] = -2.0 * g_x_xy_x_zz[i] * a_exp + 4.0 * g_xzz_xy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1644-1650)

    #pragma omp simd aligned(g_x_xy_y_xx, g_x_xy_y_xy, g_x_xy_y_xz, g_x_xy_y_yy, g_x_xy_y_yz, g_x_xy_y_zz, g_xzz_xy_y_xx, g_xzz_xy_y_xy, g_xzz_xy_y_xz, g_xzz_xy_y_yy, g_xzz_xy_y_yz, g_xzz_xy_y_zz, g_zz_0_0_0_x_xy_y_xx, g_zz_0_0_0_x_xy_y_xy, g_zz_0_0_0_x_xy_y_xz, g_zz_0_0_0_x_xy_y_yy, g_zz_0_0_0_x_xy_y_yz, g_zz_0_0_0_x_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_xy_y_xx[i] = -2.0 * g_x_xy_y_xx[i] * a_exp + 4.0 * g_xzz_xy_y_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xy_y_xy[i] = -2.0 * g_x_xy_y_xy[i] * a_exp + 4.0 * g_xzz_xy_y_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xy_y_xz[i] = -2.0 * g_x_xy_y_xz[i] * a_exp + 4.0 * g_xzz_xy_y_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xy_y_yy[i] = -2.0 * g_x_xy_y_yy[i] * a_exp + 4.0 * g_xzz_xy_y_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xy_y_yz[i] = -2.0 * g_x_xy_y_yz[i] * a_exp + 4.0 * g_xzz_xy_y_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xy_y_zz[i] = -2.0 * g_x_xy_y_zz[i] * a_exp + 4.0 * g_xzz_xy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1650-1656)

    #pragma omp simd aligned(g_x_xy_z_xx, g_x_xy_z_xy, g_x_xy_z_xz, g_x_xy_z_yy, g_x_xy_z_yz, g_x_xy_z_zz, g_xzz_xy_z_xx, g_xzz_xy_z_xy, g_xzz_xy_z_xz, g_xzz_xy_z_yy, g_xzz_xy_z_yz, g_xzz_xy_z_zz, g_zz_0_0_0_x_xy_z_xx, g_zz_0_0_0_x_xy_z_xy, g_zz_0_0_0_x_xy_z_xz, g_zz_0_0_0_x_xy_z_yy, g_zz_0_0_0_x_xy_z_yz, g_zz_0_0_0_x_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_xy_z_xx[i] = -2.0 * g_x_xy_z_xx[i] * a_exp + 4.0 * g_xzz_xy_z_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xy_z_xy[i] = -2.0 * g_x_xy_z_xy[i] * a_exp + 4.0 * g_xzz_xy_z_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xy_z_xz[i] = -2.0 * g_x_xy_z_xz[i] * a_exp + 4.0 * g_xzz_xy_z_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xy_z_yy[i] = -2.0 * g_x_xy_z_yy[i] * a_exp + 4.0 * g_xzz_xy_z_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xy_z_yz[i] = -2.0 * g_x_xy_z_yz[i] * a_exp + 4.0 * g_xzz_xy_z_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xy_z_zz[i] = -2.0 * g_x_xy_z_zz[i] * a_exp + 4.0 * g_xzz_xy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1656-1662)

    #pragma omp simd aligned(g_x_xz_x_xx, g_x_xz_x_xy, g_x_xz_x_xz, g_x_xz_x_yy, g_x_xz_x_yz, g_x_xz_x_zz, g_xzz_xz_x_xx, g_xzz_xz_x_xy, g_xzz_xz_x_xz, g_xzz_xz_x_yy, g_xzz_xz_x_yz, g_xzz_xz_x_zz, g_zz_0_0_0_x_xz_x_xx, g_zz_0_0_0_x_xz_x_xy, g_zz_0_0_0_x_xz_x_xz, g_zz_0_0_0_x_xz_x_yy, g_zz_0_0_0_x_xz_x_yz, g_zz_0_0_0_x_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_xz_x_xx[i] = -2.0 * g_x_xz_x_xx[i] * a_exp + 4.0 * g_xzz_xz_x_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xz_x_xy[i] = -2.0 * g_x_xz_x_xy[i] * a_exp + 4.0 * g_xzz_xz_x_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xz_x_xz[i] = -2.0 * g_x_xz_x_xz[i] * a_exp + 4.0 * g_xzz_xz_x_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xz_x_yy[i] = -2.0 * g_x_xz_x_yy[i] * a_exp + 4.0 * g_xzz_xz_x_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xz_x_yz[i] = -2.0 * g_x_xz_x_yz[i] * a_exp + 4.0 * g_xzz_xz_x_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xz_x_zz[i] = -2.0 * g_x_xz_x_zz[i] * a_exp + 4.0 * g_xzz_xz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1662-1668)

    #pragma omp simd aligned(g_x_xz_y_xx, g_x_xz_y_xy, g_x_xz_y_xz, g_x_xz_y_yy, g_x_xz_y_yz, g_x_xz_y_zz, g_xzz_xz_y_xx, g_xzz_xz_y_xy, g_xzz_xz_y_xz, g_xzz_xz_y_yy, g_xzz_xz_y_yz, g_xzz_xz_y_zz, g_zz_0_0_0_x_xz_y_xx, g_zz_0_0_0_x_xz_y_xy, g_zz_0_0_0_x_xz_y_xz, g_zz_0_0_0_x_xz_y_yy, g_zz_0_0_0_x_xz_y_yz, g_zz_0_0_0_x_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_xz_y_xx[i] = -2.0 * g_x_xz_y_xx[i] * a_exp + 4.0 * g_xzz_xz_y_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xz_y_xy[i] = -2.0 * g_x_xz_y_xy[i] * a_exp + 4.0 * g_xzz_xz_y_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xz_y_xz[i] = -2.0 * g_x_xz_y_xz[i] * a_exp + 4.0 * g_xzz_xz_y_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xz_y_yy[i] = -2.0 * g_x_xz_y_yy[i] * a_exp + 4.0 * g_xzz_xz_y_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xz_y_yz[i] = -2.0 * g_x_xz_y_yz[i] * a_exp + 4.0 * g_xzz_xz_y_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xz_y_zz[i] = -2.0 * g_x_xz_y_zz[i] * a_exp + 4.0 * g_xzz_xz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1668-1674)

    #pragma omp simd aligned(g_x_xz_z_xx, g_x_xz_z_xy, g_x_xz_z_xz, g_x_xz_z_yy, g_x_xz_z_yz, g_x_xz_z_zz, g_xzz_xz_z_xx, g_xzz_xz_z_xy, g_xzz_xz_z_xz, g_xzz_xz_z_yy, g_xzz_xz_z_yz, g_xzz_xz_z_zz, g_zz_0_0_0_x_xz_z_xx, g_zz_0_0_0_x_xz_z_xy, g_zz_0_0_0_x_xz_z_xz, g_zz_0_0_0_x_xz_z_yy, g_zz_0_0_0_x_xz_z_yz, g_zz_0_0_0_x_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_xz_z_xx[i] = -2.0 * g_x_xz_z_xx[i] * a_exp + 4.0 * g_xzz_xz_z_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xz_z_xy[i] = -2.0 * g_x_xz_z_xy[i] * a_exp + 4.0 * g_xzz_xz_z_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xz_z_xz[i] = -2.0 * g_x_xz_z_xz[i] * a_exp + 4.0 * g_xzz_xz_z_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xz_z_yy[i] = -2.0 * g_x_xz_z_yy[i] * a_exp + 4.0 * g_xzz_xz_z_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xz_z_yz[i] = -2.0 * g_x_xz_z_yz[i] * a_exp + 4.0 * g_xzz_xz_z_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xz_z_zz[i] = -2.0 * g_x_xz_z_zz[i] * a_exp + 4.0 * g_xzz_xz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1674-1680)

    #pragma omp simd aligned(g_x_yy_x_xx, g_x_yy_x_xy, g_x_yy_x_xz, g_x_yy_x_yy, g_x_yy_x_yz, g_x_yy_x_zz, g_xzz_yy_x_xx, g_xzz_yy_x_xy, g_xzz_yy_x_xz, g_xzz_yy_x_yy, g_xzz_yy_x_yz, g_xzz_yy_x_zz, g_zz_0_0_0_x_yy_x_xx, g_zz_0_0_0_x_yy_x_xy, g_zz_0_0_0_x_yy_x_xz, g_zz_0_0_0_x_yy_x_yy, g_zz_0_0_0_x_yy_x_yz, g_zz_0_0_0_x_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_yy_x_xx[i] = -2.0 * g_x_yy_x_xx[i] * a_exp + 4.0 * g_xzz_yy_x_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yy_x_xy[i] = -2.0 * g_x_yy_x_xy[i] * a_exp + 4.0 * g_xzz_yy_x_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yy_x_xz[i] = -2.0 * g_x_yy_x_xz[i] * a_exp + 4.0 * g_xzz_yy_x_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yy_x_yy[i] = -2.0 * g_x_yy_x_yy[i] * a_exp + 4.0 * g_xzz_yy_x_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yy_x_yz[i] = -2.0 * g_x_yy_x_yz[i] * a_exp + 4.0 * g_xzz_yy_x_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yy_x_zz[i] = -2.0 * g_x_yy_x_zz[i] * a_exp + 4.0 * g_xzz_yy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1680-1686)

    #pragma omp simd aligned(g_x_yy_y_xx, g_x_yy_y_xy, g_x_yy_y_xz, g_x_yy_y_yy, g_x_yy_y_yz, g_x_yy_y_zz, g_xzz_yy_y_xx, g_xzz_yy_y_xy, g_xzz_yy_y_xz, g_xzz_yy_y_yy, g_xzz_yy_y_yz, g_xzz_yy_y_zz, g_zz_0_0_0_x_yy_y_xx, g_zz_0_0_0_x_yy_y_xy, g_zz_0_0_0_x_yy_y_xz, g_zz_0_0_0_x_yy_y_yy, g_zz_0_0_0_x_yy_y_yz, g_zz_0_0_0_x_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_yy_y_xx[i] = -2.0 * g_x_yy_y_xx[i] * a_exp + 4.0 * g_xzz_yy_y_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yy_y_xy[i] = -2.0 * g_x_yy_y_xy[i] * a_exp + 4.0 * g_xzz_yy_y_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yy_y_xz[i] = -2.0 * g_x_yy_y_xz[i] * a_exp + 4.0 * g_xzz_yy_y_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yy_y_yy[i] = -2.0 * g_x_yy_y_yy[i] * a_exp + 4.0 * g_xzz_yy_y_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yy_y_yz[i] = -2.0 * g_x_yy_y_yz[i] * a_exp + 4.0 * g_xzz_yy_y_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yy_y_zz[i] = -2.0 * g_x_yy_y_zz[i] * a_exp + 4.0 * g_xzz_yy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1686-1692)

    #pragma omp simd aligned(g_x_yy_z_xx, g_x_yy_z_xy, g_x_yy_z_xz, g_x_yy_z_yy, g_x_yy_z_yz, g_x_yy_z_zz, g_xzz_yy_z_xx, g_xzz_yy_z_xy, g_xzz_yy_z_xz, g_xzz_yy_z_yy, g_xzz_yy_z_yz, g_xzz_yy_z_zz, g_zz_0_0_0_x_yy_z_xx, g_zz_0_0_0_x_yy_z_xy, g_zz_0_0_0_x_yy_z_xz, g_zz_0_0_0_x_yy_z_yy, g_zz_0_0_0_x_yy_z_yz, g_zz_0_0_0_x_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_yy_z_xx[i] = -2.0 * g_x_yy_z_xx[i] * a_exp + 4.0 * g_xzz_yy_z_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yy_z_xy[i] = -2.0 * g_x_yy_z_xy[i] * a_exp + 4.0 * g_xzz_yy_z_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yy_z_xz[i] = -2.0 * g_x_yy_z_xz[i] * a_exp + 4.0 * g_xzz_yy_z_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yy_z_yy[i] = -2.0 * g_x_yy_z_yy[i] * a_exp + 4.0 * g_xzz_yy_z_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yy_z_yz[i] = -2.0 * g_x_yy_z_yz[i] * a_exp + 4.0 * g_xzz_yy_z_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yy_z_zz[i] = -2.0 * g_x_yy_z_zz[i] * a_exp + 4.0 * g_xzz_yy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1692-1698)

    #pragma omp simd aligned(g_x_yz_x_xx, g_x_yz_x_xy, g_x_yz_x_xz, g_x_yz_x_yy, g_x_yz_x_yz, g_x_yz_x_zz, g_xzz_yz_x_xx, g_xzz_yz_x_xy, g_xzz_yz_x_xz, g_xzz_yz_x_yy, g_xzz_yz_x_yz, g_xzz_yz_x_zz, g_zz_0_0_0_x_yz_x_xx, g_zz_0_0_0_x_yz_x_xy, g_zz_0_0_0_x_yz_x_xz, g_zz_0_0_0_x_yz_x_yy, g_zz_0_0_0_x_yz_x_yz, g_zz_0_0_0_x_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_yz_x_xx[i] = -2.0 * g_x_yz_x_xx[i] * a_exp + 4.0 * g_xzz_yz_x_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yz_x_xy[i] = -2.0 * g_x_yz_x_xy[i] * a_exp + 4.0 * g_xzz_yz_x_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yz_x_xz[i] = -2.0 * g_x_yz_x_xz[i] * a_exp + 4.0 * g_xzz_yz_x_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yz_x_yy[i] = -2.0 * g_x_yz_x_yy[i] * a_exp + 4.0 * g_xzz_yz_x_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yz_x_yz[i] = -2.0 * g_x_yz_x_yz[i] * a_exp + 4.0 * g_xzz_yz_x_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yz_x_zz[i] = -2.0 * g_x_yz_x_zz[i] * a_exp + 4.0 * g_xzz_yz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1698-1704)

    #pragma omp simd aligned(g_x_yz_y_xx, g_x_yz_y_xy, g_x_yz_y_xz, g_x_yz_y_yy, g_x_yz_y_yz, g_x_yz_y_zz, g_xzz_yz_y_xx, g_xzz_yz_y_xy, g_xzz_yz_y_xz, g_xzz_yz_y_yy, g_xzz_yz_y_yz, g_xzz_yz_y_zz, g_zz_0_0_0_x_yz_y_xx, g_zz_0_0_0_x_yz_y_xy, g_zz_0_0_0_x_yz_y_xz, g_zz_0_0_0_x_yz_y_yy, g_zz_0_0_0_x_yz_y_yz, g_zz_0_0_0_x_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_yz_y_xx[i] = -2.0 * g_x_yz_y_xx[i] * a_exp + 4.0 * g_xzz_yz_y_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yz_y_xy[i] = -2.0 * g_x_yz_y_xy[i] * a_exp + 4.0 * g_xzz_yz_y_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yz_y_xz[i] = -2.0 * g_x_yz_y_xz[i] * a_exp + 4.0 * g_xzz_yz_y_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yz_y_yy[i] = -2.0 * g_x_yz_y_yy[i] * a_exp + 4.0 * g_xzz_yz_y_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yz_y_yz[i] = -2.0 * g_x_yz_y_yz[i] * a_exp + 4.0 * g_xzz_yz_y_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yz_y_zz[i] = -2.0 * g_x_yz_y_zz[i] * a_exp + 4.0 * g_xzz_yz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1704-1710)

    #pragma omp simd aligned(g_x_yz_z_xx, g_x_yz_z_xy, g_x_yz_z_xz, g_x_yz_z_yy, g_x_yz_z_yz, g_x_yz_z_zz, g_xzz_yz_z_xx, g_xzz_yz_z_xy, g_xzz_yz_z_xz, g_xzz_yz_z_yy, g_xzz_yz_z_yz, g_xzz_yz_z_zz, g_zz_0_0_0_x_yz_z_xx, g_zz_0_0_0_x_yz_z_xy, g_zz_0_0_0_x_yz_z_xz, g_zz_0_0_0_x_yz_z_yy, g_zz_0_0_0_x_yz_z_yz, g_zz_0_0_0_x_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_yz_z_xx[i] = -2.0 * g_x_yz_z_xx[i] * a_exp + 4.0 * g_xzz_yz_z_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yz_z_xy[i] = -2.0 * g_x_yz_z_xy[i] * a_exp + 4.0 * g_xzz_yz_z_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yz_z_xz[i] = -2.0 * g_x_yz_z_xz[i] * a_exp + 4.0 * g_xzz_yz_z_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yz_z_yy[i] = -2.0 * g_x_yz_z_yy[i] * a_exp + 4.0 * g_xzz_yz_z_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yz_z_yz[i] = -2.0 * g_x_yz_z_yz[i] * a_exp + 4.0 * g_xzz_yz_z_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yz_z_zz[i] = -2.0 * g_x_yz_z_zz[i] * a_exp + 4.0 * g_xzz_yz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1710-1716)

    #pragma omp simd aligned(g_x_zz_x_xx, g_x_zz_x_xy, g_x_zz_x_xz, g_x_zz_x_yy, g_x_zz_x_yz, g_x_zz_x_zz, g_xzz_zz_x_xx, g_xzz_zz_x_xy, g_xzz_zz_x_xz, g_xzz_zz_x_yy, g_xzz_zz_x_yz, g_xzz_zz_x_zz, g_zz_0_0_0_x_zz_x_xx, g_zz_0_0_0_x_zz_x_xy, g_zz_0_0_0_x_zz_x_xz, g_zz_0_0_0_x_zz_x_yy, g_zz_0_0_0_x_zz_x_yz, g_zz_0_0_0_x_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_zz_x_xx[i] = -2.0 * g_x_zz_x_xx[i] * a_exp + 4.0 * g_xzz_zz_x_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_zz_x_xy[i] = -2.0 * g_x_zz_x_xy[i] * a_exp + 4.0 * g_xzz_zz_x_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_zz_x_xz[i] = -2.0 * g_x_zz_x_xz[i] * a_exp + 4.0 * g_xzz_zz_x_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_zz_x_yy[i] = -2.0 * g_x_zz_x_yy[i] * a_exp + 4.0 * g_xzz_zz_x_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_zz_x_yz[i] = -2.0 * g_x_zz_x_yz[i] * a_exp + 4.0 * g_xzz_zz_x_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_zz_x_zz[i] = -2.0 * g_x_zz_x_zz[i] * a_exp + 4.0 * g_xzz_zz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1716-1722)

    #pragma omp simd aligned(g_x_zz_y_xx, g_x_zz_y_xy, g_x_zz_y_xz, g_x_zz_y_yy, g_x_zz_y_yz, g_x_zz_y_zz, g_xzz_zz_y_xx, g_xzz_zz_y_xy, g_xzz_zz_y_xz, g_xzz_zz_y_yy, g_xzz_zz_y_yz, g_xzz_zz_y_zz, g_zz_0_0_0_x_zz_y_xx, g_zz_0_0_0_x_zz_y_xy, g_zz_0_0_0_x_zz_y_xz, g_zz_0_0_0_x_zz_y_yy, g_zz_0_0_0_x_zz_y_yz, g_zz_0_0_0_x_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_zz_y_xx[i] = -2.0 * g_x_zz_y_xx[i] * a_exp + 4.0 * g_xzz_zz_y_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_zz_y_xy[i] = -2.0 * g_x_zz_y_xy[i] * a_exp + 4.0 * g_xzz_zz_y_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_zz_y_xz[i] = -2.0 * g_x_zz_y_xz[i] * a_exp + 4.0 * g_xzz_zz_y_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_zz_y_yy[i] = -2.0 * g_x_zz_y_yy[i] * a_exp + 4.0 * g_xzz_zz_y_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_zz_y_yz[i] = -2.0 * g_x_zz_y_yz[i] * a_exp + 4.0 * g_xzz_zz_y_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_zz_y_zz[i] = -2.0 * g_x_zz_y_zz[i] * a_exp + 4.0 * g_xzz_zz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1722-1728)

    #pragma omp simd aligned(g_x_zz_z_xx, g_x_zz_z_xy, g_x_zz_z_xz, g_x_zz_z_yy, g_x_zz_z_yz, g_x_zz_z_zz, g_xzz_zz_z_xx, g_xzz_zz_z_xy, g_xzz_zz_z_xz, g_xzz_zz_z_yy, g_xzz_zz_z_yz, g_xzz_zz_z_zz, g_zz_0_0_0_x_zz_z_xx, g_zz_0_0_0_x_zz_z_xy, g_zz_0_0_0_x_zz_z_xz, g_zz_0_0_0_x_zz_z_yy, g_zz_0_0_0_x_zz_z_yz, g_zz_0_0_0_x_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_zz_z_xx[i] = -2.0 * g_x_zz_z_xx[i] * a_exp + 4.0 * g_xzz_zz_z_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_x_zz_z_xy[i] = -2.0 * g_x_zz_z_xy[i] * a_exp + 4.0 * g_xzz_zz_z_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_zz_z_xz[i] = -2.0 * g_x_zz_z_xz[i] * a_exp + 4.0 * g_xzz_zz_z_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_zz_z_yy[i] = -2.0 * g_x_zz_z_yy[i] * a_exp + 4.0 * g_xzz_zz_z_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_x_zz_z_yz[i] = -2.0 * g_x_zz_z_yz[i] * a_exp + 4.0 * g_xzz_zz_z_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_x_zz_z_zz[i] = -2.0 * g_x_zz_z_zz[i] * a_exp + 4.0 * g_xzz_zz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1728-1734)

    #pragma omp simd aligned(g_y_xx_x_xx, g_y_xx_x_xy, g_y_xx_x_xz, g_y_xx_x_yy, g_y_xx_x_yz, g_y_xx_x_zz, g_yzz_xx_x_xx, g_yzz_xx_x_xy, g_yzz_xx_x_xz, g_yzz_xx_x_yy, g_yzz_xx_x_yz, g_yzz_xx_x_zz, g_zz_0_0_0_y_xx_x_xx, g_zz_0_0_0_y_xx_x_xy, g_zz_0_0_0_y_xx_x_xz, g_zz_0_0_0_y_xx_x_yy, g_zz_0_0_0_y_xx_x_yz, g_zz_0_0_0_y_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_xx_x_xx[i] = -2.0 * g_y_xx_x_xx[i] * a_exp + 4.0 * g_yzz_xx_x_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xx_x_xy[i] = -2.0 * g_y_xx_x_xy[i] * a_exp + 4.0 * g_yzz_xx_x_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xx_x_xz[i] = -2.0 * g_y_xx_x_xz[i] * a_exp + 4.0 * g_yzz_xx_x_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xx_x_yy[i] = -2.0 * g_y_xx_x_yy[i] * a_exp + 4.0 * g_yzz_xx_x_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xx_x_yz[i] = -2.0 * g_y_xx_x_yz[i] * a_exp + 4.0 * g_yzz_xx_x_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xx_x_zz[i] = -2.0 * g_y_xx_x_zz[i] * a_exp + 4.0 * g_yzz_xx_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1734-1740)

    #pragma omp simd aligned(g_y_xx_y_xx, g_y_xx_y_xy, g_y_xx_y_xz, g_y_xx_y_yy, g_y_xx_y_yz, g_y_xx_y_zz, g_yzz_xx_y_xx, g_yzz_xx_y_xy, g_yzz_xx_y_xz, g_yzz_xx_y_yy, g_yzz_xx_y_yz, g_yzz_xx_y_zz, g_zz_0_0_0_y_xx_y_xx, g_zz_0_0_0_y_xx_y_xy, g_zz_0_0_0_y_xx_y_xz, g_zz_0_0_0_y_xx_y_yy, g_zz_0_0_0_y_xx_y_yz, g_zz_0_0_0_y_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_xx_y_xx[i] = -2.0 * g_y_xx_y_xx[i] * a_exp + 4.0 * g_yzz_xx_y_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xx_y_xy[i] = -2.0 * g_y_xx_y_xy[i] * a_exp + 4.0 * g_yzz_xx_y_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xx_y_xz[i] = -2.0 * g_y_xx_y_xz[i] * a_exp + 4.0 * g_yzz_xx_y_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xx_y_yy[i] = -2.0 * g_y_xx_y_yy[i] * a_exp + 4.0 * g_yzz_xx_y_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xx_y_yz[i] = -2.0 * g_y_xx_y_yz[i] * a_exp + 4.0 * g_yzz_xx_y_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xx_y_zz[i] = -2.0 * g_y_xx_y_zz[i] * a_exp + 4.0 * g_yzz_xx_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1740-1746)

    #pragma omp simd aligned(g_y_xx_z_xx, g_y_xx_z_xy, g_y_xx_z_xz, g_y_xx_z_yy, g_y_xx_z_yz, g_y_xx_z_zz, g_yzz_xx_z_xx, g_yzz_xx_z_xy, g_yzz_xx_z_xz, g_yzz_xx_z_yy, g_yzz_xx_z_yz, g_yzz_xx_z_zz, g_zz_0_0_0_y_xx_z_xx, g_zz_0_0_0_y_xx_z_xy, g_zz_0_0_0_y_xx_z_xz, g_zz_0_0_0_y_xx_z_yy, g_zz_0_0_0_y_xx_z_yz, g_zz_0_0_0_y_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_xx_z_xx[i] = -2.0 * g_y_xx_z_xx[i] * a_exp + 4.0 * g_yzz_xx_z_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xx_z_xy[i] = -2.0 * g_y_xx_z_xy[i] * a_exp + 4.0 * g_yzz_xx_z_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xx_z_xz[i] = -2.0 * g_y_xx_z_xz[i] * a_exp + 4.0 * g_yzz_xx_z_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xx_z_yy[i] = -2.0 * g_y_xx_z_yy[i] * a_exp + 4.0 * g_yzz_xx_z_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xx_z_yz[i] = -2.0 * g_y_xx_z_yz[i] * a_exp + 4.0 * g_yzz_xx_z_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xx_z_zz[i] = -2.0 * g_y_xx_z_zz[i] * a_exp + 4.0 * g_yzz_xx_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1746-1752)

    #pragma omp simd aligned(g_y_xy_x_xx, g_y_xy_x_xy, g_y_xy_x_xz, g_y_xy_x_yy, g_y_xy_x_yz, g_y_xy_x_zz, g_yzz_xy_x_xx, g_yzz_xy_x_xy, g_yzz_xy_x_xz, g_yzz_xy_x_yy, g_yzz_xy_x_yz, g_yzz_xy_x_zz, g_zz_0_0_0_y_xy_x_xx, g_zz_0_0_0_y_xy_x_xy, g_zz_0_0_0_y_xy_x_xz, g_zz_0_0_0_y_xy_x_yy, g_zz_0_0_0_y_xy_x_yz, g_zz_0_0_0_y_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_xy_x_xx[i] = -2.0 * g_y_xy_x_xx[i] * a_exp + 4.0 * g_yzz_xy_x_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xy_x_xy[i] = -2.0 * g_y_xy_x_xy[i] * a_exp + 4.0 * g_yzz_xy_x_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xy_x_xz[i] = -2.0 * g_y_xy_x_xz[i] * a_exp + 4.0 * g_yzz_xy_x_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xy_x_yy[i] = -2.0 * g_y_xy_x_yy[i] * a_exp + 4.0 * g_yzz_xy_x_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xy_x_yz[i] = -2.0 * g_y_xy_x_yz[i] * a_exp + 4.0 * g_yzz_xy_x_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xy_x_zz[i] = -2.0 * g_y_xy_x_zz[i] * a_exp + 4.0 * g_yzz_xy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1752-1758)

    #pragma omp simd aligned(g_y_xy_y_xx, g_y_xy_y_xy, g_y_xy_y_xz, g_y_xy_y_yy, g_y_xy_y_yz, g_y_xy_y_zz, g_yzz_xy_y_xx, g_yzz_xy_y_xy, g_yzz_xy_y_xz, g_yzz_xy_y_yy, g_yzz_xy_y_yz, g_yzz_xy_y_zz, g_zz_0_0_0_y_xy_y_xx, g_zz_0_0_0_y_xy_y_xy, g_zz_0_0_0_y_xy_y_xz, g_zz_0_0_0_y_xy_y_yy, g_zz_0_0_0_y_xy_y_yz, g_zz_0_0_0_y_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_xy_y_xx[i] = -2.0 * g_y_xy_y_xx[i] * a_exp + 4.0 * g_yzz_xy_y_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xy_y_xy[i] = -2.0 * g_y_xy_y_xy[i] * a_exp + 4.0 * g_yzz_xy_y_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xy_y_xz[i] = -2.0 * g_y_xy_y_xz[i] * a_exp + 4.0 * g_yzz_xy_y_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xy_y_yy[i] = -2.0 * g_y_xy_y_yy[i] * a_exp + 4.0 * g_yzz_xy_y_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xy_y_yz[i] = -2.0 * g_y_xy_y_yz[i] * a_exp + 4.0 * g_yzz_xy_y_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xy_y_zz[i] = -2.0 * g_y_xy_y_zz[i] * a_exp + 4.0 * g_yzz_xy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1758-1764)

    #pragma omp simd aligned(g_y_xy_z_xx, g_y_xy_z_xy, g_y_xy_z_xz, g_y_xy_z_yy, g_y_xy_z_yz, g_y_xy_z_zz, g_yzz_xy_z_xx, g_yzz_xy_z_xy, g_yzz_xy_z_xz, g_yzz_xy_z_yy, g_yzz_xy_z_yz, g_yzz_xy_z_zz, g_zz_0_0_0_y_xy_z_xx, g_zz_0_0_0_y_xy_z_xy, g_zz_0_0_0_y_xy_z_xz, g_zz_0_0_0_y_xy_z_yy, g_zz_0_0_0_y_xy_z_yz, g_zz_0_0_0_y_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_xy_z_xx[i] = -2.0 * g_y_xy_z_xx[i] * a_exp + 4.0 * g_yzz_xy_z_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xy_z_xy[i] = -2.0 * g_y_xy_z_xy[i] * a_exp + 4.0 * g_yzz_xy_z_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xy_z_xz[i] = -2.0 * g_y_xy_z_xz[i] * a_exp + 4.0 * g_yzz_xy_z_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xy_z_yy[i] = -2.0 * g_y_xy_z_yy[i] * a_exp + 4.0 * g_yzz_xy_z_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xy_z_yz[i] = -2.0 * g_y_xy_z_yz[i] * a_exp + 4.0 * g_yzz_xy_z_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xy_z_zz[i] = -2.0 * g_y_xy_z_zz[i] * a_exp + 4.0 * g_yzz_xy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1764-1770)

    #pragma omp simd aligned(g_y_xz_x_xx, g_y_xz_x_xy, g_y_xz_x_xz, g_y_xz_x_yy, g_y_xz_x_yz, g_y_xz_x_zz, g_yzz_xz_x_xx, g_yzz_xz_x_xy, g_yzz_xz_x_xz, g_yzz_xz_x_yy, g_yzz_xz_x_yz, g_yzz_xz_x_zz, g_zz_0_0_0_y_xz_x_xx, g_zz_0_0_0_y_xz_x_xy, g_zz_0_0_0_y_xz_x_xz, g_zz_0_0_0_y_xz_x_yy, g_zz_0_0_0_y_xz_x_yz, g_zz_0_0_0_y_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_xz_x_xx[i] = -2.0 * g_y_xz_x_xx[i] * a_exp + 4.0 * g_yzz_xz_x_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xz_x_xy[i] = -2.0 * g_y_xz_x_xy[i] * a_exp + 4.0 * g_yzz_xz_x_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xz_x_xz[i] = -2.0 * g_y_xz_x_xz[i] * a_exp + 4.0 * g_yzz_xz_x_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xz_x_yy[i] = -2.0 * g_y_xz_x_yy[i] * a_exp + 4.0 * g_yzz_xz_x_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xz_x_yz[i] = -2.0 * g_y_xz_x_yz[i] * a_exp + 4.0 * g_yzz_xz_x_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xz_x_zz[i] = -2.0 * g_y_xz_x_zz[i] * a_exp + 4.0 * g_yzz_xz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1770-1776)

    #pragma omp simd aligned(g_y_xz_y_xx, g_y_xz_y_xy, g_y_xz_y_xz, g_y_xz_y_yy, g_y_xz_y_yz, g_y_xz_y_zz, g_yzz_xz_y_xx, g_yzz_xz_y_xy, g_yzz_xz_y_xz, g_yzz_xz_y_yy, g_yzz_xz_y_yz, g_yzz_xz_y_zz, g_zz_0_0_0_y_xz_y_xx, g_zz_0_0_0_y_xz_y_xy, g_zz_0_0_0_y_xz_y_xz, g_zz_0_0_0_y_xz_y_yy, g_zz_0_0_0_y_xz_y_yz, g_zz_0_0_0_y_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_xz_y_xx[i] = -2.0 * g_y_xz_y_xx[i] * a_exp + 4.0 * g_yzz_xz_y_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xz_y_xy[i] = -2.0 * g_y_xz_y_xy[i] * a_exp + 4.0 * g_yzz_xz_y_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xz_y_xz[i] = -2.0 * g_y_xz_y_xz[i] * a_exp + 4.0 * g_yzz_xz_y_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xz_y_yy[i] = -2.0 * g_y_xz_y_yy[i] * a_exp + 4.0 * g_yzz_xz_y_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xz_y_yz[i] = -2.0 * g_y_xz_y_yz[i] * a_exp + 4.0 * g_yzz_xz_y_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xz_y_zz[i] = -2.0 * g_y_xz_y_zz[i] * a_exp + 4.0 * g_yzz_xz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1776-1782)

    #pragma omp simd aligned(g_y_xz_z_xx, g_y_xz_z_xy, g_y_xz_z_xz, g_y_xz_z_yy, g_y_xz_z_yz, g_y_xz_z_zz, g_yzz_xz_z_xx, g_yzz_xz_z_xy, g_yzz_xz_z_xz, g_yzz_xz_z_yy, g_yzz_xz_z_yz, g_yzz_xz_z_zz, g_zz_0_0_0_y_xz_z_xx, g_zz_0_0_0_y_xz_z_xy, g_zz_0_0_0_y_xz_z_xz, g_zz_0_0_0_y_xz_z_yy, g_zz_0_0_0_y_xz_z_yz, g_zz_0_0_0_y_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_xz_z_xx[i] = -2.0 * g_y_xz_z_xx[i] * a_exp + 4.0 * g_yzz_xz_z_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xz_z_xy[i] = -2.0 * g_y_xz_z_xy[i] * a_exp + 4.0 * g_yzz_xz_z_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xz_z_xz[i] = -2.0 * g_y_xz_z_xz[i] * a_exp + 4.0 * g_yzz_xz_z_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xz_z_yy[i] = -2.0 * g_y_xz_z_yy[i] * a_exp + 4.0 * g_yzz_xz_z_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xz_z_yz[i] = -2.0 * g_y_xz_z_yz[i] * a_exp + 4.0 * g_yzz_xz_z_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xz_z_zz[i] = -2.0 * g_y_xz_z_zz[i] * a_exp + 4.0 * g_yzz_xz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1782-1788)

    #pragma omp simd aligned(g_y_yy_x_xx, g_y_yy_x_xy, g_y_yy_x_xz, g_y_yy_x_yy, g_y_yy_x_yz, g_y_yy_x_zz, g_yzz_yy_x_xx, g_yzz_yy_x_xy, g_yzz_yy_x_xz, g_yzz_yy_x_yy, g_yzz_yy_x_yz, g_yzz_yy_x_zz, g_zz_0_0_0_y_yy_x_xx, g_zz_0_0_0_y_yy_x_xy, g_zz_0_0_0_y_yy_x_xz, g_zz_0_0_0_y_yy_x_yy, g_zz_0_0_0_y_yy_x_yz, g_zz_0_0_0_y_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_yy_x_xx[i] = -2.0 * g_y_yy_x_xx[i] * a_exp + 4.0 * g_yzz_yy_x_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yy_x_xy[i] = -2.0 * g_y_yy_x_xy[i] * a_exp + 4.0 * g_yzz_yy_x_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yy_x_xz[i] = -2.0 * g_y_yy_x_xz[i] * a_exp + 4.0 * g_yzz_yy_x_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yy_x_yy[i] = -2.0 * g_y_yy_x_yy[i] * a_exp + 4.0 * g_yzz_yy_x_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yy_x_yz[i] = -2.0 * g_y_yy_x_yz[i] * a_exp + 4.0 * g_yzz_yy_x_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yy_x_zz[i] = -2.0 * g_y_yy_x_zz[i] * a_exp + 4.0 * g_yzz_yy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1788-1794)

    #pragma omp simd aligned(g_y_yy_y_xx, g_y_yy_y_xy, g_y_yy_y_xz, g_y_yy_y_yy, g_y_yy_y_yz, g_y_yy_y_zz, g_yzz_yy_y_xx, g_yzz_yy_y_xy, g_yzz_yy_y_xz, g_yzz_yy_y_yy, g_yzz_yy_y_yz, g_yzz_yy_y_zz, g_zz_0_0_0_y_yy_y_xx, g_zz_0_0_0_y_yy_y_xy, g_zz_0_0_0_y_yy_y_xz, g_zz_0_0_0_y_yy_y_yy, g_zz_0_0_0_y_yy_y_yz, g_zz_0_0_0_y_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_yy_y_xx[i] = -2.0 * g_y_yy_y_xx[i] * a_exp + 4.0 * g_yzz_yy_y_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yy_y_xy[i] = -2.0 * g_y_yy_y_xy[i] * a_exp + 4.0 * g_yzz_yy_y_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yy_y_xz[i] = -2.0 * g_y_yy_y_xz[i] * a_exp + 4.0 * g_yzz_yy_y_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yy_y_yy[i] = -2.0 * g_y_yy_y_yy[i] * a_exp + 4.0 * g_yzz_yy_y_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yy_y_yz[i] = -2.0 * g_y_yy_y_yz[i] * a_exp + 4.0 * g_yzz_yy_y_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yy_y_zz[i] = -2.0 * g_y_yy_y_zz[i] * a_exp + 4.0 * g_yzz_yy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1794-1800)

    #pragma omp simd aligned(g_y_yy_z_xx, g_y_yy_z_xy, g_y_yy_z_xz, g_y_yy_z_yy, g_y_yy_z_yz, g_y_yy_z_zz, g_yzz_yy_z_xx, g_yzz_yy_z_xy, g_yzz_yy_z_xz, g_yzz_yy_z_yy, g_yzz_yy_z_yz, g_yzz_yy_z_zz, g_zz_0_0_0_y_yy_z_xx, g_zz_0_0_0_y_yy_z_xy, g_zz_0_0_0_y_yy_z_xz, g_zz_0_0_0_y_yy_z_yy, g_zz_0_0_0_y_yy_z_yz, g_zz_0_0_0_y_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_yy_z_xx[i] = -2.0 * g_y_yy_z_xx[i] * a_exp + 4.0 * g_yzz_yy_z_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yy_z_xy[i] = -2.0 * g_y_yy_z_xy[i] * a_exp + 4.0 * g_yzz_yy_z_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yy_z_xz[i] = -2.0 * g_y_yy_z_xz[i] * a_exp + 4.0 * g_yzz_yy_z_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yy_z_yy[i] = -2.0 * g_y_yy_z_yy[i] * a_exp + 4.0 * g_yzz_yy_z_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yy_z_yz[i] = -2.0 * g_y_yy_z_yz[i] * a_exp + 4.0 * g_yzz_yy_z_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yy_z_zz[i] = -2.0 * g_y_yy_z_zz[i] * a_exp + 4.0 * g_yzz_yy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1800-1806)

    #pragma omp simd aligned(g_y_yz_x_xx, g_y_yz_x_xy, g_y_yz_x_xz, g_y_yz_x_yy, g_y_yz_x_yz, g_y_yz_x_zz, g_yzz_yz_x_xx, g_yzz_yz_x_xy, g_yzz_yz_x_xz, g_yzz_yz_x_yy, g_yzz_yz_x_yz, g_yzz_yz_x_zz, g_zz_0_0_0_y_yz_x_xx, g_zz_0_0_0_y_yz_x_xy, g_zz_0_0_0_y_yz_x_xz, g_zz_0_0_0_y_yz_x_yy, g_zz_0_0_0_y_yz_x_yz, g_zz_0_0_0_y_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_yz_x_xx[i] = -2.0 * g_y_yz_x_xx[i] * a_exp + 4.0 * g_yzz_yz_x_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yz_x_xy[i] = -2.0 * g_y_yz_x_xy[i] * a_exp + 4.0 * g_yzz_yz_x_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yz_x_xz[i] = -2.0 * g_y_yz_x_xz[i] * a_exp + 4.0 * g_yzz_yz_x_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yz_x_yy[i] = -2.0 * g_y_yz_x_yy[i] * a_exp + 4.0 * g_yzz_yz_x_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yz_x_yz[i] = -2.0 * g_y_yz_x_yz[i] * a_exp + 4.0 * g_yzz_yz_x_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yz_x_zz[i] = -2.0 * g_y_yz_x_zz[i] * a_exp + 4.0 * g_yzz_yz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1806-1812)

    #pragma omp simd aligned(g_y_yz_y_xx, g_y_yz_y_xy, g_y_yz_y_xz, g_y_yz_y_yy, g_y_yz_y_yz, g_y_yz_y_zz, g_yzz_yz_y_xx, g_yzz_yz_y_xy, g_yzz_yz_y_xz, g_yzz_yz_y_yy, g_yzz_yz_y_yz, g_yzz_yz_y_zz, g_zz_0_0_0_y_yz_y_xx, g_zz_0_0_0_y_yz_y_xy, g_zz_0_0_0_y_yz_y_xz, g_zz_0_0_0_y_yz_y_yy, g_zz_0_0_0_y_yz_y_yz, g_zz_0_0_0_y_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_yz_y_xx[i] = -2.0 * g_y_yz_y_xx[i] * a_exp + 4.0 * g_yzz_yz_y_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yz_y_xy[i] = -2.0 * g_y_yz_y_xy[i] * a_exp + 4.0 * g_yzz_yz_y_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yz_y_xz[i] = -2.0 * g_y_yz_y_xz[i] * a_exp + 4.0 * g_yzz_yz_y_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yz_y_yy[i] = -2.0 * g_y_yz_y_yy[i] * a_exp + 4.0 * g_yzz_yz_y_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yz_y_yz[i] = -2.0 * g_y_yz_y_yz[i] * a_exp + 4.0 * g_yzz_yz_y_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yz_y_zz[i] = -2.0 * g_y_yz_y_zz[i] * a_exp + 4.0 * g_yzz_yz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1812-1818)

    #pragma omp simd aligned(g_y_yz_z_xx, g_y_yz_z_xy, g_y_yz_z_xz, g_y_yz_z_yy, g_y_yz_z_yz, g_y_yz_z_zz, g_yzz_yz_z_xx, g_yzz_yz_z_xy, g_yzz_yz_z_xz, g_yzz_yz_z_yy, g_yzz_yz_z_yz, g_yzz_yz_z_zz, g_zz_0_0_0_y_yz_z_xx, g_zz_0_0_0_y_yz_z_xy, g_zz_0_0_0_y_yz_z_xz, g_zz_0_0_0_y_yz_z_yy, g_zz_0_0_0_y_yz_z_yz, g_zz_0_0_0_y_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_yz_z_xx[i] = -2.0 * g_y_yz_z_xx[i] * a_exp + 4.0 * g_yzz_yz_z_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yz_z_xy[i] = -2.0 * g_y_yz_z_xy[i] * a_exp + 4.0 * g_yzz_yz_z_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yz_z_xz[i] = -2.0 * g_y_yz_z_xz[i] * a_exp + 4.0 * g_yzz_yz_z_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yz_z_yy[i] = -2.0 * g_y_yz_z_yy[i] * a_exp + 4.0 * g_yzz_yz_z_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yz_z_yz[i] = -2.0 * g_y_yz_z_yz[i] * a_exp + 4.0 * g_yzz_yz_z_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yz_z_zz[i] = -2.0 * g_y_yz_z_zz[i] * a_exp + 4.0 * g_yzz_yz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1818-1824)

    #pragma omp simd aligned(g_y_zz_x_xx, g_y_zz_x_xy, g_y_zz_x_xz, g_y_zz_x_yy, g_y_zz_x_yz, g_y_zz_x_zz, g_yzz_zz_x_xx, g_yzz_zz_x_xy, g_yzz_zz_x_xz, g_yzz_zz_x_yy, g_yzz_zz_x_yz, g_yzz_zz_x_zz, g_zz_0_0_0_y_zz_x_xx, g_zz_0_0_0_y_zz_x_xy, g_zz_0_0_0_y_zz_x_xz, g_zz_0_0_0_y_zz_x_yy, g_zz_0_0_0_y_zz_x_yz, g_zz_0_0_0_y_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_zz_x_xx[i] = -2.0 * g_y_zz_x_xx[i] * a_exp + 4.0 * g_yzz_zz_x_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_zz_x_xy[i] = -2.0 * g_y_zz_x_xy[i] * a_exp + 4.0 * g_yzz_zz_x_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_zz_x_xz[i] = -2.0 * g_y_zz_x_xz[i] * a_exp + 4.0 * g_yzz_zz_x_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_zz_x_yy[i] = -2.0 * g_y_zz_x_yy[i] * a_exp + 4.0 * g_yzz_zz_x_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_zz_x_yz[i] = -2.0 * g_y_zz_x_yz[i] * a_exp + 4.0 * g_yzz_zz_x_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_zz_x_zz[i] = -2.0 * g_y_zz_x_zz[i] * a_exp + 4.0 * g_yzz_zz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1824-1830)

    #pragma omp simd aligned(g_y_zz_y_xx, g_y_zz_y_xy, g_y_zz_y_xz, g_y_zz_y_yy, g_y_zz_y_yz, g_y_zz_y_zz, g_yzz_zz_y_xx, g_yzz_zz_y_xy, g_yzz_zz_y_xz, g_yzz_zz_y_yy, g_yzz_zz_y_yz, g_yzz_zz_y_zz, g_zz_0_0_0_y_zz_y_xx, g_zz_0_0_0_y_zz_y_xy, g_zz_0_0_0_y_zz_y_xz, g_zz_0_0_0_y_zz_y_yy, g_zz_0_0_0_y_zz_y_yz, g_zz_0_0_0_y_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_zz_y_xx[i] = -2.0 * g_y_zz_y_xx[i] * a_exp + 4.0 * g_yzz_zz_y_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_zz_y_xy[i] = -2.0 * g_y_zz_y_xy[i] * a_exp + 4.0 * g_yzz_zz_y_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_zz_y_xz[i] = -2.0 * g_y_zz_y_xz[i] * a_exp + 4.0 * g_yzz_zz_y_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_zz_y_yy[i] = -2.0 * g_y_zz_y_yy[i] * a_exp + 4.0 * g_yzz_zz_y_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_zz_y_yz[i] = -2.0 * g_y_zz_y_yz[i] * a_exp + 4.0 * g_yzz_zz_y_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_zz_y_zz[i] = -2.0 * g_y_zz_y_zz[i] * a_exp + 4.0 * g_yzz_zz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1830-1836)

    #pragma omp simd aligned(g_y_zz_z_xx, g_y_zz_z_xy, g_y_zz_z_xz, g_y_zz_z_yy, g_y_zz_z_yz, g_y_zz_z_zz, g_yzz_zz_z_xx, g_yzz_zz_z_xy, g_yzz_zz_z_xz, g_yzz_zz_z_yy, g_yzz_zz_z_yz, g_yzz_zz_z_zz, g_zz_0_0_0_y_zz_z_xx, g_zz_0_0_0_y_zz_z_xy, g_zz_0_0_0_y_zz_z_xz, g_zz_0_0_0_y_zz_z_yy, g_zz_0_0_0_y_zz_z_yz, g_zz_0_0_0_y_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_zz_z_xx[i] = -2.0 * g_y_zz_z_xx[i] * a_exp + 4.0 * g_yzz_zz_z_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_y_zz_z_xy[i] = -2.0 * g_y_zz_z_xy[i] * a_exp + 4.0 * g_yzz_zz_z_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_zz_z_xz[i] = -2.0 * g_y_zz_z_xz[i] * a_exp + 4.0 * g_yzz_zz_z_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_zz_z_yy[i] = -2.0 * g_y_zz_z_yy[i] * a_exp + 4.0 * g_yzz_zz_z_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_y_zz_z_yz[i] = -2.0 * g_y_zz_z_yz[i] * a_exp + 4.0 * g_yzz_zz_z_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_y_zz_z_zz[i] = -2.0 * g_y_zz_z_zz[i] * a_exp + 4.0 * g_yzz_zz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1836-1842)

    #pragma omp simd aligned(g_z_xx_x_xx, g_z_xx_x_xy, g_z_xx_x_xz, g_z_xx_x_yy, g_z_xx_x_yz, g_z_xx_x_zz, g_zz_0_0_0_z_xx_x_xx, g_zz_0_0_0_z_xx_x_xy, g_zz_0_0_0_z_xx_x_xz, g_zz_0_0_0_z_xx_x_yy, g_zz_0_0_0_z_xx_x_yz, g_zz_0_0_0_z_xx_x_zz, g_zzz_xx_x_xx, g_zzz_xx_x_xy, g_zzz_xx_x_xz, g_zzz_xx_x_yy, g_zzz_xx_x_yz, g_zzz_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_xx_x_xx[i] = -6.0 * g_z_xx_x_xx[i] * a_exp + 4.0 * g_zzz_xx_x_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xx_x_xy[i] = -6.0 * g_z_xx_x_xy[i] * a_exp + 4.0 * g_zzz_xx_x_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xx_x_xz[i] = -6.0 * g_z_xx_x_xz[i] * a_exp + 4.0 * g_zzz_xx_x_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xx_x_yy[i] = -6.0 * g_z_xx_x_yy[i] * a_exp + 4.0 * g_zzz_xx_x_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xx_x_yz[i] = -6.0 * g_z_xx_x_yz[i] * a_exp + 4.0 * g_zzz_xx_x_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xx_x_zz[i] = -6.0 * g_z_xx_x_zz[i] * a_exp + 4.0 * g_zzz_xx_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1842-1848)

    #pragma omp simd aligned(g_z_xx_y_xx, g_z_xx_y_xy, g_z_xx_y_xz, g_z_xx_y_yy, g_z_xx_y_yz, g_z_xx_y_zz, g_zz_0_0_0_z_xx_y_xx, g_zz_0_0_0_z_xx_y_xy, g_zz_0_0_0_z_xx_y_xz, g_zz_0_0_0_z_xx_y_yy, g_zz_0_0_0_z_xx_y_yz, g_zz_0_0_0_z_xx_y_zz, g_zzz_xx_y_xx, g_zzz_xx_y_xy, g_zzz_xx_y_xz, g_zzz_xx_y_yy, g_zzz_xx_y_yz, g_zzz_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_xx_y_xx[i] = -6.0 * g_z_xx_y_xx[i] * a_exp + 4.0 * g_zzz_xx_y_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xx_y_xy[i] = -6.0 * g_z_xx_y_xy[i] * a_exp + 4.0 * g_zzz_xx_y_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xx_y_xz[i] = -6.0 * g_z_xx_y_xz[i] * a_exp + 4.0 * g_zzz_xx_y_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xx_y_yy[i] = -6.0 * g_z_xx_y_yy[i] * a_exp + 4.0 * g_zzz_xx_y_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xx_y_yz[i] = -6.0 * g_z_xx_y_yz[i] * a_exp + 4.0 * g_zzz_xx_y_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xx_y_zz[i] = -6.0 * g_z_xx_y_zz[i] * a_exp + 4.0 * g_zzz_xx_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1848-1854)

    #pragma omp simd aligned(g_z_xx_z_xx, g_z_xx_z_xy, g_z_xx_z_xz, g_z_xx_z_yy, g_z_xx_z_yz, g_z_xx_z_zz, g_zz_0_0_0_z_xx_z_xx, g_zz_0_0_0_z_xx_z_xy, g_zz_0_0_0_z_xx_z_xz, g_zz_0_0_0_z_xx_z_yy, g_zz_0_0_0_z_xx_z_yz, g_zz_0_0_0_z_xx_z_zz, g_zzz_xx_z_xx, g_zzz_xx_z_xy, g_zzz_xx_z_xz, g_zzz_xx_z_yy, g_zzz_xx_z_yz, g_zzz_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_xx_z_xx[i] = -6.0 * g_z_xx_z_xx[i] * a_exp + 4.0 * g_zzz_xx_z_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xx_z_xy[i] = -6.0 * g_z_xx_z_xy[i] * a_exp + 4.0 * g_zzz_xx_z_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xx_z_xz[i] = -6.0 * g_z_xx_z_xz[i] * a_exp + 4.0 * g_zzz_xx_z_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xx_z_yy[i] = -6.0 * g_z_xx_z_yy[i] * a_exp + 4.0 * g_zzz_xx_z_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xx_z_yz[i] = -6.0 * g_z_xx_z_yz[i] * a_exp + 4.0 * g_zzz_xx_z_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xx_z_zz[i] = -6.0 * g_z_xx_z_zz[i] * a_exp + 4.0 * g_zzz_xx_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1854-1860)

    #pragma omp simd aligned(g_z_xy_x_xx, g_z_xy_x_xy, g_z_xy_x_xz, g_z_xy_x_yy, g_z_xy_x_yz, g_z_xy_x_zz, g_zz_0_0_0_z_xy_x_xx, g_zz_0_0_0_z_xy_x_xy, g_zz_0_0_0_z_xy_x_xz, g_zz_0_0_0_z_xy_x_yy, g_zz_0_0_0_z_xy_x_yz, g_zz_0_0_0_z_xy_x_zz, g_zzz_xy_x_xx, g_zzz_xy_x_xy, g_zzz_xy_x_xz, g_zzz_xy_x_yy, g_zzz_xy_x_yz, g_zzz_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_xy_x_xx[i] = -6.0 * g_z_xy_x_xx[i] * a_exp + 4.0 * g_zzz_xy_x_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xy_x_xy[i] = -6.0 * g_z_xy_x_xy[i] * a_exp + 4.0 * g_zzz_xy_x_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xy_x_xz[i] = -6.0 * g_z_xy_x_xz[i] * a_exp + 4.0 * g_zzz_xy_x_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xy_x_yy[i] = -6.0 * g_z_xy_x_yy[i] * a_exp + 4.0 * g_zzz_xy_x_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xy_x_yz[i] = -6.0 * g_z_xy_x_yz[i] * a_exp + 4.0 * g_zzz_xy_x_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xy_x_zz[i] = -6.0 * g_z_xy_x_zz[i] * a_exp + 4.0 * g_zzz_xy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1860-1866)

    #pragma omp simd aligned(g_z_xy_y_xx, g_z_xy_y_xy, g_z_xy_y_xz, g_z_xy_y_yy, g_z_xy_y_yz, g_z_xy_y_zz, g_zz_0_0_0_z_xy_y_xx, g_zz_0_0_0_z_xy_y_xy, g_zz_0_0_0_z_xy_y_xz, g_zz_0_0_0_z_xy_y_yy, g_zz_0_0_0_z_xy_y_yz, g_zz_0_0_0_z_xy_y_zz, g_zzz_xy_y_xx, g_zzz_xy_y_xy, g_zzz_xy_y_xz, g_zzz_xy_y_yy, g_zzz_xy_y_yz, g_zzz_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_xy_y_xx[i] = -6.0 * g_z_xy_y_xx[i] * a_exp + 4.0 * g_zzz_xy_y_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xy_y_xy[i] = -6.0 * g_z_xy_y_xy[i] * a_exp + 4.0 * g_zzz_xy_y_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xy_y_xz[i] = -6.0 * g_z_xy_y_xz[i] * a_exp + 4.0 * g_zzz_xy_y_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xy_y_yy[i] = -6.0 * g_z_xy_y_yy[i] * a_exp + 4.0 * g_zzz_xy_y_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xy_y_yz[i] = -6.0 * g_z_xy_y_yz[i] * a_exp + 4.0 * g_zzz_xy_y_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xy_y_zz[i] = -6.0 * g_z_xy_y_zz[i] * a_exp + 4.0 * g_zzz_xy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1866-1872)

    #pragma omp simd aligned(g_z_xy_z_xx, g_z_xy_z_xy, g_z_xy_z_xz, g_z_xy_z_yy, g_z_xy_z_yz, g_z_xy_z_zz, g_zz_0_0_0_z_xy_z_xx, g_zz_0_0_0_z_xy_z_xy, g_zz_0_0_0_z_xy_z_xz, g_zz_0_0_0_z_xy_z_yy, g_zz_0_0_0_z_xy_z_yz, g_zz_0_0_0_z_xy_z_zz, g_zzz_xy_z_xx, g_zzz_xy_z_xy, g_zzz_xy_z_xz, g_zzz_xy_z_yy, g_zzz_xy_z_yz, g_zzz_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_xy_z_xx[i] = -6.0 * g_z_xy_z_xx[i] * a_exp + 4.0 * g_zzz_xy_z_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xy_z_xy[i] = -6.0 * g_z_xy_z_xy[i] * a_exp + 4.0 * g_zzz_xy_z_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xy_z_xz[i] = -6.0 * g_z_xy_z_xz[i] * a_exp + 4.0 * g_zzz_xy_z_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xy_z_yy[i] = -6.0 * g_z_xy_z_yy[i] * a_exp + 4.0 * g_zzz_xy_z_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xy_z_yz[i] = -6.0 * g_z_xy_z_yz[i] * a_exp + 4.0 * g_zzz_xy_z_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xy_z_zz[i] = -6.0 * g_z_xy_z_zz[i] * a_exp + 4.0 * g_zzz_xy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1872-1878)

    #pragma omp simd aligned(g_z_xz_x_xx, g_z_xz_x_xy, g_z_xz_x_xz, g_z_xz_x_yy, g_z_xz_x_yz, g_z_xz_x_zz, g_zz_0_0_0_z_xz_x_xx, g_zz_0_0_0_z_xz_x_xy, g_zz_0_0_0_z_xz_x_xz, g_zz_0_0_0_z_xz_x_yy, g_zz_0_0_0_z_xz_x_yz, g_zz_0_0_0_z_xz_x_zz, g_zzz_xz_x_xx, g_zzz_xz_x_xy, g_zzz_xz_x_xz, g_zzz_xz_x_yy, g_zzz_xz_x_yz, g_zzz_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_xz_x_xx[i] = -6.0 * g_z_xz_x_xx[i] * a_exp + 4.0 * g_zzz_xz_x_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xz_x_xy[i] = -6.0 * g_z_xz_x_xy[i] * a_exp + 4.0 * g_zzz_xz_x_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xz_x_xz[i] = -6.0 * g_z_xz_x_xz[i] * a_exp + 4.0 * g_zzz_xz_x_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xz_x_yy[i] = -6.0 * g_z_xz_x_yy[i] * a_exp + 4.0 * g_zzz_xz_x_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xz_x_yz[i] = -6.0 * g_z_xz_x_yz[i] * a_exp + 4.0 * g_zzz_xz_x_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xz_x_zz[i] = -6.0 * g_z_xz_x_zz[i] * a_exp + 4.0 * g_zzz_xz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1878-1884)

    #pragma omp simd aligned(g_z_xz_y_xx, g_z_xz_y_xy, g_z_xz_y_xz, g_z_xz_y_yy, g_z_xz_y_yz, g_z_xz_y_zz, g_zz_0_0_0_z_xz_y_xx, g_zz_0_0_0_z_xz_y_xy, g_zz_0_0_0_z_xz_y_xz, g_zz_0_0_0_z_xz_y_yy, g_zz_0_0_0_z_xz_y_yz, g_zz_0_0_0_z_xz_y_zz, g_zzz_xz_y_xx, g_zzz_xz_y_xy, g_zzz_xz_y_xz, g_zzz_xz_y_yy, g_zzz_xz_y_yz, g_zzz_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_xz_y_xx[i] = -6.0 * g_z_xz_y_xx[i] * a_exp + 4.0 * g_zzz_xz_y_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xz_y_xy[i] = -6.0 * g_z_xz_y_xy[i] * a_exp + 4.0 * g_zzz_xz_y_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xz_y_xz[i] = -6.0 * g_z_xz_y_xz[i] * a_exp + 4.0 * g_zzz_xz_y_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xz_y_yy[i] = -6.0 * g_z_xz_y_yy[i] * a_exp + 4.0 * g_zzz_xz_y_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xz_y_yz[i] = -6.0 * g_z_xz_y_yz[i] * a_exp + 4.0 * g_zzz_xz_y_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xz_y_zz[i] = -6.0 * g_z_xz_y_zz[i] * a_exp + 4.0 * g_zzz_xz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1884-1890)

    #pragma omp simd aligned(g_z_xz_z_xx, g_z_xz_z_xy, g_z_xz_z_xz, g_z_xz_z_yy, g_z_xz_z_yz, g_z_xz_z_zz, g_zz_0_0_0_z_xz_z_xx, g_zz_0_0_0_z_xz_z_xy, g_zz_0_0_0_z_xz_z_xz, g_zz_0_0_0_z_xz_z_yy, g_zz_0_0_0_z_xz_z_yz, g_zz_0_0_0_z_xz_z_zz, g_zzz_xz_z_xx, g_zzz_xz_z_xy, g_zzz_xz_z_xz, g_zzz_xz_z_yy, g_zzz_xz_z_yz, g_zzz_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_xz_z_xx[i] = -6.0 * g_z_xz_z_xx[i] * a_exp + 4.0 * g_zzz_xz_z_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xz_z_xy[i] = -6.0 * g_z_xz_z_xy[i] * a_exp + 4.0 * g_zzz_xz_z_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xz_z_xz[i] = -6.0 * g_z_xz_z_xz[i] * a_exp + 4.0 * g_zzz_xz_z_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xz_z_yy[i] = -6.0 * g_z_xz_z_yy[i] * a_exp + 4.0 * g_zzz_xz_z_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xz_z_yz[i] = -6.0 * g_z_xz_z_yz[i] * a_exp + 4.0 * g_zzz_xz_z_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xz_z_zz[i] = -6.0 * g_z_xz_z_zz[i] * a_exp + 4.0 * g_zzz_xz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1890-1896)

    #pragma omp simd aligned(g_z_yy_x_xx, g_z_yy_x_xy, g_z_yy_x_xz, g_z_yy_x_yy, g_z_yy_x_yz, g_z_yy_x_zz, g_zz_0_0_0_z_yy_x_xx, g_zz_0_0_0_z_yy_x_xy, g_zz_0_0_0_z_yy_x_xz, g_zz_0_0_0_z_yy_x_yy, g_zz_0_0_0_z_yy_x_yz, g_zz_0_0_0_z_yy_x_zz, g_zzz_yy_x_xx, g_zzz_yy_x_xy, g_zzz_yy_x_xz, g_zzz_yy_x_yy, g_zzz_yy_x_yz, g_zzz_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_yy_x_xx[i] = -6.0 * g_z_yy_x_xx[i] * a_exp + 4.0 * g_zzz_yy_x_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yy_x_xy[i] = -6.0 * g_z_yy_x_xy[i] * a_exp + 4.0 * g_zzz_yy_x_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yy_x_xz[i] = -6.0 * g_z_yy_x_xz[i] * a_exp + 4.0 * g_zzz_yy_x_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yy_x_yy[i] = -6.0 * g_z_yy_x_yy[i] * a_exp + 4.0 * g_zzz_yy_x_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yy_x_yz[i] = -6.0 * g_z_yy_x_yz[i] * a_exp + 4.0 * g_zzz_yy_x_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yy_x_zz[i] = -6.0 * g_z_yy_x_zz[i] * a_exp + 4.0 * g_zzz_yy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1896-1902)

    #pragma omp simd aligned(g_z_yy_y_xx, g_z_yy_y_xy, g_z_yy_y_xz, g_z_yy_y_yy, g_z_yy_y_yz, g_z_yy_y_zz, g_zz_0_0_0_z_yy_y_xx, g_zz_0_0_0_z_yy_y_xy, g_zz_0_0_0_z_yy_y_xz, g_zz_0_0_0_z_yy_y_yy, g_zz_0_0_0_z_yy_y_yz, g_zz_0_0_0_z_yy_y_zz, g_zzz_yy_y_xx, g_zzz_yy_y_xy, g_zzz_yy_y_xz, g_zzz_yy_y_yy, g_zzz_yy_y_yz, g_zzz_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_yy_y_xx[i] = -6.0 * g_z_yy_y_xx[i] * a_exp + 4.0 * g_zzz_yy_y_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yy_y_xy[i] = -6.0 * g_z_yy_y_xy[i] * a_exp + 4.0 * g_zzz_yy_y_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yy_y_xz[i] = -6.0 * g_z_yy_y_xz[i] * a_exp + 4.0 * g_zzz_yy_y_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yy_y_yy[i] = -6.0 * g_z_yy_y_yy[i] * a_exp + 4.0 * g_zzz_yy_y_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yy_y_yz[i] = -6.0 * g_z_yy_y_yz[i] * a_exp + 4.0 * g_zzz_yy_y_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yy_y_zz[i] = -6.0 * g_z_yy_y_zz[i] * a_exp + 4.0 * g_zzz_yy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1902-1908)

    #pragma omp simd aligned(g_z_yy_z_xx, g_z_yy_z_xy, g_z_yy_z_xz, g_z_yy_z_yy, g_z_yy_z_yz, g_z_yy_z_zz, g_zz_0_0_0_z_yy_z_xx, g_zz_0_0_0_z_yy_z_xy, g_zz_0_0_0_z_yy_z_xz, g_zz_0_0_0_z_yy_z_yy, g_zz_0_0_0_z_yy_z_yz, g_zz_0_0_0_z_yy_z_zz, g_zzz_yy_z_xx, g_zzz_yy_z_xy, g_zzz_yy_z_xz, g_zzz_yy_z_yy, g_zzz_yy_z_yz, g_zzz_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_yy_z_xx[i] = -6.0 * g_z_yy_z_xx[i] * a_exp + 4.0 * g_zzz_yy_z_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yy_z_xy[i] = -6.0 * g_z_yy_z_xy[i] * a_exp + 4.0 * g_zzz_yy_z_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yy_z_xz[i] = -6.0 * g_z_yy_z_xz[i] * a_exp + 4.0 * g_zzz_yy_z_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yy_z_yy[i] = -6.0 * g_z_yy_z_yy[i] * a_exp + 4.0 * g_zzz_yy_z_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yy_z_yz[i] = -6.0 * g_z_yy_z_yz[i] * a_exp + 4.0 * g_zzz_yy_z_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yy_z_zz[i] = -6.0 * g_z_yy_z_zz[i] * a_exp + 4.0 * g_zzz_yy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1908-1914)

    #pragma omp simd aligned(g_z_yz_x_xx, g_z_yz_x_xy, g_z_yz_x_xz, g_z_yz_x_yy, g_z_yz_x_yz, g_z_yz_x_zz, g_zz_0_0_0_z_yz_x_xx, g_zz_0_0_0_z_yz_x_xy, g_zz_0_0_0_z_yz_x_xz, g_zz_0_0_0_z_yz_x_yy, g_zz_0_0_0_z_yz_x_yz, g_zz_0_0_0_z_yz_x_zz, g_zzz_yz_x_xx, g_zzz_yz_x_xy, g_zzz_yz_x_xz, g_zzz_yz_x_yy, g_zzz_yz_x_yz, g_zzz_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_yz_x_xx[i] = -6.0 * g_z_yz_x_xx[i] * a_exp + 4.0 * g_zzz_yz_x_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yz_x_xy[i] = -6.0 * g_z_yz_x_xy[i] * a_exp + 4.0 * g_zzz_yz_x_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yz_x_xz[i] = -6.0 * g_z_yz_x_xz[i] * a_exp + 4.0 * g_zzz_yz_x_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yz_x_yy[i] = -6.0 * g_z_yz_x_yy[i] * a_exp + 4.0 * g_zzz_yz_x_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yz_x_yz[i] = -6.0 * g_z_yz_x_yz[i] * a_exp + 4.0 * g_zzz_yz_x_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yz_x_zz[i] = -6.0 * g_z_yz_x_zz[i] * a_exp + 4.0 * g_zzz_yz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1914-1920)

    #pragma omp simd aligned(g_z_yz_y_xx, g_z_yz_y_xy, g_z_yz_y_xz, g_z_yz_y_yy, g_z_yz_y_yz, g_z_yz_y_zz, g_zz_0_0_0_z_yz_y_xx, g_zz_0_0_0_z_yz_y_xy, g_zz_0_0_0_z_yz_y_xz, g_zz_0_0_0_z_yz_y_yy, g_zz_0_0_0_z_yz_y_yz, g_zz_0_0_0_z_yz_y_zz, g_zzz_yz_y_xx, g_zzz_yz_y_xy, g_zzz_yz_y_xz, g_zzz_yz_y_yy, g_zzz_yz_y_yz, g_zzz_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_yz_y_xx[i] = -6.0 * g_z_yz_y_xx[i] * a_exp + 4.0 * g_zzz_yz_y_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yz_y_xy[i] = -6.0 * g_z_yz_y_xy[i] * a_exp + 4.0 * g_zzz_yz_y_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yz_y_xz[i] = -6.0 * g_z_yz_y_xz[i] * a_exp + 4.0 * g_zzz_yz_y_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yz_y_yy[i] = -6.0 * g_z_yz_y_yy[i] * a_exp + 4.0 * g_zzz_yz_y_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yz_y_yz[i] = -6.0 * g_z_yz_y_yz[i] * a_exp + 4.0 * g_zzz_yz_y_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yz_y_zz[i] = -6.0 * g_z_yz_y_zz[i] * a_exp + 4.0 * g_zzz_yz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1920-1926)

    #pragma omp simd aligned(g_z_yz_z_xx, g_z_yz_z_xy, g_z_yz_z_xz, g_z_yz_z_yy, g_z_yz_z_yz, g_z_yz_z_zz, g_zz_0_0_0_z_yz_z_xx, g_zz_0_0_0_z_yz_z_xy, g_zz_0_0_0_z_yz_z_xz, g_zz_0_0_0_z_yz_z_yy, g_zz_0_0_0_z_yz_z_yz, g_zz_0_0_0_z_yz_z_zz, g_zzz_yz_z_xx, g_zzz_yz_z_xy, g_zzz_yz_z_xz, g_zzz_yz_z_yy, g_zzz_yz_z_yz, g_zzz_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_yz_z_xx[i] = -6.0 * g_z_yz_z_xx[i] * a_exp + 4.0 * g_zzz_yz_z_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yz_z_xy[i] = -6.0 * g_z_yz_z_xy[i] * a_exp + 4.0 * g_zzz_yz_z_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yz_z_xz[i] = -6.0 * g_z_yz_z_xz[i] * a_exp + 4.0 * g_zzz_yz_z_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yz_z_yy[i] = -6.0 * g_z_yz_z_yy[i] * a_exp + 4.0 * g_zzz_yz_z_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yz_z_yz[i] = -6.0 * g_z_yz_z_yz[i] * a_exp + 4.0 * g_zzz_yz_z_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yz_z_zz[i] = -6.0 * g_z_yz_z_zz[i] * a_exp + 4.0 * g_zzz_yz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (1926-1932)

    #pragma omp simd aligned(g_z_zz_x_xx, g_z_zz_x_xy, g_z_zz_x_xz, g_z_zz_x_yy, g_z_zz_x_yz, g_z_zz_x_zz, g_zz_0_0_0_z_zz_x_xx, g_zz_0_0_0_z_zz_x_xy, g_zz_0_0_0_z_zz_x_xz, g_zz_0_0_0_z_zz_x_yy, g_zz_0_0_0_z_zz_x_yz, g_zz_0_0_0_z_zz_x_zz, g_zzz_zz_x_xx, g_zzz_zz_x_xy, g_zzz_zz_x_xz, g_zzz_zz_x_yy, g_zzz_zz_x_yz, g_zzz_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_zz_x_xx[i] = -6.0 * g_z_zz_x_xx[i] * a_exp + 4.0 * g_zzz_zz_x_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_zz_x_xy[i] = -6.0 * g_z_zz_x_xy[i] * a_exp + 4.0 * g_zzz_zz_x_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_zz_x_xz[i] = -6.0 * g_z_zz_x_xz[i] * a_exp + 4.0 * g_zzz_zz_x_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_zz_x_yy[i] = -6.0 * g_z_zz_x_yy[i] * a_exp + 4.0 * g_zzz_zz_x_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_zz_x_yz[i] = -6.0 * g_z_zz_x_yz[i] * a_exp + 4.0 * g_zzz_zz_x_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_zz_x_zz[i] = -6.0 * g_z_zz_x_zz[i] * a_exp + 4.0 * g_zzz_zz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (1932-1938)

    #pragma omp simd aligned(g_z_zz_y_xx, g_z_zz_y_xy, g_z_zz_y_xz, g_z_zz_y_yy, g_z_zz_y_yz, g_z_zz_y_zz, g_zz_0_0_0_z_zz_y_xx, g_zz_0_0_0_z_zz_y_xy, g_zz_0_0_0_z_zz_y_xz, g_zz_0_0_0_z_zz_y_yy, g_zz_0_0_0_z_zz_y_yz, g_zz_0_0_0_z_zz_y_zz, g_zzz_zz_y_xx, g_zzz_zz_y_xy, g_zzz_zz_y_xz, g_zzz_zz_y_yy, g_zzz_zz_y_yz, g_zzz_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_zz_y_xx[i] = -6.0 * g_z_zz_y_xx[i] * a_exp + 4.0 * g_zzz_zz_y_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_zz_y_xy[i] = -6.0 * g_z_zz_y_xy[i] * a_exp + 4.0 * g_zzz_zz_y_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_zz_y_xz[i] = -6.0 * g_z_zz_y_xz[i] * a_exp + 4.0 * g_zzz_zz_y_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_zz_y_yy[i] = -6.0 * g_z_zz_y_yy[i] * a_exp + 4.0 * g_zzz_zz_y_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_zz_y_yz[i] = -6.0 * g_z_zz_y_yz[i] * a_exp + 4.0 * g_zzz_zz_y_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_zz_y_zz[i] = -6.0 * g_z_zz_y_zz[i] * a_exp + 4.0 * g_zzz_zz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (1938-1944)

    #pragma omp simd aligned(g_z_zz_z_xx, g_z_zz_z_xy, g_z_zz_z_xz, g_z_zz_z_yy, g_z_zz_z_yz, g_z_zz_z_zz, g_zz_0_0_0_z_zz_z_xx, g_zz_0_0_0_z_zz_z_xy, g_zz_0_0_0_z_zz_z_xz, g_zz_0_0_0_z_zz_z_yy, g_zz_0_0_0_z_zz_z_yz, g_zz_0_0_0_z_zz_z_zz, g_zzz_zz_z_xx, g_zzz_zz_z_xy, g_zzz_zz_z_xz, g_zzz_zz_z_yy, g_zzz_zz_z_yz, g_zzz_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_zz_z_xx[i] = -6.0 * g_z_zz_z_xx[i] * a_exp + 4.0 * g_zzz_zz_z_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_z_zz_z_xy[i] = -6.0 * g_z_zz_z_xy[i] * a_exp + 4.0 * g_zzz_zz_z_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_zz_z_xz[i] = -6.0 * g_z_zz_z_xz[i] * a_exp + 4.0 * g_zzz_zz_z_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_zz_z_yy[i] = -6.0 * g_z_zz_z_yy[i] * a_exp + 4.0 * g_zzz_zz_z_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_z_zz_z_yz[i] = -6.0 * g_z_zz_z_yz[i] * a_exp + 4.0 * g_zzz_zz_z_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_z_zz_z_zz[i] = -6.0 * g_z_zz_z_zz[i] * a_exp + 4.0 * g_zzz_zz_z_zz[i] * a_exp * a_exp;
    }
}

} // t4c_geom namespace


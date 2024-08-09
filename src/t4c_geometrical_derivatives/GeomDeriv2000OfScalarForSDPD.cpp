#include "GeomDeriv2000OfScalarForSDPD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom2000_sdpd_0(CSimdArray<double>& buffer_2000_sdpd,
                     const CSimdArray<double>& buffer_sdpd,
                     const CSimdArray<double>& buffer_ddpd,
                     const double a_exp) -> void
{
    const auto ndims = buffer_2000_sdpd.number_of_columns();

    /// Set up components of auxilary buffer : buffer_sdpd

    auto g_0_xx_x_xx = buffer_sdpd[0];

    auto g_0_xx_x_xy = buffer_sdpd[1];

    auto g_0_xx_x_xz = buffer_sdpd[2];

    auto g_0_xx_x_yy = buffer_sdpd[3];

    auto g_0_xx_x_yz = buffer_sdpd[4];

    auto g_0_xx_x_zz = buffer_sdpd[5];

    auto g_0_xx_y_xx = buffer_sdpd[6];

    auto g_0_xx_y_xy = buffer_sdpd[7];

    auto g_0_xx_y_xz = buffer_sdpd[8];

    auto g_0_xx_y_yy = buffer_sdpd[9];

    auto g_0_xx_y_yz = buffer_sdpd[10];

    auto g_0_xx_y_zz = buffer_sdpd[11];

    auto g_0_xx_z_xx = buffer_sdpd[12];

    auto g_0_xx_z_xy = buffer_sdpd[13];

    auto g_0_xx_z_xz = buffer_sdpd[14];

    auto g_0_xx_z_yy = buffer_sdpd[15];

    auto g_0_xx_z_yz = buffer_sdpd[16];

    auto g_0_xx_z_zz = buffer_sdpd[17];

    auto g_0_xy_x_xx = buffer_sdpd[18];

    auto g_0_xy_x_xy = buffer_sdpd[19];

    auto g_0_xy_x_xz = buffer_sdpd[20];

    auto g_0_xy_x_yy = buffer_sdpd[21];

    auto g_0_xy_x_yz = buffer_sdpd[22];

    auto g_0_xy_x_zz = buffer_sdpd[23];

    auto g_0_xy_y_xx = buffer_sdpd[24];

    auto g_0_xy_y_xy = buffer_sdpd[25];

    auto g_0_xy_y_xz = buffer_sdpd[26];

    auto g_0_xy_y_yy = buffer_sdpd[27];

    auto g_0_xy_y_yz = buffer_sdpd[28];

    auto g_0_xy_y_zz = buffer_sdpd[29];

    auto g_0_xy_z_xx = buffer_sdpd[30];

    auto g_0_xy_z_xy = buffer_sdpd[31];

    auto g_0_xy_z_xz = buffer_sdpd[32];

    auto g_0_xy_z_yy = buffer_sdpd[33];

    auto g_0_xy_z_yz = buffer_sdpd[34];

    auto g_0_xy_z_zz = buffer_sdpd[35];

    auto g_0_xz_x_xx = buffer_sdpd[36];

    auto g_0_xz_x_xy = buffer_sdpd[37];

    auto g_0_xz_x_xz = buffer_sdpd[38];

    auto g_0_xz_x_yy = buffer_sdpd[39];

    auto g_0_xz_x_yz = buffer_sdpd[40];

    auto g_0_xz_x_zz = buffer_sdpd[41];

    auto g_0_xz_y_xx = buffer_sdpd[42];

    auto g_0_xz_y_xy = buffer_sdpd[43];

    auto g_0_xz_y_xz = buffer_sdpd[44];

    auto g_0_xz_y_yy = buffer_sdpd[45];

    auto g_0_xz_y_yz = buffer_sdpd[46];

    auto g_0_xz_y_zz = buffer_sdpd[47];

    auto g_0_xz_z_xx = buffer_sdpd[48];

    auto g_0_xz_z_xy = buffer_sdpd[49];

    auto g_0_xz_z_xz = buffer_sdpd[50];

    auto g_0_xz_z_yy = buffer_sdpd[51];

    auto g_0_xz_z_yz = buffer_sdpd[52];

    auto g_0_xz_z_zz = buffer_sdpd[53];

    auto g_0_yy_x_xx = buffer_sdpd[54];

    auto g_0_yy_x_xy = buffer_sdpd[55];

    auto g_0_yy_x_xz = buffer_sdpd[56];

    auto g_0_yy_x_yy = buffer_sdpd[57];

    auto g_0_yy_x_yz = buffer_sdpd[58];

    auto g_0_yy_x_zz = buffer_sdpd[59];

    auto g_0_yy_y_xx = buffer_sdpd[60];

    auto g_0_yy_y_xy = buffer_sdpd[61];

    auto g_0_yy_y_xz = buffer_sdpd[62];

    auto g_0_yy_y_yy = buffer_sdpd[63];

    auto g_0_yy_y_yz = buffer_sdpd[64];

    auto g_0_yy_y_zz = buffer_sdpd[65];

    auto g_0_yy_z_xx = buffer_sdpd[66];

    auto g_0_yy_z_xy = buffer_sdpd[67];

    auto g_0_yy_z_xz = buffer_sdpd[68];

    auto g_0_yy_z_yy = buffer_sdpd[69];

    auto g_0_yy_z_yz = buffer_sdpd[70];

    auto g_0_yy_z_zz = buffer_sdpd[71];

    auto g_0_yz_x_xx = buffer_sdpd[72];

    auto g_0_yz_x_xy = buffer_sdpd[73];

    auto g_0_yz_x_xz = buffer_sdpd[74];

    auto g_0_yz_x_yy = buffer_sdpd[75];

    auto g_0_yz_x_yz = buffer_sdpd[76];

    auto g_0_yz_x_zz = buffer_sdpd[77];

    auto g_0_yz_y_xx = buffer_sdpd[78];

    auto g_0_yz_y_xy = buffer_sdpd[79];

    auto g_0_yz_y_xz = buffer_sdpd[80];

    auto g_0_yz_y_yy = buffer_sdpd[81];

    auto g_0_yz_y_yz = buffer_sdpd[82];

    auto g_0_yz_y_zz = buffer_sdpd[83];

    auto g_0_yz_z_xx = buffer_sdpd[84];

    auto g_0_yz_z_xy = buffer_sdpd[85];

    auto g_0_yz_z_xz = buffer_sdpd[86];

    auto g_0_yz_z_yy = buffer_sdpd[87];

    auto g_0_yz_z_yz = buffer_sdpd[88];

    auto g_0_yz_z_zz = buffer_sdpd[89];

    auto g_0_zz_x_xx = buffer_sdpd[90];

    auto g_0_zz_x_xy = buffer_sdpd[91];

    auto g_0_zz_x_xz = buffer_sdpd[92];

    auto g_0_zz_x_yy = buffer_sdpd[93];

    auto g_0_zz_x_yz = buffer_sdpd[94];

    auto g_0_zz_x_zz = buffer_sdpd[95];

    auto g_0_zz_y_xx = buffer_sdpd[96];

    auto g_0_zz_y_xy = buffer_sdpd[97];

    auto g_0_zz_y_xz = buffer_sdpd[98];

    auto g_0_zz_y_yy = buffer_sdpd[99];

    auto g_0_zz_y_yz = buffer_sdpd[100];

    auto g_0_zz_y_zz = buffer_sdpd[101];

    auto g_0_zz_z_xx = buffer_sdpd[102];

    auto g_0_zz_z_xy = buffer_sdpd[103];

    auto g_0_zz_z_xz = buffer_sdpd[104];

    auto g_0_zz_z_yy = buffer_sdpd[105];

    auto g_0_zz_z_yz = buffer_sdpd[106];

    auto g_0_zz_z_zz = buffer_sdpd[107];

    /// Set up components of auxilary buffer : buffer_ddpd

    auto g_xx_xx_x_xx = buffer_ddpd[0];

    auto g_xx_xx_x_xy = buffer_ddpd[1];

    auto g_xx_xx_x_xz = buffer_ddpd[2];

    auto g_xx_xx_x_yy = buffer_ddpd[3];

    auto g_xx_xx_x_yz = buffer_ddpd[4];

    auto g_xx_xx_x_zz = buffer_ddpd[5];

    auto g_xx_xx_y_xx = buffer_ddpd[6];

    auto g_xx_xx_y_xy = buffer_ddpd[7];

    auto g_xx_xx_y_xz = buffer_ddpd[8];

    auto g_xx_xx_y_yy = buffer_ddpd[9];

    auto g_xx_xx_y_yz = buffer_ddpd[10];

    auto g_xx_xx_y_zz = buffer_ddpd[11];

    auto g_xx_xx_z_xx = buffer_ddpd[12];

    auto g_xx_xx_z_xy = buffer_ddpd[13];

    auto g_xx_xx_z_xz = buffer_ddpd[14];

    auto g_xx_xx_z_yy = buffer_ddpd[15];

    auto g_xx_xx_z_yz = buffer_ddpd[16];

    auto g_xx_xx_z_zz = buffer_ddpd[17];

    auto g_xx_xy_x_xx = buffer_ddpd[18];

    auto g_xx_xy_x_xy = buffer_ddpd[19];

    auto g_xx_xy_x_xz = buffer_ddpd[20];

    auto g_xx_xy_x_yy = buffer_ddpd[21];

    auto g_xx_xy_x_yz = buffer_ddpd[22];

    auto g_xx_xy_x_zz = buffer_ddpd[23];

    auto g_xx_xy_y_xx = buffer_ddpd[24];

    auto g_xx_xy_y_xy = buffer_ddpd[25];

    auto g_xx_xy_y_xz = buffer_ddpd[26];

    auto g_xx_xy_y_yy = buffer_ddpd[27];

    auto g_xx_xy_y_yz = buffer_ddpd[28];

    auto g_xx_xy_y_zz = buffer_ddpd[29];

    auto g_xx_xy_z_xx = buffer_ddpd[30];

    auto g_xx_xy_z_xy = buffer_ddpd[31];

    auto g_xx_xy_z_xz = buffer_ddpd[32];

    auto g_xx_xy_z_yy = buffer_ddpd[33];

    auto g_xx_xy_z_yz = buffer_ddpd[34];

    auto g_xx_xy_z_zz = buffer_ddpd[35];

    auto g_xx_xz_x_xx = buffer_ddpd[36];

    auto g_xx_xz_x_xy = buffer_ddpd[37];

    auto g_xx_xz_x_xz = buffer_ddpd[38];

    auto g_xx_xz_x_yy = buffer_ddpd[39];

    auto g_xx_xz_x_yz = buffer_ddpd[40];

    auto g_xx_xz_x_zz = buffer_ddpd[41];

    auto g_xx_xz_y_xx = buffer_ddpd[42];

    auto g_xx_xz_y_xy = buffer_ddpd[43];

    auto g_xx_xz_y_xz = buffer_ddpd[44];

    auto g_xx_xz_y_yy = buffer_ddpd[45];

    auto g_xx_xz_y_yz = buffer_ddpd[46];

    auto g_xx_xz_y_zz = buffer_ddpd[47];

    auto g_xx_xz_z_xx = buffer_ddpd[48];

    auto g_xx_xz_z_xy = buffer_ddpd[49];

    auto g_xx_xz_z_xz = buffer_ddpd[50];

    auto g_xx_xz_z_yy = buffer_ddpd[51];

    auto g_xx_xz_z_yz = buffer_ddpd[52];

    auto g_xx_xz_z_zz = buffer_ddpd[53];

    auto g_xx_yy_x_xx = buffer_ddpd[54];

    auto g_xx_yy_x_xy = buffer_ddpd[55];

    auto g_xx_yy_x_xz = buffer_ddpd[56];

    auto g_xx_yy_x_yy = buffer_ddpd[57];

    auto g_xx_yy_x_yz = buffer_ddpd[58];

    auto g_xx_yy_x_zz = buffer_ddpd[59];

    auto g_xx_yy_y_xx = buffer_ddpd[60];

    auto g_xx_yy_y_xy = buffer_ddpd[61];

    auto g_xx_yy_y_xz = buffer_ddpd[62];

    auto g_xx_yy_y_yy = buffer_ddpd[63];

    auto g_xx_yy_y_yz = buffer_ddpd[64];

    auto g_xx_yy_y_zz = buffer_ddpd[65];

    auto g_xx_yy_z_xx = buffer_ddpd[66];

    auto g_xx_yy_z_xy = buffer_ddpd[67];

    auto g_xx_yy_z_xz = buffer_ddpd[68];

    auto g_xx_yy_z_yy = buffer_ddpd[69];

    auto g_xx_yy_z_yz = buffer_ddpd[70];

    auto g_xx_yy_z_zz = buffer_ddpd[71];

    auto g_xx_yz_x_xx = buffer_ddpd[72];

    auto g_xx_yz_x_xy = buffer_ddpd[73];

    auto g_xx_yz_x_xz = buffer_ddpd[74];

    auto g_xx_yz_x_yy = buffer_ddpd[75];

    auto g_xx_yz_x_yz = buffer_ddpd[76];

    auto g_xx_yz_x_zz = buffer_ddpd[77];

    auto g_xx_yz_y_xx = buffer_ddpd[78];

    auto g_xx_yz_y_xy = buffer_ddpd[79];

    auto g_xx_yz_y_xz = buffer_ddpd[80];

    auto g_xx_yz_y_yy = buffer_ddpd[81];

    auto g_xx_yz_y_yz = buffer_ddpd[82];

    auto g_xx_yz_y_zz = buffer_ddpd[83];

    auto g_xx_yz_z_xx = buffer_ddpd[84];

    auto g_xx_yz_z_xy = buffer_ddpd[85];

    auto g_xx_yz_z_xz = buffer_ddpd[86];

    auto g_xx_yz_z_yy = buffer_ddpd[87];

    auto g_xx_yz_z_yz = buffer_ddpd[88];

    auto g_xx_yz_z_zz = buffer_ddpd[89];

    auto g_xx_zz_x_xx = buffer_ddpd[90];

    auto g_xx_zz_x_xy = buffer_ddpd[91];

    auto g_xx_zz_x_xz = buffer_ddpd[92];

    auto g_xx_zz_x_yy = buffer_ddpd[93];

    auto g_xx_zz_x_yz = buffer_ddpd[94];

    auto g_xx_zz_x_zz = buffer_ddpd[95];

    auto g_xx_zz_y_xx = buffer_ddpd[96];

    auto g_xx_zz_y_xy = buffer_ddpd[97];

    auto g_xx_zz_y_xz = buffer_ddpd[98];

    auto g_xx_zz_y_yy = buffer_ddpd[99];

    auto g_xx_zz_y_yz = buffer_ddpd[100];

    auto g_xx_zz_y_zz = buffer_ddpd[101];

    auto g_xx_zz_z_xx = buffer_ddpd[102];

    auto g_xx_zz_z_xy = buffer_ddpd[103];

    auto g_xx_zz_z_xz = buffer_ddpd[104];

    auto g_xx_zz_z_yy = buffer_ddpd[105];

    auto g_xx_zz_z_yz = buffer_ddpd[106];

    auto g_xx_zz_z_zz = buffer_ddpd[107];

    auto g_xy_xx_x_xx = buffer_ddpd[108];

    auto g_xy_xx_x_xy = buffer_ddpd[109];

    auto g_xy_xx_x_xz = buffer_ddpd[110];

    auto g_xy_xx_x_yy = buffer_ddpd[111];

    auto g_xy_xx_x_yz = buffer_ddpd[112];

    auto g_xy_xx_x_zz = buffer_ddpd[113];

    auto g_xy_xx_y_xx = buffer_ddpd[114];

    auto g_xy_xx_y_xy = buffer_ddpd[115];

    auto g_xy_xx_y_xz = buffer_ddpd[116];

    auto g_xy_xx_y_yy = buffer_ddpd[117];

    auto g_xy_xx_y_yz = buffer_ddpd[118];

    auto g_xy_xx_y_zz = buffer_ddpd[119];

    auto g_xy_xx_z_xx = buffer_ddpd[120];

    auto g_xy_xx_z_xy = buffer_ddpd[121];

    auto g_xy_xx_z_xz = buffer_ddpd[122];

    auto g_xy_xx_z_yy = buffer_ddpd[123];

    auto g_xy_xx_z_yz = buffer_ddpd[124];

    auto g_xy_xx_z_zz = buffer_ddpd[125];

    auto g_xy_xy_x_xx = buffer_ddpd[126];

    auto g_xy_xy_x_xy = buffer_ddpd[127];

    auto g_xy_xy_x_xz = buffer_ddpd[128];

    auto g_xy_xy_x_yy = buffer_ddpd[129];

    auto g_xy_xy_x_yz = buffer_ddpd[130];

    auto g_xy_xy_x_zz = buffer_ddpd[131];

    auto g_xy_xy_y_xx = buffer_ddpd[132];

    auto g_xy_xy_y_xy = buffer_ddpd[133];

    auto g_xy_xy_y_xz = buffer_ddpd[134];

    auto g_xy_xy_y_yy = buffer_ddpd[135];

    auto g_xy_xy_y_yz = buffer_ddpd[136];

    auto g_xy_xy_y_zz = buffer_ddpd[137];

    auto g_xy_xy_z_xx = buffer_ddpd[138];

    auto g_xy_xy_z_xy = buffer_ddpd[139];

    auto g_xy_xy_z_xz = buffer_ddpd[140];

    auto g_xy_xy_z_yy = buffer_ddpd[141];

    auto g_xy_xy_z_yz = buffer_ddpd[142];

    auto g_xy_xy_z_zz = buffer_ddpd[143];

    auto g_xy_xz_x_xx = buffer_ddpd[144];

    auto g_xy_xz_x_xy = buffer_ddpd[145];

    auto g_xy_xz_x_xz = buffer_ddpd[146];

    auto g_xy_xz_x_yy = buffer_ddpd[147];

    auto g_xy_xz_x_yz = buffer_ddpd[148];

    auto g_xy_xz_x_zz = buffer_ddpd[149];

    auto g_xy_xz_y_xx = buffer_ddpd[150];

    auto g_xy_xz_y_xy = buffer_ddpd[151];

    auto g_xy_xz_y_xz = buffer_ddpd[152];

    auto g_xy_xz_y_yy = buffer_ddpd[153];

    auto g_xy_xz_y_yz = buffer_ddpd[154];

    auto g_xy_xz_y_zz = buffer_ddpd[155];

    auto g_xy_xz_z_xx = buffer_ddpd[156];

    auto g_xy_xz_z_xy = buffer_ddpd[157];

    auto g_xy_xz_z_xz = buffer_ddpd[158];

    auto g_xy_xz_z_yy = buffer_ddpd[159];

    auto g_xy_xz_z_yz = buffer_ddpd[160];

    auto g_xy_xz_z_zz = buffer_ddpd[161];

    auto g_xy_yy_x_xx = buffer_ddpd[162];

    auto g_xy_yy_x_xy = buffer_ddpd[163];

    auto g_xy_yy_x_xz = buffer_ddpd[164];

    auto g_xy_yy_x_yy = buffer_ddpd[165];

    auto g_xy_yy_x_yz = buffer_ddpd[166];

    auto g_xy_yy_x_zz = buffer_ddpd[167];

    auto g_xy_yy_y_xx = buffer_ddpd[168];

    auto g_xy_yy_y_xy = buffer_ddpd[169];

    auto g_xy_yy_y_xz = buffer_ddpd[170];

    auto g_xy_yy_y_yy = buffer_ddpd[171];

    auto g_xy_yy_y_yz = buffer_ddpd[172];

    auto g_xy_yy_y_zz = buffer_ddpd[173];

    auto g_xy_yy_z_xx = buffer_ddpd[174];

    auto g_xy_yy_z_xy = buffer_ddpd[175];

    auto g_xy_yy_z_xz = buffer_ddpd[176];

    auto g_xy_yy_z_yy = buffer_ddpd[177];

    auto g_xy_yy_z_yz = buffer_ddpd[178];

    auto g_xy_yy_z_zz = buffer_ddpd[179];

    auto g_xy_yz_x_xx = buffer_ddpd[180];

    auto g_xy_yz_x_xy = buffer_ddpd[181];

    auto g_xy_yz_x_xz = buffer_ddpd[182];

    auto g_xy_yz_x_yy = buffer_ddpd[183];

    auto g_xy_yz_x_yz = buffer_ddpd[184];

    auto g_xy_yz_x_zz = buffer_ddpd[185];

    auto g_xy_yz_y_xx = buffer_ddpd[186];

    auto g_xy_yz_y_xy = buffer_ddpd[187];

    auto g_xy_yz_y_xz = buffer_ddpd[188];

    auto g_xy_yz_y_yy = buffer_ddpd[189];

    auto g_xy_yz_y_yz = buffer_ddpd[190];

    auto g_xy_yz_y_zz = buffer_ddpd[191];

    auto g_xy_yz_z_xx = buffer_ddpd[192];

    auto g_xy_yz_z_xy = buffer_ddpd[193];

    auto g_xy_yz_z_xz = buffer_ddpd[194];

    auto g_xy_yz_z_yy = buffer_ddpd[195];

    auto g_xy_yz_z_yz = buffer_ddpd[196];

    auto g_xy_yz_z_zz = buffer_ddpd[197];

    auto g_xy_zz_x_xx = buffer_ddpd[198];

    auto g_xy_zz_x_xy = buffer_ddpd[199];

    auto g_xy_zz_x_xz = buffer_ddpd[200];

    auto g_xy_zz_x_yy = buffer_ddpd[201];

    auto g_xy_zz_x_yz = buffer_ddpd[202];

    auto g_xy_zz_x_zz = buffer_ddpd[203];

    auto g_xy_zz_y_xx = buffer_ddpd[204];

    auto g_xy_zz_y_xy = buffer_ddpd[205];

    auto g_xy_zz_y_xz = buffer_ddpd[206];

    auto g_xy_zz_y_yy = buffer_ddpd[207];

    auto g_xy_zz_y_yz = buffer_ddpd[208];

    auto g_xy_zz_y_zz = buffer_ddpd[209];

    auto g_xy_zz_z_xx = buffer_ddpd[210];

    auto g_xy_zz_z_xy = buffer_ddpd[211];

    auto g_xy_zz_z_xz = buffer_ddpd[212];

    auto g_xy_zz_z_yy = buffer_ddpd[213];

    auto g_xy_zz_z_yz = buffer_ddpd[214];

    auto g_xy_zz_z_zz = buffer_ddpd[215];

    auto g_xz_xx_x_xx = buffer_ddpd[216];

    auto g_xz_xx_x_xy = buffer_ddpd[217];

    auto g_xz_xx_x_xz = buffer_ddpd[218];

    auto g_xz_xx_x_yy = buffer_ddpd[219];

    auto g_xz_xx_x_yz = buffer_ddpd[220];

    auto g_xz_xx_x_zz = buffer_ddpd[221];

    auto g_xz_xx_y_xx = buffer_ddpd[222];

    auto g_xz_xx_y_xy = buffer_ddpd[223];

    auto g_xz_xx_y_xz = buffer_ddpd[224];

    auto g_xz_xx_y_yy = buffer_ddpd[225];

    auto g_xz_xx_y_yz = buffer_ddpd[226];

    auto g_xz_xx_y_zz = buffer_ddpd[227];

    auto g_xz_xx_z_xx = buffer_ddpd[228];

    auto g_xz_xx_z_xy = buffer_ddpd[229];

    auto g_xz_xx_z_xz = buffer_ddpd[230];

    auto g_xz_xx_z_yy = buffer_ddpd[231];

    auto g_xz_xx_z_yz = buffer_ddpd[232];

    auto g_xz_xx_z_zz = buffer_ddpd[233];

    auto g_xz_xy_x_xx = buffer_ddpd[234];

    auto g_xz_xy_x_xy = buffer_ddpd[235];

    auto g_xz_xy_x_xz = buffer_ddpd[236];

    auto g_xz_xy_x_yy = buffer_ddpd[237];

    auto g_xz_xy_x_yz = buffer_ddpd[238];

    auto g_xz_xy_x_zz = buffer_ddpd[239];

    auto g_xz_xy_y_xx = buffer_ddpd[240];

    auto g_xz_xy_y_xy = buffer_ddpd[241];

    auto g_xz_xy_y_xz = buffer_ddpd[242];

    auto g_xz_xy_y_yy = buffer_ddpd[243];

    auto g_xz_xy_y_yz = buffer_ddpd[244];

    auto g_xz_xy_y_zz = buffer_ddpd[245];

    auto g_xz_xy_z_xx = buffer_ddpd[246];

    auto g_xz_xy_z_xy = buffer_ddpd[247];

    auto g_xz_xy_z_xz = buffer_ddpd[248];

    auto g_xz_xy_z_yy = buffer_ddpd[249];

    auto g_xz_xy_z_yz = buffer_ddpd[250];

    auto g_xz_xy_z_zz = buffer_ddpd[251];

    auto g_xz_xz_x_xx = buffer_ddpd[252];

    auto g_xz_xz_x_xy = buffer_ddpd[253];

    auto g_xz_xz_x_xz = buffer_ddpd[254];

    auto g_xz_xz_x_yy = buffer_ddpd[255];

    auto g_xz_xz_x_yz = buffer_ddpd[256];

    auto g_xz_xz_x_zz = buffer_ddpd[257];

    auto g_xz_xz_y_xx = buffer_ddpd[258];

    auto g_xz_xz_y_xy = buffer_ddpd[259];

    auto g_xz_xz_y_xz = buffer_ddpd[260];

    auto g_xz_xz_y_yy = buffer_ddpd[261];

    auto g_xz_xz_y_yz = buffer_ddpd[262];

    auto g_xz_xz_y_zz = buffer_ddpd[263];

    auto g_xz_xz_z_xx = buffer_ddpd[264];

    auto g_xz_xz_z_xy = buffer_ddpd[265];

    auto g_xz_xz_z_xz = buffer_ddpd[266];

    auto g_xz_xz_z_yy = buffer_ddpd[267];

    auto g_xz_xz_z_yz = buffer_ddpd[268];

    auto g_xz_xz_z_zz = buffer_ddpd[269];

    auto g_xz_yy_x_xx = buffer_ddpd[270];

    auto g_xz_yy_x_xy = buffer_ddpd[271];

    auto g_xz_yy_x_xz = buffer_ddpd[272];

    auto g_xz_yy_x_yy = buffer_ddpd[273];

    auto g_xz_yy_x_yz = buffer_ddpd[274];

    auto g_xz_yy_x_zz = buffer_ddpd[275];

    auto g_xz_yy_y_xx = buffer_ddpd[276];

    auto g_xz_yy_y_xy = buffer_ddpd[277];

    auto g_xz_yy_y_xz = buffer_ddpd[278];

    auto g_xz_yy_y_yy = buffer_ddpd[279];

    auto g_xz_yy_y_yz = buffer_ddpd[280];

    auto g_xz_yy_y_zz = buffer_ddpd[281];

    auto g_xz_yy_z_xx = buffer_ddpd[282];

    auto g_xz_yy_z_xy = buffer_ddpd[283];

    auto g_xz_yy_z_xz = buffer_ddpd[284];

    auto g_xz_yy_z_yy = buffer_ddpd[285];

    auto g_xz_yy_z_yz = buffer_ddpd[286];

    auto g_xz_yy_z_zz = buffer_ddpd[287];

    auto g_xz_yz_x_xx = buffer_ddpd[288];

    auto g_xz_yz_x_xy = buffer_ddpd[289];

    auto g_xz_yz_x_xz = buffer_ddpd[290];

    auto g_xz_yz_x_yy = buffer_ddpd[291];

    auto g_xz_yz_x_yz = buffer_ddpd[292];

    auto g_xz_yz_x_zz = buffer_ddpd[293];

    auto g_xz_yz_y_xx = buffer_ddpd[294];

    auto g_xz_yz_y_xy = buffer_ddpd[295];

    auto g_xz_yz_y_xz = buffer_ddpd[296];

    auto g_xz_yz_y_yy = buffer_ddpd[297];

    auto g_xz_yz_y_yz = buffer_ddpd[298];

    auto g_xz_yz_y_zz = buffer_ddpd[299];

    auto g_xz_yz_z_xx = buffer_ddpd[300];

    auto g_xz_yz_z_xy = buffer_ddpd[301];

    auto g_xz_yz_z_xz = buffer_ddpd[302];

    auto g_xz_yz_z_yy = buffer_ddpd[303];

    auto g_xz_yz_z_yz = buffer_ddpd[304];

    auto g_xz_yz_z_zz = buffer_ddpd[305];

    auto g_xz_zz_x_xx = buffer_ddpd[306];

    auto g_xz_zz_x_xy = buffer_ddpd[307];

    auto g_xz_zz_x_xz = buffer_ddpd[308];

    auto g_xz_zz_x_yy = buffer_ddpd[309];

    auto g_xz_zz_x_yz = buffer_ddpd[310];

    auto g_xz_zz_x_zz = buffer_ddpd[311];

    auto g_xz_zz_y_xx = buffer_ddpd[312];

    auto g_xz_zz_y_xy = buffer_ddpd[313];

    auto g_xz_zz_y_xz = buffer_ddpd[314];

    auto g_xz_zz_y_yy = buffer_ddpd[315];

    auto g_xz_zz_y_yz = buffer_ddpd[316];

    auto g_xz_zz_y_zz = buffer_ddpd[317];

    auto g_xz_zz_z_xx = buffer_ddpd[318];

    auto g_xz_zz_z_xy = buffer_ddpd[319];

    auto g_xz_zz_z_xz = buffer_ddpd[320];

    auto g_xz_zz_z_yy = buffer_ddpd[321];

    auto g_xz_zz_z_yz = buffer_ddpd[322];

    auto g_xz_zz_z_zz = buffer_ddpd[323];

    auto g_yy_xx_x_xx = buffer_ddpd[324];

    auto g_yy_xx_x_xy = buffer_ddpd[325];

    auto g_yy_xx_x_xz = buffer_ddpd[326];

    auto g_yy_xx_x_yy = buffer_ddpd[327];

    auto g_yy_xx_x_yz = buffer_ddpd[328];

    auto g_yy_xx_x_zz = buffer_ddpd[329];

    auto g_yy_xx_y_xx = buffer_ddpd[330];

    auto g_yy_xx_y_xy = buffer_ddpd[331];

    auto g_yy_xx_y_xz = buffer_ddpd[332];

    auto g_yy_xx_y_yy = buffer_ddpd[333];

    auto g_yy_xx_y_yz = buffer_ddpd[334];

    auto g_yy_xx_y_zz = buffer_ddpd[335];

    auto g_yy_xx_z_xx = buffer_ddpd[336];

    auto g_yy_xx_z_xy = buffer_ddpd[337];

    auto g_yy_xx_z_xz = buffer_ddpd[338];

    auto g_yy_xx_z_yy = buffer_ddpd[339];

    auto g_yy_xx_z_yz = buffer_ddpd[340];

    auto g_yy_xx_z_zz = buffer_ddpd[341];

    auto g_yy_xy_x_xx = buffer_ddpd[342];

    auto g_yy_xy_x_xy = buffer_ddpd[343];

    auto g_yy_xy_x_xz = buffer_ddpd[344];

    auto g_yy_xy_x_yy = buffer_ddpd[345];

    auto g_yy_xy_x_yz = buffer_ddpd[346];

    auto g_yy_xy_x_zz = buffer_ddpd[347];

    auto g_yy_xy_y_xx = buffer_ddpd[348];

    auto g_yy_xy_y_xy = buffer_ddpd[349];

    auto g_yy_xy_y_xz = buffer_ddpd[350];

    auto g_yy_xy_y_yy = buffer_ddpd[351];

    auto g_yy_xy_y_yz = buffer_ddpd[352];

    auto g_yy_xy_y_zz = buffer_ddpd[353];

    auto g_yy_xy_z_xx = buffer_ddpd[354];

    auto g_yy_xy_z_xy = buffer_ddpd[355];

    auto g_yy_xy_z_xz = buffer_ddpd[356];

    auto g_yy_xy_z_yy = buffer_ddpd[357];

    auto g_yy_xy_z_yz = buffer_ddpd[358];

    auto g_yy_xy_z_zz = buffer_ddpd[359];

    auto g_yy_xz_x_xx = buffer_ddpd[360];

    auto g_yy_xz_x_xy = buffer_ddpd[361];

    auto g_yy_xz_x_xz = buffer_ddpd[362];

    auto g_yy_xz_x_yy = buffer_ddpd[363];

    auto g_yy_xz_x_yz = buffer_ddpd[364];

    auto g_yy_xz_x_zz = buffer_ddpd[365];

    auto g_yy_xz_y_xx = buffer_ddpd[366];

    auto g_yy_xz_y_xy = buffer_ddpd[367];

    auto g_yy_xz_y_xz = buffer_ddpd[368];

    auto g_yy_xz_y_yy = buffer_ddpd[369];

    auto g_yy_xz_y_yz = buffer_ddpd[370];

    auto g_yy_xz_y_zz = buffer_ddpd[371];

    auto g_yy_xz_z_xx = buffer_ddpd[372];

    auto g_yy_xz_z_xy = buffer_ddpd[373];

    auto g_yy_xz_z_xz = buffer_ddpd[374];

    auto g_yy_xz_z_yy = buffer_ddpd[375];

    auto g_yy_xz_z_yz = buffer_ddpd[376];

    auto g_yy_xz_z_zz = buffer_ddpd[377];

    auto g_yy_yy_x_xx = buffer_ddpd[378];

    auto g_yy_yy_x_xy = buffer_ddpd[379];

    auto g_yy_yy_x_xz = buffer_ddpd[380];

    auto g_yy_yy_x_yy = buffer_ddpd[381];

    auto g_yy_yy_x_yz = buffer_ddpd[382];

    auto g_yy_yy_x_zz = buffer_ddpd[383];

    auto g_yy_yy_y_xx = buffer_ddpd[384];

    auto g_yy_yy_y_xy = buffer_ddpd[385];

    auto g_yy_yy_y_xz = buffer_ddpd[386];

    auto g_yy_yy_y_yy = buffer_ddpd[387];

    auto g_yy_yy_y_yz = buffer_ddpd[388];

    auto g_yy_yy_y_zz = buffer_ddpd[389];

    auto g_yy_yy_z_xx = buffer_ddpd[390];

    auto g_yy_yy_z_xy = buffer_ddpd[391];

    auto g_yy_yy_z_xz = buffer_ddpd[392];

    auto g_yy_yy_z_yy = buffer_ddpd[393];

    auto g_yy_yy_z_yz = buffer_ddpd[394];

    auto g_yy_yy_z_zz = buffer_ddpd[395];

    auto g_yy_yz_x_xx = buffer_ddpd[396];

    auto g_yy_yz_x_xy = buffer_ddpd[397];

    auto g_yy_yz_x_xz = buffer_ddpd[398];

    auto g_yy_yz_x_yy = buffer_ddpd[399];

    auto g_yy_yz_x_yz = buffer_ddpd[400];

    auto g_yy_yz_x_zz = buffer_ddpd[401];

    auto g_yy_yz_y_xx = buffer_ddpd[402];

    auto g_yy_yz_y_xy = buffer_ddpd[403];

    auto g_yy_yz_y_xz = buffer_ddpd[404];

    auto g_yy_yz_y_yy = buffer_ddpd[405];

    auto g_yy_yz_y_yz = buffer_ddpd[406];

    auto g_yy_yz_y_zz = buffer_ddpd[407];

    auto g_yy_yz_z_xx = buffer_ddpd[408];

    auto g_yy_yz_z_xy = buffer_ddpd[409];

    auto g_yy_yz_z_xz = buffer_ddpd[410];

    auto g_yy_yz_z_yy = buffer_ddpd[411];

    auto g_yy_yz_z_yz = buffer_ddpd[412];

    auto g_yy_yz_z_zz = buffer_ddpd[413];

    auto g_yy_zz_x_xx = buffer_ddpd[414];

    auto g_yy_zz_x_xy = buffer_ddpd[415];

    auto g_yy_zz_x_xz = buffer_ddpd[416];

    auto g_yy_zz_x_yy = buffer_ddpd[417];

    auto g_yy_zz_x_yz = buffer_ddpd[418];

    auto g_yy_zz_x_zz = buffer_ddpd[419];

    auto g_yy_zz_y_xx = buffer_ddpd[420];

    auto g_yy_zz_y_xy = buffer_ddpd[421];

    auto g_yy_zz_y_xz = buffer_ddpd[422];

    auto g_yy_zz_y_yy = buffer_ddpd[423];

    auto g_yy_zz_y_yz = buffer_ddpd[424];

    auto g_yy_zz_y_zz = buffer_ddpd[425];

    auto g_yy_zz_z_xx = buffer_ddpd[426];

    auto g_yy_zz_z_xy = buffer_ddpd[427];

    auto g_yy_zz_z_xz = buffer_ddpd[428];

    auto g_yy_zz_z_yy = buffer_ddpd[429];

    auto g_yy_zz_z_yz = buffer_ddpd[430];

    auto g_yy_zz_z_zz = buffer_ddpd[431];

    auto g_yz_xx_x_xx = buffer_ddpd[432];

    auto g_yz_xx_x_xy = buffer_ddpd[433];

    auto g_yz_xx_x_xz = buffer_ddpd[434];

    auto g_yz_xx_x_yy = buffer_ddpd[435];

    auto g_yz_xx_x_yz = buffer_ddpd[436];

    auto g_yz_xx_x_zz = buffer_ddpd[437];

    auto g_yz_xx_y_xx = buffer_ddpd[438];

    auto g_yz_xx_y_xy = buffer_ddpd[439];

    auto g_yz_xx_y_xz = buffer_ddpd[440];

    auto g_yz_xx_y_yy = buffer_ddpd[441];

    auto g_yz_xx_y_yz = buffer_ddpd[442];

    auto g_yz_xx_y_zz = buffer_ddpd[443];

    auto g_yz_xx_z_xx = buffer_ddpd[444];

    auto g_yz_xx_z_xy = buffer_ddpd[445];

    auto g_yz_xx_z_xz = buffer_ddpd[446];

    auto g_yz_xx_z_yy = buffer_ddpd[447];

    auto g_yz_xx_z_yz = buffer_ddpd[448];

    auto g_yz_xx_z_zz = buffer_ddpd[449];

    auto g_yz_xy_x_xx = buffer_ddpd[450];

    auto g_yz_xy_x_xy = buffer_ddpd[451];

    auto g_yz_xy_x_xz = buffer_ddpd[452];

    auto g_yz_xy_x_yy = buffer_ddpd[453];

    auto g_yz_xy_x_yz = buffer_ddpd[454];

    auto g_yz_xy_x_zz = buffer_ddpd[455];

    auto g_yz_xy_y_xx = buffer_ddpd[456];

    auto g_yz_xy_y_xy = buffer_ddpd[457];

    auto g_yz_xy_y_xz = buffer_ddpd[458];

    auto g_yz_xy_y_yy = buffer_ddpd[459];

    auto g_yz_xy_y_yz = buffer_ddpd[460];

    auto g_yz_xy_y_zz = buffer_ddpd[461];

    auto g_yz_xy_z_xx = buffer_ddpd[462];

    auto g_yz_xy_z_xy = buffer_ddpd[463];

    auto g_yz_xy_z_xz = buffer_ddpd[464];

    auto g_yz_xy_z_yy = buffer_ddpd[465];

    auto g_yz_xy_z_yz = buffer_ddpd[466];

    auto g_yz_xy_z_zz = buffer_ddpd[467];

    auto g_yz_xz_x_xx = buffer_ddpd[468];

    auto g_yz_xz_x_xy = buffer_ddpd[469];

    auto g_yz_xz_x_xz = buffer_ddpd[470];

    auto g_yz_xz_x_yy = buffer_ddpd[471];

    auto g_yz_xz_x_yz = buffer_ddpd[472];

    auto g_yz_xz_x_zz = buffer_ddpd[473];

    auto g_yz_xz_y_xx = buffer_ddpd[474];

    auto g_yz_xz_y_xy = buffer_ddpd[475];

    auto g_yz_xz_y_xz = buffer_ddpd[476];

    auto g_yz_xz_y_yy = buffer_ddpd[477];

    auto g_yz_xz_y_yz = buffer_ddpd[478];

    auto g_yz_xz_y_zz = buffer_ddpd[479];

    auto g_yz_xz_z_xx = buffer_ddpd[480];

    auto g_yz_xz_z_xy = buffer_ddpd[481];

    auto g_yz_xz_z_xz = buffer_ddpd[482];

    auto g_yz_xz_z_yy = buffer_ddpd[483];

    auto g_yz_xz_z_yz = buffer_ddpd[484];

    auto g_yz_xz_z_zz = buffer_ddpd[485];

    auto g_yz_yy_x_xx = buffer_ddpd[486];

    auto g_yz_yy_x_xy = buffer_ddpd[487];

    auto g_yz_yy_x_xz = buffer_ddpd[488];

    auto g_yz_yy_x_yy = buffer_ddpd[489];

    auto g_yz_yy_x_yz = buffer_ddpd[490];

    auto g_yz_yy_x_zz = buffer_ddpd[491];

    auto g_yz_yy_y_xx = buffer_ddpd[492];

    auto g_yz_yy_y_xy = buffer_ddpd[493];

    auto g_yz_yy_y_xz = buffer_ddpd[494];

    auto g_yz_yy_y_yy = buffer_ddpd[495];

    auto g_yz_yy_y_yz = buffer_ddpd[496];

    auto g_yz_yy_y_zz = buffer_ddpd[497];

    auto g_yz_yy_z_xx = buffer_ddpd[498];

    auto g_yz_yy_z_xy = buffer_ddpd[499];

    auto g_yz_yy_z_xz = buffer_ddpd[500];

    auto g_yz_yy_z_yy = buffer_ddpd[501];

    auto g_yz_yy_z_yz = buffer_ddpd[502];

    auto g_yz_yy_z_zz = buffer_ddpd[503];

    auto g_yz_yz_x_xx = buffer_ddpd[504];

    auto g_yz_yz_x_xy = buffer_ddpd[505];

    auto g_yz_yz_x_xz = buffer_ddpd[506];

    auto g_yz_yz_x_yy = buffer_ddpd[507];

    auto g_yz_yz_x_yz = buffer_ddpd[508];

    auto g_yz_yz_x_zz = buffer_ddpd[509];

    auto g_yz_yz_y_xx = buffer_ddpd[510];

    auto g_yz_yz_y_xy = buffer_ddpd[511];

    auto g_yz_yz_y_xz = buffer_ddpd[512];

    auto g_yz_yz_y_yy = buffer_ddpd[513];

    auto g_yz_yz_y_yz = buffer_ddpd[514];

    auto g_yz_yz_y_zz = buffer_ddpd[515];

    auto g_yz_yz_z_xx = buffer_ddpd[516];

    auto g_yz_yz_z_xy = buffer_ddpd[517];

    auto g_yz_yz_z_xz = buffer_ddpd[518];

    auto g_yz_yz_z_yy = buffer_ddpd[519];

    auto g_yz_yz_z_yz = buffer_ddpd[520];

    auto g_yz_yz_z_zz = buffer_ddpd[521];

    auto g_yz_zz_x_xx = buffer_ddpd[522];

    auto g_yz_zz_x_xy = buffer_ddpd[523];

    auto g_yz_zz_x_xz = buffer_ddpd[524];

    auto g_yz_zz_x_yy = buffer_ddpd[525];

    auto g_yz_zz_x_yz = buffer_ddpd[526];

    auto g_yz_zz_x_zz = buffer_ddpd[527];

    auto g_yz_zz_y_xx = buffer_ddpd[528];

    auto g_yz_zz_y_xy = buffer_ddpd[529];

    auto g_yz_zz_y_xz = buffer_ddpd[530];

    auto g_yz_zz_y_yy = buffer_ddpd[531];

    auto g_yz_zz_y_yz = buffer_ddpd[532];

    auto g_yz_zz_y_zz = buffer_ddpd[533];

    auto g_yz_zz_z_xx = buffer_ddpd[534];

    auto g_yz_zz_z_xy = buffer_ddpd[535];

    auto g_yz_zz_z_xz = buffer_ddpd[536];

    auto g_yz_zz_z_yy = buffer_ddpd[537];

    auto g_yz_zz_z_yz = buffer_ddpd[538];

    auto g_yz_zz_z_zz = buffer_ddpd[539];

    auto g_zz_xx_x_xx = buffer_ddpd[540];

    auto g_zz_xx_x_xy = buffer_ddpd[541];

    auto g_zz_xx_x_xz = buffer_ddpd[542];

    auto g_zz_xx_x_yy = buffer_ddpd[543];

    auto g_zz_xx_x_yz = buffer_ddpd[544];

    auto g_zz_xx_x_zz = buffer_ddpd[545];

    auto g_zz_xx_y_xx = buffer_ddpd[546];

    auto g_zz_xx_y_xy = buffer_ddpd[547];

    auto g_zz_xx_y_xz = buffer_ddpd[548];

    auto g_zz_xx_y_yy = buffer_ddpd[549];

    auto g_zz_xx_y_yz = buffer_ddpd[550];

    auto g_zz_xx_y_zz = buffer_ddpd[551];

    auto g_zz_xx_z_xx = buffer_ddpd[552];

    auto g_zz_xx_z_xy = buffer_ddpd[553];

    auto g_zz_xx_z_xz = buffer_ddpd[554];

    auto g_zz_xx_z_yy = buffer_ddpd[555];

    auto g_zz_xx_z_yz = buffer_ddpd[556];

    auto g_zz_xx_z_zz = buffer_ddpd[557];

    auto g_zz_xy_x_xx = buffer_ddpd[558];

    auto g_zz_xy_x_xy = buffer_ddpd[559];

    auto g_zz_xy_x_xz = buffer_ddpd[560];

    auto g_zz_xy_x_yy = buffer_ddpd[561];

    auto g_zz_xy_x_yz = buffer_ddpd[562];

    auto g_zz_xy_x_zz = buffer_ddpd[563];

    auto g_zz_xy_y_xx = buffer_ddpd[564];

    auto g_zz_xy_y_xy = buffer_ddpd[565];

    auto g_zz_xy_y_xz = buffer_ddpd[566];

    auto g_zz_xy_y_yy = buffer_ddpd[567];

    auto g_zz_xy_y_yz = buffer_ddpd[568];

    auto g_zz_xy_y_zz = buffer_ddpd[569];

    auto g_zz_xy_z_xx = buffer_ddpd[570];

    auto g_zz_xy_z_xy = buffer_ddpd[571];

    auto g_zz_xy_z_xz = buffer_ddpd[572];

    auto g_zz_xy_z_yy = buffer_ddpd[573];

    auto g_zz_xy_z_yz = buffer_ddpd[574];

    auto g_zz_xy_z_zz = buffer_ddpd[575];

    auto g_zz_xz_x_xx = buffer_ddpd[576];

    auto g_zz_xz_x_xy = buffer_ddpd[577];

    auto g_zz_xz_x_xz = buffer_ddpd[578];

    auto g_zz_xz_x_yy = buffer_ddpd[579];

    auto g_zz_xz_x_yz = buffer_ddpd[580];

    auto g_zz_xz_x_zz = buffer_ddpd[581];

    auto g_zz_xz_y_xx = buffer_ddpd[582];

    auto g_zz_xz_y_xy = buffer_ddpd[583];

    auto g_zz_xz_y_xz = buffer_ddpd[584];

    auto g_zz_xz_y_yy = buffer_ddpd[585];

    auto g_zz_xz_y_yz = buffer_ddpd[586];

    auto g_zz_xz_y_zz = buffer_ddpd[587];

    auto g_zz_xz_z_xx = buffer_ddpd[588];

    auto g_zz_xz_z_xy = buffer_ddpd[589];

    auto g_zz_xz_z_xz = buffer_ddpd[590];

    auto g_zz_xz_z_yy = buffer_ddpd[591];

    auto g_zz_xz_z_yz = buffer_ddpd[592];

    auto g_zz_xz_z_zz = buffer_ddpd[593];

    auto g_zz_yy_x_xx = buffer_ddpd[594];

    auto g_zz_yy_x_xy = buffer_ddpd[595];

    auto g_zz_yy_x_xz = buffer_ddpd[596];

    auto g_zz_yy_x_yy = buffer_ddpd[597];

    auto g_zz_yy_x_yz = buffer_ddpd[598];

    auto g_zz_yy_x_zz = buffer_ddpd[599];

    auto g_zz_yy_y_xx = buffer_ddpd[600];

    auto g_zz_yy_y_xy = buffer_ddpd[601];

    auto g_zz_yy_y_xz = buffer_ddpd[602];

    auto g_zz_yy_y_yy = buffer_ddpd[603];

    auto g_zz_yy_y_yz = buffer_ddpd[604];

    auto g_zz_yy_y_zz = buffer_ddpd[605];

    auto g_zz_yy_z_xx = buffer_ddpd[606];

    auto g_zz_yy_z_xy = buffer_ddpd[607];

    auto g_zz_yy_z_xz = buffer_ddpd[608];

    auto g_zz_yy_z_yy = buffer_ddpd[609];

    auto g_zz_yy_z_yz = buffer_ddpd[610];

    auto g_zz_yy_z_zz = buffer_ddpd[611];

    auto g_zz_yz_x_xx = buffer_ddpd[612];

    auto g_zz_yz_x_xy = buffer_ddpd[613];

    auto g_zz_yz_x_xz = buffer_ddpd[614];

    auto g_zz_yz_x_yy = buffer_ddpd[615];

    auto g_zz_yz_x_yz = buffer_ddpd[616];

    auto g_zz_yz_x_zz = buffer_ddpd[617];

    auto g_zz_yz_y_xx = buffer_ddpd[618];

    auto g_zz_yz_y_xy = buffer_ddpd[619];

    auto g_zz_yz_y_xz = buffer_ddpd[620];

    auto g_zz_yz_y_yy = buffer_ddpd[621];

    auto g_zz_yz_y_yz = buffer_ddpd[622];

    auto g_zz_yz_y_zz = buffer_ddpd[623];

    auto g_zz_yz_z_xx = buffer_ddpd[624];

    auto g_zz_yz_z_xy = buffer_ddpd[625];

    auto g_zz_yz_z_xz = buffer_ddpd[626];

    auto g_zz_yz_z_yy = buffer_ddpd[627];

    auto g_zz_yz_z_yz = buffer_ddpd[628];

    auto g_zz_yz_z_zz = buffer_ddpd[629];

    auto g_zz_zz_x_xx = buffer_ddpd[630];

    auto g_zz_zz_x_xy = buffer_ddpd[631];

    auto g_zz_zz_x_xz = buffer_ddpd[632];

    auto g_zz_zz_x_yy = buffer_ddpd[633];

    auto g_zz_zz_x_yz = buffer_ddpd[634];

    auto g_zz_zz_x_zz = buffer_ddpd[635];

    auto g_zz_zz_y_xx = buffer_ddpd[636];

    auto g_zz_zz_y_xy = buffer_ddpd[637];

    auto g_zz_zz_y_xz = buffer_ddpd[638];

    auto g_zz_zz_y_yy = buffer_ddpd[639];

    auto g_zz_zz_y_yz = buffer_ddpd[640];

    auto g_zz_zz_y_zz = buffer_ddpd[641];

    auto g_zz_zz_z_xx = buffer_ddpd[642];

    auto g_zz_zz_z_xy = buffer_ddpd[643];

    auto g_zz_zz_z_xz = buffer_ddpd[644];

    auto g_zz_zz_z_yy = buffer_ddpd[645];

    auto g_zz_zz_z_yz = buffer_ddpd[646];

    auto g_zz_zz_z_zz = buffer_ddpd[647];

    /// Set up components of integrals buffer : buffer_2000_sdpd

    auto g_xx_0_0_0_0_xx_x_xx = buffer_2000_sdpd[0];

    auto g_xx_0_0_0_0_xx_x_xy = buffer_2000_sdpd[1];

    auto g_xx_0_0_0_0_xx_x_xz = buffer_2000_sdpd[2];

    auto g_xx_0_0_0_0_xx_x_yy = buffer_2000_sdpd[3];

    auto g_xx_0_0_0_0_xx_x_yz = buffer_2000_sdpd[4];

    auto g_xx_0_0_0_0_xx_x_zz = buffer_2000_sdpd[5];

    auto g_xx_0_0_0_0_xx_y_xx = buffer_2000_sdpd[6];

    auto g_xx_0_0_0_0_xx_y_xy = buffer_2000_sdpd[7];

    auto g_xx_0_0_0_0_xx_y_xz = buffer_2000_sdpd[8];

    auto g_xx_0_0_0_0_xx_y_yy = buffer_2000_sdpd[9];

    auto g_xx_0_0_0_0_xx_y_yz = buffer_2000_sdpd[10];

    auto g_xx_0_0_0_0_xx_y_zz = buffer_2000_sdpd[11];

    auto g_xx_0_0_0_0_xx_z_xx = buffer_2000_sdpd[12];

    auto g_xx_0_0_0_0_xx_z_xy = buffer_2000_sdpd[13];

    auto g_xx_0_0_0_0_xx_z_xz = buffer_2000_sdpd[14];

    auto g_xx_0_0_0_0_xx_z_yy = buffer_2000_sdpd[15];

    auto g_xx_0_0_0_0_xx_z_yz = buffer_2000_sdpd[16];

    auto g_xx_0_0_0_0_xx_z_zz = buffer_2000_sdpd[17];

    auto g_xx_0_0_0_0_xy_x_xx = buffer_2000_sdpd[18];

    auto g_xx_0_0_0_0_xy_x_xy = buffer_2000_sdpd[19];

    auto g_xx_0_0_0_0_xy_x_xz = buffer_2000_sdpd[20];

    auto g_xx_0_0_0_0_xy_x_yy = buffer_2000_sdpd[21];

    auto g_xx_0_0_0_0_xy_x_yz = buffer_2000_sdpd[22];

    auto g_xx_0_0_0_0_xy_x_zz = buffer_2000_sdpd[23];

    auto g_xx_0_0_0_0_xy_y_xx = buffer_2000_sdpd[24];

    auto g_xx_0_0_0_0_xy_y_xy = buffer_2000_sdpd[25];

    auto g_xx_0_0_0_0_xy_y_xz = buffer_2000_sdpd[26];

    auto g_xx_0_0_0_0_xy_y_yy = buffer_2000_sdpd[27];

    auto g_xx_0_0_0_0_xy_y_yz = buffer_2000_sdpd[28];

    auto g_xx_0_0_0_0_xy_y_zz = buffer_2000_sdpd[29];

    auto g_xx_0_0_0_0_xy_z_xx = buffer_2000_sdpd[30];

    auto g_xx_0_0_0_0_xy_z_xy = buffer_2000_sdpd[31];

    auto g_xx_0_0_0_0_xy_z_xz = buffer_2000_sdpd[32];

    auto g_xx_0_0_0_0_xy_z_yy = buffer_2000_sdpd[33];

    auto g_xx_0_0_0_0_xy_z_yz = buffer_2000_sdpd[34];

    auto g_xx_0_0_0_0_xy_z_zz = buffer_2000_sdpd[35];

    auto g_xx_0_0_0_0_xz_x_xx = buffer_2000_sdpd[36];

    auto g_xx_0_0_0_0_xz_x_xy = buffer_2000_sdpd[37];

    auto g_xx_0_0_0_0_xz_x_xz = buffer_2000_sdpd[38];

    auto g_xx_0_0_0_0_xz_x_yy = buffer_2000_sdpd[39];

    auto g_xx_0_0_0_0_xz_x_yz = buffer_2000_sdpd[40];

    auto g_xx_0_0_0_0_xz_x_zz = buffer_2000_sdpd[41];

    auto g_xx_0_0_0_0_xz_y_xx = buffer_2000_sdpd[42];

    auto g_xx_0_0_0_0_xz_y_xy = buffer_2000_sdpd[43];

    auto g_xx_0_0_0_0_xz_y_xz = buffer_2000_sdpd[44];

    auto g_xx_0_0_0_0_xz_y_yy = buffer_2000_sdpd[45];

    auto g_xx_0_0_0_0_xz_y_yz = buffer_2000_sdpd[46];

    auto g_xx_0_0_0_0_xz_y_zz = buffer_2000_sdpd[47];

    auto g_xx_0_0_0_0_xz_z_xx = buffer_2000_sdpd[48];

    auto g_xx_0_0_0_0_xz_z_xy = buffer_2000_sdpd[49];

    auto g_xx_0_0_0_0_xz_z_xz = buffer_2000_sdpd[50];

    auto g_xx_0_0_0_0_xz_z_yy = buffer_2000_sdpd[51];

    auto g_xx_0_0_0_0_xz_z_yz = buffer_2000_sdpd[52];

    auto g_xx_0_0_0_0_xz_z_zz = buffer_2000_sdpd[53];

    auto g_xx_0_0_0_0_yy_x_xx = buffer_2000_sdpd[54];

    auto g_xx_0_0_0_0_yy_x_xy = buffer_2000_sdpd[55];

    auto g_xx_0_0_0_0_yy_x_xz = buffer_2000_sdpd[56];

    auto g_xx_0_0_0_0_yy_x_yy = buffer_2000_sdpd[57];

    auto g_xx_0_0_0_0_yy_x_yz = buffer_2000_sdpd[58];

    auto g_xx_0_0_0_0_yy_x_zz = buffer_2000_sdpd[59];

    auto g_xx_0_0_0_0_yy_y_xx = buffer_2000_sdpd[60];

    auto g_xx_0_0_0_0_yy_y_xy = buffer_2000_sdpd[61];

    auto g_xx_0_0_0_0_yy_y_xz = buffer_2000_sdpd[62];

    auto g_xx_0_0_0_0_yy_y_yy = buffer_2000_sdpd[63];

    auto g_xx_0_0_0_0_yy_y_yz = buffer_2000_sdpd[64];

    auto g_xx_0_0_0_0_yy_y_zz = buffer_2000_sdpd[65];

    auto g_xx_0_0_0_0_yy_z_xx = buffer_2000_sdpd[66];

    auto g_xx_0_0_0_0_yy_z_xy = buffer_2000_sdpd[67];

    auto g_xx_0_0_0_0_yy_z_xz = buffer_2000_sdpd[68];

    auto g_xx_0_0_0_0_yy_z_yy = buffer_2000_sdpd[69];

    auto g_xx_0_0_0_0_yy_z_yz = buffer_2000_sdpd[70];

    auto g_xx_0_0_0_0_yy_z_zz = buffer_2000_sdpd[71];

    auto g_xx_0_0_0_0_yz_x_xx = buffer_2000_sdpd[72];

    auto g_xx_0_0_0_0_yz_x_xy = buffer_2000_sdpd[73];

    auto g_xx_0_0_0_0_yz_x_xz = buffer_2000_sdpd[74];

    auto g_xx_0_0_0_0_yz_x_yy = buffer_2000_sdpd[75];

    auto g_xx_0_0_0_0_yz_x_yz = buffer_2000_sdpd[76];

    auto g_xx_0_0_0_0_yz_x_zz = buffer_2000_sdpd[77];

    auto g_xx_0_0_0_0_yz_y_xx = buffer_2000_sdpd[78];

    auto g_xx_0_0_0_0_yz_y_xy = buffer_2000_sdpd[79];

    auto g_xx_0_0_0_0_yz_y_xz = buffer_2000_sdpd[80];

    auto g_xx_0_0_0_0_yz_y_yy = buffer_2000_sdpd[81];

    auto g_xx_0_0_0_0_yz_y_yz = buffer_2000_sdpd[82];

    auto g_xx_0_0_0_0_yz_y_zz = buffer_2000_sdpd[83];

    auto g_xx_0_0_0_0_yz_z_xx = buffer_2000_sdpd[84];

    auto g_xx_0_0_0_0_yz_z_xy = buffer_2000_sdpd[85];

    auto g_xx_0_0_0_0_yz_z_xz = buffer_2000_sdpd[86];

    auto g_xx_0_0_0_0_yz_z_yy = buffer_2000_sdpd[87];

    auto g_xx_0_0_0_0_yz_z_yz = buffer_2000_sdpd[88];

    auto g_xx_0_0_0_0_yz_z_zz = buffer_2000_sdpd[89];

    auto g_xx_0_0_0_0_zz_x_xx = buffer_2000_sdpd[90];

    auto g_xx_0_0_0_0_zz_x_xy = buffer_2000_sdpd[91];

    auto g_xx_0_0_0_0_zz_x_xz = buffer_2000_sdpd[92];

    auto g_xx_0_0_0_0_zz_x_yy = buffer_2000_sdpd[93];

    auto g_xx_0_0_0_0_zz_x_yz = buffer_2000_sdpd[94];

    auto g_xx_0_0_0_0_zz_x_zz = buffer_2000_sdpd[95];

    auto g_xx_0_0_0_0_zz_y_xx = buffer_2000_sdpd[96];

    auto g_xx_0_0_0_0_zz_y_xy = buffer_2000_sdpd[97];

    auto g_xx_0_0_0_0_zz_y_xz = buffer_2000_sdpd[98];

    auto g_xx_0_0_0_0_zz_y_yy = buffer_2000_sdpd[99];

    auto g_xx_0_0_0_0_zz_y_yz = buffer_2000_sdpd[100];

    auto g_xx_0_0_0_0_zz_y_zz = buffer_2000_sdpd[101];

    auto g_xx_0_0_0_0_zz_z_xx = buffer_2000_sdpd[102];

    auto g_xx_0_0_0_0_zz_z_xy = buffer_2000_sdpd[103];

    auto g_xx_0_0_0_0_zz_z_xz = buffer_2000_sdpd[104];

    auto g_xx_0_0_0_0_zz_z_yy = buffer_2000_sdpd[105];

    auto g_xx_0_0_0_0_zz_z_yz = buffer_2000_sdpd[106];

    auto g_xx_0_0_0_0_zz_z_zz = buffer_2000_sdpd[107];

    auto g_xy_0_0_0_0_xx_x_xx = buffer_2000_sdpd[108];

    auto g_xy_0_0_0_0_xx_x_xy = buffer_2000_sdpd[109];

    auto g_xy_0_0_0_0_xx_x_xz = buffer_2000_sdpd[110];

    auto g_xy_0_0_0_0_xx_x_yy = buffer_2000_sdpd[111];

    auto g_xy_0_0_0_0_xx_x_yz = buffer_2000_sdpd[112];

    auto g_xy_0_0_0_0_xx_x_zz = buffer_2000_sdpd[113];

    auto g_xy_0_0_0_0_xx_y_xx = buffer_2000_sdpd[114];

    auto g_xy_0_0_0_0_xx_y_xy = buffer_2000_sdpd[115];

    auto g_xy_0_0_0_0_xx_y_xz = buffer_2000_sdpd[116];

    auto g_xy_0_0_0_0_xx_y_yy = buffer_2000_sdpd[117];

    auto g_xy_0_0_0_0_xx_y_yz = buffer_2000_sdpd[118];

    auto g_xy_0_0_0_0_xx_y_zz = buffer_2000_sdpd[119];

    auto g_xy_0_0_0_0_xx_z_xx = buffer_2000_sdpd[120];

    auto g_xy_0_0_0_0_xx_z_xy = buffer_2000_sdpd[121];

    auto g_xy_0_0_0_0_xx_z_xz = buffer_2000_sdpd[122];

    auto g_xy_0_0_0_0_xx_z_yy = buffer_2000_sdpd[123];

    auto g_xy_0_0_0_0_xx_z_yz = buffer_2000_sdpd[124];

    auto g_xy_0_0_0_0_xx_z_zz = buffer_2000_sdpd[125];

    auto g_xy_0_0_0_0_xy_x_xx = buffer_2000_sdpd[126];

    auto g_xy_0_0_0_0_xy_x_xy = buffer_2000_sdpd[127];

    auto g_xy_0_0_0_0_xy_x_xz = buffer_2000_sdpd[128];

    auto g_xy_0_0_0_0_xy_x_yy = buffer_2000_sdpd[129];

    auto g_xy_0_0_0_0_xy_x_yz = buffer_2000_sdpd[130];

    auto g_xy_0_0_0_0_xy_x_zz = buffer_2000_sdpd[131];

    auto g_xy_0_0_0_0_xy_y_xx = buffer_2000_sdpd[132];

    auto g_xy_0_0_0_0_xy_y_xy = buffer_2000_sdpd[133];

    auto g_xy_0_0_0_0_xy_y_xz = buffer_2000_sdpd[134];

    auto g_xy_0_0_0_0_xy_y_yy = buffer_2000_sdpd[135];

    auto g_xy_0_0_0_0_xy_y_yz = buffer_2000_sdpd[136];

    auto g_xy_0_0_0_0_xy_y_zz = buffer_2000_sdpd[137];

    auto g_xy_0_0_0_0_xy_z_xx = buffer_2000_sdpd[138];

    auto g_xy_0_0_0_0_xy_z_xy = buffer_2000_sdpd[139];

    auto g_xy_0_0_0_0_xy_z_xz = buffer_2000_sdpd[140];

    auto g_xy_0_0_0_0_xy_z_yy = buffer_2000_sdpd[141];

    auto g_xy_0_0_0_0_xy_z_yz = buffer_2000_sdpd[142];

    auto g_xy_0_0_0_0_xy_z_zz = buffer_2000_sdpd[143];

    auto g_xy_0_0_0_0_xz_x_xx = buffer_2000_sdpd[144];

    auto g_xy_0_0_0_0_xz_x_xy = buffer_2000_sdpd[145];

    auto g_xy_0_0_0_0_xz_x_xz = buffer_2000_sdpd[146];

    auto g_xy_0_0_0_0_xz_x_yy = buffer_2000_sdpd[147];

    auto g_xy_0_0_0_0_xz_x_yz = buffer_2000_sdpd[148];

    auto g_xy_0_0_0_0_xz_x_zz = buffer_2000_sdpd[149];

    auto g_xy_0_0_0_0_xz_y_xx = buffer_2000_sdpd[150];

    auto g_xy_0_0_0_0_xz_y_xy = buffer_2000_sdpd[151];

    auto g_xy_0_0_0_0_xz_y_xz = buffer_2000_sdpd[152];

    auto g_xy_0_0_0_0_xz_y_yy = buffer_2000_sdpd[153];

    auto g_xy_0_0_0_0_xz_y_yz = buffer_2000_sdpd[154];

    auto g_xy_0_0_0_0_xz_y_zz = buffer_2000_sdpd[155];

    auto g_xy_0_0_0_0_xz_z_xx = buffer_2000_sdpd[156];

    auto g_xy_0_0_0_0_xz_z_xy = buffer_2000_sdpd[157];

    auto g_xy_0_0_0_0_xz_z_xz = buffer_2000_sdpd[158];

    auto g_xy_0_0_0_0_xz_z_yy = buffer_2000_sdpd[159];

    auto g_xy_0_0_0_0_xz_z_yz = buffer_2000_sdpd[160];

    auto g_xy_0_0_0_0_xz_z_zz = buffer_2000_sdpd[161];

    auto g_xy_0_0_0_0_yy_x_xx = buffer_2000_sdpd[162];

    auto g_xy_0_0_0_0_yy_x_xy = buffer_2000_sdpd[163];

    auto g_xy_0_0_0_0_yy_x_xz = buffer_2000_sdpd[164];

    auto g_xy_0_0_0_0_yy_x_yy = buffer_2000_sdpd[165];

    auto g_xy_0_0_0_0_yy_x_yz = buffer_2000_sdpd[166];

    auto g_xy_0_0_0_0_yy_x_zz = buffer_2000_sdpd[167];

    auto g_xy_0_0_0_0_yy_y_xx = buffer_2000_sdpd[168];

    auto g_xy_0_0_0_0_yy_y_xy = buffer_2000_sdpd[169];

    auto g_xy_0_0_0_0_yy_y_xz = buffer_2000_sdpd[170];

    auto g_xy_0_0_0_0_yy_y_yy = buffer_2000_sdpd[171];

    auto g_xy_0_0_0_0_yy_y_yz = buffer_2000_sdpd[172];

    auto g_xy_0_0_0_0_yy_y_zz = buffer_2000_sdpd[173];

    auto g_xy_0_0_0_0_yy_z_xx = buffer_2000_sdpd[174];

    auto g_xy_0_0_0_0_yy_z_xy = buffer_2000_sdpd[175];

    auto g_xy_0_0_0_0_yy_z_xz = buffer_2000_sdpd[176];

    auto g_xy_0_0_0_0_yy_z_yy = buffer_2000_sdpd[177];

    auto g_xy_0_0_0_0_yy_z_yz = buffer_2000_sdpd[178];

    auto g_xy_0_0_0_0_yy_z_zz = buffer_2000_sdpd[179];

    auto g_xy_0_0_0_0_yz_x_xx = buffer_2000_sdpd[180];

    auto g_xy_0_0_0_0_yz_x_xy = buffer_2000_sdpd[181];

    auto g_xy_0_0_0_0_yz_x_xz = buffer_2000_sdpd[182];

    auto g_xy_0_0_0_0_yz_x_yy = buffer_2000_sdpd[183];

    auto g_xy_0_0_0_0_yz_x_yz = buffer_2000_sdpd[184];

    auto g_xy_0_0_0_0_yz_x_zz = buffer_2000_sdpd[185];

    auto g_xy_0_0_0_0_yz_y_xx = buffer_2000_sdpd[186];

    auto g_xy_0_0_0_0_yz_y_xy = buffer_2000_sdpd[187];

    auto g_xy_0_0_0_0_yz_y_xz = buffer_2000_sdpd[188];

    auto g_xy_0_0_0_0_yz_y_yy = buffer_2000_sdpd[189];

    auto g_xy_0_0_0_0_yz_y_yz = buffer_2000_sdpd[190];

    auto g_xy_0_0_0_0_yz_y_zz = buffer_2000_sdpd[191];

    auto g_xy_0_0_0_0_yz_z_xx = buffer_2000_sdpd[192];

    auto g_xy_0_0_0_0_yz_z_xy = buffer_2000_sdpd[193];

    auto g_xy_0_0_0_0_yz_z_xz = buffer_2000_sdpd[194];

    auto g_xy_0_0_0_0_yz_z_yy = buffer_2000_sdpd[195];

    auto g_xy_0_0_0_0_yz_z_yz = buffer_2000_sdpd[196];

    auto g_xy_0_0_0_0_yz_z_zz = buffer_2000_sdpd[197];

    auto g_xy_0_0_0_0_zz_x_xx = buffer_2000_sdpd[198];

    auto g_xy_0_0_0_0_zz_x_xy = buffer_2000_sdpd[199];

    auto g_xy_0_0_0_0_zz_x_xz = buffer_2000_sdpd[200];

    auto g_xy_0_0_0_0_zz_x_yy = buffer_2000_sdpd[201];

    auto g_xy_0_0_0_0_zz_x_yz = buffer_2000_sdpd[202];

    auto g_xy_0_0_0_0_zz_x_zz = buffer_2000_sdpd[203];

    auto g_xy_0_0_0_0_zz_y_xx = buffer_2000_sdpd[204];

    auto g_xy_0_0_0_0_zz_y_xy = buffer_2000_sdpd[205];

    auto g_xy_0_0_0_0_zz_y_xz = buffer_2000_sdpd[206];

    auto g_xy_0_0_0_0_zz_y_yy = buffer_2000_sdpd[207];

    auto g_xy_0_0_0_0_zz_y_yz = buffer_2000_sdpd[208];

    auto g_xy_0_0_0_0_zz_y_zz = buffer_2000_sdpd[209];

    auto g_xy_0_0_0_0_zz_z_xx = buffer_2000_sdpd[210];

    auto g_xy_0_0_0_0_zz_z_xy = buffer_2000_sdpd[211];

    auto g_xy_0_0_0_0_zz_z_xz = buffer_2000_sdpd[212];

    auto g_xy_0_0_0_0_zz_z_yy = buffer_2000_sdpd[213];

    auto g_xy_0_0_0_0_zz_z_yz = buffer_2000_sdpd[214];

    auto g_xy_0_0_0_0_zz_z_zz = buffer_2000_sdpd[215];

    auto g_xz_0_0_0_0_xx_x_xx = buffer_2000_sdpd[216];

    auto g_xz_0_0_0_0_xx_x_xy = buffer_2000_sdpd[217];

    auto g_xz_0_0_0_0_xx_x_xz = buffer_2000_sdpd[218];

    auto g_xz_0_0_0_0_xx_x_yy = buffer_2000_sdpd[219];

    auto g_xz_0_0_0_0_xx_x_yz = buffer_2000_sdpd[220];

    auto g_xz_0_0_0_0_xx_x_zz = buffer_2000_sdpd[221];

    auto g_xz_0_0_0_0_xx_y_xx = buffer_2000_sdpd[222];

    auto g_xz_0_0_0_0_xx_y_xy = buffer_2000_sdpd[223];

    auto g_xz_0_0_0_0_xx_y_xz = buffer_2000_sdpd[224];

    auto g_xz_0_0_0_0_xx_y_yy = buffer_2000_sdpd[225];

    auto g_xz_0_0_0_0_xx_y_yz = buffer_2000_sdpd[226];

    auto g_xz_0_0_0_0_xx_y_zz = buffer_2000_sdpd[227];

    auto g_xz_0_0_0_0_xx_z_xx = buffer_2000_sdpd[228];

    auto g_xz_0_0_0_0_xx_z_xy = buffer_2000_sdpd[229];

    auto g_xz_0_0_0_0_xx_z_xz = buffer_2000_sdpd[230];

    auto g_xz_0_0_0_0_xx_z_yy = buffer_2000_sdpd[231];

    auto g_xz_0_0_0_0_xx_z_yz = buffer_2000_sdpd[232];

    auto g_xz_0_0_0_0_xx_z_zz = buffer_2000_sdpd[233];

    auto g_xz_0_0_0_0_xy_x_xx = buffer_2000_sdpd[234];

    auto g_xz_0_0_0_0_xy_x_xy = buffer_2000_sdpd[235];

    auto g_xz_0_0_0_0_xy_x_xz = buffer_2000_sdpd[236];

    auto g_xz_0_0_0_0_xy_x_yy = buffer_2000_sdpd[237];

    auto g_xz_0_0_0_0_xy_x_yz = buffer_2000_sdpd[238];

    auto g_xz_0_0_0_0_xy_x_zz = buffer_2000_sdpd[239];

    auto g_xz_0_0_0_0_xy_y_xx = buffer_2000_sdpd[240];

    auto g_xz_0_0_0_0_xy_y_xy = buffer_2000_sdpd[241];

    auto g_xz_0_0_0_0_xy_y_xz = buffer_2000_sdpd[242];

    auto g_xz_0_0_0_0_xy_y_yy = buffer_2000_sdpd[243];

    auto g_xz_0_0_0_0_xy_y_yz = buffer_2000_sdpd[244];

    auto g_xz_0_0_0_0_xy_y_zz = buffer_2000_sdpd[245];

    auto g_xz_0_0_0_0_xy_z_xx = buffer_2000_sdpd[246];

    auto g_xz_0_0_0_0_xy_z_xy = buffer_2000_sdpd[247];

    auto g_xz_0_0_0_0_xy_z_xz = buffer_2000_sdpd[248];

    auto g_xz_0_0_0_0_xy_z_yy = buffer_2000_sdpd[249];

    auto g_xz_0_0_0_0_xy_z_yz = buffer_2000_sdpd[250];

    auto g_xz_0_0_0_0_xy_z_zz = buffer_2000_sdpd[251];

    auto g_xz_0_0_0_0_xz_x_xx = buffer_2000_sdpd[252];

    auto g_xz_0_0_0_0_xz_x_xy = buffer_2000_sdpd[253];

    auto g_xz_0_0_0_0_xz_x_xz = buffer_2000_sdpd[254];

    auto g_xz_0_0_0_0_xz_x_yy = buffer_2000_sdpd[255];

    auto g_xz_0_0_0_0_xz_x_yz = buffer_2000_sdpd[256];

    auto g_xz_0_0_0_0_xz_x_zz = buffer_2000_sdpd[257];

    auto g_xz_0_0_0_0_xz_y_xx = buffer_2000_sdpd[258];

    auto g_xz_0_0_0_0_xz_y_xy = buffer_2000_sdpd[259];

    auto g_xz_0_0_0_0_xz_y_xz = buffer_2000_sdpd[260];

    auto g_xz_0_0_0_0_xz_y_yy = buffer_2000_sdpd[261];

    auto g_xz_0_0_0_0_xz_y_yz = buffer_2000_sdpd[262];

    auto g_xz_0_0_0_0_xz_y_zz = buffer_2000_sdpd[263];

    auto g_xz_0_0_0_0_xz_z_xx = buffer_2000_sdpd[264];

    auto g_xz_0_0_0_0_xz_z_xy = buffer_2000_sdpd[265];

    auto g_xz_0_0_0_0_xz_z_xz = buffer_2000_sdpd[266];

    auto g_xz_0_0_0_0_xz_z_yy = buffer_2000_sdpd[267];

    auto g_xz_0_0_0_0_xz_z_yz = buffer_2000_sdpd[268];

    auto g_xz_0_0_0_0_xz_z_zz = buffer_2000_sdpd[269];

    auto g_xz_0_0_0_0_yy_x_xx = buffer_2000_sdpd[270];

    auto g_xz_0_0_0_0_yy_x_xy = buffer_2000_sdpd[271];

    auto g_xz_0_0_0_0_yy_x_xz = buffer_2000_sdpd[272];

    auto g_xz_0_0_0_0_yy_x_yy = buffer_2000_sdpd[273];

    auto g_xz_0_0_0_0_yy_x_yz = buffer_2000_sdpd[274];

    auto g_xz_0_0_0_0_yy_x_zz = buffer_2000_sdpd[275];

    auto g_xz_0_0_0_0_yy_y_xx = buffer_2000_sdpd[276];

    auto g_xz_0_0_0_0_yy_y_xy = buffer_2000_sdpd[277];

    auto g_xz_0_0_0_0_yy_y_xz = buffer_2000_sdpd[278];

    auto g_xz_0_0_0_0_yy_y_yy = buffer_2000_sdpd[279];

    auto g_xz_0_0_0_0_yy_y_yz = buffer_2000_sdpd[280];

    auto g_xz_0_0_0_0_yy_y_zz = buffer_2000_sdpd[281];

    auto g_xz_0_0_0_0_yy_z_xx = buffer_2000_sdpd[282];

    auto g_xz_0_0_0_0_yy_z_xy = buffer_2000_sdpd[283];

    auto g_xz_0_0_0_0_yy_z_xz = buffer_2000_sdpd[284];

    auto g_xz_0_0_0_0_yy_z_yy = buffer_2000_sdpd[285];

    auto g_xz_0_0_0_0_yy_z_yz = buffer_2000_sdpd[286];

    auto g_xz_0_0_0_0_yy_z_zz = buffer_2000_sdpd[287];

    auto g_xz_0_0_0_0_yz_x_xx = buffer_2000_sdpd[288];

    auto g_xz_0_0_0_0_yz_x_xy = buffer_2000_sdpd[289];

    auto g_xz_0_0_0_0_yz_x_xz = buffer_2000_sdpd[290];

    auto g_xz_0_0_0_0_yz_x_yy = buffer_2000_sdpd[291];

    auto g_xz_0_0_0_0_yz_x_yz = buffer_2000_sdpd[292];

    auto g_xz_0_0_0_0_yz_x_zz = buffer_2000_sdpd[293];

    auto g_xz_0_0_0_0_yz_y_xx = buffer_2000_sdpd[294];

    auto g_xz_0_0_0_0_yz_y_xy = buffer_2000_sdpd[295];

    auto g_xz_0_0_0_0_yz_y_xz = buffer_2000_sdpd[296];

    auto g_xz_0_0_0_0_yz_y_yy = buffer_2000_sdpd[297];

    auto g_xz_0_0_0_0_yz_y_yz = buffer_2000_sdpd[298];

    auto g_xz_0_0_0_0_yz_y_zz = buffer_2000_sdpd[299];

    auto g_xz_0_0_0_0_yz_z_xx = buffer_2000_sdpd[300];

    auto g_xz_0_0_0_0_yz_z_xy = buffer_2000_sdpd[301];

    auto g_xz_0_0_0_0_yz_z_xz = buffer_2000_sdpd[302];

    auto g_xz_0_0_0_0_yz_z_yy = buffer_2000_sdpd[303];

    auto g_xz_0_0_0_0_yz_z_yz = buffer_2000_sdpd[304];

    auto g_xz_0_0_0_0_yz_z_zz = buffer_2000_sdpd[305];

    auto g_xz_0_0_0_0_zz_x_xx = buffer_2000_sdpd[306];

    auto g_xz_0_0_0_0_zz_x_xy = buffer_2000_sdpd[307];

    auto g_xz_0_0_0_0_zz_x_xz = buffer_2000_sdpd[308];

    auto g_xz_0_0_0_0_zz_x_yy = buffer_2000_sdpd[309];

    auto g_xz_0_0_0_0_zz_x_yz = buffer_2000_sdpd[310];

    auto g_xz_0_0_0_0_zz_x_zz = buffer_2000_sdpd[311];

    auto g_xz_0_0_0_0_zz_y_xx = buffer_2000_sdpd[312];

    auto g_xz_0_0_0_0_zz_y_xy = buffer_2000_sdpd[313];

    auto g_xz_0_0_0_0_zz_y_xz = buffer_2000_sdpd[314];

    auto g_xz_0_0_0_0_zz_y_yy = buffer_2000_sdpd[315];

    auto g_xz_0_0_0_0_zz_y_yz = buffer_2000_sdpd[316];

    auto g_xz_0_0_0_0_zz_y_zz = buffer_2000_sdpd[317];

    auto g_xz_0_0_0_0_zz_z_xx = buffer_2000_sdpd[318];

    auto g_xz_0_0_0_0_zz_z_xy = buffer_2000_sdpd[319];

    auto g_xz_0_0_0_0_zz_z_xz = buffer_2000_sdpd[320];

    auto g_xz_0_0_0_0_zz_z_yy = buffer_2000_sdpd[321];

    auto g_xz_0_0_0_0_zz_z_yz = buffer_2000_sdpd[322];

    auto g_xz_0_0_0_0_zz_z_zz = buffer_2000_sdpd[323];

    auto g_yy_0_0_0_0_xx_x_xx = buffer_2000_sdpd[324];

    auto g_yy_0_0_0_0_xx_x_xy = buffer_2000_sdpd[325];

    auto g_yy_0_0_0_0_xx_x_xz = buffer_2000_sdpd[326];

    auto g_yy_0_0_0_0_xx_x_yy = buffer_2000_sdpd[327];

    auto g_yy_0_0_0_0_xx_x_yz = buffer_2000_sdpd[328];

    auto g_yy_0_0_0_0_xx_x_zz = buffer_2000_sdpd[329];

    auto g_yy_0_0_0_0_xx_y_xx = buffer_2000_sdpd[330];

    auto g_yy_0_0_0_0_xx_y_xy = buffer_2000_sdpd[331];

    auto g_yy_0_0_0_0_xx_y_xz = buffer_2000_sdpd[332];

    auto g_yy_0_0_0_0_xx_y_yy = buffer_2000_sdpd[333];

    auto g_yy_0_0_0_0_xx_y_yz = buffer_2000_sdpd[334];

    auto g_yy_0_0_0_0_xx_y_zz = buffer_2000_sdpd[335];

    auto g_yy_0_0_0_0_xx_z_xx = buffer_2000_sdpd[336];

    auto g_yy_0_0_0_0_xx_z_xy = buffer_2000_sdpd[337];

    auto g_yy_0_0_0_0_xx_z_xz = buffer_2000_sdpd[338];

    auto g_yy_0_0_0_0_xx_z_yy = buffer_2000_sdpd[339];

    auto g_yy_0_0_0_0_xx_z_yz = buffer_2000_sdpd[340];

    auto g_yy_0_0_0_0_xx_z_zz = buffer_2000_sdpd[341];

    auto g_yy_0_0_0_0_xy_x_xx = buffer_2000_sdpd[342];

    auto g_yy_0_0_0_0_xy_x_xy = buffer_2000_sdpd[343];

    auto g_yy_0_0_0_0_xy_x_xz = buffer_2000_sdpd[344];

    auto g_yy_0_0_0_0_xy_x_yy = buffer_2000_sdpd[345];

    auto g_yy_0_0_0_0_xy_x_yz = buffer_2000_sdpd[346];

    auto g_yy_0_0_0_0_xy_x_zz = buffer_2000_sdpd[347];

    auto g_yy_0_0_0_0_xy_y_xx = buffer_2000_sdpd[348];

    auto g_yy_0_0_0_0_xy_y_xy = buffer_2000_sdpd[349];

    auto g_yy_0_0_0_0_xy_y_xz = buffer_2000_sdpd[350];

    auto g_yy_0_0_0_0_xy_y_yy = buffer_2000_sdpd[351];

    auto g_yy_0_0_0_0_xy_y_yz = buffer_2000_sdpd[352];

    auto g_yy_0_0_0_0_xy_y_zz = buffer_2000_sdpd[353];

    auto g_yy_0_0_0_0_xy_z_xx = buffer_2000_sdpd[354];

    auto g_yy_0_0_0_0_xy_z_xy = buffer_2000_sdpd[355];

    auto g_yy_0_0_0_0_xy_z_xz = buffer_2000_sdpd[356];

    auto g_yy_0_0_0_0_xy_z_yy = buffer_2000_sdpd[357];

    auto g_yy_0_0_0_0_xy_z_yz = buffer_2000_sdpd[358];

    auto g_yy_0_0_0_0_xy_z_zz = buffer_2000_sdpd[359];

    auto g_yy_0_0_0_0_xz_x_xx = buffer_2000_sdpd[360];

    auto g_yy_0_0_0_0_xz_x_xy = buffer_2000_sdpd[361];

    auto g_yy_0_0_0_0_xz_x_xz = buffer_2000_sdpd[362];

    auto g_yy_0_0_0_0_xz_x_yy = buffer_2000_sdpd[363];

    auto g_yy_0_0_0_0_xz_x_yz = buffer_2000_sdpd[364];

    auto g_yy_0_0_0_0_xz_x_zz = buffer_2000_sdpd[365];

    auto g_yy_0_0_0_0_xz_y_xx = buffer_2000_sdpd[366];

    auto g_yy_0_0_0_0_xz_y_xy = buffer_2000_sdpd[367];

    auto g_yy_0_0_0_0_xz_y_xz = buffer_2000_sdpd[368];

    auto g_yy_0_0_0_0_xz_y_yy = buffer_2000_sdpd[369];

    auto g_yy_0_0_0_0_xz_y_yz = buffer_2000_sdpd[370];

    auto g_yy_0_0_0_0_xz_y_zz = buffer_2000_sdpd[371];

    auto g_yy_0_0_0_0_xz_z_xx = buffer_2000_sdpd[372];

    auto g_yy_0_0_0_0_xz_z_xy = buffer_2000_sdpd[373];

    auto g_yy_0_0_0_0_xz_z_xz = buffer_2000_sdpd[374];

    auto g_yy_0_0_0_0_xz_z_yy = buffer_2000_sdpd[375];

    auto g_yy_0_0_0_0_xz_z_yz = buffer_2000_sdpd[376];

    auto g_yy_0_0_0_0_xz_z_zz = buffer_2000_sdpd[377];

    auto g_yy_0_0_0_0_yy_x_xx = buffer_2000_sdpd[378];

    auto g_yy_0_0_0_0_yy_x_xy = buffer_2000_sdpd[379];

    auto g_yy_0_0_0_0_yy_x_xz = buffer_2000_sdpd[380];

    auto g_yy_0_0_0_0_yy_x_yy = buffer_2000_sdpd[381];

    auto g_yy_0_0_0_0_yy_x_yz = buffer_2000_sdpd[382];

    auto g_yy_0_0_0_0_yy_x_zz = buffer_2000_sdpd[383];

    auto g_yy_0_0_0_0_yy_y_xx = buffer_2000_sdpd[384];

    auto g_yy_0_0_0_0_yy_y_xy = buffer_2000_sdpd[385];

    auto g_yy_0_0_0_0_yy_y_xz = buffer_2000_sdpd[386];

    auto g_yy_0_0_0_0_yy_y_yy = buffer_2000_sdpd[387];

    auto g_yy_0_0_0_0_yy_y_yz = buffer_2000_sdpd[388];

    auto g_yy_0_0_0_0_yy_y_zz = buffer_2000_sdpd[389];

    auto g_yy_0_0_0_0_yy_z_xx = buffer_2000_sdpd[390];

    auto g_yy_0_0_0_0_yy_z_xy = buffer_2000_sdpd[391];

    auto g_yy_0_0_0_0_yy_z_xz = buffer_2000_sdpd[392];

    auto g_yy_0_0_0_0_yy_z_yy = buffer_2000_sdpd[393];

    auto g_yy_0_0_0_0_yy_z_yz = buffer_2000_sdpd[394];

    auto g_yy_0_0_0_0_yy_z_zz = buffer_2000_sdpd[395];

    auto g_yy_0_0_0_0_yz_x_xx = buffer_2000_sdpd[396];

    auto g_yy_0_0_0_0_yz_x_xy = buffer_2000_sdpd[397];

    auto g_yy_0_0_0_0_yz_x_xz = buffer_2000_sdpd[398];

    auto g_yy_0_0_0_0_yz_x_yy = buffer_2000_sdpd[399];

    auto g_yy_0_0_0_0_yz_x_yz = buffer_2000_sdpd[400];

    auto g_yy_0_0_0_0_yz_x_zz = buffer_2000_sdpd[401];

    auto g_yy_0_0_0_0_yz_y_xx = buffer_2000_sdpd[402];

    auto g_yy_0_0_0_0_yz_y_xy = buffer_2000_sdpd[403];

    auto g_yy_0_0_0_0_yz_y_xz = buffer_2000_sdpd[404];

    auto g_yy_0_0_0_0_yz_y_yy = buffer_2000_sdpd[405];

    auto g_yy_0_0_0_0_yz_y_yz = buffer_2000_sdpd[406];

    auto g_yy_0_0_0_0_yz_y_zz = buffer_2000_sdpd[407];

    auto g_yy_0_0_0_0_yz_z_xx = buffer_2000_sdpd[408];

    auto g_yy_0_0_0_0_yz_z_xy = buffer_2000_sdpd[409];

    auto g_yy_0_0_0_0_yz_z_xz = buffer_2000_sdpd[410];

    auto g_yy_0_0_0_0_yz_z_yy = buffer_2000_sdpd[411];

    auto g_yy_0_0_0_0_yz_z_yz = buffer_2000_sdpd[412];

    auto g_yy_0_0_0_0_yz_z_zz = buffer_2000_sdpd[413];

    auto g_yy_0_0_0_0_zz_x_xx = buffer_2000_sdpd[414];

    auto g_yy_0_0_0_0_zz_x_xy = buffer_2000_sdpd[415];

    auto g_yy_0_0_0_0_zz_x_xz = buffer_2000_sdpd[416];

    auto g_yy_0_0_0_0_zz_x_yy = buffer_2000_sdpd[417];

    auto g_yy_0_0_0_0_zz_x_yz = buffer_2000_sdpd[418];

    auto g_yy_0_0_0_0_zz_x_zz = buffer_2000_sdpd[419];

    auto g_yy_0_0_0_0_zz_y_xx = buffer_2000_sdpd[420];

    auto g_yy_0_0_0_0_zz_y_xy = buffer_2000_sdpd[421];

    auto g_yy_0_0_0_0_zz_y_xz = buffer_2000_sdpd[422];

    auto g_yy_0_0_0_0_zz_y_yy = buffer_2000_sdpd[423];

    auto g_yy_0_0_0_0_zz_y_yz = buffer_2000_sdpd[424];

    auto g_yy_0_0_0_0_zz_y_zz = buffer_2000_sdpd[425];

    auto g_yy_0_0_0_0_zz_z_xx = buffer_2000_sdpd[426];

    auto g_yy_0_0_0_0_zz_z_xy = buffer_2000_sdpd[427];

    auto g_yy_0_0_0_0_zz_z_xz = buffer_2000_sdpd[428];

    auto g_yy_0_0_0_0_zz_z_yy = buffer_2000_sdpd[429];

    auto g_yy_0_0_0_0_zz_z_yz = buffer_2000_sdpd[430];

    auto g_yy_0_0_0_0_zz_z_zz = buffer_2000_sdpd[431];

    auto g_yz_0_0_0_0_xx_x_xx = buffer_2000_sdpd[432];

    auto g_yz_0_0_0_0_xx_x_xy = buffer_2000_sdpd[433];

    auto g_yz_0_0_0_0_xx_x_xz = buffer_2000_sdpd[434];

    auto g_yz_0_0_0_0_xx_x_yy = buffer_2000_sdpd[435];

    auto g_yz_0_0_0_0_xx_x_yz = buffer_2000_sdpd[436];

    auto g_yz_0_0_0_0_xx_x_zz = buffer_2000_sdpd[437];

    auto g_yz_0_0_0_0_xx_y_xx = buffer_2000_sdpd[438];

    auto g_yz_0_0_0_0_xx_y_xy = buffer_2000_sdpd[439];

    auto g_yz_0_0_0_0_xx_y_xz = buffer_2000_sdpd[440];

    auto g_yz_0_0_0_0_xx_y_yy = buffer_2000_sdpd[441];

    auto g_yz_0_0_0_0_xx_y_yz = buffer_2000_sdpd[442];

    auto g_yz_0_0_0_0_xx_y_zz = buffer_2000_sdpd[443];

    auto g_yz_0_0_0_0_xx_z_xx = buffer_2000_sdpd[444];

    auto g_yz_0_0_0_0_xx_z_xy = buffer_2000_sdpd[445];

    auto g_yz_0_0_0_0_xx_z_xz = buffer_2000_sdpd[446];

    auto g_yz_0_0_0_0_xx_z_yy = buffer_2000_sdpd[447];

    auto g_yz_0_0_0_0_xx_z_yz = buffer_2000_sdpd[448];

    auto g_yz_0_0_0_0_xx_z_zz = buffer_2000_sdpd[449];

    auto g_yz_0_0_0_0_xy_x_xx = buffer_2000_sdpd[450];

    auto g_yz_0_0_0_0_xy_x_xy = buffer_2000_sdpd[451];

    auto g_yz_0_0_0_0_xy_x_xz = buffer_2000_sdpd[452];

    auto g_yz_0_0_0_0_xy_x_yy = buffer_2000_sdpd[453];

    auto g_yz_0_0_0_0_xy_x_yz = buffer_2000_sdpd[454];

    auto g_yz_0_0_0_0_xy_x_zz = buffer_2000_sdpd[455];

    auto g_yz_0_0_0_0_xy_y_xx = buffer_2000_sdpd[456];

    auto g_yz_0_0_0_0_xy_y_xy = buffer_2000_sdpd[457];

    auto g_yz_0_0_0_0_xy_y_xz = buffer_2000_sdpd[458];

    auto g_yz_0_0_0_0_xy_y_yy = buffer_2000_sdpd[459];

    auto g_yz_0_0_0_0_xy_y_yz = buffer_2000_sdpd[460];

    auto g_yz_0_0_0_0_xy_y_zz = buffer_2000_sdpd[461];

    auto g_yz_0_0_0_0_xy_z_xx = buffer_2000_sdpd[462];

    auto g_yz_0_0_0_0_xy_z_xy = buffer_2000_sdpd[463];

    auto g_yz_0_0_0_0_xy_z_xz = buffer_2000_sdpd[464];

    auto g_yz_0_0_0_0_xy_z_yy = buffer_2000_sdpd[465];

    auto g_yz_0_0_0_0_xy_z_yz = buffer_2000_sdpd[466];

    auto g_yz_0_0_0_0_xy_z_zz = buffer_2000_sdpd[467];

    auto g_yz_0_0_0_0_xz_x_xx = buffer_2000_sdpd[468];

    auto g_yz_0_0_0_0_xz_x_xy = buffer_2000_sdpd[469];

    auto g_yz_0_0_0_0_xz_x_xz = buffer_2000_sdpd[470];

    auto g_yz_0_0_0_0_xz_x_yy = buffer_2000_sdpd[471];

    auto g_yz_0_0_0_0_xz_x_yz = buffer_2000_sdpd[472];

    auto g_yz_0_0_0_0_xz_x_zz = buffer_2000_sdpd[473];

    auto g_yz_0_0_0_0_xz_y_xx = buffer_2000_sdpd[474];

    auto g_yz_0_0_0_0_xz_y_xy = buffer_2000_sdpd[475];

    auto g_yz_0_0_0_0_xz_y_xz = buffer_2000_sdpd[476];

    auto g_yz_0_0_0_0_xz_y_yy = buffer_2000_sdpd[477];

    auto g_yz_0_0_0_0_xz_y_yz = buffer_2000_sdpd[478];

    auto g_yz_0_0_0_0_xz_y_zz = buffer_2000_sdpd[479];

    auto g_yz_0_0_0_0_xz_z_xx = buffer_2000_sdpd[480];

    auto g_yz_0_0_0_0_xz_z_xy = buffer_2000_sdpd[481];

    auto g_yz_0_0_0_0_xz_z_xz = buffer_2000_sdpd[482];

    auto g_yz_0_0_0_0_xz_z_yy = buffer_2000_sdpd[483];

    auto g_yz_0_0_0_0_xz_z_yz = buffer_2000_sdpd[484];

    auto g_yz_0_0_0_0_xz_z_zz = buffer_2000_sdpd[485];

    auto g_yz_0_0_0_0_yy_x_xx = buffer_2000_sdpd[486];

    auto g_yz_0_0_0_0_yy_x_xy = buffer_2000_sdpd[487];

    auto g_yz_0_0_0_0_yy_x_xz = buffer_2000_sdpd[488];

    auto g_yz_0_0_0_0_yy_x_yy = buffer_2000_sdpd[489];

    auto g_yz_0_0_0_0_yy_x_yz = buffer_2000_sdpd[490];

    auto g_yz_0_0_0_0_yy_x_zz = buffer_2000_sdpd[491];

    auto g_yz_0_0_0_0_yy_y_xx = buffer_2000_sdpd[492];

    auto g_yz_0_0_0_0_yy_y_xy = buffer_2000_sdpd[493];

    auto g_yz_0_0_0_0_yy_y_xz = buffer_2000_sdpd[494];

    auto g_yz_0_0_0_0_yy_y_yy = buffer_2000_sdpd[495];

    auto g_yz_0_0_0_0_yy_y_yz = buffer_2000_sdpd[496];

    auto g_yz_0_0_0_0_yy_y_zz = buffer_2000_sdpd[497];

    auto g_yz_0_0_0_0_yy_z_xx = buffer_2000_sdpd[498];

    auto g_yz_0_0_0_0_yy_z_xy = buffer_2000_sdpd[499];

    auto g_yz_0_0_0_0_yy_z_xz = buffer_2000_sdpd[500];

    auto g_yz_0_0_0_0_yy_z_yy = buffer_2000_sdpd[501];

    auto g_yz_0_0_0_0_yy_z_yz = buffer_2000_sdpd[502];

    auto g_yz_0_0_0_0_yy_z_zz = buffer_2000_sdpd[503];

    auto g_yz_0_0_0_0_yz_x_xx = buffer_2000_sdpd[504];

    auto g_yz_0_0_0_0_yz_x_xy = buffer_2000_sdpd[505];

    auto g_yz_0_0_0_0_yz_x_xz = buffer_2000_sdpd[506];

    auto g_yz_0_0_0_0_yz_x_yy = buffer_2000_sdpd[507];

    auto g_yz_0_0_0_0_yz_x_yz = buffer_2000_sdpd[508];

    auto g_yz_0_0_0_0_yz_x_zz = buffer_2000_sdpd[509];

    auto g_yz_0_0_0_0_yz_y_xx = buffer_2000_sdpd[510];

    auto g_yz_0_0_0_0_yz_y_xy = buffer_2000_sdpd[511];

    auto g_yz_0_0_0_0_yz_y_xz = buffer_2000_sdpd[512];

    auto g_yz_0_0_0_0_yz_y_yy = buffer_2000_sdpd[513];

    auto g_yz_0_0_0_0_yz_y_yz = buffer_2000_sdpd[514];

    auto g_yz_0_0_0_0_yz_y_zz = buffer_2000_sdpd[515];

    auto g_yz_0_0_0_0_yz_z_xx = buffer_2000_sdpd[516];

    auto g_yz_0_0_0_0_yz_z_xy = buffer_2000_sdpd[517];

    auto g_yz_0_0_0_0_yz_z_xz = buffer_2000_sdpd[518];

    auto g_yz_0_0_0_0_yz_z_yy = buffer_2000_sdpd[519];

    auto g_yz_0_0_0_0_yz_z_yz = buffer_2000_sdpd[520];

    auto g_yz_0_0_0_0_yz_z_zz = buffer_2000_sdpd[521];

    auto g_yz_0_0_0_0_zz_x_xx = buffer_2000_sdpd[522];

    auto g_yz_0_0_0_0_zz_x_xy = buffer_2000_sdpd[523];

    auto g_yz_0_0_0_0_zz_x_xz = buffer_2000_sdpd[524];

    auto g_yz_0_0_0_0_zz_x_yy = buffer_2000_sdpd[525];

    auto g_yz_0_0_0_0_zz_x_yz = buffer_2000_sdpd[526];

    auto g_yz_0_0_0_0_zz_x_zz = buffer_2000_sdpd[527];

    auto g_yz_0_0_0_0_zz_y_xx = buffer_2000_sdpd[528];

    auto g_yz_0_0_0_0_zz_y_xy = buffer_2000_sdpd[529];

    auto g_yz_0_0_0_0_zz_y_xz = buffer_2000_sdpd[530];

    auto g_yz_0_0_0_0_zz_y_yy = buffer_2000_sdpd[531];

    auto g_yz_0_0_0_0_zz_y_yz = buffer_2000_sdpd[532];

    auto g_yz_0_0_0_0_zz_y_zz = buffer_2000_sdpd[533];

    auto g_yz_0_0_0_0_zz_z_xx = buffer_2000_sdpd[534];

    auto g_yz_0_0_0_0_zz_z_xy = buffer_2000_sdpd[535];

    auto g_yz_0_0_0_0_zz_z_xz = buffer_2000_sdpd[536];

    auto g_yz_0_0_0_0_zz_z_yy = buffer_2000_sdpd[537];

    auto g_yz_0_0_0_0_zz_z_yz = buffer_2000_sdpd[538];

    auto g_yz_0_0_0_0_zz_z_zz = buffer_2000_sdpd[539];

    auto g_zz_0_0_0_0_xx_x_xx = buffer_2000_sdpd[540];

    auto g_zz_0_0_0_0_xx_x_xy = buffer_2000_sdpd[541];

    auto g_zz_0_0_0_0_xx_x_xz = buffer_2000_sdpd[542];

    auto g_zz_0_0_0_0_xx_x_yy = buffer_2000_sdpd[543];

    auto g_zz_0_0_0_0_xx_x_yz = buffer_2000_sdpd[544];

    auto g_zz_0_0_0_0_xx_x_zz = buffer_2000_sdpd[545];

    auto g_zz_0_0_0_0_xx_y_xx = buffer_2000_sdpd[546];

    auto g_zz_0_0_0_0_xx_y_xy = buffer_2000_sdpd[547];

    auto g_zz_0_0_0_0_xx_y_xz = buffer_2000_sdpd[548];

    auto g_zz_0_0_0_0_xx_y_yy = buffer_2000_sdpd[549];

    auto g_zz_0_0_0_0_xx_y_yz = buffer_2000_sdpd[550];

    auto g_zz_0_0_0_0_xx_y_zz = buffer_2000_sdpd[551];

    auto g_zz_0_0_0_0_xx_z_xx = buffer_2000_sdpd[552];

    auto g_zz_0_0_0_0_xx_z_xy = buffer_2000_sdpd[553];

    auto g_zz_0_0_0_0_xx_z_xz = buffer_2000_sdpd[554];

    auto g_zz_0_0_0_0_xx_z_yy = buffer_2000_sdpd[555];

    auto g_zz_0_0_0_0_xx_z_yz = buffer_2000_sdpd[556];

    auto g_zz_0_0_0_0_xx_z_zz = buffer_2000_sdpd[557];

    auto g_zz_0_0_0_0_xy_x_xx = buffer_2000_sdpd[558];

    auto g_zz_0_0_0_0_xy_x_xy = buffer_2000_sdpd[559];

    auto g_zz_0_0_0_0_xy_x_xz = buffer_2000_sdpd[560];

    auto g_zz_0_0_0_0_xy_x_yy = buffer_2000_sdpd[561];

    auto g_zz_0_0_0_0_xy_x_yz = buffer_2000_sdpd[562];

    auto g_zz_0_0_0_0_xy_x_zz = buffer_2000_sdpd[563];

    auto g_zz_0_0_0_0_xy_y_xx = buffer_2000_sdpd[564];

    auto g_zz_0_0_0_0_xy_y_xy = buffer_2000_sdpd[565];

    auto g_zz_0_0_0_0_xy_y_xz = buffer_2000_sdpd[566];

    auto g_zz_0_0_0_0_xy_y_yy = buffer_2000_sdpd[567];

    auto g_zz_0_0_0_0_xy_y_yz = buffer_2000_sdpd[568];

    auto g_zz_0_0_0_0_xy_y_zz = buffer_2000_sdpd[569];

    auto g_zz_0_0_0_0_xy_z_xx = buffer_2000_sdpd[570];

    auto g_zz_0_0_0_0_xy_z_xy = buffer_2000_sdpd[571];

    auto g_zz_0_0_0_0_xy_z_xz = buffer_2000_sdpd[572];

    auto g_zz_0_0_0_0_xy_z_yy = buffer_2000_sdpd[573];

    auto g_zz_0_0_0_0_xy_z_yz = buffer_2000_sdpd[574];

    auto g_zz_0_0_0_0_xy_z_zz = buffer_2000_sdpd[575];

    auto g_zz_0_0_0_0_xz_x_xx = buffer_2000_sdpd[576];

    auto g_zz_0_0_0_0_xz_x_xy = buffer_2000_sdpd[577];

    auto g_zz_0_0_0_0_xz_x_xz = buffer_2000_sdpd[578];

    auto g_zz_0_0_0_0_xz_x_yy = buffer_2000_sdpd[579];

    auto g_zz_0_0_0_0_xz_x_yz = buffer_2000_sdpd[580];

    auto g_zz_0_0_0_0_xz_x_zz = buffer_2000_sdpd[581];

    auto g_zz_0_0_0_0_xz_y_xx = buffer_2000_sdpd[582];

    auto g_zz_0_0_0_0_xz_y_xy = buffer_2000_sdpd[583];

    auto g_zz_0_0_0_0_xz_y_xz = buffer_2000_sdpd[584];

    auto g_zz_0_0_0_0_xz_y_yy = buffer_2000_sdpd[585];

    auto g_zz_0_0_0_0_xz_y_yz = buffer_2000_sdpd[586];

    auto g_zz_0_0_0_0_xz_y_zz = buffer_2000_sdpd[587];

    auto g_zz_0_0_0_0_xz_z_xx = buffer_2000_sdpd[588];

    auto g_zz_0_0_0_0_xz_z_xy = buffer_2000_sdpd[589];

    auto g_zz_0_0_0_0_xz_z_xz = buffer_2000_sdpd[590];

    auto g_zz_0_0_0_0_xz_z_yy = buffer_2000_sdpd[591];

    auto g_zz_0_0_0_0_xz_z_yz = buffer_2000_sdpd[592];

    auto g_zz_0_0_0_0_xz_z_zz = buffer_2000_sdpd[593];

    auto g_zz_0_0_0_0_yy_x_xx = buffer_2000_sdpd[594];

    auto g_zz_0_0_0_0_yy_x_xy = buffer_2000_sdpd[595];

    auto g_zz_0_0_0_0_yy_x_xz = buffer_2000_sdpd[596];

    auto g_zz_0_0_0_0_yy_x_yy = buffer_2000_sdpd[597];

    auto g_zz_0_0_0_0_yy_x_yz = buffer_2000_sdpd[598];

    auto g_zz_0_0_0_0_yy_x_zz = buffer_2000_sdpd[599];

    auto g_zz_0_0_0_0_yy_y_xx = buffer_2000_sdpd[600];

    auto g_zz_0_0_0_0_yy_y_xy = buffer_2000_sdpd[601];

    auto g_zz_0_0_0_0_yy_y_xz = buffer_2000_sdpd[602];

    auto g_zz_0_0_0_0_yy_y_yy = buffer_2000_sdpd[603];

    auto g_zz_0_0_0_0_yy_y_yz = buffer_2000_sdpd[604];

    auto g_zz_0_0_0_0_yy_y_zz = buffer_2000_sdpd[605];

    auto g_zz_0_0_0_0_yy_z_xx = buffer_2000_sdpd[606];

    auto g_zz_0_0_0_0_yy_z_xy = buffer_2000_sdpd[607];

    auto g_zz_0_0_0_0_yy_z_xz = buffer_2000_sdpd[608];

    auto g_zz_0_0_0_0_yy_z_yy = buffer_2000_sdpd[609];

    auto g_zz_0_0_0_0_yy_z_yz = buffer_2000_sdpd[610];

    auto g_zz_0_0_0_0_yy_z_zz = buffer_2000_sdpd[611];

    auto g_zz_0_0_0_0_yz_x_xx = buffer_2000_sdpd[612];

    auto g_zz_0_0_0_0_yz_x_xy = buffer_2000_sdpd[613];

    auto g_zz_0_0_0_0_yz_x_xz = buffer_2000_sdpd[614];

    auto g_zz_0_0_0_0_yz_x_yy = buffer_2000_sdpd[615];

    auto g_zz_0_0_0_0_yz_x_yz = buffer_2000_sdpd[616];

    auto g_zz_0_0_0_0_yz_x_zz = buffer_2000_sdpd[617];

    auto g_zz_0_0_0_0_yz_y_xx = buffer_2000_sdpd[618];

    auto g_zz_0_0_0_0_yz_y_xy = buffer_2000_sdpd[619];

    auto g_zz_0_0_0_0_yz_y_xz = buffer_2000_sdpd[620];

    auto g_zz_0_0_0_0_yz_y_yy = buffer_2000_sdpd[621];

    auto g_zz_0_0_0_0_yz_y_yz = buffer_2000_sdpd[622];

    auto g_zz_0_0_0_0_yz_y_zz = buffer_2000_sdpd[623];

    auto g_zz_0_0_0_0_yz_z_xx = buffer_2000_sdpd[624];

    auto g_zz_0_0_0_0_yz_z_xy = buffer_2000_sdpd[625];

    auto g_zz_0_0_0_0_yz_z_xz = buffer_2000_sdpd[626];

    auto g_zz_0_0_0_0_yz_z_yy = buffer_2000_sdpd[627];

    auto g_zz_0_0_0_0_yz_z_yz = buffer_2000_sdpd[628];

    auto g_zz_0_0_0_0_yz_z_zz = buffer_2000_sdpd[629];

    auto g_zz_0_0_0_0_zz_x_xx = buffer_2000_sdpd[630];

    auto g_zz_0_0_0_0_zz_x_xy = buffer_2000_sdpd[631];

    auto g_zz_0_0_0_0_zz_x_xz = buffer_2000_sdpd[632];

    auto g_zz_0_0_0_0_zz_x_yy = buffer_2000_sdpd[633];

    auto g_zz_0_0_0_0_zz_x_yz = buffer_2000_sdpd[634];

    auto g_zz_0_0_0_0_zz_x_zz = buffer_2000_sdpd[635];

    auto g_zz_0_0_0_0_zz_y_xx = buffer_2000_sdpd[636];

    auto g_zz_0_0_0_0_zz_y_xy = buffer_2000_sdpd[637];

    auto g_zz_0_0_0_0_zz_y_xz = buffer_2000_sdpd[638];

    auto g_zz_0_0_0_0_zz_y_yy = buffer_2000_sdpd[639];

    auto g_zz_0_0_0_0_zz_y_yz = buffer_2000_sdpd[640];

    auto g_zz_0_0_0_0_zz_y_zz = buffer_2000_sdpd[641];

    auto g_zz_0_0_0_0_zz_z_xx = buffer_2000_sdpd[642];

    auto g_zz_0_0_0_0_zz_z_xy = buffer_2000_sdpd[643];

    auto g_zz_0_0_0_0_zz_z_xz = buffer_2000_sdpd[644];

    auto g_zz_0_0_0_0_zz_z_yy = buffer_2000_sdpd[645];

    auto g_zz_0_0_0_0_zz_z_yz = buffer_2000_sdpd[646];

    auto g_zz_0_0_0_0_zz_z_zz = buffer_2000_sdpd[647];

    // integrals block (0-6)

    #pragma omp simd aligned(g_0_xx_x_xx, g_0_xx_x_xy, g_0_xx_x_xz, g_0_xx_x_yy, g_0_xx_x_yz, g_0_xx_x_zz, g_xx_0_0_0_0_xx_x_xx, g_xx_0_0_0_0_xx_x_xy, g_xx_0_0_0_0_xx_x_xz, g_xx_0_0_0_0_xx_x_yy, g_xx_0_0_0_0_xx_x_yz, g_xx_0_0_0_0_xx_x_zz, g_xx_xx_x_xx, g_xx_xx_x_xy, g_xx_xx_x_xz, g_xx_xx_x_yy, g_xx_xx_x_yz, g_xx_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_xx_x_xx[i] = -2.0 * g_0_xx_x_xx[i] * a_exp + 4.0 * g_xx_xx_x_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_x_xy[i] = -2.0 * g_0_xx_x_xy[i] * a_exp + 4.0 * g_xx_xx_x_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_x_xz[i] = -2.0 * g_0_xx_x_xz[i] * a_exp + 4.0 * g_xx_xx_x_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_x_yy[i] = -2.0 * g_0_xx_x_yy[i] * a_exp + 4.0 * g_xx_xx_x_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_x_yz[i] = -2.0 * g_0_xx_x_yz[i] * a_exp + 4.0 * g_xx_xx_x_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_x_zz[i] = -2.0 * g_0_xx_x_zz[i] * a_exp + 4.0 * g_xx_xx_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_0_xx_y_xx, g_0_xx_y_xy, g_0_xx_y_xz, g_0_xx_y_yy, g_0_xx_y_yz, g_0_xx_y_zz, g_xx_0_0_0_0_xx_y_xx, g_xx_0_0_0_0_xx_y_xy, g_xx_0_0_0_0_xx_y_xz, g_xx_0_0_0_0_xx_y_yy, g_xx_0_0_0_0_xx_y_yz, g_xx_0_0_0_0_xx_y_zz, g_xx_xx_y_xx, g_xx_xx_y_xy, g_xx_xx_y_xz, g_xx_xx_y_yy, g_xx_xx_y_yz, g_xx_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_xx_y_xx[i] = -2.0 * g_0_xx_y_xx[i] * a_exp + 4.0 * g_xx_xx_y_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_y_xy[i] = -2.0 * g_0_xx_y_xy[i] * a_exp + 4.0 * g_xx_xx_y_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_y_xz[i] = -2.0 * g_0_xx_y_xz[i] * a_exp + 4.0 * g_xx_xx_y_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_y_yy[i] = -2.0 * g_0_xx_y_yy[i] * a_exp + 4.0 * g_xx_xx_y_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_y_yz[i] = -2.0 * g_0_xx_y_yz[i] * a_exp + 4.0 * g_xx_xx_y_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_y_zz[i] = -2.0 * g_0_xx_y_zz[i] * a_exp + 4.0 * g_xx_xx_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_0_xx_z_xx, g_0_xx_z_xy, g_0_xx_z_xz, g_0_xx_z_yy, g_0_xx_z_yz, g_0_xx_z_zz, g_xx_0_0_0_0_xx_z_xx, g_xx_0_0_0_0_xx_z_xy, g_xx_0_0_0_0_xx_z_xz, g_xx_0_0_0_0_xx_z_yy, g_xx_0_0_0_0_xx_z_yz, g_xx_0_0_0_0_xx_z_zz, g_xx_xx_z_xx, g_xx_xx_z_xy, g_xx_xx_z_xz, g_xx_xx_z_yy, g_xx_xx_z_yz, g_xx_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_xx_z_xx[i] = -2.0 * g_0_xx_z_xx[i] * a_exp + 4.0 * g_xx_xx_z_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_z_xy[i] = -2.0 * g_0_xx_z_xy[i] * a_exp + 4.0 * g_xx_xx_z_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_z_xz[i] = -2.0 * g_0_xx_z_xz[i] * a_exp + 4.0 * g_xx_xx_z_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_z_yy[i] = -2.0 * g_0_xx_z_yy[i] * a_exp + 4.0 * g_xx_xx_z_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_z_yz[i] = -2.0 * g_0_xx_z_yz[i] * a_exp + 4.0 * g_xx_xx_z_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_z_zz[i] = -2.0 * g_0_xx_z_zz[i] * a_exp + 4.0 * g_xx_xx_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_0_xy_x_xx, g_0_xy_x_xy, g_0_xy_x_xz, g_0_xy_x_yy, g_0_xy_x_yz, g_0_xy_x_zz, g_xx_0_0_0_0_xy_x_xx, g_xx_0_0_0_0_xy_x_xy, g_xx_0_0_0_0_xy_x_xz, g_xx_0_0_0_0_xy_x_yy, g_xx_0_0_0_0_xy_x_yz, g_xx_0_0_0_0_xy_x_zz, g_xx_xy_x_xx, g_xx_xy_x_xy, g_xx_xy_x_xz, g_xx_xy_x_yy, g_xx_xy_x_yz, g_xx_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_xy_x_xx[i] = -2.0 * g_0_xy_x_xx[i] * a_exp + 4.0 * g_xx_xy_x_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_x_xy[i] = -2.0 * g_0_xy_x_xy[i] * a_exp + 4.0 * g_xx_xy_x_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_x_xz[i] = -2.0 * g_0_xy_x_xz[i] * a_exp + 4.0 * g_xx_xy_x_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_x_yy[i] = -2.0 * g_0_xy_x_yy[i] * a_exp + 4.0 * g_xx_xy_x_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_x_yz[i] = -2.0 * g_0_xy_x_yz[i] * a_exp + 4.0 * g_xx_xy_x_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_x_zz[i] = -2.0 * g_0_xy_x_zz[i] * a_exp + 4.0 * g_xx_xy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_0_xy_y_xx, g_0_xy_y_xy, g_0_xy_y_xz, g_0_xy_y_yy, g_0_xy_y_yz, g_0_xy_y_zz, g_xx_0_0_0_0_xy_y_xx, g_xx_0_0_0_0_xy_y_xy, g_xx_0_0_0_0_xy_y_xz, g_xx_0_0_0_0_xy_y_yy, g_xx_0_0_0_0_xy_y_yz, g_xx_0_0_0_0_xy_y_zz, g_xx_xy_y_xx, g_xx_xy_y_xy, g_xx_xy_y_xz, g_xx_xy_y_yy, g_xx_xy_y_yz, g_xx_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_xy_y_xx[i] = -2.0 * g_0_xy_y_xx[i] * a_exp + 4.0 * g_xx_xy_y_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_y_xy[i] = -2.0 * g_0_xy_y_xy[i] * a_exp + 4.0 * g_xx_xy_y_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_y_xz[i] = -2.0 * g_0_xy_y_xz[i] * a_exp + 4.0 * g_xx_xy_y_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_y_yy[i] = -2.0 * g_0_xy_y_yy[i] * a_exp + 4.0 * g_xx_xy_y_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_y_yz[i] = -2.0 * g_0_xy_y_yz[i] * a_exp + 4.0 * g_xx_xy_y_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_y_zz[i] = -2.0 * g_0_xy_y_zz[i] * a_exp + 4.0 * g_xx_xy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_0_xy_z_xx, g_0_xy_z_xy, g_0_xy_z_xz, g_0_xy_z_yy, g_0_xy_z_yz, g_0_xy_z_zz, g_xx_0_0_0_0_xy_z_xx, g_xx_0_0_0_0_xy_z_xy, g_xx_0_0_0_0_xy_z_xz, g_xx_0_0_0_0_xy_z_yy, g_xx_0_0_0_0_xy_z_yz, g_xx_0_0_0_0_xy_z_zz, g_xx_xy_z_xx, g_xx_xy_z_xy, g_xx_xy_z_xz, g_xx_xy_z_yy, g_xx_xy_z_yz, g_xx_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_xy_z_xx[i] = -2.0 * g_0_xy_z_xx[i] * a_exp + 4.0 * g_xx_xy_z_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_z_xy[i] = -2.0 * g_0_xy_z_xy[i] * a_exp + 4.0 * g_xx_xy_z_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_z_xz[i] = -2.0 * g_0_xy_z_xz[i] * a_exp + 4.0 * g_xx_xy_z_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_z_yy[i] = -2.0 * g_0_xy_z_yy[i] * a_exp + 4.0 * g_xx_xy_z_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_z_yz[i] = -2.0 * g_0_xy_z_yz[i] * a_exp + 4.0 * g_xx_xy_z_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_z_zz[i] = -2.0 * g_0_xy_z_zz[i] * a_exp + 4.0 * g_xx_xy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_0_xz_x_xx, g_0_xz_x_xy, g_0_xz_x_xz, g_0_xz_x_yy, g_0_xz_x_yz, g_0_xz_x_zz, g_xx_0_0_0_0_xz_x_xx, g_xx_0_0_0_0_xz_x_xy, g_xx_0_0_0_0_xz_x_xz, g_xx_0_0_0_0_xz_x_yy, g_xx_0_0_0_0_xz_x_yz, g_xx_0_0_0_0_xz_x_zz, g_xx_xz_x_xx, g_xx_xz_x_xy, g_xx_xz_x_xz, g_xx_xz_x_yy, g_xx_xz_x_yz, g_xx_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_xz_x_xx[i] = -2.0 * g_0_xz_x_xx[i] * a_exp + 4.0 * g_xx_xz_x_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_x_xy[i] = -2.0 * g_0_xz_x_xy[i] * a_exp + 4.0 * g_xx_xz_x_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_x_xz[i] = -2.0 * g_0_xz_x_xz[i] * a_exp + 4.0 * g_xx_xz_x_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_x_yy[i] = -2.0 * g_0_xz_x_yy[i] * a_exp + 4.0 * g_xx_xz_x_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_x_yz[i] = -2.0 * g_0_xz_x_yz[i] * a_exp + 4.0 * g_xx_xz_x_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_x_zz[i] = -2.0 * g_0_xz_x_zz[i] * a_exp + 4.0 * g_xx_xz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_0_xz_y_xx, g_0_xz_y_xy, g_0_xz_y_xz, g_0_xz_y_yy, g_0_xz_y_yz, g_0_xz_y_zz, g_xx_0_0_0_0_xz_y_xx, g_xx_0_0_0_0_xz_y_xy, g_xx_0_0_0_0_xz_y_xz, g_xx_0_0_0_0_xz_y_yy, g_xx_0_0_0_0_xz_y_yz, g_xx_0_0_0_0_xz_y_zz, g_xx_xz_y_xx, g_xx_xz_y_xy, g_xx_xz_y_xz, g_xx_xz_y_yy, g_xx_xz_y_yz, g_xx_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_xz_y_xx[i] = -2.0 * g_0_xz_y_xx[i] * a_exp + 4.0 * g_xx_xz_y_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_y_xy[i] = -2.0 * g_0_xz_y_xy[i] * a_exp + 4.0 * g_xx_xz_y_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_y_xz[i] = -2.0 * g_0_xz_y_xz[i] * a_exp + 4.0 * g_xx_xz_y_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_y_yy[i] = -2.0 * g_0_xz_y_yy[i] * a_exp + 4.0 * g_xx_xz_y_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_y_yz[i] = -2.0 * g_0_xz_y_yz[i] * a_exp + 4.0 * g_xx_xz_y_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_y_zz[i] = -2.0 * g_0_xz_y_zz[i] * a_exp + 4.0 * g_xx_xz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_0_xz_z_xx, g_0_xz_z_xy, g_0_xz_z_xz, g_0_xz_z_yy, g_0_xz_z_yz, g_0_xz_z_zz, g_xx_0_0_0_0_xz_z_xx, g_xx_0_0_0_0_xz_z_xy, g_xx_0_0_0_0_xz_z_xz, g_xx_0_0_0_0_xz_z_yy, g_xx_0_0_0_0_xz_z_yz, g_xx_0_0_0_0_xz_z_zz, g_xx_xz_z_xx, g_xx_xz_z_xy, g_xx_xz_z_xz, g_xx_xz_z_yy, g_xx_xz_z_yz, g_xx_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_xz_z_xx[i] = -2.0 * g_0_xz_z_xx[i] * a_exp + 4.0 * g_xx_xz_z_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_z_xy[i] = -2.0 * g_0_xz_z_xy[i] * a_exp + 4.0 * g_xx_xz_z_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_z_xz[i] = -2.0 * g_0_xz_z_xz[i] * a_exp + 4.0 * g_xx_xz_z_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_z_yy[i] = -2.0 * g_0_xz_z_yy[i] * a_exp + 4.0 * g_xx_xz_z_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_z_yz[i] = -2.0 * g_0_xz_z_yz[i] * a_exp + 4.0 * g_xx_xz_z_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_z_zz[i] = -2.0 * g_0_xz_z_zz[i] * a_exp + 4.0 * g_xx_xz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_0_yy_x_xx, g_0_yy_x_xy, g_0_yy_x_xz, g_0_yy_x_yy, g_0_yy_x_yz, g_0_yy_x_zz, g_xx_0_0_0_0_yy_x_xx, g_xx_0_0_0_0_yy_x_xy, g_xx_0_0_0_0_yy_x_xz, g_xx_0_0_0_0_yy_x_yy, g_xx_0_0_0_0_yy_x_yz, g_xx_0_0_0_0_yy_x_zz, g_xx_yy_x_xx, g_xx_yy_x_xy, g_xx_yy_x_xz, g_xx_yy_x_yy, g_xx_yy_x_yz, g_xx_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_yy_x_xx[i] = -2.0 * g_0_yy_x_xx[i] * a_exp + 4.0 * g_xx_yy_x_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_x_xy[i] = -2.0 * g_0_yy_x_xy[i] * a_exp + 4.0 * g_xx_yy_x_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_x_xz[i] = -2.0 * g_0_yy_x_xz[i] * a_exp + 4.0 * g_xx_yy_x_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_x_yy[i] = -2.0 * g_0_yy_x_yy[i] * a_exp + 4.0 * g_xx_yy_x_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_x_yz[i] = -2.0 * g_0_yy_x_yz[i] * a_exp + 4.0 * g_xx_yy_x_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_x_zz[i] = -2.0 * g_0_yy_x_zz[i] * a_exp + 4.0 * g_xx_yy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_0_yy_y_xx, g_0_yy_y_xy, g_0_yy_y_xz, g_0_yy_y_yy, g_0_yy_y_yz, g_0_yy_y_zz, g_xx_0_0_0_0_yy_y_xx, g_xx_0_0_0_0_yy_y_xy, g_xx_0_0_0_0_yy_y_xz, g_xx_0_0_0_0_yy_y_yy, g_xx_0_0_0_0_yy_y_yz, g_xx_0_0_0_0_yy_y_zz, g_xx_yy_y_xx, g_xx_yy_y_xy, g_xx_yy_y_xz, g_xx_yy_y_yy, g_xx_yy_y_yz, g_xx_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_yy_y_xx[i] = -2.0 * g_0_yy_y_xx[i] * a_exp + 4.0 * g_xx_yy_y_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_y_xy[i] = -2.0 * g_0_yy_y_xy[i] * a_exp + 4.0 * g_xx_yy_y_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_y_xz[i] = -2.0 * g_0_yy_y_xz[i] * a_exp + 4.0 * g_xx_yy_y_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_y_yy[i] = -2.0 * g_0_yy_y_yy[i] * a_exp + 4.0 * g_xx_yy_y_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_y_yz[i] = -2.0 * g_0_yy_y_yz[i] * a_exp + 4.0 * g_xx_yy_y_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_y_zz[i] = -2.0 * g_0_yy_y_zz[i] * a_exp + 4.0 * g_xx_yy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_0_yy_z_xx, g_0_yy_z_xy, g_0_yy_z_xz, g_0_yy_z_yy, g_0_yy_z_yz, g_0_yy_z_zz, g_xx_0_0_0_0_yy_z_xx, g_xx_0_0_0_0_yy_z_xy, g_xx_0_0_0_0_yy_z_xz, g_xx_0_0_0_0_yy_z_yy, g_xx_0_0_0_0_yy_z_yz, g_xx_0_0_0_0_yy_z_zz, g_xx_yy_z_xx, g_xx_yy_z_xy, g_xx_yy_z_xz, g_xx_yy_z_yy, g_xx_yy_z_yz, g_xx_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_yy_z_xx[i] = -2.0 * g_0_yy_z_xx[i] * a_exp + 4.0 * g_xx_yy_z_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_z_xy[i] = -2.0 * g_0_yy_z_xy[i] * a_exp + 4.0 * g_xx_yy_z_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_z_xz[i] = -2.0 * g_0_yy_z_xz[i] * a_exp + 4.0 * g_xx_yy_z_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_z_yy[i] = -2.0 * g_0_yy_z_yy[i] * a_exp + 4.0 * g_xx_yy_z_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_z_yz[i] = -2.0 * g_0_yy_z_yz[i] * a_exp + 4.0 * g_xx_yy_z_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_z_zz[i] = -2.0 * g_0_yy_z_zz[i] * a_exp + 4.0 * g_xx_yy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_0_yz_x_xx, g_0_yz_x_xy, g_0_yz_x_xz, g_0_yz_x_yy, g_0_yz_x_yz, g_0_yz_x_zz, g_xx_0_0_0_0_yz_x_xx, g_xx_0_0_0_0_yz_x_xy, g_xx_0_0_0_0_yz_x_xz, g_xx_0_0_0_0_yz_x_yy, g_xx_0_0_0_0_yz_x_yz, g_xx_0_0_0_0_yz_x_zz, g_xx_yz_x_xx, g_xx_yz_x_xy, g_xx_yz_x_xz, g_xx_yz_x_yy, g_xx_yz_x_yz, g_xx_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_yz_x_xx[i] = -2.0 * g_0_yz_x_xx[i] * a_exp + 4.0 * g_xx_yz_x_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_x_xy[i] = -2.0 * g_0_yz_x_xy[i] * a_exp + 4.0 * g_xx_yz_x_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_x_xz[i] = -2.0 * g_0_yz_x_xz[i] * a_exp + 4.0 * g_xx_yz_x_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_x_yy[i] = -2.0 * g_0_yz_x_yy[i] * a_exp + 4.0 * g_xx_yz_x_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_x_yz[i] = -2.0 * g_0_yz_x_yz[i] * a_exp + 4.0 * g_xx_yz_x_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_x_zz[i] = -2.0 * g_0_yz_x_zz[i] * a_exp + 4.0 * g_xx_yz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_0_yz_y_xx, g_0_yz_y_xy, g_0_yz_y_xz, g_0_yz_y_yy, g_0_yz_y_yz, g_0_yz_y_zz, g_xx_0_0_0_0_yz_y_xx, g_xx_0_0_0_0_yz_y_xy, g_xx_0_0_0_0_yz_y_xz, g_xx_0_0_0_0_yz_y_yy, g_xx_0_0_0_0_yz_y_yz, g_xx_0_0_0_0_yz_y_zz, g_xx_yz_y_xx, g_xx_yz_y_xy, g_xx_yz_y_xz, g_xx_yz_y_yy, g_xx_yz_y_yz, g_xx_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_yz_y_xx[i] = -2.0 * g_0_yz_y_xx[i] * a_exp + 4.0 * g_xx_yz_y_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_y_xy[i] = -2.0 * g_0_yz_y_xy[i] * a_exp + 4.0 * g_xx_yz_y_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_y_xz[i] = -2.0 * g_0_yz_y_xz[i] * a_exp + 4.0 * g_xx_yz_y_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_y_yy[i] = -2.0 * g_0_yz_y_yy[i] * a_exp + 4.0 * g_xx_yz_y_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_y_yz[i] = -2.0 * g_0_yz_y_yz[i] * a_exp + 4.0 * g_xx_yz_y_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_y_zz[i] = -2.0 * g_0_yz_y_zz[i] * a_exp + 4.0 * g_xx_yz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_0_yz_z_xx, g_0_yz_z_xy, g_0_yz_z_xz, g_0_yz_z_yy, g_0_yz_z_yz, g_0_yz_z_zz, g_xx_0_0_0_0_yz_z_xx, g_xx_0_0_0_0_yz_z_xy, g_xx_0_0_0_0_yz_z_xz, g_xx_0_0_0_0_yz_z_yy, g_xx_0_0_0_0_yz_z_yz, g_xx_0_0_0_0_yz_z_zz, g_xx_yz_z_xx, g_xx_yz_z_xy, g_xx_yz_z_xz, g_xx_yz_z_yy, g_xx_yz_z_yz, g_xx_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_yz_z_xx[i] = -2.0 * g_0_yz_z_xx[i] * a_exp + 4.0 * g_xx_yz_z_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_z_xy[i] = -2.0 * g_0_yz_z_xy[i] * a_exp + 4.0 * g_xx_yz_z_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_z_xz[i] = -2.0 * g_0_yz_z_xz[i] * a_exp + 4.0 * g_xx_yz_z_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_z_yy[i] = -2.0 * g_0_yz_z_yy[i] * a_exp + 4.0 * g_xx_yz_z_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_z_yz[i] = -2.0 * g_0_yz_z_yz[i] * a_exp + 4.0 * g_xx_yz_z_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_z_zz[i] = -2.0 * g_0_yz_z_zz[i] * a_exp + 4.0 * g_xx_yz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_0_zz_x_xx, g_0_zz_x_xy, g_0_zz_x_xz, g_0_zz_x_yy, g_0_zz_x_yz, g_0_zz_x_zz, g_xx_0_0_0_0_zz_x_xx, g_xx_0_0_0_0_zz_x_xy, g_xx_0_0_0_0_zz_x_xz, g_xx_0_0_0_0_zz_x_yy, g_xx_0_0_0_0_zz_x_yz, g_xx_0_0_0_0_zz_x_zz, g_xx_zz_x_xx, g_xx_zz_x_xy, g_xx_zz_x_xz, g_xx_zz_x_yy, g_xx_zz_x_yz, g_xx_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_zz_x_xx[i] = -2.0 * g_0_zz_x_xx[i] * a_exp + 4.0 * g_xx_zz_x_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_x_xy[i] = -2.0 * g_0_zz_x_xy[i] * a_exp + 4.0 * g_xx_zz_x_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_x_xz[i] = -2.0 * g_0_zz_x_xz[i] * a_exp + 4.0 * g_xx_zz_x_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_x_yy[i] = -2.0 * g_0_zz_x_yy[i] * a_exp + 4.0 * g_xx_zz_x_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_x_yz[i] = -2.0 * g_0_zz_x_yz[i] * a_exp + 4.0 * g_xx_zz_x_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_x_zz[i] = -2.0 * g_0_zz_x_zz[i] * a_exp + 4.0 * g_xx_zz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_0_zz_y_xx, g_0_zz_y_xy, g_0_zz_y_xz, g_0_zz_y_yy, g_0_zz_y_yz, g_0_zz_y_zz, g_xx_0_0_0_0_zz_y_xx, g_xx_0_0_0_0_zz_y_xy, g_xx_0_0_0_0_zz_y_xz, g_xx_0_0_0_0_zz_y_yy, g_xx_0_0_0_0_zz_y_yz, g_xx_0_0_0_0_zz_y_zz, g_xx_zz_y_xx, g_xx_zz_y_xy, g_xx_zz_y_xz, g_xx_zz_y_yy, g_xx_zz_y_yz, g_xx_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_zz_y_xx[i] = -2.0 * g_0_zz_y_xx[i] * a_exp + 4.0 * g_xx_zz_y_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_y_xy[i] = -2.0 * g_0_zz_y_xy[i] * a_exp + 4.0 * g_xx_zz_y_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_y_xz[i] = -2.0 * g_0_zz_y_xz[i] * a_exp + 4.0 * g_xx_zz_y_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_y_yy[i] = -2.0 * g_0_zz_y_yy[i] * a_exp + 4.0 * g_xx_zz_y_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_y_yz[i] = -2.0 * g_0_zz_y_yz[i] * a_exp + 4.0 * g_xx_zz_y_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_y_zz[i] = -2.0 * g_0_zz_y_zz[i] * a_exp + 4.0 * g_xx_zz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_0_zz_z_xx, g_0_zz_z_xy, g_0_zz_z_xz, g_0_zz_z_yy, g_0_zz_z_yz, g_0_zz_z_zz, g_xx_0_0_0_0_zz_z_xx, g_xx_0_0_0_0_zz_z_xy, g_xx_0_0_0_0_zz_z_xz, g_xx_0_0_0_0_zz_z_yy, g_xx_0_0_0_0_zz_z_yz, g_xx_0_0_0_0_zz_z_zz, g_xx_zz_z_xx, g_xx_zz_z_xy, g_xx_zz_z_xz, g_xx_zz_z_yy, g_xx_zz_z_yz, g_xx_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_zz_z_xx[i] = -2.0 * g_0_zz_z_xx[i] * a_exp + 4.0 * g_xx_zz_z_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_z_xy[i] = -2.0 * g_0_zz_z_xy[i] * a_exp + 4.0 * g_xx_zz_z_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_z_xz[i] = -2.0 * g_0_zz_z_xz[i] * a_exp + 4.0 * g_xx_zz_z_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_z_yy[i] = -2.0 * g_0_zz_z_yy[i] * a_exp + 4.0 * g_xx_zz_z_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_z_yz[i] = -2.0 * g_0_zz_z_yz[i] * a_exp + 4.0 * g_xx_zz_z_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_z_zz[i] = -2.0 * g_0_zz_z_zz[i] * a_exp + 4.0 * g_xx_zz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_xy_0_0_0_0_xx_x_xx, g_xy_0_0_0_0_xx_x_xy, g_xy_0_0_0_0_xx_x_xz, g_xy_0_0_0_0_xx_x_yy, g_xy_0_0_0_0_xx_x_yz, g_xy_0_0_0_0_xx_x_zz, g_xy_xx_x_xx, g_xy_xx_x_xy, g_xy_xx_x_xz, g_xy_xx_x_yy, g_xy_xx_x_yz, g_xy_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_xx_x_xx[i] = 4.0 * g_xy_xx_x_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_x_xy[i] = 4.0 * g_xy_xx_x_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_x_xz[i] = 4.0 * g_xy_xx_x_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_x_yy[i] = 4.0 * g_xy_xx_x_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_x_yz[i] = 4.0 * g_xy_xx_x_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_x_zz[i] = 4.0 * g_xy_xx_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_xy_0_0_0_0_xx_y_xx, g_xy_0_0_0_0_xx_y_xy, g_xy_0_0_0_0_xx_y_xz, g_xy_0_0_0_0_xx_y_yy, g_xy_0_0_0_0_xx_y_yz, g_xy_0_0_0_0_xx_y_zz, g_xy_xx_y_xx, g_xy_xx_y_xy, g_xy_xx_y_xz, g_xy_xx_y_yy, g_xy_xx_y_yz, g_xy_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_xx_y_xx[i] = 4.0 * g_xy_xx_y_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_y_xy[i] = 4.0 * g_xy_xx_y_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_y_xz[i] = 4.0 * g_xy_xx_y_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_y_yy[i] = 4.0 * g_xy_xx_y_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_y_yz[i] = 4.0 * g_xy_xx_y_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_y_zz[i] = 4.0 * g_xy_xx_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_xy_0_0_0_0_xx_z_xx, g_xy_0_0_0_0_xx_z_xy, g_xy_0_0_0_0_xx_z_xz, g_xy_0_0_0_0_xx_z_yy, g_xy_0_0_0_0_xx_z_yz, g_xy_0_0_0_0_xx_z_zz, g_xy_xx_z_xx, g_xy_xx_z_xy, g_xy_xx_z_xz, g_xy_xx_z_yy, g_xy_xx_z_yz, g_xy_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_xx_z_xx[i] = 4.0 * g_xy_xx_z_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_z_xy[i] = 4.0 * g_xy_xx_z_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_z_xz[i] = 4.0 * g_xy_xx_z_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_z_yy[i] = 4.0 * g_xy_xx_z_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_z_yz[i] = 4.0 * g_xy_xx_z_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_z_zz[i] = 4.0 * g_xy_xx_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_xy_0_0_0_0_xy_x_xx, g_xy_0_0_0_0_xy_x_xy, g_xy_0_0_0_0_xy_x_xz, g_xy_0_0_0_0_xy_x_yy, g_xy_0_0_0_0_xy_x_yz, g_xy_0_0_0_0_xy_x_zz, g_xy_xy_x_xx, g_xy_xy_x_xy, g_xy_xy_x_xz, g_xy_xy_x_yy, g_xy_xy_x_yz, g_xy_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_xy_x_xx[i] = 4.0 * g_xy_xy_x_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_x_xy[i] = 4.0 * g_xy_xy_x_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_x_xz[i] = 4.0 * g_xy_xy_x_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_x_yy[i] = 4.0 * g_xy_xy_x_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_x_yz[i] = 4.0 * g_xy_xy_x_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_x_zz[i] = 4.0 * g_xy_xy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_xy_0_0_0_0_xy_y_xx, g_xy_0_0_0_0_xy_y_xy, g_xy_0_0_0_0_xy_y_xz, g_xy_0_0_0_0_xy_y_yy, g_xy_0_0_0_0_xy_y_yz, g_xy_0_0_0_0_xy_y_zz, g_xy_xy_y_xx, g_xy_xy_y_xy, g_xy_xy_y_xz, g_xy_xy_y_yy, g_xy_xy_y_yz, g_xy_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_xy_y_xx[i] = 4.0 * g_xy_xy_y_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_y_xy[i] = 4.0 * g_xy_xy_y_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_y_xz[i] = 4.0 * g_xy_xy_y_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_y_yy[i] = 4.0 * g_xy_xy_y_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_y_yz[i] = 4.0 * g_xy_xy_y_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_y_zz[i] = 4.0 * g_xy_xy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_xy_0_0_0_0_xy_z_xx, g_xy_0_0_0_0_xy_z_xy, g_xy_0_0_0_0_xy_z_xz, g_xy_0_0_0_0_xy_z_yy, g_xy_0_0_0_0_xy_z_yz, g_xy_0_0_0_0_xy_z_zz, g_xy_xy_z_xx, g_xy_xy_z_xy, g_xy_xy_z_xz, g_xy_xy_z_yy, g_xy_xy_z_yz, g_xy_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_xy_z_xx[i] = 4.0 * g_xy_xy_z_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_z_xy[i] = 4.0 * g_xy_xy_z_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_z_xz[i] = 4.0 * g_xy_xy_z_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_z_yy[i] = 4.0 * g_xy_xy_z_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_z_yz[i] = 4.0 * g_xy_xy_z_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_z_zz[i] = 4.0 * g_xy_xy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_xy_0_0_0_0_xz_x_xx, g_xy_0_0_0_0_xz_x_xy, g_xy_0_0_0_0_xz_x_xz, g_xy_0_0_0_0_xz_x_yy, g_xy_0_0_0_0_xz_x_yz, g_xy_0_0_0_0_xz_x_zz, g_xy_xz_x_xx, g_xy_xz_x_xy, g_xy_xz_x_xz, g_xy_xz_x_yy, g_xy_xz_x_yz, g_xy_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_xz_x_xx[i] = 4.0 * g_xy_xz_x_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_x_xy[i] = 4.0 * g_xy_xz_x_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_x_xz[i] = 4.0 * g_xy_xz_x_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_x_yy[i] = 4.0 * g_xy_xz_x_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_x_yz[i] = 4.0 * g_xy_xz_x_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_x_zz[i] = 4.0 * g_xy_xz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_xy_0_0_0_0_xz_y_xx, g_xy_0_0_0_0_xz_y_xy, g_xy_0_0_0_0_xz_y_xz, g_xy_0_0_0_0_xz_y_yy, g_xy_0_0_0_0_xz_y_yz, g_xy_0_0_0_0_xz_y_zz, g_xy_xz_y_xx, g_xy_xz_y_xy, g_xy_xz_y_xz, g_xy_xz_y_yy, g_xy_xz_y_yz, g_xy_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_xz_y_xx[i] = 4.0 * g_xy_xz_y_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_y_xy[i] = 4.0 * g_xy_xz_y_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_y_xz[i] = 4.0 * g_xy_xz_y_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_y_yy[i] = 4.0 * g_xy_xz_y_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_y_yz[i] = 4.0 * g_xy_xz_y_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_y_zz[i] = 4.0 * g_xy_xz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_xy_0_0_0_0_xz_z_xx, g_xy_0_0_0_0_xz_z_xy, g_xy_0_0_0_0_xz_z_xz, g_xy_0_0_0_0_xz_z_yy, g_xy_0_0_0_0_xz_z_yz, g_xy_0_0_0_0_xz_z_zz, g_xy_xz_z_xx, g_xy_xz_z_xy, g_xy_xz_z_xz, g_xy_xz_z_yy, g_xy_xz_z_yz, g_xy_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_xz_z_xx[i] = 4.0 * g_xy_xz_z_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_z_xy[i] = 4.0 * g_xy_xz_z_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_z_xz[i] = 4.0 * g_xy_xz_z_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_z_yy[i] = 4.0 * g_xy_xz_z_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_z_yz[i] = 4.0 * g_xy_xz_z_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_z_zz[i] = 4.0 * g_xy_xz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_xy_0_0_0_0_yy_x_xx, g_xy_0_0_0_0_yy_x_xy, g_xy_0_0_0_0_yy_x_xz, g_xy_0_0_0_0_yy_x_yy, g_xy_0_0_0_0_yy_x_yz, g_xy_0_0_0_0_yy_x_zz, g_xy_yy_x_xx, g_xy_yy_x_xy, g_xy_yy_x_xz, g_xy_yy_x_yy, g_xy_yy_x_yz, g_xy_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_yy_x_xx[i] = 4.0 * g_xy_yy_x_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_x_xy[i] = 4.0 * g_xy_yy_x_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_x_xz[i] = 4.0 * g_xy_yy_x_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_x_yy[i] = 4.0 * g_xy_yy_x_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_x_yz[i] = 4.0 * g_xy_yy_x_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_x_zz[i] = 4.0 * g_xy_yy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_xy_0_0_0_0_yy_y_xx, g_xy_0_0_0_0_yy_y_xy, g_xy_0_0_0_0_yy_y_xz, g_xy_0_0_0_0_yy_y_yy, g_xy_0_0_0_0_yy_y_yz, g_xy_0_0_0_0_yy_y_zz, g_xy_yy_y_xx, g_xy_yy_y_xy, g_xy_yy_y_xz, g_xy_yy_y_yy, g_xy_yy_y_yz, g_xy_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_yy_y_xx[i] = 4.0 * g_xy_yy_y_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_y_xy[i] = 4.0 * g_xy_yy_y_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_y_xz[i] = 4.0 * g_xy_yy_y_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_y_yy[i] = 4.0 * g_xy_yy_y_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_y_yz[i] = 4.0 * g_xy_yy_y_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_y_zz[i] = 4.0 * g_xy_yy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_xy_0_0_0_0_yy_z_xx, g_xy_0_0_0_0_yy_z_xy, g_xy_0_0_0_0_yy_z_xz, g_xy_0_0_0_0_yy_z_yy, g_xy_0_0_0_0_yy_z_yz, g_xy_0_0_0_0_yy_z_zz, g_xy_yy_z_xx, g_xy_yy_z_xy, g_xy_yy_z_xz, g_xy_yy_z_yy, g_xy_yy_z_yz, g_xy_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_yy_z_xx[i] = 4.0 * g_xy_yy_z_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_z_xy[i] = 4.0 * g_xy_yy_z_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_z_xz[i] = 4.0 * g_xy_yy_z_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_z_yy[i] = 4.0 * g_xy_yy_z_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_z_yz[i] = 4.0 * g_xy_yy_z_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_z_zz[i] = 4.0 * g_xy_yy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_xy_0_0_0_0_yz_x_xx, g_xy_0_0_0_0_yz_x_xy, g_xy_0_0_0_0_yz_x_xz, g_xy_0_0_0_0_yz_x_yy, g_xy_0_0_0_0_yz_x_yz, g_xy_0_0_0_0_yz_x_zz, g_xy_yz_x_xx, g_xy_yz_x_xy, g_xy_yz_x_xz, g_xy_yz_x_yy, g_xy_yz_x_yz, g_xy_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_yz_x_xx[i] = 4.0 * g_xy_yz_x_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_x_xy[i] = 4.0 * g_xy_yz_x_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_x_xz[i] = 4.0 * g_xy_yz_x_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_x_yy[i] = 4.0 * g_xy_yz_x_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_x_yz[i] = 4.0 * g_xy_yz_x_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_x_zz[i] = 4.0 * g_xy_yz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_xy_0_0_0_0_yz_y_xx, g_xy_0_0_0_0_yz_y_xy, g_xy_0_0_0_0_yz_y_xz, g_xy_0_0_0_0_yz_y_yy, g_xy_0_0_0_0_yz_y_yz, g_xy_0_0_0_0_yz_y_zz, g_xy_yz_y_xx, g_xy_yz_y_xy, g_xy_yz_y_xz, g_xy_yz_y_yy, g_xy_yz_y_yz, g_xy_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_yz_y_xx[i] = 4.0 * g_xy_yz_y_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_y_xy[i] = 4.0 * g_xy_yz_y_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_y_xz[i] = 4.0 * g_xy_yz_y_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_y_yy[i] = 4.0 * g_xy_yz_y_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_y_yz[i] = 4.0 * g_xy_yz_y_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_y_zz[i] = 4.0 * g_xy_yz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_xy_0_0_0_0_yz_z_xx, g_xy_0_0_0_0_yz_z_xy, g_xy_0_0_0_0_yz_z_xz, g_xy_0_0_0_0_yz_z_yy, g_xy_0_0_0_0_yz_z_yz, g_xy_0_0_0_0_yz_z_zz, g_xy_yz_z_xx, g_xy_yz_z_xy, g_xy_yz_z_xz, g_xy_yz_z_yy, g_xy_yz_z_yz, g_xy_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_yz_z_xx[i] = 4.0 * g_xy_yz_z_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_z_xy[i] = 4.0 * g_xy_yz_z_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_z_xz[i] = 4.0 * g_xy_yz_z_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_z_yy[i] = 4.0 * g_xy_yz_z_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_z_yz[i] = 4.0 * g_xy_yz_z_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_z_zz[i] = 4.0 * g_xy_yz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_xy_0_0_0_0_zz_x_xx, g_xy_0_0_0_0_zz_x_xy, g_xy_0_0_0_0_zz_x_xz, g_xy_0_0_0_0_zz_x_yy, g_xy_0_0_0_0_zz_x_yz, g_xy_0_0_0_0_zz_x_zz, g_xy_zz_x_xx, g_xy_zz_x_xy, g_xy_zz_x_xz, g_xy_zz_x_yy, g_xy_zz_x_yz, g_xy_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_zz_x_xx[i] = 4.0 * g_xy_zz_x_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_x_xy[i] = 4.0 * g_xy_zz_x_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_x_xz[i] = 4.0 * g_xy_zz_x_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_x_yy[i] = 4.0 * g_xy_zz_x_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_x_yz[i] = 4.0 * g_xy_zz_x_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_x_zz[i] = 4.0 * g_xy_zz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_xy_0_0_0_0_zz_y_xx, g_xy_0_0_0_0_zz_y_xy, g_xy_0_0_0_0_zz_y_xz, g_xy_0_0_0_0_zz_y_yy, g_xy_0_0_0_0_zz_y_yz, g_xy_0_0_0_0_zz_y_zz, g_xy_zz_y_xx, g_xy_zz_y_xy, g_xy_zz_y_xz, g_xy_zz_y_yy, g_xy_zz_y_yz, g_xy_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_zz_y_xx[i] = 4.0 * g_xy_zz_y_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_y_xy[i] = 4.0 * g_xy_zz_y_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_y_xz[i] = 4.0 * g_xy_zz_y_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_y_yy[i] = 4.0 * g_xy_zz_y_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_y_yz[i] = 4.0 * g_xy_zz_y_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_y_zz[i] = 4.0 * g_xy_zz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_xy_0_0_0_0_zz_z_xx, g_xy_0_0_0_0_zz_z_xy, g_xy_0_0_0_0_zz_z_xz, g_xy_0_0_0_0_zz_z_yy, g_xy_0_0_0_0_zz_z_yz, g_xy_0_0_0_0_zz_z_zz, g_xy_zz_z_xx, g_xy_zz_z_xy, g_xy_zz_z_xz, g_xy_zz_z_yy, g_xy_zz_z_yz, g_xy_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_zz_z_xx[i] = 4.0 * g_xy_zz_z_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_z_xy[i] = 4.0 * g_xy_zz_z_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_z_xz[i] = 4.0 * g_xy_zz_z_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_z_yy[i] = 4.0 * g_xy_zz_z_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_z_yz[i] = 4.0 * g_xy_zz_z_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_z_zz[i] = 4.0 * g_xy_zz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_xz_0_0_0_0_xx_x_xx, g_xz_0_0_0_0_xx_x_xy, g_xz_0_0_0_0_xx_x_xz, g_xz_0_0_0_0_xx_x_yy, g_xz_0_0_0_0_xx_x_yz, g_xz_0_0_0_0_xx_x_zz, g_xz_xx_x_xx, g_xz_xx_x_xy, g_xz_xx_x_xz, g_xz_xx_x_yy, g_xz_xx_x_yz, g_xz_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_xx_x_xx[i] = 4.0 * g_xz_xx_x_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_x_xy[i] = 4.0 * g_xz_xx_x_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_x_xz[i] = 4.0 * g_xz_xx_x_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_x_yy[i] = 4.0 * g_xz_xx_x_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_x_yz[i] = 4.0 * g_xz_xx_x_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_x_zz[i] = 4.0 * g_xz_xx_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_xz_0_0_0_0_xx_y_xx, g_xz_0_0_0_0_xx_y_xy, g_xz_0_0_0_0_xx_y_xz, g_xz_0_0_0_0_xx_y_yy, g_xz_0_0_0_0_xx_y_yz, g_xz_0_0_0_0_xx_y_zz, g_xz_xx_y_xx, g_xz_xx_y_xy, g_xz_xx_y_xz, g_xz_xx_y_yy, g_xz_xx_y_yz, g_xz_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_xx_y_xx[i] = 4.0 * g_xz_xx_y_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_y_xy[i] = 4.0 * g_xz_xx_y_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_y_xz[i] = 4.0 * g_xz_xx_y_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_y_yy[i] = 4.0 * g_xz_xx_y_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_y_yz[i] = 4.0 * g_xz_xx_y_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_y_zz[i] = 4.0 * g_xz_xx_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_xz_0_0_0_0_xx_z_xx, g_xz_0_0_0_0_xx_z_xy, g_xz_0_0_0_0_xx_z_xz, g_xz_0_0_0_0_xx_z_yy, g_xz_0_0_0_0_xx_z_yz, g_xz_0_0_0_0_xx_z_zz, g_xz_xx_z_xx, g_xz_xx_z_xy, g_xz_xx_z_xz, g_xz_xx_z_yy, g_xz_xx_z_yz, g_xz_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_xx_z_xx[i] = 4.0 * g_xz_xx_z_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_z_xy[i] = 4.0 * g_xz_xx_z_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_z_xz[i] = 4.0 * g_xz_xx_z_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_z_yy[i] = 4.0 * g_xz_xx_z_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_z_yz[i] = 4.0 * g_xz_xx_z_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_z_zz[i] = 4.0 * g_xz_xx_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_xz_0_0_0_0_xy_x_xx, g_xz_0_0_0_0_xy_x_xy, g_xz_0_0_0_0_xy_x_xz, g_xz_0_0_0_0_xy_x_yy, g_xz_0_0_0_0_xy_x_yz, g_xz_0_0_0_0_xy_x_zz, g_xz_xy_x_xx, g_xz_xy_x_xy, g_xz_xy_x_xz, g_xz_xy_x_yy, g_xz_xy_x_yz, g_xz_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_xy_x_xx[i] = 4.0 * g_xz_xy_x_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_x_xy[i] = 4.0 * g_xz_xy_x_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_x_xz[i] = 4.0 * g_xz_xy_x_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_x_yy[i] = 4.0 * g_xz_xy_x_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_x_yz[i] = 4.0 * g_xz_xy_x_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_x_zz[i] = 4.0 * g_xz_xy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_xz_0_0_0_0_xy_y_xx, g_xz_0_0_0_0_xy_y_xy, g_xz_0_0_0_0_xy_y_xz, g_xz_0_0_0_0_xy_y_yy, g_xz_0_0_0_0_xy_y_yz, g_xz_0_0_0_0_xy_y_zz, g_xz_xy_y_xx, g_xz_xy_y_xy, g_xz_xy_y_xz, g_xz_xy_y_yy, g_xz_xy_y_yz, g_xz_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_xy_y_xx[i] = 4.0 * g_xz_xy_y_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_y_xy[i] = 4.0 * g_xz_xy_y_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_y_xz[i] = 4.0 * g_xz_xy_y_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_y_yy[i] = 4.0 * g_xz_xy_y_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_y_yz[i] = 4.0 * g_xz_xy_y_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_y_zz[i] = 4.0 * g_xz_xy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_xz_0_0_0_0_xy_z_xx, g_xz_0_0_0_0_xy_z_xy, g_xz_0_0_0_0_xy_z_xz, g_xz_0_0_0_0_xy_z_yy, g_xz_0_0_0_0_xy_z_yz, g_xz_0_0_0_0_xy_z_zz, g_xz_xy_z_xx, g_xz_xy_z_xy, g_xz_xy_z_xz, g_xz_xy_z_yy, g_xz_xy_z_yz, g_xz_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_xy_z_xx[i] = 4.0 * g_xz_xy_z_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_z_xy[i] = 4.0 * g_xz_xy_z_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_z_xz[i] = 4.0 * g_xz_xy_z_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_z_yy[i] = 4.0 * g_xz_xy_z_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_z_yz[i] = 4.0 * g_xz_xy_z_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_z_zz[i] = 4.0 * g_xz_xy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_xz_0_0_0_0_xz_x_xx, g_xz_0_0_0_0_xz_x_xy, g_xz_0_0_0_0_xz_x_xz, g_xz_0_0_0_0_xz_x_yy, g_xz_0_0_0_0_xz_x_yz, g_xz_0_0_0_0_xz_x_zz, g_xz_xz_x_xx, g_xz_xz_x_xy, g_xz_xz_x_xz, g_xz_xz_x_yy, g_xz_xz_x_yz, g_xz_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_xz_x_xx[i] = 4.0 * g_xz_xz_x_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_x_xy[i] = 4.0 * g_xz_xz_x_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_x_xz[i] = 4.0 * g_xz_xz_x_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_x_yy[i] = 4.0 * g_xz_xz_x_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_x_yz[i] = 4.0 * g_xz_xz_x_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_x_zz[i] = 4.0 * g_xz_xz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_xz_0_0_0_0_xz_y_xx, g_xz_0_0_0_0_xz_y_xy, g_xz_0_0_0_0_xz_y_xz, g_xz_0_0_0_0_xz_y_yy, g_xz_0_0_0_0_xz_y_yz, g_xz_0_0_0_0_xz_y_zz, g_xz_xz_y_xx, g_xz_xz_y_xy, g_xz_xz_y_xz, g_xz_xz_y_yy, g_xz_xz_y_yz, g_xz_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_xz_y_xx[i] = 4.0 * g_xz_xz_y_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_y_xy[i] = 4.0 * g_xz_xz_y_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_y_xz[i] = 4.0 * g_xz_xz_y_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_y_yy[i] = 4.0 * g_xz_xz_y_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_y_yz[i] = 4.0 * g_xz_xz_y_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_y_zz[i] = 4.0 * g_xz_xz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_xz_0_0_0_0_xz_z_xx, g_xz_0_0_0_0_xz_z_xy, g_xz_0_0_0_0_xz_z_xz, g_xz_0_0_0_0_xz_z_yy, g_xz_0_0_0_0_xz_z_yz, g_xz_0_0_0_0_xz_z_zz, g_xz_xz_z_xx, g_xz_xz_z_xy, g_xz_xz_z_xz, g_xz_xz_z_yy, g_xz_xz_z_yz, g_xz_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_xz_z_xx[i] = 4.0 * g_xz_xz_z_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_z_xy[i] = 4.0 * g_xz_xz_z_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_z_xz[i] = 4.0 * g_xz_xz_z_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_z_yy[i] = 4.0 * g_xz_xz_z_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_z_yz[i] = 4.0 * g_xz_xz_z_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_z_zz[i] = 4.0 * g_xz_xz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_xz_0_0_0_0_yy_x_xx, g_xz_0_0_0_0_yy_x_xy, g_xz_0_0_0_0_yy_x_xz, g_xz_0_0_0_0_yy_x_yy, g_xz_0_0_0_0_yy_x_yz, g_xz_0_0_0_0_yy_x_zz, g_xz_yy_x_xx, g_xz_yy_x_xy, g_xz_yy_x_xz, g_xz_yy_x_yy, g_xz_yy_x_yz, g_xz_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_yy_x_xx[i] = 4.0 * g_xz_yy_x_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_x_xy[i] = 4.0 * g_xz_yy_x_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_x_xz[i] = 4.0 * g_xz_yy_x_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_x_yy[i] = 4.0 * g_xz_yy_x_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_x_yz[i] = 4.0 * g_xz_yy_x_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_x_zz[i] = 4.0 * g_xz_yy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_xz_0_0_0_0_yy_y_xx, g_xz_0_0_0_0_yy_y_xy, g_xz_0_0_0_0_yy_y_xz, g_xz_0_0_0_0_yy_y_yy, g_xz_0_0_0_0_yy_y_yz, g_xz_0_0_0_0_yy_y_zz, g_xz_yy_y_xx, g_xz_yy_y_xy, g_xz_yy_y_xz, g_xz_yy_y_yy, g_xz_yy_y_yz, g_xz_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_yy_y_xx[i] = 4.0 * g_xz_yy_y_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_y_xy[i] = 4.0 * g_xz_yy_y_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_y_xz[i] = 4.0 * g_xz_yy_y_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_y_yy[i] = 4.0 * g_xz_yy_y_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_y_yz[i] = 4.0 * g_xz_yy_y_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_y_zz[i] = 4.0 * g_xz_yy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_xz_0_0_0_0_yy_z_xx, g_xz_0_0_0_0_yy_z_xy, g_xz_0_0_0_0_yy_z_xz, g_xz_0_0_0_0_yy_z_yy, g_xz_0_0_0_0_yy_z_yz, g_xz_0_0_0_0_yy_z_zz, g_xz_yy_z_xx, g_xz_yy_z_xy, g_xz_yy_z_xz, g_xz_yy_z_yy, g_xz_yy_z_yz, g_xz_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_yy_z_xx[i] = 4.0 * g_xz_yy_z_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_z_xy[i] = 4.0 * g_xz_yy_z_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_z_xz[i] = 4.0 * g_xz_yy_z_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_z_yy[i] = 4.0 * g_xz_yy_z_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_z_yz[i] = 4.0 * g_xz_yy_z_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_z_zz[i] = 4.0 * g_xz_yy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_xz_0_0_0_0_yz_x_xx, g_xz_0_0_0_0_yz_x_xy, g_xz_0_0_0_0_yz_x_xz, g_xz_0_0_0_0_yz_x_yy, g_xz_0_0_0_0_yz_x_yz, g_xz_0_0_0_0_yz_x_zz, g_xz_yz_x_xx, g_xz_yz_x_xy, g_xz_yz_x_xz, g_xz_yz_x_yy, g_xz_yz_x_yz, g_xz_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_yz_x_xx[i] = 4.0 * g_xz_yz_x_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_x_xy[i] = 4.0 * g_xz_yz_x_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_x_xz[i] = 4.0 * g_xz_yz_x_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_x_yy[i] = 4.0 * g_xz_yz_x_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_x_yz[i] = 4.0 * g_xz_yz_x_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_x_zz[i] = 4.0 * g_xz_yz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_xz_0_0_0_0_yz_y_xx, g_xz_0_0_0_0_yz_y_xy, g_xz_0_0_0_0_yz_y_xz, g_xz_0_0_0_0_yz_y_yy, g_xz_0_0_0_0_yz_y_yz, g_xz_0_0_0_0_yz_y_zz, g_xz_yz_y_xx, g_xz_yz_y_xy, g_xz_yz_y_xz, g_xz_yz_y_yy, g_xz_yz_y_yz, g_xz_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_yz_y_xx[i] = 4.0 * g_xz_yz_y_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_y_xy[i] = 4.0 * g_xz_yz_y_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_y_xz[i] = 4.0 * g_xz_yz_y_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_y_yy[i] = 4.0 * g_xz_yz_y_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_y_yz[i] = 4.0 * g_xz_yz_y_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_y_zz[i] = 4.0 * g_xz_yz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_xz_0_0_0_0_yz_z_xx, g_xz_0_0_0_0_yz_z_xy, g_xz_0_0_0_0_yz_z_xz, g_xz_0_0_0_0_yz_z_yy, g_xz_0_0_0_0_yz_z_yz, g_xz_0_0_0_0_yz_z_zz, g_xz_yz_z_xx, g_xz_yz_z_xy, g_xz_yz_z_xz, g_xz_yz_z_yy, g_xz_yz_z_yz, g_xz_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_yz_z_xx[i] = 4.0 * g_xz_yz_z_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_z_xy[i] = 4.0 * g_xz_yz_z_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_z_xz[i] = 4.0 * g_xz_yz_z_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_z_yy[i] = 4.0 * g_xz_yz_z_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_z_yz[i] = 4.0 * g_xz_yz_z_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_z_zz[i] = 4.0 * g_xz_yz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_xz_0_0_0_0_zz_x_xx, g_xz_0_0_0_0_zz_x_xy, g_xz_0_0_0_0_zz_x_xz, g_xz_0_0_0_0_zz_x_yy, g_xz_0_0_0_0_zz_x_yz, g_xz_0_0_0_0_zz_x_zz, g_xz_zz_x_xx, g_xz_zz_x_xy, g_xz_zz_x_xz, g_xz_zz_x_yy, g_xz_zz_x_yz, g_xz_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_zz_x_xx[i] = 4.0 * g_xz_zz_x_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_x_xy[i] = 4.0 * g_xz_zz_x_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_x_xz[i] = 4.0 * g_xz_zz_x_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_x_yy[i] = 4.0 * g_xz_zz_x_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_x_yz[i] = 4.0 * g_xz_zz_x_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_x_zz[i] = 4.0 * g_xz_zz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_xz_0_0_0_0_zz_y_xx, g_xz_0_0_0_0_zz_y_xy, g_xz_0_0_0_0_zz_y_xz, g_xz_0_0_0_0_zz_y_yy, g_xz_0_0_0_0_zz_y_yz, g_xz_0_0_0_0_zz_y_zz, g_xz_zz_y_xx, g_xz_zz_y_xy, g_xz_zz_y_xz, g_xz_zz_y_yy, g_xz_zz_y_yz, g_xz_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_zz_y_xx[i] = 4.0 * g_xz_zz_y_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_y_xy[i] = 4.0 * g_xz_zz_y_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_y_xz[i] = 4.0 * g_xz_zz_y_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_y_yy[i] = 4.0 * g_xz_zz_y_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_y_yz[i] = 4.0 * g_xz_zz_y_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_y_zz[i] = 4.0 * g_xz_zz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_xz_0_0_0_0_zz_z_xx, g_xz_0_0_0_0_zz_z_xy, g_xz_0_0_0_0_zz_z_xz, g_xz_0_0_0_0_zz_z_yy, g_xz_0_0_0_0_zz_z_yz, g_xz_0_0_0_0_zz_z_zz, g_xz_zz_z_xx, g_xz_zz_z_xy, g_xz_zz_z_xz, g_xz_zz_z_yy, g_xz_zz_z_yz, g_xz_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_zz_z_xx[i] = 4.0 * g_xz_zz_z_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_z_xy[i] = 4.0 * g_xz_zz_z_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_z_xz[i] = 4.0 * g_xz_zz_z_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_z_yy[i] = 4.0 * g_xz_zz_z_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_z_yz[i] = 4.0 * g_xz_zz_z_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_z_zz[i] = 4.0 * g_xz_zz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (324-330)

    #pragma omp simd aligned(g_0_xx_x_xx, g_0_xx_x_xy, g_0_xx_x_xz, g_0_xx_x_yy, g_0_xx_x_yz, g_0_xx_x_zz, g_yy_0_0_0_0_xx_x_xx, g_yy_0_0_0_0_xx_x_xy, g_yy_0_0_0_0_xx_x_xz, g_yy_0_0_0_0_xx_x_yy, g_yy_0_0_0_0_xx_x_yz, g_yy_0_0_0_0_xx_x_zz, g_yy_xx_x_xx, g_yy_xx_x_xy, g_yy_xx_x_xz, g_yy_xx_x_yy, g_yy_xx_x_yz, g_yy_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_xx_x_xx[i] = -2.0 * g_0_xx_x_xx[i] * a_exp + 4.0 * g_yy_xx_x_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_x_xy[i] = -2.0 * g_0_xx_x_xy[i] * a_exp + 4.0 * g_yy_xx_x_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_x_xz[i] = -2.0 * g_0_xx_x_xz[i] * a_exp + 4.0 * g_yy_xx_x_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_x_yy[i] = -2.0 * g_0_xx_x_yy[i] * a_exp + 4.0 * g_yy_xx_x_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_x_yz[i] = -2.0 * g_0_xx_x_yz[i] * a_exp + 4.0 * g_yy_xx_x_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_x_zz[i] = -2.0 * g_0_xx_x_zz[i] * a_exp + 4.0 * g_yy_xx_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (330-336)

    #pragma omp simd aligned(g_0_xx_y_xx, g_0_xx_y_xy, g_0_xx_y_xz, g_0_xx_y_yy, g_0_xx_y_yz, g_0_xx_y_zz, g_yy_0_0_0_0_xx_y_xx, g_yy_0_0_0_0_xx_y_xy, g_yy_0_0_0_0_xx_y_xz, g_yy_0_0_0_0_xx_y_yy, g_yy_0_0_0_0_xx_y_yz, g_yy_0_0_0_0_xx_y_zz, g_yy_xx_y_xx, g_yy_xx_y_xy, g_yy_xx_y_xz, g_yy_xx_y_yy, g_yy_xx_y_yz, g_yy_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_xx_y_xx[i] = -2.0 * g_0_xx_y_xx[i] * a_exp + 4.0 * g_yy_xx_y_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_y_xy[i] = -2.0 * g_0_xx_y_xy[i] * a_exp + 4.0 * g_yy_xx_y_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_y_xz[i] = -2.0 * g_0_xx_y_xz[i] * a_exp + 4.0 * g_yy_xx_y_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_y_yy[i] = -2.0 * g_0_xx_y_yy[i] * a_exp + 4.0 * g_yy_xx_y_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_y_yz[i] = -2.0 * g_0_xx_y_yz[i] * a_exp + 4.0 * g_yy_xx_y_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_y_zz[i] = -2.0 * g_0_xx_y_zz[i] * a_exp + 4.0 * g_yy_xx_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (336-342)

    #pragma omp simd aligned(g_0_xx_z_xx, g_0_xx_z_xy, g_0_xx_z_xz, g_0_xx_z_yy, g_0_xx_z_yz, g_0_xx_z_zz, g_yy_0_0_0_0_xx_z_xx, g_yy_0_0_0_0_xx_z_xy, g_yy_0_0_0_0_xx_z_xz, g_yy_0_0_0_0_xx_z_yy, g_yy_0_0_0_0_xx_z_yz, g_yy_0_0_0_0_xx_z_zz, g_yy_xx_z_xx, g_yy_xx_z_xy, g_yy_xx_z_xz, g_yy_xx_z_yy, g_yy_xx_z_yz, g_yy_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_xx_z_xx[i] = -2.0 * g_0_xx_z_xx[i] * a_exp + 4.0 * g_yy_xx_z_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_z_xy[i] = -2.0 * g_0_xx_z_xy[i] * a_exp + 4.0 * g_yy_xx_z_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_z_xz[i] = -2.0 * g_0_xx_z_xz[i] * a_exp + 4.0 * g_yy_xx_z_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_z_yy[i] = -2.0 * g_0_xx_z_yy[i] * a_exp + 4.0 * g_yy_xx_z_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_z_yz[i] = -2.0 * g_0_xx_z_yz[i] * a_exp + 4.0 * g_yy_xx_z_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_z_zz[i] = -2.0 * g_0_xx_z_zz[i] * a_exp + 4.0 * g_yy_xx_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (342-348)

    #pragma omp simd aligned(g_0_xy_x_xx, g_0_xy_x_xy, g_0_xy_x_xz, g_0_xy_x_yy, g_0_xy_x_yz, g_0_xy_x_zz, g_yy_0_0_0_0_xy_x_xx, g_yy_0_0_0_0_xy_x_xy, g_yy_0_0_0_0_xy_x_xz, g_yy_0_0_0_0_xy_x_yy, g_yy_0_0_0_0_xy_x_yz, g_yy_0_0_0_0_xy_x_zz, g_yy_xy_x_xx, g_yy_xy_x_xy, g_yy_xy_x_xz, g_yy_xy_x_yy, g_yy_xy_x_yz, g_yy_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_xy_x_xx[i] = -2.0 * g_0_xy_x_xx[i] * a_exp + 4.0 * g_yy_xy_x_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_x_xy[i] = -2.0 * g_0_xy_x_xy[i] * a_exp + 4.0 * g_yy_xy_x_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_x_xz[i] = -2.0 * g_0_xy_x_xz[i] * a_exp + 4.0 * g_yy_xy_x_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_x_yy[i] = -2.0 * g_0_xy_x_yy[i] * a_exp + 4.0 * g_yy_xy_x_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_x_yz[i] = -2.0 * g_0_xy_x_yz[i] * a_exp + 4.0 * g_yy_xy_x_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_x_zz[i] = -2.0 * g_0_xy_x_zz[i] * a_exp + 4.0 * g_yy_xy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (348-354)

    #pragma omp simd aligned(g_0_xy_y_xx, g_0_xy_y_xy, g_0_xy_y_xz, g_0_xy_y_yy, g_0_xy_y_yz, g_0_xy_y_zz, g_yy_0_0_0_0_xy_y_xx, g_yy_0_0_0_0_xy_y_xy, g_yy_0_0_0_0_xy_y_xz, g_yy_0_0_0_0_xy_y_yy, g_yy_0_0_0_0_xy_y_yz, g_yy_0_0_0_0_xy_y_zz, g_yy_xy_y_xx, g_yy_xy_y_xy, g_yy_xy_y_xz, g_yy_xy_y_yy, g_yy_xy_y_yz, g_yy_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_xy_y_xx[i] = -2.0 * g_0_xy_y_xx[i] * a_exp + 4.0 * g_yy_xy_y_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_y_xy[i] = -2.0 * g_0_xy_y_xy[i] * a_exp + 4.0 * g_yy_xy_y_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_y_xz[i] = -2.0 * g_0_xy_y_xz[i] * a_exp + 4.0 * g_yy_xy_y_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_y_yy[i] = -2.0 * g_0_xy_y_yy[i] * a_exp + 4.0 * g_yy_xy_y_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_y_yz[i] = -2.0 * g_0_xy_y_yz[i] * a_exp + 4.0 * g_yy_xy_y_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_y_zz[i] = -2.0 * g_0_xy_y_zz[i] * a_exp + 4.0 * g_yy_xy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (354-360)

    #pragma omp simd aligned(g_0_xy_z_xx, g_0_xy_z_xy, g_0_xy_z_xz, g_0_xy_z_yy, g_0_xy_z_yz, g_0_xy_z_zz, g_yy_0_0_0_0_xy_z_xx, g_yy_0_0_0_0_xy_z_xy, g_yy_0_0_0_0_xy_z_xz, g_yy_0_0_0_0_xy_z_yy, g_yy_0_0_0_0_xy_z_yz, g_yy_0_0_0_0_xy_z_zz, g_yy_xy_z_xx, g_yy_xy_z_xy, g_yy_xy_z_xz, g_yy_xy_z_yy, g_yy_xy_z_yz, g_yy_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_xy_z_xx[i] = -2.0 * g_0_xy_z_xx[i] * a_exp + 4.0 * g_yy_xy_z_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_z_xy[i] = -2.0 * g_0_xy_z_xy[i] * a_exp + 4.0 * g_yy_xy_z_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_z_xz[i] = -2.0 * g_0_xy_z_xz[i] * a_exp + 4.0 * g_yy_xy_z_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_z_yy[i] = -2.0 * g_0_xy_z_yy[i] * a_exp + 4.0 * g_yy_xy_z_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_z_yz[i] = -2.0 * g_0_xy_z_yz[i] * a_exp + 4.0 * g_yy_xy_z_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_z_zz[i] = -2.0 * g_0_xy_z_zz[i] * a_exp + 4.0 * g_yy_xy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (360-366)

    #pragma omp simd aligned(g_0_xz_x_xx, g_0_xz_x_xy, g_0_xz_x_xz, g_0_xz_x_yy, g_0_xz_x_yz, g_0_xz_x_zz, g_yy_0_0_0_0_xz_x_xx, g_yy_0_0_0_0_xz_x_xy, g_yy_0_0_0_0_xz_x_xz, g_yy_0_0_0_0_xz_x_yy, g_yy_0_0_0_0_xz_x_yz, g_yy_0_0_0_0_xz_x_zz, g_yy_xz_x_xx, g_yy_xz_x_xy, g_yy_xz_x_xz, g_yy_xz_x_yy, g_yy_xz_x_yz, g_yy_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_xz_x_xx[i] = -2.0 * g_0_xz_x_xx[i] * a_exp + 4.0 * g_yy_xz_x_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_x_xy[i] = -2.0 * g_0_xz_x_xy[i] * a_exp + 4.0 * g_yy_xz_x_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_x_xz[i] = -2.0 * g_0_xz_x_xz[i] * a_exp + 4.0 * g_yy_xz_x_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_x_yy[i] = -2.0 * g_0_xz_x_yy[i] * a_exp + 4.0 * g_yy_xz_x_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_x_yz[i] = -2.0 * g_0_xz_x_yz[i] * a_exp + 4.0 * g_yy_xz_x_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_x_zz[i] = -2.0 * g_0_xz_x_zz[i] * a_exp + 4.0 * g_yy_xz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (366-372)

    #pragma omp simd aligned(g_0_xz_y_xx, g_0_xz_y_xy, g_0_xz_y_xz, g_0_xz_y_yy, g_0_xz_y_yz, g_0_xz_y_zz, g_yy_0_0_0_0_xz_y_xx, g_yy_0_0_0_0_xz_y_xy, g_yy_0_0_0_0_xz_y_xz, g_yy_0_0_0_0_xz_y_yy, g_yy_0_0_0_0_xz_y_yz, g_yy_0_0_0_0_xz_y_zz, g_yy_xz_y_xx, g_yy_xz_y_xy, g_yy_xz_y_xz, g_yy_xz_y_yy, g_yy_xz_y_yz, g_yy_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_xz_y_xx[i] = -2.0 * g_0_xz_y_xx[i] * a_exp + 4.0 * g_yy_xz_y_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_y_xy[i] = -2.0 * g_0_xz_y_xy[i] * a_exp + 4.0 * g_yy_xz_y_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_y_xz[i] = -2.0 * g_0_xz_y_xz[i] * a_exp + 4.0 * g_yy_xz_y_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_y_yy[i] = -2.0 * g_0_xz_y_yy[i] * a_exp + 4.0 * g_yy_xz_y_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_y_yz[i] = -2.0 * g_0_xz_y_yz[i] * a_exp + 4.0 * g_yy_xz_y_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_y_zz[i] = -2.0 * g_0_xz_y_zz[i] * a_exp + 4.0 * g_yy_xz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (372-378)

    #pragma omp simd aligned(g_0_xz_z_xx, g_0_xz_z_xy, g_0_xz_z_xz, g_0_xz_z_yy, g_0_xz_z_yz, g_0_xz_z_zz, g_yy_0_0_0_0_xz_z_xx, g_yy_0_0_0_0_xz_z_xy, g_yy_0_0_0_0_xz_z_xz, g_yy_0_0_0_0_xz_z_yy, g_yy_0_0_0_0_xz_z_yz, g_yy_0_0_0_0_xz_z_zz, g_yy_xz_z_xx, g_yy_xz_z_xy, g_yy_xz_z_xz, g_yy_xz_z_yy, g_yy_xz_z_yz, g_yy_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_xz_z_xx[i] = -2.0 * g_0_xz_z_xx[i] * a_exp + 4.0 * g_yy_xz_z_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_z_xy[i] = -2.0 * g_0_xz_z_xy[i] * a_exp + 4.0 * g_yy_xz_z_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_z_xz[i] = -2.0 * g_0_xz_z_xz[i] * a_exp + 4.0 * g_yy_xz_z_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_z_yy[i] = -2.0 * g_0_xz_z_yy[i] * a_exp + 4.0 * g_yy_xz_z_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_z_yz[i] = -2.0 * g_0_xz_z_yz[i] * a_exp + 4.0 * g_yy_xz_z_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_z_zz[i] = -2.0 * g_0_xz_z_zz[i] * a_exp + 4.0 * g_yy_xz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (378-384)

    #pragma omp simd aligned(g_0_yy_x_xx, g_0_yy_x_xy, g_0_yy_x_xz, g_0_yy_x_yy, g_0_yy_x_yz, g_0_yy_x_zz, g_yy_0_0_0_0_yy_x_xx, g_yy_0_0_0_0_yy_x_xy, g_yy_0_0_0_0_yy_x_xz, g_yy_0_0_0_0_yy_x_yy, g_yy_0_0_0_0_yy_x_yz, g_yy_0_0_0_0_yy_x_zz, g_yy_yy_x_xx, g_yy_yy_x_xy, g_yy_yy_x_xz, g_yy_yy_x_yy, g_yy_yy_x_yz, g_yy_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_yy_x_xx[i] = -2.0 * g_0_yy_x_xx[i] * a_exp + 4.0 * g_yy_yy_x_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_x_xy[i] = -2.0 * g_0_yy_x_xy[i] * a_exp + 4.0 * g_yy_yy_x_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_x_xz[i] = -2.0 * g_0_yy_x_xz[i] * a_exp + 4.0 * g_yy_yy_x_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_x_yy[i] = -2.0 * g_0_yy_x_yy[i] * a_exp + 4.0 * g_yy_yy_x_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_x_yz[i] = -2.0 * g_0_yy_x_yz[i] * a_exp + 4.0 * g_yy_yy_x_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_x_zz[i] = -2.0 * g_0_yy_x_zz[i] * a_exp + 4.0 * g_yy_yy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (384-390)

    #pragma omp simd aligned(g_0_yy_y_xx, g_0_yy_y_xy, g_0_yy_y_xz, g_0_yy_y_yy, g_0_yy_y_yz, g_0_yy_y_zz, g_yy_0_0_0_0_yy_y_xx, g_yy_0_0_0_0_yy_y_xy, g_yy_0_0_0_0_yy_y_xz, g_yy_0_0_0_0_yy_y_yy, g_yy_0_0_0_0_yy_y_yz, g_yy_0_0_0_0_yy_y_zz, g_yy_yy_y_xx, g_yy_yy_y_xy, g_yy_yy_y_xz, g_yy_yy_y_yy, g_yy_yy_y_yz, g_yy_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_yy_y_xx[i] = -2.0 * g_0_yy_y_xx[i] * a_exp + 4.0 * g_yy_yy_y_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_y_xy[i] = -2.0 * g_0_yy_y_xy[i] * a_exp + 4.0 * g_yy_yy_y_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_y_xz[i] = -2.0 * g_0_yy_y_xz[i] * a_exp + 4.0 * g_yy_yy_y_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_y_yy[i] = -2.0 * g_0_yy_y_yy[i] * a_exp + 4.0 * g_yy_yy_y_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_y_yz[i] = -2.0 * g_0_yy_y_yz[i] * a_exp + 4.0 * g_yy_yy_y_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_y_zz[i] = -2.0 * g_0_yy_y_zz[i] * a_exp + 4.0 * g_yy_yy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (390-396)

    #pragma omp simd aligned(g_0_yy_z_xx, g_0_yy_z_xy, g_0_yy_z_xz, g_0_yy_z_yy, g_0_yy_z_yz, g_0_yy_z_zz, g_yy_0_0_0_0_yy_z_xx, g_yy_0_0_0_0_yy_z_xy, g_yy_0_0_0_0_yy_z_xz, g_yy_0_0_0_0_yy_z_yy, g_yy_0_0_0_0_yy_z_yz, g_yy_0_0_0_0_yy_z_zz, g_yy_yy_z_xx, g_yy_yy_z_xy, g_yy_yy_z_xz, g_yy_yy_z_yy, g_yy_yy_z_yz, g_yy_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_yy_z_xx[i] = -2.0 * g_0_yy_z_xx[i] * a_exp + 4.0 * g_yy_yy_z_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_z_xy[i] = -2.0 * g_0_yy_z_xy[i] * a_exp + 4.0 * g_yy_yy_z_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_z_xz[i] = -2.0 * g_0_yy_z_xz[i] * a_exp + 4.0 * g_yy_yy_z_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_z_yy[i] = -2.0 * g_0_yy_z_yy[i] * a_exp + 4.0 * g_yy_yy_z_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_z_yz[i] = -2.0 * g_0_yy_z_yz[i] * a_exp + 4.0 * g_yy_yy_z_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_z_zz[i] = -2.0 * g_0_yy_z_zz[i] * a_exp + 4.0 * g_yy_yy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (396-402)

    #pragma omp simd aligned(g_0_yz_x_xx, g_0_yz_x_xy, g_0_yz_x_xz, g_0_yz_x_yy, g_0_yz_x_yz, g_0_yz_x_zz, g_yy_0_0_0_0_yz_x_xx, g_yy_0_0_0_0_yz_x_xy, g_yy_0_0_0_0_yz_x_xz, g_yy_0_0_0_0_yz_x_yy, g_yy_0_0_0_0_yz_x_yz, g_yy_0_0_0_0_yz_x_zz, g_yy_yz_x_xx, g_yy_yz_x_xy, g_yy_yz_x_xz, g_yy_yz_x_yy, g_yy_yz_x_yz, g_yy_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_yz_x_xx[i] = -2.0 * g_0_yz_x_xx[i] * a_exp + 4.0 * g_yy_yz_x_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_x_xy[i] = -2.0 * g_0_yz_x_xy[i] * a_exp + 4.0 * g_yy_yz_x_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_x_xz[i] = -2.0 * g_0_yz_x_xz[i] * a_exp + 4.0 * g_yy_yz_x_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_x_yy[i] = -2.0 * g_0_yz_x_yy[i] * a_exp + 4.0 * g_yy_yz_x_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_x_yz[i] = -2.0 * g_0_yz_x_yz[i] * a_exp + 4.0 * g_yy_yz_x_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_x_zz[i] = -2.0 * g_0_yz_x_zz[i] * a_exp + 4.0 * g_yy_yz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (402-408)

    #pragma omp simd aligned(g_0_yz_y_xx, g_0_yz_y_xy, g_0_yz_y_xz, g_0_yz_y_yy, g_0_yz_y_yz, g_0_yz_y_zz, g_yy_0_0_0_0_yz_y_xx, g_yy_0_0_0_0_yz_y_xy, g_yy_0_0_0_0_yz_y_xz, g_yy_0_0_0_0_yz_y_yy, g_yy_0_0_0_0_yz_y_yz, g_yy_0_0_0_0_yz_y_zz, g_yy_yz_y_xx, g_yy_yz_y_xy, g_yy_yz_y_xz, g_yy_yz_y_yy, g_yy_yz_y_yz, g_yy_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_yz_y_xx[i] = -2.0 * g_0_yz_y_xx[i] * a_exp + 4.0 * g_yy_yz_y_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_y_xy[i] = -2.0 * g_0_yz_y_xy[i] * a_exp + 4.0 * g_yy_yz_y_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_y_xz[i] = -2.0 * g_0_yz_y_xz[i] * a_exp + 4.0 * g_yy_yz_y_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_y_yy[i] = -2.0 * g_0_yz_y_yy[i] * a_exp + 4.0 * g_yy_yz_y_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_y_yz[i] = -2.0 * g_0_yz_y_yz[i] * a_exp + 4.0 * g_yy_yz_y_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_y_zz[i] = -2.0 * g_0_yz_y_zz[i] * a_exp + 4.0 * g_yy_yz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (408-414)

    #pragma omp simd aligned(g_0_yz_z_xx, g_0_yz_z_xy, g_0_yz_z_xz, g_0_yz_z_yy, g_0_yz_z_yz, g_0_yz_z_zz, g_yy_0_0_0_0_yz_z_xx, g_yy_0_0_0_0_yz_z_xy, g_yy_0_0_0_0_yz_z_xz, g_yy_0_0_0_0_yz_z_yy, g_yy_0_0_0_0_yz_z_yz, g_yy_0_0_0_0_yz_z_zz, g_yy_yz_z_xx, g_yy_yz_z_xy, g_yy_yz_z_xz, g_yy_yz_z_yy, g_yy_yz_z_yz, g_yy_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_yz_z_xx[i] = -2.0 * g_0_yz_z_xx[i] * a_exp + 4.0 * g_yy_yz_z_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_z_xy[i] = -2.0 * g_0_yz_z_xy[i] * a_exp + 4.0 * g_yy_yz_z_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_z_xz[i] = -2.0 * g_0_yz_z_xz[i] * a_exp + 4.0 * g_yy_yz_z_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_z_yy[i] = -2.0 * g_0_yz_z_yy[i] * a_exp + 4.0 * g_yy_yz_z_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_z_yz[i] = -2.0 * g_0_yz_z_yz[i] * a_exp + 4.0 * g_yy_yz_z_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_z_zz[i] = -2.0 * g_0_yz_z_zz[i] * a_exp + 4.0 * g_yy_yz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (414-420)

    #pragma omp simd aligned(g_0_zz_x_xx, g_0_zz_x_xy, g_0_zz_x_xz, g_0_zz_x_yy, g_0_zz_x_yz, g_0_zz_x_zz, g_yy_0_0_0_0_zz_x_xx, g_yy_0_0_0_0_zz_x_xy, g_yy_0_0_0_0_zz_x_xz, g_yy_0_0_0_0_zz_x_yy, g_yy_0_0_0_0_zz_x_yz, g_yy_0_0_0_0_zz_x_zz, g_yy_zz_x_xx, g_yy_zz_x_xy, g_yy_zz_x_xz, g_yy_zz_x_yy, g_yy_zz_x_yz, g_yy_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_zz_x_xx[i] = -2.0 * g_0_zz_x_xx[i] * a_exp + 4.0 * g_yy_zz_x_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_x_xy[i] = -2.0 * g_0_zz_x_xy[i] * a_exp + 4.0 * g_yy_zz_x_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_x_xz[i] = -2.0 * g_0_zz_x_xz[i] * a_exp + 4.0 * g_yy_zz_x_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_x_yy[i] = -2.0 * g_0_zz_x_yy[i] * a_exp + 4.0 * g_yy_zz_x_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_x_yz[i] = -2.0 * g_0_zz_x_yz[i] * a_exp + 4.0 * g_yy_zz_x_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_x_zz[i] = -2.0 * g_0_zz_x_zz[i] * a_exp + 4.0 * g_yy_zz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (420-426)

    #pragma omp simd aligned(g_0_zz_y_xx, g_0_zz_y_xy, g_0_zz_y_xz, g_0_zz_y_yy, g_0_zz_y_yz, g_0_zz_y_zz, g_yy_0_0_0_0_zz_y_xx, g_yy_0_0_0_0_zz_y_xy, g_yy_0_0_0_0_zz_y_xz, g_yy_0_0_0_0_zz_y_yy, g_yy_0_0_0_0_zz_y_yz, g_yy_0_0_0_0_zz_y_zz, g_yy_zz_y_xx, g_yy_zz_y_xy, g_yy_zz_y_xz, g_yy_zz_y_yy, g_yy_zz_y_yz, g_yy_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_zz_y_xx[i] = -2.0 * g_0_zz_y_xx[i] * a_exp + 4.0 * g_yy_zz_y_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_y_xy[i] = -2.0 * g_0_zz_y_xy[i] * a_exp + 4.0 * g_yy_zz_y_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_y_xz[i] = -2.0 * g_0_zz_y_xz[i] * a_exp + 4.0 * g_yy_zz_y_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_y_yy[i] = -2.0 * g_0_zz_y_yy[i] * a_exp + 4.0 * g_yy_zz_y_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_y_yz[i] = -2.0 * g_0_zz_y_yz[i] * a_exp + 4.0 * g_yy_zz_y_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_y_zz[i] = -2.0 * g_0_zz_y_zz[i] * a_exp + 4.0 * g_yy_zz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (426-432)

    #pragma omp simd aligned(g_0_zz_z_xx, g_0_zz_z_xy, g_0_zz_z_xz, g_0_zz_z_yy, g_0_zz_z_yz, g_0_zz_z_zz, g_yy_0_0_0_0_zz_z_xx, g_yy_0_0_0_0_zz_z_xy, g_yy_0_0_0_0_zz_z_xz, g_yy_0_0_0_0_zz_z_yy, g_yy_0_0_0_0_zz_z_yz, g_yy_0_0_0_0_zz_z_zz, g_yy_zz_z_xx, g_yy_zz_z_xy, g_yy_zz_z_xz, g_yy_zz_z_yy, g_yy_zz_z_yz, g_yy_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_zz_z_xx[i] = -2.0 * g_0_zz_z_xx[i] * a_exp + 4.0 * g_yy_zz_z_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_z_xy[i] = -2.0 * g_0_zz_z_xy[i] * a_exp + 4.0 * g_yy_zz_z_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_z_xz[i] = -2.0 * g_0_zz_z_xz[i] * a_exp + 4.0 * g_yy_zz_z_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_z_yy[i] = -2.0 * g_0_zz_z_yy[i] * a_exp + 4.0 * g_yy_zz_z_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_z_yz[i] = -2.0 * g_0_zz_z_yz[i] * a_exp + 4.0 * g_yy_zz_z_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_z_zz[i] = -2.0 * g_0_zz_z_zz[i] * a_exp + 4.0 * g_yy_zz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (432-438)

    #pragma omp simd aligned(g_yz_0_0_0_0_xx_x_xx, g_yz_0_0_0_0_xx_x_xy, g_yz_0_0_0_0_xx_x_xz, g_yz_0_0_0_0_xx_x_yy, g_yz_0_0_0_0_xx_x_yz, g_yz_0_0_0_0_xx_x_zz, g_yz_xx_x_xx, g_yz_xx_x_xy, g_yz_xx_x_xz, g_yz_xx_x_yy, g_yz_xx_x_yz, g_yz_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_xx_x_xx[i] = 4.0 * g_yz_xx_x_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_x_xy[i] = 4.0 * g_yz_xx_x_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_x_xz[i] = 4.0 * g_yz_xx_x_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_x_yy[i] = 4.0 * g_yz_xx_x_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_x_yz[i] = 4.0 * g_yz_xx_x_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_x_zz[i] = 4.0 * g_yz_xx_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (438-444)

    #pragma omp simd aligned(g_yz_0_0_0_0_xx_y_xx, g_yz_0_0_0_0_xx_y_xy, g_yz_0_0_0_0_xx_y_xz, g_yz_0_0_0_0_xx_y_yy, g_yz_0_0_0_0_xx_y_yz, g_yz_0_0_0_0_xx_y_zz, g_yz_xx_y_xx, g_yz_xx_y_xy, g_yz_xx_y_xz, g_yz_xx_y_yy, g_yz_xx_y_yz, g_yz_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_xx_y_xx[i] = 4.0 * g_yz_xx_y_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_y_xy[i] = 4.0 * g_yz_xx_y_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_y_xz[i] = 4.0 * g_yz_xx_y_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_y_yy[i] = 4.0 * g_yz_xx_y_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_y_yz[i] = 4.0 * g_yz_xx_y_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_y_zz[i] = 4.0 * g_yz_xx_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (444-450)

    #pragma omp simd aligned(g_yz_0_0_0_0_xx_z_xx, g_yz_0_0_0_0_xx_z_xy, g_yz_0_0_0_0_xx_z_xz, g_yz_0_0_0_0_xx_z_yy, g_yz_0_0_0_0_xx_z_yz, g_yz_0_0_0_0_xx_z_zz, g_yz_xx_z_xx, g_yz_xx_z_xy, g_yz_xx_z_xz, g_yz_xx_z_yy, g_yz_xx_z_yz, g_yz_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_xx_z_xx[i] = 4.0 * g_yz_xx_z_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_z_xy[i] = 4.0 * g_yz_xx_z_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_z_xz[i] = 4.0 * g_yz_xx_z_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_z_yy[i] = 4.0 * g_yz_xx_z_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_z_yz[i] = 4.0 * g_yz_xx_z_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_z_zz[i] = 4.0 * g_yz_xx_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (450-456)

    #pragma omp simd aligned(g_yz_0_0_0_0_xy_x_xx, g_yz_0_0_0_0_xy_x_xy, g_yz_0_0_0_0_xy_x_xz, g_yz_0_0_0_0_xy_x_yy, g_yz_0_0_0_0_xy_x_yz, g_yz_0_0_0_0_xy_x_zz, g_yz_xy_x_xx, g_yz_xy_x_xy, g_yz_xy_x_xz, g_yz_xy_x_yy, g_yz_xy_x_yz, g_yz_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_xy_x_xx[i] = 4.0 * g_yz_xy_x_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_x_xy[i] = 4.0 * g_yz_xy_x_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_x_xz[i] = 4.0 * g_yz_xy_x_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_x_yy[i] = 4.0 * g_yz_xy_x_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_x_yz[i] = 4.0 * g_yz_xy_x_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_x_zz[i] = 4.0 * g_yz_xy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (456-462)

    #pragma omp simd aligned(g_yz_0_0_0_0_xy_y_xx, g_yz_0_0_0_0_xy_y_xy, g_yz_0_0_0_0_xy_y_xz, g_yz_0_0_0_0_xy_y_yy, g_yz_0_0_0_0_xy_y_yz, g_yz_0_0_0_0_xy_y_zz, g_yz_xy_y_xx, g_yz_xy_y_xy, g_yz_xy_y_xz, g_yz_xy_y_yy, g_yz_xy_y_yz, g_yz_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_xy_y_xx[i] = 4.0 * g_yz_xy_y_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_y_xy[i] = 4.0 * g_yz_xy_y_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_y_xz[i] = 4.0 * g_yz_xy_y_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_y_yy[i] = 4.0 * g_yz_xy_y_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_y_yz[i] = 4.0 * g_yz_xy_y_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_y_zz[i] = 4.0 * g_yz_xy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (462-468)

    #pragma omp simd aligned(g_yz_0_0_0_0_xy_z_xx, g_yz_0_0_0_0_xy_z_xy, g_yz_0_0_0_0_xy_z_xz, g_yz_0_0_0_0_xy_z_yy, g_yz_0_0_0_0_xy_z_yz, g_yz_0_0_0_0_xy_z_zz, g_yz_xy_z_xx, g_yz_xy_z_xy, g_yz_xy_z_xz, g_yz_xy_z_yy, g_yz_xy_z_yz, g_yz_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_xy_z_xx[i] = 4.0 * g_yz_xy_z_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_z_xy[i] = 4.0 * g_yz_xy_z_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_z_xz[i] = 4.0 * g_yz_xy_z_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_z_yy[i] = 4.0 * g_yz_xy_z_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_z_yz[i] = 4.0 * g_yz_xy_z_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_z_zz[i] = 4.0 * g_yz_xy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (468-474)

    #pragma omp simd aligned(g_yz_0_0_0_0_xz_x_xx, g_yz_0_0_0_0_xz_x_xy, g_yz_0_0_0_0_xz_x_xz, g_yz_0_0_0_0_xz_x_yy, g_yz_0_0_0_0_xz_x_yz, g_yz_0_0_0_0_xz_x_zz, g_yz_xz_x_xx, g_yz_xz_x_xy, g_yz_xz_x_xz, g_yz_xz_x_yy, g_yz_xz_x_yz, g_yz_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_xz_x_xx[i] = 4.0 * g_yz_xz_x_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_x_xy[i] = 4.0 * g_yz_xz_x_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_x_xz[i] = 4.0 * g_yz_xz_x_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_x_yy[i] = 4.0 * g_yz_xz_x_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_x_yz[i] = 4.0 * g_yz_xz_x_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_x_zz[i] = 4.0 * g_yz_xz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (474-480)

    #pragma omp simd aligned(g_yz_0_0_0_0_xz_y_xx, g_yz_0_0_0_0_xz_y_xy, g_yz_0_0_0_0_xz_y_xz, g_yz_0_0_0_0_xz_y_yy, g_yz_0_0_0_0_xz_y_yz, g_yz_0_0_0_0_xz_y_zz, g_yz_xz_y_xx, g_yz_xz_y_xy, g_yz_xz_y_xz, g_yz_xz_y_yy, g_yz_xz_y_yz, g_yz_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_xz_y_xx[i] = 4.0 * g_yz_xz_y_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_y_xy[i] = 4.0 * g_yz_xz_y_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_y_xz[i] = 4.0 * g_yz_xz_y_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_y_yy[i] = 4.0 * g_yz_xz_y_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_y_yz[i] = 4.0 * g_yz_xz_y_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_y_zz[i] = 4.0 * g_yz_xz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (480-486)

    #pragma omp simd aligned(g_yz_0_0_0_0_xz_z_xx, g_yz_0_0_0_0_xz_z_xy, g_yz_0_0_0_0_xz_z_xz, g_yz_0_0_0_0_xz_z_yy, g_yz_0_0_0_0_xz_z_yz, g_yz_0_0_0_0_xz_z_zz, g_yz_xz_z_xx, g_yz_xz_z_xy, g_yz_xz_z_xz, g_yz_xz_z_yy, g_yz_xz_z_yz, g_yz_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_xz_z_xx[i] = 4.0 * g_yz_xz_z_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_z_xy[i] = 4.0 * g_yz_xz_z_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_z_xz[i] = 4.0 * g_yz_xz_z_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_z_yy[i] = 4.0 * g_yz_xz_z_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_z_yz[i] = 4.0 * g_yz_xz_z_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_z_zz[i] = 4.0 * g_yz_xz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (486-492)

    #pragma omp simd aligned(g_yz_0_0_0_0_yy_x_xx, g_yz_0_0_0_0_yy_x_xy, g_yz_0_0_0_0_yy_x_xz, g_yz_0_0_0_0_yy_x_yy, g_yz_0_0_0_0_yy_x_yz, g_yz_0_0_0_0_yy_x_zz, g_yz_yy_x_xx, g_yz_yy_x_xy, g_yz_yy_x_xz, g_yz_yy_x_yy, g_yz_yy_x_yz, g_yz_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_yy_x_xx[i] = 4.0 * g_yz_yy_x_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_x_xy[i] = 4.0 * g_yz_yy_x_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_x_xz[i] = 4.0 * g_yz_yy_x_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_x_yy[i] = 4.0 * g_yz_yy_x_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_x_yz[i] = 4.0 * g_yz_yy_x_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_x_zz[i] = 4.0 * g_yz_yy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (492-498)

    #pragma omp simd aligned(g_yz_0_0_0_0_yy_y_xx, g_yz_0_0_0_0_yy_y_xy, g_yz_0_0_0_0_yy_y_xz, g_yz_0_0_0_0_yy_y_yy, g_yz_0_0_0_0_yy_y_yz, g_yz_0_0_0_0_yy_y_zz, g_yz_yy_y_xx, g_yz_yy_y_xy, g_yz_yy_y_xz, g_yz_yy_y_yy, g_yz_yy_y_yz, g_yz_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_yy_y_xx[i] = 4.0 * g_yz_yy_y_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_y_xy[i] = 4.0 * g_yz_yy_y_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_y_xz[i] = 4.0 * g_yz_yy_y_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_y_yy[i] = 4.0 * g_yz_yy_y_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_y_yz[i] = 4.0 * g_yz_yy_y_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_y_zz[i] = 4.0 * g_yz_yy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (498-504)

    #pragma omp simd aligned(g_yz_0_0_0_0_yy_z_xx, g_yz_0_0_0_0_yy_z_xy, g_yz_0_0_0_0_yy_z_xz, g_yz_0_0_0_0_yy_z_yy, g_yz_0_0_0_0_yy_z_yz, g_yz_0_0_0_0_yy_z_zz, g_yz_yy_z_xx, g_yz_yy_z_xy, g_yz_yy_z_xz, g_yz_yy_z_yy, g_yz_yy_z_yz, g_yz_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_yy_z_xx[i] = 4.0 * g_yz_yy_z_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_z_xy[i] = 4.0 * g_yz_yy_z_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_z_xz[i] = 4.0 * g_yz_yy_z_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_z_yy[i] = 4.0 * g_yz_yy_z_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_z_yz[i] = 4.0 * g_yz_yy_z_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_z_zz[i] = 4.0 * g_yz_yy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (504-510)

    #pragma omp simd aligned(g_yz_0_0_0_0_yz_x_xx, g_yz_0_0_0_0_yz_x_xy, g_yz_0_0_0_0_yz_x_xz, g_yz_0_0_0_0_yz_x_yy, g_yz_0_0_0_0_yz_x_yz, g_yz_0_0_0_0_yz_x_zz, g_yz_yz_x_xx, g_yz_yz_x_xy, g_yz_yz_x_xz, g_yz_yz_x_yy, g_yz_yz_x_yz, g_yz_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_yz_x_xx[i] = 4.0 * g_yz_yz_x_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_x_xy[i] = 4.0 * g_yz_yz_x_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_x_xz[i] = 4.0 * g_yz_yz_x_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_x_yy[i] = 4.0 * g_yz_yz_x_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_x_yz[i] = 4.0 * g_yz_yz_x_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_x_zz[i] = 4.0 * g_yz_yz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (510-516)

    #pragma omp simd aligned(g_yz_0_0_0_0_yz_y_xx, g_yz_0_0_0_0_yz_y_xy, g_yz_0_0_0_0_yz_y_xz, g_yz_0_0_0_0_yz_y_yy, g_yz_0_0_0_0_yz_y_yz, g_yz_0_0_0_0_yz_y_zz, g_yz_yz_y_xx, g_yz_yz_y_xy, g_yz_yz_y_xz, g_yz_yz_y_yy, g_yz_yz_y_yz, g_yz_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_yz_y_xx[i] = 4.0 * g_yz_yz_y_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_y_xy[i] = 4.0 * g_yz_yz_y_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_y_xz[i] = 4.0 * g_yz_yz_y_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_y_yy[i] = 4.0 * g_yz_yz_y_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_y_yz[i] = 4.0 * g_yz_yz_y_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_y_zz[i] = 4.0 * g_yz_yz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (516-522)

    #pragma omp simd aligned(g_yz_0_0_0_0_yz_z_xx, g_yz_0_0_0_0_yz_z_xy, g_yz_0_0_0_0_yz_z_xz, g_yz_0_0_0_0_yz_z_yy, g_yz_0_0_0_0_yz_z_yz, g_yz_0_0_0_0_yz_z_zz, g_yz_yz_z_xx, g_yz_yz_z_xy, g_yz_yz_z_xz, g_yz_yz_z_yy, g_yz_yz_z_yz, g_yz_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_yz_z_xx[i] = 4.0 * g_yz_yz_z_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_z_xy[i] = 4.0 * g_yz_yz_z_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_z_xz[i] = 4.0 * g_yz_yz_z_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_z_yy[i] = 4.0 * g_yz_yz_z_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_z_yz[i] = 4.0 * g_yz_yz_z_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_z_zz[i] = 4.0 * g_yz_yz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (522-528)

    #pragma omp simd aligned(g_yz_0_0_0_0_zz_x_xx, g_yz_0_0_0_0_zz_x_xy, g_yz_0_0_0_0_zz_x_xz, g_yz_0_0_0_0_zz_x_yy, g_yz_0_0_0_0_zz_x_yz, g_yz_0_0_0_0_zz_x_zz, g_yz_zz_x_xx, g_yz_zz_x_xy, g_yz_zz_x_xz, g_yz_zz_x_yy, g_yz_zz_x_yz, g_yz_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_zz_x_xx[i] = 4.0 * g_yz_zz_x_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_x_xy[i] = 4.0 * g_yz_zz_x_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_x_xz[i] = 4.0 * g_yz_zz_x_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_x_yy[i] = 4.0 * g_yz_zz_x_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_x_yz[i] = 4.0 * g_yz_zz_x_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_x_zz[i] = 4.0 * g_yz_zz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (528-534)

    #pragma omp simd aligned(g_yz_0_0_0_0_zz_y_xx, g_yz_0_0_0_0_zz_y_xy, g_yz_0_0_0_0_zz_y_xz, g_yz_0_0_0_0_zz_y_yy, g_yz_0_0_0_0_zz_y_yz, g_yz_0_0_0_0_zz_y_zz, g_yz_zz_y_xx, g_yz_zz_y_xy, g_yz_zz_y_xz, g_yz_zz_y_yy, g_yz_zz_y_yz, g_yz_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_zz_y_xx[i] = 4.0 * g_yz_zz_y_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_y_xy[i] = 4.0 * g_yz_zz_y_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_y_xz[i] = 4.0 * g_yz_zz_y_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_y_yy[i] = 4.0 * g_yz_zz_y_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_y_yz[i] = 4.0 * g_yz_zz_y_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_y_zz[i] = 4.0 * g_yz_zz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (534-540)

    #pragma omp simd aligned(g_yz_0_0_0_0_zz_z_xx, g_yz_0_0_0_0_zz_z_xy, g_yz_0_0_0_0_zz_z_xz, g_yz_0_0_0_0_zz_z_yy, g_yz_0_0_0_0_zz_z_yz, g_yz_0_0_0_0_zz_z_zz, g_yz_zz_z_xx, g_yz_zz_z_xy, g_yz_zz_z_xz, g_yz_zz_z_yy, g_yz_zz_z_yz, g_yz_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_zz_z_xx[i] = 4.0 * g_yz_zz_z_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_z_xy[i] = 4.0 * g_yz_zz_z_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_z_xz[i] = 4.0 * g_yz_zz_z_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_z_yy[i] = 4.0 * g_yz_zz_z_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_z_yz[i] = 4.0 * g_yz_zz_z_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_z_zz[i] = 4.0 * g_yz_zz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (540-546)

    #pragma omp simd aligned(g_0_xx_x_xx, g_0_xx_x_xy, g_0_xx_x_xz, g_0_xx_x_yy, g_0_xx_x_yz, g_0_xx_x_zz, g_zz_0_0_0_0_xx_x_xx, g_zz_0_0_0_0_xx_x_xy, g_zz_0_0_0_0_xx_x_xz, g_zz_0_0_0_0_xx_x_yy, g_zz_0_0_0_0_xx_x_yz, g_zz_0_0_0_0_xx_x_zz, g_zz_xx_x_xx, g_zz_xx_x_xy, g_zz_xx_x_xz, g_zz_xx_x_yy, g_zz_xx_x_yz, g_zz_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_xx_x_xx[i] = -2.0 * g_0_xx_x_xx[i] * a_exp + 4.0 * g_zz_xx_x_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_x_xy[i] = -2.0 * g_0_xx_x_xy[i] * a_exp + 4.0 * g_zz_xx_x_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_x_xz[i] = -2.0 * g_0_xx_x_xz[i] * a_exp + 4.0 * g_zz_xx_x_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_x_yy[i] = -2.0 * g_0_xx_x_yy[i] * a_exp + 4.0 * g_zz_xx_x_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_x_yz[i] = -2.0 * g_0_xx_x_yz[i] * a_exp + 4.0 * g_zz_xx_x_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_x_zz[i] = -2.0 * g_0_xx_x_zz[i] * a_exp + 4.0 * g_zz_xx_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (546-552)

    #pragma omp simd aligned(g_0_xx_y_xx, g_0_xx_y_xy, g_0_xx_y_xz, g_0_xx_y_yy, g_0_xx_y_yz, g_0_xx_y_zz, g_zz_0_0_0_0_xx_y_xx, g_zz_0_0_0_0_xx_y_xy, g_zz_0_0_0_0_xx_y_xz, g_zz_0_0_0_0_xx_y_yy, g_zz_0_0_0_0_xx_y_yz, g_zz_0_0_0_0_xx_y_zz, g_zz_xx_y_xx, g_zz_xx_y_xy, g_zz_xx_y_xz, g_zz_xx_y_yy, g_zz_xx_y_yz, g_zz_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_xx_y_xx[i] = -2.0 * g_0_xx_y_xx[i] * a_exp + 4.0 * g_zz_xx_y_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_y_xy[i] = -2.0 * g_0_xx_y_xy[i] * a_exp + 4.0 * g_zz_xx_y_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_y_xz[i] = -2.0 * g_0_xx_y_xz[i] * a_exp + 4.0 * g_zz_xx_y_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_y_yy[i] = -2.0 * g_0_xx_y_yy[i] * a_exp + 4.0 * g_zz_xx_y_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_y_yz[i] = -2.0 * g_0_xx_y_yz[i] * a_exp + 4.0 * g_zz_xx_y_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_y_zz[i] = -2.0 * g_0_xx_y_zz[i] * a_exp + 4.0 * g_zz_xx_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (552-558)

    #pragma omp simd aligned(g_0_xx_z_xx, g_0_xx_z_xy, g_0_xx_z_xz, g_0_xx_z_yy, g_0_xx_z_yz, g_0_xx_z_zz, g_zz_0_0_0_0_xx_z_xx, g_zz_0_0_0_0_xx_z_xy, g_zz_0_0_0_0_xx_z_xz, g_zz_0_0_0_0_xx_z_yy, g_zz_0_0_0_0_xx_z_yz, g_zz_0_0_0_0_xx_z_zz, g_zz_xx_z_xx, g_zz_xx_z_xy, g_zz_xx_z_xz, g_zz_xx_z_yy, g_zz_xx_z_yz, g_zz_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_xx_z_xx[i] = -2.0 * g_0_xx_z_xx[i] * a_exp + 4.0 * g_zz_xx_z_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_z_xy[i] = -2.0 * g_0_xx_z_xy[i] * a_exp + 4.0 * g_zz_xx_z_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_z_xz[i] = -2.0 * g_0_xx_z_xz[i] * a_exp + 4.0 * g_zz_xx_z_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_z_yy[i] = -2.0 * g_0_xx_z_yy[i] * a_exp + 4.0 * g_zz_xx_z_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_z_yz[i] = -2.0 * g_0_xx_z_yz[i] * a_exp + 4.0 * g_zz_xx_z_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_z_zz[i] = -2.0 * g_0_xx_z_zz[i] * a_exp + 4.0 * g_zz_xx_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (558-564)

    #pragma omp simd aligned(g_0_xy_x_xx, g_0_xy_x_xy, g_0_xy_x_xz, g_0_xy_x_yy, g_0_xy_x_yz, g_0_xy_x_zz, g_zz_0_0_0_0_xy_x_xx, g_zz_0_0_0_0_xy_x_xy, g_zz_0_0_0_0_xy_x_xz, g_zz_0_0_0_0_xy_x_yy, g_zz_0_0_0_0_xy_x_yz, g_zz_0_0_0_0_xy_x_zz, g_zz_xy_x_xx, g_zz_xy_x_xy, g_zz_xy_x_xz, g_zz_xy_x_yy, g_zz_xy_x_yz, g_zz_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_xy_x_xx[i] = -2.0 * g_0_xy_x_xx[i] * a_exp + 4.0 * g_zz_xy_x_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_x_xy[i] = -2.0 * g_0_xy_x_xy[i] * a_exp + 4.0 * g_zz_xy_x_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_x_xz[i] = -2.0 * g_0_xy_x_xz[i] * a_exp + 4.0 * g_zz_xy_x_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_x_yy[i] = -2.0 * g_0_xy_x_yy[i] * a_exp + 4.0 * g_zz_xy_x_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_x_yz[i] = -2.0 * g_0_xy_x_yz[i] * a_exp + 4.0 * g_zz_xy_x_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_x_zz[i] = -2.0 * g_0_xy_x_zz[i] * a_exp + 4.0 * g_zz_xy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (564-570)

    #pragma omp simd aligned(g_0_xy_y_xx, g_0_xy_y_xy, g_0_xy_y_xz, g_0_xy_y_yy, g_0_xy_y_yz, g_0_xy_y_zz, g_zz_0_0_0_0_xy_y_xx, g_zz_0_0_0_0_xy_y_xy, g_zz_0_0_0_0_xy_y_xz, g_zz_0_0_0_0_xy_y_yy, g_zz_0_0_0_0_xy_y_yz, g_zz_0_0_0_0_xy_y_zz, g_zz_xy_y_xx, g_zz_xy_y_xy, g_zz_xy_y_xz, g_zz_xy_y_yy, g_zz_xy_y_yz, g_zz_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_xy_y_xx[i] = -2.0 * g_0_xy_y_xx[i] * a_exp + 4.0 * g_zz_xy_y_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_y_xy[i] = -2.0 * g_0_xy_y_xy[i] * a_exp + 4.0 * g_zz_xy_y_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_y_xz[i] = -2.0 * g_0_xy_y_xz[i] * a_exp + 4.0 * g_zz_xy_y_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_y_yy[i] = -2.0 * g_0_xy_y_yy[i] * a_exp + 4.0 * g_zz_xy_y_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_y_yz[i] = -2.0 * g_0_xy_y_yz[i] * a_exp + 4.0 * g_zz_xy_y_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_y_zz[i] = -2.0 * g_0_xy_y_zz[i] * a_exp + 4.0 * g_zz_xy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (570-576)

    #pragma omp simd aligned(g_0_xy_z_xx, g_0_xy_z_xy, g_0_xy_z_xz, g_0_xy_z_yy, g_0_xy_z_yz, g_0_xy_z_zz, g_zz_0_0_0_0_xy_z_xx, g_zz_0_0_0_0_xy_z_xy, g_zz_0_0_0_0_xy_z_xz, g_zz_0_0_0_0_xy_z_yy, g_zz_0_0_0_0_xy_z_yz, g_zz_0_0_0_0_xy_z_zz, g_zz_xy_z_xx, g_zz_xy_z_xy, g_zz_xy_z_xz, g_zz_xy_z_yy, g_zz_xy_z_yz, g_zz_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_xy_z_xx[i] = -2.0 * g_0_xy_z_xx[i] * a_exp + 4.0 * g_zz_xy_z_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_z_xy[i] = -2.0 * g_0_xy_z_xy[i] * a_exp + 4.0 * g_zz_xy_z_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_z_xz[i] = -2.0 * g_0_xy_z_xz[i] * a_exp + 4.0 * g_zz_xy_z_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_z_yy[i] = -2.0 * g_0_xy_z_yy[i] * a_exp + 4.0 * g_zz_xy_z_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_z_yz[i] = -2.0 * g_0_xy_z_yz[i] * a_exp + 4.0 * g_zz_xy_z_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_z_zz[i] = -2.0 * g_0_xy_z_zz[i] * a_exp + 4.0 * g_zz_xy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (576-582)

    #pragma omp simd aligned(g_0_xz_x_xx, g_0_xz_x_xy, g_0_xz_x_xz, g_0_xz_x_yy, g_0_xz_x_yz, g_0_xz_x_zz, g_zz_0_0_0_0_xz_x_xx, g_zz_0_0_0_0_xz_x_xy, g_zz_0_0_0_0_xz_x_xz, g_zz_0_0_0_0_xz_x_yy, g_zz_0_0_0_0_xz_x_yz, g_zz_0_0_0_0_xz_x_zz, g_zz_xz_x_xx, g_zz_xz_x_xy, g_zz_xz_x_xz, g_zz_xz_x_yy, g_zz_xz_x_yz, g_zz_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_xz_x_xx[i] = -2.0 * g_0_xz_x_xx[i] * a_exp + 4.0 * g_zz_xz_x_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_x_xy[i] = -2.0 * g_0_xz_x_xy[i] * a_exp + 4.0 * g_zz_xz_x_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_x_xz[i] = -2.0 * g_0_xz_x_xz[i] * a_exp + 4.0 * g_zz_xz_x_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_x_yy[i] = -2.0 * g_0_xz_x_yy[i] * a_exp + 4.0 * g_zz_xz_x_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_x_yz[i] = -2.0 * g_0_xz_x_yz[i] * a_exp + 4.0 * g_zz_xz_x_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_x_zz[i] = -2.0 * g_0_xz_x_zz[i] * a_exp + 4.0 * g_zz_xz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (582-588)

    #pragma omp simd aligned(g_0_xz_y_xx, g_0_xz_y_xy, g_0_xz_y_xz, g_0_xz_y_yy, g_0_xz_y_yz, g_0_xz_y_zz, g_zz_0_0_0_0_xz_y_xx, g_zz_0_0_0_0_xz_y_xy, g_zz_0_0_0_0_xz_y_xz, g_zz_0_0_0_0_xz_y_yy, g_zz_0_0_0_0_xz_y_yz, g_zz_0_0_0_0_xz_y_zz, g_zz_xz_y_xx, g_zz_xz_y_xy, g_zz_xz_y_xz, g_zz_xz_y_yy, g_zz_xz_y_yz, g_zz_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_xz_y_xx[i] = -2.0 * g_0_xz_y_xx[i] * a_exp + 4.0 * g_zz_xz_y_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_y_xy[i] = -2.0 * g_0_xz_y_xy[i] * a_exp + 4.0 * g_zz_xz_y_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_y_xz[i] = -2.0 * g_0_xz_y_xz[i] * a_exp + 4.0 * g_zz_xz_y_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_y_yy[i] = -2.0 * g_0_xz_y_yy[i] * a_exp + 4.0 * g_zz_xz_y_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_y_yz[i] = -2.0 * g_0_xz_y_yz[i] * a_exp + 4.0 * g_zz_xz_y_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_y_zz[i] = -2.0 * g_0_xz_y_zz[i] * a_exp + 4.0 * g_zz_xz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (588-594)

    #pragma omp simd aligned(g_0_xz_z_xx, g_0_xz_z_xy, g_0_xz_z_xz, g_0_xz_z_yy, g_0_xz_z_yz, g_0_xz_z_zz, g_zz_0_0_0_0_xz_z_xx, g_zz_0_0_0_0_xz_z_xy, g_zz_0_0_0_0_xz_z_xz, g_zz_0_0_0_0_xz_z_yy, g_zz_0_0_0_0_xz_z_yz, g_zz_0_0_0_0_xz_z_zz, g_zz_xz_z_xx, g_zz_xz_z_xy, g_zz_xz_z_xz, g_zz_xz_z_yy, g_zz_xz_z_yz, g_zz_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_xz_z_xx[i] = -2.0 * g_0_xz_z_xx[i] * a_exp + 4.0 * g_zz_xz_z_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_z_xy[i] = -2.0 * g_0_xz_z_xy[i] * a_exp + 4.0 * g_zz_xz_z_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_z_xz[i] = -2.0 * g_0_xz_z_xz[i] * a_exp + 4.0 * g_zz_xz_z_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_z_yy[i] = -2.0 * g_0_xz_z_yy[i] * a_exp + 4.0 * g_zz_xz_z_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_z_yz[i] = -2.0 * g_0_xz_z_yz[i] * a_exp + 4.0 * g_zz_xz_z_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_z_zz[i] = -2.0 * g_0_xz_z_zz[i] * a_exp + 4.0 * g_zz_xz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (594-600)

    #pragma omp simd aligned(g_0_yy_x_xx, g_0_yy_x_xy, g_0_yy_x_xz, g_0_yy_x_yy, g_0_yy_x_yz, g_0_yy_x_zz, g_zz_0_0_0_0_yy_x_xx, g_zz_0_0_0_0_yy_x_xy, g_zz_0_0_0_0_yy_x_xz, g_zz_0_0_0_0_yy_x_yy, g_zz_0_0_0_0_yy_x_yz, g_zz_0_0_0_0_yy_x_zz, g_zz_yy_x_xx, g_zz_yy_x_xy, g_zz_yy_x_xz, g_zz_yy_x_yy, g_zz_yy_x_yz, g_zz_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_yy_x_xx[i] = -2.0 * g_0_yy_x_xx[i] * a_exp + 4.0 * g_zz_yy_x_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_x_xy[i] = -2.0 * g_0_yy_x_xy[i] * a_exp + 4.0 * g_zz_yy_x_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_x_xz[i] = -2.0 * g_0_yy_x_xz[i] * a_exp + 4.0 * g_zz_yy_x_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_x_yy[i] = -2.0 * g_0_yy_x_yy[i] * a_exp + 4.0 * g_zz_yy_x_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_x_yz[i] = -2.0 * g_0_yy_x_yz[i] * a_exp + 4.0 * g_zz_yy_x_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_x_zz[i] = -2.0 * g_0_yy_x_zz[i] * a_exp + 4.0 * g_zz_yy_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (600-606)

    #pragma omp simd aligned(g_0_yy_y_xx, g_0_yy_y_xy, g_0_yy_y_xz, g_0_yy_y_yy, g_0_yy_y_yz, g_0_yy_y_zz, g_zz_0_0_0_0_yy_y_xx, g_zz_0_0_0_0_yy_y_xy, g_zz_0_0_0_0_yy_y_xz, g_zz_0_0_0_0_yy_y_yy, g_zz_0_0_0_0_yy_y_yz, g_zz_0_0_0_0_yy_y_zz, g_zz_yy_y_xx, g_zz_yy_y_xy, g_zz_yy_y_xz, g_zz_yy_y_yy, g_zz_yy_y_yz, g_zz_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_yy_y_xx[i] = -2.0 * g_0_yy_y_xx[i] * a_exp + 4.0 * g_zz_yy_y_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_y_xy[i] = -2.0 * g_0_yy_y_xy[i] * a_exp + 4.0 * g_zz_yy_y_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_y_xz[i] = -2.0 * g_0_yy_y_xz[i] * a_exp + 4.0 * g_zz_yy_y_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_y_yy[i] = -2.0 * g_0_yy_y_yy[i] * a_exp + 4.0 * g_zz_yy_y_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_y_yz[i] = -2.0 * g_0_yy_y_yz[i] * a_exp + 4.0 * g_zz_yy_y_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_y_zz[i] = -2.0 * g_0_yy_y_zz[i] * a_exp + 4.0 * g_zz_yy_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (606-612)

    #pragma omp simd aligned(g_0_yy_z_xx, g_0_yy_z_xy, g_0_yy_z_xz, g_0_yy_z_yy, g_0_yy_z_yz, g_0_yy_z_zz, g_zz_0_0_0_0_yy_z_xx, g_zz_0_0_0_0_yy_z_xy, g_zz_0_0_0_0_yy_z_xz, g_zz_0_0_0_0_yy_z_yy, g_zz_0_0_0_0_yy_z_yz, g_zz_0_0_0_0_yy_z_zz, g_zz_yy_z_xx, g_zz_yy_z_xy, g_zz_yy_z_xz, g_zz_yy_z_yy, g_zz_yy_z_yz, g_zz_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_yy_z_xx[i] = -2.0 * g_0_yy_z_xx[i] * a_exp + 4.0 * g_zz_yy_z_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_z_xy[i] = -2.0 * g_0_yy_z_xy[i] * a_exp + 4.0 * g_zz_yy_z_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_z_xz[i] = -2.0 * g_0_yy_z_xz[i] * a_exp + 4.0 * g_zz_yy_z_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_z_yy[i] = -2.0 * g_0_yy_z_yy[i] * a_exp + 4.0 * g_zz_yy_z_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_z_yz[i] = -2.0 * g_0_yy_z_yz[i] * a_exp + 4.0 * g_zz_yy_z_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_z_zz[i] = -2.0 * g_0_yy_z_zz[i] * a_exp + 4.0 * g_zz_yy_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (612-618)

    #pragma omp simd aligned(g_0_yz_x_xx, g_0_yz_x_xy, g_0_yz_x_xz, g_0_yz_x_yy, g_0_yz_x_yz, g_0_yz_x_zz, g_zz_0_0_0_0_yz_x_xx, g_zz_0_0_0_0_yz_x_xy, g_zz_0_0_0_0_yz_x_xz, g_zz_0_0_0_0_yz_x_yy, g_zz_0_0_0_0_yz_x_yz, g_zz_0_0_0_0_yz_x_zz, g_zz_yz_x_xx, g_zz_yz_x_xy, g_zz_yz_x_xz, g_zz_yz_x_yy, g_zz_yz_x_yz, g_zz_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_yz_x_xx[i] = -2.0 * g_0_yz_x_xx[i] * a_exp + 4.0 * g_zz_yz_x_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_x_xy[i] = -2.0 * g_0_yz_x_xy[i] * a_exp + 4.0 * g_zz_yz_x_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_x_xz[i] = -2.0 * g_0_yz_x_xz[i] * a_exp + 4.0 * g_zz_yz_x_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_x_yy[i] = -2.0 * g_0_yz_x_yy[i] * a_exp + 4.0 * g_zz_yz_x_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_x_yz[i] = -2.0 * g_0_yz_x_yz[i] * a_exp + 4.0 * g_zz_yz_x_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_x_zz[i] = -2.0 * g_0_yz_x_zz[i] * a_exp + 4.0 * g_zz_yz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (618-624)

    #pragma omp simd aligned(g_0_yz_y_xx, g_0_yz_y_xy, g_0_yz_y_xz, g_0_yz_y_yy, g_0_yz_y_yz, g_0_yz_y_zz, g_zz_0_0_0_0_yz_y_xx, g_zz_0_0_0_0_yz_y_xy, g_zz_0_0_0_0_yz_y_xz, g_zz_0_0_0_0_yz_y_yy, g_zz_0_0_0_0_yz_y_yz, g_zz_0_0_0_0_yz_y_zz, g_zz_yz_y_xx, g_zz_yz_y_xy, g_zz_yz_y_xz, g_zz_yz_y_yy, g_zz_yz_y_yz, g_zz_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_yz_y_xx[i] = -2.0 * g_0_yz_y_xx[i] * a_exp + 4.0 * g_zz_yz_y_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_y_xy[i] = -2.0 * g_0_yz_y_xy[i] * a_exp + 4.0 * g_zz_yz_y_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_y_xz[i] = -2.0 * g_0_yz_y_xz[i] * a_exp + 4.0 * g_zz_yz_y_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_y_yy[i] = -2.0 * g_0_yz_y_yy[i] * a_exp + 4.0 * g_zz_yz_y_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_y_yz[i] = -2.0 * g_0_yz_y_yz[i] * a_exp + 4.0 * g_zz_yz_y_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_y_zz[i] = -2.0 * g_0_yz_y_zz[i] * a_exp + 4.0 * g_zz_yz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (624-630)

    #pragma omp simd aligned(g_0_yz_z_xx, g_0_yz_z_xy, g_0_yz_z_xz, g_0_yz_z_yy, g_0_yz_z_yz, g_0_yz_z_zz, g_zz_0_0_0_0_yz_z_xx, g_zz_0_0_0_0_yz_z_xy, g_zz_0_0_0_0_yz_z_xz, g_zz_0_0_0_0_yz_z_yy, g_zz_0_0_0_0_yz_z_yz, g_zz_0_0_0_0_yz_z_zz, g_zz_yz_z_xx, g_zz_yz_z_xy, g_zz_yz_z_xz, g_zz_yz_z_yy, g_zz_yz_z_yz, g_zz_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_yz_z_xx[i] = -2.0 * g_0_yz_z_xx[i] * a_exp + 4.0 * g_zz_yz_z_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_z_xy[i] = -2.0 * g_0_yz_z_xy[i] * a_exp + 4.0 * g_zz_yz_z_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_z_xz[i] = -2.0 * g_0_yz_z_xz[i] * a_exp + 4.0 * g_zz_yz_z_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_z_yy[i] = -2.0 * g_0_yz_z_yy[i] * a_exp + 4.0 * g_zz_yz_z_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_z_yz[i] = -2.0 * g_0_yz_z_yz[i] * a_exp + 4.0 * g_zz_yz_z_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_z_zz[i] = -2.0 * g_0_yz_z_zz[i] * a_exp + 4.0 * g_zz_yz_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (630-636)

    #pragma omp simd aligned(g_0_zz_x_xx, g_0_zz_x_xy, g_0_zz_x_xz, g_0_zz_x_yy, g_0_zz_x_yz, g_0_zz_x_zz, g_zz_0_0_0_0_zz_x_xx, g_zz_0_0_0_0_zz_x_xy, g_zz_0_0_0_0_zz_x_xz, g_zz_0_0_0_0_zz_x_yy, g_zz_0_0_0_0_zz_x_yz, g_zz_0_0_0_0_zz_x_zz, g_zz_zz_x_xx, g_zz_zz_x_xy, g_zz_zz_x_xz, g_zz_zz_x_yy, g_zz_zz_x_yz, g_zz_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_zz_x_xx[i] = -2.0 * g_0_zz_x_xx[i] * a_exp + 4.0 * g_zz_zz_x_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_x_xy[i] = -2.0 * g_0_zz_x_xy[i] * a_exp + 4.0 * g_zz_zz_x_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_x_xz[i] = -2.0 * g_0_zz_x_xz[i] * a_exp + 4.0 * g_zz_zz_x_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_x_yy[i] = -2.0 * g_0_zz_x_yy[i] * a_exp + 4.0 * g_zz_zz_x_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_x_yz[i] = -2.0 * g_0_zz_x_yz[i] * a_exp + 4.0 * g_zz_zz_x_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_x_zz[i] = -2.0 * g_0_zz_x_zz[i] * a_exp + 4.0 * g_zz_zz_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (636-642)

    #pragma omp simd aligned(g_0_zz_y_xx, g_0_zz_y_xy, g_0_zz_y_xz, g_0_zz_y_yy, g_0_zz_y_yz, g_0_zz_y_zz, g_zz_0_0_0_0_zz_y_xx, g_zz_0_0_0_0_zz_y_xy, g_zz_0_0_0_0_zz_y_xz, g_zz_0_0_0_0_zz_y_yy, g_zz_0_0_0_0_zz_y_yz, g_zz_0_0_0_0_zz_y_zz, g_zz_zz_y_xx, g_zz_zz_y_xy, g_zz_zz_y_xz, g_zz_zz_y_yy, g_zz_zz_y_yz, g_zz_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_zz_y_xx[i] = -2.0 * g_0_zz_y_xx[i] * a_exp + 4.0 * g_zz_zz_y_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_y_xy[i] = -2.0 * g_0_zz_y_xy[i] * a_exp + 4.0 * g_zz_zz_y_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_y_xz[i] = -2.0 * g_0_zz_y_xz[i] * a_exp + 4.0 * g_zz_zz_y_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_y_yy[i] = -2.0 * g_0_zz_y_yy[i] * a_exp + 4.0 * g_zz_zz_y_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_y_yz[i] = -2.0 * g_0_zz_y_yz[i] * a_exp + 4.0 * g_zz_zz_y_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_y_zz[i] = -2.0 * g_0_zz_y_zz[i] * a_exp + 4.0 * g_zz_zz_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (642-648)

    #pragma omp simd aligned(g_0_zz_z_xx, g_0_zz_z_xy, g_0_zz_z_xz, g_0_zz_z_yy, g_0_zz_z_yz, g_0_zz_z_zz, g_zz_0_0_0_0_zz_z_xx, g_zz_0_0_0_0_zz_z_xy, g_zz_0_0_0_0_zz_z_xz, g_zz_0_0_0_0_zz_z_yy, g_zz_0_0_0_0_zz_z_yz, g_zz_0_0_0_0_zz_z_zz, g_zz_zz_z_xx, g_zz_zz_z_xy, g_zz_zz_z_xz, g_zz_zz_z_yy, g_zz_zz_z_yz, g_zz_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_zz_z_xx[i] = -2.0 * g_0_zz_z_xx[i] * a_exp + 4.0 * g_zz_zz_z_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_z_xy[i] = -2.0 * g_0_zz_z_xy[i] * a_exp + 4.0 * g_zz_zz_z_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_z_xz[i] = -2.0 * g_0_zz_z_xz[i] * a_exp + 4.0 * g_zz_zz_z_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_z_yy[i] = -2.0 * g_0_zz_z_yy[i] * a_exp + 4.0 * g_zz_zz_z_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_z_yz[i] = -2.0 * g_0_zz_z_yz[i] * a_exp + 4.0 * g_zz_zz_z_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_z_zz[i] = -2.0 * g_0_zz_z_zz[i] * a_exp + 4.0 * g_zz_zz_z_zz[i] * a_exp * a_exp;
    }
}

} // t4c_geom namespace


#include "GeomDeriv1010OfScalarForPDSD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1010_pdsd_0(CSimdArray<double>& buffer_1010_pdsd,
                     const CSimdArray<double>& buffer_sdpd,
                     const CSimdArray<double>& buffer_ddpd,
                     const double a_exp,
                     const double* c_exps) -> void
{
    const auto ndims = buffer_1010_pdsd.number_of_columns();

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

    /// Set up components of integrals buffer : buffer_1010_pdsd

    auto g_x_0_x_0_x_xx_0_xx = buffer_1010_pdsd[0];

    auto g_x_0_x_0_x_xx_0_xy = buffer_1010_pdsd[1];

    auto g_x_0_x_0_x_xx_0_xz = buffer_1010_pdsd[2];

    auto g_x_0_x_0_x_xx_0_yy = buffer_1010_pdsd[3];

    auto g_x_0_x_0_x_xx_0_yz = buffer_1010_pdsd[4];

    auto g_x_0_x_0_x_xx_0_zz = buffer_1010_pdsd[5];

    auto g_x_0_x_0_x_xy_0_xx = buffer_1010_pdsd[6];

    auto g_x_0_x_0_x_xy_0_xy = buffer_1010_pdsd[7];

    auto g_x_0_x_0_x_xy_0_xz = buffer_1010_pdsd[8];

    auto g_x_0_x_0_x_xy_0_yy = buffer_1010_pdsd[9];

    auto g_x_0_x_0_x_xy_0_yz = buffer_1010_pdsd[10];

    auto g_x_0_x_0_x_xy_0_zz = buffer_1010_pdsd[11];

    auto g_x_0_x_0_x_xz_0_xx = buffer_1010_pdsd[12];

    auto g_x_0_x_0_x_xz_0_xy = buffer_1010_pdsd[13];

    auto g_x_0_x_0_x_xz_0_xz = buffer_1010_pdsd[14];

    auto g_x_0_x_0_x_xz_0_yy = buffer_1010_pdsd[15];

    auto g_x_0_x_0_x_xz_0_yz = buffer_1010_pdsd[16];

    auto g_x_0_x_0_x_xz_0_zz = buffer_1010_pdsd[17];

    auto g_x_0_x_0_x_yy_0_xx = buffer_1010_pdsd[18];

    auto g_x_0_x_0_x_yy_0_xy = buffer_1010_pdsd[19];

    auto g_x_0_x_0_x_yy_0_xz = buffer_1010_pdsd[20];

    auto g_x_0_x_0_x_yy_0_yy = buffer_1010_pdsd[21];

    auto g_x_0_x_0_x_yy_0_yz = buffer_1010_pdsd[22];

    auto g_x_0_x_0_x_yy_0_zz = buffer_1010_pdsd[23];

    auto g_x_0_x_0_x_yz_0_xx = buffer_1010_pdsd[24];

    auto g_x_0_x_0_x_yz_0_xy = buffer_1010_pdsd[25];

    auto g_x_0_x_0_x_yz_0_xz = buffer_1010_pdsd[26];

    auto g_x_0_x_0_x_yz_0_yy = buffer_1010_pdsd[27];

    auto g_x_0_x_0_x_yz_0_yz = buffer_1010_pdsd[28];

    auto g_x_0_x_0_x_yz_0_zz = buffer_1010_pdsd[29];

    auto g_x_0_x_0_x_zz_0_xx = buffer_1010_pdsd[30];

    auto g_x_0_x_0_x_zz_0_xy = buffer_1010_pdsd[31];

    auto g_x_0_x_0_x_zz_0_xz = buffer_1010_pdsd[32];

    auto g_x_0_x_0_x_zz_0_yy = buffer_1010_pdsd[33];

    auto g_x_0_x_0_x_zz_0_yz = buffer_1010_pdsd[34];

    auto g_x_0_x_0_x_zz_0_zz = buffer_1010_pdsd[35];

    auto g_x_0_x_0_y_xx_0_xx = buffer_1010_pdsd[36];

    auto g_x_0_x_0_y_xx_0_xy = buffer_1010_pdsd[37];

    auto g_x_0_x_0_y_xx_0_xz = buffer_1010_pdsd[38];

    auto g_x_0_x_0_y_xx_0_yy = buffer_1010_pdsd[39];

    auto g_x_0_x_0_y_xx_0_yz = buffer_1010_pdsd[40];

    auto g_x_0_x_0_y_xx_0_zz = buffer_1010_pdsd[41];

    auto g_x_0_x_0_y_xy_0_xx = buffer_1010_pdsd[42];

    auto g_x_0_x_0_y_xy_0_xy = buffer_1010_pdsd[43];

    auto g_x_0_x_0_y_xy_0_xz = buffer_1010_pdsd[44];

    auto g_x_0_x_0_y_xy_0_yy = buffer_1010_pdsd[45];

    auto g_x_0_x_0_y_xy_0_yz = buffer_1010_pdsd[46];

    auto g_x_0_x_0_y_xy_0_zz = buffer_1010_pdsd[47];

    auto g_x_0_x_0_y_xz_0_xx = buffer_1010_pdsd[48];

    auto g_x_0_x_0_y_xz_0_xy = buffer_1010_pdsd[49];

    auto g_x_0_x_0_y_xz_0_xz = buffer_1010_pdsd[50];

    auto g_x_0_x_0_y_xz_0_yy = buffer_1010_pdsd[51];

    auto g_x_0_x_0_y_xz_0_yz = buffer_1010_pdsd[52];

    auto g_x_0_x_0_y_xz_0_zz = buffer_1010_pdsd[53];

    auto g_x_0_x_0_y_yy_0_xx = buffer_1010_pdsd[54];

    auto g_x_0_x_0_y_yy_0_xy = buffer_1010_pdsd[55];

    auto g_x_0_x_0_y_yy_0_xz = buffer_1010_pdsd[56];

    auto g_x_0_x_0_y_yy_0_yy = buffer_1010_pdsd[57];

    auto g_x_0_x_0_y_yy_0_yz = buffer_1010_pdsd[58];

    auto g_x_0_x_0_y_yy_0_zz = buffer_1010_pdsd[59];

    auto g_x_0_x_0_y_yz_0_xx = buffer_1010_pdsd[60];

    auto g_x_0_x_0_y_yz_0_xy = buffer_1010_pdsd[61];

    auto g_x_0_x_0_y_yz_0_xz = buffer_1010_pdsd[62];

    auto g_x_0_x_0_y_yz_0_yy = buffer_1010_pdsd[63];

    auto g_x_0_x_0_y_yz_0_yz = buffer_1010_pdsd[64];

    auto g_x_0_x_0_y_yz_0_zz = buffer_1010_pdsd[65];

    auto g_x_0_x_0_y_zz_0_xx = buffer_1010_pdsd[66];

    auto g_x_0_x_0_y_zz_0_xy = buffer_1010_pdsd[67];

    auto g_x_0_x_0_y_zz_0_xz = buffer_1010_pdsd[68];

    auto g_x_0_x_0_y_zz_0_yy = buffer_1010_pdsd[69];

    auto g_x_0_x_0_y_zz_0_yz = buffer_1010_pdsd[70];

    auto g_x_0_x_0_y_zz_0_zz = buffer_1010_pdsd[71];

    auto g_x_0_x_0_z_xx_0_xx = buffer_1010_pdsd[72];

    auto g_x_0_x_0_z_xx_0_xy = buffer_1010_pdsd[73];

    auto g_x_0_x_0_z_xx_0_xz = buffer_1010_pdsd[74];

    auto g_x_0_x_0_z_xx_0_yy = buffer_1010_pdsd[75];

    auto g_x_0_x_0_z_xx_0_yz = buffer_1010_pdsd[76];

    auto g_x_0_x_0_z_xx_0_zz = buffer_1010_pdsd[77];

    auto g_x_0_x_0_z_xy_0_xx = buffer_1010_pdsd[78];

    auto g_x_0_x_0_z_xy_0_xy = buffer_1010_pdsd[79];

    auto g_x_0_x_0_z_xy_0_xz = buffer_1010_pdsd[80];

    auto g_x_0_x_0_z_xy_0_yy = buffer_1010_pdsd[81];

    auto g_x_0_x_0_z_xy_0_yz = buffer_1010_pdsd[82];

    auto g_x_0_x_0_z_xy_0_zz = buffer_1010_pdsd[83];

    auto g_x_0_x_0_z_xz_0_xx = buffer_1010_pdsd[84];

    auto g_x_0_x_0_z_xz_0_xy = buffer_1010_pdsd[85];

    auto g_x_0_x_0_z_xz_0_xz = buffer_1010_pdsd[86];

    auto g_x_0_x_0_z_xz_0_yy = buffer_1010_pdsd[87];

    auto g_x_0_x_0_z_xz_0_yz = buffer_1010_pdsd[88];

    auto g_x_0_x_0_z_xz_0_zz = buffer_1010_pdsd[89];

    auto g_x_0_x_0_z_yy_0_xx = buffer_1010_pdsd[90];

    auto g_x_0_x_0_z_yy_0_xy = buffer_1010_pdsd[91];

    auto g_x_0_x_0_z_yy_0_xz = buffer_1010_pdsd[92];

    auto g_x_0_x_0_z_yy_0_yy = buffer_1010_pdsd[93];

    auto g_x_0_x_0_z_yy_0_yz = buffer_1010_pdsd[94];

    auto g_x_0_x_0_z_yy_0_zz = buffer_1010_pdsd[95];

    auto g_x_0_x_0_z_yz_0_xx = buffer_1010_pdsd[96];

    auto g_x_0_x_0_z_yz_0_xy = buffer_1010_pdsd[97];

    auto g_x_0_x_0_z_yz_0_xz = buffer_1010_pdsd[98];

    auto g_x_0_x_0_z_yz_0_yy = buffer_1010_pdsd[99];

    auto g_x_0_x_0_z_yz_0_yz = buffer_1010_pdsd[100];

    auto g_x_0_x_0_z_yz_0_zz = buffer_1010_pdsd[101];

    auto g_x_0_x_0_z_zz_0_xx = buffer_1010_pdsd[102];

    auto g_x_0_x_0_z_zz_0_xy = buffer_1010_pdsd[103];

    auto g_x_0_x_0_z_zz_0_xz = buffer_1010_pdsd[104];

    auto g_x_0_x_0_z_zz_0_yy = buffer_1010_pdsd[105];

    auto g_x_0_x_0_z_zz_0_yz = buffer_1010_pdsd[106];

    auto g_x_0_x_0_z_zz_0_zz = buffer_1010_pdsd[107];

    auto g_x_0_y_0_x_xx_0_xx = buffer_1010_pdsd[108];

    auto g_x_0_y_0_x_xx_0_xy = buffer_1010_pdsd[109];

    auto g_x_0_y_0_x_xx_0_xz = buffer_1010_pdsd[110];

    auto g_x_0_y_0_x_xx_0_yy = buffer_1010_pdsd[111];

    auto g_x_0_y_0_x_xx_0_yz = buffer_1010_pdsd[112];

    auto g_x_0_y_0_x_xx_0_zz = buffer_1010_pdsd[113];

    auto g_x_0_y_0_x_xy_0_xx = buffer_1010_pdsd[114];

    auto g_x_0_y_0_x_xy_0_xy = buffer_1010_pdsd[115];

    auto g_x_0_y_0_x_xy_0_xz = buffer_1010_pdsd[116];

    auto g_x_0_y_0_x_xy_0_yy = buffer_1010_pdsd[117];

    auto g_x_0_y_0_x_xy_0_yz = buffer_1010_pdsd[118];

    auto g_x_0_y_0_x_xy_0_zz = buffer_1010_pdsd[119];

    auto g_x_0_y_0_x_xz_0_xx = buffer_1010_pdsd[120];

    auto g_x_0_y_0_x_xz_0_xy = buffer_1010_pdsd[121];

    auto g_x_0_y_0_x_xz_0_xz = buffer_1010_pdsd[122];

    auto g_x_0_y_0_x_xz_0_yy = buffer_1010_pdsd[123];

    auto g_x_0_y_0_x_xz_0_yz = buffer_1010_pdsd[124];

    auto g_x_0_y_0_x_xz_0_zz = buffer_1010_pdsd[125];

    auto g_x_0_y_0_x_yy_0_xx = buffer_1010_pdsd[126];

    auto g_x_0_y_0_x_yy_0_xy = buffer_1010_pdsd[127];

    auto g_x_0_y_0_x_yy_0_xz = buffer_1010_pdsd[128];

    auto g_x_0_y_0_x_yy_0_yy = buffer_1010_pdsd[129];

    auto g_x_0_y_0_x_yy_0_yz = buffer_1010_pdsd[130];

    auto g_x_0_y_0_x_yy_0_zz = buffer_1010_pdsd[131];

    auto g_x_0_y_0_x_yz_0_xx = buffer_1010_pdsd[132];

    auto g_x_0_y_0_x_yz_0_xy = buffer_1010_pdsd[133];

    auto g_x_0_y_0_x_yz_0_xz = buffer_1010_pdsd[134];

    auto g_x_0_y_0_x_yz_0_yy = buffer_1010_pdsd[135];

    auto g_x_0_y_0_x_yz_0_yz = buffer_1010_pdsd[136];

    auto g_x_0_y_0_x_yz_0_zz = buffer_1010_pdsd[137];

    auto g_x_0_y_0_x_zz_0_xx = buffer_1010_pdsd[138];

    auto g_x_0_y_0_x_zz_0_xy = buffer_1010_pdsd[139];

    auto g_x_0_y_0_x_zz_0_xz = buffer_1010_pdsd[140];

    auto g_x_0_y_0_x_zz_0_yy = buffer_1010_pdsd[141];

    auto g_x_0_y_0_x_zz_0_yz = buffer_1010_pdsd[142];

    auto g_x_0_y_0_x_zz_0_zz = buffer_1010_pdsd[143];

    auto g_x_0_y_0_y_xx_0_xx = buffer_1010_pdsd[144];

    auto g_x_0_y_0_y_xx_0_xy = buffer_1010_pdsd[145];

    auto g_x_0_y_0_y_xx_0_xz = buffer_1010_pdsd[146];

    auto g_x_0_y_0_y_xx_0_yy = buffer_1010_pdsd[147];

    auto g_x_0_y_0_y_xx_0_yz = buffer_1010_pdsd[148];

    auto g_x_0_y_0_y_xx_0_zz = buffer_1010_pdsd[149];

    auto g_x_0_y_0_y_xy_0_xx = buffer_1010_pdsd[150];

    auto g_x_0_y_0_y_xy_0_xy = buffer_1010_pdsd[151];

    auto g_x_0_y_0_y_xy_0_xz = buffer_1010_pdsd[152];

    auto g_x_0_y_0_y_xy_0_yy = buffer_1010_pdsd[153];

    auto g_x_0_y_0_y_xy_0_yz = buffer_1010_pdsd[154];

    auto g_x_0_y_0_y_xy_0_zz = buffer_1010_pdsd[155];

    auto g_x_0_y_0_y_xz_0_xx = buffer_1010_pdsd[156];

    auto g_x_0_y_0_y_xz_0_xy = buffer_1010_pdsd[157];

    auto g_x_0_y_0_y_xz_0_xz = buffer_1010_pdsd[158];

    auto g_x_0_y_0_y_xz_0_yy = buffer_1010_pdsd[159];

    auto g_x_0_y_0_y_xz_0_yz = buffer_1010_pdsd[160];

    auto g_x_0_y_0_y_xz_0_zz = buffer_1010_pdsd[161];

    auto g_x_0_y_0_y_yy_0_xx = buffer_1010_pdsd[162];

    auto g_x_0_y_0_y_yy_0_xy = buffer_1010_pdsd[163];

    auto g_x_0_y_0_y_yy_0_xz = buffer_1010_pdsd[164];

    auto g_x_0_y_0_y_yy_0_yy = buffer_1010_pdsd[165];

    auto g_x_0_y_0_y_yy_0_yz = buffer_1010_pdsd[166];

    auto g_x_0_y_0_y_yy_0_zz = buffer_1010_pdsd[167];

    auto g_x_0_y_0_y_yz_0_xx = buffer_1010_pdsd[168];

    auto g_x_0_y_0_y_yz_0_xy = buffer_1010_pdsd[169];

    auto g_x_0_y_0_y_yz_0_xz = buffer_1010_pdsd[170];

    auto g_x_0_y_0_y_yz_0_yy = buffer_1010_pdsd[171];

    auto g_x_0_y_0_y_yz_0_yz = buffer_1010_pdsd[172];

    auto g_x_0_y_0_y_yz_0_zz = buffer_1010_pdsd[173];

    auto g_x_0_y_0_y_zz_0_xx = buffer_1010_pdsd[174];

    auto g_x_0_y_0_y_zz_0_xy = buffer_1010_pdsd[175];

    auto g_x_0_y_0_y_zz_0_xz = buffer_1010_pdsd[176];

    auto g_x_0_y_0_y_zz_0_yy = buffer_1010_pdsd[177];

    auto g_x_0_y_0_y_zz_0_yz = buffer_1010_pdsd[178];

    auto g_x_0_y_0_y_zz_0_zz = buffer_1010_pdsd[179];

    auto g_x_0_y_0_z_xx_0_xx = buffer_1010_pdsd[180];

    auto g_x_0_y_0_z_xx_0_xy = buffer_1010_pdsd[181];

    auto g_x_0_y_0_z_xx_0_xz = buffer_1010_pdsd[182];

    auto g_x_0_y_0_z_xx_0_yy = buffer_1010_pdsd[183];

    auto g_x_0_y_0_z_xx_0_yz = buffer_1010_pdsd[184];

    auto g_x_0_y_0_z_xx_0_zz = buffer_1010_pdsd[185];

    auto g_x_0_y_0_z_xy_0_xx = buffer_1010_pdsd[186];

    auto g_x_0_y_0_z_xy_0_xy = buffer_1010_pdsd[187];

    auto g_x_0_y_0_z_xy_0_xz = buffer_1010_pdsd[188];

    auto g_x_0_y_0_z_xy_0_yy = buffer_1010_pdsd[189];

    auto g_x_0_y_0_z_xy_0_yz = buffer_1010_pdsd[190];

    auto g_x_0_y_0_z_xy_0_zz = buffer_1010_pdsd[191];

    auto g_x_0_y_0_z_xz_0_xx = buffer_1010_pdsd[192];

    auto g_x_0_y_0_z_xz_0_xy = buffer_1010_pdsd[193];

    auto g_x_0_y_0_z_xz_0_xz = buffer_1010_pdsd[194];

    auto g_x_0_y_0_z_xz_0_yy = buffer_1010_pdsd[195];

    auto g_x_0_y_0_z_xz_0_yz = buffer_1010_pdsd[196];

    auto g_x_0_y_0_z_xz_0_zz = buffer_1010_pdsd[197];

    auto g_x_0_y_0_z_yy_0_xx = buffer_1010_pdsd[198];

    auto g_x_0_y_0_z_yy_0_xy = buffer_1010_pdsd[199];

    auto g_x_0_y_0_z_yy_0_xz = buffer_1010_pdsd[200];

    auto g_x_0_y_0_z_yy_0_yy = buffer_1010_pdsd[201];

    auto g_x_0_y_0_z_yy_0_yz = buffer_1010_pdsd[202];

    auto g_x_0_y_0_z_yy_0_zz = buffer_1010_pdsd[203];

    auto g_x_0_y_0_z_yz_0_xx = buffer_1010_pdsd[204];

    auto g_x_0_y_0_z_yz_0_xy = buffer_1010_pdsd[205];

    auto g_x_0_y_0_z_yz_0_xz = buffer_1010_pdsd[206];

    auto g_x_0_y_0_z_yz_0_yy = buffer_1010_pdsd[207];

    auto g_x_0_y_0_z_yz_0_yz = buffer_1010_pdsd[208];

    auto g_x_0_y_0_z_yz_0_zz = buffer_1010_pdsd[209];

    auto g_x_0_y_0_z_zz_0_xx = buffer_1010_pdsd[210];

    auto g_x_0_y_0_z_zz_0_xy = buffer_1010_pdsd[211];

    auto g_x_0_y_0_z_zz_0_xz = buffer_1010_pdsd[212];

    auto g_x_0_y_0_z_zz_0_yy = buffer_1010_pdsd[213];

    auto g_x_0_y_0_z_zz_0_yz = buffer_1010_pdsd[214];

    auto g_x_0_y_0_z_zz_0_zz = buffer_1010_pdsd[215];

    auto g_x_0_z_0_x_xx_0_xx = buffer_1010_pdsd[216];

    auto g_x_0_z_0_x_xx_0_xy = buffer_1010_pdsd[217];

    auto g_x_0_z_0_x_xx_0_xz = buffer_1010_pdsd[218];

    auto g_x_0_z_0_x_xx_0_yy = buffer_1010_pdsd[219];

    auto g_x_0_z_0_x_xx_0_yz = buffer_1010_pdsd[220];

    auto g_x_0_z_0_x_xx_0_zz = buffer_1010_pdsd[221];

    auto g_x_0_z_0_x_xy_0_xx = buffer_1010_pdsd[222];

    auto g_x_0_z_0_x_xy_0_xy = buffer_1010_pdsd[223];

    auto g_x_0_z_0_x_xy_0_xz = buffer_1010_pdsd[224];

    auto g_x_0_z_0_x_xy_0_yy = buffer_1010_pdsd[225];

    auto g_x_0_z_0_x_xy_0_yz = buffer_1010_pdsd[226];

    auto g_x_0_z_0_x_xy_0_zz = buffer_1010_pdsd[227];

    auto g_x_0_z_0_x_xz_0_xx = buffer_1010_pdsd[228];

    auto g_x_0_z_0_x_xz_0_xy = buffer_1010_pdsd[229];

    auto g_x_0_z_0_x_xz_0_xz = buffer_1010_pdsd[230];

    auto g_x_0_z_0_x_xz_0_yy = buffer_1010_pdsd[231];

    auto g_x_0_z_0_x_xz_0_yz = buffer_1010_pdsd[232];

    auto g_x_0_z_0_x_xz_0_zz = buffer_1010_pdsd[233];

    auto g_x_0_z_0_x_yy_0_xx = buffer_1010_pdsd[234];

    auto g_x_0_z_0_x_yy_0_xy = buffer_1010_pdsd[235];

    auto g_x_0_z_0_x_yy_0_xz = buffer_1010_pdsd[236];

    auto g_x_0_z_0_x_yy_0_yy = buffer_1010_pdsd[237];

    auto g_x_0_z_0_x_yy_0_yz = buffer_1010_pdsd[238];

    auto g_x_0_z_0_x_yy_0_zz = buffer_1010_pdsd[239];

    auto g_x_0_z_0_x_yz_0_xx = buffer_1010_pdsd[240];

    auto g_x_0_z_0_x_yz_0_xy = buffer_1010_pdsd[241];

    auto g_x_0_z_0_x_yz_0_xz = buffer_1010_pdsd[242];

    auto g_x_0_z_0_x_yz_0_yy = buffer_1010_pdsd[243];

    auto g_x_0_z_0_x_yz_0_yz = buffer_1010_pdsd[244];

    auto g_x_0_z_0_x_yz_0_zz = buffer_1010_pdsd[245];

    auto g_x_0_z_0_x_zz_0_xx = buffer_1010_pdsd[246];

    auto g_x_0_z_0_x_zz_0_xy = buffer_1010_pdsd[247];

    auto g_x_0_z_0_x_zz_0_xz = buffer_1010_pdsd[248];

    auto g_x_0_z_0_x_zz_0_yy = buffer_1010_pdsd[249];

    auto g_x_0_z_0_x_zz_0_yz = buffer_1010_pdsd[250];

    auto g_x_0_z_0_x_zz_0_zz = buffer_1010_pdsd[251];

    auto g_x_0_z_0_y_xx_0_xx = buffer_1010_pdsd[252];

    auto g_x_0_z_0_y_xx_0_xy = buffer_1010_pdsd[253];

    auto g_x_0_z_0_y_xx_0_xz = buffer_1010_pdsd[254];

    auto g_x_0_z_0_y_xx_0_yy = buffer_1010_pdsd[255];

    auto g_x_0_z_0_y_xx_0_yz = buffer_1010_pdsd[256];

    auto g_x_0_z_0_y_xx_0_zz = buffer_1010_pdsd[257];

    auto g_x_0_z_0_y_xy_0_xx = buffer_1010_pdsd[258];

    auto g_x_0_z_0_y_xy_0_xy = buffer_1010_pdsd[259];

    auto g_x_0_z_0_y_xy_0_xz = buffer_1010_pdsd[260];

    auto g_x_0_z_0_y_xy_0_yy = buffer_1010_pdsd[261];

    auto g_x_0_z_0_y_xy_0_yz = buffer_1010_pdsd[262];

    auto g_x_0_z_0_y_xy_0_zz = buffer_1010_pdsd[263];

    auto g_x_0_z_0_y_xz_0_xx = buffer_1010_pdsd[264];

    auto g_x_0_z_0_y_xz_0_xy = buffer_1010_pdsd[265];

    auto g_x_0_z_0_y_xz_0_xz = buffer_1010_pdsd[266];

    auto g_x_0_z_0_y_xz_0_yy = buffer_1010_pdsd[267];

    auto g_x_0_z_0_y_xz_0_yz = buffer_1010_pdsd[268];

    auto g_x_0_z_0_y_xz_0_zz = buffer_1010_pdsd[269];

    auto g_x_0_z_0_y_yy_0_xx = buffer_1010_pdsd[270];

    auto g_x_0_z_0_y_yy_0_xy = buffer_1010_pdsd[271];

    auto g_x_0_z_0_y_yy_0_xz = buffer_1010_pdsd[272];

    auto g_x_0_z_0_y_yy_0_yy = buffer_1010_pdsd[273];

    auto g_x_0_z_0_y_yy_0_yz = buffer_1010_pdsd[274];

    auto g_x_0_z_0_y_yy_0_zz = buffer_1010_pdsd[275];

    auto g_x_0_z_0_y_yz_0_xx = buffer_1010_pdsd[276];

    auto g_x_0_z_0_y_yz_0_xy = buffer_1010_pdsd[277];

    auto g_x_0_z_0_y_yz_0_xz = buffer_1010_pdsd[278];

    auto g_x_0_z_0_y_yz_0_yy = buffer_1010_pdsd[279];

    auto g_x_0_z_0_y_yz_0_yz = buffer_1010_pdsd[280];

    auto g_x_0_z_0_y_yz_0_zz = buffer_1010_pdsd[281];

    auto g_x_0_z_0_y_zz_0_xx = buffer_1010_pdsd[282];

    auto g_x_0_z_0_y_zz_0_xy = buffer_1010_pdsd[283];

    auto g_x_0_z_0_y_zz_0_xz = buffer_1010_pdsd[284];

    auto g_x_0_z_0_y_zz_0_yy = buffer_1010_pdsd[285];

    auto g_x_0_z_0_y_zz_0_yz = buffer_1010_pdsd[286];

    auto g_x_0_z_0_y_zz_0_zz = buffer_1010_pdsd[287];

    auto g_x_0_z_0_z_xx_0_xx = buffer_1010_pdsd[288];

    auto g_x_0_z_0_z_xx_0_xy = buffer_1010_pdsd[289];

    auto g_x_0_z_0_z_xx_0_xz = buffer_1010_pdsd[290];

    auto g_x_0_z_0_z_xx_0_yy = buffer_1010_pdsd[291];

    auto g_x_0_z_0_z_xx_0_yz = buffer_1010_pdsd[292];

    auto g_x_0_z_0_z_xx_0_zz = buffer_1010_pdsd[293];

    auto g_x_0_z_0_z_xy_0_xx = buffer_1010_pdsd[294];

    auto g_x_0_z_0_z_xy_0_xy = buffer_1010_pdsd[295];

    auto g_x_0_z_0_z_xy_0_xz = buffer_1010_pdsd[296];

    auto g_x_0_z_0_z_xy_0_yy = buffer_1010_pdsd[297];

    auto g_x_0_z_0_z_xy_0_yz = buffer_1010_pdsd[298];

    auto g_x_0_z_0_z_xy_0_zz = buffer_1010_pdsd[299];

    auto g_x_0_z_0_z_xz_0_xx = buffer_1010_pdsd[300];

    auto g_x_0_z_0_z_xz_0_xy = buffer_1010_pdsd[301];

    auto g_x_0_z_0_z_xz_0_xz = buffer_1010_pdsd[302];

    auto g_x_0_z_0_z_xz_0_yy = buffer_1010_pdsd[303];

    auto g_x_0_z_0_z_xz_0_yz = buffer_1010_pdsd[304];

    auto g_x_0_z_0_z_xz_0_zz = buffer_1010_pdsd[305];

    auto g_x_0_z_0_z_yy_0_xx = buffer_1010_pdsd[306];

    auto g_x_0_z_0_z_yy_0_xy = buffer_1010_pdsd[307];

    auto g_x_0_z_0_z_yy_0_xz = buffer_1010_pdsd[308];

    auto g_x_0_z_0_z_yy_0_yy = buffer_1010_pdsd[309];

    auto g_x_0_z_0_z_yy_0_yz = buffer_1010_pdsd[310];

    auto g_x_0_z_0_z_yy_0_zz = buffer_1010_pdsd[311];

    auto g_x_0_z_0_z_yz_0_xx = buffer_1010_pdsd[312];

    auto g_x_0_z_0_z_yz_0_xy = buffer_1010_pdsd[313];

    auto g_x_0_z_0_z_yz_0_xz = buffer_1010_pdsd[314];

    auto g_x_0_z_0_z_yz_0_yy = buffer_1010_pdsd[315];

    auto g_x_0_z_0_z_yz_0_yz = buffer_1010_pdsd[316];

    auto g_x_0_z_0_z_yz_0_zz = buffer_1010_pdsd[317];

    auto g_x_0_z_0_z_zz_0_xx = buffer_1010_pdsd[318];

    auto g_x_0_z_0_z_zz_0_xy = buffer_1010_pdsd[319];

    auto g_x_0_z_0_z_zz_0_xz = buffer_1010_pdsd[320];

    auto g_x_0_z_0_z_zz_0_yy = buffer_1010_pdsd[321];

    auto g_x_0_z_0_z_zz_0_yz = buffer_1010_pdsd[322];

    auto g_x_0_z_0_z_zz_0_zz = buffer_1010_pdsd[323];

    auto g_y_0_x_0_x_xx_0_xx = buffer_1010_pdsd[324];

    auto g_y_0_x_0_x_xx_0_xy = buffer_1010_pdsd[325];

    auto g_y_0_x_0_x_xx_0_xz = buffer_1010_pdsd[326];

    auto g_y_0_x_0_x_xx_0_yy = buffer_1010_pdsd[327];

    auto g_y_0_x_0_x_xx_0_yz = buffer_1010_pdsd[328];

    auto g_y_0_x_0_x_xx_0_zz = buffer_1010_pdsd[329];

    auto g_y_0_x_0_x_xy_0_xx = buffer_1010_pdsd[330];

    auto g_y_0_x_0_x_xy_0_xy = buffer_1010_pdsd[331];

    auto g_y_0_x_0_x_xy_0_xz = buffer_1010_pdsd[332];

    auto g_y_0_x_0_x_xy_0_yy = buffer_1010_pdsd[333];

    auto g_y_0_x_0_x_xy_0_yz = buffer_1010_pdsd[334];

    auto g_y_0_x_0_x_xy_0_zz = buffer_1010_pdsd[335];

    auto g_y_0_x_0_x_xz_0_xx = buffer_1010_pdsd[336];

    auto g_y_0_x_0_x_xz_0_xy = buffer_1010_pdsd[337];

    auto g_y_0_x_0_x_xz_0_xz = buffer_1010_pdsd[338];

    auto g_y_0_x_0_x_xz_0_yy = buffer_1010_pdsd[339];

    auto g_y_0_x_0_x_xz_0_yz = buffer_1010_pdsd[340];

    auto g_y_0_x_0_x_xz_0_zz = buffer_1010_pdsd[341];

    auto g_y_0_x_0_x_yy_0_xx = buffer_1010_pdsd[342];

    auto g_y_0_x_0_x_yy_0_xy = buffer_1010_pdsd[343];

    auto g_y_0_x_0_x_yy_0_xz = buffer_1010_pdsd[344];

    auto g_y_0_x_0_x_yy_0_yy = buffer_1010_pdsd[345];

    auto g_y_0_x_0_x_yy_0_yz = buffer_1010_pdsd[346];

    auto g_y_0_x_0_x_yy_0_zz = buffer_1010_pdsd[347];

    auto g_y_0_x_0_x_yz_0_xx = buffer_1010_pdsd[348];

    auto g_y_0_x_0_x_yz_0_xy = buffer_1010_pdsd[349];

    auto g_y_0_x_0_x_yz_0_xz = buffer_1010_pdsd[350];

    auto g_y_0_x_0_x_yz_0_yy = buffer_1010_pdsd[351];

    auto g_y_0_x_0_x_yz_0_yz = buffer_1010_pdsd[352];

    auto g_y_0_x_0_x_yz_0_zz = buffer_1010_pdsd[353];

    auto g_y_0_x_0_x_zz_0_xx = buffer_1010_pdsd[354];

    auto g_y_0_x_0_x_zz_0_xy = buffer_1010_pdsd[355];

    auto g_y_0_x_0_x_zz_0_xz = buffer_1010_pdsd[356];

    auto g_y_0_x_0_x_zz_0_yy = buffer_1010_pdsd[357];

    auto g_y_0_x_0_x_zz_0_yz = buffer_1010_pdsd[358];

    auto g_y_0_x_0_x_zz_0_zz = buffer_1010_pdsd[359];

    auto g_y_0_x_0_y_xx_0_xx = buffer_1010_pdsd[360];

    auto g_y_0_x_0_y_xx_0_xy = buffer_1010_pdsd[361];

    auto g_y_0_x_0_y_xx_0_xz = buffer_1010_pdsd[362];

    auto g_y_0_x_0_y_xx_0_yy = buffer_1010_pdsd[363];

    auto g_y_0_x_0_y_xx_0_yz = buffer_1010_pdsd[364];

    auto g_y_0_x_0_y_xx_0_zz = buffer_1010_pdsd[365];

    auto g_y_0_x_0_y_xy_0_xx = buffer_1010_pdsd[366];

    auto g_y_0_x_0_y_xy_0_xy = buffer_1010_pdsd[367];

    auto g_y_0_x_0_y_xy_0_xz = buffer_1010_pdsd[368];

    auto g_y_0_x_0_y_xy_0_yy = buffer_1010_pdsd[369];

    auto g_y_0_x_0_y_xy_0_yz = buffer_1010_pdsd[370];

    auto g_y_0_x_0_y_xy_0_zz = buffer_1010_pdsd[371];

    auto g_y_0_x_0_y_xz_0_xx = buffer_1010_pdsd[372];

    auto g_y_0_x_0_y_xz_0_xy = buffer_1010_pdsd[373];

    auto g_y_0_x_0_y_xz_0_xz = buffer_1010_pdsd[374];

    auto g_y_0_x_0_y_xz_0_yy = buffer_1010_pdsd[375];

    auto g_y_0_x_0_y_xz_0_yz = buffer_1010_pdsd[376];

    auto g_y_0_x_0_y_xz_0_zz = buffer_1010_pdsd[377];

    auto g_y_0_x_0_y_yy_0_xx = buffer_1010_pdsd[378];

    auto g_y_0_x_0_y_yy_0_xy = buffer_1010_pdsd[379];

    auto g_y_0_x_0_y_yy_0_xz = buffer_1010_pdsd[380];

    auto g_y_0_x_0_y_yy_0_yy = buffer_1010_pdsd[381];

    auto g_y_0_x_0_y_yy_0_yz = buffer_1010_pdsd[382];

    auto g_y_0_x_0_y_yy_0_zz = buffer_1010_pdsd[383];

    auto g_y_0_x_0_y_yz_0_xx = buffer_1010_pdsd[384];

    auto g_y_0_x_0_y_yz_0_xy = buffer_1010_pdsd[385];

    auto g_y_0_x_0_y_yz_0_xz = buffer_1010_pdsd[386];

    auto g_y_0_x_0_y_yz_0_yy = buffer_1010_pdsd[387];

    auto g_y_0_x_0_y_yz_0_yz = buffer_1010_pdsd[388];

    auto g_y_0_x_0_y_yz_0_zz = buffer_1010_pdsd[389];

    auto g_y_0_x_0_y_zz_0_xx = buffer_1010_pdsd[390];

    auto g_y_0_x_0_y_zz_0_xy = buffer_1010_pdsd[391];

    auto g_y_0_x_0_y_zz_0_xz = buffer_1010_pdsd[392];

    auto g_y_0_x_0_y_zz_0_yy = buffer_1010_pdsd[393];

    auto g_y_0_x_0_y_zz_0_yz = buffer_1010_pdsd[394];

    auto g_y_0_x_0_y_zz_0_zz = buffer_1010_pdsd[395];

    auto g_y_0_x_0_z_xx_0_xx = buffer_1010_pdsd[396];

    auto g_y_0_x_0_z_xx_0_xy = buffer_1010_pdsd[397];

    auto g_y_0_x_0_z_xx_0_xz = buffer_1010_pdsd[398];

    auto g_y_0_x_0_z_xx_0_yy = buffer_1010_pdsd[399];

    auto g_y_0_x_0_z_xx_0_yz = buffer_1010_pdsd[400];

    auto g_y_0_x_0_z_xx_0_zz = buffer_1010_pdsd[401];

    auto g_y_0_x_0_z_xy_0_xx = buffer_1010_pdsd[402];

    auto g_y_0_x_0_z_xy_0_xy = buffer_1010_pdsd[403];

    auto g_y_0_x_0_z_xy_0_xz = buffer_1010_pdsd[404];

    auto g_y_0_x_0_z_xy_0_yy = buffer_1010_pdsd[405];

    auto g_y_0_x_0_z_xy_0_yz = buffer_1010_pdsd[406];

    auto g_y_0_x_0_z_xy_0_zz = buffer_1010_pdsd[407];

    auto g_y_0_x_0_z_xz_0_xx = buffer_1010_pdsd[408];

    auto g_y_0_x_0_z_xz_0_xy = buffer_1010_pdsd[409];

    auto g_y_0_x_0_z_xz_0_xz = buffer_1010_pdsd[410];

    auto g_y_0_x_0_z_xz_0_yy = buffer_1010_pdsd[411];

    auto g_y_0_x_0_z_xz_0_yz = buffer_1010_pdsd[412];

    auto g_y_0_x_0_z_xz_0_zz = buffer_1010_pdsd[413];

    auto g_y_0_x_0_z_yy_0_xx = buffer_1010_pdsd[414];

    auto g_y_0_x_0_z_yy_0_xy = buffer_1010_pdsd[415];

    auto g_y_0_x_0_z_yy_0_xz = buffer_1010_pdsd[416];

    auto g_y_0_x_0_z_yy_0_yy = buffer_1010_pdsd[417];

    auto g_y_0_x_0_z_yy_0_yz = buffer_1010_pdsd[418];

    auto g_y_0_x_0_z_yy_0_zz = buffer_1010_pdsd[419];

    auto g_y_0_x_0_z_yz_0_xx = buffer_1010_pdsd[420];

    auto g_y_0_x_0_z_yz_0_xy = buffer_1010_pdsd[421];

    auto g_y_0_x_0_z_yz_0_xz = buffer_1010_pdsd[422];

    auto g_y_0_x_0_z_yz_0_yy = buffer_1010_pdsd[423];

    auto g_y_0_x_0_z_yz_0_yz = buffer_1010_pdsd[424];

    auto g_y_0_x_0_z_yz_0_zz = buffer_1010_pdsd[425];

    auto g_y_0_x_0_z_zz_0_xx = buffer_1010_pdsd[426];

    auto g_y_0_x_0_z_zz_0_xy = buffer_1010_pdsd[427];

    auto g_y_0_x_0_z_zz_0_xz = buffer_1010_pdsd[428];

    auto g_y_0_x_0_z_zz_0_yy = buffer_1010_pdsd[429];

    auto g_y_0_x_0_z_zz_0_yz = buffer_1010_pdsd[430];

    auto g_y_0_x_0_z_zz_0_zz = buffer_1010_pdsd[431];

    auto g_y_0_y_0_x_xx_0_xx = buffer_1010_pdsd[432];

    auto g_y_0_y_0_x_xx_0_xy = buffer_1010_pdsd[433];

    auto g_y_0_y_0_x_xx_0_xz = buffer_1010_pdsd[434];

    auto g_y_0_y_0_x_xx_0_yy = buffer_1010_pdsd[435];

    auto g_y_0_y_0_x_xx_0_yz = buffer_1010_pdsd[436];

    auto g_y_0_y_0_x_xx_0_zz = buffer_1010_pdsd[437];

    auto g_y_0_y_0_x_xy_0_xx = buffer_1010_pdsd[438];

    auto g_y_0_y_0_x_xy_0_xy = buffer_1010_pdsd[439];

    auto g_y_0_y_0_x_xy_0_xz = buffer_1010_pdsd[440];

    auto g_y_0_y_0_x_xy_0_yy = buffer_1010_pdsd[441];

    auto g_y_0_y_0_x_xy_0_yz = buffer_1010_pdsd[442];

    auto g_y_0_y_0_x_xy_0_zz = buffer_1010_pdsd[443];

    auto g_y_0_y_0_x_xz_0_xx = buffer_1010_pdsd[444];

    auto g_y_0_y_0_x_xz_0_xy = buffer_1010_pdsd[445];

    auto g_y_0_y_0_x_xz_0_xz = buffer_1010_pdsd[446];

    auto g_y_0_y_0_x_xz_0_yy = buffer_1010_pdsd[447];

    auto g_y_0_y_0_x_xz_0_yz = buffer_1010_pdsd[448];

    auto g_y_0_y_0_x_xz_0_zz = buffer_1010_pdsd[449];

    auto g_y_0_y_0_x_yy_0_xx = buffer_1010_pdsd[450];

    auto g_y_0_y_0_x_yy_0_xy = buffer_1010_pdsd[451];

    auto g_y_0_y_0_x_yy_0_xz = buffer_1010_pdsd[452];

    auto g_y_0_y_0_x_yy_0_yy = buffer_1010_pdsd[453];

    auto g_y_0_y_0_x_yy_0_yz = buffer_1010_pdsd[454];

    auto g_y_0_y_0_x_yy_0_zz = buffer_1010_pdsd[455];

    auto g_y_0_y_0_x_yz_0_xx = buffer_1010_pdsd[456];

    auto g_y_0_y_0_x_yz_0_xy = buffer_1010_pdsd[457];

    auto g_y_0_y_0_x_yz_0_xz = buffer_1010_pdsd[458];

    auto g_y_0_y_0_x_yz_0_yy = buffer_1010_pdsd[459];

    auto g_y_0_y_0_x_yz_0_yz = buffer_1010_pdsd[460];

    auto g_y_0_y_0_x_yz_0_zz = buffer_1010_pdsd[461];

    auto g_y_0_y_0_x_zz_0_xx = buffer_1010_pdsd[462];

    auto g_y_0_y_0_x_zz_0_xy = buffer_1010_pdsd[463];

    auto g_y_0_y_0_x_zz_0_xz = buffer_1010_pdsd[464];

    auto g_y_0_y_0_x_zz_0_yy = buffer_1010_pdsd[465];

    auto g_y_0_y_0_x_zz_0_yz = buffer_1010_pdsd[466];

    auto g_y_0_y_0_x_zz_0_zz = buffer_1010_pdsd[467];

    auto g_y_0_y_0_y_xx_0_xx = buffer_1010_pdsd[468];

    auto g_y_0_y_0_y_xx_0_xy = buffer_1010_pdsd[469];

    auto g_y_0_y_0_y_xx_0_xz = buffer_1010_pdsd[470];

    auto g_y_0_y_0_y_xx_0_yy = buffer_1010_pdsd[471];

    auto g_y_0_y_0_y_xx_0_yz = buffer_1010_pdsd[472];

    auto g_y_0_y_0_y_xx_0_zz = buffer_1010_pdsd[473];

    auto g_y_0_y_0_y_xy_0_xx = buffer_1010_pdsd[474];

    auto g_y_0_y_0_y_xy_0_xy = buffer_1010_pdsd[475];

    auto g_y_0_y_0_y_xy_0_xz = buffer_1010_pdsd[476];

    auto g_y_0_y_0_y_xy_0_yy = buffer_1010_pdsd[477];

    auto g_y_0_y_0_y_xy_0_yz = buffer_1010_pdsd[478];

    auto g_y_0_y_0_y_xy_0_zz = buffer_1010_pdsd[479];

    auto g_y_0_y_0_y_xz_0_xx = buffer_1010_pdsd[480];

    auto g_y_0_y_0_y_xz_0_xy = buffer_1010_pdsd[481];

    auto g_y_0_y_0_y_xz_0_xz = buffer_1010_pdsd[482];

    auto g_y_0_y_0_y_xz_0_yy = buffer_1010_pdsd[483];

    auto g_y_0_y_0_y_xz_0_yz = buffer_1010_pdsd[484];

    auto g_y_0_y_0_y_xz_0_zz = buffer_1010_pdsd[485];

    auto g_y_0_y_0_y_yy_0_xx = buffer_1010_pdsd[486];

    auto g_y_0_y_0_y_yy_0_xy = buffer_1010_pdsd[487];

    auto g_y_0_y_0_y_yy_0_xz = buffer_1010_pdsd[488];

    auto g_y_0_y_0_y_yy_0_yy = buffer_1010_pdsd[489];

    auto g_y_0_y_0_y_yy_0_yz = buffer_1010_pdsd[490];

    auto g_y_0_y_0_y_yy_0_zz = buffer_1010_pdsd[491];

    auto g_y_0_y_0_y_yz_0_xx = buffer_1010_pdsd[492];

    auto g_y_0_y_0_y_yz_0_xy = buffer_1010_pdsd[493];

    auto g_y_0_y_0_y_yz_0_xz = buffer_1010_pdsd[494];

    auto g_y_0_y_0_y_yz_0_yy = buffer_1010_pdsd[495];

    auto g_y_0_y_0_y_yz_0_yz = buffer_1010_pdsd[496];

    auto g_y_0_y_0_y_yz_0_zz = buffer_1010_pdsd[497];

    auto g_y_0_y_0_y_zz_0_xx = buffer_1010_pdsd[498];

    auto g_y_0_y_0_y_zz_0_xy = buffer_1010_pdsd[499];

    auto g_y_0_y_0_y_zz_0_xz = buffer_1010_pdsd[500];

    auto g_y_0_y_0_y_zz_0_yy = buffer_1010_pdsd[501];

    auto g_y_0_y_0_y_zz_0_yz = buffer_1010_pdsd[502];

    auto g_y_0_y_0_y_zz_0_zz = buffer_1010_pdsd[503];

    auto g_y_0_y_0_z_xx_0_xx = buffer_1010_pdsd[504];

    auto g_y_0_y_0_z_xx_0_xy = buffer_1010_pdsd[505];

    auto g_y_0_y_0_z_xx_0_xz = buffer_1010_pdsd[506];

    auto g_y_0_y_0_z_xx_0_yy = buffer_1010_pdsd[507];

    auto g_y_0_y_0_z_xx_0_yz = buffer_1010_pdsd[508];

    auto g_y_0_y_0_z_xx_0_zz = buffer_1010_pdsd[509];

    auto g_y_0_y_0_z_xy_0_xx = buffer_1010_pdsd[510];

    auto g_y_0_y_0_z_xy_0_xy = buffer_1010_pdsd[511];

    auto g_y_0_y_0_z_xy_0_xz = buffer_1010_pdsd[512];

    auto g_y_0_y_0_z_xy_0_yy = buffer_1010_pdsd[513];

    auto g_y_0_y_0_z_xy_0_yz = buffer_1010_pdsd[514];

    auto g_y_0_y_0_z_xy_0_zz = buffer_1010_pdsd[515];

    auto g_y_0_y_0_z_xz_0_xx = buffer_1010_pdsd[516];

    auto g_y_0_y_0_z_xz_0_xy = buffer_1010_pdsd[517];

    auto g_y_0_y_0_z_xz_0_xz = buffer_1010_pdsd[518];

    auto g_y_0_y_0_z_xz_0_yy = buffer_1010_pdsd[519];

    auto g_y_0_y_0_z_xz_0_yz = buffer_1010_pdsd[520];

    auto g_y_0_y_0_z_xz_0_zz = buffer_1010_pdsd[521];

    auto g_y_0_y_0_z_yy_0_xx = buffer_1010_pdsd[522];

    auto g_y_0_y_0_z_yy_0_xy = buffer_1010_pdsd[523];

    auto g_y_0_y_0_z_yy_0_xz = buffer_1010_pdsd[524];

    auto g_y_0_y_0_z_yy_0_yy = buffer_1010_pdsd[525];

    auto g_y_0_y_0_z_yy_0_yz = buffer_1010_pdsd[526];

    auto g_y_0_y_0_z_yy_0_zz = buffer_1010_pdsd[527];

    auto g_y_0_y_0_z_yz_0_xx = buffer_1010_pdsd[528];

    auto g_y_0_y_0_z_yz_0_xy = buffer_1010_pdsd[529];

    auto g_y_0_y_0_z_yz_0_xz = buffer_1010_pdsd[530];

    auto g_y_0_y_0_z_yz_0_yy = buffer_1010_pdsd[531];

    auto g_y_0_y_0_z_yz_0_yz = buffer_1010_pdsd[532];

    auto g_y_0_y_0_z_yz_0_zz = buffer_1010_pdsd[533];

    auto g_y_0_y_0_z_zz_0_xx = buffer_1010_pdsd[534];

    auto g_y_0_y_0_z_zz_0_xy = buffer_1010_pdsd[535];

    auto g_y_0_y_0_z_zz_0_xz = buffer_1010_pdsd[536];

    auto g_y_0_y_0_z_zz_0_yy = buffer_1010_pdsd[537];

    auto g_y_0_y_0_z_zz_0_yz = buffer_1010_pdsd[538];

    auto g_y_0_y_0_z_zz_0_zz = buffer_1010_pdsd[539];

    auto g_y_0_z_0_x_xx_0_xx = buffer_1010_pdsd[540];

    auto g_y_0_z_0_x_xx_0_xy = buffer_1010_pdsd[541];

    auto g_y_0_z_0_x_xx_0_xz = buffer_1010_pdsd[542];

    auto g_y_0_z_0_x_xx_0_yy = buffer_1010_pdsd[543];

    auto g_y_0_z_0_x_xx_0_yz = buffer_1010_pdsd[544];

    auto g_y_0_z_0_x_xx_0_zz = buffer_1010_pdsd[545];

    auto g_y_0_z_0_x_xy_0_xx = buffer_1010_pdsd[546];

    auto g_y_0_z_0_x_xy_0_xy = buffer_1010_pdsd[547];

    auto g_y_0_z_0_x_xy_0_xz = buffer_1010_pdsd[548];

    auto g_y_0_z_0_x_xy_0_yy = buffer_1010_pdsd[549];

    auto g_y_0_z_0_x_xy_0_yz = buffer_1010_pdsd[550];

    auto g_y_0_z_0_x_xy_0_zz = buffer_1010_pdsd[551];

    auto g_y_0_z_0_x_xz_0_xx = buffer_1010_pdsd[552];

    auto g_y_0_z_0_x_xz_0_xy = buffer_1010_pdsd[553];

    auto g_y_0_z_0_x_xz_0_xz = buffer_1010_pdsd[554];

    auto g_y_0_z_0_x_xz_0_yy = buffer_1010_pdsd[555];

    auto g_y_0_z_0_x_xz_0_yz = buffer_1010_pdsd[556];

    auto g_y_0_z_0_x_xz_0_zz = buffer_1010_pdsd[557];

    auto g_y_0_z_0_x_yy_0_xx = buffer_1010_pdsd[558];

    auto g_y_0_z_0_x_yy_0_xy = buffer_1010_pdsd[559];

    auto g_y_0_z_0_x_yy_0_xz = buffer_1010_pdsd[560];

    auto g_y_0_z_0_x_yy_0_yy = buffer_1010_pdsd[561];

    auto g_y_0_z_0_x_yy_0_yz = buffer_1010_pdsd[562];

    auto g_y_0_z_0_x_yy_0_zz = buffer_1010_pdsd[563];

    auto g_y_0_z_0_x_yz_0_xx = buffer_1010_pdsd[564];

    auto g_y_0_z_0_x_yz_0_xy = buffer_1010_pdsd[565];

    auto g_y_0_z_0_x_yz_0_xz = buffer_1010_pdsd[566];

    auto g_y_0_z_0_x_yz_0_yy = buffer_1010_pdsd[567];

    auto g_y_0_z_0_x_yz_0_yz = buffer_1010_pdsd[568];

    auto g_y_0_z_0_x_yz_0_zz = buffer_1010_pdsd[569];

    auto g_y_0_z_0_x_zz_0_xx = buffer_1010_pdsd[570];

    auto g_y_0_z_0_x_zz_0_xy = buffer_1010_pdsd[571];

    auto g_y_0_z_0_x_zz_0_xz = buffer_1010_pdsd[572];

    auto g_y_0_z_0_x_zz_0_yy = buffer_1010_pdsd[573];

    auto g_y_0_z_0_x_zz_0_yz = buffer_1010_pdsd[574];

    auto g_y_0_z_0_x_zz_0_zz = buffer_1010_pdsd[575];

    auto g_y_0_z_0_y_xx_0_xx = buffer_1010_pdsd[576];

    auto g_y_0_z_0_y_xx_0_xy = buffer_1010_pdsd[577];

    auto g_y_0_z_0_y_xx_0_xz = buffer_1010_pdsd[578];

    auto g_y_0_z_0_y_xx_0_yy = buffer_1010_pdsd[579];

    auto g_y_0_z_0_y_xx_0_yz = buffer_1010_pdsd[580];

    auto g_y_0_z_0_y_xx_0_zz = buffer_1010_pdsd[581];

    auto g_y_0_z_0_y_xy_0_xx = buffer_1010_pdsd[582];

    auto g_y_0_z_0_y_xy_0_xy = buffer_1010_pdsd[583];

    auto g_y_0_z_0_y_xy_0_xz = buffer_1010_pdsd[584];

    auto g_y_0_z_0_y_xy_0_yy = buffer_1010_pdsd[585];

    auto g_y_0_z_0_y_xy_0_yz = buffer_1010_pdsd[586];

    auto g_y_0_z_0_y_xy_0_zz = buffer_1010_pdsd[587];

    auto g_y_0_z_0_y_xz_0_xx = buffer_1010_pdsd[588];

    auto g_y_0_z_0_y_xz_0_xy = buffer_1010_pdsd[589];

    auto g_y_0_z_0_y_xz_0_xz = buffer_1010_pdsd[590];

    auto g_y_0_z_0_y_xz_0_yy = buffer_1010_pdsd[591];

    auto g_y_0_z_0_y_xz_0_yz = buffer_1010_pdsd[592];

    auto g_y_0_z_0_y_xz_0_zz = buffer_1010_pdsd[593];

    auto g_y_0_z_0_y_yy_0_xx = buffer_1010_pdsd[594];

    auto g_y_0_z_0_y_yy_0_xy = buffer_1010_pdsd[595];

    auto g_y_0_z_0_y_yy_0_xz = buffer_1010_pdsd[596];

    auto g_y_0_z_0_y_yy_0_yy = buffer_1010_pdsd[597];

    auto g_y_0_z_0_y_yy_0_yz = buffer_1010_pdsd[598];

    auto g_y_0_z_0_y_yy_0_zz = buffer_1010_pdsd[599];

    auto g_y_0_z_0_y_yz_0_xx = buffer_1010_pdsd[600];

    auto g_y_0_z_0_y_yz_0_xy = buffer_1010_pdsd[601];

    auto g_y_0_z_0_y_yz_0_xz = buffer_1010_pdsd[602];

    auto g_y_0_z_0_y_yz_0_yy = buffer_1010_pdsd[603];

    auto g_y_0_z_0_y_yz_0_yz = buffer_1010_pdsd[604];

    auto g_y_0_z_0_y_yz_0_zz = buffer_1010_pdsd[605];

    auto g_y_0_z_0_y_zz_0_xx = buffer_1010_pdsd[606];

    auto g_y_0_z_0_y_zz_0_xy = buffer_1010_pdsd[607];

    auto g_y_0_z_0_y_zz_0_xz = buffer_1010_pdsd[608];

    auto g_y_0_z_0_y_zz_0_yy = buffer_1010_pdsd[609];

    auto g_y_0_z_0_y_zz_0_yz = buffer_1010_pdsd[610];

    auto g_y_0_z_0_y_zz_0_zz = buffer_1010_pdsd[611];

    auto g_y_0_z_0_z_xx_0_xx = buffer_1010_pdsd[612];

    auto g_y_0_z_0_z_xx_0_xy = buffer_1010_pdsd[613];

    auto g_y_0_z_0_z_xx_0_xz = buffer_1010_pdsd[614];

    auto g_y_0_z_0_z_xx_0_yy = buffer_1010_pdsd[615];

    auto g_y_0_z_0_z_xx_0_yz = buffer_1010_pdsd[616];

    auto g_y_0_z_0_z_xx_0_zz = buffer_1010_pdsd[617];

    auto g_y_0_z_0_z_xy_0_xx = buffer_1010_pdsd[618];

    auto g_y_0_z_0_z_xy_0_xy = buffer_1010_pdsd[619];

    auto g_y_0_z_0_z_xy_0_xz = buffer_1010_pdsd[620];

    auto g_y_0_z_0_z_xy_0_yy = buffer_1010_pdsd[621];

    auto g_y_0_z_0_z_xy_0_yz = buffer_1010_pdsd[622];

    auto g_y_0_z_0_z_xy_0_zz = buffer_1010_pdsd[623];

    auto g_y_0_z_0_z_xz_0_xx = buffer_1010_pdsd[624];

    auto g_y_0_z_0_z_xz_0_xy = buffer_1010_pdsd[625];

    auto g_y_0_z_0_z_xz_0_xz = buffer_1010_pdsd[626];

    auto g_y_0_z_0_z_xz_0_yy = buffer_1010_pdsd[627];

    auto g_y_0_z_0_z_xz_0_yz = buffer_1010_pdsd[628];

    auto g_y_0_z_0_z_xz_0_zz = buffer_1010_pdsd[629];

    auto g_y_0_z_0_z_yy_0_xx = buffer_1010_pdsd[630];

    auto g_y_0_z_0_z_yy_0_xy = buffer_1010_pdsd[631];

    auto g_y_0_z_0_z_yy_0_xz = buffer_1010_pdsd[632];

    auto g_y_0_z_0_z_yy_0_yy = buffer_1010_pdsd[633];

    auto g_y_0_z_0_z_yy_0_yz = buffer_1010_pdsd[634];

    auto g_y_0_z_0_z_yy_0_zz = buffer_1010_pdsd[635];

    auto g_y_0_z_0_z_yz_0_xx = buffer_1010_pdsd[636];

    auto g_y_0_z_0_z_yz_0_xy = buffer_1010_pdsd[637];

    auto g_y_0_z_0_z_yz_0_xz = buffer_1010_pdsd[638];

    auto g_y_0_z_0_z_yz_0_yy = buffer_1010_pdsd[639];

    auto g_y_0_z_0_z_yz_0_yz = buffer_1010_pdsd[640];

    auto g_y_0_z_0_z_yz_0_zz = buffer_1010_pdsd[641];

    auto g_y_0_z_0_z_zz_0_xx = buffer_1010_pdsd[642];

    auto g_y_0_z_0_z_zz_0_xy = buffer_1010_pdsd[643];

    auto g_y_0_z_0_z_zz_0_xz = buffer_1010_pdsd[644];

    auto g_y_0_z_0_z_zz_0_yy = buffer_1010_pdsd[645];

    auto g_y_0_z_0_z_zz_0_yz = buffer_1010_pdsd[646];

    auto g_y_0_z_0_z_zz_0_zz = buffer_1010_pdsd[647];

    auto g_z_0_x_0_x_xx_0_xx = buffer_1010_pdsd[648];

    auto g_z_0_x_0_x_xx_0_xy = buffer_1010_pdsd[649];

    auto g_z_0_x_0_x_xx_0_xz = buffer_1010_pdsd[650];

    auto g_z_0_x_0_x_xx_0_yy = buffer_1010_pdsd[651];

    auto g_z_0_x_0_x_xx_0_yz = buffer_1010_pdsd[652];

    auto g_z_0_x_0_x_xx_0_zz = buffer_1010_pdsd[653];

    auto g_z_0_x_0_x_xy_0_xx = buffer_1010_pdsd[654];

    auto g_z_0_x_0_x_xy_0_xy = buffer_1010_pdsd[655];

    auto g_z_0_x_0_x_xy_0_xz = buffer_1010_pdsd[656];

    auto g_z_0_x_0_x_xy_0_yy = buffer_1010_pdsd[657];

    auto g_z_0_x_0_x_xy_0_yz = buffer_1010_pdsd[658];

    auto g_z_0_x_0_x_xy_0_zz = buffer_1010_pdsd[659];

    auto g_z_0_x_0_x_xz_0_xx = buffer_1010_pdsd[660];

    auto g_z_0_x_0_x_xz_0_xy = buffer_1010_pdsd[661];

    auto g_z_0_x_0_x_xz_0_xz = buffer_1010_pdsd[662];

    auto g_z_0_x_0_x_xz_0_yy = buffer_1010_pdsd[663];

    auto g_z_0_x_0_x_xz_0_yz = buffer_1010_pdsd[664];

    auto g_z_0_x_0_x_xz_0_zz = buffer_1010_pdsd[665];

    auto g_z_0_x_0_x_yy_0_xx = buffer_1010_pdsd[666];

    auto g_z_0_x_0_x_yy_0_xy = buffer_1010_pdsd[667];

    auto g_z_0_x_0_x_yy_0_xz = buffer_1010_pdsd[668];

    auto g_z_0_x_0_x_yy_0_yy = buffer_1010_pdsd[669];

    auto g_z_0_x_0_x_yy_0_yz = buffer_1010_pdsd[670];

    auto g_z_0_x_0_x_yy_0_zz = buffer_1010_pdsd[671];

    auto g_z_0_x_0_x_yz_0_xx = buffer_1010_pdsd[672];

    auto g_z_0_x_0_x_yz_0_xy = buffer_1010_pdsd[673];

    auto g_z_0_x_0_x_yz_0_xz = buffer_1010_pdsd[674];

    auto g_z_0_x_0_x_yz_0_yy = buffer_1010_pdsd[675];

    auto g_z_0_x_0_x_yz_0_yz = buffer_1010_pdsd[676];

    auto g_z_0_x_0_x_yz_0_zz = buffer_1010_pdsd[677];

    auto g_z_0_x_0_x_zz_0_xx = buffer_1010_pdsd[678];

    auto g_z_0_x_0_x_zz_0_xy = buffer_1010_pdsd[679];

    auto g_z_0_x_0_x_zz_0_xz = buffer_1010_pdsd[680];

    auto g_z_0_x_0_x_zz_0_yy = buffer_1010_pdsd[681];

    auto g_z_0_x_0_x_zz_0_yz = buffer_1010_pdsd[682];

    auto g_z_0_x_0_x_zz_0_zz = buffer_1010_pdsd[683];

    auto g_z_0_x_0_y_xx_0_xx = buffer_1010_pdsd[684];

    auto g_z_0_x_0_y_xx_0_xy = buffer_1010_pdsd[685];

    auto g_z_0_x_0_y_xx_0_xz = buffer_1010_pdsd[686];

    auto g_z_0_x_0_y_xx_0_yy = buffer_1010_pdsd[687];

    auto g_z_0_x_0_y_xx_0_yz = buffer_1010_pdsd[688];

    auto g_z_0_x_0_y_xx_0_zz = buffer_1010_pdsd[689];

    auto g_z_0_x_0_y_xy_0_xx = buffer_1010_pdsd[690];

    auto g_z_0_x_0_y_xy_0_xy = buffer_1010_pdsd[691];

    auto g_z_0_x_0_y_xy_0_xz = buffer_1010_pdsd[692];

    auto g_z_0_x_0_y_xy_0_yy = buffer_1010_pdsd[693];

    auto g_z_0_x_0_y_xy_0_yz = buffer_1010_pdsd[694];

    auto g_z_0_x_0_y_xy_0_zz = buffer_1010_pdsd[695];

    auto g_z_0_x_0_y_xz_0_xx = buffer_1010_pdsd[696];

    auto g_z_0_x_0_y_xz_0_xy = buffer_1010_pdsd[697];

    auto g_z_0_x_0_y_xz_0_xz = buffer_1010_pdsd[698];

    auto g_z_0_x_0_y_xz_0_yy = buffer_1010_pdsd[699];

    auto g_z_0_x_0_y_xz_0_yz = buffer_1010_pdsd[700];

    auto g_z_0_x_0_y_xz_0_zz = buffer_1010_pdsd[701];

    auto g_z_0_x_0_y_yy_0_xx = buffer_1010_pdsd[702];

    auto g_z_0_x_0_y_yy_0_xy = buffer_1010_pdsd[703];

    auto g_z_0_x_0_y_yy_0_xz = buffer_1010_pdsd[704];

    auto g_z_0_x_0_y_yy_0_yy = buffer_1010_pdsd[705];

    auto g_z_0_x_0_y_yy_0_yz = buffer_1010_pdsd[706];

    auto g_z_0_x_0_y_yy_0_zz = buffer_1010_pdsd[707];

    auto g_z_0_x_0_y_yz_0_xx = buffer_1010_pdsd[708];

    auto g_z_0_x_0_y_yz_0_xy = buffer_1010_pdsd[709];

    auto g_z_0_x_0_y_yz_0_xz = buffer_1010_pdsd[710];

    auto g_z_0_x_0_y_yz_0_yy = buffer_1010_pdsd[711];

    auto g_z_0_x_0_y_yz_0_yz = buffer_1010_pdsd[712];

    auto g_z_0_x_0_y_yz_0_zz = buffer_1010_pdsd[713];

    auto g_z_0_x_0_y_zz_0_xx = buffer_1010_pdsd[714];

    auto g_z_0_x_0_y_zz_0_xy = buffer_1010_pdsd[715];

    auto g_z_0_x_0_y_zz_0_xz = buffer_1010_pdsd[716];

    auto g_z_0_x_0_y_zz_0_yy = buffer_1010_pdsd[717];

    auto g_z_0_x_0_y_zz_0_yz = buffer_1010_pdsd[718];

    auto g_z_0_x_0_y_zz_0_zz = buffer_1010_pdsd[719];

    auto g_z_0_x_0_z_xx_0_xx = buffer_1010_pdsd[720];

    auto g_z_0_x_0_z_xx_0_xy = buffer_1010_pdsd[721];

    auto g_z_0_x_0_z_xx_0_xz = buffer_1010_pdsd[722];

    auto g_z_0_x_0_z_xx_0_yy = buffer_1010_pdsd[723];

    auto g_z_0_x_0_z_xx_0_yz = buffer_1010_pdsd[724];

    auto g_z_0_x_0_z_xx_0_zz = buffer_1010_pdsd[725];

    auto g_z_0_x_0_z_xy_0_xx = buffer_1010_pdsd[726];

    auto g_z_0_x_0_z_xy_0_xy = buffer_1010_pdsd[727];

    auto g_z_0_x_0_z_xy_0_xz = buffer_1010_pdsd[728];

    auto g_z_0_x_0_z_xy_0_yy = buffer_1010_pdsd[729];

    auto g_z_0_x_0_z_xy_0_yz = buffer_1010_pdsd[730];

    auto g_z_0_x_0_z_xy_0_zz = buffer_1010_pdsd[731];

    auto g_z_0_x_0_z_xz_0_xx = buffer_1010_pdsd[732];

    auto g_z_0_x_0_z_xz_0_xy = buffer_1010_pdsd[733];

    auto g_z_0_x_0_z_xz_0_xz = buffer_1010_pdsd[734];

    auto g_z_0_x_0_z_xz_0_yy = buffer_1010_pdsd[735];

    auto g_z_0_x_0_z_xz_0_yz = buffer_1010_pdsd[736];

    auto g_z_0_x_0_z_xz_0_zz = buffer_1010_pdsd[737];

    auto g_z_0_x_0_z_yy_0_xx = buffer_1010_pdsd[738];

    auto g_z_0_x_0_z_yy_0_xy = buffer_1010_pdsd[739];

    auto g_z_0_x_0_z_yy_0_xz = buffer_1010_pdsd[740];

    auto g_z_0_x_0_z_yy_0_yy = buffer_1010_pdsd[741];

    auto g_z_0_x_0_z_yy_0_yz = buffer_1010_pdsd[742];

    auto g_z_0_x_0_z_yy_0_zz = buffer_1010_pdsd[743];

    auto g_z_0_x_0_z_yz_0_xx = buffer_1010_pdsd[744];

    auto g_z_0_x_0_z_yz_0_xy = buffer_1010_pdsd[745];

    auto g_z_0_x_0_z_yz_0_xz = buffer_1010_pdsd[746];

    auto g_z_0_x_0_z_yz_0_yy = buffer_1010_pdsd[747];

    auto g_z_0_x_0_z_yz_0_yz = buffer_1010_pdsd[748];

    auto g_z_0_x_0_z_yz_0_zz = buffer_1010_pdsd[749];

    auto g_z_0_x_0_z_zz_0_xx = buffer_1010_pdsd[750];

    auto g_z_0_x_0_z_zz_0_xy = buffer_1010_pdsd[751];

    auto g_z_0_x_0_z_zz_0_xz = buffer_1010_pdsd[752];

    auto g_z_0_x_0_z_zz_0_yy = buffer_1010_pdsd[753];

    auto g_z_0_x_0_z_zz_0_yz = buffer_1010_pdsd[754];

    auto g_z_0_x_0_z_zz_0_zz = buffer_1010_pdsd[755];

    auto g_z_0_y_0_x_xx_0_xx = buffer_1010_pdsd[756];

    auto g_z_0_y_0_x_xx_0_xy = buffer_1010_pdsd[757];

    auto g_z_0_y_0_x_xx_0_xz = buffer_1010_pdsd[758];

    auto g_z_0_y_0_x_xx_0_yy = buffer_1010_pdsd[759];

    auto g_z_0_y_0_x_xx_0_yz = buffer_1010_pdsd[760];

    auto g_z_0_y_0_x_xx_0_zz = buffer_1010_pdsd[761];

    auto g_z_0_y_0_x_xy_0_xx = buffer_1010_pdsd[762];

    auto g_z_0_y_0_x_xy_0_xy = buffer_1010_pdsd[763];

    auto g_z_0_y_0_x_xy_0_xz = buffer_1010_pdsd[764];

    auto g_z_0_y_0_x_xy_0_yy = buffer_1010_pdsd[765];

    auto g_z_0_y_0_x_xy_0_yz = buffer_1010_pdsd[766];

    auto g_z_0_y_0_x_xy_0_zz = buffer_1010_pdsd[767];

    auto g_z_0_y_0_x_xz_0_xx = buffer_1010_pdsd[768];

    auto g_z_0_y_0_x_xz_0_xy = buffer_1010_pdsd[769];

    auto g_z_0_y_0_x_xz_0_xz = buffer_1010_pdsd[770];

    auto g_z_0_y_0_x_xz_0_yy = buffer_1010_pdsd[771];

    auto g_z_0_y_0_x_xz_0_yz = buffer_1010_pdsd[772];

    auto g_z_0_y_0_x_xz_0_zz = buffer_1010_pdsd[773];

    auto g_z_0_y_0_x_yy_0_xx = buffer_1010_pdsd[774];

    auto g_z_0_y_0_x_yy_0_xy = buffer_1010_pdsd[775];

    auto g_z_0_y_0_x_yy_0_xz = buffer_1010_pdsd[776];

    auto g_z_0_y_0_x_yy_0_yy = buffer_1010_pdsd[777];

    auto g_z_0_y_0_x_yy_0_yz = buffer_1010_pdsd[778];

    auto g_z_0_y_0_x_yy_0_zz = buffer_1010_pdsd[779];

    auto g_z_0_y_0_x_yz_0_xx = buffer_1010_pdsd[780];

    auto g_z_0_y_0_x_yz_0_xy = buffer_1010_pdsd[781];

    auto g_z_0_y_0_x_yz_0_xz = buffer_1010_pdsd[782];

    auto g_z_0_y_0_x_yz_0_yy = buffer_1010_pdsd[783];

    auto g_z_0_y_0_x_yz_0_yz = buffer_1010_pdsd[784];

    auto g_z_0_y_0_x_yz_0_zz = buffer_1010_pdsd[785];

    auto g_z_0_y_0_x_zz_0_xx = buffer_1010_pdsd[786];

    auto g_z_0_y_0_x_zz_0_xy = buffer_1010_pdsd[787];

    auto g_z_0_y_0_x_zz_0_xz = buffer_1010_pdsd[788];

    auto g_z_0_y_0_x_zz_0_yy = buffer_1010_pdsd[789];

    auto g_z_0_y_0_x_zz_0_yz = buffer_1010_pdsd[790];

    auto g_z_0_y_0_x_zz_0_zz = buffer_1010_pdsd[791];

    auto g_z_0_y_0_y_xx_0_xx = buffer_1010_pdsd[792];

    auto g_z_0_y_0_y_xx_0_xy = buffer_1010_pdsd[793];

    auto g_z_0_y_0_y_xx_0_xz = buffer_1010_pdsd[794];

    auto g_z_0_y_0_y_xx_0_yy = buffer_1010_pdsd[795];

    auto g_z_0_y_0_y_xx_0_yz = buffer_1010_pdsd[796];

    auto g_z_0_y_0_y_xx_0_zz = buffer_1010_pdsd[797];

    auto g_z_0_y_0_y_xy_0_xx = buffer_1010_pdsd[798];

    auto g_z_0_y_0_y_xy_0_xy = buffer_1010_pdsd[799];

    auto g_z_0_y_0_y_xy_0_xz = buffer_1010_pdsd[800];

    auto g_z_0_y_0_y_xy_0_yy = buffer_1010_pdsd[801];

    auto g_z_0_y_0_y_xy_0_yz = buffer_1010_pdsd[802];

    auto g_z_0_y_0_y_xy_0_zz = buffer_1010_pdsd[803];

    auto g_z_0_y_0_y_xz_0_xx = buffer_1010_pdsd[804];

    auto g_z_0_y_0_y_xz_0_xy = buffer_1010_pdsd[805];

    auto g_z_0_y_0_y_xz_0_xz = buffer_1010_pdsd[806];

    auto g_z_0_y_0_y_xz_0_yy = buffer_1010_pdsd[807];

    auto g_z_0_y_0_y_xz_0_yz = buffer_1010_pdsd[808];

    auto g_z_0_y_0_y_xz_0_zz = buffer_1010_pdsd[809];

    auto g_z_0_y_0_y_yy_0_xx = buffer_1010_pdsd[810];

    auto g_z_0_y_0_y_yy_0_xy = buffer_1010_pdsd[811];

    auto g_z_0_y_0_y_yy_0_xz = buffer_1010_pdsd[812];

    auto g_z_0_y_0_y_yy_0_yy = buffer_1010_pdsd[813];

    auto g_z_0_y_0_y_yy_0_yz = buffer_1010_pdsd[814];

    auto g_z_0_y_0_y_yy_0_zz = buffer_1010_pdsd[815];

    auto g_z_0_y_0_y_yz_0_xx = buffer_1010_pdsd[816];

    auto g_z_0_y_0_y_yz_0_xy = buffer_1010_pdsd[817];

    auto g_z_0_y_0_y_yz_0_xz = buffer_1010_pdsd[818];

    auto g_z_0_y_0_y_yz_0_yy = buffer_1010_pdsd[819];

    auto g_z_0_y_0_y_yz_0_yz = buffer_1010_pdsd[820];

    auto g_z_0_y_0_y_yz_0_zz = buffer_1010_pdsd[821];

    auto g_z_0_y_0_y_zz_0_xx = buffer_1010_pdsd[822];

    auto g_z_0_y_0_y_zz_0_xy = buffer_1010_pdsd[823];

    auto g_z_0_y_0_y_zz_0_xz = buffer_1010_pdsd[824];

    auto g_z_0_y_0_y_zz_0_yy = buffer_1010_pdsd[825];

    auto g_z_0_y_0_y_zz_0_yz = buffer_1010_pdsd[826];

    auto g_z_0_y_0_y_zz_0_zz = buffer_1010_pdsd[827];

    auto g_z_0_y_0_z_xx_0_xx = buffer_1010_pdsd[828];

    auto g_z_0_y_0_z_xx_0_xy = buffer_1010_pdsd[829];

    auto g_z_0_y_0_z_xx_0_xz = buffer_1010_pdsd[830];

    auto g_z_0_y_0_z_xx_0_yy = buffer_1010_pdsd[831];

    auto g_z_0_y_0_z_xx_0_yz = buffer_1010_pdsd[832];

    auto g_z_0_y_0_z_xx_0_zz = buffer_1010_pdsd[833];

    auto g_z_0_y_0_z_xy_0_xx = buffer_1010_pdsd[834];

    auto g_z_0_y_0_z_xy_0_xy = buffer_1010_pdsd[835];

    auto g_z_0_y_0_z_xy_0_xz = buffer_1010_pdsd[836];

    auto g_z_0_y_0_z_xy_0_yy = buffer_1010_pdsd[837];

    auto g_z_0_y_0_z_xy_0_yz = buffer_1010_pdsd[838];

    auto g_z_0_y_0_z_xy_0_zz = buffer_1010_pdsd[839];

    auto g_z_0_y_0_z_xz_0_xx = buffer_1010_pdsd[840];

    auto g_z_0_y_0_z_xz_0_xy = buffer_1010_pdsd[841];

    auto g_z_0_y_0_z_xz_0_xz = buffer_1010_pdsd[842];

    auto g_z_0_y_0_z_xz_0_yy = buffer_1010_pdsd[843];

    auto g_z_0_y_0_z_xz_0_yz = buffer_1010_pdsd[844];

    auto g_z_0_y_0_z_xz_0_zz = buffer_1010_pdsd[845];

    auto g_z_0_y_0_z_yy_0_xx = buffer_1010_pdsd[846];

    auto g_z_0_y_0_z_yy_0_xy = buffer_1010_pdsd[847];

    auto g_z_0_y_0_z_yy_0_xz = buffer_1010_pdsd[848];

    auto g_z_0_y_0_z_yy_0_yy = buffer_1010_pdsd[849];

    auto g_z_0_y_0_z_yy_0_yz = buffer_1010_pdsd[850];

    auto g_z_0_y_0_z_yy_0_zz = buffer_1010_pdsd[851];

    auto g_z_0_y_0_z_yz_0_xx = buffer_1010_pdsd[852];

    auto g_z_0_y_0_z_yz_0_xy = buffer_1010_pdsd[853];

    auto g_z_0_y_0_z_yz_0_xz = buffer_1010_pdsd[854];

    auto g_z_0_y_0_z_yz_0_yy = buffer_1010_pdsd[855];

    auto g_z_0_y_0_z_yz_0_yz = buffer_1010_pdsd[856];

    auto g_z_0_y_0_z_yz_0_zz = buffer_1010_pdsd[857];

    auto g_z_0_y_0_z_zz_0_xx = buffer_1010_pdsd[858];

    auto g_z_0_y_0_z_zz_0_xy = buffer_1010_pdsd[859];

    auto g_z_0_y_0_z_zz_0_xz = buffer_1010_pdsd[860];

    auto g_z_0_y_0_z_zz_0_yy = buffer_1010_pdsd[861];

    auto g_z_0_y_0_z_zz_0_yz = buffer_1010_pdsd[862];

    auto g_z_0_y_0_z_zz_0_zz = buffer_1010_pdsd[863];

    auto g_z_0_z_0_x_xx_0_xx = buffer_1010_pdsd[864];

    auto g_z_0_z_0_x_xx_0_xy = buffer_1010_pdsd[865];

    auto g_z_0_z_0_x_xx_0_xz = buffer_1010_pdsd[866];

    auto g_z_0_z_0_x_xx_0_yy = buffer_1010_pdsd[867];

    auto g_z_0_z_0_x_xx_0_yz = buffer_1010_pdsd[868];

    auto g_z_0_z_0_x_xx_0_zz = buffer_1010_pdsd[869];

    auto g_z_0_z_0_x_xy_0_xx = buffer_1010_pdsd[870];

    auto g_z_0_z_0_x_xy_0_xy = buffer_1010_pdsd[871];

    auto g_z_0_z_0_x_xy_0_xz = buffer_1010_pdsd[872];

    auto g_z_0_z_0_x_xy_0_yy = buffer_1010_pdsd[873];

    auto g_z_0_z_0_x_xy_0_yz = buffer_1010_pdsd[874];

    auto g_z_0_z_0_x_xy_0_zz = buffer_1010_pdsd[875];

    auto g_z_0_z_0_x_xz_0_xx = buffer_1010_pdsd[876];

    auto g_z_0_z_0_x_xz_0_xy = buffer_1010_pdsd[877];

    auto g_z_0_z_0_x_xz_0_xz = buffer_1010_pdsd[878];

    auto g_z_0_z_0_x_xz_0_yy = buffer_1010_pdsd[879];

    auto g_z_0_z_0_x_xz_0_yz = buffer_1010_pdsd[880];

    auto g_z_0_z_0_x_xz_0_zz = buffer_1010_pdsd[881];

    auto g_z_0_z_0_x_yy_0_xx = buffer_1010_pdsd[882];

    auto g_z_0_z_0_x_yy_0_xy = buffer_1010_pdsd[883];

    auto g_z_0_z_0_x_yy_0_xz = buffer_1010_pdsd[884];

    auto g_z_0_z_0_x_yy_0_yy = buffer_1010_pdsd[885];

    auto g_z_0_z_0_x_yy_0_yz = buffer_1010_pdsd[886];

    auto g_z_0_z_0_x_yy_0_zz = buffer_1010_pdsd[887];

    auto g_z_0_z_0_x_yz_0_xx = buffer_1010_pdsd[888];

    auto g_z_0_z_0_x_yz_0_xy = buffer_1010_pdsd[889];

    auto g_z_0_z_0_x_yz_0_xz = buffer_1010_pdsd[890];

    auto g_z_0_z_0_x_yz_0_yy = buffer_1010_pdsd[891];

    auto g_z_0_z_0_x_yz_0_yz = buffer_1010_pdsd[892];

    auto g_z_0_z_0_x_yz_0_zz = buffer_1010_pdsd[893];

    auto g_z_0_z_0_x_zz_0_xx = buffer_1010_pdsd[894];

    auto g_z_0_z_0_x_zz_0_xy = buffer_1010_pdsd[895];

    auto g_z_0_z_0_x_zz_0_xz = buffer_1010_pdsd[896];

    auto g_z_0_z_0_x_zz_0_yy = buffer_1010_pdsd[897];

    auto g_z_0_z_0_x_zz_0_yz = buffer_1010_pdsd[898];

    auto g_z_0_z_0_x_zz_0_zz = buffer_1010_pdsd[899];

    auto g_z_0_z_0_y_xx_0_xx = buffer_1010_pdsd[900];

    auto g_z_0_z_0_y_xx_0_xy = buffer_1010_pdsd[901];

    auto g_z_0_z_0_y_xx_0_xz = buffer_1010_pdsd[902];

    auto g_z_0_z_0_y_xx_0_yy = buffer_1010_pdsd[903];

    auto g_z_0_z_0_y_xx_0_yz = buffer_1010_pdsd[904];

    auto g_z_0_z_0_y_xx_0_zz = buffer_1010_pdsd[905];

    auto g_z_0_z_0_y_xy_0_xx = buffer_1010_pdsd[906];

    auto g_z_0_z_0_y_xy_0_xy = buffer_1010_pdsd[907];

    auto g_z_0_z_0_y_xy_0_xz = buffer_1010_pdsd[908];

    auto g_z_0_z_0_y_xy_0_yy = buffer_1010_pdsd[909];

    auto g_z_0_z_0_y_xy_0_yz = buffer_1010_pdsd[910];

    auto g_z_0_z_0_y_xy_0_zz = buffer_1010_pdsd[911];

    auto g_z_0_z_0_y_xz_0_xx = buffer_1010_pdsd[912];

    auto g_z_0_z_0_y_xz_0_xy = buffer_1010_pdsd[913];

    auto g_z_0_z_0_y_xz_0_xz = buffer_1010_pdsd[914];

    auto g_z_0_z_0_y_xz_0_yy = buffer_1010_pdsd[915];

    auto g_z_0_z_0_y_xz_0_yz = buffer_1010_pdsd[916];

    auto g_z_0_z_0_y_xz_0_zz = buffer_1010_pdsd[917];

    auto g_z_0_z_0_y_yy_0_xx = buffer_1010_pdsd[918];

    auto g_z_0_z_0_y_yy_0_xy = buffer_1010_pdsd[919];

    auto g_z_0_z_0_y_yy_0_xz = buffer_1010_pdsd[920];

    auto g_z_0_z_0_y_yy_0_yy = buffer_1010_pdsd[921];

    auto g_z_0_z_0_y_yy_0_yz = buffer_1010_pdsd[922];

    auto g_z_0_z_0_y_yy_0_zz = buffer_1010_pdsd[923];

    auto g_z_0_z_0_y_yz_0_xx = buffer_1010_pdsd[924];

    auto g_z_0_z_0_y_yz_0_xy = buffer_1010_pdsd[925];

    auto g_z_0_z_0_y_yz_0_xz = buffer_1010_pdsd[926];

    auto g_z_0_z_0_y_yz_0_yy = buffer_1010_pdsd[927];

    auto g_z_0_z_0_y_yz_0_yz = buffer_1010_pdsd[928];

    auto g_z_0_z_0_y_yz_0_zz = buffer_1010_pdsd[929];

    auto g_z_0_z_0_y_zz_0_xx = buffer_1010_pdsd[930];

    auto g_z_0_z_0_y_zz_0_xy = buffer_1010_pdsd[931];

    auto g_z_0_z_0_y_zz_0_xz = buffer_1010_pdsd[932];

    auto g_z_0_z_0_y_zz_0_yy = buffer_1010_pdsd[933];

    auto g_z_0_z_0_y_zz_0_yz = buffer_1010_pdsd[934];

    auto g_z_0_z_0_y_zz_0_zz = buffer_1010_pdsd[935];

    auto g_z_0_z_0_z_xx_0_xx = buffer_1010_pdsd[936];

    auto g_z_0_z_0_z_xx_0_xy = buffer_1010_pdsd[937];

    auto g_z_0_z_0_z_xx_0_xz = buffer_1010_pdsd[938];

    auto g_z_0_z_0_z_xx_0_yy = buffer_1010_pdsd[939];

    auto g_z_0_z_0_z_xx_0_yz = buffer_1010_pdsd[940];

    auto g_z_0_z_0_z_xx_0_zz = buffer_1010_pdsd[941];

    auto g_z_0_z_0_z_xy_0_xx = buffer_1010_pdsd[942];

    auto g_z_0_z_0_z_xy_0_xy = buffer_1010_pdsd[943];

    auto g_z_0_z_0_z_xy_0_xz = buffer_1010_pdsd[944];

    auto g_z_0_z_0_z_xy_0_yy = buffer_1010_pdsd[945];

    auto g_z_0_z_0_z_xy_0_yz = buffer_1010_pdsd[946];

    auto g_z_0_z_0_z_xy_0_zz = buffer_1010_pdsd[947];

    auto g_z_0_z_0_z_xz_0_xx = buffer_1010_pdsd[948];

    auto g_z_0_z_0_z_xz_0_xy = buffer_1010_pdsd[949];

    auto g_z_0_z_0_z_xz_0_xz = buffer_1010_pdsd[950];

    auto g_z_0_z_0_z_xz_0_yy = buffer_1010_pdsd[951];

    auto g_z_0_z_0_z_xz_0_yz = buffer_1010_pdsd[952];

    auto g_z_0_z_0_z_xz_0_zz = buffer_1010_pdsd[953];

    auto g_z_0_z_0_z_yy_0_xx = buffer_1010_pdsd[954];

    auto g_z_0_z_0_z_yy_0_xy = buffer_1010_pdsd[955];

    auto g_z_0_z_0_z_yy_0_xz = buffer_1010_pdsd[956];

    auto g_z_0_z_0_z_yy_0_yy = buffer_1010_pdsd[957];

    auto g_z_0_z_0_z_yy_0_yz = buffer_1010_pdsd[958];

    auto g_z_0_z_0_z_yy_0_zz = buffer_1010_pdsd[959];

    auto g_z_0_z_0_z_yz_0_xx = buffer_1010_pdsd[960];

    auto g_z_0_z_0_z_yz_0_xy = buffer_1010_pdsd[961];

    auto g_z_0_z_0_z_yz_0_xz = buffer_1010_pdsd[962];

    auto g_z_0_z_0_z_yz_0_yy = buffer_1010_pdsd[963];

    auto g_z_0_z_0_z_yz_0_yz = buffer_1010_pdsd[964];

    auto g_z_0_z_0_z_yz_0_zz = buffer_1010_pdsd[965];

    auto g_z_0_z_0_z_zz_0_xx = buffer_1010_pdsd[966];

    auto g_z_0_z_0_z_zz_0_xy = buffer_1010_pdsd[967];

    auto g_z_0_z_0_z_zz_0_xz = buffer_1010_pdsd[968];

    auto g_z_0_z_0_z_zz_0_yy = buffer_1010_pdsd[969];

    auto g_z_0_z_0_z_zz_0_yz = buffer_1010_pdsd[970];

    auto g_z_0_z_0_z_zz_0_zz = buffer_1010_pdsd[971];

    // integrals block (0-6)

    #pragma omp simd aligned(g_0_xx_x_xx, g_0_xx_x_xy, g_0_xx_x_xz, g_0_xx_x_yy, g_0_xx_x_yz, g_0_xx_x_zz, g_x_0_x_0_x_xx_0_xx, g_x_0_x_0_x_xx_0_xy, g_x_0_x_0_x_xx_0_xz, g_x_0_x_0_x_xx_0_yy, g_x_0_x_0_x_xx_0_yz, g_x_0_x_0_x_xx_0_zz, g_xx_xx_x_xx, g_xx_xx_x_xy, g_xx_xx_x_xz, g_xx_xx_x_yy, g_xx_xx_x_yz, g_xx_xx_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_xx_0_xx[i] = -2.0 * g_0_xx_x_xx[i] * c_exps[i] + 4.0 * g_xx_xx_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xx_0_xy[i] = -2.0 * g_0_xx_x_xy[i] * c_exps[i] + 4.0 * g_xx_xx_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xx_0_xz[i] = -2.0 * g_0_xx_x_xz[i] * c_exps[i] + 4.0 * g_xx_xx_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xx_0_yy[i] = -2.0 * g_0_xx_x_yy[i] * c_exps[i] + 4.0 * g_xx_xx_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xx_0_yz[i] = -2.0 * g_0_xx_x_yz[i] * c_exps[i] + 4.0 * g_xx_xx_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xx_0_zz[i] = -2.0 * g_0_xx_x_zz[i] * c_exps[i] + 4.0 * g_xx_xx_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_0_xy_x_xx, g_0_xy_x_xy, g_0_xy_x_xz, g_0_xy_x_yy, g_0_xy_x_yz, g_0_xy_x_zz, g_x_0_x_0_x_xy_0_xx, g_x_0_x_0_x_xy_0_xy, g_x_0_x_0_x_xy_0_xz, g_x_0_x_0_x_xy_0_yy, g_x_0_x_0_x_xy_0_yz, g_x_0_x_0_x_xy_0_zz, g_xx_xy_x_xx, g_xx_xy_x_xy, g_xx_xy_x_xz, g_xx_xy_x_yy, g_xx_xy_x_yz, g_xx_xy_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_xy_0_xx[i] = -2.0 * g_0_xy_x_xx[i] * c_exps[i] + 4.0 * g_xx_xy_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xy_0_xy[i] = -2.0 * g_0_xy_x_xy[i] * c_exps[i] + 4.0 * g_xx_xy_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xy_0_xz[i] = -2.0 * g_0_xy_x_xz[i] * c_exps[i] + 4.0 * g_xx_xy_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xy_0_yy[i] = -2.0 * g_0_xy_x_yy[i] * c_exps[i] + 4.0 * g_xx_xy_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xy_0_yz[i] = -2.0 * g_0_xy_x_yz[i] * c_exps[i] + 4.0 * g_xx_xy_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xy_0_zz[i] = -2.0 * g_0_xy_x_zz[i] * c_exps[i] + 4.0 * g_xx_xy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_0_xz_x_xx, g_0_xz_x_xy, g_0_xz_x_xz, g_0_xz_x_yy, g_0_xz_x_yz, g_0_xz_x_zz, g_x_0_x_0_x_xz_0_xx, g_x_0_x_0_x_xz_0_xy, g_x_0_x_0_x_xz_0_xz, g_x_0_x_0_x_xz_0_yy, g_x_0_x_0_x_xz_0_yz, g_x_0_x_0_x_xz_0_zz, g_xx_xz_x_xx, g_xx_xz_x_xy, g_xx_xz_x_xz, g_xx_xz_x_yy, g_xx_xz_x_yz, g_xx_xz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_xz_0_xx[i] = -2.0 * g_0_xz_x_xx[i] * c_exps[i] + 4.0 * g_xx_xz_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xz_0_xy[i] = -2.0 * g_0_xz_x_xy[i] * c_exps[i] + 4.0 * g_xx_xz_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xz_0_xz[i] = -2.0 * g_0_xz_x_xz[i] * c_exps[i] + 4.0 * g_xx_xz_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xz_0_yy[i] = -2.0 * g_0_xz_x_yy[i] * c_exps[i] + 4.0 * g_xx_xz_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xz_0_yz[i] = -2.0 * g_0_xz_x_yz[i] * c_exps[i] + 4.0 * g_xx_xz_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xz_0_zz[i] = -2.0 * g_0_xz_x_zz[i] * c_exps[i] + 4.0 * g_xx_xz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_0_yy_x_xx, g_0_yy_x_xy, g_0_yy_x_xz, g_0_yy_x_yy, g_0_yy_x_yz, g_0_yy_x_zz, g_x_0_x_0_x_yy_0_xx, g_x_0_x_0_x_yy_0_xy, g_x_0_x_0_x_yy_0_xz, g_x_0_x_0_x_yy_0_yy, g_x_0_x_0_x_yy_0_yz, g_x_0_x_0_x_yy_0_zz, g_xx_yy_x_xx, g_xx_yy_x_xy, g_xx_yy_x_xz, g_xx_yy_x_yy, g_xx_yy_x_yz, g_xx_yy_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_yy_0_xx[i] = -2.0 * g_0_yy_x_xx[i] * c_exps[i] + 4.0 * g_xx_yy_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_yy_0_xy[i] = -2.0 * g_0_yy_x_xy[i] * c_exps[i] + 4.0 * g_xx_yy_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_yy_0_xz[i] = -2.0 * g_0_yy_x_xz[i] * c_exps[i] + 4.0 * g_xx_yy_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_yy_0_yy[i] = -2.0 * g_0_yy_x_yy[i] * c_exps[i] + 4.0 * g_xx_yy_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_yy_0_yz[i] = -2.0 * g_0_yy_x_yz[i] * c_exps[i] + 4.0 * g_xx_yy_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_yy_0_zz[i] = -2.0 * g_0_yy_x_zz[i] * c_exps[i] + 4.0 * g_xx_yy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_0_yz_x_xx, g_0_yz_x_xy, g_0_yz_x_xz, g_0_yz_x_yy, g_0_yz_x_yz, g_0_yz_x_zz, g_x_0_x_0_x_yz_0_xx, g_x_0_x_0_x_yz_0_xy, g_x_0_x_0_x_yz_0_xz, g_x_0_x_0_x_yz_0_yy, g_x_0_x_0_x_yz_0_yz, g_x_0_x_0_x_yz_0_zz, g_xx_yz_x_xx, g_xx_yz_x_xy, g_xx_yz_x_xz, g_xx_yz_x_yy, g_xx_yz_x_yz, g_xx_yz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_yz_0_xx[i] = -2.0 * g_0_yz_x_xx[i] * c_exps[i] + 4.0 * g_xx_yz_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_yz_0_xy[i] = -2.0 * g_0_yz_x_xy[i] * c_exps[i] + 4.0 * g_xx_yz_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_yz_0_xz[i] = -2.0 * g_0_yz_x_xz[i] * c_exps[i] + 4.0 * g_xx_yz_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_yz_0_yy[i] = -2.0 * g_0_yz_x_yy[i] * c_exps[i] + 4.0 * g_xx_yz_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_yz_0_yz[i] = -2.0 * g_0_yz_x_yz[i] * c_exps[i] + 4.0 * g_xx_yz_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_yz_0_zz[i] = -2.0 * g_0_yz_x_zz[i] * c_exps[i] + 4.0 * g_xx_yz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_0_zz_x_xx, g_0_zz_x_xy, g_0_zz_x_xz, g_0_zz_x_yy, g_0_zz_x_yz, g_0_zz_x_zz, g_x_0_x_0_x_zz_0_xx, g_x_0_x_0_x_zz_0_xy, g_x_0_x_0_x_zz_0_xz, g_x_0_x_0_x_zz_0_yy, g_x_0_x_0_x_zz_0_yz, g_x_0_x_0_x_zz_0_zz, g_xx_zz_x_xx, g_xx_zz_x_xy, g_xx_zz_x_xz, g_xx_zz_x_yy, g_xx_zz_x_yz, g_xx_zz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_zz_0_xx[i] = -2.0 * g_0_zz_x_xx[i] * c_exps[i] + 4.0 * g_xx_zz_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_zz_0_xy[i] = -2.0 * g_0_zz_x_xy[i] * c_exps[i] + 4.0 * g_xx_zz_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_zz_0_xz[i] = -2.0 * g_0_zz_x_xz[i] * c_exps[i] + 4.0 * g_xx_zz_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_zz_0_yy[i] = -2.0 * g_0_zz_x_yy[i] * c_exps[i] + 4.0 * g_xx_zz_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_zz_0_yz[i] = -2.0 * g_0_zz_x_yz[i] * c_exps[i] + 4.0 * g_xx_zz_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_zz_0_zz[i] = -2.0 * g_0_zz_x_zz[i] * c_exps[i] + 4.0 * g_xx_zz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_x_0_x_0_y_xx_0_xx, g_x_0_x_0_y_xx_0_xy, g_x_0_x_0_y_xx_0_xz, g_x_0_x_0_y_xx_0_yy, g_x_0_x_0_y_xx_0_yz, g_x_0_x_0_y_xx_0_zz, g_xy_xx_x_xx, g_xy_xx_x_xy, g_xy_xx_x_xz, g_xy_xx_x_yy, g_xy_xx_x_yz, g_xy_xx_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_xx_0_xx[i] = 4.0 * g_xy_xx_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xx_0_xy[i] = 4.0 * g_xy_xx_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xx_0_xz[i] = 4.0 * g_xy_xx_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xx_0_yy[i] = 4.0 * g_xy_xx_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xx_0_yz[i] = 4.0 * g_xy_xx_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xx_0_zz[i] = 4.0 * g_xy_xx_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_x_0_x_0_y_xy_0_xx, g_x_0_x_0_y_xy_0_xy, g_x_0_x_0_y_xy_0_xz, g_x_0_x_0_y_xy_0_yy, g_x_0_x_0_y_xy_0_yz, g_x_0_x_0_y_xy_0_zz, g_xy_xy_x_xx, g_xy_xy_x_xy, g_xy_xy_x_xz, g_xy_xy_x_yy, g_xy_xy_x_yz, g_xy_xy_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_xy_0_xx[i] = 4.0 * g_xy_xy_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xy_0_xy[i] = 4.0 * g_xy_xy_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xy_0_xz[i] = 4.0 * g_xy_xy_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xy_0_yy[i] = 4.0 * g_xy_xy_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xy_0_yz[i] = 4.0 * g_xy_xy_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xy_0_zz[i] = 4.0 * g_xy_xy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_x_0_x_0_y_xz_0_xx, g_x_0_x_0_y_xz_0_xy, g_x_0_x_0_y_xz_0_xz, g_x_0_x_0_y_xz_0_yy, g_x_0_x_0_y_xz_0_yz, g_x_0_x_0_y_xz_0_zz, g_xy_xz_x_xx, g_xy_xz_x_xy, g_xy_xz_x_xz, g_xy_xz_x_yy, g_xy_xz_x_yz, g_xy_xz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_xz_0_xx[i] = 4.0 * g_xy_xz_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xz_0_xy[i] = 4.0 * g_xy_xz_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xz_0_xz[i] = 4.0 * g_xy_xz_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xz_0_yy[i] = 4.0 * g_xy_xz_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xz_0_yz[i] = 4.0 * g_xy_xz_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xz_0_zz[i] = 4.0 * g_xy_xz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_x_0_x_0_y_yy_0_xx, g_x_0_x_0_y_yy_0_xy, g_x_0_x_0_y_yy_0_xz, g_x_0_x_0_y_yy_0_yy, g_x_0_x_0_y_yy_0_yz, g_x_0_x_0_y_yy_0_zz, g_xy_yy_x_xx, g_xy_yy_x_xy, g_xy_yy_x_xz, g_xy_yy_x_yy, g_xy_yy_x_yz, g_xy_yy_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_yy_0_xx[i] = 4.0 * g_xy_yy_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_yy_0_xy[i] = 4.0 * g_xy_yy_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_yy_0_xz[i] = 4.0 * g_xy_yy_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_yy_0_yy[i] = 4.0 * g_xy_yy_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_yy_0_yz[i] = 4.0 * g_xy_yy_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_yy_0_zz[i] = 4.0 * g_xy_yy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_x_0_x_0_y_yz_0_xx, g_x_0_x_0_y_yz_0_xy, g_x_0_x_0_y_yz_0_xz, g_x_0_x_0_y_yz_0_yy, g_x_0_x_0_y_yz_0_yz, g_x_0_x_0_y_yz_0_zz, g_xy_yz_x_xx, g_xy_yz_x_xy, g_xy_yz_x_xz, g_xy_yz_x_yy, g_xy_yz_x_yz, g_xy_yz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_yz_0_xx[i] = 4.0 * g_xy_yz_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_yz_0_xy[i] = 4.0 * g_xy_yz_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_yz_0_xz[i] = 4.0 * g_xy_yz_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_yz_0_yy[i] = 4.0 * g_xy_yz_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_yz_0_yz[i] = 4.0 * g_xy_yz_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_yz_0_zz[i] = 4.0 * g_xy_yz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_x_0_x_0_y_zz_0_xx, g_x_0_x_0_y_zz_0_xy, g_x_0_x_0_y_zz_0_xz, g_x_0_x_0_y_zz_0_yy, g_x_0_x_0_y_zz_0_yz, g_x_0_x_0_y_zz_0_zz, g_xy_zz_x_xx, g_xy_zz_x_xy, g_xy_zz_x_xz, g_xy_zz_x_yy, g_xy_zz_x_yz, g_xy_zz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_zz_0_xx[i] = 4.0 * g_xy_zz_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_zz_0_xy[i] = 4.0 * g_xy_zz_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_zz_0_xz[i] = 4.0 * g_xy_zz_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_zz_0_yy[i] = 4.0 * g_xy_zz_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_zz_0_yz[i] = 4.0 * g_xy_zz_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_zz_0_zz[i] = 4.0 * g_xy_zz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_x_0_x_0_z_xx_0_xx, g_x_0_x_0_z_xx_0_xy, g_x_0_x_0_z_xx_0_xz, g_x_0_x_0_z_xx_0_yy, g_x_0_x_0_z_xx_0_yz, g_x_0_x_0_z_xx_0_zz, g_xz_xx_x_xx, g_xz_xx_x_xy, g_xz_xx_x_xz, g_xz_xx_x_yy, g_xz_xx_x_yz, g_xz_xx_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_xx_0_xx[i] = 4.0 * g_xz_xx_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xx_0_xy[i] = 4.0 * g_xz_xx_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xx_0_xz[i] = 4.0 * g_xz_xx_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xx_0_yy[i] = 4.0 * g_xz_xx_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xx_0_yz[i] = 4.0 * g_xz_xx_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xx_0_zz[i] = 4.0 * g_xz_xx_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_x_0_x_0_z_xy_0_xx, g_x_0_x_0_z_xy_0_xy, g_x_0_x_0_z_xy_0_xz, g_x_0_x_0_z_xy_0_yy, g_x_0_x_0_z_xy_0_yz, g_x_0_x_0_z_xy_0_zz, g_xz_xy_x_xx, g_xz_xy_x_xy, g_xz_xy_x_xz, g_xz_xy_x_yy, g_xz_xy_x_yz, g_xz_xy_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_xy_0_xx[i] = 4.0 * g_xz_xy_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xy_0_xy[i] = 4.0 * g_xz_xy_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xy_0_xz[i] = 4.0 * g_xz_xy_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xy_0_yy[i] = 4.0 * g_xz_xy_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xy_0_yz[i] = 4.0 * g_xz_xy_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xy_0_zz[i] = 4.0 * g_xz_xy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_x_0_x_0_z_xz_0_xx, g_x_0_x_0_z_xz_0_xy, g_x_0_x_0_z_xz_0_xz, g_x_0_x_0_z_xz_0_yy, g_x_0_x_0_z_xz_0_yz, g_x_0_x_0_z_xz_0_zz, g_xz_xz_x_xx, g_xz_xz_x_xy, g_xz_xz_x_xz, g_xz_xz_x_yy, g_xz_xz_x_yz, g_xz_xz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_xz_0_xx[i] = 4.0 * g_xz_xz_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xz_0_xy[i] = 4.0 * g_xz_xz_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xz_0_xz[i] = 4.0 * g_xz_xz_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xz_0_yy[i] = 4.0 * g_xz_xz_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xz_0_yz[i] = 4.0 * g_xz_xz_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xz_0_zz[i] = 4.0 * g_xz_xz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_x_0_x_0_z_yy_0_xx, g_x_0_x_0_z_yy_0_xy, g_x_0_x_0_z_yy_0_xz, g_x_0_x_0_z_yy_0_yy, g_x_0_x_0_z_yy_0_yz, g_x_0_x_0_z_yy_0_zz, g_xz_yy_x_xx, g_xz_yy_x_xy, g_xz_yy_x_xz, g_xz_yy_x_yy, g_xz_yy_x_yz, g_xz_yy_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_yy_0_xx[i] = 4.0 * g_xz_yy_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_yy_0_xy[i] = 4.0 * g_xz_yy_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_yy_0_xz[i] = 4.0 * g_xz_yy_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_yy_0_yy[i] = 4.0 * g_xz_yy_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_yy_0_yz[i] = 4.0 * g_xz_yy_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_yy_0_zz[i] = 4.0 * g_xz_yy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_x_0_x_0_z_yz_0_xx, g_x_0_x_0_z_yz_0_xy, g_x_0_x_0_z_yz_0_xz, g_x_0_x_0_z_yz_0_yy, g_x_0_x_0_z_yz_0_yz, g_x_0_x_0_z_yz_0_zz, g_xz_yz_x_xx, g_xz_yz_x_xy, g_xz_yz_x_xz, g_xz_yz_x_yy, g_xz_yz_x_yz, g_xz_yz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_yz_0_xx[i] = 4.0 * g_xz_yz_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_yz_0_xy[i] = 4.0 * g_xz_yz_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_yz_0_xz[i] = 4.0 * g_xz_yz_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_yz_0_yy[i] = 4.0 * g_xz_yz_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_yz_0_yz[i] = 4.0 * g_xz_yz_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_yz_0_zz[i] = 4.0 * g_xz_yz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_x_0_x_0_z_zz_0_xx, g_x_0_x_0_z_zz_0_xy, g_x_0_x_0_z_zz_0_xz, g_x_0_x_0_z_zz_0_yy, g_x_0_x_0_z_zz_0_yz, g_x_0_x_0_z_zz_0_zz, g_xz_zz_x_xx, g_xz_zz_x_xy, g_xz_zz_x_xz, g_xz_zz_x_yy, g_xz_zz_x_yz, g_xz_zz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_zz_0_xx[i] = 4.0 * g_xz_zz_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_zz_0_xy[i] = 4.0 * g_xz_zz_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_zz_0_xz[i] = 4.0 * g_xz_zz_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_zz_0_yy[i] = 4.0 * g_xz_zz_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_zz_0_yz[i] = 4.0 * g_xz_zz_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_zz_0_zz[i] = 4.0 * g_xz_zz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_0_xx_y_xx, g_0_xx_y_xy, g_0_xx_y_xz, g_0_xx_y_yy, g_0_xx_y_yz, g_0_xx_y_zz, g_x_0_y_0_x_xx_0_xx, g_x_0_y_0_x_xx_0_xy, g_x_0_y_0_x_xx_0_xz, g_x_0_y_0_x_xx_0_yy, g_x_0_y_0_x_xx_0_yz, g_x_0_y_0_x_xx_0_zz, g_xx_xx_y_xx, g_xx_xx_y_xy, g_xx_xx_y_xz, g_xx_xx_y_yy, g_xx_xx_y_yz, g_xx_xx_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_xx_0_xx[i] = -2.0 * g_0_xx_y_xx[i] * c_exps[i] + 4.0 * g_xx_xx_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xx_0_xy[i] = -2.0 * g_0_xx_y_xy[i] * c_exps[i] + 4.0 * g_xx_xx_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xx_0_xz[i] = -2.0 * g_0_xx_y_xz[i] * c_exps[i] + 4.0 * g_xx_xx_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xx_0_yy[i] = -2.0 * g_0_xx_y_yy[i] * c_exps[i] + 4.0 * g_xx_xx_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xx_0_yz[i] = -2.0 * g_0_xx_y_yz[i] * c_exps[i] + 4.0 * g_xx_xx_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xx_0_zz[i] = -2.0 * g_0_xx_y_zz[i] * c_exps[i] + 4.0 * g_xx_xx_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_0_xy_y_xx, g_0_xy_y_xy, g_0_xy_y_xz, g_0_xy_y_yy, g_0_xy_y_yz, g_0_xy_y_zz, g_x_0_y_0_x_xy_0_xx, g_x_0_y_0_x_xy_0_xy, g_x_0_y_0_x_xy_0_xz, g_x_0_y_0_x_xy_0_yy, g_x_0_y_0_x_xy_0_yz, g_x_0_y_0_x_xy_0_zz, g_xx_xy_y_xx, g_xx_xy_y_xy, g_xx_xy_y_xz, g_xx_xy_y_yy, g_xx_xy_y_yz, g_xx_xy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_xy_0_xx[i] = -2.0 * g_0_xy_y_xx[i] * c_exps[i] + 4.0 * g_xx_xy_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xy_0_xy[i] = -2.0 * g_0_xy_y_xy[i] * c_exps[i] + 4.0 * g_xx_xy_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xy_0_xz[i] = -2.0 * g_0_xy_y_xz[i] * c_exps[i] + 4.0 * g_xx_xy_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xy_0_yy[i] = -2.0 * g_0_xy_y_yy[i] * c_exps[i] + 4.0 * g_xx_xy_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xy_0_yz[i] = -2.0 * g_0_xy_y_yz[i] * c_exps[i] + 4.0 * g_xx_xy_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xy_0_zz[i] = -2.0 * g_0_xy_y_zz[i] * c_exps[i] + 4.0 * g_xx_xy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_0_xz_y_xx, g_0_xz_y_xy, g_0_xz_y_xz, g_0_xz_y_yy, g_0_xz_y_yz, g_0_xz_y_zz, g_x_0_y_0_x_xz_0_xx, g_x_0_y_0_x_xz_0_xy, g_x_0_y_0_x_xz_0_xz, g_x_0_y_0_x_xz_0_yy, g_x_0_y_0_x_xz_0_yz, g_x_0_y_0_x_xz_0_zz, g_xx_xz_y_xx, g_xx_xz_y_xy, g_xx_xz_y_xz, g_xx_xz_y_yy, g_xx_xz_y_yz, g_xx_xz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_xz_0_xx[i] = -2.0 * g_0_xz_y_xx[i] * c_exps[i] + 4.0 * g_xx_xz_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xz_0_xy[i] = -2.0 * g_0_xz_y_xy[i] * c_exps[i] + 4.0 * g_xx_xz_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xz_0_xz[i] = -2.0 * g_0_xz_y_xz[i] * c_exps[i] + 4.0 * g_xx_xz_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xz_0_yy[i] = -2.0 * g_0_xz_y_yy[i] * c_exps[i] + 4.0 * g_xx_xz_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xz_0_yz[i] = -2.0 * g_0_xz_y_yz[i] * c_exps[i] + 4.0 * g_xx_xz_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xz_0_zz[i] = -2.0 * g_0_xz_y_zz[i] * c_exps[i] + 4.0 * g_xx_xz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_0_yy_y_xx, g_0_yy_y_xy, g_0_yy_y_xz, g_0_yy_y_yy, g_0_yy_y_yz, g_0_yy_y_zz, g_x_0_y_0_x_yy_0_xx, g_x_0_y_0_x_yy_0_xy, g_x_0_y_0_x_yy_0_xz, g_x_0_y_0_x_yy_0_yy, g_x_0_y_0_x_yy_0_yz, g_x_0_y_0_x_yy_0_zz, g_xx_yy_y_xx, g_xx_yy_y_xy, g_xx_yy_y_xz, g_xx_yy_y_yy, g_xx_yy_y_yz, g_xx_yy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_yy_0_xx[i] = -2.0 * g_0_yy_y_xx[i] * c_exps[i] + 4.0 * g_xx_yy_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_yy_0_xy[i] = -2.0 * g_0_yy_y_xy[i] * c_exps[i] + 4.0 * g_xx_yy_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_yy_0_xz[i] = -2.0 * g_0_yy_y_xz[i] * c_exps[i] + 4.0 * g_xx_yy_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_yy_0_yy[i] = -2.0 * g_0_yy_y_yy[i] * c_exps[i] + 4.0 * g_xx_yy_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_yy_0_yz[i] = -2.0 * g_0_yy_y_yz[i] * c_exps[i] + 4.0 * g_xx_yy_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_yy_0_zz[i] = -2.0 * g_0_yy_y_zz[i] * c_exps[i] + 4.0 * g_xx_yy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_0_yz_y_xx, g_0_yz_y_xy, g_0_yz_y_xz, g_0_yz_y_yy, g_0_yz_y_yz, g_0_yz_y_zz, g_x_0_y_0_x_yz_0_xx, g_x_0_y_0_x_yz_0_xy, g_x_0_y_0_x_yz_0_xz, g_x_0_y_0_x_yz_0_yy, g_x_0_y_0_x_yz_0_yz, g_x_0_y_0_x_yz_0_zz, g_xx_yz_y_xx, g_xx_yz_y_xy, g_xx_yz_y_xz, g_xx_yz_y_yy, g_xx_yz_y_yz, g_xx_yz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_yz_0_xx[i] = -2.0 * g_0_yz_y_xx[i] * c_exps[i] + 4.0 * g_xx_yz_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_yz_0_xy[i] = -2.0 * g_0_yz_y_xy[i] * c_exps[i] + 4.0 * g_xx_yz_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_yz_0_xz[i] = -2.0 * g_0_yz_y_xz[i] * c_exps[i] + 4.0 * g_xx_yz_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_yz_0_yy[i] = -2.0 * g_0_yz_y_yy[i] * c_exps[i] + 4.0 * g_xx_yz_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_yz_0_yz[i] = -2.0 * g_0_yz_y_yz[i] * c_exps[i] + 4.0 * g_xx_yz_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_yz_0_zz[i] = -2.0 * g_0_yz_y_zz[i] * c_exps[i] + 4.0 * g_xx_yz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_0_zz_y_xx, g_0_zz_y_xy, g_0_zz_y_xz, g_0_zz_y_yy, g_0_zz_y_yz, g_0_zz_y_zz, g_x_0_y_0_x_zz_0_xx, g_x_0_y_0_x_zz_0_xy, g_x_0_y_0_x_zz_0_xz, g_x_0_y_0_x_zz_0_yy, g_x_0_y_0_x_zz_0_yz, g_x_0_y_0_x_zz_0_zz, g_xx_zz_y_xx, g_xx_zz_y_xy, g_xx_zz_y_xz, g_xx_zz_y_yy, g_xx_zz_y_yz, g_xx_zz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_zz_0_xx[i] = -2.0 * g_0_zz_y_xx[i] * c_exps[i] + 4.0 * g_xx_zz_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_zz_0_xy[i] = -2.0 * g_0_zz_y_xy[i] * c_exps[i] + 4.0 * g_xx_zz_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_zz_0_xz[i] = -2.0 * g_0_zz_y_xz[i] * c_exps[i] + 4.0 * g_xx_zz_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_zz_0_yy[i] = -2.0 * g_0_zz_y_yy[i] * c_exps[i] + 4.0 * g_xx_zz_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_zz_0_yz[i] = -2.0 * g_0_zz_y_yz[i] * c_exps[i] + 4.0 * g_xx_zz_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_zz_0_zz[i] = -2.0 * g_0_zz_y_zz[i] * c_exps[i] + 4.0 * g_xx_zz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_x_0_y_0_y_xx_0_xx, g_x_0_y_0_y_xx_0_xy, g_x_0_y_0_y_xx_0_xz, g_x_0_y_0_y_xx_0_yy, g_x_0_y_0_y_xx_0_yz, g_x_0_y_0_y_xx_0_zz, g_xy_xx_y_xx, g_xy_xx_y_xy, g_xy_xx_y_xz, g_xy_xx_y_yy, g_xy_xx_y_yz, g_xy_xx_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_xx_0_xx[i] = 4.0 * g_xy_xx_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xx_0_xy[i] = 4.0 * g_xy_xx_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xx_0_xz[i] = 4.0 * g_xy_xx_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xx_0_yy[i] = 4.0 * g_xy_xx_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xx_0_yz[i] = 4.0 * g_xy_xx_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xx_0_zz[i] = 4.0 * g_xy_xx_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_x_0_y_0_y_xy_0_xx, g_x_0_y_0_y_xy_0_xy, g_x_0_y_0_y_xy_0_xz, g_x_0_y_0_y_xy_0_yy, g_x_0_y_0_y_xy_0_yz, g_x_0_y_0_y_xy_0_zz, g_xy_xy_y_xx, g_xy_xy_y_xy, g_xy_xy_y_xz, g_xy_xy_y_yy, g_xy_xy_y_yz, g_xy_xy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_xy_0_xx[i] = 4.0 * g_xy_xy_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xy_0_xy[i] = 4.0 * g_xy_xy_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xy_0_xz[i] = 4.0 * g_xy_xy_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xy_0_yy[i] = 4.0 * g_xy_xy_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xy_0_yz[i] = 4.0 * g_xy_xy_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xy_0_zz[i] = 4.0 * g_xy_xy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_x_0_y_0_y_xz_0_xx, g_x_0_y_0_y_xz_0_xy, g_x_0_y_0_y_xz_0_xz, g_x_0_y_0_y_xz_0_yy, g_x_0_y_0_y_xz_0_yz, g_x_0_y_0_y_xz_0_zz, g_xy_xz_y_xx, g_xy_xz_y_xy, g_xy_xz_y_xz, g_xy_xz_y_yy, g_xy_xz_y_yz, g_xy_xz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_xz_0_xx[i] = 4.0 * g_xy_xz_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xz_0_xy[i] = 4.0 * g_xy_xz_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xz_0_xz[i] = 4.0 * g_xy_xz_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xz_0_yy[i] = 4.0 * g_xy_xz_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xz_0_yz[i] = 4.0 * g_xy_xz_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xz_0_zz[i] = 4.0 * g_xy_xz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_x_0_y_0_y_yy_0_xx, g_x_0_y_0_y_yy_0_xy, g_x_0_y_0_y_yy_0_xz, g_x_0_y_0_y_yy_0_yy, g_x_0_y_0_y_yy_0_yz, g_x_0_y_0_y_yy_0_zz, g_xy_yy_y_xx, g_xy_yy_y_xy, g_xy_yy_y_xz, g_xy_yy_y_yy, g_xy_yy_y_yz, g_xy_yy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_yy_0_xx[i] = 4.0 * g_xy_yy_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_yy_0_xy[i] = 4.0 * g_xy_yy_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_yy_0_xz[i] = 4.0 * g_xy_yy_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_yy_0_yy[i] = 4.0 * g_xy_yy_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_yy_0_yz[i] = 4.0 * g_xy_yy_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_yy_0_zz[i] = 4.0 * g_xy_yy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_x_0_y_0_y_yz_0_xx, g_x_0_y_0_y_yz_0_xy, g_x_0_y_0_y_yz_0_xz, g_x_0_y_0_y_yz_0_yy, g_x_0_y_0_y_yz_0_yz, g_x_0_y_0_y_yz_0_zz, g_xy_yz_y_xx, g_xy_yz_y_xy, g_xy_yz_y_xz, g_xy_yz_y_yy, g_xy_yz_y_yz, g_xy_yz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_yz_0_xx[i] = 4.0 * g_xy_yz_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_yz_0_xy[i] = 4.0 * g_xy_yz_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_yz_0_xz[i] = 4.0 * g_xy_yz_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_yz_0_yy[i] = 4.0 * g_xy_yz_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_yz_0_yz[i] = 4.0 * g_xy_yz_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_yz_0_zz[i] = 4.0 * g_xy_yz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_x_0_y_0_y_zz_0_xx, g_x_0_y_0_y_zz_0_xy, g_x_0_y_0_y_zz_0_xz, g_x_0_y_0_y_zz_0_yy, g_x_0_y_0_y_zz_0_yz, g_x_0_y_0_y_zz_0_zz, g_xy_zz_y_xx, g_xy_zz_y_xy, g_xy_zz_y_xz, g_xy_zz_y_yy, g_xy_zz_y_yz, g_xy_zz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_zz_0_xx[i] = 4.0 * g_xy_zz_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_zz_0_xy[i] = 4.0 * g_xy_zz_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_zz_0_xz[i] = 4.0 * g_xy_zz_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_zz_0_yy[i] = 4.0 * g_xy_zz_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_zz_0_yz[i] = 4.0 * g_xy_zz_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_zz_0_zz[i] = 4.0 * g_xy_zz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_x_0_y_0_z_xx_0_xx, g_x_0_y_0_z_xx_0_xy, g_x_0_y_0_z_xx_0_xz, g_x_0_y_0_z_xx_0_yy, g_x_0_y_0_z_xx_0_yz, g_x_0_y_0_z_xx_0_zz, g_xz_xx_y_xx, g_xz_xx_y_xy, g_xz_xx_y_xz, g_xz_xx_y_yy, g_xz_xx_y_yz, g_xz_xx_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_xx_0_xx[i] = 4.0 * g_xz_xx_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xx_0_xy[i] = 4.0 * g_xz_xx_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xx_0_xz[i] = 4.0 * g_xz_xx_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xx_0_yy[i] = 4.0 * g_xz_xx_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xx_0_yz[i] = 4.0 * g_xz_xx_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xx_0_zz[i] = 4.0 * g_xz_xx_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_x_0_y_0_z_xy_0_xx, g_x_0_y_0_z_xy_0_xy, g_x_0_y_0_z_xy_0_xz, g_x_0_y_0_z_xy_0_yy, g_x_0_y_0_z_xy_0_yz, g_x_0_y_0_z_xy_0_zz, g_xz_xy_y_xx, g_xz_xy_y_xy, g_xz_xy_y_xz, g_xz_xy_y_yy, g_xz_xy_y_yz, g_xz_xy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_xy_0_xx[i] = 4.0 * g_xz_xy_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xy_0_xy[i] = 4.0 * g_xz_xy_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xy_0_xz[i] = 4.0 * g_xz_xy_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xy_0_yy[i] = 4.0 * g_xz_xy_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xy_0_yz[i] = 4.0 * g_xz_xy_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xy_0_zz[i] = 4.0 * g_xz_xy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_x_0_y_0_z_xz_0_xx, g_x_0_y_0_z_xz_0_xy, g_x_0_y_0_z_xz_0_xz, g_x_0_y_0_z_xz_0_yy, g_x_0_y_0_z_xz_0_yz, g_x_0_y_0_z_xz_0_zz, g_xz_xz_y_xx, g_xz_xz_y_xy, g_xz_xz_y_xz, g_xz_xz_y_yy, g_xz_xz_y_yz, g_xz_xz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_xz_0_xx[i] = 4.0 * g_xz_xz_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xz_0_xy[i] = 4.0 * g_xz_xz_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xz_0_xz[i] = 4.0 * g_xz_xz_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xz_0_yy[i] = 4.0 * g_xz_xz_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xz_0_yz[i] = 4.0 * g_xz_xz_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xz_0_zz[i] = 4.0 * g_xz_xz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_x_0_y_0_z_yy_0_xx, g_x_0_y_0_z_yy_0_xy, g_x_0_y_0_z_yy_0_xz, g_x_0_y_0_z_yy_0_yy, g_x_0_y_0_z_yy_0_yz, g_x_0_y_0_z_yy_0_zz, g_xz_yy_y_xx, g_xz_yy_y_xy, g_xz_yy_y_xz, g_xz_yy_y_yy, g_xz_yy_y_yz, g_xz_yy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_yy_0_xx[i] = 4.0 * g_xz_yy_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_yy_0_xy[i] = 4.0 * g_xz_yy_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_yy_0_xz[i] = 4.0 * g_xz_yy_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_yy_0_yy[i] = 4.0 * g_xz_yy_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_yy_0_yz[i] = 4.0 * g_xz_yy_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_yy_0_zz[i] = 4.0 * g_xz_yy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_x_0_y_0_z_yz_0_xx, g_x_0_y_0_z_yz_0_xy, g_x_0_y_0_z_yz_0_xz, g_x_0_y_0_z_yz_0_yy, g_x_0_y_0_z_yz_0_yz, g_x_0_y_0_z_yz_0_zz, g_xz_yz_y_xx, g_xz_yz_y_xy, g_xz_yz_y_xz, g_xz_yz_y_yy, g_xz_yz_y_yz, g_xz_yz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_yz_0_xx[i] = 4.0 * g_xz_yz_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_yz_0_xy[i] = 4.0 * g_xz_yz_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_yz_0_xz[i] = 4.0 * g_xz_yz_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_yz_0_yy[i] = 4.0 * g_xz_yz_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_yz_0_yz[i] = 4.0 * g_xz_yz_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_yz_0_zz[i] = 4.0 * g_xz_yz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_x_0_y_0_z_zz_0_xx, g_x_0_y_0_z_zz_0_xy, g_x_0_y_0_z_zz_0_xz, g_x_0_y_0_z_zz_0_yy, g_x_0_y_0_z_zz_0_yz, g_x_0_y_0_z_zz_0_zz, g_xz_zz_y_xx, g_xz_zz_y_xy, g_xz_zz_y_xz, g_xz_zz_y_yy, g_xz_zz_y_yz, g_xz_zz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_zz_0_xx[i] = 4.0 * g_xz_zz_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_zz_0_xy[i] = 4.0 * g_xz_zz_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_zz_0_xz[i] = 4.0 * g_xz_zz_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_zz_0_yy[i] = 4.0 * g_xz_zz_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_zz_0_yz[i] = 4.0 * g_xz_zz_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_zz_0_zz[i] = 4.0 * g_xz_zz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_0_xx_z_xx, g_0_xx_z_xy, g_0_xx_z_xz, g_0_xx_z_yy, g_0_xx_z_yz, g_0_xx_z_zz, g_x_0_z_0_x_xx_0_xx, g_x_0_z_0_x_xx_0_xy, g_x_0_z_0_x_xx_0_xz, g_x_0_z_0_x_xx_0_yy, g_x_0_z_0_x_xx_0_yz, g_x_0_z_0_x_xx_0_zz, g_xx_xx_z_xx, g_xx_xx_z_xy, g_xx_xx_z_xz, g_xx_xx_z_yy, g_xx_xx_z_yz, g_xx_xx_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_xx_0_xx[i] = -2.0 * g_0_xx_z_xx[i] * c_exps[i] + 4.0 * g_xx_xx_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xx_0_xy[i] = -2.0 * g_0_xx_z_xy[i] * c_exps[i] + 4.0 * g_xx_xx_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xx_0_xz[i] = -2.0 * g_0_xx_z_xz[i] * c_exps[i] + 4.0 * g_xx_xx_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xx_0_yy[i] = -2.0 * g_0_xx_z_yy[i] * c_exps[i] + 4.0 * g_xx_xx_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xx_0_yz[i] = -2.0 * g_0_xx_z_yz[i] * c_exps[i] + 4.0 * g_xx_xx_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xx_0_zz[i] = -2.0 * g_0_xx_z_zz[i] * c_exps[i] + 4.0 * g_xx_xx_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_0_xy_z_xx, g_0_xy_z_xy, g_0_xy_z_xz, g_0_xy_z_yy, g_0_xy_z_yz, g_0_xy_z_zz, g_x_0_z_0_x_xy_0_xx, g_x_0_z_0_x_xy_0_xy, g_x_0_z_0_x_xy_0_xz, g_x_0_z_0_x_xy_0_yy, g_x_0_z_0_x_xy_0_yz, g_x_0_z_0_x_xy_0_zz, g_xx_xy_z_xx, g_xx_xy_z_xy, g_xx_xy_z_xz, g_xx_xy_z_yy, g_xx_xy_z_yz, g_xx_xy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_xy_0_xx[i] = -2.0 * g_0_xy_z_xx[i] * c_exps[i] + 4.0 * g_xx_xy_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xy_0_xy[i] = -2.0 * g_0_xy_z_xy[i] * c_exps[i] + 4.0 * g_xx_xy_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xy_0_xz[i] = -2.0 * g_0_xy_z_xz[i] * c_exps[i] + 4.0 * g_xx_xy_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xy_0_yy[i] = -2.0 * g_0_xy_z_yy[i] * c_exps[i] + 4.0 * g_xx_xy_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xy_0_yz[i] = -2.0 * g_0_xy_z_yz[i] * c_exps[i] + 4.0 * g_xx_xy_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xy_0_zz[i] = -2.0 * g_0_xy_z_zz[i] * c_exps[i] + 4.0 * g_xx_xy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_0_xz_z_xx, g_0_xz_z_xy, g_0_xz_z_xz, g_0_xz_z_yy, g_0_xz_z_yz, g_0_xz_z_zz, g_x_0_z_0_x_xz_0_xx, g_x_0_z_0_x_xz_0_xy, g_x_0_z_0_x_xz_0_xz, g_x_0_z_0_x_xz_0_yy, g_x_0_z_0_x_xz_0_yz, g_x_0_z_0_x_xz_0_zz, g_xx_xz_z_xx, g_xx_xz_z_xy, g_xx_xz_z_xz, g_xx_xz_z_yy, g_xx_xz_z_yz, g_xx_xz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_xz_0_xx[i] = -2.0 * g_0_xz_z_xx[i] * c_exps[i] + 4.0 * g_xx_xz_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xz_0_xy[i] = -2.0 * g_0_xz_z_xy[i] * c_exps[i] + 4.0 * g_xx_xz_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xz_0_xz[i] = -2.0 * g_0_xz_z_xz[i] * c_exps[i] + 4.0 * g_xx_xz_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xz_0_yy[i] = -2.0 * g_0_xz_z_yy[i] * c_exps[i] + 4.0 * g_xx_xz_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xz_0_yz[i] = -2.0 * g_0_xz_z_yz[i] * c_exps[i] + 4.0 * g_xx_xz_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xz_0_zz[i] = -2.0 * g_0_xz_z_zz[i] * c_exps[i] + 4.0 * g_xx_xz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_0_yy_z_xx, g_0_yy_z_xy, g_0_yy_z_xz, g_0_yy_z_yy, g_0_yy_z_yz, g_0_yy_z_zz, g_x_0_z_0_x_yy_0_xx, g_x_0_z_0_x_yy_0_xy, g_x_0_z_0_x_yy_0_xz, g_x_0_z_0_x_yy_0_yy, g_x_0_z_0_x_yy_0_yz, g_x_0_z_0_x_yy_0_zz, g_xx_yy_z_xx, g_xx_yy_z_xy, g_xx_yy_z_xz, g_xx_yy_z_yy, g_xx_yy_z_yz, g_xx_yy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_yy_0_xx[i] = -2.0 * g_0_yy_z_xx[i] * c_exps[i] + 4.0 * g_xx_yy_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_yy_0_xy[i] = -2.0 * g_0_yy_z_xy[i] * c_exps[i] + 4.0 * g_xx_yy_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_yy_0_xz[i] = -2.0 * g_0_yy_z_xz[i] * c_exps[i] + 4.0 * g_xx_yy_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_yy_0_yy[i] = -2.0 * g_0_yy_z_yy[i] * c_exps[i] + 4.0 * g_xx_yy_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_yy_0_yz[i] = -2.0 * g_0_yy_z_yz[i] * c_exps[i] + 4.0 * g_xx_yy_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_yy_0_zz[i] = -2.0 * g_0_yy_z_zz[i] * c_exps[i] + 4.0 * g_xx_yy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_0_yz_z_xx, g_0_yz_z_xy, g_0_yz_z_xz, g_0_yz_z_yy, g_0_yz_z_yz, g_0_yz_z_zz, g_x_0_z_0_x_yz_0_xx, g_x_0_z_0_x_yz_0_xy, g_x_0_z_0_x_yz_0_xz, g_x_0_z_0_x_yz_0_yy, g_x_0_z_0_x_yz_0_yz, g_x_0_z_0_x_yz_0_zz, g_xx_yz_z_xx, g_xx_yz_z_xy, g_xx_yz_z_xz, g_xx_yz_z_yy, g_xx_yz_z_yz, g_xx_yz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_yz_0_xx[i] = -2.0 * g_0_yz_z_xx[i] * c_exps[i] + 4.0 * g_xx_yz_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_yz_0_xy[i] = -2.0 * g_0_yz_z_xy[i] * c_exps[i] + 4.0 * g_xx_yz_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_yz_0_xz[i] = -2.0 * g_0_yz_z_xz[i] * c_exps[i] + 4.0 * g_xx_yz_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_yz_0_yy[i] = -2.0 * g_0_yz_z_yy[i] * c_exps[i] + 4.0 * g_xx_yz_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_yz_0_yz[i] = -2.0 * g_0_yz_z_yz[i] * c_exps[i] + 4.0 * g_xx_yz_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_yz_0_zz[i] = -2.0 * g_0_yz_z_zz[i] * c_exps[i] + 4.0 * g_xx_yz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_0_zz_z_xx, g_0_zz_z_xy, g_0_zz_z_xz, g_0_zz_z_yy, g_0_zz_z_yz, g_0_zz_z_zz, g_x_0_z_0_x_zz_0_xx, g_x_0_z_0_x_zz_0_xy, g_x_0_z_0_x_zz_0_xz, g_x_0_z_0_x_zz_0_yy, g_x_0_z_0_x_zz_0_yz, g_x_0_z_0_x_zz_0_zz, g_xx_zz_z_xx, g_xx_zz_z_xy, g_xx_zz_z_xz, g_xx_zz_z_yy, g_xx_zz_z_yz, g_xx_zz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_zz_0_xx[i] = -2.0 * g_0_zz_z_xx[i] * c_exps[i] + 4.0 * g_xx_zz_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_zz_0_xy[i] = -2.0 * g_0_zz_z_xy[i] * c_exps[i] + 4.0 * g_xx_zz_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_zz_0_xz[i] = -2.0 * g_0_zz_z_xz[i] * c_exps[i] + 4.0 * g_xx_zz_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_zz_0_yy[i] = -2.0 * g_0_zz_z_yy[i] * c_exps[i] + 4.0 * g_xx_zz_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_zz_0_yz[i] = -2.0 * g_0_zz_z_yz[i] * c_exps[i] + 4.0 * g_xx_zz_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_zz_0_zz[i] = -2.0 * g_0_zz_z_zz[i] * c_exps[i] + 4.0 * g_xx_zz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_x_0_z_0_y_xx_0_xx, g_x_0_z_0_y_xx_0_xy, g_x_0_z_0_y_xx_0_xz, g_x_0_z_0_y_xx_0_yy, g_x_0_z_0_y_xx_0_yz, g_x_0_z_0_y_xx_0_zz, g_xy_xx_z_xx, g_xy_xx_z_xy, g_xy_xx_z_xz, g_xy_xx_z_yy, g_xy_xx_z_yz, g_xy_xx_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_xx_0_xx[i] = 4.0 * g_xy_xx_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xx_0_xy[i] = 4.0 * g_xy_xx_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xx_0_xz[i] = 4.0 * g_xy_xx_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xx_0_yy[i] = 4.0 * g_xy_xx_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xx_0_yz[i] = 4.0 * g_xy_xx_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xx_0_zz[i] = 4.0 * g_xy_xx_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_x_0_z_0_y_xy_0_xx, g_x_0_z_0_y_xy_0_xy, g_x_0_z_0_y_xy_0_xz, g_x_0_z_0_y_xy_0_yy, g_x_0_z_0_y_xy_0_yz, g_x_0_z_0_y_xy_0_zz, g_xy_xy_z_xx, g_xy_xy_z_xy, g_xy_xy_z_xz, g_xy_xy_z_yy, g_xy_xy_z_yz, g_xy_xy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_xy_0_xx[i] = 4.0 * g_xy_xy_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xy_0_xy[i] = 4.0 * g_xy_xy_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xy_0_xz[i] = 4.0 * g_xy_xy_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xy_0_yy[i] = 4.0 * g_xy_xy_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xy_0_yz[i] = 4.0 * g_xy_xy_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xy_0_zz[i] = 4.0 * g_xy_xy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_x_0_z_0_y_xz_0_xx, g_x_0_z_0_y_xz_0_xy, g_x_0_z_0_y_xz_0_xz, g_x_0_z_0_y_xz_0_yy, g_x_0_z_0_y_xz_0_yz, g_x_0_z_0_y_xz_0_zz, g_xy_xz_z_xx, g_xy_xz_z_xy, g_xy_xz_z_xz, g_xy_xz_z_yy, g_xy_xz_z_yz, g_xy_xz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_xz_0_xx[i] = 4.0 * g_xy_xz_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xz_0_xy[i] = 4.0 * g_xy_xz_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xz_0_xz[i] = 4.0 * g_xy_xz_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xz_0_yy[i] = 4.0 * g_xy_xz_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xz_0_yz[i] = 4.0 * g_xy_xz_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xz_0_zz[i] = 4.0 * g_xy_xz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_x_0_z_0_y_yy_0_xx, g_x_0_z_0_y_yy_0_xy, g_x_0_z_0_y_yy_0_xz, g_x_0_z_0_y_yy_0_yy, g_x_0_z_0_y_yy_0_yz, g_x_0_z_0_y_yy_0_zz, g_xy_yy_z_xx, g_xy_yy_z_xy, g_xy_yy_z_xz, g_xy_yy_z_yy, g_xy_yy_z_yz, g_xy_yy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_yy_0_xx[i] = 4.0 * g_xy_yy_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_yy_0_xy[i] = 4.0 * g_xy_yy_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_yy_0_xz[i] = 4.0 * g_xy_yy_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_yy_0_yy[i] = 4.0 * g_xy_yy_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_yy_0_yz[i] = 4.0 * g_xy_yy_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_yy_0_zz[i] = 4.0 * g_xy_yy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_x_0_z_0_y_yz_0_xx, g_x_0_z_0_y_yz_0_xy, g_x_0_z_0_y_yz_0_xz, g_x_0_z_0_y_yz_0_yy, g_x_0_z_0_y_yz_0_yz, g_x_0_z_0_y_yz_0_zz, g_xy_yz_z_xx, g_xy_yz_z_xy, g_xy_yz_z_xz, g_xy_yz_z_yy, g_xy_yz_z_yz, g_xy_yz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_yz_0_xx[i] = 4.0 * g_xy_yz_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_yz_0_xy[i] = 4.0 * g_xy_yz_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_yz_0_xz[i] = 4.0 * g_xy_yz_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_yz_0_yy[i] = 4.0 * g_xy_yz_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_yz_0_yz[i] = 4.0 * g_xy_yz_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_yz_0_zz[i] = 4.0 * g_xy_yz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_x_0_z_0_y_zz_0_xx, g_x_0_z_0_y_zz_0_xy, g_x_0_z_0_y_zz_0_xz, g_x_0_z_0_y_zz_0_yy, g_x_0_z_0_y_zz_0_yz, g_x_0_z_0_y_zz_0_zz, g_xy_zz_z_xx, g_xy_zz_z_xy, g_xy_zz_z_xz, g_xy_zz_z_yy, g_xy_zz_z_yz, g_xy_zz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_zz_0_xx[i] = 4.0 * g_xy_zz_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_zz_0_xy[i] = 4.0 * g_xy_zz_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_zz_0_xz[i] = 4.0 * g_xy_zz_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_zz_0_yy[i] = 4.0 * g_xy_zz_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_zz_0_yz[i] = 4.0 * g_xy_zz_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_zz_0_zz[i] = 4.0 * g_xy_zz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_x_0_z_0_z_xx_0_xx, g_x_0_z_0_z_xx_0_xy, g_x_0_z_0_z_xx_0_xz, g_x_0_z_0_z_xx_0_yy, g_x_0_z_0_z_xx_0_yz, g_x_0_z_0_z_xx_0_zz, g_xz_xx_z_xx, g_xz_xx_z_xy, g_xz_xx_z_xz, g_xz_xx_z_yy, g_xz_xx_z_yz, g_xz_xx_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_xx_0_xx[i] = 4.0 * g_xz_xx_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xx_0_xy[i] = 4.0 * g_xz_xx_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xx_0_xz[i] = 4.0 * g_xz_xx_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xx_0_yy[i] = 4.0 * g_xz_xx_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xx_0_yz[i] = 4.0 * g_xz_xx_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xx_0_zz[i] = 4.0 * g_xz_xx_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_x_0_z_0_z_xy_0_xx, g_x_0_z_0_z_xy_0_xy, g_x_0_z_0_z_xy_0_xz, g_x_0_z_0_z_xy_0_yy, g_x_0_z_0_z_xy_0_yz, g_x_0_z_0_z_xy_0_zz, g_xz_xy_z_xx, g_xz_xy_z_xy, g_xz_xy_z_xz, g_xz_xy_z_yy, g_xz_xy_z_yz, g_xz_xy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_xy_0_xx[i] = 4.0 * g_xz_xy_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xy_0_xy[i] = 4.0 * g_xz_xy_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xy_0_xz[i] = 4.0 * g_xz_xy_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xy_0_yy[i] = 4.0 * g_xz_xy_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xy_0_yz[i] = 4.0 * g_xz_xy_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xy_0_zz[i] = 4.0 * g_xz_xy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_x_0_z_0_z_xz_0_xx, g_x_0_z_0_z_xz_0_xy, g_x_0_z_0_z_xz_0_xz, g_x_0_z_0_z_xz_0_yy, g_x_0_z_0_z_xz_0_yz, g_x_0_z_0_z_xz_0_zz, g_xz_xz_z_xx, g_xz_xz_z_xy, g_xz_xz_z_xz, g_xz_xz_z_yy, g_xz_xz_z_yz, g_xz_xz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_xz_0_xx[i] = 4.0 * g_xz_xz_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xz_0_xy[i] = 4.0 * g_xz_xz_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xz_0_xz[i] = 4.0 * g_xz_xz_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xz_0_yy[i] = 4.0 * g_xz_xz_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xz_0_yz[i] = 4.0 * g_xz_xz_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xz_0_zz[i] = 4.0 * g_xz_xz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_x_0_z_0_z_yy_0_xx, g_x_0_z_0_z_yy_0_xy, g_x_0_z_0_z_yy_0_xz, g_x_0_z_0_z_yy_0_yy, g_x_0_z_0_z_yy_0_yz, g_x_0_z_0_z_yy_0_zz, g_xz_yy_z_xx, g_xz_yy_z_xy, g_xz_yy_z_xz, g_xz_yy_z_yy, g_xz_yy_z_yz, g_xz_yy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_yy_0_xx[i] = 4.0 * g_xz_yy_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_yy_0_xy[i] = 4.0 * g_xz_yy_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_yy_0_xz[i] = 4.0 * g_xz_yy_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_yy_0_yy[i] = 4.0 * g_xz_yy_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_yy_0_yz[i] = 4.0 * g_xz_yy_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_yy_0_zz[i] = 4.0 * g_xz_yy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_x_0_z_0_z_yz_0_xx, g_x_0_z_0_z_yz_0_xy, g_x_0_z_0_z_yz_0_xz, g_x_0_z_0_z_yz_0_yy, g_x_0_z_0_z_yz_0_yz, g_x_0_z_0_z_yz_0_zz, g_xz_yz_z_xx, g_xz_yz_z_xy, g_xz_yz_z_xz, g_xz_yz_z_yy, g_xz_yz_z_yz, g_xz_yz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_yz_0_xx[i] = 4.0 * g_xz_yz_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_yz_0_xy[i] = 4.0 * g_xz_yz_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_yz_0_xz[i] = 4.0 * g_xz_yz_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_yz_0_yy[i] = 4.0 * g_xz_yz_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_yz_0_yz[i] = 4.0 * g_xz_yz_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_yz_0_zz[i] = 4.0 * g_xz_yz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_x_0_z_0_z_zz_0_xx, g_x_0_z_0_z_zz_0_xy, g_x_0_z_0_z_zz_0_xz, g_x_0_z_0_z_zz_0_yy, g_x_0_z_0_z_zz_0_yz, g_x_0_z_0_z_zz_0_zz, g_xz_zz_z_xx, g_xz_zz_z_xy, g_xz_zz_z_xz, g_xz_zz_z_yy, g_xz_zz_z_yz, g_xz_zz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_zz_0_xx[i] = 4.0 * g_xz_zz_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_zz_0_xy[i] = 4.0 * g_xz_zz_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_zz_0_xz[i] = 4.0 * g_xz_zz_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_zz_0_yy[i] = 4.0 * g_xz_zz_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_zz_0_yz[i] = 4.0 * g_xz_zz_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_zz_0_zz[i] = 4.0 * g_xz_zz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (324-330)

    #pragma omp simd aligned(g_xy_xx_x_xx, g_xy_xx_x_xy, g_xy_xx_x_xz, g_xy_xx_x_yy, g_xy_xx_x_yz, g_xy_xx_x_zz, g_y_0_x_0_x_xx_0_xx, g_y_0_x_0_x_xx_0_xy, g_y_0_x_0_x_xx_0_xz, g_y_0_x_0_x_xx_0_yy, g_y_0_x_0_x_xx_0_yz, g_y_0_x_0_x_xx_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_xx_0_xx[i] = 4.0 * g_xy_xx_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xx_0_xy[i] = 4.0 * g_xy_xx_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xx_0_xz[i] = 4.0 * g_xy_xx_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xx_0_yy[i] = 4.0 * g_xy_xx_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xx_0_yz[i] = 4.0 * g_xy_xx_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xx_0_zz[i] = 4.0 * g_xy_xx_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (330-336)

    #pragma omp simd aligned(g_xy_xy_x_xx, g_xy_xy_x_xy, g_xy_xy_x_xz, g_xy_xy_x_yy, g_xy_xy_x_yz, g_xy_xy_x_zz, g_y_0_x_0_x_xy_0_xx, g_y_0_x_0_x_xy_0_xy, g_y_0_x_0_x_xy_0_xz, g_y_0_x_0_x_xy_0_yy, g_y_0_x_0_x_xy_0_yz, g_y_0_x_0_x_xy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_xy_0_xx[i] = 4.0 * g_xy_xy_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xy_0_xy[i] = 4.0 * g_xy_xy_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xy_0_xz[i] = 4.0 * g_xy_xy_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xy_0_yy[i] = 4.0 * g_xy_xy_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xy_0_yz[i] = 4.0 * g_xy_xy_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xy_0_zz[i] = 4.0 * g_xy_xy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (336-342)

    #pragma omp simd aligned(g_xy_xz_x_xx, g_xy_xz_x_xy, g_xy_xz_x_xz, g_xy_xz_x_yy, g_xy_xz_x_yz, g_xy_xz_x_zz, g_y_0_x_0_x_xz_0_xx, g_y_0_x_0_x_xz_0_xy, g_y_0_x_0_x_xz_0_xz, g_y_0_x_0_x_xz_0_yy, g_y_0_x_0_x_xz_0_yz, g_y_0_x_0_x_xz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_xz_0_xx[i] = 4.0 * g_xy_xz_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xz_0_xy[i] = 4.0 * g_xy_xz_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xz_0_xz[i] = 4.0 * g_xy_xz_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xz_0_yy[i] = 4.0 * g_xy_xz_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xz_0_yz[i] = 4.0 * g_xy_xz_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xz_0_zz[i] = 4.0 * g_xy_xz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (342-348)

    #pragma omp simd aligned(g_xy_yy_x_xx, g_xy_yy_x_xy, g_xy_yy_x_xz, g_xy_yy_x_yy, g_xy_yy_x_yz, g_xy_yy_x_zz, g_y_0_x_0_x_yy_0_xx, g_y_0_x_0_x_yy_0_xy, g_y_0_x_0_x_yy_0_xz, g_y_0_x_0_x_yy_0_yy, g_y_0_x_0_x_yy_0_yz, g_y_0_x_0_x_yy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_yy_0_xx[i] = 4.0 * g_xy_yy_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_yy_0_xy[i] = 4.0 * g_xy_yy_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_yy_0_xz[i] = 4.0 * g_xy_yy_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_yy_0_yy[i] = 4.0 * g_xy_yy_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_yy_0_yz[i] = 4.0 * g_xy_yy_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_yy_0_zz[i] = 4.0 * g_xy_yy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (348-354)

    #pragma omp simd aligned(g_xy_yz_x_xx, g_xy_yz_x_xy, g_xy_yz_x_xz, g_xy_yz_x_yy, g_xy_yz_x_yz, g_xy_yz_x_zz, g_y_0_x_0_x_yz_0_xx, g_y_0_x_0_x_yz_0_xy, g_y_0_x_0_x_yz_0_xz, g_y_0_x_0_x_yz_0_yy, g_y_0_x_0_x_yz_0_yz, g_y_0_x_0_x_yz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_yz_0_xx[i] = 4.0 * g_xy_yz_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_yz_0_xy[i] = 4.0 * g_xy_yz_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_yz_0_xz[i] = 4.0 * g_xy_yz_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_yz_0_yy[i] = 4.0 * g_xy_yz_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_yz_0_yz[i] = 4.0 * g_xy_yz_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_yz_0_zz[i] = 4.0 * g_xy_yz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (354-360)

    #pragma omp simd aligned(g_xy_zz_x_xx, g_xy_zz_x_xy, g_xy_zz_x_xz, g_xy_zz_x_yy, g_xy_zz_x_yz, g_xy_zz_x_zz, g_y_0_x_0_x_zz_0_xx, g_y_0_x_0_x_zz_0_xy, g_y_0_x_0_x_zz_0_xz, g_y_0_x_0_x_zz_0_yy, g_y_0_x_0_x_zz_0_yz, g_y_0_x_0_x_zz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_zz_0_xx[i] = 4.0 * g_xy_zz_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_zz_0_xy[i] = 4.0 * g_xy_zz_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_zz_0_xz[i] = 4.0 * g_xy_zz_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_zz_0_yy[i] = 4.0 * g_xy_zz_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_zz_0_yz[i] = 4.0 * g_xy_zz_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_zz_0_zz[i] = 4.0 * g_xy_zz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (360-366)

    #pragma omp simd aligned(g_0_xx_x_xx, g_0_xx_x_xy, g_0_xx_x_xz, g_0_xx_x_yy, g_0_xx_x_yz, g_0_xx_x_zz, g_y_0_x_0_y_xx_0_xx, g_y_0_x_0_y_xx_0_xy, g_y_0_x_0_y_xx_0_xz, g_y_0_x_0_y_xx_0_yy, g_y_0_x_0_y_xx_0_yz, g_y_0_x_0_y_xx_0_zz, g_yy_xx_x_xx, g_yy_xx_x_xy, g_yy_xx_x_xz, g_yy_xx_x_yy, g_yy_xx_x_yz, g_yy_xx_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_xx_0_xx[i] = -2.0 * g_0_xx_x_xx[i] * c_exps[i] + 4.0 * g_yy_xx_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xx_0_xy[i] = -2.0 * g_0_xx_x_xy[i] * c_exps[i] + 4.0 * g_yy_xx_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xx_0_xz[i] = -2.0 * g_0_xx_x_xz[i] * c_exps[i] + 4.0 * g_yy_xx_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xx_0_yy[i] = -2.0 * g_0_xx_x_yy[i] * c_exps[i] + 4.0 * g_yy_xx_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xx_0_yz[i] = -2.0 * g_0_xx_x_yz[i] * c_exps[i] + 4.0 * g_yy_xx_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xx_0_zz[i] = -2.0 * g_0_xx_x_zz[i] * c_exps[i] + 4.0 * g_yy_xx_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (366-372)

    #pragma omp simd aligned(g_0_xy_x_xx, g_0_xy_x_xy, g_0_xy_x_xz, g_0_xy_x_yy, g_0_xy_x_yz, g_0_xy_x_zz, g_y_0_x_0_y_xy_0_xx, g_y_0_x_0_y_xy_0_xy, g_y_0_x_0_y_xy_0_xz, g_y_0_x_0_y_xy_0_yy, g_y_0_x_0_y_xy_0_yz, g_y_0_x_0_y_xy_0_zz, g_yy_xy_x_xx, g_yy_xy_x_xy, g_yy_xy_x_xz, g_yy_xy_x_yy, g_yy_xy_x_yz, g_yy_xy_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_xy_0_xx[i] = -2.0 * g_0_xy_x_xx[i] * c_exps[i] + 4.0 * g_yy_xy_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xy_0_xy[i] = -2.0 * g_0_xy_x_xy[i] * c_exps[i] + 4.0 * g_yy_xy_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xy_0_xz[i] = -2.0 * g_0_xy_x_xz[i] * c_exps[i] + 4.0 * g_yy_xy_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xy_0_yy[i] = -2.0 * g_0_xy_x_yy[i] * c_exps[i] + 4.0 * g_yy_xy_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xy_0_yz[i] = -2.0 * g_0_xy_x_yz[i] * c_exps[i] + 4.0 * g_yy_xy_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xy_0_zz[i] = -2.0 * g_0_xy_x_zz[i] * c_exps[i] + 4.0 * g_yy_xy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (372-378)

    #pragma omp simd aligned(g_0_xz_x_xx, g_0_xz_x_xy, g_0_xz_x_xz, g_0_xz_x_yy, g_0_xz_x_yz, g_0_xz_x_zz, g_y_0_x_0_y_xz_0_xx, g_y_0_x_0_y_xz_0_xy, g_y_0_x_0_y_xz_0_xz, g_y_0_x_0_y_xz_0_yy, g_y_0_x_0_y_xz_0_yz, g_y_0_x_0_y_xz_0_zz, g_yy_xz_x_xx, g_yy_xz_x_xy, g_yy_xz_x_xz, g_yy_xz_x_yy, g_yy_xz_x_yz, g_yy_xz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_xz_0_xx[i] = -2.0 * g_0_xz_x_xx[i] * c_exps[i] + 4.0 * g_yy_xz_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xz_0_xy[i] = -2.0 * g_0_xz_x_xy[i] * c_exps[i] + 4.0 * g_yy_xz_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xz_0_xz[i] = -2.0 * g_0_xz_x_xz[i] * c_exps[i] + 4.0 * g_yy_xz_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xz_0_yy[i] = -2.0 * g_0_xz_x_yy[i] * c_exps[i] + 4.0 * g_yy_xz_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xz_0_yz[i] = -2.0 * g_0_xz_x_yz[i] * c_exps[i] + 4.0 * g_yy_xz_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xz_0_zz[i] = -2.0 * g_0_xz_x_zz[i] * c_exps[i] + 4.0 * g_yy_xz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (378-384)

    #pragma omp simd aligned(g_0_yy_x_xx, g_0_yy_x_xy, g_0_yy_x_xz, g_0_yy_x_yy, g_0_yy_x_yz, g_0_yy_x_zz, g_y_0_x_0_y_yy_0_xx, g_y_0_x_0_y_yy_0_xy, g_y_0_x_0_y_yy_0_xz, g_y_0_x_0_y_yy_0_yy, g_y_0_x_0_y_yy_0_yz, g_y_0_x_0_y_yy_0_zz, g_yy_yy_x_xx, g_yy_yy_x_xy, g_yy_yy_x_xz, g_yy_yy_x_yy, g_yy_yy_x_yz, g_yy_yy_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_yy_0_xx[i] = -2.0 * g_0_yy_x_xx[i] * c_exps[i] + 4.0 * g_yy_yy_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_yy_0_xy[i] = -2.0 * g_0_yy_x_xy[i] * c_exps[i] + 4.0 * g_yy_yy_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_yy_0_xz[i] = -2.0 * g_0_yy_x_xz[i] * c_exps[i] + 4.0 * g_yy_yy_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_yy_0_yy[i] = -2.0 * g_0_yy_x_yy[i] * c_exps[i] + 4.0 * g_yy_yy_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_yy_0_yz[i] = -2.0 * g_0_yy_x_yz[i] * c_exps[i] + 4.0 * g_yy_yy_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_yy_0_zz[i] = -2.0 * g_0_yy_x_zz[i] * c_exps[i] + 4.0 * g_yy_yy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (384-390)

    #pragma omp simd aligned(g_0_yz_x_xx, g_0_yz_x_xy, g_0_yz_x_xz, g_0_yz_x_yy, g_0_yz_x_yz, g_0_yz_x_zz, g_y_0_x_0_y_yz_0_xx, g_y_0_x_0_y_yz_0_xy, g_y_0_x_0_y_yz_0_xz, g_y_0_x_0_y_yz_0_yy, g_y_0_x_0_y_yz_0_yz, g_y_0_x_0_y_yz_0_zz, g_yy_yz_x_xx, g_yy_yz_x_xy, g_yy_yz_x_xz, g_yy_yz_x_yy, g_yy_yz_x_yz, g_yy_yz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_yz_0_xx[i] = -2.0 * g_0_yz_x_xx[i] * c_exps[i] + 4.0 * g_yy_yz_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_yz_0_xy[i] = -2.0 * g_0_yz_x_xy[i] * c_exps[i] + 4.0 * g_yy_yz_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_yz_0_xz[i] = -2.0 * g_0_yz_x_xz[i] * c_exps[i] + 4.0 * g_yy_yz_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_yz_0_yy[i] = -2.0 * g_0_yz_x_yy[i] * c_exps[i] + 4.0 * g_yy_yz_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_yz_0_yz[i] = -2.0 * g_0_yz_x_yz[i] * c_exps[i] + 4.0 * g_yy_yz_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_yz_0_zz[i] = -2.0 * g_0_yz_x_zz[i] * c_exps[i] + 4.0 * g_yy_yz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (390-396)

    #pragma omp simd aligned(g_0_zz_x_xx, g_0_zz_x_xy, g_0_zz_x_xz, g_0_zz_x_yy, g_0_zz_x_yz, g_0_zz_x_zz, g_y_0_x_0_y_zz_0_xx, g_y_0_x_0_y_zz_0_xy, g_y_0_x_0_y_zz_0_xz, g_y_0_x_0_y_zz_0_yy, g_y_0_x_0_y_zz_0_yz, g_y_0_x_0_y_zz_0_zz, g_yy_zz_x_xx, g_yy_zz_x_xy, g_yy_zz_x_xz, g_yy_zz_x_yy, g_yy_zz_x_yz, g_yy_zz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_zz_0_xx[i] = -2.0 * g_0_zz_x_xx[i] * c_exps[i] + 4.0 * g_yy_zz_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_zz_0_xy[i] = -2.0 * g_0_zz_x_xy[i] * c_exps[i] + 4.0 * g_yy_zz_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_zz_0_xz[i] = -2.0 * g_0_zz_x_xz[i] * c_exps[i] + 4.0 * g_yy_zz_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_zz_0_yy[i] = -2.0 * g_0_zz_x_yy[i] * c_exps[i] + 4.0 * g_yy_zz_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_zz_0_yz[i] = -2.0 * g_0_zz_x_yz[i] * c_exps[i] + 4.0 * g_yy_zz_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_zz_0_zz[i] = -2.0 * g_0_zz_x_zz[i] * c_exps[i] + 4.0 * g_yy_zz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (396-402)

    #pragma omp simd aligned(g_y_0_x_0_z_xx_0_xx, g_y_0_x_0_z_xx_0_xy, g_y_0_x_0_z_xx_0_xz, g_y_0_x_0_z_xx_0_yy, g_y_0_x_0_z_xx_0_yz, g_y_0_x_0_z_xx_0_zz, g_yz_xx_x_xx, g_yz_xx_x_xy, g_yz_xx_x_xz, g_yz_xx_x_yy, g_yz_xx_x_yz, g_yz_xx_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_xx_0_xx[i] = 4.0 * g_yz_xx_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xx_0_xy[i] = 4.0 * g_yz_xx_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xx_0_xz[i] = 4.0 * g_yz_xx_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xx_0_yy[i] = 4.0 * g_yz_xx_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xx_0_yz[i] = 4.0 * g_yz_xx_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xx_0_zz[i] = 4.0 * g_yz_xx_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (402-408)

    #pragma omp simd aligned(g_y_0_x_0_z_xy_0_xx, g_y_0_x_0_z_xy_0_xy, g_y_0_x_0_z_xy_0_xz, g_y_0_x_0_z_xy_0_yy, g_y_0_x_0_z_xy_0_yz, g_y_0_x_0_z_xy_0_zz, g_yz_xy_x_xx, g_yz_xy_x_xy, g_yz_xy_x_xz, g_yz_xy_x_yy, g_yz_xy_x_yz, g_yz_xy_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_xy_0_xx[i] = 4.0 * g_yz_xy_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xy_0_xy[i] = 4.0 * g_yz_xy_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xy_0_xz[i] = 4.0 * g_yz_xy_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xy_0_yy[i] = 4.0 * g_yz_xy_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xy_0_yz[i] = 4.0 * g_yz_xy_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xy_0_zz[i] = 4.0 * g_yz_xy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (408-414)

    #pragma omp simd aligned(g_y_0_x_0_z_xz_0_xx, g_y_0_x_0_z_xz_0_xy, g_y_0_x_0_z_xz_0_xz, g_y_0_x_0_z_xz_0_yy, g_y_0_x_0_z_xz_0_yz, g_y_0_x_0_z_xz_0_zz, g_yz_xz_x_xx, g_yz_xz_x_xy, g_yz_xz_x_xz, g_yz_xz_x_yy, g_yz_xz_x_yz, g_yz_xz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_xz_0_xx[i] = 4.0 * g_yz_xz_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xz_0_xy[i] = 4.0 * g_yz_xz_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xz_0_xz[i] = 4.0 * g_yz_xz_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xz_0_yy[i] = 4.0 * g_yz_xz_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xz_0_yz[i] = 4.0 * g_yz_xz_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xz_0_zz[i] = 4.0 * g_yz_xz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (414-420)

    #pragma omp simd aligned(g_y_0_x_0_z_yy_0_xx, g_y_0_x_0_z_yy_0_xy, g_y_0_x_0_z_yy_0_xz, g_y_0_x_0_z_yy_0_yy, g_y_0_x_0_z_yy_0_yz, g_y_0_x_0_z_yy_0_zz, g_yz_yy_x_xx, g_yz_yy_x_xy, g_yz_yy_x_xz, g_yz_yy_x_yy, g_yz_yy_x_yz, g_yz_yy_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_yy_0_xx[i] = 4.0 * g_yz_yy_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_yy_0_xy[i] = 4.0 * g_yz_yy_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_yy_0_xz[i] = 4.0 * g_yz_yy_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_yy_0_yy[i] = 4.0 * g_yz_yy_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_yy_0_yz[i] = 4.0 * g_yz_yy_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_yy_0_zz[i] = 4.0 * g_yz_yy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (420-426)

    #pragma omp simd aligned(g_y_0_x_0_z_yz_0_xx, g_y_0_x_0_z_yz_0_xy, g_y_0_x_0_z_yz_0_xz, g_y_0_x_0_z_yz_0_yy, g_y_0_x_0_z_yz_0_yz, g_y_0_x_0_z_yz_0_zz, g_yz_yz_x_xx, g_yz_yz_x_xy, g_yz_yz_x_xz, g_yz_yz_x_yy, g_yz_yz_x_yz, g_yz_yz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_yz_0_xx[i] = 4.0 * g_yz_yz_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_yz_0_xy[i] = 4.0 * g_yz_yz_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_yz_0_xz[i] = 4.0 * g_yz_yz_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_yz_0_yy[i] = 4.0 * g_yz_yz_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_yz_0_yz[i] = 4.0 * g_yz_yz_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_yz_0_zz[i] = 4.0 * g_yz_yz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (426-432)

    #pragma omp simd aligned(g_y_0_x_0_z_zz_0_xx, g_y_0_x_0_z_zz_0_xy, g_y_0_x_0_z_zz_0_xz, g_y_0_x_0_z_zz_0_yy, g_y_0_x_0_z_zz_0_yz, g_y_0_x_0_z_zz_0_zz, g_yz_zz_x_xx, g_yz_zz_x_xy, g_yz_zz_x_xz, g_yz_zz_x_yy, g_yz_zz_x_yz, g_yz_zz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_zz_0_xx[i] = 4.0 * g_yz_zz_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_zz_0_xy[i] = 4.0 * g_yz_zz_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_zz_0_xz[i] = 4.0 * g_yz_zz_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_zz_0_yy[i] = 4.0 * g_yz_zz_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_zz_0_yz[i] = 4.0 * g_yz_zz_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_zz_0_zz[i] = 4.0 * g_yz_zz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (432-438)

    #pragma omp simd aligned(g_xy_xx_y_xx, g_xy_xx_y_xy, g_xy_xx_y_xz, g_xy_xx_y_yy, g_xy_xx_y_yz, g_xy_xx_y_zz, g_y_0_y_0_x_xx_0_xx, g_y_0_y_0_x_xx_0_xy, g_y_0_y_0_x_xx_0_xz, g_y_0_y_0_x_xx_0_yy, g_y_0_y_0_x_xx_0_yz, g_y_0_y_0_x_xx_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_xx_0_xx[i] = 4.0 * g_xy_xx_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xx_0_xy[i] = 4.0 * g_xy_xx_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xx_0_xz[i] = 4.0 * g_xy_xx_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xx_0_yy[i] = 4.0 * g_xy_xx_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xx_0_yz[i] = 4.0 * g_xy_xx_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xx_0_zz[i] = 4.0 * g_xy_xx_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (438-444)

    #pragma omp simd aligned(g_xy_xy_y_xx, g_xy_xy_y_xy, g_xy_xy_y_xz, g_xy_xy_y_yy, g_xy_xy_y_yz, g_xy_xy_y_zz, g_y_0_y_0_x_xy_0_xx, g_y_0_y_0_x_xy_0_xy, g_y_0_y_0_x_xy_0_xz, g_y_0_y_0_x_xy_0_yy, g_y_0_y_0_x_xy_0_yz, g_y_0_y_0_x_xy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_xy_0_xx[i] = 4.0 * g_xy_xy_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xy_0_xy[i] = 4.0 * g_xy_xy_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xy_0_xz[i] = 4.0 * g_xy_xy_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xy_0_yy[i] = 4.0 * g_xy_xy_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xy_0_yz[i] = 4.0 * g_xy_xy_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xy_0_zz[i] = 4.0 * g_xy_xy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (444-450)

    #pragma omp simd aligned(g_xy_xz_y_xx, g_xy_xz_y_xy, g_xy_xz_y_xz, g_xy_xz_y_yy, g_xy_xz_y_yz, g_xy_xz_y_zz, g_y_0_y_0_x_xz_0_xx, g_y_0_y_0_x_xz_0_xy, g_y_0_y_0_x_xz_0_xz, g_y_0_y_0_x_xz_0_yy, g_y_0_y_0_x_xz_0_yz, g_y_0_y_0_x_xz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_xz_0_xx[i] = 4.0 * g_xy_xz_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xz_0_xy[i] = 4.0 * g_xy_xz_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xz_0_xz[i] = 4.0 * g_xy_xz_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xz_0_yy[i] = 4.0 * g_xy_xz_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xz_0_yz[i] = 4.0 * g_xy_xz_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xz_0_zz[i] = 4.0 * g_xy_xz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (450-456)

    #pragma omp simd aligned(g_xy_yy_y_xx, g_xy_yy_y_xy, g_xy_yy_y_xz, g_xy_yy_y_yy, g_xy_yy_y_yz, g_xy_yy_y_zz, g_y_0_y_0_x_yy_0_xx, g_y_0_y_0_x_yy_0_xy, g_y_0_y_0_x_yy_0_xz, g_y_0_y_0_x_yy_0_yy, g_y_0_y_0_x_yy_0_yz, g_y_0_y_0_x_yy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_yy_0_xx[i] = 4.0 * g_xy_yy_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_yy_0_xy[i] = 4.0 * g_xy_yy_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_yy_0_xz[i] = 4.0 * g_xy_yy_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_yy_0_yy[i] = 4.0 * g_xy_yy_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_yy_0_yz[i] = 4.0 * g_xy_yy_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_yy_0_zz[i] = 4.0 * g_xy_yy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (456-462)

    #pragma omp simd aligned(g_xy_yz_y_xx, g_xy_yz_y_xy, g_xy_yz_y_xz, g_xy_yz_y_yy, g_xy_yz_y_yz, g_xy_yz_y_zz, g_y_0_y_0_x_yz_0_xx, g_y_0_y_0_x_yz_0_xy, g_y_0_y_0_x_yz_0_xz, g_y_0_y_0_x_yz_0_yy, g_y_0_y_0_x_yz_0_yz, g_y_0_y_0_x_yz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_yz_0_xx[i] = 4.0 * g_xy_yz_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_yz_0_xy[i] = 4.0 * g_xy_yz_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_yz_0_xz[i] = 4.0 * g_xy_yz_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_yz_0_yy[i] = 4.0 * g_xy_yz_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_yz_0_yz[i] = 4.0 * g_xy_yz_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_yz_0_zz[i] = 4.0 * g_xy_yz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (462-468)

    #pragma omp simd aligned(g_xy_zz_y_xx, g_xy_zz_y_xy, g_xy_zz_y_xz, g_xy_zz_y_yy, g_xy_zz_y_yz, g_xy_zz_y_zz, g_y_0_y_0_x_zz_0_xx, g_y_0_y_0_x_zz_0_xy, g_y_0_y_0_x_zz_0_xz, g_y_0_y_0_x_zz_0_yy, g_y_0_y_0_x_zz_0_yz, g_y_0_y_0_x_zz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_zz_0_xx[i] = 4.0 * g_xy_zz_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_zz_0_xy[i] = 4.0 * g_xy_zz_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_zz_0_xz[i] = 4.0 * g_xy_zz_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_zz_0_yy[i] = 4.0 * g_xy_zz_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_zz_0_yz[i] = 4.0 * g_xy_zz_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_zz_0_zz[i] = 4.0 * g_xy_zz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (468-474)

    #pragma omp simd aligned(g_0_xx_y_xx, g_0_xx_y_xy, g_0_xx_y_xz, g_0_xx_y_yy, g_0_xx_y_yz, g_0_xx_y_zz, g_y_0_y_0_y_xx_0_xx, g_y_0_y_0_y_xx_0_xy, g_y_0_y_0_y_xx_0_xz, g_y_0_y_0_y_xx_0_yy, g_y_0_y_0_y_xx_0_yz, g_y_0_y_0_y_xx_0_zz, g_yy_xx_y_xx, g_yy_xx_y_xy, g_yy_xx_y_xz, g_yy_xx_y_yy, g_yy_xx_y_yz, g_yy_xx_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_xx_0_xx[i] = -2.0 * g_0_xx_y_xx[i] * c_exps[i] + 4.0 * g_yy_xx_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xx_0_xy[i] = -2.0 * g_0_xx_y_xy[i] * c_exps[i] + 4.0 * g_yy_xx_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xx_0_xz[i] = -2.0 * g_0_xx_y_xz[i] * c_exps[i] + 4.0 * g_yy_xx_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xx_0_yy[i] = -2.0 * g_0_xx_y_yy[i] * c_exps[i] + 4.0 * g_yy_xx_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xx_0_yz[i] = -2.0 * g_0_xx_y_yz[i] * c_exps[i] + 4.0 * g_yy_xx_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xx_0_zz[i] = -2.0 * g_0_xx_y_zz[i] * c_exps[i] + 4.0 * g_yy_xx_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (474-480)

    #pragma omp simd aligned(g_0_xy_y_xx, g_0_xy_y_xy, g_0_xy_y_xz, g_0_xy_y_yy, g_0_xy_y_yz, g_0_xy_y_zz, g_y_0_y_0_y_xy_0_xx, g_y_0_y_0_y_xy_0_xy, g_y_0_y_0_y_xy_0_xz, g_y_0_y_0_y_xy_0_yy, g_y_0_y_0_y_xy_0_yz, g_y_0_y_0_y_xy_0_zz, g_yy_xy_y_xx, g_yy_xy_y_xy, g_yy_xy_y_xz, g_yy_xy_y_yy, g_yy_xy_y_yz, g_yy_xy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_xy_0_xx[i] = -2.0 * g_0_xy_y_xx[i] * c_exps[i] + 4.0 * g_yy_xy_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xy_0_xy[i] = -2.0 * g_0_xy_y_xy[i] * c_exps[i] + 4.0 * g_yy_xy_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xy_0_xz[i] = -2.0 * g_0_xy_y_xz[i] * c_exps[i] + 4.0 * g_yy_xy_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xy_0_yy[i] = -2.0 * g_0_xy_y_yy[i] * c_exps[i] + 4.0 * g_yy_xy_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xy_0_yz[i] = -2.0 * g_0_xy_y_yz[i] * c_exps[i] + 4.0 * g_yy_xy_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xy_0_zz[i] = -2.0 * g_0_xy_y_zz[i] * c_exps[i] + 4.0 * g_yy_xy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (480-486)

    #pragma omp simd aligned(g_0_xz_y_xx, g_0_xz_y_xy, g_0_xz_y_xz, g_0_xz_y_yy, g_0_xz_y_yz, g_0_xz_y_zz, g_y_0_y_0_y_xz_0_xx, g_y_0_y_0_y_xz_0_xy, g_y_0_y_0_y_xz_0_xz, g_y_0_y_0_y_xz_0_yy, g_y_0_y_0_y_xz_0_yz, g_y_0_y_0_y_xz_0_zz, g_yy_xz_y_xx, g_yy_xz_y_xy, g_yy_xz_y_xz, g_yy_xz_y_yy, g_yy_xz_y_yz, g_yy_xz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_xz_0_xx[i] = -2.0 * g_0_xz_y_xx[i] * c_exps[i] + 4.0 * g_yy_xz_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xz_0_xy[i] = -2.0 * g_0_xz_y_xy[i] * c_exps[i] + 4.0 * g_yy_xz_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xz_0_xz[i] = -2.0 * g_0_xz_y_xz[i] * c_exps[i] + 4.0 * g_yy_xz_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xz_0_yy[i] = -2.0 * g_0_xz_y_yy[i] * c_exps[i] + 4.0 * g_yy_xz_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xz_0_yz[i] = -2.0 * g_0_xz_y_yz[i] * c_exps[i] + 4.0 * g_yy_xz_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xz_0_zz[i] = -2.0 * g_0_xz_y_zz[i] * c_exps[i] + 4.0 * g_yy_xz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (486-492)

    #pragma omp simd aligned(g_0_yy_y_xx, g_0_yy_y_xy, g_0_yy_y_xz, g_0_yy_y_yy, g_0_yy_y_yz, g_0_yy_y_zz, g_y_0_y_0_y_yy_0_xx, g_y_0_y_0_y_yy_0_xy, g_y_0_y_0_y_yy_0_xz, g_y_0_y_0_y_yy_0_yy, g_y_0_y_0_y_yy_0_yz, g_y_0_y_0_y_yy_0_zz, g_yy_yy_y_xx, g_yy_yy_y_xy, g_yy_yy_y_xz, g_yy_yy_y_yy, g_yy_yy_y_yz, g_yy_yy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_yy_0_xx[i] = -2.0 * g_0_yy_y_xx[i] * c_exps[i] + 4.0 * g_yy_yy_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_yy_0_xy[i] = -2.0 * g_0_yy_y_xy[i] * c_exps[i] + 4.0 * g_yy_yy_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_yy_0_xz[i] = -2.0 * g_0_yy_y_xz[i] * c_exps[i] + 4.0 * g_yy_yy_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_yy_0_yy[i] = -2.0 * g_0_yy_y_yy[i] * c_exps[i] + 4.0 * g_yy_yy_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_yy_0_yz[i] = -2.0 * g_0_yy_y_yz[i] * c_exps[i] + 4.0 * g_yy_yy_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_yy_0_zz[i] = -2.0 * g_0_yy_y_zz[i] * c_exps[i] + 4.0 * g_yy_yy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (492-498)

    #pragma omp simd aligned(g_0_yz_y_xx, g_0_yz_y_xy, g_0_yz_y_xz, g_0_yz_y_yy, g_0_yz_y_yz, g_0_yz_y_zz, g_y_0_y_0_y_yz_0_xx, g_y_0_y_0_y_yz_0_xy, g_y_0_y_0_y_yz_0_xz, g_y_0_y_0_y_yz_0_yy, g_y_0_y_0_y_yz_0_yz, g_y_0_y_0_y_yz_0_zz, g_yy_yz_y_xx, g_yy_yz_y_xy, g_yy_yz_y_xz, g_yy_yz_y_yy, g_yy_yz_y_yz, g_yy_yz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_yz_0_xx[i] = -2.0 * g_0_yz_y_xx[i] * c_exps[i] + 4.0 * g_yy_yz_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_yz_0_xy[i] = -2.0 * g_0_yz_y_xy[i] * c_exps[i] + 4.0 * g_yy_yz_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_yz_0_xz[i] = -2.0 * g_0_yz_y_xz[i] * c_exps[i] + 4.0 * g_yy_yz_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_yz_0_yy[i] = -2.0 * g_0_yz_y_yy[i] * c_exps[i] + 4.0 * g_yy_yz_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_yz_0_yz[i] = -2.0 * g_0_yz_y_yz[i] * c_exps[i] + 4.0 * g_yy_yz_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_yz_0_zz[i] = -2.0 * g_0_yz_y_zz[i] * c_exps[i] + 4.0 * g_yy_yz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (498-504)

    #pragma omp simd aligned(g_0_zz_y_xx, g_0_zz_y_xy, g_0_zz_y_xz, g_0_zz_y_yy, g_0_zz_y_yz, g_0_zz_y_zz, g_y_0_y_0_y_zz_0_xx, g_y_0_y_0_y_zz_0_xy, g_y_0_y_0_y_zz_0_xz, g_y_0_y_0_y_zz_0_yy, g_y_0_y_0_y_zz_0_yz, g_y_0_y_0_y_zz_0_zz, g_yy_zz_y_xx, g_yy_zz_y_xy, g_yy_zz_y_xz, g_yy_zz_y_yy, g_yy_zz_y_yz, g_yy_zz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_zz_0_xx[i] = -2.0 * g_0_zz_y_xx[i] * c_exps[i] + 4.0 * g_yy_zz_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_zz_0_xy[i] = -2.0 * g_0_zz_y_xy[i] * c_exps[i] + 4.0 * g_yy_zz_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_zz_0_xz[i] = -2.0 * g_0_zz_y_xz[i] * c_exps[i] + 4.0 * g_yy_zz_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_zz_0_yy[i] = -2.0 * g_0_zz_y_yy[i] * c_exps[i] + 4.0 * g_yy_zz_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_zz_0_yz[i] = -2.0 * g_0_zz_y_yz[i] * c_exps[i] + 4.0 * g_yy_zz_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_zz_0_zz[i] = -2.0 * g_0_zz_y_zz[i] * c_exps[i] + 4.0 * g_yy_zz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (504-510)

    #pragma omp simd aligned(g_y_0_y_0_z_xx_0_xx, g_y_0_y_0_z_xx_0_xy, g_y_0_y_0_z_xx_0_xz, g_y_0_y_0_z_xx_0_yy, g_y_0_y_0_z_xx_0_yz, g_y_0_y_0_z_xx_0_zz, g_yz_xx_y_xx, g_yz_xx_y_xy, g_yz_xx_y_xz, g_yz_xx_y_yy, g_yz_xx_y_yz, g_yz_xx_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_xx_0_xx[i] = 4.0 * g_yz_xx_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xx_0_xy[i] = 4.0 * g_yz_xx_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xx_0_xz[i] = 4.0 * g_yz_xx_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xx_0_yy[i] = 4.0 * g_yz_xx_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xx_0_yz[i] = 4.0 * g_yz_xx_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xx_0_zz[i] = 4.0 * g_yz_xx_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (510-516)

    #pragma omp simd aligned(g_y_0_y_0_z_xy_0_xx, g_y_0_y_0_z_xy_0_xy, g_y_0_y_0_z_xy_0_xz, g_y_0_y_0_z_xy_0_yy, g_y_0_y_0_z_xy_0_yz, g_y_0_y_0_z_xy_0_zz, g_yz_xy_y_xx, g_yz_xy_y_xy, g_yz_xy_y_xz, g_yz_xy_y_yy, g_yz_xy_y_yz, g_yz_xy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_xy_0_xx[i] = 4.0 * g_yz_xy_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xy_0_xy[i] = 4.0 * g_yz_xy_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xy_0_xz[i] = 4.0 * g_yz_xy_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xy_0_yy[i] = 4.0 * g_yz_xy_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xy_0_yz[i] = 4.0 * g_yz_xy_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xy_0_zz[i] = 4.0 * g_yz_xy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (516-522)

    #pragma omp simd aligned(g_y_0_y_0_z_xz_0_xx, g_y_0_y_0_z_xz_0_xy, g_y_0_y_0_z_xz_0_xz, g_y_0_y_0_z_xz_0_yy, g_y_0_y_0_z_xz_0_yz, g_y_0_y_0_z_xz_0_zz, g_yz_xz_y_xx, g_yz_xz_y_xy, g_yz_xz_y_xz, g_yz_xz_y_yy, g_yz_xz_y_yz, g_yz_xz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_xz_0_xx[i] = 4.0 * g_yz_xz_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xz_0_xy[i] = 4.0 * g_yz_xz_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xz_0_xz[i] = 4.0 * g_yz_xz_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xz_0_yy[i] = 4.0 * g_yz_xz_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xz_0_yz[i] = 4.0 * g_yz_xz_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xz_0_zz[i] = 4.0 * g_yz_xz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (522-528)

    #pragma omp simd aligned(g_y_0_y_0_z_yy_0_xx, g_y_0_y_0_z_yy_0_xy, g_y_0_y_0_z_yy_0_xz, g_y_0_y_0_z_yy_0_yy, g_y_0_y_0_z_yy_0_yz, g_y_0_y_0_z_yy_0_zz, g_yz_yy_y_xx, g_yz_yy_y_xy, g_yz_yy_y_xz, g_yz_yy_y_yy, g_yz_yy_y_yz, g_yz_yy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_yy_0_xx[i] = 4.0 * g_yz_yy_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_yy_0_xy[i] = 4.0 * g_yz_yy_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_yy_0_xz[i] = 4.0 * g_yz_yy_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_yy_0_yy[i] = 4.0 * g_yz_yy_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_yy_0_yz[i] = 4.0 * g_yz_yy_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_yy_0_zz[i] = 4.0 * g_yz_yy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (528-534)

    #pragma omp simd aligned(g_y_0_y_0_z_yz_0_xx, g_y_0_y_0_z_yz_0_xy, g_y_0_y_0_z_yz_0_xz, g_y_0_y_0_z_yz_0_yy, g_y_0_y_0_z_yz_0_yz, g_y_0_y_0_z_yz_0_zz, g_yz_yz_y_xx, g_yz_yz_y_xy, g_yz_yz_y_xz, g_yz_yz_y_yy, g_yz_yz_y_yz, g_yz_yz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_yz_0_xx[i] = 4.0 * g_yz_yz_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_yz_0_xy[i] = 4.0 * g_yz_yz_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_yz_0_xz[i] = 4.0 * g_yz_yz_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_yz_0_yy[i] = 4.0 * g_yz_yz_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_yz_0_yz[i] = 4.0 * g_yz_yz_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_yz_0_zz[i] = 4.0 * g_yz_yz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (534-540)

    #pragma omp simd aligned(g_y_0_y_0_z_zz_0_xx, g_y_0_y_0_z_zz_0_xy, g_y_0_y_0_z_zz_0_xz, g_y_0_y_0_z_zz_0_yy, g_y_0_y_0_z_zz_0_yz, g_y_0_y_0_z_zz_0_zz, g_yz_zz_y_xx, g_yz_zz_y_xy, g_yz_zz_y_xz, g_yz_zz_y_yy, g_yz_zz_y_yz, g_yz_zz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_zz_0_xx[i] = 4.0 * g_yz_zz_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_zz_0_xy[i] = 4.0 * g_yz_zz_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_zz_0_xz[i] = 4.0 * g_yz_zz_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_zz_0_yy[i] = 4.0 * g_yz_zz_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_zz_0_yz[i] = 4.0 * g_yz_zz_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_zz_0_zz[i] = 4.0 * g_yz_zz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (540-546)

    #pragma omp simd aligned(g_xy_xx_z_xx, g_xy_xx_z_xy, g_xy_xx_z_xz, g_xy_xx_z_yy, g_xy_xx_z_yz, g_xy_xx_z_zz, g_y_0_z_0_x_xx_0_xx, g_y_0_z_0_x_xx_0_xy, g_y_0_z_0_x_xx_0_xz, g_y_0_z_0_x_xx_0_yy, g_y_0_z_0_x_xx_0_yz, g_y_0_z_0_x_xx_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_xx_0_xx[i] = 4.0 * g_xy_xx_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xx_0_xy[i] = 4.0 * g_xy_xx_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xx_0_xz[i] = 4.0 * g_xy_xx_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xx_0_yy[i] = 4.0 * g_xy_xx_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xx_0_yz[i] = 4.0 * g_xy_xx_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xx_0_zz[i] = 4.0 * g_xy_xx_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (546-552)

    #pragma omp simd aligned(g_xy_xy_z_xx, g_xy_xy_z_xy, g_xy_xy_z_xz, g_xy_xy_z_yy, g_xy_xy_z_yz, g_xy_xy_z_zz, g_y_0_z_0_x_xy_0_xx, g_y_0_z_0_x_xy_0_xy, g_y_0_z_0_x_xy_0_xz, g_y_0_z_0_x_xy_0_yy, g_y_0_z_0_x_xy_0_yz, g_y_0_z_0_x_xy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_xy_0_xx[i] = 4.0 * g_xy_xy_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xy_0_xy[i] = 4.0 * g_xy_xy_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xy_0_xz[i] = 4.0 * g_xy_xy_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xy_0_yy[i] = 4.0 * g_xy_xy_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xy_0_yz[i] = 4.0 * g_xy_xy_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xy_0_zz[i] = 4.0 * g_xy_xy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (552-558)

    #pragma omp simd aligned(g_xy_xz_z_xx, g_xy_xz_z_xy, g_xy_xz_z_xz, g_xy_xz_z_yy, g_xy_xz_z_yz, g_xy_xz_z_zz, g_y_0_z_0_x_xz_0_xx, g_y_0_z_0_x_xz_0_xy, g_y_0_z_0_x_xz_0_xz, g_y_0_z_0_x_xz_0_yy, g_y_0_z_0_x_xz_0_yz, g_y_0_z_0_x_xz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_xz_0_xx[i] = 4.0 * g_xy_xz_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xz_0_xy[i] = 4.0 * g_xy_xz_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xz_0_xz[i] = 4.0 * g_xy_xz_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xz_0_yy[i] = 4.0 * g_xy_xz_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xz_0_yz[i] = 4.0 * g_xy_xz_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xz_0_zz[i] = 4.0 * g_xy_xz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (558-564)

    #pragma omp simd aligned(g_xy_yy_z_xx, g_xy_yy_z_xy, g_xy_yy_z_xz, g_xy_yy_z_yy, g_xy_yy_z_yz, g_xy_yy_z_zz, g_y_0_z_0_x_yy_0_xx, g_y_0_z_0_x_yy_0_xy, g_y_0_z_0_x_yy_0_xz, g_y_0_z_0_x_yy_0_yy, g_y_0_z_0_x_yy_0_yz, g_y_0_z_0_x_yy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_yy_0_xx[i] = 4.0 * g_xy_yy_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_yy_0_xy[i] = 4.0 * g_xy_yy_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_yy_0_xz[i] = 4.0 * g_xy_yy_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_yy_0_yy[i] = 4.0 * g_xy_yy_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_yy_0_yz[i] = 4.0 * g_xy_yy_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_yy_0_zz[i] = 4.0 * g_xy_yy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (564-570)

    #pragma omp simd aligned(g_xy_yz_z_xx, g_xy_yz_z_xy, g_xy_yz_z_xz, g_xy_yz_z_yy, g_xy_yz_z_yz, g_xy_yz_z_zz, g_y_0_z_0_x_yz_0_xx, g_y_0_z_0_x_yz_0_xy, g_y_0_z_0_x_yz_0_xz, g_y_0_z_0_x_yz_0_yy, g_y_0_z_0_x_yz_0_yz, g_y_0_z_0_x_yz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_yz_0_xx[i] = 4.0 * g_xy_yz_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_yz_0_xy[i] = 4.0 * g_xy_yz_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_yz_0_xz[i] = 4.0 * g_xy_yz_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_yz_0_yy[i] = 4.0 * g_xy_yz_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_yz_0_yz[i] = 4.0 * g_xy_yz_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_yz_0_zz[i] = 4.0 * g_xy_yz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (570-576)

    #pragma omp simd aligned(g_xy_zz_z_xx, g_xy_zz_z_xy, g_xy_zz_z_xz, g_xy_zz_z_yy, g_xy_zz_z_yz, g_xy_zz_z_zz, g_y_0_z_0_x_zz_0_xx, g_y_0_z_0_x_zz_0_xy, g_y_0_z_0_x_zz_0_xz, g_y_0_z_0_x_zz_0_yy, g_y_0_z_0_x_zz_0_yz, g_y_0_z_0_x_zz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_zz_0_xx[i] = 4.0 * g_xy_zz_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_zz_0_xy[i] = 4.0 * g_xy_zz_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_zz_0_xz[i] = 4.0 * g_xy_zz_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_zz_0_yy[i] = 4.0 * g_xy_zz_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_zz_0_yz[i] = 4.0 * g_xy_zz_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_zz_0_zz[i] = 4.0 * g_xy_zz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (576-582)

    #pragma omp simd aligned(g_0_xx_z_xx, g_0_xx_z_xy, g_0_xx_z_xz, g_0_xx_z_yy, g_0_xx_z_yz, g_0_xx_z_zz, g_y_0_z_0_y_xx_0_xx, g_y_0_z_0_y_xx_0_xy, g_y_0_z_0_y_xx_0_xz, g_y_0_z_0_y_xx_0_yy, g_y_0_z_0_y_xx_0_yz, g_y_0_z_0_y_xx_0_zz, g_yy_xx_z_xx, g_yy_xx_z_xy, g_yy_xx_z_xz, g_yy_xx_z_yy, g_yy_xx_z_yz, g_yy_xx_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_xx_0_xx[i] = -2.0 * g_0_xx_z_xx[i] * c_exps[i] + 4.0 * g_yy_xx_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xx_0_xy[i] = -2.0 * g_0_xx_z_xy[i] * c_exps[i] + 4.0 * g_yy_xx_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xx_0_xz[i] = -2.0 * g_0_xx_z_xz[i] * c_exps[i] + 4.0 * g_yy_xx_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xx_0_yy[i] = -2.0 * g_0_xx_z_yy[i] * c_exps[i] + 4.0 * g_yy_xx_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xx_0_yz[i] = -2.0 * g_0_xx_z_yz[i] * c_exps[i] + 4.0 * g_yy_xx_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xx_0_zz[i] = -2.0 * g_0_xx_z_zz[i] * c_exps[i] + 4.0 * g_yy_xx_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (582-588)

    #pragma omp simd aligned(g_0_xy_z_xx, g_0_xy_z_xy, g_0_xy_z_xz, g_0_xy_z_yy, g_0_xy_z_yz, g_0_xy_z_zz, g_y_0_z_0_y_xy_0_xx, g_y_0_z_0_y_xy_0_xy, g_y_0_z_0_y_xy_0_xz, g_y_0_z_0_y_xy_0_yy, g_y_0_z_0_y_xy_0_yz, g_y_0_z_0_y_xy_0_zz, g_yy_xy_z_xx, g_yy_xy_z_xy, g_yy_xy_z_xz, g_yy_xy_z_yy, g_yy_xy_z_yz, g_yy_xy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_xy_0_xx[i] = -2.0 * g_0_xy_z_xx[i] * c_exps[i] + 4.0 * g_yy_xy_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xy_0_xy[i] = -2.0 * g_0_xy_z_xy[i] * c_exps[i] + 4.0 * g_yy_xy_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xy_0_xz[i] = -2.0 * g_0_xy_z_xz[i] * c_exps[i] + 4.0 * g_yy_xy_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xy_0_yy[i] = -2.0 * g_0_xy_z_yy[i] * c_exps[i] + 4.0 * g_yy_xy_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xy_0_yz[i] = -2.0 * g_0_xy_z_yz[i] * c_exps[i] + 4.0 * g_yy_xy_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xy_0_zz[i] = -2.0 * g_0_xy_z_zz[i] * c_exps[i] + 4.0 * g_yy_xy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (588-594)

    #pragma omp simd aligned(g_0_xz_z_xx, g_0_xz_z_xy, g_0_xz_z_xz, g_0_xz_z_yy, g_0_xz_z_yz, g_0_xz_z_zz, g_y_0_z_0_y_xz_0_xx, g_y_0_z_0_y_xz_0_xy, g_y_0_z_0_y_xz_0_xz, g_y_0_z_0_y_xz_0_yy, g_y_0_z_0_y_xz_0_yz, g_y_0_z_0_y_xz_0_zz, g_yy_xz_z_xx, g_yy_xz_z_xy, g_yy_xz_z_xz, g_yy_xz_z_yy, g_yy_xz_z_yz, g_yy_xz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_xz_0_xx[i] = -2.0 * g_0_xz_z_xx[i] * c_exps[i] + 4.0 * g_yy_xz_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xz_0_xy[i] = -2.0 * g_0_xz_z_xy[i] * c_exps[i] + 4.0 * g_yy_xz_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xz_0_xz[i] = -2.0 * g_0_xz_z_xz[i] * c_exps[i] + 4.0 * g_yy_xz_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xz_0_yy[i] = -2.0 * g_0_xz_z_yy[i] * c_exps[i] + 4.0 * g_yy_xz_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xz_0_yz[i] = -2.0 * g_0_xz_z_yz[i] * c_exps[i] + 4.0 * g_yy_xz_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xz_0_zz[i] = -2.0 * g_0_xz_z_zz[i] * c_exps[i] + 4.0 * g_yy_xz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (594-600)

    #pragma omp simd aligned(g_0_yy_z_xx, g_0_yy_z_xy, g_0_yy_z_xz, g_0_yy_z_yy, g_0_yy_z_yz, g_0_yy_z_zz, g_y_0_z_0_y_yy_0_xx, g_y_0_z_0_y_yy_0_xy, g_y_0_z_0_y_yy_0_xz, g_y_0_z_0_y_yy_0_yy, g_y_0_z_0_y_yy_0_yz, g_y_0_z_0_y_yy_0_zz, g_yy_yy_z_xx, g_yy_yy_z_xy, g_yy_yy_z_xz, g_yy_yy_z_yy, g_yy_yy_z_yz, g_yy_yy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_yy_0_xx[i] = -2.0 * g_0_yy_z_xx[i] * c_exps[i] + 4.0 * g_yy_yy_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_yy_0_xy[i] = -2.0 * g_0_yy_z_xy[i] * c_exps[i] + 4.0 * g_yy_yy_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_yy_0_xz[i] = -2.0 * g_0_yy_z_xz[i] * c_exps[i] + 4.0 * g_yy_yy_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_yy_0_yy[i] = -2.0 * g_0_yy_z_yy[i] * c_exps[i] + 4.0 * g_yy_yy_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_yy_0_yz[i] = -2.0 * g_0_yy_z_yz[i] * c_exps[i] + 4.0 * g_yy_yy_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_yy_0_zz[i] = -2.0 * g_0_yy_z_zz[i] * c_exps[i] + 4.0 * g_yy_yy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (600-606)

    #pragma omp simd aligned(g_0_yz_z_xx, g_0_yz_z_xy, g_0_yz_z_xz, g_0_yz_z_yy, g_0_yz_z_yz, g_0_yz_z_zz, g_y_0_z_0_y_yz_0_xx, g_y_0_z_0_y_yz_0_xy, g_y_0_z_0_y_yz_0_xz, g_y_0_z_0_y_yz_0_yy, g_y_0_z_0_y_yz_0_yz, g_y_0_z_0_y_yz_0_zz, g_yy_yz_z_xx, g_yy_yz_z_xy, g_yy_yz_z_xz, g_yy_yz_z_yy, g_yy_yz_z_yz, g_yy_yz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_yz_0_xx[i] = -2.0 * g_0_yz_z_xx[i] * c_exps[i] + 4.0 * g_yy_yz_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_yz_0_xy[i] = -2.0 * g_0_yz_z_xy[i] * c_exps[i] + 4.0 * g_yy_yz_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_yz_0_xz[i] = -2.0 * g_0_yz_z_xz[i] * c_exps[i] + 4.0 * g_yy_yz_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_yz_0_yy[i] = -2.0 * g_0_yz_z_yy[i] * c_exps[i] + 4.0 * g_yy_yz_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_yz_0_yz[i] = -2.0 * g_0_yz_z_yz[i] * c_exps[i] + 4.0 * g_yy_yz_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_yz_0_zz[i] = -2.0 * g_0_yz_z_zz[i] * c_exps[i] + 4.0 * g_yy_yz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (606-612)

    #pragma omp simd aligned(g_0_zz_z_xx, g_0_zz_z_xy, g_0_zz_z_xz, g_0_zz_z_yy, g_0_zz_z_yz, g_0_zz_z_zz, g_y_0_z_0_y_zz_0_xx, g_y_0_z_0_y_zz_0_xy, g_y_0_z_0_y_zz_0_xz, g_y_0_z_0_y_zz_0_yy, g_y_0_z_0_y_zz_0_yz, g_y_0_z_0_y_zz_0_zz, g_yy_zz_z_xx, g_yy_zz_z_xy, g_yy_zz_z_xz, g_yy_zz_z_yy, g_yy_zz_z_yz, g_yy_zz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_zz_0_xx[i] = -2.0 * g_0_zz_z_xx[i] * c_exps[i] + 4.0 * g_yy_zz_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_zz_0_xy[i] = -2.0 * g_0_zz_z_xy[i] * c_exps[i] + 4.0 * g_yy_zz_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_zz_0_xz[i] = -2.0 * g_0_zz_z_xz[i] * c_exps[i] + 4.0 * g_yy_zz_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_zz_0_yy[i] = -2.0 * g_0_zz_z_yy[i] * c_exps[i] + 4.0 * g_yy_zz_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_zz_0_yz[i] = -2.0 * g_0_zz_z_yz[i] * c_exps[i] + 4.0 * g_yy_zz_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_zz_0_zz[i] = -2.0 * g_0_zz_z_zz[i] * c_exps[i] + 4.0 * g_yy_zz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (612-618)

    #pragma omp simd aligned(g_y_0_z_0_z_xx_0_xx, g_y_0_z_0_z_xx_0_xy, g_y_0_z_0_z_xx_0_xz, g_y_0_z_0_z_xx_0_yy, g_y_0_z_0_z_xx_0_yz, g_y_0_z_0_z_xx_0_zz, g_yz_xx_z_xx, g_yz_xx_z_xy, g_yz_xx_z_xz, g_yz_xx_z_yy, g_yz_xx_z_yz, g_yz_xx_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_xx_0_xx[i] = 4.0 * g_yz_xx_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xx_0_xy[i] = 4.0 * g_yz_xx_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xx_0_xz[i] = 4.0 * g_yz_xx_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xx_0_yy[i] = 4.0 * g_yz_xx_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xx_0_yz[i] = 4.0 * g_yz_xx_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xx_0_zz[i] = 4.0 * g_yz_xx_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (618-624)

    #pragma omp simd aligned(g_y_0_z_0_z_xy_0_xx, g_y_0_z_0_z_xy_0_xy, g_y_0_z_0_z_xy_0_xz, g_y_0_z_0_z_xy_0_yy, g_y_0_z_0_z_xy_0_yz, g_y_0_z_0_z_xy_0_zz, g_yz_xy_z_xx, g_yz_xy_z_xy, g_yz_xy_z_xz, g_yz_xy_z_yy, g_yz_xy_z_yz, g_yz_xy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_xy_0_xx[i] = 4.0 * g_yz_xy_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xy_0_xy[i] = 4.0 * g_yz_xy_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xy_0_xz[i] = 4.0 * g_yz_xy_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xy_0_yy[i] = 4.0 * g_yz_xy_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xy_0_yz[i] = 4.0 * g_yz_xy_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xy_0_zz[i] = 4.0 * g_yz_xy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (624-630)

    #pragma omp simd aligned(g_y_0_z_0_z_xz_0_xx, g_y_0_z_0_z_xz_0_xy, g_y_0_z_0_z_xz_0_xz, g_y_0_z_0_z_xz_0_yy, g_y_0_z_0_z_xz_0_yz, g_y_0_z_0_z_xz_0_zz, g_yz_xz_z_xx, g_yz_xz_z_xy, g_yz_xz_z_xz, g_yz_xz_z_yy, g_yz_xz_z_yz, g_yz_xz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_xz_0_xx[i] = 4.0 * g_yz_xz_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xz_0_xy[i] = 4.0 * g_yz_xz_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xz_0_xz[i] = 4.0 * g_yz_xz_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xz_0_yy[i] = 4.0 * g_yz_xz_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xz_0_yz[i] = 4.0 * g_yz_xz_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xz_0_zz[i] = 4.0 * g_yz_xz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (630-636)

    #pragma omp simd aligned(g_y_0_z_0_z_yy_0_xx, g_y_0_z_0_z_yy_0_xy, g_y_0_z_0_z_yy_0_xz, g_y_0_z_0_z_yy_0_yy, g_y_0_z_0_z_yy_0_yz, g_y_0_z_0_z_yy_0_zz, g_yz_yy_z_xx, g_yz_yy_z_xy, g_yz_yy_z_xz, g_yz_yy_z_yy, g_yz_yy_z_yz, g_yz_yy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_yy_0_xx[i] = 4.0 * g_yz_yy_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_yy_0_xy[i] = 4.0 * g_yz_yy_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_yy_0_xz[i] = 4.0 * g_yz_yy_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_yy_0_yy[i] = 4.0 * g_yz_yy_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_yy_0_yz[i] = 4.0 * g_yz_yy_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_yy_0_zz[i] = 4.0 * g_yz_yy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (636-642)

    #pragma omp simd aligned(g_y_0_z_0_z_yz_0_xx, g_y_0_z_0_z_yz_0_xy, g_y_0_z_0_z_yz_0_xz, g_y_0_z_0_z_yz_0_yy, g_y_0_z_0_z_yz_0_yz, g_y_0_z_0_z_yz_0_zz, g_yz_yz_z_xx, g_yz_yz_z_xy, g_yz_yz_z_xz, g_yz_yz_z_yy, g_yz_yz_z_yz, g_yz_yz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_yz_0_xx[i] = 4.0 * g_yz_yz_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_yz_0_xy[i] = 4.0 * g_yz_yz_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_yz_0_xz[i] = 4.0 * g_yz_yz_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_yz_0_yy[i] = 4.0 * g_yz_yz_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_yz_0_yz[i] = 4.0 * g_yz_yz_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_yz_0_zz[i] = 4.0 * g_yz_yz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (642-648)

    #pragma omp simd aligned(g_y_0_z_0_z_zz_0_xx, g_y_0_z_0_z_zz_0_xy, g_y_0_z_0_z_zz_0_xz, g_y_0_z_0_z_zz_0_yy, g_y_0_z_0_z_zz_0_yz, g_y_0_z_0_z_zz_0_zz, g_yz_zz_z_xx, g_yz_zz_z_xy, g_yz_zz_z_xz, g_yz_zz_z_yy, g_yz_zz_z_yz, g_yz_zz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_zz_0_xx[i] = 4.0 * g_yz_zz_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_zz_0_xy[i] = 4.0 * g_yz_zz_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_zz_0_xz[i] = 4.0 * g_yz_zz_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_zz_0_yy[i] = 4.0 * g_yz_zz_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_zz_0_yz[i] = 4.0 * g_yz_zz_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_zz_0_zz[i] = 4.0 * g_yz_zz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (648-654)

    #pragma omp simd aligned(g_xz_xx_x_xx, g_xz_xx_x_xy, g_xz_xx_x_xz, g_xz_xx_x_yy, g_xz_xx_x_yz, g_xz_xx_x_zz, g_z_0_x_0_x_xx_0_xx, g_z_0_x_0_x_xx_0_xy, g_z_0_x_0_x_xx_0_xz, g_z_0_x_0_x_xx_0_yy, g_z_0_x_0_x_xx_0_yz, g_z_0_x_0_x_xx_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_xx_0_xx[i] = 4.0 * g_xz_xx_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xx_0_xy[i] = 4.0 * g_xz_xx_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xx_0_xz[i] = 4.0 * g_xz_xx_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xx_0_yy[i] = 4.0 * g_xz_xx_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xx_0_yz[i] = 4.0 * g_xz_xx_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xx_0_zz[i] = 4.0 * g_xz_xx_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (654-660)

    #pragma omp simd aligned(g_xz_xy_x_xx, g_xz_xy_x_xy, g_xz_xy_x_xz, g_xz_xy_x_yy, g_xz_xy_x_yz, g_xz_xy_x_zz, g_z_0_x_0_x_xy_0_xx, g_z_0_x_0_x_xy_0_xy, g_z_0_x_0_x_xy_0_xz, g_z_0_x_0_x_xy_0_yy, g_z_0_x_0_x_xy_0_yz, g_z_0_x_0_x_xy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_xy_0_xx[i] = 4.0 * g_xz_xy_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xy_0_xy[i] = 4.0 * g_xz_xy_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xy_0_xz[i] = 4.0 * g_xz_xy_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xy_0_yy[i] = 4.0 * g_xz_xy_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xy_0_yz[i] = 4.0 * g_xz_xy_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xy_0_zz[i] = 4.0 * g_xz_xy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (660-666)

    #pragma omp simd aligned(g_xz_xz_x_xx, g_xz_xz_x_xy, g_xz_xz_x_xz, g_xz_xz_x_yy, g_xz_xz_x_yz, g_xz_xz_x_zz, g_z_0_x_0_x_xz_0_xx, g_z_0_x_0_x_xz_0_xy, g_z_0_x_0_x_xz_0_xz, g_z_0_x_0_x_xz_0_yy, g_z_0_x_0_x_xz_0_yz, g_z_0_x_0_x_xz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_xz_0_xx[i] = 4.0 * g_xz_xz_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xz_0_xy[i] = 4.0 * g_xz_xz_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xz_0_xz[i] = 4.0 * g_xz_xz_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xz_0_yy[i] = 4.0 * g_xz_xz_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xz_0_yz[i] = 4.0 * g_xz_xz_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xz_0_zz[i] = 4.0 * g_xz_xz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (666-672)

    #pragma omp simd aligned(g_xz_yy_x_xx, g_xz_yy_x_xy, g_xz_yy_x_xz, g_xz_yy_x_yy, g_xz_yy_x_yz, g_xz_yy_x_zz, g_z_0_x_0_x_yy_0_xx, g_z_0_x_0_x_yy_0_xy, g_z_0_x_0_x_yy_0_xz, g_z_0_x_0_x_yy_0_yy, g_z_0_x_0_x_yy_0_yz, g_z_0_x_0_x_yy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_yy_0_xx[i] = 4.0 * g_xz_yy_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_yy_0_xy[i] = 4.0 * g_xz_yy_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_yy_0_xz[i] = 4.0 * g_xz_yy_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_yy_0_yy[i] = 4.0 * g_xz_yy_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_yy_0_yz[i] = 4.0 * g_xz_yy_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_yy_0_zz[i] = 4.0 * g_xz_yy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (672-678)

    #pragma omp simd aligned(g_xz_yz_x_xx, g_xz_yz_x_xy, g_xz_yz_x_xz, g_xz_yz_x_yy, g_xz_yz_x_yz, g_xz_yz_x_zz, g_z_0_x_0_x_yz_0_xx, g_z_0_x_0_x_yz_0_xy, g_z_0_x_0_x_yz_0_xz, g_z_0_x_0_x_yz_0_yy, g_z_0_x_0_x_yz_0_yz, g_z_0_x_0_x_yz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_yz_0_xx[i] = 4.0 * g_xz_yz_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_yz_0_xy[i] = 4.0 * g_xz_yz_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_yz_0_xz[i] = 4.0 * g_xz_yz_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_yz_0_yy[i] = 4.0 * g_xz_yz_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_yz_0_yz[i] = 4.0 * g_xz_yz_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_yz_0_zz[i] = 4.0 * g_xz_yz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (678-684)

    #pragma omp simd aligned(g_xz_zz_x_xx, g_xz_zz_x_xy, g_xz_zz_x_xz, g_xz_zz_x_yy, g_xz_zz_x_yz, g_xz_zz_x_zz, g_z_0_x_0_x_zz_0_xx, g_z_0_x_0_x_zz_0_xy, g_z_0_x_0_x_zz_0_xz, g_z_0_x_0_x_zz_0_yy, g_z_0_x_0_x_zz_0_yz, g_z_0_x_0_x_zz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_zz_0_xx[i] = 4.0 * g_xz_zz_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_zz_0_xy[i] = 4.0 * g_xz_zz_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_zz_0_xz[i] = 4.0 * g_xz_zz_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_zz_0_yy[i] = 4.0 * g_xz_zz_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_zz_0_yz[i] = 4.0 * g_xz_zz_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_zz_0_zz[i] = 4.0 * g_xz_zz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (684-690)

    #pragma omp simd aligned(g_yz_xx_x_xx, g_yz_xx_x_xy, g_yz_xx_x_xz, g_yz_xx_x_yy, g_yz_xx_x_yz, g_yz_xx_x_zz, g_z_0_x_0_y_xx_0_xx, g_z_0_x_0_y_xx_0_xy, g_z_0_x_0_y_xx_0_xz, g_z_0_x_0_y_xx_0_yy, g_z_0_x_0_y_xx_0_yz, g_z_0_x_0_y_xx_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_xx_0_xx[i] = 4.0 * g_yz_xx_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xx_0_xy[i] = 4.0 * g_yz_xx_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xx_0_xz[i] = 4.0 * g_yz_xx_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xx_0_yy[i] = 4.0 * g_yz_xx_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xx_0_yz[i] = 4.0 * g_yz_xx_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xx_0_zz[i] = 4.0 * g_yz_xx_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (690-696)

    #pragma omp simd aligned(g_yz_xy_x_xx, g_yz_xy_x_xy, g_yz_xy_x_xz, g_yz_xy_x_yy, g_yz_xy_x_yz, g_yz_xy_x_zz, g_z_0_x_0_y_xy_0_xx, g_z_0_x_0_y_xy_0_xy, g_z_0_x_0_y_xy_0_xz, g_z_0_x_0_y_xy_0_yy, g_z_0_x_0_y_xy_0_yz, g_z_0_x_0_y_xy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_xy_0_xx[i] = 4.0 * g_yz_xy_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xy_0_xy[i] = 4.0 * g_yz_xy_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xy_0_xz[i] = 4.0 * g_yz_xy_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xy_0_yy[i] = 4.0 * g_yz_xy_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xy_0_yz[i] = 4.0 * g_yz_xy_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xy_0_zz[i] = 4.0 * g_yz_xy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (696-702)

    #pragma omp simd aligned(g_yz_xz_x_xx, g_yz_xz_x_xy, g_yz_xz_x_xz, g_yz_xz_x_yy, g_yz_xz_x_yz, g_yz_xz_x_zz, g_z_0_x_0_y_xz_0_xx, g_z_0_x_0_y_xz_0_xy, g_z_0_x_0_y_xz_0_xz, g_z_0_x_0_y_xz_0_yy, g_z_0_x_0_y_xz_0_yz, g_z_0_x_0_y_xz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_xz_0_xx[i] = 4.0 * g_yz_xz_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xz_0_xy[i] = 4.0 * g_yz_xz_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xz_0_xz[i] = 4.0 * g_yz_xz_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xz_0_yy[i] = 4.0 * g_yz_xz_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xz_0_yz[i] = 4.0 * g_yz_xz_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xz_0_zz[i] = 4.0 * g_yz_xz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (702-708)

    #pragma omp simd aligned(g_yz_yy_x_xx, g_yz_yy_x_xy, g_yz_yy_x_xz, g_yz_yy_x_yy, g_yz_yy_x_yz, g_yz_yy_x_zz, g_z_0_x_0_y_yy_0_xx, g_z_0_x_0_y_yy_0_xy, g_z_0_x_0_y_yy_0_xz, g_z_0_x_0_y_yy_0_yy, g_z_0_x_0_y_yy_0_yz, g_z_0_x_0_y_yy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_yy_0_xx[i] = 4.0 * g_yz_yy_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_yy_0_xy[i] = 4.0 * g_yz_yy_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_yy_0_xz[i] = 4.0 * g_yz_yy_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_yy_0_yy[i] = 4.0 * g_yz_yy_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_yy_0_yz[i] = 4.0 * g_yz_yy_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_yy_0_zz[i] = 4.0 * g_yz_yy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (708-714)

    #pragma omp simd aligned(g_yz_yz_x_xx, g_yz_yz_x_xy, g_yz_yz_x_xz, g_yz_yz_x_yy, g_yz_yz_x_yz, g_yz_yz_x_zz, g_z_0_x_0_y_yz_0_xx, g_z_0_x_0_y_yz_0_xy, g_z_0_x_0_y_yz_0_xz, g_z_0_x_0_y_yz_0_yy, g_z_0_x_0_y_yz_0_yz, g_z_0_x_0_y_yz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_yz_0_xx[i] = 4.0 * g_yz_yz_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_yz_0_xy[i] = 4.0 * g_yz_yz_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_yz_0_xz[i] = 4.0 * g_yz_yz_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_yz_0_yy[i] = 4.0 * g_yz_yz_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_yz_0_yz[i] = 4.0 * g_yz_yz_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_yz_0_zz[i] = 4.0 * g_yz_yz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (714-720)

    #pragma omp simd aligned(g_yz_zz_x_xx, g_yz_zz_x_xy, g_yz_zz_x_xz, g_yz_zz_x_yy, g_yz_zz_x_yz, g_yz_zz_x_zz, g_z_0_x_0_y_zz_0_xx, g_z_0_x_0_y_zz_0_xy, g_z_0_x_0_y_zz_0_xz, g_z_0_x_0_y_zz_0_yy, g_z_0_x_0_y_zz_0_yz, g_z_0_x_0_y_zz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_zz_0_xx[i] = 4.0 * g_yz_zz_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_zz_0_xy[i] = 4.0 * g_yz_zz_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_zz_0_xz[i] = 4.0 * g_yz_zz_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_zz_0_yy[i] = 4.0 * g_yz_zz_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_zz_0_yz[i] = 4.0 * g_yz_zz_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_zz_0_zz[i] = 4.0 * g_yz_zz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (720-726)

    #pragma omp simd aligned(g_0_xx_x_xx, g_0_xx_x_xy, g_0_xx_x_xz, g_0_xx_x_yy, g_0_xx_x_yz, g_0_xx_x_zz, g_z_0_x_0_z_xx_0_xx, g_z_0_x_0_z_xx_0_xy, g_z_0_x_0_z_xx_0_xz, g_z_0_x_0_z_xx_0_yy, g_z_0_x_0_z_xx_0_yz, g_z_0_x_0_z_xx_0_zz, g_zz_xx_x_xx, g_zz_xx_x_xy, g_zz_xx_x_xz, g_zz_xx_x_yy, g_zz_xx_x_yz, g_zz_xx_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_xx_0_xx[i] = -2.0 * g_0_xx_x_xx[i] * c_exps[i] + 4.0 * g_zz_xx_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xx_0_xy[i] = -2.0 * g_0_xx_x_xy[i] * c_exps[i] + 4.0 * g_zz_xx_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xx_0_xz[i] = -2.0 * g_0_xx_x_xz[i] * c_exps[i] + 4.0 * g_zz_xx_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xx_0_yy[i] = -2.0 * g_0_xx_x_yy[i] * c_exps[i] + 4.0 * g_zz_xx_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xx_0_yz[i] = -2.0 * g_0_xx_x_yz[i] * c_exps[i] + 4.0 * g_zz_xx_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xx_0_zz[i] = -2.0 * g_0_xx_x_zz[i] * c_exps[i] + 4.0 * g_zz_xx_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (726-732)

    #pragma omp simd aligned(g_0_xy_x_xx, g_0_xy_x_xy, g_0_xy_x_xz, g_0_xy_x_yy, g_0_xy_x_yz, g_0_xy_x_zz, g_z_0_x_0_z_xy_0_xx, g_z_0_x_0_z_xy_0_xy, g_z_0_x_0_z_xy_0_xz, g_z_0_x_0_z_xy_0_yy, g_z_0_x_0_z_xy_0_yz, g_z_0_x_0_z_xy_0_zz, g_zz_xy_x_xx, g_zz_xy_x_xy, g_zz_xy_x_xz, g_zz_xy_x_yy, g_zz_xy_x_yz, g_zz_xy_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_xy_0_xx[i] = -2.0 * g_0_xy_x_xx[i] * c_exps[i] + 4.0 * g_zz_xy_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xy_0_xy[i] = -2.0 * g_0_xy_x_xy[i] * c_exps[i] + 4.0 * g_zz_xy_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xy_0_xz[i] = -2.0 * g_0_xy_x_xz[i] * c_exps[i] + 4.0 * g_zz_xy_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xy_0_yy[i] = -2.0 * g_0_xy_x_yy[i] * c_exps[i] + 4.0 * g_zz_xy_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xy_0_yz[i] = -2.0 * g_0_xy_x_yz[i] * c_exps[i] + 4.0 * g_zz_xy_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xy_0_zz[i] = -2.0 * g_0_xy_x_zz[i] * c_exps[i] + 4.0 * g_zz_xy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (732-738)

    #pragma omp simd aligned(g_0_xz_x_xx, g_0_xz_x_xy, g_0_xz_x_xz, g_0_xz_x_yy, g_0_xz_x_yz, g_0_xz_x_zz, g_z_0_x_0_z_xz_0_xx, g_z_0_x_0_z_xz_0_xy, g_z_0_x_0_z_xz_0_xz, g_z_0_x_0_z_xz_0_yy, g_z_0_x_0_z_xz_0_yz, g_z_0_x_0_z_xz_0_zz, g_zz_xz_x_xx, g_zz_xz_x_xy, g_zz_xz_x_xz, g_zz_xz_x_yy, g_zz_xz_x_yz, g_zz_xz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_xz_0_xx[i] = -2.0 * g_0_xz_x_xx[i] * c_exps[i] + 4.0 * g_zz_xz_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xz_0_xy[i] = -2.0 * g_0_xz_x_xy[i] * c_exps[i] + 4.0 * g_zz_xz_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xz_0_xz[i] = -2.0 * g_0_xz_x_xz[i] * c_exps[i] + 4.0 * g_zz_xz_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xz_0_yy[i] = -2.0 * g_0_xz_x_yy[i] * c_exps[i] + 4.0 * g_zz_xz_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xz_0_yz[i] = -2.0 * g_0_xz_x_yz[i] * c_exps[i] + 4.0 * g_zz_xz_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xz_0_zz[i] = -2.0 * g_0_xz_x_zz[i] * c_exps[i] + 4.0 * g_zz_xz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (738-744)

    #pragma omp simd aligned(g_0_yy_x_xx, g_0_yy_x_xy, g_0_yy_x_xz, g_0_yy_x_yy, g_0_yy_x_yz, g_0_yy_x_zz, g_z_0_x_0_z_yy_0_xx, g_z_0_x_0_z_yy_0_xy, g_z_0_x_0_z_yy_0_xz, g_z_0_x_0_z_yy_0_yy, g_z_0_x_0_z_yy_0_yz, g_z_0_x_0_z_yy_0_zz, g_zz_yy_x_xx, g_zz_yy_x_xy, g_zz_yy_x_xz, g_zz_yy_x_yy, g_zz_yy_x_yz, g_zz_yy_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_yy_0_xx[i] = -2.0 * g_0_yy_x_xx[i] * c_exps[i] + 4.0 * g_zz_yy_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_yy_0_xy[i] = -2.0 * g_0_yy_x_xy[i] * c_exps[i] + 4.0 * g_zz_yy_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_yy_0_xz[i] = -2.0 * g_0_yy_x_xz[i] * c_exps[i] + 4.0 * g_zz_yy_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_yy_0_yy[i] = -2.0 * g_0_yy_x_yy[i] * c_exps[i] + 4.0 * g_zz_yy_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_yy_0_yz[i] = -2.0 * g_0_yy_x_yz[i] * c_exps[i] + 4.0 * g_zz_yy_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_yy_0_zz[i] = -2.0 * g_0_yy_x_zz[i] * c_exps[i] + 4.0 * g_zz_yy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (744-750)

    #pragma omp simd aligned(g_0_yz_x_xx, g_0_yz_x_xy, g_0_yz_x_xz, g_0_yz_x_yy, g_0_yz_x_yz, g_0_yz_x_zz, g_z_0_x_0_z_yz_0_xx, g_z_0_x_0_z_yz_0_xy, g_z_0_x_0_z_yz_0_xz, g_z_0_x_0_z_yz_0_yy, g_z_0_x_0_z_yz_0_yz, g_z_0_x_0_z_yz_0_zz, g_zz_yz_x_xx, g_zz_yz_x_xy, g_zz_yz_x_xz, g_zz_yz_x_yy, g_zz_yz_x_yz, g_zz_yz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_yz_0_xx[i] = -2.0 * g_0_yz_x_xx[i] * c_exps[i] + 4.0 * g_zz_yz_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_yz_0_xy[i] = -2.0 * g_0_yz_x_xy[i] * c_exps[i] + 4.0 * g_zz_yz_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_yz_0_xz[i] = -2.0 * g_0_yz_x_xz[i] * c_exps[i] + 4.0 * g_zz_yz_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_yz_0_yy[i] = -2.0 * g_0_yz_x_yy[i] * c_exps[i] + 4.0 * g_zz_yz_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_yz_0_yz[i] = -2.0 * g_0_yz_x_yz[i] * c_exps[i] + 4.0 * g_zz_yz_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_yz_0_zz[i] = -2.0 * g_0_yz_x_zz[i] * c_exps[i] + 4.0 * g_zz_yz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (750-756)

    #pragma omp simd aligned(g_0_zz_x_xx, g_0_zz_x_xy, g_0_zz_x_xz, g_0_zz_x_yy, g_0_zz_x_yz, g_0_zz_x_zz, g_z_0_x_0_z_zz_0_xx, g_z_0_x_0_z_zz_0_xy, g_z_0_x_0_z_zz_0_xz, g_z_0_x_0_z_zz_0_yy, g_z_0_x_0_z_zz_0_yz, g_z_0_x_0_z_zz_0_zz, g_zz_zz_x_xx, g_zz_zz_x_xy, g_zz_zz_x_xz, g_zz_zz_x_yy, g_zz_zz_x_yz, g_zz_zz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_zz_0_xx[i] = -2.0 * g_0_zz_x_xx[i] * c_exps[i] + 4.0 * g_zz_zz_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_zz_0_xy[i] = -2.0 * g_0_zz_x_xy[i] * c_exps[i] + 4.0 * g_zz_zz_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_zz_0_xz[i] = -2.0 * g_0_zz_x_xz[i] * c_exps[i] + 4.0 * g_zz_zz_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_zz_0_yy[i] = -2.0 * g_0_zz_x_yy[i] * c_exps[i] + 4.0 * g_zz_zz_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_zz_0_yz[i] = -2.0 * g_0_zz_x_yz[i] * c_exps[i] + 4.0 * g_zz_zz_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_zz_0_zz[i] = -2.0 * g_0_zz_x_zz[i] * c_exps[i] + 4.0 * g_zz_zz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (756-762)

    #pragma omp simd aligned(g_xz_xx_y_xx, g_xz_xx_y_xy, g_xz_xx_y_xz, g_xz_xx_y_yy, g_xz_xx_y_yz, g_xz_xx_y_zz, g_z_0_y_0_x_xx_0_xx, g_z_0_y_0_x_xx_0_xy, g_z_0_y_0_x_xx_0_xz, g_z_0_y_0_x_xx_0_yy, g_z_0_y_0_x_xx_0_yz, g_z_0_y_0_x_xx_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_xx_0_xx[i] = 4.0 * g_xz_xx_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xx_0_xy[i] = 4.0 * g_xz_xx_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xx_0_xz[i] = 4.0 * g_xz_xx_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xx_0_yy[i] = 4.0 * g_xz_xx_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xx_0_yz[i] = 4.0 * g_xz_xx_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xx_0_zz[i] = 4.0 * g_xz_xx_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (762-768)

    #pragma omp simd aligned(g_xz_xy_y_xx, g_xz_xy_y_xy, g_xz_xy_y_xz, g_xz_xy_y_yy, g_xz_xy_y_yz, g_xz_xy_y_zz, g_z_0_y_0_x_xy_0_xx, g_z_0_y_0_x_xy_0_xy, g_z_0_y_0_x_xy_0_xz, g_z_0_y_0_x_xy_0_yy, g_z_0_y_0_x_xy_0_yz, g_z_0_y_0_x_xy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_xy_0_xx[i] = 4.0 * g_xz_xy_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xy_0_xy[i] = 4.0 * g_xz_xy_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xy_0_xz[i] = 4.0 * g_xz_xy_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xy_0_yy[i] = 4.0 * g_xz_xy_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xy_0_yz[i] = 4.0 * g_xz_xy_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xy_0_zz[i] = 4.0 * g_xz_xy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (768-774)

    #pragma omp simd aligned(g_xz_xz_y_xx, g_xz_xz_y_xy, g_xz_xz_y_xz, g_xz_xz_y_yy, g_xz_xz_y_yz, g_xz_xz_y_zz, g_z_0_y_0_x_xz_0_xx, g_z_0_y_0_x_xz_0_xy, g_z_0_y_0_x_xz_0_xz, g_z_0_y_0_x_xz_0_yy, g_z_0_y_0_x_xz_0_yz, g_z_0_y_0_x_xz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_xz_0_xx[i] = 4.0 * g_xz_xz_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xz_0_xy[i] = 4.0 * g_xz_xz_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xz_0_xz[i] = 4.0 * g_xz_xz_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xz_0_yy[i] = 4.0 * g_xz_xz_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xz_0_yz[i] = 4.0 * g_xz_xz_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xz_0_zz[i] = 4.0 * g_xz_xz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (774-780)

    #pragma omp simd aligned(g_xz_yy_y_xx, g_xz_yy_y_xy, g_xz_yy_y_xz, g_xz_yy_y_yy, g_xz_yy_y_yz, g_xz_yy_y_zz, g_z_0_y_0_x_yy_0_xx, g_z_0_y_0_x_yy_0_xy, g_z_0_y_0_x_yy_0_xz, g_z_0_y_0_x_yy_0_yy, g_z_0_y_0_x_yy_0_yz, g_z_0_y_0_x_yy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_yy_0_xx[i] = 4.0 * g_xz_yy_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_yy_0_xy[i] = 4.0 * g_xz_yy_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_yy_0_xz[i] = 4.0 * g_xz_yy_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_yy_0_yy[i] = 4.0 * g_xz_yy_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_yy_0_yz[i] = 4.0 * g_xz_yy_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_yy_0_zz[i] = 4.0 * g_xz_yy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (780-786)

    #pragma omp simd aligned(g_xz_yz_y_xx, g_xz_yz_y_xy, g_xz_yz_y_xz, g_xz_yz_y_yy, g_xz_yz_y_yz, g_xz_yz_y_zz, g_z_0_y_0_x_yz_0_xx, g_z_0_y_0_x_yz_0_xy, g_z_0_y_0_x_yz_0_xz, g_z_0_y_0_x_yz_0_yy, g_z_0_y_0_x_yz_0_yz, g_z_0_y_0_x_yz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_yz_0_xx[i] = 4.0 * g_xz_yz_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_yz_0_xy[i] = 4.0 * g_xz_yz_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_yz_0_xz[i] = 4.0 * g_xz_yz_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_yz_0_yy[i] = 4.0 * g_xz_yz_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_yz_0_yz[i] = 4.0 * g_xz_yz_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_yz_0_zz[i] = 4.0 * g_xz_yz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (786-792)

    #pragma omp simd aligned(g_xz_zz_y_xx, g_xz_zz_y_xy, g_xz_zz_y_xz, g_xz_zz_y_yy, g_xz_zz_y_yz, g_xz_zz_y_zz, g_z_0_y_0_x_zz_0_xx, g_z_0_y_0_x_zz_0_xy, g_z_0_y_0_x_zz_0_xz, g_z_0_y_0_x_zz_0_yy, g_z_0_y_0_x_zz_0_yz, g_z_0_y_0_x_zz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_zz_0_xx[i] = 4.0 * g_xz_zz_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_zz_0_xy[i] = 4.0 * g_xz_zz_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_zz_0_xz[i] = 4.0 * g_xz_zz_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_zz_0_yy[i] = 4.0 * g_xz_zz_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_zz_0_yz[i] = 4.0 * g_xz_zz_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_zz_0_zz[i] = 4.0 * g_xz_zz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (792-798)

    #pragma omp simd aligned(g_yz_xx_y_xx, g_yz_xx_y_xy, g_yz_xx_y_xz, g_yz_xx_y_yy, g_yz_xx_y_yz, g_yz_xx_y_zz, g_z_0_y_0_y_xx_0_xx, g_z_0_y_0_y_xx_0_xy, g_z_0_y_0_y_xx_0_xz, g_z_0_y_0_y_xx_0_yy, g_z_0_y_0_y_xx_0_yz, g_z_0_y_0_y_xx_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_xx_0_xx[i] = 4.0 * g_yz_xx_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xx_0_xy[i] = 4.0 * g_yz_xx_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xx_0_xz[i] = 4.0 * g_yz_xx_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xx_0_yy[i] = 4.0 * g_yz_xx_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xx_0_yz[i] = 4.0 * g_yz_xx_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xx_0_zz[i] = 4.0 * g_yz_xx_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (798-804)

    #pragma omp simd aligned(g_yz_xy_y_xx, g_yz_xy_y_xy, g_yz_xy_y_xz, g_yz_xy_y_yy, g_yz_xy_y_yz, g_yz_xy_y_zz, g_z_0_y_0_y_xy_0_xx, g_z_0_y_0_y_xy_0_xy, g_z_0_y_0_y_xy_0_xz, g_z_0_y_0_y_xy_0_yy, g_z_0_y_0_y_xy_0_yz, g_z_0_y_0_y_xy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_xy_0_xx[i] = 4.0 * g_yz_xy_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xy_0_xy[i] = 4.0 * g_yz_xy_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xy_0_xz[i] = 4.0 * g_yz_xy_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xy_0_yy[i] = 4.0 * g_yz_xy_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xy_0_yz[i] = 4.0 * g_yz_xy_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xy_0_zz[i] = 4.0 * g_yz_xy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (804-810)

    #pragma omp simd aligned(g_yz_xz_y_xx, g_yz_xz_y_xy, g_yz_xz_y_xz, g_yz_xz_y_yy, g_yz_xz_y_yz, g_yz_xz_y_zz, g_z_0_y_0_y_xz_0_xx, g_z_0_y_0_y_xz_0_xy, g_z_0_y_0_y_xz_0_xz, g_z_0_y_0_y_xz_0_yy, g_z_0_y_0_y_xz_0_yz, g_z_0_y_0_y_xz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_xz_0_xx[i] = 4.0 * g_yz_xz_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xz_0_xy[i] = 4.0 * g_yz_xz_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xz_0_xz[i] = 4.0 * g_yz_xz_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xz_0_yy[i] = 4.0 * g_yz_xz_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xz_0_yz[i] = 4.0 * g_yz_xz_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xz_0_zz[i] = 4.0 * g_yz_xz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (810-816)

    #pragma omp simd aligned(g_yz_yy_y_xx, g_yz_yy_y_xy, g_yz_yy_y_xz, g_yz_yy_y_yy, g_yz_yy_y_yz, g_yz_yy_y_zz, g_z_0_y_0_y_yy_0_xx, g_z_0_y_0_y_yy_0_xy, g_z_0_y_0_y_yy_0_xz, g_z_0_y_0_y_yy_0_yy, g_z_0_y_0_y_yy_0_yz, g_z_0_y_0_y_yy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_yy_0_xx[i] = 4.0 * g_yz_yy_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_yy_0_xy[i] = 4.0 * g_yz_yy_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_yy_0_xz[i] = 4.0 * g_yz_yy_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_yy_0_yy[i] = 4.0 * g_yz_yy_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_yy_0_yz[i] = 4.0 * g_yz_yy_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_yy_0_zz[i] = 4.0 * g_yz_yy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (816-822)

    #pragma omp simd aligned(g_yz_yz_y_xx, g_yz_yz_y_xy, g_yz_yz_y_xz, g_yz_yz_y_yy, g_yz_yz_y_yz, g_yz_yz_y_zz, g_z_0_y_0_y_yz_0_xx, g_z_0_y_0_y_yz_0_xy, g_z_0_y_0_y_yz_0_xz, g_z_0_y_0_y_yz_0_yy, g_z_0_y_0_y_yz_0_yz, g_z_0_y_0_y_yz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_yz_0_xx[i] = 4.0 * g_yz_yz_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_yz_0_xy[i] = 4.0 * g_yz_yz_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_yz_0_xz[i] = 4.0 * g_yz_yz_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_yz_0_yy[i] = 4.0 * g_yz_yz_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_yz_0_yz[i] = 4.0 * g_yz_yz_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_yz_0_zz[i] = 4.0 * g_yz_yz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (822-828)

    #pragma omp simd aligned(g_yz_zz_y_xx, g_yz_zz_y_xy, g_yz_zz_y_xz, g_yz_zz_y_yy, g_yz_zz_y_yz, g_yz_zz_y_zz, g_z_0_y_0_y_zz_0_xx, g_z_0_y_0_y_zz_0_xy, g_z_0_y_0_y_zz_0_xz, g_z_0_y_0_y_zz_0_yy, g_z_0_y_0_y_zz_0_yz, g_z_0_y_0_y_zz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_zz_0_xx[i] = 4.0 * g_yz_zz_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_zz_0_xy[i] = 4.0 * g_yz_zz_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_zz_0_xz[i] = 4.0 * g_yz_zz_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_zz_0_yy[i] = 4.0 * g_yz_zz_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_zz_0_yz[i] = 4.0 * g_yz_zz_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_zz_0_zz[i] = 4.0 * g_yz_zz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (828-834)

    #pragma omp simd aligned(g_0_xx_y_xx, g_0_xx_y_xy, g_0_xx_y_xz, g_0_xx_y_yy, g_0_xx_y_yz, g_0_xx_y_zz, g_z_0_y_0_z_xx_0_xx, g_z_0_y_0_z_xx_0_xy, g_z_0_y_0_z_xx_0_xz, g_z_0_y_0_z_xx_0_yy, g_z_0_y_0_z_xx_0_yz, g_z_0_y_0_z_xx_0_zz, g_zz_xx_y_xx, g_zz_xx_y_xy, g_zz_xx_y_xz, g_zz_xx_y_yy, g_zz_xx_y_yz, g_zz_xx_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_xx_0_xx[i] = -2.0 * g_0_xx_y_xx[i] * c_exps[i] + 4.0 * g_zz_xx_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xx_0_xy[i] = -2.0 * g_0_xx_y_xy[i] * c_exps[i] + 4.0 * g_zz_xx_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xx_0_xz[i] = -2.0 * g_0_xx_y_xz[i] * c_exps[i] + 4.0 * g_zz_xx_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xx_0_yy[i] = -2.0 * g_0_xx_y_yy[i] * c_exps[i] + 4.0 * g_zz_xx_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xx_0_yz[i] = -2.0 * g_0_xx_y_yz[i] * c_exps[i] + 4.0 * g_zz_xx_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xx_0_zz[i] = -2.0 * g_0_xx_y_zz[i] * c_exps[i] + 4.0 * g_zz_xx_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (834-840)

    #pragma omp simd aligned(g_0_xy_y_xx, g_0_xy_y_xy, g_0_xy_y_xz, g_0_xy_y_yy, g_0_xy_y_yz, g_0_xy_y_zz, g_z_0_y_0_z_xy_0_xx, g_z_0_y_0_z_xy_0_xy, g_z_0_y_0_z_xy_0_xz, g_z_0_y_0_z_xy_0_yy, g_z_0_y_0_z_xy_0_yz, g_z_0_y_0_z_xy_0_zz, g_zz_xy_y_xx, g_zz_xy_y_xy, g_zz_xy_y_xz, g_zz_xy_y_yy, g_zz_xy_y_yz, g_zz_xy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_xy_0_xx[i] = -2.0 * g_0_xy_y_xx[i] * c_exps[i] + 4.0 * g_zz_xy_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xy_0_xy[i] = -2.0 * g_0_xy_y_xy[i] * c_exps[i] + 4.0 * g_zz_xy_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xy_0_xz[i] = -2.0 * g_0_xy_y_xz[i] * c_exps[i] + 4.0 * g_zz_xy_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xy_0_yy[i] = -2.0 * g_0_xy_y_yy[i] * c_exps[i] + 4.0 * g_zz_xy_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xy_0_yz[i] = -2.0 * g_0_xy_y_yz[i] * c_exps[i] + 4.0 * g_zz_xy_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xy_0_zz[i] = -2.0 * g_0_xy_y_zz[i] * c_exps[i] + 4.0 * g_zz_xy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (840-846)

    #pragma omp simd aligned(g_0_xz_y_xx, g_0_xz_y_xy, g_0_xz_y_xz, g_0_xz_y_yy, g_0_xz_y_yz, g_0_xz_y_zz, g_z_0_y_0_z_xz_0_xx, g_z_0_y_0_z_xz_0_xy, g_z_0_y_0_z_xz_0_xz, g_z_0_y_0_z_xz_0_yy, g_z_0_y_0_z_xz_0_yz, g_z_0_y_0_z_xz_0_zz, g_zz_xz_y_xx, g_zz_xz_y_xy, g_zz_xz_y_xz, g_zz_xz_y_yy, g_zz_xz_y_yz, g_zz_xz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_xz_0_xx[i] = -2.0 * g_0_xz_y_xx[i] * c_exps[i] + 4.0 * g_zz_xz_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xz_0_xy[i] = -2.0 * g_0_xz_y_xy[i] * c_exps[i] + 4.0 * g_zz_xz_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xz_0_xz[i] = -2.0 * g_0_xz_y_xz[i] * c_exps[i] + 4.0 * g_zz_xz_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xz_0_yy[i] = -2.0 * g_0_xz_y_yy[i] * c_exps[i] + 4.0 * g_zz_xz_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xz_0_yz[i] = -2.0 * g_0_xz_y_yz[i] * c_exps[i] + 4.0 * g_zz_xz_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xz_0_zz[i] = -2.0 * g_0_xz_y_zz[i] * c_exps[i] + 4.0 * g_zz_xz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (846-852)

    #pragma omp simd aligned(g_0_yy_y_xx, g_0_yy_y_xy, g_0_yy_y_xz, g_0_yy_y_yy, g_0_yy_y_yz, g_0_yy_y_zz, g_z_0_y_0_z_yy_0_xx, g_z_0_y_0_z_yy_0_xy, g_z_0_y_0_z_yy_0_xz, g_z_0_y_0_z_yy_0_yy, g_z_0_y_0_z_yy_0_yz, g_z_0_y_0_z_yy_0_zz, g_zz_yy_y_xx, g_zz_yy_y_xy, g_zz_yy_y_xz, g_zz_yy_y_yy, g_zz_yy_y_yz, g_zz_yy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_yy_0_xx[i] = -2.0 * g_0_yy_y_xx[i] * c_exps[i] + 4.0 * g_zz_yy_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_yy_0_xy[i] = -2.0 * g_0_yy_y_xy[i] * c_exps[i] + 4.0 * g_zz_yy_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_yy_0_xz[i] = -2.0 * g_0_yy_y_xz[i] * c_exps[i] + 4.0 * g_zz_yy_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_yy_0_yy[i] = -2.0 * g_0_yy_y_yy[i] * c_exps[i] + 4.0 * g_zz_yy_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_yy_0_yz[i] = -2.0 * g_0_yy_y_yz[i] * c_exps[i] + 4.0 * g_zz_yy_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_yy_0_zz[i] = -2.0 * g_0_yy_y_zz[i] * c_exps[i] + 4.0 * g_zz_yy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (852-858)

    #pragma omp simd aligned(g_0_yz_y_xx, g_0_yz_y_xy, g_0_yz_y_xz, g_0_yz_y_yy, g_0_yz_y_yz, g_0_yz_y_zz, g_z_0_y_0_z_yz_0_xx, g_z_0_y_0_z_yz_0_xy, g_z_0_y_0_z_yz_0_xz, g_z_0_y_0_z_yz_0_yy, g_z_0_y_0_z_yz_0_yz, g_z_0_y_0_z_yz_0_zz, g_zz_yz_y_xx, g_zz_yz_y_xy, g_zz_yz_y_xz, g_zz_yz_y_yy, g_zz_yz_y_yz, g_zz_yz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_yz_0_xx[i] = -2.0 * g_0_yz_y_xx[i] * c_exps[i] + 4.0 * g_zz_yz_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_yz_0_xy[i] = -2.0 * g_0_yz_y_xy[i] * c_exps[i] + 4.0 * g_zz_yz_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_yz_0_xz[i] = -2.0 * g_0_yz_y_xz[i] * c_exps[i] + 4.0 * g_zz_yz_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_yz_0_yy[i] = -2.0 * g_0_yz_y_yy[i] * c_exps[i] + 4.0 * g_zz_yz_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_yz_0_yz[i] = -2.0 * g_0_yz_y_yz[i] * c_exps[i] + 4.0 * g_zz_yz_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_yz_0_zz[i] = -2.0 * g_0_yz_y_zz[i] * c_exps[i] + 4.0 * g_zz_yz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (858-864)

    #pragma omp simd aligned(g_0_zz_y_xx, g_0_zz_y_xy, g_0_zz_y_xz, g_0_zz_y_yy, g_0_zz_y_yz, g_0_zz_y_zz, g_z_0_y_0_z_zz_0_xx, g_z_0_y_0_z_zz_0_xy, g_z_0_y_0_z_zz_0_xz, g_z_0_y_0_z_zz_0_yy, g_z_0_y_0_z_zz_0_yz, g_z_0_y_0_z_zz_0_zz, g_zz_zz_y_xx, g_zz_zz_y_xy, g_zz_zz_y_xz, g_zz_zz_y_yy, g_zz_zz_y_yz, g_zz_zz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_zz_0_xx[i] = -2.0 * g_0_zz_y_xx[i] * c_exps[i] + 4.0 * g_zz_zz_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_zz_0_xy[i] = -2.0 * g_0_zz_y_xy[i] * c_exps[i] + 4.0 * g_zz_zz_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_zz_0_xz[i] = -2.0 * g_0_zz_y_xz[i] * c_exps[i] + 4.0 * g_zz_zz_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_zz_0_yy[i] = -2.0 * g_0_zz_y_yy[i] * c_exps[i] + 4.0 * g_zz_zz_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_zz_0_yz[i] = -2.0 * g_0_zz_y_yz[i] * c_exps[i] + 4.0 * g_zz_zz_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_zz_0_zz[i] = -2.0 * g_0_zz_y_zz[i] * c_exps[i] + 4.0 * g_zz_zz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (864-870)

    #pragma omp simd aligned(g_xz_xx_z_xx, g_xz_xx_z_xy, g_xz_xx_z_xz, g_xz_xx_z_yy, g_xz_xx_z_yz, g_xz_xx_z_zz, g_z_0_z_0_x_xx_0_xx, g_z_0_z_0_x_xx_0_xy, g_z_0_z_0_x_xx_0_xz, g_z_0_z_0_x_xx_0_yy, g_z_0_z_0_x_xx_0_yz, g_z_0_z_0_x_xx_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_xx_0_xx[i] = 4.0 * g_xz_xx_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xx_0_xy[i] = 4.0 * g_xz_xx_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xx_0_xz[i] = 4.0 * g_xz_xx_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xx_0_yy[i] = 4.0 * g_xz_xx_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xx_0_yz[i] = 4.0 * g_xz_xx_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xx_0_zz[i] = 4.0 * g_xz_xx_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (870-876)

    #pragma omp simd aligned(g_xz_xy_z_xx, g_xz_xy_z_xy, g_xz_xy_z_xz, g_xz_xy_z_yy, g_xz_xy_z_yz, g_xz_xy_z_zz, g_z_0_z_0_x_xy_0_xx, g_z_0_z_0_x_xy_0_xy, g_z_0_z_0_x_xy_0_xz, g_z_0_z_0_x_xy_0_yy, g_z_0_z_0_x_xy_0_yz, g_z_0_z_0_x_xy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_xy_0_xx[i] = 4.0 * g_xz_xy_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xy_0_xy[i] = 4.0 * g_xz_xy_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xy_0_xz[i] = 4.0 * g_xz_xy_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xy_0_yy[i] = 4.0 * g_xz_xy_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xy_0_yz[i] = 4.0 * g_xz_xy_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xy_0_zz[i] = 4.0 * g_xz_xy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (876-882)

    #pragma omp simd aligned(g_xz_xz_z_xx, g_xz_xz_z_xy, g_xz_xz_z_xz, g_xz_xz_z_yy, g_xz_xz_z_yz, g_xz_xz_z_zz, g_z_0_z_0_x_xz_0_xx, g_z_0_z_0_x_xz_0_xy, g_z_0_z_0_x_xz_0_xz, g_z_0_z_0_x_xz_0_yy, g_z_0_z_0_x_xz_0_yz, g_z_0_z_0_x_xz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_xz_0_xx[i] = 4.0 * g_xz_xz_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xz_0_xy[i] = 4.0 * g_xz_xz_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xz_0_xz[i] = 4.0 * g_xz_xz_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xz_0_yy[i] = 4.0 * g_xz_xz_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xz_0_yz[i] = 4.0 * g_xz_xz_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xz_0_zz[i] = 4.0 * g_xz_xz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (882-888)

    #pragma omp simd aligned(g_xz_yy_z_xx, g_xz_yy_z_xy, g_xz_yy_z_xz, g_xz_yy_z_yy, g_xz_yy_z_yz, g_xz_yy_z_zz, g_z_0_z_0_x_yy_0_xx, g_z_0_z_0_x_yy_0_xy, g_z_0_z_0_x_yy_0_xz, g_z_0_z_0_x_yy_0_yy, g_z_0_z_0_x_yy_0_yz, g_z_0_z_0_x_yy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_yy_0_xx[i] = 4.0 * g_xz_yy_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_yy_0_xy[i] = 4.0 * g_xz_yy_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_yy_0_xz[i] = 4.0 * g_xz_yy_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_yy_0_yy[i] = 4.0 * g_xz_yy_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_yy_0_yz[i] = 4.0 * g_xz_yy_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_yy_0_zz[i] = 4.0 * g_xz_yy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (888-894)

    #pragma omp simd aligned(g_xz_yz_z_xx, g_xz_yz_z_xy, g_xz_yz_z_xz, g_xz_yz_z_yy, g_xz_yz_z_yz, g_xz_yz_z_zz, g_z_0_z_0_x_yz_0_xx, g_z_0_z_0_x_yz_0_xy, g_z_0_z_0_x_yz_0_xz, g_z_0_z_0_x_yz_0_yy, g_z_0_z_0_x_yz_0_yz, g_z_0_z_0_x_yz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_yz_0_xx[i] = 4.0 * g_xz_yz_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_yz_0_xy[i] = 4.0 * g_xz_yz_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_yz_0_xz[i] = 4.0 * g_xz_yz_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_yz_0_yy[i] = 4.0 * g_xz_yz_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_yz_0_yz[i] = 4.0 * g_xz_yz_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_yz_0_zz[i] = 4.0 * g_xz_yz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (894-900)

    #pragma omp simd aligned(g_xz_zz_z_xx, g_xz_zz_z_xy, g_xz_zz_z_xz, g_xz_zz_z_yy, g_xz_zz_z_yz, g_xz_zz_z_zz, g_z_0_z_0_x_zz_0_xx, g_z_0_z_0_x_zz_0_xy, g_z_0_z_0_x_zz_0_xz, g_z_0_z_0_x_zz_0_yy, g_z_0_z_0_x_zz_0_yz, g_z_0_z_0_x_zz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_zz_0_xx[i] = 4.0 * g_xz_zz_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_zz_0_xy[i] = 4.0 * g_xz_zz_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_zz_0_xz[i] = 4.0 * g_xz_zz_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_zz_0_yy[i] = 4.0 * g_xz_zz_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_zz_0_yz[i] = 4.0 * g_xz_zz_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_zz_0_zz[i] = 4.0 * g_xz_zz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (900-906)

    #pragma omp simd aligned(g_yz_xx_z_xx, g_yz_xx_z_xy, g_yz_xx_z_xz, g_yz_xx_z_yy, g_yz_xx_z_yz, g_yz_xx_z_zz, g_z_0_z_0_y_xx_0_xx, g_z_0_z_0_y_xx_0_xy, g_z_0_z_0_y_xx_0_xz, g_z_0_z_0_y_xx_0_yy, g_z_0_z_0_y_xx_0_yz, g_z_0_z_0_y_xx_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_xx_0_xx[i] = 4.0 * g_yz_xx_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xx_0_xy[i] = 4.0 * g_yz_xx_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xx_0_xz[i] = 4.0 * g_yz_xx_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xx_0_yy[i] = 4.0 * g_yz_xx_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xx_0_yz[i] = 4.0 * g_yz_xx_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xx_0_zz[i] = 4.0 * g_yz_xx_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (906-912)

    #pragma omp simd aligned(g_yz_xy_z_xx, g_yz_xy_z_xy, g_yz_xy_z_xz, g_yz_xy_z_yy, g_yz_xy_z_yz, g_yz_xy_z_zz, g_z_0_z_0_y_xy_0_xx, g_z_0_z_0_y_xy_0_xy, g_z_0_z_0_y_xy_0_xz, g_z_0_z_0_y_xy_0_yy, g_z_0_z_0_y_xy_0_yz, g_z_0_z_0_y_xy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_xy_0_xx[i] = 4.0 * g_yz_xy_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xy_0_xy[i] = 4.0 * g_yz_xy_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xy_0_xz[i] = 4.0 * g_yz_xy_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xy_0_yy[i] = 4.0 * g_yz_xy_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xy_0_yz[i] = 4.0 * g_yz_xy_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xy_0_zz[i] = 4.0 * g_yz_xy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (912-918)

    #pragma omp simd aligned(g_yz_xz_z_xx, g_yz_xz_z_xy, g_yz_xz_z_xz, g_yz_xz_z_yy, g_yz_xz_z_yz, g_yz_xz_z_zz, g_z_0_z_0_y_xz_0_xx, g_z_0_z_0_y_xz_0_xy, g_z_0_z_0_y_xz_0_xz, g_z_0_z_0_y_xz_0_yy, g_z_0_z_0_y_xz_0_yz, g_z_0_z_0_y_xz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_xz_0_xx[i] = 4.0 * g_yz_xz_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xz_0_xy[i] = 4.0 * g_yz_xz_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xz_0_xz[i] = 4.0 * g_yz_xz_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xz_0_yy[i] = 4.0 * g_yz_xz_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xz_0_yz[i] = 4.0 * g_yz_xz_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xz_0_zz[i] = 4.0 * g_yz_xz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (918-924)

    #pragma omp simd aligned(g_yz_yy_z_xx, g_yz_yy_z_xy, g_yz_yy_z_xz, g_yz_yy_z_yy, g_yz_yy_z_yz, g_yz_yy_z_zz, g_z_0_z_0_y_yy_0_xx, g_z_0_z_0_y_yy_0_xy, g_z_0_z_0_y_yy_0_xz, g_z_0_z_0_y_yy_0_yy, g_z_0_z_0_y_yy_0_yz, g_z_0_z_0_y_yy_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_yy_0_xx[i] = 4.0 * g_yz_yy_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_yy_0_xy[i] = 4.0 * g_yz_yy_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_yy_0_xz[i] = 4.0 * g_yz_yy_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_yy_0_yy[i] = 4.0 * g_yz_yy_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_yy_0_yz[i] = 4.0 * g_yz_yy_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_yy_0_zz[i] = 4.0 * g_yz_yy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (924-930)

    #pragma omp simd aligned(g_yz_yz_z_xx, g_yz_yz_z_xy, g_yz_yz_z_xz, g_yz_yz_z_yy, g_yz_yz_z_yz, g_yz_yz_z_zz, g_z_0_z_0_y_yz_0_xx, g_z_0_z_0_y_yz_0_xy, g_z_0_z_0_y_yz_0_xz, g_z_0_z_0_y_yz_0_yy, g_z_0_z_0_y_yz_0_yz, g_z_0_z_0_y_yz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_yz_0_xx[i] = 4.0 * g_yz_yz_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_yz_0_xy[i] = 4.0 * g_yz_yz_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_yz_0_xz[i] = 4.0 * g_yz_yz_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_yz_0_yy[i] = 4.0 * g_yz_yz_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_yz_0_yz[i] = 4.0 * g_yz_yz_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_yz_0_zz[i] = 4.0 * g_yz_yz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (930-936)

    #pragma omp simd aligned(g_yz_zz_z_xx, g_yz_zz_z_xy, g_yz_zz_z_xz, g_yz_zz_z_yy, g_yz_zz_z_yz, g_yz_zz_z_zz, g_z_0_z_0_y_zz_0_xx, g_z_0_z_0_y_zz_0_xy, g_z_0_z_0_y_zz_0_xz, g_z_0_z_0_y_zz_0_yy, g_z_0_z_0_y_zz_0_yz, g_z_0_z_0_y_zz_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_zz_0_xx[i] = 4.0 * g_yz_zz_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_zz_0_xy[i] = 4.0 * g_yz_zz_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_zz_0_xz[i] = 4.0 * g_yz_zz_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_zz_0_yy[i] = 4.0 * g_yz_zz_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_zz_0_yz[i] = 4.0 * g_yz_zz_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_zz_0_zz[i] = 4.0 * g_yz_zz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (936-942)

    #pragma omp simd aligned(g_0_xx_z_xx, g_0_xx_z_xy, g_0_xx_z_xz, g_0_xx_z_yy, g_0_xx_z_yz, g_0_xx_z_zz, g_z_0_z_0_z_xx_0_xx, g_z_0_z_0_z_xx_0_xy, g_z_0_z_0_z_xx_0_xz, g_z_0_z_0_z_xx_0_yy, g_z_0_z_0_z_xx_0_yz, g_z_0_z_0_z_xx_0_zz, g_zz_xx_z_xx, g_zz_xx_z_xy, g_zz_xx_z_xz, g_zz_xx_z_yy, g_zz_xx_z_yz, g_zz_xx_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_xx_0_xx[i] = -2.0 * g_0_xx_z_xx[i] * c_exps[i] + 4.0 * g_zz_xx_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xx_0_xy[i] = -2.0 * g_0_xx_z_xy[i] * c_exps[i] + 4.0 * g_zz_xx_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xx_0_xz[i] = -2.0 * g_0_xx_z_xz[i] * c_exps[i] + 4.0 * g_zz_xx_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xx_0_yy[i] = -2.0 * g_0_xx_z_yy[i] * c_exps[i] + 4.0 * g_zz_xx_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xx_0_yz[i] = -2.0 * g_0_xx_z_yz[i] * c_exps[i] + 4.0 * g_zz_xx_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xx_0_zz[i] = -2.0 * g_0_xx_z_zz[i] * c_exps[i] + 4.0 * g_zz_xx_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (942-948)

    #pragma omp simd aligned(g_0_xy_z_xx, g_0_xy_z_xy, g_0_xy_z_xz, g_0_xy_z_yy, g_0_xy_z_yz, g_0_xy_z_zz, g_z_0_z_0_z_xy_0_xx, g_z_0_z_0_z_xy_0_xy, g_z_0_z_0_z_xy_0_xz, g_z_0_z_0_z_xy_0_yy, g_z_0_z_0_z_xy_0_yz, g_z_0_z_0_z_xy_0_zz, g_zz_xy_z_xx, g_zz_xy_z_xy, g_zz_xy_z_xz, g_zz_xy_z_yy, g_zz_xy_z_yz, g_zz_xy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_xy_0_xx[i] = -2.0 * g_0_xy_z_xx[i] * c_exps[i] + 4.0 * g_zz_xy_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xy_0_xy[i] = -2.0 * g_0_xy_z_xy[i] * c_exps[i] + 4.0 * g_zz_xy_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xy_0_xz[i] = -2.0 * g_0_xy_z_xz[i] * c_exps[i] + 4.0 * g_zz_xy_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xy_0_yy[i] = -2.0 * g_0_xy_z_yy[i] * c_exps[i] + 4.0 * g_zz_xy_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xy_0_yz[i] = -2.0 * g_0_xy_z_yz[i] * c_exps[i] + 4.0 * g_zz_xy_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xy_0_zz[i] = -2.0 * g_0_xy_z_zz[i] * c_exps[i] + 4.0 * g_zz_xy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (948-954)

    #pragma omp simd aligned(g_0_xz_z_xx, g_0_xz_z_xy, g_0_xz_z_xz, g_0_xz_z_yy, g_0_xz_z_yz, g_0_xz_z_zz, g_z_0_z_0_z_xz_0_xx, g_z_0_z_0_z_xz_0_xy, g_z_0_z_0_z_xz_0_xz, g_z_0_z_0_z_xz_0_yy, g_z_0_z_0_z_xz_0_yz, g_z_0_z_0_z_xz_0_zz, g_zz_xz_z_xx, g_zz_xz_z_xy, g_zz_xz_z_xz, g_zz_xz_z_yy, g_zz_xz_z_yz, g_zz_xz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_xz_0_xx[i] = -2.0 * g_0_xz_z_xx[i] * c_exps[i] + 4.0 * g_zz_xz_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xz_0_xy[i] = -2.0 * g_0_xz_z_xy[i] * c_exps[i] + 4.0 * g_zz_xz_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xz_0_xz[i] = -2.0 * g_0_xz_z_xz[i] * c_exps[i] + 4.0 * g_zz_xz_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xz_0_yy[i] = -2.0 * g_0_xz_z_yy[i] * c_exps[i] + 4.0 * g_zz_xz_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xz_0_yz[i] = -2.0 * g_0_xz_z_yz[i] * c_exps[i] + 4.0 * g_zz_xz_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xz_0_zz[i] = -2.0 * g_0_xz_z_zz[i] * c_exps[i] + 4.0 * g_zz_xz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (954-960)

    #pragma omp simd aligned(g_0_yy_z_xx, g_0_yy_z_xy, g_0_yy_z_xz, g_0_yy_z_yy, g_0_yy_z_yz, g_0_yy_z_zz, g_z_0_z_0_z_yy_0_xx, g_z_0_z_0_z_yy_0_xy, g_z_0_z_0_z_yy_0_xz, g_z_0_z_0_z_yy_0_yy, g_z_0_z_0_z_yy_0_yz, g_z_0_z_0_z_yy_0_zz, g_zz_yy_z_xx, g_zz_yy_z_xy, g_zz_yy_z_xz, g_zz_yy_z_yy, g_zz_yy_z_yz, g_zz_yy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_yy_0_xx[i] = -2.0 * g_0_yy_z_xx[i] * c_exps[i] + 4.0 * g_zz_yy_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_yy_0_xy[i] = -2.0 * g_0_yy_z_xy[i] * c_exps[i] + 4.0 * g_zz_yy_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_yy_0_xz[i] = -2.0 * g_0_yy_z_xz[i] * c_exps[i] + 4.0 * g_zz_yy_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_yy_0_yy[i] = -2.0 * g_0_yy_z_yy[i] * c_exps[i] + 4.0 * g_zz_yy_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_yy_0_yz[i] = -2.0 * g_0_yy_z_yz[i] * c_exps[i] + 4.0 * g_zz_yy_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_yy_0_zz[i] = -2.0 * g_0_yy_z_zz[i] * c_exps[i] + 4.0 * g_zz_yy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (960-966)

    #pragma omp simd aligned(g_0_yz_z_xx, g_0_yz_z_xy, g_0_yz_z_xz, g_0_yz_z_yy, g_0_yz_z_yz, g_0_yz_z_zz, g_z_0_z_0_z_yz_0_xx, g_z_0_z_0_z_yz_0_xy, g_z_0_z_0_z_yz_0_xz, g_z_0_z_0_z_yz_0_yy, g_z_0_z_0_z_yz_0_yz, g_z_0_z_0_z_yz_0_zz, g_zz_yz_z_xx, g_zz_yz_z_xy, g_zz_yz_z_xz, g_zz_yz_z_yy, g_zz_yz_z_yz, g_zz_yz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_yz_0_xx[i] = -2.0 * g_0_yz_z_xx[i] * c_exps[i] + 4.0 * g_zz_yz_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_yz_0_xy[i] = -2.0 * g_0_yz_z_xy[i] * c_exps[i] + 4.0 * g_zz_yz_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_yz_0_xz[i] = -2.0 * g_0_yz_z_xz[i] * c_exps[i] + 4.0 * g_zz_yz_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_yz_0_yy[i] = -2.0 * g_0_yz_z_yy[i] * c_exps[i] + 4.0 * g_zz_yz_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_yz_0_yz[i] = -2.0 * g_0_yz_z_yz[i] * c_exps[i] + 4.0 * g_zz_yz_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_yz_0_zz[i] = -2.0 * g_0_yz_z_zz[i] * c_exps[i] + 4.0 * g_zz_yz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (966-972)

    #pragma omp simd aligned(g_0_zz_z_xx, g_0_zz_z_xy, g_0_zz_z_xz, g_0_zz_z_yy, g_0_zz_z_yz, g_0_zz_z_zz, g_z_0_z_0_z_zz_0_xx, g_z_0_z_0_z_zz_0_xy, g_z_0_z_0_z_zz_0_xz, g_z_0_z_0_z_zz_0_yy, g_z_0_z_0_z_zz_0_yz, g_z_0_z_0_z_zz_0_zz, g_zz_zz_z_xx, g_zz_zz_z_xy, g_zz_zz_z_xz, g_zz_zz_z_yy, g_zz_zz_z_yz, g_zz_zz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_zz_0_xx[i] = -2.0 * g_0_zz_z_xx[i] * c_exps[i] + 4.0 * g_zz_zz_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_zz_0_xy[i] = -2.0 * g_0_zz_z_xy[i] * c_exps[i] + 4.0 * g_zz_zz_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_zz_0_xz[i] = -2.0 * g_0_zz_z_xz[i] * c_exps[i] + 4.0 * g_zz_zz_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_zz_0_yy[i] = -2.0 * g_0_zz_z_yy[i] * c_exps[i] + 4.0 * g_zz_zz_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_zz_0_yz[i] = -2.0 * g_0_zz_z_yz[i] * c_exps[i] + 4.0 * g_zz_zz_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_zz_0_zz[i] = -2.0 * g_0_zz_z_zz[i] * c_exps[i] + 4.0 * g_zz_zz_z_zz[i] * a_exp * c_exps[i];
    }
}

} // t4c_geom namespace


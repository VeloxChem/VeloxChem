#include "GeomDeriv1000OfScalarForPDPD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1000_pdpd_0(CSimdArray<double>& buffer_1000_pdpd,
                     const CSimdArray<double>& buffer_sdpd,
                     const CSimdArray<double>& buffer_ddpd,
                     const double a_exp) -> void
{
    const auto ndims = buffer_1000_pdpd.number_of_columns();

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

    /// Set up components of integrals buffer : buffer_1000_pdpd

    auto g_x_0_0_0_x_xx_x_xx = buffer_1000_pdpd[0];

    auto g_x_0_0_0_x_xx_x_xy = buffer_1000_pdpd[1];

    auto g_x_0_0_0_x_xx_x_xz = buffer_1000_pdpd[2];

    auto g_x_0_0_0_x_xx_x_yy = buffer_1000_pdpd[3];

    auto g_x_0_0_0_x_xx_x_yz = buffer_1000_pdpd[4];

    auto g_x_0_0_0_x_xx_x_zz = buffer_1000_pdpd[5];

    auto g_x_0_0_0_x_xx_y_xx = buffer_1000_pdpd[6];

    auto g_x_0_0_0_x_xx_y_xy = buffer_1000_pdpd[7];

    auto g_x_0_0_0_x_xx_y_xz = buffer_1000_pdpd[8];

    auto g_x_0_0_0_x_xx_y_yy = buffer_1000_pdpd[9];

    auto g_x_0_0_0_x_xx_y_yz = buffer_1000_pdpd[10];

    auto g_x_0_0_0_x_xx_y_zz = buffer_1000_pdpd[11];

    auto g_x_0_0_0_x_xx_z_xx = buffer_1000_pdpd[12];

    auto g_x_0_0_0_x_xx_z_xy = buffer_1000_pdpd[13];

    auto g_x_0_0_0_x_xx_z_xz = buffer_1000_pdpd[14];

    auto g_x_0_0_0_x_xx_z_yy = buffer_1000_pdpd[15];

    auto g_x_0_0_0_x_xx_z_yz = buffer_1000_pdpd[16];

    auto g_x_0_0_0_x_xx_z_zz = buffer_1000_pdpd[17];

    auto g_x_0_0_0_x_xy_x_xx = buffer_1000_pdpd[18];

    auto g_x_0_0_0_x_xy_x_xy = buffer_1000_pdpd[19];

    auto g_x_0_0_0_x_xy_x_xz = buffer_1000_pdpd[20];

    auto g_x_0_0_0_x_xy_x_yy = buffer_1000_pdpd[21];

    auto g_x_0_0_0_x_xy_x_yz = buffer_1000_pdpd[22];

    auto g_x_0_0_0_x_xy_x_zz = buffer_1000_pdpd[23];

    auto g_x_0_0_0_x_xy_y_xx = buffer_1000_pdpd[24];

    auto g_x_0_0_0_x_xy_y_xy = buffer_1000_pdpd[25];

    auto g_x_0_0_0_x_xy_y_xz = buffer_1000_pdpd[26];

    auto g_x_0_0_0_x_xy_y_yy = buffer_1000_pdpd[27];

    auto g_x_0_0_0_x_xy_y_yz = buffer_1000_pdpd[28];

    auto g_x_0_0_0_x_xy_y_zz = buffer_1000_pdpd[29];

    auto g_x_0_0_0_x_xy_z_xx = buffer_1000_pdpd[30];

    auto g_x_0_0_0_x_xy_z_xy = buffer_1000_pdpd[31];

    auto g_x_0_0_0_x_xy_z_xz = buffer_1000_pdpd[32];

    auto g_x_0_0_0_x_xy_z_yy = buffer_1000_pdpd[33];

    auto g_x_0_0_0_x_xy_z_yz = buffer_1000_pdpd[34];

    auto g_x_0_0_0_x_xy_z_zz = buffer_1000_pdpd[35];

    auto g_x_0_0_0_x_xz_x_xx = buffer_1000_pdpd[36];

    auto g_x_0_0_0_x_xz_x_xy = buffer_1000_pdpd[37];

    auto g_x_0_0_0_x_xz_x_xz = buffer_1000_pdpd[38];

    auto g_x_0_0_0_x_xz_x_yy = buffer_1000_pdpd[39];

    auto g_x_0_0_0_x_xz_x_yz = buffer_1000_pdpd[40];

    auto g_x_0_0_0_x_xz_x_zz = buffer_1000_pdpd[41];

    auto g_x_0_0_0_x_xz_y_xx = buffer_1000_pdpd[42];

    auto g_x_0_0_0_x_xz_y_xy = buffer_1000_pdpd[43];

    auto g_x_0_0_0_x_xz_y_xz = buffer_1000_pdpd[44];

    auto g_x_0_0_0_x_xz_y_yy = buffer_1000_pdpd[45];

    auto g_x_0_0_0_x_xz_y_yz = buffer_1000_pdpd[46];

    auto g_x_0_0_0_x_xz_y_zz = buffer_1000_pdpd[47];

    auto g_x_0_0_0_x_xz_z_xx = buffer_1000_pdpd[48];

    auto g_x_0_0_0_x_xz_z_xy = buffer_1000_pdpd[49];

    auto g_x_0_0_0_x_xz_z_xz = buffer_1000_pdpd[50];

    auto g_x_0_0_0_x_xz_z_yy = buffer_1000_pdpd[51];

    auto g_x_0_0_0_x_xz_z_yz = buffer_1000_pdpd[52];

    auto g_x_0_0_0_x_xz_z_zz = buffer_1000_pdpd[53];

    auto g_x_0_0_0_x_yy_x_xx = buffer_1000_pdpd[54];

    auto g_x_0_0_0_x_yy_x_xy = buffer_1000_pdpd[55];

    auto g_x_0_0_0_x_yy_x_xz = buffer_1000_pdpd[56];

    auto g_x_0_0_0_x_yy_x_yy = buffer_1000_pdpd[57];

    auto g_x_0_0_0_x_yy_x_yz = buffer_1000_pdpd[58];

    auto g_x_0_0_0_x_yy_x_zz = buffer_1000_pdpd[59];

    auto g_x_0_0_0_x_yy_y_xx = buffer_1000_pdpd[60];

    auto g_x_0_0_0_x_yy_y_xy = buffer_1000_pdpd[61];

    auto g_x_0_0_0_x_yy_y_xz = buffer_1000_pdpd[62];

    auto g_x_0_0_0_x_yy_y_yy = buffer_1000_pdpd[63];

    auto g_x_0_0_0_x_yy_y_yz = buffer_1000_pdpd[64];

    auto g_x_0_0_0_x_yy_y_zz = buffer_1000_pdpd[65];

    auto g_x_0_0_0_x_yy_z_xx = buffer_1000_pdpd[66];

    auto g_x_0_0_0_x_yy_z_xy = buffer_1000_pdpd[67];

    auto g_x_0_0_0_x_yy_z_xz = buffer_1000_pdpd[68];

    auto g_x_0_0_0_x_yy_z_yy = buffer_1000_pdpd[69];

    auto g_x_0_0_0_x_yy_z_yz = buffer_1000_pdpd[70];

    auto g_x_0_0_0_x_yy_z_zz = buffer_1000_pdpd[71];

    auto g_x_0_0_0_x_yz_x_xx = buffer_1000_pdpd[72];

    auto g_x_0_0_0_x_yz_x_xy = buffer_1000_pdpd[73];

    auto g_x_0_0_0_x_yz_x_xz = buffer_1000_pdpd[74];

    auto g_x_0_0_0_x_yz_x_yy = buffer_1000_pdpd[75];

    auto g_x_0_0_0_x_yz_x_yz = buffer_1000_pdpd[76];

    auto g_x_0_0_0_x_yz_x_zz = buffer_1000_pdpd[77];

    auto g_x_0_0_0_x_yz_y_xx = buffer_1000_pdpd[78];

    auto g_x_0_0_0_x_yz_y_xy = buffer_1000_pdpd[79];

    auto g_x_0_0_0_x_yz_y_xz = buffer_1000_pdpd[80];

    auto g_x_0_0_0_x_yz_y_yy = buffer_1000_pdpd[81];

    auto g_x_0_0_0_x_yz_y_yz = buffer_1000_pdpd[82];

    auto g_x_0_0_0_x_yz_y_zz = buffer_1000_pdpd[83];

    auto g_x_0_0_0_x_yz_z_xx = buffer_1000_pdpd[84];

    auto g_x_0_0_0_x_yz_z_xy = buffer_1000_pdpd[85];

    auto g_x_0_0_0_x_yz_z_xz = buffer_1000_pdpd[86];

    auto g_x_0_0_0_x_yz_z_yy = buffer_1000_pdpd[87];

    auto g_x_0_0_0_x_yz_z_yz = buffer_1000_pdpd[88];

    auto g_x_0_0_0_x_yz_z_zz = buffer_1000_pdpd[89];

    auto g_x_0_0_0_x_zz_x_xx = buffer_1000_pdpd[90];

    auto g_x_0_0_0_x_zz_x_xy = buffer_1000_pdpd[91];

    auto g_x_0_0_0_x_zz_x_xz = buffer_1000_pdpd[92];

    auto g_x_0_0_0_x_zz_x_yy = buffer_1000_pdpd[93];

    auto g_x_0_0_0_x_zz_x_yz = buffer_1000_pdpd[94];

    auto g_x_0_0_0_x_zz_x_zz = buffer_1000_pdpd[95];

    auto g_x_0_0_0_x_zz_y_xx = buffer_1000_pdpd[96];

    auto g_x_0_0_0_x_zz_y_xy = buffer_1000_pdpd[97];

    auto g_x_0_0_0_x_zz_y_xz = buffer_1000_pdpd[98];

    auto g_x_0_0_0_x_zz_y_yy = buffer_1000_pdpd[99];

    auto g_x_0_0_0_x_zz_y_yz = buffer_1000_pdpd[100];

    auto g_x_0_0_0_x_zz_y_zz = buffer_1000_pdpd[101];

    auto g_x_0_0_0_x_zz_z_xx = buffer_1000_pdpd[102];

    auto g_x_0_0_0_x_zz_z_xy = buffer_1000_pdpd[103];

    auto g_x_0_0_0_x_zz_z_xz = buffer_1000_pdpd[104];

    auto g_x_0_0_0_x_zz_z_yy = buffer_1000_pdpd[105];

    auto g_x_0_0_0_x_zz_z_yz = buffer_1000_pdpd[106];

    auto g_x_0_0_0_x_zz_z_zz = buffer_1000_pdpd[107];

    auto g_x_0_0_0_y_xx_x_xx = buffer_1000_pdpd[108];

    auto g_x_0_0_0_y_xx_x_xy = buffer_1000_pdpd[109];

    auto g_x_0_0_0_y_xx_x_xz = buffer_1000_pdpd[110];

    auto g_x_0_0_0_y_xx_x_yy = buffer_1000_pdpd[111];

    auto g_x_0_0_0_y_xx_x_yz = buffer_1000_pdpd[112];

    auto g_x_0_0_0_y_xx_x_zz = buffer_1000_pdpd[113];

    auto g_x_0_0_0_y_xx_y_xx = buffer_1000_pdpd[114];

    auto g_x_0_0_0_y_xx_y_xy = buffer_1000_pdpd[115];

    auto g_x_0_0_0_y_xx_y_xz = buffer_1000_pdpd[116];

    auto g_x_0_0_0_y_xx_y_yy = buffer_1000_pdpd[117];

    auto g_x_0_0_0_y_xx_y_yz = buffer_1000_pdpd[118];

    auto g_x_0_0_0_y_xx_y_zz = buffer_1000_pdpd[119];

    auto g_x_0_0_0_y_xx_z_xx = buffer_1000_pdpd[120];

    auto g_x_0_0_0_y_xx_z_xy = buffer_1000_pdpd[121];

    auto g_x_0_0_0_y_xx_z_xz = buffer_1000_pdpd[122];

    auto g_x_0_0_0_y_xx_z_yy = buffer_1000_pdpd[123];

    auto g_x_0_0_0_y_xx_z_yz = buffer_1000_pdpd[124];

    auto g_x_0_0_0_y_xx_z_zz = buffer_1000_pdpd[125];

    auto g_x_0_0_0_y_xy_x_xx = buffer_1000_pdpd[126];

    auto g_x_0_0_0_y_xy_x_xy = buffer_1000_pdpd[127];

    auto g_x_0_0_0_y_xy_x_xz = buffer_1000_pdpd[128];

    auto g_x_0_0_0_y_xy_x_yy = buffer_1000_pdpd[129];

    auto g_x_0_0_0_y_xy_x_yz = buffer_1000_pdpd[130];

    auto g_x_0_0_0_y_xy_x_zz = buffer_1000_pdpd[131];

    auto g_x_0_0_0_y_xy_y_xx = buffer_1000_pdpd[132];

    auto g_x_0_0_0_y_xy_y_xy = buffer_1000_pdpd[133];

    auto g_x_0_0_0_y_xy_y_xz = buffer_1000_pdpd[134];

    auto g_x_0_0_0_y_xy_y_yy = buffer_1000_pdpd[135];

    auto g_x_0_0_0_y_xy_y_yz = buffer_1000_pdpd[136];

    auto g_x_0_0_0_y_xy_y_zz = buffer_1000_pdpd[137];

    auto g_x_0_0_0_y_xy_z_xx = buffer_1000_pdpd[138];

    auto g_x_0_0_0_y_xy_z_xy = buffer_1000_pdpd[139];

    auto g_x_0_0_0_y_xy_z_xz = buffer_1000_pdpd[140];

    auto g_x_0_0_0_y_xy_z_yy = buffer_1000_pdpd[141];

    auto g_x_0_0_0_y_xy_z_yz = buffer_1000_pdpd[142];

    auto g_x_0_0_0_y_xy_z_zz = buffer_1000_pdpd[143];

    auto g_x_0_0_0_y_xz_x_xx = buffer_1000_pdpd[144];

    auto g_x_0_0_0_y_xz_x_xy = buffer_1000_pdpd[145];

    auto g_x_0_0_0_y_xz_x_xz = buffer_1000_pdpd[146];

    auto g_x_0_0_0_y_xz_x_yy = buffer_1000_pdpd[147];

    auto g_x_0_0_0_y_xz_x_yz = buffer_1000_pdpd[148];

    auto g_x_0_0_0_y_xz_x_zz = buffer_1000_pdpd[149];

    auto g_x_0_0_0_y_xz_y_xx = buffer_1000_pdpd[150];

    auto g_x_0_0_0_y_xz_y_xy = buffer_1000_pdpd[151];

    auto g_x_0_0_0_y_xz_y_xz = buffer_1000_pdpd[152];

    auto g_x_0_0_0_y_xz_y_yy = buffer_1000_pdpd[153];

    auto g_x_0_0_0_y_xz_y_yz = buffer_1000_pdpd[154];

    auto g_x_0_0_0_y_xz_y_zz = buffer_1000_pdpd[155];

    auto g_x_0_0_0_y_xz_z_xx = buffer_1000_pdpd[156];

    auto g_x_0_0_0_y_xz_z_xy = buffer_1000_pdpd[157];

    auto g_x_0_0_0_y_xz_z_xz = buffer_1000_pdpd[158];

    auto g_x_0_0_0_y_xz_z_yy = buffer_1000_pdpd[159];

    auto g_x_0_0_0_y_xz_z_yz = buffer_1000_pdpd[160];

    auto g_x_0_0_0_y_xz_z_zz = buffer_1000_pdpd[161];

    auto g_x_0_0_0_y_yy_x_xx = buffer_1000_pdpd[162];

    auto g_x_0_0_0_y_yy_x_xy = buffer_1000_pdpd[163];

    auto g_x_0_0_0_y_yy_x_xz = buffer_1000_pdpd[164];

    auto g_x_0_0_0_y_yy_x_yy = buffer_1000_pdpd[165];

    auto g_x_0_0_0_y_yy_x_yz = buffer_1000_pdpd[166];

    auto g_x_0_0_0_y_yy_x_zz = buffer_1000_pdpd[167];

    auto g_x_0_0_0_y_yy_y_xx = buffer_1000_pdpd[168];

    auto g_x_0_0_0_y_yy_y_xy = buffer_1000_pdpd[169];

    auto g_x_0_0_0_y_yy_y_xz = buffer_1000_pdpd[170];

    auto g_x_0_0_0_y_yy_y_yy = buffer_1000_pdpd[171];

    auto g_x_0_0_0_y_yy_y_yz = buffer_1000_pdpd[172];

    auto g_x_0_0_0_y_yy_y_zz = buffer_1000_pdpd[173];

    auto g_x_0_0_0_y_yy_z_xx = buffer_1000_pdpd[174];

    auto g_x_0_0_0_y_yy_z_xy = buffer_1000_pdpd[175];

    auto g_x_0_0_0_y_yy_z_xz = buffer_1000_pdpd[176];

    auto g_x_0_0_0_y_yy_z_yy = buffer_1000_pdpd[177];

    auto g_x_0_0_0_y_yy_z_yz = buffer_1000_pdpd[178];

    auto g_x_0_0_0_y_yy_z_zz = buffer_1000_pdpd[179];

    auto g_x_0_0_0_y_yz_x_xx = buffer_1000_pdpd[180];

    auto g_x_0_0_0_y_yz_x_xy = buffer_1000_pdpd[181];

    auto g_x_0_0_0_y_yz_x_xz = buffer_1000_pdpd[182];

    auto g_x_0_0_0_y_yz_x_yy = buffer_1000_pdpd[183];

    auto g_x_0_0_0_y_yz_x_yz = buffer_1000_pdpd[184];

    auto g_x_0_0_0_y_yz_x_zz = buffer_1000_pdpd[185];

    auto g_x_0_0_0_y_yz_y_xx = buffer_1000_pdpd[186];

    auto g_x_0_0_0_y_yz_y_xy = buffer_1000_pdpd[187];

    auto g_x_0_0_0_y_yz_y_xz = buffer_1000_pdpd[188];

    auto g_x_0_0_0_y_yz_y_yy = buffer_1000_pdpd[189];

    auto g_x_0_0_0_y_yz_y_yz = buffer_1000_pdpd[190];

    auto g_x_0_0_0_y_yz_y_zz = buffer_1000_pdpd[191];

    auto g_x_0_0_0_y_yz_z_xx = buffer_1000_pdpd[192];

    auto g_x_0_0_0_y_yz_z_xy = buffer_1000_pdpd[193];

    auto g_x_0_0_0_y_yz_z_xz = buffer_1000_pdpd[194];

    auto g_x_0_0_0_y_yz_z_yy = buffer_1000_pdpd[195];

    auto g_x_0_0_0_y_yz_z_yz = buffer_1000_pdpd[196];

    auto g_x_0_0_0_y_yz_z_zz = buffer_1000_pdpd[197];

    auto g_x_0_0_0_y_zz_x_xx = buffer_1000_pdpd[198];

    auto g_x_0_0_0_y_zz_x_xy = buffer_1000_pdpd[199];

    auto g_x_0_0_0_y_zz_x_xz = buffer_1000_pdpd[200];

    auto g_x_0_0_0_y_zz_x_yy = buffer_1000_pdpd[201];

    auto g_x_0_0_0_y_zz_x_yz = buffer_1000_pdpd[202];

    auto g_x_0_0_0_y_zz_x_zz = buffer_1000_pdpd[203];

    auto g_x_0_0_0_y_zz_y_xx = buffer_1000_pdpd[204];

    auto g_x_0_0_0_y_zz_y_xy = buffer_1000_pdpd[205];

    auto g_x_0_0_0_y_zz_y_xz = buffer_1000_pdpd[206];

    auto g_x_0_0_0_y_zz_y_yy = buffer_1000_pdpd[207];

    auto g_x_0_0_0_y_zz_y_yz = buffer_1000_pdpd[208];

    auto g_x_0_0_0_y_zz_y_zz = buffer_1000_pdpd[209];

    auto g_x_0_0_0_y_zz_z_xx = buffer_1000_pdpd[210];

    auto g_x_0_0_0_y_zz_z_xy = buffer_1000_pdpd[211];

    auto g_x_0_0_0_y_zz_z_xz = buffer_1000_pdpd[212];

    auto g_x_0_0_0_y_zz_z_yy = buffer_1000_pdpd[213];

    auto g_x_0_0_0_y_zz_z_yz = buffer_1000_pdpd[214];

    auto g_x_0_0_0_y_zz_z_zz = buffer_1000_pdpd[215];

    auto g_x_0_0_0_z_xx_x_xx = buffer_1000_pdpd[216];

    auto g_x_0_0_0_z_xx_x_xy = buffer_1000_pdpd[217];

    auto g_x_0_0_0_z_xx_x_xz = buffer_1000_pdpd[218];

    auto g_x_0_0_0_z_xx_x_yy = buffer_1000_pdpd[219];

    auto g_x_0_0_0_z_xx_x_yz = buffer_1000_pdpd[220];

    auto g_x_0_0_0_z_xx_x_zz = buffer_1000_pdpd[221];

    auto g_x_0_0_0_z_xx_y_xx = buffer_1000_pdpd[222];

    auto g_x_0_0_0_z_xx_y_xy = buffer_1000_pdpd[223];

    auto g_x_0_0_0_z_xx_y_xz = buffer_1000_pdpd[224];

    auto g_x_0_0_0_z_xx_y_yy = buffer_1000_pdpd[225];

    auto g_x_0_0_0_z_xx_y_yz = buffer_1000_pdpd[226];

    auto g_x_0_0_0_z_xx_y_zz = buffer_1000_pdpd[227];

    auto g_x_0_0_0_z_xx_z_xx = buffer_1000_pdpd[228];

    auto g_x_0_0_0_z_xx_z_xy = buffer_1000_pdpd[229];

    auto g_x_0_0_0_z_xx_z_xz = buffer_1000_pdpd[230];

    auto g_x_0_0_0_z_xx_z_yy = buffer_1000_pdpd[231];

    auto g_x_0_0_0_z_xx_z_yz = buffer_1000_pdpd[232];

    auto g_x_0_0_0_z_xx_z_zz = buffer_1000_pdpd[233];

    auto g_x_0_0_0_z_xy_x_xx = buffer_1000_pdpd[234];

    auto g_x_0_0_0_z_xy_x_xy = buffer_1000_pdpd[235];

    auto g_x_0_0_0_z_xy_x_xz = buffer_1000_pdpd[236];

    auto g_x_0_0_0_z_xy_x_yy = buffer_1000_pdpd[237];

    auto g_x_0_0_0_z_xy_x_yz = buffer_1000_pdpd[238];

    auto g_x_0_0_0_z_xy_x_zz = buffer_1000_pdpd[239];

    auto g_x_0_0_0_z_xy_y_xx = buffer_1000_pdpd[240];

    auto g_x_0_0_0_z_xy_y_xy = buffer_1000_pdpd[241];

    auto g_x_0_0_0_z_xy_y_xz = buffer_1000_pdpd[242];

    auto g_x_0_0_0_z_xy_y_yy = buffer_1000_pdpd[243];

    auto g_x_0_0_0_z_xy_y_yz = buffer_1000_pdpd[244];

    auto g_x_0_0_0_z_xy_y_zz = buffer_1000_pdpd[245];

    auto g_x_0_0_0_z_xy_z_xx = buffer_1000_pdpd[246];

    auto g_x_0_0_0_z_xy_z_xy = buffer_1000_pdpd[247];

    auto g_x_0_0_0_z_xy_z_xz = buffer_1000_pdpd[248];

    auto g_x_0_0_0_z_xy_z_yy = buffer_1000_pdpd[249];

    auto g_x_0_0_0_z_xy_z_yz = buffer_1000_pdpd[250];

    auto g_x_0_0_0_z_xy_z_zz = buffer_1000_pdpd[251];

    auto g_x_0_0_0_z_xz_x_xx = buffer_1000_pdpd[252];

    auto g_x_0_0_0_z_xz_x_xy = buffer_1000_pdpd[253];

    auto g_x_0_0_0_z_xz_x_xz = buffer_1000_pdpd[254];

    auto g_x_0_0_0_z_xz_x_yy = buffer_1000_pdpd[255];

    auto g_x_0_0_0_z_xz_x_yz = buffer_1000_pdpd[256];

    auto g_x_0_0_0_z_xz_x_zz = buffer_1000_pdpd[257];

    auto g_x_0_0_0_z_xz_y_xx = buffer_1000_pdpd[258];

    auto g_x_0_0_0_z_xz_y_xy = buffer_1000_pdpd[259];

    auto g_x_0_0_0_z_xz_y_xz = buffer_1000_pdpd[260];

    auto g_x_0_0_0_z_xz_y_yy = buffer_1000_pdpd[261];

    auto g_x_0_0_0_z_xz_y_yz = buffer_1000_pdpd[262];

    auto g_x_0_0_0_z_xz_y_zz = buffer_1000_pdpd[263];

    auto g_x_0_0_0_z_xz_z_xx = buffer_1000_pdpd[264];

    auto g_x_0_0_0_z_xz_z_xy = buffer_1000_pdpd[265];

    auto g_x_0_0_0_z_xz_z_xz = buffer_1000_pdpd[266];

    auto g_x_0_0_0_z_xz_z_yy = buffer_1000_pdpd[267];

    auto g_x_0_0_0_z_xz_z_yz = buffer_1000_pdpd[268];

    auto g_x_0_0_0_z_xz_z_zz = buffer_1000_pdpd[269];

    auto g_x_0_0_0_z_yy_x_xx = buffer_1000_pdpd[270];

    auto g_x_0_0_0_z_yy_x_xy = buffer_1000_pdpd[271];

    auto g_x_0_0_0_z_yy_x_xz = buffer_1000_pdpd[272];

    auto g_x_0_0_0_z_yy_x_yy = buffer_1000_pdpd[273];

    auto g_x_0_0_0_z_yy_x_yz = buffer_1000_pdpd[274];

    auto g_x_0_0_0_z_yy_x_zz = buffer_1000_pdpd[275];

    auto g_x_0_0_0_z_yy_y_xx = buffer_1000_pdpd[276];

    auto g_x_0_0_0_z_yy_y_xy = buffer_1000_pdpd[277];

    auto g_x_0_0_0_z_yy_y_xz = buffer_1000_pdpd[278];

    auto g_x_0_0_0_z_yy_y_yy = buffer_1000_pdpd[279];

    auto g_x_0_0_0_z_yy_y_yz = buffer_1000_pdpd[280];

    auto g_x_0_0_0_z_yy_y_zz = buffer_1000_pdpd[281];

    auto g_x_0_0_0_z_yy_z_xx = buffer_1000_pdpd[282];

    auto g_x_0_0_0_z_yy_z_xy = buffer_1000_pdpd[283];

    auto g_x_0_0_0_z_yy_z_xz = buffer_1000_pdpd[284];

    auto g_x_0_0_0_z_yy_z_yy = buffer_1000_pdpd[285];

    auto g_x_0_0_0_z_yy_z_yz = buffer_1000_pdpd[286];

    auto g_x_0_0_0_z_yy_z_zz = buffer_1000_pdpd[287];

    auto g_x_0_0_0_z_yz_x_xx = buffer_1000_pdpd[288];

    auto g_x_0_0_0_z_yz_x_xy = buffer_1000_pdpd[289];

    auto g_x_0_0_0_z_yz_x_xz = buffer_1000_pdpd[290];

    auto g_x_0_0_0_z_yz_x_yy = buffer_1000_pdpd[291];

    auto g_x_0_0_0_z_yz_x_yz = buffer_1000_pdpd[292];

    auto g_x_0_0_0_z_yz_x_zz = buffer_1000_pdpd[293];

    auto g_x_0_0_0_z_yz_y_xx = buffer_1000_pdpd[294];

    auto g_x_0_0_0_z_yz_y_xy = buffer_1000_pdpd[295];

    auto g_x_0_0_0_z_yz_y_xz = buffer_1000_pdpd[296];

    auto g_x_0_0_0_z_yz_y_yy = buffer_1000_pdpd[297];

    auto g_x_0_0_0_z_yz_y_yz = buffer_1000_pdpd[298];

    auto g_x_0_0_0_z_yz_y_zz = buffer_1000_pdpd[299];

    auto g_x_0_0_0_z_yz_z_xx = buffer_1000_pdpd[300];

    auto g_x_0_0_0_z_yz_z_xy = buffer_1000_pdpd[301];

    auto g_x_0_0_0_z_yz_z_xz = buffer_1000_pdpd[302];

    auto g_x_0_0_0_z_yz_z_yy = buffer_1000_pdpd[303];

    auto g_x_0_0_0_z_yz_z_yz = buffer_1000_pdpd[304];

    auto g_x_0_0_0_z_yz_z_zz = buffer_1000_pdpd[305];

    auto g_x_0_0_0_z_zz_x_xx = buffer_1000_pdpd[306];

    auto g_x_0_0_0_z_zz_x_xy = buffer_1000_pdpd[307];

    auto g_x_0_0_0_z_zz_x_xz = buffer_1000_pdpd[308];

    auto g_x_0_0_0_z_zz_x_yy = buffer_1000_pdpd[309];

    auto g_x_0_0_0_z_zz_x_yz = buffer_1000_pdpd[310];

    auto g_x_0_0_0_z_zz_x_zz = buffer_1000_pdpd[311];

    auto g_x_0_0_0_z_zz_y_xx = buffer_1000_pdpd[312];

    auto g_x_0_0_0_z_zz_y_xy = buffer_1000_pdpd[313];

    auto g_x_0_0_0_z_zz_y_xz = buffer_1000_pdpd[314];

    auto g_x_0_0_0_z_zz_y_yy = buffer_1000_pdpd[315];

    auto g_x_0_0_0_z_zz_y_yz = buffer_1000_pdpd[316];

    auto g_x_0_0_0_z_zz_y_zz = buffer_1000_pdpd[317];

    auto g_x_0_0_0_z_zz_z_xx = buffer_1000_pdpd[318];

    auto g_x_0_0_0_z_zz_z_xy = buffer_1000_pdpd[319];

    auto g_x_0_0_0_z_zz_z_xz = buffer_1000_pdpd[320];

    auto g_x_0_0_0_z_zz_z_yy = buffer_1000_pdpd[321];

    auto g_x_0_0_0_z_zz_z_yz = buffer_1000_pdpd[322];

    auto g_x_0_0_0_z_zz_z_zz = buffer_1000_pdpd[323];

    auto g_y_0_0_0_x_xx_x_xx = buffer_1000_pdpd[324];

    auto g_y_0_0_0_x_xx_x_xy = buffer_1000_pdpd[325];

    auto g_y_0_0_0_x_xx_x_xz = buffer_1000_pdpd[326];

    auto g_y_0_0_0_x_xx_x_yy = buffer_1000_pdpd[327];

    auto g_y_0_0_0_x_xx_x_yz = buffer_1000_pdpd[328];

    auto g_y_0_0_0_x_xx_x_zz = buffer_1000_pdpd[329];

    auto g_y_0_0_0_x_xx_y_xx = buffer_1000_pdpd[330];

    auto g_y_0_0_0_x_xx_y_xy = buffer_1000_pdpd[331];

    auto g_y_0_0_0_x_xx_y_xz = buffer_1000_pdpd[332];

    auto g_y_0_0_0_x_xx_y_yy = buffer_1000_pdpd[333];

    auto g_y_0_0_0_x_xx_y_yz = buffer_1000_pdpd[334];

    auto g_y_0_0_0_x_xx_y_zz = buffer_1000_pdpd[335];

    auto g_y_0_0_0_x_xx_z_xx = buffer_1000_pdpd[336];

    auto g_y_0_0_0_x_xx_z_xy = buffer_1000_pdpd[337];

    auto g_y_0_0_0_x_xx_z_xz = buffer_1000_pdpd[338];

    auto g_y_0_0_0_x_xx_z_yy = buffer_1000_pdpd[339];

    auto g_y_0_0_0_x_xx_z_yz = buffer_1000_pdpd[340];

    auto g_y_0_0_0_x_xx_z_zz = buffer_1000_pdpd[341];

    auto g_y_0_0_0_x_xy_x_xx = buffer_1000_pdpd[342];

    auto g_y_0_0_0_x_xy_x_xy = buffer_1000_pdpd[343];

    auto g_y_0_0_0_x_xy_x_xz = buffer_1000_pdpd[344];

    auto g_y_0_0_0_x_xy_x_yy = buffer_1000_pdpd[345];

    auto g_y_0_0_0_x_xy_x_yz = buffer_1000_pdpd[346];

    auto g_y_0_0_0_x_xy_x_zz = buffer_1000_pdpd[347];

    auto g_y_0_0_0_x_xy_y_xx = buffer_1000_pdpd[348];

    auto g_y_0_0_0_x_xy_y_xy = buffer_1000_pdpd[349];

    auto g_y_0_0_0_x_xy_y_xz = buffer_1000_pdpd[350];

    auto g_y_0_0_0_x_xy_y_yy = buffer_1000_pdpd[351];

    auto g_y_0_0_0_x_xy_y_yz = buffer_1000_pdpd[352];

    auto g_y_0_0_0_x_xy_y_zz = buffer_1000_pdpd[353];

    auto g_y_0_0_0_x_xy_z_xx = buffer_1000_pdpd[354];

    auto g_y_0_0_0_x_xy_z_xy = buffer_1000_pdpd[355];

    auto g_y_0_0_0_x_xy_z_xz = buffer_1000_pdpd[356];

    auto g_y_0_0_0_x_xy_z_yy = buffer_1000_pdpd[357];

    auto g_y_0_0_0_x_xy_z_yz = buffer_1000_pdpd[358];

    auto g_y_0_0_0_x_xy_z_zz = buffer_1000_pdpd[359];

    auto g_y_0_0_0_x_xz_x_xx = buffer_1000_pdpd[360];

    auto g_y_0_0_0_x_xz_x_xy = buffer_1000_pdpd[361];

    auto g_y_0_0_0_x_xz_x_xz = buffer_1000_pdpd[362];

    auto g_y_0_0_0_x_xz_x_yy = buffer_1000_pdpd[363];

    auto g_y_0_0_0_x_xz_x_yz = buffer_1000_pdpd[364];

    auto g_y_0_0_0_x_xz_x_zz = buffer_1000_pdpd[365];

    auto g_y_0_0_0_x_xz_y_xx = buffer_1000_pdpd[366];

    auto g_y_0_0_0_x_xz_y_xy = buffer_1000_pdpd[367];

    auto g_y_0_0_0_x_xz_y_xz = buffer_1000_pdpd[368];

    auto g_y_0_0_0_x_xz_y_yy = buffer_1000_pdpd[369];

    auto g_y_0_0_0_x_xz_y_yz = buffer_1000_pdpd[370];

    auto g_y_0_0_0_x_xz_y_zz = buffer_1000_pdpd[371];

    auto g_y_0_0_0_x_xz_z_xx = buffer_1000_pdpd[372];

    auto g_y_0_0_0_x_xz_z_xy = buffer_1000_pdpd[373];

    auto g_y_0_0_0_x_xz_z_xz = buffer_1000_pdpd[374];

    auto g_y_0_0_0_x_xz_z_yy = buffer_1000_pdpd[375];

    auto g_y_0_0_0_x_xz_z_yz = buffer_1000_pdpd[376];

    auto g_y_0_0_0_x_xz_z_zz = buffer_1000_pdpd[377];

    auto g_y_0_0_0_x_yy_x_xx = buffer_1000_pdpd[378];

    auto g_y_0_0_0_x_yy_x_xy = buffer_1000_pdpd[379];

    auto g_y_0_0_0_x_yy_x_xz = buffer_1000_pdpd[380];

    auto g_y_0_0_0_x_yy_x_yy = buffer_1000_pdpd[381];

    auto g_y_0_0_0_x_yy_x_yz = buffer_1000_pdpd[382];

    auto g_y_0_0_0_x_yy_x_zz = buffer_1000_pdpd[383];

    auto g_y_0_0_0_x_yy_y_xx = buffer_1000_pdpd[384];

    auto g_y_0_0_0_x_yy_y_xy = buffer_1000_pdpd[385];

    auto g_y_0_0_0_x_yy_y_xz = buffer_1000_pdpd[386];

    auto g_y_0_0_0_x_yy_y_yy = buffer_1000_pdpd[387];

    auto g_y_0_0_0_x_yy_y_yz = buffer_1000_pdpd[388];

    auto g_y_0_0_0_x_yy_y_zz = buffer_1000_pdpd[389];

    auto g_y_0_0_0_x_yy_z_xx = buffer_1000_pdpd[390];

    auto g_y_0_0_0_x_yy_z_xy = buffer_1000_pdpd[391];

    auto g_y_0_0_0_x_yy_z_xz = buffer_1000_pdpd[392];

    auto g_y_0_0_0_x_yy_z_yy = buffer_1000_pdpd[393];

    auto g_y_0_0_0_x_yy_z_yz = buffer_1000_pdpd[394];

    auto g_y_0_0_0_x_yy_z_zz = buffer_1000_pdpd[395];

    auto g_y_0_0_0_x_yz_x_xx = buffer_1000_pdpd[396];

    auto g_y_0_0_0_x_yz_x_xy = buffer_1000_pdpd[397];

    auto g_y_0_0_0_x_yz_x_xz = buffer_1000_pdpd[398];

    auto g_y_0_0_0_x_yz_x_yy = buffer_1000_pdpd[399];

    auto g_y_0_0_0_x_yz_x_yz = buffer_1000_pdpd[400];

    auto g_y_0_0_0_x_yz_x_zz = buffer_1000_pdpd[401];

    auto g_y_0_0_0_x_yz_y_xx = buffer_1000_pdpd[402];

    auto g_y_0_0_0_x_yz_y_xy = buffer_1000_pdpd[403];

    auto g_y_0_0_0_x_yz_y_xz = buffer_1000_pdpd[404];

    auto g_y_0_0_0_x_yz_y_yy = buffer_1000_pdpd[405];

    auto g_y_0_0_0_x_yz_y_yz = buffer_1000_pdpd[406];

    auto g_y_0_0_0_x_yz_y_zz = buffer_1000_pdpd[407];

    auto g_y_0_0_0_x_yz_z_xx = buffer_1000_pdpd[408];

    auto g_y_0_0_0_x_yz_z_xy = buffer_1000_pdpd[409];

    auto g_y_0_0_0_x_yz_z_xz = buffer_1000_pdpd[410];

    auto g_y_0_0_0_x_yz_z_yy = buffer_1000_pdpd[411];

    auto g_y_0_0_0_x_yz_z_yz = buffer_1000_pdpd[412];

    auto g_y_0_0_0_x_yz_z_zz = buffer_1000_pdpd[413];

    auto g_y_0_0_0_x_zz_x_xx = buffer_1000_pdpd[414];

    auto g_y_0_0_0_x_zz_x_xy = buffer_1000_pdpd[415];

    auto g_y_0_0_0_x_zz_x_xz = buffer_1000_pdpd[416];

    auto g_y_0_0_0_x_zz_x_yy = buffer_1000_pdpd[417];

    auto g_y_0_0_0_x_zz_x_yz = buffer_1000_pdpd[418];

    auto g_y_0_0_0_x_zz_x_zz = buffer_1000_pdpd[419];

    auto g_y_0_0_0_x_zz_y_xx = buffer_1000_pdpd[420];

    auto g_y_0_0_0_x_zz_y_xy = buffer_1000_pdpd[421];

    auto g_y_0_0_0_x_zz_y_xz = buffer_1000_pdpd[422];

    auto g_y_0_0_0_x_zz_y_yy = buffer_1000_pdpd[423];

    auto g_y_0_0_0_x_zz_y_yz = buffer_1000_pdpd[424];

    auto g_y_0_0_0_x_zz_y_zz = buffer_1000_pdpd[425];

    auto g_y_0_0_0_x_zz_z_xx = buffer_1000_pdpd[426];

    auto g_y_0_0_0_x_zz_z_xy = buffer_1000_pdpd[427];

    auto g_y_0_0_0_x_zz_z_xz = buffer_1000_pdpd[428];

    auto g_y_0_0_0_x_zz_z_yy = buffer_1000_pdpd[429];

    auto g_y_0_0_0_x_zz_z_yz = buffer_1000_pdpd[430];

    auto g_y_0_0_0_x_zz_z_zz = buffer_1000_pdpd[431];

    auto g_y_0_0_0_y_xx_x_xx = buffer_1000_pdpd[432];

    auto g_y_0_0_0_y_xx_x_xy = buffer_1000_pdpd[433];

    auto g_y_0_0_0_y_xx_x_xz = buffer_1000_pdpd[434];

    auto g_y_0_0_0_y_xx_x_yy = buffer_1000_pdpd[435];

    auto g_y_0_0_0_y_xx_x_yz = buffer_1000_pdpd[436];

    auto g_y_0_0_0_y_xx_x_zz = buffer_1000_pdpd[437];

    auto g_y_0_0_0_y_xx_y_xx = buffer_1000_pdpd[438];

    auto g_y_0_0_0_y_xx_y_xy = buffer_1000_pdpd[439];

    auto g_y_0_0_0_y_xx_y_xz = buffer_1000_pdpd[440];

    auto g_y_0_0_0_y_xx_y_yy = buffer_1000_pdpd[441];

    auto g_y_0_0_0_y_xx_y_yz = buffer_1000_pdpd[442];

    auto g_y_0_0_0_y_xx_y_zz = buffer_1000_pdpd[443];

    auto g_y_0_0_0_y_xx_z_xx = buffer_1000_pdpd[444];

    auto g_y_0_0_0_y_xx_z_xy = buffer_1000_pdpd[445];

    auto g_y_0_0_0_y_xx_z_xz = buffer_1000_pdpd[446];

    auto g_y_0_0_0_y_xx_z_yy = buffer_1000_pdpd[447];

    auto g_y_0_0_0_y_xx_z_yz = buffer_1000_pdpd[448];

    auto g_y_0_0_0_y_xx_z_zz = buffer_1000_pdpd[449];

    auto g_y_0_0_0_y_xy_x_xx = buffer_1000_pdpd[450];

    auto g_y_0_0_0_y_xy_x_xy = buffer_1000_pdpd[451];

    auto g_y_0_0_0_y_xy_x_xz = buffer_1000_pdpd[452];

    auto g_y_0_0_0_y_xy_x_yy = buffer_1000_pdpd[453];

    auto g_y_0_0_0_y_xy_x_yz = buffer_1000_pdpd[454];

    auto g_y_0_0_0_y_xy_x_zz = buffer_1000_pdpd[455];

    auto g_y_0_0_0_y_xy_y_xx = buffer_1000_pdpd[456];

    auto g_y_0_0_0_y_xy_y_xy = buffer_1000_pdpd[457];

    auto g_y_0_0_0_y_xy_y_xz = buffer_1000_pdpd[458];

    auto g_y_0_0_0_y_xy_y_yy = buffer_1000_pdpd[459];

    auto g_y_0_0_0_y_xy_y_yz = buffer_1000_pdpd[460];

    auto g_y_0_0_0_y_xy_y_zz = buffer_1000_pdpd[461];

    auto g_y_0_0_0_y_xy_z_xx = buffer_1000_pdpd[462];

    auto g_y_0_0_0_y_xy_z_xy = buffer_1000_pdpd[463];

    auto g_y_0_0_0_y_xy_z_xz = buffer_1000_pdpd[464];

    auto g_y_0_0_0_y_xy_z_yy = buffer_1000_pdpd[465];

    auto g_y_0_0_0_y_xy_z_yz = buffer_1000_pdpd[466];

    auto g_y_0_0_0_y_xy_z_zz = buffer_1000_pdpd[467];

    auto g_y_0_0_0_y_xz_x_xx = buffer_1000_pdpd[468];

    auto g_y_0_0_0_y_xz_x_xy = buffer_1000_pdpd[469];

    auto g_y_0_0_0_y_xz_x_xz = buffer_1000_pdpd[470];

    auto g_y_0_0_0_y_xz_x_yy = buffer_1000_pdpd[471];

    auto g_y_0_0_0_y_xz_x_yz = buffer_1000_pdpd[472];

    auto g_y_0_0_0_y_xz_x_zz = buffer_1000_pdpd[473];

    auto g_y_0_0_0_y_xz_y_xx = buffer_1000_pdpd[474];

    auto g_y_0_0_0_y_xz_y_xy = buffer_1000_pdpd[475];

    auto g_y_0_0_0_y_xz_y_xz = buffer_1000_pdpd[476];

    auto g_y_0_0_0_y_xz_y_yy = buffer_1000_pdpd[477];

    auto g_y_0_0_0_y_xz_y_yz = buffer_1000_pdpd[478];

    auto g_y_0_0_0_y_xz_y_zz = buffer_1000_pdpd[479];

    auto g_y_0_0_0_y_xz_z_xx = buffer_1000_pdpd[480];

    auto g_y_0_0_0_y_xz_z_xy = buffer_1000_pdpd[481];

    auto g_y_0_0_0_y_xz_z_xz = buffer_1000_pdpd[482];

    auto g_y_0_0_0_y_xz_z_yy = buffer_1000_pdpd[483];

    auto g_y_0_0_0_y_xz_z_yz = buffer_1000_pdpd[484];

    auto g_y_0_0_0_y_xz_z_zz = buffer_1000_pdpd[485];

    auto g_y_0_0_0_y_yy_x_xx = buffer_1000_pdpd[486];

    auto g_y_0_0_0_y_yy_x_xy = buffer_1000_pdpd[487];

    auto g_y_0_0_0_y_yy_x_xz = buffer_1000_pdpd[488];

    auto g_y_0_0_0_y_yy_x_yy = buffer_1000_pdpd[489];

    auto g_y_0_0_0_y_yy_x_yz = buffer_1000_pdpd[490];

    auto g_y_0_0_0_y_yy_x_zz = buffer_1000_pdpd[491];

    auto g_y_0_0_0_y_yy_y_xx = buffer_1000_pdpd[492];

    auto g_y_0_0_0_y_yy_y_xy = buffer_1000_pdpd[493];

    auto g_y_0_0_0_y_yy_y_xz = buffer_1000_pdpd[494];

    auto g_y_0_0_0_y_yy_y_yy = buffer_1000_pdpd[495];

    auto g_y_0_0_0_y_yy_y_yz = buffer_1000_pdpd[496];

    auto g_y_0_0_0_y_yy_y_zz = buffer_1000_pdpd[497];

    auto g_y_0_0_0_y_yy_z_xx = buffer_1000_pdpd[498];

    auto g_y_0_0_0_y_yy_z_xy = buffer_1000_pdpd[499];

    auto g_y_0_0_0_y_yy_z_xz = buffer_1000_pdpd[500];

    auto g_y_0_0_0_y_yy_z_yy = buffer_1000_pdpd[501];

    auto g_y_0_0_0_y_yy_z_yz = buffer_1000_pdpd[502];

    auto g_y_0_0_0_y_yy_z_zz = buffer_1000_pdpd[503];

    auto g_y_0_0_0_y_yz_x_xx = buffer_1000_pdpd[504];

    auto g_y_0_0_0_y_yz_x_xy = buffer_1000_pdpd[505];

    auto g_y_0_0_0_y_yz_x_xz = buffer_1000_pdpd[506];

    auto g_y_0_0_0_y_yz_x_yy = buffer_1000_pdpd[507];

    auto g_y_0_0_0_y_yz_x_yz = buffer_1000_pdpd[508];

    auto g_y_0_0_0_y_yz_x_zz = buffer_1000_pdpd[509];

    auto g_y_0_0_0_y_yz_y_xx = buffer_1000_pdpd[510];

    auto g_y_0_0_0_y_yz_y_xy = buffer_1000_pdpd[511];

    auto g_y_0_0_0_y_yz_y_xz = buffer_1000_pdpd[512];

    auto g_y_0_0_0_y_yz_y_yy = buffer_1000_pdpd[513];

    auto g_y_0_0_0_y_yz_y_yz = buffer_1000_pdpd[514];

    auto g_y_0_0_0_y_yz_y_zz = buffer_1000_pdpd[515];

    auto g_y_0_0_0_y_yz_z_xx = buffer_1000_pdpd[516];

    auto g_y_0_0_0_y_yz_z_xy = buffer_1000_pdpd[517];

    auto g_y_0_0_0_y_yz_z_xz = buffer_1000_pdpd[518];

    auto g_y_0_0_0_y_yz_z_yy = buffer_1000_pdpd[519];

    auto g_y_0_0_0_y_yz_z_yz = buffer_1000_pdpd[520];

    auto g_y_0_0_0_y_yz_z_zz = buffer_1000_pdpd[521];

    auto g_y_0_0_0_y_zz_x_xx = buffer_1000_pdpd[522];

    auto g_y_0_0_0_y_zz_x_xy = buffer_1000_pdpd[523];

    auto g_y_0_0_0_y_zz_x_xz = buffer_1000_pdpd[524];

    auto g_y_0_0_0_y_zz_x_yy = buffer_1000_pdpd[525];

    auto g_y_0_0_0_y_zz_x_yz = buffer_1000_pdpd[526];

    auto g_y_0_0_0_y_zz_x_zz = buffer_1000_pdpd[527];

    auto g_y_0_0_0_y_zz_y_xx = buffer_1000_pdpd[528];

    auto g_y_0_0_0_y_zz_y_xy = buffer_1000_pdpd[529];

    auto g_y_0_0_0_y_zz_y_xz = buffer_1000_pdpd[530];

    auto g_y_0_0_0_y_zz_y_yy = buffer_1000_pdpd[531];

    auto g_y_0_0_0_y_zz_y_yz = buffer_1000_pdpd[532];

    auto g_y_0_0_0_y_zz_y_zz = buffer_1000_pdpd[533];

    auto g_y_0_0_0_y_zz_z_xx = buffer_1000_pdpd[534];

    auto g_y_0_0_0_y_zz_z_xy = buffer_1000_pdpd[535];

    auto g_y_0_0_0_y_zz_z_xz = buffer_1000_pdpd[536];

    auto g_y_0_0_0_y_zz_z_yy = buffer_1000_pdpd[537];

    auto g_y_0_0_0_y_zz_z_yz = buffer_1000_pdpd[538];

    auto g_y_0_0_0_y_zz_z_zz = buffer_1000_pdpd[539];

    auto g_y_0_0_0_z_xx_x_xx = buffer_1000_pdpd[540];

    auto g_y_0_0_0_z_xx_x_xy = buffer_1000_pdpd[541];

    auto g_y_0_0_0_z_xx_x_xz = buffer_1000_pdpd[542];

    auto g_y_0_0_0_z_xx_x_yy = buffer_1000_pdpd[543];

    auto g_y_0_0_0_z_xx_x_yz = buffer_1000_pdpd[544];

    auto g_y_0_0_0_z_xx_x_zz = buffer_1000_pdpd[545];

    auto g_y_0_0_0_z_xx_y_xx = buffer_1000_pdpd[546];

    auto g_y_0_0_0_z_xx_y_xy = buffer_1000_pdpd[547];

    auto g_y_0_0_0_z_xx_y_xz = buffer_1000_pdpd[548];

    auto g_y_0_0_0_z_xx_y_yy = buffer_1000_pdpd[549];

    auto g_y_0_0_0_z_xx_y_yz = buffer_1000_pdpd[550];

    auto g_y_0_0_0_z_xx_y_zz = buffer_1000_pdpd[551];

    auto g_y_0_0_0_z_xx_z_xx = buffer_1000_pdpd[552];

    auto g_y_0_0_0_z_xx_z_xy = buffer_1000_pdpd[553];

    auto g_y_0_0_0_z_xx_z_xz = buffer_1000_pdpd[554];

    auto g_y_0_0_0_z_xx_z_yy = buffer_1000_pdpd[555];

    auto g_y_0_0_0_z_xx_z_yz = buffer_1000_pdpd[556];

    auto g_y_0_0_0_z_xx_z_zz = buffer_1000_pdpd[557];

    auto g_y_0_0_0_z_xy_x_xx = buffer_1000_pdpd[558];

    auto g_y_0_0_0_z_xy_x_xy = buffer_1000_pdpd[559];

    auto g_y_0_0_0_z_xy_x_xz = buffer_1000_pdpd[560];

    auto g_y_0_0_0_z_xy_x_yy = buffer_1000_pdpd[561];

    auto g_y_0_0_0_z_xy_x_yz = buffer_1000_pdpd[562];

    auto g_y_0_0_0_z_xy_x_zz = buffer_1000_pdpd[563];

    auto g_y_0_0_0_z_xy_y_xx = buffer_1000_pdpd[564];

    auto g_y_0_0_0_z_xy_y_xy = buffer_1000_pdpd[565];

    auto g_y_0_0_0_z_xy_y_xz = buffer_1000_pdpd[566];

    auto g_y_0_0_0_z_xy_y_yy = buffer_1000_pdpd[567];

    auto g_y_0_0_0_z_xy_y_yz = buffer_1000_pdpd[568];

    auto g_y_0_0_0_z_xy_y_zz = buffer_1000_pdpd[569];

    auto g_y_0_0_0_z_xy_z_xx = buffer_1000_pdpd[570];

    auto g_y_0_0_0_z_xy_z_xy = buffer_1000_pdpd[571];

    auto g_y_0_0_0_z_xy_z_xz = buffer_1000_pdpd[572];

    auto g_y_0_0_0_z_xy_z_yy = buffer_1000_pdpd[573];

    auto g_y_0_0_0_z_xy_z_yz = buffer_1000_pdpd[574];

    auto g_y_0_0_0_z_xy_z_zz = buffer_1000_pdpd[575];

    auto g_y_0_0_0_z_xz_x_xx = buffer_1000_pdpd[576];

    auto g_y_0_0_0_z_xz_x_xy = buffer_1000_pdpd[577];

    auto g_y_0_0_0_z_xz_x_xz = buffer_1000_pdpd[578];

    auto g_y_0_0_0_z_xz_x_yy = buffer_1000_pdpd[579];

    auto g_y_0_0_0_z_xz_x_yz = buffer_1000_pdpd[580];

    auto g_y_0_0_0_z_xz_x_zz = buffer_1000_pdpd[581];

    auto g_y_0_0_0_z_xz_y_xx = buffer_1000_pdpd[582];

    auto g_y_0_0_0_z_xz_y_xy = buffer_1000_pdpd[583];

    auto g_y_0_0_0_z_xz_y_xz = buffer_1000_pdpd[584];

    auto g_y_0_0_0_z_xz_y_yy = buffer_1000_pdpd[585];

    auto g_y_0_0_0_z_xz_y_yz = buffer_1000_pdpd[586];

    auto g_y_0_0_0_z_xz_y_zz = buffer_1000_pdpd[587];

    auto g_y_0_0_0_z_xz_z_xx = buffer_1000_pdpd[588];

    auto g_y_0_0_0_z_xz_z_xy = buffer_1000_pdpd[589];

    auto g_y_0_0_0_z_xz_z_xz = buffer_1000_pdpd[590];

    auto g_y_0_0_0_z_xz_z_yy = buffer_1000_pdpd[591];

    auto g_y_0_0_0_z_xz_z_yz = buffer_1000_pdpd[592];

    auto g_y_0_0_0_z_xz_z_zz = buffer_1000_pdpd[593];

    auto g_y_0_0_0_z_yy_x_xx = buffer_1000_pdpd[594];

    auto g_y_0_0_0_z_yy_x_xy = buffer_1000_pdpd[595];

    auto g_y_0_0_0_z_yy_x_xz = buffer_1000_pdpd[596];

    auto g_y_0_0_0_z_yy_x_yy = buffer_1000_pdpd[597];

    auto g_y_0_0_0_z_yy_x_yz = buffer_1000_pdpd[598];

    auto g_y_0_0_0_z_yy_x_zz = buffer_1000_pdpd[599];

    auto g_y_0_0_0_z_yy_y_xx = buffer_1000_pdpd[600];

    auto g_y_0_0_0_z_yy_y_xy = buffer_1000_pdpd[601];

    auto g_y_0_0_0_z_yy_y_xz = buffer_1000_pdpd[602];

    auto g_y_0_0_0_z_yy_y_yy = buffer_1000_pdpd[603];

    auto g_y_0_0_0_z_yy_y_yz = buffer_1000_pdpd[604];

    auto g_y_0_0_0_z_yy_y_zz = buffer_1000_pdpd[605];

    auto g_y_0_0_0_z_yy_z_xx = buffer_1000_pdpd[606];

    auto g_y_0_0_0_z_yy_z_xy = buffer_1000_pdpd[607];

    auto g_y_0_0_0_z_yy_z_xz = buffer_1000_pdpd[608];

    auto g_y_0_0_0_z_yy_z_yy = buffer_1000_pdpd[609];

    auto g_y_0_0_0_z_yy_z_yz = buffer_1000_pdpd[610];

    auto g_y_0_0_0_z_yy_z_zz = buffer_1000_pdpd[611];

    auto g_y_0_0_0_z_yz_x_xx = buffer_1000_pdpd[612];

    auto g_y_0_0_0_z_yz_x_xy = buffer_1000_pdpd[613];

    auto g_y_0_0_0_z_yz_x_xz = buffer_1000_pdpd[614];

    auto g_y_0_0_0_z_yz_x_yy = buffer_1000_pdpd[615];

    auto g_y_0_0_0_z_yz_x_yz = buffer_1000_pdpd[616];

    auto g_y_0_0_0_z_yz_x_zz = buffer_1000_pdpd[617];

    auto g_y_0_0_0_z_yz_y_xx = buffer_1000_pdpd[618];

    auto g_y_0_0_0_z_yz_y_xy = buffer_1000_pdpd[619];

    auto g_y_0_0_0_z_yz_y_xz = buffer_1000_pdpd[620];

    auto g_y_0_0_0_z_yz_y_yy = buffer_1000_pdpd[621];

    auto g_y_0_0_0_z_yz_y_yz = buffer_1000_pdpd[622];

    auto g_y_0_0_0_z_yz_y_zz = buffer_1000_pdpd[623];

    auto g_y_0_0_0_z_yz_z_xx = buffer_1000_pdpd[624];

    auto g_y_0_0_0_z_yz_z_xy = buffer_1000_pdpd[625];

    auto g_y_0_0_0_z_yz_z_xz = buffer_1000_pdpd[626];

    auto g_y_0_0_0_z_yz_z_yy = buffer_1000_pdpd[627];

    auto g_y_0_0_0_z_yz_z_yz = buffer_1000_pdpd[628];

    auto g_y_0_0_0_z_yz_z_zz = buffer_1000_pdpd[629];

    auto g_y_0_0_0_z_zz_x_xx = buffer_1000_pdpd[630];

    auto g_y_0_0_0_z_zz_x_xy = buffer_1000_pdpd[631];

    auto g_y_0_0_0_z_zz_x_xz = buffer_1000_pdpd[632];

    auto g_y_0_0_0_z_zz_x_yy = buffer_1000_pdpd[633];

    auto g_y_0_0_0_z_zz_x_yz = buffer_1000_pdpd[634];

    auto g_y_0_0_0_z_zz_x_zz = buffer_1000_pdpd[635];

    auto g_y_0_0_0_z_zz_y_xx = buffer_1000_pdpd[636];

    auto g_y_0_0_0_z_zz_y_xy = buffer_1000_pdpd[637];

    auto g_y_0_0_0_z_zz_y_xz = buffer_1000_pdpd[638];

    auto g_y_0_0_0_z_zz_y_yy = buffer_1000_pdpd[639];

    auto g_y_0_0_0_z_zz_y_yz = buffer_1000_pdpd[640];

    auto g_y_0_0_0_z_zz_y_zz = buffer_1000_pdpd[641];

    auto g_y_0_0_0_z_zz_z_xx = buffer_1000_pdpd[642];

    auto g_y_0_0_0_z_zz_z_xy = buffer_1000_pdpd[643];

    auto g_y_0_0_0_z_zz_z_xz = buffer_1000_pdpd[644];

    auto g_y_0_0_0_z_zz_z_yy = buffer_1000_pdpd[645];

    auto g_y_0_0_0_z_zz_z_yz = buffer_1000_pdpd[646];

    auto g_y_0_0_0_z_zz_z_zz = buffer_1000_pdpd[647];

    auto g_z_0_0_0_x_xx_x_xx = buffer_1000_pdpd[648];

    auto g_z_0_0_0_x_xx_x_xy = buffer_1000_pdpd[649];

    auto g_z_0_0_0_x_xx_x_xz = buffer_1000_pdpd[650];

    auto g_z_0_0_0_x_xx_x_yy = buffer_1000_pdpd[651];

    auto g_z_0_0_0_x_xx_x_yz = buffer_1000_pdpd[652];

    auto g_z_0_0_0_x_xx_x_zz = buffer_1000_pdpd[653];

    auto g_z_0_0_0_x_xx_y_xx = buffer_1000_pdpd[654];

    auto g_z_0_0_0_x_xx_y_xy = buffer_1000_pdpd[655];

    auto g_z_0_0_0_x_xx_y_xz = buffer_1000_pdpd[656];

    auto g_z_0_0_0_x_xx_y_yy = buffer_1000_pdpd[657];

    auto g_z_0_0_0_x_xx_y_yz = buffer_1000_pdpd[658];

    auto g_z_0_0_0_x_xx_y_zz = buffer_1000_pdpd[659];

    auto g_z_0_0_0_x_xx_z_xx = buffer_1000_pdpd[660];

    auto g_z_0_0_0_x_xx_z_xy = buffer_1000_pdpd[661];

    auto g_z_0_0_0_x_xx_z_xz = buffer_1000_pdpd[662];

    auto g_z_0_0_0_x_xx_z_yy = buffer_1000_pdpd[663];

    auto g_z_0_0_0_x_xx_z_yz = buffer_1000_pdpd[664];

    auto g_z_0_0_0_x_xx_z_zz = buffer_1000_pdpd[665];

    auto g_z_0_0_0_x_xy_x_xx = buffer_1000_pdpd[666];

    auto g_z_0_0_0_x_xy_x_xy = buffer_1000_pdpd[667];

    auto g_z_0_0_0_x_xy_x_xz = buffer_1000_pdpd[668];

    auto g_z_0_0_0_x_xy_x_yy = buffer_1000_pdpd[669];

    auto g_z_0_0_0_x_xy_x_yz = buffer_1000_pdpd[670];

    auto g_z_0_0_0_x_xy_x_zz = buffer_1000_pdpd[671];

    auto g_z_0_0_0_x_xy_y_xx = buffer_1000_pdpd[672];

    auto g_z_0_0_0_x_xy_y_xy = buffer_1000_pdpd[673];

    auto g_z_0_0_0_x_xy_y_xz = buffer_1000_pdpd[674];

    auto g_z_0_0_0_x_xy_y_yy = buffer_1000_pdpd[675];

    auto g_z_0_0_0_x_xy_y_yz = buffer_1000_pdpd[676];

    auto g_z_0_0_0_x_xy_y_zz = buffer_1000_pdpd[677];

    auto g_z_0_0_0_x_xy_z_xx = buffer_1000_pdpd[678];

    auto g_z_0_0_0_x_xy_z_xy = buffer_1000_pdpd[679];

    auto g_z_0_0_0_x_xy_z_xz = buffer_1000_pdpd[680];

    auto g_z_0_0_0_x_xy_z_yy = buffer_1000_pdpd[681];

    auto g_z_0_0_0_x_xy_z_yz = buffer_1000_pdpd[682];

    auto g_z_0_0_0_x_xy_z_zz = buffer_1000_pdpd[683];

    auto g_z_0_0_0_x_xz_x_xx = buffer_1000_pdpd[684];

    auto g_z_0_0_0_x_xz_x_xy = buffer_1000_pdpd[685];

    auto g_z_0_0_0_x_xz_x_xz = buffer_1000_pdpd[686];

    auto g_z_0_0_0_x_xz_x_yy = buffer_1000_pdpd[687];

    auto g_z_0_0_0_x_xz_x_yz = buffer_1000_pdpd[688];

    auto g_z_0_0_0_x_xz_x_zz = buffer_1000_pdpd[689];

    auto g_z_0_0_0_x_xz_y_xx = buffer_1000_pdpd[690];

    auto g_z_0_0_0_x_xz_y_xy = buffer_1000_pdpd[691];

    auto g_z_0_0_0_x_xz_y_xz = buffer_1000_pdpd[692];

    auto g_z_0_0_0_x_xz_y_yy = buffer_1000_pdpd[693];

    auto g_z_0_0_0_x_xz_y_yz = buffer_1000_pdpd[694];

    auto g_z_0_0_0_x_xz_y_zz = buffer_1000_pdpd[695];

    auto g_z_0_0_0_x_xz_z_xx = buffer_1000_pdpd[696];

    auto g_z_0_0_0_x_xz_z_xy = buffer_1000_pdpd[697];

    auto g_z_0_0_0_x_xz_z_xz = buffer_1000_pdpd[698];

    auto g_z_0_0_0_x_xz_z_yy = buffer_1000_pdpd[699];

    auto g_z_0_0_0_x_xz_z_yz = buffer_1000_pdpd[700];

    auto g_z_0_0_0_x_xz_z_zz = buffer_1000_pdpd[701];

    auto g_z_0_0_0_x_yy_x_xx = buffer_1000_pdpd[702];

    auto g_z_0_0_0_x_yy_x_xy = buffer_1000_pdpd[703];

    auto g_z_0_0_0_x_yy_x_xz = buffer_1000_pdpd[704];

    auto g_z_0_0_0_x_yy_x_yy = buffer_1000_pdpd[705];

    auto g_z_0_0_0_x_yy_x_yz = buffer_1000_pdpd[706];

    auto g_z_0_0_0_x_yy_x_zz = buffer_1000_pdpd[707];

    auto g_z_0_0_0_x_yy_y_xx = buffer_1000_pdpd[708];

    auto g_z_0_0_0_x_yy_y_xy = buffer_1000_pdpd[709];

    auto g_z_0_0_0_x_yy_y_xz = buffer_1000_pdpd[710];

    auto g_z_0_0_0_x_yy_y_yy = buffer_1000_pdpd[711];

    auto g_z_0_0_0_x_yy_y_yz = buffer_1000_pdpd[712];

    auto g_z_0_0_0_x_yy_y_zz = buffer_1000_pdpd[713];

    auto g_z_0_0_0_x_yy_z_xx = buffer_1000_pdpd[714];

    auto g_z_0_0_0_x_yy_z_xy = buffer_1000_pdpd[715];

    auto g_z_0_0_0_x_yy_z_xz = buffer_1000_pdpd[716];

    auto g_z_0_0_0_x_yy_z_yy = buffer_1000_pdpd[717];

    auto g_z_0_0_0_x_yy_z_yz = buffer_1000_pdpd[718];

    auto g_z_0_0_0_x_yy_z_zz = buffer_1000_pdpd[719];

    auto g_z_0_0_0_x_yz_x_xx = buffer_1000_pdpd[720];

    auto g_z_0_0_0_x_yz_x_xy = buffer_1000_pdpd[721];

    auto g_z_0_0_0_x_yz_x_xz = buffer_1000_pdpd[722];

    auto g_z_0_0_0_x_yz_x_yy = buffer_1000_pdpd[723];

    auto g_z_0_0_0_x_yz_x_yz = buffer_1000_pdpd[724];

    auto g_z_0_0_0_x_yz_x_zz = buffer_1000_pdpd[725];

    auto g_z_0_0_0_x_yz_y_xx = buffer_1000_pdpd[726];

    auto g_z_0_0_0_x_yz_y_xy = buffer_1000_pdpd[727];

    auto g_z_0_0_0_x_yz_y_xz = buffer_1000_pdpd[728];

    auto g_z_0_0_0_x_yz_y_yy = buffer_1000_pdpd[729];

    auto g_z_0_0_0_x_yz_y_yz = buffer_1000_pdpd[730];

    auto g_z_0_0_0_x_yz_y_zz = buffer_1000_pdpd[731];

    auto g_z_0_0_0_x_yz_z_xx = buffer_1000_pdpd[732];

    auto g_z_0_0_0_x_yz_z_xy = buffer_1000_pdpd[733];

    auto g_z_0_0_0_x_yz_z_xz = buffer_1000_pdpd[734];

    auto g_z_0_0_0_x_yz_z_yy = buffer_1000_pdpd[735];

    auto g_z_0_0_0_x_yz_z_yz = buffer_1000_pdpd[736];

    auto g_z_0_0_0_x_yz_z_zz = buffer_1000_pdpd[737];

    auto g_z_0_0_0_x_zz_x_xx = buffer_1000_pdpd[738];

    auto g_z_0_0_0_x_zz_x_xy = buffer_1000_pdpd[739];

    auto g_z_0_0_0_x_zz_x_xz = buffer_1000_pdpd[740];

    auto g_z_0_0_0_x_zz_x_yy = buffer_1000_pdpd[741];

    auto g_z_0_0_0_x_zz_x_yz = buffer_1000_pdpd[742];

    auto g_z_0_0_0_x_zz_x_zz = buffer_1000_pdpd[743];

    auto g_z_0_0_0_x_zz_y_xx = buffer_1000_pdpd[744];

    auto g_z_0_0_0_x_zz_y_xy = buffer_1000_pdpd[745];

    auto g_z_0_0_0_x_zz_y_xz = buffer_1000_pdpd[746];

    auto g_z_0_0_0_x_zz_y_yy = buffer_1000_pdpd[747];

    auto g_z_0_0_0_x_zz_y_yz = buffer_1000_pdpd[748];

    auto g_z_0_0_0_x_zz_y_zz = buffer_1000_pdpd[749];

    auto g_z_0_0_0_x_zz_z_xx = buffer_1000_pdpd[750];

    auto g_z_0_0_0_x_zz_z_xy = buffer_1000_pdpd[751];

    auto g_z_0_0_0_x_zz_z_xz = buffer_1000_pdpd[752];

    auto g_z_0_0_0_x_zz_z_yy = buffer_1000_pdpd[753];

    auto g_z_0_0_0_x_zz_z_yz = buffer_1000_pdpd[754];

    auto g_z_0_0_0_x_zz_z_zz = buffer_1000_pdpd[755];

    auto g_z_0_0_0_y_xx_x_xx = buffer_1000_pdpd[756];

    auto g_z_0_0_0_y_xx_x_xy = buffer_1000_pdpd[757];

    auto g_z_0_0_0_y_xx_x_xz = buffer_1000_pdpd[758];

    auto g_z_0_0_0_y_xx_x_yy = buffer_1000_pdpd[759];

    auto g_z_0_0_0_y_xx_x_yz = buffer_1000_pdpd[760];

    auto g_z_0_0_0_y_xx_x_zz = buffer_1000_pdpd[761];

    auto g_z_0_0_0_y_xx_y_xx = buffer_1000_pdpd[762];

    auto g_z_0_0_0_y_xx_y_xy = buffer_1000_pdpd[763];

    auto g_z_0_0_0_y_xx_y_xz = buffer_1000_pdpd[764];

    auto g_z_0_0_0_y_xx_y_yy = buffer_1000_pdpd[765];

    auto g_z_0_0_0_y_xx_y_yz = buffer_1000_pdpd[766];

    auto g_z_0_0_0_y_xx_y_zz = buffer_1000_pdpd[767];

    auto g_z_0_0_0_y_xx_z_xx = buffer_1000_pdpd[768];

    auto g_z_0_0_0_y_xx_z_xy = buffer_1000_pdpd[769];

    auto g_z_0_0_0_y_xx_z_xz = buffer_1000_pdpd[770];

    auto g_z_0_0_0_y_xx_z_yy = buffer_1000_pdpd[771];

    auto g_z_0_0_0_y_xx_z_yz = buffer_1000_pdpd[772];

    auto g_z_0_0_0_y_xx_z_zz = buffer_1000_pdpd[773];

    auto g_z_0_0_0_y_xy_x_xx = buffer_1000_pdpd[774];

    auto g_z_0_0_0_y_xy_x_xy = buffer_1000_pdpd[775];

    auto g_z_0_0_0_y_xy_x_xz = buffer_1000_pdpd[776];

    auto g_z_0_0_0_y_xy_x_yy = buffer_1000_pdpd[777];

    auto g_z_0_0_0_y_xy_x_yz = buffer_1000_pdpd[778];

    auto g_z_0_0_0_y_xy_x_zz = buffer_1000_pdpd[779];

    auto g_z_0_0_0_y_xy_y_xx = buffer_1000_pdpd[780];

    auto g_z_0_0_0_y_xy_y_xy = buffer_1000_pdpd[781];

    auto g_z_0_0_0_y_xy_y_xz = buffer_1000_pdpd[782];

    auto g_z_0_0_0_y_xy_y_yy = buffer_1000_pdpd[783];

    auto g_z_0_0_0_y_xy_y_yz = buffer_1000_pdpd[784];

    auto g_z_0_0_0_y_xy_y_zz = buffer_1000_pdpd[785];

    auto g_z_0_0_0_y_xy_z_xx = buffer_1000_pdpd[786];

    auto g_z_0_0_0_y_xy_z_xy = buffer_1000_pdpd[787];

    auto g_z_0_0_0_y_xy_z_xz = buffer_1000_pdpd[788];

    auto g_z_0_0_0_y_xy_z_yy = buffer_1000_pdpd[789];

    auto g_z_0_0_0_y_xy_z_yz = buffer_1000_pdpd[790];

    auto g_z_0_0_0_y_xy_z_zz = buffer_1000_pdpd[791];

    auto g_z_0_0_0_y_xz_x_xx = buffer_1000_pdpd[792];

    auto g_z_0_0_0_y_xz_x_xy = buffer_1000_pdpd[793];

    auto g_z_0_0_0_y_xz_x_xz = buffer_1000_pdpd[794];

    auto g_z_0_0_0_y_xz_x_yy = buffer_1000_pdpd[795];

    auto g_z_0_0_0_y_xz_x_yz = buffer_1000_pdpd[796];

    auto g_z_0_0_0_y_xz_x_zz = buffer_1000_pdpd[797];

    auto g_z_0_0_0_y_xz_y_xx = buffer_1000_pdpd[798];

    auto g_z_0_0_0_y_xz_y_xy = buffer_1000_pdpd[799];

    auto g_z_0_0_0_y_xz_y_xz = buffer_1000_pdpd[800];

    auto g_z_0_0_0_y_xz_y_yy = buffer_1000_pdpd[801];

    auto g_z_0_0_0_y_xz_y_yz = buffer_1000_pdpd[802];

    auto g_z_0_0_0_y_xz_y_zz = buffer_1000_pdpd[803];

    auto g_z_0_0_0_y_xz_z_xx = buffer_1000_pdpd[804];

    auto g_z_0_0_0_y_xz_z_xy = buffer_1000_pdpd[805];

    auto g_z_0_0_0_y_xz_z_xz = buffer_1000_pdpd[806];

    auto g_z_0_0_0_y_xz_z_yy = buffer_1000_pdpd[807];

    auto g_z_0_0_0_y_xz_z_yz = buffer_1000_pdpd[808];

    auto g_z_0_0_0_y_xz_z_zz = buffer_1000_pdpd[809];

    auto g_z_0_0_0_y_yy_x_xx = buffer_1000_pdpd[810];

    auto g_z_0_0_0_y_yy_x_xy = buffer_1000_pdpd[811];

    auto g_z_0_0_0_y_yy_x_xz = buffer_1000_pdpd[812];

    auto g_z_0_0_0_y_yy_x_yy = buffer_1000_pdpd[813];

    auto g_z_0_0_0_y_yy_x_yz = buffer_1000_pdpd[814];

    auto g_z_0_0_0_y_yy_x_zz = buffer_1000_pdpd[815];

    auto g_z_0_0_0_y_yy_y_xx = buffer_1000_pdpd[816];

    auto g_z_0_0_0_y_yy_y_xy = buffer_1000_pdpd[817];

    auto g_z_0_0_0_y_yy_y_xz = buffer_1000_pdpd[818];

    auto g_z_0_0_0_y_yy_y_yy = buffer_1000_pdpd[819];

    auto g_z_0_0_0_y_yy_y_yz = buffer_1000_pdpd[820];

    auto g_z_0_0_0_y_yy_y_zz = buffer_1000_pdpd[821];

    auto g_z_0_0_0_y_yy_z_xx = buffer_1000_pdpd[822];

    auto g_z_0_0_0_y_yy_z_xy = buffer_1000_pdpd[823];

    auto g_z_0_0_0_y_yy_z_xz = buffer_1000_pdpd[824];

    auto g_z_0_0_0_y_yy_z_yy = buffer_1000_pdpd[825];

    auto g_z_0_0_0_y_yy_z_yz = buffer_1000_pdpd[826];

    auto g_z_0_0_0_y_yy_z_zz = buffer_1000_pdpd[827];

    auto g_z_0_0_0_y_yz_x_xx = buffer_1000_pdpd[828];

    auto g_z_0_0_0_y_yz_x_xy = buffer_1000_pdpd[829];

    auto g_z_0_0_0_y_yz_x_xz = buffer_1000_pdpd[830];

    auto g_z_0_0_0_y_yz_x_yy = buffer_1000_pdpd[831];

    auto g_z_0_0_0_y_yz_x_yz = buffer_1000_pdpd[832];

    auto g_z_0_0_0_y_yz_x_zz = buffer_1000_pdpd[833];

    auto g_z_0_0_0_y_yz_y_xx = buffer_1000_pdpd[834];

    auto g_z_0_0_0_y_yz_y_xy = buffer_1000_pdpd[835];

    auto g_z_0_0_0_y_yz_y_xz = buffer_1000_pdpd[836];

    auto g_z_0_0_0_y_yz_y_yy = buffer_1000_pdpd[837];

    auto g_z_0_0_0_y_yz_y_yz = buffer_1000_pdpd[838];

    auto g_z_0_0_0_y_yz_y_zz = buffer_1000_pdpd[839];

    auto g_z_0_0_0_y_yz_z_xx = buffer_1000_pdpd[840];

    auto g_z_0_0_0_y_yz_z_xy = buffer_1000_pdpd[841];

    auto g_z_0_0_0_y_yz_z_xz = buffer_1000_pdpd[842];

    auto g_z_0_0_0_y_yz_z_yy = buffer_1000_pdpd[843];

    auto g_z_0_0_0_y_yz_z_yz = buffer_1000_pdpd[844];

    auto g_z_0_0_0_y_yz_z_zz = buffer_1000_pdpd[845];

    auto g_z_0_0_0_y_zz_x_xx = buffer_1000_pdpd[846];

    auto g_z_0_0_0_y_zz_x_xy = buffer_1000_pdpd[847];

    auto g_z_0_0_0_y_zz_x_xz = buffer_1000_pdpd[848];

    auto g_z_0_0_0_y_zz_x_yy = buffer_1000_pdpd[849];

    auto g_z_0_0_0_y_zz_x_yz = buffer_1000_pdpd[850];

    auto g_z_0_0_0_y_zz_x_zz = buffer_1000_pdpd[851];

    auto g_z_0_0_0_y_zz_y_xx = buffer_1000_pdpd[852];

    auto g_z_0_0_0_y_zz_y_xy = buffer_1000_pdpd[853];

    auto g_z_0_0_0_y_zz_y_xz = buffer_1000_pdpd[854];

    auto g_z_0_0_0_y_zz_y_yy = buffer_1000_pdpd[855];

    auto g_z_0_0_0_y_zz_y_yz = buffer_1000_pdpd[856];

    auto g_z_0_0_0_y_zz_y_zz = buffer_1000_pdpd[857];

    auto g_z_0_0_0_y_zz_z_xx = buffer_1000_pdpd[858];

    auto g_z_0_0_0_y_zz_z_xy = buffer_1000_pdpd[859];

    auto g_z_0_0_0_y_zz_z_xz = buffer_1000_pdpd[860];

    auto g_z_0_0_0_y_zz_z_yy = buffer_1000_pdpd[861];

    auto g_z_0_0_0_y_zz_z_yz = buffer_1000_pdpd[862];

    auto g_z_0_0_0_y_zz_z_zz = buffer_1000_pdpd[863];

    auto g_z_0_0_0_z_xx_x_xx = buffer_1000_pdpd[864];

    auto g_z_0_0_0_z_xx_x_xy = buffer_1000_pdpd[865];

    auto g_z_0_0_0_z_xx_x_xz = buffer_1000_pdpd[866];

    auto g_z_0_0_0_z_xx_x_yy = buffer_1000_pdpd[867];

    auto g_z_0_0_0_z_xx_x_yz = buffer_1000_pdpd[868];

    auto g_z_0_0_0_z_xx_x_zz = buffer_1000_pdpd[869];

    auto g_z_0_0_0_z_xx_y_xx = buffer_1000_pdpd[870];

    auto g_z_0_0_0_z_xx_y_xy = buffer_1000_pdpd[871];

    auto g_z_0_0_0_z_xx_y_xz = buffer_1000_pdpd[872];

    auto g_z_0_0_0_z_xx_y_yy = buffer_1000_pdpd[873];

    auto g_z_0_0_0_z_xx_y_yz = buffer_1000_pdpd[874];

    auto g_z_0_0_0_z_xx_y_zz = buffer_1000_pdpd[875];

    auto g_z_0_0_0_z_xx_z_xx = buffer_1000_pdpd[876];

    auto g_z_0_0_0_z_xx_z_xy = buffer_1000_pdpd[877];

    auto g_z_0_0_0_z_xx_z_xz = buffer_1000_pdpd[878];

    auto g_z_0_0_0_z_xx_z_yy = buffer_1000_pdpd[879];

    auto g_z_0_0_0_z_xx_z_yz = buffer_1000_pdpd[880];

    auto g_z_0_0_0_z_xx_z_zz = buffer_1000_pdpd[881];

    auto g_z_0_0_0_z_xy_x_xx = buffer_1000_pdpd[882];

    auto g_z_0_0_0_z_xy_x_xy = buffer_1000_pdpd[883];

    auto g_z_0_0_0_z_xy_x_xz = buffer_1000_pdpd[884];

    auto g_z_0_0_0_z_xy_x_yy = buffer_1000_pdpd[885];

    auto g_z_0_0_0_z_xy_x_yz = buffer_1000_pdpd[886];

    auto g_z_0_0_0_z_xy_x_zz = buffer_1000_pdpd[887];

    auto g_z_0_0_0_z_xy_y_xx = buffer_1000_pdpd[888];

    auto g_z_0_0_0_z_xy_y_xy = buffer_1000_pdpd[889];

    auto g_z_0_0_0_z_xy_y_xz = buffer_1000_pdpd[890];

    auto g_z_0_0_0_z_xy_y_yy = buffer_1000_pdpd[891];

    auto g_z_0_0_0_z_xy_y_yz = buffer_1000_pdpd[892];

    auto g_z_0_0_0_z_xy_y_zz = buffer_1000_pdpd[893];

    auto g_z_0_0_0_z_xy_z_xx = buffer_1000_pdpd[894];

    auto g_z_0_0_0_z_xy_z_xy = buffer_1000_pdpd[895];

    auto g_z_0_0_0_z_xy_z_xz = buffer_1000_pdpd[896];

    auto g_z_0_0_0_z_xy_z_yy = buffer_1000_pdpd[897];

    auto g_z_0_0_0_z_xy_z_yz = buffer_1000_pdpd[898];

    auto g_z_0_0_0_z_xy_z_zz = buffer_1000_pdpd[899];

    auto g_z_0_0_0_z_xz_x_xx = buffer_1000_pdpd[900];

    auto g_z_0_0_0_z_xz_x_xy = buffer_1000_pdpd[901];

    auto g_z_0_0_0_z_xz_x_xz = buffer_1000_pdpd[902];

    auto g_z_0_0_0_z_xz_x_yy = buffer_1000_pdpd[903];

    auto g_z_0_0_0_z_xz_x_yz = buffer_1000_pdpd[904];

    auto g_z_0_0_0_z_xz_x_zz = buffer_1000_pdpd[905];

    auto g_z_0_0_0_z_xz_y_xx = buffer_1000_pdpd[906];

    auto g_z_0_0_0_z_xz_y_xy = buffer_1000_pdpd[907];

    auto g_z_0_0_0_z_xz_y_xz = buffer_1000_pdpd[908];

    auto g_z_0_0_0_z_xz_y_yy = buffer_1000_pdpd[909];

    auto g_z_0_0_0_z_xz_y_yz = buffer_1000_pdpd[910];

    auto g_z_0_0_0_z_xz_y_zz = buffer_1000_pdpd[911];

    auto g_z_0_0_0_z_xz_z_xx = buffer_1000_pdpd[912];

    auto g_z_0_0_0_z_xz_z_xy = buffer_1000_pdpd[913];

    auto g_z_0_0_0_z_xz_z_xz = buffer_1000_pdpd[914];

    auto g_z_0_0_0_z_xz_z_yy = buffer_1000_pdpd[915];

    auto g_z_0_0_0_z_xz_z_yz = buffer_1000_pdpd[916];

    auto g_z_0_0_0_z_xz_z_zz = buffer_1000_pdpd[917];

    auto g_z_0_0_0_z_yy_x_xx = buffer_1000_pdpd[918];

    auto g_z_0_0_0_z_yy_x_xy = buffer_1000_pdpd[919];

    auto g_z_0_0_0_z_yy_x_xz = buffer_1000_pdpd[920];

    auto g_z_0_0_0_z_yy_x_yy = buffer_1000_pdpd[921];

    auto g_z_0_0_0_z_yy_x_yz = buffer_1000_pdpd[922];

    auto g_z_0_0_0_z_yy_x_zz = buffer_1000_pdpd[923];

    auto g_z_0_0_0_z_yy_y_xx = buffer_1000_pdpd[924];

    auto g_z_0_0_0_z_yy_y_xy = buffer_1000_pdpd[925];

    auto g_z_0_0_0_z_yy_y_xz = buffer_1000_pdpd[926];

    auto g_z_0_0_0_z_yy_y_yy = buffer_1000_pdpd[927];

    auto g_z_0_0_0_z_yy_y_yz = buffer_1000_pdpd[928];

    auto g_z_0_0_0_z_yy_y_zz = buffer_1000_pdpd[929];

    auto g_z_0_0_0_z_yy_z_xx = buffer_1000_pdpd[930];

    auto g_z_0_0_0_z_yy_z_xy = buffer_1000_pdpd[931];

    auto g_z_0_0_0_z_yy_z_xz = buffer_1000_pdpd[932];

    auto g_z_0_0_0_z_yy_z_yy = buffer_1000_pdpd[933];

    auto g_z_0_0_0_z_yy_z_yz = buffer_1000_pdpd[934];

    auto g_z_0_0_0_z_yy_z_zz = buffer_1000_pdpd[935];

    auto g_z_0_0_0_z_yz_x_xx = buffer_1000_pdpd[936];

    auto g_z_0_0_0_z_yz_x_xy = buffer_1000_pdpd[937];

    auto g_z_0_0_0_z_yz_x_xz = buffer_1000_pdpd[938];

    auto g_z_0_0_0_z_yz_x_yy = buffer_1000_pdpd[939];

    auto g_z_0_0_0_z_yz_x_yz = buffer_1000_pdpd[940];

    auto g_z_0_0_0_z_yz_x_zz = buffer_1000_pdpd[941];

    auto g_z_0_0_0_z_yz_y_xx = buffer_1000_pdpd[942];

    auto g_z_0_0_0_z_yz_y_xy = buffer_1000_pdpd[943];

    auto g_z_0_0_0_z_yz_y_xz = buffer_1000_pdpd[944];

    auto g_z_0_0_0_z_yz_y_yy = buffer_1000_pdpd[945];

    auto g_z_0_0_0_z_yz_y_yz = buffer_1000_pdpd[946];

    auto g_z_0_0_0_z_yz_y_zz = buffer_1000_pdpd[947];

    auto g_z_0_0_0_z_yz_z_xx = buffer_1000_pdpd[948];

    auto g_z_0_0_0_z_yz_z_xy = buffer_1000_pdpd[949];

    auto g_z_0_0_0_z_yz_z_xz = buffer_1000_pdpd[950];

    auto g_z_0_0_0_z_yz_z_yy = buffer_1000_pdpd[951];

    auto g_z_0_0_0_z_yz_z_yz = buffer_1000_pdpd[952];

    auto g_z_0_0_0_z_yz_z_zz = buffer_1000_pdpd[953];

    auto g_z_0_0_0_z_zz_x_xx = buffer_1000_pdpd[954];

    auto g_z_0_0_0_z_zz_x_xy = buffer_1000_pdpd[955];

    auto g_z_0_0_0_z_zz_x_xz = buffer_1000_pdpd[956];

    auto g_z_0_0_0_z_zz_x_yy = buffer_1000_pdpd[957];

    auto g_z_0_0_0_z_zz_x_yz = buffer_1000_pdpd[958];

    auto g_z_0_0_0_z_zz_x_zz = buffer_1000_pdpd[959];

    auto g_z_0_0_0_z_zz_y_xx = buffer_1000_pdpd[960];

    auto g_z_0_0_0_z_zz_y_xy = buffer_1000_pdpd[961];

    auto g_z_0_0_0_z_zz_y_xz = buffer_1000_pdpd[962];

    auto g_z_0_0_0_z_zz_y_yy = buffer_1000_pdpd[963];

    auto g_z_0_0_0_z_zz_y_yz = buffer_1000_pdpd[964];

    auto g_z_0_0_0_z_zz_y_zz = buffer_1000_pdpd[965];

    auto g_z_0_0_0_z_zz_z_xx = buffer_1000_pdpd[966];

    auto g_z_0_0_0_z_zz_z_xy = buffer_1000_pdpd[967];

    auto g_z_0_0_0_z_zz_z_xz = buffer_1000_pdpd[968];

    auto g_z_0_0_0_z_zz_z_yy = buffer_1000_pdpd[969];

    auto g_z_0_0_0_z_zz_z_yz = buffer_1000_pdpd[970];

    auto g_z_0_0_0_z_zz_z_zz = buffer_1000_pdpd[971];

    // integrals block (0-6)

    #pragma omp simd aligned(g_0_xx_x_xx, g_0_xx_x_xy, g_0_xx_x_xz, g_0_xx_x_yy, g_0_xx_x_yz, g_0_xx_x_zz, g_x_0_0_0_x_xx_x_xx, g_x_0_0_0_x_xx_x_xy, g_x_0_0_0_x_xx_x_xz, g_x_0_0_0_x_xx_x_yy, g_x_0_0_0_x_xx_x_yz, g_x_0_0_0_x_xx_x_zz, g_xx_xx_x_xx, g_xx_xx_x_xy, g_xx_xx_x_xz, g_xx_xx_x_yy, g_xx_xx_x_yz, g_xx_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xx_x_xx[i] = -g_0_xx_x_xx[i] + 2.0 * g_xx_xx_x_xx[i] * a_exp;

        g_x_0_0_0_x_xx_x_xy[i] = -g_0_xx_x_xy[i] + 2.0 * g_xx_xx_x_xy[i] * a_exp;

        g_x_0_0_0_x_xx_x_xz[i] = -g_0_xx_x_xz[i] + 2.0 * g_xx_xx_x_xz[i] * a_exp;

        g_x_0_0_0_x_xx_x_yy[i] = -g_0_xx_x_yy[i] + 2.0 * g_xx_xx_x_yy[i] * a_exp;

        g_x_0_0_0_x_xx_x_yz[i] = -g_0_xx_x_yz[i] + 2.0 * g_xx_xx_x_yz[i] * a_exp;

        g_x_0_0_0_x_xx_x_zz[i] = -g_0_xx_x_zz[i] + 2.0 * g_xx_xx_x_zz[i] * a_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_0_xx_y_xx, g_0_xx_y_xy, g_0_xx_y_xz, g_0_xx_y_yy, g_0_xx_y_yz, g_0_xx_y_zz, g_x_0_0_0_x_xx_y_xx, g_x_0_0_0_x_xx_y_xy, g_x_0_0_0_x_xx_y_xz, g_x_0_0_0_x_xx_y_yy, g_x_0_0_0_x_xx_y_yz, g_x_0_0_0_x_xx_y_zz, g_xx_xx_y_xx, g_xx_xx_y_xy, g_xx_xx_y_xz, g_xx_xx_y_yy, g_xx_xx_y_yz, g_xx_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xx_y_xx[i] = -g_0_xx_y_xx[i] + 2.0 * g_xx_xx_y_xx[i] * a_exp;

        g_x_0_0_0_x_xx_y_xy[i] = -g_0_xx_y_xy[i] + 2.0 * g_xx_xx_y_xy[i] * a_exp;

        g_x_0_0_0_x_xx_y_xz[i] = -g_0_xx_y_xz[i] + 2.0 * g_xx_xx_y_xz[i] * a_exp;

        g_x_0_0_0_x_xx_y_yy[i] = -g_0_xx_y_yy[i] + 2.0 * g_xx_xx_y_yy[i] * a_exp;

        g_x_0_0_0_x_xx_y_yz[i] = -g_0_xx_y_yz[i] + 2.0 * g_xx_xx_y_yz[i] * a_exp;

        g_x_0_0_0_x_xx_y_zz[i] = -g_0_xx_y_zz[i] + 2.0 * g_xx_xx_y_zz[i] * a_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_0_xx_z_xx, g_0_xx_z_xy, g_0_xx_z_xz, g_0_xx_z_yy, g_0_xx_z_yz, g_0_xx_z_zz, g_x_0_0_0_x_xx_z_xx, g_x_0_0_0_x_xx_z_xy, g_x_0_0_0_x_xx_z_xz, g_x_0_0_0_x_xx_z_yy, g_x_0_0_0_x_xx_z_yz, g_x_0_0_0_x_xx_z_zz, g_xx_xx_z_xx, g_xx_xx_z_xy, g_xx_xx_z_xz, g_xx_xx_z_yy, g_xx_xx_z_yz, g_xx_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xx_z_xx[i] = -g_0_xx_z_xx[i] + 2.0 * g_xx_xx_z_xx[i] * a_exp;

        g_x_0_0_0_x_xx_z_xy[i] = -g_0_xx_z_xy[i] + 2.0 * g_xx_xx_z_xy[i] * a_exp;

        g_x_0_0_0_x_xx_z_xz[i] = -g_0_xx_z_xz[i] + 2.0 * g_xx_xx_z_xz[i] * a_exp;

        g_x_0_0_0_x_xx_z_yy[i] = -g_0_xx_z_yy[i] + 2.0 * g_xx_xx_z_yy[i] * a_exp;

        g_x_0_0_0_x_xx_z_yz[i] = -g_0_xx_z_yz[i] + 2.0 * g_xx_xx_z_yz[i] * a_exp;

        g_x_0_0_0_x_xx_z_zz[i] = -g_0_xx_z_zz[i] + 2.0 * g_xx_xx_z_zz[i] * a_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_0_xy_x_xx, g_0_xy_x_xy, g_0_xy_x_xz, g_0_xy_x_yy, g_0_xy_x_yz, g_0_xy_x_zz, g_x_0_0_0_x_xy_x_xx, g_x_0_0_0_x_xy_x_xy, g_x_0_0_0_x_xy_x_xz, g_x_0_0_0_x_xy_x_yy, g_x_0_0_0_x_xy_x_yz, g_x_0_0_0_x_xy_x_zz, g_xx_xy_x_xx, g_xx_xy_x_xy, g_xx_xy_x_xz, g_xx_xy_x_yy, g_xx_xy_x_yz, g_xx_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xy_x_xx[i] = -g_0_xy_x_xx[i] + 2.0 * g_xx_xy_x_xx[i] * a_exp;

        g_x_0_0_0_x_xy_x_xy[i] = -g_0_xy_x_xy[i] + 2.0 * g_xx_xy_x_xy[i] * a_exp;

        g_x_0_0_0_x_xy_x_xz[i] = -g_0_xy_x_xz[i] + 2.0 * g_xx_xy_x_xz[i] * a_exp;

        g_x_0_0_0_x_xy_x_yy[i] = -g_0_xy_x_yy[i] + 2.0 * g_xx_xy_x_yy[i] * a_exp;

        g_x_0_0_0_x_xy_x_yz[i] = -g_0_xy_x_yz[i] + 2.0 * g_xx_xy_x_yz[i] * a_exp;

        g_x_0_0_0_x_xy_x_zz[i] = -g_0_xy_x_zz[i] + 2.0 * g_xx_xy_x_zz[i] * a_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_0_xy_y_xx, g_0_xy_y_xy, g_0_xy_y_xz, g_0_xy_y_yy, g_0_xy_y_yz, g_0_xy_y_zz, g_x_0_0_0_x_xy_y_xx, g_x_0_0_0_x_xy_y_xy, g_x_0_0_0_x_xy_y_xz, g_x_0_0_0_x_xy_y_yy, g_x_0_0_0_x_xy_y_yz, g_x_0_0_0_x_xy_y_zz, g_xx_xy_y_xx, g_xx_xy_y_xy, g_xx_xy_y_xz, g_xx_xy_y_yy, g_xx_xy_y_yz, g_xx_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xy_y_xx[i] = -g_0_xy_y_xx[i] + 2.0 * g_xx_xy_y_xx[i] * a_exp;

        g_x_0_0_0_x_xy_y_xy[i] = -g_0_xy_y_xy[i] + 2.0 * g_xx_xy_y_xy[i] * a_exp;

        g_x_0_0_0_x_xy_y_xz[i] = -g_0_xy_y_xz[i] + 2.0 * g_xx_xy_y_xz[i] * a_exp;

        g_x_0_0_0_x_xy_y_yy[i] = -g_0_xy_y_yy[i] + 2.0 * g_xx_xy_y_yy[i] * a_exp;

        g_x_0_0_0_x_xy_y_yz[i] = -g_0_xy_y_yz[i] + 2.0 * g_xx_xy_y_yz[i] * a_exp;

        g_x_0_0_0_x_xy_y_zz[i] = -g_0_xy_y_zz[i] + 2.0 * g_xx_xy_y_zz[i] * a_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_0_xy_z_xx, g_0_xy_z_xy, g_0_xy_z_xz, g_0_xy_z_yy, g_0_xy_z_yz, g_0_xy_z_zz, g_x_0_0_0_x_xy_z_xx, g_x_0_0_0_x_xy_z_xy, g_x_0_0_0_x_xy_z_xz, g_x_0_0_0_x_xy_z_yy, g_x_0_0_0_x_xy_z_yz, g_x_0_0_0_x_xy_z_zz, g_xx_xy_z_xx, g_xx_xy_z_xy, g_xx_xy_z_xz, g_xx_xy_z_yy, g_xx_xy_z_yz, g_xx_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xy_z_xx[i] = -g_0_xy_z_xx[i] + 2.0 * g_xx_xy_z_xx[i] * a_exp;

        g_x_0_0_0_x_xy_z_xy[i] = -g_0_xy_z_xy[i] + 2.0 * g_xx_xy_z_xy[i] * a_exp;

        g_x_0_0_0_x_xy_z_xz[i] = -g_0_xy_z_xz[i] + 2.0 * g_xx_xy_z_xz[i] * a_exp;

        g_x_0_0_0_x_xy_z_yy[i] = -g_0_xy_z_yy[i] + 2.0 * g_xx_xy_z_yy[i] * a_exp;

        g_x_0_0_0_x_xy_z_yz[i] = -g_0_xy_z_yz[i] + 2.0 * g_xx_xy_z_yz[i] * a_exp;

        g_x_0_0_0_x_xy_z_zz[i] = -g_0_xy_z_zz[i] + 2.0 * g_xx_xy_z_zz[i] * a_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_0_xz_x_xx, g_0_xz_x_xy, g_0_xz_x_xz, g_0_xz_x_yy, g_0_xz_x_yz, g_0_xz_x_zz, g_x_0_0_0_x_xz_x_xx, g_x_0_0_0_x_xz_x_xy, g_x_0_0_0_x_xz_x_xz, g_x_0_0_0_x_xz_x_yy, g_x_0_0_0_x_xz_x_yz, g_x_0_0_0_x_xz_x_zz, g_xx_xz_x_xx, g_xx_xz_x_xy, g_xx_xz_x_xz, g_xx_xz_x_yy, g_xx_xz_x_yz, g_xx_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xz_x_xx[i] = -g_0_xz_x_xx[i] + 2.0 * g_xx_xz_x_xx[i] * a_exp;

        g_x_0_0_0_x_xz_x_xy[i] = -g_0_xz_x_xy[i] + 2.0 * g_xx_xz_x_xy[i] * a_exp;

        g_x_0_0_0_x_xz_x_xz[i] = -g_0_xz_x_xz[i] + 2.0 * g_xx_xz_x_xz[i] * a_exp;

        g_x_0_0_0_x_xz_x_yy[i] = -g_0_xz_x_yy[i] + 2.0 * g_xx_xz_x_yy[i] * a_exp;

        g_x_0_0_0_x_xz_x_yz[i] = -g_0_xz_x_yz[i] + 2.0 * g_xx_xz_x_yz[i] * a_exp;

        g_x_0_0_0_x_xz_x_zz[i] = -g_0_xz_x_zz[i] + 2.0 * g_xx_xz_x_zz[i] * a_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_0_xz_y_xx, g_0_xz_y_xy, g_0_xz_y_xz, g_0_xz_y_yy, g_0_xz_y_yz, g_0_xz_y_zz, g_x_0_0_0_x_xz_y_xx, g_x_0_0_0_x_xz_y_xy, g_x_0_0_0_x_xz_y_xz, g_x_0_0_0_x_xz_y_yy, g_x_0_0_0_x_xz_y_yz, g_x_0_0_0_x_xz_y_zz, g_xx_xz_y_xx, g_xx_xz_y_xy, g_xx_xz_y_xz, g_xx_xz_y_yy, g_xx_xz_y_yz, g_xx_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xz_y_xx[i] = -g_0_xz_y_xx[i] + 2.0 * g_xx_xz_y_xx[i] * a_exp;

        g_x_0_0_0_x_xz_y_xy[i] = -g_0_xz_y_xy[i] + 2.0 * g_xx_xz_y_xy[i] * a_exp;

        g_x_0_0_0_x_xz_y_xz[i] = -g_0_xz_y_xz[i] + 2.0 * g_xx_xz_y_xz[i] * a_exp;

        g_x_0_0_0_x_xz_y_yy[i] = -g_0_xz_y_yy[i] + 2.0 * g_xx_xz_y_yy[i] * a_exp;

        g_x_0_0_0_x_xz_y_yz[i] = -g_0_xz_y_yz[i] + 2.0 * g_xx_xz_y_yz[i] * a_exp;

        g_x_0_0_0_x_xz_y_zz[i] = -g_0_xz_y_zz[i] + 2.0 * g_xx_xz_y_zz[i] * a_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_0_xz_z_xx, g_0_xz_z_xy, g_0_xz_z_xz, g_0_xz_z_yy, g_0_xz_z_yz, g_0_xz_z_zz, g_x_0_0_0_x_xz_z_xx, g_x_0_0_0_x_xz_z_xy, g_x_0_0_0_x_xz_z_xz, g_x_0_0_0_x_xz_z_yy, g_x_0_0_0_x_xz_z_yz, g_x_0_0_0_x_xz_z_zz, g_xx_xz_z_xx, g_xx_xz_z_xy, g_xx_xz_z_xz, g_xx_xz_z_yy, g_xx_xz_z_yz, g_xx_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xz_z_xx[i] = -g_0_xz_z_xx[i] + 2.0 * g_xx_xz_z_xx[i] * a_exp;

        g_x_0_0_0_x_xz_z_xy[i] = -g_0_xz_z_xy[i] + 2.0 * g_xx_xz_z_xy[i] * a_exp;

        g_x_0_0_0_x_xz_z_xz[i] = -g_0_xz_z_xz[i] + 2.0 * g_xx_xz_z_xz[i] * a_exp;

        g_x_0_0_0_x_xz_z_yy[i] = -g_0_xz_z_yy[i] + 2.0 * g_xx_xz_z_yy[i] * a_exp;

        g_x_0_0_0_x_xz_z_yz[i] = -g_0_xz_z_yz[i] + 2.0 * g_xx_xz_z_yz[i] * a_exp;

        g_x_0_0_0_x_xz_z_zz[i] = -g_0_xz_z_zz[i] + 2.0 * g_xx_xz_z_zz[i] * a_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_0_yy_x_xx, g_0_yy_x_xy, g_0_yy_x_xz, g_0_yy_x_yy, g_0_yy_x_yz, g_0_yy_x_zz, g_x_0_0_0_x_yy_x_xx, g_x_0_0_0_x_yy_x_xy, g_x_0_0_0_x_yy_x_xz, g_x_0_0_0_x_yy_x_yy, g_x_0_0_0_x_yy_x_yz, g_x_0_0_0_x_yy_x_zz, g_xx_yy_x_xx, g_xx_yy_x_xy, g_xx_yy_x_xz, g_xx_yy_x_yy, g_xx_yy_x_yz, g_xx_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_yy_x_xx[i] = -g_0_yy_x_xx[i] + 2.0 * g_xx_yy_x_xx[i] * a_exp;

        g_x_0_0_0_x_yy_x_xy[i] = -g_0_yy_x_xy[i] + 2.0 * g_xx_yy_x_xy[i] * a_exp;

        g_x_0_0_0_x_yy_x_xz[i] = -g_0_yy_x_xz[i] + 2.0 * g_xx_yy_x_xz[i] * a_exp;

        g_x_0_0_0_x_yy_x_yy[i] = -g_0_yy_x_yy[i] + 2.0 * g_xx_yy_x_yy[i] * a_exp;

        g_x_0_0_0_x_yy_x_yz[i] = -g_0_yy_x_yz[i] + 2.0 * g_xx_yy_x_yz[i] * a_exp;

        g_x_0_0_0_x_yy_x_zz[i] = -g_0_yy_x_zz[i] + 2.0 * g_xx_yy_x_zz[i] * a_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_0_yy_y_xx, g_0_yy_y_xy, g_0_yy_y_xz, g_0_yy_y_yy, g_0_yy_y_yz, g_0_yy_y_zz, g_x_0_0_0_x_yy_y_xx, g_x_0_0_0_x_yy_y_xy, g_x_0_0_0_x_yy_y_xz, g_x_0_0_0_x_yy_y_yy, g_x_0_0_0_x_yy_y_yz, g_x_0_0_0_x_yy_y_zz, g_xx_yy_y_xx, g_xx_yy_y_xy, g_xx_yy_y_xz, g_xx_yy_y_yy, g_xx_yy_y_yz, g_xx_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_yy_y_xx[i] = -g_0_yy_y_xx[i] + 2.0 * g_xx_yy_y_xx[i] * a_exp;

        g_x_0_0_0_x_yy_y_xy[i] = -g_0_yy_y_xy[i] + 2.0 * g_xx_yy_y_xy[i] * a_exp;

        g_x_0_0_0_x_yy_y_xz[i] = -g_0_yy_y_xz[i] + 2.0 * g_xx_yy_y_xz[i] * a_exp;

        g_x_0_0_0_x_yy_y_yy[i] = -g_0_yy_y_yy[i] + 2.0 * g_xx_yy_y_yy[i] * a_exp;

        g_x_0_0_0_x_yy_y_yz[i] = -g_0_yy_y_yz[i] + 2.0 * g_xx_yy_y_yz[i] * a_exp;

        g_x_0_0_0_x_yy_y_zz[i] = -g_0_yy_y_zz[i] + 2.0 * g_xx_yy_y_zz[i] * a_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_0_yy_z_xx, g_0_yy_z_xy, g_0_yy_z_xz, g_0_yy_z_yy, g_0_yy_z_yz, g_0_yy_z_zz, g_x_0_0_0_x_yy_z_xx, g_x_0_0_0_x_yy_z_xy, g_x_0_0_0_x_yy_z_xz, g_x_0_0_0_x_yy_z_yy, g_x_0_0_0_x_yy_z_yz, g_x_0_0_0_x_yy_z_zz, g_xx_yy_z_xx, g_xx_yy_z_xy, g_xx_yy_z_xz, g_xx_yy_z_yy, g_xx_yy_z_yz, g_xx_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_yy_z_xx[i] = -g_0_yy_z_xx[i] + 2.0 * g_xx_yy_z_xx[i] * a_exp;

        g_x_0_0_0_x_yy_z_xy[i] = -g_0_yy_z_xy[i] + 2.0 * g_xx_yy_z_xy[i] * a_exp;

        g_x_0_0_0_x_yy_z_xz[i] = -g_0_yy_z_xz[i] + 2.0 * g_xx_yy_z_xz[i] * a_exp;

        g_x_0_0_0_x_yy_z_yy[i] = -g_0_yy_z_yy[i] + 2.0 * g_xx_yy_z_yy[i] * a_exp;

        g_x_0_0_0_x_yy_z_yz[i] = -g_0_yy_z_yz[i] + 2.0 * g_xx_yy_z_yz[i] * a_exp;

        g_x_0_0_0_x_yy_z_zz[i] = -g_0_yy_z_zz[i] + 2.0 * g_xx_yy_z_zz[i] * a_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_0_yz_x_xx, g_0_yz_x_xy, g_0_yz_x_xz, g_0_yz_x_yy, g_0_yz_x_yz, g_0_yz_x_zz, g_x_0_0_0_x_yz_x_xx, g_x_0_0_0_x_yz_x_xy, g_x_0_0_0_x_yz_x_xz, g_x_0_0_0_x_yz_x_yy, g_x_0_0_0_x_yz_x_yz, g_x_0_0_0_x_yz_x_zz, g_xx_yz_x_xx, g_xx_yz_x_xy, g_xx_yz_x_xz, g_xx_yz_x_yy, g_xx_yz_x_yz, g_xx_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_yz_x_xx[i] = -g_0_yz_x_xx[i] + 2.0 * g_xx_yz_x_xx[i] * a_exp;

        g_x_0_0_0_x_yz_x_xy[i] = -g_0_yz_x_xy[i] + 2.0 * g_xx_yz_x_xy[i] * a_exp;

        g_x_0_0_0_x_yz_x_xz[i] = -g_0_yz_x_xz[i] + 2.0 * g_xx_yz_x_xz[i] * a_exp;

        g_x_0_0_0_x_yz_x_yy[i] = -g_0_yz_x_yy[i] + 2.0 * g_xx_yz_x_yy[i] * a_exp;

        g_x_0_0_0_x_yz_x_yz[i] = -g_0_yz_x_yz[i] + 2.0 * g_xx_yz_x_yz[i] * a_exp;

        g_x_0_0_0_x_yz_x_zz[i] = -g_0_yz_x_zz[i] + 2.0 * g_xx_yz_x_zz[i] * a_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_0_yz_y_xx, g_0_yz_y_xy, g_0_yz_y_xz, g_0_yz_y_yy, g_0_yz_y_yz, g_0_yz_y_zz, g_x_0_0_0_x_yz_y_xx, g_x_0_0_0_x_yz_y_xy, g_x_0_0_0_x_yz_y_xz, g_x_0_0_0_x_yz_y_yy, g_x_0_0_0_x_yz_y_yz, g_x_0_0_0_x_yz_y_zz, g_xx_yz_y_xx, g_xx_yz_y_xy, g_xx_yz_y_xz, g_xx_yz_y_yy, g_xx_yz_y_yz, g_xx_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_yz_y_xx[i] = -g_0_yz_y_xx[i] + 2.0 * g_xx_yz_y_xx[i] * a_exp;

        g_x_0_0_0_x_yz_y_xy[i] = -g_0_yz_y_xy[i] + 2.0 * g_xx_yz_y_xy[i] * a_exp;

        g_x_0_0_0_x_yz_y_xz[i] = -g_0_yz_y_xz[i] + 2.0 * g_xx_yz_y_xz[i] * a_exp;

        g_x_0_0_0_x_yz_y_yy[i] = -g_0_yz_y_yy[i] + 2.0 * g_xx_yz_y_yy[i] * a_exp;

        g_x_0_0_0_x_yz_y_yz[i] = -g_0_yz_y_yz[i] + 2.0 * g_xx_yz_y_yz[i] * a_exp;

        g_x_0_0_0_x_yz_y_zz[i] = -g_0_yz_y_zz[i] + 2.0 * g_xx_yz_y_zz[i] * a_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_0_yz_z_xx, g_0_yz_z_xy, g_0_yz_z_xz, g_0_yz_z_yy, g_0_yz_z_yz, g_0_yz_z_zz, g_x_0_0_0_x_yz_z_xx, g_x_0_0_0_x_yz_z_xy, g_x_0_0_0_x_yz_z_xz, g_x_0_0_0_x_yz_z_yy, g_x_0_0_0_x_yz_z_yz, g_x_0_0_0_x_yz_z_zz, g_xx_yz_z_xx, g_xx_yz_z_xy, g_xx_yz_z_xz, g_xx_yz_z_yy, g_xx_yz_z_yz, g_xx_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_yz_z_xx[i] = -g_0_yz_z_xx[i] + 2.0 * g_xx_yz_z_xx[i] * a_exp;

        g_x_0_0_0_x_yz_z_xy[i] = -g_0_yz_z_xy[i] + 2.0 * g_xx_yz_z_xy[i] * a_exp;

        g_x_0_0_0_x_yz_z_xz[i] = -g_0_yz_z_xz[i] + 2.0 * g_xx_yz_z_xz[i] * a_exp;

        g_x_0_0_0_x_yz_z_yy[i] = -g_0_yz_z_yy[i] + 2.0 * g_xx_yz_z_yy[i] * a_exp;

        g_x_0_0_0_x_yz_z_yz[i] = -g_0_yz_z_yz[i] + 2.0 * g_xx_yz_z_yz[i] * a_exp;

        g_x_0_0_0_x_yz_z_zz[i] = -g_0_yz_z_zz[i] + 2.0 * g_xx_yz_z_zz[i] * a_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_0_zz_x_xx, g_0_zz_x_xy, g_0_zz_x_xz, g_0_zz_x_yy, g_0_zz_x_yz, g_0_zz_x_zz, g_x_0_0_0_x_zz_x_xx, g_x_0_0_0_x_zz_x_xy, g_x_0_0_0_x_zz_x_xz, g_x_0_0_0_x_zz_x_yy, g_x_0_0_0_x_zz_x_yz, g_x_0_0_0_x_zz_x_zz, g_xx_zz_x_xx, g_xx_zz_x_xy, g_xx_zz_x_xz, g_xx_zz_x_yy, g_xx_zz_x_yz, g_xx_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_zz_x_xx[i] = -g_0_zz_x_xx[i] + 2.0 * g_xx_zz_x_xx[i] * a_exp;

        g_x_0_0_0_x_zz_x_xy[i] = -g_0_zz_x_xy[i] + 2.0 * g_xx_zz_x_xy[i] * a_exp;

        g_x_0_0_0_x_zz_x_xz[i] = -g_0_zz_x_xz[i] + 2.0 * g_xx_zz_x_xz[i] * a_exp;

        g_x_0_0_0_x_zz_x_yy[i] = -g_0_zz_x_yy[i] + 2.0 * g_xx_zz_x_yy[i] * a_exp;

        g_x_0_0_0_x_zz_x_yz[i] = -g_0_zz_x_yz[i] + 2.0 * g_xx_zz_x_yz[i] * a_exp;

        g_x_0_0_0_x_zz_x_zz[i] = -g_0_zz_x_zz[i] + 2.0 * g_xx_zz_x_zz[i] * a_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_0_zz_y_xx, g_0_zz_y_xy, g_0_zz_y_xz, g_0_zz_y_yy, g_0_zz_y_yz, g_0_zz_y_zz, g_x_0_0_0_x_zz_y_xx, g_x_0_0_0_x_zz_y_xy, g_x_0_0_0_x_zz_y_xz, g_x_0_0_0_x_zz_y_yy, g_x_0_0_0_x_zz_y_yz, g_x_0_0_0_x_zz_y_zz, g_xx_zz_y_xx, g_xx_zz_y_xy, g_xx_zz_y_xz, g_xx_zz_y_yy, g_xx_zz_y_yz, g_xx_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_zz_y_xx[i] = -g_0_zz_y_xx[i] + 2.0 * g_xx_zz_y_xx[i] * a_exp;

        g_x_0_0_0_x_zz_y_xy[i] = -g_0_zz_y_xy[i] + 2.0 * g_xx_zz_y_xy[i] * a_exp;

        g_x_0_0_0_x_zz_y_xz[i] = -g_0_zz_y_xz[i] + 2.0 * g_xx_zz_y_xz[i] * a_exp;

        g_x_0_0_0_x_zz_y_yy[i] = -g_0_zz_y_yy[i] + 2.0 * g_xx_zz_y_yy[i] * a_exp;

        g_x_0_0_0_x_zz_y_yz[i] = -g_0_zz_y_yz[i] + 2.0 * g_xx_zz_y_yz[i] * a_exp;

        g_x_0_0_0_x_zz_y_zz[i] = -g_0_zz_y_zz[i] + 2.0 * g_xx_zz_y_zz[i] * a_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_0_zz_z_xx, g_0_zz_z_xy, g_0_zz_z_xz, g_0_zz_z_yy, g_0_zz_z_yz, g_0_zz_z_zz, g_x_0_0_0_x_zz_z_xx, g_x_0_0_0_x_zz_z_xy, g_x_0_0_0_x_zz_z_xz, g_x_0_0_0_x_zz_z_yy, g_x_0_0_0_x_zz_z_yz, g_x_0_0_0_x_zz_z_zz, g_xx_zz_z_xx, g_xx_zz_z_xy, g_xx_zz_z_xz, g_xx_zz_z_yy, g_xx_zz_z_yz, g_xx_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_zz_z_xx[i] = -g_0_zz_z_xx[i] + 2.0 * g_xx_zz_z_xx[i] * a_exp;

        g_x_0_0_0_x_zz_z_xy[i] = -g_0_zz_z_xy[i] + 2.0 * g_xx_zz_z_xy[i] * a_exp;

        g_x_0_0_0_x_zz_z_xz[i] = -g_0_zz_z_xz[i] + 2.0 * g_xx_zz_z_xz[i] * a_exp;

        g_x_0_0_0_x_zz_z_yy[i] = -g_0_zz_z_yy[i] + 2.0 * g_xx_zz_z_yy[i] * a_exp;

        g_x_0_0_0_x_zz_z_yz[i] = -g_0_zz_z_yz[i] + 2.0 * g_xx_zz_z_yz[i] * a_exp;

        g_x_0_0_0_x_zz_z_zz[i] = -g_0_zz_z_zz[i] + 2.0 * g_xx_zz_z_zz[i] * a_exp;
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_x_0_0_0_y_xx_x_xx, g_x_0_0_0_y_xx_x_xy, g_x_0_0_0_y_xx_x_xz, g_x_0_0_0_y_xx_x_yy, g_x_0_0_0_y_xx_x_yz, g_x_0_0_0_y_xx_x_zz, g_xy_xx_x_xx, g_xy_xx_x_xy, g_xy_xx_x_xz, g_xy_xx_x_yy, g_xy_xx_x_yz, g_xy_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xx_x_xx[i] = 2.0 * g_xy_xx_x_xx[i] * a_exp;

        g_x_0_0_0_y_xx_x_xy[i] = 2.0 * g_xy_xx_x_xy[i] * a_exp;

        g_x_0_0_0_y_xx_x_xz[i] = 2.0 * g_xy_xx_x_xz[i] * a_exp;

        g_x_0_0_0_y_xx_x_yy[i] = 2.0 * g_xy_xx_x_yy[i] * a_exp;

        g_x_0_0_0_y_xx_x_yz[i] = 2.0 * g_xy_xx_x_yz[i] * a_exp;

        g_x_0_0_0_y_xx_x_zz[i] = 2.0 * g_xy_xx_x_zz[i] * a_exp;
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_x_0_0_0_y_xx_y_xx, g_x_0_0_0_y_xx_y_xy, g_x_0_0_0_y_xx_y_xz, g_x_0_0_0_y_xx_y_yy, g_x_0_0_0_y_xx_y_yz, g_x_0_0_0_y_xx_y_zz, g_xy_xx_y_xx, g_xy_xx_y_xy, g_xy_xx_y_xz, g_xy_xx_y_yy, g_xy_xx_y_yz, g_xy_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xx_y_xx[i] = 2.0 * g_xy_xx_y_xx[i] * a_exp;

        g_x_0_0_0_y_xx_y_xy[i] = 2.0 * g_xy_xx_y_xy[i] * a_exp;

        g_x_0_0_0_y_xx_y_xz[i] = 2.0 * g_xy_xx_y_xz[i] * a_exp;

        g_x_0_0_0_y_xx_y_yy[i] = 2.0 * g_xy_xx_y_yy[i] * a_exp;

        g_x_0_0_0_y_xx_y_yz[i] = 2.0 * g_xy_xx_y_yz[i] * a_exp;

        g_x_0_0_0_y_xx_y_zz[i] = 2.0 * g_xy_xx_y_zz[i] * a_exp;
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_x_0_0_0_y_xx_z_xx, g_x_0_0_0_y_xx_z_xy, g_x_0_0_0_y_xx_z_xz, g_x_0_0_0_y_xx_z_yy, g_x_0_0_0_y_xx_z_yz, g_x_0_0_0_y_xx_z_zz, g_xy_xx_z_xx, g_xy_xx_z_xy, g_xy_xx_z_xz, g_xy_xx_z_yy, g_xy_xx_z_yz, g_xy_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xx_z_xx[i] = 2.0 * g_xy_xx_z_xx[i] * a_exp;

        g_x_0_0_0_y_xx_z_xy[i] = 2.0 * g_xy_xx_z_xy[i] * a_exp;

        g_x_0_0_0_y_xx_z_xz[i] = 2.0 * g_xy_xx_z_xz[i] * a_exp;

        g_x_0_0_0_y_xx_z_yy[i] = 2.0 * g_xy_xx_z_yy[i] * a_exp;

        g_x_0_0_0_y_xx_z_yz[i] = 2.0 * g_xy_xx_z_yz[i] * a_exp;

        g_x_0_0_0_y_xx_z_zz[i] = 2.0 * g_xy_xx_z_zz[i] * a_exp;
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_x_0_0_0_y_xy_x_xx, g_x_0_0_0_y_xy_x_xy, g_x_0_0_0_y_xy_x_xz, g_x_0_0_0_y_xy_x_yy, g_x_0_0_0_y_xy_x_yz, g_x_0_0_0_y_xy_x_zz, g_xy_xy_x_xx, g_xy_xy_x_xy, g_xy_xy_x_xz, g_xy_xy_x_yy, g_xy_xy_x_yz, g_xy_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xy_x_xx[i] = 2.0 * g_xy_xy_x_xx[i] * a_exp;

        g_x_0_0_0_y_xy_x_xy[i] = 2.0 * g_xy_xy_x_xy[i] * a_exp;

        g_x_0_0_0_y_xy_x_xz[i] = 2.0 * g_xy_xy_x_xz[i] * a_exp;

        g_x_0_0_0_y_xy_x_yy[i] = 2.0 * g_xy_xy_x_yy[i] * a_exp;

        g_x_0_0_0_y_xy_x_yz[i] = 2.0 * g_xy_xy_x_yz[i] * a_exp;

        g_x_0_0_0_y_xy_x_zz[i] = 2.0 * g_xy_xy_x_zz[i] * a_exp;
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_x_0_0_0_y_xy_y_xx, g_x_0_0_0_y_xy_y_xy, g_x_0_0_0_y_xy_y_xz, g_x_0_0_0_y_xy_y_yy, g_x_0_0_0_y_xy_y_yz, g_x_0_0_0_y_xy_y_zz, g_xy_xy_y_xx, g_xy_xy_y_xy, g_xy_xy_y_xz, g_xy_xy_y_yy, g_xy_xy_y_yz, g_xy_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xy_y_xx[i] = 2.0 * g_xy_xy_y_xx[i] * a_exp;

        g_x_0_0_0_y_xy_y_xy[i] = 2.0 * g_xy_xy_y_xy[i] * a_exp;

        g_x_0_0_0_y_xy_y_xz[i] = 2.0 * g_xy_xy_y_xz[i] * a_exp;

        g_x_0_0_0_y_xy_y_yy[i] = 2.0 * g_xy_xy_y_yy[i] * a_exp;

        g_x_0_0_0_y_xy_y_yz[i] = 2.0 * g_xy_xy_y_yz[i] * a_exp;

        g_x_0_0_0_y_xy_y_zz[i] = 2.0 * g_xy_xy_y_zz[i] * a_exp;
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_x_0_0_0_y_xy_z_xx, g_x_0_0_0_y_xy_z_xy, g_x_0_0_0_y_xy_z_xz, g_x_0_0_0_y_xy_z_yy, g_x_0_0_0_y_xy_z_yz, g_x_0_0_0_y_xy_z_zz, g_xy_xy_z_xx, g_xy_xy_z_xy, g_xy_xy_z_xz, g_xy_xy_z_yy, g_xy_xy_z_yz, g_xy_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xy_z_xx[i] = 2.0 * g_xy_xy_z_xx[i] * a_exp;

        g_x_0_0_0_y_xy_z_xy[i] = 2.0 * g_xy_xy_z_xy[i] * a_exp;

        g_x_0_0_0_y_xy_z_xz[i] = 2.0 * g_xy_xy_z_xz[i] * a_exp;

        g_x_0_0_0_y_xy_z_yy[i] = 2.0 * g_xy_xy_z_yy[i] * a_exp;

        g_x_0_0_0_y_xy_z_yz[i] = 2.0 * g_xy_xy_z_yz[i] * a_exp;

        g_x_0_0_0_y_xy_z_zz[i] = 2.0 * g_xy_xy_z_zz[i] * a_exp;
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_x_0_0_0_y_xz_x_xx, g_x_0_0_0_y_xz_x_xy, g_x_0_0_0_y_xz_x_xz, g_x_0_0_0_y_xz_x_yy, g_x_0_0_0_y_xz_x_yz, g_x_0_0_0_y_xz_x_zz, g_xy_xz_x_xx, g_xy_xz_x_xy, g_xy_xz_x_xz, g_xy_xz_x_yy, g_xy_xz_x_yz, g_xy_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xz_x_xx[i] = 2.0 * g_xy_xz_x_xx[i] * a_exp;

        g_x_0_0_0_y_xz_x_xy[i] = 2.0 * g_xy_xz_x_xy[i] * a_exp;

        g_x_0_0_0_y_xz_x_xz[i] = 2.0 * g_xy_xz_x_xz[i] * a_exp;

        g_x_0_0_0_y_xz_x_yy[i] = 2.0 * g_xy_xz_x_yy[i] * a_exp;

        g_x_0_0_0_y_xz_x_yz[i] = 2.0 * g_xy_xz_x_yz[i] * a_exp;

        g_x_0_0_0_y_xz_x_zz[i] = 2.0 * g_xy_xz_x_zz[i] * a_exp;
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_x_0_0_0_y_xz_y_xx, g_x_0_0_0_y_xz_y_xy, g_x_0_0_0_y_xz_y_xz, g_x_0_0_0_y_xz_y_yy, g_x_0_0_0_y_xz_y_yz, g_x_0_0_0_y_xz_y_zz, g_xy_xz_y_xx, g_xy_xz_y_xy, g_xy_xz_y_xz, g_xy_xz_y_yy, g_xy_xz_y_yz, g_xy_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xz_y_xx[i] = 2.0 * g_xy_xz_y_xx[i] * a_exp;

        g_x_0_0_0_y_xz_y_xy[i] = 2.0 * g_xy_xz_y_xy[i] * a_exp;

        g_x_0_0_0_y_xz_y_xz[i] = 2.0 * g_xy_xz_y_xz[i] * a_exp;

        g_x_0_0_0_y_xz_y_yy[i] = 2.0 * g_xy_xz_y_yy[i] * a_exp;

        g_x_0_0_0_y_xz_y_yz[i] = 2.0 * g_xy_xz_y_yz[i] * a_exp;

        g_x_0_0_0_y_xz_y_zz[i] = 2.0 * g_xy_xz_y_zz[i] * a_exp;
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_x_0_0_0_y_xz_z_xx, g_x_0_0_0_y_xz_z_xy, g_x_0_0_0_y_xz_z_xz, g_x_0_0_0_y_xz_z_yy, g_x_0_0_0_y_xz_z_yz, g_x_0_0_0_y_xz_z_zz, g_xy_xz_z_xx, g_xy_xz_z_xy, g_xy_xz_z_xz, g_xy_xz_z_yy, g_xy_xz_z_yz, g_xy_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xz_z_xx[i] = 2.0 * g_xy_xz_z_xx[i] * a_exp;

        g_x_0_0_0_y_xz_z_xy[i] = 2.0 * g_xy_xz_z_xy[i] * a_exp;

        g_x_0_0_0_y_xz_z_xz[i] = 2.0 * g_xy_xz_z_xz[i] * a_exp;

        g_x_0_0_0_y_xz_z_yy[i] = 2.0 * g_xy_xz_z_yy[i] * a_exp;

        g_x_0_0_0_y_xz_z_yz[i] = 2.0 * g_xy_xz_z_yz[i] * a_exp;

        g_x_0_0_0_y_xz_z_zz[i] = 2.0 * g_xy_xz_z_zz[i] * a_exp;
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_x_0_0_0_y_yy_x_xx, g_x_0_0_0_y_yy_x_xy, g_x_0_0_0_y_yy_x_xz, g_x_0_0_0_y_yy_x_yy, g_x_0_0_0_y_yy_x_yz, g_x_0_0_0_y_yy_x_zz, g_xy_yy_x_xx, g_xy_yy_x_xy, g_xy_yy_x_xz, g_xy_yy_x_yy, g_xy_yy_x_yz, g_xy_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_yy_x_xx[i] = 2.0 * g_xy_yy_x_xx[i] * a_exp;

        g_x_0_0_0_y_yy_x_xy[i] = 2.0 * g_xy_yy_x_xy[i] * a_exp;

        g_x_0_0_0_y_yy_x_xz[i] = 2.0 * g_xy_yy_x_xz[i] * a_exp;

        g_x_0_0_0_y_yy_x_yy[i] = 2.0 * g_xy_yy_x_yy[i] * a_exp;

        g_x_0_0_0_y_yy_x_yz[i] = 2.0 * g_xy_yy_x_yz[i] * a_exp;

        g_x_0_0_0_y_yy_x_zz[i] = 2.0 * g_xy_yy_x_zz[i] * a_exp;
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_x_0_0_0_y_yy_y_xx, g_x_0_0_0_y_yy_y_xy, g_x_0_0_0_y_yy_y_xz, g_x_0_0_0_y_yy_y_yy, g_x_0_0_0_y_yy_y_yz, g_x_0_0_0_y_yy_y_zz, g_xy_yy_y_xx, g_xy_yy_y_xy, g_xy_yy_y_xz, g_xy_yy_y_yy, g_xy_yy_y_yz, g_xy_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_yy_y_xx[i] = 2.0 * g_xy_yy_y_xx[i] * a_exp;

        g_x_0_0_0_y_yy_y_xy[i] = 2.0 * g_xy_yy_y_xy[i] * a_exp;

        g_x_0_0_0_y_yy_y_xz[i] = 2.0 * g_xy_yy_y_xz[i] * a_exp;

        g_x_0_0_0_y_yy_y_yy[i] = 2.0 * g_xy_yy_y_yy[i] * a_exp;

        g_x_0_0_0_y_yy_y_yz[i] = 2.0 * g_xy_yy_y_yz[i] * a_exp;

        g_x_0_0_0_y_yy_y_zz[i] = 2.0 * g_xy_yy_y_zz[i] * a_exp;
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_x_0_0_0_y_yy_z_xx, g_x_0_0_0_y_yy_z_xy, g_x_0_0_0_y_yy_z_xz, g_x_0_0_0_y_yy_z_yy, g_x_0_0_0_y_yy_z_yz, g_x_0_0_0_y_yy_z_zz, g_xy_yy_z_xx, g_xy_yy_z_xy, g_xy_yy_z_xz, g_xy_yy_z_yy, g_xy_yy_z_yz, g_xy_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_yy_z_xx[i] = 2.0 * g_xy_yy_z_xx[i] * a_exp;

        g_x_0_0_0_y_yy_z_xy[i] = 2.0 * g_xy_yy_z_xy[i] * a_exp;

        g_x_0_0_0_y_yy_z_xz[i] = 2.0 * g_xy_yy_z_xz[i] * a_exp;

        g_x_0_0_0_y_yy_z_yy[i] = 2.0 * g_xy_yy_z_yy[i] * a_exp;

        g_x_0_0_0_y_yy_z_yz[i] = 2.0 * g_xy_yy_z_yz[i] * a_exp;

        g_x_0_0_0_y_yy_z_zz[i] = 2.0 * g_xy_yy_z_zz[i] * a_exp;
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_x_0_0_0_y_yz_x_xx, g_x_0_0_0_y_yz_x_xy, g_x_0_0_0_y_yz_x_xz, g_x_0_0_0_y_yz_x_yy, g_x_0_0_0_y_yz_x_yz, g_x_0_0_0_y_yz_x_zz, g_xy_yz_x_xx, g_xy_yz_x_xy, g_xy_yz_x_xz, g_xy_yz_x_yy, g_xy_yz_x_yz, g_xy_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_yz_x_xx[i] = 2.0 * g_xy_yz_x_xx[i] * a_exp;

        g_x_0_0_0_y_yz_x_xy[i] = 2.0 * g_xy_yz_x_xy[i] * a_exp;

        g_x_0_0_0_y_yz_x_xz[i] = 2.0 * g_xy_yz_x_xz[i] * a_exp;

        g_x_0_0_0_y_yz_x_yy[i] = 2.0 * g_xy_yz_x_yy[i] * a_exp;

        g_x_0_0_0_y_yz_x_yz[i] = 2.0 * g_xy_yz_x_yz[i] * a_exp;

        g_x_0_0_0_y_yz_x_zz[i] = 2.0 * g_xy_yz_x_zz[i] * a_exp;
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_x_0_0_0_y_yz_y_xx, g_x_0_0_0_y_yz_y_xy, g_x_0_0_0_y_yz_y_xz, g_x_0_0_0_y_yz_y_yy, g_x_0_0_0_y_yz_y_yz, g_x_0_0_0_y_yz_y_zz, g_xy_yz_y_xx, g_xy_yz_y_xy, g_xy_yz_y_xz, g_xy_yz_y_yy, g_xy_yz_y_yz, g_xy_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_yz_y_xx[i] = 2.0 * g_xy_yz_y_xx[i] * a_exp;

        g_x_0_0_0_y_yz_y_xy[i] = 2.0 * g_xy_yz_y_xy[i] * a_exp;

        g_x_0_0_0_y_yz_y_xz[i] = 2.0 * g_xy_yz_y_xz[i] * a_exp;

        g_x_0_0_0_y_yz_y_yy[i] = 2.0 * g_xy_yz_y_yy[i] * a_exp;

        g_x_0_0_0_y_yz_y_yz[i] = 2.0 * g_xy_yz_y_yz[i] * a_exp;

        g_x_0_0_0_y_yz_y_zz[i] = 2.0 * g_xy_yz_y_zz[i] * a_exp;
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_x_0_0_0_y_yz_z_xx, g_x_0_0_0_y_yz_z_xy, g_x_0_0_0_y_yz_z_xz, g_x_0_0_0_y_yz_z_yy, g_x_0_0_0_y_yz_z_yz, g_x_0_0_0_y_yz_z_zz, g_xy_yz_z_xx, g_xy_yz_z_xy, g_xy_yz_z_xz, g_xy_yz_z_yy, g_xy_yz_z_yz, g_xy_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_yz_z_xx[i] = 2.0 * g_xy_yz_z_xx[i] * a_exp;

        g_x_0_0_0_y_yz_z_xy[i] = 2.0 * g_xy_yz_z_xy[i] * a_exp;

        g_x_0_0_0_y_yz_z_xz[i] = 2.0 * g_xy_yz_z_xz[i] * a_exp;

        g_x_0_0_0_y_yz_z_yy[i] = 2.0 * g_xy_yz_z_yy[i] * a_exp;

        g_x_0_0_0_y_yz_z_yz[i] = 2.0 * g_xy_yz_z_yz[i] * a_exp;

        g_x_0_0_0_y_yz_z_zz[i] = 2.0 * g_xy_yz_z_zz[i] * a_exp;
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_x_0_0_0_y_zz_x_xx, g_x_0_0_0_y_zz_x_xy, g_x_0_0_0_y_zz_x_xz, g_x_0_0_0_y_zz_x_yy, g_x_0_0_0_y_zz_x_yz, g_x_0_0_0_y_zz_x_zz, g_xy_zz_x_xx, g_xy_zz_x_xy, g_xy_zz_x_xz, g_xy_zz_x_yy, g_xy_zz_x_yz, g_xy_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_zz_x_xx[i] = 2.0 * g_xy_zz_x_xx[i] * a_exp;

        g_x_0_0_0_y_zz_x_xy[i] = 2.0 * g_xy_zz_x_xy[i] * a_exp;

        g_x_0_0_0_y_zz_x_xz[i] = 2.0 * g_xy_zz_x_xz[i] * a_exp;

        g_x_0_0_0_y_zz_x_yy[i] = 2.0 * g_xy_zz_x_yy[i] * a_exp;

        g_x_0_0_0_y_zz_x_yz[i] = 2.0 * g_xy_zz_x_yz[i] * a_exp;

        g_x_0_0_0_y_zz_x_zz[i] = 2.0 * g_xy_zz_x_zz[i] * a_exp;
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_x_0_0_0_y_zz_y_xx, g_x_0_0_0_y_zz_y_xy, g_x_0_0_0_y_zz_y_xz, g_x_0_0_0_y_zz_y_yy, g_x_0_0_0_y_zz_y_yz, g_x_0_0_0_y_zz_y_zz, g_xy_zz_y_xx, g_xy_zz_y_xy, g_xy_zz_y_xz, g_xy_zz_y_yy, g_xy_zz_y_yz, g_xy_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_zz_y_xx[i] = 2.0 * g_xy_zz_y_xx[i] * a_exp;

        g_x_0_0_0_y_zz_y_xy[i] = 2.0 * g_xy_zz_y_xy[i] * a_exp;

        g_x_0_0_0_y_zz_y_xz[i] = 2.0 * g_xy_zz_y_xz[i] * a_exp;

        g_x_0_0_0_y_zz_y_yy[i] = 2.0 * g_xy_zz_y_yy[i] * a_exp;

        g_x_0_0_0_y_zz_y_yz[i] = 2.0 * g_xy_zz_y_yz[i] * a_exp;

        g_x_0_0_0_y_zz_y_zz[i] = 2.0 * g_xy_zz_y_zz[i] * a_exp;
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_x_0_0_0_y_zz_z_xx, g_x_0_0_0_y_zz_z_xy, g_x_0_0_0_y_zz_z_xz, g_x_0_0_0_y_zz_z_yy, g_x_0_0_0_y_zz_z_yz, g_x_0_0_0_y_zz_z_zz, g_xy_zz_z_xx, g_xy_zz_z_xy, g_xy_zz_z_xz, g_xy_zz_z_yy, g_xy_zz_z_yz, g_xy_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_zz_z_xx[i] = 2.0 * g_xy_zz_z_xx[i] * a_exp;

        g_x_0_0_0_y_zz_z_xy[i] = 2.0 * g_xy_zz_z_xy[i] * a_exp;

        g_x_0_0_0_y_zz_z_xz[i] = 2.0 * g_xy_zz_z_xz[i] * a_exp;

        g_x_0_0_0_y_zz_z_yy[i] = 2.0 * g_xy_zz_z_yy[i] * a_exp;

        g_x_0_0_0_y_zz_z_yz[i] = 2.0 * g_xy_zz_z_yz[i] * a_exp;

        g_x_0_0_0_y_zz_z_zz[i] = 2.0 * g_xy_zz_z_zz[i] * a_exp;
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_x_0_0_0_z_xx_x_xx, g_x_0_0_0_z_xx_x_xy, g_x_0_0_0_z_xx_x_xz, g_x_0_0_0_z_xx_x_yy, g_x_0_0_0_z_xx_x_yz, g_x_0_0_0_z_xx_x_zz, g_xz_xx_x_xx, g_xz_xx_x_xy, g_xz_xx_x_xz, g_xz_xx_x_yy, g_xz_xx_x_yz, g_xz_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xx_x_xx[i] = 2.0 * g_xz_xx_x_xx[i] * a_exp;

        g_x_0_0_0_z_xx_x_xy[i] = 2.0 * g_xz_xx_x_xy[i] * a_exp;

        g_x_0_0_0_z_xx_x_xz[i] = 2.0 * g_xz_xx_x_xz[i] * a_exp;

        g_x_0_0_0_z_xx_x_yy[i] = 2.0 * g_xz_xx_x_yy[i] * a_exp;

        g_x_0_0_0_z_xx_x_yz[i] = 2.0 * g_xz_xx_x_yz[i] * a_exp;

        g_x_0_0_0_z_xx_x_zz[i] = 2.0 * g_xz_xx_x_zz[i] * a_exp;
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_x_0_0_0_z_xx_y_xx, g_x_0_0_0_z_xx_y_xy, g_x_0_0_0_z_xx_y_xz, g_x_0_0_0_z_xx_y_yy, g_x_0_0_0_z_xx_y_yz, g_x_0_0_0_z_xx_y_zz, g_xz_xx_y_xx, g_xz_xx_y_xy, g_xz_xx_y_xz, g_xz_xx_y_yy, g_xz_xx_y_yz, g_xz_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xx_y_xx[i] = 2.0 * g_xz_xx_y_xx[i] * a_exp;

        g_x_0_0_0_z_xx_y_xy[i] = 2.0 * g_xz_xx_y_xy[i] * a_exp;

        g_x_0_0_0_z_xx_y_xz[i] = 2.0 * g_xz_xx_y_xz[i] * a_exp;

        g_x_0_0_0_z_xx_y_yy[i] = 2.0 * g_xz_xx_y_yy[i] * a_exp;

        g_x_0_0_0_z_xx_y_yz[i] = 2.0 * g_xz_xx_y_yz[i] * a_exp;

        g_x_0_0_0_z_xx_y_zz[i] = 2.0 * g_xz_xx_y_zz[i] * a_exp;
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_x_0_0_0_z_xx_z_xx, g_x_0_0_0_z_xx_z_xy, g_x_0_0_0_z_xx_z_xz, g_x_0_0_0_z_xx_z_yy, g_x_0_0_0_z_xx_z_yz, g_x_0_0_0_z_xx_z_zz, g_xz_xx_z_xx, g_xz_xx_z_xy, g_xz_xx_z_xz, g_xz_xx_z_yy, g_xz_xx_z_yz, g_xz_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xx_z_xx[i] = 2.0 * g_xz_xx_z_xx[i] * a_exp;

        g_x_0_0_0_z_xx_z_xy[i] = 2.0 * g_xz_xx_z_xy[i] * a_exp;

        g_x_0_0_0_z_xx_z_xz[i] = 2.0 * g_xz_xx_z_xz[i] * a_exp;

        g_x_0_0_0_z_xx_z_yy[i] = 2.0 * g_xz_xx_z_yy[i] * a_exp;

        g_x_0_0_0_z_xx_z_yz[i] = 2.0 * g_xz_xx_z_yz[i] * a_exp;

        g_x_0_0_0_z_xx_z_zz[i] = 2.0 * g_xz_xx_z_zz[i] * a_exp;
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_x_0_0_0_z_xy_x_xx, g_x_0_0_0_z_xy_x_xy, g_x_0_0_0_z_xy_x_xz, g_x_0_0_0_z_xy_x_yy, g_x_0_0_0_z_xy_x_yz, g_x_0_0_0_z_xy_x_zz, g_xz_xy_x_xx, g_xz_xy_x_xy, g_xz_xy_x_xz, g_xz_xy_x_yy, g_xz_xy_x_yz, g_xz_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xy_x_xx[i] = 2.0 * g_xz_xy_x_xx[i] * a_exp;

        g_x_0_0_0_z_xy_x_xy[i] = 2.0 * g_xz_xy_x_xy[i] * a_exp;

        g_x_0_0_0_z_xy_x_xz[i] = 2.0 * g_xz_xy_x_xz[i] * a_exp;

        g_x_0_0_0_z_xy_x_yy[i] = 2.0 * g_xz_xy_x_yy[i] * a_exp;

        g_x_0_0_0_z_xy_x_yz[i] = 2.0 * g_xz_xy_x_yz[i] * a_exp;

        g_x_0_0_0_z_xy_x_zz[i] = 2.0 * g_xz_xy_x_zz[i] * a_exp;
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_x_0_0_0_z_xy_y_xx, g_x_0_0_0_z_xy_y_xy, g_x_0_0_0_z_xy_y_xz, g_x_0_0_0_z_xy_y_yy, g_x_0_0_0_z_xy_y_yz, g_x_0_0_0_z_xy_y_zz, g_xz_xy_y_xx, g_xz_xy_y_xy, g_xz_xy_y_xz, g_xz_xy_y_yy, g_xz_xy_y_yz, g_xz_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xy_y_xx[i] = 2.0 * g_xz_xy_y_xx[i] * a_exp;

        g_x_0_0_0_z_xy_y_xy[i] = 2.0 * g_xz_xy_y_xy[i] * a_exp;

        g_x_0_0_0_z_xy_y_xz[i] = 2.0 * g_xz_xy_y_xz[i] * a_exp;

        g_x_0_0_0_z_xy_y_yy[i] = 2.0 * g_xz_xy_y_yy[i] * a_exp;

        g_x_0_0_0_z_xy_y_yz[i] = 2.0 * g_xz_xy_y_yz[i] * a_exp;

        g_x_0_0_0_z_xy_y_zz[i] = 2.0 * g_xz_xy_y_zz[i] * a_exp;
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_x_0_0_0_z_xy_z_xx, g_x_0_0_0_z_xy_z_xy, g_x_0_0_0_z_xy_z_xz, g_x_0_0_0_z_xy_z_yy, g_x_0_0_0_z_xy_z_yz, g_x_0_0_0_z_xy_z_zz, g_xz_xy_z_xx, g_xz_xy_z_xy, g_xz_xy_z_xz, g_xz_xy_z_yy, g_xz_xy_z_yz, g_xz_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xy_z_xx[i] = 2.0 * g_xz_xy_z_xx[i] * a_exp;

        g_x_0_0_0_z_xy_z_xy[i] = 2.0 * g_xz_xy_z_xy[i] * a_exp;

        g_x_0_0_0_z_xy_z_xz[i] = 2.0 * g_xz_xy_z_xz[i] * a_exp;

        g_x_0_0_0_z_xy_z_yy[i] = 2.0 * g_xz_xy_z_yy[i] * a_exp;

        g_x_0_0_0_z_xy_z_yz[i] = 2.0 * g_xz_xy_z_yz[i] * a_exp;

        g_x_0_0_0_z_xy_z_zz[i] = 2.0 * g_xz_xy_z_zz[i] * a_exp;
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_x_0_0_0_z_xz_x_xx, g_x_0_0_0_z_xz_x_xy, g_x_0_0_0_z_xz_x_xz, g_x_0_0_0_z_xz_x_yy, g_x_0_0_0_z_xz_x_yz, g_x_0_0_0_z_xz_x_zz, g_xz_xz_x_xx, g_xz_xz_x_xy, g_xz_xz_x_xz, g_xz_xz_x_yy, g_xz_xz_x_yz, g_xz_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xz_x_xx[i] = 2.0 * g_xz_xz_x_xx[i] * a_exp;

        g_x_0_0_0_z_xz_x_xy[i] = 2.0 * g_xz_xz_x_xy[i] * a_exp;

        g_x_0_0_0_z_xz_x_xz[i] = 2.0 * g_xz_xz_x_xz[i] * a_exp;

        g_x_0_0_0_z_xz_x_yy[i] = 2.0 * g_xz_xz_x_yy[i] * a_exp;

        g_x_0_0_0_z_xz_x_yz[i] = 2.0 * g_xz_xz_x_yz[i] * a_exp;

        g_x_0_0_0_z_xz_x_zz[i] = 2.0 * g_xz_xz_x_zz[i] * a_exp;
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_x_0_0_0_z_xz_y_xx, g_x_0_0_0_z_xz_y_xy, g_x_0_0_0_z_xz_y_xz, g_x_0_0_0_z_xz_y_yy, g_x_0_0_0_z_xz_y_yz, g_x_0_0_0_z_xz_y_zz, g_xz_xz_y_xx, g_xz_xz_y_xy, g_xz_xz_y_xz, g_xz_xz_y_yy, g_xz_xz_y_yz, g_xz_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xz_y_xx[i] = 2.0 * g_xz_xz_y_xx[i] * a_exp;

        g_x_0_0_0_z_xz_y_xy[i] = 2.0 * g_xz_xz_y_xy[i] * a_exp;

        g_x_0_0_0_z_xz_y_xz[i] = 2.0 * g_xz_xz_y_xz[i] * a_exp;

        g_x_0_0_0_z_xz_y_yy[i] = 2.0 * g_xz_xz_y_yy[i] * a_exp;

        g_x_0_0_0_z_xz_y_yz[i] = 2.0 * g_xz_xz_y_yz[i] * a_exp;

        g_x_0_0_0_z_xz_y_zz[i] = 2.0 * g_xz_xz_y_zz[i] * a_exp;
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_x_0_0_0_z_xz_z_xx, g_x_0_0_0_z_xz_z_xy, g_x_0_0_0_z_xz_z_xz, g_x_0_0_0_z_xz_z_yy, g_x_0_0_0_z_xz_z_yz, g_x_0_0_0_z_xz_z_zz, g_xz_xz_z_xx, g_xz_xz_z_xy, g_xz_xz_z_xz, g_xz_xz_z_yy, g_xz_xz_z_yz, g_xz_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xz_z_xx[i] = 2.0 * g_xz_xz_z_xx[i] * a_exp;

        g_x_0_0_0_z_xz_z_xy[i] = 2.0 * g_xz_xz_z_xy[i] * a_exp;

        g_x_0_0_0_z_xz_z_xz[i] = 2.0 * g_xz_xz_z_xz[i] * a_exp;

        g_x_0_0_0_z_xz_z_yy[i] = 2.0 * g_xz_xz_z_yy[i] * a_exp;

        g_x_0_0_0_z_xz_z_yz[i] = 2.0 * g_xz_xz_z_yz[i] * a_exp;

        g_x_0_0_0_z_xz_z_zz[i] = 2.0 * g_xz_xz_z_zz[i] * a_exp;
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_x_0_0_0_z_yy_x_xx, g_x_0_0_0_z_yy_x_xy, g_x_0_0_0_z_yy_x_xz, g_x_0_0_0_z_yy_x_yy, g_x_0_0_0_z_yy_x_yz, g_x_0_0_0_z_yy_x_zz, g_xz_yy_x_xx, g_xz_yy_x_xy, g_xz_yy_x_xz, g_xz_yy_x_yy, g_xz_yy_x_yz, g_xz_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_yy_x_xx[i] = 2.0 * g_xz_yy_x_xx[i] * a_exp;

        g_x_0_0_0_z_yy_x_xy[i] = 2.0 * g_xz_yy_x_xy[i] * a_exp;

        g_x_0_0_0_z_yy_x_xz[i] = 2.0 * g_xz_yy_x_xz[i] * a_exp;

        g_x_0_0_0_z_yy_x_yy[i] = 2.0 * g_xz_yy_x_yy[i] * a_exp;

        g_x_0_0_0_z_yy_x_yz[i] = 2.0 * g_xz_yy_x_yz[i] * a_exp;

        g_x_0_0_0_z_yy_x_zz[i] = 2.0 * g_xz_yy_x_zz[i] * a_exp;
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_x_0_0_0_z_yy_y_xx, g_x_0_0_0_z_yy_y_xy, g_x_0_0_0_z_yy_y_xz, g_x_0_0_0_z_yy_y_yy, g_x_0_0_0_z_yy_y_yz, g_x_0_0_0_z_yy_y_zz, g_xz_yy_y_xx, g_xz_yy_y_xy, g_xz_yy_y_xz, g_xz_yy_y_yy, g_xz_yy_y_yz, g_xz_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_yy_y_xx[i] = 2.0 * g_xz_yy_y_xx[i] * a_exp;

        g_x_0_0_0_z_yy_y_xy[i] = 2.0 * g_xz_yy_y_xy[i] * a_exp;

        g_x_0_0_0_z_yy_y_xz[i] = 2.0 * g_xz_yy_y_xz[i] * a_exp;

        g_x_0_0_0_z_yy_y_yy[i] = 2.0 * g_xz_yy_y_yy[i] * a_exp;

        g_x_0_0_0_z_yy_y_yz[i] = 2.0 * g_xz_yy_y_yz[i] * a_exp;

        g_x_0_0_0_z_yy_y_zz[i] = 2.0 * g_xz_yy_y_zz[i] * a_exp;
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_x_0_0_0_z_yy_z_xx, g_x_0_0_0_z_yy_z_xy, g_x_0_0_0_z_yy_z_xz, g_x_0_0_0_z_yy_z_yy, g_x_0_0_0_z_yy_z_yz, g_x_0_0_0_z_yy_z_zz, g_xz_yy_z_xx, g_xz_yy_z_xy, g_xz_yy_z_xz, g_xz_yy_z_yy, g_xz_yy_z_yz, g_xz_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_yy_z_xx[i] = 2.0 * g_xz_yy_z_xx[i] * a_exp;

        g_x_0_0_0_z_yy_z_xy[i] = 2.0 * g_xz_yy_z_xy[i] * a_exp;

        g_x_0_0_0_z_yy_z_xz[i] = 2.0 * g_xz_yy_z_xz[i] * a_exp;

        g_x_0_0_0_z_yy_z_yy[i] = 2.0 * g_xz_yy_z_yy[i] * a_exp;

        g_x_0_0_0_z_yy_z_yz[i] = 2.0 * g_xz_yy_z_yz[i] * a_exp;

        g_x_0_0_0_z_yy_z_zz[i] = 2.0 * g_xz_yy_z_zz[i] * a_exp;
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_x_0_0_0_z_yz_x_xx, g_x_0_0_0_z_yz_x_xy, g_x_0_0_0_z_yz_x_xz, g_x_0_0_0_z_yz_x_yy, g_x_0_0_0_z_yz_x_yz, g_x_0_0_0_z_yz_x_zz, g_xz_yz_x_xx, g_xz_yz_x_xy, g_xz_yz_x_xz, g_xz_yz_x_yy, g_xz_yz_x_yz, g_xz_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_yz_x_xx[i] = 2.0 * g_xz_yz_x_xx[i] * a_exp;

        g_x_0_0_0_z_yz_x_xy[i] = 2.0 * g_xz_yz_x_xy[i] * a_exp;

        g_x_0_0_0_z_yz_x_xz[i] = 2.0 * g_xz_yz_x_xz[i] * a_exp;

        g_x_0_0_0_z_yz_x_yy[i] = 2.0 * g_xz_yz_x_yy[i] * a_exp;

        g_x_0_0_0_z_yz_x_yz[i] = 2.0 * g_xz_yz_x_yz[i] * a_exp;

        g_x_0_0_0_z_yz_x_zz[i] = 2.0 * g_xz_yz_x_zz[i] * a_exp;
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_x_0_0_0_z_yz_y_xx, g_x_0_0_0_z_yz_y_xy, g_x_0_0_0_z_yz_y_xz, g_x_0_0_0_z_yz_y_yy, g_x_0_0_0_z_yz_y_yz, g_x_0_0_0_z_yz_y_zz, g_xz_yz_y_xx, g_xz_yz_y_xy, g_xz_yz_y_xz, g_xz_yz_y_yy, g_xz_yz_y_yz, g_xz_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_yz_y_xx[i] = 2.0 * g_xz_yz_y_xx[i] * a_exp;

        g_x_0_0_0_z_yz_y_xy[i] = 2.0 * g_xz_yz_y_xy[i] * a_exp;

        g_x_0_0_0_z_yz_y_xz[i] = 2.0 * g_xz_yz_y_xz[i] * a_exp;

        g_x_0_0_0_z_yz_y_yy[i] = 2.0 * g_xz_yz_y_yy[i] * a_exp;

        g_x_0_0_0_z_yz_y_yz[i] = 2.0 * g_xz_yz_y_yz[i] * a_exp;

        g_x_0_0_0_z_yz_y_zz[i] = 2.0 * g_xz_yz_y_zz[i] * a_exp;
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_x_0_0_0_z_yz_z_xx, g_x_0_0_0_z_yz_z_xy, g_x_0_0_0_z_yz_z_xz, g_x_0_0_0_z_yz_z_yy, g_x_0_0_0_z_yz_z_yz, g_x_0_0_0_z_yz_z_zz, g_xz_yz_z_xx, g_xz_yz_z_xy, g_xz_yz_z_xz, g_xz_yz_z_yy, g_xz_yz_z_yz, g_xz_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_yz_z_xx[i] = 2.0 * g_xz_yz_z_xx[i] * a_exp;

        g_x_0_0_0_z_yz_z_xy[i] = 2.0 * g_xz_yz_z_xy[i] * a_exp;

        g_x_0_0_0_z_yz_z_xz[i] = 2.0 * g_xz_yz_z_xz[i] * a_exp;

        g_x_0_0_0_z_yz_z_yy[i] = 2.0 * g_xz_yz_z_yy[i] * a_exp;

        g_x_0_0_0_z_yz_z_yz[i] = 2.0 * g_xz_yz_z_yz[i] * a_exp;

        g_x_0_0_0_z_yz_z_zz[i] = 2.0 * g_xz_yz_z_zz[i] * a_exp;
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_x_0_0_0_z_zz_x_xx, g_x_0_0_0_z_zz_x_xy, g_x_0_0_0_z_zz_x_xz, g_x_0_0_0_z_zz_x_yy, g_x_0_0_0_z_zz_x_yz, g_x_0_0_0_z_zz_x_zz, g_xz_zz_x_xx, g_xz_zz_x_xy, g_xz_zz_x_xz, g_xz_zz_x_yy, g_xz_zz_x_yz, g_xz_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_zz_x_xx[i] = 2.0 * g_xz_zz_x_xx[i] * a_exp;

        g_x_0_0_0_z_zz_x_xy[i] = 2.0 * g_xz_zz_x_xy[i] * a_exp;

        g_x_0_0_0_z_zz_x_xz[i] = 2.0 * g_xz_zz_x_xz[i] * a_exp;

        g_x_0_0_0_z_zz_x_yy[i] = 2.0 * g_xz_zz_x_yy[i] * a_exp;

        g_x_0_0_0_z_zz_x_yz[i] = 2.0 * g_xz_zz_x_yz[i] * a_exp;

        g_x_0_0_0_z_zz_x_zz[i] = 2.0 * g_xz_zz_x_zz[i] * a_exp;
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_x_0_0_0_z_zz_y_xx, g_x_0_0_0_z_zz_y_xy, g_x_0_0_0_z_zz_y_xz, g_x_0_0_0_z_zz_y_yy, g_x_0_0_0_z_zz_y_yz, g_x_0_0_0_z_zz_y_zz, g_xz_zz_y_xx, g_xz_zz_y_xy, g_xz_zz_y_xz, g_xz_zz_y_yy, g_xz_zz_y_yz, g_xz_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_zz_y_xx[i] = 2.0 * g_xz_zz_y_xx[i] * a_exp;

        g_x_0_0_0_z_zz_y_xy[i] = 2.0 * g_xz_zz_y_xy[i] * a_exp;

        g_x_0_0_0_z_zz_y_xz[i] = 2.0 * g_xz_zz_y_xz[i] * a_exp;

        g_x_0_0_0_z_zz_y_yy[i] = 2.0 * g_xz_zz_y_yy[i] * a_exp;

        g_x_0_0_0_z_zz_y_yz[i] = 2.0 * g_xz_zz_y_yz[i] * a_exp;

        g_x_0_0_0_z_zz_y_zz[i] = 2.0 * g_xz_zz_y_zz[i] * a_exp;
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_x_0_0_0_z_zz_z_xx, g_x_0_0_0_z_zz_z_xy, g_x_0_0_0_z_zz_z_xz, g_x_0_0_0_z_zz_z_yy, g_x_0_0_0_z_zz_z_yz, g_x_0_0_0_z_zz_z_zz, g_xz_zz_z_xx, g_xz_zz_z_xy, g_xz_zz_z_xz, g_xz_zz_z_yy, g_xz_zz_z_yz, g_xz_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_zz_z_xx[i] = 2.0 * g_xz_zz_z_xx[i] * a_exp;

        g_x_0_0_0_z_zz_z_xy[i] = 2.0 * g_xz_zz_z_xy[i] * a_exp;

        g_x_0_0_0_z_zz_z_xz[i] = 2.0 * g_xz_zz_z_xz[i] * a_exp;

        g_x_0_0_0_z_zz_z_yy[i] = 2.0 * g_xz_zz_z_yy[i] * a_exp;

        g_x_0_0_0_z_zz_z_yz[i] = 2.0 * g_xz_zz_z_yz[i] * a_exp;

        g_x_0_0_0_z_zz_z_zz[i] = 2.0 * g_xz_zz_z_zz[i] * a_exp;
    }
    // integrals block (324-330)

    #pragma omp simd aligned(g_xy_xx_x_xx, g_xy_xx_x_xy, g_xy_xx_x_xz, g_xy_xx_x_yy, g_xy_xx_x_yz, g_xy_xx_x_zz, g_y_0_0_0_x_xx_x_xx, g_y_0_0_0_x_xx_x_xy, g_y_0_0_0_x_xx_x_xz, g_y_0_0_0_x_xx_x_yy, g_y_0_0_0_x_xx_x_yz, g_y_0_0_0_x_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xx_x_xx[i] = 2.0 * g_xy_xx_x_xx[i] * a_exp;

        g_y_0_0_0_x_xx_x_xy[i] = 2.0 * g_xy_xx_x_xy[i] * a_exp;

        g_y_0_0_0_x_xx_x_xz[i] = 2.0 * g_xy_xx_x_xz[i] * a_exp;

        g_y_0_0_0_x_xx_x_yy[i] = 2.0 * g_xy_xx_x_yy[i] * a_exp;

        g_y_0_0_0_x_xx_x_yz[i] = 2.0 * g_xy_xx_x_yz[i] * a_exp;

        g_y_0_0_0_x_xx_x_zz[i] = 2.0 * g_xy_xx_x_zz[i] * a_exp;
    }
    // integrals block (330-336)

    #pragma omp simd aligned(g_xy_xx_y_xx, g_xy_xx_y_xy, g_xy_xx_y_xz, g_xy_xx_y_yy, g_xy_xx_y_yz, g_xy_xx_y_zz, g_y_0_0_0_x_xx_y_xx, g_y_0_0_0_x_xx_y_xy, g_y_0_0_0_x_xx_y_xz, g_y_0_0_0_x_xx_y_yy, g_y_0_0_0_x_xx_y_yz, g_y_0_0_0_x_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xx_y_xx[i] = 2.0 * g_xy_xx_y_xx[i] * a_exp;

        g_y_0_0_0_x_xx_y_xy[i] = 2.0 * g_xy_xx_y_xy[i] * a_exp;

        g_y_0_0_0_x_xx_y_xz[i] = 2.0 * g_xy_xx_y_xz[i] * a_exp;

        g_y_0_0_0_x_xx_y_yy[i] = 2.0 * g_xy_xx_y_yy[i] * a_exp;

        g_y_0_0_0_x_xx_y_yz[i] = 2.0 * g_xy_xx_y_yz[i] * a_exp;

        g_y_0_0_0_x_xx_y_zz[i] = 2.0 * g_xy_xx_y_zz[i] * a_exp;
    }
    // integrals block (336-342)

    #pragma omp simd aligned(g_xy_xx_z_xx, g_xy_xx_z_xy, g_xy_xx_z_xz, g_xy_xx_z_yy, g_xy_xx_z_yz, g_xy_xx_z_zz, g_y_0_0_0_x_xx_z_xx, g_y_0_0_0_x_xx_z_xy, g_y_0_0_0_x_xx_z_xz, g_y_0_0_0_x_xx_z_yy, g_y_0_0_0_x_xx_z_yz, g_y_0_0_0_x_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xx_z_xx[i] = 2.0 * g_xy_xx_z_xx[i] * a_exp;

        g_y_0_0_0_x_xx_z_xy[i] = 2.0 * g_xy_xx_z_xy[i] * a_exp;

        g_y_0_0_0_x_xx_z_xz[i] = 2.0 * g_xy_xx_z_xz[i] * a_exp;

        g_y_0_0_0_x_xx_z_yy[i] = 2.0 * g_xy_xx_z_yy[i] * a_exp;

        g_y_0_0_0_x_xx_z_yz[i] = 2.0 * g_xy_xx_z_yz[i] * a_exp;

        g_y_0_0_0_x_xx_z_zz[i] = 2.0 * g_xy_xx_z_zz[i] * a_exp;
    }
    // integrals block (342-348)

    #pragma omp simd aligned(g_xy_xy_x_xx, g_xy_xy_x_xy, g_xy_xy_x_xz, g_xy_xy_x_yy, g_xy_xy_x_yz, g_xy_xy_x_zz, g_y_0_0_0_x_xy_x_xx, g_y_0_0_0_x_xy_x_xy, g_y_0_0_0_x_xy_x_xz, g_y_0_0_0_x_xy_x_yy, g_y_0_0_0_x_xy_x_yz, g_y_0_0_0_x_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xy_x_xx[i] = 2.0 * g_xy_xy_x_xx[i] * a_exp;

        g_y_0_0_0_x_xy_x_xy[i] = 2.0 * g_xy_xy_x_xy[i] * a_exp;

        g_y_0_0_0_x_xy_x_xz[i] = 2.0 * g_xy_xy_x_xz[i] * a_exp;

        g_y_0_0_0_x_xy_x_yy[i] = 2.0 * g_xy_xy_x_yy[i] * a_exp;

        g_y_0_0_0_x_xy_x_yz[i] = 2.0 * g_xy_xy_x_yz[i] * a_exp;

        g_y_0_0_0_x_xy_x_zz[i] = 2.0 * g_xy_xy_x_zz[i] * a_exp;
    }
    // integrals block (348-354)

    #pragma omp simd aligned(g_xy_xy_y_xx, g_xy_xy_y_xy, g_xy_xy_y_xz, g_xy_xy_y_yy, g_xy_xy_y_yz, g_xy_xy_y_zz, g_y_0_0_0_x_xy_y_xx, g_y_0_0_0_x_xy_y_xy, g_y_0_0_0_x_xy_y_xz, g_y_0_0_0_x_xy_y_yy, g_y_0_0_0_x_xy_y_yz, g_y_0_0_0_x_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xy_y_xx[i] = 2.0 * g_xy_xy_y_xx[i] * a_exp;

        g_y_0_0_0_x_xy_y_xy[i] = 2.0 * g_xy_xy_y_xy[i] * a_exp;

        g_y_0_0_0_x_xy_y_xz[i] = 2.0 * g_xy_xy_y_xz[i] * a_exp;

        g_y_0_0_0_x_xy_y_yy[i] = 2.0 * g_xy_xy_y_yy[i] * a_exp;

        g_y_0_0_0_x_xy_y_yz[i] = 2.0 * g_xy_xy_y_yz[i] * a_exp;

        g_y_0_0_0_x_xy_y_zz[i] = 2.0 * g_xy_xy_y_zz[i] * a_exp;
    }
    // integrals block (354-360)

    #pragma omp simd aligned(g_xy_xy_z_xx, g_xy_xy_z_xy, g_xy_xy_z_xz, g_xy_xy_z_yy, g_xy_xy_z_yz, g_xy_xy_z_zz, g_y_0_0_0_x_xy_z_xx, g_y_0_0_0_x_xy_z_xy, g_y_0_0_0_x_xy_z_xz, g_y_0_0_0_x_xy_z_yy, g_y_0_0_0_x_xy_z_yz, g_y_0_0_0_x_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xy_z_xx[i] = 2.0 * g_xy_xy_z_xx[i] * a_exp;

        g_y_0_0_0_x_xy_z_xy[i] = 2.0 * g_xy_xy_z_xy[i] * a_exp;

        g_y_0_0_0_x_xy_z_xz[i] = 2.0 * g_xy_xy_z_xz[i] * a_exp;

        g_y_0_0_0_x_xy_z_yy[i] = 2.0 * g_xy_xy_z_yy[i] * a_exp;

        g_y_0_0_0_x_xy_z_yz[i] = 2.0 * g_xy_xy_z_yz[i] * a_exp;

        g_y_0_0_0_x_xy_z_zz[i] = 2.0 * g_xy_xy_z_zz[i] * a_exp;
    }
    // integrals block (360-366)

    #pragma omp simd aligned(g_xy_xz_x_xx, g_xy_xz_x_xy, g_xy_xz_x_xz, g_xy_xz_x_yy, g_xy_xz_x_yz, g_xy_xz_x_zz, g_y_0_0_0_x_xz_x_xx, g_y_0_0_0_x_xz_x_xy, g_y_0_0_0_x_xz_x_xz, g_y_0_0_0_x_xz_x_yy, g_y_0_0_0_x_xz_x_yz, g_y_0_0_0_x_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xz_x_xx[i] = 2.0 * g_xy_xz_x_xx[i] * a_exp;

        g_y_0_0_0_x_xz_x_xy[i] = 2.0 * g_xy_xz_x_xy[i] * a_exp;

        g_y_0_0_0_x_xz_x_xz[i] = 2.0 * g_xy_xz_x_xz[i] * a_exp;

        g_y_0_0_0_x_xz_x_yy[i] = 2.0 * g_xy_xz_x_yy[i] * a_exp;

        g_y_0_0_0_x_xz_x_yz[i] = 2.0 * g_xy_xz_x_yz[i] * a_exp;

        g_y_0_0_0_x_xz_x_zz[i] = 2.0 * g_xy_xz_x_zz[i] * a_exp;
    }
    // integrals block (366-372)

    #pragma omp simd aligned(g_xy_xz_y_xx, g_xy_xz_y_xy, g_xy_xz_y_xz, g_xy_xz_y_yy, g_xy_xz_y_yz, g_xy_xz_y_zz, g_y_0_0_0_x_xz_y_xx, g_y_0_0_0_x_xz_y_xy, g_y_0_0_0_x_xz_y_xz, g_y_0_0_0_x_xz_y_yy, g_y_0_0_0_x_xz_y_yz, g_y_0_0_0_x_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xz_y_xx[i] = 2.0 * g_xy_xz_y_xx[i] * a_exp;

        g_y_0_0_0_x_xz_y_xy[i] = 2.0 * g_xy_xz_y_xy[i] * a_exp;

        g_y_0_0_0_x_xz_y_xz[i] = 2.0 * g_xy_xz_y_xz[i] * a_exp;

        g_y_0_0_0_x_xz_y_yy[i] = 2.0 * g_xy_xz_y_yy[i] * a_exp;

        g_y_0_0_0_x_xz_y_yz[i] = 2.0 * g_xy_xz_y_yz[i] * a_exp;

        g_y_0_0_0_x_xz_y_zz[i] = 2.0 * g_xy_xz_y_zz[i] * a_exp;
    }
    // integrals block (372-378)

    #pragma omp simd aligned(g_xy_xz_z_xx, g_xy_xz_z_xy, g_xy_xz_z_xz, g_xy_xz_z_yy, g_xy_xz_z_yz, g_xy_xz_z_zz, g_y_0_0_0_x_xz_z_xx, g_y_0_0_0_x_xz_z_xy, g_y_0_0_0_x_xz_z_xz, g_y_0_0_0_x_xz_z_yy, g_y_0_0_0_x_xz_z_yz, g_y_0_0_0_x_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xz_z_xx[i] = 2.0 * g_xy_xz_z_xx[i] * a_exp;

        g_y_0_0_0_x_xz_z_xy[i] = 2.0 * g_xy_xz_z_xy[i] * a_exp;

        g_y_0_0_0_x_xz_z_xz[i] = 2.0 * g_xy_xz_z_xz[i] * a_exp;

        g_y_0_0_0_x_xz_z_yy[i] = 2.0 * g_xy_xz_z_yy[i] * a_exp;

        g_y_0_0_0_x_xz_z_yz[i] = 2.0 * g_xy_xz_z_yz[i] * a_exp;

        g_y_0_0_0_x_xz_z_zz[i] = 2.0 * g_xy_xz_z_zz[i] * a_exp;
    }
    // integrals block (378-384)

    #pragma omp simd aligned(g_xy_yy_x_xx, g_xy_yy_x_xy, g_xy_yy_x_xz, g_xy_yy_x_yy, g_xy_yy_x_yz, g_xy_yy_x_zz, g_y_0_0_0_x_yy_x_xx, g_y_0_0_0_x_yy_x_xy, g_y_0_0_0_x_yy_x_xz, g_y_0_0_0_x_yy_x_yy, g_y_0_0_0_x_yy_x_yz, g_y_0_0_0_x_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_yy_x_xx[i] = 2.0 * g_xy_yy_x_xx[i] * a_exp;

        g_y_0_0_0_x_yy_x_xy[i] = 2.0 * g_xy_yy_x_xy[i] * a_exp;

        g_y_0_0_0_x_yy_x_xz[i] = 2.0 * g_xy_yy_x_xz[i] * a_exp;

        g_y_0_0_0_x_yy_x_yy[i] = 2.0 * g_xy_yy_x_yy[i] * a_exp;

        g_y_0_0_0_x_yy_x_yz[i] = 2.0 * g_xy_yy_x_yz[i] * a_exp;

        g_y_0_0_0_x_yy_x_zz[i] = 2.0 * g_xy_yy_x_zz[i] * a_exp;
    }
    // integrals block (384-390)

    #pragma omp simd aligned(g_xy_yy_y_xx, g_xy_yy_y_xy, g_xy_yy_y_xz, g_xy_yy_y_yy, g_xy_yy_y_yz, g_xy_yy_y_zz, g_y_0_0_0_x_yy_y_xx, g_y_0_0_0_x_yy_y_xy, g_y_0_0_0_x_yy_y_xz, g_y_0_0_0_x_yy_y_yy, g_y_0_0_0_x_yy_y_yz, g_y_0_0_0_x_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_yy_y_xx[i] = 2.0 * g_xy_yy_y_xx[i] * a_exp;

        g_y_0_0_0_x_yy_y_xy[i] = 2.0 * g_xy_yy_y_xy[i] * a_exp;

        g_y_0_0_0_x_yy_y_xz[i] = 2.0 * g_xy_yy_y_xz[i] * a_exp;

        g_y_0_0_0_x_yy_y_yy[i] = 2.0 * g_xy_yy_y_yy[i] * a_exp;

        g_y_0_0_0_x_yy_y_yz[i] = 2.0 * g_xy_yy_y_yz[i] * a_exp;

        g_y_0_0_0_x_yy_y_zz[i] = 2.0 * g_xy_yy_y_zz[i] * a_exp;
    }
    // integrals block (390-396)

    #pragma omp simd aligned(g_xy_yy_z_xx, g_xy_yy_z_xy, g_xy_yy_z_xz, g_xy_yy_z_yy, g_xy_yy_z_yz, g_xy_yy_z_zz, g_y_0_0_0_x_yy_z_xx, g_y_0_0_0_x_yy_z_xy, g_y_0_0_0_x_yy_z_xz, g_y_0_0_0_x_yy_z_yy, g_y_0_0_0_x_yy_z_yz, g_y_0_0_0_x_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_yy_z_xx[i] = 2.0 * g_xy_yy_z_xx[i] * a_exp;

        g_y_0_0_0_x_yy_z_xy[i] = 2.0 * g_xy_yy_z_xy[i] * a_exp;

        g_y_0_0_0_x_yy_z_xz[i] = 2.0 * g_xy_yy_z_xz[i] * a_exp;

        g_y_0_0_0_x_yy_z_yy[i] = 2.0 * g_xy_yy_z_yy[i] * a_exp;

        g_y_0_0_0_x_yy_z_yz[i] = 2.0 * g_xy_yy_z_yz[i] * a_exp;

        g_y_0_0_0_x_yy_z_zz[i] = 2.0 * g_xy_yy_z_zz[i] * a_exp;
    }
    // integrals block (396-402)

    #pragma omp simd aligned(g_xy_yz_x_xx, g_xy_yz_x_xy, g_xy_yz_x_xz, g_xy_yz_x_yy, g_xy_yz_x_yz, g_xy_yz_x_zz, g_y_0_0_0_x_yz_x_xx, g_y_0_0_0_x_yz_x_xy, g_y_0_0_0_x_yz_x_xz, g_y_0_0_0_x_yz_x_yy, g_y_0_0_0_x_yz_x_yz, g_y_0_0_0_x_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_yz_x_xx[i] = 2.0 * g_xy_yz_x_xx[i] * a_exp;

        g_y_0_0_0_x_yz_x_xy[i] = 2.0 * g_xy_yz_x_xy[i] * a_exp;

        g_y_0_0_0_x_yz_x_xz[i] = 2.0 * g_xy_yz_x_xz[i] * a_exp;

        g_y_0_0_0_x_yz_x_yy[i] = 2.0 * g_xy_yz_x_yy[i] * a_exp;

        g_y_0_0_0_x_yz_x_yz[i] = 2.0 * g_xy_yz_x_yz[i] * a_exp;

        g_y_0_0_0_x_yz_x_zz[i] = 2.0 * g_xy_yz_x_zz[i] * a_exp;
    }
    // integrals block (402-408)

    #pragma omp simd aligned(g_xy_yz_y_xx, g_xy_yz_y_xy, g_xy_yz_y_xz, g_xy_yz_y_yy, g_xy_yz_y_yz, g_xy_yz_y_zz, g_y_0_0_0_x_yz_y_xx, g_y_0_0_0_x_yz_y_xy, g_y_0_0_0_x_yz_y_xz, g_y_0_0_0_x_yz_y_yy, g_y_0_0_0_x_yz_y_yz, g_y_0_0_0_x_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_yz_y_xx[i] = 2.0 * g_xy_yz_y_xx[i] * a_exp;

        g_y_0_0_0_x_yz_y_xy[i] = 2.0 * g_xy_yz_y_xy[i] * a_exp;

        g_y_0_0_0_x_yz_y_xz[i] = 2.0 * g_xy_yz_y_xz[i] * a_exp;

        g_y_0_0_0_x_yz_y_yy[i] = 2.0 * g_xy_yz_y_yy[i] * a_exp;

        g_y_0_0_0_x_yz_y_yz[i] = 2.0 * g_xy_yz_y_yz[i] * a_exp;

        g_y_0_0_0_x_yz_y_zz[i] = 2.0 * g_xy_yz_y_zz[i] * a_exp;
    }
    // integrals block (408-414)

    #pragma omp simd aligned(g_xy_yz_z_xx, g_xy_yz_z_xy, g_xy_yz_z_xz, g_xy_yz_z_yy, g_xy_yz_z_yz, g_xy_yz_z_zz, g_y_0_0_0_x_yz_z_xx, g_y_0_0_0_x_yz_z_xy, g_y_0_0_0_x_yz_z_xz, g_y_0_0_0_x_yz_z_yy, g_y_0_0_0_x_yz_z_yz, g_y_0_0_0_x_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_yz_z_xx[i] = 2.0 * g_xy_yz_z_xx[i] * a_exp;

        g_y_0_0_0_x_yz_z_xy[i] = 2.0 * g_xy_yz_z_xy[i] * a_exp;

        g_y_0_0_0_x_yz_z_xz[i] = 2.0 * g_xy_yz_z_xz[i] * a_exp;

        g_y_0_0_0_x_yz_z_yy[i] = 2.0 * g_xy_yz_z_yy[i] * a_exp;

        g_y_0_0_0_x_yz_z_yz[i] = 2.0 * g_xy_yz_z_yz[i] * a_exp;

        g_y_0_0_0_x_yz_z_zz[i] = 2.0 * g_xy_yz_z_zz[i] * a_exp;
    }
    // integrals block (414-420)

    #pragma omp simd aligned(g_xy_zz_x_xx, g_xy_zz_x_xy, g_xy_zz_x_xz, g_xy_zz_x_yy, g_xy_zz_x_yz, g_xy_zz_x_zz, g_y_0_0_0_x_zz_x_xx, g_y_0_0_0_x_zz_x_xy, g_y_0_0_0_x_zz_x_xz, g_y_0_0_0_x_zz_x_yy, g_y_0_0_0_x_zz_x_yz, g_y_0_0_0_x_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_zz_x_xx[i] = 2.0 * g_xy_zz_x_xx[i] * a_exp;

        g_y_0_0_0_x_zz_x_xy[i] = 2.0 * g_xy_zz_x_xy[i] * a_exp;

        g_y_0_0_0_x_zz_x_xz[i] = 2.0 * g_xy_zz_x_xz[i] * a_exp;

        g_y_0_0_0_x_zz_x_yy[i] = 2.0 * g_xy_zz_x_yy[i] * a_exp;

        g_y_0_0_0_x_zz_x_yz[i] = 2.0 * g_xy_zz_x_yz[i] * a_exp;

        g_y_0_0_0_x_zz_x_zz[i] = 2.0 * g_xy_zz_x_zz[i] * a_exp;
    }
    // integrals block (420-426)

    #pragma omp simd aligned(g_xy_zz_y_xx, g_xy_zz_y_xy, g_xy_zz_y_xz, g_xy_zz_y_yy, g_xy_zz_y_yz, g_xy_zz_y_zz, g_y_0_0_0_x_zz_y_xx, g_y_0_0_0_x_zz_y_xy, g_y_0_0_0_x_zz_y_xz, g_y_0_0_0_x_zz_y_yy, g_y_0_0_0_x_zz_y_yz, g_y_0_0_0_x_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_zz_y_xx[i] = 2.0 * g_xy_zz_y_xx[i] * a_exp;

        g_y_0_0_0_x_zz_y_xy[i] = 2.0 * g_xy_zz_y_xy[i] * a_exp;

        g_y_0_0_0_x_zz_y_xz[i] = 2.0 * g_xy_zz_y_xz[i] * a_exp;

        g_y_0_0_0_x_zz_y_yy[i] = 2.0 * g_xy_zz_y_yy[i] * a_exp;

        g_y_0_0_0_x_zz_y_yz[i] = 2.0 * g_xy_zz_y_yz[i] * a_exp;

        g_y_0_0_0_x_zz_y_zz[i] = 2.0 * g_xy_zz_y_zz[i] * a_exp;
    }
    // integrals block (426-432)

    #pragma omp simd aligned(g_xy_zz_z_xx, g_xy_zz_z_xy, g_xy_zz_z_xz, g_xy_zz_z_yy, g_xy_zz_z_yz, g_xy_zz_z_zz, g_y_0_0_0_x_zz_z_xx, g_y_0_0_0_x_zz_z_xy, g_y_0_0_0_x_zz_z_xz, g_y_0_0_0_x_zz_z_yy, g_y_0_0_0_x_zz_z_yz, g_y_0_0_0_x_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_zz_z_xx[i] = 2.0 * g_xy_zz_z_xx[i] * a_exp;

        g_y_0_0_0_x_zz_z_xy[i] = 2.0 * g_xy_zz_z_xy[i] * a_exp;

        g_y_0_0_0_x_zz_z_xz[i] = 2.0 * g_xy_zz_z_xz[i] * a_exp;

        g_y_0_0_0_x_zz_z_yy[i] = 2.0 * g_xy_zz_z_yy[i] * a_exp;

        g_y_0_0_0_x_zz_z_yz[i] = 2.0 * g_xy_zz_z_yz[i] * a_exp;

        g_y_0_0_0_x_zz_z_zz[i] = 2.0 * g_xy_zz_z_zz[i] * a_exp;
    }
    // integrals block (432-438)

    #pragma omp simd aligned(g_0_xx_x_xx, g_0_xx_x_xy, g_0_xx_x_xz, g_0_xx_x_yy, g_0_xx_x_yz, g_0_xx_x_zz, g_y_0_0_0_y_xx_x_xx, g_y_0_0_0_y_xx_x_xy, g_y_0_0_0_y_xx_x_xz, g_y_0_0_0_y_xx_x_yy, g_y_0_0_0_y_xx_x_yz, g_y_0_0_0_y_xx_x_zz, g_yy_xx_x_xx, g_yy_xx_x_xy, g_yy_xx_x_xz, g_yy_xx_x_yy, g_yy_xx_x_yz, g_yy_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xx_x_xx[i] = -g_0_xx_x_xx[i] + 2.0 * g_yy_xx_x_xx[i] * a_exp;

        g_y_0_0_0_y_xx_x_xy[i] = -g_0_xx_x_xy[i] + 2.0 * g_yy_xx_x_xy[i] * a_exp;

        g_y_0_0_0_y_xx_x_xz[i] = -g_0_xx_x_xz[i] + 2.0 * g_yy_xx_x_xz[i] * a_exp;

        g_y_0_0_0_y_xx_x_yy[i] = -g_0_xx_x_yy[i] + 2.0 * g_yy_xx_x_yy[i] * a_exp;

        g_y_0_0_0_y_xx_x_yz[i] = -g_0_xx_x_yz[i] + 2.0 * g_yy_xx_x_yz[i] * a_exp;

        g_y_0_0_0_y_xx_x_zz[i] = -g_0_xx_x_zz[i] + 2.0 * g_yy_xx_x_zz[i] * a_exp;
    }
    // integrals block (438-444)

    #pragma omp simd aligned(g_0_xx_y_xx, g_0_xx_y_xy, g_0_xx_y_xz, g_0_xx_y_yy, g_0_xx_y_yz, g_0_xx_y_zz, g_y_0_0_0_y_xx_y_xx, g_y_0_0_0_y_xx_y_xy, g_y_0_0_0_y_xx_y_xz, g_y_0_0_0_y_xx_y_yy, g_y_0_0_0_y_xx_y_yz, g_y_0_0_0_y_xx_y_zz, g_yy_xx_y_xx, g_yy_xx_y_xy, g_yy_xx_y_xz, g_yy_xx_y_yy, g_yy_xx_y_yz, g_yy_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xx_y_xx[i] = -g_0_xx_y_xx[i] + 2.0 * g_yy_xx_y_xx[i] * a_exp;

        g_y_0_0_0_y_xx_y_xy[i] = -g_0_xx_y_xy[i] + 2.0 * g_yy_xx_y_xy[i] * a_exp;

        g_y_0_0_0_y_xx_y_xz[i] = -g_0_xx_y_xz[i] + 2.0 * g_yy_xx_y_xz[i] * a_exp;

        g_y_0_0_0_y_xx_y_yy[i] = -g_0_xx_y_yy[i] + 2.0 * g_yy_xx_y_yy[i] * a_exp;

        g_y_0_0_0_y_xx_y_yz[i] = -g_0_xx_y_yz[i] + 2.0 * g_yy_xx_y_yz[i] * a_exp;

        g_y_0_0_0_y_xx_y_zz[i] = -g_0_xx_y_zz[i] + 2.0 * g_yy_xx_y_zz[i] * a_exp;
    }
    // integrals block (444-450)

    #pragma omp simd aligned(g_0_xx_z_xx, g_0_xx_z_xy, g_0_xx_z_xz, g_0_xx_z_yy, g_0_xx_z_yz, g_0_xx_z_zz, g_y_0_0_0_y_xx_z_xx, g_y_0_0_0_y_xx_z_xy, g_y_0_0_0_y_xx_z_xz, g_y_0_0_0_y_xx_z_yy, g_y_0_0_0_y_xx_z_yz, g_y_0_0_0_y_xx_z_zz, g_yy_xx_z_xx, g_yy_xx_z_xy, g_yy_xx_z_xz, g_yy_xx_z_yy, g_yy_xx_z_yz, g_yy_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xx_z_xx[i] = -g_0_xx_z_xx[i] + 2.0 * g_yy_xx_z_xx[i] * a_exp;

        g_y_0_0_0_y_xx_z_xy[i] = -g_0_xx_z_xy[i] + 2.0 * g_yy_xx_z_xy[i] * a_exp;

        g_y_0_0_0_y_xx_z_xz[i] = -g_0_xx_z_xz[i] + 2.0 * g_yy_xx_z_xz[i] * a_exp;

        g_y_0_0_0_y_xx_z_yy[i] = -g_0_xx_z_yy[i] + 2.0 * g_yy_xx_z_yy[i] * a_exp;

        g_y_0_0_0_y_xx_z_yz[i] = -g_0_xx_z_yz[i] + 2.0 * g_yy_xx_z_yz[i] * a_exp;

        g_y_0_0_0_y_xx_z_zz[i] = -g_0_xx_z_zz[i] + 2.0 * g_yy_xx_z_zz[i] * a_exp;
    }
    // integrals block (450-456)

    #pragma omp simd aligned(g_0_xy_x_xx, g_0_xy_x_xy, g_0_xy_x_xz, g_0_xy_x_yy, g_0_xy_x_yz, g_0_xy_x_zz, g_y_0_0_0_y_xy_x_xx, g_y_0_0_0_y_xy_x_xy, g_y_0_0_0_y_xy_x_xz, g_y_0_0_0_y_xy_x_yy, g_y_0_0_0_y_xy_x_yz, g_y_0_0_0_y_xy_x_zz, g_yy_xy_x_xx, g_yy_xy_x_xy, g_yy_xy_x_xz, g_yy_xy_x_yy, g_yy_xy_x_yz, g_yy_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xy_x_xx[i] = -g_0_xy_x_xx[i] + 2.0 * g_yy_xy_x_xx[i] * a_exp;

        g_y_0_0_0_y_xy_x_xy[i] = -g_0_xy_x_xy[i] + 2.0 * g_yy_xy_x_xy[i] * a_exp;

        g_y_0_0_0_y_xy_x_xz[i] = -g_0_xy_x_xz[i] + 2.0 * g_yy_xy_x_xz[i] * a_exp;

        g_y_0_0_0_y_xy_x_yy[i] = -g_0_xy_x_yy[i] + 2.0 * g_yy_xy_x_yy[i] * a_exp;

        g_y_0_0_0_y_xy_x_yz[i] = -g_0_xy_x_yz[i] + 2.0 * g_yy_xy_x_yz[i] * a_exp;

        g_y_0_0_0_y_xy_x_zz[i] = -g_0_xy_x_zz[i] + 2.0 * g_yy_xy_x_zz[i] * a_exp;
    }
    // integrals block (456-462)

    #pragma omp simd aligned(g_0_xy_y_xx, g_0_xy_y_xy, g_0_xy_y_xz, g_0_xy_y_yy, g_0_xy_y_yz, g_0_xy_y_zz, g_y_0_0_0_y_xy_y_xx, g_y_0_0_0_y_xy_y_xy, g_y_0_0_0_y_xy_y_xz, g_y_0_0_0_y_xy_y_yy, g_y_0_0_0_y_xy_y_yz, g_y_0_0_0_y_xy_y_zz, g_yy_xy_y_xx, g_yy_xy_y_xy, g_yy_xy_y_xz, g_yy_xy_y_yy, g_yy_xy_y_yz, g_yy_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xy_y_xx[i] = -g_0_xy_y_xx[i] + 2.0 * g_yy_xy_y_xx[i] * a_exp;

        g_y_0_0_0_y_xy_y_xy[i] = -g_0_xy_y_xy[i] + 2.0 * g_yy_xy_y_xy[i] * a_exp;

        g_y_0_0_0_y_xy_y_xz[i] = -g_0_xy_y_xz[i] + 2.0 * g_yy_xy_y_xz[i] * a_exp;

        g_y_0_0_0_y_xy_y_yy[i] = -g_0_xy_y_yy[i] + 2.0 * g_yy_xy_y_yy[i] * a_exp;

        g_y_0_0_0_y_xy_y_yz[i] = -g_0_xy_y_yz[i] + 2.0 * g_yy_xy_y_yz[i] * a_exp;

        g_y_0_0_0_y_xy_y_zz[i] = -g_0_xy_y_zz[i] + 2.0 * g_yy_xy_y_zz[i] * a_exp;
    }
    // integrals block (462-468)

    #pragma omp simd aligned(g_0_xy_z_xx, g_0_xy_z_xy, g_0_xy_z_xz, g_0_xy_z_yy, g_0_xy_z_yz, g_0_xy_z_zz, g_y_0_0_0_y_xy_z_xx, g_y_0_0_0_y_xy_z_xy, g_y_0_0_0_y_xy_z_xz, g_y_0_0_0_y_xy_z_yy, g_y_0_0_0_y_xy_z_yz, g_y_0_0_0_y_xy_z_zz, g_yy_xy_z_xx, g_yy_xy_z_xy, g_yy_xy_z_xz, g_yy_xy_z_yy, g_yy_xy_z_yz, g_yy_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xy_z_xx[i] = -g_0_xy_z_xx[i] + 2.0 * g_yy_xy_z_xx[i] * a_exp;

        g_y_0_0_0_y_xy_z_xy[i] = -g_0_xy_z_xy[i] + 2.0 * g_yy_xy_z_xy[i] * a_exp;

        g_y_0_0_0_y_xy_z_xz[i] = -g_0_xy_z_xz[i] + 2.0 * g_yy_xy_z_xz[i] * a_exp;

        g_y_0_0_0_y_xy_z_yy[i] = -g_0_xy_z_yy[i] + 2.0 * g_yy_xy_z_yy[i] * a_exp;

        g_y_0_0_0_y_xy_z_yz[i] = -g_0_xy_z_yz[i] + 2.0 * g_yy_xy_z_yz[i] * a_exp;

        g_y_0_0_0_y_xy_z_zz[i] = -g_0_xy_z_zz[i] + 2.0 * g_yy_xy_z_zz[i] * a_exp;
    }
    // integrals block (468-474)

    #pragma omp simd aligned(g_0_xz_x_xx, g_0_xz_x_xy, g_0_xz_x_xz, g_0_xz_x_yy, g_0_xz_x_yz, g_0_xz_x_zz, g_y_0_0_0_y_xz_x_xx, g_y_0_0_0_y_xz_x_xy, g_y_0_0_0_y_xz_x_xz, g_y_0_0_0_y_xz_x_yy, g_y_0_0_0_y_xz_x_yz, g_y_0_0_0_y_xz_x_zz, g_yy_xz_x_xx, g_yy_xz_x_xy, g_yy_xz_x_xz, g_yy_xz_x_yy, g_yy_xz_x_yz, g_yy_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xz_x_xx[i] = -g_0_xz_x_xx[i] + 2.0 * g_yy_xz_x_xx[i] * a_exp;

        g_y_0_0_0_y_xz_x_xy[i] = -g_0_xz_x_xy[i] + 2.0 * g_yy_xz_x_xy[i] * a_exp;

        g_y_0_0_0_y_xz_x_xz[i] = -g_0_xz_x_xz[i] + 2.0 * g_yy_xz_x_xz[i] * a_exp;

        g_y_0_0_0_y_xz_x_yy[i] = -g_0_xz_x_yy[i] + 2.0 * g_yy_xz_x_yy[i] * a_exp;

        g_y_0_0_0_y_xz_x_yz[i] = -g_0_xz_x_yz[i] + 2.0 * g_yy_xz_x_yz[i] * a_exp;

        g_y_0_0_0_y_xz_x_zz[i] = -g_0_xz_x_zz[i] + 2.0 * g_yy_xz_x_zz[i] * a_exp;
    }
    // integrals block (474-480)

    #pragma omp simd aligned(g_0_xz_y_xx, g_0_xz_y_xy, g_0_xz_y_xz, g_0_xz_y_yy, g_0_xz_y_yz, g_0_xz_y_zz, g_y_0_0_0_y_xz_y_xx, g_y_0_0_0_y_xz_y_xy, g_y_0_0_0_y_xz_y_xz, g_y_0_0_0_y_xz_y_yy, g_y_0_0_0_y_xz_y_yz, g_y_0_0_0_y_xz_y_zz, g_yy_xz_y_xx, g_yy_xz_y_xy, g_yy_xz_y_xz, g_yy_xz_y_yy, g_yy_xz_y_yz, g_yy_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xz_y_xx[i] = -g_0_xz_y_xx[i] + 2.0 * g_yy_xz_y_xx[i] * a_exp;

        g_y_0_0_0_y_xz_y_xy[i] = -g_0_xz_y_xy[i] + 2.0 * g_yy_xz_y_xy[i] * a_exp;

        g_y_0_0_0_y_xz_y_xz[i] = -g_0_xz_y_xz[i] + 2.0 * g_yy_xz_y_xz[i] * a_exp;

        g_y_0_0_0_y_xz_y_yy[i] = -g_0_xz_y_yy[i] + 2.0 * g_yy_xz_y_yy[i] * a_exp;

        g_y_0_0_0_y_xz_y_yz[i] = -g_0_xz_y_yz[i] + 2.0 * g_yy_xz_y_yz[i] * a_exp;

        g_y_0_0_0_y_xz_y_zz[i] = -g_0_xz_y_zz[i] + 2.0 * g_yy_xz_y_zz[i] * a_exp;
    }
    // integrals block (480-486)

    #pragma omp simd aligned(g_0_xz_z_xx, g_0_xz_z_xy, g_0_xz_z_xz, g_0_xz_z_yy, g_0_xz_z_yz, g_0_xz_z_zz, g_y_0_0_0_y_xz_z_xx, g_y_0_0_0_y_xz_z_xy, g_y_0_0_0_y_xz_z_xz, g_y_0_0_0_y_xz_z_yy, g_y_0_0_0_y_xz_z_yz, g_y_0_0_0_y_xz_z_zz, g_yy_xz_z_xx, g_yy_xz_z_xy, g_yy_xz_z_xz, g_yy_xz_z_yy, g_yy_xz_z_yz, g_yy_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xz_z_xx[i] = -g_0_xz_z_xx[i] + 2.0 * g_yy_xz_z_xx[i] * a_exp;

        g_y_0_0_0_y_xz_z_xy[i] = -g_0_xz_z_xy[i] + 2.0 * g_yy_xz_z_xy[i] * a_exp;

        g_y_0_0_0_y_xz_z_xz[i] = -g_0_xz_z_xz[i] + 2.0 * g_yy_xz_z_xz[i] * a_exp;

        g_y_0_0_0_y_xz_z_yy[i] = -g_0_xz_z_yy[i] + 2.0 * g_yy_xz_z_yy[i] * a_exp;

        g_y_0_0_0_y_xz_z_yz[i] = -g_0_xz_z_yz[i] + 2.0 * g_yy_xz_z_yz[i] * a_exp;

        g_y_0_0_0_y_xz_z_zz[i] = -g_0_xz_z_zz[i] + 2.0 * g_yy_xz_z_zz[i] * a_exp;
    }
    // integrals block (486-492)

    #pragma omp simd aligned(g_0_yy_x_xx, g_0_yy_x_xy, g_0_yy_x_xz, g_0_yy_x_yy, g_0_yy_x_yz, g_0_yy_x_zz, g_y_0_0_0_y_yy_x_xx, g_y_0_0_0_y_yy_x_xy, g_y_0_0_0_y_yy_x_xz, g_y_0_0_0_y_yy_x_yy, g_y_0_0_0_y_yy_x_yz, g_y_0_0_0_y_yy_x_zz, g_yy_yy_x_xx, g_yy_yy_x_xy, g_yy_yy_x_xz, g_yy_yy_x_yy, g_yy_yy_x_yz, g_yy_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_yy_x_xx[i] = -g_0_yy_x_xx[i] + 2.0 * g_yy_yy_x_xx[i] * a_exp;

        g_y_0_0_0_y_yy_x_xy[i] = -g_0_yy_x_xy[i] + 2.0 * g_yy_yy_x_xy[i] * a_exp;

        g_y_0_0_0_y_yy_x_xz[i] = -g_0_yy_x_xz[i] + 2.0 * g_yy_yy_x_xz[i] * a_exp;

        g_y_0_0_0_y_yy_x_yy[i] = -g_0_yy_x_yy[i] + 2.0 * g_yy_yy_x_yy[i] * a_exp;

        g_y_0_0_0_y_yy_x_yz[i] = -g_0_yy_x_yz[i] + 2.0 * g_yy_yy_x_yz[i] * a_exp;

        g_y_0_0_0_y_yy_x_zz[i] = -g_0_yy_x_zz[i] + 2.0 * g_yy_yy_x_zz[i] * a_exp;
    }
    // integrals block (492-498)

    #pragma omp simd aligned(g_0_yy_y_xx, g_0_yy_y_xy, g_0_yy_y_xz, g_0_yy_y_yy, g_0_yy_y_yz, g_0_yy_y_zz, g_y_0_0_0_y_yy_y_xx, g_y_0_0_0_y_yy_y_xy, g_y_0_0_0_y_yy_y_xz, g_y_0_0_0_y_yy_y_yy, g_y_0_0_0_y_yy_y_yz, g_y_0_0_0_y_yy_y_zz, g_yy_yy_y_xx, g_yy_yy_y_xy, g_yy_yy_y_xz, g_yy_yy_y_yy, g_yy_yy_y_yz, g_yy_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_yy_y_xx[i] = -g_0_yy_y_xx[i] + 2.0 * g_yy_yy_y_xx[i] * a_exp;

        g_y_0_0_0_y_yy_y_xy[i] = -g_0_yy_y_xy[i] + 2.0 * g_yy_yy_y_xy[i] * a_exp;

        g_y_0_0_0_y_yy_y_xz[i] = -g_0_yy_y_xz[i] + 2.0 * g_yy_yy_y_xz[i] * a_exp;

        g_y_0_0_0_y_yy_y_yy[i] = -g_0_yy_y_yy[i] + 2.0 * g_yy_yy_y_yy[i] * a_exp;

        g_y_0_0_0_y_yy_y_yz[i] = -g_0_yy_y_yz[i] + 2.0 * g_yy_yy_y_yz[i] * a_exp;

        g_y_0_0_0_y_yy_y_zz[i] = -g_0_yy_y_zz[i] + 2.0 * g_yy_yy_y_zz[i] * a_exp;
    }
    // integrals block (498-504)

    #pragma omp simd aligned(g_0_yy_z_xx, g_0_yy_z_xy, g_0_yy_z_xz, g_0_yy_z_yy, g_0_yy_z_yz, g_0_yy_z_zz, g_y_0_0_0_y_yy_z_xx, g_y_0_0_0_y_yy_z_xy, g_y_0_0_0_y_yy_z_xz, g_y_0_0_0_y_yy_z_yy, g_y_0_0_0_y_yy_z_yz, g_y_0_0_0_y_yy_z_zz, g_yy_yy_z_xx, g_yy_yy_z_xy, g_yy_yy_z_xz, g_yy_yy_z_yy, g_yy_yy_z_yz, g_yy_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_yy_z_xx[i] = -g_0_yy_z_xx[i] + 2.0 * g_yy_yy_z_xx[i] * a_exp;

        g_y_0_0_0_y_yy_z_xy[i] = -g_0_yy_z_xy[i] + 2.0 * g_yy_yy_z_xy[i] * a_exp;

        g_y_0_0_0_y_yy_z_xz[i] = -g_0_yy_z_xz[i] + 2.0 * g_yy_yy_z_xz[i] * a_exp;

        g_y_0_0_0_y_yy_z_yy[i] = -g_0_yy_z_yy[i] + 2.0 * g_yy_yy_z_yy[i] * a_exp;

        g_y_0_0_0_y_yy_z_yz[i] = -g_0_yy_z_yz[i] + 2.0 * g_yy_yy_z_yz[i] * a_exp;

        g_y_0_0_0_y_yy_z_zz[i] = -g_0_yy_z_zz[i] + 2.0 * g_yy_yy_z_zz[i] * a_exp;
    }
    // integrals block (504-510)

    #pragma omp simd aligned(g_0_yz_x_xx, g_0_yz_x_xy, g_0_yz_x_xz, g_0_yz_x_yy, g_0_yz_x_yz, g_0_yz_x_zz, g_y_0_0_0_y_yz_x_xx, g_y_0_0_0_y_yz_x_xy, g_y_0_0_0_y_yz_x_xz, g_y_0_0_0_y_yz_x_yy, g_y_0_0_0_y_yz_x_yz, g_y_0_0_0_y_yz_x_zz, g_yy_yz_x_xx, g_yy_yz_x_xy, g_yy_yz_x_xz, g_yy_yz_x_yy, g_yy_yz_x_yz, g_yy_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_yz_x_xx[i] = -g_0_yz_x_xx[i] + 2.0 * g_yy_yz_x_xx[i] * a_exp;

        g_y_0_0_0_y_yz_x_xy[i] = -g_0_yz_x_xy[i] + 2.0 * g_yy_yz_x_xy[i] * a_exp;

        g_y_0_0_0_y_yz_x_xz[i] = -g_0_yz_x_xz[i] + 2.0 * g_yy_yz_x_xz[i] * a_exp;

        g_y_0_0_0_y_yz_x_yy[i] = -g_0_yz_x_yy[i] + 2.0 * g_yy_yz_x_yy[i] * a_exp;

        g_y_0_0_0_y_yz_x_yz[i] = -g_0_yz_x_yz[i] + 2.0 * g_yy_yz_x_yz[i] * a_exp;

        g_y_0_0_0_y_yz_x_zz[i] = -g_0_yz_x_zz[i] + 2.0 * g_yy_yz_x_zz[i] * a_exp;
    }
    // integrals block (510-516)

    #pragma omp simd aligned(g_0_yz_y_xx, g_0_yz_y_xy, g_0_yz_y_xz, g_0_yz_y_yy, g_0_yz_y_yz, g_0_yz_y_zz, g_y_0_0_0_y_yz_y_xx, g_y_0_0_0_y_yz_y_xy, g_y_0_0_0_y_yz_y_xz, g_y_0_0_0_y_yz_y_yy, g_y_0_0_0_y_yz_y_yz, g_y_0_0_0_y_yz_y_zz, g_yy_yz_y_xx, g_yy_yz_y_xy, g_yy_yz_y_xz, g_yy_yz_y_yy, g_yy_yz_y_yz, g_yy_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_yz_y_xx[i] = -g_0_yz_y_xx[i] + 2.0 * g_yy_yz_y_xx[i] * a_exp;

        g_y_0_0_0_y_yz_y_xy[i] = -g_0_yz_y_xy[i] + 2.0 * g_yy_yz_y_xy[i] * a_exp;

        g_y_0_0_0_y_yz_y_xz[i] = -g_0_yz_y_xz[i] + 2.0 * g_yy_yz_y_xz[i] * a_exp;

        g_y_0_0_0_y_yz_y_yy[i] = -g_0_yz_y_yy[i] + 2.0 * g_yy_yz_y_yy[i] * a_exp;

        g_y_0_0_0_y_yz_y_yz[i] = -g_0_yz_y_yz[i] + 2.0 * g_yy_yz_y_yz[i] * a_exp;

        g_y_0_0_0_y_yz_y_zz[i] = -g_0_yz_y_zz[i] + 2.0 * g_yy_yz_y_zz[i] * a_exp;
    }
    // integrals block (516-522)

    #pragma omp simd aligned(g_0_yz_z_xx, g_0_yz_z_xy, g_0_yz_z_xz, g_0_yz_z_yy, g_0_yz_z_yz, g_0_yz_z_zz, g_y_0_0_0_y_yz_z_xx, g_y_0_0_0_y_yz_z_xy, g_y_0_0_0_y_yz_z_xz, g_y_0_0_0_y_yz_z_yy, g_y_0_0_0_y_yz_z_yz, g_y_0_0_0_y_yz_z_zz, g_yy_yz_z_xx, g_yy_yz_z_xy, g_yy_yz_z_xz, g_yy_yz_z_yy, g_yy_yz_z_yz, g_yy_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_yz_z_xx[i] = -g_0_yz_z_xx[i] + 2.0 * g_yy_yz_z_xx[i] * a_exp;

        g_y_0_0_0_y_yz_z_xy[i] = -g_0_yz_z_xy[i] + 2.0 * g_yy_yz_z_xy[i] * a_exp;

        g_y_0_0_0_y_yz_z_xz[i] = -g_0_yz_z_xz[i] + 2.0 * g_yy_yz_z_xz[i] * a_exp;

        g_y_0_0_0_y_yz_z_yy[i] = -g_0_yz_z_yy[i] + 2.0 * g_yy_yz_z_yy[i] * a_exp;

        g_y_0_0_0_y_yz_z_yz[i] = -g_0_yz_z_yz[i] + 2.0 * g_yy_yz_z_yz[i] * a_exp;

        g_y_0_0_0_y_yz_z_zz[i] = -g_0_yz_z_zz[i] + 2.0 * g_yy_yz_z_zz[i] * a_exp;
    }
    // integrals block (522-528)

    #pragma omp simd aligned(g_0_zz_x_xx, g_0_zz_x_xy, g_0_zz_x_xz, g_0_zz_x_yy, g_0_zz_x_yz, g_0_zz_x_zz, g_y_0_0_0_y_zz_x_xx, g_y_0_0_0_y_zz_x_xy, g_y_0_0_0_y_zz_x_xz, g_y_0_0_0_y_zz_x_yy, g_y_0_0_0_y_zz_x_yz, g_y_0_0_0_y_zz_x_zz, g_yy_zz_x_xx, g_yy_zz_x_xy, g_yy_zz_x_xz, g_yy_zz_x_yy, g_yy_zz_x_yz, g_yy_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_zz_x_xx[i] = -g_0_zz_x_xx[i] + 2.0 * g_yy_zz_x_xx[i] * a_exp;

        g_y_0_0_0_y_zz_x_xy[i] = -g_0_zz_x_xy[i] + 2.0 * g_yy_zz_x_xy[i] * a_exp;

        g_y_0_0_0_y_zz_x_xz[i] = -g_0_zz_x_xz[i] + 2.0 * g_yy_zz_x_xz[i] * a_exp;

        g_y_0_0_0_y_zz_x_yy[i] = -g_0_zz_x_yy[i] + 2.0 * g_yy_zz_x_yy[i] * a_exp;

        g_y_0_0_0_y_zz_x_yz[i] = -g_0_zz_x_yz[i] + 2.0 * g_yy_zz_x_yz[i] * a_exp;

        g_y_0_0_0_y_zz_x_zz[i] = -g_0_zz_x_zz[i] + 2.0 * g_yy_zz_x_zz[i] * a_exp;
    }
    // integrals block (528-534)

    #pragma omp simd aligned(g_0_zz_y_xx, g_0_zz_y_xy, g_0_zz_y_xz, g_0_zz_y_yy, g_0_zz_y_yz, g_0_zz_y_zz, g_y_0_0_0_y_zz_y_xx, g_y_0_0_0_y_zz_y_xy, g_y_0_0_0_y_zz_y_xz, g_y_0_0_0_y_zz_y_yy, g_y_0_0_0_y_zz_y_yz, g_y_0_0_0_y_zz_y_zz, g_yy_zz_y_xx, g_yy_zz_y_xy, g_yy_zz_y_xz, g_yy_zz_y_yy, g_yy_zz_y_yz, g_yy_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_zz_y_xx[i] = -g_0_zz_y_xx[i] + 2.0 * g_yy_zz_y_xx[i] * a_exp;

        g_y_0_0_0_y_zz_y_xy[i] = -g_0_zz_y_xy[i] + 2.0 * g_yy_zz_y_xy[i] * a_exp;

        g_y_0_0_0_y_zz_y_xz[i] = -g_0_zz_y_xz[i] + 2.0 * g_yy_zz_y_xz[i] * a_exp;

        g_y_0_0_0_y_zz_y_yy[i] = -g_0_zz_y_yy[i] + 2.0 * g_yy_zz_y_yy[i] * a_exp;

        g_y_0_0_0_y_zz_y_yz[i] = -g_0_zz_y_yz[i] + 2.0 * g_yy_zz_y_yz[i] * a_exp;

        g_y_0_0_0_y_zz_y_zz[i] = -g_0_zz_y_zz[i] + 2.0 * g_yy_zz_y_zz[i] * a_exp;
    }
    // integrals block (534-540)

    #pragma omp simd aligned(g_0_zz_z_xx, g_0_zz_z_xy, g_0_zz_z_xz, g_0_zz_z_yy, g_0_zz_z_yz, g_0_zz_z_zz, g_y_0_0_0_y_zz_z_xx, g_y_0_0_0_y_zz_z_xy, g_y_0_0_0_y_zz_z_xz, g_y_0_0_0_y_zz_z_yy, g_y_0_0_0_y_zz_z_yz, g_y_0_0_0_y_zz_z_zz, g_yy_zz_z_xx, g_yy_zz_z_xy, g_yy_zz_z_xz, g_yy_zz_z_yy, g_yy_zz_z_yz, g_yy_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_zz_z_xx[i] = -g_0_zz_z_xx[i] + 2.0 * g_yy_zz_z_xx[i] * a_exp;

        g_y_0_0_0_y_zz_z_xy[i] = -g_0_zz_z_xy[i] + 2.0 * g_yy_zz_z_xy[i] * a_exp;

        g_y_0_0_0_y_zz_z_xz[i] = -g_0_zz_z_xz[i] + 2.0 * g_yy_zz_z_xz[i] * a_exp;

        g_y_0_0_0_y_zz_z_yy[i] = -g_0_zz_z_yy[i] + 2.0 * g_yy_zz_z_yy[i] * a_exp;

        g_y_0_0_0_y_zz_z_yz[i] = -g_0_zz_z_yz[i] + 2.0 * g_yy_zz_z_yz[i] * a_exp;

        g_y_0_0_0_y_zz_z_zz[i] = -g_0_zz_z_zz[i] + 2.0 * g_yy_zz_z_zz[i] * a_exp;
    }
    // integrals block (540-546)

    #pragma omp simd aligned(g_y_0_0_0_z_xx_x_xx, g_y_0_0_0_z_xx_x_xy, g_y_0_0_0_z_xx_x_xz, g_y_0_0_0_z_xx_x_yy, g_y_0_0_0_z_xx_x_yz, g_y_0_0_0_z_xx_x_zz, g_yz_xx_x_xx, g_yz_xx_x_xy, g_yz_xx_x_xz, g_yz_xx_x_yy, g_yz_xx_x_yz, g_yz_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xx_x_xx[i] = 2.0 * g_yz_xx_x_xx[i] * a_exp;

        g_y_0_0_0_z_xx_x_xy[i] = 2.0 * g_yz_xx_x_xy[i] * a_exp;

        g_y_0_0_0_z_xx_x_xz[i] = 2.0 * g_yz_xx_x_xz[i] * a_exp;

        g_y_0_0_0_z_xx_x_yy[i] = 2.0 * g_yz_xx_x_yy[i] * a_exp;

        g_y_0_0_0_z_xx_x_yz[i] = 2.0 * g_yz_xx_x_yz[i] * a_exp;

        g_y_0_0_0_z_xx_x_zz[i] = 2.0 * g_yz_xx_x_zz[i] * a_exp;
    }
    // integrals block (546-552)

    #pragma omp simd aligned(g_y_0_0_0_z_xx_y_xx, g_y_0_0_0_z_xx_y_xy, g_y_0_0_0_z_xx_y_xz, g_y_0_0_0_z_xx_y_yy, g_y_0_0_0_z_xx_y_yz, g_y_0_0_0_z_xx_y_zz, g_yz_xx_y_xx, g_yz_xx_y_xy, g_yz_xx_y_xz, g_yz_xx_y_yy, g_yz_xx_y_yz, g_yz_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xx_y_xx[i] = 2.0 * g_yz_xx_y_xx[i] * a_exp;

        g_y_0_0_0_z_xx_y_xy[i] = 2.0 * g_yz_xx_y_xy[i] * a_exp;

        g_y_0_0_0_z_xx_y_xz[i] = 2.0 * g_yz_xx_y_xz[i] * a_exp;

        g_y_0_0_0_z_xx_y_yy[i] = 2.0 * g_yz_xx_y_yy[i] * a_exp;

        g_y_0_0_0_z_xx_y_yz[i] = 2.0 * g_yz_xx_y_yz[i] * a_exp;

        g_y_0_0_0_z_xx_y_zz[i] = 2.0 * g_yz_xx_y_zz[i] * a_exp;
    }
    // integrals block (552-558)

    #pragma omp simd aligned(g_y_0_0_0_z_xx_z_xx, g_y_0_0_0_z_xx_z_xy, g_y_0_0_0_z_xx_z_xz, g_y_0_0_0_z_xx_z_yy, g_y_0_0_0_z_xx_z_yz, g_y_0_0_0_z_xx_z_zz, g_yz_xx_z_xx, g_yz_xx_z_xy, g_yz_xx_z_xz, g_yz_xx_z_yy, g_yz_xx_z_yz, g_yz_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xx_z_xx[i] = 2.0 * g_yz_xx_z_xx[i] * a_exp;

        g_y_0_0_0_z_xx_z_xy[i] = 2.0 * g_yz_xx_z_xy[i] * a_exp;

        g_y_0_0_0_z_xx_z_xz[i] = 2.0 * g_yz_xx_z_xz[i] * a_exp;

        g_y_0_0_0_z_xx_z_yy[i] = 2.0 * g_yz_xx_z_yy[i] * a_exp;

        g_y_0_0_0_z_xx_z_yz[i] = 2.0 * g_yz_xx_z_yz[i] * a_exp;

        g_y_0_0_0_z_xx_z_zz[i] = 2.0 * g_yz_xx_z_zz[i] * a_exp;
    }
    // integrals block (558-564)

    #pragma omp simd aligned(g_y_0_0_0_z_xy_x_xx, g_y_0_0_0_z_xy_x_xy, g_y_0_0_0_z_xy_x_xz, g_y_0_0_0_z_xy_x_yy, g_y_0_0_0_z_xy_x_yz, g_y_0_0_0_z_xy_x_zz, g_yz_xy_x_xx, g_yz_xy_x_xy, g_yz_xy_x_xz, g_yz_xy_x_yy, g_yz_xy_x_yz, g_yz_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xy_x_xx[i] = 2.0 * g_yz_xy_x_xx[i] * a_exp;

        g_y_0_0_0_z_xy_x_xy[i] = 2.0 * g_yz_xy_x_xy[i] * a_exp;

        g_y_0_0_0_z_xy_x_xz[i] = 2.0 * g_yz_xy_x_xz[i] * a_exp;

        g_y_0_0_0_z_xy_x_yy[i] = 2.0 * g_yz_xy_x_yy[i] * a_exp;

        g_y_0_0_0_z_xy_x_yz[i] = 2.0 * g_yz_xy_x_yz[i] * a_exp;

        g_y_0_0_0_z_xy_x_zz[i] = 2.0 * g_yz_xy_x_zz[i] * a_exp;
    }
    // integrals block (564-570)

    #pragma omp simd aligned(g_y_0_0_0_z_xy_y_xx, g_y_0_0_0_z_xy_y_xy, g_y_0_0_0_z_xy_y_xz, g_y_0_0_0_z_xy_y_yy, g_y_0_0_0_z_xy_y_yz, g_y_0_0_0_z_xy_y_zz, g_yz_xy_y_xx, g_yz_xy_y_xy, g_yz_xy_y_xz, g_yz_xy_y_yy, g_yz_xy_y_yz, g_yz_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xy_y_xx[i] = 2.0 * g_yz_xy_y_xx[i] * a_exp;

        g_y_0_0_0_z_xy_y_xy[i] = 2.0 * g_yz_xy_y_xy[i] * a_exp;

        g_y_0_0_0_z_xy_y_xz[i] = 2.0 * g_yz_xy_y_xz[i] * a_exp;

        g_y_0_0_0_z_xy_y_yy[i] = 2.0 * g_yz_xy_y_yy[i] * a_exp;

        g_y_0_0_0_z_xy_y_yz[i] = 2.0 * g_yz_xy_y_yz[i] * a_exp;

        g_y_0_0_0_z_xy_y_zz[i] = 2.0 * g_yz_xy_y_zz[i] * a_exp;
    }
    // integrals block (570-576)

    #pragma omp simd aligned(g_y_0_0_0_z_xy_z_xx, g_y_0_0_0_z_xy_z_xy, g_y_0_0_0_z_xy_z_xz, g_y_0_0_0_z_xy_z_yy, g_y_0_0_0_z_xy_z_yz, g_y_0_0_0_z_xy_z_zz, g_yz_xy_z_xx, g_yz_xy_z_xy, g_yz_xy_z_xz, g_yz_xy_z_yy, g_yz_xy_z_yz, g_yz_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xy_z_xx[i] = 2.0 * g_yz_xy_z_xx[i] * a_exp;

        g_y_0_0_0_z_xy_z_xy[i] = 2.0 * g_yz_xy_z_xy[i] * a_exp;

        g_y_0_0_0_z_xy_z_xz[i] = 2.0 * g_yz_xy_z_xz[i] * a_exp;

        g_y_0_0_0_z_xy_z_yy[i] = 2.0 * g_yz_xy_z_yy[i] * a_exp;

        g_y_0_0_0_z_xy_z_yz[i] = 2.0 * g_yz_xy_z_yz[i] * a_exp;

        g_y_0_0_0_z_xy_z_zz[i] = 2.0 * g_yz_xy_z_zz[i] * a_exp;
    }
    // integrals block (576-582)

    #pragma omp simd aligned(g_y_0_0_0_z_xz_x_xx, g_y_0_0_0_z_xz_x_xy, g_y_0_0_0_z_xz_x_xz, g_y_0_0_0_z_xz_x_yy, g_y_0_0_0_z_xz_x_yz, g_y_0_0_0_z_xz_x_zz, g_yz_xz_x_xx, g_yz_xz_x_xy, g_yz_xz_x_xz, g_yz_xz_x_yy, g_yz_xz_x_yz, g_yz_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xz_x_xx[i] = 2.0 * g_yz_xz_x_xx[i] * a_exp;

        g_y_0_0_0_z_xz_x_xy[i] = 2.0 * g_yz_xz_x_xy[i] * a_exp;

        g_y_0_0_0_z_xz_x_xz[i] = 2.0 * g_yz_xz_x_xz[i] * a_exp;

        g_y_0_0_0_z_xz_x_yy[i] = 2.0 * g_yz_xz_x_yy[i] * a_exp;

        g_y_0_0_0_z_xz_x_yz[i] = 2.0 * g_yz_xz_x_yz[i] * a_exp;

        g_y_0_0_0_z_xz_x_zz[i] = 2.0 * g_yz_xz_x_zz[i] * a_exp;
    }
    // integrals block (582-588)

    #pragma omp simd aligned(g_y_0_0_0_z_xz_y_xx, g_y_0_0_0_z_xz_y_xy, g_y_0_0_0_z_xz_y_xz, g_y_0_0_0_z_xz_y_yy, g_y_0_0_0_z_xz_y_yz, g_y_0_0_0_z_xz_y_zz, g_yz_xz_y_xx, g_yz_xz_y_xy, g_yz_xz_y_xz, g_yz_xz_y_yy, g_yz_xz_y_yz, g_yz_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xz_y_xx[i] = 2.0 * g_yz_xz_y_xx[i] * a_exp;

        g_y_0_0_0_z_xz_y_xy[i] = 2.0 * g_yz_xz_y_xy[i] * a_exp;

        g_y_0_0_0_z_xz_y_xz[i] = 2.0 * g_yz_xz_y_xz[i] * a_exp;

        g_y_0_0_0_z_xz_y_yy[i] = 2.0 * g_yz_xz_y_yy[i] * a_exp;

        g_y_0_0_0_z_xz_y_yz[i] = 2.0 * g_yz_xz_y_yz[i] * a_exp;

        g_y_0_0_0_z_xz_y_zz[i] = 2.0 * g_yz_xz_y_zz[i] * a_exp;
    }
    // integrals block (588-594)

    #pragma omp simd aligned(g_y_0_0_0_z_xz_z_xx, g_y_0_0_0_z_xz_z_xy, g_y_0_0_0_z_xz_z_xz, g_y_0_0_0_z_xz_z_yy, g_y_0_0_0_z_xz_z_yz, g_y_0_0_0_z_xz_z_zz, g_yz_xz_z_xx, g_yz_xz_z_xy, g_yz_xz_z_xz, g_yz_xz_z_yy, g_yz_xz_z_yz, g_yz_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xz_z_xx[i] = 2.0 * g_yz_xz_z_xx[i] * a_exp;

        g_y_0_0_0_z_xz_z_xy[i] = 2.0 * g_yz_xz_z_xy[i] * a_exp;

        g_y_0_0_0_z_xz_z_xz[i] = 2.0 * g_yz_xz_z_xz[i] * a_exp;

        g_y_0_0_0_z_xz_z_yy[i] = 2.0 * g_yz_xz_z_yy[i] * a_exp;

        g_y_0_0_0_z_xz_z_yz[i] = 2.0 * g_yz_xz_z_yz[i] * a_exp;

        g_y_0_0_0_z_xz_z_zz[i] = 2.0 * g_yz_xz_z_zz[i] * a_exp;
    }
    // integrals block (594-600)

    #pragma omp simd aligned(g_y_0_0_0_z_yy_x_xx, g_y_0_0_0_z_yy_x_xy, g_y_0_0_0_z_yy_x_xz, g_y_0_0_0_z_yy_x_yy, g_y_0_0_0_z_yy_x_yz, g_y_0_0_0_z_yy_x_zz, g_yz_yy_x_xx, g_yz_yy_x_xy, g_yz_yy_x_xz, g_yz_yy_x_yy, g_yz_yy_x_yz, g_yz_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_yy_x_xx[i] = 2.0 * g_yz_yy_x_xx[i] * a_exp;

        g_y_0_0_0_z_yy_x_xy[i] = 2.0 * g_yz_yy_x_xy[i] * a_exp;

        g_y_0_0_0_z_yy_x_xz[i] = 2.0 * g_yz_yy_x_xz[i] * a_exp;

        g_y_0_0_0_z_yy_x_yy[i] = 2.0 * g_yz_yy_x_yy[i] * a_exp;

        g_y_0_0_0_z_yy_x_yz[i] = 2.0 * g_yz_yy_x_yz[i] * a_exp;

        g_y_0_0_0_z_yy_x_zz[i] = 2.0 * g_yz_yy_x_zz[i] * a_exp;
    }
    // integrals block (600-606)

    #pragma omp simd aligned(g_y_0_0_0_z_yy_y_xx, g_y_0_0_0_z_yy_y_xy, g_y_0_0_0_z_yy_y_xz, g_y_0_0_0_z_yy_y_yy, g_y_0_0_0_z_yy_y_yz, g_y_0_0_0_z_yy_y_zz, g_yz_yy_y_xx, g_yz_yy_y_xy, g_yz_yy_y_xz, g_yz_yy_y_yy, g_yz_yy_y_yz, g_yz_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_yy_y_xx[i] = 2.0 * g_yz_yy_y_xx[i] * a_exp;

        g_y_0_0_0_z_yy_y_xy[i] = 2.0 * g_yz_yy_y_xy[i] * a_exp;

        g_y_0_0_0_z_yy_y_xz[i] = 2.0 * g_yz_yy_y_xz[i] * a_exp;

        g_y_0_0_0_z_yy_y_yy[i] = 2.0 * g_yz_yy_y_yy[i] * a_exp;

        g_y_0_0_0_z_yy_y_yz[i] = 2.0 * g_yz_yy_y_yz[i] * a_exp;

        g_y_0_0_0_z_yy_y_zz[i] = 2.0 * g_yz_yy_y_zz[i] * a_exp;
    }
    // integrals block (606-612)

    #pragma omp simd aligned(g_y_0_0_0_z_yy_z_xx, g_y_0_0_0_z_yy_z_xy, g_y_0_0_0_z_yy_z_xz, g_y_0_0_0_z_yy_z_yy, g_y_0_0_0_z_yy_z_yz, g_y_0_0_0_z_yy_z_zz, g_yz_yy_z_xx, g_yz_yy_z_xy, g_yz_yy_z_xz, g_yz_yy_z_yy, g_yz_yy_z_yz, g_yz_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_yy_z_xx[i] = 2.0 * g_yz_yy_z_xx[i] * a_exp;

        g_y_0_0_0_z_yy_z_xy[i] = 2.0 * g_yz_yy_z_xy[i] * a_exp;

        g_y_0_0_0_z_yy_z_xz[i] = 2.0 * g_yz_yy_z_xz[i] * a_exp;

        g_y_0_0_0_z_yy_z_yy[i] = 2.0 * g_yz_yy_z_yy[i] * a_exp;

        g_y_0_0_0_z_yy_z_yz[i] = 2.0 * g_yz_yy_z_yz[i] * a_exp;

        g_y_0_0_0_z_yy_z_zz[i] = 2.0 * g_yz_yy_z_zz[i] * a_exp;
    }
    // integrals block (612-618)

    #pragma omp simd aligned(g_y_0_0_0_z_yz_x_xx, g_y_0_0_0_z_yz_x_xy, g_y_0_0_0_z_yz_x_xz, g_y_0_0_0_z_yz_x_yy, g_y_0_0_0_z_yz_x_yz, g_y_0_0_0_z_yz_x_zz, g_yz_yz_x_xx, g_yz_yz_x_xy, g_yz_yz_x_xz, g_yz_yz_x_yy, g_yz_yz_x_yz, g_yz_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_yz_x_xx[i] = 2.0 * g_yz_yz_x_xx[i] * a_exp;

        g_y_0_0_0_z_yz_x_xy[i] = 2.0 * g_yz_yz_x_xy[i] * a_exp;

        g_y_0_0_0_z_yz_x_xz[i] = 2.0 * g_yz_yz_x_xz[i] * a_exp;

        g_y_0_0_0_z_yz_x_yy[i] = 2.0 * g_yz_yz_x_yy[i] * a_exp;

        g_y_0_0_0_z_yz_x_yz[i] = 2.0 * g_yz_yz_x_yz[i] * a_exp;

        g_y_0_0_0_z_yz_x_zz[i] = 2.0 * g_yz_yz_x_zz[i] * a_exp;
    }
    // integrals block (618-624)

    #pragma omp simd aligned(g_y_0_0_0_z_yz_y_xx, g_y_0_0_0_z_yz_y_xy, g_y_0_0_0_z_yz_y_xz, g_y_0_0_0_z_yz_y_yy, g_y_0_0_0_z_yz_y_yz, g_y_0_0_0_z_yz_y_zz, g_yz_yz_y_xx, g_yz_yz_y_xy, g_yz_yz_y_xz, g_yz_yz_y_yy, g_yz_yz_y_yz, g_yz_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_yz_y_xx[i] = 2.0 * g_yz_yz_y_xx[i] * a_exp;

        g_y_0_0_0_z_yz_y_xy[i] = 2.0 * g_yz_yz_y_xy[i] * a_exp;

        g_y_0_0_0_z_yz_y_xz[i] = 2.0 * g_yz_yz_y_xz[i] * a_exp;

        g_y_0_0_0_z_yz_y_yy[i] = 2.0 * g_yz_yz_y_yy[i] * a_exp;

        g_y_0_0_0_z_yz_y_yz[i] = 2.0 * g_yz_yz_y_yz[i] * a_exp;

        g_y_0_0_0_z_yz_y_zz[i] = 2.0 * g_yz_yz_y_zz[i] * a_exp;
    }
    // integrals block (624-630)

    #pragma omp simd aligned(g_y_0_0_0_z_yz_z_xx, g_y_0_0_0_z_yz_z_xy, g_y_0_0_0_z_yz_z_xz, g_y_0_0_0_z_yz_z_yy, g_y_0_0_0_z_yz_z_yz, g_y_0_0_0_z_yz_z_zz, g_yz_yz_z_xx, g_yz_yz_z_xy, g_yz_yz_z_xz, g_yz_yz_z_yy, g_yz_yz_z_yz, g_yz_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_yz_z_xx[i] = 2.0 * g_yz_yz_z_xx[i] * a_exp;

        g_y_0_0_0_z_yz_z_xy[i] = 2.0 * g_yz_yz_z_xy[i] * a_exp;

        g_y_0_0_0_z_yz_z_xz[i] = 2.0 * g_yz_yz_z_xz[i] * a_exp;

        g_y_0_0_0_z_yz_z_yy[i] = 2.0 * g_yz_yz_z_yy[i] * a_exp;

        g_y_0_0_0_z_yz_z_yz[i] = 2.0 * g_yz_yz_z_yz[i] * a_exp;

        g_y_0_0_0_z_yz_z_zz[i] = 2.0 * g_yz_yz_z_zz[i] * a_exp;
    }
    // integrals block (630-636)

    #pragma omp simd aligned(g_y_0_0_0_z_zz_x_xx, g_y_0_0_0_z_zz_x_xy, g_y_0_0_0_z_zz_x_xz, g_y_0_0_0_z_zz_x_yy, g_y_0_0_0_z_zz_x_yz, g_y_0_0_0_z_zz_x_zz, g_yz_zz_x_xx, g_yz_zz_x_xy, g_yz_zz_x_xz, g_yz_zz_x_yy, g_yz_zz_x_yz, g_yz_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_zz_x_xx[i] = 2.0 * g_yz_zz_x_xx[i] * a_exp;

        g_y_0_0_0_z_zz_x_xy[i] = 2.0 * g_yz_zz_x_xy[i] * a_exp;

        g_y_0_0_0_z_zz_x_xz[i] = 2.0 * g_yz_zz_x_xz[i] * a_exp;

        g_y_0_0_0_z_zz_x_yy[i] = 2.0 * g_yz_zz_x_yy[i] * a_exp;

        g_y_0_0_0_z_zz_x_yz[i] = 2.0 * g_yz_zz_x_yz[i] * a_exp;

        g_y_0_0_0_z_zz_x_zz[i] = 2.0 * g_yz_zz_x_zz[i] * a_exp;
    }
    // integrals block (636-642)

    #pragma omp simd aligned(g_y_0_0_0_z_zz_y_xx, g_y_0_0_0_z_zz_y_xy, g_y_0_0_0_z_zz_y_xz, g_y_0_0_0_z_zz_y_yy, g_y_0_0_0_z_zz_y_yz, g_y_0_0_0_z_zz_y_zz, g_yz_zz_y_xx, g_yz_zz_y_xy, g_yz_zz_y_xz, g_yz_zz_y_yy, g_yz_zz_y_yz, g_yz_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_zz_y_xx[i] = 2.0 * g_yz_zz_y_xx[i] * a_exp;

        g_y_0_0_0_z_zz_y_xy[i] = 2.0 * g_yz_zz_y_xy[i] * a_exp;

        g_y_0_0_0_z_zz_y_xz[i] = 2.0 * g_yz_zz_y_xz[i] * a_exp;

        g_y_0_0_0_z_zz_y_yy[i] = 2.0 * g_yz_zz_y_yy[i] * a_exp;

        g_y_0_0_0_z_zz_y_yz[i] = 2.0 * g_yz_zz_y_yz[i] * a_exp;

        g_y_0_0_0_z_zz_y_zz[i] = 2.0 * g_yz_zz_y_zz[i] * a_exp;
    }
    // integrals block (642-648)

    #pragma omp simd aligned(g_y_0_0_0_z_zz_z_xx, g_y_0_0_0_z_zz_z_xy, g_y_0_0_0_z_zz_z_xz, g_y_0_0_0_z_zz_z_yy, g_y_0_0_0_z_zz_z_yz, g_y_0_0_0_z_zz_z_zz, g_yz_zz_z_xx, g_yz_zz_z_xy, g_yz_zz_z_xz, g_yz_zz_z_yy, g_yz_zz_z_yz, g_yz_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_zz_z_xx[i] = 2.0 * g_yz_zz_z_xx[i] * a_exp;

        g_y_0_0_0_z_zz_z_xy[i] = 2.0 * g_yz_zz_z_xy[i] * a_exp;

        g_y_0_0_0_z_zz_z_xz[i] = 2.0 * g_yz_zz_z_xz[i] * a_exp;

        g_y_0_0_0_z_zz_z_yy[i] = 2.0 * g_yz_zz_z_yy[i] * a_exp;

        g_y_0_0_0_z_zz_z_yz[i] = 2.0 * g_yz_zz_z_yz[i] * a_exp;

        g_y_0_0_0_z_zz_z_zz[i] = 2.0 * g_yz_zz_z_zz[i] * a_exp;
    }
    // integrals block (648-654)

    #pragma omp simd aligned(g_xz_xx_x_xx, g_xz_xx_x_xy, g_xz_xx_x_xz, g_xz_xx_x_yy, g_xz_xx_x_yz, g_xz_xx_x_zz, g_z_0_0_0_x_xx_x_xx, g_z_0_0_0_x_xx_x_xy, g_z_0_0_0_x_xx_x_xz, g_z_0_0_0_x_xx_x_yy, g_z_0_0_0_x_xx_x_yz, g_z_0_0_0_x_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xx_x_xx[i] = 2.0 * g_xz_xx_x_xx[i] * a_exp;

        g_z_0_0_0_x_xx_x_xy[i] = 2.0 * g_xz_xx_x_xy[i] * a_exp;

        g_z_0_0_0_x_xx_x_xz[i] = 2.0 * g_xz_xx_x_xz[i] * a_exp;

        g_z_0_0_0_x_xx_x_yy[i] = 2.0 * g_xz_xx_x_yy[i] * a_exp;

        g_z_0_0_0_x_xx_x_yz[i] = 2.0 * g_xz_xx_x_yz[i] * a_exp;

        g_z_0_0_0_x_xx_x_zz[i] = 2.0 * g_xz_xx_x_zz[i] * a_exp;
    }
    // integrals block (654-660)

    #pragma omp simd aligned(g_xz_xx_y_xx, g_xz_xx_y_xy, g_xz_xx_y_xz, g_xz_xx_y_yy, g_xz_xx_y_yz, g_xz_xx_y_zz, g_z_0_0_0_x_xx_y_xx, g_z_0_0_0_x_xx_y_xy, g_z_0_0_0_x_xx_y_xz, g_z_0_0_0_x_xx_y_yy, g_z_0_0_0_x_xx_y_yz, g_z_0_0_0_x_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xx_y_xx[i] = 2.0 * g_xz_xx_y_xx[i] * a_exp;

        g_z_0_0_0_x_xx_y_xy[i] = 2.0 * g_xz_xx_y_xy[i] * a_exp;

        g_z_0_0_0_x_xx_y_xz[i] = 2.0 * g_xz_xx_y_xz[i] * a_exp;

        g_z_0_0_0_x_xx_y_yy[i] = 2.0 * g_xz_xx_y_yy[i] * a_exp;

        g_z_0_0_0_x_xx_y_yz[i] = 2.0 * g_xz_xx_y_yz[i] * a_exp;

        g_z_0_0_0_x_xx_y_zz[i] = 2.0 * g_xz_xx_y_zz[i] * a_exp;
    }
    // integrals block (660-666)

    #pragma omp simd aligned(g_xz_xx_z_xx, g_xz_xx_z_xy, g_xz_xx_z_xz, g_xz_xx_z_yy, g_xz_xx_z_yz, g_xz_xx_z_zz, g_z_0_0_0_x_xx_z_xx, g_z_0_0_0_x_xx_z_xy, g_z_0_0_0_x_xx_z_xz, g_z_0_0_0_x_xx_z_yy, g_z_0_0_0_x_xx_z_yz, g_z_0_0_0_x_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xx_z_xx[i] = 2.0 * g_xz_xx_z_xx[i] * a_exp;

        g_z_0_0_0_x_xx_z_xy[i] = 2.0 * g_xz_xx_z_xy[i] * a_exp;

        g_z_0_0_0_x_xx_z_xz[i] = 2.0 * g_xz_xx_z_xz[i] * a_exp;

        g_z_0_0_0_x_xx_z_yy[i] = 2.0 * g_xz_xx_z_yy[i] * a_exp;

        g_z_0_0_0_x_xx_z_yz[i] = 2.0 * g_xz_xx_z_yz[i] * a_exp;

        g_z_0_0_0_x_xx_z_zz[i] = 2.0 * g_xz_xx_z_zz[i] * a_exp;
    }
    // integrals block (666-672)

    #pragma omp simd aligned(g_xz_xy_x_xx, g_xz_xy_x_xy, g_xz_xy_x_xz, g_xz_xy_x_yy, g_xz_xy_x_yz, g_xz_xy_x_zz, g_z_0_0_0_x_xy_x_xx, g_z_0_0_0_x_xy_x_xy, g_z_0_0_0_x_xy_x_xz, g_z_0_0_0_x_xy_x_yy, g_z_0_0_0_x_xy_x_yz, g_z_0_0_0_x_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xy_x_xx[i] = 2.0 * g_xz_xy_x_xx[i] * a_exp;

        g_z_0_0_0_x_xy_x_xy[i] = 2.0 * g_xz_xy_x_xy[i] * a_exp;

        g_z_0_0_0_x_xy_x_xz[i] = 2.0 * g_xz_xy_x_xz[i] * a_exp;

        g_z_0_0_0_x_xy_x_yy[i] = 2.0 * g_xz_xy_x_yy[i] * a_exp;

        g_z_0_0_0_x_xy_x_yz[i] = 2.0 * g_xz_xy_x_yz[i] * a_exp;

        g_z_0_0_0_x_xy_x_zz[i] = 2.0 * g_xz_xy_x_zz[i] * a_exp;
    }
    // integrals block (672-678)

    #pragma omp simd aligned(g_xz_xy_y_xx, g_xz_xy_y_xy, g_xz_xy_y_xz, g_xz_xy_y_yy, g_xz_xy_y_yz, g_xz_xy_y_zz, g_z_0_0_0_x_xy_y_xx, g_z_0_0_0_x_xy_y_xy, g_z_0_0_0_x_xy_y_xz, g_z_0_0_0_x_xy_y_yy, g_z_0_0_0_x_xy_y_yz, g_z_0_0_0_x_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xy_y_xx[i] = 2.0 * g_xz_xy_y_xx[i] * a_exp;

        g_z_0_0_0_x_xy_y_xy[i] = 2.0 * g_xz_xy_y_xy[i] * a_exp;

        g_z_0_0_0_x_xy_y_xz[i] = 2.0 * g_xz_xy_y_xz[i] * a_exp;

        g_z_0_0_0_x_xy_y_yy[i] = 2.0 * g_xz_xy_y_yy[i] * a_exp;

        g_z_0_0_0_x_xy_y_yz[i] = 2.0 * g_xz_xy_y_yz[i] * a_exp;

        g_z_0_0_0_x_xy_y_zz[i] = 2.0 * g_xz_xy_y_zz[i] * a_exp;
    }
    // integrals block (678-684)

    #pragma omp simd aligned(g_xz_xy_z_xx, g_xz_xy_z_xy, g_xz_xy_z_xz, g_xz_xy_z_yy, g_xz_xy_z_yz, g_xz_xy_z_zz, g_z_0_0_0_x_xy_z_xx, g_z_0_0_0_x_xy_z_xy, g_z_0_0_0_x_xy_z_xz, g_z_0_0_0_x_xy_z_yy, g_z_0_0_0_x_xy_z_yz, g_z_0_0_0_x_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xy_z_xx[i] = 2.0 * g_xz_xy_z_xx[i] * a_exp;

        g_z_0_0_0_x_xy_z_xy[i] = 2.0 * g_xz_xy_z_xy[i] * a_exp;

        g_z_0_0_0_x_xy_z_xz[i] = 2.0 * g_xz_xy_z_xz[i] * a_exp;

        g_z_0_0_0_x_xy_z_yy[i] = 2.0 * g_xz_xy_z_yy[i] * a_exp;

        g_z_0_0_0_x_xy_z_yz[i] = 2.0 * g_xz_xy_z_yz[i] * a_exp;

        g_z_0_0_0_x_xy_z_zz[i] = 2.0 * g_xz_xy_z_zz[i] * a_exp;
    }
    // integrals block (684-690)

    #pragma omp simd aligned(g_xz_xz_x_xx, g_xz_xz_x_xy, g_xz_xz_x_xz, g_xz_xz_x_yy, g_xz_xz_x_yz, g_xz_xz_x_zz, g_z_0_0_0_x_xz_x_xx, g_z_0_0_0_x_xz_x_xy, g_z_0_0_0_x_xz_x_xz, g_z_0_0_0_x_xz_x_yy, g_z_0_0_0_x_xz_x_yz, g_z_0_0_0_x_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xz_x_xx[i] = 2.0 * g_xz_xz_x_xx[i] * a_exp;

        g_z_0_0_0_x_xz_x_xy[i] = 2.0 * g_xz_xz_x_xy[i] * a_exp;

        g_z_0_0_0_x_xz_x_xz[i] = 2.0 * g_xz_xz_x_xz[i] * a_exp;

        g_z_0_0_0_x_xz_x_yy[i] = 2.0 * g_xz_xz_x_yy[i] * a_exp;

        g_z_0_0_0_x_xz_x_yz[i] = 2.0 * g_xz_xz_x_yz[i] * a_exp;

        g_z_0_0_0_x_xz_x_zz[i] = 2.0 * g_xz_xz_x_zz[i] * a_exp;
    }
    // integrals block (690-696)

    #pragma omp simd aligned(g_xz_xz_y_xx, g_xz_xz_y_xy, g_xz_xz_y_xz, g_xz_xz_y_yy, g_xz_xz_y_yz, g_xz_xz_y_zz, g_z_0_0_0_x_xz_y_xx, g_z_0_0_0_x_xz_y_xy, g_z_0_0_0_x_xz_y_xz, g_z_0_0_0_x_xz_y_yy, g_z_0_0_0_x_xz_y_yz, g_z_0_0_0_x_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xz_y_xx[i] = 2.0 * g_xz_xz_y_xx[i] * a_exp;

        g_z_0_0_0_x_xz_y_xy[i] = 2.0 * g_xz_xz_y_xy[i] * a_exp;

        g_z_0_0_0_x_xz_y_xz[i] = 2.0 * g_xz_xz_y_xz[i] * a_exp;

        g_z_0_0_0_x_xz_y_yy[i] = 2.0 * g_xz_xz_y_yy[i] * a_exp;

        g_z_0_0_0_x_xz_y_yz[i] = 2.0 * g_xz_xz_y_yz[i] * a_exp;

        g_z_0_0_0_x_xz_y_zz[i] = 2.0 * g_xz_xz_y_zz[i] * a_exp;
    }
    // integrals block (696-702)

    #pragma omp simd aligned(g_xz_xz_z_xx, g_xz_xz_z_xy, g_xz_xz_z_xz, g_xz_xz_z_yy, g_xz_xz_z_yz, g_xz_xz_z_zz, g_z_0_0_0_x_xz_z_xx, g_z_0_0_0_x_xz_z_xy, g_z_0_0_0_x_xz_z_xz, g_z_0_0_0_x_xz_z_yy, g_z_0_0_0_x_xz_z_yz, g_z_0_0_0_x_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xz_z_xx[i] = 2.0 * g_xz_xz_z_xx[i] * a_exp;

        g_z_0_0_0_x_xz_z_xy[i] = 2.0 * g_xz_xz_z_xy[i] * a_exp;

        g_z_0_0_0_x_xz_z_xz[i] = 2.0 * g_xz_xz_z_xz[i] * a_exp;

        g_z_0_0_0_x_xz_z_yy[i] = 2.0 * g_xz_xz_z_yy[i] * a_exp;

        g_z_0_0_0_x_xz_z_yz[i] = 2.0 * g_xz_xz_z_yz[i] * a_exp;

        g_z_0_0_0_x_xz_z_zz[i] = 2.0 * g_xz_xz_z_zz[i] * a_exp;
    }
    // integrals block (702-708)

    #pragma omp simd aligned(g_xz_yy_x_xx, g_xz_yy_x_xy, g_xz_yy_x_xz, g_xz_yy_x_yy, g_xz_yy_x_yz, g_xz_yy_x_zz, g_z_0_0_0_x_yy_x_xx, g_z_0_0_0_x_yy_x_xy, g_z_0_0_0_x_yy_x_xz, g_z_0_0_0_x_yy_x_yy, g_z_0_0_0_x_yy_x_yz, g_z_0_0_0_x_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_yy_x_xx[i] = 2.0 * g_xz_yy_x_xx[i] * a_exp;

        g_z_0_0_0_x_yy_x_xy[i] = 2.0 * g_xz_yy_x_xy[i] * a_exp;

        g_z_0_0_0_x_yy_x_xz[i] = 2.0 * g_xz_yy_x_xz[i] * a_exp;

        g_z_0_0_0_x_yy_x_yy[i] = 2.0 * g_xz_yy_x_yy[i] * a_exp;

        g_z_0_0_0_x_yy_x_yz[i] = 2.0 * g_xz_yy_x_yz[i] * a_exp;

        g_z_0_0_0_x_yy_x_zz[i] = 2.0 * g_xz_yy_x_zz[i] * a_exp;
    }
    // integrals block (708-714)

    #pragma omp simd aligned(g_xz_yy_y_xx, g_xz_yy_y_xy, g_xz_yy_y_xz, g_xz_yy_y_yy, g_xz_yy_y_yz, g_xz_yy_y_zz, g_z_0_0_0_x_yy_y_xx, g_z_0_0_0_x_yy_y_xy, g_z_0_0_0_x_yy_y_xz, g_z_0_0_0_x_yy_y_yy, g_z_0_0_0_x_yy_y_yz, g_z_0_0_0_x_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_yy_y_xx[i] = 2.0 * g_xz_yy_y_xx[i] * a_exp;

        g_z_0_0_0_x_yy_y_xy[i] = 2.0 * g_xz_yy_y_xy[i] * a_exp;

        g_z_0_0_0_x_yy_y_xz[i] = 2.0 * g_xz_yy_y_xz[i] * a_exp;

        g_z_0_0_0_x_yy_y_yy[i] = 2.0 * g_xz_yy_y_yy[i] * a_exp;

        g_z_0_0_0_x_yy_y_yz[i] = 2.0 * g_xz_yy_y_yz[i] * a_exp;

        g_z_0_0_0_x_yy_y_zz[i] = 2.0 * g_xz_yy_y_zz[i] * a_exp;
    }
    // integrals block (714-720)

    #pragma omp simd aligned(g_xz_yy_z_xx, g_xz_yy_z_xy, g_xz_yy_z_xz, g_xz_yy_z_yy, g_xz_yy_z_yz, g_xz_yy_z_zz, g_z_0_0_0_x_yy_z_xx, g_z_0_0_0_x_yy_z_xy, g_z_0_0_0_x_yy_z_xz, g_z_0_0_0_x_yy_z_yy, g_z_0_0_0_x_yy_z_yz, g_z_0_0_0_x_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_yy_z_xx[i] = 2.0 * g_xz_yy_z_xx[i] * a_exp;

        g_z_0_0_0_x_yy_z_xy[i] = 2.0 * g_xz_yy_z_xy[i] * a_exp;

        g_z_0_0_0_x_yy_z_xz[i] = 2.0 * g_xz_yy_z_xz[i] * a_exp;

        g_z_0_0_0_x_yy_z_yy[i] = 2.0 * g_xz_yy_z_yy[i] * a_exp;

        g_z_0_0_0_x_yy_z_yz[i] = 2.0 * g_xz_yy_z_yz[i] * a_exp;

        g_z_0_0_0_x_yy_z_zz[i] = 2.0 * g_xz_yy_z_zz[i] * a_exp;
    }
    // integrals block (720-726)

    #pragma omp simd aligned(g_xz_yz_x_xx, g_xz_yz_x_xy, g_xz_yz_x_xz, g_xz_yz_x_yy, g_xz_yz_x_yz, g_xz_yz_x_zz, g_z_0_0_0_x_yz_x_xx, g_z_0_0_0_x_yz_x_xy, g_z_0_0_0_x_yz_x_xz, g_z_0_0_0_x_yz_x_yy, g_z_0_0_0_x_yz_x_yz, g_z_0_0_0_x_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_yz_x_xx[i] = 2.0 * g_xz_yz_x_xx[i] * a_exp;

        g_z_0_0_0_x_yz_x_xy[i] = 2.0 * g_xz_yz_x_xy[i] * a_exp;

        g_z_0_0_0_x_yz_x_xz[i] = 2.0 * g_xz_yz_x_xz[i] * a_exp;

        g_z_0_0_0_x_yz_x_yy[i] = 2.0 * g_xz_yz_x_yy[i] * a_exp;

        g_z_0_0_0_x_yz_x_yz[i] = 2.0 * g_xz_yz_x_yz[i] * a_exp;

        g_z_0_0_0_x_yz_x_zz[i] = 2.0 * g_xz_yz_x_zz[i] * a_exp;
    }
    // integrals block (726-732)

    #pragma omp simd aligned(g_xz_yz_y_xx, g_xz_yz_y_xy, g_xz_yz_y_xz, g_xz_yz_y_yy, g_xz_yz_y_yz, g_xz_yz_y_zz, g_z_0_0_0_x_yz_y_xx, g_z_0_0_0_x_yz_y_xy, g_z_0_0_0_x_yz_y_xz, g_z_0_0_0_x_yz_y_yy, g_z_0_0_0_x_yz_y_yz, g_z_0_0_0_x_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_yz_y_xx[i] = 2.0 * g_xz_yz_y_xx[i] * a_exp;

        g_z_0_0_0_x_yz_y_xy[i] = 2.0 * g_xz_yz_y_xy[i] * a_exp;

        g_z_0_0_0_x_yz_y_xz[i] = 2.0 * g_xz_yz_y_xz[i] * a_exp;

        g_z_0_0_0_x_yz_y_yy[i] = 2.0 * g_xz_yz_y_yy[i] * a_exp;

        g_z_0_0_0_x_yz_y_yz[i] = 2.0 * g_xz_yz_y_yz[i] * a_exp;

        g_z_0_0_0_x_yz_y_zz[i] = 2.0 * g_xz_yz_y_zz[i] * a_exp;
    }
    // integrals block (732-738)

    #pragma omp simd aligned(g_xz_yz_z_xx, g_xz_yz_z_xy, g_xz_yz_z_xz, g_xz_yz_z_yy, g_xz_yz_z_yz, g_xz_yz_z_zz, g_z_0_0_0_x_yz_z_xx, g_z_0_0_0_x_yz_z_xy, g_z_0_0_0_x_yz_z_xz, g_z_0_0_0_x_yz_z_yy, g_z_0_0_0_x_yz_z_yz, g_z_0_0_0_x_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_yz_z_xx[i] = 2.0 * g_xz_yz_z_xx[i] * a_exp;

        g_z_0_0_0_x_yz_z_xy[i] = 2.0 * g_xz_yz_z_xy[i] * a_exp;

        g_z_0_0_0_x_yz_z_xz[i] = 2.0 * g_xz_yz_z_xz[i] * a_exp;

        g_z_0_0_0_x_yz_z_yy[i] = 2.0 * g_xz_yz_z_yy[i] * a_exp;

        g_z_0_0_0_x_yz_z_yz[i] = 2.0 * g_xz_yz_z_yz[i] * a_exp;

        g_z_0_0_0_x_yz_z_zz[i] = 2.0 * g_xz_yz_z_zz[i] * a_exp;
    }
    // integrals block (738-744)

    #pragma omp simd aligned(g_xz_zz_x_xx, g_xz_zz_x_xy, g_xz_zz_x_xz, g_xz_zz_x_yy, g_xz_zz_x_yz, g_xz_zz_x_zz, g_z_0_0_0_x_zz_x_xx, g_z_0_0_0_x_zz_x_xy, g_z_0_0_0_x_zz_x_xz, g_z_0_0_0_x_zz_x_yy, g_z_0_0_0_x_zz_x_yz, g_z_0_0_0_x_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_zz_x_xx[i] = 2.0 * g_xz_zz_x_xx[i] * a_exp;

        g_z_0_0_0_x_zz_x_xy[i] = 2.0 * g_xz_zz_x_xy[i] * a_exp;

        g_z_0_0_0_x_zz_x_xz[i] = 2.0 * g_xz_zz_x_xz[i] * a_exp;

        g_z_0_0_0_x_zz_x_yy[i] = 2.0 * g_xz_zz_x_yy[i] * a_exp;

        g_z_0_0_0_x_zz_x_yz[i] = 2.0 * g_xz_zz_x_yz[i] * a_exp;

        g_z_0_0_0_x_zz_x_zz[i] = 2.0 * g_xz_zz_x_zz[i] * a_exp;
    }
    // integrals block (744-750)

    #pragma omp simd aligned(g_xz_zz_y_xx, g_xz_zz_y_xy, g_xz_zz_y_xz, g_xz_zz_y_yy, g_xz_zz_y_yz, g_xz_zz_y_zz, g_z_0_0_0_x_zz_y_xx, g_z_0_0_0_x_zz_y_xy, g_z_0_0_0_x_zz_y_xz, g_z_0_0_0_x_zz_y_yy, g_z_0_0_0_x_zz_y_yz, g_z_0_0_0_x_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_zz_y_xx[i] = 2.0 * g_xz_zz_y_xx[i] * a_exp;

        g_z_0_0_0_x_zz_y_xy[i] = 2.0 * g_xz_zz_y_xy[i] * a_exp;

        g_z_0_0_0_x_zz_y_xz[i] = 2.0 * g_xz_zz_y_xz[i] * a_exp;

        g_z_0_0_0_x_zz_y_yy[i] = 2.0 * g_xz_zz_y_yy[i] * a_exp;

        g_z_0_0_0_x_zz_y_yz[i] = 2.0 * g_xz_zz_y_yz[i] * a_exp;

        g_z_0_0_0_x_zz_y_zz[i] = 2.0 * g_xz_zz_y_zz[i] * a_exp;
    }
    // integrals block (750-756)

    #pragma omp simd aligned(g_xz_zz_z_xx, g_xz_zz_z_xy, g_xz_zz_z_xz, g_xz_zz_z_yy, g_xz_zz_z_yz, g_xz_zz_z_zz, g_z_0_0_0_x_zz_z_xx, g_z_0_0_0_x_zz_z_xy, g_z_0_0_0_x_zz_z_xz, g_z_0_0_0_x_zz_z_yy, g_z_0_0_0_x_zz_z_yz, g_z_0_0_0_x_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_zz_z_xx[i] = 2.0 * g_xz_zz_z_xx[i] * a_exp;

        g_z_0_0_0_x_zz_z_xy[i] = 2.0 * g_xz_zz_z_xy[i] * a_exp;

        g_z_0_0_0_x_zz_z_xz[i] = 2.0 * g_xz_zz_z_xz[i] * a_exp;

        g_z_0_0_0_x_zz_z_yy[i] = 2.0 * g_xz_zz_z_yy[i] * a_exp;

        g_z_0_0_0_x_zz_z_yz[i] = 2.0 * g_xz_zz_z_yz[i] * a_exp;

        g_z_0_0_0_x_zz_z_zz[i] = 2.0 * g_xz_zz_z_zz[i] * a_exp;
    }
    // integrals block (756-762)

    #pragma omp simd aligned(g_yz_xx_x_xx, g_yz_xx_x_xy, g_yz_xx_x_xz, g_yz_xx_x_yy, g_yz_xx_x_yz, g_yz_xx_x_zz, g_z_0_0_0_y_xx_x_xx, g_z_0_0_0_y_xx_x_xy, g_z_0_0_0_y_xx_x_xz, g_z_0_0_0_y_xx_x_yy, g_z_0_0_0_y_xx_x_yz, g_z_0_0_0_y_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xx_x_xx[i] = 2.0 * g_yz_xx_x_xx[i] * a_exp;

        g_z_0_0_0_y_xx_x_xy[i] = 2.0 * g_yz_xx_x_xy[i] * a_exp;

        g_z_0_0_0_y_xx_x_xz[i] = 2.0 * g_yz_xx_x_xz[i] * a_exp;

        g_z_0_0_0_y_xx_x_yy[i] = 2.0 * g_yz_xx_x_yy[i] * a_exp;

        g_z_0_0_0_y_xx_x_yz[i] = 2.0 * g_yz_xx_x_yz[i] * a_exp;

        g_z_0_0_0_y_xx_x_zz[i] = 2.0 * g_yz_xx_x_zz[i] * a_exp;
    }
    // integrals block (762-768)

    #pragma omp simd aligned(g_yz_xx_y_xx, g_yz_xx_y_xy, g_yz_xx_y_xz, g_yz_xx_y_yy, g_yz_xx_y_yz, g_yz_xx_y_zz, g_z_0_0_0_y_xx_y_xx, g_z_0_0_0_y_xx_y_xy, g_z_0_0_0_y_xx_y_xz, g_z_0_0_0_y_xx_y_yy, g_z_0_0_0_y_xx_y_yz, g_z_0_0_0_y_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xx_y_xx[i] = 2.0 * g_yz_xx_y_xx[i] * a_exp;

        g_z_0_0_0_y_xx_y_xy[i] = 2.0 * g_yz_xx_y_xy[i] * a_exp;

        g_z_0_0_0_y_xx_y_xz[i] = 2.0 * g_yz_xx_y_xz[i] * a_exp;

        g_z_0_0_0_y_xx_y_yy[i] = 2.0 * g_yz_xx_y_yy[i] * a_exp;

        g_z_0_0_0_y_xx_y_yz[i] = 2.0 * g_yz_xx_y_yz[i] * a_exp;

        g_z_0_0_0_y_xx_y_zz[i] = 2.0 * g_yz_xx_y_zz[i] * a_exp;
    }
    // integrals block (768-774)

    #pragma omp simd aligned(g_yz_xx_z_xx, g_yz_xx_z_xy, g_yz_xx_z_xz, g_yz_xx_z_yy, g_yz_xx_z_yz, g_yz_xx_z_zz, g_z_0_0_0_y_xx_z_xx, g_z_0_0_0_y_xx_z_xy, g_z_0_0_0_y_xx_z_xz, g_z_0_0_0_y_xx_z_yy, g_z_0_0_0_y_xx_z_yz, g_z_0_0_0_y_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xx_z_xx[i] = 2.0 * g_yz_xx_z_xx[i] * a_exp;

        g_z_0_0_0_y_xx_z_xy[i] = 2.0 * g_yz_xx_z_xy[i] * a_exp;

        g_z_0_0_0_y_xx_z_xz[i] = 2.0 * g_yz_xx_z_xz[i] * a_exp;

        g_z_0_0_0_y_xx_z_yy[i] = 2.0 * g_yz_xx_z_yy[i] * a_exp;

        g_z_0_0_0_y_xx_z_yz[i] = 2.0 * g_yz_xx_z_yz[i] * a_exp;

        g_z_0_0_0_y_xx_z_zz[i] = 2.0 * g_yz_xx_z_zz[i] * a_exp;
    }
    // integrals block (774-780)

    #pragma omp simd aligned(g_yz_xy_x_xx, g_yz_xy_x_xy, g_yz_xy_x_xz, g_yz_xy_x_yy, g_yz_xy_x_yz, g_yz_xy_x_zz, g_z_0_0_0_y_xy_x_xx, g_z_0_0_0_y_xy_x_xy, g_z_0_0_0_y_xy_x_xz, g_z_0_0_0_y_xy_x_yy, g_z_0_0_0_y_xy_x_yz, g_z_0_0_0_y_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xy_x_xx[i] = 2.0 * g_yz_xy_x_xx[i] * a_exp;

        g_z_0_0_0_y_xy_x_xy[i] = 2.0 * g_yz_xy_x_xy[i] * a_exp;

        g_z_0_0_0_y_xy_x_xz[i] = 2.0 * g_yz_xy_x_xz[i] * a_exp;

        g_z_0_0_0_y_xy_x_yy[i] = 2.0 * g_yz_xy_x_yy[i] * a_exp;

        g_z_0_0_0_y_xy_x_yz[i] = 2.0 * g_yz_xy_x_yz[i] * a_exp;

        g_z_0_0_0_y_xy_x_zz[i] = 2.0 * g_yz_xy_x_zz[i] * a_exp;
    }
    // integrals block (780-786)

    #pragma omp simd aligned(g_yz_xy_y_xx, g_yz_xy_y_xy, g_yz_xy_y_xz, g_yz_xy_y_yy, g_yz_xy_y_yz, g_yz_xy_y_zz, g_z_0_0_0_y_xy_y_xx, g_z_0_0_0_y_xy_y_xy, g_z_0_0_0_y_xy_y_xz, g_z_0_0_0_y_xy_y_yy, g_z_0_0_0_y_xy_y_yz, g_z_0_0_0_y_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xy_y_xx[i] = 2.0 * g_yz_xy_y_xx[i] * a_exp;

        g_z_0_0_0_y_xy_y_xy[i] = 2.0 * g_yz_xy_y_xy[i] * a_exp;

        g_z_0_0_0_y_xy_y_xz[i] = 2.0 * g_yz_xy_y_xz[i] * a_exp;

        g_z_0_0_0_y_xy_y_yy[i] = 2.0 * g_yz_xy_y_yy[i] * a_exp;

        g_z_0_0_0_y_xy_y_yz[i] = 2.0 * g_yz_xy_y_yz[i] * a_exp;

        g_z_0_0_0_y_xy_y_zz[i] = 2.0 * g_yz_xy_y_zz[i] * a_exp;
    }
    // integrals block (786-792)

    #pragma omp simd aligned(g_yz_xy_z_xx, g_yz_xy_z_xy, g_yz_xy_z_xz, g_yz_xy_z_yy, g_yz_xy_z_yz, g_yz_xy_z_zz, g_z_0_0_0_y_xy_z_xx, g_z_0_0_0_y_xy_z_xy, g_z_0_0_0_y_xy_z_xz, g_z_0_0_0_y_xy_z_yy, g_z_0_0_0_y_xy_z_yz, g_z_0_0_0_y_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xy_z_xx[i] = 2.0 * g_yz_xy_z_xx[i] * a_exp;

        g_z_0_0_0_y_xy_z_xy[i] = 2.0 * g_yz_xy_z_xy[i] * a_exp;

        g_z_0_0_0_y_xy_z_xz[i] = 2.0 * g_yz_xy_z_xz[i] * a_exp;

        g_z_0_0_0_y_xy_z_yy[i] = 2.0 * g_yz_xy_z_yy[i] * a_exp;

        g_z_0_0_0_y_xy_z_yz[i] = 2.0 * g_yz_xy_z_yz[i] * a_exp;

        g_z_0_0_0_y_xy_z_zz[i] = 2.0 * g_yz_xy_z_zz[i] * a_exp;
    }
    // integrals block (792-798)

    #pragma omp simd aligned(g_yz_xz_x_xx, g_yz_xz_x_xy, g_yz_xz_x_xz, g_yz_xz_x_yy, g_yz_xz_x_yz, g_yz_xz_x_zz, g_z_0_0_0_y_xz_x_xx, g_z_0_0_0_y_xz_x_xy, g_z_0_0_0_y_xz_x_xz, g_z_0_0_0_y_xz_x_yy, g_z_0_0_0_y_xz_x_yz, g_z_0_0_0_y_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xz_x_xx[i] = 2.0 * g_yz_xz_x_xx[i] * a_exp;

        g_z_0_0_0_y_xz_x_xy[i] = 2.0 * g_yz_xz_x_xy[i] * a_exp;

        g_z_0_0_0_y_xz_x_xz[i] = 2.0 * g_yz_xz_x_xz[i] * a_exp;

        g_z_0_0_0_y_xz_x_yy[i] = 2.0 * g_yz_xz_x_yy[i] * a_exp;

        g_z_0_0_0_y_xz_x_yz[i] = 2.0 * g_yz_xz_x_yz[i] * a_exp;

        g_z_0_0_0_y_xz_x_zz[i] = 2.0 * g_yz_xz_x_zz[i] * a_exp;
    }
    // integrals block (798-804)

    #pragma omp simd aligned(g_yz_xz_y_xx, g_yz_xz_y_xy, g_yz_xz_y_xz, g_yz_xz_y_yy, g_yz_xz_y_yz, g_yz_xz_y_zz, g_z_0_0_0_y_xz_y_xx, g_z_0_0_0_y_xz_y_xy, g_z_0_0_0_y_xz_y_xz, g_z_0_0_0_y_xz_y_yy, g_z_0_0_0_y_xz_y_yz, g_z_0_0_0_y_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xz_y_xx[i] = 2.0 * g_yz_xz_y_xx[i] * a_exp;

        g_z_0_0_0_y_xz_y_xy[i] = 2.0 * g_yz_xz_y_xy[i] * a_exp;

        g_z_0_0_0_y_xz_y_xz[i] = 2.0 * g_yz_xz_y_xz[i] * a_exp;

        g_z_0_0_0_y_xz_y_yy[i] = 2.0 * g_yz_xz_y_yy[i] * a_exp;

        g_z_0_0_0_y_xz_y_yz[i] = 2.0 * g_yz_xz_y_yz[i] * a_exp;

        g_z_0_0_0_y_xz_y_zz[i] = 2.0 * g_yz_xz_y_zz[i] * a_exp;
    }
    // integrals block (804-810)

    #pragma omp simd aligned(g_yz_xz_z_xx, g_yz_xz_z_xy, g_yz_xz_z_xz, g_yz_xz_z_yy, g_yz_xz_z_yz, g_yz_xz_z_zz, g_z_0_0_0_y_xz_z_xx, g_z_0_0_0_y_xz_z_xy, g_z_0_0_0_y_xz_z_xz, g_z_0_0_0_y_xz_z_yy, g_z_0_0_0_y_xz_z_yz, g_z_0_0_0_y_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xz_z_xx[i] = 2.0 * g_yz_xz_z_xx[i] * a_exp;

        g_z_0_0_0_y_xz_z_xy[i] = 2.0 * g_yz_xz_z_xy[i] * a_exp;

        g_z_0_0_0_y_xz_z_xz[i] = 2.0 * g_yz_xz_z_xz[i] * a_exp;

        g_z_0_0_0_y_xz_z_yy[i] = 2.0 * g_yz_xz_z_yy[i] * a_exp;

        g_z_0_0_0_y_xz_z_yz[i] = 2.0 * g_yz_xz_z_yz[i] * a_exp;

        g_z_0_0_0_y_xz_z_zz[i] = 2.0 * g_yz_xz_z_zz[i] * a_exp;
    }
    // integrals block (810-816)

    #pragma omp simd aligned(g_yz_yy_x_xx, g_yz_yy_x_xy, g_yz_yy_x_xz, g_yz_yy_x_yy, g_yz_yy_x_yz, g_yz_yy_x_zz, g_z_0_0_0_y_yy_x_xx, g_z_0_0_0_y_yy_x_xy, g_z_0_0_0_y_yy_x_xz, g_z_0_0_0_y_yy_x_yy, g_z_0_0_0_y_yy_x_yz, g_z_0_0_0_y_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_yy_x_xx[i] = 2.0 * g_yz_yy_x_xx[i] * a_exp;

        g_z_0_0_0_y_yy_x_xy[i] = 2.0 * g_yz_yy_x_xy[i] * a_exp;

        g_z_0_0_0_y_yy_x_xz[i] = 2.0 * g_yz_yy_x_xz[i] * a_exp;

        g_z_0_0_0_y_yy_x_yy[i] = 2.0 * g_yz_yy_x_yy[i] * a_exp;

        g_z_0_0_0_y_yy_x_yz[i] = 2.0 * g_yz_yy_x_yz[i] * a_exp;

        g_z_0_0_0_y_yy_x_zz[i] = 2.0 * g_yz_yy_x_zz[i] * a_exp;
    }
    // integrals block (816-822)

    #pragma omp simd aligned(g_yz_yy_y_xx, g_yz_yy_y_xy, g_yz_yy_y_xz, g_yz_yy_y_yy, g_yz_yy_y_yz, g_yz_yy_y_zz, g_z_0_0_0_y_yy_y_xx, g_z_0_0_0_y_yy_y_xy, g_z_0_0_0_y_yy_y_xz, g_z_0_0_0_y_yy_y_yy, g_z_0_0_0_y_yy_y_yz, g_z_0_0_0_y_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_yy_y_xx[i] = 2.0 * g_yz_yy_y_xx[i] * a_exp;

        g_z_0_0_0_y_yy_y_xy[i] = 2.0 * g_yz_yy_y_xy[i] * a_exp;

        g_z_0_0_0_y_yy_y_xz[i] = 2.0 * g_yz_yy_y_xz[i] * a_exp;

        g_z_0_0_0_y_yy_y_yy[i] = 2.0 * g_yz_yy_y_yy[i] * a_exp;

        g_z_0_0_0_y_yy_y_yz[i] = 2.0 * g_yz_yy_y_yz[i] * a_exp;

        g_z_0_0_0_y_yy_y_zz[i] = 2.0 * g_yz_yy_y_zz[i] * a_exp;
    }
    // integrals block (822-828)

    #pragma omp simd aligned(g_yz_yy_z_xx, g_yz_yy_z_xy, g_yz_yy_z_xz, g_yz_yy_z_yy, g_yz_yy_z_yz, g_yz_yy_z_zz, g_z_0_0_0_y_yy_z_xx, g_z_0_0_0_y_yy_z_xy, g_z_0_0_0_y_yy_z_xz, g_z_0_0_0_y_yy_z_yy, g_z_0_0_0_y_yy_z_yz, g_z_0_0_0_y_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_yy_z_xx[i] = 2.0 * g_yz_yy_z_xx[i] * a_exp;

        g_z_0_0_0_y_yy_z_xy[i] = 2.0 * g_yz_yy_z_xy[i] * a_exp;

        g_z_0_0_0_y_yy_z_xz[i] = 2.0 * g_yz_yy_z_xz[i] * a_exp;

        g_z_0_0_0_y_yy_z_yy[i] = 2.0 * g_yz_yy_z_yy[i] * a_exp;

        g_z_0_0_0_y_yy_z_yz[i] = 2.0 * g_yz_yy_z_yz[i] * a_exp;

        g_z_0_0_0_y_yy_z_zz[i] = 2.0 * g_yz_yy_z_zz[i] * a_exp;
    }
    // integrals block (828-834)

    #pragma omp simd aligned(g_yz_yz_x_xx, g_yz_yz_x_xy, g_yz_yz_x_xz, g_yz_yz_x_yy, g_yz_yz_x_yz, g_yz_yz_x_zz, g_z_0_0_0_y_yz_x_xx, g_z_0_0_0_y_yz_x_xy, g_z_0_0_0_y_yz_x_xz, g_z_0_0_0_y_yz_x_yy, g_z_0_0_0_y_yz_x_yz, g_z_0_0_0_y_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_yz_x_xx[i] = 2.0 * g_yz_yz_x_xx[i] * a_exp;

        g_z_0_0_0_y_yz_x_xy[i] = 2.0 * g_yz_yz_x_xy[i] * a_exp;

        g_z_0_0_0_y_yz_x_xz[i] = 2.0 * g_yz_yz_x_xz[i] * a_exp;

        g_z_0_0_0_y_yz_x_yy[i] = 2.0 * g_yz_yz_x_yy[i] * a_exp;

        g_z_0_0_0_y_yz_x_yz[i] = 2.0 * g_yz_yz_x_yz[i] * a_exp;

        g_z_0_0_0_y_yz_x_zz[i] = 2.0 * g_yz_yz_x_zz[i] * a_exp;
    }
    // integrals block (834-840)

    #pragma omp simd aligned(g_yz_yz_y_xx, g_yz_yz_y_xy, g_yz_yz_y_xz, g_yz_yz_y_yy, g_yz_yz_y_yz, g_yz_yz_y_zz, g_z_0_0_0_y_yz_y_xx, g_z_0_0_0_y_yz_y_xy, g_z_0_0_0_y_yz_y_xz, g_z_0_0_0_y_yz_y_yy, g_z_0_0_0_y_yz_y_yz, g_z_0_0_0_y_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_yz_y_xx[i] = 2.0 * g_yz_yz_y_xx[i] * a_exp;

        g_z_0_0_0_y_yz_y_xy[i] = 2.0 * g_yz_yz_y_xy[i] * a_exp;

        g_z_0_0_0_y_yz_y_xz[i] = 2.0 * g_yz_yz_y_xz[i] * a_exp;

        g_z_0_0_0_y_yz_y_yy[i] = 2.0 * g_yz_yz_y_yy[i] * a_exp;

        g_z_0_0_0_y_yz_y_yz[i] = 2.0 * g_yz_yz_y_yz[i] * a_exp;

        g_z_0_0_0_y_yz_y_zz[i] = 2.0 * g_yz_yz_y_zz[i] * a_exp;
    }
    // integrals block (840-846)

    #pragma omp simd aligned(g_yz_yz_z_xx, g_yz_yz_z_xy, g_yz_yz_z_xz, g_yz_yz_z_yy, g_yz_yz_z_yz, g_yz_yz_z_zz, g_z_0_0_0_y_yz_z_xx, g_z_0_0_0_y_yz_z_xy, g_z_0_0_0_y_yz_z_xz, g_z_0_0_0_y_yz_z_yy, g_z_0_0_0_y_yz_z_yz, g_z_0_0_0_y_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_yz_z_xx[i] = 2.0 * g_yz_yz_z_xx[i] * a_exp;

        g_z_0_0_0_y_yz_z_xy[i] = 2.0 * g_yz_yz_z_xy[i] * a_exp;

        g_z_0_0_0_y_yz_z_xz[i] = 2.0 * g_yz_yz_z_xz[i] * a_exp;

        g_z_0_0_0_y_yz_z_yy[i] = 2.0 * g_yz_yz_z_yy[i] * a_exp;

        g_z_0_0_0_y_yz_z_yz[i] = 2.0 * g_yz_yz_z_yz[i] * a_exp;

        g_z_0_0_0_y_yz_z_zz[i] = 2.0 * g_yz_yz_z_zz[i] * a_exp;
    }
    // integrals block (846-852)

    #pragma omp simd aligned(g_yz_zz_x_xx, g_yz_zz_x_xy, g_yz_zz_x_xz, g_yz_zz_x_yy, g_yz_zz_x_yz, g_yz_zz_x_zz, g_z_0_0_0_y_zz_x_xx, g_z_0_0_0_y_zz_x_xy, g_z_0_0_0_y_zz_x_xz, g_z_0_0_0_y_zz_x_yy, g_z_0_0_0_y_zz_x_yz, g_z_0_0_0_y_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_zz_x_xx[i] = 2.0 * g_yz_zz_x_xx[i] * a_exp;

        g_z_0_0_0_y_zz_x_xy[i] = 2.0 * g_yz_zz_x_xy[i] * a_exp;

        g_z_0_0_0_y_zz_x_xz[i] = 2.0 * g_yz_zz_x_xz[i] * a_exp;

        g_z_0_0_0_y_zz_x_yy[i] = 2.0 * g_yz_zz_x_yy[i] * a_exp;

        g_z_0_0_0_y_zz_x_yz[i] = 2.0 * g_yz_zz_x_yz[i] * a_exp;

        g_z_0_0_0_y_zz_x_zz[i] = 2.0 * g_yz_zz_x_zz[i] * a_exp;
    }
    // integrals block (852-858)

    #pragma omp simd aligned(g_yz_zz_y_xx, g_yz_zz_y_xy, g_yz_zz_y_xz, g_yz_zz_y_yy, g_yz_zz_y_yz, g_yz_zz_y_zz, g_z_0_0_0_y_zz_y_xx, g_z_0_0_0_y_zz_y_xy, g_z_0_0_0_y_zz_y_xz, g_z_0_0_0_y_zz_y_yy, g_z_0_0_0_y_zz_y_yz, g_z_0_0_0_y_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_zz_y_xx[i] = 2.0 * g_yz_zz_y_xx[i] * a_exp;

        g_z_0_0_0_y_zz_y_xy[i] = 2.0 * g_yz_zz_y_xy[i] * a_exp;

        g_z_0_0_0_y_zz_y_xz[i] = 2.0 * g_yz_zz_y_xz[i] * a_exp;

        g_z_0_0_0_y_zz_y_yy[i] = 2.0 * g_yz_zz_y_yy[i] * a_exp;

        g_z_0_0_0_y_zz_y_yz[i] = 2.0 * g_yz_zz_y_yz[i] * a_exp;

        g_z_0_0_0_y_zz_y_zz[i] = 2.0 * g_yz_zz_y_zz[i] * a_exp;
    }
    // integrals block (858-864)

    #pragma omp simd aligned(g_yz_zz_z_xx, g_yz_zz_z_xy, g_yz_zz_z_xz, g_yz_zz_z_yy, g_yz_zz_z_yz, g_yz_zz_z_zz, g_z_0_0_0_y_zz_z_xx, g_z_0_0_0_y_zz_z_xy, g_z_0_0_0_y_zz_z_xz, g_z_0_0_0_y_zz_z_yy, g_z_0_0_0_y_zz_z_yz, g_z_0_0_0_y_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_zz_z_xx[i] = 2.0 * g_yz_zz_z_xx[i] * a_exp;

        g_z_0_0_0_y_zz_z_xy[i] = 2.0 * g_yz_zz_z_xy[i] * a_exp;

        g_z_0_0_0_y_zz_z_xz[i] = 2.0 * g_yz_zz_z_xz[i] * a_exp;

        g_z_0_0_0_y_zz_z_yy[i] = 2.0 * g_yz_zz_z_yy[i] * a_exp;

        g_z_0_0_0_y_zz_z_yz[i] = 2.0 * g_yz_zz_z_yz[i] * a_exp;

        g_z_0_0_0_y_zz_z_zz[i] = 2.0 * g_yz_zz_z_zz[i] * a_exp;
    }
    // integrals block (864-870)

    #pragma omp simd aligned(g_0_xx_x_xx, g_0_xx_x_xy, g_0_xx_x_xz, g_0_xx_x_yy, g_0_xx_x_yz, g_0_xx_x_zz, g_z_0_0_0_z_xx_x_xx, g_z_0_0_0_z_xx_x_xy, g_z_0_0_0_z_xx_x_xz, g_z_0_0_0_z_xx_x_yy, g_z_0_0_0_z_xx_x_yz, g_z_0_0_0_z_xx_x_zz, g_zz_xx_x_xx, g_zz_xx_x_xy, g_zz_xx_x_xz, g_zz_xx_x_yy, g_zz_xx_x_yz, g_zz_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xx_x_xx[i] = -g_0_xx_x_xx[i] + 2.0 * g_zz_xx_x_xx[i] * a_exp;

        g_z_0_0_0_z_xx_x_xy[i] = -g_0_xx_x_xy[i] + 2.0 * g_zz_xx_x_xy[i] * a_exp;

        g_z_0_0_0_z_xx_x_xz[i] = -g_0_xx_x_xz[i] + 2.0 * g_zz_xx_x_xz[i] * a_exp;

        g_z_0_0_0_z_xx_x_yy[i] = -g_0_xx_x_yy[i] + 2.0 * g_zz_xx_x_yy[i] * a_exp;

        g_z_0_0_0_z_xx_x_yz[i] = -g_0_xx_x_yz[i] + 2.0 * g_zz_xx_x_yz[i] * a_exp;

        g_z_0_0_0_z_xx_x_zz[i] = -g_0_xx_x_zz[i] + 2.0 * g_zz_xx_x_zz[i] * a_exp;
    }
    // integrals block (870-876)

    #pragma omp simd aligned(g_0_xx_y_xx, g_0_xx_y_xy, g_0_xx_y_xz, g_0_xx_y_yy, g_0_xx_y_yz, g_0_xx_y_zz, g_z_0_0_0_z_xx_y_xx, g_z_0_0_0_z_xx_y_xy, g_z_0_0_0_z_xx_y_xz, g_z_0_0_0_z_xx_y_yy, g_z_0_0_0_z_xx_y_yz, g_z_0_0_0_z_xx_y_zz, g_zz_xx_y_xx, g_zz_xx_y_xy, g_zz_xx_y_xz, g_zz_xx_y_yy, g_zz_xx_y_yz, g_zz_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xx_y_xx[i] = -g_0_xx_y_xx[i] + 2.0 * g_zz_xx_y_xx[i] * a_exp;

        g_z_0_0_0_z_xx_y_xy[i] = -g_0_xx_y_xy[i] + 2.0 * g_zz_xx_y_xy[i] * a_exp;

        g_z_0_0_0_z_xx_y_xz[i] = -g_0_xx_y_xz[i] + 2.0 * g_zz_xx_y_xz[i] * a_exp;

        g_z_0_0_0_z_xx_y_yy[i] = -g_0_xx_y_yy[i] + 2.0 * g_zz_xx_y_yy[i] * a_exp;

        g_z_0_0_0_z_xx_y_yz[i] = -g_0_xx_y_yz[i] + 2.0 * g_zz_xx_y_yz[i] * a_exp;

        g_z_0_0_0_z_xx_y_zz[i] = -g_0_xx_y_zz[i] + 2.0 * g_zz_xx_y_zz[i] * a_exp;
    }
    // integrals block (876-882)

    #pragma omp simd aligned(g_0_xx_z_xx, g_0_xx_z_xy, g_0_xx_z_xz, g_0_xx_z_yy, g_0_xx_z_yz, g_0_xx_z_zz, g_z_0_0_0_z_xx_z_xx, g_z_0_0_0_z_xx_z_xy, g_z_0_0_0_z_xx_z_xz, g_z_0_0_0_z_xx_z_yy, g_z_0_0_0_z_xx_z_yz, g_z_0_0_0_z_xx_z_zz, g_zz_xx_z_xx, g_zz_xx_z_xy, g_zz_xx_z_xz, g_zz_xx_z_yy, g_zz_xx_z_yz, g_zz_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xx_z_xx[i] = -g_0_xx_z_xx[i] + 2.0 * g_zz_xx_z_xx[i] * a_exp;

        g_z_0_0_0_z_xx_z_xy[i] = -g_0_xx_z_xy[i] + 2.0 * g_zz_xx_z_xy[i] * a_exp;

        g_z_0_0_0_z_xx_z_xz[i] = -g_0_xx_z_xz[i] + 2.0 * g_zz_xx_z_xz[i] * a_exp;

        g_z_0_0_0_z_xx_z_yy[i] = -g_0_xx_z_yy[i] + 2.0 * g_zz_xx_z_yy[i] * a_exp;

        g_z_0_0_0_z_xx_z_yz[i] = -g_0_xx_z_yz[i] + 2.0 * g_zz_xx_z_yz[i] * a_exp;

        g_z_0_0_0_z_xx_z_zz[i] = -g_0_xx_z_zz[i] + 2.0 * g_zz_xx_z_zz[i] * a_exp;
    }
    // integrals block (882-888)

    #pragma omp simd aligned(g_0_xy_x_xx, g_0_xy_x_xy, g_0_xy_x_xz, g_0_xy_x_yy, g_0_xy_x_yz, g_0_xy_x_zz, g_z_0_0_0_z_xy_x_xx, g_z_0_0_0_z_xy_x_xy, g_z_0_0_0_z_xy_x_xz, g_z_0_0_0_z_xy_x_yy, g_z_0_0_0_z_xy_x_yz, g_z_0_0_0_z_xy_x_zz, g_zz_xy_x_xx, g_zz_xy_x_xy, g_zz_xy_x_xz, g_zz_xy_x_yy, g_zz_xy_x_yz, g_zz_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xy_x_xx[i] = -g_0_xy_x_xx[i] + 2.0 * g_zz_xy_x_xx[i] * a_exp;

        g_z_0_0_0_z_xy_x_xy[i] = -g_0_xy_x_xy[i] + 2.0 * g_zz_xy_x_xy[i] * a_exp;

        g_z_0_0_0_z_xy_x_xz[i] = -g_0_xy_x_xz[i] + 2.0 * g_zz_xy_x_xz[i] * a_exp;

        g_z_0_0_0_z_xy_x_yy[i] = -g_0_xy_x_yy[i] + 2.0 * g_zz_xy_x_yy[i] * a_exp;

        g_z_0_0_0_z_xy_x_yz[i] = -g_0_xy_x_yz[i] + 2.0 * g_zz_xy_x_yz[i] * a_exp;

        g_z_0_0_0_z_xy_x_zz[i] = -g_0_xy_x_zz[i] + 2.0 * g_zz_xy_x_zz[i] * a_exp;
    }
    // integrals block (888-894)

    #pragma omp simd aligned(g_0_xy_y_xx, g_0_xy_y_xy, g_0_xy_y_xz, g_0_xy_y_yy, g_0_xy_y_yz, g_0_xy_y_zz, g_z_0_0_0_z_xy_y_xx, g_z_0_0_0_z_xy_y_xy, g_z_0_0_0_z_xy_y_xz, g_z_0_0_0_z_xy_y_yy, g_z_0_0_0_z_xy_y_yz, g_z_0_0_0_z_xy_y_zz, g_zz_xy_y_xx, g_zz_xy_y_xy, g_zz_xy_y_xz, g_zz_xy_y_yy, g_zz_xy_y_yz, g_zz_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xy_y_xx[i] = -g_0_xy_y_xx[i] + 2.0 * g_zz_xy_y_xx[i] * a_exp;

        g_z_0_0_0_z_xy_y_xy[i] = -g_0_xy_y_xy[i] + 2.0 * g_zz_xy_y_xy[i] * a_exp;

        g_z_0_0_0_z_xy_y_xz[i] = -g_0_xy_y_xz[i] + 2.0 * g_zz_xy_y_xz[i] * a_exp;

        g_z_0_0_0_z_xy_y_yy[i] = -g_0_xy_y_yy[i] + 2.0 * g_zz_xy_y_yy[i] * a_exp;

        g_z_0_0_0_z_xy_y_yz[i] = -g_0_xy_y_yz[i] + 2.0 * g_zz_xy_y_yz[i] * a_exp;

        g_z_0_0_0_z_xy_y_zz[i] = -g_0_xy_y_zz[i] + 2.0 * g_zz_xy_y_zz[i] * a_exp;
    }
    // integrals block (894-900)

    #pragma omp simd aligned(g_0_xy_z_xx, g_0_xy_z_xy, g_0_xy_z_xz, g_0_xy_z_yy, g_0_xy_z_yz, g_0_xy_z_zz, g_z_0_0_0_z_xy_z_xx, g_z_0_0_0_z_xy_z_xy, g_z_0_0_0_z_xy_z_xz, g_z_0_0_0_z_xy_z_yy, g_z_0_0_0_z_xy_z_yz, g_z_0_0_0_z_xy_z_zz, g_zz_xy_z_xx, g_zz_xy_z_xy, g_zz_xy_z_xz, g_zz_xy_z_yy, g_zz_xy_z_yz, g_zz_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xy_z_xx[i] = -g_0_xy_z_xx[i] + 2.0 * g_zz_xy_z_xx[i] * a_exp;

        g_z_0_0_0_z_xy_z_xy[i] = -g_0_xy_z_xy[i] + 2.0 * g_zz_xy_z_xy[i] * a_exp;

        g_z_0_0_0_z_xy_z_xz[i] = -g_0_xy_z_xz[i] + 2.0 * g_zz_xy_z_xz[i] * a_exp;

        g_z_0_0_0_z_xy_z_yy[i] = -g_0_xy_z_yy[i] + 2.0 * g_zz_xy_z_yy[i] * a_exp;

        g_z_0_0_0_z_xy_z_yz[i] = -g_0_xy_z_yz[i] + 2.0 * g_zz_xy_z_yz[i] * a_exp;

        g_z_0_0_0_z_xy_z_zz[i] = -g_0_xy_z_zz[i] + 2.0 * g_zz_xy_z_zz[i] * a_exp;
    }
    // integrals block (900-906)

    #pragma omp simd aligned(g_0_xz_x_xx, g_0_xz_x_xy, g_0_xz_x_xz, g_0_xz_x_yy, g_0_xz_x_yz, g_0_xz_x_zz, g_z_0_0_0_z_xz_x_xx, g_z_0_0_0_z_xz_x_xy, g_z_0_0_0_z_xz_x_xz, g_z_0_0_0_z_xz_x_yy, g_z_0_0_0_z_xz_x_yz, g_z_0_0_0_z_xz_x_zz, g_zz_xz_x_xx, g_zz_xz_x_xy, g_zz_xz_x_xz, g_zz_xz_x_yy, g_zz_xz_x_yz, g_zz_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xz_x_xx[i] = -g_0_xz_x_xx[i] + 2.0 * g_zz_xz_x_xx[i] * a_exp;

        g_z_0_0_0_z_xz_x_xy[i] = -g_0_xz_x_xy[i] + 2.0 * g_zz_xz_x_xy[i] * a_exp;

        g_z_0_0_0_z_xz_x_xz[i] = -g_0_xz_x_xz[i] + 2.0 * g_zz_xz_x_xz[i] * a_exp;

        g_z_0_0_0_z_xz_x_yy[i] = -g_0_xz_x_yy[i] + 2.0 * g_zz_xz_x_yy[i] * a_exp;

        g_z_0_0_0_z_xz_x_yz[i] = -g_0_xz_x_yz[i] + 2.0 * g_zz_xz_x_yz[i] * a_exp;

        g_z_0_0_0_z_xz_x_zz[i] = -g_0_xz_x_zz[i] + 2.0 * g_zz_xz_x_zz[i] * a_exp;
    }
    // integrals block (906-912)

    #pragma omp simd aligned(g_0_xz_y_xx, g_0_xz_y_xy, g_0_xz_y_xz, g_0_xz_y_yy, g_0_xz_y_yz, g_0_xz_y_zz, g_z_0_0_0_z_xz_y_xx, g_z_0_0_0_z_xz_y_xy, g_z_0_0_0_z_xz_y_xz, g_z_0_0_0_z_xz_y_yy, g_z_0_0_0_z_xz_y_yz, g_z_0_0_0_z_xz_y_zz, g_zz_xz_y_xx, g_zz_xz_y_xy, g_zz_xz_y_xz, g_zz_xz_y_yy, g_zz_xz_y_yz, g_zz_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xz_y_xx[i] = -g_0_xz_y_xx[i] + 2.0 * g_zz_xz_y_xx[i] * a_exp;

        g_z_0_0_0_z_xz_y_xy[i] = -g_0_xz_y_xy[i] + 2.0 * g_zz_xz_y_xy[i] * a_exp;

        g_z_0_0_0_z_xz_y_xz[i] = -g_0_xz_y_xz[i] + 2.0 * g_zz_xz_y_xz[i] * a_exp;

        g_z_0_0_0_z_xz_y_yy[i] = -g_0_xz_y_yy[i] + 2.0 * g_zz_xz_y_yy[i] * a_exp;

        g_z_0_0_0_z_xz_y_yz[i] = -g_0_xz_y_yz[i] + 2.0 * g_zz_xz_y_yz[i] * a_exp;

        g_z_0_0_0_z_xz_y_zz[i] = -g_0_xz_y_zz[i] + 2.0 * g_zz_xz_y_zz[i] * a_exp;
    }
    // integrals block (912-918)

    #pragma omp simd aligned(g_0_xz_z_xx, g_0_xz_z_xy, g_0_xz_z_xz, g_0_xz_z_yy, g_0_xz_z_yz, g_0_xz_z_zz, g_z_0_0_0_z_xz_z_xx, g_z_0_0_0_z_xz_z_xy, g_z_0_0_0_z_xz_z_xz, g_z_0_0_0_z_xz_z_yy, g_z_0_0_0_z_xz_z_yz, g_z_0_0_0_z_xz_z_zz, g_zz_xz_z_xx, g_zz_xz_z_xy, g_zz_xz_z_xz, g_zz_xz_z_yy, g_zz_xz_z_yz, g_zz_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xz_z_xx[i] = -g_0_xz_z_xx[i] + 2.0 * g_zz_xz_z_xx[i] * a_exp;

        g_z_0_0_0_z_xz_z_xy[i] = -g_0_xz_z_xy[i] + 2.0 * g_zz_xz_z_xy[i] * a_exp;

        g_z_0_0_0_z_xz_z_xz[i] = -g_0_xz_z_xz[i] + 2.0 * g_zz_xz_z_xz[i] * a_exp;

        g_z_0_0_0_z_xz_z_yy[i] = -g_0_xz_z_yy[i] + 2.0 * g_zz_xz_z_yy[i] * a_exp;

        g_z_0_0_0_z_xz_z_yz[i] = -g_0_xz_z_yz[i] + 2.0 * g_zz_xz_z_yz[i] * a_exp;

        g_z_0_0_0_z_xz_z_zz[i] = -g_0_xz_z_zz[i] + 2.0 * g_zz_xz_z_zz[i] * a_exp;
    }
    // integrals block (918-924)

    #pragma omp simd aligned(g_0_yy_x_xx, g_0_yy_x_xy, g_0_yy_x_xz, g_0_yy_x_yy, g_0_yy_x_yz, g_0_yy_x_zz, g_z_0_0_0_z_yy_x_xx, g_z_0_0_0_z_yy_x_xy, g_z_0_0_0_z_yy_x_xz, g_z_0_0_0_z_yy_x_yy, g_z_0_0_0_z_yy_x_yz, g_z_0_0_0_z_yy_x_zz, g_zz_yy_x_xx, g_zz_yy_x_xy, g_zz_yy_x_xz, g_zz_yy_x_yy, g_zz_yy_x_yz, g_zz_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_yy_x_xx[i] = -g_0_yy_x_xx[i] + 2.0 * g_zz_yy_x_xx[i] * a_exp;

        g_z_0_0_0_z_yy_x_xy[i] = -g_0_yy_x_xy[i] + 2.0 * g_zz_yy_x_xy[i] * a_exp;

        g_z_0_0_0_z_yy_x_xz[i] = -g_0_yy_x_xz[i] + 2.0 * g_zz_yy_x_xz[i] * a_exp;

        g_z_0_0_0_z_yy_x_yy[i] = -g_0_yy_x_yy[i] + 2.0 * g_zz_yy_x_yy[i] * a_exp;

        g_z_0_0_0_z_yy_x_yz[i] = -g_0_yy_x_yz[i] + 2.0 * g_zz_yy_x_yz[i] * a_exp;

        g_z_0_0_0_z_yy_x_zz[i] = -g_0_yy_x_zz[i] + 2.0 * g_zz_yy_x_zz[i] * a_exp;
    }
    // integrals block (924-930)

    #pragma omp simd aligned(g_0_yy_y_xx, g_0_yy_y_xy, g_0_yy_y_xz, g_0_yy_y_yy, g_0_yy_y_yz, g_0_yy_y_zz, g_z_0_0_0_z_yy_y_xx, g_z_0_0_0_z_yy_y_xy, g_z_0_0_0_z_yy_y_xz, g_z_0_0_0_z_yy_y_yy, g_z_0_0_0_z_yy_y_yz, g_z_0_0_0_z_yy_y_zz, g_zz_yy_y_xx, g_zz_yy_y_xy, g_zz_yy_y_xz, g_zz_yy_y_yy, g_zz_yy_y_yz, g_zz_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_yy_y_xx[i] = -g_0_yy_y_xx[i] + 2.0 * g_zz_yy_y_xx[i] * a_exp;

        g_z_0_0_0_z_yy_y_xy[i] = -g_0_yy_y_xy[i] + 2.0 * g_zz_yy_y_xy[i] * a_exp;

        g_z_0_0_0_z_yy_y_xz[i] = -g_0_yy_y_xz[i] + 2.0 * g_zz_yy_y_xz[i] * a_exp;

        g_z_0_0_0_z_yy_y_yy[i] = -g_0_yy_y_yy[i] + 2.0 * g_zz_yy_y_yy[i] * a_exp;

        g_z_0_0_0_z_yy_y_yz[i] = -g_0_yy_y_yz[i] + 2.0 * g_zz_yy_y_yz[i] * a_exp;

        g_z_0_0_0_z_yy_y_zz[i] = -g_0_yy_y_zz[i] + 2.0 * g_zz_yy_y_zz[i] * a_exp;
    }
    // integrals block (930-936)

    #pragma omp simd aligned(g_0_yy_z_xx, g_0_yy_z_xy, g_0_yy_z_xz, g_0_yy_z_yy, g_0_yy_z_yz, g_0_yy_z_zz, g_z_0_0_0_z_yy_z_xx, g_z_0_0_0_z_yy_z_xy, g_z_0_0_0_z_yy_z_xz, g_z_0_0_0_z_yy_z_yy, g_z_0_0_0_z_yy_z_yz, g_z_0_0_0_z_yy_z_zz, g_zz_yy_z_xx, g_zz_yy_z_xy, g_zz_yy_z_xz, g_zz_yy_z_yy, g_zz_yy_z_yz, g_zz_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_yy_z_xx[i] = -g_0_yy_z_xx[i] + 2.0 * g_zz_yy_z_xx[i] * a_exp;

        g_z_0_0_0_z_yy_z_xy[i] = -g_0_yy_z_xy[i] + 2.0 * g_zz_yy_z_xy[i] * a_exp;

        g_z_0_0_0_z_yy_z_xz[i] = -g_0_yy_z_xz[i] + 2.0 * g_zz_yy_z_xz[i] * a_exp;

        g_z_0_0_0_z_yy_z_yy[i] = -g_0_yy_z_yy[i] + 2.0 * g_zz_yy_z_yy[i] * a_exp;

        g_z_0_0_0_z_yy_z_yz[i] = -g_0_yy_z_yz[i] + 2.0 * g_zz_yy_z_yz[i] * a_exp;

        g_z_0_0_0_z_yy_z_zz[i] = -g_0_yy_z_zz[i] + 2.0 * g_zz_yy_z_zz[i] * a_exp;
    }
    // integrals block (936-942)

    #pragma omp simd aligned(g_0_yz_x_xx, g_0_yz_x_xy, g_0_yz_x_xz, g_0_yz_x_yy, g_0_yz_x_yz, g_0_yz_x_zz, g_z_0_0_0_z_yz_x_xx, g_z_0_0_0_z_yz_x_xy, g_z_0_0_0_z_yz_x_xz, g_z_0_0_0_z_yz_x_yy, g_z_0_0_0_z_yz_x_yz, g_z_0_0_0_z_yz_x_zz, g_zz_yz_x_xx, g_zz_yz_x_xy, g_zz_yz_x_xz, g_zz_yz_x_yy, g_zz_yz_x_yz, g_zz_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_yz_x_xx[i] = -g_0_yz_x_xx[i] + 2.0 * g_zz_yz_x_xx[i] * a_exp;

        g_z_0_0_0_z_yz_x_xy[i] = -g_0_yz_x_xy[i] + 2.0 * g_zz_yz_x_xy[i] * a_exp;

        g_z_0_0_0_z_yz_x_xz[i] = -g_0_yz_x_xz[i] + 2.0 * g_zz_yz_x_xz[i] * a_exp;

        g_z_0_0_0_z_yz_x_yy[i] = -g_0_yz_x_yy[i] + 2.0 * g_zz_yz_x_yy[i] * a_exp;

        g_z_0_0_0_z_yz_x_yz[i] = -g_0_yz_x_yz[i] + 2.0 * g_zz_yz_x_yz[i] * a_exp;

        g_z_0_0_0_z_yz_x_zz[i] = -g_0_yz_x_zz[i] + 2.0 * g_zz_yz_x_zz[i] * a_exp;
    }
    // integrals block (942-948)

    #pragma omp simd aligned(g_0_yz_y_xx, g_0_yz_y_xy, g_0_yz_y_xz, g_0_yz_y_yy, g_0_yz_y_yz, g_0_yz_y_zz, g_z_0_0_0_z_yz_y_xx, g_z_0_0_0_z_yz_y_xy, g_z_0_0_0_z_yz_y_xz, g_z_0_0_0_z_yz_y_yy, g_z_0_0_0_z_yz_y_yz, g_z_0_0_0_z_yz_y_zz, g_zz_yz_y_xx, g_zz_yz_y_xy, g_zz_yz_y_xz, g_zz_yz_y_yy, g_zz_yz_y_yz, g_zz_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_yz_y_xx[i] = -g_0_yz_y_xx[i] + 2.0 * g_zz_yz_y_xx[i] * a_exp;

        g_z_0_0_0_z_yz_y_xy[i] = -g_0_yz_y_xy[i] + 2.0 * g_zz_yz_y_xy[i] * a_exp;

        g_z_0_0_0_z_yz_y_xz[i] = -g_0_yz_y_xz[i] + 2.0 * g_zz_yz_y_xz[i] * a_exp;

        g_z_0_0_0_z_yz_y_yy[i] = -g_0_yz_y_yy[i] + 2.0 * g_zz_yz_y_yy[i] * a_exp;

        g_z_0_0_0_z_yz_y_yz[i] = -g_0_yz_y_yz[i] + 2.0 * g_zz_yz_y_yz[i] * a_exp;

        g_z_0_0_0_z_yz_y_zz[i] = -g_0_yz_y_zz[i] + 2.0 * g_zz_yz_y_zz[i] * a_exp;
    }
    // integrals block (948-954)

    #pragma omp simd aligned(g_0_yz_z_xx, g_0_yz_z_xy, g_0_yz_z_xz, g_0_yz_z_yy, g_0_yz_z_yz, g_0_yz_z_zz, g_z_0_0_0_z_yz_z_xx, g_z_0_0_0_z_yz_z_xy, g_z_0_0_0_z_yz_z_xz, g_z_0_0_0_z_yz_z_yy, g_z_0_0_0_z_yz_z_yz, g_z_0_0_0_z_yz_z_zz, g_zz_yz_z_xx, g_zz_yz_z_xy, g_zz_yz_z_xz, g_zz_yz_z_yy, g_zz_yz_z_yz, g_zz_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_yz_z_xx[i] = -g_0_yz_z_xx[i] + 2.0 * g_zz_yz_z_xx[i] * a_exp;

        g_z_0_0_0_z_yz_z_xy[i] = -g_0_yz_z_xy[i] + 2.0 * g_zz_yz_z_xy[i] * a_exp;

        g_z_0_0_0_z_yz_z_xz[i] = -g_0_yz_z_xz[i] + 2.0 * g_zz_yz_z_xz[i] * a_exp;

        g_z_0_0_0_z_yz_z_yy[i] = -g_0_yz_z_yy[i] + 2.0 * g_zz_yz_z_yy[i] * a_exp;

        g_z_0_0_0_z_yz_z_yz[i] = -g_0_yz_z_yz[i] + 2.0 * g_zz_yz_z_yz[i] * a_exp;

        g_z_0_0_0_z_yz_z_zz[i] = -g_0_yz_z_zz[i] + 2.0 * g_zz_yz_z_zz[i] * a_exp;
    }
    // integrals block (954-960)

    #pragma omp simd aligned(g_0_zz_x_xx, g_0_zz_x_xy, g_0_zz_x_xz, g_0_zz_x_yy, g_0_zz_x_yz, g_0_zz_x_zz, g_z_0_0_0_z_zz_x_xx, g_z_0_0_0_z_zz_x_xy, g_z_0_0_0_z_zz_x_xz, g_z_0_0_0_z_zz_x_yy, g_z_0_0_0_z_zz_x_yz, g_z_0_0_0_z_zz_x_zz, g_zz_zz_x_xx, g_zz_zz_x_xy, g_zz_zz_x_xz, g_zz_zz_x_yy, g_zz_zz_x_yz, g_zz_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_zz_x_xx[i] = -g_0_zz_x_xx[i] + 2.0 * g_zz_zz_x_xx[i] * a_exp;

        g_z_0_0_0_z_zz_x_xy[i] = -g_0_zz_x_xy[i] + 2.0 * g_zz_zz_x_xy[i] * a_exp;

        g_z_0_0_0_z_zz_x_xz[i] = -g_0_zz_x_xz[i] + 2.0 * g_zz_zz_x_xz[i] * a_exp;

        g_z_0_0_0_z_zz_x_yy[i] = -g_0_zz_x_yy[i] + 2.0 * g_zz_zz_x_yy[i] * a_exp;

        g_z_0_0_0_z_zz_x_yz[i] = -g_0_zz_x_yz[i] + 2.0 * g_zz_zz_x_yz[i] * a_exp;

        g_z_0_0_0_z_zz_x_zz[i] = -g_0_zz_x_zz[i] + 2.0 * g_zz_zz_x_zz[i] * a_exp;
    }
    // integrals block (960-966)

    #pragma omp simd aligned(g_0_zz_y_xx, g_0_zz_y_xy, g_0_zz_y_xz, g_0_zz_y_yy, g_0_zz_y_yz, g_0_zz_y_zz, g_z_0_0_0_z_zz_y_xx, g_z_0_0_0_z_zz_y_xy, g_z_0_0_0_z_zz_y_xz, g_z_0_0_0_z_zz_y_yy, g_z_0_0_0_z_zz_y_yz, g_z_0_0_0_z_zz_y_zz, g_zz_zz_y_xx, g_zz_zz_y_xy, g_zz_zz_y_xz, g_zz_zz_y_yy, g_zz_zz_y_yz, g_zz_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_zz_y_xx[i] = -g_0_zz_y_xx[i] + 2.0 * g_zz_zz_y_xx[i] * a_exp;

        g_z_0_0_0_z_zz_y_xy[i] = -g_0_zz_y_xy[i] + 2.0 * g_zz_zz_y_xy[i] * a_exp;

        g_z_0_0_0_z_zz_y_xz[i] = -g_0_zz_y_xz[i] + 2.0 * g_zz_zz_y_xz[i] * a_exp;

        g_z_0_0_0_z_zz_y_yy[i] = -g_0_zz_y_yy[i] + 2.0 * g_zz_zz_y_yy[i] * a_exp;

        g_z_0_0_0_z_zz_y_yz[i] = -g_0_zz_y_yz[i] + 2.0 * g_zz_zz_y_yz[i] * a_exp;

        g_z_0_0_0_z_zz_y_zz[i] = -g_0_zz_y_zz[i] + 2.0 * g_zz_zz_y_zz[i] * a_exp;
    }
    // integrals block (966-972)

    #pragma omp simd aligned(g_0_zz_z_xx, g_0_zz_z_xy, g_0_zz_z_xz, g_0_zz_z_yy, g_0_zz_z_yz, g_0_zz_z_zz, g_z_0_0_0_z_zz_z_xx, g_z_0_0_0_z_zz_z_xy, g_z_0_0_0_z_zz_z_xz, g_z_0_0_0_z_zz_z_yy, g_z_0_0_0_z_zz_z_yz, g_z_0_0_0_z_zz_z_zz, g_zz_zz_z_xx, g_zz_zz_z_xy, g_zz_zz_z_xz, g_zz_zz_z_yy, g_zz_zz_z_yz, g_zz_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_zz_z_xx[i] = -g_0_zz_z_xx[i] + 2.0 * g_zz_zz_z_xx[i] * a_exp;

        g_z_0_0_0_z_zz_z_xy[i] = -g_0_zz_z_xy[i] + 2.0 * g_zz_zz_z_xy[i] * a_exp;

        g_z_0_0_0_z_zz_z_xz[i] = -g_0_zz_z_xz[i] + 2.0 * g_zz_zz_z_xz[i] * a_exp;

        g_z_0_0_0_z_zz_z_yy[i] = -g_0_zz_z_yy[i] + 2.0 * g_zz_zz_z_yy[i] * a_exp;

        g_z_0_0_0_z_zz_z_yz[i] = -g_0_zz_z_yz[i] + 2.0 * g_zz_zz_z_yz[i] * a_exp;

        g_z_0_0_0_z_zz_z_zz[i] = -g_0_zz_z_zz[i] + 2.0 * g_zz_zz_z_zz[i] * a_exp;
    }
}

} // t4c_geom namespace


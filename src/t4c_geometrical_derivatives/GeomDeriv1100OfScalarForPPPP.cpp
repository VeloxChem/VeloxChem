#include "GeomDeriv1100OfScalarForPPPP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1100_pppp_0(CSimdArray<double>& buffer_1100_pppp,
                     const CSimdArray<double>& buffer_sspp,
                     const CSimdArray<double>& buffer_sdpp,
                     const CSimdArray<double>& buffer_dspp,
                     const CSimdArray<double>& buffer_ddpp,
                     const double a_exp,
                     const double b_exp) -> void
{
    const auto ndims = buffer_1100_pppp.number_of_columns();

    /// Set up components of auxilary buffer : buffer_sspp

    auto g_0_0_x_x = buffer_sspp[0];

    auto g_0_0_x_y = buffer_sspp[1];

    auto g_0_0_x_z = buffer_sspp[2];

    auto g_0_0_y_x = buffer_sspp[3];

    auto g_0_0_y_y = buffer_sspp[4];

    auto g_0_0_y_z = buffer_sspp[5];

    auto g_0_0_z_x = buffer_sspp[6];

    auto g_0_0_z_y = buffer_sspp[7];

    auto g_0_0_z_z = buffer_sspp[8];

    /// Set up components of auxilary buffer : buffer_sdpp

    auto g_0_xx_x_x = buffer_sdpp[0];

    auto g_0_xx_x_y = buffer_sdpp[1];

    auto g_0_xx_x_z = buffer_sdpp[2];

    auto g_0_xx_y_x = buffer_sdpp[3];

    auto g_0_xx_y_y = buffer_sdpp[4];

    auto g_0_xx_y_z = buffer_sdpp[5];

    auto g_0_xx_z_x = buffer_sdpp[6];

    auto g_0_xx_z_y = buffer_sdpp[7];

    auto g_0_xx_z_z = buffer_sdpp[8];

    auto g_0_xy_x_x = buffer_sdpp[9];

    auto g_0_xy_x_y = buffer_sdpp[10];

    auto g_0_xy_x_z = buffer_sdpp[11];

    auto g_0_xy_y_x = buffer_sdpp[12];

    auto g_0_xy_y_y = buffer_sdpp[13];

    auto g_0_xy_y_z = buffer_sdpp[14];

    auto g_0_xy_z_x = buffer_sdpp[15];

    auto g_0_xy_z_y = buffer_sdpp[16];

    auto g_0_xy_z_z = buffer_sdpp[17];

    auto g_0_xz_x_x = buffer_sdpp[18];

    auto g_0_xz_x_y = buffer_sdpp[19];

    auto g_0_xz_x_z = buffer_sdpp[20];

    auto g_0_xz_y_x = buffer_sdpp[21];

    auto g_0_xz_y_y = buffer_sdpp[22];

    auto g_0_xz_y_z = buffer_sdpp[23];

    auto g_0_xz_z_x = buffer_sdpp[24];

    auto g_0_xz_z_y = buffer_sdpp[25];

    auto g_0_xz_z_z = buffer_sdpp[26];

    auto g_0_yy_x_x = buffer_sdpp[27];

    auto g_0_yy_x_y = buffer_sdpp[28];

    auto g_0_yy_x_z = buffer_sdpp[29];

    auto g_0_yy_y_x = buffer_sdpp[30];

    auto g_0_yy_y_y = buffer_sdpp[31];

    auto g_0_yy_y_z = buffer_sdpp[32];

    auto g_0_yy_z_x = buffer_sdpp[33];

    auto g_0_yy_z_y = buffer_sdpp[34];

    auto g_0_yy_z_z = buffer_sdpp[35];

    auto g_0_yz_x_x = buffer_sdpp[36];

    auto g_0_yz_x_y = buffer_sdpp[37];

    auto g_0_yz_x_z = buffer_sdpp[38];

    auto g_0_yz_y_x = buffer_sdpp[39];

    auto g_0_yz_y_y = buffer_sdpp[40];

    auto g_0_yz_y_z = buffer_sdpp[41];

    auto g_0_yz_z_x = buffer_sdpp[42];

    auto g_0_yz_z_y = buffer_sdpp[43];

    auto g_0_yz_z_z = buffer_sdpp[44];

    auto g_0_zz_x_x = buffer_sdpp[45];

    auto g_0_zz_x_y = buffer_sdpp[46];

    auto g_0_zz_x_z = buffer_sdpp[47];

    auto g_0_zz_y_x = buffer_sdpp[48];

    auto g_0_zz_y_y = buffer_sdpp[49];

    auto g_0_zz_y_z = buffer_sdpp[50];

    auto g_0_zz_z_x = buffer_sdpp[51];

    auto g_0_zz_z_y = buffer_sdpp[52];

    auto g_0_zz_z_z = buffer_sdpp[53];

    /// Set up components of auxilary buffer : buffer_dspp

    auto g_xx_0_x_x = buffer_dspp[0];

    auto g_xx_0_x_y = buffer_dspp[1];

    auto g_xx_0_x_z = buffer_dspp[2];

    auto g_xx_0_y_x = buffer_dspp[3];

    auto g_xx_0_y_y = buffer_dspp[4];

    auto g_xx_0_y_z = buffer_dspp[5];

    auto g_xx_0_z_x = buffer_dspp[6];

    auto g_xx_0_z_y = buffer_dspp[7];

    auto g_xx_0_z_z = buffer_dspp[8];

    auto g_xy_0_x_x = buffer_dspp[9];

    auto g_xy_0_x_y = buffer_dspp[10];

    auto g_xy_0_x_z = buffer_dspp[11];

    auto g_xy_0_y_x = buffer_dspp[12];

    auto g_xy_0_y_y = buffer_dspp[13];

    auto g_xy_0_y_z = buffer_dspp[14];

    auto g_xy_0_z_x = buffer_dspp[15];

    auto g_xy_0_z_y = buffer_dspp[16];

    auto g_xy_0_z_z = buffer_dspp[17];

    auto g_xz_0_x_x = buffer_dspp[18];

    auto g_xz_0_x_y = buffer_dspp[19];

    auto g_xz_0_x_z = buffer_dspp[20];

    auto g_xz_0_y_x = buffer_dspp[21];

    auto g_xz_0_y_y = buffer_dspp[22];

    auto g_xz_0_y_z = buffer_dspp[23];

    auto g_xz_0_z_x = buffer_dspp[24];

    auto g_xz_0_z_y = buffer_dspp[25];

    auto g_xz_0_z_z = buffer_dspp[26];

    auto g_yy_0_x_x = buffer_dspp[27];

    auto g_yy_0_x_y = buffer_dspp[28];

    auto g_yy_0_x_z = buffer_dspp[29];

    auto g_yy_0_y_x = buffer_dspp[30];

    auto g_yy_0_y_y = buffer_dspp[31];

    auto g_yy_0_y_z = buffer_dspp[32];

    auto g_yy_0_z_x = buffer_dspp[33];

    auto g_yy_0_z_y = buffer_dspp[34];

    auto g_yy_0_z_z = buffer_dspp[35];

    auto g_yz_0_x_x = buffer_dspp[36];

    auto g_yz_0_x_y = buffer_dspp[37];

    auto g_yz_0_x_z = buffer_dspp[38];

    auto g_yz_0_y_x = buffer_dspp[39];

    auto g_yz_0_y_y = buffer_dspp[40];

    auto g_yz_0_y_z = buffer_dspp[41];

    auto g_yz_0_z_x = buffer_dspp[42];

    auto g_yz_0_z_y = buffer_dspp[43];

    auto g_yz_0_z_z = buffer_dspp[44];

    auto g_zz_0_x_x = buffer_dspp[45];

    auto g_zz_0_x_y = buffer_dspp[46];

    auto g_zz_0_x_z = buffer_dspp[47];

    auto g_zz_0_y_x = buffer_dspp[48];

    auto g_zz_0_y_y = buffer_dspp[49];

    auto g_zz_0_y_z = buffer_dspp[50];

    auto g_zz_0_z_x = buffer_dspp[51];

    auto g_zz_0_z_y = buffer_dspp[52];

    auto g_zz_0_z_z = buffer_dspp[53];

    /// Set up components of auxilary buffer : buffer_ddpp

    auto g_xx_xx_x_x = buffer_ddpp[0];

    auto g_xx_xx_x_y = buffer_ddpp[1];

    auto g_xx_xx_x_z = buffer_ddpp[2];

    auto g_xx_xx_y_x = buffer_ddpp[3];

    auto g_xx_xx_y_y = buffer_ddpp[4];

    auto g_xx_xx_y_z = buffer_ddpp[5];

    auto g_xx_xx_z_x = buffer_ddpp[6];

    auto g_xx_xx_z_y = buffer_ddpp[7];

    auto g_xx_xx_z_z = buffer_ddpp[8];

    auto g_xx_xy_x_x = buffer_ddpp[9];

    auto g_xx_xy_x_y = buffer_ddpp[10];

    auto g_xx_xy_x_z = buffer_ddpp[11];

    auto g_xx_xy_y_x = buffer_ddpp[12];

    auto g_xx_xy_y_y = buffer_ddpp[13];

    auto g_xx_xy_y_z = buffer_ddpp[14];

    auto g_xx_xy_z_x = buffer_ddpp[15];

    auto g_xx_xy_z_y = buffer_ddpp[16];

    auto g_xx_xy_z_z = buffer_ddpp[17];

    auto g_xx_xz_x_x = buffer_ddpp[18];

    auto g_xx_xz_x_y = buffer_ddpp[19];

    auto g_xx_xz_x_z = buffer_ddpp[20];

    auto g_xx_xz_y_x = buffer_ddpp[21];

    auto g_xx_xz_y_y = buffer_ddpp[22];

    auto g_xx_xz_y_z = buffer_ddpp[23];

    auto g_xx_xz_z_x = buffer_ddpp[24];

    auto g_xx_xz_z_y = buffer_ddpp[25];

    auto g_xx_xz_z_z = buffer_ddpp[26];

    auto g_xx_yy_x_x = buffer_ddpp[27];

    auto g_xx_yy_x_y = buffer_ddpp[28];

    auto g_xx_yy_x_z = buffer_ddpp[29];

    auto g_xx_yy_y_x = buffer_ddpp[30];

    auto g_xx_yy_y_y = buffer_ddpp[31];

    auto g_xx_yy_y_z = buffer_ddpp[32];

    auto g_xx_yy_z_x = buffer_ddpp[33];

    auto g_xx_yy_z_y = buffer_ddpp[34];

    auto g_xx_yy_z_z = buffer_ddpp[35];

    auto g_xx_yz_x_x = buffer_ddpp[36];

    auto g_xx_yz_x_y = buffer_ddpp[37];

    auto g_xx_yz_x_z = buffer_ddpp[38];

    auto g_xx_yz_y_x = buffer_ddpp[39];

    auto g_xx_yz_y_y = buffer_ddpp[40];

    auto g_xx_yz_y_z = buffer_ddpp[41];

    auto g_xx_yz_z_x = buffer_ddpp[42];

    auto g_xx_yz_z_y = buffer_ddpp[43];

    auto g_xx_yz_z_z = buffer_ddpp[44];

    auto g_xx_zz_x_x = buffer_ddpp[45];

    auto g_xx_zz_x_y = buffer_ddpp[46];

    auto g_xx_zz_x_z = buffer_ddpp[47];

    auto g_xx_zz_y_x = buffer_ddpp[48];

    auto g_xx_zz_y_y = buffer_ddpp[49];

    auto g_xx_zz_y_z = buffer_ddpp[50];

    auto g_xx_zz_z_x = buffer_ddpp[51];

    auto g_xx_zz_z_y = buffer_ddpp[52];

    auto g_xx_zz_z_z = buffer_ddpp[53];

    auto g_xy_xx_x_x = buffer_ddpp[54];

    auto g_xy_xx_x_y = buffer_ddpp[55];

    auto g_xy_xx_x_z = buffer_ddpp[56];

    auto g_xy_xx_y_x = buffer_ddpp[57];

    auto g_xy_xx_y_y = buffer_ddpp[58];

    auto g_xy_xx_y_z = buffer_ddpp[59];

    auto g_xy_xx_z_x = buffer_ddpp[60];

    auto g_xy_xx_z_y = buffer_ddpp[61];

    auto g_xy_xx_z_z = buffer_ddpp[62];

    auto g_xy_xy_x_x = buffer_ddpp[63];

    auto g_xy_xy_x_y = buffer_ddpp[64];

    auto g_xy_xy_x_z = buffer_ddpp[65];

    auto g_xy_xy_y_x = buffer_ddpp[66];

    auto g_xy_xy_y_y = buffer_ddpp[67];

    auto g_xy_xy_y_z = buffer_ddpp[68];

    auto g_xy_xy_z_x = buffer_ddpp[69];

    auto g_xy_xy_z_y = buffer_ddpp[70];

    auto g_xy_xy_z_z = buffer_ddpp[71];

    auto g_xy_xz_x_x = buffer_ddpp[72];

    auto g_xy_xz_x_y = buffer_ddpp[73];

    auto g_xy_xz_x_z = buffer_ddpp[74];

    auto g_xy_xz_y_x = buffer_ddpp[75];

    auto g_xy_xz_y_y = buffer_ddpp[76];

    auto g_xy_xz_y_z = buffer_ddpp[77];

    auto g_xy_xz_z_x = buffer_ddpp[78];

    auto g_xy_xz_z_y = buffer_ddpp[79];

    auto g_xy_xz_z_z = buffer_ddpp[80];

    auto g_xy_yy_x_x = buffer_ddpp[81];

    auto g_xy_yy_x_y = buffer_ddpp[82];

    auto g_xy_yy_x_z = buffer_ddpp[83];

    auto g_xy_yy_y_x = buffer_ddpp[84];

    auto g_xy_yy_y_y = buffer_ddpp[85];

    auto g_xy_yy_y_z = buffer_ddpp[86];

    auto g_xy_yy_z_x = buffer_ddpp[87];

    auto g_xy_yy_z_y = buffer_ddpp[88];

    auto g_xy_yy_z_z = buffer_ddpp[89];

    auto g_xy_yz_x_x = buffer_ddpp[90];

    auto g_xy_yz_x_y = buffer_ddpp[91];

    auto g_xy_yz_x_z = buffer_ddpp[92];

    auto g_xy_yz_y_x = buffer_ddpp[93];

    auto g_xy_yz_y_y = buffer_ddpp[94];

    auto g_xy_yz_y_z = buffer_ddpp[95];

    auto g_xy_yz_z_x = buffer_ddpp[96];

    auto g_xy_yz_z_y = buffer_ddpp[97];

    auto g_xy_yz_z_z = buffer_ddpp[98];

    auto g_xy_zz_x_x = buffer_ddpp[99];

    auto g_xy_zz_x_y = buffer_ddpp[100];

    auto g_xy_zz_x_z = buffer_ddpp[101];

    auto g_xy_zz_y_x = buffer_ddpp[102];

    auto g_xy_zz_y_y = buffer_ddpp[103];

    auto g_xy_zz_y_z = buffer_ddpp[104];

    auto g_xy_zz_z_x = buffer_ddpp[105];

    auto g_xy_zz_z_y = buffer_ddpp[106];

    auto g_xy_zz_z_z = buffer_ddpp[107];

    auto g_xz_xx_x_x = buffer_ddpp[108];

    auto g_xz_xx_x_y = buffer_ddpp[109];

    auto g_xz_xx_x_z = buffer_ddpp[110];

    auto g_xz_xx_y_x = buffer_ddpp[111];

    auto g_xz_xx_y_y = buffer_ddpp[112];

    auto g_xz_xx_y_z = buffer_ddpp[113];

    auto g_xz_xx_z_x = buffer_ddpp[114];

    auto g_xz_xx_z_y = buffer_ddpp[115];

    auto g_xz_xx_z_z = buffer_ddpp[116];

    auto g_xz_xy_x_x = buffer_ddpp[117];

    auto g_xz_xy_x_y = buffer_ddpp[118];

    auto g_xz_xy_x_z = buffer_ddpp[119];

    auto g_xz_xy_y_x = buffer_ddpp[120];

    auto g_xz_xy_y_y = buffer_ddpp[121];

    auto g_xz_xy_y_z = buffer_ddpp[122];

    auto g_xz_xy_z_x = buffer_ddpp[123];

    auto g_xz_xy_z_y = buffer_ddpp[124];

    auto g_xz_xy_z_z = buffer_ddpp[125];

    auto g_xz_xz_x_x = buffer_ddpp[126];

    auto g_xz_xz_x_y = buffer_ddpp[127];

    auto g_xz_xz_x_z = buffer_ddpp[128];

    auto g_xz_xz_y_x = buffer_ddpp[129];

    auto g_xz_xz_y_y = buffer_ddpp[130];

    auto g_xz_xz_y_z = buffer_ddpp[131];

    auto g_xz_xz_z_x = buffer_ddpp[132];

    auto g_xz_xz_z_y = buffer_ddpp[133];

    auto g_xz_xz_z_z = buffer_ddpp[134];

    auto g_xz_yy_x_x = buffer_ddpp[135];

    auto g_xz_yy_x_y = buffer_ddpp[136];

    auto g_xz_yy_x_z = buffer_ddpp[137];

    auto g_xz_yy_y_x = buffer_ddpp[138];

    auto g_xz_yy_y_y = buffer_ddpp[139];

    auto g_xz_yy_y_z = buffer_ddpp[140];

    auto g_xz_yy_z_x = buffer_ddpp[141];

    auto g_xz_yy_z_y = buffer_ddpp[142];

    auto g_xz_yy_z_z = buffer_ddpp[143];

    auto g_xz_yz_x_x = buffer_ddpp[144];

    auto g_xz_yz_x_y = buffer_ddpp[145];

    auto g_xz_yz_x_z = buffer_ddpp[146];

    auto g_xz_yz_y_x = buffer_ddpp[147];

    auto g_xz_yz_y_y = buffer_ddpp[148];

    auto g_xz_yz_y_z = buffer_ddpp[149];

    auto g_xz_yz_z_x = buffer_ddpp[150];

    auto g_xz_yz_z_y = buffer_ddpp[151];

    auto g_xz_yz_z_z = buffer_ddpp[152];

    auto g_xz_zz_x_x = buffer_ddpp[153];

    auto g_xz_zz_x_y = buffer_ddpp[154];

    auto g_xz_zz_x_z = buffer_ddpp[155];

    auto g_xz_zz_y_x = buffer_ddpp[156];

    auto g_xz_zz_y_y = buffer_ddpp[157];

    auto g_xz_zz_y_z = buffer_ddpp[158];

    auto g_xz_zz_z_x = buffer_ddpp[159];

    auto g_xz_zz_z_y = buffer_ddpp[160];

    auto g_xz_zz_z_z = buffer_ddpp[161];

    auto g_yy_xx_x_x = buffer_ddpp[162];

    auto g_yy_xx_x_y = buffer_ddpp[163];

    auto g_yy_xx_x_z = buffer_ddpp[164];

    auto g_yy_xx_y_x = buffer_ddpp[165];

    auto g_yy_xx_y_y = buffer_ddpp[166];

    auto g_yy_xx_y_z = buffer_ddpp[167];

    auto g_yy_xx_z_x = buffer_ddpp[168];

    auto g_yy_xx_z_y = buffer_ddpp[169];

    auto g_yy_xx_z_z = buffer_ddpp[170];

    auto g_yy_xy_x_x = buffer_ddpp[171];

    auto g_yy_xy_x_y = buffer_ddpp[172];

    auto g_yy_xy_x_z = buffer_ddpp[173];

    auto g_yy_xy_y_x = buffer_ddpp[174];

    auto g_yy_xy_y_y = buffer_ddpp[175];

    auto g_yy_xy_y_z = buffer_ddpp[176];

    auto g_yy_xy_z_x = buffer_ddpp[177];

    auto g_yy_xy_z_y = buffer_ddpp[178];

    auto g_yy_xy_z_z = buffer_ddpp[179];

    auto g_yy_xz_x_x = buffer_ddpp[180];

    auto g_yy_xz_x_y = buffer_ddpp[181];

    auto g_yy_xz_x_z = buffer_ddpp[182];

    auto g_yy_xz_y_x = buffer_ddpp[183];

    auto g_yy_xz_y_y = buffer_ddpp[184];

    auto g_yy_xz_y_z = buffer_ddpp[185];

    auto g_yy_xz_z_x = buffer_ddpp[186];

    auto g_yy_xz_z_y = buffer_ddpp[187];

    auto g_yy_xz_z_z = buffer_ddpp[188];

    auto g_yy_yy_x_x = buffer_ddpp[189];

    auto g_yy_yy_x_y = buffer_ddpp[190];

    auto g_yy_yy_x_z = buffer_ddpp[191];

    auto g_yy_yy_y_x = buffer_ddpp[192];

    auto g_yy_yy_y_y = buffer_ddpp[193];

    auto g_yy_yy_y_z = buffer_ddpp[194];

    auto g_yy_yy_z_x = buffer_ddpp[195];

    auto g_yy_yy_z_y = buffer_ddpp[196];

    auto g_yy_yy_z_z = buffer_ddpp[197];

    auto g_yy_yz_x_x = buffer_ddpp[198];

    auto g_yy_yz_x_y = buffer_ddpp[199];

    auto g_yy_yz_x_z = buffer_ddpp[200];

    auto g_yy_yz_y_x = buffer_ddpp[201];

    auto g_yy_yz_y_y = buffer_ddpp[202];

    auto g_yy_yz_y_z = buffer_ddpp[203];

    auto g_yy_yz_z_x = buffer_ddpp[204];

    auto g_yy_yz_z_y = buffer_ddpp[205];

    auto g_yy_yz_z_z = buffer_ddpp[206];

    auto g_yy_zz_x_x = buffer_ddpp[207];

    auto g_yy_zz_x_y = buffer_ddpp[208];

    auto g_yy_zz_x_z = buffer_ddpp[209];

    auto g_yy_zz_y_x = buffer_ddpp[210];

    auto g_yy_zz_y_y = buffer_ddpp[211];

    auto g_yy_zz_y_z = buffer_ddpp[212];

    auto g_yy_zz_z_x = buffer_ddpp[213];

    auto g_yy_zz_z_y = buffer_ddpp[214];

    auto g_yy_zz_z_z = buffer_ddpp[215];

    auto g_yz_xx_x_x = buffer_ddpp[216];

    auto g_yz_xx_x_y = buffer_ddpp[217];

    auto g_yz_xx_x_z = buffer_ddpp[218];

    auto g_yz_xx_y_x = buffer_ddpp[219];

    auto g_yz_xx_y_y = buffer_ddpp[220];

    auto g_yz_xx_y_z = buffer_ddpp[221];

    auto g_yz_xx_z_x = buffer_ddpp[222];

    auto g_yz_xx_z_y = buffer_ddpp[223];

    auto g_yz_xx_z_z = buffer_ddpp[224];

    auto g_yz_xy_x_x = buffer_ddpp[225];

    auto g_yz_xy_x_y = buffer_ddpp[226];

    auto g_yz_xy_x_z = buffer_ddpp[227];

    auto g_yz_xy_y_x = buffer_ddpp[228];

    auto g_yz_xy_y_y = buffer_ddpp[229];

    auto g_yz_xy_y_z = buffer_ddpp[230];

    auto g_yz_xy_z_x = buffer_ddpp[231];

    auto g_yz_xy_z_y = buffer_ddpp[232];

    auto g_yz_xy_z_z = buffer_ddpp[233];

    auto g_yz_xz_x_x = buffer_ddpp[234];

    auto g_yz_xz_x_y = buffer_ddpp[235];

    auto g_yz_xz_x_z = buffer_ddpp[236];

    auto g_yz_xz_y_x = buffer_ddpp[237];

    auto g_yz_xz_y_y = buffer_ddpp[238];

    auto g_yz_xz_y_z = buffer_ddpp[239];

    auto g_yz_xz_z_x = buffer_ddpp[240];

    auto g_yz_xz_z_y = buffer_ddpp[241];

    auto g_yz_xz_z_z = buffer_ddpp[242];

    auto g_yz_yy_x_x = buffer_ddpp[243];

    auto g_yz_yy_x_y = buffer_ddpp[244];

    auto g_yz_yy_x_z = buffer_ddpp[245];

    auto g_yz_yy_y_x = buffer_ddpp[246];

    auto g_yz_yy_y_y = buffer_ddpp[247];

    auto g_yz_yy_y_z = buffer_ddpp[248];

    auto g_yz_yy_z_x = buffer_ddpp[249];

    auto g_yz_yy_z_y = buffer_ddpp[250];

    auto g_yz_yy_z_z = buffer_ddpp[251];

    auto g_yz_yz_x_x = buffer_ddpp[252];

    auto g_yz_yz_x_y = buffer_ddpp[253];

    auto g_yz_yz_x_z = buffer_ddpp[254];

    auto g_yz_yz_y_x = buffer_ddpp[255];

    auto g_yz_yz_y_y = buffer_ddpp[256];

    auto g_yz_yz_y_z = buffer_ddpp[257];

    auto g_yz_yz_z_x = buffer_ddpp[258];

    auto g_yz_yz_z_y = buffer_ddpp[259];

    auto g_yz_yz_z_z = buffer_ddpp[260];

    auto g_yz_zz_x_x = buffer_ddpp[261];

    auto g_yz_zz_x_y = buffer_ddpp[262];

    auto g_yz_zz_x_z = buffer_ddpp[263];

    auto g_yz_zz_y_x = buffer_ddpp[264];

    auto g_yz_zz_y_y = buffer_ddpp[265];

    auto g_yz_zz_y_z = buffer_ddpp[266];

    auto g_yz_zz_z_x = buffer_ddpp[267];

    auto g_yz_zz_z_y = buffer_ddpp[268];

    auto g_yz_zz_z_z = buffer_ddpp[269];

    auto g_zz_xx_x_x = buffer_ddpp[270];

    auto g_zz_xx_x_y = buffer_ddpp[271];

    auto g_zz_xx_x_z = buffer_ddpp[272];

    auto g_zz_xx_y_x = buffer_ddpp[273];

    auto g_zz_xx_y_y = buffer_ddpp[274];

    auto g_zz_xx_y_z = buffer_ddpp[275];

    auto g_zz_xx_z_x = buffer_ddpp[276];

    auto g_zz_xx_z_y = buffer_ddpp[277];

    auto g_zz_xx_z_z = buffer_ddpp[278];

    auto g_zz_xy_x_x = buffer_ddpp[279];

    auto g_zz_xy_x_y = buffer_ddpp[280];

    auto g_zz_xy_x_z = buffer_ddpp[281];

    auto g_zz_xy_y_x = buffer_ddpp[282];

    auto g_zz_xy_y_y = buffer_ddpp[283];

    auto g_zz_xy_y_z = buffer_ddpp[284];

    auto g_zz_xy_z_x = buffer_ddpp[285];

    auto g_zz_xy_z_y = buffer_ddpp[286];

    auto g_zz_xy_z_z = buffer_ddpp[287];

    auto g_zz_xz_x_x = buffer_ddpp[288];

    auto g_zz_xz_x_y = buffer_ddpp[289];

    auto g_zz_xz_x_z = buffer_ddpp[290];

    auto g_zz_xz_y_x = buffer_ddpp[291];

    auto g_zz_xz_y_y = buffer_ddpp[292];

    auto g_zz_xz_y_z = buffer_ddpp[293];

    auto g_zz_xz_z_x = buffer_ddpp[294];

    auto g_zz_xz_z_y = buffer_ddpp[295];

    auto g_zz_xz_z_z = buffer_ddpp[296];

    auto g_zz_yy_x_x = buffer_ddpp[297];

    auto g_zz_yy_x_y = buffer_ddpp[298];

    auto g_zz_yy_x_z = buffer_ddpp[299];

    auto g_zz_yy_y_x = buffer_ddpp[300];

    auto g_zz_yy_y_y = buffer_ddpp[301];

    auto g_zz_yy_y_z = buffer_ddpp[302];

    auto g_zz_yy_z_x = buffer_ddpp[303];

    auto g_zz_yy_z_y = buffer_ddpp[304];

    auto g_zz_yy_z_z = buffer_ddpp[305];

    auto g_zz_yz_x_x = buffer_ddpp[306];

    auto g_zz_yz_x_y = buffer_ddpp[307];

    auto g_zz_yz_x_z = buffer_ddpp[308];

    auto g_zz_yz_y_x = buffer_ddpp[309];

    auto g_zz_yz_y_y = buffer_ddpp[310];

    auto g_zz_yz_y_z = buffer_ddpp[311];

    auto g_zz_yz_z_x = buffer_ddpp[312];

    auto g_zz_yz_z_y = buffer_ddpp[313];

    auto g_zz_yz_z_z = buffer_ddpp[314];

    auto g_zz_zz_x_x = buffer_ddpp[315];

    auto g_zz_zz_x_y = buffer_ddpp[316];

    auto g_zz_zz_x_z = buffer_ddpp[317];

    auto g_zz_zz_y_x = buffer_ddpp[318];

    auto g_zz_zz_y_y = buffer_ddpp[319];

    auto g_zz_zz_y_z = buffer_ddpp[320];

    auto g_zz_zz_z_x = buffer_ddpp[321];

    auto g_zz_zz_z_y = buffer_ddpp[322];

    auto g_zz_zz_z_z = buffer_ddpp[323];

    /// Set up components of integrals buffer : buffer_1100_pppp

    auto g_x_x_0_0_x_x_x_x = buffer_1100_pppp[0];

    auto g_x_x_0_0_x_x_x_y = buffer_1100_pppp[1];

    auto g_x_x_0_0_x_x_x_z = buffer_1100_pppp[2];

    auto g_x_x_0_0_x_x_y_x = buffer_1100_pppp[3];

    auto g_x_x_0_0_x_x_y_y = buffer_1100_pppp[4];

    auto g_x_x_0_0_x_x_y_z = buffer_1100_pppp[5];

    auto g_x_x_0_0_x_x_z_x = buffer_1100_pppp[6];

    auto g_x_x_0_0_x_x_z_y = buffer_1100_pppp[7];

    auto g_x_x_0_0_x_x_z_z = buffer_1100_pppp[8];

    auto g_x_x_0_0_x_y_x_x = buffer_1100_pppp[9];

    auto g_x_x_0_0_x_y_x_y = buffer_1100_pppp[10];

    auto g_x_x_0_0_x_y_x_z = buffer_1100_pppp[11];

    auto g_x_x_0_0_x_y_y_x = buffer_1100_pppp[12];

    auto g_x_x_0_0_x_y_y_y = buffer_1100_pppp[13];

    auto g_x_x_0_0_x_y_y_z = buffer_1100_pppp[14];

    auto g_x_x_0_0_x_y_z_x = buffer_1100_pppp[15];

    auto g_x_x_0_0_x_y_z_y = buffer_1100_pppp[16];

    auto g_x_x_0_0_x_y_z_z = buffer_1100_pppp[17];

    auto g_x_x_0_0_x_z_x_x = buffer_1100_pppp[18];

    auto g_x_x_0_0_x_z_x_y = buffer_1100_pppp[19];

    auto g_x_x_0_0_x_z_x_z = buffer_1100_pppp[20];

    auto g_x_x_0_0_x_z_y_x = buffer_1100_pppp[21];

    auto g_x_x_0_0_x_z_y_y = buffer_1100_pppp[22];

    auto g_x_x_0_0_x_z_y_z = buffer_1100_pppp[23];

    auto g_x_x_0_0_x_z_z_x = buffer_1100_pppp[24];

    auto g_x_x_0_0_x_z_z_y = buffer_1100_pppp[25];

    auto g_x_x_0_0_x_z_z_z = buffer_1100_pppp[26];

    auto g_x_x_0_0_y_x_x_x = buffer_1100_pppp[27];

    auto g_x_x_0_0_y_x_x_y = buffer_1100_pppp[28];

    auto g_x_x_0_0_y_x_x_z = buffer_1100_pppp[29];

    auto g_x_x_0_0_y_x_y_x = buffer_1100_pppp[30];

    auto g_x_x_0_0_y_x_y_y = buffer_1100_pppp[31];

    auto g_x_x_0_0_y_x_y_z = buffer_1100_pppp[32];

    auto g_x_x_0_0_y_x_z_x = buffer_1100_pppp[33];

    auto g_x_x_0_0_y_x_z_y = buffer_1100_pppp[34];

    auto g_x_x_0_0_y_x_z_z = buffer_1100_pppp[35];

    auto g_x_x_0_0_y_y_x_x = buffer_1100_pppp[36];

    auto g_x_x_0_0_y_y_x_y = buffer_1100_pppp[37];

    auto g_x_x_0_0_y_y_x_z = buffer_1100_pppp[38];

    auto g_x_x_0_0_y_y_y_x = buffer_1100_pppp[39];

    auto g_x_x_0_0_y_y_y_y = buffer_1100_pppp[40];

    auto g_x_x_0_0_y_y_y_z = buffer_1100_pppp[41];

    auto g_x_x_0_0_y_y_z_x = buffer_1100_pppp[42];

    auto g_x_x_0_0_y_y_z_y = buffer_1100_pppp[43];

    auto g_x_x_0_0_y_y_z_z = buffer_1100_pppp[44];

    auto g_x_x_0_0_y_z_x_x = buffer_1100_pppp[45];

    auto g_x_x_0_0_y_z_x_y = buffer_1100_pppp[46];

    auto g_x_x_0_0_y_z_x_z = buffer_1100_pppp[47];

    auto g_x_x_0_0_y_z_y_x = buffer_1100_pppp[48];

    auto g_x_x_0_0_y_z_y_y = buffer_1100_pppp[49];

    auto g_x_x_0_0_y_z_y_z = buffer_1100_pppp[50];

    auto g_x_x_0_0_y_z_z_x = buffer_1100_pppp[51];

    auto g_x_x_0_0_y_z_z_y = buffer_1100_pppp[52];

    auto g_x_x_0_0_y_z_z_z = buffer_1100_pppp[53];

    auto g_x_x_0_0_z_x_x_x = buffer_1100_pppp[54];

    auto g_x_x_0_0_z_x_x_y = buffer_1100_pppp[55];

    auto g_x_x_0_0_z_x_x_z = buffer_1100_pppp[56];

    auto g_x_x_0_0_z_x_y_x = buffer_1100_pppp[57];

    auto g_x_x_0_0_z_x_y_y = buffer_1100_pppp[58];

    auto g_x_x_0_0_z_x_y_z = buffer_1100_pppp[59];

    auto g_x_x_0_0_z_x_z_x = buffer_1100_pppp[60];

    auto g_x_x_0_0_z_x_z_y = buffer_1100_pppp[61];

    auto g_x_x_0_0_z_x_z_z = buffer_1100_pppp[62];

    auto g_x_x_0_0_z_y_x_x = buffer_1100_pppp[63];

    auto g_x_x_0_0_z_y_x_y = buffer_1100_pppp[64];

    auto g_x_x_0_0_z_y_x_z = buffer_1100_pppp[65];

    auto g_x_x_0_0_z_y_y_x = buffer_1100_pppp[66];

    auto g_x_x_0_0_z_y_y_y = buffer_1100_pppp[67];

    auto g_x_x_0_0_z_y_y_z = buffer_1100_pppp[68];

    auto g_x_x_0_0_z_y_z_x = buffer_1100_pppp[69];

    auto g_x_x_0_0_z_y_z_y = buffer_1100_pppp[70];

    auto g_x_x_0_0_z_y_z_z = buffer_1100_pppp[71];

    auto g_x_x_0_0_z_z_x_x = buffer_1100_pppp[72];

    auto g_x_x_0_0_z_z_x_y = buffer_1100_pppp[73];

    auto g_x_x_0_0_z_z_x_z = buffer_1100_pppp[74];

    auto g_x_x_0_0_z_z_y_x = buffer_1100_pppp[75];

    auto g_x_x_0_0_z_z_y_y = buffer_1100_pppp[76];

    auto g_x_x_0_0_z_z_y_z = buffer_1100_pppp[77];

    auto g_x_x_0_0_z_z_z_x = buffer_1100_pppp[78];

    auto g_x_x_0_0_z_z_z_y = buffer_1100_pppp[79];

    auto g_x_x_0_0_z_z_z_z = buffer_1100_pppp[80];

    auto g_x_y_0_0_x_x_x_x = buffer_1100_pppp[81];

    auto g_x_y_0_0_x_x_x_y = buffer_1100_pppp[82];

    auto g_x_y_0_0_x_x_x_z = buffer_1100_pppp[83];

    auto g_x_y_0_0_x_x_y_x = buffer_1100_pppp[84];

    auto g_x_y_0_0_x_x_y_y = buffer_1100_pppp[85];

    auto g_x_y_0_0_x_x_y_z = buffer_1100_pppp[86];

    auto g_x_y_0_0_x_x_z_x = buffer_1100_pppp[87];

    auto g_x_y_0_0_x_x_z_y = buffer_1100_pppp[88];

    auto g_x_y_0_0_x_x_z_z = buffer_1100_pppp[89];

    auto g_x_y_0_0_x_y_x_x = buffer_1100_pppp[90];

    auto g_x_y_0_0_x_y_x_y = buffer_1100_pppp[91];

    auto g_x_y_0_0_x_y_x_z = buffer_1100_pppp[92];

    auto g_x_y_0_0_x_y_y_x = buffer_1100_pppp[93];

    auto g_x_y_0_0_x_y_y_y = buffer_1100_pppp[94];

    auto g_x_y_0_0_x_y_y_z = buffer_1100_pppp[95];

    auto g_x_y_0_0_x_y_z_x = buffer_1100_pppp[96];

    auto g_x_y_0_0_x_y_z_y = buffer_1100_pppp[97];

    auto g_x_y_0_0_x_y_z_z = buffer_1100_pppp[98];

    auto g_x_y_0_0_x_z_x_x = buffer_1100_pppp[99];

    auto g_x_y_0_0_x_z_x_y = buffer_1100_pppp[100];

    auto g_x_y_0_0_x_z_x_z = buffer_1100_pppp[101];

    auto g_x_y_0_0_x_z_y_x = buffer_1100_pppp[102];

    auto g_x_y_0_0_x_z_y_y = buffer_1100_pppp[103];

    auto g_x_y_0_0_x_z_y_z = buffer_1100_pppp[104];

    auto g_x_y_0_0_x_z_z_x = buffer_1100_pppp[105];

    auto g_x_y_0_0_x_z_z_y = buffer_1100_pppp[106];

    auto g_x_y_0_0_x_z_z_z = buffer_1100_pppp[107];

    auto g_x_y_0_0_y_x_x_x = buffer_1100_pppp[108];

    auto g_x_y_0_0_y_x_x_y = buffer_1100_pppp[109];

    auto g_x_y_0_0_y_x_x_z = buffer_1100_pppp[110];

    auto g_x_y_0_0_y_x_y_x = buffer_1100_pppp[111];

    auto g_x_y_0_0_y_x_y_y = buffer_1100_pppp[112];

    auto g_x_y_0_0_y_x_y_z = buffer_1100_pppp[113];

    auto g_x_y_0_0_y_x_z_x = buffer_1100_pppp[114];

    auto g_x_y_0_0_y_x_z_y = buffer_1100_pppp[115];

    auto g_x_y_0_0_y_x_z_z = buffer_1100_pppp[116];

    auto g_x_y_0_0_y_y_x_x = buffer_1100_pppp[117];

    auto g_x_y_0_0_y_y_x_y = buffer_1100_pppp[118];

    auto g_x_y_0_0_y_y_x_z = buffer_1100_pppp[119];

    auto g_x_y_0_0_y_y_y_x = buffer_1100_pppp[120];

    auto g_x_y_0_0_y_y_y_y = buffer_1100_pppp[121];

    auto g_x_y_0_0_y_y_y_z = buffer_1100_pppp[122];

    auto g_x_y_0_0_y_y_z_x = buffer_1100_pppp[123];

    auto g_x_y_0_0_y_y_z_y = buffer_1100_pppp[124];

    auto g_x_y_0_0_y_y_z_z = buffer_1100_pppp[125];

    auto g_x_y_0_0_y_z_x_x = buffer_1100_pppp[126];

    auto g_x_y_0_0_y_z_x_y = buffer_1100_pppp[127];

    auto g_x_y_0_0_y_z_x_z = buffer_1100_pppp[128];

    auto g_x_y_0_0_y_z_y_x = buffer_1100_pppp[129];

    auto g_x_y_0_0_y_z_y_y = buffer_1100_pppp[130];

    auto g_x_y_0_0_y_z_y_z = buffer_1100_pppp[131];

    auto g_x_y_0_0_y_z_z_x = buffer_1100_pppp[132];

    auto g_x_y_0_0_y_z_z_y = buffer_1100_pppp[133];

    auto g_x_y_0_0_y_z_z_z = buffer_1100_pppp[134];

    auto g_x_y_0_0_z_x_x_x = buffer_1100_pppp[135];

    auto g_x_y_0_0_z_x_x_y = buffer_1100_pppp[136];

    auto g_x_y_0_0_z_x_x_z = buffer_1100_pppp[137];

    auto g_x_y_0_0_z_x_y_x = buffer_1100_pppp[138];

    auto g_x_y_0_0_z_x_y_y = buffer_1100_pppp[139];

    auto g_x_y_0_0_z_x_y_z = buffer_1100_pppp[140];

    auto g_x_y_0_0_z_x_z_x = buffer_1100_pppp[141];

    auto g_x_y_0_0_z_x_z_y = buffer_1100_pppp[142];

    auto g_x_y_0_0_z_x_z_z = buffer_1100_pppp[143];

    auto g_x_y_0_0_z_y_x_x = buffer_1100_pppp[144];

    auto g_x_y_0_0_z_y_x_y = buffer_1100_pppp[145];

    auto g_x_y_0_0_z_y_x_z = buffer_1100_pppp[146];

    auto g_x_y_0_0_z_y_y_x = buffer_1100_pppp[147];

    auto g_x_y_0_0_z_y_y_y = buffer_1100_pppp[148];

    auto g_x_y_0_0_z_y_y_z = buffer_1100_pppp[149];

    auto g_x_y_0_0_z_y_z_x = buffer_1100_pppp[150];

    auto g_x_y_0_0_z_y_z_y = buffer_1100_pppp[151];

    auto g_x_y_0_0_z_y_z_z = buffer_1100_pppp[152];

    auto g_x_y_0_0_z_z_x_x = buffer_1100_pppp[153];

    auto g_x_y_0_0_z_z_x_y = buffer_1100_pppp[154];

    auto g_x_y_0_0_z_z_x_z = buffer_1100_pppp[155];

    auto g_x_y_0_0_z_z_y_x = buffer_1100_pppp[156];

    auto g_x_y_0_0_z_z_y_y = buffer_1100_pppp[157];

    auto g_x_y_0_0_z_z_y_z = buffer_1100_pppp[158];

    auto g_x_y_0_0_z_z_z_x = buffer_1100_pppp[159];

    auto g_x_y_0_0_z_z_z_y = buffer_1100_pppp[160];

    auto g_x_y_0_0_z_z_z_z = buffer_1100_pppp[161];

    auto g_x_z_0_0_x_x_x_x = buffer_1100_pppp[162];

    auto g_x_z_0_0_x_x_x_y = buffer_1100_pppp[163];

    auto g_x_z_0_0_x_x_x_z = buffer_1100_pppp[164];

    auto g_x_z_0_0_x_x_y_x = buffer_1100_pppp[165];

    auto g_x_z_0_0_x_x_y_y = buffer_1100_pppp[166];

    auto g_x_z_0_0_x_x_y_z = buffer_1100_pppp[167];

    auto g_x_z_0_0_x_x_z_x = buffer_1100_pppp[168];

    auto g_x_z_0_0_x_x_z_y = buffer_1100_pppp[169];

    auto g_x_z_0_0_x_x_z_z = buffer_1100_pppp[170];

    auto g_x_z_0_0_x_y_x_x = buffer_1100_pppp[171];

    auto g_x_z_0_0_x_y_x_y = buffer_1100_pppp[172];

    auto g_x_z_0_0_x_y_x_z = buffer_1100_pppp[173];

    auto g_x_z_0_0_x_y_y_x = buffer_1100_pppp[174];

    auto g_x_z_0_0_x_y_y_y = buffer_1100_pppp[175];

    auto g_x_z_0_0_x_y_y_z = buffer_1100_pppp[176];

    auto g_x_z_0_0_x_y_z_x = buffer_1100_pppp[177];

    auto g_x_z_0_0_x_y_z_y = buffer_1100_pppp[178];

    auto g_x_z_0_0_x_y_z_z = buffer_1100_pppp[179];

    auto g_x_z_0_0_x_z_x_x = buffer_1100_pppp[180];

    auto g_x_z_0_0_x_z_x_y = buffer_1100_pppp[181];

    auto g_x_z_0_0_x_z_x_z = buffer_1100_pppp[182];

    auto g_x_z_0_0_x_z_y_x = buffer_1100_pppp[183];

    auto g_x_z_0_0_x_z_y_y = buffer_1100_pppp[184];

    auto g_x_z_0_0_x_z_y_z = buffer_1100_pppp[185];

    auto g_x_z_0_0_x_z_z_x = buffer_1100_pppp[186];

    auto g_x_z_0_0_x_z_z_y = buffer_1100_pppp[187];

    auto g_x_z_0_0_x_z_z_z = buffer_1100_pppp[188];

    auto g_x_z_0_0_y_x_x_x = buffer_1100_pppp[189];

    auto g_x_z_0_0_y_x_x_y = buffer_1100_pppp[190];

    auto g_x_z_0_0_y_x_x_z = buffer_1100_pppp[191];

    auto g_x_z_0_0_y_x_y_x = buffer_1100_pppp[192];

    auto g_x_z_0_0_y_x_y_y = buffer_1100_pppp[193];

    auto g_x_z_0_0_y_x_y_z = buffer_1100_pppp[194];

    auto g_x_z_0_0_y_x_z_x = buffer_1100_pppp[195];

    auto g_x_z_0_0_y_x_z_y = buffer_1100_pppp[196];

    auto g_x_z_0_0_y_x_z_z = buffer_1100_pppp[197];

    auto g_x_z_0_0_y_y_x_x = buffer_1100_pppp[198];

    auto g_x_z_0_0_y_y_x_y = buffer_1100_pppp[199];

    auto g_x_z_0_0_y_y_x_z = buffer_1100_pppp[200];

    auto g_x_z_0_0_y_y_y_x = buffer_1100_pppp[201];

    auto g_x_z_0_0_y_y_y_y = buffer_1100_pppp[202];

    auto g_x_z_0_0_y_y_y_z = buffer_1100_pppp[203];

    auto g_x_z_0_0_y_y_z_x = buffer_1100_pppp[204];

    auto g_x_z_0_0_y_y_z_y = buffer_1100_pppp[205];

    auto g_x_z_0_0_y_y_z_z = buffer_1100_pppp[206];

    auto g_x_z_0_0_y_z_x_x = buffer_1100_pppp[207];

    auto g_x_z_0_0_y_z_x_y = buffer_1100_pppp[208];

    auto g_x_z_0_0_y_z_x_z = buffer_1100_pppp[209];

    auto g_x_z_0_0_y_z_y_x = buffer_1100_pppp[210];

    auto g_x_z_0_0_y_z_y_y = buffer_1100_pppp[211];

    auto g_x_z_0_0_y_z_y_z = buffer_1100_pppp[212];

    auto g_x_z_0_0_y_z_z_x = buffer_1100_pppp[213];

    auto g_x_z_0_0_y_z_z_y = buffer_1100_pppp[214];

    auto g_x_z_0_0_y_z_z_z = buffer_1100_pppp[215];

    auto g_x_z_0_0_z_x_x_x = buffer_1100_pppp[216];

    auto g_x_z_0_0_z_x_x_y = buffer_1100_pppp[217];

    auto g_x_z_0_0_z_x_x_z = buffer_1100_pppp[218];

    auto g_x_z_0_0_z_x_y_x = buffer_1100_pppp[219];

    auto g_x_z_0_0_z_x_y_y = buffer_1100_pppp[220];

    auto g_x_z_0_0_z_x_y_z = buffer_1100_pppp[221];

    auto g_x_z_0_0_z_x_z_x = buffer_1100_pppp[222];

    auto g_x_z_0_0_z_x_z_y = buffer_1100_pppp[223];

    auto g_x_z_0_0_z_x_z_z = buffer_1100_pppp[224];

    auto g_x_z_0_0_z_y_x_x = buffer_1100_pppp[225];

    auto g_x_z_0_0_z_y_x_y = buffer_1100_pppp[226];

    auto g_x_z_0_0_z_y_x_z = buffer_1100_pppp[227];

    auto g_x_z_0_0_z_y_y_x = buffer_1100_pppp[228];

    auto g_x_z_0_0_z_y_y_y = buffer_1100_pppp[229];

    auto g_x_z_0_0_z_y_y_z = buffer_1100_pppp[230];

    auto g_x_z_0_0_z_y_z_x = buffer_1100_pppp[231];

    auto g_x_z_0_0_z_y_z_y = buffer_1100_pppp[232];

    auto g_x_z_0_0_z_y_z_z = buffer_1100_pppp[233];

    auto g_x_z_0_0_z_z_x_x = buffer_1100_pppp[234];

    auto g_x_z_0_0_z_z_x_y = buffer_1100_pppp[235];

    auto g_x_z_0_0_z_z_x_z = buffer_1100_pppp[236];

    auto g_x_z_0_0_z_z_y_x = buffer_1100_pppp[237];

    auto g_x_z_0_0_z_z_y_y = buffer_1100_pppp[238];

    auto g_x_z_0_0_z_z_y_z = buffer_1100_pppp[239];

    auto g_x_z_0_0_z_z_z_x = buffer_1100_pppp[240];

    auto g_x_z_0_0_z_z_z_y = buffer_1100_pppp[241];

    auto g_x_z_0_0_z_z_z_z = buffer_1100_pppp[242];

    auto g_y_x_0_0_x_x_x_x = buffer_1100_pppp[243];

    auto g_y_x_0_0_x_x_x_y = buffer_1100_pppp[244];

    auto g_y_x_0_0_x_x_x_z = buffer_1100_pppp[245];

    auto g_y_x_0_0_x_x_y_x = buffer_1100_pppp[246];

    auto g_y_x_0_0_x_x_y_y = buffer_1100_pppp[247];

    auto g_y_x_0_0_x_x_y_z = buffer_1100_pppp[248];

    auto g_y_x_0_0_x_x_z_x = buffer_1100_pppp[249];

    auto g_y_x_0_0_x_x_z_y = buffer_1100_pppp[250];

    auto g_y_x_0_0_x_x_z_z = buffer_1100_pppp[251];

    auto g_y_x_0_0_x_y_x_x = buffer_1100_pppp[252];

    auto g_y_x_0_0_x_y_x_y = buffer_1100_pppp[253];

    auto g_y_x_0_0_x_y_x_z = buffer_1100_pppp[254];

    auto g_y_x_0_0_x_y_y_x = buffer_1100_pppp[255];

    auto g_y_x_0_0_x_y_y_y = buffer_1100_pppp[256];

    auto g_y_x_0_0_x_y_y_z = buffer_1100_pppp[257];

    auto g_y_x_0_0_x_y_z_x = buffer_1100_pppp[258];

    auto g_y_x_0_0_x_y_z_y = buffer_1100_pppp[259];

    auto g_y_x_0_0_x_y_z_z = buffer_1100_pppp[260];

    auto g_y_x_0_0_x_z_x_x = buffer_1100_pppp[261];

    auto g_y_x_0_0_x_z_x_y = buffer_1100_pppp[262];

    auto g_y_x_0_0_x_z_x_z = buffer_1100_pppp[263];

    auto g_y_x_0_0_x_z_y_x = buffer_1100_pppp[264];

    auto g_y_x_0_0_x_z_y_y = buffer_1100_pppp[265];

    auto g_y_x_0_0_x_z_y_z = buffer_1100_pppp[266];

    auto g_y_x_0_0_x_z_z_x = buffer_1100_pppp[267];

    auto g_y_x_0_0_x_z_z_y = buffer_1100_pppp[268];

    auto g_y_x_0_0_x_z_z_z = buffer_1100_pppp[269];

    auto g_y_x_0_0_y_x_x_x = buffer_1100_pppp[270];

    auto g_y_x_0_0_y_x_x_y = buffer_1100_pppp[271];

    auto g_y_x_0_0_y_x_x_z = buffer_1100_pppp[272];

    auto g_y_x_0_0_y_x_y_x = buffer_1100_pppp[273];

    auto g_y_x_0_0_y_x_y_y = buffer_1100_pppp[274];

    auto g_y_x_0_0_y_x_y_z = buffer_1100_pppp[275];

    auto g_y_x_0_0_y_x_z_x = buffer_1100_pppp[276];

    auto g_y_x_0_0_y_x_z_y = buffer_1100_pppp[277];

    auto g_y_x_0_0_y_x_z_z = buffer_1100_pppp[278];

    auto g_y_x_0_0_y_y_x_x = buffer_1100_pppp[279];

    auto g_y_x_0_0_y_y_x_y = buffer_1100_pppp[280];

    auto g_y_x_0_0_y_y_x_z = buffer_1100_pppp[281];

    auto g_y_x_0_0_y_y_y_x = buffer_1100_pppp[282];

    auto g_y_x_0_0_y_y_y_y = buffer_1100_pppp[283];

    auto g_y_x_0_0_y_y_y_z = buffer_1100_pppp[284];

    auto g_y_x_0_0_y_y_z_x = buffer_1100_pppp[285];

    auto g_y_x_0_0_y_y_z_y = buffer_1100_pppp[286];

    auto g_y_x_0_0_y_y_z_z = buffer_1100_pppp[287];

    auto g_y_x_0_0_y_z_x_x = buffer_1100_pppp[288];

    auto g_y_x_0_0_y_z_x_y = buffer_1100_pppp[289];

    auto g_y_x_0_0_y_z_x_z = buffer_1100_pppp[290];

    auto g_y_x_0_0_y_z_y_x = buffer_1100_pppp[291];

    auto g_y_x_0_0_y_z_y_y = buffer_1100_pppp[292];

    auto g_y_x_0_0_y_z_y_z = buffer_1100_pppp[293];

    auto g_y_x_0_0_y_z_z_x = buffer_1100_pppp[294];

    auto g_y_x_0_0_y_z_z_y = buffer_1100_pppp[295];

    auto g_y_x_0_0_y_z_z_z = buffer_1100_pppp[296];

    auto g_y_x_0_0_z_x_x_x = buffer_1100_pppp[297];

    auto g_y_x_0_0_z_x_x_y = buffer_1100_pppp[298];

    auto g_y_x_0_0_z_x_x_z = buffer_1100_pppp[299];

    auto g_y_x_0_0_z_x_y_x = buffer_1100_pppp[300];

    auto g_y_x_0_0_z_x_y_y = buffer_1100_pppp[301];

    auto g_y_x_0_0_z_x_y_z = buffer_1100_pppp[302];

    auto g_y_x_0_0_z_x_z_x = buffer_1100_pppp[303];

    auto g_y_x_0_0_z_x_z_y = buffer_1100_pppp[304];

    auto g_y_x_0_0_z_x_z_z = buffer_1100_pppp[305];

    auto g_y_x_0_0_z_y_x_x = buffer_1100_pppp[306];

    auto g_y_x_0_0_z_y_x_y = buffer_1100_pppp[307];

    auto g_y_x_0_0_z_y_x_z = buffer_1100_pppp[308];

    auto g_y_x_0_0_z_y_y_x = buffer_1100_pppp[309];

    auto g_y_x_0_0_z_y_y_y = buffer_1100_pppp[310];

    auto g_y_x_0_0_z_y_y_z = buffer_1100_pppp[311];

    auto g_y_x_0_0_z_y_z_x = buffer_1100_pppp[312];

    auto g_y_x_0_0_z_y_z_y = buffer_1100_pppp[313];

    auto g_y_x_0_0_z_y_z_z = buffer_1100_pppp[314];

    auto g_y_x_0_0_z_z_x_x = buffer_1100_pppp[315];

    auto g_y_x_0_0_z_z_x_y = buffer_1100_pppp[316];

    auto g_y_x_0_0_z_z_x_z = buffer_1100_pppp[317];

    auto g_y_x_0_0_z_z_y_x = buffer_1100_pppp[318];

    auto g_y_x_0_0_z_z_y_y = buffer_1100_pppp[319];

    auto g_y_x_0_0_z_z_y_z = buffer_1100_pppp[320];

    auto g_y_x_0_0_z_z_z_x = buffer_1100_pppp[321];

    auto g_y_x_0_0_z_z_z_y = buffer_1100_pppp[322];

    auto g_y_x_0_0_z_z_z_z = buffer_1100_pppp[323];

    auto g_y_y_0_0_x_x_x_x = buffer_1100_pppp[324];

    auto g_y_y_0_0_x_x_x_y = buffer_1100_pppp[325];

    auto g_y_y_0_0_x_x_x_z = buffer_1100_pppp[326];

    auto g_y_y_0_0_x_x_y_x = buffer_1100_pppp[327];

    auto g_y_y_0_0_x_x_y_y = buffer_1100_pppp[328];

    auto g_y_y_0_0_x_x_y_z = buffer_1100_pppp[329];

    auto g_y_y_0_0_x_x_z_x = buffer_1100_pppp[330];

    auto g_y_y_0_0_x_x_z_y = buffer_1100_pppp[331];

    auto g_y_y_0_0_x_x_z_z = buffer_1100_pppp[332];

    auto g_y_y_0_0_x_y_x_x = buffer_1100_pppp[333];

    auto g_y_y_0_0_x_y_x_y = buffer_1100_pppp[334];

    auto g_y_y_0_0_x_y_x_z = buffer_1100_pppp[335];

    auto g_y_y_0_0_x_y_y_x = buffer_1100_pppp[336];

    auto g_y_y_0_0_x_y_y_y = buffer_1100_pppp[337];

    auto g_y_y_0_0_x_y_y_z = buffer_1100_pppp[338];

    auto g_y_y_0_0_x_y_z_x = buffer_1100_pppp[339];

    auto g_y_y_0_0_x_y_z_y = buffer_1100_pppp[340];

    auto g_y_y_0_0_x_y_z_z = buffer_1100_pppp[341];

    auto g_y_y_0_0_x_z_x_x = buffer_1100_pppp[342];

    auto g_y_y_0_0_x_z_x_y = buffer_1100_pppp[343];

    auto g_y_y_0_0_x_z_x_z = buffer_1100_pppp[344];

    auto g_y_y_0_0_x_z_y_x = buffer_1100_pppp[345];

    auto g_y_y_0_0_x_z_y_y = buffer_1100_pppp[346];

    auto g_y_y_0_0_x_z_y_z = buffer_1100_pppp[347];

    auto g_y_y_0_0_x_z_z_x = buffer_1100_pppp[348];

    auto g_y_y_0_0_x_z_z_y = buffer_1100_pppp[349];

    auto g_y_y_0_0_x_z_z_z = buffer_1100_pppp[350];

    auto g_y_y_0_0_y_x_x_x = buffer_1100_pppp[351];

    auto g_y_y_0_0_y_x_x_y = buffer_1100_pppp[352];

    auto g_y_y_0_0_y_x_x_z = buffer_1100_pppp[353];

    auto g_y_y_0_0_y_x_y_x = buffer_1100_pppp[354];

    auto g_y_y_0_0_y_x_y_y = buffer_1100_pppp[355];

    auto g_y_y_0_0_y_x_y_z = buffer_1100_pppp[356];

    auto g_y_y_0_0_y_x_z_x = buffer_1100_pppp[357];

    auto g_y_y_0_0_y_x_z_y = buffer_1100_pppp[358];

    auto g_y_y_0_0_y_x_z_z = buffer_1100_pppp[359];

    auto g_y_y_0_0_y_y_x_x = buffer_1100_pppp[360];

    auto g_y_y_0_0_y_y_x_y = buffer_1100_pppp[361];

    auto g_y_y_0_0_y_y_x_z = buffer_1100_pppp[362];

    auto g_y_y_0_0_y_y_y_x = buffer_1100_pppp[363];

    auto g_y_y_0_0_y_y_y_y = buffer_1100_pppp[364];

    auto g_y_y_0_0_y_y_y_z = buffer_1100_pppp[365];

    auto g_y_y_0_0_y_y_z_x = buffer_1100_pppp[366];

    auto g_y_y_0_0_y_y_z_y = buffer_1100_pppp[367];

    auto g_y_y_0_0_y_y_z_z = buffer_1100_pppp[368];

    auto g_y_y_0_0_y_z_x_x = buffer_1100_pppp[369];

    auto g_y_y_0_0_y_z_x_y = buffer_1100_pppp[370];

    auto g_y_y_0_0_y_z_x_z = buffer_1100_pppp[371];

    auto g_y_y_0_0_y_z_y_x = buffer_1100_pppp[372];

    auto g_y_y_0_0_y_z_y_y = buffer_1100_pppp[373];

    auto g_y_y_0_0_y_z_y_z = buffer_1100_pppp[374];

    auto g_y_y_0_0_y_z_z_x = buffer_1100_pppp[375];

    auto g_y_y_0_0_y_z_z_y = buffer_1100_pppp[376];

    auto g_y_y_0_0_y_z_z_z = buffer_1100_pppp[377];

    auto g_y_y_0_0_z_x_x_x = buffer_1100_pppp[378];

    auto g_y_y_0_0_z_x_x_y = buffer_1100_pppp[379];

    auto g_y_y_0_0_z_x_x_z = buffer_1100_pppp[380];

    auto g_y_y_0_0_z_x_y_x = buffer_1100_pppp[381];

    auto g_y_y_0_0_z_x_y_y = buffer_1100_pppp[382];

    auto g_y_y_0_0_z_x_y_z = buffer_1100_pppp[383];

    auto g_y_y_0_0_z_x_z_x = buffer_1100_pppp[384];

    auto g_y_y_0_0_z_x_z_y = buffer_1100_pppp[385];

    auto g_y_y_0_0_z_x_z_z = buffer_1100_pppp[386];

    auto g_y_y_0_0_z_y_x_x = buffer_1100_pppp[387];

    auto g_y_y_0_0_z_y_x_y = buffer_1100_pppp[388];

    auto g_y_y_0_0_z_y_x_z = buffer_1100_pppp[389];

    auto g_y_y_0_0_z_y_y_x = buffer_1100_pppp[390];

    auto g_y_y_0_0_z_y_y_y = buffer_1100_pppp[391];

    auto g_y_y_0_0_z_y_y_z = buffer_1100_pppp[392];

    auto g_y_y_0_0_z_y_z_x = buffer_1100_pppp[393];

    auto g_y_y_0_0_z_y_z_y = buffer_1100_pppp[394];

    auto g_y_y_0_0_z_y_z_z = buffer_1100_pppp[395];

    auto g_y_y_0_0_z_z_x_x = buffer_1100_pppp[396];

    auto g_y_y_0_0_z_z_x_y = buffer_1100_pppp[397];

    auto g_y_y_0_0_z_z_x_z = buffer_1100_pppp[398];

    auto g_y_y_0_0_z_z_y_x = buffer_1100_pppp[399];

    auto g_y_y_0_0_z_z_y_y = buffer_1100_pppp[400];

    auto g_y_y_0_0_z_z_y_z = buffer_1100_pppp[401];

    auto g_y_y_0_0_z_z_z_x = buffer_1100_pppp[402];

    auto g_y_y_0_0_z_z_z_y = buffer_1100_pppp[403];

    auto g_y_y_0_0_z_z_z_z = buffer_1100_pppp[404];

    auto g_y_z_0_0_x_x_x_x = buffer_1100_pppp[405];

    auto g_y_z_0_0_x_x_x_y = buffer_1100_pppp[406];

    auto g_y_z_0_0_x_x_x_z = buffer_1100_pppp[407];

    auto g_y_z_0_0_x_x_y_x = buffer_1100_pppp[408];

    auto g_y_z_0_0_x_x_y_y = buffer_1100_pppp[409];

    auto g_y_z_0_0_x_x_y_z = buffer_1100_pppp[410];

    auto g_y_z_0_0_x_x_z_x = buffer_1100_pppp[411];

    auto g_y_z_0_0_x_x_z_y = buffer_1100_pppp[412];

    auto g_y_z_0_0_x_x_z_z = buffer_1100_pppp[413];

    auto g_y_z_0_0_x_y_x_x = buffer_1100_pppp[414];

    auto g_y_z_0_0_x_y_x_y = buffer_1100_pppp[415];

    auto g_y_z_0_0_x_y_x_z = buffer_1100_pppp[416];

    auto g_y_z_0_0_x_y_y_x = buffer_1100_pppp[417];

    auto g_y_z_0_0_x_y_y_y = buffer_1100_pppp[418];

    auto g_y_z_0_0_x_y_y_z = buffer_1100_pppp[419];

    auto g_y_z_0_0_x_y_z_x = buffer_1100_pppp[420];

    auto g_y_z_0_0_x_y_z_y = buffer_1100_pppp[421];

    auto g_y_z_0_0_x_y_z_z = buffer_1100_pppp[422];

    auto g_y_z_0_0_x_z_x_x = buffer_1100_pppp[423];

    auto g_y_z_0_0_x_z_x_y = buffer_1100_pppp[424];

    auto g_y_z_0_0_x_z_x_z = buffer_1100_pppp[425];

    auto g_y_z_0_0_x_z_y_x = buffer_1100_pppp[426];

    auto g_y_z_0_0_x_z_y_y = buffer_1100_pppp[427];

    auto g_y_z_0_0_x_z_y_z = buffer_1100_pppp[428];

    auto g_y_z_0_0_x_z_z_x = buffer_1100_pppp[429];

    auto g_y_z_0_0_x_z_z_y = buffer_1100_pppp[430];

    auto g_y_z_0_0_x_z_z_z = buffer_1100_pppp[431];

    auto g_y_z_0_0_y_x_x_x = buffer_1100_pppp[432];

    auto g_y_z_0_0_y_x_x_y = buffer_1100_pppp[433];

    auto g_y_z_0_0_y_x_x_z = buffer_1100_pppp[434];

    auto g_y_z_0_0_y_x_y_x = buffer_1100_pppp[435];

    auto g_y_z_0_0_y_x_y_y = buffer_1100_pppp[436];

    auto g_y_z_0_0_y_x_y_z = buffer_1100_pppp[437];

    auto g_y_z_0_0_y_x_z_x = buffer_1100_pppp[438];

    auto g_y_z_0_0_y_x_z_y = buffer_1100_pppp[439];

    auto g_y_z_0_0_y_x_z_z = buffer_1100_pppp[440];

    auto g_y_z_0_0_y_y_x_x = buffer_1100_pppp[441];

    auto g_y_z_0_0_y_y_x_y = buffer_1100_pppp[442];

    auto g_y_z_0_0_y_y_x_z = buffer_1100_pppp[443];

    auto g_y_z_0_0_y_y_y_x = buffer_1100_pppp[444];

    auto g_y_z_0_0_y_y_y_y = buffer_1100_pppp[445];

    auto g_y_z_0_0_y_y_y_z = buffer_1100_pppp[446];

    auto g_y_z_0_0_y_y_z_x = buffer_1100_pppp[447];

    auto g_y_z_0_0_y_y_z_y = buffer_1100_pppp[448];

    auto g_y_z_0_0_y_y_z_z = buffer_1100_pppp[449];

    auto g_y_z_0_0_y_z_x_x = buffer_1100_pppp[450];

    auto g_y_z_0_0_y_z_x_y = buffer_1100_pppp[451];

    auto g_y_z_0_0_y_z_x_z = buffer_1100_pppp[452];

    auto g_y_z_0_0_y_z_y_x = buffer_1100_pppp[453];

    auto g_y_z_0_0_y_z_y_y = buffer_1100_pppp[454];

    auto g_y_z_0_0_y_z_y_z = buffer_1100_pppp[455];

    auto g_y_z_0_0_y_z_z_x = buffer_1100_pppp[456];

    auto g_y_z_0_0_y_z_z_y = buffer_1100_pppp[457];

    auto g_y_z_0_0_y_z_z_z = buffer_1100_pppp[458];

    auto g_y_z_0_0_z_x_x_x = buffer_1100_pppp[459];

    auto g_y_z_0_0_z_x_x_y = buffer_1100_pppp[460];

    auto g_y_z_0_0_z_x_x_z = buffer_1100_pppp[461];

    auto g_y_z_0_0_z_x_y_x = buffer_1100_pppp[462];

    auto g_y_z_0_0_z_x_y_y = buffer_1100_pppp[463];

    auto g_y_z_0_0_z_x_y_z = buffer_1100_pppp[464];

    auto g_y_z_0_0_z_x_z_x = buffer_1100_pppp[465];

    auto g_y_z_0_0_z_x_z_y = buffer_1100_pppp[466];

    auto g_y_z_0_0_z_x_z_z = buffer_1100_pppp[467];

    auto g_y_z_0_0_z_y_x_x = buffer_1100_pppp[468];

    auto g_y_z_0_0_z_y_x_y = buffer_1100_pppp[469];

    auto g_y_z_0_0_z_y_x_z = buffer_1100_pppp[470];

    auto g_y_z_0_0_z_y_y_x = buffer_1100_pppp[471];

    auto g_y_z_0_0_z_y_y_y = buffer_1100_pppp[472];

    auto g_y_z_0_0_z_y_y_z = buffer_1100_pppp[473];

    auto g_y_z_0_0_z_y_z_x = buffer_1100_pppp[474];

    auto g_y_z_0_0_z_y_z_y = buffer_1100_pppp[475];

    auto g_y_z_0_0_z_y_z_z = buffer_1100_pppp[476];

    auto g_y_z_0_0_z_z_x_x = buffer_1100_pppp[477];

    auto g_y_z_0_0_z_z_x_y = buffer_1100_pppp[478];

    auto g_y_z_0_0_z_z_x_z = buffer_1100_pppp[479];

    auto g_y_z_0_0_z_z_y_x = buffer_1100_pppp[480];

    auto g_y_z_0_0_z_z_y_y = buffer_1100_pppp[481];

    auto g_y_z_0_0_z_z_y_z = buffer_1100_pppp[482];

    auto g_y_z_0_0_z_z_z_x = buffer_1100_pppp[483];

    auto g_y_z_0_0_z_z_z_y = buffer_1100_pppp[484];

    auto g_y_z_0_0_z_z_z_z = buffer_1100_pppp[485];

    auto g_z_x_0_0_x_x_x_x = buffer_1100_pppp[486];

    auto g_z_x_0_0_x_x_x_y = buffer_1100_pppp[487];

    auto g_z_x_0_0_x_x_x_z = buffer_1100_pppp[488];

    auto g_z_x_0_0_x_x_y_x = buffer_1100_pppp[489];

    auto g_z_x_0_0_x_x_y_y = buffer_1100_pppp[490];

    auto g_z_x_0_0_x_x_y_z = buffer_1100_pppp[491];

    auto g_z_x_0_0_x_x_z_x = buffer_1100_pppp[492];

    auto g_z_x_0_0_x_x_z_y = buffer_1100_pppp[493];

    auto g_z_x_0_0_x_x_z_z = buffer_1100_pppp[494];

    auto g_z_x_0_0_x_y_x_x = buffer_1100_pppp[495];

    auto g_z_x_0_0_x_y_x_y = buffer_1100_pppp[496];

    auto g_z_x_0_0_x_y_x_z = buffer_1100_pppp[497];

    auto g_z_x_0_0_x_y_y_x = buffer_1100_pppp[498];

    auto g_z_x_0_0_x_y_y_y = buffer_1100_pppp[499];

    auto g_z_x_0_0_x_y_y_z = buffer_1100_pppp[500];

    auto g_z_x_0_0_x_y_z_x = buffer_1100_pppp[501];

    auto g_z_x_0_0_x_y_z_y = buffer_1100_pppp[502];

    auto g_z_x_0_0_x_y_z_z = buffer_1100_pppp[503];

    auto g_z_x_0_0_x_z_x_x = buffer_1100_pppp[504];

    auto g_z_x_0_0_x_z_x_y = buffer_1100_pppp[505];

    auto g_z_x_0_0_x_z_x_z = buffer_1100_pppp[506];

    auto g_z_x_0_0_x_z_y_x = buffer_1100_pppp[507];

    auto g_z_x_0_0_x_z_y_y = buffer_1100_pppp[508];

    auto g_z_x_0_0_x_z_y_z = buffer_1100_pppp[509];

    auto g_z_x_0_0_x_z_z_x = buffer_1100_pppp[510];

    auto g_z_x_0_0_x_z_z_y = buffer_1100_pppp[511];

    auto g_z_x_0_0_x_z_z_z = buffer_1100_pppp[512];

    auto g_z_x_0_0_y_x_x_x = buffer_1100_pppp[513];

    auto g_z_x_0_0_y_x_x_y = buffer_1100_pppp[514];

    auto g_z_x_0_0_y_x_x_z = buffer_1100_pppp[515];

    auto g_z_x_0_0_y_x_y_x = buffer_1100_pppp[516];

    auto g_z_x_0_0_y_x_y_y = buffer_1100_pppp[517];

    auto g_z_x_0_0_y_x_y_z = buffer_1100_pppp[518];

    auto g_z_x_0_0_y_x_z_x = buffer_1100_pppp[519];

    auto g_z_x_0_0_y_x_z_y = buffer_1100_pppp[520];

    auto g_z_x_0_0_y_x_z_z = buffer_1100_pppp[521];

    auto g_z_x_0_0_y_y_x_x = buffer_1100_pppp[522];

    auto g_z_x_0_0_y_y_x_y = buffer_1100_pppp[523];

    auto g_z_x_0_0_y_y_x_z = buffer_1100_pppp[524];

    auto g_z_x_0_0_y_y_y_x = buffer_1100_pppp[525];

    auto g_z_x_0_0_y_y_y_y = buffer_1100_pppp[526];

    auto g_z_x_0_0_y_y_y_z = buffer_1100_pppp[527];

    auto g_z_x_0_0_y_y_z_x = buffer_1100_pppp[528];

    auto g_z_x_0_0_y_y_z_y = buffer_1100_pppp[529];

    auto g_z_x_0_0_y_y_z_z = buffer_1100_pppp[530];

    auto g_z_x_0_0_y_z_x_x = buffer_1100_pppp[531];

    auto g_z_x_0_0_y_z_x_y = buffer_1100_pppp[532];

    auto g_z_x_0_0_y_z_x_z = buffer_1100_pppp[533];

    auto g_z_x_0_0_y_z_y_x = buffer_1100_pppp[534];

    auto g_z_x_0_0_y_z_y_y = buffer_1100_pppp[535];

    auto g_z_x_0_0_y_z_y_z = buffer_1100_pppp[536];

    auto g_z_x_0_0_y_z_z_x = buffer_1100_pppp[537];

    auto g_z_x_0_0_y_z_z_y = buffer_1100_pppp[538];

    auto g_z_x_0_0_y_z_z_z = buffer_1100_pppp[539];

    auto g_z_x_0_0_z_x_x_x = buffer_1100_pppp[540];

    auto g_z_x_0_0_z_x_x_y = buffer_1100_pppp[541];

    auto g_z_x_0_0_z_x_x_z = buffer_1100_pppp[542];

    auto g_z_x_0_0_z_x_y_x = buffer_1100_pppp[543];

    auto g_z_x_0_0_z_x_y_y = buffer_1100_pppp[544];

    auto g_z_x_0_0_z_x_y_z = buffer_1100_pppp[545];

    auto g_z_x_0_0_z_x_z_x = buffer_1100_pppp[546];

    auto g_z_x_0_0_z_x_z_y = buffer_1100_pppp[547];

    auto g_z_x_0_0_z_x_z_z = buffer_1100_pppp[548];

    auto g_z_x_0_0_z_y_x_x = buffer_1100_pppp[549];

    auto g_z_x_0_0_z_y_x_y = buffer_1100_pppp[550];

    auto g_z_x_0_0_z_y_x_z = buffer_1100_pppp[551];

    auto g_z_x_0_0_z_y_y_x = buffer_1100_pppp[552];

    auto g_z_x_0_0_z_y_y_y = buffer_1100_pppp[553];

    auto g_z_x_0_0_z_y_y_z = buffer_1100_pppp[554];

    auto g_z_x_0_0_z_y_z_x = buffer_1100_pppp[555];

    auto g_z_x_0_0_z_y_z_y = buffer_1100_pppp[556];

    auto g_z_x_0_0_z_y_z_z = buffer_1100_pppp[557];

    auto g_z_x_0_0_z_z_x_x = buffer_1100_pppp[558];

    auto g_z_x_0_0_z_z_x_y = buffer_1100_pppp[559];

    auto g_z_x_0_0_z_z_x_z = buffer_1100_pppp[560];

    auto g_z_x_0_0_z_z_y_x = buffer_1100_pppp[561];

    auto g_z_x_0_0_z_z_y_y = buffer_1100_pppp[562];

    auto g_z_x_0_0_z_z_y_z = buffer_1100_pppp[563];

    auto g_z_x_0_0_z_z_z_x = buffer_1100_pppp[564];

    auto g_z_x_0_0_z_z_z_y = buffer_1100_pppp[565];

    auto g_z_x_0_0_z_z_z_z = buffer_1100_pppp[566];

    auto g_z_y_0_0_x_x_x_x = buffer_1100_pppp[567];

    auto g_z_y_0_0_x_x_x_y = buffer_1100_pppp[568];

    auto g_z_y_0_0_x_x_x_z = buffer_1100_pppp[569];

    auto g_z_y_0_0_x_x_y_x = buffer_1100_pppp[570];

    auto g_z_y_0_0_x_x_y_y = buffer_1100_pppp[571];

    auto g_z_y_0_0_x_x_y_z = buffer_1100_pppp[572];

    auto g_z_y_0_0_x_x_z_x = buffer_1100_pppp[573];

    auto g_z_y_0_0_x_x_z_y = buffer_1100_pppp[574];

    auto g_z_y_0_0_x_x_z_z = buffer_1100_pppp[575];

    auto g_z_y_0_0_x_y_x_x = buffer_1100_pppp[576];

    auto g_z_y_0_0_x_y_x_y = buffer_1100_pppp[577];

    auto g_z_y_0_0_x_y_x_z = buffer_1100_pppp[578];

    auto g_z_y_0_0_x_y_y_x = buffer_1100_pppp[579];

    auto g_z_y_0_0_x_y_y_y = buffer_1100_pppp[580];

    auto g_z_y_0_0_x_y_y_z = buffer_1100_pppp[581];

    auto g_z_y_0_0_x_y_z_x = buffer_1100_pppp[582];

    auto g_z_y_0_0_x_y_z_y = buffer_1100_pppp[583];

    auto g_z_y_0_0_x_y_z_z = buffer_1100_pppp[584];

    auto g_z_y_0_0_x_z_x_x = buffer_1100_pppp[585];

    auto g_z_y_0_0_x_z_x_y = buffer_1100_pppp[586];

    auto g_z_y_0_0_x_z_x_z = buffer_1100_pppp[587];

    auto g_z_y_0_0_x_z_y_x = buffer_1100_pppp[588];

    auto g_z_y_0_0_x_z_y_y = buffer_1100_pppp[589];

    auto g_z_y_0_0_x_z_y_z = buffer_1100_pppp[590];

    auto g_z_y_0_0_x_z_z_x = buffer_1100_pppp[591];

    auto g_z_y_0_0_x_z_z_y = buffer_1100_pppp[592];

    auto g_z_y_0_0_x_z_z_z = buffer_1100_pppp[593];

    auto g_z_y_0_0_y_x_x_x = buffer_1100_pppp[594];

    auto g_z_y_0_0_y_x_x_y = buffer_1100_pppp[595];

    auto g_z_y_0_0_y_x_x_z = buffer_1100_pppp[596];

    auto g_z_y_0_0_y_x_y_x = buffer_1100_pppp[597];

    auto g_z_y_0_0_y_x_y_y = buffer_1100_pppp[598];

    auto g_z_y_0_0_y_x_y_z = buffer_1100_pppp[599];

    auto g_z_y_0_0_y_x_z_x = buffer_1100_pppp[600];

    auto g_z_y_0_0_y_x_z_y = buffer_1100_pppp[601];

    auto g_z_y_0_0_y_x_z_z = buffer_1100_pppp[602];

    auto g_z_y_0_0_y_y_x_x = buffer_1100_pppp[603];

    auto g_z_y_0_0_y_y_x_y = buffer_1100_pppp[604];

    auto g_z_y_0_0_y_y_x_z = buffer_1100_pppp[605];

    auto g_z_y_0_0_y_y_y_x = buffer_1100_pppp[606];

    auto g_z_y_0_0_y_y_y_y = buffer_1100_pppp[607];

    auto g_z_y_0_0_y_y_y_z = buffer_1100_pppp[608];

    auto g_z_y_0_0_y_y_z_x = buffer_1100_pppp[609];

    auto g_z_y_0_0_y_y_z_y = buffer_1100_pppp[610];

    auto g_z_y_0_0_y_y_z_z = buffer_1100_pppp[611];

    auto g_z_y_0_0_y_z_x_x = buffer_1100_pppp[612];

    auto g_z_y_0_0_y_z_x_y = buffer_1100_pppp[613];

    auto g_z_y_0_0_y_z_x_z = buffer_1100_pppp[614];

    auto g_z_y_0_0_y_z_y_x = buffer_1100_pppp[615];

    auto g_z_y_0_0_y_z_y_y = buffer_1100_pppp[616];

    auto g_z_y_0_0_y_z_y_z = buffer_1100_pppp[617];

    auto g_z_y_0_0_y_z_z_x = buffer_1100_pppp[618];

    auto g_z_y_0_0_y_z_z_y = buffer_1100_pppp[619];

    auto g_z_y_0_0_y_z_z_z = buffer_1100_pppp[620];

    auto g_z_y_0_0_z_x_x_x = buffer_1100_pppp[621];

    auto g_z_y_0_0_z_x_x_y = buffer_1100_pppp[622];

    auto g_z_y_0_0_z_x_x_z = buffer_1100_pppp[623];

    auto g_z_y_0_0_z_x_y_x = buffer_1100_pppp[624];

    auto g_z_y_0_0_z_x_y_y = buffer_1100_pppp[625];

    auto g_z_y_0_0_z_x_y_z = buffer_1100_pppp[626];

    auto g_z_y_0_0_z_x_z_x = buffer_1100_pppp[627];

    auto g_z_y_0_0_z_x_z_y = buffer_1100_pppp[628];

    auto g_z_y_0_0_z_x_z_z = buffer_1100_pppp[629];

    auto g_z_y_0_0_z_y_x_x = buffer_1100_pppp[630];

    auto g_z_y_0_0_z_y_x_y = buffer_1100_pppp[631];

    auto g_z_y_0_0_z_y_x_z = buffer_1100_pppp[632];

    auto g_z_y_0_0_z_y_y_x = buffer_1100_pppp[633];

    auto g_z_y_0_0_z_y_y_y = buffer_1100_pppp[634];

    auto g_z_y_0_0_z_y_y_z = buffer_1100_pppp[635];

    auto g_z_y_0_0_z_y_z_x = buffer_1100_pppp[636];

    auto g_z_y_0_0_z_y_z_y = buffer_1100_pppp[637];

    auto g_z_y_0_0_z_y_z_z = buffer_1100_pppp[638];

    auto g_z_y_0_0_z_z_x_x = buffer_1100_pppp[639];

    auto g_z_y_0_0_z_z_x_y = buffer_1100_pppp[640];

    auto g_z_y_0_0_z_z_x_z = buffer_1100_pppp[641];

    auto g_z_y_0_0_z_z_y_x = buffer_1100_pppp[642];

    auto g_z_y_0_0_z_z_y_y = buffer_1100_pppp[643];

    auto g_z_y_0_0_z_z_y_z = buffer_1100_pppp[644];

    auto g_z_y_0_0_z_z_z_x = buffer_1100_pppp[645];

    auto g_z_y_0_0_z_z_z_y = buffer_1100_pppp[646];

    auto g_z_y_0_0_z_z_z_z = buffer_1100_pppp[647];

    auto g_z_z_0_0_x_x_x_x = buffer_1100_pppp[648];

    auto g_z_z_0_0_x_x_x_y = buffer_1100_pppp[649];

    auto g_z_z_0_0_x_x_x_z = buffer_1100_pppp[650];

    auto g_z_z_0_0_x_x_y_x = buffer_1100_pppp[651];

    auto g_z_z_0_0_x_x_y_y = buffer_1100_pppp[652];

    auto g_z_z_0_0_x_x_y_z = buffer_1100_pppp[653];

    auto g_z_z_0_0_x_x_z_x = buffer_1100_pppp[654];

    auto g_z_z_0_0_x_x_z_y = buffer_1100_pppp[655];

    auto g_z_z_0_0_x_x_z_z = buffer_1100_pppp[656];

    auto g_z_z_0_0_x_y_x_x = buffer_1100_pppp[657];

    auto g_z_z_0_0_x_y_x_y = buffer_1100_pppp[658];

    auto g_z_z_0_0_x_y_x_z = buffer_1100_pppp[659];

    auto g_z_z_0_0_x_y_y_x = buffer_1100_pppp[660];

    auto g_z_z_0_0_x_y_y_y = buffer_1100_pppp[661];

    auto g_z_z_0_0_x_y_y_z = buffer_1100_pppp[662];

    auto g_z_z_0_0_x_y_z_x = buffer_1100_pppp[663];

    auto g_z_z_0_0_x_y_z_y = buffer_1100_pppp[664];

    auto g_z_z_0_0_x_y_z_z = buffer_1100_pppp[665];

    auto g_z_z_0_0_x_z_x_x = buffer_1100_pppp[666];

    auto g_z_z_0_0_x_z_x_y = buffer_1100_pppp[667];

    auto g_z_z_0_0_x_z_x_z = buffer_1100_pppp[668];

    auto g_z_z_0_0_x_z_y_x = buffer_1100_pppp[669];

    auto g_z_z_0_0_x_z_y_y = buffer_1100_pppp[670];

    auto g_z_z_0_0_x_z_y_z = buffer_1100_pppp[671];

    auto g_z_z_0_0_x_z_z_x = buffer_1100_pppp[672];

    auto g_z_z_0_0_x_z_z_y = buffer_1100_pppp[673];

    auto g_z_z_0_0_x_z_z_z = buffer_1100_pppp[674];

    auto g_z_z_0_0_y_x_x_x = buffer_1100_pppp[675];

    auto g_z_z_0_0_y_x_x_y = buffer_1100_pppp[676];

    auto g_z_z_0_0_y_x_x_z = buffer_1100_pppp[677];

    auto g_z_z_0_0_y_x_y_x = buffer_1100_pppp[678];

    auto g_z_z_0_0_y_x_y_y = buffer_1100_pppp[679];

    auto g_z_z_0_0_y_x_y_z = buffer_1100_pppp[680];

    auto g_z_z_0_0_y_x_z_x = buffer_1100_pppp[681];

    auto g_z_z_0_0_y_x_z_y = buffer_1100_pppp[682];

    auto g_z_z_0_0_y_x_z_z = buffer_1100_pppp[683];

    auto g_z_z_0_0_y_y_x_x = buffer_1100_pppp[684];

    auto g_z_z_0_0_y_y_x_y = buffer_1100_pppp[685];

    auto g_z_z_0_0_y_y_x_z = buffer_1100_pppp[686];

    auto g_z_z_0_0_y_y_y_x = buffer_1100_pppp[687];

    auto g_z_z_0_0_y_y_y_y = buffer_1100_pppp[688];

    auto g_z_z_0_0_y_y_y_z = buffer_1100_pppp[689];

    auto g_z_z_0_0_y_y_z_x = buffer_1100_pppp[690];

    auto g_z_z_0_0_y_y_z_y = buffer_1100_pppp[691];

    auto g_z_z_0_0_y_y_z_z = buffer_1100_pppp[692];

    auto g_z_z_0_0_y_z_x_x = buffer_1100_pppp[693];

    auto g_z_z_0_0_y_z_x_y = buffer_1100_pppp[694];

    auto g_z_z_0_0_y_z_x_z = buffer_1100_pppp[695];

    auto g_z_z_0_0_y_z_y_x = buffer_1100_pppp[696];

    auto g_z_z_0_0_y_z_y_y = buffer_1100_pppp[697];

    auto g_z_z_0_0_y_z_y_z = buffer_1100_pppp[698];

    auto g_z_z_0_0_y_z_z_x = buffer_1100_pppp[699];

    auto g_z_z_0_0_y_z_z_y = buffer_1100_pppp[700];

    auto g_z_z_0_0_y_z_z_z = buffer_1100_pppp[701];

    auto g_z_z_0_0_z_x_x_x = buffer_1100_pppp[702];

    auto g_z_z_0_0_z_x_x_y = buffer_1100_pppp[703];

    auto g_z_z_0_0_z_x_x_z = buffer_1100_pppp[704];

    auto g_z_z_0_0_z_x_y_x = buffer_1100_pppp[705];

    auto g_z_z_0_0_z_x_y_y = buffer_1100_pppp[706];

    auto g_z_z_0_0_z_x_y_z = buffer_1100_pppp[707];

    auto g_z_z_0_0_z_x_z_x = buffer_1100_pppp[708];

    auto g_z_z_0_0_z_x_z_y = buffer_1100_pppp[709];

    auto g_z_z_0_0_z_x_z_z = buffer_1100_pppp[710];

    auto g_z_z_0_0_z_y_x_x = buffer_1100_pppp[711];

    auto g_z_z_0_0_z_y_x_y = buffer_1100_pppp[712];

    auto g_z_z_0_0_z_y_x_z = buffer_1100_pppp[713];

    auto g_z_z_0_0_z_y_y_x = buffer_1100_pppp[714];

    auto g_z_z_0_0_z_y_y_y = buffer_1100_pppp[715];

    auto g_z_z_0_0_z_y_y_z = buffer_1100_pppp[716];

    auto g_z_z_0_0_z_y_z_x = buffer_1100_pppp[717];

    auto g_z_z_0_0_z_y_z_y = buffer_1100_pppp[718];

    auto g_z_z_0_0_z_y_z_z = buffer_1100_pppp[719];

    auto g_z_z_0_0_z_z_x_x = buffer_1100_pppp[720];

    auto g_z_z_0_0_z_z_x_y = buffer_1100_pppp[721];

    auto g_z_z_0_0_z_z_x_z = buffer_1100_pppp[722];

    auto g_z_z_0_0_z_z_y_x = buffer_1100_pppp[723];

    auto g_z_z_0_0_z_z_y_y = buffer_1100_pppp[724];

    auto g_z_z_0_0_z_z_y_z = buffer_1100_pppp[725];

    auto g_z_z_0_0_z_z_z_x = buffer_1100_pppp[726];

    auto g_z_z_0_0_z_z_z_y = buffer_1100_pppp[727];

    auto g_z_z_0_0_z_z_z_z = buffer_1100_pppp[728];

    // integrals block (0-3)

    #pragma omp simd aligned(g_0_0_x_x, g_0_0_x_y, g_0_0_x_z, g_0_xx_x_x, g_0_xx_x_y, g_0_xx_x_z, g_x_x_0_0_x_x_x_x, g_x_x_0_0_x_x_x_y, g_x_x_0_0_x_x_x_z, g_xx_0_x_x, g_xx_0_x_y, g_xx_0_x_z, g_xx_xx_x_x, g_xx_xx_x_y, g_xx_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_x_x_x[i] = g_0_0_x_x[i] - 2.0 * g_0_xx_x_x[i] * b_exp - 2.0 * g_xx_0_x_x[i] * a_exp + 4.0 * g_xx_xx_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_x_y[i] = g_0_0_x_y[i] - 2.0 * g_0_xx_x_y[i] * b_exp - 2.0 * g_xx_0_x_y[i] * a_exp + 4.0 * g_xx_xx_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_x_z[i] = g_0_0_x_z[i] - 2.0 * g_0_xx_x_z[i] * b_exp - 2.0 * g_xx_0_x_z[i] * a_exp + 4.0 * g_xx_xx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_0_0_y_x, g_0_0_y_y, g_0_0_y_z, g_0_xx_y_x, g_0_xx_y_y, g_0_xx_y_z, g_x_x_0_0_x_x_y_x, g_x_x_0_0_x_x_y_y, g_x_x_0_0_x_x_y_z, g_xx_0_y_x, g_xx_0_y_y, g_xx_0_y_z, g_xx_xx_y_x, g_xx_xx_y_y, g_xx_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_x_y_x[i] = g_0_0_y_x[i] - 2.0 * g_0_xx_y_x[i] * b_exp - 2.0 * g_xx_0_y_x[i] * a_exp + 4.0 * g_xx_xx_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_y_y[i] = g_0_0_y_y[i] - 2.0 * g_0_xx_y_y[i] * b_exp - 2.0 * g_xx_0_y_y[i] * a_exp + 4.0 * g_xx_xx_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_y_z[i] = g_0_0_y_z[i] - 2.0 * g_0_xx_y_z[i] * b_exp - 2.0 * g_xx_0_y_z[i] * a_exp + 4.0 * g_xx_xx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_0_0_z_x, g_0_0_z_y, g_0_0_z_z, g_0_xx_z_x, g_0_xx_z_y, g_0_xx_z_z, g_x_x_0_0_x_x_z_x, g_x_x_0_0_x_x_z_y, g_x_x_0_0_x_x_z_z, g_xx_0_z_x, g_xx_0_z_y, g_xx_0_z_z, g_xx_xx_z_x, g_xx_xx_z_y, g_xx_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_x_z_x[i] = g_0_0_z_x[i] - 2.0 * g_0_xx_z_x[i] * b_exp - 2.0 * g_xx_0_z_x[i] * a_exp + 4.0 * g_xx_xx_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_z_y[i] = g_0_0_z_y[i] - 2.0 * g_0_xx_z_y[i] * b_exp - 2.0 * g_xx_0_z_y[i] * a_exp + 4.0 * g_xx_xx_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_z_z[i] = g_0_0_z_z[i] - 2.0 * g_0_xx_z_z[i] * b_exp - 2.0 * g_xx_0_z_z[i] * a_exp + 4.0 * g_xx_xx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_0_xy_x_x, g_0_xy_x_y, g_0_xy_x_z, g_x_x_0_0_x_y_x_x, g_x_x_0_0_x_y_x_y, g_x_x_0_0_x_y_x_z, g_xx_xy_x_x, g_xx_xy_x_y, g_xx_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_y_x_x[i] = -2.0 * g_0_xy_x_x[i] * b_exp + 4.0 * g_xx_xy_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_x_y[i] = -2.0 * g_0_xy_x_y[i] * b_exp + 4.0 * g_xx_xy_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_x_z[i] = -2.0 * g_0_xy_x_z[i] * b_exp + 4.0 * g_xx_xy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_0_xy_y_x, g_0_xy_y_y, g_0_xy_y_z, g_x_x_0_0_x_y_y_x, g_x_x_0_0_x_y_y_y, g_x_x_0_0_x_y_y_z, g_xx_xy_y_x, g_xx_xy_y_y, g_xx_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_y_y_x[i] = -2.0 * g_0_xy_y_x[i] * b_exp + 4.0 * g_xx_xy_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_y_y[i] = -2.0 * g_0_xy_y_y[i] * b_exp + 4.0 * g_xx_xy_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_y_z[i] = -2.0 * g_0_xy_y_z[i] * b_exp + 4.0 * g_xx_xy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_0_xy_z_x, g_0_xy_z_y, g_0_xy_z_z, g_x_x_0_0_x_y_z_x, g_x_x_0_0_x_y_z_y, g_x_x_0_0_x_y_z_z, g_xx_xy_z_x, g_xx_xy_z_y, g_xx_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_y_z_x[i] = -2.0 * g_0_xy_z_x[i] * b_exp + 4.0 * g_xx_xy_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_z_y[i] = -2.0 * g_0_xy_z_y[i] * b_exp + 4.0 * g_xx_xy_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_z_z[i] = -2.0 * g_0_xy_z_z[i] * b_exp + 4.0 * g_xx_xy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_0_xz_x_x, g_0_xz_x_y, g_0_xz_x_z, g_x_x_0_0_x_z_x_x, g_x_x_0_0_x_z_x_y, g_x_x_0_0_x_z_x_z, g_xx_xz_x_x, g_xx_xz_x_y, g_xx_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_z_x_x[i] = -2.0 * g_0_xz_x_x[i] * b_exp + 4.0 * g_xx_xz_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_x_y[i] = -2.0 * g_0_xz_x_y[i] * b_exp + 4.0 * g_xx_xz_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_x_z[i] = -2.0 * g_0_xz_x_z[i] * b_exp + 4.0 * g_xx_xz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_0_xz_y_x, g_0_xz_y_y, g_0_xz_y_z, g_x_x_0_0_x_z_y_x, g_x_x_0_0_x_z_y_y, g_x_x_0_0_x_z_y_z, g_xx_xz_y_x, g_xx_xz_y_y, g_xx_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_z_y_x[i] = -2.0 * g_0_xz_y_x[i] * b_exp + 4.0 * g_xx_xz_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_y_y[i] = -2.0 * g_0_xz_y_y[i] * b_exp + 4.0 * g_xx_xz_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_y_z[i] = -2.0 * g_0_xz_y_z[i] * b_exp + 4.0 * g_xx_xz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_0_xz_z_x, g_0_xz_z_y, g_0_xz_z_z, g_x_x_0_0_x_z_z_x, g_x_x_0_0_x_z_z_y, g_x_x_0_0_x_z_z_z, g_xx_xz_z_x, g_xx_xz_z_y, g_xx_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_z_z_x[i] = -2.0 * g_0_xz_z_x[i] * b_exp + 4.0 * g_xx_xz_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_z_y[i] = -2.0 * g_0_xz_z_y[i] * b_exp + 4.0 * g_xx_xz_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_z_z[i] = -2.0 * g_0_xz_z_z[i] * b_exp + 4.0 * g_xx_xz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (27-30)

    #pragma omp simd aligned(g_x_x_0_0_y_x_x_x, g_x_x_0_0_y_x_x_y, g_x_x_0_0_y_x_x_z, g_xy_0_x_x, g_xy_0_x_y, g_xy_0_x_z, g_xy_xx_x_x, g_xy_xx_x_y, g_xy_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_x_x_x[i] = -2.0 * g_xy_0_x_x[i] * a_exp + 4.0 * g_xy_xx_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_x_y[i] = -2.0 * g_xy_0_x_y[i] * a_exp + 4.0 * g_xy_xx_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_x_z[i] = -2.0 * g_xy_0_x_z[i] * a_exp + 4.0 * g_xy_xx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (30-33)

    #pragma omp simd aligned(g_x_x_0_0_y_x_y_x, g_x_x_0_0_y_x_y_y, g_x_x_0_0_y_x_y_z, g_xy_0_y_x, g_xy_0_y_y, g_xy_0_y_z, g_xy_xx_y_x, g_xy_xx_y_y, g_xy_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_x_y_x[i] = -2.0 * g_xy_0_y_x[i] * a_exp + 4.0 * g_xy_xx_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_y_y[i] = -2.0 * g_xy_0_y_y[i] * a_exp + 4.0 * g_xy_xx_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_y_z[i] = -2.0 * g_xy_0_y_z[i] * a_exp + 4.0 * g_xy_xx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (33-36)

    #pragma omp simd aligned(g_x_x_0_0_y_x_z_x, g_x_x_0_0_y_x_z_y, g_x_x_0_0_y_x_z_z, g_xy_0_z_x, g_xy_0_z_y, g_xy_0_z_z, g_xy_xx_z_x, g_xy_xx_z_y, g_xy_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_x_z_x[i] = -2.0 * g_xy_0_z_x[i] * a_exp + 4.0 * g_xy_xx_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_z_y[i] = -2.0 * g_xy_0_z_y[i] * a_exp + 4.0 * g_xy_xx_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_z_z[i] = -2.0 * g_xy_0_z_z[i] * a_exp + 4.0 * g_xy_xx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (36-39)

    #pragma omp simd aligned(g_x_x_0_0_y_y_x_x, g_x_x_0_0_y_y_x_y, g_x_x_0_0_y_y_x_z, g_xy_xy_x_x, g_xy_xy_x_y, g_xy_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_y_x_x[i] = 4.0 * g_xy_xy_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_x_y[i] = 4.0 * g_xy_xy_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_x_z[i] = 4.0 * g_xy_xy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (39-42)

    #pragma omp simd aligned(g_x_x_0_0_y_y_y_x, g_x_x_0_0_y_y_y_y, g_x_x_0_0_y_y_y_z, g_xy_xy_y_x, g_xy_xy_y_y, g_xy_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_y_y_x[i] = 4.0 * g_xy_xy_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_y_y[i] = 4.0 * g_xy_xy_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_y_z[i] = 4.0 * g_xy_xy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (42-45)

    #pragma omp simd aligned(g_x_x_0_0_y_y_z_x, g_x_x_0_0_y_y_z_y, g_x_x_0_0_y_y_z_z, g_xy_xy_z_x, g_xy_xy_z_y, g_xy_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_y_z_x[i] = 4.0 * g_xy_xy_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_z_y[i] = 4.0 * g_xy_xy_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_z_z[i] = 4.0 * g_xy_xy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (45-48)

    #pragma omp simd aligned(g_x_x_0_0_y_z_x_x, g_x_x_0_0_y_z_x_y, g_x_x_0_0_y_z_x_z, g_xy_xz_x_x, g_xy_xz_x_y, g_xy_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_z_x_x[i] = 4.0 * g_xy_xz_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_x_y[i] = 4.0 * g_xy_xz_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_x_z[i] = 4.0 * g_xy_xz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (48-51)

    #pragma omp simd aligned(g_x_x_0_0_y_z_y_x, g_x_x_0_0_y_z_y_y, g_x_x_0_0_y_z_y_z, g_xy_xz_y_x, g_xy_xz_y_y, g_xy_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_z_y_x[i] = 4.0 * g_xy_xz_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_y_y[i] = 4.0 * g_xy_xz_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_y_z[i] = 4.0 * g_xy_xz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (51-54)

    #pragma omp simd aligned(g_x_x_0_0_y_z_z_x, g_x_x_0_0_y_z_z_y, g_x_x_0_0_y_z_z_z, g_xy_xz_z_x, g_xy_xz_z_y, g_xy_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_z_z_x[i] = 4.0 * g_xy_xz_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_z_y[i] = 4.0 * g_xy_xz_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_z_z[i] = 4.0 * g_xy_xz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (54-57)

    #pragma omp simd aligned(g_x_x_0_0_z_x_x_x, g_x_x_0_0_z_x_x_y, g_x_x_0_0_z_x_x_z, g_xz_0_x_x, g_xz_0_x_y, g_xz_0_x_z, g_xz_xx_x_x, g_xz_xx_x_y, g_xz_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_x_x_x[i] = -2.0 * g_xz_0_x_x[i] * a_exp + 4.0 * g_xz_xx_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_x_y[i] = -2.0 * g_xz_0_x_y[i] * a_exp + 4.0 * g_xz_xx_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_x_z[i] = -2.0 * g_xz_0_x_z[i] * a_exp + 4.0 * g_xz_xx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (57-60)

    #pragma omp simd aligned(g_x_x_0_0_z_x_y_x, g_x_x_0_0_z_x_y_y, g_x_x_0_0_z_x_y_z, g_xz_0_y_x, g_xz_0_y_y, g_xz_0_y_z, g_xz_xx_y_x, g_xz_xx_y_y, g_xz_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_x_y_x[i] = -2.0 * g_xz_0_y_x[i] * a_exp + 4.0 * g_xz_xx_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_y_y[i] = -2.0 * g_xz_0_y_y[i] * a_exp + 4.0 * g_xz_xx_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_y_z[i] = -2.0 * g_xz_0_y_z[i] * a_exp + 4.0 * g_xz_xx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (60-63)

    #pragma omp simd aligned(g_x_x_0_0_z_x_z_x, g_x_x_0_0_z_x_z_y, g_x_x_0_0_z_x_z_z, g_xz_0_z_x, g_xz_0_z_y, g_xz_0_z_z, g_xz_xx_z_x, g_xz_xx_z_y, g_xz_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_x_z_x[i] = -2.0 * g_xz_0_z_x[i] * a_exp + 4.0 * g_xz_xx_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_z_y[i] = -2.0 * g_xz_0_z_y[i] * a_exp + 4.0 * g_xz_xx_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_z_z[i] = -2.0 * g_xz_0_z_z[i] * a_exp + 4.0 * g_xz_xx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (63-66)

    #pragma omp simd aligned(g_x_x_0_0_z_y_x_x, g_x_x_0_0_z_y_x_y, g_x_x_0_0_z_y_x_z, g_xz_xy_x_x, g_xz_xy_x_y, g_xz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_y_x_x[i] = 4.0 * g_xz_xy_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_x_y[i] = 4.0 * g_xz_xy_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_x_z[i] = 4.0 * g_xz_xy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (66-69)

    #pragma omp simd aligned(g_x_x_0_0_z_y_y_x, g_x_x_0_0_z_y_y_y, g_x_x_0_0_z_y_y_z, g_xz_xy_y_x, g_xz_xy_y_y, g_xz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_y_y_x[i] = 4.0 * g_xz_xy_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_y_y[i] = 4.0 * g_xz_xy_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_y_z[i] = 4.0 * g_xz_xy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (69-72)

    #pragma omp simd aligned(g_x_x_0_0_z_y_z_x, g_x_x_0_0_z_y_z_y, g_x_x_0_0_z_y_z_z, g_xz_xy_z_x, g_xz_xy_z_y, g_xz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_y_z_x[i] = 4.0 * g_xz_xy_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_z_y[i] = 4.0 * g_xz_xy_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_z_z[i] = 4.0 * g_xz_xy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (72-75)

    #pragma omp simd aligned(g_x_x_0_0_z_z_x_x, g_x_x_0_0_z_z_x_y, g_x_x_0_0_z_z_x_z, g_xz_xz_x_x, g_xz_xz_x_y, g_xz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_z_x_x[i] = 4.0 * g_xz_xz_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_x_y[i] = 4.0 * g_xz_xz_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_x_z[i] = 4.0 * g_xz_xz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (75-78)

    #pragma omp simd aligned(g_x_x_0_0_z_z_y_x, g_x_x_0_0_z_z_y_y, g_x_x_0_0_z_z_y_z, g_xz_xz_y_x, g_xz_xz_y_y, g_xz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_z_y_x[i] = 4.0 * g_xz_xz_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_y_y[i] = 4.0 * g_xz_xz_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_y_z[i] = 4.0 * g_xz_xz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (78-81)

    #pragma omp simd aligned(g_x_x_0_0_z_z_z_x, g_x_x_0_0_z_z_z_y, g_x_x_0_0_z_z_z_z, g_xz_xz_z_x, g_xz_xz_z_y, g_xz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_z_z_x[i] = 4.0 * g_xz_xz_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_z_y[i] = 4.0 * g_xz_xz_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_z_z[i] = 4.0 * g_xz_xz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (81-84)

    #pragma omp simd aligned(g_0_xy_x_x, g_0_xy_x_y, g_0_xy_x_z, g_x_y_0_0_x_x_x_x, g_x_y_0_0_x_x_x_y, g_x_y_0_0_x_x_x_z, g_xx_xy_x_x, g_xx_xy_x_y, g_xx_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_x_x_x[i] = -2.0 * g_0_xy_x_x[i] * b_exp + 4.0 * g_xx_xy_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_x_y[i] = -2.0 * g_0_xy_x_y[i] * b_exp + 4.0 * g_xx_xy_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_x_z[i] = -2.0 * g_0_xy_x_z[i] * b_exp + 4.0 * g_xx_xy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (84-87)

    #pragma omp simd aligned(g_0_xy_y_x, g_0_xy_y_y, g_0_xy_y_z, g_x_y_0_0_x_x_y_x, g_x_y_0_0_x_x_y_y, g_x_y_0_0_x_x_y_z, g_xx_xy_y_x, g_xx_xy_y_y, g_xx_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_x_y_x[i] = -2.0 * g_0_xy_y_x[i] * b_exp + 4.0 * g_xx_xy_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_y_y[i] = -2.0 * g_0_xy_y_y[i] * b_exp + 4.0 * g_xx_xy_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_y_z[i] = -2.0 * g_0_xy_y_z[i] * b_exp + 4.0 * g_xx_xy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (87-90)

    #pragma omp simd aligned(g_0_xy_z_x, g_0_xy_z_y, g_0_xy_z_z, g_x_y_0_0_x_x_z_x, g_x_y_0_0_x_x_z_y, g_x_y_0_0_x_x_z_z, g_xx_xy_z_x, g_xx_xy_z_y, g_xx_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_x_z_x[i] = -2.0 * g_0_xy_z_x[i] * b_exp + 4.0 * g_xx_xy_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_z_y[i] = -2.0 * g_0_xy_z_y[i] * b_exp + 4.0 * g_xx_xy_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_z_z[i] = -2.0 * g_0_xy_z_z[i] * b_exp + 4.0 * g_xx_xy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (90-93)

    #pragma omp simd aligned(g_0_0_x_x, g_0_0_x_y, g_0_0_x_z, g_0_yy_x_x, g_0_yy_x_y, g_0_yy_x_z, g_x_y_0_0_x_y_x_x, g_x_y_0_0_x_y_x_y, g_x_y_0_0_x_y_x_z, g_xx_0_x_x, g_xx_0_x_y, g_xx_0_x_z, g_xx_yy_x_x, g_xx_yy_x_y, g_xx_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_y_x_x[i] = g_0_0_x_x[i] - 2.0 * g_0_yy_x_x[i] * b_exp - 2.0 * g_xx_0_x_x[i] * a_exp + 4.0 * g_xx_yy_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_x_y[i] = g_0_0_x_y[i] - 2.0 * g_0_yy_x_y[i] * b_exp - 2.0 * g_xx_0_x_y[i] * a_exp + 4.0 * g_xx_yy_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_x_z[i] = g_0_0_x_z[i] - 2.0 * g_0_yy_x_z[i] * b_exp - 2.0 * g_xx_0_x_z[i] * a_exp + 4.0 * g_xx_yy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (93-96)

    #pragma omp simd aligned(g_0_0_y_x, g_0_0_y_y, g_0_0_y_z, g_0_yy_y_x, g_0_yy_y_y, g_0_yy_y_z, g_x_y_0_0_x_y_y_x, g_x_y_0_0_x_y_y_y, g_x_y_0_0_x_y_y_z, g_xx_0_y_x, g_xx_0_y_y, g_xx_0_y_z, g_xx_yy_y_x, g_xx_yy_y_y, g_xx_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_y_y_x[i] = g_0_0_y_x[i] - 2.0 * g_0_yy_y_x[i] * b_exp - 2.0 * g_xx_0_y_x[i] * a_exp + 4.0 * g_xx_yy_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_y_y[i] = g_0_0_y_y[i] - 2.0 * g_0_yy_y_y[i] * b_exp - 2.0 * g_xx_0_y_y[i] * a_exp + 4.0 * g_xx_yy_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_y_z[i] = g_0_0_y_z[i] - 2.0 * g_0_yy_y_z[i] * b_exp - 2.0 * g_xx_0_y_z[i] * a_exp + 4.0 * g_xx_yy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (96-99)

    #pragma omp simd aligned(g_0_0_z_x, g_0_0_z_y, g_0_0_z_z, g_0_yy_z_x, g_0_yy_z_y, g_0_yy_z_z, g_x_y_0_0_x_y_z_x, g_x_y_0_0_x_y_z_y, g_x_y_0_0_x_y_z_z, g_xx_0_z_x, g_xx_0_z_y, g_xx_0_z_z, g_xx_yy_z_x, g_xx_yy_z_y, g_xx_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_y_z_x[i] = g_0_0_z_x[i] - 2.0 * g_0_yy_z_x[i] * b_exp - 2.0 * g_xx_0_z_x[i] * a_exp + 4.0 * g_xx_yy_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_z_y[i] = g_0_0_z_y[i] - 2.0 * g_0_yy_z_y[i] * b_exp - 2.0 * g_xx_0_z_y[i] * a_exp + 4.0 * g_xx_yy_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_z_z[i] = g_0_0_z_z[i] - 2.0 * g_0_yy_z_z[i] * b_exp - 2.0 * g_xx_0_z_z[i] * a_exp + 4.0 * g_xx_yy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (99-102)

    #pragma omp simd aligned(g_0_yz_x_x, g_0_yz_x_y, g_0_yz_x_z, g_x_y_0_0_x_z_x_x, g_x_y_0_0_x_z_x_y, g_x_y_0_0_x_z_x_z, g_xx_yz_x_x, g_xx_yz_x_y, g_xx_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_z_x_x[i] = -2.0 * g_0_yz_x_x[i] * b_exp + 4.0 * g_xx_yz_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_x_y[i] = -2.0 * g_0_yz_x_y[i] * b_exp + 4.0 * g_xx_yz_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_x_z[i] = -2.0 * g_0_yz_x_z[i] * b_exp + 4.0 * g_xx_yz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (102-105)

    #pragma omp simd aligned(g_0_yz_y_x, g_0_yz_y_y, g_0_yz_y_z, g_x_y_0_0_x_z_y_x, g_x_y_0_0_x_z_y_y, g_x_y_0_0_x_z_y_z, g_xx_yz_y_x, g_xx_yz_y_y, g_xx_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_z_y_x[i] = -2.0 * g_0_yz_y_x[i] * b_exp + 4.0 * g_xx_yz_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_y_y[i] = -2.0 * g_0_yz_y_y[i] * b_exp + 4.0 * g_xx_yz_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_y_z[i] = -2.0 * g_0_yz_y_z[i] * b_exp + 4.0 * g_xx_yz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (105-108)

    #pragma omp simd aligned(g_0_yz_z_x, g_0_yz_z_y, g_0_yz_z_z, g_x_y_0_0_x_z_z_x, g_x_y_0_0_x_z_z_y, g_x_y_0_0_x_z_z_z, g_xx_yz_z_x, g_xx_yz_z_y, g_xx_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_z_z_x[i] = -2.0 * g_0_yz_z_x[i] * b_exp + 4.0 * g_xx_yz_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_z_y[i] = -2.0 * g_0_yz_z_y[i] * b_exp + 4.0 * g_xx_yz_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_z_z[i] = -2.0 * g_0_yz_z_z[i] * b_exp + 4.0 * g_xx_yz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (108-111)

    #pragma omp simd aligned(g_x_y_0_0_y_x_x_x, g_x_y_0_0_y_x_x_y, g_x_y_0_0_y_x_x_z, g_xy_xy_x_x, g_xy_xy_x_y, g_xy_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_x_x_x[i] = 4.0 * g_xy_xy_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_x_y[i] = 4.0 * g_xy_xy_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_x_z[i] = 4.0 * g_xy_xy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (111-114)

    #pragma omp simd aligned(g_x_y_0_0_y_x_y_x, g_x_y_0_0_y_x_y_y, g_x_y_0_0_y_x_y_z, g_xy_xy_y_x, g_xy_xy_y_y, g_xy_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_x_y_x[i] = 4.0 * g_xy_xy_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_y_y[i] = 4.0 * g_xy_xy_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_y_z[i] = 4.0 * g_xy_xy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (114-117)

    #pragma omp simd aligned(g_x_y_0_0_y_x_z_x, g_x_y_0_0_y_x_z_y, g_x_y_0_0_y_x_z_z, g_xy_xy_z_x, g_xy_xy_z_y, g_xy_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_x_z_x[i] = 4.0 * g_xy_xy_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_z_y[i] = 4.0 * g_xy_xy_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_z_z[i] = 4.0 * g_xy_xy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (117-120)

    #pragma omp simd aligned(g_x_y_0_0_y_y_x_x, g_x_y_0_0_y_y_x_y, g_x_y_0_0_y_y_x_z, g_xy_0_x_x, g_xy_0_x_y, g_xy_0_x_z, g_xy_yy_x_x, g_xy_yy_x_y, g_xy_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_y_x_x[i] = -2.0 * g_xy_0_x_x[i] * a_exp + 4.0 * g_xy_yy_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_x_y[i] = -2.0 * g_xy_0_x_y[i] * a_exp + 4.0 * g_xy_yy_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_x_z[i] = -2.0 * g_xy_0_x_z[i] * a_exp + 4.0 * g_xy_yy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (120-123)

    #pragma omp simd aligned(g_x_y_0_0_y_y_y_x, g_x_y_0_0_y_y_y_y, g_x_y_0_0_y_y_y_z, g_xy_0_y_x, g_xy_0_y_y, g_xy_0_y_z, g_xy_yy_y_x, g_xy_yy_y_y, g_xy_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_y_y_x[i] = -2.0 * g_xy_0_y_x[i] * a_exp + 4.0 * g_xy_yy_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_y_y[i] = -2.0 * g_xy_0_y_y[i] * a_exp + 4.0 * g_xy_yy_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_y_z[i] = -2.0 * g_xy_0_y_z[i] * a_exp + 4.0 * g_xy_yy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (123-126)

    #pragma omp simd aligned(g_x_y_0_0_y_y_z_x, g_x_y_0_0_y_y_z_y, g_x_y_0_0_y_y_z_z, g_xy_0_z_x, g_xy_0_z_y, g_xy_0_z_z, g_xy_yy_z_x, g_xy_yy_z_y, g_xy_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_y_z_x[i] = -2.0 * g_xy_0_z_x[i] * a_exp + 4.0 * g_xy_yy_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_z_y[i] = -2.0 * g_xy_0_z_y[i] * a_exp + 4.0 * g_xy_yy_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_z_z[i] = -2.0 * g_xy_0_z_z[i] * a_exp + 4.0 * g_xy_yy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (126-129)

    #pragma omp simd aligned(g_x_y_0_0_y_z_x_x, g_x_y_0_0_y_z_x_y, g_x_y_0_0_y_z_x_z, g_xy_yz_x_x, g_xy_yz_x_y, g_xy_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_z_x_x[i] = 4.0 * g_xy_yz_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_x_y[i] = 4.0 * g_xy_yz_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_x_z[i] = 4.0 * g_xy_yz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (129-132)

    #pragma omp simd aligned(g_x_y_0_0_y_z_y_x, g_x_y_0_0_y_z_y_y, g_x_y_0_0_y_z_y_z, g_xy_yz_y_x, g_xy_yz_y_y, g_xy_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_z_y_x[i] = 4.0 * g_xy_yz_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_y_y[i] = 4.0 * g_xy_yz_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_y_z[i] = 4.0 * g_xy_yz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (132-135)

    #pragma omp simd aligned(g_x_y_0_0_y_z_z_x, g_x_y_0_0_y_z_z_y, g_x_y_0_0_y_z_z_z, g_xy_yz_z_x, g_xy_yz_z_y, g_xy_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_z_z_x[i] = 4.0 * g_xy_yz_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_z_y[i] = 4.0 * g_xy_yz_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_z_z[i] = 4.0 * g_xy_yz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (135-138)

    #pragma omp simd aligned(g_x_y_0_0_z_x_x_x, g_x_y_0_0_z_x_x_y, g_x_y_0_0_z_x_x_z, g_xz_xy_x_x, g_xz_xy_x_y, g_xz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_x_x_x[i] = 4.0 * g_xz_xy_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_x_y[i] = 4.0 * g_xz_xy_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_x_z[i] = 4.0 * g_xz_xy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (138-141)

    #pragma omp simd aligned(g_x_y_0_0_z_x_y_x, g_x_y_0_0_z_x_y_y, g_x_y_0_0_z_x_y_z, g_xz_xy_y_x, g_xz_xy_y_y, g_xz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_x_y_x[i] = 4.0 * g_xz_xy_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_y_y[i] = 4.0 * g_xz_xy_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_y_z[i] = 4.0 * g_xz_xy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (141-144)

    #pragma omp simd aligned(g_x_y_0_0_z_x_z_x, g_x_y_0_0_z_x_z_y, g_x_y_0_0_z_x_z_z, g_xz_xy_z_x, g_xz_xy_z_y, g_xz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_x_z_x[i] = 4.0 * g_xz_xy_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_z_y[i] = 4.0 * g_xz_xy_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_z_z[i] = 4.0 * g_xz_xy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (144-147)

    #pragma omp simd aligned(g_x_y_0_0_z_y_x_x, g_x_y_0_0_z_y_x_y, g_x_y_0_0_z_y_x_z, g_xz_0_x_x, g_xz_0_x_y, g_xz_0_x_z, g_xz_yy_x_x, g_xz_yy_x_y, g_xz_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_y_x_x[i] = -2.0 * g_xz_0_x_x[i] * a_exp + 4.0 * g_xz_yy_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_x_y[i] = -2.0 * g_xz_0_x_y[i] * a_exp + 4.0 * g_xz_yy_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_x_z[i] = -2.0 * g_xz_0_x_z[i] * a_exp + 4.0 * g_xz_yy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (147-150)

    #pragma omp simd aligned(g_x_y_0_0_z_y_y_x, g_x_y_0_0_z_y_y_y, g_x_y_0_0_z_y_y_z, g_xz_0_y_x, g_xz_0_y_y, g_xz_0_y_z, g_xz_yy_y_x, g_xz_yy_y_y, g_xz_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_y_y_x[i] = -2.0 * g_xz_0_y_x[i] * a_exp + 4.0 * g_xz_yy_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_y_y[i] = -2.0 * g_xz_0_y_y[i] * a_exp + 4.0 * g_xz_yy_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_y_z[i] = -2.0 * g_xz_0_y_z[i] * a_exp + 4.0 * g_xz_yy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (150-153)

    #pragma omp simd aligned(g_x_y_0_0_z_y_z_x, g_x_y_0_0_z_y_z_y, g_x_y_0_0_z_y_z_z, g_xz_0_z_x, g_xz_0_z_y, g_xz_0_z_z, g_xz_yy_z_x, g_xz_yy_z_y, g_xz_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_y_z_x[i] = -2.0 * g_xz_0_z_x[i] * a_exp + 4.0 * g_xz_yy_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_z_y[i] = -2.0 * g_xz_0_z_y[i] * a_exp + 4.0 * g_xz_yy_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_z_z[i] = -2.0 * g_xz_0_z_z[i] * a_exp + 4.0 * g_xz_yy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (153-156)

    #pragma omp simd aligned(g_x_y_0_0_z_z_x_x, g_x_y_0_0_z_z_x_y, g_x_y_0_0_z_z_x_z, g_xz_yz_x_x, g_xz_yz_x_y, g_xz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_z_x_x[i] = 4.0 * g_xz_yz_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_x_y[i] = 4.0 * g_xz_yz_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_x_z[i] = 4.0 * g_xz_yz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (156-159)

    #pragma omp simd aligned(g_x_y_0_0_z_z_y_x, g_x_y_0_0_z_z_y_y, g_x_y_0_0_z_z_y_z, g_xz_yz_y_x, g_xz_yz_y_y, g_xz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_z_y_x[i] = 4.0 * g_xz_yz_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_y_y[i] = 4.0 * g_xz_yz_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_y_z[i] = 4.0 * g_xz_yz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (159-162)

    #pragma omp simd aligned(g_x_y_0_0_z_z_z_x, g_x_y_0_0_z_z_z_y, g_x_y_0_0_z_z_z_z, g_xz_yz_z_x, g_xz_yz_z_y, g_xz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_z_z_x[i] = 4.0 * g_xz_yz_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_z_y[i] = 4.0 * g_xz_yz_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_z_z[i] = 4.0 * g_xz_yz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (162-165)

    #pragma omp simd aligned(g_0_xz_x_x, g_0_xz_x_y, g_0_xz_x_z, g_x_z_0_0_x_x_x_x, g_x_z_0_0_x_x_x_y, g_x_z_0_0_x_x_x_z, g_xx_xz_x_x, g_xx_xz_x_y, g_xx_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_x_x_x[i] = -2.0 * g_0_xz_x_x[i] * b_exp + 4.0 * g_xx_xz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_x_y[i] = -2.0 * g_0_xz_x_y[i] * b_exp + 4.0 * g_xx_xz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_x_z[i] = -2.0 * g_0_xz_x_z[i] * b_exp + 4.0 * g_xx_xz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (165-168)

    #pragma omp simd aligned(g_0_xz_y_x, g_0_xz_y_y, g_0_xz_y_z, g_x_z_0_0_x_x_y_x, g_x_z_0_0_x_x_y_y, g_x_z_0_0_x_x_y_z, g_xx_xz_y_x, g_xx_xz_y_y, g_xx_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_x_y_x[i] = -2.0 * g_0_xz_y_x[i] * b_exp + 4.0 * g_xx_xz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_y_y[i] = -2.0 * g_0_xz_y_y[i] * b_exp + 4.0 * g_xx_xz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_y_z[i] = -2.0 * g_0_xz_y_z[i] * b_exp + 4.0 * g_xx_xz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (168-171)

    #pragma omp simd aligned(g_0_xz_z_x, g_0_xz_z_y, g_0_xz_z_z, g_x_z_0_0_x_x_z_x, g_x_z_0_0_x_x_z_y, g_x_z_0_0_x_x_z_z, g_xx_xz_z_x, g_xx_xz_z_y, g_xx_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_x_z_x[i] = -2.0 * g_0_xz_z_x[i] * b_exp + 4.0 * g_xx_xz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_z_y[i] = -2.0 * g_0_xz_z_y[i] * b_exp + 4.0 * g_xx_xz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_z_z[i] = -2.0 * g_0_xz_z_z[i] * b_exp + 4.0 * g_xx_xz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (171-174)

    #pragma omp simd aligned(g_0_yz_x_x, g_0_yz_x_y, g_0_yz_x_z, g_x_z_0_0_x_y_x_x, g_x_z_0_0_x_y_x_y, g_x_z_0_0_x_y_x_z, g_xx_yz_x_x, g_xx_yz_x_y, g_xx_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_y_x_x[i] = -2.0 * g_0_yz_x_x[i] * b_exp + 4.0 * g_xx_yz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_x_y[i] = -2.0 * g_0_yz_x_y[i] * b_exp + 4.0 * g_xx_yz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_x_z[i] = -2.0 * g_0_yz_x_z[i] * b_exp + 4.0 * g_xx_yz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (174-177)

    #pragma omp simd aligned(g_0_yz_y_x, g_0_yz_y_y, g_0_yz_y_z, g_x_z_0_0_x_y_y_x, g_x_z_0_0_x_y_y_y, g_x_z_0_0_x_y_y_z, g_xx_yz_y_x, g_xx_yz_y_y, g_xx_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_y_y_x[i] = -2.0 * g_0_yz_y_x[i] * b_exp + 4.0 * g_xx_yz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_y_y[i] = -2.0 * g_0_yz_y_y[i] * b_exp + 4.0 * g_xx_yz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_y_z[i] = -2.0 * g_0_yz_y_z[i] * b_exp + 4.0 * g_xx_yz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (177-180)

    #pragma omp simd aligned(g_0_yz_z_x, g_0_yz_z_y, g_0_yz_z_z, g_x_z_0_0_x_y_z_x, g_x_z_0_0_x_y_z_y, g_x_z_0_0_x_y_z_z, g_xx_yz_z_x, g_xx_yz_z_y, g_xx_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_y_z_x[i] = -2.0 * g_0_yz_z_x[i] * b_exp + 4.0 * g_xx_yz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_z_y[i] = -2.0 * g_0_yz_z_y[i] * b_exp + 4.0 * g_xx_yz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_z_z[i] = -2.0 * g_0_yz_z_z[i] * b_exp + 4.0 * g_xx_yz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (180-183)

    #pragma omp simd aligned(g_0_0_x_x, g_0_0_x_y, g_0_0_x_z, g_0_zz_x_x, g_0_zz_x_y, g_0_zz_x_z, g_x_z_0_0_x_z_x_x, g_x_z_0_0_x_z_x_y, g_x_z_0_0_x_z_x_z, g_xx_0_x_x, g_xx_0_x_y, g_xx_0_x_z, g_xx_zz_x_x, g_xx_zz_x_y, g_xx_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_z_x_x[i] = g_0_0_x_x[i] - 2.0 * g_0_zz_x_x[i] * b_exp - 2.0 * g_xx_0_x_x[i] * a_exp + 4.0 * g_xx_zz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_x_y[i] = g_0_0_x_y[i] - 2.0 * g_0_zz_x_y[i] * b_exp - 2.0 * g_xx_0_x_y[i] * a_exp + 4.0 * g_xx_zz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_x_z[i] = g_0_0_x_z[i] - 2.0 * g_0_zz_x_z[i] * b_exp - 2.0 * g_xx_0_x_z[i] * a_exp + 4.0 * g_xx_zz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (183-186)

    #pragma omp simd aligned(g_0_0_y_x, g_0_0_y_y, g_0_0_y_z, g_0_zz_y_x, g_0_zz_y_y, g_0_zz_y_z, g_x_z_0_0_x_z_y_x, g_x_z_0_0_x_z_y_y, g_x_z_0_0_x_z_y_z, g_xx_0_y_x, g_xx_0_y_y, g_xx_0_y_z, g_xx_zz_y_x, g_xx_zz_y_y, g_xx_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_z_y_x[i] = g_0_0_y_x[i] - 2.0 * g_0_zz_y_x[i] * b_exp - 2.0 * g_xx_0_y_x[i] * a_exp + 4.0 * g_xx_zz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_y_y[i] = g_0_0_y_y[i] - 2.0 * g_0_zz_y_y[i] * b_exp - 2.0 * g_xx_0_y_y[i] * a_exp + 4.0 * g_xx_zz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_y_z[i] = g_0_0_y_z[i] - 2.0 * g_0_zz_y_z[i] * b_exp - 2.0 * g_xx_0_y_z[i] * a_exp + 4.0 * g_xx_zz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (186-189)

    #pragma omp simd aligned(g_0_0_z_x, g_0_0_z_y, g_0_0_z_z, g_0_zz_z_x, g_0_zz_z_y, g_0_zz_z_z, g_x_z_0_0_x_z_z_x, g_x_z_0_0_x_z_z_y, g_x_z_0_0_x_z_z_z, g_xx_0_z_x, g_xx_0_z_y, g_xx_0_z_z, g_xx_zz_z_x, g_xx_zz_z_y, g_xx_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_z_z_x[i] = g_0_0_z_x[i] - 2.0 * g_0_zz_z_x[i] * b_exp - 2.0 * g_xx_0_z_x[i] * a_exp + 4.0 * g_xx_zz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_z_y[i] = g_0_0_z_y[i] - 2.0 * g_0_zz_z_y[i] * b_exp - 2.0 * g_xx_0_z_y[i] * a_exp + 4.0 * g_xx_zz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_z_z[i] = g_0_0_z_z[i] - 2.0 * g_0_zz_z_z[i] * b_exp - 2.0 * g_xx_0_z_z[i] * a_exp + 4.0 * g_xx_zz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (189-192)

    #pragma omp simd aligned(g_x_z_0_0_y_x_x_x, g_x_z_0_0_y_x_x_y, g_x_z_0_0_y_x_x_z, g_xy_xz_x_x, g_xy_xz_x_y, g_xy_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_x_x_x[i] = 4.0 * g_xy_xz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_x_y[i] = 4.0 * g_xy_xz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_x_z[i] = 4.0 * g_xy_xz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (192-195)

    #pragma omp simd aligned(g_x_z_0_0_y_x_y_x, g_x_z_0_0_y_x_y_y, g_x_z_0_0_y_x_y_z, g_xy_xz_y_x, g_xy_xz_y_y, g_xy_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_x_y_x[i] = 4.0 * g_xy_xz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_y_y[i] = 4.0 * g_xy_xz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_y_z[i] = 4.0 * g_xy_xz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (195-198)

    #pragma omp simd aligned(g_x_z_0_0_y_x_z_x, g_x_z_0_0_y_x_z_y, g_x_z_0_0_y_x_z_z, g_xy_xz_z_x, g_xy_xz_z_y, g_xy_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_x_z_x[i] = 4.0 * g_xy_xz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_z_y[i] = 4.0 * g_xy_xz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_z_z[i] = 4.0 * g_xy_xz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (198-201)

    #pragma omp simd aligned(g_x_z_0_0_y_y_x_x, g_x_z_0_0_y_y_x_y, g_x_z_0_0_y_y_x_z, g_xy_yz_x_x, g_xy_yz_x_y, g_xy_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_y_x_x[i] = 4.0 * g_xy_yz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_x_y[i] = 4.0 * g_xy_yz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_x_z[i] = 4.0 * g_xy_yz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (201-204)

    #pragma omp simd aligned(g_x_z_0_0_y_y_y_x, g_x_z_0_0_y_y_y_y, g_x_z_0_0_y_y_y_z, g_xy_yz_y_x, g_xy_yz_y_y, g_xy_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_y_y_x[i] = 4.0 * g_xy_yz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_y_y[i] = 4.0 * g_xy_yz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_y_z[i] = 4.0 * g_xy_yz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (204-207)

    #pragma omp simd aligned(g_x_z_0_0_y_y_z_x, g_x_z_0_0_y_y_z_y, g_x_z_0_0_y_y_z_z, g_xy_yz_z_x, g_xy_yz_z_y, g_xy_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_y_z_x[i] = 4.0 * g_xy_yz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_z_y[i] = 4.0 * g_xy_yz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_z_z[i] = 4.0 * g_xy_yz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (207-210)

    #pragma omp simd aligned(g_x_z_0_0_y_z_x_x, g_x_z_0_0_y_z_x_y, g_x_z_0_0_y_z_x_z, g_xy_0_x_x, g_xy_0_x_y, g_xy_0_x_z, g_xy_zz_x_x, g_xy_zz_x_y, g_xy_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_z_x_x[i] = -2.0 * g_xy_0_x_x[i] * a_exp + 4.0 * g_xy_zz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_x_y[i] = -2.0 * g_xy_0_x_y[i] * a_exp + 4.0 * g_xy_zz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_x_z[i] = -2.0 * g_xy_0_x_z[i] * a_exp + 4.0 * g_xy_zz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (210-213)

    #pragma omp simd aligned(g_x_z_0_0_y_z_y_x, g_x_z_0_0_y_z_y_y, g_x_z_0_0_y_z_y_z, g_xy_0_y_x, g_xy_0_y_y, g_xy_0_y_z, g_xy_zz_y_x, g_xy_zz_y_y, g_xy_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_z_y_x[i] = -2.0 * g_xy_0_y_x[i] * a_exp + 4.0 * g_xy_zz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_y_y[i] = -2.0 * g_xy_0_y_y[i] * a_exp + 4.0 * g_xy_zz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_y_z[i] = -2.0 * g_xy_0_y_z[i] * a_exp + 4.0 * g_xy_zz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (213-216)

    #pragma omp simd aligned(g_x_z_0_0_y_z_z_x, g_x_z_0_0_y_z_z_y, g_x_z_0_0_y_z_z_z, g_xy_0_z_x, g_xy_0_z_y, g_xy_0_z_z, g_xy_zz_z_x, g_xy_zz_z_y, g_xy_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_z_z_x[i] = -2.0 * g_xy_0_z_x[i] * a_exp + 4.0 * g_xy_zz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_z_y[i] = -2.0 * g_xy_0_z_y[i] * a_exp + 4.0 * g_xy_zz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_z_z[i] = -2.0 * g_xy_0_z_z[i] * a_exp + 4.0 * g_xy_zz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (216-219)

    #pragma omp simd aligned(g_x_z_0_0_z_x_x_x, g_x_z_0_0_z_x_x_y, g_x_z_0_0_z_x_x_z, g_xz_xz_x_x, g_xz_xz_x_y, g_xz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_x_x_x[i] = 4.0 * g_xz_xz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_x_y[i] = 4.0 * g_xz_xz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_x_z[i] = 4.0 * g_xz_xz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (219-222)

    #pragma omp simd aligned(g_x_z_0_0_z_x_y_x, g_x_z_0_0_z_x_y_y, g_x_z_0_0_z_x_y_z, g_xz_xz_y_x, g_xz_xz_y_y, g_xz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_x_y_x[i] = 4.0 * g_xz_xz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_y_y[i] = 4.0 * g_xz_xz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_y_z[i] = 4.0 * g_xz_xz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (222-225)

    #pragma omp simd aligned(g_x_z_0_0_z_x_z_x, g_x_z_0_0_z_x_z_y, g_x_z_0_0_z_x_z_z, g_xz_xz_z_x, g_xz_xz_z_y, g_xz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_x_z_x[i] = 4.0 * g_xz_xz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_z_y[i] = 4.0 * g_xz_xz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_z_z[i] = 4.0 * g_xz_xz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (225-228)

    #pragma omp simd aligned(g_x_z_0_0_z_y_x_x, g_x_z_0_0_z_y_x_y, g_x_z_0_0_z_y_x_z, g_xz_yz_x_x, g_xz_yz_x_y, g_xz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_y_x_x[i] = 4.0 * g_xz_yz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_x_y[i] = 4.0 * g_xz_yz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_x_z[i] = 4.0 * g_xz_yz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (228-231)

    #pragma omp simd aligned(g_x_z_0_0_z_y_y_x, g_x_z_0_0_z_y_y_y, g_x_z_0_0_z_y_y_z, g_xz_yz_y_x, g_xz_yz_y_y, g_xz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_y_y_x[i] = 4.0 * g_xz_yz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_y_y[i] = 4.0 * g_xz_yz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_y_z[i] = 4.0 * g_xz_yz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (231-234)

    #pragma omp simd aligned(g_x_z_0_0_z_y_z_x, g_x_z_0_0_z_y_z_y, g_x_z_0_0_z_y_z_z, g_xz_yz_z_x, g_xz_yz_z_y, g_xz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_y_z_x[i] = 4.0 * g_xz_yz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_z_y[i] = 4.0 * g_xz_yz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_z_z[i] = 4.0 * g_xz_yz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (234-237)

    #pragma omp simd aligned(g_x_z_0_0_z_z_x_x, g_x_z_0_0_z_z_x_y, g_x_z_0_0_z_z_x_z, g_xz_0_x_x, g_xz_0_x_y, g_xz_0_x_z, g_xz_zz_x_x, g_xz_zz_x_y, g_xz_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_z_x_x[i] = -2.0 * g_xz_0_x_x[i] * a_exp + 4.0 * g_xz_zz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_x_y[i] = -2.0 * g_xz_0_x_y[i] * a_exp + 4.0 * g_xz_zz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_x_z[i] = -2.0 * g_xz_0_x_z[i] * a_exp + 4.0 * g_xz_zz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (237-240)

    #pragma omp simd aligned(g_x_z_0_0_z_z_y_x, g_x_z_0_0_z_z_y_y, g_x_z_0_0_z_z_y_z, g_xz_0_y_x, g_xz_0_y_y, g_xz_0_y_z, g_xz_zz_y_x, g_xz_zz_y_y, g_xz_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_z_y_x[i] = -2.0 * g_xz_0_y_x[i] * a_exp + 4.0 * g_xz_zz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_y_y[i] = -2.0 * g_xz_0_y_y[i] * a_exp + 4.0 * g_xz_zz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_y_z[i] = -2.0 * g_xz_0_y_z[i] * a_exp + 4.0 * g_xz_zz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (240-243)

    #pragma omp simd aligned(g_x_z_0_0_z_z_z_x, g_x_z_0_0_z_z_z_y, g_x_z_0_0_z_z_z_z, g_xz_0_z_x, g_xz_0_z_y, g_xz_0_z_z, g_xz_zz_z_x, g_xz_zz_z_y, g_xz_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_z_z_x[i] = -2.0 * g_xz_0_z_x[i] * a_exp + 4.0 * g_xz_zz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_z_y[i] = -2.0 * g_xz_0_z_y[i] * a_exp + 4.0 * g_xz_zz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_z_z[i] = -2.0 * g_xz_0_z_z[i] * a_exp + 4.0 * g_xz_zz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (243-246)

    #pragma omp simd aligned(g_xy_0_x_x, g_xy_0_x_y, g_xy_0_x_z, g_xy_xx_x_x, g_xy_xx_x_y, g_xy_xx_x_z, g_y_x_0_0_x_x_x_x, g_y_x_0_0_x_x_x_y, g_y_x_0_0_x_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_x_x_x[i] = -2.0 * g_xy_0_x_x[i] * a_exp + 4.0 * g_xy_xx_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_x_y[i] = -2.0 * g_xy_0_x_y[i] * a_exp + 4.0 * g_xy_xx_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_x_z[i] = -2.0 * g_xy_0_x_z[i] * a_exp + 4.0 * g_xy_xx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (246-249)

    #pragma omp simd aligned(g_xy_0_y_x, g_xy_0_y_y, g_xy_0_y_z, g_xy_xx_y_x, g_xy_xx_y_y, g_xy_xx_y_z, g_y_x_0_0_x_x_y_x, g_y_x_0_0_x_x_y_y, g_y_x_0_0_x_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_x_y_x[i] = -2.0 * g_xy_0_y_x[i] * a_exp + 4.0 * g_xy_xx_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_y_y[i] = -2.0 * g_xy_0_y_y[i] * a_exp + 4.0 * g_xy_xx_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_y_z[i] = -2.0 * g_xy_0_y_z[i] * a_exp + 4.0 * g_xy_xx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (249-252)

    #pragma omp simd aligned(g_xy_0_z_x, g_xy_0_z_y, g_xy_0_z_z, g_xy_xx_z_x, g_xy_xx_z_y, g_xy_xx_z_z, g_y_x_0_0_x_x_z_x, g_y_x_0_0_x_x_z_y, g_y_x_0_0_x_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_x_z_x[i] = -2.0 * g_xy_0_z_x[i] * a_exp + 4.0 * g_xy_xx_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_z_y[i] = -2.0 * g_xy_0_z_y[i] * a_exp + 4.0 * g_xy_xx_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_z_z[i] = -2.0 * g_xy_0_z_z[i] * a_exp + 4.0 * g_xy_xx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (252-255)

    #pragma omp simd aligned(g_xy_xy_x_x, g_xy_xy_x_y, g_xy_xy_x_z, g_y_x_0_0_x_y_x_x, g_y_x_0_0_x_y_x_y, g_y_x_0_0_x_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_y_x_x[i] = 4.0 * g_xy_xy_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_x_y[i] = 4.0 * g_xy_xy_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_x_z[i] = 4.0 * g_xy_xy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (255-258)

    #pragma omp simd aligned(g_xy_xy_y_x, g_xy_xy_y_y, g_xy_xy_y_z, g_y_x_0_0_x_y_y_x, g_y_x_0_0_x_y_y_y, g_y_x_0_0_x_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_y_y_x[i] = 4.0 * g_xy_xy_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_y_y[i] = 4.0 * g_xy_xy_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_y_z[i] = 4.0 * g_xy_xy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (258-261)

    #pragma omp simd aligned(g_xy_xy_z_x, g_xy_xy_z_y, g_xy_xy_z_z, g_y_x_0_0_x_y_z_x, g_y_x_0_0_x_y_z_y, g_y_x_0_0_x_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_y_z_x[i] = 4.0 * g_xy_xy_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_z_y[i] = 4.0 * g_xy_xy_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_z_z[i] = 4.0 * g_xy_xy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (261-264)

    #pragma omp simd aligned(g_xy_xz_x_x, g_xy_xz_x_y, g_xy_xz_x_z, g_y_x_0_0_x_z_x_x, g_y_x_0_0_x_z_x_y, g_y_x_0_0_x_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_z_x_x[i] = 4.0 * g_xy_xz_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_x_y[i] = 4.0 * g_xy_xz_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_x_z[i] = 4.0 * g_xy_xz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (264-267)

    #pragma omp simd aligned(g_xy_xz_y_x, g_xy_xz_y_y, g_xy_xz_y_z, g_y_x_0_0_x_z_y_x, g_y_x_0_0_x_z_y_y, g_y_x_0_0_x_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_z_y_x[i] = 4.0 * g_xy_xz_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_y_y[i] = 4.0 * g_xy_xz_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_y_z[i] = 4.0 * g_xy_xz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (267-270)

    #pragma omp simd aligned(g_xy_xz_z_x, g_xy_xz_z_y, g_xy_xz_z_z, g_y_x_0_0_x_z_z_x, g_y_x_0_0_x_z_z_y, g_y_x_0_0_x_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_z_z_x[i] = 4.0 * g_xy_xz_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_z_y[i] = 4.0 * g_xy_xz_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_z_z[i] = 4.0 * g_xy_xz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (270-273)

    #pragma omp simd aligned(g_0_0_x_x, g_0_0_x_y, g_0_0_x_z, g_0_xx_x_x, g_0_xx_x_y, g_0_xx_x_z, g_y_x_0_0_y_x_x_x, g_y_x_0_0_y_x_x_y, g_y_x_0_0_y_x_x_z, g_yy_0_x_x, g_yy_0_x_y, g_yy_0_x_z, g_yy_xx_x_x, g_yy_xx_x_y, g_yy_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_x_x_x[i] = g_0_0_x_x[i] - 2.0 * g_0_xx_x_x[i] * b_exp - 2.0 * g_yy_0_x_x[i] * a_exp + 4.0 * g_yy_xx_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_x_y[i] = g_0_0_x_y[i] - 2.0 * g_0_xx_x_y[i] * b_exp - 2.0 * g_yy_0_x_y[i] * a_exp + 4.0 * g_yy_xx_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_x_z[i] = g_0_0_x_z[i] - 2.0 * g_0_xx_x_z[i] * b_exp - 2.0 * g_yy_0_x_z[i] * a_exp + 4.0 * g_yy_xx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (273-276)

    #pragma omp simd aligned(g_0_0_y_x, g_0_0_y_y, g_0_0_y_z, g_0_xx_y_x, g_0_xx_y_y, g_0_xx_y_z, g_y_x_0_0_y_x_y_x, g_y_x_0_0_y_x_y_y, g_y_x_0_0_y_x_y_z, g_yy_0_y_x, g_yy_0_y_y, g_yy_0_y_z, g_yy_xx_y_x, g_yy_xx_y_y, g_yy_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_x_y_x[i] = g_0_0_y_x[i] - 2.0 * g_0_xx_y_x[i] * b_exp - 2.0 * g_yy_0_y_x[i] * a_exp + 4.0 * g_yy_xx_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_y_y[i] = g_0_0_y_y[i] - 2.0 * g_0_xx_y_y[i] * b_exp - 2.0 * g_yy_0_y_y[i] * a_exp + 4.0 * g_yy_xx_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_y_z[i] = g_0_0_y_z[i] - 2.0 * g_0_xx_y_z[i] * b_exp - 2.0 * g_yy_0_y_z[i] * a_exp + 4.0 * g_yy_xx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (276-279)

    #pragma omp simd aligned(g_0_0_z_x, g_0_0_z_y, g_0_0_z_z, g_0_xx_z_x, g_0_xx_z_y, g_0_xx_z_z, g_y_x_0_0_y_x_z_x, g_y_x_0_0_y_x_z_y, g_y_x_0_0_y_x_z_z, g_yy_0_z_x, g_yy_0_z_y, g_yy_0_z_z, g_yy_xx_z_x, g_yy_xx_z_y, g_yy_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_x_z_x[i] = g_0_0_z_x[i] - 2.0 * g_0_xx_z_x[i] * b_exp - 2.0 * g_yy_0_z_x[i] * a_exp + 4.0 * g_yy_xx_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_z_y[i] = g_0_0_z_y[i] - 2.0 * g_0_xx_z_y[i] * b_exp - 2.0 * g_yy_0_z_y[i] * a_exp + 4.0 * g_yy_xx_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_z_z[i] = g_0_0_z_z[i] - 2.0 * g_0_xx_z_z[i] * b_exp - 2.0 * g_yy_0_z_z[i] * a_exp + 4.0 * g_yy_xx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (279-282)

    #pragma omp simd aligned(g_0_xy_x_x, g_0_xy_x_y, g_0_xy_x_z, g_y_x_0_0_y_y_x_x, g_y_x_0_0_y_y_x_y, g_y_x_0_0_y_y_x_z, g_yy_xy_x_x, g_yy_xy_x_y, g_yy_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_y_x_x[i] = -2.0 * g_0_xy_x_x[i] * b_exp + 4.0 * g_yy_xy_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_x_y[i] = -2.0 * g_0_xy_x_y[i] * b_exp + 4.0 * g_yy_xy_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_x_z[i] = -2.0 * g_0_xy_x_z[i] * b_exp + 4.0 * g_yy_xy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (282-285)

    #pragma omp simd aligned(g_0_xy_y_x, g_0_xy_y_y, g_0_xy_y_z, g_y_x_0_0_y_y_y_x, g_y_x_0_0_y_y_y_y, g_y_x_0_0_y_y_y_z, g_yy_xy_y_x, g_yy_xy_y_y, g_yy_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_y_y_x[i] = -2.0 * g_0_xy_y_x[i] * b_exp + 4.0 * g_yy_xy_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_y_y[i] = -2.0 * g_0_xy_y_y[i] * b_exp + 4.0 * g_yy_xy_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_y_z[i] = -2.0 * g_0_xy_y_z[i] * b_exp + 4.0 * g_yy_xy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (285-288)

    #pragma omp simd aligned(g_0_xy_z_x, g_0_xy_z_y, g_0_xy_z_z, g_y_x_0_0_y_y_z_x, g_y_x_0_0_y_y_z_y, g_y_x_0_0_y_y_z_z, g_yy_xy_z_x, g_yy_xy_z_y, g_yy_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_y_z_x[i] = -2.0 * g_0_xy_z_x[i] * b_exp + 4.0 * g_yy_xy_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_z_y[i] = -2.0 * g_0_xy_z_y[i] * b_exp + 4.0 * g_yy_xy_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_z_z[i] = -2.0 * g_0_xy_z_z[i] * b_exp + 4.0 * g_yy_xy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (288-291)

    #pragma omp simd aligned(g_0_xz_x_x, g_0_xz_x_y, g_0_xz_x_z, g_y_x_0_0_y_z_x_x, g_y_x_0_0_y_z_x_y, g_y_x_0_0_y_z_x_z, g_yy_xz_x_x, g_yy_xz_x_y, g_yy_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_z_x_x[i] = -2.0 * g_0_xz_x_x[i] * b_exp + 4.0 * g_yy_xz_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_x_y[i] = -2.0 * g_0_xz_x_y[i] * b_exp + 4.0 * g_yy_xz_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_x_z[i] = -2.0 * g_0_xz_x_z[i] * b_exp + 4.0 * g_yy_xz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (291-294)

    #pragma omp simd aligned(g_0_xz_y_x, g_0_xz_y_y, g_0_xz_y_z, g_y_x_0_0_y_z_y_x, g_y_x_0_0_y_z_y_y, g_y_x_0_0_y_z_y_z, g_yy_xz_y_x, g_yy_xz_y_y, g_yy_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_z_y_x[i] = -2.0 * g_0_xz_y_x[i] * b_exp + 4.0 * g_yy_xz_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_y_y[i] = -2.0 * g_0_xz_y_y[i] * b_exp + 4.0 * g_yy_xz_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_y_z[i] = -2.0 * g_0_xz_y_z[i] * b_exp + 4.0 * g_yy_xz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (294-297)

    #pragma omp simd aligned(g_0_xz_z_x, g_0_xz_z_y, g_0_xz_z_z, g_y_x_0_0_y_z_z_x, g_y_x_0_0_y_z_z_y, g_y_x_0_0_y_z_z_z, g_yy_xz_z_x, g_yy_xz_z_y, g_yy_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_z_z_x[i] = -2.0 * g_0_xz_z_x[i] * b_exp + 4.0 * g_yy_xz_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_z_y[i] = -2.0 * g_0_xz_z_y[i] * b_exp + 4.0 * g_yy_xz_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_z_z[i] = -2.0 * g_0_xz_z_z[i] * b_exp + 4.0 * g_yy_xz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (297-300)

    #pragma omp simd aligned(g_y_x_0_0_z_x_x_x, g_y_x_0_0_z_x_x_y, g_y_x_0_0_z_x_x_z, g_yz_0_x_x, g_yz_0_x_y, g_yz_0_x_z, g_yz_xx_x_x, g_yz_xx_x_y, g_yz_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_x_x_x[i] = -2.0 * g_yz_0_x_x[i] * a_exp + 4.0 * g_yz_xx_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_x_y[i] = -2.0 * g_yz_0_x_y[i] * a_exp + 4.0 * g_yz_xx_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_x_z[i] = -2.0 * g_yz_0_x_z[i] * a_exp + 4.0 * g_yz_xx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (300-303)

    #pragma omp simd aligned(g_y_x_0_0_z_x_y_x, g_y_x_0_0_z_x_y_y, g_y_x_0_0_z_x_y_z, g_yz_0_y_x, g_yz_0_y_y, g_yz_0_y_z, g_yz_xx_y_x, g_yz_xx_y_y, g_yz_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_x_y_x[i] = -2.0 * g_yz_0_y_x[i] * a_exp + 4.0 * g_yz_xx_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_y_y[i] = -2.0 * g_yz_0_y_y[i] * a_exp + 4.0 * g_yz_xx_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_y_z[i] = -2.0 * g_yz_0_y_z[i] * a_exp + 4.0 * g_yz_xx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (303-306)

    #pragma omp simd aligned(g_y_x_0_0_z_x_z_x, g_y_x_0_0_z_x_z_y, g_y_x_0_0_z_x_z_z, g_yz_0_z_x, g_yz_0_z_y, g_yz_0_z_z, g_yz_xx_z_x, g_yz_xx_z_y, g_yz_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_x_z_x[i] = -2.0 * g_yz_0_z_x[i] * a_exp + 4.0 * g_yz_xx_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_z_y[i] = -2.0 * g_yz_0_z_y[i] * a_exp + 4.0 * g_yz_xx_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_z_z[i] = -2.0 * g_yz_0_z_z[i] * a_exp + 4.0 * g_yz_xx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (306-309)

    #pragma omp simd aligned(g_y_x_0_0_z_y_x_x, g_y_x_0_0_z_y_x_y, g_y_x_0_0_z_y_x_z, g_yz_xy_x_x, g_yz_xy_x_y, g_yz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_y_x_x[i] = 4.0 * g_yz_xy_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_x_y[i] = 4.0 * g_yz_xy_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_x_z[i] = 4.0 * g_yz_xy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (309-312)

    #pragma omp simd aligned(g_y_x_0_0_z_y_y_x, g_y_x_0_0_z_y_y_y, g_y_x_0_0_z_y_y_z, g_yz_xy_y_x, g_yz_xy_y_y, g_yz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_y_y_x[i] = 4.0 * g_yz_xy_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_y_y[i] = 4.0 * g_yz_xy_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_y_z[i] = 4.0 * g_yz_xy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (312-315)

    #pragma omp simd aligned(g_y_x_0_0_z_y_z_x, g_y_x_0_0_z_y_z_y, g_y_x_0_0_z_y_z_z, g_yz_xy_z_x, g_yz_xy_z_y, g_yz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_y_z_x[i] = 4.0 * g_yz_xy_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_z_y[i] = 4.0 * g_yz_xy_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_z_z[i] = 4.0 * g_yz_xy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (315-318)

    #pragma omp simd aligned(g_y_x_0_0_z_z_x_x, g_y_x_0_0_z_z_x_y, g_y_x_0_0_z_z_x_z, g_yz_xz_x_x, g_yz_xz_x_y, g_yz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_z_x_x[i] = 4.0 * g_yz_xz_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_x_y[i] = 4.0 * g_yz_xz_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_x_z[i] = 4.0 * g_yz_xz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (318-321)

    #pragma omp simd aligned(g_y_x_0_0_z_z_y_x, g_y_x_0_0_z_z_y_y, g_y_x_0_0_z_z_y_z, g_yz_xz_y_x, g_yz_xz_y_y, g_yz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_z_y_x[i] = 4.0 * g_yz_xz_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_y_y[i] = 4.0 * g_yz_xz_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_y_z[i] = 4.0 * g_yz_xz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (321-324)

    #pragma omp simd aligned(g_y_x_0_0_z_z_z_x, g_y_x_0_0_z_z_z_y, g_y_x_0_0_z_z_z_z, g_yz_xz_z_x, g_yz_xz_z_y, g_yz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_z_z_x[i] = 4.0 * g_yz_xz_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_z_y[i] = 4.0 * g_yz_xz_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_z_z[i] = 4.0 * g_yz_xz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (324-327)

    #pragma omp simd aligned(g_xy_xy_x_x, g_xy_xy_x_y, g_xy_xy_x_z, g_y_y_0_0_x_x_x_x, g_y_y_0_0_x_x_x_y, g_y_y_0_0_x_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_x_x_x[i] = 4.0 * g_xy_xy_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_x_y[i] = 4.0 * g_xy_xy_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_x_z[i] = 4.0 * g_xy_xy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (327-330)

    #pragma omp simd aligned(g_xy_xy_y_x, g_xy_xy_y_y, g_xy_xy_y_z, g_y_y_0_0_x_x_y_x, g_y_y_0_0_x_x_y_y, g_y_y_0_0_x_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_x_y_x[i] = 4.0 * g_xy_xy_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_y_y[i] = 4.0 * g_xy_xy_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_y_z[i] = 4.0 * g_xy_xy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (330-333)

    #pragma omp simd aligned(g_xy_xy_z_x, g_xy_xy_z_y, g_xy_xy_z_z, g_y_y_0_0_x_x_z_x, g_y_y_0_0_x_x_z_y, g_y_y_0_0_x_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_x_z_x[i] = 4.0 * g_xy_xy_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_z_y[i] = 4.0 * g_xy_xy_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_z_z[i] = 4.0 * g_xy_xy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (333-336)

    #pragma omp simd aligned(g_xy_0_x_x, g_xy_0_x_y, g_xy_0_x_z, g_xy_yy_x_x, g_xy_yy_x_y, g_xy_yy_x_z, g_y_y_0_0_x_y_x_x, g_y_y_0_0_x_y_x_y, g_y_y_0_0_x_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_y_x_x[i] = -2.0 * g_xy_0_x_x[i] * a_exp + 4.0 * g_xy_yy_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_x_y[i] = -2.0 * g_xy_0_x_y[i] * a_exp + 4.0 * g_xy_yy_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_x_z[i] = -2.0 * g_xy_0_x_z[i] * a_exp + 4.0 * g_xy_yy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (336-339)

    #pragma omp simd aligned(g_xy_0_y_x, g_xy_0_y_y, g_xy_0_y_z, g_xy_yy_y_x, g_xy_yy_y_y, g_xy_yy_y_z, g_y_y_0_0_x_y_y_x, g_y_y_0_0_x_y_y_y, g_y_y_0_0_x_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_y_y_x[i] = -2.0 * g_xy_0_y_x[i] * a_exp + 4.0 * g_xy_yy_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_y_y[i] = -2.0 * g_xy_0_y_y[i] * a_exp + 4.0 * g_xy_yy_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_y_z[i] = -2.0 * g_xy_0_y_z[i] * a_exp + 4.0 * g_xy_yy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (339-342)

    #pragma omp simd aligned(g_xy_0_z_x, g_xy_0_z_y, g_xy_0_z_z, g_xy_yy_z_x, g_xy_yy_z_y, g_xy_yy_z_z, g_y_y_0_0_x_y_z_x, g_y_y_0_0_x_y_z_y, g_y_y_0_0_x_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_y_z_x[i] = -2.0 * g_xy_0_z_x[i] * a_exp + 4.0 * g_xy_yy_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_z_y[i] = -2.0 * g_xy_0_z_y[i] * a_exp + 4.0 * g_xy_yy_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_z_z[i] = -2.0 * g_xy_0_z_z[i] * a_exp + 4.0 * g_xy_yy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (342-345)

    #pragma omp simd aligned(g_xy_yz_x_x, g_xy_yz_x_y, g_xy_yz_x_z, g_y_y_0_0_x_z_x_x, g_y_y_0_0_x_z_x_y, g_y_y_0_0_x_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_z_x_x[i] = 4.0 * g_xy_yz_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_x_y[i] = 4.0 * g_xy_yz_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_x_z[i] = 4.0 * g_xy_yz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (345-348)

    #pragma omp simd aligned(g_xy_yz_y_x, g_xy_yz_y_y, g_xy_yz_y_z, g_y_y_0_0_x_z_y_x, g_y_y_0_0_x_z_y_y, g_y_y_0_0_x_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_z_y_x[i] = 4.0 * g_xy_yz_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_y_y[i] = 4.0 * g_xy_yz_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_y_z[i] = 4.0 * g_xy_yz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (348-351)

    #pragma omp simd aligned(g_xy_yz_z_x, g_xy_yz_z_y, g_xy_yz_z_z, g_y_y_0_0_x_z_z_x, g_y_y_0_0_x_z_z_y, g_y_y_0_0_x_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_z_z_x[i] = 4.0 * g_xy_yz_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_z_y[i] = 4.0 * g_xy_yz_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_z_z[i] = 4.0 * g_xy_yz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (351-354)

    #pragma omp simd aligned(g_0_xy_x_x, g_0_xy_x_y, g_0_xy_x_z, g_y_y_0_0_y_x_x_x, g_y_y_0_0_y_x_x_y, g_y_y_0_0_y_x_x_z, g_yy_xy_x_x, g_yy_xy_x_y, g_yy_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_x_x_x[i] = -2.0 * g_0_xy_x_x[i] * b_exp + 4.0 * g_yy_xy_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_x_y[i] = -2.0 * g_0_xy_x_y[i] * b_exp + 4.0 * g_yy_xy_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_x_z[i] = -2.0 * g_0_xy_x_z[i] * b_exp + 4.0 * g_yy_xy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (354-357)

    #pragma omp simd aligned(g_0_xy_y_x, g_0_xy_y_y, g_0_xy_y_z, g_y_y_0_0_y_x_y_x, g_y_y_0_0_y_x_y_y, g_y_y_0_0_y_x_y_z, g_yy_xy_y_x, g_yy_xy_y_y, g_yy_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_x_y_x[i] = -2.0 * g_0_xy_y_x[i] * b_exp + 4.0 * g_yy_xy_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_y_y[i] = -2.0 * g_0_xy_y_y[i] * b_exp + 4.0 * g_yy_xy_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_y_z[i] = -2.0 * g_0_xy_y_z[i] * b_exp + 4.0 * g_yy_xy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (357-360)

    #pragma omp simd aligned(g_0_xy_z_x, g_0_xy_z_y, g_0_xy_z_z, g_y_y_0_0_y_x_z_x, g_y_y_0_0_y_x_z_y, g_y_y_0_0_y_x_z_z, g_yy_xy_z_x, g_yy_xy_z_y, g_yy_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_x_z_x[i] = -2.0 * g_0_xy_z_x[i] * b_exp + 4.0 * g_yy_xy_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_z_y[i] = -2.0 * g_0_xy_z_y[i] * b_exp + 4.0 * g_yy_xy_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_z_z[i] = -2.0 * g_0_xy_z_z[i] * b_exp + 4.0 * g_yy_xy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (360-363)

    #pragma omp simd aligned(g_0_0_x_x, g_0_0_x_y, g_0_0_x_z, g_0_yy_x_x, g_0_yy_x_y, g_0_yy_x_z, g_y_y_0_0_y_y_x_x, g_y_y_0_0_y_y_x_y, g_y_y_0_0_y_y_x_z, g_yy_0_x_x, g_yy_0_x_y, g_yy_0_x_z, g_yy_yy_x_x, g_yy_yy_x_y, g_yy_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_y_x_x[i] = g_0_0_x_x[i] - 2.0 * g_0_yy_x_x[i] * b_exp - 2.0 * g_yy_0_x_x[i] * a_exp + 4.0 * g_yy_yy_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_x_y[i] = g_0_0_x_y[i] - 2.0 * g_0_yy_x_y[i] * b_exp - 2.0 * g_yy_0_x_y[i] * a_exp + 4.0 * g_yy_yy_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_x_z[i] = g_0_0_x_z[i] - 2.0 * g_0_yy_x_z[i] * b_exp - 2.0 * g_yy_0_x_z[i] * a_exp + 4.0 * g_yy_yy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (363-366)

    #pragma omp simd aligned(g_0_0_y_x, g_0_0_y_y, g_0_0_y_z, g_0_yy_y_x, g_0_yy_y_y, g_0_yy_y_z, g_y_y_0_0_y_y_y_x, g_y_y_0_0_y_y_y_y, g_y_y_0_0_y_y_y_z, g_yy_0_y_x, g_yy_0_y_y, g_yy_0_y_z, g_yy_yy_y_x, g_yy_yy_y_y, g_yy_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_y_y_x[i] = g_0_0_y_x[i] - 2.0 * g_0_yy_y_x[i] * b_exp - 2.0 * g_yy_0_y_x[i] * a_exp + 4.0 * g_yy_yy_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_y_y[i] = g_0_0_y_y[i] - 2.0 * g_0_yy_y_y[i] * b_exp - 2.0 * g_yy_0_y_y[i] * a_exp + 4.0 * g_yy_yy_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_y_z[i] = g_0_0_y_z[i] - 2.0 * g_0_yy_y_z[i] * b_exp - 2.0 * g_yy_0_y_z[i] * a_exp + 4.0 * g_yy_yy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (366-369)

    #pragma omp simd aligned(g_0_0_z_x, g_0_0_z_y, g_0_0_z_z, g_0_yy_z_x, g_0_yy_z_y, g_0_yy_z_z, g_y_y_0_0_y_y_z_x, g_y_y_0_0_y_y_z_y, g_y_y_0_0_y_y_z_z, g_yy_0_z_x, g_yy_0_z_y, g_yy_0_z_z, g_yy_yy_z_x, g_yy_yy_z_y, g_yy_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_y_z_x[i] = g_0_0_z_x[i] - 2.0 * g_0_yy_z_x[i] * b_exp - 2.0 * g_yy_0_z_x[i] * a_exp + 4.0 * g_yy_yy_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_z_y[i] = g_0_0_z_y[i] - 2.0 * g_0_yy_z_y[i] * b_exp - 2.0 * g_yy_0_z_y[i] * a_exp + 4.0 * g_yy_yy_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_z_z[i] = g_0_0_z_z[i] - 2.0 * g_0_yy_z_z[i] * b_exp - 2.0 * g_yy_0_z_z[i] * a_exp + 4.0 * g_yy_yy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (369-372)

    #pragma omp simd aligned(g_0_yz_x_x, g_0_yz_x_y, g_0_yz_x_z, g_y_y_0_0_y_z_x_x, g_y_y_0_0_y_z_x_y, g_y_y_0_0_y_z_x_z, g_yy_yz_x_x, g_yy_yz_x_y, g_yy_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_z_x_x[i] = -2.0 * g_0_yz_x_x[i] * b_exp + 4.0 * g_yy_yz_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_x_y[i] = -2.0 * g_0_yz_x_y[i] * b_exp + 4.0 * g_yy_yz_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_x_z[i] = -2.0 * g_0_yz_x_z[i] * b_exp + 4.0 * g_yy_yz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (372-375)

    #pragma omp simd aligned(g_0_yz_y_x, g_0_yz_y_y, g_0_yz_y_z, g_y_y_0_0_y_z_y_x, g_y_y_0_0_y_z_y_y, g_y_y_0_0_y_z_y_z, g_yy_yz_y_x, g_yy_yz_y_y, g_yy_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_z_y_x[i] = -2.0 * g_0_yz_y_x[i] * b_exp + 4.0 * g_yy_yz_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_y_y[i] = -2.0 * g_0_yz_y_y[i] * b_exp + 4.0 * g_yy_yz_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_y_z[i] = -2.0 * g_0_yz_y_z[i] * b_exp + 4.0 * g_yy_yz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (375-378)

    #pragma omp simd aligned(g_0_yz_z_x, g_0_yz_z_y, g_0_yz_z_z, g_y_y_0_0_y_z_z_x, g_y_y_0_0_y_z_z_y, g_y_y_0_0_y_z_z_z, g_yy_yz_z_x, g_yy_yz_z_y, g_yy_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_z_z_x[i] = -2.0 * g_0_yz_z_x[i] * b_exp + 4.0 * g_yy_yz_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_z_y[i] = -2.0 * g_0_yz_z_y[i] * b_exp + 4.0 * g_yy_yz_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_z_z[i] = -2.0 * g_0_yz_z_z[i] * b_exp + 4.0 * g_yy_yz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (378-381)

    #pragma omp simd aligned(g_y_y_0_0_z_x_x_x, g_y_y_0_0_z_x_x_y, g_y_y_0_0_z_x_x_z, g_yz_xy_x_x, g_yz_xy_x_y, g_yz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_x_x_x[i] = 4.0 * g_yz_xy_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_x_y[i] = 4.0 * g_yz_xy_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_x_z[i] = 4.0 * g_yz_xy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (381-384)

    #pragma omp simd aligned(g_y_y_0_0_z_x_y_x, g_y_y_0_0_z_x_y_y, g_y_y_0_0_z_x_y_z, g_yz_xy_y_x, g_yz_xy_y_y, g_yz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_x_y_x[i] = 4.0 * g_yz_xy_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_y_y[i] = 4.0 * g_yz_xy_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_y_z[i] = 4.0 * g_yz_xy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (384-387)

    #pragma omp simd aligned(g_y_y_0_0_z_x_z_x, g_y_y_0_0_z_x_z_y, g_y_y_0_0_z_x_z_z, g_yz_xy_z_x, g_yz_xy_z_y, g_yz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_x_z_x[i] = 4.0 * g_yz_xy_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_z_y[i] = 4.0 * g_yz_xy_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_z_z[i] = 4.0 * g_yz_xy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (387-390)

    #pragma omp simd aligned(g_y_y_0_0_z_y_x_x, g_y_y_0_0_z_y_x_y, g_y_y_0_0_z_y_x_z, g_yz_0_x_x, g_yz_0_x_y, g_yz_0_x_z, g_yz_yy_x_x, g_yz_yy_x_y, g_yz_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_y_x_x[i] = -2.0 * g_yz_0_x_x[i] * a_exp + 4.0 * g_yz_yy_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_x_y[i] = -2.0 * g_yz_0_x_y[i] * a_exp + 4.0 * g_yz_yy_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_x_z[i] = -2.0 * g_yz_0_x_z[i] * a_exp + 4.0 * g_yz_yy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (390-393)

    #pragma omp simd aligned(g_y_y_0_0_z_y_y_x, g_y_y_0_0_z_y_y_y, g_y_y_0_0_z_y_y_z, g_yz_0_y_x, g_yz_0_y_y, g_yz_0_y_z, g_yz_yy_y_x, g_yz_yy_y_y, g_yz_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_y_y_x[i] = -2.0 * g_yz_0_y_x[i] * a_exp + 4.0 * g_yz_yy_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_y_y[i] = -2.0 * g_yz_0_y_y[i] * a_exp + 4.0 * g_yz_yy_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_y_z[i] = -2.0 * g_yz_0_y_z[i] * a_exp + 4.0 * g_yz_yy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (393-396)

    #pragma omp simd aligned(g_y_y_0_0_z_y_z_x, g_y_y_0_0_z_y_z_y, g_y_y_0_0_z_y_z_z, g_yz_0_z_x, g_yz_0_z_y, g_yz_0_z_z, g_yz_yy_z_x, g_yz_yy_z_y, g_yz_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_y_z_x[i] = -2.0 * g_yz_0_z_x[i] * a_exp + 4.0 * g_yz_yy_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_z_y[i] = -2.0 * g_yz_0_z_y[i] * a_exp + 4.0 * g_yz_yy_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_z_z[i] = -2.0 * g_yz_0_z_z[i] * a_exp + 4.0 * g_yz_yy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (396-399)

    #pragma omp simd aligned(g_y_y_0_0_z_z_x_x, g_y_y_0_0_z_z_x_y, g_y_y_0_0_z_z_x_z, g_yz_yz_x_x, g_yz_yz_x_y, g_yz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_z_x_x[i] = 4.0 * g_yz_yz_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_x_y[i] = 4.0 * g_yz_yz_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_x_z[i] = 4.0 * g_yz_yz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (399-402)

    #pragma omp simd aligned(g_y_y_0_0_z_z_y_x, g_y_y_0_0_z_z_y_y, g_y_y_0_0_z_z_y_z, g_yz_yz_y_x, g_yz_yz_y_y, g_yz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_z_y_x[i] = 4.0 * g_yz_yz_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_y_y[i] = 4.0 * g_yz_yz_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_y_z[i] = 4.0 * g_yz_yz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (402-405)

    #pragma omp simd aligned(g_y_y_0_0_z_z_z_x, g_y_y_0_0_z_z_z_y, g_y_y_0_0_z_z_z_z, g_yz_yz_z_x, g_yz_yz_z_y, g_yz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_z_z_x[i] = 4.0 * g_yz_yz_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_z_y[i] = 4.0 * g_yz_yz_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_z_z[i] = 4.0 * g_yz_yz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (405-408)

    #pragma omp simd aligned(g_xy_xz_x_x, g_xy_xz_x_y, g_xy_xz_x_z, g_y_z_0_0_x_x_x_x, g_y_z_0_0_x_x_x_y, g_y_z_0_0_x_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_x_x_x[i] = 4.0 * g_xy_xz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_x_y[i] = 4.0 * g_xy_xz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_x_z[i] = 4.0 * g_xy_xz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (408-411)

    #pragma omp simd aligned(g_xy_xz_y_x, g_xy_xz_y_y, g_xy_xz_y_z, g_y_z_0_0_x_x_y_x, g_y_z_0_0_x_x_y_y, g_y_z_0_0_x_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_x_y_x[i] = 4.0 * g_xy_xz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_y_y[i] = 4.0 * g_xy_xz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_y_z[i] = 4.0 * g_xy_xz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (411-414)

    #pragma omp simd aligned(g_xy_xz_z_x, g_xy_xz_z_y, g_xy_xz_z_z, g_y_z_0_0_x_x_z_x, g_y_z_0_0_x_x_z_y, g_y_z_0_0_x_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_x_z_x[i] = 4.0 * g_xy_xz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_z_y[i] = 4.0 * g_xy_xz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_z_z[i] = 4.0 * g_xy_xz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (414-417)

    #pragma omp simd aligned(g_xy_yz_x_x, g_xy_yz_x_y, g_xy_yz_x_z, g_y_z_0_0_x_y_x_x, g_y_z_0_0_x_y_x_y, g_y_z_0_0_x_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_y_x_x[i] = 4.0 * g_xy_yz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_x_y[i] = 4.0 * g_xy_yz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_x_z[i] = 4.0 * g_xy_yz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (417-420)

    #pragma omp simd aligned(g_xy_yz_y_x, g_xy_yz_y_y, g_xy_yz_y_z, g_y_z_0_0_x_y_y_x, g_y_z_0_0_x_y_y_y, g_y_z_0_0_x_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_y_y_x[i] = 4.0 * g_xy_yz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_y_y[i] = 4.0 * g_xy_yz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_y_z[i] = 4.0 * g_xy_yz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (420-423)

    #pragma omp simd aligned(g_xy_yz_z_x, g_xy_yz_z_y, g_xy_yz_z_z, g_y_z_0_0_x_y_z_x, g_y_z_0_0_x_y_z_y, g_y_z_0_0_x_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_y_z_x[i] = 4.0 * g_xy_yz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_z_y[i] = 4.0 * g_xy_yz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_z_z[i] = 4.0 * g_xy_yz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (423-426)

    #pragma omp simd aligned(g_xy_0_x_x, g_xy_0_x_y, g_xy_0_x_z, g_xy_zz_x_x, g_xy_zz_x_y, g_xy_zz_x_z, g_y_z_0_0_x_z_x_x, g_y_z_0_0_x_z_x_y, g_y_z_0_0_x_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_z_x_x[i] = -2.0 * g_xy_0_x_x[i] * a_exp + 4.0 * g_xy_zz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_x_y[i] = -2.0 * g_xy_0_x_y[i] * a_exp + 4.0 * g_xy_zz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_x_z[i] = -2.0 * g_xy_0_x_z[i] * a_exp + 4.0 * g_xy_zz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (426-429)

    #pragma omp simd aligned(g_xy_0_y_x, g_xy_0_y_y, g_xy_0_y_z, g_xy_zz_y_x, g_xy_zz_y_y, g_xy_zz_y_z, g_y_z_0_0_x_z_y_x, g_y_z_0_0_x_z_y_y, g_y_z_0_0_x_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_z_y_x[i] = -2.0 * g_xy_0_y_x[i] * a_exp + 4.0 * g_xy_zz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_y_y[i] = -2.0 * g_xy_0_y_y[i] * a_exp + 4.0 * g_xy_zz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_y_z[i] = -2.0 * g_xy_0_y_z[i] * a_exp + 4.0 * g_xy_zz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (429-432)

    #pragma omp simd aligned(g_xy_0_z_x, g_xy_0_z_y, g_xy_0_z_z, g_xy_zz_z_x, g_xy_zz_z_y, g_xy_zz_z_z, g_y_z_0_0_x_z_z_x, g_y_z_0_0_x_z_z_y, g_y_z_0_0_x_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_z_z_x[i] = -2.0 * g_xy_0_z_x[i] * a_exp + 4.0 * g_xy_zz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_z_y[i] = -2.0 * g_xy_0_z_y[i] * a_exp + 4.0 * g_xy_zz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_z_z[i] = -2.0 * g_xy_0_z_z[i] * a_exp + 4.0 * g_xy_zz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (432-435)

    #pragma omp simd aligned(g_0_xz_x_x, g_0_xz_x_y, g_0_xz_x_z, g_y_z_0_0_y_x_x_x, g_y_z_0_0_y_x_x_y, g_y_z_0_0_y_x_x_z, g_yy_xz_x_x, g_yy_xz_x_y, g_yy_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_x_x_x[i] = -2.0 * g_0_xz_x_x[i] * b_exp + 4.0 * g_yy_xz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_x_y[i] = -2.0 * g_0_xz_x_y[i] * b_exp + 4.0 * g_yy_xz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_x_z[i] = -2.0 * g_0_xz_x_z[i] * b_exp + 4.0 * g_yy_xz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (435-438)

    #pragma omp simd aligned(g_0_xz_y_x, g_0_xz_y_y, g_0_xz_y_z, g_y_z_0_0_y_x_y_x, g_y_z_0_0_y_x_y_y, g_y_z_0_0_y_x_y_z, g_yy_xz_y_x, g_yy_xz_y_y, g_yy_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_x_y_x[i] = -2.0 * g_0_xz_y_x[i] * b_exp + 4.0 * g_yy_xz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_y_y[i] = -2.0 * g_0_xz_y_y[i] * b_exp + 4.0 * g_yy_xz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_y_z[i] = -2.0 * g_0_xz_y_z[i] * b_exp + 4.0 * g_yy_xz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (438-441)

    #pragma omp simd aligned(g_0_xz_z_x, g_0_xz_z_y, g_0_xz_z_z, g_y_z_0_0_y_x_z_x, g_y_z_0_0_y_x_z_y, g_y_z_0_0_y_x_z_z, g_yy_xz_z_x, g_yy_xz_z_y, g_yy_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_x_z_x[i] = -2.0 * g_0_xz_z_x[i] * b_exp + 4.0 * g_yy_xz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_z_y[i] = -2.0 * g_0_xz_z_y[i] * b_exp + 4.0 * g_yy_xz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_z_z[i] = -2.0 * g_0_xz_z_z[i] * b_exp + 4.0 * g_yy_xz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (441-444)

    #pragma omp simd aligned(g_0_yz_x_x, g_0_yz_x_y, g_0_yz_x_z, g_y_z_0_0_y_y_x_x, g_y_z_0_0_y_y_x_y, g_y_z_0_0_y_y_x_z, g_yy_yz_x_x, g_yy_yz_x_y, g_yy_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_y_x_x[i] = -2.0 * g_0_yz_x_x[i] * b_exp + 4.0 * g_yy_yz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_x_y[i] = -2.0 * g_0_yz_x_y[i] * b_exp + 4.0 * g_yy_yz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_x_z[i] = -2.0 * g_0_yz_x_z[i] * b_exp + 4.0 * g_yy_yz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (444-447)

    #pragma omp simd aligned(g_0_yz_y_x, g_0_yz_y_y, g_0_yz_y_z, g_y_z_0_0_y_y_y_x, g_y_z_0_0_y_y_y_y, g_y_z_0_0_y_y_y_z, g_yy_yz_y_x, g_yy_yz_y_y, g_yy_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_y_y_x[i] = -2.0 * g_0_yz_y_x[i] * b_exp + 4.0 * g_yy_yz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_y_y[i] = -2.0 * g_0_yz_y_y[i] * b_exp + 4.0 * g_yy_yz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_y_z[i] = -2.0 * g_0_yz_y_z[i] * b_exp + 4.0 * g_yy_yz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (447-450)

    #pragma omp simd aligned(g_0_yz_z_x, g_0_yz_z_y, g_0_yz_z_z, g_y_z_0_0_y_y_z_x, g_y_z_0_0_y_y_z_y, g_y_z_0_0_y_y_z_z, g_yy_yz_z_x, g_yy_yz_z_y, g_yy_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_y_z_x[i] = -2.0 * g_0_yz_z_x[i] * b_exp + 4.0 * g_yy_yz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_z_y[i] = -2.0 * g_0_yz_z_y[i] * b_exp + 4.0 * g_yy_yz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_z_z[i] = -2.0 * g_0_yz_z_z[i] * b_exp + 4.0 * g_yy_yz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (450-453)

    #pragma omp simd aligned(g_0_0_x_x, g_0_0_x_y, g_0_0_x_z, g_0_zz_x_x, g_0_zz_x_y, g_0_zz_x_z, g_y_z_0_0_y_z_x_x, g_y_z_0_0_y_z_x_y, g_y_z_0_0_y_z_x_z, g_yy_0_x_x, g_yy_0_x_y, g_yy_0_x_z, g_yy_zz_x_x, g_yy_zz_x_y, g_yy_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_z_x_x[i] = g_0_0_x_x[i] - 2.0 * g_0_zz_x_x[i] * b_exp - 2.0 * g_yy_0_x_x[i] * a_exp + 4.0 * g_yy_zz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_x_y[i] = g_0_0_x_y[i] - 2.0 * g_0_zz_x_y[i] * b_exp - 2.0 * g_yy_0_x_y[i] * a_exp + 4.0 * g_yy_zz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_x_z[i] = g_0_0_x_z[i] - 2.0 * g_0_zz_x_z[i] * b_exp - 2.0 * g_yy_0_x_z[i] * a_exp + 4.0 * g_yy_zz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (453-456)

    #pragma omp simd aligned(g_0_0_y_x, g_0_0_y_y, g_0_0_y_z, g_0_zz_y_x, g_0_zz_y_y, g_0_zz_y_z, g_y_z_0_0_y_z_y_x, g_y_z_0_0_y_z_y_y, g_y_z_0_0_y_z_y_z, g_yy_0_y_x, g_yy_0_y_y, g_yy_0_y_z, g_yy_zz_y_x, g_yy_zz_y_y, g_yy_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_z_y_x[i] = g_0_0_y_x[i] - 2.0 * g_0_zz_y_x[i] * b_exp - 2.0 * g_yy_0_y_x[i] * a_exp + 4.0 * g_yy_zz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_y_y[i] = g_0_0_y_y[i] - 2.0 * g_0_zz_y_y[i] * b_exp - 2.0 * g_yy_0_y_y[i] * a_exp + 4.0 * g_yy_zz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_y_z[i] = g_0_0_y_z[i] - 2.0 * g_0_zz_y_z[i] * b_exp - 2.0 * g_yy_0_y_z[i] * a_exp + 4.0 * g_yy_zz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (456-459)

    #pragma omp simd aligned(g_0_0_z_x, g_0_0_z_y, g_0_0_z_z, g_0_zz_z_x, g_0_zz_z_y, g_0_zz_z_z, g_y_z_0_0_y_z_z_x, g_y_z_0_0_y_z_z_y, g_y_z_0_0_y_z_z_z, g_yy_0_z_x, g_yy_0_z_y, g_yy_0_z_z, g_yy_zz_z_x, g_yy_zz_z_y, g_yy_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_z_z_x[i] = g_0_0_z_x[i] - 2.0 * g_0_zz_z_x[i] * b_exp - 2.0 * g_yy_0_z_x[i] * a_exp + 4.0 * g_yy_zz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_z_y[i] = g_0_0_z_y[i] - 2.0 * g_0_zz_z_y[i] * b_exp - 2.0 * g_yy_0_z_y[i] * a_exp + 4.0 * g_yy_zz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_z_z[i] = g_0_0_z_z[i] - 2.0 * g_0_zz_z_z[i] * b_exp - 2.0 * g_yy_0_z_z[i] * a_exp + 4.0 * g_yy_zz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (459-462)

    #pragma omp simd aligned(g_y_z_0_0_z_x_x_x, g_y_z_0_0_z_x_x_y, g_y_z_0_0_z_x_x_z, g_yz_xz_x_x, g_yz_xz_x_y, g_yz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_x_x_x[i] = 4.0 * g_yz_xz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_x_y[i] = 4.0 * g_yz_xz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_x_z[i] = 4.0 * g_yz_xz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (462-465)

    #pragma omp simd aligned(g_y_z_0_0_z_x_y_x, g_y_z_0_0_z_x_y_y, g_y_z_0_0_z_x_y_z, g_yz_xz_y_x, g_yz_xz_y_y, g_yz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_x_y_x[i] = 4.0 * g_yz_xz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_y_y[i] = 4.0 * g_yz_xz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_y_z[i] = 4.0 * g_yz_xz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (465-468)

    #pragma omp simd aligned(g_y_z_0_0_z_x_z_x, g_y_z_0_0_z_x_z_y, g_y_z_0_0_z_x_z_z, g_yz_xz_z_x, g_yz_xz_z_y, g_yz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_x_z_x[i] = 4.0 * g_yz_xz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_z_y[i] = 4.0 * g_yz_xz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_z_z[i] = 4.0 * g_yz_xz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (468-471)

    #pragma omp simd aligned(g_y_z_0_0_z_y_x_x, g_y_z_0_0_z_y_x_y, g_y_z_0_0_z_y_x_z, g_yz_yz_x_x, g_yz_yz_x_y, g_yz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_y_x_x[i] = 4.0 * g_yz_yz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_x_y[i] = 4.0 * g_yz_yz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_x_z[i] = 4.0 * g_yz_yz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (471-474)

    #pragma omp simd aligned(g_y_z_0_0_z_y_y_x, g_y_z_0_0_z_y_y_y, g_y_z_0_0_z_y_y_z, g_yz_yz_y_x, g_yz_yz_y_y, g_yz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_y_y_x[i] = 4.0 * g_yz_yz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_y_y[i] = 4.0 * g_yz_yz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_y_z[i] = 4.0 * g_yz_yz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (474-477)

    #pragma omp simd aligned(g_y_z_0_0_z_y_z_x, g_y_z_0_0_z_y_z_y, g_y_z_0_0_z_y_z_z, g_yz_yz_z_x, g_yz_yz_z_y, g_yz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_y_z_x[i] = 4.0 * g_yz_yz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_z_y[i] = 4.0 * g_yz_yz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_z_z[i] = 4.0 * g_yz_yz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (477-480)

    #pragma omp simd aligned(g_y_z_0_0_z_z_x_x, g_y_z_0_0_z_z_x_y, g_y_z_0_0_z_z_x_z, g_yz_0_x_x, g_yz_0_x_y, g_yz_0_x_z, g_yz_zz_x_x, g_yz_zz_x_y, g_yz_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_z_x_x[i] = -2.0 * g_yz_0_x_x[i] * a_exp + 4.0 * g_yz_zz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_x_y[i] = -2.0 * g_yz_0_x_y[i] * a_exp + 4.0 * g_yz_zz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_x_z[i] = -2.0 * g_yz_0_x_z[i] * a_exp + 4.0 * g_yz_zz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (480-483)

    #pragma omp simd aligned(g_y_z_0_0_z_z_y_x, g_y_z_0_0_z_z_y_y, g_y_z_0_0_z_z_y_z, g_yz_0_y_x, g_yz_0_y_y, g_yz_0_y_z, g_yz_zz_y_x, g_yz_zz_y_y, g_yz_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_z_y_x[i] = -2.0 * g_yz_0_y_x[i] * a_exp + 4.0 * g_yz_zz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_y_y[i] = -2.0 * g_yz_0_y_y[i] * a_exp + 4.0 * g_yz_zz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_y_z[i] = -2.0 * g_yz_0_y_z[i] * a_exp + 4.0 * g_yz_zz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (483-486)

    #pragma omp simd aligned(g_y_z_0_0_z_z_z_x, g_y_z_0_0_z_z_z_y, g_y_z_0_0_z_z_z_z, g_yz_0_z_x, g_yz_0_z_y, g_yz_0_z_z, g_yz_zz_z_x, g_yz_zz_z_y, g_yz_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_z_z_x[i] = -2.0 * g_yz_0_z_x[i] * a_exp + 4.0 * g_yz_zz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_z_y[i] = -2.0 * g_yz_0_z_y[i] * a_exp + 4.0 * g_yz_zz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_z_z[i] = -2.0 * g_yz_0_z_z[i] * a_exp + 4.0 * g_yz_zz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (486-489)

    #pragma omp simd aligned(g_xz_0_x_x, g_xz_0_x_y, g_xz_0_x_z, g_xz_xx_x_x, g_xz_xx_x_y, g_xz_xx_x_z, g_z_x_0_0_x_x_x_x, g_z_x_0_0_x_x_x_y, g_z_x_0_0_x_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_x_x_x[i] = -2.0 * g_xz_0_x_x[i] * a_exp + 4.0 * g_xz_xx_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_x_y[i] = -2.0 * g_xz_0_x_y[i] * a_exp + 4.0 * g_xz_xx_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_x_z[i] = -2.0 * g_xz_0_x_z[i] * a_exp + 4.0 * g_xz_xx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (489-492)

    #pragma omp simd aligned(g_xz_0_y_x, g_xz_0_y_y, g_xz_0_y_z, g_xz_xx_y_x, g_xz_xx_y_y, g_xz_xx_y_z, g_z_x_0_0_x_x_y_x, g_z_x_0_0_x_x_y_y, g_z_x_0_0_x_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_x_y_x[i] = -2.0 * g_xz_0_y_x[i] * a_exp + 4.0 * g_xz_xx_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_y_y[i] = -2.0 * g_xz_0_y_y[i] * a_exp + 4.0 * g_xz_xx_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_y_z[i] = -2.0 * g_xz_0_y_z[i] * a_exp + 4.0 * g_xz_xx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (492-495)

    #pragma omp simd aligned(g_xz_0_z_x, g_xz_0_z_y, g_xz_0_z_z, g_xz_xx_z_x, g_xz_xx_z_y, g_xz_xx_z_z, g_z_x_0_0_x_x_z_x, g_z_x_0_0_x_x_z_y, g_z_x_0_0_x_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_x_z_x[i] = -2.0 * g_xz_0_z_x[i] * a_exp + 4.0 * g_xz_xx_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_z_y[i] = -2.0 * g_xz_0_z_y[i] * a_exp + 4.0 * g_xz_xx_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_z_z[i] = -2.0 * g_xz_0_z_z[i] * a_exp + 4.0 * g_xz_xx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (495-498)

    #pragma omp simd aligned(g_xz_xy_x_x, g_xz_xy_x_y, g_xz_xy_x_z, g_z_x_0_0_x_y_x_x, g_z_x_0_0_x_y_x_y, g_z_x_0_0_x_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_y_x_x[i] = 4.0 * g_xz_xy_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_x_y[i] = 4.0 * g_xz_xy_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_x_z[i] = 4.0 * g_xz_xy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (498-501)

    #pragma omp simd aligned(g_xz_xy_y_x, g_xz_xy_y_y, g_xz_xy_y_z, g_z_x_0_0_x_y_y_x, g_z_x_0_0_x_y_y_y, g_z_x_0_0_x_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_y_y_x[i] = 4.0 * g_xz_xy_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_y_y[i] = 4.0 * g_xz_xy_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_y_z[i] = 4.0 * g_xz_xy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (501-504)

    #pragma omp simd aligned(g_xz_xy_z_x, g_xz_xy_z_y, g_xz_xy_z_z, g_z_x_0_0_x_y_z_x, g_z_x_0_0_x_y_z_y, g_z_x_0_0_x_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_y_z_x[i] = 4.0 * g_xz_xy_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_z_y[i] = 4.0 * g_xz_xy_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_z_z[i] = 4.0 * g_xz_xy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (504-507)

    #pragma omp simd aligned(g_xz_xz_x_x, g_xz_xz_x_y, g_xz_xz_x_z, g_z_x_0_0_x_z_x_x, g_z_x_0_0_x_z_x_y, g_z_x_0_0_x_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_z_x_x[i] = 4.0 * g_xz_xz_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_x_y[i] = 4.0 * g_xz_xz_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_x_z[i] = 4.0 * g_xz_xz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (507-510)

    #pragma omp simd aligned(g_xz_xz_y_x, g_xz_xz_y_y, g_xz_xz_y_z, g_z_x_0_0_x_z_y_x, g_z_x_0_0_x_z_y_y, g_z_x_0_0_x_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_z_y_x[i] = 4.0 * g_xz_xz_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_y_y[i] = 4.0 * g_xz_xz_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_y_z[i] = 4.0 * g_xz_xz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (510-513)

    #pragma omp simd aligned(g_xz_xz_z_x, g_xz_xz_z_y, g_xz_xz_z_z, g_z_x_0_0_x_z_z_x, g_z_x_0_0_x_z_z_y, g_z_x_0_0_x_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_z_z_x[i] = 4.0 * g_xz_xz_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_z_y[i] = 4.0 * g_xz_xz_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_z_z[i] = 4.0 * g_xz_xz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (513-516)

    #pragma omp simd aligned(g_yz_0_x_x, g_yz_0_x_y, g_yz_0_x_z, g_yz_xx_x_x, g_yz_xx_x_y, g_yz_xx_x_z, g_z_x_0_0_y_x_x_x, g_z_x_0_0_y_x_x_y, g_z_x_0_0_y_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_x_x_x[i] = -2.0 * g_yz_0_x_x[i] * a_exp + 4.0 * g_yz_xx_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_x_y[i] = -2.0 * g_yz_0_x_y[i] * a_exp + 4.0 * g_yz_xx_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_x_z[i] = -2.0 * g_yz_0_x_z[i] * a_exp + 4.0 * g_yz_xx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (516-519)

    #pragma omp simd aligned(g_yz_0_y_x, g_yz_0_y_y, g_yz_0_y_z, g_yz_xx_y_x, g_yz_xx_y_y, g_yz_xx_y_z, g_z_x_0_0_y_x_y_x, g_z_x_0_0_y_x_y_y, g_z_x_0_0_y_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_x_y_x[i] = -2.0 * g_yz_0_y_x[i] * a_exp + 4.0 * g_yz_xx_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_y_y[i] = -2.0 * g_yz_0_y_y[i] * a_exp + 4.0 * g_yz_xx_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_y_z[i] = -2.0 * g_yz_0_y_z[i] * a_exp + 4.0 * g_yz_xx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (519-522)

    #pragma omp simd aligned(g_yz_0_z_x, g_yz_0_z_y, g_yz_0_z_z, g_yz_xx_z_x, g_yz_xx_z_y, g_yz_xx_z_z, g_z_x_0_0_y_x_z_x, g_z_x_0_0_y_x_z_y, g_z_x_0_0_y_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_x_z_x[i] = -2.0 * g_yz_0_z_x[i] * a_exp + 4.0 * g_yz_xx_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_z_y[i] = -2.0 * g_yz_0_z_y[i] * a_exp + 4.0 * g_yz_xx_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_z_z[i] = -2.0 * g_yz_0_z_z[i] * a_exp + 4.0 * g_yz_xx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (522-525)

    #pragma omp simd aligned(g_yz_xy_x_x, g_yz_xy_x_y, g_yz_xy_x_z, g_z_x_0_0_y_y_x_x, g_z_x_0_0_y_y_x_y, g_z_x_0_0_y_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_y_x_x[i] = 4.0 * g_yz_xy_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_x_y[i] = 4.0 * g_yz_xy_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_x_z[i] = 4.0 * g_yz_xy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (525-528)

    #pragma omp simd aligned(g_yz_xy_y_x, g_yz_xy_y_y, g_yz_xy_y_z, g_z_x_0_0_y_y_y_x, g_z_x_0_0_y_y_y_y, g_z_x_0_0_y_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_y_y_x[i] = 4.0 * g_yz_xy_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_y_y[i] = 4.0 * g_yz_xy_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_y_z[i] = 4.0 * g_yz_xy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (528-531)

    #pragma omp simd aligned(g_yz_xy_z_x, g_yz_xy_z_y, g_yz_xy_z_z, g_z_x_0_0_y_y_z_x, g_z_x_0_0_y_y_z_y, g_z_x_0_0_y_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_y_z_x[i] = 4.0 * g_yz_xy_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_z_y[i] = 4.0 * g_yz_xy_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_z_z[i] = 4.0 * g_yz_xy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (531-534)

    #pragma omp simd aligned(g_yz_xz_x_x, g_yz_xz_x_y, g_yz_xz_x_z, g_z_x_0_0_y_z_x_x, g_z_x_0_0_y_z_x_y, g_z_x_0_0_y_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_z_x_x[i] = 4.0 * g_yz_xz_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_x_y[i] = 4.0 * g_yz_xz_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_x_z[i] = 4.0 * g_yz_xz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (534-537)

    #pragma omp simd aligned(g_yz_xz_y_x, g_yz_xz_y_y, g_yz_xz_y_z, g_z_x_0_0_y_z_y_x, g_z_x_0_0_y_z_y_y, g_z_x_0_0_y_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_z_y_x[i] = 4.0 * g_yz_xz_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_y_y[i] = 4.0 * g_yz_xz_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_y_z[i] = 4.0 * g_yz_xz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (537-540)

    #pragma omp simd aligned(g_yz_xz_z_x, g_yz_xz_z_y, g_yz_xz_z_z, g_z_x_0_0_y_z_z_x, g_z_x_0_0_y_z_z_y, g_z_x_0_0_y_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_z_z_x[i] = 4.0 * g_yz_xz_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_z_y[i] = 4.0 * g_yz_xz_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_z_z[i] = 4.0 * g_yz_xz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (540-543)

    #pragma omp simd aligned(g_0_0_x_x, g_0_0_x_y, g_0_0_x_z, g_0_xx_x_x, g_0_xx_x_y, g_0_xx_x_z, g_z_x_0_0_z_x_x_x, g_z_x_0_0_z_x_x_y, g_z_x_0_0_z_x_x_z, g_zz_0_x_x, g_zz_0_x_y, g_zz_0_x_z, g_zz_xx_x_x, g_zz_xx_x_y, g_zz_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_x_x_x[i] = g_0_0_x_x[i] - 2.0 * g_0_xx_x_x[i] * b_exp - 2.0 * g_zz_0_x_x[i] * a_exp + 4.0 * g_zz_xx_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_x_y[i] = g_0_0_x_y[i] - 2.0 * g_0_xx_x_y[i] * b_exp - 2.0 * g_zz_0_x_y[i] * a_exp + 4.0 * g_zz_xx_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_x_z[i] = g_0_0_x_z[i] - 2.0 * g_0_xx_x_z[i] * b_exp - 2.0 * g_zz_0_x_z[i] * a_exp + 4.0 * g_zz_xx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (543-546)

    #pragma omp simd aligned(g_0_0_y_x, g_0_0_y_y, g_0_0_y_z, g_0_xx_y_x, g_0_xx_y_y, g_0_xx_y_z, g_z_x_0_0_z_x_y_x, g_z_x_0_0_z_x_y_y, g_z_x_0_0_z_x_y_z, g_zz_0_y_x, g_zz_0_y_y, g_zz_0_y_z, g_zz_xx_y_x, g_zz_xx_y_y, g_zz_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_x_y_x[i] = g_0_0_y_x[i] - 2.0 * g_0_xx_y_x[i] * b_exp - 2.0 * g_zz_0_y_x[i] * a_exp + 4.0 * g_zz_xx_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_y_y[i] = g_0_0_y_y[i] - 2.0 * g_0_xx_y_y[i] * b_exp - 2.0 * g_zz_0_y_y[i] * a_exp + 4.0 * g_zz_xx_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_y_z[i] = g_0_0_y_z[i] - 2.0 * g_0_xx_y_z[i] * b_exp - 2.0 * g_zz_0_y_z[i] * a_exp + 4.0 * g_zz_xx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (546-549)

    #pragma omp simd aligned(g_0_0_z_x, g_0_0_z_y, g_0_0_z_z, g_0_xx_z_x, g_0_xx_z_y, g_0_xx_z_z, g_z_x_0_0_z_x_z_x, g_z_x_0_0_z_x_z_y, g_z_x_0_0_z_x_z_z, g_zz_0_z_x, g_zz_0_z_y, g_zz_0_z_z, g_zz_xx_z_x, g_zz_xx_z_y, g_zz_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_x_z_x[i] = g_0_0_z_x[i] - 2.0 * g_0_xx_z_x[i] * b_exp - 2.0 * g_zz_0_z_x[i] * a_exp + 4.0 * g_zz_xx_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_z_y[i] = g_0_0_z_y[i] - 2.0 * g_0_xx_z_y[i] * b_exp - 2.0 * g_zz_0_z_y[i] * a_exp + 4.0 * g_zz_xx_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_z_z[i] = g_0_0_z_z[i] - 2.0 * g_0_xx_z_z[i] * b_exp - 2.0 * g_zz_0_z_z[i] * a_exp + 4.0 * g_zz_xx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (549-552)

    #pragma omp simd aligned(g_0_xy_x_x, g_0_xy_x_y, g_0_xy_x_z, g_z_x_0_0_z_y_x_x, g_z_x_0_0_z_y_x_y, g_z_x_0_0_z_y_x_z, g_zz_xy_x_x, g_zz_xy_x_y, g_zz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_y_x_x[i] = -2.0 * g_0_xy_x_x[i] * b_exp + 4.0 * g_zz_xy_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_x_y[i] = -2.0 * g_0_xy_x_y[i] * b_exp + 4.0 * g_zz_xy_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_x_z[i] = -2.0 * g_0_xy_x_z[i] * b_exp + 4.0 * g_zz_xy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (552-555)

    #pragma omp simd aligned(g_0_xy_y_x, g_0_xy_y_y, g_0_xy_y_z, g_z_x_0_0_z_y_y_x, g_z_x_0_0_z_y_y_y, g_z_x_0_0_z_y_y_z, g_zz_xy_y_x, g_zz_xy_y_y, g_zz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_y_y_x[i] = -2.0 * g_0_xy_y_x[i] * b_exp + 4.0 * g_zz_xy_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_y_y[i] = -2.0 * g_0_xy_y_y[i] * b_exp + 4.0 * g_zz_xy_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_y_z[i] = -2.0 * g_0_xy_y_z[i] * b_exp + 4.0 * g_zz_xy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (555-558)

    #pragma omp simd aligned(g_0_xy_z_x, g_0_xy_z_y, g_0_xy_z_z, g_z_x_0_0_z_y_z_x, g_z_x_0_0_z_y_z_y, g_z_x_0_0_z_y_z_z, g_zz_xy_z_x, g_zz_xy_z_y, g_zz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_y_z_x[i] = -2.0 * g_0_xy_z_x[i] * b_exp + 4.0 * g_zz_xy_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_z_y[i] = -2.0 * g_0_xy_z_y[i] * b_exp + 4.0 * g_zz_xy_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_z_z[i] = -2.0 * g_0_xy_z_z[i] * b_exp + 4.0 * g_zz_xy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (558-561)

    #pragma omp simd aligned(g_0_xz_x_x, g_0_xz_x_y, g_0_xz_x_z, g_z_x_0_0_z_z_x_x, g_z_x_0_0_z_z_x_y, g_z_x_0_0_z_z_x_z, g_zz_xz_x_x, g_zz_xz_x_y, g_zz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_z_x_x[i] = -2.0 * g_0_xz_x_x[i] * b_exp + 4.0 * g_zz_xz_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_x_y[i] = -2.0 * g_0_xz_x_y[i] * b_exp + 4.0 * g_zz_xz_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_x_z[i] = -2.0 * g_0_xz_x_z[i] * b_exp + 4.0 * g_zz_xz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (561-564)

    #pragma omp simd aligned(g_0_xz_y_x, g_0_xz_y_y, g_0_xz_y_z, g_z_x_0_0_z_z_y_x, g_z_x_0_0_z_z_y_y, g_z_x_0_0_z_z_y_z, g_zz_xz_y_x, g_zz_xz_y_y, g_zz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_z_y_x[i] = -2.0 * g_0_xz_y_x[i] * b_exp + 4.0 * g_zz_xz_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_y_y[i] = -2.0 * g_0_xz_y_y[i] * b_exp + 4.0 * g_zz_xz_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_y_z[i] = -2.0 * g_0_xz_y_z[i] * b_exp + 4.0 * g_zz_xz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (564-567)

    #pragma omp simd aligned(g_0_xz_z_x, g_0_xz_z_y, g_0_xz_z_z, g_z_x_0_0_z_z_z_x, g_z_x_0_0_z_z_z_y, g_z_x_0_0_z_z_z_z, g_zz_xz_z_x, g_zz_xz_z_y, g_zz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_z_z_x[i] = -2.0 * g_0_xz_z_x[i] * b_exp + 4.0 * g_zz_xz_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_z_y[i] = -2.0 * g_0_xz_z_y[i] * b_exp + 4.0 * g_zz_xz_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_z_z[i] = -2.0 * g_0_xz_z_z[i] * b_exp + 4.0 * g_zz_xz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (567-570)

    #pragma omp simd aligned(g_xz_xy_x_x, g_xz_xy_x_y, g_xz_xy_x_z, g_z_y_0_0_x_x_x_x, g_z_y_0_0_x_x_x_y, g_z_y_0_0_x_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_x_x_x[i] = 4.0 * g_xz_xy_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_x_y[i] = 4.0 * g_xz_xy_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_x_z[i] = 4.0 * g_xz_xy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (570-573)

    #pragma omp simd aligned(g_xz_xy_y_x, g_xz_xy_y_y, g_xz_xy_y_z, g_z_y_0_0_x_x_y_x, g_z_y_0_0_x_x_y_y, g_z_y_0_0_x_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_x_y_x[i] = 4.0 * g_xz_xy_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_y_y[i] = 4.0 * g_xz_xy_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_y_z[i] = 4.0 * g_xz_xy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (573-576)

    #pragma omp simd aligned(g_xz_xy_z_x, g_xz_xy_z_y, g_xz_xy_z_z, g_z_y_0_0_x_x_z_x, g_z_y_0_0_x_x_z_y, g_z_y_0_0_x_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_x_z_x[i] = 4.0 * g_xz_xy_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_z_y[i] = 4.0 * g_xz_xy_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_z_z[i] = 4.0 * g_xz_xy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (576-579)

    #pragma omp simd aligned(g_xz_0_x_x, g_xz_0_x_y, g_xz_0_x_z, g_xz_yy_x_x, g_xz_yy_x_y, g_xz_yy_x_z, g_z_y_0_0_x_y_x_x, g_z_y_0_0_x_y_x_y, g_z_y_0_0_x_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_y_x_x[i] = -2.0 * g_xz_0_x_x[i] * a_exp + 4.0 * g_xz_yy_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_x_y[i] = -2.0 * g_xz_0_x_y[i] * a_exp + 4.0 * g_xz_yy_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_x_z[i] = -2.0 * g_xz_0_x_z[i] * a_exp + 4.0 * g_xz_yy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (579-582)

    #pragma omp simd aligned(g_xz_0_y_x, g_xz_0_y_y, g_xz_0_y_z, g_xz_yy_y_x, g_xz_yy_y_y, g_xz_yy_y_z, g_z_y_0_0_x_y_y_x, g_z_y_0_0_x_y_y_y, g_z_y_0_0_x_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_y_y_x[i] = -2.0 * g_xz_0_y_x[i] * a_exp + 4.0 * g_xz_yy_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_y_y[i] = -2.0 * g_xz_0_y_y[i] * a_exp + 4.0 * g_xz_yy_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_y_z[i] = -2.0 * g_xz_0_y_z[i] * a_exp + 4.0 * g_xz_yy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (582-585)

    #pragma omp simd aligned(g_xz_0_z_x, g_xz_0_z_y, g_xz_0_z_z, g_xz_yy_z_x, g_xz_yy_z_y, g_xz_yy_z_z, g_z_y_0_0_x_y_z_x, g_z_y_0_0_x_y_z_y, g_z_y_0_0_x_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_y_z_x[i] = -2.0 * g_xz_0_z_x[i] * a_exp + 4.0 * g_xz_yy_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_z_y[i] = -2.0 * g_xz_0_z_y[i] * a_exp + 4.0 * g_xz_yy_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_z_z[i] = -2.0 * g_xz_0_z_z[i] * a_exp + 4.0 * g_xz_yy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (585-588)

    #pragma omp simd aligned(g_xz_yz_x_x, g_xz_yz_x_y, g_xz_yz_x_z, g_z_y_0_0_x_z_x_x, g_z_y_0_0_x_z_x_y, g_z_y_0_0_x_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_z_x_x[i] = 4.0 * g_xz_yz_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_x_y[i] = 4.0 * g_xz_yz_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_x_z[i] = 4.0 * g_xz_yz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (588-591)

    #pragma omp simd aligned(g_xz_yz_y_x, g_xz_yz_y_y, g_xz_yz_y_z, g_z_y_0_0_x_z_y_x, g_z_y_0_0_x_z_y_y, g_z_y_0_0_x_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_z_y_x[i] = 4.0 * g_xz_yz_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_y_y[i] = 4.0 * g_xz_yz_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_y_z[i] = 4.0 * g_xz_yz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (591-594)

    #pragma omp simd aligned(g_xz_yz_z_x, g_xz_yz_z_y, g_xz_yz_z_z, g_z_y_0_0_x_z_z_x, g_z_y_0_0_x_z_z_y, g_z_y_0_0_x_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_z_z_x[i] = 4.0 * g_xz_yz_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_z_y[i] = 4.0 * g_xz_yz_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_z_z[i] = 4.0 * g_xz_yz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (594-597)

    #pragma omp simd aligned(g_yz_xy_x_x, g_yz_xy_x_y, g_yz_xy_x_z, g_z_y_0_0_y_x_x_x, g_z_y_0_0_y_x_x_y, g_z_y_0_0_y_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_x_x_x[i] = 4.0 * g_yz_xy_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_x_y[i] = 4.0 * g_yz_xy_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_x_z[i] = 4.0 * g_yz_xy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (597-600)

    #pragma omp simd aligned(g_yz_xy_y_x, g_yz_xy_y_y, g_yz_xy_y_z, g_z_y_0_0_y_x_y_x, g_z_y_0_0_y_x_y_y, g_z_y_0_0_y_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_x_y_x[i] = 4.0 * g_yz_xy_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_y_y[i] = 4.0 * g_yz_xy_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_y_z[i] = 4.0 * g_yz_xy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (600-603)

    #pragma omp simd aligned(g_yz_xy_z_x, g_yz_xy_z_y, g_yz_xy_z_z, g_z_y_0_0_y_x_z_x, g_z_y_0_0_y_x_z_y, g_z_y_0_0_y_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_x_z_x[i] = 4.0 * g_yz_xy_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_z_y[i] = 4.0 * g_yz_xy_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_z_z[i] = 4.0 * g_yz_xy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (603-606)

    #pragma omp simd aligned(g_yz_0_x_x, g_yz_0_x_y, g_yz_0_x_z, g_yz_yy_x_x, g_yz_yy_x_y, g_yz_yy_x_z, g_z_y_0_0_y_y_x_x, g_z_y_0_0_y_y_x_y, g_z_y_0_0_y_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_y_x_x[i] = -2.0 * g_yz_0_x_x[i] * a_exp + 4.0 * g_yz_yy_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_x_y[i] = -2.0 * g_yz_0_x_y[i] * a_exp + 4.0 * g_yz_yy_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_x_z[i] = -2.0 * g_yz_0_x_z[i] * a_exp + 4.0 * g_yz_yy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (606-609)

    #pragma omp simd aligned(g_yz_0_y_x, g_yz_0_y_y, g_yz_0_y_z, g_yz_yy_y_x, g_yz_yy_y_y, g_yz_yy_y_z, g_z_y_0_0_y_y_y_x, g_z_y_0_0_y_y_y_y, g_z_y_0_0_y_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_y_y_x[i] = -2.0 * g_yz_0_y_x[i] * a_exp + 4.0 * g_yz_yy_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_y_y[i] = -2.0 * g_yz_0_y_y[i] * a_exp + 4.0 * g_yz_yy_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_y_z[i] = -2.0 * g_yz_0_y_z[i] * a_exp + 4.0 * g_yz_yy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (609-612)

    #pragma omp simd aligned(g_yz_0_z_x, g_yz_0_z_y, g_yz_0_z_z, g_yz_yy_z_x, g_yz_yy_z_y, g_yz_yy_z_z, g_z_y_0_0_y_y_z_x, g_z_y_0_0_y_y_z_y, g_z_y_0_0_y_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_y_z_x[i] = -2.0 * g_yz_0_z_x[i] * a_exp + 4.0 * g_yz_yy_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_z_y[i] = -2.0 * g_yz_0_z_y[i] * a_exp + 4.0 * g_yz_yy_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_z_z[i] = -2.0 * g_yz_0_z_z[i] * a_exp + 4.0 * g_yz_yy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (612-615)

    #pragma omp simd aligned(g_yz_yz_x_x, g_yz_yz_x_y, g_yz_yz_x_z, g_z_y_0_0_y_z_x_x, g_z_y_0_0_y_z_x_y, g_z_y_0_0_y_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_z_x_x[i] = 4.0 * g_yz_yz_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_x_y[i] = 4.0 * g_yz_yz_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_x_z[i] = 4.0 * g_yz_yz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (615-618)

    #pragma omp simd aligned(g_yz_yz_y_x, g_yz_yz_y_y, g_yz_yz_y_z, g_z_y_0_0_y_z_y_x, g_z_y_0_0_y_z_y_y, g_z_y_0_0_y_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_z_y_x[i] = 4.0 * g_yz_yz_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_y_y[i] = 4.0 * g_yz_yz_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_y_z[i] = 4.0 * g_yz_yz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (618-621)

    #pragma omp simd aligned(g_yz_yz_z_x, g_yz_yz_z_y, g_yz_yz_z_z, g_z_y_0_0_y_z_z_x, g_z_y_0_0_y_z_z_y, g_z_y_0_0_y_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_z_z_x[i] = 4.0 * g_yz_yz_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_z_y[i] = 4.0 * g_yz_yz_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_z_z[i] = 4.0 * g_yz_yz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (621-624)

    #pragma omp simd aligned(g_0_xy_x_x, g_0_xy_x_y, g_0_xy_x_z, g_z_y_0_0_z_x_x_x, g_z_y_0_0_z_x_x_y, g_z_y_0_0_z_x_x_z, g_zz_xy_x_x, g_zz_xy_x_y, g_zz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_x_x_x[i] = -2.0 * g_0_xy_x_x[i] * b_exp + 4.0 * g_zz_xy_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_x_y[i] = -2.0 * g_0_xy_x_y[i] * b_exp + 4.0 * g_zz_xy_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_x_z[i] = -2.0 * g_0_xy_x_z[i] * b_exp + 4.0 * g_zz_xy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (624-627)

    #pragma omp simd aligned(g_0_xy_y_x, g_0_xy_y_y, g_0_xy_y_z, g_z_y_0_0_z_x_y_x, g_z_y_0_0_z_x_y_y, g_z_y_0_0_z_x_y_z, g_zz_xy_y_x, g_zz_xy_y_y, g_zz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_x_y_x[i] = -2.0 * g_0_xy_y_x[i] * b_exp + 4.0 * g_zz_xy_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_y_y[i] = -2.0 * g_0_xy_y_y[i] * b_exp + 4.0 * g_zz_xy_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_y_z[i] = -2.0 * g_0_xy_y_z[i] * b_exp + 4.0 * g_zz_xy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (627-630)

    #pragma omp simd aligned(g_0_xy_z_x, g_0_xy_z_y, g_0_xy_z_z, g_z_y_0_0_z_x_z_x, g_z_y_0_0_z_x_z_y, g_z_y_0_0_z_x_z_z, g_zz_xy_z_x, g_zz_xy_z_y, g_zz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_x_z_x[i] = -2.0 * g_0_xy_z_x[i] * b_exp + 4.0 * g_zz_xy_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_z_y[i] = -2.0 * g_0_xy_z_y[i] * b_exp + 4.0 * g_zz_xy_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_z_z[i] = -2.0 * g_0_xy_z_z[i] * b_exp + 4.0 * g_zz_xy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (630-633)

    #pragma omp simd aligned(g_0_0_x_x, g_0_0_x_y, g_0_0_x_z, g_0_yy_x_x, g_0_yy_x_y, g_0_yy_x_z, g_z_y_0_0_z_y_x_x, g_z_y_0_0_z_y_x_y, g_z_y_0_0_z_y_x_z, g_zz_0_x_x, g_zz_0_x_y, g_zz_0_x_z, g_zz_yy_x_x, g_zz_yy_x_y, g_zz_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_y_x_x[i] = g_0_0_x_x[i] - 2.0 * g_0_yy_x_x[i] * b_exp - 2.0 * g_zz_0_x_x[i] * a_exp + 4.0 * g_zz_yy_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_x_y[i] = g_0_0_x_y[i] - 2.0 * g_0_yy_x_y[i] * b_exp - 2.0 * g_zz_0_x_y[i] * a_exp + 4.0 * g_zz_yy_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_x_z[i] = g_0_0_x_z[i] - 2.0 * g_0_yy_x_z[i] * b_exp - 2.0 * g_zz_0_x_z[i] * a_exp + 4.0 * g_zz_yy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (633-636)

    #pragma omp simd aligned(g_0_0_y_x, g_0_0_y_y, g_0_0_y_z, g_0_yy_y_x, g_0_yy_y_y, g_0_yy_y_z, g_z_y_0_0_z_y_y_x, g_z_y_0_0_z_y_y_y, g_z_y_0_0_z_y_y_z, g_zz_0_y_x, g_zz_0_y_y, g_zz_0_y_z, g_zz_yy_y_x, g_zz_yy_y_y, g_zz_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_y_y_x[i] = g_0_0_y_x[i] - 2.0 * g_0_yy_y_x[i] * b_exp - 2.0 * g_zz_0_y_x[i] * a_exp + 4.0 * g_zz_yy_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_y_y[i] = g_0_0_y_y[i] - 2.0 * g_0_yy_y_y[i] * b_exp - 2.0 * g_zz_0_y_y[i] * a_exp + 4.0 * g_zz_yy_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_y_z[i] = g_0_0_y_z[i] - 2.0 * g_0_yy_y_z[i] * b_exp - 2.0 * g_zz_0_y_z[i] * a_exp + 4.0 * g_zz_yy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (636-639)

    #pragma omp simd aligned(g_0_0_z_x, g_0_0_z_y, g_0_0_z_z, g_0_yy_z_x, g_0_yy_z_y, g_0_yy_z_z, g_z_y_0_0_z_y_z_x, g_z_y_0_0_z_y_z_y, g_z_y_0_0_z_y_z_z, g_zz_0_z_x, g_zz_0_z_y, g_zz_0_z_z, g_zz_yy_z_x, g_zz_yy_z_y, g_zz_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_y_z_x[i] = g_0_0_z_x[i] - 2.0 * g_0_yy_z_x[i] * b_exp - 2.0 * g_zz_0_z_x[i] * a_exp + 4.0 * g_zz_yy_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_z_y[i] = g_0_0_z_y[i] - 2.0 * g_0_yy_z_y[i] * b_exp - 2.0 * g_zz_0_z_y[i] * a_exp + 4.0 * g_zz_yy_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_z_z[i] = g_0_0_z_z[i] - 2.0 * g_0_yy_z_z[i] * b_exp - 2.0 * g_zz_0_z_z[i] * a_exp + 4.0 * g_zz_yy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (639-642)

    #pragma omp simd aligned(g_0_yz_x_x, g_0_yz_x_y, g_0_yz_x_z, g_z_y_0_0_z_z_x_x, g_z_y_0_0_z_z_x_y, g_z_y_0_0_z_z_x_z, g_zz_yz_x_x, g_zz_yz_x_y, g_zz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_z_x_x[i] = -2.0 * g_0_yz_x_x[i] * b_exp + 4.0 * g_zz_yz_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_x_y[i] = -2.0 * g_0_yz_x_y[i] * b_exp + 4.0 * g_zz_yz_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_x_z[i] = -2.0 * g_0_yz_x_z[i] * b_exp + 4.0 * g_zz_yz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (642-645)

    #pragma omp simd aligned(g_0_yz_y_x, g_0_yz_y_y, g_0_yz_y_z, g_z_y_0_0_z_z_y_x, g_z_y_0_0_z_z_y_y, g_z_y_0_0_z_z_y_z, g_zz_yz_y_x, g_zz_yz_y_y, g_zz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_z_y_x[i] = -2.0 * g_0_yz_y_x[i] * b_exp + 4.0 * g_zz_yz_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_y_y[i] = -2.0 * g_0_yz_y_y[i] * b_exp + 4.0 * g_zz_yz_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_y_z[i] = -2.0 * g_0_yz_y_z[i] * b_exp + 4.0 * g_zz_yz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (645-648)

    #pragma omp simd aligned(g_0_yz_z_x, g_0_yz_z_y, g_0_yz_z_z, g_z_y_0_0_z_z_z_x, g_z_y_0_0_z_z_z_y, g_z_y_0_0_z_z_z_z, g_zz_yz_z_x, g_zz_yz_z_y, g_zz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_z_z_x[i] = -2.0 * g_0_yz_z_x[i] * b_exp + 4.0 * g_zz_yz_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_z_y[i] = -2.0 * g_0_yz_z_y[i] * b_exp + 4.0 * g_zz_yz_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_z_z[i] = -2.0 * g_0_yz_z_z[i] * b_exp + 4.0 * g_zz_yz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (648-651)

    #pragma omp simd aligned(g_xz_xz_x_x, g_xz_xz_x_y, g_xz_xz_x_z, g_z_z_0_0_x_x_x_x, g_z_z_0_0_x_x_x_y, g_z_z_0_0_x_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_x_x_x[i] = 4.0 * g_xz_xz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_x_y[i] = 4.0 * g_xz_xz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_x_z[i] = 4.0 * g_xz_xz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (651-654)

    #pragma omp simd aligned(g_xz_xz_y_x, g_xz_xz_y_y, g_xz_xz_y_z, g_z_z_0_0_x_x_y_x, g_z_z_0_0_x_x_y_y, g_z_z_0_0_x_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_x_y_x[i] = 4.0 * g_xz_xz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_y_y[i] = 4.0 * g_xz_xz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_y_z[i] = 4.0 * g_xz_xz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (654-657)

    #pragma omp simd aligned(g_xz_xz_z_x, g_xz_xz_z_y, g_xz_xz_z_z, g_z_z_0_0_x_x_z_x, g_z_z_0_0_x_x_z_y, g_z_z_0_0_x_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_x_z_x[i] = 4.0 * g_xz_xz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_z_y[i] = 4.0 * g_xz_xz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_z_z[i] = 4.0 * g_xz_xz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (657-660)

    #pragma omp simd aligned(g_xz_yz_x_x, g_xz_yz_x_y, g_xz_yz_x_z, g_z_z_0_0_x_y_x_x, g_z_z_0_0_x_y_x_y, g_z_z_0_0_x_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_y_x_x[i] = 4.0 * g_xz_yz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_x_y[i] = 4.0 * g_xz_yz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_x_z[i] = 4.0 * g_xz_yz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (660-663)

    #pragma omp simd aligned(g_xz_yz_y_x, g_xz_yz_y_y, g_xz_yz_y_z, g_z_z_0_0_x_y_y_x, g_z_z_0_0_x_y_y_y, g_z_z_0_0_x_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_y_y_x[i] = 4.0 * g_xz_yz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_y_y[i] = 4.0 * g_xz_yz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_y_z[i] = 4.0 * g_xz_yz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (663-666)

    #pragma omp simd aligned(g_xz_yz_z_x, g_xz_yz_z_y, g_xz_yz_z_z, g_z_z_0_0_x_y_z_x, g_z_z_0_0_x_y_z_y, g_z_z_0_0_x_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_y_z_x[i] = 4.0 * g_xz_yz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_z_y[i] = 4.0 * g_xz_yz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_z_z[i] = 4.0 * g_xz_yz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (666-669)

    #pragma omp simd aligned(g_xz_0_x_x, g_xz_0_x_y, g_xz_0_x_z, g_xz_zz_x_x, g_xz_zz_x_y, g_xz_zz_x_z, g_z_z_0_0_x_z_x_x, g_z_z_0_0_x_z_x_y, g_z_z_0_0_x_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_z_x_x[i] = -2.0 * g_xz_0_x_x[i] * a_exp + 4.0 * g_xz_zz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_x_y[i] = -2.0 * g_xz_0_x_y[i] * a_exp + 4.0 * g_xz_zz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_x_z[i] = -2.0 * g_xz_0_x_z[i] * a_exp + 4.0 * g_xz_zz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (669-672)

    #pragma omp simd aligned(g_xz_0_y_x, g_xz_0_y_y, g_xz_0_y_z, g_xz_zz_y_x, g_xz_zz_y_y, g_xz_zz_y_z, g_z_z_0_0_x_z_y_x, g_z_z_0_0_x_z_y_y, g_z_z_0_0_x_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_z_y_x[i] = -2.0 * g_xz_0_y_x[i] * a_exp + 4.0 * g_xz_zz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_y_y[i] = -2.0 * g_xz_0_y_y[i] * a_exp + 4.0 * g_xz_zz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_y_z[i] = -2.0 * g_xz_0_y_z[i] * a_exp + 4.0 * g_xz_zz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (672-675)

    #pragma omp simd aligned(g_xz_0_z_x, g_xz_0_z_y, g_xz_0_z_z, g_xz_zz_z_x, g_xz_zz_z_y, g_xz_zz_z_z, g_z_z_0_0_x_z_z_x, g_z_z_0_0_x_z_z_y, g_z_z_0_0_x_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_z_z_x[i] = -2.0 * g_xz_0_z_x[i] * a_exp + 4.0 * g_xz_zz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_z_y[i] = -2.0 * g_xz_0_z_y[i] * a_exp + 4.0 * g_xz_zz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_z_z[i] = -2.0 * g_xz_0_z_z[i] * a_exp + 4.0 * g_xz_zz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (675-678)

    #pragma omp simd aligned(g_yz_xz_x_x, g_yz_xz_x_y, g_yz_xz_x_z, g_z_z_0_0_y_x_x_x, g_z_z_0_0_y_x_x_y, g_z_z_0_0_y_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_x_x_x[i] = 4.0 * g_yz_xz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_x_y[i] = 4.0 * g_yz_xz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_x_z[i] = 4.0 * g_yz_xz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (678-681)

    #pragma omp simd aligned(g_yz_xz_y_x, g_yz_xz_y_y, g_yz_xz_y_z, g_z_z_0_0_y_x_y_x, g_z_z_0_0_y_x_y_y, g_z_z_0_0_y_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_x_y_x[i] = 4.0 * g_yz_xz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_y_y[i] = 4.0 * g_yz_xz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_y_z[i] = 4.0 * g_yz_xz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (681-684)

    #pragma omp simd aligned(g_yz_xz_z_x, g_yz_xz_z_y, g_yz_xz_z_z, g_z_z_0_0_y_x_z_x, g_z_z_0_0_y_x_z_y, g_z_z_0_0_y_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_x_z_x[i] = 4.0 * g_yz_xz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_z_y[i] = 4.0 * g_yz_xz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_z_z[i] = 4.0 * g_yz_xz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (684-687)

    #pragma omp simd aligned(g_yz_yz_x_x, g_yz_yz_x_y, g_yz_yz_x_z, g_z_z_0_0_y_y_x_x, g_z_z_0_0_y_y_x_y, g_z_z_0_0_y_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_y_x_x[i] = 4.0 * g_yz_yz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_x_y[i] = 4.0 * g_yz_yz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_x_z[i] = 4.0 * g_yz_yz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (687-690)

    #pragma omp simd aligned(g_yz_yz_y_x, g_yz_yz_y_y, g_yz_yz_y_z, g_z_z_0_0_y_y_y_x, g_z_z_0_0_y_y_y_y, g_z_z_0_0_y_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_y_y_x[i] = 4.0 * g_yz_yz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_y_y[i] = 4.0 * g_yz_yz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_y_z[i] = 4.0 * g_yz_yz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (690-693)

    #pragma omp simd aligned(g_yz_yz_z_x, g_yz_yz_z_y, g_yz_yz_z_z, g_z_z_0_0_y_y_z_x, g_z_z_0_0_y_y_z_y, g_z_z_0_0_y_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_y_z_x[i] = 4.0 * g_yz_yz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_z_y[i] = 4.0 * g_yz_yz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_z_z[i] = 4.0 * g_yz_yz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (693-696)

    #pragma omp simd aligned(g_yz_0_x_x, g_yz_0_x_y, g_yz_0_x_z, g_yz_zz_x_x, g_yz_zz_x_y, g_yz_zz_x_z, g_z_z_0_0_y_z_x_x, g_z_z_0_0_y_z_x_y, g_z_z_0_0_y_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_z_x_x[i] = -2.0 * g_yz_0_x_x[i] * a_exp + 4.0 * g_yz_zz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_x_y[i] = -2.0 * g_yz_0_x_y[i] * a_exp + 4.0 * g_yz_zz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_x_z[i] = -2.0 * g_yz_0_x_z[i] * a_exp + 4.0 * g_yz_zz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (696-699)

    #pragma omp simd aligned(g_yz_0_y_x, g_yz_0_y_y, g_yz_0_y_z, g_yz_zz_y_x, g_yz_zz_y_y, g_yz_zz_y_z, g_z_z_0_0_y_z_y_x, g_z_z_0_0_y_z_y_y, g_z_z_0_0_y_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_z_y_x[i] = -2.0 * g_yz_0_y_x[i] * a_exp + 4.0 * g_yz_zz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_y_y[i] = -2.0 * g_yz_0_y_y[i] * a_exp + 4.0 * g_yz_zz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_y_z[i] = -2.0 * g_yz_0_y_z[i] * a_exp + 4.0 * g_yz_zz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (699-702)

    #pragma omp simd aligned(g_yz_0_z_x, g_yz_0_z_y, g_yz_0_z_z, g_yz_zz_z_x, g_yz_zz_z_y, g_yz_zz_z_z, g_z_z_0_0_y_z_z_x, g_z_z_0_0_y_z_z_y, g_z_z_0_0_y_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_z_z_x[i] = -2.0 * g_yz_0_z_x[i] * a_exp + 4.0 * g_yz_zz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_z_y[i] = -2.0 * g_yz_0_z_y[i] * a_exp + 4.0 * g_yz_zz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_z_z[i] = -2.0 * g_yz_0_z_z[i] * a_exp + 4.0 * g_yz_zz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (702-705)

    #pragma omp simd aligned(g_0_xz_x_x, g_0_xz_x_y, g_0_xz_x_z, g_z_z_0_0_z_x_x_x, g_z_z_0_0_z_x_x_y, g_z_z_0_0_z_x_x_z, g_zz_xz_x_x, g_zz_xz_x_y, g_zz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_x_x_x[i] = -2.0 * g_0_xz_x_x[i] * b_exp + 4.0 * g_zz_xz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_x_y[i] = -2.0 * g_0_xz_x_y[i] * b_exp + 4.0 * g_zz_xz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_x_z[i] = -2.0 * g_0_xz_x_z[i] * b_exp + 4.0 * g_zz_xz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (705-708)

    #pragma omp simd aligned(g_0_xz_y_x, g_0_xz_y_y, g_0_xz_y_z, g_z_z_0_0_z_x_y_x, g_z_z_0_0_z_x_y_y, g_z_z_0_0_z_x_y_z, g_zz_xz_y_x, g_zz_xz_y_y, g_zz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_x_y_x[i] = -2.0 * g_0_xz_y_x[i] * b_exp + 4.0 * g_zz_xz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_y_y[i] = -2.0 * g_0_xz_y_y[i] * b_exp + 4.0 * g_zz_xz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_y_z[i] = -2.0 * g_0_xz_y_z[i] * b_exp + 4.0 * g_zz_xz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (708-711)

    #pragma omp simd aligned(g_0_xz_z_x, g_0_xz_z_y, g_0_xz_z_z, g_z_z_0_0_z_x_z_x, g_z_z_0_0_z_x_z_y, g_z_z_0_0_z_x_z_z, g_zz_xz_z_x, g_zz_xz_z_y, g_zz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_x_z_x[i] = -2.0 * g_0_xz_z_x[i] * b_exp + 4.0 * g_zz_xz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_z_y[i] = -2.0 * g_0_xz_z_y[i] * b_exp + 4.0 * g_zz_xz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_z_z[i] = -2.0 * g_0_xz_z_z[i] * b_exp + 4.0 * g_zz_xz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (711-714)

    #pragma omp simd aligned(g_0_yz_x_x, g_0_yz_x_y, g_0_yz_x_z, g_z_z_0_0_z_y_x_x, g_z_z_0_0_z_y_x_y, g_z_z_0_0_z_y_x_z, g_zz_yz_x_x, g_zz_yz_x_y, g_zz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_y_x_x[i] = -2.0 * g_0_yz_x_x[i] * b_exp + 4.0 * g_zz_yz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_x_y[i] = -2.0 * g_0_yz_x_y[i] * b_exp + 4.0 * g_zz_yz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_x_z[i] = -2.0 * g_0_yz_x_z[i] * b_exp + 4.0 * g_zz_yz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (714-717)

    #pragma omp simd aligned(g_0_yz_y_x, g_0_yz_y_y, g_0_yz_y_z, g_z_z_0_0_z_y_y_x, g_z_z_0_0_z_y_y_y, g_z_z_0_0_z_y_y_z, g_zz_yz_y_x, g_zz_yz_y_y, g_zz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_y_y_x[i] = -2.0 * g_0_yz_y_x[i] * b_exp + 4.0 * g_zz_yz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_y_y[i] = -2.0 * g_0_yz_y_y[i] * b_exp + 4.0 * g_zz_yz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_y_z[i] = -2.0 * g_0_yz_y_z[i] * b_exp + 4.0 * g_zz_yz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (717-720)

    #pragma omp simd aligned(g_0_yz_z_x, g_0_yz_z_y, g_0_yz_z_z, g_z_z_0_0_z_y_z_x, g_z_z_0_0_z_y_z_y, g_z_z_0_0_z_y_z_z, g_zz_yz_z_x, g_zz_yz_z_y, g_zz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_y_z_x[i] = -2.0 * g_0_yz_z_x[i] * b_exp + 4.0 * g_zz_yz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_z_y[i] = -2.0 * g_0_yz_z_y[i] * b_exp + 4.0 * g_zz_yz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_z_z[i] = -2.0 * g_0_yz_z_z[i] * b_exp + 4.0 * g_zz_yz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (720-723)

    #pragma omp simd aligned(g_0_0_x_x, g_0_0_x_y, g_0_0_x_z, g_0_zz_x_x, g_0_zz_x_y, g_0_zz_x_z, g_z_z_0_0_z_z_x_x, g_z_z_0_0_z_z_x_y, g_z_z_0_0_z_z_x_z, g_zz_0_x_x, g_zz_0_x_y, g_zz_0_x_z, g_zz_zz_x_x, g_zz_zz_x_y, g_zz_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_z_x_x[i] = g_0_0_x_x[i] - 2.0 * g_0_zz_x_x[i] * b_exp - 2.0 * g_zz_0_x_x[i] * a_exp + 4.0 * g_zz_zz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_x_y[i] = g_0_0_x_y[i] - 2.0 * g_0_zz_x_y[i] * b_exp - 2.0 * g_zz_0_x_y[i] * a_exp + 4.0 * g_zz_zz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_x_z[i] = g_0_0_x_z[i] - 2.0 * g_0_zz_x_z[i] * b_exp - 2.0 * g_zz_0_x_z[i] * a_exp + 4.0 * g_zz_zz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (723-726)

    #pragma omp simd aligned(g_0_0_y_x, g_0_0_y_y, g_0_0_y_z, g_0_zz_y_x, g_0_zz_y_y, g_0_zz_y_z, g_z_z_0_0_z_z_y_x, g_z_z_0_0_z_z_y_y, g_z_z_0_0_z_z_y_z, g_zz_0_y_x, g_zz_0_y_y, g_zz_0_y_z, g_zz_zz_y_x, g_zz_zz_y_y, g_zz_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_z_y_x[i] = g_0_0_y_x[i] - 2.0 * g_0_zz_y_x[i] * b_exp - 2.0 * g_zz_0_y_x[i] * a_exp + 4.0 * g_zz_zz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_y_y[i] = g_0_0_y_y[i] - 2.0 * g_0_zz_y_y[i] * b_exp - 2.0 * g_zz_0_y_y[i] * a_exp + 4.0 * g_zz_zz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_y_z[i] = g_0_0_y_z[i] - 2.0 * g_0_zz_y_z[i] * b_exp - 2.0 * g_zz_0_y_z[i] * a_exp + 4.0 * g_zz_zz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (726-729)

    #pragma omp simd aligned(g_0_0_z_x, g_0_0_z_y, g_0_0_z_z, g_0_zz_z_x, g_0_zz_z_y, g_0_zz_z_z, g_z_z_0_0_z_z_z_x, g_z_z_0_0_z_z_z_y, g_z_z_0_0_z_z_z_z, g_zz_0_z_x, g_zz_0_z_y, g_zz_0_z_z, g_zz_zz_z_x, g_zz_zz_z_y, g_zz_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_z_z_x[i] = g_0_0_z_x[i] - 2.0 * g_0_zz_z_x[i] * b_exp - 2.0 * g_zz_0_z_x[i] * a_exp + 4.0 * g_zz_zz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_z_y[i] = g_0_0_z_y[i] - 2.0 * g_0_zz_z_y[i] * b_exp - 2.0 * g_zz_0_z_y[i] * a_exp + 4.0 * g_zz_zz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_z_z[i] = g_0_0_z_z[i] - 2.0 * g_0_zz_z_z[i] * b_exp - 2.0 * g_zz_0_z_z[i] * a_exp + 4.0 * g_zz_zz_z_z[i] * a_exp * b_exp;
    }
}

} // t4c_geom namespace


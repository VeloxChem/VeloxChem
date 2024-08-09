#include "GeomDeriv1010OfScalarForPPPP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1010_pppp_0(CSimdArray<double>& buffer_1010_pppp,
                     const CSimdArray<double>& buffer_spsp,
                     const CSimdArray<double>& buffer_spdp,
                     const CSimdArray<double>& buffer_dpsp,
                     const CSimdArray<double>& buffer_dpdp,
                     const double a_exp,
                     const double* c_exps) -> void
{
    const auto ndims = buffer_1010_pppp.number_of_columns();

    /// Set up components of auxilary buffer : buffer_spsp

    auto g_0_x_0_x = buffer_spsp[0];

    auto g_0_x_0_y = buffer_spsp[1];

    auto g_0_x_0_z = buffer_spsp[2];

    auto g_0_y_0_x = buffer_spsp[3];

    auto g_0_y_0_y = buffer_spsp[4];

    auto g_0_y_0_z = buffer_spsp[5];

    auto g_0_z_0_x = buffer_spsp[6];

    auto g_0_z_0_y = buffer_spsp[7];

    auto g_0_z_0_z = buffer_spsp[8];

    /// Set up components of auxilary buffer : buffer_spdp

    auto g_0_x_xx_x = buffer_spdp[0];

    auto g_0_x_xx_y = buffer_spdp[1];

    auto g_0_x_xx_z = buffer_spdp[2];

    auto g_0_x_xy_x = buffer_spdp[3];

    auto g_0_x_xy_y = buffer_spdp[4];

    auto g_0_x_xy_z = buffer_spdp[5];

    auto g_0_x_xz_x = buffer_spdp[6];

    auto g_0_x_xz_y = buffer_spdp[7];

    auto g_0_x_xz_z = buffer_spdp[8];

    auto g_0_x_yy_x = buffer_spdp[9];

    auto g_0_x_yy_y = buffer_spdp[10];

    auto g_0_x_yy_z = buffer_spdp[11];

    auto g_0_x_yz_x = buffer_spdp[12];

    auto g_0_x_yz_y = buffer_spdp[13];

    auto g_0_x_yz_z = buffer_spdp[14];

    auto g_0_x_zz_x = buffer_spdp[15];

    auto g_0_x_zz_y = buffer_spdp[16];

    auto g_0_x_zz_z = buffer_spdp[17];

    auto g_0_y_xx_x = buffer_spdp[18];

    auto g_0_y_xx_y = buffer_spdp[19];

    auto g_0_y_xx_z = buffer_spdp[20];

    auto g_0_y_xy_x = buffer_spdp[21];

    auto g_0_y_xy_y = buffer_spdp[22];

    auto g_0_y_xy_z = buffer_spdp[23];

    auto g_0_y_xz_x = buffer_spdp[24];

    auto g_0_y_xz_y = buffer_spdp[25];

    auto g_0_y_xz_z = buffer_spdp[26];

    auto g_0_y_yy_x = buffer_spdp[27];

    auto g_0_y_yy_y = buffer_spdp[28];

    auto g_0_y_yy_z = buffer_spdp[29];

    auto g_0_y_yz_x = buffer_spdp[30];

    auto g_0_y_yz_y = buffer_spdp[31];

    auto g_0_y_yz_z = buffer_spdp[32];

    auto g_0_y_zz_x = buffer_spdp[33];

    auto g_0_y_zz_y = buffer_spdp[34];

    auto g_0_y_zz_z = buffer_spdp[35];

    auto g_0_z_xx_x = buffer_spdp[36];

    auto g_0_z_xx_y = buffer_spdp[37];

    auto g_0_z_xx_z = buffer_spdp[38];

    auto g_0_z_xy_x = buffer_spdp[39];

    auto g_0_z_xy_y = buffer_spdp[40];

    auto g_0_z_xy_z = buffer_spdp[41];

    auto g_0_z_xz_x = buffer_spdp[42];

    auto g_0_z_xz_y = buffer_spdp[43];

    auto g_0_z_xz_z = buffer_spdp[44];

    auto g_0_z_yy_x = buffer_spdp[45];

    auto g_0_z_yy_y = buffer_spdp[46];

    auto g_0_z_yy_z = buffer_spdp[47];

    auto g_0_z_yz_x = buffer_spdp[48];

    auto g_0_z_yz_y = buffer_spdp[49];

    auto g_0_z_yz_z = buffer_spdp[50];

    auto g_0_z_zz_x = buffer_spdp[51];

    auto g_0_z_zz_y = buffer_spdp[52];

    auto g_0_z_zz_z = buffer_spdp[53];

    /// Set up components of auxilary buffer : buffer_dpsp

    auto g_xx_x_0_x = buffer_dpsp[0];

    auto g_xx_x_0_y = buffer_dpsp[1];

    auto g_xx_x_0_z = buffer_dpsp[2];

    auto g_xx_y_0_x = buffer_dpsp[3];

    auto g_xx_y_0_y = buffer_dpsp[4];

    auto g_xx_y_0_z = buffer_dpsp[5];

    auto g_xx_z_0_x = buffer_dpsp[6];

    auto g_xx_z_0_y = buffer_dpsp[7];

    auto g_xx_z_0_z = buffer_dpsp[8];

    auto g_xy_x_0_x = buffer_dpsp[9];

    auto g_xy_x_0_y = buffer_dpsp[10];

    auto g_xy_x_0_z = buffer_dpsp[11];

    auto g_xy_y_0_x = buffer_dpsp[12];

    auto g_xy_y_0_y = buffer_dpsp[13];

    auto g_xy_y_0_z = buffer_dpsp[14];

    auto g_xy_z_0_x = buffer_dpsp[15];

    auto g_xy_z_0_y = buffer_dpsp[16];

    auto g_xy_z_0_z = buffer_dpsp[17];

    auto g_xz_x_0_x = buffer_dpsp[18];

    auto g_xz_x_0_y = buffer_dpsp[19];

    auto g_xz_x_0_z = buffer_dpsp[20];

    auto g_xz_y_0_x = buffer_dpsp[21];

    auto g_xz_y_0_y = buffer_dpsp[22];

    auto g_xz_y_0_z = buffer_dpsp[23];

    auto g_xz_z_0_x = buffer_dpsp[24];

    auto g_xz_z_0_y = buffer_dpsp[25];

    auto g_xz_z_0_z = buffer_dpsp[26];

    auto g_yy_x_0_x = buffer_dpsp[27];

    auto g_yy_x_0_y = buffer_dpsp[28];

    auto g_yy_x_0_z = buffer_dpsp[29];

    auto g_yy_y_0_x = buffer_dpsp[30];

    auto g_yy_y_0_y = buffer_dpsp[31];

    auto g_yy_y_0_z = buffer_dpsp[32];

    auto g_yy_z_0_x = buffer_dpsp[33];

    auto g_yy_z_0_y = buffer_dpsp[34];

    auto g_yy_z_0_z = buffer_dpsp[35];

    auto g_yz_x_0_x = buffer_dpsp[36];

    auto g_yz_x_0_y = buffer_dpsp[37];

    auto g_yz_x_0_z = buffer_dpsp[38];

    auto g_yz_y_0_x = buffer_dpsp[39];

    auto g_yz_y_0_y = buffer_dpsp[40];

    auto g_yz_y_0_z = buffer_dpsp[41];

    auto g_yz_z_0_x = buffer_dpsp[42];

    auto g_yz_z_0_y = buffer_dpsp[43];

    auto g_yz_z_0_z = buffer_dpsp[44];

    auto g_zz_x_0_x = buffer_dpsp[45];

    auto g_zz_x_0_y = buffer_dpsp[46];

    auto g_zz_x_0_z = buffer_dpsp[47];

    auto g_zz_y_0_x = buffer_dpsp[48];

    auto g_zz_y_0_y = buffer_dpsp[49];

    auto g_zz_y_0_z = buffer_dpsp[50];

    auto g_zz_z_0_x = buffer_dpsp[51];

    auto g_zz_z_0_y = buffer_dpsp[52];

    auto g_zz_z_0_z = buffer_dpsp[53];

    /// Set up components of auxilary buffer : buffer_dpdp

    auto g_xx_x_xx_x = buffer_dpdp[0];

    auto g_xx_x_xx_y = buffer_dpdp[1];

    auto g_xx_x_xx_z = buffer_dpdp[2];

    auto g_xx_x_xy_x = buffer_dpdp[3];

    auto g_xx_x_xy_y = buffer_dpdp[4];

    auto g_xx_x_xy_z = buffer_dpdp[5];

    auto g_xx_x_xz_x = buffer_dpdp[6];

    auto g_xx_x_xz_y = buffer_dpdp[7];

    auto g_xx_x_xz_z = buffer_dpdp[8];

    auto g_xx_x_yy_x = buffer_dpdp[9];

    auto g_xx_x_yy_y = buffer_dpdp[10];

    auto g_xx_x_yy_z = buffer_dpdp[11];

    auto g_xx_x_yz_x = buffer_dpdp[12];

    auto g_xx_x_yz_y = buffer_dpdp[13];

    auto g_xx_x_yz_z = buffer_dpdp[14];

    auto g_xx_x_zz_x = buffer_dpdp[15];

    auto g_xx_x_zz_y = buffer_dpdp[16];

    auto g_xx_x_zz_z = buffer_dpdp[17];

    auto g_xx_y_xx_x = buffer_dpdp[18];

    auto g_xx_y_xx_y = buffer_dpdp[19];

    auto g_xx_y_xx_z = buffer_dpdp[20];

    auto g_xx_y_xy_x = buffer_dpdp[21];

    auto g_xx_y_xy_y = buffer_dpdp[22];

    auto g_xx_y_xy_z = buffer_dpdp[23];

    auto g_xx_y_xz_x = buffer_dpdp[24];

    auto g_xx_y_xz_y = buffer_dpdp[25];

    auto g_xx_y_xz_z = buffer_dpdp[26];

    auto g_xx_y_yy_x = buffer_dpdp[27];

    auto g_xx_y_yy_y = buffer_dpdp[28];

    auto g_xx_y_yy_z = buffer_dpdp[29];

    auto g_xx_y_yz_x = buffer_dpdp[30];

    auto g_xx_y_yz_y = buffer_dpdp[31];

    auto g_xx_y_yz_z = buffer_dpdp[32];

    auto g_xx_y_zz_x = buffer_dpdp[33];

    auto g_xx_y_zz_y = buffer_dpdp[34];

    auto g_xx_y_zz_z = buffer_dpdp[35];

    auto g_xx_z_xx_x = buffer_dpdp[36];

    auto g_xx_z_xx_y = buffer_dpdp[37];

    auto g_xx_z_xx_z = buffer_dpdp[38];

    auto g_xx_z_xy_x = buffer_dpdp[39];

    auto g_xx_z_xy_y = buffer_dpdp[40];

    auto g_xx_z_xy_z = buffer_dpdp[41];

    auto g_xx_z_xz_x = buffer_dpdp[42];

    auto g_xx_z_xz_y = buffer_dpdp[43];

    auto g_xx_z_xz_z = buffer_dpdp[44];

    auto g_xx_z_yy_x = buffer_dpdp[45];

    auto g_xx_z_yy_y = buffer_dpdp[46];

    auto g_xx_z_yy_z = buffer_dpdp[47];

    auto g_xx_z_yz_x = buffer_dpdp[48];

    auto g_xx_z_yz_y = buffer_dpdp[49];

    auto g_xx_z_yz_z = buffer_dpdp[50];

    auto g_xx_z_zz_x = buffer_dpdp[51];

    auto g_xx_z_zz_y = buffer_dpdp[52];

    auto g_xx_z_zz_z = buffer_dpdp[53];

    auto g_xy_x_xx_x = buffer_dpdp[54];

    auto g_xy_x_xx_y = buffer_dpdp[55];

    auto g_xy_x_xx_z = buffer_dpdp[56];

    auto g_xy_x_xy_x = buffer_dpdp[57];

    auto g_xy_x_xy_y = buffer_dpdp[58];

    auto g_xy_x_xy_z = buffer_dpdp[59];

    auto g_xy_x_xz_x = buffer_dpdp[60];

    auto g_xy_x_xz_y = buffer_dpdp[61];

    auto g_xy_x_xz_z = buffer_dpdp[62];

    auto g_xy_x_yy_x = buffer_dpdp[63];

    auto g_xy_x_yy_y = buffer_dpdp[64];

    auto g_xy_x_yy_z = buffer_dpdp[65];

    auto g_xy_x_yz_x = buffer_dpdp[66];

    auto g_xy_x_yz_y = buffer_dpdp[67];

    auto g_xy_x_yz_z = buffer_dpdp[68];

    auto g_xy_x_zz_x = buffer_dpdp[69];

    auto g_xy_x_zz_y = buffer_dpdp[70];

    auto g_xy_x_zz_z = buffer_dpdp[71];

    auto g_xy_y_xx_x = buffer_dpdp[72];

    auto g_xy_y_xx_y = buffer_dpdp[73];

    auto g_xy_y_xx_z = buffer_dpdp[74];

    auto g_xy_y_xy_x = buffer_dpdp[75];

    auto g_xy_y_xy_y = buffer_dpdp[76];

    auto g_xy_y_xy_z = buffer_dpdp[77];

    auto g_xy_y_xz_x = buffer_dpdp[78];

    auto g_xy_y_xz_y = buffer_dpdp[79];

    auto g_xy_y_xz_z = buffer_dpdp[80];

    auto g_xy_y_yy_x = buffer_dpdp[81];

    auto g_xy_y_yy_y = buffer_dpdp[82];

    auto g_xy_y_yy_z = buffer_dpdp[83];

    auto g_xy_y_yz_x = buffer_dpdp[84];

    auto g_xy_y_yz_y = buffer_dpdp[85];

    auto g_xy_y_yz_z = buffer_dpdp[86];

    auto g_xy_y_zz_x = buffer_dpdp[87];

    auto g_xy_y_zz_y = buffer_dpdp[88];

    auto g_xy_y_zz_z = buffer_dpdp[89];

    auto g_xy_z_xx_x = buffer_dpdp[90];

    auto g_xy_z_xx_y = buffer_dpdp[91];

    auto g_xy_z_xx_z = buffer_dpdp[92];

    auto g_xy_z_xy_x = buffer_dpdp[93];

    auto g_xy_z_xy_y = buffer_dpdp[94];

    auto g_xy_z_xy_z = buffer_dpdp[95];

    auto g_xy_z_xz_x = buffer_dpdp[96];

    auto g_xy_z_xz_y = buffer_dpdp[97];

    auto g_xy_z_xz_z = buffer_dpdp[98];

    auto g_xy_z_yy_x = buffer_dpdp[99];

    auto g_xy_z_yy_y = buffer_dpdp[100];

    auto g_xy_z_yy_z = buffer_dpdp[101];

    auto g_xy_z_yz_x = buffer_dpdp[102];

    auto g_xy_z_yz_y = buffer_dpdp[103];

    auto g_xy_z_yz_z = buffer_dpdp[104];

    auto g_xy_z_zz_x = buffer_dpdp[105];

    auto g_xy_z_zz_y = buffer_dpdp[106];

    auto g_xy_z_zz_z = buffer_dpdp[107];

    auto g_xz_x_xx_x = buffer_dpdp[108];

    auto g_xz_x_xx_y = buffer_dpdp[109];

    auto g_xz_x_xx_z = buffer_dpdp[110];

    auto g_xz_x_xy_x = buffer_dpdp[111];

    auto g_xz_x_xy_y = buffer_dpdp[112];

    auto g_xz_x_xy_z = buffer_dpdp[113];

    auto g_xz_x_xz_x = buffer_dpdp[114];

    auto g_xz_x_xz_y = buffer_dpdp[115];

    auto g_xz_x_xz_z = buffer_dpdp[116];

    auto g_xz_x_yy_x = buffer_dpdp[117];

    auto g_xz_x_yy_y = buffer_dpdp[118];

    auto g_xz_x_yy_z = buffer_dpdp[119];

    auto g_xz_x_yz_x = buffer_dpdp[120];

    auto g_xz_x_yz_y = buffer_dpdp[121];

    auto g_xz_x_yz_z = buffer_dpdp[122];

    auto g_xz_x_zz_x = buffer_dpdp[123];

    auto g_xz_x_zz_y = buffer_dpdp[124];

    auto g_xz_x_zz_z = buffer_dpdp[125];

    auto g_xz_y_xx_x = buffer_dpdp[126];

    auto g_xz_y_xx_y = buffer_dpdp[127];

    auto g_xz_y_xx_z = buffer_dpdp[128];

    auto g_xz_y_xy_x = buffer_dpdp[129];

    auto g_xz_y_xy_y = buffer_dpdp[130];

    auto g_xz_y_xy_z = buffer_dpdp[131];

    auto g_xz_y_xz_x = buffer_dpdp[132];

    auto g_xz_y_xz_y = buffer_dpdp[133];

    auto g_xz_y_xz_z = buffer_dpdp[134];

    auto g_xz_y_yy_x = buffer_dpdp[135];

    auto g_xz_y_yy_y = buffer_dpdp[136];

    auto g_xz_y_yy_z = buffer_dpdp[137];

    auto g_xz_y_yz_x = buffer_dpdp[138];

    auto g_xz_y_yz_y = buffer_dpdp[139];

    auto g_xz_y_yz_z = buffer_dpdp[140];

    auto g_xz_y_zz_x = buffer_dpdp[141];

    auto g_xz_y_zz_y = buffer_dpdp[142];

    auto g_xz_y_zz_z = buffer_dpdp[143];

    auto g_xz_z_xx_x = buffer_dpdp[144];

    auto g_xz_z_xx_y = buffer_dpdp[145];

    auto g_xz_z_xx_z = buffer_dpdp[146];

    auto g_xz_z_xy_x = buffer_dpdp[147];

    auto g_xz_z_xy_y = buffer_dpdp[148];

    auto g_xz_z_xy_z = buffer_dpdp[149];

    auto g_xz_z_xz_x = buffer_dpdp[150];

    auto g_xz_z_xz_y = buffer_dpdp[151];

    auto g_xz_z_xz_z = buffer_dpdp[152];

    auto g_xz_z_yy_x = buffer_dpdp[153];

    auto g_xz_z_yy_y = buffer_dpdp[154];

    auto g_xz_z_yy_z = buffer_dpdp[155];

    auto g_xz_z_yz_x = buffer_dpdp[156];

    auto g_xz_z_yz_y = buffer_dpdp[157];

    auto g_xz_z_yz_z = buffer_dpdp[158];

    auto g_xz_z_zz_x = buffer_dpdp[159];

    auto g_xz_z_zz_y = buffer_dpdp[160];

    auto g_xz_z_zz_z = buffer_dpdp[161];

    auto g_yy_x_xx_x = buffer_dpdp[162];

    auto g_yy_x_xx_y = buffer_dpdp[163];

    auto g_yy_x_xx_z = buffer_dpdp[164];

    auto g_yy_x_xy_x = buffer_dpdp[165];

    auto g_yy_x_xy_y = buffer_dpdp[166];

    auto g_yy_x_xy_z = buffer_dpdp[167];

    auto g_yy_x_xz_x = buffer_dpdp[168];

    auto g_yy_x_xz_y = buffer_dpdp[169];

    auto g_yy_x_xz_z = buffer_dpdp[170];

    auto g_yy_x_yy_x = buffer_dpdp[171];

    auto g_yy_x_yy_y = buffer_dpdp[172];

    auto g_yy_x_yy_z = buffer_dpdp[173];

    auto g_yy_x_yz_x = buffer_dpdp[174];

    auto g_yy_x_yz_y = buffer_dpdp[175];

    auto g_yy_x_yz_z = buffer_dpdp[176];

    auto g_yy_x_zz_x = buffer_dpdp[177];

    auto g_yy_x_zz_y = buffer_dpdp[178];

    auto g_yy_x_zz_z = buffer_dpdp[179];

    auto g_yy_y_xx_x = buffer_dpdp[180];

    auto g_yy_y_xx_y = buffer_dpdp[181];

    auto g_yy_y_xx_z = buffer_dpdp[182];

    auto g_yy_y_xy_x = buffer_dpdp[183];

    auto g_yy_y_xy_y = buffer_dpdp[184];

    auto g_yy_y_xy_z = buffer_dpdp[185];

    auto g_yy_y_xz_x = buffer_dpdp[186];

    auto g_yy_y_xz_y = buffer_dpdp[187];

    auto g_yy_y_xz_z = buffer_dpdp[188];

    auto g_yy_y_yy_x = buffer_dpdp[189];

    auto g_yy_y_yy_y = buffer_dpdp[190];

    auto g_yy_y_yy_z = buffer_dpdp[191];

    auto g_yy_y_yz_x = buffer_dpdp[192];

    auto g_yy_y_yz_y = buffer_dpdp[193];

    auto g_yy_y_yz_z = buffer_dpdp[194];

    auto g_yy_y_zz_x = buffer_dpdp[195];

    auto g_yy_y_zz_y = buffer_dpdp[196];

    auto g_yy_y_zz_z = buffer_dpdp[197];

    auto g_yy_z_xx_x = buffer_dpdp[198];

    auto g_yy_z_xx_y = buffer_dpdp[199];

    auto g_yy_z_xx_z = buffer_dpdp[200];

    auto g_yy_z_xy_x = buffer_dpdp[201];

    auto g_yy_z_xy_y = buffer_dpdp[202];

    auto g_yy_z_xy_z = buffer_dpdp[203];

    auto g_yy_z_xz_x = buffer_dpdp[204];

    auto g_yy_z_xz_y = buffer_dpdp[205];

    auto g_yy_z_xz_z = buffer_dpdp[206];

    auto g_yy_z_yy_x = buffer_dpdp[207];

    auto g_yy_z_yy_y = buffer_dpdp[208];

    auto g_yy_z_yy_z = buffer_dpdp[209];

    auto g_yy_z_yz_x = buffer_dpdp[210];

    auto g_yy_z_yz_y = buffer_dpdp[211];

    auto g_yy_z_yz_z = buffer_dpdp[212];

    auto g_yy_z_zz_x = buffer_dpdp[213];

    auto g_yy_z_zz_y = buffer_dpdp[214];

    auto g_yy_z_zz_z = buffer_dpdp[215];

    auto g_yz_x_xx_x = buffer_dpdp[216];

    auto g_yz_x_xx_y = buffer_dpdp[217];

    auto g_yz_x_xx_z = buffer_dpdp[218];

    auto g_yz_x_xy_x = buffer_dpdp[219];

    auto g_yz_x_xy_y = buffer_dpdp[220];

    auto g_yz_x_xy_z = buffer_dpdp[221];

    auto g_yz_x_xz_x = buffer_dpdp[222];

    auto g_yz_x_xz_y = buffer_dpdp[223];

    auto g_yz_x_xz_z = buffer_dpdp[224];

    auto g_yz_x_yy_x = buffer_dpdp[225];

    auto g_yz_x_yy_y = buffer_dpdp[226];

    auto g_yz_x_yy_z = buffer_dpdp[227];

    auto g_yz_x_yz_x = buffer_dpdp[228];

    auto g_yz_x_yz_y = buffer_dpdp[229];

    auto g_yz_x_yz_z = buffer_dpdp[230];

    auto g_yz_x_zz_x = buffer_dpdp[231];

    auto g_yz_x_zz_y = buffer_dpdp[232];

    auto g_yz_x_zz_z = buffer_dpdp[233];

    auto g_yz_y_xx_x = buffer_dpdp[234];

    auto g_yz_y_xx_y = buffer_dpdp[235];

    auto g_yz_y_xx_z = buffer_dpdp[236];

    auto g_yz_y_xy_x = buffer_dpdp[237];

    auto g_yz_y_xy_y = buffer_dpdp[238];

    auto g_yz_y_xy_z = buffer_dpdp[239];

    auto g_yz_y_xz_x = buffer_dpdp[240];

    auto g_yz_y_xz_y = buffer_dpdp[241];

    auto g_yz_y_xz_z = buffer_dpdp[242];

    auto g_yz_y_yy_x = buffer_dpdp[243];

    auto g_yz_y_yy_y = buffer_dpdp[244];

    auto g_yz_y_yy_z = buffer_dpdp[245];

    auto g_yz_y_yz_x = buffer_dpdp[246];

    auto g_yz_y_yz_y = buffer_dpdp[247];

    auto g_yz_y_yz_z = buffer_dpdp[248];

    auto g_yz_y_zz_x = buffer_dpdp[249];

    auto g_yz_y_zz_y = buffer_dpdp[250];

    auto g_yz_y_zz_z = buffer_dpdp[251];

    auto g_yz_z_xx_x = buffer_dpdp[252];

    auto g_yz_z_xx_y = buffer_dpdp[253];

    auto g_yz_z_xx_z = buffer_dpdp[254];

    auto g_yz_z_xy_x = buffer_dpdp[255];

    auto g_yz_z_xy_y = buffer_dpdp[256];

    auto g_yz_z_xy_z = buffer_dpdp[257];

    auto g_yz_z_xz_x = buffer_dpdp[258];

    auto g_yz_z_xz_y = buffer_dpdp[259];

    auto g_yz_z_xz_z = buffer_dpdp[260];

    auto g_yz_z_yy_x = buffer_dpdp[261];

    auto g_yz_z_yy_y = buffer_dpdp[262];

    auto g_yz_z_yy_z = buffer_dpdp[263];

    auto g_yz_z_yz_x = buffer_dpdp[264];

    auto g_yz_z_yz_y = buffer_dpdp[265];

    auto g_yz_z_yz_z = buffer_dpdp[266];

    auto g_yz_z_zz_x = buffer_dpdp[267];

    auto g_yz_z_zz_y = buffer_dpdp[268];

    auto g_yz_z_zz_z = buffer_dpdp[269];

    auto g_zz_x_xx_x = buffer_dpdp[270];

    auto g_zz_x_xx_y = buffer_dpdp[271];

    auto g_zz_x_xx_z = buffer_dpdp[272];

    auto g_zz_x_xy_x = buffer_dpdp[273];

    auto g_zz_x_xy_y = buffer_dpdp[274];

    auto g_zz_x_xy_z = buffer_dpdp[275];

    auto g_zz_x_xz_x = buffer_dpdp[276];

    auto g_zz_x_xz_y = buffer_dpdp[277];

    auto g_zz_x_xz_z = buffer_dpdp[278];

    auto g_zz_x_yy_x = buffer_dpdp[279];

    auto g_zz_x_yy_y = buffer_dpdp[280];

    auto g_zz_x_yy_z = buffer_dpdp[281];

    auto g_zz_x_yz_x = buffer_dpdp[282];

    auto g_zz_x_yz_y = buffer_dpdp[283];

    auto g_zz_x_yz_z = buffer_dpdp[284];

    auto g_zz_x_zz_x = buffer_dpdp[285];

    auto g_zz_x_zz_y = buffer_dpdp[286];

    auto g_zz_x_zz_z = buffer_dpdp[287];

    auto g_zz_y_xx_x = buffer_dpdp[288];

    auto g_zz_y_xx_y = buffer_dpdp[289];

    auto g_zz_y_xx_z = buffer_dpdp[290];

    auto g_zz_y_xy_x = buffer_dpdp[291];

    auto g_zz_y_xy_y = buffer_dpdp[292];

    auto g_zz_y_xy_z = buffer_dpdp[293];

    auto g_zz_y_xz_x = buffer_dpdp[294];

    auto g_zz_y_xz_y = buffer_dpdp[295];

    auto g_zz_y_xz_z = buffer_dpdp[296];

    auto g_zz_y_yy_x = buffer_dpdp[297];

    auto g_zz_y_yy_y = buffer_dpdp[298];

    auto g_zz_y_yy_z = buffer_dpdp[299];

    auto g_zz_y_yz_x = buffer_dpdp[300];

    auto g_zz_y_yz_y = buffer_dpdp[301];

    auto g_zz_y_yz_z = buffer_dpdp[302];

    auto g_zz_y_zz_x = buffer_dpdp[303];

    auto g_zz_y_zz_y = buffer_dpdp[304];

    auto g_zz_y_zz_z = buffer_dpdp[305];

    auto g_zz_z_xx_x = buffer_dpdp[306];

    auto g_zz_z_xx_y = buffer_dpdp[307];

    auto g_zz_z_xx_z = buffer_dpdp[308];

    auto g_zz_z_xy_x = buffer_dpdp[309];

    auto g_zz_z_xy_y = buffer_dpdp[310];

    auto g_zz_z_xy_z = buffer_dpdp[311];

    auto g_zz_z_xz_x = buffer_dpdp[312];

    auto g_zz_z_xz_y = buffer_dpdp[313];

    auto g_zz_z_xz_z = buffer_dpdp[314];

    auto g_zz_z_yy_x = buffer_dpdp[315];

    auto g_zz_z_yy_y = buffer_dpdp[316];

    auto g_zz_z_yy_z = buffer_dpdp[317];

    auto g_zz_z_yz_x = buffer_dpdp[318];

    auto g_zz_z_yz_y = buffer_dpdp[319];

    auto g_zz_z_yz_z = buffer_dpdp[320];

    auto g_zz_z_zz_x = buffer_dpdp[321];

    auto g_zz_z_zz_y = buffer_dpdp[322];

    auto g_zz_z_zz_z = buffer_dpdp[323];

    /// Set up components of integrals buffer : buffer_1010_pppp

    auto g_x_0_x_0_x_x_x_x = buffer_1010_pppp[0];

    auto g_x_0_x_0_x_x_x_y = buffer_1010_pppp[1];

    auto g_x_0_x_0_x_x_x_z = buffer_1010_pppp[2];

    auto g_x_0_x_0_x_x_y_x = buffer_1010_pppp[3];

    auto g_x_0_x_0_x_x_y_y = buffer_1010_pppp[4];

    auto g_x_0_x_0_x_x_y_z = buffer_1010_pppp[5];

    auto g_x_0_x_0_x_x_z_x = buffer_1010_pppp[6];

    auto g_x_0_x_0_x_x_z_y = buffer_1010_pppp[7];

    auto g_x_0_x_0_x_x_z_z = buffer_1010_pppp[8];

    auto g_x_0_x_0_x_y_x_x = buffer_1010_pppp[9];

    auto g_x_0_x_0_x_y_x_y = buffer_1010_pppp[10];

    auto g_x_0_x_0_x_y_x_z = buffer_1010_pppp[11];

    auto g_x_0_x_0_x_y_y_x = buffer_1010_pppp[12];

    auto g_x_0_x_0_x_y_y_y = buffer_1010_pppp[13];

    auto g_x_0_x_0_x_y_y_z = buffer_1010_pppp[14];

    auto g_x_0_x_0_x_y_z_x = buffer_1010_pppp[15];

    auto g_x_0_x_0_x_y_z_y = buffer_1010_pppp[16];

    auto g_x_0_x_0_x_y_z_z = buffer_1010_pppp[17];

    auto g_x_0_x_0_x_z_x_x = buffer_1010_pppp[18];

    auto g_x_0_x_0_x_z_x_y = buffer_1010_pppp[19];

    auto g_x_0_x_0_x_z_x_z = buffer_1010_pppp[20];

    auto g_x_0_x_0_x_z_y_x = buffer_1010_pppp[21];

    auto g_x_0_x_0_x_z_y_y = buffer_1010_pppp[22];

    auto g_x_0_x_0_x_z_y_z = buffer_1010_pppp[23];

    auto g_x_0_x_0_x_z_z_x = buffer_1010_pppp[24];

    auto g_x_0_x_0_x_z_z_y = buffer_1010_pppp[25];

    auto g_x_0_x_0_x_z_z_z = buffer_1010_pppp[26];

    auto g_x_0_x_0_y_x_x_x = buffer_1010_pppp[27];

    auto g_x_0_x_0_y_x_x_y = buffer_1010_pppp[28];

    auto g_x_0_x_0_y_x_x_z = buffer_1010_pppp[29];

    auto g_x_0_x_0_y_x_y_x = buffer_1010_pppp[30];

    auto g_x_0_x_0_y_x_y_y = buffer_1010_pppp[31];

    auto g_x_0_x_0_y_x_y_z = buffer_1010_pppp[32];

    auto g_x_0_x_0_y_x_z_x = buffer_1010_pppp[33];

    auto g_x_0_x_0_y_x_z_y = buffer_1010_pppp[34];

    auto g_x_0_x_0_y_x_z_z = buffer_1010_pppp[35];

    auto g_x_0_x_0_y_y_x_x = buffer_1010_pppp[36];

    auto g_x_0_x_0_y_y_x_y = buffer_1010_pppp[37];

    auto g_x_0_x_0_y_y_x_z = buffer_1010_pppp[38];

    auto g_x_0_x_0_y_y_y_x = buffer_1010_pppp[39];

    auto g_x_0_x_0_y_y_y_y = buffer_1010_pppp[40];

    auto g_x_0_x_0_y_y_y_z = buffer_1010_pppp[41];

    auto g_x_0_x_0_y_y_z_x = buffer_1010_pppp[42];

    auto g_x_0_x_0_y_y_z_y = buffer_1010_pppp[43];

    auto g_x_0_x_0_y_y_z_z = buffer_1010_pppp[44];

    auto g_x_0_x_0_y_z_x_x = buffer_1010_pppp[45];

    auto g_x_0_x_0_y_z_x_y = buffer_1010_pppp[46];

    auto g_x_0_x_0_y_z_x_z = buffer_1010_pppp[47];

    auto g_x_0_x_0_y_z_y_x = buffer_1010_pppp[48];

    auto g_x_0_x_0_y_z_y_y = buffer_1010_pppp[49];

    auto g_x_0_x_0_y_z_y_z = buffer_1010_pppp[50];

    auto g_x_0_x_0_y_z_z_x = buffer_1010_pppp[51];

    auto g_x_0_x_0_y_z_z_y = buffer_1010_pppp[52];

    auto g_x_0_x_0_y_z_z_z = buffer_1010_pppp[53];

    auto g_x_0_x_0_z_x_x_x = buffer_1010_pppp[54];

    auto g_x_0_x_0_z_x_x_y = buffer_1010_pppp[55];

    auto g_x_0_x_0_z_x_x_z = buffer_1010_pppp[56];

    auto g_x_0_x_0_z_x_y_x = buffer_1010_pppp[57];

    auto g_x_0_x_0_z_x_y_y = buffer_1010_pppp[58];

    auto g_x_0_x_0_z_x_y_z = buffer_1010_pppp[59];

    auto g_x_0_x_0_z_x_z_x = buffer_1010_pppp[60];

    auto g_x_0_x_0_z_x_z_y = buffer_1010_pppp[61];

    auto g_x_0_x_0_z_x_z_z = buffer_1010_pppp[62];

    auto g_x_0_x_0_z_y_x_x = buffer_1010_pppp[63];

    auto g_x_0_x_0_z_y_x_y = buffer_1010_pppp[64];

    auto g_x_0_x_0_z_y_x_z = buffer_1010_pppp[65];

    auto g_x_0_x_0_z_y_y_x = buffer_1010_pppp[66];

    auto g_x_0_x_0_z_y_y_y = buffer_1010_pppp[67];

    auto g_x_0_x_0_z_y_y_z = buffer_1010_pppp[68];

    auto g_x_0_x_0_z_y_z_x = buffer_1010_pppp[69];

    auto g_x_0_x_0_z_y_z_y = buffer_1010_pppp[70];

    auto g_x_0_x_0_z_y_z_z = buffer_1010_pppp[71];

    auto g_x_0_x_0_z_z_x_x = buffer_1010_pppp[72];

    auto g_x_0_x_0_z_z_x_y = buffer_1010_pppp[73];

    auto g_x_0_x_0_z_z_x_z = buffer_1010_pppp[74];

    auto g_x_0_x_0_z_z_y_x = buffer_1010_pppp[75];

    auto g_x_0_x_0_z_z_y_y = buffer_1010_pppp[76];

    auto g_x_0_x_0_z_z_y_z = buffer_1010_pppp[77];

    auto g_x_0_x_0_z_z_z_x = buffer_1010_pppp[78];

    auto g_x_0_x_0_z_z_z_y = buffer_1010_pppp[79];

    auto g_x_0_x_0_z_z_z_z = buffer_1010_pppp[80];

    auto g_x_0_y_0_x_x_x_x = buffer_1010_pppp[81];

    auto g_x_0_y_0_x_x_x_y = buffer_1010_pppp[82];

    auto g_x_0_y_0_x_x_x_z = buffer_1010_pppp[83];

    auto g_x_0_y_0_x_x_y_x = buffer_1010_pppp[84];

    auto g_x_0_y_0_x_x_y_y = buffer_1010_pppp[85];

    auto g_x_0_y_0_x_x_y_z = buffer_1010_pppp[86];

    auto g_x_0_y_0_x_x_z_x = buffer_1010_pppp[87];

    auto g_x_0_y_0_x_x_z_y = buffer_1010_pppp[88];

    auto g_x_0_y_0_x_x_z_z = buffer_1010_pppp[89];

    auto g_x_0_y_0_x_y_x_x = buffer_1010_pppp[90];

    auto g_x_0_y_0_x_y_x_y = buffer_1010_pppp[91];

    auto g_x_0_y_0_x_y_x_z = buffer_1010_pppp[92];

    auto g_x_0_y_0_x_y_y_x = buffer_1010_pppp[93];

    auto g_x_0_y_0_x_y_y_y = buffer_1010_pppp[94];

    auto g_x_0_y_0_x_y_y_z = buffer_1010_pppp[95];

    auto g_x_0_y_0_x_y_z_x = buffer_1010_pppp[96];

    auto g_x_0_y_0_x_y_z_y = buffer_1010_pppp[97];

    auto g_x_0_y_0_x_y_z_z = buffer_1010_pppp[98];

    auto g_x_0_y_0_x_z_x_x = buffer_1010_pppp[99];

    auto g_x_0_y_0_x_z_x_y = buffer_1010_pppp[100];

    auto g_x_0_y_0_x_z_x_z = buffer_1010_pppp[101];

    auto g_x_0_y_0_x_z_y_x = buffer_1010_pppp[102];

    auto g_x_0_y_0_x_z_y_y = buffer_1010_pppp[103];

    auto g_x_0_y_0_x_z_y_z = buffer_1010_pppp[104];

    auto g_x_0_y_0_x_z_z_x = buffer_1010_pppp[105];

    auto g_x_0_y_0_x_z_z_y = buffer_1010_pppp[106];

    auto g_x_0_y_0_x_z_z_z = buffer_1010_pppp[107];

    auto g_x_0_y_0_y_x_x_x = buffer_1010_pppp[108];

    auto g_x_0_y_0_y_x_x_y = buffer_1010_pppp[109];

    auto g_x_0_y_0_y_x_x_z = buffer_1010_pppp[110];

    auto g_x_0_y_0_y_x_y_x = buffer_1010_pppp[111];

    auto g_x_0_y_0_y_x_y_y = buffer_1010_pppp[112];

    auto g_x_0_y_0_y_x_y_z = buffer_1010_pppp[113];

    auto g_x_0_y_0_y_x_z_x = buffer_1010_pppp[114];

    auto g_x_0_y_0_y_x_z_y = buffer_1010_pppp[115];

    auto g_x_0_y_0_y_x_z_z = buffer_1010_pppp[116];

    auto g_x_0_y_0_y_y_x_x = buffer_1010_pppp[117];

    auto g_x_0_y_0_y_y_x_y = buffer_1010_pppp[118];

    auto g_x_0_y_0_y_y_x_z = buffer_1010_pppp[119];

    auto g_x_0_y_0_y_y_y_x = buffer_1010_pppp[120];

    auto g_x_0_y_0_y_y_y_y = buffer_1010_pppp[121];

    auto g_x_0_y_0_y_y_y_z = buffer_1010_pppp[122];

    auto g_x_0_y_0_y_y_z_x = buffer_1010_pppp[123];

    auto g_x_0_y_0_y_y_z_y = buffer_1010_pppp[124];

    auto g_x_0_y_0_y_y_z_z = buffer_1010_pppp[125];

    auto g_x_0_y_0_y_z_x_x = buffer_1010_pppp[126];

    auto g_x_0_y_0_y_z_x_y = buffer_1010_pppp[127];

    auto g_x_0_y_0_y_z_x_z = buffer_1010_pppp[128];

    auto g_x_0_y_0_y_z_y_x = buffer_1010_pppp[129];

    auto g_x_0_y_0_y_z_y_y = buffer_1010_pppp[130];

    auto g_x_0_y_0_y_z_y_z = buffer_1010_pppp[131];

    auto g_x_0_y_0_y_z_z_x = buffer_1010_pppp[132];

    auto g_x_0_y_0_y_z_z_y = buffer_1010_pppp[133];

    auto g_x_0_y_0_y_z_z_z = buffer_1010_pppp[134];

    auto g_x_0_y_0_z_x_x_x = buffer_1010_pppp[135];

    auto g_x_0_y_0_z_x_x_y = buffer_1010_pppp[136];

    auto g_x_0_y_0_z_x_x_z = buffer_1010_pppp[137];

    auto g_x_0_y_0_z_x_y_x = buffer_1010_pppp[138];

    auto g_x_0_y_0_z_x_y_y = buffer_1010_pppp[139];

    auto g_x_0_y_0_z_x_y_z = buffer_1010_pppp[140];

    auto g_x_0_y_0_z_x_z_x = buffer_1010_pppp[141];

    auto g_x_0_y_0_z_x_z_y = buffer_1010_pppp[142];

    auto g_x_0_y_0_z_x_z_z = buffer_1010_pppp[143];

    auto g_x_0_y_0_z_y_x_x = buffer_1010_pppp[144];

    auto g_x_0_y_0_z_y_x_y = buffer_1010_pppp[145];

    auto g_x_0_y_0_z_y_x_z = buffer_1010_pppp[146];

    auto g_x_0_y_0_z_y_y_x = buffer_1010_pppp[147];

    auto g_x_0_y_0_z_y_y_y = buffer_1010_pppp[148];

    auto g_x_0_y_0_z_y_y_z = buffer_1010_pppp[149];

    auto g_x_0_y_0_z_y_z_x = buffer_1010_pppp[150];

    auto g_x_0_y_0_z_y_z_y = buffer_1010_pppp[151];

    auto g_x_0_y_0_z_y_z_z = buffer_1010_pppp[152];

    auto g_x_0_y_0_z_z_x_x = buffer_1010_pppp[153];

    auto g_x_0_y_0_z_z_x_y = buffer_1010_pppp[154];

    auto g_x_0_y_0_z_z_x_z = buffer_1010_pppp[155];

    auto g_x_0_y_0_z_z_y_x = buffer_1010_pppp[156];

    auto g_x_0_y_0_z_z_y_y = buffer_1010_pppp[157];

    auto g_x_0_y_0_z_z_y_z = buffer_1010_pppp[158];

    auto g_x_0_y_0_z_z_z_x = buffer_1010_pppp[159];

    auto g_x_0_y_0_z_z_z_y = buffer_1010_pppp[160];

    auto g_x_0_y_0_z_z_z_z = buffer_1010_pppp[161];

    auto g_x_0_z_0_x_x_x_x = buffer_1010_pppp[162];

    auto g_x_0_z_0_x_x_x_y = buffer_1010_pppp[163];

    auto g_x_0_z_0_x_x_x_z = buffer_1010_pppp[164];

    auto g_x_0_z_0_x_x_y_x = buffer_1010_pppp[165];

    auto g_x_0_z_0_x_x_y_y = buffer_1010_pppp[166];

    auto g_x_0_z_0_x_x_y_z = buffer_1010_pppp[167];

    auto g_x_0_z_0_x_x_z_x = buffer_1010_pppp[168];

    auto g_x_0_z_0_x_x_z_y = buffer_1010_pppp[169];

    auto g_x_0_z_0_x_x_z_z = buffer_1010_pppp[170];

    auto g_x_0_z_0_x_y_x_x = buffer_1010_pppp[171];

    auto g_x_0_z_0_x_y_x_y = buffer_1010_pppp[172];

    auto g_x_0_z_0_x_y_x_z = buffer_1010_pppp[173];

    auto g_x_0_z_0_x_y_y_x = buffer_1010_pppp[174];

    auto g_x_0_z_0_x_y_y_y = buffer_1010_pppp[175];

    auto g_x_0_z_0_x_y_y_z = buffer_1010_pppp[176];

    auto g_x_0_z_0_x_y_z_x = buffer_1010_pppp[177];

    auto g_x_0_z_0_x_y_z_y = buffer_1010_pppp[178];

    auto g_x_0_z_0_x_y_z_z = buffer_1010_pppp[179];

    auto g_x_0_z_0_x_z_x_x = buffer_1010_pppp[180];

    auto g_x_0_z_0_x_z_x_y = buffer_1010_pppp[181];

    auto g_x_0_z_0_x_z_x_z = buffer_1010_pppp[182];

    auto g_x_0_z_0_x_z_y_x = buffer_1010_pppp[183];

    auto g_x_0_z_0_x_z_y_y = buffer_1010_pppp[184];

    auto g_x_0_z_0_x_z_y_z = buffer_1010_pppp[185];

    auto g_x_0_z_0_x_z_z_x = buffer_1010_pppp[186];

    auto g_x_0_z_0_x_z_z_y = buffer_1010_pppp[187];

    auto g_x_0_z_0_x_z_z_z = buffer_1010_pppp[188];

    auto g_x_0_z_0_y_x_x_x = buffer_1010_pppp[189];

    auto g_x_0_z_0_y_x_x_y = buffer_1010_pppp[190];

    auto g_x_0_z_0_y_x_x_z = buffer_1010_pppp[191];

    auto g_x_0_z_0_y_x_y_x = buffer_1010_pppp[192];

    auto g_x_0_z_0_y_x_y_y = buffer_1010_pppp[193];

    auto g_x_0_z_0_y_x_y_z = buffer_1010_pppp[194];

    auto g_x_0_z_0_y_x_z_x = buffer_1010_pppp[195];

    auto g_x_0_z_0_y_x_z_y = buffer_1010_pppp[196];

    auto g_x_0_z_0_y_x_z_z = buffer_1010_pppp[197];

    auto g_x_0_z_0_y_y_x_x = buffer_1010_pppp[198];

    auto g_x_0_z_0_y_y_x_y = buffer_1010_pppp[199];

    auto g_x_0_z_0_y_y_x_z = buffer_1010_pppp[200];

    auto g_x_0_z_0_y_y_y_x = buffer_1010_pppp[201];

    auto g_x_0_z_0_y_y_y_y = buffer_1010_pppp[202];

    auto g_x_0_z_0_y_y_y_z = buffer_1010_pppp[203];

    auto g_x_0_z_0_y_y_z_x = buffer_1010_pppp[204];

    auto g_x_0_z_0_y_y_z_y = buffer_1010_pppp[205];

    auto g_x_0_z_0_y_y_z_z = buffer_1010_pppp[206];

    auto g_x_0_z_0_y_z_x_x = buffer_1010_pppp[207];

    auto g_x_0_z_0_y_z_x_y = buffer_1010_pppp[208];

    auto g_x_0_z_0_y_z_x_z = buffer_1010_pppp[209];

    auto g_x_0_z_0_y_z_y_x = buffer_1010_pppp[210];

    auto g_x_0_z_0_y_z_y_y = buffer_1010_pppp[211];

    auto g_x_0_z_0_y_z_y_z = buffer_1010_pppp[212];

    auto g_x_0_z_0_y_z_z_x = buffer_1010_pppp[213];

    auto g_x_0_z_0_y_z_z_y = buffer_1010_pppp[214];

    auto g_x_0_z_0_y_z_z_z = buffer_1010_pppp[215];

    auto g_x_0_z_0_z_x_x_x = buffer_1010_pppp[216];

    auto g_x_0_z_0_z_x_x_y = buffer_1010_pppp[217];

    auto g_x_0_z_0_z_x_x_z = buffer_1010_pppp[218];

    auto g_x_0_z_0_z_x_y_x = buffer_1010_pppp[219];

    auto g_x_0_z_0_z_x_y_y = buffer_1010_pppp[220];

    auto g_x_0_z_0_z_x_y_z = buffer_1010_pppp[221];

    auto g_x_0_z_0_z_x_z_x = buffer_1010_pppp[222];

    auto g_x_0_z_0_z_x_z_y = buffer_1010_pppp[223];

    auto g_x_0_z_0_z_x_z_z = buffer_1010_pppp[224];

    auto g_x_0_z_0_z_y_x_x = buffer_1010_pppp[225];

    auto g_x_0_z_0_z_y_x_y = buffer_1010_pppp[226];

    auto g_x_0_z_0_z_y_x_z = buffer_1010_pppp[227];

    auto g_x_0_z_0_z_y_y_x = buffer_1010_pppp[228];

    auto g_x_0_z_0_z_y_y_y = buffer_1010_pppp[229];

    auto g_x_0_z_0_z_y_y_z = buffer_1010_pppp[230];

    auto g_x_0_z_0_z_y_z_x = buffer_1010_pppp[231];

    auto g_x_0_z_0_z_y_z_y = buffer_1010_pppp[232];

    auto g_x_0_z_0_z_y_z_z = buffer_1010_pppp[233];

    auto g_x_0_z_0_z_z_x_x = buffer_1010_pppp[234];

    auto g_x_0_z_0_z_z_x_y = buffer_1010_pppp[235];

    auto g_x_0_z_0_z_z_x_z = buffer_1010_pppp[236];

    auto g_x_0_z_0_z_z_y_x = buffer_1010_pppp[237];

    auto g_x_0_z_0_z_z_y_y = buffer_1010_pppp[238];

    auto g_x_0_z_0_z_z_y_z = buffer_1010_pppp[239];

    auto g_x_0_z_0_z_z_z_x = buffer_1010_pppp[240];

    auto g_x_0_z_0_z_z_z_y = buffer_1010_pppp[241];

    auto g_x_0_z_0_z_z_z_z = buffer_1010_pppp[242];

    auto g_y_0_x_0_x_x_x_x = buffer_1010_pppp[243];

    auto g_y_0_x_0_x_x_x_y = buffer_1010_pppp[244];

    auto g_y_0_x_0_x_x_x_z = buffer_1010_pppp[245];

    auto g_y_0_x_0_x_x_y_x = buffer_1010_pppp[246];

    auto g_y_0_x_0_x_x_y_y = buffer_1010_pppp[247];

    auto g_y_0_x_0_x_x_y_z = buffer_1010_pppp[248];

    auto g_y_0_x_0_x_x_z_x = buffer_1010_pppp[249];

    auto g_y_0_x_0_x_x_z_y = buffer_1010_pppp[250];

    auto g_y_0_x_0_x_x_z_z = buffer_1010_pppp[251];

    auto g_y_0_x_0_x_y_x_x = buffer_1010_pppp[252];

    auto g_y_0_x_0_x_y_x_y = buffer_1010_pppp[253];

    auto g_y_0_x_0_x_y_x_z = buffer_1010_pppp[254];

    auto g_y_0_x_0_x_y_y_x = buffer_1010_pppp[255];

    auto g_y_0_x_0_x_y_y_y = buffer_1010_pppp[256];

    auto g_y_0_x_0_x_y_y_z = buffer_1010_pppp[257];

    auto g_y_0_x_0_x_y_z_x = buffer_1010_pppp[258];

    auto g_y_0_x_0_x_y_z_y = buffer_1010_pppp[259];

    auto g_y_0_x_0_x_y_z_z = buffer_1010_pppp[260];

    auto g_y_0_x_0_x_z_x_x = buffer_1010_pppp[261];

    auto g_y_0_x_0_x_z_x_y = buffer_1010_pppp[262];

    auto g_y_0_x_0_x_z_x_z = buffer_1010_pppp[263];

    auto g_y_0_x_0_x_z_y_x = buffer_1010_pppp[264];

    auto g_y_0_x_0_x_z_y_y = buffer_1010_pppp[265];

    auto g_y_0_x_0_x_z_y_z = buffer_1010_pppp[266];

    auto g_y_0_x_0_x_z_z_x = buffer_1010_pppp[267];

    auto g_y_0_x_0_x_z_z_y = buffer_1010_pppp[268];

    auto g_y_0_x_0_x_z_z_z = buffer_1010_pppp[269];

    auto g_y_0_x_0_y_x_x_x = buffer_1010_pppp[270];

    auto g_y_0_x_0_y_x_x_y = buffer_1010_pppp[271];

    auto g_y_0_x_0_y_x_x_z = buffer_1010_pppp[272];

    auto g_y_0_x_0_y_x_y_x = buffer_1010_pppp[273];

    auto g_y_0_x_0_y_x_y_y = buffer_1010_pppp[274];

    auto g_y_0_x_0_y_x_y_z = buffer_1010_pppp[275];

    auto g_y_0_x_0_y_x_z_x = buffer_1010_pppp[276];

    auto g_y_0_x_0_y_x_z_y = buffer_1010_pppp[277];

    auto g_y_0_x_0_y_x_z_z = buffer_1010_pppp[278];

    auto g_y_0_x_0_y_y_x_x = buffer_1010_pppp[279];

    auto g_y_0_x_0_y_y_x_y = buffer_1010_pppp[280];

    auto g_y_0_x_0_y_y_x_z = buffer_1010_pppp[281];

    auto g_y_0_x_0_y_y_y_x = buffer_1010_pppp[282];

    auto g_y_0_x_0_y_y_y_y = buffer_1010_pppp[283];

    auto g_y_0_x_0_y_y_y_z = buffer_1010_pppp[284];

    auto g_y_0_x_0_y_y_z_x = buffer_1010_pppp[285];

    auto g_y_0_x_0_y_y_z_y = buffer_1010_pppp[286];

    auto g_y_0_x_0_y_y_z_z = buffer_1010_pppp[287];

    auto g_y_0_x_0_y_z_x_x = buffer_1010_pppp[288];

    auto g_y_0_x_0_y_z_x_y = buffer_1010_pppp[289];

    auto g_y_0_x_0_y_z_x_z = buffer_1010_pppp[290];

    auto g_y_0_x_0_y_z_y_x = buffer_1010_pppp[291];

    auto g_y_0_x_0_y_z_y_y = buffer_1010_pppp[292];

    auto g_y_0_x_0_y_z_y_z = buffer_1010_pppp[293];

    auto g_y_0_x_0_y_z_z_x = buffer_1010_pppp[294];

    auto g_y_0_x_0_y_z_z_y = buffer_1010_pppp[295];

    auto g_y_0_x_0_y_z_z_z = buffer_1010_pppp[296];

    auto g_y_0_x_0_z_x_x_x = buffer_1010_pppp[297];

    auto g_y_0_x_0_z_x_x_y = buffer_1010_pppp[298];

    auto g_y_0_x_0_z_x_x_z = buffer_1010_pppp[299];

    auto g_y_0_x_0_z_x_y_x = buffer_1010_pppp[300];

    auto g_y_0_x_0_z_x_y_y = buffer_1010_pppp[301];

    auto g_y_0_x_0_z_x_y_z = buffer_1010_pppp[302];

    auto g_y_0_x_0_z_x_z_x = buffer_1010_pppp[303];

    auto g_y_0_x_0_z_x_z_y = buffer_1010_pppp[304];

    auto g_y_0_x_0_z_x_z_z = buffer_1010_pppp[305];

    auto g_y_0_x_0_z_y_x_x = buffer_1010_pppp[306];

    auto g_y_0_x_0_z_y_x_y = buffer_1010_pppp[307];

    auto g_y_0_x_0_z_y_x_z = buffer_1010_pppp[308];

    auto g_y_0_x_0_z_y_y_x = buffer_1010_pppp[309];

    auto g_y_0_x_0_z_y_y_y = buffer_1010_pppp[310];

    auto g_y_0_x_0_z_y_y_z = buffer_1010_pppp[311];

    auto g_y_0_x_0_z_y_z_x = buffer_1010_pppp[312];

    auto g_y_0_x_0_z_y_z_y = buffer_1010_pppp[313];

    auto g_y_0_x_0_z_y_z_z = buffer_1010_pppp[314];

    auto g_y_0_x_0_z_z_x_x = buffer_1010_pppp[315];

    auto g_y_0_x_0_z_z_x_y = buffer_1010_pppp[316];

    auto g_y_0_x_0_z_z_x_z = buffer_1010_pppp[317];

    auto g_y_0_x_0_z_z_y_x = buffer_1010_pppp[318];

    auto g_y_0_x_0_z_z_y_y = buffer_1010_pppp[319];

    auto g_y_0_x_0_z_z_y_z = buffer_1010_pppp[320];

    auto g_y_0_x_0_z_z_z_x = buffer_1010_pppp[321];

    auto g_y_0_x_0_z_z_z_y = buffer_1010_pppp[322];

    auto g_y_0_x_0_z_z_z_z = buffer_1010_pppp[323];

    auto g_y_0_y_0_x_x_x_x = buffer_1010_pppp[324];

    auto g_y_0_y_0_x_x_x_y = buffer_1010_pppp[325];

    auto g_y_0_y_0_x_x_x_z = buffer_1010_pppp[326];

    auto g_y_0_y_0_x_x_y_x = buffer_1010_pppp[327];

    auto g_y_0_y_0_x_x_y_y = buffer_1010_pppp[328];

    auto g_y_0_y_0_x_x_y_z = buffer_1010_pppp[329];

    auto g_y_0_y_0_x_x_z_x = buffer_1010_pppp[330];

    auto g_y_0_y_0_x_x_z_y = buffer_1010_pppp[331];

    auto g_y_0_y_0_x_x_z_z = buffer_1010_pppp[332];

    auto g_y_0_y_0_x_y_x_x = buffer_1010_pppp[333];

    auto g_y_0_y_0_x_y_x_y = buffer_1010_pppp[334];

    auto g_y_0_y_0_x_y_x_z = buffer_1010_pppp[335];

    auto g_y_0_y_0_x_y_y_x = buffer_1010_pppp[336];

    auto g_y_0_y_0_x_y_y_y = buffer_1010_pppp[337];

    auto g_y_0_y_0_x_y_y_z = buffer_1010_pppp[338];

    auto g_y_0_y_0_x_y_z_x = buffer_1010_pppp[339];

    auto g_y_0_y_0_x_y_z_y = buffer_1010_pppp[340];

    auto g_y_0_y_0_x_y_z_z = buffer_1010_pppp[341];

    auto g_y_0_y_0_x_z_x_x = buffer_1010_pppp[342];

    auto g_y_0_y_0_x_z_x_y = buffer_1010_pppp[343];

    auto g_y_0_y_0_x_z_x_z = buffer_1010_pppp[344];

    auto g_y_0_y_0_x_z_y_x = buffer_1010_pppp[345];

    auto g_y_0_y_0_x_z_y_y = buffer_1010_pppp[346];

    auto g_y_0_y_0_x_z_y_z = buffer_1010_pppp[347];

    auto g_y_0_y_0_x_z_z_x = buffer_1010_pppp[348];

    auto g_y_0_y_0_x_z_z_y = buffer_1010_pppp[349];

    auto g_y_0_y_0_x_z_z_z = buffer_1010_pppp[350];

    auto g_y_0_y_0_y_x_x_x = buffer_1010_pppp[351];

    auto g_y_0_y_0_y_x_x_y = buffer_1010_pppp[352];

    auto g_y_0_y_0_y_x_x_z = buffer_1010_pppp[353];

    auto g_y_0_y_0_y_x_y_x = buffer_1010_pppp[354];

    auto g_y_0_y_0_y_x_y_y = buffer_1010_pppp[355];

    auto g_y_0_y_0_y_x_y_z = buffer_1010_pppp[356];

    auto g_y_0_y_0_y_x_z_x = buffer_1010_pppp[357];

    auto g_y_0_y_0_y_x_z_y = buffer_1010_pppp[358];

    auto g_y_0_y_0_y_x_z_z = buffer_1010_pppp[359];

    auto g_y_0_y_0_y_y_x_x = buffer_1010_pppp[360];

    auto g_y_0_y_0_y_y_x_y = buffer_1010_pppp[361];

    auto g_y_0_y_0_y_y_x_z = buffer_1010_pppp[362];

    auto g_y_0_y_0_y_y_y_x = buffer_1010_pppp[363];

    auto g_y_0_y_0_y_y_y_y = buffer_1010_pppp[364];

    auto g_y_0_y_0_y_y_y_z = buffer_1010_pppp[365];

    auto g_y_0_y_0_y_y_z_x = buffer_1010_pppp[366];

    auto g_y_0_y_0_y_y_z_y = buffer_1010_pppp[367];

    auto g_y_0_y_0_y_y_z_z = buffer_1010_pppp[368];

    auto g_y_0_y_0_y_z_x_x = buffer_1010_pppp[369];

    auto g_y_0_y_0_y_z_x_y = buffer_1010_pppp[370];

    auto g_y_0_y_0_y_z_x_z = buffer_1010_pppp[371];

    auto g_y_0_y_0_y_z_y_x = buffer_1010_pppp[372];

    auto g_y_0_y_0_y_z_y_y = buffer_1010_pppp[373];

    auto g_y_0_y_0_y_z_y_z = buffer_1010_pppp[374];

    auto g_y_0_y_0_y_z_z_x = buffer_1010_pppp[375];

    auto g_y_0_y_0_y_z_z_y = buffer_1010_pppp[376];

    auto g_y_0_y_0_y_z_z_z = buffer_1010_pppp[377];

    auto g_y_0_y_0_z_x_x_x = buffer_1010_pppp[378];

    auto g_y_0_y_0_z_x_x_y = buffer_1010_pppp[379];

    auto g_y_0_y_0_z_x_x_z = buffer_1010_pppp[380];

    auto g_y_0_y_0_z_x_y_x = buffer_1010_pppp[381];

    auto g_y_0_y_0_z_x_y_y = buffer_1010_pppp[382];

    auto g_y_0_y_0_z_x_y_z = buffer_1010_pppp[383];

    auto g_y_0_y_0_z_x_z_x = buffer_1010_pppp[384];

    auto g_y_0_y_0_z_x_z_y = buffer_1010_pppp[385];

    auto g_y_0_y_0_z_x_z_z = buffer_1010_pppp[386];

    auto g_y_0_y_0_z_y_x_x = buffer_1010_pppp[387];

    auto g_y_0_y_0_z_y_x_y = buffer_1010_pppp[388];

    auto g_y_0_y_0_z_y_x_z = buffer_1010_pppp[389];

    auto g_y_0_y_0_z_y_y_x = buffer_1010_pppp[390];

    auto g_y_0_y_0_z_y_y_y = buffer_1010_pppp[391];

    auto g_y_0_y_0_z_y_y_z = buffer_1010_pppp[392];

    auto g_y_0_y_0_z_y_z_x = buffer_1010_pppp[393];

    auto g_y_0_y_0_z_y_z_y = buffer_1010_pppp[394];

    auto g_y_0_y_0_z_y_z_z = buffer_1010_pppp[395];

    auto g_y_0_y_0_z_z_x_x = buffer_1010_pppp[396];

    auto g_y_0_y_0_z_z_x_y = buffer_1010_pppp[397];

    auto g_y_0_y_0_z_z_x_z = buffer_1010_pppp[398];

    auto g_y_0_y_0_z_z_y_x = buffer_1010_pppp[399];

    auto g_y_0_y_0_z_z_y_y = buffer_1010_pppp[400];

    auto g_y_0_y_0_z_z_y_z = buffer_1010_pppp[401];

    auto g_y_0_y_0_z_z_z_x = buffer_1010_pppp[402];

    auto g_y_0_y_0_z_z_z_y = buffer_1010_pppp[403];

    auto g_y_0_y_0_z_z_z_z = buffer_1010_pppp[404];

    auto g_y_0_z_0_x_x_x_x = buffer_1010_pppp[405];

    auto g_y_0_z_0_x_x_x_y = buffer_1010_pppp[406];

    auto g_y_0_z_0_x_x_x_z = buffer_1010_pppp[407];

    auto g_y_0_z_0_x_x_y_x = buffer_1010_pppp[408];

    auto g_y_0_z_0_x_x_y_y = buffer_1010_pppp[409];

    auto g_y_0_z_0_x_x_y_z = buffer_1010_pppp[410];

    auto g_y_0_z_0_x_x_z_x = buffer_1010_pppp[411];

    auto g_y_0_z_0_x_x_z_y = buffer_1010_pppp[412];

    auto g_y_0_z_0_x_x_z_z = buffer_1010_pppp[413];

    auto g_y_0_z_0_x_y_x_x = buffer_1010_pppp[414];

    auto g_y_0_z_0_x_y_x_y = buffer_1010_pppp[415];

    auto g_y_0_z_0_x_y_x_z = buffer_1010_pppp[416];

    auto g_y_0_z_0_x_y_y_x = buffer_1010_pppp[417];

    auto g_y_0_z_0_x_y_y_y = buffer_1010_pppp[418];

    auto g_y_0_z_0_x_y_y_z = buffer_1010_pppp[419];

    auto g_y_0_z_0_x_y_z_x = buffer_1010_pppp[420];

    auto g_y_0_z_0_x_y_z_y = buffer_1010_pppp[421];

    auto g_y_0_z_0_x_y_z_z = buffer_1010_pppp[422];

    auto g_y_0_z_0_x_z_x_x = buffer_1010_pppp[423];

    auto g_y_0_z_0_x_z_x_y = buffer_1010_pppp[424];

    auto g_y_0_z_0_x_z_x_z = buffer_1010_pppp[425];

    auto g_y_0_z_0_x_z_y_x = buffer_1010_pppp[426];

    auto g_y_0_z_0_x_z_y_y = buffer_1010_pppp[427];

    auto g_y_0_z_0_x_z_y_z = buffer_1010_pppp[428];

    auto g_y_0_z_0_x_z_z_x = buffer_1010_pppp[429];

    auto g_y_0_z_0_x_z_z_y = buffer_1010_pppp[430];

    auto g_y_0_z_0_x_z_z_z = buffer_1010_pppp[431];

    auto g_y_0_z_0_y_x_x_x = buffer_1010_pppp[432];

    auto g_y_0_z_0_y_x_x_y = buffer_1010_pppp[433];

    auto g_y_0_z_0_y_x_x_z = buffer_1010_pppp[434];

    auto g_y_0_z_0_y_x_y_x = buffer_1010_pppp[435];

    auto g_y_0_z_0_y_x_y_y = buffer_1010_pppp[436];

    auto g_y_0_z_0_y_x_y_z = buffer_1010_pppp[437];

    auto g_y_0_z_0_y_x_z_x = buffer_1010_pppp[438];

    auto g_y_0_z_0_y_x_z_y = buffer_1010_pppp[439];

    auto g_y_0_z_0_y_x_z_z = buffer_1010_pppp[440];

    auto g_y_0_z_0_y_y_x_x = buffer_1010_pppp[441];

    auto g_y_0_z_0_y_y_x_y = buffer_1010_pppp[442];

    auto g_y_0_z_0_y_y_x_z = buffer_1010_pppp[443];

    auto g_y_0_z_0_y_y_y_x = buffer_1010_pppp[444];

    auto g_y_0_z_0_y_y_y_y = buffer_1010_pppp[445];

    auto g_y_0_z_0_y_y_y_z = buffer_1010_pppp[446];

    auto g_y_0_z_0_y_y_z_x = buffer_1010_pppp[447];

    auto g_y_0_z_0_y_y_z_y = buffer_1010_pppp[448];

    auto g_y_0_z_0_y_y_z_z = buffer_1010_pppp[449];

    auto g_y_0_z_0_y_z_x_x = buffer_1010_pppp[450];

    auto g_y_0_z_0_y_z_x_y = buffer_1010_pppp[451];

    auto g_y_0_z_0_y_z_x_z = buffer_1010_pppp[452];

    auto g_y_0_z_0_y_z_y_x = buffer_1010_pppp[453];

    auto g_y_0_z_0_y_z_y_y = buffer_1010_pppp[454];

    auto g_y_0_z_0_y_z_y_z = buffer_1010_pppp[455];

    auto g_y_0_z_0_y_z_z_x = buffer_1010_pppp[456];

    auto g_y_0_z_0_y_z_z_y = buffer_1010_pppp[457];

    auto g_y_0_z_0_y_z_z_z = buffer_1010_pppp[458];

    auto g_y_0_z_0_z_x_x_x = buffer_1010_pppp[459];

    auto g_y_0_z_0_z_x_x_y = buffer_1010_pppp[460];

    auto g_y_0_z_0_z_x_x_z = buffer_1010_pppp[461];

    auto g_y_0_z_0_z_x_y_x = buffer_1010_pppp[462];

    auto g_y_0_z_0_z_x_y_y = buffer_1010_pppp[463];

    auto g_y_0_z_0_z_x_y_z = buffer_1010_pppp[464];

    auto g_y_0_z_0_z_x_z_x = buffer_1010_pppp[465];

    auto g_y_0_z_0_z_x_z_y = buffer_1010_pppp[466];

    auto g_y_0_z_0_z_x_z_z = buffer_1010_pppp[467];

    auto g_y_0_z_0_z_y_x_x = buffer_1010_pppp[468];

    auto g_y_0_z_0_z_y_x_y = buffer_1010_pppp[469];

    auto g_y_0_z_0_z_y_x_z = buffer_1010_pppp[470];

    auto g_y_0_z_0_z_y_y_x = buffer_1010_pppp[471];

    auto g_y_0_z_0_z_y_y_y = buffer_1010_pppp[472];

    auto g_y_0_z_0_z_y_y_z = buffer_1010_pppp[473];

    auto g_y_0_z_0_z_y_z_x = buffer_1010_pppp[474];

    auto g_y_0_z_0_z_y_z_y = buffer_1010_pppp[475];

    auto g_y_0_z_0_z_y_z_z = buffer_1010_pppp[476];

    auto g_y_0_z_0_z_z_x_x = buffer_1010_pppp[477];

    auto g_y_0_z_0_z_z_x_y = buffer_1010_pppp[478];

    auto g_y_0_z_0_z_z_x_z = buffer_1010_pppp[479];

    auto g_y_0_z_0_z_z_y_x = buffer_1010_pppp[480];

    auto g_y_0_z_0_z_z_y_y = buffer_1010_pppp[481];

    auto g_y_0_z_0_z_z_y_z = buffer_1010_pppp[482];

    auto g_y_0_z_0_z_z_z_x = buffer_1010_pppp[483];

    auto g_y_0_z_0_z_z_z_y = buffer_1010_pppp[484];

    auto g_y_0_z_0_z_z_z_z = buffer_1010_pppp[485];

    auto g_z_0_x_0_x_x_x_x = buffer_1010_pppp[486];

    auto g_z_0_x_0_x_x_x_y = buffer_1010_pppp[487];

    auto g_z_0_x_0_x_x_x_z = buffer_1010_pppp[488];

    auto g_z_0_x_0_x_x_y_x = buffer_1010_pppp[489];

    auto g_z_0_x_0_x_x_y_y = buffer_1010_pppp[490];

    auto g_z_0_x_0_x_x_y_z = buffer_1010_pppp[491];

    auto g_z_0_x_0_x_x_z_x = buffer_1010_pppp[492];

    auto g_z_0_x_0_x_x_z_y = buffer_1010_pppp[493];

    auto g_z_0_x_0_x_x_z_z = buffer_1010_pppp[494];

    auto g_z_0_x_0_x_y_x_x = buffer_1010_pppp[495];

    auto g_z_0_x_0_x_y_x_y = buffer_1010_pppp[496];

    auto g_z_0_x_0_x_y_x_z = buffer_1010_pppp[497];

    auto g_z_0_x_0_x_y_y_x = buffer_1010_pppp[498];

    auto g_z_0_x_0_x_y_y_y = buffer_1010_pppp[499];

    auto g_z_0_x_0_x_y_y_z = buffer_1010_pppp[500];

    auto g_z_0_x_0_x_y_z_x = buffer_1010_pppp[501];

    auto g_z_0_x_0_x_y_z_y = buffer_1010_pppp[502];

    auto g_z_0_x_0_x_y_z_z = buffer_1010_pppp[503];

    auto g_z_0_x_0_x_z_x_x = buffer_1010_pppp[504];

    auto g_z_0_x_0_x_z_x_y = buffer_1010_pppp[505];

    auto g_z_0_x_0_x_z_x_z = buffer_1010_pppp[506];

    auto g_z_0_x_0_x_z_y_x = buffer_1010_pppp[507];

    auto g_z_0_x_0_x_z_y_y = buffer_1010_pppp[508];

    auto g_z_0_x_0_x_z_y_z = buffer_1010_pppp[509];

    auto g_z_0_x_0_x_z_z_x = buffer_1010_pppp[510];

    auto g_z_0_x_0_x_z_z_y = buffer_1010_pppp[511];

    auto g_z_0_x_0_x_z_z_z = buffer_1010_pppp[512];

    auto g_z_0_x_0_y_x_x_x = buffer_1010_pppp[513];

    auto g_z_0_x_0_y_x_x_y = buffer_1010_pppp[514];

    auto g_z_0_x_0_y_x_x_z = buffer_1010_pppp[515];

    auto g_z_0_x_0_y_x_y_x = buffer_1010_pppp[516];

    auto g_z_0_x_0_y_x_y_y = buffer_1010_pppp[517];

    auto g_z_0_x_0_y_x_y_z = buffer_1010_pppp[518];

    auto g_z_0_x_0_y_x_z_x = buffer_1010_pppp[519];

    auto g_z_0_x_0_y_x_z_y = buffer_1010_pppp[520];

    auto g_z_0_x_0_y_x_z_z = buffer_1010_pppp[521];

    auto g_z_0_x_0_y_y_x_x = buffer_1010_pppp[522];

    auto g_z_0_x_0_y_y_x_y = buffer_1010_pppp[523];

    auto g_z_0_x_0_y_y_x_z = buffer_1010_pppp[524];

    auto g_z_0_x_0_y_y_y_x = buffer_1010_pppp[525];

    auto g_z_0_x_0_y_y_y_y = buffer_1010_pppp[526];

    auto g_z_0_x_0_y_y_y_z = buffer_1010_pppp[527];

    auto g_z_0_x_0_y_y_z_x = buffer_1010_pppp[528];

    auto g_z_0_x_0_y_y_z_y = buffer_1010_pppp[529];

    auto g_z_0_x_0_y_y_z_z = buffer_1010_pppp[530];

    auto g_z_0_x_0_y_z_x_x = buffer_1010_pppp[531];

    auto g_z_0_x_0_y_z_x_y = buffer_1010_pppp[532];

    auto g_z_0_x_0_y_z_x_z = buffer_1010_pppp[533];

    auto g_z_0_x_0_y_z_y_x = buffer_1010_pppp[534];

    auto g_z_0_x_0_y_z_y_y = buffer_1010_pppp[535];

    auto g_z_0_x_0_y_z_y_z = buffer_1010_pppp[536];

    auto g_z_0_x_0_y_z_z_x = buffer_1010_pppp[537];

    auto g_z_0_x_0_y_z_z_y = buffer_1010_pppp[538];

    auto g_z_0_x_0_y_z_z_z = buffer_1010_pppp[539];

    auto g_z_0_x_0_z_x_x_x = buffer_1010_pppp[540];

    auto g_z_0_x_0_z_x_x_y = buffer_1010_pppp[541];

    auto g_z_0_x_0_z_x_x_z = buffer_1010_pppp[542];

    auto g_z_0_x_0_z_x_y_x = buffer_1010_pppp[543];

    auto g_z_0_x_0_z_x_y_y = buffer_1010_pppp[544];

    auto g_z_0_x_0_z_x_y_z = buffer_1010_pppp[545];

    auto g_z_0_x_0_z_x_z_x = buffer_1010_pppp[546];

    auto g_z_0_x_0_z_x_z_y = buffer_1010_pppp[547];

    auto g_z_0_x_0_z_x_z_z = buffer_1010_pppp[548];

    auto g_z_0_x_0_z_y_x_x = buffer_1010_pppp[549];

    auto g_z_0_x_0_z_y_x_y = buffer_1010_pppp[550];

    auto g_z_0_x_0_z_y_x_z = buffer_1010_pppp[551];

    auto g_z_0_x_0_z_y_y_x = buffer_1010_pppp[552];

    auto g_z_0_x_0_z_y_y_y = buffer_1010_pppp[553];

    auto g_z_0_x_0_z_y_y_z = buffer_1010_pppp[554];

    auto g_z_0_x_0_z_y_z_x = buffer_1010_pppp[555];

    auto g_z_0_x_0_z_y_z_y = buffer_1010_pppp[556];

    auto g_z_0_x_0_z_y_z_z = buffer_1010_pppp[557];

    auto g_z_0_x_0_z_z_x_x = buffer_1010_pppp[558];

    auto g_z_0_x_0_z_z_x_y = buffer_1010_pppp[559];

    auto g_z_0_x_0_z_z_x_z = buffer_1010_pppp[560];

    auto g_z_0_x_0_z_z_y_x = buffer_1010_pppp[561];

    auto g_z_0_x_0_z_z_y_y = buffer_1010_pppp[562];

    auto g_z_0_x_0_z_z_y_z = buffer_1010_pppp[563];

    auto g_z_0_x_0_z_z_z_x = buffer_1010_pppp[564];

    auto g_z_0_x_0_z_z_z_y = buffer_1010_pppp[565];

    auto g_z_0_x_0_z_z_z_z = buffer_1010_pppp[566];

    auto g_z_0_y_0_x_x_x_x = buffer_1010_pppp[567];

    auto g_z_0_y_0_x_x_x_y = buffer_1010_pppp[568];

    auto g_z_0_y_0_x_x_x_z = buffer_1010_pppp[569];

    auto g_z_0_y_0_x_x_y_x = buffer_1010_pppp[570];

    auto g_z_0_y_0_x_x_y_y = buffer_1010_pppp[571];

    auto g_z_0_y_0_x_x_y_z = buffer_1010_pppp[572];

    auto g_z_0_y_0_x_x_z_x = buffer_1010_pppp[573];

    auto g_z_0_y_0_x_x_z_y = buffer_1010_pppp[574];

    auto g_z_0_y_0_x_x_z_z = buffer_1010_pppp[575];

    auto g_z_0_y_0_x_y_x_x = buffer_1010_pppp[576];

    auto g_z_0_y_0_x_y_x_y = buffer_1010_pppp[577];

    auto g_z_0_y_0_x_y_x_z = buffer_1010_pppp[578];

    auto g_z_0_y_0_x_y_y_x = buffer_1010_pppp[579];

    auto g_z_0_y_0_x_y_y_y = buffer_1010_pppp[580];

    auto g_z_0_y_0_x_y_y_z = buffer_1010_pppp[581];

    auto g_z_0_y_0_x_y_z_x = buffer_1010_pppp[582];

    auto g_z_0_y_0_x_y_z_y = buffer_1010_pppp[583];

    auto g_z_0_y_0_x_y_z_z = buffer_1010_pppp[584];

    auto g_z_0_y_0_x_z_x_x = buffer_1010_pppp[585];

    auto g_z_0_y_0_x_z_x_y = buffer_1010_pppp[586];

    auto g_z_0_y_0_x_z_x_z = buffer_1010_pppp[587];

    auto g_z_0_y_0_x_z_y_x = buffer_1010_pppp[588];

    auto g_z_0_y_0_x_z_y_y = buffer_1010_pppp[589];

    auto g_z_0_y_0_x_z_y_z = buffer_1010_pppp[590];

    auto g_z_0_y_0_x_z_z_x = buffer_1010_pppp[591];

    auto g_z_0_y_0_x_z_z_y = buffer_1010_pppp[592];

    auto g_z_0_y_0_x_z_z_z = buffer_1010_pppp[593];

    auto g_z_0_y_0_y_x_x_x = buffer_1010_pppp[594];

    auto g_z_0_y_0_y_x_x_y = buffer_1010_pppp[595];

    auto g_z_0_y_0_y_x_x_z = buffer_1010_pppp[596];

    auto g_z_0_y_0_y_x_y_x = buffer_1010_pppp[597];

    auto g_z_0_y_0_y_x_y_y = buffer_1010_pppp[598];

    auto g_z_0_y_0_y_x_y_z = buffer_1010_pppp[599];

    auto g_z_0_y_0_y_x_z_x = buffer_1010_pppp[600];

    auto g_z_0_y_0_y_x_z_y = buffer_1010_pppp[601];

    auto g_z_0_y_0_y_x_z_z = buffer_1010_pppp[602];

    auto g_z_0_y_0_y_y_x_x = buffer_1010_pppp[603];

    auto g_z_0_y_0_y_y_x_y = buffer_1010_pppp[604];

    auto g_z_0_y_0_y_y_x_z = buffer_1010_pppp[605];

    auto g_z_0_y_0_y_y_y_x = buffer_1010_pppp[606];

    auto g_z_0_y_0_y_y_y_y = buffer_1010_pppp[607];

    auto g_z_0_y_0_y_y_y_z = buffer_1010_pppp[608];

    auto g_z_0_y_0_y_y_z_x = buffer_1010_pppp[609];

    auto g_z_0_y_0_y_y_z_y = buffer_1010_pppp[610];

    auto g_z_0_y_0_y_y_z_z = buffer_1010_pppp[611];

    auto g_z_0_y_0_y_z_x_x = buffer_1010_pppp[612];

    auto g_z_0_y_0_y_z_x_y = buffer_1010_pppp[613];

    auto g_z_0_y_0_y_z_x_z = buffer_1010_pppp[614];

    auto g_z_0_y_0_y_z_y_x = buffer_1010_pppp[615];

    auto g_z_0_y_0_y_z_y_y = buffer_1010_pppp[616];

    auto g_z_0_y_0_y_z_y_z = buffer_1010_pppp[617];

    auto g_z_0_y_0_y_z_z_x = buffer_1010_pppp[618];

    auto g_z_0_y_0_y_z_z_y = buffer_1010_pppp[619];

    auto g_z_0_y_0_y_z_z_z = buffer_1010_pppp[620];

    auto g_z_0_y_0_z_x_x_x = buffer_1010_pppp[621];

    auto g_z_0_y_0_z_x_x_y = buffer_1010_pppp[622];

    auto g_z_0_y_0_z_x_x_z = buffer_1010_pppp[623];

    auto g_z_0_y_0_z_x_y_x = buffer_1010_pppp[624];

    auto g_z_0_y_0_z_x_y_y = buffer_1010_pppp[625];

    auto g_z_0_y_0_z_x_y_z = buffer_1010_pppp[626];

    auto g_z_0_y_0_z_x_z_x = buffer_1010_pppp[627];

    auto g_z_0_y_0_z_x_z_y = buffer_1010_pppp[628];

    auto g_z_0_y_0_z_x_z_z = buffer_1010_pppp[629];

    auto g_z_0_y_0_z_y_x_x = buffer_1010_pppp[630];

    auto g_z_0_y_0_z_y_x_y = buffer_1010_pppp[631];

    auto g_z_0_y_0_z_y_x_z = buffer_1010_pppp[632];

    auto g_z_0_y_0_z_y_y_x = buffer_1010_pppp[633];

    auto g_z_0_y_0_z_y_y_y = buffer_1010_pppp[634];

    auto g_z_0_y_0_z_y_y_z = buffer_1010_pppp[635];

    auto g_z_0_y_0_z_y_z_x = buffer_1010_pppp[636];

    auto g_z_0_y_0_z_y_z_y = buffer_1010_pppp[637];

    auto g_z_0_y_0_z_y_z_z = buffer_1010_pppp[638];

    auto g_z_0_y_0_z_z_x_x = buffer_1010_pppp[639];

    auto g_z_0_y_0_z_z_x_y = buffer_1010_pppp[640];

    auto g_z_0_y_0_z_z_x_z = buffer_1010_pppp[641];

    auto g_z_0_y_0_z_z_y_x = buffer_1010_pppp[642];

    auto g_z_0_y_0_z_z_y_y = buffer_1010_pppp[643];

    auto g_z_0_y_0_z_z_y_z = buffer_1010_pppp[644];

    auto g_z_0_y_0_z_z_z_x = buffer_1010_pppp[645];

    auto g_z_0_y_0_z_z_z_y = buffer_1010_pppp[646];

    auto g_z_0_y_0_z_z_z_z = buffer_1010_pppp[647];

    auto g_z_0_z_0_x_x_x_x = buffer_1010_pppp[648];

    auto g_z_0_z_0_x_x_x_y = buffer_1010_pppp[649];

    auto g_z_0_z_0_x_x_x_z = buffer_1010_pppp[650];

    auto g_z_0_z_0_x_x_y_x = buffer_1010_pppp[651];

    auto g_z_0_z_0_x_x_y_y = buffer_1010_pppp[652];

    auto g_z_0_z_0_x_x_y_z = buffer_1010_pppp[653];

    auto g_z_0_z_0_x_x_z_x = buffer_1010_pppp[654];

    auto g_z_0_z_0_x_x_z_y = buffer_1010_pppp[655];

    auto g_z_0_z_0_x_x_z_z = buffer_1010_pppp[656];

    auto g_z_0_z_0_x_y_x_x = buffer_1010_pppp[657];

    auto g_z_0_z_0_x_y_x_y = buffer_1010_pppp[658];

    auto g_z_0_z_0_x_y_x_z = buffer_1010_pppp[659];

    auto g_z_0_z_0_x_y_y_x = buffer_1010_pppp[660];

    auto g_z_0_z_0_x_y_y_y = buffer_1010_pppp[661];

    auto g_z_0_z_0_x_y_y_z = buffer_1010_pppp[662];

    auto g_z_0_z_0_x_y_z_x = buffer_1010_pppp[663];

    auto g_z_0_z_0_x_y_z_y = buffer_1010_pppp[664];

    auto g_z_0_z_0_x_y_z_z = buffer_1010_pppp[665];

    auto g_z_0_z_0_x_z_x_x = buffer_1010_pppp[666];

    auto g_z_0_z_0_x_z_x_y = buffer_1010_pppp[667];

    auto g_z_0_z_0_x_z_x_z = buffer_1010_pppp[668];

    auto g_z_0_z_0_x_z_y_x = buffer_1010_pppp[669];

    auto g_z_0_z_0_x_z_y_y = buffer_1010_pppp[670];

    auto g_z_0_z_0_x_z_y_z = buffer_1010_pppp[671];

    auto g_z_0_z_0_x_z_z_x = buffer_1010_pppp[672];

    auto g_z_0_z_0_x_z_z_y = buffer_1010_pppp[673];

    auto g_z_0_z_0_x_z_z_z = buffer_1010_pppp[674];

    auto g_z_0_z_0_y_x_x_x = buffer_1010_pppp[675];

    auto g_z_0_z_0_y_x_x_y = buffer_1010_pppp[676];

    auto g_z_0_z_0_y_x_x_z = buffer_1010_pppp[677];

    auto g_z_0_z_0_y_x_y_x = buffer_1010_pppp[678];

    auto g_z_0_z_0_y_x_y_y = buffer_1010_pppp[679];

    auto g_z_0_z_0_y_x_y_z = buffer_1010_pppp[680];

    auto g_z_0_z_0_y_x_z_x = buffer_1010_pppp[681];

    auto g_z_0_z_0_y_x_z_y = buffer_1010_pppp[682];

    auto g_z_0_z_0_y_x_z_z = buffer_1010_pppp[683];

    auto g_z_0_z_0_y_y_x_x = buffer_1010_pppp[684];

    auto g_z_0_z_0_y_y_x_y = buffer_1010_pppp[685];

    auto g_z_0_z_0_y_y_x_z = buffer_1010_pppp[686];

    auto g_z_0_z_0_y_y_y_x = buffer_1010_pppp[687];

    auto g_z_0_z_0_y_y_y_y = buffer_1010_pppp[688];

    auto g_z_0_z_0_y_y_y_z = buffer_1010_pppp[689];

    auto g_z_0_z_0_y_y_z_x = buffer_1010_pppp[690];

    auto g_z_0_z_0_y_y_z_y = buffer_1010_pppp[691];

    auto g_z_0_z_0_y_y_z_z = buffer_1010_pppp[692];

    auto g_z_0_z_0_y_z_x_x = buffer_1010_pppp[693];

    auto g_z_0_z_0_y_z_x_y = buffer_1010_pppp[694];

    auto g_z_0_z_0_y_z_x_z = buffer_1010_pppp[695];

    auto g_z_0_z_0_y_z_y_x = buffer_1010_pppp[696];

    auto g_z_0_z_0_y_z_y_y = buffer_1010_pppp[697];

    auto g_z_0_z_0_y_z_y_z = buffer_1010_pppp[698];

    auto g_z_0_z_0_y_z_z_x = buffer_1010_pppp[699];

    auto g_z_0_z_0_y_z_z_y = buffer_1010_pppp[700];

    auto g_z_0_z_0_y_z_z_z = buffer_1010_pppp[701];

    auto g_z_0_z_0_z_x_x_x = buffer_1010_pppp[702];

    auto g_z_0_z_0_z_x_x_y = buffer_1010_pppp[703];

    auto g_z_0_z_0_z_x_x_z = buffer_1010_pppp[704];

    auto g_z_0_z_0_z_x_y_x = buffer_1010_pppp[705];

    auto g_z_0_z_0_z_x_y_y = buffer_1010_pppp[706];

    auto g_z_0_z_0_z_x_y_z = buffer_1010_pppp[707];

    auto g_z_0_z_0_z_x_z_x = buffer_1010_pppp[708];

    auto g_z_0_z_0_z_x_z_y = buffer_1010_pppp[709];

    auto g_z_0_z_0_z_x_z_z = buffer_1010_pppp[710];

    auto g_z_0_z_0_z_y_x_x = buffer_1010_pppp[711];

    auto g_z_0_z_0_z_y_x_y = buffer_1010_pppp[712];

    auto g_z_0_z_0_z_y_x_z = buffer_1010_pppp[713];

    auto g_z_0_z_0_z_y_y_x = buffer_1010_pppp[714];

    auto g_z_0_z_0_z_y_y_y = buffer_1010_pppp[715];

    auto g_z_0_z_0_z_y_y_z = buffer_1010_pppp[716];

    auto g_z_0_z_0_z_y_z_x = buffer_1010_pppp[717];

    auto g_z_0_z_0_z_y_z_y = buffer_1010_pppp[718];

    auto g_z_0_z_0_z_y_z_z = buffer_1010_pppp[719];

    auto g_z_0_z_0_z_z_x_x = buffer_1010_pppp[720];

    auto g_z_0_z_0_z_z_x_y = buffer_1010_pppp[721];

    auto g_z_0_z_0_z_z_x_z = buffer_1010_pppp[722];

    auto g_z_0_z_0_z_z_y_x = buffer_1010_pppp[723];

    auto g_z_0_z_0_z_z_y_y = buffer_1010_pppp[724];

    auto g_z_0_z_0_z_z_y_z = buffer_1010_pppp[725];

    auto g_z_0_z_0_z_z_z_x = buffer_1010_pppp[726];

    auto g_z_0_z_0_z_z_z_y = buffer_1010_pppp[727];

    auto g_z_0_z_0_z_z_z_z = buffer_1010_pppp[728];

    // integrals block (0-3)

    #pragma omp simd aligned(g_0_x_0_x, g_0_x_0_y, g_0_x_0_z, g_0_x_xx_x, g_0_x_xx_y, g_0_x_xx_z, g_x_0_x_0_x_x_x_x, g_x_0_x_0_x_x_x_y, g_x_0_x_0_x_x_x_z, g_xx_x_0_x, g_xx_x_0_y, g_xx_x_0_z, g_xx_x_xx_x, g_xx_x_xx_y, g_xx_x_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_x_x_x[i] = g_0_x_0_x[i] - 2.0 * g_0_x_xx_x[i] * c_exps[i] - 2.0 * g_xx_x_0_x[i] * a_exp + 4.0 * g_xx_x_xx_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_x_y[i] = g_0_x_0_y[i] - 2.0 * g_0_x_xx_y[i] * c_exps[i] - 2.0 * g_xx_x_0_y[i] * a_exp + 4.0 * g_xx_x_xx_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_x_z[i] = g_0_x_0_z[i] - 2.0 * g_0_x_xx_z[i] * c_exps[i] - 2.0 * g_xx_x_0_z[i] * a_exp + 4.0 * g_xx_x_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_0_x_xy_x, g_0_x_xy_y, g_0_x_xy_z, g_x_0_x_0_x_x_y_x, g_x_0_x_0_x_x_y_y, g_x_0_x_0_x_x_y_z, g_xx_x_xy_x, g_xx_x_xy_y, g_xx_x_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_x_y_x[i] = -2.0 * g_0_x_xy_x[i] * c_exps[i] + 4.0 * g_xx_x_xy_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_y_y[i] = -2.0 * g_0_x_xy_y[i] * c_exps[i] + 4.0 * g_xx_x_xy_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_y_z[i] = -2.0 * g_0_x_xy_z[i] * c_exps[i] + 4.0 * g_xx_x_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_0_x_xz_x, g_0_x_xz_y, g_0_x_xz_z, g_x_0_x_0_x_x_z_x, g_x_0_x_0_x_x_z_y, g_x_0_x_0_x_x_z_z, g_xx_x_xz_x, g_xx_x_xz_y, g_xx_x_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_x_z_x[i] = -2.0 * g_0_x_xz_x[i] * c_exps[i] + 4.0 * g_xx_x_xz_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_z_y[i] = -2.0 * g_0_x_xz_y[i] * c_exps[i] + 4.0 * g_xx_x_xz_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_z_z[i] = -2.0 * g_0_x_xz_z[i] * c_exps[i] + 4.0 * g_xx_x_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_0_y_0_x, g_0_y_0_y, g_0_y_0_z, g_0_y_xx_x, g_0_y_xx_y, g_0_y_xx_z, g_x_0_x_0_x_y_x_x, g_x_0_x_0_x_y_x_y, g_x_0_x_0_x_y_x_z, g_xx_y_0_x, g_xx_y_0_y, g_xx_y_0_z, g_xx_y_xx_x, g_xx_y_xx_y, g_xx_y_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_y_x_x[i] = g_0_y_0_x[i] - 2.0 * g_0_y_xx_x[i] * c_exps[i] - 2.0 * g_xx_y_0_x[i] * a_exp + 4.0 * g_xx_y_xx_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_x_y[i] = g_0_y_0_y[i] - 2.0 * g_0_y_xx_y[i] * c_exps[i] - 2.0 * g_xx_y_0_y[i] * a_exp + 4.0 * g_xx_y_xx_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_x_z[i] = g_0_y_0_z[i] - 2.0 * g_0_y_xx_z[i] * c_exps[i] - 2.0 * g_xx_y_0_z[i] * a_exp + 4.0 * g_xx_y_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_0_y_xy_x, g_0_y_xy_y, g_0_y_xy_z, g_x_0_x_0_x_y_y_x, g_x_0_x_0_x_y_y_y, g_x_0_x_0_x_y_y_z, g_xx_y_xy_x, g_xx_y_xy_y, g_xx_y_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_y_y_x[i] = -2.0 * g_0_y_xy_x[i] * c_exps[i] + 4.0 * g_xx_y_xy_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_y_y[i] = -2.0 * g_0_y_xy_y[i] * c_exps[i] + 4.0 * g_xx_y_xy_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_y_z[i] = -2.0 * g_0_y_xy_z[i] * c_exps[i] + 4.0 * g_xx_y_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_0_y_xz_x, g_0_y_xz_y, g_0_y_xz_z, g_x_0_x_0_x_y_z_x, g_x_0_x_0_x_y_z_y, g_x_0_x_0_x_y_z_z, g_xx_y_xz_x, g_xx_y_xz_y, g_xx_y_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_y_z_x[i] = -2.0 * g_0_y_xz_x[i] * c_exps[i] + 4.0 * g_xx_y_xz_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_z_y[i] = -2.0 * g_0_y_xz_y[i] * c_exps[i] + 4.0 * g_xx_y_xz_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_z_z[i] = -2.0 * g_0_y_xz_z[i] * c_exps[i] + 4.0 * g_xx_y_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_0_z_0_x, g_0_z_0_y, g_0_z_0_z, g_0_z_xx_x, g_0_z_xx_y, g_0_z_xx_z, g_x_0_x_0_x_z_x_x, g_x_0_x_0_x_z_x_y, g_x_0_x_0_x_z_x_z, g_xx_z_0_x, g_xx_z_0_y, g_xx_z_0_z, g_xx_z_xx_x, g_xx_z_xx_y, g_xx_z_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_z_x_x[i] = g_0_z_0_x[i] - 2.0 * g_0_z_xx_x[i] * c_exps[i] - 2.0 * g_xx_z_0_x[i] * a_exp + 4.0 * g_xx_z_xx_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_x_y[i] = g_0_z_0_y[i] - 2.0 * g_0_z_xx_y[i] * c_exps[i] - 2.0 * g_xx_z_0_y[i] * a_exp + 4.0 * g_xx_z_xx_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_x_z[i] = g_0_z_0_z[i] - 2.0 * g_0_z_xx_z[i] * c_exps[i] - 2.0 * g_xx_z_0_z[i] * a_exp + 4.0 * g_xx_z_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_0_z_xy_x, g_0_z_xy_y, g_0_z_xy_z, g_x_0_x_0_x_z_y_x, g_x_0_x_0_x_z_y_y, g_x_0_x_0_x_z_y_z, g_xx_z_xy_x, g_xx_z_xy_y, g_xx_z_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_z_y_x[i] = -2.0 * g_0_z_xy_x[i] * c_exps[i] + 4.0 * g_xx_z_xy_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_y_y[i] = -2.0 * g_0_z_xy_y[i] * c_exps[i] + 4.0 * g_xx_z_xy_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_y_z[i] = -2.0 * g_0_z_xy_z[i] * c_exps[i] + 4.0 * g_xx_z_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_0_z_xz_x, g_0_z_xz_y, g_0_z_xz_z, g_x_0_x_0_x_z_z_x, g_x_0_x_0_x_z_z_y, g_x_0_x_0_x_z_z_z, g_xx_z_xz_x, g_xx_z_xz_y, g_xx_z_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_z_z_x[i] = -2.0 * g_0_z_xz_x[i] * c_exps[i] + 4.0 * g_xx_z_xz_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_z_y[i] = -2.0 * g_0_z_xz_y[i] * c_exps[i] + 4.0 * g_xx_z_xz_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_z_z[i] = -2.0 * g_0_z_xz_z[i] * c_exps[i] + 4.0 * g_xx_z_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (27-30)

    #pragma omp simd aligned(g_x_0_x_0_y_x_x_x, g_x_0_x_0_y_x_x_y, g_x_0_x_0_y_x_x_z, g_xy_x_0_x, g_xy_x_0_y, g_xy_x_0_z, g_xy_x_xx_x, g_xy_x_xx_y, g_xy_x_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_x_x_x[i] = -2.0 * g_xy_x_0_x[i] * a_exp + 4.0 * g_xy_x_xx_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_x_y[i] = -2.0 * g_xy_x_0_y[i] * a_exp + 4.0 * g_xy_x_xx_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_x_z[i] = -2.0 * g_xy_x_0_z[i] * a_exp + 4.0 * g_xy_x_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (30-33)

    #pragma omp simd aligned(g_x_0_x_0_y_x_y_x, g_x_0_x_0_y_x_y_y, g_x_0_x_0_y_x_y_z, g_xy_x_xy_x, g_xy_x_xy_y, g_xy_x_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_x_y_x[i] = 4.0 * g_xy_x_xy_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_y_y[i] = 4.0 * g_xy_x_xy_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_y_z[i] = 4.0 * g_xy_x_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (33-36)

    #pragma omp simd aligned(g_x_0_x_0_y_x_z_x, g_x_0_x_0_y_x_z_y, g_x_0_x_0_y_x_z_z, g_xy_x_xz_x, g_xy_x_xz_y, g_xy_x_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_x_z_x[i] = 4.0 * g_xy_x_xz_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_z_y[i] = 4.0 * g_xy_x_xz_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_z_z[i] = 4.0 * g_xy_x_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (36-39)

    #pragma omp simd aligned(g_x_0_x_0_y_y_x_x, g_x_0_x_0_y_y_x_y, g_x_0_x_0_y_y_x_z, g_xy_y_0_x, g_xy_y_0_y, g_xy_y_0_z, g_xy_y_xx_x, g_xy_y_xx_y, g_xy_y_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_y_x_x[i] = -2.0 * g_xy_y_0_x[i] * a_exp + 4.0 * g_xy_y_xx_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_x_y[i] = -2.0 * g_xy_y_0_y[i] * a_exp + 4.0 * g_xy_y_xx_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_x_z[i] = -2.0 * g_xy_y_0_z[i] * a_exp + 4.0 * g_xy_y_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (39-42)

    #pragma omp simd aligned(g_x_0_x_0_y_y_y_x, g_x_0_x_0_y_y_y_y, g_x_0_x_0_y_y_y_z, g_xy_y_xy_x, g_xy_y_xy_y, g_xy_y_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_y_y_x[i] = 4.0 * g_xy_y_xy_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_y_y[i] = 4.0 * g_xy_y_xy_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_y_z[i] = 4.0 * g_xy_y_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (42-45)

    #pragma omp simd aligned(g_x_0_x_0_y_y_z_x, g_x_0_x_0_y_y_z_y, g_x_0_x_0_y_y_z_z, g_xy_y_xz_x, g_xy_y_xz_y, g_xy_y_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_y_z_x[i] = 4.0 * g_xy_y_xz_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_z_y[i] = 4.0 * g_xy_y_xz_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_z_z[i] = 4.0 * g_xy_y_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (45-48)

    #pragma omp simd aligned(g_x_0_x_0_y_z_x_x, g_x_0_x_0_y_z_x_y, g_x_0_x_0_y_z_x_z, g_xy_z_0_x, g_xy_z_0_y, g_xy_z_0_z, g_xy_z_xx_x, g_xy_z_xx_y, g_xy_z_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_z_x_x[i] = -2.0 * g_xy_z_0_x[i] * a_exp + 4.0 * g_xy_z_xx_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_x_y[i] = -2.0 * g_xy_z_0_y[i] * a_exp + 4.0 * g_xy_z_xx_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_x_z[i] = -2.0 * g_xy_z_0_z[i] * a_exp + 4.0 * g_xy_z_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (48-51)

    #pragma omp simd aligned(g_x_0_x_0_y_z_y_x, g_x_0_x_0_y_z_y_y, g_x_0_x_0_y_z_y_z, g_xy_z_xy_x, g_xy_z_xy_y, g_xy_z_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_z_y_x[i] = 4.0 * g_xy_z_xy_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_y_y[i] = 4.0 * g_xy_z_xy_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_y_z[i] = 4.0 * g_xy_z_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (51-54)

    #pragma omp simd aligned(g_x_0_x_0_y_z_z_x, g_x_0_x_0_y_z_z_y, g_x_0_x_0_y_z_z_z, g_xy_z_xz_x, g_xy_z_xz_y, g_xy_z_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_z_z_x[i] = 4.0 * g_xy_z_xz_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_z_y[i] = 4.0 * g_xy_z_xz_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_z_z[i] = 4.0 * g_xy_z_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (54-57)

    #pragma omp simd aligned(g_x_0_x_0_z_x_x_x, g_x_0_x_0_z_x_x_y, g_x_0_x_0_z_x_x_z, g_xz_x_0_x, g_xz_x_0_y, g_xz_x_0_z, g_xz_x_xx_x, g_xz_x_xx_y, g_xz_x_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_x_x_x[i] = -2.0 * g_xz_x_0_x[i] * a_exp + 4.0 * g_xz_x_xx_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_x_y[i] = -2.0 * g_xz_x_0_y[i] * a_exp + 4.0 * g_xz_x_xx_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_x_z[i] = -2.0 * g_xz_x_0_z[i] * a_exp + 4.0 * g_xz_x_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (57-60)

    #pragma omp simd aligned(g_x_0_x_0_z_x_y_x, g_x_0_x_0_z_x_y_y, g_x_0_x_0_z_x_y_z, g_xz_x_xy_x, g_xz_x_xy_y, g_xz_x_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_x_y_x[i] = 4.0 * g_xz_x_xy_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_y_y[i] = 4.0 * g_xz_x_xy_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_y_z[i] = 4.0 * g_xz_x_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (60-63)

    #pragma omp simd aligned(g_x_0_x_0_z_x_z_x, g_x_0_x_0_z_x_z_y, g_x_0_x_0_z_x_z_z, g_xz_x_xz_x, g_xz_x_xz_y, g_xz_x_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_x_z_x[i] = 4.0 * g_xz_x_xz_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_z_y[i] = 4.0 * g_xz_x_xz_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_z_z[i] = 4.0 * g_xz_x_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (63-66)

    #pragma omp simd aligned(g_x_0_x_0_z_y_x_x, g_x_0_x_0_z_y_x_y, g_x_0_x_0_z_y_x_z, g_xz_y_0_x, g_xz_y_0_y, g_xz_y_0_z, g_xz_y_xx_x, g_xz_y_xx_y, g_xz_y_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_y_x_x[i] = -2.0 * g_xz_y_0_x[i] * a_exp + 4.0 * g_xz_y_xx_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_x_y[i] = -2.0 * g_xz_y_0_y[i] * a_exp + 4.0 * g_xz_y_xx_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_x_z[i] = -2.0 * g_xz_y_0_z[i] * a_exp + 4.0 * g_xz_y_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (66-69)

    #pragma omp simd aligned(g_x_0_x_0_z_y_y_x, g_x_0_x_0_z_y_y_y, g_x_0_x_0_z_y_y_z, g_xz_y_xy_x, g_xz_y_xy_y, g_xz_y_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_y_y_x[i] = 4.0 * g_xz_y_xy_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_y_y[i] = 4.0 * g_xz_y_xy_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_y_z[i] = 4.0 * g_xz_y_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (69-72)

    #pragma omp simd aligned(g_x_0_x_0_z_y_z_x, g_x_0_x_0_z_y_z_y, g_x_0_x_0_z_y_z_z, g_xz_y_xz_x, g_xz_y_xz_y, g_xz_y_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_y_z_x[i] = 4.0 * g_xz_y_xz_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_z_y[i] = 4.0 * g_xz_y_xz_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_z_z[i] = 4.0 * g_xz_y_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (72-75)

    #pragma omp simd aligned(g_x_0_x_0_z_z_x_x, g_x_0_x_0_z_z_x_y, g_x_0_x_0_z_z_x_z, g_xz_z_0_x, g_xz_z_0_y, g_xz_z_0_z, g_xz_z_xx_x, g_xz_z_xx_y, g_xz_z_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_z_x_x[i] = -2.0 * g_xz_z_0_x[i] * a_exp + 4.0 * g_xz_z_xx_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_x_y[i] = -2.0 * g_xz_z_0_y[i] * a_exp + 4.0 * g_xz_z_xx_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_x_z[i] = -2.0 * g_xz_z_0_z[i] * a_exp + 4.0 * g_xz_z_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (75-78)

    #pragma omp simd aligned(g_x_0_x_0_z_z_y_x, g_x_0_x_0_z_z_y_y, g_x_0_x_0_z_z_y_z, g_xz_z_xy_x, g_xz_z_xy_y, g_xz_z_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_z_y_x[i] = 4.0 * g_xz_z_xy_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_y_y[i] = 4.0 * g_xz_z_xy_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_y_z[i] = 4.0 * g_xz_z_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (78-81)

    #pragma omp simd aligned(g_x_0_x_0_z_z_z_x, g_x_0_x_0_z_z_z_y, g_x_0_x_0_z_z_z_z, g_xz_z_xz_x, g_xz_z_xz_y, g_xz_z_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_z_z_x[i] = 4.0 * g_xz_z_xz_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_z_y[i] = 4.0 * g_xz_z_xz_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_z_z[i] = 4.0 * g_xz_z_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (81-84)

    #pragma omp simd aligned(g_0_x_xy_x, g_0_x_xy_y, g_0_x_xy_z, g_x_0_y_0_x_x_x_x, g_x_0_y_0_x_x_x_y, g_x_0_y_0_x_x_x_z, g_xx_x_xy_x, g_xx_x_xy_y, g_xx_x_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_x_x_x[i] = -2.0 * g_0_x_xy_x[i] * c_exps[i] + 4.0 * g_xx_x_xy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_x_y[i] = -2.0 * g_0_x_xy_y[i] * c_exps[i] + 4.0 * g_xx_x_xy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_x_z[i] = -2.0 * g_0_x_xy_z[i] * c_exps[i] + 4.0 * g_xx_x_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (84-87)

    #pragma omp simd aligned(g_0_x_0_x, g_0_x_0_y, g_0_x_0_z, g_0_x_yy_x, g_0_x_yy_y, g_0_x_yy_z, g_x_0_y_0_x_x_y_x, g_x_0_y_0_x_x_y_y, g_x_0_y_0_x_x_y_z, g_xx_x_0_x, g_xx_x_0_y, g_xx_x_0_z, g_xx_x_yy_x, g_xx_x_yy_y, g_xx_x_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_x_y_x[i] = g_0_x_0_x[i] - 2.0 * g_0_x_yy_x[i] * c_exps[i] - 2.0 * g_xx_x_0_x[i] * a_exp + 4.0 * g_xx_x_yy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_y_y[i] = g_0_x_0_y[i] - 2.0 * g_0_x_yy_y[i] * c_exps[i] - 2.0 * g_xx_x_0_y[i] * a_exp + 4.0 * g_xx_x_yy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_y_z[i] = g_0_x_0_z[i] - 2.0 * g_0_x_yy_z[i] * c_exps[i] - 2.0 * g_xx_x_0_z[i] * a_exp + 4.0 * g_xx_x_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (87-90)

    #pragma omp simd aligned(g_0_x_yz_x, g_0_x_yz_y, g_0_x_yz_z, g_x_0_y_0_x_x_z_x, g_x_0_y_0_x_x_z_y, g_x_0_y_0_x_x_z_z, g_xx_x_yz_x, g_xx_x_yz_y, g_xx_x_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_x_z_x[i] = -2.0 * g_0_x_yz_x[i] * c_exps[i] + 4.0 * g_xx_x_yz_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_z_y[i] = -2.0 * g_0_x_yz_y[i] * c_exps[i] + 4.0 * g_xx_x_yz_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_z_z[i] = -2.0 * g_0_x_yz_z[i] * c_exps[i] + 4.0 * g_xx_x_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (90-93)

    #pragma omp simd aligned(g_0_y_xy_x, g_0_y_xy_y, g_0_y_xy_z, g_x_0_y_0_x_y_x_x, g_x_0_y_0_x_y_x_y, g_x_0_y_0_x_y_x_z, g_xx_y_xy_x, g_xx_y_xy_y, g_xx_y_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_y_x_x[i] = -2.0 * g_0_y_xy_x[i] * c_exps[i] + 4.0 * g_xx_y_xy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_x_y[i] = -2.0 * g_0_y_xy_y[i] * c_exps[i] + 4.0 * g_xx_y_xy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_x_z[i] = -2.0 * g_0_y_xy_z[i] * c_exps[i] + 4.0 * g_xx_y_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (93-96)

    #pragma omp simd aligned(g_0_y_0_x, g_0_y_0_y, g_0_y_0_z, g_0_y_yy_x, g_0_y_yy_y, g_0_y_yy_z, g_x_0_y_0_x_y_y_x, g_x_0_y_0_x_y_y_y, g_x_0_y_0_x_y_y_z, g_xx_y_0_x, g_xx_y_0_y, g_xx_y_0_z, g_xx_y_yy_x, g_xx_y_yy_y, g_xx_y_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_y_y_x[i] = g_0_y_0_x[i] - 2.0 * g_0_y_yy_x[i] * c_exps[i] - 2.0 * g_xx_y_0_x[i] * a_exp + 4.0 * g_xx_y_yy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_y_y[i] = g_0_y_0_y[i] - 2.0 * g_0_y_yy_y[i] * c_exps[i] - 2.0 * g_xx_y_0_y[i] * a_exp + 4.0 * g_xx_y_yy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_y_z[i] = g_0_y_0_z[i] - 2.0 * g_0_y_yy_z[i] * c_exps[i] - 2.0 * g_xx_y_0_z[i] * a_exp + 4.0 * g_xx_y_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (96-99)

    #pragma omp simd aligned(g_0_y_yz_x, g_0_y_yz_y, g_0_y_yz_z, g_x_0_y_0_x_y_z_x, g_x_0_y_0_x_y_z_y, g_x_0_y_0_x_y_z_z, g_xx_y_yz_x, g_xx_y_yz_y, g_xx_y_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_y_z_x[i] = -2.0 * g_0_y_yz_x[i] * c_exps[i] + 4.0 * g_xx_y_yz_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_z_y[i] = -2.0 * g_0_y_yz_y[i] * c_exps[i] + 4.0 * g_xx_y_yz_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_z_z[i] = -2.0 * g_0_y_yz_z[i] * c_exps[i] + 4.0 * g_xx_y_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (99-102)

    #pragma omp simd aligned(g_0_z_xy_x, g_0_z_xy_y, g_0_z_xy_z, g_x_0_y_0_x_z_x_x, g_x_0_y_0_x_z_x_y, g_x_0_y_0_x_z_x_z, g_xx_z_xy_x, g_xx_z_xy_y, g_xx_z_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_z_x_x[i] = -2.0 * g_0_z_xy_x[i] * c_exps[i] + 4.0 * g_xx_z_xy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_x_y[i] = -2.0 * g_0_z_xy_y[i] * c_exps[i] + 4.0 * g_xx_z_xy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_x_z[i] = -2.0 * g_0_z_xy_z[i] * c_exps[i] + 4.0 * g_xx_z_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (102-105)

    #pragma omp simd aligned(g_0_z_0_x, g_0_z_0_y, g_0_z_0_z, g_0_z_yy_x, g_0_z_yy_y, g_0_z_yy_z, g_x_0_y_0_x_z_y_x, g_x_0_y_0_x_z_y_y, g_x_0_y_0_x_z_y_z, g_xx_z_0_x, g_xx_z_0_y, g_xx_z_0_z, g_xx_z_yy_x, g_xx_z_yy_y, g_xx_z_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_z_y_x[i] = g_0_z_0_x[i] - 2.0 * g_0_z_yy_x[i] * c_exps[i] - 2.0 * g_xx_z_0_x[i] * a_exp + 4.0 * g_xx_z_yy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_y_y[i] = g_0_z_0_y[i] - 2.0 * g_0_z_yy_y[i] * c_exps[i] - 2.0 * g_xx_z_0_y[i] * a_exp + 4.0 * g_xx_z_yy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_y_z[i] = g_0_z_0_z[i] - 2.0 * g_0_z_yy_z[i] * c_exps[i] - 2.0 * g_xx_z_0_z[i] * a_exp + 4.0 * g_xx_z_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (105-108)

    #pragma omp simd aligned(g_0_z_yz_x, g_0_z_yz_y, g_0_z_yz_z, g_x_0_y_0_x_z_z_x, g_x_0_y_0_x_z_z_y, g_x_0_y_0_x_z_z_z, g_xx_z_yz_x, g_xx_z_yz_y, g_xx_z_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_z_z_x[i] = -2.0 * g_0_z_yz_x[i] * c_exps[i] + 4.0 * g_xx_z_yz_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_z_y[i] = -2.0 * g_0_z_yz_y[i] * c_exps[i] + 4.0 * g_xx_z_yz_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_z_z[i] = -2.0 * g_0_z_yz_z[i] * c_exps[i] + 4.0 * g_xx_z_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (108-111)

    #pragma omp simd aligned(g_x_0_y_0_y_x_x_x, g_x_0_y_0_y_x_x_y, g_x_0_y_0_y_x_x_z, g_xy_x_xy_x, g_xy_x_xy_y, g_xy_x_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_x_x_x[i] = 4.0 * g_xy_x_xy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_x_y[i] = 4.0 * g_xy_x_xy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_x_z[i] = 4.0 * g_xy_x_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (111-114)

    #pragma omp simd aligned(g_x_0_y_0_y_x_y_x, g_x_0_y_0_y_x_y_y, g_x_0_y_0_y_x_y_z, g_xy_x_0_x, g_xy_x_0_y, g_xy_x_0_z, g_xy_x_yy_x, g_xy_x_yy_y, g_xy_x_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_x_y_x[i] = -2.0 * g_xy_x_0_x[i] * a_exp + 4.0 * g_xy_x_yy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_y_y[i] = -2.0 * g_xy_x_0_y[i] * a_exp + 4.0 * g_xy_x_yy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_y_z[i] = -2.0 * g_xy_x_0_z[i] * a_exp + 4.0 * g_xy_x_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (114-117)

    #pragma omp simd aligned(g_x_0_y_0_y_x_z_x, g_x_0_y_0_y_x_z_y, g_x_0_y_0_y_x_z_z, g_xy_x_yz_x, g_xy_x_yz_y, g_xy_x_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_x_z_x[i] = 4.0 * g_xy_x_yz_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_z_y[i] = 4.0 * g_xy_x_yz_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_z_z[i] = 4.0 * g_xy_x_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (117-120)

    #pragma omp simd aligned(g_x_0_y_0_y_y_x_x, g_x_0_y_0_y_y_x_y, g_x_0_y_0_y_y_x_z, g_xy_y_xy_x, g_xy_y_xy_y, g_xy_y_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_y_x_x[i] = 4.0 * g_xy_y_xy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_x_y[i] = 4.0 * g_xy_y_xy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_x_z[i] = 4.0 * g_xy_y_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (120-123)

    #pragma omp simd aligned(g_x_0_y_0_y_y_y_x, g_x_0_y_0_y_y_y_y, g_x_0_y_0_y_y_y_z, g_xy_y_0_x, g_xy_y_0_y, g_xy_y_0_z, g_xy_y_yy_x, g_xy_y_yy_y, g_xy_y_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_y_y_x[i] = -2.0 * g_xy_y_0_x[i] * a_exp + 4.0 * g_xy_y_yy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_y_y[i] = -2.0 * g_xy_y_0_y[i] * a_exp + 4.0 * g_xy_y_yy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_y_z[i] = -2.0 * g_xy_y_0_z[i] * a_exp + 4.0 * g_xy_y_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (123-126)

    #pragma omp simd aligned(g_x_0_y_0_y_y_z_x, g_x_0_y_0_y_y_z_y, g_x_0_y_0_y_y_z_z, g_xy_y_yz_x, g_xy_y_yz_y, g_xy_y_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_y_z_x[i] = 4.0 * g_xy_y_yz_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_z_y[i] = 4.0 * g_xy_y_yz_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_z_z[i] = 4.0 * g_xy_y_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (126-129)

    #pragma omp simd aligned(g_x_0_y_0_y_z_x_x, g_x_0_y_0_y_z_x_y, g_x_0_y_0_y_z_x_z, g_xy_z_xy_x, g_xy_z_xy_y, g_xy_z_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_z_x_x[i] = 4.0 * g_xy_z_xy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_x_y[i] = 4.0 * g_xy_z_xy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_x_z[i] = 4.0 * g_xy_z_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (129-132)

    #pragma omp simd aligned(g_x_0_y_0_y_z_y_x, g_x_0_y_0_y_z_y_y, g_x_0_y_0_y_z_y_z, g_xy_z_0_x, g_xy_z_0_y, g_xy_z_0_z, g_xy_z_yy_x, g_xy_z_yy_y, g_xy_z_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_z_y_x[i] = -2.0 * g_xy_z_0_x[i] * a_exp + 4.0 * g_xy_z_yy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_y_y[i] = -2.0 * g_xy_z_0_y[i] * a_exp + 4.0 * g_xy_z_yy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_y_z[i] = -2.0 * g_xy_z_0_z[i] * a_exp + 4.0 * g_xy_z_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (132-135)

    #pragma omp simd aligned(g_x_0_y_0_y_z_z_x, g_x_0_y_0_y_z_z_y, g_x_0_y_0_y_z_z_z, g_xy_z_yz_x, g_xy_z_yz_y, g_xy_z_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_z_z_x[i] = 4.0 * g_xy_z_yz_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_z_y[i] = 4.0 * g_xy_z_yz_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_z_z[i] = 4.0 * g_xy_z_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (135-138)

    #pragma omp simd aligned(g_x_0_y_0_z_x_x_x, g_x_0_y_0_z_x_x_y, g_x_0_y_0_z_x_x_z, g_xz_x_xy_x, g_xz_x_xy_y, g_xz_x_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_x_x_x[i] = 4.0 * g_xz_x_xy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_x_y[i] = 4.0 * g_xz_x_xy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_x_z[i] = 4.0 * g_xz_x_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (138-141)

    #pragma omp simd aligned(g_x_0_y_0_z_x_y_x, g_x_0_y_0_z_x_y_y, g_x_0_y_0_z_x_y_z, g_xz_x_0_x, g_xz_x_0_y, g_xz_x_0_z, g_xz_x_yy_x, g_xz_x_yy_y, g_xz_x_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_x_y_x[i] = -2.0 * g_xz_x_0_x[i] * a_exp + 4.0 * g_xz_x_yy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_y_y[i] = -2.0 * g_xz_x_0_y[i] * a_exp + 4.0 * g_xz_x_yy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_y_z[i] = -2.0 * g_xz_x_0_z[i] * a_exp + 4.0 * g_xz_x_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (141-144)

    #pragma omp simd aligned(g_x_0_y_0_z_x_z_x, g_x_0_y_0_z_x_z_y, g_x_0_y_0_z_x_z_z, g_xz_x_yz_x, g_xz_x_yz_y, g_xz_x_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_x_z_x[i] = 4.0 * g_xz_x_yz_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_z_y[i] = 4.0 * g_xz_x_yz_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_z_z[i] = 4.0 * g_xz_x_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (144-147)

    #pragma omp simd aligned(g_x_0_y_0_z_y_x_x, g_x_0_y_0_z_y_x_y, g_x_0_y_0_z_y_x_z, g_xz_y_xy_x, g_xz_y_xy_y, g_xz_y_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_y_x_x[i] = 4.0 * g_xz_y_xy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_x_y[i] = 4.0 * g_xz_y_xy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_x_z[i] = 4.0 * g_xz_y_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (147-150)

    #pragma omp simd aligned(g_x_0_y_0_z_y_y_x, g_x_0_y_0_z_y_y_y, g_x_0_y_0_z_y_y_z, g_xz_y_0_x, g_xz_y_0_y, g_xz_y_0_z, g_xz_y_yy_x, g_xz_y_yy_y, g_xz_y_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_y_y_x[i] = -2.0 * g_xz_y_0_x[i] * a_exp + 4.0 * g_xz_y_yy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_y_y[i] = -2.0 * g_xz_y_0_y[i] * a_exp + 4.0 * g_xz_y_yy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_y_z[i] = -2.0 * g_xz_y_0_z[i] * a_exp + 4.0 * g_xz_y_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (150-153)

    #pragma omp simd aligned(g_x_0_y_0_z_y_z_x, g_x_0_y_0_z_y_z_y, g_x_0_y_0_z_y_z_z, g_xz_y_yz_x, g_xz_y_yz_y, g_xz_y_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_y_z_x[i] = 4.0 * g_xz_y_yz_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_z_y[i] = 4.0 * g_xz_y_yz_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_z_z[i] = 4.0 * g_xz_y_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (153-156)

    #pragma omp simd aligned(g_x_0_y_0_z_z_x_x, g_x_0_y_0_z_z_x_y, g_x_0_y_0_z_z_x_z, g_xz_z_xy_x, g_xz_z_xy_y, g_xz_z_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_z_x_x[i] = 4.0 * g_xz_z_xy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_x_y[i] = 4.0 * g_xz_z_xy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_x_z[i] = 4.0 * g_xz_z_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (156-159)

    #pragma omp simd aligned(g_x_0_y_0_z_z_y_x, g_x_0_y_0_z_z_y_y, g_x_0_y_0_z_z_y_z, g_xz_z_0_x, g_xz_z_0_y, g_xz_z_0_z, g_xz_z_yy_x, g_xz_z_yy_y, g_xz_z_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_z_y_x[i] = -2.0 * g_xz_z_0_x[i] * a_exp + 4.0 * g_xz_z_yy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_y_y[i] = -2.0 * g_xz_z_0_y[i] * a_exp + 4.0 * g_xz_z_yy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_y_z[i] = -2.0 * g_xz_z_0_z[i] * a_exp + 4.0 * g_xz_z_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (159-162)

    #pragma omp simd aligned(g_x_0_y_0_z_z_z_x, g_x_0_y_0_z_z_z_y, g_x_0_y_0_z_z_z_z, g_xz_z_yz_x, g_xz_z_yz_y, g_xz_z_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_z_z_x[i] = 4.0 * g_xz_z_yz_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_z_y[i] = 4.0 * g_xz_z_yz_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_z_z[i] = 4.0 * g_xz_z_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (162-165)

    #pragma omp simd aligned(g_0_x_xz_x, g_0_x_xz_y, g_0_x_xz_z, g_x_0_z_0_x_x_x_x, g_x_0_z_0_x_x_x_y, g_x_0_z_0_x_x_x_z, g_xx_x_xz_x, g_xx_x_xz_y, g_xx_x_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_x_x_x[i] = -2.0 * g_0_x_xz_x[i] * c_exps[i] + 4.0 * g_xx_x_xz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_x_y[i] = -2.0 * g_0_x_xz_y[i] * c_exps[i] + 4.0 * g_xx_x_xz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_x_z[i] = -2.0 * g_0_x_xz_z[i] * c_exps[i] + 4.0 * g_xx_x_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (165-168)

    #pragma omp simd aligned(g_0_x_yz_x, g_0_x_yz_y, g_0_x_yz_z, g_x_0_z_0_x_x_y_x, g_x_0_z_0_x_x_y_y, g_x_0_z_0_x_x_y_z, g_xx_x_yz_x, g_xx_x_yz_y, g_xx_x_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_x_y_x[i] = -2.0 * g_0_x_yz_x[i] * c_exps[i] + 4.0 * g_xx_x_yz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_y_y[i] = -2.0 * g_0_x_yz_y[i] * c_exps[i] + 4.0 * g_xx_x_yz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_y_z[i] = -2.0 * g_0_x_yz_z[i] * c_exps[i] + 4.0 * g_xx_x_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (168-171)

    #pragma omp simd aligned(g_0_x_0_x, g_0_x_0_y, g_0_x_0_z, g_0_x_zz_x, g_0_x_zz_y, g_0_x_zz_z, g_x_0_z_0_x_x_z_x, g_x_0_z_0_x_x_z_y, g_x_0_z_0_x_x_z_z, g_xx_x_0_x, g_xx_x_0_y, g_xx_x_0_z, g_xx_x_zz_x, g_xx_x_zz_y, g_xx_x_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_x_z_x[i] = g_0_x_0_x[i] - 2.0 * g_0_x_zz_x[i] * c_exps[i] - 2.0 * g_xx_x_0_x[i] * a_exp + 4.0 * g_xx_x_zz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_z_y[i] = g_0_x_0_y[i] - 2.0 * g_0_x_zz_y[i] * c_exps[i] - 2.0 * g_xx_x_0_y[i] * a_exp + 4.0 * g_xx_x_zz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_z_z[i] = g_0_x_0_z[i] - 2.0 * g_0_x_zz_z[i] * c_exps[i] - 2.0 * g_xx_x_0_z[i] * a_exp + 4.0 * g_xx_x_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (171-174)

    #pragma omp simd aligned(g_0_y_xz_x, g_0_y_xz_y, g_0_y_xz_z, g_x_0_z_0_x_y_x_x, g_x_0_z_0_x_y_x_y, g_x_0_z_0_x_y_x_z, g_xx_y_xz_x, g_xx_y_xz_y, g_xx_y_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_y_x_x[i] = -2.0 * g_0_y_xz_x[i] * c_exps[i] + 4.0 * g_xx_y_xz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_x_y[i] = -2.0 * g_0_y_xz_y[i] * c_exps[i] + 4.0 * g_xx_y_xz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_x_z[i] = -2.0 * g_0_y_xz_z[i] * c_exps[i] + 4.0 * g_xx_y_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (174-177)

    #pragma omp simd aligned(g_0_y_yz_x, g_0_y_yz_y, g_0_y_yz_z, g_x_0_z_0_x_y_y_x, g_x_0_z_0_x_y_y_y, g_x_0_z_0_x_y_y_z, g_xx_y_yz_x, g_xx_y_yz_y, g_xx_y_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_y_y_x[i] = -2.0 * g_0_y_yz_x[i] * c_exps[i] + 4.0 * g_xx_y_yz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_y_y[i] = -2.0 * g_0_y_yz_y[i] * c_exps[i] + 4.0 * g_xx_y_yz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_y_z[i] = -2.0 * g_0_y_yz_z[i] * c_exps[i] + 4.0 * g_xx_y_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (177-180)

    #pragma omp simd aligned(g_0_y_0_x, g_0_y_0_y, g_0_y_0_z, g_0_y_zz_x, g_0_y_zz_y, g_0_y_zz_z, g_x_0_z_0_x_y_z_x, g_x_0_z_0_x_y_z_y, g_x_0_z_0_x_y_z_z, g_xx_y_0_x, g_xx_y_0_y, g_xx_y_0_z, g_xx_y_zz_x, g_xx_y_zz_y, g_xx_y_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_y_z_x[i] = g_0_y_0_x[i] - 2.0 * g_0_y_zz_x[i] * c_exps[i] - 2.0 * g_xx_y_0_x[i] * a_exp + 4.0 * g_xx_y_zz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_z_y[i] = g_0_y_0_y[i] - 2.0 * g_0_y_zz_y[i] * c_exps[i] - 2.0 * g_xx_y_0_y[i] * a_exp + 4.0 * g_xx_y_zz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_z_z[i] = g_0_y_0_z[i] - 2.0 * g_0_y_zz_z[i] * c_exps[i] - 2.0 * g_xx_y_0_z[i] * a_exp + 4.0 * g_xx_y_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (180-183)

    #pragma omp simd aligned(g_0_z_xz_x, g_0_z_xz_y, g_0_z_xz_z, g_x_0_z_0_x_z_x_x, g_x_0_z_0_x_z_x_y, g_x_0_z_0_x_z_x_z, g_xx_z_xz_x, g_xx_z_xz_y, g_xx_z_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_z_x_x[i] = -2.0 * g_0_z_xz_x[i] * c_exps[i] + 4.0 * g_xx_z_xz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_x_y[i] = -2.0 * g_0_z_xz_y[i] * c_exps[i] + 4.0 * g_xx_z_xz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_x_z[i] = -2.0 * g_0_z_xz_z[i] * c_exps[i] + 4.0 * g_xx_z_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (183-186)

    #pragma omp simd aligned(g_0_z_yz_x, g_0_z_yz_y, g_0_z_yz_z, g_x_0_z_0_x_z_y_x, g_x_0_z_0_x_z_y_y, g_x_0_z_0_x_z_y_z, g_xx_z_yz_x, g_xx_z_yz_y, g_xx_z_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_z_y_x[i] = -2.0 * g_0_z_yz_x[i] * c_exps[i] + 4.0 * g_xx_z_yz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_y_y[i] = -2.0 * g_0_z_yz_y[i] * c_exps[i] + 4.0 * g_xx_z_yz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_y_z[i] = -2.0 * g_0_z_yz_z[i] * c_exps[i] + 4.0 * g_xx_z_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (186-189)

    #pragma omp simd aligned(g_0_z_0_x, g_0_z_0_y, g_0_z_0_z, g_0_z_zz_x, g_0_z_zz_y, g_0_z_zz_z, g_x_0_z_0_x_z_z_x, g_x_0_z_0_x_z_z_y, g_x_0_z_0_x_z_z_z, g_xx_z_0_x, g_xx_z_0_y, g_xx_z_0_z, g_xx_z_zz_x, g_xx_z_zz_y, g_xx_z_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_z_z_x[i] = g_0_z_0_x[i] - 2.0 * g_0_z_zz_x[i] * c_exps[i] - 2.0 * g_xx_z_0_x[i] * a_exp + 4.0 * g_xx_z_zz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_z_y[i] = g_0_z_0_y[i] - 2.0 * g_0_z_zz_y[i] * c_exps[i] - 2.0 * g_xx_z_0_y[i] * a_exp + 4.0 * g_xx_z_zz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_z_z[i] = g_0_z_0_z[i] - 2.0 * g_0_z_zz_z[i] * c_exps[i] - 2.0 * g_xx_z_0_z[i] * a_exp + 4.0 * g_xx_z_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (189-192)

    #pragma omp simd aligned(g_x_0_z_0_y_x_x_x, g_x_0_z_0_y_x_x_y, g_x_0_z_0_y_x_x_z, g_xy_x_xz_x, g_xy_x_xz_y, g_xy_x_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_x_x_x[i] = 4.0 * g_xy_x_xz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_x_y[i] = 4.0 * g_xy_x_xz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_x_z[i] = 4.0 * g_xy_x_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (192-195)

    #pragma omp simd aligned(g_x_0_z_0_y_x_y_x, g_x_0_z_0_y_x_y_y, g_x_0_z_0_y_x_y_z, g_xy_x_yz_x, g_xy_x_yz_y, g_xy_x_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_x_y_x[i] = 4.0 * g_xy_x_yz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_y_y[i] = 4.0 * g_xy_x_yz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_y_z[i] = 4.0 * g_xy_x_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (195-198)

    #pragma omp simd aligned(g_x_0_z_0_y_x_z_x, g_x_0_z_0_y_x_z_y, g_x_0_z_0_y_x_z_z, g_xy_x_0_x, g_xy_x_0_y, g_xy_x_0_z, g_xy_x_zz_x, g_xy_x_zz_y, g_xy_x_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_x_z_x[i] = -2.0 * g_xy_x_0_x[i] * a_exp + 4.0 * g_xy_x_zz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_z_y[i] = -2.0 * g_xy_x_0_y[i] * a_exp + 4.0 * g_xy_x_zz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_z_z[i] = -2.0 * g_xy_x_0_z[i] * a_exp + 4.0 * g_xy_x_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (198-201)

    #pragma omp simd aligned(g_x_0_z_0_y_y_x_x, g_x_0_z_0_y_y_x_y, g_x_0_z_0_y_y_x_z, g_xy_y_xz_x, g_xy_y_xz_y, g_xy_y_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_y_x_x[i] = 4.0 * g_xy_y_xz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_x_y[i] = 4.0 * g_xy_y_xz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_x_z[i] = 4.0 * g_xy_y_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (201-204)

    #pragma omp simd aligned(g_x_0_z_0_y_y_y_x, g_x_0_z_0_y_y_y_y, g_x_0_z_0_y_y_y_z, g_xy_y_yz_x, g_xy_y_yz_y, g_xy_y_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_y_y_x[i] = 4.0 * g_xy_y_yz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_y_y[i] = 4.0 * g_xy_y_yz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_y_z[i] = 4.0 * g_xy_y_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (204-207)

    #pragma omp simd aligned(g_x_0_z_0_y_y_z_x, g_x_0_z_0_y_y_z_y, g_x_0_z_0_y_y_z_z, g_xy_y_0_x, g_xy_y_0_y, g_xy_y_0_z, g_xy_y_zz_x, g_xy_y_zz_y, g_xy_y_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_y_z_x[i] = -2.0 * g_xy_y_0_x[i] * a_exp + 4.0 * g_xy_y_zz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_z_y[i] = -2.0 * g_xy_y_0_y[i] * a_exp + 4.0 * g_xy_y_zz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_z_z[i] = -2.0 * g_xy_y_0_z[i] * a_exp + 4.0 * g_xy_y_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (207-210)

    #pragma omp simd aligned(g_x_0_z_0_y_z_x_x, g_x_0_z_0_y_z_x_y, g_x_0_z_0_y_z_x_z, g_xy_z_xz_x, g_xy_z_xz_y, g_xy_z_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_z_x_x[i] = 4.0 * g_xy_z_xz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_x_y[i] = 4.0 * g_xy_z_xz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_x_z[i] = 4.0 * g_xy_z_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (210-213)

    #pragma omp simd aligned(g_x_0_z_0_y_z_y_x, g_x_0_z_0_y_z_y_y, g_x_0_z_0_y_z_y_z, g_xy_z_yz_x, g_xy_z_yz_y, g_xy_z_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_z_y_x[i] = 4.0 * g_xy_z_yz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_y_y[i] = 4.0 * g_xy_z_yz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_y_z[i] = 4.0 * g_xy_z_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (213-216)

    #pragma omp simd aligned(g_x_0_z_0_y_z_z_x, g_x_0_z_0_y_z_z_y, g_x_0_z_0_y_z_z_z, g_xy_z_0_x, g_xy_z_0_y, g_xy_z_0_z, g_xy_z_zz_x, g_xy_z_zz_y, g_xy_z_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_z_z_x[i] = -2.0 * g_xy_z_0_x[i] * a_exp + 4.0 * g_xy_z_zz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_z_y[i] = -2.0 * g_xy_z_0_y[i] * a_exp + 4.0 * g_xy_z_zz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_z_z[i] = -2.0 * g_xy_z_0_z[i] * a_exp + 4.0 * g_xy_z_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (216-219)

    #pragma omp simd aligned(g_x_0_z_0_z_x_x_x, g_x_0_z_0_z_x_x_y, g_x_0_z_0_z_x_x_z, g_xz_x_xz_x, g_xz_x_xz_y, g_xz_x_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_x_x_x[i] = 4.0 * g_xz_x_xz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_x_y[i] = 4.0 * g_xz_x_xz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_x_z[i] = 4.0 * g_xz_x_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (219-222)

    #pragma omp simd aligned(g_x_0_z_0_z_x_y_x, g_x_0_z_0_z_x_y_y, g_x_0_z_0_z_x_y_z, g_xz_x_yz_x, g_xz_x_yz_y, g_xz_x_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_x_y_x[i] = 4.0 * g_xz_x_yz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_y_y[i] = 4.0 * g_xz_x_yz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_y_z[i] = 4.0 * g_xz_x_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (222-225)

    #pragma omp simd aligned(g_x_0_z_0_z_x_z_x, g_x_0_z_0_z_x_z_y, g_x_0_z_0_z_x_z_z, g_xz_x_0_x, g_xz_x_0_y, g_xz_x_0_z, g_xz_x_zz_x, g_xz_x_zz_y, g_xz_x_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_x_z_x[i] = -2.0 * g_xz_x_0_x[i] * a_exp + 4.0 * g_xz_x_zz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_z_y[i] = -2.0 * g_xz_x_0_y[i] * a_exp + 4.0 * g_xz_x_zz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_z_z[i] = -2.0 * g_xz_x_0_z[i] * a_exp + 4.0 * g_xz_x_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (225-228)

    #pragma omp simd aligned(g_x_0_z_0_z_y_x_x, g_x_0_z_0_z_y_x_y, g_x_0_z_0_z_y_x_z, g_xz_y_xz_x, g_xz_y_xz_y, g_xz_y_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_y_x_x[i] = 4.0 * g_xz_y_xz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_x_y[i] = 4.0 * g_xz_y_xz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_x_z[i] = 4.0 * g_xz_y_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (228-231)

    #pragma omp simd aligned(g_x_0_z_0_z_y_y_x, g_x_0_z_0_z_y_y_y, g_x_0_z_0_z_y_y_z, g_xz_y_yz_x, g_xz_y_yz_y, g_xz_y_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_y_y_x[i] = 4.0 * g_xz_y_yz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_y_y[i] = 4.0 * g_xz_y_yz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_y_z[i] = 4.0 * g_xz_y_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (231-234)

    #pragma omp simd aligned(g_x_0_z_0_z_y_z_x, g_x_0_z_0_z_y_z_y, g_x_0_z_0_z_y_z_z, g_xz_y_0_x, g_xz_y_0_y, g_xz_y_0_z, g_xz_y_zz_x, g_xz_y_zz_y, g_xz_y_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_y_z_x[i] = -2.0 * g_xz_y_0_x[i] * a_exp + 4.0 * g_xz_y_zz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_z_y[i] = -2.0 * g_xz_y_0_y[i] * a_exp + 4.0 * g_xz_y_zz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_z_z[i] = -2.0 * g_xz_y_0_z[i] * a_exp + 4.0 * g_xz_y_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (234-237)

    #pragma omp simd aligned(g_x_0_z_0_z_z_x_x, g_x_0_z_0_z_z_x_y, g_x_0_z_0_z_z_x_z, g_xz_z_xz_x, g_xz_z_xz_y, g_xz_z_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_z_x_x[i] = 4.0 * g_xz_z_xz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_x_y[i] = 4.0 * g_xz_z_xz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_x_z[i] = 4.0 * g_xz_z_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (237-240)

    #pragma omp simd aligned(g_x_0_z_0_z_z_y_x, g_x_0_z_0_z_z_y_y, g_x_0_z_0_z_z_y_z, g_xz_z_yz_x, g_xz_z_yz_y, g_xz_z_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_z_y_x[i] = 4.0 * g_xz_z_yz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_y_y[i] = 4.0 * g_xz_z_yz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_y_z[i] = 4.0 * g_xz_z_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (240-243)

    #pragma omp simd aligned(g_x_0_z_0_z_z_z_x, g_x_0_z_0_z_z_z_y, g_x_0_z_0_z_z_z_z, g_xz_z_0_x, g_xz_z_0_y, g_xz_z_0_z, g_xz_z_zz_x, g_xz_z_zz_y, g_xz_z_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_z_z_x[i] = -2.0 * g_xz_z_0_x[i] * a_exp + 4.0 * g_xz_z_zz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_z_y[i] = -2.0 * g_xz_z_0_y[i] * a_exp + 4.0 * g_xz_z_zz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_z_z[i] = -2.0 * g_xz_z_0_z[i] * a_exp + 4.0 * g_xz_z_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (243-246)

    #pragma omp simd aligned(g_xy_x_0_x, g_xy_x_0_y, g_xy_x_0_z, g_xy_x_xx_x, g_xy_x_xx_y, g_xy_x_xx_z, g_y_0_x_0_x_x_x_x, g_y_0_x_0_x_x_x_y, g_y_0_x_0_x_x_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_x_x_x[i] = -2.0 * g_xy_x_0_x[i] * a_exp + 4.0 * g_xy_x_xx_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_x_y[i] = -2.0 * g_xy_x_0_y[i] * a_exp + 4.0 * g_xy_x_xx_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_x_z[i] = -2.0 * g_xy_x_0_z[i] * a_exp + 4.0 * g_xy_x_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (246-249)

    #pragma omp simd aligned(g_xy_x_xy_x, g_xy_x_xy_y, g_xy_x_xy_z, g_y_0_x_0_x_x_y_x, g_y_0_x_0_x_x_y_y, g_y_0_x_0_x_x_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_x_y_x[i] = 4.0 * g_xy_x_xy_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_y_y[i] = 4.0 * g_xy_x_xy_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_y_z[i] = 4.0 * g_xy_x_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (249-252)

    #pragma omp simd aligned(g_xy_x_xz_x, g_xy_x_xz_y, g_xy_x_xz_z, g_y_0_x_0_x_x_z_x, g_y_0_x_0_x_x_z_y, g_y_0_x_0_x_x_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_x_z_x[i] = 4.0 * g_xy_x_xz_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_z_y[i] = 4.0 * g_xy_x_xz_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_z_z[i] = 4.0 * g_xy_x_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (252-255)

    #pragma omp simd aligned(g_xy_y_0_x, g_xy_y_0_y, g_xy_y_0_z, g_xy_y_xx_x, g_xy_y_xx_y, g_xy_y_xx_z, g_y_0_x_0_x_y_x_x, g_y_0_x_0_x_y_x_y, g_y_0_x_0_x_y_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_y_x_x[i] = -2.0 * g_xy_y_0_x[i] * a_exp + 4.0 * g_xy_y_xx_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_x_y[i] = -2.0 * g_xy_y_0_y[i] * a_exp + 4.0 * g_xy_y_xx_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_x_z[i] = -2.0 * g_xy_y_0_z[i] * a_exp + 4.0 * g_xy_y_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (255-258)

    #pragma omp simd aligned(g_xy_y_xy_x, g_xy_y_xy_y, g_xy_y_xy_z, g_y_0_x_0_x_y_y_x, g_y_0_x_0_x_y_y_y, g_y_0_x_0_x_y_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_y_y_x[i] = 4.0 * g_xy_y_xy_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_y_y[i] = 4.0 * g_xy_y_xy_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_y_z[i] = 4.0 * g_xy_y_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (258-261)

    #pragma omp simd aligned(g_xy_y_xz_x, g_xy_y_xz_y, g_xy_y_xz_z, g_y_0_x_0_x_y_z_x, g_y_0_x_0_x_y_z_y, g_y_0_x_0_x_y_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_y_z_x[i] = 4.0 * g_xy_y_xz_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_z_y[i] = 4.0 * g_xy_y_xz_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_z_z[i] = 4.0 * g_xy_y_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (261-264)

    #pragma omp simd aligned(g_xy_z_0_x, g_xy_z_0_y, g_xy_z_0_z, g_xy_z_xx_x, g_xy_z_xx_y, g_xy_z_xx_z, g_y_0_x_0_x_z_x_x, g_y_0_x_0_x_z_x_y, g_y_0_x_0_x_z_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_z_x_x[i] = -2.0 * g_xy_z_0_x[i] * a_exp + 4.0 * g_xy_z_xx_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_x_y[i] = -2.0 * g_xy_z_0_y[i] * a_exp + 4.0 * g_xy_z_xx_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_x_z[i] = -2.0 * g_xy_z_0_z[i] * a_exp + 4.0 * g_xy_z_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (264-267)

    #pragma omp simd aligned(g_xy_z_xy_x, g_xy_z_xy_y, g_xy_z_xy_z, g_y_0_x_0_x_z_y_x, g_y_0_x_0_x_z_y_y, g_y_0_x_0_x_z_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_z_y_x[i] = 4.0 * g_xy_z_xy_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_y_y[i] = 4.0 * g_xy_z_xy_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_y_z[i] = 4.0 * g_xy_z_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (267-270)

    #pragma omp simd aligned(g_xy_z_xz_x, g_xy_z_xz_y, g_xy_z_xz_z, g_y_0_x_0_x_z_z_x, g_y_0_x_0_x_z_z_y, g_y_0_x_0_x_z_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_z_z_x[i] = 4.0 * g_xy_z_xz_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_z_y[i] = 4.0 * g_xy_z_xz_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_z_z[i] = 4.0 * g_xy_z_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (270-273)

    #pragma omp simd aligned(g_0_x_0_x, g_0_x_0_y, g_0_x_0_z, g_0_x_xx_x, g_0_x_xx_y, g_0_x_xx_z, g_y_0_x_0_y_x_x_x, g_y_0_x_0_y_x_x_y, g_y_0_x_0_y_x_x_z, g_yy_x_0_x, g_yy_x_0_y, g_yy_x_0_z, g_yy_x_xx_x, g_yy_x_xx_y, g_yy_x_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_x_x_x[i] = g_0_x_0_x[i] - 2.0 * g_0_x_xx_x[i] * c_exps[i] - 2.0 * g_yy_x_0_x[i] * a_exp + 4.0 * g_yy_x_xx_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_x_y[i] = g_0_x_0_y[i] - 2.0 * g_0_x_xx_y[i] * c_exps[i] - 2.0 * g_yy_x_0_y[i] * a_exp + 4.0 * g_yy_x_xx_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_x_z[i] = g_0_x_0_z[i] - 2.0 * g_0_x_xx_z[i] * c_exps[i] - 2.0 * g_yy_x_0_z[i] * a_exp + 4.0 * g_yy_x_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (273-276)

    #pragma omp simd aligned(g_0_x_xy_x, g_0_x_xy_y, g_0_x_xy_z, g_y_0_x_0_y_x_y_x, g_y_0_x_0_y_x_y_y, g_y_0_x_0_y_x_y_z, g_yy_x_xy_x, g_yy_x_xy_y, g_yy_x_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_x_y_x[i] = -2.0 * g_0_x_xy_x[i] * c_exps[i] + 4.0 * g_yy_x_xy_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_y_y[i] = -2.0 * g_0_x_xy_y[i] * c_exps[i] + 4.0 * g_yy_x_xy_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_y_z[i] = -2.0 * g_0_x_xy_z[i] * c_exps[i] + 4.0 * g_yy_x_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (276-279)

    #pragma omp simd aligned(g_0_x_xz_x, g_0_x_xz_y, g_0_x_xz_z, g_y_0_x_0_y_x_z_x, g_y_0_x_0_y_x_z_y, g_y_0_x_0_y_x_z_z, g_yy_x_xz_x, g_yy_x_xz_y, g_yy_x_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_x_z_x[i] = -2.0 * g_0_x_xz_x[i] * c_exps[i] + 4.0 * g_yy_x_xz_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_z_y[i] = -2.0 * g_0_x_xz_y[i] * c_exps[i] + 4.0 * g_yy_x_xz_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_z_z[i] = -2.0 * g_0_x_xz_z[i] * c_exps[i] + 4.0 * g_yy_x_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (279-282)

    #pragma omp simd aligned(g_0_y_0_x, g_0_y_0_y, g_0_y_0_z, g_0_y_xx_x, g_0_y_xx_y, g_0_y_xx_z, g_y_0_x_0_y_y_x_x, g_y_0_x_0_y_y_x_y, g_y_0_x_0_y_y_x_z, g_yy_y_0_x, g_yy_y_0_y, g_yy_y_0_z, g_yy_y_xx_x, g_yy_y_xx_y, g_yy_y_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_y_x_x[i] = g_0_y_0_x[i] - 2.0 * g_0_y_xx_x[i] * c_exps[i] - 2.0 * g_yy_y_0_x[i] * a_exp + 4.0 * g_yy_y_xx_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_x_y[i] = g_0_y_0_y[i] - 2.0 * g_0_y_xx_y[i] * c_exps[i] - 2.0 * g_yy_y_0_y[i] * a_exp + 4.0 * g_yy_y_xx_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_x_z[i] = g_0_y_0_z[i] - 2.0 * g_0_y_xx_z[i] * c_exps[i] - 2.0 * g_yy_y_0_z[i] * a_exp + 4.0 * g_yy_y_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (282-285)

    #pragma omp simd aligned(g_0_y_xy_x, g_0_y_xy_y, g_0_y_xy_z, g_y_0_x_0_y_y_y_x, g_y_0_x_0_y_y_y_y, g_y_0_x_0_y_y_y_z, g_yy_y_xy_x, g_yy_y_xy_y, g_yy_y_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_y_y_x[i] = -2.0 * g_0_y_xy_x[i] * c_exps[i] + 4.0 * g_yy_y_xy_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_y_y[i] = -2.0 * g_0_y_xy_y[i] * c_exps[i] + 4.0 * g_yy_y_xy_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_y_z[i] = -2.0 * g_0_y_xy_z[i] * c_exps[i] + 4.0 * g_yy_y_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (285-288)

    #pragma omp simd aligned(g_0_y_xz_x, g_0_y_xz_y, g_0_y_xz_z, g_y_0_x_0_y_y_z_x, g_y_0_x_0_y_y_z_y, g_y_0_x_0_y_y_z_z, g_yy_y_xz_x, g_yy_y_xz_y, g_yy_y_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_y_z_x[i] = -2.0 * g_0_y_xz_x[i] * c_exps[i] + 4.0 * g_yy_y_xz_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_z_y[i] = -2.0 * g_0_y_xz_y[i] * c_exps[i] + 4.0 * g_yy_y_xz_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_z_z[i] = -2.0 * g_0_y_xz_z[i] * c_exps[i] + 4.0 * g_yy_y_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (288-291)

    #pragma omp simd aligned(g_0_z_0_x, g_0_z_0_y, g_0_z_0_z, g_0_z_xx_x, g_0_z_xx_y, g_0_z_xx_z, g_y_0_x_0_y_z_x_x, g_y_0_x_0_y_z_x_y, g_y_0_x_0_y_z_x_z, g_yy_z_0_x, g_yy_z_0_y, g_yy_z_0_z, g_yy_z_xx_x, g_yy_z_xx_y, g_yy_z_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_z_x_x[i] = g_0_z_0_x[i] - 2.0 * g_0_z_xx_x[i] * c_exps[i] - 2.0 * g_yy_z_0_x[i] * a_exp + 4.0 * g_yy_z_xx_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_x_y[i] = g_0_z_0_y[i] - 2.0 * g_0_z_xx_y[i] * c_exps[i] - 2.0 * g_yy_z_0_y[i] * a_exp + 4.0 * g_yy_z_xx_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_x_z[i] = g_0_z_0_z[i] - 2.0 * g_0_z_xx_z[i] * c_exps[i] - 2.0 * g_yy_z_0_z[i] * a_exp + 4.0 * g_yy_z_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (291-294)

    #pragma omp simd aligned(g_0_z_xy_x, g_0_z_xy_y, g_0_z_xy_z, g_y_0_x_0_y_z_y_x, g_y_0_x_0_y_z_y_y, g_y_0_x_0_y_z_y_z, g_yy_z_xy_x, g_yy_z_xy_y, g_yy_z_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_z_y_x[i] = -2.0 * g_0_z_xy_x[i] * c_exps[i] + 4.0 * g_yy_z_xy_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_y_y[i] = -2.0 * g_0_z_xy_y[i] * c_exps[i] + 4.0 * g_yy_z_xy_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_y_z[i] = -2.0 * g_0_z_xy_z[i] * c_exps[i] + 4.0 * g_yy_z_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (294-297)

    #pragma omp simd aligned(g_0_z_xz_x, g_0_z_xz_y, g_0_z_xz_z, g_y_0_x_0_y_z_z_x, g_y_0_x_0_y_z_z_y, g_y_0_x_0_y_z_z_z, g_yy_z_xz_x, g_yy_z_xz_y, g_yy_z_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_z_z_x[i] = -2.0 * g_0_z_xz_x[i] * c_exps[i] + 4.0 * g_yy_z_xz_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_z_y[i] = -2.0 * g_0_z_xz_y[i] * c_exps[i] + 4.0 * g_yy_z_xz_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_z_z[i] = -2.0 * g_0_z_xz_z[i] * c_exps[i] + 4.0 * g_yy_z_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (297-300)

    #pragma omp simd aligned(g_y_0_x_0_z_x_x_x, g_y_0_x_0_z_x_x_y, g_y_0_x_0_z_x_x_z, g_yz_x_0_x, g_yz_x_0_y, g_yz_x_0_z, g_yz_x_xx_x, g_yz_x_xx_y, g_yz_x_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_x_x_x[i] = -2.0 * g_yz_x_0_x[i] * a_exp + 4.0 * g_yz_x_xx_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_x_y[i] = -2.0 * g_yz_x_0_y[i] * a_exp + 4.0 * g_yz_x_xx_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_x_z[i] = -2.0 * g_yz_x_0_z[i] * a_exp + 4.0 * g_yz_x_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (300-303)

    #pragma omp simd aligned(g_y_0_x_0_z_x_y_x, g_y_0_x_0_z_x_y_y, g_y_0_x_0_z_x_y_z, g_yz_x_xy_x, g_yz_x_xy_y, g_yz_x_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_x_y_x[i] = 4.0 * g_yz_x_xy_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_y_y[i] = 4.0 * g_yz_x_xy_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_y_z[i] = 4.0 * g_yz_x_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (303-306)

    #pragma omp simd aligned(g_y_0_x_0_z_x_z_x, g_y_0_x_0_z_x_z_y, g_y_0_x_0_z_x_z_z, g_yz_x_xz_x, g_yz_x_xz_y, g_yz_x_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_x_z_x[i] = 4.0 * g_yz_x_xz_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_z_y[i] = 4.0 * g_yz_x_xz_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_z_z[i] = 4.0 * g_yz_x_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (306-309)

    #pragma omp simd aligned(g_y_0_x_0_z_y_x_x, g_y_0_x_0_z_y_x_y, g_y_0_x_0_z_y_x_z, g_yz_y_0_x, g_yz_y_0_y, g_yz_y_0_z, g_yz_y_xx_x, g_yz_y_xx_y, g_yz_y_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_y_x_x[i] = -2.0 * g_yz_y_0_x[i] * a_exp + 4.0 * g_yz_y_xx_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_x_y[i] = -2.0 * g_yz_y_0_y[i] * a_exp + 4.0 * g_yz_y_xx_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_x_z[i] = -2.0 * g_yz_y_0_z[i] * a_exp + 4.0 * g_yz_y_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (309-312)

    #pragma omp simd aligned(g_y_0_x_0_z_y_y_x, g_y_0_x_0_z_y_y_y, g_y_0_x_0_z_y_y_z, g_yz_y_xy_x, g_yz_y_xy_y, g_yz_y_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_y_y_x[i] = 4.0 * g_yz_y_xy_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_y_y[i] = 4.0 * g_yz_y_xy_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_y_z[i] = 4.0 * g_yz_y_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (312-315)

    #pragma omp simd aligned(g_y_0_x_0_z_y_z_x, g_y_0_x_0_z_y_z_y, g_y_0_x_0_z_y_z_z, g_yz_y_xz_x, g_yz_y_xz_y, g_yz_y_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_y_z_x[i] = 4.0 * g_yz_y_xz_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_z_y[i] = 4.0 * g_yz_y_xz_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_z_z[i] = 4.0 * g_yz_y_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (315-318)

    #pragma omp simd aligned(g_y_0_x_0_z_z_x_x, g_y_0_x_0_z_z_x_y, g_y_0_x_0_z_z_x_z, g_yz_z_0_x, g_yz_z_0_y, g_yz_z_0_z, g_yz_z_xx_x, g_yz_z_xx_y, g_yz_z_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_z_x_x[i] = -2.0 * g_yz_z_0_x[i] * a_exp + 4.0 * g_yz_z_xx_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_x_y[i] = -2.0 * g_yz_z_0_y[i] * a_exp + 4.0 * g_yz_z_xx_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_x_z[i] = -2.0 * g_yz_z_0_z[i] * a_exp + 4.0 * g_yz_z_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (318-321)

    #pragma omp simd aligned(g_y_0_x_0_z_z_y_x, g_y_0_x_0_z_z_y_y, g_y_0_x_0_z_z_y_z, g_yz_z_xy_x, g_yz_z_xy_y, g_yz_z_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_z_y_x[i] = 4.0 * g_yz_z_xy_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_y_y[i] = 4.0 * g_yz_z_xy_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_y_z[i] = 4.0 * g_yz_z_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (321-324)

    #pragma omp simd aligned(g_y_0_x_0_z_z_z_x, g_y_0_x_0_z_z_z_y, g_y_0_x_0_z_z_z_z, g_yz_z_xz_x, g_yz_z_xz_y, g_yz_z_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_z_z_x[i] = 4.0 * g_yz_z_xz_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_z_y[i] = 4.0 * g_yz_z_xz_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_z_z[i] = 4.0 * g_yz_z_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (324-327)

    #pragma omp simd aligned(g_xy_x_xy_x, g_xy_x_xy_y, g_xy_x_xy_z, g_y_0_y_0_x_x_x_x, g_y_0_y_0_x_x_x_y, g_y_0_y_0_x_x_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_x_x_x[i] = 4.0 * g_xy_x_xy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_x_y[i] = 4.0 * g_xy_x_xy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_x_z[i] = 4.0 * g_xy_x_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (327-330)

    #pragma omp simd aligned(g_xy_x_0_x, g_xy_x_0_y, g_xy_x_0_z, g_xy_x_yy_x, g_xy_x_yy_y, g_xy_x_yy_z, g_y_0_y_0_x_x_y_x, g_y_0_y_0_x_x_y_y, g_y_0_y_0_x_x_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_x_y_x[i] = -2.0 * g_xy_x_0_x[i] * a_exp + 4.0 * g_xy_x_yy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_y_y[i] = -2.0 * g_xy_x_0_y[i] * a_exp + 4.0 * g_xy_x_yy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_y_z[i] = -2.0 * g_xy_x_0_z[i] * a_exp + 4.0 * g_xy_x_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (330-333)

    #pragma omp simd aligned(g_xy_x_yz_x, g_xy_x_yz_y, g_xy_x_yz_z, g_y_0_y_0_x_x_z_x, g_y_0_y_0_x_x_z_y, g_y_0_y_0_x_x_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_x_z_x[i] = 4.0 * g_xy_x_yz_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_z_y[i] = 4.0 * g_xy_x_yz_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_z_z[i] = 4.0 * g_xy_x_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (333-336)

    #pragma omp simd aligned(g_xy_y_xy_x, g_xy_y_xy_y, g_xy_y_xy_z, g_y_0_y_0_x_y_x_x, g_y_0_y_0_x_y_x_y, g_y_0_y_0_x_y_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_y_x_x[i] = 4.0 * g_xy_y_xy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_x_y[i] = 4.0 * g_xy_y_xy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_x_z[i] = 4.0 * g_xy_y_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (336-339)

    #pragma omp simd aligned(g_xy_y_0_x, g_xy_y_0_y, g_xy_y_0_z, g_xy_y_yy_x, g_xy_y_yy_y, g_xy_y_yy_z, g_y_0_y_0_x_y_y_x, g_y_0_y_0_x_y_y_y, g_y_0_y_0_x_y_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_y_y_x[i] = -2.0 * g_xy_y_0_x[i] * a_exp + 4.0 * g_xy_y_yy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_y_y[i] = -2.0 * g_xy_y_0_y[i] * a_exp + 4.0 * g_xy_y_yy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_y_z[i] = -2.0 * g_xy_y_0_z[i] * a_exp + 4.0 * g_xy_y_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (339-342)

    #pragma omp simd aligned(g_xy_y_yz_x, g_xy_y_yz_y, g_xy_y_yz_z, g_y_0_y_0_x_y_z_x, g_y_0_y_0_x_y_z_y, g_y_0_y_0_x_y_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_y_z_x[i] = 4.0 * g_xy_y_yz_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_z_y[i] = 4.0 * g_xy_y_yz_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_z_z[i] = 4.0 * g_xy_y_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (342-345)

    #pragma omp simd aligned(g_xy_z_xy_x, g_xy_z_xy_y, g_xy_z_xy_z, g_y_0_y_0_x_z_x_x, g_y_0_y_0_x_z_x_y, g_y_0_y_0_x_z_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_z_x_x[i] = 4.0 * g_xy_z_xy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_x_y[i] = 4.0 * g_xy_z_xy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_x_z[i] = 4.0 * g_xy_z_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (345-348)

    #pragma omp simd aligned(g_xy_z_0_x, g_xy_z_0_y, g_xy_z_0_z, g_xy_z_yy_x, g_xy_z_yy_y, g_xy_z_yy_z, g_y_0_y_0_x_z_y_x, g_y_0_y_0_x_z_y_y, g_y_0_y_0_x_z_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_z_y_x[i] = -2.0 * g_xy_z_0_x[i] * a_exp + 4.0 * g_xy_z_yy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_y_y[i] = -2.0 * g_xy_z_0_y[i] * a_exp + 4.0 * g_xy_z_yy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_y_z[i] = -2.0 * g_xy_z_0_z[i] * a_exp + 4.0 * g_xy_z_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (348-351)

    #pragma omp simd aligned(g_xy_z_yz_x, g_xy_z_yz_y, g_xy_z_yz_z, g_y_0_y_0_x_z_z_x, g_y_0_y_0_x_z_z_y, g_y_0_y_0_x_z_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_z_z_x[i] = 4.0 * g_xy_z_yz_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_z_y[i] = 4.0 * g_xy_z_yz_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_z_z[i] = 4.0 * g_xy_z_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (351-354)

    #pragma omp simd aligned(g_0_x_xy_x, g_0_x_xy_y, g_0_x_xy_z, g_y_0_y_0_y_x_x_x, g_y_0_y_0_y_x_x_y, g_y_0_y_0_y_x_x_z, g_yy_x_xy_x, g_yy_x_xy_y, g_yy_x_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_x_x_x[i] = -2.0 * g_0_x_xy_x[i] * c_exps[i] + 4.0 * g_yy_x_xy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_x_y[i] = -2.0 * g_0_x_xy_y[i] * c_exps[i] + 4.0 * g_yy_x_xy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_x_z[i] = -2.0 * g_0_x_xy_z[i] * c_exps[i] + 4.0 * g_yy_x_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (354-357)

    #pragma omp simd aligned(g_0_x_0_x, g_0_x_0_y, g_0_x_0_z, g_0_x_yy_x, g_0_x_yy_y, g_0_x_yy_z, g_y_0_y_0_y_x_y_x, g_y_0_y_0_y_x_y_y, g_y_0_y_0_y_x_y_z, g_yy_x_0_x, g_yy_x_0_y, g_yy_x_0_z, g_yy_x_yy_x, g_yy_x_yy_y, g_yy_x_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_x_y_x[i] = g_0_x_0_x[i] - 2.0 * g_0_x_yy_x[i] * c_exps[i] - 2.0 * g_yy_x_0_x[i] * a_exp + 4.0 * g_yy_x_yy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_y_y[i] = g_0_x_0_y[i] - 2.0 * g_0_x_yy_y[i] * c_exps[i] - 2.0 * g_yy_x_0_y[i] * a_exp + 4.0 * g_yy_x_yy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_y_z[i] = g_0_x_0_z[i] - 2.0 * g_0_x_yy_z[i] * c_exps[i] - 2.0 * g_yy_x_0_z[i] * a_exp + 4.0 * g_yy_x_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (357-360)

    #pragma omp simd aligned(g_0_x_yz_x, g_0_x_yz_y, g_0_x_yz_z, g_y_0_y_0_y_x_z_x, g_y_0_y_0_y_x_z_y, g_y_0_y_0_y_x_z_z, g_yy_x_yz_x, g_yy_x_yz_y, g_yy_x_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_x_z_x[i] = -2.0 * g_0_x_yz_x[i] * c_exps[i] + 4.0 * g_yy_x_yz_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_z_y[i] = -2.0 * g_0_x_yz_y[i] * c_exps[i] + 4.0 * g_yy_x_yz_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_z_z[i] = -2.0 * g_0_x_yz_z[i] * c_exps[i] + 4.0 * g_yy_x_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (360-363)

    #pragma omp simd aligned(g_0_y_xy_x, g_0_y_xy_y, g_0_y_xy_z, g_y_0_y_0_y_y_x_x, g_y_0_y_0_y_y_x_y, g_y_0_y_0_y_y_x_z, g_yy_y_xy_x, g_yy_y_xy_y, g_yy_y_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_y_x_x[i] = -2.0 * g_0_y_xy_x[i] * c_exps[i] + 4.0 * g_yy_y_xy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_x_y[i] = -2.0 * g_0_y_xy_y[i] * c_exps[i] + 4.0 * g_yy_y_xy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_x_z[i] = -2.0 * g_0_y_xy_z[i] * c_exps[i] + 4.0 * g_yy_y_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (363-366)

    #pragma omp simd aligned(g_0_y_0_x, g_0_y_0_y, g_0_y_0_z, g_0_y_yy_x, g_0_y_yy_y, g_0_y_yy_z, g_y_0_y_0_y_y_y_x, g_y_0_y_0_y_y_y_y, g_y_0_y_0_y_y_y_z, g_yy_y_0_x, g_yy_y_0_y, g_yy_y_0_z, g_yy_y_yy_x, g_yy_y_yy_y, g_yy_y_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_y_y_x[i] = g_0_y_0_x[i] - 2.0 * g_0_y_yy_x[i] * c_exps[i] - 2.0 * g_yy_y_0_x[i] * a_exp + 4.0 * g_yy_y_yy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_y_y[i] = g_0_y_0_y[i] - 2.0 * g_0_y_yy_y[i] * c_exps[i] - 2.0 * g_yy_y_0_y[i] * a_exp + 4.0 * g_yy_y_yy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_y_z[i] = g_0_y_0_z[i] - 2.0 * g_0_y_yy_z[i] * c_exps[i] - 2.0 * g_yy_y_0_z[i] * a_exp + 4.0 * g_yy_y_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (366-369)

    #pragma omp simd aligned(g_0_y_yz_x, g_0_y_yz_y, g_0_y_yz_z, g_y_0_y_0_y_y_z_x, g_y_0_y_0_y_y_z_y, g_y_0_y_0_y_y_z_z, g_yy_y_yz_x, g_yy_y_yz_y, g_yy_y_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_y_z_x[i] = -2.0 * g_0_y_yz_x[i] * c_exps[i] + 4.0 * g_yy_y_yz_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_z_y[i] = -2.0 * g_0_y_yz_y[i] * c_exps[i] + 4.0 * g_yy_y_yz_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_z_z[i] = -2.0 * g_0_y_yz_z[i] * c_exps[i] + 4.0 * g_yy_y_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (369-372)

    #pragma omp simd aligned(g_0_z_xy_x, g_0_z_xy_y, g_0_z_xy_z, g_y_0_y_0_y_z_x_x, g_y_0_y_0_y_z_x_y, g_y_0_y_0_y_z_x_z, g_yy_z_xy_x, g_yy_z_xy_y, g_yy_z_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_z_x_x[i] = -2.0 * g_0_z_xy_x[i] * c_exps[i] + 4.0 * g_yy_z_xy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_x_y[i] = -2.0 * g_0_z_xy_y[i] * c_exps[i] + 4.0 * g_yy_z_xy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_x_z[i] = -2.0 * g_0_z_xy_z[i] * c_exps[i] + 4.0 * g_yy_z_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (372-375)

    #pragma omp simd aligned(g_0_z_0_x, g_0_z_0_y, g_0_z_0_z, g_0_z_yy_x, g_0_z_yy_y, g_0_z_yy_z, g_y_0_y_0_y_z_y_x, g_y_0_y_0_y_z_y_y, g_y_0_y_0_y_z_y_z, g_yy_z_0_x, g_yy_z_0_y, g_yy_z_0_z, g_yy_z_yy_x, g_yy_z_yy_y, g_yy_z_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_z_y_x[i] = g_0_z_0_x[i] - 2.0 * g_0_z_yy_x[i] * c_exps[i] - 2.0 * g_yy_z_0_x[i] * a_exp + 4.0 * g_yy_z_yy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_y_y[i] = g_0_z_0_y[i] - 2.0 * g_0_z_yy_y[i] * c_exps[i] - 2.0 * g_yy_z_0_y[i] * a_exp + 4.0 * g_yy_z_yy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_y_z[i] = g_0_z_0_z[i] - 2.0 * g_0_z_yy_z[i] * c_exps[i] - 2.0 * g_yy_z_0_z[i] * a_exp + 4.0 * g_yy_z_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (375-378)

    #pragma omp simd aligned(g_0_z_yz_x, g_0_z_yz_y, g_0_z_yz_z, g_y_0_y_0_y_z_z_x, g_y_0_y_0_y_z_z_y, g_y_0_y_0_y_z_z_z, g_yy_z_yz_x, g_yy_z_yz_y, g_yy_z_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_z_z_x[i] = -2.0 * g_0_z_yz_x[i] * c_exps[i] + 4.0 * g_yy_z_yz_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_z_y[i] = -2.0 * g_0_z_yz_y[i] * c_exps[i] + 4.0 * g_yy_z_yz_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_z_z[i] = -2.0 * g_0_z_yz_z[i] * c_exps[i] + 4.0 * g_yy_z_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (378-381)

    #pragma omp simd aligned(g_y_0_y_0_z_x_x_x, g_y_0_y_0_z_x_x_y, g_y_0_y_0_z_x_x_z, g_yz_x_xy_x, g_yz_x_xy_y, g_yz_x_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_x_x_x[i] = 4.0 * g_yz_x_xy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_x_y[i] = 4.0 * g_yz_x_xy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_x_z[i] = 4.0 * g_yz_x_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (381-384)

    #pragma omp simd aligned(g_y_0_y_0_z_x_y_x, g_y_0_y_0_z_x_y_y, g_y_0_y_0_z_x_y_z, g_yz_x_0_x, g_yz_x_0_y, g_yz_x_0_z, g_yz_x_yy_x, g_yz_x_yy_y, g_yz_x_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_x_y_x[i] = -2.0 * g_yz_x_0_x[i] * a_exp + 4.0 * g_yz_x_yy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_y_y[i] = -2.0 * g_yz_x_0_y[i] * a_exp + 4.0 * g_yz_x_yy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_y_z[i] = -2.0 * g_yz_x_0_z[i] * a_exp + 4.0 * g_yz_x_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (384-387)

    #pragma omp simd aligned(g_y_0_y_0_z_x_z_x, g_y_0_y_0_z_x_z_y, g_y_0_y_0_z_x_z_z, g_yz_x_yz_x, g_yz_x_yz_y, g_yz_x_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_x_z_x[i] = 4.0 * g_yz_x_yz_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_z_y[i] = 4.0 * g_yz_x_yz_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_z_z[i] = 4.0 * g_yz_x_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (387-390)

    #pragma omp simd aligned(g_y_0_y_0_z_y_x_x, g_y_0_y_0_z_y_x_y, g_y_0_y_0_z_y_x_z, g_yz_y_xy_x, g_yz_y_xy_y, g_yz_y_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_y_x_x[i] = 4.0 * g_yz_y_xy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_x_y[i] = 4.0 * g_yz_y_xy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_x_z[i] = 4.0 * g_yz_y_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (390-393)

    #pragma omp simd aligned(g_y_0_y_0_z_y_y_x, g_y_0_y_0_z_y_y_y, g_y_0_y_0_z_y_y_z, g_yz_y_0_x, g_yz_y_0_y, g_yz_y_0_z, g_yz_y_yy_x, g_yz_y_yy_y, g_yz_y_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_y_y_x[i] = -2.0 * g_yz_y_0_x[i] * a_exp + 4.0 * g_yz_y_yy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_y_y[i] = -2.0 * g_yz_y_0_y[i] * a_exp + 4.0 * g_yz_y_yy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_y_z[i] = -2.0 * g_yz_y_0_z[i] * a_exp + 4.0 * g_yz_y_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (393-396)

    #pragma omp simd aligned(g_y_0_y_0_z_y_z_x, g_y_0_y_0_z_y_z_y, g_y_0_y_0_z_y_z_z, g_yz_y_yz_x, g_yz_y_yz_y, g_yz_y_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_y_z_x[i] = 4.0 * g_yz_y_yz_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_z_y[i] = 4.0 * g_yz_y_yz_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_z_z[i] = 4.0 * g_yz_y_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (396-399)

    #pragma omp simd aligned(g_y_0_y_0_z_z_x_x, g_y_0_y_0_z_z_x_y, g_y_0_y_0_z_z_x_z, g_yz_z_xy_x, g_yz_z_xy_y, g_yz_z_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_z_x_x[i] = 4.0 * g_yz_z_xy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_x_y[i] = 4.0 * g_yz_z_xy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_x_z[i] = 4.0 * g_yz_z_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (399-402)

    #pragma omp simd aligned(g_y_0_y_0_z_z_y_x, g_y_0_y_0_z_z_y_y, g_y_0_y_0_z_z_y_z, g_yz_z_0_x, g_yz_z_0_y, g_yz_z_0_z, g_yz_z_yy_x, g_yz_z_yy_y, g_yz_z_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_z_y_x[i] = -2.0 * g_yz_z_0_x[i] * a_exp + 4.0 * g_yz_z_yy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_y_y[i] = -2.0 * g_yz_z_0_y[i] * a_exp + 4.0 * g_yz_z_yy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_y_z[i] = -2.0 * g_yz_z_0_z[i] * a_exp + 4.0 * g_yz_z_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (402-405)

    #pragma omp simd aligned(g_y_0_y_0_z_z_z_x, g_y_0_y_0_z_z_z_y, g_y_0_y_0_z_z_z_z, g_yz_z_yz_x, g_yz_z_yz_y, g_yz_z_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_z_z_x[i] = 4.0 * g_yz_z_yz_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_z_y[i] = 4.0 * g_yz_z_yz_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_z_z[i] = 4.0 * g_yz_z_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (405-408)

    #pragma omp simd aligned(g_xy_x_xz_x, g_xy_x_xz_y, g_xy_x_xz_z, g_y_0_z_0_x_x_x_x, g_y_0_z_0_x_x_x_y, g_y_0_z_0_x_x_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_x_x_x[i] = 4.0 * g_xy_x_xz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_x_y[i] = 4.0 * g_xy_x_xz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_x_z[i] = 4.0 * g_xy_x_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (408-411)

    #pragma omp simd aligned(g_xy_x_yz_x, g_xy_x_yz_y, g_xy_x_yz_z, g_y_0_z_0_x_x_y_x, g_y_0_z_0_x_x_y_y, g_y_0_z_0_x_x_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_x_y_x[i] = 4.0 * g_xy_x_yz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_y_y[i] = 4.0 * g_xy_x_yz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_y_z[i] = 4.0 * g_xy_x_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (411-414)

    #pragma omp simd aligned(g_xy_x_0_x, g_xy_x_0_y, g_xy_x_0_z, g_xy_x_zz_x, g_xy_x_zz_y, g_xy_x_zz_z, g_y_0_z_0_x_x_z_x, g_y_0_z_0_x_x_z_y, g_y_0_z_0_x_x_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_x_z_x[i] = -2.0 * g_xy_x_0_x[i] * a_exp + 4.0 * g_xy_x_zz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_z_y[i] = -2.0 * g_xy_x_0_y[i] * a_exp + 4.0 * g_xy_x_zz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_z_z[i] = -2.0 * g_xy_x_0_z[i] * a_exp + 4.0 * g_xy_x_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (414-417)

    #pragma omp simd aligned(g_xy_y_xz_x, g_xy_y_xz_y, g_xy_y_xz_z, g_y_0_z_0_x_y_x_x, g_y_0_z_0_x_y_x_y, g_y_0_z_0_x_y_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_y_x_x[i] = 4.0 * g_xy_y_xz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_x_y[i] = 4.0 * g_xy_y_xz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_x_z[i] = 4.0 * g_xy_y_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (417-420)

    #pragma omp simd aligned(g_xy_y_yz_x, g_xy_y_yz_y, g_xy_y_yz_z, g_y_0_z_0_x_y_y_x, g_y_0_z_0_x_y_y_y, g_y_0_z_0_x_y_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_y_y_x[i] = 4.0 * g_xy_y_yz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_y_y[i] = 4.0 * g_xy_y_yz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_y_z[i] = 4.0 * g_xy_y_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (420-423)

    #pragma omp simd aligned(g_xy_y_0_x, g_xy_y_0_y, g_xy_y_0_z, g_xy_y_zz_x, g_xy_y_zz_y, g_xy_y_zz_z, g_y_0_z_0_x_y_z_x, g_y_0_z_0_x_y_z_y, g_y_0_z_0_x_y_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_y_z_x[i] = -2.0 * g_xy_y_0_x[i] * a_exp + 4.0 * g_xy_y_zz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_z_y[i] = -2.0 * g_xy_y_0_y[i] * a_exp + 4.0 * g_xy_y_zz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_z_z[i] = -2.0 * g_xy_y_0_z[i] * a_exp + 4.0 * g_xy_y_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (423-426)

    #pragma omp simd aligned(g_xy_z_xz_x, g_xy_z_xz_y, g_xy_z_xz_z, g_y_0_z_0_x_z_x_x, g_y_0_z_0_x_z_x_y, g_y_0_z_0_x_z_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_z_x_x[i] = 4.0 * g_xy_z_xz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_x_y[i] = 4.0 * g_xy_z_xz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_x_z[i] = 4.0 * g_xy_z_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (426-429)

    #pragma omp simd aligned(g_xy_z_yz_x, g_xy_z_yz_y, g_xy_z_yz_z, g_y_0_z_0_x_z_y_x, g_y_0_z_0_x_z_y_y, g_y_0_z_0_x_z_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_z_y_x[i] = 4.0 * g_xy_z_yz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_y_y[i] = 4.0 * g_xy_z_yz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_y_z[i] = 4.0 * g_xy_z_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (429-432)

    #pragma omp simd aligned(g_xy_z_0_x, g_xy_z_0_y, g_xy_z_0_z, g_xy_z_zz_x, g_xy_z_zz_y, g_xy_z_zz_z, g_y_0_z_0_x_z_z_x, g_y_0_z_0_x_z_z_y, g_y_0_z_0_x_z_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_z_z_x[i] = -2.0 * g_xy_z_0_x[i] * a_exp + 4.0 * g_xy_z_zz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_z_y[i] = -2.0 * g_xy_z_0_y[i] * a_exp + 4.0 * g_xy_z_zz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_z_z[i] = -2.0 * g_xy_z_0_z[i] * a_exp + 4.0 * g_xy_z_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (432-435)

    #pragma omp simd aligned(g_0_x_xz_x, g_0_x_xz_y, g_0_x_xz_z, g_y_0_z_0_y_x_x_x, g_y_0_z_0_y_x_x_y, g_y_0_z_0_y_x_x_z, g_yy_x_xz_x, g_yy_x_xz_y, g_yy_x_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_x_x_x[i] = -2.0 * g_0_x_xz_x[i] * c_exps[i] + 4.0 * g_yy_x_xz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_x_y[i] = -2.0 * g_0_x_xz_y[i] * c_exps[i] + 4.0 * g_yy_x_xz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_x_z[i] = -2.0 * g_0_x_xz_z[i] * c_exps[i] + 4.0 * g_yy_x_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (435-438)

    #pragma omp simd aligned(g_0_x_yz_x, g_0_x_yz_y, g_0_x_yz_z, g_y_0_z_0_y_x_y_x, g_y_0_z_0_y_x_y_y, g_y_0_z_0_y_x_y_z, g_yy_x_yz_x, g_yy_x_yz_y, g_yy_x_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_x_y_x[i] = -2.0 * g_0_x_yz_x[i] * c_exps[i] + 4.0 * g_yy_x_yz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_y_y[i] = -2.0 * g_0_x_yz_y[i] * c_exps[i] + 4.0 * g_yy_x_yz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_y_z[i] = -2.0 * g_0_x_yz_z[i] * c_exps[i] + 4.0 * g_yy_x_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (438-441)

    #pragma omp simd aligned(g_0_x_0_x, g_0_x_0_y, g_0_x_0_z, g_0_x_zz_x, g_0_x_zz_y, g_0_x_zz_z, g_y_0_z_0_y_x_z_x, g_y_0_z_0_y_x_z_y, g_y_0_z_0_y_x_z_z, g_yy_x_0_x, g_yy_x_0_y, g_yy_x_0_z, g_yy_x_zz_x, g_yy_x_zz_y, g_yy_x_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_x_z_x[i] = g_0_x_0_x[i] - 2.0 * g_0_x_zz_x[i] * c_exps[i] - 2.0 * g_yy_x_0_x[i] * a_exp + 4.0 * g_yy_x_zz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_z_y[i] = g_0_x_0_y[i] - 2.0 * g_0_x_zz_y[i] * c_exps[i] - 2.0 * g_yy_x_0_y[i] * a_exp + 4.0 * g_yy_x_zz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_z_z[i] = g_0_x_0_z[i] - 2.0 * g_0_x_zz_z[i] * c_exps[i] - 2.0 * g_yy_x_0_z[i] * a_exp + 4.0 * g_yy_x_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (441-444)

    #pragma omp simd aligned(g_0_y_xz_x, g_0_y_xz_y, g_0_y_xz_z, g_y_0_z_0_y_y_x_x, g_y_0_z_0_y_y_x_y, g_y_0_z_0_y_y_x_z, g_yy_y_xz_x, g_yy_y_xz_y, g_yy_y_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_y_x_x[i] = -2.0 * g_0_y_xz_x[i] * c_exps[i] + 4.0 * g_yy_y_xz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_x_y[i] = -2.0 * g_0_y_xz_y[i] * c_exps[i] + 4.0 * g_yy_y_xz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_x_z[i] = -2.0 * g_0_y_xz_z[i] * c_exps[i] + 4.0 * g_yy_y_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (444-447)

    #pragma omp simd aligned(g_0_y_yz_x, g_0_y_yz_y, g_0_y_yz_z, g_y_0_z_0_y_y_y_x, g_y_0_z_0_y_y_y_y, g_y_0_z_0_y_y_y_z, g_yy_y_yz_x, g_yy_y_yz_y, g_yy_y_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_y_y_x[i] = -2.0 * g_0_y_yz_x[i] * c_exps[i] + 4.0 * g_yy_y_yz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_y_y[i] = -2.0 * g_0_y_yz_y[i] * c_exps[i] + 4.0 * g_yy_y_yz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_y_z[i] = -2.0 * g_0_y_yz_z[i] * c_exps[i] + 4.0 * g_yy_y_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (447-450)

    #pragma omp simd aligned(g_0_y_0_x, g_0_y_0_y, g_0_y_0_z, g_0_y_zz_x, g_0_y_zz_y, g_0_y_zz_z, g_y_0_z_0_y_y_z_x, g_y_0_z_0_y_y_z_y, g_y_0_z_0_y_y_z_z, g_yy_y_0_x, g_yy_y_0_y, g_yy_y_0_z, g_yy_y_zz_x, g_yy_y_zz_y, g_yy_y_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_y_z_x[i] = g_0_y_0_x[i] - 2.0 * g_0_y_zz_x[i] * c_exps[i] - 2.0 * g_yy_y_0_x[i] * a_exp + 4.0 * g_yy_y_zz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_z_y[i] = g_0_y_0_y[i] - 2.0 * g_0_y_zz_y[i] * c_exps[i] - 2.0 * g_yy_y_0_y[i] * a_exp + 4.0 * g_yy_y_zz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_z_z[i] = g_0_y_0_z[i] - 2.0 * g_0_y_zz_z[i] * c_exps[i] - 2.0 * g_yy_y_0_z[i] * a_exp + 4.0 * g_yy_y_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (450-453)

    #pragma omp simd aligned(g_0_z_xz_x, g_0_z_xz_y, g_0_z_xz_z, g_y_0_z_0_y_z_x_x, g_y_0_z_0_y_z_x_y, g_y_0_z_0_y_z_x_z, g_yy_z_xz_x, g_yy_z_xz_y, g_yy_z_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_z_x_x[i] = -2.0 * g_0_z_xz_x[i] * c_exps[i] + 4.0 * g_yy_z_xz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_x_y[i] = -2.0 * g_0_z_xz_y[i] * c_exps[i] + 4.0 * g_yy_z_xz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_x_z[i] = -2.0 * g_0_z_xz_z[i] * c_exps[i] + 4.0 * g_yy_z_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (453-456)

    #pragma omp simd aligned(g_0_z_yz_x, g_0_z_yz_y, g_0_z_yz_z, g_y_0_z_0_y_z_y_x, g_y_0_z_0_y_z_y_y, g_y_0_z_0_y_z_y_z, g_yy_z_yz_x, g_yy_z_yz_y, g_yy_z_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_z_y_x[i] = -2.0 * g_0_z_yz_x[i] * c_exps[i] + 4.0 * g_yy_z_yz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_y_y[i] = -2.0 * g_0_z_yz_y[i] * c_exps[i] + 4.0 * g_yy_z_yz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_y_z[i] = -2.0 * g_0_z_yz_z[i] * c_exps[i] + 4.0 * g_yy_z_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (456-459)

    #pragma omp simd aligned(g_0_z_0_x, g_0_z_0_y, g_0_z_0_z, g_0_z_zz_x, g_0_z_zz_y, g_0_z_zz_z, g_y_0_z_0_y_z_z_x, g_y_0_z_0_y_z_z_y, g_y_0_z_0_y_z_z_z, g_yy_z_0_x, g_yy_z_0_y, g_yy_z_0_z, g_yy_z_zz_x, g_yy_z_zz_y, g_yy_z_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_z_z_x[i] = g_0_z_0_x[i] - 2.0 * g_0_z_zz_x[i] * c_exps[i] - 2.0 * g_yy_z_0_x[i] * a_exp + 4.0 * g_yy_z_zz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_z_y[i] = g_0_z_0_y[i] - 2.0 * g_0_z_zz_y[i] * c_exps[i] - 2.0 * g_yy_z_0_y[i] * a_exp + 4.0 * g_yy_z_zz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_z_z[i] = g_0_z_0_z[i] - 2.0 * g_0_z_zz_z[i] * c_exps[i] - 2.0 * g_yy_z_0_z[i] * a_exp + 4.0 * g_yy_z_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (459-462)

    #pragma omp simd aligned(g_y_0_z_0_z_x_x_x, g_y_0_z_0_z_x_x_y, g_y_0_z_0_z_x_x_z, g_yz_x_xz_x, g_yz_x_xz_y, g_yz_x_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_x_x_x[i] = 4.0 * g_yz_x_xz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_x_y[i] = 4.0 * g_yz_x_xz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_x_z[i] = 4.0 * g_yz_x_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (462-465)

    #pragma omp simd aligned(g_y_0_z_0_z_x_y_x, g_y_0_z_0_z_x_y_y, g_y_0_z_0_z_x_y_z, g_yz_x_yz_x, g_yz_x_yz_y, g_yz_x_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_x_y_x[i] = 4.0 * g_yz_x_yz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_y_y[i] = 4.0 * g_yz_x_yz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_y_z[i] = 4.0 * g_yz_x_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (465-468)

    #pragma omp simd aligned(g_y_0_z_0_z_x_z_x, g_y_0_z_0_z_x_z_y, g_y_0_z_0_z_x_z_z, g_yz_x_0_x, g_yz_x_0_y, g_yz_x_0_z, g_yz_x_zz_x, g_yz_x_zz_y, g_yz_x_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_x_z_x[i] = -2.0 * g_yz_x_0_x[i] * a_exp + 4.0 * g_yz_x_zz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_z_y[i] = -2.0 * g_yz_x_0_y[i] * a_exp + 4.0 * g_yz_x_zz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_z_z[i] = -2.0 * g_yz_x_0_z[i] * a_exp + 4.0 * g_yz_x_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (468-471)

    #pragma omp simd aligned(g_y_0_z_0_z_y_x_x, g_y_0_z_0_z_y_x_y, g_y_0_z_0_z_y_x_z, g_yz_y_xz_x, g_yz_y_xz_y, g_yz_y_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_y_x_x[i] = 4.0 * g_yz_y_xz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_x_y[i] = 4.0 * g_yz_y_xz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_x_z[i] = 4.0 * g_yz_y_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (471-474)

    #pragma omp simd aligned(g_y_0_z_0_z_y_y_x, g_y_0_z_0_z_y_y_y, g_y_0_z_0_z_y_y_z, g_yz_y_yz_x, g_yz_y_yz_y, g_yz_y_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_y_y_x[i] = 4.0 * g_yz_y_yz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_y_y[i] = 4.0 * g_yz_y_yz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_y_z[i] = 4.0 * g_yz_y_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (474-477)

    #pragma omp simd aligned(g_y_0_z_0_z_y_z_x, g_y_0_z_0_z_y_z_y, g_y_0_z_0_z_y_z_z, g_yz_y_0_x, g_yz_y_0_y, g_yz_y_0_z, g_yz_y_zz_x, g_yz_y_zz_y, g_yz_y_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_y_z_x[i] = -2.0 * g_yz_y_0_x[i] * a_exp + 4.0 * g_yz_y_zz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_z_y[i] = -2.0 * g_yz_y_0_y[i] * a_exp + 4.0 * g_yz_y_zz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_z_z[i] = -2.0 * g_yz_y_0_z[i] * a_exp + 4.0 * g_yz_y_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (477-480)

    #pragma omp simd aligned(g_y_0_z_0_z_z_x_x, g_y_0_z_0_z_z_x_y, g_y_0_z_0_z_z_x_z, g_yz_z_xz_x, g_yz_z_xz_y, g_yz_z_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_z_x_x[i] = 4.0 * g_yz_z_xz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_x_y[i] = 4.0 * g_yz_z_xz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_x_z[i] = 4.0 * g_yz_z_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (480-483)

    #pragma omp simd aligned(g_y_0_z_0_z_z_y_x, g_y_0_z_0_z_z_y_y, g_y_0_z_0_z_z_y_z, g_yz_z_yz_x, g_yz_z_yz_y, g_yz_z_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_z_y_x[i] = 4.0 * g_yz_z_yz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_y_y[i] = 4.0 * g_yz_z_yz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_y_z[i] = 4.0 * g_yz_z_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (483-486)

    #pragma omp simd aligned(g_y_0_z_0_z_z_z_x, g_y_0_z_0_z_z_z_y, g_y_0_z_0_z_z_z_z, g_yz_z_0_x, g_yz_z_0_y, g_yz_z_0_z, g_yz_z_zz_x, g_yz_z_zz_y, g_yz_z_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_z_z_x[i] = -2.0 * g_yz_z_0_x[i] * a_exp + 4.0 * g_yz_z_zz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_z_y[i] = -2.0 * g_yz_z_0_y[i] * a_exp + 4.0 * g_yz_z_zz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_z_z[i] = -2.0 * g_yz_z_0_z[i] * a_exp + 4.0 * g_yz_z_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (486-489)

    #pragma omp simd aligned(g_xz_x_0_x, g_xz_x_0_y, g_xz_x_0_z, g_xz_x_xx_x, g_xz_x_xx_y, g_xz_x_xx_z, g_z_0_x_0_x_x_x_x, g_z_0_x_0_x_x_x_y, g_z_0_x_0_x_x_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_x_x_x[i] = -2.0 * g_xz_x_0_x[i] * a_exp + 4.0 * g_xz_x_xx_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_x_y[i] = -2.0 * g_xz_x_0_y[i] * a_exp + 4.0 * g_xz_x_xx_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_x_z[i] = -2.0 * g_xz_x_0_z[i] * a_exp + 4.0 * g_xz_x_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (489-492)

    #pragma omp simd aligned(g_xz_x_xy_x, g_xz_x_xy_y, g_xz_x_xy_z, g_z_0_x_0_x_x_y_x, g_z_0_x_0_x_x_y_y, g_z_0_x_0_x_x_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_x_y_x[i] = 4.0 * g_xz_x_xy_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_y_y[i] = 4.0 * g_xz_x_xy_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_y_z[i] = 4.0 * g_xz_x_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (492-495)

    #pragma omp simd aligned(g_xz_x_xz_x, g_xz_x_xz_y, g_xz_x_xz_z, g_z_0_x_0_x_x_z_x, g_z_0_x_0_x_x_z_y, g_z_0_x_0_x_x_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_x_z_x[i] = 4.0 * g_xz_x_xz_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_z_y[i] = 4.0 * g_xz_x_xz_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_z_z[i] = 4.0 * g_xz_x_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (495-498)

    #pragma omp simd aligned(g_xz_y_0_x, g_xz_y_0_y, g_xz_y_0_z, g_xz_y_xx_x, g_xz_y_xx_y, g_xz_y_xx_z, g_z_0_x_0_x_y_x_x, g_z_0_x_0_x_y_x_y, g_z_0_x_0_x_y_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_y_x_x[i] = -2.0 * g_xz_y_0_x[i] * a_exp + 4.0 * g_xz_y_xx_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_x_y[i] = -2.0 * g_xz_y_0_y[i] * a_exp + 4.0 * g_xz_y_xx_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_x_z[i] = -2.0 * g_xz_y_0_z[i] * a_exp + 4.0 * g_xz_y_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (498-501)

    #pragma omp simd aligned(g_xz_y_xy_x, g_xz_y_xy_y, g_xz_y_xy_z, g_z_0_x_0_x_y_y_x, g_z_0_x_0_x_y_y_y, g_z_0_x_0_x_y_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_y_y_x[i] = 4.0 * g_xz_y_xy_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_y_y[i] = 4.0 * g_xz_y_xy_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_y_z[i] = 4.0 * g_xz_y_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (501-504)

    #pragma omp simd aligned(g_xz_y_xz_x, g_xz_y_xz_y, g_xz_y_xz_z, g_z_0_x_0_x_y_z_x, g_z_0_x_0_x_y_z_y, g_z_0_x_0_x_y_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_y_z_x[i] = 4.0 * g_xz_y_xz_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_z_y[i] = 4.0 * g_xz_y_xz_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_z_z[i] = 4.0 * g_xz_y_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (504-507)

    #pragma omp simd aligned(g_xz_z_0_x, g_xz_z_0_y, g_xz_z_0_z, g_xz_z_xx_x, g_xz_z_xx_y, g_xz_z_xx_z, g_z_0_x_0_x_z_x_x, g_z_0_x_0_x_z_x_y, g_z_0_x_0_x_z_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_z_x_x[i] = -2.0 * g_xz_z_0_x[i] * a_exp + 4.0 * g_xz_z_xx_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_x_y[i] = -2.0 * g_xz_z_0_y[i] * a_exp + 4.0 * g_xz_z_xx_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_x_z[i] = -2.0 * g_xz_z_0_z[i] * a_exp + 4.0 * g_xz_z_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (507-510)

    #pragma omp simd aligned(g_xz_z_xy_x, g_xz_z_xy_y, g_xz_z_xy_z, g_z_0_x_0_x_z_y_x, g_z_0_x_0_x_z_y_y, g_z_0_x_0_x_z_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_z_y_x[i] = 4.0 * g_xz_z_xy_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_y_y[i] = 4.0 * g_xz_z_xy_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_y_z[i] = 4.0 * g_xz_z_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (510-513)

    #pragma omp simd aligned(g_xz_z_xz_x, g_xz_z_xz_y, g_xz_z_xz_z, g_z_0_x_0_x_z_z_x, g_z_0_x_0_x_z_z_y, g_z_0_x_0_x_z_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_z_z_x[i] = 4.0 * g_xz_z_xz_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_z_y[i] = 4.0 * g_xz_z_xz_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_z_z[i] = 4.0 * g_xz_z_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (513-516)

    #pragma omp simd aligned(g_yz_x_0_x, g_yz_x_0_y, g_yz_x_0_z, g_yz_x_xx_x, g_yz_x_xx_y, g_yz_x_xx_z, g_z_0_x_0_y_x_x_x, g_z_0_x_0_y_x_x_y, g_z_0_x_0_y_x_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_x_x_x[i] = -2.0 * g_yz_x_0_x[i] * a_exp + 4.0 * g_yz_x_xx_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_x_y[i] = -2.0 * g_yz_x_0_y[i] * a_exp + 4.0 * g_yz_x_xx_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_x_z[i] = -2.0 * g_yz_x_0_z[i] * a_exp + 4.0 * g_yz_x_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (516-519)

    #pragma omp simd aligned(g_yz_x_xy_x, g_yz_x_xy_y, g_yz_x_xy_z, g_z_0_x_0_y_x_y_x, g_z_0_x_0_y_x_y_y, g_z_0_x_0_y_x_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_x_y_x[i] = 4.0 * g_yz_x_xy_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_y_y[i] = 4.0 * g_yz_x_xy_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_y_z[i] = 4.0 * g_yz_x_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (519-522)

    #pragma omp simd aligned(g_yz_x_xz_x, g_yz_x_xz_y, g_yz_x_xz_z, g_z_0_x_0_y_x_z_x, g_z_0_x_0_y_x_z_y, g_z_0_x_0_y_x_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_x_z_x[i] = 4.0 * g_yz_x_xz_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_z_y[i] = 4.0 * g_yz_x_xz_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_z_z[i] = 4.0 * g_yz_x_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (522-525)

    #pragma omp simd aligned(g_yz_y_0_x, g_yz_y_0_y, g_yz_y_0_z, g_yz_y_xx_x, g_yz_y_xx_y, g_yz_y_xx_z, g_z_0_x_0_y_y_x_x, g_z_0_x_0_y_y_x_y, g_z_0_x_0_y_y_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_y_x_x[i] = -2.0 * g_yz_y_0_x[i] * a_exp + 4.0 * g_yz_y_xx_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_x_y[i] = -2.0 * g_yz_y_0_y[i] * a_exp + 4.0 * g_yz_y_xx_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_x_z[i] = -2.0 * g_yz_y_0_z[i] * a_exp + 4.0 * g_yz_y_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (525-528)

    #pragma omp simd aligned(g_yz_y_xy_x, g_yz_y_xy_y, g_yz_y_xy_z, g_z_0_x_0_y_y_y_x, g_z_0_x_0_y_y_y_y, g_z_0_x_0_y_y_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_y_y_x[i] = 4.0 * g_yz_y_xy_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_y_y[i] = 4.0 * g_yz_y_xy_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_y_z[i] = 4.0 * g_yz_y_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (528-531)

    #pragma omp simd aligned(g_yz_y_xz_x, g_yz_y_xz_y, g_yz_y_xz_z, g_z_0_x_0_y_y_z_x, g_z_0_x_0_y_y_z_y, g_z_0_x_0_y_y_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_y_z_x[i] = 4.0 * g_yz_y_xz_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_z_y[i] = 4.0 * g_yz_y_xz_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_z_z[i] = 4.0 * g_yz_y_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (531-534)

    #pragma omp simd aligned(g_yz_z_0_x, g_yz_z_0_y, g_yz_z_0_z, g_yz_z_xx_x, g_yz_z_xx_y, g_yz_z_xx_z, g_z_0_x_0_y_z_x_x, g_z_0_x_0_y_z_x_y, g_z_0_x_0_y_z_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_z_x_x[i] = -2.0 * g_yz_z_0_x[i] * a_exp + 4.0 * g_yz_z_xx_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_x_y[i] = -2.0 * g_yz_z_0_y[i] * a_exp + 4.0 * g_yz_z_xx_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_x_z[i] = -2.0 * g_yz_z_0_z[i] * a_exp + 4.0 * g_yz_z_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (534-537)

    #pragma omp simd aligned(g_yz_z_xy_x, g_yz_z_xy_y, g_yz_z_xy_z, g_z_0_x_0_y_z_y_x, g_z_0_x_0_y_z_y_y, g_z_0_x_0_y_z_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_z_y_x[i] = 4.0 * g_yz_z_xy_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_y_y[i] = 4.0 * g_yz_z_xy_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_y_z[i] = 4.0 * g_yz_z_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (537-540)

    #pragma omp simd aligned(g_yz_z_xz_x, g_yz_z_xz_y, g_yz_z_xz_z, g_z_0_x_0_y_z_z_x, g_z_0_x_0_y_z_z_y, g_z_0_x_0_y_z_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_z_z_x[i] = 4.0 * g_yz_z_xz_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_z_y[i] = 4.0 * g_yz_z_xz_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_z_z[i] = 4.0 * g_yz_z_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (540-543)

    #pragma omp simd aligned(g_0_x_0_x, g_0_x_0_y, g_0_x_0_z, g_0_x_xx_x, g_0_x_xx_y, g_0_x_xx_z, g_z_0_x_0_z_x_x_x, g_z_0_x_0_z_x_x_y, g_z_0_x_0_z_x_x_z, g_zz_x_0_x, g_zz_x_0_y, g_zz_x_0_z, g_zz_x_xx_x, g_zz_x_xx_y, g_zz_x_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_x_x_x[i] = g_0_x_0_x[i] - 2.0 * g_0_x_xx_x[i] * c_exps[i] - 2.0 * g_zz_x_0_x[i] * a_exp + 4.0 * g_zz_x_xx_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_x_y[i] = g_0_x_0_y[i] - 2.0 * g_0_x_xx_y[i] * c_exps[i] - 2.0 * g_zz_x_0_y[i] * a_exp + 4.0 * g_zz_x_xx_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_x_z[i] = g_0_x_0_z[i] - 2.0 * g_0_x_xx_z[i] * c_exps[i] - 2.0 * g_zz_x_0_z[i] * a_exp + 4.0 * g_zz_x_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (543-546)

    #pragma omp simd aligned(g_0_x_xy_x, g_0_x_xy_y, g_0_x_xy_z, g_z_0_x_0_z_x_y_x, g_z_0_x_0_z_x_y_y, g_z_0_x_0_z_x_y_z, g_zz_x_xy_x, g_zz_x_xy_y, g_zz_x_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_x_y_x[i] = -2.0 * g_0_x_xy_x[i] * c_exps[i] + 4.0 * g_zz_x_xy_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_y_y[i] = -2.0 * g_0_x_xy_y[i] * c_exps[i] + 4.0 * g_zz_x_xy_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_y_z[i] = -2.0 * g_0_x_xy_z[i] * c_exps[i] + 4.0 * g_zz_x_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (546-549)

    #pragma omp simd aligned(g_0_x_xz_x, g_0_x_xz_y, g_0_x_xz_z, g_z_0_x_0_z_x_z_x, g_z_0_x_0_z_x_z_y, g_z_0_x_0_z_x_z_z, g_zz_x_xz_x, g_zz_x_xz_y, g_zz_x_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_x_z_x[i] = -2.0 * g_0_x_xz_x[i] * c_exps[i] + 4.0 * g_zz_x_xz_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_z_y[i] = -2.0 * g_0_x_xz_y[i] * c_exps[i] + 4.0 * g_zz_x_xz_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_z_z[i] = -2.0 * g_0_x_xz_z[i] * c_exps[i] + 4.0 * g_zz_x_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (549-552)

    #pragma omp simd aligned(g_0_y_0_x, g_0_y_0_y, g_0_y_0_z, g_0_y_xx_x, g_0_y_xx_y, g_0_y_xx_z, g_z_0_x_0_z_y_x_x, g_z_0_x_0_z_y_x_y, g_z_0_x_0_z_y_x_z, g_zz_y_0_x, g_zz_y_0_y, g_zz_y_0_z, g_zz_y_xx_x, g_zz_y_xx_y, g_zz_y_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_y_x_x[i] = g_0_y_0_x[i] - 2.0 * g_0_y_xx_x[i] * c_exps[i] - 2.0 * g_zz_y_0_x[i] * a_exp + 4.0 * g_zz_y_xx_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_x_y[i] = g_0_y_0_y[i] - 2.0 * g_0_y_xx_y[i] * c_exps[i] - 2.0 * g_zz_y_0_y[i] * a_exp + 4.0 * g_zz_y_xx_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_x_z[i] = g_0_y_0_z[i] - 2.0 * g_0_y_xx_z[i] * c_exps[i] - 2.0 * g_zz_y_0_z[i] * a_exp + 4.0 * g_zz_y_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (552-555)

    #pragma omp simd aligned(g_0_y_xy_x, g_0_y_xy_y, g_0_y_xy_z, g_z_0_x_0_z_y_y_x, g_z_0_x_0_z_y_y_y, g_z_0_x_0_z_y_y_z, g_zz_y_xy_x, g_zz_y_xy_y, g_zz_y_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_y_y_x[i] = -2.0 * g_0_y_xy_x[i] * c_exps[i] + 4.0 * g_zz_y_xy_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_y_y[i] = -2.0 * g_0_y_xy_y[i] * c_exps[i] + 4.0 * g_zz_y_xy_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_y_z[i] = -2.0 * g_0_y_xy_z[i] * c_exps[i] + 4.0 * g_zz_y_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (555-558)

    #pragma omp simd aligned(g_0_y_xz_x, g_0_y_xz_y, g_0_y_xz_z, g_z_0_x_0_z_y_z_x, g_z_0_x_0_z_y_z_y, g_z_0_x_0_z_y_z_z, g_zz_y_xz_x, g_zz_y_xz_y, g_zz_y_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_y_z_x[i] = -2.0 * g_0_y_xz_x[i] * c_exps[i] + 4.0 * g_zz_y_xz_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_z_y[i] = -2.0 * g_0_y_xz_y[i] * c_exps[i] + 4.0 * g_zz_y_xz_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_z_z[i] = -2.0 * g_0_y_xz_z[i] * c_exps[i] + 4.0 * g_zz_y_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (558-561)

    #pragma omp simd aligned(g_0_z_0_x, g_0_z_0_y, g_0_z_0_z, g_0_z_xx_x, g_0_z_xx_y, g_0_z_xx_z, g_z_0_x_0_z_z_x_x, g_z_0_x_0_z_z_x_y, g_z_0_x_0_z_z_x_z, g_zz_z_0_x, g_zz_z_0_y, g_zz_z_0_z, g_zz_z_xx_x, g_zz_z_xx_y, g_zz_z_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_z_x_x[i] = g_0_z_0_x[i] - 2.0 * g_0_z_xx_x[i] * c_exps[i] - 2.0 * g_zz_z_0_x[i] * a_exp + 4.0 * g_zz_z_xx_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_x_y[i] = g_0_z_0_y[i] - 2.0 * g_0_z_xx_y[i] * c_exps[i] - 2.0 * g_zz_z_0_y[i] * a_exp + 4.0 * g_zz_z_xx_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_x_z[i] = g_0_z_0_z[i] - 2.0 * g_0_z_xx_z[i] * c_exps[i] - 2.0 * g_zz_z_0_z[i] * a_exp + 4.0 * g_zz_z_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (561-564)

    #pragma omp simd aligned(g_0_z_xy_x, g_0_z_xy_y, g_0_z_xy_z, g_z_0_x_0_z_z_y_x, g_z_0_x_0_z_z_y_y, g_z_0_x_0_z_z_y_z, g_zz_z_xy_x, g_zz_z_xy_y, g_zz_z_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_z_y_x[i] = -2.0 * g_0_z_xy_x[i] * c_exps[i] + 4.0 * g_zz_z_xy_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_y_y[i] = -2.0 * g_0_z_xy_y[i] * c_exps[i] + 4.0 * g_zz_z_xy_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_y_z[i] = -2.0 * g_0_z_xy_z[i] * c_exps[i] + 4.0 * g_zz_z_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (564-567)

    #pragma omp simd aligned(g_0_z_xz_x, g_0_z_xz_y, g_0_z_xz_z, g_z_0_x_0_z_z_z_x, g_z_0_x_0_z_z_z_y, g_z_0_x_0_z_z_z_z, g_zz_z_xz_x, g_zz_z_xz_y, g_zz_z_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_z_z_x[i] = -2.0 * g_0_z_xz_x[i] * c_exps[i] + 4.0 * g_zz_z_xz_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_z_y[i] = -2.0 * g_0_z_xz_y[i] * c_exps[i] + 4.0 * g_zz_z_xz_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_z_z[i] = -2.0 * g_0_z_xz_z[i] * c_exps[i] + 4.0 * g_zz_z_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (567-570)

    #pragma omp simd aligned(g_xz_x_xy_x, g_xz_x_xy_y, g_xz_x_xy_z, g_z_0_y_0_x_x_x_x, g_z_0_y_0_x_x_x_y, g_z_0_y_0_x_x_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_x_x_x[i] = 4.0 * g_xz_x_xy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_x_y[i] = 4.0 * g_xz_x_xy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_x_z[i] = 4.0 * g_xz_x_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (570-573)

    #pragma omp simd aligned(g_xz_x_0_x, g_xz_x_0_y, g_xz_x_0_z, g_xz_x_yy_x, g_xz_x_yy_y, g_xz_x_yy_z, g_z_0_y_0_x_x_y_x, g_z_0_y_0_x_x_y_y, g_z_0_y_0_x_x_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_x_y_x[i] = -2.0 * g_xz_x_0_x[i] * a_exp + 4.0 * g_xz_x_yy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_y_y[i] = -2.0 * g_xz_x_0_y[i] * a_exp + 4.0 * g_xz_x_yy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_y_z[i] = -2.0 * g_xz_x_0_z[i] * a_exp + 4.0 * g_xz_x_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (573-576)

    #pragma omp simd aligned(g_xz_x_yz_x, g_xz_x_yz_y, g_xz_x_yz_z, g_z_0_y_0_x_x_z_x, g_z_0_y_0_x_x_z_y, g_z_0_y_0_x_x_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_x_z_x[i] = 4.0 * g_xz_x_yz_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_z_y[i] = 4.0 * g_xz_x_yz_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_z_z[i] = 4.0 * g_xz_x_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (576-579)

    #pragma omp simd aligned(g_xz_y_xy_x, g_xz_y_xy_y, g_xz_y_xy_z, g_z_0_y_0_x_y_x_x, g_z_0_y_0_x_y_x_y, g_z_0_y_0_x_y_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_y_x_x[i] = 4.0 * g_xz_y_xy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_x_y[i] = 4.0 * g_xz_y_xy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_x_z[i] = 4.0 * g_xz_y_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (579-582)

    #pragma omp simd aligned(g_xz_y_0_x, g_xz_y_0_y, g_xz_y_0_z, g_xz_y_yy_x, g_xz_y_yy_y, g_xz_y_yy_z, g_z_0_y_0_x_y_y_x, g_z_0_y_0_x_y_y_y, g_z_0_y_0_x_y_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_y_y_x[i] = -2.0 * g_xz_y_0_x[i] * a_exp + 4.0 * g_xz_y_yy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_y_y[i] = -2.0 * g_xz_y_0_y[i] * a_exp + 4.0 * g_xz_y_yy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_y_z[i] = -2.0 * g_xz_y_0_z[i] * a_exp + 4.0 * g_xz_y_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (582-585)

    #pragma omp simd aligned(g_xz_y_yz_x, g_xz_y_yz_y, g_xz_y_yz_z, g_z_0_y_0_x_y_z_x, g_z_0_y_0_x_y_z_y, g_z_0_y_0_x_y_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_y_z_x[i] = 4.0 * g_xz_y_yz_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_z_y[i] = 4.0 * g_xz_y_yz_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_z_z[i] = 4.0 * g_xz_y_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (585-588)

    #pragma omp simd aligned(g_xz_z_xy_x, g_xz_z_xy_y, g_xz_z_xy_z, g_z_0_y_0_x_z_x_x, g_z_0_y_0_x_z_x_y, g_z_0_y_0_x_z_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_z_x_x[i] = 4.0 * g_xz_z_xy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_x_y[i] = 4.0 * g_xz_z_xy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_x_z[i] = 4.0 * g_xz_z_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (588-591)

    #pragma omp simd aligned(g_xz_z_0_x, g_xz_z_0_y, g_xz_z_0_z, g_xz_z_yy_x, g_xz_z_yy_y, g_xz_z_yy_z, g_z_0_y_0_x_z_y_x, g_z_0_y_0_x_z_y_y, g_z_0_y_0_x_z_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_z_y_x[i] = -2.0 * g_xz_z_0_x[i] * a_exp + 4.0 * g_xz_z_yy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_y_y[i] = -2.0 * g_xz_z_0_y[i] * a_exp + 4.0 * g_xz_z_yy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_y_z[i] = -2.0 * g_xz_z_0_z[i] * a_exp + 4.0 * g_xz_z_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (591-594)

    #pragma omp simd aligned(g_xz_z_yz_x, g_xz_z_yz_y, g_xz_z_yz_z, g_z_0_y_0_x_z_z_x, g_z_0_y_0_x_z_z_y, g_z_0_y_0_x_z_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_z_z_x[i] = 4.0 * g_xz_z_yz_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_z_y[i] = 4.0 * g_xz_z_yz_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_z_z[i] = 4.0 * g_xz_z_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (594-597)

    #pragma omp simd aligned(g_yz_x_xy_x, g_yz_x_xy_y, g_yz_x_xy_z, g_z_0_y_0_y_x_x_x, g_z_0_y_0_y_x_x_y, g_z_0_y_0_y_x_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_x_x_x[i] = 4.0 * g_yz_x_xy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_x_y[i] = 4.0 * g_yz_x_xy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_x_z[i] = 4.0 * g_yz_x_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (597-600)

    #pragma omp simd aligned(g_yz_x_0_x, g_yz_x_0_y, g_yz_x_0_z, g_yz_x_yy_x, g_yz_x_yy_y, g_yz_x_yy_z, g_z_0_y_0_y_x_y_x, g_z_0_y_0_y_x_y_y, g_z_0_y_0_y_x_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_x_y_x[i] = -2.0 * g_yz_x_0_x[i] * a_exp + 4.0 * g_yz_x_yy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_y_y[i] = -2.0 * g_yz_x_0_y[i] * a_exp + 4.0 * g_yz_x_yy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_y_z[i] = -2.0 * g_yz_x_0_z[i] * a_exp + 4.0 * g_yz_x_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (600-603)

    #pragma omp simd aligned(g_yz_x_yz_x, g_yz_x_yz_y, g_yz_x_yz_z, g_z_0_y_0_y_x_z_x, g_z_0_y_0_y_x_z_y, g_z_0_y_0_y_x_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_x_z_x[i] = 4.0 * g_yz_x_yz_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_z_y[i] = 4.0 * g_yz_x_yz_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_z_z[i] = 4.0 * g_yz_x_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (603-606)

    #pragma omp simd aligned(g_yz_y_xy_x, g_yz_y_xy_y, g_yz_y_xy_z, g_z_0_y_0_y_y_x_x, g_z_0_y_0_y_y_x_y, g_z_0_y_0_y_y_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_y_x_x[i] = 4.0 * g_yz_y_xy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_x_y[i] = 4.0 * g_yz_y_xy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_x_z[i] = 4.0 * g_yz_y_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (606-609)

    #pragma omp simd aligned(g_yz_y_0_x, g_yz_y_0_y, g_yz_y_0_z, g_yz_y_yy_x, g_yz_y_yy_y, g_yz_y_yy_z, g_z_0_y_0_y_y_y_x, g_z_0_y_0_y_y_y_y, g_z_0_y_0_y_y_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_y_y_x[i] = -2.0 * g_yz_y_0_x[i] * a_exp + 4.0 * g_yz_y_yy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_y_y[i] = -2.0 * g_yz_y_0_y[i] * a_exp + 4.0 * g_yz_y_yy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_y_z[i] = -2.0 * g_yz_y_0_z[i] * a_exp + 4.0 * g_yz_y_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (609-612)

    #pragma omp simd aligned(g_yz_y_yz_x, g_yz_y_yz_y, g_yz_y_yz_z, g_z_0_y_0_y_y_z_x, g_z_0_y_0_y_y_z_y, g_z_0_y_0_y_y_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_y_z_x[i] = 4.0 * g_yz_y_yz_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_z_y[i] = 4.0 * g_yz_y_yz_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_z_z[i] = 4.0 * g_yz_y_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (612-615)

    #pragma omp simd aligned(g_yz_z_xy_x, g_yz_z_xy_y, g_yz_z_xy_z, g_z_0_y_0_y_z_x_x, g_z_0_y_0_y_z_x_y, g_z_0_y_0_y_z_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_z_x_x[i] = 4.0 * g_yz_z_xy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_x_y[i] = 4.0 * g_yz_z_xy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_x_z[i] = 4.0 * g_yz_z_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (615-618)

    #pragma omp simd aligned(g_yz_z_0_x, g_yz_z_0_y, g_yz_z_0_z, g_yz_z_yy_x, g_yz_z_yy_y, g_yz_z_yy_z, g_z_0_y_0_y_z_y_x, g_z_0_y_0_y_z_y_y, g_z_0_y_0_y_z_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_z_y_x[i] = -2.0 * g_yz_z_0_x[i] * a_exp + 4.0 * g_yz_z_yy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_y_y[i] = -2.0 * g_yz_z_0_y[i] * a_exp + 4.0 * g_yz_z_yy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_y_z[i] = -2.0 * g_yz_z_0_z[i] * a_exp + 4.0 * g_yz_z_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (618-621)

    #pragma omp simd aligned(g_yz_z_yz_x, g_yz_z_yz_y, g_yz_z_yz_z, g_z_0_y_0_y_z_z_x, g_z_0_y_0_y_z_z_y, g_z_0_y_0_y_z_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_z_z_x[i] = 4.0 * g_yz_z_yz_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_z_y[i] = 4.0 * g_yz_z_yz_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_z_z[i] = 4.0 * g_yz_z_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (621-624)

    #pragma omp simd aligned(g_0_x_xy_x, g_0_x_xy_y, g_0_x_xy_z, g_z_0_y_0_z_x_x_x, g_z_0_y_0_z_x_x_y, g_z_0_y_0_z_x_x_z, g_zz_x_xy_x, g_zz_x_xy_y, g_zz_x_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_x_x_x[i] = -2.0 * g_0_x_xy_x[i] * c_exps[i] + 4.0 * g_zz_x_xy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_x_y[i] = -2.0 * g_0_x_xy_y[i] * c_exps[i] + 4.0 * g_zz_x_xy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_x_z[i] = -2.0 * g_0_x_xy_z[i] * c_exps[i] + 4.0 * g_zz_x_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (624-627)

    #pragma omp simd aligned(g_0_x_0_x, g_0_x_0_y, g_0_x_0_z, g_0_x_yy_x, g_0_x_yy_y, g_0_x_yy_z, g_z_0_y_0_z_x_y_x, g_z_0_y_0_z_x_y_y, g_z_0_y_0_z_x_y_z, g_zz_x_0_x, g_zz_x_0_y, g_zz_x_0_z, g_zz_x_yy_x, g_zz_x_yy_y, g_zz_x_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_x_y_x[i] = g_0_x_0_x[i] - 2.0 * g_0_x_yy_x[i] * c_exps[i] - 2.0 * g_zz_x_0_x[i] * a_exp + 4.0 * g_zz_x_yy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_y_y[i] = g_0_x_0_y[i] - 2.0 * g_0_x_yy_y[i] * c_exps[i] - 2.0 * g_zz_x_0_y[i] * a_exp + 4.0 * g_zz_x_yy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_y_z[i] = g_0_x_0_z[i] - 2.0 * g_0_x_yy_z[i] * c_exps[i] - 2.0 * g_zz_x_0_z[i] * a_exp + 4.0 * g_zz_x_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (627-630)

    #pragma omp simd aligned(g_0_x_yz_x, g_0_x_yz_y, g_0_x_yz_z, g_z_0_y_0_z_x_z_x, g_z_0_y_0_z_x_z_y, g_z_0_y_0_z_x_z_z, g_zz_x_yz_x, g_zz_x_yz_y, g_zz_x_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_x_z_x[i] = -2.0 * g_0_x_yz_x[i] * c_exps[i] + 4.0 * g_zz_x_yz_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_z_y[i] = -2.0 * g_0_x_yz_y[i] * c_exps[i] + 4.0 * g_zz_x_yz_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_z_z[i] = -2.0 * g_0_x_yz_z[i] * c_exps[i] + 4.0 * g_zz_x_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (630-633)

    #pragma omp simd aligned(g_0_y_xy_x, g_0_y_xy_y, g_0_y_xy_z, g_z_0_y_0_z_y_x_x, g_z_0_y_0_z_y_x_y, g_z_0_y_0_z_y_x_z, g_zz_y_xy_x, g_zz_y_xy_y, g_zz_y_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_y_x_x[i] = -2.0 * g_0_y_xy_x[i] * c_exps[i] + 4.0 * g_zz_y_xy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_x_y[i] = -2.0 * g_0_y_xy_y[i] * c_exps[i] + 4.0 * g_zz_y_xy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_x_z[i] = -2.0 * g_0_y_xy_z[i] * c_exps[i] + 4.0 * g_zz_y_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (633-636)

    #pragma omp simd aligned(g_0_y_0_x, g_0_y_0_y, g_0_y_0_z, g_0_y_yy_x, g_0_y_yy_y, g_0_y_yy_z, g_z_0_y_0_z_y_y_x, g_z_0_y_0_z_y_y_y, g_z_0_y_0_z_y_y_z, g_zz_y_0_x, g_zz_y_0_y, g_zz_y_0_z, g_zz_y_yy_x, g_zz_y_yy_y, g_zz_y_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_y_y_x[i] = g_0_y_0_x[i] - 2.0 * g_0_y_yy_x[i] * c_exps[i] - 2.0 * g_zz_y_0_x[i] * a_exp + 4.0 * g_zz_y_yy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_y_y[i] = g_0_y_0_y[i] - 2.0 * g_0_y_yy_y[i] * c_exps[i] - 2.0 * g_zz_y_0_y[i] * a_exp + 4.0 * g_zz_y_yy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_y_z[i] = g_0_y_0_z[i] - 2.0 * g_0_y_yy_z[i] * c_exps[i] - 2.0 * g_zz_y_0_z[i] * a_exp + 4.0 * g_zz_y_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (636-639)

    #pragma omp simd aligned(g_0_y_yz_x, g_0_y_yz_y, g_0_y_yz_z, g_z_0_y_0_z_y_z_x, g_z_0_y_0_z_y_z_y, g_z_0_y_0_z_y_z_z, g_zz_y_yz_x, g_zz_y_yz_y, g_zz_y_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_y_z_x[i] = -2.0 * g_0_y_yz_x[i] * c_exps[i] + 4.0 * g_zz_y_yz_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_z_y[i] = -2.0 * g_0_y_yz_y[i] * c_exps[i] + 4.0 * g_zz_y_yz_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_z_z[i] = -2.0 * g_0_y_yz_z[i] * c_exps[i] + 4.0 * g_zz_y_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (639-642)

    #pragma omp simd aligned(g_0_z_xy_x, g_0_z_xy_y, g_0_z_xy_z, g_z_0_y_0_z_z_x_x, g_z_0_y_0_z_z_x_y, g_z_0_y_0_z_z_x_z, g_zz_z_xy_x, g_zz_z_xy_y, g_zz_z_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_z_x_x[i] = -2.0 * g_0_z_xy_x[i] * c_exps[i] + 4.0 * g_zz_z_xy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_x_y[i] = -2.0 * g_0_z_xy_y[i] * c_exps[i] + 4.0 * g_zz_z_xy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_x_z[i] = -2.0 * g_0_z_xy_z[i] * c_exps[i] + 4.0 * g_zz_z_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (642-645)

    #pragma omp simd aligned(g_0_z_0_x, g_0_z_0_y, g_0_z_0_z, g_0_z_yy_x, g_0_z_yy_y, g_0_z_yy_z, g_z_0_y_0_z_z_y_x, g_z_0_y_0_z_z_y_y, g_z_0_y_0_z_z_y_z, g_zz_z_0_x, g_zz_z_0_y, g_zz_z_0_z, g_zz_z_yy_x, g_zz_z_yy_y, g_zz_z_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_z_y_x[i] = g_0_z_0_x[i] - 2.0 * g_0_z_yy_x[i] * c_exps[i] - 2.0 * g_zz_z_0_x[i] * a_exp + 4.0 * g_zz_z_yy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_y_y[i] = g_0_z_0_y[i] - 2.0 * g_0_z_yy_y[i] * c_exps[i] - 2.0 * g_zz_z_0_y[i] * a_exp + 4.0 * g_zz_z_yy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_y_z[i] = g_0_z_0_z[i] - 2.0 * g_0_z_yy_z[i] * c_exps[i] - 2.0 * g_zz_z_0_z[i] * a_exp + 4.0 * g_zz_z_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (645-648)

    #pragma omp simd aligned(g_0_z_yz_x, g_0_z_yz_y, g_0_z_yz_z, g_z_0_y_0_z_z_z_x, g_z_0_y_0_z_z_z_y, g_z_0_y_0_z_z_z_z, g_zz_z_yz_x, g_zz_z_yz_y, g_zz_z_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_z_z_x[i] = -2.0 * g_0_z_yz_x[i] * c_exps[i] + 4.0 * g_zz_z_yz_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_z_y[i] = -2.0 * g_0_z_yz_y[i] * c_exps[i] + 4.0 * g_zz_z_yz_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_z_z[i] = -2.0 * g_0_z_yz_z[i] * c_exps[i] + 4.0 * g_zz_z_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (648-651)

    #pragma omp simd aligned(g_xz_x_xz_x, g_xz_x_xz_y, g_xz_x_xz_z, g_z_0_z_0_x_x_x_x, g_z_0_z_0_x_x_x_y, g_z_0_z_0_x_x_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_x_x_x[i] = 4.0 * g_xz_x_xz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_x_y[i] = 4.0 * g_xz_x_xz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_x_z[i] = 4.0 * g_xz_x_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (651-654)

    #pragma omp simd aligned(g_xz_x_yz_x, g_xz_x_yz_y, g_xz_x_yz_z, g_z_0_z_0_x_x_y_x, g_z_0_z_0_x_x_y_y, g_z_0_z_0_x_x_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_x_y_x[i] = 4.0 * g_xz_x_yz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_y_y[i] = 4.0 * g_xz_x_yz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_y_z[i] = 4.0 * g_xz_x_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (654-657)

    #pragma omp simd aligned(g_xz_x_0_x, g_xz_x_0_y, g_xz_x_0_z, g_xz_x_zz_x, g_xz_x_zz_y, g_xz_x_zz_z, g_z_0_z_0_x_x_z_x, g_z_0_z_0_x_x_z_y, g_z_0_z_0_x_x_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_x_z_x[i] = -2.0 * g_xz_x_0_x[i] * a_exp + 4.0 * g_xz_x_zz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_z_y[i] = -2.0 * g_xz_x_0_y[i] * a_exp + 4.0 * g_xz_x_zz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_z_z[i] = -2.0 * g_xz_x_0_z[i] * a_exp + 4.0 * g_xz_x_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (657-660)

    #pragma omp simd aligned(g_xz_y_xz_x, g_xz_y_xz_y, g_xz_y_xz_z, g_z_0_z_0_x_y_x_x, g_z_0_z_0_x_y_x_y, g_z_0_z_0_x_y_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_y_x_x[i] = 4.0 * g_xz_y_xz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_x_y[i] = 4.0 * g_xz_y_xz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_x_z[i] = 4.0 * g_xz_y_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (660-663)

    #pragma omp simd aligned(g_xz_y_yz_x, g_xz_y_yz_y, g_xz_y_yz_z, g_z_0_z_0_x_y_y_x, g_z_0_z_0_x_y_y_y, g_z_0_z_0_x_y_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_y_y_x[i] = 4.0 * g_xz_y_yz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_y_y[i] = 4.0 * g_xz_y_yz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_y_z[i] = 4.0 * g_xz_y_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (663-666)

    #pragma omp simd aligned(g_xz_y_0_x, g_xz_y_0_y, g_xz_y_0_z, g_xz_y_zz_x, g_xz_y_zz_y, g_xz_y_zz_z, g_z_0_z_0_x_y_z_x, g_z_0_z_0_x_y_z_y, g_z_0_z_0_x_y_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_y_z_x[i] = -2.0 * g_xz_y_0_x[i] * a_exp + 4.0 * g_xz_y_zz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_z_y[i] = -2.0 * g_xz_y_0_y[i] * a_exp + 4.0 * g_xz_y_zz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_z_z[i] = -2.0 * g_xz_y_0_z[i] * a_exp + 4.0 * g_xz_y_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (666-669)

    #pragma omp simd aligned(g_xz_z_xz_x, g_xz_z_xz_y, g_xz_z_xz_z, g_z_0_z_0_x_z_x_x, g_z_0_z_0_x_z_x_y, g_z_0_z_0_x_z_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_z_x_x[i] = 4.0 * g_xz_z_xz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_x_y[i] = 4.0 * g_xz_z_xz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_x_z[i] = 4.0 * g_xz_z_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (669-672)

    #pragma omp simd aligned(g_xz_z_yz_x, g_xz_z_yz_y, g_xz_z_yz_z, g_z_0_z_0_x_z_y_x, g_z_0_z_0_x_z_y_y, g_z_0_z_0_x_z_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_z_y_x[i] = 4.0 * g_xz_z_yz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_y_y[i] = 4.0 * g_xz_z_yz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_y_z[i] = 4.0 * g_xz_z_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (672-675)

    #pragma omp simd aligned(g_xz_z_0_x, g_xz_z_0_y, g_xz_z_0_z, g_xz_z_zz_x, g_xz_z_zz_y, g_xz_z_zz_z, g_z_0_z_0_x_z_z_x, g_z_0_z_0_x_z_z_y, g_z_0_z_0_x_z_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_z_z_x[i] = -2.0 * g_xz_z_0_x[i] * a_exp + 4.0 * g_xz_z_zz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_z_y[i] = -2.0 * g_xz_z_0_y[i] * a_exp + 4.0 * g_xz_z_zz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_z_z[i] = -2.0 * g_xz_z_0_z[i] * a_exp + 4.0 * g_xz_z_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (675-678)

    #pragma omp simd aligned(g_yz_x_xz_x, g_yz_x_xz_y, g_yz_x_xz_z, g_z_0_z_0_y_x_x_x, g_z_0_z_0_y_x_x_y, g_z_0_z_0_y_x_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_x_x_x[i] = 4.0 * g_yz_x_xz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_x_y[i] = 4.0 * g_yz_x_xz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_x_z[i] = 4.0 * g_yz_x_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (678-681)

    #pragma omp simd aligned(g_yz_x_yz_x, g_yz_x_yz_y, g_yz_x_yz_z, g_z_0_z_0_y_x_y_x, g_z_0_z_0_y_x_y_y, g_z_0_z_0_y_x_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_x_y_x[i] = 4.0 * g_yz_x_yz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_y_y[i] = 4.0 * g_yz_x_yz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_y_z[i] = 4.0 * g_yz_x_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (681-684)

    #pragma omp simd aligned(g_yz_x_0_x, g_yz_x_0_y, g_yz_x_0_z, g_yz_x_zz_x, g_yz_x_zz_y, g_yz_x_zz_z, g_z_0_z_0_y_x_z_x, g_z_0_z_0_y_x_z_y, g_z_0_z_0_y_x_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_x_z_x[i] = -2.0 * g_yz_x_0_x[i] * a_exp + 4.0 * g_yz_x_zz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_z_y[i] = -2.0 * g_yz_x_0_y[i] * a_exp + 4.0 * g_yz_x_zz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_z_z[i] = -2.0 * g_yz_x_0_z[i] * a_exp + 4.0 * g_yz_x_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (684-687)

    #pragma omp simd aligned(g_yz_y_xz_x, g_yz_y_xz_y, g_yz_y_xz_z, g_z_0_z_0_y_y_x_x, g_z_0_z_0_y_y_x_y, g_z_0_z_0_y_y_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_y_x_x[i] = 4.0 * g_yz_y_xz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_x_y[i] = 4.0 * g_yz_y_xz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_x_z[i] = 4.0 * g_yz_y_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (687-690)

    #pragma omp simd aligned(g_yz_y_yz_x, g_yz_y_yz_y, g_yz_y_yz_z, g_z_0_z_0_y_y_y_x, g_z_0_z_0_y_y_y_y, g_z_0_z_0_y_y_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_y_y_x[i] = 4.0 * g_yz_y_yz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_y_y[i] = 4.0 * g_yz_y_yz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_y_z[i] = 4.0 * g_yz_y_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (690-693)

    #pragma omp simd aligned(g_yz_y_0_x, g_yz_y_0_y, g_yz_y_0_z, g_yz_y_zz_x, g_yz_y_zz_y, g_yz_y_zz_z, g_z_0_z_0_y_y_z_x, g_z_0_z_0_y_y_z_y, g_z_0_z_0_y_y_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_y_z_x[i] = -2.0 * g_yz_y_0_x[i] * a_exp + 4.0 * g_yz_y_zz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_z_y[i] = -2.0 * g_yz_y_0_y[i] * a_exp + 4.0 * g_yz_y_zz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_z_z[i] = -2.0 * g_yz_y_0_z[i] * a_exp + 4.0 * g_yz_y_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (693-696)

    #pragma omp simd aligned(g_yz_z_xz_x, g_yz_z_xz_y, g_yz_z_xz_z, g_z_0_z_0_y_z_x_x, g_z_0_z_0_y_z_x_y, g_z_0_z_0_y_z_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_z_x_x[i] = 4.0 * g_yz_z_xz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_x_y[i] = 4.0 * g_yz_z_xz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_x_z[i] = 4.0 * g_yz_z_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (696-699)

    #pragma omp simd aligned(g_yz_z_yz_x, g_yz_z_yz_y, g_yz_z_yz_z, g_z_0_z_0_y_z_y_x, g_z_0_z_0_y_z_y_y, g_z_0_z_0_y_z_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_z_y_x[i] = 4.0 * g_yz_z_yz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_y_y[i] = 4.0 * g_yz_z_yz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_y_z[i] = 4.0 * g_yz_z_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (699-702)

    #pragma omp simd aligned(g_yz_z_0_x, g_yz_z_0_y, g_yz_z_0_z, g_yz_z_zz_x, g_yz_z_zz_y, g_yz_z_zz_z, g_z_0_z_0_y_z_z_x, g_z_0_z_0_y_z_z_y, g_z_0_z_0_y_z_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_z_z_x[i] = -2.0 * g_yz_z_0_x[i] * a_exp + 4.0 * g_yz_z_zz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_z_y[i] = -2.0 * g_yz_z_0_y[i] * a_exp + 4.0 * g_yz_z_zz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_z_z[i] = -2.0 * g_yz_z_0_z[i] * a_exp + 4.0 * g_yz_z_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (702-705)

    #pragma omp simd aligned(g_0_x_xz_x, g_0_x_xz_y, g_0_x_xz_z, g_z_0_z_0_z_x_x_x, g_z_0_z_0_z_x_x_y, g_z_0_z_0_z_x_x_z, g_zz_x_xz_x, g_zz_x_xz_y, g_zz_x_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_x_x_x[i] = -2.0 * g_0_x_xz_x[i] * c_exps[i] + 4.0 * g_zz_x_xz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_x_y[i] = -2.0 * g_0_x_xz_y[i] * c_exps[i] + 4.0 * g_zz_x_xz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_x_z[i] = -2.0 * g_0_x_xz_z[i] * c_exps[i] + 4.0 * g_zz_x_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (705-708)

    #pragma omp simd aligned(g_0_x_yz_x, g_0_x_yz_y, g_0_x_yz_z, g_z_0_z_0_z_x_y_x, g_z_0_z_0_z_x_y_y, g_z_0_z_0_z_x_y_z, g_zz_x_yz_x, g_zz_x_yz_y, g_zz_x_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_x_y_x[i] = -2.0 * g_0_x_yz_x[i] * c_exps[i] + 4.0 * g_zz_x_yz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_y_y[i] = -2.0 * g_0_x_yz_y[i] * c_exps[i] + 4.0 * g_zz_x_yz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_y_z[i] = -2.0 * g_0_x_yz_z[i] * c_exps[i] + 4.0 * g_zz_x_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (708-711)

    #pragma omp simd aligned(g_0_x_0_x, g_0_x_0_y, g_0_x_0_z, g_0_x_zz_x, g_0_x_zz_y, g_0_x_zz_z, g_z_0_z_0_z_x_z_x, g_z_0_z_0_z_x_z_y, g_z_0_z_0_z_x_z_z, g_zz_x_0_x, g_zz_x_0_y, g_zz_x_0_z, g_zz_x_zz_x, g_zz_x_zz_y, g_zz_x_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_x_z_x[i] = g_0_x_0_x[i] - 2.0 * g_0_x_zz_x[i] * c_exps[i] - 2.0 * g_zz_x_0_x[i] * a_exp + 4.0 * g_zz_x_zz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_z_y[i] = g_0_x_0_y[i] - 2.0 * g_0_x_zz_y[i] * c_exps[i] - 2.0 * g_zz_x_0_y[i] * a_exp + 4.0 * g_zz_x_zz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_z_z[i] = g_0_x_0_z[i] - 2.0 * g_0_x_zz_z[i] * c_exps[i] - 2.0 * g_zz_x_0_z[i] * a_exp + 4.0 * g_zz_x_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (711-714)

    #pragma omp simd aligned(g_0_y_xz_x, g_0_y_xz_y, g_0_y_xz_z, g_z_0_z_0_z_y_x_x, g_z_0_z_0_z_y_x_y, g_z_0_z_0_z_y_x_z, g_zz_y_xz_x, g_zz_y_xz_y, g_zz_y_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_y_x_x[i] = -2.0 * g_0_y_xz_x[i] * c_exps[i] + 4.0 * g_zz_y_xz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_x_y[i] = -2.0 * g_0_y_xz_y[i] * c_exps[i] + 4.0 * g_zz_y_xz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_x_z[i] = -2.0 * g_0_y_xz_z[i] * c_exps[i] + 4.0 * g_zz_y_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (714-717)

    #pragma omp simd aligned(g_0_y_yz_x, g_0_y_yz_y, g_0_y_yz_z, g_z_0_z_0_z_y_y_x, g_z_0_z_0_z_y_y_y, g_z_0_z_0_z_y_y_z, g_zz_y_yz_x, g_zz_y_yz_y, g_zz_y_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_y_y_x[i] = -2.0 * g_0_y_yz_x[i] * c_exps[i] + 4.0 * g_zz_y_yz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_y_y[i] = -2.0 * g_0_y_yz_y[i] * c_exps[i] + 4.0 * g_zz_y_yz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_y_z[i] = -2.0 * g_0_y_yz_z[i] * c_exps[i] + 4.0 * g_zz_y_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (717-720)

    #pragma omp simd aligned(g_0_y_0_x, g_0_y_0_y, g_0_y_0_z, g_0_y_zz_x, g_0_y_zz_y, g_0_y_zz_z, g_z_0_z_0_z_y_z_x, g_z_0_z_0_z_y_z_y, g_z_0_z_0_z_y_z_z, g_zz_y_0_x, g_zz_y_0_y, g_zz_y_0_z, g_zz_y_zz_x, g_zz_y_zz_y, g_zz_y_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_y_z_x[i] = g_0_y_0_x[i] - 2.0 * g_0_y_zz_x[i] * c_exps[i] - 2.0 * g_zz_y_0_x[i] * a_exp + 4.0 * g_zz_y_zz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_z_y[i] = g_0_y_0_y[i] - 2.0 * g_0_y_zz_y[i] * c_exps[i] - 2.0 * g_zz_y_0_y[i] * a_exp + 4.0 * g_zz_y_zz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_z_z[i] = g_0_y_0_z[i] - 2.0 * g_0_y_zz_z[i] * c_exps[i] - 2.0 * g_zz_y_0_z[i] * a_exp + 4.0 * g_zz_y_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (720-723)

    #pragma omp simd aligned(g_0_z_xz_x, g_0_z_xz_y, g_0_z_xz_z, g_z_0_z_0_z_z_x_x, g_z_0_z_0_z_z_x_y, g_z_0_z_0_z_z_x_z, g_zz_z_xz_x, g_zz_z_xz_y, g_zz_z_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_z_x_x[i] = -2.0 * g_0_z_xz_x[i] * c_exps[i] + 4.0 * g_zz_z_xz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_x_y[i] = -2.0 * g_0_z_xz_y[i] * c_exps[i] + 4.0 * g_zz_z_xz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_x_z[i] = -2.0 * g_0_z_xz_z[i] * c_exps[i] + 4.0 * g_zz_z_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (723-726)

    #pragma omp simd aligned(g_0_z_yz_x, g_0_z_yz_y, g_0_z_yz_z, g_z_0_z_0_z_z_y_x, g_z_0_z_0_z_z_y_y, g_z_0_z_0_z_z_y_z, g_zz_z_yz_x, g_zz_z_yz_y, g_zz_z_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_z_y_x[i] = -2.0 * g_0_z_yz_x[i] * c_exps[i] + 4.0 * g_zz_z_yz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_y_y[i] = -2.0 * g_0_z_yz_y[i] * c_exps[i] + 4.0 * g_zz_z_yz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_y_z[i] = -2.0 * g_0_z_yz_z[i] * c_exps[i] + 4.0 * g_zz_z_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (726-729)

    #pragma omp simd aligned(g_0_z_0_x, g_0_z_0_y, g_0_z_0_z, g_0_z_zz_x, g_0_z_zz_y, g_0_z_zz_z, g_z_0_z_0_z_z_z_x, g_z_0_z_0_z_z_z_y, g_z_0_z_0_z_z_z_z, g_zz_z_0_x, g_zz_z_0_y, g_zz_z_0_z, g_zz_z_zz_x, g_zz_z_zz_y, g_zz_z_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_z_z_x[i] = g_0_z_0_x[i] - 2.0 * g_0_z_zz_x[i] * c_exps[i] - 2.0 * g_zz_z_0_x[i] * a_exp + 4.0 * g_zz_z_zz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_z_y[i] = g_0_z_0_y[i] - 2.0 * g_0_z_zz_y[i] * c_exps[i] - 2.0 * g_zz_z_0_y[i] * a_exp + 4.0 * g_zz_z_zz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_z_z[i] = g_0_z_0_z[i] - 2.0 * g_0_z_zz_z[i] * c_exps[i] - 2.0 * g_zz_z_0_z[i] * a_exp + 4.0 * g_zz_z_zz_z[i] * a_exp * c_exps[i];
    }
}

} // t4c_geom namespace


#include "GeomDeriv1010OfScalarForPDSS.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1010_pdss_0(CSimdArray<double>& buffer_1010_pdss,
                     const CSimdArray<double>& buffer_sdps,
                     const CSimdArray<double>& buffer_ddps,
                     const double a_exp,
                     const double* c_exps) -> void
{
    const auto ndims = buffer_1010_pdss.number_of_columns();

    /// Set up components of auxilary buffer : buffer_sdps

    auto g_0_xx_x_0 = buffer_sdps[0];

    auto g_0_xx_y_0 = buffer_sdps[1];

    auto g_0_xx_z_0 = buffer_sdps[2];

    auto g_0_xy_x_0 = buffer_sdps[3];

    auto g_0_xy_y_0 = buffer_sdps[4];

    auto g_0_xy_z_0 = buffer_sdps[5];

    auto g_0_xz_x_0 = buffer_sdps[6];

    auto g_0_xz_y_0 = buffer_sdps[7];

    auto g_0_xz_z_0 = buffer_sdps[8];

    auto g_0_yy_x_0 = buffer_sdps[9];

    auto g_0_yy_y_0 = buffer_sdps[10];

    auto g_0_yy_z_0 = buffer_sdps[11];

    auto g_0_yz_x_0 = buffer_sdps[12];

    auto g_0_yz_y_0 = buffer_sdps[13];

    auto g_0_yz_z_0 = buffer_sdps[14];

    auto g_0_zz_x_0 = buffer_sdps[15];

    auto g_0_zz_y_0 = buffer_sdps[16];

    auto g_0_zz_z_0 = buffer_sdps[17];

    /// Set up components of auxilary buffer : buffer_ddps

    auto g_xx_xx_x_0 = buffer_ddps[0];

    auto g_xx_xx_y_0 = buffer_ddps[1];

    auto g_xx_xx_z_0 = buffer_ddps[2];

    auto g_xx_xy_x_0 = buffer_ddps[3];

    auto g_xx_xy_y_0 = buffer_ddps[4];

    auto g_xx_xy_z_0 = buffer_ddps[5];

    auto g_xx_xz_x_0 = buffer_ddps[6];

    auto g_xx_xz_y_0 = buffer_ddps[7];

    auto g_xx_xz_z_0 = buffer_ddps[8];

    auto g_xx_yy_x_0 = buffer_ddps[9];

    auto g_xx_yy_y_0 = buffer_ddps[10];

    auto g_xx_yy_z_0 = buffer_ddps[11];

    auto g_xx_yz_x_0 = buffer_ddps[12];

    auto g_xx_yz_y_0 = buffer_ddps[13];

    auto g_xx_yz_z_0 = buffer_ddps[14];

    auto g_xx_zz_x_0 = buffer_ddps[15];

    auto g_xx_zz_y_0 = buffer_ddps[16];

    auto g_xx_zz_z_0 = buffer_ddps[17];

    auto g_xy_xx_x_0 = buffer_ddps[18];

    auto g_xy_xx_y_0 = buffer_ddps[19];

    auto g_xy_xx_z_0 = buffer_ddps[20];

    auto g_xy_xy_x_0 = buffer_ddps[21];

    auto g_xy_xy_y_0 = buffer_ddps[22];

    auto g_xy_xy_z_0 = buffer_ddps[23];

    auto g_xy_xz_x_0 = buffer_ddps[24];

    auto g_xy_xz_y_0 = buffer_ddps[25];

    auto g_xy_xz_z_0 = buffer_ddps[26];

    auto g_xy_yy_x_0 = buffer_ddps[27];

    auto g_xy_yy_y_0 = buffer_ddps[28];

    auto g_xy_yy_z_0 = buffer_ddps[29];

    auto g_xy_yz_x_0 = buffer_ddps[30];

    auto g_xy_yz_y_0 = buffer_ddps[31];

    auto g_xy_yz_z_0 = buffer_ddps[32];

    auto g_xy_zz_x_0 = buffer_ddps[33];

    auto g_xy_zz_y_0 = buffer_ddps[34];

    auto g_xy_zz_z_0 = buffer_ddps[35];

    auto g_xz_xx_x_0 = buffer_ddps[36];

    auto g_xz_xx_y_0 = buffer_ddps[37];

    auto g_xz_xx_z_0 = buffer_ddps[38];

    auto g_xz_xy_x_0 = buffer_ddps[39];

    auto g_xz_xy_y_0 = buffer_ddps[40];

    auto g_xz_xy_z_0 = buffer_ddps[41];

    auto g_xz_xz_x_0 = buffer_ddps[42];

    auto g_xz_xz_y_0 = buffer_ddps[43];

    auto g_xz_xz_z_0 = buffer_ddps[44];

    auto g_xz_yy_x_0 = buffer_ddps[45];

    auto g_xz_yy_y_0 = buffer_ddps[46];

    auto g_xz_yy_z_0 = buffer_ddps[47];

    auto g_xz_yz_x_0 = buffer_ddps[48];

    auto g_xz_yz_y_0 = buffer_ddps[49];

    auto g_xz_yz_z_0 = buffer_ddps[50];

    auto g_xz_zz_x_0 = buffer_ddps[51];

    auto g_xz_zz_y_0 = buffer_ddps[52];

    auto g_xz_zz_z_0 = buffer_ddps[53];

    auto g_yy_xx_x_0 = buffer_ddps[54];

    auto g_yy_xx_y_0 = buffer_ddps[55];

    auto g_yy_xx_z_0 = buffer_ddps[56];

    auto g_yy_xy_x_0 = buffer_ddps[57];

    auto g_yy_xy_y_0 = buffer_ddps[58];

    auto g_yy_xy_z_0 = buffer_ddps[59];

    auto g_yy_xz_x_0 = buffer_ddps[60];

    auto g_yy_xz_y_0 = buffer_ddps[61];

    auto g_yy_xz_z_0 = buffer_ddps[62];

    auto g_yy_yy_x_0 = buffer_ddps[63];

    auto g_yy_yy_y_0 = buffer_ddps[64];

    auto g_yy_yy_z_0 = buffer_ddps[65];

    auto g_yy_yz_x_0 = buffer_ddps[66];

    auto g_yy_yz_y_0 = buffer_ddps[67];

    auto g_yy_yz_z_0 = buffer_ddps[68];

    auto g_yy_zz_x_0 = buffer_ddps[69];

    auto g_yy_zz_y_0 = buffer_ddps[70];

    auto g_yy_zz_z_0 = buffer_ddps[71];

    auto g_yz_xx_x_0 = buffer_ddps[72];

    auto g_yz_xx_y_0 = buffer_ddps[73];

    auto g_yz_xx_z_0 = buffer_ddps[74];

    auto g_yz_xy_x_0 = buffer_ddps[75];

    auto g_yz_xy_y_0 = buffer_ddps[76];

    auto g_yz_xy_z_0 = buffer_ddps[77];

    auto g_yz_xz_x_0 = buffer_ddps[78];

    auto g_yz_xz_y_0 = buffer_ddps[79];

    auto g_yz_xz_z_0 = buffer_ddps[80];

    auto g_yz_yy_x_0 = buffer_ddps[81];

    auto g_yz_yy_y_0 = buffer_ddps[82];

    auto g_yz_yy_z_0 = buffer_ddps[83];

    auto g_yz_yz_x_0 = buffer_ddps[84];

    auto g_yz_yz_y_0 = buffer_ddps[85];

    auto g_yz_yz_z_0 = buffer_ddps[86];

    auto g_yz_zz_x_0 = buffer_ddps[87];

    auto g_yz_zz_y_0 = buffer_ddps[88];

    auto g_yz_zz_z_0 = buffer_ddps[89];

    auto g_zz_xx_x_0 = buffer_ddps[90];

    auto g_zz_xx_y_0 = buffer_ddps[91];

    auto g_zz_xx_z_0 = buffer_ddps[92];

    auto g_zz_xy_x_0 = buffer_ddps[93];

    auto g_zz_xy_y_0 = buffer_ddps[94];

    auto g_zz_xy_z_0 = buffer_ddps[95];

    auto g_zz_xz_x_0 = buffer_ddps[96];

    auto g_zz_xz_y_0 = buffer_ddps[97];

    auto g_zz_xz_z_0 = buffer_ddps[98];

    auto g_zz_yy_x_0 = buffer_ddps[99];

    auto g_zz_yy_y_0 = buffer_ddps[100];

    auto g_zz_yy_z_0 = buffer_ddps[101];

    auto g_zz_yz_x_0 = buffer_ddps[102];

    auto g_zz_yz_y_0 = buffer_ddps[103];

    auto g_zz_yz_z_0 = buffer_ddps[104];

    auto g_zz_zz_x_0 = buffer_ddps[105];

    auto g_zz_zz_y_0 = buffer_ddps[106];

    auto g_zz_zz_z_0 = buffer_ddps[107];

    /// Set up components of integrals buffer : buffer_1010_pdss

    auto g_x_0_x_0_x_xx_0_0 = buffer_1010_pdss[0];

    auto g_x_0_x_0_x_xy_0_0 = buffer_1010_pdss[1];

    auto g_x_0_x_0_x_xz_0_0 = buffer_1010_pdss[2];

    auto g_x_0_x_0_x_yy_0_0 = buffer_1010_pdss[3];

    auto g_x_0_x_0_x_yz_0_0 = buffer_1010_pdss[4];

    auto g_x_0_x_0_x_zz_0_0 = buffer_1010_pdss[5];

    auto g_x_0_x_0_y_xx_0_0 = buffer_1010_pdss[6];

    auto g_x_0_x_0_y_xy_0_0 = buffer_1010_pdss[7];

    auto g_x_0_x_0_y_xz_0_0 = buffer_1010_pdss[8];

    auto g_x_0_x_0_y_yy_0_0 = buffer_1010_pdss[9];

    auto g_x_0_x_0_y_yz_0_0 = buffer_1010_pdss[10];

    auto g_x_0_x_0_y_zz_0_0 = buffer_1010_pdss[11];

    auto g_x_0_x_0_z_xx_0_0 = buffer_1010_pdss[12];

    auto g_x_0_x_0_z_xy_0_0 = buffer_1010_pdss[13];

    auto g_x_0_x_0_z_xz_0_0 = buffer_1010_pdss[14];

    auto g_x_0_x_0_z_yy_0_0 = buffer_1010_pdss[15];

    auto g_x_0_x_0_z_yz_0_0 = buffer_1010_pdss[16];

    auto g_x_0_x_0_z_zz_0_0 = buffer_1010_pdss[17];

    auto g_x_0_y_0_x_xx_0_0 = buffer_1010_pdss[18];

    auto g_x_0_y_0_x_xy_0_0 = buffer_1010_pdss[19];

    auto g_x_0_y_0_x_xz_0_0 = buffer_1010_pdss[20];

    auto g_x_0_y_0_x_yy_0_0 = buffer_1010_pdss[21];

    auto g_x_0_y_0_x_yz_0_0 = buffer_1010_pdss[22];

    auto g_x_0_y_0_x_zz_0_0 = buffer_1010_pdss[23];

    auto g_x_0_y_0_y_xx_0_0 = buffer_1010_pdss[24];

    auto g_x_0_y_0_y_xy_0_0 = buffer_1010_pdss[25];

    auto g_x_0_y_0_y_xz_0_0 = buffer_1010_pdss[26];

    auto g_x_0_y_0_y_yy_0_0 = buffer_1010_pdss[27];

    auto g_x_0_y_0_y_yz_0_0 = buffer_1010_pdss[28];

    auto g_x_0_y_0_y_zz_0_0 = buffer_1010_pdss[29];

    auto g_x_0_y_0_z_xx_0_0 = buffer_1010_pdss[30];

    auto g_x_0_y_0_z_xy_0_0 = buffer_1010_pdss[31];

    auto g_x_0_y_0_z_xz_0_0 = buffer_1010_pdss[32];

    auto g_x_0_y_0_z_yy_0_0 = buffer_1010_pdss[33];

    auto g_x_0_y_0_z_yz_0_0 = buffer_1010_pdss[34];

    auto g_x_0_y_0_z_zz_0_0 = buffer_1010_pdss[35];

    auto g_x_0_z_0_x_xx_0_0 = buffer_1010_pdss[36];

    auto g_x_0_z_0_x_xy_0_0 = buffer_1010_pdss[37];

    auto g_x_0_z_0_x_xz_0_0 = buffer_1010_pdss[38];

    auto g_x_0_z_0_x_yy_0_0 = buffer_1010_pdss[39];

    auto g_x_0_z_0_x_yz_0_0 = buffer_1010_pdss[40];

    auto g_x_0_z_0_x_zz_0_0 = buffer_1010_pdss[41];

    auto g_x_0_z_0_y_xx_0_0 = buffer_1010_pdss[42];

    auto g_x_0_z_0_y_xy_0_0 = buffer_1010_pdss[43];

    auto g_x_0_z_0_y_xz_0_0 = buffer_1010_pdss[44];

    auto g_x_0_z_0_y_yy_0_0 = buffer_1010_pdss[45];

    auto g_x_0_z_0_y_yz_0_0 = buffer_1010_pdss[46];

    auto g_x_0_z_0_y_zz_0_0 = buffer_1010_pdss[47];

    auto g_x_0_z_0_z_xx_0_0 = buffer_1010_pdss[48];

    auto g_x_0_z_0_z_xy_0_0 = buffer_1010_pdss[49];

    auto g_x_0_z_0_z_xz_0_0 = buffer_1010_pdss[50];

    auto g_x_0_z_0_z_yy_0_0 = buffer_1010_pdss[51];

    auto g_x_0_z_0_z_yz_0_0 = buffer_1010_pdss[52];

    auto g_x_0_z_0_z_zz_0_0 = buffer_1010_pdss[53];

    auto g_y_0_x_0_x_xx_0_0 = buffer_1010_pdss[54];

    auto g_y_0_x_0_x_xy_0_0 = buffer_1010_pdss[55];

    auto g_y_0_x_0_x_xz_0_0 = buffer_1010_pdss[56];

    auto g_y_0_x_0_x_yy_0_0 = buffer_1010_pdss[57];

    auto g_y_0_x_0_x_yz_0_0 = buffer_1010_pdss[58];

    auto g_y_0_x_0_x_zz_0_0 = buffer_1010_pdss[59];

    auto g_y_0_x_0_y_xx_0_0 = buffer_1010_pdss[60];

    auto g_y_0_x_0_y_xy_0_0 = buffer_1010_pdss[61];

    auto g_y_0_x_0_y_xz_0_0 = buffer_1010_pdss[62];

    auto g_y_0_x_0_y_yy_0_0 = buffer_1010_pdss[63];

    auto g_y_0_x_0_y_yz_0_0 = buffer_1010_pdss[64];

    auto g_y_0_x_0_y_zz_0_0 = buffer_1010_pdss[65];

    auto g_y_0_x_0_z_xx_0_0 = buffer_1010_pdss[66];

    auto g_y_0_x_0_z_xy_0_0 = buffer_1010_pdss[67];

    auto g_y_0_x_0_z_xz_0_0 = buffer_1010_pdss[68];

    auto g_y_0_x_0_z_yy_0_0 = buffer_1010_pdss[69];

    auto g_y_0_x_0_z_yz_0_0 = buffer_1010_pdss[70];

    auto g_y_0_x_0_z_zz_0_0 = buffer_1010_pdss[71];

    auto g_y_0_y_0_x_xx_0_0 = buffer_1010_pdss[72];

    auto g_y_0_y_0_x_xy_0_0 = buffer_1010_pdss[73];

    auto g_y_0_y_0_x_xz_0_0 = buffer_1010_pdss[74];

    auto g_y_0_y_0_x_yy_0_0 = buffer_1010_pdss[75];

    auto g_y_0_y_0_x_yz_0_0 = buffer_1010_pdss[76];

    auto g_y_0_y_0_x_zz_0_0 = buffer_1010_pdss[77];

    auto g_y_0_y_0_y_xx_0_0 = buffer_1010_pdss[78];

    auto g_y_0_y_0_y_xy_0_0 = buffer_1010_pdss[79];

    auto g_y_0_y_0_y_xz_0_0 = buffer_1010_pdss[80];

    auto g_y_0_y_0_y_yy_0_0 = buffer_1010_pdss[81];

    auto g_y_0_y_0_y_yz_0_0 = buffer_1010_pdss[82];

    auto g_y_0_y_0_y_zz_0_0 = buffer_1010_pdss[83];

    auto g_y_0_y_0_z_xx_0_0 = buffer_1010_pdss[84];

    auto g_y_0_y_0_z_xy_0_0 = buffer_1010_pdss[85];

    auto g_y_0_y_0_z_xz_0_0 = buffer_1010_pdss[86];

    auto g_y_0_y_0_z_yy_0_0 = buffer_1010_pdss[87];

    auto g_y_0_y_0_z_yz_0_0 = buffer_1010_pdss[88];

    auto g_y_0_y_0_z_zz_0_0 = buffer_1010_pdss[89];

    auto g_y_0_z_0_x_xx_0_0 = buffer_1010_pdss[90];

    auto g_y_0_z_0_x_xy_0_0 = buffer_1010_pdss[91];

    auto g_y_0_z_0_x_xz_0_0 = buffer_1010_pdss[92];

    auto g_y_0_z_0_x_yy_0_0 = buffer_1010_pdss[93];

    auto g_y_0_z_0_x_yz_0_0 = buffer_1010_pdss[94];

    auto g_y_0_z_0_x_zz_0_0 = buffer_1010_pdss[95];

    auto g_y_0_z_0_y_xx_0_0 = buffer_1010_pdss[96];

    auto g_y_0_z_0_y_xy_0_0 = buffer_1010_pdss[97];

    auto g_y_0_z_0_y_xz_0_0 = buffer_1010_pdss[98];

    auto g_y_0_z_0_y_yy_0_0 = buffer_1010_pdss[99];

    auto g_y_0_z_0_y_yz_0_0 = buffer_1010_pdss[100];

    auto g_y_0_z_0_y_zz_0_0 = buffer_1010_pdss[101];

    auto g_y_0_z_0_z_xx_0_0 = buffer_1010_pdss[102];

    auto g_y_0_z_0_z_xy_0_0 = buffer_1010_pdss[103];

    auto g_y_0_z_0_z_xz_0_0 = buffer_1010_pdss[104];

    auto g_y_0_z_0_z_yy_0_0 = buffer_1010_pdss[105];

    auto g_y_0_z_0_z_yz_0_0 = buffer_1010_pdss[106];

    auto g_y_0_z_0_z_zz_0_0 = buffer_1010_pdss[107];

    auto g_z_0_x_0_x_xx_0_0 = buffer_1010_pdss[108];

    auto g_z_0_x_0_x_xy_0_0 = buffer_1010_pdss[109];

    auto g_z_0_x_0_x_xz_0_0 = buffer_1010_pdss[110];

    auto g_z_0_x_0_x_yy_0_0 = buffer_1010_pdss[111];

    auto g_z_0_x_0_x_yz_0_0 = buffer_1010_pdss[112];

    auto g_z_0_x_0_x_zz_0_0 = buffer_1010_pdss[113];

    auto g_z_0_x_0_y_xx_0_0 = buffer_1010_pdss[114];

    auto g_z_0_x_0_y_xy_0_0 = buffer_1010_pdss[115];

    auto g_z_0_x_0_y_xz_0_0 = buffer_1010_pdss[116];

    auto g_z_0_x_0_y_yy_0_0 = buffer_1010_pdss[117];

    auto g_z_0_x_0_y_yz_0_0 = buffer_1010_pdss[118];

    auto g_z_0_x_0_y_zz_0_0 = buffer_1010_pdss[119];

    auto g_z_0_x_0_z_xx_0_0 = buffer_1010_pdss[120];

    auto g_z_0_x_0_z_xy_0_0 = buffer_1010_pdss[121];

    auto g_z_0_x_0_z_xz_0_0 = buffer_1010_pdss[122];

    auto g_z_0_x_0_z_yy_0_0 = buffer_1010_pdss[123];

    auto g_z_0_x_0_z_yz_0_0 = buffer_1010_pdss[124];

    auto g_z_0_x_0_z_zz_0_0 = buffer_1010_pdss[125];

    auto g_z_0_y_0_x_xx_0_0 = buffer_1010_pdss[126];

    auto g_z_0_y_0_x_xy_0_0 = buffer_1010_pdss[127];

    auto g_z_0_y_0_x_xz_0_0 = buffer_1010_pdss[128];

    auto g_z_0_y_0_x_yy_0_0 = buffer_1010_pdss[129];

    auto g_z_0_y_0_x_yz_0_0 = buffer_1010_pdss[130];

    auto g_z_0_y_0_x_zz_0_0 = buffer_1010_pdss[131];

    auto g_z_0_y_0_y_xx_0_0 = buffer_1010_pdss[132];

    auto g_z_0_y_0_y_xy_0_0 = buffer_1010_pdss[133];

    auto g_z_0_y_0_y_xz_0_0 = buffer_1010_pdss[134];

    auto g_z_0_y_0_y_yy_0_0 = buffer_1010_pdss[135];

    auto g_z_0_y_0_y_yz_0_0 = buffer_1010_pdss[136];

    auto g_z_0_y_0_y_zz_0_0 = buffer_1010_pdss[137];

    auto g_z_0_y_0_z_xx_0_0 = buffer_1010_pdss[138];

    auto g_z_0_y_0_z_xy_0_0 = buffer_1010_pdss[139];

    auto g_z_0_y_0_z_xz_0_0 = buffer_1010_pdss[140];

    auto g_z_0_y_0_z_yy_0_0 = buffer_1010_pdss[141];

    auto g_z_0_y_0_z_yz_0_0 = buffer_1010_pdss[142];

    auto g_z_0_y_0_z_zz_0_0 = buffer_1010_pdss[143];

    auto g_z_0_z_0_x_xx_0_0 = buffer_1010_pdss[144];

    auto g_z_0_z_0_x_xy_0_0 = buffer_1010_pdss[145];

    auto g_z_0_z_0_x_xz_0_0 = buffer_1010_pdss[146];

    auto g_z_0_z_0_x_yy_0_0 = buffer_1010_pdss[147];

    auto g_z_0_z_0_x_yz_0_0 = buffer_1010_pdss[148];

    auto g_z_0_z_0_x_zz_0_0 = buffer_1010_pdss[149];

    auto g_z_0_z_0_y_xx_0_0 = buffer_1010_pdss[150];

    auto g_z_0_z_0_y_xy_0_0 = buffer_1010_pdss[151];

    auto g_z_0_z_0_y_xz_0_0 = buffer_1010_pdss[152];

    auto g_z_0_z_0_y_yy_0_0 = buffer_1010_pdss[153];

    auto g_z_0_z_0_y_yz_0_0 = buffer_1010_pdss[154];

    auto g_z_0_z_0_y_zz_0_0 = buffer_1010_pdss[155];

    auto g_z_0_z_0_z_xx_0_0 = buffer_1010_pdss[156];

    auto g_z_0_z_0_z_xy_0_0 = buffer_1010_pdss[157];

    auto g_z_0_z_0_z_xz_0_0 = buffer_1010_pdss[158];

    auto g_z_0_z_0_z_yy_0_0 = buffer_1010_pdss[159];

    auto g_z_0_z_0_z_yz_0_0 = buffer_1010_pdss[160];

    auto g_z_0_z_0_z_zz_0_0 = buffer_1010_pdss[161];

    // integrals block (0-6)

    #pragma omp simd aligned(g_0_xx_x_0, g_0_xy_x_0, g_0_xz_x_0, g_0_yy_x_0, g_0_yz_x_0, g_0_zz_x_0, g_x_0_x_0_x_xx_0_0, g_x_0_x_0_x_xy_0_0, g_x_0_x_0_x_xz_0_0, g_x_0_x_0_x_yy_0_0, g_x_0_x_0_x_yz_0_0, g_x_0_x_0_x_zz_0_0, g_xx_xx_x_0, g_xx_xy_x_0, g_xx_xz_x_0, g_xx_yy_x_0, g_xx_yz_x_0, g_xx_zz_x_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_xx_0_0[i] = -2.0 * g_0_xx_x_0[i] * c_exps[i] + 4.0 * g_xx_xx_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xy_0_0[i] = -2.0 * g_0_xy_x_0[i] * c_exps[i] + 4.0 * g_xx_xy_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xz_0_0[i] = -2.0 * g_0_xz_x_0[i] * c_exps[i] + 4.0 * g_xx_xz_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_yy_0_0[i] = -2.0 * g_0_yy_x_0[i] * c_exps[i] + 4.0 * g_xx_yy_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_yz_0_0[i] = -2.0 * g_0_yz_x_0[i] * c_exps[i] + 4.0 * g_xx_yz_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_zz_0_0[i] = -2.0 * g_0_zz_x_0[i] * c_exps[i] + 4.0 * g_xx_zz_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_0_x_0_y_xx_0_0, g_x_0_x_0_y_xy_0_0, g_x_0_x_0_y_xz_0_0, g_x_0_x_0_y_yy_0_0, g_x_0_x_0_y_yz_0_0, g_x_0_x_0_y_zz_0_0, g_xy_xx_x_0, g_xy_xy_x_0, g_xy_xz_x_0, g_xy_yy_x_0, g_xy_yz_x_0, g_xy_zz_x_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_xx_0_0[i] = 4.0 * g_xy_xx_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xy_0_0[i] = 4.0 * g_xy_xy_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xz_0_0[i] = 4.0 * g_xy_xz_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_yy_0_0[i] = 4.0 * g_xy_yy_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_yz_0_0[i] = 4.0 * g_xy_yz_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_zz_0_0[i] = 4.0 * g_xy_zz_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_0_x_0_z_xx_0_0, g_x_0_x_0_z_xy_0_0, g_x_0_x_0_z_xz_0_0, g_x_0_x_0_z_yy_0_0, g_x_0_x_0_z_yz_0_0, g_x_0_x_0_z_zz_0_0, g_xz_xx_x_0, g_xz_xy_x_0, g_xz_xz_x_0, g_xz_yy_x_0, g_xz_yz_x_0, g_xz_zz_x_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_xx_0_0[i] = 4.0 * g_xz_xx_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xy_0_0[i] = 4.0 * g_xz_xy_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xz_0_0[i] = 4.0 * g_xz_xz_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_yy_0_0[i] = 4.0 * g_xz_yy_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_yz_0_0[i] = 4.0 * g_xz_yz_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_zz_0_0[i] = 4.0 * g_xz_zz_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_0_xx_y_0, g_0_xy_y_0, g_0_xz_y_0, g_0_yy_y_0, g_0_yz_y_0, g_0_zz_y_0, g_x_0_y_0_x_xx_0_0, g_x_0_y_0_x_xy_0_0, g_x_0_y_0_x_xz_0_0, g_x_0_y_0_x_yy_0_0, g_x_0_y_0_x_yz_0_0, g_x_0_y_0_x_zz_0_0, g_xx_xx_y_0, g_xx_xy_y_0, g_xx_xz_y_0, g_xx_yy_y_0, g_xx_yz_y_0, g_xx_zz_y_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_xx_0_0[i] = -2.0 * g_0_xx_y_0[i] * c_exps[i] + 4.0 * g_xx_xx_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xy_0_0[i] = -2.0 * g_0_xy_y_0[i] * c_exps[i] + 4.0 * g_xx_xy_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xz_0_0[i] = -2.0 * g_0_xz_y_0[i] * c_exps[i] + 4.0 * g_xx_xz_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_yy_0_0[i] = -2.0 * g_0_yy_y_0[i] * c_exps[i] + 4.0 * g_xx_yy_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_yz_0_0[i] = -2.0 * g_0_yz_y_0[i] * c_exps[i] + 4.0 * g_xx_yz_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_zz_0_0[i] = -2.0 * g_0_zz_y_0[i] * c_exps[i] + 4.0 * g_xx_zz_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_x_0_y_0_y_xx_0_0, g_x_0_y_0_y_xy_0_0, g_x_0_y_0_y_xz_0_0, g_x_0_y_0_y_yy_0_0, g_x_0_y_0_y_yz_0_0, g_x_0_y_0_y_zz_0_0, g_xy_xx_y_0, g_xy_xy_y_0, g_xy_xz_y_0, g_xy_yy_y_0, g_xy_yz_y_0, g_xy_zz_y_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_xx_0_0[i] = 4.0 * g_xy_xx_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xy_0_0[i] = 4.0 * g_xy_xy_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xz_0_0[i] = 4.0 * g_xy_xz_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_yy_0_0[i] = 4.0 * g_xy_yy_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_yz_0_0[i] = 4.0 * g_xy_yz_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_zz_0_0[i] = 4.0 * g_xy_zz_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_x_0_y_0_z_xx_0_0, g_x_0_y_0_z_xy_0_0, g_x_0_y_0_z_xz_0_0, g_x_0_y_0_z_yy_0_0, g_x_0_y_0_z_yz_0_0, g_x_0_y_0_z_zz_0_0, g_xz_xx_y_0, g_xz_xy_y_0, g_xz_xz_y_0, g_xz_yy_y_0, g_xz_yz_y_0, g_xz_zz_y_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_xx_0_0[i] = 4.0 * g_xz_xx_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xy_0_0[i] = 4.0 * g_xz_xy_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xz_0_0[i] = 4.0 * g_xz_xz_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_yy_0_0[i] = 4.0 * g_xz_yy_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_yz_0_0[i] = 4.0 * g_xz_yz_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_zz_0_0[i] = 4.0 * g_xz_zz_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_0_xx_z_0, g_0_xy_z_0, g_0_xz_z_0, g_0_yy_z_0, g_0_yz_z_0, g_0_zz_z_0, g_x_0_z_0_x_xx_0_0, g_x_0_z_0_x_xy_0_0, g_x_0_z_0_x_xz_0_0, g_x_0_z_0_x_yy_0_0, g_x_0_z_0_x_yz_0_0, g_x_0_z_0_x_zz_0_0, g_xx_xx_z_0, g_xx_xy_z_0, g_xx_xz_z_0, g_xx_yy_z_0, g_xx_yz_z_0, g_xx_zz_z_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_xx_0_0[i] = -2.0 * g_0_xx_z_0[i] * c_exps[i] + 4.0 * g_xx_xx_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xy_0_0[i] = -2.0 * g_0_xy_z_0[i] * c_exps[i] + 4.0 * g_xx_xy_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xz_0_0[i] = -2.0 * g_0_xz_z_0[i] * c_exps[i] + 4.0 * g_xx_xz_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_yy_0_0[i] = -2.0 * g_0_yy_z_0[i] * c_exps[i] + 4.0 * g_xx_yy_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_yz_0_0[i] = -2.0 * g_0_yz_z_0[i] * c_exps[i] + 4.0 * g_xx_yz_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_zz_0_0[i] = -2.0 * g_0_zz_z_0[i] * c_exps[i] + 4.0 * g_xx_zz_z_0[i] * a_exp * c_exps[i];
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_x_0_z_0_y_xx_0_0, g_x_0_z_0_y_xy_0_0, g_x_0_z_0_y_xz_0_0, g_x_0_z_0_y_yy_0_0, g_x_0_z_0_y_yz_0_0, g_x_0_z_0_y_zz_0_0, g_xy_xx_z_0, g_xy_xy_z_0, g_xy_xz_z_0, g_xy_yy_z_0, g_xy_yz_z_0, g_xy_zz_z_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_xx_0_0[i] = 4.0 * g_xy_xx_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xy_0_0[i] = 4.0 * g_xy_xy_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xz_0_0[i] = 4.0 * g_xy_xz_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_yy_0_0[i] = 4.0 * g_xy_yy_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_yz_0_0[i] = 4.0 * g_xy_yz_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_zz_0_0[i] = 4.0 * g_xy_zz_z_0[i] * a_exp * c_exps[i];
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_x_0_z_0_z_xx_0_0, g_x_0_z_0_z_xy_0_0, g_x_0_z_0_z_xz_0_0, g_x_0_z_0_z_yy_0_0, g_x_0_z_0_z_yz_0_0, g_x_0_z_0_z_zz_0_0, g_xz_xx_z_0, g_xz_xy_z_0, g_xz_xz_z_0, g_xz_yy_z_0, g_xz_yz_z_0, g_xz_zz_z_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_xx_0_0[i] = 4.0 * g_xz_xx_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xy_0_0[i] = 4.0 * g_xz_xy_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xz_0_0[i] = 4.0 * g_xz_xz_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_yy_0_0[i] = 4.0 * g_xz_yy_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_yz_0_0[i] = 4.0 * g_xz_yz_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_zz_0_0[i] = 4.0 * g_xz_zz_z_0[i] * a_exp * c_exps[i];
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_xy_xx_x_0, g_xy_xy_x_0, g_xy_xz_x_0, g_xy_yy_x_0, g_xy_yz_x_0, g_xy_zz_x_0, g_y_0_x_0_x_xx_0_0, g_y_0_x_0_x_xy_0_0, g_y_0_x_0_x_xz_0_0, g_y_0_x_0_x_yy_0_0, g_y_0_x_0_x_yz_0_0, g_y_0_x_0_x_zz_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_xx_0_0[i] = 4.0 * g_xy_xx_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xy_0_0[i] = 4.0 * g_xy_xy_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xz_0_0[i] = 4.0 * g_xy_xz_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_yy_0_0[i] = 4.0 * g_xy_yy_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_yz_0_0[i] = 4.0 * g_xy_yz_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_zz_0_0[i] = 4.0 * g_xy_zz_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_0_xx_x_0, g_0_xy_x_0, g_0_xz_x_0, g_0_yy_x_0, g_0_yz_x_0, g_0_zz_x_0, g_y_0_x_0_y_xx_0_0, g_y_0_x_0_y_xy_0_0, g_y_0_x_0_y_xz_0_0, g_y_0_x_0_y_yy_0_0, g_y_0_x_0_y_yz_0_0, g_y_0_x_0_y_zz_0_0, g_yy_xx_x_0, g_yy_xy_x_0, g_yy_xz_x_0, g_yy_yy_x_0, g_yy_yz_x_0, g_yy_zz_x_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_xx_0_0[i] = -2.0 * g_0_xx_x_0[i] * c_exps[i] + 4.0 * g_yy_xx_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xy_0_0[i] = -2.0 * g_0_xy_x_0[i] * c_exps[i] + 4.0 * g_yy_xy_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xz_0_0[i] = -2.0 * g_0_xz_x_0[i] * c_exps[i] + 4.0 * g_yy_xz_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_yy_0_0[i] = -2.0 * g_0_yy_x_0[i] * c_exps[i] + 4.0 * g_yy_yy_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_yz_0_0[i] = -2.0 * g_0_yz_x_0[i] * c_exps[i] + 4.0 * g_yy_yz_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_zz_0_0[i] = -2.0 * g_0_zz_x_0[i] * c_exps[i] + 4.0 * g_yy_zz_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_y_0_x_0_z_xx_0_0, g_y_0_x_0_z_xy_0_0, g_y_0_x_0_z_xz_0_0, g_y_0_x_0_z_yy_0_0, g_y_0_x_0_z_yz_0_0, g_y_0_x_0_z_zz_0_0, g_yz_xx_x_0, g_yz_xy_x_0, g_yz_xz_x_0, g_yz_yy_x_0, g_yz_yz_x_0, g_yz_zz_x_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_xx_0_0[i] = 4.0 * g_yz_xx_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xy_0_0[i] = 4.0 * g_yz_xy_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xz_0_0[i] = 4.0 * g_yz_xz_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_yy_0_0[i] = 4.0 * g_yz_yy_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_yz_0_0[i] = 4.0 * g_yz_yz_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_zz_0_0[i] = 4.0 * g_yz_zz_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_xy_xx_y_0, g_xy_xy_y_0, g_xy_xz_y_0, g_xy_yy_y_0, g_xy_yz_y_0, g_xy_zz_y_0, g_y_0_y_0_x_xx_0_0, g_y_0_y_0_x_xy_0_0, g_y_0_y_0_x_xz_0_0, g_y_0_y_0_x_yy_0_0, g_y_0_y_0_x_yz_0_0, g_y_0_y_0_x_zz_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_xx_0_0[i] = 4.0 * g_xy_xx_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xy_0_0[i] = 4.0 * g_xy_xy_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xz_0_0[i] = 4.0 * g_xy_xz_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_yy_0_0[i] = 4.0 * g_xy_yy_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_yz_0_0[i] = 4.0 * g_xy_yz_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_zz_0_0[i] = 4.0 * g_xy_zz_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_0_xx_y_0, g_0_xy_y_0, g_0_xz_y_0, g_0_yy_y_0, g_0_yz_y_0, g_0_zz_y_0, g_y_0_y_0_y_xx_0_0, g_y_0_y_0_y_xy_0_0, g_y_0_y_0_y_xz_0_0, g_y_0_y_0_y_yy_0_0, g_y_0_y_0_y_yz_0_0, g_y_0_y_0_y_zz_0_0, g_yy_xx_y_0, g_yy_xy_y_0, g_yy_xz_y_0, g_yy_yy_y_0, g_yy_yz_y_0, g_yy_zz_y_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_xx_0_0[i] = -2.0 * g_0_xx_y_0[i] * c_exps[i] + 4.0 * g_yy_xx_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xy_0_0[i] = -2.0 * g_0_xy_y_0[i] * c_exps[i] + 4.0 * g_yy_xy_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xz_0_0[i] = -2.0 * g_0_xz_y_0[i] * c_exps[i] + 4.0 * g_yy_xz_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_yy_0_0[i] = -2.0 * g_0_yy_y_0[i] * c_exps[i] + 4.0 * g_yy_yy_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_yz_0_0[i] = -2.0 * g_0_yz_y_0[i] * c_exps[i] + 4.0 * g_yy_yz_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_zz_0_0[i] = -2.0 * g_0_zz_y_0[i] * c_exps[i] + 4.0 * g_yy_zz_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_y_0_y_0_z_xx_0_0, g_y_0_y_0_z_xy_0_0, g_y_0_y_0_z_xz_0_0, g_y_0_y_0_z_yy_0_0, g_y_0_y_0_z_yz_0_0, g_y_0_y_0_z_zz_0_0, g_yz_xx_y_0, g_yz_xy_y_0, g_yz_xz_y_0, g_yz_yy_y_0, g_yz_yz_y_0, g_yz_zz_y_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_xx_0_0[i] = 4.0 * g_yz_xx_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xy_0_0[i] = 4.0 * g_yz_xy_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xz_0_0[i] = 4.0 * g_yz_xz_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_yy_0_0[i] = 4.0 * g_yz_yy_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_yz_0_0[i] = 4.0 * g_yz_yz_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_zz_0_0[i] = 4.0 * g_yz_zz_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_xy_xx_z_0, g_xy_xy_z_0, g_xy_xz_z_0, g_xy_yy_z_0, g_xy_yz_z_0, g_xy_zz_z_0, g_y_0_z_0_x_xx_0_0, g_y_0_z_0_x_xy_0_0, g_y_0_z_0_x_xz_0_0, g_y_0_z_0_x_yy_0_0, g_y_0_z_0_x_yz_0_0, g_y_0_z_0_x_zz_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_xx_0_0[i] = 4.0 * g_xy_xx_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xy_0_0[i] = 4.0 * g_xy_xy_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xz_0_0[i] = 4.0 * g_xy_xz_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_yy_0_0[i] = 4.0 * g_xy_yy_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_yz_0_0[i] = 4.0 * g_xy_yz_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_zz_0_0[i] = 4.0 * g_xy_zz_z_0[i] * a_exp * c_exps[i];
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_0_xx_z_0, g_0_xy_z_0, g_0_xz_z_0, g_0_yy_z_0, g_0_yz_z_0, g_0_zz_z_0, g_y_0_z_0_y_xx_0_0, g_y_0_z_0_y_xy_0_0, g_y_0_z_0_y_xz_0_0, g_y_0_z_0_y_yy_0_0, g_y_0_z_0_y_yz_0_0, g_y_0_z_0_y_zz_0_0, g_yy_xx_z_0, g_yy_xy_z_0, g_yy_xz_z_0, g_yy_yy_z_0, g_yy_yz_z_0, g_yy_zz_z_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_xx_0_0[i] = -2.0 * g_0_xx_z_0[i] * c_exps[i] + 4.0 * g_yy_xx_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xy_0_0[i] = -2.0 * g_0_xy_z_0[i] * c_exps[i] + 4.0 * g_yy_xy_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xz_0_0[i] = -2.0 * g_0_xz_z_0[i] * c_exps[i] + 4.0 * g_yy_xz_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_yy_0_0[i] = -2.0 * g_0_yy_z_0[i] * c_exps[i] + 4.0 * g_yy_yy_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_yz_0_0[i] = -2.0 * g_0_yz_z_0[i] * c_exps[i] + 4.0 * g_yy_yz_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_zz_0_0[i] = -2.0 * g_0_zz_z_0[i] * c_exps[i] + 4.0 * g_yy_zz_z_0[i] * a_exp * c_exps[i];
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_y_0_z_0_z_xx_0_0, g_y_0_z_0_z_xy_0_0, g_y_0_z_0_z_xz_0_0, g_y_0_z_0_z_yy_0_0, g_y_0_z_0_z_yz_0_0, g_y_0_z_0_z_zz_0_0, g_yz_xx_z_0, g_yz_xy_z_0, g_yz_xz_z_0, g_yz_yy_z_0, g_yz_yz_z_0, g_yz_zz_z_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_xx_0_0[i] = 4.0 * g_yz_xx_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xy_0_0[i] = 4.0 * g_yz_xy_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xz_0_0[i] = 4.0 * g_yz_xz_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_yy_0_0[i] = 4.0 * g_yz_yy_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_yz_0_0[i] = 4.0 * g_yz_yz_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_zz_0_0[i] = 4.0 * g_yz_zz_z_0[i] * a_exp * c_exps[i];
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_xz_xx_x_0, g_xz_xy_x_0, g_xz_xz_x_0, g_xz_yy_x_0, g_xz_yz_x_0, g_xz_zz_x_0, g_z_0_x_0_x_xx_0_0, g_z_0_x_0_x_xy_0_0, g_z_0_x_0_x_xz_0_0, g_z_0_x_0_x_yy_0_0, g_z_0_x_0_x_yz_0_0, g_z_0_x_0_x_zz_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_xx_0_0[i] = 4.0 * g_xz_xx_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xy_0_0[i] = 4.0 * g_xz_xy_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xz_0_0[i] = 4.0 * g_xz_xz_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_yy_0_0[i] = 4.0 * g_xz_yy_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_yz_0_0[i] = 4.0 * g_xz_yz_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_zz_0_0[i] = 4.0 * g_xz_zz_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_yz_xx_x_0, g_yz_xy_x_0, g_yz_xz_x_0, g_yz_yy_x_0, g_yz_yz_x_0, g_yz_zz_x_0, g_z_0_x_0_y_xx_0_0, g_z_0_x_0_y_xy_0_0, g_z_0_x_0_y_xz_0_0, g_z_0_x_0_y_yy_0_0, g_z_0_x_0_y_yz_0_0, g_z_0_x_0_y_zz_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_xx_0_0[i] = 4.0 * g_yz_xx_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xy_0_0[i] = 4.0 * g_yz_xy_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xz_0_0[i] = 4.0 * g_yz_xz_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_yy_0_0[i] = 4.0 * g_yz_yy_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_yz_0_0[i] = 4.0 * g_yz_yz_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_zz_0_0[i] = 4.0 * g_yz_zz_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_0_xx_x_0, g_0_xy_x_0, g_0_xz_x_0, g_0_yy_x_0, g_0_yz_x_0, g_0_zz_x_0, g_z_0_x_0_z_xx_0_0, g_z_0_x_0_z_xy_0_0, g_z_0_x_0_z_xz_0_0, g_z_0_x_0_z_yy_0_0, g_z_0_x_0_z_yz_0_0, g_z_0_x_0_z_zz_0_0, g_zz_xx_x_0, g_zz_xy_x_0, g_zz_xz_x_0, g_zz_yy_x_0, g_zz_yz_x_0, g_zz_zz_x_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_xx_0_0[i] = -2.0 * g_0_xx_x_0[i] * c_exps[i] + 4.0 * g_zz_xx_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xy_0_0[i] = -2.0 * g_0_xy_x_0[i] * c_exps[i] + 4.0 * g_zz_xy_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xz_0_0[i] = -2.0 * g_0_xz_x_0[i] * c_exps[i] + 4.0 * g_zz_xz_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_yy_0_0[i] = -2.0 * g_0_yy_x_0[i] * c_exps[i] + 4.0 * g_zz_yy_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_yz_0_0[i] = -2.0 * g_0_yz_x_0[i] * c_exps[i] + 4.0 * g_zz_yz_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_zz_0_0[i] = -2.0 * g_0_zz_x_0[i] * c_exps[i] + 4.0 * g_zz_zz_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_xz_xx_y_0, g_xz_xy_y_0, g_xz_xz_y_0, g_xz_yy_y_0, g_xz_yz_y_0, g_xz_zz_y_0, g_z_0_y_0_x_xx_0_0, g_z_0_y_0_x_xy_0_0, g_z_0_y_0_x_xz_0_0, g_z_0_y_0_x_yy_0_0, g_z_0_y_0_x_yz_0_0, g_z_0_y_0_x_zz_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_xx_0_0[i] = 4.0 * g_xz_xx_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xy_0_0[i] = 4.0 * g_xz_xy_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xz_0_0[i] = 4.0 * g_xz_xz_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_yy_0_0[i] = 4.0 * g_xz_yy_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_yz_0_0[i] = 4.0 * g_xz_yz_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_zz_0_0[i] = 4.0 * g_xz_zz_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_yz_xx_y_0, g_yz_xy_y_0, g_yz_xz_y_0, g_yz_yy_y_0, g_yz_yz_y_0, g_yz_zz_y_0, g_z_0_y_0_y_xx_0_0, g_z_0_y_0_y_xy_0_0, g_z_0_y_0_y_xz_0_0, g_z_0_y_0_y_yy_0_0, g_z_0_y_0_y_yz_0_0, g_z_0_y_0_y_zz_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_xx_0_0[i] = 4.0 * g_yz_xx_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xy_0_0[i] = 4.0 * g_yz_xy_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xz_0_0[i] = 4.0 * g_yz_xz_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_yy_0_0[i] = 4.0 * g_yz_yy_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_yz_0_0[i] = 4.0 * g_yz_yz_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_zz_0_0[i] = 4.0 * g_yz_zz_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_0_xx_y_0, g_0_xy_y_0, g_0_xz_y_0, g_0_yy_y_0, g_0_yz_y_0, g_0_zz_y_0, g_z_0_y_0_z_xx_0_0, g_z_0_y_0_z_xy_0_0, g_z_0_y_0_z_xz_0_0, g_z_0_y_0_z_yy_0_0, g_z_0_y_0_z_yz_0_0, g_z_0_y_0_z_zz_0_0, g_zz_xx_y_0, g_zz_xy_y_0, g_zz_xz_y_0, g_zz_yy_y_0, g_zz_yz_y_0, g_zz_zz_y_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_xx_0_0[i] = -2.0 * g_0_xx_y_0[i] * c_exps[i] + 4.0 * g_zz_xx_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xy_0_0[i] = -2.0 * g_0_xy_y_0[i] * c_exps[i] + 4.0 * g_zz_xy_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xz_0_0[i] = -2.0 * g_0_xz_y_0[i] * c_exps[i] + 4.0 * g_zz_xz_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_yy_0_0[i] = -2.0 * g_0_yy_y_0[i] * c_exps[i] + 4.0 * g_zz_yy_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_yz_0_0[i] = -2.0 * g_0_yz_y_0[i] * c_exps[i] + 4.0 * g_zz_yz_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_zz_0_0[i] = -2.0 * g_0_zz_y_0[i] * c_exps[i] + 4.0 * g_zz_zz_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_xz_xx_z_0, g_xz_xy_z_0, g_xz_xz_z_0, g_xz_yy_z_0, g_xz_yz_z_0, g_xz_zz_z_0, g_z_0_z_0_x_xx_0_0, g_z_0_z_0_x_xy_0_0, g_z_0_z_0_x_xz_0_0, g_z_0_z_0_x_yy_0_0, g_z_0_z_0_x_yz_0_0, g_z_0_z_0_x_zz_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_xx_0_0[i] = 4.0 * g_xz_xx_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xy_0_0[i] = 4.0 * g_xz_xy_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xz_0_0[i] = 4.0 * g_xz_xz_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_yy_0_0[i] = 4.0 * g_xz_yy_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_yz_0_0[i] = 4.0 * g_xz_yz_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_zz_0_0[i] = 4.0 * g_xz_zz_z_0[i] * a_exp * c_exps[i];
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_yz_xx_z_0, g_yz_xy_z_0, g_yz_xz_z_0, g_yz_yy_z_0, g_yz_yz_z_0, g_yz_zz_z_0, g_z_0_z_0_y_xx_0_0, g_z_0_z_0_y_xy_0_0, g_z_0_z_0_y_xz_0_0, g_z_0_z_0_y_yy_0_0, g_z_0_z_0_y_yz_0_0, g_z_0_z_0_y_zz_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_xx_0_0[i] = 4.0 * g_yz_xx_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xy_0_0[i] = 4.0 * g_yz_xy_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xz_0_0[i] = 4.0 * g_yz_xz_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_yy_0_0[i] = 4.0 * g_yz_yy_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_yz_0_0[i] = 4.0 * g_yz_yz_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_zz_0_0[i] = 4.0 * g_yz_zz_z_0[i] * a_exp * c_exps[i];
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_0_xx_z_0, g_0_xy_z_0, g_0_xz_z_0, g_0_yy_z_0, g_0_yz_z_0, g_0_zz_z_0, g_z_0_z_0_z_xx_0_0, g_z_0_z_0_z_xy_0_0, g_z_0_z_0_z_xz_0_0, g_z_0_z_0_z_yy_0_0, g_z_0_z_0_z_yz_0_0, g_z_0_z_0_z_zz_0_0, g_zz_xx_z_0, g_zz_xy_z_0, g_zz_xz_z_0, g_zz_yy_z_0, g_zz_yz_z_0, g_zz_zz_z_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_xx_0_0[i] = -2.0 * g_0_xx_z_0[i] * c_exps[i] + 4.0 * g_zz_xx_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xy_0_0[i] = -2.0 * g_0_xy_z_0[i] * c_exps[i] + 4.0 * g_zz_xy_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xz_0_0[i] = -2.0 * g_0_xz_z_0[i] * c_exps[i] + 4.0 * g_zz_xz_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_yy_0_0[i] = -2.0 * g_0_yy_z_0[i] * c_exps[i] + 4.0 * g_zz_yy_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_yz_0_0[i] = -2.0 * g_0_yz_z_0[i] * c_exps[i] + 4.0 * g_zz_yz_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_zz_0_0[i] = -2.0 * g_0_zz_z_0[i] * c_exps[i] + 4.0 * g_zz_zz_z_0[i] * a_exp * c_exps[i];
    }
}

} // t4c_geom namespace


#include "GeomDeriv1100OfScalarForSSPD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1100_sspd_0(CSimdArray<double>& buffer_1100_sspd,
                     const CSimdArray<double>& buffer_pppd,
                     const double a_exp,
                     const double b_exp) -> void
{
    const auto ndims = buffer_1100_sspd.number_of_columns();

    /// Set up components of auxilary buffer : buffer_pppd

    auto g_x_x_x_xx = buffer_pppd[0];

    auto g_x_x_x_xy = buffer_pppd[1];

    auto g_x_x_x_xz = buffer_pppd[2];

    auto g_x_x_x_yy = buffer_pppd[3];

    auto g_x_x_x_yz = buffer_pppd[4];

    auto g_x_x_x_zz = buffer_pppd[5];

    auto g_x_x_y_xx = buffer_pppd[6];

    auto g_x_x_y_xy = buffer_pppd[7];

    auto g_x_x_y_xz = buffer_pppd[8];

    auto g_x_x_y_yy = buffer_pppd[9];

    auto g_x_x_y_yz = buffer_pppd[10];

    auto g_x_x_y_zz = buffer_pppd[11];

    auto g_x_x_z_xx = buffer_pppd[12];

    auto g_x_x_z_xy = buffer_pppd[13];

    auto g_x_x_z_xz = buffer_pppd[14];

    auto g_x_x_z_yy = buffer_pppd[15];

    auto g_x_x_z_yz = buffer_pppd[16];

    auto g_x_x_z_zz = buffer_pppd[17];

    auto g_x_y_x_xx = buffer_pppd[18];

    auto g_x_y_x_xy = buffer_pppd[19];

    auto g_x_y_x_xz = buffer_pppd[20];

    auto g_x_y_x_yy = buffer_pppd[21];

    auto g_x_y_x_yz = buffer_pppd[22];

    auto g_x_y_x_zz = buffer_pppd[23];

    auto g_x_y_y_xx = buffer_pppd[24];

    auto g_x_y_y_xy = buffer_pppd[25];

    auto g_x_y_y_xz = buffer_pppd[26];

    auto g_x_y_y_yy = buffer_pppd[27];

    auto g_x_y_y_yz = buffer_pppd[28];

    auto g_x_y_y_zz = buffer_pppd[29];

    auto g_x_y_z_xx = buffer_pppd[30];

    auto g_x_y_z_xy = buffer_pppd[31];

    auto g_x_y_z_xz = buffer_pppd[32];

    auto g_x_y_z_yy = buffer_pppd[33];

    auto g_x_y_z_yz = buffer_pppd[34];

    auto g_x_y_z_zz = buffer_pppd[35];

    auto g_x_z_x_xx = buffer_pppd[36];

    auto g_x_z_x_xy = buffer_pppd[37];

    auto g_x_z_x_xz = buffer_pppd[38];

    auto g_x_z_x_yy = buffer_pppd[39];

    auto g_x_z_x_yz = buffer_pppd[40];

    auto g_x_z_x_zz = buffer_pppd[41];

    auto g_x_z_y_xx = buffer_pppd[42];

    auto g_x_z_y_xy = buffer_pppd[43];

    auto g_x_z_y_xz = buffer_pppd[44];

    auto g_x_z_y_yy = buffer_pppd[45];

    auto g_x_z_y_yz = buffer_pppd[46];

    auto g_x_z_y_zz = buffer_pppd[47];

    auto g_x_z_z_xx = buffer_pppd[48];

    auto g_x_z_z_xy = buffer_pppd[49];

    auto g_x_z_z_xz = buffer_pppd[50];

    auto g_x_z_z_yy = buffer_pppd[51];

    auto g_x_z_z_yz = buffer_pppd[52];

    auto g_x_z_z_zz = buffer_pppd[53];

    auto g_y_x_x_xx = buffer_pppd[54];

    auto g_y_x_x_xy = buffer_pppd[55];

    auto g_y_x_x_xz = buffer_pppd[56];

    auto g_y_x_x_yy = buffer_pppd[57];

    auto g_y_x_x_yz = buffer_pppd[58];

    auto g_y_x_x_zz = buffer_pppd[59];

    auto g_y_x_y_xx = buffer_pppd[60];

    auto g_y_x_y_xy = buffer_pppd[61];

    auto g_y_x_y_xz = buffer_pppd[62];

    auto g_y_x_y_yy = buffer_pppd[63];

    auto g_y_x_y_yz = buffer_pppd[64];

    auto g_y_x_y_zz = buffer_pppd[65];

    auto g_y_x_z_xx = buffer_pppd[66];

    auto g_y_x_z_xy = buffer_pppd[67];

    auto g_y_x_z_xz = buffer_pppd[68];

    auto g_y_x_z_yy = buffer_pppd[69];

    auto g_y_x_z_yz = buffer_pppd[70];

    auto g_y_x_z_zz = buffer_pppd[71];

    auto g_y_y_x_xx = buffer_pppd[72];

    auto g_y_y_x_xy = buffer_pppd[73];

    auto g_y_y_x_xz = buffer_pppd[74];

    auto g_y_y_x_yy = buffer_pppd[75];

    auto g_y_y_x_yz = buffer_pppd[76];

    auto g_y_y_x_zz = buffer_pppd[77];

    auto g_y_y_y_xx = buffer_pppd[78];

    auto g_y_y_y_xy = buffer_pppd[79];

    auto g_y_y_y_xz = buffer_pppd[80];

    auto g_y_y_y_yy = buffer_pppd[81];

    auto g_y_y_y_yz = buffer_pppd[82];

    auto g_y_y_y_zz = buffer_pppd[83];

    auto g_y_y_z_xx = buffer_pppd[84];

    auto g_y_y_z_xy = buffer_pppd[85];

    auto g_y_y_z_xz = buffer_pppd[86];

    auto g_y_y_z_yy = buffer_pppd[87];

    auto g_y_y_z_yz = buffer_pppd[88];

    auto g_y_y_z_zz = buffer_pppd[89];

    auto g_y_z_x_xx = buffer_pppd[90];

    auto g_y_z_x_xy = buffer_pppd[91];

    auto g_y_z_x_xz = buffer_pppd[92];

    auto g_y_z_x_yy = buffer_pppd[93];

    auto g_y_z_x_yz = buffer_pppd[94];

    auto g_y_z_x_zz = buffer_pppd[95];

    auto g_y_z_y_xx = buffer_pppd[96];

    auto g_y_z_y_xy = buffer_pppd[97];

    auto g_y_z_y_xz = buffer_pppd[98];

    auto g_y_z_y_yy = buffer_pppd[99];

    auto g_y_z_y_yz = buffer_pppd[100];

    auto g_y_z_y_zz = buffer_pppd[101];

    auto g_y_z_z_xx = buffer_pppd[102];

    auto g_y_z_z_xy = buffer_pppd[103];

    auto g_y_z_z_xz = buffer_pppd[104];

    auto g_y_z_z_yy = buffer_pppd[105];

    auto g_y_z_z_yz = buffer_pppd[106];

    auto g_y_z_z_zz = buffer_pppd[107];

    auto g_z_x_x_xx = buffer_pppd[108];

    auto g_z_x_x_xy = buffer_pppd[109];

    auto g_z_x_x_xz = buffer_pppd[110];

    auto g_z_x_x_yy = buffer_pppd[111];

    auto g_z_x_x_yz = buffer_pppd[112];

    auto g_z_x_x_zz = buffer_pppd[113];

    auto g_z_x_y_xx = buffer_pppd[114];

    auto g_z_x_y_xy = buffer_pppd[115];

    auto g_z_x_y_xz = buffer_pppd[116];

    auto g_z_x_y_yy = buffer_pppd[117];

    auto g_z_x_y_yz = buffer_pppd[118];

    auto g_z_x_y_zz = buffer_pppd[119];

    auto g_z_x_z_xx = buffer_pppd[120];

    auto g_z_x_z_xy = buffer_pppd[121];

    auto g_z_x_z_xz = buffer_pppd[122];

    auto g_z_x_z_yy = buffer_pppd[123];

    auto g_z_x_z_yz = buffer_pppd[124];

    auto g_z_x_z_zz = buffer_pppd[125];

    auto g_z_y_x_xx = buffer_pppd[126];

    auto g_z_y_x_xy = buffer_pppd[127];

    auto g_z_y_x_xz = buffer_pppd[128];

    auto g_z_y_x_yy = buffer_pppd[129];

    auto g_z_y_x_yz = buffer_pppd[130];

    auto g_z_y_x_zz = buffer_pppd[131];

    auto g_z_y_y_xx = buffer_pppd[132];

    auto g_z_y_y_xy = buffer_pppd[133];

    auto g_z_y_y_xz = buffer_pppd[134];

    auto g_z_y_y_yy = buffer_pppd[135];

    auto g_z_y_y_yz = buffer_pppd[136];

    auto g_z_y_y_zz = buffer_pppd[137];

    auto g_z_y_z_xx = buffer_pppd[138];

    auto g_z_y_z_xy = buffer_pppd[139];

    auto g_z_y_z_xz = buffer_pppd[140];

    auto g_z_y_z_yy = buffer_pppd[141];

    auto g_z_y_z_yz = buffer_pppd[142];

    auto g_z_y_z_zz = buffer_pppd[143];

    auto g_z_z_x_xx = buffer_pppd[144];

    auto g_z_z_x_xy = buffer_pppd[145];

    auto g_z_z_x_xz = buffer_pppd[146];

    auto g_z_z_x_yy = buffer_pppd[147];

    auto g_z_z_x_yz = buffer_pppd[148];

    auto g_z_z_x_zz = buffer_pppd[149];

    auto g_z_z_y_xx = buffer_pppd[150];

    auto g_z_z_y_xy = buffer_pppd[151];

    auto g_z_z_y_xz = buffer_pppd[152];

    auto g_z_z_y_yy = buffer_pppd[153];

    auto g_z_z_y_yz = buffer_pppd[154];

    auto g_z_z_y_zz = buffer_pppd[155];

    auto g_z_z_z_xx = buffer_pppd[156];

    auto g_z_z_z_xy = buffer_pppd[157];

    auto g_z_z_z_xz = buffer_pppd[158];

    auto g_z_z_z_yy = buffer_pppd[159];

    auto g_z_z_z_yz = buffer_pppd[160];

    auto g_z_z_z_zz = buffer_pppd[161];

    /// Set up components of integrals buffer : buffer_1100_sspd

    auto g_x_x_0_0_0_0_x_xx = buffer_1100_sspd[0];

    auto g_x_x_0_0_0_0_x_xy = buffer_1100_sspd[1];

    auto g_x_x_0_0_0_0_x_xz = buffer_1100_sspd[2];

    auto g_x_x_0_0_0_0_x_yy = buffer_1100_sspd[3];

    auto g_x_x_0_0_0_0_x_yz = buffer_1100_sspd[4];

    auto g_x_x_0_0_0_0_x_zz = buffer_1100_sspd[5];

    auto g_x_x_0_0_0_0_y_xx = buffer_1100_sspd[6];

    auto g_x_x_0_0_0_0_y_xy = buffer_1100_sspd[7];

    auto g_x_x_0_0_0_0_y_xz = buffer_1100_sspd[8];

    auto g_x_x_0_0_0_0_y_yy = buffer_1100_sspd[9];

    auto g_x_x_0_0_0_0_y_yz = buffer_1100_sspd[10];

    auto g_x_x_0_0_0_0_y_zz = buffer_1100_sspd[11];

    auto g_x_x_0_0_0_0_z_xx = buffer_1100_sspd[12];

    auto g_x_x_0_0_0_0_z_xy = buffer_1100_sspd[13];

    auto g_x_x_0_0_0_0_z_xz = buffer_1100_sspd[14];

    auto g_x_x_0_0_0_0_z_yy = buffer_1100_sspd[15];

    auto g_x_x_0_0_0_0_z_yz = buffer_1100_sspd[16];

    auto g_x_x_0_0_0_0_z_zz = buffer_1100_sspd[17];

    auto g_x_y_0_0_0_0_x_xx = buffer_1100_sspd[18];

    auto g_x_y_0_0_0_0_x_xy = buffer_1100_sspd[19];

    auto g_x_y_0_0_0_0_x_xz = buffer_1100_sspd[20];

    auto g_x_y_0_0_0_0_x_yy = buffer_1100_sspd[21];

    auto g_x_y_0_0_0_0_x_yz = buffer_1100_sspd[22];

    auto g_x_y_0_0_0_0_x_zz = buffer_1100_sspd[23];

    auto g_x_y_0_0_0_0_y_xx = buffer_1100_sspd[24];

    auto g_x_y_0_0_0_0_y_xy = buffer_1100_sspd[25];

    auto g_x_y_0_0_0_0_y_xz = buffer_1100_sspd[26];

    auto g_x_y_0_0_0_0_y_yy = buffer_1100_sspd[27];

    auto g_x_y_0_0_0_0_y_yz = buffer_1100_sspd[28];

    auto g_x_y_0_0_0_0_y_zz = buffer_1100_sspd[29];

    auto g_x_y_0_0_0_0_z_xx = buffer_1100_sspd[30];

    auto g_x_y_0_0_0_0_z_xy = buffer_1100_sspd[31];

    auto g_x_y_0_0_0_0_z_xz = buffer_1100_sspd[32];

    auto g_x_y_0_0_0_0_z_yy = buffer_1100_sspd[33];

    auto g_x_y_0_0_0_0_z_yz = buffer_1100_sspd[34];

    auto g_x_y_0_0_0_0_z_zz = buffer_1100_sspd[35];

    auto g_x_z_0_0_0_0_x_xx = buffer_1100_sspd[36];

    auto g_x_z_0_0_0_0_x_xy = buffer_1100_sspd[37];

    auto g_x_z_0_0_0_0_x_xz = buffer_1100_sspd[38];

    auto g_x_z_0_0_0_0_x_yy = buffer_1100_sspd[39];

    auto g_x_z_0_0_0_0_x_yz = buffer_1100_sspd[40];

    auto g_x_z_0_0_0_0_x_zz = buffer_1100_sspd[41];

    auto g_x_z_0_0_0_0_y_xx = buffer_1100_sspd[42];

    auto g_x_z_0_0_0_0_y_xy = buffer_1100_sspd[43];

    auto g_x_z_0_0_0_0_y_xz = buffer_1100_sspd[44];

    auto g_x_z_0_0_0_0_y_yy = buffer_1100_sspd[45];

    auto g_x_z_0_0_0_0_y_yz = buffer_1100_sspd[46];

    auto g_x_z_0_0_0_0_y_zz = buffer_1100_sspd[47];

    auto g_x_z_0_0_0_0_z_xx = buffer_1100_sspd[48];

    auto g_x_z_0_0_0_0_z_xy = buffer_1100_sspd[49];

    auto g_x_z_0_0_0_0_z_xz = buffer_1100_sspd[50];

    auto g_x_z_0_0_0_0_z_yy = buffer_1100_sspd[51];

    auto g_x_z_0_0_0_0_z_yz = buffer_1100_sspd[52];

    auto g_x_z_0_0_0_0_z_zz = buffer_1100_sspd[53];

    auto g_y_x_0_0_0_0_x_xx = buffer_1100_sspd[54];

    auto g_y_x_0_0_0_0_x_xy = buffer_1100_sspd[55];

    auto g_y_x_0_0_0_0_x_xz = buffer_1100_sspd[56];

    auto g_y_x_0_0_0_0_x_yy = buffer_1100_sspd[57];

    auto g_y_x_0_0_0_0_x_yz = buffer_1100_sspd[58];

    auto g_y_x_0_0_0_0_x_zz = buffer_1100_sspd[59];

    auto g_y_x_0_0_0_0_y_xx = buffer_1100_sspd[60];

    auto g_y_x_0_0_0_0_y_xy = buffer_1100_sspd[61];

    auto g_y_x_0_0_0_0_y_xz = buffer_1100_sspd[62];

    auto g_y_x_0_0_0_0_y_yy = buffer_1100_sspd[63];

    auto g_y_x_0_0_0_0_y_yz = buffer_1100_sspd[64];

    auto g_y_x_0_0_0_0_y_zz = buffer_1100_sspd[65];

    auto g_y_x_0_0_0_0_z_xx = buffer_1100_sspd[66];

    auto g_y_x_0_0_0_0_z_xy = buffer_1100_sspd[67];

    auto g_y_x_0_0_0_0_z_xz = buffer_1100_sspd[68];

    auto g_y_x_0_0_0_0_z_yy = buffer_1100_sspd[69];

    auto g_y_x_0_0_0_0_z_yz = buffer_1100_sspd[70];

    auto g_y_x_0_0_0_0_z_zz = buffer_1100_sspd[71];

    auto g_y_y_0_0_0_0_x_xx = buffer_1100_sspd[72];

    auto g_y_y_0_0_0_0_x_xy = buffer_1100_sspd[73];

    auto g_y_y_0_0_0_0_x_xz = buffer_1100_sspd[74];

    auto g_y_y_0_0_0_0_x_yy = buffer_1100_sspd[75];

    auto g_y_y_0_0_0_0_x_yz = buffer_1100_sspd[76];

    auto g_y_y_0_0_0_0_x_zz = buffer_1100_sspd[77];

    auto g_y_y_0_0_0_0_y_xx = buffer_1100_sspd[78];

    auto g_y_y_0_0_0_0_y_xy = buffer_1100_sspd[79];

    auto g_y_y_0_0_0_0_y_xz = buffer_1100_sspd[80];

    auto g_y_y_0_0_0_0_y_yy = buffer_1100_sspd[81];

    auto g_y_y_0_0_0_0_y_yz = buffer_1100_sspd[82];

    auto g_y_y_0_0_0_0_y_zz = buffer_1100_sspd[83];

    auto g_y_y_0_0_0_0_z_xx = buffer_1100_sspd[84];

    auto g_y_y_0_0_0_0_z_xy = buffer_1100_sspd[85];

    auto g_y_y_0_0_0_0_z_xz = buffer_1100_sspd[86];

    auto g_y_y_0_0_0_0_z_yy = buffer_1100_sspd[87];

    auto g_y_y_0_0_0_0_z_yz = buffer_1100_sspd[88];

    auto g_y_y_0_0_0_0_z_zz = buffer_1100_sspd[89];

    auto g_y_z_0_0_0_0_x_xx = buffer_1100_sspd[90];

    auto g_y_z_0_0_0_0_x_xy = buffer_1100_sspd[91];

    auto g_y_z_0_0_0_0_x_xz = buffer_1100_sspd[92];

    auto g_y_z_0_0_0_0_x_yy = buffer_1100_sspd[93];

    auto g_y_z_0_0_0_0_x_yz = buffer_1100_sspd[94];

    auto g_y_z_0_0_0_0_x_zz = buffer_1100_sspd[95];

    auto g_y_z_0_0_0_0_y_xx = buffer_1100_sspd[96];

    auto g_y_z_0_0_0_0_y_xy = buffer_1100_sspd[97];

    auto g_y_z_0_0_0_0_y_xz = buffer_1100_sspd[98];

    auto g_y_z_0_0_0_0_y_yy = buffer_1100_sspd[99];

    auto g_y_z_0_0_0_0_y_yz = buffer_1100_sspd[100];

    auto g_y_z_0_0_0_0_y_zz = buffer_1100_sspd[101];

    auto g_y_z_0_0_0_0_z_xx = buffer_1100_sspd[102];

    auto g_y_z_0_0_0_0_z_xy = buffer_1100_sspd[103];

    auto g_y_z_0_0_0_0_z_xz = buffer_1100_sspd[104];

    auto g_y_z_0_0_0_0_z_yy = buffer_1100_sspd[105];

    auto g_y_z_0_0_0_0_z_yz = buffer_1100_sspd[106];

    auto g_y_z_0_0_0_0_z_zz = buffer_1100_sspd[107];

    auto g_z_x_0_0_0_0_x_xx = buffer_1100_sspd[108];

    auto g_z_x_0_0_0_0_x_xy = buffer_1100_sspd[109];

    auto g_z_x_0_0_0_0_x_xz = buffer_1100_sspd[110];

    auto g_z_x_0_0_0_0_x_yy = buffer_1100_sspd[111];

    auto g_z_x_0_0_0_0_x_yz = buffer_1100_sspd[112];

    auto g_z_x_0_0_0_0_x_zz = buffer_1100_sspd[113];

    auto g_z_x_0_0_0_0_y_xx = buffer_1100_sspd[114];

    auto g_z_x_0_0_0_0_y_xy = buffer_1100_sspd[115];

    auto g_z_x_0_0_0_0_y_xz = buffer_1100_sspd[116];

    auto g_z_x_0_0_0_0_y_yy = buffer_1100_sspd[117];

    auto g_z_x_0_0_0_0_y_yz = buffer_1100_sspd[118];

    auto g_z_x_0_0_0_0_y_zz = buffer_1100_sspd[119];

    auto g_z_x_0_0_0_0_z_xx = buffer_1100_sspd[120];

    auto g_z_x_0_0_0_0_z_xy = buffer_1100_sspd[121];

    auto g_z_x_0_0_0_0_z_xz = buffer_1100_sspd[122];

    auto g_z_x_0_0_0_0_z_yy = buffer_1100_sspd[123];

    auto g_z_x_0_0_0_0_z_yz = buffer_1100_sspd[124];

    auto g_z_x_0_0_0_0_z_zz = buffer_1100_sspd[125];

    auto g_z_y_0_0_0_0_x_xx = buffer_1100_sspd[126];

    auto g_z_y_0_0_0_0_x_xy = buffer_1100_sspd[127];

    auto g_z_y_0_0_0_0_x_xz = buffer_1100_sspd[128];

    auto g_z_y_0_0_0_0_x_yy = buffer_1100_sspd[129];

    auto g_z_y_0_0_0_0_x_yz = buffer_1100_sspd[130];

    auto g_z_y_0_0_0_0_x_zz = buffer_1100_sspd[131];

    auto g_z_y_0_0_0_0_y_xx = buffer_1100_sspd[132];

    auto g_z_y_0_0_0_0_y_xy = buffer_1100_sspd[133];

    auto g_z_y_0_0_0_0_y_xz = buffer_1100_sspd[134];

    auto g_z_y_0_0_0_0_y_yy = buffer_1100_sspd[135];

    auto g_z_y_0_0_0_0_y_yz = buffer_1100_sspd[136];

    auto g_z_y_0_0_0_0_y_zz = buffer_1100_sspd[137];

    auto g_z_y_0_0_0_0_z_xx = buffer_1100_sspd[138];

    auto g_z_y_0_0_0_0_z_xy = buffer_1100_sspd[139];

    auto g_z_y_0_0_0_0_z_xz = buffer_1100_sspd[140];

    auto g_z_y_0_0_0_0_z_yy = buffer_1100_sspd[141];

    auto g_z_y_0_0_0_0_z_yz = buffer_1100_sspd[142];

    auto g_z_y_0_0_0_0_z_zz = buffer_1100_sspd[143];

    auto g_z_z_0_0_0_0_x_xx = buffer_1100_sspd[144];

    auto g_z_z_0_0_0_0_x_xy = buffer_1100_sspd[145];

    auto g_z_z_0_0_0_0_x_xz = buffer_1100_sspd[146];

    auto g_z_z_0_0_0_0_x_yy = buffer_1100_sspd[147];

    auto g_z_z_0_0_0_0_x_yz = buffer_1100_sspd[148];

    auto g_z_z_0_0_0_0_x_zz = buffer_1100_sspd[149];

    auto g_z_z_0_0_0_0_y_xx = buffer_1100_sspd[150];

    auto g_z_z_0_0_0_0_y_xy = buffer_1100_sspd[151];

    auto g_z_z_0_0_0_0_y_xz = buffer_1100_sspd[152];

    auto g_z_z_0_0_0_0_y_yy = buffer_1100_sspd[153];

    auto g_z_z_0_0_0_0_y_yz = buffer_1100_sspd[154];

    auto g_z_z_0_0_0_0_y_zz = buffer_1100_sspd[155];

    auto g_z_z_0_0_0_0_z_xx = buffer_1100_sspd[156];

    auto g_z_z_0_0_0_0_z_xy = buffer_1100_sspd[157];

    auto g_z_z_0_0_0_0_z_xz = buffer_1100_sspd[158];

    auto g_z_z_0_0_0_0_z_yy = buffer_1100_sspd[159];

    auto g_z_z_0_0_0_0_z_yz = buffer_1100_sspd[160];

    auto g_z_z_0_0_0_0_z_zz = buffer_1100_sspd[161];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_x_0_0_0_0_x_xx, g_x_x_0_0_0_0_x_xy, g_x_x_0_0_0_0_x_xz, g_x_x_0_0_0_0_x_yy, g_x_x_0_0_0_0_x_yz, g_x_x_0_0_0_0_x_zz, g_x_x_x_xx, g_x_x_x_xy, g_x_x_x_xz, g_x_x_x_yy, g_x_x_x_yz, g_x_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_0_x_xx[i] = 4.0 * g_x_x_x_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_x_xy[i] = 4.0 * g_x_x_x_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_x_xz[i] = 4.0 * g_x_x_x_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_x_yy[i] = 4.0 * g_x_x_x_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_x_yz[i] = 4.0 * g_x_x_x_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_x_zz[i] = 4.0 * g_x_x_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_x_0_0_0_0_y_xx, g_x_x_0_0_0_0_y_xy, g_x_x_0_0_0_0_y_xz, g_x_x_0_0_0_0_y_yy, g_x_x_0_0_0_0_y_yz, g_x_x_0_0_0_0_y_zz, g_x_x_y_xx, g_x_x_y_xy, g_x_x_y_xz, g_x_x_y_yy, g_x_x_y_yz, g_x_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_0_y_xx[i] = 4.0 * g_x_x_y_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_y_xy[i] = 4.0 * g_x_x_y_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_y_xz[i] = 4.0 * g_x_x_y_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_y_yy[i] = 4.0 * g_x_x_y_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_y_yz[i] = 4.0 * g_x_x_y_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_y_zz[i] = 4.0 * g_x_x_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_x_0_0_0_0_z_xx, g_x_x_0_0_0_0_z_xy, g_x_x_0_0_0_0_z_xz, g_x_x_0_0_0_0_z_yy, g_x_x_0_0_0_0_z_yz, g_x_x_0_0_0_0_z_zz, g_x_x_z_xx, g_x_x_z_xy, g_x_x_z_xz, g_x_x_z_yy, g_x_x_z_yz, g_x_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_0_z_xx[i] = 4.0 * g_x_x_z_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_z_xy[i] = 4.0 * g_x_x_z_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_z_xz[i] = 4.0 * g_x_x_z_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_z_yy[i] = 4.0 * g_x_x_z_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_z_yz[i] = 4.0 * g_x_x_z_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_z_zz[i] = 4.0 * g_x_x_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_x_y_0_0_0_0_x_xx, g_x_y_0_0_0_0_x_xy, g_x_y_0_0_0_0_x_xz, g_x_y_0_0_0_0_x_yy, g_x_y_0_0_0_0_x_yz, g_x_y_0_0_0_0_x_zz, g_x_y_x_xx, g_x_y_x_xy, g_x_y_x_xz, g_x_y_x_yy, g_x_y_x_yz, g_x_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_0_x_xx[i] = 4.0 * g_x_y_x_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_x_xy[i] = 4.0 * g_x_y_x_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_x_xz[i] = 4.0 * g_x_y_x_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_x_yy[i] = 4.0 * g_x_y_x_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_x_yz[i] = 4.0 * g_x_y_x_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_x_zz[i] = 4.0 * g_x_y_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_x_y_0_0_0_0_y_xx, g_x_y_0_0_0_0_y_xy, g_x_y_0_0_0_0_y_xz, g_x_y_0_0_0_0_y_yy, g_x_y_0_0_0_0_y_yz, g_x_y_0_0_0_0_y_zz, g_x_y_y_xx, g_x_y_y_xy, g_x_y_y_xz, g_x_y_y_yy, g_x_y_y_yz, g_x_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_0_y_xx[i] = 4.0 * g_x_y_y_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_y_xy[i] = 4.0 * g_x_y_y_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_y_xz[i] = 4.0 * g_x_y_y_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_y_yy[i] = 4.0 * g_x_y_y_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_y_yz[i] = 4.0 * g_x_y_y_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_y_zz[i] = 4.0 * g_x_y_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_x_y_0_0_0_0_z_xx, g_x_y_0_0_0_0_z_xy, g_x_y_0_0_0_0_z_xz, g_x_y_0_0_0_0_z_yy, g_x_y_0_0_0_0_z_yz, g_x_y_0_0_0_0_z_zz, g_x_y_z_xx, g_x_y_z_xy, g_x_y_z_xz, g_x_y_z_yy, g_x_y_z_yz, g_x_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_0_z_xx[i] = 4.0 * g_x_y_z_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_z_xy[i] = 4.0 * g_x_y_z_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_z_xz[i] = 4.0 * g_x_y_z_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_z_yy[i] = 4.0 * g_x_y_z_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_z_yz[i] = 4.0 * g_x_y_z_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_z_zz[i] = 4.0 * g_x_y_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_x_z_0_0_0_0_x_xx, g_x_z_0_0_0_0_x_xy, g_x_z_0_0_0_0_x_xz, g_x_z_0_0_0_0_x_yy, g_x_z_0_0_0_0_x_yz, g_x_z_0_0_0_0_x_zz, g_x_z_x_xx, g_x_z_x_xy, g_x_z_x_xz, g_x_z_x_yy, g_x_z_x_yz, g_x_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_0_x_xx[i] = 4.0 * g_x_z_x_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_x_xy[i] = 4.0 * g_x_z_x_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_x_xz[i] = 4.0 * g_x_z_x_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_x_yy[i] = 4.0 * g_x_z_x_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_x_yz[i] = 4.0 * g_x_z_x_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_x_zz[i] = 4.0 * g_x_z_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_x_z_0_0_0_0_y_xx, g_x_z_0_0_0_0_y_xy, g_x_z_0_0_0_0_y_xz, g_x_z_0_0_0_0_y_yy, g_x_z_0_0_0_0_y_yz, g_x_z_0_0_0_0_y_zz, g_x_z_y_xx, g_x_z_y_xy, g_x_z_y_xz, g_x_z_y_yy, g_x_z_y_yz, g_x_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_0_y_xx[i] = 4.0 * g_x_z_y_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_y_xy[i] = 4.0 * g_x_z_y_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_y_xz[i] = 4.0 * g_x_z_y_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_y_yy[i] = 4.0 * g_x_z_y_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_y_yz[i] = 4.0 * g_x_z_y_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_y_zz[i] = 4.0 * g_x_z_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_x_z_0_0_0_0_z_xx, g_x_z_0_0_0_0_z_xy, g_x_z_0_0_0_0_z_xz, g_x_z_0_0_0_0_z_yy, g_x_z_0_0_0_0_z_yz, g_x_z_0_0_0_0_z_zz, g_x_z_z_xx, g_x_z_z_xy, g_x_z_z_xz, g_x_z_z_yy, g_x_z_z_yz, g_x_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_0_z_xx[i] = 4.0 * g_x_z_z_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_z_xy[i] = 4.0 * g_x_z_z_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_z_xz[i] = 4.0 * g_x_z_z_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_z_yy[i] = 4.0 * g_x_z_z_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_z_yz[i] = 4.0 * g_x_z_z_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_z_zz[i] = 4.0 * g_x_z_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_y_x_0_0_0_0_x_xx, g_y_x_0_0_0_0_x_xy, g_y_x_0_0_0_0_x_xz, g_y_x_0_0_0_0_x_yy, g_y_x_0_0_0_0_x_yz, g_y_x_0_0_0_0_x_zz, g_y_x_x_xx, g_y_x_x_xy, g_y_x_x_xz, g_y_x_x_yy, g_y_x_x_yz, g_y_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_0_x_xx[i] = 4.0 * g_y_x_x_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_x_xy[i] = 4.0 * g_y_x_x_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_x_xz[i] = 4.0 * g_y_x_x_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_x_yy[i] = 4.0 * g_y_x_x_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_x_yz[i] = 4.0 * g_y_x_x_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_x_zz[i] = 4.0 * g_y_x_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_y_x_0_0_0_0_y_xx, g_y_x_0_0_0_0_y_xy, g_y_x_0_0_0_0_y_xz, g_y_x_0_0_0_0_y_yy, g_y_x_0_0_0_0_y_yz, g_y_x_0_0_0_0_y_zz, g_y_x_y_xx, g_y_x_y_xy, g_y_x_y_xz, g_y_x_y_yy, g_y_x_y_yz, g_y_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_0_y_xx[i] = 4.0 * g_y_x_y_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_y_xy[i] = 4.0 * g_y_x_y_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_y_xz[i] = 4.0 * g_y_x_y_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_y_yy[i] = 4.0 * g_y_x_y_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_y_yz[i] = 4.0 * g_y_x_y_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_y_zz[i] = 4.0 * g_y_x_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_y_x_0_0_0_0_z_xx, g_y_x_0_0_0_0_z_xy, g_y_x_0_0_0_0_z_xz, g_y_x_0_0_0_0_z_yy, g_y_x_0_0_0_0_z_yz, g_y_x_0_0_0_0_z_zz, g_y_x_z_xx, g_y_x_z_xy, g_y_x_z_xz, g_y_x_z_yy, g_y_x_z_yz, g_y_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_0_z_xx[i] = 4.0 * g_y_x_z_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_z_xy[i] = 4.0 * g_y_x_z_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_z_xz[i] = 4.0 * g_y_x_z_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_z_yy[i] = 4.0 * g_y_x_z_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_z_yz[i] = 4.0 * g_y_x_z_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_z_zz[i] = 4.0 * g_y_x_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_y_y_0_0_0_0_x_xx, g_y_y_0_0_0_0_x_xy, g_y_y_0_0_0_0_x_xz, g_y_y_0_0_0_0_x_yy, g_y_y_0_0_0_0_x_yz, g_y_y_0_0_0_0_x_zz, g_y_y_x_xx, g_y_y_x_xy, g_y_y_x_xz, g_y_y_x_yy, g_y_y_x_yz, g_y_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_0_x_xx[i] = 4.0 * g_y_y_x_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_x_xy[i] = 4.0 * g_y_y_x_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_x_xz[i] = 4.0 * g_y_y_x_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_x_yy[i] = 4.0 * g_y_y_x_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_x_yz[i] = 4.0 * g_y_y_x_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_x_zz[i] = 4.0 * g_y_y_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_y_y_0_0_0_0_y_xx, g_y_y_0_0_0_0_y_xy, g_y_y_0_0_0_0_y_xz, g_y_y_0_0_0_0_y_yy, g_y_y_0_0_0_0_y_yz, g_y_y_0_0_0_0_y_zz, g_y_y_y_xx, g_y_y_y_xy, g_y_y_y_xz, g_y_y_y_yy, g_y_y_y_yz, g_y_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_0_y_xx[i] = 4.0 * g_y_y_y_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_y_xy[i] = 4.0 * g_y_y_y_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_y_xz[i] = 4.0 * g_y_y_y_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_y_yy[i] = 4.0 * g_y_y_y_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_y_yz[i] = 4.0 * g_y_y_y_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_y_zz[i] = 4.0 * g_y_y_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_y_y_0_0_0_0_z_xx, g_y_y_0_0_0_0_z_xy, g_y_y_0_0_0_0_z_xz, g_y_y_0_0_0_0_z_yy, g_y_y_0_0_0_0_z_yz, g_y_y_0_0_0_0_z_zz, g_y_y_z_xx, g_y_y_z_xy, g_y_y_z_xz, g_y_y_z_yy, g_y_y_z_yz, g_y_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_0_z_xx[i] = 4.0 * g_y_y_z_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_z_xy[i] = 4.0 * g_y_y_z_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_z_xz[i] = 4.0 * g_y_y_z_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_z_yy[i] = 4.0 * g_y_y_z_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_z_yz[i] = 4.0 * g_y_y_z_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_z_zz[i] = 4.0 * g_y_y_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_y_z_0_0_0_0_x_xx, g_y_z_0_0_0_0_x_xy, g_y_z_0_0_0_0_x_xz, g_y_z_0_0_0_0_x_yy, g_y_z_0_0_0_0_x_yz, g_y_z_0_0_0_0_x_zz, g_y_z_x_xx, g_y_z_x_xy, g_y_z_x_xz, g_y_z_x_yy, g_y_z_x_yz, g_y_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_0_x_xx[i] = 4.0 * g_y_z_x_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_x_xy[i] = 4.0 * g_y_z_x_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_x_xz[i] = 4.0 * g_y_z_x_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_x_yy[i] = 4.0 * g_y_z_x_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_x_yz[i] = 4.0 * g_y_z_x_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_x_zz[i] = 4.0 * g_y_z_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_y_z_0_0_0_0_y_xx, g_y_z_0_0_0_0_y_xy, g_y_z_0_0_0_0_y_xz, g_y_z_0_0_0_0_y_yy, g_y_z_0_0_0_0_y_yz, g_y_z_0_0_0_0_y_zz, g_y_z_y_xx, g_y_z_y_xy, g_y_z_y_xz, g_y_z_y_yy, g_y_z_y_yz, g_y_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_0_y_xx[i] = 4.0 * g_y_z_y_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_y_xy[i] = 4.0 * g_y_z_y_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_y_xz[i] = 4.0 * g_y_z_y_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_y_yy[i] = 4.0 * g_y_z_y_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_y_yz[i] = 4.0 * g_y_z_y_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_y_zz[i] = 4.0 * g_y_z_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_y_z_0_0_0_0_z_xx, g_y_z_0_0_0_0_z_xy, g_y_z_0_0_0_0_z_xz, g_y_z_0_0_0_0_z_yy, g_y_z_0_0_0_0_z_yz, g_y_z_0_0_0_0_z_zz, g_y_z_z_xx, g_y_z_z_xy, g_y_z_z_xz, g_y_z_z_yy, g_y_z_z_yz, g_y_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_0_z_xx[i] = 4.0 * g_y_z_z_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_z_xy[i] = 4.0 * g_y_z_z_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_z_xz[i] = 4.0 * g_y_z_z_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_z_yy[i] = 4.0 * g_y_z_z_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_z_yz[i] = 4.0 * g_y_z_z_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_z_zz[i] = 4.0 * g_y_z_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_z_x_0_0_0_0_x_xx, g_z_x_0_0_0_0_x_xy, g_z_x_0_0_0_0_x_xz, g_z_x_0_0_0_0_x_yy, g_z_x_0_0_0_0_x_yz, g_z_x_0_0_0_0_x_zz, g_z_x_x_xx, g_z_x_x_xy, g_z_x_x_xz, g_z_x_x_yy, g_z_x_x_yz, g_z_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_0_x_xx[i] = 4.0 * g_z_x_x_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_x_xy[i] = 4.0 * g_z_x_x_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_x_xz[i] = 4.0 * g_z_x_x_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_x_yy[i] = 4.0 * g_z_x_x_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_x_yz[i] = 4.0 * g_z_x_x_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_x_zz[i] = 4.0 * g_z_x_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_z_x_0_0_0_0_y_xx, g_z_x_0_0_0_0_y_xy, g_z_x_0_0_0_0_y_xz, g_z_x_0_0_0_0_y_yy, g_z_x_0_0_0_0_y_yz, g_z_x_0_0_0_0_y_zz, g_z_x_y_xx, g_z_x_y_xy, g_z_x_y_xz, g_z_x_y_yy, g_z_x_y_yz, g_z_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_0_y_xx[i] = 4.0 * g_z_x_y_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_y_xy[i] = 4.0 * g_z_x_y_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_y_xz[i] = 4.0 * g_z_x_y_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_y_yy[i] = 4.0 * g_z_x_y_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_y_yz[i] = 4.0 * g_z_x_y_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_y_zz[i] = 4.0 * g_z_x_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_z_x_0_0_0_0_z_xx, g_z_x_0_0_0_0_z_xy, g_z_x_0_0_0_0_z_xz, g_z_x_0_0_0_0_z_yy, g_z_x_0_0_0_0_z_yz, g_z_x_0_0_0_0_z_zz, g_z_x_z_xx, g_z_x_z_xy, g_z_x_z_xz, g_z_x_z_yy, g_z_x_z_yz, g_z_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_0_z_xx[i] = 4.0 * g_z_x_z_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_z_xy[i] = 4.0 * g_z_x_z_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_z_xz[i] = 4.0 * g_z_x_z_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_z_yy[i] = 4.0 * g_z_x_z_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_z_yz[i] = 4.0 * g_z_x_z_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_z_zz[i] = 4.0 * g_z_x_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_z_y_0_0_0_0_x_xx, g_z_y_0_0_0_0_x_xy, g_z_y_0_0_0_0_x_xz, g_z_y_0_0_0_0_x_yy, g_z_y_0_0_0_0_x_yz, g_z_y_0_0_0_0_x_zz, g_z_y_x_xx, g_z_y_x_xy, g_z_y_x_xz, g_z_y_x_yy, g_z_y_x_yz, g_z_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_0_x_xx[i] = 4.0 * g_z_y_x_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_x_xy[i] = 4.0 * g_z_y_x_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_x_xz[i] = 4.0 * g_z_y_x_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_x_yy[i] = 4.0 * g_z_y_x_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_x_yz[i] = 4.0 * g_z_y_x_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_x_zz[i] = 4.0 * g_z_y_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_z_y_0_0_0_0_y_xx, g_z_y_0_0_0_0_y_xy, g_z_y_0_0_0_0_y_xz, g_z_y_0_0_0_0_y_yy, g_z_y_0_0_0_0_y_yz, g_z_y_0_0_0_0_y_zz, g_z_y_y_xx, g_z_y_y_xy, g_z_y_y_xz, g_z_y_y_yy, g_z_y_y_yz, g_z_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_0_y_xx[i] = 4.0 * g_z_y_y_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_y_xy[i] = 4.0 * g_z_y_y_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_y_xz[i] = 4.0 * g_z_y_y_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_y_yy[i] = 4.0 * g_z_y_y_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_y_yz[i] = 4.0 * g_z_y_y_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_y_zz[i] = 4.0 * g_z_y_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_z_y_0_0_0_0_z_xx, g_z_y_0_0_0_0_z_xy, g_z_y_0_0_0_0_z_xz, g_z_y_0_0_0_0_z_yy, g_z_y_0_0_0_0_z_yz, g_z_y_0_0_0_0_z_zz, g_z_y_z_xx, g_z_y_z_xy, g_z_y_z_xz, g_z_y_z_yy, g_z_y_z_yz, g_z_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_0_z_xx[i] = 4.0 * g_z_y_z_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_z_xy[i] = 4.0 * g_z_y_z_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_z_xz[i] = 4.0 * g_z_y_z_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_z_yy[i] = 4.0 * g_z_y_z_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_z_yz[i] = 4.0 * g_z_y_z_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_z_zz[i] = 4.0 * g_z_y_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_z_z_0_0_0_0_x_xx, g_z_z_0_0_0_0_x_xy, g_z_z_0_0_0_0_x_xz, g_z_z_0_0_0_0_x_yy, g_z_z_0_0_0_0_x_yz, g_z_z_0_0_0_0_x_zz, g_z_z_x_xx, g_z_z_x_xy, g_z_z_x_xz, g_z_z_x_yy, g_z_z_x_yz, g_z_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_0_x_xx[i] = 4.0 * g_z_z_x_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_x_xy[i] = 4.0 * g_z_z_x_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_x_xz[i] = 4.0 * g_z_z_x_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_x_yy[i] = 4.0 * g_z_z_x_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_x_yz[i] = 4.0 * g_z_z_x_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_x_zz[i] = 4.0 * g_z_z_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_z_z_0_0_0_0_y_xx, g_z_z_0_0_0_0_y_xy, g_z_z_0_0_0_0_y_xz, g_z_z_0_0_0_0_y_yy, g_z_z_0_0_0_0_y_yz, g_z_z_0_0_0_0_y_zz, g_z_z_y_xx, g_z_z_y_xy, g_z_z_y_xz, g_z_z_y_yy, g_z_z_y_yz, g_z_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_0_y_xx[i] = 4.0 * g_z_z_y_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_y_xy[i] = 4.0 * g_z_z_y_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_y_xz[i] = 4.0 * g_z_z_y_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_y_yy[i] = 4.0 * g_z_z_y_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_y_yz[i] = 4.0 * g_z_z_y_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_y_zz[i] = 4.0 * g_z_z_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_z_z_0_0_0_0_z_xx, g_z_z_0_0_0_0_z_xy, g_z_z_0_0_0_0_z_xz, g_z_z_0_0_0_0_z_yy, g_z_z_0_0_0_0_z_yz, g_z_z_0_0_0_0_z_zz, g_z_z_z_xx, g_z_z_z_xy, g_z_z_z_xz, g_z_z_z_yy, g_z_z_z_yz, g_z_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_0_z_xx[i] = 4.0 * g_z_z_z_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_z_xy[i] = 4.0 * g_z_z_z_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_z_xz[i] = 4.0 * g_z_z_z_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_z_yy[i] = 4.0 * g_z_z_z_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_z_yz[i] = 4.0 * g_z_z_z_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_z_zz[i] = 4.0 * g_z_z_z_zz[i] * a_exp * b_exp;
    }
}

} // t4c_geom namespace


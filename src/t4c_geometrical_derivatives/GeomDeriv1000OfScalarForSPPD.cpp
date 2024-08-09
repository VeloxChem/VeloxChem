#include "GeomDeriv1000OfScalarForSPPD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1000_sppd_0(CSimdArray<double>& buffer_1000_sppd,
                     const CSimdArray<double>& buffer_pppd,
                     const double a_exp) -> void
{
    const auto ndims = buffer_1000_sppd.number_of_columns();

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

    /// Set up components of integrals buffer : buffer_1000_sppd

    auto g_x_0_0_0_0_x_x_xx = buffer_1000_sppd[0];

    auto g_x_0_0_0_0_x_x_xy = buffer_1000_sppd[1];

    auto g_x_0_0_0_0_x_x_xz = buffer_1000_sppd[2];

    auto g_x_0_0_0_0_x_x_yy = buffer_1000_sppd[3];

    auto g_x_0_0_0_0_x_x_yz = buffer_1000_sppd[4];

    auto g_x_0_0_0_0_x_x_zz = buffer_1000_sppd[5];

    auto g_x_0_0_0_0_x_y_xx = buffer_1000_sppd[6];

    auto g_x_0_0_0_0_x_y_xy = buffer_1000_sppd[7];

    auto g_x_0_0_0_0_x_y_xz = buffer_1000_sppd[8];

    auto g_x_0_0_0_0_x_y_yy = buffer_1000_sppd[9];

    auto g_x_0_0_0_0_x_y_yz = buffer_1000_sppd[10];

    auto g_x_0_0_0_0_x_y_zz = buffer_1000_sppd[11];

    auto g_x_0_0_0_0_x_z_xx = buffer_1000_sppd[12];

    auto g_x_0_0_0_0_x_z_xy = buffer_1000_sppd[13];

    auto g_x_0_0_0_0_x_z_xz = buffer_1000_sppd[14];

    auto g_x_0_0_0_0_x_z_yy = buffer_1000_sppd[15];

    auto g_x_0_0_0_0_x_z_yz = buffer_1000_sppd[16];

    auto g_x_0_0_0_0_x_z_zz = buffer_1000_sppd[17];

    auto g_x_0_0_0_0_y_x_xx = buffer_1000_sppd[18];

    auto g_x_0_0_0_0_y_x_xy = buffer_1000_sppd[19];

    auto g_x_0_0_0_0_y_x_xz = buffer_1000_sppd[20];

    auto g_x_0_0_0_0_y_x_yy = buffer_1000_sppd[21];

    auto g_x_0_0_0_0_y_x_yz = buffer_1000_sppd[22];

    auto g_x_0_0_0_0_y_x_zz = buffer_1000_sppd[23];

    auto g_x_0_0_0_0_y_y_xx = buffer_1000_sppd[24];

    auto g_x_0_0_0_0_y_y_xy = buffer_1000_sppd[25];

    auto g_x_0_0_0_0_y_y_xz = buffer_1000_sppd[26];

    auto g_x_0_0_0_0_y_y_yy = buffer_1000_sppd[27];

    auto g_x_0_0_0_0_y_y_yz = buffer_1000_sppd[28];

    auto g_x_0_0_0_0_y_y_zz = buffer_1000_sppd[29];

    auto g_x_0_0_0_0_y_z_xx = buffer_1000_sppd[30];

    auto g_x_0_0_0_0_y_z_xy = buffer_1000_sppd[31];

    auto g_x_0_0_0_0_y_z_xz = buffer_1000_sppd[32];

    auto g_x_0_0_0_0_y_z_yy = buffer_1000_sppd[33];

    auto g_x_0_0_0_0_y_z_yz = buffer_1000_sppd[34];

    auto g_x_0_0_0_0_y_z_zz = buffer_1000_sppd[35];

    auto g_x_0_0_0_0_z_x_xx = buffer_1000_sppd[36];

    auto g_x_0_0_0_0_z_x_xy = buffer_1000_sppd[37];

    auto g_x_0_0_0_0_z_x_xz = buffer_1000_sppd[38];

    auto g_x_0_0_0_0_z_x_yy = buffer_1000_sppd[39];

    auto g_x_0_0_0_0_z_x_yz = buffer_1000_sppd[40];

    auto g_x_0_0_0_0_z_x_zz = buffer_1000_sppd[41];

    auto g_x_0_0_0_0_z_y_xx = buffer_1000_sppd[42];

    auto g_x_0_0_0_0_z_y_xy = buffer_1000_sppd[43];

    auto g_x_0_0_0_0_z_y_xz = buffer_1000_sppd[44];

    auto g_x_0_0_0_0_z_y_yy = buffer_1000_sppd[45];

    auto g_x_0_0_0_0_z_y_yz = buffer_1000_sppd[46];

    auto g_x_0_0_0_0_z_y_zz = buffer_1000_sppd[47];

    auto g_x_0_0_0_0_z_z_xx = buffer_1000_sppd[48];

    auto g_x_0_0_0_0_z_z_xy = buffer_1000_sppd[49];

    auto g_x_0_0_0_0_z_z_xz = buffer_1000_sppd[50];

    auto g_x_0_0_0_0_z_z_yy = buffer_1000_sppd[51];

    auto g_x_0_0_0_0_z_z_yz = buffer_1000_sppd[52];

    auto g_x_0_0_0_0_z_z_zz = buffer_1000_sppd[53];

    auto g_y_0_0_0_0_x_x_xx = buffer_1000_sppd[54];

    auto g_y_0_0_0_0_x_x_xy = buffer_1000_sppd[55];

    auto g_y_0_0_0_0_x_x_xz = buffer_1000_sppd[56];

    auto g_y_0_0_0_0_x_x_yy = buffer_1000_sppd[57];

    auto g_y_0_0_0_0_x_x_yz = buffer_1000_sppd[58];

    auto g_y_0_0_0_0_x_x_zz = buffer_1000_sppd[59];

    auto g_y_0_0_0_0_x_y_xx = buffer_1000_sppd[60];

    auto g_y_0_0_0_0_x_y_xy = buffer_1000_sppd[61];

    auto g_y_0_0_0_0_x_y_xz = buffer_1000_sppd[62];

    auto g_y_0_0_0_0_x_y_yy = buffer_1000_sppd[63];

    auto g_y_0_0_0_0_x_y_yz = buffer_1000_sppd[64];

    auto g_y_0_0_0_0_x_y_zz = buffer_1000_sppd[65];

    auto g_y_0_0_0_0_x_z_xx = buffer_1000_sppd[66];

    auto g_y_0_0_0_0_x_z_xy = buffer_1000_sppd[67];

    auto g_y_0_0_0_0_x_z_xz = buffer_1000_sppd[68];

    auto g_y_0_0_0_0_x_z_yy = buffer_1000_sppd[69];

    auto g_y_0_0_0_0_x_z_yz = buffer_1000_sppd[70];

    auto g_y_0_0_0_0_x_z_zz = buffer_1000_sppd[71];

    auto g_y_0_0_0_0_y_x_xx = buffer_1000_sppd[72];

    auto g_y_0_0_0_0_y_x_xy = buffer_1000_sppd[73];

    auto g_y_0_0_0_0_y_x_xz = buffer_1000_sppd[74];

    auto g_y_0_0_0_0_y_x_yy = buffer_1000_sppd[75];

    auto g_y_0_0_0_0_y_x_yz = buffer_1000_sppd[76];

    auto g_y_0_0_0_0_y_x_zz = buffer_1000_sppd[77];

    auto g_y_0_0_0_0_y_y_xx = buffer_1000_sppd[78];

    auto g_y_0_0_0_0_y_y_xy = buffer_1000_sppd[79];

    auto g_y_0_0_0_0_y_y_xz = buffer_1000_sppd[80];

    auto g_y_0_0_0_0_y_y_yy = buffer_1000_sppd[81];

    auto g_y_0_0_0_0_y_y_yz = buffer_1000_sppd[82];

    auto g_y_0_0_0_0_y_y_zz = buffer_1000_sppd[83];

    auto g_y_0_0_0_0_y_z_xx = buffer_1000_sppd[84];

    auto g_y_0_0_0_0_y_z_xy = buffer_1000_sppd[85];

    auto g_y_0_0_0_0_y_z_xz = buffer_1000_sppd[86];

    auto g_y_0_0_0_0_y_z_yy = buffer_1000_sppd[87];

    auto g_y_0_0_0_0_y_z_yz = buffer_1000_sppd[88];

    auto g_y_0_0_0_0_y_z_zz = buffer_1000_sppd[89];

    auto g_y_0_0_0_0_z_x_xx = buffer_1000_sppd[90];

    auto g_y_0_0_0_0_z_x_xy = buffer_1000_sppd[91];

    auto g_y_0_0_0_0_z_x_xz = buffer_1000_sppd[92];

    auto g_y_0_0_0_0_z_x_yy = buffer_1000_sppd[93];

    auto g_y_0_0_0_0_z_x_yz = buffer_1000_sppd[94];

    auto g_y_0_0_0_0_z_x_zz = buffer_1000_sppd[95];

    auto g_y_0_0_0_0_z_y_xx = buffer_1000_sppd[96];

    auto g_y_0_0_0_0_z_y_xy = buffer_1000_sppd[97];

    auto g_y_0_0_0_0_z_y_xz = buffer_1000_sppd[98];

    auto g_y_0_0_0_0_z_y_yy = buffer_1000_sppd[99];

    auto g_y_0_0_0_0_z_y_yz = buffer_1000_sppd[100];

    auto g_y_0_0_0_0_z_y_zz = buffer_1000_sppd[101];

    auto g_y_0_0_0_0_z_z_xx = buffer_1000_sppd[102];

    auto g_y_0_0_0_0_z_z_xy = buffer_1000_sppd[103];

    auto g_y_0_0_0_0_z_z_xz = buffer_1000_sppd[104];

    auto g_y_0_0_0_0_z_z_yy = buffer_1000_sppd[105];

    auto g_y_0_0_0_0_z_z_yz = buffer_1000_sppd[106];

    auto g_y_0_0_0_0_z_z_zz = buffer_1000_sppd[107];

    auto g_z_0_0_0_0_x_x_xx = buffer_1000_sppd[108];

    auto g_z_0_0_0_0_x_x_xy = buffer_1000_sppd[109];

    auto g_z_0_0_0_0_x_x_xz = buffer_1000_sppd[110];

    auto g_z_0_0_0_0_x_x_yy = buffer_1000_sppd[111];

    auto g_z_0_0_0_0_x_x_yz = buffer_1000_sppd[112];

    auto g_z_0_0_0_0_x_x_zz = buffer_1000_sppd[113];

    auto g_z_0_0_0_0_x_y_xx = buffer_1000_sppd[114];

    auto g_z_0_0_0_0_x_y_xy = buffer_1000_sppd[115];

    auto g_z_0_0_0_0_x_y_xz = buffer_1000_sppd[116];

    auto g_z_0_0_0_0_x_y_yy = buffer_1000_sppd[117];

    auto g_z_0_0_0_0_x_y_yz = buffer_1000_sppd[118];

    auto g_z_0_0_0_0_x_y_zz = buffer_1000_sppd[119];

    auto g_z_0_0_0_0_x_z_xx = buffer_1000_sppd[120];

    auto g_z_0_0_0_0_x_z_xy = buffer_1000_sppd[121];

    auto g_z_0_0_0_0_x_z_xz = buffer_1000_sppd[122];

    auto g_z_0_0_0_0_x_z_yy = buffer_1000_sppd[123];

    auto g_z_0_0_0_0_x_z_yz = buffer_1000_sppd[124];

    auto g_z_0_0_0_0_x_z_zz = buffer_1000_sppd[125];

    auto g_z_0_0_0_0_y_x_xx = buffer_1000_sppd[126];

    auto g_z_0_0_0_0_y_x_xy = buffer_1000_sppd[127];

    auto g_z_0_0_0_0_y_x_xz = buffer_1000_sppd[128];

    auto g_z_0_0_0_0_y_x_yy = buffer_1000_sppd[129];

    auto g_z_0_0_0_0_y_x_yz = buffer_1000_sppd[130];

    auto g_z_0_0_0_0_y_x_zz = buffer_1000_sppd[131];

    auto g_z_0_0_0_0_y_y_xx = buffer_1000_sppd[132];

    auto g_z_0_0_0_0_y_y_xy = buffer_1000_sppd[133];

    auto g_z_0_0_0_0_y_y_xz = buffer_1000_sppd[134];

    auto g_z_0_0_0_0_y_y_yy = buffer_1000_sppd[135];

    auto g_z_0_0_0_0_y_y_yz = buffer_1000_sppd[136];

    auto g_z_0_0_0_0_y_y_zz = buffer_1000_sppd[137];

    auto g_z_0_0_0_0_y_z_xx = buffer_1000_sppd[138];

    auto g_z_0_0_0_0_y_z_xy = buffer_1000_sppd[139];

    auto g_z_0_0_0_0_y_z_xz = buffer_1000_sppd[140];

    auto g_z_0_0_0_0_y_z_yy = buffer_1000_sppd[141];

    auto g_z_0_0_0_0_y_z_yz = buffer_1000_sppd[142];

    auto g_z_0_0_0_0_y_z_zz = buffer_1000_sppd[143];

    auto g_z_0_0_0_0_z_x_xx = buffer_1000_sppd[144];

    auto g_z_0_0_0_0_z_x_xy = buffer_1000_sppd[145];

    auto g_z_0_0_0_0_z_x_xz = buffer_1000_sppd[146];

    auto g_z_0_0_0_0_z_x_yy = buffer_1000_sppd[147];

    auto g_z_0_0_0_0_z_x_yz = buffer_1000_sppd[148];

    auto g_z_0_0_0_0_z_x_zz = buffer_1000_sppd[149];

    auto g_z_0_0_0_0_z_y_xx = buffer_1000_sppd[150];

    auto g_z_0_0_0_0_z_y_xy = buffer_1000_sppd[151];

    auto g_z_0_0_0_0_z_y_xz = buffer_1000_sppd[152];

    auto g_z_0_0_0_0_z_y_yy = buffer_1000_sppd[153];

    auto g_z_0_0_0_0_z_y_yz = buffer_1000_sppd[154];

    auto g_z_0_0_0_0_z_y_zz = buffer_1000_sppd[155];

    auto g_z_0_0_0_0_z_z_xx = buffer_1000_sppd[156];

    auto g_z_0_0_0_0_z_z_xy = buffer_1000_sppd[157];

    auto g_z_0_0_0_0_z_z_xz = buffer_1000_sppd[158];

    auto g_z_0_0_0_0_z_z_yy = buffer_1000_sppd[159];

    auto g_z_0_0_0_0_z_z_yz = buffer_1000_sppd[160];

    auto g_z_0_0_0_0_z_z_zz = buffer_1000_sppd[161];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_0_0_0_0_x_x_xx, g_x_0_0_0_0_x_x_xy, g_x_0_0_0_0_x_x_xz, g_x_0_0_0_0_x_x_yy, g_x_0_0_0_0_x_x_yz, g_x_0_0_0_0_x_x_zz, g_x_x_x_xx, g_x_x_x_xy, g_x_x_x_xz, g_x_x_x_yy, g_x_x_x_yz, g_x_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_x_x_xx[i] = 2.0 * g_x_x_x_xx[i] * a_exp;

        g_x_0_0_0_0_x_x_xy[i] = 2.0 * g_x_x_x_xy[i] * a_exp;

        g_x_0_0_0_0_x_x_xz[i] = 2.0 * g_x_x_x_xz[i] * a_exp;

        g_x_0_0_0_0_x_x_yy[i] = 2.0 * g_x_x_x_yy[i] * a_exp;

        g_x_0_0_0_0_x_x_yz[i] = 2.0 * g_x_x_x_yz[i] * a_exp;

        g_x_0_0_0_0_x_x_zz[i] = 2.0 * g_x_x_x_zz[i] * a_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_0_0_0_0_x_y_xx, g_x_0_0_0_0_x_y_xy, g_x_0_0_0_0_x_y_xz, g_x_0_0_0_0_x_y_yy, g_x_0_0_0_0_x_y_yz, g_x_0_0_0_0_x_y_zz, g_x_x_y_xx, g_x_x_y_xy, g_x_x_y_xz, g_x_x_y_yy, g_x_x_y_yz, g_x_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_x_y_xx[i] = 2.0 * g_x_x_y_xx[i] * a_exp;

        g_x_0_0_0_0_x_y_xy[i] = 2.0 * g_x_x_y_xy[i] * a_exp;

        g_x_0_0_0_0_x_y_xz[i] = 2.0 * g_x_x_y_xz[i] * a_exp;

        g_x_0_0_0_0_x_y_yy[i] = 2.0 * g_x_x_y_yy[i] * a_exp;

        g_x_0_0_0_0_x_y_yz[i] = 2.0 * g_x_x_y_yz[i] * a_exp;

        g_x_0_0_0_0_x_y_zz[i] = 2.0 * g_x_x_y_zz[i] * a_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_0_0_0_0_x_z_xx, g_x_0_0_0_0_x_z_xy, g_x_0_0_0_0_x_z_xz, g_x_0_0_0_0_x_z_yy, g_x_0_0_0_0_x_z_yz, g_x_0_0_0_0_x_z_zz, g_x_x_z_xx, g_x_x_z_xy, g_x_x_z_xz, g_x_x_z_yy, g_x_x_z_yz, g_x_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_x_z_xx[i] = 2.0 * g_x_x_z_xx[i] * a_exp;

        g_x_0_0_0_0_x_z_xy[i] = 2.0 * g_x_x_z_xy[i] * a_exp;

        g_x_0_0_0_0_x_z_xz[i] = 2.0 * g_x_x_z_xz[i] * a_exp;

        g_x_0_0_0_0_x_z_yy[i] = 2.0 * g_x_x_z_yy[i] * a_exp;

        g_x_0_0_0_0_x_z_yz[i] = 2.0 * g_x_x_z_yz[i] * a_exp;

        g_x_0_0_0_0_x_z_zz[i] = 2.0 * g_x_x_z_zz[i] * a_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_x_0_0_0_0_y_x_xx, g_x_0_0_0_0_y_x_xy, g_x_0_0_0_0_y_x_xz, g_x_0_0_0_0_y_x_yy, g_x_0_0_0_0_y_x_yz, g_x_0_0_0_0_y_x_zz, g_x_y_x_xx, g_x_y_x_xy, g_x_y_x_xz, g_x_y_x_yy, g_x_y_x_yz, g_x_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_y_x_xx[i] = 2.0 * g_x_y_x_xx[i] * a_exp;

        g_x_0_0_0_0_y_x_xy[i] = 2.0 * g_x_y_x_xy[i] * a_exp;

        g_x_0_0_0_0_y_x_xz[i] = 2.0 * g_x_y_x_xz[i] * a_exp;

        g_x_0_0_0_0_y_x_yy[i] = 2.0 * g_x_y_x_yy[i] * a_exp;

        g_x_0_0_0_0_y_x_yz[i] = 2.0 * g_x_y_x_yz[i] * a_exp;

        g_x_0_0_0_0_y_x_zz[i] = 2.0 * g_x_y_x_zz[i] * a_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_x_0_0_0_0_y_y_xx, g_x_0_0_0_0_y_y_xy, g_x_0_0_0_0_y_y_xz, g_x_0_0_0_0_y_y_yy, g_x_0_0_0_0_y_y_yz, g_x_0_0_0_0_y_y_zz, g_x_y_y_xx, g_x_y_y_xy, g_x_y_y_xz, g_x_y_y_yy, g_x_y_y_yz, g_x_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_y_y_xx[i] = 2.0 * g_x_y_y_xx[i] * a_exp;

        g_x_0_0_0_0_y_y_xy[i] = 2.0 * g_x_y_y_xy[i] * a_exp;

        g_x_0_0_0_0_y_y_xz[i] = 2.0 * g_x_y_y_xz[i] * a_exp;

        g_x_0_0_0_0_y_y_yy[i] = 2.0 * g_x_y_y_yy[i] * a_exp;

        g_x_0_0_0_0_y_y_yz[i] = 2.0 * g_x_y_y_yz[i] * a_exp;

        g_x_0_0_0_0_y_y_zz[i] = 2.0 * g_x_y_y_zz[i] * a_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_x_0_0_0_0_y_z_xx, g_x_0_0_0_0_y_z_xy, g_x_0_0_0_0_y_z_xz, g_x_0_0_0_0_y_z_yy, g_x_0_0_0_0_y_z_yz, g_x_0_0_0_0_y_z_zz, g_x_y_z_xx, g_x_y_z_xy, g_x_y_z_xz, g_x_y_z_yy, g_x_y_z_yz, g_x_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_y_z_xx[i] = 2.0 * g_x_y_z_xx[i] * a_exp;

        g_x_0_0_0_0_y_z_xy[i] = 2.0 * g_x_y_z_xy[i] * a_exp;

        g_x_0_0_0_0_y_z_xz[i] = 2.0 * g_x_y_z_xz[i] * a_exp;

        g_x_0_0_0_0_y_z_yy[i] = 2.0 * g_x_y_z_yy[i] * a_exp;

        g_x_0_0_0_0_y_z_yz[i] = 2.0 * g_x_y_z_yz[i] * a_exp;

        g_x_0_0_0_0_y_z_zz[i] = 2.0 * g_x_y_z_zz[i] * a_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_x_0_0_0_0_z_x_xx, g_x_0_0_0_0_z_x_xy, g_x_0_0_0_0_z_x_xz, g_x_0_0_0_0_z_x_yy, g_x_0_0_0_0_z_x_yz, g_x_0_0_0_0_z_x_zz, g_x_z_x_xx, g_x_z_x_xy, g_x_z_x_xz, g_x_z_x_yy, g_x_z_x_yz, g_x_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_z_x_xx[i] = 2.0 * g_x_z_x_xx[i] * a_exp;

        g_x_0_0_0_0_z_x_xy[i] = 2.0 * g_x_z_x_xy[i] * a_exp;

        g_x_0_0_0_0_z_x_xz[i] = 2.0 * g_x_z_x_xz[i] * a_exp;

        g_x_0_0_0_0_z_x_yy[i] = 2.0 * g_x_z_x_yy[i] * a_exp;

        g_x_0_0_0_0_z_x_yz[i] = 2.0 * g_x_z_x_yz[i] * a_exp;

        g_x_0_0_0_0_z_x_zz[i] = 2.0 * g_x_z_x_zz[i] * a_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_x_0_0_0_0_z_y_xx, g_x_0_0_0_0_z_y_xy, g_x_0_0_0_0_z_y_xz, g_x_0_0_0_0_z_y_yy, g_x_0_0_0_0_z_y_yz, g_x_0_0_0_0_z_y_zz, g_x_z_y_xx, g_x_z_y_xy, g_x_z_y_xz, g_x_z_y_yy, g_x_z_y_yz, g_x_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_z_y_xx[i] = 2.0 * g_x_z_y_xx[i] * a_exp;

        g_x_0_0_0_0_z_y_xy[i] = 2.0 * g_x_z_y_xy[i] * a_exp;

        g_x_0_0_0_0_z_y_xz[i] = 2.0 * g_x_z_y_xz[i] * a_exp;

        g_x_0_0_0_0_z_y_yy[i] = 2.0 * g_x_z_y_yy[i] * a_exp;

        g_x_0_0_0_0_z_y_yz[i] = 2.0 * g_x_z_y_yz[i] * a_exp;

        g_x_0_0_0_0_z_y_zz[i] = 2.0 * g_x_z_y_zz[i] * a_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_x_0_0_0_0_z_z_xx, g_x_0_0_0_0_z_z_xy, g_x_0_0_0_0_z_z_xz, g_x_0_0_0_0_z_z_yy, g_x_0_0_0_0_z_z_yz, g_x_0_0_0_0_z_z_zz, g_x_z_z_xx, g_x_z_z_xy, g_x_z_z_xz, g_x_z_z_yy, g_x_z_z_yz, g_x_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_z_z_xx[i] = 2.0 * g_x_z_z_xx[i] * a_exp;

        g_x_0_0_0_0_z_z_xy[i] = 2.0 * g_x_z_z_xy[i] * a_exp;

        g_x_0_0_0_0_z_z_xz[i] = 2.0 * g_x_z_z_xz[i] * a_exp;

        g_x_0_0_0_0_z_z_yy[i] = 2.0 * g_x_z_z_yy[i] * a_exp;

        g_x_0_0_0_0_z_z_yz[i] = 2.0 * g_x_z_z_yz[i] * a_exp;

        g_x_0_0_0_0_z_z_zz[i] = 2.0 * g_x_z_z_zz[i] * a_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_y_0_0_0_0_x_x_xx, g_y_0_0_0_0_x_x_xy, g_y_0_0_0_0_x_x_xz, g_y_0_0_0_0_x_x_yy, g_y_0_0_0_0_x_x_yz, g_y_0_0_0_0_x_x_zz, g_y_x_x_xx, g_y_x_x_xy, g_y_x_x_xz, g_y_x_x_yy, g_y_x_x_yz, g_y_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_x_x_xx[i] = 2.0 * g_y_x_x_xx[i] * a_exp;

        g_y_0_0_0_0_x_x_xy[i] = 2.0 * g_y_x_x_xy[i] * a_exp;

        g_y_0_0_0_0_x_x_xz[i] = 2.0 * g_y_x_x_xz[i] * a_exp;

        g_y_0_0_0_0_x_x_yy[i] = 2.0 * g_y_x_x_yy[i] * a_exp;

        g_y_0_0_0_0_x_x_yz[i] = 2.0 * g_y_x_x_yz[i] * a_exp;

        g_y_0_0_0_0_x_x_zz[i] = 2.0 * g_y_x_x_zz[i] * a_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_y_0_0_0_0_x_y_xx, g_y_0_0_0_0_x_y_xy, g_y_0_0_0_0_x_y_xz, g_y_0_0_0_0_x_y_yy, g_y_0_0_0_0_x_y_yz, g_y_0_0_0_0_x_y_zz, g_y_x_y_xx, g_y_x_y_xy, g_y_x_y_xz, g_y_x_y_yy, g_y_x_y_yz, g_y_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_x_y_xx[i] = 2.0 * g_y_x_y_xx[i] * a_exp;

        g_y_0_0_0_0_x_y_xy[i] = 2.0 * g_y_x_y_xy[i] * a_exp;

        g_y_0_0_0_0_x_y_xz[i] = 2.0 * g_y_x_y_xz[i] * a_exp;

        g_y_0_0_0_0_x_y_yy[i] = 2.0 * g_y_x_y_yy[i] * a_exp;

        g_y_0_0_0_0_x_y_yz[i] = 2.0 * g_y_x_y_yz[i] * a_exp;

        g_y_0_0_0_0_x_y_zz[i] = 2.0 * g_y_x_y_zz[i] * a_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_y_0_0_0_0_x_z_xx, g_y_0_0_0_0_x_z_xy, g_y_0_0_0_0_x_z_xz, g_y_0_0_0_0_x_z_yy, g_y_0_0_0_0_x_z_yz, g_y_0_0_0_0_x_z_zz, g_y_x_z_xx, g_y_x_z_xy, g_y_x_z_xz, g_y_x_z_yy, g_y_x_z_yz, g_y_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_x_z_xx[i] = 2.0 * g_y_x_z_xx[i] * a_exp;

        g_y_0_0_0_0_x_z_xy[i] = 2.0 * g_y_x_z_xy[i] * a_exp;

        g_y_0_0_0_0_x_z_xz[i] = 2.0 * g_y_x_z_xz[i] * a_exp;

        g_y_0_0_0_0_x_z_yy[i] = 2.0 * g_y_x_z_yy[i] * a_exp;

        g_y_0_0_0_0_x_z_yz[i] = 2.0 * g_y_x_z_yz[i] * a_exp;

        g_y_0_0_0_0_x_z_zz[i] = 2.0 * g_y_x_z_zz[i] * a_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_y_0_0_0_0_y_x_xx, g_y_0_0_0_0_y_x_xy, g_y_0_0_0_0_y_x_xz, g_y_0_0_0_0_y_x_yy, g_y_0_0_0_0_y_x_yz, g_y_0_0_0_0_y_x_zz, g_y_y_x_xx, g_y_y_x_xy, g_y_y_x_xz, g_y_y_x_yy, g_y_y_x_yz, g_y_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_y_x_xx[i] = 2.0 * g_y_y_x_xx[i] * a_exp;

        g_y_0_0_0_0_y_x_xy[i] = 2.0 * g_y_y_x_xy[i] * a_exp;

        g_y_0_0_0_0_y_x_xz[i] = 2.0 * g_y_y_x_xz[i] * a_exp;

        g_y_0_0_0_0_y_x_yy[i] = 2.0 * g_y_y_x_yy[i] * a_exp;

        g_y_0_0_0_0_y_x_yz[i] = 2.0 * g_y_y_x_yz[i] * a_exp;

        g_y_0_0_0_0_y_x_zz[i] = 2.0 * g_y_y_x_zz[i] * a_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_y_0_0_0_0_y_y_xx, g_y_0_0_0_0_y_y_xy, g_y_0_0_0_0_y_y_xz, g_y_0_0_0_0_y_y_yy, g_y_0_0_0_0_y_y_yz, g_y_0_0_0_0_y_y_zz, g_y_y_y_xx, g_y_y_y_xy, g_y_y_y_xz, g_y_y_y_yy, g_y_y_y_yz, g_y_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_y_y_xx[i] = 2.0 * g_y_y_y_xx[i] * a_exp;

        g_y_0_0_0_0_y_y_xy[i] = 2.0 * g_y_y_y_xy[i] * a_exp;

        g_y_0_0_0_0_y_y_xz[i] = 2.0 * g_y_y_y_xz[i] * a_exp;

        g_y_0_0_0_0_y_y_yy[i] = 2.0 * g_y_y_y_yy[i] * a_exp;

        g_y_0_0_0_0_y_y_yz[i] = 2.0 * g_y_y_y_yz[i] * a_exp;

        g_y_0_0_0_0_y_y_zz[i] = 2.0 * g_y_y_y_zz[i] * a_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_y_0_0_0_0_y_z_xx, g_y_0_0_0_0_y_z_xy, g_y_0_0_0_0_y_z_xz, g_y_0_0_0_0_y_z_yy, g_y_0_0_0_0_y_z_yz, g_y_0_0_0_0_y_z_zz, g_y_y_z_xx, g_y_y_z_xy, g_y_y_z_xz, g_y_y_z_yy, g_y_y_z_yz, g_y_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_y_z_xx[i] = 2.0 * g_y_y_z_xx[i] * a_exp;

        g_y_0_0_0_0_y_z_xy[i] = 2.0 * g_y_y_z_xy[i] * a_exp;

        g_y_0_0_0_0_y_z_xz[i] = 2.0 * g_y_y_z_xz[i] * a_exp;

        g_y_0_0_0_0_y_z_yy[i] = 2.0 * g_y_y_z_yy[i] * a_exp;

        g_y_0_0_0_0_y_z_yz[i] = 2.0 * g_y_y_z_yz[i] * a_exp;

        g_y_0_0_0_0_y_z_zz[i] = 2.0 * g_y_y_z_zz[i] * a_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_y_0_0_0_0_z_x_xx, g_y_0_0_0_0_z_x_xy, g_y_0_0_0_0_z_x_xz, g_y_0_0_0_0_z_x_yy, g_y_0_0_0_0_z_x_yz, g_y_0_0_0_0_z_x_zz, g_y_z_x_xx, g_y_z_x_xy, g_y_z_x_xz, g_y_z_x_yy, g_y_z_x_yz, g_y_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_z_x_xx[i] = 2.0 * g_y_z_x_xx[i] * a_exp;

        g_y_0_0_0_0_z_x_xy[i] = 2.0 * g_y_z_x_xy[i] * a_exp;

        g_y_0_0_0_0_z_x_xz[i] = 2.0 * g_y_z_x_xz[i] * a_exp;

        g_y_0_0_0_0_z_x_yy[i] = 2.0 * g_y_z_x_yy[i] * a_exp;

        g_y_0_0_0_0_z_x_yz[i] = 2.0 * g_y_z_x_yz[i] * a_exp;

        g_y_0_0_0_0_z_x_zz[i] = 2.0 * g_y_z_x_zz[i] * a_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_y_0_0_0_0_z_y_xx, g_y_0_0_0_0_z_y_xy, g_y_0_0_0_0_z_y_xz, g_y_0_0_0_0_z_y_yy, g_y_0_0_0_0_z_y_yz, g_y_0_0_0_0_z_y_zz, g_y_z_y_xx, g_y_z_y_xy, g_y_z_y_xz, g_y_z_y_yy, g_y_z_y_yz, g_y_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_z_y_xx[i] = 2.0 * g_y_z_y_xx[i] * a_exp;

        g_y_0_0_0_0_z_y_xy[i] = 2.0 * g_y_z_y_xy[i] * a_exp;

        g_y_0_0_0_0_z_y_xz[i] = 2.0 * g_y_z_y_xz[i] * a_exp;

        g_y_0_0_0_0_z_y_yy[i] = 2.0 * g_y_z_y_yy[i] * a_exp;

        g_y_0_0_0_0_z_y_yz[i] = 2.0 * g_y_z_y_yz[i] * a_exp;

        g_y_0_0_0_0_z_y_zz[i] = 2.0 * g_y_z_y_zz[i] * a_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_y_0_0_0_0_z_z_xx, g_y_0_0_0_0_z_z_xy, g_y_0_0_0_0_z_z_xz, g_y_0_0_0_0_z_z_yy, g_y_0_0_0_0_z_z_yz, g_y_0_0_0_0_z_z_zz, g_y_z_z_xx, g_y_z_z_xy, g_y_z_z_xz, g_y_z_z_yy, g_y_z_z_yz, g_y_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_z_z_xx[i] = 2.0 * g_y_z_z_xx[i] * a_exp;

        g_y_0_0_0_0_z_z_xy[i] = 2.0 * g_y_z_z_xy[i] * a_exp;

        g_y_0_0_0_0_z_z_xz[i] = 2.0 * g_y_z_z_xz[i] * a_exp;

        g_y_0_0_0_0_z_z_yy[i] = 2.0 * g_y_z_z_yy[i] * a_exp;

        g_y_0_0_0_0_z_z_yz[i] = 2.0 * g_y_z_z_yz[i] * a_exp;

        g_y_0_0_0_0_z_z_zz[i] = 2.0 * g_y_z_z_zz[i] * a_exp;
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_z_0_0_0_0_x_x_xx, g_z_0_0_0_0_x_x_xy, g_z_0_0_0_0_x_x_xz, g_z_0_0_0_0_x_x_yy, g_z_0_0_0_0_x_x_yz, g_z_0_0_0_0_x_x_zz, g_z_x_x_xx, g_z_x_x_xy, g_z_x_x_xz, g_z_x_x_yy, g_z_x_x_yz, g_z_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_x_x_xx[i] = 2.0 * g_z_x_x_xx[i] * a_exp;

        g_z_0_0_0_0_x_x_xy[i] = 2.0 * g_z_x_x_xy[i] * a_exp;

        g_z_0_0_0_0_x_x_xz[i] = 2.0 * g_z_x_x_xz[i] * a_exp;

        g_z_0_0_0_0_x_x_yy[i] = 2.0 * g_z_x_x_yy[i] * a_exp;

        g_z_0_0_0_0_x_x_yz[i] = 2.0 * g_z_x_x_yz[i] * a_exp;

        g_z_0_0_0_0_x_x_zz[i] = 2.0 * g_z_x_x_zz[i] * a_exp;
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_z_0_0_0_0_x_y_xx, g_z_0_0_0_0_x_y_xy, g_z_0_0_0_0_x_y_xz, g_z_0_0_0_0_x_y_yy, g_z_0_0_0_0_x_y_yz, g_z_0_0_0_0_x_y_zz, g_z_x_y_xx, g_z_x_y_xy, g_z_x_y_xz, g_z_x_y_yy, g_z_x_y_yz, g_z_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_x_y_xx[i] = 2.0 * g_z_x_y_xx[i] * a_exp;

        g_z_0_0_0_0_x_y_xy[i] = 2.0 * g_z_x_y_xy[i] * a_exp;

        g_z_0_0_0_0_x_y_xz[i] = 2.0 * g_z_x_y_xz[i] * a_exp;

        g_z_0_0_0_0_x_y_yy[i] = 2.0 * g_z_x_y_yy[i] * a_exp;

        g_z_0_0_0_0_x_y_yz[i] = 2.0 * g_z_x_y_yz[i] * a_exp;

        g_z_0_0_0_0_x_y_zz[i] = 2.0 * g_z_x_y_zz[i] * a_exp;
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_z_0_0_0_0_x_z_xx, g_z_0_0_0_0_x_z_xy, g_z_0_0_0_0_x_z_xz, g_z_0_0_0_0_x_z_yy, g_z_0_0_0_0_x_z_yz, g_z_0_0_0_0_x_z_zz, g_z_x_z_xx, g_z_x_z_xy, g_z_x_z_xz, g_z_x_z_yy, g_z_x_z_yz, g_z_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_x_z_xx[i] = 2.0 * g_z_x_z_xx[i] * a_exp;

        g_z_0_0_0_0_x_z_xy[i] = 2.0 * g_z_x_z_xy[i] * a_exp;

        g_z_0_0_0_0_x_z_xz[i] = 2.0 * g_z_x_z_xz[i] * a_exp;

        g_z_0_0_0_0_x_z_yy[i] = 2.0 * g_z_x_z_yy[i] * a_exp;

        g_z_0_0_0_0_x_z_yz[i] = 2.0 * g_z_x_z_yz[i] * a_exp;

        g_z_0_0_0_0_x_z_zz[i] = 2.0 * g_z_x_z_zz[i] * a_exp;
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_z_0_0_0_0_y_x_xx, g_z_0_0_0_0_y_x_xy, g_z_0_0_0_0_y_x_xz, g_z_0_0_0_0_y_x_yy, g_z_0_0_0_0_y_x_yz, g_z_0_0_0_0_y_x_zz, g_z_y_x_xx, g_z_y_x_xy, g_z_y_x_xz, g_z_y_x_yy, g_z_y_x_yz, g_z_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_y_x_xx[i] = 2.0 * g_z_y_x_xx[i] * a_exp;

        g_z_0_0_0_0_y_x_xy[i] = 2.0 * g_z_y_x_xy[i] * a_exp;

        g_z_0_0_0_0_y_x_xz[i] = 2.0 * g_z_y_x_xz[i] * a_exp;

        g_z_0_0_0_0_y_x_yy[i] = 2.0 * g_z_y_x_yy[i] * a_exp;

        g_z_0_0_0_0_y_x_yz[i] = 2.0 * g_z_y_x_yz[i] * a_exp;

        g_z_0_0_0_0_y_x_zz[i] = 2.0 * g_z_y_x_zz[i] * a_exp;
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_z_0_0_0_0_y_y_xx, g_z_0_0_0_0_y_y_xy, g_z_0_0_0_0_y_y_xz, g_z_0_0_0_0_y_y_yy, g_z_0_0_0_0_y_y_yz, g_z_0_0_0_0_y_y_zz, g_z_y_y_xx, g_z_y_y_xy, g_z_y_y_xz, g_z_y_y_yy, g_z_y_y_yz, g_z_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_y_y_xx[i] = 2.0 * g_z_y_y_xx[i] * a_exp;

        g_z_0_0_0_0_y_y_xy[i] = 2.0 * g_z_y_y_xy[i] * a_exp;

        g_z_0_0_0_0_y_y_xz[i] = 2.0 * g_z_y_y_xz[i] * a_exp;

        g_z_0_0_0_0_y_y_yy[i] = 2.0 * g_z_y_y_yy[i] * a_exp;

        g_z_0_0_0_0_y_y_yz[i] = 2.0 * g_z_y_y_yz[i] * a_exp;

        g_z_0_0_0_0_y_y_zz[i] = 2.0 * g_z_y_y_zz[i] * a_exp;
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_z_0_0_0_0_y_z_xx, g_z_0_0_0_0_y_z_xy, g_z_0_0_0_0_y_z_xz, g_z_0_0_0_0_y_z_yy, g_z_0_0_0_0_y_z_yz, g_z_0_0_0_0_y_z_zz, g_z_y_z_xx, g_z_y_z_xy, g_z_y_z_xz, g_z_y_z_yy, g_z_y_z_yz, g_z_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_y_z_xx[i] = 2.0 * g_z_y_z_xx[i] * a_exp;

        g_z_0_0_0_0_y_z_xy[i] = 2.0 * g_z_y_z_xy[i] * a_exp;

        g_z_0_0_0_0_y_z_xz[i] = 2.0 * g_z_y_z_xz[i] * a_exp;

        g_z_0_0_0_0_y_z_yy[i] = 2.0 * g_z_y_z_yy[i] * a_exp;

        g_z_0_0_0_0_y_z_yz[i] = 2.0 * g_z_y_z_yz[i] * a_exp;

        g_z_0_0_0_0_y_z_zz[i] = 2.0 * g_z_y_z_zz[i] * a_exp;
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_z_0_0_0_0_z_x_xx, g_z_0_0_0_0_z_x_xy, g_z_0_0_0_0_z_x_xz, g_z_0_0_0_0_z_x_yy, g_z_0_0_0_0_z_x_yz, g_z_0_0_0_0_z_x_zz, g_z_z_x_xx, g_z_z_x_xy, g_z_z_x_xz, g_z_z_x_yy, g_z_z_x_yz, g_z_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_z_x_xx[i] = 2.0 * g_z_z_x_xx[i] * a_exp;

        g_z_0_0_0_0_z_x_xy[i] = 2.0 * g_z_z_x_xy[i] * a_exp;

        g_z_0_0_0_0_z_x_xz[i] = 2.0 * g_z_z_x_xz[i] * a_exp;

        g_z_0_0_0_0_z_x_yy[i] = 2.0 * g_z_z_x_yy[i] * a_exp;

        g_z_0_0_0_0_z_x_yz[i] = 2.0 * g_z_z_x_yz[i] * a_exp;

        g_z_0_0_0_0_z_x_zz[i] = 2.0 * g_z_z_x_zz[i] * a_exp;
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_z_0_0_0_0_z_y_xx, g_z_0_0_0_0_z_y_xy, g_z_0_0_0_0_z_y_xz, g_z_0_0_0_0_z_y_yy, g_z_0_0_0_0_z_y_yz, g_z_0_0_0_0_z_y_zz, g_z_z_y_xx, g_z_z_y_xy, g_z_z_y_xz, g_z_z_y_yy, g_z_z_y_yz, g_z_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_z_y_xx[i] = 2.0 * g_z_z_y_xx[i] * a_exp;

        g_z_0_0_0_0_z_y_xy[i] = 2.0 * g_z_z_y_xy[i] * a_exp;

        g_z_0_0_0_0_z_y_xz[i] = 2.0 * g_z_z_y_xz[i] * a_exp;

        g_z_0_0_0_0_z_y_yy[i] = 2.0 * g_z_z_y_yy[i] * a_exp;

        g_z_0_0_0_0_z_y_yz[i] = 2.0 * g_z_z_y_yz[i] * a_exp;

        g_z_0_0_0_0_z_y_zz[i] = 2.0 * g_z_z_y_zz[i] * a_exp;
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_z_0_0_0_0_z_z_xx, g_z_0_0_0_0_z_z_xy, g_z_0_0_0_0_z_z_xz, g_z_0_0_0_0_z_z_yy, g_z_0_0_0_0_z_z_yz, g_z_0_0_0_0_z_z_zz, g_z_z_z_xx, g_z_z_z_xy, g_z_z_z_xz, g_z_z_z_yy, g_z_z_z_yz, g_z_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_z_z_xx[i] = 2.0 * g_z_z_z_xx[i] * a_exp;

        g_z_0_0_0_0_z_z_xy[i] = 2.0 * g_z_z_z_xy[i] * a_exp;

        g_z_0_0_0_0_z_z_xz[i] = 2.0 * g_z_z_z_xz[i] * a_exp;

        g_z_0_0_0_0_z_z_yy[i] = 2.0 * g_z_z_z_yy[i] * a_exp;

        g_z_0_0_0_0_z_z_yz[i] = 2.0 * g_z_z_z_yz[i] * a_exp;

        g_z_0_0_0_0_z_z_zz[i] = 2.0 * g_z_z_z_zz[i] * a_exp;
    }
}

} // t4c_geom namespace


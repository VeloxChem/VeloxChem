#include "GeomDeriv1010OfScalarForSPSD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1010_spsd_0(CSimdArray<double>& buffer_1010_spsd,
                     const CSimdArray<double>& buffer_pppd,
                     const double a_exp,
                     const double* c_exps) -> void
{
    const auto ndims = buffer_1010_spsd.number_of_columns();

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

    /// Set up components of integrals buffer : buffer_1010_spsd

    auto g_x_0_x_0_0_x_0_xx = buffer_1010_spsd[0];

    auto g_x_0_x_0_0_x_0_xy = buffer_1010_spsd[1];

    auto g_x_0_x_0_0_x_0_xz = buffer_1010_spsd[2];

    auto g_x_0_x_0_0_x_0_yy = buffer_1010_spsd[3];

    auto g_x_0_x_0_0_x_0_yz = buffer_1010_spsd[4];

    auto g_x_0_x_0_0_x_0_zz = buffer_1010_spsd[5];

    auto g_x_0_x_0_0_y_0_xx = buffer_1010_spsd[6];

    auto g_x_0_x_0_0_y_0_xy = buffer_1010_spsd[7];

    auto g_x_0_x_0_0_y_0_xz = buffer_1010_spsd[8];

    auto g_x_0_x_0_0_y_0_yy = buffer_1010_spsd[9];

    auto g_x_0_x_0_0_y_0_yz = buffer_1010_spsd[10];

    auto g_x_0_x_0_0_y_0_zz = buffer_1010_spsd[11];

    auto g_x_0_x_0_0_z_0_xx = buffer_1010_spsd[12];

    auto g_x_0_x_0_0_z_0_xy = buffer_1010_spsd[13];

    auto g_x_0_x_0_0_z_0_xz = buffer_1010_spsd[14];

    auto g_x_0_x_0_0_z_0_yy = buffer_1010_spsd[15];

    auto g_x_0_x_0_0_z_0_yz = buffer_1010_spsd[16];

    auto g_x_0_x_0_0_z_0_zz = buffer_1010_spsd[17];

    auto g_x_0_y_0_0_x_0_xx = buffer_1010_spsd[18];

    auto g_x_0_y_0_0_x_0_xy = buffer_1010_spsd[19];

    auto g_x_0_y_0_0_x_0_xz = buffer_1010_spsd[20];

    auto g_x_0_y_0_0_x_0_yy = buffer_1010_spsd[21];

    auto g_x_0_y_0_0_x_0_yz = buffer_1010_spsd[22];

    auto g_x_0_y_0_0_x_0_zz = buffer_1010_spsd[23];

    auto g_x_0_y_0_0_y_0_xx = buffer_1010_spsd[24];

    auto g_x_0_y_0_0_y_0_xy = buffer_1010_spsd[25];

    auto g_x_0_y_0_0_y_0_xz = buffer_1010_spsd[26];

    auto g_x_0_y_0_0_y_0_yy = buffer_1010_spsd[27];

    auto g_x_0_y_0_0_y_0_yz = buffer_1010_spsd[28];

    auto g_x_0_y_0_0_y_0_zz = buffer_1010_spsd[29];

    auto g_x_0_y_0_0_z_0_xx = buffer_1010_spsd[30];

    auto g_x_0_y_0_0_z_0_xy = buffer_1010_spsd[31];

    auto g_x_0_y_0_0_z_0_xz = buffer_1010_spsd[32];

    auto g_x_0_y_0_0_z_0_yy = buffer_1010_spsd[33];

    auto g_x_0_y_0_0_z_0_yz = buffer_1010_spsd[34];

    auto g_x_0_y_0_0_z_0_zz = buffer_1010_spsd[35];

    auto g_x_0_z_0_0_x_0_xx = buffer_1010_spsd[36];

    auto g_x_0_z_0_0_x_0_xy = buffer_1010_spsd[37];

    auto g_x_0_z_0_0_x_0_xz = buffer_1010_spsd[38];

    auto g_x_0_z_0_0_x_0_yy = buffer_1010_spsd[39];

    auto g_x_0_z_0_0_x_0_yz = buffer_1010_spsd[40];

    auto g_x_0_z_0_0_x_0_zz = buffer_1010_spsd[41];

    auto g_x_0_z_0_0_y_0_xx = buffer_1010_spsd[42];

    auto g_x_0_z_0_0_y_0_xy = buffer_1010_spsd[43];

    auto g_x_0_z_0_0_y_0_xz = buffer_1010_spsd[44];

    auto g_x_0_z_0_0_y_0_yy = buffer_1010_spsd[45];

    auto g_x_0_z_0_0_y_0_yz = buffer_1010_spsd[46];

    auto g_x_0_z_0_0_y_0_zz = buffer_1010_spsd[47];

    auto g_x_0_z_0_0_z_0_xx = buffer_1010_spsd[48];

    auto g_x_0_z_0_0_z_0_xy = buffer_1010_spsd[49];

    auto g_x_0_z_0_0_z_0_xz = buffer_1010_spsd[50];

    auto g_x_0_z_0_0_z_0_yy = buffer_1010_spsd[51];

    auto g_x_0_z_0_0_z_0_yz = buffer_1010_spsd[52];

    auto g_x_0_z_0_0_z_0_zz = buffer_1010_spsd[53];

    auto g_y_0_x_0_0_x_0_xx = buffer_1010_spsd[54];

    auto g_y_0_x_0_0_x_0_xy = buffer_1010_spsd[55];

    auto g_y_0_x_0_0_x_0_xz = buffer_1010_spsd[56];

    auto g_y_0_x_0_0_x_0_yy = buffer_1010_spsd[57];

    auto g_y_0_x_0_0_x_0_yz = buffer_1010_spsd[58];

    auto g_y_0_x_0_0_x_0_zz = buffer_1010_spsd[59];

    auto g_y_0_x_0_0_y_0_xx = buffer_1010_spsd[60];

    auto g_y_0_x_0_0_y_0_xy = buffer_1010_spsd[61];

    auto g_y_0_x_0_0_y_0_xz = buffer_1010_spsd[62];

    auto g_y_0_x_0_0_y_0_yy = buffer_1010_spsd[63];

    auto g_y_0_x_0_0_y_0_yz = buffer_1010_spsd[64];

    auto g_y_0_x_0_0_y_0_zz = buffer_1010_spsd[65];

    auto g_y_0_x_0_0_z_0_xx = buffer_1010_spsd[66];

    auto g_y_0_x_0_0_z_0_xy = buffer_1010_spsd[67];

    auto g_y_0_x_0_0_z_0_xz = buffer_1010_spsd[68];

    auto g_y_0_x_0_0_z_0_yy = buffer_1010_spsd[69];

    auto g_y_0_x_0_0_z_0_yz = buffer_1010_spsd[70];

    auto g_y_0_x_0_0_z_0_zz = buffer_1010_spsd[71];

    auto g_y_0_y_0_0_x_0_xx = buffer_1010_spsd[72];

    auto g_y_0_y_0_0_x_0_xy = buffer_1010_spsd[73];

    auto g_y_0_y_0_0_x_0_xz = buffer_1010_spsd[74];

    auto g_y_0_y_0_0_x_0_yy = buffer_1010_spsd[75];

    auto g_y_0_y_0_0_x_0_yz = buffer_1010_spsd[76];

    auto g_y_0_y_0_0_x_0_zz = buffer_1010_spsd[77];

    auto g_y_0_y_0_0_y_0_xx = buffer_1010_spsd[78];

    auto g_y_0_y_0_0_y_0_xy = buffer_1010_spsd[79];

    auto g_y_0_y_0_0_y_0_xz = buffer_1010_spsd[80];

    auto g_y_0_y_0_0_y_0_yy = buffer_1010_spsd[81];

    auto g_y_0_y_0_0_y_0_yz = buffer_1010_spsd[82];

    auto g_y_0_y_0_0_y_0_zz = buffer_1010_spsd[83];

    auto g_y_0_y_0_0_z_0_xx = buffer_1010_spsd[84];

    auto g_y_0_y_0_0_z_0_xy = buffer_1010_spsd[85];

    auto g_y_0_y_0_0_z_0_xz = buffer_1010_spsd[86];

    auto g_y_0_y_0_0_z_0_yy = buffer_1010_spsd[87];

    auto g_y_0_y_0_0_z_0_yz = buffer_1010_spsd[88];

    auto g_y_0_y_0_0_z_0_zz = buffer_1010_spsd[89];

    auto g_y_0_z_0_0_x_0_xx = buffer_1010_spsd[90];

    auto g_y_0_z_0_0_x_0_xy = buffer_1010_spsd[91];

    auto g_y_0_z_0_0_x_0_xz = buffer_1010_spsd[92];

    auto g_y_0_z_0_0_x_0_yy = buffer_1010_spsd[93];

    auto g_y_0_z_0_0_x_0_yz = buffer_1010_spsd[94];

    auto g_y_0_z_0_0_x_0_zz = buffer_1010_spsd[95];

    auto g_y_0_z_0_0_y_0_xx = buffer_1010_spsd[96];

    auto g_y_0_z_0_0_y_0_xy = buffer_1010_spsd[97];

    auto g_y_0_z_0_0_y_0_xz = buffer_1010_spsd[98];

    auto g_y_0_z_0_0_y_0_yy = buffer_1010_spsd[99];

    auto g_y_0_z_0_0_y_0_yz = buffer_1010_spsd[100];

    auto g_y_0_z_0_0_y_0_zz = buffer_1010_spsd[101];

    auto g_y_0_z_0_0_z_0_xx = buffer_1010_spsd[102];

    auto g_y_0_z_0_0_z_0_xy = buffer_1010_spsd[103];

    auto g_y_0_z_0_0_z_0_xz = buffer_1010_spsd[104];

    auto g_y_0_z_0_0_z_0_yy = buffer_1010_spsd[105];

    auto g_y_0_z_0_0_z_0_yz = buffer_1010_spsd[106];

    auto g_y_0_z_0_0_z_0_zz = buffer_1010_spsd[107];

    auto g_z_0_x_0_0_x_0_xx = buffer_1010_spsd[108];

    auto g_z_0_x_0_0_x_0_xy = buffer_1010_spsd[109];

    auto g_z_0_x_0_0_x_0_xz = buffer_1010_spsd[110];

    auto g_z_0_x_0_0_x_0_yy = buffer_1010_spsd[111];

    auto g_z_0_x_0_0_x_0_yz = buffer_1010_spsd[112];

    auto g_z_0_x_0_0_x_0_zz = buffer_1010_spsd[113];

    auto g_z_0_x_0_0_y_0_xx = buffer_1010_spsd[114];

    auto g_z_0_x_0_0_y_0_xy = buffer_1010_spsd[115];

    auto g_z_0_x_0_0_y_0_xz = buffer_1010_spsd[116];

    auto g_z_0_x_0_0_y_0_yy = buffer_1010_spsd[117];

    auto g_z_0_x_0_0_y_0_yz = buffer_1010_spsd[118];

    auto g_z_0_x_0_0_y_0_zz = buffer_1010_spsd[119];

    auto g_z_0_x_0_0_z_0_xx = buffer_1010_spsd[120];

    auto g_z_0_x_0_0_z_0_xy = buffer_1010_spsd[121];

    auto g_z_0_x_0_0_z_0_xz = buffer_1010_spsd[122];

    auto g_z_0_x_0_0_z_0_yy = buffer_1010_spsd[123];

    auto g_z_0_x_0_0_z_0_yz = buffer_1010_spsd[124];

    auto g_z_0_x_0_0_z_0_zz = buffer_1010_spsd[125];

    auto g_z_0_y_0_0_x_0_xx = buffer_1010_spsd[126];

    auto g_z_0_y_0_0_x_0_xy = buffer_1010_spsd[127];

    auto g_z_0_y_0_0_x_0_xz = buffer_1010_spsd[128];

    auto g_z_0_y_0_0_x_0_yy = buffer_1010_spsd[129];

    auto g_z_0_y_0_0_x_0_yz = buffer_1010_spsd[130];

    auto g_z_0_y_0_0_x_0_zz = buffer_1010_spsd[131];

    auto g_z_0_y_0_0_y_0_xx = buffer_1010_spsd[132];

    auto g_z_0_y_0_0_y_0_xy = buffer_1010_spsd[133];

    auto g_z_0_y_0_0_y_0_xz = buffer_1010_spsd[134];

    auto g_z_0_y_0_0_y_0_yy = buffer_1010_spsd[135];

    auto g_z_0_y_0_0_y_0_yz = buffer_1010_spsd[136];

    auto g_z_0_y_0_0_y_0_zz = buffer_1010_spsd[137];

    auto g_z_0_y_0_0_z_0_xx = buffer_1010_spsd[138];

    auto g_z_0_y_0_0_z_0_xy = buffer_1010_spsd[139];

    auto g_z_0_y_0_0_z_0_xz = buffer_1010_spsd[140];

    auto g_z_0_y_0_0_z_0_yy = buffer_1010_spsd[141];

    auto g_z_0_y_0_0_z_0_yz = buffer_1010_spsd[142];

    auto g_z_0_y_0_0_z_0_zz = buffer_1010_spsd[143];

    auto g_z_0_z_0_0_x_0_xx = buffer_1010_spsd[144];

    auto g_z_0_z_0_0_x_0_xy = buffer_1010_spsd[145];

    auto g_z_0_z_0_0_x_0_xz = buffer_1010_spsd[146];

    auto g_z_0_z_0_0_x_0_yy = buffer_1010_spsd[147];

    auto g_z_0_z_0_0_x_0_yz = buffer_1010_spsd[148];

    auto g_z_0_z_0_0_x_0_zz = buffer_1010_spsd[149];

    auto g_z_0_z_0_0_y_0_xx = buffer_1010_spsd[150];

    auto g_z_0_z_0_0_y_0_xy = buffer_1010_spsd[151];

    auto g_z_0_z_0_0_y_0_xz = buffer_1010_spsd[152];

    auto g_z_0_z_0_0_y_0_yy = buffer_1010_spsd[153];

    auto g_z_0_z_0_0_y_0_yz = buffer_1010_spsd[154];

    auto g_z_0_z_0_0_y_0_zz = buffer_1010_spsd[155];

    auto g_z_0_z_0_0_z_0_xx = buffer_1010_spsd[156];

    auto g_z_0_z_0_0_z_0_xy = buffer_1010_spsd[157];

    auto g_z_0_z_0_0_z_0_xz = buffer_1010_spsd[158];

    auto g_z_0_z_0_0_z_0_yy = buffer_1010_spsd[159];

    auto g_z_0_z_0_0_z_0_yz = buffer_1010_spsd[160];

    auto g_z_0_z_0_0_z_0_zz = buffer_1010_spsd[161];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_0_x_0_0_x_0_xx, g_x_0_x_0_0_x_0_xy, g_x_0_x_0_0_x_0_xz, g_x_0_x_0_0_x_0_yy, g_x_0_x_0_0_x_0_yz, g_x_0_x_0_0_x_0_zz, g_x_x_x_xx, g_x_x_x_xy, g_x_x_x_xz, g_x_x_x_yy, g_x_x_x_yz, g_x_x_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_x_0_xx[i] = 4.0 * g_x_x_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_0_xy[i] = 4.0 * g_x_x_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_0_xz[i] = 4.0 * g_x_x_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_0_yy[i] = 4.0 * g_x_x_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_0_yz[i] = 4.0 * g_x_x_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_0_zz[i] = 4.0 * g_x_x_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_0_x_0_0_y_0_xx, g_x_0_x_0_0_y_0_xy, g_x_0_x_0_0_y_0_xz, g_x_0_x_0_0_y_0_yy, g_x_0_x_0_0_y_0_yz, g_x_0_x_0_0_y_0_zz, g_x_y_x_xx, g_x_y_x_xy, g_x_y_x_xz, g_x_y_x_yy, g_x_y_x_yz, g_x_y_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_y_0_xx[i] = 4.0 * g_x_y_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_0_xy[i] = 4.0 * g_x_y_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_0_xz[i] = 4.0 * g_x_y_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_0_yy[i] = 4.0 * g_x_y_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_0_yz[i] = 4.0 * g_x_y_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_0_zz[i] = 4.0 * g_x_y_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_0_x_0_0_z_0_xx, g_x_0_x_0_0_z_0_xy, g_x_0_x_0_0_z_0_xz, g_x_0_x_0_0_z_0_yy, g_x_0_x_0_0_z_0_yz, g_x_0_x_0_0_z_0_zz, g_x_z_x_xx, g_x_z_x_xy, g_x_z_x_xz, g_x_z_x_yy, g_x_z_x_yz, g_x_z_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_z_0_xx[i] = 4.0 * g_x_z_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_0_xy[i] = 4.0 * g_x_z_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_0_xz[i] = 4.0 * g_x_z_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_0_yy[i] = 4.0 * g_x_z_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_0_yz[i] = 4.0 * g_x_z_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_0_zz[i] = 4.0 * g_x_z_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_x_0_y_0_0_x_0_xx, g_x_0_y_0_0_x_0_xy, g_x_0_y_0_0_x_0_xz, g_x_0_y_0_0_x_0_yy, g_x_0_y_0_0_x_0_yz, g_x_0_y_0_0_x_0_zz, g_x_x_y_xx, g_x_x_y_xy, g_x_x_y_xz, g_x_x_y_yy, g_x_x_y_yz, g_x_x_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_x_0_xx[i] = 4.0 * g_x_x_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_0_xy[i] = 4.0 * g_x_x_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_0_xz[i] = 4.0 * g_x_x_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_0_yy[i] = 4.0 * g_x_x_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_0_yz[i] = 4.0 * g_x_x_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_0_zz[i] = 4.0 * g_x_x_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_x_0_y_0_0_y_0_xx, g_x_0_y_0_0_y_0_xy, g_x_0_y_0_0_y_0_xz, g_x_0_y_0_0_y_0_yy, g_x_0_y_0_0_y_0_yz, g_x_0_y_0_0_y_0_zz, g_x_y_y_xx, g_x_y_y_xy, g_x_y_y_xz, g_x_y_y_yy, g_x_y_y_yz, g_x_y_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_y_0_xx[i] = 4.0 * g_x_y_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_0_xy[i] = 4.0 * g_x_y_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_0_xz[i] = 4.0 * g_x_y_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_0_yy[i] = 4.0 * g_x_y_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_0_yz[i] = 4.0 * g_x_y_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_0_zz[i] = 4.0 * g_x_y_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_x_0_y_0_0_z_0_xx, g_x_0_y_0_0_z_0_xy, g_x_0_y_0_0_z_0_xz, g_x_0_y_0_0_z_0_yy, g_x_0_y_0_0_z_0_yz, g_x_0_y_0_0_z_0_zz, g_x_z_y_xx, g_x_z_y_xy, g_x_z_y_xz, g_x_z_y_yy, g_x_z_y_yz, g_x_z_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_z_0_xx[i] = 4.0 * g_x_z_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_0_xy[i] = 4.0 * g_x_z_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_0_xz[i] = 4.0 * g_x_z_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_0_yy[i] = 4.0 * g_x_z_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_0_yz[i] = 4.0 * g_x_z_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_0_zz[i] = 4.0 * g_x_z_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_x_0_z_0_0_x_0_xx, g_x_0_z_0_0_x_0_xy, g_x_0_z_0_0_x_0_xz, g_x_0_z_0_0_x_0_yy, g_x_0_z_0_0_x_0_yz, g_x_0_z_0_0_x_0_zz, g_x_x_z_xx, g_x_x_z_xy, g_x_x_z_xz, g_x_x_z_yy, g_x_x_z_yz, g_x_x_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_x_0_xx[i] = 4.0 * g_x_x_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_0_xy[i] = 4.0 * g_x_x_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_0_xz[i] = 4.0 * g_x_x_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_0_yy[i] = 4.0 * g_x_x_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_0_yz[i] = 4.0 * g_x_x_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_0_zz[i] = 4.0 * g_x_x_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_x_0_z_0_0_y_0_xx, g_x_0_z_0_0_y_0_xy, g_x_0_z_0_0_y_0_xz, g_x_0_z_0_0_y_0_yy, g_x_0_z_0_0_y_0_yz, g_x_0_z_0_0_y_0_zz, g_x_y_z_xx, g_x_y_z_xy, g_x_y_z_xz, g_x_y_z_yy, g_x_y_z_yz, g_x_y_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_y_0_xx[i] = 4.0 * g_x_y_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_0_xy[i] = 4.0 * g_x_y_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_0_xz[i] = 4.0 * g_x_y_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_0_yy[i] = 4.0 * g_x_y_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_0_yz[i] = 4.0 * g_x_y_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_0_zz[i] = 4.0 * g_x_y_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_x_0_z_0_0_z_0_xx, g_x_0_z_0_0_z_0_xy, g_x_0_z_0_0_z_0_xz, g_x_0_z_0_0_z_0_yy, g_x_0_z_0_0_z_0_yz, g_x_0_z_0_0_z_0_zz, g_x_z_z_xx, g_x_z_z_xy, g_x_z_z_xz, g_x_z_z_yy, g_x_z_z_yz, g_x_z_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_z_0_xx[i] = 4.0 * g_x_z_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_0_xy[i] = 4.0 * g_x_z_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_0_xz[i] = 4.0 * g_x_z_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_0_yy[i] = 4.0 * g_x_z_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_0_yz[i] = 4.0 * g_x_z_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_0_zz[i] = 4.0 * g_x_z_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_y_0_x_0_0_x_0_xx, g_y_0_x_0_0_x_0_xy, g_y_0_x_0_0_x_0_xz, g_y_0_x_0_0_x_0_yy, g_y_0_x_0_0_x_0_yz, g_y_0_x_0_0_x_0_zz, g_y_x_x_xx, g_y_x_x_xy, g_y_x_x_xz, g_y_x_x_yy, g_y_x_x_yz, g_y_x_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_x_0_xx[i] = 4.0 * g_y_x_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_0_xy[i] = 4.0 * g_y_x_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_0_xz[i] = 4.0 * g_y_x_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_0_yy[i] = 4.0 * g_y_x_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_0_yz[i] = 4.0 * g_y_x_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_0_zz[i] = 4.0 * g_y_x_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_y_0_x_0_0_y_0_xx, g_y_0_x_0_0_y_0_xy, g_y_0_x_0_0_y_0_xz, g_y_0_x_0_0_y_0_yy, g_y_0_x_0_0_y_0_yz, g_y_0_x_0_0_y_0_zz, g_y_y_x_xx, g_y_y_x_xy, g_y_y_x_xz, g_y_y_x_yy, g_y_y_x_yz, g_y_y_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_y_0_xx[i] = 4.0 * g_y_y_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_0_xy[i] = 4.0 * g_y_y_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_0_xz[i] = 4.0 * g_y_y_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_0_yy[i] = 4.0 * g_y_y_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_0_yz[i] = 4.0 * g_y_y_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_0_zz[i] = 4.0 * g_y_y_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_y_0_x_0_0_z_0_xx, g_y_0_x_0_0_z_0_xy, g_y_0_x_0_0_z_0_xz, g_y_0_x_0_0_z_0_yy, g_y_0_x_0_0_z_0_yz, g_y_0_x_0_0_z_0_zz, g_y_z_x_xx, g_y_z_x_xy, g_y_z_x_xz, g_y_z_x_yy, g_y_z_x_yz, g_y_z_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_z_0_xx[i] = 4.0 * g_y_z_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_0_xy[i] = 4.0 * g_y_z_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_0_xz[i] = 4.0 * g_y_z_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_0_yy[i] = 4.0 * g_y_z_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_0_yz[i] = 4.0 * g_y_z_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_0_zz[i] = 4.0 * g_y_z_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_y_0_y_0_0_x_0_xx, g_y_0_y_0_0_x_0_xy, g_y_0_y_0_0_x_0_xz, g_y_0_y_0_0_x_0_yy, g_y_0_y_0_0_x_0_yz, g_y_0_y_0_0_x_0_zz, g_y_x_y_xx, g_y_x_y_xy, g_y_x_y_xz, g_y_x_y_yy, g_y_x_y_yz, g_y_x_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_x_0_xx[i] = 4.0 * g_y_x_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_0_xy[i] = 4.0 * g_y_x_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_0_xz[i] = 4.0 * g_y_x_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_0_yy[i] = 4.0 * g_y_x_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_0_yz[i] = 4.0 * g_y_x_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_0_zz[i] = 4.0 * g_y_x_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_y_0_y_0_0_y_0_xx, g_y_0_y_0_0_y_0_xy, g_y_0_y_0_0_y_0_xz, g_y_0_y_0_0_y_0_yy, g_y_0_y_0_0_y_0_yz, g_y_0_y_0_0_y_0_zz, g_y_y_y_xx, g_y_y_y_xy, g_y_y_y_xz, g_y_y_y_yy, g_y_y_y_yz, g_y_y_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_y_0_xx[i] = 4.0 * g_y_y_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_0_xy[i] = 4.0 * g_y_y_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_0_xz[i] = 4.0 * g_y_y_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_0_yy[i] = 4.0 * g_y_y_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_0_yz[i] = 4.0 * g_y_y_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_0_zz[i] = 4.0 * g_y_y_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_y_0_y_0_0_z_0_xx, g_y_0_y_0_0_z_0_xy, g_y_0_y_0_0_z_0_xz, g_y_0_y_0_0_z_0_yy, g_y_0_y_0_0_z_0_yz, g_y_0_y_0_0_z_0_zz, g_y_z_y_xx, g_y_z_y_xy, g_y_z_y_xz, g_y_z_y_yy, g_y_z_y_yz, g_y_z_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_z_0_xx[i] = 4.0 * g_y_z_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_0_xy[i] = 4.0 * g_y_z_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_0_xz[i] = 4.0 * g_y_z_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_0_yy[i] = 4.0 * g_y_z_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_0_yz[i] = 4.0 * g_y_z_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_0_zz[i] = 4.0 * g_y_z_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_y_0_z_0_0_x_0_xx, g_y_0_z_0_0_x_0_xy, g_y_0_z_0_0_x_0_xz, g_y_0_z_0_0_x_0_yy, g_y_0_z_0_0_x_0_yz, g_y_0_z_0_0_x_0_zz, g_y_x_z_xx, g_y_x_z_xy, g_y_x_z_xz, g_y_x_z_yy, g_y_x_z_yz, g_y_x_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_x_0_xx[i] = 4.0 * g_y_x_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_0_xy[i] = 4.0 * g_y_x_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_0_xz[i] = 4.0 * g_y_x_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_0_yy[i] = 4.0 * g_y_x_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_0_yz[i] = 4.0 * g_y_x_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_0_zz[i] = 4.0 * g_y_x_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_y_0_z_0_0_y_0_xx, g_y_0_z_0_0_y_0_xy, g_y_0_z_0_0_y_0_xz, g_y_0_z_0_0_y_0_yy, g_y_0_z_0_0_y_0_yz, g_y_0_z_0_0_y_0_zz, g_y_y_z_xx, g_y_y_z_xy, g_y_y_z_xz, g_y_y_z_yy, g_y_y_z_yz, g_y_y_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_y_0_xx[i] = 4.0 * g_y_y_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_0_xy[i] = 4.0 * g_y_y_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_0_xz[i] = 4.0 * g_y_y_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_0_yy[i] = 4.0 * g_y_y_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_0_yz[i] = 4.0 * g_y_y_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_0_zz[i] = 4.0 * g_y_y_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_y_0_z_0_0_z_0_xx, g_y_0_z_0_0_z_0_xy, g_y_0_z_0_0_z_0_xz, g_y_0_z_0_0_z_0_yy, g_y_0_z_0_0_z_0_yz, g_y_0_z_0_0_z_0_zz, g_y_z_z_xx, g_y_z_z_xy, g_y_z_z_xz, g_y_z_z_yy, g_y_z_z_yz, g_y_z_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_z_0_xx[i] = 4.0 * g_y_z_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_0_xy[i] = 4.0 * g_y_z_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_0_xz[i] = 4.0 * g_y_z_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_0_yy[i] = 4.0 * g_y_z_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_0_yz[i] = 4.0 * g_y_z_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_0_zz[i] = 4.0 * g_y_z_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_z_0_x_0_0_x_0_xx, g_z_0_x_0_0_x_0_xy, g_z_0_x_0_0_x_0_xz, g_z_0_x_0_0_x_0_yy, g_z_0_x_0_0_x_0_yz, g_z_0_x_0_0_x_0_zz, g_z_x_x_xx, g_z_x_x_xy, g_z_x_x_xz, g_z_x_x_yy, g_z_x_x_yz, g_z_x_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_x_0_xx[i] = 4.0 * g_z_x_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_0_xy[i] = 4.0 * g_z_x_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_0_xz[i] = 4.0 * g_z_x_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_0_yy[i] = 4.0 * g_z_x_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_0_yz[i] = 4.0 * g_z_x_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_0_zz[i] = 4.0 * g_z_x_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_z_0_x_0_0_y_0_xx, g_z_0_x_0_0_y_0_xy, g_z_0_x_0_0_y_0_xz, g_z_0_x_0_0_y_0_yy, g_z_0_x_0_0_y_0_yz, g_z_0_x_0_0_y_0_zz, g_z_y_x_xx, g_z_y_x_xy, g_z_y_x_xz, g_z_y_x_yy, g_z_y_x_yz, g_z_y_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_y_0_xx[i] = 4.0 * g_z_y_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_0_xy[i] = 4.0 * g_z_y_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_0_xz[i] = 4.0 * g_z_y_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_0_yy[i] = 4.0 * g_z_y_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_0_yz[i] = 4.0 * g_z_y_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_0_zz[i] = 4.0 * g_z_y_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_z_0_x_0_0_z_0_xx, g_z_0_x_0_0_z_0_xy, g_z_0_x_0_0_z_0_xz, g_z_0_x_0_0_z_0_yy, g_z_0_x_0_0_z_0_yz, g_z_0_x_0_0_z_0_zz, g_z_z_x_xx, g_z_z_x_xy, g_z_z_x_xz, g_z_z_x_yy, g_z_z_x_yz, g_z_z_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_z_0_xx[i] = 4.0 * g_z_z_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_0_xy[i] = 4.0 * g_z_z_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_0_xz[i] = 4.0 * g_z_z_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_0_yy[i] = 4.0 * g_z_z_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_0_yz[i] = 4.0 * g_z_z_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_0_zz[i] = 4.0 * g_z_z_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_z_0_y_0_0_x_0_xx, g_z_0_y_0_0_x_0_xy, g_z_0_y_0_0_x_0_xz, g_z_0_y_0_0_x_0_yy, g_z_0_y_0_0_x_0_yz, g_z_0_y_0_0_x_0_zz, g_z_x_y_xx, g_z_x_y_xy, g_z_x_y_xz, g_z_x_y_yy, g_z_x_y_yz, g_z_x_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_x_0_xx[i] = 4.0 * g_z_x_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_0_xy[i] = 4.0 * g_z_x_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_0_xz[i] = 4.0 * g_z_x_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_0_yy[i] = 4.0 * g_z_x_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_0_yz[i] = 4.0 * g_z_x_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_0_zz[i] = 4.0 * g_z_x_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_z_0_y_0_0_y_0_xx, g_z_0_y_0_0_y_0_xy, g_z_0_y_0_0_y_0_xz, g_z_0_y_0_0_y_0_yy, g_z_0_y_0_0_y_0_yz, g_z_0_y_0_0_y_0_zz, g_z_y_y_xx, g_z_y_y_xy, g_z_y_y_xz, g_z_y_y_yy, g_z_y_y_yz, g_z_y_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_y_0_xx[i] = 4.0 * g_z_y_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_0_xy[i] = 4.0 * g_z_y_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_0_xz[i] = 4.0 * g_z_y_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_0_yy[i] = 4.0 * g_z_y_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_0_yz[i] = 4.0 * g_z_y_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_0_zz[i] = 4.0 * g_z_y_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_z_0_y_0_0_z_0_xx, g_z_0_y_0_0_z_0_xy, g_z_0_y_0_0_z_0_xz, g_z_0_y_0_0_z_0_yy, g_z_0_y_0_0_z_0_yz, g_z_0_y_0_0_z_0_zz, g_z_z_y_xx, g_z_z_y_xy, g_z_z_y_xz, g_z_z_y_yy, g_z_z_y_yz, g_z_z_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_z_0_xx[i] = 4.0 * g_z_z_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_0_xy[i] = 4.0 * g_z_z_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_0_xz[i] = 4.0 * g_z_z_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_0_yy[i] = 4.0 * g_z_z_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_0_yz[i] = 4.0 * g_z_z_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_0_zz[i] = 4.0 * g_z_z_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_z_0_z_0_0_x_0_xx, g_z_0_z_0_0_x_0_xy, g_z_0_z_0_0_x_0_xz, g_z_0_z_0_0_x_0_yy, g_z_0_z_0_0_x_0_yz, g_z_0_z_0_0_x_0_zz, g_z_x_z_xx, g_z_x_z_xy, g_z_x_z_xz, g_z_x_z_yy, g_z_x_z_yz, g_z_x_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_x_0_xx[i] = 4.0 * g_z_x_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_0_xy[i] = 4.0 * g_z_x_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_0_xz[i] = 4.0 * g_z_x_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_0_yy[i] = 4.0 * g_z_x_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_0_yz[i] = 4.0 * g_z_x_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_0_zz[i] = 4.0 * g_z_x_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_z_0_z_0_0_y_0_xx, g_z_0_z_0_0_y_0_xy, g_z_0_z_0_0_y_0_xz, g_z_0_z_0_0_y_0_yy, g_z_0_z_0_0_y_0_yz, g_z_0_z_0_0_y_0_zz, g_z_y_z_xx, g_z_y_z_xy, g_z_y_z_xz, g_z_y_z_yy, g_z_y_z_yz, g_z_y_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_y_0_xx[i] = 4.0 * g_z_y_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_0_xy[i] = 4.0 * g_z_y_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_0_xz[i] = 4.0 * g_z_y_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_0_yy[i] = 4.0 * g_z_y_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_0_yz[i] = 4.0 * g_z_y_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_0_zz[i] = 4.0 * g_z_y_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_z_0_z_0_0_z_0_xx, g_z_0_z_0_0_z_0_xy, g_z_0_z_0_0_z_0_xz, g_z_0_z_0_0_z_0_yy, g_z_0_z_0_0_z_0_yz, g_z_0_z_0_0_z_0_zz, g_z_z_z_xx, g_z_z_z_xy, g_z_z_z_xz, g_z_z_z_yy, g_z_z_z_yz, g_z_z_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_z_0_xx[i] = 4.0 * g_z_z_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_0_xy[i] = 4.0 * g_z_z_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_0_xz[i] = 4.0 * g_z_z_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_0_yy[i] = 4.0 * g_z_z_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_0_yz[i] = 4.0 * g_z_z_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_0_zz[i] = 4.0 * g_z_z_z_zz[i] * a_exp * c_exps[i];
    }
}

} // t4c_geom namespace


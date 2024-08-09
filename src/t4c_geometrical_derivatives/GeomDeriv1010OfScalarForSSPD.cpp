#include "GeomDeriv1010OfScalarForSSPD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1010_sspd_0(CSimdArray<double>& buffer_1010_sspd,
                     const CSimdArray<double>& buffer_pssd,
                     const CSimdArray<double>& buffer_psdd,
                     const double a_exp,
                     const double* c_exps) -> void
{
    const auto ndims = buffer_1010_sspd.number_of_columns();

    /// Set up components of auxilary buffer : buffer_pssd

    auto g_x_0_0_xx = buffer_pssd[0];

    auto g_x_0_0_xy = buffer_pssd[1];

    auto g_x_0_0_xz = buffer_pssd[2];

    auto g_x_0_0_yy = buffer_pssd[3];

    auto g_x_0_0_yz = buffer_pssd[4];

    auto g_x_0_0_zz = buffer_pssd[5];

    auto g_y_0_0_xx = buffer_pssd[6];

    auto g_y_0_0_xy = buffer_pssd[7];

    auto g_y_0_0_xz = buffer_pssd[8];

    auto g_y_0_0_yy = buffer_pssd[9];

    auto g_y_0_0_yz = buffer_pssd[10];

    auto g_y_0_0_zz = buffer_pssd[11];

    auto g_z_0_0_xx = buffer_pssd[12];

    auto g_z_0_0_xy = buffer_pssd[13];

    auto g_z_0_0_xz = buffer_pssd[14];

    auto g_z_0_0_yy = buffer_pssd[15];

    auto g_z_0_0_yz = buffer_pssd[16];

    auto g_z_0_0_zz = buffer_pssd[17];

    /// Set up components of auxilary buffer : buffer_psdd

    auto g_x_0_xx_xx = buffer_psdd[0];

    auto g_x_0_xx_xy = buffer_psdd[1];

    auto g_x_0_xx_xz = buffer_psdd[2];

    auto g_x_0_xx_yy = buffer_psdd[3];

    auto g_x_0_xx_yz = buffer_psdd[4];

    auto g_x_0_xx_zz = buffer_psdd[5];

    auto g_x_0_xy_xx = buffer_psdd[6];

    auto g_x_0_xy_xy = buffer_psdd[7];

    auto g_x_0_xy_xz = buffer_psdd[8];

    auto g_x_0_xy_yy = buffer_psdd[9];

    auto g_x_0_xy_yz = buffer_psdd[10];

    auto g_x_0_xy_zz = buffer_psdd[11];

    auto g_x_0_xz_xx = buffer_psdd[12];

    auto g_x_0_xz_xy = buffer_psdd[13];

    auto g_x_0_xz_xz = buffer_psdd[14];

    auto g_x_0_xz_yy = buffer_psdd[15];

    auto g_x_0_xz_yz = buffer_psdd[16];

    auto g_x_0_xz_zz = buffer_psdd[17];

    auto g_x_0_yy_xx = buffer_psdd[18];

    auto g_x_0_yy_xy = buffer_psdd[19];

    auto g_x_0_yy_xz = buffer_psdd[20];

    auto g_x_0_yy_yy = buffer_psdd[21];

    auto g_x_0_yy_yz = buffer_psdd[22];

    auto g_x_0_yy_zz = buffer_psdd[23];

    auto g_x_0_yz_xx = buffer_psdd[24];

    auto g_x_0_yz_xy = buffer_psdd[25];

    auto g_x_0_yz_xz = buffer_psdd[26];

    auto g_x_0_yz_yy = buffer_psdd[27];

    auto g_x_0_yz_yz = buffer_psdd[28];

    auto g_x_0_yz_zz = buffer_psdd[29];

    auto g_x_0_zz_xx = buffer_psdd[30];

    auto g_x_0_zz_xy = buffer_psdd[31];

    auto g_x_0_zz_xz = buffer_psdd[32];

    auto g_x_0_zz_yy = buffer_psdd[33];

    auto g_x_0_zz_yz = buffer_psdd[34];

    auto g_x_0_zz_zz = buffer_psdd[35];

    auto g_y_0_xx_xx = buffer_psdd[36];

    auto g_y_0_xx_xy = buffer_psdd[37];

    auto g_y_0_xx_xz = buffer_psdd[38];

    auto g_y_0_xx_yy = buffer_psdd[39];

    auto g_y_0_xx_yz = buffer_psdd[40];

    auto g_y_0_xx_zz = buffer_psdd[41];

    auto g_y_0_xy_xx = buffer_psdd[42];

    auto g_y_0_xy_xy = buffer_psdd[43];

    auto g_y_0_xy_xz = buffer_psdd[44];

    auto g_y_0_xy_yy = buffer_psdd[45];

    auto g_y_0_xy_yz = buffer_psdd[46];

    auto g_y_0_xy_zz = buffer_psdd[47];

    auto g_y_0_xz_xx = buffer_psdd[48];

    auto g_y_0_xz_xy = buffer_psdd[49];

    auto g_y_0_xz_xz = buffer_psdd[50];

    auto g_y_0_xz_yy = buffer_psdd[51];

    auto g_y_0_xz_yz = buffer_psdd[52];

    auto g_y_0_xz_zz = buffer_psdd[53];

    auto g_y_0_yy_xx = buffer_psdd[54];

    auto g_y_0_yy_xy = buffer_psdd[55];

    auto g_y_0_yy_xz = buffer_psdd[56];

    auto g_y_0_yy_yy = buffer_psdd[57];

    auto g_y_0_yy_yz = buffer_psdd[58];

    auto g_y_0_yy_zz = buffer_psdd[59];

    auto g_y_0_yz_xx = buffer_psdd[60];

    auto g_y_0_yz_xy = buffer_psdd[61];

    auto g_y_0_yz_xz = buffer_psdd[62];

    auto g_y_0_yz_yy = buffer_psdd[63];

    auto g_y_0_yz_yz = buffer_psdd[64];

    auto g_y_0_yz_zz = buffer_psdd[65];

    auto g_y_0_zz_xx = buffer_psdd[66];

    auto g_y_0_zz_xy = buffer_psdd[67];

    auto g_y_0_zz_xz = buffer_psdd[68];

    auto g_y_0_zz_yy = buffer_psdd[69];

    auto g_y_0_zz_yz = buffer_psdd[70];

    auto g_y_0_zz_zz = buffer_psdd[71];

    auto g_z_0_xx_xx = buffer_psdd[72];

    auto g_z_0_xx_xy = buffer_psdd[73];

    auto g_z_0_xx_xz = buffer_psdd[74];

    auto g_z_0_xx_yy = buffer_psdd[75];

    auto g_z_0_xx_yz = buffer_psdd[76];

    auto g_z_0_xx_zz = buffer_psdd[77];

    auto g_z_0_xy_xx = buffer_psdd[78];

    auto g_z_0_xy_xy = buffer_psdd[79];

    auto g_z_0_xy_xz = buffer_psdd[80];

    auto g_z_0_xy_yy = buffer_psdd[81];

    auto g_z_0_xy_yz = buffer_psdd[82];

    auto g_z_0_xy_zz = buffer_psdd[83];

    auto g_z_0_xz_xx = buffer_psdd[84];

    auto g_z_0_xz_xy = buffer_psdd[85];

    auto g_z_0_xz_xz = buffer_psdd[86];

    auto g_z_0_xz_yy = buffer_psdd[87];

    auto g_z_0_xz_yz = buffer_psdd[88];

    auto g_z_0_xz_zz = buffer_psdd[89];

    auto g_z_0_yy_xx = buffer_psdd[90];

    auto g_z_0_yy_xy = buffer_psdd[91];

    auto g_z_0_yy_xz = buffer_psdd[92];

    auto g_z_0_yy_yy = buffer_psdd[93];

    auto g_z_0_yy_yz = buffer_psdd[94];

    auto g_z_0_yy_zz = buffer_psdd[95];

    auto g_z_0_yz_xx = buffer_psdd[96];

    auto g_z_0_yz_xy = buffer_psdd[97];

    auto g_z_0_yz_xz = buffer_psdd[98];

    auto g_z_0_yz_yy = buffer_psdd[99];

    auto g_z_0_yz_yz = buffer_psdd[100];

    auto g_z_0_yz_zz = buffer_psdd[101];

    auto g_z_0_zz_xx = buffer_psdd[102];

    auto g_z_0_zz_xy = buffer_psdd[103];

    auto g_z_0_zz_xz = buffer_psdd[104];

    auto g_z_0_zz_yy = buffer_psdd[105];

    auto g_z_0_zz_yz = buffer_psdd[106];

    auto g_z_0_zz_zz = buffer_psdd[107];

    /// Set up components of integrals buffer : buffer_1010_sspd

    auto g_x_0_x_0_0_0_x_xx = buffer_1010_sspd[0];

    auto g_x_0_x_0_0_0_x_xy = buffer_1010_sspd[1];

    auto g_x_0_x_0_0_0_x_xz = buffer_1010_sspd[2];

    auto g_x_0_x_0_0_0_x_yy = buffer_1010_sspd[3];

    auto g_x_0_x_0_0_0_x_yz = buffer_1010_sspd[4];

    auto g_x_0_x_0_0_0_x_zz = buffer_1010_sspd[5];

    auto g_x_0_x_0_0_0_y_xx = buffer_1010_sspd[6];

    auto g_x_0_x_0_0_0_y_xy = buffer_1010_sspd[7];

    auto g_x_0_x_0_0_0_y_xz = buffer_1010_sspd[8];

    auto g_x_0_x_0_0_0_y_yy = buffer_1010_sspd[9];

    auto g_x_0_x_0_0_0_y_yz = buffer_1010_sspd[10];

    auto g_x_0_x_0_0_0_y_zz = buffer_1010_sspd[11];

    auto g_x_0_x_0_0_0_z_xx = buffer_1010_sspd[12];

    auto g_x_0_x_0_0_0_z_xy = buffer_1010_sspd[13];

    auto g_x_0_x_0_0_0_z_xz = buffer_1010_sspd[14];

    auto g_x_0_x_0_0_0_z_yy = buffer_1010_sspd[15];

    auto g_x_0_x_0_0_0_z_yz = buffer_1010_sspd[16];

    auto g_x_0_x_0_0_0_z_zz = buffer_1010_sspd[17];

    auto g_x_0_y_0_0_0_x_xx = buffer_1010_sspd[18];

    auto g_x_0_y_0_0_0_x_xy = buffer_1010_sspd[19];

    auto g_x_0_y_0_0_0_x_xz = buffer_1010_sspd[20];

    auto g_x_0_y_0_0_0_x_yy = buffer_1010_sspd[21];

    auto g_x_0_y_0_0_0_x_yz = buffer_1010_sspd[22];

    auto g_x_0_y_0_0_0_x_zz = buffer_1010_sspd[23];

    auto g_x_0_y_0_0_0_y_xx = buffer_1010_sspd[24];

    auto g_x_0_y_0_0_0_y_xy = buffer_1010_sspd[25];

    auto g_x_0_y_0_0_0_y_xz = buffer_1010_sspd[26];

    auto g_x_0_y_0_0_0_y_yy = buffer_1010_sspd[27];

    auto g_x_0_y_0_0_0_y_yz = buffer_1010_sspd[28];

    auto g_x_0_y_0_0_0_y_zz = buffer_1010_sspd[29];

    auto g_x_0_y_0_0_0_z_xx = buffer_1010_sspd[30];

    auto g_x_0_y_0_0_0_z_xy = buffer_1010_sspd[31];

    auto g_x_0_y_0_0_0_z_xz = buffer_1010_sspd[32];

    auto g_x_0_y_0_0_0_z_yy = buffer_1010_sspd[33];

    auto g_x_0_y_0_0_0_z_yz = buffer_1010_sspd[34];

    auto g_x_0_y_0_0_0_z_zz = buffer_1010_sspd[35];

    auto g_x_0_z_0_0_0_x_xx = buffer_1010_sspd[36];

    auto g_x_0_z_0_0_0_x_xy = buffer_1010_sspd[37];

    auto g_x_0_z_0_0_0_x_xz = buffer_1010_sspd[38];

    auto g_x_0_z_0_0_0_x_yy = buffer_1010_sspd[39];

    auto g_x_0_z_0_0_0_x_yz = buffer_1010_sspd[40];

    auto g_x_0_z_0_0_0_x_zz = buffer_1010_sspd[41];

    auto g_x_0_z_0_0_0_y_xx = buffer_1010_sspd[42];

    auto g_x_0_z_0_0_0_y_xy = buffer_1010_sspd[43];

    auto g_x_0_z_0_0_0_y_xz = buffer_1010_sspd[44];

    auto g_x_0_z_0_0_0_y_yy = buffer_1010_sspd[45];

    auto g_x_0_z_0_0_0_y_yz = buffer_1010_sspd[46];

    auto g_x_0_z_0_0_0_y_zz = buffer_1010_sspd[47];

    auto g_x_0_z_0_0_0_z_xx = buffer_1010_sspd[48];

    auto g_x_0_z_0_0_0_z_xy = buffer_1010_sspd[49];

    auto g_x_0_z_0_0_0_z_xz = buffer_1010_sspd[50];

    auto g_x_0_z_0_0_0_z_yy = buffer_1010_sspd[51];

    auto g_x_0_z_0_0_0_z_yz = buffer_1010_sspd[52];

    auto g_x_0_z_0_0_0_z_zz = buffer_1010_sspd[53];

    auto g_y_0_x_0_0_0_x_xx = buffer_1010_sspd[54];

    auto g_y_0_x_0_0_0_x_xy = buffer_1010_sspd[55];

    auto g_y_0_x_0_0_0_x_xz = buffer_1010_sspd[56];

    auto g_y_0_x_0_0_0_x_yy = buffer_1010_sspd[57];

    auto g_y_0_x_0_0_0_x_yz = buffer_1010_sspd[58];

    auto g_y_0_x_0_0_0_x_zz = buffer_1010_sspd[59];

    auto g_y_0_x_0_0_0_y_xx = buffer_1010_sspd[60];

    auto g_y_0_x_0_0_0_y_xy = buffer_1010_sspd[61];

    auto g_y_0_x_0_0_0_y_xz = buffer_1010_sspd[62];

    auto g_y_0_x_0_0_0_y_yy = buffer_1010_sspd[63];

    auto g_y_0_x_0_0_0_y_yz = buffer_1010_sspd[64];

    auto g_y_0_x_0_0_0_y_zz = buffer_1010_sspd[65];

    auto g_y_0_x_0_0_0_z_xx = buffer_1010_sspd[66];

    auto g_y_0_x_0_0_0_z_xy = buffer_1010_sspd[67];

    auto g_y_0_x_0_0_0_z_xz = buffer_1010_sspd[68];

    auto g_y_0_x_0_0_0_z_yy = buffer_1010_sspd[69];

    auto g_y_0_x_0_0_0_z_yz = buffer_1010_sspd[70];

    auto g_y_0_x_0_0_0_z_zz = buffer_1010_sspd[71];

    auto g_y_0_y_0_0_0_x_xx = buffer_1010_sspd[72];

    auto g_y_0_y_0_0_0_x_xy = buffer_1010_sspd[73];

    auto g_y_0_y_0_0_0_x_xz = buffer_1010_sspd[74];

    auto g_y_0_y_0_0_0_x_yy = buffer_1010_sspd[75];

    auto g_y_0_y_0_0_0_x_yz = buffer_1010_sspd[76];

    auto g_y_0_y_0_0_0_x_zz = buffer_1010_sspd[77];

    auto g_y_0_y_0_0_0_y_xx = buffer_1010_sspd[78];

    auto g_y_0_y_0_0_0_y_xy = buffer_1010_sspd[79];

    auto g_y_0_y_0_0_0_y_xz = buffer_1010_sspd[80];

    auto g_y_0_y_0_0_0_y_yy = buffer_1010_sspd[81];

    auto g_y_0_y_0_0_0_y_yz = buffer_1010_sspd[82];

    auto g_y_0_y_0_0_0_y_zz = buffer_1010_sspd[83];

    auto g_y_0_y_0_0_0_z_xx = buffer_1010_sspd[84];

    auto g_y_0_y_0_0_0_z_xy = buffer_1010_sspd[85];

    auto g_y_0_y_0_0_0_z_xz = buffer_1010_sspd[86];

    auto g_y_0_y_0_0_0_z_yy = buffer_1010_sspd[87];

    auto g_y_0_y_0_0_0_z_yz = buffer_1010_sspd[88];

    auto g_y_0_y_0_0_0_z_zz = buffer_1010_sspd[89];

    auto g_y_0_z_0_0_0_x_xx = buffer_1010_sspd[90];

    auto g_y_0_z_0_0_0_x_xy = buffer_1010_sspd[91];

    auto g_y_0_z_0_0_0_x_xz = buffer_1010_sspd[92];

    auto g_y_0_z_0_0_0_x_yy = buffer_1010_sspd[93];

    auto g_y_0_z_0_0_0_x_yz = buffer_1010_sspd[94];

    auto g_y_0_z_0_0_0_x_zz = buffer_1010_sspd[95];

    auto g_y_0_z_0_0_0_y_xx = buffer_1010_sspd[96];

    auto g_y_0_z_0_0_0_y_xy = buffer_1010_sspd[97];

    auto g_y_0_z_0_0_0_y_xz = buffer_1010_sspd[98];

    auto g_y_0_z_0_0_0_y_yy = buffer_1010_sspd[99];

    auto g_y_0_z_0_0_0_y_yz = buffer_1010_sspd[100];

    auto g_y_0_z_0_0_0_y_zz = buffer_1010_sspd[101];

    auto g_y_0_z_0_0_0_z_xx = buffer_1010_sspd[102];

    auto g_y_0_z_0_0_0_z_xy = buffer_1010_sspd[103];

    auto g_y_0_z_0_0_0_z_xz = buffer_1010_sspd[104];

    auto g_y_0_z_0_0_0_z_yy = buffer_1010_sspd[105];

    auto g_y_0_z_0_0_0_z_yz = buffer_1010_sspd[106];

    auto g_y_0_z_0_0_0_z_zz = buffer_1010_sspd[107];

    auto g_z_0_x_0_0_0_x_xx = buffer_1010_sspd[108];

    auto g_z_0_x_0_0_0_x_xy = buffer_1010_sspd[109];

    auto g_z_0_x_0_0_0_x_xz = buffer_1010_sspd[110];

    auto g_z_0_x_0_0_0_x_yy = buffer_1010_sspd[111];

    auto g_z_0_x_0_0_0_x_yz = buffer_1010_sspd[112];

    auto g_z_0_x_0_0_0_x_zz = buffer_1010_sspd[113];

    auto g_z_0_x_0_0_0_y_xx = buffer_1010_sspd[114];

    auto g_z_0_x_0_0_0_y_xy = buffer_1010_sspd[115];

    auto g_z_0_x_0_0_0_y_xz = buffer_1010_sspd[116];

    auto g_z_0_x_0_0_0_y_yy = buffer_1010_sspd[117];

    auto g_z_0_x_0_0_0_y_yz = buffer_1010_sspd[118];

    auto g_z_0_x_0_0_0_y_zz = buffer_1010_sspd[119];

    auto g_z_0_x_0_0_0_z_xx = buffer_1010_sspd[120];

    auto g_z_0_x_0_0_0_z_xy = buffer_1010_sspd[121];

    auto g_z_0_x_0_0_0_z_xz = buffer_1010_sspd[122];

    auto g_z_0_x_0_0_0_z_yy = buffer_1010_sspd[123];

    auto g_z_0_x_0_0_0_z_yz = buffer_1010_sspd[124];

    auto g_z_0_x_0_0_0_z_zz = buffer_1010_sspd[125];

    auto g_z_0_y_0_0_0_x_xx = buffer_1010_sspd[126];

    auto g_z_0_y_0_0_0_x_xy = buffer_1010_sspd[127];

    auto g_z_0_y_0_0_0_x_xz = buffer_1010_sspd[128];

    auto g_z_0_y_0_0_0_x_yy = buffer_1010_sspd[129];

    auto g_z_0_y_0_0_0_x_yz = buffer_1010_sspd[130];

    auto g_z_0_y_0_0_0_x_zz = buffer_1010_sspd[131];

    auto g_z_0_y_0_0_0_y_xx = buffer_1010_sspd[132];

    auto g_z_0_y_0_0_0_y_xy = buffer_1010_sspd[133];

    auto g_z_0_y_0_0_0_y_xz = buffer_1010_sspd[134];

    auto g_z_0_y_0_0_0_y_yy = buffer_1010_sspd[135];

    auto g_z_0_y_0_0_0_y_yz = buffer_1010_sspd[136];

    auto g_z_0_y_0_0_0_y_zz = buffer_1010_sspd[137];

    auto g_z_0_y_0_0_0_z_xx = buffer_1010_sspd[138];

    auto g_z_0_y_0_0_0_z_xy = buffer_1010_sspd[139];

    auto g_z_0_y_0_0_0_z_xz = buffer_1010_sspd[140];

    auto g_z_0_y_0_0_0_z_yy = buffer_1010_sspd[141];

    auto g_z_0_y_0_0_0_z_yz = buffer_1010_sspd[142];

    auto g_z_0_y_0_0_0_z_zz = buffer_1010_sspd[143];

    auto g_z_0_z_0_0_0_x_xx = buffer_1010_sspd[144];

    auto g_z_0_z_0_0_0_x_xy = buffer_1010_sspd[145];

    auto g_z_0_z_0_0_0_x_xz = buffer_1010_sspd[146];

    auto g_z_0_z_0_0_0_x_yy = buffer_1010_sspd[147];

    auto g_z_0_z_0_0_0_x_yz = buffer_1010_sspd[148];

    auto g_z_0_z_0_0_0_x_zz = buffer_1010_sspd[149];

    auto g_z_0_z_0_0_0_y_xx = buffer_1010_sspd[150];

    auto g_z_0_z_0_0_0_y_xy = buffer_1010_sspd[151];

    auto g_z_0_z_0_0_0_y_xz = buffer_1010_sspd[152];

    auto g_z_0_z_0_0_0_y_yy = buffer_1010_sspd[153];

    auto g_z_0_z_0_0_0_y_yz = buffer_1010_sspd[154];

    auto g_z_0_z_0_0_0_y_zz = buffer_1010_sspd[155];

    auto g_z_0_z_0_0_0_z_xx = buffer_1010_sspd[156];

    auto g_z_0_z_0_0_0_z_xy = buffer_1010_sspd[157];

    auto g_z_0_z_0_0_0_z_xz = buffer_1010_sspd[158];

    auto g_z_0_z_0_0_0_z_yy = buffer_1010_sspd[159];

    auto g_z_0_z_0_0_0_z_yz = buffer_1010_sspd[160];

    auto g_z_0_z_0_0_0_z_zz = buffer_1010_sspd[161];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_0_0_xx, g_x_0_0_xy, g_x_0_0_xz, g_x_0_0_yy, g_x_0_0_yz, g_x_0_0_zz, g_x_0_x_0_0_0_x_xx, g_x_0_x_0_0_0_x_xy, g_x_0_x_0_0_0_x_xz, g_x_0_x_0_0_0_x_yy, g_x_0_x_0_0_0_x_yz, g_x_0_x_0_0_0_x_zz, g_x_0_xx_xx, g_x_0_xx_xy, g_x_0_xx_xz, g_x_0_xx_yy, g_x_0_xx_yz, g_x_0_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_0_x_xx[i] = -2.0 * g_x_0_0_xx[i] * a_exp + 4.0 * g_x_0_xx_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_x_xy[i] = -2.0 * g_x_0_0_xy[i] * a_exp + 4.0 * g_x_0_xx_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_x_xz[i] = -2.0 * g_x_0_0_xz[i] * a_exp + 4.0 * g_x_0_xx_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_x_yy[i] = -2.0 * g_x_0_0_yy[i] * a_exp + 4.0 * g_x_0_xx_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_x_yz[i] = -2.0 * g_x_0_0_yz[i] * a_exp + 4.0 * g_x_0_xx_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_x_zz[i] = -2.0 * g_x_0_0_zz[i] * a_exp + 4.0 * g_x_0_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_0_x_0_0_0_y_xx, g_x_0_x_0_0_0_y_xy, g_x_0_x_0_0_0_y_xz, g_x_0_x_0_0_0_y_yy, g_x_0_x_0_0_0_y_yz, g_x_0_x_0_0_0_y_zz, g_x_0_xy_xx, g_x_0_xy_xy, g_x_0_xy_xz, g_x_0_xy_yy, g_x_0_xy_yz, g_x_0_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_0_y_xx[i] = 4.0 * g_x_0_xy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_y_xy[i] = 4.0 * g_x_0_xy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_y_xz[i] = 4.0 * g_x_0_xy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_y_yy[i] = 4.0 * g_x_0_xy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_y_yz[i] = 4.0 * g_x_0_xy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_y_zz[i] = 4.0 * g_x_0_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_0_x_0_0_0_z_xx, g_x_0_x_0_0_0_z_xy, g_x_0_x_0_0_0_z_xz, g_x_0_x_0_0_0_z_yy, g_x_0_x_0_0_0_z_yz, g_x_0_x_0_0_0_z_zz, g_x_0_xz_xx, g_x_0_xz_xy, g_x_0_xz_xz, g_x_0_xz_yy, g_x_0_xz_yz, g_x_0_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_0_z_xx[i] = 4.0 * g_x_0_xz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_z_xy[i] = 4.0 * g_x_0_xz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_z_xz[i] = 4.0 * g_x_0_xz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_z_yy[i] = 4.0 * g_x_0_xz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_z_yz[i] = 4.0 * g_x_0_xz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_z_zz[i] = 4.0 * g_x_0_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_x_0_xy_xx, g_x_0_xy_xy, g_x_0_xy_xz, g_x_0_xy_yy, g_x_0_xy_yz, g_x_0_xy_zz, g_x_0_y_0_0_0_x_xx, g_x_0_y_0_0_0_x_xy, g_x_0_y_0_0_0_x_xz, g_x_0_y_0_0_0_x_yy, g_x_0_y_0_0_0_x_yz, g_x_0_y_0_0_0_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_0_x_xx[i] = 4.0 * g_x_0_xy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_x_xy[i] = 4.0 * g_x_0_xy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_x_xz[i] = 4.0 * g_x_0_xy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_x_yy[i] = 4.0 * g_x_0_xy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_x_yz[i] = 4.0 * g_x_0_xy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_x_zz[i] = 4.0 * g_x_0_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_x_0_0_xx, g_x_0_0_xy, g_x_0_0_xz, g_x_0_0_yy, g_x_0_0_yz, g_x_0_0_zz, g_x_0_y_0_0_0_y_xx, g_x_0_y_0_0_0_y_xy, g_x_0_y_0_0_0_y_xz, g_x_0_y_0_0_0_y_yy, g_x_0_y_0_0_0_y_yz, g_x_0_y_0_0_0_y_zz, g_x_0_yy_xx, g_x_0_yy_xy, g_x_0_yy_xz, g_x_0_yy_yy, g_x_0_yy_yz, g_x_0_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_0_y_xx[i] = -2.0 * g_x_0_0_xx[i] * a_exp + 4.0 * g_x_0_yy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_y_xy[i] = -2.0 * g_x_0_0_xy[i] * a_exp + 4.0 * g_x_0_yy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_y_xz[i] = -2.0 * g_x_0_0_xz[i] * a_exp + 4.0 * g_x_0_yy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_y_yy[i] = -2.0 * g_x_0_0_yy[i] * a_exp + 4.0 * g_x_0_yy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_y_yz[i] = -2.0 * g_x_0_0_yz[i] * a_exp + 4.0 * g_x_0_yy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_y_zz[i] = -2.0 * g_x_0_0_zz[i] * a_exp + 4.0 * g_x_0_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_x_0_y_0_0_0_z_xx, g_x_0_y_0_0_0_z_xy, g_x_0_y_0_0_0_z_xz, g_x_0_y_0_0_0_z_yy, g_x_0_y_0_0_0_z_yz, g_x_0_y_0_0_0_z_zz, g_x_0_yz_xx, g_x_0_yz_xy, g_x_0_yz_xz, g_x_0_yz_yy, g_x_0_yz_yz, g_x_0_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_0_z_xx[i] = 4.0 * g_x_0_yz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_z_xy[i] = 4.0 * g_x_0_yz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_z_xz[i] = 4.0 * g_x_0_yz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_z_yy[i] = 4.0 * g_x_0_yz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_z_yz[i] = 4.0 * g_x_0_yz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_z_zz[i] = 4.0 * g_x_0_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_x_0_xz_xx, g_x_0_xz_xy, g_x_0_xz_xz, g_x_0_xz_yy, g_x_0_xz_yz, g_x_0_xz_zz, g_x_0_z_0_0_0_x_xx, g_x_0_z_0_0_0_x_xy, g_x_0_z_0_0_0_x_xz, g_x_0_z_0_0_0_x_yy, g_x_0_z_0_0_0_x_yz, g_x_0_z_0_0_0_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_0_x_xx[i] = 4.0 * g_x_0_xz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_x_xy[i] = 4.0 * g_x_0_xz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_x_xz[i] = 4.0 * g_x_0_xz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_x_yy[i] = 4.0 * g_x_0_xz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_x_yz[i] = 4.0 * g_x_0_xz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_x_zz[i] = 4.0 * g_x_0_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_x_0_yz_xx, g_x_0_yz_xy, g_x_0_yz_xz, g_x_0_yz_yy, g_x_0_yz_yz, g_x_0_yz_zz, g_x_0_z_0_0_0_y_xx, g_x_0_z_0_0_0_y_xy, g_x_0_z_0_0_0_y_xz, g_x_0_z_0_0_0_y_yy, g_x_0_z_0_0_0_y_yz, g_x_0_z_0_0_0_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_0_y_xx[i] = 4.0 * g_x_0_yz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_y_xy[i] = 4.0 * g_x_0_yz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_y_xz[i] = 4.0 * g_x_0_yz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_y_yy[i] = 4.0 * g_x_0_yz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_y_yz[i] = 4.0 * g_x_0_yz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_y_zz[i] = 4.0 * g_x_0_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_x_0_0_xx, g_x_0_0_xy, g_x_0_0_xz, g_x_0_0_yy, g_x_0_0_yz, g_x_0_0_zz, g_x_0_z_0_0_0_z_xx, g_x_0_z_0_0_0_z_xy, g_x_0_z_0_0_0_z_xz, g_x_0_z_0_0_0_z_yy, g_x_0_z_0_0_0_z_yz, g_x_0_z_0_0_0_z_zz, g_x_0_zz_xx, g_x_0_zz_xy, g_x_0_zz_xz, g_x_0_zz_yy, g_x_0_zz_yz, g_x_0_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_0_z_xx[i] = -2.0 * g_x_0_0_xx[i] * a_exp + 4.0 * g_x_0_zz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_z_xy[i] = -2.0 * g_x_0_0_xy[i] * a_exp + 4.0 * g_x_0_zz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_z_xz[i] = -2.0 * g_x_0_0_xz[i] * a_exp + 4.0 * g_x_0_zz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_z_yy[i] = -2.0 * g_x_0_0_yy[i] * a_exp + 4.0 * g_x_0_zz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_z_yz[i] = -2.0 * g_x_0_0_yz[i] * a_exp + 4.0 * g_x_0_zz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_z_zz[i] = -2.0 * g_x_0_0_zz[i] * a_exp + 4.0 * g_x_0_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_y_0_0_xx, g_y_0_0_xy, g_y_0_0_xz, g_y_0_0_yy, g_y_0_0_yz, g_y_0_0_zz, g_y_0_x_0_0_0_x_xx, g_y_0_x_0_0_0_x_xy, g_y_0_x_0_0_0_x_xz, g_y_0_x_0_0_0_x_yy, g_y_0_x_0_0_0_x_yz, g_y_0_x_0_0_0_x_zz, g_y_0_xx_xx, g_y_0_xx_xy, g_y_0_xx_xz, g_y_0_xx_yy, g_y_0_xx_yz, g_y_0_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_0_x_xx[i] = -2.0 * g_y_0_0_xx[i] * a_exp + 4.0 * g_y_0_xx_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_x_xy[i] = -2.0 * g_y_0_0_xy[i] * a_exp + 4.0 * g_y_0_xx_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_x_xz[i] = -2.0 * g_y_0_0_xz[i] * a_exp + 4.0 * g_y_0_xx_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_x_yy[i] = -2.0 * g_y_0_0_yy[i] * a_exp + 4.0 * g_y_0_xx_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_x_yz[i] = -2.0 * g_y_0_0_yz[i] * a_exp + 4.0 * g_y_0_xx_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_x_zz[i] = -2.0 * g_y_0_0_zz[i] * a_exp + 4.0 * g_y_0_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_y_0_x_0_0_0_y_xx, g_y_0_x_0_0_0_y_xy, g_y_0_x_0_0_0_y_xz, g_y_0_x_0_0_0_y_yy, g_y_0_x_0_0_0_y_yz, g_y_0_x_0_0_0_y_zz, g_y_0_xy_xx, g_y_0_xy_xy, g_y_0_xy_xz, g_y_0_xy_yy, g_y_0_xy_yz, g_y_0_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_0_y_xx[i] = 4.0 * g_y_0_xy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_y_xy[i] = 4.0 * g_y_0_xy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_y_xz[i] = 4.0 * g_y_0_xy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_y_yy[i] = 4.0 * g_y_0_xy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_y_yz[i] = 4.0 * g_y_0_xy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_y_zz[i] = 4.0 * g_y_0_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_y_0_x_0_0_0_z_xx, g_y_0_x_0_0_0_z_xy, g_y_0_x_0_0_0_z_xz, g_y_0_x_0_0_0_z_yy, g_y_0_x_0_0_0_z_yz, g_y_0_x_0_0_0_z_zz, g_y_0_xz_xx, g_y_0_xz_xy, g_y_0_xz_xz, g_y_0_xz_yy, g_y_0_xz_yz, g_y_0_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_0_z_xx[i] = 4.0 * g_y_0_xz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_z_xy[i] = 4.0 * g_y_0_xz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_z_xz[i] = 4.0 * g_y_0_xz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_z_yy[i] = 4.0 * g_y_0_xz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_z_yz[i] = 4.0 * g_y_0_xz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_z_zz[i] = 4.0 * g_y_0_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_y_0_xy_xx, g_y_0_xy_xy, g_y_0_xy_xz, g_y_0_xy_yy, g_y_0_xy_yz, g_y_0_xy_zz, g_y_0_y_0_0_0_x_xx, g_y_0_y_0_0_0_x_xy, g_y_0_y_0_0_0_x_xz, g_y_0_y_0_0_0_x_yy, g_y_0_y_0_0_0_x_yz, g_y_0_y_0_0_0_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_0_x_xx[i] = 4.0 * g_y_0_xy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_x_xy[i] = 4.0 * g_y_0_xy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_x_xz[i] = 4.0 * g_y_0_xy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_x_yy[i] = 4.0 * g_y_0_xy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_x_yz[i] = 4.0 * g_y_0_xy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_x_zz[i] = 4.0 * g_y_0_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_y_0_0_xx, g_y_0_0_xy, g_y_0_0_xz, g_y_0_0_yy, g_y_0_0_yz, g_y_0_0_zz, g_y_0_y_0_0_0_y_xx, g_y_0_y_0_0_0_y_xy, g_y_0_y_0_0_0_y_xz, g_y_0_y_0_0_0_y_yy, g_y_0_y_0_0_0_y_yz, g_y_0_y_0_0_0_y_zz, g_y_0_yy_xx, g_y_0_yy_xy, g_y_0_yy_xz, g_y_0_yy_yy, g_y_0_yy_yz, g_y_0_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_0_y_xx[i] = -2.0 * g_y_0_0_xx[i] * a_exp + 4.0 * g_y_0_yy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_y_xy[i] = -2.0 * g_y_0_0_xy[i] * a_exp + 4.0 * g_y_0_yy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_y_xz[i] = -2.0 * g_y_0_0_xz[i] * a_exp + 4.0 * g_y_0_yy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_y_yy[i] = -2.0 * g_y_0_0_yy[i] * a_exp + 4.0 * g_y_0_yy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_y_yz[i] = -2.0 * g_y_0_0_yz[i] * a_exp + 4.0 * g_y_0_yy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_y_zz[i] = -2.0 * g_y_0_0_zz[i] * a_exp + 4.0 * g_y_0_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_y_0_y_0_0_0_z_xx, g_y_0_y_0_0_0_z_xy, g_y_0_y_0_0_0_z_xz, g_y_0_y_0_0_0_z_yy, g_y_0_y_0_0_0_z_yz, g_y_0_y_0_0_0_z_zz, g_y_0_yz_xx, g_y_0_yz_xy, g_y_0_yz_xz, g_y_0_yz_yy, g_y_0_yz_yz, g_y_0_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_0_z_xx[i] = 4.0 * g_y_0_yz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_z_xy[i] = 4.0 * g_y_0_yz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_z_xz[i] = 4.0 * g_y_0_yz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_z_yy[i] = 4.0 * g_y_0_yz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_z_yz[i] = 4.0 * g_y_0_yz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_z_zz[i] = 4.0 * g_y_0_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_y_0_xz_xx, g_y_0_xz_xy, g_y_0_xz_xz, g_y_0_xz_yy, g_y_0_xz_yz, g_y_0_xz_zz, g_y_0_z_0_0_0_x_xx, g_y_0_z_0_0_0_x_xy, g_y_0_z_0_0_0_x_xz, g_y_0_z_0_0_0_x_yy, g_y_0_z_0_0_0_x_yz, g_y_0_z_0_0_0_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_0_x_xx[i] = 4.0 * g_y_0_xz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_x_xy[i] = 4.0 * g_y_0_xz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_x_xz[i] = 4.0 * g_y_0_xz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_x_yy[i] = 4.0 * g_y_0_xz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_x_yz[i] = 4.0 * g_y_0_xz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_x_zz[i] = 4.0 * g_y_0_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_y_0_yz_xx, g_y_0_yz_xy, g_y_0_yz_xz, g_y_0_yz_yy, g_y_0_yz_yz, g_y_0_yz_zz, g_y_0_z_0_0_0_y_xx, g_y_0_z_0_0_0_y_xy, g_y_0_z_0_0_0_y_xz, g_y_0_z_0_0_0_y_yy, g_y_0_z_0_0_0_y_yz, g_y_0_z_0_0_0_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_0_y_xx[i] = 4.0 * g_y_0_yz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_y_xy[i] = 4.0 * g_y_0_yz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_y_xz[i] = 4.0 * g_y_0_yz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_y_yy[i] = 4.0 * g_y_0_yz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_y_yz[i] = 4.0 * g_y_0_yz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_y_zz[i] = 4.0 * g_y_0_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_y_0_0_xx, g_y_0_0_xy, g_y_0_0_xz, g_y_0_0_yy, g_y_0_0_yz, g_y_0_0_zz, g_y_0_z_0_0_0_z_xx, g_y_0_z_0_0_0_z_xy, g_y_0_z_0_0_0_z_xz, g_y_0_z_0_0_0_z_yy, g_y_0_z_0_0_0_z_yz, g_y_0_z_0_0_0_z_zz, g_y_0_zz_xx, g_y_0_zz_xy, g_y_0_zz_xz, g_y_0_zz_yy, g_y_0_zz_yz, g_y_0_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_0_z_xx[i] = -2.0 * g_y_0_0_xx[i] * a_exp + 4.0 * g_y_0_zz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_z_xy[i] = -2.0 * g_y_0_0_xy[i] * a_exp + 4.0 * g_y_0_zz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_z_xz[i] = -2.0 * g_y_0_0_xz[i] * a_exp + 4.0 * g_y_0_zz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_z_yy[i] = -2.0 * g_y_0_0_yy[i] * a_exp + 4.0 * g_y_0_zz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_z_yz[i] = -2.0 * g_y_0_0_yz[i] * a_exp + 4.0 * g_y_0_zz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_z_zz[i] = -2.0 * g_y_0_0_zz[i] * a_exp + 4.0 * g_y_0_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_z_0_0_xx, g_z_0_0_xy, g_z_0_0_xz, g_z_0_0_yy, g_z_0_0_yz, g_z_0_0_zz, g_z_0_x_0_0_0_x_xx, g_z_0_x_0_0_0_x_xy, g_z_0_x_0_0_0_x_xz, g_z_0_x_0_0_0_x_yy, g_z_0_x_0_0_0_x_yz, g_z_0_x_0_0_0_x_zz, g_z_0_xx_xx, g_z_0_xx_xy, g_z_0_xx_xz, g_z_0_xx_yy, g_z_0_xx_yz, g_z_0_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_0_x_xx[i] = -2.0 * g_z_0_0_xx[i] * a_exp + 4.0 * g_z_0_xx_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_x_xy[i] = -2.0 * g_z_0_0_xy[i] * a_exp + 4.0 * g_z_0_xx_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_x_xz[i] = -2.0 * g_z_0_0_xz[i] * a_exp + 4.0 * g_z_0_xx_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_x_yy[i] = -2.0 * g_z_0_0_yy[i] * a_exp + 4.0 * g_z_0_xx_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_x_yz[i] = -2.0 * g_z_0_0_yz[i] * a_exp + 4.0 * g_z_0_xx_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_x_zz[i] = -2.0 * g_z_0_0_zz[i] * a_exp + 4.0 * g_z_0_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_z_0_x_0_0_0_y_xx, g_z_0_x_0_0_0_y_xy, g_z_0_x_0_0_0_y_xz, g_z_0_x_0_0_0_y_yy, g_z_0_x_0_0_0_y_yz, g_z_0_x_0_0_0_y_zz, g_z_0_xy_xx, g_z_0_xy_xy, g_z_0_xy_xz, g_z_0_xy_yy, g_z_0_xy_yz, g_z_0_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_0_y_xx[i] = 4.0 * g_z_0_xy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_y_xy[i] = 4.0 * g_z_0_xy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_y_xz[i] = 4.0 * g_z_0_xy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_y_yy[i] = 4.0 * g_z_0_xy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_y_yz[i] = 4.0 * g_z_0_xy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_y_zz[i] = 4.0 * g_z_0_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_z_0_x_0_0_0_z_xx, g_z_0_x_0_0_0_z_xy, g_z_0_x_0_0_0_z_xz, g_z_0_x_0_0_0_z_yy, g_z_0_x_0_0_0_z_yz, g_z_0_x_0_0_0_z_zz, g_z_0_xz_xx, g_z_0_xz_xy, g_z_0_xz_xz, g_z_0_xz_yy, g_z_0_xz_yz, g_z_0_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_0_z_xx[i] = 4.0 * g_z_0_xz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_z_xy[i] = 4.0 * g_z_0_xz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_z_xz[i] = 4.0 * g_z_0_xz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_z_yy[i] = 4.0 * g_z_0_xz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_z_yz[i] = 4.0 * g_z_0_xz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_z_zz[i] = 4.0 * g_z_0_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_z_0_xy_xx, g_z_0_xy_xy, g_z_0_xy_xz, g_z_0_xy_yy, g_z_0_xy_yz, g_z_0_xy_zz, g_z_0_y_0_0_0_x_xx, g_z_0_y_0_0_0_x_xy, g_z_0_y_0_0_0_x_xz, g_z_0_y_0_0_0_x_yy, g_z_0_y_0_0_0_x_yz, g_z_0_y_0_0_0_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_0_x_xx[i] = 4.0 * g_z_0_xy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_x_xy[i] = 4.0 * g_z_0_xy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_x_xz[i] = 4.0 * g_z_0_xy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_x_yy[i] = 4.0 * g_z_0_xy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_x_yz[i] = 4.0 * g_z_0_xy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_x_zz[i] = 4.0 * g_z_0_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_z_0_0_xx, g_z_0_0_xy, g_z_0_0_xz, g_z_0_0_yy, g_z_0_0_yz, g_z_0_0_zz, g_z_0_y_0_0_0_y_xx, g_z_0_y_0_0_0_y_xy, g_z_0_y_0_0_0_y_xz, g_z_0_y_0_0_0_y_yy, g_z_0_y_0_0_0_y_yz, g_z_0_y_0_0_0_y_zz, g_z_0_yy_xx, g_z_0_yy_xy, g_z_0_yy_xz, g_z_0_yy_yy, g_z_0_yy_yz, g_z_0_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_0_y_xx[i] = -2.0 * g_z_0_0_xx[i] * a_exp + 4.0 * g_z_0_yy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_y_xy[i] = -2.0 * g_z_0_0_xy[i] * a_exp + 4.0 * g_z_0_yy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_y_xz[i] = -2.0 * g_z_0_0_xz[i] * a_exp + 4.0 * g_z_0_yy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_y_yy[i] = -2.0 * g_z_0_0_yy[i] * a_exp + 4.0 * g_z_0_yy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_y_yz[i] = -2.0 * g_z_0_0_yz[i] * a_exp + 4.0 * g_z_0_yy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_y_zz[i] = -2.0 * g_z_0_0_zz[i] * a_exp + 4.0 * g_z_0_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_z_0_y_0_0_0_z_xx, g_z_0_y_0_0_0_z_xy, g_z_0_y_0_0_0_z_xz, g_z_0_y_0_0_0_z_yy, g_z_0_y_0_0_0_z_yz, g_z_0_y_0_0_0_z_zz, g_z_0_yz_xx, g_z_0_yz_xy, g_z_0_yz_xz, g_z_0_yz_yy, g_z_0_yz_yz, g_z_0_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_0_z_xx[i] = 4.0 * g_z_0_yz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_z_xy[i] = 4.0 * g_z_0_yz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_z_xz[i] = 4.0 * g_z_0_yz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_z_yy[i] = 4.0 * g_z_0_yz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_z_yz[i] = 4.0 * g_z_0_yz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_z_zz[i] = 4.0 * g_z_0_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_z_0_xz_xx, g_z_0_xz_xy, g_z_0_xz_xz, g_z_0_xz_yy, g_z_0_xz_yz, g_z_0_xz_zz, g_z_0_z_0_0_0_x_xx, g_z_0_z_0_0_0_x_xy, g_z_0_z_0_0_0_x_xz, g_z_0_z_0_0_0_x_yy, g_z_0_z_0_0_0_x_yz, g_z_0_z_0_0_0_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_0_x_xx[i] = 4.0 * g_z_0_xz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_x_xy[i] = 4.0 * g_z_0_xz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_x_xz[i] = 4.0 * g_z_0_xz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_x_yy[i] = 4.0 * g_z_0_xz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_x_yz[i] = 4.0 * g_z_0_xz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_x_zz[i] = 4.0 * g_z_0_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_z_0_yz_xx, g_z_0_yz_xy, g_z_0_yz_xz, g_z_0_yz_yy, g_z_0_yz_yz, g_z_0_yz_zz, g_z_0_z_0_0_0_y_xx, g_z_0_z_0_0_0_y_xy, g_z_0_z_0_0_0_y_xz, g_z_0_z_0_0_0_y_yy, g_z_0_z_0_0_0_y_yz, g_z_0_z_0_0_0_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_0_y_xx[i] = 4.0 * g_z_0_yz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_y_xy[i] = 4.0 * g_z_0_yz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_y_xz[i] = 4.0 * g_z_0_yz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_y_yy[i] = 4.0 * g_z_0_yz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_y_yz[i] = 4.0 * g_z_0_yz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_y_zz[i] = 4.0 * g_z_0_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_z_0_0_xx, g_z_0_0_xy, g_z_0_0_xz, g_z_0_0_yy, g_z_0_0_yz, g_z_0_0_zz, g_z_0_z_0_0_0_z_xx, g_z_0_z_0_0_0_z_xy, g_z_0_z_0_0_0_z_xz, g_z_0_z_0_0_0_z_yy, g_z_0_z_0_0_0_z_yz, g_z_0_z_0_0_0_z_zz, g_z_0_zz_xx, g_z_0_zz_xy, g_z_0_zz_xz, g_z_0_zz_yy, g_z_0_zz_yz, g_z_0_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_0_z_xx[i] = -2.0 * g_z_0_0_xx[i] * a_exp + 4.0 * g_z_0_zz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_z_xy[i] = -2.0 * g_z_0_0_xy[i] * a_exp + 4.0 * g_z_0_zz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_z_xz[i] = -2.0 * g_z_0_0_xz[i] * a_exp + 4.0 * g_z_0_zz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_z_yy[i] = -2.0 * g_z_0_0_yy[i] * a_exp + 4.0 * g_z_0_zz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_z_yz[i] = -2.0 * g_z_0_0_yz[i] * a_exp + 4.0 * g_z_0_zz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_z_zz[i] = -2.0 * g_z_0_0_zz[i] * a_exp + 4.0 * g_z_0_zz_zz[i] * a_exp * c_exps[i];
    }
}

} // t4c_geom namespace


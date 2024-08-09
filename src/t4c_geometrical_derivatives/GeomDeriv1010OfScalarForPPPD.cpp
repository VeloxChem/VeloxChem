#include "GeomDeriv1010OfScalarForPPPD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1010_pppd_0(CSimdArray<double>& buffer_1010_pppd,
                     const CSimdArray<double>& buffer_spsd,
                     const CSimdArray<double>& buffer_spdd,
                     const CSimdArray<double>& buffer_dpsd,
                     const CSimdArray<double>& buffer_dpdd,
                     const double a_exp,
                     const double* c_exps) -> void
{
    const auto ndims = buffer_1010_pppd.number_of_columns();

    /// Set up components of auxilary buffer : buffer_spsd

    auto g_0_x_0_xx = buffer_spsd[0];

    auto g_0_x_0_xy = buffer_spsd[1];

    auto g_0_x_0_xz = buffer_spsd[2];

    auto g_0_x_0_yy = buffer_spsd[3];

    auto g_0_x_0_yz = buffer_spsd[4];

    auto g_0_x_0_zz = buffer_spsd[5];

    auto g_0_y_0_xx = buffer_spsd[6];

    auto g_0_y_0_xy = buffer_spsd[7];

    auto g_0_y_0_xz = buffer_spsd[8];

    auto g_0_y_0_yy = buffer_spsd[9];

    auto g_0_y_0_yz = buffer_spsd[10];

    auto g_0_y_0_zz = buffer_spsd[11];

    auto g_0_z_0_xx = buffer_spsd[12];

    auto g_0_z_0_xy = buffer_spsd[13];

    auto g_0_z_0_xz = buffer_spsd[14];

    auto g_0_z_0_yy = buffer_spsd[15];

    auto g_0_z_0_yz = buffer_spsd[16];

    auto g_0_z_0_zz = buffer_spsd[17];

    /// Set up components of auxilary buffer : buffer_spdd

    auto g_0_x_xx_xx = buffer_spdd[0];

    auto g_0_x_xx_xy = buffer_spdd[1];

    auto g_0_x_xx_xz = buffer_spdd[2];

    auto g_0_x_xx_yy = buffer_spdd[3];

    auto g_0_x_xx_yz = buffer_spdd[4];

    auto g_0_x_xx_zz = buffer_spdd[5];

    auto g_0_x_xy_xx = buffer_spdd[6];

    auto g_0_x_xy_xy = buffer_spdd[7];

    auto g_0_x_xy_xz = buffer_spdd[8];

    auto g_0_x_xy_yy = buffer_spdd[9];

    auto g_0_x_xy_yz = buffer_spdd[10];

    auto g_0_x_xy_zz = buffer_spdd[11];

    auto g_0_x_xz_xx = buffer_spdd[12];

    auto g_0_x_xz_xy = buffer_spdd[13];

    auto g_0_x_xz_xz = buffer_spdd[14];

    auto g_0_x_xz_yy = buffer_spdd[15];

    auto g_0_x_xz_yz = buffer_spdd[16];

    auto g_0_x_xz_zz = buffer_spdd[17];

    auto g_0_x_yy_xx = buffer_spdd[18];

    auto g_0_x_yy_xy = buffer_spdd[19];

    auto g_0_x_yy_xz = buffer_spdd[20];

    auto g_0_x_yy_yy = buffer_spdd[21];

    auto g_0_x_yy_yz = buffer_spdd[22];

    auto g_0_x_yy_zz = buffer_spdd[23];

    auto g_0_x_yz_xx = buffer_spdd[24];

    auto g_0_x_yz_xy = buffer_spdd[25];

    auto g_0_x_yz_xz = buffer_spdd[26];

    auto g_0_x_yz_yy = buffer_spdd[27];

    auto g_0_x_yz_yz = buffer_spdd[28];

    auto g_0_x_yz_zz = buffer_spdd[29];

    auto g_0_x_zz_xx = buffer_spdd[30];

    auto g_0_x_zz_xy = buffer_spdd[31];

    auto g_0_x_zz_xz = buffer_spdd[32];

    auto g_0_x_zz_yy = buffer_spdd[33];

    auto g_0_x_zz_yz = buffer_spdd[34];

    auto g_0_x_zz_zz = buffer_spdd[35];

    auto g_0_y_xx_xx = buffer_spdd[36];

    auto g_0_y_xx_xy = buffer_spdd[37];

    auto g_0_y_xx_xz = buffer_spdd[38];

    auto g_0_y_xx_yy = buffer_spdd[39];

    auto g_0_y_xx_yz = buffer_spdd[40];

    auto g_0_y_xx_zz = buffer_spdd[41];

    auto g_0_y_xy_xx = buffer_spdd[42];

    auto g_0_y_xy_xy = buffer_spdd[43];

    auto g_0_y_xy_xz = buffer_spdd[44];

    auto g_0_y_xy_yy = buffer_spdd[45];

    auto g_0_y_xy_yz = buffer_spdd[46];

    auto g_0_y_xy_zz = buffer_spdd[47];

    auto g_0_y_xz_xx = buffer_spdd[48];

    auto g_0_y_xz_xy = buffer_spdd[49];

    auto g_0_y_xz_xz = buffer_spdd[50];

    auto g_0_y_xz_yy = buffer_spdd[51];

    auto g_0_y_xz_yz = buffer_spdd[52];

    auto g_0_y_xz_zz = buffer_spdd[53];

    auto g_0_y_yy_xx = buffer_spdd[54];

    auto g_0_y_yy_xy = buffer_spdd[55];

    auto g_0_y_yy_xz = buffer_spdd[56];

    auto g_0_y_yy_yy = buffer_spdd[57];

    auto g_0_y_yy_yz = buffer_spdd[58];

    auto g_0_y_yy_zz = buffer_spdd[59];

    auto g_0_y_yz_xx = buffer_spdd[60];

    auto g_0_y_yz_xy = buffer_spdd[61];

    auto g_0_y_yz_xz = buffer_spdd[62];

    auto g_0_y_yz_yy = buffer_spdd[63];

    auto g_0_y_yz_yz = buffer_spdd[64];

    auto g_0_y_yz_zz = buffer_spdd[65];

    auto g_0_y_zz_xx = buffer_spdd[66];

    auto g_0_y_zz_xy = buffer_spdd[67];

    auto g_0_y_zz_xz = buffer_spdd[68];

    auto g_0_y_zz_yy = buffer_spdd[69];

    auto g_0_y_zz_yz = buffer_spdd[70];

    auto g_0_y_zz_zz = buffer_spdd[71];

    auto g_0_z_xx_xx = buffer_spdd[72];

    auto g_0_z_xx_xy = buffer_spdd[73];

    auto g_0_z_xx_xz = buffer_spdd[74];

    auto g_0_z_xx_yy = buffer_spdd[75];

    auto g_0_z_xx_yz = buffer_spdd[76];

    auto g_0_z_xx_zz = buffer_spdd[77];

    auto g_0_z_xy_xx = buffer_spdd[78];

    auto g_0_z_xy_xy = buffer_spdd[79];

    auto g_0_z_xy_xz = buffer_spdd[80];

    auto g_0_z_xy_yy = buffer_spdd[81];

    auto g_0_z_xy_yz = buffer_spdd[82];

    auto g_0_z_xy_zz = buffer_spdd[83];

    auto g_0_z_xz_xx = buffer_spdd[84];

    auto g_0_z_xz_xy = buffer_spdd[85];

    auto g_0_z_xz_xz = buffer_spdd[86];

    auto g_0_z_xz_yy = buffer_spdd[87];

    auto g_0_z_xz_yz = buffer_spdd[88];

    auto g_0_z_xz_zz = buffer_spdd[89];

    auto g_0_z_yy_xx = buffer_spdd[90];

    auto g_0_z_yy_xy = buffer_spdd[91];

    auto g_0_z_yy_xz = buffer_spdd[92];

    auto g_0_z_yy_yy = buffer_spdd[93];

    auto g_0_z_yy_yz = buffer_spdd[94];

    auto g_0_z_yy_zz = buffer_spdd[95];

    auto g_0_z_yz_xx = buffer_spdd[96];

    auto g_0_z_yz_xy = buffer_spdd[97];

    auto g_0_z_yz_xz = buffer_spdd[98];

    auto g_0_z_yz_yy = buffer_spdd[99];

    auto g_0_z_yz_yz = buffer_spdd[100];

    auto g_0_z_yz_zz = buffer_spdd[101];

    auto g_0_z_zz_xx = buffer_spdd[102];

    auto g_0_z_zz_xy = buffer_spdd[103];

    auto g_0_z_zz_xz = buffer_spdd[104];

    auto g_0_z_zz_yy = buffer_spdd[105];

    auto g_0_z_zz_yz = buffer_spdd[106];

    auto g_0_z_zz_zz = buffer_spdd[107];

    /// Set up components of auxilary buffer : buffer_dpsd

    auto g_xx_x_0_xx = buffer_dpsd[0];

    auto g_xx_x_0_xy = buffer_dpsd[1];

    auto g_xx_x_0_xz = buffer_dpsd[2];

    auto g_xx_x_0_yy = buffer_dpsd[3];

    auto g_xx_x_0_yz = buffer_dpsd[4];

    auto g_xx_x_0_zz = buffer_dpsd[5];

    auto g_xx_y_0_xx = buffer_dpsd[6];

    auto g_xx_y_0_xy = buffer_dpsd[7];

    auto g_xx_y_0_xz = buffer_dpsd[8];

    auto g_xx_y_0_yy = buffer_dpsd[9];

    auto g_xx_y_0_yz = buffer_dpsd[10];

    auto g_xx_y_0_zz = buffer_dpsd[11];

    auto g_xx_z_0_xx = buffer_dpsd[12];

    auto g_xx_z_0_xy = buffer_dpsd[13];

    auto g_xx_z_0_xz = buffer_dpsd[14];

    auto g_xx_z_0_yy = buffer_dpsd[15];

    auto g_xx_z_0_yz = buffer_dpsd[16];

    auto g_xx_z_0_zz = buffer_dpsd[17];

    auto g_xy_x_0_xx = buffer_dpsd[18];

    auto g_xy_x_0_xy = buffer_dpsd[19];

    auto g_xy_x_0_xz = buffer_dpsd[20];

    auto g_xy_x_0_yy = buffer_dpsd[21];

    auto g_xy_x_0_yz = buffer_dpsd[22];

    auto g_xy_x_0_zz = buffer_dpsd[23];

    auto g_xy_y_0_xx = buffer_dpsd[24];

    auto g_xy_y_0_xy = buffer_dpsd[25];

    auto g_xy_y_0_xz = buffer_dpsd[26];

    auto g_xy_y_0_yy = buffer_dpsd[27];

    auto g_xy_y_0_yz = buffer_dpsd[28];

    auto g_xy_y_0_zz = buffer_dpsd[29];

    auto g_xy_z_0_xx = buffer_dpsd[30];

    auto g_xy_z_0_xy = buffer_dpsd[31];

    auto g_xy_z_0_xz = buffer_dpsd[32];

    auto g_xy_z_0_yy = buffer_dpsd[33];

    auto g_xy_z_0_yz = buffer_dpsd[34];

    auto g_xy_z_0_zz = buffer_dpsd[35];

    auto g_xz_x_0_xx = buffer_dpsd[36];

    auto g_xz_x_0_xy = buffer_dpsd[37];

    auto g_xz_x_0_xz = buffer_dpsd[38];

    auto g_xz_x_0_yy = buffer_dpsd[39];

    auto g_xz_x_0_yz = buffer_dpsd[40];

    auto g_xz_x_0_zz = buffer_dpsd[41];

    auto g_xz_y_0_xx = buffer_dpsd[42];

    auto g_xz_y_0_xy = buffer_dpsd[43];

    auto g_xz_y_0_xz = buffer_dpsd[44];

    auto g_xz_y_0_yy = buffer_dpsd[45];

    auto g_xz_y_0_yz = buffer_dpsd[46];

    auto g_xz_y_0_zz = buffer_dpsd[47];

    auto g_xz_z_0_xx = buffer_dpsd[48];

    auto g_xz_z_0_xy = buffer_dpsd[49];

    auto g_xz_z_0_xz = buffer_dpsd[50];

    auto g_xz_z_0_yy = buffer_dpsd[51];

    auto g_xz_z_0_yz = buffer_dpsd[52];

    auto g_xz_z_0_zz = buffer_dpsd[53];

    auto g_yy_x_0_xx = buffer_dpsd[54];

    auto g_yy_x_0_xy = buffer_dpsd[55];

    auto g_yy_x_0_xz = buffer_dpsd[56];

    auto g_yy_x_0_yy = buffer_dpsd[57];

    auto g_yy_x_0_yz = buffer_dpsd[58];

    auto g_yy_x_0_zz = buffer_dpsd[59];

    auto g_yy_y_0_xx = buffer_dpsd[60];

    auto g_yy_y_0_xy = buffer_dpsd[61];

    auto g_yy_y_0_xz = buffer_dpsd[62];

    auto g_yy_y_0_yy = buffer_dpsd[63];

    auto g_yy_y_0_yz = buffer_dpsd[64];

    auto g_yy_y_0_zz = buffer_dpsd[65];

    auto g_yy_z_0_xx = buffer_dpsd[66];

    auto g_yy_z_0_xy = buffer_dpsd[67];

    auto g_yy_z_0_xz = buffer_dpsd[68];

    auto g_yy_z_0_yy = buffer_dpsd[69];

    auto g_yy_z_0_yz = buffer_dpsd[70];

    auto g_yy_z_0_zz = buffer_dpsd[71];

    auto g_yz_x_0_xx = buffer_dpsd[72];

    auto g_yz_x_0_xy = buffer_dpsd[73];

    auto g_yz_x_0_xz = buffer_dpsd[74];

    auto g_yz_x_0_yy = buffer_dpsd[75];

    auto g_yz_x_0_yz = buffer_dpsd[76];

    auto g_yz_x_0_zz = buffer_dpsd[77];

    auto g_yz_y_0_xx = buffer_dpsd[78];

    auto g_yz_y_0_xy = buffer_dpsd[79];

    auto g_yz_y_0_xz = buffer_dpsd[80];

    auto g_yz_y_0_yy = buffer_dpsd[81];

    auto g_yz_y_0_yz = buffer_dpsd[82];

    auto g_yz_y_0_zz = buffer_dpsd[83];

    auto g_yz_z_0_xx = buffer_dpsd[84];

    auto g_yz_z_0_xy = buffer_dpsd[85];

    auto g_yz_z_0_xz = buffer_dpsd[86];

    auto g_yz_z_0_yy = buffer_dpsd[87];

    auto g_yz_z_0_yz = buffer_dpsd[88];

    auto g_yz_z_0_zz = buffer_dpsd[89];

    auto g_zz_x_0_xx = buffer_dpsd[90];

    auto g_zz_x_0_xy = buffer_dpsd[91];

    auto g_zz_x_0_xz = buffer_dpsd[92];

    auto g_zz_x_0_yy = buffer_dpsd[93];

    auto g_zz_x_0_yz = buffer_dpsd[94];

    auto g_zz_x_0_zz = buffer_dpsd[95];

    auto g_zz_y_0_xx = buffer_dpsd[96];

    auto g_zz_y_0_xy = buffer_dpsd[97];

    auto g_zz_y_0_xz = buffer_dpsd[98];

    auto g_zz_y_0_yy = buffer_dpsd[99];

    auto g_zz_y_0_yz = buffer_dpsd[100];

    auto g_zz_y_0_zz = buffer_dpsd[101];

    auto g_zz_z_0_xx = buffer_dpsd[102];

    auto g_zz_z_0_xy = buffer_dpsd[103];

    auto g_zz_z_0_xz = buffer_dpsd[104];

    auto g_zz_z_0_yy = buffer_dpsd[105];

    auto g_zz_z_0_yz = buffer_dpsd[106];

    auto g_zz_z_0_zz = buffer_dpsd[107];

    /// Set up components of auxilary buffer : buffer_dpdd

    auto g_xx_x_xx_xx = buffer_dpdd[0];

    auto g_xx_x_xx_xy = buffer_dpdd[1];

    auto g_xx_x_xx_xz = buffer_dpdd[2];

    auto g_xx_x_xx_yy = buffer_dpdd[3];

    auto g_xx_x_xx_yz = buffer_dpdd[4];

    auto g_xx_x_xx_zz = buffer_dpdd[5];

    auto g_xx_x_xy_xx = buffer_dpdd[6];

    auto g_xx_x_xy_xy = buffer_dpdd[7];

    auto g_xx_x_xy_xz = buffer_dpdd[8];

    auto g_xx_x_xy_yy = buffer_dpdd[9];

    auto g_xx_x_xy_yz = buffer_dpdd[10];

    auto g_xx_x_xy_zz = buffer_dpdd[11];

    auto g_xx_x_xz_xx = buffer_dpdd[12];

    auto g_xx_x_xz_xy = buffer_dpdd[13];

    auto g_xx_x_xz_xz = buffer_dpdd[14];

    auto g_xx_x_xz_yy = buffer_dpdd[15];

    auto g_xx_x_xz_yz = buffer_dpdd[16];

    auto g_xx_x_xz_zz = buffer_dpdd[17];

    auto g_xx_x_yy_xx = buffer_dpdd[18];

    auto g_xx_x_yy_xy = buffer_dpdd[19];

    auto g_xx_x_yy_xz = buffer_dpdd[20];

    auto g_xx_x_yy_yy = buffer_dpdd[21];

    auto g_xx_x_yy_yz = buffer_dpdd[22];

    auto g_xx_x_yy_zz = buffer_dpdd[23];

    auto g_xx_x_yz_xx = buffer_dpdd[24];

    auto g_xx_x_yz_xy = buffer_dpdd[25];

    auto g_xx_x_yz_xz = buffer_dpdd[26];

    auto g_xx_x_yz_yy = buffer_dpdd[27];

    auto g_xx_x_yz_yz = buffer_dpdd[28];

    auto g_xx_x_yz_zz = buffer_dpdd[29];

    auto g_xx_x_zz_xx = buffer_dpdd[30];

    auto g_xx_x_zz_xy = buffer_dpdd[31];

    auto g_xx_x_zz_xz = buffer_dpdd[32];

    auto g_xx_x_zz_yy = buffer_dpdd[33];

    auto g_xx_x_zz_yz = buffer_dpdd[34];

    auto g_xx_x_zz_zz = buffer_dpdd[35];

    auto g_xx_y_xx_xx = buffer_dpdd[36];

    auto g_xx_y_xx_xy = buffer_dpdd[37];

    auto g_xx_y_xx_xz = buffer_dpdd[38];

    auto g_xx_y_xx_yy = buffer_dpdd[39];

    auto g_xx_y_xx_yz = buffer_dpdd[40];

    auto g_xx_y_xx_zz = buffer_dpdd[41];

    auto g_xx_y_xy_xx = buffer_dpdd[42];

    auto g_xx_y_xy_xy = buffer_dpdd[43];

    auto g_xx_y_xy_xz = buffer_dpdd[44];

    auto g_xx_y_xy_yy = buffer_dpdd[45];

    auto g_xx_y_xy_yz = buffer_dpdd[46];

    auto g_xx_y_xy_zz = buffer_dpdd[47];

    auto g_xx_y_xz_xx = buffer_dpdd[48];

    auto g_xx_y_xz_xy = buffer_dpdd[49];

    auto g_xx_y_xz_xz = buffer_dpdd[50];

    auto g_xx_y_xz_yy = buffer_dpdd[51];

    auto g_xx_y_xz_yz = buffer_dpdd[52];

    auto g_xx_y_xz_zz = buffer_dpdd[53];

    auto g_xx_y_yy_xx = buffer_dpdd[54];

    auto g_xx_y_yy_xy = buffer_dpdd[55];

    auto g_xx_y_yy_xz = buffer_dpdd[56];

    auto g_xx_y_yy_yy = buffer_dpdd[57];

    auto g_xx_y_yy_yz = buffer_dpdd[58];

    auto g_xx_y_yy_zz = buffer_dpdd[59];

    auto g_xx_y_yz_xx = buffer_dpdd[60];

    auto g_xx_y_yz_xy = buffer_dpdd[61];

    auto g_xx_y_yz_xz = buffer_dpdd[62];

    auto g_xx_y_yz_yy = buffer_dpdd[63];

    auto g_xx_y_yz_yz = buffer_dpdd[64];

    auto g_xx_y_yz_zz = buffer_dpdd[65];

    auto g_xx_y_zz_xx = buffer_dpdd[66];

    auto g_xx_y_zz_xy = buffer_dpdd[67];

    auto g_xx_y_zz_xz = buffer_dpdd[68];

    auto g_xx_y_zz_yy = buffer_dpdd[69];

    auto g_xx_y_zz_yz = buffer_dpdd[70];

    auto g_xx_y_zz_zz = buffer_dpdd[71];

    auto g_xx_z_xx_xx = buffer_dpdd[72];

    auto g_xx_z_xx_xy = buffer_dpdd[73];

    auto g_xx_z_xx_xz = buffer_dpdd[74];

    auto g_xx_z_xx_yy = buffer_dpdd[75];

    auto g_xx_z_xx_yz = buffer_dpdd[76];

    auto g_xx_z_xx_zz = buffer_dpdd[77];

    auto g_xx_z_xy_xx = buffer_dpdd[78];

    auto g_xx_z_xy_xy = buffer_dpdd[79];

    auto g_xx_z_xy_xz = buffer_dpdd[80];

    auto g_xx_z_xy_yy = buffer_dpdd[81];

    auto g_xx_z_xy_yz = buffer_dpdd[82];

    auto g_xx_z_xy_zz = buffer_dpdd[83];

    auto g_xx_z_xz_xx = buffer_dpdd[84];

    auto g_xx_z_xz_xy = buffer_dpdd[85];

    auto g_xx_z_xz_xz = buffer_dpdd[86];

    auto g_xx_z_xz_yy = buffer_dpdd[87];

    auto g_xx_z_xz_yz = buffer_dpdd[88];

    auto g_xx_z_xz_zz = buffer_dpdd[89];

    auto g_xx_z_yy_xx = buffer_dpdd[90];

    auto g_xx_z_yy_xy = buffer_dpdd[91];

    auto g_xx_z_yy_xz = buffer_dpdd[92];

    auto g_xx_z_yy_yy = buffer_dpdd[93];

    auto g_xx_z_yy_yz = buffer_dpdd[94];

    auto g_xx_z_yy_zz = buffer_dpdd[95];

    auto g_xx_z_yz_xx = buffer_dpdd[96];

    auto g_xx_z_yz_xy = buffer_dpdd[97];

    auto g_xx_z_yz_xz = buffer_dpdd[98];

    auto g_xx_z_yz_yy = buffer_dpdd[99];

    auto g_xx_z_yz_yz = buffer_dpdd[100];

    auto g_xx_z_yz_zz = buffer_dpdd[101];

    auto g_xx_z_zz_xx = buffer_dpdd[102];

    auto g_xx_z_zz_xy = buffer_dpdd[103];

    auto g_xx_z_zz_xz = buffer_dpdd[104];

    auto g_xx_z_zz_yy = buffer_dpdd[105];

    auto g_xx_z_zz_yz = buffer_dpdd[106];

    auto g_xx_z_zz_zz = buffer_dpdd[107];

    auto g_xy_x_xx_xx = buffer_dpdd[108];

    auto g_xy_x_xx_xy = buffer_dpdd[109];

    auto g_xy_x_xx_xz = buffer_dpdd[110];

    auto g_xy_x_xx_yy = buffer_dpdd[111];

    auto g_xy_x_xx_yz = buffer_dpdd[112];

    auto g_xy_x_xx_zz = buffer_dpdd[113];

    auto g_xy_x_xy_xx = buffer_dpdd[114];

    auto g_xy_x_xy_xy = buffer_dpdd[115];

    auto g_xy_x_xy_xz = buffer_dpdd[116];

    auto g_xy_x_xy_yy = buffer_dpdd[117];

    auto g_xy_x_xy_yz = buffer_dpdd[118];

    auto g_xy_x_xy_zz = buffer_dpdd[119];

    auto g_xy_x_xz_xx = buffer_dpdd[120];

    auto g_xy_x_xz_xy = buffer_dpdd[121];

    auto g_xy_x_xz_xz = buffer_dpdd[122];

    auto g_xy_x_xz_yy = buffer_dpdd[123];

    auto g_xy_x_xz_yz = buffer_dpdd[124];

    auto g_xy_x_xz_zz = buffer_dpdd[125];

    auto g_xy_x_yy_xx = buffer_dpdd[126];

    auto g_xy_x_yy_xy = buffer_dpdd[127];

    auto g_xy_x_yy_xz = buffer_dpdd[128];

    auto g_xy_x_yy_yy = buffer_dpdd[129];

    auto g_xy_x_yy_yz = buffer_dpdd[130];

    auto g_xy_x_yy_zz = buffer_dpdd[131];

    auto g_xy_x_yz_xx = buffer_dpdd[132];

    auto g_xy_x_yz_xy = buffer_dpdd[133];

    auto g_xy_x_yz_xz = buffer_dpdd[134];

    auto g_xy_x_yz_yy = buffer_dpdd[135];

    auto g_xy_x_yz_yz = buffer_dpdd[136];

    auto g_xy_x_yz_zz = buffer_dpdd[137];

    auto g_xy_x_zz_xx = buffer_dpdd[138];

    auto g_xy_x_zz_xy = buffer_dpdd[139];

    auto g_xy_x_zz_xz = buffer_dpdd[140];

    auto g_xy_x_zz_yy = buffer_dpdd[141];

    auto g_xy_x_zz_yz = buffer_dpdd[142];

    auto g_xy_x_zz_zz = buffer_dpdd[143];

    auto g_xy_y_xx_xx = buffer_dpdd[144];

    auto g_xy_y_xx_xy = buffer_dpdd[145];

    auto g_xy_y_xx_xz = buffer_dpdd[146];

    auto g_xy_y_xx_yy = buffer_dpdd[147];

    auto g_xy_y_xx_yz = buffer_dpdd[148];

    auto g_xy_y_xx_zz = buffer_dpdd[149];

    auto g_xy_y_xy_xx = buffer_dpdd[150];

    auto g_xy_y_xy_xy = buffer_dpdd[151];

    auto g_xy_y_xy_xz = buffer_dpdd[152];

    auto g_xy_y_xy_yy = buffer_dpdd[153];

    auto g_xy_y_xy_yz = buffer_dpdd[154];

    auto g_xy_y_xy_zz = buffer_dpdd[155];

    auto g_xy_y_xz_xx = buffer_dpdd[156];

    auto g_xy_y_xz_xy = buffer_dpdd[157];

    auto g_xy_y_xz_xz = buffer_dpdd[158];

    auto g_xy_y_xz_yy = buffer_dpdd[159];

    auto g_xy_y_xz_yz = buffer_dpdd[160];

    auto g_xy_y_xz_zz = buffer_dpdd[161];

    auto g_xy_y_yy_xx = buffer_dpdd[162];

    auto g_xy_y_yy_xy = buffer_dpdd[163];

    auto g_xy_y_yy_xz = buffer_dpdd[164];

    auto g_xy_y_yy_yy = buffer_dpdd[165];

    auto g_xy_y_yy_yz = buffer_dpdd[166];

    auto g_xy_y_yy_zz = buffer_dpdd[167];

    auto g_xy_y_yz_xx = buffer_dpdd[168];

    auto g_xy_y_yz_xy = buffer_dpdd[169];

    auto g_xy_y_yz_xz = buffer_dpdd[170];

    auto g_xy_y_yz_yy = buffer_dpdd[171];

    auto g_xy_y_yz_yz = buffer_dpdd[172];

    auto g_xy_y_yz_zz = buffer_dpdd[173];

    auto g_xy_y_zz_xx = buffer_dpdd[174];

    auto g_xy_y_zz_xy = buffer_dpdd[175];

    auto g_xy_y_zz_xz = buffer_dpdd[176];

    auto g_xy_y_zz_yy = buffer_dpdd[177];

    auto g_xy_y_zz_yz = buffer_dpdd[178];

    auto g_xy_y_zz_zz = buffer_dpdd[179];

    auto g_xy_z_xx_xx = buffer_dpdd[180];

    auto g_xy_z_xx_xy = buffer_dpdd[181];

    auto g_xy_z_xx_xz = buffer_dpdd[182];

    auto g_xy_z_xx_yy = buffer_dpdd[183];

    auto g_xy_z_xx_yz = buffer_dpdd[184];

    auto g_xy_z_xx_zz = buffer_dpdd[185];

    auto g_xy_z_xy_xx = buffer_dpdd[186];

    auto g_xy_z_xy_xy = buffer_dpdd[187];

    auto g_xy_z_xy_xz = buffer_dpdd[188];

    auto g_xy_z_xy_yy = buffer_dpdd[189];

    auto g_xy_z_xy_yz = buffer_dpdd[190];

    auto g_xy_z_xy_zz = buffer_dpdd[191];

    auto g_xy_z_xz_xx = buffer_dpdd[192];

    auto g_xy_z_xz_xy = buffer_dpdd[193];

    auto g_xy_z_xz_xz = buffer_dpdd[194];

    auto g_xy_z_xz_yy = buffer_dpdd[195];

    auto g_xy_z_xz_yz = buffer_dpdd[196];

    auto g_xy_z_xz_zz = buffer_dpdd[197];

    auto g_xy_z_yy_xx = buffer_dpdd[198];

    auto g_xy_z_yy_xy = buffer_dpdd[199];

    auto g_xy_z_yy_xz = buffer_dpdd[200];

    auto g_xy_z_yy_yy = buffer_dpdd[201];

    auto g_xy_z_yy_yz = buffer_dpdd[202];

    auto g_xy_z_yy_zz = buffer_dpdd[203];

    auto g_xy_z_yz_xx = buffer_dpdd[204];

    auto g_xy_z_yz_xy = buffer_dpdd[205];

    auto g_xy_z_yz_xz = buffer_dpdd[206];

    auto g_xy_z_yz_yy = buffer_dpdd[207];

    auto g_xy_z_yz_yz = buffer_dpdd[208];

    auto g_xy_z_yz_zz = buffer_dpdd[209];

    auto g_xy_z_zz_xx = buffer_dpdd[210];

    auto g_xy_z_zz_xy = buffer_dpdd[211];

    auto g_xy_z_zz_xz = buffer_dpdd[212];

    auto g_xy_z_zz_yy = buffer_dpdd[213];

    auto g_xy_z_zz_yz = buffer_dpdd[214];

    auto g_xy_z_zz_zz = buffer_dpdd[215];

    auto g_xz_x_xx_xx = buffer_dpdd[216];

    auto g_xz_x_xx_xy = buffer_dpdd[217];

    auto g_xz_x_xx_xz = buffer_dpdd[218];

    auto g_xz_x_xx_yy = buffer_dpdd[219];

    auto g_xz_x_xx_yz = buffer_dpdd[220];

    auto g_xz_x_xx_zz = buffer_dpdd[221];

    auto g_xz_x_xy_xx = buffer_dpdd[222];

    auto g_xz_x_xy_xy = buffer_dpdd[223];

    auto g_xz_x_xy_xz = buffer_dpdd[224];

    auto g_xz_x_xy_yy = buffer_dpdd[225];

    auto g_xz_x_xy_yz = buffer_dpdd[226];

    auto g_xz_x_xy_zz = buffer_dpdd[227];

    auto g_xz_x_xz_xx = buffer_dpdd[228];

    auto g_xz_x_xz_xy = buffer_dpdd[229];

    auto g_xz_x_xz_xz = buffer_dpdd[230];

    auto g_xz_x_xz_yy = buffer_dpdd[231];

    auto g_xz_x_xz_yz = buffer_dpdd[232];

    auto g_xz_x_xz_zz = buffer_dpdd[233];

    auto g_xz_x_yy_xx = buffer_dpdd[234];

    auto g_xz_x_yy_xy = buffer_dpdd[235];

    auto g_xz_x_yy_xz = buffer_dpdd[236];

    auto g_xz_x_yy_yy = buffer_dpdd[237];

    auto g_xz_x_yy_yz = buffer_dpdd[238];

    auto g_xz_x_yy_zz = buffer_dpdd[239];

    auto g_xz_x_yz_xx = buffer_dpdd[240];

    auto g_xz_x_yz_xy = buffer_dpdd[241];

    auto g_xz_x_yz_xz = buffer_dpdd[242];

    auto g_xz_x_yz_yy = buffer_dpdd[243];

    auto g_xz_x_yz_yz = buffer_dpdd[244];

    auto g_xz_x_yz_zz = buffer_dpdd[245];

    auto g_xz_x_zz_xx = buffer_dpdd[246];

    auto g_xz_x_zz_xy = buffer_dpdd[247];

    auto g_xz_x_zz_xz = buffer_dpdd[248];

    auto g_xz_x_zz_yy = buffer_dpdd[249];

    auto g_xz_x_zz_yz = buffer_dpdd[250];

    auto g_xz_x_zz_zz = buffer_dpdd[251];

    auto g_xz_y_xx_xx = buffer_dpdd[252];

    auto g_xz_y_xx_xy = buffer_dpdd[253];

    auto g_xz_y_xx_xz = buffer_dpdd[254];

    auto g_xz_y_xx_yy = buffer_dpdd[255];

    auto g_xz_y_xx_yz = buffer_dpdd[256];

    auto g_xz_y_xx_zz = buffer_dpdd[257];

    auto g_xz_y_xy_xx = buffer_dpdd[258];

    auto g_xz_y_xy_xy = buffer_dpdd[259];

    auto g_xz_y_xy_xz = buffer_dpdd[260];

    auto g_xz_y_xy_yy = buffer_dpdd[261];

    auto g_xz_y_xy_yz = buffer_dpdd[262];

    auto g_xz_y_xy_zz = buffer_dpdd[263];

    auto g_xz_y_xz_xx = buffer_dpdd[264];

    auto g_xz_y_xz_xy = buffer_dpdd[265];

    auto g_xz_y_xz_xz = buffer_dpdd[266];

    auto g_xz_y_xz_yy = buffer_dpdd[267];

    auto g_xz_y_xz_yz = buffer_dpdd[268];

    auto g_xz_y_xz_zz = buffer_dpdd[269];

    auto g_xz_y_yy_xx = buffer_dpdd[270];

    auto g_xz_y_yy_xy = buffer_dpdd[271];

    auto g_xz_y_yy_xz = buffer_dpdd[272];

    auto g_xz_y_yy_yy = buffer_dpdd[273];

    auto g_xz_y_yy_yz = buffer_dpdd[274];

    auto g_xz_y_yy_zz = buffer_dpdd[275];

    auto g_xz_y_yz_xx = buffer_dpdd[276];

    auto g_xz_y_yz_xy = buffer_dpdd[277];

    auto g_xz_y_yz_xz = buffer_dpdd[278];

    auto g_xz_y_yz_yy = buffer_dpdd[279];

    auto g_xz_y_yz_yz = buffer_dpdd[280];

    auto g_xz_y_yz_zz = buffer_dpdd[281];

    auto g_xz_y_zz_xx = buffer_dpdd[282];

    auto g_xz_y_zz_xy = buffer_dpdd[283];

    auto g_xz_y_zz_xz = buffer_dpdd[284];

    auto g_xz_y_zz_yy = buffer_dpdd[285];

    auto g_xz_y_zz_yz = buffer_dpdd[286];

    auto g_xz_y_zz_zz = buffer_dpdd[287];

    auto g_xz_z_xx_xx = buffer_dpdd[288];

    auto g_xz_z_xx_xy = buffer_dpdd[289];

    auto g_xz_z_xx_xz = buffer_dpdd[290];

    auto g_xz_z_xx_yy = buffer_dpdd[291];

    auto g_xz_z_xx_yz = buffer_dpdd[292];

    auto g_xz_z_xx_zz = buffer_dpdd[293];

    auto g_xz_z_xy_xx = buffer_dpdd[294];

    auto g_xz_z_xy_xy = buffer_dpdd[295];

    auto g_xz_z_xy_xz = buffer_dpdd[296];

    auto g_xz_z_xy_yy = buffer_dpdd[297];

    auto g_xz_z_xy_yz = buffer_dpdd[298];

    auto g_xz_z_xy_zz = buffer_dpdd[299];

    auto g_xz_z_xz_xx = buffer_dpdd[300];

    auto g_xz_z_xz_xy = buffer_dpdd[301];

    auto g_xz_z_xz_xz = buffer_dpdd[302];

    auto g_xz_z_xz_yy = buffer_dpdd[303];

    auto g_xz_z_xz_yz = buffer_dpdd[304];

    auto g_xz_z_xz_zz = buffer_dpdd[305];

    auto g_xz_z_yy_xx = buffer_dpdd[306];

    auto g_xz_z_yy_xy = buffer_dpdd[307];

    auto g_xz_z_yy_xz = buffer_dpdd[308];

    auto g_xz_z_yy_yy = buffer_dpdd[309];

    auto g_xz_z_yy_yz = buffer_dpdd[310];

    auto g_xz_z_yy_zz = buffer_dpdd[311];

    auto g_xz_z_yz_xx = buffer_dpdd[312];

    auto g_xz_z_yz_xy = buffer_dpdd[313];

    auto g_xz_z_yz_xz = buffer_dpdd[314];

    auto g_xz_z_yz_yy = buffer_dpdd[315];

    auto g_xz_z_yz_yz = buffer_dpdd[316];

    auto g_xz_z_yz_zz = buffer_dpdd[317];

    auto g_xz_z_zz_xx = buffer_dpdd[318];

    auto g_xz_z_zz_xy = buffer_dpdd[319];

    auto g_xz_z_zz_xz = buffer_dpdd[320];

    auto g_xz_z_zz_yy = buffer_dpdd[321];

    auto g_xz_z_zz_yz = buffer_dpdd[322];

    auto g_xz_z_zz_zz = buffer_dpdd[323];

    auto g_yy_x_xx_xx = buffer_dpdd[324];

    auto g_yy_x_xx_xy = buffer_dpdd[325];

    auto g_yy_x_xx_xz = buffer_dpdd[326];

    auto g_yy_x_xx_yy = buffer_dpdd[327];

    auto g_yy_x_xx_yz = buffer_dpdd[328];

    auto g_yy_x_xx_zz = buffer_dpdd[329];

    auto g_yy_x_xy_xx = buffer_dpdd[330];

    auto g_yy_x_xy_xy = buffer_dpdd[331];

    auto g_yy_x_xy_xz = buffer_dpdd[332];

    auto g_yy_x_xy_yy = buffer_dpdd[333];

    auto g_yy_x_xy_yz = buffer_dpdd[334];

    auto g_yy_x_xy_zz = buffer_dpdd[335];

    auto g_yy_x_xz_xx = buffer_dpdd[336];

    auto g_yy_x_xz_xy = buffer_dpdd[337];

    auto g_yy_x_xz_xz = buffer_dpdd[338];

    auto g_yy_x_xz_yy = buffer_dpdd[339];

    auto g_yy_x_xz_yz = buffer_dpdd[340];

    auto g_yy_x_xz_zz = buffer_dpdd[341];

    auto g_yy_x_yy_xx = buffer_dpdd[342];

    auto g_yy_x_yy_xy = buffer_dpdd[343];

    auto g_yy_x_yy_xz = buffer_dpdd[344];

    auto g_yy_x_yy_yy = buffer_dpdd[345];

    auto g_yy_x_yy_yz = buffer_dpdd[346];

    auto g_yy_x_yy_zz = buffer_dpdd[347];

    auto g_yy_x_yz_xx = buffer_dpdd[348];

    auto g_yy_x_yz_xy = buffer_dpdd[349];

    auto g_yy_x_yz_xz = buffer_dpdd[350];

    auto g_yy_x_yz_yy = buffer_dpdd[351];

    auto g_yy_x_yz_yz = buffer_dpdd[352];

    auto g_yy_x_yz_zz = buffer_dpdd[353];

    auto g_yy_x_zz_xx = buffer_dpdd[354];

    auto g_yy_x_zz_xy = buffer_dpdd[355];

    auto g_yy_x_zz_xz = buffer_dpdd[356];

    auto g_yy_x_zz_yy = buffer_dpdd[357];

    auto g_yy_x_zz_yz = buffer_dpdd[358];

    auto g_yy_x_zz_zz = buffer_dpdd[359];

    auto g_yy_y_xx_xx = buffer_dpdd[360];

    auto g_yy_y_xx_xy = buffer_dpdd[361];

    auto g_yy_y_xx_xz = buffer_dpdd[362];

    auto g_yy_y_xx_yy = buffer_dpdd[363];

    auto g_yy_y_xx_yz = buffer_dpdd[364];

    auto g_yy_y_xx_zz = buffer_dpdd[365];

    auto g_yy_y_xy_xx = buffer_dpdd[366];

    auto g_yy_y_xy_xy = buffer_dpdd[367];

    auto g_yy_y_xy_xz = buffer_dpdd[368];

    auto g_yy_y_xy_yy = buffer_dpdd[369];

    auto g_yy_y_xy_yz = buffer_dpdd[370];

    auto g_yy_y_xy_zz = buffer_dpdd[371];

    auto g_yy_y_xz_xx = buffer_dpdd[372];

    auto g_yy_y_xz_xy = buffer_dpdd[373];

    auto g_yy_y_xz_xz = buffer_dpdd[374];

    auto g_yy_y_xz_yy = buffer_dpdd[375];

    auto g_yy_y_xz_yz = buffer_dpdd[376];

    auto g_yy_y_xz_zz = buffer_dpdd[377];

    auto g_yy_y_yy_xx = buffer_dpdd[378];

    auto g_yy_y_yy_xy = buffer_dpdd[379];

    auto g_yy_y_yy_xz = buffer_dpdd[380];

    auto g_yy_y_yy_yy = buffer_dpdd[381];

    auto g_yy_y_yy_yz = buffer_dpdd[382];

    auto g_yy_y_yy_zz = buffer_dpdd[383];

    auto g_yy_y_yz_xx = buffer_dpdd[384];

    auto g_yy_y_yz_xy = buffer_dpdd[385];

    auto g_yy_y_yz_xz = buffer_dpdd[386];

    auto g_yy_y_yz_yy = buffer_dpdd[387];

    auto g_yy_y_yz_yz = buffer_dpdd[388];

    auto g_yy_y_yz_zz = buffer_dpdd[389];

    auto g_yy_y_zz_xx = buffer_dpdd[390];

    auto g_yy_y_zz_xy = buffer_dpdd[391];

    auto g_yy_y_zz_xz = buffer_dpdd[392];

    auto g_yy_y_zz_yy = buffer_dpdd[393];

    auto g_yy_y_zz_yz = buffer_dpdd[394];

    auto g_yy_y_zz_zz = buffer_dpdd[395];

    auto g_yy_z_xx_xx = buffer_dpdd[396];

    auto g_yy_z_xx_xy = buffer_dpdd[397];

    auto g_yy_z_xx_xz = buffer_dpdd[398];

    auto g_yy_z_xx_yy = buffer_dpdd[399];

    auto g_yy_z_xx_yz = buffer_dpdd[400];

    auto g_yy_z_xx_zz = buffer_dpdd[401];

    auto g_yy_z_xy_xx = buffer_dpdd[402];

    auto g_yy_z_xy_xy = buffer_dpdd[403];

    auto g_yy_z_xy_xz = buffer_dpdd[404];

    auto g_yy_z_xy_yy = buffer_dpdd[405];

    auto g_yy_z_xy_yz = buffer_dpdd[406];

    auto g_yy_z_xy_zz = buffer_dpdd[407];

    auto g_yy_z_xz_xx = buffer_dpdd[408];

    auto g_yy_z_xz_xy = buffer_dpdd[409];

    auto g_yy_z_xz_xz = buffer_dpdd[410];

    auto g_yy_z_xz_yy = buffer_dpdd[411];

    auto g_yy_z_xz_yz = buffer_dpdd[412];

    auto g_yy_z_xz_zz = buffer_dpdd[413];

    auto g_yy_z_yy_xx = buffer_dpdd[414];

    auto g_yy_z_yy_xy = buffer_dpdd[415];

    auto g_yy_z_yy_xz = buffer_dpdd[416];

    auto g_yy_z_yy_yy = buffer_dpdd[417];

    auto g_yy_z_yy_yz = buffer_dpdd[418];

    auto g_yy_z_yy_zz = buffer_dpdd[419];

    auto g_yy_z_yz_xx = buffer_dpdd[420];

    auto g_yy_z_yz_xy = buffer_dpdd[421];

    auto g_yy_z_yz_xz = buffer_dpdd[422];

    auto g_yy_z_yz_yy = buffer_dpdd[423];

    auto g_yy_z_yz_yz = buffer_dpdd[424];

    auto g_yy_z_yz_zz = buffer_dpdd[425];

    auto g_yy_z_zz_xx = buffer_dpdd[426];

    auto g_yy_z_zz_xy = buffer_dpdd[427];

    auto g_yy_z_zz_xz = buffer_dpdd[428];

    auto g_yy_z_zz_yy = buffer_dpdd[429];

    auto g_yy_z_zz_yz = buffer_dpdd[430];

    auto g_yy_z_zz_zz = buffer_dpdd[431];

    auto g_yz_x_xx_xx = buffer_dpdd[432];

    auto g_yz_x_xx_xy = buffer_dpdd[433];

    auto g_yz_x_xx_xz = buffer_dpdd[434];

    auto g_yz_x_xx_yy = buffer_dpdd[435];

    auto g_yz_x_xx_yz = buffer_dpdd[436];

    auto g_yz_x_xx_zz = buffer_dpdd[437];

    auto g_yz_x_xy_xx = buffer_dpdd[438];

    auto g_yz_x_xy_xy = buffer_dpdd[439];

    auto g_yz_x_xy_xz = buffer_dpdd[440];

    auto g_yz_x_xy_yy = buffer_dpdd[441];

    auto g_yz_x_xy_yz = buffer_dpdd[442];

    auto g_yz_x_xy_zz = buffer_dpdd[443];

    auto g_yz_x_xz_xx = buffer_dpdd[444];

    auto g_yz_x_xz_xy = buffer_dpdd[445];

    auto g_yz_x_xz_xz = buffer_dpdd[446];

    auto g_yz_x_xz_yy = buffer_dpdd[447];

    auto g_yz_x_xz_yz = buffer_dpdd[448];

    auto g_yz_x_xz_zz = buffer_dpdd[449];

    auto g_yz_x_yy_xx = buffer_dpdd[450];

    auto g_yz_x_yy_xy = buffer_dpdd[451];

    auto g_yz_x_yy_xz = buffer_dpdd[452];

    auto g_yz_x_yy_yy = buffer_dpdd[453];

    auto g_yz_x_yy_yz = buffer_dpdd[454];

    auto g_yz_x_yy_zz = buffer_dpdd[455];

    auto g_yz_x_yz_xx = buffer_dpdd[456];

    auto g_yz_x_yz_xy = buffer_dpdd[457];

    auto g_yz_x_yz_xz = buffer_dpdd[458];

    auto g_yz_x_yz_yy = buffer_dpdd[459];

    auto g_yz_x_yz_yz = buffer_dpdd[460];

    auto g_yz_x_yz_zz = buffer_dpdd[461];

    auto g_yz_x_zz_xx = buffer_dpdd[462];

    auto g_yz_x_zz_xy = buffer_dpdd[463];

    auto g_yz_x_zz_xz = buffer_dpdd[464];

    auto g_yz_x_zz_yy = buffer_dpdd[465];

    auto g_yz_x_zz_yz = buffer_dpdd[466];

    auto g_yz_x_zz_zz = buffer_dpdd[467];

    auto g_yz_y_xx_xx = buffer_dpdd[468];

    auto g_yz_y_xx_xy = buffer_dpdd[469];

    auto g_yz_y_xx_xz = buffer_dpdd[470];

    auto g_yz_y_xx_yy = buffer_dpdd[471];

    auto g_yz_y_xx_yz = buffer_dpdd[472];

    auto g_yz_y_xx_zz = buffer_dpdd[473];

    auto g_yz_y_xy_xx = buffer_dpdd[474];

    auto g_yz_y_xy_xy = buffer_dpdd[475];

    auto g_yz_y_xy_xz = buffer_dpdd[476];

    auto g_yz_y_xy_yy = buffer_dpdd[477];

    auto g_yz_y_xy_yz = buffer_dpdd[478];

    auto g_yz_y_xy_zz = buffer_dpdd[479];

    auto g_yz_y_xz_xx = buffer_dpdd[480];

    auto g_yz_y_xz_xy = buffer_dpdd[481];

    auto g_yz_y_xz_xz = buffer_dpdd[482];

    auto g_yz_y_xz_yy = buffer_dpdd[483];

    auto g_yz_y_xz_yz = buffer_dpdd[484];

    auto g_yz_y_xz_zz = buffer_dpdd[485];

    auto g_yz_y_yy_xx = buffer_dpdd[486];

    auto g_yz_y_yy_xy = buffer_dpdd[487];

    auto g_yz_y_yy_xz = buffer_dpdd[488];

    auto g_yz_y_yy_yy = buffer_dpdd[489];

    auto g_yz_y_yy_yz = buffer_dpdd[490];

    auto g_yz_y_yy_zz = buffer_dpdd[491];

    auto g_yz_y_yz_xx = buffer_dpdd[492];

    auto g_yz_y_yz_xy = buffer_dpdd[493];

    auto g_yz_y_yz_xz = buffer_dpdd[494];

    auto g_yz_y_yz_yy = buffer_dpdd[495];

    auto g_yz_y_yz_yz = buffer_dpdd[496];

    auto g_yz_y_yz_zz = buffer_dpdd[497];

    auto g_yz_y_zz_xx = buffer_dpdd[498];

    auto g_yz_y_zz_xy = buffer_dpdd[499];

    auto g_yz_y_zz_xz = buffer_dpdd[500];

    auto g_yz_y_zz_yy = buffer_dpdd[501];

    auto g_yz_y_zz_yz = buffer_dpdd[502];

    auto g_yz_y_zz_zz = buffer_dpdd[503];

    auto g_yz_z_xx_xx = buffer_dpdd[504];

    auto g_yz_z_xx_xy = buffer_dpdd[505];

    auto g_yz_z_xx_xz = buffer_dpdd[506];

    auto g_yz_z_xx_yy = buffer_dpdd[507];

    auto g_yz_z_xx_yz = buffer_dpdd[508];

    auto g_yz_z_xx_zz = buffer_dpdd[509];

    auto g_yz_z_xy_xx = buffer_dpdd[510];

    auto g_yz_z_xy_xy = buffer_dpdd[511];

    auto g_yz_z_xy_xz = buffer_dpdd[512];

    auto g_yz_z_xy_yy = buffer_dpdd[513];

    auto g_yz_z_xy_yz = buffer_dpdd[514];

    auto g_yz_z_xy_zz = buffer_dpdd[515];

    auto g_yz_z_xz_xx = buffer_dpdd[516];

    auto g_yz_z_xz_xy = buffer_dpdd[517];

    auto g_yz_z_xz_xz = buffer_dpdd[518];

    auto g_yz_z_xz_yy = buffer_dpdd[519];

    auto g_yz_z_xz_yz = buffer_dpdd[520];

    auto g_yz_z_xz_zz = buffer_dpdd[521];

    auto g_yz_z_yy_xx = buffer_dpdd[522];

    auto g_yz_z_yy_xy = buffer_dpdd[523];

    auto g_yz_z_yy_xz = buffer_dpdd[524];

    auto g_yz_z_yy_yy = buffer_dpdd[525];

    auto g_yz_z_yy_yz = buffer_dpdd[526];

    auto g_yz_z_yy_zz = buffer_dpdd[527];

    auto g_yz_z_yz_xx = buffer_dpdd[528];

    auto g_yz_z_yz_xy = buffer_dpdd[529];

    auto g_yz_z_yz_xz = buffer_dpdd[530];

    auto g_yz_z_yz_yy = buffer_dpdd[531];

    auto g_yz_z_yz_yz = buffer_dpdd[532];

    auto g_yz_z_yz_zz = buffer_dpdd[533];

    auto g_yz_z_zz_xx = buffer_dpdd[534];

    auto g_yz_z_zz_xy = buffer_dpdd[535];

    auto g_yz_z_zz_xz = buffer_dpdd[536];

    auto g_yz_z_zz_yy = buffer_dpdd[537];

    auto g_yz_z_zz_yz = buffer_dpdd[538];

    auto g_yz_z_zz_zz = buffer_dpdd[539];

    auto g_zz_x_xx_xx = buffer_dpdd[540];

    auto g_zz_x_xx_xy = buffer_dpdd[541];

    auto g_zz_x_xx_xz = buffer_dpdd[542];

    auto g_zz_x_xx_yy = buffer_dpdd[543];

    auto g_zz_x_xx_yz = buffer_dpdd[544];

    auto g_zz_x_xx_zz = buffer_dpdd[545];

    auto g_zz_x_xy_xx = buffer_dpdd[546];

    auto g_zz_x_xy_xy = buffer_dpdd[547];

    auto g_zz_x_xy_xz = buffer_dpdd[548];

    auto g_zz_x_xy_yy = buffer_dpdd[549];

    auto g_zz_x_xy_yz = buffer_dpdd[550];

    auto g_zz_x_xy_zz = buffer_dpdd[551];

    auto g_zz_x_xz_xx = buffer_dpdd[552];

    auto g_zz_x_xz_xy = buffer_dpdd[553];

    auto g_zz_x_xz_xz = buffer_dpdd[554];

    auto g_zz_x_xz_yy = buffer_dpdd[555];

    auto g_zz_x_xz_yz = buffer_dpdd[556];

    auto g_zz_x_xz_zz = buffer_dpdd[557];

    auto g_zz_x_yy_xx = buffer_dpdd[558];

    auto g_zz_x_yy_xy = buffer_dpdd[559];

    auto g_zz_x_yy_xz = buffer_dpdd[560];

    auto g_zz_x_yy_yy = buffer_dpdd[561];

    auto g_zz_x_yy_yz = buffer_dpdd[562];

    auto g_zz_x_yy_zz = buffer_dpdd[563];

    auto g_zz_x_yz_xx = buffer_dpdd[564];

    auto g_zz_x_yz_xy = buffer_dpdd[565];

    auto g_zz_x_yz_xz = buffer_dpdd[566];

    auto g_zz_x_yz_yy = buffer_dpdd[567];

    auto g_zz_x_yz_yz = buffer_dpdd[568];

    auto g_zz_x_yz_zz = buffer_dpdd[569];

    auto g_zz_x_zz_xx = buffer_dpdd[570];

    auto g_zz_x_zz_xy = buffer_dpdd[571];

    auto g_zz_x_zz_xz = buffer_dpdd[572];

    auto g_zz_x_zz_yy = buffer_dpdd[573];

    auto g_zz_x_zz_yz = buffer_dpdd[574];

    auto g_zz_x_zz_zz = buffer_dpdd[575];

    auto g_zz_y_xx_xx = buffer_dpdd[576];

    auto g_zz_y_xx_xy = buffer_dpdd[577];

    auto g_zz_y_xx_xz = buffer_dpdd[578];

    auto g_zz_y_xx_yy = buffer_dpdd[579];

    auto g_zz_y_xx_yz = buffer_dpdd[580];

    auto g_zz_y_xx_zz = buffer_dpdd[581];

    auto g_zz_y_xy_xx = buffer_dpdd[582];

    auto g_zz_y_xy_xy = buffer_dpdd[583];

    auto g_zz_y_xy_xz = buffer_dpdd[584];

    auto g_zz_y_xy_yy = buffer_dpdd[585];

    auto g_zz_y_xy_yz = buffer_dpdd[586];

    auto g_zz_y_xy_zz = buffer_dpdd[587];

    auto g_zz_y_xz_xx = buffer_dpdd[588];

    auto g_zz_y_xz_xy = buffer_dpdd[589];

    auto g_zz_y_xz_xz = buffer_dpdd[590];

    auto g_zz_y_xz_yy = buffer_dpdd[591];

    auto g_zz_y_xz_yz = buffer_dpdd[592];

    auto g_zz_y_xz_zz = buffer_dpdd[593];

    auto g_zz_y_yy_xx = buffer_dpdd[594];

    auto g_zz_y_yy_xy = buffer_dpdd[595];

    auto g_zz_y_yy_xz = buffer_dpdd[596];

    auto g_zz_y_yy_yy = buffer_dpdd[597];

    auto g_zz_y_yy_yz = buffer_dpdd[598];

    auto g_zz_y_yy_zz = buffer_dpdd[599];

    auto g_zz_y_yz_xx = buffer_dpdd[600];

    auto g_zz_y_yz_xy = buffer_dpdd[601];

    auto g_zz_y_yz_xz = buffer_dpdd[602];

    auto g_zz_y_yz_yy = buffer_dpdd[603];

    auto g_zz_y_yz_yz = buffer_dpdd[604];

    auto g_zz_y_yz_zz = buffer_dpdd[605];

    auto g_zz_y_zz_xx = buffer_dpdd[606];

    auto g_zz_y_zz_xy = buffer_dpdd[607];

    auto g_zz_y_zz_xz = buffer_dpdd[608];

    auto g_zz_y_zz_yy = buffer_dpdd[609];

    auto g_zz_y_zz_yz = buffer_dpdd[610];

    auto g_zz_y_zz_zz = buffer_dpdd[611];

    auto g_zz_z_xx_xx = buffer_dpdd[612];

    auto g_zz_z_xx_xy = buffer_dpdd[613];

    auto g_zz_z_xx_xz = buffer_dpdd[614];

    auto g_zz_z_xx_yy = buffer_dpdd[615];

    auto g_zz_z_xx_yz = buffer_dpdd[616];

    auto g_zz_z_xx_zz = buffer_dpdd[617];

    auto g_zz_z_xy_xx = buffer_dpdd[618];

    auto g_zz_z_xy_xy = buffer_dpdd[619];

    auto g_zz_z_xy_xz = buffer_dpdd[620];

    auto g_zz_z_xy_yy = buffer_dpdd[621];

    auto g_zz_z_xy_yz = buffer_dpdd[622];

    auto g_zz_z_xy_zz = buffer_dpdd[623];

    auto g_zz_z_xz_xx = buffer_dpdd[624];

    auto g_zz_z_xz_xy = buffer_dpdd[625];

    auto g_zz_z_xz_xz = buffer_dpdd[626];

    auto g_zz_z_xz_yy = buffer_dpdd[627];

    auto g_zz_z_xz_yz = buffer_dpdd[628];

    auto g_zz_z_xz_zz = buffer_dpdd[629];

    auto g_zz_z_yy_xx = buffer_dpdd[630];

    auto g_zz_z_yy_xy = buffer_dpdd[631];

    auto g_zz_z_yy_xz = buffer_dpdd[632];

    auto g_zz_z_yy_yy = buffer_dpdd[633];

    auto g_zz_z_yy_yz = buffer_dpdd[634];

    auto g_zz_z_yy_zz = buffer_dpdd[635];

    auto g_zz_z_yz_xx = buffer_dpdd[636];

    auto g_zz_z_yz_xy = buffer_dpdd[637];

    auto g_zz_z_yz_xz = buffer_dpdd[638];

    auto g_zz_z_yz_yy = buffer_dpdd[639];

    auto g_zz_z_yz_yz = buffer_dpdd[640];

    auto g_zz_z_yz_zz = buffer_dpdd[641];

    auto g_zz_z_zz_xx = buffer_dpdd[642];

    auto g_zz_z_zz_xy = buffer_dpdd[643];

    auto g_zz_z_zz_xz = buffer_dpdd[644];

    auto g_zz_z_zz_yy = buffer_dpdd[645];

    auto g_zz_z_zz_yz = buffer_dpdd[646];

    auto g_zz_z_zz_zz = buffer_dpdd[647];

    /// Set up components of integrals buffer : buffer_1010_pppd

    auto g_x_0_x_0_x_x_x_xx = buffer_1010_pppd[0];

    auto g_x_0_x_0_x_x_x_xy = buffer_1010_pppd[1];

    auto g_x_0_x_0_x_x_x_xz = buffer_1010_pppd[2];

    auto g_x_0_x_0_x_x_x_yy = buffer_1010_pppd[3];

    auto g_x_0_x_0_x_x_x_yz = buffer_1010_pppd[4];

    auto g_x_0_x_0_x_x_x_zz = buffer_1010_pppd[5];

    auto g_x_0_x_0_x_x_y_xx = buffer_1010_pppd[6];

    auto g_x_0_x_0_x_x_y_xy = buffer_1010_pppd[7];

    auto g_x_0_x_0_x_x_y_xz = buffer_1010_pppd[8];

    auto g_x_0_x_0_x_x_y_yy = buffer_1010_pppd[9];

    auto g_x_0_x_0_x_x_y_yz = buffer_1010_pppd[10];

    auto g_x_0_x_0_x_x_y_zz = buffer_1010_pppd[11];

    auto g_x_0_x_0_x_x_z_xx = buffer_1010_pppd[12];

    auto g_x_0_x_0_x_x_z_xy = buffer_1010_pppd[13];

    auto g_x_0_x_0_x_x_z_xz = buffer_1010_pppd[14];

    auto g_x_0_x_0_x_x_z_yy = buffer_1010_pppd[15];

    auto g_x_0_x_0_x_x_z_yz = buffer_1010_pppd[16];

    auto g_x_0_x_0_x_x_z_zz = buffer_1010_pppd[17];

    auto g_x_0_x_0_x_y_x_xx = buffer_1010_pppd[18];

    auto g_x_0_x_0_x_y_x_xy = buffer_1010_pppd[19];

    auto g_x_0_x_0_x_y_x_xz = buffer_1010_pppd[20];

    auto g_x_0_x_0_x_y_x_yy = buffer_1010_pppd[21];

    auto g_x_0_x_0_x_y_x_yz = buffer_1010_pppd[22];

    auto g_x_0_x_0_x_y_x_zz = buffer_1010_pppd[23];

    auto g_x_0_x_0_x_y_y_xx = buffer_1010_pppd[24];

    auto g_x_0_x_0_x_y_y_xy = buffer_1010_pppd[25];

    auto g_x_0_x_0_x_y_y_xz = buffer_1010_pppd[26];

    auto g_x_0_x_0_x_y_y_yy = buffer_1010_pppd[27];

    auto g_x_0_x_0_x_y_y_yz = buffer_1010_pppd[28];

    auto g_x_0_x_0_x_y_y_zz = buffer_1010_pppd[29];

    auto g_x_0_x_0_x_y_z_xx = buffer_1010_pppd[30];

    auto g_x_0_x_0_x_y_z_xy = buffer_1010_pppd[31];

    auto g_x_0_x_0_x_y_z_xz = buffer_1010_pppd[32];

    auto g_x_0_x_0_x_y_z_yy = buffer_1010_pppd[33];

    auto g_x_0_x_0_x_y_z_yz = buffer_1010_pppd[34];

    auto g_x_0_x_0_x_y_z_zz = buffer_1010_pppd[35];

    auto g_x_0_x_0_x_z_x_xx = buffer_1010_pppd[36];

    auto g_x_0_x_0_x_z_x_xy = buffer_1010_pppd[37];

    auto g_x_0_x_0_x_z_x_xz = buffer_1010_pppd[38];

    auto g_x_0_x_0_x_z_x_yy = buffer_1010_pppd[39];

    auto g_x_0_x_0_x_z_x_yz = buffer_1010_pppd[40];

    auto g_x_0_x_0_x_z_x_zz = buffer_1010_pppd[41];

    auto g_x_0_x_0_x_z_y_xx = buffer_1010_pppd[42];

    auto g_x_0_x_0_x_z_y_xy = buffer_1010_pppd[43];

    auto g_x_0_x_0_x_z_y_xz = buffer_1010_pppd[44];

    auto g_x_0_x_0_x_z_y_yy = buffer_1010_pppd[45];

    auto g_x_0_x_0_x_z_y_yz = buffer_1010_pppd[46];

    auto g_x_0_x_0_x_z_y_zz = buffer_1010_pppd[47];

    auto g_x_0_x_0_x_z_z_xx = buffer_1010_pppd[48];

    auto g_x_0_x_0_x_z_z_xy = buffer_1010_pppd[49];

    auto g_x_0_x_0_x_z_z_xz = buffer_1010_pppd[50];

    auto g_x_0_x_0_x_z_z_yy = buffer_1010_pppd[51];

    auto g_x_0_x_0_x_z_z_yz = buffer_1010_pppd[52];

    auto g_x_0_x_0_x_z_z_zz = buffer_1010_pppd[53];

    auto g_x_0_x_0_y_x_x_xx = buffer_1010_pppd[54];

    auto g_x_0_x_0_y_x_x_xy = buffer_1010_pppd[55];

    auto g_x_0_x_0_y_x_x_xz = buffer_1010_pppd[56];

    auto g_x_0_x_0_y_x_x_yy = buffer_1010_pppd[57];

    auto g_x_0_x_0_y_x_x_yz = buffer_1010_pppd[58];

    auto g_x_0_x_0_y_x_x_zz = buffer_1010_pppd[59];

    auto g_x_0_x_0_y_x_y_xx = buffer_1010_pppd[60];

    auto g_x_0_x_0_y_x_y_xy = buffer_1010_pppd[61];

    auto g_x_0_x_0_y_x_y_xz = buffer_1010_pppd[62];

    auto g_x_0_x_0_y_x_y_yy = buffer_1010_pppd[63];

    auto g_x_0_x_0_y_x_y_yz = buffer_1010_pppd[64];

    auto g_x_0_x_0_y_x_y_zz = buffer_1010_pppd[65];

    auto g_x_0_x_0_y_x_z_xx = buffer_1010_pppd[66];

    auto g_x_0_x_0_y_x_z_xy = buffer_1010_pppd[67];

    auto g_x_0_x_0_y_x_z_xz = buffer_1010_pppd[68];

    auto g_x_0_x_0_y_x_z_yy = buffer_1010_pppd[69];

    auto g_x_0_x_0_y_x_z_yz = buffer_1010_pppd[70];

    auto g_x_0_x_0_y_x_z_zz = buffer_1010_pppd[71];

    auto g_x_0_x_0_y_y_x_xx = buffer_1010_pppd[72];

    auto g_x_0_x_0_y_y_x_xy = buffer_1010_pppd[73];

    auto g_x_0_x_0_y_y_x_xz = buffer_1010_pppd[74];

    auto g_x_0_x_0_y_y_x_yy = buffer_1010_pppd[75];

    auto g_x_0_x_0_y_y_x_yz = buffer_1010_pppd[76];

    auto g_x_0_x_0_y_y_x_zz = buffer_1010_pppd[77];

    auto g_x_0_x_0_y_y_y_xx = buffer_1010_pppd[78];

    auto g_x_0_x_0_y_y_y_xy = buffer_1010_pppd[79];

    auto g_x_0_x_0_y_y_y_xz = buffer_1010_pppd[80];

    auto g_x_0_x_0_y_y_y_yy = buffer_1010_pppd[81];

    auto g_x_0_x_0_y_y_y_yz = buffer_1010_pppd[82];

    auto g_x_0_x_0_y_y_y_zz = buffer_1010_pppd[83];

    auto g_x_0_x_0_y_y_z_xx = buffer_1010_pppd[84];

    auto g_x_0_x_0_y_y_z_xy = buffer_1010_pppd[85];

    auto g_x_0_x_0_y_y_z_xz = buffer_1010_pppd[86];

    auto g_x_0_x_0_y_y_z_yy = buffer_1010_pppd[87];

    auto g_x_0_x_0_y_y_z_yz = buffer_1010_pppd[88];

    auto g_x_0_x_0_y_y_z_zz = buffer_1010_pppd[89];

    auto g_x_0_x_0_y_z_x_xx = buffer_1010_pppd[90];

    auto g_x_0_x_0_y_z_x_xy = buffer_1010_pppd[91];

    auto g_x_0_x_0_y_z_x_xz = buffer_1010_pppd[92];

    auto g_x_0_x_0_y_z_x_yy = buffer_1010_pppd[93];

    auto g_x_0_x_0_y_z_x_yz = buffer_1010_pppd[94];

    auto g_x_0_x_0_y_z_x_zz = buffer_1010_pppd[95];

    auto g_x_0_x_0_y_z_y_xx = buffer_1010_pppd[96];

    auto g_x_0_x_0_y_z_y_xy = buffer_1010_pppd[97];

    auto g_x_0_x_0_y_z_y_xz = buffer_1010_pppd[98];

    auto g_x_0_x_0_y_z_y_yy = buffer_1010_pppd[99];

    auto g_x_0_x_0_y_z_y_yz = buffer_1010_pppd[100];

    auto g_x_0_x_0_y_z_y_zz = buffer_1010_pppd[101];

    auto g_x_0_x_0_y_z_z_xx = buffer_1010_pppd[102];

    auto g_x_0_x_0_y_z_z_xy = buffer_1010_pppd[103];

    auto g_x_0_x_0_y_z_z_xz = buffer_1010_pppd[104];

    auto g_x_0_x_0_y_z_z_yy = buffer_1010_pppd[105];

    auto g_x_0_x_0_y_z_z_yz = buffer_1010_pppd[106];

    auto g_x_0_x_0_y_z_z_zz = buffer_1010_pppd[107];

    auto g_x_0_x_0_z_x_x_xx = buffer_1010_pppd[108];

    auto g_x_0_x_0_z_x_x_xy = buffer_1010_pppd[109];

    auto g_x_0_x_0_z_x_x_xz = buffer_1010_pppd[110];

    auto g_x_0_x_0_z_x_x_yy = buffer_1010_pppd[111];

    auto g_x_0_x_0_z_x_x_yz = buffer_1010_pppd[112];

    auto g_x_0_x_0_z_x_x_zz = buffer_1010_pppd[113];

    auto g_x_0_x_0_z_x_y_xx = buffer_1010_pppd[114];

    auto g_x_0_x_0_z_x_y_xy = buffer_1010_pppd[115];

    auto g_x_0_x_0_z_x_y_xz = buffer_1010_pppd[116];

    auto g_x_0_x_0_z_x_y_yy = buffer_1010_pppd[117];

    auto g_x_0_x_0_z_x_y_yz = buffer_1010_pppd[118];

    auto g_x_0_x_0_z_x_y_zz = buffer_1010_pppd[119];

    auto g_x_0_x_0_z_x_z_xx = buffer_1010_pppd[120];

    auto g_x_0_x_0_z_x_z_xy = buffer_1010_pppd[121];

    auto g_x_0_x_0_z_x_z_xz = buffer_1010_pppd[122];

    auto g_x_0_x_0_z_x_z_yy = buffer_1010_pppd[123];

    auto g_x_0_x_0_z_x_z_yz = buffer_1010_pppd[124];

    auto g_x_0_x_0_z_x_z_zz = buffer_1010_pppd[125];

    auto g_x_0_x_0_z_y_x_xx = buffer_1010_pppd[126];

    auto g_x_0_x_0_z_y_x_xy = buffer_1010_pppd[127];

    auto g_x_0_x_0_z_y_x_xz = buffer_1010_pppd[128];

    auto g_x_0_x_0_z_y_x_yy = buffer_1010_pppd[129];

    auto g_x_0_x_0_z_y_x_yz = buffer_1010_pppd[130];

    auto g_x_0_x_0_z_y_x_zz = buffer_1010_pppd[131];

    auto g_x_0_x_0_z_y_y_xx = buffer_1010_pppd[132];

    auto g_x_0_x_0_z_y_y_xy = buffer_1010_pppd[133];

    auto g_x_0_x_0_z_y_y_xz = buffer_1010_pppd[134];

    auto g_x_0_x_0_z_y_y_yy = buffer_1010_pppd[135];

    auto g_x_0_x_0_z_y_y_yz = buffer_1010_pppd[136];

    auto g_x_0_x_0_z_y_y_zz = buffer_1010_pppd[137];

    auto g_x_0_x_0_z_y_z_xx = buffer_1010_pppd[138];

    auto g_x_0_x_0_z_y_z_xy = buffer_1010_pppd[139];

    auto g_x_0_x_0_z_y_z_xz = buffer_1010_pppd[140];

    auto g_x_0_x_0_z_y_z_yy = buffer_1010_pppd[141];

    auto g_x_0_x_0_z_y_z_yz = buffer_1010_pppd[142];

    auto g_x_0_x_0_z_y_z_zz = buffer_1010_pppd[143];

    auto g_x_0_x_0_z_z_x_xx = buffer_1010_pppd[144];

    auto g_x_0_x_0_z_z_x_xy = buffer_1010_pppd[145];

    auto g_x_0_x_0_z_z_x_xz = buffer_1010_pppd[146];

    auto g_x_0_x_0_z_z_x_yy = buffer_1010_pppd[147];

    auto g_x_0_x_0_z_z_x_yz = buffer_1010_pppd[148];

    auto g_x_0_x_0_z_z_x_zz = buffer_1010_pppd[149];

    auto g_x_0_x_0_z_z_y_xx = buffer_1010_pppd[150];

    auto g_x_0_x_0_z_z_y_xy = buffer_1010_pppd[151];

    auto g_x_0_x_0_z_z_y_xz = buffer_1010_pppd[152];

    auto g_x_0_x_0_z_z_y_yy = buffer_1010_pppd[153];

    auto g_x_0_x_0_z_z_y_yz = buffer_1010_pppd[154];

    auto g_x_0_x_0_z_z_y_zz = buffer_1010_pppd[155];

    auto g_x_0_x_0_z_z_z_xx = buffer_1010_pppd[156];

    auto g_x_0_x_0_z_z_z_xy = buffer_1010_pppd[157];

    auto g_x_0_x_0_z_z_z_xz = buffer_1010_pppd[158];

    auto g_x_0_x_0_z_z_z_yy = buffer_1010_pppd[159];

    auto g_x_0_x_0_z_z_z_yz = buffer_1010_pppd[160];

    auto g_x_0_x_0_z_z_z_zz = buffer_1010_pppd[161];

    auto g_x_0_y_0_x_x_x_xx = buffer_1010_pppd[162];

    auto g_x_0_y_0_x_x_x_xy = buffer_1010_pppd[163];

    auto g_x_0_y_0_x_x_x_xz = buffer_1010_pppd[164];

    auto g_x_0_y_0_x_x_x_yy = buffer_1010_pppd[165];

    auto g_x_0_y_0_x_x_x_yz = buffer_1010_pppd[166];

    auto g_x_0_y_0_x_x_x_zz = buffer_1010_pppd[167];

    auto g_x_0_y_0_x_x_y_xx = buffer_1010_pppd[168];

    auto g_x_0_y_0_x_x_y_xy = buffer_1010_pppd[169];

    auto g_x_0_y_0_x_x_y_xz = buffer_1010_pppd[170];

    auto g_x_0_y_0_x_x_y_yy = buffer_1010_pppd[171];

    auto g_x_0_y_0_x_x_y_yz = buffer_1010_pppd[172];

    auto g_x_0_y_0_x_x_y_zz = buffer_1010_pppd[173];

    auto g_x_0_y_0_x_x_z_xx = buffer_1010_pppd[174];

    auto g_x_0_y_0_x_x_z_xy = buffer_1010_pppd[175];

    auto g_x_0_y_0_x_x_z_xz = buffer_1010_pppd[176];

    auto g_x_0_y_0_x_x_z_yy = buffer_1010_pppd[177];

    auto g_x_0_y_0_x_x_z_yz = buffer_1010_pppd[178];

    auto g_x_0_y_0_x_x_z_zz = buffer_1010_pppd[179];

    auto g_x_0_y_0_x_y_x_xx = buffer_1010_pppd[180];

    auto g_x_0_y_0_x_y_x_xy = buffer_1010_pppd[181];

    auto g_x_0_y_0_x_y_x_xz = buffer_1010_pppd[182];

    auto g_x_0_y_0_x_y_x_yy = buffer_1010_pppd[183];

    auto g_x_0_y_0_x_y_x_yz = buffer_1010_pppd[184];

    auto g_x_0_y_0_x_y_x_zz = buffer_1010_pppd[185];

    auto g_x_0_y_0_x_y_y_xx = buffer_1010_pppd[186];

    auto g_x_0_y_0_x_y_y_xy = buffer_1010_pppd[187];

    auto g_x_0_y_0_x_y_y_xz = buffer_1010_pppd[188];

    auto g_x_0_y_0_x_y_y_yy = buffer_1010_pppd[189];

    auto g_x_0_y_0_x_y_y_yz = buffer_1010_pppd[190];

    auto g_x_0_y_0_x_y_y_zz = buffer_1010_pppd[191];

    auto g_x_0_y_0_x_y_z_xx = buffer_1010_pppd[192];

    auto g_x_0_y_0_x_y_z_xy = buffer_1010_pppd[193];

    auto g_x_0_y_0_x_y_z_xz = buffer_1010_pppd[194];

    auto g_x_0_y_0_x_y_z_yy = buffer_1010_pppd[195];

    auto g_x_0_y_0_x_y_z_yz = buffer_1010_pppd[196];

    auto g_x_0_y_0_x_y_z_zz = buffer_1010_pppd[197];

    auto g_x_0_y_0_x_z_x_xx = buffer_1010_pppd[198];

    auto g_x_0_y_0_x_z_x_xy = buffer_1010_pppd[199];

    auto g_x_0_y_0_x_z_x_xz = buffer_1010_pppd[200];

    auto g_x_0_y_0_x_z_x_yy = buffer_1010_pppd[201];

    auto g_x_0_y_0_x_z_x_yz = buffer_1010_pppd[202];

    auto g_x_0_y_0_x_z_x_zz = buffer_1010_pppd[203];

    auto g_x_0_y_0_x_z_y_xx = buffer_1010_pppd[204];

    auto g_x_0_y_0_x_z_y_xy = buffer_1010_pppd[205];

    auto g_x_0_y_0_x_z_y_xz = buffer_1010_pppd[206];

    auto g_x_0_y_0_x_z_y_yy = buffer_1010_pppd[207];

    auto g_x_0_y_0_x_z_y_yz = buffer_1010_pppd[208];

    auto g_x_0_y_0_x_z_y_zz = buffer_1010_pppd[209];

    auto g_x_0_y_0_x_z_z_xx = buffer_1010_pppd[210];

    auto g_x_0_y_0_x_z_z_xy = buffer_1010_pppd[211];

    auto g_x_0_y_0_x_z_z_xz = buffer_1010_pppd[212];

    auto g_x_0_y_0_x_z_z_yy = buffer_1010_pppd[213];

    auto g_x_0_y_0_x_z_z_yz = buffer_1010_pppd[214];

    auto g_x_0_y_0_x_z_z_zz = buffer_1010_pppd[215];

    auto g_x_0_y_0_y_x_x_xx = buffer_1010_pppd[216];

    auto g_x_0_y_0_y_x_x_xy = buffer_1010_pppd[217];

    auto g_x_0_y_0_y_x_x_xz = buffer_1010_pppd[218];

    auto g_x_0_y_0_y_x_x_yy = buffer_1010_pppd[219];

    auto g_x_0_y_0_y_x_x_yz = buffer_1010_pppd[220];

    auto g_x_0_y_0_y_x_x_zz = buffer_1010_pppd[221];

    auto g_x_0_y_0_y_x_y_xx = buffer_1010_pppd[222];

    auto g_x_0_y_0_y_x_y_xy = buffer_1010_pppd[223];

    auto g_x_0_y_0_y_x_y_xz = buffer_1010_pppd[224];

    auto g_x_0_y_0_y_x_y_yy = buffer_1010_pppd[225];

    auto g_x_0_y_0_y_x_y_yz = buffer_1010_pppd[226];

    auto g_x_0_y_0_y_x_y_zz = buffer_1010_pppd[227];

    auto g_x_0_y_0_y_x_z_xx = buffer_1010_pppd[228];

    auto g_x_0_y_0_y_x_z_xy = buffer_1010_pppd[229];

    auto g_x_0_y_0_y_x_z_xz = buffer_1010_pppd[230];

    auto g_x_0_y_0_y_x_z_yy = buffer_1010_pppd[231];

    auto g_x_0_y_0_y_x_z_yz = buffer_1010_pppd[232];

    auto g_x_0_y_0_y_x_z_zz = buffer_1010_pppd[233];

    auto g_x_0_y_0_y_y_x_xx = buffer_1010_pppd[234];

    auto g_x_0_y_0_y_y_x_xy = buffer_1010_pppd[235];

    auto g_x_0_y_0_y_y_x_xz = buffer_1010_pppd[236];

    auto g_x_0_y_0_y_y_x_yy = buffer_1010_pppd[237];

    auto g_x_0_y_0_y_y_x_yz = buffer_1010_pppd[238];

    auto g_x_0_y_0_y_y_x_zz = buffer_1010_pppd[239];

    auto g_x_0_y_0_y_y_y_xx = buffer_1010_pppd[240];

    auto g_x_0_y_0_y_y_y_xy = buffer_1010_pppd[241];

    auto g_x_0_y_0_y_y_y_xz = buffer_1010_pppd[242];

    auto g_x_0_y_0_y_y_y_yy = buffer_1010_pppd[243];

    auto g_x_0_y_0_y_y_y_yz = buffer_1010_pppd[244];

    auto g_x_0_y_0_y_y_y_zz = buffer_1010_pppd[245];

    auto g_x_0_y_0_y_y_z_xx = buffer_1010_pppd[246];

    auto g_x_0_y_0_y_y_z_xy = buffer_1010_pppd[247];

    auto g_x_0_y_0_y_y_z_xz = buffer_1010_pppd[248];

    auto g_x_0_y_0_y_y_z_yy = buffer_1010_pppd[249];

    auto g_x_0_y_0_y_y_z_yz = buffer_1010_pppd[250];

    auto g_x_0_y_0_y_y_z_zz = buffer_1010_pppd[251];

    auto g_x_0_y_0_y_z_x_xx = buffer_1010_pppd[252];

    auto g_x_0_y_0_y_z_x_xy = buffer_1010_pppd[253];

    auto g_x_0_y_0_y_z_x_xz = buffer_1010_pppd[254];

    auto g_x_0_y_0_y_z_x_yy = buffer_1010_pppd[255];

    auto g_x_0_y_0_y_z_x_yz = buffer_1010_pppd[256];

    auto g_x_0_y_0_y_z_x_zz = buffer_1010_pppd[257];

    auto g_x_0_y_0_y_z_y_xx = buffer_1010_pppd[258];

    auto g_x_0_y_0_y_z_y_xy = buffer_1010_pppd[259];

    auto g_x_0_y_0_y_z_y_xz = buffer_1010_pppd[260];

    auto g_x_0_y_0_y_z_y_yy = buffer_1010_pppd[261];

    auto g_x_0_y_0_y_z_y_yz = buffer_1010_pppd[262];

    auto g_x_0_y_0_y_z_y_zz = buffer_1010_pppd[263];

    auto g_x_0_y_0_y_z_z_xx = buffer_1010_pppd[264];

    auto g_x_0_y_0_y_z_z_xy = buffer_1010_pppd[265];

    auto g_x_0_y_0_y_z_z_xz = buffer_1010_pppd[266];

    auto g_x_0_y_0_y_z_z_yy = buffer_1010_pppd[267];

    auto g_x_0_y_0_y_z_z_yz = buffer_1010_pppd[268];

    auto g_x_0_y_0_y_z_z_zz = buffer_1010_pppd[269];

    auto g_x_0_y_0_z_x_x_xx = buffer_1010_pppd[270];

    auto g_x_0_y_0_z_x_x_xy = buffer_1010_pppd[271];

    auto g_x_0_y_0_z_x_x_xz = buffer_1010_pppd[272];

    auto g_x_0_y_0_z_x_x_yy = buffer_1010_pppd[273];

    auto g_x_0_y_0_z_x_x_yz = buffer_1010_pppd[274];

    auto g_x_0_y_0_z_x_x_zz = buffer_1010_pppd[275];

    auto g_x_0_y_0_z_x_y_xx = buffer_1010_pppd[276];

    auto g_x_0_y_0_z_x_y_xy = buffer_1010_pppd[277];

    auto g_x_0_y_0_z_x_y_xz = buffer_1010_pppd[278];

    auto g_x_0_y_0_z_x_y_yy = buffer_1010_pppd[279];

    auto g_x_0_y_0_z_x_y_yz = buffer_1010_pppd[280];

    auto g_x_0_y_0_z_x_y_zz = buffer_1010_pppd[281];

    auto g_x_0_y_0_z_x_z_xx = buffer_1010_pppd[282];

    auto g_x_0_y_0_z_x_z_xy = buffer_1010_pppd[283];

    auto g_x_0_y_0_z_x_z_xz = buffer_1010_pppd[284];

    auto g_x_0_y_0_z_x_z_yy = buffer_1010_pppd[285];

    auto g_x_0_y_0_z_x_z_yz = buffer_1010_pppd[286];

    auto g_x_0_y_0_z_x_z_zz = buffer_1010_pppd[287];

    auto g_x_0_y_0_z_y_x_xx = buffer_1010_pppd[288];

    auto g_x_0_y_0_z_y_x_xy = buffer_1010_pppd[289];

    auto g_x_0_y_0_z_y_x_xz = buffer_1010_pppd[290];

    auto g_x_0_y_0_z_y_x_yy = buffer_1010_pppd[291];

    auto g_x_0_y_0_z_y_x_yz = buffer_1010_pppd[292];

    auto g_x_0_y_0_z_y_x_zz = buffer_1010_pppd[293];

    auto g_x_0_y_0_z_y_y_xx = buffer_1010_pppd[294];

    auto g_x_0_y_0_z_y_y_xy = buffer_1010_pppd[295];

    auto g_x_0_y_0_z_y_y_xz = buffer_1010_pppd[296];

    auto g_x_0_y_0_z_y_y_yy = buffer_1010_pppd[297];

    auto g_x_0_y_0_z_y_y_yz = buffer_1010_pppd[298];

    auto g_x_0_y_0_z_y_y_zz = buffer_1010_pppd[299];

    auto g_x_0_y_0_z_y_z_xx = buffer_1010_pppd[300];

    auto g_x_0_y_0_z_y_z_xy = buffer_1010_pppd[301];

    auto g_x_0_y_0_z_y_z_xz = buffer_1010_pppd[302];

    auto g_x_0_y_0_z_y_z_yy = buffer_1010_pppd[303];

    auto g_x_0_y_0_z_y_z_yz = buffer_1010_pppd[304];

    auto g_x_0_y_0_z_y_z_zz = buffer_1010_pppd[305];

    auto g_x_0_y_0_z_z_x_xx = buffer_1010_pppd[306];

    auto g_x_0_y_0_z_z_x_xy = buffer_1010_pppd[307];

    auto g_x_0_y_0_z_z_x_xz = buffer_1010_pppd[308];

    auto g_x_0_y_0_z_z_x_yy = buffer_1010_pppd[309];

    auto g_x_0_y_0_z_z_x_yz = buffer_1010_pppd[310];

    auto g_x_0_y_0_z_z_x_zz = buffer_1010_pppd[311];

    auto g_x_0_y_0_z_z_y_xx = buffer_1010_pppd[312];

    auto g_x_0_y_0_z_z_y_xy = buffer_1010_pppd[313];

    auto g_x_0_y_0_z_z_y_xz = buffer_1010_pppd[314];

    auto g_x_0_y_0_z_z_y_yy = buffer_1010_pppd[315];

    auto g_x_0_y_0_z_z_y_yz = buffer_1010_pppd[316];

    auto g_x_0_y_0_z_z_y_zz = buffer_1010_pppd[317];

    auto g_x_0_y_0_z_z_z_xx = buffer_1010_pppd[318];

    auto g_x_0_y_0_z_z_z_xy = buffer_1010_pppd[319];

    auto g_x_0_y_0_z_z_z_xz = buffer_1010_pppd[320];

    auto g_x_0_y_0_z_z_z_yy = buffer_1010_pppd[321];

    auto g_x_0_y_0_z_z_z_yz = buffer_1010_pppd[322];

    auto g_x_0_y_0_z_z_z_zz = buffer_1010_pppd[323];

    auto g_x_0_z_0_x_x_x_xx = buffer_1010_pppd[324];

    auto g_x_0_z_0_x_x_x_xy = buffer_1010_pppd[325];

    auto g_x_0_z_0_x_x_x_xz = buffer_1010_pppd[326];

    auto g_x_0_z_0_x_x_x_yy = buffer_1010_pppd[327];

    auto g_x_0_z_0_x_x_x_yz = buffer_1010_pppd[328];

    auto g_x_0_z_0_x_x_x_zz = buffer_1010_pppd[329];

    auto g_x_0_z_0_x_x_y_xx = buffer_1010_pppd[330];

    auto g_x_0_z_0_x_x_y_xy = buffer_1010_pppd[331];

    auto g_x_0_z_0_x_x_y_xz = buffer_1010_pppd[332];

    auto g_x_0_z_0_x_x_y_yy = buffer_1010_pppd[333];

    auto g_x_0_z_0_x_x_y_yz = buffer_1010_pppd[334];

    auto g_x_0_z_0_x_x_y_zz = buffer_1010_pppd[335];

    auto g_x_0_z_0_x_x_z_xx = buffer_1010_pppd[336];

    auto g_x_0_z_0_x_x_z_xy = buffer_1010_pppd[337];

    auto g_x_0_z_0_x_x_z_xz = buffer_1010_pppd[338];

    auto g_x_0_z_0_x_x_z_yy = buffer_1010_pppd[339];

    auto g_x_0_z_0_x_x_z_yz = buffer_1010_pppd[340];

    auto g_x_0_z_0_x_x_z_zz = buffer_1010_pppd[341];

    auto g_x_0_z_0_x_y_x_xx = buffer_1010_pppd[342];

    auto g_x_0_z_0_x_y_x_xy = buffer_1010_pppd[343];

    auto g_x_0_z_0_x_y_x_xz = buffer_1010_pppd[344];

    auto g_x_0_z_0_x_y_x_yy = buffer_1010_pppd[345];

    auto g_x_0_z_0_x_y_x_yz = buffer_1010_pppd[346];

    auto g_x_0_z_0_x_y_x_zz = buffer_1010_pppd[347];

    auto g_x_0_z_0_x_y_y_xx = buffer_1010_pppd[348];

    auto g_x_0_z_0_x_y_y_xy = buffer_1010_pppd[349];

    auto g_x_0_z_0_x_y_y_xz = buffer_1010_pppd[350];

    auto g_x_0_z_0_x_y_y_yy = buffer_1010_pppd[351];

    auto g_x_0_z_0_x_y_y_yz = buffer_1010_pppd[352];

    auto g_x_0_z_0_x_y_y_zz = buffer_1010_pppd[353];

    auto g_x_0_z_0_x_y_z_xx = buffer_1010_pppd[354];

    auto g_x_0_z_0_x_y_z_xy = buffer_1010_pppd[355];

    auto g_x_0_z_0_x_y_z_xz = buffer_1010_pppd[356];

    auto g_x_0_z_0_x_y_z_yy = buffer_1010_pppd[357];

    auto g_x_0_z_0_x_y_z_yz = buffer_1010_pppd[358];

    auto g_x_0_z_0_x_y_z_zz = buffer_1010_pppd[359];

    auto g_x_0_z_0_x_z_x_xx = buffer_1010_pppd[360];

    auto g_x_0_z_0_x_z_x_xy = buffer_1010_pppd[361];

    auto g_x_0_z_0_x_z_x_xz = buffer_1010_pppd[362];

    auto g_x_0_z_0_x_z_x_yy = buffer_1010_pppd[363];

    auto g_x_0_z_0_x_z_x_yz = buffer_1010_pppd[364];

    auto g_x_0_z_0_x_z_x_zz = buffer_1010_pppd[365];

    auto g_x_0_z_0_x_z_y_xx = buffer_1010_pppd[366];

    auto g_x_0_z_0_x_z_y_xy = buffer_1010_pppd[367];

    auto g_x_0_z_0_x_z_y_xz = buffer_1010_pppd[368];

    auto g_x_0_z_0_x_z_y_yy = buffer_1010_pppd[369];

    auto g_x_0_z_0_x_z_y_yz = buffer_1010_pppd[370];

    auto g_x_0_z_0_x_z_y_zz = buffer_1010_pppd[371];

    auto g_x_0_z_0_x_z_z_xx = buffer_1010_pppd[372];

    auto g_x_0_z_0_x_z_z_xy = buffer_1010_pppd[373];

    auto g_x_0_z_0_x_z_z_xz = buffer_1010_pppd[374];

    auto g_x_0_z_0_x_z_z_yy = buffer_1010_pppd[375];

    auto g_x_0_z_0_x_z_z_yz = buffer_1010_pppd[376];

    auto g_x_0_z_0_x_z_z_zz = buffer_1010_pppd[377];

    auto g_x_0_z_0_y_x_x_xx = buffer_1010_pppd[378];

    auto g_x_0_z_0_y_x_x_xy = buffer_1010_pppd[379];

    auto g_x_0_z_0_y_x_x_xz = buffer_1010_pppd[380];

    auto g_x_0_z_0_y_x_x_yy = buffer_1010_pppd[381];

    auto g_x_0_z_0_y_x_x_yz = buffer_1010_pppd[382];

    auto g_x_0_z_0_y_x_x_zz = buffer_1010_pppd[383];

    auto g_x_0_z_0_y_x_y_xx = buffer_1010_pppd[384];

    auto g_x_0_z_0_y_x_y_xy = buffer_1010_pppd[385];

    auto g_x_0_z_0_y_x_y_xz = buffer_1010_pppd[386];

    auto g_x_0_z_0_y_x_y_yy = buffer_1010_pppd[387];

    auto g_x_0_z_0_y_x_y_yz = buffer_1010_pppd[388];

    auto g_x_0_z_0_y_x_y_zz = buffer_1010_pppd[389];

    auto g_x_0_z_0_y_x_z_xx = buffer_1010_pppd[390];

    auto g_x_0_z_0_y_x_z_xy = buffer_1010_pppd[391];

    auto g_x_0_z_0_y_x_z_xz = buffer_1010_pppd[392];

    auto g_x_0_z_0_y_x_z_yy = buffer_1010_pppd[393];

    auto g_x_0_z_0_y_x_z_yz = buffer_1010_pppd[394];

    auto g_x_0_z_0_y_x_z_zz = buffer_1010_pppd[395];

    auto g_x_0_z_0_y_y_x_xx = buffer_1010_pppd[396];

    auto g_x_0_z_0_y_y_x_xy = buffer_1010_pppd[397];

    auto g_x_0_z_0_y_y_x_xz = buffer_1010_pppd[398];

    auto g_x_0_z_0_y_y_x_yy = buffer_1010_pppd[399];

    auto g_x_0_z_0_y_y_x_yz = buffer_1010_pppd[400];

    auto g_x_0_z_0_y_y_x_zz = buffer_1010_pppd[401];

    auto g_x_0_z_0_y_y_y_xx = buffer_1010_pppd[402];

    auto g_x_0_z_0_y_y_y_xy = buffer_1010_pppd[403];

    auto g_x_0_z_0_y_y_y_xz = buffer_1010_pppd[404];

    auto g_x_0_z_0_y_y_y_yy = buffer_1010_pppd[405];

    auto g_x_0_z_0_y_y_y_yz = buffer_1010_pppd[406];

    auto g_x_0_z_0_y_y_y_zz = buffer_1010_pppd[407];

    auto g_x_0_z_0_y_y_z_xx = buffer_1010_pppd[408];

    auto g_x_0_z_0_y_y_z_xy = buffer_1010_pppd[409];

    auto g_x_0_z_0_y_y_z_xz = buffer_1010_pppd[410];

    auto g_x_0_z_0_y_y_z_yy = buffer_1010_pppd[411];

    auto g_x_0_z_0_y_y_z_yz = buffer_1010_pppd[412];

    auto g_x_0_z_0_y_y_z_zz = buffer_1010_pppd[413];

    auto g_x_0_z_0_y_z_x_xx = buffer_1010_pppd[414];

    auto g_x_0_z_0_y_z_x_xy = buffer_1010_pppd[415];

    auto g_x_0_z_0_y_z_x_xz = buffer_1010_pppd[416];

    auto g_x_0_z_0_y_z_x_yy = buffer_1010_pppd[417];

    auto g_x_0_z_0_y_z_x_yz = buffer_1010_pppd[418];

    auto g_x_0_z_0_y_z_x_zz = buffer_1010_pppd[419];

    auto g_x_0_z_0_y_z_y_xx = buffer_1010_pppd[420];

    auto g_x_0_z_0_y_z_y_xy = buffer_1010_pppd[421];

    auto g_x_0_z_0_y_z_y_xz = buffer_1010_pppd[422];

    auto g_x_0_z_0_y_z_y_yy = buffer_1010_pppd[423];

    auto g_x_0_z_0_y_z_y_yz = buffer_1010_pppd[424];

    auto g_x_0_z_0_y_z_y_zz = buffer_1010_pppd[425];

    auto g_x_0_z_0_y_z_z_xx = buffer_1010_pppd[426];

    auto g_x_0_z_0_y_z_z_xy = buffer_1010_pppd[427];

    auto g_x_0_z_0_y_z_z_xz = buffer_1010_pppd[428];

    auto g_x_0_z_0_y_z_z_yy = buffer_1010_pppd[429];

    auto g_x_0_z_0_y_z_z_yz = buffer_1010_pppd[430];

    auto g_x_0_z_0_y_z_z_zz = buffer_1010_pppd[431];

    auto g_x_0_z_0_z_x_x_xx = buffer_1010_pppd[432];

    auto g_x_0_z_0_z_x_x_xy = buffer_1010_pppd[433];

    auto g_x_0_z_0_z_x_x_xz = buffer_1010_pppd[434];

    auto g_x_0_z_0_z_x_x_yy = buffer_1010_pppd[435];

    auto g_x_0_z_0_z_x_x_yz = buffer_1010_pppd[436];

    auto g_x_0_z_0_z_x_x_zz = buffer_1010_pppd[437];

    auto g_x_0_z_0_z_x_y_xx = buffer_1010_pppd[438];

    auto g_x_0_z_0_z_x_y_xy = buffer_1010_pppd[439];

    auto g_x_0_z_0_z_x_y_xz = buffer_1010_pppd[440];

    auto g_x_0_z_0_z_x_y_yy = buffer_1010_pppd[441];

    auto g_x_0_z_0_z_x_y_yz = buffer_1010_pppd[442];

    auto g_x_0_z_0_z_x_y_zz = buffer_1010_pppd[443];

    auto g_x_0_z_0_z_x_z_xx = buffer_1010_pppd[444];

    auto g_x_0_z_0_z_x_z_xy = buffer_1010_pppd[445];

    auto g_x_0_z_0_z_x_z_xz = buffer_1010_pppd[446];

    auto g_x_0_z_0_z_x_z_yy = buffer_1010_pppd[447];

    auto g_x_0_z_0_z_x_z_yz = buffer_1010_pppd[448];

    auto g_x_0_z_0_z_x_z_zz = buffer_1010_pppd[449];

    auto g_x_0_z_0_z_y_x_xx = buffer_1010_pppd[450];

    auto g_x_0_z_0_z_y_x_xy = buffer_1010_pppd[451];

    auto g_x_0_z_0_z_y_x_xz = buffer_1010_pppd[452];

    auto g_x_0_z_0_z_y_x_yy = buffer_1010_pppd[453];

    auto g_x_0_z_0_z_y_x_yz = buffer_1010_pppd[454];

    auto g_x_0_z_0_z_y_x_zz = buffer_1010_pppd[455];

    auto g_x_0_z_0_z_y_y_xx = buffer_1010_pppd[456];

    auto g_x_0_z_0_z_y_y_xy = buffer_1010_pppd[457];

    auto g_x_0_z_0_z_y_y_xz = buffer_1010_pppd[458];

    auto g_x_0_z_0_z_y_y_yy = buffer_1010_pppd[459];

    auto g_x_0_z_0_z_y_y_yz = buffer_1010_pppd[460];

    auto g_x_0_z_0_z_y_y_zz = buffer_1010_pppd[461];

    auto g_x_0_z_0_z_y_z_xx = buffer_1010_pppd[462];

    auto g_x_0_z_0_z_y_z_xy = buffer_1010_pppd[463];

    auto g_x_0_z_0_z_y_z_xz = buffer_1010_pppd[464];

    auto g_x_0_z_0_z_y_z_yy = buffer_1010_pppd[465];

    auto g_x_0_z_0_z_y_z_yz = buffer_1010_pppd[466];

    auto g_x_0_z_0_z_y_z_zz = buffer_1010_pppd[467];

    auto g_x_0_z_0_z_z_x_xx = buffer_1010_pppd[468];

    auto g_x_0_z_0_z_z_x_xy = buffer_1010_pppd[469];

    auto g_x_0_z_0_z_z_x_xz = buffer_1010_pppd[470];

    auto g_x_0_z_0_z_z_x_yy = buffer_1010_pppd[471];

    auto g_x_0_z_0_z_z_x_yz = buffer_1010_pppd[472];

    auto g_x_0_z_0_z_z_x_zz = buffer_1010_pppd[473];

    auto g_x_0_z_0_z_z_y_xx = buffer_1010_pppd[474];

    auto g_x_0_z_0_z_z_y_xy = buffer_1010_pppd[475];

    auto g_x_0_z_0_z_z_y_xz = buffer_1010_pppd[476];

    auto g_x_0_z_0_z_z_y_yy = buffer_1010_pppd[477];

    auto g_x_0_z_0_z_z_y_yz = buffer_1010_pppd[478];

    auto g_x_0_z_0_z_z_y_zz = buffer_1010_pppd[479];

    auto g_x_0_z_0_z_z_z_xx = buffer_1010_pppd[480];

    auto g_x_0_z_0_z_z_z_xy = buffer_1010_pppd[481];

    auto g_x_0_z_0_z_z_z_xz = buffer_1010_pppd[482];

    auto g_x_0_z_0_z_z_z_yy = buffer_1010_pppd[483];

    auto g_x_0_z_0_z_z_z_yz = buffer_1010_pppd[484];

    auto g_x_0_z_0_z_z_z_zz = buffer_1010_pppd[485];

    auto g_y_0_x_0_x_x_x_xx = buffer_1010_pppd[486];

    auto g_y_0_x_0_x_x_x_xy = buffer_1010_pppd[487];

    auto g_y_0_x_0_x_x_x_xz = buffer_1010_pppd[488];

    auto g_y_0_x_0_x_x_x_yy = buffer_1010_pppd[489];

    auto g_y_0_x_0_x_x_x_yz = buffer_1010_pppd[490];

    auto g_y_0_x_0_x_x_x_zz = buffer_1010_pppd[491];

    auto g_y_0_x_0_x_x_y_xx = buffer_1010_pppd[492];

    auto g_y_0_x_0_x_x_y_xy = buffer_1010_pppd[493];

    auto g_y_0_x_0_x_x_y_xz = buffer_1010_pppd[494];

    auto g_y_0_x_0_x_x_y_yy = buffer_1010_pppd[495];

    auto g_y_0_x_0_x_x_y_yz = buffer_1010_pppd[496];

    auto g_y_0_x_0_x_x_y_zz = buffer_1010_pppd[497];

    auto g_y_0_x_0_x_x_z_xx = buffer_1010_pppd[498];

    auto g_y_0_x_0_x_x_z_xy = buffer_1010_pppd[499];

    auto g_y_0_x_0_x_x_z_xz = buffer_1010_pppd[500];

    auto g_y_0_x_0_x_x_z_yy = buffer_1010_pppd[501];

    auto g_y_0_x_0_x_x_z_yz = buffer_1010_pppd[502];

    auto g_y_0_x_0_x_x_z_zz = buffer_1010_pppd[503];

    auto g_y_0_x_0_x_y_x_xx = buffer_1010_pppd[504];

    auto g_y_0_x_0_x_y_x_xy = buffer_1010_pppd[505];

    auto g_y_0_x_0_x_y_x_xz = buffer_1010_pppd[506];

    auto g_y_0_x_0_x_y_x_yy = buffer_1010_pppd[507];

    auto g_y_0_x_0_x_y_x_yz = buffer_1010_pppd[508];

    auto g_y_0_x_0_x_y_x_zz = buffer_1010_pppd[509];

    auto g_y_0_x_0_x_y_y_xx = buffer_1010_pppd[510];

    auto g_y_0_x_0_x_y_y_xy = buffer_1010_pppd[511];

    auto g_y_0_x_0_x_y_y_xz = buffer_1010_pppd[512];

    auto g_y_0_x_0_x_y_y_yy = buffer_1010_pppd[513];

    auto g_y_0_x_0_x_y_y_yz = buffer_1010_pppd[514];

    auto g_y_0_x_0_x_y_y_zz = buffer_1010_pppd[515];

    auto g_y_0_x_0_x_y_z_xx = buffer_1010_pppd[516];

    auto g_y_0_x_0_x_y_z_xy = buffer_1010_pppd[517];

    auto g_y_0_x_0_x_y_z_xz = buffer_1010_pppd[518];

    auto g_y_0_x_0_x_y_z_yy = buffer_1010_pppd[519];

    auto g_y_0_x_0_x_y_z_yz = buffer_1010_pppd[520];

    auto g_y_0_x_0_x_y_z_zz = buffer_1010_pppd[521];

    auto g_y_0_x_0_x_z_x_xx = buffer_1010_pppd[522];

    auto g_y_0_x_0_x_z_x_xy = buffer_1010_pppd[523];

    auto g_y_0_x_0_x_z_x_xz = buffer_1010_pppd[524];

    auto g_y_0_x_0_x_z_x_yy = buffer_1010_pppd[525];

    auto g_y_0_x_0_x_z_x_yz = buffer_1010_pppd[526];

    auto g_y_0_x_0_x_z_x_zz = buffer_1010_pppd[527];

    auto g_y_0_x_0_x_z_y_xx = buffer_1010_pppd[528];

    auto g_y_0_x_0_x_z_y_xy = buffer_1010_pppd[529];

    auto g_y_0_x_0_x_z_y_xz = buffer_1010_pppd[530];

    auto g_y_0_x_0_x_z_y_yy = buffer_1010_pppd[531];

    auto g_y_0_x_0_x_z_y_yz = buffer_1010_pppd[532];

    auto g_y_0_x_0_x_z_y_zz = buffer_1010_pppd[533];

    auto g_y_0_x_0_x_z_z_xx = buffer_1010_pppd[534];

    auto g_y_0_x_0_x_z_z_xy = buffer_1010_pppd[535];

    auto g_y_0_x_0_x_z_z_xz = buffer_1010_pppd[536];

    auto g_y_0_x_0_x_z_z_yy = buffer_1010_pppd[537];

    auto g_y_0_x_0_x_z_z_yz = buffer_1010_pppd[538];

    auto g_y_0_x_0_x_z_z_zz = buffer_1010_pppd[539];

    auto g_y_0_x_0_y_x_x_xx = buffer_1010_pppd[540];

    auto g_y_0_x_0_y_x_x_xy = buffer_1010_pppd[541];

    auto g_y_0_x_0_y_x_x_xz = buffer_1010_pppd[542];

    auto g_y_0_x_0_y_x_x_yy = buffer_1010_pppd[543];

    auto g_y_0_x_0_y_x_x_yz = buffer_1010_pppd[544];

    auto g_y_0_x_0_y_x_x_zz = buffer_1010_pppd[545];

    auto g_y_0_x_0_y_x_y_xx = buffer_1010_pppd[546];

    auto g_y_0_x_0_y_x_y_xy = buffer_1010_pppd[547];

    auto g_y_0_x_0_y_x_y_xz = buffer_1010_pppd[548];

    auto g_y_0_x_0_y_x_y_yy = buffer_1010_pppd[549];

    auto g_y_0_x_0_y_x_y_yz = buffer_1010_pppd[550];

    auto g_y_0_x_0_y_x_y_zz = buffer_1010_pppd[551];

    auto g_y_0_x_0_y_x_z_xx = buffer_1010_pppd[552];

    auto g_y_0_x_0_y_x_z_xy = buffer_1010_pppd[553];

    auto g_y_0_x_0_y_x_z_xz = buffer_1010_pppd[554];

    auto g_y_0_x_0_y_x_z_yy = buffer_1010_pppd[555];

    auto g_y_0_x_0_y_x_z_yz = buffer_1010_pppd[556];

    auto g_y_0_x_0_y_x_z_zz = buffer_1010_pppd[557];

    auto g_y_0_x_0_y_y_x_xx = buffer_1010_pppd[558];

    auto g_y_0_x_0_y_y_x_xy = buffer_1010_pppd[559];

    auto g_y_0_x_0_y_y_x_xz = buffer_1010_pppd[560];

    auto g_y_0_x_0_y_y_x_yy = buffer_1010_pppd[561];

    auto g_y_0_x_0_y_y_x_yz = buffer_1010_pppd[562];

    auto g_y_0_x_0_y_y_x_zz = buffer_1010_pppd[563];

    auto g_y_0_x_0_y_y_y_xx = buffer_1010_pppd[564];

    auto g_y_0_x_0_y_y_y_xy = buffer_1010_pppd[565];

    auto g_y_0_x_0_y_y_y_xz = buffer_1010_pppd[566];

    auto g_y_0_x_0_y_y_y_yy = buffer_1010_pppd[567];

    auto g_y_0_x_0_y_y_y_yz = buffer_1010_pppd[568];

    auto g_y_0_x_0_y_y_y_zz = buffer_1010_pppd[569];

    auto g_y_0_x_0_y_y_z_xx = buffer_1010_pppd[570];

    auto g_y_0_x_0_y_y_z_xy = buffer_1010_pppd[571];

    auto g_y_0_x_0_y_y_z_xz = buffer_1010_pppd[572];

    auto g_y_0_x_0_y_y_z_yy = buffer_1010_pppd[573];

    auto g_y_0_x_0_y_y_z_yz = buffer_1010_pppd[574];

    auto g_y_0_x_0_y_y_z_zz = buffer_1010_pppd[575];

    auto g_y_0_x_0_y_z_x_xx = buffer_1010_pppd[576];

    auto g_y_0_x_0_y_z_x_xy = buffer_1010_pppd[577];

    auto g_y_0_x_0_y_z_x_xz = buffer_1010_pppd[578];

    auto g_y_0_x_0_y_z_x_yy = buffer_1010_pppd[579];

    auto g_y_0_x_0_y_z_x_yz = buffer_1010_pppd[580];

    auto g_y_0_x_0_y_z_x_zz = buffer_1010_pppd[581];

    auto g_y_0_x_0_y_z_y_xx = buffer_1010_pppd[582];

    auto g_y_0_x_0_y_z_y_xy = buffer_1010_pppd[583];

    auto g_y_0_x_0_y_z_y_xz = buffer_1010_pppd[584];

    auto g_y_0_x_0_y_z_y_yy = buffer_1010_pppd[585];

    auto g_y_0_x_0_y_z_y_yz = buffer_1010_pppd[586];

    auto g_y_0_x_0_y_z_y_zz = buffer_1010_pppd[587];

    auto g_y_0_x_0_y_z_z_xx = buffer_1010_pppd[588];

    auto g_y_0_x_0_y_z_z_xy = buffer_1010_pppd[589];

    auto g_y_0_x_0_y_z_z_xz = buffer_1010_pppd[590];

    auto g_y_0_x_0_y_z_z_yy = buffer_1010_pppd[591];

    auto g_y_0_x_0_y_z_z_yz = buffer_1010_pppd[592];

    auto g_y_0_x_0_y_z_z_zz = buffer_1010_pppd[593];

    auto g_y_0_x_0_z_x_x_xx = buffer_1010_pppd[594];

    auto g_y_0_x_0_z_x_x_xy = buffer_1010_pppd[595];

    auto g_y_0_x_0_z_x_x_xz = buffer_1010_pppd[596];

    auto g_y_0_x_0_z_x_x_yy = buffer_1010_pppd[597];

    auto g_y_0_x_0_z_x_x_yz = buffer_1010_pppd[598];

    auto g_y_0_x_0_z_x_x_zz = buffer_1010_pppd[599];

    auto g_y_0_x_0_z_x_y_xx = buffer_1010_pppd[600];

    auto g_y_0_x_0_z_x_y_xy = buffer_1010_pppd[601];

    auto g_y_0_x_0_z_x_y_xz = buffer_1010_pppd[602];

    auto g_y_0_x_0_z_x_y_yy = buffer_1010_pppd[603];

    auto g_y_0_x_0_z_x_y_yz = buffer_1010_pppd[604];

    auto g_y_0_x_0_z_x_y_zz = buffer_1010_pppd[605];

    auto g_y_0_x_0_z_x_z_xx = buffer_1010_pppd[606];

    auto g_y_0_x_0_z_x_z_xy = buffer_1010_pppd[607];

    auto g_y_0_x_0_z_x_z_xz = buffer_1010_pppd[608];

    auto g_y_0_x_0_z_x_z_yy = buffer_1010_pppd[609];

    auto g_y_0_x_0_z_x_z_yz = buffer_1010_pppd[610];

    auto g_y_0_x_0_z_x_z_zz = buffer_1010_pppd[611];

    auto g_y_0_x_0_z_y_x_xx = buffer_1010_pppd[612];

    auto g_y_0_x_0_z_y_x_xy = buffer_1010_pppd[613];

    auto g_y_0_x_0_z_y_x_xz = buffer_1010_pppd[614];

    auto g_y_0_x_0_z_y_x_yy = buffer_1010_pppd[615];

    auto g_y_0_x_0_z_y_x_yz = buffer_1010_pppd[616];

    auto g_y_0_x_0_z_y_x_zz = buffer_1010_pppd[617];

    auto g_y_0_x_0_z_y_y_xx = buffer_1010_pppd[618];

    auto g_y_0_x_0_z_y_y_xy = buffer_1010_pppd[619];

    auto g_y_0_x_0_z_y_y_xz = buffer_1010_pppd[620];

    auto g_y_0_x_0_z_y_y_yy = buffer_1010_pppd[621];

    auto g_y_0_x_0_z_y_y_yz = buffer_1010_pppd[622];

    auto g_y_0_x_0_z_y_y_zz = buffer_1010_pppd[623];

    auto g_y_0_x_0_z_y_z_xx = buffer_1010_pppd[624];

    auto g_y_0_x_0_z_y_z_xy = buffer_1010_pppd[625];

    auto g_y_0_x_0_z_y_z_xz = buffer_1010_pppd[626];

    auto g_y_0_x_0_z_y_z_yy = buffer_1010_pppd[627];

    auto g_y_0_x_0_z_y_z_yz = buffer_1010_pppd[628];

    auto g_y_0_x_0_z_y_z_zz = buffer_1010_pppd[629];

    auto g_y_0_x_0_z_z_x_xx = buffer_1010_pppd[630];

    auto g_y_0_x_0_z_z_x_xy = buffer_1010_pppd[631];

    auto g_y_0_x_0_z_z_x_xz = buffer_1010_pppd[632];

    auto g_y_0_x_0_z_z_x_yy = buffer_1010_pppd[633];

    auto g_y_0_x_0_z_z_x_yz = buffer_1010_pppd[634];

    auto g_y_0_x_0_z_z_x_zz = buffer_1010_pppd[635];

    auto g_y_0_x_0_z_z_y_xx = buffer_1010_pppd[636];

    auto g_y_0_x_0_z_z_y_xy = buffer_1010_pppd[637];

    auto g_y_0_x_0_z_z_y_xz = buffer_1010_pppd[638];

    auto g_y_0_x_0_z_z_y_yy = buffer_1010_pppd[639];

    auto g_y_0_x_0_z_z_y_yz = buffer_1010_pppd[640];

    auto g_y_0_x_0_z_z_y_zz = buffer_1010_pppd[641];

    auto g_y_0_x_0_z_z_z_xx = buffer_1010_pppd[642];

    auto g_y_0_x_0_z_z_z_xy = buffer_1010_pppd[643];

    auto g_y_0_x_0_z_z_z_xz = buffer_1010_pppd[644];

    auto g_y_0_x_0_z_z_z_yy = buffer_1010_pppd[645];

    auto g_y_0_x_0_z_z_z_yz = buffer_1010_pppd[646];

    auto g_y_0_x_0_z_z_z_zz = buffer_1010_pppd[647];

    auto g_y_0_y_0_x_x_x_xx = buffer_1010_pppd[648];

    auto g_y_0_y_0_x_x_x_xy = buffer_1010_pppd[649];

    auto g_y_0_y_0_x_x_x_xz = buffer_1010_pppd[650];

    auto g_y_0_y_0_x_x_x_yy = buffer_1010_pppd[651];

    auto g_y_0_y_0_x_x_x_yz = buffer_1010_pppd[652];

    auto g_y_0_y_0_x_x_x_zz = buffer_1010_pppd[653];

    auto g_y_0_y_0_x_x_y_xx = buffer_1010_pppd[654];

    auto g_y_0_y_0_x_x_y_xy = buffer_1010_pppd[655];

    auto g_y_0_y_0_x_x_y_xz = buffer_1010_pppd[656];

    auto g_y_0_y_0_x_x_y_yy = buffer_1010_pppd[657];

    auto g_y_0_y_0_x_x_y_yz = buffer_1010_pppd[658];

    auto g_y_0_y_0_x_x_y_zz = buffer_1010_pppd[659];

    auto g_y_0_y_0_x_x_z_xx = buffer_1010_pppd[660];

    auto g_y_0_y_0_x_x_z_xy = buffer_1010_pppd[661];

    auto g_y_0_y_0_x_x_z_xz = buffer_1010_pppd[662];

    auto g_y_0_y_0_x_x_z_yy = buffer_1010_pppd[663];

    auto g_y_0_y_0_x_x_z_yz = buffer_1010_pppd[664];

    auto g_y_0_y_0_x_x_z_zz = buffer_1010_pppd[665];

    auto g_y_0_y_0_x_y_x_xx = buffer_1010_pppd[666];

    auto g_y_0_y_0_x_y_x_xy = buffer_1010_pppd[667];

    auto g_y_0_y_0_x_y_x_xz = buffer_1010_pppd[668];

    auto g_y_0_y_0_x_y_x_yy = buffer_1010_pppd[669];

    auto g_y_0_y_0_x_y_x_yz = buffer_1010_pppd[670];

    auto g_y_0_y_0_x_y_x_zz = buffer_1010_pppd[671];

    auto g_y_0_y_0_x_y_y_xx = buffer_1010_pppd[672];

    auto g_y_0_y_0_x_y_y_xy = buffer_1010_pppd[673];

    auto g_y_0_y_0_x_y_y_xz = buffer_1010_pppd[674];

    auto g_y_0_y_0_x_y_y_yy = buffer_1010_pppd[675];

    auto g_y_0_y_0_x_y_y_yz = buffer_1010_pppd[676];

    auto g_y_0_y_0_x_y_y_zz = buffer_1010_pppd[677];

    auto g_y_0_y_0_x_y_z_xx = buffer_1010_pppd[678];

    auto g_y_0_y_0_x_y_z_xy = buffer_1010_pppd[679];

    auto g_y_0_y_0_x_y_z_xz = buffer_1010_pppd[680];

    auto g_y_0_y_0_x_y_z_yy = buffer_1010_pppd[681];

    auto g_y_0_y_0_x_y_z_yz = buffer_1010_pppd[682];

    auto g_y_0_y_0_x_y_z_zz = buffer_1010_pppd[683];

    auto g_y_0_y_0_x_z_x_xx = buffer_1010_pppd[684];

    auto g_y_0_y_0_x_z_x_xy = buffer_1010_pppd[685];

    auto g_y_0_y_0_x_z_x_xz = buffer_1010_pppd[686];

    auto g_y_0_y_0_x_z_x_yy = buffer_1010_pppd[687];

    auto g_y_0_y_0_x_z_x_yz = buffer_1010_pppd[688];

    auto g_y_0_y_0_x_z_x_zz = buffer_1010_pppd[689];

    auto g_y_0_y_0_x_z_y_xx = buffer_1010_pppd[690];

    auto g_y_0_y_0_x_z_y_xy = buffer_1010_pppd[691];

    auto g_y_0_y_0_x_z_y_xz = buffer_1010_pppd[692];

    auto g_y_0_y_0_x_z_y_yy = buffer_1010_pppd[693];

    auto g_y_0_y_0_x_z_y_yz = buffer_1010_pppd[694];

    auto g_y_0_y_0_x_z_y_zz = buffer_1010_pppd[695];

    auto g_y_0_y_0_x_z_z_xx = buffer_1010_pppd[696];

    auto g_y_0_y_0_x_z_z_xy = buffer_1010_pppd[697];

    auto g_y_0_y_0_x_z_z_xz = buffer_1010_pppd[698];

    auto g_y_0_y_0_x_z_z_yy = buffer_1010_pppd[699];

    auto g_y_0_y_0_x_z_z_yz = buffer_1010_pppd[700];

    auto g_y_0_y_0_x_z_z_zz = buffer_1010_pppd[701];

    auto g_y_0_y_0_y_x_x_xx = buffer_1010_pppd[702];

    auto g_y_0_y_0_y_x_x_xy = buffer_1010_pppd[703];

    auto g_y_0_y_0_y_x_x_xz = buffer_1010_pppd[704];

    auto g_y_0_y_0_y_x_x_yy = buffer_1010_pppd[705];

    auto g_y_0_y_0_y_x_x_yz = buffer_1010_pppd[706];

    auto g_y_0_y_0_y_x_x_zz = buffer_1010_pppd[707];

    auto g_y_0_y_0_y_x_y_xx = buffer_1010_pppd[708];

    auto g_y_0_y_0_y_x_y_xy = buffer_1010_pppd[709];

    auto g_y_0_y_0_y_x_y_xz = buffer_1010_pppd[710];

    auto g_y_0_y_0_y_x_y_yy = buffer_1010_pppd[711];

    auto g_y_0_y_0_y_x_y_yz = buffer_1010_pppd[712];

    auto g_y_0_y_0_y_x_y_zz = buffer_1010_pppd[713];

    auto g_y_0_y_0_y_x_z_xx = buffer_1010_pppd[714];

    auto g_y_0_y_0_y_x_z_xy = buffer_1010_pppd[715];

    auto g_y_0_y_0_y_x_z_xz = buffer_1010_pppd[716];

    auto g_y_0_y_0_y_x_z_yy = buffer_1010_pppd[717];

    auto g_y_0_y_0_y_x_z_yz = buffer_1010_pppd[718];

    auto g_y_0_y_0_y_x_z_zz = buffer_1010_pppd[719];

    auto g_y_0_y_0_y_y_x_xx = buffer_1010_pppd[720];

    auto g_y_0_y_0_y_y_x_xy = buffer_1010_pppd[721];

    auto g_y_0_y_0_y_y_x_xz = buffer_1010_pppd[722];

    auto g_y_0_y_0_y_y_x_yy = buffer_1010_pppd[723];

    auto g_y_0_y_0_y_y_x_yz = buffer_1010_pppd[724];

    auto g_y_0_y_0_y_y_x_zz = buffer_1010_pppd[725];

    auto g_y_0_y_0_y_y_y_xx = buffer_1010_pppd[726];

    auto g_y_0_y_0_y_y_y_xy = buffer_1010_pppd[727];

    auto g_y_0_y_0_y_y_y_xz = buffer_1010_pppd[728];

    auto g_y_0_y_0_y_y_y_yy = buffer_1010_pppd[729];

    auto g_y_0_y_0_y_y_y_yz = buffer_1010_pppd[730];

    auto g_y_0_y_0_y_y_y_zz = buffer_1010_pppd[731];

    auto g_y_0_y_0_y_y_z_xx = buffer_1010_pppd[732];

    auto g_y_0_y_0_y_y_z_xy = buffer_1010_pppd[733];

    auto g_y_0_y_0_y_y_z_xz = buffer_1010_pppd[734];

    auto g_y_0_y_0_y_y_z_yy = buffer_1010_pppd[735];

    auto g_y_0_y_0_y_y_z_yz = buffer_1010_pppd[736];

    auto g_y_0_y_0_y_y_z_zz = buffer_1010_pppd[737];

    auto g_y_0_y_0_y_z_x_xx = buffer_1010_pppd[738];

    auto g_y_0_y_0_y_z_x_xy = buffer_1010_pppd[739];

    auto g_y_0_y_0_y_z_x_xz = buffer_1010_pppd[740];

    auto g_y_0_y_0_y_z_x_yy = buffer_1010_pppd[741];

    auto g_y_0_y_0_y_z_x_yz = buffer_1010_pppd[742];

    auto g_y_0_y_0_y_z_x_zz = buffer_1010_pppd[743];

    auto g_y_0_y_0_y_z_y_xx = buffer_1010_pppd[744];

    auto g_y_0_y_0_y_z_y_xy = buffer_1010_pppd[745];

    auto g_y_0_y_0_y_z_y_xz = buffer_1010_pppd[746];

    auto g_y_0_y_0_y_z_y_yy = buffer_1010_pppd[747];

    auto g_y_0_y_0_y_z_y_yz = buffer_1010_pppd[748];

    auto g_y_0_y_0_y_z_y_zz = buffer_1010_pppd[749];

    auto g_y_0_y_0_y_z_z_xx = buffer_1010_pppd[750];

    auto g_y_0_y_0_y_z_z_xy = buffer_1010_pppd[751];

    auto g_y_0_y_0_y_z_z_xz = buffer_1010_pppd[752];

    auto g_y_0_y_0_y_z_z_yy = buffer_1010_pppd[753];

    auto g_y_0_y_0_y_z_z_yz = buffer_1010_pppd[754];

    auto g_y_0_y_0_y_z_z_zz = buffer_1010_pppd[755];

    auto g_y_0_y_0_z_x_x_xx = buffer_1010_pppd[756];

    auto g_y_0_y_0_z_x_x_xy = buffer_1010_pppd[757];

    auto g_y_0_y_0_z_x_x_xz = buffer_1010_pppd[758];

    auto g_y_0_y_0_z_x_x_yy = buffer_1010_pppd[759];

    auto g_y_0_y_0_z_x_x_yz = buffer_1010_pppd[760];

    auto g_y_0_y_0_z_x_x_zz = buffer_1010_pppd[761];

    auto g_y_0_y_0_z_x_y_xx = buffer_1010_pppd[762];

    auto g_y_0_y_0_z_x_y_xy = buffer_1010_pppd[763];

    auto g_y_0_y_0_z_x_y_xz = buffer_1010_pppd[764];

    auto g_y_0_y_0_z_x_y_yy = buffer_1010_pppd[765];

    auto g_y_0_y_0_z_x_y_yz = buffer_1010_pppd[766];

    auto g_y_0_y_0_z_x_y_zz = buffer_1010_pppd[767];

    auto g_y_0_y_0_z_x_z_xx = buffer_1010_pppd[768];

    auto g_y_0_y_0_z_x_z_xy = buffer_1010_pppd[769];

    auto g_y_0_y_0_z_x_z_xz = buffer_1010_pppd[770];

    auto g_y_0_y_0_z_x_z_yy = buffer_1010_pppd[771];

    auto g_y_0_y_0_z_x_z_yz = buffer_1010_pppd[772];

    auto g_y_0_y_0_z_x_z_zz = buffer_1010_pppd[773];

    auto g_y_0_y_0_z_y_x_xx = buffer_1010_pppd[774];

    auto g_y_0_y_0_z_y_x_xy = buffer_1010_pppd[775];

    auto g_y_0_y_0_z_y_x_xz = buffer_1010_pppd[776];

    auto g_y_0_y_0_z_y_x_yy = buffer_1010_pppd[777];

    auto g_y_0_y_0_z_y_x_yz = buffer_1010_pppd[778];

    auto g_y_0_y_0_z_y_x_zz = buffer_1010_pppd[779];

    auto g_y_0_y_0_z_y_y_xx = buffer_1010_pppd[780];

    auto g_y_0_y_0_z_y_y_xy = buffer_1010_pppd[781];

    auto g_y_0_y_0_z_y_y_xz = buffer_1010_pppd[782];

    auto g_y_0_y_0_z_y_y_yy = buffer_1010_pppd[783];

    auto g_y_0_y_0_z_y_y_yz = buffer_1010_pppd[784];

    auto g_y_0_y_0_z_y_y_zz = buffer_1010_pppd[785];

    auto g_y_0_y_0_z_y_z_xx = buffer_1010_pppd[786];

    auto g_y_0_y_0_z_y_z_xy = buffer_1010_pppd[787];

    auto g_y_0_y_0_z_y_z_xz = buffer_1010_pppd[788];

    auto g_y_0_y_0_z_y_z_yy = buffer_1010_pppd[789];

    auto g_y_0_y_0_z_y_z_yz = buffer_1010_pppd[790];

    auto g_y_0_y_0_z_y_z_zz = buffer_1010_pppd[791];

    auto g_y_0_y_0_z_z_x_xx = buffer_1010_pppd[792];

    auto g_y_0_y_0_z_z_x_xy = buffer_1010_pppd[793];

    auto g_y_0_y_0_z_z_x_xz = buffer_1010_pppd[794];

    auto g_y_0_y_0_z_z_x_yy = buffer_1010_pppd[795];

    auto g_y_0_y_0_z_z_x_yz = buffer_1010_pppd[796];

    auto g_y_0_y_0_z_z_x_zz = buffer_1010_pppd[797];

    auto g_y_0_y_0_z_z_y_xx = buffer_1010_pppd[798];

    auto g_y_0_y_0_z_z_y_xy = buffer_1010_pppd[799];

    auto g_y_0_y_0_z_z_y_xz = buffer_1010_pppd[800];

    auto g_y_0_y_0_z_z_y_yy = buffer_1010_pppd[801];

    auto g_y_0_y_0_z_z_y_yz = buffer_1010_pppd[802];

    auto g_y_0_y_0_z_z_y_zz = buffer_1010_pppd[803];

    auto g_y_0_y_0_z_z_z_xx = buffer_1010_pppd[804];

    auto g_y_0_y_0_z_z_z_xy = buffer_1010_pppd[805];

    auto g_y_0_y_0_z_z_z_xz = buffer_1010_pppd[806];

    auto g_y_0_y_0_z_z_z_yy = buffer_1010_pppd[807];

    auto g_y_0_y_0_z_z_z_yz = buffer_1010_pppd[808];

    auto g_y_0_y_0_z_z_z_zz = buffer_1010_pppd[809];

    auto g_y_0_z_0_x_x_x_xx = buffer_1010_pppd[810];

    auto g_y_0_z_0_x_x_x_xy = buffer_1010_pppd[811];

    auto g_y_0_z_0_x_x_x_xz = buffer_1010_pppd[812];

    auto g_y_0_z_0_x_x_x_yy = buffer_1010_pppd[813];

    auto g_y_0_z_0_x_x_x_yz = buffer_1010_pppd[814];

    auto g_y_0_z_0_x_x_x_zz = buffer_1010_pppd[815];

    auto g_y_0_z_0_x_x_y_xx = buffer_1010_pppd[816];

    auto g_y_0_z_0_x_x_y_xy = buffer_1010_pppd[817];

    auto g_y_0_z_0_x_x_y_xz = buffer_1010_pppd[818];

    auto g_y_0_z_0_x_x_y_yy = buffer_1010_pppd[819];

    auto g_y_0_z_0_x_x_y_yz = buffer_1010_pppd[820];

    auto g_y_0_z_0_x_x_y_zz = buffer_1010_pppd[821];

    auto g_y_0_z_0_x_x_z_xx = buffer_1010_pppd[822];

    auto g_y_0_z_0_x_x_z_xy = buffer_1010_pppd[823];

    auto g_y_0_z_0_x_x_z_xz = buffer_1010_pppd[824];

    auto g_y_0_z_0_x_x_z_yy = buffer_1010_pppd[825];

    auto g_y_0_z_0_x_x_z_yz = buffer_1010_pppd[826];

    auto g_y_0_z_0_x_x_z_zz = buffer_1010_pppd[827];

    auto g_y_0_z_0_x_y_x_xx = buffer_1010_pppd[828];

    auto g_y_0_z_0_x_y_x_xy = buffer_1010_pppd[829];

    auto g_y_0_z_0_x_y_x_xz = buffer_1010_pppd[830];

    auto g_y_0_z_0_x_y_x_yy = buffer_1010_pppd[831];

    auto g_y_0_z_0_x_y_x_yz = buffer_1010_pppd[832];

    auto g_y_0_z_0_x_y_x_zz = buffer_1010_pppd[833];

    auto g_y_0_z_0_x_y_y_xx = buffer_1010_pppd[834];

    auto g_y_0_z_0_x_y_y_xy = buffer_1010_pppd[835];

    auto g_y_0_z_0_x_y_y_xz = buffer_1010_pppd[836];

    auto g_y_0_z_0_x_y_y_yy = buffer_1010_pppd[837];

    auto g_y_0_z_0_x_y_y_yz = buffer_1010_pppd[838];

    auto g_y_0_z_0_x_y_y_zz = buffer_1010_pppd[839];

    auto g_y_0_z_0_x_y_z_xx = buffer_1010_pppd[840];

    auto g_y_0_z_0_x_y_z_xy = buffer_1010_pppd[841];

    auto g_y_0_z_0_x_y_z_xz = buffer_1010_pppd[842];

    auto g_y_0_z_0_x_y_z_yy = buffer_1010_pppd[843];

    auto g_y_0_z_0_x_y_z_yz = buffer_1010_pppd[844];

    auto g_y_0_z_0_x_y_z_zz = buffer_1010_pppd[845];

    auto g_y_0_z_0_x_z_x_xx = buffer_1010_pppd[846];

    auto g_y_0_z_0_x_z_x_xy = buffer_1010_pppd[847];

    auto g_y_0_z_0_x_z_x_xz = buffer_1010_pppd[848];

    auto g_y_0_z_0_x_z_x_yy = buffer_1010_pppd[849];

    auto g_y_0_z_0_x_z_x_yz = buffer_1010_pppd[850];

    auto g_y_0_z_0_x_z_x_zz = buffer_1010_pppd[851];

    auto g_y_0_z_0_x_z_y_xx = buffer_1010_pppd[852];

    auto g_y_0_z_0_x_z_y_xy = buffer_1010_pppd[853];

    auto g_y_0_z_0_x_z_y_xz = buffer_1010_pppd[854];

    auto g_y_0_z_0_x_z_y_yy = buffer_1010_pppd[855];

    auto g_y_0_z_0_x_z_y_yz = buffer_1010_pppd[856];

    auto g_y_0_z_0_x_z_y_zz = buffer_1010_pppd[857];

    auto g_y_0_z_0_x_z_z_xx = buffer_1010_pppd[858];

    auto g_y_0_z_0_x_z_z_xy = buffer_1010_pppd[859];

    auto g_y_0_z_0_x_z_z_xz = buffer_1010_pppd[860];

    auto g_y_0_z_0_x_z_z_yy = buffer_1010_pppd[861];

    auto g_y_0_z_0_x_z_z_yz = buffer_1010_pppd[862];

    auto g_y_0_z_0_x_z_z_zz = buffer_1010_pppd[863];

    auto g_y_0_z_0_y_x_x_xx = buffer_1010_pppd[864];

    auto g_y_0_z_0_y_x_x_xy = buffer_1010_pppd[865];

    auto g_y_0_z_0_y_x_x_xz = buffer_1010_pppd[866];

    auto g_y_0_z_0_y_x_x_yy = buffer_1010_pppd[867];

    auto g_y_0_z_0_y_x_x_yz = buffer_1010_pppd[868];

    auto g_y_0_z_0_y_x_x_zz = buffer_1010_pppd[869];

    auto g_y_0_z_0_y_x_y_xx = buffer_1010_pppd[870];

    auto g_y_0_z_0_y_x_y_xy = buffer_1010_pppd[871];

    auto g_y_0_z_0_y_x_y_xz = buffer_1010_pppd[872];

    auto g_y_0_z_0_y_x_y_yy = buffer_1010_pppd[873];

    auto g_y_0_z_0_y_x_y_yz = buffer_1010_pppd[874];

    auto g_y_0_z_0_y_x_y_zz = buffer_1010_pppd[875];

    auto g_y_0_z_0_y_x_z_xx = buffer_1010_pppd[876];

    auto g_y_0_z_0_y_x_z_xy = buffer_1010_pppd[877];

    auto g_y_0_z_0_y_x_z_xz = buffer_1010_pppd[878];

    auto g_y_0_z_0_y_x_z_yy = buffer_1010_pppd[879];

    auto g_y_0_z_0_y_x_z_yz = buffer_1010_pppd[880];

    auto g_y_0_z_0_y_x_z_zz = buffer_1010_pppd[881];

    auto g_y_0_z_0_y_y_x_xx = buffer_1010_pppd[882];

    auto g_y_0_z_0_y_y_x_xy = buffer_1010_pppd[883];

    auto g_y_0_z_0_y_y_x_xz = buffer_1010_pppd[884];

    auto g_y_0_z_0_y_y_x_yy = buffer_1010_pppd[885];

    auto g_y_0_z_0_y_y_x_yz = buffer_1010_pppd[886];

    auto g_y_0_z_0_y_y_x_zz = buffer_1010_pppd[887];

    auto g_y_0_z_0_y_y_y_xx = buffer_1010_pppd[888];

    auto g_y_0_z_0_y_y_y_xy = buffer_1010_pppd[889];

    auto g_y_0_z_0_y_y_y_xz = buffer_1010_pppd[890];

    auto g_y_0_z_0_y_y_y_yy = buffer_1010_pppd[891];

    auto g_y_0_z_0_y_y_y_yz = buffer_1010_pppd[892];

    auto g_y_0_z_0_y_y_y_zz = buffer_1010_pppd[893];

    auto g_y_0_z_0_y_y_z_xx = buffer_1010_pppd[894];

    auto g_y_0_z_0_y_y_z_xy = buffer_1010_pppd[895];

    auto g_y_0_z_0_y_y_z_xz = buffer_1010_pppd[896];

    auto g_y_0_z_0_y_y_z_yy = buffer_1010_pppd[897];

    auto g_y_0_z_0_y_y_z_yz = buffer_1010_pppd[898];

    auto g_y_0_z_0_y_y_z_zz = buffer_1010_pppd[899];

    auto g_y_0_z_0_y_z_x_xx = buffer_1010_pppd[900];

    auto g_y_0_z_0_y_z_x_xy = buffer_1010_pppd[901];

    auto g_y_0_z_0_y_z_x_xz = buffer_1010_pppd[902];

    auto g_y_0_z_0_y_z_x_yy = buffer_1010_pppd[903];

    auto g_y_0_z_0_y_z_x_yz = buffer_1010_pppd[904];

    auto g_y_0_z_0_y_z_x_zz = buffer_1010_pppd[905];

    auto g_y_0_z_0_y_z_y_xx = buffer_1010_pppd[906];

    auto g_y_0_z_0_y_z_y_xy = buffer_1010_pppd[907];

    auto g_y_0_z_0_y_z_y_xz = buffer_1010_pppd[908];

    auto g_y_0_z_0_y_z_y_yy = buffer_1010_pppd[909];

    auto g_y_0_z_0_y_z_y_yz = buffer_1010_pppd[910];

    auto g_y_0_z_0_y_z_y_zz = buffer_1010_pppd[911];

    auto g_y_0_z_0_y_z_z_xx = buffer_1010_pppd[912];

    auto g_y_0_z_0_y_z_z_xy = buffer_1010_pppd[913];

    auto g_y_0_z_0_y_z_z_xz = buffer_1010_pppd[914];

    auto g_y_0_z_0_y_z_z_yy = buffer_1010_pppd[915];

    auto g_y_0_z_0_y_z_z_yz = buffer_1010_pppd[916];

    auto g_y_0_z_0_y_z_z_zz = buffer_1010_pppd[917];

    auto g_y_0_z_0_z_x_x_xx = buffer_1010_pppd[918];

    auto g_y_0_z_0_z_x_x_xy = buffer_1010_pppd[919];

    auto g_y_0_z_0_z_x_x_xz = buffer_1010_pppd[920];

    auto g_y_0_z_0_z_x_x_yy = buffer_1010_pppd[921];

    auto g_y_0_z_0_z_x_x_yz = buffer_1010_pppd[922];

    auto g_y_0_z_0_z_x_x_zz = buffer_1010_pppd[923];

    auto g_y_0_z_0_z_x_y_xx = buffer_1010_pppd[924];

    auto g_y_0_z_0_z_x_y_xy = buffer_1010_pppd[925];

    auto g_y_0_z_0_z_x_y_xz = buffer_1010_pppd[926];

    auto g_y_0_z_0_z_x_y_yy = buffer_1010_pppd[927];

    auto g_y_0_z_0_z_x_y_yz = buffer_1010_pppd[928];

    auto g_y_0_z_0_z_x_y_zz = buffer_1010_pppd[929];

    auto g_y_0_z_0_z_x_z_xx = buffer_1010_pppd[930];

    auto g_y_0_z_0_z_x_z_xy = buffer_1010_pppd[931];

    auto g_y_0_z_0_z_x_z_xz = buffer_1010_pppd[932];

    auto g_y_0_z_0_z_x_z_yy = buffer_1010_pppd[933];

    auto g_y_0_z_0_z_x_z_yz = buffer_1010_pppd[934];

    auto g_y_0_z_0_z_x_z_zz = buffer_1010_pppd[935];

    auto g_y_0_z_0_z_y_x_xx = buffer_1010_pppd[936];

    auto g_y_0_z_0_z_y_x_xy = buffer_1010_pppd[937];

    auto g_y_0_z_0_z_y_x_xz = buffer_1010_pppd[938];

    auto g_y_0_z_0_z_y_x_yy = buffer_1010_pppd[939];

    auto g_y_0_z_0_z_y_x_yz = buffer_1010_pppd[940];

    auto g_y_0_z_0_z_y_x_zz = buffer_1010_pppd[941];

    auto g_y_0_z_0_z_y_y_xx = buffer_1010_pppd[942];

    auto g_y_0_z_0_z_y_y_xy = buffer_1010_pppd[943];

    auto g_y_0_z_0_z_y_y_xz = buffer_1010_pppd[944];

    auto g_y_0_z_0_z_y_y_yy = buffer_1010_pppd[945];

    auto g_y_0_z_0_z_y_y_yz = buffer_1010_pppd[946];

    auto g_y_0_z_0_z_y_y_zz = buffer_1010_pppd[947];

    auto g_y_0_z_0_z_y_z_xx = buffer_1010_pppd[948];

    auto g_y_0_z_0_z_y_z_xy = buffer_1010_pppd[949];

    auto g_y_0_z_0_z_y_z_xz = buffer_1010_pppd[950];

    auto g_y_0_z_0_z_y_z_yy = buffer_1010_pppd[951];

    auto g_y_0_z_0_z_y_z_yz = buffer_1010_pppd[952];

    auto g_y_0_z_0_z_y_z_zz = buffer_1010_pppd[953];

    auto g_y_0_z_0_z_z_x_xx = buffer_1010_pppd[954];

    auto g_y_0_z_0_z_z_x_xy = buffer_1010_pppd[955];

    auto g_y_0_z_0_z_z_x_xz = buffer_1010_pppd[956];

    auto g_y_0_z_0_z_z_x_yy = buffer_1010_pppd[957];

    auto g_y_0_z_0_z_z_x_yz = buffer_1010_pppd[958];

    auto g_y_0_z_0_z_z_x_zz = buffer_1010_pppd[959];

    auto g_y_0_z_0_z_z_y_xx = buffer_1010_pppd[960];

    auto g_y_0_z_0_z_z_y_xy = buffer_1010_pppd[961];

    auto g_y_0_z_0_z_z_y_xz = buffer_1010_pppd[962];

    auto g_y_0_z_0_z_z_y_yy = buffer_1010_pppd[963];

    auto g_y_0_z_0_z_z_y_yz = buffer_1010_pppd[964];

    auto g_y_0_z_0_z_z_y_zz = buffer_1010_pppd[965];

    auto g_y_0_z_0_z_z_z_xx = buffer_1010_pppd[966];

    auto g_y_0_z_0_z_z_z_xy = buffer_1010_pppd[967];

    auto g_y_0_z_0_z_z_z_xz = buffer_1010_pppd[968];

    auto g_y_0_z_0_z_z_z_yy = buffer_1010_pppd[969];

    auto g_y_0_z_0_z_z_z_yz = buffer_1010_pppd[970];

    auto g_y_0_z_0_z_z_z_zz = buffer_1010_pppd[971];

    auto g_z_0_x_0_x_x_x_xx = buffer_1010_pppd[972];

    auto g_z_0_x_0_x_x_x_xy = buffer_1010_pppd[973];

    auto g_z_0_x_0_x_x_x_xz = buffer_1010_pppd[974];

    auto g_z_0_x_0_x_x_x_yy = buffer_1010_pppd[975];

    auto g_z_0_x_0_x_x_x_yz = buffer_1010_pppd[976];

    auto g_z_0_x_0_x_x_x_zz = buffer_1010_pppd[977];

    auto g_z_0_x_0_x_x_y_xx = buffer_1010_pppd[978];

    auto g_z_0_x_0_x_x_y_xy = buffer_1010_pppd[979];

    auto g_z_0_x_0_x_x_y_xz = buffer_1010_pppd[980];

    auto g_z_0_x_0_x_x_y_yy = buffer_1010_pppd[981];

    auto g_z_0_x_0_x_x_y_yz = buffer_1010_pppd[982];

    auto g_z_0_x_0_x_x_y_zz = buffer_1010_pppd[983];

    auto g_z_0_x_0_x_x_z_xx = buffer_1010_pppd[984];

    auto g_z_0_x_0_x_x_z_xy = buffer_1010_pppd[985];

    auto g_z_0_x_0_x_x_z_xz = buffer_1010_pppd[986];

    auto g_z_0_x_0_x_x_z_yy = buffer_1010_pppd[987];

    auto g_z_0_x_0_x_x_z_yz = buffer_1010_pppd[988];

    auto g_z_0_x_0_x_x_z_zz = buffer_1010_pppd[989];

    auto g_z_0_x_0_x_y_x_xx = buffer_1010_pppd[990];

    auto g_z_0_x_0_x_y_x_xy = buffer_1010_pppd[991];

    auto g_z_0_x_0_x_y_x_xz = buffer_1010_pppd[992];

    auto g_z_0_x_0_x_y_x_yy = buffer_1010_pppd[993];

    auto g_z_0_x_0_x_y_x_yz = buffer_1010_pppd[994];

    auto g_z_0_x_0_x_y_x_zz = buffer_1010_pppd[995];

    auto g_z_0_x_0_x_y_y_xx = buffer_1010_pppd[996];

    auto g_z_0_x_0_x_y_y_xy = buffer_1010_pppd[997];

    auto g_z_0_x_0_x_y_y_xz = buffer_1010_pppd[998];

    auto g_z_0_x_0_x_y_y_yy = buffer_1010_pppd[999];

    auto g_z_0_x_0_x_y_y_yz = buffer_1010_pppd[1000];

    auto g_z_0_x_0_x_y_y_zz = buffer_1010_pppd[1001];

    auto g_z_0_x_0_x_y_z_xx = buffer_1010_pppd[1002];

    auto g_z_0_x_0_x_y_z_xy = buffer_1010_pppd[1003];

    auto g_z_0_x_0_x_y_z_xz = buffer_1010_pppd[1004];

    auto g_z_0_x_0_x_y_z_yy = buffer_1010_pppd[1005];

    auto g_z_0_x_0_x_y_z_yz = buffer_1010_pppd[1006];

    auto g_z_0_x_0_x_y_z_zz = buffer_1010_pppd[1007];

    auto g_z_0_x_0_x_z_x_xx = buffer_1010_pppd[1008];

    auto g_z_0_x_0_x_z_x_xy = buffer_1010_pppd[1009];

    auto g_z_0_x_0_x_z_x_xz = buffer_1010_pppd[1010];

    auto g_z_0_x_0_x_z_x_yy = buffer_1010_pppd[1011];

    auto g_z_0_x_0_x_z_x_yz = buffer_1010_pppd[1012];

    auto g_z_0_x_0_x_z_x_zz = buffer_1010_pppd[1013];

    auto g_z_0_x_0_x_z_y_xx = buffer_1010_pppd[1014];

    auto g_z_0_x_0_x_z_y_xy = buffer_1010_pppd[1015];

    auto g_z_0_x_0_x_z_y_xz = buffer_1010_pppd[1016];

    auto g_z_0_x_0_x_z_y_yy = buffer_1010_pppd[1017];

    auto g_z_0_x_0_x_z_y_yz = buffer_1010_pppd[1018];

    auto g_z_0_x_0_x_z_y_zz = buffer_1010_pppd[1019];

    auto g_z_0_x_0_x_z_z_xx = buffer_1010_pppd[1020];

    auto g_z_0_x_0_x_z_z_xy = buffer_1010_pppd[1021];

    auto g_z_0_x_0_x_z_z_xz = buffer_1010_pppd[1022];

    auto g_z_0_x_0_x_z_z_yy = buffer_1010_pppd[1023];

    auto g_z_0_x_0_x_z_z_yz = buffer_1010_pppd[1024];

    auto g_z_0_x_0_x_z_z_zz = buffer_1010_pppd[1025];

    auto g_z_0_x_0_y_x_x_xx = buffer_1010_pppd[1026];

    auto g_z_0_x_0_y_x_x_xy = buffer_1010_pppd[1027];

    auto g_z_0_x_0_y_x_x_xz = buffer_1010_pppd[1028];

    auto g_z_0_x_0_y_x_x_yy = buffer_1010_pppd[1029];

    auto g_z_0_x_0_y_x_x_yz = buffer_1010_pppd[1030];

    auto g_z_0_x_0_y_x_x_zz = buffer_1010_pppd[1031];

    auto g_z_0_x_0_y_x_y_xx = buffer_1010_pppd[1032];

    auto g_z_0_x_0_y_x_y_xy = buffer_1010_pppd[1033];

    auto g_z_0_x_0_y_x_y_xz = buffer_1010_pppd[1034];

    auto g_z_0_x_0_y_x_y_yy = buffer_1010_pppd[1035];

    auto g_z_0_x_0_y_x_y_yz = buffer_1010_pppd[1036];

    auto g_z_0_x_0_y_x_y_zz = buffer_1010_pppd[1037];

    auto g_z_0_x_0_y_x_z_xx = buffer_1010_pppd[1038];

    auto g_z_0_x_0_y_x_z_xy = buffer_1010_pppd[1039];

    auto g_z_0_x_0_y_x_z_xz = buffer_1010_pppd[1040];

    auto g_z_0_x_0_y_x_z_yy = buffer_1010_pppd[1041];

    auto g_z_0_x_0_y_x_z_yz = buffer_1010_pppd[1042];

    auto g_z_0_x_0_y_x_z_zz = buffer_1010_pppd[1043];

    auto g_z_0_x_0_y_y_x_xx = buffer_1010_pppd[1044];

    auto g_z_0_x_0_y_y_x_xy = buffer_1010_pppd[1045];

    auto g_z_0_x_0_y_y_x_xz = buffer_1010_pppd[1046];

    auto g_z_0_x_0_y_y_x_yy = buffer_1010_pppd[1047];

    auto g_z_0_x_0_y_y_x_yz = buffer_1010_pppd[1048];

    auto g_z_0_x_0_y_y_x_zz = buffer_1010_pppd[1049];

    auto g_z_0_x_0_y_y_y_xx = buffer_1010_pppd[1050];

    auto g_z_0_x_0_y_y_y_xy = buffer_1010_pppd[1051];

    auto g_z_0_x_0_y_y_y_xz = buffer_1010_pppd[1052];

    auto g_z_0_x_0_y_y_y_yy = buffer_1010_pppd[1053];

    auto g_z_0_x_0_y_y_y_yz = buffer_1010_pppd[1054];

    auto g_z_0_x_0_y_y_y_zz = buffer_1010_pppd[1055];

    auto g_z_0_x_0_y_y_z_xx = buffer_1010_pppd[1056];

    auto g_z_0_x_0_y_y_z_xy = buffer_1010_pppd[1057];

    auto g_z_0_x_0_y_y_z_xz = buffer_1010_pppd[1058];

    auto g_z_0_x_0_y_y_z_yy = buffer_1010_pppd[1059];

    auto g_z_0_x_0_y_y_z_yz = buffer_1010_pppd[1060];

    auto g_z_0_x_0_y_y_z_zz = buffer_1010_pppd[1061];

    auto g_z_0_x_0_y_z_x_xx = buffer_1010_pppd[1062];

    auto g_z_0_x_0_y_z_x_xy = buffer_1010_pppd[1063];

    auto g_z_0_x_0_y_z_x_xz = buffer_1010_pppd[1064];

    auto g_z_0_x_0_y_z_x_yy = buffer_1010_pppd[1065];

    auto g_z_0_x_0_y_z_x_yz = buffer_1010_pppd[1066];

    auto g_z_0_x_0_y_z_x_zz = buffer_1010_pppd[1067];

    auto g_z_0_x_0_y_z_y_xx = buffer_1010_pppd[1068];

    auto g_z_0_x_0_y_z_y_xy = buffer_1010_pppd[1069];

    auto g_z_0_x_0_y_z_y_xz = buffer_1010_pppd[1070];

    auto g_z_0_x_0_y_z_y_yy = buffer_1010_pppd[1071];

    auto g_z_0_x_0_y_z_y_yz = buffer_1010_pppd[1072];

    auto g_z_0_x_0_y_z_y_zz = buffer_1010_pppd[1073];

    auto g_z_0_x_0_y_z_z_xx = buffer_1010_pppd[1074];

    auto g_z_0_x_0_y_z_z_xy = buffer_1010_pppd[1075];

    auto g_z_0_x_0_y_z_z_xz = buffer_1010_pppd[1076];

    auto g_z_0_x_0_y_z_z_yy = buffer_1010_pppd[1077];

    auto g_z_0_x_0_y_z_z_yz = buffer_1010_pppd[1078];

    auto g_z_0_x_0_y_z_z_zz = buffer_1010_pppd[1079];

    auto g_z_0_x_0_z_x_x_xx = buffer_1010_pppd[1080];

    auto g_z_0_x_0_z_x_x_xy = buffer_1010_pppd[1081];

    auto g_z_0_x_0_z_x_x_xz = buffer_1010_pppd[1082];

    auto g_z_0_x_0_z_x_x_yy = buffer_1010_pppd[1083];

    auto g_z_0_x_0_z_x_x_yz = buffer_1010_pppd[1084];

    auto g_z_0_x_0_z_x_x_zz = buffer_1010_pppd[1085];

    auto g_z_0_x_0_z_x_y_xx = buffer_1010_pppd[1086];

    auto g_z_0_x_0_z_x_y_xy = buffer_1010_pppd[1087];

    auto g_z_0_x_0_z_x_y_xz = buffer_1010_pppd[1088];

    auto g_z_0_x_0_z_x_y_yy = buffer_1010_pppd[1089];

    auto g_z_0_x_0_z_x_y_yz = buffer_1010_pppd[1090];

    auto g_z_0_x_0_z_x_y_zz = buffer_1010_pppd[1091];

    auto g_z_0_x_0_z_x_z_xx = buffer_1010_pppd[1092];

    auto g_z_0_x_0_z_x_z_xy = buffer_1010_pppd[1093];

    auto g_z_0_x_0_z_x_z_xz = buffer_1010_pppd[1094];

    auto g_z_0_x_0_z_x_z_yy = buffer_1010_pppd[1095];

    auto g_z_0_x_0_z_x_z_yz = buffer_1010_pppd[1096];

    auto g_z_0_x_0_z_x_z_zz = buffer_1010_pppd[1097];

    auto g_z_0_x_0_z_y_x_xx = buffer_1010_pppd[1098];

    auto g_z_0_x_0_z_y_x_xy = buffer_1010_pppd[1099];

    auto g_z_0_x_0_z_y_x_xz = buffer_1010_pppd[1100];

    auto g_z_0_x_0_z_y_x_yy = buffer_1010_pppd[1101];

    auto g_z_0_x_0_z_y_x_yz = buffer_1010_pppd[1102];

    auto g_z_0_x_0_z_y_x_zz = buffer_1010_pppd[1103];

    auto g_z_0_x_0_z_y_y_xx = buffer_1010_pppd[1104];

    auto g_z_0_x_0_z_y_y_xy = buffer_1010_pppd[1105];

    auto g_z_0_x_0_z_y_y_xz = buffer_1010_pppd[1106];

    auto g_z_0_x_0_z_y_y_yy = buffer_1010_pppd[1107];

    auto g_z_0_x_0_z_y_y_yz = buffer_1010_pppd[1108];

    auto g_z_0_x_0_z_y_y_zz = buffer_1010_pppd[1109];

    auto g_z_0_x_0_z_y_z_xx = buffer_1010_pppd[1110];

    auto g_z_0_x_0_z_y_z_xy = buffer_1010_pppd[1111];

    auto g_z_0_x_0_z_y_z_xz = buffer_1010_pppd[1112];

    auto g_z_0_x_0_z_y_z_yy = buffer_1010_pppd[1113];

    auto g_z_0_x_0_z_y_z_yz = buffer_1010_pppd[1114];

    auto g_z_0_x_0_z_y_z_zz = buffer_1010_pppd[1115];

    auto g_z_0_x_0_z_z_x_xx = buffer_1010_pppd[1116];

    auto g_z_0_x_0_z_z_x_xy = buffer_1010_pppd[1117];

    auto g_z_0_x_0_z_z_x_xz = buffer_1010_pppd[1118];

    auto g_z_0_x_0_z_z_x_yy = buffer_1010_pppd[1119];

    auto g_z_0_x_0_z_z_x_yz = buffer_1010_pppd[1120];

    auto g_z_0_x_0_z_z_x_zz = buffer_1010_pppd[1121];

    auto g_z_0_x_0_z_z_y_xx = buffer_1010_pppd[1122];

    auto g_z_0_x_0_z_z_y_xy = buffer_1010_pppd[1123];

    auto g_z_0_x_0_z_z_y_xz = buffer_1010_pppd[1124];

    auto g_z_0_x_0_z_z_y_yy = buffer_1010_pppd[1125];

    auto g_z_0_x_0_z_z_y_yz = buffer_1010_pppd[1126];

    auto g_z_0_x_0_z_z_y_zz = buffer_1010_pppd[1127];

    auto g_z_0_x_0_z_z_z_xx = buffer_1010_pppd[1128];

    auto g_z_0_x_0_z_z_z_xy = buffer_1010_pppd[1129];

    auto g_z_0_x_0_z_z_z_xz = buffer_1010_pppd[1130];

    auto g_z_0_x_0_z_z_z_yy = buffer_1010_pppd[1131];

    auto g_z_0_x_0_z_z_z_yz = buffer_1010_pppd[1132];

    auto g_z_0_x_0_z_z_z_zz = buffer_1010_pppd[1133];

    auto g_z_0_y_0_x_x_x_xx = buffer_1010_pppd[1134];

    auto g_z_0_y_0_x_x_x_xy = buffer_1010_pppd[1135];

    auto g_z_0_y_0_x_x_x_xz = buffer_1010_pppd[1136];

    auto g_z_0_y_0_x_x_x_yy = buffer_1010_pppd[1137];

    auto g_z_0_y_0_x_x_x_yz = buffer_1010_pppd[1138];

    auto g_z_0_y_0_x_x_x_zz = buffer_1010_pppd[1139];

    auto g_z_0_y_0_x_x_y_xx = buffer_1010_pppd[1140];

    auto g_z_0_y_0_x_x_y_xy = buffer_1010_pppd[1141];

    auto g_z_0_y_0_x_x_y_xz = buffer_1010_pppd[1142];

    auto g_z_0_y_0_x_x_y_yy = buffer_1010_pppd[1143];

    auto g_z_0_y_0_x_x_y_yz = buffer_1010_pppd[1144];

    auto g_z_0_y_0_x_x_y_zz = buffer_1010_pppd[1145];

    auto g_z_0_y_0_x_x_z_xx = buffer_1010_pppd[1146];

    auto g_z_0_y_0_x_x_z_xy = buffer_1010_pppd[1147];

    auto g_z_0_y_0_x_x_z_xz = buffer_1010_pppd[1148];

    auto g_z_0_y_0_x_x_z_yy = buffer_1010_pppd[1149];

    auto g_z_0_y_0_x_x_z_yz = buffer_1010_pppd[1150];

    auto g_z_0_y_0_x_x_z_zz = buffer_1010_pppd[1151];

    auto g_z_0_y_0_x_y_x_xx = buffer_1010_pppd[1152];

    auto g_z_0_y_0_x_y_x_xy = buffer_1010_pppd[1153];

    auto g_z_0_y_0_x_y_x_xz = buffer_1010_pppd[1154];

    auto g_z_0_y_0_x_y_x_yy = buffer_1010_pppd[1155];

    auto g_z_0_y_0_x_y_x_yz = buffer_1010_pppd[1156];

    auto g_z_0_y_0_x_y_x_zz = buffer_1010_pppd[1157];

    auto g_z_0_y_0_x_y_y_xx = buffer_1010_pppd[1158];

    auto g_z_0_y_0_x_y_y_xy = buffer_1010_pppd[1159];

    auto g_z_0_y_0_x_y_y_xz = buffer_1010_pppd[1160];

    auto g_z_0_y_0_x_y_y_yy = buffer_1010_pppd[1161];

    auto g_z_0_y_0_x_y_y_yz = buffer_1010_pppd[1162];

    auto g_z_0_y_0_x_y_y_zz = buffer_1010_pppd[1163];

    auto g_z_0_y_0_x_y_z_xx = buffer_1010_pppd[1164];

    auto g_z_0_y_0_x_y_z_xy = buffer_1010_pppd[1165];

    auto g_z_0_y_0_x_y_z_xz = buffer_1010_pppd[1166];

    auto g_z_0_y_0_x_y_z_yy = buffer_1010_pppd[1167];

    auto g_z_0_y_0_x_y_z_yz = buffer_1010_pppd[1168];

    auto g_z_0_y_0_x_y_z_zz = buffer_1010_pppd[1169];

    auto g_z_0_y_0_x_z_x_xx = buffer_1010_pppd[1170];

    auto g_z_0_y_0_x_z_x_xy = buffer_1010_pppd[1171];

    auto g_z_0_y_0_x_z_x_xz = buffer_1010_pppd[1172];

    auto g_z_0_y_0_x_z_x_yy = buffer_1010_pppd[1173];

    auto g_z_0_y_0_x_z_x_yz = buffer_1010_pppd[1174];

    auto g_z_0_y_0_x_z_x_zz = buffer_1010_pppd[1175];

    auto g_z_0_y_0_x_z_y_xx = buffer_1010_pppd[1176];

    auto g_z_0_y_0_x_z_y_xy = buffer_1010_pppd[1177];

    auto g_z_0_y_0_x_z_y_xz = buffer_1010_pppd[1178];

    auto g_z_0_y_0_x_z_y_yy = buffer_1010_pppd[1179];

    auto g_z_0_y_0_x_z_y_yz = buffer_1010_pppd[1180];

    auto g_z_0_y_0_x_z_y_zz = buffer_1010_pppd[1181];

    auto g_z_0_y_0_x_z_z_xx = buffer_1010_pppd[1182];

    auto g_z_0_y_0_x_z_z_xy = buffer_1010_pppd[1183];

    auto g_z_0_y_0_x_z_z_xz = buffer_1010_pppd[1184];

    auto g_z_0_y_0_x_z_z_yy = buffer_1010_pppd[1185];

    auto g_z_0_y_0_x_z_z_yz = buffer_1010_pppd[1186];

    auto g_z_0_y_0_x_z_z_zz = buffer_1010_pppd[1187];

    auto g_z_0_y_0_y_x_x_xx = buffer_1010_pppd[1188];

    auto g_z_0_y_0_y_x_x_xy = buffer_1010_pppd[1189];

    auto g_z_0_y_0_y_x_x_xz = buffer_1010_pppd[1190];

    auto g_z_0_y_0_y_x_x_yy = buffer_1010_pppd[1191];

    auto g_z_0_y_0_y_x_x_yz = buffer_1010_pppd[1192];

    auto g_z_0_y_0_y_x_x_zz = buffer_1010_pppd[1193];

    auto g_z_0_y_0_y_x_y_xx = buffer_1010_pppd[1194];

    auto g_z_0_y_0_y_x_y_xy = buffer_1010_pppd[1195];

    auto g_z_0_y_0_y_x_y_xz = buffer_1010_pppd[1196];

    auto g_z_0_y_0_y_x_y_yy = buffer_1010_pppd[1197];

    auto g_z_0_y_0_y_x_y_yz = buffer_1010_pppd[1198];

    auto g_z_0_y_0_y_x_y_zz = buffer_1010_pppd[1199];

    auto g_z_0_y_0_y_x_z_xx = buffer_1010_pppd[1200];

    auto g_z_0_y_0_y_x_z_xy = buffer_1010_pppd[1201];

    auto g_z_0_y_0_y_x_z_xz = buffer_1010_pppd[1202];

    auto g_z_0_y_0_y_x_z_yy = buffer_1010_pppd[1203];

    auto g_z_0_y_0_y_x_z_yz = buffer_1010_pppd[1204];

    auto g_z_0_y_0_y_x_z_zz = buffer_1010_pppd[1205];

    auto g_z_0_y_0_y_y_x_xx = buffer_1010_pppd[1206];

    auto g_z_0_y_0_y_y_x_xy = buffer_1010_pppd[1207];

    auto g_z_0_y_0_y_y_x_xz = buffer_1010_pppd[1208];

    auto g_z_0_y_0_y_y_x_yy = buffer_1010_pppd[1209];

    auto g_z_0_y_0_y_y_x_yz = buffer_1010_pppd[1210];

    auto g_z_0_y_0_y_y_x_zz = buffer_1010_pppd[1211];

    auto g_z_0_y_0_y_y_y_xx = buffer_1010_pppd[1212];

    auto g_z_0_y_0_y_y_y_xy = buffer_1010_pppd[1213];

    auto g_z_0_y_0_y_y_y_xz = buffer_1010_pppd[1214];

    auto g_z_0_y_0_y_y_y_yy = buffer_1010_pppd[1215];

    auto g_z_0_y_0_y_y_y_yz = buffer_1010_pppd[1216];

    auto g_z_0_y_0_y_y_y_zz = buffer_1010_pppd[1217];

    auto g_z_0_y_0_y_y_z_xx = buffer_1010_pppd[1218];

    auto g_z_0_y_0_y_y_z_xy = buffer_1010_pppd[1219];

    auto g_z_0_y_0_y_y_z_xz = buffer_1010_pppd[1220];

    auto g_z_0_y_0_y_y_z_yy = buffer_1010_pppd[1221];

    auto g_z_0_y_0_y_y_z_yz = buffer_1010_pppd[1222];

    auto g_z_0_y_0_y_y_z_zz = buffer_1010_pppd[1223];

    auto g_z_0_y_0_y_z_x_xx = buffer_1010_pppd[1224];

    auto g_z_0_y_0_y_z_x_xy = buffer_1010_pppd[1225];

    auto g_z_0_y_0_y_z_x_xz = buffer_1010_pppd[1226];

    auto g_z_0_y_0_y_z_x_yy = buffer_1010_pppd[1227];

    auto g_z_0_y_0_y_z_x_yz = buffer_1010_pppd[1228];

    auto g_z_0_y_0_y_z_x_zz = buffer_1010_pppd[1229];

    auto g_z_0_y_0_y_z_y_xx = buffer_1010_pppd[1230];

    auto g_z_0_y_0_y_z_y_xy = buffer_1010_pppd[1231];

    auto g_z_0_y_0_y_z_y_xz = buffer_1010_pppd[1232];

    auto g_z_0_y_0_y_z_y_yy = buffer_1010_pppd[1233];

    auto g_z_0_y_0_y_z_y_yz = buffer_1010_pppd[1234];

    auto g_z_0_y_0_y_z_y_zz = buffer_1010_pppd[1235];

    auto g_z_0_y_0_y_z_z_xx = buffer_1010_pppd[1236];

    auto g_z_0_y_0_y_z_z_xy = buffer_1010_pppd[1237];

    auto g_z_0_y_0_y_z_z_xz = buffer_1010_pppd[1238];

    auto g_z_0_y_0_y_z_z_yy = buffer_1010_pppd[1239];

    auto g_z_0_y_0_y_z_z_yz = buffer_1010_pppd[1240];

    auto g_z_0_y_0_y_z_z_zz = buffer_1010_pppd[1241];

    auto g_z_0_y_0_z_x_x_xx = buffer_1010_pppd[1242];

    auto g_z_0_y_0_z_x_x_xy = buffer_1010_pppd[1243];

    auto g_z_0_y_0_z_x_x_xz = buffer_1010_pppd[1244];

    auto g_z_0_y_0_z_x_x_yy = buffer_1010_pppd[1245];

    auto g_z_0_y_0_z_x_x_yz = buffer_1010_pppd[1246];

    auto g_z_0_y_0_z_x_x_zz = buffer_1010_pppd[1247];

    auto g_z_0_y_0_z_x_y_xx = buffer_1010_pppd[1248];

    auto g_z_0_y_0_z_x_y_xy = buffer_1010_pppd[1249];

    auto g_z_0_y_0_z_x_y_xz = buffer_1010_pppd[1250];

    auto g_z_0_y_0_z_x_y_yy = buffer_1010_pppd[1251];

    auto g_z_0_y_0_z_x_y_yz = buffer_1010_pppd[1252];

    auto g_z_0_y_0_z_x_y_zz = buffer_1010_pppd[1253];

    auto g_z_0_y_0_z_x_z_xx = buffer_1010_pppd[1254];

    auto g_z_0_y_0_z_x_z_xy = buffer_1010_pppd[1255];

    auto g_z_0_y_0_z_x_z_xz = buffer_1010_pppd[1256];

    auto g_z_0_y_0_z_x_z_yy = buffer_1010_pppd[1257];

    auto g_z_0_y_0_z_x_z_yz = buffer_1010_pppd[1258];

    auto g_z_0_y_0_z_x_z_zz = buffer_1010_pppd[1259];

    auto g_z_0_y_0_z_y_x_xx = buffer_1010_pppd[1260];

    auto g_z_0_y_0_z_y_x_xy = buffer_1010_pppd[1261];

    auto g_z_0_y_0_z_y_x_xz = buffer_1010_pppd[1262];

    auto g_z_0_y_0_z_y_x_yy = buffer_1010_pppd[1263];

    auto g_z_0_y_0_z_y_x_yz = buffer_1010_pppd[1264];

    auto g_z_0_y_0_z_y_x_zz = buffer_1010_pppd[1265];

    auto g_z_0_y_0_z_y_y_xx = buffer_1010_pppd[1266];

    auto g_z_0_y_0_z_y_y_xy = buffer_1010_pppd[1267];

    auto g_z_0_y_0_z_y_y_xz = buffer_1010_pppd[1268];

    auto g_z_0_y_0_z_y_y_yy = buffer_1010_pppd[1269];

    auto g_z_0_y_0_z_y_y_yz = buffer_1010_pppd[1270];

    auto g_z_0_y_0_z_y_y_zz = buffer_1010_pppd[1271];

    auto g_z_0_y_0_z_y_z_xx = buffer_1010_pppd[1272];

    auto g_z_0_y_0_z_y_z_xy = buffer_1010_pppd[1273];

    auto g_z_0_y_0_z_y_z_xz = buffer_1010_pppd[1274];

    auto g_z_0_y_0_z_y_z_yy = buffer_1010_pppd[1275];

    auto g_z_0_y_0_z_y_z_yz = buffer_1010_pppd[1276];

    auto g_z_0_y_0_z_y_z_zz = buffer_1010_pppd[1277];

    auto g_z_0_y_0_z_z_x_xx = buffer_1010_pppd[1278];

    auto g_z_0_y_0_z_z_x_xy = buffer_1010_pppd[1279];

    auto g_z_0_y_0_z_z_x_xz = buffer_1010_pppd[1280];

    auto g_z_0_y_0_z_z_x_yy = buffer_1010_pppd[1281];

    auto g_z_0_y_0_z_z_x_yz = buffer_1010_pppd[1282];

    auto g_z_0_y_0_z_z_x_zz = buffer_1010_pppd[1283];

    auto g_z_0_y_0_z_z_y_xx = buffer_1010_pppd[1284];

    auto g_z_0_y_0_z_z_y_xy = buffer_1010_pppd[1285];

    auto g_z_0_y_0_z_z_y_xz = buffer_1010_pppd[1286];

    auto g_z_0_y_0_z_z_y_yy = buffer_1010_pppd[1287];

    auto g_z_0_y_0_z_z_y_yz = buffer_1010_pppd[1288];

    auto g_z_0_y_0_z_z_y_zz = buffer_1010_pppd[1289];

    auto g_z_0_y_0_z_z_z_xx = buffer_1010_pppd[1290];

    auto g_z_0_y_0_z_z_z_xy = buffer_1010_pppd[1291];

    auto g_z_0_y_0_z_z_z_xz = buffer_1010_pppd[1292];

    auto g_z_0_y_0_z_z_z_yy = buffer_1010_pppd[1293];

    auto g_z_0_y_0_z_z_z_yz = buffer_1010_pppd[1294];

    auto g_z_0_y_0_z_z_z_zz = buffer_1010_pppd[1295];

    auto g_z_0_z_0_x_x_x_xx = buffer_1010_pppd[1296];

    auto g_z_0_z_0_x_x_x_xy = buffer_1010_pppd[1297];

    auto g_z_0_z_0_x_x_x_xz = buffer_1010_pppd[1298];

    auto g_z_0_z_0_x_x_x_yy = buffer_1010_pppd[1299];

    auto g_z_0_z_0_x_x_x_yz = buffer_1010_pppd[1300];

    auto g_z_0_z_0_x_x_x_zz = buffer_1010_pppd[1301];

    auto g_z_0_z_0_x_x_y_xx = buffer_1010_pppd[1302];

    auto g_z_0_z_0_x_x_y_xy = buffer_1010_pppd[1303];

    auto g_z_0_z_0_x_x_y_xz = buffer_1010_pppd[1304];

    auto g_z_0_z_0_x_x_y_yy = buffer_1010_pppd[1305];

    auto g_z_0_z_0_x_x_y_yz = buffer_1010_pppd[1306];

    auto g_z_0_z_0_x_x_y_zz = buffer_1010_pppd[1307];

    auto g_z_0_z_0_x_x_z_xx = buffer_1010_pppd[1308];

    auto g_z_0_z_0_x_x_z_xy = buffer_1010_pppd[1309];

    auto g_z_0_z_0_x_x_z_xz = buffer_1010_pppd[1310];

    auto g_z_0_z_0_x_x_z_yy = buffer_1010_pppd[1311];

    auto g_z_0_z_0_x_x_z_yz = buffer_1010_pppd[1312];

    auto g_z_0_z_0_x_x_z_zz = buffer_1010_pppd[1313];

    auto g_z_0_z_0_x_y_x_xx = buffer_1010_pppd[1314];

    auto g_z_0_z_0_x_y_x_xy = buffer_1010_pppd[1315];

    auto g_z_0_z_0_x_y_x_xz = buffer_1010_pppd[1316];

    auto g_z_0_z_0_x_y_x_yy = buffer_1010_pppd[1317];

    auto g_z_0_z_0_x_y_x_yz = buffer_1010_pppd[1318];

    auto g_z_0_z_0_x_y_x_zz = buffer_1010_pppd[1319];

    auto g_z_0_z_0_x_y_y_xx = buffer_1010_pppd[1320];

    auto g_z_0_z_0_x_y_y_xy = buffer_1010_pppd[1321];

    auto g_z_0_z_0_x_y_y_xz = buffer_1010_pppd[1322];

    auto g_z_0_z_0_x_y_y_yy = buffer_1010_pppd[1323];

    auto g_z_0_z_0_x_y_y_yz = buffer_1010_pppd[1324];

    auto g_z_0_z_0_x_y_y_zz = buffer_1010_pppd[1325];

    auto g_z_0_z_0_x_y_z_xx = buffer_1010_pppd[1326];

    auto g_z_0_z_0_x_y_z_xy = buffer_1010_pppd[1327];

    auto g_z_0_z_0_x_y_z_xz = buffer_1010_pppd[1328];

    auto g_z_0_z_0_x_y_z_yy = buffer_1010_pppd[1329];

    auto g_z_0_z_0_x_y_z_yz = buffer_1010_pppd[1330];

    auto g_z_0_z_0_x_y_z_zz = buffer_1010_pppd[1331];

    auto g_z_0_z_0_x_z_x_xx = buffer_1010_pppd[1332];

    auto g_z_0_z_0_x_z_x_xy = buffer_1010_pppd[1333];

    auto g_z_0_z_0_x_z_x_xz = buffer_1010_pppd[1334];

    auto g_z_0_z_0_x_z_x_yy = buffer_1010_pppd[1335];

    auto g_z_0_z_0_x_z_x_yz = buffer_1010_pppd[1336];

    auto g_z_0_z_0_x_z_x_zz = buffer_1010_pppd[1337];

    auto g_z_0_z_0_x_z_y_xx = buffer_1010_pppd[1338];

    auto g_z_0_z_0_x_z_y_xy = buffer_1010_pppd[1339];

    auto g_z_0_z_0_x_z_y_xz = buffer_1010_pppd[1340];

    auto g_z_0_z_0_x_z_y_yy = buffer_1010_pppd[1341];

    auto g_z_0_z_0_x_z_y_yz = buffer_1010_pppd[1342];

    auto g_z_0_z_0_x_z_y_zz = buffer_1010_pppd[1343];

    auto g_z_0_z_0_x_z_z_xx = buffer_1010_pppd[1344];

    auto g_z_0_z_0_x_z_z_xy = buffer_1010_pppd[1345];

    auto g_z_0_z_0_x_z_z_xz = buffer_1010_pppd[1346];

    auto g_z_0_z_0_x_z_z_yy = buffer_1010_pppd[1347];

    auto g_z_0_z_0_x_z_z_yz = buffer_1010_pppd[1348];

    auto g_z_0_z_0_x_z_z_zz = buffer_1010_pppd[1349];

    auto g_z_0_z_0_y_x_x_xx = buffer_1010_pppd[1350];

    auto g_z_0_z_0_y_x_x_xy = buffer_1010_pppd[1351];

    auto g_z_0_z_0_y_x_x_xz = buffer_1010_pppd[1352];

    auto g_z_0_z_0_y_x_x_yy = buffer_1010_pppd[1353];

    auto g_z_0_z_0_y_x_x_yz = buffer_1010_pppd[1354];

    auto g_z_0_z_0_y_x_x_zz = buffer_1010_pppd[1355];

    auto g_z_0_z_0_y_x_y_xx = buffer_1010_pppd[1356];

    auto g_z_0_z_0_y_x_y_xy = buffer_1010_pppd[1357];

    auto g_z_0_z_0_y_x_y_xz = buffer_1010_pppd[1358];

    auto g_z_0_z_0_y_x_y_yy = buffer_1010_pppd[1359];

    auto g_z_0_z_0_y_x_y_yz = buffer_1010_pppd[1360];

    auto g_z_0_z_0_y_x_y_zz = buffer_1010_pppd[1361];

    auto g_z_0_z_0_y_x_z_xx = buffer_1010_pppd[1362];

    auto g_z_0_z_0_y_x_z_xy = buffer_1010_pppd[1363];

    auto g_z_0_z_0_y_x_z_xz = buffer_1010_pppd[1364];

    auto g_z_0_z_0_y_x_z_yy = buffer_1010_pppd[1365];

    auto g_z_0_z_0_y_x_z_yz = buffer_1010_pppd[1366];

    auto g_z_0_z_0_y_x_z_zz = buffer_1010_pppd[1367];

    auto g_z_0_z_0_y_y_x_xx = buffer_1010_pppd[1368];

    auto g_z_0_z_0_y_y_x_xy = buffer_1010_pppd[1369];

    auto g_z_0_z_0_y_y_x_xz = buffer_1010_pppd[1370];

    auto g_z_0_z_0_y_y_x_yy = buffer_1010_pppd[1371];

    auto g_z_0_z_0_y_y_x_yz = buffer_1010_pppd[1372];

    auto g_z_0_z_0_y_y_x_zz = buffer_1010_pppd[1373];

    auto g_z_0_z_0_y_y_y_xx = buffer_1010_pppd[1374];

    auto g_z_0_z_0_y_y_y_xy = buffer_1010_pppd[1375];

    auto g_z_0_z_0_y_y_y_xz = buffer_1010_pppd[1376];

    auto g_z_0_z_0_y_y_y_yy = buffer_1010_pppd[1377];

    auto g_z_0_z_0_y_y_y_yz = buffer_1010_pppd[1378];

    auto g_z_0_z_0_y_y_y_zz = buffer_1010_pppd[1379];

    auto g_z_0_z_0_y_y_z_xx = buffer_1010_pppd[1380];

    auto g_z_0_z_0_y_y_z_xy = buffer_1010_pppd[1381];

    auto g_z_0_z_0_y_y_z_xz = buffer_1010_pppd[1382];

    auto g_z_0_z_0_y_y_z_yy = buffer_1010_pppd[1383];

    auto g_z_0_z_0_y_y_z_yz = buffer_1010_pppd[1384];

    auto g_z_0_z_0_y_y_z_zz = buffer_1010_pppd[1385];

    auto g_z_0_z_0_y_z_x_xx = buffer_1010_pppd[1386];

    auto g_z_0_z_0_y_z_x_xy = buffer_1010_pppd[1387];

    auto g_z_0_z_0_y_z_x_xz = buffer_1010_pppd[1388];

    auto g_z_0_z_0_y_z_x_yy = buffer_1010_pppd[1389];

    auto g_z_0_z_0_y_z_x_yz = buffer_1010_pppd[1390];

    auto g_z_0_z_0_y_z_x_zz = buffer_1010_pppd[1391];

    auto g_z_0_z_0_y_z_y_xx = buffer_1010_pppd[1392];

    auto g_z_0_z_0_y_z_y_xy = buffer_1010_pppd[1393];

    auto g_z_0_z_0_y_z_y_xz = buffer_1010_pppd[1394];

    auto g_z_0_z_0_y_z_y_yy = buffer_1010_pppd[1395];

    auto g_z_0_z_0_y_z_y_yz = buffer_1010_pppd[1396];

    auto g_z_0_z_0_y_z_y_zz = buffer_1010_pppd[1397];

    auto g_z_0_z_0_y_z_z_xx = buffer_1010_pppd[1398];

    auto g_z_0_z_0_y_z_z_xy = buffer_1010_pppd[1399];

    auto g_z_0_z_0_y_z_z_xz = buffer_1010_pppd[1400];

    auto g_z_0_z_0_y_z_z_yy = buffer_1010_pppd[1401];

    auto g_z_0_z_0_y_z_z_yz = buffer_1010_pppd[1402];

    auto g_z_0_z_0_y_z_z_zz = buffer_1010_pppd[1403];

    auto g_z_0_z_0_z_x_x_xx = buffer_1010_pppd[1404];

    auto g_z_0_z_0_z_x_x_xy = buffer_1010_pppd[1405];

    auto g_z_0_z_0_z_x_x_xz = buffer_1010_pppd[1406];

    auto g_z_0_z_0_z_x_x_yy = buffer_1010_pppd[1407];

    auto g_z_0_z_0_z_x_x_yz = buffer_1010_pppd[1408];

    auto g_z_0_z_0_z_x_x_zz = buffer_1010_pppd[1409];

    auto g_z_0_z_0_z_x_y_xx = buffer_1010_pppd[1410];

    auto g_z_0_z_0_z_x_y_xy = buffer_1010_pppd[1411];

    auto g_z_0_z_0_z_x_y_xz = buffer_1010_pppd[1412];

    auto g_z_0_z_0_z_x_y_yy = buffer_1010_pppd[1413];

    auto g_z_0_z_0_z_x_y_yz = buffer_1010_pppd[1414];

    auto g_z_0_z_0_z_x_y_zz = buffer_1010_pppd[1415];

    auto g_z_0_z_0_z_x_z_xx = buffer_1010_pppd[1416];

    auto g_z_0_z_0_z_x_z_xy = buffer_1010_pppd[1417];

    auto g_z_0_z_0_z_x_z_xz = buffer_1010_pppd[1418];

    auto g_z_0_z_0_z_x_z_yy = buffer_1010_pppd[1419];

    auto g_z_0_z_0_z_x_z_yz = buffer_1010_pppd[1420];

    auto g_z_0_z_0_z_x_z_zz = buffer_1010_pppd[1421];

    auto g_z_0_z_0_z_y_x_xx = buffer_1010_pppd[1422];

    auto g_z_0_z_0_z_y_x_xy = buffer_1010_pppd[1423];

    auto g_z_0_z_0_z_y_x_xz = buffer_1010_pppd[1424];

    auto g_z_0_z_0_z_y_x_yy = buffer_1010_pppd[1425];

    auto g_z_0_z_0_z_y_x_yz = buffer_1010_pppd[1426];

    auto g_z_0_z_0_z_y_x_zz = buffer_1010_pppd[1427];

    auto g_z_0_z_0_z_y_y_xx = buffer_1010_pppd[1428];

    auto g_z_0_z_0_z_y_y_xy = buffer_1010_pppd[1429];

    auto g_z_0_z_0_z_y_y_xz = buffer_1010_pppd[1430];

    auto g_z_0_z_0_z_y_y_yy = buffer_1010_pppd[1431];

    auto g_z_0_z_0_z_y_y_yz = buffer_1010_pppd[1432];

    auto g_z_0_z_0_z_y_y_zz = buffer_1010_pppd[1433];

    auto g_z_0_z_0_z_y_z_xx = buffer_1010_pppd[1434];

    auto g_z_0_z_0_z_y_z_xy = buffer_1010_pppd[1435];

    auto g_z_0_z_0_z_y_z_xz = buffer_1010_pppd[1436];

    auto g_z_0_z_0_z_y_z_yy = buffer_1010_pppd[1437];

    auto g_z_0_z_0_z_y_z_yz = buffer_1010_pppd[1438];

    auto g_z_0_z_0_z_y_z_zz = buffer_1010_pppd[1439];

    auto g_z_0_z_0_z_z_x_xx = buffer_1010_pppd[1440];

    auto g_z_0_z_0_z_z_x_xy = buffer_1010_pppd[1441];

    auto g_z_0_z_0_z_z_x_xz = buffer_1010_pppd[1442];

    auto g_z_0_z_0_z_z_x_yy = buffer_1010_pppd[1443];

    auto g_z_0_z_0_z_z_x_yz = buffer_1010_pppd[1444];

    auto g_z_0_z_0_z_z_x_zz = buffer_1010_pppd[1445];

    auto g_z_0_z_0_z_z_y_xx = buffer_1010_pppd[1446];

    auto g_z_0_z_0_z_z_y_xy = buffer_1010_pppd[1447];

    auto g_z_0_z_0_z_z_y_xz = buffer_1010_pppd[1448];

    auto g_z_0_z_0_z_z_y_yy = buffer_1010_pppd[1449];

    auto g_z_0_z_0_z_z_y_yz = buffer_1010_pppd[1450];

    auto g_z_0_z_0_z_z_y_zz = buffer_1010_pppd[1451];

    auto g_z_0_z_0_z_z_z_xx = buffer_1010_pppd[1452];

    auto g_z_0_z_0_z_z_z_xy = buffer_1010_pppd[1453];

    auto g_z_0_z_0_z_z_z_xz = buffer_1010_pppd[1454];

    auto g_z_0_z_0_z_z_z_yy = buffer_1010_pppd[1455];

    auto g_z_0_z_0_z_z_z_yz = buffer_1010_pppd[1456];

    auto g_z_0_z_0_z_z_z_zz = buffer_1010_pppd[1457];

    // integrals block (0-6)

    #pragma omp simd aligned(g_0_x_0_xx, g_0_x_0_xy, g_0_x_0_xz, g_0_x_0_yy, g_0_x_0_yz, g_0_x_0_zz, g_0_x_xx_xx, g_0_x_xx_xy, g_0_x_xx_xz, g_0_x_xx_yy, g_0_x_xx_yz, g_0_x_xx_zz, g_x_0_x_0_x_x_x_xx, g_x_0_x_0_x_x_x_xy, g_x_0_x_0_x_x_x_xz, g_x_0_x_0_x_x_x_yy, g_x_0_x_0_x_x_x_yz, g_x_0_x_0_x_x_x_zz, g_xx_x_0_xx, g_xx_x_0_xy, g_xx_x_0_xz, g_xx_x_0_yy, g_xx_x_0_yz, g_xx_x_0_zz, g_xx_x_xx_xx, g_xx_x_xx_xy, g_xx_x_xx_xz, g_xx_x_xx_yy, g_xx_x_xx_yz, g_xx_x_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_x_x_xx[i] = g_0_x_0_xx[i] - 2.0 * g_0_x_xx_xx[i] * c_exps[i] - 2.0 * g_xx_x_0_xx[i] * a_exp + 4.0 * g_xx_x_xx_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_x_xy[i] = g_0_x_0_xy[i] - 2.0 * g_0_x_xx_xy[i] * c_exps[i] - 2.0 * g_xx_x_0_xy[i] * a_exp + 4.0 * g_xx_x_xx_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_x_xz[i] = g_0_x_0_xz[i] - 2.0 * g_0_x_xx_xz[i] * c_exps[i] - 2.0 * g_xx_x_0_xz[i] * a_exp + 4.0 * g_xx_x_xx_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_x_yy[i] = g_0_x_0_yy[i] - 2.0 * g_0_x_xx_yy[i] * c_exps[i] - 2.0 * g_xx_x_0_yy[i] * a_exp + 4.0 * g_xx_x_xx_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_x_yz[i] = g_0_x_0_yz[i] - 2.0 * g_0_x_xx_yz[i] * c_exps[i] - 2.0 * g_xx_x_0_yz[i] * a_exp + 4.0 * g_xx_x_xx_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_x_zz[i] = g_0_x_0_zz[i] - 2.0 * g_0_x_xx_zz[i] * c_exps[i] - 2.0 * g_xx_x_0_zz[i] * a_exp + 4.0 * g_xx_x_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_0_x_xy_xx, g_0_x_xy_xy, g_0_x_xy_xz, g_0_x_xy_yy, g_0_x_xy_yz, g_0_x_xy_zz, g_x_0_x_0_x_x_y_xx, g_x_0_x_0_x_x_y_xy, g_x_0_x_0_x_x_y_xz, g_x_0_x_0_x_x_y_yy, g_x_0_x_0_x_x_y_yz, g_x_0_x_0_x_x_y_zz, g_xx_x_xy_xx, g_xx_x_xy_xy, g_xx_x_xy_xz, g_xx_x_xy_yy, g_xx_x_xy_yz, g_xx_x_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_x_y_xx[i] = -2.0 * g_0_x_xy_xx[i] * c_exps[i] + 4.0 * g_xx_x_xy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_y_xy[i] = -2.0 * g_0_x_xy_xy[i] * c_exps[i] + 4.0 * g_xx_x_xy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_y_xz[i] = -2.0 * g_0_x_xy_xz[i] * c_exps[i] + 4.0 * g_xx_x_xy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_y_yy[i] = -2.0 * g_0_x_xy_yy[i] * c_exps[i] + 4.0 * g_xx_x_xy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_y_yz[i] = -2.0 * g_0_x_xy_yz[i] * c_exps[i] + 4.0 * g_xx_x_xy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_y_zz[i] = -2.0 * g_0_x_xy_zz[i] * c_exps[i] + 4.0 * g_xx_x_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_0_x_xz_xx, g_0_x_xz_xy, g_0_x_xz_xz, g_0_x_xz_yy, g_0_x_xz_yz, g_0_x_xz_zz, g_x_0_x_0_x_x_z_xx, g_x_0_x_0_x_x_z_xy, g_x_0_x_0_x_x_z_xz, g_x_0_x_0_x_x_z_yy, g_x_0_x_0_x_x_z_yz, g_x_0_x_0_x_x_z_zz, g_xx_x_xz_xx, g_xx_x_xz_xy, g_xx_x_xz_xz, g_xx_x_xz_yy, g_xx_x_xz_yz, g_xx_x_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_x_z_xx[i] = -2.0 * g_0_x_xz_xx[i] * c_exps[i] + 4.0 * g_xx_x_xz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_z_xy[i] = -2.0 * g_0_x_xz_xy[i] * c_exps[i] + 4.0 * g_xx_x_xz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_z_xz[i] = -2.0 * g_0_x_xz_xz[i] * c_exps[i] + 4.0 * g_xx_x_xz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_z_yy[i] = -2.0 * g_0_x_xz_yy[i] * c_exps[i] + 4.0 * g_xx_x_xz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_z_yz[i] = -2.0 * g_0_x_xz_yz[i] * c_exps[i] + 4.0 * g_xx_x_xz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_z_zz[i] = -2.0 * g_0_x_xz_zz[i] * c_exps[i] + 4.0 * g_xx_x_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_0_y_0_xx, g_0_y_0_xy, g_0_y_0_xz, g_0_y_0_yy, g_0_y_0_yz, g_0_y_0_zz, g_0_y_xx_xx, g_0_y_xx_xy, g_0_y_xx_xz, g_0_y_xx_yy, g_0_y_xx_yz, g_0_y_xx_zz, g_x_0_x_0_x_y_x_xx, g_x_0_x_0_x_y_x_xy, g_x_0_x_0_x_y_x_xz, g_x_0_x_0_x_y_x_yy, g_x_0_x_0_x_y_x_yz, g_x_0_x_0_x_y_x_zz, g_xx_y_0_xx, g_xx_y_0_xy, g_xx_y_0_xz, g_xx_y_0_yy, g_xx_y_0_yz, g_xx_y_0_zz, g_xx_y_xx_xx, g_xx_y_xx_xy, g_xx_y_xx_xz, g_xx_y_xx_yy, g_xx_y_xx_yz, g_xx_y_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_y_x_xx[i] = g_0_y_0_xx[i] - 2.0 * g_0_y_xx_xx[i] * c_exps[i] - 2.0 * g_xx_y_0_xx[i] * a_exp + 4.0 * g_xx_y_xx_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_x_xy[i] = g_0_y_0_xy[i] - 2.0 * g_0_y_xx_xy[i] * c_exps[i] - 2.0 * g_xx_y_0_xy[i] * a_exp + 4.0 * g_xx_y_xx_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_x_xz[i] = g_0_y_0_xz[i] - 2.0 * g_0_y_xx_xz[i] * c_exps[i] - 2.0 * g_xx_y_0_xz[i] * a_exp + 4.0 * g_xx_y_xx_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_x_yy[i] = g_0_y_0_yy[i] - 2.0 * g_0_y_xx_yy[i] * c_exps[i] - 2.0 * g_xx_y_0_yy[i] * a_exp + 4.0 * g_xx_y_xx_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_x_yz[i] = g_0_y_0_yz[i] - 2.0 * g_0_y_xx_yz[i] * c_exps[i] - 2.0 * g_xx_y_0_yz[i] * a_exp + 4.0 * g_xx_y_xx_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_x_zz[i] = g_0_y_0_zz[i] - 2.0 * g_0_y_xx_zz[i] * c_exps[i] - 2.0 * g_xx_y_0_zz[i] * a_exp + 4.0 * g_xx_y_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_0_y_xy_xx, g_0_y_xy_xy, g_0_y_xy_xz, g_0_y_xy_yy, g_0_y_xy_yz, g_0_y_xy_zz, g_x_0_x_0_x_y_y_xx, g_x_0_x_0_x_y_y_xy, g_x_0_x_0_x_y_y_xz, g_x_0_x_0_x_y_y_yy, g_x_0_x_0_x_y_y_yz, g_x_0_x_0_x_y_y_zz, g_xx_y_xy_xx, g_xx_y_xy_xy, g_xx_y_xy_xz, g_xx_y_xy_yy, g_xx_y_xy_yz, g_xx_y_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_y_y_xx[i] = -2.0 * g_0_y_xy_xx[i] * c_exps[i] + 4.0 * g_xx_y_xy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_y_xy[i] = -2.0 * g_0_y_xy_xy[i] * c_exps[i] + 4.0 * g_xx_y_xy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_y_xz[i] = -2.0 * g_0_y_xy_xz[i] * c_exps[i] + 4.0 * g_xx_y_xy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_y_yy[i] = -2.0 * g_0_y_xy_yy[i] * c_exps[i] + 4.0 * g_xx_y_xy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_y_yz[i] = -2.0 * g_0_y_xy_yz[i] * c_exps[i] + 4.0 * g_xx_y_xy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_y_zz[i] = -2.0 * g_0_y_xy_zz[i] * c_exps[i] + 4.0 * g_xx_y_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_0_y_xz_xx, g_0_y_xz_xy, g_0_y_xz_xz, g_0_y_xz_yy, g_0_y_xz_yz, g_0_y_xz_zz, g_x_0_x_0_x_y_z_xx, g_x_0_x_0_x_y_z_xy, g_x_0_x_0_x_y_z_xz, g_x_0_x_0_x_y_z_yy, g_x_0_x_0_x_y_z_yz, g_x_0_x_0_x_y_z_zz, g_xx_y_xz_xx, g_xx_y_xz_xy, g_xx_y_xz_xz, g_xx_y_xz_yy, g_xx_y_xz_yz, g_xx_y_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_y_z_xx[i] = -2.0 * g_0_y_xz_xx[i] * c_exps[i] + 4.0 * g_xx_y_xz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_z_xy[i] = -2.0 * g_0_y_xz_xy[i] * c_exps[i] + 4.0 * g_xx_y_xz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_z_xz[i] = -2.0 * g_0_y_xz_xz[i] * c_exps[i] + 4.0 * g_xx_y_xz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_z_yy[i] = -2.0 * g_0_y_xz_yy[i] * c_exps[i] + 4.0 * g_xx_y_xz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_z_yz[i] = -2.0 * g_0_y_xz_yz[i] * c_exps[i] + 4.0 * g_xx_y_xz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_z_zz[i] = -2.0 * g_0_y_xz_zz[i] * c_exps[i] + 4.0 * g_xx_y_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_0_z_0_xx, g_0_z_0_xy, g_0_z_0_xz, g_0_z_0_yy, g_0_z_0_yz, g_0_z_0_zz, g_0_z_xx_xx, g_0_z_xx_xy, g_0_z_xx_xz, g_0_z_xx_yy, g_0_z_xx_yz, g_0_z_xx_zz, g_x_0_x_0_x_z_x_xx, g_x_0_x_0_x_z_x_xy, g_x_0_x_0_x_z_x_xz, g_x_0_x_0_x_z_x_yy, g_x_0_x_0_x_z_x_yz, g_x_0_x_0_x_z_x_zz, g_xx_z_0_xx, g_xx_z_0_xy, g_xx_z_0_xz, g_xx_z_0_yy, g_xx_z_0_yz, g_xx_z_0_zz, g_xx_z_xx_xx, g_xx_z_xx_xy, g_xx_z_xx_xz, g_xx_z_xx_yy, g_xx_z_xx_yz, g_xx_z_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_z_x_xx[i] = g_0_z_0_xx[i] - 2.0 * g_0_z_xx_xx[i] * c_exps[i] - 2.0 * g_xx_z_0_xx[i] * a_exp + 4.0 * g_xx_z_xx_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_x_xy[i] = g_0_z_0_xy[i] - 2.0 * g_0_z_xx_xy[i] * c_exps[i] - 2.0 * g_xx_z_0_xy[i] * a_exp + 4.0 * g_xx_z_xx_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_x_xz[i] = g_0_z_0_xz[i] - 2.0 * g_0_z_xx_xz[i] * c_exps[i] - 2.0 * g_xx_z_0_xz[i] * a_exp + 4.0 * g_xx_z_xx_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_x_yy[i] = g_0_z_0_yy[i] - 2.0 * g_0_z_xx_yy[i] * c_exps[i] - 2.0 * g_xx_z_0_yy[i] * a_exp + 4.0 * g_xx_z_xx_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_x_yz[i] = g_0_z_0_yz[i] - 2.0 * g_0_z_xx_yz[i] * c_exps[i] - 2.0 * g_xx_z_0_yz[i] * a_exp + 4.0 * g_xx_z_xx_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_x_zz[i] = g_0_z_0_zz[i] - 2.0 * g_0_z_xx_zz[i] * c_exps[i] - 2.0 * g_xx_z_0_zz[i] * a_exp + 4.0 * g_xx_z_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_0_z_xy_xx, g_0_z_xy_xy, g_0_z_xy_xz, g_0_z_xy_yy, g_0_z_xy_yz, g_0_z_xy_zz, g_x_0_x_0_x_z_y_xx, g_x_0_x_0_x_z_y_xy, g_x_0_x_0_x_z_y_xz, g_x_0_x_0_x_z_y_yy, g_x_0_x_0_x_z_y_yz, g_x_0_x_0_x_z_y_zz, g_xx_z_xy_xx, g_xx_z_xy_xy, g_xx_z_xy_xz, g_xx_z_xy_yy, g_xx_z_xy_yz, g_xx_z_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_z_y_xx[i] = -2.0 * g_0_z_xy_xx[i] * c_exps[i] + 4.0 * g_xx_z_xy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_y_xy[i] = -2.0 * g_0_z_xy_xy[i] * c_exps[i] + 4.0 * g_xx_z_xy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_y_xz[i] = -2.0 * g_0_z_xy_xz[i] * c_exps[i] + 4.0 * g_xx_z_xy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_y_yy[i] = -2.0 * g_0_z_xy_yy[i] * c_exps[i] + 4.0 * g_xx_z_xy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_y_yz[i] = -2.0 * g_0_z_xy_yz[i] * c_exps[i] + 4.0 * g_xx_z_xy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_y_zz[i] = -2.0 * g_0_z_xy_zz[i] * c_exps[i] + 4.0 * g_xx_z_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_0_z_xz_xx, g_0_z_xz_xy, g_0_z_xz_xz, g_0_z_xz_yy, g_0_z_xz_yz, g_0_z_xz_zz, g_x_0_x_0_x_z_z_xx, g_x_0_x_0_x_z_z_xy, g_x_0_x_0_x_z_z_xz, g_x_0_x_0_x_z_z_yy, g_x_0_x_0_x_z_z_yz, g_x_0_x_0_x_z_z_zz, g_xx_z_xz_xx, g_xx_z_xz_xy, g_xx_z_xz_xz, g_xx_z_xz_yy, g_xx_z_xz_yz, g_xx_z_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_z_z_xx[i] = -2.0 * g_0_z_xz_xx[i] * c_exps[i] + 4.0 * g_xx_z_xz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_z_xy[i] = -2.0 * g_0_z_xz_xy[i] * c_exps[i] + 4.0 * g_xx_z_xz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_z_xz[i] = -2.0 * g_0_z_xz_xz[i] * c_exps[i] + 4.0 * g_xx_z_xz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_z_yy[i] = -2.0 * g_0_z_xz_yy[i] * c_exps[i] + 4.0 * g_xx_z_xz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_z_yz[i] = -2.0 * g_0_z_xz_yz[i] * c_exps[i] + 4.0 * g_xx_z_xz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_z_zz[i] = -2.0 * g_0_z_xz_zz[i] * c_exps[i] + 4.0 * g_xx_z_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_x_0_x_0_y_x_x_xx, g_x_0_x_0_y_x_x_xy, g_x_0_x_0_y_x_x_xz, g_x_0_x_0_y_x_x_yy, g_x_0_x_0_y_x_x_yz, g_x_0_x_0_y_x_x_zz, g_xy_x_0_xx, g_xy_x_0_xy, g_xy_x_0_xz, g_xy_x_0_yy, g_xy_x_0_yz, g_xy_x_0_zz, g_xy_x_xx_xx, g_xy_x_xx_xy, g_xy_x_xx_xz, g_xy_x_xx_yy, g_xy_x_xx_yz, g_xy_x_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_x_x_xx[i] = -2.0 * g_xy_x_0_xx[i] * a_exp + 4.0 * g_xy_x_xx_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_x_xy[i] = -2.0 * g_xy_x_0_xy[i] * a_exp + 4.0 * g_xy_x_xx_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_x_xz[i] = -2.0 * g_xy_x_0_xz[i] * a_exp + 4.0 * g_xy_x_xx_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_x_yy[i] = -2.0 * g_xy_x_0_yy[i] * a_exp + 4.0 * g_xy_x_xx_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_x_yz[i] = -2.0 * g_xy_x_0_yz[i] * a_exp + 4.0 * g_xy_x_xx_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_x_zz[i] = -2.0 * g_xy_x_0_zz[i] * a_exp + 4.0 * g_xy_x_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_x_0_x_0_y_x_y_xx, g_x_0_x_0_y_x_y_xy, g_x_0_x_0_y_x_y_xz, g_x_0_x_0_y_x_y_yy, g_x_0_x_0_y_x_y_yz, g_x_0_x_0_y_x_y_zz, g_xy_x_xy_xx, g_xy_x_xy_xy, g_xy_x_xy_xz, g_xy_x_xy_yy, g_xy_x_xy_yz, g_xy_x_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_x_y_xx[i] = 4.0 * g_xy_x_xy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_y_xy[i] = 4.0 * g_xy_x_xy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_y_xz[i] = 4.0 * g_xy_x_xy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_y_yy[i] = 4.0 * g_xy_x_xy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_y_yz[i] = 4.0 * g_xy_x_xy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_y_zz[i] = 4.0 * g_xy_x_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_x_0_x_0_y_x_z_xx, g_x_0_x_0_y_x_z_xy, g_x_0_x_0_y_x_z_xz, g_x_0_x_0_y_x_z_yy, g_x_0_x_0_y_x_z_yz, g_x_0_x_0_y_x_z_zz, g_xy_x_xz_xx, g_xy_x_xz_xy, g_xy_x_xz_xz, g_xy_x_xz_yy, g_xy_x_xz_yz, g_xy_x_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_x_z_xx[i] = 4.0 * g_xy_x_xz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_z_xy[i] = 4.0 * g_xy_x_xz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_z_xz[i] = 4.0 * g_xy_x_xz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_z_yy[i] = 4.0 * g_xy_x_xz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_z_yz[i] = 4.0 * g_xy_x_xz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_z_zz[i] = 4.0 * g_xy_x_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_x_0_x_0_y_y_x_xx, g_x_0_x_0_y_y_x_xy, g_x_0_x_0_y_y_x_xz, g_x_0_x_0_y_y_x_yy, g_x_0_x_0_y_y_x_yz, g_x_0_x_0_y_y_x_zz, g_xy_y_0_xx, g_xy_y_0_xy, g_xy_y_0_xz, g_xy_y_0_yy, g_xy_y_0_yz, g_xy_y_0_zz, g_xy_y_xx_xx, g_xy_y_xx_xy, g_xy_y_xx_xz, g_xy_y_xx_yy, g_xy_y_xx_yz, g_xy_y_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_y_x_xx[i] = -2.0 * g_xy_y_0_xx[i] * a_exp + 4.0 * g_xy_y_xx_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_x_xy[i] = -2.0 * g_xy_y_0_xy[i] * a_exp + 4.0 * g_xy_y_xx_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_x_xz[i] = -2.0 * g_xy_y_0_xz[i] * a_exp + 4.0 * g_xy_y_xx_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_x_yy[i] = -2.0 * g_xy_y_0_yy[i] * a_exp + 4.0 * g_xy_y_xx_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_x_yz[i] = -2.0 * g_xy_y_0_yz[i] * a_exp + 4.0 * g_xy_y_xx_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_x_zz[i] = -2.0 * g_xy_y_0_zz[i] * a_exp + 4.0 * g_xy_y_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_x_0_x_0_y_y_y_xx, g_x_0_x_0_y_y_y_xy, g_x_0_x_0_y_y_y_xz, g_x_0_x_0_y_y_y_yy, g_x_0_x_0_y_y_y_yz, g_x_0_x_0_y_y_y_zz, g_xy_y_xy_xx, g_xy_y_xy_xy, g_xy_y_xy_xz, g_xy_y_xy_yy, g_xy_y_xy_yz, g_xy_y_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_y_y_xx[i] = 4.0 * g_xy_y_xy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_y_xy[i] = 4.0 * g_xy_y_xy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_y_xz[i] = 4.0 * g_xy_y_xy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_y_yy[i] = 4.0 * g_xy_y_xy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_y_yz[i] = 4.0 * g_xy_y_xy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_y_zz[i] = 4.0 * g_xy_y_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_x_0_x_0_y_y_z_xx, g_x_0_x_0_y_y_z_xy, g_x_0_x_0_y_y_z_xz, g_x_0_x_0_y_y_z_yy, g_x_0_x_0_y_y_z_yz, g_x_0_x_0_y_y_z_zz, g_xy_y_xz_xx, g_xy_y_xz_xy, g_xy_y_xz_xz, g_xy_y_xz_yy, g_xy_y_xz_yz, g_xy_y_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_y_z_xx[i] = 4.0 * g_xy_y_xz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_z_xy[i] = 4.0 * g_xy_y_xz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_z_xz[i] = 4.0 * g_xy_y_xz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_z_yy[i] = 4.0 * g_xy_y_xz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_z_yz[i] = 4.0 * g_xy_y_xz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_z_zz[i] = 4.0 * g_xy_y_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_x_0_x_0_y_z_x_xx, g_x_0_x_0_y_z_x_xy, g_x_0_x_0_y_z_x_xz, g_x_0_x_0_y_z_x_yy, g_x_0_x_0_y_z_x_yz, g_x_0_x_0_y_z_x_zz, g_xy_z_0_xx, g_xy_z_0_xy, g_xy_z_0_xz, g_xy_z_0_yy, g_xy_z_0_yz, g_xy_z_0_zz, g_xy_z_xx_xx, g_xy_z_xx_xy, g_xy_z_xx_xz, g_xy_z_xx_yy, g_xy_z_xx_yz, g_xy_z_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_z_x_xx[i] = -2.0 * g_xy_z_0_xx[i] * a_exp + 4.0 * g_xy_z_xx_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_x_xy[i] = -2.0 * g_xy_z_0_xy[i] * a_exp + 4.0 * g_xy_z_xx_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_x_xz[i] = -2.0 * g_xy_z_0_xz[i] * a_exp + 4.0 * g_xy_z_xx_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_x_yy[i] = -2.0 * g_xy_z_0_yy[i] * a_exp + 4.0 * g_xy_z_xx_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_x_yz[i] = -2.0 * g_xy_z_0_yz[i] * a_exp + 4.0 * g_xy_z_xx_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_x_zz[i] = -2.0 * g_xy_z_0_zz[i] * a_exp + 4.0 * g_xy_z_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_x_0_x_0_y_z_y_xx, g_x_0_x_0_y_z_y_xy, g_x_0_x_0_y_z_y_xz, g_x_0_x_0_y_z_y_yy, g_x_0_x_0_y_z_y_yz, g_x_0_x_0_y_z_y_zz, g_xy_z_xy_xx, g_xy_z_xy_xy, g_xy_z_xy_xz, g_xy_z_xy_yy, g_xy_z_xy_yz, g_xy_z_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_z_y_xx[i] = 4.0 * g_xy_z_xy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_y_xy[i] = 4.0 * g_xy_z_xy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_y_xz[i] = 4.0 * g_xy_z_xy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_y_yy[i] = 4.0 * g_xy_z_xy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_y_yz[i] = 4.0 * g_xy_z_xy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_y_zz[i] = 4.0 * g_xy_z_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_x_0_x_0_y_z_z_xx, g_x_0_x_0_y_z_z_xy, g_x_0_x_0_y_z_z_xz, g_x_0_x_0_y_z_z_yy, g_x_0_x_0_y_z_z_yz, g_x_0_x_0_y_z_z_zz, g_xy_z_xz_xx, g_xy_z_xz_xy, g_xy_z_xz_xz, g_xy_z_xz_yy, g_xy_z_xz_yz, g_xy_z_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_z_z_xx[i] = 4.0 * g_xy_z_xz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_z_xy[i] = 4.0 * g_xy_z_xz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_z_xz[i] = 4.0 * g_xy_z_xz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_z_yy[i] = 4.0 * g_xy_z_xz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_z_yz[i] = 4.0 * g_xy_z_xz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_z_zz[i] = 4.0 * g_xy_z_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_x_0_x_0_z_x_x_xx, g_x_0_x_0_z_x_x_xy, g_x_0_x_0_z_x_x_xz, g_x_0_x_0_z_x_x_yy, g_x_0_x_0_z_x_x_yz, g_x_0_x_0_z_x_x_zz, g_xz_x_0_xx, g_xz_x_0_xy, g_xz_x_0_xz, g_xz_x_0_yy, g_xz_x_0_yz, g_xz_x_0_zz, g_xz_x_xx_xx, g_xz_x_xx_xy, g_xz_x_xx_xz, g_xz_x_xx_yy, g_xz_x_xx_yz, g_xz_x_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_x_x_xx[i] = -2.0 * g_xz_x_0_xx[i] * a_exp + 4.0 * g_xz_x_xx_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_x_xy[i] = -2.0 * g_xz_x_0_xy[i] * a_exp + 4.0 * g_xz_x_xx_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_x_xz[i] = -2.0 * g_xz_x_0_xz[i] * a_exp + 4.0 * g_xz_x_xx_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_x_yy[i] = -2.0 * g_xz_x_0_yy[i] * a_exp + 4.0 * g_xz_x_xx_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_x_yz[i] = -2.0 * g_xz_x_0_yz[i] * a_exp + 4.0 * g_xz_x_xx_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_x_zz[i] = -2.0 * g_xz_x_0_zz[i] * a_exp + 4.0 * g_xz_x_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_x_0_x_0_z_x_y_xx, g_x_0_x_0_z_x_y_xy, g_x_0_x_0_z_x_y_xz, g_x_0_x_0_z_x_y_yy, g_x_0_x_0_z_x_y_yz, g_x_0_x_0_z_x_y_zz, g_xz_x_xy_xx, g_xz_x_xy_xy, g_xz_x_xy_xz, g_xz_x_xy_yy, g_xz_x_xy_yz, g_xz_x_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_x_y_xx[i] = 4.0 * g_xz_x_xy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_y_xy[i] = 4.0 * g_xz_x_xy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_y_xz[i] = 4.0 * g_xz_x_xy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_y_yy[i] = 4.0 * g_xz_x_xy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_y_yz[i] = 4.0 * g_xz_x_xy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_y_zz[i] = 4.0 * g_xz_x_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_x_0_x_0_z_x_z_xx, g_x_0_x_0_z_x_z_xy, g_x_0_x_0_z_x_z_xz, g_x_0_x_0_z_x_z_yy, g_x_0_x_0_z_x_z_yz, g_x_0_x_0_z_x_z_zz, g_xz_x_xz_xx, g_xz_x_xz_xy, g_xz_x_xz_xz, g_xz_x_xz_yy, g_xz_x_xz_yz, g_xz_x_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_x_z_xx[i] = 4.0 * g_xz_x_xz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_z_xy[i] = 4.0 * g_xz_x_xz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_z_xz[i] = 4.0 * g_xz_x_xz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_z_yy[i] = 4.0 * g_xz_x_xz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_z_yz[i] = 4.0 * g_xz_x_xz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_z_zz[i] = 4.0 * g_xz_x_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_x_0_x_0_z_y_x_xx, g_x_0_x_0_z_y_x_xy, g_x_0_x_0_z_y_x_xz, g_x_0_x_0_z_y_x_yy, g_x_0_x_0_z_y_x_yz, g_x_0_x_0_z_y_x_zz, g_xz_y_0_xx, g_xz_y_0_xy, g_xz_y_0_xz, g_xz_y_0_yy, g_xz_y_0_yz, g_xz_y_0_zz, g_xz_y_xx_xx, g_xz_y_xx_xy, g_xz_y_xx_xz, g_xz_y_xx_yy, g_xz_y_xx_yz, g_xz_y_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_y_x_xx[i] = -2.0 * g_xz_y_0_xx[i] * a_exp + 4.0 * g_xz_y_xx_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_x_xy[i] = -2.0 * g_xz_y_0_xy[i] * a_exp + 4.0 * g_xz_y_xx_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_x_xz[i] = -2.0 * g_xz_y_0_xz[i] * a_exp + 4.0 * g_xz_y_xx_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_x_yy[i] = -2.0 * g_xz_y_0_yy[i] * a_exp + 4.0 * g_xz_y_xx_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_x_yz[i] = -2.0 * g_xz_y_0_yz[i] * a_exp + 4.0 * g_xz_y_xx_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_x_zz[i] = -2.0 * g_xz_y_0_zz[i] * a_exp + 4.0 * g_xz_y_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_x_0_x_0_z_y_y_xx, g_x_0_x_0_z_y_y_xy, g_x_0_x_0_z_y_y_xz, g_x_0_x_0_z_y_y_yy, g_x_0_x_0_z_y_y_yz, g_x_0_x_0_z_y_y_zz, g_xz_y_xy_xx, g_xz_y_xy_xy, g_xz_y_xy_xz, g_xz_y_xy_yy, g_xz_y_xy_yz, g_xz_y_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_y_y_xx[i] = 4.0 * g_xz_y_xy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_y_xy[i] = 4.0 * g_xz_y_xy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_y_xz[i] = 4.0 * g_xz_y_xy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_y_yy[i] = 4.0 * g_xz_y_xy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_y_yz[i] = 4.0 * g_xz_y_xy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_y_zz[i] = 4.0 * g_xz_y_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_x_0_x_0_z_y_z_xx, g_x_0_x_0_z_y_z_xy, g_x_0_x_0_z_y_z_xz, g_x_0_x_0_z_y_z_yy, g_x_0_x_0_z_y_z_yz, g_x_0_x_0_z_y_z_zz, g_xz_y_xz_xx, g_xz_y_xz_xy, g_xz_y_xz_xz, g_xz_y_xz_yy, g_xz_y_xz_yz, g_xz_y_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_y_z_xx[i] = 4.0 * g_xz_y_xz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_z_xy[i] = 4.0 * g_xz_y_xz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_z_xz[i] = 4.0 * g_xz_y_xz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_z_yy[i] = 4.0 * g_xz_y_xz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_z_yz[i] = 4.0 * g_xz_y_xz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_z_zz[i] = 4.0 * g_xz_y_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_x_0_x_0_z_z_x_xx, g_x_0_x_0_z_z_x_xy, g_x_0_x_0_z_z_x_xz, g_x_0_x_0_z_z_x_yy, g_x_0_x_0_z_z_x_yz, g_x_0_x_0_z_z_x_zz, g_xz_z_0_xx, g_xz_z_0_xy, g_xz_z_0_xz, g_xz_z_0_yy, g_xz_z_0_yz, g_xz_z_0_zz, g_xz_z_xx_xx, g_xz_z_xx_xy, g_xz_z_xx_xz, g_xz_z_xx_yy, g_xz_z_xx_yz, g_xz_z_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_z_x_xx[i] = -2.0 * g_xz_z_0_xx[i] * a_exp + 4.0 * g_xz_z_xx_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_x_xy[i] = -2.0 * g_xz_z_0_xy[i] * a_exp + 4.0 * g_xz_z_xx_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_x_xz[i] = -2.0 * g_xz_z_0_xz[i] * a_exp + 4.0 * g_xz_z_xx_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_x_yy[i] = -2.0 * g_xz_z_0_yy[i] * a_exp + 4.0 * g_xz_z_xx_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_x_yz[i] = -2.0 * g_xz_z_0_yz[i] * a_exp + 4.0 * g_xz_z_xx_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_x_zz[i] = -2.0 * g_xz_z_0_zz[i] * a_exp + 4.0 * g_xz_z_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_x_0_x_0_z_z_y_xx, g_x_0_x_0_z_z_y_xy, g_x_0_x_0_z_z_y_xz, g_x_0_x_0_z_z_y_yy, g_x_0_x_0_z_z_y_yz, g_x_0_x_0_z_z_y_zz, g_xz_z_xy_xx, g_xz_z_xy_xy, g_xz_z_xy_xz, g_xz_z_xy_yy, g_xz_z_xy_yz, g_xz_z_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_z_y_xx[i] = 4.0 * g_xz_z_xy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_y_xy[i] = 4.0 * g_xz_z_xy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_y_xz[i] = 4.0 * g_xz_z_xy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_y_yy[i] = 4.0 * g_xz_z_xy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_y_yz[i] = 4.0 * g_xz_z_xy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_y_zz[i] = 4.0 * g_xz_z_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_x_0_x_0_z_z_z_xx, g_x_0_x_0_z_z_z_xy, g_x_0_x_0_z_z_z_xz, g_x_0_x_0_z_z_z_yy, g_x_0_x_0_z_z_z_yz, g_x_0_x_0_z_z_z_zz, g_xz_z_xz_xx, g_xz_z_xz_xy, g_xz_z_xz_xz, g_xz_z_xz_yy, g_xz_z_xz_yz, g_xz_z_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_z_z_xx[i] = 4.0 * g_xz_z_xz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_z_xy[i] = 4.0 * g_xz_z_xz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_z_xz[i] = 4.0 * g_xz_z_xz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_z_yy[i] = 4.0 * g_xz_z_xz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_z_yz[i] = 4.0 * g_xz_z_xz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_z_zz[i] = 4.0 * g_xz_z_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_0_x_xy_xx, g_0_x_xy_xy, g_0_x_xy_xz, g_0_x_xy_yy, g_0_x_xy_yz, g_0_x_xy_zz, g_x_0_y_0_x_x_x_xx, g_x_0_y_0_x_x_x_xy, g_x_0_y_0_x_x_x_xz, g_x_0_y_0_x_x_x_yy, g_x_0_y_0_x_x_x_yz, g_x_0_y_0_x_x_x_zz, g_xx_x_xy_xx, g_xx_x_xy_xy, g_xx_x_xy_xz, g_xx_x_xy_yy, g_xx_x_xy_yz, g_xx_x_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_x_x_xx[i] = -2.0 * g_0_x_xy_xx[i] * c_exps[i] + 4.0 * g_xx_x_xy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_x_xy[i] = -2.0 * g_0_x_xy_xy[i] * c_exps[i] + 4.0 * g_xx_x_xy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_x_xz[i] = -2.0 * g_0_x_xy_xz[i] * c_exps[i] + 4.0 * g_xx_x_xy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_x_yy[i] = -2.0 * g_0_x_xy_yy[i] * c_exps[i] + 4.0 * g_xx_x_xy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_x_yz[i] = -2.0 * g_0_x_xy_yz[i] * c_exps[i] + 4.0 * g_xx_x_xy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_x_zz[i] = -2.0 * g_0_x_xy_zz[i] * c_exps[i] + 4.0 * g_xx_x_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_0_x_0_xx, g_0_x_0_xy, g_0_x_0_xz, g_0_x_0_yy, g_0_x_0_yz, g_0_x_0_zz, g_0_x_yy_xx, g_0_x_yy_xy, g_0_x_yy_xz, g_0_x_yy_yy, g_0_x_yy_yz, g_0_x_yy_zz, g_x_0_y_0_x_x_y_xx, g_x_0_y_0_x_x_y_xy, g_x_0_y_0_x_x_y_xz, g_x_0_y_0_x_x_y_yy, g_x_0_y_0_x_x_y_yz, g_x_0_y_0_x_x_y_zz, g_xx_x_0_xx, g_xx_x_0_xy, g_xx_x_0_xz, g_xx_x_0_yy, g_xx_x_0_yz, g_xx_x_0_zz, g_xx_x_yy_xx, g_xx_x_yy_xy, g_xx_x_yy_xz, g_xx_x_yy_yy, g_xx_x_yy_yz, g_xx_x_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_x_y_xx[i] = g_0_x_0_xx[i] - 2.0 * g_0_x_yy_xx[i] * c_exps[i] - 2.0 * g_xx_x_0_xx[i] * a_exp + 4.0 * g_xx_x_yy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_y_xy[i] = g_0_x_0_xy[i] - 2.0 * g_0_x_yy_xy[i] * c_exps[i] - 2.0 * g_xx_x_0_xy[i] * a_exp + 4.0 * g_xx_x_yy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_y_xz[i] = g_0_x_0_xz[i] - 2.0 * g_0_x_yy_xz[i] * c_exps[i] - 2.0 * g_xx_x_0_xz[i] * a_exp + 4.0 * g_xx_x_yy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_y_yy[i] = g_0_x_0_yy[i] - 2.0 * g_0_x_yy_yy[i] * c_exps[i] - 2.0 * g_xx_x_0_yy[i] * a_exp + 4.0 * g_xx_x_yy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_y_yz[i] = g_0_x_0_yz[i] - 2.0 * g_0_x_yy_yz[i] * c_exps[i] - 2.0 * g_xx_x_0_yz[i] * a_exp + 4.0 * g_xx_x_yy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_y_zz[i] = g_0_x_0_zz[i] - 2.0 * g_0_x_yy_zz[i] * c_exps[i] - 2.0 * g_xx_x_0_zz[i] * a_exp + 4.0 * g_xx_x_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_0_x_yz_xx, g_0_x_yz_xy, g_0_x_yz_xz, g_0_x_yz_yy, g_0_x_yz_yz, g_0_x_yz_zz, g_x_0_y_0_x_x_z_xx, g_x_0_y_0_x_x_z_xy, g_x_0_y_0_x_x_z_xz, g_x_0_y_0_x_x_z_yy, g_x_0_y_0_x_x_z_yz, g_x_0_y_0_x_x_z_zz, g_xx_x_yz_xx, g_xx_x_yz_xy, g_xx_x_yz_xz, g_xx_x_yz_yy, g_xx_x_yz_yz, g_xx_x_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_x_z_xx[i] = -2.0 * g_0_x_yz_xx[i] * c_exps[i] + 4.0 * g_xx_x_yz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_z_xy[i] = -2.0 * g_0_x_yz_xy[i] * c_exps[i] + 4.0 * g_xx_x_yz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_z_xz[i] = -2.0 * g_0_x_yz_xz[i] * c_exps[i] + 4.0 * g_xx_x_yz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_z_yy[i] = -2.0 * g_0_x_yz_yy[i] * c_exps[i] + 4.0 * g_xx_x_yz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_z_yz[i] = -2.0 * g_0_x_yz_yz[i] * c_exps[i] + 4.0 * g_xx_x_yz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_z_zz[i] = -2.0 * g_0_x_yz_zz[i] * c_exps[i] + 4.0 * g_xx_x_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_0_y_xy_xx, g_0_y_xy_xy, g_0_y_xy_xz, g_0_y_xy_yy, g_0_y_xy_yz, g_0_y_xy_zz, g_x_0_y_0_x_y_x_xx, g_x_0_y_0_x_y_x_xy, g_x_0_y_0_x_y_x_xz, g_x_0_y_0_x_y_x_yy, g_x_0_y_0_x_y_x_yz, g_x_0_y_0_x_y_x_zz, g_xx_y_xy_xx, g_xx_y_xy_xy, g_xx_y_xy_xz, g_xx_y_xy_yy, g_xx_y_xy_yz, g_xx_y_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_y_x_xx[i] = -2.0 * g_0_y_xy_xx[i] * c_exps[i] + 4.0 * g_xx_y_xy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_x_xy[i] = -2.0 * g_0_y_xy_xy[i] * c_exps[i] + 4.0 * g_xx_y_xy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_x_xz[i] = -2.0 * g_0_y_xy_xz[i] * c_exps[i] + 4.0 * g_xx_y_xy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_x_yy[i] = -2.0 * g_0_y_xy_yy[i] * c_exps[i] + 4.0 * g_xx_y_xy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_x_yz[i] = -2.0 * g_0_y_xy_yz[i] * c_exps[i] + 4.0 * g_xx_y_xy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_x_zz[i] = -2.0 * g_0_y_xy_zz[i] * c_exps[i] + 4.0 * g_xx_y_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_0_y_0_xx, g_0_y_0_xy, g_0_y_0_xz, g_0_y_0_yy, g_0_y_0_yz, g_0_y_0_zz, g_0_y_yy_xx, g_0_y_yy_xy, g_0_y_yy_xz, g_0_y_yy_yy, g_0_y_yy_yz, g_0_y_yy_zz, g_x_0_y_0_x_y_y_xx, g_x_0_y_0_x_y_y_xy, g_x_0_y_0_x_y_y_xz, g_x_0_y_0_x_y_y_yy, g_x_0_y_0_x_y_y_yz, g_x_0_y_0_x_y_y_zz, g_xx_y_0_xx, g_xx_y_0_xy, g_xx_y_0_xz, g_xx_y_0_yy, g_xx_y_0_yz, g_xx_y_0_zz, g_xx_y_yy_xx, g_xx_y_yy_xy, g_xx_y_yy_xz, g_xx_y_yy_yy, g_xx_y_yy_yz, g_xx_y_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_y_y_xx[i] = g_0_y_0_xx[i] - 2.0 * g_0_y_yy_xx[i] * c_exps[i] - 2.0 * g_xx_y_0_xx[i] * a_exp + 4.0 * g_xx_y_yy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_y_xy[i] = g_0_y_0_xy[i] - 2.0 * g_0_y_yy_xy[i] * c_exps[i] - 2.0 * g_xx_y_0_xy[i] * a_exp + 4.0 * g_xx_y_yy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_y_xz[i] = g_0_y_0_xz[i] - 2.0 * g_0_y_yy_xz[i] * c_exps[i] - 2.0 * g_xx_y_0_xz[i] * a_exp + 4.0 * g_xx_y_yy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_y_yy[i] = g_0_y_0_yy[i] - 2.0 * g_0_y_yy_yy[i] * c_exps[i] - 2.0 * g_xx_y_0_yy[i] * a_exp + 4.0 * g_xx_y_yy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_y_yz[i] = g_0_y_0_yz[i] - 2.0 * g_0_y_yy_yz[i] * c_exps[i] - 2.0 * g_xx_y_0_yz[i] * a_exp + 4.0 * g_xx_y_yy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_y_zz[i] = g_0_y_0_zz[i] - 2.0 * g_0_y_yy_zz[i] * c_exps[i] - 2.0 * g_xx_y_0_zz[i] * a_exp + 4.0 * g_xx_y_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_0_y_yz_xx, g_0_y_yz_xy, g_0_y_yz_xz, g_0_y_yz_yy, g_0_y_yz_yz, g_0_y_yz_zz, g_x_0_y_0_x_y_z_xx, g_x_0_y_0_x_y_z_xy, g_x_0_y_0_x_y_z_xz, g_x_0_y_0_x_y_z_yy, g_x_0_y_0_x_y_z_yz, g_x_0_y_0_x_y_z_zz, g_xx_y_yz_xx, g_xx_y_yz_xy, g_xx_y_yz_xz, g_xx_y_yz_yy, g_xx_y_yz_yz, g_xx_y_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_y_z_xx[i] = -2.0 * g_0_y_yz_xx[i] * c_exps[i] + 4.0 * g_xx_y_yz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_z_xy[i] = -2.0 * g_0_y_yz_xy[i] * c_exps[i] + 4.0 * g_xx_y_yz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_z_xz[i] = -2.0 * g_0_y_yz_xz[i] * c_exps[i] + 4.0 * g_xx_y_yz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_z_yy[i] = -2.0 * g_0_y_yz_yy[i] * c_exps[i] + 4.0 * g_xx_y_yz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_z_yz[i] = -2.0 * g_0_y_yz_yz[i] * c_exps[i] + 4.0 * g_xx_y_yz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_z_zz[i] = -2.0 * g_0_y_yz_zz[i] * c_exps[i] + 4.0 * g_xx_y_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_0_z_xy_xx, g_0_z_xy_xy, g_0_z_xy_xz, g_0_z_xy_yy, g_0_z_xy_yz, g_0_z_xy_zz, g_x_0_y_0_x_z_x_xx, g_x_0_y_0_x_z_x_xy, g_x_0_y_0_x_z_x_xz, g_x_0_y_0_x_z_x_yy, g_x_0_y_0_x_z_x_yz, g_x_0_y_0_x_z_x_zz, g_xx_z_xy_xx, g_xx_z_xy_xy, g_xx_z_xy_xz, g_xx_z_xy_yy, g_xx_z_xy_yz, g_xx_z_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_z_x_xx[i] = -2.0 * g_0_z_xy_xx[i] * c_exps[i] + 4.0 * g_xx_z_xy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_x_xy[i] = -2.0 * g_0_z_xy_xy[i] * c_exps[i] + 4.0 * g_xx_z_xy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_x_xz[i] = -2.0 * g_0_z_xy_xz[i] * c_exps[i] + 4.0 * g_xx_z_xy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_x_yy[i] = -2.0 * g_0_z_xy_yy[i] * c_exps[i] + 4.0 * g_xx_z_xy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_x_yz[i] = -2.0 * g_0_z_xy_yz[i] * c_exps[i] + 4.0 * g_xx_z_xy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_x_zz[i] = -2.0 * g_0_z_xy_zz[i] * c_exps[i] + 4.0 * g_xx_z_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_0_z_0_xx, g_0_z_0_xy, g_0_z_0_xz, g_0_z_0_yy, g_0_z_0_yz, g_0_z_0_zz, g_0_z_yy_xx, g_0_z_yy_xy, g_0_z_yy_xz, g_0_z_yy_yy, g_0_z_yy_yz, g_0_z_yy_zz, g_x_0_y_0_x_z_y_xx, g_x_0_y_0_x_z_y_xy, g_x_0_y_0_x_z_y_xz, g_x_0_y_0_x_z_y_yy, g_x_0_y_0_x_z_y_yz, g_x_0_y_0_x_z_y_zz, g_xx_z_0_xx, g_xx_z_0_xy, g_xx_z_0_xz, g_xx_z_0_yy, g_xx_z_0_yz, g_xx_z_0_zz, g_xx_z_yy_xx, g_xx_z_yy_xy, g_xx_z_yy_xz, g_xx_z_yy_yy, g_xx_z_yy_yz, g_xx_z_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_z_y_xx[i] = g_0_z_0_xx[i] - 2.0 * g_0_z_yy_xx[i] * c_exps[i] - 2.0 * g_xx_z_0_xx[i] * a_exp + 4.0 * g_xx_z_yy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_y_xy[i] = g_0_z_0_xy[i] - 2.0 * g_0_z_yy_xy[i] * c_exps[i] - 2.0 * g_xx_z_0_xy[i] * a_exp + 4.0 * g_xx_z_yy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_y_xz[i] = g_0_z_0_xz[i] - 2.0 * g_0_z_yy_xz[i] * c_exps[i] - 2.0 * g_xx_z_0_xz[i] * a_exp + 4.0 * g_xx_z_yy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_y_yy[i] = g_0_z_0_yy[i] - 2.0 * g_0_z_yy_yy[i] * c_exps[i] - 2.0 * g_xx_z_0_yy[i] * a_exp + 4.0 * g_xx_z_yy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_y_yz[i] = g_0_z_0_yz[i] - 2.0 * g_0_z_yy_yz[i] * c_exps[i] - 2.0 * g_xx_z_0_yz[i] * a_exp + 4.0 * g_xx_z_yy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_y_zz[i] = g_0_z_0_zz[i] - 2.0 * g_0_z_yy_zz[i] * c_exps[i] - 2.0 * g_xx_z_0_zz[i] * a_exp + 4.0 * g_xx_z_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_0_z_yz_xx, g_0_z_yz_xy, g_0_z_yz_xz, g_0_z_yz_yy, g_0_z_yz_yz, g_0_z_yz_zz, g_x_0_y_0_x_z_z_xx, g_x_0_y_0_x_z_z_xy, g_x_0_y_0_x_z_z_xz, g_x_0_y_0_x_z_z_yy, g_x_0_y_0_x_z_z_yz, g_x_0_y_0_x_z_z_zz, g_xx_z_yz_xx, g_xx_z_yz_xy, g_xx_z_yz_xz, g_xx_z_yz_yy, g_xx_z_yz_yz, g_xx_z_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_z_z_xx[i] = -2.0 * g_0_z_yz_xx[i] * c_exps[i] + 4.0 * g_xx_z_yz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_z_xy[i] = -2.0 * g_0_z_yz_xy[i] * c_exps[i] + 4.0 * g_xx_z_yz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_z_xz[i] = -2.0 * g_0_z_yz_xz[i] * c_exps[i] + 4.0 * g_xx_z_yz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_z_yy[i] = -2.0 * g_0_z_yz_yy[i] * c_exps[i] + 4.0 * g_xx_z_yz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_z_yz[i] = -2.0 * g_0_z_yz_yz[i] * c_exps[i] + 4.0 * g_xx_z_yz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_z_zz[i] = -2.0 * g_0_z_yz_zz[i] * c_exps[i] + 4.0 * g_xx_z_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_x_0_y_0_y_x_x_xx, g_x_0_y_0_y_x_x_xy, g_x_0_y_0_y_x_x_xz, g_x_0_y_0_y_x_x_yy, g_x_0_y_0_y_x_x_yz, g_x_0_y_0_y_x_x_zz, g_xy_x_xy_xx, g_xy_x_xy_xy, g_xy_x_xy_xz, g_xy_x_xy_yy, g_xy_x_xy_yz, g_xy_x_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_x_x_xx[i] = 4.0 * g_xy_x_xy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_x_xy[i] = 4.0 * g_xy_x_xy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_x_xz[i] = 4.0 * g_xy_x_xy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_x_yy[i] = 4.0 * g_xy_x_xy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_x_yz[i] = 4.0 * g_xy_x_xy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_x_zz[i] = 4.0 * g_xy_x_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_x_0_y_0_y_x_y_xx, g_x_0_y_0_y_x_y_xy, g_x_0_y_0_y_x_y_xz, g_x_0_y_0_y_x_y_yy, g_x_0_y_0_y_x_y_yz, g_x_0_y_0_y_x_y_zz, g_xy_x_0_xx, g_xy_x_0_xy, g_xy_x_0_xz, g_xy_x_0_yy, g_xy_x_0_yz, g_xy_x_0_zz, g_xy_x_yy_xx, g_xy_x_yy_xy, g_xy_x_yy_xz, g_xy_x_yy_yy, g_xy_x_yy_yz, g_xy_x_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_x_y_xx[i] = -2.0 * g_xy_x_0_xx[i] * a_exp + 4.0 * g_xy_x_yy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_y_xy[i] = -2.0 * g_xy_x_0_xy[i] * a_exp + 4.0 * g_xy_x_yy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_y_xz[i] = -2.0 * g_xy_x_0_xz[i] * a_exp + 4.0 * g_xy_x_yy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_y_yy[i] = -2.0 * g_xy_x_0_yy[i] * a_exp + 4.0 * g_xy_x_yy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_y_yz[i] = -2.0 * g_xy_x_0_yz[i] * a_exp + 4.0 * g_xy_x_yy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_y_zz[i] = -2.0 * g_xy_x_0_zz[i] * a_exp + 4.0 * g_xy_x_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_x_0_y_0_y_x_z_xx, g_x_0_y_0_y_x_z_xy, g_x_0_y_0_y_x_z_xz, g_x_0_y_0_y_x_z_yy, g_x_0_y_0_y_x_z_yz, g_x_0_y_0_y_x_z_zz, g_xy_x_yz_xx, g_xy_x_yz_xy, g_xy_x_yz_xz, g_xy_x_yz_yy, g_xy_x_yz_yz, g_xy_x_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_x_z_xx[i] = 4.0 * g_xy_x_yz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_z_xy[i] = 4.0 * g_xy_x_yz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_z_xz[i] = 4.0 * g_xy_x_yz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_z_yy[i] = 4.0 * g_xy_x_yz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_z_yz[i] = 4.0 * g_xy_x_yz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_z_zz[i] = 4.0 * g_xy_x_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_x_0_y_0_y_y_x_xx, g_x_0_y_0_y_y_x_xy, g_x_0_y_0_y_y_x_xz, g_x_0_y_0_y_y_x_yy, g_x_0_y_0_y_y_x_yz, g_x_0_y_0_y_y_x_zz, g_xy_y_xy_xx, g_xy_y_xy_xy, g_xy_y_xy_xz, g_xy_y_xy_yy, g_xy_y_xy_yz, g_xy_y_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_y_x_xx[i] = 4.0 * g_xy_y_xy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_x_xy[i] = 4.0 * g_xy_y_xy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_x_xz[i] = 4.0 * g_xy_y_xy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_x_yy[i] = 4.0 * g_xy_y_xy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_x_yz[i] = 4.0 * g_xy_y_xy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_x_zz[i] = 4.0 * g_xy_y_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_x_0_y_0_y_y_y_xx, g_x_0_y_0_y_y_y_xy, g_x_0_y_0_y_y_y_xz, g_x_0_y_0_y_y_y_yy, g_x_0_y_0_y_y_y_yz, g_x_0_y_0_y_y_y_zz, g_xy_y_0_xx, g_xy_y_0_xy, g_xy_y_0_xz, g_xy_y_0_yy, g_xy_y_0_yz, g_xy_y_0_zz, g_xy_y_yy_xx, g_xy_y_yy_xy, g_xy_y_yy_xz, g_xy_y_yy_yy, g_xy_y_yy_yz, g_xy_y_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_y_y_xx[i] = -2.0 * g_xy_y_0_xx[i] * a_exp + 4.0 * g_xy_y_yy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_y_xy[i] = -2.0 * g_xy_y_0_xy[i] * a_exp + 4.0 * g_xy_y_yy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_y_xz[i] = -2.0 * g_xy_y_0_xz[i] * a_exp + 4.0 * g_xy_y_yy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_y_yy[i] = -2.0 * g_xy_y_0_yy[i] * a_exp + 4.0 * g_xy_y_yy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_y_yz[i] = -2.0 * g_xy_y_0_yz[i] * a_exp + 4.0 * g_xy_y_yy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_y_zz[i] = -2.0 * g_xy_y_0_zz[i] * a_exp + 4.0 * g_xy_y_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_x_0_y_0_y_y_z_xx, g_x_0_y_0_y_y_z_xy, g_x_0_y_0_y_y_z_xz, g_x_0_y_0_y_y_z_yy, g_x_0_y_0_y_y_z_yz, g_x_0_y_0_y_y_z_zz, g_xy_y_yz_xx, g_xy_y_yz_xy, g_xy_y_yz_xz, g_xy_y_yz_yy, g_xy_y_yz_yz, g_xy_y_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_y_z_xx[i] = 4.0 * g_xy_y_yz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_z_xy[i] = 4.0 * g_xy_y_yz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_z_xz[i] = 4.0 * g_xy_y_yz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_z_yy[i] = 4.0 * g_xy_y_yz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_z_yz[i] = 4.0 * g_xy_y_yz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_z_zz[i] = 4.0 * g_xy_y_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_x_0_y_0_y_z_x_xx, g_x_0_y_0_y_z_x_xy, g_x_0_y_0_y_z_x_xz, g_x_0_y_0_y_z_x_yy, g_x_0_y_0_y_z_x_yz, g_x_0_y_0_y_z_x_zz, g_xy_z_xy_xx, g_xy_z_xy_xy, g_xy_z_xy_xz, g_xy_z_xy_yy, g_xy_z_xy_yz, g_xy_z_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_z_x_xx[i] = 4.0 * g_xy_z_xy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_x_xy[i] = 4.0 * g_xy_z_xy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_x_xz[i] = 4.0 * g_xy_z_xy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_x_yy[i] = 4.0 * g_xy_z_xy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_x_yz[i] = 4.0 * g_xy_z_xy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_x_zz[i] = 4.0 * g_xy_z_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_x_0_y_0_y_z_y_xx, g_x_0_y_0_y_z_y_xy, g_x_0_y_0_y_z_y_xz, g_x_0_y_0_y_z_y_yy, g_x_0_y_0_y_z_y_yz, g_x_0_y_0_y_z_y_zz, g_xy_z_0_xx, g_xy_z_0_xy, g_xy_z_0_xz, g_xy_z_0_yy, g_xy_z_0_yz, g_xy_z_0_zz, g_xy_z_yy_xx, g_xy_z_yy_xy, g_xy_z_yy_xz, g_xy_z_yy_yy, g_xy_z_yy_yz, g_xy_z_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_z_y_xx[i] = -2.0 * g_xy_z_0_xx[i] * a_exp + 4.0 * g_xy_z_yy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_y_xy[i] = -2.0 * g_xy_z_0_xy[i] * a_exp + 4.0 * g_xy_z_yy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_y_xz[i] = -2.0 * g_xy_z_0_xz[i] * a_exp + 4.0 * g_xy_z_yy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_y_yy[i] = -2.0 * g_xy_z_0_yy[i] * a_exp + 4.0 * g_xy_z_yy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_y_yz[i] = -2.0 * g_xy_z_0_yz[i] * a_exp + 4.0 * g_xy_z_yy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_y_zz[i] = -2.0 * g_xy_z_0_zz[i] * a_exp + 4.0 * g_xy_z_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_x_0_y_0_y_z_z_xx, g_x_0_y_0_y_z_z_xy, g_x_0_y_0_y_z_z_xz, g_x_0_y_0_y_z_z_yy, g_x_0_y_0_y_z_z_yz, g_x_0_y_0_y_z_z_zz, g_xy_z_yz_xx, g_xy_z_yz_xy, g_xy_z_yz_xz, g_xy_z_yz_yy, g_xy_z_yz_yz, g_xy_z_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_z_z_xx[i] = 4.0 * g_xy_z_yz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_z_xy[i] = 4.0 * g_xy_z_yz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_z_xz[i] = 4.0 * g_xy_z_yz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_z_yy[i] = 4.0 * g_xy_z_yz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_z_yz[i] = 4.0 * g_xy_z_yz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_z_zz[i] = 4.0 * g_xy_z_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_x_0_y_0_z_x_x_xx, g_x_0_y_0_z_x_x_xy, g_x_0_y_0_z_x_x_xz, g_x_0_y_0_z_x_x_yy, g_x_0_y_0_z_x_x_yz, g_x_0_y_0_z_x_x_zz, g_xz_x_xy_xx, g_xz_x_xy_xy, g_xz_x_xy_xz, g_xz_x_xy_yy, g_xz_x_xy_yz, g_xz_x_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_x_x_xx[i] = 4.0 * g_xz_x_xy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_x_xy[i] = 4.0 * g_xz_x_xy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_x_xz[i] = 4.0 * g_xz_x_xy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_x_yy[i] = 4.0 * g_xz_x_xy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_x_yz[i] = 4.0 * g_xz_x_xy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_x_zz[i] = 4.0 * g_xz_x_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_x_0_y_0_z_x_y_xx, g_x_0_y_0_z_x_y_xy, g_x_0_y_0_z_x_y_xz, g_x_0_y_0_z_x_y_yy, g_x_0_y_0_z_x_y_yz, g_x_0_y_0_z_x_y_zz, g_xz_x_0_xx, g_xz_x_0_xy, g_xz_x_0_xz, g_xz_x_0_yy, g_xz_x_0_yz, g_xz_x_0_zz, g_xz_x_yy_xx, g_xz_x_yy_xy, g_xz_x_yy_xz, g_xz_x_yy_yy, g_xz_x_yy_yz, g_xz_x_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_x_y_xx[i] = -2.0 * g_xz_x_0_xx[i] * a_exp + 4.0 * g_xz_x_yy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_y_xy[i] = -2.0 * g_xz_x_0_xy[i] * a_exp + 4.0 * g_xz_x_yy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_y_xz[i] = -2.0 * g_xz_x_0_xz[i] * a_exp + 4.0 * g_xz_x_yy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_y_yy[i] = -2.0 * g_xz_x_0_yy[i] * a_exp + 4.0 * g_xz_x_yy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_y_yz[i] = -2.0 * g_xz_x_0_yz[i] * a_exp + 4.0 * g_xz_x_yy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_y_zz[i] = -2.0 * g_xz_x_0_zz[i] * a_exp + 4.0 * g_xz_x_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_x_0_y_0_z_x_z_xx, g_x_0_y_0_z_x_z_xy, g_x_0_y_0_z_x_z_xz, g_x_0_y_0_z_x_z_yy, g_x_0_y_0_z_x_z_yz, g_x_0_y_0_z_x_z_zz, g_xz_x_yz_xx, g_xz_x_yz_xy, g_xz_x_yz_xz, g_xz_x_yz_yy, g_xz_x_yz_yz, g_xz_x_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_x_z_xx[i] = 4.0 * g_xz_x_yz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_z_xy[i] = 4.0 * g_xz_x_yz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_z_xz[i] = 4.0 * g_xz_x_yz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_z_yy[i] = 4.0 * g_xz_x_yz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_z_yz[i] = 4.0 * g_xz_x_yz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_z_zz[i] = 4.0 * g_xz_x_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_x_0_y_0_z_y_x_xx, g_x_0_y_0_z_y_x_xy, g_x_0_y_0_z_y_x_xz, g_x_0_y_0_z_y_x_yy, g_x_0_y_0_z_y_x_yz, g_x_0_y_0_z_y_x_zz, g_xz_y_xy_xx, g_xz_y_xy_xy, g_xz_y_xy_xz, g_xz_y_xy_yy, g_xz_y_xy_yz, g_xz_y_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_y_x_xx[i] = 4.0 * g_xz_y_xy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_x_xy[i] = 4.0 * g_xz_y_xy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_x_xz[i] = 4.0 * g_xz_y_xy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_x_yy[i] = 4.0 * g_xz_y_xy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_x_yz[i] = 4.0 * g_xz_y_xy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_x_zz[i] = 4.0 * g_xz_y_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_x_0_y_0_z_y_y_xx, g_x_0_y_0_z_y_y_xy, g_x_0_y_0_z_y_y_xz, g_x_0_y_0_z_y_y_yy, g_x_0_y_0_z_y_y_yz, g_x_0_y_0_z_y_y_zz, g_xz_y_0_xx, g_xz_y_0_xy, g_xz_y_0_xz, g_xz_y_0_yy, g_xz_y_0_yz, g_xz_y_0_zz, g_xz_y_yy_xx, g_xz_y_yy_xy, g_xz_y_yy_xz, g_xz_y_yy_yy, g_xz_y_yy_yz, g_xz_y_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_y_y_xx[i] = -2.0 * g_xz_y_0_xx[i] * a_exp + 4.0 * g_xz_y_yy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_y_xy[i] = -2.0 * g_xz_y_0_xy[i] * a_exp + 4.0 * g_xz_y_yy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_y_xz[i] = -2.0 * g_xz_y_0_xz[i] * a_exp + 4.0 * g_xz_y_yy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_y_yy[i] = -2.0 * g_xz_y_0_yy[i] * a_exp + 4.0 * g_xz_y_yy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_y_yz[i] = -2.0 * g_xz_y_0_yz[i] * a_exp + 4.0 * g_xz_y_yy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_y_zz[i] = -2.0 * g_xz_y_0_zz[i] * a_exp + 4.0 * g_xz_y_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_x_0_y_0_z_y_z_xx, g_x_0_y_0_z_y_z_xy, g_x_0_y_0_z_y_z_xz, g_x_0_y_0_z_y_z_yy, g_x_0_y_0_z_y_z_yz, g_x_0_y_0_z_y_z_zz, g_xz_y_yz_xx, g_xz_y_yz_xy, g_xz_y_yz_xz, g_xz_y_yz_yy, g_xz_y_yz_yz, g_xz_y_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_y_z_xx[i] = 4.0 * g_xz_y_yz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_z_xy[i] = 4.0 * g_xz_y_yz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_z_xz[i] = 4.0 * g_xz_y_yz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_z_yy[i] = 4.0 * g_xz_y_yz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_z_yz[i] = 4.0 * g_xz_y_yz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_z_zz[i] = 4.0 * g_xz_y_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_x_0_y_0_z_z_x_xx, g_x_0_y_0_z_z_x_xy, g_x_0_y_0_z_z_x_xz, g_x_0_y_0_z_z_x_yy, g_x_0_y_0_z_z_x_yz, g_x_0_y_0_z_z_x_zz, g_xz_z_xy_xx, g_xz_z_xy_xy, g_xz_z_xy_xz, g_xz_z_xy_yy, g_xz_z_xy_yz, g_xz_z_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_z_x_xx[i] = 4.0 * g_xz_z_xy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_x_xy[i] = 4.0 * g_xz_z_xy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_x_xz[i] = 4.0 * g_xz_z_xy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_x_yy[i] = 4.0 * g_xz_z_xy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_x_yz[i] = 4.0 * g_xz_z_xy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_x_zz[i] = 4.0 * g_xz_z_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_x_0_y_0_z_z_y_xx, g_x_0_y_0_z_z_y_xy, g_x_0_y_0_z_z_y_xz, g_x_0_y_0_z_z_y_yy, g_x_0_y_0_z_z_y_yz, g_x_0_y_0_z_z_y_zz, g_xz_z_0_xx, g_xz_z_0_xy, g_xz_z_0_xz, g_xz_z_0_yy, g_xz_z_0_yz, g_xz_z_0_zz, g_xz_z_yy_xx, g_xz_z_yy_xy, g_xz_z_yy_xz, g_xz_z_yy_yy, g_xz_z_yy_yz, g_xz_z_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_z_y_xx[i] = -2.0 * g_xz_z_0_xx[i] * a_exp + 4.0 * g_xz_z_yy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_y_xy[i] = -2.0 * g_xz_z_0_xy[i] * a_exp + 4.0 * g_xz_z_yy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_y_xz[i] = -2.0 * g_xz_z_0_xz[i] * a_exp + 4.0 * g_xz_z_yy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_y_yy[i] = -2.0 * g_xz_z_0_yy[i] * a_exp + 4.0 * g_xz_z_yy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_y_yz[i] = -2.0 * g_xz_z_0_yz[i] * a_exp + 4.0 * g_xz_z_yy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_y_zz[i] = -2.0 * g_xz_z_0_zz[i] * a_exp + 4.0 * g_xz_z_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_x_0_y_0_z_z_z_xx, g_x_0_y_0_z_z_z_xy, g_x_0_y_0_z_z_z_xz, g_x_0_y_0_z_z_z_yy, g_x_0_y_0_z_z_z_yz, g_x_0_y_0_z_z_z_zz, g_xz_z_yz_xx, g_xz_z_yz_xy, g_xz_z_yz_xz, g_xz_z_yz_yy, g_xz_z_yz_yz, g_xz_z_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_z_z_xx[i] = 4.0 * g_xz_z_yz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_z_xy[i] = 4.0 * g_xz_z_yz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_z_xz[i] = 4.0 * g_xz_z_yz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_z_yy[i] = 4.0 * g_xz_z_yz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_z_yz[i] = 4.0 * g_xz_z_yz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_z_zz[i] = 4.0 * g_xz_z_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (324-330)

    #pragma omp simd aligned(g_0_x_xz_xx, g_0_x_xz_xy, g_0_x_xz_xz, g_0_x_xz_yy, g_0_x_xz_yz, g_0_x_xz_zz, g_x_0_z_0_x_x_x_xx, g_x_0_z_0_x_x_x_xy, g_x_0_z_0_x_x_x_xz, g_x_0_z_0_x_x_x_yy, g_x_0_z_0_x_x_x_yz, g_x_0_z_0_x_x_x_zz, g_xx_x_xz_xx, g_xx_x_xz_xy, g_xx_x_xz_xz, g_xx_x_xz_yy, g_xx_x_xz_yz, g_xx_x_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_x_x_xx[i] = -2.0 * g_0_x_xz_xx[i] * c_exps[i] + 4.0 * g_xx_x_xz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_x_xy[i] = -2.0 * g_0_x_xz_xy[i] * c_exps[i] + 4.0 * g_xx_x_xz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_x_xz[i] = -2.0 * g_0_x_xz_xz[i] * c_exps[i] + 4.0 * g_xx_x_xz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_x_yy[i] = -2.0 * g_0_x_xz_yy[i] * c_exps[i] + 4.0 * g_xx_x_xz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_x_yz[i] = -2.0 * g_0_x_xz_yz[i] * c_exps[i] + 4.0 * g_xx_x_xz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_x_zz[i] = -2.0 * g_0_x_xz_zz[i] * c_exps[i] + 4.0 * g_xx_x_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (330-336)

    #pragma omp simd aligned(g_0_x_yz_xx, g_0_x_yz_xy, g_0_x_yz_xz, g_0_x_yz_yy, g_0_x_yz_yz, g_0_x_yz_zz, g_x_0_z_0_x_x_y_xx, g_x_0_z_0_x_x_y_xy, g_x_0_z_0_x_x_y_xz, g_x_0_z_0_x_x_y_yy, g_x_0_z_0_x_x_y_yz, g_x_0_z_0_x_x_y_zz, g_xx_x_yz_xx, g_xx_x_yz_xy, g_xx_x_yz_xz, g_xx_x_yz_yy, g_xx_x_yz_yz, g_xx_x_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_x_y_xx[i] = -2.0 * g_0_x_yz_xx[i] * c_exps[i] + 4.0 * g_xx_x_yz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_y_xy[i] = -2.0 * g_0_x_yz_xy[i] * c_exps[i] + 4.0 * g_xx_x_yz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_y_xz[i] = -2.0 * g_0_x_yz_xz[i] * c_exps[i] + 4.0 * g_xx_x_yz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_y_yy[i] = -2.0 * g_0_x_yz_yy[i] * c_exps[i] + 4.0 * g_xx_x_yz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_y_yz[i] = -2.0 * g_0_x_yz_yz[i] * c_exps[i] + 4.0 * g_xx_x_yz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_y_zz[i] = -2.0 * g_0_x_yz_zz[i] * c_exps[i] + 4.0 * g_xx_x_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (336-342)

    #pragma omp simd aligned(g_0_x_0_xx, g_0_x_0_xy, g_0_x_0_xz, g_0_x_0_yy, g_0_x_0_yz, g_0_x_0_zz, g_0_x_zz_xx, g_0_x_zz_xy, g_0_x_zz_xz, g_0_x_zz_yy, g_0_x_zz_yz, g_0_x_zz_zz, g_x_0_z_0_x_x_z_xx, g_x_0_z_0_x_x_z_xy, g_x_0_z_0_x_x_z_xz, g_x_0_z_0_x_x_z_yy, g_x_0_z_0_x_x_z_yz, g_x_0_z_0_x_x_z_zz, g_xx_x_0_xx, g_xx_x_0_xy, g_xx_x_0_xz, g_xx_x_0_yy, g_xx_x_0_yz, g_xx_x_0_zz, g_xx_x_zz_xx, g_xx_x_zz_xy, g_xx_x_zz_xz, g_xx_x_zz_yy, g_xx_x_zz_yz, g_xx_x_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_x_z_xx[i] = g_0_x_0_xx[i] - 2.0 * g_0_x_zz_xx[i] * c_exps[i] - 2.0 * g_xx_x_0_xx[i] * a_exp + 4.0 * g_xx_x_zz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_z_xy[i] = g_0_x_0_xy[i] - 2.0 * g_0_x_zz_xy[i] * c_exps[i] - 2.0 * g_xx_x_0_xy[i] * a_exp + 4.0 * g_xx_x_zz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_z_xz[i] = g_0_x_0_xz[i] - 2.0 * g_0_x_zz_xz[i] * c_exps[i] - 2.0 * g_xx_x_0_xz[i] * a_exp + 4.0 * g_xx_x_zz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_z_yy[i] = g_0_x_0_yy[i] - 2.0 * g_0_x_zz_yy[i] * c_exps[i] - 2.0 * g_xx_x_0_yy[i] * a_exp + 4.0 * g_xx_x_zz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_z_yz[i] = g_0_x_0_yz[i] - 2.0 * g_0_x_zz_yz[i] * c_exps[i] - 2.0 * g_xx_x_0_yz[i] * a_exp + 4.0 * g_xx_x_zz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_z_zz[i] = g_0_x_0_zz[i] - 2.0 * g_0_x_zz_zz[i] * c_exps[i] - 2.0 * g_xx_x_0_zz[i] * a_exp + 4.0 * g_xx_x_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (342-348)

    #pragma omp simd aligned(g_0_y_xz_xx, g_0_y_xz_xy, g_0_y_xz_xz, g_0_y_xz_yy, g_0_y_xz_yz, g_0_y_xz_zz, g_x_0_z_0_x_y_x_xx, g_x_0_z_0_x_y_x_xy, g_x_0_z_0_x_y_x_xz, g_x_0_z_0_x_y_x_yy, g_x_0_z_0_x_y_x_yz, g_x_0_z_0_x_y_x_zz, g_xx_y_xz_xx, g_xx_y_xz_xy, g_xx_y_xz_xz, g_xx_y_xz_yy, g_xx_y_xz_yz, g_xx_y_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_y_x_xx[i] = -2.0 * g_0_y_xz_xx[i] * c_exps[i] + 4.0 * g_xx_y_xz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_x_xy[i] = -2.0 * g_0_y_xz_xy[i] * c_exps[i] + 4.0 * g_xx_y_xz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_x_xz[i] = -2.0 * g_0_y_xz_xz[i] * c_exps[i] + 4.0 * g_xx_y_xz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_x_yy[i] = -2.0 * g_0_y_xz_yy[i] * c_exps[i] + 4.0 * g_xx_y_xz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_x_yz[i] = -2.0 * g_0_y_xz_yz[i] * c_exps[i] + 4.0 * g_xx_y_xz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_x_zz[i] = -2.0 * g_0_y_xz_zz[i] * c_exps[i] + 4.0 * g_xx_y_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (348-354)

    #pragma omp simd aligned(g_0_y_yz_xx, g_0_y_yz_xy, g_0_y_yz_xz, g_0_y_yz_yy, g_0_y_yz_yz, g_0_y_yz_zz, g_x_0_z_0_x_y_y_xx, g_x_0_z_0_x_y_y_xy, g_x_0_z_0_x_y_y_xz, g_x_0_z_0_x_y_y_yy, g_x_0_z_0_x_y_y_yz, g_x_0_z_0_x_y_y_zz, g_xx_y_yz_xx, g_xx_y_yz_xy, g_xx_y_yz_xz, g_xx_y_yz_yy, g_xx_y_yz_yz, g_xx_y_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_y_y_xx[i] = -2.0 * g_0_y_yz_xx[i] * c_exps[i] + 4.0 * g_xx_y_yz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_y_xy[i] = -2.0 * g_0_y_yz_xy[i] * c_exps[i] + 4.0 * g_xx_y_yz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_y_xz[i] = -2.0 * g_0_y_yz_xz[i] * c_exps[i] + 4.0 * g_xx_y_yz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_y_yy[i] = -2.0 * g_0_y_yz_yy[i] * c_exps[i] + 4.0 * g_xx_y_yz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_y_yz[i] = -2.0 * g_0_y_yz_yz[i] * c_exps[i] + 4.0 * g_xx_y_yz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_y_zz[i] = -2.0 * g_0_y_yz_zz[i] * c_exps[i] + 4.0 * g_xx_y_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (354-360)

    #pragma omp simd aligned(g_0_y_0_xx, g_0_y_0_xy, g_0_y_0_xz, g_0_y_0_yy, g_0_y_0_yz, g_0_y_0_zz, g_0_y_zz_xx, g_0_y_zz_xy, g_0_y_zz_xz, g_0_y_zz_yy, g_0_y_zz_yz, g_0_y_zz_zz, g_x_0_z_0_x_y_z_xx, g_x_0_z_0_x_y_z_xy, g_x_0_z_0_x_y_z_xz, g_x_0_z_0_x_y_z_yy, g_x_0_z_0_x_y_z_yz, g_x_0_z_0_x_y_z_zz, g_xx_y_0_xx, g_xx_y_0_xy, g_xx_y_0_xz, g_xx_y_0_yy, g_xx_y_0_yz, g_xx_y_0_zz, g_xx_y_zz_xx, g_xx_y_zz_xy, g_xx_y_zz_xz, g_xx_y_zz_yy, g_xx_y_zz_yz, g_xx_y_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_y_z_xx[i] = g_0_y_0_xx[i] - 2.0 * g_0_y_zz_xx[i] * c_exps[i] - 2.0 * g_xx_y_0_xx[i] * a_exp + 4.0 * g_xx_y_zz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_z_xy[i] = g_0_y_0_xy[i] - 2.0 * g_0_y_zz_xy[i] * c_exps[i] - 2.0 * g_xx_y_0_xy[i] * a_exp + 4.0 * g_xx_y_zz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_z_xz[i] = g_0_y_0_xz[i] - 2.0 * g_0_y_zz_xz[i] * c_exps[i] - 2.0 * g_xx_y_0_xz[i] * a_exp + 4.0 * g_xx_y_zz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_z_yy[i] = g_0_y_0_yy[i] - 2.0 * g_0_y_zz_yy[i] * c_exps[i] - 2.0 * g_xx_y_0_yy[i] * a_exp + 4.0 * g_xx_y_zz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_z_yz[i] = g_0_y_0_yz[i] - 2.0 * g_0_y_zz_yz[i] * c_exps[i] - 2.0 * g_xx_y_0_yz[i] * a_exp + 4.0 * g_xx_y_zz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_z_zz[i] = g_0_y_0_zz[i] - 2.0 * g_0_y_zz_zz[i] * c_exps[i] - 2.0 * g_xx_y_0_zz[i] * a_exp + 4.0 * g_xx_y_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (360-366)

    #pragma omp simd aligned(g_0_z_xz_xx, g_0_z_xz_xy, g_0_z_xz_xz, g_0_z_xz_yy, g_0_z_xz_yz, g_0_z_xz_zz, g_x_0_z_0_x_z_x_xx, g_x_0_z_0_x_z_x_xy, g_x_0_z_0_x_z_x_xz, g_x_0_z_0_x_z_x_yy, g_x_0_z_0_x_z_x_yz, g_x_0_z_0_x_z_x_zz, g_xx_z_xz_xx, g_xx_z_xz_xy, g_xx_z_xz_xz, g_xx_z_xz_yy, g_xx_z_xz_yz, g_xx_z_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_z_x_xx[i] = -2.0 * g_0_z_xz_xx[i] * c_exps[i] + 4.0 * g_xx_z_xz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_x_xy[i] = -2.0 * g_0_z_xz_xy[i] * c_exps[i] + 4.0 * g_xx_z_xz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_x_xz[i] = -2.0 * g_0_z_xz_xz[i] * c_exps[i] + 4.0 * g_xx_z_xz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_x_yy[i] = -2.0 * g_0_z_xz_yy[i] * c_exps[i] + 4.0 * g_xx_z_xz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_x_yz[i] = -2.0 * g_0_z_xz_yz[i] * c_exps[i] + 4.0 * g_xx_z_xz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_x_zz[i] = -2.0 * g_0_z_xz_zz[i] * c_exps[i] + 4.0 * g_xx_z_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (366-372)

    #pragma omp simd aligned(g_0_z_yz_xx, g_0_z_yz_xy, g_0_z_yz_xz, g_0_z_yz_yy, g_0_z_yz_yz, g_0_z_yz_zz, g_x_0_z_0_x_z_y_xx, g_x_0_z_0_x_z_y_xy, g_x_0_z_0_x_z_y_xz, g_x_0_z_0_x_z_y_yy, g_x_0_z_0_x_z_y_yz, g_x_0_z_0_x_z_y_zz, g_xx_z_yz_xx, g_xx_z_yz_xy, g_xx_z_yz_xz, g_xx_z_yz_yy, g_xx_z_yz_yz, g_xx_z_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_z_y_xx[i] = -2.0 * g_0_z_yz_xx[i] * c_exps[i] + 4.0 * g_xx_z_yz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_y_xy[i] = -2.0 * g_0_z_yz_xy[i] * c_exps[i] + 4.0 * g_xx_z_yz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_y_xz[i] = -2.0 * g_0_z_yz_xz[i] * c_exps[i] + 4.0 * g_xx_z_yz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_y_yy[i] = -2.0 * g_0_z_yz_yy[i] * c_exps[i] + 4.0 * g_xx_z_yz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_y_yz[i] = -2.0 * g_0_z_yz_yz[i] * c_exps[i] + 4.0 * g_xx_z_yz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_y_zz[i] = -2.0 * g_0_z_yz_zz[i] * c_exps[i] + 4.0 * g_xx_z_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (372-378)

    #pragma omp simd aligned(g_0_z_0_xx, g_0_z_0_xy, g_0_z_0_xz, g_0_z_0_yy, g_0_z_0_yz, g_0_z_0_zz, g_0_z_zz_xx, g_0_z_zz_xy, g_0_z_zz_xz, g_0_z_zz_yy, g_0_z_zz_yz, g_0_z_zz_zz, g_x_0_z_0_x_z_z_xx, g_x_0_z_0_x_z_z_xy, g_x_0_z_0_x_z_z_xz, g_x_0_z_0_x_z_z_yy, g_x_0_z_0_x_z_z_yz, g_x_0_z_0_x_z_z_zz, g_xx_z_0_xx, g_xx_z_0_xy, g_xx_z_0_xz, g_xx_z_0_yy, g_xx_z_0_yz, g_xx_z_0_zz, g_xx_z_zz_xx, g_xx_z_zz_xy, g_xx_z_zz_xz, g_xx_z_zz_yy, g_xx_z_zz_yz, g_xx_z_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_z_z_xx[i] = g_0_z_0_xx[i] - 2.0 * g_0_z_zz_xx[i] * c_exps[i] - 2.0 * g_xx_z_0_xx[i] * a_exp + 4.0 * g_xx_z_zz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_z_xy[i] = g_0_z_0_xy[i] - 2.0 * g_0_z_zz_xy[i] * c_exps[i] - 2.0 * g_xx_z_0_xy[i] * a_exp + 4.0 * g_xx_z_zz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_z_xz[i] = g_0_z_0_xz[i] - 2.0 * g_0_z_zz_xz[i] * c_exps[i] - 2.0 * g_xx_z_0_xz[i] * a_exp + 4.0 * g_xx_z_zz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_z_yy[i] = g_0_z_0_yy[i] - 2.0 * g_0_z_zz_yy[i] * c_exps[i] - 2.0 * g_xx_z_0_yy[i] * a_exp + 4.0 * g_xx_z_zz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_z_yz[i] = g_0_z_0_yz[i] - 2.0 * g_0_z_zz_yz[i] * c_exps[i] - 2.0 * g_xx_z_0_yz[i] * a_exp + 4.0 * g_xx_z_zz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_z_zz[i] = g_0_z_0_zz[i] - 2.0 * g_0_z_zz_zz[i] * c_exps[i] - 2.0 * g_xx_z_0_zz[i] * a_exp + 4.0 * g_xx_z_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (378-384)

    #pragma omp simd aligned(g_x_0_z_0_y_x_x_xx, g_x_0_z_0_y_x_x_xy, g_x_0_z_0_y_x_x_xz, g_x_0_z_0_y_x_x_yy, g_x_0_z_0_y_x_x_yz, g_x_0_z_0_y_x_x_zz, g_xy_x_xz_xx, g_xy_x_xz_xy, g_xy_x_xz_xz, g_xy_x_xz_yy, g_xy_x_xz_yz, g_xy_x_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_x_x_xx[i] = 4.0 * g_xy_x_xz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_x_xy[i] = 4.0 * g_xy_x_xz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_x_xz[i] = 4.0 * g_xy_x_xz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_x_yy[i] = 4.0 * g_xy_x_xz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_x_yz[i] = 4.0 * g_xy_x_xz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_x_zz[i] = 4.0 * g_xy_x_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (384-390)

    #pragma omp simd aligned(g_x_0_z_0_y_x_y_xx, g_x_0_z_0_y_x_y_xy, g_x_0_z_0_y_x_y_xz, g_x_0_z_0_y_x_y_yy, g_x_0_z_0_y_x_y_yz, g_x_0_z_0_y_x_y_zz, g_xy_x_yz_xx, g_xy_x_yz_xy, g_xy_x_yz_xz, g_xy_x_yz_yy, g_xy_x_yz_yz, g_xy_x_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_x_y_xx[i] = 4.0 * g_xy_x_yz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_y_xy[i] = 4.0 * g_xy_x_yz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_y_xz[i] = 4.0 * g_xy_x_yz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_y_yy[i] = 4.0 * g_xy_x_yz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_y_yz[i] = 4.0 * g_xy_x_yz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_y_zz[i] = 4.0 * g_xy_x_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (390-396)

    #pragma omp simd aligned(g_x_0_z_0_y_x_z_xx, g_x_0_z_0_y_x_z_xy, g_x_0_z_0_y_x_z_xz, g_x_0_z_0_y_x_z_yy, g_x_0_z_0_y_x_z_yz, g_x_0_z_0_y_x_z_zz, g_xy_x_0_xx, g_xy_x_0_xy, g_xy_x_0_xz, g_xy_x_0_yy, g_xy_x_0_yz, g_xy_x_0_zz, g_xy_x_zz_xx, g_xy_x_zz_xy, g_xy_x_zz_xz, g_xy_x_zz_yy, g_xy_x_zz_yz, g_xy_x_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_x_z_xx[i] = -2.0 * g_xy_x_0_xx[i] * a_exp + 4.0 * g_xy_x_zz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_z_xy[i] = -2.0 * g_xy_x_0_xy[i] * a_exp + 4.0 * g_xy_x_zz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_z_xz[i] = -2.0 * g_xy_x_0_xz[i] * a_exp + 4.0 * g_xy_x_zz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_z_yy[i] = -2.0 * g_xy_x_0_yy[i] * a_exp + 4.0 * g_xy_x_zz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_z_yz[i] = -2.0 * g_xy_x_0_yz[i] * a_exp + 4.0 * g_xy_x_zz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_z_zz[i] = -2.0 * g_xy_x_0_zz[i] * a_exp + 4.0 * g_xy_x_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (396-402)

    #pragma omp simd aligned(g_x_0_z_0_y_y_x_xx, g_x_0_z_0_y_y_x_xy, g_x_0_z_0_y_y_x_xz, g_x_0_z_0_y_y_x_yy, g_x_0_z_0_y_y_x_yz, g_x_0_z_0_y_y_x_zz, g_xy_y_xz_xx, g_xy_y_xz_xy, g_xy_y_xz_xz, g_xy_y_xz_yy, g_xy_y_xz_yz, g_xy_y_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_y_x_xx[i] = 4.0 * g_xy_y_xz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_x_xy[i] = 4.0 * g_xy_y_xz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_x_xz[i] = 4.0 * g_xy_y_xz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_x_yy[i] = 4.0 * g_xy_y_xz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_x_yz[i] = 4.0 * g_xy_y_xz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_x_zz[i] = 4.0 * g_xy_y_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (402-408)

    #pragma omp simd aligned(g_x_0_z_0_y_y_y_xx, g_x_0_z_0_y_y_y_xy, g_x_0_z_0_y_y_y_xz, g_x_0_z_0_y_y_y_yy, g_x_0_z_0_y_y_y_yz, g_x_0_z_0_y_y_y_zz, g_xy_y_yz_xx, g_xy_y_yz_xy, g_xy_y_yz_xz, g_xy_y_yz_yy, g_xy_y_yz_yz, g_xy_y_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_y_y_xx[i] = 4.0 * g_xy_y_yz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_y_xy[i] = 4.0 * g_xy_y_yz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_y_xz[i] = 4.0 * g_xy_y_yz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_y_yy[i] = 4.0 * g_xy_y_yz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_y_yz[i] = 4.0 * g_xy_y_yz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_y_zz[i] = 4.0 * g_xy_y_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (408-414)

    #pragma omp simd aligned(g_x_0_z_0_y_y_z_xx, g_x_0_z_0_y_y_z_xy, g_x_0_z_0_y_y_z_xz, g_x_0_z_0_y_y_z_yy, g_x_0_z_0_y_y_z_yz, g_x_0_z_0_y_y_z_zz, g_xy_y_0_xx, g_xy_y_0_xy, g_xy_y_0_xz, g_xy_y_0_yy, g_xy_y_0_yz, g_xy_y_0_zz, g_xy_y_zz_xx, g_xy_y_zz_xy, g_xy_y_zz_xz, g_xy_y_zz_yy, g_xy_y_zz_yz, g_xy_y_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_y_z_xx[i] = -2.0 * g_xy_y_0_xx[i] * a_exp + 4.0 * g_xy_y_zz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_z_xy[i] = -2.0 * g_xy_y_0_xy[i] * a_exp + 4.0 * g_xy_y_zz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_z_xz[i] = -2.0 * g_xy_y_0_xz[i] * a_exp + 4.0 * g_xy_y_zz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_z_yy[i] = -2.0 * g_xy_y_0_yy[i] * a_exp + 4.0 * g_xy_y_zz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_z_yz[i] = -2.0 * g_xy_y_0_yz[i] * a_exp + 4.0 * g_xy_y_zz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_z_zz[i] = -2.0 * g_xy_y_0_zz[i] * a_exp + 4.0 * g_xy_y_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (414-420)

    #pragma omp simd aligned(g_x_0_z_0_y_z_x_xx, g_x_0_z_0_y_z_x_xy, g_x_0_z_0_y_z_x_xz, g_x_0_z_0_y_z_x_yy, g_x_0_z_0_y_z_x_yz, g_x_0_z_0_y_z_x_zz, g_xy_z_xz_xx, g_xy_z_xz_xy, g_xy_z_xz_xz, g_xy_z_xz_yy, g_xy_z_xz_yz, g_xy_z_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_z_x_xx[i] = 4.0 * g_xy_z_xz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_x_xy[i] = 4.0 * g_xy_z_xz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_x_xz[i] = 4.0 * g_xy_z_xz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_x_yy[i] = 4.0 * g_xy_z_xz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_x_yz[i] = 4.0 * g_xy_z_xz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_x_zz[i] = 4.0 * g_xy_z_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (420-426)

    #pragma omp simd aligned(g_x_0_z_0_y_z_y_xx, g_x_0_z_0_y_z_y_xy, g_x_0_z_0_y_z_y_xz, g_x_0_z_0_y_z_y_yy, g_x_0_z_0_y_z_y_yz, g_x_0_z_0_y_z_y_zz, g_xy_z_yz_xx, g_xy_z_yz_xy, g_xy_z_yz_xz, g_xy_z_yz_yy, g_xy_z_yz_yz, g_xy_z_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_z_y_xx[i] = 4.0 * g_xy_z_yz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_y_xy[i] = 4.0 * g_xy_z_yz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_y_xz[i] = 4.0 * g_xy_z_yz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_y_yy[i] = 4.0 * g_xy_z_yz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_y_yz[i] = 4.0 * g_xy_z_yz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_y_zz[i] = 4.0 * g_xy_z_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (426-432)

    #pragma omp simd aligned(g_x_0_z_0_y_z_z_xx, g_x_0_z_0_y_z_z_xy, g_x_0_z_0_y_z_z_xz, g_x_0_z_0_y_z_z_yy, g_x_0_z_0_y_z_z_yz, g_x_0_z_0_y_z_z_zz, g_xy_z_0_xx, g_xy_z_0_xy, g_xy_z_0_xz, g_xy_z_0_yy, g_xy_z_0_yz, g_xy_z_0_zz, g_xy_z_zz_xx, g_xy_z_zz_xy, g_xy_z_zz_xz, g_xy_z_zz_yy, g_xy_z_zz_yz, g_xy_z_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_z_z_xx[i] = -2.0 * g_xy_z_0_xx[i] * a_exp + 4.0 * g_xy_z_zz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_z_xy[i] = -2.0 * g_xy_z_0_xy[i] * a_exp + 4.0 * g_xy_z_zz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_z_xz[i] = -2.0 * g_xy_z_0_xz[i] * a_exp + 4.0 * g_xy_z_zz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_z_yy[i] = -2.0 * g_xy_z_0_yy[i] * a_exp + 4.0 * g_xy_z_zz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_z_yz[i] = -2.0 * g_xy_z_0_yz[i] * a_exp + 4.0 * g_xy_z_zz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_z_zz[i] = -2.0 * g_xy_z_0_zz[i] * a_exp + 4.0 * g_xy_z_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (432-438)

    #pragma omp simd aligned(g_x_0_z_0_z_x_x_xx, g_x_0_z_0_z_x_x_xy, g_x_0_z_0_z_x_x_xz, g_x_0_z_0_z_x_x_yy, g_x_0_z_0_z_x_x_yz, g_x_0_z_0_z_x_x_zz, g_xz_x_xz_xx, g_xz_x_xz_xy, g_xz_x_xz_xz, g_xz_x_xz_yy, g_xz_x_xz_yz, g_xz_x_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_x_x_xx[i] = 4.0 * g_xz_x_xz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_x_xy[i] = 4.0 * g_xz_x_xz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_x_xz[i] = 4.0 * g_xz_x_xz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_x_yy[i] = 4.0 * g_xz_x_xz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_x_yz[i] = 4.0 * g_xz_x_xz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_x_zz[i] = 4.0 * g_xz_x_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (438-444)

    #pragma omp simd aligned(g_x_0_z_0_z_x_y_xx, g_x_0_z_0_z_x_y_xy, g_x_0_z_0_z_x_y_xz, g_x_0_z_0_z_x_y_yy, g_x_0_z_0_z_x_y_yz, g_x_0_z_0_z_x_y_zz, g_xz_x_yz_xx, g_xz_x_yz_xy, g_xz_x_yz_xz, g_xz_x_yz_yy, g_xz_x_yz_yz, g_xz_x_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_x_y_xx[i] = 4.0 * g_xz_x_yz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_y_xy[i] = 4.0 * g_xz_x_yz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_y_xz[i] = 4.0 * g_xz_x_yz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_y_yy[i] = 4.0 * g_xz_x_yz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_y_yz[i] = 4.0 * g_xz_x_yz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_y_zz[i] = 4.0 * g_xz_x_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (444-450)

    #pragma omp simd aligned(g_x_0_z_0_z_x_z_xx, g_x_0_z_0_z_x_z_xy, g_x_0_z_0_z_x_z_xz, g_x_0_z_0_z_x_z_yy, g_x_0_z_0_z_x_z_yz, g_x_0_z_0_z_x_z_zz, g_xz_x_0_xx, g_xz_x_0_xy, g_xz_x_0_xz, g_xz_x_0_yy, g_xz_x_0_yz, g_xz_x_0_zz, g_xz_x_zz_xx, g_xz_x_zz_xy, g_xz_x_zz_xz, g_xz_x_zz_yy, g_xz_x_zz_yz, g_xz_x_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_x_z_xx[i] = -2.0 * g_xz_x_0_xx[i] * a_exp + 4.0 * g_xz_x_zz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_z_xy[i] = -2.0 * g_xz_x_0_xy[i] * a_exp + 4.0 * g_xz_x_zz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_z_xz[i] = -2.0 * g_xz_x_0_xz[i] * a_exp + 4.0 * g_xz_x_zz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_z_yy[i] = -2.0 * g_xz_x_0_yy[i] * a_exp + 4.0 * g_xz_x_zz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_z_yz[i] = -2.0 * g_xz_x_0_yz[i] * a_exp + 4.0 * g_xz_x_zz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_z_zz[i] = -2.0 * g_xz_x_0_zz[i] * a_exp + 4.0 * g_xz_x_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (450-456)

    #pragma omp simd aligned(g_x_0_z_0_z_y_x_xx, g_x_0_z_0_z_y_x_xy, g_x_0_z_0_z_y_x_xz, g_x_0_z_0_z_y_x_yy, g_x_0_z_0_z_y_x_yz, g_x_0_z_0_z_y_x_zz, g_xz_y_xz_xx, g_xz_y_xz_xy, g_xz_y_xz_xz, g_xz_y_xz_yy, g_xz_y_xz_yz, g_xz_y_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_y_x_xx[i] = 4.0 * g_xz_y_xz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_x_xy[i] = 4.0 * g_xz_y_xz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_x_xz[i] = 4.0 * g_xz_y_xz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_x_yy[i] = 4.0 * g_xz_y_xz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_x_yz[i] = 4.0 * g_xz_y_xz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_x_zz[i] = 4.0 * g_xz_y_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (456-462)

    #pragma omp simd aligned(g_x_0_z_0_z_y_y_xx, g_x_0_z_0_z_y_y_xy, g_x_0_z_0_z_y_y_xz, g_x_0_z_0_z_y_y_yy, g_x_0_z_0_z_y_y_yz, g_x_0_z_0_z_y_y_zz, g_xz_y_yz_xx, g_xz_y_yz_xy, g_xz_y_yz_xz, g_xz_y_yz_yy, g_xz_y_yz_yz, g_xz_y_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_y_y_xx[i] = 4.0 * g_xz_y_yz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_y_xy[i] = 4.0 * g_xz_y_yz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_y_xz[i] = 4.0 * g_xz_y_yz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_y_yy[i] = 4.0 * g_xz_y_yz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_y_yz[i] = 4.0 * g_xz_y_yz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_y_zz[i] = 4.0 * g_xz_y_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (462-468)

    #pragma omp simd aligned(g_x_0_z_0_z_y_z_xx, g_x_0_z_0_z_y_z_xy, g_x_0_z_0_z_y_z_xz, g_x_0_z_0_z_y_z_yy, g_x_0_z_0_z_y_z_yz, g_x_0_z_0_z_y_z_zz, g_xz_y_0_xx, g_xz_y_0_xy, g_xz_y_0_xz, g_xz_y_0_yy, g_xz_y_0_yz, g_xz_y_0_zz, g_xz_y_zz_xx, g_xz_y_zz_xy, g_xz_y_zz_xz, g_xz_y_zz_yy, g_xz_y_zz_yz, g_xz_y_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_y_z_xx[i] = -2.0 * g_xz_y_0_xx[i] * a_exp + 4.0 * g_xz_y_zz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_z_xy[i] = -2.0 * g_xz_y_0_xy[i] * a_exp + 4.0 * g_xz_y_zz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_z_xz[i] = -2.0 * g_xz_y_0_xz[i] * a_exp + 4.0 * g_xz_y_zz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_z_yy[i] = -2.0 * g_xz_y_0_yy[i] * a_exp + 4.0 * g_xz_y_zz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_z_yz[i] = -2.0 * g_xz_y_0_yz[i] * a_exp + 4.0 * g_xz_y_zz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_z_zz[i] = -2.0 * g_xz_y_0_zz[i] * a_exp + 4.0 * g_xz_y_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (468-474)

    #pragma omp simd aligned(g_x_0_z_0_z_z_x_xx, g_x_0_z_0_z_z_x_xy, g_x_0_z_0_z_z_x_xz, g_x_0_z_0_z_z_x_yy, g_x_0_z_0_z_z_x_yz, g_x_0_z_0_z_z_x_zz, g_xz_z_xz_xx, g_xz_z_xz_xy, g_xz_z_xz_xz, g_xz_z_xz_yy, g_xz_z_xz_yz, g_xz_z_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_z_x_xx[i] = 4.0 * g_xz_z_xz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_x_xy[i] = 4.0 * g_xz_z_xz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_x_xz[i] = 4.0 * g_xz_z_xz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_x_yy[i] = 4.0 * g_xz_z_xz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_x_yz[i] = 4.0 * g_xz_z_xz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_x_zz[i] = 4.0 * g_xz_z_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (474-480)

    #pragma omp simd aligned(g_x_0_z_0_z_z_y_xx, g_x_0_z_0_z_z_y_xy, g_x_0_z_0_z_z_y_xz, g_x_0_z_0_z_z_y_yy, g_x_0_z_0_z_z_y_yz, g_x_0_z_0_z_z_y_zz, g_xz_z_yz_xx, g_xz_z_yz_xy, g_xz_z_yz_xz, g_xz_z_yz_yy, g_xz_z_yz_yz, g_xz_z_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_z_y_xx[i] = 4.0 * g_xz_z_yz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_y_xy[i] = 4.0 * g_xz_z_yz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_y_xz[i] = 4.0 * g_xz_z_yz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_y_yy[i] = 4.0 * g_xz_z_yz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_y_yz[i] = 4.0 * g_xz_z_yz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_y_zz[i] = 4.0 * g_xz_z_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (480-486)

    #pragma omp simd aligned(g_x_0_z_0_z_z_z_xx, g_x_0_z_0_z_z_z_xy, g_x_0_z_0_z_z_z_xz, g_x_0_z_0_z_z_z_yy, g_x_0_z_0_z_z_z_yz, g_x_0_z_0_z_z_z_zz, g_xz_z_0_xx, g_xz_z_0_xy, g_xz_z_0_xz, g_xz_z_0_yy, g_xz_z_0_yz, g_xz_z_0_zz, g_xz_z_zz_xx, g_xz_z_zz_xy, g_xz_z_zz_xz, g_xz_z_zz_yy, g_xz_z_zz_yz, g_xz_z_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_z_z_xx[i] = -2.0 * g_xz_z_0_xx[i] * a_exp + 4.0 * g_xz_z_zz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_z_xy[i] = -2.0 * g_xz_z_0_xy[i] * a_exp + 4.0 * g_xz_z_zz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_z_xz[i] = -2.0 * g_xz_z_0_xz[i] * a_exp + 4.0 * g_xz_z_zz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_z_yy[i] = -2.0 * g_xz_z_0_yy[i] * a_exp + 4.0 * g_xz_z_zz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_z_yz[i] = -2.0 * g_xz_z_0_yz[i] * a_exp + 4.0 * g_xz_z_zz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_z_zz[i] = -2.0 * g_xz_z_0_zz[i] * a_exp + 4.0 * g_xz_z_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (486-492)

    #pragma omp simd aligned(g_xy_x_0_xx, g_xy_x_0_xy, g_xy_x_0_xz, g_xy_x_0_yy, g_xy_x_0_yz, g_xy_x_0_zz, g_xy_x_xx_xx, g_xy_x_xx_xy, g_xy_x_xx_xz, g_xy_x_xx_yy, g_xy_x_xx_yz, g_xy_x_xx_zz, g_y_0_x_0_x_x_x_xx, g_y_0_x_0_x_x_x_xy, g_y_0_x_0_x_x_x_xz, g_y_0_x_0_x_x_x_yy, g_y_0_x_0_x_x_x_yz, g_y_0_x_0_x_x_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_x_x_xx[i] = -2.0 * g_xy_x_0_xx[i] * a_exp + 4.0 * g_xy_x_xx_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_x_xy[i] = -2.0 * g_xy_x_0_xy[i] * a_exp + 4.0 * g_xy_x_xx_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_x_xz[i] = -2.0 * g_xy_x_0_xz[i] * a_exp + 4.0 * g_xy_x_xx_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_x_yy[i] = -2.0 * g_xy_x_0_yy[i] * a_exp + 4.0 * g_xy_x_xx_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_x_yz[i] = -2.0 * g_xy_x_0_yz[i] * a_exp + 4.0 * g_xy_x_xx_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_x_zz[i] = -2.0 * g_xy_x_0_zz[i] * a_exp + 4.0 * g_xy_x_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (492-498)

    #pragma omp simd aligned(g_xy_x_xy_xx, g_xy_x_xy_xy, g_xy_x_xy_xz, g_xy_x_xy_yy, g_xy_x_xy_yz, g_xy_x_xy_zz, g_y_0_x_0_x_x_y_xx, g_y_0_x_0_x_x_y_xy, g_y_0_x_0_x_x_y_xz, g_y_0_x_0_x_x_y_yy, g_y_0_x_0_x_x_y_yz, g_y_0_x_0_x_x_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_x_y_xx[i] = 4.0 * g_xy_x_xy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_y_xy[i] = 4.0 * g_xy_x_xy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_y_xz[i] = 4.0 * g_xy_x_xy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_y_yy[i] = 4.0 * g_xy_x_xy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_y_yz[i] = 4.0 * g_xy_x_xy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_y_zz[i] = 4.0 * g_xy_x_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (498-504)

    #pragma omp simd aligned(g_xy_x_xz_xx, g_xy_x_xz_xy, g_xy_x_xz_xz, g_xy_x_xz_yy, g_xy_x_xz_yz, g_xy_x_xz_zz, g_y_0_x_0_x_x_z_xx, g_y_0_x_0_x_x_z_xy, g_y_0_x_0_x_x_z_xz, g_y_0_x_0_x_x_z_yy, g_y_0_x_0_x_x_z_yz, g_y_0_x_0_x_x_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_x_z_xx[i] = 4.0 * g_xy_x_xz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_z_xy[i] = 4.0 * g_xy_x_xz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_z_xz[i] = 4.0 * g_xy_x_xz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_z_yy[i] = 4.0 * g_xy_x_xz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_z_yz[i] = 4.0 * g_xy_x_xz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_z_zz[i] = 4.0 * g_xy_x_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (504-510)

    #pragma omp simd aligned(g_xy_y_0_xx, g_xy_y_0_xy, g_xy_y_0_xz, g_xy_y_0_yy, g_xy_y_0_yz, g_xy_y_0_zz, g_xy_y_xx_xx, g_xy_y_xx_xy, g_xy_y_xx_xz, g_xy_y_xx_yy, g_xy_y_xx_yz, g_xy_y_xx_zz, g_y_0_x_0_x_y_x_xx, g_y_0_x_0_x_y_x_xy, g_y_0_x_0_x_y_x_xz, g_y_0_x_0_x_y_x_yy, g_y_0_x_0_x_y_x_yz, g_y_0_x_0_x_y_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_y_x_xx[i] = -2.0 * g_xy_y_0_xx[i] * a_exp + 4.0 * g_xy_y_xx_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_x_xy[i] = -2.0 * g_xy_y_0_xy[i] * a_exp + 4.0 * g_xy_y_xx_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_x_xz[i] = -2.0 * g_xy_y_0_xz[i] * a_exp + 4.0 * g_xy_y_xx_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_x_yy[i] = -2.0 * g_xy_y_0_yy[i] * a_exp + 4.0 * g_xy_y_xx_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_x_yz[i] = -2.0 * g_xy_y_0_yz[i] * a_exp + 4.0 * g_xy_y_xx_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_x_zz[i] = -2.0 * g_xy_y_0_zz[i] * a_exp + 4.0 * g_xy_y_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (510-516)

    #pragma omp simd aligned(g_xy_y_xy_xx, g_xy_y_xy_xy, g_xy_y_xy_xz, g_xy_y_xy_yy, g_xy_y_xy_yz, g_xy_y_xy_zz, g_y_0_x_0_x_y_y_xx, g_y_0_x_0_x_y_y_xy, g_y_0_x_0_x_y_y_xz, g_y_0_x_0_x_y_y_yy, g_y_0_x_0_x_y_y_yz, g_y_0_x_0_x_y_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_y_y_xx[i] = 4.0 * g_xy_y_xy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_y_xy[i] = 4.0 * g_xy_y_xy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_y_xz[i] = 4.0 * g_xy_y_xy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_y_yy[i] = 4.0 * g_xy_y_xy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_y_yz[i] = 4.0 * g_xy_y_xy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_y_zz[i] = 4.0 * g_xy_y_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (516-522)

    #pragma omp simd aligned(g_xy_y_xz_xx, g_xy_y_xz_xy, g_xy_y_xz_xz, g_xy_y_xz_yy, g_xy_y_xz_yz, g_xy_y_xz_zz, g_y_0_x_0_x_y_z_xx, g_y_0_x_0_x_y_z_xy, g_y_0_x_0_x_y_z_xz, g_y_0_x_0_x_y_z_yy, g_y_0_x_0_x_y_z_yz, g_y_0_x_0_x_y_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_y_z_xx[i] = 4.0 * g_xy_y_xz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_z_xy[i] = 4.0 * g_xy_y_xz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_z_xz[i] = 4.0 * g_xy_y_xz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_z_yy[i] = 4.0 * g_xy_y_xz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_z_yz[i] = 4.0 * g_xy_y_xz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_z_zz[i] = 4.0 * g_xy_y_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (522-528)

    #pragma omp simd aligned(g_xy_z_0_xx, g_xy_z_0_xy, g_xy_z_0_xz, g_xy_z_0_yy, g_xy_z_0_yz, g_xy_z_0_zz, g_xy_z_xx_xx, g_xy_z_xx_xy, g_xy_z_xx_xz, g_xy_z_xx_yy, g_xy_z_xx_yz, g_xy_z_xx_zz, g_y_0_x_0_x_z_x_xx, g_y_0_x_0_x_z_x_xy, g_y_0_x_0_x_z_x_xz, g_y_0_x_0_x_z_x_yy, g_y_0_x_0_x_z_x_yz, g_y_0_x_0_x_z_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_z_x_xx[i] = -2.0 * g_xy_z_0_xx[i] * a_exp + 4.0 * g_xy_z_xx_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_x_xy[i] = -2.0 * g_xy_z_0_xy[i] * a_exp + 4.0 * g_xy_z_xx_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_x_xz[i] = -2.0 * g_xy_z_0_xz[i] * a_exp + 4.0 * g_xy_z_xx_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_x_yy[i] = -2.0 * g_xy_z_0_yy[i] * a_exp + 4.0 * g_xy_z_xx_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_x_yz[i] = -2.0 * g_xy_z_0_yz[i] * a_exp + 4.0 * g_xy_z_xx_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_x_zz[i] = -2.0 * g_xy_z_0_zz[i] * a_exp + 4.0 * g_xy_z_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (528-534)

    #pragma omp simd aligned(g_xy_z_xy_xx, g_xy_z_xy_xy, g_xy_z_xy_xz, g_xy_z_xy_yy, g_xy_z_xy_yz, g_xy_z_xy_zz, g_y_0_x_0_x_z_y_xx, g_y_0_x_0_x_z_y_xy, g_y_0_x_0_x_z_y_xz, g_y_0_x_0_x_z_y_yy, g_y_0_x_0_x_z_y_yz, g_y_0_x_0_x_z_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_z_y_xx[i] = 4.0 * g_xy_z_xy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_y_xy[i] = 4.0 * g_xy_z_xy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_y_xz[i] = 4.0 * g_xy_z_xy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_y_yy[i] = 4.0 * g_xy_z_xy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_y_yz[i] = 4.0 * g_xy_z_xy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_y_zz[i] = 4.0 * g_xy_z_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (534-540)

    #pragma omp simd aligned(g_xy_z_xz_xx, g_xy_z_xz_xy, g_xy_z_xz_xz, g_xy_z_xz_yy, g_xy_z_xz_yz, g_xy_z_xz_zz, g_y_0_x_0_x_z_z_xx, g_y_0_x_0_x_z_z_xy, g_y_0_x_0_x_z_z_xz, g_y_0_x_0_x_z_z_yy, g_y_0_x_0_x_z_z_yz, g_y_0_x_0_x_z_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_z_z_xx[i] = 4.0 * g_xy_z_xz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_z_xy[i] = 4.0 * g_xy_z_xz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_z_xz[i] = 4.0 * g_xy_z_xz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_z_yy[i] = 4.0 * g_xy_z_xz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_z_yz[i] = 4.0 * g_xy_z_xz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_z_zz[i] = 4.0 * g_xy_z_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (540-546)

    #pragma omp simd aligned(g_0_x_0_xx, g_0_x_0_xy, g_0_x_0_xz, g_0_x_0_yy, g_0_x_0_yz, g_0_x_0_zz, g_0_x_xx_xx, g_0_x_xx_xy, g_0_x_xx_xz, g_0_x_xx_yy, g_0_x_xx_yz, g_0_x_xx_zz, g_y_0_x_0_y_x_x_xx, g_y_0_x_0_y_x_x_xy, g_y_0_x_0_y_x_x_xz, g_y_0_x_0_y_x_x_yy, g_y_0_x_0_y_x_x_yz, g_y_0_x_0_y_x_x_zz, g_yy_x_0_xx, g_yy_x_0_xy, g_yy_x_0_xz, g_yy_x_0_yy, g_yy_x_0_yz, g_yy_x_0_zz, g_yy_x_xx_xx, g_yy_x_xx_xy, g_yy_x_xx_xz, g_yy_x_xx_yy, g_yy_x_xx_yz, g_yy_x_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_x_x_xx[i] = g_0_x_0_xx[i] - 2.0 * g_0_x_xx_xx[i] * c_exps[i] - 2.0 * g_yy_x_0_xx[i] * a_exp + 4.0 * g_yy_x_xx_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_x_xy[i] = g_0_x_0_xy[i] - 2.0 * g_0_x_xx_xy[i] * c_exps[i] - 2.0 * g_yy_x_0_xy[i] * a_exp + 4.0 * g_yy_x_xx_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_x_xz[i] = g_0_x_0_xz[i] - 2.0 * g_0_x_xx_xz[i] * c_exps[i] - 2.0 * g_yy_x_0_xz[i] * a_exp + 4.0 * g_yy_x_xx_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_x_yy[i] = g_0_x_0_yy[i] - 2.0 * g_0_x_xx_yy[i] * c_exps[i] - 2.0 * g_yy_x_0_yy[i] * a_exp + 4.0 * g_yy_x_xx_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_x_yz[i] = g_0_x_0_yz[i] - 2.0 * g_0_x_xx_yz[i] * c_exps[i] - 2.0 * g_yy_x_0_yz[i] * a_exp + 4.0 * g_yy_x_xx_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_x_zz[i] = g_0_x_0_zz[i] - 2.0 * g_0_x_xx_zz[i] * c_exps[i] - 2.0 * g_yy_x_0_zz[i] * a_exp + 4.0 * g_yy_x_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (546-552)

    #pragma omp simd aligned(g_0_x_xy_xx, g_0_x_xy_xy, g_0_x_xy_xz, g_0_x_xy_yy, g_0_x_xy_yz, g_0_x_xy_zz, g_y_0_x_0_y_x_y_xx, g_y_0_x_0_y_x_y_xy, g_y_0_x_0_y_x_y_xz, g_y_0_x_0_y_x_y_yy, g_y_0_x_0_y_x_y_yz, g_y_0_x_0_y_x_y_zz, g_yy_x_xy_xx, g_yy_x_xy_xy, g_yy_x_xy_xz, g_yy_x_xy_yy, g_yy_x_xy_yz, g_yy_x_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_x_y_xx[i] = -2.0 * g_0_x_xy_xx[i] * c_exps[i] + 4.0 * g_yy_x_xy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_y_xy[i] = -2.0 * g_0_x_xy_xy[i] * c_exps[i] + 4.0 * g_yy_x_xy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_y_xz[i] = -2.0 * g_0_x_xy_xz[i] * c_exps[i] + 4.0 * g_yy_x_xy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_y_yy[i] = -2.0 * g_0_x_xy_yy[i] * c_exps[i] + 4.0 * g_yy_x_xy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_y_yz[i] = -2.0 * g_0_x_xy_yz[i] * c_exps[i] + 4.0 * g_yy_x_xy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_y_zz[i] = -2.0 * g_0_x_xy_zz[i] * c_exps[i] + 4.0 * g_yy_x_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (552-558)

    #pragma omp simd aligned(g_0_x_xz_xx, g_0_x_xz_xy, g_0_x_xz_xz, g_0_x_xz_yy, g_0_x_xz_yz, g_0_x_xz_zz, g_y_0_x_0_y_x_z_xx, g_y_0_x_0_y_x_z_xy, g_y_0_x_0_y_x_z_xz, g_y_0_x_0_y_x_z_yy, g_y_0_x_0_y_x_z_yz, g_y_0_x_0_y_x_z_zz, g_yy_x_xz_xx, g_yy_x_xz_xy, g_yy_x_xz_xz, g_yy_x_xz_yy, g_yy_x_xz_yz, g_yy_x_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_x_z_xx[i] = -2.0 * g_0_x_xz_xx[i] * c_exps[i] + 4.0 * g_yy_x_xz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_z_xy[i] = -2.0 * g_0_x_xz_xy[i] * c_exps[i] + 4.0 * g_yy_x_xz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_z_xz[i] = -2.0 * g_0_x_xz_xz[i] * c_exps[i] + 4.0 * g_yy_x_xz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_z_yy[i] = -2.0 * g_0_x_xz_yy[i] * c_exps[i] + 4.0 * g_yy_x_xz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_z_yz[i] = -2.0 * g_0_x_xz_yz[i] * c_exps[i] + 4.0 * g_yy_x_xz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_z_zz[i] = -2.0 * g_0_x_xz_zz[i] * c_exps[i] + 4.0 * g_yy_x_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (558-564)

    #pragma omp simd aligned(g_0_y_0_xx, g_0_y_0_xy, g_0_y_0_xz, g_0_y_0_yy, g_0_y_0_yz, g_0_y_0_zz, g_0_y_xx_xx, g_0_y_xx_xy, g_0_y_xx_xz, g_0_y_xx_yy, g_0_y_xx_yz, g_0_y_xx_zz, g_y_0_x_0_y_y_x_xx, g_y_0_x_0_y_y_x_xy, g_y_0_x_0_y_y_x_xz, g_y_0_x_0_y_y_x_yy, g_y_0_x_0_y_y_x_yz, g_y_0_x_0_y_y_x_zz, g_yy_y_0_xx, g_yy_y_0_xy, g_yy_y_0_xz, g_yy_y_0_yy, g_yy_y_0_yz, g_yy_y_0_zz, g_yy_y_xx_xx, g_yy_y_xx_xy, g_yy_y_xx_xz, g_yy_y_xx_yy, g_yy_y_xx_yz, g_yy_y_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_y_x_xx[i] = g_0_y_0_xx[i] - 2.0 * g_0_y_xx_xx[i] * c_exps[i] - 2.0 * g_yy_y_0_xx[i] * a_exp + 4.0 * g_yy_y_xx_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_x_xy[i] = g_0_y_0_xy[i] - 2.0 * g_0_y_xx_xy[i] * c_exps[i] - 2.0 * g_yy_y_0_xy[i] * a_exp + 4.0 * g_yy_y_xx_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_x_xz[i] = g_0_y_0_xz[i] - 2.0 * g_0_y_xx_xz[i] * c_exps[i] - 2.0 * g_yy_y_0_xz[i] * a_exp + 4.0 * g_yy_y_xx_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_x_yy[i] = g_0_y_0_yy[i] - 2.0 * g_0_y_xx_yy[i] * c_exps[i] - 2.0 * g_yy_y_0_yy[i] * a_exp + 4.0 * g_yy_y_xx_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_x_yz[i] = g_0_y_0_yz[i] - 2.0 * g_0_y_xx_yz[i] * c_exps[i] - 2.0 * g_yy_y_0_yz[i] * a_exp + 4.0 * g_yy_y_xx_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_x_zz[i] = g_0_y_0_zz[i] - 2.0 * g_0_y_xx_zz[i] * c_exps[i] - 2.0 * g_yy_y_0_zz[i] * a_exp + 4.0 * g_yy_y_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (564-570)

    #pragma omp simd aligned(g_0_y_xy_xx, g_0_y_xy_xy, g_0_y_xy_xz, g_0_y_xy_yy, g_0_y_xy_yz, g_0_y_xy_zz, g_y_0_x_0_y_y_y_xx, g_y_0_x_0_y_y_y_xy, g_y_0_x_0_y_y_y_xz, g_y_0_x_0_y_y_y_yy, g_y_0_x_0_y_y_y_yz, g_y_0_x_0_y_y_y_zz, g_yy_y_xy_xx, g_yy_y_xy_xy, g_yy_y_xy_xz, g_yy_y_xy_yy, g_yy_y_xy_yz, g_yy_y_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_y_y_xx[i] = -2.0 * g_0_y_xy_xx[i] * c_exps[i] + 4.0 * g_yy_y_xy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_y_xy[i] = -2.0 * g_0_y_xy_xy[i] * c_exps[i] + 4.0 * g_yy_y_xy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_y_xz[i] = -2.0 * g_0_y_xy_xz[i] * c_exps[i] + 4.0 * g_yy_y_xy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_y_yy[i] = -2.0 * g_0_y_xy_yy[i] * c_exps[i] + 4.0 * g_yy_y_xy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_y_yz[i] = -2.0 * g_0_y_xy_yz[i] * c_exps[i] + 4.0 * g_yy_y_xy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_y_zz[i] = -2.0 * g_0_y_xy_zz[i] * c_exps[i] + 4.0 * g_yy_y_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (570-576)

    #pragma omp simd aligned(g_0_y_xz_xx, g_0_y_xz_xy, g_0_y_xz_xz, g_0_y_xz_yy, g_0_y_xz_yz, g_0_y_xz_zz, g_y_0_x_0_y_y_z_xx, g_y_0_x_0_y_y_z_xy, g_y_0_x_0_y_y_z_xz, g_y_0_x_0_y_y_z_yy, g_y_0_x_0_y_y_z_yz, g_y_0_x_0_y_y_z_zz, g_yy_y_xz_xx, g_yy_y_xz_xy, g_yy_y_xz_xz, g_yy_y_xz_yy, g_yy_y_xz_yz, g_yy_y_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_y_z_xx[i] = -2.0 * g_0_y_xz_xx[i] * c_exps[i] + 4.0 * g_yy_y_xz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_z_xy[i] = -2.0 * g_0_y_xz_xy[i] * c_exps[i] + 4.0 * g_yy_y_xz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_z_xz[i] = -2.0 * g_0_y_xz_xz[i] * c_exps[i] + 4.0 * g_yy_y_xz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_z_yy[i] = -2.0 * g_0_y_xz_yy[i] * c_exps[i] + 4.0 * g_yy_y_xz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_z_yz[i] = -2.0 * g_0_y_xz_yz[i] * c_exps[i] + 4.0 * g_yy_y_xz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_z_zz[i] = -2.0 * g_0_y_xz_zz[i] * c_exps[i] + 4.0 * g_yy_y_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (576-582)

    #pragma omp simd aligned(g_0_z_0_xx, g_0_z_0_xy, g_0_z_0_xz, g_0_z_0_yy, g_0_z_0_yz, g_0_z_0_zz, g_0_z_xx_xx, g_0_z_xx_xy, g_0_z_xx_xz, g_0_z_xx_yy, g_0_z_xx_yz, g_0_z_xx_zz, g_y_0_x_0_y_z_x_xx, g_y_0_x_0_y_z_x_xy, g_y_0_x_0_y_z_x_xz, g_y_0_x_0_y_z_x_yy, g_y_0_x_0_y_z_x_yz, g_y_0_x_0_y_z_x_zz, g_yy_z_0_xx, g_yy_z_0_xy, g_yy_z_0_xz, g_yy_z_0_yy, g_yy_z_0_yz, g_yy_z_0_zz, g_yy_z_xx_xx, g_yy_z_xx_xy, g_yy_z_xx_xz, g_yy_z_xx_yy, g_yy_z_xx_yz, g_yy_z_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_z_x_xx[i] = g_0_z_0_xx[i] - 2.0 * g_0_z_xx_xx[i] * c_exps[i] - 2.0 * g_yy_z_0_xx[i] * a_exp + 4.0 * g_yy_z_xx_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_x_xy[i] = g_0_z_0_xy[i] - 2.0 * g_0_z_xx_xy[i] * c_exps[i] - 2.0 * g_yy_z_0_xy[i] * a_exp + 4.0 * g_yy_z_xx_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_x_xz[i] = g_0_z_0_xz[i] - 2.0 * g_0_z_xx_xz[i] * c_exps[i] - 2.0 * g_yy_z_0_xz[i] * a_exp + 4.0 * g_yy_z_xx_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_x_yy[i] = g_0_z_0_yy[i] - 2.0 * g_0_z_xx_yy[i] * c_exps[i] - 2.0 * g_yy_z_0_yy[i] * a_exp + 4.0 * g_yy_z_xx_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_x_yz[i] = g_0_z_0_yz[i] - 2.0 * g_0_z_xx_yz[i] * c_exps[i] - 2.0 * g_yy_z_0_yz[i] * a_exp + 4.0 * g_yy_z_xx_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_x_zz[i] = g_0_z_0_zz[i] - 2.0 * g_0_z_xx_zz[i] * c_exps[i] - 2.0 * g_yy_z_0_zz[i] * a_exp + 4.0 * g_yy_z_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (582-588)

    #pragma omp simd aligned(g_0_z_xy_xx, g_0_z_xy_xy, g_0_z_xy_xz, g_0_z_xy_yy, g_0_z_xy_yz, g_0_z_xy_zz, g_y_0_x_0_y_z_y_xx, g_y_0_x_0_y_z_y_xy, g_y_0_x_0_y_z_y_xz, g_y_0_x_0_y_z_y_yy, g_y_0_x_0_y_z_y_yz, g_y_0_x_0_y_z_y_zz, g_yy_z_xy_xx, g_yy_z_xy_xy, g_yy_z_xy_xz, g_yy_z_xy_yy, g_yy_z_xy_yz, g_yy_z_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_z_y_xx[i] = -2.0 * g_0_z_xy_xx[i] * c_exps[i] + 4.0 * g_yy_z_xy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_y_xy[i] = -2.0 * g_0_z_xy_xy[i] * c_exps[i] + 4.0 * g_yy_z_xy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_y_xz[i] = -2.0 * g_0_z_xy_xz[i] * c_exps[i] + 4.0 * g_yy_z_xy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_y_yy[i] = -2.0 * g_0_z_xy_yy[i] * c_exps[i] + 4.0 * g_yy_z_xy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_y_yz[i] = -2.0 * g_0_z_xy_yz[i] * c_exps[i] + 4.0 * g_yy_z_xy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_y_zz[i] = -2.0 * g_0_z_xy_zz[i] * c_exps[i] + 4.0 * g_yy_z_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (588-594)

    #pragma omp simd aligned(g_0_z_xz_xx, g_0_z_xz_xy, g_0_z_xz_xz, g_0_z_xz_yy, g_0_z_xz_yz, g_0_z_xz_zz, g_y_0_x_0_y_z_z_xx, g_y_0_x_0_y_z_z_xy, g_y_0_x_0_y_z_z_xz, g_y_0_x_0_y_z_z_yy, g_y_0_x_0_y_z_z_yz, g_y_0_x_0_y_z_z_zz, g_yy_z_xz_xx, g_yy_z_xz_xy, g_yy_z_xz_xz, g_yy_z_xz_yy, g_yy_z_xz_yz, g_yy_z_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_z_z_xx[i] = -2.0 * g_0_z_xz_xx[i] * c_exps[i] + 4.0 * g_yy_z_xz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_z_xy[i] = -2.0 * g_0_z_xz_xy[i] * c_exps[i] + 4.0 * g_yy_z_xz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_z_xz[i] = -2.0 * g_0_z_xz_xz[i] * c_exps[i] + 4.0 * g_yy_z_xz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_z_yy[i] = -2.0 * g_0_z_xz_yy[i] * c_exps[i] + 4.0 * g_yy_z_xz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_z_yz[i] = -2.0 * g_0_z_xz_yz[i] * c_exps[i] + 4.0 * g_yy_z_xz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_z_zz[i] = -2.0 * g_0_z_xz_zz[i] * c_exps[i] + 4.0 * g_yy_z_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (594-600)

    #pragma omp simd aligned(g_y_0_x_0_z_x_x_xx, g_y_0_x_0_z_x_x_xy, g_y_0_x_0_z_x_x_xz, g_y_0_x_0_z_x_x_yy, g_y_0_x_0_z_x_x_yz, g_y_0_x_0_z_x_x_zz, g_yz_x_0_xx, g_yz_x_0_xy, g_yz_x_0_xz, g_yz_x_0_yy, g_yz_x_0_yz, g_yz_x_0_zz, g_yz_x_xx_xx, g_yz_x_xx_xy, g_yz_x_xx_xz, g_yz_x_xx_yy, g_yz_x_xx_yz, g_yz_x_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_x_x_xx[i] = -2.0 * g_yz_x_0_xx[i] * a_exp + 4.0 * g_yz_x_xx_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_x_xy[i] = -2.0 * g_yz_x_0_xy[i] * a_exp + 4.0 * g_yz_x_xx_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_x_xz[i] = -2.0 * g_yz_x_0_xz[i] * a_exp + 4.0 * g_yz_x_xx_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_x_yy[i] = -2.0 * g_yz_x_0_yy[i] * a_exp + 4.0 * g_yz_x_xx_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_x_yz[i] = -2.0 * g_yz_x_0_yz[i] * a_exp + 4.0 * g_yz_x_xx_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_x_zz[i] = -2.0 * g_yz_x_0_zz[i] * a_exp + 4.0 * g_yz_x_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (600-606)

    #pragma omp simd aligned(g_y_0_x_0_z_x_y_xx, g_y_0_x_0_z_x_y_xy, g_y_0_x_0_z_x_y_xz, g_y_0_x_0_z_x_y_yy, g_y_0_x_0_z_x_y_yz, g_y_0_x_0_z_x_y_zz, g_yz_x_xy_xx, g_yz_x_xy_xy, g_yz_x_xy_xz, g_yz_x_xy_yy, g_yz_x_xy_yz, g_yz_x_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_x_y_xx[i] = 4.0 * g_yz_x_xy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_y_xy[i] = 4.0 * g_yz_x_xy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_y_xz[i] = 4.0 * g_yz_x_xy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_y_yy[i] = 4.0 * g_yz_x_xy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_y_yz[i] = 4.0 * g_yz_x_xy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_y_zz[i] = 4.0 * g_yz_x_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (606-612)

    #pragma omp simd aligned(g_y_0_x_0_z_x_z_xx, g_y_0_x_0_z_x_z_xy, g_y_0_x_0_z_x_z_xz, g_y_0_x_0_z_x_z_yy, g_y_0_x_0_z_x_z_yz, g_y_0_x_0_z_x_z_zz, g_yz_x_xz_xx, g_yz_x_xz_xy, g_yz_x_xz_xz, g_yz_x_xz_yy, g_yz_x_xz_yz, g_yz_x_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_x_z_xx[i] = 4.0 * g_yz_x_xz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_z_xy[i] = 4.0 * g_yz_x_xz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_z_xz[i] = 4.0 * g_yz_x_xz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_z_yy[i] = 4.0 * g_yz_x_xz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_z_yz[i] = 4.0 * g_yz_x_xz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_z_zz[i] = 4.0 * g_yz_x_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (612-618)

    #pragma omp simd aligned(g_y_0_x_0_z_y_x_xx, g_y_0_x_0_z_y_x_xy, g_y_0_x_0_z_y_x_xz, g_y_0_x_0_z_y_x_yy, g_y_0_x_0_z_y_x_yz, g_y_0_x_0_z_y_x_zz, g_yz_y_0_xx, g_yz_y_0_xy, g_yz_y_0_xz, g_yz_y_0_yy, g_yz_y_0_yz, g_yz_y_0_zz, g_yz_y_xx_xx, g_yz_y_xx_xy, g_yz_y_xx_xz, g_yz_y_xx_yy, g_yz_y_xx_yz, g_yz_y_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_y_x_xx[i] = -2.0 * g_yz_y_0_xx[i] * a_exp + 4.0 * g_yz_y_xx_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_x_xy[i] = -2.0 * g_yz_y_0_xy[i] * a_exp + 4.0 * g_yz_y_xx_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_x_xz[i] = -2.0 * g_yz_y_0_xz[i] * a_exp + 4.0 * g_yz_y_xx_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_x_yy[i] = -2.0 * g_yz_y_0_yy[i] * a_exp + 4.0 * g_yz_y_xx_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_x_yz[i] = -2.0 * g_yz_y_0_yz[i] * a_exp + 4.0 * g_yz_y_xx_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_x_zz[i] = -2.0 * g_yz_y_0_zz[i] * a_exp + 4.0 * g_yz_y_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (618-624)

    #pragma omp simd aligned(g_y_0_x_0_z_y_y_xx, g_y_0_x_0_z_y_y_xy, g_y_0_x_0_z_y_y_xz, g_y_0_x_0_z_y_y_yy, g_y_0_x_0_z_y_y_yz, g_y_0_x_0_z_y_y_zz, g_yz_y_xy_xx, g_yz_y_xy_xy, g_yz_y_xy_xz, g_yz_y_xy_yy, g_yz_y_xy_yz, g_yz_y_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_y_y_xx[i] = 4.0 * g_yz_y_xy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_y_xy[i] = 4.0 * g_yz_y_xy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_y_xz[i] = 4.0 * g_yz_y_xy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_y_yy[i] = 4.0 * g_yz_y_xy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_y_yz[i] = 4.0 * g_yz_y_xy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_y_zz[i] = 4.0 * g_yz_y_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (624-630)

    #pragma omp simd aligned(g_y_0_x_0_z_y_z_xx, g_y_0_x_0_z_y_z_xy, g_y_0_x_0_z_y_z_xz, g_y_0_x_0_z_y_z_yy, g_y_0_x_0_z_y_z_yz, g_y_0_x_0_z_y_z_zz, g_yz_y_xz_xx, g_yz_y_xz_xy, g_yz_y_xz_xz, g_yz_y_xz_yy, g_yz_y_xz_yz, g_yz_y_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_y_z_xx[i] = 4.0 * g_yz_y_xz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_z_xy[i] = 4.0 * g_yz_y_xz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_z_xz[i] = 4.0 * g_yz_y_xz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_z_yy[i] = 4.0 * g_yz_y_xz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_z_yz[i] = 4.0 * g_yz_y_xz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_z_zz[i] = 4.0 * g_yz_y_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (630-636)

    #pragma omp simd aligned(g_y_0_x_0_z_z_x_xx, g_y_0_x_0_z_z_x_xy, g_y_0_x_0_z_z_x_xz, g_y_0_x_0_z_z_x_yy, g_y_0_x_0_z_z_x_yz, g_y_0_x_0_z_z_x_zz, g_yz_z_0_xx, g_yz_z_0_xy, g_yz_z_0_xz, g_yz_z_0_yy, g_yz_z_0_yz, g_yz_z_0_zz, g_yz_z_xx_xx, g_yz_z_xx_xy, g_yz_z_xx_xz, g_yz_z_xx_yy, g_yz_z_xx_yz, g_yz_z_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_z_x_xx[i] = -2.0 * g_yz_z_0_xx[i] * a_exp + 4.0 * g_yz_z_xx_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_x_xy[i] = -2.0 * g_yz_z_0_xy[i] * a_exp + 4.0 * g_yz_z_xx_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_x_xz[i] = -2.0 * g_yz_z_0_xz[i] * a_exp + 4.0 * g_yz_z_xx_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_x_yy[i] = -2.0 * g_yz_z_0_yy[i] * a_exp + 4.0 * g_yz_z_xx_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_x_yz[i] = -2.0 * g_yz_z_0_yz[i] * a_exp + 4.0 * g_yz_z_xx_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_x_zz[i] = -2.0 * g_yz_z_0_zz[i] * a_exp + 4.0 * g_yz_z_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (636-642)

    #pragma omp simd aligned(g_y_0_x_0_z_z_y_xx, g_y_0_x_0_z_z_y_xy, g_y_0_x_0_z_z_y_xz, g_y_0_x_0_z_z_y_yy, g_y_0_x_0_z_z_y_yz, g_y_0_x_0_z_z_y_zz, g_yz_z_xy_xx, g_yz_z_xy_xy, g_yz_z_xy_xz, g_yz_z_xy_yy, g_yz_z_xy_yz, g_yz_z_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_z_y_xx[i] = 4.0 * g_yz_z_xy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_y_xy[i] = 4.0 * g_yz_z_xy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_y_xz[i] = 4.0 * g_yz_z_xy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_y_yy[i] = 4.0 * g_yz_z_xy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_y_yz[i] = 4.0 * g_yz_z_xy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_y_zz[i] = 4.0 * g_yz_z_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (642-648)

    #pragma omp simd aligned(g_y_0_x_0_z_z_z_xx, g_y_0_x_0_z_z_z_xy, g_y_0_x_0_z_z_z_xz, g_y_0_x_0_z_z_z_yy, g_y_0_x_0_z_z_z_yz, g_y_0_x_0_z_z_z_zz, g_yz_z_xz_xx, g_yz_z_xz_xy, g_yz_z_xz_xz, g_yz_z_xz_yy, g_yz_z_xz_yz, g_yz_z_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_z_z_xx[i] = 4.0 * g_yz_z_xz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_z_xy[i] = 4.0 * g_yz_z_xz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_z_xz[i] = 4.0 * g_yz_z_xz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_z_yy[i] = 4.0 * g_yz_z_xz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_z_yz[i] = 4.0 * g_yz_z_xz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_z_zz[i] = 4.0 * g_yz_z_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (648-654)

    #pragma omp simd aligned(g_xy_x_xy_xx, g_xy_x_xy_xy, g_xy_x_xy_xz, g_xy_x_xy_yy, g_xy_x_xy_yz, g_xy_x_xy_zz, g_y_0_y_0_x_x_x_xx, g_y_0_y_0_x_x_x_xy, g_y_0_y_0_x_x_x_xz, g_y_0_y_0_x_x_x_yy, g_y_0_y_0_x_x_x_yz, g_y_0_y_0_x_x_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_x_x_xx[i] = 4.0 * g_xy_x_xy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_x_xy[i] = 4.0 * g_xy_x_xy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_x_xz[i] = 4.0 * g_xy_x_xy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_x_yy[i] = 4.0 * g_xy_x_xy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_x_yz[i] = 4.0 * g_xy_x_xy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_x_zz[i] = 4.0 * g_xy_x_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (654-660)

    #pragma omp simd aligned(g_xy_x_0_xx, g_xy_x_0_xy, g_xy_x_0_xz, g_xy_x_0_yy, g_xy_x_0_yz, g_xy_x_0_zz, g_xy_x_yy_xx, g_xy_x_yy_xy, g_xy_x_yy_xz, g_xy_x_yy_yy, g_xy_x_yy_yz, g_xy_x_yy_zz, g_y_0_y_0_x_x_y_xx, g_y_0_y_0_x_x_y_xy, g_y_0_y_0_x_x_y_xz, g_y_0_y_0_x_x_y_yy, g_y_0_y_0_x_x_y_yz, g_y_0_y_0_x_x_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_x_y_xx[i] = -2.0 * g_xy_x_0_xx[i] * a_exp + 4.0 * g_xy_x_yy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_y_xy[i] = -2.0 * g_xy_x_0_xy[i] * a_exp + 4.0 * g_xy_x_yy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_y_xz[i] = -2.0 * g_xy_x_0_xz[i] * a_exp + 4.0 * g_xy_x_yy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_y_yy[i] = -2.0 * g_xy_x_0_yy[i] * a_exp + 4.0 * g_xy_x_yy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_y_yz[i] = -2.0 * g_xy_x_0_yz[i] * a_exp + 4.0 * g_xy_x_yy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_y_zz[i] = -2.0 * g_xy_x_0_zz[i] * a_exp + 4.0 * g_xy_x_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (660-666)

    #pragma omp simd aligned(g_xy_x_yz_xx, g_xy_x_yz_xy, g_xy_x_yz_xz, g_xy_x_yz_yy, g_xy_x_yz_yz, g_xy_x_yz_zz, g_y_0_y_0_x_x_z_xx, g_y_0_y_0_x_x_z_xy, g_y_0_y_0_x_x_z_xz, g_y_0_y_0_x_x_z_yy, g_y_0_y_0_x_x_z_yz, g_y_0_y_0_x_x_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_x_z_xx[i] = 4.0 * g_xy_x_yz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_z_xy[i] = 4.0 * g_xy_x_yz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_z_xz[i] = 4.0 * g_xy_x_yz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_z_yy[i] = 4.0 * g_xy_x_yz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_z_yz[i] = 4.0 * g_xy_x_yz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_z_zz[i] = 4.0 * g_xy_x_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (666-672)

    #pragma omp simd aligned(g_xy_y_xy_xx, g_xy_y_xy_xy, g_xy_y_xy_xz, g_xy_y_xy_yy, g_xy_y_xy_yz, g_xy_y_xy_zz, g_y_0_y_0_x_y_x_xx, g_y_0_y_0_x_y_x_xy, g_y_0_y_0_x_y_x_xz, g_y_0_y_0_x_y_x_yy, g_y_0_y_0_x_y_x_yz, g_y_0_y_0_x_y_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_y_x_xx[i] = 4.0 * g_xy_y_xy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_x_xy[i] = 4.0 * g_xy_y_xy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_x_xz[i] = 4.0 * g_xy_y_xy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_x_yy[i] = 4.0 * g_xy_y_xy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_x_yz[i] = 4.0 * g_xy_y_xy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_x_zz[i] = 4.0 * g_xy_y_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (672-678)

    #pragma omp simd aligned(g_xy_y_0_xx, g_xy_y_0_xy, g_xy_y_0_xz, g_xy_y_0_yy, g_xy_y_0_yz, g_xy_y_0_zz, g_xy_y_yy_xx, g_xy_y_yy_xy, g_xy_y_yy_xz, g_xy_y_yy_yy, g_xy_y_yy_yz, g_xy_y_yy_zz, g_y_0_y_0_x_y_y_xx, g_y_0_y_0_x_y_y_xy, g_y_0_y_0_x_y_y_xz, g_y_0_y_0_x_y_y_yy, g_y_0_y_0_x_y_y_yz, g_y_0_y_0_x_y_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_y_y_xx[i] = -2.0 * g_xy_y_0_xx[i] * a_exp + 4.0 * g_xy_y_yy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_y_xy[i] = -2.0 * g_xy_y_0_xy[i] * a_exp + 4.0 * g_xy_y_yy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_y_xz[i] = -2.0 * g_xy_y_0_xz[i] * a_exp + 4.0 * g_xy_y_yy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_y_yy[i] = -2.0 * g_xy_y_0_yy[i] * a_exp + 4.0 * g_xy_y_yy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_y_yz[i] = -2.0 * g_xy_y_0_yz[i] * a_exp + 4.0 * g_xy_y_yy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_y_zz[i] = -2.0 * g_xy_y_0_zz[i] * a_exp + 4.0 * g_xy_y_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (678-684)

    #pragma omp simd aligned(g_xy_y_yz_xx, g_xy_y_yz_xy, g_xy_y_yz_xz, g_xy_y_yz_yy, g_xy_y_yz_yz, g_xy_y_yz_zz, g_y_0_y_0_x_y_z_xx, g_y_0_y_0_x_y_z_xy, g_y_0_y_0_x_y_z_xz, g_y_0_y_0_x_y_z_yy, g_y_0_y_0_x_y_z_yz, g_y_0_y_0_x_y_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_y_z_xx[i] = 4.0 * g_xy_y_yz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_z_xy[i] = 4.0 * g_xy_y_yz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_z_xz[i] = 4.0 * g_xy_y_yz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_z_yy[i] = 4.0 * g_xy_y_yz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_z_yz[i] = 4.0 * g_xy_y_yz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_z_zz[i] = 4.0 * g_xy_y_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (684-690)

    #pragma omp simd aligned(g_xy_z_xy_xx, g_xy_z_xy_xy, g_xy_z_xy_xz, g_xy_z_xy_yy, g_xy_z_xy_yz, g_xy_z_xy_zz, g_y_0_y_0_x_z_x_xx, g_y_0_y_0_x_z_x_xy, g_y_0_y_0_x_z_x_xz, g_y_0_y_0_x_z_x_yy, g_y_0_y_0_x_z_x_yz, g_y_0_y_0_x_z_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_z_x_xx[i] = 4.0 * g_xy_z_xy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_x_xy[i] = 4.0 * g_xy_z_xy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_x_xz[i] = 4.0 * g_xy_z_xy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_x_yy[i] = 4.0 * g_xy_z_xy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_x_yz[i] = 4.0 * g_xy_z_xy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_x_zz[i] = 4.0 * g_xy_z_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (690-696)

    #pragma omp simd aligned(g_xy_z_0_xx, g_xy_z_0_xy, g_xy_z_0_xz, g_xy_z_0_yy, g_xy_z_0_yz, g_xy_z_0_zz, g_xy_z_yy_xx, g_xy_z_yy_xy, g_xy_z_yy_xz, g_xy_z_yy_yy, g_xy_z_yy_yz, g_xy_z_yy_zz, g_y_0_y_0_x_z_y_xx, g_y_0_y_0_x_z_y_xy, g_y_0_y_0_x_z_y_xz, g_y_0_y_0_x_z_y_yy, g_y_0_y_0_x_z_y_yz, g_y_0_y_0_x_z_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_z_y_xx[i] = -2.0 * g_xy_z_0_xx[i] * a_exp + 4.0 * g_xy_z_yy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_y_xy[i] = -2.0 * g_xy_z_0_xy[i] * a_exp + 4.0 * g_xy_z_yy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_y_xz[i] = -2.0 * g_xy_z_0_xz[i] * a_exp + 4.0 * g_xy_z_yy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_y_yy[i] = -2.0 * g_xy_z_0_yy[i] * a_exp + 4.0 * g_xy_z_yy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_y_yz[i] = -2.0 * g_xy_z_0_yz[i] * a_exp + 4.0 * g_xy_z_yy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_y_zz[i] = -2.0 * g_xy_z_0_zz[i] * a_exp + 4.0 * g_xy_z_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (696-702)

    #pragma omp simd aligned(g_xy_z_yz_xx, g_xy_z_yz_xy, g_xy_z_yz_xz, g_xy_z_yz_yy, g_xy_z_yz_yz, g_xy_z_yz_zz, g_y_0_y_0_x_z_z_xx, g_y_0_y_0_x_z_z_xy, g_y_0_y_0_x_z_z_xz, g_y_0_y_0_x_z_z_yy, g_y_0_y_0_x_z_z_yz, g_y_0_y_0_x_z_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_z_z_xx[i] = 4.0 * g_xy_z_yz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_z_xy[i] = 4.0 * g_xy_z_yz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_z_xz[i] = 4.0 * g_xy_z_yz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_z_yy[i] = 4.0 * g_xy_z_yz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_z_yz[i] = 4.0 * g_xy_z_yz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_z_zz[i] = 4.0 * g_xy_z_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (702-708)

    #pragma omp simd aligned(g_0_x_xy_xx, g_0_x_xy_xy, g_0_x_xy_xz, g_0_x_xy_yy, g_0_x_xy_yz, g_0_x_xy_zz, g_y_0_y_0_y_x_x_xx, g_y_0_y_0_y_x_x_xy, g_y_0_y_0_y_x_x_xz, g_y_0_y_0_y_x_x_yy, g_y_0_y_0_y_x_x_yz, g_y_0_y_0_y_x_x_zz, g_yy_x_xy_xx, g_yy_x_xy_xy, g_yy_x_xy_xz, g_yy_x_xy_yy, g_yy_x_xy_yz, g_yy_x_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_x_x_xx[i] = -2.0 * g_0_x_xy_xx[i] * c_exps[i] + 4.0 * g_yy_x_xy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_x_xy[i] = -2.0 * g_0_x_xy_xy[i] * c_exps[i] + 4.0 * g_yy_x_xy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_x_xz[i] = -2.0 * g_0_x_xy_xz[i] * c_exps[i] + 4.0 * g_yy_x_xy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_x_yy[i] = -2.0 * g_0_x_xy_yy[i] * c_exps[i] + 4.0 * g_yy_x_xy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_x_yz[i] = -2.0 * g_0_x_xy_yz[i] * c_exps[i] + 4.0 * g_yy_x_xy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_x_zz[i] = -2.0 * g_0_x_xy_zz[i] * c_exps[i] + 4.0 * g_yy_x_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (708-714)

    #pragma omp simd aligned(g_0_x_0_xx, g_0_x_0_xy, g_0_x_0_xz, g_0_x_0_yy, g_0_x_0_yz, g_0_x_0_zz, g_0_x_yy_xx, g_0_x_yy_xy, g_0_x_yy_xz, g_0_x_yy_yy, g_0_x_yy_yz, g_0_x_yy_zz, g_y_0_y_0_y_x_y_xx, g_y_0_y_0_y_x_y_xy, g_y_0_y_0_y_x_y_xz, g_y_0_y_0_y_x_y_yy, g_y_0_y_0_y_x_y_yz, g_y_0_y_0_y_x_y_zz, g_yy_x_0_xx, g_yy_x_0_xy, g_yy_x_0_xz, g_yy_x_0_yy, g_yy_x_0_yz, g_yy_x_0_zz, g_yy_x_yy_xx, g_yy_x_yy_xy, g_yy_x_yy_xz, g_yy_x_yy_yy, g_yy_x_yy_yz, g_yy_x_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_x_y_xx[i] = g_0_x_0_xx[i] - 2.0 * g_0_x_yy_xx[i] * c_exps[i] - 2.0 * g_yy_x_0_xx[i] * a_exp + 4.0 * g_yy_x_yy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_y_xy[i] = g_0_x_0_xy[i] - 2.0 * g_0_x_yy_xy[i] * c_exps[i] - 2.0 * g_yy_x_0_xy[i] * a_exp + 4.0 * g_yy_x_yy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_y_xz[i] = g_0_x_0_xz[i] - 2.0 * g_0_x_yy_xz[i] * c_exps[i] - 2.0 * g_yy_x_0_xz[i] * a_exp + 4.0 * g_yy_x_yy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_y_yy[i] = g_0_x_0_yy[i] - 2.0 * g_0_x_yy_yy[i] * c_exps[i] - 2.0 * g_yy_x_0_yy[i] * a_exp + 4.0 * g_yy_x_yy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_y_yz[i] = g_0_x_0_yz[i] - 2.0 * g_0_x_yy_yz[i] * c_exps[i] - 2.0 * g_yy_x_0_yz[i] * a_exp + 4.0 * g_yy_x_yy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_y_zz[i] = g_0_x_0_zz[i] - 2.0 * g_0_x_yy_zz[i] * c_exps[i] - 2.0 * g_yy_x_0_zz[i] * a_exp + 4.0 * g_yy_x_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (714-720)

    #pragma omp simd aligned(g_0_x_yz_xx, g_0_x_yz_xy, g_0_x_yz_xz, g_0_x_yz_yy, g_0_x_yz_yz, g_0_x_yz_zz, g_y_0_y_0_y_x_z_xx, g_y_0_y_0_y_x_z_xy, g_y_0_y_0_y_x_z_xz, g_y_0_y_0_y_x_z_yy, g_y_0_y_0_y_x_z_yz, g_y_0_y_0_y_x_z_zz, g_yy_x_yz_xx, g_yy_x_yz_xy, g_yy_x_yz_xz, g_yy_x_yz_yy, g_yy_x_yz_yz, g_yy_x_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_x_z_xx[i] = -2.0 * g_0_x_yz_xx[i] * c_exps[i] + 4.0 * g_yy_x_yz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_z_xy[i] = -2.0 * g_0_x_yz_xy[i] * c_exps[i] + 4.0 * g_yy_x_yz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_z_xz[i] = -2.0 * g_0_x_yz_xz[i] * c_exps[i] + 4.0 * g_yy_x_yz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_z_yy[i] = -2.0 * g_0_x_yz_yy[i] * c_exps[i] + 4.0 * g_yy_x_yz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_z_yz[i] = -2.0 * g_0_x_yz_yz[i] * c_exps[i] + 4.0 * g_yy_x_yz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_z_zz[i] = -2.0 * g_0_x_yz_zz[i] * c_exps[i] + 4.0 * g_yy_x_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (720-726)

    #pragma omp simd aligned(g_0_y_xy_xx, g_0_y_xy_xy, g_0_y_xy_xz, g_0_y_xy_yy, g_0_y_xy_yz, g_0_y_xy_zz, g_y_0_y_0_y_y_x_xx, g_y_0_y_0_y_y_x_xy, g_y_0_y_0_y_y_x_xz, g_y_0_y_0_y_y_x_yy, g_y_0_y_0_y_y_x_yz, g_y_0_y_0_y_y_x_zz, g_yy_y_xy_xx, g_yy_y_xy_xy, g_yy_y_xy_xz, g_yy_y_xy_yy, g_yy_y_xy_yz, g_yy_y_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_y_x_xx[i] = -2.0 * g_0_y_xy_xx[i] * c_exps[i] + 4.0 * g_yy_y_xy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_x_xy[i] = -2.0 * g_0_y_xy_xy[i] * c_exps[i] + 4.0 * g_yy_y_xy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_x_xz[i] = -2.0 * g_0_y_xy_xz[i] * c_exps[i] + 4.0 * g_yy_y_xy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_x_yy[i] = -2.0 * g_0_y_xy_yy[i] * c_exps[i] + 4.0 * g_yy_y_xy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_x_yz[i] = -2.0 * g_0_y_xy_yz[i] * c_exps[i] + 4.0 * g_yy_y_xy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_x_zz[i] = -2.0 * g_0_y_xy_zz[i] * c_exps[i] + 4.0 * g_yy_y_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (726-732)

    #pragma omp simd aligned(g_0_y_0_xx, g_0_y_0_xy, g_0_y_0_xz, g_0_y_0_yy, g_0_y_0_yz, g_0_y_0_zz, g_0_y_yy_xx, g_0_y_yy_xy, g_0_y_yy_xz, g_0_y_yy_yy, g_0_y_yy_yz, g_0_y_yy_zz, g_y_0_y_0_y_y_y_xx, g_y_0_y_0_y_y_y_xy, g_y_0_y_0_y_y_y_xz, g_y_0_y_0_y_y_y_yy, g_y_0_y_0_y_y_y_yz, g_y_0_y_0_y_y_y_zz, g_yy_y_0_xx, g_yy_y_0_xy, g_yy_y_0_xz, g_yy_y_0_yy, g_yy_y_0_yz, g_yy_y_0_zz, g_yy_y_yy_xx, g_yy_y_yy_xy, g_yy_y_yy_xz, g_yy_y_yy_yy, g_yy_y_yy_yz, g_yy_y_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_y_y_xx[i] = g_0_y_0_xx[i] - 2.0 * g_0_y_yy_xx[i] * c_exps[i] - 2.0 * g_yy_y_0_xx[i] * a_exp + 4.0 * g_yy_y_yy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_y_xy[i] = g_0_y_0_xy[i] - 2.0 * g_0_y_yy_xy[i] * c_exps[i] - 2.0 * g_yy_y_0_xy[i] * a_exp + 4.0 * g_yy_y_yy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_y_xz[i] = g_0_y_0_xz[i] - 2.0 * g_0_y_yy_xz[i] * c_exps[i] - 2.0 * g_yy_y_0_xz[i] * a_exp + 4.0 * g_yy_y_yy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_y_yy[i] = g_0_y_0_yy[i] - 2.0 * g_0_y_yy_yy[i] * c_exps[i] - 2.0 * g_yy_y_0_yy[i] * a_exp + 4.0 * g_yy_y_yy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_y_yz[i] = g_0_y_0_yz[i] - 2.0 * g_0_y_yy_yz[i] * c_exps[i] - 2.0 * g_yy_y_0_yz[i] * a_exp + 4.0 * g_yy_y_yy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_y_zz[i] = g_0_y_0_zz[i] - 2.0 * g_0_y_yy_zz[i] * c_exps[i] - 2.0 * g_yy_y_0_zz[i] * a_exp + 4.0 * g_yy_y_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (732-738)

    #pragma omp simd aligned(g_0_y_yz_xx, g_0_y_yz_xy, g_0_y_yz_xz, g_0_y_yz_yy, g_0_y_yz_yz, g_0_y_yz_zz, g_y_0_y_0_y_y_z_xx, g_y_0_y_0_y_y_z_xy, g_y_0_y_0_y_y_z_xz, g_y_0_y_0_y_y_z_yy, g_y_0_y_0_y_y_z_yz, g_y_0_y_0_y_y_z_zz, g_yy_y_yz_xx, g_yy_y_yz_xy, g_yy_y_yz_xz, g_yy_y_yz_yy, g_yy_y_yz_yz, g_yy_y_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_y_z_xx[i] = -2.0 * g_0_y_yz_xx[i] * c_exps[i] + 4.0 * g_yy_y_yz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_z_xy[i] = -2.0 * g_0_y_yz_xy[i] * c_exps[i] + 4.0 * g_yy_y_yz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_z_xz[i] = -2.0 * g_0_y_yz_xz[i] * c_exps[i] + 4.0 * g_yy_y_yz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_z_yy[i] = -2.0 * g_0_y_yz_yy[i] * c_exps[i] + 4.0 * g_yy_y_yz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_z_yz[i] = -2.0 * g_0_y_yz_yz[i] * c_exps[i] + 4.0 * g_yy_y_yz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_z_zz[i] = -2.0 * g_0_y_yz_zz[i] * c_exps[i] + 4.0 * g_yy_y_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (738-744)

    #pragma omp simd aligned(g_0_z_xy_xx, g_0_z_xy_xy, g_0_z_xy_xz, g_0_z_xy_yy, g_0_z_xy_yz, g_0_z_xy_zz, g_y_0_y_0_y_z_x_xx, g_y_0_y_0_y_z_x_xy, g_y_0_y_0_y_z_x_xz, g_y_0_y_0_y_z_x_yy, g_y_0_y_0_y_z_x_yz, g_y_0_y_0_y_z_x_zz, g_yy_z_xy_xx, g_yy_z_xy_xy, g_yy_z_xy_xz, g_yy_z_xy_yy, g_yy_z_xy_yz, g_yy_z_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_z_x_xx[i] = -2.0 * g_0_z_xy_xx[i] * c_exps[i] + 4.0 * g_yy_z_xy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_x_xy[i] = -2.0 * g_0_z_xy_xy[i] * c_exps[i] + 4.0 * g_yy_z_xy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_x_xz[i] = -2.0 * g_0_z_xy_xz[i] * c_exps[i] + 4.0 * g_yy_z_xy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_x_yy[i] = -2.0 * g_0_z_xy_yy[i] * c_exps[i] + 4.0 * g_yy_z_xy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_x_yz[i] = -2.0 * g_0_z_xy_yz[i] * c_exps[i] + 4.0 * g_yy_z_xy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_x_zz[i] = -2.0 * g_0_z_xy_zz[i] * c_exps[i] + 4.0 * g_yy_z_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (744-750)

    #pragma omp simd aligned(g_0_z_0_xx, g_0_z_0_xy, g_0_z_0_xz, g_0_z_0_yy, g_0_z_0_yz, g_0_z_0_zz, g_0_z_yy_xx, g_0_z_yy_xy, g_0_z_yy_xz, g_0_z_yy_yy, g_0_z_yy_yz, g_0_z_yy_zz, g_y_0_y_0_y_z_y_xx, g_y_0_y_0_y_z_y_xy, g_y_0_y_0_y_z_y_xz, g_y_0_y_0_y_z_y_yy, g_y_0_y_0_y_z_y_yz, g_y_0_y_0_y_z_y_zz, g_yy_z_0_xx, g_yy_z_0_xy, g_yy_z_0_xz, g_yy_z_0_yy, g_yy_z_0_yz, g_yy_z_0_zz, g_yy_z_yy_xx, g_yy_z_yy_xy, g_yy_z_yy_xz, g_yy_z_yy_yy, g_yy_z_yy_yz, g_yy_z_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_z_y_xx[i] = g_0_z_0_xx[i] - 2.0 * g_0_z_yy_xx[i] * c_exps[i] - 2.0 * g_yy_z_0_xx[i] * a_exp + 4.0 * g_yy_z_yy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_y_xy[i] = g_0_z_0_xy[i] - 2.0 * g_0_z_yy_xy[i] * c_exps[i] - 2.0 * g_yy_z_0_xy[i] * a_exp + 4.0 * g_yy_z_yy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_y_xz[i] = g_0_z_0_xz[i] - 2.0 * g_0_z_yy_xz[i] * c_exps[i] - 2.0 * g_yy_z_0_xz[i] * a_exp + 4.0 * g_yy_z_yy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_y_yy[i] = g_0_z_0_yy[i] - 2.0 * g_0_z_yy_yy[i] * c_exps[i] - 2.0 * g_yy_z_0_yy[i] * a_exp + 4.0 * g_yy_z_yy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_y_yz[i] = g_0_z_0_yz[i] - 2.0 * g_0_z_yy_yz[i] * c_exps[i] - 2.0 * g_yy_z_0_yz[i] * a_exp + 4.0 * g_yy_z_yy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_y_zz[i] = g_0_z_0_zz[i] - 2.0 * g_0_z_yy_zz[i] * c_exps[i] - 2.0 * g_yy_z_0_zz[i] * a_exp + 4.0 * g_yy_z_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (750-756)

    #pragma omp simd aligned(g_0_z_yz_xx, g_0_z_yz_xy, g_0_z_yz_xz, g_0_z_yz_yy, g_0_z_yz_yz, g_0_z_yz_zz, g_y_0_y_0_y_z_z_xx, g_y_0_y_0_y_z_z_xy, g_y_0_y_0_y_z_z_xz, g_y_0_y_0_y_z_z_yy, g_y_0_y_0_y_z_z_yz, g_y_0_y_0_y_z_z_zz, g_yy_z_yz_xx, g_yy_z_yz_xy, g_yy_z_yz_xz, g_yy_z_yz_yy, g_yy_z_yz_yz, g_yy_z_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_z_z_xx[i] = -2.0 * g_0_z_yz_xx[i] * c_exps[i] + 4.0 * g_yy_z_yz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_z_xy[i] = -2.0 * g_0_z_yz_xy[i] * c_exps[i] + 4.0 * g_yy_z_yz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_z_xz[i] = -2.0 * g_0_z_yz_xz[i] * c_exps[i] + 4.0 * g_yy_z_yz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_z_yy[i] = -2.0 * g_0_z_yz_yy[i] * c_exps[i] + 4.0 * g_yy_z_yz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_z_yz[i] = -2.0 * g_0_z_yz_yz[i] * c_exps[i] + 4.0 * g_yy_z_yz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_z_zz[i] = -2.0 * g_0_z_yz_zz[i] * c_exps[i] + 4.0 * g_yy_z_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (756-762)

    #pragma omp simd aligned(g_y_0_y_0_z_x_x_xx, g_y_0_y_0_z_x_x_xy, g_y_0_y_0_z_x_x_xz, g_y_0_y_0_z_x_x_yy, g_y_0_y_0_z_x_x_yz, g_y_0_y_0_z_x_x_zz, g_yz_x_xy_xx, g_yz_x_xy_xy, g_yz_x_xy_xz, g_yz_x_xy_yy, g_yz_x_xy_yz, g_yz_x_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_x_x_xx[i] = 4.0 * g_yz_x_xy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_x_xy[i] = 4.0 * g_yz_x_xy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_x_xz[i] = 4.0 * g_yz_x_xy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_x_yy[i] = 4.0 * g_yz_x_xy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_x_yz[i] = 4.0 * g_yz_x_xy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_x_zz[i] = 4.0 * g_yz_x_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (762-768)

    #pragma omp simd aligned(g_y_0_y_0_z_x_y_xx, g_y_0_y_0_z_x_y_xy, g_y_0_y_0_z_x_y_xz, g_y_0_y_0_z_x_y_yy, g_y_0_y_0_z_x_y_yz, g_y_0_y_0_z_x_y_zz, g_yz_x_0_xx, g_yz_x_0_xy, g_yz_x_0_xz, g_yz_x_0_yy, g_yz_x_0_yz, g_yz_x_0_zz, g_yz_x_yy_xx, g_yz_x_yy_xy, g_yz_x_yy_xz, g_yz_x_yy_yy, g_yz_x_yy_yz, g_yz_x_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_x_y_xx[i] = -2.0 * g_yz_x_0_xx[i] * a_exp + 4.0 * g_yz_x_yy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_y_xy[i] = -2.0 * g_yz_x_0_xy[i] * a_exp + 4.0 * g_yz_x_yy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_y_xz[i] = -2.0 * g_yz_x_0_xz[i] * a_exp + 4.0 * g_yz_x_yy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_y_yy[i] = -2.0 * g_yz_x_0_yy[i] * a_exp + 4.0 * g_yz_x_yy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_y_yz[i] = -2.0 * g_yz_x_0_yz[i] * a_exp + 4.0 * g_yz_x_yy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_y_zz[i] = -2.0 * g_yz_x_0_zz[i] * a_exp + 4.0 * g_yz_x_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (768-774)

    #pragma omp simd aligned(g_y_0_y_0_z_x_z_xx, g_y_0_y_0_z_x_z_xy, g_y_0_y_0_z_x_z_xz, g_y_0_y_0_z_x_z_yy, g_y_0_y_0_z_x_z_yz, g_y_0_y_0_z_x_z_zz, g_yz_x_yz_xx, g_yz_x_yz_xy, g_yz_x_yz_xz, g_yz_x_yz_yy, g_yz_x_yz_yz, g_yz_x_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_x_z_xx[i] = 4.0 * g_yz_x_yz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_z_xy[i] = 4.0 * g_yz_x_yz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_z_xz[i] = 4.0 * g_yz_x_yz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_z_yy[i] = 4.0 * g_yz_x_yz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_z_yz[i] = 4.0 * g_yz_x_yz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_z_zz[i] = 4.0 * g_yz_x_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (774-780)

    #pragma omp simd aligned(g_y_0_y_0_z_y_x_xx, g_y_0_y_0_z_y_x_xy, g_y_0_y_0_z_y_x_xz, g_y_0_y_0_z_y_x_yy, g_y_0_y_0_z_y_x_yz, g_y_0_y_0_z_y_x_zz, g_yz_y_xy_xx, g_yz_y_xy_xy, g_yz_y_xy_xz, g_yz_y_xy_yy, g_yz_y_xy_yz, g_yz_y_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_y_x_xx[i] = 4.0 * g_yz_y_xy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_x_xy[i] = 4.0 * g_yz_y_xy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_x_xz[i] = 4.0 * g_yz_y_xy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_x_yy[i] = 4.0 * g_yz_y_xy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_x_yz[i] = 4.0 * g_yz_y_xy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_x_zz[i] = 4.0 * g_yz_y_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (780-786)

    #pragma omp simd aligned(g_y_0_y_0_z_y_y_xx, g_y_0_y_0_z_y_y_xy, g_y_0_y_0_z_y_y_xz, g_y_0_y_0_z_y_y_yy, g_y_0_y_0_z_y_y_yz, g_y_0_y_0_z_y_y_zz, g_yz_y_0_xx, g_yz_y_0_xy, g_yz_y_0_xz, g_yz_y_0_yy, g_yz_y_0_yz, g_yz_y_0_zz, g_yz_y_yy_xx, g_yz_y_yy_xy, g_yz_y_yy_xz, g_yz_y_yy_yy, g_yz_y_yy_yz, g_yz_y_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_y_y_xx[i] = -2.0 * g_yz_y_0_xx[i] * a_exp + 4.0 * g_yz_y_yy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_y_xy[i] = -2.0 * g_yz_y_0_xy[i] * a_exp + 4.0 * g_yz_y_yy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_y_xz[i] = -2.0 * g_yz_y_0_xz[i] * a_exp + 4.0 * g_yz_y_yy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_y_yy[i] = -2.0 * g_yz_y_0_yy[i] * a_exp + 4.0 * g_yz_y_yy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_y_yz[i] = -2.0 * g_yz_y_0_yz[i] * a_exp + 4.0 * g_yz_y_yy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_y_zz[i] = -2.0 * g_yz_y_0_zz[i] * a_exp + 4.0 * g_yz_y_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (786-792)

    #pragma omp simd aligned(g_y_0_y_0_z_y_z_xx, g_y_0_y_0_z_y_z_xy, g_y_0_y_0_z_y_z_xz, g_y_0_y_0_z_y_z_yy, g_y_0_y_0_z_y_z_yz, g_y_0_y_0_z_y_z_zz, g_yz_y_yz_xx, g_yz_y_yz_xy, g_yz_y_yz_xz, g_yz_y_yz_yy, g_yz_y_yz_yz, g_yz_y_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_y_z_xx[i] = 4.0 * g_yz_y_yz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_z_xy[i] = 4.0 * g_yz_y_yz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_z_xz[i] = 4.0 * g_yz_y_yz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_z_yy[i] = 4.0 * g_yz_y_yz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_z_yz[i] = 4.0 * g_yz_y_yz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_z_zz[i] = 4.0 * g_yz_y_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (792-798)

    #pragma omp simd aligned(g_y_0_y_0_z_z_x_xx, g_y_0_y_0_z_z_x_xy, g_y_0_y_0_z_z_x_xz, g_y_0_y_0_z_z_x_yy, g_y_0_y_0_z_z_x_yz, g_y_0_y_0_z_z_x_zz, g_yz_z_xy_xx, g_yz_z_xy_xy, g_yz_z_xy_xz, g_yz_z_xy_yy, g_yz_z_xy_yz, g_yz_z_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_z_x_xx[i] = 4.0 * g_yz_z_xy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_x_xy[i] = 4.0 * g_yz_z_xy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_x_xz[i] = 4.0 * g_yz_z_xy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_x_yy[i] = 4.0 * g_yz_z_xy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_x_yz[i] = 4.0 * g_yz_z_xy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_x_zz[i] = 4.0 * g_yz_z_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (798-804)

    #pragma omp simd aligned(g_y_0_y_0_z_z_y_xx, g_y_0_y_0_z_z_y_xy, g_y_0_y_0_z_z_y_xz, g_y_0_y_0_z_z_y_yy, g_y_0_y_0_z_z_y_yz, g_y_0_y_0_z_z_y_zz, g_yz_z_0_xx, g_yz_z_0_xy, g_yz_z_0_xz, g_yz_z_0_yy, g_yz_z_0_yz, g_yz_z_0_zz, g_yz_z_yy_xx, g_yz_z_yy_xy, g_yz_z_yy_xz, g_yz_z_yy_yy, g_yz_z_yy_yz, g_yz_z_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_z_y_xx[i] = -2.0 * g_yz_z_0_xx[i] * a_exp + 4.0 * g_yz_z_yy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_y_xy[i] = -2.0 * g_yz_z_0_xy[i] * a_exp + 4.0 * g_yz_z_yy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_y_xz[i] = -2.0 * g_yz_z_0_xz[i] * a_exp + 4.0 * g_yz_z_yy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_y_yy[i] = -2.0 * g_yz_z_0_yy[i] * a_exp + 4.0 * g_yz_z_yy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_y_yz[i] = -2.0 * g_yz_z_0_yz[i] * a_exp + 4.0 * g_yz_z_yy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_y_zz[i] = -2.0 * g_yz_z_0_zz[i] * a_exp + 4.0 * g_yz_z_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (804-810)

    #pragma omp simd aligned(g_y_0_y_0_z_z_z_xx, g_y_0_y_0_z_z_z_xy, g_y_0_y_0_z_z_z_xz, g_y_0_y_0_z_z_z_yy, g_y_0_y_0_z_z_z_yz, g_y_0_y_0_z_z_z_zz, g_yz_z_yz_xx, g_yz_z_yz_xy, g_yz_z_yz_xz, g_yz_z_yz_yy, g_yz_z_yz_yz, g_yz_z_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_z_z_xx[i] = 4.0 * g_yz_z_yz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_z_xy[i] = 4.0 * g_yz_z_yz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_z_xz[i] = 4.0 * g_yz_z_yz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_z_yy[i] = 4.0 * g_yz_z_yz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_z_yz[i] = 4.0 * g_yz_z_yz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_z_zz[i] = 4.0 * g_yz_z_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (810-816)

    #pragma omp simd aligned(g_xy_x_xz_xx, g_xy_x_xz_xy, g_xy_x_xz_xz, g_xy_x_xz_yy, g_xy_x_xz_yz, g_xy_x_xz_zz, g_y_0_z_0_x_x_x_xx, g_y_0_z_0_x_x_x_xy, g_y_0_z_0_x_x_x_xz, g_y_0_z_0_x_x_x_yy, g_y_0_z_0_x_x_x_yz, g_y_0_z_0_x_x_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_x_x_xx[i] = 4.0 * g_xy_x_xz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_x_xy[i] = 4.0 * g_xy_x_xz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_x_xz[i] = 4.0 * g_xy_x_xz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_x_yy[i] = 4.0 * g_xy_x_xz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_x_yz[i] = 4.0 * g_xy_x_xz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_x_zz[i] = 4.0 * g_xy_x_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (816-822)

    #pragma omp simd aligned(g_xy_x_yz_xx, g_xy_x_yz_xy, g_xy_x_yz_xz, g_xy_x_yz_yy, g_xy_x_yz_yz, g_xy_x_yz_zz, g_y_0_z_0_x_x_y_xx, g_y_0_z_0_x_x_y_xy, g_y_0_z_0_x_x_y_xz, g_y_0_z_0_x_x_y_yy, g_y_0_z_0_x_x_y_yz, g_y_0_z_0_x_x_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_x_y_xx[i] = 4.0 * g_xy_x_yz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_y_xy[i] = 4.0 * g_xy_x_yz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_y_xz[i] = 4.0 * g_xy_x_yz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_y_yy[i] = 4.0 * g_xy_x_yz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_y_yz[i] = 4.0 * g_xy_x_yz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_y_zz[i] = 4.0 * g_xy_x_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (822-828)

    #pragma omp simd aligned(g_xy_x_0_xx, g_xy_x_0_xy, g_xy_x_0_xz, g_xy_x_0_yy, g_xy_x_0_yz, g_xy_x_0_zz, g_xy_x_zz_xx, g_xy_x_zz_xy, g_xy_x_zz_xz, g_xy_x_zz_yy, g_xy_x_zz_yz, g_xy_x_zz_zz, g_y_0_z_0_x_x_z_xx, g_y_0_z_0_x_x_z_xy, g_y_0_z_0_x_x_z_xz, g_y_0_z_0_x_x_z_yy, g_y_0_z_0_x_x_z_yz, g_y_0_z_0_x_x_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_x_z_xx[i] = -2.0 * g_xy_x_0_xx[i] * a_exp + 4.0 * g_xy_x_zz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_z_xy[i] = -2.0 * g_xy_x_0_xy[i] * a_exp + 4.0 * g_xy_x_zz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_z_xz[i] = -2.0 * g_xy_x_0_xz[i] * a_exp + 4.0 * g_xy_x_zz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_z_yy[i] = -2.0 * g_xy_x_0_yy[i] * a_exp + 4.0 * g_xy_x_zz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_z_yz[i] = -2.0 * g_xy_x_0_yz[i] * a_exp + 4.0 * g_xy_x_zz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_z_zz[i] = -2.0 * g_xy_x_0_zz[i] * a_exp + 4.0 * g_xy_x_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (828-834)

    #pragma omp simd aligned(g_xy_y_xz_xx, g_xy_y_xz_xy, g_xy_y_xz_xz, g_xy_y_xz_yy, g_xy_y_xz_yz, g_xy_y_xz_zz, g_y_0_z_0_x_y_x_xx, g_y_0_z_0_x_y_x_xy, g_y_0_z_0_x_y_x_xz, g_y_0_z_0_x_y_x_yy, g_y_0_z_0_x_y_x_yz, g_y_0_z_0_x_y_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_y_x_xx[i] = 4.0 * g_xy_y_xz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_x_xy[i] = 4.0 * g_xy_y_xz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_x_xz[i] = 4.0 * g_xy_y_xz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_x_yy[i] = 4.0 * g_xy_y_xz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_x_yz[i] = 4.0 * g_xy_y_xz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_x_zz[i] = 4.0 * g_xy_y_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (834-840)

    #pragma omp simd aligned(g_xy_y_yz_xx, g_xy_y_yz_xy, g_xy_y_yz_xz, g_xy_y_yz_yy, g_xy_y_yz_yz, g_xy_y_yz_zz, g_y_0_z_0_x_y_y_xx, g_y_0_z_0_x_y_y_xy, g_y_0_z_0_x_y_y_xz, g_y_0_z_0_x_y_y_yy, g_y_0_z_0_x_y_y_yz, g_y_0_z_0_x_y_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_y_y_xx[i] = 4.0 * g_xy_y_yz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_y_xy[i] = 4.0 * g_xy_y_yz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_y_xz[i] = 4.0 * g_xy_y_yz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_y_yy[i] = 4.0 * g_xy_y_yz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_y_yz[i] = 4.0 * g_xy_y_yz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_y_zz[i] = 4.0 * g_xy_y_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (840-846)

    #pragma omp simd aligned(g_xy_y_0_xx, g_xy_y_0_xy, g_xy_y_0_xz, g_xy_y_0_yy, g_xy_y_0_yz, g_xy_y_0_zz, g_xy_y_zz_xx, g_xy_y_zz_xy, g_xy_y_zz_xz, g_xy_y_zz_yy, g_xy_y_zz_yz, g_xy_y_zz_zz, g_y_0_z_0_x_y_z_xx, g_y_0_z_0_x_y_z_xy, g_y_0_z_0_x_y_z_xz, g_y_0_z_0_x_y_z_yy, g_y_0_z_0_x_y_z_yz, g_y_0_z_0_x_y_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_y_z_xx[i] = -2.0 * g_xy_y_0_xx[i] * a_exp + 4.0 * g_xy_y_zz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_z_xy[i] = -2.0 * g_xy_y_0_xy[i] * a_exp + 4.0 * g_xy_y_zz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_z_xz[i] = -2.0 * g_xy_y_0_xz[i] * a_exp + 4.0 * g_xy_y_zz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_z_yy[i] = -2.0 * g_xy_y_0_yy[i] * a_exp + 4.0 * g_xy_y_zz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_z_yz[i] = -2.0 * g_xy_y_0_yz[i] * a_exp + 4.0 * g_xy_y_zz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_z_zz[i] = -2.0 * g_xy_y_0_zz[i] * a_exp + 4.0 * g_xy_y_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (846-852)

    #pragma omp simd aligned(g_xy_z_xz_xx, g_xy_z_xz_xy, g_xy_z_xz_xz, g_xy_z_xz_yy, g_xy_z_xz_yz, g_xy_z_xz_zz, g_y_0_z_0_x_z_x_xx, g_y_0_z_0_x_z_x_xy, g_y_0_z_0_x_z_x_xz, g_y_0_z_0_x_z_x_yy, g_y_0_z_0_x_z_x_yz, g_y_0_z_0_x_z_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_z_x_xx[i] = 4.0 * g_xy_z_xz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_x_xy[i] = 4.0 * g_xy_z_xz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_x_xz[i] = 4.0 * g_xy_z_xz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_x_yy[i] = 4.0 * g_xy_z_xz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_x_yz[i] = 4.0 * g_xy_z_xz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_x_zz[i] = 4.0 * g_xy_z_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (852-858)

    #pragma omp simd aligned(g_xy_z_yz_xx, g_xy_z_yz_xy, g_xy_z_yz_xz, g_xy_z_yz_yy, g_xy_z_yz_yz, g_xy_z_yz_zz, g_y_0_z_0_x_z_y_xx, g_y_0_z_0_x_z_y_xy, g_y_0_z_0_x_z_y_xz, g_y_0_z_0_x_z_y_yy, g_y_0_z_0_x_z_y_yz, g_y_0_z_0_x_z_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_z_y_xx[i] = 4.0 * g_xy_z_yz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_y_xy[i] = 4.0 * g_xy_z_yz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_y_xz[i] = 4.0 * g_xy_z_yz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_y_yy[i] = 4.0 * g_xy_z_yz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_y_yz[i] = 4.0 * g_xy_z_yz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_y_zz[i] = 4.0 * g_xy_z_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (858-864)

    #pragma omp simd aligned(g_xy_z_0_xx, g_xy_z_0_xy, g_xy_z_0_xz, g_xy_z_0_yy, g_xy_z_0_yz, g_xy_z_0_zz, g_xy_z_zz_xx, g_xy_z_zz_xy, g_xy_z_zz_xz, g_xy_z_zz_yy, g_xy_z_zz_yz, g_xy_z_zz_zz, g_y_0_z_0_x_z_z_xx, g_y_0_z_0_x_z_z_xy, g_y_0_z_0_x_z_z_xz, g_y_0_z_0_x_z_z_yy, g_y_0_z_0_x_z_z_yz, g_y_0_z_0_x_z_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_z_z_xx[i] = -2.0 * g_xy_z_0_xx[i] * a_exp + 4.0 * g_xy_z_zz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_z_xy[i] = -2.0 * g_xy_z_0_xy[i] * a_exp + 4.0 * g_xy_z_zz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_z_xz[i] = -2.0 * g_xy_z_0_xz[i] * a_exp + 4.0 * g_xy_z_zz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_z_yy[i] = -2.0 * g_xy_z_0_yy[i] * a_exp + 4.0 * g_xy_z_zz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_z_yz[i] = -2.0 * g_xy_z_0_yz[i] * a_exp + 4.0 * g_xy_z_zz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_z_zz[i] = -2.0 * g_xy_z_0_zz[i] * a_exp + 4.0 * g_xy_z_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (864-870)

    #pragma omp simd aligned(g_0_x_xz_xx, g_0_x_xz_xy, g_0_x_xz_xz, g_0_x_xz_yy, g_0_x_xz_yz, g_0_x_xz_zz, g_y_0_z_0_y_x_x_xx, g_y_0_z_0_y_x_x_xy, g_y_0_z_0_y_x_x_xz, g_y_0_z_0_y_x_x_yy, g_y_0_z_0_y_x_x_yz, g_y_0_z_0_y_x_x_zz, g_yy_x_xz_xx, g_yy_x_xz_xy, g_yy_x_xz_xz, g_yy_x_xz_yy, g_yy_x_xz_yz, g_yy_x_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_x_x_xx[i] = -2.0 * g_0_x_xz_xx[i] * c_exps[i] + 4.0 * g_yy_x_xz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_x_xy[i] = -2.0 * g_0_x_xz_xy[i] * c_exps[i] + 4.0 * g_yy_x_xz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_x_xz[i] = -2.0 * g_0_x_xz_xz[i] * c_exps[i] + 4.0 * g_yy_x_xz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_x_yy[i] = -2.0 * g_0_x_xz_yy[i] * c_exps[i] + 4.0 * g_yy_x_xz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_x_yz[i] = -2.0 * g_0_x_xz_yz[i] * c_exps[i] + 4.0 * g_yy_x_xz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_x_zz[i] = -2.0 * g_0_x_xz_zz[i] * c_exps[i] + 4.0 * g_yy_x_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (870-876)

    #pragma omp simd aligned(g_0_x_yz_xx, g_0_x_yz_xy, g_0_x_yz_xz, g_0_x_yz_yy, g_0_x_yz_yz, g_0_x_yz_zz, g_y_0_z_0_y_x_y_xx, g_y_0_z_0_y_x_y_xy, g_y_0_z_0_y_x_y_xz, g_y_0_z_0_y_x_y_yy, g_y_0_z_0_y_x_y_yz, g_y_0_z_0_y_x_y_zz, g_yy_x_yz_xx, g_yy_x_yz_xy, g_yy_x_yz_xz, g_yy_x_yz_yy, g_yy_x_yz_yz, g_yy_x_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_x_y_xx[i] = -2.0 * g_0_x_yz_xx[i] * c_exps[i] + 4.0 * g_yy_x_yz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_y_xy[i] = -2.0 * g_0_x_yz_xy[i] * c_exps[i] + 4.0 * g_yy_x_yz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_y_xz[i] = -2.0 * g_0_x_yz_xz[i] * c_exps[i] + 4.0 * g_yy_x_yz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_y_yy[i] = -2.0 * g_0_x_yz_yy[i] * c_exps[i] + 4.0 * g_yy_x_yz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_y_yz[i] = -2.0 * g_0_x_yz_yz[i] * c_exps[i] + 4.0 * g_yy_x_yz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_y_zz[i] = -2.0 * g_0_x_yz_zz[i] * c_exps[i] + 4.0 * g_yy_x_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (876-882)

    #pragma omp simd aligned(g_0_x_0_xx, g_0_x_0_xy, g_0_x_0_xz, g_0_x_0_yy, g_0_x_0_yz, g_0_x_0_zz, g_0_x_zz_xx, g_0_x_zz_xy, g_0_x_zz_xz, g_0_x_zz_yy, g_0_x_zz_yz, g_0_x_zz_zz, g_y_0_z_0_y_x_z_xx, g_y_0_z_0_y_x_z_xy, g_y_0_z_0_y_x_z_xz, g_y_0_z_0_y_x_z_yy, g_y_0_z_0_y_x_z_yz, g_y_0_z_0_y_x_z_zz, g_yy_x_0_xx, g_yy_x_0_xy, g_yy_x_0_xz, g_yy_x_0_yy, g_yy_x_0_yz, g_yy_x_0_zz, g_yy_x_zz_xx, g_yy_x_zz_xy, g_yy_x_zz_xz, g_yy_x_zz_yy, g_yy_x_zz_yz, g_yy_x_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_x_z_xx[i] = g_0_x_0_xx[i] - 2.0 * g_0_x_zz_xx[i] * c_exps[i] - 2.0 * g_yy_x_0_xx[i] * a_exp + 4.0 * g_yy_x_zz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_z_xy[i] = g_0_x_0_xy[i] - 2.0 * g_0_x_zz_xy[i] * c_exps[i] - 2.0 * g_yy_x_0_xy[i] * a_exp + 4.0 * g_yy_x_zz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_z_xz[i] = g_0_x_0_xz[i] - 2.0 * g_0_x_zz_xz[i] * c_exps[i] - 2.0 * g_yy_x_0_xz[i] * a_exp + 4.0 * g_yy_x_zz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_z_yy[i] = g_0_x_0_yy[i] - 2.0 * g_0_x_zz_yy[i] * c_exps[i] - 2.0 * g_yy_x_0_yy[i] * a_exp + 4.0 * g_yy_x_zz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_z_yz[i] = g_0_x_0_yz[i] - 2.0 * g_0_x_zz_yz[i] * c_exps[i] - 2.0 * g_yy_x_0_yz[i] * a_exp + 4.0 * g_yy_x_zz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_z_zz[i] = g_0_x_0_zz[i] - 2.0 * g_0_x_zz_zz[i] * c_exps[i] - 2.0 * g_yy_x_0_zz[i] * a_exp + 4.0 * g_yy_x_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (882-888)

    #pragma omp simd aligned(g_0_y_xz_xx, g_0_y_xz_xy, g_0_y_xz_xz, g_0_y_xz_yy, g_0_y_xz_yz, g_0_y_xz_zz, g_y_0_z_0_y_y_x_xx, g_y_0_z_0_y_y_x_xy, g_y_0_z_0_y_y_x_xz, g_y_0_z_0_y_y_x_yy, g_y_0_z_0_y_y_x_yz, g_y_0_z_0_y_y_x_zz, g_yy_y_xz_xx, g_yy_y_xz_xy, g_yy_y_xz_xz, g_yy_y_xz_yy, g_yy_y_xz_yz, g_yy_y_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_y_x_xx[i] = -2.0 * g_0_y_xz_xx[i] * c_exps[i] + 4.0 * g_yy_y_xz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_x_xy[i] = -2.0 * g_0_y_xz_xy[i] * c_exps[i] + 4.0 * g_yy_y_xz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_x_xz[i] = -2.0 * g_0_y_xz_xz[i] * c_exps[i] + 4.0 * g_yy_y_xz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_x_yy[i] = -2.0 * g_0_y_xz_yy[i] * c_exps[i] + 4.0 * g_yy_y_xz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_x_yz[i] = -2.0 * g_0_y_xz_yz[i] * c_exps[i] + 4.0 * g_yy_y_xz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_x_zz[i] = -2.0 * g_0_y_xz_zz[i] * c_exps[i] + 4.0 * g_yy_y_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (888-894)

    #pragma omp simd aligned(g_0_y_yz_xx, g_0_y_yz_xy, g_0_y_yz_xz, g_0_y_yz_yy, g_0_y_yz_yz, g_0_y_yz_zz, g_y_0_z_0_y_y_y_xx, g_y_0_z_0_y_y_y_xy, g_y_0_z_0_y_y_y_xz, g_y_0_z_0_y_y_y_yy, g_y_0_z_0_y_y_y_yz, g_y_0_z_0_y_y_y_zz, g_yy_y_yz_xx, g_yy_y_yz_xy, g_yy_y_yz_xz, g_yy_y_yz_yy, g_yy_y_yz_yz, g_yy_y_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_y_y_xx[i] = -2.0 * g_0_y_yz_xx[i] * c_exps[i] + 4.0 * g_yy_y_yz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_y_xy[i] = -2.0 * g_0_y_yz_xy[i] * c_exps[i] + 4.0 * g_yy_y_yz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_y_xz[i] = -2.0 * g_0_y_yz_xz[i] * c_exps[i] + 4.0 * g_yy_y_yz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_y_yy[i] = -2.0 * g_0_y_yz_yy[i] * c_exps[i] + 4.0 * g_yy_y_yz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_y_yz[i] = -2.0 * g_0_y_yz_yz[i] * c_exps[i] + 4.0 * g_yy_y_yz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_y_zz[i] = -2.0 * g_0_y_yz_zz[i] * c_exps[i] + 4.0 * g_yy_y_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (894-900)

    #pragma omp simd aligned(g_0_y_0_xx, g_0_y_0_xy, g_0_y_0_xz, g_0_y_0_yy, g_0_y_0_yz, g_0_y_0_zz, g_0_y_zz_xx, g_0_y_zz_xy, g_0_y_zz_xz, g_0_y_zz_yy, g_0_y_zz_yz, g_0_y_zz_zz, g_y_0_z_0_y_y_z_xx, g_y_0_z_0_y_y_z_xy, g_y_0_z_0_y_y_z_xz, g_y_0_z_0_y_y_z_yy, g_y_0_z_0_y_y_z_yz, g_y_0_z_0_y_y_z_zz, g_yy_y_0_xx, g_yy_y_0_xy, g_yy_y_0_xz, g_yy_y_0_yy, g_yy_y_0_yz, g_yy_y_0_zz, g_yy_y_zz_xx, g_yy_y_zz_xy, g_yy_y_zz_xz, g_yy_y_zz_yy, g_yy_y_zz_yz, g_yy_y_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_y_z_xx[i] = g_0_y_0_xx[i] - 2.0 * g_0_y_zz_xx[i] * c_exps[i] - 2.0 * g_yy_y_0_xx[i] * a_exp + 4.0 * g_yy_y_zz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_z_xy[i] = g_0_y_0_xy[i] - 2.0 * g_0_y_zz_xy[i] * c_exps[i] - 2.0 * g_yy_y_0_xy[i] * a_exp + 4.0 * g_yy_y_zz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_z_xz[i] = g_0_y_0_xz[i] - 2.0 * g_0_y_zz_xz[i] * c_exps[i] - 2.0 * g_yy_y_0_xz[i] * a_exp + 4.0 * g_yy_y_zz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_z_yy[i] = g_0_y_0_yy[i] - 2.0 * g_0_y_zz_yy[i] * c_exps[i] - 2.0 * g_yy_y_0_yy[i] * a_exp + 4.0 * g_yy_y_zz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_z_yz[i] = g_0_y_0_yz[i] - 2.0 * g_0_y_zz_yz[i] * c_exps[i] - 2.0 * g_yy_y_0_yz[i] * a_exp + 4.0 * g_yy_y_zz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_z_zz[i] = g_0_y_0_zz[i] - 2.0 * g_0_y_zz_zz[i] * c_exps[i] - 2.0 * g_yy_y_0_zz[i] * a_exp + 4.0 * g_yy_y_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (900-906)

    #pragma omp simd aligned(g_0_z_xz_xx, g_0_z_xz_xy, g_0_z_xz_xz, g_0_z_xz_yy, g_0_z_xz_yz, g_0_z_xz_zz, g_y_0_z_0_y_z_x_xx, g_y_0_z_0_y_z_x_xy, g_y_0_z_0_y_z_x_xz, g_y_0_z_0_y_z_x_yy, g_y_0_z_0_y_z_x_yz, g_y_0_z_0_y_z_x_zz, g_yy_z_xz_xx, g_yy_z_xz_xy, g_yy_z_xz_xz, g_yy_z_xz_yy, g_yy_z_xz_yz, g_yy_z_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_z_x_xx[i] = -2.0 * g_0_z_xz_xx[i] * c_exps[i] + 4.0 * g_yy_z_xz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_x_xy[i] = -2.0 * g_0_z_xz_xy[i] * c_exps[i] + 4.0 * g_yy_z_xz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_x_xz[i] = -2.0 * g_0_z_xz_xz[i] * c_exps[i] + 4.0 * g_yy_z_xz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_x_yy[i] = -2.0 * g_0_z_xz_yy[i] * c_exps[i] + 4.0 * g_yy_z_xz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_x_yz[i] = -2.0 * g_0_z_xz_yz[i] * c_exps[i] + 4.0 * g_yy_z_xz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_x_zz[i] = -2.0 * g_0_z_xz_zz[i] * c_exps[i] + 4.0 * g_yy_z_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (906-912)

    #pragma omp simd aligned(g_0_z_yz_xx, g_0_z_yz_xy, g_0_z_yz_xz, g_0_z_yz_yy, g_0_z_yz_yz, g_0_z_yz_zz, g_y_0_z_0_y_z_y_xx, g_y_0_z_0_y_z_y_xy, g_y_0_z_0_y_z_y_xz, g_y_0_z_0_y_z_y_yy, g_y_0_z_0_y_z_y_yz, g_y_0_z_0_y_z_y_zz, g_yy_z_yz_xx, g_yy_z_yz_xy, g_yy_z_yz_xz, g_yy_z_yz_yy, g_yy_z_yz_yz, g_yy_z_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_z_y_xx[i] = -2.0 * g_0_z_yz_xx[i] * c_exps[i] + 4.0 * g_yy_z_yz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_y_xy[i] = -2.0 * g_0_z_yz_xy[i] * c_exps[i] + 4.0 * g_yy_z_yz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_y_xz[i] = -2.0 * g_0_z_yz_xz[i] * c_exps[i] + 4.0 * g_yy_z_yz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_y_yy[i] = -2.0 * g_0_z_yz_yy[i] * c_exps[i] + 4.0 * g_yy_z_yz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_y_yz[i] = -2.0 * g_0_z_yz_yz[i] * c_exps[i] + 4.0 * g_yy_z_yz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_y_zz[i] = -2.0 * g_0_z_yz_zz[i] * c_exps[i] + 4.0 * g_yy_z_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (912-918)

    #pragma omp simd aligned(g_0_z_0_xx, g_0_z_0_xy, g_0_z_0_xz, g_0_z_0_yy, g_0_z_0_yz, g_0_z_0_zz, g_0_z_zz_xx, g_0_z_zz_xy, g_0_z_zz_xz, g_0_z_zz_yy, g_0_z_zz_yz, g_0_z_zz_zz, g_y_0_z_0_y_z_z_xx, g_y_0_z_0_y_z_z_xy, g_y_0_z_0_y_z_z_xz, g_y_0_z_0_y_z_z_yy, g_y_0_z_0_y_z_z_yz, g_y_0_z_0_y_z_z_zz, g_yy_z_0_xx, g_yy_z_0_xy, g_yy_z_0_xz, g_yy_z_0_yy, g_yy_z_0_yz, g_yy_z_0_zz, g_yy_z_zz_xx, g_yy_z_zz_xy, g_yy_z_zz_xz, g_yy_z_zz_yy, g_yy_z_zz_yz, g_yy_z_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_z_z_xx[i] = g_0_z_0_xx[i] - 2.0 * g_0_z_zz_xx[i] * c_exps[i] - 2.0 * g_yy_z_0_xx[i] * a_exp + 4.0 * g_yy_z_zz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_z_xy[i] = g_0_z_0_xy[i] - 2.0 * g_0_z_zz_xy[i] * c_exps[i] - 2.0 * g_yy_z_0_xy[i] * a_exp + 4.0 * g_yy_z_zz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_z_xz[i] = g_0_z_0_xz[i] - 2.0 * g_0_z_zz_xz[i] * c_exps[i] - 2.0 * g_yy_z_0_xz[i] * a_exp + 4.0 * g_yy_z_zz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_z_yy[i] = g_0_z_0_yy[i] - 2.0 * g_0_z_zz_yy[i] * c_exps[i] - 2.0 * g_yy_z_0_yy[i] * a_exp + 4.0 * g_yy_z_zz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_z_yz[i] = g_0_z_0_yz[i] - 2.0 * g_0_z_zz_yz[i] * c_exps[i] - 2.0 * g_yy_z_0_yz[i] * a_exp + 4.0 * g_yy_z_zz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_z_zz[i] = g_0_z_0_zz[i] - 2.0 * g_0_z_zz_zz[i] * c_exps[i] - 2.0 * g_yy_z_0_zz[i] * a_exp + 4.0 * g_yy_z_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (918-924)

    #pragma omp simd aligned(g_y_0_z_0_z_x_x_xx, g_y_0_z_0_z_x_x_xy, g_y_0_z_0_z_x_x_xz, g_y_0_z_0_z_x_x_yy, g_y_0_z_0_z_x_x_yz, g_y_0_z_0_z_x_x_zz, g_yz_x_xz_xx, g_yz_x_xz_xy, g_yz_x_xz_xz, g_yz_x_xz_yy, g_yz_x_xz_yz, g_yz_x_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_x_x_xx[i] = 4.0 * g_yz_x_xz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_x_xy[i] = 4.0 * g_yz_x_xz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_x_xz[i] = 4.0 * g_yz_x_xz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_x_yy[i] = 4.0 * g_yz_x_xz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_x_yz[i] = 4.0 * g_yz_x_xz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_x_zz[i] = 4.0 * g_yz_x_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (924-930)

    #pragma omp simd aligned(g_y_0_z_0_z_x_y_xx, g_y_0_z_0_z_x_y_xy, g_y_0_z_0_z_x_y_xz, g_y_0_z_0_z_x_y_yy, g_y_0_z_0_z_x_y_yz, g_y_0_z_0_z_x_y_zz, g_yz_x_yz_xx, g_yz_x_yz_xy, g_yz_x_yz_xz, g_yz_x_yz_yy, g_yz_x_yz_yz, g_yz_x_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_x_y_xx[i] = 4.0 * g_yz_x_yz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_y_xy[i] = 4.0 * g_yz_x_yz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_y_xz[i] = 4.0 * g_yz_x_yz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_y_yy[i] = 4.0 * g_yz_x_yz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_y_yz[i] = 4.0 * g_yz_x_yz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_y_zz[i] = 4.0 * g_yz_x_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (930-936)

    #pragma omp simd aligned(g_y_0_z_0_z_x_z_xx, g_y_0_z_0_z_x_z_xy, g_y_0_z_0_z_x_z_xz, g_y_0_z_0_z_x_z_yy, g_y_0_z_0_z_x_z_yz, g_y_0_z_0_z_x_z_zz, g_yz_x_0_xx, g_yz_x_0_xy, g_yz_x_0_xz, g_yz_x_0_yy, g_yz_x_0_yz, g_yz_x_0_zz, g_yz_x_zz_xx, g_yz_x_zz_xy, g_yz_x_zz_xz, g_yz_x_zz_yy, g_yz_x_zz_yz, g_yz_x_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_x_z_xx[i] = -2.0 * g_yz_x_0_xx[i] * a_exp + 4.0 * g_yz_x_zz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_z_xy[i] = -2.0 * g_yz_x_0_xy[i] * a_exp + 4.0 * g_yz_x_zz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_z_xz[i] = -2.0 * g_yz_x_0_xz[i] * a_exp + 4.0 * g_yz_x_zz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_z_yy[i] = -2.0 * g_yz_x_0_yy[i] * a_exp + 4.0 * g_yz_x_zz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_z_yz[i] = -2.0 * g_yz_x_0_yz[i] * a_exp + 4.0 * g_yz_x_zz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_z_zz[i] = -2.0 * g_yz_x_0_zz[i] * a_exp + 4.0 * g_yz_x_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (936-942)

    #pragma omp simd aligned(g_y_0_z_0_z_y_x_xx, g_y_0_z_0_z_y_x_xy, g_y_0_z_0_z_y_x_xz, g_y_0_z_0_z_y_x_yy, g_y_0_z_0_z_y_x_yz, g_y_0_z_0_z_y_x_zz, g_yz_y_xz_xx, g_yz_y_xz_xy, g_yz_y_xz_xz, g_yz_y_xz_yy, g_yz_y_xz_yz, g_yz_y_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_y_x_xx[i] = 4.0 * g_yz_y_xz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_x_xy[i] = 4.0 * g_yz_y_xz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_x_xz[i] = 4.0 * g_yz_y_xz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_x_yy[i] = 4.0 * g_yz_y_xz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_x_yz[i] = 4.0 * g_yz_y_xz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_x_zz[i] = 4.0 * g_yz_y_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (942-948)

    #pragma omp simd aligned(g_y_0_z_0_z_y_y_xx, g_y_0_z_0_z_y_y_xy, g_y_0_z_0_z_y_y_xz, g_y_0_z_0_z_y_y_yy, g_y_0_z_0_z_y_y_yz, g_y_0_z_0_z_y_y_zz, g_yz_y_yz_xx, g_yz_y_yz_xy, g_yz_y_yz_xz, g_yz_y_yz_yy, g_yz_y_yz_yz, g_yz_y_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_y_y_xx[i] = 4.0 * g_yz_y_yz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_y_xy[i] = 4.0 * g_yz_y_yz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_y_xz[i] = 4.0 * g_yz_y_yz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_y_yy[i] = 4.0 * g_yz_y_yz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_y_yz[i] = 4.0 * g_yz_y_yz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_y_zz[i] = 4.0 * g_yz_y_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (948-954)

    #pragma omp simd aligned(g_y_0_z_0_z_y_z_xx, g_y_0_z_0_z_y_z_xy, g_y_0_z_0_z_y_z_xz, g_y_0_z_0_z_y_z_yy, g_y_0_z_0_z_y_z_yz, g_y_0_z_0_z_y_z_zz, g_yz_y_0_xx, g_yz_y_0_xy, g_yz_y_0_xz, g_yz_y_0_yy, g_yz_y_0_yz, g_yz_y_0_zz, g_yz_y_zz_xx, g_yz_y_zz_xy, g_yz_y_zz_xz, g_yz_y_zz_yy, g_yz_y_zz_yz, g_yz_y_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_y_z_xx[i] = -2.0 * g_yz_y_0_xx[i] * a_exp + 4.0 * g_yz_y_zz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_z_xy[i] = -2.0 * g_yz_y_0_xy[i] * a_exp + 4.0 * g_yz_y_zz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_z_xz[i] = -2.0 * g_yz_y_0_xz[i] * a_exp + 4.0 * g_yz_y_zz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_z_yy[i] = -2.0 * g_yz_y_0_yy[i] * a_exp + 4.0 * g_yz_y_zz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_z_yz[i] = -2.0 * g_yz_y_0_yz[i] * a_exp + 4.0 * g_yz_y_zz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_z_zz[i] = -2.0 * g_yz_y_0_zz[i] * a_exp + 4.0 * g_yz_y_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (954-960)

    #pragma omp simd aligned(g_y_0_z_0_z_z_x_xx, g_y_0_z_0_z_z_x_xy, g_y_0_z_0_z_z_x_xz, g_y_0_z_0_z_z_x_yy, g_y_0_z_0_z_z_x_yz, g_y_0_z_0_z_z_x_zz, g_yz_z_xz_xx, g_yz_z_xz_xy, g_yz_z_xz_xz, g_yz_z_xz_yy, g_yz_z_xz_yz, g_yz_z_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_z_x_xx[i] = 4.0 * g_yz_z_xz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_x_xy[i] = 4.0 * g_yz_z_xz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_x_xz[i] = 4.0 * g_yz_z_xz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_x_yy[i] = 4.0 * g_yz_z_xz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_x_yz[i] = 4.0 * g_yz_z_xz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_x_zz[i] = 4.0 * g_yz_z_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (960-966)

    #pragma omp simd aligned(g_y_0_z_0_z_z_y_xx, g_y_0_z_0_z_z_y_xy, g_y_0_z_0_z_z_y_xz, g_y_0_z_0_z_z_y_yy, g_y_0_z_0_z_z_y_yz, g_y_0_z_0_z_z_y_zz, g_yz_z_yz_xx, g_yz_z_yz_xy, g_yz_z_yz_xz, g_yz_z_yz_yy, g_yz_z_yz_yz, g_yz_z_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_z_y_xx[i] = 4.0 * g_yz_z_yz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_y_xy[i] = 4.0 * g_yz_z_yz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_y_xz[i] = 4.0 * g_yz_z_yz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_y_yy[i] = 4.0 * g_yz_z_yz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_y_yz[i] = 4.0 * g_yz_z_yz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_y_zz[i] = 4.0 * g_yz_z_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (966-972)

    #pragma omp simd aligned(g_y_0_z_0_z_z_z_xx, g_y_0_z_0_z_z_z_xy, g_y_0_z_0_z_z_z_xz, g_y_0_z_0_z_z_z_yy, g_y_0_z_0_z_z_z_yz, g_y_0_z_0_z_z_z_zz, g_yz_z_0_xx, g_yz_z_0_xy, g_yz_z_0_xz, g_yz_z_0_yy, g_yz_z_0_yz, g_yz_z_0_zz, g_yz_z_zz_xx, g_yz_z_zz_xy, g_yz_z_zz_xz, g_yz_z_zz_yy, g_yz_z_zz_yz, g_yz_z_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_z_z_xx[i] = -2.0 * g_yz_z_0_xx[i] * a_exp + 4.0 * g_yz_z_zz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_z_xy[i] = -2.0 * g_yz_z_0_xy[i] * a_exp + 4.0 * g_yz_z_zz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_z_xz[i] = -2.0 * g_yz_z_0_xz[i] * a_exp + 4.0 * g_yz_z_zz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_z_yy[i] = -2.0 * g_yz_z_0_yy[i] * a_exp + 4.0 * g_yz_z_zz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_z_yz[i] = -2.0 * g_yz_z_0_yz[i] * a_exp + 4.0 * g_yz_z_zz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_z_zz[i] = -2.0 * g_yz_z_0_zz[i] * a_exp + 4.0 * g_yz_z_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (972-978)

    #pragma omp simd aligned(g_xz_x_0_xx, g_xz_x_0_xy, g_xz_x_0_xz, g_xz_x_0_yy, g_xz_x_0_yz, g_xz_x_0_zz, g_xz_x_xx_xx, g_xz_x_xx_xy, g_xz_x_xx_xz, g_xz_x_xx_yy, g_xz_x_xx_yz, g_xz_x_xx_zz, g_z_0_x_0_x_x_x_xx, g_z_0_x_0_x_x_x_xy, g_z_0_x_0_x_x_x_xz, g_z_0_x_0_x_x_x_yy, g_z_0_x_0_x_x_x_yz, g_z_0_x_0_x_x_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_x_x_xx[i] = -2.0 * g_xz_x_0_xx[i] * a_exp + 4.0 * g_xz_x_xx_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_x_xy[i] = -2.0 * g_xz_x_0_xy[i] * a_exp + 4.0 * g_xz_x_xx_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_x_xz[i] = -2.0 * g_xz_x_0_xz[i] * a_exp + 4.0 * g_xz_x_xx_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_x_yy[i] = -2.0 * g_xz_x_0_yy[i] * a_exp + 4.0 * g_xz_x_xx_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_x_yz[i] = -2.0 * g_xz_x_0_yz[i] * a_exp + 4.0 * g_xz_x_xx_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_x_zz[i] = -2.0 * g_xz_x_0_zz[i] * a_exp + 4.0 * g_xz_x_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (978-984)

    #pragma omp simd aligned(g_xz_x_xy_xx, g_xz_x_xy_xy, g_xz_x_xy_xz, g_xz_x_xy_yy, g_xz_x_xy_yz, g_xz_x_xy_zz, g_z_0_x_0_x_x_y_xx, g_z_0_x_0_x_x_y_xy, g_z_0_x_0_x_x_y_xz, g_z_0_x_0_x_x_y_yy, g_z_0_x_0_x_x_y_yz, g_z_0_x_0_x_x_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_x_y_xx[i] = 4.0 * g_xz_x_xy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_y_xy[i] = 4.0 * g_xz_x_xy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_y_xz[i] = 4.0 * g_xz_x_xy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_y_yy[i] = 4.0 * g_xz_x_xy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_y_yz[i] = 4.0 * g_xz_x_xy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_y_zz[i] = 4.0 * g_xz_x_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (984-990)

    #pragma omp simd aligned(g_xz_x_xz_xx, g_xz_x_xz_xy, g_xz_x_xz_xz, g_xz_x_xz_yy, g_xz_x_xz_yz, g_xz_x_xz_zz, g_z_0_x_0_x_x_z_xx, g_z_0_x_0_x_x_z_xy, g_z_0_x_0_x_x_z_xz, g_z_0_x_0_x_x_z_yy, g_z_0_x_0_x_x_z_yz, g_z_0_x_0_x_x_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_x_z_xx[i] = 4.0 * g_xz_x_xz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_z_xy[i] = 4.0 * g_xz_x_xz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_z_xz[i] = 4.0 * g_xz_x_xz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_z_yy[i] = 4.0 * g_xz_x_xz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_z_yz[i] = 4.0 * g_xz_x_xz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_z_zz[i] = 4.0 * g_xz_x_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (990-996)

    #pragma omp simd aligned(g_xz_y_0_xx, g_xz_y_0_xy, g_xz_y_0_xz, g_xz_y_0_yy, g_xz_y_0_yz, g_xz_y_0_zz, g_xz_y_xx_xx, g_xz_y_xx_xy, g_xz_y_xx_xz, g_xz_y_xx_yy, g_xz_y_xx_yz, g_xz_y_xx_zz, g_z_0_x_0_x_y_x_xx, g_z_0_x_0_x_y_x_xy, g_z_0_x_0_x_y_x_xz, g_z_0_x_0_x_y_x_yy, g_z_0_x_0_x_y_x_yz, g_z_0_x_0_x_y_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_y_x_xx[i] = -2.0 * g_xz_y_0_xx[i] * a_exp + 4.0 * g_xz_y_xx_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_x_xy[i] = -2.0 * g_xz_y_0_xy[i] * a_exp + 4.0 * g_xz_y_xx_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_x_xz[i] = -2.0 * g_xz_y_0_xz[i] * a_exp + 4.0 * g_xz_y_xx_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_x_yy[i] = -2.0 * g_xz_y_0_yy[i] * a_exp + 4.0 * g_xz_y_xx_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_x_yz[i] = -2.0 * g_xz_y_0_yz[i] * a_exp + 4.0 * g_xz_y_xx_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_x_zz[i] = -2.0 * g_xz_y_0_zz[i] * a_exp + 4.0 * g_xz_y_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (996-1002)

    #pragma omp simd aligned(g_xz_y_xy_xx, g_xz_y_xy_xy, g_xz_y_xy_xz, g_xz_y_xy_yy, g_xz_y_xy_yz, g_xz_y_xy_zz, g_z_0_x_0_x_y_y_xx, g_z_0_x_0_x_y_y_xy, g_z_0_x_0_x_y_y_xz, g_z_0_x_0_x_y_y_yy, g_z_0_x_0_x_y_y_yz, g_z_0_x_0_x_y_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_y_y_xx[i] = 4.0 * g_xz_y_xy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_y_xy[i] = 4.0 * g_xz_y_xy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_y_xz[i] = 4.0 * g_xz_y_xy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_y_yy[i] = 4.0 * g_xz_y_xy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_y_yz[i] = 4.0 * g_xz_y_xy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_y_zz[i] = 4.0 * g_xz_y_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1002-1008)

    #pragma omp simd aligned(g_xz_y_xz_xx, g_xz_y_xz_xy, g_xz_y_xz_xz, g_xz_y_xz_yy, g_xz_y_xz_yz, g_xz_y_xz_zz, g_z_0_x_0_x_y_z_xx, g_z_0_x_0_x_y_z_xy, g_z_0_x_0_x_y_z_xz, g_z_0_x_0_x_y_z_yy, g_z_0_x_0_x_y_z_yz, g_z_0_x_0_x_y_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_y_z_xx[i] = 4.0 * g_xz_y_xz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_z_xy[i] = 4.0 * g_xz_y_xz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_z_xz[i] = 4.0 * g_xz_y_xz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_z_yy[i] = 4.0 * g_xz_y_xz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_z_yz[i] = 4.0 * g_xz_y_xz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_z_zz[i] = 4.0 * g_xz_y_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1008-1014)

    #pragma omp simd aligned(g_xz_z_0_xx, g_xz_z_0_xy, g_xz_z_0_xz, g_xz_z_0_yy, g_xz_z_0_yz, g_xz_z_0_zz, g_xz_z_xx_xx, g_xz_z_xx_xy, g_xz_z_xx_xz, g_xz_z_xx_yy, g_xz_z_xx_yz, g_xz_z_xx_zz, g_z_0_x_0_x_z_x_xx, g_z_0_x_0_x_z_x_xy, g_z_0_x_0_x_z_x_xz, g_z_0_x_0_x_z_x_yy, g_z_0_x_0_x_z_x_yz, g_z_0_x_0_x_z_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_z_x_xx[i] = -2.0 * g_xz_z_0_xx[i] * a_exp + 4.0 * g_xz_z_xx_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_x_xy[i] = -2.0 * g_xz_z_0_xy[i] * a_exp + 4.0 * g_xz_z_xx_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_x_xz[i] = -2.0 * g_xz_z_0_xz[i] * a_exp + 4.0 * g_xz_z_xx_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_x_yy[i] = -2.0 * g_xz_z_0_yy[i] * a_exp + 4.0 * g_xz_z_xx_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_x_yz[i] = -2.0 * g_xz_z_0_yz[i] * a_exp + 4.0 * g_xz_z_xx_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_x_zz[i] = -2.0 * g_xz_z_0_zz[i] * a_exp + 4.0 * g_xz_z_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1014-1020)

    #pragma omp simd aligned(g_xz_z_xy_xx, g_xz_z_xy_xy, g_xz_z_xy_xz, g_xz_z_xy_yy, g_xz_z_xy_yz, g_xz_z_xy_zz, g_z_0_x_0_x_z_y_xx, g_z_0_x_0_x_z_y_xy, g_z_0_x_0_x_z_y_xz, g_z_0_x_0_x_z_y_yy, g_z_0_x_0_x_z_y_yz, g_z_0_x_0_x_z_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_z_y_xx[i] = 4.0 * g_xz_z_xy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_y_xy[i] = 4.0 * g_xz_z_xy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_y_xz[i] = 4.0 * g_xz_z_xy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_y_yy[i] = 4.0 * g_xz_z_xy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_y_yz[i] = 4.0 * g_xz_z_xy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_y_zz[i] = 4.0 * g_xz_z_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1020-1026)

    #pragma omp simd aligned(g_xz_z_xz_xx, g_xz_z_xz_xy, g_xz_z_xz_xz, g_xz_z_xz_yy, g_xz_z_xz_yz, g_xz_z_xz_zz, g_z_0_x_0_x_z_z_xx, g_z_0_x_0_x_z_z_xy, g_z_0_x_0_x_z_z_xz, g_z_0_x_0_x_z_z_yy, g_z_0_x_0_x_z_z_yz, g_z_0_x_0_x_z_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_z_z_xx[i] = 4.0 * g_xz_z_xz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_z_xy[i] = 4.0 * g_xz_z_xz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_z_xz[i] = 4.0 * g_xz_z_xz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_z_yy[i] = 4.0 * g_xz_z_xz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_z_yz[i] = 4.0 * g_xz_z_xz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_z_zz[i] = 4.0 * g_xz_z_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1026-1032)

    #pragma omp simd aligned(g_yz_x_0_xx, g_yz_x_0_xy, g_yz_x_0_xz, g_yz_x_0_yy, g_yz_x_0_yz, g_yz_x_0_zz, g_yz_x_xx_xx, g_yz_x_xx_xy, g_yz_x_xx_xz, g_yz_x_xx_yy, g_yz_x_xx_yz, g_yz_x_xx_zz, g_z_0_x_0_y_x_x_xx, g_z_0_x_0_y_x_x_xy, g_z_0_x_0_y_x_x_xz, g_z_0_x_0_y_x_x_yy, g_z_0_x_0_y_x_x_yz, g_z_0_x_0_y_x_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_x_x_xx[i] = -2.0 * g_yz_x_0_xx[i] * a_exp + 4.0 * g_yz_x_xx_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_x_xy[i] = -2.0 * g_yz_x_0_xy[i] * a_exp + 4.0 * g_yz_x_xx_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_x_xz[i] = -2.0 * g_yz_x_0_xz[i] * a_exp + 4.0 * g_yz_x_xx_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_x_yy[i] = -2.0 * g_yz_x_0_yy[i] * a_exp + 4.0 * g_yz_x_xx_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_x_yz[i] = -2.0 * g_yz_x_0_yz[i] * a_exp + 4.0 * g_yz_x_xx_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_x_zz[i] = -2.0 * g_yz_x_0_zz[i] * a_exp + 4.0 * g_yz_x_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1032-1038)

    #pragma omp simd aligned(g_yz_x_xy_xx, g_yz_x_xy_xy, g_yz_x_xy_xz, g_yz_x_xy_yy, g_yz_x_xy_yz, g_yz_x_xy_zz, g_z_0_x_0_y_x_y_xx, g_z_0_x_0_y_x_y_xy, g_z_0_x_0_y_x_y_xz, g_z_0_x_0_y_x_y_yy, g_z_0_x_0_y_x_y_yz, g_z_0_x_0_y_x_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_x_y_xx[i] = 4.0 * g_yz_x_xy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_y_xy[i] = 4.0 * g_yz_x_xy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_y_xz[i] = 4.0 * g_yz_x_xy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_y_yy[i] = 4.0 * g_yz_x_xy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_y_yz[i] = 4.0 * g_yz_x_xy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_y_zz[i] = 4.0 * g_yz_x_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1038-1044)

    #pragma omp simd aligned(g_yz_x_xz_xx, g_yz_x_xz_xy, g_yz_x_xz_xz, g_yz_x_xz_yy, g_yz_x_xz_yz, g_yz_x_xz_zz, g_z_0_x_0_y_x_z_xx, g_z_0_x_0_y_x_z_xy, g_z_0_x_0_y_x_z_xz, g_z_0_x_0_y_x_z_yy, g_z_0_x_0_y_x_z_yz, g_z_0_x_0_y_x_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_x_z_xx[i] = 4.0 * g_yz_x_xz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_z_xy[i] = 4.0 * g_yz_x_xz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_z_xz[i] = 4.0 * g_yz_x_xz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_z_yy[i] = 4.0 * g_yz_x_xz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_z_yz[i] = 4.0 * g_yz_x_xz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_z_zz[i] = 4.0 * g_yz_x_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1044-1050)

    #pragma omp simd aligned(g_yz_y_0_xx, g_yz_y_0_xy, g_yz_y_0_xz, g_yz_y_0_yy, g_yz_y_0_yz, g_yz_y_0_zz, g_yz_y_xx_xx, g_yz_y_xx_xy, g_yz_y_xx_xz, g_yz_y_xx_yy, g_yz_y_xx_yz, g_yz_y_xx_zz, g_z_0_x_0_y_y_x_xx, g_z_0_x_0_y_y_x_xy, g_z_0_x_0_y_y_x_xz, g_z_0_x_0_y_y_x_yy, g_z_0_x_0_y_y_x_yz, g_z_0_x_0_y_y_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_y_x_xx[i] = -2.0 * g_yz_y_0_xx[i] * a_exp + 4.0 * g_yz_y_xx_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_x_xy[i] = -2.0 * g_yz_y_0_xy[i] * a_exp + 4.0 * g_yz_y_xx_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_x_xz[i] = -2.0 * g_yz_y_0_xz[i] * a_exp + 4.0 * g_yz_y_xx_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_x_yy[i] = -2.0 * g_yz_y_0_yy[i] * a_exp + 4.0 * g_yz_y_xx_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_x_yz[i] = -2.0 * g_yz_y_0_yz[i] * a_exp + 4.0 * g_yz_y_xx_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_x_zz[i] = -2.0 * g_yz_y_0_zz[i] * a_exp + 4.0 * g_yz_y_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1050-1056)

    #pragma omp simd aligned(g_yz_y_xy_xx, g_yz_y_xy_xy, g_yz_y_xy_xz, g_yz_y_xy_yy, g_yz_y_xy_yz, g_yz_y_xy_zz, g_z_0_x_0_y_y_y_xx, g_z_0_x_0_y_y_y_xy, g_z_0_x_0_y_y_y_xz, g_z_0_x_0_y_y_y_yy, g_z_0_x_0_y_y_y_yz, g_z_0_x_0_y_y_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_y_y_xx[i] = 4.0 * g_yz_y_xy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_y_xy[i] = 4.0 * g_yz_y_xy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_y_xz[i] = 4.0 * g_yz_y_xy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_y_yy[i] = 4.0 * g_yz_y_xy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_y_yz[i] = 4.0 * g_yz_y_xy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_y_zz[i] = 4.0 * g_yz_y_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1056-1062)

    #pragma omp simd aligned(g_yz_y_xz_xx, g_yz_y_xz_xy, g_yz_y_xz_xz, g_yz_y_xz_yy, g_yz_y_xz_yz, g_yz_y_xz_zz, g_z_0_x_0_y_y_z_xx, g_z_0_x_0_y_y_z_xy, g_z_0_x_0_y_y_z_xz, g_z_0_x_0_y_y_z_yy, g_z_0_x_0_y_y_z_yz, g_z_0_x_0_y_y_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_y_z_xx[i] = 4.0 * g_yz_y_xz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_z_xy[i] = 4.0 * g_yz_y_xz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_z_xz[i] = 4.0 * g_yz_y_xz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_z_yy[i] = 4.0 * g_yz_y_xz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_z_yz[i] = 4.0 * g_yz_y_xz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_z_zz[i] = 4.0 * g_yz_y_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1062-1068)

    #pragma omp simd aligned(g_yz_z_0_xx, g_yz_z_0_xy, g_yz_z_0_xz, g_yz_z_0_yy, g_yz_z_0_yz, g_yz_z_0_zz, g_yz_z_xx_xx, g_yz_z_xx_xy, g_yz_z_xx_xz, g_yz_z_xx_yy, g_yz_z_xx_yz, g_yz_z_xx_zz, g_z_0_x_0_y_z_x_xx, g_z_0_x_0_y_z_x_xy, g_z_0_x_0_y_z_x_xz, g_z_0_x_0_y_z_x_yy, g_z_0_x_0_y_z_x_yz, g_z_0_x_0_y_z_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_z_x_xx[i] = -2.0 * g_yz_z_0_xx[i] * a_exp + 4.0 * g_yz_z_xx_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_x_xy[i] = -2.0 * g_yz_z_0_xy[i] * a_exp + 4.0 * g_yz_z_xx_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_x_xz[i] = -2.0 * g_yz_z_0_xz[i] * a_exp + 4.0 * g_yz_z_xx_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_x_yy[i] = -2.0 * g_yz_z_0_yy[i] * a_exp + 4.0 * g_yz_z_xx_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_x_yz[i] = -2.0 * g_yz_z_0_yz[i] * a_exp + 4.0 * g_yz_z_xx_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_x_zz[i] = -2.0 * g_yz_z_0_zz[i] * a_exp + 4.0 * g_yz_z_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1068-1074)

    #pragma omp simd aligned(g_yz_z_xy_xx, g_yz_z_xy_xy, g_yz_z_xy_xz, g_yz_z_xy_yy, g_yz_z_xy_yz, g_yz_z_xy_zz, g_z_0_x_0_y_z_y_xx, g_z_0_x_0_y_z_y_xy, g_z_0_x_0_y_z_y_xz, g_z_0_x_0_y_z_y_yy, g_z_0_x_0_y_z_y_yz, g_z_0_x_0_y_z_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_z_y_xx[i] = 4.0 * g_yz_z_xy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_y_xy[i] = 4.0 * g_yz_z_xy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_y_xz[i] = 4.0 * g_yz_z_xy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_y_yy[i] = 4.0 * g_yz_z_xy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_y_yz[i] = 4.0 * g_yz_z_xy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_y_zz[i] = 4.0 * g_yz_z_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1074-1080)

    #pragma omp simd aligned(g_yz_z_xz_xx, g_yz_z_xz_xy, g_yz_z_xz_xz, g_yz_z_xz_yy, g_yz_z_xz_yz, g_yz_z_xz_zz, g_z_0_x_0_y_z_z_xx, g_z_0_x_0_y_z_z_xy, g_z_0_x_0_y_z_z_xz, g_z_0_x_0_y_z_z_yy, g_z_0_x_0_y_z_z_yz, g_z_0_x_0_y_z_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_z_z_xx[i] = 4.0 * g_yz_z_xz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_z_xy[i] = 4.0 * g_yz_z_xz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_z_xz[i] = 4.0 * g_yz_z_xz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_z_yy[i] = 4.0 * g_yz_z_xz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_z_yz[i] = 4.0 * g_yz_z_xz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_z_zz[i] = 4.0 * g_yz_z_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1080-1086)

    #pragma omp simd aligned(g_0_x_0_xx, g_0_x_0_xy, g_0_x_0_xz, g_0_x_0_yy, g_0_x_0_yz, g_0_x_0_zz, g_0_x_xx_xx, g_0_x_xx_xy, g_0_x_xx_xz, g_0_x_xx_yy, g_0_x_xx_yz, g_0_x_xx_zz, g_z_0_x_0_z_x_x_xx, g_z_0_x_0_z_x_x_xy, g_z_0_x_0_z_x_x_xz, g_z_0_x_0_z_x_x_yy, g_z_0_x_0_z_x_x_yz, g_z_0_x_0_z_x_x_zz, g_zz_x_0_xx, g_zz_x_0_xy, g_zz_x_0_xz, g_zz_x_0_yy, g_zz_x_0_yz, g_zz_x_0_zz, g_zz_x_xx_xx, g_zz_x_xx_xy, g_zz_x_xx_xz, g_zz_x_xx_yy, g_zz_x_xx_yz, g_zz_x_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_x_x_xx[i] = g_0_x_0_xx[i] - 2.0 * g_0_x_xx_xx[i] * c_exps[i] - 2.0 * g_zz_x_0_xx[i] * a_exp + 4.0 * g_zz_x_xx_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_x_xy[i] = g_0_x_0_xy[i] - 2.0 * g_0_x_xx_xy[i] * c_exps[i] - 2.0 * g_zz_x_0_xy[i] * a_exp + 4.0 * g_zz_x_xx_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_x_xz[i] = g_0_x_0_xz[i] - 2.0 * g_0_x_xx_xz[i] * c_exps[i] - 2.0 * g_zz_x_0_xz[i] * a_exp + 4.0 * g_zz_x_xx_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_x_yy[i] = g_0_x_0_yy[i] - 2.0 * g_0_x_xx_yy[i] * c_exps[i] - 2.0 * g_zz_x_0_yy[i] * a_exp + 4.0 * g_zz_x_xx_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_x_yz[i] = g_0_x_0_yz[i] - 2.0 * g_0_x_xx_yz[i] * c_exps[i] - 2.0 * g_zz_x_0_yz[i] * a_exp + 4.0 * g_zz_x_xx_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_x_zz[i] = g_0_x_0_zz[i] - 2.0 * g_0_x_xx_zz[i] * c_exps[i] - 2.0 * g_zz_x_0_zz[i] * a_exp + 4.0 * g_zz_x_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1086-1092)

    #pragma omp simd aligned(g_0_x_xy_xx, g_0_x_xy_xy, g_0_x_xy_xz, g_0_x_xy_yy, g_0_x_xy_yz, g_0_x_xy_zz, g_z_0_x_0_z_x_y_xx, g_z_0_x_0_z_x_y_xy, g_z_0_x_0_z_x_y_xz, g_z_0_x_0_z_x_y_yy, g_z_0_x_0_z_x_y_yz, g_z_0_x_0_z_x_y_zz, g_zz_x_xy_xx, g_zz_x_xy_xy, g_zz_x_xy_xz, g_zz_x_xy_yy, g_zz_x_xy_yz, g_zz_x_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_x_y_xx[i] = -2.0 * g_0_x_xy_xx[i] * c_exps[i] + 4.0 * g_zz_x_xy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_y_xy[i] = -2.0 * g_0_x_xy_xy[i] * c_exps[i] + 4.0 * g_zz_x_xy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_y_xz[i] = -2.0 * g_0_x_xy_xz[i] * c_exps[i] + 4.0 * g_zz_x_xy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_y_yy[i] = -2.0 * g_0_x_xy_yy[i] * c_exps[i] + 4.0 * g_zz_x_xy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_y_yz[i] = -2.0 * g_0_x_xy_yz[i] * c_exps[i] + 4.0 * g_zz_x_xy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_y_zz[i] = -2.0 * g_0_x_xy_zz[i] * c_exps[i] + 4.0 * g_zz_x_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1092-1098)

    #pragma omp simd aligned(g_0_x_xz_xx, g_0_x_xz_xy, g_0_x_xz_xz, g_0_x_xz_yy, g_0_x_xz_yz, g_0_x_xz_zz, g_z_0_x_0_z_x_z_xx, g_z_0_x_0_z_x_z_xy, g_z_0_x_0_z_x_z_xz, g_z_0_x_0_z_x_z_yy, g_z_0_x_0_z_x_z_yz, g_z_0_x_0_z_x_z_zz, g_zz_x_xz_xx, g_zz_x_xz_xy, g_zz_x_xz_xz, g_zz_x_xz_yy, g_zz_x_xz_yz, g_zz_x_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_x_z_xx[i] = -2.0 * g_0_x_xz_xx[i] * c_exps[i] + 4.0 * g_zz_x_xz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_z_xy[i] = -2.0 * g_0_x_xz_xy[i] * c_exps[i] + 4.0 * g_zz_x_xz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_z_xz[i] = -2.0 * g_0_x_xz_xz[i] * c_exps[i] + 4.0 * g_zz_x_xz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_z_yy[i] = -2.0 * g_0_x_xz_yy[i] * c_exps[i] + 4.0 * g_zz_x_xz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_z_yz[i] = -2.0 * g_0_x_xz_yz[i] * c_exps[i] + 4.0 * g_zz_x_xz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_z_zz[i] = -2.0 * g_0_x_xz_zz[i] * c_exps[i] + 4.0 * g_zz_x_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1098-1104)

    #pragma omp simd aligned(g_0_y_0_xx, g_0_y_0_xy, g_0_y_0_xz, g_0_y_0_yy, g_0_y_0_yz, g_0_y_0_zz, g_0_y_xx_xx, g_0_y_xx_xy, g_0_y_xx_xz, g_0_y_xx_yy, g_0_y_xx_yz, g_0_y_xx_zz, g_z_0_x_0_z_y_x_xx, g_z_0_x_0_z_y_x_xy, g_z_0_x_0_z_y_x_xz, g_z_0_x_0_z_y_x_yy, g_z_0_x_0_z_y_x_yz, g_z_0_x_0_z_y_x_zz, g_zz_y_0_xx, g_zz_y_0_xy, g_zz_y_0_xz, g_zz_y_0_yy, g_zz_y_0_yz, g_zz_y_0_zz, g_zz_y_xx_xx, g_zz_y_xx_xy, g_zz_y_xx_xz, g_zz_y_xx_yy, g_zz_y_xx_yz, g_zz_y_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_y_x_xx[i] = g_0_y_0_xx[i] - 2.0 * g_0_y_xx_xx[i] * c_exps[i] - 2.0 * g_zz_y_0_xx[i] * a_exp + 4.0 * g_zz_y_xx_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_x_xy[i] = g_0_y_0_xy[i] - 2.0 * g_0_y_xx_xy[i] * c_exps[i] - 2.0 * g_zz_y_0_xy[i] * a_exp + 4.0 * g_zz_y_xx_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_x_xz[i] = g_0_y_0_xz[i] - 2.0 * g_0_y_xx_xz[i] * c_exps[i] - 2.0 * g_zz_y_0_xz[i] * a_exp + 4.0 * g_zz_y_xx_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_x_yy[i] = g_0_y_0_yy[i] - 2.0 * g_0_y_xx_yy[i] * c_exps[i] - 2.0 * g_zz_y_0_yy[i] * a_exp + 4.0 * g_zz_y_xx_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_x_yz[i] = g_0_y_0_yz[i] - 2.0 * g_0_y_xx_yz[i] * c_exps[i] - 2.0 * g_zz_y_0_yz[i] * a_exp + 4.0 * g_zz_y_xx_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_x_zz[i] = g_0_y_0_zz[i] - 2.0 * g_0_y_xx_zz[i] * c_exps[i] - 2.0 * g_zz_y_0_zz[i] * a_exp + 4.0 * g_zz_y_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1104-1110)

    #pragma omp simd aligned(g_0_y_xy_xx, g_0_y_xy_xy, g_0_y_xy_xz, g_0_y_xy_yy, g_0_y_xy_yz, g_0_y_xy_zz, g_z_0_x_0_z_y_y_xx, g_z_0_x_0_z_y_y_xy, g_z_0_x_0_z_y_y_xz, g_z_0_x_0_z_y_y_yy, g_z_0_x_0_z_y_y_yz, g_z_0_x_0_z_y_y_zz, g_zz_y_xy_xx, g_zz_y_xy_xy, g_zz_y_xy_xz, g_zz_y_xy_yy, g_zz_y_xy_yz, g_zz_y_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_y_y_xx[i] = -2.0 * g_0_y_xy_xx[i] * c_exps[i] + 4.0 * g_zz_y_xy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_y_xy[i] = -2.0 * g_0_y_xy_xy[i] * c_exps[i] + 4.0 * g_zz_y_xy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_y_xz[i] = -2.0 * g_0_y_xy_xz[i] * c_exps[i] + 4.0 * g_zz_y_xy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_y_yy[i] = -2.0 * g_0_y_xy_yy[i] * c_exps[i] + 4.0 * g_zz_y_xy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_y_yz[i] = -2.0 * g_0_y_xy_yz[i] * c_exps[i] + 4.0 * g_zz_y_xy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_y_zz[i] = -2.0 * g_0_y_xy_zz[i] * c_exps[i] + 4.0 * g_zz_y_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1110-1116)

    #pragma omp simd aligned(g_0_y_xz_xx, g_0_y_xz_xy, g_0_y_xz_xz, g_0_y_xz_yy, g_0_y_xz_yz, g_0_y_xz_zz, g_z_0_x_0_z_y_z_xx, g_z_0_x_0_z_y_z_xy, g_z_0_x_0_z_y_z_xz, g_z_0_x_0_z_y_z_yy, g_z_0_x_0_z_y_z_yz, g_z_0_x_0_z_y_z_zz, g_zz_y_xz_xx, g_zz_y_xz_xy, g_zz_y_xz_xz, g_zz_y_xz_yy, g_zz_y_xz_yz, g_zz_y_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_y_z_xx[i] = -2.0 * g_0_y_xz_xx[i] * c_exps[i] + 4.0 * g_zz_y_xz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_z_xy[i] = -2.0 * g_0_y_xz_xy[i] * c_exps[i] + 4.0 * g_zz_y_xz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_z_xz[i] = -2.0 * g_0_y_xz_xz[i] * c_exps[i] + 4.0 * g_zz_y_xz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_z_yy[i] = -2.0 * g_0_y_xz_yy[i] * c_exps[i] + 4.0 * g_zz_y_xz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_z_yz[i] = -2.0 * g_0_y_xz_yz[i] * c_exps[i] + 4.0 * g_zz_y_xz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_z_zz[i] = -2.0 * g_0_y_xz_zz[i] * c_exps[i] + 4.0 * g_zz_y_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1116-1122)

    #pragma omp simd aligned(g_0_z_0_xx, g_0_z_0_xy, g_0_z_0_xz, g_0_z_0_yy, g_0_z_0_yz, g_0_z_0_zz, g_0_z_xx_xx, g_0_z_xx_xy, g_0_z_xx_xz, g_0_z_xx_yy, g_0_z_xx_yz, g_0_z_xx_zz, g_z_0_x_0_z_z_x_xx, g_z_0_x_0_z_z_x_xy, g_z_0_x_0_z_z_x_xz, g_z_0_x_0_z_z_x_yy, g_z_0_x_0_z_z_x_yz, g_z_0_x_0_z_z_x_zz, g_zz_z_0_xx, g_zz_z_0_xy, g_zz_z_0_xz, g_zz_z_0_yy, g_zz_z_0_yz, g_zz_z_0_zz, g_zz_z_xx_xx, g_zz_z_xx_xy, g_zz_z_xx_xz, g_zz_z_xx_yy, g_zz_z_xx_yz, g_zz_z_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_z_x_xx[i] = g_0_z_0_xx[i] - 2.0 * g_0_z_xx_xx[i] * c_exps[i] - 2.0 * g_zz_z_0_xx[i] * a_exp + 4.0 * g_zz_z_xx_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_x_xy[i] = g_0_z_0_xy[i] - 2.0 * g_0_z_xx_xy[i] * c_exps[i] - 2.0 * g_zz_z_0_xy[i] * a_exp + 4.0 * g_zz_z_xx_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_x_xz[i] = g_0_z_0_xz[i] - 2.0 * g_0_z_xx_xz[i] * c_exps[i] - 2.0 * g_zz_z_0_xz[i] * a_exp + 4.0 * g_zz_z_xx_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_x_yy[i] = g_0_z_0_yy[i] - 2.0 * g_0_z_xx_yy[i] * c_exps[i] - 2.0 * g_zz_z_0_yy[i] * a_exp + 4.0 * g_zz_z_xx_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_x_yz[i] = g_0_z_0_yz[i] - 2.0 * g_0_z_xx_yz[i] * c_exps[i] - 2.0 * g_zz_z_0_yz[i] * a_exp + 4.0 * g_zz_z_xx_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_x_zz[i] = g_0_z_0_zz[i] - 2.0 * g_0_z_xx_zz[i] * c_exps[i] - 2.0 * g_zz_z_0_zz[i] * a_exp + 4.0 * g_zz_z_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1122-1128)

    #pragma omp simd aligned(g_0_z_xy_xx, g_0_z_xy_xy, g_0_z_xy_xz, g_0_z_xy_yy, g_0_z_xy_yz, g_0_z_xy_zz, g_z_0_x_0_z_z_y_xx, g_z_0_x_0_z_z_y_xy, g_z_0_x_0_z_z_y_xz, g_z_0_x_0_z_z_y_yy, g_z_0_x_0_z_z_y_yz, g_z_0_x_0_z_z_y_zz, g_zz_z_xy_xx, g_zz_z_xy_xy, g_zz_z_xy_xz, g_zz_z_xy_yy, g_zz_z_xy_yz, g_zz_z_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_z_y_xx[i] = -2.0 * g_0_z_xy_xx[i] * c_exps[i] + 4.0 * g_zz_z_xy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_y_xy[i] = -2.0 * g_0_z_xy_xy[i] * c_exps[i] + 4.0 * g_zz_z_xy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_y_xz[i] = -2.0 * g_0_z_xy_xz[i] * c_exps[i] + 4.0 * g_zz_z_xy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_y_yy[i] = -2.0 * g_0_z_xy_yy[i] * c_exps[i] + 4.0 * g_zz_z_xy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_y_yz[i] = -2.0 * g_0_z_xy_yz[i] * c_exps[i] + 4.0 * g_zz_z_xy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_y_zz[i] = -2.0 * g_0_z_xy_zz[i] * c_exps[i] + 4.0 * g_zz_z_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1128-1134)

    #pragma omp simd aligned(g_0_z_xz_xx, g_0_z_xz_xy, g_0_z_xz_xz, g_0_z_xz_yy, g_0_z_xz_yz, g_0_z_xz_zz, g_z_0_x_0_z_z_z_xx, g_z_0_x_0_z_z_z_xy, g_z_0_x_0_z_z_z_xz, g_z_0_x_0_z_z_z_yy, g_z_0_x_0_z_z_z_yz, g_z_0_x_0_z_z_z_zz, g_zz_z_xz_xx, g_zz_z_xz_xy, g_zz_z_xz_xz, g_zz_z_xz_yy, g_zz_z_xz_yz, g_zz_z_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_z_z_xx[i] = -2.0 * g_0_z_xz_xx[i] * c_exps[i] + 4.0 * g_zz_z_xz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_z_xy[i] = -2.0 * g_0_z_xz_xy[i] * c_exps[i] + 4.0 * g_zz_z_xz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_z_xz[i] = -2.0 * g_0_z_xz_xz[i] * c_exps[i] + 4.0 * g_zz_z_xz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_z_yy[i] = -2.0 * g_0_z_xz_yy[i] * c_exps[i] + 4.0 * g_zz_z_xz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_z_yz[i] = -2.0 * g_0_z_xz_yz[i] * c_exps[i] + 4.0 * g_zz_z_xz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_z_zz[i] = -2.0 * g_0_z_xz_zz[i] * c_exps[i] + 4.0 * g_zz_z_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1134-1140)

    #pragma omp simd aligned(g_xz_x_xy_xx, g_xz_x_xy_xy, g_xz_x_xy_xz, g_xz_x_xy_yy, g_xz_x_xy_yz, g_xz_x_xy_zz, g_z_0_y_0_x_x_x_xx, g_z_0_y_0_x_x_x_xy, g_z_0_y_0_x_x_x_xz, g_z_0_y_0_x_x_x_yy, g_z_0_y_0_x_x_x_yz, g_z_0_y_0_x_x_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_x_x_xx[i] = 4.0 * g_xz_x_xy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_x_xy[i] = 4.0 * g_xz_x_xy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_x_xz[i] = 4.0 * g_xz_x_xy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_x_yy[i] = 4.0 * g_xz_x_xy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_x_yz[i] = 4.0 * g_xz_x_xy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_x_zz[i] = 4.0 * g_xz_x_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1140-1146)

    #pragma omp simd aligned(g_xz_x_0_xx, g_xz_x_0_xy, g_xz_x_0_xz, g_xz_x_0_yy, g_xz_x_0_yz, g_xz_x_0_zz, g_xz_x_yy_xx, g_xz_x_yy_xy, g_xz_x_yy_xz, g_xz_x_yy_yy, g_xz_x_yy_yz, g_xz_x_yy_zz, g_z_0_y_0_x_x_y_xx, g_z_0_y_0_x_x_y_xy, g_z_0_y_0_x_x_y_xz, g_z_0_y_0_x_x_y_yy, g_z_0_y_0_x_x_y_yz, g_z_0_y_0_x_x_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_x_y_xx[i] = -2.0 * g_xz_x_0_xx[i] * a_exp + 4.0 * g_xz_x_yy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_y_xy[i] = -2.0 * g_xz_x_0_xy[i] * a_exp + 4.0 * g_xz_x_yy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_y_xz[i] = -2.0 * g_xz_x_0_xz[i] * a_exp + 4.0 * g_xz_x_yy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_y_yy[i] = -2.0 * g_xz_x_0_yy[i] * a_exp + 4.0 * g_xz_x_yy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_y_yz[i] = -2.0 * g_xz_x_0_yz[i] * a_exp + 4.0 * g_xz_x_yy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_y_zz[i] = -2.0 * g_xz_x_0_zz[i] * a_exp + 4.0 * g_xz_x_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1146-1152)

    #pragma omp simd aligned(g_xz_x_yz_xx, g_xz_x_yz_xy, g_xz_x_yz_xz, g_xz_x_yz_yy, g_xz_x_yz_yz, g_xz_x_yz_zz, g_z_0_y_0_x_x_z_xx, g_z_0_y_0_x_x_z_xy, g_z_0_y_0_x_x_z_xz, g_z_0_y_0_x_x_z_yy, g_z_0_y_0_x_x_z_yz, g_z_0_y_0_x_x_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_x_z_xx[i] = 4.0 * g_xz_x_yz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_z_xy[i] = 4.0 * g_xz_x_yz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_z_xz[i] = 4.0 * g_xz_x_yz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_z_yy[i] = 4.0 * g_xz_x_yz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_z_yz[i] = 4.0 * g_xz_x_yz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_z_zz[i] = 4.0 * g_xz_x_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1152-1158)

    #pragma omp simd aligned(g_xz_y_xy_xx, g_xz_y_xy_xy, g_xz_y_xy_xz, g_xz_y_xy_yy, g_xz_y_xy_yz, g_xz_y_xy_zz, g_z_0_y_0_x_y_x_xx, g_z_0_y_0_x_y_x_xy, g_z_0_y_0_x_y_x_xz, g_z_0_y_0_x_y_x_yy, g_z_0_y_0_x_y_x_yz, g_z_0_y_0_x_y_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_y_x_xx[i] = 4.0 * g_xz_y_xy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_x_xy[i] = 4.0 * g_xz_y_xy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_x_xz[i] = 4.0 * g_xz_y_xy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_x_yy[i] = 4.0 * g_xz_y_xy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_x_yz[i] = 4.0 * g_xz_y_xy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_x_zz[i] = 4.0 * g_xz_y_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1158-1164)

    #pragma omp simd aligned(g_xz_y_0_xx, g_xz_y_0_xy, g_xz_y_0_xz, g_xz_y_0_yy, g_xz_y_0_yz, g_xz_y_0_zz, g_xz_y_yy_xx, g_xz_y_yy_xy, g_xz_y_yy_xz, g_xz_y_yy_yy, g_xz_y_yy_yz, g_xz_y_yy_zz, g_z_0_y_0_x_y_y_xx, g_z_0_y_0_x_y_y_xy, g_z_0_y_0_x_y_y_xz, g_z_0_y_0_x_y_y_yy, g_z_0_y_0_x_y_y_yz, g_z_0_y_0_x_y_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_y_y_xx[i] = -2.0 * g_xz_y_0_xx[i] * a_exp + 4.0 * g_xz_y_yy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_y_xy[i] = -2.0 * g_xz_y_0_xy[i] * a_exp + 4.0 * g_xz_y_yy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_y_xz[i] = -2.0 * g_xz_y_0_xz[i] * a_exp + 4.0 * g_xz_y_yy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_y_yy[i] = -2.0 * g_xz_y_0_yy[i] * a_exp + 4.0 * g_xz_y_yy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_y_yz[i] = -2.0 * g_xz_y_0_yz[i] * a_exp + 4.0 * g_xz_y_yy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_y_zz[i] = -2.0 * g_xz_y_0_zz[i] * a_exp + 4.0 * g_xz_y_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1164-1170)

    #pragma omp simd aligned(g_xz_y_yz_xx, g_xz_y_yz_xy, g_xz_y_yz_xz, g_xz_y_yz_yy, g_xz_y_yz_yz, g_xz_y_yz_zz, g_z_0_y_0_x_y_z_xx, g_z_0_y_0_x_y_z_xy, g_z_0_y_0_x_y_z_xz, g_z_0_y_0_x_y_z_yy, g_z_0_y_0_x_y_z_yz, g_z_0_y_0_x_y_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_y_z_xx[i] = 4.0 * g_xz_y_yz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_z_xy[i] = 4.0 * g_xz_y_yz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_z_xz[i] = 4.0 * g_xz_y_yz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_z_yy[i] = 4.0 * g_xz_y_yz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_z_yz[i] = 4.0 * g_xz_y_yz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_z_zz[i] = 4.0 * g_xz_y_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1170-1176)

    #pragma omp simd aligned(g_xz_z_xy_xx, g_xz_z_xy_xy, g_xz_z_xy_xz, g_xz_z_xy_yy, g_xz_z_xy_yz, g_xz_z_xy_zz, g_z_0_y_0_x_z_x_xx, g_z_0_y_0_x_z_x_xy, g_z_0_y_0_x_z_x_xz, g_z_0_y_0_x_z_x_yy, g_z_0_y_0_x_z_x_yz, g_z_0_y_0_x_z_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_z_x_xx[i] = 4.0 * g_xz_z_xy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_x_xy[i] = 4.0 * g_xz_z_xy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_x_xz[i] = 4.0 * g_xz_z_xy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_x_yy[i] = 4.0 * g_xz_z_xy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_x_yz[i] = 4.0 * g_xz_z_xy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_x_zz[i] = 4.0 * g_xz_z_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1176-1182)

    #pragma omp simd aligned(g_xz_z_0_xx, g_xz_z_0_xy, g_xz_z_0_xz, g_xz_z_0_yy, g_xz_z_0_yz, g_xz_z_0_zz, g_xz_z_yy_xx, g_xz_z_yy_xy, g_xz_z_yy_xz, g_xz_z_yy_yy, g_xz_z_yy_yz, g_xz_z_yy_zz, g_z_0_y_0_x_z_y_xx, g_z_0_y_0_x_z_y_xy, g_z_0_y_0_x_z_y_xz, g_z_0_y_0_x_z_y_yy, g_z_0_y_0_x_z_y_yz, g_z_0_y_0_x_z_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_z_y_xx[i] = -2.0 * g_xz_z_0_xx[i] * a_exp + 4.0 * g_xz_z_yy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_y_xy[i] = -2.0 * g_xz_z_0_xy[i] * a_exp + 4.0 * g_xz_z_yy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_y_xz[i] = -2.0 * g_xz_z_0_xz[i] * a_exp + 4.0 * g_xz_z_yy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_y_yy[i] = -2.0 * g_xz_z_0_yy[i] * a_exp + 4.0 * g_xz_z_yy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_y_yz[i] = -2.0 * g_xz_z_0_yz[i] * a_exp + 4.0 * g_xz_z_yy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_y_zz[i] = -2.0 * g_xz_z_0_zz[i] * a_exp + 4.0 * g_xz_z_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1182-1188)

    #pragma omp simd aligned(g_xz_z_yz_xx, g_xz_z_yz_xy, g_xz_z_yz_xz, g_xz_z_yz_yy, g_xz_z_yz_yz, g_xz_z_yz_zz, g_z_0_y_0_x_z_z_xx, g_z_0_y_0_x_z_z_xy, g_z_0_y_0_x_z_z_xz, g_z_0_y_0_x_z_z_yy, g_z_0_y_0_x_z_z_yz, g_z_0_y_0_x_z_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_z_z_xx[i] = 4.0 * g_xz_z_yz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_z_xy[i] = 4.0 * g_xz_z_yz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_z_xz[i] = 4.0 * g_xz_z_yz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_z_yy[i] = 4.0 * g_xz_z_yz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_z_yz[i] = 4.0 * g_xz_z_yz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_z_zz[i] = 4.0 * g_xz_z_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1188-1194)

    #pragma omp simd aligned(g_yz_x_xy_xx, g_yz_x_xy_xy, g_yz_x_xy_xz, g_yz_x_xy_yy, g_yz_x_xy_yz, g_yz_x_xy_zz, g_z_0_y_0_y_x_x_xx, g_z_0_y_0_y_x_x_xy, g_z_0_y_0_y_x_x_xz, g_z_0_y_0_y_x_x_yy, g_z_0_y_0_y_x_x_yz, g_z_0_y_0_y_x_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_x_x_xx[i] = 4.0 * g_yz_x_xy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_x_xy[i] = 4.0 * g_yz_x_xy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_x_xz[i] = 4.0 * g_yz_x_xy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_x_yy[i] = 4.0 * g_yz_x_xy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_x_yz[i] = 4.0 * g_yz_x_xy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_x_zz[i] = 4.0 * g_yz_x_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1194-1200)

    #pragma omp simd aligned(g_yz_x_0_xx, g_yz_x_0_xy, g_yz_x_0_xz, g_yz_x_0_yy, g_yz_x_0_yz, g_yz_x_0_zz, g_yz_x_yy_xx, g_yz_x_yy_xy, g_yz_x_yy_xz, g_yz_x_yy_yy, g_yz_x_yy_yz, g_yz_x_yy_zz, g_z_0_y_0_y_x_y_xx, g_z_0_y_0_y_x_y_xy, g_z_0_y_0_y_x_y_xz, g_z_0_y_0_y_x_y_yy, g_z_0_y_0_y_x_y_yz, g_z_0_y_0_y_x_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_x_y_xx[i] = -2.0 * g_yz_x_0_xx[i] * a_exp + 4.0 * g_yz_x_yy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_y_xy[i] = -2.0 * g_yz_x_0_xy[i] * a_exp + 4.0 * g_yz_x_yy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_y_xz[i] = -2.0 * g_yz_x_0_xz[i] * a_exp + 4.0 * g_yz_x_yy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_y_yy[i] = -2.0 * g_yz_x_0_yy[i] * a_exp + 4.0 * g_yz_x_yy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_y_yz[i] = -2.0 * g_yz_x_0_yz[i] * a_exp + 4.0 * g_yz_x_yy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_y_zz[i] = -2.0 * g_yz_x_0_zz[i] * a_exp + 4.0 * g_yz_x_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1200-1206)

    #pragma omp simd aligned(g_yz_x_yz_xx, g_yz_x_yz_xy, g_yz_x_yz_xz, g_yz_x_yz_yy, g_yz_x_yz_yz, g_yz_x_yz_zz, g_z_0_y_0_y_x_z_xx, g_z_0_y_0_y_x_z_xy, g_z_0_y_0_y_x_z_xz, g_z_0_y_0_y_x_z_yy, g_z_0_y_0_y_x_z_yz, g_z_0_y_0_y_x_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_x_z_xx[i] = 4.0 * g_yz_x_yz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_z_xy[i] = 4.0 * g_yz_x_yz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_z_xz[i] = 4.0 * g_yz_x_yz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_z_yy[i] = 4.0 * g_yz_x_yz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_z_yz[i] = 4.0 * g_yz_x_yz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_z_zz[i] = 4.0 * g_yz_x_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1206-1212)

    #pragma omp simd aligned(g_yz_y_xy_xx, g_yz_y_xy_xy, g_yz_y_xy_xz, g_yz_y_xy_yy, g_yz_y_xy_yz, g_yz_y_xy_zz, g_z_0_y_0_y_y_x_xx, g_z_0_y_0_y_y_x_xy, g_z_0_y_0_y_y_x_xz, g_z_0_y_0_y_y_x_yy, g_z_0_y_0_y_y_x_yz, g_z_0_y_0_y_y_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_y_x_xx[i] = 4.0 * g_yz_y_xy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_x_xy[i] = 4.0 * g_yz_y_xy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_x_xz[i] = 4.0 * g_yz_y_xy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_x_yy[i] = 4.0 * g_yz_y_xy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_x_yz[i] = 4.0 * g_yz_y_xy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_x_zz[i] = 4.0 * g_yz_y_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1212-1218)

    #pragma omp simd aligned(g_yz_y_0_xx, g_yz_y_0_xy, g_yz_y_0_xz, g_yz_y_0_yy, g_yz_y_0_yz, g_yz_y_0_zz, g_yz_y_yy_xx, g_yz_y_yy_xy, g_yz_y_yy_xz, g_yz_y_yy_yy, g_yz_y_yy_yz, g_yz_y_yy_zz, g_z_0_y_0_y_y_y_xx, g_z_0_y_0_y_y_y_xy, g_z_0_y_0_y_y_y_xz, g_z_0_y_0_y_y_y_yy, g_z_0_y_0_y_y_y_yz, g_z_0_y_0_y_y_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_y_y_xx[i] = -2.0 * g_yz_y_0_xx[i] * a_exp + 4.0 * g_yz_y_yy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_y_xy[i] = -2.0 * g_yz_y_0_xy[i] * a_exp + 4.0 * g_yz_y_yy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_y_xz[i] = -2.0 * g_yz_y_0_xz[i] * a_exp + 4.0 * g_yz_y_yy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_y_yy[i] = -2.0 * g_yz_y_0_yy[i] * a_exp + 4.0 * g_yz_y_yy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_y_yz[i] = -2.0 * g_yz_y_0_yz[i] * a_exp + 4.0 * g_yz_y_yy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_y_zz[i] = -2.0 * g_yz_y_0_zz[i] * a_exp + 4.0 * g_yz_y_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1218-1224)

    #pragma omp simd aligned(g_yz_y_yz_xx, g_yz_y_yz_xy, g_yz_y_yz_xz, g_yz_y_yz_yy, g_yz_y_yz_yz, g_yz_y_yz_zz, g_z_0_y_0_y_y_z_xx, g_z_0_y_0_y_y_z_xy, g_z_0_y_0_y_y_z_xz, g_z_0_y_0_y_y_z_yy, g_z_0_y_0_y_y_z_yz, g_z_0_y_0_y_y_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_y_z_xx[i] = 4.0 * g_yz_y_yz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_z_xy[i] = 4.0 * g_yz_y_yz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_z_xz[i] = 4.0 * g_yz_y_yz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_z_yy[i] = 4.0 * g_yz_y_yz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_z_yz[i] = 4.0 * g_yz_y_yz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_z_zz[i] = 4.0 * g_yz_y_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1224-1230)

    #pragma omp simd aligned(g_yz_z_xy_xx, g_yz_z_xy_xy, g_yz_z_xy_xz, g_yz_z_xy_yy, g_yz_z_xy_yz, g_yz_z_xy_zz, g_z_0_y_0_y_z_x_xx, g_z_0_y_0_y_z_x_xy, g_z_0_y_0_y_z_x_xz, g_z_0_y_0_y_z_x_yy, g_z_0_y_0_y_z_x_yz, g_z_0_y_0_y_z_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_z_x_xx[i] = 4.0 * g_yz_z_xy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_x_xy[i] = 4.0 * g_yz_z_xy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_x_xz[i] = 4.0 * g_yz_z_xy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_x_yy[i] = 4.0 * g_yz_z_xy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_x_yz[i] = 4.0 * g_yz_z_xy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_x_zz[i] = 4.0 * g_yz_z_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1230-1236)

    #pragma omp simd aligned(g_yz_z_0_xx, g_yz_z_0_xy, g_yz_z_0_xz, g_yz_z_0_yy, g_yz_z_0_yz, g_yz_z_0_zz, g_yz_z_yy_xx, g_yz_z_yy_xy, g_yz_z_yy_xz, g_yz_z_yy_yy, g_yz_z_yy_yz, g_yz_z_yy_zz, g_z_0_y_0_y_z_y_xx, g_z_0_y_0_y_z_y_xy, g_z_0_y_0_y_z_y_xz, g_z_0_y_0_y_z_y_yy, g_z_0_y_0_y_z_y_yz, g_z_0_y_0_y_z_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_z_y_xx[i] = -2.0 * g_yz_z_0_xx[i] * a_exp + 4.0 * g_yz_z_yy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_y_xy[i] = -2.0 * g_yz_z_0_xy[i] * a_exp + 4.0 * g_yz_z_yy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_y_xz[i] = -2.0 * g_yz_z_0_xz[i] * a_exp + 4.0 * g_yz_z_yy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_y_yy[i] = -2.0 * g_yz_z_0_yy[i] * a_exp + 4.0 * g_yz_z_yy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_y_yz[i] = -2.0 * g_yz_z_0_yz[i] * a_exp + 4.0 * g_yz_z_yy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_y_zz[i] = -2.0 * g_yz_z_0_zz[i] * a_exp + 4.0 * g_yz_z_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1236-1242)

    #pragma omp simd aligned(g_yz_z_yz_xx, g_yz_z_yz_xy, g_yz_z_yz_xz, g_yz_z_yz_yy, g_yz_z_yz_yz, g_yz_z_yz_zz, g_z_0_y_0_y_z_z_xx, g_z_0_y_0_y_z_z_xy, g_z_0_y_0_y_z_z_xz, g_z_0_y_0_y_z_z_yy, g_z_0_y_0_y_z_z_yz, g_z_0_y_0_y_z_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_z_z_xx[i] = 4.0 * g_yz_z_yz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_z_xy[i] = 4.0 * g_yz_z_yz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_z_xz[i] = 4.0 * g_yz_z_yz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_z_yy[i] = 4.0 * g_yz_z_yz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_z_yz[i] = 4.0 * g_yz_z_yz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_z_zz[i] = 4.0 * g_yz_z_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1242-1248)

    #pragma omp simd aligned(g_0_x_xy_xx, g_0_x_xy_xy, g_0_x_xy_xz, g_0_x_xy_yy, g_0_x_xy_yz, g_0_x_xy_zz, g_z_0_y_0_z_x_x_xx, g_z_0_y_0_z_x_x_xy, g_z_0_y_0_z_x_x_xz, g_z_0_y_0_z_x_x_yy, g_z_0_y_0_z_x_x_yz, g_z_0_y_0_z_x_x_zz, g_zz_x_xy_xx, g_zz_x_xy_xy, g_zz_x_xy_xz, g_zz_x_xy_yy, g_zz_x_xy_yz, g_zz_x_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_x_x_xx[i] = -2.0 * g_0_x_xy_xx[i] * c_exps[i] + 4.0 * g_zz_x_xy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_x_xy[i] = -2.0 * g_0_x_xy_xy[i] * c_exps[i] + 4.0 * g_zz_x_xy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_x_xz[i] = -2.0 * g_0_x_xy_xz[i] * c_exps[i] + 4.0 * g_zz_x_xy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_x_yy[i] = -2.0 * g_0_x_xy_yy[i] * c_exps[i] + 4.0 * g_zz_x_xy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_x_yz[i] = -2.0 * g_0_x_xy_yz[i] * c_exps[i] + 4.0 * g_zz_x_xy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_x_zz[i] = -2.0 * g_0_x_xy_zz[i] * c_exps[i] + 4.0 * g_zz_x_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1248-1254)

    #pragma omp simd aligned(g_0_x_0_xx, g_0_x_0_xy, g_0_x_0_xz, g_0_x_0_yy, g_0_x_0_yz, g_0_x_0_zz, g_0_x_yy_xx, g_0_x_yy_xy, g_0_x_yy_xz, g_0_x_yy_yy, g_0_x_yy_yz, g_0_x_yy_zz, g_z_0_y_0_z_x_y_xx, g_z_0_y_0_z_x_y_xy, g_z_0_y_0_z_x_y_xz, g_z_0_y_0_z_x_y_yy, g_z_0_y_0_z_x_y_yz, g_z_0_y_0_z_x_y_zz, g_zz_x_0_xx, g_zz_x_0_xy, g_zz_x_0_xz, g_zz_x_0_yy, g_zz_x_0_yz, g_zz_x_0_zz, g_zz_x_yy_xx, g_zz_x_yy_xy, g_zz_x_yy_xz, g_zz_x_yy_yy, g_zz_x_yy_yz, g_zz_x_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_x_y_xx[i] = g_0_x_0_xx[i] - 2.0 * g_0_x_yy_xx[i] * c_exps[i] - 2.0 * g_zz_x_0_xx[i] * a_exp + 4.0 * g_zz_x_yy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_y_xy[i] = g_0_x_0_xy[i] - 2.0 * g_0_x_yy_xy[i] * c_exps[i] - 2.0 * g_zz_x_0_xy[i] * a_exp + 4.0 * g_zz_x_yy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_y_xz[i] = g_0_x_0_xz[i] - 2.0 * g_0_x_yy_xz[i] * c_exps[i] - 2.0 * g_zz_x_0_xz[i] * a_exp + 4.0 * g_zz_x_yy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_y_yy[i] = g_0_x_0_yy[i] - 2.0 * g_0_x_yy_yy[i] * c_exps[i] - 2.0 * g_zz_x_0_yy[i] * a_exp + 4.0 * g_zz_x_yy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_y_yz[i] = g_0_x_0_yz[i] - 2.0 * g_0_x_yy_yz[i] * c_exps[i] - 2.0 * g_zz_x_0_yz[i] * a_exp + 4.0 * g_zz_x_yy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_y_zz[i] = g_0_x_0_zz[i] - 2.0 * g_0_x_yy_zz[i] * c_exps[i] - 2.0 * g_zz_x_0_zz[i] * a_exp + 4.0 * g_zz_x_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1254-1260)

    #pragma omp simd aligned(g_0_x_yz_xx, g_0_x_yz_xy, g_0_x_yz_xz, g_0_x_yz_yy, g_0_x_yz_yz, g_0_x_yz_zz, g_z_0_y_0_z_x_z_xx, g_z_0_y_0_z_x_z_xy, g_z_0_y_0_z_x_z_xz, g_z_0_y_0_z_x_z_yy, g_z_0_y_0_z_x_z_yz, g_z_0_y_0_z_x_z_zz, g_zz_x_yz_xx, g_zz_x_yz_xy, g_zz_x_yz_xz, g_zz_x_yz_yy, g_zz_x_yz_yz, g_zz_x_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_x_z_xx[i] = -2.0 * g_0_x_yz_xx[i] * c_exps[i] + 4.0 * g_zz_x_yz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_z_xy[i] = -2.0 * g_0_x_yz_xy[i] * c_exps[i] + 4.0 * g_zz_x_yz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_z_xz[i] = -2.0 * g_0_x_yz_xz[i] * c_exps[i] + 4.0 * g_zz_x_yz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_z_yy[i] = -2.0 * g_0_x_yz_yy[i] * c_exps[i] + 4.0 * g_zz_x_yz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_z_yz[i] = -2.0 * g_0_x_yz_yz[i] * c_exps[i] + 4.0 * g_zz_x_yz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_z_zz[i] = -2.0 * g_0_x_yz_zz[i] * c_exps[i] + 4.0 * g_zz_x_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1260-1266)

    #pragma omp simd aligned(g_0_y_xy_xx, g_0_y_xy_xy, g_0_y_xy_xz, g_0_y_xy_yy, g_0_y_xy_yz, g_0_y_xy_zz, g_z_0_y_0_z_y_x_xx, g_z_0_y_0_z_y_x_xy, g_z_0_y_0_z_y_x_xz, g_z_0_y_0_z_y_x_yy, g_z_0_y_0_z_y_x_yz, g_z_0_y_0_z_y_x_zz, g_zz_y_xy_xx, g_zz_y_xy_xy, g_zz_y_xy_xz, g_zz_y_xy_yy, g_zz_y_xy_yz, g_zz_y_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_y_x_xx[i] = -2.0 * g_0_y_xy_xx[i] * c_exps[i] + 4.0 * g_zz_y_xy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_x_xy[i] = -2.0 * g_0_y_xy_xy[i] * c_exps[i] + 4.0 * g_zz_y_xy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_x_xz[i] = -2.0 * g_0_y_xy_xz[i] * c_exps[i] + 4.0 * g_zz_y_xy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_x_yy[i] = -2.0 * g_0_y_xy_yy[i] * c_exps[i] + 4.0 * g_zz_y_xy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_x_yz[i] = -2.0 * g_0_y_xy_yz[i] * c_exps[i] + 4.0 * g_zz_y_xy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_x_zz[i] = -2.0 * g_0_y_xy_zz[i] * c_exps[i] + 4.0 * g_zz_y_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1266-1272)

    #pragma omp simd aligned(g_0_y_0_xx, g_0_y_0_xy, g_0_y_0_xz, g_0_y_0_yy, g_0_y_0_yz, g_0_y_0_zz, g_0_y_yy_xx, g_0_y_yy_xy, g_0_y_yy_xz, g_0_y_yy_yy, g_0_y_yy_yz, g_0_y_yy_zz, g_z_0_y_0_z_y_y_xx, g_z_0_y_0_z_y_y_xy, g_z_0_y_0_z_y_y_xz, g_z_0_y_0_z_y_y_yy, g_z_0_y_0_z_y_y_yz, g_z_0_y_0_z_y_y_zz, g_zz_y_0_xx, g_zz_y_0_xy, g_zz_y_0_xz, g_zz_y_0_yy, g_zz_y_0_yz, g_zz_y_0_zz, g_zz_y_yy_xx, g_zz_y_yy_xy, g_zz_y_yy_xz, g_zz_y_yy_yy, g_zz_y_yy_yz, g_zz_y_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_y_y_xx[i] = g_0_y_0_xx[i] - 2.0 * g_0_y_yy_xx[i] * c_exps[i] - 2.0 * g_zz_y_0_xx[i] * a_exp + 4.0 * g_zz_y_yy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_y_xy[i] = g_0_y_0_xy[i] - 2.0 * g_0_y_yy_xy[i] * c_exps[i] - 2.0 * g_zz_y_0_xy[i] * a_exp + 4.0 * g_zz_y_yy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_y_xz[i] = g_0_y_0_xz[i] - 2.0 * g_0_y_yy_xz[i] * c_exps[i] - 2.0 * g_zz_y_0_xz[i] * a_exp + 4.0 * g_zz_y_yy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_y_yy[i] = g_0_y_0_yy[i] - 2.0 * g_0_y_yy_yy[i] * c_exps[i] - 2.0 * g_zz_y_0_yy[i] * a_exp + 4.0 * g_zz_y_yy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_y_yz[i] = g_0_y_0_yz[i] - 2.0 * g_0_y_yy_yz[i] * c_exps[i] - 2.0 * g_zz_y_0_yz[i] * a_exp + 4.0 * g_zz_y_yy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_y_zz[i] = g_0_y_0_zz[i] - 2.0 * g_0_y_yy_zz[i] * c_exps[i] - 2.0 * g_zz_y_0_zz[i] * a_exp + 4.0 * g_zz_y_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1272-1278)

    #pragma omp simd aligned(g_0_y_yz_xx, g_0_y_yz_xy, g_0_y_yz_xz, g_0_y_yz_yy, g_0_y_yz_yz, g_0_y_yz_zz, g_z_0_y_0_z_y_z_xx, g_z_0_y_0_z_y_z_xy, g_z_0_y_0_z_y_z_xz, g_z_0_y_0_z_y_z_yy, g_z_0_y_0_z_y_z_yz, g_z_0_y_0_z_y_z_zz, g_zz_y_yz_xx, g_zz_y_yz_xy, g_zz_y_yz_xz, g_zz_y_yz_yy, g_zz_y_yz_yz, g_zz_y_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_y_z_xx[i] = -2.0 * g_0_y_yz_xx[i] * c_exps[i] + 4.0 * g_zz_y_yz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_z_xy[i] = -2.0 * g_0_y_yz_xy[i] * c_exps[i] + 4.0 * g_zz_y_yz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_z_xz[i] = -2.0 * g_0_y_yz_xz[i] * c_exps[i] + 4.0 * g_zz_y_yz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_z_yy[i] = -2.0 * g_0_y_yz_yy[i] * c_exps[i] + 4.0 * g_zz_y_yz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_z_yz[i] = -2.0 * g_0_y_yz_yz[i] * c_exps[i] + 4.0 * g_zz_y_yz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_z_zz[i] = -2.0 * g_0_y_yz_zz[i] * c_exps[i] + 4.0 * g_zz_y_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1278-1284)

    #pragma omp simd aligned(g_0_z_xy_xx, g_0_z_xy_xy, g_0_z_xy_xz, g_0_z_xy_yy, g_0_z_xy_yz, g_0_z_xy_zz, g_z_0_y_0_z_z_x_xx, g_z_0_y_0_z_z_x_xy, g_z_0_y_0_z_z_x_xz, g_z_0_y_0_z_z_x_yy, g_z_0_y_0_z_z_x_yz, g_z_0_y_0_z_z_x_zz, g_zz_z_xy_xx, g_zz_z_xy_xy, g_zz_z_xy_xz, g_zz_z_xy_yy, g_zz_z_xy_yz, g_zz_z_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_z_x_xx[i] = -2.0 * g_0_z_xy_xx[i] * c_exps[i] + 4.0 * g_zz_z_xy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_x_xy[i] = -2.0 * g_0_z_xy_xy[i] * c_exps[i] + 4.0 * g_zz_z_xy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_x_xz[i] = -2.0 * g_0_z_xy_xz[i] * c_exps[i] + 4.0 * g_zz_z_xy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_x_yy[i] = -2.0 * g_0_z_xy_yy[i] * c_exps[i] + 4.0 * g_zz_z_xy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_x_yz[i] = -2.0 * g_0_z_xy_yz[i] * c_exps[i] + 4.0 * g_zz_z_xy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_x_zz[i] = -2.0 * g_0_z_xy_zz[i] * c_exps[i] + 4.0 * g_zz_z_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1284-1290)

    #pragma omp simd aligned(g_0_z_0_xx, g_0_z_0_xy, g_0_z_0_xz, g_0_z_0_yy, g_0_z_0_yz, g_0_z_0_zz, g_0_z_yy_xx, g_0_z_yy_xy, g_0_z_yy_xz, g_0_z_yy_yy, g_0_z_yy_yz, g_0_z_yy_zz, g_z_0_y_0_z_z_y_xx, g_z_0_y_0_z_z_y_xy, g_z_0_y_0_z_z_y_xz, g_z_0_y_0_z_z_y_yy, g_z_0_y_0_z_z_y_yz, g_z_0_y_0_z_z_y_zz, g_zz_z_0_xx, g_zz_z_0_xy, g_zz_z_0_xz, g_zz_z_0_yy, g_zz_z_0_yz, g_zz_z_0_zz, g_zz_z_yy_xx, g_zz_z_yy_xy, g_zz_z_yy_xz, g_zz_z_yy_yy, g_zz_z_yy_yz, g_zz_z_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_z_y_xx[i] = g_0_z_0_xx[i] - 2.0 * g_0_z_yy_xx[i] * c_exps[i] - 2.0 * g_zz_z_0_xx[i] * a_exp + 4.0 * g_zz_z_yy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_y_xy[i] = g_0_z_0_xy[i] - 2.0 * g_0_z_yy_xy[i] * c_exps[i] - 2.0 * g_zz_z_0_xy[i] * a_exp + 4.0 * g_zz_z_yy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_y_xz[i] = g_0_z_0_xz[i] - 2.0 * g_0_z_yy_xz[i] * c_exps[i] - 2.0 * g_zz_z_0_xz[i] * a_exp + 4.0 * g_zz_z_yy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_y_yy[i] = g_0_z_0_yy[i] - 2.0 * g_0_z_yy_yy[i] * c_exps[i] - 2.0 * g_zz_z_0_yy[i] * a_exp + 4.0 * g_zz_z_yy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_y_yz[i] = g_0_z_0_yz[i] - 2.0 * g_0_z_yy_yz[i] * c_exps[i] - 2.0 * g_zz_z_0_yz[i] * a_exp + 4.0 * g_zz_z_yy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_y_zz[i] = g_0_z_0_zz[i] - 2.0 * g_0_z_yy_zz[i] * c_exps[i] - 2.0 * g_zz_z_0_zz[i] * a_exp + 4.0 * g_zz_z_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1290-1296)

    #pragma omp simd aligned(g_0_z_yz_xx, g_0_z_yz_xy, g_0_z_yz_xz, g_0_z_yz_yy, g_0_z_yz_yz, g_0_z_yz_zz, g_z_0_y_0_z_z_z_xx, g_z_0_y_0_z_z_z_xy, g_z_0_y_0_z_z_z_xz, g_z_0_y_0_z_z_z_yy, g_z_0_y_0_z_z_z_yz, g_z_0_y_0_z_z_z_zz, g_zz_z_yz_xx, g_zz_z_yz_xy, g_zz_z_yz_xz, g_zz_z_yz_yy, g_zz_z_yz_yz, g_zz_z_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_z_z_xx[i] = -2.0 * g_0_z_yz_xx[i] * c_exps[i] + 4.0 * g_zz_z_yz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_z_xy[i] = -2.0 * g_0_z_yz_xy[i] * c_exps[i] + 4.0 * g_zz_z_yz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_z_xz[i] = -2.0 * g_0_z_yz_xz[i] * c_exps[i] + 4.0 * g_zz_z_yz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_z_yy[i] = -2.0 * g_0_z_yz_yy[i] * c_exps[i] + 4.0 * g_zz_z_yz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_z_yz[i] = -2.0 * g_0_z_yz_yz[i] * c_exps[i] + 4.0 * g_zz_z_yz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_z_zz[i] = -2.0 * g_0_z_yz_zz[i] * c_exps[i] + 4.0 * g_zz_z_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1296-1302)

    #pragma omp simd aligned(g_xz_x_xz_xx, g_xz_x_xz_xy, g_xz_x_xz_xz, g_xz_x_xz_yy, g_xz_x_xz_yz, g_xz_x_xz_zz, g_z_0_z_0_x_x_x_xx, g_z_0_z_0_x_x_x_xy, g_z_0_z_0_x_x_x_xz, g_z_0_z_0_x_x_x_yy, g_z_0_z_0_x_x_x_yz, g_z_0_z_0_x_x_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_x_x_xx[i] = 4.0 * g_xz_x_xz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_x_xy[i] = 4.0 * g_xz_x_xz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_x_xz[i] = 4.0 * g_xz_x_xz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_x_yy[i] = 4.0 * g_xz_x_xz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_x_yz[i] = 4.0 * g_xz_x_xz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_x_zz[i] = 4.0 * g_xz_x_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1302-1308)

    #pragma omp simd aligned(g_xz_x_yz_xx, g_xz_x_yz_xy, g_xz_x_yz_xz, g_xz_x_yz_yy, g_xz_x_yz_yz, g_xz_x_yz_zz, g_z_0_z_0_x_x_y_xx, g_z_0_z_0_x_x_y_xy, g_z_0_z_0_x_x_y_xz, g_z_0_z_0_x_x_y_yy, g_z_0_z_0_x_x_y_yz, g_z_0_z_0_x_x_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_x_y_xx[i] = 4.0 * g_xz_x_yz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_y_xy[i] = 4.0 * g_xz_x_yz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_y_xz[i] = 4.0 * g_xz_x_yz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_y_yy[i] = 4.0 * g_xz_x_yz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_y_yz[i] = 4.0 * g_xz_x_yz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_y_zz[i] = 4.0 * g_xz_x_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1308-1314)

    #pragma omp simd aligned(g_xz_x_0_xx, g_xz_x_0_xy, g_xz_x_0_xz, g_xz_x_0_yy, g_xz_x_0_yz, g_xz_x_0_zz, g_xz_x_zz_xx, g_xz_x_zz_xy, g_xz_x_zz_xz, g_xz_x_zz_yy, g_xz_x_zz_yz, g_xz_x_zz_zz, g_z_0_z_0_x_x_z_xx, g_z_0_z_0_x_x_z_xy, g_z_0_z_0_x_x_z_xz, g_z_0_z_0_x_x_z_yy, g_z_0_z_0_x_x_z_yz, g_z_0_z_0_x_x_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_x_z_xx[i] = -2.0 * g_xz_x_0_xx[i] * a_exp + 4.0 * g_xz_x_zz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_z_xy[i] = -2.0 * g_xz_x_0_xy[i] * a_exp + 4.0 * g_xz_x_zz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_z_xz[i] = -2.0 * g_xz_x_0_xz[i] * a_exp + 4.0 * g_xz_x_zz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_z_yy[i] = -2.0 * g_xz_x_0_yy[i] * a_exp + 4.0 * g_xz_x_zz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_z_yz[i] = -2.0 * g_xz_x_0_yz[i] * a_exp + 4.0 * g_xz_x_zz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_z_zz[i] = -2.0 * g_xz_x_0_zz[i] * a_exp + 4.0 * g_xz_x_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1314-1320)

    #pragma omp simd aligned(g_xz_y_xz_xx, g_xz_y_xz_xy, g_xz_y_xz_xz, g_xz_y_xz_yy, g_xz_y_xz_yz, g_xz_y_xz_zz, g_z_0_z_0_x_y_x_xx, g_z_0_z_0_x_y_x_xy, g_z_0_z_0_x_y_x_xz, g_z_0_z_0_x_y_x_yy, g_z_0_z_0_x_y_x_yz, g_z_0_z_0_x_y_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_y_x_xx[i] = 4.0 * g_xz_y_xz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_x_xy[i] = 4.0 * g_xz_y_xz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_x_xz[i] = 4.0 * g_xz_y_xz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_x_yy[i] = 4.0 * g_xz_y_xz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_x_yz[i] = 4.0 * g_xz_y_xz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_x_zz[i] = 4.0 * g_xz_y_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1320-1326)

    #pragma omp simd aligned(g_xz_y_yz_xx, g_xz_y_yz_xy, g_xz_y_yz_xz, g_xz_y_yz_yy, g_xz_y_yz_yz, g_xz_y_yz_zz, g_z_0_z_0_x_y_y_xx, g_z_0_z_0_x_y_y_xy, g_z_0_z_0_x_y_y_xz, g_z_0_z_0_x_y_y_yy, g_z_0_z_0_x_y_y_yz, g_z_0_z_0_x_y_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_y_y_xx[i] = 4.0 * g_xz_y_yz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_y_xy[i] = 4.0 * g_xz_y_yz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_y_xz[i] = 4.0 * g_xz_y_yz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_y_yy[i] = 4.0 * g_xz_y_yz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_y_yz[i] = 4.0 * g_xz_y_yz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_y_zz[i] = 4.0 * g_xz_y_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1326-1332)

    #pragma omp simd aligned(g_xz_y_0_xx, g_xz_y_0_xy, g_xz_y_0_xz, g_xz_y_0_yy, g_xz_y_0_yz, g_xz_y_0_zz, g_xz_y_zz_xx, g_xz_y_zz_xy, g_xz_y_zz_xz, g_xz_y_zz_yy, g_xz_y_zz_yz, g_xz_y_zz_zz, g_z_0_z_0_x_y_z_xx, g_z_0_z_0_x_y_z_xy, g_z_0_z_0_x_y_z_xz, g_z_0_z_0_x_y_z_yy, g_z_0_z_0_x_y_z_yz, g_z_0_z_0_x_y_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_y_z_xx[i] = -2.0 * g_xz_y_0_xx[i] * a_exp + 4.0 * g_xz_y_zz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_z_xy[i] = -2.0 * g_xz_y_0_xy[i] * a_exp + 4.0 * g_xz_y_zz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_z_xz[i] = -2.0 * g_xz_y_0_xz[i] * a_exp + 4.0 * g_xz_y_zz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_z_yy[i] = -2.0 * g_xz_y_0_yy[i] * a_exp + 4.0 * g_xz_y_zz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_z_yz[i] = -2.0 * g_xz_y_0_yz[i] * a_exp + 4.0 * g_xz_y_zz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_z_zz[i] = -2.0 * g_xz_y_0_zz[i] * a_exp + 4.0 * g_xz_y_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1332-1338)

    #pragma omp simd aligned(g_xz_z_xz_xx, g_xz_z_xz_xy, g_xz_z_xz_xz, g_xz_z_xz_yy, g_xz_z_xz_yz, g_xz_z_xz_zz, g_z_0_z_0_x_z_x_xx, g_z_0_z_0_x_z_x_xy, g_z_0_z_0_x_z_x_xz, g_z_0_z_0_x_z_x_yy, g_z_0_z_0_x_z_x_yz, g_z_0_z_0_x_z_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_z_x_xx[i] = 4.0 * g_xz_z_xz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_x_xy[i] = 4.0 * g_xz_z_xz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_x_xz[i] = 4.0 * g_xz_z_xz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_x_yy[i] = 4.0 * g_xz_z_xz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_x_yz[i] = 4.0 * g_xz_z_xz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_x_zz[i] = 4.0 * g_xz_z_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1338-1344)

    #pragma omp simd aligned(g_xz_z_yz_xx, g_xz_z_yz_xy, g_xz_z_yz_xz, g_xz_z_yz_yy, g_xz_z_yz_yz, g_xz_z_yz_zz, g_z_0_z_0_x_z_y_xx, g_z_0_z_0_x_z_y_xy, g_z_0_z_0_x_z_y_xz, g_z_0_z_0_x_z_y_yy, g_z_0_z_0_x_z_y_yz, g_z_0_z_0_x_z_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_z_y_xx[i] = 4.0 * g_xz_z_yz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_y_xy[i] = 4.0 * g_xz_z_yz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_y_xz[i] = 4.0 * g_xz_z_yz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_y_yy[i] = 4.0 * g_xz_z_yz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_y_yz[i] = 4.0 * g_xz_z_yz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_y_zz[i] = 4.0 * g_xz_z_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1344-1350)

    #pragma omp simd aligned(g_xz_z_0_xx, g_xz_z_0_xy, g_xz_z_0_xz, g_xz_z_0_yy, g_xz_z_0_yz, g_xz_z_0_zz, g_xz_z_zz_xx, g_xz_z_zz_xy, g_xz_z_zz_xz, g_xz_z_zz_yy, g_xz_z_zz_yz, g_xz_z_zz_zz, g_z_0_z_0_x_z_z_xx, g_z_0_z_0_x_z_z_xy, g_z_0_z_0_x_z_z_xz, g_z_0_z_0_x_z_z_yy, g_z_0_z_0_x_z_z_yz, g_z_0_z_0_x_z_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_z_z_xx[i] = -2.0 * g_xz_z_0_xx[i] * a_exp + 4.0 * g_xz_z_zz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_z_xy[i] = -2.0 * g_xz_z_0_xy[i] * a_exp + 4.0 * g_xz_z_zz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_z_xz[i] = -2.0 * g_xz_z_0_xz[i] * a_exp + 4.0 * g_xz_z_zz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_z_yy[i] = -2.0 * g_xz_z_0_yy[i] * a_exp + 4.0 * g_xz_z_zz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_z_yz[i] = -2.0 * g_xz_z_0_yz[i] * a_exp + 4.0 * g_xz_z_zz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_z_zz[i] = -2.0 * g_xz_z_0_zz[i] * a_exp + 4.0 * g_xz_z_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1350-1356)

    #pragma omp simd aligned(g_yz_x_xz_xx, g_yz_x_xz_xy, g_yz_x_xz_xz, g_yz_x_xz_yy, g_yz_x_xz_yz, g_yz_x_xz_zz, g_z_0_z_0_y_x_x_xx, g_z_0_z_0_y_x_x_xy, g_z_0_z_0_y_x_x_xz, g_z_0_z_0_y_x_x_yy, g_z_0_z_0_y_x_x_yz, g_z_0_z_0_y_x_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_x_x_xx[i] = 4.0 * g_yz_x_xz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_x_xy[i] = 4.0 * g_yz_x_xz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_x_xz[i] = 4.0 * g_yz_x_xz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_x_yy[i] = 4.0 * g_yz_x_xz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_x_yz[i] = 4.0 * g_yz_x_xz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_x_zz[i] = 4.0 * g_yz_x_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1356-1362)

    #pragma omp simd aligned(g_yz_x_yz_xx, g_yz_x_yz_xy, g_yz_x_yz_xz, g_yz_x_yz_yy, g_yz_x_yz_yz, g_yz_x_yz_zz, g_z_0_z_0_y_x_y_xx, g_z_0_z_0_y_x_y_xy, g_z_0_z_0_y_x_y_xz, g_z_0_z_0_y_x_y_yy, g_z_0_z_0_y_x_y_yz, g_z_0_z_0_y_x_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_x_y_xx[i] = 4.0 * g_yz_x_yz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_y_xy[i] = 4.0 * g_yz_x_yz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_y_xz[i] = 4.0 * g_yz_x_yz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_y_yy[i] = 4.0 * g_yz_x_yz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_y_yz[i] = 4.0 * g_yz_x_yz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_y_zz[i] = 4.0 * g_yz_x_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1362-1368)

    #pragma omp simd aligned(g_yz_x_0_xx, g_yz_x_0_xy, g_yz_x_0_xz, g_yz_x_0_yy, g_yz_x_0_yz, g_yz_x_0_zz, g_yz_x_zz_xx, g_yz_x_zz_xy, g_yz_x_zz_xz, g_yz_x_zz_yy, g_yz_x_zz_yz, g_yz_x_zz_zz, g_z_0_z_0_y_x_z_xx, g_z_0_z_0_y_x_z_xy, g_z_0_z_0_y_x_z_xz, g_z_0_z_0_y_x_z_yy, g_z_0_z_0_y_x_z_yz, g_z_0_z_0_y_x_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_x_z_xx[i] = -2.0 * g_yz_x_0_xx[i] * a_exp + 4.0 * g_yz_x_zz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_z_xy[i] = -2.0 * g_yz_x_0_xy[i] * a_exp + 4.0 * g_yz_x_zz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_z_xz[i] = -2.0 * g_yz_x_0_xz[i] * a_exp + 4.0 * g_yz_x_zz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_z_yy[i] = -2.0 * g_yz_x_0_yy[i] * a_exp + 4.0 * g_yz_x_zz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_z_yz[i] = -2.0 * g_yz_x_0_yz[i] * a_exp + 4.0 * g_yz_x_zz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_z_zz[i] = -2.0 * g_yz_x_0_zz[i] * a_exp + 4.0 * g_yz_x_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1368-1374)

    #pragma omp simd aligned(g_yz_y_xz_xx, g_yz_y_xz_xy, g_yz_y_xz_xz, g_yz_y_xz_yy, g_yz_y_xz_yz, g_yz_y_xz_zz, g_z_0_z_0_y_y_x_xx, g_z_0_z_0_y_y_x_xy, g_z_0_z_0_y_y_x_xz, g_z_0_z_0_y_y_x_yy, g_z_0_z_0_y_y_x_yz, g_z_0_z_0_y_y_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_y_x_xx[i] = 4.0 * g_yz_y_xz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_x_xy[i] = 4.0 * g_yz_y_xz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_x_xz[i] = 4.0 * g_yz_y_xz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_x_yy[i] = 4.0 * g_yz_y_xz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_x_yz[i] = 4.0 * g_yz_y_xz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_x_zz[i] = 4.0 * g_yz_y_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1374-1380)

    #pragma omp simd aligned(g_yz_y_yz_xx, g_yz_y_yz_xy, g_yz_y_yz_xz, g_yz_y_yz_yy, g_yz_y_yz_yz, g_yz_y_yz_zz, g_z_0_z_0_y_y_y_xx, g_z_0_z_0_y_y_y_xy, g_z_0_z_0_y_y_y_xz, g_z_0_z_0_y_y_y_yy, g_z_0_z_0_y_y_y_yz, g_z_0_z_0_y_y_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_y_y_xx[i] = 4.0 * g_yz_y_yz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_y_xy[i] = 4.0 * g_yz_y_yz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_y_xz[i] = 4.0 * g_yz_y_yz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_y_yy[i] = 4.0 * g_yz_y_yz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_y_yz[i] = 4.0 * g_yz_y_yz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_y_zz[i] = 4.0 * g_yz_y_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1380-1386)

    #pragma omp simd aligned(g_yz_y_0_xx, g_yz_y_0_xy, g_yz_y_0_xz, g_yz_y_0_yy, g_yz_y_0_yz, g_yz_y_0_zz, g_yz_y_zz_xx, g_yz_y_zz_xy, g_yz_y_zz_xz, g_yz_y_zz_yy, g_yz_y_zz_yz, g_yz_y_zz_zz, g_z_0_z_0_y_y_z_xx, g_z_0_z_0_y_y_z_xy, g_z_0_z_0_y_y_z_xz, g_z_0_z_0_y_y_z_yy, g_z_0_z_0_y_y_z_yz, g_z_0_z_0_y_y_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_y_z_xx[i] = -2.0 * g_yz_y_0_xx[i] * a_exp + 4.0 * g_yz_y_zz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_z_xy[i] = -2.0 * g_yz_y_0_xy[i] * a_exp + 4.0 * g_yz_y_zz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_z_xz[i] = -2.0 * g_yz_y_0_xz[i] * a_exp + 4.0 * g_yz_y_zz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_z_yy[i] = -2.0 * g_yz_y_0_yy[i] * a_exp + 4.0 * g_yz_y_zz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_z_yz[i] = -2.0 * g_yz_y_0_yz[i] * a_exp + 4.0 * g_yz_y_zz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_z_zz[i] = -2.0 * g_yz_y_0_zz[i] * a_exp + 4.0 * g_yz_y_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1386-1392)

    #pragma omp simd aligned(g_yz_z_xz_xx, g_yz_z_xz_xy, g_yz_z_xz_xz, g_yz_z_xz_yy, g_yz_z_xz_yz, g_yz_z_xz_zz, g_z_0_z_0_y_z_x_xx, g_z_0_z_0_y_z_x_xy, g_z_0_z_0_y_z_x_xz, g_z_0_z_0_y_z_x_yy, g_z_0_z_0_y_z_x_yz, g_z_0_z_0_y_z_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_z_x_xx[i] = 4.0 * g_yz_z_xz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_x_xy[i] = 4.0 * g_yz_z_xz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_x_xz[i] = 4.0 * g_yz_z_xz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_x_yy[i] = 4.0 * g_yz_z_xz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_x_yz[i] = 4.0 * g_yz_z_xz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_x_zz[i] = 4.0 * g_yz_z_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1392-1398)

    #pragma omp simd aligned(g_yz_z_yz_xx, g_yz_z_yz_xy, g_yz_z_yz_xz, g_yz_z_yz_yy, g_yz_z_yz_yz, g_yz_z_yz_zz, g_z_0_z_0_y_z_y_xx, g_z_0_z_0_y_z_y_xy, g_z_0_z_0_y_z_y_xz, g_z_0_z_0_y_z_y_yy, g_z_0_z_0_y_z_y_yz, g_z_0_z_0_y_z_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_z_y_xx[i] = 4.0 * g_yz_z_yz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_y_xy[i] = 4.0 * g_yz_z_yz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_y_xz[i] = 4.0 * g_yz_z_yz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_y_yy[i] = 4.0 * g_yz_z_yz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_y_yz[i] = 4.0 * g_yz_z_yz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_y_zz[i] = 4.0 * g_yz_z_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1398-1404)

    #pragma omp simd aligned(g_yz_z_0_xx, g_yz_z_0_xy, g_yz_z_0_xz, g_yz_z_0_yy, g_yz_z_0_yz, g_yz_z_0_zz, g_yz_z_zz_xx, g_yz_z_zz_xy, g_yz_z_zz_xz, g_yz_z_zz_yy, g_yz_z_zz_yz, g_yz_z_zz_zz, g_z_0_z_0_y_z_z_xx, g_z_0_z_0_y_z_z_xy, g_z_0_z_0_y_z_z_xz, g_z_0_z_0_y_z_z_yy, g_z_0_z_0_y_z_z_yz, g_z_0_z_0_y_z_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_z_z_xx[i] = -2.0 * g_yz_z_0_xx[i] * a_exp + 4.0 * g_yz_z_zz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_z_xy[i] = -2.0 * g_yz_z_0_xy[i] * a_exp + 4.0 * g_yz_z_zz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_z_xz[i] = -2.0 * g_yz_z_0_xz[i] * a_exp + 4.0 * g_yz_z_zz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_z_yy[i] = -2.0 * g_yz_z_0_yy[i] * a_exp + 4.0 * g_yz_z_zz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_z_yz[i] = -2.0 * g_yz_z_0_yz[i] * a_exp + 4.0 * g_yz_z_zz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_z_zz[i] = -2.0 * g_yz_z_0_zz[i] * a_exp + 4.0 * g_yz_z_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1404-1410)

    #pragma omp simd aligned(g_0_x_xz_xx, g_0_x_xz_xy, g_0_x_xz_xz, g_0_x_xz_yy, g_0_x_xz_yz, g_0_x_xz_zz, g_z_0_z_0_z_x_x_xx, g_z_0_z_0_z_x_x_xy, g_z_0_z_0_z_x_x_xz, g_z_0_z_0_z_x_x_yy, g_z_0_z_0_z_x_x_yz, g_z_0_z_0_z_x_x_zz, g_zz_x_xz_xx, g_zz_x_xz_xy, g_zz_x_xz_xz, g_zz_x_xz_yy, g_zz_x_xz_yz, g_zz_x_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_x_x_xx[i] = -2.0 * g_0_x_xz_xx[i] * c_exps[i] + 4.0 * g_zz_x_xz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_x_xy[i] = -2.0 * g_0_x_xz_xy[i] * c_exps[i] + 4.0 * g_zz_x_xz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_x_xz[i] = -2.0 * g_0_x_xz_xz[i] * c_exps[i] + 4.0 * g_zz_x_xz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_x_yy[i] = -2.0 * g_0_x_xz_yy[i] * c_exps[i] + 4.0 * g_zz_x_xz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_x_yz[i] = -2.0 * g_0_x_xz_yz[i] * c_exps[i] + 4.0 * g_zz_x_xz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_x_zz[i] = -2.0 * g_0_x_xz_zz[i] * c_exps[i] + 4.0 * g_zz_x_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1410-1416)

    #pragma omp simd aligned(g_0_x_yz_xx, g_0_x_yz_xy, g_0_x_yz_xz, g_0_x_yz_yy, g_0_x_yz_yz, g_0_x_yz_zz, g_z_0_z_0_z_x_y_xx, g_z_0_z_0_z_x_y_xy, g_z_0_z_0_z_x_y_xz, g_z_0_z_0_z_x_y_yy, g_z_0_z_0_z_x_y_yz, g_z_0_z_0_z_x_y_zz, g_zz_x_yz_xx, g_zz_x_yz_xy, g_zz_x_yz_xz, g_zz_x_yz_yy, g_zz_x_yz_yz, g_zz_x_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_x_y_xx[i] = -2.0 * g_0_x_yz_xx[i] * c_exps[i] + 4.0 * g_zz_x_yz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_y_xy[i] = -2.0 * g_0_x_yz_xy[i] * c_exps[i] + 4.0 * g_zz_x_yz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_y_xz[i] = -2.0 * g_0_x_yz_xz[i] * c_exps[i] + 4.0 * g_zz_x_yz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_y_yy[i] = -2.0 * g_0_x_yz_yy[i] * c_exps[i] + 4.0 * g_zz_x_yz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_y_yz[i] = -2.0 * g_0_x_yz_yz[i] * c_exps[i] + 4.0 * g_zz_x_yz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_y_zz[i] = -2.0 * g_0_x_yz_zz[i] * c_exps[i] + 4.0 * g_zz_x_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1416-1422)

    #pragma omp simd aligned(g_0_x_0_xx, g_0_x_0_xy, g_0_x_0_xz, g_0_x_0_yy, g_0_x_0_yz, g_0_x_0_zz, g_0_x_zz_xx, g_0_x_zz_xy, g_0_x_zz_xz, g_0_x_zz_yy, g_0_x_zz_yz, g_0_x_zz_zz, g_z_0_z_0_z_x_z_xx, g_z_0_z_0_z_x_z_xy, g_z_0_z_0_z_x_z_xz, g_z_0_z_0_z_x_z_yy, g_z_0_z_0_z_x_z_yz, g_z_0_z_0_z_x_z_zz, g_zz_x_0_xx, g_zz_x_0_xy, g_zz_x_0_xz, g_zz_x_0_yy, g_zz_x_0_yz, g_zz_x_0_zz, g_zz_x_zz_xx, g_zz_x_zz_xy, g_zz_x_zz_xz, g_zz_x_zz_yy, g_zz_x_zz_yz, g_zz_x_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_x_z_xx[i] = g_0_x_0_xx[i] - 2.0 * g_0_x_zz_xx[i] * c_exps[i] - 2.0 * g_zz_x_0_xx[i] * a_exp + 4.0 * g_zz_x_zz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_z_xy[i] = g_0_x_0_xy[i] - 2.0 * g_0_x_zz_xy[i] * c_exps[i] - 2.0 * g_zz_x_0_xy[i] * a_exp + 4.0 * g_zz_x_zz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_z_xz[i] = g_0_x_0_xz[i] - 2.0 * g_0_x_zz_xz[i] * c_exps[i] - 2.0 * g_zz_x_0_xz[i] * a_exp + 4.0 * g_zz_x_zz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_z_yy[i] = g_0_x_0_yy[i] - 2.0 * g_0_x_zz_yy[i] * c_exps[i] - 2.0 * g_zz_x_0_yy[i] * a_exp + 4.0 * g_zz_x_zz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_z_yz[i] = g_0_x_0_yz[i] - 2.0 * g_0_x_zz_yz[i] * c_exps[i] - 2.0 * g_zz_x_0_yz[i] * a_exp + 4.0 * g_zz_x_zz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_z_zz[i] = g_0_x_0_zz[i] - 2.0 * g_0_x_zz_zz[i] * c_exps[i] - 2.0 * g_zz_x_0_zz[i] * a_exp + 4.0 * g_zz_x_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1422-1428)

    #pragma omp simd aligned(g_0_y_xz_xx, g_0_y_xz_xy, g_0_y_xz_xz, g_0_y_xz_yy, g_0_y_xz_yz, g_0_y_xz_zz, g_z_0_z_0_z_y_x_xx, g_z_0_z_0_z_y_x_xy, g_z_0_z_0_z_y_x_xz, g_z_0_z_0_z_y_x_yy, g_z_0_z_0_z_y_x_yz, g_z_0_z_0_z_y_x_zz, g_zz_y_xz_xx, g_zz_y_xz_xy, g_zz_y_xz_xz, g_zz_y_xz_yy, g_zz_y_xz_yz, g_zz_y_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_y_x_xx[i] = -2.0 * g_0_y_xz_xx[i] * c_exps[i] + 4.0 * g_zz_y_xz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_x_xy[i] = -2.0 * g_0_y_xz_xy[i] * c_exps[i] + 4.0 * g_zz_y_xz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_x_xz[i] = -2.0 * g_0_y_xz_xz[i] * c_exps[i] + 4.0 * g_zz_y_xz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_x_yy[i] = -2.0 * g_0_y_xz_yy[i] * c_exps[i] + 4.0 * g_zz_y_xz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_x_yz[i] = -2.0 * g_0_y_xz_yz[i] * c_exps[i] + 4.0 * g_zz_y_xz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_x_zz[i] = -2.0 * g_0_y_xz_zz[i] * c_exps[i] + 4.0 * g_zz_y_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1428-1434)

    #pragma omp simd aligned(g_0_y_yz_xx, g_0_y_yz_xy, g_0_y_yz_xz, g_0_y_yz_yy, g_0_y_yz_yz, g_0_y_yz_zz, g_z_0_z_0_z_y_y_xx, g_z_0_z_0_z_y_y_xy, g_z_0_z_0_z_y_y_xz, g_z_0_z_0_z_y_y_yy, g_z_0_z_0_z_y_y_yz, g_z_0_z_0_z_y_y_zz, g_zz_y_yz_xx, g_zz_y_yz_xy, g_zz_y_yz_xz, g_zz_y_yz_yy, g_zz_y_yz_yz, g_zz_y_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_y_y_xx[i] = -2.0 * g_0_y_yz_xx[i] * c_exps[i] + 4.0 * g_zz_y_yz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_y_xy[i] = -2.0 * g_0_y_yz_xy[i] * c_exps[i] + 4.0 * g_zz_y_yz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_y_xz[i] = -2.0 * g_0_y_yz_xz[i] * c_exps[i] + 4.0 * g_zz_y_yz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_y_yy[i] = -2.0 * g_0_y_yz_yy[i] * c_exps[i] + 4.0 * g_zz_y_yz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_y_yz[i] = -2.0 * g_0_y_yz_yz[i] * c_exps[i] + 4.0 * g_zz_y_yz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_y_zz[i] = -2.0 * g_0_y_yz_zz[i] * c_exps[i] + 4.0 * g_zz_y_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1434-1440)

    #pragma omp simd aligned(g_0_y_0_xx, g_0_y_0_xy, g_0_y_0_xz, g_0_y_0_yy, g_0_y_0_yz, g_0_y_0_zz, g_0_y_zz_xx, g_0_y_zz_xy, g_0_y_zz_xz, g_0_y_zz_yy, g_0_y_zz_yz, g_0_y_zz_zz, g_z_0_z_0_z_y_z_xx, g_z_0_z_0_z_y_z_xy, g_z_0_z_0_z_y_z_xz, g_z_0_z_0_z_y_z_yy, g_z_0_z_0_z_y_z_yz, g_z_0_z_0_z_y_z_zz, g_zz_y_0_xx, g_zz_y_0_xy, g_zz_y_0_xz, g_zz_y_0_yy, g_zz_y_0_yz, g_zz_y_0_zz, g_zz_y_zz_xx, g_zz_y_zz_xy, g_zz_y_zz_xz, g_zz_y_zz_yy, g_zz_y_zz_yz, g_zz_y_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_y_z_xx[i] = g_0_y_0_xx[i] - 2.0 * g_0_y_zz_xx[i] * c_exps[i] - 2.0 * g_zz_y_0_xx[i] * a_exp + 4.0 * g_zz_y_zz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_z_xy[i] = g_0_y_0_xy[i] - 2.0 * g_0_y_zz_xy[i] * c_exps[i] - 2.0 * g_zz_y_0_xy[i] * a_exp + 4.0 * g_zz_y_zz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_z_xz[i] = g_0_y_0_xz[i] - 2.0 * g_0_y_zz_xz[i] * c_exps[i] - 2.0 * g_zz_y_0_xz[i] * a_exp + 4.0 * g_zz_y_zz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_z_yy[i] = g_0_y_0_yy[i] - 2.0 * g_0_y_zz_yy[i] * c_exps[i] - 2.0 * g_zz_y_0_yy[i] * a_exp + 4.0 * g_zz_y_zz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_z_yz[i] = g_0_y_0_yz[i] - 2.0 * g_0_y_zz_yz[i] * c_exps[i] - 2.0 * g_zz_y_0_yz[i] * a_exp + 4.0 * g_zz_y_zz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_z_zz[i] = g_0_y_0_zz[i] - 2.0 * g_0_y_zz_zz[i] * c_exps[i] - 2.0 * g_zz_y_0_zz[i] * a_exp + 4.0 * g_zz_y_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1440-1446)

    #pragma omp simd aligned(g_0_z_xz_xx, g_0_z_xz_xy, g_0_z_xz_xz, g_0_z_xz_yy, g_0_z_xz_yz, g_0_z_xz_zz, g_z_0_z_0_z_z_x_xx, g_z_0_z_0_z_z_x_xy, g_z_0_z_0_z_z_x_xz, g_z_0_z_0_z_z_x_yy, g_z_0_z_0_z_z_x_yz, g_z_0_z_0_z_z_x_zz, g_zz_z_xz_xx, g_zz_z_xz_xy, g_zz_z_xz_xz, g_zz_z_xz_yy, g_zz_z_xz_yz, g_zz_z_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_z_x_xx[i] = -2.0 * g_0_z_xz_xx[i] * c_exps[i] + 4.0 * g_zz_z_xz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_x_xy[i] = -2.0 * g_0_z_xz_xy[i] * c_exps[i] + 4.0 * g_zz_z_xz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_x_xz[i] = -2.0 * g_0_z_xz_xz[i] * c_exps[i] + 4.0 * g_zz_z_xz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_x_yy[i] = -2.0 * g_0_z_xz_yy[i] * c_exps[i] + 4.0 * g_zz_z_xz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_x_yz[i] = -2.0 * g_0_z_xz_yz[i] * c_exps[i] + 4.0 * g_zz_z_xz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_x_zz[i] = -2.0 * g_0_z_xz_zz[i] * c_exps[i] + 4.0 * g_zz_z_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1446-1452)

    #pragma omp simd aligned(g_0_z_yz_xx, g_0_z_yz_xy, g_0_z_yz_xz, g_0_z_yz_yy, g_0_z_yz_yz, g_0_z_yz_zz, g_z_0_z_0_z_z_y_xx, g_z_0_z_0_z_z_y_xy, g_z_0_z_0_z_z_y_xz, g_z_0_z_0_z_z_y_yy, g_z_0_z_0_z_z_y_yz, g_z_0_z_0_z_z_y_zz, g_zz_z_yz_xx, g_zz_z_yz_xy, g_zz_z_yz_xz, g_zz_z_yz_yy, g_zz_z_yz_yz, g_zz_z_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_z_y_xx[i] = -2.0 * g_0_z_yz_xx[i] * c_exps[i] + 4.0 * g_zz_z_yz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_y_xy[i] = -2.0 * g_0_z_yz_xy[i] * c_exps[i] + 4.0 * g_zz_z_yz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_y_xz[i] = -2.0 * g_0_z_yz_xz[i] * c_exps[i] + 4.0 * g_zz_z_yz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_y_yy[i] = -2.0 * g_0_z_yz_yy[i] * c_exps[i] + 4.0 * g_zz_z_yz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_y_yz[i] = -2.0 * g_0_z_yz_yz[i] * c_exps[i] + 4.0 * g_zz_z_yz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_y_zz[i] = -2.0 * g_0_z_yz_zz[i] * c_exps[i] + 4.0 * g_zz_z_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (1452-1458)

    #pragma omp simd aligned(g_0_z_0_xx, g_0_z_0_xy, g_0_z_0_xz, g_0_z_0_yy, g_0_z_0_yz, g_0_z_0_zz, g_0_z_zz_xx, g_0_z_zz_xy, g_0_z_zz_xz, g_0_z_zz_yy, g_0_z_zz_yz, g_0_z_zz_zz, g_z_0_z_0_z_z_z_xx, g_z_0_z_0_z_z_z_xy, g_z_0_z_0_z_z_z_xz, g_z_0_z_0_z_z_z_yy, g_z_0_z_0_z_z_z_yz, g_z_0_z_0_z_z_z_zz, g_zz_z_0_xx, g_zz_z_0_xy, g_zz_z_0_xz, g_zz_z_0_yy, g_zz_z_0_yz, g_zz_z_0_zz, g_zz_z_zz_xx, g_zz_z_zz_xy, g_zz_z_zz_xz, g_zz_z_zz_yy, g_zz_z_zz_yz, g_zz_z_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_z_z_xx[i] = g_0_z_0_xx[i] - 2.0 * g_0_z_zz_xx[i] * c_exps[i] - 2.0 * g_zz_z_0_xx[i] * a_exp + 4.0 * g_zz_z_zz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_z_xy[i] = g_0_z_0_xy[i] - 2.0 * g_0_z_zz_xy[i] * c_exps[i] - 2.0 * g_zz_z_0_xy[i] * a_exp + 4.0 * g_zz_z_zz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_z_xz[i] = g_0_z_0_xz[i] - 2.0 * g_0_z_zz_xz[i] * c_exps[i] - 2.0 * g_zz_z_0_xz[i] * a_exp + 4.0 * g_zz_z_zz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_z_yy[i] = g_0_z_0_yy[i] - 2.0 * g_0_z_zz_yy[i] * c_exps[i] - 2.0 * g_zz_z_0_yy[i] * a_exp + 4.0 * g_zz_z_zz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_z_yz[i] = g_0_z_0_yz[i] - 2.0 * g_0_z_zz_yz[i] * c_exps[i] - 2.0 * g_zz_z_0_yz[i] * a_exp + 4.0 * g_zz_z_zz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_z_zz[i] = g_0_z_0_zz[i] - 2.0 * g_0_z_zz_zz[i] * c_exps[i] - 2.0 * g_zz_z_0_zz[i] * a_exp + 4.0 * g_zz_z_zz_zz[i] * a_exp * c_exps[i];
    }
}

} // t4c_geom namespace


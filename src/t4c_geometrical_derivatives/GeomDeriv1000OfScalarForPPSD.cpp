#include "GeomDeriv1000OfScalarForPPSD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1000_ppsd_0(CSimdArray<double>& buffer_1000_ppsd,
                     const CSimdArray<double>& buffer_spsd,
                     const CSimdArray<double>& buffer_dpsd,
                     const double a_exp) -> void
{
    const auto ndims = buffer_1000_ppsd.number_of_columns();

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

    /// Set up components of integrals buffer : buffer_1000_ppsd

    auto g_x_0_0_0_x_x_0_xx = buffer_1000_ppsd[0];

    auto g_x_0_0_0_x_x_0_xy = buffer_1000_ppsd[1];

    auto g_x_0_0_0_x_x_0_xz = buffer_1000_ppsd[2];

    auto g_x_0_0_0_x_x_0_yy = buffer_1000_ppsd[3];

    auto g_x_0_0_0_x_x_0_yz = buffer_1000_ppsd[4];

    auto g_x_0_0_0_x_x_0_zz = buffer_1000_ppsd[5];

    auto g_x_0_0_0_x_y_0_xx = buffer_1000_ppsd[6];

    auto g_x_0_0_0_x_y_0_xy = buffer_1000_ppsd[7];

    auto g_x_0_0_0_x_y_0_xz = buffer_1000_ppsd[8];

    auto g_x_0_0_0_x_y_0_yy = buffer_1000_ppsd[9];

    auto g_x_0_0_0_x_y_0_yz = buffer_1000_ppsd[10];

    auto g_x_0_0_0_x_y_0_zz = buffer_1000_ppsd[11];

    auto g_x_0_0_0_x_z_0_xx = buffer_1000_ppsd[12];

    auto g_x_0_0_0_x_z_0_xy = buffer_1000_ppsd[13];

    auto g_x_0_0_0_x_z_0_xz = buffer_1000_ppsd[14];

    auto g_x_0_0_0_x_z_0_yy = buffer_1000_ppsd[15];

    auto g_x_0_0_0_x_z_0_yz = buffer_1000_ppsd[16];

    auto g_x_0_0_0_x_z_0_zz = buffer_1000_ppsd[17];

    auto g_x_0_0_0_y_x_0_xx = buffer_1000_ppsd[18];

    auto g_x_0_0_0_y_x_0_xy = buffer_1000_ppsd[19];

    auto g_x_0_0_0_y_x_0_xz = buffer_1000_ppsd[20];

    auto g_x_0_0_0_y_x_0_yy = buffer_1000_ppsd[21];

    auto g_x_0_0_0_y_x_0_yz = buffer_1000_ppsd[22];

    auto g_x_0_0_0_y_x_0_zz = buffer_1000_ppsd[23];

    auto g_x_0_0_0_y_y_0_xx = buffer_1000_ppsd[24];

    auto g_x_0_0_0_y_y_0_xy = buffer_1000_ppsd[25];

    auto g_x_0_0_0_y_y_0_xz = buffer_1000_ppsd[26];

    auto g_x_0_0_0_y_y_0_yy = buffer_1000_ppsd[27];

    auto g_x_0_0_0_y_y_0_yz = buffer_1000_ppsd[28];

    auto g_x_0_0_0_y_y_0_zz = buffer_1000_ppsd[29];

    auto g_x_0_0_0_y_z_0_xx = buffer_1000_ppsd[30];

    auto g_x_0_0_0_y_z_0_xy = buffer_1000_ppsd[31];

    auto g_x_0_0_0_y_z_0_xz = buffer_1000_ppsd[32];

    auto g_x_0_0_0_y_z_0_yy = buffer_1000_ppsd[33];

    auto g_x_0_0_0_y_z_0_yz = buffer_1000_ppsd[34];

    auto g_x_0_0_0_y_z_0_zz = buffer_1000_ppsd[35];

    auto g_x_0_0_0_z_x_0_xx = buffer_1000_ppsd[36];

    auto g_x_0_0_0_z_x_0_xy = buffer_1000_ppsd[37];

    auto g_x_0_0_0_z_x_0_xz = buffer_1000_ppsd[38];

    auto g_x_0_0_0_z_x_0_yy = buffer_1000_ppsd[39];

    auto g_x_0_0_0_z_x_0_yz = buffer_1000_ppsd[40];

    auto g_x_0_0_0_z_x_0_zz = buffer_1000_ppsd[41];

    auto g_x_0_0_0_z_y_0_xx = buffer_1000_ppsd[42];

    auto g_x_0_0_0_z_y_0_xy = buffer_1000_ppsd[43];

    auto g_x_0_0_0_z_y_0_xz = buffer_1000_ppsd[44];

    auto g_x_0_0_0_z_y_0_yy = buffer_1000_ppsd[45];

    auto g_x_0_0_0_z_y_0_yz = buffer_1000_ppsd[46];

    auto g_x_0_0_0_z_y_0_zz = buffer_1000_ppsd[47];

    auto g_x_0_0_0_z_z_0_xx = buffer_1000_ppsd[48];

    auto g_x_0_0_0_z_z_0_xy = buffer_1000_ppsd[49];

    auto g_x_0_0_0_z_z_0_xz = buffer_1000_ppsd[50];

    auto g_x_0_0_0_z_z_0_yy = buffer_1000_ppsd[51];

    auto g_x_0_0_0_z_z_0_yz = buffer_1000_ppsd[52];

    auto g_x_0_0_0_z_z_0_zz = buffer_1000_ppsd[53];

    auto g_y_0_0_0_x_x_0_xx = buffer_1000_ppsd[54];

    auto g_y_0_0_0_x_x_0_xy = buffer_1000_ppsd[55];

    auto g_y_0_0_0_x_x_0_xz = buffer_1000_ppsd[56];

    auto g_y_0_0_0_x_x_0_yy = buffer_1000_ppsd[57];

    auto g_y_0_0_0_x_x_0_yz = buffer_1000_ppsd[58];

    auto g_y_0_0_0_x_x_0_zz = buffer_1000_ppsd[59];

    auto g_y_0_0_0_x_y_0_xx = buffer_1000_ppsd[60];

    auto g_y_0_0_0_x_y_0_xy = buffer_1000_ppsd[61];

    auto g_y_0_0_0_x_y_0_xz = buffer_1000_ppsd[62];

    auto g_y_0_0_0_x_y_0_yy = buffer_1000_ppsd[63];

    auto g_y_0_0_0_x_y_0_yz = buffer_1000_ppsd[64];

    auto g_y_0_0_0_x_y_0_zz = buffer_1000_ppsd[65];

    auto g_y_0_0_0_x_z_0_xx = buffer_1000_ppsd[66];

    auto g_y_0_0_0_x_z_0_xy = buffer_1000_ppsd[67];

    auto g_y_0_0_0_x_z_0_xz = buffer_1000_ppsd[68];

    auto g_y_0_0_0_x_z_0_yy = buffer_1000_ppsd[69];

    auto g_y_0_0_0_x_z_0_yz = buffer_1000_ppsd[70];

    auto g_y_0_0_0_x_z_0_zz = buffer_1000_ppsd[71];

    auto g_y_0_0_0_y_x_0_xx = buffer_1000_ppsd[72];

    auto g_y_0_0_0_y_x_0_xy = buffer_1000_ppsd[73];

    auto g_y_0_0_0_y_x_0_xz = buffer_1000_ppsd[74];

    auto g_y_0_0_0_y_x_0_yy = buffer_1000_ppsd[75];

    auto g_y_0_0_0_y_x_0_yz = buffer_1000_ppsd[76];

    auto g_y_0_0_0_y_x_0_zz = buffer_1000_ppsd[77];

    auto g_y_0_0_0_y_y_0_xx = buffer_1000_ppsd[78];

    auto g_y_0_0_0_y_y_0_xy = buffer_1000_ppsd[79];

    auto g_y_0_0_0_y_y_0_xz = buffer_1000_ppsd[80];

    auto g_y_0_0_0_y_y_0_yy = buffer_1000_ppsd[81];

    auto g_y_0_0_0_y_y_0_yz = buffer_1000_ppsd[82];

    auto g_y_0_0_0_y_y_0_zz = buffer_1000_ppsd[83];

    auto g_y_0_0_0_y_z_0_xx = buffer_1000_ppsd[84];

    auto g_y_0_0_0_y_z_0_xy = buffer_1000_ppsd[85];

    auto g_y_0_0_0_y_z_0_xz = buffer_1000_ppsd[86];

    auto g_y_0_0_0_y_z_0_yy = buffer_1000_ppsd[87];

    auto g_y_0_0_0_y_z_0_yz = buffer_1000_ppsd[88];

    auto g_y_0_0_0_y_z_0_zz = buffer_1000_ppsd[89];

    auto g_y_0_0_0_z_x_0_xx = buffer_1000_ppsd[90];

    auto g_y_0_0_0_z_x_0_xy = buffer_1000_ppsd[91];

    auto g_y_0_0_0_z_x_0_xz = buffer_1000_ppsd[92];

    auto g_y_0_0_0_z_x_0_yy = buffer_1000_ppsd[93];

    auto g_y_0_0_0_z_x_0_yz = buffer_1000_ppsd[94];

    auto g_y_0_0_0_z_x_0_zz = buffer_1000_ppsd[95];

    auto g_y_0_0_0_z_y_0_xx = buffer_1000_ppsd[96];

    auto g_y_0_0_0_z_y_0_xy = buffer_1000_ppsd[97];

    auto g_y_0_0_0_z_y_0_xz = buffer_1000_ppsd[98];

    auto g_y_0_0_0_z_y_0_yy = buffer_1000_ppsd[99];

    auto g_y_0_0_0_z_y_0_yz = buffer_1000_ppsd[100];

    auto g_y_0_0_0_z_y_0_zz = buffer_1000_ppsd[101];

    auto g_y_0_0_0_z_z_0_xx = buffer_1000_ppsd[102];

    auto g_y_0_0_0_z_z_0_xy = buffer_1000_ppsd[103];

    auto g_y_0_0_0_z_z_0_xz = buffer_1000_ppsd[104];

    auto g_y_0_0_0_z_z_0_yy = buffer_1000_ppsd[105];

    auto g_y_0_0_0_z_z_0_yz = buffer_1000_ppsd[106];

    auto g_y_0_0_0_z_z_0_zz = buffer_1000_ppsd[107];

    auto g_z_0_0_0_x_x_0_xx = buffer_1000_ppsd[108];

    auto g_z_0_0_0_x_x_0_xy = buffer_1000_ppsd[109];

    auto g_z_0_0_0_x_x_0_xz = buffer_1000_ppsd[110];

    auto g_z_0_0_0_x_x_0_yy = buffer_1000_ppsd[111];

    auto g_z_0_0_0_x_x_0_yz = buffer_1000_ppsd[112];

    auto g_z_0_0_0_x_x_0_zz = buffer_1000_ppsd[113];

    auto g_z_0_0_0_x_y_0_xx = buffer_1000_ppsd[114];

    auto g_z_0_0_0_x_y_0_xy = buffer_1000_ppsd[115];

    auto g_z_0_0_0_x_y_0_xz = buffer_1000_ppsd[116];

    auto g_z_0_0_0_x_y_0_yy = buffer_1000_ppsd[117];

    auto g_z_0_0_0_x_y_0_yz = buffer_1000_ppsd[118];

    auto g_z_0_0_0_x_y_0_zz = buffer_1000_ppsd[119];

    auto g_z_0_0_0_x_z_0_xx = buffer_1000_ppsd[120];

    auto g_z_0_0_0_x_z_0_xy = buffer_1000_ppsd[121];

    auto g_z_0_0_0_x_z_0_xz = buffer_1000_ppsd[122];

    auto g_z_0_0_0_x_z_0_yy = buffer_1000_ppsd[123];

    auto g_z_0_0_0_x_z_0_yz = buffer_1000_ppsd[124];

    auto g_z_0_0_0_x_z_0_zz = buffer_1000_ppsd[125];

    auto g_z_0_0_0_y_x_0_xx = buffer_1000_ppsd[126];

    auto g_z_0_0_0_y_x_0_xy = buffer_1000_ppsd[127];

    auto g_z_0_0_0_y_x_0_xz = buffer_1000_ppsd[128];

    auto g_z_0_0_0_y_x_0_yy = buffer_1000_ppsd[129];

    auto g_z_0_0_0_y_x_0_yz = buffer_1000_ppsd[130];

    auto g_z_0_0_0_y_x_0_zz = buffer_1000_ppsd[131];

    auto g_z_0_0_0_y_y_0_xx = buffer_1000_ppsd[132];

    auto g_z_0_0_0_y_y_0_xy = buffer_1000_ppsd[133];

    auto g_z_0_0_0_y_y_0_xz = buffer_1000_ppsd[134];

    auto g_z_0_0_0_y_y_0_yy = buffer_1000_ppsd[135];

    auto g_z_0_0_0_y_y_0_yz = buffer_1000_ppsd[136];

    auto g_z_0_0_0_y_y_0_zz = buffer_1000_ppsd[137];

    auto g_z_0_0_0_y_z_0_xx = buffer_1000_ppsd[138];

    auto g_z_0_0_0_y_z_0_xy = buffer_1000_ppsd[139];

    auto g_z_0_0_0_y_z_0_xz = buffer_1000_ppsd[140];

    auto g_z_0_0_0_y_z_0_yy = buffer_1000_ppsd[141];

    auto g_z_0_0_0_y_z_0_yz = buffer_1000_ppsd[142];

    auto g_z_0_0_0_y_z_0_zz = buffer_1000_ppsd[143];

    auto g_z_0_0_0_z_x_0_xx = buffer_1000_ppsd[144];

    auto g_z_0_0_0_z_x_0_xy = buffer_1000_ppsd[145];

    auto g_z_0_0_0_z_x_0_xz = buffer_1000_ppsd[146];

    auto g_z_0_0_0_z_x_0_yy = buffer_1000_ppsd[147];

    auto g_z_0_0_0_z_x_0_yz = buffer_1000_ppsd[148];

    auto g_z_0_0_0_z_x_0_zz = buffer_1000_ppsd[149];

    auto g_z_0_0_0_z_y_0_xx = buffer_1000_ppsd[150];

    auto g_z_0_0_0_z_y_0_xy = buffer_1000_ppsd[151];

    auto g_z_0_0_0_z_y_0_xz = buffer_1000_ppsd[152];

    auto g_z_0_0_0_z_y_0_yy = buffer_1000_ppsd[153];

    auto g_z_0_0_0_z_y_0_yz = buffer_1000_ppsd[154];

    auto g_z_0_0_0_z_y_0_zz = buffer_1000_ppsd[155];

    auto g_z_0_0_0_z_z_0_xx = buffer_1000_ppsd[156];

    auto g_z_0_0_0_z_z_0_xy = buffer_1000_ppsd[157];

    auto g_z_0_0_0_z_z_0_xz = buffer_1000_ppsd[158];

    auto g_z_0_0_0_z_z_0_yy = buffer_1000_ppsd[159];

    auto g_z_0_0_0_z_z_0_yz = buffer_1000_ppsd[160];

    auto g_z_0_0_0_z_z_0_zz = buffer_1000_ppsd[161];

    // integrals block (0-6)

    #pragma omp simd aligned(g_0_x_0_xx, g_0_x_0_xy, g_0_x_0_xz, g_0_x_0_yy, g_0_x_0_yz, g_0_x_0_zz, g_x_0_0_0_x_x_0_xx, g_x_0_0_0_x_x_0_xy, g_x_0_0_0_x_x_0_xz, g_x_0_0_0_x_x_0_yy, g_x_0_0_0_x_x_0_yz, g_x_0_0_0_x_x_0_zz, g_xx_x_0_xx, g_xx_x_0_xy, g_xx_x_0_xz, g_xx_x_0_yy, g_xx_x_0_yz, g_xx_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_x_0_xx[i] = -g_0_x_0_xx[i] + 2.0 * g_xx_x_0_xx[i] * a_exp;

        g_x_0_0_0_x_x_0_xy[i] = -g_0_x_0_xy[i] + 2.0 * g_xx_x_0_xy[i] * a_exp;

        g_x_0_0_0_x_x_0_xz[i] = -g_0_x_0_xz[i] + 2.0 * g_xx_x_0_xz[i] * a_exp;

        g_x_0_0_0_x_x_0_yy[i] = -g_0_x_0_yy[i] + 2.0 * g_xx_x_0_yy[i] * a_exp;

        g_x_0_0_0_x_x_0_yz[i] = -g_0_x_0_yz[i] + 2.0 * g_xx_x_0_yz[i] * a_exp;

        g_x_0_0_0_x_x_0_zz[i] = -g_0_x_0_zz[i] + 2.0 * g_xx_x_0_zz[i] * a_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_0_y_0_xx, g_0_y_0_xy, g_0_y_0_xz, g_0_y_0_yy, g_0_y_0_yz, g_0_y_0_zz, g_x_0_0_0_x_y_0_xx, g_x_0_0_0_x_y_0_xy, g_x_0_0_0_x_y_0_xz, g_x_0_0_0_x_y_0_yy, g_x_0_0_0_x_y_0_yz, g_x_0_0_0_x_y_0_zz, g_xx_y_0_xx, g_xx_y_0_xy, g_xx_y_0_xz, g_xx_y_0_yy, g_xx_y_0_yz, g_xx_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_y_0_xx[i] = -g_0_y_0_xx[i] + 2.0 * g_xx_y_0_xx[i] * a_exp;

        g_x_0_0_0_x_y_0_xy[i] = -g_0_y_0_xy[i] + 2.0 * g_xx_y_0_xy[i] * a_exp;

        g_x_0_0_0_x_y_0_xz[i] = -g_0_y_0_xz[i] + 2.0 * g_xx_y_0_xz[i] * a_exp;

        g_x_0_0_0_x_y_0_yy[i] = -g_0_y_0_yy[i] + 2.0 * g_xx_y_0_yy[i] * a_exp;

        g_x_0_0_0_x_y_0_yz[i] = -g_0_y_0_yz[i] + 2.0 * g_xx_y_0_yz[i] * a_exp;

        g_x_0_0_0_x_y_0_zz[i] = -g_0_y_0_zz[i] + 2.0 * g_xx_y_0_zz[i] * a_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_0_z_0_xx, g_0_z_0_xy, g_0_z_0_xz, g_0_z_0_yy, g_0_z_0_yz, g_0_z_0_zz, g_x_0_0_0_x_z_0_xx, g_x_0_0_0_x_z_0_xy, g_x_0_0_0_x_z_0_xz, g_x_0_0_0_x_z_0_yy, g_x_0_0_0_x_z_0_yz, g_x_0_0_0_x_z_0_zz, g_xx_z_0_xx, g_xx_z_0_xy, g_xx_z_0_xz, g_xx_z_0_yy, g_xx_z_0_yz, g_xx_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_z_0_xx[i] = -g_0_z_0_xx[i] + 2.0 * g_xx_z_0_xx[i] * a_exp;

        g_x_0_0_0_x_z_0_xy[i] = -g_0_z_0_xy[i] + 2.0 * g_xx_z_0_xy[i] * a_exp;

        g_x_0_0_0_x_z_0_xz[i] = -g_0_z_0_xz[i] + 2.0 * g_xx_z_0_xz[i] * a_exp;

        g_x_0_0_0_x_z_0_yy[i] = -g_0_z_0_yy[i] + 2.0 * g_xx_z_0_yy[i] * a_exp;

        g_x_0_0_0_x_z_0_yz[i] = -g_0_z_0_yz[i] + 2.0 * g_xx_z_0_yz[i] * a_exp;

        g_x_0_0_0_x_z_0_zz[i] = -g_0_z_0_zz[i] + 2.0 * g_xx_z_0_zz[i] * a_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_x_0_0_0_y_x_0_xx, g_x_0_0_0_y_x_0_xy, g_x_0_0_0_y_x_0_xz, g_x_0_0_0_y_x_0_yy, g_x_0_0_0_y_x_0_yz, g_x_0_0_0_y_x_0_zz, g_xy_x_0_xx, g_xy_x_0_xy, g_xy_x_0_xz, g_xy_x_0_yy, g_xy_x_0_yz, g_xy_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_x_0_xx[i] = 2.0 * g_xy_x_0_xx[i] * a_exp;

        g_x_0_0_0_y_x_0_xy[i] = 2.0 * g_xy_x_0_xy[i] * a_exp;

        g_x_0_0_0_y_x_0_xz[i] = 2.0 * g_xy_x_0_xz[i] * a_exp;

        g_x_0_0_0_y_x_0_yy[i] = 2.0 * g_xy_x_0_yy[i] * a_exp;

        g_x_0_0_0_y_x_0_yz[i] = 2.0 * g_xy_x_0_yz[i] * a_exp;

        g_x_0_0_0_y_x_0_zz[i] = 2.0 * g_xy_x_0_zz[i] * a_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_x_0_0_0_y_y_0_xx, g_x_0_0_0_y_y_0_xy, g_x_0_0_0_y_y_0_xz, g_x_0_0_0_y_y_0_yy, g_x_0_0_0_y_y_0_yz, g_x_0_0_0_y_y_0_zz, g_xy_y_0_xx, g_xy_y_0_xy, g_xy_y_0_xz, g_xy_y_0_yy, g_xy_y_0_yz, g_xy_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_y_0_xx[i] = 2.0 * g_xy_y_0_xx[i] * a_exp;

        g_x_0_0_0_y_y_0_xy[i] = 2.0 * g_xy_y_0_xy[i] * a_exp;

        g_x_0_0_0_y_y_0_xz[i] = 2.0 * g_xy_y_0_xz[i] * a_exp;

        g_x_0_0_0_y_y_0_yy[i] = 2.0 * g_xy_y_0_yy[i] * a_exp;

        g_x_0_0_0_y_y_0_yz[i] = 2.0 * g_xy_y_0_yz[i] * a_exp;

        g_x_0_0_0_y_y_0_zz[i] = 2.0 * g_xy_y_0_zz[i] * a_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_x_0_0_0_y_z_0_xx, g_x_0_0_0_y_z_0_xy, g_x_0_0_0_y_z_0_xz, g_x_0_0_0_y_z_0_yy, g_x_0_0_0_y_z_0_yz, g_x_0_0_0_y_z_0_zz, g_xy_z_0_xx, g_xy_z_0_xy, g_xy_z_0_xz, g_xy_z_0_yy, g_xy_z_0_yz, g_xy_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_z_0_xx[i] = 2.0 * g_xy_z_0_xx[i] * a_exp;

        g_x_0_0_0_y_z_0_xy[i] = 2.0 * g_xy_z_0_xy[i] * a_exp;

        g_x_0_0_0_y_z_0_xz[i] = 2.0 * g_xy_z_0_xz[i] * a_exp;

        g_x_0_0_0_y_z_0_yy[i] = 2.0 * g_xy_z_0_yy[i] * a_exp;

        g_x_0_0_0_y_z_0_yz[i] = 2.0 * g_xy_z_0_yz[i] * a_exp;

        g_x_0_0_0_y_z_0_zz[i] = 2.0 * g_xy_z_0_zz[i] * a_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_x_0_0_0_z_x_0_xx, g_x_0_0_0_z_x_0_xy, g_x_0_0_0_z_x_0_xz, g_x_0_0_0_z_x_0_yy, g_x_0_0_0_z_x_0_yz, g_x_0_0_0_z_x_0_zz, g_xz_x_0_xx, g_xz_x_0_xy, g_xz_x_0_xz, g_xz_x_0_yy, g_xz_x_0_yz, g_xz_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_x_0_xx[i] = 2.0 * g_xz_x_0_xx[i] * a_exp;

        g_x_0_0_0_z_x_0_xy[i] = 2.0 * g_xz_x_0_xy[i] * a_exp;

        g_x_0_0_0_z_x_0_xz[i] = 2.0 * g_xz_x_0_xz[i] * a_exp;

        g_x_0_0_0_z_x_0_yy[i] = 2.0 * g_xz_x_0_yy[i] * a_exp;

        g_x_0_0_0_z_x_0_yz[i] = 2.0 * g_xz_x_0_yz[i] * a_exp;

        g_x_0_0_0_z_x_0_zz[i] = 2.0 * g_xz_x_0_zz[i] * a_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_x_0_0_0_z_y_0_xx, g_x_0_0_0_z_y_0_xy, g_x_0_0_0_z_y_0_xz, g_x_0_0_0_z_y_0_yy, g_x_0_0_0_z_y_0_yz, g_x_0_0_0_z_y_0_zz, g_xz_y_0_xx, g_xz_y_0_xy, g_xz_y_0_xz, g_xz_y_0_yy, g_xz_y_0_yz, g_xz_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_y_0_xx[i] = 2.0 * g_xz_y_0_xx[i] * a_exp;

        g_x_0_0_0_z_y_0_xy[i] = 2.0 * g_xz_y_0_xy[i] * a_exp;

        g_x_0_0_0_z_y_0_xz[i] = 2.0 * g_xz_y_0_xz[i] * a_exp;

        g_x_0_0_0_z_y_0_yy[i] = 2.0 * g_xz_y_0_yy[i] * a_exp;

        g_x_0_0_0_z_y_0_yz[i] = 2.0 * g_xz_y_0_yz[i] * a_exp;

        g_x_0_0_0_z_y_0_zz[i] = 2.0 * g_xz_y_0_zz[i] * a_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_x_0_0_0_z_z_0_xx, g_x_0_0_0_z_z_0_xy, g_x_0_0_0_z_z_0_xz, g_x_0_0_0_z_z_0_yy, g_x_0_0_0_z_z_0_yz, g_x_0_0_0_z_z_0_zz, g_xz_z_0_xx, g_xz_z_0_xy, g_xz_z_0_xz, g_xz_z_0_yy, g_xz_z_0_yz, g_xz_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_z_0_xx[i] = 2.0 * g_xz_z_0_xx[i] * a_exp;

        g_x_0_0_0_z_z_0_xy[i] = 2.0 * g_xz_z_0_xy[i] * a_exp;

        g_x_0_0_0_z_z_0_xz[i] = 2.0 * g_xz_z_0_xz[i] * a_exp;

        g_x_0_0_0_z_z_0_yy[i] = 2.0 * g_xz_z_0_yy[i] * a_exp;

        g_x_0_0_0_z_z_0_yz[i] = 2.0 * g_xz_z_0_yz[i] * a_exp;

        g_x_0_0_0_z_z_0_zz[i] = 2.0 * g_xz_z_0_zz[i] * a_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_xy_x_0_xx, g_xy_x_0_xy, g_xy_x_0_xz, g_xy_x_0_yy, g_xy_x_0_yz, g_xy_x_0_zz, g_y_0_0_0_x_x_0_xx, g_y_0_0_0_x_x_0_xy, g_y_0_0_0_x_x_0_xz, g_y_0_0_0_x_x_0_yy, g_y_0_0_0_x_x_0_yz, g_y_0_0_0_x_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_x_0_xx[i] = 2.0 * g_xy_x_0_xx[i] * a_exp;

        g_y_0_0_0_x_x_0_xy[i] = 2.0 * g_xy_x_0_xy[i] * a_exp;

        g_y_0_0_0_x_x_0_xz[i] = 2.0 * g_xy_x_0_xz[i] * a_exp;

        g_y_0_0_0_x_x_0_yy[i] = 2.0 * g_xy_x_0_yy[i] * a_exp;

        g_y_0_0_0_x_x_0_yz[i] = 2.0 * g_xy_x_0_yz[i] * a_exp;

        g_y_0_0_0_x_x_0_zz[i] = 2.0 * g_xy_x_0_zz[i] * a_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_xy_y_0_xx, g_xy_y_0_xy, g_xy_y_0_xz, g_xy_y_0_yy, g_xy_y_0_yz, g_xy_y_0_zz, g_y_0_0_0_x_y_0_xx, g_y_0_0_0_x_y_0_xy, g_y_0_0_0_x_y_0_xz, g_y_0_0_0_x_y_0_yy, g_y_0_0_0_x_y_0_yz, g_y_0_0_0_x_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_y_0_xx[i] = 2.0 * g_xy_y_0_xx[i] * a_exp;

        g_y_0_0_0_x_y_0_xy[i] = 2.0 * g_xy_y_0_xy[i] * a_exp;

        g_y_0_0_0_x_y_0_xz[i] = 2.0 * g_xy_y_0_xz[i] * a_exp;

        g_y_0_0_0_x_y_0_yy[i] = 2.0 * g_xy_y_0_yy[i] * a_exp;

        g_y_0_0_0_x_y_0_yz[i] = 2.0 * g_xy_y_0_yz[i] * a_exp;

        g_y_0_0_0_x_y_0_zz[i] = 2.0 * g_xy_y_0_zz[i] * a_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_xy_z_0_xx, g_xy_z_0_xy, g_xy_z_0_xz, g_xy_z_0_yy, g_xy_z_0_yz, g_xy_z_0_zz, g_y_0_0_0_x_z_0_xx, g_y_0_0_0_x_z_0_xy, g_y_0_0_0_x_z_0_xz, g_y_0_0_0_x_z_0_yy, g_y_0_0_0_x_z_0_yz, g_y_0_0_0_x_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_z_0_xx[i] = 2.0 * g_xy_z_0_xx[i] * a_exp;

        g_y_0_0_0_x_z_0_xy[i] = 2.0 * g_xy_z_0_xy[i] * a_exp;

        g_y_0_0_0_x_z_0_xz[i] = 2.0 * g_xy_z_0_xz[i] * a_exp;

        g_y_0_0_0_x_z_0_yy[i] = 2.0 * g_xy_z_0_yy[i] * a_exp;

        g_y_0_0_0_x_z_0_yz[i] = 2.0 * g_xy_z_0_yz[i] * a_exp;

        g_y_0_0_0_x_z_0_zz[i] = 2.0 * g_xy_z_0_zz[i] * a_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_0_x_0_xx, g_0_x_0_xy, g_0_x_0_xz, g_0_x_0_yy, g_0_x_0_yz, g_0_x_0_zz, g_y_0_0_0_y_x_0_xx, g_y_0_0_0_y_x_0_xy, g_y_0_0_0_y_x_0_xz, g_y_0_0_0_y_x_0_yy, g_y_0_0_0_y_x_0_yz, g_y_0_0_0_y_x_0_zz, g_yy_x_0_xx, g_yy_x_0_xy, g_yy_x_0_xz, g_yy_x_0_yy, g_yy_x_0_yz, g_yy_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_x_0_xx[i] = -g_0_x_0_xx[i] + 2.0 * g_yy_x_0_xx[i] * a_exp;

        g_y_0_0_0_y_x_0_xy[i] = -g_0_x_0_xy[i] + 2.0 * g_yy_x_0_xy[i] * a_exp;

        g_y_0_0_0_y_x_0_xz[i] = -g_0_x_0_xz[i] + 2.0 * g_yy_x_0_xz[i] * a_exp;

        g_y_0_0_0_y_x_0_yy[i] = -g_0_x_0_yy[i] + 2.0 * g_yy_x_0_yy[i] * a_exp;

        g_y_0_0_0_y_x_0_yz[i] = -g_0_x_0_yz[i] + 2.0 * g_yy_x_0_yz[i] * a_exp;

        g_y_0_0_0_y_x_0_zz[i] = -g_0_x_0_zz[i] + 2.0 * g_yy_x_0_zz[i] * a_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_0_y_0_xx, g_0_y_0_xy, g_0_y_0_xz, g_0_y_0_yy, g_0_y_0_yz, g_0_y_0_zz, g_y_0_0_0_y_y_0_xx, g_y_0_0_0_y_y_0_xy, g_y_0_0_0_y_y_0_xz, g_y_0_0_0_y_y_0_yy, g_y_0_0_0_y_y_0_yz, g_y_0_0_0_y_y_0_zz, g_yy_y_0_xx, g_yy_y_0_xy, g_yy_y_0_xz, g_yy_y_0_yy, g_yy_y_0_yz, g_yy_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_y_0_xx[i] = -g_0_y_0_xx[i] + 2.0 * g_yy_y_0_xx[i] * a_exp;

        g_y_0_0_0_y_y_0_xy[i] = -g_0_y_0_xy[i] + 2.0 * g_yy_y_0_xy[i] * a_exp;

        g_y_0_0_0_y_y_0_xz[i] = -g_0_y_0_xz[i] + 2.0 * g_yy_y_0_xz[i] * a_exp;

        g_y_0_0_0_y_y_0_yy[i] = -g_0_y_0_yy[i] + 2.0 * g_yy_y_0_yy[i] * a_exp;

        g_y_0_0_0_y_y_0_yz[i] = -g_0_y_0_yz[i] + 2.0 * g_yy_y_0_yz[i] * a_exp;

        g_y_0_0_0_y_y_0_zz[i] = -g_0_y_0_zz[i] + 2.0 * g_yy_y_0_zz[i] * a_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_0_z_0_xx, g_0_z_0_xy, g_0_z_0_xz, g_0_z_0_yy, g_0_z_0_yz, g_0_z_0_zz, g_y_0_0_0_y_z_0_xx, g_y_0_0_0_y_z_0_xy, g_y_0_0_0_y_z_0_xz, g_y_0_0_0_y_z_0_yy, g_y_0_0_0_y_z_0_yz, g_y_0_0_0_y_z_0_zz, g_yy_z_0_xx, g_yy_z_0_xy, g_yy_z_0_xz, g_yy_z_0_yy, g_yy_z_0_yz, g_yy_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_z_0_xx[i] = -g_0_z_0_xx[i] + 2.0 * g_yy_z_0_xx[i] * a_exp;

        g_y_0_0_0_y_z_0_xy[i] = -g_0_z_0_xy[i] + 2.0 * g_yy_z_0_xy[i] * a_exp;

        g_y_0_0_0_y_z_0_xz[i] = -g_0_z_0_xz[i] + 2.0 * g_yy_z_0_xz[i] * a_exp;

        g_y_0_0_0_y_z_0_yy[i] = -g_0_z_0_yy[i] + 2.0 * g_yy_z_0_yy[i] * a_exp;

        g_y_0_0_0_y_z_0_yz[i] = -g_0_z_0_yz[i] + 2.0 * g_yy_z_0_yz[i] * a_exp;

        g_y_0_0_0_y_z_0_zz[i] = -g_0_z_0_zz[i] + 2.0 * g_yy_z_0_zz[i] * a_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_y_0_0_0_z_x_0_xx, g_y_0_0_0_z_x_0_xy, g_y_0_0_0_z_x_0_xz, g_y_0_0_0_z_x_0_yy, g_y_0_0_0_z_x_0_yz, g_y_0_0_0_z_x_0_zz, g_yz_x_0_xx, g_yz_x_0_xy, g_yz_x_0_xz, g_yz_x_0_yy, g_yz_x_0_yz, g_yz_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_x_0_xx[i] = 2.0 * g_yz_x_0_xx[i] * a_exp;

        g_y_0_0_0_z_x_0_xy[i] = 2.0 * g_yz_x_0_xy[i] * a_exp;

        g_y_0_0_0_z_x_0_xz[i] = 2.0 * g_yz_x_0_xz[i] * a_exp;

        g_y_0_0_0_z_x_0_yy[i] = 2.0 * g_yz_x_0_yy[i] * a_exp;

        g_y_0_0_0_z_x_0_yz[i] = 2.0 * g_yz_x_0_yz[i] * a_exp;

        g_y_0_0_0_z_x_0_zz[i] = 2.0 * g_yz_x_0_zz[i] * a_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_y_0_0_0_z_y_0_xx, g_y_0_0_0_z_y_0_xy, g_y_0_0_0_z_y_0_xz, g_y_0_0_0_z_y_0_yy, g_y_0_0_0_z_y_0_yz, g_y_0_0_0_z_y_0_zz, g_yz_y_0_xx, g_yz_y_0_xy, g_yz_y_0_xz, g_yz_y_0_yy, g_yz_y_0_yz, g_yz_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_y_0_xx[i] = 2.0 * g_yz_y_0_xx[i] * a_exp;

        g_y_0_0_0_z_y_0_xy[i] = 2.0 * g_yz_y_0_xy[i] * a_exp;

        g_y_0_0_0_z_y_0_xz[i] = 2.0 * g_yz_y_0_xz[i] * a_exp;

        g_y_0_0_0_z_y_0_yy[i] = 2.0 * g_yz_y_0_yy[i] * a_exp;

        g_y_0_0_0_z_y_0_yz[i] = 2.0 * g_yz_y_0_yz[i] * a_exp;

        g_y_0_0_0_z_y_0_zz[i] = 2.0 * g_yz_y_0_zz[i] * a_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_y_0_0_0_z_z_0_xx, g_y_0_0_0_z_z_0_xy, g_y_0_0_0_z_z_0_xz, g_y_0_0_0_z_z_0_yy, g_y_0_0_0_z_z_0_yz, g_y_0_0_0_z_z_0_zz, g_yz_z_0_xx, g_yz_z_0_xy, g_yz_z_0_xz, g_yz_z_0_yy, g_yz_z_0_yz, g_yz_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_z_0_xx[i] = 2.0 * g_yz_z_0_xx[i] * a_exp;

        g_y_0_0_0_z_z_0_xy[i] = 2.0 * g_yz_z_0_xy[i] * a_exp;

        g_y_0_0_0_z_z_0_xz[i] = 2.0 * g_yz_z_0_xz[i] * a_exp;

        g_y_0_0_0_z_z_0_yy[i] = 2.0 * g_yz_z_0_yy[i] * a_exp;

        g_y_0_0_0_z_z_0_yz[i] = 2.0 * g_yz_z_0_yz[i] * a_exp;

        g_y_0_0_0_z_z_0_zz[i] = 2.0 * g_yz_z_0_zz[i] * a_exp;
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_xz_x_0_xx, g_xz_x_0_xy, g_xz_x_0_xz, g_xz_x_0_yy, g_xz_x_0_yz, g_xz_x_0_zz, g_z_0_0_0_x_x_0_xx, g_z_0_0_0_x_x_0_xy, g_z_0_0_0_x_x_0_xz, g_z_0_0_0_x_x_0_yy, g_z_0_0_0_x_x_0_yz, g_z_0_0_0_x_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_x_0_xx[i] = 2.0 * g_xz_x_0_xx[i] * a_exp;

        g_z_0_0_0_x_x_0_xy[i] = 2.0 * g_xz_x_0_xy[i] * a_exp;

        g_z_0_0_0_x_x_0_xz[i] = 2.0 * g_xz_x_0_xz[i] * a_exp;

        g_z_0_0_0_x_x_0_yy[i] = 2.0 * g_xz_x_0_yy[i] * a_exp;

        g_z_0_0_0_x_x_0_yz[i] = 2.0 * g_xz_x_0_yz[i] * a_exp;

        g_z_0_0_0_x_x_0_zz[i] = 2.0 * g_xz_x_0_zz[i] * a_exp;
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_xz_y_0_xx, g_xz_y_0_xy, g_xz_y_0_xz, g_xz_y_0_yy, g_xz_y_0_yz, g_xz_y_0_zz, g_z_0_0_0_x_y_0_xx, g_z_0_0_0_x_y_0_xy, g_z_0_0_0_x_y_0_xz, g_z_0_0_0_x_y_0_yy, g_z_0_0_0_x_y_0_yz, g_z_0_0_0_x_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_y_0_xx[i] = 2.0 * g_xz_y_0_xx[i] * a_exp;

        g_z_0_0_0_x_y_0_xy[i] = 2.0 * g_xz_y_0_xy[i] * a_exp;

        g_z_0_0_0_x_y_0_xz[i] = 2.0 * g_xz_y_0_xz[i] * a_exp;

        g_z_0_0_0_x_y_0_yy[i] = 2.0 * g_xz_y_0_yy[i] * a_exp;

        g_z_0_0_0_x_y_0_yz[i] = 2.0 * g_xz_y_0_yz[i] * a_exp;

        g_z_0_0_0_x_y_0_zz[i] = 2.0 * g_xz_y_0_zz[i] * a_exp;
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_xz_z_0_xx, g_xz_z_0_xy, g_xz_z_0_xz, g_xz_z_0_yy, g_xz_z_0_yz, g_xz_z_0_zz, g_z_0_0_0_x_z_0_xx, g_z_0_0_0_x_z_0_xy, g_z_0_0_0_x_z_0_xz, g_z_0_0_0_x_z_0_yy, g_z_0_0_0_x_z_0_yz, g_z_0_0_0_x_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_z_0_xx[i] = 2.0 * g_xz_z_0_xx[i] * a_exp;

        g_z_0_0_0_x_z_0_xy[i] = 2.0 * g_xz_z_0_xy[i] * a_exp;

        g_z_0_0_0_x_z_0_xz[i] = 2.0 * g_xz_z_0_xz[i] * a_exp;

        g_z_0_0_0_x_z_0_yy[i] = 2.0 * g_xz_z_0_yy[i] * a_exp;

        g_z_0_0_0_x_z_0_yz[i] = 2.0 * g_xz_z_0_yz[i] * a_exp;

        g_z_0_0_0_x_z_0_zz[i] = 2.0 * g_xz_z_0_zz[i] * a_exp;
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_yz_x_0_xx, g_yz_x_0_xy, g_yz_x_0_xz, g_yz_x_0_yy, g_yz_x_0_yz, g_yz_x_0_zz, g_z_0_0_0_y_x_0_xx, g_z_0_0_0_y_x_0_xy, g_z_0_0_0_y_x_0_xz, g_z_0_0_0_y_x_0_yy, g_z_0_0_0_y_x_0_yz, g_z_0_0_0_y_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_x_0_xx[i] = 2.0 * g_yz_x_0_xx[i] * a_exp;

        g_z_0_0_0_y_x_0_xy[i] = 2.0 * g_yz_x_0_xy[i] * a_exp;

        g_z_0_0_0_y_x_0_xz[i] = 2.0 * g_yz_x_0_xz[i] * a_exp;

        g_z_0_0_0_y_x_0_yy[i] = 2.0 * g_yz_x_0_yy[i] * a_exp;

        g_z_0_0_0_y_x_0_yz[i] = 2.0 * g_yz_x_0_yz[i] * a_exp;

        g_z_0_0_0_y_x_0_zz[i] = 2.0 * g_yz_x_0_zz[i] * a_exp;
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_yz_y_0_xx, g_yz_y_0_xy, g_yz_y_0_xz, g_yz_y_0_yy, g_yz_y_0_yz, g_yz_y_0_zz, g_z_0_0_0_y_y_0_xx, g_z_0_0_0_y_y_0_xy, g_z_0_0_0_y_y_0_xz, g_z_0_0_0_y_y_0_yy, g_z_0_0_0_y_y_0_yz, g_z_0_0_0_y_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_y_0_xx[i] = 2.0 * g_yz_y_0_xx[i] * a_exp;

        g_z_0_0_0_y_y_0_xy[i] = 2.0 * g_yz_y_0_xy[i] * a_exp;

        g_z_0_0_0_y_y_0_xz[i] = 2.0 * g_yz_y_0_xz[i] * a_exp;

        g_z_0_0_0_y_y_0_yy[i] = 2.0 * g_yz_y_0_yy[i] * a_exp;

        g_z_0_0_0_y_y_0_yz[i] = 2.0 * g_yz_y_0_yz[i] * a_exp;

        g_z_0_0_0_y_y_0_zz[i] = 2.0 * g_yz_y_0_zz[i] * a_exp;
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_yz_z_0_xx, g_yz_z_0_xy, g_yz_z_0_xz, g_yz_z_0_yy, g_yz_z_0_yz, g_yz_z_0_zz, g_z_0_0_0_y_z_0_xx, g_z_0_0_0_y_z_0_xy, g_z_0_0_0_y_z_0_xz, g_z_0_0_0_y_z_0_yy, g_z_0_0_0_y_z_0_yz, g_z_0_0_0_y_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_z_0_xx[i] = 2.0 * g_yz_z_0_xx[i] * a_exp;

        g_z_0_0_0_y_z_0_xy[i] = 2.0 * g_yz_z_0_xy[i] * a_exp;

        g_z_0_0_0_y_z_0_xz[i] = 2.0 * g_yz_z_0_xz[i] * a_exp;

        g_z_0_0_0_y_z_0_yy[i] = 2.0 * g_yz_z_0_yy[i] * a_exp;

        g_z_0_0_0_y_z_0_yz[i] = 2.0 * g_yz_z_0_yz[i] * a_exp;

        g_z_0_0_0_y_z_0_zz[i] = 2.0 * g_yz_z_0_zz[i] * a_exp;
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_0_x_0_xx, g_0_x_0_xy, g_0_x_0_xz, g_0_x_0_yy, g_0_x_0_yz, g_0_x_0_zz, g_z_0_0_0_z_x_0_xx, g_z_0_0_0_z_x_0_xy, g_z_0_0_0_z_x_0_xz, g_z_0_0_0_z_x_0_yy, g_z_0_0_0_z_x_0_yz, g_z_0_0_0_z_x_0_zz, g_zz_x_0_xx, g_zz_x_0_xy, g_zz_x_0_xz, g_zz_x_0_yy, g_zz_x_0_yz, g_zz_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_x_0_xx[i] = -g_0_x_0_xx[i] + 2.0 * g_zz_x_0_xx[i] * a_exp;

        g_z_0_0_0_z_x_0_xy[i] = -g_0_x_0_xy[i] + 2.0 * g_zz_x_0_xy[i] * a_exp;

        g_z_0_0_0_z_x_0_xz[i] = -g_0_x_0_xz[i] + 2.0 * g_zz_x_0_xz[i] * a_exp;

        g_z_0_0_0_z_x_0_yy[i] = -g_0_x_0_yy[i] + 2.0 * g_zz_x_0_yy[i] * a_exp;

        g_z_0_0_0_z_x_0_yz[i] = -g_0_x_0_yz[i] + 2.0 * g_zz_x_0_yz[i] * a_exp;

        g_z_0_0_0_z_x_0_zz[i] = -g_0_x_0_zz[i] + 2.0 * g_zz_x_0_zz[i] * a_exp;
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_0_y_0_xx, g_0_y_0_xy, g_0_y_0_xz, g_0_y_0_yy, g_0_y_0_yz, g_0_y_0_zz, g_z_0_0_0_z_y_0_xx, g_z_0_0_0_z_y_0_xy, g_z_0_0_0_z_y_0_xz, g_z_0_0_0_z_y_0_yy, g_z_0_0_0_z_y_0_yz, g_z_0_0_0_z_y_0_zz, g_zz_y_0_xx, g_zz_y_0_xy, g_zz_y_0_xz, g_zz_y_0_yy, g_zz_y_0_yz, g_zz_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_y_0_xx[i] = -g_0_y_0_xx[i] + 2.0 * g_zz_y_0_xx[i] * a_exp;

        g_z_0_0_0_z_y_0_xy[i] = -g_0_y_0_xy[i] + 2.0 * g_zz_y_0_xy[i] * a_exp;

        g_z_0_0_0_z_y_0_xz[i] = -g_0_y_0_xz[i] + 2.0 * g_zz_y_0_xz[i] * a_exp;

        g_z_0_0_0_z_y_0_yy[i] = -g_0_y_0_yy[i] + 2.0 * g_zz_y_0_yy[i] * a_exp;

        g_z_0_0_0_z_y_0_yz[i] = -g_0_y_0_yz[i] + 2.0 * g_zz_y_0_yz[i] * a_exp;

        g_z_0_0_0_z_y_0_zz[i] = -g_0_y_0_zz[i] + 2.0 * g_zz_y_0_zz[i] * a_exp;
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_0_z_0_xx, g_0_z_0_xy, g_0_z_0_xz, g_0_z_0_yy, g_0_z_0_yz, g_0_z_0_zz, g_z_0_0_0_z_z_0_xx, g_z_0_0_0_z_z_0_xy, g_z_0_0_0_z_z_0_xz, g_z_0_0_0_z_z_0_yy, g_z_0_0_0_z_z_0_yz, g_z_0_0_0_z_z_0_zz, g_zz_z_0_xx, g_zz_z_0_xy, g_zz_z_0_xz, g_zz_z_0_yy, g_zz_z_0_yz, g_zz_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_z_0_xx[i] = -g_0_z_0_xx[i] + 2.0 * g_zz_z_0_xx[i] * a_exp;

        g_z_0_0_0_z_z_0_xy[i] = -g_0_z_0_xy[i] + 2.0 * g_zz_z_0_xy[i] * a_exp;

        g_z_0_0_0_z_z_0_xz[i] = -g_0_z_0_xz[i] + 2.0 * g_zz_z_0_xz[i] * a_exp;

        g_z_0_0_0_z_z_0_yy[i] = -g_0_z_0_yy[i] + 2.0 * g_zz_z_0_yy[i] * a_exp;

        g_z_0_0_0_z_z_0_yz[i] = -g_0_z_0_yz[i] + 2.0 * g_zz_z_0_yz[i] * a_exp;

        g_z_0_0_0_z_z_0_zz[i] = -g_0_z_0_zz[i] + 2.0 * g_zz_z_0_zz[i] * a_exp;
    }
}

} // t4c_geom namespace


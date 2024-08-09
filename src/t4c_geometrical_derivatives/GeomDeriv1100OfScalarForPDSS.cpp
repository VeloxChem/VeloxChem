#include "GeomDeriv1100OfScalarForPDSS.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1100_pdss_0(CSimdArray<double>& buffer_1100_pdss,
                     const CSimdArray<double>& buffer_spss,
                     const CSimdArray<double>& buffer_sfss,
                     const CSimdArray<double>& buffer_dpss,
                     const CSimdArray<double>& buffer_dfss,
                     const double a_exp,
                     const double b_exp) -> void
{
    const auto ndims = buffer_1100_pdss.number_of_columns();

    /// Set up components of auxilary buffer : buffer_spss

    auto g_0_x_0_0 = buffer_spss[0];

    auto g_0_y_0_0 = buffer_spss[1];

    auto g_0_z_0_0 = buffer_spss[2];

    /// Set up components of auxilary buffer : buffer_sfss

    auto g_0_xxx_0_0 = buffer_sfss[0];

    auto g_0_xxy_0_0 = buffer_sfss[1];

    auto g_0_xxz_0_0 = buffer_sfss[2];

    auto g_0_xyy_0_0 = buffer_sfss[3];

    auto g_0_xyz_0_0 = buffer_sfss[4];

    auto g_0_xzz_0_0 = buffer_sfss[5];

    auto g_0_yyy_0_0 = buffer_sfss[6];

    auto g_0_yyz_0_0 = buffer_sfss[7];

    auto g_0_yzz_0_0 = buffer_sfss[8];

    auto g_0_zzz_0_0 = buffer_sfss[9];

    /// Set up components of auxilary buffer : buffer_dpss

    auto g_xx_x_0_0 = buffer_dpss[0];

    auto g_xx_y_0_0 = buffer_dpss[1];

    auto g_xx_z_0_0 = buffer_dpss[2];

    auto g_xy_x_0_0 = buffer_dpss[3];

    auto g_xy_y_0_0 = buffer_dpss[4];

    auto g_xy_z_0_0 = buffer_dpss[5];

    auto g_xz_x_0_0 = buffer_dpss[6];

    auto g_xz_y_0_0 = buffer_dpss[7];

    auto g_xz_z_0_0 = buffer_dpss[8];

    auto g_yy_x_0_0 = buffer_dpss[9];

    auto g_yy_y_0_0 = buffer_dpss[10];

    auto g_yy_z_0_0 = buffer_dpss[11];

    auto g_yz_x_0_0 = buffer_dpss[12];

    auto g_yz_y_0_0 = buffer_dpss[13];

    auto g_yz_z_0_0 = buffer_dpss[14];

    auto g_zz_x_0_0 = buffer_dpss[15];

    auto g_zz_y_0_0 = buffer_dpss[16];

    auto g_zz_z_0_0 = buffer_dpss[17];

    /// Set up components of auxilary buffer : buffer_dfss

    auto g_xx_xxx_0_0 = buffer_dfss[0];

    auto g_xx_xxy_0_0 = buffer_dfss[1];

    auto g_xx_xxz_0_0 = buffer_dfss[2];

    auto g_xx_xyy_0_0 = buffer_dfss[3];

    auto g_xx_xyz_0_0 = buffer_dfss[4];

    auto g_xx_xzz_0_0 = buffer_dfss[5];

    auto g_xx_yyy_0_0 = buffer_dfss[6];

    auto g_xx_yyz_0_0 = buffer_dfss[7];

    auto g_xx_yzz_0_0 = buffer_dfss[8];

    auto g_xx_zzz_0_0 = buffer_dfss[9];

    auto g_xy_xxx_0_0 = buffer_dfss[10];

    auto g_xy_xxy_0_0 = buffer_dfss[11];

    auto g_xy_xxz_0_0 = buffer_dfss[12];

    auto g_xy_xyy_0_0 = buffer_dfss[13];

    auto g_xy_xyz_0_0 = buffer_dfss[14];

    auto g_xy_xzz_0_0 = buffer_dfss[15];

    auto g_xy_yyy_0_0 = buffer_dfss[16];

    auto g_xy_yyz_0_0 = buffer_dfss[17];

    auto g_xy_yzz_0_0 = buffer_dfss[18];

    auto g_xy_zzz_0_0 = buffer_dfss[19];

    auto g_xz_xxx_0_0 = buffer_dfss[20];

    auto g_xz_xxy_0_0 = buffer_dfss[21];

    auto g_xz_xxz_0_0 = buffer_dfss[22];

    auto g_xz_xyy_0_0 = buffer_dfss[23];

    auto g_xz_xyz_0_0 = buffer_dfss[24];

    auto g_xz_xzz_0_0 = buffer_dfss[25];

    auto g_xz_yyy_0_0 = buffer_dfss[26];

    auto g_xz_yyz_0_0 = buffer_dfss[27];

    auto g_xz_yzz_0_0 = buffer_dfss[28];

    auto g_xz_zzz_0_0 = buffer_dfss[29];

    auto g_yy_xxx_0_0 = buffer_dfss[30];

    auto g_yy_xxy_0_0 = buffer_dfss[31];

    auto g_yy_xxz_0_0 = buffer_dfss[32];

    auto g_yy_xyy_0_0 = buffer_dfss[33];

    auto g_yy_xyz_0_0 = buffer_dfss[34];

    auto g_yy_xzz_0_0 = buffer_dfss[35];

    auto g_yy_yyy_0_0 = buffer_dfss[36];

    auto g_yy_yyz_0_0 = buffer_dfss[37];

    auto g_yy_yzz_0_0 = buffer_dfss[38];

    auto g_yy_zzz_0_0 = buffer_dfss[39];

    auto g_yz_xxx_0_0 = buffer_dfss[40];

    auto g_yz_xxy_0_0 = buffer_dfss[41];

    auto g_yz_xxz_0_0 = buffer_dfss[42];

    auto g_yz_xyy_0_0 = buffer_dfss[43];

    auto g_yz_xyz_0_0 = buffer_dfss[44];

    auto g_yz_xzz_0_0 = buffer_dfss[45];

    auto g_yz_yyy_0_0 = buffer_dfss[46];

    auto g_yz_yyz_0_0 = buffer_dfss[47];

    auto g_yz_yzz_0_0 = buffer_dfss[48];

    auto g_yz_zzz_0_0 = buffer_dfss[49];

    auto g_zz_xxx_0_0 = buffer_dfss[50];

    auto g_zz_xxy_0_0 = buffer_dfss[51];

    auto g_zz_xxz_0_0 = buffer_dfss[52];

    auto g_zz_xyy_0_0 = buffer_dfss[53];

    auto g_zz_xyz_0_0 = buffer_dfss[54];

    auto g_zz_xzz_0_0 = buffer_dfss[55];

    auto g_zz_yyy_0_0 = buffer_dfss[56];

    auto g_zz_yyz_0_0 = buffer_dfss[57];

    auto g_zz_yzz_0_0 = buffer_dfss[58];

    auto g_zz_zzz_0_0 = buffer_dfss[59];

    /// Set up components of integrals buffer : buffer_1100_pdss

    auto g_x_x_0_0_x_xx_0_0 = buffer_1100_pdss[0];

    auto g_x_x_0_0_x_xy_0_0 = buffer_1100_pdss[1];

    auto g_x_x_0_0_x_xz_0_0 = buffer_1100_pdss[2];

    auto g_x_x_0_0_x_yy_0_0 = buffer_1100_pdss[3];

    auto g_x_x_0_0_x_yz_0_0 = buffer_1100_pdss[4];

    auto g_x_x_0_0_x_zz_0_0 = buffer_1100_pdss[5];

    auto g_x_x_0_0_y_xx_0_0 = buffer_1100_pdss[6];

    auto g_x_x_0_0_y_xy_0_0 = buffer_1100_pdss[7];

    auto g_x_x_0_0_y_xz_0_0 = buffer_1100_pdss[8];

    auto g_x_x_0_0_y_yy_0_0 = buffer_1100_pdss[9];

    auto g_x_x_0_0_y_yz_0_0 = buffer_1100_pdss[10];

    auto g_x_x_0_0_y_zz_0_0 = buffer_1100_pdss[11];

    auto g_x_x_0_0_z_xx_0_0 = buffer_1100_pdss[12];

    auto g_x_x_0_0_z_xy_0_0 = buffer_1100_pdss[13];

    auto g_x_x_0_0_z_xz_0_0 = buffer_1100_pdss[14];

    auto g_x_x_0_0_z_yy_0_0 = buffer_1100_pdss[15];

    auto g_x_x_0_0_z_yz_0_0 = buffer_1100_pdss[16];

    auto g_x_x_0_0_z_zz_0_0 = buffer_1100_pdss[17];

    auto g_x_y_0_0_x_xx_0_0 = buffer_1100_pdss[18];

    auto g_x_y_0_0_x_xy_0_0 = buffer_1100_pdss[19];

    auto g_x_y_0_0_x_xz_0_0 = buffer_1100_pdss[20];

    auto g_x_y_0_0_x_yy_0_0 = buffer_1100_pdss[21];

    auto g_x_y_0_0_x_yz_0_0 = buffer_1100_pdss[22];

    auto g_x_y_0_0_x_zz_0_0 = buffer_1100_pdss[23];

    auto g_x_y_0_0_y_xx_0_0 = buffer_1100_pdss[24];

    auto g_x_y_0_0_y_xy_0_0 = buffer_1100_pdss[25];

    auto g_x_y_0_0_y_xz_0_0 = buffer_1100_pdss[26];

    auto g_x_y_0_0_y_yy_0_0 = buffer_1100_pdss[27];

    auto g_x_y_0_0_y_yz_0_0 = buffer_1100_pdss[28];

    auto g_x_y_0_0_y_zz_0_0 = buffer_1100_pdss[29];

    auto g_x_y_0_0_z_xx_0_0 = buffer_1100_pdss[30];

    auto g_x_y_0_0_z_xy_0_0 = buffer_1100_pdss[31];

    auto g_x_y_0_0_z_xz_0_0 = buffer_1100_pdss[32];

    auto g_x_y_0_0_z_yy_0_0 = buffer_1100_pdss[33];

    auto g_x_y_0_0_z_yz_0_0 = buffer_1100_pdss[34];

    auto g_x_y_0_0_z_zz_0_0 = buffer_1100_pdss[35];

    auto g_x_z_0_0_x_xx_0_0 = buffer_1100_pdss[36];

    auto g_x_z_0_0_x_xy_0_0 = buffer_1100_pdss[37];

    auto g_x_z_0_0_x_xz_0_0 = buffer_1100_pdss[38];

    auto g_x_z_0_0_x_yy_0_0 = buffer_1100_pdss[39];

    auto g_x_z_0_0_x_yz_0_0 = buffer_1100_pdss[40];

    auto g_x_z_0_0_x_zz_0_0 = buffer_1100_pdss[41];

    auto g_x_z_0_0_y_xx_0_0 = buffer_1100_pdss[42];

    auto g_x_z_0_0_y_xy_0_0 = buffer_1100_pdss[43];

    auto g_x_z_0_0_y_xz_0_0 = buffer_1100_pdss[44];

    auto g_x_z_0_0_y_yy_0_0 = buffer_1100_pdss[45];

    auto g_x_z_0_0_y_yz_0_0 = buffer_1100_pdss[46];

    auto g_x_z_0_0_y_zz_0_0 = buffer_1100_pdss[47];

    auto g_x_z_0_0_z_xx_0_0 = buffer_1100_pdss[48];

    auto g_x_z_0_0_z_xy_0_0 = buffer_1100_pdss[49];

    auto g_x_z_0_0_z_xz_0_0 = buffer_1100_pdss[50];

    auto g_x_z_0_0_z_yy_0_0 = buffer_1100_pdss[51];

    auto g_x_z_0_0_z_yz_0_0 = buffer_1100_pdss[52];

    auto g_x_z_0_0_z_zz_0_0 = buffer_1100_pdss[53];

    auto g_y_x_0_0_x_xx_0_0 = buffer_1100_pdss[54];

    auto g_y_x_0_0_x_xy_0_0 = buffer_1100_pdss[55];

    auto g_y_x_0_0_x_xz_0_0 = buffer_1100_pdss[56];

    auto g_y_x_0_0_x_yy_0_0 = buffer_1100_pdss[57];

    auto g_y_x_0_0_x_yz_0_0 = buffer_1100_pdss[58];

    auto g_y_x_0_0_x_zz_0_0 = buffer_1100_pdss[59];

    auto g_y_x_0_0_y_xx_0_0 = buffer_1100_pdss[60];

    auto g_y_x_0_0_y_xy_0_0 = buffer_1100_pdss[61];

    auto g_y_x_0_0_y_xz_0_0 = buffer_1100_pdss[62];

    auto g_y_x_0_0_y_yy_0_0 = buffer_1100_pdss[63];

    auto g_y_x_0_0_y_yz_0_0 = buffer_1100_pdss[64];

    auto g_y_x_0_0_y_zz_0_0 = buffer_1100_pdss[65];

    auto g_y_x_0_0_z_xx_0_0 = buffer_1100_pdss[66];

    auto g_y_x_0_0_z_xy_0_0 = buffer_1100_pdss[67];

    auto g_y_x_0_0_z_xz_0_0 = buffer_1100_pdss[68];

    auto g_y_x_0_0_z_yy_0_0 = buffer_1100_pdss[69];

    auto g_y_x_0_0_z_yz_0_0 = buffer_1100_pdss[70];

    auto g_y_x_0_0_z_zz_0_0 = buffer_1100_pdss[71];

    auto g_y_y_0_0_x_xx_0_0 = buffer_1100_pdss[72];

    auto g_y_y_0_0_x_xy_0_0 = buffer_1100_pdss[73];

    auto g_y_y_0_0_x_xz_0_0 = buffer_1100_pdss[74];

    auto g_y_y_0_0_x_yy_0_0 = buffer_1100_pdss[75];

    auto g_y_y_0_0_x_yz_0_0 = buffer_1100_pdss[76];

    auto g_y_y_0_0_x_zz_0_0 = buffer_1100_pdss[77];

    auto g_y_y_0_0_y_xx_0_0 = buffer_1100_pdss[78];

    auto g_y_y_0_0_y_xy_0_0 = buffer_1100_pdss[79];

    auto g_y_y_0_0_y_xz_0_0 = buffer_1100_pdss[80];

    auto g_y_y_0_0_y_yy_0_0 = buffer_1100_pdss[81];

    auto g_y_y_0_0_y_yz_0_0 = buffer_1100_pdss[82];

    auto g_y_y_0_0_y_zz_0_0 = buffer_1100_pdss[83];

    auto g_y_y_0_0_z_xx_0_0 = buffer_1100_pdss[84];

    auto g_y_y_0_0_z_xy_0_0 = buffer_1100_pdss[85];

    auto g_y_y_0_0_z_xz_0_0 = buffer_1100_pdss[86];

    auto g_y_y_0_0_z_yy_0_0 = buffer_1100_pdss[87];

    auto g_y_y_0_0_z_yz_0_0 = buffer_1100_pdss[88];

    auto g_y_y_0_0_z_zz_0_0 = buffer_1100_pdss[89];

    auto g_y_z_0_0_x_xx_0_0 = buffer_1100_pdss[90];

    auto g_y_z_0_0_x_xy_0_0 = buffer_1100_pdss[91];

    auto g_y_z_0_0_x_xz_0_0 = buffer_1100_pdss[92];

    auto g_y_z_0_0_x_yy_0_0 = buffer_1100_pdss[93];

    auto g_y_z_0_0_x_yz_0_0 = buffer_1100_pdss[94];

    auto g_y_z_0_0_x_zz_0_0 = buffer_1100_pdss[95];

    auto g_y_z_0_0_y_xx_0_0 = buffer_1100_pdss[96];

    auto g_y_z_0_0_y_xy_0_0 = buffer_1100_pdss[97];

    auto g_y_z_0_0_y_xz_0_0 = buffer_1100_pdss[98];

    auto g_y_z_0_0_y_yy_0_0 = buffer_1100_pdss[99];

    auto g_y_z_0_0_y_yz_0_0 = buffer_1100_pdss[100];

    auto g_y_z_0_0_y_zz_0_0 = buffer_1100_pdss[101];

    auto g_y_z_0_0_z_xx_0_0 = buffer_1100_pdss[102];

    auto g_y_z_0_0_z_xy_0_0 = buffer_1100_pdss[103];

    auto g_y_z_0_0_z_xz_0_0 = buffer_1100_pdss[104];

    auto g_y_z_0_0_z_yy_0_0 = buffer_1100_pdss[105];

    auto g_y_z_0_0_z_yz_0_0 = buffer_1100_pdss[106];

    auto g_y_z_0_0_z_zz_0_0 = buffer_1100_pdss[107];

    auto g_z_x_0_0_x_xx_0_0 = buffer_1100_pdss[108];

    auto g_z_x_0_0_x_xy_0_0 = buffer_1100_pdss[109];

    auto g_z_x_0_0_x_xz_0_0 = buffer_1100_pdss[110];

    auto g_z_x_0_0_x_yy_0_0 = buffer_1100_pdss[111];

    auto g_z_x_0_0_x_yz_0_0 = buffer_1100_pdss[112];

    auto g_z_x_0_0_x_zz_0_0 = buffer_1100_pdss[113];

    auto g_z_x_0_0_y_xx_0_0 = buffer_1100_pdss[114];

    auto g_z_x_0_0_y_xy_0_0 = buffer_1100_pdss[115];

    auto g_z_x_0_0_y_xz_0_0 = buffer_1100_pdss[116];

    auto g_z_x_0_0_y_yy_0_0 = buffer_1100_pdss[117];

    auto g_z_x_0_0_y_yz_0_0 = buffer_1100_pdss[118];

    auto g_z_x_0_0_y_zz_0_0 = buffer_1100_pdss[119];

    auto g_z_x_0_0_z_xx_0_0 = buffer_1100_pdss[120];

    auto g_z_x_0_0_z_xy_0_0 = buffer_1100_pdss[121];

    auto g_z_x_0_0_z_xz_0_0 = buffer_1100_pdss[122];

    auto g_z_x_0_0_z_yy_0_0 = buffer_1100_pdss[123];

    auto g_z_x_0_0_z_yz_0_0 = buffer_1100_pdss[124];

    auto g_z_x_0_0_z_zz_0_0 = buffer_1100_pdss[125];

    auto g_z_y_0_0_x_xx_0_0 = buffer_1100_pdss[126];

    auto g_z_y_0_0_x_xy_0_0 = buffer_1100_pdss[127];

    auto g_z_y_0_0_x_xz_0_0 = buffer_1100_pdss[128];

    auto g_z_y_0_0_x_yy_0_0 = buffer_1100_pdss[129];

    auto g_z_y_0_0_x_yz_0_0 = buffer_1100_pdss[130];

    auto g_z_y_0_0_x_zz_0_0 = buffer_1100_pdss[131];

    auto g_z_y_0_0_y_xx_0_0 = buffer_1100_pdss[132];

    auto g_z_y_0_0_y_xy_0_0 = buffer_1100_pdss[133];

    auto g_z_y_0_0_y_xz_0_0 = buffer_1100_pdss[134];

    auto g_z_y_0_0_y_yy_0_0 = buffer_1100_pdss[135];

    auto g_z_y_0_0_y_yz_0_0 = buffer_1100_pdss[136];

    auto g_z_y_0_0_y_zz_0_0 = buffer_1100_pdss[137];

    auto g_z_y_0_0_z_xx_0_0 = buffer_1100_pdss[138];

    auto g_z_y_0_0_z_xy_0_0 = buffer_1100_pdss[139];

    auto g_z_y_0_0_z_xz_0_0 = buffer_1100_pdss[140];

    auto g_z_y_0_0_z_yy_0_0 = buffer_1100_pdss[141];

    auto g_z_y_0_0_z_yz_0_0 = buffer_1100_pdss[142];

    auto g_z_y_0_0_z_zz_0_0 = buffer_1100_pdss[143];

    auto g_z_z_0_0_x_xx_0_0 = buffer_1100_pdss[144];

    auto g_z_z_0_0_x_xy_0_0 = buffer_1100_pdss[145];

    auto g_z_z_0_0_x_xz_0_0 = buffer_1100_pdss[146];

    auto g_z_z_0_0_x_yy_0_0 = buffer_1100_pdss[147];

    auto g_z_z_0_0_x_yz_0_0 = buffer_1100_pdss[148];

    auto g_z_z_0_0_x_zz_0_0 = buffer_1100_pdss[149];

    auto g_z_z_0_0_y_xx_0_0 = buffer_1100_pdss[150];

    auto g_z_z_0_0_y_xy_0_0 = buffer_1100_pdss[151];

    auto g_z_z_0_0_y_xz_0_0 = buffer_1100_pdss[152];

    auto g_z_z_0_0_y_yy_0_0 = buffer_1100_pdss[153];

    auto g_z_z_0_0_y_yz_0_0 = buffer_1100_pdss[154];

    auto g_z_z_0_0_y_zz_0_0 = buffer_1100_pdss[155];

    auto g_z_z_0_0_z_xx_0_0 = buffer_1100_pdss[156];

    auto g_z_z_0_0_z_xy_0_0 = buffer_1100_pdss[157];

    auto g_z_z_0_0_z_xz_0_0 = buffer_1100_pdss[158];

    auto g_z_z_0_0_z_yy_0_0 = buffer_1100_pdss[159];

    auto g_z_z_0_0_z_yz_0_0 = buffer_1100_pdss[160];

    auto g_z_z_0_0_z_zz_0_0 = buffer_1100_pdss[161];

    // integrals block (0-6)

    #pragma omp simd aligned(g_0_x_0_0, g_0_xxx_0_0, g_0_xxy_0_0, g_0_xxz_0_0, g_0_xyy_0_0, g_0_xyz_0_0, g_0_xzz_0_0, g_0_y_0_0, g_0_z_0_0, g_x_x_0_0_x_xx_0_0, g_x_x_0_0_x_xy_0_0, g_x_x_0_0_x_xz_0_0, g_x_x_0_0_x_yy_0_0, g_x_x_0_0_x_yz_0_0, g_x_x_0_0_x_zz_0_0, g_xx_x_0_0, g_xx_xxx_0_0, g_xx_xxy_0_0, g_xx_xxz_0_0, g_xx_xyy_0_0, g_xx_xyz_0_0, g_xx_xzz_0_0, g_xx_y_0_0, g_xx_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_xx_0_0[i] = 2.0 * g_0_x_0_0[i] - 2.0 * g_0_xxx_0_0[i] * b_exp - 4.0 * g_xx_x_0_0[i] * a_exp + 4.0 * g_xx_xxx_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_x_xy_0_0[i] = g_0_y_0_0[i] - 2.0 * g_0_xxy_0_0[i] * b_exp - 2.0 * g_xx_y_0_0[i] * a_exp + 4.0 * g_xx_xxy_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_x_xz_0_0[i] = g_0_z_0_0[i] - 2.0 * g_0_xxz_0_0[i] * b_exp - 2.0 * g_xx_z_0_0[i] * a_exp + 4.0 * g_xx_xxz_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_x_yy_0_0[i] = -2.0 * g_0_xyy_0_0[i] * b_exp + 4.0 * g_xx_xyy_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_x_yz_0_0[i] = -2.0 * g_0_xyz_0_0[i] * b_exp + 4.0 * g_xx_xyz_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_x_zz_0_0[i] = -2.0 * g_0_xzz_0_0[i] * b_exp + 4.0 * g_xx_xzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_x_0_0_y_xx_0_0, g_x_x_0_0_y_xy_0_0, g_x_x_0_0_y_xz_0_0, g_x_x_0_0_y_yy_0_0, g_x_x_0_0_y_yz_0_0, g_x_x_0_0_y_zz_0_0, g_xy_x_0_0, g_xy_xxx_0_0, g_xy_xxy_0_0, g_xy_xxz_0_0, g_xy_xyy_0_0, g_xy_xyz_0_0, g_xy_xzz_0_0, g_xy_y_0_0, g_xy_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_xx_0_0[i] = -4.0 * g_xy_x_0_0[i] * a_exp + 4.0 * g_xy_xxx_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_y_xy_0_0[i] = -2.0 * g_xy_y_0_0[i] * a_exp + 4.0 * g_xy_xxy_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_y_xz_0_0[i] = -2.0 * g_xy_z_0_0[i] * a_exp + 4.0 * g_xy_xxz_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_y_yy_0_0[i] = 4.0 * g_xy_xyy_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_y_yz_0_0[i] = 4.0 * g_xy_xyz_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_y_zz_0_0[i] = 4.0 * g_xy_xzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_x_0_0_z_xx_0_0, g_x_x_0_0_z_xy_0_0, g_x_x_0_0_z_xz_0_0, g_x_x_0_0_z_yy_0_0, g_x_x_0_0_z_yz_0_0, g_x_x_0_0_z_zz_0_0, g_xz_x_0_0, g_xz_xxx_0_0, g_xz_xxy_0_0, g_xz_xxz_0_0, g_xz_xyy_0_0, g_xz_xyz_0_0, g_xz_xzz_0_0, g_xz_y_0_0, g_xz_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_xx_0_0[i] = -4.0 * g_xz_x_0_0[i] * a_exp + 4.0 * g_xz_xxx_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_z_xy_0_0[i] = -2.0 * g_xz_y_0_0[i] * a_exp + 4.0 * g_xz_xxy_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_z_xz_0_0[i] = -2.0 * g_xz_z_0_0[i] * a_exp + 4.0 * g_xz_xxz_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_z_yy_0_0[i] = 4.0 * g_xz_xyy_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_z_yz_0_0[i] = 4.0 * g_xz_xyz_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_z_zz_0_0[i] = 4.0 * g_xz_xzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_0_x_0_0, g_0_xxy_0_0, g_0_xyy_0_0, g_0_xyz_0_0, g_0_y_0_0, g_0_yyy_0_0, g_0_yyz_0_0, g_0_yzz_0_0, g_0_z_0_0, g_x_y_0_0_x_xx_0_0, g_x_y_0_0_x_xy_0_0, g_x_y_0_0_x_xz_0_0, g_x_y_0_0_x_yy_0_0, g_x_y_0_0_x_yz_0_0, g_x_y_0_0_x_zz_0_0, g_xx_x_0_0, g_xx_xxy_0_0, g_xx_xyy_0_0, g_xx_xyz_0_0, g_xx_y_0_0, g_xx_yyy_0_0, g_xx_yyz_0_0, g_xx_yzz_0_0, g_xx_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_xx_0_0[i] = -2.0 * g_0_xxy_0_0[i] * b_exp + 4.0 * g_xx_xxy_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_x_xy_0_0[i] = g_0_x_0_0[i] - 2.0 * g_0_xyy_0_0[i] * b_exp - 2.0 * g_xx_x_0_0[i] * a_exp + 4.0 * g_xx_xyy_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_x_xz_0_0[i] = -2.0 * g_0_xyz_0_0[i] * b_exp + 4.0 * g_xx_xyz_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_x_yy_0_0[i] = 2.0 * g_0_y_0_0[i] - 2.0 * g_0_yyy_0_0[i] * b_exp - 4.0 * g_xx_y_0_0[i] * a_exp + 4.0 * g_xx_yyy_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_x_yz_0_0[i] = g_0_z_0_0[i] - 2.0 * g_0_yyz_0_0[i] * b_exp - 2.0 * g_xx_z_0_0[i] * a_exp + 4.0 * g_xx_yyz_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_x_zz_0_0[i] = -2.0 * g_0_yzz_0_0[i] * b_exp + 4.0 * g_xx_yzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_x_y_0_0_y_xx_0_0, g_x_y_0_0_y_xy_0_0, g_x_y_0_0_y_xz_0_0, g_x_y_0_0_y_yy_0_0, g_x_y_0_0_y_yz_0_0, g_x_y_0_0_y_zz_0_0, g_xy_x_0_0, g_xy_xxy_0_0, g_xy_xyy_0_0, g_xy_xyz_0_0, g_xy_y_0_0, g_xy_yyy_0_0, g_xy_yyz_0_0, g_xy_yzz_0_0, g_xy_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_xx_0_0[i] = 4.0 * g_xy_xxy_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_y_xy_0_0[i] = -2.0 * g_xy_x_0_0[i] * a_exp + 4.0 * g_xy_xyy_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_y_xz_0_0[i] = 4.0 * g_xy_xyz_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_y_yy_0_0[i] = -4.0 * g_xy_y_0_0[i] * a_exp + 4.0 * g_xy_yyy_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_y_yz_0_0[i] = -2.0 * g_xy_z_0_0[i] * a_exp + 4.0 * g_xy_yyz_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_y_zz_0_0[i] = 4.0 * g_xy_yzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_x_y_0_0_z_xx_0_0, g_x_y_0_0_z_xy_0_0, g_x_y_0_0_z_xz_0_0, g_x_y_0_0_z_yy_0_0, g_x_y_0_0_z_yz_0_0, g_x_y_0_0_z_zz_0_0, g_xz_x_0_0, g_xz_xxy_0_0, g_xz_xyy_0_0, g_xz_xyz_0_0, g_xz_y_0_0, g_xz_yyy_0_0, g_xz_yyz_0_0, g_xz_yzz_0_0, g_xz_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_xx_0_0[i] = 4.0 * g_xz_xxy_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_z_xy_0_0[i] = -2.0 * g_xz_x_0_0[i] * a_exp + 4.0 * g_xz_xyy_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_z_xz_0_0[i] = 4.0 * g_xz_xyz_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_z_yy_0_0[i] = -4.0 * g_xz_y_0_0[i] * a_exp + 4.0 * g_xz_yyy_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_z_yz_0_0[i] = -2.0 * g_xz_z_0_0[i] * a_exp + 4.0 * g_xz_yyz_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_z_zz_0_0[i] = 4.0 * g_xz_yzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_0_x_0_0, g_0_xxz_0_0, g_0_xyz_0_0, g_0_xzz_0_0, g_0_y_0_0, g_0_yyz_0_0, g_0_yzz_0_0, g_0_z_0_0, g_0_zzz_0_0, g_x_z_0_0_x_xx_0_0, g_x_z_0_0_x_xy_0_0, g_x_z_0_0_x_xz_0_0, g_x_z_0_0_x_yy_0_0, g_x_z_0_0_x_yz_0_0, g_x_z_0_0_x_zz_0_0, g_xx_x_0_0, g_xx_xxz_0_0, g_xx_xyz_0_0, g_xx_xzz_0_0, g_xx_y_0_0, g_xx_yyz_0_0, g_xx_yzz_0_0, g_xx_z_0_0, g_xx_zzz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_xx_0_0[i] = -2.0 * g_0_xxz_0_0[i] * b_exp + 4.0 * g_xx_xxz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_x_xy_0_0[i] = -2.0 * g_0_xyz_0_0[i] * b_exp + 4.0 * g_xx_xyz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_x_xz_0_0[i] = g_0_x_0_0[i] - 2.0 * g_0_xzz_0_0[i] * b_exp - 2.0 * g_xx_x_0_0[i] * a_exp + 4.0 * g_xx_xzz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_x_yy_0_0[i] = -2.0 * g_0_yyz_0_0[i] * b_exp + 4.0 * g_xx_yyz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_x_yz_0_0[i] = g_0_y_0_0[i] - 2.0 * g_0_yzz_0_0[i] * b_exp - 2.0 * g_xx_y_0_0[i] * a_exp + 4.0 * g_xx_yzz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_x_zz_0_0[i] = 2.0 * g_0_z_0_0[i] - 2.0 * g_0_zzz_0_0[i] * b_exp - 4.0 * g_xx_z_0_0[i] * a_exp + 4.0 * g_xx_zzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_x_z_0_0_y_xx_0_0, g_x_z_0_0_y_xy_0_0, g_x_z_0_0_y_xz_0_0, g_x_z_0_0_y_yy_0_0, g_x_z_0_0_y_yz_0_0, g_x_z_0_0_y_zz_0_0, g_xy_x_0_0, g_xy_xxz_0_0, g_xy_xyz_0_0, g_xy_xzz_0_0, g_xy_y_0_0, g_xy_yyz_0_0, g_xy_yzz_0_0, g_xy_z_0_0, g_xy_zzz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_xx_0_0[i] = 4.0 * g_xy_xxz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_y_xy_0_0[i] = 4.0 * g_xy_xyz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_y_xz_0_0[i] = -2.0 * g_xy_x_0_0[i] * a_exp + 4.0 * g_xy_xzz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_y_yy_0_0[i] = 4.0 * g_xy_yyz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_y_yz_0_0[i] = -2.0 * g_xy_y_0_0[i] * a_exp + 4.0 * g_xy_yzz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_y_zz_0_0[i] = -4.0 * g_xy_z_0_0[i] * a_exp + 4.0 * g_xy_zzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_x_z_0_0_z_xx_0_0, g_x_z_0_0_z_xy_0_0, g_x_z_0_0_z_xz_0_0, g_x_z_0_0_z_yy_0_0, g_x_z_0_0_z_yz_0_0, g_x_z_0_0_z_zz_0_0, g_xz_x_0_0, g_xz_xxz_0_0, g_xz_xyz_0_0, g_xz_xzz_0_0, g_xz_y_0_0, g_xz_yyz_0_0, g_xz_yzz_0_0, g_xz_z_0_0, g_xz_zzz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_xx_0_0[i] = 4.0 * g_xz_xxz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_z_xy_0_0[i] = 4.0 * g_xz_xyz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_z_xz_0_0[i] = -2.0 * g_xz_x_0_0[i] * a_exp + 4.0 * g_xz_xzz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_z_yy_0_0[i] = 4.0 * g_xz_yyz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_z_yz_0_0[i] = -2.0 * g_xz_y_0_0[i] * a_exp + 4.0 * g_xz_yzz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_z_zz_0_0[i] = -4.0 * g_xz_z_0_0[i] * a_exp + 4.0 * g_xz_zzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_xy_x_0_0, g_xy_xxx_0_0, g_xy_xxy_0_0, g_xy_xxz_0_0, g_xy_xyy_0_0, g_xy_xyz_0_0, g_xy_xzz_0_0, g_xy_y_0_0, g_xy_z_0_0, g_y_x_0_0_x_xx_0_0, g_y_x_0_0_x_xy_0_0, g_y_x_0_0_x_xz_0_0, g_y_x_0_0_x_yy_0_0, g_y_x_0_0_x_yz_0_0, g_y_x_0_0_x_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_xx_0_0[i] = -4.0 * g_xy_x_0_0[i] * a_exp + 4.0 * g_xy_xxx_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_x_xy_0_0[i] = -2.0 * g_xy_y_0_0[i] * a_exp + 4.0 * g_xy_xxy_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_x_xz_0_0[i] = -2.0 * g_xy_z_0_0[i] * a_exp + 4.0 * g_xy_xxz_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_x_yy_0_0[i] = 4.0 * g_xy_xyy_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_x_yz_0_0[i] = 4.0 * g_xy_xyz_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_x_zz_0_0[i] = 4.0 * g_xy_xzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_0_x_0_0, g_0_xxx_0_0, g_0_xxy_0_0, g_0_xxz_0_0, g_0_xyy_0_0, g_0_xyz_0_0, g_0_xzz_0_0, g_0_y_0_0, g_0_z_0_0, g_y_x_0_0_y_xx_0_0, g_y_x_0_0_y_xy_0_0, g_y_x_0_0_y_xz_0_0, g_y_x_0_0_y_yy_0_0, g_y_x_0_0_y_yz_0_0, g_y_x_0_0_y_zz_0_0, g_yy_x_0_0, g_yy_xxx_0_0, g_yy_xxy_0_0, g_yy_xxz_0_0, g_yy_xyy_0_0, g_yy_xyz_0_0, g_yy_xzz_0_0, g_yy_y_0_0, g_yy_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_xx_0_0[i] = 2.0 * g_0_x_0_0[i] - 2.0 * g_0_xxx_0_0[i] * b_exp - 4.0 * g_yy_x_0_0[i] * a_exp + 4.0 * g_yy_xxx_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_y_xy_0_0[i] = g_0_y_0_0[i] - 2.0 * g_0_xxy_0_0[i] * b_exp - 2.0 * g_yy_y_0_0[i] * a_exp + 4.0 * g_yy_xxy_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_y_xz_0_0[i] = g_0_z_0_0[i] - 2.0 * g_0_xxz_0_0[i] * b_exp - 2.0 * g_yy_z_0_0[i] * a_exp + 4.0 * g_yy_xxz_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_y_yy_0_0[i] = -2.0 * g_0_xyy_0_0[i] * b_exp + 4.0 * g_yy_xyy_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_y_yz_0_0[i] = -2.0 * g_0_xyz_0_0[i] * b_exp + 4.0 * g_yy_xyz_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_y_zz_0_0[i] = -2.0 * g_0_xzz_0_0[i] * b_exp + 4.0 * g_yy_xzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_y_x_0_0_z_xx_0_0, g_y_x_0_0_z_xy_0_0, g_y_x_0_0_z_xz_0_0, g_y_x_0_0_z_yy_0_0, g_y_x_0_0_z_yz_0_0, g_y_x_0_0_z_zz_0_0, g_yz_x_0_0, g_yz_xxx_0_0, g_yz_xxy_0_0, g_yz_xxz_0_0, g_yz_xyy_0_0, g_yz_xyz_0_0, g_yz_xzz_0_0, g_yz_y_0_0, g_yz_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_xx_0_0[i] = -4.0 * g_yz_x_0_0[i] * a_exp + 4.0 * g_yz_xxx_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_z_xy_0_0[i] = -2.0 * g_yz_y_0_0[i] * a_exp + 4.0 * g_yz_xxy_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_z_xz_0_0[i] = -2.0 * g_yz_z_0_0[i] * a_exp + 4.0 * g_yz_xxz_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_z_yy_0_0[i] = 4.0 * g_yz_xyy_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_z_yz_0_0[i] = 4.0 * g_yz_xyz_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_z_zz_0_0[i] = 4.0 * g_yz_xzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_xy_x_0_0, g_xy_xxy_0_0, g_xy_xyy_0_0, g_xy_xyz_0_0, g_xy_y_0_0, g_xy_yyy_0_0, g_xy_yyz_0_0, g_xy_yzz_0_0, g_xy_z_0_0, g_y_y_0_0_x_xx_0_0, g_y_y_0_0_x_xy_0_0, g_y_y_0_0_x_xz_0_0, g_y_y_0_0_x_yy_0_0, g_y_y_0_0_x_yz_0_0, g_y_y_0_0_x_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_xx_0_0[i] = 4.0 * g_xy_xxy_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_x_xy_0_0[i] = -2.0 * g_xy_x_0_0[i] * a_exp + 4.0 * g_xy_xyy_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_x_xz_0_0[i] = 4.0 * g_xy_xyz_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_x_yy_0_0[i] = -4.0 * g_xy_y_0_0[i] * a_exp + 4.0 * g_xy_yyy_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_x_yz_0_0[i] = -2.0 * g_xy_z_0_0[i] * a_exp + 4.0 * g_xy_yyz_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_x_zz_0_0[i] = 4.0 * g_xy_yzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_0_x_0_0, g_0_xxy_0_0, g_0_xyy_0_0, g_0_xyz_0_0, g_0_y_0_0, g_0_yyy_0_0, g_0_yyz_0_0, g_0_yzz_0_0, g_0_z_0_0, g_y_y_0_0_y_xx_0_0, g_y_y_0_0_y_xy_0_0, g_y_y_0_0_y_xz_0_0, g_y_y_0_0_y_yy_0_0, g_y_y_0_0_y_yz_0_0, g_y_y_0_0_y_zz_0_0, g_yy_x_0_0, g_yy_xxy_0_0, g_yy_xyy_0_0, g_yy_xyz_0_0, g_yy_y_0_0, g_yy_yyy_0_0, g_yy_yyz_0_0, g_yy_yzz_0_0, g_yy_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_xx_0_0[i] = -2.0 * g_0_xxy_0_0[i] * b_exp + 4.0 * g_yy_xxy_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_y_xy_0_0[i] = g_0_x_0_0[i] - 2.0 * g_0_xyy_0_0[i] * b_exp - 2.0 * g_yy_x_0_0[i] * a_exp + 4.0 * g_yy_xyy_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_y_xz_0_0[i] = -2.0 * g_0_xyz_0_0[i] * b_exp + 4.0 * g_yy_xyz_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_y_yy_0_0[i] = 2.0 * g_0_y_0_0[i] - 2.0 * g_0_yyy_0_0[i] * b_exp - 4.0 * g_yy_y_0_0[i] * a_exp + 4.0 * g_yy_yyy_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_y_yz_0_0[i] = g_0_z_0_0[i] - 2.0 * g_0_yyz_0_0[i] * b_exp - 2.0 * g_yy_z_0_0[i] * a_exp + 4.0 * g_yy_yyz_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_y_zz_0_0[i] = -2.0 * g_0_yzz_0_0[i] * b_exp + 4.0 * g_yy_yzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_y_y_0_0_z_xx_0_0, g_y_y_0_0_z_xy_0_0, g_y_y_0_0_z_xz_0_0, g_y_y_0_0_z_yy_0_0, g_y_y_0_0_z_yz_0_0, g_y_y_0_0_z_zz_0_0, g_yz_x_0_0, g_yz_xxy_0_0, g_yz_xyy_0_0, g_yz_xyz_0_0, g_yz_y_0_0, g_yz_yyy_0_0, g_yz_yyz_0_0, g_yz_yzz_0_0, g_yz_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_xx_0_0[i] = 4.0 * g_yz_xxy_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_z_xy_0_0[i] = -2.0 * g_yz_x_0_0[i] * a_exp + 4.0 * g_yz_xyy_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_z_xz_0_0[i] = 4.0 * g_yz_xyz_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_z_yy_0_0[i] = -4.0 * g_yz_y_0_0[i] * a_exp + 4.0 * g_yz_yyy_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_z_yz_0_0[i] = -2.0 * g_yz_z_0_0[i] * a_exp + 4.0 * g_yz_yyz_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_z_zz_0_0[i] = 4.0 * g_yz_yzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_xy_x_0_0, g_xy_xxz_0_0, g_xy_xyz_0_0, g_xy_xzz_0_0, g_xy_y_0_0, g_xy_yyz_0_0, g_xy_yzz_0_0, g_xy_z_0_0, g_xy_zzz_0_0, g_y_z_0_0_x_xx_0_0, g_y_z_0_0_x_xy_0_0, g_y_z_0_0_x_xz_0_0, g_y_z_0_0_x_yy_0_0, g_y_z_0_0_x_yz_0_0, g_y_z_0_0_x_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_xx_0_0[i] = 4.0 * g_xy_xxz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_x_xy_0_0[i] = 4.0 * g_xy_xyz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_x_xz_0_0[i] = -2.0 * g_xy_x_0_0[i] * a_exp + 4.0 * g_xy_xzz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_x_yy_0_0[i] = 4.0 * g_xy_yyz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_x_yz_0_0[i] = -2.0 * g_xy_y_0_0[i] * a_exp + 4.0 * g_xy_yzz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_x_zz_0_0[i] = -4.0 * g_xy_z_0_0[i] * a_exp + 4.0 * g_xy_zzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_0_x_0_0, g_0_xxz_0_0, g_0_xyz_0_0, g_0_xzz_0_0, g_0_y_0_0, g_0_yyz_0_0, g_0_yzz_0_0, g_0_z_0_0, g_0_zzz_0_0, g_y_z_0_0_y_xx_0_0, g_y_z_0_0_y_xy_0_0, g_y_z_0_0_y_xz_0_0, g_y_z_0_0_y_yy_0_0, g_y_z_0_0_y_yz_0_0, g_y_z_0_0_y_zz_0_0, g_yy_x_0_0, g_yy_xxz_0_0, g_yy_xyz_0_0, g_yy_xzz_0_0, g_yy_y_0_0, g_yy_yyz_0_0, g_yy_yzz_0_0, g_yy_z_0_0, g_yy_zzz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_xx_0_0[i] = -2.0 * g_0_xxz_0_0[i] * b_exp + 4.0 * g_yy_xxz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_y_xy_0_0[i] = -2.0 * g_0_xyz_0_0[i] * b_exp + 4.0 * g_yy_xyz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_y_xz_0_0[i] = g_0_x_0_0[i] - 2.0 * g_0_xzz_0_0[i] * b_exp - 2.0 * g_yy_x_0_0[i] * a_exp + 4.0 * g_yy_xzz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_y_yy_0_0[i] = -2.0 * g_0_yyz_0_0[i] * b_exp + 4.0 * g_yy_yyz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_y_yz_0_0[i] = g_0_y_0_0[i] - 2.0 * g_0_yzz_0_0[i] * b_exp - 2.0 * g_yy_y_0_0[i] * a_exp + 4.0 * g_yy_yzz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_y_zz_0_0[i] = 2.0 * g_0_z_0_0[i] - 2.0 * g_0_zzz_0_0[i] * b_exp - 4.0 * g_yy_z_0_0[i] * a_exp + 4.0 * g_yy_zzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_y_z_0_0_z_xx_0_0, g_y_z_0_0_z_xy_0_0, g_y_z_0_0_z_xz_0_0, g_y_z_0_0_z_yy_0_0, g_y_z_0_0_z_yz_0_0, g_y_z_0_0_z_zz_0_0, g_yz_x_0_0, g_yz_xxz_0_0, g_yz_xyz_0_0, g_yz_xzz_0_0, g_yz_y_0_0, g_yz_yyz_0_0, g_yz_yzz_0_0, g_yz_z_0_0, g_yz_zzz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_xx_0_0[i] = 4.0 * g_yz_xxz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_z_xy_0_0[i] = 4.0 * g_yz_xyz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_z_xz_0_0[i] = -2.0 * g_yz_x_0_0[i] * a_exp + 4.0 * g_yz_xzz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_z_yy_0_0[i] = 4.0 * g_yz_yyz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_z_yz_0_0[i] = -2.0 * g_yz_y_0_0[i] * a_exp + 4.0 * g_yz_yzz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_z_zz_0_0[i] = -4.0 * g_yz_z_0_0[i] * a_exp + 4.0 * g_yz_zzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_xz_x_0_0, g_xz_xxx_0_0, g_xz_xxy_0_0, g_xz_xxz_0_0, g_xz_xyy_0_0, g_xz_xyz_0_0, g_xz_xzz_0_0, g_xz_y_0_0, g_xz_z_0_0, g_z_x_0_0_x_xx_0_0, g_z_x_0_0_x_xy_0_0, g_z_x_0_0_x_xz_0_0, g_z_x_0_0_x_yy_0_0, g_z_x_0_0_x_yz_0_0, g_z_x_0_0_x_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_xx_0_0[i] = -4.0 * g_xz_x_0_0[i] * a_exp + 4.0 * g_xz_xxx_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_x_xy_0_0[i] = -2.0 * g_xz_y_0_0[i] * a_exp + 4.0 * g_xz_xxy_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_x_xz_0_0[i] = -2.0 * g_xz_z_0_0[i] * a_exp + 4.0 * g_xz_xxz_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_x_yy_0_0[i] = 4.0 * g_xz_xyy_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_x_yz_0_0[i] = 4.0 * g_xz_xyz_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_x_zz_0_0[i] = 4.0 * g_xz_xzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_yz_x_0_0, g_yz_xxx_0_0, g_yz_xxy_0_0, g_yz_xxz_0_0, g_yz_xyy_0_0, g_yz_xyz_0_0, g_yz_xzz_0_0, g_yz_y_0_0, g_yz_z_0_0, g_z_x_0_0_y_xx_0_0, g_z_x_0_0_y_xy_0_0, g_z_x_0_0_y_xz_0_0, g_z_x_0_0_y_yy_0_0, g_z_x_0_0_y_yz_0_0, g_z_x_0_0_y_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_xx_0_0[i] = -4.0 * g_yz_x_0_0[i] * a_exp + 4.0 * g_yz_xxx_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_y_xy_0_0[i] = -2.0 * g_yz_y_0_0[i] * a_exp + 4.0 * g_yz_xxy_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_y_xz_0_0[i] = -2.0 * g_yz_z_0_0[i] * a_exp + 4.0 * g_yz_xxz_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_y_yy_0_0[i] = 4.0 * g_yz_xyy_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_y_yz_0_0[i] = 4.0 * g_yz_xyz_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_y_zz_0_0[i] = 4.0 * g_yz_xzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_0_x_0_0, g_0_xxx_0_0, g_0_xxy_0_0, g_0_xxz_0_0, g_0_xyy_0_0, g_0_xyz_0_0, g_0_xzz_0_0, g_0_y_0_0, g_0_z_0_0, g_z_x_0_0_z_xx_0_0, g_z_x_0_0_z_xy_0_0, g_z_x_0_0_z_xz_0_0, g_z_x_0_0_z_yy_0_0, g_z_x_0_0_z_yz_0_0, g_z_x_0_0_z_zz_0_0, g_zz_x_0_0, g_zz_xxx_0_0, g_zz_xxy_0_0, g_zz_xxz_0_0, g_zz_xyy_0_0, g_zz_xyz_0_0, g_zz_xzz_0_0, g_zz_y_0_0, g_zz_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_xx_0_0[i] = 2.0 * g_0_x_0_0[i] - 2.0 * g_0_xxx_0_0[i] * b_exp - 4.0 * g_zz_x_0_0[i] * a_exp + 4.0 * g_zz_xxx_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_z_xy_0_0[i] = g_0_y_0_0[i] - 2.0 * g_0_xxy_0_0[i] * b_exp - 2.0 * g_zz_y_0_0[i] * a_exp + 4.0 * g_zz_xxy_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_z_xz_0_0[i] = g_0_z_0_0[i] - 2.0 * g_0_xxz_0_0[i] * b_exp - 2.0 * g_zz_z_0_0[i] * a_exp + 4.0 * g_zz_xxz_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_z_yy_0_0[i] = -2.0 * g_0_xyy_0_0[i] * b_exp + 4.0 * g_zz_xyy_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_z_yz_0_0[i] = -2.0 * g_0_xyz_0_0[i] * b_exp + 4.0 * g_zz_xyz_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_z_zz_0_0[i] = -2.0 * g_0_xzz_0_0[i] * b_exp + 4.0 * g_zz_xzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_xz_x_0_0, g_xz_xxy_0_0, g_xz_xyy_0_0, g_xz_xyz_0_0, g_xz_y_0_0, g_xz_yyy_0_0, g_xz_yyz_0_0, g_xz_yzz_0_0, g_xz_z_0_0, g_z_y_0_0_x_xx_0_0, g_z_y_0_0_x_xy_0_0, g_z_y_0_0_x_xz_0_0, g_z_y_0_0_x_yy_0_0, g_z_y_0_0_x_yz_0_0, g_z_y_0_0_x_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_xx_0_0[i] = 4.0 * g_xz_xxy_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_x_xy_0_0[i] = -2.0 * g_xz_x_0_0[i] * a_exp + 4.0 * g_xz_xyy_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_x_xz_0_0[i] = 4.0 * g_xz_xyz_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_x_yy_0_0[i] = -4.0 * g_xz_y_0_0[i] * a_exp + 4.0 * g_xz_yyy_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_x_yz_0_0[i] = -2.0 * g_xz_z_0_0[i] * a_exp + 4.0 * g_xz_yyz_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_x_zz_0_0[i] = 4.0 * g_xz_yzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_yz_x_0_0, g_yz_xxy_0_0, g_yz_xyy_0_0, g_yz_xyz_0_0, g_yz_y_0_0, g_yz_yyy_0_0, g_yz_yyz_0_0, g_yz_yzz_0_0, g_yz_z_0_0, g_z_y_0_0_y_xx_0_0, g_z_y_0_0_y_xy_0_0, g_z_y_0_0_y_xz_0_0, g_z_y_0_0_y_yy_0_0, g_z_y_0_0_y_yz_0_0, g_z_y_0_0_y_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_xx_0_0[i] = 4.0 * g_yz_xxy_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_y_xy_0_0[i] = -2.0 * g_yz_x_0_0[i] * a_exp + 4.0 * g_yz_xyy_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_y_xz_0_0[i] = 4.0 * g_yz_xyz_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_y_yy_0_0[i] = -4.0 * g_yz_y_0_0[i] * a_exp + 4.0 * g_yz_yyy_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_y_yz_0_0[i] = -2.0 * g_yz_z_0_0[i] * a_exp + 4.0 * g_yz_yyz_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_y_zz_0_0[i] = 4.0 * g_yz_yzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_0_x_0_0, g_0_xxy_0_0, g_0_xyy_0_0, g_0_xyz_0_0, g_0_y_0_0, g_0_yyy_0_0, g_0_yyz_0_0, g_0_yzz_0_0, g_0_z_0_0, g_z_y_0_0_z_xx_0_0, g_z_y_0_0_z_xy_0_0, g_z_y_0_0_z_xz_0_0, g_z_y_0_0_z_yy_0_0, g_z_y_0_0_z_yz_0_0, g_z_y_0_0_z_zz_0_0, g_zz_x_0_0, g_zz_xxy_0_0, g_zz_xyy_0_0, g_zz_xyz_0_0, g_zz_y_0_0, g_zz_yyy_0_0, g_zz_yyz_0_0, g_zz_yzz_0_0, g_zz_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_xx_0_0[i] = -2.0 * g_0_xxy_0_0[i] * b_exp + 4.0 * g_zz_xxy_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_z_xy_0_0[i] = g_0_x_0_0[i] - 2.0 * g_0_xyy_0_0[i] * b_exp - 2.0 * g_zz_x_0_0[i] * a_exp + 4.0 * g_zz_xyy_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_z_xz_0_0[i] = -2.0 * g_0_xyz_0_0[i] * b_exp + 4.0 * g_zz_xyz_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_z_yy_0_0[i] = 2.0 * g_0_y_0_0[i] - 2.0 * g_0_yyy_0_0[i] * b_exp - 4.0 * g_zz_y_0_0[i] * a_exp + 4.0 * g_zz_yyy_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_z_yz_0_0[i] = g_0_z_0_0[i] - 2.0 * g_0_yyz_0_0[i] * b_exp - 2.0 * g_zz_z_0_0[i] * a_exp + 4.0 * g_zz_yyz_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_z_zz_0_0[i] = -2.0 * g_0_yzz_0_0[i] * b_exp + 4.0 * g_zz_yzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_xz_x_0_0, g_xz_xxz_0_0, g_xz_xyz_0_0, g_xz_xzz_0_0, g_xz_y_0_0, g_xz_yyz_0_0, g_xz_yzz_0_0, g_xz_z_0_0, g_xz_zzz_0_0, g_z_z_0_0_x_xx_0_0, g_z_z_0_0_x_xy_0_0, g_z_z_0_0_x_xz_0_0, g_z_z_0_0_x_yy_0_0, g_z_z_0_0_x_yz_0_0, g_z_z_0_0_x_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_xx_0_0[i] = 4.0 * g_xz_xxz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_x_xy_0_0[i] = 4.0 * g_xz_xyz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_x_xz_0_0[i] = -2.0 * g_xz_x_0_0[i] * a_exp + 4.0 * g_xz_xzz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_x_yy_0_0[i] = 4.0 * g_xz_yyz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_x_yz_0_0[i] = -2.0 * g_xz_y_0_0[i] * a_exp + 4.0 * g_xz_yzz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_x_zz_0_0[i] = -4.0 * g_xz_z_0_0[i] * a_exp + 4.0 * g_xz_zzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_yz_x_0_0, g_yz_xxz_0_0, g_yz_xyz_0_0, g_yz_xzz_0_0, g_yz_y_0_0, g_yz_yyz_0_0, g_yz_yzz_0_0, g_yz_z_0_0, g_yz_zzz_0_0, g_z_z_0_0_y_xx_0_0, g_z_z_0_0_y_xy_0_0, g_z_z_0_0_y_xz_0_0, g_z_z_0_0_y_yy_0_0, g_z_z_0_0_y_yz_0_0, g_z_z_0_0_y_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_xx_0_0[i] = 4.0 * g_yz_xxz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_y_xy_0_0[i] = 4.0 * g_yz_xyz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_y_xz_0_0[i] = -2.0 * g_yz_x_0_0[i] * a_exp + 4.0 * g_yz_xzz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_y_yy_0_0[i] = 4.0 * g_yz_yyz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_y_yz_0_0[i] = -2.0 * g_yz_y_0_0[i] * a_exp + 4.0 * g_yz_yzz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_y_zz_0_0[i] = -4.0 * g_yz_z_0_0[i] * a_exp + 4.0 * g_yz_zzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_0_x_0_0, g_0_xxz_0_0, g_0_xyz_0_0, g_0_xzz_0_0, g_0_y_0_0, g_0_yyz_0_0, g_0_yzz_0_0, g_0_z_0_0, g_0_zzz_0_0, g_z_z_0_0_z_xx_0_0, g_z_z_0_0_z_xy_0_0, g_z_z_0_0_z_xz_0_0, g_z_z_0_0_z_yy_0_0, g_z_z_0_0_z_yz_0_0, g_z_z_0_0_z_zz_0_0, g_zz_x_0_0, g_zz_xxz_0_0, g_zz_xyz_0_0, g_zz_xzz_0_0, g_zz_y_0_0, g_zz_yyz_0_0, g_zz_yzz_0_0, g_zz_z_0_0, g_zz_zzz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_xx_0_0[i] = -2.0 * g_0_xxz_0_0[i] * b_exp + 4.0 * g_zz_xxz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_z_xy_0_0[i] = -2.0 * g_0_xyz_0_0[i] * b_exp + 4.0 * g_zz_xyz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_z_xz_0_0[i] = g_0_x_0_0[i] - 2.0 * g_0_xzz_0_0[i] * b_exp - 2.0 * g_zz_x_0_0[i] * a_exp + 4.0 * g_zz_xzz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_z_yy_0_0[i] = -2.0 * g_0_yyz_0_0[i] * b_exp + 4.0 * g_zz_yyz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_z_yz_0_0[i] = g_0_y_0_0[i] - 2.0 * g_0_yzz_0_0[i] * b_exp - 2.0 * g_zz_y_0_0[i] * a_exp + 4.0 * g_zz_yzz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_z_zz_0_0[i] = 2.0 * g_0_z_0_0[i] - 2.0 * g_0_zzz_0_0[i] * b_exp - 4.0 * g_zz_z_0_0[i] * a_exp + 4.0 * g_zz_zzz_0_0[i] * a_exp * b_exp;
    }
}

} // t4c_geom namespace


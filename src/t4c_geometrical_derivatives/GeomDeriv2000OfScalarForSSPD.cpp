#include "GeomDeriv2000OfScalarForSSPD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom2000_sspd_0(CSimdArray<double>& buffer_2000_sspd,
                     const CSimdArray<double>& buffer_sspd,
                     const CSimdArray<double>& buffer_dspd,
                     const double a_exp) -> void
{
    const auto ndims = buffer_2000_sspd.number_of_columns();

    /// Set up components of auxilary buffer : buffer_sspd

    auto g_0_0_x_xx = buffer_sspd[0];

    auto g_0_0_x_xy = buffer_sspd[1];

    auto g_0_0_x_xz = buffer_sspd[2];

    auto g_0_0_x_yy = buffer_sspd[3];

    auto g_0_0_x_yz = buffer_sspd[4];

    auto g_0_0_x_zz = buffer_sspd[5];

    auto g_0_0_y_xx = buffer_sspd[6];

    auto g_0_0_y_xy = buffer_sspd[7];

    auto g_0_0_y_xz = buffer_sspd[8];

    auto g_0_0_y_yy = buffer_sspd[9];

    auto g_0_0_y_yz = buffer_sspd[10];

    auto g_0_0_y_zz = buffer_sspd[11];

    auto g_0_0_z_xx = buffer_sspd[12];

    auto g_0_0_z_xy = buffer_sspd[13];

    auto g_0_0_z_xz = buffer_sspd[14];

    auto g_0_0_z_yy = buffer_sspd[15];

    auto g_0_0_z_yz = buffer_sspd[16];

    auto g_0_0_z_zz = buffer_sspd[17];

    /// Set up components of auxilary buffer : buffer_dspd

    auto g_xx_0_x_xx = buffer_dspd[0];

    auto g_xx_0_x_xy = buffer_dspd[1];

    auto g_xx_0_x_xz = buffer_dspd[2];

    auto g_xx_0_x_yy = buffer_dspd[3];

    auto g_xx_0_x_yz = buffer_dspd[4];

    auto g_xx_0_x_zz = buffer_dspd[5];

    auto g_xx_0_y_xx = buffer_dspd[6];

    auto g_xx_0_y_xy = buffer_dspd[7];

    auto g_xx_0_y_xz = buffer_dspd[8];

    auto g_xx_0_y_yy = buffer_dspd[9];

    auto g_xx_0_y_yz = buffer_dspd[10];

    auto g_xx_0_y_zz = buffer_dspd[11];

    auto g_xx_0_z_xx = buffer_dspd[12];

    auto g_xx_0_z_xy = buffer_dspd[13];

    auto g_xx_0_z_xz = buffer_dspd[14];

    auto g_xx_0_z_yy = buffer_dspd[15];

    auto g_xx_0_z_yz = buffer_dspd[16];

    auto g_xx_0_z_zz = buffer_dspd[17];

    auto g_xy_0_x_xx = buffer_dspd[18];

    auto g_xy_0_x_xy = buffer_dspd[19];

    auto g_xy_0_x_xz = buffer_dspd[20];

    auto g_xy_0_x_yy = buffer_dspd[21];

    auto g_xy_0_x_yz = buffer_dspd[22];

    auto g_xy_0_x_zz = buffer_dspd[23];

    auto g_xy_0_y_xx = buffer_dspd[24];

    auto g_xy_0_y_xy = buffer_dspd[25];

    auto g_xy_0_y_xz = buffer_dspd[26];

    auto g_xy_0_y_yy = buffer_dspd[27];

    auto g_xy_0_y_yz = buffer_dspd[28];

    auto g_xy_0_y_zz = buffer_dspd[29];

    auto g_xy_0_z_xx = buffer_dspd[30];

    auto g_xy_0_z_xy = buffer_dspd[31];

    auto g_xy_0_z_xz = buffer_dspd[32];

    auto g_xy_0_z_yy = buffer_dspd[33];

    auto g_xy_0_z_yz = buffer_dspd[34];

    auto g_xy_0_z_zz = buffer_dspd[35];

    auto g_xz_0_x_xx = buffer_dspd[36];

    auto g_xz_0_x_xy = buffer_dspd[37];

    auto g_xz_0_x_xz = buffer_dspd[38];

    auto g_xz_0_x_yy = buffer_dspd[39];

    auto g_xz_0_x_yz = buffer_dspd[40];

    auto g_xz_0_x_zz = buffer_dspd[41];

    auto g_xz_0_y_xx = buffer_dspd[42];

    auto g_xz_0_y_xy = buffer_dspd[43];

    auto g_xz_0_y_xz = buffer_dspd[44];

    auto g_xz_0_y_yy = buffer_dspd[45];

    auto g_xz_0_y_yz = buffer_dspd[46];

    auto g_xz_0_y_zz = buffer_dspd[47];

    auto g_xz_0_z_xx = buffer_dspd[48];

    auto g_xz_0_z_xy = buffer_dspd[49];

    auto g_xz_0_z_xz = buffer_dspd[50];

    auto g_xz_0_z_yy = buffer_dspd[51];

    auto g_xz_0_z_yz = buffer_dspd[52];

    auto g_xz_0_z_zz = buffer_dspd[53];

    auto g_yy_0_x_xx = buffer_dspd[54];

    auto g_yy_0_x_xy = buffer_dspd[55];

    auto g_yy_0_x_xz = buffer_dspd[56];

    auto g_yy_0_x_yy = buffer_dspd[57];

    auto g_yy_0_x_yz = buffer_dspd[58];

    auto g_yy_0_x_zz = buffer_dspd[59];

    auto g_yy_0_y_xx = buffer_dspd[60];

    auto g_yy_0_y_xy = buffer_dspd[61];

    auto g_yy_0_y_xz = buffer_dspd[62];

    auto g_yy_0_y_yy = buffer_dspd[63];

    auto g_yy_0_y_yz = buffer_dspd[64];

    auto g_yy_0_y_zz = buffer_dspd[65];

    auto g_yy_0_z_xx = buffer_dspd[66];

    auto g_yy_0_z_xy = buffer_dspd[67];

    auto g_yy_0_z_xz = buffer_dspd[68];

    auto g_yy_0_z_yy = buffer_dspd[69];

    auto g_yy_0_z_yz = buffer_dspd[70];

    auto g_yy_0_z_zz = buffer_dspd[71];

    auto g_yz_0_x_xx = buffer_dspd[72];

    auto g_yz_0_x_xy = buffer_dspd[73];

    auto g_yz_0_x_xz = buffer_dspd[74];

    auto g_yz_0_x_yy = buffer_dspd[75];

    auto g_yz_0_x_yz = buffer_dspd[76];

    auto g_yz_0_x_zz = buffer_dspd[77];

    auto g_yz_0_y_xx = buffer_dspd[78];

    auto g_yz_0_y_xy = buffer_dspd[79];

    auto g_yz_0_y_xz = buffer_dspd[80];

    auto g_yz_0_y_yy = buffer_dspd[81];

    auto g_yz_0_y_yz = buffer_dspd[82];

    auto g_yz_0_y_zz = buffer_dspd[83];

    auto g_yz_0_z_xx = buffer_dspd[84];

    auto g_yz_0_z_xy = buffer_dspd[85];

    auto g_yz_0_z_xz = buffer_dspd[86];

    auto g_yz_0_z_yy = buffer_dspd[87];

    auto g_yz_0_z_yz = buffer_dspd[88];

    auto g_yz_0_z_zz = buffer_dspd[89];

    auto g_zz_0_x_xx = buffer_dspd[90];

    auto g_zz_0_x_xy = buffer_dspd[91];

    auto g_zz_0_x_xz = buffer_dspd[92];

    auto g_zz_0_x_yy = buffer_dspd[93];

    auto g_zz_0_x_yz = buffer_dspd[94];

    auto g_zz_0_x_zz = buffer_dspd[95];

    auto g_zz_0_y_xx = buffer_dspd[96];

    auto g_zz_0_y_xy = buffer_dspd[97];

    auto g_zz_0_y_xz = buffer_dspd[98];

    auto g_zz_0_y_yy = buffer_dspd[99];

    auto g_zz_0_y_yz = buffer_dspd[100];

    auto g_zz_0_y_zz = buffer_dspd[101];

    auto g_zz_0_z_xx = buffer_dspd[102];

    auto g_zz_0_z_xy = buffer_dspd[103];

    auto g_zz_0_z_xz = buffer_dspd[104];

    auto g_zz_0_z_yy = buffer_dspd[105];

    auto g_zz_0_z_yz = buffer_dspd[106];

    auto g_zz_0_z_zz = buffer_dspd[107];

    /// Set up components of integrals buffer : buffer_2000_sspd

    auto g_xx_0_0_0_0_0_x_xx = buffer_2000_sspd[0];

    auto g_xx_0_0_0_0_0_x_xy = buffer_2000_sspd[1];

    auto g_xx_0_0_0_0_0_x_xz = buffer_2000_sspd[2];

    auto g_xx_0_0_0_0_0_x_yy = buffer_2000_sspd[3];

    auto g_xx_0_0_0_0_0_x_yz = buffer_2000_sspd[4];

    auto g_xx_0_0_0_0_0_x_zz = buffer_2000_sspd[5];

    auto g_xx_0_0_0_0_0_y_xx = buffer_2000_sspd[6];

    auto g_xx_0_0_0_0_0_y_xy = buffer_2000_sspd[7];

    auto g_xx_0_0_0_0_0_y_xz = buffer_2000_sspd[8];

    auto g_xx_0_0_0_0_0_y_yy = buffer_2000_sspd[9];

    auto g_xx_0_0_0_0_0_y_yz = buffer_2000_sspd[10];

    auto g_xx_0_0_0_0_0_y_zz = buffer_2000_sspd[11];

    auto g_xx_0_0_0_0_0_z_xx = buffer_2000_sspd[12];

    auto g_xx_0_0_0_0_0_z_xy = buffer_2000_sspd[13];

    auto g_xx_0_0_0_0_0_z_xz = buffer_2000_sspd[14];

    auto g_xx_0_0_0_0_0_z_yy = buffer_2000_sspd[15];

    auto g_xx_0_0_0_0_0_z_yz = buffer_2000_sspd[16];

    auto g_xx_0_0_0_0_0_z_zz = buffer_2000_sspd[17];

    auto g_xy_0_0_0_0_0_x_xx = buffer_2000_sspd[18];

    auto g_xy_0_0_0_0_0_x_xy = buffer_2000_sspd[19];

    auto g_xy_0_0_0_0_0_x_xz = buffer_2000_sspd[20];

    auto g_xy_0_0_0_0_0_x_yy = buffer_2000_sspd[21];

    auto g_xy_0_0_0_0_0_x_yz = buffer_2000_sspd[22];

    auto g_xy_0_0_0_0_0_x_zz = buffer_2000_sspd[23];

    auto g_xy_0_0_0_0_0_y_xx = buffer_2000_sspd[24];

    auto g_xy_0_0_0_0_0_y_xy = buffer_2000_sspd[25];

    auto g_xy_0_0_0_0_0_y_xz = buffer_2000_sspd[26];

    auto g_xy_0_0_0_0_0_y_yy = buffer_2000_sspd[27];

    auto g_xy_0_0_0_0_0_y_yz = buffer_2000_sspd[28];

    auto g_xy_0_0_0_0_0_y_zz = buffer_2000_sspd[29];

    auto g_xy_0_0_0_0_0_z_xx = buffer_2000_sspd[30];

    auto g_xy_0_0_0_0_0_z_xy = buffer_2000_sspd[31];

    auto g_xy_0_0_0_0_0_z_xz = buffer_2000_sspd[32];

    auto g_xy_0_0_0_0_0_z_yy = buffer_2000_sspd[33];

    auto g_xy_0_0_0_0_0_z_yz = buffer_2000_sspd[34];

    auto g_xy_0_0_0_0_0_z_zz = buffer_2000_sspd[35];

    auto g_xz_0_0_0_0_0_x_xx = buffer_2000_sspd[36];

    auto g_xz_0_0_0_0_0_x_xy = buffer_2000_sspd[37];

    auto g_xz_0_0_0_0_0_x_xz = buffer_2000_sspd[38];

    auto g_xz_0_0_0_0_0_x_yy = buffer_2000_sspd[39];

    auto g_xz_0_0_0_0_0_x_yz = buffer_2000_sspd[40];

    auto g_xz_0_0_0_0_0_x_zz = buffer_2000_sspd[41];

    auto g_xz_0_0_0_0_0_y_xx = buffer_2000_sspd[42];

    auto g_xz_0_0_0_0_0_y_xy = buffer_2000_sspd[43];

    auto g_xz_0_0_0_0_0_y_xz = buffer_2000_sspd[44];

    auto g_xz_0_0_0_0_0_y_yy = buffer_2000_sspd[45];

    auto g_xz_0_0_0_0_0_y_yz = buffer_2000_sspd[46];

    auto g_xz_0_0_0_0_0_y_zz = buffer_2000_sspd[47];

    auto g_xz_0_0_0_0_0_z_xx = buffer_2000_sspd[48];

    auto g_xz_0_0_0_0_0_z_xy = buffer_2000_sspd[49];

    auto g_xz_0_0_0_0_0_z_xz = buffer_2000_sspd[50];

    auto g_xz_0_0_0_0_0_z_yy = buffer_2000_sspd[51];

    auto g_xz_0_0_0_0_0_z_yz = buffer_2000_sspd[52];

    auto g_xz_0_0_0_0_0_z_zz = buffer_2000_sspd[53];

    auto g_yy_0_0_0_0_0_x_xx = buffer_2000_sspd[54];

    auto g_yy_0_0_0_0_0_x_xy = buffer_2000_sspd[55];

    auto g_yy_0_0_0_0_0_x_xz = buffer_2000_sspd[56];

    auto g_yy_0_0_0_0_0_x_yy = buffer_2000_sspd[57];

    auto g_yy_0_0_0_0_0_x_yz = buffer_2000_sspd[58];

    auto g_yy_0_0_0_0_0_x_zz = buffer_2000_sspd[59];

    auto g_yy_0_0_0_0_0_y_xx = buffer_2000_sspd[60];

    auto g_yy_0_0_0_0_0_y_xy = buffer_2000_sspd[61];

    auto g_yy_0_0_0_0_0_y_xz = buffer_2000_sspd[62];

    auto g_yy_0_0_0_0_0_y_yy = buffer_2000_sspd[63];

    auto g_yy_0_0_0_0_0_y_yz = buffer_2000_sspd[64];

    auto g_yy_0_0_0_0_0_y_zz = buffer_2000_sspd[65];

    auto g_yy_0_0_0_0_0_z_xx = buffer_2000_sspd[66];

    auto g_yy_0_0_0_0_0_z_xy = buffer_2000_sspd[67];

    auto g_yy_0_0_0_0_0_z_xz = buffer_2000_sspd[68];

    auto g_yy_0_0_0_0_0_z_yy = buffer_2000_sspd[69];

    auto g_yy_0_0_0_0_0_z_yz = buffer_2000_sspd[70];

    auto g_yy_0_0_0_0_0_z_zz = buffer_2000_sspd[71];

    auto g_yz_0_0_0_0_0_x_xx = buffer_2000_sspd[72];

    auto g_yz_0_0_0_0_0_x_xy = buffer_2000_sspd[73];

    auto g_yz_0_0_0_0_0_x_xz = buffer_2000_sspd[74];

    auto g_yz_0_0_0_0_0_x_yy = buffer_2000_sspd[75];

    auto g_yz_0_0_0_0_0_x_yz = buffer_2000_sspd[76];

    auto g_yz_0_0_0_0_0_x_zz = buffer_2000_sspd[77];

    auto g_yz_0_0_0_0_0_y_xx = buffer_2000_sspd[78];

    auto g_yz_0_0_0_0_0_y_xy = buffer_2000_sspd[79];

    auto g_yz_0_0_0_0_0_y_xz = buffer_2000_sspd[80];

    auto g_yz_0_0_0_0_0_y_yy = buffer_2000_sspd[81];

    auto g_yz_0_0_0_0_0_y_yz = buffer_2000_sspd[82];

    auto g_yz_0_0_0_0_0_y_zz = buffer_2000_sspd[83];

    auto g_yz_0_0_0_0_0_z_xx = buffer_2000_sspd[84];

    auto g_yz_0_0_0_0_0_z_xy = buffer_2000_sspd[85];

    auto g_yz_0_0_0_0_0_z_xz = buffer_2000_sspd[86];

    auto g_yz_0_0_0_0_0_z_yy = buffer_2000_sspd[87];

    auto g_yz_0_0_0_0_0_z_yz = buffer_2000_sspd[88];

    auto g_yz_0_0_0_0_0_z_zz = buffer_2000_sspd[89];

    auto g_zz_0_0_0_0_0_x_xx = buffer_2000_sspd[90];

    auto g_zz_0_0_0_0_0_x_xy = buffer_2000_sspd[91];

    auto g_zz_0_0_0_0_0_x_xz = buffer_2000_sspd[92];

    auto g_zz_0_0_0_0_0_x_yy = buffer_2000_sspd[93];

    auto g_zz_0_0_0_0_0_x_yz = buffer_2000_sspd[94];

    auto g_zz_0_0_0_0_0_x_zz = buffer_2000_sspd[95];

    auto g_zz_0_0_0_0_0_y_xx = buffer_2000_sspd[96];

    auto g_zz_0_0_0_0_0_y_xy = buffer_2000_sspd[97];

    auto g_zz_0_0_0_0_0_y_xz = buffer_2000_sspd[98];

    auto g_zz_0_0_0_0_0_y_yy = buffer_2000_sspd[99];

    auto g_zz_0_0_0_0_0_y_yz = buffer_2000_sspd[100];

    auto g_zz_0_0_0_0_0_y_zz = buffer_2000_sspd[101];

    auto g_zz_0_0_0_0_0_z_xx = buffer_2000_sspd[102];

    auto g_zz_0_0_0_0_0_z_xy = buffer_2000_sspd[103];

    auto g_zz_0_0_0_0_0_z_xz = buffer_2000_sspd[104];

    auto g_zz_0_0_0_0_0_z_yy = buffer_2000_sspd[105];

    auto g_zz_0_0_0_0_0_z_yz = buffer_2000_sspd[106];

    auto g_zz_0_0_0_0_0_z_zz = buffer_2000_sspd[107];

    // integrals block (0-6)

    #pragma omp simd aligned(g_0_0_x_xx, g_0_0_x_xy, g_0_0_x_xz, g_0_0_x_yy, g_0_0_x_yz, g_0_0_x_zz, g_xx_0_0_0_0_0_x_xx, g_xx_0_0_0_0_0_x_xy, g_xx_0_0_0_0_0_x_xz, g_xx_0_0_0_0_0_x_yy, g_xx_0_0_0_0_0_x_yz, g_xx_0_0_0_0_0_x_zz, g_xx_0_x_xx, g_xx_0_x_xy, g_xx_0_x_xz, g_xx_0_x_yy, g_xx_0_x_yz, g_xx_0_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_0_x_xx[i] = -2.0 * g_0_0_x_xx[i] * a_exp + 4.0 * g_xx_0_x_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_x_xy[i] = -2.0 * g_0_0_x_xy[i] * a_exp + 4.0 * g_xx_0_x_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_x_xz[i] = -2.0 * g_0_0_x_xz[i] * a_exp + 4.0 * g_xx_0_x_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_x_yy[i] = -2.0 * g_0_0_x_yy[i] * a_exp + 4.0 * g_xx_0_x_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_x_yz[i] = -2.0 * g_0_0_x_yz[i] * a_exp + 4.0 * g_xx_0_x_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_x_zz[i] = -2.0 * g_0_0_x_zz[i] * a_exp + 4.0 * g_xx_0_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_0_0_y_xx, g_0_0_y_xy, g_0_0_y_xz, g_0_0_y_yy, g_0_0_y_yz, g_0_0_y_zz, g_xx_0_0_0_0_0_y_xx, g_xx_0_0_0_0_0_y_xy, g_xx_0_0_0_0_0_y_xz, g_xx_0_0_0_0_0_y_yy, g_xx_0_0_0_0_0_y_yz, g_xx_0_0_0_0_0_y_zz, g_xx_0_y_xx, g_xx_0_y_xy, g_xx_0_y_xz, g_xx_0_y_yy, g_xx_0_y_yz, g_xx_0_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_0_y_xx[i] = -2.0 * g_0_0_y_xx[i] * a_exp + 4.0 * g_xx_0_y_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_y_xy[i] = -2.0 * g_0_0_y_xy[i] * a_exp + 4.0 * g_xx_0_y_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_y_xz[i] = -2.0 * g_0_0_y_xz[i] * a_exp + 4.0 * g_xx_0_y_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_y_yy[i] = -2.0 * g_0_0_y_yy[i] * a_exp + 4.0 * g_xx_0_y_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_y_yz[i] = -2.0 * g_0_0_y_yz[i] * a_exp + 4.0 * g_xx_0_y_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_y_zz[i] = -2.0 * g_0_0_y_zz[i] * a_exp + 4.0 * g_xx_0_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_0_0_z_xx, g_0_0_z_xy, g_0_0_z_xz, g_0_0_z_yy, g_0_0_z_yz, g_0_0_z_zz, g_xx_0_0_0_0_0_z_xx, g_xx_0_0_0_0_0_z_xy, g_xx_0_0_0_0_0_z_xz, g_xx_0_0_0_0_0_z_yy, g_xx_0_0_0_0_0_z_yz, g_xx_0_0_0_0_0_z_zz, g_xx_0_z_xx, g_xx_0_z_xy, g_xx_0_z_xz, g_xx_0_z_yy, g_xx_0_z_yz, g_xx_0_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_0_z_xx[i] = -2.0 * g_0_0_z_xx[i] * a_exp + 4.0 * g_xx_0_z_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_z_xy[i] = -2.0 * g_0_0_z_xy[i] * a_exp + 4.0 * g_xx_0_z_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_z_xz[i] = -2.0 * g_0_0_z_xz[i] * a_exp + 4.0 * g_xx_0_z_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_z_yy[i] = -2.0 * g_0_0_z_yy[i] * a_exp + 4.0 * g_xx_0_z_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_z_yz[i] = -2.0 * g_0_0_z_yz[i] * a_exp + 4.0 * g_xx_0_z_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_z_zz[i] = -2.0 * g_0_0_z_zz[i] * a_exp + 4.0 * g_xx_0_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_xy_0_0_0_0_0_x_xx, g_xy_0_0_0_0_0_x_xy, g_xy_0_0_0_0_0_x_xz, g_xy_0_0_0_0_0_x_yy, g_xy_0_0_0_0_0_x_yz, g_xy_0_0_0_0_0_x_zz, g_xy_0_x_xx, g_xy_0_x_xy, g_xy_0_x_xz, g_xy_0_x_yy, g_xy_0_x_yz, g_xy_0_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_0_x_xx[i] = 4.0 * g_xy_0_x_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_x_xy[i] = 4.0 * g_xy_0_x_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_x_xz[i] = 4.0 * g_xy_0_x_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_x_yy[i] = 4.0 * g_xy_0_x_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_x_yz[i] = 4.0 * g_xy_0_x_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_x_zz[i] = 4.0 * g_xy_0_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_xy_0_0_0_0_0_y_xx, g_xy_0_0_0_0_0_y_xy, g_xy_0_0_0_0_0_y_xz, g_xy_0_0_0_0_0_y_yy, g_xy_0_0_0_0_0_y_yz, g_xy_0_0_0_0_0_y_zz, g_xy_0_y_xx, g_xy_0_y_xy, g_xy_0_y_xz, g_xy_0_y_yy, g_xy_0_y_yz, g_xy_0_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_0_y_xx[i] = 4.0 * g_xy_0_y_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_y_xy[i] = 4.0 * g_xy_0_y_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_y_xz[i] = 4.0 * g_xy_0_y_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_y_yy[i] = 4.0 * g_xy_0_y_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_y_yz[i] = 4.0 * g_xy_0_y_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_y_zz[i] = 4.0 * g_xy_0_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_xy_0_0_0_0_0_z_xx, g_xy_0_0_0_0_0_z_xy, g_xy_0_0_0_0_0_z_xz, g_xy_0_0_0_0_0_z_yy, g_xy_0_0_0_0_0_z_yz, g_xy_0_0_0_0_0_z_zz, g_xy_0_z_xx, g_xy_0_z_xy, g_xy_0_z_xz, g_xy_0_z_yy, g_xy_0_z_yz, g_xy_0_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_0_z_xx[i] = 4.0 * g_xy_0_z_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_z_xy[i] = 4.0 * g_xy_0_z_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_z_xz[i] = 4.0 * g_xy_0_z_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_z_yy[i] = 4.0 * g_xy_0_z_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_z_yz[i] = 4.0 * g_xy_0_z_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_z_zz[i] = 4.0 * g_xy_0_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_xz_0_0_0_0_0_x_xx, g_xz_0_0_0_0_0_x_xy, g_xz_0_0_0_0_0_x_xz, g_xz_0_0_0_0_0_x_yy, g_xz_0_0_0_0_0_x_yz, g_xz_0_0_0_0_0_x_zz, g_xz_0_x_xx, g_xz_0_x_xy, g_xz_0_x_xz, g_xz_0_x_yy, g_xz_0_x_yz, g_xz_0_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_0_x_xx[i] = 4.0 * g_xz_0_x_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_x_xy[i] = 4.0 * g_xz_0_x_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_x_xz[i] = 4.0 * g_xz_0_x_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_x_yy[i] = 4.0 * g_xz_0_x_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_x_yz[i] = 4.0 * g_xz_0_x_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_x_zz[i] = 4.0 * g_xz_0_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_xz_0_0_0_0_0_y_xx, g_xz_0_0_0_0_0_y_xy, g_xz_0_0_0_0_0_y_xz, g_xz_0_0_0_0_0_y_yy, g_xz_0_0_0_0_0_y_yz, g_xz_0_0_0_0_0_y_zz, g_xz_0_y_xx, g_xz_0_y_xy, g_xz_0_y_xz, g_xz_0_y_yy, g_xz_0_y_yz, g_xz_0_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_0_y_xx[i] = 4.0 * g_xz_0_y_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_y_xy[i] = 4.0 * g_xz_0_y_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_y_xz[i] = 4.0 * g_xz_0_y_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_y_yy[i] = 4.0 * g_xz_0_y_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_y_yz[i] = 4.0 * g_xz_0_y_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_y_zz[i] = 4.0 * g_xz_0_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_xz_0_0_0_0_0_z_xx, g_xz_0_0_0_0_0_z_xy, g_xz_0_0_0_0_0_z_xz, g_xz_0_0_0_0_0_z_yy, g_xz_0_0_0_0_0_z_yz, g_xz_0_0_0_0_0_z_zz, g_xz_0_z_xx, g_xz_0_z_xy, g_xz_0_z_xz, g_xz_0_z_yy, g_xz_0_z_yz, g_xz_0_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_0_z_xx[i] = 4.0 * g_xz_0_z_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_z_xy[i] = 4.0 * g_xz_0_z_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_z_xz[i] = 4.0 * g_xz_0_z_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_z_yy[i] = 4.0 * g_xz_0_z_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_z_yz[i] = 4.0 * g_xz_0_z_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_z_zz[i] = 4.0 * g_xz_0_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_0_0_x_xx, g_0_0_x_xy, g_0_0_x_xz, g_0_0_x_yy, g_0_0_x_yz, g_0_0_x_zz, g_yy_0_0_0_0_0_x_xx, g_yy_0_0_0_0_0_x_xy, g_yy_0_0_0_0_0_x_xz, g_yy_0_0_0_0_0_x_yy, g_yy_0_0_0_0_0_x_yz, g_yy_0_0_0_0_0_x_zz, g_yy_0_x_xx, g_yy_0_x_xy, g_yy_0_x_xz, g_yy_0_x_yy, g_yy_0_x_yz, g_yy_0_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_0_x_xx[i] = -2.0 * g_0_0_x_xx[i] * a_exp + 4.0 * g_yy_0_x_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_x_xy[i] = -2.0 * g_0_0_x_xy[i] * a_exp + 4.0 * g_yy_0_x_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_x_xz[i] = -2.0 * g_0_0_x_xz[i] * a_exp + 4.0 * g_yy_0_x_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_x_yy[i] = -2.0 * g_0_0_x_yy[i] * a_exp + 4.0 * g_yy_0_x_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_x_yz[i] = -2.0 * g_0_0_x_yz[i] * a_exp + 4.0 * g_yy_0_x_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_x_zz[i] = -2.0 * g_0_0_x_zz[i] * a_exp + 4.0 * g_yy_0_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_0_0_y_xx, g_0_0_y_xy, g_0_0_y_xz, g_0_0_y_yy, g_0_0_y_yz, g_0_0_y_zz, g_yy_0_0_0_0_0_y_xx, g_yy_0_0_0_0_0_y_xy, g_yy_0_0_0_0_0_y_xz, g_yy_0_0_0_0_0_y_yy, g_yy_0_0_0_0_0_y_yz, g_yy_0_0_0_0_0_y_zz, g_yy_0_y_xx, g_yy_0_y_xy, g_yy_0_y_xz, g_yy_0_y_yy, g_yy_0_y_yz, g_yy_0_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_0_y_xx[i] = -2.0 * g_0_0_y_xx[i] * a_exp + 4.0 * g_yy_0_y_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_y_xy[i] = -2.0 * g_0_0_y_xy[i] * a_exp + 4.0 * g_yy_0_y_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_y_xz[i] = -2.0 * g_0_0_y_xz[i] * a_exp + 4.0 * g_yy_0_y_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_y_yy[i] = -2.0 * g_0_0_y_yy[i] * a_exp + 4.0 * g_yy_0_y_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_y_yz[i] = -2.0 * g_0_0_y_yz[i] * a_exp + 4.0 * g_yy_0_y_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_y_zz[i] = -2.0 * g_0_0_y_zz[i] * a_exp + 4.0 * g_yy_0_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_0_0_z_xx, g_0_0_z_xy, g_0_0_z_xz, g_0_0_z_yy, g_0_0_z_yz, g_0_0_z_zz, g_yy_0_0_0_0_0_z_xx, g_yy_0_0_0_0_0_z_xy, g_yy_0_0_0_0_0_z_xz, g_yy_0_0_0_0_0_z_yy, g_yy_0_0_0_0_0_z_yz, g_yy_0_0_0_0_0_z_zz, g_yy_0_z_xx, g_yy_0_z_xy, g_yy_0_z_xz, g_yy_0_z_yy, g_yy_0_z_yz, g_yy_0_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_0_z_xx[i] = -2.0 * g_0_0_z_xx[i] * a_exp + 4.0 * g_yy_0_z_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_z_xy[i] = -2.0 * g_0_0_z_xy[i] * a_exp + 4.0 * g_yy_0_z_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_z_xz[i] = -2.0 * g_0_0_z_xz[i] * a_exp + 4.0 * g_yy_0_z_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_z_yy[i] = -2.0 * g_0_0_z_yy[i] * a_exp + 4.0 * g_yy_0_z_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_z_yz[i] = -2.0 * g_0_0_z_yz[i] * a_exp + 4.0 * g_yy_0_z_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_z_zz[i] = -2.0 * g_0_0_z_zz[i] * a_exp + 4.0 * g_yy_0_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_yz_0_0_0_0_0_x_xx, g_yz_0_0_0_0_0_x_xy, g_yz_0_0_0_0_0_x_xz, g_yz_0_0_0_0_0_x_yy, g_yz_0_0_0_0_0_x_yz, g_yz_0_0_0_0_0_x_zz, g_yz_0_x_xx, g_yz_0_x_xy, g_yz_0_x_xz, g_yz_0_x_yy, g_yz_0_x_yz, g_yz_0_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_0_x_xx[i] = 4.0 * g_yz_0_x_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_x_xy[i] = 4.0 * g_yz_0_x_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_x_xz[i] = 4.0 * g_yz_0_x_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_x_yy[i] = 4.0 * g_yz_0_x_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_x_yz[i] = 4.0 * g_yz_0_x_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_x_zz[i] = 4.0 * g_yz_0_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_yz_0_0_0_0_0_y_xx, g_yz_0_0_0_0_0_y_xy, g_yz_0_0_0_0_0_y_xz, g_yz_0_0_0_0_0_y_yy, g_yz_0_0_0_0_0_y_yz, g_yz_0_0_0_0_0_y_zz, g_yz_0_y_xx, g_yz_0_y_xy, g_yz_0_y_xz, g_yz_0_y_yy, g_yz_0_y_yz, g_yz_0_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_0_y_xx[i] = 4.0 * g_yz_0_y_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_y_xy[i] = 4.0 * g_yz_0_y_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_y_xz[i] = 4.0 * g_yz_0_y_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_y_yy[i] = 4.0 * g_yz_0_y_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_y_yz[i] = 4.0 * g_yz_0_y_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_y_zz[i] = 4.0 * g_yz_0_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_yz_0_0_0_0_0_z_xx, g_yz_0_0_0_0_0_z_xy, g_yz_0_0_0_0_0_z_xz, g_yz_0_0_0_0_0_z_yy, g_yz_0_0_0_0_0_z_yz, g_yz_0_0_0_0_0_z_zz, g_yz_0_z_xx, g_yz_0_z_xy, g_yz_0_z_xz, g_yz_0_z_yy, g_yz_0_z_yz, g_yz_0_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_0_z_xx[i] = 4.0 * g_yz_0_z_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_z_xy[i] = 4.0 * g_yz_0_z_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_z_xz[i] = 4.0 * g_yz_0_z_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_z_yy[i] = 4.0 * g_yz_0_z_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_z_yz[i] = 4.0 * g_yz_0_z_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_z_zz[i] = 4.0 * g_yz_0_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_0_0_x_xx, g_0_0_x_xy, g_0_0_x_xz, g_0_0_x_yy, g_0_0_x_yz, g_0_0_x_zz, g_zz_0_0_0_0_0_x_xx, g_zz_0_0_0_0_0_x_xy, g_zz_0_0_0_0_0_x_xz, g_zz_0_0_0_0_0_x_yy, g_zz_0_0_0_0_0_x_yz, g_zz_0_0_0_0_0_x_zz, g_zz_0_x_xx, g_zz_0_x_xy, g_zz_0_x_xz, g_zz_0_x_yy, g_zz_0_x_yz, g_zz_0_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_0_x_xx[i] = -2.0 * g_0_0_x_xx[i] * a_exp + 4.0 * g_zz_0_x_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_x_xy[i] = -2.0 * g_0_0_x_xy[i] * a_exp + 4.0 * g_zz_0_x_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_x_xz[i] = -2.0 * g_0_0_x_xz[i] * a_exp + 4.0 * g_zz_0_x_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_x_yy[i] = -2.0 * g_0_0_x_yy[i] * a_exp + 4.0 * g_zz_0_x_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_x_yz[i] = -2.0 * g_0_0_x_yz[i] * a_exp + 4.0 * g_zz_0_x_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_x_zz[i] = -2.0 * g_0_0_x_zz[i] * a_exp + 4.0 * g_zz_0_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_0_0_y_xx, g_0_0_y_xy, g_0_0_y_xz, g_0_0_y_yy, g_0_0_y_yz, g_0_0_y_zz, g_zz_0_0_0_0_0_y_xx, g_zz_0_0_0_0_0_y_xy, g_zz_0_0_0_0_0_y_xz, g_zz_0_0_0_0_0_y_yy, g_zz_0_0_0_0_0_y_yz, g_zz_0_0_0_0_0_y_zz, g_zz_0_y_xx, g_zz_0_y_xy, g_zz_0_y_xz, g_zz_0_y_yy, g_zz_0_y_yz, g_zz_0_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_0_y_xx[i] = -2.0 * g_0_0_y_xx[i] * a_exp + 4.0 * g_zz_0_y_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_y_xy[i] = -2.0 * g_0_0_y_xy[i] * a_exp + 4.0 * g_zz_0_y_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_y_xz[i] = -2.0 * g_0_0_y_xz[i] * a_exp + 4.0 * g_zz_0_y_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_y_yy[i] = -2.0 * g_0_0_y_yy[i] * a_exp + 4.0 * g_zz_0_y_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_y_yz[i] = -2.0 * g_0_0_y_yz[i] * a_exp + 4.0 * g_zz_0_y_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_y_zz[i] = -2.0 * g_0_0_y_zz[i] * a_exp + 4.0 * g_zz_0_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_0_0_z_xx, g_0_0_z_xy, g_0_0_z_xz, g_0_0_z_yy, g_0_0_z_yz, g_0_0_z_zz, g_zz_0_0_0_0_0_z_xx, g_zz_0_0_0_0_0_z_xy, g_zz_0_0_0_0_0_z_xz, g_zz_0_0_0_0_0_z_yy, g_zz_0_0_0_0_0_z_yz, g_zz_0_0_0_0_0_z_zz, g_zz_0_z_xx, g_zz_0_z_xy, g_zz_0_z_xz, g_zz_0_z_yy, g_zz_0_z_yz, g_zz_0_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_0_z_xx[i] = -2.0 * g_0_0_z_xx[i] * a_exp + 4.0 * g_zz_0_z_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_z_xy[i] = -2.0 * g_0_0_z_xy[i] * a_exp + 4.0 * g_zz_0_z_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_z_xz[i] = -2.0 * g_0_0_z_xz[i] * a_exp + 4.0 * g_zz_0_z_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_z_yy[i] = -2.0 * g_0_0_z_yy[i] * a_exp + 4.0 * g_zz_0_z_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_z_yz[i] = -2.0 * g_0_0_z_yz[i] * a_exp + 4.0 * g_zz_0_z_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_z_zz[i] = -2.0 * g_0_0_z_zz[i] * a_exp + 4.0 * g_zz_0_z_zz[i] * a_exp * a_exp;
    }
}

} // t4c_geom namespace


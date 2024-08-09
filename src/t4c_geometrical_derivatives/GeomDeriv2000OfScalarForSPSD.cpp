#include "GeomDeriv2000OfScalarForSPSD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom2000_spsd_0(CSimdArray<double>& buffer_2000_spsd,
                     const CSimdArray<double>& buffer_spsd,
                     const CSimdArray<double>& buffer_dpsd,
                     const double a_exp) -> void
{
    const auto ndims = buffer_2000_spsd.number_of_columns();

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

    /// Set up components of integrals buffer : buffer_2000_spsd

    auto g_xx_0_0_0_0_x_0_xx = buffer_2000_spsd[0];

    auto g_xx_0_0_0_0_x_0_xy = buffer_2000_spsd[1];

    auto g_xx_0_0_0_0_x_0_xz = buffer_2000_spsd[2];

    auto g_xx_0_0_0_0_x_0_yy = buffer_2000_spsd[3];

    auto g_xx_0_0_0_0_x_0_yz = buffer_2000_spsd[4];

    auto g_xx_0_0_0_0_x_0_zz = buffer_2000_spsd[5];

    auto g_xx_0_0_0_0_y_0_xx = buffer_2000_spsd[6];

    auto g_xx_0_0_0_0_y_0_xy = buffer_2000_spsd[7];

    auto g_xx_0_0_0_0_y_0_xz = buffer_2000_spsd[8];

    auto g_xx_0_0_0_0_y_0_yy = buffer_2000_spsd[9];

    auto g_xx_0_0_0_0_y_0_yz = buffer_2000_spsd[10];

    auto g_xx_0_0_0_0_y_0_zz = buffer_2000_spsd[11];

    auto g_xx_0_0_0_0_z_0_xx = buffer_2000_spsd[12];

    auto g_xx_0_0_0_0_z_0_xy = buffer_2000_spsd[13];

    auto g_xx_0_0_0_0_z_0_xz = buffer_2000_spsd[14];

    auto g_xx_0_0_0_0_z_0_yy = buffer_2000_spsd[15];

    auto g_xx_0_0_0_0_z_0_yz = buffer_2000_spsd[16];

    auto g_xx_0_0_0_0_z_0_zz = buffer_2000_spsd[17];

    auto g_xy_0_0_0_0_x_0_xx = buffer_2000_spsd[18];

    auto g_xy_0_0_0_0_x_0_xy = buffer_2000_spsd[19];

    auto g_xy_0_0_0_0_x_0_xz = buffer_2000_spsd[20];

    auto g_xy_0_0_0_0_x_0_yy = buffer_2000_spsd[21];

    auto g_xy_0_0_0_0_x_0_yz = buffer_2000_spsd[22];

    auto g_xy_0_0_0_0_x_0_zz = buffer_2000_spsd[23];

    auto g_xy_0_0_0_0_y_0_xx = buffer_2000_spsd[24];

    auto g_xy_0_0_0_0_y_0_xy = buffer_2000_spsd[25];

    auto g_xy_0_0_0_0_y_0_xz = buffer_2000_spsd[26];

    auto g_xy_0_0_0_0_y_0_yy = buffer_2000_spsd[27];

    auto g_xy_0_0_0_0_y_0_yz = buffer_2000_spsd[28];

    auto g_xy_0_0_0_0_y_0_zz = buffer_2000_spsd[29];

    auto g_xy_0_0_0_0_z_0_xx = buffer_2000_spsd[30];

    auto g_xy_0_0_0_0_z_0_xy = buffer_2000_spsd[31];

    auto g_xy_0_0_0_0_z_0_xz = buffer_2000_spsd[32];

    auto g_xy_0_0_0_0_z_0_yy = buffer_2000_spsd[33];

    auto g_xy_0_0_0_0_z_0_yz = buffer_2000_spsd[34];

    auto g_xy_0_0_0_0_z_0_zz = buffer_2000_spsd[35];

    auto g_xz_0_0_0_0_x_0_xx = buffer_2000_spsd[36];

    auto g_xz_0_0_0_0_x_0_xy = buffer_2000_spsd[37];

    auto g_xz_0_0_0_0_x_0_xz = buffer_2000_spsd[38];

    auto g_xz_0_0_0_0_x_0_yy = buffer_2000_spsd[39];

    auto g_xz_0_0_0_0_x_0_yz = buffer_2000_spsd[40];

    auto g_xz_0_0_0_0_x_0_zz = buffer_2000_spsd[41];

    auto g_xz_0_0_0_0_y_0_xx = buffer_2000_spsd[42];

    auto g_xz_0_0_0_0_y_0_xy = buffer_2000_spsd[43];

    auto g_xz_0_0_0_0_y_0_xz = buffer_2000_spsd[44];

    auto g_xz_0_0_0_0_y_0_yy = buffer_2000_spsd[45];

    auto g_xz_0_0_0_0_y_0_yz = buffer_2000_spsd[46];

    auto g_xz_0_0_0_0_y_0_zz = buffer_2000_spsd[47];

    auto g_xz_0_0_0_0_z_0_xx = buffer_2000_spsd[48];

    auto g_xz_0_0_0_0_z_0_xy = buffer_2000_spsd[49];

    auto g_xz_0_0_0_0_z_0_xz = buffer_2000_spsd[50];

    auto g_xz_0_0_0_0_z_0_yy = buffer_2000_spsd[51];

    auto g_xz_0_0_0_0_z_0_yz = buffer_2000_spsd[52];

    auto g_xz_0_0_0_0_z_0_zz = buffer_2000_spsd[53];

    auto g_yy_0_0_0_0_x_0_xx = buffer_2000_spsd[54];

    auto g_yy_0_0_0_0_x_0_xy = buffer_2000_spsd[55];

    auto g_yy_0_0_0_0_x_0_xz = buffer_2000_spsd[56];

    auto g_yy_0_0_0_0_x_0_yy = buffer_2000_spsd[57];

    auto g_yy_0_0_0_0_x_0_yz = buffer_2000_spsd[58];

    auto g_yy_0_0_0_0_x_0_zz = buffer_2000_spsd[59];

    auto g_yy_0_0_0_0_y_0_xx = buffer_2000_spsd[60];

    auto g_yy_0_0_0_0_y_0_xy = buffer_2000_spsd[61];

    auto g_yy_0_0_0_0_y_0_xz = buffer_2000_spsd[62];

    auto g_yy_0_0_0_0_y_0_yy = buffer_2000_spsd[63];

    auto g_yy_0_0_0_0_y_0_yz = buffer_2000_spsd[64];

    auto g_yy_0_0_0_0_y_0_zz = buffer_2000_spsd[65];

    auto g_yy_0_0_0_0_z_0_xx = buffer_2000_spsd[66];

    auto g_yy_0_0_0_0_z_0_xy = buffer_2000_spsd[67];

    auto g_yy_0_0_0_0_z_0_xz = buffer_2000_spsd[68];

    auto g_yy_0_0_0_0_z_0_yy = buffer_2000_spsd[69];

    auto g_yy_0_0_0_0_z_0_yz = buffer_2000_spsd[70];

    auto g_yy_0_0_0_0_z_0_zz = buffer_2000_spsd[71];

    auto g_yz_0_0_0_0_x_0_xx = buffer_2000_spsd[72];

    auto g_yz_0_0_0_0_x_0_xy = buffer_2000_spsd[73];

    auto g_yz_0_0_0_0_x_0_xz = buffer_2000_spsd[74];

    auto g_yz_0_0_0_0_x_0_yy = buffer_2000_spsd[75];

    auto g_yz_0_0_0_0_x_0_yz = buffer_2000_spsd[76];

    auto g_yz_0_0_0_0_x_0_zz = buffer_2000_spsd[77];

    auto g_yz_0_0_0_0_y_0_xx = buffer_2000_spsd[78];

    auto g_yz_0_0_0_0_y_0_xy = buffer_2000_spsd[79];

    auto g_yz_0_0_0_0_y_0_xz = buffer_2000_spsd[80];

    auto g_yz_0_0_0_0_y_0_yy = buffer_2000_spsd[81];

    auto g_yz_0_0_0_0_y_0_yz = buffer_2000_spsd[82];

    auto g_yz_0_0_0_0_y_0_zz = buffer_2000_spsd[83];

    auto g_yz_0_0_0_0_z_0_xx = buffer_2000_spsd[84];

    auto g_yz_0_0_0_0_z_0_xy = buffer_2000_spsd[85];

    auto g_yz_0_0_0_0_z_0_xz = buffer_2000_spsd[86];

    auto g_yz_0_0_0_0_z_0_yy = buffer_2000_spsd[87];

    auto g_yz_0_0_0_0_z_0_yz = buffer_2000_spsd[88];

    auto g_yz_0_0_0_0_z_0_zz = buffer_2000_spsd[89];

    auto g_zz_0_0_0_0_x_0_xx = buffer_2000_spsd[90];

    auto g_zz_0_0_0_0_x_0_xy = buffer_2000_spsd[91];

    auto g_zz_0_0_0_0_x_0_xz = buffer_2000_spsd[92];

    auto g_zz_0_0_0_0_x_0_yy = buffer_2000_spsd[93];

    auto g_zz_0_0_0_0_x_0_yz = buffer_2000_spsd[94];

    auto g_zz_0_0_0_0_x_0_zz = buffer_2000_spsd[95];

    auto g_zz_0_0_0_0_y_0_xx = buffer_2000_spsd[96];

    auto g_zz_0_0_0_0_y_0_xy = buffer_2000_spsd[97];

    auto g_zz_0_0_0_0_y_0_xz = buffer_2000_spsd[98];

    auto g_zz_0_0_0_0_y_0_yy = buffer_2000_spsd[99];

    auto g_zz_0_0_0_0_y_0_yz = buffer_2000_spsd[100];

    auto g_zz_0_0_0_0_y_0_zz = buffer_2000_spsd[101];

    auto g_zz_0_0_0_0_z_0_xx = buffer_2000_spsd[102];

    auto g_zz_0_0_0_0_z_0_xy = buffer_2000_spsd[103];

    auto g_zz_0_0_0_0_z_0_xz = buffer_2000_spsd[104];

    auto g_zz_0_0_0_0_z_0_yy = buffer_2000_spsd[105];

    auto g_zz_0_0_0_0_z_0_yz = buffer_2000_spsd[106];

    auto g_zz_0_0_0_0_z_0_zz = buffer_2000_spsd[107];

    // integrals block (0-6)

    #pragma omp simd aligned(g_0_x_0_xx, g_0_x_0_xy, g_0_x_0_xz, g_0_x_0_yy, g_0_x_0_yz, g_0_x_0_zz, g_xx_0_0_0_0_x_0_xx, g_xx_0_0_0_0_x_0_xy, g_xx_0_0_0_0_x_0_xz, g_xx_0_0_0_0_x_0_yy, g_xx_0_0_0_0_x_0_yz, g_xx_0_0_0_0_x_0_zz, g_xx_x_0_xx, g_xx_x_0_xy, g_xx_x_0_xz, g_xx_x_0_yy, g_xx_x_0_yz, g_xx_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_x_0_xx[i] = -2.0 * g_0_x_0_xx[i] * a_exp + 4.0 * g_xx_x_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_0_xy[i] = -2.0 * g_0_x_0_xy[i] * a_exp + 4.0 * g_xx_x_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_0_xz[i] = -2.0 * g_0_x_0_xz[i] * a_exp + 4.0 * g_xx_x_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_0_yy[i] = -2.0 * g_0_x_0_yy[i] * a_exp + 4.0 * g_xx_x_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_0_yz[i] = -2.0 * g_0_x_0_yz[i] * a_exp + 4.0 * g_xx_x_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_0_zz[i] = -2.0 * g_0_x_0_zz[i] * a_exp + 4.0 * g_xx_x_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_0_y_0_xx, g_0_y_0_xy, g_0_y_0_xz, g_0_y_0_yy, g_0_y_0_yz, g_0_y_0_zz, g_xx_0_0_0_0_y_0_xx, g_xx_0_0_0_0_y_0_xy, g_xx_0_0_0_0_y_0_xz, g_xx_0_0_0_0_y_0_yy, g_xx_0_0_0_0_y_0_yz, g_xx_0_0_0_0_y_0_zz, g_xx_y_0_xx, g_xx_y_0_xy, g_xx_y_0_xz, g_xx_y_0_yy, g_xx_y_0_yz, g_xx_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_y_0_xx[i] = -2.0 * g_0_y_0_xx[i] * a_exp + 4.0 * g_xx_y_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_0_xy[i] = -2.0 * g_0_y_0_xy[i] * a_exp + 4.0 * g_xx_y_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_0_xz[i] = -2.0 * g_0_y_0_xz[i] * a_exp + 4.0 * g_xx_y_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_0_yy[i] = -2.0 * g_0_y_0_yy[i] * a_exp + 4.0 * g_xx_y_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_0_yz[i] = -2.0 * g_0_y_0_yz[i] * a_exp + 4.0 * g_xx_y_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_0_zz[i] = -2.0 * g_0_y_0_zz[i] * a_exp + 4.0 * g_xx_y_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_0_z_0_xx, g_0_z_0_xy, g_0_z_0_xz, g_0_z_0_yy, g_0_z_0_yz, g_0_z_0_zz, g_xx_0_0_0_0_z_0_xx, g_xx_0_0_0_0_z_0_xy, g_xx_0_0_0_0_z_0_xz, g_xx_0_0_0_0_z_0_yy, g_xx_0_0_0_0_z_0_yz, g_xx_0_0_0_0_z_0_zz, g_xx_z_0_xx, g_xx_z_0_xy, g_xx_z_0_xz, g_xx_z_0_yy, g_xx_z_0_yz, g_xx_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_z_0_xx[i] = -2.0 * g_0_z_0_xx[i] * a_exp + 4.0 * g_xx_z_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_0_xy[i] = -2.0 * g_0_z_0_xy[i] * a_exp + 4.0 * g_xx_z_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_0_xz[i] = -2.0 * g_0_z_0_xz[i] * a_exp + 4.0 * g_xx_z_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_0_yy[i] = -2.0 * g_0_z_0_yy[i] * a_exp + 4.0 * g_xx_z_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_0_yz[i] = -2.0 * g_0_z_0_yz[i] * a_exp + 4.0 * g_xx_z_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_0_zz[i] = -2.0 * g_0_z_0_zz[i] * a_exp + 4.0 * g_xx_z_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_xy_0_0_0_0_x_0_xx, g_xy_0_0_0_0_x_0_xy, g_xy_0_0_0_0_x_0_xz, g_xy_0_0_0_0_x_0_yy, g_xy_0_0_0_0_x_0_yz, g_xy_0_0_0_0_x_0_zz, g_xy_x_0_xx, g_xy_x_0_xy, g_xy_x_0_xz, g_xy_x_0_yy, g_xy_x_0_yz, g_xy_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_x_0_xx[i] = 4.0 * g_xy_x_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_0_xy[i] = 4.0 * g_xy_x_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_0_xz[i] = 4.0 * g_xy_x_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_0_yy[i] = 4.0 * g_xy_x_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_0_yz[i] = 4.0 * g_xy_x_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_0_zz[i] = 4.0 * g_xy_x_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_xy_0_0_0_0_y_0_xx, g_xy_0_0_0_0_y_0_xy, g_xy_0_0_0_0_y_0_xz, g_xy_0_0_0_0_y_0_yy, g_xy_0_0_0_0_y_0_yz, g_xy_0_0_0_0_y_0_zz, g_xy_y_0_xx, g_xy_y_0_xy, g_xy_y_0_xz, g_xy_y_0_yy, g_xy_y_0_yz, g_xy_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_y_0_xx[i] = 4.0 * g_xy_y_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_0_xy[i] = 4.0 * g_xy_y_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_0_xz[i] = 4.0 * g_xy_y_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_0_yy[i] = 4.0 * g_xy_y_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_0_yz[i] = 4.0 * g_xy_y_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_0_zz[i] = 4.0 * g_xy_y_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_xy_0_0_0_0_z_0_xx, g_xy_0_0_0_0_z_0_xy, g_xy_0_0_0_0_z_0_xz, g_xy_0_0_0_0_z_0_yy, g_xy_0_0_0_0_z_0_yz, g_xy_0_0_0_0_z_0_zz, g_xy_z_0_xx, g_xy_z_0_xy, g_xy_z_0_xz, g_xy_z_0_yy, g_xy_z_0_yz, g_xy_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_z_0_xx[i] = 4.0 * g_xy_z_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_0_xy[i] = 4.0 * g_xy_z_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_0_xz[i] = 4.0 * g_xy_z_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_0_yy[i] = 4.0 * g_xy_z_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_0_yz[i] = 4.0 * g_xy_z_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_0_zz[i] = 4.0 * g_xy_z_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_xz_0_0_0_0_x_0_xx, g_xz_0_0_0_0_x_0_xy, g_xz_0_0_0_0_x_0_xz, g_xz_0_0_0_0_x_0_yy, g_xz_0_0_0_0_x_0_yz, g_xz_0_0_0_0_x_0_zz, g_xz_x_0_xx, g_xz_x_0_xy, g_xz_x_0_xz, g_xz_x_0_yy, g_xz_x_0_yz, g_xz_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_x_0_xx[i] = 4.0 * g_xz_x_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_0_xy[i] = 4.0 * g_xz_x_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_0_xz[i] = 4.0 * g_xz_x_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_0_yy[i] = 4.0 * g_xz_x_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_0_yz[i] = 4.0 * g_xz_x_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_0_zz[i] = 4.0 * g_xz_x_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_xz_0_0_0_0_y_0_xx, g_xz_0_0_0_0_y_0_xy, g_xz_0_0_0_0_y_0_xz, g_xz_0_0_0_0_y_0_yy, g_xz_0_0_0_0_y_0_yz, g_xz_0_0_0_0_y_0_zz, g_xz_y_0_xx, g_xz_y_0_xy, g_xz_y_0_xz, g_xz_y_0_yy, g_xz_y_0_yz, g_xz_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_y_0_xx[i] = 4.0 * g_xz_y_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_0_xy[i] = 4.0 * g_xz_y_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_0_xz[i] = 4.0 * g_xz_y_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_0_yy[i] = 4.0 * g_xz_y_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_0_yz[i] = 4.0 * g_xz_y_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_0_zz[i] = 4.0 * g_xz_y_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_xz_0_0_0_0_z_0_xx, g_xz_0_0_0_0_z_0_xy, g_xz_0_0_0_0_z_0_xz, g_xz_0_0_0_0_z_0_yy, g_xz_0_0_0_0_z_0_yz, g_xz_0_0_0_0_z_0_zz, g_xz_z_0_xx, g_xz_z_0_xy, g_xz_z_0_xz, g_xz_z_0_yy, g_xz_z_0_yz, g_xz_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_z_0_xx[i] = 4.0 * g_xz_z_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_0_xy[i] = 4.0 * g_xz_z_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_0_xz[i] = 4.0 * g_xz_z_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_0_yy[i] = 4.0 * g_xz_z_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_0_yz[i] = 4.0 * g_xz_z_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_0_zz[i] = 4.0 * g_xz_z_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_0_x_0_xx, g_0_x_0_xy, g_0_x_0_xz, g_0_x_0_yy, g_0_x_0_yz, g_0_x_0_zz, g_yy_0_0_0_0_x_0_xx, g_yy_0_0_0_0_x_0_xy, g_yy_0_0_0_0_x_0_xz, g_yy_0_0_0_0_x_0_yy, g_yy_0_0_0_0_x_0_yz, g_yy_0_0_0_0_x_0_zz, g_yy_x_0_xx, g_yy_x_0_xy, g_yy_x_0_xz, g_yy_x_0_yy, g_yy_x_0_yz, g_yy_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_x_0_xx[i] = -2.0 * g_0_x_0_xx[i] * a_exp + 4.0 * g_yy_x_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_0_xy[i] = -2.0 * g_0_x_0_xy[i] * a_exp + 4.0 * g_yy_x_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_0_xz[i] = -2.0 * g_0_x_0_xz[i] * a_exp + 4.0 * g_yy_x_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_0_yy[i] = -2.0 * g_0_x_0_yy[i] * a_exp + 4.0 * g_yy_x_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_0_yz[i] = -2.0 * g_0_x_0_yz[i] * a_exp + 4.0 * g_yy_x_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_0_zz[i] = -2.0 * g_0_x_0_zz[i] * a_exp + 4.0 * g_yy_x_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_0_y_0_xx, g_0_y_0_xy, g_0_y_0_xz, g_0_y_0_yy, g_0_y_0_yz, g_0_y_0_zz, g_yy_0_0_0_0_y_0_xx, g_yy_0_0_0_0_y_0_xy, g_yy_0_0_0_0_y_0_xz, g_yy_0_0_0_0_y_0_yy, g_yy_0_0_0_0_y_0_yz, g_yy_0_0_0_0_y_0_zz, g_yy_y_0_xx, g_yy_y_0_xy, g_yy_y_0_xz, g_yy_y_0_yy, g_yy_y_0_yz, g_yy_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_y_0_xx[i] = -2.0 * g_0_y_0_xx[i] * a_exp + 4.0 * g_yy_y_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_0_xy[i] = -2.0 * g_0_y_0_xy[i] * a_exp + 4.0 * g_yy_y_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_0_xz[i] = -2.0 * g_0_y_0_xz[i] * a_exp + 4.0 * g_yy_y_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_0_yy[i] = -2.0 * g_0_y_0_yy[i] * a_exp + 4.0 * g_yy_y_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_0_yz[i] = -2.0 * g_0_y_0_yz[i] * a_exp + 4.0 * g_yy_y_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_0_zz[i] = -2.0 * g_0_y_0_zz[i] * a_exp + 4.0 * g_yy_y_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_0_z_0_xx, g_0_z_0_xy, g_0_z_0_xz, g_0_z_0_yy, g_0_z_0_yz, g_0_z_0_zz, g_yy_0_0_0_0_z_0_xx, g_yy_0_0_0_0_z_0_xy, g_yy_0_0_0_0_z_0_xz, g_yy_0_0_0_0_z_0_yy, g_yy_0_0_0_0_z_0_yz, g_yy_0_0_0_0_z_0_zz, g_yy_z_0_xx, g_yy_z_0_xy, g_yy_z_0_xz, g_yy_z_0_yy, g_yy_z_0_yz, g_yy_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_z_0_xx[i] = -2.0 * g_0_z_0_xx[i] * a_exp + 4.0 * g_yy_z_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_0_xy[i] = -2.0 * g_0_z_0_xy[i] * a_exp + 4.0 * g_yy_z_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_0_xz[i] = -2.0 * g_0_z_0_xz[i] * a_exp + 4.0 * g_yy_z_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_0_yy[i] = -2.0 * g_0_z_0_yy[i] * a_exp + 4.0 * g_yy_z_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_0_yz[i] = -2.0 * g_0_z_0_yz[i] * a_exp + 4.0 * g_yy_z_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_0_zz[i] = -2.0 * g_0_z_0_zz[i] * a_exp + 4.0 * g_yy_z_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_yz_0_0_0_0_x_0_xx, g_yz_0_0_0_0_x_0_xy, g_yz_0_0_0_0_x_0_xz, g_yz_0_0_0_0_x_0_yy, g_yz_0_0_0_0_x_0_yz, g_yz_0_0_0_0_x_0_zz, g_yz_x_0_xx, g_yz_x_0_xy, g_yz_x_0_xz, g_yz_x_0_yy, g_yz_x_0_yz, g_yz_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_x_0_xx[i] = 4.0 * g_yz_x_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_0_xy[i] = 4.0 * g_yz_x_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_0_xz[i] = 4.0 * g_yz_x_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_0_yy[i] = 4.0 * g_yz_x_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_0_yz[i] = 4.0 * g_yz_x_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_0_zz[i] = 4.0 * g_yz_x_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_yz_0_0_0_0_y_0_xx, g_yz_0_0_0_0_y_0_xy, g_yz_0_0_0_0_y_0_xz, g_yz_0_0_0_0_y_0_yy, g_yz_0_0_0_0_y_0_yz, g_yz_0_0_0_0_y_0_zz, g_yz_y_0_xx, g_yz_y_0_xy, g_yz_y_0_xz, g_yz_y_0_yy, g_yz_y_0_yz, g_yz_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_y_0_xx[i] = 4.0 * g_yz_y_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_0_xy[i] = 4.0 * g_yz_y_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_0_xz[i] = 4.0 * g_yz_y_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_0_yy[i] = 4.0 * g_yz_y_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_0_yz[i] = 4.0 * g_yz_y_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_0_zz[i] = 4.0 * g_yz_y_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_yz_0_0_0_0_z_0_xx, g_yz_0_0_0_0_z_0_xy, g_yz_0_0_0_0_z_0_xz, g_yz_0_0_0_0_z_0_yy, g_yz_0_0_0_0_z_0_yz, g_yz_0_0_0_0_z_0_zz, g_yz_z_0_xx, g_yz_z_0_xy, g_yz_z_0_xz, g_yz_z_0_yy, g_yz_z_0_yz, g_yz_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_z_0_xx[i] = 4.0 * g_yz_z_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_0_xy[i] = 4.0 * g_yz_z_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_0_xz[i] = 4.0 * g_yz_z_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_0_yy[i] = 4.0 * g_yz_z_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_0_yz[i] = 4.0 * g_yz_z_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_0_zz[i] = 4.0 * g_yz_z_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_0_x_0_xx, g_0_x_0_xy, g_0_x_0_xz, g_0_x_0_yy, g_0_x_0_yz, g_0_x_0_zz, g_zz_0_0_0_0_x_0_xx, g_zz_0_0_0_0_x_0_xy, g_zz_0_0_0_0_x_0_xz, g_zz_0_0_0_0_x_0_yy, g_zz_0_0_0_0_x_0_yz, g_zz_0_0_0_0_x_0_zz, g_zz_x_0_xx, g_zz_x_0_xy, g_zz_x_0_xz, g_zz_x_0_yy, g_zz_x_0_yz, g_zz_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_x_0_xx[i] = -2.0 * g_0_x_0_xx[i] * a_exp + 4.0 * g_zz_x_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_0_xy[i] = -2.0 * g_0_x_0_xy[i] * a_exp + 4.0 * g_zz_x_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_0_xz[i] = -2.0 * g_0_x_0_xz[i] * a_exp + 4.0 * g_zz_x_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_0_yy[i] = -2.0 * g_0_x_0_yy[i] * a_exp + 4.0 * g_zz_x_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_0_yz[i] = -2.0 * g_0_x_0_yz[i] * a_exp + 4.0 * g_zz_x_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_0_zz[i] = -2.0 * g_0_x_0_zz[i] * a_exp + 4.0 * g_zz_x_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_0_y_0_xx, g_0_y_0_xy, g_0_y_0_xz, g_0_y_0_yy, g_0_y_0_yz, g_0_y_0_zz, g_zz_0_0_0_0_y_0_xx, g_zz_0_0_0_0_y_0_xy, g_zz_0_0_0_0_y_0_xz, g_zz_0_0_0_0_y_0_yy, g_zz_0_0_0_0_y_0_yz, g_zz_0_0_0_0_y_0_zz, g_zz_y_0_xx, g_zz_y_0_xy, g_zz_y_0_xz, g_zz_y_0_yy, g_zz_y_0_yz, g_zz_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_y_0_xx[i] = -2.0 * g_0_y_0_xx[i] * a_exp + 4.0 * g_zz_y_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_0_xy[i] = -2.0 * g_0_y_0_xy[i] * a_exp + 4.0 * g_zz_y_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_0_xz[i] = -2.0 * g_0_y_0_xz[i] * a_exp + 4.0 * g_zz_y_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_0_yy[i] = -2.0 * g_0_y_0_yy[i] * a_exp + 4.0 * g_zz_y_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_0_yz[i] = -2.0 * g_0_y_0_yz[i] * a_exp + 4.0 * g_zz_y_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_0_zz[i] = -2.0 * g_0_y_0_zz[i] * a_exp + 4.0 * g_zz_y_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_0_z_0_xx, g_0_z_0_xy, g_0_z_0_xz, g_0_z_0_yy, g_0_z_0_yz, g_0_z_0_zz, g_zz_0_0_0_0_z_0_xx, g_zz_0_0_0_0_z_0_xy, g_zz_0_0_0_0_z_0_xz, g_zz_0_0_0_0_z_0_yy, g_zz_0_0_0_0_z_0_yz, g_zz_0_0_0_0_z_0_zz, g_zz_z_0_xx, g_zz_z_0_xy, g_zz_z_0_xz, g_zz_z_0_yy, g_zz_z_0_yz, g_zz_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_z_0_xx[i] = -2.0 * g_0_z_0_xx[i] * a_exp + 4.0 * g_zz_z_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_0_xy[i] = -2.0 * g_0_z_0_xy[i] * a_exp + 4.0 * g_zz_z_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_0_xz[i] = -2.0 * g_0_z_0_xz[i] * a_exp + 4.0 * g_zz_z_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_0_yy[i] = -2.0 * g_0_z_0_yy[i] * a_exp + 4.0 * g_zz_z_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_0_yz[i] = -2.0 * g_0_z_0_yz[i] * a_exp + 4.0 * g_zz_z_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_0_zz[i] = -2.0 * g_0_z_0_zz[i] * a_exp + 4.0 * g_zz_z_0_zz[i] * a_exp * a_exp;
    }
}

} // t4c_geom namespace


#include "GeomDeriv1000OfScalarForSSDD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1000_ssdd_0(CSimdArray<double>& buffer_1000_ssdd,
                     const CSimdArray<double>& buffer_psdd,
                     const double a_exp) -> void
{
    const auto ndims = buffer_1000_ssdd.number_of_columns();

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

    /// Set up components of integrals buffer : buffer_1000_ssdd

    auto g_x_0_0_0_0_0_xx_xx = buffer_1000_ssdd[0];

    auto g_x_0_0_0_0_0_xx_xy = buffer_1000_ssdd[1];

    auto g_x_0_0_0_0_0_xx_xz = buffer_1000_ssdd[2];

    auto g_x_0_0_0_0_0_xx_yy = buffer_1000_ssdd[3];

    auto g_x_0_0_0_0_0_xx_yz = buffer_1000_ssdd[4];

    auto g_x_0_0_0_0_0_xx_zz = buffer_1000_ssdd[5];

    auto g_x_0_0_0_0_0_xy_xx = buffer_1000_ssdd[6];

    auto g_x_0_0_0_0_0_xy_xy = buffer_1000_ssdd[7];

    auto g_x_0_0_0_0_0_xy_xz = buffer_1000_ssdd[8];

    auto g_x_0_0_0_0_0_xy_yy = buffer_1000_ssdd[9];

    auto g_x_0_0_0_0_0_xy_yz = buffer_1000_ssdd[10];

    auto g_x_0_0_0_0_0_xy_zz = buffer_1000_ssdd[11];

    auto g_x_0_0_0_0_0_xz_xx = buffer_1000_ssdd[12];

    auto g_x_0_0_0_0_0_xz_xy = buffer_1000_ssdd[13];

    auto g_x_0_0_0_0_0_xz_xz = buffer_1000_ssdd[14];

    auto g_x_0_0_0_0_0_xz_yy = buffer_1000_ssdd[15];

    auto g_x_0_0_0_0_0_xz_yz = buffer_1000_ssdd[16];

    auto g_x_0_0_0_0_0_xz_zz = buffer_1000_ssdd[17];

    auto g_x_0_0_0_0_0_yy_xx = buffer_1000_ssdd[18];

    auto g_x_0_0_0_0_0_yy_xy = buffer_1000_ssdd[19];

    auto g_x_0_0_0_0_0_yy_xz = buffer_1000_ssdd[20];

    auto g_x_0_0_0_0_0_yy_yy = buffer_1000_ssdd[21];

    auto g_x_0_0_0_0_0_yy_yz = buffer_1000_ssdd[22];

    auto g_x_0_0_0_0_0_yy_zz = buffer_1000_ssdd[23];

    auto g_x_0_0_0_0_0_yz_xx = buffer_1000_ssdd[24];

    auto g_x_0_0_0_0_0_yz_xy = buffer_1000_ssdd[25];

    auto g_x_0_0_0_0_0_yz_xz = buffer_1000_ssdd[26];

    auto g_x_0_0_0_0_0_yz_yy = buffer_1000_ssdd[27];

    auto g_x_0_0_0_0_0_yz_yz = buffer_1000_ssdd[28];

    auto g_x_0_0_0_0_0_yz_zz = buffer_1000_ssdd[29];

    auto g_x_0_0_0_0_0_zz_xx = buffer_1000_ssdd[30];

    auto g_x_0_0_0_0_0_zz_xy = buffer_1000_ssdd[31];

    auto g_x_0_0_0_0_0_zz_xz = buffer_1000_ssdd[32];

    auto g_x_0_0_0_0_0_zz_yy = buffer_1000_ssdd[33];

    auto g_x_0_0_0_0_0_zz_yz = buffer_1000_ssdd[34];

    auto g_x_0_0_0_0_0_zz_zz = buffer_1000_ssdd[35];

    auto g_y_0_0_0_0_0_xx_xx = buffer_1000_ssdd[36];

    auto g_y_0_0_0_0_0_xx_xy = buffer_1000_ssdd[37];

    auto g_y_0_0_0_0_0_xx_xz = buffer_1000_ssdd[38];

    auto g_y_0_0_0_0_0_xx_yy = buffer_1000_ssdd[39];

    auto g_y_0_0_0_0_0_xx_yz = buffer_1000_ssdd[40];

    auto g_y_0_0_0_0_0_xx_zz = buffer_1000_ssdd[41];

    auto g_y_0_0_0_0_0_xy_xx = buffer_1000_ssdd[42];

    auto g_y_0_0_0_0_0_xy_xy = buffer_1000_ssdd[43];

    auto g_y_0_0_0_0_0_xy_xz = buffer_1000_ssdd[44];

    auto g_y_0_0_0_0_0_xy_yy = buffer_1000_ssdd[45];

    auto g_y_0_0_0_0_0_xy_yz = buffer_1000_ssdd[46];

    auto g_y_0_0_0_0_0_xy_zz = buffer_1000_ssdd[47];

    auto g_y_0_0_0_0_0_xz_xx = buffer_1000_ssdd[48];

    auto g_y_0_0_0_0_0_xz_xy = buffer_1000_ssdd[49];

    auto g_y_0_0_0_0_0_xz_xz = buffer_1000_ssdd[50];

    auto g_y_0_0_0_0_0_xz_yy = buffer_1000_ssdd[51];

    auto g_y_0_0_0_0_0_xz_yz = buffer_1000_ssdd[52];

    auto g_y_0_0_0_0_0_xz_zz = buffer_1000_ssdd[53];

    auto g_y_0_0_0_0_0_yy_xx = buffer_1000_ssdd[54];

    auto g_y_0_0_0_0_0_yy_xy = buffer_1000_ssdd[55];

    auto g_y_0_0_0_0_0_yy_xz = buffer_1000_ssdd[56];

    auto g_y_0_0_0_0_0_yy_yy = buffer_1000_ssdd[57];

    auto g_y_0_0_0_0_0_yy_yz = buffer_1000_ssdd[58];

    auto g_y_0_0_0_0_0_yy_zz = buffer_1000_ssdd[59];

    auto g_y_0_0_0_0_0_yz_xx = buffer_1000_ssdd[60];

    auto g_y_0_0_0_0_0_yz_xy = buffer_1000_ssdd[61];

    auto g_y_0_0_0_0_0_yz_xz = buffer_1000_ssdd[62];

    auto g_y_0_0_0_0_0_yz_yy = buffer_1000_ssdd[63];

    auto g_y_0_0_0_0_0_yz_yz = buffer_1000_ssdd[64];

    auto g_y_0_0_0_0_0_yz_zz = buffer_1000_ssdd[65];

    auto g_y_0_0_0_0_0_zz_xx = buffer_1000_ssdd[66];

    auto g_y_0_0_0_0_0_zz_xy = buffer_1000_ssdd[67];

    auto g_y_0_0_0_0_0_zz_xz = buffer_1000_ssdd[68];

    auto g_y_0_0_0_0_0_zz_yy = buffer_1000_ssdd[69];

    auto g_y_0_0_0_0_0_zz_yz = buffer_1000_ssdd[70];

    auto g_y_0_0_0_0_0_zz_zz = buffer_1000_ssdd[71];

    auto g_z_0_0_0_0_0_xx_xx = buffer_1000_ssdd[72];

    auto g_z_0_0_0_0_0_xx_xy = buffer_1000_ssdd[73];

    auto g_z_0_0_0_0_0_xx_xz = buffer_1000_ssdd[74];

    auto g_z_0_0_0_0_0_xx_yy = buffer_1000_ssdd[75];

    auto g_z_0_0_0_0_0_xx_yz = buffer_1000_ssdd[76];

    auto g_z_0_0_0_0_0_xx_zz = buffer_1000_ssdd[77];

    auto g_z_0_0_0_0_0_xy_xx = buffer_1000_ssdd[78];

    auto g_z_0_0_0_0_0_xy_xy = buffer_1000_ssdd[79];

    auto g_z_0_0_0_0_0_xy_xz = buffer_1000_ssdd[80];

    auto g_z_0_0_0_0_0_xy_yy = buffer_1000_ssdd[81];

    auto g_z_0_0_0_0_0_xy_yz = buffer_1000_ssdd[82];

    auto g_z_0_0_0_0_0_xy_zz = buffer_1000_ssdd[83];

    auto g_z_0_0_0_0_0_xz_xx = buffer_1000_ssdd[84];

    auto g_z_0_0_0_0_0_xz_xy = buffer_1000_ssdd[85];

    auto g_z_0_0_0_0_0_xz_xz = buffer_1000_ssdd[86];

    auto g_z_0_0_0_0_0_xz_yy = buffer_1000_ssdd[87];

    auto g_z_0_0_0_0_0_xz_yz = buffer_1000_ssdd[88];

    auto g_z_0_0_0_0_0_xz_zz = buffer_1000_ssdd[89];

    auto g_z_0_0_0_0_0_yy_xx = buffer_1000_ssdd[90];

    auto g_z_0_0_0_0_0_yy_xy = buffer_1000_ssdd[91];

    auto g_z_0_0_0_0_0_yy_xz = buffer_1000_ssdd[92];

    auto g_z_0_0_0_0_0_yy_yy = buffer_1000_ssdd[93];

    auto g_z_0_0_0_0_0_yy_yz = buffer_1000_ssdd[94];

    auto g_z_0_0_0_0_0_yy_zz = buffer_1000_ssdd[95];

    auto g_z_0_0_0_0_0_yz_xx = buffer_1000_ssdd[96];

    auto g_z_0_0_0_0_0_yz_xy = buffer_1000_ssdd[97];

    auto g_z_0_0_0_0_0_yz_xz = buffer_1000_ssdd[98];

    auto g_z_0_0_0_0_0_yz_yy = buffer_1000_ssdd[99];

    auto g_z_0_0_0_0_0_yz_yz = buffer_1000_ssdd[100];

    auto g_z_0_0_0_0_0_yz_zz = buffer_1000_ssdd[101];

    auto g_z_0_0_0_0_0_zz_xx = buffer_1000_ssdd[102];

    auto g_z_0_0_0_0_0_zz_xy = buffer_1000_ssdd[103];

    auto g_z_0_0_0_0_0_zz_xz = buffer_1000_ssdd[104];

    auto g_z_0_0_0_0_0_zz_yy = buffer_1000_ssdd[105];

    auto g_z_0_0_0_0_0_zz_yz = buffer_1000_ssdd[106];

    auto g_z_0_0_0_0_0_zz_zz = buffer_1000_ssdd[107];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_0_0_0_0_0_xx_xx, g_x_0_0_0_0_0_xx_xy, g_x_0_0_0_0_0_xx_xz, g_x_0_0_0_0_0_xx_yy, g_x_0_0_0_0_0_xx_yz, g_x_0_0_0_0_0_xx_zz, g_x_0_xx_xx, g_x_0_xx_xy, g_x_0_xx_xz, g_x_0_xx_yy, g_x_0_xx_yz, g_x_0_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_0_xx_xx[i] = 2.0 * g_x_0_xx_xx[i] * a_exp;

        g_x_0_0_0_0_0_xx_xy[i] = 2.0 * g_x_0_xx_xy[i] * a_exp;

        g_x_0_0_0_0_0_xx_xz[i] = 2.0 * g_x_0_xx_xz[i] * a_exp;

        g_x_0_0_0_0_0_xx_yy[i] = 2.0 * g_x_0_xx_yy[i] * a_exp;

        g_x_0_0_0_0_0_xx_yz[i] = 2.0 * g_x_0_xx_yz[i] * a_exp;

        g_x_0_0_0_0_0_xx_zz[i] = 2.0 * g_x_0_xx_zz[i] * a_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_0_0_0_0_0_xy_xx, g_x_0_0_0_0_0_xy_xy, g_x_0_0_0_0_0_xy_xz, g_x_0_0_0_0_0_xy_yy, g_x_0_0_0_0_0_xy_yz, g_x_0_0_0_0_0_xy_zz, g_x_0_xy_xx, g_x_0_xy_xy, g_x_0_xy_xz, g_x_0_xy_yy, g_x_0_xy_yz, g_x_0_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_0_xy_xx[i] = 2.0 * g_x_0_xy_xx[i] * a_exp;

        g_x_0_0_0_0_0_xy_xy[i] = 2.0 * g_x_0_xy_xy[i] * a_exp;

        g_x_0_0_0_0_0_xy_xz[i] = 2.0 * g_x_0_xy_xz[i] * a_exp;

        g_x_0_0_0_0_0_xy_yy[i] = 2.0 * g_x_0_xy_yy[i] * a_exp;

        g_x_0_0_0_0_0_xy_yz[i] = 2.0 * g_x_0_xy_yz[i] * a_exp;

        g_x_0_0_0_0_0_xy_zz[i] = 2.0 * g_x_0_xy_zz[i] * a_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_0_0_0_0_0_xz_xx, g_x_0_0_0_0_0_xz_xy, g_x_0_0_0_0_0_xz_xz, g_x_0_0_0_0_0_xz_yy, g_x_0_0_0_0_0_xz_yz, g_x_0_0_0_0_0_xz_zz, g_x_0_xz_xx, g_x_0_xz_xy, g_x_0_xz_xz, g_x_0_xz_yy, g_x_0_xz_yz, g_x_0_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_0_xz_xx[i] = 2.0 * g_x_0_xz_xx[i] * a_exp;

        g_x_0_0_0_0_0_xz_xy[i] = 2.0 * g_x_0_xz_xy[i] * a_exp;

        g_x_0_0_0_0_0_xz_xz[i] = 2.0 * g_x_0_xz_xz[i] * a_exp;

        g_x_0_0_0_0_0_xz_yy[i] = 2.0 * g_x_0_xz_yy[i] * a_exp;

        g_x_0_0_0_0_0_xz_yz[i] = 2.0 * g_x_0_xz_yz[i] * a_exp;

        g_x_0_0_0_0_0_xz_zz[i] = 2.0 * g_x_0_xz_zz[i] * a_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_x_0_0_0_0_0_yy_xx, g_x_0_0_0_0_0_yy_xy, g_x_0_0_0_0_0_yy_xz, g_x_0_0_0_0_0_yy_yy, g_x_0_0_0_0_0_yy_yz, g_x_0_0_0_0_0_yy_zz, g_x_0_yy_xx, g_x_0_yy_xy, g_x_0_yy_xz, g_x_0_yy_yy, g_x_0_yy_yz, g_x_0_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_0_yy_xx[i] = 2.0 * g_x_0_yy_xx[i] * a_exp;

        g_x_0_0_0_0_0_yy_xy[i] = 2.0 * g_x_0_yy_xy[i] * a_exp;

        g_x_0_0_0_0_0_yy_xz[i] = 2.0 * g_x_0_yy_xz[i] * a_exp;

        g_x_0_0_0_0_0_yy_yy[i] = 2.0 * g_x_0_yy_yy[i] * a_exp;

        g_x_0_0_0_0_0_yy_yz[i] = 2.0 * g_x_0_yy_yz[i] * a_exp;

        g_x_0_0_0_0_0_yy_zz[i] = 2.0 * g_x_0_yy_zz[i] * a_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_x_0_0_0_0_0_yz_xx, g_x_0_0_0_0_0_yz_xy, g_x_0_0_0_0_0_yz_xz, g_x_0_0_0_0_0_yz_yy, g_x_0_0_0_0_0_yz_yz, g_x_0_0_0_0_0_yz_zz, g_x_0_yz_xx, g_x_0_yz_xy, g_x_0_yz_xz, g_x_0_yz_yy, g_x_0_yz_yz, g_x_0_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_0_yz_xx[i] = 2.0 * g_x_0_yz_xx[i] * a_exp;

        g_x_0_0_0_0_0_yz_xy[i] = 2.0 * g_x_0_yz_xy[i] * a_exp;

        g_x_0_0_0_0_0_yz_xz[i] = 2.0 * g_x_0_yz_xz[i] * a_exp;

        g_x_0_0_0_0_0_yz_yy[i] = 2.0 * g_x_0_yz_yy[i] * a_exp;

        g_x_0_0_0_0_0_yz_yz[i] = 2.0 * g_x_0_yz_yz[i] * a_exp;

        g_x_0_0_0_0_0_yz_zz[i] = 2.0 * g_x_0_yz_zz[i] * a_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_x_0_0_0_0_0_zz_xx, g_x_0_0_0_0_0_zz_xy, g_x_0_0_0_0_0_zz_xz, g_x_0_0_0_0_0_zz_yy, g_x_0_0_0_0_0_zz_yz, g_x_0_0_0_0_0_zz_zz, g_x_0_zz_xx, g_x_0_zz_xy, g_x_0_zz_xz, g_x_0_zz_yy, g_x_0_zz_yz, g_x_0_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_0_zz_xx[i] = 2.0 * g_x_0_zz_xx[i] * a_exp;

        g_x_0_0_0_0_0_zz_xy[i] = 2.0 * g_x_0_zz_xy[i] * a_exp;

        g_x_0_0_0_0_0_zz_xz[i] = 2.0 * g_x_0_zz_xz[i] * a_exp;

        g_x_0_0_0_0_0_zz_yy[i] = 2.0 * g_x_0_zz_yy[i] * a_exp;

        g_x_0_0_0_0_0_zz_yz[i] = 2.0 * g_x_0_zz_yz[i] * a_exp;

        g_x_0_0_0_0_0_zz_zz[i] = 2.0 * g_x_0_zz_zz[i] * a_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_y_0_0_0_0_0_xx_xx, g_y_0_0_0_0_0_xx_xy, g_y_0_0_0_0_0_xx_xz, g_y_0_0_0_0_0_xx_yy, g_y_0_0_0_0_0_xx_yz, g_y_0_0_0_0_0_xx_zz, g_y_0_xx_xx, g_y_0_xx_xy, g_y_0_xx_xz, g_y_0_xx_yy, g_y_0_xx_yz, g_y_0_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_0_xx_xx[i] = 2.0 * g_y_0_xx_xx[i] * a_exp;

        g_y_0_0_0_0_0_xx_xy[i] = 2.0 * g_y_0_xx_xy[i] * a_exp;

        g_y_0_0_0_0_0_xx_xz[i] = 2.0 * g_y_0_xx_xz[i] * a_exp;

        g_y_0_0_0_0_0_xx_yy[i] = 2.0 * g_y_0_xx_yy[i] * a_exp;

        g_y_0_0_0_0_0_xx_yz[i] = 2.0 * g_y_0_xx_yz[i] * a_exp;

        g_y_0_0_0_0_0_xx_zz[i] = 2.0 * g_y_0_xx_zz[i] * a_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_y_0_0_0_0_0_xy_xx, g_y_0_0_0_0_0_xy_xy, g_y_0_0_0_0_0_xy_xz, g_y_0_0_0_0_0_xy_yy, g_y_0_0_0_0_0_xy_yz, g_y_0_0_0_0_0_xy_zz, g_y_0_xy_xx, g_y_0_xy_xy, g_y_0_xy_xz, g_y_0_xy_yy, g_y_0_xy_yz, g_y_0_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_0_xy_xx[i] = 2.0 * g_y_0_xy_xx[i] * a_exp;

        g_y_0_0_0_0_0_xy_xy[i] = 2.0 * g_y_0_xy_xy[i] * a_exp;

        g_y_0_0_0_0_0_xy_xz[i] = 2.0 * g_y_0_xy_xz[i] * a_exp;

        g_y_0_0_0_0_0_xy_yy[i] = 2.0 * g_y_0_xy_yy[i] * a_exp;

        g_y_0_0_0_0_0_xy_yz[i] = 2.0 * g_y_0_xy_yz[i] * a_exp;

        g_y_0_0_0_0_0_xy_zz[i] = 2.0 * g_y_0_xy_zz[i] * a_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_y_0_0_0_0_0_xz_xx, g_y_0_0_0_0_0_xz_xy, g_y_0_0_0_0_0_xz_xz, g_y_0_0_0_0_0_xz_yy, g_y_0_0_0_0_0_xz_yz, g_y_0_0_0_0_0_xz_zz, g_y_0_xz_xx, g_y_0_xz_xy, g_y_0_xz_xz, g_y_0_xz_yy, g_y_0_xz_yz, g_y_0_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_0_xz_xx[i] = 2.0 * g_y_0_xz_xx[i] * a_exp;

        g_y_0_0_0_0_0_xz_xy[i] = 2.0 * g_y_0_xz_xy[i] * a_exp;

        g_y_0_0_0_0_0_xz_xz[i] = 2.0 * g_y_0_xz_xz[i] * a_exp;

        g_y_0_0_0_0_0_xz_yy[i] = 2.0 * g_y_0_xz_yy[i] * a_exp;

        g_y_0_0_0_0_0_xz_yz[i] = 2.0 * g_y_0_xz_yz[i] * a_exp;

        g_y_0_0_0_0_0_xz_zz[i] = 2.0 * g_y_0_xz_zz[i] * a_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_y_0_0_0_0_0_yy_xx, g_y_0_0_0_0_0_yy_xy, g_y_0_0_0_0_0_yy_xz, g_y_0_0_0_0_0_yy_yy, g_y_0_0_0_0_0_yy_yz, g_y_0_0_0_0_0_yy_zz, g_y_0_yy_xx, g_y_0_yy_xy, g_y_0_yy_xz, g_y_0_yy_yy, g_y_0_yy_yz, g_y_0_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_0_yy_xx[i] = 2.0 * g_y_0_yy_xx[i] * a_exp;

        g_y_0_0_0_0_0_yy_xy[i] = 2.0 * g_y_0_yy_xy[i] * a_exp;

        g_y_0_0_0_0_0_yy_xz[i] = 2.0 * g_y_0_yy_xz[i] * a_exp;

        g_y_0_0_0_0_0_yy_yy[i] = 2.0 * g_y_0_yy_yy[i] * a_exp;

        g_y_0_0_0_0_0_yy_yz[i] = 2.0 * g_y_0_yy_yz[i] * a_exp;

        g_y_0_0_0_0_0_yy_zz[i] = 2.0 * g_y_0_yy_zz[i] * a_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_y_0_0_0_0_0_yz_xx, g_y_0_0_0_0_0_yz_xy, g_y_0_0_0_0_0_yz_xz, g_y_0_0_0_0_0_yz_yy, g_y_0_0_0_0_0_yz_yz, g_y_0_0_0_0_0_yz_zz, g_y_0_yz_xx, g_y_0_yz_xy, g_y_0_yz_xz, g_y_0_yz_yy, g_y_0_yz_yz, g_y_0_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_0_yz_xx[i] = 2.0 * g_y_0_yz_xx[i] * a_exp;

        g_y_0_0_0_0_0_yz_xy[i] = 2.0 * g_y_0_yz_xy[i] * a_exp;

        g_y_0_0_0_0_0_yz_xz[i] = 2.0 * g_y_0_yz_xz[i] * a_exp;

        g_y_0_0_0_0_0_yz_yy[i] = 2.0 * g_y_0_yz_yy[i] * a_exp;

        g_y_0_0_0_0_0_yz_yz[i] = 2.0 * g_y_0_yz_yz[i] * a_exp;

        g_y_0_0_0_0_0_yz_zz[i] = 2.0 * g_y_0_yz_zz[i] * a_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_y_0_0_0_0_0_zz_xx, g_y_0_0_0_0_0_zz_xy, g_y_0_0_0_0_0_zz_xz, g_y_0_0_0_0_0_zz_yy, g_y_0_0_0_0_0_zz_yz, g_y_0_0_0_0_0_zz_zz, g_y_0_zz_xx, g_y_0_zz_xy, g_y_0_zz_xz, g_y_0_zz_yy, g_y_0_zz_yz, g_y_0_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_0_zz_xx[i] = 2.0 * g_y_0_zz_xx[i] * a_exp;

        g_y_0_0_0_0_0_zz_xy[i] = 2.0 * g_y_0_zz_xy[i] * a_exp;

        g_y_0_0_0_0_0_zz_xz[i] = 2.0 * g_y_0_zz_xz[i] * a_exp;

        g_y_0_0_0_0_0_zz_yy[i] = 2.0 * g_y_0_zz_yy[i] * a_exp;

        g_y_0_0_0_0_0_zz_yz[i] = 2.0 * g_y_0_zz_yz[i] * a_exp;

        g_y_0_0_0_0_0_zz_zz[i] = 2.0 * g_y_0_zz_zz[i] * a_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_z_0_0_0_0_0_xx_xx, g_z_0_0_0_0_0_xx_xy, g_z_0_0_0_0_0_xx_xz, g_z_0_0_0_0_0_xx_yy, g_z_0_0_0_0_0_xx_yz, g_z_0_0_0_0_0_xx_zz, g_z_0_xx_xx, g_z_0_xx_xy, g_z_0_xx_xz, g_z_0_xx_yy, g_z_0_xx_yz, g_z_0_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_0_xx_xx[i] = 2.0 * g_z_0_xx_xx[i] * a_exp;

        g_z_0_0_0_0_0_xx_xy[i] = 2.0 * g_z_0_xx_xy[i] * a_exp;

        g_z_0_0_0_0_0_xx_xz[i] = 2.0 * g_z_0_xx_xz[i] * a_exp;

        g_z_0_0_0_0_0_xx_yy[i] = 2.0 * g_z_0_xx_yy[i] * a_exp;

        g_z_0_0_0_0_0_xx_yz[i] = 2.0 * g_z_0_xx_yz[i] * a_exp;

        g_z_0_0_0_0_0_xx_zz[i] = 2.0 * g_z_0_xx_zz[i] * a_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_z_0_0_0_0_0_xy_xx, g_z_0_0_0_0_0_xy_xy, g_z_0_0_0_0_0_xy_xz, g_z_0_0_0_0_0_xy_yy, g_z_0_0_0_0_0_xy_yz, g_z_0_0_0_0_0_xy_zz, g_z_0_xy_xx, g_z_0_xy_xy, g_z_0_xy_xz, g_z_0_xy_yy, g_z_0_xy_yz, g_z_0_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_0_xy_xx[i] = 2.0 * g_z_0_xy_xx[i] * a_exp;

        g_z_0_0_0_0_0_xy_xy[i] = 2.0 * g_z_0_xy_xy[i] * a_exp;

        g_z_0_0_0_0_0_xy_xz[i] = 2.0 * g_z_0_xy_xz[i] * a_exp;

        g_z_0_0_0_0_0_xy_yy[i] = 2.0 * g_z_0_xy_yy[i] * a_exp;

        g_z_0_0_0_0_0_xy_yz[i] = 2.0 * g_z_0_xy_yz[i] * a_exp;

        g_z_0_0_0_0_0_xy_zz[i] = 2.0 * g_z_0_xy_zz[i] * a_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_z_0_0_0_0_0_xz_xx, g_z_0_0_0_0_0_xz_xy, g_z_0_0_0_0_0_xz_xz, g_z_0_0_0_0_0_xz_yy, g_z_0_0_0_0_0_xz_yz, g_z_0_0_0_0_0_xz_zz, g_z_0_xz_xx, g_z_0_xz_xy, g_z_0_xz_xz, g_z_0_xz_yy, g_z_0_xz_yz, g_z_0_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_0_xz_xx[i] = 2.0 * g_z_0_xz_xx[i] * a_exp;

        g_z_0_0_0_0_0_xz_xy[i] = 2.0 * g_z_0_xz_xy[i] * a_exp;

        g_z_0_0_0_0_0_xz_xz[i] = 2.0 * g_z_0_xz_xz[i] * a_exp;

        g_z_0_0_0_0_0_xz_yy[i] = 2.0 * g_z_0_xz_yy[i] * a_exp;

        g_z_0_0_0_0_0_xz_yz[i] = 2.0 * g_z_0_xz_yz[i] * a_exp;

        g_z_0_0_0_0_0_xz_zz[i] = 2.0 * g_z_0_xz_zz[i] * a_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_z_0_0_0_0_0_yy_xx, g_z_0_0_0_0_0_yy_xy, g_z_0_0_0_0_0_yy_xz, g_z_0_0_0_0_0_yy_yy, g_z_0_0_0_0_0_yy_yz, g_z_0_0_0_0_0_yy_zz, g_z_0_yy_xx, g_z_0_yy_xy, g_z_0_yy_xz, g_z_0_yy_yy, g_z_0_yy_yz, g_z_0_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_0_yy_xx[i] = 2.0 * g_z_0_yy_xx[i] * a_exp;

        g_z_0_0_0_0_0_yy_xy[i] = 2.0 * g_z_0_yy_xy[i] * a_exp;

        g_z_0_0_0_0_0_yy_xz[i] = 2.0 * g_z_0_yy_xz[i] * a_exp;

        g_z_0_0_0_0_0_yy_yy[i] = 2.0 * g_z_0_yy_yy[i] * a_exp;

        g_z_0_0_0_0_0_yy_yz[i] = 2.0 * g_z_0_yy_yz[i] * a_exp;

        g_z_0_0_0_0_0_yy_zz[i] = 2.0 * g_z_0_yy_zz[i] * a_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_z_0_0_0_0_0_yz_xx, g_z_0_0_0_0_0_yz_xy, g_z_0_0_0_0_0_yz_xz, g_z_0_0_0_0_0_yz_yy, g_z_0_0_0_0_0_yz_yz, g_z_0_0_0_0_0_yz_zz, g_z_0_yz_xx, g_z_0_yz_xy, g_z_0_yz_xz, g_z_0_yz_yy, g_z_0_yz_yz, g_z_0_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_0_yz_xx[i] = 2.0 * g_z_0_yz_xx[i] * a_exp;

        g_z_0_0_0_0_0_yz_xy[i] = 2.0 * g_z_0_yz_xy[i] * a_exp;

        g_z_0_0_0_0_0_yz_xz[i] = 2.0 * g_z_0_yz_xz[i] * a_exp;

        g_z_0_0_0_0_0_yz_yy[i] = 2.0 * g_z_0_yz_yy[i] * a_exp;

        g_z_0_0_0_0_0_yz_yz[i] = 2.0 * g_z_0_yz_yz[i] * a_exp;

        g_z_0_0_0_0_0_yz_zz[i] = 2.0 * g_z_0_yz_zz[i] * a_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_z_0_0_0_0_0_zz_xx, g_z_0_0_0_0_0_zz_xy, g_z_0_0_0_0_0_zz_xz, g_z_0_0_0_0_0_zz_yy, g_z_0_0_0_0_0_zz_yz, g_z_0_0_0_0_0_zz_zz, g_z_0_zz_xx, g_z_0_zz_xy, g_z_0_zz_xz, g_z_0_zz_yy, g_z_0_zz_yz, g_z_0_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_0_zz_xx[i] = 2.0 * g_z_0_zz_xx[i] * a_exp;

        g_z_0_0_0_0_0_zz_xy[i] = 2.0 * g_z_0_zz_xy[i] * a_exp;

        g_z_0_0_0_0_0_zz_xz[i] = 2.0 * g_z_0_zz_xz[i] * a_exp;

        g_z_0_0_0_0_0_zz_yy[i] = 2.0 * g_z_0_zz_yy[i] * a_exp;

        g_z_0_0_0_0_0_zz_yz[i] = 2.0 * g_z_0_zz_yz[i] * a_exp;

        g_z_0_0_0_0_0_zz_zz[i] = 2.0 * g_z_0_zz_zz[i] * a_exp;
    }
}

} // t4c_geom namespace


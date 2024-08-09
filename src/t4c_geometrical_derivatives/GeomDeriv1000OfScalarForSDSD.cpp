#include "GeomDeriv1000OfScalarForSDSD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1000_sdsd_0(CSimdArray<double>& buffer_1000_sdsd,
                     const CSimdArray<double>& buffer_pdsd,
                     const double a_exp) -> void
{
    const auto ndims = buffer_1000_sdsd.number_of_columns();

    /// Set up components of auxilary buffer : buffer_pdsd

    auto g_x_xx_0_xx = buffer_pdsd[0];

    auto g_x_xx_0_xy = buffer_pdsd[1];

    auto g_x_xx_0_xz = buffer_pdsd[2];

    auto g_x_xx_0_yy = buffer_pdsd[3];

    auto g_x_xx_0_yz = buffer_pdsd[4];

    auto g_x_xx_0_zz = buffer_pdsd[5];

    auto g_x_xy_0_xx = buffer_pdsd[6];

    auto g_x_xy_0_xy = buffer_pdsd[7];

    auto g_x_xy_0_xz = buffer_pdsd[8];

    auto g_x_xy_0_yy = buffer_pdsd[9];

    auto g_x_xy_0_yz = buffer_pdsd[10];

    auto g_x_xy_0_zz = buffer_pdsd[11];

    auto g_x_xz_0_xx = buffer_pdsd[12];

    auto g_x_xz_0_xy = buffer_pdsd[13];

    auto g_x_xz_0_xz = buffer_pdsd[14];

    auto g_x_xz_0_yy = buffer_pdsd[15];

    auto g_x_xz_0_yz = buffer_pdsd[16];

    auto g_x_xz_0_zz = buffer_pdsd[17];

    auto g_x_yy_0_xx = buffer_pdsd[18];

    auto g_x_yy_0_xy = buffer_pdsd[19];

    auto g_x_yy_0_xz = buffer_pdsd[20];

    auto g_x_yy_0_yy = buffer_pdsd[21];

    auto g_x_yy_0_yz = buffer_pdsd[22];

    auto g_x_yy_0_zz = buffer_pdsd[23];

    auto g_x_yz_0_xx = buffer_pdsd[24];

    auto g_x_yz_0_xy = buffer_pdsd[25];

    auto g_x_yz_0_xz = buffer_pdsd[26];

    auto g_x_yz_0_yy = buffer_pdsd[27];

    auto g_x_yz_0_yz = buffer_pdsd[28];

    auto g_x_yz_0_zz = buffer_pdsd[29];

    auto g_x_zz_0_xx = buffer_pdsd[30];

    auto g_x_zz_0_xy = buffer_pdsd[31];

    auto g_x_zz_0_xz = buffer_pdsd[32];

    auto g_x_zz_0_yy = buffer_pdsd[33];

    auto g_x_zz_0_yz = buffer_pdsd[34];

    auto g_x_zz_0_zz = buffer_pdsd[35];

    auto g_y_xx_0_xx = buffer_pdsd[36];

    auto g_y_xx_0_xy = buffer_pdsd[37];

    auto g_y_xx_0_xz = buffer_pdsd[38];

    auto g_y_xx_0_yy = buffer_pdsd[39];

    auto g_y_xx_0_yz = buffer_pdsd[40];

    auto g_y_xx_0_zz = buffer_pdsd[41];

    auto g_y_xy_0_xx = buffer_pdsd[42];

    auto g_y_xy_0_xy = buffer_pdsd[43];

    auto g_y_xy_0_xz = buffer_pdsd[44];

    auto g_y_xy_0_yy = buffer_pdsd[45];

    auto g_y_xy_0_yz = buffer_pdsd[46];

    auto g_y_xy_0_zz = buffer_pdsd[47];

    auto g_y_xz_0_xx = buffer_pdsd[48];

    auto g_y_xz_0_xy = buffer_pdsd[49];

    auto g_y_xz_0_xz = buffer_pdsd[50];

    auto g_y_xz_0_yy = buffer_pdsd[51];

    auto g_y_xz_0_yz = buffer_pdsd[52];

    auto g_y_xz_0_zz = buffer_pdsd[53];

    auto g_y_yy_0_xx = buffer_pdsd[54];

    auto g_y_yy_0_xy = buffer_pdsd[55];

    auto g_y_yy_0_xz = buffer_pdsd[56];

    auto g_y_yy_0_yy = buffer_pdsd[57];

    auto g_y_yy_0_yz = buffer_pdsd[58];

    auto g_y_yy_0_zz = buffer_pdsd[59];

    auto g_y_yz_0_xx = buffer_pdsd[60];

    auto g_y_yz_0_xy = buffer_pdsd[61];

    auto g_y_yz_0_xz = buffer_pdsd[62];

    auto g_y_yz_0_yy = buffer_pdsd[63];

    auto g_y_yz_0_yz = buffer_pdsd[64];

    auto g_y_yz_0_zz = buffer_pdsd[65];

    auto g_y_zz_0_xx = buffer_pdsd[66];

    auto g_y_zz_0_xy = buffer_pdsd[67];

    auto g_y_zz_0_xz = buffer_pdsd[68];

    auto g_y_zz_0_yy = buffer_pdsd[69];

    auto g_y_zz_0_yz = buffer_pdsd[70];

    auto g_y_zz_0_zz = buffer_pdsd[71];

    auto g_z_xx_0_xx = buffer_pdsd[72];

    auto g_z_xx_0_xy = buffer_pdsd[73];

    auto g_z_xx_0_xz = buffer_pdsd[74];

    auto g_z_xx_0_yy = buffer_pdsd[75];

    auto g_z_xx_0_yz = buffer_pdsd[76];

    auto g_z_xx_0_zz = buffer_pdsd[77];

    auto g_z_xy_0_xx = buffer_pdsd[78];

    auto g_z_xy_0_xy = buffer_pdsd[79];

    auto g_z_xy_0_xz = buffer_pdsd[80];

    auto g_z_xy_0_yy = buffer_pdsd[81];

    auto g_z_xy_0_yz = buffer_pdsd[82];

    auto g_z_xy_0_zz = buffer_pdsd[83];

    auto g_z_xz_0_xx = buffer_pdsd[84];

    auto g_z_xz_0_xy = buffer_pdsd[85];

    auto g_z_xz_0_xz = buffer_pdsd[86];

    auto g_z_xz_0_yy = buffer_pdsd[87];

    auto g_z_xz_0_yz = buffer_pdsd[88];

    auto g_z_xz_0_zz = buffer_pdsd[89];

    auto g_z_yy_0_xx = buffer_pdsd[90];

    auto g_z_yy_0_xy = buffer_pdsd[91];

    auto g_z_yy_0_xz = buffer_pdsd[92];

    auto g_z_yy_0_yy = buffer_pdsd[93];

    auto g_z_yy_0_yz = buffer_pdsd[94];

    auto g_z_yy_0_zz = buffer_pdsd[95];

    auto g_z_yz_0_xx = buffer_pdsd[96];

    auto g_z_yz_0_xy = buffer_pdsd[97];

    auto g_z_yz_0_xz = buffer_pdsd[98];

    auto g_z_yz_0_yy = buffer_pdsd[99];

    auto g_z_yz_0_yz = buffer_pdsd[100];

    auto g_z_yz_0_zz = buffer_pdsd[101];

    auto g_z_zz_0_xx = buffer_pdsd[102];

    auto g_z_zz_0_xy = buffer_pdsd[103];

    auto g_z_zz_0_xz = buffer_pdsd[104];

    auto g_z_zz_0_yy = buffer_pdsd[105];

    auto g_z_zz_0_yz = buffer_pdsd[106];

    auto g_z_zz_0_zz = buffer_pdsd[107];

    /// Set up components of integrals buffer : buffer_1000_sdsd

    auto g_x_0_0_0_0_xx_0_xx = buffer_1000_sdsd[0];

    auto g_x_0_0_0_0_xx_0_xy = buffer_1000_sdsd[1];

    auto g_x_0_0_0_0_xx_0_xz = buffer_1000_sdsd[2];

    auto g_x_0_0_0_0_xx_0_yy = buffer_1000_sdsd[3];

    auto g_x_0_0_0_0_xx_0_yz = buffer_1000_sdsd[4];

    auto g_x_0_0_0_0_xx_0_zz = buffer_1000_sdsd[5];

    auto g_x_0_0_0_0_xy_0_xx = buffer_1000_sdsd[6];

    auto g_x_0_0_0_0_xy_0_xy = buffer_1000_sdsd[7];

    auto g_x_0_0_0_0_xy_0_xz = buffer_1000_sdsd[8];

    auto g_x_0_0_0_0_xy_0_yy = buffer_1000_sdsd[9];

    auto g_x_0_0_0_0_xy_0_yz = buffer_1000_sdsd[10];

    auto g_x_0_0_0_0_xy_0_zz = buffer_1000_sdsd[11];

    auto g_x_0_0_0_0_xz_0_xx = buffer_1000_sdsd[12];

    auto g_x_0_0_0_0_xz_0_xy = buffer_1000_sdsd[13];

    auto g_x_0_0_0_0_xz_0_xz = buffer_1000_sdsd[14];

    auto g_x_0_0_0_0_xz_0_yy = buffer_1000_sdsd[15];

    auto g_x_0_0_0_0_xz_0_yz = buffer_1000_sdsd[16];

    auto g_x_0_0_0_0_xz_0_zz = buffer_1000_sdsd[17];

    auto g_x_0_0_0_0_yy_0_xx = buffer_1000_sdsd[18];

    auto g_x_0_0_0_0_yy_0_xy = buffer_1000_sdsd[19];

    auto g_x_0_0_0_0_yy_0_xz = buffer_1000_sdsd[20];

    auto g_x_0_0_0_0_yy_0_yy = buffer_1000_sdsd[21];

    auto g_x_0_0_0_0_yy_0_yz = buffer_1000_sdsd[22];

    auto g_x_0_0_0_0_yy_0_zz = buffer_1000_sdsd[23];

    auto g_x_0_0_0_0_yz_0_xx = buffer_1000_sdsd[24];

    auto g_x_0_0_0_0_yz_0_xy = buffer_1000_sdsd[25];

    auto g_x_0_0_0_0_yz_0_xz = buffer_1000_sdsd[26];

    auto g_x_0_0_0_0_yz_0_yy = buffer_1000_sdsd[27];

    auto g_x_0_0_0_0_yz_0_yz = buffer_1000_sdsd[28];

    auto g_x_0_0_0_0_yz_0_zz = buffer_1000_sdsd[29];

    auto g_x_0_0_0_0_zz_0_xx = buffer_1000_sdsd[30];

    auto g_x_0_0_0_0_zz_0_xy = buffer_1000_sdsd[31];

    auto g_x_0_0_0_0_zz_0_xz = buffer_1000_sdsd[32];

    auto g_x_0_0_0_0_zz_0_yy = buffer_1000_sdsd[33];

    auto g_x_0_0_0_0_zz_0_yz = buffer_1000_sdsd[34];

    auto g_x_0_0_0_0_zz_0_zz = buffer_1000_sdsd[35];

    auto g_y_0_0_0_0_xx_0_xx = buffer_1000_sdsd[36];

    auto g_y_0_0_0_0_xx_0_xy = buffer_1000_sdsd[37];

    auto g_y_0_0_0_0_xx_0_xz = buffer_1000_sdsd[38];

    auto g_y_0_0_0_0_xx_0_yy = buffer_1000_sdsd[39];

    auto g_y_0_0_0_0_xx_0_yz = buffer_1000_sdsd[40];

    auto g_y_0_0_0_0_xx_0_zz = buffer_1000_sdsd[41];

    auto g_y_0_0_0_0_xy_0_xx = buffer_1000_sdsd[42];

    auto g_y_0_0_0_0_xy_0_xy = buffer_1000_sdsd[43];

    auto g_y_0_0_0_0_xy_0_xz = buffer_1000_sdsd[44];

    auto g_y_0_0_0_0_xy_0_yy = buffer_1000_sdsd[45];

    auto g_y_0_0_0_0_xy_0_yz = buffer_1000_sdsd[46];

    auto g_y_0_0_0_0_xy_0_zz = buffer_1000_sdsd[47];

    auto g_y_0_0_0_0_xz_0_xx = buffer_1000_sdsd[48];

    auto g_y_0_0_0_0_xz_0_xy = buffer_1000_sdsd[49];

    auto g_y_0_0_0_0_xz_0_xz = buffer_1000_sdsd[50];

    auto g_y_0_0_0_0_xz_0_yy = buffer_1000_sdsd[51];

    auto g_y_0_0_0_0_xz_0_yz = buffer_1000_sdsd[52];

    auto g_y_0_0_0_0_xz_0_zz = buffer_1000_sdsd[53];

    auto g_y_0_0_0_0_yy_0_xx = buffer_1000_sdsd[54];

    auto g_y_0_0_0_0_yy_0_xy = buffer_1000_sdsd[55];

    auto g_y_0_0_0_0_yy_0_xz = buffer_1000_sdsd[56];

    auto g_y_0_0_0_0_yy_0_yy = buffer_1000_sdsd[57];

    auto g_y_0_0_0_0_yy_0_yz = buffer_1000_sdsd[58];

    auto g_y_0_0_0_0_yy_0_zz = buffer_1000_sdsd[59];

    auto g_y_0_0_0_0_yz_0_xx = buffer_1000_sdsd[60];

    auto g_y_0_0_0_0_yz_0_xy = buffer_1000_sdsd[61];

    auto g_y_0_0_0_0_yz_0_xz = buffer_1000_sdsd[62];

    auto g_y_0_0_0_0_yz_0_yy = buffer_1000_sdsd[63];

    auto g_y_0_0_0_0_yz_0_yz = buffer_1000_sdsd[64];

    auto g_y_0_0_0_0_yz_0_zz = buffer_1000_sdsd[65];

    auto g_y_0_0_0_0_zz_0_xx = buffer_1000_sdsd[66];

    auto g_y_0_0_0_0_zz_0_xy = buffer_1000_sdsd[67];

    auto g_y_0_0_0_0_zz_0_xz = buffer_1000_sdsd[68];

    auto g_y_0_0_0_0_zz_0_yy = buffer_1000_sdsd[69];

    auto g_y_0_0_0_0_zz_0_yz = buffer_1000_sdsd[70];

    auto g_y_0_0_0_0_zz_0_zz = buffer_1000_sdsd[71];

    auto g_z_0_0_0_0_xx_0_xx = buffer_1000_sdsd[72];

    auto g_z_0_0_0_0_xx_0_xy = buffer_1000_sdsd[73];

    auto g_z_0_0_0_0_xx_0_xz = buffer_1000_sdsd[74];

    auto g_z_0_0_0_0_xx_0_yy = buffer_1000_sdsd[75];

    auto g_z_0_0_0_0_xx_0_yz = buffer_1000_sdsd[76];

    auto g_z_0_0_0_0_xx_0_zz = buffer_1000_sdsd[77];

    auto g_z_0_0_0_0_xy_0_xx = buffer_1000_sdsd[78];

    auto g_z_0_0_0_0_xy_0_xy = buffer_1000_sdsd[79];

    auto g_z_0_0_0_0_xy_0_xz = buffer_1000_sdsd[80];

    auto g_z_0_0_0_0_xy_0_yy = buffer_1000_sdsd[81];

    auto g_z_0_0_0_0_xy_0_yz = buffer_1000_sdsd[82];

    auto g_z_0_0_0_0_xy_0_zz = buffer_1000_sdsd[83];

    auto g_z_0_0_0_0_xz_0_xx = buffer_1000_sdsd[84];

    auto g_z_0_0_0_0_xz_0_xy = buffer_1000_sdsd[85];

    auto g_z_0_0_0_0_xz_0_xz = buffer_1000_sdsd[86];

    auto g_z_0_0_0_0_xz_0_yy = buffer_1000_sdsd[87];

    auto g_z_0_0_0_0_xz_0_yz = buffer_1000_sdsd[88];

    auto g_z_0_0_0_0_xz_0_zz = buffer_1000_sdsd[89];

    auto g_z_0_0_0_0_yy_0_xx = buffer_1000_sdsd[90];

    auto g_z_0_0_0_0_yy_0_xy = buffer_1000_sdsd[91];

    auto g_z_0_0_0_0_yy_0_xz = buffer_1000_sdsd[92];

    auto g_z_0_0_0_0_yy_0_yy = buffer_1000_sdsd[93];

    auto g_z_0_0_0_0_yy_0_yz = buffer_1000_sdsd[94];

    auto g_z_0_0_0_0_yy_0_zz = buffer_1000_sdsd[95];

    auto g_z_0_0_0_0_yz_0_xx = buffer_1000_sdsd[96];

    auto g_z_0_0_0_0_yz_0_xy = buffer_1000_sdsd[97];

    auto g_z_0_0_0_0_yz_0_xz = buffer_1000_sdsd[98];

    auto g_z_0_0_0_0_yz_0_yy = buffer_1000_sdsd[99];

    auto g_z_0_0_0_0_yz_0_yz = buffer_1000_sdsd[100];

    auto g_z_0_0_0_0_yz_0_zz = buffer_1000_sdsd[101];

    auto g_z_0_0_0_0_zz_0_xx = buffer_1000_sdsd[102];

    auto g_z_0_0_0_0_zz_0_xy = buffer_1000_sdsd[103];

    auto g_z_0_0_0_0_zz_0_xz = buffer_1000_sdsd[104];

    auto g_z_0_0_0_0_zz_0_yy = buffer_1000_sdsd[105];

    auto g_z_0_0_0_0_zz_0_yz = buffer_1000_sdsd[106];

    auto g_z_0_0_0_0_zz_0_zz = buffer_1000_sdsd[107];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_0_0_0_0_xx_0_xx, g_x_0_0_0_0_xx_0_xy, g_x_0_0_0_0_xx_0_xz, g_x_0_0_0_0_xx_0_yy, g_x_0_0_0_0_xx_0_yz, g_x_0_0_0_0_xx_0_zz, g_x_xx_0_xx, g_x_xx_0_xy, g_x_xx_0_xz, g_x_xx_0_yy, g_x_xx_0_yz, g_x_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xx_0_xx[i] = 2.0 * g_x_xx_0_xx[i] * a_exp;

        g_x_0_0_0_0_xx_0_xy[i] = 2.0 * g_x_xx_0_xy[i] * a_exp;

        g_x_0_0_0_0_xx_0_xz[i] = 2.0 * g_x_xx_0_xz[i] * a_exp;

        g_x_0_0_0_0_xx_0_yy[i] = 2.0 * g_x_xx_0_yy[i] * a_exp;

        g_x_0_0_0_0_xx_0_yz[i] = 2.0 * g_x_xx_0_yz[i] * a_exp;

        g_x_0_0_0_0_xx_0_zz[i] = 2.0 * g_x_xx_0_zz[i] * a_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_0_0_0_0_xy_0_xx, g_x_0_0_0_0_xy_0_xy, g_x_0_0_0_0_xy_0_xz, g_x_0_0_0_0_xy_0_yy, g_x_0_0_0_0_xy_0_yz, g_x_0_0_0_0_xy_0_zz, g_x_xy_0_xx, g_x_xy_0_xy, g_x_xy_0_xz, g_x_xy_0_yy, g_x_xy_0_yz, g_x_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xy_0_xx[i] = 2.0 * g_x_xy_0_xx[i] * a_exp;

        g_x_0_0_0_0_xy_0_xy[i] = 2.0 * g_x_xy_0_xy[i] * a_exp;

        g_x_0_0_0_0_xy_0_xz[i] = 2.0 * g_x_xy_0_xz[i] * a_exp;

        g_x_0_0_0_0_xy_0_yy[i] = 2.0 * g_x_xy_0_yy[i] * a_exp;

        g_x_0_0_0_0_xy_0_yz[i] = 2.0 * g_x_xy_0_yz[i] * a_exp;

        g_x_0_0_0_0_xy_0_zz[i] = 2.0 * g_x_xy_0_zz[i] * a_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_0_0_0_0_xz_0_xx, g_x_0_0_0_0_xz_0_xy, g_x_0_0_0_0_xz_0_xz, g_x_0_0_0_0_xz_0_yy, g_x_0_0_0_0_xz_0_yz, g_x_0_0_0_0_xz_0_zz, g_x_xz_0_xx, g_x_xz_0_xy, g_x_xz_0_xz, g_x_xz_0_yy, g_x_xz_0_yz, g_x_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xz_0_xx[i] = 2.0 * g_x_xz_0_xx[i] * a_exp;

        g_x_0_0_0_0_xz_0_xy[i] = 2.0 * g_x_xz_0_xy[i] * a_exp;

        g_x_0_0_0_0_xz_0_xz[i] = 2.0 * g_x_xz_0_xz[i] * a_exp;

        g_x_0_0_0_0_xz_0_yy[i] = 2.0 * g_x_xz_0_yy[i] * a_exp;

        g_x_0_0_0_0_xz_0_yz[i] = 2.0 * g_x_xz_0_yz[i] * a_exp;

        g_x_0_0_0_0_xz_0_zz[i] = 2.0 * g_x_xz_0_zz[i] * a_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_x_0_0_0_0_yy_0_xx, g_x_0_0_0_0_yy_0_xy, g_x_0_0_0_0_yy_0_xz, g_x_0_0_0_0_yy_0_yy, g_x_0_0_0_0_yy_0_yz, g_x_0_0_0_0_yy_0_zz, g_x_yy_0_xx, g_x_yy_0_xy, g_x_yy_0_xz, g_x_yy_0_yy, g_x_yy_0_yz, g_x_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_yy_0_xx[i] = 2.0 * g_x_yy_0_xx[i] * a_exp;

        g_x_0_0_0_0_yy_0_xy[i] = 2.0 * g_x_yy_0_xy[i] * a_exp;

        g_x_0_0_0_0_yy_0_xz[i] = 2.0 * g_x_yy_0_xz[i] * a_exp;

        g_x_0_0_0_0_yy_0_yy[i] = 2.0 * g_x_yy_0_yy[i] * a_exp;

        g_x_0_0_0_0_yy_0_yz[i] = 2.0 * g_x_yy_0_yz[i] * a_exp;

        g_x_0_0_0_0_yy_0_zz[i] = 2.0 * g_x_yy_0_zz[i] * a_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_x_0_0_0_0_yz_0_xx, g_x_0_0_0_0_yz_0_xy, g_x_0_0_0_0_yz_0_xz, g_x_0_0_0_0_yz_0_yy, g_x_0_0_0_0_yz_0_yz, g_x_0_0_0_0_yz_0_zz, g_x_yz_0_xx, g_x_yz_0_xy, g_x_yz_0_xz, g_x_yz_0_yy, g_x_yz_0_yz, g_x_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_yz_0_xx[i] = 2.0 * g_x_yz_0_xx[i] * a_exp;

        g_x_0_0_0_0_yz_0_xy[i] = 2.0 * g_x_yz_0_xy[i] * a_exp;

        g_x_0_0_0_0_yz_0_xz[i] = 2.0 * g_x_yz_0_xz[i] * a_exp;

        g_x_0_0_0_0_yz_0_yy[i] = 2.0 * g_x_yz_0_yy[i] * a_exp;

        g_x_0_0_0_0_yz_0_yz[i] = 2.0 * g_x_yz_0_yz[i] * a_exp;

        g_x_0_0_0_0_yz_0_zz[i] = 2.0 * g_x_yz_0_zz[i] * a_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_x_0_0_0_0_zz_0_xx, g_x_0_0_0_0_zz_0_xy, g_x_0_0_0_0_zz_0_xz, g_x_0_0_0_0_zz_0_yy, g_x_0_0_0_0_zz_0_yz, g_x_0_0_0_0_zz_0_zz, g_x_zz_0_xx, g_x_zz_0_xy, g_x_zz_0_xz, g_x_zz_0_yy, g_x_zz_0_yz, g_x_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_zz_0_xx[i] = 2.0 * g_x_zz_0_xx[i] * a_exp;

        g_x_0_0_0_0_zz_0_xy[i] = 2.0 * g_x_zz_0_xy[i] * a_exp;

        g_x_0_0_0_0_zz_0_xz[i] = 2.0 * g_x_zz_0_xz[i] * a_exp;

        g_x_0_0_0_0_zz_0_yy[i] = 2.0 * g_x_zz_0_yy[i] * a_exp;

        g_x_0_0_0_0_zz_0_yz[i] = 2.0 * g_x_zz_0_yz[i] * a_exp;

        g_x_0_0_0_0_zz_0_zz[i] = 2.0 * g_x_zz_0_zz[i] * a_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_y_0_0_0_0_xx_0_xx, g_y_0_0_0_0_xx_0_xy, g_y_0_0_0_0_xx_0_xz, g_y_0_0_0_0_xx_0_yy, g_y_0_0_0_0_xx_0_yz, g_y_0_0_0_0_xx_0_zz, g_y_xx_0_xx, g_y_xx_0_xy, g_y_xx_0_xz, g_y_xx_0_yy, g_y_xx_0_yz, g_y_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xx_0_xx[i] = 2.0 * g_y_xx_0_xx[i] * a_exp;

        g_y_0_0_0_0_xx_0_xy[i] = 2.0 * g_y_xx_0_xy[i] * a_exp;

        g_y_0_0_0_0_xx_0_xz[i] = 2.0 * g_y_xx_0_xz[i] * a_exp;

        g_y_0_0_0_0_xx_0_yy[i] = 2.0 * g_y_xx_0_yy[i] * a_exp;

        g_y_0_0_0_0_xx_0_yz[i] = 2.0 * g_y_xx_0_yz[i] * a_exp;

        g_y_0_0_0_0_xx_0_zz[i] = 2.0 * g_y_xx_0_zz[i] * a_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_y_0_0_0_0_xy_0_xx, g_y_0_0_0_0_xy_0_xy, g_y_0_0_0_0_xy_0_xz, g_y_0_0_0_0_xy_0_yy, g_y_0_0_0_0_xy_0_yz, g_y_0_0_0_0_xy_0_zz, g_y_xy_0_xx, g_y_xy_0_xy, g_y_xy_0_xz, g_y_xy_0_yy, g_y_xy_0_yz, g_y_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xy_0_xx[i] = 2.0 * g_y_xy_0_xx[i] * a_exp;

        g_y_0_0_0_0_xy_0_xy[i] = 2.0 * g_y_xy_0_xy[i] * a_exp;

        g_y_0_0_0_0_xy_0_xz[i] = 2.0 * g_y_xy_0_xz[i] * a_exp;

        g_y_0_0_0_0_xy_0_yy[i] = 2.0 * g_y_xy_0_yy[i] * a_exp;

        g_y_0_0_0_0_xy_0_yz[i] = 2.0 * g_y_xy_0_yz[i] * a_exp;

        g_y_0_0_0_0_xy_0_zz[i] = 2.0 * g_y_xy_0_zz[i] * a_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_y_0_0_0_0_xz_0_xx, g_y_0_0_0_0_xz_0_xy, g_y_0_0_0_0_xz_0_xz, g_y_0_0_0_0_xz_0_yy, g_y_0_0_0_0_xz_0_yz, g_y_0_0_0_0_xz_0_zz, g_y_xz_0_xx, g_y_xz_0_xy, g_y_xz_0_xz, g_y_xz_0_yy, g_y_xz_0_yz, g_y_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xz_0_xx[i] = 2.0 * g_y_xz_0_xx[i] * a_exp;

        g_y_0_0_0_0_xz_0_xy[i] = 2.0 * g_y_xz_0_xy[i] * a_exp;

        g_y_0_0_0_0_xz_0_xz[i] = 2.0 * g_y_xz_0_xz[i] * a_exp;

        g_y_0_0_0_0_xz_0_yy[i] = 2.0 * g_y_xz_0_yy[i] * a_exp;

        g_y_0_0_0_0_xz_0_yz[i] = 2.0 * g_y_xz_0_yz[i] * a_exp;

        g_y_0_0_0_0_xz_0_zz[i] = 2.0 * g_y_xz_0_zz[i] * a_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_y_0_0_0_0_yy_0_xx, g_y_0_0_0_0_yy_0_xy, g_y_0_0_0_0_yy_0_xz, g_y_0_0_0_0_yy_0_yy, g_y_0_0_0_0_yy_0_yz, g_y_0_0_0_0_yy_0_zz, g_y_yy_0_xx, g_y_yy_0_xy, g_y_yy_0_xz, g_y_yy_0_yy, g_y_yy_0_yz, g_y_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_yy_0_xx[i] = 2.0 * g_y_yy_0_xx[i] * a_exp;

        g_y_0_0_0_0_yy_0_xy[i] = 2.0 * g_y_yy_0_xy[i] * a_exp;

        g_y_0_0_0_0_yy_0_xz[i] = 2.0 * g_y_yy_0_xz[i] * a_exp;

        g_y_0_0_0_0_yy_0_yy[i] = 2.0 * g_y_yy_0_yy[i] * a_exp;

        g_y_0_0_0_0_yy_0_yz[i] = 2.0 * g_y_yy_0_yz[i] * a_exp;

        g_y_0_0_0_0_yy_0_zz[i] = 2.0 * g_y_yy_0_zz[i] * a_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_y_0_0_0_0_yz_0_xx, g_y_0_0_0_0_yz_0_xy, g_y_0_0_0_0_yz_0_xz, g_y_0_0_0_0_yz_0_yy, g_y_0_0_0_0_yz_0_yz, g_y_0_0_0_0_yz_0_zz, g_y_yz_0_xx, g_y_yz_0_xy, g_y_yz_0_xz, g_y_yz_0_yy, g_y_yz_0_yz, g_y_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_yz_0_xx[i] = 2.0 * g_y_yz_0_xx[i] * a_exp;

        g_y_0_0_0_0_yz_0_xy[i] = 2.0 * g_y_yz_0_xy[i] * a_exp;

        g_y_0_0_0_0_yz_0_xz[i] = 2.0 * g_y_yz_0_xz[i] * a_exp;

        g_y_0_0_0_0_yz_0_yy[i] = 2.0 * g_y_yz_0_yy[i] * a_exp;

        g_y_0_0_0_0_yz_0_yz[i] = 2.0 * g_y_yz_0_yz[i] * a_exp;

        g_y_0_0_0_0_yz_0_zz[i] = 2.0 * g_y_yz_0_zz[i] * a_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_y_0_0_0_0_zz_0_xx, g_y_0_0_0_0_zz_0_xy, g_y_0_0_0_0_zz_0_xz, g_y_0_0_0_0_zz_0_yy, g_y_0_0_0_0_zz_0_yz, g_y_0_0_0_0_zz_0_zz, g_y_zz_0_xx, g_y_zz_0_xy, g_y_zz_0_xz, g_y_zz_0_yy, g_y_zz_0_yz, g_y_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_zz_0_xx[i] = 2.0 * g_y_zz_0_xx[i] * a_exp;

        g_y_0_0_0_0_zz_0_xy[i] = 2.0 * g_y_zz_0_xy[i] * a_exp;

        g_y_0_0_0_0_zz_0_xz[i] = 2.0 * g_y_zz_0_xz[i] * a_exp;

        g_y_0_0_0_0_zz_0_yy[i] = 2.0 * g_y_zz_0_yy[i] * a_exp;

        g_y_0_0_0_0_zz_0_yz[i] = 2.0 * g_y_zz_0_yz[i] * a_exp;

        g_y_0_0_0_0_zz_0_zz[i] = 2.0 * g_y_zz_0_zz[i] * a_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_z_0_0_0_0_xx_0_xx, g_z_0_0_0_0_xx_0_xy, g_z_0_0_0_0_xx_0_xz, g_z_0_0_0_0_xx_0_yy, g_z_0_0_0_0_xx_0_yz, g_z_0_0_0_0_xx_0_zz, g_z_xx_0_xx, g_z_xx_0_xy, g_z_xx_0_xz, g_z_xx_0_yy, g_z_xx_0_yz, g_z_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xx_0_xx[i] = 2.0 * g_z_xx_0_xx[i] * a_exp;

        g_z_0_0_0_0_xx_0_xy[i] = 2.0 * g_z_xx_0_xy[i] * a_exp;

        g_z_0_0_0_0_xx_0_xz[i] = 2.0 * g_z_xx_0_xz[i] * a_exp;

        g_z_0_0_0_0_xx_0_yy[i] = 2.0 * g_z_xx_0_yy[i] * a_exp;

        g_z_0_0_0_0_xx_0_yz[i] = 2.0 * g_z_xx_0_yz[i] * a_exp;

        g_z_0_0_0_0_xx_0_zz[i] = 2.0 * g_z_xx_0_zz[i] * a_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_z_0_0_0_0_xy_0_xx, g_z_0_0_0_0_xy_0_xy, g_z_0_0_0_0_xy_0_xz, g_z_0_0_0_0_xy_0_yy, g_z_0_0_0_0_xy_0_yz, g_z_0_0_0_0_xy_0_zz, g_z_xy_0_xx, g_z_xy_0_xy, g_z_xy_0_xz, g_z_xy_0_yy, g_z_xy_0_yz, g_z_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xy_0_xx[i] = 2.0 * g_z_xy_0_xx[i] * a_exp;

        g_z_0_0_0_0_xy_0_xy[i] = 2.0 * g_z_xy_0_xy[i] * a_exp;

        g_z_0_0_0_0_xy_0_xz[i] = 2.0 * g_z_xy_0_xz[i] * a_exp;

        g_z_0_0_0_0_xy_0_yy[i] = 2.0 * g_z_xy_0_yy[i] * a_exp;

        g_z_0_0_0_0_xy_0_yz[i] = 2.0 * g_z_xy_0_yz[i] * a_exp;

        g_z_0_0_0_0_xy_0_zz[i] = 2.0 * g_z_xy_0_zz[i] * a_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_z_0_0_0_0_xz_0_xx, g_z_0_0_0_0_xz_0_xy, g_z_0_0_0_0_xz_0_xz, g_z_0_0_0_0_xz_0_yy, g_z_0_0_0_0_xz_0_yz, g_z_0_0_0_0_xz_0_zz, g_z_xz_0_xx, g_z_xz_0_xy, g_z_xz_0_xz, g_z_xz_0_yy, g_z_xz_0_yz, g_z_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xz_0_xx[i] = 2.0 * g_z_xz_0_xx[i] * a_exp;

        g_z_0_0_0_0_xz_0_xy[i] = 2.0 * g_z_xz_0_xy[i] * a_exp;

        g_z_0_0_0_0_xz_0_xz[i] = 2.0 * g_z_xz_0_xz[i] * a_exp;

        g_z_0_0_0_0_xz_0_yy[i] = 2.0 * g_z_xz_0_yy[i] * a_exp;

        g_z_0_0_0_0_xz_0_yz[i] = 2.0 * g_z_xz_0_yz[i] * a_exp;

        g_z_0_0_0_0_xz_0_zz[i] = 2.0 * g_z_xz_0_zz[i] * a_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_z_0_0_0_0_yy_0_xx, g_z_0_0_0_0_yy_0_xy, g_z_0_0_0_0_yy_0_xz, g_z_0_0_0_0_yy_0_yy, g_z_0_0_0_0_yy_0_yz, g_z_0_0_0_0_yy_0_zz, g_z_yy_0_xx, g_z_yy_0_xy, g_z_yy_0_xz, g_z_yy_0_yy, g_z_yy_0_yz, g_z_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_yy_0_xx[i] = 2.0 * g_z_yy_0_xx[i] * a_exp;

        g_z_0_0_0_0_yy_0_xy[i] = 2.0 * g_z_yy_0_xy[i] * a_exp;

        g_z_0_0_0_0_yy_0_xz[i] = 2.0 * g_z_yy_0_xz[i] * a_exp;

        g_z_0_0_0_0_yy_0_yy[i] = 2.0 * g_z_yy_0_yy[i] * a_exp;

        g_z_0_0_0_0_yy_0_yz[i] = 2.0 * g_z_yy_0_yz[i] * a_exp;

        g_z_0_0_0_0_yy_0_zz[i] = 2.0 * g_z_yy_0_zz[i] * a_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_z_0_0_0_0_yz_0_xx, g_z_0_0_0_0_yz_0_xy, g_z_0_0_0_0_yz_0_xz, g_z_0_0_0_0_yz_0_yy, g_z_0_0_0_0_yz_0_yz, g_z_0_0_0_0_yz_0_zz, g_z_yz_0_xx, g_z_yz_0_xy, g_z_yz_0_xz, g_z_yz_0_yy, g_z_yz_0_yz, g_z_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_yz_0_xx[i] = 2.0 * g_z_yz_0_xx[i] * a_exp;

        g_z_0_0_0_0_yz_0_xy[i] = 2.0 * g_z_yz_0_xy[i] * a_exp;

        g_z_0_0_0_0_yz_0_xz[i] = 2.0 * g_z_yz_0_xz[i] * a_exp;

        g_z_0_0_0_0_yz_0_yy[i] = 2.0 * g_z_yz_0_yy[i] * a_exp;

        g_z_0_0_0_0_yz_0_yz[i] = 2.0 * g_z_yz_0_yz[i] * a_exp;

        g_z_0_0_0_0_yz_0_zz[i] = 2.0 * g_z_yz_0_zz[i] * a_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_z_0_0_0_0_zz_0_xx, g_z_0_0_0_0_zz_0_xy, g_z_0_0_0_0_zz_0_xz, g_z_0_0_0_0_zz_0_yy, g_z_0_0_0_0_zz_0_yz, g_z_0_0_0_0_zz_0_zz, g_z_zz_0_xx, g_z_zz_0_xy, g_z_zz_0_xz, g_z_zz_0_yy, g_z_zz_0_yz, g_z_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_zz_0_xx[i] = 2.0 * g_z_zz_0_xx[i] * a_exp;

        g_z_0_0_0_0_zz_0_xy[i] = 2.0 * g_z_zz_0_xy[i] * a_exp;

        g_z_0_0_0_0_zz_0_xz[i] = 2.0 * g_z_zz_0_xz[i] * a_exp;

        g_z_0_0_0_0_zz_0_yy[i] = 2.0 * g_z_zz_0_yy[i] * a_exp;

        g_z_0_0_0_0_zz_0_yz[i] = 2.0 * g_z_zz_0_yz[i] * a_exp;

        g_z_0_0_0_0_zz_0_zz[i] = 2.0 * g_z_zz_0_zz[i] * a_exp;
    }
}

} // t4c_geom namespace


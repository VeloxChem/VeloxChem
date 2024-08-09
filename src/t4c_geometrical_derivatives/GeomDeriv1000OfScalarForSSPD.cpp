#include "GeomDeriv1000OfScalarForSSPD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1000_sspd_0(CSimdArray<double>& buffer_1000_sspd,
                     const CSimdArray<double>& buffer_pspd,
                     const double a_exp) -> void
{
    const auto ndims = buffer_1000_sspd.number_of_columns();

    /// Set up components of auxilary buffer : buffer_pspd

    auto g_x_0_x_xx = buffer_pspd[0];

    auto g_x_0_x_xy = buffer_pspd[1];

    auto g_x_0_x_xz = buffer_pspd[2];

    auto g_x_0_x_yy = buffer_pspd[3];

    auto g_x_0_x_yz = buffer_pspd[4];

    auto g_x_0_x_zz = buffer_pspd[5];

    auto g_x_0_y_xx = buffer_pspd[6];

    auto g_x_0_y_xy = buffer_pspd[7];

    auto g_x_0_y_xz = buffer_pspd[8];

    auto g_x_0_y_yy = buffer_pspd[9];

    auto g_x_0_y_yz = buffer_pspd[10];

    auto g_x_0_y_zz = buffer_pspd[11];

    auto g_x_0_z_xx = buffer_pspd[12];

    auto g_x_0_z_xy = buffer_pspd[13];

    auto g_x_0_z_xz = buffer_pspd[14];

    auto g_x_0_z_yy = buffer_pspd[15];

    auto g_x_0_z_yz = buffer_pspd[16];

    auto g_x_0_z_zz = buffer_pspd[17];

    auto g_y_0_x_xx = buffer_pspd[18];

    auto g_y_0_x_xy = buffer_pspd[19];

    auto g_y_0_x_xz = buffer_pspd[20];

    auto g_y_0_x_yy = buffer_pspd[21];

    auto g_y_0_x_yz = buffer_pspd[22];

    auto g_y_0_x_zz = buffer_pspd[23];

    auto g_y_0_y_xx = buffer_pspd[24];

    auto g_y_0_y_xy = buffer_pspd[25];

    auto g_y_0_y_xz = buffer_pspd[26];

    auto g_y_0_y_yy = buffer_pspd[27];

    auto g_y_0_y_yz = buffer_pspd[28];

    auto g_y_0_y_zz = buffer_pspd[29];

    auto g_y_0_z_xx = buffer_pspd[30];

    auto g_y_0_z_xy = buffer_pspd[31];

    auto g_y_0_z_xz = buffer_pspd[32];

    auto g_y_0_z_yy = buffer_pspd[33];

    auto g_y_0_z_yz = buffer_pspd[34];

    auto g_y_0_z_zz = buffer_pspd[35];

    auto g_z_0_x_xx = buffer_pspd[36];

    auto g_z_0_x_xy = buffer_pspd[37];

    auto g_z_0_x_xz = buffer_pspd[38];

    auto g_z_0_x_yy = buffer_pspd[39];

    auto g_z_0_x_yz = buffer_pspd[40];

    auto g_z_0_x_zz = buffer_pspd[41];

    auto g_z_0_y_xx = buffer_pspd[42];

    auto g_z_0_y_xy = buffer_pspd[43];

    auto g_z_0_y_xz = buffer_pspd[44];

    auto g_z_0_y_yy = buffer_pspd[45];

    auto g_z_0_y_yz = buffer_pspd[46];

    auto g_z_0_y_zz = buffer_pspd[47];

    auto g_z_0_z_xx = buffer_pspd[48];

    auto g_z_0_z_xy = buffer_pspd[49];

    auto g_z_0_z_xz = buffer_pspd[50];

    auto g_z_0_z_yy = buffer_pspd[51];

    auto g_z_0_z_yz = buffer_pspd[52];

    auto g_z_0_z_zz = buffer_pspd[53];

    /// Set up components of integrals buffer : buffer_1000_sspd

    auto g_x_0_0_0_0_0_x_xx = buffer_1000_sspd[0];

    auto g_x_0_0_0_0_0_x_xy = buffer_1000_sspd[1];

    auto g_x_0_0_0_0_0_x_xz = buffer_1000_sspd[2];

    auto g_x_0_0_0_0_0_x_yy = buffer_1000_sspd[3];

    auto g_x_0_0_0_0_0_x_yz = buffer_1000_sspd[4];

    auto g_x_0_0_0_0_0_x_zz = buffer_1000_sspd[5];

    auto g_x_0_0_0_0_0_y_xx = buffer_1000_sspd[6];

    auto g_x_0_0_0_0_0_y_xy = buffer_1000_sspd[7];

    auto g_x_0_0_0_0_0_y_xz = buffer_1000_sspd[8];

    auto g_x_0_0_0_0_0_y_yy = buffer_1000_sspd[9];

    auto g_x_0_0_0_0_0_y_yz = buffer_1000_sspd[10];

    auto g_x_0_0_0_0_0_y_zz = buffer_1000_sspd[11];

    auto g_x_0_0_0_0_0_z_xx = buffer_1000_sspd[12];

    auto g_x_0_0_0_0_0_z_xy = buffer_1000_sspd[13];

    auto g_x_0_0_0_0_0_z_xz = buffer_1000_sspd[14];

    auto g_x_0_0_0_0_0_z_yy = buffer_1000_sspd[15];

    auto g_x_0_0_0_0_0_z_yz = buffer_1000_sspd[16];

    auto g_x_0_0_0_0_0_z_zz = buffer_1000_sspd[17];

    auto g_y_0_0_0_0_0_x_xx = buffer_1000_sspd[18];

    auto g_y_0_0_0_0_0_x_xy = buffer_1000_sspd[19];

    auto g_y_0_0_0_0_0_x_xz = buffer_1000_sspd[20];

    auto g_y_0_0_0_0_0_x_yy = buffer_1000_sspd[21];

    auto g_y_0_0_0_0_0_x_yz = buffer_1000_sspd[22];

    auto g_y_0_0_0_0_0_x_zz = buffer_1000_sspd[23];

    auto g_y_0_0_0_0_0_y_xx = buffer_1000_sspd[24];

    auto g_y_0_0_0_0_0_y_xy = buffer_1000_sspd[25];

    auto g_y_0_0_0_0_0_y_xz = buffer_1000_sspd[26];

    auto g_y_0_0_0_0_0_y_yy = buffer_1000_sspd[27];

    auto g_y_0_0_0_0_0_y_yz = buffer_1000_sspd[28];

    auto g_y_0_0_0_0_0_y_zz = buffer_1000_sspd[29];

    auto g_y_0_0_0_0_0_z_xx = buffer_1000_sspd[30];

    auto g_y_0_0_0_0_0_z_xy = buffer_1000_sspd[31];

    auto g_y_0_0_0_0_0_z_xz = buffer_1000_sspd[32];

    auto g_y_0_0_0_0_0_z_yy = buffer_1000_sspd[33];

    auto g_y_0_0_0_0_0_z_yz = buffer_1000_sspd[34];

    auto g_y_0_0_0_0_0_z_zz = buffer_1000_sspd[35];

    auto g_z_0_0_0_0_0_x_xx = buffer_1000_sspd[36];

    auto g_z_0_0_0_0_0_x_xy = buffer_1000_sspd[37];

    auto g_z_0_0_0_0_0_x_xz = buffer_1000_sspd[38];

    auto g_z_0_0_0_0_0_x_yy = buffer_1000_sspd[39];

    auto g_z_0_0_0_0_0_x_yz = buffer_1000_sspd[40];

    auto g_z_0_0_0_0_0_x_zz = buffer_1000_sspd[41];

    auto g_z_0_0_0_0_0_y_xx = buffer_1000_sspd[42];

    auto g_z_0_0_0_0_0_y_xy = buffer_1000_sspd[43];

    auto g_z_0_0_0_0_0_y_xz = buffer_1000_sspd[44];

    auto g_z_0_0_0_0_0_y_yy = buffer_1000_sspd[45];

    auto g_z_0_0_0_0_0_y_yz = buffer_1000_sspd[46];

    auto g_z_0_0_0_0_0_y_zz = buffer_1000_sspd[47];

    auto g_z_0_0_0_0_0_z_xx = buffer_1000_sspd[48];

    auto g_z_0_0_0_0_0_z_xy = buffer_1000_sspd[49];

    auto g_z_0_0_0_0_0_z_xz = buffer_1000_sspd[50];

    auto g_z_0_0_0_0_0_z_yy = buffer_1000_sspd[51];

    auto g_z_0_0_0_0_0_z_yz = buffer_1000_sspd[52];

    auto g_z_0_0_0_0_0_z_zz = buffer_1000_sspd[53];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_0_0_0_0_0_x_xx, g_x_0_0_0_0_0_x_xy, g_x_0_0_0_0_0_x_xz, g_x_0_0_0_0_0_x_yy, g_x_0_0_0_0_0_x_yz, g_x_0_0_0_0_0_x_zz, g_x_0_x_xx, g_x_0_x_xy, g_x_0_x_xz, g_x_0_x_yy, g_x_0_x_yz, g_x_0_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_0_x_xx[i] = 2.0 * g_x_0_x_xx[i] * a_exp;

        g_x_0_0_0_0_0_x_xy[i] = 2.0 * g_x_0_x_xy[i] * a_exp;

        g_x_0_0_0_0_0_x_xz[i] = 2.0 * g_x_0_x_xz[i] * a_exp;

        g_x_0_0_0_0_0_x_yy[i] = 2.0 * g_x_0_x_yy[i] * a_exp;

        g_x_0_0_0_0_0_x_yz[i] = 2.0 * g_x_0_x_yz[i] * a_exp;

        g_x_0_0_0_0_0_x_zz[i] = 2.0 * g_x_0_x_zz[i] * a_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_0_0_0_0_0_y_xx, g_x_0_0_0_0_0_y_xy, g_x_0_0_0_0_0_y_xz, g_x_0_0_0_0_0_y_yy, g_x_0_0_0_0_0_y_yz, g_x_0_0_0_0_0_y_zz, g_x_0_y_xx, g_x_0_y_xy, g_x_0_y_xz, g_x_0_y_yy, g_x_0_y_yz, g_x_0_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_0_y_xx[i] = 2.0 * g_x_0_y_xx[i] * a_exp;

        g_x_0_0_0_0_0_y_xy[i] = 2.0 * g_x_0_y_xy[i] * a_exp;

        g_x_0_0_0_0_0_y_xz[i] = 2.0 * g_x_0_y_xz[i] * a_exp;

        g_x_0_0_0_0_0_y_yy[i] = 2.0 * g_x_0_y_yy[i] * a_exp;

        g_x_0_0_0_0_0_y_yz[i] = 2.0 * g_x_0_y_yz[i] * a_exp;

        g_x_0_0_0_0_0_y_zz[i] = 2.0 * g_x_0_y_zz[i] * a_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_0_0_0_0_0_z_xx, g_x_0_0_0_0_0_z_xy, g_x_0_0_0_0_0_z_xz, g_x_0_0_0_0_0_z_yy, g_x_0_0_0_0_0_z_yz, g_x_0_0_0_0_0_z_zz, g_x_0_z_xx, g_x_0_z_xy, g_x_0_z_xz, g_x_0_z_yy, g_x_0_z_yz, g_x_0_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_0_z_xx[i] = 2.0 * g_x_0_z_xx[i] * a_exp;

        g_x_0_0_0_0_0_z_xy[i] = 2.0 * g_x_0_z_xy[i] * a_exp;

        g_x_0_0_0_0_0_z_xz[i] = 2.0 * g_x_0_z_xz[i] * a_exp;

        g_x_0_0_0_0_0_z_yy[i] = 2.0 * g_x_0_z_yy[i] * a_exp;

        g_x_0_0_0_0_0_z_yz[i] = 2.0 * g_x_0_z_yz[i] * a_exp;

        g_x_0_0_0_0_0_z_zz[i] = 2.0 * g_x_0_z_zz[i] * a_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_y_0_0_0_0_0_x_xx, g_y_0_0_0_0_0_x_xy, g_y_0_0_0_0_0_x_xz, g_y_0_0_0_0_0_x_yy, g_y_0_0_0_0_0_x_yz, g_y_0_0_0_0_0_x_zz, g_y_0_x_xx, g_y_0_x_xy, g_y_0_x_xz, g_y_0_x_yy, g_y_0_x_yz, g_y_0_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_0_x_xx[i] = 2.0 * g_y_0_x_xx[i] * a_exp;

        g_y_0_0_0_0_0_x_xy[i] = 2.0 * g_y_0_x_xy[i] * a_exp;

        g_y_0_0_0_0_0_x_xz[i] = 2.0 * g_y_0_x_xz[i] * a_exp;

        g_y_0_0_0_0_0_x_yy[i] = 2.0 * g_y_0_x_yy[i] * a_exp;

        g_y_0_0_0_0_0_x_yz[i] = 2.0 * g_y_0_x_yz[i] * a_exp;

        g_y_0_0_0_0_0_x_zz[i] = 2.0 * g_y_0_x_zz[i] * a_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_y_0_0_0_0_0_y_xx, g_y_0_0_0_0_0_y_xy, g_y_0_0_0_0_0_y_xz, g_y_0_0_0_0_0_y_yy, g_y_0_0_0_0_0_y_yz, g_y_0_0_0_0_0_y_zz, g_y_0_y_xx, g_y_0_y_xy, g_y_0_y_xz, g_y_0_y_yy, g_y_0_y_yz, g_y_0_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_0_y_xx[i] = 2.0 * g_y_0_y_xx[i] * a_exp;

        g_y_0_0_0_0_0_y_xy[i] = 2.0 * g_y_0_y_xy[i] * a_exp;

        g_y_0_0_0_0_0_y_xz[i] = 2.0 * g_y_0_y_xz[i] * a_exp;

        g_y_0_0_0_0_0_y_yy[i] = 2.0 * g_y_0_y_yy[i] * a_exp;

        g_y_0_0_0_0_0_y_yz[i] = 2.0 * g_y_0_y_yz[i] * a_exp;

        g_y_0_0_0_0_0_y_zz[i] = 2.0 * g_y_0_y_zz[i] * a_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_y_0_0_0_0_0_z_xx, g_y_0_0_0_0_0_z_xy, g_y_0_0_0_0_0_z_xz, g_y_0_0_0_0_0_z_yy, g_y_0_0_0_0_0_z_yz, g_y_0_0_0_0_0_z_zz, g_y_0_z_xx, g_y_0_z_xy, g_y_0_z_xz, g_y_0_z_yy, g_y_0_z_yz, g_y_0_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_0_z_xx[i] = 2.0 * g_y_0_z_xx[i] * a_exp;

        g_y_0_0_0_0_0_z_xy[i] = 2.0 * g_y_0_z_xy[i] * a_exp;

        g_y_0_0_0_0_0_z_xz[i] = 2.0 * g_y_0_z_xz[i] * a_exp;

        g_y_0_0_0_0_0_z_yy[i] = 2.0 * g_y_0_z_yy[i] * a_exp;

        g_y_0_0_0_0_0_z_yz[i] = 2.0 * g_y_0_z_yz[i] * a_exp;

        g_y_0_0_0_0_0_z_zz[i] = 2.0 * g_y_0_z_zz[i] * a_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_z_0_0_0_0_0_x_xx, g_z_0_0_0_0_0_x_xy, g_z_0_0_0_0_0_x_xz, g_z_0_0_0_0_0_x_yy, g_z_0_0_0_0_0_x_yz, g_z_0_0_0_0_0_x_zz, g_z_0_x_xx, g_z_0_x_xy, g_z_0_x_xz, g_z_0_x_yy, g_z_0_x_yz, g_z_0_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_0_x_xx[i] = 2.0 * g_z_0_x_xx[i] * a_exp;

        g_z_0_0_0_0_0_x_xy[i] = 2.0 * g_z_0_x_xy[i] * a_exp;

        g_z_0_0_0_0_0_x_xz[i] = 2.0 * g_z_0_x_xz[i] * a_exp;

        g_z_0_0_0_0_0_x_yy[i] = 2.0 * g_z_0_x_yy[i] * a_exp;

        g_z_0_0_0_0_0_x_yz[i] = 2.0 * g_z_0_x_yz[i] * a_exp;

        g_z_0_0_0_0_0_x_zz[i] = 2.0 * g_z_0_x_zz[i] * a_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_z_0_0_0_0_0_y_xx, g_z_0_0_0_0_0_y_xy, g_z_0_0_0_0_0_y_xz, g_z_0_0_0_0_0_y_yy, g_z_0_0_0_0_0_y_yz, g_z_0_0_0_0_0_y_zz, g_z_0_y_xx, g_z_0_y_xy, g_z_0_y_xz, g_z_0_y_yy, g_z_0_y_yz, g_z_0_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_0_y_xx[i] = 2.0 * g_z_0_y_xx[i] * a_exp;

        g_z_0_0_0_0_0_y_xy[i] = 2.0 * g_z_0_y_xy[i] * a_exp;

        g_z_0_0_0_0_0_y_xz[i] = 2.0 * g_z_0_y_xz[i] * a_exp;

        g_z_0_0_0_0_0_y_yy[i] = 2.0 * g_z_0_y_yy[i] * a_exp;

        g_z_0_0_0_0_0_y_yz[i] = 2.0 * g_z_0_y_yz[i] * a_exp;

        g_z_0_0_0_0_0_y_zz[i] = 2.0 * g_z_0_y_zz[i] * a_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_z_0_0_0_0_0_z_xx, g_z_0_0_0_0_0_z_xy, g_z_0_0_0_0_0_z_xz, g_z_0_0_0_0_0_z_yy, g_z_0_0_0_0_0_z_yz, g_z_0_0_0_0_0_z_zz, g_z_0_z_xx, g_z_0_z_xy, g_z_0_z_xz, g_z_0_z_yy, g_z_0_z_yz, g_z_0_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_0_z_xx[i] = 2.0 * g_z_0_z_xx[i] * a_exp;

        g_z_0_0_0_0_0_z_xy[i] = 2.0 * g_z_0_z_xy[i] * a_exp;

        g_z_0_0_0_0_0_z_xz[i] = 2.0 * g_z_0_z_xz[i] * a_exp;

        g_z_0_0_0_0_0_z_yy[i] = 2.0 * g_z_0_z_yy[i] * a_exp;

        g_z_0_0_0_0_0_z_yz[i] = 2.0 * g_z_0_z_yz[i] * a_exp;

        g_z_0_0_0_0_0_z_zz[i] = 2.0 * g_z_0_z_zz[i] * a_exp;
    }
}

} // t4c_geom namespace


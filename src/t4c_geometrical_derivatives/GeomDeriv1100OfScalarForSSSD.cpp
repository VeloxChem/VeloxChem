#include "GeomDeriv1100OfScalarForSSSD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1100_sssd_0(CSimdArray<double>& buffer_1100_sssd,
                     const CSimdArray<double>& buffer_ppsd,
                     const double a_exp,
                     const double b_exp) -> void
{
    const auto ndims = buffer_1100_sssd.number_of_columns();

    /// Set up components of auxilary buffer : buffer_ppsd

    auto g_x_x_0_xx = buffer_ppsd[0];

    auto g_x_x_0_xy = buffer_ppsd[1];

    auto g_x_x_0_xz = buffer_ppsd[2];

    auto g_x_x_0_yy = buffer_ppsd[3];

    auto g_x_x_0_yz = buffer_ppsd[4];

    auto g_x_x_0_zz = buffer_ppsd[5];

    auto g_x_y_0_xx = buffer_ppsd[6];

    auto g_x_y_0_xy = buffer_ppsd[7];

    auto g_x_y_0_xz = buffer_ppsd[8];

    auto g_x_y_0_yy = buffer_ppsd[9];

    auto g_x_y_0_yz = buffer_ppsd[10];

    auto g_x_y_0_zz = buffer_ppsd[11];

    auto g_x_z_0_xx = buffer_ppsd[12];

    auto g_x_z_0_xy = buffer_ppsd[13];

    auto g_x_z_0_xz = buffer_ppsd[14];

    auto g_x_z_0_yy = buffer_ppsd[15];

    auto g_x_z_0_yz = buffer_ppsd[16];

    auto g_x_z_0_zz = buffer_ppsd[17];

    auto g_y_x_0_xx = buffer_ppsd[18];

    auto g_y_x_0_xy = buffer_ppsd[19];

    auto g_y_x_0_xz = buffer_ppsd[20];

    auto g_y_x_0_yy = buffer_ppsd[21];

    auto g_y_x_0_yz = buffer_ppsd[22];

    auto g_y_x_0_zz = buffer_ppsd[23];

    auto g_y_y_0_xx = buffer_ppsd[24];

    auto g_y_y_0_xy = buffer_ppsd[25];

    auto g_y_y_0_xz = buffer_ppsd[26];

    auto g_y_y_0_yy = buffer_ppsd[27];

    auto g_y_y_0_yz = buffer_ppsd[28];

    auto g_y_y_0_zz = buffer_ppsd[29];

    auto g_y_z_0_xx = buffer_ppsd[30];

    auto g_y_z_0_xy = buffer_ppsd[31];

    auto g_y_z_0_xz = buffer_ppsd[32];

    auto g_y_z_0_yy = buffer_ppsd[33];

    auto g_y_z_0_yz = buffer_ppsd[34];

    auto g_y_z_0_zz = buffer_ppsd[35];

    auto g_z_x_0_xx = buffer_ppsd[36];

    auto g_z_x_0_xy = buffer_ppsd[37];

    auto g_z_x_0_xz = buffer_ppsd[38];

    auto g_z_x_0_yy = buffer_ppsd[39];

    auto g_z_x_0_yz = buffer_ppsd[40];

    auto g_z_x_0_zz = buffer_ppsd[41];

    auto g_z_y_0_xx = buffer_ppsd[42];

    auto g_z_y_0_xy = buffer_ppsd[43];

    auto g_z_y_0_xz = buffer_ppsd[44];

    auto g_z_y_0_yy = buffer_ppsd[45];

    auto g_z_y_0_yz = buffer_ppsd[46];

    auto g_z_y_0_zz = buffer_ppsd[47];

    auto g_z_z_0_xx = buffer_ppsd[48];

    auto g_z_z_0_xy = buffer_ppsd[49];

    auto g_z_z_0_xz = buffer_ppsd[50];

    auto g_z_z_0_yy = buffer_ppsd[51];

    auto g_z_z_0_yz = buffer_ppsd[52];

    auto g_z_z_0_zz = buffer_ppsd[53];

    /// Set up components of integrals buffer : buffer_1100_sssd

    auto g_x_x_0_0_0_0_0_xx = buffer_1100_sssd[0];

    auto g_x_x_0_0_0_0_0_xy = buffer_1100_sssd[1];

    auto g_x_x_0_0_0_0_0_xz = buffer_1100_sssd[2];

    auto g_x_x_0_0_0_0_0_yy = buffer_1100_sssd[3];

    auto g_x_x_0_0_0_0_0_yz = buffer_1100_sssd[4];

    auto g_x_x_0_0_0_0_0_zz = buffer_1100_sssd[5];

    auto g_x_y_0_0_0_0_0_xx = buffer_1100_sssd[6];

    auto g_x_y_0_0_0_0_0_xy = buffer_1100_sssd[7];

    auto g_x_y_0_0_0_0_0_xz = buffer_1100_sssd[8];

    auto g_x_y_0_0_0_0_0_yy = buffer_1100_sssd[9];

    auto g_x_y_0_0_0_0_0_yz = buffer_1100_sssd[10];

    auto g_x_y_0_0_0_0_0_zz = buffer_1100_sssd[11];

    auto g_x_z_0_0_0_0_0_xx = buffer_1100_sssd[12];

    auto g_x_z_0_0_0_0_0_xy = buffer_1100_sssd[13];

    auto g_x_z_0_0_0_0_0_xz = buffer_1100_sssd[14];

    auto g_x_z_0_0_0_0_0_yy = buffer_1100_sssd[15];

    auto g_x_z_0_0_0_0_0_yz = buffer_1100_sssd[16];

    auto g_x_z_0_0_0_0_0_zz = buffer_1100_sssd[17];

    auto g_y_x_0_0_0_0_0_xx = buffer_1100_sssd[18];

    auto g_y_x_0_0_0_0_0_xy = buffer_1100_sssd[19];

    auto g_y_x_0_0_0_0_0_xz = buffer_1100_sssd[20];

    auto g_y_x_0_0_0_0_0_yy = buffer_1100_sssd[21];

    auto g_y_x_0_0_0_0_0_yz = buffer_1100_sssd[22];

    auto g_y_x_0_0_0_0_0_zz = buffer_1100_sssd[23];

    auto g_y_y_0_0_0_0_0_xx = buffer_1100_sssd[24];

    auto g_y_y_0_0_0_0_0_xy = buffer_1100_sssd[25];

    auto g_y_y_0_0_0_0_0_xz = buffer_1100_sssd[26];

    auto g_y_y_0_0_0_0_0_yy = buffer_1100_sssd[27];

    auto g_y_y_0_0_0_0_0_yz = buffer_1100_sssd[28];

    auto g_y_y_0_0_0_0_0_zz = buffer_1100_sssd[29];

    auto g_y_z_0_0_0_0_0_xx = buffer_1100_sssd[30];

    auto g_y_z_0_0_0_0_0_xy = buffer_1100_sssd[31];

    auto g_y_z_0_0_0_0_0_xz = buffer_1100_sssd[32];

    auto g_y_z_0_0_0_0_0_yy = buffer_1100_sssd[33];

    auto g_y_z_0_0_0_0_0_yz = buffer_1100_sssd[34];

    auto g_y_z_0_0_0_0_0_zz = buffer_1100_sssd[35];

    auto g_z_x_0_0_0_0_0_xx = buffer_1100_sssd[36];

    auto g_z_x_0_0_0_0_0_xy = buffer_1100_sssd[37];

    auto g_z_x_0_0_0_0_0_xz = buffer_1100_sssd[38];

    auto g_z_x_0_0_0_0_0_yy = buffer_1100_sssd[39];

    auto g_z_x_0_0_0_0_0_yz = buffer_1100_sssd[40];

    auto g_z_x_0_0_0_0_0_zz = buffer_1100_sssd[41];

    auto g_z_y_0_0_0_0_0_xx = buffer_1100_sssd[42];

    auto g_z_y_0_0_0_0_0_xy = buffer_1100_sssd[43];

    auto g_z_y_0_0_0_0_0_xz = buffer_1100_sssd[44];

    auto g_z_y_0_0_0_0_0_yy = buffer_1100_sssd[45];

    auto g_z_y_0_0_0_0_0_yz = buffer_1100_sssd[46];

    auto g_z_y_0_0_0_0_0_zz = buffer_1100_sssd[47];

    auto g_z_z_0_0_0_0_0_xx = buffer_1100_sssd[48];

    auto g_z_z_0_0_0_0_0_xy = buffer_1100_sssd[49];

    auto g_z_z_0_0_0_0_0_xz = buffer_1100_sssd[50];

    auto g_z_z_0_0_0_0_0_yy = buffer_1100_sssd[51];

    auto g_z_z_0_0_0_0_0_yz = buffer_1100_sssd[52];

    auto g_z_z_0_0_0_0_0_zz = buffer_1100_sssd[53];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_x_0_0_0_0_0_xx, g_x_x_0_0_0_0_0_xy, g_x_x_0_0_0_0_0_xz, g_x_x_0_0_0_0_0_yy, g_x_x_0_0_0_0_0_yz, g_x_x_0_0_0_0_0_zz, g_x_x_0_xx, g_x_x_0_xy, g_x_x_0_xz, g_x_x_0_yy, g_x_x_0_yz, g_x_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_0_0_xx[i] = 4.0 * g_x_x_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_0_xy[i] = 4.0 * g_x_x_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_0_xz[i] = 4.0 * g_x_x_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_0_yy[i] = 4.0 * g_x_x_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_0_yz[i] = 4.0 * g_x_x_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_0_zz[i] = 4.0 * g_x_x_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_y_0_0_0_0_0_xx, g_x_y_0_0_0_0_0_xy, g_x_y_0_0_0_0_0_xz, g_x_y_0_0_0_0_0_yy, g_x_y_0_0_0_0_0_yz, g_x_y_0_0_0_0_0_zz, g_x_y_0_xx, g_x_y_0_xy, g_x_y_0_xz, g_x_y_0_yy, g_x_y_0_yz, g_x_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_0_0_xx[i] = 4.0 * g_x_y_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_0_xy[i] = 4.0 * g_x_y_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_0_xz[i] = 4.0 * g_x_y_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_0_yy[i] = 4.0 * g_x_y_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_0_yz[i] = 4.0 * g_x_y_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_0_zz[i] = 4.0 * g_x_y_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_z_0_0_0_0_0_xx, g_x_z_0_0_0_0_0_xy, g_x_z_0_0_0_0_0_xz, g_x_z_0_0_0_0_0_yy, g_x_z_0_0_0_0_0_yz, g_x_z_0_0_0_0_0_zz, g_x_z_0_xx, g_x_z_0_xy, g_x_z_0_xz, g_x_z_0_yy, g_x_z_0_yz, g_x_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_0_0_xx[i] = 4.0 * g_x_z_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_0_xy[i] = 4.0 * g_x_z_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_0_xz[i] = 4.0 * g_x_z_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_0_yy[i] = 4.0 * g_x_z_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_0_yz[i] = 4.0 * g_x_z_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_0_zz[i] = 4.0 * g_x_z_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_y_x_0_0_0_0_0_xx, g_y_x_0_0_0_0_0_xy, g_y_x_0_0_0_0_0_xz, g_y_x_0_0_0_0_0_yy, g_y_x_0_0_0_0_0_yz, g_y_x_0_0_0_0_0_zz, g_y_x_0_xx, g_y_x_0_xy, g_y_x_0_xz, g_y_x_0_yy, g_y_x_0_yz, g_y_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_0_0_xx[i] = 4.0 * g_y_x_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_0_xy[i] = 4.0 * g_y_x_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_0_xz[i] = 4.0 * g_y_x_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_0_yy[i] = 4.0 * g_y_x_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_0_yz[i] = 4.0 * g_y_x_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_0_zz[i] = 4.0 * g_y_x_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_y_y_0_0_0_0_0_xx, g_y_y_0_0_0_0_0_xy, g_y_y_0_0_0_0_0_xz, g_y_y_0_0_0_0_0_yy, g_y_y_0_0_0_0_0_yz, g_y_y_0_0_0_0_0_zz, g_y_y_0_xx, g_y_y_0_xy, g_y_y_0_xz, g_y_y_0_yy, g_y_y_0_yz, g_y_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_0_0_xx[i] = 4.0 * g_y_y_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_0_xy[i] = 4.0 * g_y_y_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_0_xz[i] = 4.0 * g_y_y_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_0_yy[i] = 4.0 * g_y_y_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_0_yz[i] = 4.0 * g_y_y_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_0_zz[i] = 4.0 * g_y_y_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_y_z_0_0_0_0_0_xx, g_y_z_0_0_0_0_0_xy, g_y_z_0_0_0_0_0_xz, g_y_z_0_0_0_0_0_yy, g_y_z_0_0_0_0_0_yz, g_y_z_0_0_0_0_0_zz, g_y_z_0_xx, g_y_z_0_xy, g_y_z_0_xz, g_y_z_0_yy, g_y_z_0_yz, g_y_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_0_0_xx[i] = 4.0 * g_y_z_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_0_xy[i] = 4.0 * g_y_z_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_0_xz[i] = 4.0 * g_y_z_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_0_yy[i] = 4.0 * g_y_z_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_0_yz[i] = 4.0 * g_y_z_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_0_zz[i] = 4.0 * g_y_z_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_z_x_0_0_0_0_0_xx, g_z_x_0_0_0_0_0_xy, g_z_x_0_0_0_0_0_xz, g_z_x_0_0_0_0_0_yy, g_z_x_0_0_0_0_0_yz, g_z_x_0_0_0_0_0_zz, g_z_x_0_xx, g_z_x_0_xy, g_z_x_0_xz, g_z_x_0_yy, g_z_x_0_yz, g_z_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_0_0_xx[i] = 4.0 * g_z_x_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_0_xy[i] = 4.0 * g_z_x_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_0_xz[i] = 4.0 * g_z_x_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_0_yy[i] = 4.0 * g_z_x_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_0_yz[i] = 4.0 * g_z_x_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_0_zz[i] = 4.0 * g_z_x_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_z_y_0_0_0_0_0_xx, g_z_y_0_0_0_0_0_xy, g_z_y_0_0_0_0_0_xz, g_z_y_0_0_0_0_0_yy, g_z_y_0_0_0_0_0_yz, g_z_y_0_0_0_0_0_zz, g_z_y_0_xx, g_z_y_0_xy, g_z_y_0_xz, g_z_y_0_yy, g_z_y_0_yz, g_z_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_0_0_xx[i] = 4.0 * g_z_y_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_0_xy[i] = 4.0 * g_z_y_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_0_xz[i] = 4.0 * g_z_y_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_0_yy[i] = 4.0 * g_z_y_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_0_yz[i] = 4.0 * g_z_y_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_0_zz[i] = 4.0 * g_z_y_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_z_z_0_0_0_0_0_xx, g_z_z_0_0_0_0_0_xy, g_z_z_0_0_0_0_0_xz, g_z_z_0_0_0_0_0_yy, g_z_z_0_0_0_0_0_yz, g_z_z_0_0_0_0_0_zz, g_z_z_0_xx, g_z_z_0_xy, g_z_z_0_xz, g_z_z_0_yy, g_z_z_0_yz, g_z_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_0_0_xx[i] = 4.0 * g_z_z_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_0_xy[i] = 4.0 * g_z_z_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_0_xz[i] = 4.0 * g_z_z_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_0_yy[i] = 4.0 * g_z_z_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_0_yz[i] = 4.0 * g_z_z_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_0_zz[i] = 4.0 * g_z_z_0_zz[i] * a_exp * b_exp;
    }
}

} // t4c_geom namespace


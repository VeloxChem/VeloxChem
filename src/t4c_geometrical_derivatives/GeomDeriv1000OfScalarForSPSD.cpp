#include "GeomDeriv1000OfScalarForSPSD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1000_spsd_0(CSimdArray<double>& buffer_1000_spsd,
                     const CSimdArray<double>& buffer_ppsd,
                     const double a_exp) -> void
{
    const auto ndims = buffer_1000_spsd.number_of_columns();

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

    /// Set up components of integrals buffer : buffer_1000_spsd

    auto g_x_0_0_0_0_x_0_xx = buffer_1000_spsd[0];

    auto g_x_0_0_0_0_x_0_xy = buffer_1000_spsd[1];

    auto g_x_0_0_0_0_x_0_xz = buffer_1000_spsd[2];

    auto g_x_0_0_0_0_x_0_yy = buffer_1000_spsd[3];

    auto g_x_0_0_0_0_x_0_yz = buffer_1000_spsd[4];

    auto g_x_0_0_0_0_x_0_zz = buffer_1000_spsd[5];

    auto g_x_0_0_0_0_y_0_xx = buffer_1000_spsd[6];

    auto g_x_0_0_0_0_y_0_xy = buffer_1000_spsd[7];

    auto g_x_0_0_0_0_y_0_xz = buffer_1000_spsd[8];

    auto g_x_0_0_0_0_y_0_yy = buffer_1000_spsd[9];

    auto g_x_0_0_0_0_y_0_yz = buffer_1000_spsd[10];

    auto g_x_0_0_0_0_y_0_zz = buffer_1000_spsd[11];

    auto g_x_0_0_0_0_z_0_xx = buffer_1000_spsd[12];

    auto g_x_0_0_0_0_z_0_xy = buffer_1000_spsd[13];

    auto g_x_0_0_0_0_z_0_xz = buffer_1000_spsd[14];

    auto g_x_0_0_0_0_z_0_yy = buffer_1000_spsd[15];

    auto g_x_0_0_0_0_z_0_yz = buffer_1000_spsd[16];

    auto g_x_0_0_0_0_z_0_zz = buffer_1000_spsd[17];

    auto g_y_0_0_0_0_x_0_xx = buffer_1000_spsd[18];

    auto g_y_0_0_0_0_x_0_xy = buffer_1000_spsd[19];

    auto g_y_0_0_0_0_x_0_xz = buffer_1000_spsd[20];

    auto g_y_0_0_0_0_x_0_yy = buffer_1000_spsd[21];

    auto g_y_0_0_0_0_x_0_yz = buffer_1000_spsd[22];

    auto g_y_0_0_0_0_x_0_zz = buffer_1000_spsd[23];

    auto g_y_0_0_0_0_y_0_xx = buffer_1000_spsd[24];

    auto g_y_0_0_0_0_y_0_xy = buffer_1000_spsd[25];

    auto g_y_0_0_0_0_y_0_xz = buffer_1000_spsd[26];

    auto g_y_0_0_0_0_y_0_yy = buffer_1000_spsd[27];

    auto g_y_0_0_0_0_y_0_yz = buffer_1000_spsd[28];

    auto g_y_0_0_0_0_y_0_zz = buffer_1000_spsd[29];

    auto g_y_0_0_0_0_z_0_xx = buffer_1000_spsd[30];

    auto g_y_0_0_0_0_z_0_xy = buffer_1000_spsd[31];

    auto g_y_0_0_0_0_z_0_xz = buffer_1000_spsd[32];

    auto g_y_0_0_0_0_z_0_yy = buffer_1000_spsd[33];

    auto g_y_0_0_0_0_z_0_yz = buffer_1000_spsd[34];

    auto g_y_0_0_0_0_z_0_zz = buffer_1000_spsd[35];

    auto g_z_0_0_0_0_x_0_xx = buffer_1000_spsd[36];

    auto g_z_0_0_0_0_x_0_xy = buffer_1000_spsd[37];

    auto g_z_0_0_0_0_x_0_xz = buffer_1000_spsd[38];

    auto g_z_0_0_0_0_x_0_yy = buffer_1000_spsd[39];

    auto g_z_0_0_0_0_x_0_yz = buffer_1000_spsd[40];

    auto g_z_0_0_0_0_x_0_zz = buffer_1000_spsd[41];

    auto g_z_0_0_0_0_y_0_xx = buffer_1000_spsd[42];

    auto g_z_0_0_0_0_y_0_xy = buffer_1000_spsd[43];

    auto g_z_0_0_0_0_y_0_xz = buffer_1000_spsd[44];

    auto g_z_0_0_0_0_y_0_yy = buffer_1000_spsd[45];

    auto g_z_0_0_0_0_y_0_yz = buffer_1000_spsd[46];

    auto g_z_0_0_0_0_y_0_zz = buffer_1000_spsd[47];

    auto g_z_0_0_0_0_z_0_xx = buffer_1000_spsd[48];

    auto g_z_0_0_0_0_z_0_xy = buffer_1000_spsd[49];

    auto g_z_0_0_0_0_z_0_xz = buffer_1000_spsd[50];

    auto g_z_0_0_0_0_z_0_yy = buffer_1000_spsd[51];

    auto g_z_0_0_0_0_z_0_yz = buffer_1000_spsd[52];

    auto g_z_0_0_0_0_z_0_zz = buffer_1000_spsd[53];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_0_0_0_0_x_0_xx, g_x_0_0_0_0_x_0_xy, g_x_0_0_0_0_x_0_xz, g_x_0_0_0_0_x_0_yy, g_x_0_0_0_0_x_0_yz, g_x_0_0_0_0_x_0_zz, g_x_x_0_xx, g_x_x_0_xy, g_x_x_0_xz, g_x_x_0_yy, g_x_x_0_yz, g_x_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_x_0_xx[i] = 2.0 * g_x_x_0_xx[i] * a_exp;

        g_x_0_0_0_0_x_0_xy[i] = 2.0 * g_x_x_0_xy[i] * a_exp;

        g_x_0_0_0_0_x_0_xz[i] = 2.0 * g_x_x_0_xz[i] * a_exp;

        g_x_0_0_0_0_x_0_yy[i] = 2.0 * g_x_x_0_yy[i] * a_exp;

        g_x_0_0_0_0_x_0_yz[i] = 2.0 * g_x_x_0_yz[i] * a_exp;

        g_x_0_0_0_0_x_0_zz[i] = 2.0 * g_x_x_0_zz[i] * a_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_0_0_0_0_y_0_xx, g_x_0_0_0_0_y_0_xy, g_x_0_0_0_0_y_0_xz, g_x_0_0_0_0_y_0_yy, g_x_0_0_0_0_y_0_yz, g_x_0_0_0_0_y_0_zz, g_x_y_0_xx, g_x_y_0_xy, g_x_y_0_xz, g_x_y_0_yy, g_x_y_0_yz, g_x_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_y_0_xx[i] = 2.0 * g_x_y_0_xx[i] * a_exp;

        g_x_0_0_0_0_y_0_xy[i] = 2.0 * g_x_y_0_xy[i] * a_exp;

        g_x_0_0_0_0_y_0_xz[i] = 2.0 * g_x_y_0_xz[i] * a_exp;

        g_x_0_0_0_0_y_0_yy[i] = 2.0 * g_x_y_0_yy[i] * a_exp;

        g_x_0_0_0_0_y_0_yz[i] = 2.0 * g_x_y_0_yz[i] * a_exp;

        g_x_0_0_0_0_y_0_zz[i] = 2.0 * g_x_y_0_zz[i] * a_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_0_0_0_0_z_0_xx, g_x_0_0_0_0_z_0_xy, g_x_0_0_0_0_z_0_xz, g_x_0_0_0_0_z_0_yy, g_x_0_0_0_0_z_0_yz, g_x_0_0_0_0_z_0_zz, g_x_z_0_xx, g_x_z_0_xy, g_x_z_0_xz, g_x_z_0_yy, g_x_z_0_yz, g_x_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_z_0_xx[i] = 2.0 * g_x_z_0_xx[i] * a_exp;

        g_x_0_0_0_0_z_0_xy[i] = 2.0 * g_x_z_0_xy[i] * a_exp;

        g_x_0_0_0_0_z_0_xz[i] = 2.0 * g_x_z_0_xz[i] * a_exp;

        g_x_0_0_0_0_z_0_yy[i] = 2.0 * g_x_z_0_yy[i] * a_exp;

        g_x_0_0_0_0_z_0_yz[i] = 2.0 * g_x_z_0_yz[i] * a_exp;

        g_x_0_0_0_0_z_0_zz[i] = 2.0 * g_x_z_0_zz[i] * a_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_y_0_0_0_0_x_0_xx, g_y_0_0_0_0_x_0_xy, g_y_0_0_0_0_x_0_xz, g_y_0_0_0_0_x_0_yy, g_y_0_0_0_0_x_0_yz, g_y_0_0_0_0_x_0_zz, g_y_x_0_xx, g_y_x_0_xy, g_y_x_0_xz, g_y_x_0_yy, g_y_x_0_yz, g_y_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_x_0_xx[i] = 2.0 * g_y_x_0_xx[i] * a_exp;

        g_y_0_0_0_0_x_0_xy[i] = 2.0 * g_y_x_0_xy[i] * a_exp;

        g_y_0_0_0_0_x_0_xz[i] = 2.0 * g_y_x_0_xz[i] * a_exp;

        g_y_0_0_0_0_x_0_yy[i] = 2.0 * g_y_x_0_yy[i] * a_exp;

        g_y_0_0_0_0_x_0_yz[i] = 2.0 * g_y_x_0_yz[i] * a_exp;

        g_y_0_0_0_0_x_0_zz[i] = 2.0 * g_y_x_0_zz[i] * a_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_y_0_0_0_0_y_0_xx, g_y_0_0_0_0_y_0_xy, g_y_0_0_0_0_y_0_xz, g_y_0_0_0_0_y_0_yy, g_y_0_0_0_0_y_0_yz, g_y_0_0_0_0_y_0_zz, g_y_y_0_xx, g_y_y_0_xy, g_y_y_0_xz, g_y_y_0_yy, g_y_y_0_yz, g_y_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_y_0_xx[i] = 2.0 * g_y_y_0_xx[i] * a_exp;

        g_y_0_0_0_0_y_0_xy[i] = 2.0 * g_y_y_0_xy[i] * a_exp;

        g_y_0_0_0_0_y_0_xz[i] = 2.0 * g_y_y_0_xz[i] * a_exp;

        g_y_0_0_0_0_y_0_yy[i] = 2.0 * g_y_y_0_yy[i] * a_exp;

        g_y_0_0_0_0_y_0_yz[i] = 2.0 * g_y_y_0_yz[i] * a_exp;

        g_y_0_0_0_0_y_0_zz[i] = 2.0 * g_y_y_0_zz[i] * a_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_y_0_0_0_0_z_0_xx, g_y_0_0_0_0_z_0_xy, g_y_0_0_0_0_z_0_xz, g_y_0_0_0_0_z_0_yy, g_y_0_0_0_0_z_0_yz, g_y_0_0_0_0_z_0_zz, g_y_z_0_xx, g_y_z_0_xy, g_y_z_0_xz, g_y_z_0_yy, g_y_z_0_yz, g_y_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_z_0_xx[i] = 2.0 * g_y_z_0_xx[i] * a_exp;

        g_y_0_0_0_0_z_0_xy[i] = 2.0 * g_y_z_0_xy[i] * a_exp;

        g_y_0_0_0_0_z_0_xz[i] = 2.0 * g_y_z_0_xz[i] * a_exp;

        g_y_0_0_0_0_z_0_yy[i] = 2.0 * g_y_z_0_yy[i] * a_exp;

        g_y_0_0_0_0_z_0_yz[i] = 2.0 * g_y_z_0_yz[i] * a_exp;

        g_y_0_0_0_0_z_0_zz[i] = 2.0 * g_y_z_0_zz[i] * a_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_z_0_0_0_0_x_0_xx, g_z_0_0_0_0_x_0_xy, g_z_0_0_0_0_x_0_xz, g_z_0_0_0_0_x_0_yy, g_z_0_0_0_0_x_0_yz, g_z_0_0_0_0_x_0_zz, g_z_x_0_xx, g_z_x_0_xy, g_z_x_0_xz, g_z_x_0_yy, g_z_x_0_yz, g_z_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_x_0_xx[i] = 2.0 * g_z_x_0_xx[i] * a_exp;

        g_z_0_0_0_0_x_0_xy[i] = 2.0 * g_z_x_0_xy[i] * a_exp;

        g_z_0_0_0_0_x_0_xz[i] = 2.0 * g_z_x_0_xz[i] * a_exp;

        g_z_0_0_0_0_x_0_yy[i] = 2.0 * g_z_x_0_yy[i] * a_exp;

        g_z_0_0_0_0_x_0_yz[i] = 2.0 * g_z_x_0_yz[i] * a_exp;

        g_z_0_0_0_0_x_0_zz[i] = 2.0 * g_z_x_0_zz[i] * a_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_z_0_0_0_0_y_0_xx, g_z_0_0_0_0_y_0_xy, g_z_0_0_0_0_y_0_xz, g_z_0_0_0_0_y_0_yy, g_z_0_0_0_0_y_0_yz, g_z_0_0_0_0_y_0_zz, g_z_y_0_xx, g_z_y_0_xy, g_z_y_0_xz, g_z_y_0_yy, g_z_y_0_yz, g_z_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_y_0_xx[i] = 2.0 * g_z_y_0_xx[i] * a_exp;

        g_z_0_0_0_0_y_0_xy[i] = 2.0 * g_z_y_0_xy[i] * a_exp;

        g_z_0_0_0_0_y_0_xz[i] = 2.0 * g_z_y_0_xz[i] * a_exp;

        g_z_0_0_0_0_y_0_yy[i] = 2.0 * g_z_y_0_yy[i] * a_exp;

        g_z_0_0_0_0_y_0_yz[i] = 2.0 * g_z_y_0_yz[i] * a_exp;

        g_z_0_0_0_0_y_0_zz[i] = 2.0 * g_z_y_0_zz[i] * a_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_z_0_0_0_0_z_0_xx, g_z_0_0_0_0_z_0_xy, g_z_0_0_0_0_z_0_xz, g_z_0_0_0_0_z_0_yy, g_z_0_0_0_0_z_0_yz, g_z_0_0_0_0_z_0_zz, g_z_z_0_xx, g_z_z_0_xy, g_z_z_0_xz, g_z_z_0_yy, g_z_z_0_yz, g_z_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_z_0_xx[i] = 2.0 * g_z_z_0_xx[i] * a_exp;

        g_z_0_0_0_0_z_0_xy[i] = 2.0 * g_z_z_0_xy[i] * a_exp;

        g_z_0_0_0_0_z_0_xz[i] = 2.0 * g_z_z_0_xz[i] * a_exp;

        g_z_0_0_0_0_z_0_yy[i] = 2.0 * g_z_z_0_yy[i] * a_exp;

        g_z_0_0_0_0_z_0_yz[i] = 2.0 * g_z_z_0_yz[i] * a_exp;

        g_z_0_0_0_0_z_0_zz[i] = 2.0 * g_z_z_0_zz[i] * a_exp;
    }
}

} // t4c_geom namespace


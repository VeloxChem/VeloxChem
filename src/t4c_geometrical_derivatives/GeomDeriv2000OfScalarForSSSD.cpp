#include "GeomDeriv2000OfScalarForSSSD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom2000_sssd_0(CSimdArray<double>& buffer_2000_sssd,
                     const CSimdArray<double>& buffer_sssd,
                     const CSimdArray<double>& buffer_dssd,
                     const double a_exp) -> void
{
    const auto ndims = buffer_2000_sssd.number_of_columns();

    /// Set up components of auxilary buffer : buffer_sssd

    auto g_0_0_0_xx = buffer_sssd[0];

    auto g_0_0_0_xy = buffer_sssd[1];

    auto g_0_0_0_xz = buffer_sssd[2];

    auto g_0_0_0_yy = buffer_sssd[3];

    auto g_0_0_0_yz = buffer_sssd[4];

    auto g_0_0_0_zz = buffer_sssd[5];

    /// Set up components of auxilary buffer : buffer_dssd

    auto g_xx_0_0_xx = buffer_dssd[0];

    auto g_xx_0_0_xy = buffer_dssd[1];

    auto g_xx_0_0_xz = buffer_dssd[2];

    auto g_xx_0_0_yy = buffer_dssd[3];

    auto g_xx_0_0_yz = buffer_dssd[4];

    auto g_xx_0_0_zz = buffer_dssd[5];

    auto g_xy_0_0_xx = buffer_dssd[6];

    auto g_xy_0_0_xy = buffer_dssd[7];

    auto g_xy_0_0_xz = buffer_dssd[8];

    auto g_xy_0_0_yy = buffer_dssd[9];

    auto g_xy_0_0_yz = buffer_dssd[10];

    auto g_xy_0_0_zz = buffer_dssd[11];

    auto g_xz_0_0_xx = buffer_dssd[12];

    auto g_xz_0_0_xy = buffer_dssd[13];

    auto g_xz_0_0_xz = buffer_dssd[14];

    auto g_xz_0_0_yy = buffer_dssd[15];

    auto g_xz_0_0_yz = buffer_dssd[16];

    auto g_xz_0_0_zz = buffer_dssd[17];

    auto g_yy_0_0_xx = buffer_dssd[18];

    auto g_yy_0_0_xy = buffer_dssd[19];

    auto g_yy_0_0_xz = buffer_dssd[20];

    auto g_yy_0_0_yy = buffer_dssd[21];

    auto g_yy_0_0_yz = buffer_dssd[22];

    auto g_yy_0_0_zz = buffer_dssd[23];

    auto g_yz_0_0_xx = buffer_dssd[24];

    auto g_yz_0_0_xy = buffer_dssd[25];

    auto g_yz_0_0_xz = buffer_dssd[26];

    auto g_yz_0_0_yy = buffer_dssd[27];

    auto g_yz_0_0_yz = buffer_dssd[28];

    auto g_yz_0_0_zz = buffer_dssd[29];

    auto g_zz_0_0_xx = buffer_dssd[30];

    auto g_zz_0_0_xy = buffer_dssd[31];

    auto g_zz_0_0_xz = buffer_dssd[32];

    auto g_zz_0_0_yy = buffer_dssd[33];

    auto g_zz_0_0_yz = buffer_dssd[34];

    auto g_zz_0_0_zz = buffer_dssd[35];

    /// Set up components of integrals buffer : buffer_2000_sssd

    auto g_xx_0_0_0_0_0_0_xx = buffer_2000_sssd[0];

    auto g_xx_0_0_0_0_0_0_xy = buffer_2000_sssd[1];

    auto g_xx_0_0_0_0_0_0_xz = buffer_2000_sssd[2];

    auto g_xx_0_0_0_0_0_0_yy = buffer_2000_sssd[3];

    auto g_xx_0_0_0_0_0_0_yz = buffer_2000_sssd[4];

    auto g_xx_0_0_0_0_0_0_zz = buffer_2000_sssd[5];

    auto g_xy_0_0_0_0_0_0_xx = buffer_2000_sssd[6];

    auto g_xy_0_0_0_0_0_0_xy = buffer_2000_sssd[7];

    auto g_xy_0_0_0_0_0_0_xz = buffer_2000_sssd[8];

    auto g_xy_0_0_0_0_0_0_yy = buffer_2000_sssd[9];

    auto g_xy_0_0_0_0_0_0_yz = buffer_2000_sssd[10];

    auto g_xy_0_0_0_0_0_0_zz = buffer_2000_sssd[11];

    auto g_xz_0_0_0_0_0_0_xx = buffer_2000_sssd[12];

    auto g_xz_0_0_0_0_0_0_xy = buffer_2000_sssd[13];

    auto g_xz_0_0_0_0_0_0_xz = buffer_2000_sssd[14];

    auto g_xz_0_0_0_0_0_0_yy = buffer_2000_sssd[15];

    auto g_xz_0_0_0_0_0_0_yz = buffer_2000_sssd[16];

    auto g_xz_0_0_0_0_0_0_zz = buffer_2000_sssd[17];

    auto g_yy_0_0_0_0_0_0_xx = buffer_2000_sssd[18];

    auto g_yy_0_0_0_0_0_0_xy = buffer_2000_sssd[19];

    auto g_yy_0_0_0_0_0_0_xz = buffer_2000_sssd[20];

    auto g_yy_0_0_0_0_0_0_yy = buffer_2000_sssd[21];

    auto g_yy_0_0_0_0_0_0_yz = buffer_2000_sssd[22];

    auto g_yy_0_0_0_0_0_0_zz = buffer_2000_sssd[23];

    auto g_yz_0_0_0_0_0_0_xx = buffer_2000_sssd[24];

    auto g_yz_0_0_0_0_0_0_xy = buffer_2000_sssd[25];

    auto g_yz_0_0_0_0_0_0_xz = buffer_2000_sssd[26];

    auto g_yz_0_0_0_0_0_0_yy = buffer_2000_sssd[27];

    auto g_yz_0_0_0_0_0_0_yz = buffer_2000_sssd[28];

    auto g_yz_0_0_0_0_0_0_zz = buffer_2000_sssd[29];

    auto g_zz_0_0_0_0_0_0_xx = buffer_2000_sssd[30];

    auto g_zz_0_0_0_0_0_0_xy = buffer_2000_sssd[31];

    auto g_zz_0_0_0_0_0_0_xz = buffer_2000_sssd[32];

    auto g_zz_0_0_0_0_0_0_yy = buffer_2000_sssd[33];

    auto g_zz_0_0_0_0_0_0_yz = buffer_2000_sssd[34];

    auto g_zz_0_0_0_0_0_0_zz = buffer_2000_sssd[35];

    // integrals block (0-6)

    #pragma omp simd aligned(g_0_0_0_xx, g_0_0_0_xy, g_0_0_0_xz, g_0_0_0_yy, g_0_0_0_yz, g_0_0_0_zz, g_xx_0_0_0_0_0_0_xx, g_xx_0_0_0_0_0_0_xy, g_xx_0_0_0_0_0_0_xz, g_xx_0_0_0_0_0_0_yy, g_xx_0_0_0_0_0_0_yz, g_xx_0_0_0_0_0_0_zz, g_xx_0_0_xx, g_xx_0_0_xy, g_xx_0_0_xz, g_xx_0_0_yy, g_xx_0_0_yz, g_xx_0_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_0_0_xx[i] = -2.0 * g_0_0_0_xx[i] * a_exp + 4.0 * g_xx_0_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_0_xy[i] = -2.0 * g_0_0_0_xy[i] * a_exp + 4.0 * g_xx_0_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_0_xz[i] = -2.0 * g_0_0_0_xz[i] * a_exp + 4.0 * g_xx_0_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_0_yy[i] = -2.0 * g_0_0_0_yy[i] * a_exp + 4.0 * g_xx_0_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_0_yz[i] = -2.0 * g_0_0_0_yz[i] * a_exp + 4.0 * g_xx_0_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_0_zz[i] = -2.0 * g_0_0_0_zz[i] * a_exp + 4.0 * g_xx_0_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_xy_0_0_0_0_0_0_xx, g_xy_0_0_0_0_0_0_xy, g_xy_0_0_0_0_0_0_xz, g_xy_0_0_0_0_0_0_yy, g_xy_0_0_0_0_0_0_yz, g_xy_0_0_0_0_0_0_zz, g_xy_0_0_xx, g_xy_0_0_xy, g_xy_0_0_xz, g_xy_0_0_yy, g_xy_0_0_yz, g_xy_0_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_0_0_xx[i] = 4.0 * g_xy_0_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_0_xy[i] = 4.0 * g_xy_0_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_0_xz[i] = 4.0 * g_xy_0_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_0_yy[i] = 4.0 * g_xy_0_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_0_yz[i] = 4.0 * g_xy_0_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_0_zz[i] = 4.0 * g_xy_0_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_xz_0_0_0_0_0_0_xx, g_xz_0_0_0_0_0_0_xy, g_xz_0_0_0_0_0_0_xz, g_xz_0_0_0_0_0_0_yy, g_xz_0_0_0_0_0_0_yz, g_xz_0_0_0_0_0_0_zz, g_xz_0_0_xx, g_xz_0_0_xy, g_xz_0_0_xz, g_xz_0_0_yy, g_xz_0_0_yz, g_xz_0_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_0_0_xx[i] = 4.0 * g_xz_0_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_0_xy[i] = 4.0 * g_xz_0_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_0_xz[i] = 4.0 * g_xz_0_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_0_yy[i] = 4.0 * g_xz_0_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_0_yz[i] = 4.0 * g_xz_0_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_0_zz[i] = 4.0 * g_xz_0_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_0_0_0_xx, g_0_0_0_xy, g_0_0_0_xz, g_0_0_0_yy, g_0_0_0_yz, g_0_0_0_zz, g_yy_0_0_0_0_0_0_xx, g_yy_0_0_0_0_0_0_xy, g_yy_0_0_0_0_0_0_xz, g_yy_0_0_0_0_0_0_yy, g_yy_0_0_0_0_0_0_yz, g_yy_0_0_0_0_0_0_zz, g_yy_0_0_xx, g_yy_0_0_xy, g_yy_0_0_xz, g_yy_0_0_yy, g_yy_0_0_yz, g_yy_0_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_0_0_xx[i] = -2.0 * g_0_0_0_xx[i] * a_exp + 4.0 * g_yy_0_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_0_xy[i] = -2.0 * g_0_0_0_xy[i] * a_exp + 4.0 * g_yy_0_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_0_xz[i] = -2.0 * g_0_0_0_xz[i] * a_exp + 4.0 * g_yy_0_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_0_yy[i] = -2.0 * g_0_0_0_yy[i] * a_exp + 4.0 * g_yy_0_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_0_yz[i] = -2.0 * g_0_0_0_yz[i] * a_exp + 4.0 * g_yy_0_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_0_zz[i] = -2.0 * g_0_0_0_zz[i] * a_exp + 4.0 * g_yy_0_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_yz_0_0_0_0_0_0_xx, g_yz_0_0_0_0_0_0_xy, g_yz_0_0_0_0_0_0_xz, g_yz_0_0_0_0_0_0_yy, g_yz_0_0_0_0_0_0_yz, g_yz_0_0_0_0_0_0_zz, g_yz_0_0_xx, g_yz_0_0_xy, g_yz_0_0_xz, g_yz_0_0_yy, g_yz_0_0_yz, g_yz_0_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_0_0_xx[i] = 4.0 * g_yz_0_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_0_xy[i] = 4.0 * g_yz_0_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_0_xz[i] = 4.0 * g_yz_0_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_0_yy[i] = 4.0 * g_yz_0_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_0_yz[i] = 4.0 * g_yz_0_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_0_zz[i] = 4.0 * g_yz_0_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_0_0_0_xx, g_0_0_0_xy, g_0_0_0_xz, g_0_0_0_yy, g_0_0_0_yz, g_0_0_0_zz, g_zz_0_0_0_0_0_0_xx, g_zz_0_0_0_0_0_0_xy, g_zz_0_0_0_0_0_0_xz, g_zz_0_0_0_0_0_0_yy, g_zz_0_0_0_0_0_0_yz, g_zz_0_0_0_0_0_0_zz, g_zz_0_0_xx, g_zz_0_0_xy, g_zz_0_0_xz, g_zz_0_0_yy, g_zz_0_0_yz, g_zz_0_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_0_0_xx[i] = -2.0 * g_0_0_0_xx[i] * a_exp + 4.0 * g_zz_0_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_0_xy[i] = -2.0 * g_0_0_0_xy[i] * a_exp + 4.0 * g_zz_0_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_0_xz[i] = -2.0 * g_0_0_0_xz[i] * a_exp + 4.0 * g_zz_0_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_0_yy[i] = -2.0 * g_0_0_0_yy[i] * a_exp + 4.0 * g_zz_0_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_0_yz[i] = -2.0 * g_0_0_0_yz[i] * a_exp + 4.0 * g_zz_0_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_0_zz[i] = -2.0 * g_0_0_0_zz[i] * a_exp + 4.0 * g_zz_0_0_zz[i] * a_exp * a_exp;
    }
}

} // t4c_geom namespace


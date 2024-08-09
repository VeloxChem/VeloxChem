#include "GeomDeriv2000OfScalarForSDSS.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom2000_sdss_0(CSimdArray<double>& buffer_2000_sdss,
                     const CSimdArray<double>& buffer_sdss,
                     const CSimdArray<double>& buffer_ddss,
                     const double a_exp) -> void
{
    const auto ndims = buffer_2000_sdss.number_of_columns();

    /// Set up components of auxilary buffer : buffer_sdss

    auto g_0_xx_0_0 = buffer_sdss[0];

    auto g_0_xy_0_0 = buffer_sdss[1];

    auto g_0_xz_0_0 = buffer_sdss[2];

    auto g_0_yy_0_0 = buffer_sdss[3];

    auto g_0_yz_0_0 = buffer_sdss[4];

    auto g_0_zz_0_0 = buffer_sdss[5];

    /// Set up components of auxilary buffer : buffer_ddss

    auto g_xx_xx_0_0 = buffer_ddss[0];

    auto g_xx_xy_0_0 = buffer_ddss[1];

    auto g_xx_xz_0_0 = buffer_ddss[2];

    auto g_xx_yy_0_0 = buffer_ddss[3];

    auto g_xx_yz_0_0 = buffer_ddss[4];

    auto g_xx_zz_0_0 = buffer_ddss[5];

    auto g_xy_xx_0_0 = buffer_ddss[6];

    auto g_xy_xy_0_0 = buffer_ddss[7];

    auto g_xy_xz_0_0 = buffer_ddss[8];

    auto g_xy_yy_0_0 = buffer_ddss[9];

    auto g_xy_yz_0_0 = buffer_ddss[10];

    auto g_xy_zz_0_0 = buffer_ddss[11];

    auto g_xz_xx_0_0 = buffer_ddss[12];

    auto g_xz_xy_0_0 = buffer_ddss[13];

    auto g_xz_xz_0_0 = buffer_ddss[14];

    auto g_xz_yy_0_0 = buffer_ddss[15];

    auto g_xz_yz_0_0 = buffer_ddss[16];

    auto g_xz_zz_0_0 = buffer_ddss[17];

    auto g_yy_xx_0_0 = buffer_ddss[18];

    auto g_yy_xy_0_0 = buffer_ddss[19];

    auto g_yy_xz_0_0 = buffer_ddss[20];

    auto g_yy_yy_0_0 = buffer_ddss[21];

    auto g_yy_yz_0_0 = buffer_ddss[22];

    auto g_yy_zz_0_0 = buffer_ddss[23];

    auto g_yz_xx_0_0 = buffer_ddss[24];

    auto g_yz_xy_0_0 = buffer_ddss[25];

    auto g_yz_xz_0_0 = buffer_ddss[26];

    auto g_yz_yy_0_0 = buffer_ddss[27];

    auto g_yz_yz_0_0 = buffer_ddss[28];

    auto g_yz_zz_0_0 = buffer_ddss[29];

    auto g_zz_xx_0_0 = buffer_ddss[30];

    auto g_zz_xy_0_0 = buffer_ddss[31];

    auto g_zz_xz_0_0 = buffer_ddss[32];

    auto g_zz_yy_0_0 = buffer_ddss[33];

    auto g_zz_yz_0_0 = buffer_ddss[34];

    auto g_zz_zz_0_0 = buffer_ddss[35];

    /// Set up components of integrals buffer : buffer_2000_sdss

    auto g_xx_0_0_0_0_xx_0_0 = buffer_2000_sdss[0];

    auto g_xx_0_0_0_0_xy_0_0 = buffer_2000_sdss[1];

    auto g_xx_0_0_0_0_xz_0_0 = buffer_2000_sdss[2];

    auto g_xx_0_0_0_0_yy_0_0 = buffer_2000_sdss[3];

    auto g_xx_0_0_0_0_yz_0_0 = buffer_2000_sdss[4];

    auto g_xx_0_0_0_0_zz_0_0 = buffer_2000_sdss[5];

    auto g_xy_0_0_0_0_xx_0_0 = buffer_2000_sdss[6];

    auto g_xy_0_0_0_0_xy_0_0 = buffer_2000_sdss[7];

    auto g_xy_0_0_0_0_xz_0_0 = buffer_2000_sdss[8];

    auto g_xy_0_0_0_0_yy_0_0 = buffer_2000_sdss[9];

    auto g_xy_0_0_0_0_yz_0_0 = buffer_2000_sdss[10];

    auto g_xy_0_0_0_0_zz_0_0 = buffer_2000_sdss[11];

    auto g_xz_0_0_0_0_xx_0_0 = buffer_2000_sdss[12];

    auto g_xz_0_0_0_0_xy_0_0 = buffer_2000_sdss[13];

    auto g_xz_0_0_0_0_xz_0_0 = buffer_2000_sdss[14];

    auto g_xz_0_0_0_0_yy_0_0 = buffer_2000_sdss[15];

    auto g_xz_0_0_0_0_yz_0_0 = buffer_2000_sdss[16];

    auto g_xz_0_0_0_0_zz_0_0 = buffer_2000_sdss[17];

    auto g_yy_0_0_0_0_xx_0_0 = buffer_2000_sdss[18];

    auto g_yy_0_0_0_0_xy_0_0 = buffer_2000_sdss[19];

    auto g_yy_0_0_0_0_xz_0_0 = buffer_2000_sdss[20];

    auto g_yy_0_0_0_0_yy_0_0 = buffer_2000_sdss[21];

    auto g_yy_0_0_0_0_yz_0_0 = buffer_2000_sdss[22];

    auto g_yy_0_0_0_0_zz_0_0 = buffer_2000_sdss[23];

    auto g_yz_0_0_0_0_xx_0_0 = buffer_2000_sdss[24];

    auto g_yz_0_0_0_0_xy_0_0 = buffer_2000_sdss[25];

    auto g_yz_0_0_0_0_xz_0_0 = buffer_2000_sdss[26];

    auto g_yz_0_0_0_0_yy_0_0 = buffer_2000_sdss[27];

    auto g_yz_0_0_0_0_yz_0_0 = buffer_2000_sdss[28];

    auto g_yz_0_0_0_0_zz_0_0 = buffer_2000_sdss[29];

    auto g_zz_0_0_0_0_xx_0_0 = buffer_2000_sdss[30];

    auto g_zz_0_0_0_0_xy_0_0 = buffer_2000_sdss[31];

    auto g_zz_0_0_0_0_xz_0_0 = buffer_2000_sdss[32];

    auto g_zz_0_0_0_0_yy_0_0 = buffer_2000_sdss[33];

    auto g_zz_0_0_0_0_yz_0_0 = buffer_2000_sdss[34];

    auto g_zz_0_0_0_0_zz_0_0 = buffer_2000_sdss[35];

    // integrals block (0-6)

    #pragma omp simd aligned(g_0_xx_0_0, g_0_xy_0_0, g_0_xz_0_0, g_0_yy_0_0, g_0_yz_0_0, g_0_zz_0_0, g_xx_0_0_0_0_xx_0_0, g_xx_0_0_0_0_xy_0_0, g_xx_0_0_0_0_xz_0_0, g_xx_0_0_0_0_yy_0_0, g_xx_0_0_0_0_yz_0_0, g_xx_0_0_0_0_zz_0_0, g_xx_xx_0_0, g_xx_xy_0_0, g_xx_xz_0_0, g_xx_yy_0_0, g_xx_yz_0_0, g_xx_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_xx_0_0[i] = -2.0 * g_0_xx_0_0[i] * a_exp + 4.0 * g_xx_xx_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_0_0[i] = -2.0 * g_0_xy_0_0[i] * a_exp + 4.0 * g_xx_xy_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_0_0[i] = -2.0 * g_0_xz_0_0[i] * a_exp + 4.0 * g_xx_xz_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_0_0[i] = -2.0 * g_0_yy_0_0[i] * a_exp + 4.0 * g_xx_yy_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_0_0[i] = -2.0 * g_0_yz_0_0[i] * a_exp + 4.0 * g_xx_yz_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_0_0[i] = -2.0 * g_0_zz_0_0[i] * a_exp + 4.0 * g_xx_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_xy_0_0_0_0_xx_0_0, g_xy_0_0_0_0_xy_0_0, g_xy_0_0_0_0_xz_0_0, g_xy_0_0_0_0_yy_0_0, g_xy_0_0_0_0_yz_0_0, g_xy_0_0_0_0_zz_0_0, g_xy_xx_0_0, g_xy_xy_0_0, g_xy_xz_0_0, g_xy_yy_0_0, g_xy_yz_0_0, g_xy_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_xx_0_0[i] = 4.0 * g_xy_xx_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_0_0[i] = 4.0 * g_xy_xy_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_0_0[i] = 4.0 * g_xy_xz_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_0_0[i] = 4.0 * g_xy_yy_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_0_0[i] = 4.0 * g_xy_yz_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_0_0[i] = 4.0 * g_xy_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_xz_0_0_0_0_xx_0_0, g_xz_0_0_0_0_xy_0_0, g_xz_0_0_0_0_xz_0_0, g_xz_0_0_0_0_yy_0_0, g_xz_0_0_0_0_yz_0_0, g_xz_0_0_0_0_zz_0_0, g_xz_xx_0_0, g_xz_xy_0_0, g_xz_xz_0_0, g_xz_yy_0_0, g_xz_yz_0_0, g_xz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_xx_0_0[i] = 4.0 * g_xz_xx_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_0_0[i] = 4.0 * g_xz_xy_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_0_0[i] = 4.0 * g_xz_xz_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_0_0[i] = 4.0 * g_xz_yy_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_0_0[i] = 4.0 * g_xz_yz_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_0_0[i] = 4.0 * g_xz_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_0_xx_0_0, g_0_xy_0_0, g_0_xz_0_0, g_0_yy_0_0, g_0_yz_0_0, g_0_zz_0_0, g_yy_0_0_0_0_xx_0_0, g_yy_0_0_0_0_xy_0_0, g_yy_0_0_0_0_xz_0_0, g_yy_0_0_0_0_yy_0_0, g_yy_0_0_0_0_yz_0_0, g_yy_0_0_0_0_zz_0_0, g_yy_xx_0_0, g_yy_xy_0_0, g_yy_xz_0_0, g_yy_yy_0_0, g_yy_yz_0_0, g_yy_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_xx_0_0[i] = -2.0 * g_0_xx_0_0[i] * a_exp + 4.0 * g_yy_xx_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_0_0[i] = -2.0 * g_0_xy_0_0[i] * a_exp + 4.0 * g_yy_xy_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_0_0[i] = -2.0 * g_0_xz_0_0[i] * a_exp + 4.0 * g_yy_xz_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_0_0[i] = -2.0 * g_0_yy_0_0[i] * a_exp + 4.0 * g_yy_yy_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_0_0[i] = -2.0 * g_0_yz_0_0[i] * a_exp + 4.0 * g_yy_yz_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_0_0[i] = -2.0 * g_0_zz_0_0[i] * a_exp + 4.0 * g_yy_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_yz_0_0_0_0_xx_0_0, g_yz_0_0_0_0_xy_0_0, g_yz_0_0_0_0_xz_0_0, g_yz_0_0_0_0_yy_0_0, g_yz_0_0_0_0_yz_0_0, g_yz_0_0_0_0_zz_0_0, g_yz_xx_0_0, g_yz_xy_0_0, g_yz_xz_0_0, g_yz_yy_0_0, g_yz_yz_0_0, g_yz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_xx_0_0[i] = 4.0 * g_yz_xx_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_0_0[i] = 4.0 * g_yz_xy_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_0_0[i] = 4.0 * g_yz_xz_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_0_0[i] = 4.0 * g_yz_yy_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_0_0[i] = 4.0 * g_yz_yz_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_0_0[i] = 4.0 * g_yz_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_0_xx_0_0, g_0_xy_0_0, g_0_xz_0_0, g_0_yy_0_0, g_0_yz_0_0, g_0_zz_0_0, g_zz_0_0_0_0_xx_0_0, g_zz_0_0_0_0_xy_0_0, g_zz_0_0_0_0_xz_0_0, g_zz_0_0_0_0_yy_0_0, g_zz_0_0_0_0_yz_0_0, g_zz_0_0_0_0_zz_0_0, g_zz_xx_0_0, g_zz_xy_0_0, g_zz_xz_0_0, g_zz_yy_0_0, g_zz_yz_0_0, g_zz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_xx_0_0[i] = -2.0 * g_0_xx_0_0[i] * a_exp + 4.0 * g_zz_xx_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_0_0[i] = -2.0 * g_0_xy_0_0[i] * a_exp + 4.0 * g_zz_xy_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_0_0[i] = -2.0 * g_0_xz_0_0[i] * a_exp + 4.0 * g_zz_xz_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_0_0[i] = -2.0 * g_0_yy_0_0[i] * a_exp + 4.0 * g_zz_yy_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_0_0[i] = -2.0 * g_0_yz_0_0[i] * a_exp + 4.0 * g_zz_yz_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_0_0[i] = -2.0 * g_0_zz_0_0[i] * a_exp + 4.0 * g_zz_zz_0_0[i] * a_exp * a_exp;
    }
}

} // t4c_geom namespace


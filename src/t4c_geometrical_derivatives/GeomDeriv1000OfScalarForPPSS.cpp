#include "GeomDeriv1000OfScalarForPPSS.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1000_ppss_0(CSimdArray<double>& buffer_1000_ppss,
                     const CSimdArray<double>& buffer_spss,
                     const CSimdArray<double>& buffer_dpss,
                     const double a_exp) -> void
{
    const auto ndims = buffer_1000_ppss.number_of_columns();

    /// Set up components of auxilary buffer : buffer_spss

    auto g_0_x_0_0 = buffer_spss[0];

    auto g_0_y_0_0 = buffer_spss[1];

    auto g_0_z_0_0 = buffer_spss[2];

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

    /// Set up components of integrals buffer : buffer_1000_ppss

    auto g_x_0_0_0_x_x_0_0 = buffer_1000_ppss[0];

    auto g_x_0_0_0_x_y_0_0 = buffer_1000_ppss[1];

    auto g_x_0_0_0_x_z_0_0 = buffer_1000_ppss[2];

    auto g_x_0_0_0_y_x_0_0 = buffer_1000_ppss[3];

    auto g_x_0_0_0_y_y_0_0 = buffer_1000_ppss[4];

    auto g_x_0_0_0_y_z_0_0 = buffer_1000_ppss[5];

    auto g_x_0_0_0_z_x_0_0 = buffer_1000_ppss[6];

    auto g_x_0_0_0_z_y_0_0 = buffer_1000_ppss[7];

    auto g_x_0_0_0_z_z_0_0 = buffer_1000_ppss[8];

    auto g_y_0_0_0_x_x_0_0 = buffer_1000_ppss[9];

    auto g_y_0_0_0_x_y_0_0 = buffer_1000_ppss[10];

    auto g_y_0_0_0_x_z_0_0 = buffer_1000_ppss[11];

    auto g_y_0_0_0_y_x_0_0 = buffer_1000_ppss[12];

    auto g_y_0_0_0_y_y_0_0 = buffer_1000_ppss[13];

    auto g_y_0_0_0_y_z_0_0 = buffer_1000_ppss[14];

    auto g_y_0_0_0_z_x_0_0 = buffer_1000_ppss[15];

    auto g_y_0_0_0_z_y_0_0 = buffer_1000_ppss[16];

    auto g_y_0_0_0_z_z_0_0 = buffer_1000_ppss[17];

    auto g_z_0_0_0_x_x_0_0 = buffer_1000_ppss[18];

    auto g_z_0_0_0_x_y_0_0 = buffer_1000_ppss[19];

    auto g_z_0_0_0_x_z_0_0 = buffer_1000_ppss[20];

    auto g_z_0_0_0_y_x_0_0 = buffer_1000_ppss[21];

    auto g_z_0_0_0_y_y_0_0 = buffer_1000_ppss[22];

    auto g_z_0_0_0_y_z_0_0 = buffer_1000_ppss[23];

    auto g_z_0_0_0_z_x_0_0 = buffer_1000_ppss[24];

    auto g_z_0_0_0_z_y_0_0 = buffer_1000_ppss[25];

    auto g_z_0_0_0_z_z_0_0 = buffer_1000_ppss[26];

    // integrals block (0-3)

    #pragma omp simd aligned(g_0_x_0_0, g_0_y_0_0, g_0_z_0_0, g_x_0_0_0_x_x_0_0, g_x_0_0_0_x_y_0_0, g_x_0_0_0_x_z_0_0, g_xx_x_0_0, g_xx_y_0_0, g_xx_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_x_0_0[i] = -g_0_x_0_0[i] + 2.0 * g_xx_x_0_0[i] * a_exp;

        g_x_0_0_0_x_y_0_0[i] = -g_0_y_0_0[i] + 2.0 * g_xx_y_0_0[i] * a_exp;

        g_x_0_0_0_x_z_0_0[i] = -g_0_z_0_0[i] + 2.0 * g_xx_z_0_0[i] * a_exp;
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_x_0_0_0_y_x_0_0, g_x_0_0_0_y_y_0_0, g_x_0_0_0_y_z_0_0, g_xy_x_0_0, g_xy_y_0_0, g_xy_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_x_0_0[i] = 2.0 * g_xy_x_0_0[i] * a_exp;

        g_x_0_0_0_y_y_0_0[i] = 2.0 * g_xy_y_0_0[i] * a_exp;

        g_x_0_0_0_y_z_0_0[i] = 2.0 * g_xy_z_0_0[i] * a_exp;
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_x_0_0_0_z_x_0_0, g_x_0_0_0_z_y_0_0, g_x_0_0_0_z_z_0_0, g_xz_x_0_0, g_xz_y_0_0, g_xz_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_x_0_0[i] = 2.0 * g_xz_x_0_0[i] * a_exp;

        g_x_0_0_0_z_y_0_0[i] = 2.0 * g_xz_y_0_0[i] * a_exp;

        g_x_0_0_0_z_z_0_0[i] = 2.0 * g_xz_z_0_0[i] * a_exp;
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_xy_x_0_0, g_xy_y_0_0, g_xy_z_0_0, g_y_0_0_0_x_x_0_0, g_y_0_0_0_x_y_0_0, g_y_0_0_0_x_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_x_0_0[i] = 2.0 * g_xy_x_0_0[i] * a_exp;

        g_y_0_0_0_x_y_0_0[i] = 2.0 * g_xy_y_0_0[i] * a_exp;

        g_y_0_0_0_x_z_0_0[i] = 2.0 * g_xy_z_0_0[i] * a_exp;
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_0_x_0_0, g_0_y_0_0, g_0_z_0_0, g_y_0_0_0_y_x_0_0, g_y_0_0_0_y_y_0_0, g_y_0_0_0_y_z_0_0, g_yy_x_0_0, g_yy_y_0_0, g_yy_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_x_0_0[i] = -g_0_x_0_0[i] + 2.0 * g_yy_x_0_0[i] * a_exp;

        g_y_0_0_0_y_y_0_0[i] = -g_0_y_0_0[i] + 2.0 * g_yy_y_0_0[i] * a_exp;

        g_y_0_0_0_y_z_0_0[i] = -g_0_z_0_0[i] + 2.0 * g_yy_z_0_0[i] * a_exp;
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_y_0_0_0_z_x_0_0, g_y_0_0_0_z_y_0_0, g_y_0_0_0_z_z_0_0, g_yz_x_0_0, g_yz_y_0_0, g_yz_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_x_0_0[i] = 2.0 * g_yz_x_0_0[i] * a_exp;

        g_y_0_0_0_z_y_0_0[i] = 2.0 * g_yz_y_0_0[i] * a_exp;

        g_y_0_0_0_z_z_0_0[i] = 2.0 * g_yz_z_0_0[i] * a_exp;
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_xz_x_0_0, g_xz_y_0_0, g_xz_z_0_0, g_z_0_0_0_x_x_0_0, g_z_0_0_0_x_y_0_0, g_z_0_0_0_x_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_x_0_0[i] = 2.0 * g_xz_x_0_0[i] * a_exp;

        g_z_0_0_0_x_y_0_0[i] = 2.0 * g_xz_y_0_0[i] * a_exp;

        g_z_0_0_0_x_z_0_0[i] = 2.0 * g_xz_z_0_0[i] * a_exp;
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_yz_x_0_0, g_yz_y_0_0, g_yz_z_0_0, g_z_0_0_0_y_x_0_0, g_z_0_0_0_y_y_0_0, g_z_0_0_0_y_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_x_0_0[i] = 2.0 * g_yz_x_0_0[i] * a_exp;

        g_z_0_0_0_y_y_0_0[i] = 2.0 * g_yz_y_0_0[i] * a_exp;

        g_z_0_0_0_y_z_0_0[i] = 2.0 * g_yz_z_0_0[i] * a_exp;
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_0_x_0_0, g_0_y_0_0, g_0_z_0_0, g_z_0_0_0_z_x_0_0, g_z_0_0_0_z_y_0_0, g_z_0_0_0_z_z_0_0, g_zz_x_0_0, g_zz_y_0_0, g_zz_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_x_0_0[i] = -g_0_x_0_0[i] + 2.0 * g_zz_x_0_0[i] * a_exp;

        g_z_0_0_0_z_y_0_0[i] = -g_0_y_0_0[i] + 2.0 * g_zz_y_0_0[i] * a_exp;

        g_z_0_0_0_z_z_0_0[i] = -g_0_z_0_0[i] + 2.0 * g_zz_z_0_0[i] * a_exp;
    }
}

} // t4c_geom namespace


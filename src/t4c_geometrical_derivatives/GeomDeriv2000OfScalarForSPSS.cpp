#include "GeomDeriv2000OfScalarForSPSS.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom2000_spss_0(CSimdArray<double>& buffer_2000_spss,
                     const CSimdArray<double>& buffer_spss,
                     const CSimdArray<double>& buffer_dpss,
                     const double a_exp) -> void
{
    const auto ndims = buffer_2000_spss.number_of_columns();

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

    /// Set up components of integrals buffer : buffer_2000_spss

    auto g_xx_0_0_0_0_x_0_0 = buffer_2000_spss[0];

    auto g_xx_0_0_0_0_y_0_0 = buffer_2000_spss[1];

    auto g_xx_0_0_0_0_z_0_0 = buffer_2000_spss[2];

    auto g_xy_0_0_0_0_x_0_0 = buffer_2000_spss[3];

    auto g_xy_0_0_0_0_y_0_0 = buffer_2000_spss[4];

    auto g_xy_0_0_0_0_z_0_0 = buffer_2000_spss[5];

    auto g_xz_0_0_0_0_x_0_0 = buffer_2000_spss[6];

    auto g_xz_0_0_0_0_y_0_0 = buffer_2000_spss[7];

    auto g_xz_0_0_0_0_z_0_0 = buffer_2000_spss[8];

    auto g_yy_0_0_0_0_x_0_0 = buffer_2000_spss[9];

    auto g_yy_0_0_0_0_y_0_0 = buffer_2000_spss[10];

    auto g_yy_0_0_0_0_z_0_0 = buffer_2000_spss[11];

    auto g_yz_0_0_0_0_x_0_0 = buffer_2000_spss[12];

    auto g_yz_0_0_0_0_y_0_0 = buffer_2000_spss[13];

    auto g_yz_0_0_0_0_z_0_0 = buffer_2000_spss[14];

    auto g_zz_0_0_0_0_x_0_0 = buffer_2000_spss[15];

    auto g_zz_0_0_0_0_y_0_0 = buffer_2000_spss[16];

    auto g_zz_0_0_0_0_z_0_0 = buffer_2000_spss[17];

    // integrals block (0-3)

    #pragma omp simd aligned(g_0_x_0_0, g_0_y_0_0, g_0_z_0_0, g_xx_0_0_0_0_x_0_0, g_xx_0_0_0_0_y_0_0, g_xx_0_0_0_0_z_0_0, g_xx_x_0_0, g_xx_y_0_0, g_xx_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_x_0_0[i] = -2.0 * g_0_x_0_0[i] * a_exp + 4.0 * g_xx_x_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_0_0[i] = -2.0 * g_0_y_0_0[i] * a_exp + 4.0 * g_xx_y_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_0_0[i] = -2.0 * g_0_z_0_0[i] * a_exp + 4.0 * g_xx_z_0_0[i] * a_exp * a_exp;
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_xy_0_0_0_0_x_0_0, g_xy_0_0_0_0_y_0_0, g_xy_0_0_0_0_z_0_0, g_xy_x_0_0, g_xy_y_0_0, g_xy_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_x_0_0[i] = 4.0 * g_xy_x_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_0_0[i] = 4.0 * g_xy_y_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_0_0[i] = 4.0 * g_xy_z_0_0[i] * a_exp * a_exp;
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_xz_0_0_0_0_x_0_0, g_xz_0_0_0_0_y_0_0, g_xz_0_0_0_0_z_0_0, g_xz_x_0_0, g_xz_y_0_0, g_xz_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_x_0_0[i] = 4.0 * g_xz_x_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_0_0[i] = 4.0 * g_xz_y_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_0_0[i] = 4.0 * g_xz_z_0_0[i] * a_exp * a_exp;
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_0_x_0_0, g_0_y_0_0, g_0_z_0_0, g_yy_0_0_0_0_x_0_0, g_yy_0_0_0_0_y_0_0, g_yy_0_0_0_0_z_0_0, g_yy_x_0_0, g_yy_y_0_0, g_yy_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_x_0_0[i] = -2.0 * g_0_x_0_0[i] * a_exp + 4.0 * g_yy_x_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_0_0[i] = -2.0 * g_0_y_0_0[i] * a_exp + 4.0 * g_yy_y_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_0_0[i] = -2.0 * g_0_z_0_0[i] * a_exp + 4.0 * g_yy_z_0_0[i] * a_exp * a_exp;
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_yz_0_0_0_0_x_0_0, g_yz_0_0_0_0_y_0_0, g_yz_0_0_0_0_z_0_0, g_yz_x_0_0, g_yz_y_0_0, g_yz_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_x_0_0[i] = 4.0 * g_yz_x_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_0_0[i] = 4.0 * g_yz_y_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_0_0[i] = 4.0 * g_yz_z_0_0[i] * a_exp * a_exp;
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_0_x_0_0, g_0_y_0_0, g_0_z_0_0, g_zz_0_0_0_0_x_0_0, g_zz_0_0_0_0_y_0_0, g_zz_0_0_0_0_z_0_0, g_zz_x_0_0, g_zz_y_0_0, g_zz_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_x_0_0[i] = -2.0 * g_0_x_0_0[i] * a_exp + 4.0 * g_zz_x_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_0_0[i] = -2.0 * g_0_y_0_0[i] * a_exp + 4.0 * g_zz_y_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_0_0[i] = -2.0 * g_0_z_0_0[i] * a_exp + 4.0 * g_zz_z_0_0[i] * a_exp * a_exp;
    }
}

} // t4c_geom namespace


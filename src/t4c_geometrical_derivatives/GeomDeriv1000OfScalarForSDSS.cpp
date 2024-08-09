#include "GeomDeriv1000OfScalarForSDSS.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1000_sdss_0(CSimdArray<double>& buffer_1000_sdss,
                     const CSimdArray<double>& buffer_pdss,
                     const double a_exp) -> void
{
    const auto ndims = buffer_1000_sdss.number_of_columns();

    /// Set up components of auxilary buffer : buffer_pdss

    auto g_x_xx_0_0 = buffer_pdss[0];

    auto g_x_xy_0_0 = buffer_pdss[1];

    auto g_x_xz_0_0 = buffer_pdss[2];

    auto g_x_yy_0_0 = buffer_pdss[3];

    auto g_x_yz_0_0 = buffer_pdss[4];

    auto g_x_zz_0_0 = buffer_pdss[5];

    auto g_y_xx_0_0 = buffer_pdss[6];

    auto g_y_xy_0_0 = buffer_pdss[7];

    auto g_y_xz_0_0 = buffer_pdss[8];

    auto g_y_yy_0_0 = buffer_pdss[9];

    auto g_y_yz_0_0 = buffer_pdss[10];

    auto g_y_zz_0_0 = buffer_pdss[11];

    auto g_z_xx_0_0 = buffer_pdss[12];

    auto g_z_xy_0_0 = buffer_pdss[13];

    auto g_z_xz_0_0 = buffer_pdss[14];

    auto g_z_yy_0_0 = buffer_pdss[15];

    auto g_z_yz_0_0 = buffer_pdss[16];

    auto g_z_zz_0_0 = buffer_pdss[17];

    /// Set up components of integrals buffer : buffer_1000_sdss

    auto g_x_0_0_0_0_xx_0_0 = buffer_1000_sdss[0];

    auto g_x_0_0_0_0_xy_0_0 = buffer_1000_sdss[1];

    auto g_x_0_0_0_0_xz_0_0 = buffer_1000_sdss[2];

    auto g_x_0_0_0_0_yy_0_0 = buffer_1000_sdss[3];

    auto g_x_0_0_0_0_yz_0_0 = buffer_1000_sdss[4];

    auto g_x_0_0_0_0_zz_0_0 = buffer_1000_sdss[5];

    auto g_y_0_0_0_0_xx_0_0 = buffer_1000_sdss[6];

    auto g_y_0_0_0_0_xy_0_0 = buffer_1000_sdss[7];

    auto g_y_0_0_0_0_xz_0_0 = buffer_1000_sdss[8];

    auto g_y_0_0_0_0_yy_0_0 = buffer_1000_sdss[9];

    auto g_y_0_0_0_0_yz_0_0 = buffer_1000_sdss[10];

    auto g_y_0_0_0_0_zz_0_0 = buffer_1000_sdss[11];

    auto g_z_0_0_0_0_xx_0_0 = buffer_1000_sdss[12];

    auto g_z_0_0_0_0_xy_0_0 = buffer_1000_sdss[13];

    auto g_z_0_0_0_0_xz_0_0 = buffer_1000_sdss[14];

    auto g_z_0_0_0_0_yy_0_0 = buffer_1000_sdss[15];

    auto g_z_0_0_0_0_yz_0_0 = buffer_1000_sdss[16];

    auto g_z_0_0_0_0_zz_0_0 = buffer_1000_sdss[17];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_0_0_0_0_xx_0_0, g_x_0_0_0_0_xy_0_0, g_x_0_0_0_0_xz_0_0, g_x_0_0_0_0_yy_0_0, g_x_0_0_0_0_yz_0_0, g_x_0_0_0_0_zz_0_0, g_x_xx_0_0, g_x_xy_0_0, g_x_xz_0_0, g_x_yy_0_0, g_x_yz_0_0, g_x_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xx_0_0[i] = 2.0 * g_x_xx_0_0[i] * a_exp;

        g_x_0_0_0_0_xy_0_0[i] = 2.0 * g_x_xy_0_0[i] * a_exp;

        g_x_0_0_0_0_xz_0_0[i] = 2.0 * g_x_xz_0_0[i] * a_exp;

        g_x_0_0_0_0_yy_0_0[i] = 2.0 * g_x_yy_0_0[i] * a_exp;

        g_x_0_0_0_0_yz_0_0[i] = 2.0 * g_x_yz_0_0[i] * a_exp;

        g_x_0_0_0_0_zz_0_0[i] = 2.0 * g_x_zz_0_0[i] * a_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_y_0_0_0_0_xx_0_0, g_y_0_0_0_0_xy_0_0, g_y_0_0_0_0_xz_0_0, g_y_0_0_0_0_yy_0_0, g_y_0_0_0_0_yz_0_0, g_y_0_0_0_0_zz_0_0, g_y_xx_0_0, g_y_xy_0_0, g_y_xz_0_0, g_y_yy_0_0, g_y_yz_0_0, g_y_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xx_0_0[i] = 2.0 * g_y_xx_0_0[i] * a_exp;

        g_y_0_0_0_0_xy_0_0[i] = 2.0 * g_y_xy_0_0[i] * a_exp;

        g_y_0_0_0_0_xz_0_0[i] = 2.0 * g_y_xz_0_0[i] * a_exp;

        g_y_0_0_0_0_yy_0_0[i] = 2.0 * g_y_yy_0_0[i] * a_exp;

        g_y_0_0_0_0_yz_0_0[i] = 2.0 * g_y_yz_0_0[i] * a_exp;

        g_y_0_0_0_0_zz_0_0[i] = 2.0 * g_y_zz_0_0[i] * a_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_z_0_0_0_0_xx_0_0, g_z_0_0_0_0_xy_0_0, g_z_0_0_0_0_xz_0_0, g_z_0_0_0_0_yy_0_0, g_z_0_0_0_0_yz_0_0, g_z_0_0_0_0_zz_0_0, g_z_xx_0_0, g_z_xy_0_0, g_z_xz_0_0, g_z_yy_0_0, g_z_yz_0_0, g_z_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xx_0_0[i] = 2.0 * g_z_xx_0_0[i] * a_exp;

        g_z_0_0_0_0_xy_0_0[i] = 2.0 * g_z_xy_0_0[i] * a_exp;

        g_z_0_0_0_0_xz_0_0[i] = 2.0 * g_z_xz_0_0[i] * a_exp;

        g_z_0_0_0_0_yy_0_0[i] = 2.0 * g_z_yy_0_0[i] * a_exp;

        g_z_0_0_0_0_yz_0_0[i] = 2.0 * g_z_yz_0_0[i] * a_exp;

        g_z_0_0_0_0_zz_0_0[i] = 2.0 * g_z_zz_0_0[i] * a_exp;
    }
}

} // t4c_geom namespace


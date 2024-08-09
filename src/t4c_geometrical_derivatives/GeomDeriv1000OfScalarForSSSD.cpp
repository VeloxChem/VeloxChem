#include "GeomDeriv1000OfScalarForSSSD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1000_sssd_0(CSimdArray<double>& buffer_1000_sssd,
                     const CSimdArray<double>& buffer_pssd,
                     const double a_exp) -> void
{
    const auto ndims = buffer_1000_sssd.number_of_columns();

    /// Set up components of auxilary buffer : buffer_pssd

    auto g_x_0_0_xx = buffer_pssd[0];

    auto g_x_0_0_xy = buffer_pssd[1];

    auto g_x_0_0_xz = buffer_pssd[2];

    auto g_x_0_0_yy = buffer_pssd[3];

    auto g_x_0_0_yz = buffer_pssd[4];

    auto g_x_0_0_zz = buffer_pssd[5];

    auto g_y_0_0_xx = buffer_pssd[6];

    auto g_y_0_0_xy = buffer_pssd[7];

    auto g_y_0_0_xz = buffer_pssd[8];

    auto g_y_0_0_yy = buffer_pssd[9];

    auto g_y_0_0_yz = buffer_pssd[10];

    auto g_y_0_0_zz = buffer_pssd[11];

    auto g_z_0_0_xx = buffer_pssd[12];

    auto g_z_0_0_xy = buffer_pssd[13];

    auto g_z_0_0_xz = buffer_pssd[14];

    auto g_z_0_0_yy = buffer_pssd[15];

    auto g_z_0_0_yz = buffer_pssd[16];

    auto g_z_0_0_zz = buffer_pssd[17];

    /// Set up components of integrals buffer : buffer_1000_sssd

    auto g_x_0_0_0_0_0_0_xx = buffer_1000_sssd[0];

    auto g_x_0_0_0_0_0_0_xy = buffer_1000_sssd[1];

    auto g_x_0_0_0_0_0_0_xz = buffer_1000_sssd[2];

    auto g_x_0_0_0_0_0_0_yy = buffer_1000_sssd[3];

    auto g_x_0_0_0_0_0_0_yz = buffer_1000_sssd[4];

    auto g_x_0_0_0_0_0_0_zz = buffer_1000_sssd[5];

    auto g_y_0_0_0_0_0_0_xx = buffer_1000_sssd[6];

    auto g_y_0_0_0_0_0_0_xy = buffer_1000_sssd[7];

    auto g_y_0_0_0_0_0_0_xz = buffer_1000_sssd[8];

    auto g_y_0_0_0_0_0_0_yy = buffer_1000_sssd[9];

    auto g_y_0_0_0_0_0_0_yz = buffer_1000_sssd[10];

    auto g_y_0_0_0_0_0_0_zz = buffer_1000_sssd[11];

    auto g_z_0_0_0_0_0_0_xx = buffer_1000_sssd[12];

    auto g_z_0_0_0_0_0_0_xy = buffer_1000_sssd[13];

    auto g_z_0_0_0_0_0_0_xz = buffer_1000_sssd[14];

    auto g_z_0_0_0_0_0_0_yy = buffer_1000_sssd[15];

    auto g_z_0_0_0_0_0_0_yz = buffer_1000_sssd[16];

    auto g_z_0_0_0_0_0_0_zz = buffer_1000_sssd[17];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_0_0_0_0_0_0_xx, g_x_0_0_0_0_0_0_xy, g_x_0_0_0_0_0_0_xz, g_x_0_0_0_0_0_0_yy, g_x_0_0_0_0_0_0_yz, g_x_0_0_0_0_0_0_zz, g_x_0_0_xx, g_x_0_0_xy, g_x_0_0_xz, g_x_0_0_yy, g_x_0_0_yz, g_x_0_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_0_0_xx[i] = 2.0 * g_x_0_0_xx[i] * a_exp;

        g_x_0_0_0_0_0_0_xy[i] = 2.0 * g_x_0_0_xy[i] * a_exp;

        g_x_0_0_0_0_0_0_xz[i] = 2.0 * g_x_0_0_xz[i] * a_exp;

        g_x_0_0_0_0_0_0_yy[i] = 2.0 * g_x_0_0_yy[i] * a_exp;

        g_x_0_0_0_0_0_0_yz[i] = 2.0 * g_x_0_0_yz[i] * a_exp;

        g_x_0_0_0_0_0_0_zz[i] = 2.0 * g_x_0_0_zz[i] * a_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_y_0_0_0_0_0_0_xx, g_y_0_0_0_0_0_0_xy, g_y_0_0_0_0_0_0_xz, g_y_0_0_0_0_0_0_yy, g_y_0_0_0_0_0_0_yz, g_y_0_0_0_0_0_0_zz, g_y_0_0_xx, g_y_0_0_xy, g_y_0_0_xz, g_y_0_0_yy, g_y_0_0_yz, g_y_0_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_0_0_xx[i] = 2.0 * g_y_0_0_xx[i] * a_exp;

        g_y_0_0_0_0_0_0_xy[i] = 2.0 * g_y_0_0_xy[i] * a_exp;

        g_y_0_0_0_0_0_0_xz[i] = 2.0 * g_y_0_0_xz[i] * a_exp;

        g_y_0_0_0_0_0_0_yy[i] = 2.0 * g_y_0_0_yy[i] * a_exp;

        g_y_0_0_0_0_0_0_yz[i] = 2.0 * g_y_0_0_yz[i] * a_exp;

        g_y_0_0_0_0_0_0_zz[i] = 2.0 * g_y_0_0_zz[i] * a_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_z_0_0_0_0_0_0_xx, g_z_0_0_0_0_0_0_xy, g_z_0_0_0_0_0_0_xz, g_z_0_0_0_0_0_0_yy, g_z_0_0_0_0_0_0_yz, g_z_0_0_0_0_0_0_zz, g_z_0_0_xx, g_z_0_0_xy, g_z_0_0_xz, g_z_0_0_yy, g_z_0_0_yz, g_z_0_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_0_0_xx[i] = 2.0 * g_z_0_0_xx[i] * a_exp;

        g_z_0_0_0_0_0_0_xy[i] = 2.0 * g_z_0_0_xy[i] * a_exp;

        g_z_0_0_0_0_0_0_xz[i] = 2.0 * g_z_0_0_xz[i] * a_exp;

        g_z_0_0_0_0_0_0_yy[i] = 2.0 * g_z_0_0_yy[i] * a_exp;

        g_z_0_0_0_0_0_0_yz[i] = 2.0 * g_z_0_0_yz[i] * a_exp;

        g_z_0_0_0_0_0_0_zz[i] = 2.0 * g_z_0_0_zz[i] * a_exp;
    }
}

} // t4c_geom namespace


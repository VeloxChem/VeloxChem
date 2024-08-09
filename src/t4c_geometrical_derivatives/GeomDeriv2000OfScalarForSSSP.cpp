#include "GeomDeriv2000OfScalarForSSSP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom2000_sssp_0(CSimdArray<double>& buffer_2000_sssp,
                     const CSimdArray<double>& buffer_sssp,
                     const CSimdArray<double>& buffer_dssp,
                     const double a_exp) -> void
{
    const auto ndims = buffer_2000_sssp.number_of_columns();

    /// Set up components of auxilary buffer : buffer_sssp

    auto g_0_0_0_x = buffer_sssp[0];

    auto g_0_0_0_y = buffer_sssp[1];

    auto g_0_0_0_z = buffer_sssp[2];

    /// Set up components of auxilary buffer : buffer_dssp

    auto g_xx_0_0_x = buffer_dssp[0];

    auto g_xx_0_0_y = buffer_dssp[1];

    auto g_xx_0_0_z = buffer_dssp[2];

    auto g_xy_0_0_x = buffer_dssp[3];

    auto g_xy_0_0_y = buffer_dssp[4];

    auto g_xy_0_0_z = buffer_dssp[5];

    auto g_xz_0_0_x = buffer_dssp[6];

    auto g_xz_0_0_y = buffer_dssp[7];

    auto g_xz_0_0_z = buffer_dssp[8];

    auto g_yy_0_0_x = buffer_dssp[9];

    auto g_yy_0_0_y = buffer_dssp[10];

    auto g_yy_0_0_z = buffer_dssp[11];

    auto g_yz_0_0_x = buffer_dssp[12];

    auto g_yz_0_0_y = buffer_dssp[13];

    auto g_yz_0_0_z = buffer_dssp[14];

    auto g_zz_0_0_x = buffer_dssp[15];

    auto g_zz_0_0_y = buffer_dssp[16];

    auto g_zz_0_0_z = buffer_dssp[17];

    /// Set up components of integrals buffer : buffer_2000_sssp

    auto g_xx_0_0_0_0_0_0_x = buffer_2000_sssp[0];

    auto g_xx_0_0_0_0_0_0_y = buffer_2000_sssp[1];

    auto g_xx_0_0_0_0_0_0_z = buffer_2000_sssp[2];

    auto g_xy_0_0_0_0_0_0_x = buffer_2000_sssp[3];

    auto g_xy_0_0_0_0_0_0_y = buffer_2000_sssp[4];

    auto g_xy_0_0_0_0_0_0_z = buffer_2000_sssp[5];

    auto g_xz_0_0_0_0_0_0_x = buffer_2000_sssp[6];

    auto g_xz_0_0_0_0_0_0_y = buffer_2000_sssp[7];

    auto g_xz_0_0_0_0_0_0_z = buffer_2000_sssp[8];

    auto g_yy_0_0_0_0_0_0_x = buffer_2000_sssp[9];

    auto g_yy_0_0_0_0_0_0_y = buffer_2000_sssp[10];

    auto g_yy_0_0_0_0_0_0_z = buffer_2000_sssp[11];

    auto g_yz_0_0_0_0_0_0_x = buffer_2000_sssp[12];

    auto g_yz_0_0_0_0_0_0_y = buffer_2000_sssp[13];

    auto g_yz_0_0_0_0_0_0_z = buffer_2000_sssp[14];

    auto g_zz_0_0_0_0_0_0_x = buffer_2000_sssp[15];

    auto g_zz_0_0_0_0_0_0_y = buffer_2000_sssp[16];

    auto g_zz_0_0_0_0_0_0_z = buffer_2000_sssp[17];

    // integrals block (0-3)

    #pragma omp simd aligned(g_0_0_0_x, g_0_0_0_y, g_0_0_0_z, g_xx_0_0_0_0_0_0_x, g_xx_0_0_0_0_0_0_y, g_xx_0_0_0_0_0_0_z, g_xx_0_0_x, g_xx_0_0_y, g_xx_0_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_0_0_x[i] = -2.0 * g_0_0_0_x[i] * a_exp + 4.0 * g_xx_0_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_0_y[i] = -2.0 * g_0_0_0_y[i] * a_exp + 4.0 * g_xx_0_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_0_z[i] = -2.0 * g_0_0_0_z[i] * a_exp + 4.0 * g_xx_0_0_z[i] * a_exp * a_exp;
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_xy_0_0_0_0_0_0_x, g_xy_0_0_0_0_0_0_y, g_xy_0_0_0_0_0_0_z, g_xy_0_0_x, g_xy_0_0_y, g_xy_0_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_0_0_x[i] = 4.0 * g_xy_0_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_0_y[i] = 4.0 * g_xy_0_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_0_z[i] = 4.0 * g_xy_0_0_z[i] * a_exp * a_exp;
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_xz_0_0_0_0_0_0_x, g_xz_0_0_0_0_0_0_y, g_xz_0_0_0_0_0_0_z, g_xz_0_0_x, g_xz_0_0_y, g_xz_0_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_0_0_x[i] = 4.0 * g_xz_0_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_0_y[i] = 4.0 * g_xz_0_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_0_z[i] = 4.0 * g_xz_0_0_z[i] * a_exp * a_exp;
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_0_0_0_x, g_0_0_0_y, g_0_0_0_z, g_yy_0_0_0_0_0_0_x, g_yy_0_0_0_0_0_0_y, g_yy_0_0_0_0_0_0_z, g_yy_0_0_x, g_yy_0_0_y, g_yy_0_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_0_0_x[i] = -2.0 * g_0_0_0_x[i] * a_exp + 4.0 * g_yy_0_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_0_y[i] = -2.0 * g_0_0_0_y[i] * a_exp + 4.0 * g_yy_0_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_0_z[i] = -2.0 * g_0_0_0_z[i] * a_exp + 4.0 * g_yy_0_0_z[i] * a_exp * a_exp;
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_yz_0_0_0_0_0_0_x, g_yz_0_0_0_0_0_0_y, g_yz_0_0_0_0_0_0_z, g_yz_0_0_x, g_yz_0_0_y, g_yz_0_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_0_0_x[i] = 4.0 * g_yz_0_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_0_y[i] = 4.0 * g_yz_0_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_0_z[i] = 4.0 * g_yz_0_0_z[i] * a_exp * a_exp;
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_0_0_0_x, g_0_0_0_y, g_0_0_0_z, g_zz_0_0_0_0_0_0_x, g_zz_0_0_0_0_0_0_y, g_zz_0_0_0_0_0_0_z, g_zz_0_0_x, g_zz_0_0_y, g_zz_0_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_0_0_x[i] = -2.0 * g_0_0_0_x[i] * a_exp + 4.0 * g_zz_0_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_0_y[i] = -2.0 * g_0_0_0_y[i] * a_exp + 4.0 * g_zz_0_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_0_z[i] = -2.0 * g_0_0_0_z[i] * a_exp + 4.0 * g_zz_0_0_z[i] * a_exp * a_exp;
    }
}

} // t4c_geom namespace


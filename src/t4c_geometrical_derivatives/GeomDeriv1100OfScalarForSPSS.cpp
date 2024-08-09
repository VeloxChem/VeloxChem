#include "GeomDeriv1100OfScalarForSPSS.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1100_spss_0(CSimdArray<double>& buffer_1100_spss,
                     const CSimdArray<double>& buffer_psss,
                     const CSimdArray<double>& buffer_pdss,
                     const double a_exp,
                     const double b_exp) -> void
{
    const auto ndims = buffer_1100_spss.number_of_columns();

    /// Set up components of auxilary buffer : buffer_psss

    auto g_x_0_0_0 = buffer_psss[0];

    auto g_y_0_0_0 = buffer_psss[1];

    auto g_z_0_0_0 = buffer_psss[2];

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

    /// Set up components of integrals buffer : buffer_1100_spss

    auto g_x_x_0_0_0_x_0_0 = buffer_1100_spss[0];

    auto g_x_x_0_0_0_y_0_0 = buffer_1100_spss[1];

    auto g_x_x_0_0_0_z_0_0 = buffer_1100_spss[2];

    auto g_x_y_0_0_0_x_0_0 = buffer_1100_spss[3];

    auto g_x_y_0_0_0_y_0_0 = buffer_1100_spss[4];

    auto g_x_y_0_0_0_z_0_0 = buffer_1100_spss[5];

    auto g_x_z_0_0_0_x_0_0 = buffer_1100_spss[6];

    auto g_x_z_0_0_0_y_0_0 = buffer_1100_spss[7];

    auto g_x_z_0_0_0_z_0_0 = buffer_1100_spss[8];

    auto g_y_x_0_0_0_x_0_0 = buffer_1100_spss[9];

    auto g_y_x_0_0_0_y_0_0 = buffer_1100_spss[10];

    auto g_y_x_0_0_0_z_0_0 = buffer_1100_spss[11];

    auto g_y_y_0_0_0_x_0_0 = buffer_1100_spss[12];

    auto g_y_y_0_0_0_y_0_0 = buffer_1100_spss[13];

    auto g_y_y_0_0_0_z_0_0 = buffer_1100_spss[14];

    auto g_y_z_0_0_0_x_0_0 = buffer_1100_spss[15];

    auto g_y_z_0_0_0_y_0_0 = buffer_1100_spss[16];

    auto g_y_z_0_0_0_z_0_0 = buffer_1100_spss[17];

    auto g_z_x_0_0_0_x_0_0 = buffer_1100_spss[18];

    auto g_z_x_0_0_0_y_0_0 = buffer_1100_spss[19];

    auto g_z_x_0_0_0_z_0_0 = buffer_1100_spss[20];

    auto g_z_y_0_0_0_x_0_0 = buffer_1100_spss[21];

    auto g_z_y_0_0_0_y_0_0 = buffer_1100_spss[22];

    auto g_z_y_0_0_0_z_0_0 = buffer_1100_spss[23];

    auto g_z_z_0_0_0_x_0_0 = buffer_1100_spss[24];

    auto g_z_z_0_0_0_y_0_0 = buffer_1100_spss[25];

    auto g_z_z_0_0_0_z_0_0 = buffer_1100_spss[26];

    // integrals block (0-3)

    #pragma omp simd aligned(g_x_0_0_0, g_x_x_0_0_0_x_0_0, g_x_x_0_0_0_y_0_0, g_x_x_0_0_0_z_0_0, g_x_xx_0_0, g_x_xy_0_0, g_x_xz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_x_0_0[i] = -2.0 * g_x_0_0_0[i] * a_exp + 4.0 * g_x_xx_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_0_0[i] = 4.0 * g_x_xy_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_0_0[i] = 4.0 * g_x_xz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_x_0_0_0, g_x_xy_0_0, g_x_y_0_0_0_x_0_0, g_x_y_0_0_0_y_0_0, g_x_y_0_0_0_z_0_0, g_x_yy_0_0, g_x_yz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_x_0_0[i] = 4.0 * g_x_xy_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_0_0[i] = -2.0 * g_x_0_0_0[i] * a_exp + 4.0 * g_x_yy_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_0_0[i] = 4.0 * g_x_yz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_x_0_0_0, g_x_xz_0_0, g_x_yz_0_0, g_x_z_0_0_0_x_0_0, g_x_z_0_0_0_y_0_0, g_x_z_0_0_0_z_0_0, g_x_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_x_0_0[i] = 4.0 * g_x_xz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_0_0[i] = 4.0 * g_x_yz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_0_0[i] = -2.0 * g_x_0_0_0[i] * a_exp + 4.0 * g_x_zz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_y_0_0_0, g_y_x_0_0_0_x_0_0, g_y_x_0_0_0_y_0_0, g_y_x_0_0_0_z_0_0, g_y_xx_0_0, g_y_xy_0_0, g_y_xz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_x_0_0[i] = -2.0 * g_y_0_0_0[i] * a_exp + 4.0 * g_y_xx_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_0_0[i] = 4.0 * g_y_xy_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_0_0[i] = 4.0 * g_y_xz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_y_0_0_0, g_y_xy_0_0, g_y_y_0_0_0_x_0_0, g_y_y_0_0_0_y_0_0, g_y_y_0_0_0_z_0_0, g_y_yy_0_0, g_y_yz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_x_0_0[i] = 4.0 * g_y_xy_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_0_0[i] = -2.0 * g_y_0_0_0[i] * a_exp + 4.0 * g_y_yy_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_0_0[i] = 4.0 * g_y_yz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_y_0_0_0, g_y_xz_0_0, g_y_yz_0_0, g_y_z_0_0_0_x_0_0, g_y_z_0_0_0_y_0_0, g_y_z_0_0_0_z_0_0, g_y_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_x_0_0[i] = 4.0 * g_y_xz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_0_0[i] = 4.0 * g_y_yz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_0_0[i] = -2.0 * g_y_0_0_0[i] * a_exp + 4.0 * g_y_zz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_z_0_0_0, g_z_x_0_0_0_x_0_0, g_z_x_0_0_0_y_0_0, g_z_x_0_0_0_z_0_0, g_z_xx_0_0, g_z_xy_0_0, g_z_xz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_x_0_0[i] = -2.0 * g_z_0_0_0[i] * a_exp + 4.0 * g_z_xx_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_0_0[i] = 4.0 * g_z_xy_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_0_0[i] = 4.0 * g_z_xz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_z_0_0_0, g_z_xy_0_0, g_z_y_0_0_0_x_0_0, g_z_y_0_0_0_y_0_0, g_z_y_0_0_0_z_0_0, g_z_yy_0_0, g_z_yz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_x_0_0[i] = 4.0 * g_z_xy_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_0_0[i] = -2.0 * g_z_0_0_0[i] * a_exp + 4.0 * g_z_yy_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_0_0[i] = 4.0 * g_z_yz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_z_0_0_0, g_z_xz_0_0, g_z_yz_0_0, g_z_z_0_0_0_x_0_0, g_z_z_0_0_0_y_0_0, g_z_z_0_0_0_z_0_0, g_z_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_x_0_0[i] = 4.0 * g_z_xz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_0_0[i] = 4.0 * g_z_yz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_0_0[i] = -2.0 * g_z_0_0_0[i] * a_exp + 4.0 * g_z_zz_0_0[i] * a_exp * b_exp;
    }
}

} // t4c_geom namespace


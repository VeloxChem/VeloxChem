#include "GeomDeriv1000OfScalarForPDSS.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1000_pdss_0(CSimdArray<double>& buffer_1000_pdss,
                     const CSimdArray<double>& buffer_sdss,
                     const CSimdArray<double>& buffer_ddss,
                     const double a_exp) -> void
{
    const auto ndims = buffer_1000_pdss.number_of_columns();

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

    /// Set up components of integrals buffer : buffer_1000_pdss

    auto g_x_0_0_0_x_xx_0_0 = buffer_1000_pdss[0];

    auto g_x_0_0_0_x_xy_0_0 = buffer_1000_pdss[1];

    auto g_x_0_0_0_x_xz_0_0 = buffer_1000_pdss[2];

    auto g_x_0_0_0_x_yy_0_0 = buffer_1000_pdss[3];

    auto g_x_0_0_0_x_yz_0_0 = buffer_1000_pdss[4];

    auto g_x_0_0_0_x_zz_0_0 = buffer_1000_pdss[5];

    auto g_x_0_0_0_y_xx_0_0 = buffer_1000_pdss[6];

    auto g_x_0_0_0_y_xy_0_0 = buffer_1000_pdss[7];

    auto g_x_0_0_0_y_xz_0_0 = buffer_1000_pdss[8];

    auto g_x_0_0_0_y_yy_0_0 = buffer_1000_pdss[9];

    auto g_x_0_0_0_y_yz_0_0 = buffer_1000_pdss[10];

    auto g_x_0_0_0_y_zz_0_0 = buffer_1000_pdss[11];

    auto g_x_0_0_0_z_xx_0_0 = buffer_1000_pdss[12];

    auto g_x_0_0_0_z_xy_0_0 = buffer_1000_pdss[13];

    auto g_x_0_0_0_z_xz_0_0 = buffer_1000_pdss[14];

    auto g_x_0_0_0_z_yy_0_0 = buffer_1000_pdss[15];

    auto g_x_0_0_0_z_yz_0_0 = buffer_1000_pdss[16];

    auto g_x_0_0_0_z_zz_0_0 = buffer_1000_pdss[17];

    auto g_y_0_0_0_x_xx_0_0 = buffer_1000_pdss[18];

    auto g_y_0_0_0_x_xy_0_0 = buffer_1000_pdss[19];

    auto g_y_0_0_0_x_xz_0_0 = buffer_1000_pdss[20];

    auto g_y_0_0_0_x_yy_0_0 = buffer_1000_pdss[21];

    auto g_y_0_0_0_x_yz_0_0 = buffer_1000_pdss[22];

    auto g_y_0_0_0_x_zz_0_0 = buffer_1000_pdss[23];

    auto g_y_0_0_0_y_xx_0_0 = buffer_1000_pdss[24];

    auto g_y_0_0_0_y_xy_0_0 = buffer_1000_pdss[25];

    auto g_y_0_0_0_y_xz_0_0 = buffer_1000_pdss[26];

    auto g_y_0_0_0_y_yy_0_0 = buffer_1000_pdss[27];

    auto g_y_0_0_0_y_yz_0_0 = buffer_1000_pdss[28];

    auto g_y_0_0_0_y_zz_0_0 = buffer_1000_pdss[29];

    auto g_y_0_0_0_z_xx_0_0 = buffer_1000_pdss[30];

    auto g_y_0_0_0_z_xy_0_0 = buffer_1000_pdss[31];

    auto g_y_0_0_0_z_xz_0_0 = buffer_1000_pdss[32];

    auto g_y_0_0_0_z_yy_0_0 = buffer_1000_pdss[33];

    auto g_y_0_0_0_z_yz_0_0 = buffer_1000_pdss[34];

    auto g_y_0_0_0_z_zz_0_0 = buffer_1000_pdss[35];

    auto g_z_0_0_0_x_xx_0_0 = buffer_1000_pdss[36];

    auto g_z_0_0_0_x_xy_0_0 = buffer_1000_pdss[37];

    auto g_z_0_0_0_x_xz_0_0 = buffer_1000_pdss[38];

    auto g_z_0_0_0_x_yy_0_0 = buffer_1000_pdss[39];

    auto g_z_0_0_0_x_yz_0_0 = buffer_1000_pdss[40];

    auto g_z_0_0_0_x_zz_0_0 = buffer_1000_pdss[41];

    auto g_z_0_0_0_y_xx_0_0 = buffer_1000_pdss[42];

    auto g_z_0_0_0_y_xy_0_0 = buffer_1000_pdss[43];

    auto g_z_0_0_0_y_xz_0_0 = buffer_1000_pdss[44];

    auto g_z_0_0_0_y_yy_0_0 = buffer_1000_pdss[45];

    auto g_z_0_0_0_y_yz_0_0 = buffer_1000_pdss[46];

    auto g_z_0_0_0_y_zz_0_0 = buffer_1000_pdss[47];

    auto g_z_0_0_0_z_xx_0_0 = buffer_1000_pdss[48];

    auto g_z_0_0_0_z_xy_0_0 = buffer_1000_pdss[49];

    auto g_z_0_0_0_z_xz_0_0 = buffer_1000_pdss[50];

    auto g_z_0_0_0_z_yy_0_0 = buffer_1000_pdss[51];

    auto g_z_0_0_0_z_yz_0_0 = buffer_1000_pdss[52];

    auto g_z_0_0_0_z_zz_0_0 = buffer_1000_pdss[53];

    // integrals block (0-6)

    #pragma omp simd aligned(g_0_xx_0_0, g_0_xy_0_0, g_0_xz_0_0, g_0_yy_0_0, g_0_yz_0_0, g_0_zz_0_0, g_x_0_0_0_x_xx_0_0, g_x_0_0_0_x_xy_0_0, g_x_0_0_0_x_xz_0_0, g_x_0_0_0_x_yy_0_0, g_x_0_0_0_x_yz_0_0, g_x_0_0_0_x_zz_0_0, g_xx_xx_0_0, g_xx_xy_0_0, g_xx_xz_0_0, g_xx_yy_0_0, g_xx_yz_0_0, g_xx_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xx_0_0[i] = -g_0_xx_0_0[i] + 2.0 * g_xx_xx_0_0[i] * a_exp;

        g_x_0_0_0_x_xy_0_0[i] = -g_0_xy_0_0[i] + 2.0 * g_xx_xy_0_0[i] * a_exp;

        g_x_0_0_0_x_xz_0_0[i] = -g_0_xz_0_0[i] + 2.0 * g_xx_xz_0_0[i] * a_exp;

        g_x_0_0_0_x_yy_0_0[i] = -g_0_yy_0_0[i] + 2.0 * g_xx_yy_0_0[i] * a_exp;

        g_x_0_0_0_x_yz_0_0[i] = -g_0_yz_0_0[i] + 2.0 * g_xx_yz_0_0[i] * a_exp;

        g_x_0_0_0_x_zz_0_0[i] = -g_0_zz_0_0[i] + 2.0 * g_xx_zz_0_0[i] * a_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_0_0_0_y_xx_0_0, g_x_0_0_0_y_xy_0_0, g_x_0_0_0_y_xz_0_0, g_x_0_0_0_y_yy_0_0, g_x_0_0_0_y_yz_0_0, g_x_0_0_0_y_zz_0_0, g_xy_xx_0_0, g_xy_xy_0_0, g_xy_xz_0_0, g_xy_yy_0_0, g_xy_yz_0_0, g_xy_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xx_0_0[i] = 2.0 * g_xy_xx_0_0[i] * a_exp;

        g_x_0_0_0_y_xy_0_0[i] = 2.0 * g_xy_xy_0_0[i] * a_exp;

        g_x_0_0_0_y_xz_0_0[i] = 2.0 * g_xy_xz_0_0[i] * a_exp;

        g_x_0_0_0_y_yy_0_0[i] = 2.0 * g_xy_yy_0_0[i] * a_exp;

        g_x_0_0_0_y_yz_0_0[i] = 2.0 * g_xy_yz_0_0[i] * a_exp;

        g_x_0_0_0_y_zz_0_0[i] = 2.0 * g_xy_zz_0_0[i] * a_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_0_0_0_z_xx_0_0, g_x_0_0_0_z_xy_0_0, g_x_0_0_0_z_xz_0_0, g_x_0_0_0_z_yy_0_0, g_x_0_0_0_z_yz_0_0, g_x_0_0_0_z_zz_0_0, g_xz_xx_0_0, g_xz_xy_0_0, g_xz_xz_0_0, g_xz_yy_0_0, g_xz_yz_0_0, g_xz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xx_0_0[i] = 2.0 * g_xz_xx_0_0[i] * a_exp;

        g_x_0_0_0_z_xy_0_0[i] = 2.0 * g_xz_xy_0_0[i] * a_exp;

        g_x_0_0_0_z_xz_0_0[i] = 2.0 * g_xz_xz_0_0[i] * a_exp;

        g_x_0_0_0_z_yy_0_0[i] = 2.0 * g_xz_yy_0_0[i] * a_exp;

        g_x_0_0_0_z_yz_0_0[i] = 2.0 * g_xz_yz_0_0[i] * a_exp;

        g_x_0_0_0_z_zz_0_0[i] = 2.0 * g_xz_zz_0_0[i] * a_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_xy_xx_0_0, g_xy_xy_0_0, g_xy_xz_0_0, g_xy_yy_0_0, g_xy_yz_0_0, g_xy_zz_0_0, g_y_0_0_0_x_xx_0_0, g_y_0_0_0_x_xy_0_0, g_y_0_0_0_x_xz_0_0, g_y_0_0_0_x_yy_0_0, g_y_0_0_0_x_yz_0_0, g_y_0_0_0_x_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xx_0_0[i] = 2.0 * g_xy_xx_0_0[i] * a_exp;

        g_y_0_0_0_x_xy_0_0[i] = 2.0 * g_xy_xy_0_0[i] * a_exp;

        g_y_0_0_0_x_xz_0_0[i] = 2.0 * g_xy_xz_0_0[i] * a_exp;

        g_y_0_0_0_x_yy_0_0[i] = 2.0 * g_xy_yy_0_0[i] * a_exp;

        g_y_0_0_0_x_yz_0_0[i] = 2.0 * g_xy_yz_0_0[i] * a_exp;

        g_y_0_0_0_x_zz_0_0[i] = 2.0 * g_xy_zz_0_0[i] * a_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_0_xx_0_0, g_0_xy_0_0, g_0_xz_0_0, g_0_yy_0_0, g_0_yz_0_0, g_0_zz_0_0, g_y_0_0_0_y_xx_0_0, g_y_0_0_0_y_xy_0_0, g_y_0_0_0_y_xz_0_0, g_y_0_0_0_y_yy_0_0, g_y_0_0_0_y_yz_0_0, g_y_0_0_0_y_zz_0_0, g_yy_xx_0_0, g_yy_xy_0_0, g_yy_xz_0_0, g_yy_yy_0_0, g_yy_yz_0_0, g_yy_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xx_0_0[i] = -g_0_xx_0_0[i] + 2.0 * g_yy_xx_0_0[i] * a_exp;

        g_y_0_0_0_y_xy_0_0[i] = -g_0_xy_0_0[i] + 2.0 * g_yy_xy_0_0[i] * a_exp;

        g_y_0_0_0_y_xz_0_0[i] = -g_0_xz_0_0[i] + 2.0 * g_yy_xz_0_0[i] * a_exp;

        g_y_0_0_0_y_yy_0_0[i] = -g_0_yy_0_0[i] + 2.0 * g_yy_yy_0_0[i] * a_exp;

        g_y_0_0_0_y_yz_0_0[i] = -g_0_yz_0_0[i] + 2.0 * g_yy_yz_0_0[i] * a_exp;

        g_y_0_0_0_y_zz_0_0[i] = -g_0_zz_0_0[i] + 2.0 * g_yy_zz_0_0[i] * a_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_y_0_0_0_z_xx_0_0, g_y_0_0_0_z_xy_0_0, g_y_0_0_0_z_xz_0_0, g_y_0_0_0_z_yy_0_0, g_y_0_0_0_z_yz_0_0, g_y_0_0_0_z_zz_0_0, g_yz_xx_0_0, g_yz_xy_0_0, g_yz_xz_0_0, g_yz_yy_0_0, g_yz_yz_0_0, g_yz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xx_0_0[i] = 2.0 * g_yz_xx_0_0[i] * a_exp;

        g_y_0_0_0_z_xy_0_0[i] = 2.0 * g_yz_xy_0_0[i] * a_exp;

        g_y_0_0_0_z_xz_0_0[i] = 2.0 * g_yz_xz_0_0[i] * a_exp;

        g_y_0_0_0_z_yy_0_0[i] = 2.0 * g_yz_yy_0_0[i] * a_exp;

        g_y_0_0_0_z_yz_0_0[i] = 2.0 * g_yz_yz_0_0[i] * a_exp;

        g_y_0_0_0_z_zz_0_0[i] = 2.0 * g_yz_zz_0_0[i] * a_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_xz_xx_0_0, g_xz_xy_0_0, g_xz_xz_0_0, g_xz_yy_0_0, g_xz_yz_0_0, g_xz_zz_0_0, g_z_0_0_0_x_xx_0_0, g_z_0_0_0_x_xy_0_0, g_z_0_0_0_x_xz_0_0, g_z_0_0_0_x_yy_0_0, g_z_0_0_0_x_yz_0_0, g_z_0_0_0_x_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xx_0_0[i] = 2.0 * g_xz_xx_0_0[i] * a_exp;

        g_z_0_0_0_x_xy_0_0[i] = 2.0 * g_xz_xy_0_0[i] * a_exp;

        g_z_0_0_0_x_xz_0_0[i] = 2.0 * g_xz_xz_0_0[i] * a_exp;

        g_z_0_0_0_x_yy_0_0[i] = 2.0 * g_xz_yy_0_0[i] * a_exp;

        g_z_0_0_0_x_yz_0_0[i] = 2.0 * g_xz_yz_0_0[i] * a_exp;

        g_z_0_0_0_x_zz_0_0[i] = 2.0 * g_xz_zz_0_0[i] * a_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_yz_xx_0_0, g_yz_xy_0_0, g_yz_xz_0_0, g_yz_yy_0_0, g_yz_yz_0_0, g_yz_zz_0_0, g_z_0_0_0_y_xx_0_0, g_z_0_0_0_y_xy_0_0, g_z_0_0_0_y_xz_0_0, g_z_0_0_0_y_yy_0_0, g_z_0_0_0_y_yz_0_0, g_z_0_0_0_y_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xx_0_0[i] = 2.0 * g_yz_xx_0_0[i] * a_exp;

        g_z_0_0_0_y_xy_0_0[i] = 2.0 * g_yz_xy_0_0[i] * a_exp;

        g_z_0_0_0_y_xz_0_0[i] = 2.0 * g_yz_xz_0_0[i] * a_exp;

        g_z_0_0_0_y_yy_0_0[i] = 2.0 * g_yz_yy_0_0[i] * a_exp;

        g_z_0_0_0_y_yz_0_0[i] = 2.0 * g_yz_yz_0_0[i] * a_exp;

        g_z_0_0_0_y_zz_0_0[i] = 2.0 * g_yz_zz_0_0[i] * a_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_0_xx_0_0, g_0_xy_0_0, g_0_xz_0_0, g_0_yy_0_0, g_0_yz_0_0, g_0_zz_0_0, g_z_0_0_0_z_xx_0_0, g_z_0_0_0_z_xy_0_0, g_z_0_0_0_z_xz_0_0, g_z_0_0_0_z_yy_0_0, g_z_0_0_0_z_yz_0_0, g_z_0_0_0_z_zz_0_0, g_zz_xx_0_0, g_zz_xy_0_0, g_zz_xz_0_0, g_zz_yy_0_0, g_zz_yz_0_0, g_zz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xx_0_0[i] = -g_0_xx_0_0[i] + 2.0 * g_zz_xx_0_0[i] * a_exp;

        g_z_0_0_0_z_xy_0_0[i] = -g_0_xy_0_0[i] + 2.0 * g_zz_xy_0_0[i] * a_exp;

        g_z_0_0_0_z_xz_0_0[i] = -g_0_xz_0_0[i] + 2.0 * g_zz_xz_0_0[i] * a_exp;

        g_z_0_0_0_z_yy_0_0[i] = -g_0_yy_0_0[i] + 2.0 * g_zz_yy_0_0[i] * a_exp;

        g_z_0_0_0_z_yz_0_0[i] = -g_0_yz_0_0[i] + 2.0 * g_zz_yz_0_0[i] * a_exp;

        g_z_0_0_0_z_zz_0_0[i] = -g_0_zz_0_0[i] + 2.0 * g_zz_zz_0_0[i] * a_exp;
    }
}

} // t4c_geom namespace


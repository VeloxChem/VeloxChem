#include "GeomDeriv1010OfScalarForSDSS.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1010_sdss_0(CSimdArray<double>& buffer_1010_sdss,
                     const CSimdArray<double>& buffer_pdps,
                     const double a_exp,
                     const double* c_exps) -> void
{
    const auto ndims = buffer_1010_sdss.number_of_columns();

    /// Set up components of auxilary buffer : buffer_pdps

    auto g_x_xx_x_0 = buffer_pdps[0];

    auto g_x_xx_y_0 = buffer_pdps[1];

    auto g_x_xx_z_0 = buffer_pdps[2];

    auto g_x_xy_x_0 = buffer_pdps[3];

    auto g_x_xy_y_0 = buffer_pdps[4];

    auto g_x_xy_z_0 = buffer_pdps[5];

    auto g_x_xz_x_0 = buffer_pdps[6];

    auto g_x_xz_y_0 = buffer_pdps[7];

    auto g_x_xz_z_0 = buffer_pdps[8];

    auto g_x_yy_x_0 = buffer_pdps[9];

    auto g_x_yy_y_0 = buffer_pdps[10];

    auto g_x_yy_z_0 = buffer_pdps[11];

    auto g_x_yz_x_0 = buffer_pdps[12];

    auto g_x_yz_y_0 = buffer_pdps[13];

    auto g_x_yz_z_0 = buffer_pdps[14];

    auto g_x_zz_x_0 = buffer_pdps[15];

    auto g_x_zz_y_0 = buffer_pdps[16];

    auto g_x_zz_z_0 = buffer_pdps[17];

    auto g_y_xx_x_0 = buffer_pdps[18];

    auto g_y_xx_y_0 = buffer_pdps[19];

    auto g_y_xx_z_0 = buffer_pdps[20];

    auto g_y_xy_x_0 = buffer_pdps[21];

    auto g_y_xy_y_0 = buffer_pdps[22];

    auto g_y_xy_z_0 = buffer_pdps[23];

    auto g_y_xz_x_0 = buffer_pdps[24];

    auto g_y_xz_y_0 = buffer_pdps[25];

    auto g_y_xz_z_0 = buffer_pdps[26];

    auto g_y_yy_x_0 = buffer_pdps[27];

    auto g_y_yy_y_0 = buffer_pdps[28];

    auto g_y_yy_z_0 = buffer_pdps[29];

    auto g_y_yz_x_0 = buffer_pdps[30];

    auto g_y_yz_y_0 = buffer_pdps[31];

    auto g_y_yz_z_0 = buffer_pdps[32];

    auto g_y_zz_x_0 = buffer_pdps[33];

    auto g_y_zz_y_0 = buffer_pdps[34];

    auto g_y_zz_z_0 = buffer_pdps[35];

    auto g_z_xx_x_0 = buffer_pdps[36];

    auto g_z_xx_y_0 = buffer_pdps[37];

    auto g_z_xx_z_0 = buffer_pdps[38];

    auto g_z_xy_x_0 = buffer_pdps[39];

    auto g_z_xy_y_0 = buffer_pdps[40];

    auto g_z_xy_z_0 = buffer_pdps[41];

    auto g_z_xz_x_0 = buffer_pdps[42];

    auto g_z_xz_y_0 = buffer_pdps[43];

    auto g_z_xz_z_0 = buffer_pdps[44];

    auto g_z_yy_x_0 = buffer_pdps[45];

    auto g_z_yy_y_0 = buffer_pdps[46];

    auto g_z_yy_z_0 = buffer_pdps[47];

    auto g_z_yz_x_0 = buffer_pdps[48];

    auto g_z_yz_y_0 = buffer_pdps[49];

    auto g_z_yz_z_0 = buffer_pdps[50];

    auto g_z_zz_x_0 = buffer_pdps[51];

    auto g_z_zz_y_0 = buffer_pdps[52];

    auto g_z_zz_z_0 = buffer_pdps[53];

    /// Set up components of integrals buffer : buffer_1010_sdss

    auto g_x_0_x_0_0_xx_0_0 = buffer_1010_sdss[0];

    auto g_x_0_x_0_0_xy_0_0 = buffer_1010_sdss[1];

    auto g_x_0_x_0_0_xz_0_0 = buffer_1010_sdss[2];

    auto g_x_0_x_0_0_yy_0_0 = buffer_1010_sdss[3];

    auto g_x_0_x_0_0_yz_0_0 = buffer_1010_sdss[4];

    auto g_x_0_x_0_0_zz_0_0 = buffer_1010_sdss[5];

    auto g_x_0_y_0_0_xx_0_0 = buffer_1010_sdss[6];

    auto g_x_0_y_0_0_xy_0_0 = buffer_1010_sdss[7];

    auto g_x_0_y_0_0_xz_0_0 = buffer_1010_sdss[8];

    auto g_x_0_y_0_0_yy_0_0 = buffer_1010_sdss[9];

    auto g_x_0_y_0_0_yz_0_0 = buffer_1010_sdss[10];

    auto g_x_0_y_0_0_zz_0_0 = buffer_1010_sdss[11];

    auto g_x_0_z_0_0_xx_0_0 = buffer_1010_sdss[12];

    auto g_x_0_z_0_0_xy_0_0 = buffer_1010_sdss[13];

    auto g_x_0_z_0_0_xz_0_0 = buffer_1010_sdss[14];

    auto g_x_0_z_0_0_yy_0_0 = buffer_1010_sdss[15];

    auto g_x_0_z_0_0_yz_0_0 = buffer_1010_sdss[16];

    auto g_x_0_z_0_0_zz_0_0 = buffer_1010_sdss[17];

    auto g_y_0_x_0_0_xx_0_0 = buffer_1010_sdss[18];

    auto g_y_0_x_0_0_xy_0_0 = buffer_1010_sdss[19];

    auto g_y_0_x_0_0_xz_0_0 = buffer_1010_sdss[20];

    auto g_y_0_x_0_0_yy_0_0 = buffer_1010_sdss[21];

    auto g_y_0_x_0_0_yz_0_0 = buffer_1010_sdss[22];

    auto g_y_0_x_0_0_zz_0_0 = buffer_1010_sdss[23];

    auto g_y_0_y_0_0_xx_0_0 = buffer_1010_sdss[24];

    auto g_y_0_y_0_0_xy_0_0 = buffer_1010_sdss[25];

    auto g_y_0_y_0_0_xz_0_0 = buffer_1010_sdss[26];

    auto g_y_0_y_0_0_yy_0_0 = buffer_1010_sdss[27];

    auto g_y_0_y_0_0_yz_0_0 = buffer_1010_sdss[28];

    auto g_y_0_y_0_0_zz_0_0 = buffer_1010_sdss[29];

    auto g_y_0_z_0_0_xx_0_0 = buffer_1010_sdss[30];

    auto g_y_0_z_0_0_xy_0_0 = buffer_1010_sdss[31];

    auto g_y_0_z_0_0_xz_0_0 = buffer_1010_sdss[32];

    auto g_y_0_z_0_0_yy_0_0 = buffer_1010_sdss[33];

    auto g_y_0_z_0_0_yz_0_0 = buffer_1010_sdss[34];

    auto g_y_0_z_0_0_zz_0_0 = buffer_1010_sdss[35];

    auto g_z_0_x_0_0_xx_0_0 = buffer_1010_sdss[36];

    auto g_z_0_x_0_0_xy_0_0 = buffer_1010_sdss[37];

    auto g_z_0_x_0_0_xz_0_0 = buffer_1010_sdss[38];

    auto g_z_0_x_0_0_yy_0_0 = buffer_1010_sdss[39];

    auto g_z_0_x_0_0_yz_0_0 = buffer_1010_sdss[40];

    auto g_z_0_x_0_0_zz_0_0 = buffer_1010_sdss[41];

    auto g_z_0_y_0_0_xx_0_0 = buffer_1010_sdss[42];

    auto g_z_0_y_0_0_xy_0_0 = buffer_1010_sdss[43];

    auto g_z_0_y_0_0_xz_0_0 = buffer_1010_sdss[44];

    auto g_z_0_y_0_0_yy_0_0 = buffer_1010_sdss[45];

    auto g_z_0_y_0_0_yz_0_0 = buffer_1010_sdss[46];

    auto g_z_0_y_0_0_zz_0_0 = buffer_1010_sdss[47];

    auto g_z_0_z_0_0_xx_0_0 = buffer_1010_sdss[48];

    auto g_z_0_z_0_0_xy_0_0 = buffer_1010_sdss[49];

    auto g_z_0_z_0_0_xz_0_0 = buffer_1010_sdss[50];

    auto g_z_0_z_0_0_yy_0_0 = buffer_1010_sdss[51];

    auto g_z_0_z_0_0_yz_0_0 = buffer_1010_sdss[52];

    auto g_z_0_z_0_0_zz_0_0 = buffer_1010_sdss[53];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_0_x_0_0_xx_0_0, g_x_0_x_0_0_xy_0_0, g_x_0_x_0_0_xz_0_0, g_x_0_x_0_0_yy_0_0, g_x_0_x_0_0_yz_0_0, g_x_0_x_0_0_zz_0_0, g_x_xx_x_0, g_x_xy_x_0, g_x_xz_x_0, g_x_yy_x_0, g_x_yz_x_0, g_x_zz_x_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xx_0_0[i] = 4.0 * g_x_xx_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_0_0[i] = 4.0 * g_x_xy_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_0_0[i] = 4.0 * g_x_xz_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_0_0[i] = 4.0 * g_x_yy_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_0_0[i] = 4.0 * g_x_yz_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_0_0[i] = 4.0 * g_x_zz_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_0_y_0_0_xx_0_0, g_x_0_y_0_0_xy_0_0, g_x_0_y_0_0_xz_0_0, g_x_0_y_0_0_yy_0_0, g_x_0_y_0_0_yz_0_0, g_x_0_y_0_0_zz_0_0, g_x_xx_y_0, g_x_xy_y_0, g_x_xz_y_0, g_x_yy_y_0, g_x_yz_y_0, g_x_zz_y_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xx_0_0[i] = 4.0 * g_x_xx_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_0_0[i] = 4.0 * g_x_xy_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_0_0[i] = 4.0 * g_x_xz_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_0_0[i] = 4.0 * g_x_yy_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_0_0[i] = 4.0 * g_x_yz_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_0_0[i] = 4.0 * g_x_zz_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_0_z_0_0_xx_0_0, g_x_0_z_0_0_xy_0_0, g_x_0_z_0_0_xz_0_0, g_x_0_z_0_0_yy_0_0, g_x_0_z_0_0_yz_0_0, g_x_0_z_0_0_zz_0_0, g_x_xx_z_0, g_x_xy_z_0, g_x_xz_z_0, g_x_yy_z_0, g_x_yz_z_0, g_x_zz_z_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xx_0_0[i] = 4.0 * g_x_xx_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_0_0[i] = 4.0 * g_x_xy_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_0_0[i] = 4.0 * g_x_xz_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_0_0[i] = 4.0 * g_x_yy_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_0_0[i] = 4.0 * g_x_yz_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_0_0[i] = 4.0 * g_x_zz_z_0[i] * a_exp * c_exps[i];
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_y_0_x_0_0_xx_0_0, g_y_0_x_0_0_xy_0_0, g_y_0_x_0_0_xz_0_0, g_y_0_x_0_0_yy_0_0, g_y_0_x_0_0_yz_0_0, g_y_0_x_0_0_zz_0_0, g_y_xx_x_0, g_y_xy_x_0, g_y_xz_x_0, g_y_yy_x_0, g_y_yz_x_0, g_y_zz_x_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xx_0_0[i] = 4.0 * g_y_xx_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_0_0[i] = 4.0 * g_y_xy_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_0_0[i] = 4.0 * g_y_xz_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_0_0[i] = 4.0 * g_y_yy_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_0_0[i] = 4.0 * g_y_yz_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_0_0[i] = 4.0 * g_y_zz_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_y_0_y_0_0_xx_0_0, g_y_0_y_0_0_xy_0_0, g_y_0_y_0_0_xz_0_0, g_y_0_y_0_0_yy_0_0, g_y_0_y_0_0_yz_0_0, g_y_0_y_0_0_zz_0_0, g_y_xx_y_0, g_y_xy_y_0, g_y_xz_y_0, g_y_yy_y_0, g_y_yz_y_0, g_y_zz_y_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xx_0_0[i] = 4.0 * g_y_xx_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_0_0[i] = 4.0 * g_y_xy_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_0_0[i] = 4.0 * g_y_xz_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_0_0[i] = 4.0 * g_y_yy_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_0_0[i] = 4.0 * g_y_yz_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_0_0[i] = 4.0 * g_y_zz_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_y_0_z_0_0_xx_0_0, g_y_0_z_0_0_xy_0_0, g_y_0_z_0_0_xz_0_0, g_y_0_z_0_0_yy_0_0, g_y_0_z_0_0_yz_0_0, g_y_0_z_0_0_zz_0_0, g_y_xx_z_0, g_y_xy_z_0, g_y_xz_z_0, g_y_yy_z_0, g_y_yz_z_0, g_y_zz_z_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xx_0_0[i] = 4.0 * g_y_xx_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_0_0[i] = 4.0 * g_y_xy_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_0_0[i] = 4.0 * g_y_xz_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_0_0[i] = 4.0 * g_y_yy_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_0_0[i] = 4.0 * g_y_yz_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_0_0[i] = 4.0 * g_y_zz_z_0[i] * a_exp * c_exps[i];
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_z_0_x_0_0_xx_0_0, g_z_0_x_0_0_xy_0_0, g_z_0_x_0_0_xz_0_0, g_z_0_x_0_0_yy_0_0, g_z_0_x_0_0_yz_0_0, g_z_0_x_0_0_zz_0_0, g_z_xx_x_0, g_z_xy_x_0, g_z_xz_x_0, g_z_yy_x_0, g_z_yz_x_0, g_z_zz_x_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xx_0_0[i] = 4.0 * g_z_xx_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_0_0[i] = 4.0 * g_z_xy_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_0_0[i] = 4.0 * g_z_xz_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_0_0[i] = 4.0 * g_z_yy_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_0_0[i] = 4.0 * g_z_yz_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_0_0[i] = 4.0 * g_z_zz_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_z_0_y_0_0_xx_0_0, g_z_0_y_0_0_xy_0_0, g_z_0_y_0_0_xz_0_0, g_z_0_y_0_0_yy_0_0, g_z_0_y_0_0_yz_0_0, g_z_0_y_0_0_zz_0_0, g_z_xx_y_0, g_z_xy_y_0, g_z_xz_y_0, g_z_yy_y_0, g_z_yz_y_0, g_z_zz_y_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xx_0_0[i] = 4.0 * g_z_xx_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_0_0[i] = 4.0 * g_z_xy_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_0_0[i] = 4.0 * g_z_xz_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_0_0[i] = 4.0 * g_z_yy_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_0_0[i] = 4.0 * g_z_yz_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_0_0[i] = 4.0 * g_z_zz_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_z_0_z_0_0_xx_0_0, g_z_0_z_0_0_xy_0_0, g_z_0_z_0_0_xz_0_0, g_z_0_z_0_0_yy_0_0, g_z_0_z_0_0_yz_0_0, g_z_0_z_0_0_zz_0_0, g_z_xx_z_0, g_z_xy_z_0, g_z_xz_z_0, g_z_yy_z_0, g_z_yz_z_0, g_z_zz_z_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xx_0_0[i] = 4.0 * g_z_xx_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_0_0[i] = 4.0 * g_z_xy_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_0_0[i] = 4.0 * g_z_xz_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_0_0[i] = 4.0 * g_z_yy_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_0_0[i] = 4.0 * g_z_yz_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_0_0[i] = 4.0 * g_z_zz_z_0[i] * a_exp * c_exps[i];
    }
}

} // t4c_geom namespace


#include "GeomDeriv1010OfScalarForPPSS.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1010_ppss_0(CSimdArray<double>& buffer_1010_ppss,
                     const CSimdArray<double>& buffer_spps,
                     const CSimdArray<double>& buffer_dpps,
                     const double a_exp,
                     const double* c_exps) -> void
{
    const auto ndims = buffer_1010_ppss.number_of_columns();

    /// Set up components of auxilary buffer : buffer_spps

    auto g_0_x_x_0 = buffer_spps[0];

    auto g_0_x_y_0 = buffer_spps[1];

    auto g_0_x_z_0 = buffer_spps[2];

    auto g_0_y_x_0 = buffer_spps[3];

    auto g_0_y_y_0 = buffer_spps[4];

    auto g_0_y_z_0 = buffer_spps[5];

    auto g_0_z_x_0 = buffer_spps[6];

    auto g_0_z_y_0 = buffer_spps[7];

    auto g_0_z_z_0 = buffer_spps[8];

    /// Set up components of auxilary buffer : buffer_dpps

    auto g_xx_x_x_0 = buffer_dpps[0];

    auto g_xx_x_y_0 = buffer_dpps[1];

    auto g_xx_x_z_0 = buffer_dpps[2];

    auto g_xx_y_x_0 = buffer_dpps[3];

    auto g_xx_y_y_0 = buffer_dpps[4];

    auto g_xx_y_z_0 = buffer_dpps[5];

    auto g_xx_z_x_0 = buffer_dpps[6];

    auto g_xx_z_y_0 = buffer_dpps[7];

    auto g_xx_z_z_0 = buffer_dpps[8];

    auto g_xy_x_x_0 = buffer_dpps[9];

    auto g_xy_x_y_0 = buffer_dpps[10];

    auto g_xy_x_z_0 = buffer_dpps[11];

    auto g_xy_y_x_0 = buffer_dpps[12];

    auto g_xy_y_y_0 = buffer_dpps[13];

    auto g_xy_y_z_0 = buffer_dpps[14];

    auto g_xy_z_x_0 = buffer_dpps[15];

    auto g_xy_z_y_0 = buffer_dpps[16];

    auto g_xy_z_z_0 = buffer_dpps[17];

    auto g_xz_x_x_0 = buffer_dpps[18];

    auto g_xz_x_y_0 = buffer_dpps[19];

    auto g_xz_x_z_0 = buffer_dpps[20];

    auto g_xz_y_x_0 = buffer_dpps[21];

    auto g_xz_y_y_0 = buffer_dpps[22];

    auto g_xz_y_z_0 = buffer_dpps[23];

    auto g_xz_z_x_0 = buffer_dpps[24];

    auto g_xz_z_y_0 = buffer_dpps[25];

    auto g_xz_z_z_0 = buffer_dpps[26];

    auto g_yy_x_x_0 = buffer_dpps[27];

    auto g_yy_x_y_0 = buffer_dpps[28];

    auto g_yy_x_z_0 = buffer_dpps[29];

    auto g_yy_y_x_0 = buffer_dpps[30];

    auto g_yy_y_y_0 = buffer_dpps[31];

    auto g_yy_y_z_0 = buffer_dpps[32];

    auto g_yy_z_x_0 = buffer_dpps[33];

    auto g_yy_z_y_0 = buffer_dpps[34];

    auto g_yy_z_z_0 = buffer_dpps[35];

    auto g_yz_x_x_0 = buffer_dpps[36];

    auto g_yz_x_y_0 = buffer_dpps[37];

    auto g_yz_x_z_0 = buffer_dpps[38];

    auto g_yz_y_x_0 = buffer_dpps[39];

    auto g_yz_y_y_0 = buffer_dpps[40];

    auto g_yz_y_z_0 = buffer_dpps[41];

    auto g_yz_z_x_0 = buffer_dpps[42];

    auto g_yz_z_y_0 = buffer_dpps[43];

    auto g_yz_z_z_0 = buffer_dpps[44];

    auto g_zz_x_x_0 = buffer_dpps[45];

    auto g_zz_x_y_0 = buffer_dpps[46];

    auto g_zz_x_z_0 = buffer_dpps[47];

    auto g_zz_y_x_0 = buffer_dpps[48];

    auto g_zz_y_y_0 = buffer_dpps[49];

    auto g_zz_y_z_0 = buffer_dpps[50];

    auto g_zz_z_x_0 = buffer_dpps[51];

    auto g_zz_z_y_0 = buffer_dpps[52];

    auto g_zz_z_z_0 = buffer_dpps[53];

    /// Set up components of integrals buffer : buffer_1010_ppss

    auto g_x_0_x_0_x_x_0_0 = buffer_1010_ppss[0];

    auto g_x_0_x_0_x_y_0_0 = buffer_1010_ppss[1];

    auto g_x_0_x_0_x_z_0_0 = buffer_1010_ppss[2];

    auto g_x_0_x_0_y_x_0_0 = buffer_1010_ppss[3];

    auto g_x_0_x_0_y_y_0_0 = buffer_1010_ppss[4];

    auto g_x_0_x_0_y_z_0_0 = buffer_1010_ppss[5];

    auto g_x_0_x_0_z_x_0_0 = buffer_1010_ppss[6];

    auto g_x_0_x_0_z_y_0_0 = buffer_1010_ppss[7];

    auto g_x_0_x_0_z_z_0_0 = buffer_1010_ppss[8];

    auto g_x_0_y_0_x_x_0_0 = buffer_1010_ppss[9];

    auto g_x_0_y_0_x_y_0_0 = buffer_1010_ppss[10];

    auto g_x_0_y_0_x_z_0_0 = buffer_1010_ppss[11];

    auto g_x_0_y_0_y_x_0_0 = buffer_1010_ppss[12];

    auto g_x_0_y_0_y_y_0_0 = buffer_1010_ppss[13];

    auto g_x_0_y_0_y_z_0_0 = buffer_1010_ppss[14];

    auto g_x_0_y_0_z_x_0_0 = buffer_1010_ppss[15];

    auto g_x_0_y_0_z_y_0_0 = buffer_1010_ppss[16];

    auto g_x_0_y_0_z_z_0_0 = buffer_1010_ppss[17];

    auto g_x_0_z_0_x_x_0_0 = buffer_1010_ppss[18];

    auto g_x_0_z_0_x_y_0_0 = buffer_1010_ppss[19];

    auto g_x_0_z_0_x_z_0_0 = buffer_1010_ppss[20];

    auto g_x_0_z_0_y_x_0_0 = buffer_1010_ppss[21];

    auto g_x_0_z_0_y_y_0_0 = buffer_1010_ppss[22];

    auto g_x_0_z_0_y_z_0_0 = buffer_1010_ppss[23];

    auto g_x_0_z_0_z_x_0_0 = buffer_1010_ppss[24];

    auto g_x_0_z_0_z_y_0_0 = buffer_1010_ppss[25];

    auto g_x_0_z_0_z_z_0_0 = buffer_1010_ppss[26];

    auto g_y_0_x_0_x_x_0_0 = buffer_1010_ppss[27];

    auto g_y_0_x_0_x_y_0_0 = buffer_1010_ppss[28];

    auto g_y_0_x_0_x_z_0_0 = buffer_1010_ppss[29];

    auto g_y_0_x_0_y_x_0_0 = buffer_1010_ppss[30];

    auto g_y_0_x_0_y_y_0_0 = buffer_1010_ppss[31];

    auto g_y_0_x_0_y_z_0_0 = buffer_1010_ppss[32];

    auto g_y_0_x_0_z_x_0_0 = buffer_1010_ppss[33];

    auto g_y_0_x_0_z_y_0_0 = buffer_1010_ppss[34];

    auto g_y_0_x_0_z_z_0_0 = buffer_1010_ppss[35];

    auto g_y_0_y_0_x_x_0_0 = buffer_1010_ppss[36];

    auto g_y_0_y_0_x_y_0_0 = buffer_1010_ppss[37];

    auto g_y_0_y_0_x_z_0_0 = buffer_1010_ppss[38];

    auto g_y_0_y_0_y_x_0_0 = buffer_1010_ppss[39];

    auto g_y_0_y_0_y_y_0_0 = buffer_1010_ppss[40];

    auto g_y_0_y_0_y_z_0_0 = buffer_1010_ppss[41];

    auto g_y_0_y_0_z_x_0_0 = buffer_1010_ppss[42];

    auto g_y_0_y_0_z_y_0_0 = buffer_1010_ppss[43];

    auto g_y_0_y_0_z_z_0_0 = buffer_1010_ppss[44];

    auto g_y_0_z_0_x_x_0_0 = buffer_1010_ppss[45];

    auto g_y_0_z_0_x_y_0_0 = buffer_1010_ppss[46];

    auto g_y_0_z_0_x_z_0_0 = buffer_1010_ppss[47];

    auto g_y_0_z_0_y_x_0_0 = buffer_1010_ppss[48];

    auto g_y_0_z_0_y_y_0_0 = buffer_1010_ppss[49];

    auto g_y_0_z_0_y_z_0_0 = buffer_1010_ppss[50];

    auto g_y_0_z_0_z_x_0_0 = buffer_1010_ppss[51];

    auto g_y_0_z_0_z_y_0_0 = buffer_1010_ppss[52];

    auto g_y_0_z_0_z_z_0_0 = buffer_1010_ppss[53];

    auto g_z_0_x_0_x_x_0_0 = buffer_1010_ppss[54];

    auto g_z_0_x_0_x_y_0_0 = buffer_1010_ppss[55];

    auto g_z_0_x_0_x_z_0_0 = buffer_1010_ppss[56];

    auto g_z_0_x_0_y_x_0_0 = buffer_1010_ppss[57];

    auto g_z_0_x_0_y_y_0_0 = buffer_1010_ppss[58];

    auto g_z_0_x_0_y_z_0_0 = buffer_1010_ppss[59];

    auto g_z_0_x_0_z_x_0_0 = buffer_1010_ppss[60];

    auto g_z_0_x_0_z_y_0_0 = buffer_1010_ppss[61];

    auto g_z_0_x_0_z_z_0_0 = buffer_1010_ppss[62];

    auto g_z_0_y_0_x_x_0_0 = buffer_1010_ppss[63];

    auto g_z_0_y_0_x_y_0_0 = buffer_1010_ppss[64];

    auto g_z_0_y_0_x_z_0_0 = buffer_1010_ppss[65];

    auto g_z_0_y_0_y_x_0_0 = buffer_1010_ppss[66];

    auto g_z_0_y_0_y_y_0_0 = buffer_1010_ppss[67];

    auto g_z_0_y_0_y_z_0_0 = buffer_1010_ppss[68];

    auto g_z_0_y_0_z_x_0_0 = buffer_1010_ppss[69];

    auto g_z_0_y_0_z_y_0_0 = buffer_1010_ppss[70];

    auto g_z_0_y_0_z_z_0_0 = buffer_1010_ppss[71];

    auto g_z_0_z_0_x_x_0_0 = buffer_1010_ppss[72];

    auto g_z_0_z_0_x_y_0_0 = buffer_1010_ppss[73];

    auto g_z_0_z_0_x_z_0_0 = buffer_1010_ppss[74];

    auto g_z_0_z_0_y_x_0_0 = buffer_1010_ppss[75];

    auto g_z_0_z_0_y_y_0_0 = buffer_1010_ppss[76];

    auto g_z_0_z_0_y_z_0_0 = buffer_1010_ppss[77];

    auto g_z_0_z_0_z_x_0_0 = buffer_1010_ppss[78];

    auto g_z_0_z_0_z_y_0_0 = buffer_1010_ppss[79];

    auto g_z_0_z_0_z_z_0_0 = buffer_1010_ppss[80];

    // integrals block (0-3)

    #pragma omp simd aligned(g_0_x_x_0, g_0_y_x_0, g_0_z_x_0, g_x_0_x_0_x_x_0_0, g_x_0_x_0_x_y_0_0, g_x_0_x_0_x_z_0_0, g_xx_x_x_0, g_xx_y_x_0, g_xx_z_x_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_x_0_0[i] = -2.0 * g_0_x_x_0[i] * c_exps[i] + 4.0 * g_xx_x_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_0_0[i] = -2.0 * g_0_y_x_0[i] * c_exps[i] + 4.0 * g_xx_y_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_0_0[i] = -2.0 * g_0_z_x_0[i] * c_exps[i] + 4.0 * g_xx_z_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_x_0_x_0_y_x_0_0, g_x_0_x_0_y_y_0_0, g_x_0_x_0_y_z_0_0, g_xy_x_x_0, g_xy_y_x_0, g_xy_z_x_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_x_0_0[i] = 4.0 * g_xy_x_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_0_0[i] = 4.0 * g_xy_y_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_0_0[i] = 4.0 * g_xy_z_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_x_0_x_0_z_x_0_0, g_x_0_x_0_z_y_0_0, g_x_0_x_0_z_z_0_0, g_xz_x_x_0, g_xz_y_x_0, g_xz_z_x_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_x_0_0[i] = 4.0 * g_xz_x_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_0_0[i] = 4.0 * g_xz_y_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_0_0[i] = 4.0 * g_xz_z_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_0_x_y_0, g_0_y_y_0, g_0_z_y_0, g_x_0_y_0_x_x_0_0, g_x_0_y_0_x_y_0_0, g_x_0_y_0_x_z_0_0, g_xx_x_y_0, g_xx_y_y_0, g_xx_z_y_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_x_0_0[i] = -2.0 * g_0_x_y_0[i] * c_exps[i] + 4.0 * g_xx_x_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_0_0[i] = -2.0 * g_0_y_y_0[i] * c_exps[i] + 4.0 * g_xx_y_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_0_0[i] = -2.0 * g_0_z_y_0[i] * c_exps[i] + 4.0 * g_xx_z_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_x_0_y_0_y_x_0_0, g_x_0_y_0_y_y_0_0, g_x_0_y_0_y_z_0_0, g_xy_x_y_0, g_xy_y_y_0, g_xy_z_y_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_x_0_0[i] = 4.0 * g_xy_x_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_0_0[i] = 4.0 * g_xy_y_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_0_0[i] = 4.0 * g_xy_z_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_x_0_y_0_z_x_0_0, g_x_0_y_0_z_y_0_0, g_x_0_y_0_z_z_0_0, g_xz_x_y_0, g_xz_y_y_0, g_xz_z_y_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_x_0_0[i] = 4.0 * g_xz_x_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_0_0[i] = 4.0 * g_xz_y_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_0_0[i] = 4.0 * g_xz_z_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_0_x_z_0, g_0_y_z_0, g_0_z_z_0, g_x_0_z_0_x_x_0_0, g_x_0_z_0_x_y_0_0, g_x_0_z_0_x_z_0_0, g_xx_x_z_0, g_xx_y_z_0, g_xx_z_z_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_x_0_0[i] = -2.0 * g_0_x_z_0[i] * c_exps[i] + 4.0 * g_xx_x_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_0_0[i] = -2.0 * g_0_y_z_0[i] * c_exps[i] + 4.0 * g_xx_y_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_0_0[i] = -2.0 * g_0_z_z_0[i] * c_exps[i] + 4.0 * g_xx_z_z_0[i] * a_exp * c_exps[i];
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_x_0_z_0_y_x_0_0, g_x_0_z_0_y_y_0_0, g_x_0_z_0_y_z_0_0, g_xy_x_z_0, g_xy_y_z_0, g_xy_z_z_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_x_0_0[i] = 4.0 * g_xy_x_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_0_0[i] = 4.0 * g_xy_y_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_0_0[i] = 4.0 * g_xy_z_z_0[i] * a_exp * c_exps[i];
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_x_0_z_0_z_x_0_0, g_x_0_z_0_z_y_0_0, g_x_0_z_0_z_z_0_0, g_xz_x_z_0, g_xz_y_z_0, g_xz_z_z_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_x_0_0[i] = 4.0 * g_xz_x_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_0_0[i] = 4.0 * g_xz_y_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_0_0[i] = 4.0 * g_xz_z_z_0[i] * a_exp * c_exps[i];
    }
    // integrals block (27-30)

    #pragma omp simd aligned(g_xy_x_x_0, g_xy_y_x_0, g_xy_z_x_0, g_y_0_x_0_x_x_0_0, g_y_0_x_0_x_y_0_0, g_y_0_x_0_x_z_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_x_0_0[i] = 4.0 * g_xy_x_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_0_0[i] = 4.0 * g_xy_y_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_0_0[i] = 4.0 * g_xy_z_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (30-33)

    #pragma omp simd aligned(g_0_x_x_0, g_0_y_x_0, g_0_z_x_0, g_y_0_x_0_y_x_0_0, g_y_0_x_0_y_y_0_0, g_y_0_x_0_y_z_0_0, g_yy_x_x_0, g_yy_y_x_0, g_yy_z_x_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_x_0_0[i] = -2.0 * g_0_x_x_0[i] * c_exps[i] + 4.0 * g_yy_x_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_0_0[i] = -2.0 * g_0_y_x_0[i] * c_exps[i] + 4.0 * g_yy_y_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_0_0[i] = -2.0 * g_0_z_x_0[i] * c_exps[i] + 4.0 * g_yy_z_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (33-36)

    #pragma omp simd aligned(g_y_0_x_0_z_x_0_0, g_y_0_x_0_z_y_0_0, g_y_0_x_0_z_z_0_0, g_yz_x_x_0, g_yz_y_x_0, g_yz_z_x_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_x_0_0[i] = 4.0 * g_yz_x_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_0_0[i] = 4.0 * g_yz_y_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_0_0[i] = 4.0 * g_yz_z_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (36-39)

    #pragma omp simd aligned(g_xy_x_y_0, g_xy_y_y_0, g_xy_z_y_0, g_y_0_y_0_x_x_0_0, g_y_0_y_0_x_y_0_0, g_y_0_y_0_x_z_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_x_0_0[i] = 4.0 * g_xy_x_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_0_0[i] = 4.0 * g_xy_y_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_0_0[i] = 4.0 * g_xy_z_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (39-42)

    #pragma omp simd aligned(g_0_x_y_0, g_0_y_y_0, g_0_z_y_0, g_y_0_y_0_y_x_0_0, g_y_0_y_0_y_y_0_0, g_y_0_y_0_y_z_0_0, g_yy_x_y_0, g_yy_y_y_0, g_yy_z_y_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_x_0_0[i] = -2.0 * g_0_x_y_0[i] * c_exps[i] + 4.0 * g_yy_x_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_0_0[i] = -2.0 * g_0_y_y_0[i] * c_exps[i] + 4.0 * g_yy_y_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_0_0[i] = -2.0 * g_0_z_y_0[i] * c_exps[i] + 4.0 * g_yy_z_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (42-45)

    #pragma omp simd aligned(g_y_0_y_0_z_x_0_0, g_y_0_y_0_z_y_0_0, g_y_0_y_0_z_z_0_0, g_yz_x_y_0, g_yz_y_y_0, g_yz_z_y_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_x_0_0[i] = 4.0 * g_yz_x_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_0_0[i] = 4.0 * g_yz_y_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_0_0[i] = 4.0 * g_yz_z_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (45-48)

    #pragma omp simd aligned(g_xy_x_z_0, g_xy_y_z_0, g_xy_z_z_0, g_y_0_z_0_x_x_0_0, g_y_0_z_0_x_y_0_0, g_y_0_z_0_x_z_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_x_0_0[i] = 4.0 * g_xy_x_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_0_0[i] = 4.0 * g_xy_y_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_0_0[i] = 4.0 * g_xy_z_z_0[i] * a_exp * c_exps[i];
    }
    // integrals block (48-51)

    #pragma omp simd aligned(g_0_x_z_0, g_0_y_z_0, g_0_z_z_0, g_y_0_z_0_y_x_0_0, g_y_0_z_0_y_y_0_0, g_y_0_z_0_y_z_0_0, g_yy_x_z_0, g_yy_y_z_0, g_yy_z_z_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_x_0_0[i] = -2.0 * g_0_x_z_0[i] * c_exps[i] + 4.0 * g_yy_x_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_0_0[i] = -2.0 * g_0_y_z_0[i] * c_exps[i] + 4.0 * g_yy_y_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_0_0[i] = -2.0 * g_0_z_z_0[i] * c_exps[i] + 4.0 * g_yy_z_z_0[i] * a_exp * c_exps[i];
    }
    // integrals block (51-54)

    #pragma omp simd aligned(g_y_0_z_0_z_x_0_0, g_y_0_z_0_z_y_0_0, g_y_0_z_0_z_z_0_0, g_yz_x_z_0, g_yz_y_z_0, g_yz_z_z_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_x_0_0[i] = 4.0 * g_yz_x_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_0_0[i] = 4.0 * g_yz_y_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_0_0[i] = 4.0 * g_yz_z_z_0[i] * a_exp * c_exps[i];
    }
    // integrals block (54-57)

    #pragma omp simd aligned(g_xz_x_x_0, g_xz_y_x_0, g_xz_z_x_0, g_z_0_x_0_x_x_0_0, g_z_0_x_0_x_y_0_0, g_z_0_x_0_x_z_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_x_0_0[i] = 4.0 * g_xz_x_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_0_0[i] = 4.0 * g_xz_y_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_0_0[i] = 4.0 * g_xz_z_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (57-60)

    #pragma omp simd aligned(g_yz_x_x_0, g_yz_y_x_0, g_yz_z_x_0, g_z_0_x_0_y_x_0_0, g_z_0_x_0_y_y_0_0, g_z_0_x_0_y_z_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_x_0_0[i] = 4.0 * g_yz_x_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_0_0[i] = 4.0 * g_yz_y_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_0_0[i] = 4.0 * g_yz_z_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (60-63)

    #pragma omp simd aligned(g_0_x_x_0, g_0_y_x_0, g_0_z_x_0, g_z_0_x_0_z_x_0_0, g_z_0_x_0_z_y_0_0, g_z_0_x_0_z_z_0_0, g_zz_x_x_0, g_zz_y_x_0, g_zz_z_x_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_x_0_0[i] = -2.0 * g_0_x_x_0[i] * c_exps[i] + 4.0 * g_zz_x_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_0_0[i] = -2.0 * g_0_y_x_0[i] * c_exps[i] + 4.0 * g_zz_y_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_0_0[i] = -2.0 * g_0_z_x_0[i] * c_exps[i] + 4.0 * g_zz_z_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (63-66)

    #pragma omp simd aligned(g_xz_x_y_0, g_xz_y_y_0, g_xz_z_y_0, g_z_0_y_0_x_x_0_0, g_z_0_y_0_x_y_0_0, g_z_0_y_0_x_z_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_x_0_0[i] = 4.0 * g_xz_x_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_0_0[i] = 4.0 * g_xz_y_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_0_0[i] = 4.0 * g_xz_z_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (66-69)

    #pragma omp simd aligned(g_yz_x_y_0, g_yz_y_y_0, g_yz_z_y_0, g_z_0_y_0_y_x_0_0, g_z_0_y_0_y_y_0_0, g_z_0_y_0_y_z_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_x_0_0[i] = 4.0 * g_yz_x_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_0_0[i] = 4.0 * g_yz_y_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_0_0[i] = 4.0 * g_yz_z_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (69-72)

    #pragma omp simd aligned(g_0_x_y_0, g_0_y_y_0, g_0_z_y_0, g_z_0_y_0_z_x_0_0, g_z_0_y_0_z_y_0_0, g_z_0_y_0_z_z_0_0, g_zz_x_y_0, g_zz_y_y_0, g_zz_z_y_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_x_0_0[i] = -2.0 * g_0_x_y_0[i] * c_exps[i] + 4.0 * g_zz_x_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_0_0[i] = -2.0 * g_0_y_y_0[i] * c_exps[i] + 4.0 * g_zz_y_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_0_0[i] = -2.0 * g_0_z_y_0[i] * c_exps[i] + 4.0 * g_zz_z_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (72-75)

    #pragma omp simd aligned(g_xz_x_z_0, g_xz_y_z_0, g_xz_z_z_0, g_z_0_z_0_x_x_0_0, g_z_0_z_0_x_y_0_0, g_z_0_z_0_x_z_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_x_0_0[i] = 4.0 * g_xz_x_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_0_0[i] = 4.0 * g_xz_y_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_0_0[i] = 4.0 * g_xz_z_z_0[i] * a_exp * c_exps[i];
    }
    // integrals block (75-78)

    #pragma omp simd aligned(g_yz_x_z_0, g_yz_y_z_0, g_yz_z_z_0, g_z_0_z_0_y_x_0_0, g_z_0_z_0_y_y_0_0, g_z_0_z_0_y_z_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_x_0_0[i] = 4.0 * g_yz_x_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_0_0[i] = 4.0 * g_yz_y_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_0_0[i] = 4.0 * g_yz_z_z_0[i] * a_exp * c_exps[i];
    }
    // integrals block (78-81)

    #pragma omp simd aligned(g_0_x_z_0, g_0_y_z_0, g_0_z_z_0, g_z_0_z_0_z_x_0_0, g_z_0_z_0_z_y_0_0, g_z_0_z_0_z_z_0_0, g_zz_x_z_0, g_zz_y_z_0, g_zz_z_z_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_x_0_0[i] = -2.0 * g_0_x_z_0[i] * c_exps[i] + 4.0 * g_zz_x_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_0_0[i] = -2.0 * g_0_y_z_0[i] * c_exps[i] + 4.0 * g_zz_y_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_0_0[i] = -2.0 * g_0_z_z_0[i] * c_exps[i] + 4.0 * g_zz_z_z_0[i] * a_exp * c_exps[i];
    }
}

} // t4c_geom namespace


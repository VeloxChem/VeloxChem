#include "GeomDeriv1010OfScalarForSSPP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1010_sspp_0(CSimdArray<double>& buffer_1010_sspp,
                     const CSimdArray<double>& buffer_pssp,
                     const CSimdArray<double>& buffer_psdp,
                     const double a_exp,
                     const double* c_exps) -> void
{
    const auto ndims = buffer_1010_sspp.number_of_columns();

    /// Set up components of auxilary buffer : buffer_pssp

    auto g_x_0_0_x = buffer_pssp[0];

    auto g_x_0_0_y = buffer_pssp[1];

    auto g_x_0_0_z = buffer_pssp[2];

    auto g_y_0_0_x = buffer_pssp[3];

    auto g_y_0_0_y = buffer_pssp[4];

    auto g_y_0_0_z = buffer_pssp[5];

    auto g_z_0_0_x = buffer_pssp[6];

    auto g_z_0_0_y = buffer_pssp[7];

    auto g_z_0_0_z = buffer_pssp[8];

    /// Set up components of auxilary buffer : buffer_psdp

    auto g_x_0_xx_x = buffer_psdp[0];

    auto g_x_0_xx_y = buffer_psdp[1];

    auto g_x_0_xx_z = buffer_psdp[2];

    auto g_x_0_xy_x = buffer_psdp[3];

    auto g_x_0_xy_y = buffer_psdp[4];

    auto g_x_0_xy_z = buffer_psdp[5];

    auto g_x_0_xz_x = buffer_psdp[6];

    auto g_x_0_xz_y = buffer_psdp[7];

    auto g_x_0_xz_z = buffer_psdp[8];

    auto g_x_0_yy_x = buffer_psdp[9];

    auto g_x_0_yy_y = buffer_psdp[10];

    auto g_x_0_yy_z = buffer_psdp[11];

    auto g_x_0_yz_x = buffer_psdp[12];

    auto g_x_0_yz_y = buffer_psdp[13];

    auto g_x_0_yz_z = buffer_psdp[14];

    auto g_x_0_zz_x = buffer_psdp[15];

    auto g_x_0_zz_y = buffer_psdp[16];

    auto g_x_0_zz_z = buffer_psdp[17];

    auto g_y_0_xx_x = buffer_psdp[18];

    auto g_y_0_xx_y = buffer_psdp[19];

    auto g_y_0_xx_z = buffer_psdp[20];

    auto g_y_0_xy_x = buffer_psdp[21];

    auto g_y_0_xy_y = buffer_psdp[22];

    auto g_y_0_xy_z = buffer_psdp[23];

    auto g_y_0_xz_x = buffer_psdp[24];

    auto g_y_0_xz_y = buffer_psdp[25];

    auto g_y_0_xz_z = buffer_psdp[26];

    auto g_y_0_yy_x = buffer_psdp[27];

    auto g_y_0_yy_y = buffer_psdp[28];

    auto g_y_0_yy_z = buffer_psdp[29];

    auto g_y_0_yz_x = buffer_psdp[30];

    auto g_y_0_yz_y = buffer_psdp[31];

    auto g_y_0_yz_z = buffer_psdp[32];

    auto g_y_0_zz_x = buffer_psdp[33];

    auto g_y_0_zz_y = buffer_psdp[34];

    auto g_y_0_zz_z = buffer_psdp[35];

    auto g_z_0_xx_x = buffer_psdp[36];

    auto g_z_0_xx_y = buffer_psdp[37];

    auto g_z_0_xx_z = buffer_psdp[38];

    auto g_z_0_xy_x = buffer_psdp[39];

    auto g_z_0_xy_y = buffer_psdp[40];

    auto g_z_0_xy_z = buffer_psdp[41];

    auto g_z_0_xz_x = buffer_psdp[42];

    auto g_z_0_xz_y = buffer_psdp[43];

    auto g_z_0_xz_z = buffer_psdp[44];

    auto g_z_0_yy_x = buffer_psdp[45];

    auto g_z_0_yy_y = buffer_psdp[46];

    auto g_z_0_yy_z = buffer_psdp[47];

    auto g_z_0_yz_x = buffer_psdp[48];

    auto g_z_0_yz_y = buffer_psdp[49];

    auto g_z_0_yz_z = buffer_psdp[50];

    auto g_z_0_zz_x = buffer_psdp[51];

    auto g_z_0_zz_y = buffer_psdp[52];

    auto g_z_0_zz_z = buffer_psdp[53];

    /// Set up components of integrals buffer : buffer_1010_sspp

    auto g_x_0_x_0_0_0_x_x = buffer_1010_sspp[0];

    auto g_x_0_x_0_0_0_x_y = buffer_1010_sspp[1];

    auto g_x_0_x_0_0_0_x_z = buffer_1010_sspp[2];

    auto g_x_0_x_0_0_0_y_x = buffer_1010_sspp[3];

    auto g_x_0_x_0_0_0_y_y = buffer_1010_sspp[4];

    auto g_x_0_x_0_0_0_y_z = buffer_1010_sspp[5];

    auto g_x_0_x_0_0_0_z_x = buffer_1010_sspp[6];

    auto g_x_0_x_0_0_0_z_y = buffer_1010_sspp[7];

    auto g_x_0_x_0_0_0_z_z = buffer_1010_sspp[8];

    auto g_x_0_y_0_0_0_x_x = buffer_1010_sspp[9];

    auto g_x_0_y_0_0_0_x_y = buffer_1010_sspp[10];

    auto g_x_0_y_0_0_0_x_z = buffer_1010_sspp[11];

    auto g_x_0_y_0_0_0_y_x = buffer_1010_sspp[12];

    auto g_x_0_y_0_0_0_y_y = buffer_1010_sspp[13];

    auto g_x_0_y_0_0_0_y_z = buffer_1010_sspp[14];

    auto g_x_0_y_0_0_0_z_x = buffer_1010_sspp[15];

    auto g_x_0_y_0_0_0_z_y = buffer_1010_sspp[16];

    auto g_x_0_y_0_0_0_z_z = buffer_1010_sspp[17];

    auto g_x_0_z_0_0_0_x_x = buffer_1010_sspp[18];

    auto g_x_0_z_0_0_0_x_y = buffer_1010_sspp[19];

    auto g_x_0_z_0_0_0_x_z = buffer_1010_sspp[20];

    auto g_x_0_z_0_0_0_y_x = buffer_1010_sspp[21];

    auto g_x_0_z_0_0_0_y_y = buffer_1010_sspp[22];

    auto g_x_0_z_0_0_0_y_z = buffer_1010_sspp[23];

    auto g_x_0_z_0_0_0_z_x = buffer_1010_sspp[24];

    auto g_x_0_z_0_0_0_z_y = buffer_1010_sspp[25];

    auto g_x_0_z_0_0_0_z_z = buffer_1010_sspp[26];

    auto g_y_0_x_0_0_0_x_x = buffer_1010_sspp[27];

    auto g_y_0_x_0_0_0_x_y = buffer_1010_sspp[28];

    auto g_y_0_x_0_0_0_x_z = buffer_1010_sspp[29];

    auto g_y_0_x_0_0_0_y_x = buffer_1010_sspp[30];

    auto g_y_0_x_0_0_0_y_y = buffer_1010_sspp[31];

    auto g_y_0_x_0_0_0_y_z = buffer_1010_sspp[32];

    auto g_y_0_x_0_0_0_z_x = buffer_1010_sspp[33];

    auto g_y_0_x_0_0_0_z_y = buffer_1010_sspp[34];

    auto g_y_0_x_0_0_0_z_z = buffer_1010_sspp[35];

    auto g_y_0_y_0_0_0_x_x = buffer_1010_sspp[36];

    auto g_y_0_y_0_0_0_x_y = buffer_1010_sspp[37];

    auto g_y_0_y_0_0_0_x_z = buffer_1010_sspp[38];

    auto g_y_0_y_0_0_0_y_x = buffer_1010_sspp[39];

    auto g_y_0_y_0_0_0_y_y = buffer_1010_sspp[40];

    auto g_y_0_y_0_0_0_y_z = buffer_1010_sspp[41];

    auto g_y_0_y_0_0_0_z_x = buffer_1010_sspp[42];

    auto g_y_0_y_0_0_0_z_y = buffer_1010_sspp[43];

    auto g_y_0_y_0_0_0_z_z = buffer_1010_sspp[44];

    auto g_y_0_z_0_0_0_x_x = buffer_1010_sspp[45];

    auto g_y_0_z_0_0_0_x_y = buffer_1010_sspp[46];

    auto g_y_0_z_0_0_0_x_z = buffer_1010_sspp[47];

    auto g_y_0_z_0_0_0_y_x = buffer_1010_sspp[48];

    auto g_y_0_z_0_0_0_y_y = buffer_1010_sspp[49];

    auto g_y_0_z_0_0_0_y_z = buffer_1010_sspp[50];

    auto g_y_0_z_0_0_0_z_x = buffer_1010_sspp[51];

    auto g_y_0_z_0_0_0_z_y = buffer_1010_sspp[52];

    auto g_y_0_z_0_0_0_z_z = buffer_1010_sspp[53];

    auto g_z_0_x_0_0_0_x_x = buffer_1010_sspp[54];

    auto g_z_0_x_0_0_0_x_y = buffer_1010_sspp[55];

    auto g_z_0_x_0_0_0_x_z = buffer_1010_sspp[56];

    auto g_z_0_x_0_0_0_y_x = buffer_1010_sspp[57];

    auto g_z_0_x_0_0_0_y_y = buffer_1010_sspp[58];

    auto g_z_0_x_0_0_0_y_z = buffer_1010_sspp[59];

    auto g_z_0_x_0_0_0_z_x = buffer_1010_sspp[60];

    auto g_z_0_x_0_0_0_z_y = buffer_1010_sspp[61];

    auto g_z_0_x_0_0_0_z_z = buffer_1010_sspp[62];

    auto g_z_0_y_0_0_0_x_x = buffer_1010_sspp[63];

    auto g_z_0_y_0_0_0_x_y = buffer_1010_sspp[64];

    auto g_z_0_y_0_0_0_x_z = buffer_1010_sspp[65];

    auto g_z_0_y_0_0_0_y_x = buffer_1010_sspp[66];

    auto g_z_0_y_0_0_0_y_y = buffer_1010_sspp[67];

    auto g_z_0_y_0_0_0_y_z = buffer_1010_sspp[68];

    auto g_z_0_y_0_0_0_z_x = buffer_1010_sspp[69];

    auto g_z_0_y_0_0_0_z_y = buffer_1010_sspp[70];

    auto g_z_0_y_0_0_0_z_z = buffer_1010_sspp[71];

    auto g_z_0_z_0_0_0_x_x = buffer_1010_sspp[72];

    auto g_z_0_z_0_0_0_x_y = buffer_1010_sspp[73];

    auto g_z_0_z_0_0_0_x_z = buffer_1010_sspp[74];

    auto g_z_0_z_0_0_0_y_x = buffer_1010_sspp[75];

    auto g_z_0_z_0_0_0_y_y = buffer_1010_sspp[76];

    auto g_z_0_z_0_0_0_y_z = buffer_1010_sspp[77];

    auto g_z_0_z_0_0_0_z_x = buffer_1010_sspp[78];

    auto g_z_0_z_0_0_0_z_y = buffer_1010_sspp[79];

    auto g_z_0_z_0_0_0_z_z = buffer_1010_sspp[80];

    // integrals block (0-3)

    #pragma omp simd aligned(g_x_0_0_x, g_x_0_0_y, g_x_0_0_z, g_x_0_x_0_0_0_x_x, g_x_0_x_0_0_0_x_y, g_x_0_x_0_0_0_x_z, g_x_0_xx_x, g_x_0_xx_y, g_x_0_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_0_x_x[i] = -2.0 * g_x_0_0_x[i] * a_exp + 4.0 * g_x_0_xx_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_x_y[i] = -2.0 * g_x_0_0_y[i] * a_exp + 4.0 * g_x_0_xx_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_x_z[i] = -2.0 * g_x_0_0_z[i] * a_exp + 4.0 * g_x_0_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_x_0_x_0_0_0_y_x, g_x_0_x_0_0_0_y_y, g_x_0_x_0_0_0_y_z, g_x_0_xy_x, g_x_0_xy_y, g_x_0_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_0_y_x[i] = 4.0 * g_x_0_xy_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_y_y[i] = 4.0 * g_x_0_xy_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_y_z[i] = 4.0 * g_x_0_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_x_0_x_0_0_0_z_x, g_x_0_x_0_0_0_z_y, g_x_0_x_0_0_0_z_z, g_x_0_xz_x, g_x_0_xz_y, g_x_0_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_0_z_x[i] = 4.0 * g_x_0_xz_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_z_y[i] = 4.0 * g_x_0_xz_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_z_z[i] = 4.0 * g_x_0_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_x_0_xy_x, g_x_0_xy_y, g_x_0_xy_z, g_x_0_y_0_0_0_x_x, g_x_0_y_0_0_0_x_y, g_x_0_y_0_0_0_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_0_x_x[i] = 4.0 * g_x_0_xy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_x_y[i] = 4.0 * g_x_0_xy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_x_z[i] = 4.0 * g_x_0_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_x_0_0_x, g_x_0_0_y, g_x_0_0_z, g_x_0_y_0_0_0_y_x, g_x_0_y_0_0_0_y_y, g_x_0_y_0_0_0_y_z, g_x_0_yy_x, g_x_0_yy_y, g_x_0_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_0_y_x[i] = -2.0 * g_x_0_0_x[i] * a_exp + 4.0 * g_x_0_yy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_y_y[i] = -2.0 * g_x_0_0_y[i] * a_exp + 4.0 * g_x_0_yy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_y_z[i] = -2.0 * g_x_0_0_z[i] * a_exp + 4.0 * g_x_0_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_x_0_y_0_0_0_z_x, g_x_0_y_0_0_0_z_y, g_x_0_y_0_0_0_z_z, g_x_0_yz_x, g_x_0_yz_y, g_x_0_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_0_z_x[i] = 4.0 * g_x_0_yz_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_z_y[i] = 4.0 * g_x_0_yz_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_z_z[i] = 4.0 * g_x_0_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_x_0_xz_x, g_x_0_xz_y, g_x_0_xz_z, g_x_0_z_0_0_0_x_x, g_x_0_z_0_0_0_x_y, g_x_0_z_0_0_0_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_0_x_x[i] = 4.0 * g_x_0_xz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_x_y[i] = 4.0 * g_x_0_xz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_x_z[i] = 4.0 * g_x_0_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_x_0_yz_x, g_x_0_yz_y, g_x_0_yz_z, g_x_0_z_0_0_0_y_x, g_x_0_z_0_0_0_y_y, g_x_0_z_0_0_0_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_0_y_x[i] = 4.0 * g_x_0_yz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_y_y[i] = 4.0 * g_x_0_yz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_y_z[i] = 4.0 * g_x_0_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_x_0_0_x, g_x_0_0_y, g_x_0_0_z, g_x_0_z_0_0_0_z_x, g_x_0_z_0_0_0_z_y, g_x_0_z_0_0_0_z_z, g_x_0_zz_x, g_x_0_zz_y, g_x_0_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_0_z_x[i] = -2.0 * g_x_0_0_x[i] * a_exp + 4.0 * g_x_0_zz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_z_y[i] = -2.0 * g_x_0_0_y[i] * a_exp + 4.0 * g_x_0_zz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_z_z[i] = -2.0 * g_x_0_0_z[i] * a_exp + 4.0 * g_x_0_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (27-30)

    #pragma omp simd aligned(g_y_0_0_x, g_y_0_0_y, g_y_0_0_z, g_y_0_x_0_0_0_x_x, g_y_0_x_0_0_0_x_y, g_y_0_x_0_0_0_x_z, g_y_0_xx_x, g_y_0_xx_y, g_y_0_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_0_x_x[i] = -2.0 * g_y_0_0_x[i] * a_exp + 4.0 * g_y_0_xx_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_x_y[i] = -2.0 * g_y_0_0_y[i] * a_exp + 4.0 * g_y_0_xx_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_x_z[i] = -2.0 * g_y_0_0_z[i] * a_exp + 4.0 * g_y_0_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (30-33)

    #pragma omp simd aligned(g_y_0_x_0_0_0_y_x, g_y_0_x_0_0_0_y_y, g_y_0_x_0_0_0_y_z, g_y_0_xy_x, g_y_0_xy_y, g_y_0_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_0_y_x[i] = 4.0 * g_y_0_xy_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_y_y[i] = 4.0 * g_y_0_xy_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_y_z[i] = 4.0 * g_y_0_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (33-36)

    #pragma omp simd aligned(g_y_0_x_0_0_0_z_x, g_y_0_x_0_0_0_z_y, g_y_0_x_0_0_0_z_z, g_y_0_xz_x, g_y_0_xz_y, g_y_0_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_0_z_x[i] = 4.0 * g_y_0_xz_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_z_y[i] = 4.0 * g_y_0_xz_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_z_z[i] = 4.0 * g_y_0_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (36-39)

    #pragma omp simd aligned(g_y_0_xy_x, g_y_0_xy_y, g_y_0_xy_z, g_y_0_y_0_0_0_x_x, g_y_0_y_0_0_0_x_y, g_y_0_y_0_0_0_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_0_x_x[i] = 4.0 * g_y_0_xy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_x_y[i] = 4.0 * g_y_0_xy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_x_z[i] = 4.0 * g_y_0_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (39-42)

    #pragma omp simd aligned(g_y_0_0_x, g_y_0_0_y, g_y_0_0_z, g_y_0_y_0_0_0_y_x, g_y_0_y_0_0_0_y_y, g_y_0_y_0_0_0_y_z, g_y_0_yy_x, g_y_0_yy_y, g_y_0_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_0_y_x[i] = -2.0 * g_y_0_0_x[i] * a_exp + 4.0 * g_y_0_yy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_y_y[i] = -2.0 * g_y_0_0_y[i] * a_exp + 4.0 * g_y_0_yy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_y_z[i] = -2.0 * g_y_0_0_z[i] * a_exp + 4.0 * g_y_0_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (42-45)

    #pragma omp simd aligned(g_y_0_y_0_0_0_z_x, g_y_0_y_0_0_0_z_y, g_y_0_y_0_0_0_z_z, g_y_0_yz_x, g_y_0_yz_y, g_y_0_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_0_z_x[i] = 4.0 * g_y_0_yz_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_z_y[i] = 4.0 * g_y_0_yz_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_z_z[i] = 4.0 * g_y_0_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (45-48)

    #pragma omp simd aligned(g_y_0_xz_x, g_y_0_xz_y, g_y_0_xz_z, g_y_0_z_0_0_0_x_x, g_y_0_z_0_0_0_x_y, g_y_0_z_0_0_0_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_0_x_x[i] = 4.0 * g_y_0_xz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_x_y[i] = 4.0 * g_y_0_xz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_x_z[i] = 4.0 * g_y_0_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (48-51)

    #pragma omp simd aligned(g_y_0_yz_x, g_y_0_yz_y, g_y_0_yz_z, g_y_0_z_0_0_0_y_x, g_y_0_z_0_0_0_y_y, g_y_0_z_0_0_0_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_0_y_x[i] = 4.0 * g_y_0_yz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_y_y[i] = 4.0 * g_y_0_yz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_y_z[i] = 4.0 * g_y_0_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (51-54)

    #pragma omp simd aligned(g_y_0_0_x, g_y_0_0_y, g_y_0_0_z, g_y_0_z_0_0_0_z_x, g_y_0_z_0_0_0_z_y, g_y_0_z_0_0_0_z_z, g_y_0_zz_x, g_y_0_zz_y, g_y_0_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_0_z_x[i] = -2.0 * g_y_0_0_x[i] * a_exp + 4.0 * g_y_0_zz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_z_y[i] = -2.0 * g_y_0_0_y[i] * a_exp + 4.0 * g_y_0_zz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_z_z[i] = -2.0 * g_y_0_0_z[i] * a_exp + 4.0 * g_y_0_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (54-57)

    #pragma omp simd aligned(g_z_0_0_x, g_z_0_0_y, g_z_0_0_z, g_z_0_x_0_0_0_x_x, g_z_0_x_0_0_0_x_y, g_z_0_x_0_0_0_x_z, g_z_0_xx_x, g_z_0_xx_y, g_z_0_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_0_x_x[i] = -2.0 * g_z_0_0_x[i] * a_exp + 4.0 * g_z_0_xx_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_x_y[i] = -2.0 * g_z_0_0_y[i] * a_exp + 4.0 * g_z_0_xx_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_x_z[i] = -2.0 * g_z_0_0_z[i] * a_exp + 4.0 * g_z_0_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (57-60)

    #pragma omp simd aligned(g_z_0_x_0_0_0_y_x, g_z_0_x_0_0_0_y_y, g_z_0_x_0_0_0_y_z, g_z_0_xy_x, g_z_0_xy_y, g_z_0_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_0_y_x[i] = 4.0 * g_z_0_xy_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_y_y[i] = 4.0 * g_z_0_xy_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_y_z[i] = 4.0 * g_z_0_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (60-63)

    #pragma omp simd aligned(g_z_0_x_0_0_0_z_x, g_z_0_x_0_0_0_z_y, g_z_0_x_0_0_0_z_z, g_z_0_xz_x, g_z_0_xz_y, g_z_0_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_0_z_x[i] = 4.0 * g_z_0_xz_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_z_y[i] = 4.0 * g_z_0_xz_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_z_z[i] = 4.0 * g_z_0_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (63-66)

    #pragma omp simd aligned(g_z_0_xy_x, g_z_0_xy_y, g_z_0_xy_z, g_z_0_y_0_0_0_x_x, g_z_0_y_0_0_0_x_y, g_z_0_y_0_0_0_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_0_x_x[i] = 4.0 * g_z_0_xy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_x_y[i] = 4.0 * g_z_0_xy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_x_z[i] = 4.0 * g_z_0_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (66-69)

    #pragma omp simd aligned(g_z_0_0_x, g_z_0_0_y, g_z_0_0_z, g_z_0_y_0_0_0_y_x, g_z_0_y_0_0_0_y_y, g_z_0_y_0_0_0_y_z, g_z_0_yy_x, g_z_0_yy_y, g_z_0_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_0_y_x[i] = -2.0 * g_z_0_0_x[i] * a_exp + 4.0 * g_z_0_yy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_y_y[i] = -2.0 * g_z_0_0_y[i] * a_exp + 4.0 * g_z_0_yy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_y_z[i] = -2.0 * g_z_0_0_z[i] * a_exp + 4.0 * g_z_0_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (69-72)

    #pragma omp simd aligned(g_z_0_y_0_0_0_z_x, g_z_0_y_0_0_0_z_y, g_z_0_y_0_0_0_z_z, g_z_0_yz_x, g_z_0_yz_y, g_z_0_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_0_z_x[i] = 4.0 * g_z_0_yz_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_z_y[i] = 4.0 * g_z_0_yz_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_z_z[i] = 4.0 * g_z_0_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (72-75)

    #pragma omp simd aligned(g_z_0_xz_x, g_z_0_xz_y, g_z_0_xz_z, g_z_0_z_0_0_0_x_x, g_z_0_z_0_0_0_x_y, g_z_0_z_0_0_0_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_0_x_x[i] = 4.0 * g_z_0_xz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_x_y[i] = 4.0 * g_z_0_xz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_x_z[i] = 4.0 * g_z_0_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (75-78)

    #pragma omp simd aligned(g_z_0_yz_x, g_z_0_yz_y, g_z_0_yz_z, g_z_0_z_0_0_0_y_x, g_z_0_z_0_0_0_y_y, g_z_0_z_0_0_0_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_0_y_x[i] = 4.0 * g_z_0_yz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_y_y[i] = 4.0 * g_z_0_yz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_y_z[i] = 4.0 * g_z_0_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (78-81)

    #pragma omp simd aligned(g_z_0_0_x, g_z_0_0_y, g_z_0_0_z, g_z_0_z_0_0_0_z_x, g_z_0_z_0_0_0_z_y, g_z_0_z_0_0_0_z_z, g_z_0_zz_x, g_z_0_zz_y, g_z_0_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_0_z_x[i] = -2.0 * g_z_0_0_x[i] * a_exp + 4.0 * g_z_0_zz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_z_y[i] = -2.0 * g_z_0_0_y[i] * a_exp + 4.0 * g_z_0_zz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_z_z[i] = -2.0 * g_z_0_0_z[i] * a_exp + 4.0 * g_z_0_zz_z[i] * a_exp * c_exps[i];
    }
}

} // t4c_geom namespace


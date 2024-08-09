#include "GeomDeriv1010OfScalarForSSSD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1010_sssd_0(CSimdArray<double>& buffer_1010_sssd,
                     const CSimdArray<double>& buffer_pspd,
                     const double a_exp,
                     const double* c_exps) -> void
{
    const auto ndims = buffer_1010_sssd.number_of_columns();

    /// Set up components of auxilary buffer : buffer_pspd

    auto g_x_0_x_xx = buffer_pspd[0];

    auto g_x_0_x_xy = buffer_pspd[1];

    auto g_x_0_x_xz = buffer_pspd[2];

    auto g_x_0_x_yy = buffer_pspd[3];

    auto g_x_0_x_yz = buffer_pspd[4];

    auto g_x_0_x_zz = buffer_pspd[5];

    auto g_x_0_y_xx = buffer_pspd[6];

    auto g_x_0_y_xy = buffer_pspd[7];

    auto g_x_0_y_xz = buffer_pspd[8];

    auto g_x_0_y_yy = buffer_pspd[9];

    auto g_x_0_y_yz = buffer_pspd[10];

    auto g_x_0_y_zz = buffer_pspd[11];

    auto g_x_0_z_xx = buffer_pspd[12];

    auto g_x_0_z_xy = buffer_pspd[13];

    auto g_x_0_z_xz = buffer_pspd[14];

    auto g_x_0_z_yy = buffer_pspd[15];

    auto g_x_0_z_yz = buffer_pspd[16];

    auto g_x_0_z_zz = buffer_pspd[17];

    auto g_y_0_x_xx = buffer_pspd[18];

    auto g_y_0_x_xy = buffer_pspd[19];

    auto g_y_0_x_xz = buffer_pspd[20];

    auto g_y_0_x_yy = buffer_pspd[21];

    auto g_y_0_x_yz = buffer_pspd[22];

    auto g_y_0_x_zz = buffer_pspd[23];

    auto g_y_0_y_xx = buffer_pspd[24];

    auto g_y_0_y_xy = buffer_pspd[25];

    auto g_y_0_y_xz = buffer_pspd[26];

    auto g_y_0_y_yy = buffer_pspd[27];

    auto g_y_0_y_yz = buffer_pspd[28];

    auto g_y_0_y_zz = buffer_pspd[29];

    auto g_y_0_z_xx = buffer_pspd[30];

    auto g_y_0_z_xy = buffer_pspd[31];

    auto g_y_0_z_xz = buffer_pspd[32];

    auto g_y_0_z_yy = buffer_pspd[33];

    auto g_y_0_z_yz = buffer_pspd[34];

    auto g_y_0_z_zz = buffer_pspd[35];

    auto g_z_0_x_xx = buffer_pspd[36];

    auto g_z_0_x_xy = buffer_pspd[37];

    auto g_z_0_x_xz = buffer_pspd[38];

    auto g_z_0_x_yy = buffer_pspd[39];

    auto g_z_0_x_yz = buffer_pspd[40];

    auto g_z_0_x_zz = buffer_pspd[41];

    auto g_z_0_y_xx = buffer_pspd[42];

    auto g_z_0_y_xy = buffer_pspd[43];

    auto g_z_0_y_xz = buffer_pspd[44];

    auto g_z_0_y_yy = buffer_pspd[45];

    auto g_z_0_y_yz = buffer_pspd[46];

    auto g_z_0_y_zz = buffer_pspd[47];

    auto g_z_0_z_xx = buffer_pspd[48];

    auto g_z_0_z_xy = buffer_pspd[49];

    auto g_z_0_z_xz = buffer_pspd[50];

    auto g_z_0_z_yy = buffer_pspd[51];

    auto g_z_0_z_yz = buffer_pspd[52];

    auto g_z_0_z_zz = buffer_pspd[53];

    /// Set up components of integrals buffer : buffer_1010_sssd

    auto g_x_0_x_0_0_0_0_xx = buffer_1010_sssd[0];

    auto g_x_0_x_0_0_0_0_xy = buffer_1010_sssd[1];

    auto g_x_0_x_0_0_0_0_xz = buffer_1010_sssd[2];

    auto g_x_0_x_0_0_0_0_yy = buffer_1010_sssd[3];

    auto g_x_0_x_0_0_0_0_yz = buffer_1010_sssd[4];

    auto g_x_0_x_0_0_0_0_zz = buffer_1010_sssd[5];

    auto g_x_0_y_0_0_0_0_xx = buffer_1010_sssd[6];

    auto g_x_0_y_0_0_0_0_xy = buffer_1010_sssd[7];

    auto g_x_0_y_0_0_0_0_xz = buffer_1010_sssd[8];

    auto g_x_0_y_0_0_0_0_yy = buffer_1010_sssd[9];

    auto g_x_0_y_0_0_0_0_yz = buffer_1010_sssd[10];

    auto g_x_0_y_0_0_0_0_zz = buffer_1010_sssd[11];

    auto g_x_0_z_0_0_0_0_xx = buffer_1010_sssd[12];

    auto g_x_0_z_0_0_0_0_xy = buffer_1010_sssd[13];

    auto g_x_0_z_0_0_0_0_xz = buffer_1010_sssd[14];

    auto g_x_0_z_0_0_0_0_yy = buffer_1010_sssd[15];

    auto g_x_0_z_0_0_0_0_yz = buffer_1010_sssd[16];

    auto g_x_0_z_0_0_0_0_zz = buffer_1010_sssd[17];

    auto g_y_0_x_0_0_0_0_xx = buffer_1010_sssd[18];

    auto g_y_0_x_0_0_0_0_xy = buffer_1010_sssd[19];

    auto g_y_0_x_0_0_0_0_xz = buffer_1010_sssd[20];

    auto g_y_0_x_0_0_0_0_yy = buffer_1010_sssd[21];

    auto g_y_0_x_0_0_0_0_yz = buffer_1010_sssd[22];

    auto g_y_0_x_0_0_0_0_zz = buffer_1010_sssd[23];

    auto g_y_0_y_0_0_0_0_xx = buffer_1010_sssd[24];

    auto g_y_0_y_0_0_0_0_xy = buffer_1010_sssd[25];

    auto g_y_0_y_0_0_0_0_xz = buffer_1010_sssd[26];

    auto g_y_0_y_0_0_0_0_yy = buffer_1010_sssd[27];

    auto g_y_0_y_0_0_0_0_yz = buffer_1010_sssd[28];

    auto g_y_0_y_0_0_0_0_zz = buffer_1010_sssd[29];

    auto g_y_0_z_0_0_0_0_xx = buffer_1010_sssd[30];

    auto g_y_0_z_0_0_0_0_xy = buffer_1010_sssd[31];

    auto g_y_0_z_0_0_0_0_xz = buffer_1010_sssd[32];

    auto g_y_0_z_0_0_0_0_yy = buffer_1010_sssd[33];

    auto g_y_0_z_0_0_0_0_yz = buffer_1010_sssd[34];

    auto g_y_0_z_0_0_0_0_zz = buffer_1010_sssd[35];

    auto g_z_0_x_0_0_0_0_xx = buffer_1010_sssd[36];

    auto g_z_0_x_0_0_0_0_xy = buffer_1010_sssd[37];

    auto g_z_0_x_0_0_0_0_xz = buffer_1010_sssd[38];

    auto g_z_0_x_0_0_0_0_yy = buffer_1010_sssd[39];

    auto g_z_0_x_0_0_0_0_yz = buffer_1010_sssd[40];

    auto g_z_0_x_0_0_0_0_zz = buffer_1010_sssd[41];

    auto g_z_0_y_0_0_0_0_xx = buffer_1010_sssd[42];

    auto g_z_0_y_0_0_0_0_xy = buffer_1010_sssd[43];

    auto g_z_0_y_0_0_0_0_xz = buffer_1010_sssd[44];

    auto g_z_0_y_0_0_0_0_yy = buffer_1010_sssd[45];

    auto g_z_0_y_0_0_0_0_yz = buffer_1010_sssd[46];

    auto g_z_0_y_0_0_0_0_zz = buffer_1010_sssd[47];

    auto g_z_0_z_0_0_0_0_xx = buffer_1010_sssd[48];

    auto g_z_0_z_0_0_0_0_xy = buffer_1010_sssd[49];

    auto g_z_0_z_0_0_0_0_xz = buffer_1010_sssd[50];

    auto g_z_0_z_0_0_0_0_yy = buffer_1010_sssd[51];

    auto g_z_0_z_0_0_0_0_yz = buffer_1010_sssd[52];

    auto g_z_0_z_0_0_0_0_zz = buffer_1010_sssd[53];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_0_x_0_0_0_0_xx, g_x_0_x_0_0_0_0_xy, g_x_0_x_0_0_0_0_xz, g_x_0_x_0_0_0_0_yy, g_x_0_x_0_0_0_0_yz, g_x_0_x_0_0_0_0_zz, g_x_0_x_xx, g_x_0_x_xy, g_x_0_x_xz, g_x_0_x_yy, g_x_0_x_yz, g_x_0_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_0_0_xx[i] = 4.0 * g_x_0_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_0_xy[i] = 4.0 * g_x_0_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_0_xz[i] = 4.0 * g_x_0_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_0_yy[i] = 4.0 * g_x_0_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_0_yz[i] = 4.0 * g_x_0_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_0_zz[i] = 4.0 * g_x_0_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_0_y_0_0_0_0_xx, g_x_0_y_0_0_0_0_xy, g_x_0_y_0_0_0_0_xz, g_x_0_y_0_0_0_0_yy, g_x_0_y_0_0_0_0_yz, g_x_0_y_0_0_0_0_zz, g_x_0_y_xx, g_x_0_y_xy, g_x_0_y_xz, g_x_0_y_yy, g_x_0_y_yz, g_x_0_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_0_0_xx[i] = 4.0 * g_x_0_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_0_xy[i] = 4.0 * g_x_0_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_0_xz[i] = 4.0 * g_x_0_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_0_yy[i] = 4.0 * g_x_0_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_0_yz[i] = 4.0 * g_x_0_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_0_zz[i] = 4.0 * g_x_0_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_0_z_0_0_0_0_xx, g_x_0_z_0_0_0_0_xy, g_x_0_z_0_0_0_0_xz, g_x_0_z_0_0_0_0_yy, g_x_0_z_0_0_0_0_yz, g_x_0_z_0_0_0_0_zz, g_x_0_z_xx, g_x_0_z_xy, g_x_0_z_xz, g_x_0_z_yy, g_x_0_z_yz, g_x_0_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_0_0_xx[i] = 4.0 * g_x_0_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_0_xy[i] = 4.0 * g_x_0_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_0_xz[i] = 4.0 * g_x_0_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_0_yy[i] = 4.0 * g_x_0_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_0_yz[i] = 4.0 * g_x_0_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_0_zz[i] = 4.0 * g_x_0_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_y_0_x_0_0_0_0_xx, g_y_0_x_0_0_0_0_xy, g_y_0_x_0_0_0_0_xz, g_y_0_x_0_0_0_0_yy, g_y_0_x_0_0_0_0_yz, g_y_0_x_0_0_0_0_zz, g_y_0_x_xx, g_y_0_x_xy, g_y_0_x_xz, g_y_0_x_yy, g_y_0_x_yz, g_y_0_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_0_0_xx[i] = 4.0 * g_y_0_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_0_xy[i] = 4.0 * g_y_0_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_0_xz[i] = 4.0 * g_y_0_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_0_yy[i] = 4.0 * g_y_0_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_0_yz[i] = 4.0 * g_y_0_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_0_zz[i] = 4.0 * g_y_0_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_y_0_y_0_0_0_0_xx, g_y_0_y_0_0_0_0_xy, g_y_0_y_0_0_0_0_xz, g_y_0_y_0_0_0_0_yy, g_y_0_y_0_0_0_0_yz, g_y_0_y_0_0_0_0_zz, g_y_0_y_xx, g_y_0_y_xy, g_y_0_y_xz, g_y_0_y_yy, g_y_0_y_yz, g_y_0_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_0_0_xx[i] = 4.0 * g_y_0_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_0_xy[i] = 4.0 * g_y_0_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_0_xz[i] = 4.0 * g_y_0_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_0_yy[i] = 4.0 * g_y_0_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_0_yz[i] = 4.0 * g_y_0_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_0_zz[i] = 4.0 * g_y_0_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_y_0_z_0_0_0_0_xx, g_y_0_z_0_0_0_0_xy, g_y_0_z_0_0_0_0_xz, g_y_0_z_0_0_0_0_yy, g_y_0_z_0_0_0_0_yz, g_y_0_z_0_0_0_0_zz, g_y_0_z_xx, g_y_0_z_xy, g_y_0_z_xz, g_y_0_z_yy, g_y_0_z_yz, g_y_0_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_0_0_xx[i] = 4.0 * g_y_0_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_0_xy[i] = 4.0 * g_y_0_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_0_xz[i] = 4.0 * g_y_0_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_0_yy[i] = 4.0 * g_y_0_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_0_yz[i] = 4.0 * g_y_0_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_0_zz[i] = 4.0 * g_y_0_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_z_0_x_0_0_0_0_xx, g_z_0_x_0_0_0_0_xy, g_z_0_x_0_0_0_0_xz, g_z_0_x_0_0_0_0_yy, g_z_0_x_0_0_0_0_yz, g_z_0_x_0_0_0_0_zz, g_z_0_x_xx, g_z_0_x_xy, g_z_0_x_xz, g_z_0_x_yy, g_z_0_x_yz, g_z_0_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_0_0_xx[i] = 4.0 * g_z_0_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_0_xy[i] = 4.0 * g_z_0_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_0_xz[i] = 4.0 * g_z_0_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_0_yy[i] = 4.0 * g_z_0_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_0_yz[i] = 4.0 * g_z_0_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_0_zz[i] = 4.0 * g_z_0_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_z_0_y_0_0_0_0_xx, g_z_0_y_0_0_0_0_xy, g_z_0_y_0_0_0_0_xz, g_z_0_y_0_0_0_0_yy, g_z_0_y_0_0_0_0_yz, g_z_0_y_0_0_0_0_zz, g_z_0_y_xx, g_z_0_y_xy, g_z_0_y_xz, g_z_0_y_yy, g_z_0_y_yz, g_z_0_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_0_0_xx[i] = 4.0 * g_z_0_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_0_xy[i] = 4.0 * g_z_0_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_0_xz[i] = 4.0 * g_z_0_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_0_yy[i] = 4.0 * g_z_0_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_0_yz[i] = 4.0 * g_z_0_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_0_zz[i] = 4.0 * g_z_0_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_z_0_z_0_0_0_0_xx, g_z_0_z_0_0_0_0_xy, g_z_0_z_0_0_0_0_xz, g_z_0_z_0_0_0_0_yy, g_z_0_z_0_0_0_0_yz, g_z_0_z_0_0_0_0_zz, g_z_0_z_xx, g_z_0_z_xy, g_z_0_z_xz, g_z_0_z_yy, g_z_0_z_yz, g_z_0_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_0_0_xx[i] = 4.0 * g_z_0_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_0_xy[i] = 4.0 * g_z_0_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_0_xz[i] = 4.0 * g_z_0_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_0_yy[i] = 4.0 * g_z_0_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_0_yz[i] = 4.0 * g_z_0_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_0_zz[i] = 4.0 * g_z_0_z_zz[i] * a_exp * c_exps[i];
    }
}

} // t4c_geom namespace


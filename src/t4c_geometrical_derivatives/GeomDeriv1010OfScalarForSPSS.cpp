#include "GeomDeriv1010OfScalarForSPSS.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1010_spss_0(CSimdArray<double>& buffer_1010_spss,
                     const CSimdArray<double>& buffer_ppps,
                     const double a_exp,
                     const double* c_exps) -> void
{
    const auto ndims = buffer_1010_spss.number_of_columns();

    /// Set up components of auxilary buffer : buffer_ppps

    auto g_x_x_x_0 = buffer_ppps[0];

    auto g_x_x_y_0 = buffer_ppps[1];

    auto g_x_x_z_0 = buffer_ppps[2];

    auto g_x_y_x_0 = buffer_ppps[3];

    auto g_x_y_y_0 = buffer_ppps[4];

    auto g_x_y_z_0 = buffer_ppps[5];

    auto g_x_z_x_0 = buffer_ppps[6];

    auto g_x_z_y_0 = buffer_ppps[7];

    auto g_x_z_z_0 = buffer_ppps[8];

    auto g_y_x_x_0 = buffer_ppps[9];

    auto g_y_x_y_0 = buffer_ppps[10];

    auto g_y_x_z_0 = buffer_ppps[11];

    auto g_y_y_x_0 = buffer_ppps[12];

    auto g_y_y_y_0 = buffer_ppps[13];

    auto g_y_y_z_0 = buffer_ppps[14];

    auto g_y_z_x_0 = buffer_ppps[15];

    auto g_y_z_y_0 = buffer_ppps[16];

    auto g_y_z_z_0 = buffer_ppps[17];

    auto g_z_x_x_0 = buffer_ppps[18];

    auto g_z_x_y_0 = buffer_ppps[19];

    auto g_z_x_z_0 = buffer_ppps[20];

    auto g_z_y_x_0 = buffer_ppps[21];

    auto g_z_y_y_0 = buffer_ppps[22];

    auto g_z_y_z_0 = buffer_ppps[23];

    auto g_z_z_x_0 = buffer_ppps[24];

    auto g_z_z_y_0 = buffer_ppps[25];

    auto g_z_z_z_0 = buffer_ppps[26];

    /// Set up components of integrals buffer : buffer_1010_spss

    auto g_x_0_x_0_0_x_0_0 = buffer_1010_spss[0];

    auto g_x_0_x_0_0_y_0_0 = buffer_1010_spss[1];

    auto g_x_0_x_0_0_z_0_0 = buffer_1010_spss[2];

    auto g_x_0_y_0_0_x_0_0 = buffer_1010_spss[3];

    auto g_x_0_y_0_0_y_0_0 = buffer_1010_spss[4];

    auto g_x_0_y_0_0_z_0_0 = buffer_1010_spss[5];

    auto g_x_0_z_0_0_x_0_0 = buffer_1010_spss[6];

    auto g_x_0_z_0_0_y_0_0 = buffer_1010_spss[7];

    auto g_x_0_z_0_0_z_0_0 = buffer_1010_spss[8];

    auto g_y_0_x_0_0_x_0_0 = buffer_1010_spss[9];

    auto g_y_0_x_0_0_y_0_0 = buffer_1010_spss[10];

    auto g_y_0_x_0_0_z_0_0 = buffer_1010_spss[11];

    auto g_y_0_y_0_0_x_0_0 = buffer_1010_spss[12];

    auto g_y_0_y_0_0_y_0_0 = buffer_1010_spss[13];

    auto g_y_0_y_0_0_z_0_0 = buffer_1010_spss[14];

    auto g_y_0_z_0_0_x_0_0 = buffer_1010_spss[15];

    auto g_y_0_z_0_0_y_0_0 = buffer_1010_spss[16];

    auto g_y_0_z_0_0_z_0_0 = buffer_1010_spss[17];

    auto g_z_0_x_0_0_x_0_0 = buffer_1010_spss[18];

    auto g_z_0_x_0_0_y_0_0 = buffer_1010_spss[19];

    auto g_z_0_x_0_0_z_0_0 = buffer_1010_spss[20];

    auto g_z_0_y_0_0_x_0_0 = buffer_1010_spss[21];

    auto g_z_0_y_0_0_y_0_0 = buffer_1010_spss[22];

    auto g_z_0_y_0_0_z_0_0 = buffer_1010_spss[23];

    auto g_z_0_z_0_0_x_0_0 = buffer_1010_spss[24];

    auto g_z_0_z_0_0_y_0_0 = buffer_1010_spss[25];

    auto g_z_0_z_0_0_z_0_0 = buffer_1010_spss[26];

    // integrals block (0-3)

    #pragma omp simd aligned(g_x_0_x_0_0_x_0_0, g_x_0_x_0_0_y_0_0, g_x_0_x_0_0_z_0_0, g_x_x_x_0, g_x_y_x_0, g_x_z_x_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_x_0_0[i] = 4.0 * g_x_x_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_0_0[i] = 4.0 * g_x_y_x_0[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_0_0[i] = 4.0 * g_x_z_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_x_0_y_0_0_x_0_0, g_x_0_y_0_0_y_0_0, g_x_0_y_0_0_z_0_0, g_x_x_y_0, g_x_y_y_0, g_x_z_y_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_x_0_0[i] = 4.0 * g_x_x_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_0_0[i] = 4.0 * g_x_y_y_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_0_0[i] = 4.0 * g_x_z_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_x_0_z_0_0_x_0_0, g_x_0_z_0_0_y_0_0, g_x_0_z_0_0_z_0_0, g_x_x_z_0, g_x_y_z_0, g_x_z_z_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_x_0_0[i] = 4.0 * g_x_x_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_0_0[i] = 4.0 * g_x_y_z_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_0_0[i] = 4.0 * g_x_z_z_0[i] * a_exp * c_exps[i];
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_y_0_x_0_0_x_0_0, g_y_0_x_0_0_y_0_0, g_y_0_x_0_0_z_0_0, g_y_x_x_0, g_y_y_x_0, g_y_z_x_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_x_0_0[i] = 4.0 * g_y_x_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_0_0[i] = 4.0 * g_y_y_x_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_0_0[i] = 4.0 * g_y_z_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_y_0_y_0_0_x_0_0, g_y_0_y_0_0_y_0_0, g_y_0_y_0_0_z_0_0, g_y_x_y_0, g_y_y_y_0, g_y_z_y_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_x_0_0[i] = 4.0 * g_y_x_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_0_0[i] = 4.0 * g_y_y_y_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_0_0[i] = 4.0 * g_y_z_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_y_0_z_0_0_x_0_0, g_y_0_z_0_0_y_0_0, g_y_0_z_0_0_z_0_0, g_y_x_z_0, g_y_y_z_0, g_y_z_z_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_x_0_0[i] = 4.0 * g_y_x_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_0_0[i] = 4.0 * g_y_y_z_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_0_0[i] = 4.0 * g_y_z_z_0[i] * a_exp * c_exps[i];
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_z_0_x_0_0_x_0_0, g_z_0_x_0_0_y_0_0, g_z_0_x_0_0_z_0_0, g_z_x_x_0, g_z_y_x_0, g_z_z_x_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_x_0_0[i] = 4.0 * g_z_x_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_0_0[i] = 4.0 * g_z_y_x_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_0_0[i] = 4.0 * g_z_z_x_0[i] * a_exp * c_exps[i];
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_z_0_y_0_0_x_0_0, g_z_0_y_0_0_y_0_0, g_z_0_y_0_0_z_0_0, g_z_x_y_0, g_z_y_y_0, g_z_z_y_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_x_0_0[i] = 4.0 * g_z_x_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_0_0[i] = 4.0 * g_z_y_y_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_0_0[i] = 4.0 * g_z_z_y_0[i] * a_exp * c_exps[i];
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_z_0_z_0_0_x_0_0, g_z_0_z_0_0_y_0_0, g_z_0_z_0_0_z_0_0, g_z_x_z_0, g_z_y_z_0, g_z_z_z_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_x_0_0[i] = 4.0 * g_z_x_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_0_0[i] = 4.0 * g_z_y_z_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_0_0[i] = 4.0 * g_z_z_z_0[i] * a_exp * c_exps[i];
    }
}

} // t4c_geom namespace


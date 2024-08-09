#include "GeomDeriv2000OfScalarForSSSS.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom2000_ssss_0(CSimdArray<double>& buffer_2000_ssss,
                     const CSimdArray<double>& buffer_ssss,
                     const CSimdArray<double>& buffer_dsss,
                     const double a_exp) -> void
{
    const auto ndims = buffer_2000_ssss.number_of_columns();

    /// Set up components of auxilary buffer : buffer_ssss

    auto g_0_0_0_0 = buffer_ssss[0];

    /// Set up components of auxilary buffer : buffer_dsss

    auto g_xx_0_0_0 = buffer_dsss[0];

    auto g_xy_0_0_0 = buffer_dsss[1];

    auto g_xz_0_0_0 = buffer_dsss[2];

    auto g_yy_0_0_0 = buffer_dsss[3];

    auto g_yz_0_0_0 = buffer_dsss[4];

    auto g_zz_0_0_0 = buffer_dsss[5];

    /// Set up components of integrals buffer : buffer_2000_ssss

    auto g_xx_0_0_0_0_0_0_0 = buffer_2000_ssss[0];

    auto g_xy_0_0_0_0_0_0_0 = buffer_2000_ssss[1];

    auto g_xz_0_0_0_0_0_0_0 = buffer_2000_ssss[2];

    auto g_yy_0_0_0_0_0_0_0 = buffer_2000_ssss[3];

    auto g_yz_0_0_0_0_0_0_0 = buffer_2000_ssss[4];

    auto g_zz_0_0_0_0_0_0_0 = buffer_2000_ssss[5];

    // integrals block (0-6)

    #pragma omp simd aligned(g_0_0_0_0, g_xx_0_0_0, g_xx_0_0_0_0_0_0_0, g_xy_0_0_0, g_xy_0_0_0_0_0_0_0, g_xz_0_0_0, g_xz_0_0_0_0_0_0_0, g_yy_0_0_0, g_yy_0_0_0_0_0_0_0, g_yz_0_0_0, g_yz_0_0_0_0_0_0_0, g_zz_0_0_0, g_zz_0_0_0_0_0_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_0_0_0[i] = -2.0 * g_0_0_0_0[i] * a_exp + 4.0 * g_xx_0_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_0_0[i] = 4.0 * g_xy_0_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_0_0[i] = 4.0 * g_xz_0_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_0_0[i] = -2.0 * g_0_0_0_0[i] * a_exp + 4.0 * g_yy_0_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_0_0[i] = 4.0 * g_yz_0_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_0_0[i] = -2.0 * g_0_0_0_0[i] * a_exp + 4.0 * g_zz_0_0_0[i] * a_exp * a_exp;
    }
}

} // t4c_geom namespace


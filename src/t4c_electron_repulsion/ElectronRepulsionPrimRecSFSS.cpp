#include "ElectronRepulsionPrimRecSFSS.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_sfss(CSimdArray<double>& prim_buffer_0_sfss,
                                  const CSimdArray<double>& prim_buffer_0_spss,
                                  const CSimdArray<double>& prim_buffer_1_spss,
                                  const CSimdArray<double>& prim_buffer_0_sdss,
                                  const CSimdArray<double>& prim_buffer_1_sdss,
                                  const double pb_x,
                                  const double pb_y,
                                  const double pb_z,
                                  const double* wp_x,
                                  const double* wp_y,
                                  const double* wp_z,
                                  const double a_exp,
                                  const double b_exp,
                                  const double* c_exps,
                                  const double* d_exps) -> void
{
    const auto ndims = prim_buffer_0_sfss.number_of_columns();

    /// Set up components of auxilary buffer : prim_buffer_0_spss

    auto g_0_x_0_0_0 = prim_buffer_0_spss[0];

    auto g_0_y_0_0_0 = prim_buffer_0_spss[1];

    auto g_0_z_0_0_0 = prim_buffer_0_spss[2];

    /// Set up components of auxilary buffer : prim_buffer_1_spss

    auto g_0_x_0_0_1 = prim_buffer_1_spss[0];

    auto g_0_y_0_0_1 = prim_buffer_1_spss[1];

    auto g_0_z_0_0_1 = prim_buffer_1_spss[2];

    /// Set up components of auxilary buffer : prim_buffer_0_sdss

    auto g_0_xx_0_0_0 = prim_buffer_0_sdss[0];

    auto g_0_yy_0_0_0 = prim_buffer_0_sdss[3];

    auto g_0_yz_0_0_0 = prim_buffer_0_sdss[4];

    auto g_0_zz_0_0_0 = prim_buffer_0_sdss[5];

    /// Set up components of auxilary buffer : prim_buffer_1_sdss

    auto g_0_xx_0_0_1 = prim_buffer_1_sdss[0];

    auto g_0_yy_0_0_1 = prim_buffer_1_sdss[3];

    auto g_0_yz_0_0_1 = prim_buffer_1_sdss[4];

    auto g_0_zz_0_0_1 = prim_buffer_1_sdss[5];

    /// Set up components of targeted buffer : prim_buffer_0_sfss

    auto g_0_xxx_0_0_0 = prim_buffer_0_sfss[0];

    auto g_0_xxy_0_0_0 = prim_buffer_0_sfss[1];

    auto g_0_xxz_0_0_0 = prim_buffer_0_sfss[2];

    auto g_0_xyy_0_0_0 = prim_buffer_0_sfss[3];

    auto g_0_xyz_0_0_0 = prim_buffer_0_sfss[4];

    auto g_0_xzz_0_0_0 = prim_buffer_0_sfss[5];

    auto g_0_yyy_0_0_0 = prim_buffer_0_sfss[6];

    auto g_0_yyz_0_0_0 = prim_buffer_0_sfss[7];

    auto g_0_yzz_0_0_0 = prim_buffer_0_sfss[8];

    auto g_0_zzz_0_0_0 = prim_buffer_0_sfss[9];

    #pragma omp simd aligned(g_0_x_0_0_0, g_0_x_0_0_1, g_0_xx_0_0_0, g_0_xx_0_0_1, g_0_xxx_0_0_0, g_0_xxy_0_0_0, g_0_xxz_0_0_0, g_0_xyy_0_0_0, g_0_xyz_0_0_0, g_0_xzz_0_0_0, g_0_y_0_0_0, g_0_y_0_0_1, g_0_yy_0_0_0, g_0_yy_0_0_1, g_0_yyy_0_0_0, g_0_yyz_0_0_0, g_0_yz_0_0_0, g_0_yz_0_0_1, g_0_yzz_0_0_0, g_0_z_0_0_0, g_0_z_0_0_1, g_0_zz_0_0_0, g_0_zz_0_0_1, g_0_zzz_0_0_0, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxx_0_0_0[i] = 2.0 * g_0_x_0_0_0[i] * fi_ab_0 - 2.0 * g_0_x_0_0_1[i] * fti_ab_0 + g_0_xx_0_0_0[i] * pb_x + g_0_xx_0_0_1[i] * wp_x[i];

        g_0_xxy_0_0_0[i] = g_0_xx_0_0_0[i] * pb_y + g_0_xx_0_0_1[i] * wp_y[i];

        g_0_xxz_0_0_0[i] = g_0_xx_0_0_0[i] * pb_z + g_0_xx_0_0_1[i] * wp_z[i];

        g_0_xyy_0_0_0[i] = g_0_yy_0_0_0[i] * pb_x + g_0_yy_0_0_1[i] * wp_x[i];

        g_0_xyz_0_0_0[i] = g_0_yz_0_0_0[i] * pb_x + g_0_yz_0_0_1[i] * wp_x[i];

        g_0_xzz_0_0_0[i] = g_0_zz_0_0_0[i] * pb_x + g_0_zz_0_0_1[i] * wp_x[i];

        g_0_yyy_0_0_0[i] = 2.0 * g_0_y_0_0_0[i] * fi_ab_0 - 2.0 * g_0_y_0_0_1[i] * fti_ab_0 + g_0_yy_0_0_0[i] * pb_y + g_0_yy_0_0_1[i] * wp_y[i];

        g_0_yyz_0_0_0[i] = g_0_yy_0_0_0[i] * pb_z + g_0_yy_0_0_1[i] * wp_z[i];

        g_0_yzz_0_0_0[i] = g_0_zz_0_0_0[i] * pb_y + g_0_zz_0_0_1[i] * wp_y[i];

        g_0_zzz_0_0_0[i] = 2.0 * g_0_z_0_0_0[i] * fi_ab_0 - 2.0 * g_0_z_0_0_1[i] * fti_ab_0 + g_0_zz_0_0_0[i] * pb_z + g_0_zz_0_0_1[i] * wp_z[i];
    }
}

} // erirec namespace


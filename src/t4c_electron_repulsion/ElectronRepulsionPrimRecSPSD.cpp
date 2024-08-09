#include "ElectronRepulsionPrimRecSPSD.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_spsd(CSimdArray<double>& prim_buffer_0_spsd,
                                  const CSimdArray<double>& prim_buffer_1_sssp,
                                  const CSimdArray<double>& prim_buffer_0_sssd,
                                  const CSimdArray<double>& prim_buffer_1_sssd,
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
    const auto ndims = prim_buffer_0_spsd.number_of_columns();

    /// Set up components of auxilary buffer : prim_buffer_1_sssp

    auto g_0_0_0_x_1 = prim_buffer_1_sssp[0];

    auto g_0_0_0_y_1 = prim_buffer_1_sssp[1];

    auto g_0_0_0_z_1 = prim_buffer_1_sssp[2];

    /// Set up components of auxilary buffer : prim_buffer_0_sssd

    auto g_0_0_0_xx_0 = prim_buffer_0_sssd[0];

    auto g_0_0_0_xy_0 = prim_buffer_0_sssd[1];

    auto g_0_0_0_xz_0 = prim_buffer_0_sssd[2];

    auto g_0_0_0_yy_0 = prim_buffer_0_sssd[3];

    auto g_0_0_0_yz_0 = prim_buffer_0_sssd[4];

    auto g_0_0_0_zz_0 = prim_buffer_0_sssd[5];

    /// Set up components of auxilary buffer : prim_buffer_1_sssd

    auto g_0_0_0_xx_1 = prim_buffer_1_sssd[0];

    auto g_0_0_0_xy_1 = prim_buffer_1_sssd[1];

    auto g_0_0_0_xz_1 = prim_buffer_1_sssd[2];

    auto g_0_0_0_yy_1 = prim_buffer_1_sssd[3];

    auto g_0_0_0_yz_1 = prim_buffer_1_sssd[4];

    auto g_0_0_0_zz_1 = prim_buffer_1_sssd[5];

    /// Set up 0-6 components of targeted buffer : prim_buffer_0_spsd

    auto g_0_x_0_xx_0 = prim_buffer_0_spsd[0];

    auto g_0_x_0_xy_0 = prim_buffer_0_spsd[1];

    auto g_0_x_0_xz_0 = prim_buffer_0_spsd[2];

    auto g_0_x_0_yy_0 = prim_buffer_0_spsd[3];

    auto g_0_x_0_yz_0 = prim_buffer_0_spsd[4];

    auto g_0_x_0_zz_0 = prim_buffer_0_spsd[5];

    #pragma omp simd aligned(g_0_0_0_x_1, g_0_0_0_xx_0, g_0_0_0_xx_1, g_0_0_0_xy_0, g_0_0_0_xy_1, g_0_0_0_xz_0, g_0_0_0_xz_1, g_0_0_0_y_1, g_0_0_0_yy_0, g_0_0_0_yy_1, g_0_0_0_yz_0, g_0_0_0_yz_1, g_0_0_0_z_1, g_0_0_0_zz_0, g_0_0_0_zz_1, g_0_x_0_xx_0, g_0_x_0_xy_0, g_0_x_0_xz_0, g_0_x_0_yy_0, g_0_x_0_yz_0, g_0_x_0_zz_0, wp_x  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_x_0_xx_0[i] = 2.0 * g_0_0_0_x_1[i] * fi_abcd_0 + g_0_0_0_xx_0[i] * pb_x + g_0_0_0_xx_1[i] * wp_x[i];

        g_0_x_0_xy_0[i] = g_0_0_0_y_1[i] * fi_abcd_0 + g_0_0_0_xy_0[i] * pb_x + g_0_0_0_xy_1[i] * wp_x[i];

        g_0_x_0_xz_0[i] = g_0_0_0_z_1[i] * fi_abcd_0 + g_0_0_0_xz_0[i] * pb_x + g_0_0_0_xz_1[i] * wp_x[i];

        g_0_x_0_yy_0[i] = g_0_0_0_yy_0[i] * pb_x + g_0_0_0_yy_1[i] * wp_x[i];

        g_0_x_0_yz_0[i] = g_0_0_0_yz_0[i] * pb_x + g_0_0_0_yz_1[i] * wp_x[i];

        g_0_x_0_zz_0[i] = g_0_0_0_zz_0[i] * pb_x + g_0_0_0_zz_1[i] * wp_x[i];
    }

    /// Set up 6-12 components of targeted buffer : prim_buffer_0_spsd

    auto g_0_y_0_xx_0 = prim_buffer_0_spsd[6];

    auto g_0_y_0_xy_0 = prim_buffer_0_spsd[7];

    auto g_0_y_0_xz_0 = prim_buffer_0_spsd[8];

    auto g_0_y_0_yy_0 = prim_buffer_0_spsd[9];

    auto g_0_y_0_yz_0 = prim_buffer_0_spsd[10];

    auto g_0_y_0_zz_0 = prim_buffer_0_spsd[11];

    #pragma omp simd aligned(g_0_0_0_x_1, g_0_0_0_xx_0, g_0_0_0_xx_1, g_0_0_0_xy_0, g_0_0_0_xy_1, g_0_0_0_xz_0, g_0_0_0_xz_1, g_0_0_0_y_1, g_0_0_0_yy_0, g_0_0_0_yy_1, g_0_0_0_yz_0, g_0_0_0_yz_1, g_0_0_0_z_1, g_0_0_0_zz_0, g_0_0_0_zz_1, g_0_y_0_xx_0, g_0_y_0_xy_0, g_0_y_0_xz_0, g_0_y_0_yy_0, g_0_y_0_yz_0, g_0_y_0_zz_0, wp_y  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_y_0_xx_0[i] = g_0_0_0_xx_0[i] * pb_y + g_0_0_0_xx_1[i] * wp_y[i];

        g_0_y_0_xy_0[i] = g_0_0_0_x_1[i] * fi_abcd_0 + g_0_0_0_xy_0[i] * pb_y + g_0_0_0_xy_1[i] * wp_y[i];

        g_0_y_0_xz_0[i] = g_0_0_0_xz_0[i] * pb_y + g_0_0_0_xz_1[i] * wp_y[i];

        g_0_y_0_yy_0[i] = 2.0 * g_0_0_0_y_1[i] * fi_abcd_0 + g_0_0_0_yy_0[i] * pb_y + g_0_0_0_yy_1[i] * wp_y[i];

        g_0_y_0_yz_0[i] = g_0_0_0_z_1[i] * fi_abcd_0 + g_0_0_0_yz_0[i] * pb_y + g_0_0_0_yz_1[i] * wp_y[i];

        g_0_y_0_zz_0[i] = g_0_0_0_zz_0[i] * pb_y + g_0_0_0_zz_1[i] * wp_y[i];
    }

    /// Set up 12-18 components of targeted buffer : prim_buffer_0_spsd

    auto g_0_z_0_xx_0 = prim_buffer_0_spsd[12];

    auto g_0_z_0_xy_0 = prim_buffer_0_spsd[13];

    auto g_0_z_0_xz_0 = prim_buffer_0_spsd[14];

    auto g_0_z_0_yy_0 = prim_buffer_0_spsd[15];

    auto g_0_z_0_yz_0 = prim_buffer_0_spsd[16];

    auto g_0_z_0_zz_0 = prim_buffer_0_spsd[17];

    #pragma omp simd aligned(g_0_0_0_x_1, g_0_0_0_xx_0, g_0_0_0_xx_1, g_0_0_0_xy_0, g_0_0_0_xy_1, g_0_0_0_xz_0, g_0_0_0_xz_1, g_0_0_0_y_1, g_0_0_0_yy_0, g_0_0_0_yy_1, g_0_0_0_yz_0, g_0_0_0_yz_1, g_0_0_0_z_1, g_0_0_0_zz_0, g_0_0_0_zz_1, g_0_z_0_xx_0, g_0_z_0_xy_0, g_0_z_0_xz_0, g_0_z_0_yy_0, g_0_z_0_yz_0, g_0_z_0_zz_0, wp_z  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_z_0_xx_0[i] = g_0_0_0_xx_0[i] * pb_z + g_0_0_0_xx_1[i] * wp_z[i];

        g_0_z_0_xy_0[i] = g_0_0_0_xy_0[i] * pb_z + g_0_0_0_xy_1[i] * wp_z[i];

        g_0_z_0_xz_0[i] = g_0_0_0_x_1[i] * fi_abcd_0 + g_0_0_0_xz_0[i] * pb_z + g_0_0_0_xz_1[i] * wp_z[i];

        g_0_z_0_yy_0[i] = g_0_0_0_yy_0[i] * pb_z + g_0_0_0_yy_1[i] * wp_z[i];

        g_0_z_0_yz_0[i] = g_0_0_0_y_1[i] * fi_abcd_0 + g_0_0_0_yz_0[i] * pb_z + g_0_0_0_yz_1[i] * wp_z[i];

        g_0_z_0_zz_0[i] = 2.0 * g_0_0_0_z_1[i] * fi_abcd_0 + g_0_0_0_zz_0[i] * pb_z + g_0_0_0_zz_1[i] * wp_z[i];
    }
}

} // erirec namespace


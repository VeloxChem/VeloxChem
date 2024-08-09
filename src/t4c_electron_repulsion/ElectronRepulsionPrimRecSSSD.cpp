#include "ElectronRepulsionPrimRecSSSD.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_sssd(CSimdArray<double>& prim_buffer_0_sssd,
                                  const CSimdArray<double>& prim_buffer_0_ssss,
                                  const CSimdArray<double>& prim_buffer_1_ssss,
                                  const CSimdArray<double>& prim_buffer_0_sssp,
                                  const CSimdArray<double>& prim_buffer_1_sssp,
                                  const double* qd_x,
                                  const double* qd_y,
                                  const double* qd_z,
                                  const double* wq_x,
                                  const double* wq_y,
                                  const double* wq_z,
                                  const double a_exp,
                                  const double b_exp,
                                  const double* c_exps,
                                  const double* d_exps) -> void
{
    const auto ndims = prim_buffer_0_sssd.number_of_columns();

    /// Set up components of auxilary buffer : prim_buffer_0_ssss

    auto g_0_0_0_0_0 = prim_buffer_0_ssss[0];

    /// Set up components of auxilary buffer : prim_buffer_1_ssss

    auto g_0_0_0_0_1 = prim_buffer_1_ssss[0];

    /// Set up components of auxilary buffer : prim_buffer_0_sssp

    auto g_0_0_0_x_0 = prim_buffer_0_sssp[0];

    auto g_0_0_0_y_0 = prim_buffer_0_sssp[1];

    auto g_0_0_0_z_0 = prim_buffer_0_sssp[2];

    /// Set up components of auxilary buffer : prim_buffer_1_sssp

    auto g_0_0_0_x_1 = prim_buffer_1_sssp[0];

    auto g_0_0_0_y_1 = prim_buffer_1_sssp[1];

    auto g_0_0_0_z_1 = prim_buffer_1_sssp[2];

    /// Set up components of targeted buffer : prim_buffer_0_sssd

    auto g_0_0_0_xx_0 = prim_buffer_0_sssd[0];

    auto g_0_0_0_xy_0 = prim_buffer_0_sssd[1];

    auto g_0_0_0_xz_0 = prim_buffer_0_sssd[2];

    auto g_0_0_0_yy_0 = prim_buffer_0_sssd[3];

    auto g_0_0_0_yz_0 = prim_buffer_0_sssd[4];

    auto g_0_0_0_zz_0 = prim_buffer_0_sssd[5];

    #pragma omp simd aligned(g_0_0_0_0_0, g_0_0_0_0_1, g_0_0_0_x_0, g_0_0_0_x_1, g_0_0_0_xx_0, g_0_0_0_xy_0, g_0_0_0_xz_0, g_0_0_0_y_0, g_0_0_0_y_1, g_0_0_0_yy_0, g_0_0_0_yz_0, g_0_0_0_z_0, g_0_0_0_z_1, g_0_0_0_zz_0, qd_x, qd_y, qd_z, wq_x, wq_y, wq_z  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_cd_0 = 0.5 / (c_exps[i] + d_exps[i]);

        const double fti_cd_0 =  fi_cd_0 * (a_exp + b_exp) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_0_0_xx_0[i] = g_0_0_0_0_0[i] * fi_cd_0 - g_0_0_0_0_1[i] * fti_cd_0 + g_0_0_0_x_0[i] * qd_x[i] + g_0_0_0_x_1[i] * wq_x[i];

        g_0_0_0_xy_0[i] = g_0_0_0_y_0[i] * qd_x[i] + g_0_0_0_y_1[i] * wq_x[i];

        g_0_0_0_xz_0[i] = g_0_0_0_z_0[i] * qd_x[i] + g_0_0_0_z_1[i] * wq_x[i];

        g_0_0_0_yy_0[i] = g_0_0_0_0_0[i] * fi_cd_0 - g_0_0_0_0_1[i] * fti_cd_0 + g_0_0_0_y_0[i] * qd_y[i] + g_0_0_0_y_1[i] * wq_y[i];

        g_0_0_0_yz_0[i] = g_0_0_0_z_0[i] * qd_y[i] + g_0_0_0_z_1[i] * wq_y[i];

        g_0_0_0_zz_0[i] = g_0_0_0_0_0[i] * fi_cd_0 - g_0_0_0_0_1[i] * fti_cd_0 + g_0_0_0_z_0[i] * qd_z[i] + g_0_0_0_z_1[i] * wq_z[i];
    }
}

} // erirec namespace


#include "ElectronRepulsionPrimRecSSSD.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_sssd(CSimdArray<double>& pbuffer,
                                  const size_t idx_eri_0_sssd,
                                  size_t idx_eri_0_ssss,
                                  size_t idx_eri_1_ssss,
                                  size_t idx_eri_0_sssp,
                                  size_t idx_eri_1_sssp,
                                  CSimdArray<double>& factors,
                                  const size_t idx_qd,
                                  const size_t idx_wq,
                                  const double a_exp,
                                  const double b_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto c_exps = factors.data(0);

    auto d_exps = factors.data(1);

    // Set up R(QD) distances

    auto qd_x = factors.data(idx_qd);

    auto qd_y = factors.data(idx_qd + 1);

    auto qd_z = factors.data(idx_qd + 2);

    // Set up R(WQ) distances

    auto wq_x = factors.data(idx_wq);

    auto wq_y = factors.data(idx_wq + 1);

    auto wq_z = factors.data(idx_wq + 2);

    /// Set up components of auxilary buffer : SSSS

    auto g_0_0_0_0_0 = pbuffer.data(idx_eri_0_ssss);

    /// Set up components of auxilary buffer : SSSS

    auto g_0_0_0_0_1 = pbuffer.data(idx_eri_1_ssss);

    /// Set up components of auxilary buffer : SSSP

    auto g_0_0_0_x_0 = pbuffer.data(idx_eri_0_sssp);

    auto g_0_0_0_y_0 = pbuffer.data(idx_eri_0_sssp + 1);

    auto g_0_0_0_z_0 = pbuffer.data(idx_eri_0_sssp + 2);

    /// Set up components of auxilary buffer : SSSP

    auto g_0_0_0_x_1 = pbuffer.data(idx_eri_1_sssp);

    auto g_0_0_0_y_1 = pbuffer.data(idx_eri_1_sssp + 1);

    auto g_0_0_0_z_1 = pbuffer.data(idx_eri_1_sssp + 2);

    /// Set up components of targeted buffer : SSSD

    auto g_0_0_0_xx_0 = pbuffer.data(idx_eri_0_sssd);

    auto g_0_0_0_xy_0 = pbuffer.data(idx_eri_0_sssd + 1);

    auto g_0_0_0_xz_0 = pbuffer.data(idx_eri_0_sssd + 2);

    auto g_0_0_0_yy_0 = pbuffer.data(idx_eri_0_sssd + 3);

    auto g_0_0_0_yz_0 = pbuffer.data(idx_eri_0_sssd + 4);

    auto g_0_0_0_zz_0 = pbuffer.data(idx_eri_0_sssd + 5);

    #pragma omp simd aligned(g_0_0_0_0_0, g_0_0_0_0_1, g_0_0_0_x_0, g_0_0_0_x_1, g_0_0_0_xx_0, g_0_0_0_xy_0, g_0_0_0_xz_0, g_0_0_0_y_0, g_0_0_0_y_1, g_0_0_0_yy_0, g_0_0_0_yz_0, g_0_0_0_z_0, g_0_0_0_z_1, g_0_0_0_zz_0, qd_x, qd_y, qd_z, wq_x, wq_y, wq_z  : 64)
    for (size_t i = 0; i < nelems; i++)
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


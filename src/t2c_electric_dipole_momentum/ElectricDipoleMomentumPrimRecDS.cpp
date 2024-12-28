#include "ElectricDipoleMomentumPrimRecDS.hpp"

namespace diprec {  // diprec namespace

auto
comp_prim_electric_dipole_momentum_ds(CSimdArray<double>&       pbuffer,
                                      const size_t              idx_dip_ds,
                                      const size_t              idx_dip_ss,
                                      const size_t              idx_ovl_ps,
                                      const size_t              idx_dip_ps,
                                      const CSimdArray<double>& factors,
                                      const size_t              idx_rpa,
                                      const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : SS

    auto tr_x_0_0 = pbuffer.data(idx_dip_ss);

    auto tr_y_0_0 = pbuffer.data(idx_dip_ss + 1);

    auto tr_z_0_0 = pbuffer.data(idx_dip_ss + 2);

    // Set up components of auxiliary buffer : PS

    auto ts_x_0 = pbuffer.data(idx_ovl_ps);

    auto ts_y_0 = pbuffer.data(idx_ovl_ps + 1);

    auto ts_z_0 = pbuffer.data(idx_ovl_ps + 2);

    // Set up components of auxiliary buffer : PS

    auto tr_x_x_0 = pbuffer.data(idx_dip_ps);

    auto tr_x_y_0 = pbuffer.data(idx_dip_ps + 1);

    auto tr_x_z_0 = pbuffer.data(idx_dip_ps + 2);

    auto tr_y_x_0 = pbuffer.data(idx_dip_ps + 3);

    auto tr_y_y_0 = pbuffer.data(idx_dip_ps + 4);

    auto tr_y_z_0 = pbuffer.data(idx_dip_ps + 5);

    auto tr_z_x_0 = pbuffer.data(idx_dip_ps + 6);

    auto tr_z_y_0 = pbuffer.data(idx_dip_ps + 7);

    auto tr_z_z_0 = pbuffer.data(idx_dip_ps + 8);

    // Set up components of targeted buffer : DS

    auto tr_x_xx_0 = pbuffer.data(idx_dip_ds);

    auto tr_x_xy_0 = pbuffer.data(idx_dip_ds + 1);

    auto tr_x_xz_0 = pbuffer.data(idx_dip_ds + 2);

    auto tr_x_yy_0 = pbuffer.data(idx_dip_ds + 3);

    auto tr_x_yz_0 = pbuffer.data(idx_dip_ds + 4);

    auto tr_x_zz_0 = pbuffer.data(idx_dip_ds + 5);

    auto tr_y_xx_0 = pbuffer.data(idx_dip_ds + 6);

    auto tr_y_xy_0 = pbuffer.data(idx_dip_ds + 7);

    auto tr_y_xz_0 = pbuffer.data(idx_dip_ds + 8);

    auto tr_y_yy_0 = pbuffer.data(idx_dip_ds + 9);

    auto tr_y_yz_0 = pbuffer.data(idx_dip_ds + 10);

    auto tr_y_zz_0 = pbuffer.data(idx_dip_ds + 11);

    auto tr_z_xx_0 = pbuffer.data(idx_dip_ds + 12);

    auto tr_z_xy_0 = pbuffer.data(idx_dip_ds + 13);

    auto tr_z_xz_0 = pbuffer.data(idx_dip_ds + 14);

    auto tr_z_yy_0 = pbuffer.data(idx_dip_ds + 15);

    auto tr_z_yz_0 = pbuffer.data(idx_dip_ds + 16);

    auto tr_z_zz_0 = pbuffer.data(idx_dip_ds + 17);

#pragma omp simd aligned(pa_x,          \
                             pa_y,      \
                             pa_z,      \
                             tr_x_0_0,  \
                             tr_x_x_0,  \
                             tr_x_xx_0, \
                             tr_x_xy_0, \
                             tr_x_xz_0, \
                             tr_x_y_0,  \
                             tr_x_yy_0, \
                             tr_x_yz_0, \
                             tr_x_z_0,  \
                             tr_x_zz_0, \
                             tr_y_0_0,  \
                             tr_y_x_0,  \
                             tr_y_xx_0, \
                             tr_y_xy_0, \
                             tr_y_xz_0, \
                             tr_y_y_0,  \
                             tr_y_yy_0, \
                             tr_y_yz_0, \
                             tr_y_z_0,  \
                             tr_y_zz_0, \
                             tr_z_0_0,  \
                             tr_z_x_0,  \
                             tr_z_xx_0, \
                             tr_z_xy_0, \
                             tr_z_xz_0, \
                             tr_z_y_0,  \
                             tr_z_yy_0, \
                             tr_z_yz_0, \
                             tr_z_z_0,  \
                             tr_z_zz_0, \
                             ts_x_0,    \
                             ts_y_0,    \
                             ts_z_0,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xx_0[i] = tr_x_0_0[i] * fe_0 + ts_x_0[i] * fe_0 + tr_x_x_0[i] * pa_x[i];

        tr_x_xy_0[i] = tr_x_x_0[i] * pa_y[i];

        tr_x_xz_0[i] = tr_x_x_0[i] * pa_z[i];

        tr_x_yy_0[i] = tr_x_0_0[i] * fe_0 + tr_x_y_0[i] * pa_y[i];

        tr_x_yz_0[i] = tr_x_z_0[i] * pa_y[i];

        tr_x_zz_0[i] = tr_x_0_0[i] * fe_0 + tr_x_z_0[i] * pa_z[i];

        tr_y_xx_0[i] = tr_y_0_0[i] * fe_0 + tr_y_x_0[i] * pa_x[i];

        tr_y_xy_0[i] = tr_y_y_0[i] * pa_x[i];

        tr_y_xz_0[i] = tr_y_z_0[i] * pa_x[i];

        tr_y_yy_0[i] = tr_y_0_0[i] * fe_0 + ts_y_0[i] * fe_0 + tr_y_y_0[i] * pa_y[i];

        tr_y_yz_0[i] = tr_y_y_0[i] * pa_z[i];

        tr_y_zz_0[i] = tr_y_0_0[i] * fe_0 + tr_y_z_0[i] * pa_z[i];

        tr_z_xx_0[i] = tr_z_0_0[i] * fe_0 + tr_z_x_0[i] * pa_x[i];

        tr_z_xy_0[i] = tr_z_y_0[i] * pa_x[i];

        tr_z_xz_0[i] = tr_z_z_0[i] * pa_x[i];

        tr_z_yy_0[i] = tr_z_0_0[i] * fe_0 + tr_z_y_0[i] * pa_y[i];

        tr_z_yz_0[i] = tr_z_z_0[i] * pa_y[i];

        tr_z_zz_0[i] = tr_z_0_0[i] * fe_0 + ts_z_0[i] * fe_0 + tr_z_z_0[i] * pa_z[i];
    }
}

}  // namespace diprec
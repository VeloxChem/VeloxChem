#include "KineticEnergyPrimRecDS.hpp"

namespace kinrec {  // kinrec namespace

auto
comp_prim_kinetic_energy_ds(CSimdArray<double>&       pbuffer,
                            const size_t              idx_kin_ds,
                            const size_t              idx_ovl_ss,
                            const size_t              idx_kin_ss,
                            const size_t              idx_kin_ps,
                            const size_t              idx_ovl_ds,
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

    auto ts_0_0 = pbuffer.data(idx_ovl_ss);

    // Set up components of auxiliary buffer : SS

    auto tk_0_0 = pbuffer.data(idx_kin_ss);

    // Set up components of auxiliary buffer : PS

    auto tk_x_0 = pbuffer.data(idx_kin_ps);

    auto tk_y_0 = pbuffer.data(idx_kin_ps + 1);

    auto tk_z_0 = pbuffer.data(idx_kin_ps + 2);

    // Set up components of auxiliary buffer : DS

    auto ts_xx_0 = pbuffer.data(idx_ovl_ds);

    auto ts_xy_0 = pbuffer.data(idx_ovl_ds + 1);

    auto ts_xz_0 = pbuffer.data(idx_ovl_ds + 2);

    auto ts_yy_0 = pbuffer.data(idx_ovl_ds + 3);

    auto ts_yz_0 = pbuffer.data(idx_ovl_ds + 4);

    auto ts_zz_0 = pbuffer.data(idx_ovl_ds + 5);

    // Set up components of targeted buffer : DS

    auto tk_xx_0 = pbuffer.data(idx_kin_ds);

    auto tk_xy_0 = pbuffer.data(idx_kin_ds + 1);

    auto tk_xz_0 = pbuffer.data(idx_kin_ds + 2);

    auto tk_yy_0 = pbuffer.data(idx_kin_ds + 3);

    auto tk_yz_0 = pbuffer.data(idx_kin_ds + 4);

    auto tk_zz_0 = pbuffer.data(idx_kin_ds + 5);

#pragma omp simd aligned(pa_x,        \
                             pa_y,    \
                             pa_z,    \
                             tk_0_0,  \
                             tk_x_0,  \
                             tk_xx_0, \
                             tk_xy_0, \
                             tk_xz_0, \
                             tk_y_0,  \
                             tk_yy_0, \
                             tk_yz_0, \
                             tk_z_0,  \
                             tk_zz_0, \
                             ts_0_0,  \
                             ts_xx_0, \
                             ts_xy_0, \
                             ts_xz_0, \
                             ts_yy_0, \
                             ts_yz_0, \
                             ts_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xx_0[i] = -2.0 * ts_0_0[i] * fbe_0 * fz_0 + tk_0_0[i] * fe_0 + tk_x_0[i] * pa_x[i] + 2.0 * ts_xx_0[i] * fz_0;

        tk_xy_0[i] = tk_y_0[i] * pa_x[i] + 2.0 * ts_xy_0[i] * fz_0;

        tk_xz_0[i] = tk_z_0[i] * pa_x[i] + 2.0 * ts_xz_0[i] * fz_0;

        tk_yy_0[i] = -2.0 * ts_0_0[i] * fbe_0 * fz_0 + tk_0_0[i] * fe_0 + tk_y_0[i] * pa_y[i] + 2.0 * ts_yy_0[i] * fz_0;

        tk_yz_0[i] = tk_z_0[i] * pa_y[i] + 2.0 * ts_yz_0[i] * fz_0;

        tk_zz_0[i] = -2.0 * ts_0_0[i] * fbe_0 * fz_0 + tk_0_0[i] * fe_0 + tk_z_0[i] * pa_z[i] + 2.0 * ts_zz_0[i] * fz_0;
    }
}

}  // namespace kinrec

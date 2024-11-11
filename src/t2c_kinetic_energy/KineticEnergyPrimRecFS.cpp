#include "KineticEnergyPrimRecFS.hpp"

namespace kinrec {  // kinrec namespace

auto
comp_prim_kinetic_energy_fs(CSimdArray<double>&       pbuffer,
                            const size_t              idx_kin_fs,
                            const size_t              idx_ovl_ps,
                            const size_t              idx_kin_ps,
                            const size_t              idx_kin_ds,
                            const size_t              idx_ovl_fs,
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

    // Set up components of auxiliary buffer : PS

    auto ts_x_0 = pbuffer.data(idx_ovl_ps);

    auto ts_y_0 = pbuffer.data(idx_ovl_ps + 1);

    auto ts_z_0 = pbuffer.data(idx_ovl_ps + 2);

    // Set up components of auxiliary buffer : PS

    auto tk_x_0 = pbuffer.data(idx_kin_ps);

    auto tk_y_0 = pbuffer.data(idx_kin_ps + 1);

    auto tk_z_0 = pbuffer.data(idx_kin_ps + 2);

    // Set up components of auxiliary buffer : DS

    auto tk_xx_0 = pbuffer.data(idx_kin_ds);

    auto tk_yy_0 = pbuffer.data(idx_kin_ds + 3);

    auto tk_yz_0 = pbuffer.data(idx_kin_ds + 4);

    auto tk_zz_0 = pbuffer.data(idx_kin_ds + 5);

    // Set up components of auxiliary buffer : FS

    auto ts_xxx_0 = pbuffer.data(idx_ovl_fs);

    auto ts_xxy_0 = pbuffer.data(idx_ovl_fs + 1);

    auto ts_xxz_0 = pbuffer.data(idx_ovl_fs + 2);

    auto ts_xyy_0 = pbuffer.data(idx_ovl_fs + 3);

    auto ts_xyz_0 = pbuffer.data(idx_ovl_fs + 4);

    auto ts_xzz_0 = pbuffer.data(idx_ovl_fs + 5);

    auto ts_yyy_0 = pbuffer.data(idx_ovl_fs + 6);

    auto ts_yyz_0 = pbuffer.data(idx_ovl_fs + 7);

    auto ts_yzz_0 = pbuffer.data(idx_ovl_fs + 8);

    auto ts_zzz_0 = pbuffer.data(idx_ovl_fs + 9);

    // Set up components of targeted buffer : FS

    auto tk_xxx_0 = pbuffer.data(idx_kin_fs);

    auto tk_xxy_0 = pbuffer.data(idx_kin_fs + 1);

    auto tk_xxz_0 = pbuffer.data(idx_kin_fs + 2);

    auto tk_xyy_0 = pbuffer.data(idx_kin_fs + 3);

    auto tk_xyz_0 = pbuffer.data(idx_kin_fs + 4);

    auto tk_xzz_0 = pbuffer.data(idx_kin_fs + 5);

    auto tk_yyy_0 = pbuffer.data(idx_kin_fs + 6);

    auto tk_yyz_0 = pbuffer.data(idx_kin_fs + 7);

    auto tk_yzz_0 = pbuffer.data(idx_kin_fs + 8);

    auto tk_zzz_0 = pbuffer.data(idx_kin_fs + 9);

#pragma omp simd aligned(pa_x,         \
                             pa_y,     \
                             pa_z,     \
                             tk_x_0,   \
                             tk_xx_0,  \
                             tk_xxx_0, \
                             tk_xxy_0, \
                             tk_xxz_0, \
                             tk_xyy_0, \
                             tk_xyz_0, \
                             tk_xzz_0, \
                             tk_y_0,   \
                             tk_yy_0,  \
                             tk_yyy_0, \
                             tk_yyz_0, \
                             tk_yz_0,  \
                             tk_yzz_0, \
                             tk_z_0,   \
                             tk_zz_0,  \
                             tk_zzz_0, \
                             ts_x_0,   \
                             ts_xxx_0, \
                             ts_xxy_0, \
                             ts_xxz_0, \
                             ts_xyy_0, \
                             ts_xyz_0, \
                             ts_xzz_0, \
                             ts_y_0,   \
                             ts_yyy_0, \
                             ts_yyz_0, \
                             ts_yzz_0, \
                             ts_z_0,   \
                             ts_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxx_0[i] = -4.0 * ts_x_0[i] * fbe_0 * fz_0 + 2.0 * tk_x_0[i] * fe_0 + tk_xx_0[i] * pa_x[i] + 2.0 * ts_xxx_0[i] * fz_0;

        tk_xxy_0[i] = tk_xx_0[i] * pa_y[i] + 2.0 * ts_xxy_0[i] * fz_0;

        tk_xxz_0[i] = tk_xx_0[i] * pa_z[i] + 2.0 * ts_xxz_0[i] * fz_0;

        tk_xyy_0[i] = tk_yy_0[i] * pa_x[i] + 2.0 * ts_xyy_0[i] * fz_0;

        tk_xyz_0[i] = tk_yz_0[i] * pa_x[i] + 2.0 * ts_xyz_0[i] * fz_0;

        tk_xzz_0[i] = tk_zz_0[i] * pa_x[i] + 2.0 * ts_xzz_0[i] * fz_0;

        tk_yyy_0[i] = -4.0 * ts_y_0[i] * fbe_0 * fz_0 + 2.0 * tk_y_0[i] * fe_0 + tk_yy_0[i] * pa_y[i] + 2.0 * ts_yyy_0[i] * fz_0;

        tk_yyz_0[i] = tk_yy_0[i] * pa_z[i] + 2.0 * ts_yyz_0[i] * fz_0;

        tk_yzz_0[i] = tk_zz_0[i] * pa_y[i] + 2.0 * ts_yzz_0[i] * fz_0;

        tk_zzz_0[i] = -4.0 * ts_z_0[i] * fbe_0 * fz_0 + 2.0 * tk_z_0[i] * fe_0 + tk_zz_0[i] * pa_z[i] + 2.0 * ts_zzz_0[i] * fz_0;
    }
}

}  // namespace kinrec

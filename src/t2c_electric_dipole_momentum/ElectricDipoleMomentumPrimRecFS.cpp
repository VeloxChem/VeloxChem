#include "ElectricDipoleMomentumPrimRecFS.hpp"

namespace diprec { // diprec namespace

auto
comp_prim_electric_dipole_momentum_fs(CSimdArray<double>& pbuffer, 
                                      const size_t idx_dip_fs,
                                      const size_t idx_dip_ps,
                                      const size_t idx_ovl_ds,
                                      const size_t idx_dip_ds,
                                      const CSimdArray<double>& factors,
                                      const size_t idx_rpa,
                                      const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

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

    // Set up components of auxiliary buffer : DS

    auto ts_xx_0 = pbuffer.data(idx_ovl_ds);

    auto ts_yy_0 = pbuffer.data(idx_ovl_ds + 3);

    auto ts_zz_0 = pbuffer.data(idx_ovl_ds + 5);

    // Set up components of auxiliary buffer : DS

    auto tr_x_xx_0 = pbuffer.data(idx_dip_ds);

    auto tr_x_xz_0 = pbuffer.data(idx_dip_ds + 2);

    auto tr_x_yy_0 = pbuffer.data(idx_dip_ds + 3);

    auto tr_x_zz_0 = pbuffer.data(idx_dip_ds + 5);

    auto tr_y_xx_0 = pbuffer.data(idx_dip_ds + 6);

    auto tr_y_xy_0 = pbuffer.data(idx_dip_ds + 7);

    auto tr_y_yy_0 = pbuffer.data(idx_dip_ds + 9);

    auto tr_y_yz_0 = pbuffer.data(idx_dip_ds + 10);

    auto tr_y_zz_0 = pbuffer.data(idx_dip_ds + 11);

    auto tr_z_xx_0 = pbuffer.data(idx_dip_ds + 12);

    auto tr_z_xz_0 = pbuffer.data(idx_dip_ds + 14);

    auto tr_z_yy_0 = pbuffer.data(idx_dip_ds + 15);

    auto tr_z_yz_0 = pbuffer.data(idx_dip_ds + 16);

    auto tr_z_zz_0 = pbuffer.data(idx_dip_ds + 17);

    // Set up components of targeted buffer : FS

    auto tr_x_xxx_0 = pbuffer.data(idx_dip_fs);

    auto tr_x_xxy_0 = pbuffer.data(idx_dip_fs + 1);

    auto tr_x_xxz_0 = pbuffer.data(idx_dip_fs + 2);

    auto tr_x_xyy_0 = pbuffer.data(idx_dip_fs + 3);

    auto tr_x_xyz_0 = pbuffer.data(idx_dip_fs + 4);

    auto tr_x_xzz_0 = pbuffer.data(idx_dip_fs + 5);

    auto tr_x_yyy_0 = pbuffer.data(idx_dip_fs + 6);

    auto tr_x_yyz_0 = pbuffer.data(idx_dip_fs + 7);

    auto tr_x_yzz_0 = pbuffer.data(idx_dip_fs + 8);

    auto tr_x_zzz_0 = pbuffer.data(idx_dip_fs + 9);

    auto tr_y_xxx_0 = pbuffer.data(idx_dip_fs + 10);

    auto tr_y_xxy_0 = pbuffer.data(idx_dip_fs + 11);

    auto tr_y_xxz_0 = pbuffer.data(idx_dip_fs + 12);

    auto tr_y_xyy_0 = pbuffer.data(idx_dip_fs + 13);

    auto tr_y_xyz_0 = pbuffer.data(idx_dip_fs + 14);

    auto tr_y_xzz_0 = pbuffer.data(idx_dip_fs + 15);

    auto tr_y_yyy_0 = pbuffer.data(idx_dip_fs + 16);

    auto tr_y_yyz_0 = pbuffer.data(idx_dip_fs + 17);

    auto tr_y_yzz_0 = pbuffer.data(idx_dip_fs + 18);

    auto tr_y_zzz_0 = pbuffer.data(idx_dip_fs + 19);

    auto tr_z_xxx_0 = pbuffer.data(idx_dip_fs + 20);

    auto tr_z_xxy_0 = pbuffer.data(idx_dip_fs + 21);

    auto tr_z_xxz_0 = pbuffer.data(idx_dip_fs + 22);

    auto tr_z_xyy_0 = pbuffer.data(idx_dip_fs + 23);

    auto tr_z_xyz_0 = pbuffer.data(idx_dip_fs + 24);

    auto tr_z_xzz_0 = pbuffer.data(idx_dip_fs + 25);

    auto tr_z_yyy_0 = pbuffer.data(idx_dip_fs + 26);

    auto tr_z_yyz_0 = pbuffer.data(idx_dip_fs + 27);

    auto tr_z_yzz_0 = pbuffer.data(idx_dip_fs + 28);

    auto tr_z_zzz_0 = pbuffer.data(idx_dip_fs + 29);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_x_x_0, tr_x_xx_0, tr_x_xxx_0, tr_x_xxy_0, tr_x_xxz_0, tr_x_xyy_0, tr_x_xyz_0, tr_x_xz_0, tr_x_xzz_0, tr_x_y_0, tr_x_yy_0, tr_x_yyy_0, tr_x_yyz_0, tr_x_yzz_0, tr_x_z_0, tr_x_zz_0, tr_x_zzz_0, tr_y_x_0, tr_y_xx_0, tr_y_xxx_0, tr_y_xxy_0, tr_y_xxz_0, tr_y_xy_0, tr_y_xyy_0, tr_y_xyz_0, tr_y_xzz_0, tr_y_y_0, tr_y_yy_0, tr_y_yyy_0, tr_y_yyz_0, tr_y_yz_0, tr_y_yzz_0, tr_y_z_0, tr_y_zz_0, tr_y_zzz_0, tr_z_x_0, tr_z_xx_0, tr_z_xxx_0, tr_z_xxy_0, tr_z_xxz_0, tr_z_xyy_0, tr_z_xyz_0, tr_z_xz_0, tr_z_xzz_0, tr_z_y_0, tr_z_yy_0, tr_z_yyy_0, tr_z_yyz_0, tr_z_yz_0, tr_z_yzz_0, tr_z_z_0, tr_z_zz_0, tr_z_zzz_0, ts_xx_0, ts_yy_0, ts_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxx_0[i] = 2.0 * tr_x_x_0[i] * fe_0 + ts_xx_0[i] * fe_0 + tr_x_xx_0[i] * pa_x[i];

        tr_x_xxy_0[i] = tr_x_xx_0[i] * pa_y[i];

        tr_x_xxz_0[i] = tr_x_xx_0[i] * pa_z[i];

        tr_x_xyy_0[i] = ts_yy_0[i] * fe_0 + tr_x_yy_0[i] * pa_x[i];

        tr_x_xyz_0[i] = tr_x_xz_0[i] * pa_y[i];

        tr_x_xzz_0[i] = ts_zz_0[i] * fe_0 + tr_x_zz_0[i] * pa_x[i];

        tr_x_yyy_0[i] = 2.0 * tr_x_y_0[i] * fe_0 + tr_x_yy_0[i] * pa_y[i];

        tr_x_yyz_0[i] = tr_x_yy_0[i] * pa_z[i];

        tr_x_yzz_0[i] = tr_x_zz_0[i] * pa_y[i];

        tr_x_zzz_0[i] = 2.0 * tr_x_z_0[i] * fe_0 + tr_x_zz_0[i] * pa_z[i];

        tr_y_xxx_0[i] = 2.0 * tr_y_x_0[i] * fe_0 + tr_y_xx_0[i] * pa_x[i];

        tr_y_xxy_0[i] = tr_y_y_0[i] * fe_0 + tr_y_xy_0[i] * pa_x[i];

        tr_y_xxz_0[i] = tr_y_xx_0[i] * pa_z[i];

        tr_y_xyy_0[i] = tr_y_yy_0[i] * pa_x[i];

        tr_y_xyz_0[i] = tr_y_yz_0[i] * pa_x[i];

        tr_y_xzz_0[i] = tr_y_zz_0[i] * pa_x[i];

        tr_y_yyy_0[i] = 2.0 * tr_y_y_0[i] * fe_0 + ts_yy_0[i] * fe_0 + tr_y_yy_0[i] * pa_y[i];

        tr_y_yyz_0[i] = tr_y_yy_0[i] * pa_z[i];

        tr_y_yzz_0[i] = ts_zz_0[i] * fe_0 + tr_y_zz_0[i] * pa_y[i];

        tr_y_zzz_0[i] = 2.0 * tr_y_z_0[i] * fe_0 + tr_y_zz_0[i] * pa_z[i];

        tr_z_xxx_0[i] = 2.0 * tr_z_x_0[i] * fe_0 + tr_z_xx_0[i] * pa_x[i];

        tr_z_xxy_0[i] = tr_z_xx_0[i] * pa_y[i];

        tr_z_xxz_0[i] = tr_z_z_0[i] * fe_0 + tr_z_xz_0[i] * pa_x[i];

        tr_z_xyy_0[i] = tr_z_yy_0[i] * pa_x[i];

        tr_z_xyz_0[i] = tr_z_yz_0[i] * pa_x[i];

        tr_z_xzz_0[i] = tr_z_zz_0[i] * pa_x[i];

        tr_z_yyy_0[i] = 2.0 * tr_z_y_0[i] * fe_0 + tr_z_yy_0[i] * pa_y[i];

        tr_z_yyz_0[i] = tr_z_z_0[i] * fe_0 + tr_z_yz_0[i] * pa_y[i];

        tr_z_yzz_0[i] = tr_z_zz_0[i] * pa_y[i];

        tr_z_zzz_0[i] = 2.0 * tr_z_z_0[i] * fe_0 + ts_zz_0[i] * fe_0 + tr_z_zz_0[i] * pa_z[i];
    }
}

} // diprec namespace


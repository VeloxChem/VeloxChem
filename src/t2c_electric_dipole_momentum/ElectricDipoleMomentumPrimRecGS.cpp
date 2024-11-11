#include "ElectricDipoleMomentumPrimRecGS.hpp"

namespace diprec {  // diprec namespace

auto
comp_prim_electric_dipole_momentum_gs(CSimdArray<double>&       pbuffer,
                                      const size_t              idx_dip_gs,
                                      const size_t              idx_dip_ds,
                                      const size_t              idx_ovl_fs,
                                      const size_t              idx_dip_fs,
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

    // Set up components of auxiliary buffer : DS

    auto tr_x_xx_0 = pbuffer.data(idx_dip_ds);

    auto tr_x_yy_0 = pbuffer.data(idx_dip_ds + 3);

    auto tr_x_zz_0 = pbuffer.data(idx_dip_ds + 5);

    auto tr_y_xx_0 = pbuffer.data(idx_dip_ds + 6);

    auto tr_y_xy_0 = pbuffer.data(idx_dip_ds + 7);

    auto tr_y_yy_0 = pbuffer.data(idx_dip_ds + 9);

    auto tr_y_zz_0 = pbuffer.data(idx_dip_ds + 11);

    auto tr_z_xx_0 = pbuffer.data(idx_dip_ds + 12);

    auto tr_z_xz_0 = pbuffer.data(idx_dip_ds + 14);

    auto tr_z_yy_0 = pbuffer.data(idx_dip_ds + 15);

    auto tr_z_yz_0 = pbuffer.data(idx_dip_ds + 16);

    auto tr_z_zz_0 = pbuffer.data(idx_dip_ds + 17);

    // Set up components of auxiliary buffer : FS

    auto ts_xxx_0 = pbuffer.data(idx_ovl_fs);

    auto ts_yyy_0 = pbuffer.data(idx_ovl_fs + 6);

    auto ts_zzz_0 = pbuffer.data(idx_ovl_fs + 9);

    // Set up components of auxiliary buffer : FS

    auto tr_x_xxx_0 = pbuffer.data(idx_dip_fs);

    auto tr_x_xxy_0 = pbuffer.data(idx_dip_fs + 1);

    auto tr_x_xxz_0 = pbuffer.data(idx_dip_fs + 2);

    auto tr_x_xyy_0 = pbuffer.data(idx_dip_fs + 3);

    auto tr_x_xzz_0 = pbuffer.data(idx_dip_fs + 5);

    auto tr_x_yyy_0 = pbuffer.data(idx_dip_fs + 6);

    auto tr_x_yzz_0 = pbuffer.data(idx_dip_fs + 8);

    auto tr_x_zzz_0 = pbuffer.data(idx_dip_fs + 9);

    auto tr_y_xxx_0 = pbuffer.data(idx_dip_fs + 10);

    auto tr_y_xxy_0 = pbuffer.data(idx_dip_fs + 11);

    auto tr_y_xyy_0 = pbuffer.data(idx_dip_fs + 13);

    auto tr_y_xzz_0 = pbuffer.data(idx_dip_fs + 15);

    auto tr_y_yyy_0 = pbuffer.data(idx_dip_fs + 16);

    auto tr_y_yyz_0 = pbuffer.data(idx_dip_fs + 17);

    auto tr_y_yzz_0 = pbuffer.data(idx_dip_fs + 18);

    auto tr_y_zzz_0 = pbuffer.data(idx_dip_fs + 19);

    auto tr_z_xxx_0 = pbuffer.data(idx_dip_fs + 20);

    auto tr_z_xxz_0 = pbuffer.data(idx_dip_fs + 22);

    auto tr_z_xyy_0 = pbuffer.data(idx_dip_fs + 23);

    auto tr_z_xzz_0 = pbuffer.data(idx_dip_fs + 25);

    auto tr_z_yyy_0 = pbuffer.data(idx_dip_fs + 26);

    auto tr_z_yyz_0 = pbuffer.data(idx_dip_fs + 27);

    auto tr_z_yzz_0 = pbuffer.data(idx_dip_fs + 28);

    auto tr_z_zzz_0 = pbuffer.data(idx_dip_fs + 29);

    // Set up components of targeted buffer : GS

    auto tr_x_xxxx_0 = pbuffer.data(idx_dip_gs);

    auto tr_x_xxxy_0 = pbuffer.data(idx_dip_gs + 1);

    auto tr_x_xxxz_0 = pbuffer.data(idx_dip_gs + 2);

    auto tr_x_xxyy_0 = pbuffer.data(idx_dip_gs + 3);

    auto tr_x_xxyz_0 = pbuffer.data(idx_dip_gs + 4);

    auto tr_x_xxzz_0 = pbuffer.data(idx_dip_gs + 5);

    auto tr_x_xyyy_0 = pbuffer.data(idx_dip_gs + 6);

    auto tr_x_xyyz_0 = pbuffer.data(idx_dip_gs + 7);

    auto tr_x_xyzz_0 = pbuffer.data(idx_dip_gs + 8);

    auto tr_x_xzzz_0 = pbuffer.data(idx_dip_gs + 9);

    auto tr_x_yyyy_0 = pbuffer.data(idx_dip_gs + 10);

    auto tr_x_yyyz_0 = pbuffer.data(idx_dip_gs + 11);

    auto tr_x_yyzz_0 = pbuffer.data(idx_dip_gs + 12);

    auto tr_x_yzzz_0 = pbuffer.data(idx_dip_gs + 13);

    auto tr_x_zzzz_0 = pbuffer.data(idx_dip_gs + 14);

    auto tr_y_xxxx_0 = pbuffer.data(idx_dip_gs + 15);

    auto tr_y_xxxy_0 = pbuffer.data(idx_dip_gs + 16);

    auto tr_y_xxxz_0 = pbuffer.data(idx_dip_gs + 17);

    auto tr_y_xxyy_0 = pbuffer.data(idx_dip_gs + 18);

    auto tr_y_xxyz_0 = pbuffer.data(idx_dip_gs + 19);

    auto tr_y_xxzz_0 = pbuffer.data(idx_dip_gs + 20);

    auto tr_y_xyyy_0 = pbuffer.data(idx_dip_gs + 21);

    auto tr_y_xyyz_0 = pbuffer.data(idx_dip_gs + 22);

    auto tr_y_xyzz_0 = pbuffer.data(idx_dip_gs + 23);

    auto tr_y_xzzz_0 = pbuffer.data(idx_dip_gs + 24);

    auto tr_y_yyyy_0 = pbuffer.data(idx_dip_gs + 25);

    auto tr_y_yyyz_0 = pbuffer.data(idx_dip_gs + 26);

    auto tr_y_yyzz_0 = pbuffer.data(idx_dip_gs + 27);

    auto tr_y_yzzz_0 = pbuffer.data(idx_dip_gs + 28);

    auto tr_y_zzzz_0 = pbuffer.data(idx_dip_gs + 29);

    auto tr_z_xxxx_0 = pbuffer.data(idx_dip_gs + 30);

    auto tr_z_xxxy_0 = pbuffer.data(idx_dip_gs + 31);

    auto tr_z_xxxz_0 = pbuffer.data(idx_dip_gs + 32);

    auto tr_z_xxyy_0 = pbuffer.data(idx_dip_gs + 33);

    auto tr_z_xxyz_0 = pbuffer.data(idx_dip_gs + 34);

    auto tr_z_xxzz_0 = pbuffer.data(idx_dip_gs + 35);

    auto tr_z_xyyy_0 = pbuffer.data(idx_dip_gs + 36);

    auto tr_z_xyyz_0 = pbuffer.data(idx_dip_gs + 37);

    auto tr_z_xyzz_0 = pbuffer.data(idx_dip_gs + 38);

    auto tr_z_xzzz_0 = pbuffer.data(idx_dip_gs + 39);

    auto tr_z_yyyy_0 = pbuffer.data(idx_dip_gs + 40);

    auto tr_z_yyyz_0 = pbuffer.data(idx_dip_gs + 41);

    auto tr_z_yyzz_0 = pbuffer.data(idx_dip_gs + 42);

    auto tr_z_yzzz_0 = pbuffer.data(idx_dip_gs + 43);

    auto tr_z_zzzz_0 = pbuffer.data(idx_dip_gs + 44);

#pragma omp simd aligned(pa_x,            \
                             pa_y,        \
                             pa_z,        \
                             tr_x_xx_0,   \
                             tr_x_xxx_0,  \
                             tr_x_xxxx_0, \
                             tr_x_xxxy_0, \
                             tr_x_xxxz_0, \
                             tr_x_xxy_0,  \
                             tr_x_xxyy_0, \
                             tr_x_xxyz_0, \
                             tr_x_xxz_0,  \
                             tr_x_xxzz_0, \
                             tr_x_xyy_0,  \
                             tr_x_xyyy_0, \
                             tr_x_xyyz_0, \
                             tr_x_xyzz_0, \
                             tr_x_xzz_0,  \
                             tr_x_xzzz_0, \
                             tr_x_yy_0,   \
                             tr_x_yyy_0,  \
                             tr_x_yyyy_0, \
                             tr_x_yyyz_0, \
                             tr_x_yyzz_0, \
                             tr_x_yzz_0,  \
                             tr_x_yzzz_0, \
                             tr_x_zz_0,   \
                             tr_x_zzz_0,  \
                             tr_x_zzzz_0, \
                             tr_y_xx_0,   \
                             tr_y_xxx_0,  \
                             tr_y_xxxx_0, \
                             tr_y_xxxy_0, \
                             tr_y_xxxz_0, \
                             tr_y_xxy_0,  \
                             tr_y_xxyy_0, \
                             tr_y_xxyz_0, \
                             tr_y_xxzz_0, \
                             tr_y_xy_0,   \
                             tr_y_xyy_0,  \
                             tr_y_xyyy_0, \
                             tr_y_xyyz_0, \
                             tr_y_xyzz_0, \
                             tr_y_xzz_0,  \
                             tr_y_xzzz_0, \
                             tr_y_yy_0,   \
                             tr_y_yyy_0,  \
                             tr_y_yyyy_0, \
                             tr_y_yyyz_0, \
                             tr_y_yyz_0,  \
                             tr_y_yyzz_0, \
                             tr_y_yzz_0,  \
                             tr_y_yzzz_0, \
                             tr_y_zz_0,   \
                             tr_y_zzz_0,  \
                             tr_y_zzzz_0, \
                             tr_z_xx_0,   \
                             tr_z_xxx_0,  \
                             tr_z_xxxx_0, \
                             tr_z_xxxy_0, \
                             tr_z_xxxz_0, \
                             tr_z_xxyy_0, \
                             tr_z_xxyz_0, \
                             tr_z_xxz_0,  \
                             tr_z_xxzz_0, \
                             tr_z_xyy_0,  \
                             tr_z_xyyy_0, \
                             tr_z_xyyz_0, \
                             tr_z_xyzz_0, \
                             tr_z_xz_0,   \
                             tr_z_xzz_0,  \
                             tr_z_xzzz_0, \
                             tr_z_yy_0,   \
                             tr_z_yyy_0,  \
                             tr_z_yyyy_0, \
                             tr_z_yyyz_0, \
                             tr_z_yyz_0,  \
                             tr_z_yyzz_0, \
                             tr_z_yz_0,   \
                             tr_z_yzz_0,  \
                             tr_z_yzzz_0, \
                             tr_z_zz_0,   \
                             tr_z_zzz_0,  \
                             tr_z_zzzz_0, \
                             ts_xxx_0,    \
                             ts_yyy_0,    \
                             ts_zzz_0,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxx_0[i] = 3.0 * tr_x_xx_0[i] * fe_0 + ts_xxx_0[i] * fe_0 + tr_x_xxx_0[i] * pa_x[i];

        tr_x_xxxy_0[i] = tr_x_xxx_0[i] * pa_y[i];

        tr_x_xxxz_0[i] = tr_x_xxx_0[i] * pa_z[i];

        tr_x_xxyy_0[i] = tr_x_xx_0[i] * fe_0 + tr_x_xxy_0[i] * pa_y[i];

        tr_x_xxyz_0[i] = tr_x_xxz_0[i] * pa_y[i];

        tr_x_xxzz_0[i] = tr_x_xx_0[i] * fe_0 + tr_x_xxz_0[i] * pa_z[i];

        tr_x_xyyy_0[i] = ts_yyy_0[i] * fe_0 + tr_x_yyy_0[i] * pa_x[i];

        tr_x_xyyz_0[i] = tr_x_xyy_0[i] * pa_z[i];

        tr_x_xyzz_0[i] = tr_x_xzz_0[i] * pa_y[i];

        tr_x_xzzz_0[i] = ts_zzz_0[i] * fe_0 + tr_x_zzz_0[i] * pa_x[i];

        tr_x_yyyy_0[i] = 3.0 * tr_x_yy_0[i] * fe_0 + tr_x_yyy_0[i] * pa_y[i];

        tr_x_yyyz_0[i] = tr_x_yyy_0[i] * pa_z[i];

        tr_x_yyzz_0[i] = tr_x_zz_0[i] * fe_0 + tr_x_yzz_0[i] * pa_y[i];

        tr_x_yzzz_0[i] = tr_x_zzz_0[i] * pa_y[i];

        tr_x_zzzz_0[i] = 3.0 * tr_x_zz_0[i] * fe_0 + tr_x_zzz_0[i] * pa_z[i];

        tr_y_xxxx_0[i] = 3.0 * tr_y_xx_0[i] * fe_0 + tr_y_xxx_0[i] * pa_x[i];

        tr_y_xxxy_0[i] = 2.0 * tr_y_xy_0[i] * fe_0 + tr_y_xxy_0[i] * pa_x[i];

        tr_y_xxxz_0[i] = tr_y_xxx_0[i] * pa_z[i];

        tr_y_xxyy_0[i] = tr_y_yy_0[i] * fe_0 + tr_y_xyy_0[i] * pa_x[i];

        tr_y_xxyz_0[i] = tr_y_xxy_0[i] * pa_z[i];

        tr_y_xxzz_0[i] = tr_y_zz_0[i] * fe_0 + tr_y_xzz_0[i] * pa_x[i];

        tr_y_xyyy_0[i] = tr_y_yyy_0[i] * pa_x[i];

        tr_y_xyyz_0[i] = tr_y_yyz_0[i] * pa_x[i];

        tr_y_xyzz_0[i] = tr_y_yzz_0[i] * pa_x[i];

        tr_y_xzzz_0[i] = tr_y_zzz_0[i] * pa_x[i];

        tr_y_yyyy_0[i] = 3.0 * tr_y_yy_0[i] * fe_0 + ts_yyy_0[i] * fe_0 + tr_y_yyy_0[i] * pa_y[i];

        tr_y_yyyz_0[i] = tr_y_yyy_0[i] * pa_z[i];

        tr_y_yyzz_0[i] = tr_y_yy_0[i] * fe_0 + tr_y_yyz_0[i] * pa_z[i];

        tr_y_yzzz_0[i] = ts_zzz_0[i] * fe_0 + tr_y_zzz_0[i] * pa_y[i];

        tr_y_zzzz_0[i] = 3.0 * tr_y_zz_0[i] * fe_0 + tr_y_zzz_0[i] * pa_z[i];

        tr_z_xxxx_0[i] = 3.0 * tr_z_xx_0[i] * fe_0 + tr_z_xxx_0[i] * pa_x[i];

        tr_z_xxxy_0[i] = tr_z_xxx_0[i] * pa_y[i];

        tr_z_xxxz_0[i] = 2.0 * tr_z_xz_0[i] * fe_0 + tr_z_xxz_0[i] * pa_x[i];

        tr_z_xxyy_0[i] = tr_z_yy_0[i] * fe_0 + tr_z_xyy_0[i] * pa_x[i];

        tr_z_xxyz_0[i] = tr_z_xxz_0[i] * pa_y[i];

        tr_z_xxzz_0[i] = tr_z_zz_0[i] * fe_0 + tr_z_xzz_0[i] * pa_x[i];

        tr_z_xyyy_0[i] = tr_z_yyy_0[i] * pa_x[i];

        tr_z_xyyz_0[i] = tr_z_yyz_0[i] * pa_x[i];

        tr_z_xyzz_0[i] = tr_z_yzz_0[i] * pa_x[i];

        tr_z_xzzz_0[i] = tr_z_zzz_0[i] * pa_x[i];

        tr_z_yyyy_0[i] = 3.0 * tr_z_yy_0[i] * fe_0 + tr_z_yyy_0[i] * pa_y[i];

        tr_z_yyyz_0[i] = 2.0 * tr_z_yz_0[i] * fe_0 + tr_z_yyz_0[i] * pa_y[i];

        tr_z_yyzz_0[i] = tr_z_zz_0[i] * fe_0 + tr_z_yzz_0[i] * pa_y[i];

        tr_z_yzzz_0[i] = tr_z_zzz_0[i] * pa_y[i];

        tr_z_zzzz_0[i] = 3.0 * tr_z_zz_0[i] * fe_0 + ts_zzz_0[i] * fe_0 + tr_z_zzz_0[i] * pa_z[i];
    }
}

}  // namespace diprec

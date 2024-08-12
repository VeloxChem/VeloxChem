#include "ElectricDipoleMomentumPrimRecSD.hpp"

namespace diprec {  // diprec namespace

auto
comp_prim_electric_dipole_momentum_sd(CSimdArray<double>&       pbuffer,
                                      const size_t              idx_dip_sd,
                                      const size_t              idx_dip_ss,
                                      const size_t              idx_ovl_sp,
                                      const size_t              idx_dip_sp,
                                      const CSimdArray<double>& factors,
                                      const size_t              idx_rpb,
                                      const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PB) distances

    auto pb_x = factors.data(idx_rpb);

    auto pb_y = factors.data(idx_rpb + 1);

    auto pb_z = factors.data(idx_rpb + 2);

    // Set up components of auxiliary buffer : SS

    auto tr_x_0_0 = pbuffer.data(idx_dip_ss);

    auto tr_y_0_0 = pbuffer.data(idx_dip_ss + 1);

    auto tr_z_0_0 = pbuffer.data(idx_dip_ss + 2);

    // Set up components of auxiliary buffer : SP

    auto ts_0_x = pbuffer.data(idx_ovl_sp);

    auto ts_0_y = pbuffer.data(idx_ovl_sp + 1);

    auto ts_0_z = pbuffer.data(idx_ovl_sp + 2);

    // Set up components of auxiliary buffer : SP

    auto tr_x_0_x = pbuffer.data(idx_dip_sp);

    auto tr_x_0_y = pbuffer.data(idx_dip_sp + 1);

    auto tr_x_0_z = pbuffer.data(idx_dip_sp + 2);

    auto tr_y_0_x = pbuffer.data(idx_dip_sp + 3);

    auto tr_y_0_y = pbuffer.data(idx_dip_sp + 4);

    auto tr_y_0_z = pbuffer.data(idx_dip_sp + 5);

    auto tr_z_0_x = pbuffer.data(idx_dip_sp + 6);

    auto tr_z_0_y = pbuffer.data(idx_dip_sp + 7);

    auto tr_z_0_z = pbuffer.data(idx_dip_sp + 8);

    // Set up components of targeted buffer : SD

    auto tr_x_0_xx = pbuffer.data(idx_dip_sd);

    auto tr_x_0_xy = pbuffer.data(idx_dip_sd + 1);

    auto tr_x_0_xz = pbuffer.data(idx_dip_sd + 2);

    auto tr_x_0_yy = pbuffer.data(idx_dip_sd + 3);

    auto tr_x_0_yz = pbuffer.data(idx_dip_sd + 4);

    auto tr_x_0_zz = pbuffer.data(idx_dip_sd + 5);

    auto tr_y_0_xx = pbuffer.data(idx_dip_sd + 6);

    auto tr_y_0_xy = pbuffer.data(idx_dip_sd + 7);

    auto tr_y_0_xz = pbuffer.data(idx_dip_sd + 8);

    auto tr_y_0_yy = pbuffer.data(idx_dip_sd + 9);

    auto tr_y_0_yz = pbuffer.data(idx_dip_sd + 10);

    auto tr_y_0_zz = pbuffer.data(idx_dip_sd + 11);

    auto tr_z_0_xx = pbuffer.data(idx_dip_sd + 12);

    auto tr_z_0_xy = pbuffer.data(idx_dip_sd + 13);

    auto tr_z_0_xz = pbuffer.data(idx_dip_sd + 14);

    auto tr_z_0_yy = pbuffer.data(idx_dip_sd + 15);

    auto tr_z_0_yz = pbuffer.data(idx_dip_sd + 16);

    auto tr_z_0_zz = pbuffer.data(idx_dip_sd + 17);

#pragma omp simd aligned(pb_x,          \
                             pb_y,      \
                             pb_z,      \
                             tr_x_0_0,  \
                             tr_x_0_x,  \
                             tr_x_0_xx, \
                             tr_x_0_xy, \
                             tr_x_0_xz, \
                             tr_x_0_y,  \
                             tr_x_0_yy, \
                             tr_x_0_yz, \
                             tr_x_0_z,  \
                             tr_x_0_zz, \
                             tr_y_0_0,  \
                             tr_y_0_x,  \
                             tr_y_0_xx, \
                             tr_y_0_xy, \
                             tr_y_0_xz, \
                             tr_y_0_y,  \
                             tr_y_0_yy, \
                             tr_y_0_yz, \
                             tr_y_0_z,  \
                             tr_y_0_zz, \
                             tr_z_0_0,  \
                             tr_z_0_x,  \
                             tr_z_0_xx, \
                             tr_z_0_xy, \
                             tr_z_0_xz, \
                             tr_z_0_y,  \
                             tr_z_0_yy, \
                             tr_z_0_yz, \
                             tr_z_0_z,  \
                             tr_z_0_zz, \
                             ts_0_x,    \
                             ts_0_y,    \
                             ts_0_z,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_0_xx[i] = tr_x_0_0[i] * fe_0 + ts_0_x[i] * fe_0 + tr_x_0_x[i] * pb_x[i];

        tr_x_0_xy[i] = tr_x_0_x[i] * pb_y[i];

        tr_x_0_xz[i] = tr_x_0_x[i] * pb_z[i];

        tr_x_0_yy[i] = tr_x_0_0[i] * fe_0 + tr_x_0_y[i] * pb_y[i];

        tr_x_0_yz[i] = tr_x_0_z[i] * pb_y[i];

        tr_x_0_zz[i] = tr_x_0_0[i] * fe_0 + tr_x_0_z[i] * pb_z[i];

        tr_y_0_xx[i] = tr_y_0_0[i] * fe_0 + tr_y_0_x[i] * pb_x[i];

        tr_y_0_xy[i] = tr_y_0_y[i] * pb_x[i];

        tr_y_0_xz[i] = tr_y_0_z[i] * pb_x[i];

        tr_y_0_yy[i] = tr_y_0_0[i] * fe_0 + ts_0_y[i] * fe_0 + tr_y_0_y[i] * pb_y[i];

        tr_y_0_yz[i] = tr_y_0_y[i] * pb_z[i];

        tr_y_0_zz[i] = tr_y_0_0[i] * fe_0 + tr_y_0_z[i] * pb_z[i];

        tr_z_0_xx[i] = tr_z_0_0[i] * fe_0 + tr_z_0_x[i] * pb_x[i];

        tr_z_0_xy[i] = tr_z_0_y[i] * pb_x[i];

        tr_z_0_xz[i] = tr_z_0_z[i] * pb_x[i];

        tr_z_0_yy[i] = tr_z_0_0[i] * fe_0 + tr_z_0_y[i] * pb_y[i];

        tr_z_0_yz[i] = tr_z_0_z[i] * pb_y[i];

        tr_z_0_zz[i] = tr_z_0_0[i] * fe_0 + ts_0_z[i] * fe_0 + tr_z_0_z[i] * pb_z[i];
    }
}

}  // namespace diprec

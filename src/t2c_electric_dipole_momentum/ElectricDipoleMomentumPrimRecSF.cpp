#include "ElectricDipoleMomentumPrimRecSF.hpp"

namespace diprec { // diprec namespace

auto
comp_prim_electric_dipole_momentum_sf(CSimdArray<double>& pbuffer, 
                                      const size_t idx_dip_sf,
                                      const size_t idx_dip_sp,
                                      const size_t idx_ovl_sd,
                                      const size_t idx_dip_sd,
                                      const CSimdArray<double>& factors,
                                      const size_t idx_rpb,
                                      const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PB) distances

    auto pb_x = factors.data(idx_rpb);

    auto pb_y = factors.data(idx_rpb + 1);

    auto pb_z = factors.data(idx_rpb + 2);

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

    // Set up components of auxiliary buffer : SD

    auto ts_0_xx = pbuffer.data(idx_ovl_sd);

    auto ts_0_yy = pbuffer.data(idx_ovl_sd + 3);

    auto ts_0_zz = pbuffer.data(idx_ovl_sd + 5);

    // Set up components of auxiliary buffer : SD

    auto tr_x_0_xx = pbuffer.data(idx_dip_sd);

    auto tr_x_0_xz = pbuffer.data(idx_dip_sd + 2);

    auto tr_x_0_yy = pbuffer.data(idx_dip_sd + 3);

    auto tr_x_0_zz = pbuffer.data(idx_dip_sd + 5);

    auto tr_y_0_xx = pbuffer.data(idx_dip_sd + 6);

    auto tr_y_0_xy = pbuffer.data(idx_dip_sd + 7);

    auto tr_y_0_yy = pbuffer.data(idx_dip_sd + 9);

    auto tr_y_0_yz = pbuffer.data(idx_dip_sd + 10);

    auto tr_y_0_zz = pbuffer.data(idx_dip_sd + 11);

    auto tr_z_0_xx = pbuffer.data(idx_dip_sd + 12);

    auto tr_z_0_xz = pbuffer.data(idx_dip_sd + 14);

    auto tr_z_0_yy = pbuffer.data(idx_dip_sd + 15);

    auto tr_z_0_yz = pbuffer.data(idx_dip_sd + 16);

    auto tr_z_0_zz = pbuffer.data(idx_dip_sd + 17);

    // Set up components of targeted buffer : SF

    auto tr_x_0_xxx = pbuffer.data(idx_dip_sf);

    auto tr_x_0_xxy = pbuffer.data(idx_dip_sf + 1);

    auto tr_x_0_xxz = pbuffer.data(idx_dip_sf + 2);

    auto tr_x_0_xyy = pbuffer.data(idx_dip_sf + 3);

    auto tr_x_0_xyz = pbuffer.data(idx_dip_sf + 4);

    auto tr_x_0_xzz = pbuffer.data(idx_dip_sf + 5);

    auto tr_x_0_yyy = pbuffer.data(idx_dip_sf + 6);

    auto tr_x_0_yyz = pbuffer.data(idx_dip_sf + 7);

    auto tr_x_0_yzz = pbuffer.data(idx_dip_sf + 8);

    auto tr_x_0_zzz = pbuffer.data(idx_dip_sf + 9);

    auto tr_y_0_xxx = pbuffer.data(idx_dip_sf + 10);

    auto tr_y_0_xxy = pbuffer.data(idx_dip_sf + 11);

    auto tr_y_0_xxz = pbuffer.data(idx_dip_sf + 12);

    auto tr_y_0_xyy = pbuffer.data(idx_dip_sf + 13);

    auto tr_y_0_xyz = pbuffer.data(idx_dip_sf + 14);

    auto tr_y_0_xzz = pbuffer.data(idx_dip_sf + 15);

    auto tr_y_0_yyy = pbuffer.data(idx_dip_sf + 16);

    auto tr_y_0_yyz = pbuffer.data(idx_dip_sf + 17);

    auto tr_y_0_yzz = pbuffer.data(idx_dip_sf + 18);

    auto tr_y_0_zzz = pbuffer.data(idx_dip_sf + 19);

    auto tr_z_0_xxx = pbuffer.data(idx_dip_sf + 20);

    auto tr_z_0_xxy = pbuffer.data(idx_dip_sf + 21);

    auto tr_z_0_xxz = pbuffer.data(idx_dip_sf + 22);

    auto tr_z_0_xyy = pbuffer.data(idx_dip_sf + 23);

    auto tr_z_0_xyz = pbuffer.data(idx_dip_sf + 24);

    auto tr_z_0_xzz = pbuffer.data(idx_dip_sf + 25);

    auto tr_z_0_yyy = pbuffer.data(idx_dip_sf + 26);

    auto tr_z_0_yyz = pbuffer.data(idx_dip_sf + 27);

    auto tr_z_0_yzz = pbuffer.data(idx_dip_sf + 28);

    auto tr_z_0_zzz = pbuffer.data(idx_dip_sf + 29);

    #pragma omp simd aligned(pb_x, pb_y, pb_z, tr_x_0_x, tr_x_0_xx, tr_x_0_xxx, tr_x_0_xxy, tr_x_0_xxz, tr_x_0_xyy, tr_x_0_xyz, tr_x_0_xz, tr_x_0_xzz, tr_x_0_y, tr_x_0_yy, tr_x_0_yyy, tr_x_0_yyz, tr_x_0_yzz, tr_x_0_z, tr_x_0_zz, tr_x_0_zzz, tr_y_0_x, tr_y_0_xx, tr_y_0_xxx, tr_y_0_xxy, tr_y_0_xxz, tr_y_0_xy, tr_y_0_xyy, tr_y_0_xyz, tr_y_0_xzz, tr_y_0_y, tr_y_0_yy, tr_y_0_yyy, tr_y_0_yyz, tr_y_0_yz, tr_y_0_yzz, tr_y_0_z, tr_y_0_zz, tr_y_0_zzz, tr_z_0_x, tr_z_0_xx, tr_z_0_xxx, tr_z_0_xxy, tr_z_0_xxz, tr_z_0_xyy, tr_z_0_xyz, tr_z_0_xz, tr_z_0_xzz, tr_z_0_y, tr_z_0_yy, tr_z_0_yyy, tr_z_0_yyz, tr_z_0_yz, tr_z_0_yzz, tr_z_0_z, tr_z_0_zz, tr_z_0_zzz, ts_0_xx, ts_0_yy, ts_0_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_0_xxx[i] = 2.0 * tr_x_0_x[i] * fe_0 + ts_0_xx[i] * fe_0 + tr_x_0_xx[i] * pb_x[i];

        tr_x_0_xxy[i] = tr_x_0_xx[i] * pb_y[i];

        tr_x_0_xxz[i] = tr_x_0_xx[i] * pb_z[i];

        tr_x_0_xyy[i] = ts_0_yy[i] * fe_0 + tr_x_0_yy[i] * pb_x[i];

        tr_x_0_xyz[i] = tr_x_0_xz[i] * pb_y[i];

        tr_x_0_xzz[i] = ts_0_zz[i] * fe_0 + tr_x_0_zz[i] * pb_x[i];

        tr_x_0_yyy[i] = 2.0 * tr_x_0_y[i] * fe_0 + tr_x_0_yy[i] * pb_y[i];

        tr_x_0_yyz[i] = tr_x_0_yy[i] * pb_z[i];

        tr_x_0_yzz[i] = tr_x_0_zz[i] * pb_y[i];

        tr_x_0_zzz[i] = 2.0 * tr_x_0_z[i] * fe_0 + tr_x_0_zz[i] * pb_z[i];

        tr_y_0_xxx[i] = 2.0 * tr_y_0_x[i] * fe_0 + tr_y_0_xx[i] * pb_x[i];

        tr_y_0_xxy[i] = tr_y_0_y[i] * fe_0 + tr_y_0_xy[i] * pb_x[i];

        tr_y_0_xxz[i] = tr_y_0_xx[i] * pb_z[i];

        tr_y_0_xyy[i] = tr_y_0_yy[i] * pb_x[i];

        tr_y_0_xyz[i] = tr_y_0_yz[i] * pb_x[i];

        tr_y_0_xzz[i] = tr_y_0_zz[i] * pb_x[i];

        tr_y_0_yyy[i] = 2.0 * tr_y_0_y[i] * fe_0 + ts_0_yy[i] * fe_0 + tr_y_0_yy[i] * pb_y[i];

        tr_y_0_yyz[i] = tr_y_0_yy[i] * pb_z[i];

        tr_y_0_yzz[i] = ts_0_zz[i] * fe_0 + tr_y_0_zz[i] * pb_y[i];

        tr_y_0_zzz[i] = 2.0 * tr_y_0_z[i] * fe_0 + tr_y_0_zz[i] * pb_z[i];

        tr_z_0_xxx[i] = 2.0 * tr_z_0_x[i] * fe_0 + tr_z_0_xx[i] * pb_x[i];

        tr_z_0_xxy[i] = tr_z_0_xx[i] * pb_y[i];

        tr_z_0_xxz[i] = tr_z_0_z[i] * fe_0 + tr_z_0_xz[i] * pb_x[i];

        tr_z_0_xyy[i] = tr_z_0_yy[i] * pb_x[i];

        tr_z_0_xyz[i] = tr_z_0_yz[i] * pb_x[i];

        tr_z_0_xzz[i] = tr_z_0_zz[i] * pb_x[i];

        tr_z_0_yyy[i] = 2.0 * tr_z_0_y[i] * fe_0 + tr_z_0_yy[i] * pb_y[i];

        tr_z_0_yyz[i] = tr_z_0_z[i] * fe_0 + tr_z_0_yz[i] * pb_y[i];

        tr_z_0_yzz[i] = tr_z_0_zz[i] * pb_y[i];

        tr_z_0_zzz[i] = 2.0 * tr_z_0_z[i] * fe_0 + ts_0_zz[i] * fe_0 + tr_z_0_zz[i] * pb_z[i];
    }
}

} // diprec namespace


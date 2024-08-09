#include "ElectricDipoleMomentumPrimRecPF.hpp"

namespace diprec { // diprec namespace

auto
comp_prim_electric_dipole_momentum_pf(CSimdArray<double>& pbuffer, 
                                      const size_t idx_dip_pf,
                                      const size_t idx_dip_sd,
                                      const size_t idx_ovl_sf,
                                      const size_t idx_dip_sf,
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

    // Set up components of auxiliary buffer : SD

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

    // Set up components of auxiliary buffer : SF

    auto ts_0_xxx = pbuffer.data(idx_ovl_sf);

    auto ts_0_xxy = pbuffer.data(idx_ovl_sf + 1);

    auto ts_0_xxz = pbuffer.data(idx_ovl_sf + 2);

    auto ts_0_xyy = pbuffer.data(idx_ovl_sf + 3);

    auto ts_0_xyz = pbuffer.data(idx_ovl_sf + 4);

    auto ts_0_xzz = pbuffer.data(idx_ovl_sf + 5);

    auto ts_0_yyy = pbuffer.data(idx_ovl_sf + 6);

    auto ts_0_yyz = pbuffer.data(idx_ovl_sf + 7);

    auto ts_0_yzz = pbuffer.data(idx_ovl_sf + 8);

    auto ts_0_zzz = pbuffer.data(idx_ovl_sf + 9);

    // Set up components of auxiliary buffer : SF

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

    // Set up 0-10 components of targeted buffer : PF

    auto tr_x_x_xxx = pbuffer.data(idx_dip_pf);

    auto tr_x_x_xxy = pbuffer.data(idx_dip_pf + 1);

    auto tr_x_x_xxz = pbuffer.data(idx_dip_pf + 2);

    auto tr_x_x_xyy = pbuffer.data(idx_dip_pf + 3);

    auto tr_x_x_xyz = pbuffer.data(idx_dip_pf + 4);

    auto tr_x_x_xzz = pbuffer.data(idx_dip_pf + 5);

    auto tr_x_x_yyy = pbuffer.data(idx_dip_pf + 6);

    auto tr_x_x_yyz = pbuffer.data(idx_dip_pf + 7);

    auto tr_x_x_yzz = pbuffer.data(idx_dip_pf + 8);

    auto tr_x_x_zzz = pbuffer.data(idx_dip_pf + 9);

    #pragma omp simd aligned(pa_x, tr_x_0_xx, tr_x_0_xxx, tr_x_0_xxy, tr_x_0_xxz, tr_x_0_xy, tr_x_0_xyy, tr_x_0_xyz, tr_x_0_xz, tr_x_0_xzz, tr_x_0_yy, tr_x_0_yyy, tr_x_0_yyz, tr_x_0_yz, tr_x_0_yzz, tr_x_0_zz, tr_x_0_zzz, tr_x_x_xxx, tr_x_x_xxy, tr_x_x_xxz, tr_x_x_xyy, tr_x_x_xyz, tr_x_x_xzz, tr_x_x_yyy, tr_x_x_yyz, tr_x_x_yzz, tr_x_x_zzz, ts_0_xxx, ts_0_xxy, ts_0_xxz, ts_0_xyy, ts_0_xyz, ts_0_xzz, ts_0_yyy, ts_0_yyz, ts_0_yzz, ts_0_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_x_xxx[i] = 3.0 * tr_x_0_xx[i] * fe_0 + ts_0_xxx[i] * fe_0 + tr_x_0_xxx[i] * pa_x[i];

        tr_x_x_xxy[i] = 2.0 * tr_x_0_xy[i] * fe_0 + ts_0_xxy[i] * fe_0 + tr_x_0_xxy[i] * pa_x[i];

        tr_x_x_xxz[i] = 2.0 * tr_x_0_xz[i] * fe_0 + ts_0_xxz[i] * fe_0 + tr_x_0_xxz[i] * pa_x[i];

        tr_x_x_xyy[i] = tr_x_0_yy[i] * fe_0 + ts_0_xyy[i] * fe_0 + tr_x_0_xyy[i] * pa_x[i];

        tr_x_x_xyz[i] = tr_x_0_yz[i] * fe_0 + ts_0_xyz[i] * fe_0 + tr_x_0_xyz[i] * pa_x[i];

        tr_x_x_xzz[i] = tr_x_0_zz[i] * fe_0 + ts_0_xzz[i] * fe_0 + tr_x_0_xzz[i] * pa_x[i];

        tr_x_x_yyy[i] = ts_0_yyy[i] * fe_0 + tr_x_0_yyy[i] * pa_x[i];

        tr_x_x_yyz[i] = ts_0_yyz[i] * fe_0 + tr_x_0_yyz[i] * pa_x[i];

        tr_x_x_yzz[i] = ts_0_yzz[i] * fe_0 + tr_x_0_yzz[i] * pa_x[i];

        tr_x_x_zzz[i] = ts_0_zzz[i] * fe_0 + tr_x_0_zzz[i] * pa_x[i];
    }

    // Set up 10-20 components of targeted buffer : PF

    auto tr_x_y_xxx = pbuffer.data(idx_dip_pf + 10);

    auto tr_x_y_xxy = pbuffer.data(idx_dip_pf + 11);

    auto tr_x_y_xxz = pbuffer.data(idx_dip_pf + 12);

    auto tr_x_y_xyy = pbuffer.data(idx_dip_pf + 13);

    auto tr_x_y_xyz = pbuffer.data(idx_dip_pf + 14);

    auto tr_x_y_xzz = pbuffer.data(idx_dip_pf + 15);

    auto tr_x_y_yyy = pbuffer.data(idx_dip_pf + 16);

    auto tr_x_y_yyz = pbuffer.data(idx_dip_pf + 17);

    auto tr_x_y_yzz = pbuffer.data(idx_dip_pf + 18);

    auto tr_x_y_zzz = pbuffer.data(idx_dip_pf + 19);

    #pragma omp simd aligned(pa_y, tr_x_0_xx, tr_x_0_xxx, tr_x_0_xxy, tr_x_0_xxz, tr_x_0_xy, tr_x_0_xyy, tr_x_0_xyz, tr_x_0_xz, tr_x_0_xzz, tr_x_0_yy, tr_x_0_yyy, tr_x_0_yyz, tr_x_0_yz, tr_x_0_yzz, tr_x_0_zz, tr_x_0_zzz, tr_x_y_xxx, tr_x_y_xxy, tr_x_y_xxz, tr_x_y_xyy, tr_x_y_xyz, tr_x_y_xzz, tr_x_y_yyy, tr_x_y_yyz, tr_x_y_yzz, tr_x_y_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_y_xxx[i] = tr_x_0_xxx[i] * pa_y[i];

        tr_x_y_xxy[i] = tr_x_0_xx[i] * fe_0 + tr_x_0_xxy[i] * pa_y[i];

        tr_x_y_xxz[i] = tr_x_0_xxz[i] * pa_y[i];

        tr_x_y_xyy[i] = 2.0 * tr_x_0_xy[i] * fe_0 + tr_x_0_xyy[i] * pa_y[i];

        tr_x_y_xyz[i] = tr_x_0_xz[i] * fe_0 + tr_x_0_xyz[i] * pa_y[i];

        tr_x_y_xzz[i] = tr_x_0_xzz[i] * pa_y[i];

        tr_x_y_yyy[i] = 3.0 * tr_x_0_yy[i] * fe_0 + tr_x_0_yyy[i] * pa_y[i];

        tr_x_y_yyz[i] = 2.0 * tr_x_0_yz[i] * fe_0 + tr_x_0_yyz[i] * pa_y[i];

        tr_x_y_yzz[i] = tr_x_0_zz[i] * fe_0 + tr_x_0_yzz[i] * pa_y[i];

        tr_x_y_zzz[i] = tr_x_0_zzz[i] * pa_y[i];
    }

    // Set up 20-30 components of targeted buffer : PF

    auto tr_x_z_xxx = pbuffer.data(idx_dip_pf + 20);

    auto tr_x_z_xxy = pbuffer.data(idx_dip_pf + 21);

    auto tr_x_z_xxz = pbuffer.data(idx_dip_pf + 22);

    auto tr_x_z_xyy = pbuffer.data(idx_dip_pf + 23);

    auto tr_x_z_xyz = pbuffer.data(idx_dip_pf + 24);

    auto tr_x_z_xzz = pbuffer.data(idx_dip_pf + 25);

    auto tr_x_z_yyy = pbuffer.data(idx_dip_pf + 26);

    auto tr_x_z_yyz = pbuffer.data(idx_dip_pf + 27);

    auto tr_x_z_yzz = pbuffer.data(idx_dip_pf + 28);

    auto tr_x_z_zzz = pbuffer.data(idx_dip_pf + 29);

    #pragma omp simd aligned(pa_z, tr_x_0_xx, tr_x_0_xxx, tr_x_0_xxy, tr_x_0_xxz, tr_x_0_xy, tr_x_0_xyy, tr_x_0_xyz, tr_x_0_xz, tr_x_0_xzz, tr_x_0_yy, tr_x_0_yyy, tr_x_0_yyz, tr_x_0_yz, tr_x_0_yzz, tr_x_0_zz, tr_x_0_zzz, tr_x_z_xxx, tr_x_z_xxy, tr_x_z_xxz, tr_x_z_xyy, tr_x_z_xyz, tr_x_z_xzz, tr_x_z_yyy, tr_x_z_yyz, tr_x_z_yzz, tr_x_z_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_z_xxx[i] = tr_x_0_xxx[i] * pa_z[i];

        tr_x_z_xxy[i] = tr_x_0_xxy[i] * pa_z[i];

        tr_x_z_xxz[i] = tr_x_0_xx[i] * fe_0 + tr_x_0_xxz[i] * pa_z[i];

        tr_x_z_xyy[i] = tr_x_0_xyy[i] * pa_z[i];

        tr_x_z_xyz[i] = tr_x_0_xy[i] * fe_0 + tr_x_0_xyz[i] * pa_z[i];

        tr_x_z_xzz[i] = 2.0 * tr_x_0_xz[i] * fe_0 + tr_x_0_xzz[i] * pa_z[i];

        tr_x_z_yyy[i] = tr_x_0_yyy[i] * pa_z[i];

        tr_x_z_yyz[i] = tr_x_0_yy[i] * fe_0 + tr_x_0_yyz[i] * pa_z[i];

        tr_x_z_yzz[i] = 2.0 * tr_x_0_yz[i] * fe_0 + tr_x_0_yzz[i] * pa_z[i];

        tr_x_z_zzz[i] = 3.0 * tr_x_0_zz[i] * fe_0 + tr_x_0_zzz[i] * pa_z[i];
    }

    // Set up 30-40 components of targeted buffer : PF

    auto tr_y_x_xxx = pbuffer.data(idx_dip_pf + 30);

    auto tr_y_x_xxy = pbuffer.data(idx_dip_pf + 31);

    auto tr_y_x_xxz = pbuffer.data(idx_dip_pf + 32);

    auto tr_y_x_xyy = pbuffer.data(idx_dip_pf + 33);

    auto tr_y_x_xyz = pbuffer.data(idx_dip_pf + 34);

    auto tr_y_x_xzz = pbuffer.data(idx_dip_pf + 35);

    auto tr_y_x_yyy = pbuffer.data(idx_dip_pf + 36);

    auto tr_y_x_yyz = pbuffer.data(idx_dip_pf + 37);

    auto tr_y_x_yzz = pbuffer.data(idx_dip_pf + 38);

    auto tr_y_x_zzz = pbuffer.data(idx_dip_pf + 39);

    #pragma omp simd aligned(pa_x, tr_y_0_xx, tr_y_0_xxx, tr_y_0_xxy, tr_y_0_xxz, tr_y_0_xy, tr_y_0_xyy, tr_y_0_xyz, tr_y_0_xz, tr_y_0_xzz, tr_y_0_yy, tr_y_0_yyy, tr_y_0_yyz, tr_y_0_yz, tr_y_0_yzz, tr_y_0_zz, tr_y_0_zzz, tr_y_x_xxx, tr_y_x_xxy, tr_y_x_xxz, tr_y_x_xyy, tr_y_x_xyz, tr_y_x_xzz, tr_y_x_yyy, tr_y_x_yyz, tr_y_x_yzz, tr_y_x_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_x_xxx[i] = 3.0 * tr_y_0_xx[i] * fe_0 + tr_y_0_xxx[i] * pa_x[i];

        tr_y_x_xxy[i] = 2.0 * tr_y_0_xy[i] * fe_0 + tr_y_0_xxy[i] * pa_x[i];

        tr_y_x_xxz[i] = 2.0 * tr_y_0_xz[i] * fe_0 + tr_y_0_xxz[i] * pa_x[i];

        tr_y_x_xyy[i] = tr_y_0_yy[i] * fe_0 + tr_y_0_xyy[i] * pa_x[i];

        tr_y_x_xyz[i] = tr_y_0_yz[i] * fe_0 + tr_y_0_xyz[i] * pa_x[i];

        tr_y_x_xzz[i] = tr_y_0_zz[i] * fe_0 + tr_y_0_xzz[i] * pa_x[i];

        tr_y_x_yyy[i] = tr_y_0_yyy[i] * pa_x[i];

        tr_y_x_yyz[i] = tr_y_0_yyz[i] * pa_x[i];

        tr_y_x_yzz[i] = tr_y_0_yzz[i] * pa_x[i];

        tr_y_x_zzz[i] = tr_y_0_zzz[i] * pa_x[i];
    }

    // Set up 40-50 components of targeted buffer : PF

    auto tr_y_y_xxx = pbuffer.data(idx_dip_pf + 40);

    auto tr_y_y_xxy = pbuffer.data(idx_dip_pf + 41);

    auto tr_y_y_xxz = pbuffer.data(idx_dip_pf + 42);

    auto tr_y_y_xyy = pbuffer.data(idx_dip_pf + 43);

    auto tr_y_y_xyz = pbuffer.data(idx_dip_pf + 44);

    auto tr_y_y_xzz = pbuffer.data(idx_dip_pf + 45);

    auto tr_y_y_yyy = pbuffer.data(idx_dip_pf + 46);

    auto tr_y_y_yyz = pbuffer.data(idx_dip_pf + 47);

    auto tr_y_y_yzz = pbuffer.data(idx_dip_pf + 48);

    auto tr_y_y_zzz = pbuffer.data(idx_dip_pf + 49);

    #pragma omp simd aligned(pa_y, tr_y_0_xx, tr_y_0_xxx, tr_y_0_xxy, tr_y_0_xxz, tr_y_0_xy, tr_y_0_xyy, tr_y_0_xyz, tr_y_0_xz, tr_y_0_xzz, tr_y_0_yy, tr_y_0_yyy, tr_y_0_yyz, tr_y_0_yz, tr_y_0_yzz, tr_y_0_zz, tr_y_0_zzz, tr_y_y_xxx, tr_y_y_xxy, tr_y_y_xxz, tr_y_y_xyy, tr_y_y_xyz, tr_y_y_xzz, tr_y_y_yyy, tr_y_y_yyz, tr_y_y_yzz, tr_y_y_zzz, ts_0_xxx, ts_0_xxy, ts_0_xxz, ts_0_xyy, ts_0_xyz, ts_0_xzz, ts_0_yyy, ts_0_yyz, ts_0_yzz, ts_0_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_y_xxx[i] = ts_0_xxx[i] * fe_0 + tr_y_0_xxx[i] * pa_y[i];

        tr_y_y_xxy[i] = tr_y_0_xx[i] * fe_0 + ts_0_xxy[i] * fe_0 + tr_y_0_xxy[i] * pa_y[i];

        tr_y_y_xxz[i] = ts_0_xxz[i] * fe_0 + tr_y_0_xxz[i] * pa_y[i];

        tr_y_y_xyy[i] = 2.0 * tr_y_0_xy[i] * fe_0 + ts_0_xyy[i] * fe_0 + tr_y_0_xyy[i] * pa_y[i];

        tr_y_y_xyz[i] = tr_y_0_xz[i] * fe_0 + ts_0_xyz[i] * fe_0 + tr_y_0_xyz[i] * pa_y[i];

        tr_y_y_xzz[i] = ts_0_xzz[i] * fe_0 + tr_y_0_xzz[i] * pa_y[i];

        tr_y_y_yyy[i] = 3.0 * tr_y_0_yy[i] * fe_0 + ts_0_yyy[i] * fe_0 + tr_y_0_yyy[i] * pa_y[i];

        tr_y_y_yyz[i] = 2.0 * tr_y_0_yz[i] * fe_0 + ts_0_yyz[i] * fe_0 + tr_y_0_yyz[i] * pa_y[i];

        tr_y_y_yzz[i] = tr_y_0_zz[i] * fe_0 + ts_0_yzz[i] * fe_0 + tr_y_0_yzz[i] * pa_y[i];

        tr_y_y_zzz[i] = ts_0_zzz[i] * fe_0 + tr_y_0_zzz[i] * pa_y[i];
    }

    // Set up 50-60 components of targeted buffer : PF

    auto tr_y_z_xxx = pbuffer.data(idx_dip_pf + 50);

    auto tr_y_z_xxy = pbuffer.data(idx_dip_pf + 51);

    auto tr_y_z_xxz = pbuffer.data(idx_dip_pf + 52);

    auto tr_y_z_xyy = pbuffer.data(idx_dip_pf + 53);

    auto tr_y_z_xyz = pbuffer.data(idx_dip_pf + 54);

    auto tr_y_z_xzz = pbuffer.data(idx_dip_pf + 55);

    auto tr_y_z_yyy = pbuffer.data(idx_dip_pf + 56);

    auto tr_y_z_yyz = pbuffer.data(idx_dip_pf + 57);

    auto tr_y_z_yzz = pbuffer.data(idx_dip_pf + 58);

    auto tr_y_z_zzz = pbuffer.data(idx_dip_pf + 59);

    #pragma omp simd aligned(pa_z, tr_y_0_xx, tr_y_0_xxx, tr_y_0_xxy, tr_y_0_xxz, tr_y_0_xy, tr_y_0_xyy, tr_y_0_xyz, tr_y_0_xz, tr_y_0_xzz, tr_y_0_yy, tr_y_0_yyy, tr_y_0_yyz, tr_y_0_yz, tr_y_0_yzz, tr_y_0_zz, tr_y_0_zzz, tr_y_z_xxx, tr_y_z_xxy, tr_y_z_xxz, tr_y_z_xyy, tr_y_z_xyz, tr_y_z_xzz, tr_y_z_yyy, tr_y_z_yyz, tr_y_z_yzz, tr_y_z_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_z_xxx[i] = tr_y_0_xxx[i] * pa_z[i];

        tr_y_z_xxy[i] = tr_y_0_xxy[i] * pa_z[i];

        tr_y_z_xxz[i] = tr_y_0_xx[i] * fe_0 + tr_y_0_xxz[i] * pa_z[i];

        tr_y_z_xyy[i] = tr_y_0_xyy[i] * pa_z[i];

        tr_y_z_xyz[i] = tr_y_0_xy[i] * fe_0 + tr_y_0_xyz[i] * pa_z[i];

        tr_y_z_xzz[i] = 2.0 * tr_y_0_xz[i] * fe_0 + tr_y_0_xzz[i] * pa_z[i];

        tr_y_z_yyy[i] = tr_y_0_yyy[i] * pa_z[i];

        tr_y_z_yyz[i] = tr_y_0_yy[i] * fe_0 + tr_y_0_yyz[i] * pa_z[i];

        tr_y_z_yzz[i] = 2.0 * tr_y_0_yz[i] * fe_0 + tr_y_0_yzz[i] * pa_z[i];

        tr_y_z_zzz[i] = 3.0 * tr_y_0_zz[i] * fe_0 + tr_y_0_zzz[i] * pa_z[i];
    }

    // Set up 60-70 components of targeted buffer : PF

    auto tr_z_x_xxx = pbuffer.data(idx_dip_pf + 60);

    auto tr_z_x_xxy = pbuffer.data(idx_dip_pf + 61);

    auto tr_z_x_xxz = pbuffer.data(idx_dip_pf + 62);

    auto tr_z_x_xyy = pbuffer.data(idx_dip_pf + 63);

    auto tr_z_x_xyz = pbuffer.data(idx_dip_pf + 64);

    auto tr_z_x_xzz = pbuffer.data(idx_dip_pf + 65);

    auto tr_z_x_yyy = pbuffer.data(idx_dip_pf + 66);

    auto tr_z_x_yyz = pbuffer.data(idx_dip_pf + 67);

    auto tr_z_x_yzz = pbuffer.data(idx_dip_pf + 68);

    auto tr_z_x_zzz = pbuffer.data(idx_dip_pf + 69);

    #pragma omp simd aligned(pa_x, tr_z_0_xx, tr_z_0_xxx, tr_z_0_xxy, tr_z_0_xxz, tr_z_0_xy, tr_z_0_xyy, tr_z_0_xyz, tr_z_0_xz, tr_z_0_xzz, tr_z_0_yy, tr_z_0_yyy, tr_z_0_yyz, tr_z_0_yz, tr_z_0_yzz, tr_z_0_zz, tr_z_0_zzz, tr_z_x_xxx, tr_z_x_xxy, tr_z_x_xxz, tr_z_x_xyy, tr_z_x_xyz, tr_z_x_xzz, tr_z_x_yyy, tr_z_x_yyz, tr_z_x_yzz, tr_z_x_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_x_xxx[i] = 3.0 * tr_z_0_xx[i] * fe_0 + tr_z_0_xxx[i] * pa_x[i];

        tr_z_x_xxy[i] = 2.0 * tr_z_0_xy[i] * fe_0 + tr_z_0_xxy[i] * pa_x[i];

        tr_z_x_xxz[i] = 2.0 * tr_z_0_xz[i] * fe_0 + tr_z_0_xxz[i] * pa_x[i];

        tr_z_x_xyy[i] = tr_z_0_yy[i] * fe_0 + tr_z_0_xyy[i] * pa_x[i];

        tr_z_x_xyz[i] = tr_z_0_yz[i] * fe_0 + tr_z_0_xyz[i] * pa_x[i];

        tr_z_x_xzz[i] = tr_z_0_zz[i] * fe_0 + tr_z_0_xzz[i] * pa_x[i];

        tr_z_x_yyy[i] = tr_z_0_yyy[i] * pa_x[i];

        tr_z_x_yyz[i] = tr_z_0_yyz[i] * pa_x[i];

        tr_z_x_yzz[i] = tr_z_0_yzz[i] * pa_x[i];

        tr_z_x_zzz[i] = tr_z_0_zzz[i] * pa_x[i];
    }

    // Set up 70-80 components of targeted buffer : PF

    auto tr_z_y_xxx = pbuffer.data(idx_dip_pf + 70);

    auto tr_z_y_xxy = pbuffer.data(idx_dip_pf + 71);

    auto tr_z_y_xxz = pbuffer.data(idx_dip_pf + 72);

    auto tr_z_y_xyy = pbuffer.data(idx_dip_pf + 73);

    auto tr_z_y_xyz = pbuffer.data(idx_dip_pf + 74);

    auto tr_z_y_xzz = pbuffer.data(idx_dip_pf + 75);

    auto tr_z_y_yyy = pbuffer.data(idx_dip_pf + 76);

    auto tr_z_y_yyz = pbuffer.data(idx_dip_pf + 77);

    auto tr_z_y_yzz = pbuffer.data(idx_dip_pf + 78);

    auto tr_z_y_zzz = pbuffer.data(idx_dip_pf + 79);

    #pragma omp simd aligned(pa_y, tr_z_0_xx, tr_z_0_xxx, tr_z_0_xxy, tr_z_0_xxz, tr_z_0_xy, tr_z_0_xyy, tr_z_0_xyz, tr_z_0_xz, tr_z_0_xzz, tr_z_0_yy, tr_z_0_yyy, tr_z_0_yyz, tr_z_0_yz, tr_z_0_yzz, tr_z_0_zz, tr_z_0_zzz, tr_z_y_xxx, tr_z_y_xxy, tr_z_y_xxz, tr_z_y_xyy, tr_z_y_xyz, tr_z_y_xzz, tr_z_y_yyy, tr_z_y_yyz, tr_z_y_yzz, tr_z_y_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_y_xxx[i] = tr_z_0_xxx[i] * pa_y[i];

        tr_z_y_xxy[i] = tr_z_0_xx[i] * fe_0 + tr_z_0_xxy[i] * pa_y[i];

        tr_z_y_xxz[i] = tr_z_0_xxz[i] * pa_y[i];

        tr_z_y_xyy[i] = 2.0 * tr_z_0_xy[i] * fe_0 + tr_z_0_xyy[i] * pa_y[i];

        tr_z_y_xyz[i] = tr_z_0_xz[i] * fe_0 + tr_z_0_xyz[i] * pa_y[i];

        tr_z_y_xzz[i] = tr_z_0_xzz[i] * pa_y[i];

        tr_z_y_yyy[i] = 3.0 * tr_z_0_yy[i] * fe_0 + tr_z_0_yyy[i] * pa_y[i];

        tr_z_y_yyz[i] = 2.0 * tr_z_0_yz[i] * fe_0 + tr_z_0_yyz[i] * pa_y[i];

        tr_z_y_yzz[i] = tr_z_0_zz[i] * fe_0 + tr_z_0_yzz[i] * pa_y[i];

        tr_z_y_zzz[i] = tr_z_0_zzz[i] * pa_y[i];
    }

    // Set up 80-90 components of targeted buffer : PF

    auto tr_z_z_xxx = pbuffer.data(idx_dip_pf + 80);

    auto tr_z_z_xxy = pbuffer.data(idx_dip_pf + 81);

    auto tr_z_z_xxz = pbuffer.data(idx_dip_pf + 82);

    auto tr_z_z_xyy = pbuffer.data(idx_dip_pf + 83);

    auto tr_z_z_xyz = pbuffer.data(idx_dip_pf + 84);

    auto tr_z_z_xzz = pbuffer.data(idx_dip_pf + 85);

    auto tr_z_z_yyy = pbuffer.data(idx_dip_pf + 86);

    auto tr_z_z_yyz = pbuffer.data(idx_dip_pf + 87);

    auto tr_z_z_yzz = pbuffer.data(idx_dip_pf + 88);

    auto tr_z_z_zzz = pbuffer.data(idx_dip_pf + 89);

    #pragma omp simd aligned(pa_z, tr_z_0_xx, tr_z_0_xxx, tr_z_0_xxy, tr_z_0_xxz, tr_z_0_xy, tr_z_0_xyy, tr_z_0_xyz, tr_z_0_xz, tr_z_0_xzz, tr_z_0_yy, tr_z_0_yyy, tr_z_0_yyz, tr_z_0_yz, tr_z_0_yzz, tr_z_0_zz, tr_z_0_zzz, tr_z_z_xxx, tr_z_z_xxy, tr_z_z_xxz, tr_z_z_xyy, tr_z_z_xyz, tr_z_z_xzz, tr_z_z_yyy, tr_z_z_yyz, tr_z_z_yzz, tr_z_z_zzz, ts_0_xxx, ts_0_xxy, ts_0_xxz, ts_0_xyy, ts_0_xyz, ts_0_xzz, ts_0_yyy, ts_0_yyz, ts_0_yzz, ts_0_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_z_xxx[i] = ts_0_xxx[i] * fe_0 + tr_z_0_xxx[i] * pa_z[i];

        tr_z_z_xxy[i] = ts_0_xxy[i] * fe_0 + tr_z_0_xxy[i] * pa_z[i];

        tr_z_z_xxz[i] = tr_z_0_xx[i] * fe_0 + ts_0_xxz[i] * fe_0 + tr_z_0_xxz[i] * pa_z[i];

        tr_z_z_xyy[i] = ts_0_xyy[i] * fe_0 + tr_z_0_xyy[i] * pa_z[i];

        tr_z_z_xyz[i] = tr_z_0_xy[i] * fe_0 + ts_0_xyz[i] * fe_0 + tr_z_0_xyz[i] * pa_z[i];

        tr_z_z_xzz[i] = 2.0 * tr_z_0_xz[i] * fe_0 + ts_0_xzz[i] * fe_0 + tr_z_0_xzz[i] * pa_z[i];

        tr_z_z_yyy[i] = ts_0_yyy[i] * fe_0 + tr_z_0_yyy[i] * pa_z[i];

        tr_z_z_yyz[i] = tr_z_0_yy[i] * fe_0 + ts_0_yyz[i] * fe_0 + tr_z_0_yyz[i] * pa_z[i];

        tr_z_z_yzz[i] = 2.0 * tr_z_0_yz[i] * fe_0 + ts_0_yzz[i] * fe_0 + tr_z_0_yzz[i] * pa_z[i];

        tr_z_z_zzz[i] = 3.0 * tr_z_0_zz[i] * fe_0 + ts_0_zzz[i] * fe_0 + tr_z_0_zzz[i] * pa_z[i];
    }

}

} // diprec namespace


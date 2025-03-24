#include "ThreeCenterOverlapPrimRecPF.hpp"

namespace t3ovlrec { // t3ovlrec namespace

auto
comp_prim_overlap_pf(CSimdArray<double>& pbuffer, 
                     const size_t idx_pf,
                     const size_t idx_sd,
                     const size_t idx_sf,
                     const CSimdArray<double>& factors,
                     const size_t idx_rga,
                     const double a_exp,
                     const double c_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(GA) distances

    auto ga_x = factors.data(idx_rga);

    auto ga_y = factors.data(idx_rga + 1);

    auto ga_z = factors.data(idx_rga + 2);

    // Set up components of auxiliary buffer : SD

    auto ts_0_xx = pbuffer.data(idx_sd);

    auto ts_0_xy = pbuffer.data(idx_sd + 1);

    auto ts_0_xz = pbuffer.data(idx_sd + 2);

    auto ts_0_yy = pbuffer.data(idx_sd + 3);

    auto ts_0_yz = pbuffer.data(idx_sd + 4);

    auto ts_0_zz = pbuffer.data(idx_sd + 5);

    // Set up components of auxiliary buffer : SF

    auto ts_0_xxx = pbuffer.data(idx_sf);

    auto ts_0_xxy = pbuffer.data(idx_sf + 1);

    auto ts_0_xxz = pbuffer.data(idx_sf + 2);

    auto ts_0_xyy = pbuffer.data(idx_sf + 3);

    auto ts_0_xyz = pbuffer.data(idx_sf + 4);

    auto ts_0_xzz = pbuffer.data(idx_sf + 5);

    auto ts_0_yyy = pbuffer.data(idx_sf + 6);

    auto ts_0_yyz = pbuffer.data(idx_sf + 7);

    auto ts_0_yzz = pbuffer.data(idx_sf + 8);

    auto ts_0_zzz = pbuffer.data(idx_sf + 9);

    // Set up 0-10 components of targeted buffer : PF

    auto ts_x_xxx = pbuffer.data(idx_pf);

    auto ts_x_xxy = pbuffer.data(idx_pf + 1);

    auto ts_x_xxz = pbuffer.data(idx_pf + 2);

    auto ts_x_xyy = pbuffer.data(idx_pf + 3);

    auto ts_x_xyz = pbuffer.data(idx_pf + 4);

    auto ts_x_xzz = pbuffer.data(idx_pf + 5);

    auto ts_x_yyy = pbuffer.data(idx_pf + 6);

    auto ts_x_yyz = pbuffer.data(idx_pf + 7);

    auto ts_x_yzz = pbuffer.data(idx_pf + 8);

    auto ts_x_zzz = pbuffer.data(idx_pf + 9);

    #pragma omp simd aligned(ga_x, ts_0_xx, ts_0_xxx, ts_0_xxy, ts_0_xxz, ts_0_xy, ts_0_xyy, ts_0_xyz, ts_0_xz, ts_0_xzz, ts_0_yy, ts_0_yyy, ts_0_yyz, ts_0_yz, ts_0_yzz, ts_0_zz, ts_0_zzz, ts_x_xxx, ts_x_xxy, ts_x_xxz, ts_x_xyy, ts_x_xyz, ts_x_xzz, ts_x_yyy, ts_x_yyz, ts_x_yzz, ts_x_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_x_xxx[i] = 3.0 * ts_0_xx[i] * gfe_0 + ts_0_xxx[i] * ga_x[i];

        ts_x_xxy[i] = 2.0 * ts_0_xy[i] * gfe_0 + ts_0_xxy[i] * ga_x[i];

        ts_x_xxz[i] = 2.0 * ts_0_xz[i] * gfe_0 + ts_0_xxz[i] * ga_x[i];

        ts_x_xyy[i] = ts_0_yy[i] * gfe_0 + ts_0_xyy[i] * ga_x[i];

        ts_x_xyz[i] = ts_0_yz[i] * gfe_0 + ts_0_xyz[i] * ga_x[i];

        ts_x_xzz[i] = ts_0_zz[i] * gfe_0 + ts_0_xzz[i] * ga_x[i];

        ts_x_yyy[i] = ts_0_yyy[i] * ga_x[i];

        ts_x_yyz[i] = ts_0_yyz[i] * ga_x[i];

        ts_x_yzz[i] = ts_0_yzz[i] * ga_x[i];

        ts_x_zzz[i] = ts_0_zzz[i] * ga_x[i];
    }

    // Set up 10-20 components of targeted buffer : PF

    auto ts_y_xxx = pbuffer.data(idx_pf + 10);

    auto ts_y_xxy = pbuffer.data(idx_pf + 11);

    auto ts_y_xxz = pbuffer.data(idx_pf + 12);

    auto ts_y_xyy = pbuffer.data(idx_pf + 13);

    auto ts_y_xyz = pbuffer.data(idx_pf + 14);

    auto ts_y_xzz = pbuffer.data(idx_pf + 15);

    auto ts_y_yyy = pbuffer.data(idx_pf + 16);

    auto ts_y_yyz = pbuffer.data(idx_pf + 17);

    auto ts_y_yzz = pbuffer.data(idx_pf + 18);

    auto ts_y_zzz = pbuffer.data(idx_pf + 19);

    #pragma omp simd aligned(ga_y, ts_0_xx, ts_0_xxx, ts_0_xxy, ts_0_xxz, ts_0_xy, ts_0_xyy, ts_0_xyz, ts_0_xz, ts_0_xzz, ts_0_yy, ts_0_yyy, ts_0_yyz, ts_0_yz, ts_0_yzz, ts_0_zz, ts_0_zzz, ts_y_xxx, ts_y_xxy, ts_y_xxz, ts_y_xyy, ts_y_xyz, ts_y_xzz, ts_y_yyy, ts_y_yyz, ts_y_yzz, ts_y_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_y_xxx[i] = ts_0_xxx[i] * ga_y[i];

        ts_y_xxy[i] = ts_0_xx[i] * gfe_0 + ts_0_xxy[i] * ga_y[i];

        ts_y_xxz[i] = ts_0_xxz[i] * ga_y[i];

        ts_y_xyy[i] = 2.0 * ts_0_xy[i] * gfe_0 + ts_0_xyy[i] * ga_y[i];

        ts_y_xyz[i] = ts_0_xz[i] * gfe_0 + ts_0_xyz[i] * ga_y[i];

        ts_y_xzz[i] = ts_0_xzz[i] * ga_y[i];

        ts_y_yyy[i] = 3.0 * ts_0_yy[i] * gfe_0 + ts_0_yyy[i] * ga_y[i];

        ts_y_yyz[i] = 2.0 * ts_0_yz[i] * gfe_0 + ts_0_yyz[i] * ga_y[i];

        ts_y_yzz[i] = ts_0_zz[i] * gfe_0 + ts_0_yzz[i] * ga_y[i];

        ts_y_zzz[i] = ts_0_zzz[i] * ga_y[i];
    }

    // Set up 20-30 components of targeted buffer : PF

    auto ts_z_xxx = pbuffer.data(idx_pf + 20);

    auto ts_z_xxy = pbuffer.data(idx_pf + 21);

    auto ts_z_xxz = pbuffer.data(idx_pf + 22);

    auto ts_z_xyy = pbuffer.data(idx_pf + 23);

    auto ts_z_xyz = pbuffer.data(idx_pf + 24);

    auto ts_z_xzz = pbuffer.data(idx_pf + 25);

    auto ts_z_yyy = pbuffer.data(idx_pf + 26);

    auto ts_z_yyz = pbuffer.data(idx_pf + 27);

    auto ts_z_yzz = pbuffer.data(idx_pf + 28);

    auto ts_z_zzz = pbuffer.data(idx_pf + 29);

    #pragma omp simd aligned(ga_z, ts_0_xx, ts_0_xxx, ts_0_xxy, ts_0_xxz, ts_0_xy, ts_0_xyy, ts_0_xyz, ts_0_xz, ts_0_xzz, ts_0_yy, ts_0_yyy, ts_0_yyz, ts_0_yz, ts_0_yzz, ts_0_zz, ts_0_zzz, ts_z_xxx, ts_z_xxy, ts_z_xxz, ts_z_xyy, ts_z_xyz, ts_z_xzz, ts_z_yyy, ts_z_yyz, ts_z_yzz, ts_z_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_z_xxx[i] = ts_0_xxx[i] * ga_z[i];

        ts_z_xxy[i] = ts_0_xxy[i] * ga_z[i];

        ts_z_xxz[i] = ts_0_xx[i] * gfe_0 + ts_0_xxz[i] * ga_z[i];

        ts_z_xyy[i] = ts_0_xyy[i] * ga_z[i];

        ts_z_xyz[i] = ts_0_xy[i] * gfe_0 + ts_0_xyz[i] * ga_z[i];

        ts_z_xzz[i] = 2.0 * ts_0_xz[i] * gfe_0 + ts_0_xzz[i] * ga_z[i];

        ts_z_yyy[i] = ts_0_yyy[i] * ga_z[i];

        ts_z_yyz[i] = ts_0_yy[i] * gfe_0 + ts_0_yyz[i] * ga_z[i];

        ts_z_yzz[i] = 2.0 * ts_0_yz[i] * gfe_0 + ts_0_yzz[i] * ga_z[i];

        ts_z_zzz[i] = 3.0 * ts_0_zz[i] * gfe_0 + ts_0_zzz[i] * ga_z[i];
    }

}

} // t3ovlrec namespace


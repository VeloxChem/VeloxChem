#include "T2CHrrABRecDD.hpp"

namespace t2chrr { // t2chrr namespace

auto
comp_hrr_dd(CSimdArray<double>& cbuffer, 
            const size_t idx_dd,
            const size_t idx_pd,
            const size_t idx_pf,
            const CSimdArray<double>& factors) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    // Set up R(AB) distances

    auto ab_x = factors.data(3);

    auto ab_y = factors.data(4);

    auto ab_z = factors.data(5);

    // Set up components of auxiliary buffer : PD

    auto t_x_xx = cbuffer.data(idx_pd);

    auto t_x_xy = cbuffer.data(idx_pd + 1);

    auto t_x_xz = cbuffer.data(idx_pd + 2);

    auto t_x_yy = cbuffer.data(idx_pd + 3);

    auto t_x_yz = cbuffer.data(idx_pd + 4);

    auto t_x_zz = cbuffer.data(idx_pd + 5);

    auto t_y_xx = cbuffer.data(idx_pd + 6);

    auto t_y_xy = cbuffer.data(idx_pd + 7);

    auto t_y_xz = cbuffer.data(idx_pd + 8);

    auto t_y_yy = cbuffer.data(idx_pd + 9);

    auto t_y_yz = cbuffer.data(idx_pd + 10);

    auto t_y_zz = cbuffer.data(idx_pd + 11);

    auto t_z_xx = cbuffer.data(idx_pd + 12);

    auto t_z_xy = cbuffer.data(idx_pd + 13);

    auto t_z_xz = cbuffer.data(idx_pd + 14);

    auto t_z_yy = cbuffer.data(idx_pd + 15);

    auto t_z_yz = cbuffer.data(idx_pd + 16);

    auto t_z_zz = cbuffer.data(idx_pd + 17);

    // Set up components of auxiliary buffer : PF

    auto t_x_xxx = cbuffer.data(idx_pf);

    auto t_x_xxy = cbuffer.data(idx_pf + 1);

    auto t_x_xxz = cbuffer.data(idx_pf + 2);

    auto t_x_xyy = cbuffer.data(idx_pf + 3);

    auto t_x_xyz = cbuffer.data(idx_pf + 4);

    auto t_x_xzz = cbuffer.data(idx_pf + 5);

    auto t_y_xxx = cbuffer.data(idx_pf + 10);

    auto t_y_xxy = cbuffer.data(idx_pf + 11);

    auto t_y_xxz = cbuffer.data(idx_pf + 12);

    auto t_y_xyy = cbuffer.data(idx_pf + 13);

    auto t_y_xyz = cbuffer.data(idx_pf + 14);

    auto t_y_xzz = cbuffer.data(idx_pf + 15);

    auto t_y_yyy = cbuffer.data(idx_pf + 16);

    auto t_y_yyz = cbuffer.data(idx_pf + 17);

    auto t_y_yzz = cbuffer.data(idx_pf + 18);

    auto t_z_xxx = cbuffer.data(idx_pf + 20);

    auto t_z_xxy = cbuffer.data(idx_pf + 21);

    auto t_z_xxz = cbuffer.data(idx_pf + 22);

    auto t_z_xyy = cbuffer.data(idx_pf + 23);

    auto t_z_xyz = cbuffer.data(idx_pf + 24);

    auto t_z_xzz = cbuffer.data(idx_pf + 25);

    auto t_z_yyy = cbuffer.data(idx_pf + 26);

    auto t_z_yyz = cbuffer.data(idx_pf + 27);

    auto t_z_yzz = cbuffer.data(idx_pf + 28);

    auto t_z_zzz = cbuffer.data(idx_pf + 29);

    // Set up components of targeted buffer : DD

    auto t_xx_xx = cbuffer.data(idx_dd);

    auto t_xx_xy = cbuffer.data(idx_dd + 1);

    auto t_xx_xz = cbuffer.data(idx_dd + 2);

    auto t_xx_yy = cbuffer.data(idx_dd + 3);

    auto t_xx_yz = cbuffer.data(idx_dd + 4);

    auto t_xx_zz = cbuffer.data(idx_dd + 5);

    auto t_xy_xx = cbuffer.data(idx_dd + 6);

    auto t_xy_xy = cbuffer.data(idx_dd + 7);

    auto t_xy_xz = cbuffer.data(idx_dd + 8);

    auto t_xy_yy = cbuffer.data(idx_dd + 9);

    auto t_xy_yz = cbuffer.data(idx_dd + 10);

    auto t_xy_zz = cbuffer.data(idx_dd + 11);

    auto t_xz_xx = cbuffer.data(idx_dd + 12);

    auto t_xz_xy = cbuffer.data(idx_dd + 13);

    auto t_xz_xz = cbuffer.data(idx_dd + 14);

    auto t_xz_yy = cbuffer.data(idx_dd + 15);

    auto t_xz_yz = cbuffer.data(idx_dd + 16);

    auto t_xz_zz = cbuffer.data(idx_dd + 17);

    auto t_yy_xx = cbuffer.data(idx_dd + 18);

    auto t_yy_xy = cbuffer.data(idx_dd + 19);

    auto t_yy_xz = cbuffer.data(idx_dd + 20);

    auto t_yy_yy = cbuffer.data(idx_dd + 21);

    auto t_yy_yz = cbuffer.data(idx_dd + 22);

    auto t_yy_zz = cbuffer.data(idx_dd + 23);

    auto t_yz_xx = cbuffer.data(idx_dd + 24);

    auto t_yz_xy = cbuffer.data(idx_dd + 25);

    auto t_yz_xz = cbuffer.data(idx_dd + 26);

    auto t_yz_yy = cbuffer.data(idx_dd + 27);

    auto t_yz_yz = cbuffer.data(idx_dd + 28);

    auto t_yz_zz = cbuffer.data(idx_dd + 29);

    auto t_zz_xx = cbuffer.data(idx_dd + 30);

    auto t_zz_xy = cbuffer.data(idx_dd + 31);

    auto t_zz_xz = cbuffer.data(idx_dd + 32);

    auto t_zz_yy = cbuffer.data(idx_dd + 33);

    auto t_zz_yz = cbuffer.data(idx_dd + 34);

    auto t_zz_zz = cbuffer.data(idx_dd + 35);

    #pragma omp simd aligned(ab_x, ab_y, ab_z, t_x_xx, t_x_xxx, t_x_xxy, t_x_xxz, t_x_xy, t_x_xyy, t_x_xyz, t_x_xz, t_x_xzz, t_x_yy, t_x_yz, t_x_zz, t_xx_xx, t_xx_xy, t_xx_xz, t_xx_yy, t_xx_yz, t_xx_zz, t_xy_xx, t_xy_xy, t_xy_xz, t_xy_yy, t_xy_yz, t_xy_zz, t_xz_xx, t_xz_xy, t_xz_xz, t_xz_yy, t_xz_yz, t_xz_zz, t_y_xx, t_y_xxx, t_y_xxy, t_y_xxz, t_y_xy, t_y_xyy, t_y_xyz, t_y_xz, t_y_xzz, t_y_yy, t_y_yyy, t_y_yyz, t_y_yz, t_y_yzz, t_y_zz, t_yy_xx, t_yy_xy, t_yy_xz, t_yy_yy, t_yy_yz, t_yy_zz, t_yz_xx, t_yz_xy, t_yz_xz, t_yz_yy, t_yz_yz, t_yz_zz, t_z_xx, t_z_xxx, t_z_xxy, t_z_xxz, t_z_xy, t_z_xyy, t_z_xyz, t_z_xz, t_z_xzz, t_z_yy, t_z_yyy, t_z_yyz, t_z_yz, t_z_yzz, t_z_zz, t_z_zzz, t_zz_xx, t_zz_xy, t_zz_xz, t_zz_yy, t_zz_yz, t_zz_zz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        t_xx_xx[i] = -t_x_xx[i] * ab_x[i] + t_x_xxx[i];

        t_xx_xy[i] = -t_x_xy[i] * ab_x[i] + t_x_xxy[i];

        t_xx_xz[i] = -t_x_xz[i] * ab_x[i] + t_x_xxz[i];

        t_xx_yy[i] = -t_x_yy[i] * ab_x[i] + t_x_xyy[i];

        t_xx_yz[i] = -t_x_yz[i] * ab_x[i] + t_x_xyz[i];

        t_xx_zz[i] = -t_x_zz[i] * ab_x[i] + t_x_xzz[i];

        t_xy_xx[i] = -t_y_xx[i] * ab_x[i] + t_y_xxx[i];

        t_xy_xy[i] = -t_y_xy[i] * ab_x[i] + t_y_xxy[i];

        t_xy_xz[i] = -t_y_xz[i] * ab_x[i] + t_y_xxz[i];

        t_xy_yy[i] = -t_y_yy[i] * ab_x[i] + t_y_xyy[i];

        t_xy_yz[i] = -t_y_yz[i] * ab_x[i] + t_y_xyz[i];

        t_xy_zz[i] = -t_y_zz[i] * ab_x[i] + t_y_xzz[i];

        t_xz_xx[i] = -t_z_xx[i] * ab_x[i] + t_z_xxx[i];

        t_xz_xy[i] = -t_z_xy[i] * ab_x[i] + t_z_xxy[i];

        t_xz_xz[i] = -t_z_xz[i] * ab_x[i] + t_z_xxz[i];

        t_xz_yy[i] = -t_z_yy[i] * ab_x[i] + t_z_xyy[i];

        t_xz_yz[i] = -t_z_yz[i] * ab_x[i] + t_z_xyz[i];

        t_xz_zz[i] = -t_z_zz[i] * ab_x[i] + t_z_xzz[i];

        t_yy_xx[i] = -t_y_xx[i] * ab_y[i] + t_y_xxy[i];

        t_yy_xy[i] = -t_y_xy[i] * ab_y[i] + t_y_xyy[i];

        t_yy_xz[i] = -t_y_xz[i] * ab_y[i] + t_y_xyz[i];

        t_yy_yy[i] = -t_y_yy[i] * ab_y[i] + t_y_yyy[i];

        t_yy_yz[i] = -t_y_yz[i] * ab_y[i] + t_y_yyz[i];

        t_yy_zz[i] = -t_y_zz[i] * ab_y[i] + t_y_yzz[i];

        t_yz_xx[i] = -t_z_xx[i] * ab_y[i] + t_z_xxy[i];

        t_yz_xy[i] = -t_z_xy[i] * ab_y[i] + t_z_xyy[i];

        t_yz_xz[i] = -t_z_xz[i] * ab_y[i] + t_z_xyz[i];

        t_yz_yy[i] = -t_z_yy[i] * ab_y[i] + t_z_yyy[i];

        t_yz_yz[i] = -t_z_yz[i] * ab_y[i] + t_z_yyz[i];

        t_yz_zz[i] = -t_z_zz[i] * ab_y[i] + t_z_yzz[i];

        t_zz_xx[i] = -t_z_xx[i] * ab_z[i] + t_z_xxz[i];

        t_zz_xy[i] = -t_z_xy[i] * ab_z[i] + t_z_xyz[i];

        t_zz_xz[i] = -t_z_xz[i] * ab_z[i] + t_z_xzz[i];

        t_zz_yy[i] = -t_z_yy[i] * ab_z[i] + t_z_yyz[i];

        t_zz_yz[i] = -t_z_yz[i] * ab_z[i] + t_z_yzz[i];

        t_zz_zz[i] = -t_z_zz[i] * ab_z[i] + t_z_zzz[i];
    }
}

} // t2chrr namespace


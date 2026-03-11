#include "T2CHrrABRecPD.hpp"

namespace t2chrr { // t2chrr namespace

auto
comp_hrr_pd(CSimdArray<double>& cbuffer, 
            const size_t idx_pd,
            const size_t idx_sd,
            const size_t idx_sf,
            const CSimdArray<double>& factors) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    // Set up R(AB) distances

    auto ab_x = factors.data(3);

    auto ab_y = factors.data(4);

    auto ab_z = factors.data(5);

    // Set up components of auxiliary buffer : SD

    auto t_0_xx = cbuffer.data(idx_sd);

    auto t_0_xy = cbuffer.data(idx_sd + 1);

    auto t_0_xz = cbuffer.data(idx_sd + 2);

    auto t_0_yy = cbuffer.data(idx_sd + 3);

    auto t_0_yz = cbuffer.data(idx_sd + 4);

    auto t_0_zz = cbuffer.data(idx_sd + 5);

    // Set up components of auxiliary buffer : SF

    auto t_0_xxx = cbuffer.data(idx_sf);

    auto t_0_xxy = cbuffer.data(idx_sf + 1);

    auto t_0_xxz = cbuffer.data(idx_sf + 2);

    auto t_0_xyy = cbuffer.data(idx_sf + 3);

    auto t_0_xyz = cbuffer.data(idx_sf + 4);

    auto t_0_xzz = cbuffer.data(idx_sf + 5);

    auto t_0_yyy = cbuffer.data(idx_sf + 6);

    auto t_0_yyz = cbuffer.data(idx_sf + 7);

    auto t_0_yzz = cbuffer.data(idx_sf + 8);

    auto t_0_zzz = cbuffer.data(idx_sf + 9);

    // Set up components of targeted buffer : PD

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

    #pragma omp simd aligned(ab_x, ab_y, ab_z, t_0_xx, t_0_xxx, t_0_xxy, t_0_xxz, t_0_xy, t_0_xyy, t_0_xyz, t_0_xz, t_0_xzz, t_0_yy, t_0_yyy, t_0_yyz, t_0_yz, t_0_yzz, t_0_zz, t_0_zzz, t_x_xx, t_x_xy, t_x_xz, t_x_yy, t_x_yz, t_x_zz, t_y_xx, t_y_xy, t_y_xz, t_y_yy, t_y_yz, t_y_zz, t_z_xx, t_z_xy, t_z_xz, t_z_yy, t_z_yz, t_z_zz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        t_x_xx[i] = -t_0_xx[i] * ab_x[i] + t_0_xxx[i];

        t_x_xy[i] = -t_0_xy[i] * ab_x[i] + t_0_xxy[i];

        t_x_xz[i] = -t_0_xz[i] * ab_x[i] + t_0_xxz[i];

        t_x_yy[i] = -t_0_yy[i] * ab_x[i] + t_0_xyy[i];

        t_x_yz[i] = -t_0_yz[i] * ab_x[i] + t_0_xyz[i];

        t_x_zz[i] = -t_0_zz[i] * ab_x[i] + t_0_xzz[i];

        t_y_xx[i] = -t_0_xx[i] * ab_y[i] + t_0_xxy[i];

        t_y_xy[i] = -t_0_xy[i] * ab_y[i] + t_0_xyy[i];

        t_y_xz[i] = -t_0_xz[i] * ab_y[i] + t_0_xyz[i];

        t_y_yy[i] = -t_0_yy[i] * ab_y[i] + t_0_yyy[i];

        t_y_yz[i] = -t_0_yz[i] * ab_y[i] + t_0_yyz[i];

        t_y_zz[i] = -t_0_zz[i] * ab_y[i] + t_0_yzz[i];

        t_z_xx[i] = -t_0_xx[i] * ab_z[i] + t_0_xxz[i];

        t_z_xy[i] = -t_0_xy[i] * ab_z[i] + t_0_xyz[i];

        t_z_xz[i] = -t_0_xz[i] * ab_z[i] + t_0_xzz[i];

        t_z_yy[i] = -t_0_yy[i] * ab_z[i] + t_0_yyz[i];

        t_z_yz[i] = -t_0_yz[i] * ab_z[i] + t_0_yzz[i];

        t_z_zz[i] = -t_0_zz[i] * ab_z[i] + t_0_zzz[i];
    }
}

} // t2chrr namespace


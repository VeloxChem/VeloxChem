#include "T2CHrrABRecPF.hpp"

namespace t2chrr { // t2chrr namespace

auto
comp_hrr_pf(CSimdArray<double>& cbuffer, 
            const size_t idx_pf,
            const size_t idx_sf,
            const size_t idx_sg,
            const CSimdArray<double>& factors) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    // Set up R(AB) distances

    auto ab_x = factors.data(3);

    auto ab_y = factors.data(4);

    auto ab_z = factors.data(5);

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

    // Set up components of auxiliary buffer : SG

    auto t_0_xxxx = cbuffer.data(idx_sg);

    auto t_0_xxxy = cbuffer.data(idx_sg + 1);

    auto t_0_xxxz = cbuffer.data(idx_sg + 2);

    auto t_0_xxyy = cbuffer.data(idx_sg + 3);

    auto t_0_xxyz = cbuffer.data(idx_sg + 4);

    auto t_0_xxzz = cbuffer.data(idx_sg + 5);

    auto t_0_xyyy = cbuffer.data(idx_sg + 6);

    auto t_0_xyyz = cbuffer.data(idx_sg + 7);

    auto t_0_xyzz = cbuffer.data(idx_sg + 8);

    auto t_0_xzzz = cbuffer.data(idx_sg + 9);

    auto t_0_yyyy = cbuffer.data(idx_sg + 10);

    auto t_0_yyyz = cbuffer.data(idx_sg + 11);

    auto t_0_yyzz = cbuffer.data(idx_sg + 12);

    auto t_0_yzzz = cbuffer.data(idx_sg + 13);

    auto t_0_zzzz = cbuffer.data(idx_sg + 14);

    // Set up components of targeted buffer : PF

    auto t_x_xxx = cbuffer.data(idx_pf);

    auto t_x_xxy = cbuffer.data(idx_pf + 1);

    auto t_x_xxz = cbuffer.data(idx_pf + 2);

    auto t_x_xyy = cbuffer.data(idx_pf + 3);

    auto t_x_xyz = cbuffer.data(idx_pf + 4);

    auto t_x_xzz = cbuffer.data(idx_pf + 5);

    auto t_x_yyy = cbuffer.data(idx_pf + 6);

    auto t_x_yyz = cbuffer.data(idx_pf + 7);

    auto t_x_yzz = cbuffer.data(idx_pf + 8);

    auto t_x_zzz = cbuffer.data(idx_pf + 9);

    auto t_y_xxx = cbuffer.data(idx_pf + 10);

    auto t_y_xxy = cbuffer.data(idx_pf + 11);

    auto t_y_xxz = cbuffer.data(idx_pf + 12);

    auto t_y_xyy = cbuffer.data(idx_pf + 13);

    auto t_y_xyz = cbuffer.data(idx_pf + 14);

    auto t_y_xzz = cbuffer.data(idx_pf + 15);

    auto t_y_yyy = cbuffer.data(idx_pf + 16);

    auto t_y_yyz = cbuffer.data(idx_pf + 17);

    auto t_y_yzz = cbuffer.data(idx_pf + 18);

    auto t_y_zzz = cbuffer.data(idx_pf + 19);

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

    #pragma omp simd aligned(ab_x, ab_y, ab_z, t_0_xxx, t_0_xxxx, t_0_xxxy, t_0_xxxz, t_0_xxy, t_0_xxyy, t_0_xxyz, t_0_xxz, t_0_xxzz, t_0_xyy, t_0_xyyy, t_0_xyyz, t_0_xyz, t_0_xyzz, t_0_xzz, t_0_xzzz, t_0_yyy, t_0_yyyy, t_0_yyyz, t_0_yyz, t_0_yyzz, t_0_yzz, t_0_yzzz, t_0_zzz, t_0_zzzz, t_x_xxx, t_x_xxy, t_x_xxz, t_x_xyy, t_x_xyz, t_x_xzz, t_x_yyy, t_x_yyz, t_x_yzz, t_x_zzz, t_y_xxx, t_y_xxy, t_y_xxz, t_y_xyy, t_y_xyz, t_y_xzz, t_y_yyy, t_y_yyz, t_y_yzz, t_y_zzz, t_z_xxx, t_z_xxy, t_z_xxz, t_z_xyy, t_z_xyz, t_z_xzz, t_z_yyy, t_z_yyz, t_z_yzz, t_z_zzz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        t_x_xxx[i] = -t_0_xxx[i] * ab_x[i] + t_0_xxxx[i];

        t_x_xxy[i] = -t_0_xxy[i] * ab_x[i] + t_0_xxxy[i];

        t_x_xxz[i] = -t_0_xxz[i] * ab_x[i] + t_0_xxxz[i];

        t_x_xyy[i] = -t_0_xyy[i] * ab_x[i] + t_0_xxyy[i];

        t_x_xyz[i] = -t_0_xyz[i] * ab_x[i] + t_0_xxyz[i];

        t_x_xzz[i] = -t_0_xzz[i] * ab_x[i] + t_0_xxzz[i];

        t_x_yyy[i] = -t_0_yyy[i] * ab_x[i] + t_0_xyyy[i];

        t_x_yyz[i] = -t_0_yyz[i] * ab_x[i] + t_0_xyyz[i];

        t_x_yzz[i] = -t_0_yzz[i] * ab_x[i] + t_0_xyzz[i];

        t_x_zzz[i] = -t_0_zzz[i] * ab_x[i] + t_0_xzzz[i];

        t_y_xxx[i] = -t_0_xxx[i] * ab_y[i] + t_0_xxxy[i];

        t_y_xxy[i] = -t_0_xxy[i] * ab_y[i] + t_0_xxyy[i];

        t_y_xxz[i] = -t_0_xxz[i] * ab_y[i] + t_0_xxyz[i];

        t_y_xyy[i] = -t_0_xyy[i] * ab_y[i] + t_0_xyyy[i];

        t_y_xyz[i] = -t_0_xyz[i] * ab_y[i] + t_0_xyyz[i];

        t_y_xzz[i] = -t_0_xzz[i] * ab_y[i] + t_0_xyzz[i];

        t_y_yyy[i] = -t_0_yyy[i] * ab_y[i] + t_0_yyyy[i];

        t_y_yyz[i] = -t_0_yyz[i] * ab_y[i] + t_0_yyyz[i];

        t_y_yzz[i] = -t_0_yzz[i] * ab_y[i] + t_0_yyzz[i];

        t_y_zzz[i] = -t_0_zzz[i] * ab_y[i] + t_0_yzzz[i];

        t_z_xxx[i] = -t_0_xxx[i] * ab_z[i] + t_0_xxxz[i];

        t_z_xxy[i] = -t_0_xxy[i] * ab_z[i] + t_0_xxyz[i];

        t_z_xxz[i] = -t_0_xxz[i] * ab_z[i] + t_0_xxzz[i];

        t_z_xyy[i] = -t_0_xyy[i] * ab_z[i] + t_0_xyyz[i];

        t_z_xyz[i] = -t_0_xyz[i] * ab_z[i] + t_0_xyzz[i];

        t_z_xzz[i] = -t_0_xzz[i] * ab_z[i] + t_0_xzzz[i];

        t_z_yyy[i] = -t_0_yyy[i] * ab_z[i] + t_0_yyyz[i];

        t_z_yyz[i] = -t_0_yyz[i] * ab_z[i] + t_0_yyzz[i];

        t_z_yzz[i] = -t_0_yzz[i] * ab_z[i] + t_0_yzzz[i];

        t_z_zzz[i] = -t_0_zzz[i] * ab_z[i] + t_0_zzzz[i];
    }
}

} // t2chrr namespace


#include "T2CHrrABRecDF.hpp"

namespace t2chrr { // t2chrr namespace

auto
comp_hrr_df(CSimdArray<double>& cbuffer, 
            const size_t idx_df,
            const size_t idx_pf,
            const size_t idx_pg,
            const CSimdArray<double>& factors) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    // Set up R(AB) distances

    auto ab_x = factors.data(3);

    auto ab_y = factors.data(4);

    auto ab_z = factors.data(5);

    // Set up components of auxiliary buffer : PF

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

    // Set up components of auxiliary buffer : PG

    auto t_x_xxxx = cbuffer.data(idx_pg);

    auto t_x_xxxy = cbuffer.data(idx_pg + 1);

    auto t_x_xxxz = cbuffer.data(idx_pg + 2);

    auto t_x_xxyy = cbuffer.data(idx_pg + 3);

    auto t_x_xxyz = cbuffer.data(idx_pg + 4);

    auto t_x_xxzz = cbuffer.data(idx_pg + 5);

    auto t_x_xyyy = cbuffer.data(idx_pg + 6);

    auto t_x_xyyz = cbuffer.data(idx_pg + 7);

    auto t_x_xyzz = cbuffer.data(idx_pg + 8);

    auto t_x_xzzz = cbuffer.data(idx_pg + 9);

    auto t_y_xxxx = cbuffer.data(idx_pg + 15);

    auto t_y_xxxy = cbuffer.data(idx_pg + 16);

    auto t_y_xxxz = cbuffer.data(idx_pg + 17);

    auto t_y_xxyy = cbuffer.data(idx_pg + 18);

    auto t_y_xxyz = cbuffer.data(idx_pg + 19);

    auto t_y_xxzz = cbuffer.data(idx_pg + 20);

    auto t_y_xyyy = cbuffer.data(idx_pg + 21);

    auto t_y_xyyz = cbuffer.data(idx_pg + 22);

    auto t_y_xyzz = cbuffer.data(idx_pg + 23);

    auto t_y_xzzz = cbuffer.data(idx_pg + 24);

    auto t_y_yyyy = cbuffer.data(idx_pg + 25);

    auto t_y_yyyz = cbuffer.data(idx_pg + 26);

    auto t_y_yyzz = cbuffer.data(idx_pg + 27);

    auto t_y_yzzz = cbuffer.data(idx_pg + 28);

    auto t_z_xxxx = cbuffer.data(idx_pg + 30);

    auto t_z_xxxy = cbuffer.data(idx_pg + 31);

    auto t_z_xxxz = cbuffer.data(idx_pg + 32);

    auto t_z_xxyy = cbuffer.data(idx_pg + 33);

    auto t_z_xxyz = cbuffer.data(idx_pg + 34);

    auto t_z_xxzz = cbuffer.data(idx_pg + 35);

    auto t_z_xyyy = cbuffer.data(idx_pg + 36);

    auto t_z_xyyz = cbuffer.data(idx_pg + 37);

    auto t_z_xyzz = cbuffer.data(idx_pg + 38);

    auto t_z_xzzz = cbuffer.data(idx_pg + 39);

    auto t_z_yyyy = cbuffer.data(idx_pg + 40);

    auto t_z_yyyz = cbuffer.data(idx_pg + 41);

    auto t_z_yyzz = cbuffer.data(idx_pg + 42);

    auto t_z_yzzz = cbuffer.data(idx_pg + 43);

    auto t_z_zzzz = cbuffer.data(idx_pg + 44);

    // Set up components of targeted buffer : DF

    auto t_xx_xxx = cbuffer.data(idx_df);

    auto t_xx_xxy = cbuffer.data(idx_df + 1);

    auto t_xx_xxz = cbuffer.data(idx_df + 2);

    auto t_xx_xyy = cbuffer.data(idx_df + 3);

    auto t_xx_xyz = cbuffer.data(idx_df + 4);

    auto t_xx_xzz = cbuffer.data(idx_df + 5);

    auto t_xx_yyy = cbuffer.data(idx_df + 6);

    auto t_xx_yyz = cbuffer.data(idx_df + 7);

    auto t_xx_yzz = cbuffer.data(idx_df + 8);

    auto t_xx_zzz = cbuffer.data(idx_df + 9);

    auto t_xy_xxx = cbuffer.data(idx_df + 10);

    auto t_xy_xxy = cbuffer.data(idx_df + 11);

    auto t_xy_xxz = cbuffer.data(idx_df + 12);

    auto t_xy_xyy = cbuffer.data(idx_df + 13);

    auto t_xy_xyz = cbuffer.data(idx_df + 14);

    auto t_xy_xzz = cbuffer.data(idx_df + 15);

    auto t_xy_yyy = cbuffer.data(idx_df + 16);

    auto t_xy_yyz = cbuffer.data(idx_df + 17);

    auto t_xy_yzz = cbuffer.data(idx_df + 18);

    auto t_xy_zzz = cbuffer.data(idx_df + 19);

    auto t_xz_xxx = cbuffer.data(idx_df + 20);

    auto t_xz_xxy = cbuffer.data(idx_df + 21);

    auto t_xz_xxz = cbuffer.data(idx_df + 22);

    auto t_xz_xyy = cbuffer.data(idx_df + 23);

    auto t_xz_xyz = cbuffer.data(idx_df + 24);

    auto t_xz_xzz = cbuffer.data(idx_df + 25);

    auto t_xz_yyy = cbuffer.data(idx_df + 26);

    auto t_xz_yyz = cbuffer.data(idx_df + 27);

    auto t_xz_yzz = cbuffer.data(idx_df + 28);

    auto t_xz_zzz = cbuffer.data(idx_df + 29);

    auto t_yy_xxx = cbuffer.data(idx_df + 30);

    auto t_yy_xxy = cbuffer.data(idx_df + 31);

    auto t_yy_xxz = cbuffer.data(idx_df + 32);

    auto t_yy_xyy = cbuffer.data(idx_df + 33);

    auto t_yy_xyz = cbuffer.data(idx_df + 34);

    auto t_yy_xzz = cbuffer.data(idx_df + 35);

    auto t_yy_yyy = cbuffer.data(idx_df + 36);

    auto t_yy_yyz = cbuffer.data(idx_df + 37);

    auto t_yy_yzz = cbuffer.data(idx_df + 38);

    auto t_yy_zzz = cbuffer.data(idx_df + 39);

    auto t_yz_xxx = cbuffer.data(idx_df + 40);

    auto t_yz_xxy = cbuffer.data(idx_df + 41);

    auto t_yz_xxz = cbuffer.data(idx_df + 42);

    auto t_yz_xyy = cbuffer.data(idx_df + 43);

    auto t_yz_xyz = cbuffer.data(idx_df + 44);

    auto t_yz_xzz = cbuffer.data(idx_df + 45);

    auto t_yz_yyy = cbuffer.data(idx_df + 46);

    auto t_yz_yyz = cbuffer.data(idx_df + 47);

    auto t_yz_yzz = cbuffer.data(idx_df + 48);

    auto t_yz_zzz = cbuffer.data(idx_df + 49);

    auto t_zz_xxx = cbuffer.data(idx_df + 50);

    auto t_zz_xxy = cbuffer.data(idx_df + 51);

    auto t_zz_xxz = cbuffer.data(idx_df + 52);

    auto t_zz_xyy = cbuffer.data(idx_df + 53);

    auto t_zz_xyz = cbuffer.data(idx_df + 54);

    auto t_zz_xzz = cbuffer.data(idx_df + 55);

    auto t_zz_yyy = cbuffer.data(idx_df + 56);

    auto t_zz_yyz = cbuffer.data(idx_df + 57);

    auto t_zz_yzz = cbuffer.data(idx_df + 58);

    auto t_zz_zzz = cbuffer.data(idx_df + 59);

    #pragma omp simd aligned(ab_x, ab_y, ab_z, t_x_xxx, t_x_xxxx, t_x_xxxy, t_x_xxxz, t_x_xxy, t_x_xxyy, t_x_xxyz, t_x_xxz, t_x_xxzz, t_x_xyy, t_x_xyyy, t_x_xyyz, t_x_xyz, t_x_xyzz, t_x_xzz, t_x_xzzz, t_x_yyy, t_x_yyz, t_x_yzz, t_x_zzz, t_xx_xxx, t_xx_xxy, t_xx_xxz, t_xx_xyy, t_xx_xyz, t_xx_xzz, t_xx_yyy, t_xx_yyz, t_xx_yzz, t_xx_zzz, t_xy_xxx, t_xy_xxy, t_xy_xxz, t_xy_xyy, t_xy_xyz, t_xy_xzz, t_xy_yyy, t_xy_yyz, t_xy_yzz, t_xy_zzz, t_xz_xxx, t_xz_xxy, t_xz_xxz, t_xz_xyy, t_xz_xyz, t_xz_xzz, t_xz_yyy, t_xz_yyz, t_xz_yzz, t_xz_zzz, t_y_xxx, t_y_xxxx, t_y_xxxy, t_y_xxxz, t_y_xxy, t_y_xxyy, t_y_xxyz, t_y_xxz, t_y_xxzz, t_y_xyy, t_y_xyyy, t_y_xyyz, t_y_xyz, t_y_xyzz, t_y_xzz, t_y_xzzz, t_y_yyy, t_y_yyyy, t_y_yyyz, t_y_yyz, t_y_yyzz, t_y_yzz, t_y_yzzz, t_y_zzz, t_yy_xxx, t_yy_xxy, t_yy_xxz, t_yy_xyy, t_yy_xyz, t_yy_xzz, t_yy_yyy, t_yy_yyz, t_yy_yzz, t_yy_zzz, t_yz_xxx, t_yz_xxy, t_yz_xxz, t_yz_xyy, t_yz_xyz, t_yz_xzz, t_yz_yyy, t_yz_yyz, t_yz_yzz, t_yz_zzz, t_z_xxx, t_z_xxxx, t_z_xxxy, t_z_xxxz, t_z_xxy, t_z_xxyy, t_z_xxyz, t_z_xxz, t_z_xxzz, t_z_xyy, t_z_xyyy, t_z_xyyz, t_z_xyz, t_z_xyzz, t_z_xzz, t_z_xzzz, t_z_yyy, t_z_yyyy, t_z_yyyz, t_z_yyz, t_z_yyzz, t_z_yzz, t_z_yzzz, t_z_zzz, t_z_zzzz, t_zz_xxx, t_zz_xxy, t_zz_xxz, t_zz_xyy, t_zz_xyz, t_zz_xzz, t_zz_yyy, t_zz_yyz, t_zz_yzz, t_zz_zzz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        t_xx_xxx[i] = -t_x_xxx[i] * ab_x[i] + t_x_xxxx[i];

        t_xx_xxy[i] = -t_x_xxy[i] * ab_x[i] + t_x_xxxy[i];

        t_xx_xxz[i] = -t_x_xxz[i] * ab_x[i] + t_x_xxxz[i];

        t_xx_xyy[i] = -t_x_xyy[i] * ab_x[i] + t_x_xxyy[i];

        t_xx_xyz[i] = -t_x_xyz[i] * ab_x[i] + t_x_xxyz[i];

        t_xx_xzz[i] = -t_x_xzz[i] * ab_x[i] + t_x_xxzz[i];

        t_xx_yyy[i] = -t_x_yyy[i] * ab_x[i] + t_x_xyyy[i];

        t_xx_yyz[i] = -t_x_yyz[i] * ab_x[i] + t_x_xyyz[i];

        t_xx_yzz[i] = -t_x_yzz[i] * ab_x[i] + t_x_xyzz[i];

        t_xx_zzz[i] = -t_x_zzz[i] * ab_x[i] + t_x_xzzz[i];

        t_xy_xxx[i] = -t_y_xxx[i] * ab_x[i] + t_y_xxxx[i];

        t_xy_xxy[i] = -t_y_xxy[i] * ab_x[i] + t_y_xxxy[i];

        t_xy_xxz[i] = -t_y_xxz[i] * ab_x[i] + t_y_xxxz[i];

        t_xy_xyy[i] = -t_y_xyy[i] * ab_x[i] + t_y_xxyy[i];

        t_xy_xyz[i] = -t_y_xyz[i] * ab_x[i] + t_y_xxyz[i];

        t_xy_xzz[i] = -t_y_xzz[i] * ab_x[i] + t_y_xxzz[i];

        t_xy_yyy[i] = -t_y_yyy[i] * ab_x[i] + t_y_xyyy[i];

        t_xy_yyz[i] = -t_y_yyz[i] * ab_x[i] + t_y_xyyz[i];

        t_xy_yzz[i] = -t_y_yzz[i] * ab_x[i] + t_y_xyzz[i];

        t_xy_zzz[i] = -t_y_zzz[i] * ab_x[i] + t_y_xzzz[i];

        t_xz_xxx[i] = -t_z_xxx[i] * ab_x[i] + t_z_xxxx[i];

        t_xz_xxy[i] = -t_z_xxy[i] * ab_x[i] + t_z_xxxy[i];

        t_xz_xxz[i] = -t_z_xxz[i] * ab_x[i] + t_z_xxxz[i];

        t_xz_xyy[i] = -t_z_xyy[i] * ab_x[i] + t_z_xxyy[i];

        t_xz_xyz[i] = -t_z_xyz[i] * ab_x[i] + t_z_xxyz[i];

        t_xz_xzz[i] = -t_z_xzz[i] * ab_x[i] + t_z_xxzz[i];

        t_xz_yyy[i] = -t_z_yyy[i] * ab_x[i] + t_z_xyyy[i];

        t_xz_yyz[i] = -t_z_yyz[i] * ab_x[i] + t_z_xyyz[i];

        t_xz_yzz[i] = -t_z_yzz[i] * ab_x[i] + t_z_xyzz[i];

        t_xz_zzz[i] = -t_z_zzz[i] * ab_x[i] + t_z_xzzz[i];

        t_yy_xxx[i] = -t_y_xxx[i] * ab_y[i] + t_y_xxxy[i];

        t_yy_xxy[i] = -t_y_xxy[i] * ab_y[i] + t_y_xxyy[i];

        t_yy_xxz[i] = -t_y_xxz[i] * ab_y[i] + t_y_xxyz[i];

        t_yy_xyy[i] = -t_y_xyy[i] * ab_y[i] + t_y_xyyy[i];

        t_yy_xyz[i] = -t_y_xyz[i] * ab_y[i] + t_y_xyyz[i];

        t_yy_xzz[i] = -t_y_xzz[i] * ab_y[i] + t_y_xyzz[i];

        t_yy_yyy[i] = -t_y_yyy[i] * ab_y[i] + t_y_yyyy[i];

        t_yy_yyz[i] = -t_y_yyz[i] * ab_y[i] + t_y_yyyz[i];

        t_yy_yzz[i] = -t_y_yzz[i] * ab_y[i] + t_y_yyzz[i];

        t_yy_zzz[i] = -t_y_zzz[i] * ab_y[i] + t_y_yzzz[i];

        t_yz_xxx[i] = -t_z_xxx[i] * ab_y[i] + t_z_xxxy[i];

        t_yz_xxy[i] = -t_z_xxy[i] * ab_y[i] + t_z_xxyy[i];

        t_yz_xxz[i] = -t_z_xxz[i] * ab_y[i] + t_z_xxyz[i];

        t_yz_xyy[i] = -t_z_xyy[i] * ab_y[i] + t_z_xyyy[i];

        t_yz_xyz[i] = -t_z_xyz[i] * ab_y[i] + t_z_xyyz[i];

        t_yz_xzz[i] = -t_z_xzz[i] * ab_y[i] + t_z_xyzz[i];

        t_yz_yyy[i] = -t_z_yyy[i] * ab_y[i] + t_z_yyyy[i];

        t_yz_yyz[i] = -t_z_yyz[i] * ab_y[i] + t_z_yyyz[i];

        t_yz_yzz[i] = -t_z_yzz[i] * ab_y[i] + t_z_yyzz[i];

        t_yz_zzz[i] = -t_z_zzz[i] * ab_y[i] + t_z_yzzz[i];

        t_zz_xxx[i] = -t_z_xxx[i] * ab_z[i] + t_z_xxxz[i];

        t_zz_xxy[i] = -t_z_xxy[i] * ab_z[i] + t_z_xxyz[i];

        t_zz_xxz[i] = -t_z_xxz[i] * ab_z[i] + t_z_xxzz[i];

        t_zz_xyy[i] = -t_z_xyy[i] * ab_z[i] + t_z_xyyz[i];

        t_zz_xyz[i] = -t_z_xyz[i] * ab_z[i] + t_z_xyzz[i];

        t_zz_xzz[i] = -t_z_xzz[i] * ab_z[i] + t_z_xzzz[i];

        t_zz_yyy[i] = -t_z_yyy[i] * ab_z[i] + t_z_yyyz[i];

        t_zz_yyz[i] = -t_z_yyz[i] * ab_z[i] + t_z_yyzz[i];

        t_zz_yzz[i] = -t_z_yzz[i] * ab_z[i] + t_z_yzzz[i];

        t_zz_zzz[i] = -t_z_zzz[i] * ab_z[i] + t_z_zzzz[i];
    }
}

} // t2chrr namespace


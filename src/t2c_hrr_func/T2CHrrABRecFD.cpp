#include "T2CHrrABRecFD.hpp"

namespace t2chrr { // t2chrr namespace

auto
comp_hrr_fd(CSimdArray<double>& cbuffer, 
            const size_t idx_fd,
            const size_t idx_fp,
            const size_t idx_gp,
            const CSimdArray<double>& factors) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    // Set up R(AB) distances

    auto ab_x = factors.data(3);

    auto ab_y = factors.data(4);

    auto ab_z = factors.data(5);

    // Set up components of auxiliary buffer : FP

    auto t_xxx_x = cbuffer.data(idx_fp);

    auto t_xxx_y = cbuffer.data(idx_fp + 1);

    auto t_xxx_z = cbuffer.data(idx_fp + 2);

    auto t_xxy_x = cbuffer.data(idx_fp + 3);

    auto t_xxy_y = cbuffer.data(idx_fp + 4);

    auto t_xxy_z = cbuffer.data(idx_fp + 5);

    auto t_xxz_x = cbuffer.data(idx_fp + 6);

    auto t_xxz_y = cbuffer.data(idx_fp + 7);

    auto t_xxz_z = cbuffer.data(idx_fp + 8);

    auto t_xyy_x = cbuffer.data(idx_fp + 9);

    auto t_xyy_y = cbuffer.data(idx_fp + 10);

    auto t_xyy_z = cbuffer.data(idx_fp + 11);

    auto t_xyz_x = cbuffer.data(idx_fp + 12);

    auto t_xyz_y = cbuffer.data(idx_fp + 13);

    auto t_xyz_z = cbuffer.data(idx_fp + 14);

    auto t_xzz_x = cbuffer.data(idx_fp + 15);

    auto t_xzz_y = cbuffer.data(idx_fp + 16);

    auto t_xzz_z = cbuffer.data(idx_fp + 17);

    auto t_yyy_x = cbuffer.data(idx_fp + 18);

    auto t_yyy_y = cbuffer.data(idx_fp + 19);

    auto t_yyy_z = cbuffer.data(idx_fp + 20);

    auto t_yyz_x = cbuffer.data(idx_fp + 21);

    auto t_yyz_y = cbuffer.data(idx_fp + 22);

    auto t_yyz_z = cbuffer.data(idx_fp + 23);

    auto t_yzz_x = cbuffer.data(idx_fp + 24);

    auto t_yzz_y = cbuffer.data(idx_fp + 25);

    auto t_yzz_z = cbuffer.data(idx_fp + 26);

    auto t_zzz_x = cbuffer.data(idx_fp + 27);

    auto t_zzz_y = cbuffer.data(idx_fp + 28);

    auto t_zzz_z = cbuffer.data(idx_fp + 29);

    // Set up components of auxiliary buffer : GP

    auto t_xxxx_x = cbuffer.data(idx_gp);

    auto t_xxxx_y = cbuffer.data(idx_gp + 1);

    auto t_xxxx_z = cbuffer.data(idx_gp + 2);

    auto t_xxxy_x = cbuffer.data(idx_gp + 3);

    auto t_xxxy_y = cbuffer.data(idx_gp + 4);

    auto t_xxxy_z = cbuffer.data(idx_gp + 5);

    auto t_xxxz_x = cbuffer.data(idx_gp + 6);

    auto t_xxxz_y = cbuffer.data(idx_gp + 7);

    auto t_xxxz_z = cbuffer.data(idx_gp + 8);

    auto t_xxyy_x = cbuffer.data(idx_gp + 9);

    auto t_xxyy_y = cbuffer.data(idx_gp + 10);

    auto t_xxyy_z = cbuffer.data(idx_gp + 11);

    auto t_xxyz_x = cbuffer.data(idx_gp + 12);

    auto t_xxyz_y = cbuffer.data(idx_gp + 13);

    auto t_xxyz_z = cbuffer.data(idx_gp + 14);

    auto t_xxzz_x = cbuffer.data(idx_gp + 15);

    auto t_xxzz_y = cbuffer.data(idx_gp + 16);

    auto t_xxzz_z = cbuffer.data(idx_gp + 17);

    auto t_xyyy_x = cbuffer.data(idx_gp + 18);

    auto t_xyyy_y = cbuffer.data(idx_gp + 19);

    auto t_xyyy_z = cbuffer.data(idx_gp + 20);

    auto t_xyyz_x = cbuffer.data(idx_gp + 21);

    auto t_xyyz_y = cbuffer.data(idx_gp + 22);

    auto t_xyyz_z = cbuffer.data(idx_gp + 23);

    auto t_xyzz_x = cbuffer.data(idx_gp + 24);

    auto t_xyzz_y = cbuffer.data(idx_gp + 25);

    auto t_xyzz_z = cbuffer.data(idx_gp + 26);

    auto t_xzzz_x = cbuffer.data(idx_gp + 27);

    auto t_xzzz_y = cbuffer.data(idx_gp + 28);

    auto t_xzzz_z = cbuffer.data(idx_gp + 29);

    auto t_yyyy_y = cbuffer.data(idx_gp + 31);

    auto t_yyyy_z = cbuffer.data(idx_gp + 32);

    auto t_yyyz_y = cbuffer.data(idx_gp + 34);

    auto t_yyyz_z = cbuffer.data(idx_gp + 35);

    auto t_yyzz_y = cbuffer.data(idx_gp + 37);

    auto t_yyzz_z = cbuffer.data(idx_gp + 38);

    auto t_yzzz_y = cbuffer.data(idx_gp + 40);

    auto t_yzzz_z = cbuffer.data(idx_gp + 41);

    auto t_zzzz_z = cbuffer.data(idx_gp + 44);

    // Set up components of targeted buffer : FD

    auto t_xxx_xx = cbuffer.data(idx_fd);

    auto t_xxx_xy = cbuffer.data(idx_fd + 1);

    auto t_xxx_xz = cbuffer.data(idx_fd + 2);

    auto t_xxx_yy = cbuffer.data(idx_fd + 3);

    auto t_xxx_yz = cbuffer.data(idx_fd + 4);

    auto t_xxx_zz = cbuffer.data(idx_fd + 5);

    auto t_xxy_xx = cbuffer.data(idx_fd + 6);

    auto t_xxy_xy = cbuffer.data(idx_fd + 7);

    auto t_xxy_xz = cbuffer.data(idx_fd + 8);

    auto t_xxy_yy = cbuffer.data(idx_fd + 9);

    auto t_xxy_yz = cbuffer.data(idx_fd + 10);

    auto t_xxy_zz = cbuffer.data(idx_fd + 11);

    auto t_xxz_xx = cbuffer.data(idx_fd + 12);

    auto t_xxz_xy = cbuffer.data(idx_fd + 13);

    auto t_xxz_xz = cbuffer.data(idx_fd + 14);

    auto t_xxz_yy = cbuffer.data(idx_fd + 15);

    auto t_xxz_yz = cbuffer.data(idx_fd + 16);

    auto t_xxz_zz = cbuffer.data(idx_fd + 17);

    auto t_xyy_xx = cbuffer.data(idx_fd + 18);

    auto t_xyy_xy = cbuffer.data(idx_fd + 19);

    auto t_xyy_xz = cbuffer.data(idx_fd + 20);

    auto t_xyy_yy = cbuffer.data(idx_fd + 21);

    auto t_xyy_yz = cbuffer.data(idx_fd + 22);

    auto t_xyy_zz = cbuffer.data(idx_fd + 23);

    auto t_xyz_xx = cbuffer.data(idx_fd + 24);

    auto t_xyz_xy = cbuffer.data(idx_fd + 25);

    auto t_xyz_xz = cbuffer.data(idx_fd + 26);

    auto t_xyz_yy = cbuffer.data(idx_fd + 27);

    auto t_xyz_yz = cbuffer.data(idx_fd + 28);

    auto t_xyz_zz = cbuffer.data(idx_fd + 29);

    auto t_xzz_xx = cbuffer.data(idx_fd + 30);

    auto t_xzz_xy = cbuffer.data(idx_fd + 31);

    auto t_xzz_xz = cbuffer.data(idx_fd + 32);

    auto t_xzz_yy = cbuffer.data(idx_fd + 33);

    auto t_xzz_yz = cbuffer.data(idx_fd + 34);

    auto t_xzz_zz = cbuffer.data(idx_fd + 35);

    auto t_yyy_xx = cbuffer.data(idx_fd + 36);

    auto t_yyy_xy = cbuffer.data(idx_fd + 37);

    auto t_yyy_xz = cbuffer.data(idx_fd + 38);

    auto t_yyy_yy = cbuffer.data(idx_fd + 39);

    auto t_yyy_yz = cbuffer.data(idx_fd + 40);

    auto t_yyy_zz = cbuffer.data(idx_fd + 41);

    auto t_yyz_xx = cbuffer.data(idx_fd + 42);

    auto t_yyz_xy = cbuffer.data(idx_fd + 43);

    auto t_yyz_xz = cbuffer.data(idx_fd + 44);

    auto t_yyz_yy = cbuffer.data(idx_fd + 45);

    auto t_yyz_yz = cbuffer.data(idx_fd + 46);

    auto t_yyz_zz = cbuffer.data(idx_fd + 47);

    auto t_yzz_xx = cbuffer.data(idx_fd + 48);

    auto t_yzz_xy = cbuffer.data(idx_fd + 49);

    auto t_yzz_xz = cbuffer.data(idx_fd + 50);

    auto t_yzz_yy = cbuffer.data(idx_fd + 51);

    auto t_yzz_yz = cbuffer.data(idx_fd + 52);

    auto t_yzz_zz = cbuffer.data(idx_fd + 53);

    auto t_zzz_xx = cbuffer.data(idx_fd + 54);

    auto t_zzz_xy = cbuffer.data(idx_fd + 55);

    auto t_zzz_xz = cbuffer.data(idx_fd + 56);

    auto t_zzz_yy = cbuffer.data(idx_fd + 57);

    auto t_zzz_yz = cbuffer.data(idx_fd + 58);

    auto t_zzz_zz = cbuffer.data(idx_fd + 59);

    #pragma omp simd aligned(ab_x, ab_y, ab_z, t_xxx_x, t_xxx_xx, t_xxx_xy, t_xxx_xz, t_xxx_y, t_xxx_yy, t_xxx_yz, t_xxx_z, t_xxx_zz, t_xxxx_x, t_xxxx_y, t_xxxx_z, t_xxxy_x, t_xxxy_y, t_xxxy_z, t_xxxz_x, t_xxxz_y, t_xxxz_z, t_xxy_x, t_xxy_xx, t_xxy_xy, t_xxy_xz, t_xxy_y, t_xxy_yy, t_xxy_yz, t_xxy_z, t_xxy_zz, t_xxyy_x, t_xxyy_y, t_xxyy_z, t_xxyz_x, t_xxyz_y, t_xxyz_z, t_xxz_x, t_xxz_xx, t_xxz_xy, t_xxz_xz, t_xxz_y, t_xxz_yy, t_xxz_yz, t_xxz_z, t_xxz_zz, t_xxzz_x, t_xxzz_y, t_xxzz_z, t_xyy_x, t_xyy_xx, t_xyy_xy, t_xyy_xz, t_xyy_y, t_xyy_yy, t_xyy_yz, t_xyy_z, t_xyy_zz, t_xyyy_x, t_xyyy_y, t_xyyy_z, t_xyyz_x, t_xyyz_y, t_xyyz_z, t_xyz_x, t_xyz_xx, t_xyz_xy, t_xyz_xz, t_xyz_y, t_xyz_yy, t_xyz_yz, t_xyz_z, t_xyz_zz, t_xyzz_x, t_xyzz_y, t_xyzz_z, t_xzz_x, t_xzz_xx, t_xzz_xy, t_xzz_xz, t_xzz_y, t_xzz_yy, t_xzz_yz, t_xzz_z, t_xzz_zz, t_xzzz_x, t_xzzz_y, t_xzzz_z, t_yyy_x, t_yyy_xx, t_yyy_xy, t_yyy_xz, t_yyy_y, t_yyy_yy, t_yyy_yz, t_yyy_z, t_yyy_zz, t_yyyy_y, t_yyyy_z, t_yyyz_y, t_yyyz_z, t_yyz_x, t_yyz_xx, t_yyz_xy, t_yyz_xz, t_yyz_y, t_yyz_yy, t_yyz_yz, t_yyz_z, t_yyz_zz, t_yyzz_y, t_yyzz_z, t_yzz_x, t_yzz_xx, t_yzz_xy, t_yzz_xz, t_yzz_y, t_yzz_yy, t_yzz_yz, t_yzz_z, t_yzz_zz, t_yzzz_y, t_yzzz_z, t_zzz_x, t_zzz_xx, t_zzz_xy, t_zzz_xz, t_zzz_y, t_zzz_yy, t_zzz_yz, t_zzz_z, t_zzz_zz, t_zzzz_z  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        t_xxx_xx[i] = t_xxx_x[i] * ab_x[i] + t_xxxx_x[i];

        t_xxx_xy[i] = t_xxx_y[i] * ab_x[i] + t_xxxx_y[i];

        t_xxx_xz[i] = t_xxx_z[i] * ab_x[i] + t_xxxx_z[i];

        t_xxx_yy[i] = t_xxx_y[i] * ab_y[i] + t_xxxy_y[i];

        t_xxx_yz[i] = t_xxx_z[i] * ab_y[i] + t_xxxy_z[i];

        t_xxx_zz[i] = t_xxx_z[i] * ab_z[i] + t_xxxz_z[i];

        t_xxy_xx[i] = t_xxy_x[i] * ab_x[i] + t_xxxy_x[i];

        t_xxy_xy[i] = t_xxy_y[i] * ab_x[i] + t_xxxy_y[i];

        t_xxy_xz[i] = t_xxy_z[i] * ab_x[i] + t_xxxy_z[i];

        t_xxy_yy[i] = t_xxy_y[i] * ab_y[i] + t_xxyy_y[i];

        t_xxy_yz[i] = t_xxy_z[i] * ab_y[i] + t_xxyy_z[i];

        t_xxy_zz[i] = t_xxy_z[i] * ab_z[i] + t_xxyz_z[i];

        t_xxz_xx[i] = t_xxz_x[i] * ab_x[i] + t_xxxz_x[i];

        t_xxz_xy[i] = t_xxz_y[i] * ab_x[i] + t_xxxz_y[i];

        t_xxz_xz[i] = t_xxz_z[i] * ab_x[i] + t_xxxz_z[i];

        t_xxz_yy[i] = t_xxz_y[i] * ab_y[i] + t_xxyz_y[i];

        t_xxz_yz[i] = t_xxz_z[i] * ab_y[i] + t_xxyz_z[i];

        t_xxz_zz[i] = t_xxz_z[i] * ab_z[i] + t_xxzz_z[i];

        t_xyy_xx[i] = t_xyy_x[i] * ab_x[i] + t_xxyy_x[i];

        t_xyy_xy[i] = t_xyy_y[i] * ab_x[i] + t_xxyy_y[i];

        t_xyy_xz[i] = t_xyy_z[i] * ab_x[i] + t_xxyy_z[i];

        t_xyy_yy[i] = t_xyy_y[i] * ab_y[i] + t_xyyy_y[i];

        t_xyy_yz[i] = t_xyy_z[i] * ab_y[i] + t_xyyy_z[i];

        t_xyy_zz[i] = t_xyy_z[i] * ab_z[i] + t_xyyz_z[i];

        t_xyz_xx[i] = t_xyz_x[i] * ab_x[i] + t_xxyz_x[i];

        t_xyz_xy[i] = t_xyz_y[i] * ab_x[i] + t_xxyz_y[i];

        t_xyz_xz[i] = t_xyz_z[i] * ab_x[i] + t_xxyz_z[i];

        t_xyz_yy[i] = t_xyz_y[i] * ab_y[i] + t_xyyz_y[i];

        t_xyz_yz[i] = t_xyz_z[i] * ab_y[i] + t_xyyz_z[i];

        t_xyz_zz[i] = t_xyz_z[i] * ab_z[i] + t_xyzz_z[i];

        t_xzz_xx[i] = t_xzz_x[i] * ab_x[i] + t_xxzz_x[i];

        t_xzz_xy[i] = t_xzz_y[i] * ab_x[i] + t_xxzz_y[i];

        t_xzz_xz[i] = t_xzz_z[i] * ab_x[i] + t_xxzz_z[i];

        t_xzz_yy[i] = t_xzz_y[i] * ab_y[i] + t_xyzz_y[i];

        t_xzz_yz[i] = t_xzz_z[i] * ab_y[i] + t_xyzz_z[i];

        t_xzz_zz[i] = t_xzz_z[i] * ab_z[i] + t_xzzz_z[i];

        t_yyy_xx[i] = t_yyy_x[i] * ab_x[i] + t_xyyy_x[i];

        t_yyy_xy[i] = t_yyy_y[i] * ab_x[i] + t_xyyy_y[i];

        t_yyy_xz[i] = t_yyy_z[i] * ab_x[i] + t_xyyy_z[i];

        t_yyy_yy[i] = t_yyy_y[i] * ab_y[i] + t_yyyy_y[i];

        t_yyy_yz[i] = t_yyy_z[i] * ab_y[i] + t_yyyy_z[i];

        t_yyy_zz[i] = t_yyy_z[i] * ab_z[i] + t_yyyz_z[i];

        t_yyz_xx[i] = t_yyz_x[i] * ab_x[i] + t_xyyz_x[i];

        t_yyz_xy[i] = t_yyz_y[i] * ab_x[i] + t_xyyz_y[i];

        t_yyz_xz[i] = t_yyz_z[i] * ab_x[i] + t_xyyz_z[i];

        t_yyz_yy[i] = t_yyz_y[i] * ab_y[i] + t_yyyz_y[i];

        t_yyz_yz[i] = t_yyz_z[i] * ab_y[i] + t_yyyz_z[i];

        t_yyz_zz[i] = t_yyz_z[i] * ab_z[i] + t_yyzz_z[i];

        t_yzz_xx[i] = t_yzz_x[i] * ab_x[i] + t_xyzz_x[i];

        t_yzz_xy[i] = t_yzz_y[i] * ab_x[i] + t_xyzz_y[i];

        t_yzz_xz[i] = t_yzz_z[i] * ab_x[i] + t_xyzz_z[i];

        t_yzz_yy[i] = t_yzz_y[i] * ab_y[i] + t_yyzz_y[i];

        t_yzz_yz[i] = t_yzz_z[i] * ab_y[i] + t_yyzz_z[i];

        t_yzz_zz[i] = t_yzz_z[i] * ab_z[i] + t_yzzz_z[i];

        t_zzz_xx[i] = t_zzz_x[i] * ab_x[i] + t_xzzz_x[i];

        t_zzz_xy[i] = t_zzz_y[i] * ab_x[i] + t_xzzz_y[i];

        t_zzz_xz[i] = t_zzz_z[i] * ab_x[i] + t_xzzz_z[i];

        t_zzz_yy[i] = t_zzz_y[i] * ab_y[i] + t_yzzz_y[i];

        t_zzz_yz[i] = t_zzz_z[i] * ab_y[i] + t_yzzz_z[i];

        t_zzz_zz[i] = t_zzz_z[i] * ab_z[i] + t_zzzz_z[i];
    }
}

} // t2chrr namespace


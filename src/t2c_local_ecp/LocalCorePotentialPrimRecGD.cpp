#include "LocalCorePotentialPrimRecGD.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_gd(CSimdArray<double>& pbuffer, 
                                  const size_t idx_gd,
                                  const size_t idx_dd,
                                  const size_t idx_fp,
                                  const size_t idx_fd,
                                  const CSimdArray<double>& factors,
                                  const size_t idx_ra,
                                  const size_t idx_zeta) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up R(RA) distances

    auto ra_x = factors.data(idx_ra);

    auto ra_y = factors.data(idx_ra + 1);

    auto ra_z = factors.data(idx_ra + 2);

    // Set up inverted 1/2xi

    auto fxi = factors.data(idx_zeta);

    // Set up components of auxiliary buffer : DD

    auto tg_xx_xx = pbuffer.data(idx_dd);

    auto tg_xx_xy = pbuffer.data(idx_dd + 1);

    auto tg_xx_xz = pbuffer.data(idx_dd + 2);

    auto tg_xx_yy = pbuffer.data(idx_dd + 3);

    auto tg_xx_yz = pbuffer.data(idx_dd + 4);

    auto tg_xx_zz = pbuffer.data(idx_dd + 5);

    auto tg_xy_yy = pbuffer.data(idx_dd + 9);

    auto tg_xy_yz = pbuffer.data(idx_dd + 10);

    auto tg_xz_yz = pbuffer.data(idx_dd + 16);

    auto tg_xz_zz = pbuffer.data(idx_dd + 17);

    auto tg_yy_xx = pbuffer.data(idx_dd + 18);

    auto tg_yy_xy = pbuffer.data(idx_dd + 19);

    auto tg_yy_xz = pbuffer.data(idx_dd + 20);

    auto tg_yy_yy = pbuffer.data(idx_dd + 21);

    auto tg_yy_yz = pbuffer.data(idx_dd + 22);

    auto tg_yy_zz = pbuffer.data(idx_dd + 23);

    auto tg_yz_xz = pbuffer.data(idx_dd + 26);

    auto tg_yz_yz = pbuffer.data(idx_dd + 28);

    auto tg_yz_zz = pbuffer.data(idx_dd + 29);

    auto tg_zz_xx = pbuffer.data(idx_dd + 30);

    auto tg_zz_xy = pbuffer.data(idx_dd + 31);

    auto tg_zz_xz = pbuffer.data(idx_dd + 32);

    auto tg_zz_yy = pbuffer.data(idx_dd + 33);

    auto tg_zz_yz = pbuffer.data(idx_dd + 34);

    auto tg_zz_zz = pbuffer.data(idx_dd + 35);

    // Set up components of auxiliary buffer : FP

    auto tg_xxx_x = pbuffer.data(idx_fp);

    auto tg_xxx_y = pbuffer.data(idx_fp + 1);

    auto tg_xxx_z = pbuffer.data(idx_fp + 2);

    auto tg_xyy_y = pbuffer.data(idx_fp + 10);

    auto tg_xzz_z = pbuffer.data(idx_fp + 17);

    auto tg_yyy_x = pbuffer.data(idx_fp + 18);

    auto tg_yyy_y = pbuffer.data(idx_fp + 19);

    auto tg_yyy_z = pbuffer.data(idx_fp + 20);

    auto tg_yyz_z = pbuffer.data(idx_fp + 23);

    auto tg_yzz_y = pbuffer.data(idx_fp + 25);

    auto tg_yzz_z = pbuffer.data(idx_fp + 26);

    auto tg_zzz_x = pbuffer.data(idx_fp + 27);

    auto tg_zzz_y = pbuffer.data(idx_fp + 28);

    auto tg_zzz_z = pbuffer.data(idx_fp + 29);

    // Set up components of auxiliary buffer : FD

    auto tg_xxx_xx = pbuffer.data(idx_fd);

    auto tg_xxx_xy = pbuffer.data(idx_fd + 1);

    auto tg_xxx_xz = pbuffer.data(idx_fd + 2);

    auto tg_xxx_yy = pbuffer.data(idx_fd + 3);

    auto tg_xxx_yz = pbuffer.data(idx_fd + 4);

    auto tg_xxx_zz = pbuffer.data(idx_fd + 5);

    auto tg_xxy_xx = pbuffer.data(idx_fd + 6);

    auto tg_xxy_xy = pbuffer.data(idx_fd + 7);

    auto tg_xxy_xz = pbuffer.data(idx_fd + 8);

    auto tg_xxy_yy = pbuffer.data(idx_fd + 9);

    auto tg_xxy_yz = pbuffer.data(idx_fd + 10);

    auto tg_xxz_xx = pbuffer.data(idx_fd + 12);

    auto tg_xxz_xy = pbuffer.data(idx_fd + 13);

    auto tg_xxz_xz = pbuffer.data(idx_fd + 14);

    auto tg_xxz_yz = pbuffer.data(idx_fd + 16);

    auto tg_xxz_zz = pbuffer.data(idx_fd + 17);

    auto tg_xyy_xx = pbuffer.data(idx_fd + 18);

    auto tg_xyy_xy = pbuffer.data(idx_fd + 19);

    auto tg_xyy_yy = pbuffer.data(idx_fd + 21);

    auto tg_xyy_yz = pbuffer.data(idx_fd + 22);

    auto tg_xyy_zz = pbuffer.data(idx_fd + 23);

    auto tg_xyz_yz = pbuffer.data(idx_fd + 28);

    auto tg_xzz_xx = pbuffer.data(idx_fd + 30);

    auto tg_xzz_xz = pbuffer.data(idx_fd + 32);

    auto tg_xzz_yy = pbuffer.data(idx_fd + 33);

    auto tg_xzz_yz = pbuffer.data(idx_fd + 34);

    auto tg_xzz_zz = pbuffer.data(idx_fd + 35);

    auto tg_yyy_xx = pbuffer.data(idx_fd + 36);

    auto tg_yyy_xy = pbuffer.data(idx_fd + 37);

    auto tg_yyy_xz = pbuffer.data(idx_fd + 38);

    auto tg_yyy_yy = pbuffer.data(idx_fd + 39);

    auto tg_yyy_yz = pbuffer.data(idx_fd + 40);

    auto tg_yyy_zz = pbuffer.data(idx_fd + 41);

    auto tg_yyz_xy = pbuffer.data(idx_fd + 43);

    auto tg_yyz_xz = pbuffer.data(idx_fd + 44);

    auto tg_yyz_yy = pbuffer.data(idx_fd + 45);

    auto tg_yyz_yz = pbuffer.data(idx_fd + 46);

    auto tg_yyz_zz = pbuffer.data(idx_fd + 47);

    auto tg_yzz_xx = pbuffer.data(idx_fd + 48);

    auto tg_yzz_xy = pbuffer.data(idx_fd + 49);

    auto tg_yzz_xz = pbuffer.data(idx_fd + 50);

    auto tg_yzz_yy = pbuffer.data(idx_fd + 51);

    auto tg_yzz_yz = pbuffer.data(idx_fd + 52);

    auto tg_yzz_zz = pbuffer.data(idx_fd + 53);

    auto tg_zzz_xx = pbuffer.data(idx_fd + 54);

    auto tg_zzz_xy = pbuffer.data(idx_fd + 55);

    auto tg_zzz_xz = pbuffer.data(idx_fd + 56);

    auto tg_zzz_yy = pbuffer.data(idx_fd + 57);

    auto tg_zzz_yz = pbuffer.data(idx_fd + 58);

    auto tg_zzz_zz = pbuffer.data(idx_fd + 59);

    // Set up components of targeted buffer : GD

    auto tg_xxxx_xx = pbuffer.data(idx_gd);

    auto tg_xxxx_xy = pbuffer.data(idx_gd + 1);

    auto tg_xxxx_xz = pbuffer.data(idx_gd + 2);

    auto tg_xxxx_yy = pbuffer.data(idx_gd + 3);

    auto tg_xxxx_yz = pbuffer.data(idx_gd + 4);

    auto tg_xxxx_zz = pbuffer.data(idx_gd + 5);

    auto tg_xxxy_xx = pbuffer.data(idx_gd + 6);

    auto tg_xxxy_xy = pbuffer.data(idx_gd + 7);

    auto tg_xxxy_xz = pbuffer.data(idx_gd + 8);

    auto tg_xxxy_yy = pbuffer.data(idx_gd + 9);

    auto tg_xxxy_yz = pbuffer.data(idx_gd + 10);

    auto tg_xxxy_zz = pbuffer.data(idx_gd + 11);

    auto tg_xxxz_xx = pbuffer.data(idx_gd + 12);

    auto tg_xxxz_xy = pbuffer.data(idx_gd + 13);

    auto tg_xxxz_xz = pbuffer.data(idx_gd + 14);

    auto tg_xxxz_yy = pbuffer.data(idx_gd + 15);

    auto tg_xxxz_yz = pbuffer.data(idx_gd + 16);

    auto tg_xxxz_zz = pbuffer.data(idx_gd + 17);

    auto tg_xxyy_xx = pbuffer.data(idx_gd + 18);

    auto tg_xxyy_xy = pbuffer.data(idx_gd + 19);

    auto tg_xxyy_xz = pbuffer.data(idx_gd + 20);

    auto tg_xxyy_yy = pbuffer.data(idx_gd + 21);

    auto tg_xxyy_yz = pbuffer.data(idx_gd + 22);

    auto tg_xxyy_zz = pbuffer.data(idx_gd + 23);

    auto tg_xxyz_xx = pbuffer.data(idx_gd + 24);

    auto tg_xxyz_xy = pbuffer.data(idx_gd + 25);

    auto tg_xxyz_xz = pbuffer.data(idx_gd + 26);

    auto tg_xxyz_yy = pbuffer.data(idx_gd + 27);

    auto tg_xxyz_yz = pbuffer.data(idx_gd + 28);

    auto tg_xxyz_zz = pbuffer.data(idx_gd + 29);

    auto tg_xxzz_xx = pbuffer.data(idx_gd + 30);

    auto tg_xxzz_xy = pbuffer.data(idx_gd + 31);

    auto tg_xxzz_xz = pbuffer.data(idx_gd + 32);

    auto tg_xxzz_yy = pbuffer.data(idx_gd + 33);

    auto tg_xxzz_yz = pbuffer.data(idx_gd + 34);

    auto tg_xxzz_zz = pbuffer.data(idx_gd + 35);

    auto tg_xyyy_xx = pbuffer.data(idx_gd + 36);

    auto tg_xyyy_xy = pbuffer.data(idx_gd + 37);

    auto tg_xyyy_xz = pbuffer.data(idx_gd + 38);

    auto tg_xyyy_yy = pbuffer.data(idx_gd + 39);

    auto tg_xyyy_yz = pbuffer.data(idx_gd + 40);

    auto tg_xyyy_zz = pbuffer.data(idx_gd + 41);

    auto tg_xyyz_xx = pbuffer.data(idx_gd + 42);

    auto tg_xyyz_xy = pbuffer.data(idx_gd + 43);

    auto tg_xyyz_xz = pbuffer.data(idx_gd + 44);

    auto tg_xyyz_yy = pbuffer.data(idx_gd + 45);

    auto tg_xyyz_yz = pbuffer.data(idx_gd + 46);

    auto tg_xyyz_zz = pbuffer.data(idx_gd + 47);

    auto tg_xyzz_xx = pbuffer.data(idx_gd + 48);

    auto tg_xyzz_xy = pbuffer.data(idx_gd + 49);

    auto tg_xyzz_xz = pbuffer.data(idx_gd + 50);

    auto tg_xyzz_yy = pbuffer.data(idx_gd + 51);

    auto tg_xyzz_yz = pbuffer.data(idx_gd + 52);

    auto tg_xyzz_zz = pbuffer.data(idx_gd + 53);

    auto tg_xzzz_xx = pbuffer.data(idx_gd + 54);

    auto tg_xzzz_xy = pbuffer.data(idx_gd + 55);

    auto tg_xzzz_xz = pbuffer.data(idx_gd + 56);

    auto tg_xzzz_yy = pbuffer.data(idx_gd + 57);

    auto tg_xzzz_yz = pbuffer.data(idx_gd + 58);

    auto tg_xzzz_zz = pbuffer.data(idx_gd + 59);

    auto tg_yyyy_xx = pbuffer.data(idx_gd + 60);

    auto tg_yyyy_xy = pbuffer.data(idx_gd + 61);

    auto tg_yyyy_xz = pbuffer.data(idx_gd + 62);

    auto tg_yyyy_yy = pbuffer.data(idx_gd + 63);

    auto tg_yyyy_yz = pbuffer.data(idx_gd + 64);

    auto tg_yyyy_zz = pbuffer.data(idx_gd + 65);

    auto tg_yyyz_xx = pbuffer.data(idx_gd + 66);

    auto tg_yyyz_xy = pbuffer.data(idx_gd + 67);

    auto tg_yyyz_xz = pbuffer.data(idx_gd + 68);

    auto tg_yyyz_yy = pbuffer.data(idx_gd + 69);

    auto tg_yyyz_yz = pbuffer.data(idx_gd + 70);

    auto tg_yyyz_zz = pbuffer.data(idx_gd + 71);

    auto tg_yyzz_xx = pbuffer.data(idx_gd + 72);

    auto tg_yyzz_xy = pbuffer.data(idx_gd + 73);

    auto tg_yyzz_xz = pbuffer.data(idx_gd + 74);

    auto tg_yyzz_yy = pbuffer.data(idx_gd + 75);

    auto tg_yyzz_yz = pbuffer.data(idx_gd + 76);

    auto tg_yyzz_zz = pbuffer.data(idx_gd + 77);

    auto tg_yzzz_xx = pbuffer.data(idx_gd + 78);

    auto tg_yzzz_xy = pbuffer.data(idx_gd + 79);

    auto tg_yzzz_xz = pbuffer.data(idx_gd + 80);

    auto tg_yzzz_yy = pbuffer.data(idx_gd + 81);

    auto tg_yzzz_yz = pbuffer.data(idx_gd + 82);

    auto tg_yzzz_zz = pbuffer.data(idx_gd + 83);

    auto tg_zzzz_xx = pbuffer.data(idx_gd + 84);

    auto tg_zzzz_xy = pbuffer.data(idx_gd + 85);

    auto tg_zzzz_xz = pbuffer.data(idx_gd + 86);

    auto tg_zzzz_yy = pbuffer.data(idx_gd + 87);

    auto tg_zzzz_yz = pbuffer.data(idx_gd + 88);

    auto tg_zzzz_zz = pbuffer.data(idx_gd + 89);

    #pragma omp simd aligned(fxi, ra_x, ra_y, ra_z, tg_xx_xx, tg_xx_xy, tg_xx_xz, tg_xx_yy, tg_xx_yz, tg_xx_zz, tg_xxx_x, tg_xxx_xx, tg_xxx_xy, tg_xxx_xz, tg_xxx_y, tg_xxx_yy, tg_xxx_yz, tg_xxx_z, tg_xxx_zz, tg_xxxx_xx, tg_xxxx_xy, tg_xxxx_xz, tg_xxxx_yy, tg_xxxx_yz, tg_xxxx_zz, tg_xxxy_xx, tg_xxxy_xy, tg_xxxy_xz, tg_xxxy_yy, tg_xxxy_yz, tg_xxxy_zz, tg_xxxz_xx, tg_xxxz_xy, tg_xxxz_xz, tg_xxxz_yy, tg_xxxz_yz, tg_xxxz_zz, tg_xxy_xx, tg_xxy_xy, tg_xxy_xz, tg_xxy_yy, tg_xxy_yz, tg_xxyy_xx, tg_xxyy_xy, tg_xxyy_xz, tg_xxyy_yy, tg_xxyy_yz, tg_xxyy_zz, tg_xxyz_xx, tg_xxyz_xy, tg_xxyz_xz, tg_xxyz_yy, tg_xxyz_yz, tg_xxyz_zz, tg_xxz_xx, tg_xxz_xy, tg_xxz_xz, tg_xxz_yz, tg_xxz_zz, tg_xxzz_xx, tg_xxzz_xy, tg_xxzz_xz, tg_xxzz_yy, tg_xxzz_yz, tg_xxzz_zz, tg_xy_yy, tg_xy_yz, tg_xyy_xx, tg_xyy_xy, tg_xyy_y, tg_xyy_yy, tg_xyy_yz, tg_xyy_zz, tg_xyyy_xx, tg_xyyy_xy, tg_xyyy_xz, tg_xyyy_yy, tg_xyyy_yz, tg_xyyy_zz, tg_xyyz_xx, tg_xyyz_xy, tg_xyyz_xz, tg_xyyz_yy, tg_xyyz_yz, tg_xyyz_zz, tg_xyz_yz, tg_xyzz_xx, tg_xyzz_xy, tg_xyzz_xz, tg_xyzz_yy, tg_xyzz_yz, tg_xyzz_zz, tg_xz_yz, tg_xz_zz, tg_xzz_xx, tg_xzz_xz, tg_xzz_yy, tg_xzz_yz, tg_xzz_z, tg_xzz_zz, tg_xzzz_xx, tg_xzzz_xy, tg_xzzz_xz, tg_xzzz_yy, tg_xzzz_yz, tg_xzzz_zz, tg_yy_xx, tg_yy_xy, tg_yy_xz, tg_yy_yy, tg_yy_yz, tg_yy_zz, tg_yyy_x, tg_yyy_xx, tg_yyy_xy, tg_yyy_xz, tg_yyy_y, tg_yyy_yy, tg_yyy_yz, tg_yyy_z, tg_yyy_zz, tg_yyyy_xx, tg_yyyy_xy, tg_yyyy_xz, tg_yyyy_yy, tg_yyyy_yz, tg_yyyy_zz, tg_yyyz_xx, tg_yyyz_xy, tg_yyyz_xz, tg_yyyz_yy, tg_yyyz_yz, tg_yyyz_zz, tg_yyz_xy, tg_yyz_xz, tg_yyz_yy, tg_yyz_yz, tg_yyz_z, tg_yyz_zz, tg_yyzz_xx, tg_yyzz_xy, tg_yyzz_xz, tg_yyzz_yy, tg_yyzz_yz, tg_yyzz_zz, tg_yz_xz, tg_yz_yz, tg_yz_zz, tg_yzz_xx, tg_yzz_xy, tg_yzz_xz, tg_yzz_y, tg_yzz_yy, tg_yzz_yz, tg_yzz_z, tg_yzz_zz, tg_yzzz_xx, tg_yzzz_xy, tg_yzzz_xz, tg_yzzz_yy, tg_yzzz_yz, tg_yzzz_zz, tg_zz_xx, tg_zz_xy, tg_zz_xz, tg_zz_yy, tg_zz_yz, tg_zz_zz, tg_zzz_x, tg_zzz_xx, tg_zzz_xy, tg_zzz_xz, tg_zzz_y, tg_zzz_yy, tg_zzz_yz, tg_zzz_z, tg_zzz_zz, tg_zzzz_xx, tg_zzzz_xy, tg_zzzz_xz, tg_zzzz_yy, tg_zzzz_yz, tg_zzzz_zz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_xxxx_xx[i] = 3.0 * tg_xx_xx[i] * fxi[i] + 2.0 * tg_xxx_x[i] * fxi[i] + tg_xxx_xx[i] * ra_x[i];

        tg_xxxx_xy[i] = 3.0 * tg_xx_xy[i] * fxi[i] + tg_xxx_y[i] * fxi[i] + tg_xxx_xy[i] * ra_x[i];

        tg_xxxx_xz[i] = 3.0 * tg_xx_xz[i] * fxi[i] + tg_xxx_z[i] * fxi[i] + tg_xxx_xz[i] * ra_x[i];

        tg_xxxx_yy[i] = 3.0 * tg_xx_yy[i] * fxi[i] + tg_xxx_yy[i] * ra_x[i];

        tg_xxxx_yz[i] = 3.0 * tg_xx_yz[i] * fxi[i] + tg_xxx_yz[i] * ra_x[i];

        tg_xxxx_zz[i] = 3.0 * tg_xx_zz[i] * fxi[i] + tg_xxx_zz[i] * ra_x[i];

        tg_xxxy_xx[i] = tg_xxx_xx[i] * ra_y[i];

        tg_xxxy_xy[i] = tg_xxx_x[i] * fxi[i] + tg_xxx_xy[i] * ra_y[i];

        tg_xxxy_xz[i] = tg_xxx_xz[i] * ra_y[i];

        tg_xxxy_yy[i] = 2.0 * tg_xy_yy[i] * fxi[i] + tg_xxy_yy[i] * ra_x[i];

        tg_xxxy_yz[i] = 2.0 * tg_xy_yz[i] * fxi[i] + tg_xxy_yz[i] * ra_x[i];

        tg_xxxy_zz[i] = tg_xxx_zz[i] * ra_y[i];

        tg_xxxz_xx[i] = tg_xxx_xx[i] * ra_z[i];

        tg_xxxz_xy[i] = tg_xxx_xy[i] * ra_z[i];

        tg_xxxz_xz[i] = tg_xxx_x[i] * fxi[i] + tg_xxx_xz[i] * ra_z[i];

        tg_xxxz_yy[i] = tg_xxx_yy[i] * ra_z[i];

        tg_xxxz_yz[i] = 2.0 * tg_xz_yz[i] * fxi[i] + tg_xxz_yz[i] * ra_x[i];

        tg_xxxz_zz[i] = 2.0 * tg_xz_zz[i] * fxi[i] + tg_xxz_zz[i] * ra_x[i];

        tg_xxyy_xx[i] = tg_xx_xx[i] * fxi[i] + tg_xxy_xx[i] * ra_y[i];

        tg_xxyy_xy[i] = tg_yy_xy[i] * fxi[i] + tg_xyy_y[i] * fxi[i] + tg_xyy_xy[i] * ra_x[i];

        tg_xxyy_xz[i] = tg_xx_xz[i] * fxi[i] + tg_xxy_xz[i] * ra_y[i];

        tg_xxyy_yy[i] = tg_yy_yy[i] * fxi[i] + tg_xyy_yy[i] * ra_x[i];

        tg_xxyy_yz[i] = tg_yy_yz[i] * fxi[i] + tg_xyy_yz[i] * ra_x[i];

        tg_xxyy_zz[i] = tg_yy_zz[i] * fxi[i] + tg_xyy_zz[i] * ra_x[i];

        tg_xxyz_xx[i] = tg_xxz_xx[i] * ra_y[i];

        tg_xxyz_xy[i] = tg_xxy_xy[i] * ra_z[i];

        tg_xxyz_xz[i] = tg_xxz_xz[i] * ra_y[i];

        tg_xxyz_yy[i] = tg_xxy_yy[i] * ra_z[i];

        tg_xxyz_yz[i] = tg_yz_yz[i] * fxi[i] + tg_xyz_yz[i] * ra_x[i];

        tg_xxyz_zz[i] = tg_xxz_zz[i] * ra_y[i];

        tg_xxzz_xx[i] = tg_xx_xx[i] * fxi[i] + tg_xxz_xx[i] * ra_z[i];

        tg_xxzz_xy[i] = tg_xx_xy[i] * fxi[i] + tg_xxz_xy[i] * ra_z[i];

        tg_xxzz_xz[i] = tg_zz_xz[i] * fxi[i] + tg_xzz_z[i] * fxi[i] + tg_xzz_xz[i] * ra_x[i];

        tg_xxzz_yy[i] = tg_zz_yy[i] * fxi[i] + tg_xzz_yy[i] * ra_x[i];

        tg_xxzz_yz[i] = tg_zz_yz[i] * fxi[i] + tg_xzz_yz[i] * ra_x[i];

        tg_xxzz_zz[i] = tg_zz_zz[i] * fxi[i] + tg_xzz_zz[i] * ra_x[i];

        tg_xyyy_xx[i] = 2.0 * tg_yyy_x[i] * fxi[i] + tg_yyy_xx[i] * ra_x[i];

        tg_xyyy_xy[i] = tg_yyy_y[i] * fxi[i] + tg_yyy_xy[i] * ra_x[i];

        tg_xyyy_xz[i] = tg_yyy_z[i] * fxi[i] + tg_yyy_xz[i] * ra_x[i];

        tg_xyyy_yy[i] = tg_yyy_yy[i] * ra_x[i];

        tg_xyyy_yz[i] = tg_yyy_yz[i] * ra_x[i];

        tg_xyyy_zz[i] = tg_yyy_zz[i] * ra_x[i];

        tg_xyyz_xx[i] = tg_xyy_xx[i] * ra_z[i];

        tg_xyyz_xy[i] = tg_xyy_xy[i] * ra_z[i];

        tg_xyyz_xz[i] = tg_yyz_z[i] * fxi[i] + tg_yyz_xz[i] * ra_x[i];

        tg_xyyz_yy[i] = tg_yyz_yy[i] * ra_x[i];

        tg_xyyz_yz[i] = tg_yyz_yz[i] * ra_x[i];

        tg_xyyz_zz[i] = tg_yyz_zz[i] * ra_x[i];

        tg_xyzz_xx[i] = tg_xzz_xx[i] * ra_y[i];

        tg_xyzz_xy[i] = tg_yzz_y[i] * fxi[i] + tg_yzz_xy[i] * ra_x[i];

        tg_xyzz_xz[i] = tg_xzz_xz[i] * ra_y[i];

        tg_xyzz_yy[i] = tg_yzz_yy[i] * ra_x[i];

        tg_xyzz_yz[i] = tg_yzz_yz[i] * ra_x[i];

        tg_xyzz_zz[i] = tg_yzz_zz[i] * ra_x[i];

        tg_xzzz_xx[i] = 2.0 * tg_zzz_x[i] * fxi[i] + tg_zzz_xx[i] * ra_x[i];

        tg_xzzz_xy[i] = tg_zzz_y[i] * fxi[i] + tg_zzz_xy[i] * ra_x[i];

        tg_xzzz_xz[i] = tg_zzz_z[i] * fxi[i] + tg_zzz_xz[i] * ra_x[i];

        tg_xzzz_yy[i] = tg_zzz_yy[i] * ra_x[i];

        tg_xzzz_yz[i] = tg_zzz_yz[i] * ra_x[i];

        tg_xzzz_zz[i] = tg_zzz_zz[i] * ra_x[i];

        tg_yyyy_xx[i] = 3.0 * tg_yy_xx[i] * fxi[i] + tg_yyy_xx[i] * ra_y[i];

        tg_yyyy_xy[i] = 3.0 * tg_yy_xy[i] * fxi[i] + tg_yyy_x[i] * fxi[i] + tg_yyy_xy[i] * ra_y[i];

        tg_yyyy_xz[i] = 3.0 * tg_yy_xz[i] * fxi[i] + tg_yyy_xz[i] * ra_y[i];

        tg_yyyy_yy[i] = 3.0 * tg_yy_yy[i] * fxi[i] + 2.0 * tg_yyy_y[i] * fxi[i] + tg_yyy_yy[i] * ra_y[i];

        tg_yyyy_yz[i] = 3.0 * tg_yy_yz[i] * fxi[i] + tg_yyy_z[i] * fxi[i] + tg_yyy_yz[i] * ra_y[i];

        tg_yyyy_zz[i] = 3.0 * tg_yy_zz[i] * fxi[i] + tg_yyy_zz[i] * ra_y[i];

        tg_yyyz_xx[i] = tg_yyy_xx[i] * ra_z[i];

        tg_yyyz_xy[i] = tg_yyy_xy[i] * ra_z[i];

        tg_yyyz_xz[i] = 2.0 * tg_yz_xz[i] * fxi[i] + tg_yyz_xz[i] * ra_y[i];

        tg_yyyz_yy[i] = tg_yyy_yy[i] * ra_z[i];

        tg_yyyz_yz[i] = tg_yyy_y[i] * fxi[i] + tg_yyy_yz[i] * ra_z[i];

        tg_yyyz_zz[i] = 2.0 * tg_yz_zz[i] * fxi[i] + tg_yyz_zz[i] * ra_y[i];

        tg_yyzz_xx[i] = tg_zz_xx[i] * fxi[i] + tg_yzz_xx[i] * ra_y[i];

        tg_yyzz_xy[i] = tg_yy_xy[i] * fxi[i] + tg_yyz_xy[i] * ra_z[i];

        tg_yyzz_xz[i] = tg_zz_xz[i] * fxi[i] + tg_yzz_xz[i] * ra_y[i];

        tg_yyzz_yy[i] = tg_yy_yy[i] * fxi[i] + tg_yyz_yy[i] * ra_z[i];

        tg_yyzz_yz[i] = tg_zz_yz[i] * fxi[i] + tg_yzz_z[i] * fxi[i] + tg_yzz_yz[i] * ra_y[i];

        tg_yyzz_zz[i] = tg_zz_zz[i] * fxi[i] + tg_yzz_zz[i] * ra_y[i];

        tg_yzzz_xx[i] = tg_zzz_xx[i] * ra_y[i];

        tg_yzzz_xy[i] = tg_zzz_x[i] * fxi[i] + tg_zzz_xy[i] * ra_y[i];

        tg_yzzz_xz[i] = tg_zzz_xz[i] * ra_y[i];

        tg_yzzz_yy[i] = 2.0 * tg_zzz_y[i] * fxi[i] + tg_zzz_yy[i] * ra_y[i];

        tg_yzzz_yz[i] = tg_zzz_z[i] * fxi[i] + tg_zzz_yz[i] * ra_y[i];

        tg_yzzz_zz[i] = tg_zzz_zz[i] * ra_y[i];

        tg_zzzz_xx[i] = 3.0 * tg_zz_xx[i] * fxi[i] + tg_zzz_xx[i] * ra_z[i];

        tg_zzzz_xy[i] = 3.0 * tg_zz_xy[i] * fxi[i] + tg_zzz_xy[i] * ra_z[i];

        tg_zzzz_xz[i] = 3.0 * tg_zz_xz[i] * fxi[i] + tg_zzz_x[i] * fxi[i] + tg_zzz_xz[i] * ra_z[i];

        tg_zzzz_yy[i] = 3.0 * tg_zz_yy[i] * fxi[i] + tg_zzz_yy[i] * ra_z[i];

        tg_zzzz_yz[i] = 3.0 * tg_zz_yz[i] * fxi[i] + tg_zzz_y[i] * fxi[i] + tg_zzz_yz[i] * ra_z[i];

        tg_zzzz_zz[i] = 3.0 * tg_zz_zz[i] * fxi[i] + 2.0 * tg_zzz_z[i] * fxi[i] + tg_zzz_zz[i] * ra_z[i];
    }
}

} // t2lecp namespace


#include "LocalCorePotentialPrimRecFF.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_ff(CSimdArray<double>& pbuffer, 
                                  const size_t idx_ff,
                                  const size_t idx_pf,
                                  const size_t idx_dd,
                                  const size_t idx_df,
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

    // Set up components of auxiliary buffer : PF

    auto tg_x_xxx = pbuffer.data(idx_pf);

    auto tg_x_xxy = pbuffer.data(idx_pf + 1);

    auto tg_x_xxz = pbuffer.data(idx_pf + 2);

    auto tg_x_xyy = pbuffer.data(idx_pf + 3);

    auto tg_x_xyz = pbuffer.data(idx_pf + 4);

    auto tg_x_xzz = pbuffer.data(idx_pf + 5);

    auto tg_x_yyy = pbuffer.data(idx_pf + 6);

    auto tg_x_yyz = pbuffer.data(idx_pf + 7);

    auto tg_x_yzz = pbuffer.data(idx_pf + 8);

    auto tg_x_zzz = pbuffer.data(idx_pf + 9);

    auto tg_y_xxx = pbuffer.data(idx_pf + 10);

    auto tg_y_xxy = pbuffer.data(idx_pf + 11);

    auto tg_y_xxz = pbuffer.data(idx_pf + 12);

    auto tg_y_xyy = pbuffer.data(idx_pf + 13);

    auto tg_y_xyz = pbuffer.data(idx_pf + 14);

    auto tg_y_xzz = pbuffer.data(idx_pf + 15);

    auto tg_y_yyy = pbuffer.data(idx_pf + 16);

    auto tg_y_yyz = pbuffer.data(idx_pf + 17);

    auto tg_y_yzz = pbuffer.data(idx_pf + 18);

    auto tg_y_zzz = pbuffer.data(idx_pf + 19);

    auto tg_z_xxx = pbuffer.data(idx_pf + 20);

    auto tg_z_xxy = pbuffer.data(idx_pf + 21);

    auto tg_z_xxz = pbuffer.data(idx_pf + 22);

    auto tg_z_xyy = pbuffer.data(idx_pf + 23);

    auto tg_z_xyz = pbuffer.data(idx_pf + 24);

    auto tg_z_xzz = pbuffer.data(idx_pf + 25);

    auto tg_z_yyy = pbuffer.data(idx_pf + 26);

    auto tg_z_yyz = pbuffer.data(idx_pf + 27);

    auto tg_z_yzz = pbuffer.data(idx_pf + 28);

    auto tg_z_zzz = pbuffer.data(idx_pf + 29);

    // Set up components of auxiliary buffer : DD

    auto tg_xx_xx = pbuffer.data(idx_dd);

    auto tg_xx_xy = pbuffer.data(idx_dd + 1);

    auto tg_xx_xz = pbuffer.data(idx_dd + 2);

    auto tg_xx_yy = pbuffer.data(idx_dd + 3);

    auto tg_xx_yz = pbuffer.data(idx_dd + 4);

    auto tg_xx_zz = pbuffer.data(idx_dd + 5);

    auto tg_yy_xx = pbuffer.data(idx_dd + 18);

    auto tg_yy_xy = pbuffer.data(idx_dd + 19);

    auto tg_yy_xz = pbuffer.data(idx_dd + 20);

    auto tg_yy_yy = pbuffer.data(idx_dd + 21);

    auto tg_yy_yz = pbuffer.data(idx_dd + 22);

    auto tg_yy_zz = pbuffer.data(idx_dd + 23);

    auto tg_yz_yz = pbuffer.data(idx_dd + 28);

    auto tg_zz_xx = pbuffer.data(idx_dd + 30);

    auto tg_zz_xy = pbuffer.data(idx_dd + 31);

    auto tg_zz_xz = pbuffer.data(idx_dd + 32);

    auto tg_zz_yy = pbuffer.data(idx_dd + 33);

    auto tg_zz_yz = pbuffer.data(idx_dd + 34);

    auto tg_zz_zz = pbuffer.data(idx_dd + 35);

    // Set up components of auxiliary buffer : DF

    auto tg_xx_xxx = pbuffer.data(idx_df);

    auto tg_xx_xxy = pbuffer.data(idx_df + 1);

    auto tg_xx_xxz = pbuffer.data(idx_df + 2);

    auto tg_xx_xyy = pbuffer.data(idx_df + 3);

    auto tg_xx_xyz = pbuffer.data(idx_df + 4);

    auto tg_xx_xzz = pbuffer.data(idx_df + 5);

    auto tg_xx_yyy = pbuffer.data(idx_df + 6);

    auto tg_xx_yyz = pbuffer.data(idx_df + 7);

    auto tg_xx_yzz = pbuffer.data(idx_df + 8);

    auto tg_xx_zzz = pbuffer.data(idx_df + 9);

    auto tg_xy_xxy = pbuffer.data(idx_df + 11);

    auto tg_xy_xyy = pbuffer.data(idx_df + 13);

    auto tg_xy_yyy = pbuffer.data(idx_df + 16);

    auto tg_xy_yyz = pbuffer.data(idx_df + 17);

    auto tg_xy_yzz = pbuffer.data(idx_df + 18);

    auto tg_xz_xxx = pbuffer.data(idx_df + 20);

    auto tg_xz_xxz = pbuffer.data(idx_df + 22);

    auto tg_xz_xzz = pbuffer.data(idx_df + 25);

    auto tg_xz_yyz = pbuffer.data(idx_df + 27);

    auto tg_xz_yzz = pbuffer.data(idx_df + 28);

    auto tg_xz_zzz = pbuffer.data(idx_df + 29);

    auto tg_yy_xxx = pbuffer.data(idx_df + 30);

    auto tg_yy_xxy = pbuffer.data(idx_df + 31);

    auto tg_yy_xxz = pbuffer.data(idx_df + 32);

    auto tg_yy_xyy = pbuffer.data(idx_df + 33);

    auto tg_yy_xyz = pbuffer.data(idx_df + 34);

    auto tg_yy_xzz = pbuffer.data(idx_df + 35);

    auto tg_yy_yyy = pbuffer.data(idx_df + 36);

    auto tg_yy_yyz = pbuffer.data(idx_df + 37);

    auto tg_yy_yzz = pbuffer.data(idx_df + 38);

    auto tg_yy_zzz = pbuffer.data(idx_df + 39);

    auto tg_yz_xxz = pbuffer.data(idx_df + 42);

    auto tg_yz_xyz = pbuffer.data(idx_df + 44);

    auto tg_yz_xzz = pbuffer.data(idx_df + 45);

    auto tg_yz_yyy = pbuffer.data(idx_df + 46);

    auto tg_yz_yyz = pbuffer.data(idx_df + 47);

    auto tg_yz_yzz = pbuffer.data(idx_df + 48);

    auto tg_yz_zzz = pbuffer.data(idx_df + 49);

    auto tg_zz_xxx = pbuffer.data(idx_df + 50);

    auto tg_zz_xxy = pbuffer.data(idx_df + 51);

    auto tg_zz_xxz = pbuffer.data(idx_df + 52);

    auto tg_zz_xyy = pbuffer.data(idx_df + 53);

    auto tg_zz_xyz = pbuffer.data(idx_df + 54);

    auto tg_zz_xzz = pbuffer.data(idx_df + 55);

    auto tg_zz_yyy = pbuffer.data(idx_df + 56);

    auto tg_zz_yyz = pbuffer.data(idx_df + 57);

    auto tg_zz_yzz = pbuffer.data(idx_df + 58);

    auto tg_zz_zzz = pbuffer.data(idx_df + 59);

    // Set up components of targeted buffer : FF

    auto tg_xxx_xxx = pbuffer.data(idx_ff);

    auto tg_xxx_xxy = pbuffer.data(idx_ff + 1);

    auto tg_xxx_xxz = pbuffer.data(idx_ff + 2);

    auto tg_xxx_xyy = pbuffer.data(idx_ff + 3);

    auto tg_xxx_xyz = pbuffer.data(idx_ff + 4);

    auto tg_xxx_xzz = pbuffer.data(idx_ff + 5);

    auto tg_xxx_yyy = pbuffer.data(idx_ff + 6);

    auto tg_xxx_yyz = pbuffer.data(idx_ff + 7);

    auto tg_xxx_yzz = pbuffer.data(idx_ff + 8);

    auto tg_xxx_zzz = pbuffer.data(idx_ff + 9);

    auto tg_xxy_xxx = pbuffer.data(idx_ff + 10);

    auto tg_xxy_xxy = pbuffer.data(idx_ff + 11);

    auto tg_xxy_xxz = pbuffer.data(idx_ff + 12);

    auto tg_xxy_xyy = pbuffer.data(idx_ff + 13);

    auto tg_xxy_xyz = pbuffer.data(idx_ff + 14);

    auto tg_xxy_xzz = pbuffer.data(idx_ff + 15);

    auto tg_xxy_yyy = pbuffer.data(idx_ff + 16);

    auto tg_xxy_yyz = pbuffer.data(idx_ff + 17);

    auto tg_xxy_yzz = pbuffer.data(idx_ff + 18);

    auto tg_xxy_zzz = pbuffer.data(idx_ff + 19);

    auto tg_xxz_xxx = pbuffer.data(idx_ff + 20);

    auto tg_xxz_xxy = pbuffer.data(idx_ff + 21);

    auto tg_xxz_xxz = pbuffer.data(idx_ff + 22);

    auto tg_xxz_xyy = pbuffer.data(idx_ff + 23);

    auto tg_xxz_xyz = pbuffer.data(idx_ff + 24);

    auto tg_xxz_xzz = pbuffer.data(idx_ff + 25);

    auto tg_xxz_yyy = pbuffer.data(idx_ff + 26);

    auto tg_xxz_yyz = pbuffer.data(idx_ff + 27);

    auto tg_xxz_yzz = pbuffer.data(idx_ff + 28);

    auto tg_xxz_zzz = pbuffer.data(idx_ff + 29);

    auto tg_xyy_xxx = pbuffer.data(idx_ff + 30);

    auto tg_xyy_xxy = pbuffer.data(idx_ff + 31);

    auto tg_xyy_xxz = pbuffer.data(idx_ff + 32);

    auto tg_xyy_xyy = pbuffer.data(idx_ff + 33);

    auto tg_xyy_xyz = pbuffer.data(idx_ff + 34);

    auto tg_xyy_xzz = pbuffer.data(idx_ff + 35);

    auto tg_xyy_yyy = pbuffer.data(idx_ff + 36);

    auto tg_xyy_yyz = pbuffer.data(idx_ff + 37);

    auto tg_xyy_yzz = pbuffer.data(idx_ff + 38);

    auto tg_xyy_zzz = pbuffer.data(idx_ff + 39);

    auto tg_xyz_xxx = pbuffer.data(idx_ff + 40);

    auto tg_xyz_xxy = pbuffer.data(idx_ff + 41);

    auto tg_xyz_xxz = pbuffer.data(idx_ff + 42);

    auto tg_xyz_xyy = pbuffer.data(idx_ff + 43);

    auto tg_xyz_xyz = pbuffer.data(idx_ff + 44);

    auto tg_xyz_xzz = pbuffer.data(idx_ff + 45);

    auto tg_xyz_yyy = pbuffer.data(idx_ff + 46);

    auto tg_xyz_yyz = pbuffer.data(idx_ff + 47);

    auto tg_xyz_yzz = pbuffer.data(idx_ff + 48);

    auto tg_xyz_zzz = pbuffer.data(idx_ff + 49);

    auto tg_xzz_xxx = pbuffer.data(idx_ff + 50);

    auto tg_xzz_xxy = pbuffer.data(idx_ff + 51);

    auto tg_xzz_xxz = pbuffer.data(idx_ff + 52);

    auto tg_xzz_xyy = pbuffer.data(idx_ff + 53);

    auto tg_xzz_xyz = pbuffer.data(idx_ff + 54);

    auto tg_xzz_xzz = pbuffer.data(idx_ff + 55);

    auto tg_xzz_yyy = pbuffer.data(idx_ff + 56);

    auto tg_xzz_yyz = pbuffer.data(idx_ff + 57);

    auto tg_xzz_yzz = pbuffer.data(idx_ff + 58);

    auto tg_xzz_zzz = pbuffer.data(idx_ff + 59);

    auto tg_yyy_xxx = pbuffer.data(idx_ff + 60);

    auto tg_yyy_xxy = pbuffer.data(idx_ff + 61);

    auto tg_yyy_xxz = pbuffer.data(idx_ff + 62);

    auto tg_yyy_xyy = pbuffer.data(idx_ff + 63);

    auto tg_yyy_xyz = pbuffer.data(idx_ff + 64);

    auto tg_yyy_xzz = pbuffer.data(idx_ff + 65);

    auto tg_yyy_yyy = pbuffer.data(idx_ff + 66);

    auto tg_yyy_yyz = pbuffer.data(idx_ff + 67);

    auto tg_yyy_yzz = pbuffer.data(idx_ff + 68);

    auto tg_yyy_zzz = pbuffer.data(idx_ff + 69);

    auto tg_yyz_xxx = pbuffer.data(idx_ff + 70);

    auto tg_yyz_xxy = pbuffer.data(idx_ff + 71);

    auto tg_yyz_xxz = pbuffer.data(idx_ff + 72);

    auto tg_yyz_xyy = pbuffer.data(idx_ff + 73);

    auto tg_yyz_xyz = pbuffer.data(idx_ff + 74);

    auto tg_yyz_xzz = pbuffer.data(idx_ff + 75);

    auto tg_yyz_yyy = pbuffer.data(idx_ff + 76);

    auto tg_yyz_yyz = pbuffer.data(idx_ff + 77);

    auto tg_yyz_yzz = pbuffer.data(idx_ff + 78);

    auto tg_yyz_zzz = pbuffer.data(idx_ff + 79);

    auto tg_yzz_xxx = pbuffer.data(idx_ff + 80);

    auto tg_yzz_xxy = pbuffer.data(idx_ff + 81);

    auto tg_yzz_xxz = pbuffer.data(idx_ff + 82);

    auto tg_yzz_xyy = pbuffer.data(idx_ff + 83);

    auto tg_yzz_xyz = pbuffer.data(idx_ff + 84);

    auto tg_yzz_xzz = pbuffer.data(idx_ff + 85);

    auto tg_yzz_yyy = pbuffer.data(idx_ff + 86);

    auto tg_yzz_yyz = pbuffer.data(idx_ff + 87);

    auto tg_yzz_yzz = pbuffer.data(idx_ff + 88);

    auto tg_yzz_zzz = pbuffer.data(idx_ff + 89);

    auto tg_zzz_xxx = pbuffer.data(idx_ff + 90);

    auto tg_zzz_xxy = pbuffer.data(idx_ff + 91);

    auto tg_zzz_xxz = pbuffer.data(idx_ff + 92);

    auto tg_zzz_xyy = pbuffer.data(idx_ff + 93);

    auto tg_zzz_xyz = pbuffer.data(idx_ff + 94);

    auto tg_zzz_xzz = pbuffer.data(idx_ff + 95);

    auto tg_zzz_yyy = pbuffer.data(idx_ff + 96);

    auto tg_zzz_yyz = pbuffer.data(idx_ff + 97);

    auto tg_zzz_yzz = pbuffer.data(idx_ff + 98);

    auto tg_zzz_zzz = pbuffer.data(idx_ff + 99);

    #pragma omp simd aligned(fxi, ra_x, ra_y, ra_z, tg_x_xxx, tg_x_xxy, tg_x_xxz, tg_x_xyy, tg_x_xyz, tg_x_xzz, tg_x_yyy, tg_x_yyz, tg_x_yzz, tg_x_zzz, tg_xx_xx, tg_xx_xxx, tg_xx_xxy, tg_xx_xxz, tg_xx_xy, tg_xx_xyy, tg_xx_xyz, tg_xx_xz, tg_xx_xzz, tg_xx_yy, tg_xx_yyy, tg_xx_yyz, tg_xx_yz, tg_xx_yzz, tg_xx_zz, tg_xx_zzz, tg_xxx_xxx, tg_xxx_xxy, tg_xxx_xxz, tg_xxx_xyy, tg_xxx_xyz, tg_xxx_xzz, tg_xxx_yyy, tg_xxx_yyz, tg_xxx_yzz, tg_xxx_zzz, tg_xxy_xxx, tg_xxy_xxy, tg_xxy_xxz, tg_xxy_xyy, tg_xxy_xyz, tg_xxy_xzz, tg_xxy_yyy, tg_xxy_yyz, tg_xxy_yzz, tg_xxy_zzz, tg_xxz_xxx, tg_xxz_xxy, tg_xxz_xxz, tg_xxz_xyy, tg_xxz_xyz, tg_xxz_xzz, tg_xxz_yyy, tg_xxz_yyz, tg_xxz_yzz, tg_xxz_zzz, tg_xy_xxy, tg_xy_xyy, tg_xy_yyy, tg_xy_yyz, tg_xy_yzz, tg_xyy_xxx, tg_xyy_xxy, tg_xyy_xxz, tg_xyy_xyy, tg_xyy_xyz, tg_xyy_xzz, tg_xyy_yyy, tg_xyy_yyz, tg_xyy_yzz, tg_xyy_zzz, tg_xyz_xxx, tg_xyz_xxy, tg_xyz_xxz, tg_xyz_xyy, tg_xyz_xyz, tg_xyz_xzz, tg_xyz_yyy, tg_xyz_yyz, tg_xyz_yzz, tg_xyz_zzz, tg_xz_xxx, tg_xz_xxz, tg_xz_xzz, tg_xz_yyz, tg_xz_yzz, tg_xz_zzz, tg_xzz_xxx, tg_xzz_xxy, tg_xzz_xxz, tg_xzz_xyy, tg_xzz_xyz, tg_xzz_xzz, tg_xzz_yyy, tg_xzz_yyz, tg_xzz_yzz, tg_xzz_zzz, tg_y_xxx, tg_y_xxy, tg_y_xxz, tg_y_xyy, tg_y_xyz, tg_y_xzz, tg_y_yyy, tg_y_yyz, tg_y_yzz, tg_y_zzz, tg_yy_xx, tg_yy_xxx, tg_yy_xxy, tg_yy_xxz, tg_yy_xy, tg_yy_xyy, tg_yy_xyz, tg_yy_xz, tg_yy_xzz, tg_yy_yy, tg_yy_yyy, tg_yy_yyz, tg_yy_yz, tg_yy_yzz, tg_yy_zz, tg_yy_zzz, tg_yyy_xxx, tg_yyy_xxy, tg_yyy_xxz, tg_yyy_xyy, tg_yyy_xyz, tg_yyy_xzz, tg_yyy_yyy, tg_yyy_yyz, tg_yyy_yzz, tg_yyy_zzz, tg_yyz_xxx, tg_yyz_xxy, tg_yyz_xxz, tg_yyz_xyy, tg_yyz_xyz, tg_yyz_xzz, tg_yyz_yyy, tg_yyz_yyz, tg_yyz_yzz, tg_yyz_zzz, tg_yz_xxz, tg_yz_xyz, tg_yz_xzz, tg_yz_yyy, tg_yz_yyz, tg_yz_yz, tg_yz_yzz, tg_yz_zzz, tg_yzz_xxx, tg_yzz_xxy, tg_yzz_xxz, tg_yzz_xyy, tg_yzz_xyz, tg_yzz_xzz, tg_yzz_yyy, tg_yzz_yyz, tg_yzz_yzz, tg_yzz_zzz, tg_z_xxx, tg_z_xxy, tg_z_xxz, tg_z_xyy, tg_z_xyz, tg_z_xzz, tg_z_yyy, tg_z_yyz, tg_z_yzz, tg_z_zzz, tg_zz_xx, tg_zz_xxx, tg_zz_xxy, tg_zz_xxz, tg_zz_xy, tg_zz_xyy, tg_zz_xyz, tg_zz_xz, tg_zz_xzz, tg_zz_yy, tg_zz_yyy, tg_zz_yyz, tg_zz_yz, tg_zz_yzz, tg_zz_zz, tg_zz_zzz, tg_zzz_xxx, tg_zzz_xxy, tg_zzz_xxz, tg_zzz_xyy, tg_zzz_xyz, tg_zzz_xzz, tg_zzz_yyy, tg_zzz_yyz, tg_zzz_yzz, tg_zzz_zzz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_xxx_xxx[i] = 2.0 * tg_x_xxx[i] * fxi[i] + 3.0 * tg_xx_xx[i] * fxi[i] + tg_xx_xxx[i] * ra_x[i];

        tg_xxx_xxy[i] = 2.0 * tg_x_xxy[i] * fxi[i] + 2.0 * tg_xx_xy[i] * fxi[i] + tg_xx_xxy[i] * ra_x[i];

        tg_xxx_xxz[i] = 2.0 * tg_x_xxz[i] * fxi[i] + 2.0 * tg_xx_xz[i] * fxi[i] + tg_xx_xxz[i] * ra_x[i];

        tg_xxx_xyy[i] = 2.0 * tg_x_xyy[i] * fxi[i] + tg_xx_yy[i] * fxi[i] + tg_xx_xyy[i] * ra_x[i];

        tg_xxx_xyz[i] = 2.0 * tg_x_xyz[i] * fxi[i] + tg_xx_yz[i] * fxi[i] + tg_xx_xyz[i] * ra_x[i];

        tg_xxx_xzz[i] = 2.0 * tg_x_xzz[i] * fxi[i] + tg_xx_zz[i] * fxi[i] + tg_xx_xzz[i] * ra_x[i];

        tg_xxx_yyy[i] = 2.0 * tg_x_yyy[i] * fxi[i] + tg_xx_yyy[i] * ra_x[i];

        tg_xxx_yyz[i] = 2.0 * tg_x_yyz[i] * fxi[i] + tg_xx_yyz[i] * ra_x[i];

        tg_xxx_yzz[i] = 2.0 * tg_x_yzz[i] * fxi[i] + tg_xx_yzz[i] * ra_x[i];

        tg_xxx_zzz[i] = 2.0 * tg_x_zzz[i] * fxi[i] + tg_xx_zzz[i] * ra_x[i];

        tg_xxy_xxx[i] = tg_xx_xxx[i] * ra_y[i];

        tg_xxy_xxy[i] = tg_xx_xx[i] * fxi[i] + tg_xx_xxy[i] * ra_y[i];

        tg_xxy_xxz[i] = tg_xx_xxz[i] * ra_y[i];

        tg_xxy_xyy[i] = 2.0 * tg_xx_xy[i] * fxi[i] + tg_xx_xyy[i] * ra_y[i];

        tg_xxy_xyz[i] = tg_xx_xz[i] * fxi[i] + tg_xx_xyz[i] * ra_y[i];

        tg_xxy_xzz[i] = tg_xx_xzz[i] * ra_y[i];

        tg_xxy_yyy[i] = tg_y_yyy[i] * fxi[i] + tg_xy_yyy[i] * ra_x[i];

        tg_xxy_yyz[i] = tg_y_yyz[i] * fxi[i] + tg_xy_yyz[i] * ra_x[i];

        tg_xxy_yzz[i] = tg_y_yzz[i] * fxi[i] + tg_xy_yzz[i] * ra_x[i];

        tg_xxy_zzz[i] = tg_xx_zzz[i] * ra_y[i];

        tg_xxz_xxx[i] = tg_xx_xxx[i] * ra_z[i];

        tg_xxz_xxy[i] = tg_xx_xxy[i] * ra_z[i];

        tg_xxz_xxz[i] = tg_xx_xx[i] * fxi[i] + tg_xx_xxz[i] * ra_z[i];

        tg_xxz_xyy[i] = tg_xx_xyy[i] * ra_z[i];

        tg_xxz_xyz[i] = tg_xx_xy[i] * fxi[i] + tg_xx_xyz[i] * ra_z[i];

        tg_xxz_xzz[i] = 2.0 * tg_xx_xz[i] * fxi[i] + tg_xx_xzz[i] * ra_z[i];

        tg_xxz_yyy[i] = tg_xx_yyy[i] * ra_z[i];

        tg_xxz_yyz[i] = tg_z_yyz[i] * fxi[i] + tg_xz_yyz[i] * ra_x[i];

        tg_xxz_yzz[i] = tg_z_yzz[i] * fxi[i] + tg_xz_yzz[i] * ra_x[i];

        tg_xxz_zzz[i] = tg_z_zzz[i] * fxi[i] + tg_xz_zzz[i] * ra_x[i];

        tg_xyy_xxx[i] = 3.0 * tg_yy_xx[i] * fxi[i] + tg_yy_xxx[i] * ra_x[i];

        tg_xyy_xxy[i] = 2.0 * tg_yy_xy[i] * fxi[i] + tg_yy_xxy[i] * ra_x[i];

        tg_xyy_xxz[i] = 2.0 * tg_yy_xz[i] * fxi[i] + tg_yy_xxz[i] * ra_x[i];

        tg_xyy_xyy[i] = tg_yy_yy[i] * fxi[i] + tg_yy_xyy[i] * ra_x[i];

        tg_xyy_xyz[i] = tg_yy_yz[i] * fxi[i] + tg_yy_xyz[i] * ra_x[i];

        tg_xyy_xzz[i] = tg_yy_zz[i] * fxi[i] + tg_yy_xzz[i] * ra_x[i];

        tg_xyy_yyy[i] = tg_yy_yyy[i] * ra_x[i];

        tg_xyy_yyz[i] = tg_yy_yyz[i] * ra_x[i];

        tg_xyy_yzz[i] = tg_yy_yzz[i] * ra_x[i];

        tg_xyy_zzz[i] = tg_yy_zzz[i] * ra_x[i];

        tg_xyz_xxx[i] = tg_xz_xxx[i] * ra_y[i];

        tg_xyz_xxy[i] = tg_xy_xxy[i] * ra_z[i];

        tg_xyz_xxz[i] = tg_xz_xxz[i] * ra_y[i];

        tg_xyz_xyy[i] = tg_xy_xyy[i] * ra_z[i];

        tg_xyz_xyz[i] = tg_yz_yz[i] * fxi[i] + tg_yz_xyz[i] * ra_x[i];

        tg_xyz_xzz[i] = tg_xz_xzz[i] * ra_y[i];

        tg_xyz_yyy[i] = tg_yz_yyy[i] * ra_x[i];

        tg_xyz_yyz[i] = tg_yz_yyz[i] * ra_x[i];

        tg_xyz_yzz[i] = tg_yz_yzz[i] * ra_x[i];

        tg_xyz_zzz[i] = tg_yz_zzz[i] * ra_x[i];

        tg_xzz_xxx[i] = 3.0 * tg_zz_xx[i] * fxi[i] + tg_zz_xxx[i] * ra_x[i];

        tg_xzz_xxy[i] = 2.0 * tg_zz_xy[i] * fxi[i] + tg_zz_xxy[i] * ra_x[i];

        tg_xzz_xxz[i] = 2.0 * tg_zz_xz[i] * fxi[i] + tg_zz_xxz[i] * ra_x[i];

        tg_xzz_xyy[i] = tg_zz_yy[i] * fxi[i] + tg_zz_xyy[i] * ra_x[i];

        tg_xzz_xyz[i] = tg_zz_yz[i] * fxi[i] + tg_zz_xyz[i] * ra_x[i];

        tg_xzz_xzz[i] = tg_zz_zz[i] * fxi[i] + tg_zz_xzz[i] * ra_x[i];

        tg_xzz_yyy[i] = tg_zz_yyy[i] * ra_x[i];

        tg_xzz_yyz[i] = tg_zz_yyz[i] * ra_x[i];

        tg_xzz_yzz[i] = tg_zz_yzz[i] * ra_x[i];

        tg_xzz_zzz[i] = tg_zz_zzz[i] * ra_x[i];

        tg_yyy_xxx[i] = 2.0 * tg_y_xxx[i] * fxi[i] + tg_yy_xxx[i] * ra_y[i];

        tg_yyy_xxy[i] = 2.0 * tg_y_xxy[i] * fxi[i] + tg_yy_xx[i] * fxi[i] + tg_yy_xxy[i] * ra_y[i];

        tg_yyy_xxz[i] = 2.0 * tg_y_xxz[i] * fxi[i] + tg_yy_xxz[i] * ra_y[i];

        tg_yyy_xyy[i] = 2.0 * tg_y_xyy[i] * fxi[i] + 2.0 * tg_yy_xy[i] * fxi[i] + tg_yy_xyy[i] * ra_y[i];

        tg_yyy_xyz[i] = 2.0 * tg_y_xyz[i] * fxi[i] + tg_yy_xz[i] * fxi[i] + tg_yy_xyz[i] * ra_y[i];

        tg_yyy_xzz[i] = 2.0 * tg_y_xzz[i] * fxi[i] + tg_yy_xzz[i] * ra_y[i];

        tg_yyy_yyy[i] = 2.0 * tg_y_yyy[i] * fxi[i] + 3.0 * tg_yy_yy[i] * fxi[i] + tg_yy_yyy[i] * ra_y[i];

        tg_yyy_yyz[i] = 2.0 * tg_y_yyz[i] * fxi[i] + 2.0 * tg_yy_yz[i] * fxi[i] + tg_yy_yyz[i] * ra_y[i];

        tg_yyy_yzz[i] = 2.0 * tg_y_yzz[i] * fxi[i] + tg_yy_zz[i] * fxi[i] + tg_yy_yzz[i] * ra_y[i];

        tg_yyy_zzz[i] = 2.0 * tg_y_zzz[i] * fxi[i] + tg_yy_zzz[i] * ra_y[i];

        tg_yyz_xxx[i] = tg_yy_xxx[i] * ra_z[i];

        tg_yyz_xxy[i] = tg_yy_xxy[i] * ra_z[i];

        tg_yyz_xxz[i] = tg_z_xxz[i] * fxi[i] + tg_yz_xxz[i] * ra_y[i];

        tg_yyz_xyy[i] = tg_yy_xyy[i] * ra_z[i];

        tg_yyz_xyz[i] = tg_yy_xy[i] * fxi[i] + tg_yy_xyz[i] * ra_z[i];

        tg_yyz_xzz[i] = tg_z_xzz[i] * fxi[i] + tg_yz_xzz[i] * ra_y[i];

        tg_yyz_yyy[i] = tg_yy_yyy[i] * ra_z[i];

        tg_yyz_yyz[i] = tg_yy_yy[i] * fxi[i] + tg_yy_yyz[i] * ra_z[i];

        tg_yyz_yzz[i] = 2.0 * tg_yy_yz[i] * fxi[i] + tg_yy_yzz[i] * ra_z[i];

        tg_yyz_zzz[i] = tg_z_zzz[i] * fxi[i] + tg_yz_zzz[i] * ra_y[i];

        tg_yzz_xxx[i] = tg_zz_xxx[i] * ra_y[i];

        tg_yzz_xxy[i] = tg_zz_xx[i] * fxi[i] + tg_zz_xxy[i] * ra_y[i];

        tg_yzz_xxz[i] = tg_zz_xxz[i] * ra_y[i];

        tg_yzz_xyy[i] = 2.0 * tg_zz_xy[i] * fxi[i] + tg_zz_xyy[i] * ra_y[i];

        tg_yzz_xyz[i] = tg_zz_xz[i] * fxi[i] + tg_zz_xyz[i] * ra_y[i];

        tg_yzz_xzz[i] = tg_zz_xzz[i] * ra_y[i];

        tg_yzz_yyy[i] = 3.0 * tg_zz_yy[i] * fxi[i] + tg_zz_yyy[i] * ra_y[i];

        tg_yzz_yyz[i] = 2.0 * tg_zz_yz[i] * fxi[i] + tg_zz_yyz[i] * ra_y[i];

        tg_yzz_yzz[i] = tg_zz_zz[i] * fxi[i] + tg_zz_yzz[i] * ra_y[i];

        tg_yzz_zzz[i] = tg_zz_zzz[i] * ra_y[i];

        tg_zzz_xxx[i] = 2.0 * tg_z_xxx[i] * fxi[i] + tg_zz_xxx[i] * ra_z[i];

        tg_zzz_xxy[i] = 2.0 * tg_z_xxy[i] * fxi[i] + tg_zz_xxy[i] * ra_z[i];

        tg_zzz_xxz[i] = 2.0 * tg_z_xxz[i] * fxi[i] + tg_zz_xx[i] * fxi[i] + tg_zz_xxz[i] * ra_z[i];

        tg_zzz_xyy[i] = 2.0 * tg_z_xyy[i] * fxi[i] + tg_zz_xyy[i] * ra_z[i];

        tg_zzz_xyz[i] = 2.0 * tg_z_xyz[i] * fxi[i] + tg_zz_xy[i] * fxi[i] + tg_zz_xyz[i] * ra_z[i];

        tg_zzz_xzz[i] = 2.0 * tg_z_xzz[i] * fxi[i] + 2.0 * tg_zz_xz[i] * fxi[i] + tg_zz_xzz[i] * ra_z[i];

        tg_zzz_yyy[i] = 2.0 * tg_z_yyy[i] * fxi[i] + tg_zz_yyy[i] * ra_z[i];

        tg_zzz_yyz[i] = 2.0 * tg_z_yyz[i] * fxi[i] + tg_zz_yy[i] * fxi[i] + tg_zz_yyz[i] * ra_z[i];

        tg_zzz_yzz[i] = 2.0 * tg_z_yzz[i] * fxi[i] + 2.0 * tg_zz_yz[i] * fxi[i] + tg_zz_yzz[i] * ra_z[i];

        tg_zzz_zzz[i] = 2.0 * tg_z_zzz[i] * fxi[i] + 3.0 * tg_zz_zz[i] * fxi[i] + tg_zz_zzz[i] * ra_z[i];
    }
}

} // t2lecp namespace


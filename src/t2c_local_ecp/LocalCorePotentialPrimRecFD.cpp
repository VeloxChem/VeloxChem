#include "LocalCorePotentialPrimRecFD.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_fd(CSimdArray<double>& pbuffer, 
                                  const size_t idx_fd,
                                  const size_t idx_pd,
                                  const size_t idx_dp,
                                  const size_t idx_dd,
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

    // Set up components of auxiliary buffer : PD

    auto tg_x_xx = pbuffer.data(idx_pd);

    auto tg_x_xy = pbuffer.data(idx_pd + 1);

    auto tg_x_xz = pbuffer.data(idx_pd + 2);

    auto tg_x_yy = pbuffer.data(idx_pd + 3);

    auto tg_x_yz = pbuffer.data(idx_pd + 4);

    auto tg_x_zz = pbuffer.data(idx_pd + 5);

    auto tg_y_xx = pbuffer.data(idx_pd + 6);

    auto tg_y_xy = pbuffer.data(idx_pd + 7);

    auto tg_y_xz = pbuffer.data(idx_pd + 8);

    auto tg_y_yy = pbuffer.data(idx_pd + 9);

    auto tg_y_yz = pbuffer.data(idx_pd + 10);

    auto tg_y_zz = pbuffer.data(idx_pd + 11);

    auto tg_z_xx = pbuffer.data(idx_pd + 12);

    auto tg_z_xy = pbuffer.data(idx_pd + 13);

    auto tg_z_xz = pbuffer.data(idx_pd + 14);

    auto tg_z_yy = pbuffer.data(idx_pd + 15);

    auto tg_z_yz = pbuffer.data(idx_pd + 16);

    auto tg_z_zz = pbuffer.data(idx_pd + 17);

    // Set up components of auxiliary buffer : DP

    auto tg_xx_x = pbuffer.data(idx_dp);

    auto tg_xx_y = pbuffer.data(idx_dp + 1);

    auto tg_xx_z = pbuffer.data(idx_dp + 2);

    auto tg_yy_x = pbuffer.data(idx_dp + 9);

    auto tg_yy_y = pbuffer.data(idx_dp + 10);

    auto tg_yy_z = pbuffer.data(idx_dp + 11);

    auto tg_zz_x = pbuffer.data(idx_dp + 15);

    auto tg_zz_y = pbuffer.data(idx_dp + 16);

    auto tg_zz_z = pbuffer.data(idx_dp + 17);

    // Set up components of auxiliary buffer : DD

    auto tg_xx_xx = pbuffer.data(idx_dd);

    auto tg_xx_xy = pbuffer.data(idx_dd + 1);

    auto tg_xx_xz = pbuffer.data(idx_dd + 2);

    auto tg_xx_yy = pbuffer.data(idx_dd + 3);

    auto tg_xx_yz = pbuffer.data(idx_dd + 4);

    auto tg_xx_zz = pbuffer.data(idx_dd + 5);

    auto tg_xy_xy = pbuffer.data(idx_dd + 7);

    auto tg_xy_yy = pbuffer.data(idx_dd + 9);

    auto tg_xy_yz = pbuffer.data(idx_dd + 10);

    auto tg_xz_xx = pbuffer.data(idx_dd + 12);

    auto tg_xz_xz = pbuffer.data(idx_dd + 14);

    auto tg_xz_yz = pbuffer.data(idx_dd + 16);

    auto tg_xz_zz = pbuffer.data(idx_dd + 17);

    auto tg_yy_xx = pbuffer.data(idx_dd + 18);

    auto tg_yy_xy = pbuffer.data(idx_dd + 19);

    auto tg_yy_xz = pbuffer.data(idx_dd + 20);

    auto tg_yy_yy = pbuffer.data(idx_dd + 21);

    auto tg_yy_yz = pbuffer.data(idx_dd + 22);

    auto tg_yy_zz = pbuffer.data(idx_dd + 23);

    auto tg_yz_xz = pbuffer.data(idx_dd + 26);

    auto tg_yz_yy = pbuffer.data(idx_dd + 27);

    auto tg_yz_yz = pbuffer.data(idx_dd + 28);

    auto tg_yz_zz = pbuffer.data(idx_dd + 29);

    auto tg_zz_xx = pbuffer.data(idx_dd + 30);

    auto tg_zz_xy = pbuffer.data(idx_dd + 31);

    auto tg_zz_xz = pbuffer.data(idx_dd + 32);

    auto tg_zz_yy = pbuffer.data(idx_dd + 33);

    auto tg_zz_yz = pbuffer.data(idx_dd + 34);

    auto tg_zz_zz = pbuffer.data(idx_dd + 35);

    // Set up components of targeted buffer : FD

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

    auto tg_xxy_zz = pbuffer.data(idx_fd + 11);

    auto tg_xxz_xx = pbuffer.data(idx_fd + 12);

    auto tg_xxz_xy = pbuffer.data(idx_fd + 13);

    auto tg_xxz_xz = pbuffer.data(idx_fd + 14);

    auto tg_xxz_yy = pbuffer.data(idx_fd + 15);

    auto tg_xxz_yz = pbuffer.data(idx_fd + 16);

    auto tg_xxz_zz = pbuffer.data(idx_fd + 17);

    auto tg_xyy_xx = pbuffer.data(idx_fd + 18);

    auto tg_xyy_xy = pbuffer.data(idx_fd + 19);

    auto tg_xyy_xz = pbuffer.data(idx_fd + 20);

    auto tg_xyy_yy = pbuffer.data(idx_fd + 21);

    auto tg_xyy_yz = pbuffer.data(idx_fd + 22);

    auto tg_xyy_zz = pbuffer.data(idx_fd + 23);

    auto tg_xyz_xx = pbuffer.data(idx_fd + 24);

    auto tg_xyz_xy = pbuffer.data(idx_fd + 25);

    auto tg_xyz_xz = pbuffer.data(idx_fd + 26);

    auto tg_xyz_yy = pbuffer.data(idx_fd + 27);

    auto tg_xyz_yz = pbuffer.data(idx_fd + 28);

    auto tg_xyz_zz = pbuffer.data(idx_fd + 29);

    auto tg_xzz_xx = pbuffer.data(idx_fd + 30);

    auto tg_xzz_xy = pbuffer.data(idx_fd + 31);

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

    auto tg_yyz_xx = pbuffer.data(idx_fd + 42);

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

    #pragma omp simd aligned(fxi, ra_x, ra_y, ra_z, tg_x_xx, tg_x_xy, tg_x_xz, tg_x_yy, tg_x_yz, tg_x_zz, tg_xx_x, tg_xx_xx, tg_xx_xy, tg_xx_xz, tg_xx_y, tg_xx_yy, tg_xx_yz, tg_xx_z, tg_xx_zz, tg_xxx_xx, tg_xxx_xy, tg_xxx_xz, tg_xxx_yy, tg_xxx_yz, tg_xxx_zz, tg_xxy_xx, tg_xxy_xy, tg_xxy_xz, tg_xxy_yy, tg_xxy_yz, tg_xxy_zz, tg_xxz_xx, tg_xxz_xy, tg_xxz_xz, tg_xxz_yy, tg_xxz_yz, tg_xxz_zz, tg_xy_xy, tg_xy_yy, tg_xy_yz, tg_xyy_xx, tg_xyy_xy, tg_xyy_xz, tg_xyy_yy, tg_xyy_yz, tg_xyy_zz, tg_xyz_xx, tg_xyz_xy, tg_xyz_xz, tg_xyz_yy, tg_xyz_yz, tg_xyz_zz, tg_xz_xx, tg_xz_xz, tg_xz_yz, tg_xz_zz, tg_xzz_xx, tg_xzz_xy, tg_xzz_xz, tg_xzz_yy, tg_xzz_yz, tg_xzz_zz, tg_y_xx, tg_y_xy, tg_y_xz, tg_y_yy, tg_y_yz, tg_y_zz, tg_yy_x, tg_yy_xx, tg_yy_xy, tg_yy_xz, tg_yy_y, tg_yy_yy, tg_yy_yz, tg_yy_z, tg_yy_zz, tg_yyy_xx, tg_yyy_xy, tg_yyy_xz, tg_yyy_yy, tg_yyy_yz, tg_yyy_zz, tg_yyz_xx, tg_yyz_xy, tg_yyz_xz, tg_yyz_yy, tg_yyz_yz, tg_yyz_zz, tg_yz_xz, tg_yz_yy, tg_yz_yz, tg_yz_zz, tg_yzz_xx, tg_yzz_xy, tg_yzz_xz, tg_yzz_yy, tg_yzz_yz, tg_yzz_zz, tg_z_xx, tg_z_xy, tg_z_xz, tg_z_yy, tg_z_yz, tg_z_zz, tg_zz_x, tg_zz_xx, tg_zz_xy, tg_zz_xz, tg_zz_y, tg_zz_yy, tg_zz_yz, tg_zz_z, tg_zz_zz, tg_zzz_xx, tg_zzz_xy, tg_zzz_xz, tg_zzz_yy, tg_zzz_yz, tg_zzz_zz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_xxx_xx[i] = 2.0 * tg_x_xx[i] * fxi[i] + 2.0 * tg_xx_x[i] * fxi[i] + tg_xx_xx[i] * ra_x[i];

        tg_xxx_xy[i] = 2.0 * tg_x_xy[i] * fxi[i] + tg_xx_y[i] * fxi[i] + tg_xx_xy[i] * ra_x[i];

        tg_xxx_xz[i] = 2.0 * tg_x_xz[i] * fxi[i] + tg_xx_z[i] * fxi[i] + tg_xx_xz[i] * ra_x[i];

        tg_xxx_yy[i] = 2.0 * tg_x_yy[i] * fxi[i] + tg_xx_yy[i] * ra_x[i];

        tg_xxx_yz[i] = 2.0 * tg_x_yz[i] * fxi[i] + tg_xx_yz[i] * ra_x[i];

        tg_xxx_zz[i] = 2.0 * tg_x_zz[i] * fxi[i] + tg_xx_zz[i] * ra_x[i];

        tg_xxy_xx[i] = tg_xx_xx[i] * ra_y[i];

        tg_xxy_xy[i] = tg_xx_x[i] * fxi[i] + tg_xx_xy[i] * ra_y[i];

        tg_xxy_xz[i] = tg_xx_xz[i] * ra_y[i];

        tg_xxy_yy[i] = tg_y_yy[i] * fxi[i] + tg_xy_yy[i] * ra_x[i];

        tg_xxy_yz[i] = tg_y_yz[i] * fxi[i] + tg_xy_yz[i] * ra_x[i];

        tg_xxy_zz[i] = tg_xx_zz[i] * ra_y[i];

        tg_xxz_xx[i] = tg_xx_xx[i] * ra_z[i];

        tg_xxz_xy[i] = tg_xx_xy[i] * ra_z[i];

        tg_xxz_xz[i] = tg_xx_x[i] * fxi[i] + tg_xx_xz[i] * ra_z[i];

        tg_xxz_yy[i] = tg_xx_yy[i] * ra_z[i];

        tg_xxz_yz[i] = tg_z_yz[i] * fxi[i] + tg_xz_yz[i] * ra_x[i];

        tg_xxz_zz[i] = tg_z_zz[i] * fxi[i] + tg_xz_zz[i] * ra_x[i];

        tg_xyy_xx[i] = 2.0 * tg_yy_x[i] * fxi[i] + tg_yy_xx[i] * ra_x[i];

        tg_xyy_xy[i] = tg_yy_y[i] * fxi[i] + tg_yy_xy[i] * ra_x[i];

        tg_xyy_xz[i] = tg_yy_z[i] * fxi[i] + tg_yy_xz[i] * ra_x[i];

        tg_xyy_yy[i] = tg_yy_yy[i] * ra_x[i];

        tg_xyy_yz[i] = tg_yy_yz[i] * ra_x[i];

        tg_xyy_zz[i] = tg_yy_zz[i] * ra_x[i];

        tg_xyz_xx[i] = tg_xz_xx[i] * ra_y[i];

        tg_xyz_xy[i] = tg_xy_xy[i] * ra_z[i];

        tg_xyz_xz[i] = tg_xz_xz[i] * ra_y[i];

        tg_xyz_yy[i] = tg_yz_yy[i] * ra_x[i];

        tg_xyz_yz[i] = tg_yz_yz[i] * ra_x[i];

        tg_xyz_zz[i] = tg_yz_zz[i] * ra_x[i];

        tg_xzz_xx[i] = 2.0 * tg_zz_x[i] * fxi[i] + tg_zz_xx[i] * ra_x[i];

        tg_xzz_xy[i] = tg_zz_y[i] * fxi[i] + tg_zz_xy[i] * ra_x[i];

        tg_xzz_xz[i] = tg_zz_z[i] * fxi[i] + tg_zz_xz[i] * ra_x[i];

        tg_xzz_yy[i] = tg_zz_yy[i] * ra_x[i];

        tg_xzz_yz[i] = tg_zz_yz[i] * ra_x[i];

        tg_xzz_zz[i] = tg_zz_zz[i] * ra_x[i];

        tg_yyy_xx[i] = 2.0 * tg_y_xx[i] * fxi[i] + tg_yy_xx[i] * ra_y[i];

        tg_yyy_xy[i] = 2.0 * tg_y_xy[i] * fxi[i] + tg_yy_x[i] * fxi[i] + tg_yy_xy[i] * ra_y[i];

        tg_yyy_xz[i] = 2.0 * tg_y_xz[i] * fxi[i] + tg_yy_xz[i] * ra_y[i];

        tg_yyy_yy[i] = 2.0 * tg_y_yy[i] * fxi[i] + 2.0 * tg_yy_y[i] * fxi[i] + tg_yy_yy[i] * ra_y[i];

        tg_yyy_yz[i] = 2.0 * tg_y_yz[i] * fxi[i] + tg_yy_z[i] * fxi[i] + tg_yy_yz[i] * ra_y[i];

        tg_yyy_zz[i] = 2.0 * tg_y_zz[i] * fxi[i] + tg_yy_zz[i] * ra_y[i];

        tg_yyz_xx[i] = tg_yy_xx[i] * ra_z[i];

        tg_yyz_xy[i] = tg_yy_xy[i] * ra_z[i];

        tg_yyz_xz[i] = tg_z_xz[i] * fxi[i] + tg_yz_xz[i] * ra_y[i];

        tg_yyz_yy[i] = tg_yy_yy[i] * ra_z[i];

        tg_yyz_yz[i] = tg_yy_y[i] * fxi[i] + tg_yy_yz[i] * ra_z[i];

        tg_yyz_zz[i] = tg_z_zz[i] * fxi[i] + tg_yz_zz[i] * ra_y[i];

        tg_yzz_xx[i] = tg_zz_xx[i] * ra_y[i];

        tg_yzz_xy[i] = tg_zz_x[i] * fxi[i] + tg_zz_xy[i] * ra_y[i];

        tg_yzz_xz[i] = tg_zz_xz[i] * ra_y[i];

        tg_yzz_yy[i] = 2.0 * tg_zz_y[i] * fxi[i] + tg_zz_yy[i] * ra_y[i];

        tg_yzz_yz[i] = tg_zz_z[i] * fxi[i] + tg_zz_yz[i] * ra_y[i];

        tg_yzz_zz[i] = tg_zz_zz[i] * ra_y[i];

        tg_zzz_xx[i] = 2.0 * tg_z_xx[i] * fxi[i] + tg_zz_xx[i] * ra_z[i];

        tg_zzz_xy[i] = 2.0 * tg_z_xy[i] * fxi[i] + tg_zz_xy[i] * ra_z[i];

        tg_zzz_xz[i] = 2.0 * tg_z_xz[i] * fxi[i] + tg_zz_x[i] * fxi[i] + tg_zz_xz[i] * ra_z[i];

        tg_zzz_yy[i] = 2.0 * tg_z_yy[i] * fxi[i] + tg_zz_yy[i] * ra_z[i];

        tg_zzz_yz[i] = 2.0 * tg_z_yz[i] * fxi[i] + tg_zz_y[i] * fxi[i] + tg_zz_yz[i] * ra_z[i];

        tg_zzz_zz[i] = 2.0 * tg_z_zz[i] * fxi[i] + 2.0 * tg_zz_z[i] * fxi[i] + tg_zz_zz[i] * ra_z[i];
    }
}

} // t2lecp namespace


#include "LocalCorePotentialPrimRecDF.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_df(CSimdArray<double>& pbuffer, 
                                  const size_t idx_df,
                                  const size_t idx_sf,
                                  const size_t idx_pd,
                                  const size_t idx_pf,
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

    // Set up components of auxiliary buffer : SF

    auto tg_0_xxx = pbuffer.data(idx_sf);

    auto tg_0_xxy = pbuffer.data(idx_sf + 1);

    auto tg_0_xxz = pbuffer.data(idx_sf + 2);

    auto tg_0_xyy = pbuffer.data(idx_sf + 3);

    auto tg_0_xyz = pbuffer.data(idx_sf + 4);

    auto tg_0_xzz = pbuffer.data(idx_sf + 5);

    auto tg_0_yyy = pbuffer.data(idx_sf + 6);

    auto tg_0_yyz = pbuffer.data(idx_sf + 7);

    auto tg_0_yzz = pbuffer.data(idx_sf + 8);

    auto tg_0_zzz = pbuffer.data(idx_sf + 9);

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

    // Set up components of targeted buffer : DF

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

    auto tg_xy_xxx = pbuffer.data(idx_df + 10);

    auto tg_xy_xxy = pbuffer.data(idx_df + 11);

    auto tg_xy_xxz = pbuffer.data(idx_df + 12);

    auto tg_xy_xyy = pbuffer.data(idx_df + 13);

    auto tg_xy_xyz = pbuffer.data(idx_df + 14);

    auto tg_xy_xzz = pbuffer.data(idx_df + 15);

    auto tg_xy_yyy = pbuffer.data(idx_df + 16);

    auto tg_xy_yyz = pbuffer.data(idx_df + 17);

    auto tg_xy_yzz = pbuffer.data(idx_df + 18);

    auto tg_xy_zzz = pbuffer.data(idx_df + 19);

    auto tg_xz_xxx = pbuffer.data(idx_df + 20);

    auto tg_xz_xxy = pbuffer.data(idx_df + 21);

    auto tg_xz_xxz = pbuffer.data(idx_df + 22);

    auto tg_xz_xyy = pbuffer.data(idx_df + 23);

    auto tg_xz_xyz = pbuffer.data(idx_df + 24);

    auto tg_xz_xzz = pbuffer.data(idx_df + 25);

    auto tg_xz_yyy = pbuffer.data(idx_df + 26);

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

    auto tg_yz_xxx = pbuffer.data(idx_df + 40);

    auto tg_yz_xxy = pbuffer.data(idx_df + 41);

    auto tg_yz_xxz = pbuffer.data(idx_df + 42);

    auto tg_yz_xyy = pbuffer.data(idx_df + 43);

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

    #pragma omp simd aligned(fxi, ra_x, ra_y, ra_z, tg_0_xxx, tg_0_xxy, tg_0_xxz, tg_0_xyy, tg_0_xyz, tg_0_xzz, tg_0_yyy, tg_0_yyz, tg_0_yzz, tg_0_zzz, tg_x_xx, tg_x_xxx, tg_x_xxy, tg_x_xxz, tg_x_xy, tg_x_xyy, tg_x_xyz, tg_x_xz, tg_x_xzz, tg_x_yy, tg_x_yyy, tg_x_yyz, tg_x_yz, tg_x_yzz, tg_x_zz, tg_x_zzz, tg_xx_xxx, tg_xx_xxy, tg_xx_xxz, tg_xx_xyy, tg_xx_xyz, tg_xx_xzz, tg_xx_yyy, tg_xx_yyz, tg_xx_yzz, tg_xx_zzz, tg_xy_xxx, tg_xy_xxy, tg_xy_xxz, tg_xy_xyy, tg_xy_xyz, tg_xy_xzz, tg_xy_yyy, tg_xy_yyz, tg_xy_yzz, tg_xy_zzz, tg_xz_xxx, tg_xz_xxy, tg_xz_xxz, tg_xz_xyy, tg_xz_xyz, tg_xz_xzz, tg_xz_yyy, tg_xz_yyz, tg_xz_yzz, tg_xz_zzz, tg_y_xx, tg_y_xxx, tg_y_xxy, tg_y_xxz, tg_y_xy, tg_y_xyy, tg_y_xyz, tg_y_xz, tg_y_xzz, tg_y_yy, tg_y_yyy, tg_y_yyz, tg_y_yz, tg_y_yzz, tg_y_zz, tg_y_zzz, tg_yy_xxx, tg_yy_xxy, tg_yy_xxz, tg_yy_xyy, tg_yy_xyz, tg_yy_xzz, tg_yy_yyy, tg_yy_yyz, tg_yy_yzz, tg_yy_zzz, tg_yz_xxx, tg_yz_xxy, tg_yz_xxz, tg_yz_xyy, tg_yz_xyz, tg_yz_xzz, tg_yz_yyy, tg_yz_yyz, tg_yz_yzz, tg_yz_zzz, tg_z_xx, tg_z_xxx, tg_z_xxy, tg_z_xxz, tg_z_xy, tg_z_xyy, tg_z_xyz, tg_z_xz, tg_z_xzz, tg_z_yy, tg_z_yyy, tg_z_yyz, tg_z_yz, tg_z_yzz, tg_z_zz, tg_z_zzz, tg_zz_xxx, tg_zz_xxy, tg_zz_xxz, tg_zz_xyy, tg_zz_xyz, tg_zz_xzz, tg_zz_yyy, tg_zz_yyz, tg_zz_yzz, tg_zz_zzz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_xx_xxx[i] = tg_0_xxx[i] * fxi[i] + 3.0 * tg_x_xx[i] * fxi[i] + tg_x_xxx[i] * ra_x[i];

        tg_xx_xxy[i] = tg_0_xxy[i] * fxi[i] + 2.0 * tg_x_xy[i] * fxi[i] + tg_x_xxy[i] * ra_x[i];

        tg_xx_xxz[i] = tg_0_xxz[i] * fxi[i] + 2.0 * tg_x_xz[i] * fxi[i] + tg_x_xxz[i] * ra_x[i];

        tg_xx_xyy[i] = tg_0_xyy[i] * fxi[i] + tg_x_yy[i] * fxi[i] + tg_x_xyy[i] * ra_x[i];

        tg_xx_xyz[i] = tg_0_xyz[i] * fxi[i] + tg_x_yz[i] * fxi[i] + tg_x_xyz[i] * ra_x[i];

        tg_xx_xzz[i] = tg_0_xzz[i] * fxi[i] + tg_x_zz[i] * fxi[i] + tg_x_xzz[i] * ra_x[i];

        tg_xx_yyy[i] = tg_0_yyy[i] * fxi[i] + tg_x_yyy[i] * ra_x[i];

        tg_xx_yyz[i] = tg_0_yyz[i] * fxi[i] + tg_x_yyz[i] * ra_x[i];

        tg_xx_yzz[i] = tg_0_yzz[i] * fxi[i] + tg_x_yzz[i] * ra_x[i];

        tg_xx_zzz[i] = tg_0_zzz[i] * fxi[i] + tg_x_zzz[i] * ra_x[i];

        tg_xy_xxx[i] = tg_x_xxx[i] * ra_y[i];

        tg_xy_xxy[i] = 2.0 * tg_y_xy[i] * fxi[i] + tg_y_xxy[i] * ra_x[i];

        tg_xy_xxz[i] = tg_x_xxz[i] * ra_y[i];

        tg_xy_xyy[i] = tg_y_yy[i] * fxi[i] + tg_y_xyy[i] * ra_x[i];

        tg_xy_xyz[i] = tg_y_yz[i] * fxi[i] + tg_y_xyz[i] * ra_x[i];

        tg_xy_xzz[i] = tg_x_xzz[i] * ra_y[i];

        tg_xy_yyy[i] = tg_y_yyy[i] * ra_x[i];

        tg_xy_yyz[i] = tg_y_yyz[i] * ra_x[i];

        tg_xy_yzz[i] = tg_y_yzz[i] * ra_x[i];

        tg_xy_zzz[i] = tg_y_zzz[i] * ra_x[i];

        tg_xz_xxx[i] = tg_x_xxx[i] * ra_z[i];

        tg_xz_xxy[i] = tg_x_xxy[i] * ra_z[i];

        tg_xz_xxz[i] = 2.0 * tg_z_xz[i] * fxi[i] + tg_z_xxz[i] * ra_x[i];

        tg_xz_xyy[i] = tg_x_xyy[i] * ra_z[i];

        tg_xz_xyz[i] = tg_z_yz[i] * fxi[i] + tg_z_xyz[i] * ra_x[i];

        tg_xz_xzz[i] = tg_z_zz[i] * fxi[i] + tg_z_xzz[i] * ra_x[i];

        tg_xz_yyy[i] = tg_z_yyy[i] * ra_x[i];

        tg_xz_yyz[i] = tg_z_yyz[i] * ra_x[i];

        tg_xz_yzz[i] = tg_z_yzz[i] * ra_x[i];

        tg_xz_zzz[i] = tg_z_zzz[i] * ra_x[i];

        tg_yy_xxx[i] = tg_0_xxx[i] * fxi[i] + tg_y_xxx[i] * ra_y[i];

        tg_yy_xxy[i] = tg_0_xxy[i] * fxi[i] + tg_y_xx[i] * fxi[i] + tg_y_xxy[i] * ra_y[i];

        tg_yy_xxz[i] = tg_0_xxz[i] * fxi[i] + tg_y_xxz[i] * ra_y[i];

        tg_yy_xyy[i] = tg_0_xyy[i] * fxi[i] + 2.0 * tg_y_xy[i] * fxi[i] + tg_y_xyy[i] * ra_y[i];

        tg_yy_xyz[i] = tg_0_xyz[i] * fxi[i] + tg_y_xz[i] * fxi[i] + tg_y_xyz[i] * ra_y[i];

        tg_yy_xzz[i] = tg_0_xzz[i] * fxi[i] + tg_y_xzz[i] * ra_y[i];

        tg_yy_yyy[i] = tg_0_yyy[i] * fxi[i] + 3.0 * tg_y_yy[i] * fxi[i] + tg_y_yyy[i] * ra_y[i];

        tg_yy_yyz[i] = tg_0_yyz[i] * fxi[i] + 2.0 * tg_y_yz[i] * fxi[i] + tg_y_yyz[i] * ra_y[i];

        tg_yy_yzz[i] = tg_0_yzz[i] * fxi[i] + tg_y_zz[i] * fxi[i] + tg_y_yzz[i] * ra_y[i];

        tg_yy_zzz[i] = tg_0_zzz[i] * fxi[i] + tg_y_zzz[i] * ra_y[i];

        tg_yz_xxx[i] = tg_z_xxx[i] * ra_y[i];

        tg_yz_xxy[i] = tg_y_xxy[i] * ra_z[i];

        tg_yz_xxz[i] = tg_z_xxz[i] * ra_y[i];

        tg_yz_xyy[i] = tg_y_xyy[i] * ra_z[i];

        tg_yz_xyz[i] = tg_z_xz[i] * fxi[i] + tg_z_xyz[i] * ra_y[i];

        tg_yz_xzz[i] = tg_z_xzz[i] * ra_y[i];

        tg_yz_yyy[i] = tg_y_yyy[i] * ra_z[i];

        tg_yz_yyz[i] = 2.0 * tg_z_yz[i] * fxi[i] + tg_z_yyz[i] * ra_y[i];

        tg_yz_yzz[i] = tg_z_zz[i] * fxi[i] + tg_z_yzz[i] * ra_y[i];

        tg_yz_zzz[i] = tg_z_zzz[i] * ra_y[i];

        tg_zz_xxx[i] = tg_0_xxx[i] * fxi[i] + tg_z_xxx[i] * ra_z[i];

        tg_zz_xxy[i] = tg_0_xxy[i] * fxi[i] + tg_z_xxy[i] * ra_z[i];

        tg_zz_xxz[i] = tg_0_xxz[i] * fxi[i] + tg_z_xx[i] * fxi[i] + tg_z_xxz[i] * ra_z[i];

        tg_zz_xyy[i] = tg_0_xyy[i] * fxi[i] + tg_z_xyy[i] * ra_z[i];

        tg_zz_xyz[i] = tg_0_xyz[i] * fxi[i] + tg_z_xy[i] * fxi[i] + tg_z_xyz[i] * ra_z[i];

        tg_zz_xzz[i] = tg_0_xzz[i] * fxi[i] + 2.0 * tg_z_xz[i] * fxi[i] + tg_z_xzz[i] * ra_z[i];

        tg_zz_yyy[i] = tg_0_yyy[i] * fxi[i] + tg_z_yyy[i] * ra_z[i];

        tg_zz_yyz[i] = tg_0_yyz[i] * fxi[i] + tg_z_yy[i] * fxi[i] + tg_z_yyz[i] * ra_z[i];

        tg_zz_yzz[i] = tg_0_yzz[i] * fxi[i] + 2.0 * tg_z_yz[i] * fxi[i] + tg_z_yzz[i] * ra_z[i];

        tg_zz_zzz[i] = tg_0_zzz[i] * fxi[i] + 3.0 * tg_z_zz[i] * fxi[i] + tg_z_zzz[i] * ra_z[i];
    }
}

} // t2lecp namespace


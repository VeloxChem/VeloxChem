#include "LocalCorePotentialPrimRecDD.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_dd(CSimdArray<double>& pbuffer, 
                                  const size_t idx_dd,
                                  const size_t idx_sd,
                                  const size_t idx_pp,
                                  const size_t idx_pd,
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

    // Set up components of auxiliary buffer : SD

    auto tg_0_xx = pbuffer.data(idx_sd);

    auto tg_0_xy = pbuffer.data(idx_sd + 1);

    auto tg_0_xz = pbuffer.data(idx_sd + 2);

    auto tg_0_yy = pbuffer.data(idx_sd + 3);

    auto tg_0_yz = pbuffer.data(idx_sd + 4);

    auto tg_0_zz = pbuffer.data(idx_sd + 5);

    // Set up components of auxiliary buffer : PP

    auto tg_x_x = pbuffer.data(idx_pp);

    auto tg_x_y = pbuffer.data(idx_pp + 1);

    auto tg_x_z = pbuffer.data(idx_pp + 2);

    auto tg_y_x = pbuffer.data(idx_pp + 3);

    auto tg_y_y = pbuffer.data(idx_pp + 4);

    auto tg_y_z = pbuffer.data(idx_pp + 5);

    auto tg_z_x = pbuffer.data(idx_pp + 6);

    auto tg_z_y = pbuffer.data(idx_pp + 7);

    auto tg_z_z = pbuffer.data(idx_pp + 8);

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

    // Set up components of targeted buffer : DD

    auto tg_xx_xx = pbuffer.data(idx_dd);

    auto tg_xx_xy = pbuffer.data(idx_dd + 1);

    auto tg_xx_xz = pbuffer.data(idx_dd + 2);

    auto tg_xx_yy = pbuffer.data(idx_dd + 3);

    auto tg_xx_yz = pbuffer.data(idx_dd + 4);

    auto tg_xx_zz = pbuffer.data(idx_dd + 5);

    auto tg_xy_xx = pbuffer.data(idx_dd + 6);

    auto tg_xy_xy = pbuffer.data(idx_dd + 7);

    auto tg_xy_xz = pbuffer.data(idx_dd + 8);

    auto tg_xy_yy = pbuffer.data(idx_dd + 9);

    auto tg_xy_yz = pbuffer.data(idx_dd + 10);

    auto tg_xy_zz = pbuffer.data(idx_dd + 11);

    auto tg_xz_xx = pbuffer.data(idx_dd + 12);

    auto tg_xz_xy = pbuffer.data(idx_dd + 13);

    auto tg_xz_xz = pbuffer.data(idx_dd + 14);

    auto tg_xz_yy = pbuffer.data(idx_dd + 15);

    auto tg_xz_yz = pbuffer.data(idx_dd + 16);

    auto tg_xz_zz = pbuffer.data(idx_dd + 17);

    auto tg_yy_xx = pbuffer.data(idx_dd + 18);

    auto tg_yy_xy = pbuffer.data(idx_dd + 19);

    auto tg_yy_xz = pbuffer.data(idx_dd + 20);

    auto tg_yy_yy = pbuffer.data(idx_dd + 21);

    auto tg_yy_yz = pbuffer.data(idx_dd + 22);

    auto tg_yy_zz = pbuffer.data(idx_dd + 23);

    auto tg_yz_xx = pbuffer.data(idx_dd + 24);

    auto tg_yz_xy = pbuffer.data(idx_dd + 25);

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

    #pragma omp simd aligned(fxi, ra_x, ra_y, ra_z, tg_0_xx, tg_0_xy, tg_0_xz, tg_0_yy, tg_0_yz, tg_0_zz, tg_x_x, tg_x_xx, tg_x_xy, tg_x_xz, tg_x_y, tg_x_yy, tg_x_yz, tg_x_z, tg_x_zz, tg_xx_xx, tg_xx_xy, tg_xx_xz, tg_xx_yy, tg_xx_yz, tg_xx_zz, tg_xy_xx, tg_xy_xy, tg_xy_xz, tg_xy_yy, tg_xy_yz, tg_xy_zz, tg_xz_xx, tg_xz_xy, tg_xz_xz, tg_xz_yy, tg_xz_yz, tg_xz_zz, tg_y_x, tg_y_xx, tg_y_xy, tg_y_xz, tg_y_y, tg_y_yy, tg_y_yz, tg_y_z, tg_y_zz, tg_yy_xx, tg_yy_xy, tg_yy_xz, tg_yy_yy, tg_yy_yz, tg_yy_zz, tg_yz_xx, tg_yz_xy, tg_yz_xz, tg_yz_yy, tg_yz_yz, tg_yz_zz, tg_z_x, tg_z_xx, tg_z_xy, tg_z_xz, tg_z_y, tg_z_yy, tg_z_yz, tg_z_z, tg_z_zz, tg_zz_xx, tg_zz_xy, tg_zz_xz, tg_zz_yy, tg_zz_yz, tg_zz_zz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_xx_xx[i] = tg_0_xx[i] * fxi[i] + 2.0 * tg_x_x[i] * fxi[i] + tg_x_xx[i] * ra_x[i];

        tg_xx_xy[i] = tg_0_xy[i] * fxi[i] + tg_x_y[i] * fxi[i] + tg_x_xy[i] * ra_x[i];

        tg_xx_xz[i] = tg_0_xz[i] * fxi[i] + tg_x_z[i] * fxi[i] + tg_x_xz[i] * ra_x[i];

        tg_xx_yy[i] = tg_0_yy[i] * fxi[i] + tg_x_yy[i] * ra_x[i];

        tg_xx_yz[i] = tg_0_yz[i] * fxi[i] + tg_x_yz[i] * ra_x[i];

        tg_xx_zz[i] = tg_0_zz[i] * fxi[i] + tg_x_zz[i] * ra_x[i];

        tg_xy_xx[i] = tg_x_xx[i] * ra_y[i];

        tg_xy_xy[i] = tg_y_y[i] * fxi[i] + tg_y_xy[i] * ra_x[i];

        tg_xy_xz[i] = tg_x_xz[i] * ra_y[i];

        tg_xy_yy[i] = tg_y_yy[i] * ra_x[i];

        tg_xy_yz[i] = tg_y_yz[i] * ra_x[i];

        tg_xy_zz[i] = tg_y_zz[i] * ra_x[i];

        tg_xz_xx[i] = tg_x_xx[i] * ra_z[i];

        tg_xz_xy[i] = tg_x_xy[i] * ra_z[i];

        tg_xz_xz[i] = tg_z_z[i] * fxi[i] + tg_z_xz[i] * ra_x[i];

        tg_xz_yy[i] = tg_z_yy[i] * ra_x[i];

        tg_xz_yz[i] = tg_z_yz[i] * ra_x[i];

        tg_xz_zz[i] = tg_z_zz[i] * ra_x[i];

        tg_yy_xx[i] = tg_0_xx[i] * fxi[i] + tg_y_xx[i] * ra_y[i];

        tg_yy_xy[i] = tg_0_xy[i] * fxi[i] + tg_y_x[i] * fxi[i] + tg_y_xy[i] * ra_y[i];

        tg_yy_xz[i] = tg_0_xz[i] * fxi[i] + tg_y_xz[i] * ra_y[i];

        tg_yy_yy[i] = tg_0_yy[i] * fxi[i] + 2.0 * tg_y_y[i] * fxi[i] + tg_y_yy[i] * ra_y[i];

        tg_yy_yz[i] = tg_0_yz[i] * fxi[i] + tg_y_z[i] * fxi[i] + tg_y_yz[i] * ra_y[i];

        tg_yy_zz[i] = tg_0_zz[i] * fxi[i] + tg_y_zz[i] * ra_y[i];

        tg_yz_xx[i] = tg_z_xx[i] * ra_y[i];

        tg_yz_xy[i] = tg_y_xy[i] * ra_z[i];

        tg_yz_xz[i] = tg_z_xz[i] * ra_y[i];

        tg_yz_yy[i] = tg_y_yy[i] * ra_z[i];

        tg_yz_yz[i] = tg_z_z[i] * fxi[i] + tg_z_yz[i] * ra_y[i];

        tg_yz_zz[i] = tg_z_zz[i] * ra_y[i];

        tg_zz_xx[i] = tg_0_xx[i] * fxi[i] + tg_z_xx[i] * ra_z[i];

        tg_zz_xy[i] = tg_0_xy[i] * fxi[i] + tg_z_xy[i] * ra_z[i];

        tg_zz_xz[i] = tg_0_xz[i] * fxi[i] + tg_z_x[i] * fxi[i] + tg_z_xz[i] * ra_z[i];

        tg_zz_yy[i] = tg_0_yy[i] * fxi[i] + tg_z_yy[i] * ra_z[i];

        tg_zz_yz[i] = tg_0_yz[i] * fxi[i] + tg_z_y[i] * fxi[i] + tg_z_yz[i] * ra_z[i];

        tg_zz_zz[i] = tg_0_zz[i] * fxi[i] + 2.0 * tg_z_z[i] * fxi[i] + tg_z_zz[i] * ra_z[i];
    }
}

} // t2lecp namespace


#include "LocalCorePotentialPrimRecPF.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_pf(CSimdArray<double>& pbuffer, 
                                  const size_t idx_pf,
                                  const size_t idx_sd,
                                  const size_t idx_sf,
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

    // Set up components of targeted buffer : PF

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

    #pragma omp simd aligned(fxi, ra_x, ra_y, ra_z, tg_0_xx, tg_0_xxx, tg_0_xxy, tg_0_xxz, tg_0_xy, tg_0_xyy, tg_0_xyz, tg_0_xz, tg_0_xzz, tg_0_yy, tg_0_yyy, tg_0_yyz, tg_0_yz, tg_0_yzz, tg_0_zz, tg_0_zzz, tg_x_xxx, tg_x_xxy, tg_x_xxz, tg_x_xyy, tg_x_xyz, tg_x_xzz, tg_x_yyy, tg_x_yyz, tg_x_yzz, tg_x_zzz, tg_y_xxx, tg_y_xxy, tg_y_xxz, tg_y_xyy, tg_y_xyz, tg_y_xzz, tg_y_yyy, tg_y_yyz, tg_y_yzz, tg_y_zzz, tg_z_xxx, tg_z_xxy, tg_z_xxz, tg_z_xyy, tg_z_xyz, tg_z_xzz, tg_z_yyy, tg_z_yyz, tg_z_yzz, tg_z_zzz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_x_xxx[i] = 3.0 * tg_0_xx[i] * fxi[i] + tg_0_xxx[i] * ra_x[i];

        tg_x_xxy[i] = 2.0 * tg_0_xy[i] * fxi[i] + tg_0_xxy[i] * ra_x[i];

        tg_x_xxz[i] = 2.0 * tg_0_xz[i] * fxi[i] + tg_0_xxz[i] * ra_x[i];

        tg_x_xyy[i] = tg_0_yy[i] * fxi[i] + tg_0_xyy[i] * ra_x[i];

        tg_x_xyz[i] = tg_0_yz[i] * fxi[i] + tg_0_xyz[i] * ra_x[i];

        tg_x_xzz[i] = tg_0_zz[i] * fxi[i] + tg_0_xzz[i] * ra_x[i];

        tg_x_yyy[i] = tg_0_yyy[i] * ra_x[i];

        tg_x_yyz[i] = tg_0_yyz[i] * ra_x[i];

        tg_x_yzz[i] = tg_0_yzz[i] * ra_x[i];

        tg_x_zzz[i] = tg_0_zzz[i] * ra_x[i];

        tg_y_xxx[i] = tg_0_xxx[i] * ra_y[i];

        tg_y_xxy[i] = tg_0_xx[i] * fxi[i] + tg_0_xxy[i] * ra_y[i];

        tg_y_xxz[i] = tg_0_xxz[i] * ra_y[i];

        tg_y_xyy[i] = 2.0 * tg_0_xy[i] * fxi[i] + tg_0_xyy[i] * ra_y[i];

        tg_y_xyz[i] = tg_0_xz[i] * fxi[i] + tg_0_xyz[i] * ra_y[i];

        tg_y_xzz[i] = tg_0_xzz[i] * ra_y[i];

        tg_y_yyy[i] = 3.0 * tg_0_yy[i] * fxi[i] + tg_0_yyy[i] * ra_y[i];

        tg_y_yyz[i] = 2.0 * tg_0_yz[i] * fxi[i] + tg_0_yyz[i] * ra_y[i];

        tg_y_yzz[i] = tg_0_zz[i] * fxi[i] + tg_0_yzz[i] * ra_y[i];

        tg_y_zzz[i] = tg_0_zzz[i] * ra_y[i];

        tg_z_xxx[i] = tg_0_xxx[i] * ra_z[i];

        tg_z_xxy[i] = tg_0_xxy[i] * ra_z[i];

        tg_z_xxz[i] = tg_0_xx[i] * fxi[i] + tg_0_xxz[i] * ra_z[i];

        tg_z_xyy[i] = tg_0_xyy[i] * ra_z[i];

        tg_z_xyz[i] = tg_0_xy[i] * fxi[i] + tg_0_xyz[i] * ra_z[i];

        tg_z_xzz[i] = 2.0 * tg_0_xz[i] * fxi[i] + tg_0_xzz[i] * ra_z[i];

        tg_z_yyy[i] = tg_0_yyy[i] * ra_z[i];

        tg_z_yyz[i] = tg_0_yy[i] * fxi[i] + tg_0_yyz[i] * ra_z[i];

        tg_z_yzz[i] = 2.0 * tg_0_yz[i] * fxi[i] + tg_0_yzz[i] * ra_z[i];

        tg_z_zzz[i] = 3.0 * tg_0_zz[i] * fxi[i] + tg_0_zzz[i] * ra_z[i];
    }
}

} // t2lecp namespace


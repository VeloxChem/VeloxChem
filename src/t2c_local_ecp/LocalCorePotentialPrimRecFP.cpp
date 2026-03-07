#include "LocalCorePotentialPrimRecFP.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_fp(CSimdArray<double>& pbuffer, 
                                  const size_t idx_fp,
                                  const size_t idx_pp,
                                  const size_t idx_ds,
                                  const size_t idx_dp,
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

    // Set up components of auxiliary buffer : DS

    auto tg_xx_0 = pbuffer.data(idx_ds);

    auto tg_yy_0 = pbuffer.data(idx_ds + 3);

    auto tg_zz_0 = pbuffer.data(idx_ds + 5);

    // Set up components of auxiliary buffer : DP

    auto tg_xx_x = pbuffer.data(idx_dp);

    auto tg_xx_y = pbuffer.data(idx_dp + 1);

    auto tg_xx_z = pbuffer.data(idx_dp + 2);

    auto tg_xy_y = pbuffer.data(idx_dp + 4);

    auto tg_xz_x = pbuffer.data(idx_dp + 6);

    auto tg_xz_z = pbuffer.data(idx_dp + 8);

    auto tg_yy_x = pbuffer.data(idx_dp + 9);

    auto tg_yy_y = pbuffer.data(idx_dp + 10);

    auto tg_yy_z = pbuffer.data(idx_dp + 11);

    auto tg_yz_y = pbuffer.data(idx_dp + 13);

    auto tg_yz_z = pbuffer.data(idx_dp + 14);

    auto tg_zz_x = pbuffer.data(idx_dp + 15);

    auto tg_zz_y = pbuffer.data(idx_dp + 16);

    auto tg_zz_z = pbuffer.data(idx_dp + 17);

    // Set up components of targeted buffer : FP

    auto tg_xxx_x = pbuffer.data(idx_fp);

    auto tg_xxx_y = pbuffer.data(idx_fp + 1);

    auto tg_xxx_z = pbuffer.data(idx_fp + 2);

    auto tg_xxy_x = pbuffer.data(idx_fp + 3);

    auto tg_xxy_y = pbuffer.data(idx_fp + 4);

    auto tg_xxy_z = pbuffer.data(idx_fp + 5);

    auto tg_xxz_x = pbuffer.data(idx_fp + 6);

    auto tg_xxz_y = pbuffer.data(idx_fp + 7);

    auto tg_xxz_z = pbuffer.data(idx_fp + 8);

    auto tg_xyy_x = pbuffer.data(idx_fp + 9);

    auto tg_xyy_y = pbuffer.data(idx_fp + 10);

    auto tg_xyy_z = pbuffer.data(idx_fp + 11);

    auto tg_xyz_x = pbuffer.data(idx_fp + 12);

    auto tg_xyz_y = pbuffer.data(idx_fp + 13);

    auto tg_xyz_z = pbuffer.data(idx_fp + 14);

    auto tg_xzz_x = pbuffer.data(idx_fp + 15);

    auto tg_xzz_y = pbuffer.data(idx_fp + 16);

    auto tg_xzz_z = pbuffer.data(idx_fp + 17);

    auto tg_yyy_x = pbuffer.data(idx_fp + 18);

    auto tg_yyy_y = pbuffer.data(idx_fp + 19);

    auto tg_yyy_z = pbuffer.data(idx_fp + 20);

    auto tg_yyz_x = pbuffer.data(idx_fp + 21);

    auto tg_yyz_y = pbuffer.data(idx_fp + 22);

    auto tg_yyz_z = pbuffer.data(idx_fp + 23);

    auto tg_yzz_x = pbuffer.data(idx_fp + 24);

    auto tg_yzz_y = pbuffer.data(idx_fp + 25);

    auto tg_yzz_z = pbuffer.data(idx_fp + 26);

    auto tg_zzz_x = pbuffer.data(idx_fp + 27);

    auto tg_zzz_y = pbuffer.data(idx_fp + 28);

    auto tg_zzz_z = pbuffer.data(idx_fp + 29);

    #pragma omp simd aligned(fxi, ra_x, ra_y, ra_z, tg_x_x, tg_x_y, tg_x_z, tg_xx_0, tg_xx_x, tg_xx_y, tg_xx_z, tg_xxx_x, tg_xxx_y, tg_xxx_z, tg_xxy_x, tg_xxy_y, tg_xxy_z, tg_xxz_x, tg_xxz_y, tg_xxz_z, tg_xy_y, tg_xyy_x, tg_xyy_y, tg_xyy_z, tg_xyz_x, tg_xyz_y, tg_xyz_z, tg_xz_x, tg_xz_z, tg_xzz_x, tg_xzz_y, tg_xzz_z, tg_y_x, tg_y_y, tg_y_z, tg_yy_0, tg_yy_x, tg_yy_y, tg_yy_z, tg_yyy_x, tg_yyy_y, tg_yyy_z, tg_yyz_x, tg_yyz_y, tg_yyz_z, tg_yz_y, tg_yz_z, tg_yzz_x, tg_yzz_y, tg_yzz_z, tg_z_x, tg_z_y, tg_z_z, tg_zz_0, tg_zz_x, tg_zz_y, tg_zz_z, tg_zzz_x, tg_zzz_y, tg_zzz_z  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_xxx_x[i] = 2.0 * tg_x_x[i] * fxi[i] + tg_xx_0[i] * fxi[i] + tg_xx_x[i] * ra_x[i];

        tg_xxx_y[i] = 2.0 * tg_x_y[i] * fxi[i] + tg_xx_y[i] * ra_x[i];

        tg_xxx_z[i] = 2.0 * tg_x_z[i] * fxi[i] + tg_xx_z[i] * ra_x[i];

        tg_xxy_x[i] = tg_xx_x[i] * ra_y[i];

        tg_xxy_y[i] = tg_y_y[i] * fxi[i] + tg_xy_y[i] * ra_x[i];

        tg_xxy_z[i] = tg_xx_z[i] * ra_y[i];

        tg_xxz_x[i] = tg_xx_x[i] * ra_z[i];

        tg_xxz_y[i] = tg_xx_y[i] * ra_z[i];

        tg_xxz_z[i] = tg_z_z[i] * fxi[i] + tg_xz_z[i] * ra_x[i];

        tg_xyy_x[i] = tg_yy_0[i] * fxi[i] + tg_yy_x[i] * ra_x[i];

        tg_xyy_y[i] = tg_yy_y[i] * ra_x[i];

        tg_xyy_z[i] = tg_yy_z[i] * ra_x[i];

        tg_xyz_x[i] = tg_xz_x[i] * ra_y[i];

        tg_xyz_y[i] = tg_yz_y[i] * ra_x[i];

        tg_xyz_z[i] = tg_yz_z[i] * ra_x[i];

        tg_xzz_x[i] = tg_zz_0[i] * fxi[i] + tg_zz_x[i] * ra_x[i];

        tg_xzz_y[i] = tg_zz_y[i] * ra_x[i];

        tg_xzz_z[i] = tg_zz_z[i] * ra_x[i];

        tg_yyy_x[i] = 2.0 * tg_y_x[i] * fxi[i] + tg_yy_x[i] * ra_y[i];

        tg_yyy_y[i] = 2.0 * tg_y_y[i] * fxi[i] + tg_yy_0[i] * fxi[i] + tg_yy_y[i] * ra_y[i];

        tg_yyy_z[i] = 2.0 * tg_y_z[i] * fxi[i] + tg_yy_z[i] * ra_y[i];

        tg_yyz_x[i] = tg_yy_x[i] * ra_z[i];

        tg_yyz_y[i] = tg_yy_y[i] * ra_z[i];

        tg_yyz_z[i] = tg_z_z[i] * fxi[i] + tg_yz_z[i] * ra_y[i];

        tg_yzz_x[i] = tg_zz_x[i] * ra_y[i];

        tg_yzz_y[i] = tg_zz_0[i] * fxi[i] + tg_zz_y[i] * ra_y[i];

        tg_yzz_z[i] = tg_zz_z[i] * ra_y[i];

        tg_zzz_x[i] = 2.0 * tg_z_x[i] * fxi[i] + tg_zz_x[i] * ra_z[i];

        tg_zzz_y[i] = 2.0 * tg_z_y[i] * fxi[i] + tg_zz_y[i] * ra_z[i];

        tg_zzz_z[i] = 2.0 * tg_z_z[i] * fxi[i] + tg_zz_0[i] * fxi[i] + tg_zz_z[i] * ra_z[i];
    }
}

} // t2lecp namespace


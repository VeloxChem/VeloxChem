#include "LocalCorePotentialPrimRecDP.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_dp(CSimdArray<double>& pbuffer, 
                                  const size_t idx_dp,
                                  const size_t idx_sp,
                                  const size_t idx_ps,
                                  const size_t idx_pp,
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

    // Set up components of auxiliary buffer : SP

    auto tg_0_x = pbuffer.data(idx_sp);

    auto tg_0_y = pbuffer.data(idx_sp + 1);

    auto tg_0_z = pbuffer.data(idx_sp + 2);

    // Set up components of auxiliary buffer : PS

    auto tg_x_0 = pbuffer.data(idx_ps);

    auto tg_y_0 = pbuffer.data(idx_ps + 1);

    auto tg_z_0 = pbuffer.data(idx_ps + 2);

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

    // Set up components of targeted buffer : DP

    auto tg_xx_x = pbuffer.data(idx_dp);

    auto tg_xx_y = pbuffer.data(idx_dp + 1);

    auto tg_xx_z = pbuffer.data(idx_dp + 2);

    auto tg_xy_x = pbuffer.data(idx_dp + 3);

    auto tg_xy_y = pbuffer.data(idx_dp + 4);

    auto tg_xy_z = pbuffer.data(idx_dp + 5);

    auto tg_xz_x = pbuffer.data(idx_dp + 6);

    auto tg_xz_y = pbuffer.data(idx_dp + 7);

    auto tg_xz_z = pbuffer.data(idx_dp + 8);

    auto tg_yy_x = pbuffer.data(idx_dp + 9);

    auto tg_yy_y = pbuffer.data(idx_dp + 10);

    auto tg_yy_z = pbuffer.data(idx_dp + 11);

    auto tg_yz_x = pbuffer.data(idx_dp + 12);

    auto tg_yz_y = pbuffer.data(idx_dp + 13);

    auto tg_yz_z = pbuffer.data(idx_dp + 14);

    auto tg_zz_x = pbuffer.data(idx_dp + 15);

    auto tg_zz_y = pbuffer.data(idx_dp + 16);

    auto tg_zz_z = pbuffer.data(idx_dp + 17);

    #pragma omp simd aligned(fxi, ra_x, ra_y, ra_z, tg_0_x, tg_0_y, tg_0_z, tg_x_0, tg_x_x, tg_x_y, tg_x_z, tg_xx_x, tg_xx_y, tg_xx_z, tg_xy_x, tg_xy_y, tg_xy_z, tg_xz_x, tg_xz_y, tg_xz_z, tg_y_0, tg_y_x, tg_y_y, tg_y_z, tg_yy_x, tg_yy_y, tg_yy_z, tg_yz_x, tg_yz_y, tg_yz_z, tg_z_0, tg_z_x, tg_z_y, tg_z_z, tg_zz_x, tg_zz_y, tg_zz_z  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_xx_x[i] = tg_0_x[i] * fxi[i] + tg_x_0[i] * fxi[i] + tg_x_x[i] * ra_x[i];

        tg_xx_y[i] = tg_0_y[i] * fxi[i] + tg_x_y[i] * ra_x[i];

        tg_xx_z[i] = tg_0_z[i] * fxi[i] + tg_x_z[i] * ra_x[i];

        tg_xy_x[i] = tg_x_x[i] * ra_y[i];

        tg_xy_y[i] = tg_y_y[i] * ra_x[i];

        tg_xy_z[i] = tg_y_z[i] * ra_x[i];

        tg_xz_x[i] = tg_x_x[i] * ra_z[i];

        tg_xz_y[i] = tg_z_y[i] * ra_x[i];

        tg_xz_z[i] = tg_z_z[i] * ra_x[i];

        tg_yy_x[i] = tg_0_x[i] * fxi[i] + tg_y_x[i] * ra_y[i];

        tg_yy_y[i] = tg_0_y[i] * fxi[i] + tg_y_0[i] * fxi[i] + tg_y_y[i] * ra_y[i];

        tg_yy_z[i] = tg_0_z[i] * fxi[i] + tg_y_z[i] * ra_y[i];

        tg_yz_x[i] = tg_z_x[i] * ra_y[i];

        tg_yz_y[i] = tg_y_y[i] * ra_z[i];

        tg_yz_z[i] = tg_z_z[i] * ra_y[i];

        tg_zz_x[i] = tg_0_x[i] * fxi[i] + tg_z_x[i] * ra_z[i];

        tg_zz_y[i] = tg_0_y[i] * fxi[i] + tg_z_y[i] * ra_z[i];

        tg_zz_z[i] = tg_0_z[i] * fxi[i] + tg_z_0[i] * fxi[i] + tg_z_z[i] * ra_z[i];
    }
}

} // t2lecp namespace


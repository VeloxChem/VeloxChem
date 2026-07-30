#include "LocalCorePotentialPrimRecPD.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_pd(CSimdArray<double>& pbuffer, 
                                  const size_t idx_pd,
                                  const size_t idx_sp,
                                  const size_t idx_sd,
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

    // Set up components of auxiliary buffer : SD

    auto tg_0_xx = pbuffer.data(idx_sd);

    auto tg_0_xy = pbuffer.data(idx_sd + 1);

    auto tg_0_xz = pbuffer.data(idx_sd + 2);

    auto tg_0_yy = pbuffer.data(idx_sd + 3);

    auto tg_0_yz = pbuffer.data(idx_sd + 4);

    auto tg_0_zz = pbuffer.data(idx_sd + 5);

    // Set up components of targeted buffer : PD

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

    #pragma omp simd aligned(fxi, ra_x, ra_y, ra_z, tg_0_x, tg_0_xx, tg_0_xy, tg_0_xz, tg_0_y, tg_0_yy, tg_0_yz, tg_0_z, tg_0_zz, tg_x_xx, tg_x_xy, tg_x_xz, tg_x_yy, tg_x_yz, tg_x_zz, tg_y_xx, tg_y_xy, tg_y_xz, tg_y_yy, tg_y_yz, tg_y_zz, tg_z_xx, tg_z_xy, tg_z_xz, tg_z_yy, tg_z_yz, tg_z_zz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_x_xx[i] = 2.0 * tg_0_x[i] * fxi[i] + tg_0_xx[i] * ra_x[i];

        tg_x_xy[i] = tg_0_y[i] * fxi[i] + tg_0_xy[i] * ra_x[i];

        tg_x_xz[i] = tg_0_z[i] * fxi[i] + tg_0_xz[i] * ra_x[i];

        tg_x_yy[i] = tg_0_yy[i] * ra_x[i];

        tg_x_yz[i] = tg_0_yz[i] * ra_x[i];

        tg_x_zz[i] = tg_0_zz[i] * ra_x[i];

        tg_y_xx[i] = tg_0_xx[i] * ra_y[i];

        tg_y_xy[i] = tg_0_x[i] * fxi[i] + tg_0_xy[i] * ra_y[i];

        tg_y_xz[i] = tg_0_xz[i] * ra_y[i];

        tg_y_yy[i] = 2.0 * tg_0_y[i] * fxi[i] + tg_0_yy[i] * ra_y[i];

        tg_y_yz[i] = tg_0_z[i] * fxi[i] + tg_0_yz[i] * ra_y[i];

        tg_y_zz[i] = tg_0_zz[i] * ra_y[i];

        tg_z_xx[i] = tg_0_xx[i] * ra_z[i];

        tg_z_xy[i] = tg_0_xy[i] * ra_z[i];

        tg_z_xz[i] = tg_0_x[i] * fxi[i] + tg_0_xz[i] * ra_z[i];

        tg_z_yy[i] = tg_0_yy[i] * ra_z[i];

        tg_z_yz[i] = tg_0_y[i] * fxi[i] + tg_0_yz[i] * ra_z[i];

        tg_z_zz[i] = 2.0 * tg_0_z[i] * fxi[i] + tg_0_zz[i] * ra_z[i];
    }
}

} // t2lecp namespace


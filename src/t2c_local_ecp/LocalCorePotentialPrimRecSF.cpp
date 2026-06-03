#include "LocalCorePotentialPrimRecSF.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_sf(CSimdArray<double>& pbuffer, 
                                  const size_t idx_sf,
                                  const size_t idx_sp,
                                  const size_t idx_sd,
                                  const CSimdArray<double>& factors,
                                  const size_t idx_rb,
                                  const size_t idx_zeta) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up R(RB) distances

    auto rb_x = factors.data(idx_rb);

    auto rb_y = factors.data(idx_rb + 1);

    auto rb_z = factors.data(idx_rb + 2);

    // Set up inverted 1/2xi

    auto fxi = factors.data(idx_zeta);

    // Set up components of auxiliary buffer : SP

    auto tg_0_x = pbuffer.data(idx_sp);

    auto tg_0_y = pbuffer.data(idx_sp + 1);

    auto tg_0_z = pbuffer.data(idx_sp + 2);

    // Set up components of auxiliary buffer : SD

    auto tg_0_xx = pbuffer.data(idx_sd);

    auto tg_0_yy = pbuffer.data(idx_sd + 3);

    auto tg_0_yz = pbuffer.data(idx_sd + 4);

    auto tg_0_zz = pbuffer.data(idx_sd + 5);

    // Set up components of targeted buffer : SF

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

    #pragma omp simd aligned(fxi, rb_x, rb_y, rb_z, tg_0_x, tg_0_xx, tg_0_xxx, tg_0_xxy, tg_0_xxz, tg_0_xyy, tg_0_xyz, tg_0_xzz, tg_0_y, tg_0_yy, tg_0_yyy, tg_0_yyz, tg_0_yz, tg_0_yzz, tg_0_z, tg_0_zz, tg_0_zzz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_0_xxx[i] = 2.0 * tg_0_x[i] * fxi[i] + tg_0_xx[i] * rb_x[i];

        tg_0_xxy[i] = tg_0_xx[i] * rb_y[i];

        tg_0_xxz[i] = tg_0_xx[i] * rb_z[i];

        tg_0_xyy[i] = tg_0_yy[i] * rb_x[i];

        tg_0_xyz[i] = tg_0_yz[i] * rb_x[i];

        tg_0_xzz[i] = tg_0_zz[i] * rb_x[i];

        tg_0_yyy[i] = 2.0 * tg_0_y[i] * fxi[i] + tg_0_yy[i] * rb_y[i];

        tg_0_yyz[i] = tg_0_yy[i] * rb_z[i];

        tg_0_yzz[i] = tg_0_zz[i] * rb_y[i];

        tg_0_zzz[i] = 2.0 * tg_0_z[i] * fxi[i] + tg_0_zz[i] * rb_z[i];
    }
}

} // t2lecp namespace


#include "LocalCorePotentialPrimRecSD.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_sd(CSimdArray<double>& pbuffer, 
                                  const size_t idx_sd,
                                  const size_t idx_ss,
                                  const size_t idx_sp,
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

    // Set up components of auxiliary buffer : SS

    auto tg_0_0 = pbuffer.data(idx_ss);

    // Set up components of auxiliary buffer : SP

    auto tg_0_x = pbuffer.data(idx_sp);

    auto tg_0_y = pbuffer.data(idx_sp + 1);

    auto tg_0_z = pbuffer.data(idx_sp + 2);

    // Set up components of targeted buffer : SD

    auto tg_0_xx = pbuffer.data(idx_sd);

    auto tg_0_xy = pbuffer.data(idx_sd + 1);

    auto tg_0_xz = pbuffer.data(idx_sd + 2);

    auto tg_0_yy = pbuffer.data(idx_sd + 3);

    auto tg_0_yz = pbuffer.data(idx_sd + 4);

    auto tg_0_zz = pbuffer.data(idx_sd + 5);

    #pragma omp simd aligned(fxi, rb_x, rb_y, rb_z, tg_0_0, tg_0_x, tg_0_xx, tg_0_xy, tg_0_xz, tg_0_y, tg_0_yy, tg_0_yz, tg_0_z, tg_0_zz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_0_xx[i] = tg_0_0[i] * fxi[i] + tg_0_x[i] * rb_x[i];

        tg_0_xy[i] = tg_0_y[i] * rb_x[i];

        tg_0_xz[i] = tg_0_z[i] * rb_x[i];

        tg_0_yy[i] = tg_0_0[i] * fxi[i] + tg_0_y[i] * rb_y[i];

        tg_0_yz[i] = tg_0_z[i] * rb_y[i];

        tg_0_zz[i] = tg_0_0[i] * fxi[i] + tg_0_z[i] * rb_z[i];
    }
}

} // t2lecp namespace


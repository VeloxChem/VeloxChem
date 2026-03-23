#include "LocalCorePotentialPrimRecSG.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_sg(CSimdArray<double>& pbuffer, 
                                  const size_t idx_sg,
                                  const size_t idx_sd,
                                  const size_t idx_sf,
                                  const CSimdArray<double>& factors) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up R(RB) distances

    auto rb_x = factors.data(8);

    auto rb_y = factors.data(9);

    auto rb_z = factors.data(10);

    // Set up inverted 1/2xi

    auto fxi = factors.data(11);

    // Set up components of auxiliary buffer : SD

    auto tg_0_xx = pbuffer.data(idx_sd);

    auto tg_0_yy = pbuffer.data(idx_sd + 3);

    auto tg_0_zz = pbuffer.data(idx_sd + 5);

    // Set up components of auxiliary buffer : SF

    auto tg_0_xxx = pbuffer.data(idx_sf);

    auto tg_0_xxz = pbuffer.data(idx_sf + 2);

    auto tg_0_xyy = pbuffer.data(idx_sf + 3);

    auto tg_0_xzz = pbuffer.data(idx_sf + 5);

    auto tg_0_yyy = pbuffer.data(idx_sf + 6);

    auto tg_0_yyz = pbuffer.data(idx_sf + 7);

    auto tg_0_yzz = pbuffer.data(idx_sf + 8);

    auto tg_0_zzz = pbuffer.data(idx_sf + 9);

    // Set up components of targeted buffer : SG

    auto tg_0_xxxx = pbuffer.data(idx_sg);

    auto tg_0_xxxy = pbuffer.data(idx_sg + 1);

    auto tg_0_xxxz = pbuffer.data(idx_sg + 2);

    auto tg_0_xxyy = pbuffer.data(idx_sg + 3);

    auto tg_0_xxyz = pbuffer.data(idx_sg + 4);

    auto tg_0_xxzz = pbuffer.data(idx_sg + 5);

    auto tg_0_xyyy = pbuffer.data(idx_sg + 6);

    auto tg_0_xyyz = pbuffer.data(idx_sg + 7);

    auto tg_0_xyzz = pbuffer.data(idx_sg + 8);

    auto tg_0_xzzz = pbuffer.data(idx_sg + 9);

    auto tg_0_yyyy = pbuffer.data(idx_sg + 10);

    auto tg_0_yyyz = pbuffer.data(idx_sg + 11);

    auto tg_0_yyzz = pbuffer.data(idx_sg + 12);

    auto tg_0_yzzz = pbuffer.data(idx_sg + 13);

    auto tg_0_zzzz = pbuffer.data(idx_sg + 14);

    #pragma omp simd aligned(fxi, rb_x, rb_y, rb_z, tg_0_xx, tg_0_xxx, tg_0_xxxx, tg_0_xxxy, tg_0_xxxz, tg_0_xxyy, tg_0_xxyz, tg_0_xxz, tg_0_xxzz, tg_0_xyy, tg_0_xyyy, tg_0_xyyz, tg_0_xyzz, tg_0_xzz, tg_0_xzzz, tg_0_yy, tg_0_yyy, tg_0_yyyy, tg_0_yyyz, tg_0_yyz, tg_0_yyzz, tg_0_yzz, tg_0_yzzz, tg_0_zz, tg_0_zzz, tg_0_zzzz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_0_xxxx[i] = 3.0 * tg_0_xx[i] * fxi[i] + tg_0_xxx[i] * rb_x[i];

        tg_0_xxxy[i] = tg_0_xxx[i] * rb_y[i];

        tg_0_xxxz[i] = tg_0_xxx[i] * rb_z[i];

        tg_0_xxyy[i] = tg_0_yy[i] * fxi[i] + tg_0_xyy[i] * rb_x[i];

        tg_0_xxyz[i] = tg_0_xxz[i] * rb_y[i];

        tg_0_xxzz[i] = tg_0_zz[i] * fxi[i] + tg_0_xzz[i] * rb_x[i];

        tg_0_xyyy[i] = tg_0_yyy[i] * rb_x[i];

        tg_0_xyyz[i] = tg_0_yyz[i] * rb_x[i];

        tg_0_xyzz[i] = tg_0_yzz[i] * rb_x[i];

        tg_0_xzzz[i] = tg_0_zzz[i] * rb_x[i];

        tg_0_yyyy[i] = 3.0 * tg_0_yy[i] * fxi[i] + tg_0_yyy[i] * rb_y[i];

        tg_0_yyyz[i] = tg_0_yyy[i] * rb_z[i];

        tg_0_yyzz[i] = tg_0_zz[i] * fxi[i] + tg_0_yzz[i] * rb_y[i];

        tg_0_yzzz[i] = tg_0_zzz[i] * rb_y[i];

        tg_0_zzzz[i] = 3.0 * tg_0_zz[i] * fxi[i] + tg_0_zzz[i] * rb_z[i];
    }
}

} // t2lecp namespace


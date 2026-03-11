#include "LocalCorePotentialPrimRecSH.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_sh(CSimdArray<double>& pbuffer, 
                                  const size_t idx_sh,
                                  const size_t idx_sf,
                                  const size_t idx_sg,
                                  const CSimdArray<double>& factors) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up R(RB) distances

    auto rb_x = factors.data(8);

    auto rb_y = factors.data(9);

    auto rb_z = factors.data(10);

    // Set up inverted 1/2xi

    auto fxi = factors.data(11);

    // Set up components of auxiliary buffer : SF

    auto tg_0_xxx = pbuffer.data(idx_sf);

    auto tg_0_xyy = pbuffer.data(idx_sf + 3);

    auto tg_0_xzz = pbuffer.data(idx_sf + 5);

    auto tg_0_yyy = pbuffer.data(idx_sf + 6);

    auto tg_0_yzz = pbuffer.data(idx_sf + 8);

    auto tg_0_zzz = pbuffer.data(idx_sf + 9);

    // Set up components of auxiliary buffer : SG

    auto tg_0_xxxx = pbuffer.data(idx_sg);

    auto tg_0_xxxz = pbuffer.data(idx_sg + 2);

    auto tg_0_xxyy = pbuffer.data(idx_sg + 3);

    auto tg_0_xxzz = pbuffer.data(idx_sg + 5);

    auto tg_0_xyyy = pbuffer.data(idx_sg + 6);

    auto tg_0_xzzz = pbuffer.data(idx_sg + 9);

    auto tg_0_yyyy = pbuffer.data(idx_sg + 10);

    auto tg_0_yyyz = pbuffer.data(idx_sg + 11);

    auto tg_0_yyzz = pbuffer.data(idx_sg + 12);

    auto tg_0_yzzz = pbuffer.data(idx_sg + 13);

    auto tg_0_zzzz = pbuffer.data(idx_sg + 14);

    // Set up components of targeted buffer : SH

    auto tg_0_xxxxx = pbuffer.data(idx_sh);

    auto tg_0_xxxxy = pbuffer.data(idx_sh + 1);

    auto tg_0_xxxxz = pbuffer.data(idx_sh + 2);

    auto tg_0_xxxyy = pbuffer.data(idx_sh + 3);

    auto tg_0_xxxyz = pbuffer.data(idx_sh + 4);

    auto tg_0_xxxzz = pbuffer.data(idx_sh + 5);

    auto tg_0_xxyyy = pbuffer.data(idx_sh + 6);

    auto tg_0_xxyyz = pbuffer.data(idx_sh + 7);

    auto tg_0_xxyzz = pbuffer.data(idx_sh + 8);

    auto tg_0_xxzzz = pbuffer.data(idx_sh + 9);

    auto tg_0_xyyyy = pbuffer.data(idx_sh + 10);

    auto tg_0_xyyyz = pbuffer.data(idx_sh + 11);

    auto tg_0_xyyzz = pbuffer.data(idx_sh + 12);

    auto tg_0_xyzzz = pbuffer.data(idx_sh + 13);

    auto tg_0_xzzzz = pbuffer.data(idx_sh + 14);

    auto tg_0_yyyyy = pbuffer.data(idx_sh + 15);

    auto tg_0_yyyyz = pbuffer.data(idx_sh + 16);

    auto tg_0_yyyzz = pbuffer.data(idx_sh + 17);

    auto tg_0_yyzzz = pbuffer.data(idx_sh + 18);

    auto tg_0_yzzzz = pbuffer.data(idx_sh + 19);

    auto tg_0_zzzzz = pbuffer.data(idx_sh + 20);

    #pragma omp simd aligned(fxi, rb_x, rb_y, rb_z, tg_0_xxx, tg_0_xxxx, tg_0_xxxxx, tg_0_xxxxy, tg_0_xxxxz, tg_0_xxxyy, tg_0_xxxyz, tg_0_xxxz, tg_0_xxxzz, tg_0_xxyy, tg_0_xxyyy, tg_0_xxyyz, tg_0_xxyzz, tg_0_xxzz, tg_0_xxzzz, tg_0_xyy, tg_0_xyyy, tg_0_xyyyy, tg_0_xyyyz, tg_0_xyyzz, tg_0_xyzzz, tg_0_xzz, tg_0_xzzz, tg_0_xzzzz, tg_0_yyy, tg_0_yyyy, tg_0_yyyyy, tg_0_yyyyz, tg_0_yyyz, tg_0_yyyzz, tg_0_yyzz, tg_0_yyzzz, tg_0_yzz, tg_0_yzzz, tg_0_yzzzz, tg_0_zzz, tg_0_zzzz, tg_0_zzzzz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_0_xxxxx[i] = 4.0 * tg_0_xxx[i] * fxi[i] + tg_0_xxxx[i] * rb_x[i];

        tg_0_xxxxy[i] = tg_0_xxxx[i] * rb_y[i];

        tg_0_xxxxz[i] = tg_0_xxxx[i] * rb_z[i];

        tg_0_xxxyy[i] = 2.0 * tg_0_xyy[i] * fxi[i] + tg_0_xxyy[i] * rb_x[i];

        tg_0_xxxyz[i] = tg_0_xxxz[i] * rb_y[i];

        tg_0_xxxzz[i] = 2.0 * tg_0_xzz[i] * fxi[i] + tg_0_xxzz[i] * rb_x[i];

        tg_0_xxyyy[i] = tg_0_yyy[i] * fxi[i] + tg_0_xyyy[i] * rb_x[i];

        tg_0_xxyyz[i] = tg_0_xxyy[i] * rb_z[i];

        tg_0_xxyzz[i] = tg_0_xxzz[i] * rb_y[i];

        tg_0_xxzzz[i] = tg_0_zzz[i] * fxi[i] + tg_0_xzzz[i] * rb_x[i];

        tg_0_xyyyy[i] = tg_0_yyyy[i] * rb_x[i];

        tg_0_xyyyz[i] = tg_0_yyyz[i] * rb_x[i];

        tg_0_xyyzz[i] = tg_0_yyzz[i] * rb_x[i];

        tg_0_xyzzz[i] = tg_0_yzzz[i] * rb_x[i];

        tg_0_xzzzz[i] = tg_0_zzzz[i] * rb_x[i];

        tg_0_yyyyy[i] = 4.0 * tg_0_yyy[i] * fxi[i] + tg_0_yyyy[i] * rb_y[i];

        tg_0_yyyyz[i] = tg_0_yyyy[i] * rb_z[i];

        tg_0_yyyzz[i] = 2.0 * tg_0_yzz[i] * fxi[i] + tg_0_yyzz[i] * rb_y[i];

        tg_0_yyzzz[i] = tg_0_zzz[i] * fxi[i] + tg_0_yzzz[i] * rb_y[i];

        tg_0_yzzzz[i] = tg_0_zzzz[i] * rb_y[i];

        tg_0_zzzzz[i] = 4.0 * tg_0_zzz[i] * fxi[i] + tg_0_zzzz[i] * rb_z[i];
    }
}

} // t2lecp namespace


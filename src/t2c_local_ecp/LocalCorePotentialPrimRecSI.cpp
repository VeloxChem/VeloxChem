#include "LocalCorePotentialPrimRecSI.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_si(CSimdArray<double>& pbuffer, 
                                  const size_t idx_si,
                                  const size_t idx_sg,
                                  const size_t idx_sh,
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

    // Set up components of auxiliary buffer : SG

    auto tg_0_xxxx = pbuffer.data(idx_sg);

    auto tg_0_xxyy = pbuffer.data(idx_sg + 3);

    auto tg_0_xxzz = pbuffer.data(idx_sg + 5);

    auto tg_0_xyyy = pbuffer.data(idx_sg + 6);

    auto tg_0_xzzz = pbuffer.data(idx_sg + 9);

    auto tg_0_yyyy = pbuffer.data(idx_sg + 10);

    auto tg_0_yyzz = pbuffer.data(idx_sg + 12);

    auto tg_0_yzzz = pbuffer.data(idx_sg + 13);

    auto tg_0_zzzz = pbuffer.data(idx_sg + 14);

    // Set up components of auxiliary buffer : SH

    auto tg_0_xxxxx = pbuffer.data(idx_sh);

    auto tg_0_xxxxz = pbuffer.data(idx_sh + 2);

    auto tg_0_xxxyy = pbuffer.data(idx_sh + 3);

    auto tg_0_xxxzz = pbuffer.data(idx_sh + 5);

    auto tg_0_xxyyy = pbuffer.data(idx_sh + 6);

    auto tg_0_xxzzz = pbuffer.data(idx_sh + 9);

    auto tg_0_xyyyy = pbuffer.data(idx_sh + 10);

    auto tg_0_xyyzz = pbuffer.data(idx_sh + 12);

    auto tg_0_xzzzz = pbuffer.data(idx_sh + 14);

    auto tg_0_yyyyy = pbuffer.data(idx_sh + 15);

    auto tg_0_yyyyz = pbuffer.data(idx_sh + 16);

    auto tg_0_yyyzz = pbuffer.data(idx_sh + 17);

    auto tg_0_yyzzz = pbuffer.data(idx_sh + 18);

    auto tg_0_yzzzz = pbuffer.data(idx_sh + 19);

    auto tg_0_zzzzz = pbuffer.data(idx_sh + 20);

    // Set up components of targeted buffer : SI

    auto tg_0_xxxxxx = pbuffer.data(idx_si);

    auto tg_0_xxxxxy = pbuffer.data(idx_si + 1);

    auto tg_0_xxxxxz = pbuffer.data(idx_si + 2);

    auto tg_0_xxxxyy = pbuffer.data(idx_si + 3);

    auto tg_0_xxxxyz = pbuffer.data(idx_si + 4);

    auto tg_0_xxxxzz = pbuffer.data(idx_si + 5);

    auto tg_0_xxxyyy = pbuffer.data(idx_si + 6);

    auto tg_0_xxxyyz = pbuffer.data(idx_si + 7);

    auto tg_0_xxxyzz = pbuffer.data(idx_si + 8);

    auto tg_0_xxxzzz = pbuffer.data(idx_si + 9);

    auto tg_0_xxyyyy = pbuffer.data(idx_si + 10);

    auto tg_0_xxyyyz = pbuffer.data(idx_si + 11);

    auto tg_0_xxyyzz = pbuffer.data(idx_si + 12);

    auto tg_0_xxyzzz = pbuffer.data(idx_si + 13);

    auto tg_0_xxzzzz = pbuffer.data(idx_si + 14);

    auto tg_0_xyyyyy = pbuffer.data(idx_si + 15);

    auto tg_0_xyyyyz = pbuffer.data(idx_si + 16);

    auto tg_0_xyyyzz = pbuffer.data(idx_si + 17);

    auto tg_0_xyyzzz = pbuffer.data(idx_si + 18);

    auto tg_0_xyzzzz = pbuffer.data(idx_si + 19);

    auto tg_0_xzzzzz = pbuffer.data(idx_si + 20);

    auto tg_0_yyyyyy = pbuffer.data(idx_si + 21);

    auto tg_0_yyyyyz = pbuffer.data(idx_si + 22);

    auto tg_0_yyyyzz = pbuffer.data(idx_si + 23);

    auto tg_0_yyyzzz = pbuffer.data(idx_si + 24);

    auto tg_0_yyzzzz = pbuffer.data(idx_si + 25);

    auto tg_0_yzzzzz = pbuffer.data(idx_si + 26);

    auto tg_0_zzzzzz = pbuffer.data(idx_si + 27);

    #pragma omp simd aligned(fxi, rb_x, rb_y, rb_z, tg_0_xxxx, tg_0_xxxxx, tg_0_xxxxxx, tg_0_xxxxxy, tg_0_xxxxxz, tg_0_xxxxyy, tg_0_xxxxyz, tg_0_xxxxz, tg_0_xxxxzz, tg_0_xxxyy, tg_0_xxxyyy, tg_0_xxxyyz, tg_0_xxxyzz, tg_0_xxxzz, tg_0_xxxzzz, tg_0_xxyy, tg_0_xxyyy, tg_0_xxyyyy, tg_0_xxyyyz, tg_0_xxyyzz, tg_0_xxyzzz, tg_0_xxzz, tg_0_xxzzz, tg_0_xxzzzz, tg_0_xyyy, tg_0_xyyyy, tg_0_xyyyyy, tg_0_xyyyyz, tg_0_xyyyzz, tg_0_xyyzz, tg_0_xyyzzz, tg_0_xyzzzz, tg_0_xzzz, tg_0_xzzzz, tg_0_xzzzzz, tg_0_yyyy, tg_0_yyyyy, tg_0_yyyyyy, tg_0_yyyyyz, tg_0_yyyyz, tg_0_yyyyzz, tg_0_yyyzz, tg_0_yyyzzz, tg_0_yyzz, tg_0_yyzzz, tg_0_yyzzzz, tg_0_yzzz, tg_0_yzzzz, tg_0_yzzzzz, tg_0_zzzz, tg_0_zzzzz, tg_0_zzzzzz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_0_xxxxxx[i] = 5.0 * tg_0_xxxx[i] * fxi[i] + tg_0_xxxxx[i] * rb_x[i];

        tg_0_xxxxxy[i] = tg_0_xxxxx[i] * rb_y[i];

        tg_0_xxxxxz[i] = tg_0_xxxxx[i] * rb_z[i];

        tg_0_xxxxyy[i] = 3.0 * tg_0_xxyy[i] * fxi[i] + tg_0_xxxyy[i] * rb_x[i];

        tg_0_xxxxyz[i] = tg_0_xxxxz[i] * rb_y[i];

        tg_0_xxxxzz[i] = 3.0 * tg_0_xxzz[i] * fxi[i] + tg_0_xxxzz[i] * rb_x[i];

        tg_0_xxxyyy[i] = 2.0 * tg_0_xyyy[i] * fxi[i] + tg_0_xxyyy[i] * rb_x[i];

        tg_0_xxxyyz[i] = tg_0_xxxyy[i] * rb_z[i];

        tg_0_xxxyzz[i] = tg_0_xxxzz[i] * rb_y[i];

        tg_0_xxxzzz[i] = 2.0 * tg_0_xzzz[i] * fxi[i] + tg_0_xxzzz[i] * rb_x[i];

        tg_0_xxyyyy[i] = tg_0_yyyy[i] * fxi[i] + tg_0_xyyyy[i] * rb_x[i];

        tg_0_xxyyyz[i] = tg_0_xxyyy[i] * rb_z[i];

        tg_0_xxyyzz[i] = tg_0_yyzz[i] * fxi[i] + tg_0_xyyzz[i] * rb_x[i];

        tg_0_xxyzzz[i] = tg_0_xxzzz[i] * rb_y[i];

        tg_0_xxzzzz[i] = tg_0_zzzz[i] * fxi[i] + tg_0_xzzzz[i] * rb_x[i];

        tg_0_xyyyyy[i] = tg_0_yyyyy[i] * rb_x[i];

        tg_0_xyyyyz[i] = tg_0_yyyyz[i] * rb_x[i];

        tg_0_xyyyzz[i] = tg_0_yyyzz[i] * rb_x[i];

        tg_0_xyyzzz[i] = tg_0_yyzzz[i] * rb_x[i];

        tg_0_xyzzzz[i] = tg_0_yzzzz[i] * rb_x[i];

        tg_0_xzzzzz[i] = tg_0_zzzzz[i] * rb_x[i];

        tg_0_yyyyyy[i] = 5.0 * tg_0_yyyy[i] * fxi[i] + tg_0_yyyyy[i] * rb_y[i];

        tg_0_yyyyyz[i] = tg_0_yyyyy[i] * rb_z[i];

        tg_0_yyyyzz[i] = 3.0 * tg_0_yyzz[i] * fxi[i] + tg_0_yyyzz[i] * rb_y[i];

        tg_0_yyyzzz[i] = 2.0 * tg_0_yzzz[i] * fxi[i] + tg_0_yyzzz[i] * rb_y[i];

        tg_0_yyzzzz[i] = tg_0_zzzz[i] * fxi[i] + tg_0_yzzzz[i] * rb_y[i];

        tg_0_yzzzzz[i] = tg_0_zzzzz[i] * rb_y[i];

        tg_0_zzzzzz[i] = 5.0 * tg_0_zzzz[i] * fxi[i] + tg_0_zzzzz[i] * rb_z[i];
    }
}

} // t2lecp namespace


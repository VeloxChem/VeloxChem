#include "LocalCorePotentialPrimRecHS.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_hs(CSimdArray<double>& pbuffer, 
                                  const size_t idx_hs,
                                  const size_t idx_fs,
                                  const size_t idx_gs,
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

    // Set up components of auxiliary buffer : FS

    auto tg_xxx_0 = pbuffer.data(idx_fs);

    auto tg_xyy_0 = pbuffer.data(idx_fs + 3);

    auto tg_xzz_0 = pbuffer.data(idx_fs + 5);

    auto tg_yyy_0 = pbuffer.data(idx_fs + 6);

    auto tg_yzz_0 = pbuffer.data(idx_fs + 8);

    auto tg_zzz_0 = pbuffer.data(idx_fs + 9);

    // Set up components of auxiliary buffer : GS

    auto tg_xxxx_0 = pbuffer.data(idx_gs);

    auto tg_xxxz_0 = pbuffer.data(idx_gs + 2);

    auto tg_xxyy_0 = pbuffer.data(idx_gs + 3);

    auto tg_xxzz_0 = pbuffer.data(idx_gs + 5);

    auto tg_xyyy_0 = pbuffer.data(idx_gs + 6);

    auto tg_xzzz_0 = pbuffer.data(idx_gs + 9);

    auto tg_yyyy_0 = pbuffer.data(idx_gs + 10);

    auto tg_yyyz_0 = pbuffer.data(idx_gs + 11);

    auto tg_yyzz_0 = pbuffer.data(idx_gs + 12);

    auto tg_yzzz_0 = pbuffer.data(idx_gs + 13);

    auto tg_zzzz_0 = pbuffer.data(idx_gs + 14);

    // Set up components of targeted buffer : HS

    auto tg_xxxxx_0 = pbuffer.data(idx_hs);

    auto tg_xxxxy_0 = pbuffer.data(idx_hs + 1);

    auto tg_xxxxz_0 = pbuffer.data(idx_hs + 2);

    auto tg_xxxyy_0 = pbuffer.data(idx_hs + 3);

    auto tg_xxxyz_0 = pbuffer.data(idx_hs + 4);

    auto tg_xxxzz_0 = pbuffer.data(idx_hs + 5);

    auto tg_xxyyy_0 = pbuffer.data(idx_hs + 6);

    auto tg_xxyyz_0 = pbuffer.data(idx_hs + 7);

    auto tg_xxyzz_0 = pbuffer.data(idx_hs + 8);

    auto tg_xxzzz_0 = pbuffer.data(idx_hs + 9);

    auto tg_xyyyy_0 = pbuffer.data(idx_hs + 10);

    auto tg_xyyyz_0 = pbuffer.data(idx_hs + 11);

    auto tg_xyyzz_0 = pbuffer.data(idx_hs + 12);

    auto tg_xyzzz_0 = pbuffer.data(idx_hs + 13);

    auto tg_xzzzz_0 = pbuffer.data(idx_hs + 14);

    auto tg_yyyyy_0 = pbuffer.data(idx_hs + 15);

    auto tg_yyyyz_0 = pbuffer.data(idx_hs + 16);

    auto tg_yyyzz_0 = pbuffer.data(idx_hs + 17);

    auto tg_yyzzz_0 = pbuffer.data(idx_hs + 18);

    auto tg_yzzzz_0 = pbuffer.data(idx_hs + 19);

    auto tg_zzzzz_0 = pbuffer.data(idx_hs + 20);

    #pragma omp simd aligned(fxi, ra_x, ra_y, ra_z, tg_xxx_0, tg_xxxx_0, tg_xxxxx_0, tg_xxxxy_0, tg_xxxxz_0, tg_xxxyy_0, tg_xxxyz_0, tg_xxxz_0, tg_xxxzz_0, tg_xxyy_0, tg_xxyyy_0, tg_xxyyz_0, tg_xxyzz_0, tg_xxzz_0, tg_xxzzz_0, tg_xyy_0, tg_xyyy_0, tg_xyyyy_0, tg_xyyyz_0, tg_xyyzz_0, tg_xyzzz_0, tg_xzz_0, tg_xzzz_0, tg_xzzzz_0, tg_yyy_0, tg_yyyy_0, tg_yyyyy_0, tg_yyyyz_0, tg_yyyz_0, tg_yyyzz_0, tg_yyzz_0, tg_yyzzz_0, tg_yzz_0, tg_yzzz_0, tg_yzzzz_0, tg_zzz_0, tg_zzzz_0, tg_zzzzz_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_xxxxx_0[i] = 4.0 * tg_xxx_0[i] * fxi[i] + tg_xxxx_0[i] * ra_x[i];

        tg_xxxxy_0[i] = tg_xxxx_0[i] * ra_y[i];

        tg_xxxxz_0[i] = tg_xxxx_0[i] * ra_z[i];

        tg_xxxyy_0[i] = 2.0 * tg_xyy_0[i] * fxi[i] + tg_xxyy_0[i] * ra_x[i];

        tg_xxxyz_0[i] = tg_xxxz_0[i] * ra_y[i];

        tg_xxxzz_0[i] = 2.0 * tg_xzz_0[i] * fxi[i] + tg_xxzz_0[i] * ra_x[i];

        tg_xxyyy_0[i] = tg_yyy_0[i] * fxi[i] + tg_xyyy_0[i] * ra_x[i];

        tg_xxyyz_0[i] = tg_xxyy_0[i] * ra_z[i];

        tg_xxyzz_0[i] = tg_xxzz_0[i] * ra_y[i];

        tg_xxzzz_0[i] = tg_zzz_0[i] * fxi[i] + tg_xzzz_0[i] * ra_x[i];

        tg_xyyyy_0[i] = tg_yyyy_0[i] * ra_x[i];

        tg_xyyyz_0[i] = tg_yyyz_0[i] * ra_x[i];

        tg_xyyzz_0[i] = tg_yyzz_0[i] * ra_x[i];

        tg_xyzzz_0[i] = tg_yzzz_0[i] * ra_x[i];

        tg_xzzzz_0[i] = tg_zzzz_0[i] * ra_x[i];

        tg_yyyyy_0[i] = 4.0 * tg_yyy_0[i] * fxi[i] + tg_yyyy_0[i] * ra_y[i];

        tg_yyyyz_0[i] = tg_yyyy_0[i] * ra_z[i];

        tg_yyyzz_0[i] = 2.0 * tg_yzz_0[i] * fxi[i] + tg_yyzz_0[i] * ra_y[i];

        tg_yyzzz_0[i] = tg_zzz_0[i] * fxi[i] + tg_yzzz_0[i] * ra_y[i];

        tg_yzzzz_0[i] = tg_zzzz_0[i] * ra_y[i];

        tg_zzzzz_0[i] = 4.0 * tg_zzz_0[i] * fxi[i] + tg_zzzz_0[i] * ra_z[i];
    }
}

} // t2lecp namespace


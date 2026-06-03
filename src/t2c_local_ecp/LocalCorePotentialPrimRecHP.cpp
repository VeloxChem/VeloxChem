#include "LocalCorePotentialPrimRecHP.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_hp(CSimdArray<double>& pbuffer, 
                                  const size_t idx_hp,
                                  const size_t idx_fp,
                                  const size_t idx_gs,
                                  const size_t idx_gp,
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

    // Set up components of auxiliary buffer : FP

    auto tg_xxx_x = pbuffer.data(idx_fp);

    auto tg_xxx_y = pbuffer.data(idx_fp + 1);

    auto tg_xxx_z = pbuffer.data(idx_fp + 2);

    auto tg_xxy_x = pbuffer.data(idx_fp + 3);

    auto tg_xxy_y = pbuffer.data(idx_fp + 4);

    auto tg_xxz_x = pbuffer.data(idx_fp + 6);

    auto tg_xxz_z = pbuffer.data(idx_fp + 8);

    auto tg_xyy_y = pbuffer.data(idx_fp + 10);

    auto tg_xyy_z = pbuffer.data(idx_fp + 11);

    auto tg_xzz_y = pbuffer.data(idx_fp + 16);

    auto tg_xzz_z = pbuffer.data(idx_fp + 17);

    auto tg_yyy_x = pbuffer.data(idx_fp + 18);

    auto tg_yyy_y = pbuffer.data(idx_fp + 19);

    auto tg_yyy_z = pbuffer.data(idx_fp + 20);

    auto tg_yyz_y = pbuffer.data(idx_fp + 22);

    auto tg_yyz_z = pbuffer.data(idx_fp + 23);

    auto tg_yzz_x = pbuffer.data(idx_fp + 24);

    auto tg_yzz_y = pbuffer.data(idx_fp + 25);

    auto tg_yzz_z = pbuffer.data(idx_fp + 26);

    auto tg_zzz_x = pbuffer.data(idx_fp + 27);

    auto tg_zzz_y = pbuffer.data(idx_fp + 28);

    auto tg_zzz_z = pbuffer.data(idx_fp + 29);

    // Set up components of auxiliary buffer : GS

    auto tg_xxxx_0 = pbuffer.data(idx_gs);

    auto tg_yyyy_0 = pbuffer.data(idx_gs + 10);

    auto tg_yyzz_0 = pbuffer.data(idx_gs + 12);

    auto tg_zzzz_0 = pbuffer.data(idx_gs + 14);

    // Set up components of auxiliary buffer : GP

    auto tg_xxxx_x = pbuffer.data(idx_gp);

    auto tg_xxxx_y = pbuffer.data(idx_gp + 1);

    auto tg_xxxx_z = pbuffer.data(idx_gp + 2);

    auto tg_xxxy_x = pbuffer.data(idx_gp + 3);

    auto tg_xxxy_y = pbuffer.data(idx_gp + 4);

    auto tg_xxxz_x = pbuffer.data(idx_gp + 6);

    auto tg_xxxz_z = pbuffer.data(idx_gp + 8);

    auto tg_xxyy_x = pbuffer.data(idx_gp + 9);

    auto tg_xxyy_y = pbuffer.data(idx_gp + 10);

    auto tg_xxyy_z = pbuffer.data(idx_gp + 11);

    auto tg_xxzz_x = pbuffer.data(idx_gp + 15);

    auto tg_xxzz_y = pbuffer.data(idx_gp + 16);

    auto tg_xxzz_z = pbuffer.data(idx_gp + 17);

    auto tg_xyyy_x = pbuffer.data(idx_gp + 18);

    auto tg_xyyy_y = pbuffer.data(idx_gp + 19);

    auto tg_xyyy_z = pbuffer.data(idx_gp + 20);

    auto tg_xyyz_z = pbuffer.data(idx_gp + 23);

    auto tg_xyzz_y = pbuffer.data(idx_gp + 25);

    auto tg_xzzz_x = pbuffer.data(idx_gp + 27);

    auto tg_xzzz_y = pbuffer.data(idx_gp + 28);

    auto tg_xzzz_z = pbuffer.data(idx_gp + 29);

    auto tg_yyyy_x = pbuffer.data(idx_gp + 30);

    auto tg_yyyy_y = pbuffer.data(idx_gp + 31);

    auto tg_yyyy_z = pbuffer.data(idx_gp + 32);

    auto tg_yyyz_y = pbuffer.data(idx_gp + 34);

    auto tg_yyyz_z = pbuffer.data(idx_gp + 35);

    auto tg_yyzz_x = pbuffer.data(idx_gp + 36);

    auto tg_yyzz_y = pbuffer.data(idx_gp + 37);

    auto tg_yyzz_z = pbuffer.data(idx_gp + 38);

    auto tg_yzzz_x = pbuffer.data(idx_gp + 39);

    auto tg_yzzz_y = pbuffer.data(idx_gp + 40);

    auto tg_yzzz_z = pbuffer.data(idx_gp + 41);

    auto tg_zzzz_x = pbuffer.data(idx_gp + 42);

    auto tg_zzzz_y = pbuffer.data(idx_gp + 43);

    auto tg_zzzz_z = pbuffer.data(idx_gp + 44);

    // Set up components of targeted buffer : HP

    auto tg_xxxxx_x = pbuffer.data(idx_hp);

    auto tg_xxxxx_y = pbuffer.data(idx_hp + 1);

    auto tg_xxxxx_z = pbuffer.data(idx_hp + 2);

    auto tg_xxxxy_x = pbuffer.data(idx_hp + 3);

    auto tg_xxxxy_y = pbuffer.data(idx_hp + 4);

    auto tg_xxxxy_z = pbuffer.data(idx_hp + 5);

    auto tg_xxxxz_x = pbuffer.data(idx_hp + 6);

    auto tg_xxxxz_y = pbuffer.data(idx_hp + 7);

    auto tg_xxxxz_z = pbuffer.data(idx_hp + 8);

    auto tg_xxxyy_x = pbuffer.data(idx_hp + 9);

    auto tg_xxxyy_y = pbuffer.data(idx_hp + 10);

    auto tg_xxxyy_z = pbuffer.data(idx_hp + 11);

    auto tg_xxxyz_x = pbuffer.data(idx_hp + 12);

    auto tg_xxxyz_y = pbuffer.data(idx_hp + 13);

    auto tg_xxxyz_z = pbuffer.data(idx_hp + 14);

    auto tg_xxxzz_x = pbuffer.data(idx_hp + 15);

    auto tg_xxxzz_y = pbuffer.data(idx_hp + 16);

    auto tg_xxxzz_z = pbuffer.data(idx_hp + 17);

    auto tg_xxyyy_x = pbuffer.data(idx_hp + 18);

    auto tg_xxyyy_y = pbuffer.data(idx_hp + 19);

    auto tg_xxyyy_z = pbuffer.data(idx_hp + 20);

    auto tg_xxyyz_x = pbuffer.data(idx_hp + 21);

    auto tg_xxyyz_y = pbuffer.data(idx_hp + 22);

    auto tg_xxyyz_z = pbuffer.data(idx_hp + 23);

    auto tg_xxyzz_x = pbuffer.data(idx_hp + 24);

    auto tg_xxyzz_y = pbuffer.data(idx_hp + 25);

    auto tg_xxyzz_z = pbuffer.data(idx_hp + 26);

    auto tg_xxzzz_x = pbuffer.data(idx_hp + 27);

    auto tg_xxzzz_y = pbuffer.data(idx_hp + 28);

    auto tg_xxzzz_z = pbuffer.data(idx_hp + 29);

    auto tg_xyyyy_x = pbuffer.data(idx_hp + 30);

    auto tg_xyyyy_y = pbuffer.data(idx_hp + 31);

    auto tg_xyyyy_z = pbuffer.data(idx_hp + 32);

    auto tg_xyyyz_x = pbuffer.data(idx_hp + 33);

    auto tg_xyyyz_y = pbuffer.data(idx_hp + 34);

    auto tg_xyyyz_z = pbuffer.data(idx_hp + 35);

    auto tg_xyyzz_x = pbuffer.data(idx_hp + 36);

    auto tg_xyyzz_y = pbuffer.data(idx_hp + 37);

    auto tg_xyyzz_z = pbuffer.data(idx_hp + 38);

    auto tg_xyzzz_x = pbuffer.data(idx_hp + 39);

    auto tg_xyzzz_y = pbuffer.data(idx_hp + 40);

    auto tg_xyzzz_z = pbuffer.data(idx_hp + 41);

    auto tg_xzzzz_x = pbuffer.data(idx_hp + 42);

    auto tg_xzzzz_y = pbuffer.data(idx_hp + 43);

    auto tg_xzzzz_z = pbuffer.data(idx_hp + 44);

    auto tg_yyyyy_x = pbuffer.data(idx_hp + 45);

    auto tg_yyyyy_y = pbuffer.data(idx_hp + 46);

    auto tg_yyyyy_z = pbuffer.data(idx_hp + 47);

    auto tg_yyyyz_x = pbuffer.data(idx_hp + 48);

    auto tg_yyyyz_y = pbuffer.data(idx_hp + 49);

    auto tg_yyyyz_z = pbuffer.data(idx_hp + 50);

    auto tg_yyyzz_x = pbuffer.data(idx_hp + 51);

    auto tg_yyyzz_y = pbuffer.data(idx_hp + 52);

    auto tg_yyyzz_z = pbuffer.data(idx_hp + 53);

    auto tg_yyzzz_x = pbuffer.data(idx_hp + 54);

    auto tg_yyzzz_y = pbuffer.data(idx_hp + 55);

    auto tg_yyzzz_z = pbuffer.data(idx_hp + 56);

    auto tg_yzzzz_x = pbuffer.data(idx_hp + 57);

    auto tg_yzzzz_y = pbuffer.data(idx_hp + 58);

    auto tg_yzzzz_z = pbuffer.data(idx_hp + 59);

    auto tg_zzzzz_x = pbuffer.data(idx_hp + 60);

    auto tg_zzzzz_y = pbuffer.data(idx_hp + 61);

    auto tg_zzzzz_z = pbuffer.data(idx_hp + 62);

    #pragma omp simd aligned(fxi, ra_x, ra_y, ra_z, tg_xxx_x, tg_xxx_y, tg_xxx_z, tg_xxxx_0, tg_xxxx_x, tg_xxxx_y, tg_xxxx_z, tg_xxxxx_x, tg_xxxxx_y, tg_xxxxx_z, tg_xxxxy_x, tg_xxxxy_y, tg_xxxxy_z, tg_xxxxz_x, tg_xxxxz_y, tg_xxxxz_z, tg_xxxy_x, tg_xxxy_y, tg_xxxyy_x, tg_xxxyy_y, tg_xxxyy_z, tg_xxxyz_x, tg_xxxyz_y, tg_xxxyz_z, tg_xxxz_x, tg_xxxz_z, tg_xxxzz_x, tg_xxxzz_y, tg_xxxzz_z, tg_xxy_x, tg_xxy_y, tg_xxyy_x, tg_xxyy_y, tg_xxyy_z, tg_xxyyy_x, tg_xxyyy_y, tg_xxyyy_z, tg_xxyyz_x, tg_xxyyz_y, tg_xxyyz_z, tg_xxyzz_x, tg_xxyzz_y, tg_xxyzz_z, tg_xxz_x, tg_xxz_z, tg_xxzz_x, tg_xxzz_y, tg_xxzz_z, tg_xxzzz_x, tg_xxzzz_y, tg_xxzzz_z, tg_xyy_y, tg_xyy_z, tg_xyyy_x, tg_xyyy_y, tg_xyyy_z, tg_xyyyy_x, tg_xyyyy_y, tg_xyyyy_z, tg_xyyyz_x, tg_xyyyz_y, tg_xyyyz_z, tg_xyyz_z, tg_xyyzz_x, tg_xyyzz_y, tg_xyyzz_z, tg_xyzz_y, tg_xyzzz_x, tg_xyzzz_y, tg_xyzzz_z, tg_xzz_y, tg_xzz_z, tg_xzzz_x, tg_xzzz_y, tg_xzzz_z, tg_xzzzz_x, tg_xzzzz_y, tg_xzzzz_z, tg_yyy_x, tg_yyy_y, tg_yyy_z, tg_yyyy_0, tg_yyyy_x, tg_yyyy_y, tg_yyyy_z, tg_yyyyy_x, tg_yyyyy_y, tg_yyyyy_z, tg_yyyyz_x, tg_yyyyz_y, tg_yyyyz_z, tg_yyyz_y, tg_yyyz_z, tg_yyyzz_x, tg_yyyzz_y, tg_yyyzz_z, tg_yyz_y, tg_yyz_z, tg_yyzz_0, tg_yyzz_x, tg_yyzz_y, tg_yyzz_z, tg_yyzzz_x, tg_yyzzz_y, tg_yyzzz_z, tg_yzz_x, tg_yzz_y, tg_yzz_z, tg_yzzz_x, tg_yzzz_y, tg_yzzz_z, tg_yzzzz_x, tg_yzzzz_y, tg_yzzzz_z, tg_zzz_x, tg_zzz_y, tg_zzz_z, tg_zzzz_0, tg_zzzz_x, tg_zzzz_y, tg_zzzz_z, tg_zzzzz_x, tg_zzzzz_y, tg_zzzzz_z  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_xxxxx_x[i] = 4.0 * tg_xxx_x[i] * fxi[i] + tg_xxxx_0[i] * fxi[i] + tg_xxxx_x[i] * ra_x[i];

        tg_xxxxx_y[i] = 4.0 * tg_xxx_y[i] * fxi[i] + tg_xxxx_y[i] * ra_x[i];

        tg_xxxxx_z[i] = 4.0 * tg_xxx_z[i] * fxi[i] + tg_xxxx_z[i] * ra_x[i];

        tg_xxxxy_x[i] = tg_xxxx_x[i] * ra_y[i];

        tg_xxxxy_y[i] = 3.0 * tg_xxy_y[i] * fxi[i] + tg_xxxy_y[i] * ra_x[i];

        tg_xxxxy_z[i] = tg_xxxx_z[i] * ra_y[i];

        tg_xxxxz_x[i] = tg_xxxx_x[i] * ra_z[i];

        tg_xxxxz_y[i] = tg_xxxx_y[i] * ra_z[i];

        tg_xxxxz_z[i] = 3.0 * tg_xxz_z[i] * fxi[i] + tg_xxxz_z[i] * ra_x[i];

        tg_xxxyy_x[i] = tg_xxx_x[i] * fxi[i] + tg_xxxy_x[i] * ra_y[i];

        tg_xxxyy_y[i] = 2.0 * tg_xyy_y[i] * fxi[i] + tg_xxyy_y[i] * ra_x[i];

        tg_xxxyy_z[i] = 2.0 * tg_xyy_z[i] * fxi[i] + tg_xxyy_z[i] * ra_x[i];

        tg_xxxyz_x[i] = tg_xxxz_x[i] * ra_y[i];

        tg_xxxyz_y[i] = tg_xxxy_y[i] * ra_z[i];

        tg_xxxyz_z[i] = tg_xxxz_z[i] * ra_y[i];

        tg_xxxzz_x[i] = tg_xxx_x[i] * fxi[i] + tg_xxxz_x[i] * ra_z[i];

        tg_xxxzz_y[i] = 2.0 * tg_xzz_y[i] * fxi[i] + tg_xxzz_y[i] * ra_x[i];

        tg_xxxzz_z[i] = 2.0 * tg_xzz_z[i] * fxi[i] + tg_xxzz_z[i] * ra_x[i];

        tg_xxyyy_x[i] = 2.0 * tg_xxy_x[i] * fxi[i] + tg_xxyy_x[i] * ra_y[i];

        tg_xxyyy_y[i] = tg_yyy_y[i] * fxi[i] + tg_xyyy_y[i] * ra_x[i];

        tg_xxyyy_z[i] = tg_yyy_z[i] * fxi[i] + tg_xyyy_z[i] * ra_x[i];

        tg_xxyyz_x[i] = tg_xxyy_x[i] * ra_z[i];

        tg_xxyyz_y[i] = tg_xxyy_y[i] * ra_z[i];

        tg_xxyyz_z[i] = tg_yyz_z[i] * fxi[i] + tg_xyyz_z[i] * ra_x[i];

        tg_xxyzz_x[i] = tg_xxzz_x[i] * ra_y[i];

        tg_xxyzz_y[i] = tg_yzz_y[i] * fxi[i] + tg_xyzz_y[i] * ra_x[i];

        tg_xxyzz_z[i] = tg_xxzz_z[i] * ra_y[i];

        tg_xxzzz_x[i] = 2.0 * tg_xxz_x[i] * fxi[i] + tg_xxzz_x[i] * ra_z[i];

        tg_xxzzz_y[i] = tg_zzz_y[i] * fxi[i] + tg_xzzz_y[i] * ra_x[i];

        tg_xxzzz_z[i] = tg_zzz_z[i] * fxi[i] + tg_xzzz_z[i] * ra_x[i];

        tg_xyyyy_x[i] = tg_yyyy_0[i] * fxi[i] + tg_yyyy_x[i] * ra_x[i];

        tg_xyyyy_y[i] = tg_yyyy_y[i] * ra_x[i];

        tg_xyyyy_z[i] = tg_yyyy_z[i] * ra_x[i];

        tg_xyyyz_x[i] = tg_xyyy_x[i] * ra_z[i];

        tg_xyyyz_y[i] = tg_yyyz_y[i] * ra_x[i];

        tg_xyyyz_z[i] = tg_yyyz_z[i] * ra_x[i];

        tg_xyyzz_x[i] = tg_yyzz_0[i] * fxi[i] + tg_yyzz_x[i] * ra_x[i];

        tg_xyyzz_y[i] = tg_yyzz_y[i] * ra_x[i];

        tg_xyyzz_z[i] = tg_yyzz_z[i] * ra_x[i];

        tg_xyzzz_x[i] = tg_xzzz_x[i] * ra_y[i];

        tg_xyzzz_y[i] = tg_yzzz_y[i] * ra_x[i];

        tg_xyzzz_z[i] = tg_yzzz_z[i] * ra_x[i];

        tg_xzzzz_x[i] = tg_zzzz_0[i] * fxi[i] + tg_zzzz_x[i] * ra_x[i];

        tg_xzzzz_y[i] = tg_zzzz_y[i] * ra_x[i];

        tg_xzzzz_z[i] = tg_zzzz_z[i] * ra_x[i];

        tg_yyyyy_x[i] = 4.0 * tg_yyy_x[i] * fxi[i] + tg_yyyy_x[i] * ra_y[i];

        tg_yyyyy_y[i] = 4.0 * tg_yyy_y[i] * fxi[i] + tg_yyyy_0[i] * fxi[i] + tg_yyyy_y[i] * ra_y[i];

        tg_yyyyy_z[i] = 4.0 * tg_yyy_z[i] * fxi[i] + tg_yyyy_z[i] * ra_y[i];

        tg_yyyyz_x[i] = tg_yyyy_x[i] * ra_z[i];

        tg_yyyyz_y[i] = tg_yyyy_y[i] * ra_z[i];

        tg_yyyyz_z[i] = 3.0 * tg_yyz_z[i] * fxi[i] + tg_yyyz_z[i] * ra_y[i];

        tg_yyyzz_x[i] = 2.0 * tg_yzz_x[i] * fxi[i] + tg_yyzz_x[i] * ra_y[i];

        tg_yyyzz_y[i] = tg_yyy_y[i] * fxi[i] + tg_yyyz_y[i] * ra_z[i];

        tg_yyyzz_z[i] = 2.0 * tg_yzz_z[i] * fxi[i] + tg_yyzz_z[i] * ra_y[i];

        tg_yyzzz_x[i] = tg_zzz_x[i] * fxi[i] + tg_yzzz_x[i] * ra_y[i];

        tg_yyzzz_y[i] = 2.0 * tg_yyz_y[i] * fxi[i] + tg_yyzz_y[i] * ra_z[i];

        tg_yyzzz_z[i] = tg_zzz_z[i] * fxi[i] + tg_yzzz_z[i] * ra_y[i];

        tg_yzzzz_x[i] = tg_zzzz_x[i] * ra_y[i];

        tg_yzzzz_y[i] = tg_zzzz_0[i] * fxi[i] + tg_zzzz_y[i] * ra_y[i];

        tg_yzzzz_z[i] = tg_zzzz_z[i] * ra_y[i];

        tg_zzzzz_x[i] = 4.0 * tg_zzz_x[i] * fxi[i] + tg_zzzz_x[i] * ra_z[i];

        tg_zzzzz_y[i] = 4.0 * tg_zzz_y[i] * fxi[i] + tg_zzzz_y[i] * ra_z[i];

        tg_zzzzz_z[i] = 4.0 * tg_zzz_z[i] * fxi[i] + tg_zzzz_0[i] * fxi[i] + tg_zzzz_z[i] * ra_z[i];
    }
}

} // t2lecp namespace


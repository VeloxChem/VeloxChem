#include "LocalCorePotentialPrimRecIP.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_ip(CSimdArray<double>& pbuffer, 
                                  const size_t idx_ip,
                                  const size_t idx_gp,
                                  const size_t idx_hs,
                                  const size_t idx_hp,
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

    auto tg_xyyy_y = pbuffer.data(idx_gp + 19);

    auto tg_xyyy_z = pbuffer.data(idx_gp + 20);

    auto tg_xyyz_z = pbuffer.data(idx_gp + 23);

    auto tg_xyzz_y = pbuffer.data(idx_gp + 25);

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

    // Set up components of auxiliary buffer : HS

    auto tg_xxxxx_0 = pbuffer.data(idx_hs);

    auto tg_yyyyy_0 = pbuffer.data(idx_hs + 15);

    auto tg_yyyzz_0 = pbuffer.data(idx_hs + 17);

    auto tg_yyzzz_0 = pbuffer.data(idx_hs + 18);

    auto tg_zzzzz_0 = pbuffer.data(idx_hs + 20);

    // Set up components of auxiliary buffer : HP

    auto tg_xxxxx_x = pbuffer.data(idx_hp);

    auto tg_xxxxx_y = pbuffer.data(idx_hp + 1);

    auto tg_xxxxx_z = pbuffer.data(idx_hp + 2);

    auto tg_xxxxy_x = pbuffer.data(idx_hp + 3);

    auto tg_xxxxy_y = pbuffer.data(idx_hp + 4);

    auto tg_xxxxz_x = pbuffer.data(idx_hp + 6);

    auto tg_xxxxz_z = pbuffer.data(idx_hp + 8);

    auto tg_xxxyy_x = pbuffer.data(idx_hp + 9);

    auto tg_xxxyy_y = pbuffer.data(idx_hp + 10);

    auto tg_xxxyy_z = pbuffer.data(idx_hp + 11);

    auto tg_xxxzz_x = pbuffer.data(idx_hp + 15);

    auto tg_xxxzz_y = pbuffer.data(idx_hp + 16);

    auto tg_xxxzz_z = pbuffer.data(idx_hp + 17);

    auto tg_xxyyy_x = pbuffer.data(idx_hp + 18);

    auto tg_xxyyy_y = pbuffer.data(idx_hp + 19);

    auto tg_xxyyy_z = pbuffer.data(idx_hp + 20);

    auto tg_xxyyz_z = pbuffer.data(idx_hp + 23);

    auto tg_xxyzz_x = pbuffer.data(idx_hp + 24);

    auto tg_xxyzz_y = pbuffer.data(idx_hp + 25);

    auto tg_xxzzz_x = pbuffer.data(idx_hp + 27);

    auto tg_xxzzz_y = pbuffer.data(idx_hp + 28);

    auto tg_xxzzz_z = pbuffer.data(idx_hp + 29);

    auto tg_xyyyy_x = pbuffer.data(idx_hp + 30);

    auto tg_xyyyy_y = pbuffer.data(idx_hp + 31);

    auto tg_xyyyy_z = pbuffer.data(idx_hp + 32);

    auto tg_xyyyz_z = pbuffer.data(idx_hp + 35);

    auto tg_xyyzz_y = pbuffer.data(idx_hp + 37);

    auto tg_xyyzz_z = pbuffer.data(idx_hp + 38);

    auto tg_xyzzz_y = pbuffer.data(idx_hp + 40);

    auto tg_xzzzz_x = pbuffer.data(idx_hp + 42);

    auto tg_xzzzz_y = pbuffer.data(idx_hp + 43);

    auto tg_xzzzz_z = pbuffer.data(idx_hp + 44);

    auto tg_yyyyy_x = pbuffer.data(idx_hp + 45);

    auto tg_yyyyy_y = pbuffer.data(idx_hp + 46);

    auto tg_yyyyy_z = pbuffer.data(idx_hp + 47);

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

    // Set up components of targeted buffer : IP

    auto tg_xxxxxx_x = pbuffer.data(idx_ip);

    auto tg_xxxxxx_y = pbuffer.data(idx_ip + 1);

    auto tg_xxxxxx_z = pbuffer.data(idx_ip + 2);

    auto tg_xxxxxy_x = pbuffer.data(idx_ip + 3);

    auto tg_xxxxxy_y = pbuffer.data(idx_ip + 4);

    auto tg_xxxxxy_z = pbuffer.data(idx_ip + 5);

    auto tg_xxxxxz_x = pbuffer.data(idx_ip + 6);

    auto tg_xxxxxz_y = pbuffer.data(idx_ip + 7);

    auto tg_xxxxxz_z = pbuffer.data(idx_ip + 8);

    auto tg_xxxxyy_x = pbuffer.data(idx_ip + 9);

    auto tg_xxxxyy_y = pbuffer.data(idx_ip + 10);

    auto tg_xxxxyy_z = pbuffer.data(idx_ip + 11);

    auto tg_xxxxyz_x = pbuffer.data(idx_ip + 12);

    auto tg_xxxxyz_y = pbuffer.data(idx_ip + 13);

    auto tg_xxxxyz_z = pbuffer.data(idx_ip + 14);

    auto tg_xxxxzz_x = pbuffer.data(idx_ip + 15);

    auto tg_xxxxzz_y = pbuffer.data(idx_ip + 16);

    auto tg_xxxxzz_z = pbuffer.data(idx_ip + 17);

    auto tg_xxxyyy_x = pbuffer.data(idx_ip + 18);

    auto tg_xxxyyy_y = pbuffer.data(idx_ip + 19);

    auto tg_xxxyyy_z = pbuffer.data(idx_ip + 20);

    auto tg_xxxyyz_x = pbuffer.data(idx_ip + 21);

    auto tg_xxxyyz_y = pbuffer.data(idx_ip + 22);

    auto tg_xxxyyz_z = pbuffer.data(idx_ip + 23);

    auto tg_xxxyzz_x = pbuffer.data(idx_ip + 24);

    auto tg_xxxyzz_y = pbuffer.data(idx_ip + 25);

    auto tg_xxxyzz_z = pbuffer.data(idx_ip + 26);

    auto tg_xxxzzz_x = pbuffer.data(idx_ip + 27);

    auto tg_xxxzzz_y = pbuffer.data(idx_ip + 28);

    auto tg_xxxzzz_z = pbuffer.data(idx_ip + 29);

    auto tg_xxyyyy_x = pbuffer.data(idx_ip + 30);

    auto tg_xxyyyy_y = pbuffer.data(idx_ip + 31);

    auto tg_xxyyyy_z = pbuffer.data(idx_ip + 32);

    auto tg_xxyyyz_x = pbuffer.data(idx_ip + 33);

    auto tg_xxyyyz_y = pbuffer.data(idx_ip + 34);

    auto tg_xxyyyz_z = pbuffer.data(idx_ip + 35);

    auto tg_xxyyzz_x = pbuffer.data(idx_ip + 36);

    auto tg_xxyyzz_y = pbuffer.data(idx_ip + 37);

    auto tg_xxyyzz_z = pbuffer.data(idx_ip + 38);

    auto tg_xxyzzz_x = pbuffer.data(idx_ip + 39);

    auto tg_xxyzzz_y = pbuffer.data(idx_ip + 40);

    auto tg_xxyzzz_z = pbuffer.data(idx_ip + 41);

    auto tg_xxzzzz_x = pbuffer.data(idx_ip + 42);

    auto tg_xxzzzz_y = pbuffer.data(idx_ip + 43);

    auto tg_xxzzzz_z = pbuffer.data(idx_ip + 44);

    auto tg_xyyyyy_x = pbuffer.data(idx_ip + 45);

    auto tg_xyyyyy_y = pbuffer.data(idx_ip + 46);

    auto tg_xyyyyy_z = pbuffer.data(idx_ip + 47);

    auto tg_xyyyyz_x = pbuffer.data(idx_ip + 48);

    auto tg_xyyyyz_y = pbuffer.data(idx_ip + 49);

    auto tg_xyyyyz_z = pbuffer.data(idx_ip + 50);

    auto tg_xyyyzz_x = pbuffer.data(idx_ip + 51);

    auto tg_xyyyzz_y = pbuffer.data(idx_ip + 52);

    auto tg_xyyyzz_z = pbuffer.data(idx_ip + 53);

    auto tg_xyyzzz_x = pbuffer.data(idx_ip + 54);

    auto tg_xyyzzz_y = pbuffer.data(idx_ip + 55);

    auto tg_xyyzzz_z = pbuffer.data(idx_ip + 56);

    auto tg_xyzzzz_x = pbuffer.data(idx_ip + 57);

    auto tg_xyzzzz_y = pbuffer.data(idx_ip + 58);

    auto tg_xyzzzz_z = pbuffer.data(idx_ip + 59);

    auto tg_xzzzzz_x = pbuffer.data(idx_ip + 60);

    auto tg_xzzzzz_y = pbuffer.data(idx_ip + 61);

    auto tg_xzzzzz_z = pbuffer.data(idx_ip + 62);

    auto tg_yyyyyy_x = pbuffer.data(idx_ip + 63);

    auto tg_yyyyyy_y = pbuffer.data(idx_ip + 64);

    auto tg_yyyyyy_z = pbuffer.data(idx_ip + 65);

    auto tg_yyyyyz_x = pbuffer.data(idx_ip + 66);

    auto tg_yyyyyz_y = pbuffer.data(idx_ip + 67);

    auto tg_yyyyyz_z = pbuffer.data(idx_ip + 68);

    auto tg_yyyyzz_x = pbuffer.data(idx_ip + 69);

    auto tg_yyyyzz_y = pbuffer.data(idx_ip + 70);

    auto tg_yyyyzz_z = pbuffer.data(idx_ip + 71);

    auto tg_yyyzzz_x = pbuffer.data(idx_ip + 72);

    auto tg_yyyzzz_y = pbuffer.data(idx_ip + 73);

    auto tg_yyyzzz_z = pbuffer.data(idx_ip + 74);

    auto tg_yyzzzz_x = pbuffer.data(idx_ip + 75);

    auto tg_yyzzzz_y = pbuffer.data(idx_ip + 76);

    auto tg_yyzzzz_z = pbuffer.data(idx_ip + 77);

    auto tg_yzzzzz_x = pbuffer.data(idx_ip + 78);

    auto tg_yzzzzz_y = pbuffer.data(idx_ip + 79);

    auto tg_yzzzzz_z = pbuffer.data(idx_ip + 80);

    auto tg_zzzzzz_x = pbuffer.data(idx_ip + 81);

    auto tg_zzzzzz_y = pbuffer.data(idx_ip + 82);

    auto tg_zzzzzz_z = pbuffer.data(idx_ip + 83);

    #pragma omp simd aligned(fxi, ra_x, ra_y, ra_z, tg_xxxx_x, tg_xxxx_y, tg_xxxx_z, tg_xxxxx_0, tg_xxxxx_x, tg_xxxxx_y, tg_xxxxx_z, tg_xxxxxx_x, tg_xxxxxx_y, tg_xxxxxx_z, tg_xxxxxy_x, tg_xxxxxy_y, tg_xxxxxy_z, tg_xxxxxz_x, tg_xxxxxz_y, tg_xxxxxz_z, tg_xxxxy_x, tg_xxxxy_y, tg_xxxxyy_x, tg_xxxxyy_y, tg_xxxxyy_z, tg_xxxxyz_x, tg_xxxxyz_y, tg_xxxxyz_z, tg_xxxxz_x, tg_xxxxz_z, tg_xxxxzz_x, tg_xxxxzz_y, tg_xxxxzz_z, tg_xxxy_x, tg_xxxy_y, tg_xxxyy_x, tg_xxxyy_y, tg_xxxyy_z, tg_xxxyyy_x, tg_xxxyyy_y, tg_xxxyyy_z, tg_xxxyyz_x, tg_xxxyyz_y, tg_xxxyyz_z, tg_xxxyzz_x, tg_xxxyzz_y, tg_xxxyzz_z, tg_xxxz_x, tg_xxxz_z, tg_xxxzz_x, tg_xxxzz_y, tg_xxxzz_z, tg_xxxzzz_x, tg_xxxzzz_y, tg_xxxzzz_z, tg_xxyy_x, tg_xxyy_y, tg_xxyy_z, tg_xxyyy_x, tg_xxyyy_y, tg_xxyyy_z, tg_xxyyyy_x, tg_xxyyyy_y, tg_xxyyyy_z, tg_xxyyyz_x, tg_xxyyyz_y, tg_xxyyyz_z, tg_xxyyz_z, tg_xxyyzz_x, tg_xxyyzz_y, tg_xxyyzz_z, tg_xxyzz_x, tg_xxyzz_y, tg_xxyzzz_x, tg_xxyzzz_y, tg_xxyzzz_z, tg_xxzz_x, tg_xxzz_y, tg_xxzz_z, tg_xxzzz_x, tg_xxzzz_y, tg_xxzzz_z, tg_xxzzzz_x, tg_xxzzzz_y, tg_xxzzzz_z, tg_xyyy_y, tg_xyyy_z, tg_xyyyy_x, tg_xyyyy_y, tg_xyyyy_z, tg_xyyyyy_x, tg_xyyyyy_y, tg_xyyyyy_z, tg_xyyyyz_x, tg_xyyyyz_y, tg_xyyyyz_z, tg_xyyyz_z, tg_xyyyzz_x, tg_xyyyzz_y, tg_xyyyzz_z, tg_xyyz_z, tg_xyyzz_y, tg_xyyzz_z, tg_xyyzzz_x, tg_xyyzzz_y, tg_xyyzzz_z, tg_xyzz_y, tg_xyzzz_y, tg_xyzzzz_x, tg_xyzzzz_y, tg_xyzzzz_z, tg_xzzz_y, tg_xzzz_z, tg_xzzzz_x, tg_xzzzz_y, tg_xzzzz_z, tg_xzzzzz_x, tg_xzzzzz_y, tg_xzzzzz_z, tg_yyyy_x, tg_yyyy_y, tg_yyyy_z, tg_yyyyy_0, tg_yyyyy_x, tg_yyyyy_y, tg_yyyyy_z, tg_yyyyyy_x, tg_yyyyyy_y, tg_yyyyyy_z, tg_yyyyyz_x, tg_yyyyyz_y, tg_yyyyyz_z, tg_yyyyz_y, tg_yyyyz_z, tg_yyyyzz_x, tg_yyyyzz_y, tg_yyyyzz_z, tg_yyyz_y, tg_yyyz_z, tg_yyyzz_0, tg_yyyzz_x, tg_yyyzz_y, tg_yyyzz_z, tg_yyyzzz_x, tg_yyyzzz_y, tg_yyyzzz_z, tg_yyzz_x, tg_yyzz_y, tg_yyzz_z, tg_yyzzz_0, tg_yyzzz_x, tg_yyzzz_y, tg_yyzzz_z, tg_yyzzzz_x, tg_yyzzzz_y, tg_yyzzzz_z, tg_yzzz_x, tg_yzzz_y, tg_yzzz_z, tg_yzzzz_x, tg_yzzzz_y, tg_yzzzz_z, tg_yzzzzz_x, tg_yzzzzz_y, tg_yzzzzz_z, tg_zzzz_x, tg_zzzz_y, tg_zzzz_z, tg_zzzzz_0, tg_zzzzz_x, tg_zzzzz_y, tg_zzzzz_z, tg_zzzzzz_x, tg_zzzzzz_y, tg_zzzzzz_z  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_xxxxxx_x[i] = 5.0 * tg_xxxx_x[i] * fxi[i] + tg_xxxxx_0[i] * fxi[i] + tg_xxxxx_x[i] * ra_x[i];

        tg_xxxxxx_y[i] = 5.0 * tg_xxxx_y[i] * fxi[i] + tg_xxxxx_y[i] * ra_x[i];

        tg_xxxxxx_z[i] = 5.0 * tg_xxxx_z[i] * fxi[i] + tg_xxxxx_z[i] * ra_x[i];

        tg_xxxxxy_x[i] = tg_xxxxx_x[i] * ra_y[i];

        tg_xxxxxy_y[i] = 4.0 * tg_xxxy_y[i] * fxi[i] + tg_xxxxy_y[i] * ra_x[i];

        tg_xxxxxy_z[i] = tg_xxxxx_z[i] * ra_y[i];

        tg_xxxxxz_x[i] = tg_xxxxx_x[i] * ra_z[i];

        tg_xxxxxz_y[i] = tg_xxxxx_y[i] * ra_z[i];

        tg_xxxxxz_z[i] = 4.0 * tg_xxxz_z[i] * fxi[i] + tg_xxxxz_z[i] * ra_x[i];

        tg_xxxxyy_x[i] = tg_xxxx_x[i] * fxi[i] + tg_xxxxy_x[i] * ra_y[i];

        tg_xxxxyy_y[i] = 3.0 * tg_xxyy_y[i] * fxi[i] + tg_xxxyy_y[i] * ra_x[i];

        tg_xxxxyy_z[i] = 3.0 * tg_xxyy_z[i] * fxi[i] + tg_xxxyy_z[i] * ra_x[i];

        tg_xxxxyz_x[i] = tg_xxxxz_x[i] * ra_y[i];

        tg_xxxxyz_y[i] = tg_xxxxy_y[i] * ra_z[i];

        tg_xxxxyz_z[i] = tg_xxxxz_z[i] * ra_y[i];

        tg_xxxxzz_x[i] = tg_xxxx_x[i] * fxi[i] + tg_xxxxz_x[i] * ra_z[i];

        tg_xxxxzz_y[i] = 3.0 * tg_xxzz_y[i] * fxi[i] + tg_xxxzz_y[i] * ra_x[i];

        tg_xxxxzz_z[i] = 3.0 * tg_xxzz_z[i] * fxi[i] + tg_xxxzz_z[i] * ra_x[i];

        tg_xxxyyy_x[i] = 2.0 * tg_xxxy_x[i] * fxi[i] + tg_xxxyy_x[i] * ra_y[i];

        tg_xxxyyy_y[i] = 2.0 * tg_xyyy_y[i] * fxi[i] + tg_xxyyy_y[i] * ra_x[i];

        tg_xxxyyy_z[i] = 2.0 * tg_xyyy_z[i] * fxi[i] + tg_xxyyy_z[i] * ra_x[i];

        tg_xxxyyz_x[i] = tg_xxxyy_x[i] * ra_z[i];

        tg_xxxyyz_y[i] = tg_xxxyy_y[i] * ra_z[i];

        tg_xxxyyz_z[i] = 2.0 * tg_xyyz_z[i] * fxi[i] + tg_xxyyz_z[i] * ra_x[i];

        tg_xxxyzz_x[i] = tg_xxxzz_x[i] * ra_y[i];

        tg_xxxyzz_y[i] = 2.0 * tg_xyzz_y[i] * fxi[i] + tg_xxyzz_y[i] * ra_x[i];

        tg_xxxyzz_z[i] = tg_xxxzz_z[i] * ra_y[i];

        tg_xxxzzz_x[i] = 2.0 * tg_xxxz_x[i] * fxi[i] + tg_xxxzz_x[i] * ra_z[i];

        tg_xxxzzz_y[i] = 2.0 * tg_xzzz_y[i] * fxi[i] + tg_xxzzz_y[i] * ra_x[i];

        tg_xxxzzz_z[i] = 2.0 * tg_xzzz_z[i] * fxi[i] + tg_xxzzz_z[i] * ra_x[i];

        tg_xxyyyy_x[i] = 3.0 * tg_xxyy_x[i] * fxi[i] + tg_xxyyy_x[i] * ra_y[i];

        tg_xxyyyy_y[i] = tg_yyyy_y[i] * fxi[i] + tg_xyyyy_y[i] * ra_x[i];

        tg_xxyyyy_z[i] = tg_yyyy_z[i] * fxi[i] + tg_xyyyy_z[i] * ra_x[i];

        tg_xxyyyz_x[i] = tg_xxyyy_x[i] * ra_z[i];

        tg_xxyyyz_y[i] = tg_xxyyy_y[i] * ra_z[i];

        tg_xxyyyz_z[i] = tg_yyyz_z[i] * fxi[i] + tg_xyyyz_z[i] * ra_x[i];

        tg_xxyyzz_x[i] = tg_xxzz_x[i] * fxi[i] + tg_xxyzz_x[i] * ra_y[i];

        tg_xxyyzz_y[i] = tg_yyzz_y[i] * fxi[i] + tg_xyyzz_y[i] * ra_x[i];

        tg_xxyyzz_z[i] = tg_yyzz_z[i] * fxi[i] + tg_xyyzz_z[i] * ra_x[i];

        tg_xxyzzz_x[i] = tg_xxzzz_x[i] * ra_y[i];

        tg_xxyzzz_y[i] = tg_yzzz_y[i] * fxi[i] + tg_xyzzz_y[i] * ra_x[i];

        tg_xxyzzz_z[i] = tg_xxzzz_z[i] * ra_y[i];

        tg_xxzzzz_x[i] = 3.0 * tg_xxzz_x[i] * fxi[i] + tg_xxzzz_x[i] * ra_z[i];

        tg_xxzzzz_y[i] = tg_zzzz_y[i] * fxi[i] + tg_xzzzz_y[i] * ra_x[i];

        tg_xxzzzz_z[i] = tg_zzzz_z[i] * fxi[i] + tg_xzzzz_z[i] * ra_x[i];

        tg_xyyyyy_x[i] = tg_yyyyy_0[i] * fxi[i] + tg_yyyyy_x[i] * ra_x[i];

        tg_xyyyyy_y[i] = tg_yyyyy_y[i] * ra_x[i];

        tg_xyyyyy_z[i] = tg_yyyyy_z[i] * ra_x[i];

        tg_xyyyyz_x[i] = tg_xyyyy_x[i] * ra_z[i];

        tg_xyyyyz_y[i] = tg_yyyyz_y[i] * ra_x[i];

        tg_xyyyyz_z[i] = tg_yyyyz_z[i] * ra_x[i];

        tg_xyyyzz_x[i] = tg_yyyzz_0[i] * fxi[i] + tg_yyyzz_x[i] * ra_x[i];

        tg_xyyyzz_y[i] = tg_yyyzz_y[i] * ra_x[i];

        tg_xyyyzz_z[i] = tg_yyyzz_z[i] * ra_x[i];

        tg_xyyzzz_x[i] = tg_yyzzz_0[i] * fxi[i] + tg_yyzzz_x[i] * ra_x[i];

        tg_xyyzzz_y[i] = tg_yyzzz_y[i] * ra_x[i];

        tg_xyyzzz_z[i] = tg_yyzzz_z[i] * ra_x[i];

        tg_xyzzzz_x[i] = tg_xzzzz_x[i] * ra_y[i];

        tg_xyzzzz_y[i] = tg_yzzzz_y[i] * ra_x[i];

        tg_xyzzzz_z[i] = tg_yzzzz_z[i] * ra_x[i];

        tg_xzzzzz_x[i] = tg_zzzzz_0[i] * fxi[i] + tg_zzzzz_x[i] * ra_x[i];

        tg_xzzzzz_y[i] = tg_zzzzz_y[i] * ra_x[i];

        tg_xzzzzz_z[i] = tg_zzzzz_z[i] * ra_x[i];

        tg_yyyyyy_x[i] = 5.0 * tg_yyyy_x[i] * fxi[i] + tg_yyyyy_x[i] * ra_y[i];

        tg_yyyyyy_y[i] = 5.0 * tg_yyyy_y[i] * fxi[i] + tg_yyyyy_0[i] * fxi[i] + tg_yyyyy_y[i] * ra_y[i];

        tg_yyyyyy_z[i] = 5.0 * tg_yyyy_z[i] * fxi[i] + tg_yyyyy_z[i] * ra_y[i];

        tg_yyyyyz_x[i] = tg_yyyyy_x[i] * ra_z[i];

        tg_yyyyyz_y[i] = tg_yyyyy_y[i] * ra_z[i];

        tg_yyyyyz_z[i] = 4.0 * tg_yyyz_z[i] * fxi[i] + tg_yyyyz_z[i] * ra_y[i];

        tg_yyyyzz_x[i] = 3.0 * tg_yyzz_x[i] * fxi[i] + tg_yyyzz_x[i] * ra_y[i];

        tg_yyyyzz_y[i] = tg_yyyy_y[i] * fxi[i] + tg_yyyyz_y[i] * ra_z[i];

        tg_yyyyzz_z[i] = 3.0 * tg_yyzz_z[i] * fxi[i] + tg_yyyzz_z[i] * ra_y[i];

        tg_yyyzzz_x[i] = 2.0 * tg_yzzz_x[i] * fxi[i] + tg_yyzzz_x[i] * ra_y[i];

        tg_yyyzzz_y[i] = 2.0 * tg_yyyz_y[i] * fxi[i] + tg_yyyzz_y[i] * ra_z[i];

        tg_yyyzzz_z[i] = 2.0 * tg_yzzz_z[i] * fxi[i] + tg_yyzzz_z[i] * ra_y[i];

        tg_yyzzzz_x[i] = tg_zzzz_x[i] * fxi[i] + tg_yzzzz_x[i] * ra_y[i];

        tg_yyzzzz_y[i] = 3.0 * tg_yyzz_y[i] * fxi[i] + tg_yyzzz_y[i] * ra_z[i];

        tg_yyzzzz_z[i] = tg_zzzz_z[i] * fxi[i] + tg_yzzzz_z[i] * ra_y[i];

        tg_yzzzzz_x[i] = tg_zzzzz_x[i] * ra_y[i];

        tg_yzzzzz_y[i] = tg_zzzzz_0[i] * fxi[i] + tg_zzzzz_y[i] * ra_y[i];

        tg_yzzzzz_z[i] = tg_zzzzz_z[i] * ra_y[i];

        tg_zzzzzz_x[i] = 5.0 * tg_zzzz_x[i] * fxi[i] + tg_zzzzz_x[i] * ra_z[i];

        tg_zzzzzz_y[i] = 5.0 * tg_zzzz_y[i] * fxi[i] + tg_zzzzz_y[i] * ra_z[i];

        tg_zzzzzz_z[i] = 5.0 * tg_zzzz_z[i] * fxi[i] + tg_zzzzz_0[i] * fxi[i] + tg_zzzzz_z[i] * ra_z[i];
    }
}

} // t2lecp namespace


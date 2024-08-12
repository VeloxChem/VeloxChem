#include "ElectricDipoleMomentumPrimRecIS.hpp"

namespace diprec { // diprec namespace

auto
comp_prim_electric_dipole_momentum_is(CSimdArray<double>& pbuffer, 
                                      const size_t idx_dip_is,
                                      const size_t idx_dip_gs,
                                      const size_t idx_ovl_hs,
                                      const size_t idx_dip_hs,
                                      const CSimdArray<double>& factors,
                                      const size_t idx_rpa,
                                      const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : GS

    auto tr_x_xxxx_0 = pbuffer.data(idx_dip_gs);

    auto tr_x_xxxy_0 = pbuffer.data(idx_dip_gs + 1);

    auto tr_x_xxxz_0 = pbuffer.data(idx_dip_gs + 2);

    auto tr_x_xxyy_0 = pbuffer.data(idx_dip_gs + 3);

    auto tr_x_xxzz_0 = pbuffer.data(idx_dip_gs + 5);

    auto tr_x_yyyy_0 = pbuffer.data(idx_dip_gs + 10);

    auto tr_x_yyzz_0 = pbuffer.data(idx_dip_gs + 12);

    auto tr_x_yzzz_0 = pbuffer.data(idx_dip_gs + 13);

    auto tr_x_zzzz_0 = pbuffer.data(idx_dip_gs + 14);

    auto tr_y_xxxx_0 = pbuffer.data(idx_dip_gs + 15);

    auto tr_y_xxxy_0 = pbuffer.data(idx_dip_gs + 16);

    auto tr_y_xxyy_0 = pbuffer.data(idx_dip_gs + 18);

    auto tr_y_xxzz_0 = pbuffer.data(idx_dip_gs + 20);

    auto tr_y_xyyy_0 = pbuffer.data(idx_dip_gs + 21);

    auto tr_y_xyzz_0 = pbuffer.data(idx_dip_gs + 23);

    auto tr_y_xzzz_0 = pbuffer.data(idx_dip_gs + 24);

    auto tr_y_yyyy_0 = pbuffer.data(idx_dip_gs + 25);

    auto tr_y_yyyz_0 = pbuffer.data(idx_dip_gs + 26);

    auto tr_y_yyzz_0 = pbuffer.data(idx_dip_gs + 27);

    auto tr_y_yzzz_0 = pbuffer.data(idx_dip_gs + 28);

    auto tr_y_zzzz_0 = pbuffer.data(idx_dip_gs + 29);

    auto tr_z_xxxx_0 = pbuffer.data(idx_dip_gs + 30);

    auto tr_z_xxxz_0 = pbuffer.data(idx_dip_gs + 32);

    auto tr_z_xxyy_0 = pbuffer.data(idx_dip_gs + 33);

    auto tr_z_xxzz_0 = pbuffer.data(idx_dip_gs + 35);

    auto tr_z_xyyy_0 = pbuffer.data(idx_dip_gs + 36);

    auto tr_z_xyyz_0 = pbuffer.data(idx_dip_gs + 37);

    auto tr_z_xzzz_0 = pbuffer.data(idx_dip_gs + 39);

    auto tr_z_yyyy_0 = pbuffer.data(idx_dip_gs + 40);

    auto tr_z_yyyz_0 = pbuffer.data(idx_dip_gs + 41);

    auto tr_z_yyzz_0 = pbuffer.data(idx_dip_gs + 42);

    auto tr_z_yzzz_0 = pbuffer.data(idx_dip_gs + 43);

    auto tr_z_zzzz_0 = pbuffer.data(idx_dip_gs + 44);

    // Set up components of auxiliary buffer : HS

    auto ts_xxxxx_0 = pbuffer.data(idx_ovl_hs);

    auto ts_yyyyy_0 = pbuffer.data(idx_ovl_hs + 15);

    auto ts_yyyzz_0 = pbuffer.data(idx_ovl_hs + 17);

    auto ts_yyzzz_0 = pbuffer.data(idx_ovl_hs + 18);

    auto ts_zzzzz_0 = pbuffer.data(idx_ovl_hs + 20);

    // Set up components of auxiliary buffer : HS

    auto tr_x_xxxxx_0 = pbuffer.data(idx_dip_hs);

    auto tr_x_xxxxy_0 = pbuffer.data(idx_dip_hs + 1);

    auto tr_x_xxxxz_0 = pbuffer.data(idx_dip_hs + 2);

    auto tr_x_xxxyy_0 = pbuffer.data(idx_dip_hs + 3);

    auto tr_x_xxxzz_0 = pbuffer.data(idx_dip_hs + 5);

    auto tr_x_xxyyy_0 = pbuffer.data(idx_dip_hs + 6);

    auto tr_x_xxyzz_0 = pbuffer.data(idx_dip_hs + 8);

    auto tr_x_xxzzz_0 = pbuffer.data(idx_dip_hs + 9);

    auto tr_x_xyyyy_0 = pbuffer.data(idx_dip_hs + 10);

    auto tr_x_xzzzz_0 = pbuffer.data(idx_dip_hs + 14);

    auto tr_x_yyyyy_0 = pbuffer.data(idx_dip_hs + 15);

    auto tr_x_yyyzz_0 = pbuffer.data(idx_dip_hs + 17);

    auto tr_x_yyzzz_0 = pbuffer.data(idx_dip_hs + 18);

    auto tr_x_yzzzz_0 = pbuffer.data(idx_dip_hs + 19);

    auto tr_x_zzzzz_0 = pbuffer.data(idx_dip_hs + 20);

    auto tr_y_xxxxx_0 = pbuffer.data(idx_dip_hs + 21);

    auto tr_y_xxxxy_0 = pbuffer.data(idx_dip_hs + 22);

    auto tr_y_xxxyy_0 = pbuffer.data(idx_dip_hs + 24);

    auto tr_y_xxxzz_0 = pbuffer.data(idx_dip_hs + 26);

    auto tr_y_xxyyy_0 = pbuffer.data(idx_dip_hs + 27);

    auto tr_y_xxyzz_0 = pbuffer.data(idx_dip_hs + 29);

    auto tr_y_xxzzz_0 = pbuffer.data(idx_dip_hs + 30);

    auto tr_y_xyyyy_0 = pbuffer.data(idx_dip_hs + 31);

    auto tr_y_xyyzz_0 = pbuffer.data(idx_dip_hs + 33);

    auto tr_y_xyzzz_0 = pbuffer.data(idx_dip_hs + 34);

    auto tr_y_xzzzz_0 = pbuffer.data(idx_dip_hs + 35);

    auto tr_y_yyyyy_0 = pbuffer.data(idx_dip_hs + 36);

    auto tr_y_yyyyz_0 = pbuffer.data(idx_dip_hs + 37);

    auto tr_y_yyyzz_0 = pbuffer.data(idx_dip_hs + 38);

    auto tr_y_yyzzz_0 = pbuffer.data(idx_dip_hs + 39);

    auto tr_y_yzzzz_0 = pbuffer.data(idx_dip_hs + 40);

    auto tr_y_zzzzz_0 = pbuffer.data(idx_dip_hs + 41);

    auto tr_z_xxxxx_0 = pbuffer.data(idx_dip_hs + 42);

    auto tr_z_xxxxz_0 = pbuffer.data(idx_dip_hs + 44);

    auto tr_z_xxxyy_0 = pbuffer.data(idx_dip_hs + 45);

    auto tr_z_xxxzz_0 = pbuffer.data(idx_dip_hs + 47);

    auto tr_z_xxyyy_0 = pbuffer.data(idx_dip_hs + 48);

    auto tr_z_xxyyz_0 = pbuffer.data(idx_dip_hs + 49);

    auto tr_z_xxzzz_0 = pbuffer.data(idx_dip_hs + 51);

    auto tr_z_xyyyy_0 = pbuffer.data(idx_dip_hs + 52);

    auto tr_z_xyyyz_0 = pbuffer.data(idx_dip_hs + 53);

    auto tr_z_xyyzz_0 = pbuffer.data(idx_dip_hs + 54);

    auto tr_z_xzzzz_0 = pbuffer.data(idx_dip_hs + 56);

    auto tr_z_yyyyy_0 = pbuffer.data(idx_dip_hs + 57);

    auto tr_z_yyyyz_0 = pbuffer.data(idx_dip_hs + 58);

    auto tr_z_yyyzz_0 = pbuffer.data(idx_dip_hs + 59);

    auto tr_z_yyzzz_0 = pbuffer.data(idx_dip_hs + 60);

    auto tr_z_yzzzz_0 = pbuffer.data(idx_dip_hs + 61);

    auto tr_z_zzzzz_0 = pbuffer.data(idx_dip_hs + 62);

    // Set up components of targeted buffer : IS

    auto tr_x_xxxxxx_0 = pbuffer.data(idx_dip_is);

    auto tr_x_xxxxxy_0 = pbuffer.data(idx_dip_is + 1);

    auto tr_x_xxxxxz_0 = pbuffer.data(idx_dip_is + 2);

    auto tr_x_xxxxyy_0 = pbuffer.data(idx_dip_is + 3);

    auto tr_x_xxxxyz_0 = pbuffer.data(idx_dip_is + 4);

    auto tr_x_xxxxzz_0 = pbuffer.data(idx_dip_is + 5);

    auto tr_x_xxxyyy_0 = pbuffer.data(idx_dip_is + 6);

    auto tr_x_xxxyyz_0 = pbuffer.data(idx_dip_is + 7);

    auto tr_x_xxxyzz_0 = pbuffer.data(idx_dip_is + 8);

    auto tr_x_xxxzzz_0 = pbuffer.data(idx_dip_is + 9);

    auto tr_x_xxyyyy_0 = pbuffer.data(idx_dip_is + 10);

    auto tr_x_xxyyyz_0 = pbuffer.data(idx_dip_is + 11);

    auto tr_x_xxyyzz_0 = pbuffer.data(idx_dip_is + 12);

    auto tr_x_xxyzzz_0 = pbuffer.data(idx_dip_is + 13);

    auto tr_x_xxzzzz_0 = pbuffer.data(idx_dip_is + 14);

    auto tr_x_xyyyyy_0 = pbuffer.data(idx_dip_is + 15);

    auto tr_x_xyyyyz_0 = pbuffer.data(idx_dip_is + 16);

    auto tr_x_xyyyzz_0 = pbuffer.data(idx_dip_is + 17);

    auto tr_x_xyyzzz_0 = pbuffer.data(idx_dip_is + 18);

    auto tr_x_xyzzzz_0 = pbuffer.data(idx_dip_is + 19);

    auto tr_x_xzzzzz_0 = pbuffer.data(idx_dip_is + 20);

    auto tr_x_yyyyyy_0 = pbuffer.data(idx_dip_is + 21);

    auto tr_x_yyyyyz_0 = pbuffer.data(idx_dip_is + 22);

    auto tr_x_yyyyzz_0 = pbuffer.data(idx_dip_is + 23);

    auto tr_x_yyyzzz_0 = pbuffer.data(idx_dip_is + 24);

    auto tr_x_yyzzzz_0 = pbuffer.data(idx_dip_is + 25);

    auto tr_x_yzzzzz_0 = pbuffer.data(idx_dip_is + 26);

    auto tr_x_zzzzzz_0 = pbuffer.data(idx_dip_is + 27);

    auto tr_y_xxxxxx_0 = pbuffer.data(idx_dip_is + 28);

    auto tr_y_xxxxxy_0 = pbuffer.data(idx_dip_is + 29);

    auto tr_y_xxxxxz_0 = pbuffer.data(idx_dip_is + 30);

    auto tr_y_xxxxyy_0 = pbuffer.data(idx_dip_is + 31);

    auto tr_y_xxxxyz_0 = pbuffer.data(idx_dip_is + 32);

    auto tr_y_xxxxzz_0 = pbuffer.data(idx_dip_is + 33);

    auto tr_y_xxxyyy_0 = pbuffer.data(idx_dip_is + 34);

    auto tr_y_xxxyyz_0 = pbuffer.data(idx_dip_is + 35);

    auto tr_y_xxxyzz_0 = pbuffer.data(idx_dip_is + 36);

    auto tr_y_xxxzzz_0 = pbuffer.data(idx_dip_is + 37);

    auto tr_y_xxyyyy_0 = pbuffer.data(idx_dip_is + 38);

    auto tr_y_xxyyyz_0 = pbuffer.data(idx_dip_is + 39);

    auto tr_y_xxyyzz_0 = pbuffer.data(idx_dip_is + 40);

    auto tr_y_xxyzzz_0 = pbuffer.data(idx_dip_is + 41);

    auto tr_y_xxzzzz_0 = pbuffer.data(idx_dip_is + 42);

    auto tr_y_xyyyyy_0 = pbuffer.data(idx_dip_is + 43);

    auto tr_y_xyyyyz_0 = pbuffer.data(idx_dip_is + 44);

    auto tr_y_xyyyzz_0 = pbuffer.data(idx_dip_is + 45);

    auto tr_y_xyyzzz_0 = pbuffer.data(idx_dip_is + 46);

    auto tr_y_xyzzzz_0 = pbuffer.data(idx_dip_is + 47);

    auto tr_y_xzzzzz_0 = pbuffer.data(idx_dip_is + 48);

    auto tr_y_yyyyyy_0 = pbuffer.data(idx_dip_is + 49);

    auto tr_y_yyyyyz_0 = pbuffer.data(idx_dip_is + 50);

    auto tr_y_yyyyzz_0 = pbuffer.data(idx_dip_is + 51);

    auto tr_y_yyyzzz_0 = pbuffer.data(idx_dip_is + 52);

    auto tr_y_yyzzzz_0 = pbuffer.data(idx_dip_is + 53);

    auto tr_y_yzzzzz_0 = pbuffer.data(idx_dip_is + 54);

    auto tr_y_zzzzzz_0 = pbuffer.data(idx_dip_is + 55);

    auto tr_z_xxxxxx_0 = pbuffer.data(idx_dip_is + 56);

    auto tr_z_xxxxxy_0 = pbuffer.data(idx_dip_is + 57);

    auto tr_z_xxxxxz_0 = pbuffer.data(idx_dip_is + 58);

    auto tr_z_xxxxyy_0 = pbuffer.data(idx_dip_is + 59);

    auto tr_z_xxxxyz_0 = pbuffer.data(idx_dip_is + 60);

    auto tr_z_xxxxzz_0 = pbuffer.data(idx_dip_is + 61);

    auto tr_z_xxxyyy_0 = pbuffer.data(idx_dip_is + 62);

    auto tr_z_xxxyyz_0 = pbuffer.data(idx_dip_is + 63);

    auto tr_z_xxxyzz_0 = pbuffer.data(idx_dip_is + 64);

    auto tr_z_xxxzzz_0 = pbuffer.data(idx_dip_is + 65);

    auto tr_z_xxyyyy_0 = pbuffer.data(idx_dip_is + 66);

    auto tr_z_xxyyyz_0 = pbuffer.data(idx_dip_is + 67);

    auto tr_z_xxyyzz_0 = pbuffer.data(idx_dip_is + 68);

    auto tr_z_xxyzzz_0 = pbuffer.data(idx_dip_is + 69);

    auto tr_z_xxzzzz_0 = pbuffer.data(idx_dip_is + 70);

    auto tr_z_xyyyyy_0 = pbuffer.data(idx_dip_is + 71);

    auto tr_z_xyyyyz_0 = pbuffer.data(idx_dip_is + 72);

    auto tr_z_xyyyzz_0 = pbuffer.data(idx_dip_is + 73);

    auto tr_z_xyyzzz_0 = pbuffer.data(idx_dip_is + 74);

    auto tr_z_xyzzzz_0 = pbuffer.data(idx_dip_is + 75);

    auto tr_z_xzzzzz_0 = pbuffer.data(idx_dip_is + 76);

    auto tr_z_yyyyyy_0 = pbuffer.data(idx_dip_is + 77);

    auto tr_z_yyyyyz_0 = pbuffer.data(idx_dip_is + 78);

    auto tr_z_yyyyzz_0 = pbuffer.data(idx_dip_is + 79);

    auto tr_z_yyyzzz_0 = pbuffer.data(idx_dip_is + 80);

    auto tr_z_yyzzzz_0 = pbuffer.data(idx_dip_is + 81);

    auto tr_z_yzzzzz_0 = pbuffer.data(idx_dip_is + 82);

    auto tr_z_zzzzzz_0 = pbuffer.data(idx_dip_is + 83);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_x_xxxx_0, tr_x_xxxxx_0, tr_x_xxxxxx_0, tr_x_xxxxxy_0, tr_x_xxxxxz_0, tr_x_xxxxy_0, tr_x_xxxxyy_0, tr_x_xxxxyz_0, tr_x_xxxxz_0, tr_x_xxxxzz_0, tr_x_xxxy_0, tr_x_xxxyy_0, tr_x_xxxyyy_0, tr_x_xxxyyz_0, tr_x_xxxyzz_0, tr_x_xxxz_0, tr_x_xxxzz_0, tr_x_xxxzzz_0, tr_x_xxyy_0, tr_x_xxyyy_0, tr_x_xxyyyy_0, tr_x_xxyyyz_0, tr_x_xxyyzz_0, tr_x_xxyzz_0, tr_x_xxyzzz_0, tr_x_xxzz_0, tr_x_xxzzz_0, tr_x_xxzzzz_0, tr_x_xyyyy_0, tr_x_xyyyyy_0, tr_x_xyyyyz_0, tr_x_xyyyzz_0, tr_x_xyyzzz_0, tr_x_xyzzzz_0, tr_x_xzzzz_0, tr_x_xzzzzz_0, tr_x_yyyy_0, tr_x_yyyyy_0, tr_x_yyyyyy_0, tr_x_yyyyyz_0, tr_x_yyyyzz_0, tr_x_yyyzz_0, tr_x_yyyzzz_0, tr_x_yyzz_0, tr_x_yyzzz_0, tr_x_yyzzzz_0, tr_x_yzzz_0, tr_x_yzzzz_0, tr_x_yzzzzz_0, tr_x_zzzz_0, tr_x_zzzzz_0, tr_x_zzzzzz_0, tr_y_xxxx_0, tr_y_xxxxx_0, tr_y_xxxxxx_0, tr_y_xxxxxy_0, tr_y_xxxxxz_0, tr_y_xxxxy_0, tr_y_xxxxyy_0, tr_y_xxxxyz_0, tr_y_xxxxzz_0, tr_y_xxxy_0, tr_y_xxxyy_0, tr_y_xxxyyy_0, tr_y_xxxyyz_0, tr_y_xxxyzz_0, tr_y_xxxzz_0, tr_y_xxxzzz_0, tr_y_xxyy_0, tr_y_xxyyy_0, tr_y_xxyyyy_0, tr_y_xxyyyz_0, tr_y_xxyyzz_0, tr_y_xxyzz_0, tr_y_xxyzzz_0, tr_y_xxzz_0, tr_y_xxzzz_0, tr_y_xxzzzz_0, tr_y_xyyy_0, tr_y_xyyyy_0, tr_y_xyyyyy_0, tr_y_xyyyyz_0, tr_y_xyyyzz_0, tr_y_xyyzz_0, tr_y_xyyzzz_0, tr_y_xyzz_0, tr_y_xyzzz_0, tr_y_xyzzzz_0, tr_y_xzzz_0, tr_y_xzzzz_0, tr_y_xzzzzz_0, tr_y_yyyy_0, tr_y_yyyyy_0, tr_y_yyyyyy_0, tr_y_yyyyyz_0, tr_y_yyyyz_0, tr_y_yyyyzz_0, tr_y_yyyz_0, tr_y_yyyzz_0, tr_y_yyyzzz_0, tr_y_yyzz_0, tr_y_yyzzz_0, tr_y_yyzzzz_0, tr_y_yzzz_0, tr_y_yzzzz_0, tr_y_yzzzzz_0, tr_y_zzzz_0, tr_y_zzzzz_0, tr_y_zzzzzz_0, tr_z_xxxx_0, tr_z_xxxxx_0, tr_z_xxxxxx_0, tr_z_xxxxxy_0, tr_z_xxxxxz_0, tr_z_xxxxyy_0, tr_z_xxxxyz_0, tr_z_xxxxz_0, tr_z_xxxxzz_0, tr_z_xxxyy_0, tr_z_xxxyyy_0, tr_z_xxxyyz_0, tr_z_xxxyzz_0, tr_z_xxxz_0, tr_z_xxxzz_0, tr_z_xxxzzz_0, tr_z_xxyy_0, tr_z_xxyyy_0, tr_z_xxyyyy_0, tr_z_xxyyyz_0, tr_z_xxyyz_0, tr_z_xxyyzz_0, tr_z_xxyzzz_0, tr_z_xxzz_0, tr_z_xxzzz_0, tr_z_xxzzzz_0, tr_z_xyyy_0, tr_z_xyyyy_0, tr_z_xyyyyy_0, tr_z_xyyyyz_0, tr_z_xyyyz_0, tr_z_xyyyzz_0, tr_z_xyyz_0, tr_z_xyyzz_0, tr_z_xyyzzz_0, tr_z_xyzzzz_0, tr_z_xzzz_0, tr_z_xzzzz_0, tr_z_xzzzzz_0, tr_z_yyyy_0, tr_z_yyyyy_0, tr_z_yyyyyy_0, tr_z_yyyyyz_0, tr_z_yyyyz_0, tr_z_yyyyzz_0, tr_z_yyyz_0, tr_z_yyyzz_0, tr_z_yyyzzz_0, tr_z_yyzz_0, tr_z_yyzzz_0, tr_z_yyzzzz_0, tr_z_yzzz_0, tr_z_yzzzz_0, tr_z_yzzzzz_0, tr_z_zzzz_0, tr_z_zzzzz_0, tr_z_zzzzzz_0, ts_xxxxx_0, ts_yyyyy_0, ts_yyyzz_0, ts_yyzzz_0, ts_zzzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxxx_0[i] = 5.0 * tr_x_xxxx_0[i] * fe_0 + ts_xxxxx_0[i] * fe_0 + tr_x_xxxxx_0[i] * pa_x[i];

        tr_x_xxxxxy_0[i] = tr_x_xxxxx_0[i] * pa_y[i];

        tr_x_xxxxxz_0[i] = tr_x_xxxxx_0[i] * pa_z[i];

        tr_x_xxxxyy_0[i] = tr_x_xxxx_0[i] * fe_0 + tr_x_xxxxy_0[i] * pa_y[i];

        tr_x_xxxxyz_0[i] = tr_x_xxxxz_0[i] * pa_y[i];

        tr_x_xxxxzz_0[i] = tr_x_xxxx_0[i] * fe_0 + tr_x_xxxxz_0[i] * pa_z[i];

        tr_x_xxxyyy_0[i] = 2.0 * tr_x_xxxy_0[i] * fe_0 + tr_x_xxxyy_0[i] * pa_y[i];

        tr_x_xxxyyz_0[i] = tr_x_xxxyy_0[i] * pa_z[i];

        tr_x_xxxyzz_0[i] = tr_x_xxxzz_0[i] * pa_y[i];

        tr_x_xxxzzz_0[i] = 2.0 * tr_x_xxxz_0[i] * fe_0 + tr_x_xxxzz_0[i] * pa_z[i];

        tr_x_xxyyyy_0[i] = 3.0 * tr_x_xxyy_0[i] * fe_0 + tr_x_xxyyy_0[i] * pa_y[i];

        tr_x_xxyyyz_0[i] = tr_x_xxyyy_0[i] * pa_z[i];

        tr_x_xxyyzz_0[i] = tr_x_xxzz_0[i] * fe_0 + tr_x_xxyzz_0[i] * pa_y[i];

        tr_x_xxyzzz_0[i] = tr_x_xxzzz_0[i] * pa_y[i];

        tr_x_xxzzzz_0[i] = 3.0 * tr_x_xxzz_0[i] * fe_0 + tr_x_xxzzz_0[i] * pa_z[i];

        tr_x_xyyyyy_0[i] = ts_yyyyy_0[i] * fe_0 + tr_x_yyyyy_0[i] * pa_x[i];

        tr_x_xyyyyz_0[i] = tr_x_xyyyy_0[i] * pa_z[i];

        tr_x_xyyyzz_0[i] = ts_yyyzz_0[i] * fe_0 + tr_x_yyyzz_0[i] * pa_x[i];

        tr_x_xyyzzz_0[i] = ts_yyzzz_0[i] * fe_0 + tr_x_yyzzz_0[i] * pa_x[i];

        tr_x_xyzzzz_0[i] = tr_x_xzzzz_0[i] * pa_y[i];

        tr_x_xzzzzz_0[i] = ts_zzzzz_0[i] * fe_0 + tr_x_zzzzz_0[i] * pa_x[i];

        tr_x_yyyyyy_0[i] = 5.0 * tr_x_yyyy_0[i] * fe_0 + tr_x_yyyyy_0[i] * pa_y[i];

        tr_x_yyyyyz_0[i] = tr_x_yyyyy_0[i] * pa_z[i];

        tr_x_yyyyzz_0[i] = 3.0 * tr_x_yyzz_0[i] * fe_0 + tr_x_yyyzz_0[i] * pa_y[i];

        tr_x_yyyzzz_0[i] = 2.0 * tr_x_yzzz_0[i] * fe_0 + tr_x_yyzzz_0[i] * pa_y[i];

        tr_x_yyzzzz_0[i] = tr_x_zzzz_0[i] * fe_0 + tr_x_yzzzz_0[i] * pa_y[i];

        tr_x_yzzzzz_0[i] = tr_x_zzzzz_0[i] * pa_y[i];

        tr_x_zzzzzz_0[i] = 5.0 * tr_x_zzzz_0[i] * fe_0 + tr_x_zzzzz_0[i] * pa_z[i];

        tr_y_xxxxxx_0[i] = 5.0 * tr_y_xxxx_0[i] * fe_0 + tr_y_xxxxx_0[i] * pa_x[i];

        tr_y_xxxxxy_0[i] = 4.0 * tr_y_xxxy_0[i] * fe_0 + tr_y_xxxxy_0[i] * pa_x[i];

        tr_y_xxxxxz_0[i] = tr_y_xxxxx_0[i] * pa_z[i];

        tr_y_xxxxyy_0[i] = 3.0 * tr_y_xxyy_0[i] * fe_0 + tr_y_xxxyy_0[i] * pa_x[i];

        tr_y_xxxxyz_0[i] = tr_y_xxxxy_0[i] * pa_z[i];

        tr_y_xxxxzz_0[i] = 3.0 * tr_y_xxzz_0[i] * fe_0 + tr_y_xxxzz_0[i] * pa_x[i];

        tr_y_xxxyyy_0[i] = 2.0 * tr_y_xyyy_0[i] * fe_0 + tr_y_xxyyy_0[i] * pa_x[i];

        tr_y_xxxyyz_0[i] = tr_y_xxxyy_0[i] * pa_z[i];

        tr_y_xxxyzz_0[i] = 2.0 * tr_y_xyzz_0[i] * fe_0 + tr_y_xxyzz_0[i] * pa_x[i];

        tr_y_xxxzzz_0[i] = 2.0 * tr_y_xzzz_0[i] * fe_0 + tr_y_xxzzz_0[i] * pa_x[i];

        tr_y_xxyyyy_0[i] = tr_y_yyyy_0[i] * fe_0 + tr_y_xyyyy_0[i] * pa_x[i];

        tr_y_xxyyyz_0[i] = tr_y_xxyyy_0[i] * pa_z[i];

        tr_y_xxyyzz_0[i] = tr_y_yyzz_0[i] * fe_0 + tr_y_xyyzz_0[i] * pa_x[i];

        tr_y_xxyzzz_0[i] = tr_y_yzzz_0[i] * fe_0 + tr_y_xyzzz_0[i] * pa_x[i];

        tr_y_xxzzzz_0[i] = tr_y_zzzz_0[i] * fe_0 + tr_y_xzzzz_0[i] * pa_x[i];

        tr_y_xyyyyy_0[i] = tr_y_yyyyy_0[i] * pa_x[i];

        tr_y_xyyyyz_0[i] = tr_y_yyyyz_0[i] * pa_x[i];

        tr_y_xyyyzz_0[i] = tr_y_yyyzz_0[i] * pa_x[i];

        tr_y_xyyzzz_0[i] = tr_y_yyzzz_0[i] * pa_x[i];

        tr_y_xyzzzz_0[i] = tr_y_yzzzz_0[i] * pa_x[i];

        tr_y_xzzzzz_0[i] = tr_y_zzzzz_0[i] * pa_x[i];

        tr_y_yyyyyy_0[i] = 5.0 * tr_y_yyyy_0[i] * fe_0 + ts_yyyyy_0[i] * fe_0 + tr_y_yyyyy_0[i] * pa_y[i];

        tr_y_yyyyyz_0[i] = tr_y_yyyyy_0[i] * pa_z[i];

        tr_y_yyyyzz_0[i] = tr_y_yyyy_0[i] * fe_0 + tr_y_yyyyz_0[i] * pa_z[i];

        tr_y_yyyzzz_0[i] = 2.0 * tr_y_yyyz_0[i] * fe_0 + tr_y_yyyzz_0[i] * pa_z[i];

        tr_y_yyzzzz_0[i] = 3.0 * tr_y_yyzz_0[i] * fe_0 + tr_y_yyzzz_0[i] * pa_z[i];

        tr_y_yzzzzz_0[i] = ts_zzzzz_0[i] * fe_0 + tr_y_zzzzz_0[i] * pa_y[i];

        tr_y_zzzzzz_0[i] = 5.0 * tr_y_zzzz_0[i] * fe_0 + tr_y_zzzzz_0[i] * pa_z[i];

        tr_z_xxxxxx_0[i] = 5.0 * tr_z_xxxx_0[i] * fe_0 + tr_z_xxxxx_0[i] * pa_x[i];

        tr_z_xxxxxy_0[i] = tr_z_xxxxx_0[i] * pa_y[i];

        tr_z_xxxxxz_0[i] = 4.0 * tr_z_xxxz_0[i] * fe_0 + tr_z_xxxxz_0[i] * pa_x[i];

        tr_z_xxxxyy_0[i] = 3.0 * tr_z_xxyy_0[i] * fe_0 + tr_z_xxxyy_0[i] * pa_x[i];

        tr_z_xxxxyz_0[i] = tr_z_xxxxz_0[i] * pa_y[i];

        tr_z_xxxxzz_0[i] = 3.0 * tr_z_xxzz_0[i] * fe_0 + tr_z_xxxzz_0[i] * pa_x[i];

        tr_z_xxxyyy_0[i] = 2.0 * tr_z_xyyy_0[i] * fe_0 + tr_z_xxyyy_0[i] * pa_x[i];

        tr_z_xxxyyz_0[i] = 2.0 * tr_z_xyyz_0[i] * fe_0 + tr_z_xxyyz_0[i] * pa_x[i];

        tr_z_xxxyzz_0[i] = tr_z_xxxzz_0[i] * pa_y[i];

        tr_z_xxxzzz_0[i] = 2.0 * tr_z_xzzz_0[i] * fe_0 + tr_z_xxzzz_0[i] * pa_x[i];

        tr_z_xxyyyy_0[i] = tr_z_yyyy_0[i] * fe_0 + tr_z_xyyyy_0[i] * pa_x[i];

        tr_z_xxyyyz_0[i] = tr_z_yyyz_0[i] * fe_0 + tr_z_xyyyz_0[i] * pa_x[i];

        tr_z_xxyyzz_0[i] = tr_z_yyzz_0[i] * fe_0 + tr_z_xyyzz_0[i] * pa_x[i];

        tr_z_xxyzzz_0[i] = tr_z_xxzzz_0[i] * pa_y[i];

        tr_z_xxzzzz_0[i] = tr_z_zzzz_0[i] * fe_0 + tr_z_xzzzz_0[i] * pa_x[i];

        tr_z_xyyyyy_0[i] = tr_z_yyyyy_0[i] * pa_x[i];

        tr_z_xyyyyz_0[i] = tr_z_yyyyz_0[i] * pa_x[i];

        tr_z_xyyyzz_0[i] = tr_z_yyyzz_0[i] * pa_x[i];

        tr_z_xyyzzz_0[i] = tr_z_yyzzz_0[i] * pa_x[i];

        tr_z_xyzzzz_0[i] = tr_z_yzzzz_0[i] * pa_x[i];

        tr_z_xzzzzz_0[i] = tr_z_zzzzz_0[i] * pa_x[i];

        tr_z_yyyyyy_0[i] = 5.0 * tr_z_yyyy_0[i] * fe_0 + tr_z_yyyyy_0[i] * pa_y[i];

        tr_z_yyyyyz_0[i] = 4.0 * tr_z_yyyz_0[i] * fe_0 + tr_z_yyyyz_0[i] * pa_y[i];

        tr_z_yyyyzz_0[i] = 3.0 * tr_z_yyzz_0[i] * fe_0 + tr_z_yyyzz_0[i] * pa_y[i];

        tr_z_yyyzzz_0[i] = 2.0 * tr_z_yzzz_0[i] * fe_0 + tr_z_yyzzz_0[i] * pa_y[i];

        tr_z_yyzzzz_0[i] = tr_z_zzzz_0[i] * fe_0 + tr_z_yzzzz_0[i] * pa_y[i];

        tr_z_yzzzzz_0[i] = tr_z_zzzzz_0[i] * pa_y[i];

        tr_z_zzzzzz_0[i] = 5.0 * tr_z_zzzz_0[i] * fe_0 + ts_zzzzz_0[i] * fe_0 + tr_z_zzzzz_0[i] * pa_z[i];
    }
}

} // diprec namespace


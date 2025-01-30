#include "TwoCenterElectronRepulsionPrimRecDI.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_di(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_di,
                                const size_t idx_eri_0_si,
                                const size_t idx_eri_1_si,
                                const size_t idx_eri_1_ph,
                                const size_t idx_eri_1_pi,
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

    // Set up components of auxiliary buffer : SI

    auto g_0_xxxxxx_0 = pbuffer.data(idx_eri_0_si);

    auto g_0_xxxxxy_0 = pbuffer.data(idx_eri_0_si + 1);

    auto g_0_xxxxxz_0 = pbuffer.data(idx_eri_0_si + 2);

    auto g_0_xxxxyy_0 = pbuffer.data(idx_eri_0_si + 3);

    auto g_0_xxxxyz_0 = pbuffer.data(idx_eri_0_si + 4);

    auto g_0_xxxxzz_0 = pbuffer.data(idx_eri_0_si + 5);

    auto g_0_xxxyyy_0 = pbuffer.data(idx_eri_0_si + 6);

    auto g_0_xxxyyz_0 = pbuffer.data(idx_eri_0_si + 7);

    auto g_0_xxxyzz_0 = pbuffer.data(idx_eri_0_si + 8);

    auto g_0_xxxzzz_0 = pbuffer.data(idx_eri_0_si + 9);

    auto g_0_xxyyyy_0 = pbuffer.data(idx_eri_0_si + 10);

    auto g_0_xxyyyz_0 = pbuffer.data(idx_eri_0_si + 11);

    auto g_0_xxyyzz_0 = pbuffer.data(idx_eri_0_si + 12);

    auto g_0_xxyzzz_0 = pbuffer.data(idx_eri_0_si + 13);

    auto g_0_xxzzzz_0 = pbuffer.data(idx_eri_0_si + 14);

    auto g_0_xyyyyy_0 = pbuffer.data(idx_eri_0_si + 15);

    auto g_0_xyyyyz_0 = pbuffer.data(idx_eri_0_si + 16);

    auto g_0_xyyyzz_0 = pbuffer.data(idx_eri_0_si + 17);

    auto g_0_xyyzzz_0 = pbuffer.data(idx_eri_0_si + 18);

    auto g_0_xyzzzz_0 = pbuffer.data(idx_eri_0_si + 19);

    auto g_0_xzzzzz_0 = pbuffer.data(idx_eri_0_si + 20);

    auto g_0_yyyyyy_0 = pbuffer.data(idx_eri_0_si + 21);

    auto g_0_yyyyyz_0 = pbuffer.data(idx_eri_0_si + 22);

    auto g_0_yyyyzz_0 = pbuffer.data(idx_eri_0_si + 23);

    auto g_0_yyyzzz_0 = pbuffer.data(idx_eri_0_si + 24);

    auto g_0_yyzzzz_0 = pbuffer.data(idx_eri_0_si + 25);

    auto g_0_yzzzzz_0 = pbuffer.data(idx_eri_0_si + 26);

    auto g_0_zzzzzz_0 = pbuffer.data(idx_eri_0_si + 27);

    // Set up components of auxiliary buffer : SI

    auto g_0_xxxxxx_1 = pbuffer.data(idx_eri_1_si);

    auto g_0_xxxxxy_1 = pbuffer.data(idx_eri_1_si + 1);

    auto g_0_xxxxxz_1 = pbuffer.data(idx_eri_1_si + 2);

    auto g_0_xxxxyy_1 = pbuffer.data(idx_eri_1_si + 3);

    auto g_0_xxxxyz_1 = pbuffer.data(idx_eri_1_si + 4);

    auto g_0_xxxxzz_1 = pbuffer.data(idx_eri_1_si + 5);

    auto g_0_xxxyyy_1 = pbuffer.data(idx_eri_1_si + 6);

    auto g_0_xxxyyz_1 = pbuffer.data(idx_eri_1_si + 7);

    auto g_0_xxxyzz_1 = pbuffer.data(idx_eri_1_si + 8);

    auto g_0_xxxzzz_1 = pbuffer.data(idx_eri_1_si + 9);

    auto g_0_xxyyyy_1 = pbuffer.data(idx_eri_1_si + 10);

    auto g_0_xxyyyz_1 = pbuffer.data(idx_eri_1_si + 11);

    auto g_0_xxyyzz_1 = pbuffer.data(idx_eri_1_si + 12);

    auto g_0_xxyzzz_1 = pbuffer.data(idx_eri_1_si + 13);

    auto g_0_xxzzzz_1 = pbuffer.data(idx_eri_1_si + 14);

    auto g_0_xyyyyy_1 = pbuffer.data(idx_eri_1_si + 15);

    auto g_0_xyyyyz_1 = pbuffer.data(idx_eri_1_si + 16);

    auto g_0_xyyyzz_1 = pbuffer.data(idx_eri_1_si + 17);

    auto g_0_xyyzzz_1 = pbuffer.data(idx_eri_1_si + 18);

    auto g_0_xyzzzz_1 = pbuffer.data(idx_eri_1_si + 19);

    auto g_0_xzzzzz_1 = pbuffer.data(idx_eri_1_si + 20);

    auto g_0_yyyyyy_1 = pbuffer.data(idx_eri_1_si + 21);

    auto g_0_yyyyyz_1 = pbuffer.data(idx_eri_1_si + 22);

    auto g_0_yyyyzz_1 = pbuffer.data(idx_eri_1_si + 23);

    auto g_0_yyyzzz_1 = pbuffer.data(idx_eri_1_si + 24);

    auto g_0_yyzzzz_1 = pbuffer.data(idx_eri_1_si + 25);

    auto g_0_yzzzzz_1 = pbuffer.data(idx_eri_1_si + 26);

    auto g_0_zzzzzz_1 = pbuffer.data(idx_eri_1_si + 27);

    // Set up components of auxiliary buffer : PH

    auto g_x_xxxxx_1 = pbuffer.data(idx_eri_1_ph);

    auto g_x_xxxxy_1 = pbuffer.data(idx_eri_1_ph + 1);

    auto g_x_xxxxz_1 = pbuffer.data(idx_eri_1_ph + 2);

    auto g_x_xxxyy_1 = pbuffer.data(idx_eri_1_ph + 3);

    auto g_x_xxxyz_1 = pbuffer.data(idx_eri_1_ph + 4);

    auto g_x_xxxzz_1 = pbuffer.data(idx_eri_1_ph + 5);

    auto g_x_xxyyy_1 = pbuffer.data(idx_eri_1_ph + 6);

    auto g_x_xxyyz_1 = pbuffer.data(idx_eri_1_ph + 7);

    auto g_x_xxyzz_1 = pbuffer.data(idx_eri_1_ph + 8);

    auto g_x_xxzzz_1 = pbuffer.data(idx_eri_1_ph + 9);

    auto g_x_xyyyy_1 = pbuffer.data(idx_eri_1_ph + 10);

    auto g_x_xyyyz_1 = pbuffer.data(idx_eri_1_ph + 11);

    auto g_x_xyyzz_1 = pbuffer.data(idx_eri_1_ph + 12);

    auto g_x_xyzzz_1 = pbuffer.data(idx_eri_1_ph + 13);

    auto g_x_xzzzz_1 = pbuffer.data(idx_eri_1_ph + 14);

    auto g_x_yyyyy_1 = pbuffer.data(idx_eri_1_ph + 15);

    auto g_x_yyyyz_1 = pbuffer.data(idx_eri_1_ph + 16);

    auto g_x_yyyzz_1 = pbuffer.data(idx_eri_1_ph + 17);

    auto g_x_yyzzz_1 = pbuffer.data(idx_eri_1_ph + 18);

    auto g_x_yzzzz_1 = pbuffer.data(idx_eri_1_ph + 19);

    auto g_x_zzzzz_1 = pbuffer.data(idx_eri_1_ph + 20);

    auto g_y_xxxxx_1 = pbuffer.data(idx_eri_1_ph + 21);

    auto g_y_xxxxy_1 = pbuffer.data(idx_eri_1_ph + 22);

    auto g_y_xxxxz_1 = pbuffer.data(idx_eri_1_ph + 23);

    auto g_y_xxxyy_1 = pbuffer.data(idx_eri_1_ph + 24);

    auto g_y_xxxyz_1 = pbuffer.data(idx_eri_1_ph + 25);

    auto g_y_xxxzz_1 = pbuffer.data(idx_eri_1_ph + 26);

    auto g_y_xxyyy_1 = pbuffer.data(idx_eri_1_ph + 27);

    auto g_y_xxyyz_1 = pbuffer.data(idx_eri_1_ph + 28);

    auto g_y_xxyzz_1 = pbuffer.data(idx_eri_1_ph + 29);

    auto g_y_xxzzz_1 = pbuffer.data(idx_eri_1_ph + 30);

    auto g_y_xyyyy_1 = pbuffer.data(idx_eri_1_ph + 31);

    auto g_y_xyyyz_1 = pbuffer.data(idx_eri_1_ph + 32);

    auto g_y_xyyzz_1 = pbuffer.data(idx_eri_1_ph + 33);

    auto g_y_xyzzz_1 = pbuffer.data(idx_eri_1_ph + 34);

    auto g_y_xzzzz_1 = pbuffer.data(idx_eri_1_ph + 35);

    auto g_y_yyyyy_1 = pbuffer.data(idx_eri_1_ph + 36);

    auto g_y_yyyyz_1 = pbuffer.data(idx_eri_1_ph + 37);

    auto g_y_yyyzz_1 = pbuffer.data(idx_eri_1_ph + 38);

    auto g_y_yyzzz_1 = pbuffer.data(idx_eri_1_ph + 39);

    auto g_y_yzzzz_1 = pbuffer.data(idx_eri_1_ph + 40);

    auto g_y_zzzzz_1 = pbuffer.data(idx_eri_1_ph + 41);

    auto g_z_xxxxx_1 = pbuffer.data(idx_eri_1_ph + 42);

    auto g_z_xxxxy_1 = pbuffer.data(idx_eri_1_ph + 43);

    auto g_z_xxxxz_1 = pbuffer.data(idx_eri_1_ph + 44);

    auto g_z_xxxyy_1 = pbuffer.data(idx_eri_1_ph + 45);

    auto g_z_xxxyz_1 = pbuffer.data(idx_eri_1_ph + 46);

    auto g_z_xxxzz_1 = pbuffer.data(idx_eri_1_ph + 47);

    auto g_z_xxyyy_1 = pbuffer.data(idx_eri_1_ph + 48);

    auto g_z_xxyyz_1 = pbuffer.data(idx_eri_1_ph + 49);

    auto g_z_xxyzz_1 = pbuffer.data(idx_eri_1_ph + 50);

    auto g_z_xxzzz_1 = pbuffer.data(idx_eri_1_ph + 51);

    auto g_z_xyyyy_1 = pbuffer.data(idx_eri_1_ph + 52);

    auto g_z_xyyyz_1 = pbuffer.data(idx_eri_1_ph + 53);

    auto g_z_xyyzz_1 = pbuffer.data(idx_eri_1_ph + 54);

    auto g_z_xyzzz_1 = pbuffer.data(idx_eri_1_ph + 55);

    auto g_z_xzzzz_1 = pbuffer.data(idx_eri_1_ph + 56);

    auto g_z_yyyyy_1 = pbuffer.data(idx_eri_1_ph + 57);

    auto g_z_yyyyz_1 = pbuffer.data(idx_eri_1_ph + 58);

    auto g_z_yyyzz_1 = pbuffer.data(idx_eri_1_ph + 59);

    auto g_z_yyzzz_1 = pbuffer.data(idx_eri_1_ph + 60);

    auto g_z_yzzzz_1 = pbuffer.data(idx_eri_1_ph + 61);

    auto g_z_zzzzz_1 = pbuffer.data(idx_eri_1_ph + 62);

    // Set up components of auxiliary buffer : PI

    auto g_x_xxxxxx_1 = pbuffer.data(idx_eri_1_pi);

    auto g_x_xxxxxy_1 = pbuffer.data(idx_eri_1_pi + 1);

    auto g_x_xxxxxz_1 = pbuffer.data(idx_eri_1_pi + 2);

    auto g_x_xxxxyy_1 = pbuffer.data(idx_eri_1_pi + 3);

    auto g_x_xxxxyz_1 = pbuffer.data(idx_eri_1_pi + 4);

    auto g_x_xxxxzz_1 = pbuffer.data(idx_eri_1_pi + 5);

    auto g_x_xxxyyy_1 = pbuffer.data(idx_eri_1_pi + 6);

    auto g_x_xxxyyz_1 = pbuffer.data(idx_eri_1_pi + 7);

    auto g_x_xxxyzz_1 = pbuffer.data(idx_eri_1_pi + 8);

    auto g_x_xxxzzz_1 = pbuffer.data(idx_eri_1_pi + 9);

    auto g_x_xxyyyy_1 = pbuffer.data(idx_eri_1_pi + 10);

    auto g_x_xxyyyz_1 = pbuffer.data(idx_eri_1_pi + 11);

    auto g_x_xxyyzz_1 = pbuffer.data(idx_eri_1_pi + 12);

    auto g_x_xxyzzz_1 = pbuffer.data(idx_eri_1_pi + 13);

    auto g_x_xxzzzz_1 = pbuffer.data(idx_eri_1_pi + 14);

    auto g_x_xyyyyy_1 = pbuffer.data(idx_eri_1_pi + 15);

    auto g_x_xyyyyz_1 = pbuffer.data(idx_eri_1_pi + 16);

    auto g_x_xyyyzz_1 = pbuffer.data(idx_eri_1_pi + 17);

    auto g_x_xyyzzz_1 = pbuffer.data(idx_eri_1_pi + 18);

    auto g_x_xyzzzz_1 = pbuffer.data(idx_eri_1_pi + 19);

    auto g_x_xzzzzz_1 = pbuffer.data(idx_eri_1_pi + 20);

    auto g_x_yyyyyy_1 = pbuffer.data(idx_eri_1_pi + 21);

    auto g_x_yyyyyz_1 = pbuffer.data(idx_eri_1_pi + 22);

    auto g_x_yyyyzz_1 = pbuffer.data(idx_eri_1_pi + 23);

    auto g_x_yyyzzz_1 = pbuffer.data(idx_eri_1_pi + 24);

    auto g_x_yyzzzz_1 = pbuffer.data(idx_eri_1_pi + 25);

    auto g_x_yzzzzz_1 = pbuffer.data(idx_eri_1_pi + 26);

    auto g_x_zzzzzz_1 = pbuffer.data(idx_eri_1_pi + 27);

    auto g_y_xxxxxx_1 = pbuffer.data(idx_eri_1_pi + 28);

    auto g_y_xxxxxy_1 = pbuffer.data(idx_eri_1_pi + 29);

    auto g_y_xxxxxz_1 = pbuffer.data(idx_eri_1_pi + 30);

    auto g_y_xxxxyy_1 = pbuffer.data(idx_eri_1_pi + 31);

    auto g_y_xxxxyz_1 = pbuffer.data(idx_eri_1_pi + 32);

    auto g_y_xxxxzz_1 = pbuffer.data(idx_eri_1_pi + 33);

    auto g_y_xxxyyy_1 = pbuffer.data(idx_eri_1_pi + 34);

    auto g_y_xxxyyz_1 = pbuffer.data(idx_eri_1_pi + 35);

    auto g_y_xxxyzz_1 = pbuffer.data(idx_eri_1_pi + 36);

    auto g_y_xxxzzz_1 = pbuffer.data(idx_eri_1_pi + 37);

    auto g_y_xxyyyy_1 = pbuffer.data(idx_eri_1_pi + 38);

    auto g_y_xxyyyz_1 = pbuffer.data(idx_eri_1_pi + 39);

    auto g_y_xxyyzz_1 = pbuffer.data(idx_eri_1_pi + 40);

    auto g_y_xxyzzz_1 = pbuffer.data(idx_eri_1_pi + 41);

    auto g_y_xxzzzz_1 = pbuffer.data(idx_eri_1_pi + 42);

    auto g_y_xyyyyy_1 = pbuffer.data(idx_eri_1_pi + 43);

    auto g_y_xyyyyz_1 = pbuffer.data(idx_eri_1_pi + 44);

    auto g_y_xyyyzz_1 = pbuffer.data(idx_eri_1_pi + 45);

    auto g_y_xyyzzz_1 = pbuffer.data(idx_eri_1_pi + 46);

    auto g_y_xyzzzz_1 = pbuffer.data(idx_eri_1_pi + 47);

    auto g_y_xzzzzz_1 = pbuffer.data(idx_eri_1_pi + 48);

    auto g_y_yyyyyy_1 = pbuffer.data(idx_eri_1_pi + 49);

    auto g_y_yyyyyz_1 = pbuffer.data(idx_eri_1_pi + 50);

    auto g_y_yyyyzz_1 = pbuffer.data(idx_eri_1_pi + 51);

    auto g_y_yyyzzz_1 = pbuffer.data(idx_eri_1_pi + 52);

    auto g_y_yyzzzz_1 = pbuffer.data(idx_eri_1_pi + 53);

    auto g_y_yzzzzz_1 = pbuffer.data(idx_eri_1_pi + 54);

    auto g_y_zzzzzz_1 = pbuffer.data(idx_eri_1_pi + 55);

    auto g_z_xxxxxx_1 = pbuffer.data(idx_eri_1_pi + 56);

    auto g_z_xxxxxy_1 = pbuffer.data(idx_eri_1_pi + 57);

    auto g_z_xxxxxz_1 = pbuffer.data(idx_eri_1_pi + 58);

    auto g_z_xxxxyy_1 = pbuffer.data(idx_eri_1_pi + 59);

    auto g_z_xxxxyz_1 = pbuffer.data(idx_eri_1_pi + 60);

    auto g_z_xxxxzz_1 = pbuffer.data(idx_eri_1_pi + 61);

    auto g_z_xxxyyy_1 = pbuffer.data(idx_eri_1_pi + 62);

    auto g_z_xxxyyz_1 = pbuffer.data(idx_eri_1_pi + 63);

    auto g_z_xxxyzz_1 = pbuffer.data(idx_eri_1_pi + 64);

    auto g_z_xxxzzz_1 = pbuffer.data(idx_eri_1_pi + 65);

    auto g_z_xxyyyy_1 = pbuffer.data(idx_eri_1_pi + 66);

    auto g_z_xxyyyz_1 = pbuffer.data(idx_eri_1_pi + 67);

    auto g_z_xxyyzz_1 = pbuffer.data(idx_eri_1_pi + 68);

    auto g_z_xxyzzz_1 = pbuffer.data(idx_eri_1_pi + 69);

    auto g_z_xxzzzz_1 = pbuffer.data(idx_eri_1_pi + 70);

    auto g_z_xyyyyy_1 = pbuffer.data(idx_eri_1_pi + 71);

    auto g_z_xyyyyz_1 = pbuffer.data(idx_eri_1_pi + 72);

    auto g_z_xyyyzz_1 = pbuffer.data(idx_eri_1_pi + 73);

    auto g_z_xyyzzz_1 = pbuffer.data(idx_eri_1_pi + 74);

    auto g_z_xyzzzz_1 = pbuffer.data(idx_eri_1_pi + 75);

    auto g_z_xzzzzz_1 = pbuffer.data(idx_eri_1_pi + 76);

    auto g_z_yyyyyy_1 = pbuffer.data(idx_eri_1_pi + 77);

    auto g_z_yyyyyz_1 = pbuffer.data(idx_eri_1_pi + 78);

    auto g_z_yyyyzz_1 = pbuffer.data(idx_eri_1_pi + 79);

    auto g_z_yyyzzz_1 = pbuffer.data(idx_eri_1_pi + 80);

    auto g_z_yyzzzz_1 = pbuffer.data(idx_eri_1_pi + 81);

    auto g_z_yzzzzz_1 = pbuffer.data(idx_eri_1_pi + 82);

    auto g_z_zzzzzz_1 = pbuffer.data(idx_eri_1_pi + 83);

    // Set up 0-28 components of targeted buffer : DI

    auto g_xx_xxxxxx_0 = pbuffer.data(idx_eri_0_di);

    auto g_xx_xxxxxy_0 = pbuffer.data(idx_eri_0_di + 1);

    auto g_xx_xxxxxz_0 = pbuffer.data(idx_eri_0_di + 2);

    auto g_xx_xxxxyy_0 = pbuffer.data(idx_eri_0_di + 3);

    auto g_xx_xxxxyz_0 = pbuffer.data(idx_eri_0_di + 4);

    auto g_xx_xxxxzz_0 = pbuffer.data(idx_eri_0_di + 5);

    auto g_xx_xxxyyy_0 = pbuffer.data(idx_eri_0_di + 6);

    auto g_xx_xxxyyz_0 = pbuffer.data(idx_eri_0_di + 7);

    auto g_xx_xxxyzz_0 = pbuffer.data(idx_eri_0_di + 8);

    auto g_xx_xxxzzz_0 = pbuffer.data(idx_eri_0_di + 9);

    auto g_xx_xxyyyy_0 = pbuffer.data(idx_eri_0_di + 10);

    auto g_xx_xxyyyz_0 = pbuffer.data(idx_eri_0_di + 11);

    auto g_xx_xxyyzz_0 = pbuffer.data(idx_eri_0_di + 12);

    auto g_xx_xxyzzz_0 = pbuffer.data(idx_eri_0_di + 13);

    auto g_xx_xxzzzz_0 = pbuffer.data(idx_eri_0_di + 14);

    auto g_xx_xyyyyy_0 = pbuffer.data(idx_eri_0_di + 15);

    auto g_xx_xyyyyz_0 = pbuffer.data(idx_eri_0_di + 16);

    auto g_xx_xyyyzz_0 = pbuffer.data(idx_eri_0_di + 17);

    auto g_xx_xyyzzz_0 = pbuffer.data(idx_eri_0_di + 18);

    auto g_xx_xyzzzz_0 = pbuffer.data(idx_eri_0_di + 19);

    auto g_xx_xzzzzz_0 = pbuffer.data(idx_eri_0_di + 20);

    auto g_xx_yyyyyy_0 = pbuffer.data(idx_eri_0_di + 21);

    auto g_xx_yyyyyz_0 = pbuffer.data(idx_eri_0_di + 22);

    auto g_xx_yyyyzz_0 = pbuffer.data(idx_eri_0_di + 23);

    auto g_xx_yyyzzz_0 = pbuffer.data(idx_eri_0_di + 24);

    auto g_xx_yyzzzz_0 = pbuffer.data(idx_eri_0_di + 25);

    auto g_xx_yzzzzz_0 = pbuffer.data(idx_eri_0_di + 26);

    auto g_xx_zzzzzz_0 = pbuffer.data(idx_eri_0_di + 27);

    #pragma omp simd aligned(g_0_xxxxxx_0, g_0_xxxxxx_1, g_0_xxxxxy_0, g_0_xxxxxy_1, g_0_xxxxxz_0, g_0_xxxxxz_1, g_0_xxxxyy_0, g_0_xxxxyy_1, g_0_xxxxyz_0, g_0_xxxxyz_1, g_0_xxxxzz_0, g_0_xxxxzz_1, g_0_xxxyyy_0, g_0_xxxyyy_1, g_0_xxxyyz_0, g_0_xxxyyz_1, g_0_xxxyzz_0, g_0_xxxyzz_1, g_0_xxxzzz_0, g_0_xxxzzz_1, g_0_xxyyyy_0, g_0_xxyyyy_1, g_0_xxyyyz_0, g_0_xxyyyz_1, g_0_xxyyzz_0, g_0_xxyyzz_1, g_0_xxyzzz_0, g_0_xxyzzz_1, g_0_xxzzzz_0, g_0_xxzzzz_1, g_0_xyyyyy_0, g_0_xyyyyy_1, g_0_xyyyyz_0, g_0_xyyyyz_1, g_0_xyyyzz_0, g_0_xyyyzz_1, g_0_xyyzzz_0, g_0_xyyzzz_1, g_0_xyzzzz_0, g_0_xyzzzz_1, g_0_xzzzzz_0, g_0_xzzzzz_1, g_0_yyyyyy_0, g_0_yyyyyy_1, g_0_yyyyyz_0, g_0_yyyyyz_1, g_0_yyyyzz_0, g_0_yyyyzz_1, g_0_yyyzzz_0, g_0_yyyzzz_1, g_0_yyzzzz_0, g_0_yyzzzz_1, g_0_yzzzzz_0, g_0_yzzzzz_1, g_0_zzzzzz_0, g_0_zzzzzz_1, g_x_xxxxx_1, g_x_xxxxxx_1, g_x_xxxxxy_1, g_x_xxxxxz_1, g_x_xxxxy_1, g_x_xxxxyy_1, g_x_xxxxyz_1, g_x_xxxxz_1, g_x_xxxxzz_1, g_x_xxxyy_1, g_x_xxxyyy_1, g_x_xxxyyz_1, g_x_xxxyz_1, g_x_xxxyzz_1, g_x_xxxzz_1, g_x_xxxzzz_1, g_x_xxyyy_1, g_x_xxyyyy_1, g_x_xxyyyz_1, g_x_xxyyz_1, g_x_xxyyzz_1, g_x_xxyzz_1, g_x_xxyzzz_1, g_x_xxzzz_1, g_x_xxzzzz_1, g_x_xyyyy_1, g_x_xyyyyy_1, g_x_xyyyyz_1, g_x_xyyyz_1, g_x_xyyyzz_1, g_x_xyyzz_1, g_x_xyyzzz_1, g_x_xyzzz_1, g_x_xyzzzz_1, g_x_xzzzz_1, g_x_xzzzzz_1, g_x_yyyyy_1, g_x_yyyyyy_1, g_x_yyyyyz_1, g_x_yyyyz_1, g_x_yyyyzz_1, g_x_yyyzz_1, g_x_yyyzzz_1, g_x_yyzzz_1, g_x_yyzzzz_1, g_x_yzzzz_1, g_x_yzzzzz_1, g_x_zzzzz_1, g_x_zzzzzz_1, g_xx_xxxxxx_0, g_xx_xxxxxy_0, g_xx_xxxxxz_0, g_xx_xxxxyy_0, g_xx_xxxxyz_0, g_xx_xxxxzz_0, g_xx_xxxyyy_0, g_xx_xxxyyz_0, g_xx_xxxyzz_0, g_xx_xxxzzz_0, g_xx_xxyyyy_0, g_xx_xxyyyz_0, g_xx_xxyyzz_0, g_xx_xxyzzz_0, g_xx_xxzzzz_0, g_xx_xyyyyy_0, g_xx_xyyyyz_0, g_xx_xyyyzz_0, g_xx_xyyzzz_0, g_xx_xyzzzz_0, g_xx_xzzzzz_0, g_xx_yyyyyy_0, g_xx_yyyyyz_0, g_xx_yyyyzz_0, g_xx_yyyzzz_0, g_xx_yyzzzz_0, g_xx_yzzzzz_0, g_xx_zzzzzz_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xx_xxxxxx_0[i] = g_0_xxxxxx_0[i] * fbe_0 - g_0_xxxxxx_1[i] * fz_be_0 + 6.0 * g_x_xxxxx_1[i] * fe_0 + g_x_xxxxxx_1[i] * pa_x[i];

        g_xx_xxxxxy_0[i] = g_0_xxxxxy_0[i] * fbe_0 - g_0_xxxxxy_1[i] * fz_be_0 + 5.0 * g_x_xxxxy_1[i] * fe_0 + g_x_xxxxxy_1[i] * pa_x[i];

        g_xx_xxxxxz_0[i] = g_0_xxxxxz_0[i] * fbe_0 - g_0_xxxxxz_1[i] * fz_be_0 + 5.0 * g_x_xxxxz_1[i] * fe_0 + g_x_xxxxxz_1[i] * pa_x[i];

        g_xx_xxxxyy_0[i] = g_0_xxxxyy_0[i] * fbe_0 - g_0_xxxxyy_1[i] * fz_be_0 + 4.0 * g_x_xxxyy_1[i] * fe_0 + g_x_xxxxyy_1[i] * pa_x[i];

        g_xx_xxxxyz_0[i] = g_0_xxxxyz_0[i] * fbe_0 - g_0_xxxxyz_1[i] * fz_be_0 + 4.0 * g_x_xxxyz_1[i] * fe_0 + g_x_xxxxyz_1[i] * pa_x[i];

        g_xx_xxxxzz_0[i] = g_0_xxxxzz_0[i] * fbe_0 - g_0_xxxxzz_1[i] * fz_be_0 + 4.0 * g_x_xxxzz_1[i] * fe_0 + g_x_xxxxzz_1[i] * pa_x[i];

        g_xx_xxxyyy_0[i] = g_0_xxxyyy_0[i] * fbe_0 - g_0_xxxyyy_1[i] * fz_be_0 + 3.0 * g_x_xxyyy_1[i] * fe_0 + g_x_xxxyyy_1[i] * pa_x[i];

        g_xx_xxxyyz_0[i] = g_0_xxxyyz_0[i] * fbe_0 - g_0_xxxyyz_1[i] * fz_be_0 + 3.0 * g_x_xxyyz_1[i] * fe_0 + g_x_xxxyyz_1[i] * pa_x[i];

        g_xx_xxxyzz_0[i] = g_0_xxxyzz_0[i] * fbe_0 - g_0_xxxyzz_1[i] * fz_be_0 + 3.0 * g_x_xxyzz_1[i] * fe_0 + g_x_xxxyzz_1[i] * pa_x[i];

        g_xx_xxxzzz_0[i] = g_0_xxxzzz_0[i] * fbe_0 - g_0_xxxzzz_1[i] * fz_be_0 + 3.0 * g_x_xxzzz_1[i] * fe_0 + g_x_xxxzzz_1[i] * pa_x[i];

        g_xx_xxyyyy_0[i] = g_0_xxyyyy_0[i] * fbe_0 - g_0_xxyyyy_1[i] * fz_be_0 + 2.0 * g_x_xyyyy_1[i] * fe_0 + g_x_xxyyyy_1[i] * pa_x[i];

        g_xx_xxyyyz_0[i] = g_0_xxyyyz_0[i] * fbe_0 - g_0_xxyyyz_1[i] * fz_be_0 + 2.0 * g_x_xyyyz_1[i] * fe_0 + g_x_xxyyyz_1[i] * pa_x[i];

        g_xx_xxyyzz_0[i] = g_0_xxyyzz_0[i] * fbe_0 - g_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_x_xyyzz_1[i] * fe_0 + g_x_xxyyzz_1[i] * pa_x[i];

        g_xx_xxyzzz_0[i] = g_0_xxyzzz_0[i] * fbe_0 - g_0_xxyzzz_1[i] * fz_be_0 + 2.0 * g_x_xyzzz_1[i] * fe_0 + g_x_xxyzzz_1[i] * pa_x[i];

        g_xx_xxzzzz_0[i] = g_0_xxzzzz_0[i] * fbe_0 - g_0_xxzzzz_1[i] * fz_be_0 + 2.0 * g_x_xzzzz_1[i] * fe_0 + g_x_xxzzzz_1[i] * pa_x[i];

        g_xx_xyyyyy_0[i] = g_0_xyyyyy_0[i] * fbe_0 - g_0_xyyyyy_1[i] * fz_be_0 + g_x_yyyyy_1[i] * fe_0 + g_x_xyyyyy_1[i] * pa_x[i];

        g_xx_xyyyyz_0[i] = g_0_xyyyyz_0[i] * fbe_0 - g_0_xyyyyz_1[i] * fz_be_0 + g_x_yyyyz_1[i] * fe_0 + g_x_xyyyyz_1[i] * pa_x[i];

        g_xx_xyyyzz_0[i] = g_0_xyyyzz_0[i] * fbe_0 - g_0_xyyyzz_1[i] * fz_be_0 + g_x_yyyzz_1[i] * fe_0 + g_x_xyyyzz_1[i] * pa_x[i];

        g_xx_xyyzzz_0[i] = g_0_xyyzzz_0[i] * fbe_0 - g_0_xyyzzz_1[i] * fz_be_0 + g_x_yyzzz_1[i] * fe_0 + g_x_xyyzzz_1[i] * pa_x[i];

        g_xx_xyzzzz_0[i] = g_0_xyzzzz_0[i] * fbe_0 - g_0_xyzzzz_1[i] * fz_be_0 + g_x_yzzzz_1[i] * fe_0 + g_x_xyzzzz_1[i] * pa_x[i];

        g_xx_xzzzzz_0[i] = g_0_xzzzzz_0[i] * fbe_0 - g_0_xzzzzz_1[i] * fz_be_0 + g_x_zzzzz_1[i] * fe_0 + g_x_xzzzzz_1[i] * pa_x[i];

        g_xx_yyyyyy_0[i] = g_0_yyyyyy_0[i] * fbe_0 - g_0_yyyyyy_1[i] * fz_be_0 + g_x_yyyyyy_1[i] * pa_x[i];

        g_xx_yyyyyz_0[i] = g_0_yyyyyz_0[i] * fbe_0 - g_0_yyyyyz_1[i] * fz_be_0 + g_x_yyyyyz_1[i] * pa_x[i];

        g_xx_yyyyzz_0[i] = g_0_yyyyzz_0[i] * fbe_0 - g_0_yyyyzz_1[i] * fz_be_0 + g_x_yyyyzz_1[i] * pa_x[i];

        g_xx_yyyzzz_0[i] = g_0_yyyzzz_0[i] * fbe_0 - g_0_yyyzzz_1[i] * fz_be_0 + g_x_yyyzzz_1[i] * pa_x[i];

        g_xx_yyzzzz_0[i] = g_0_yyzzzz_0[i] * fbe_0 - g_0_yyzzzz_1[i] * fz_be_0 + g_x_yyzzzz_1[i] * pa_x[i];

        g_xx_yzzzzz_0[i] = g_0_yzzzzz_0[i] * fbe_0 - g_0_yzzzzz_1[i] * fz_be_0 + g_x_yzzzzz_1[i] * pa_x[i];

        g_xx_zzzzzz_0[i] = g_0_zzzzzz_0[i] * fbe_0 - g_0_zzzzzz_1[i] * fz_be_0 + g_x_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 28-56 components of targeted buffer : DI

    auto g_xy_xxxxxx_0 = pbuffer.data(idx_eri_0_di + 28);

    auto g_xy_xxxxxy_0 = pbuffer.data(idx_eri_0_di + 29);

    auto g_xy_xxxxxz_0 = pbuffer.data(idx_eri_0_di + 30);

    auto g_xy_xxxxyy_0 = pbuffer.data(idx_eri_0_di + 31);

    auto g_xy_xxxxyz_0 = pbuffer.data(idx_eri_0_di + 32);

    auto g_xy_xxxxzz_0 = pbuffer.data(idx_eri_0_di + 33);

    auto g_xy_xxxyyy_0 = pbuffer.data(idx_eri_0_di + 34);

    auto g_xy_xxxyyz_0 = pbuffer.data(idx_eri_0_di + 35);

    auto g_xy_xxxyzz_0 = pbuffer.data(idx_eri_0_di + 36);

    auto g_xy_xxxzzz_0 = pbuffer.data(idx_eri_0_di + 37);

    auto g_xy_xxyyyy_0 = pbuffer.data(idx_eri_0_di + 38);

    auto g_xy_xxyyyz_0 = pbuffer.data(idx_eri_0_di + 39);

    auto g_xy_xxyyzz_0 = pbuffer.data(idx_eri_0_di + 40);

    auto g_xy_xxyzzz_0 = pbuffer.data(idx_eri_0_di + 41);

    auto g_xy_xxzzzz_0 = pbuffer.data(idx_eri_0_di + 42);

    auto g_xy_xyyyyy_0 = pbuffer.data(idx_eri_0_di + 43);

    auto g_xy_xyyyyz_0 = pbuffer.data(idx_eri_0_di + 44);

    auto g_xy_xyyyzz_0 = pbuffer.data(idx_eri_0_di + 45);

    auto g_xy_xyyzzz_0 = pbuffer.data(idx_eri_0_di + 46);

    auto g_xy_xyzzzz_0 = pbuffer.data(idx_eri_0_di + 47);

    auto g_xy_xzzzzz_0 = pbuffer.data(idx_eri_0_di + 48);

    auto g_xy_yyyyyy_0 = pbuffer.data(idx_eri_0_di + 49);

    auto g_xy_yyyyyz_0 = pbuffer.data(idx_eri_0_di + 50);

    auto g_xy_yyyyzz_0 = pbuffer.data(idx_eri_0_di + 51);

    auto g_xy_yyyzzz_0 = pbuffer.data(idx_eri_0_di + 52);

    auto g_xy_yyzzzz_0 = pbuffer.data(idx_eri_0_di + 53);

    auto g_xy_yzzzzz_0 = pbuffer.data(idx_eri_0_di + 54);

    auto g_xy_zzzzzz_0 = pbuffer.data(idx_eri_0_di + 55);

    #pragma omp simd aligned(g_x_xxxxxx_1, g_x_xxxxxz_1, g_x_xxxxzz_1, g_x_xxxzzz_1, g_x_xxzzzz_1, g_x_xzzzzz_1, g_xy_xxxxxx_0, g_xy_xxxxxy_0, g_xy_xxxxxz_0, g_xy_xxxxyy_0, g_xy_xxxxyz_0, g_xy_xxxxzz_0, g_xy_xxxyyy_0, g_xy_xxxyyz_0, g_xy_xxxyzz_0, g_xy_xxxzzz_0, g_xy_xxyyyy_0, g_xy_xxyyyz_0, g_xy_xxyyzz_0, g_xy_xxyzzz_0, g_xy_xxzzzz_0, g_xy_xyyyyy_0, g_xy_xyyyyz_0, g_xy_xyyyzz_0, g_xy_xyyzzz_0, g_xy_xyzzzz_0, g_xy_xzzzzz_0, g_xy_yyyyyy_0, g_xy_yyyyyz_0, g_xy_yyyyzz_0, g_xy_yyyzzz_0, g_xy_yyzzzz_0, g_xy_yzzzzz_0, g_xy_zzzzzz_0, g_y_xxxxxy_1, g_y_xxxxy_1, g_y_xxxxyy_1, g_y_xxxxyz_1, g_y_xxxyy_1, g_y_xxxyyy_1, g_y_xxxyyz_1, g_y_xxxyz_1, g_y_xxxyzz_1, g_y_xxyyy_1, g_y_xxyyyy_1, g_y_xxyyyz_1, g_y_xxyyz_1, g_y_xxyyzz_1, g_y_xxyzz_1, g_y_xxyzzz_1, g_y_xyyyy_1, g_y_xyyyyy_1, g_y_xyyyyz_1, g_y_xyyyz_1, g_y_xyyyzz_1, g_y_xyyzz_1, g_y_xyyzzz_1, g_y_xyzzz_1, g_y_xyzzzz_1, g_y_yyyyy_1, g_y_yyyyyy_1, g_y_yyyyyz_1, g_y_yyyyz_1, g_y_yyyyzz_1, g_y_yyyzz_1, g_y_yyyzzz_1, g_y_yyzzz_1, g_y_yyzzzz_1, g_y_yzzzz_1, g_y_yzzzzz_1, g_y_zzzzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xy_xxxxxx_0[i] = g_x_xxxxxx_1[i] * pa_y[i];

        g_xy_xxxxxy_0[i] = 5.0 * g_y_xxxxy_1[i] * fe_0 + g_y_xxxxxy_1[i] * pa_x[i];

        g_xy_xxxxxz_0[i] = g_x_xxxxxz_1[i] * pa_y[i];

        g_xy_xxxxyy_0[i] = 4.0 * g_y_xxxyy_1[i] * fe_0 + g_y_xxxxyy_1[i] * pa_x[i];

        g_xy_xxxxyz_0[i] = 4.0 * g_y_xxxyz_1[i] * fe_0 + g_y_xxxxyz_1[i] * pa_x[i];

        g_xy_xxxxzz_0[i] = g_x_xxxxzz_1[i] * pa_y[i];

        g_xy_xxxyyy_0[i] = 3.0 * g_y_xxyyy_1[i] * fe_0 + g_y_xxxyyy_1[i] * pa_x[i];

        g_xy_xxxyyz_0[i] = 3.0 * g_y_xxyyz_1[i] * fe_0 + g_y_xxxyyz_1[i] * pa_x[i];

        g_xy_xxxyzz_0[i] = 3.0 * g_y_xxyzz_1[i] * fe_0 + g_y_xxxyzz_1[i] * pa_x[i];

        g_xy_xxxzzz_0[i] = g_x_xxxzzz_1[i] * pa_y[i];

        g_xy_xxyyyy_0[i] = 2.0 * g_y_xyyyy_1[i] * fe_0 + g_y_xxyyyy_1[i] * pa_x[i];

        g_xy_xxyyyz_0[i] = 2.0 * g_y_xyyyz_1[i] * fe_0 + g_y_xxyyyz_1[i] * pa_x[i];

        g_xy_xxyyzz_0[i] = 2.0 * g_y_xyyzz_1[i] * fe_0 + g_y_xxyyzz_1[i] * pa_x[i];

        g_xy_xxyzzz_0[i] = 2.0 * g_y_xyzzz_1[i] * fe_0 + g_y_xxyzzz_1[i] * pa_x[i];

        g_xy_xxzzzz_0[i] = g_x_xxzzzz_1[i] * pa_y[i];

        g_xy_xyyyyy_0[i] = g_y_yyyyy_1[i] * fe_0 + g_y_xyyyyy_1[i] * pa_x[i];

        g_xy_xyyyyz_0[i] = g_y_yyyyz_1[i] * fe_0 + g_y_xyyyyz_1[i] * pa_x[i];

        g_xy_xyyyzz_0[i] = g_y_yyyzz_1[i] * fe_0 + g_y_xyyyzz_1[i] * pa_x[i];

        g_xy_xyyzzz_0[i] = g_y_yyzzz_1[i] * fe_0 + g_y_xyyzzz_1[i] * pa_x[i];

        g_xy_xyzzzz_0[i] = g_y_yzzzz_1[i] * fe_0 + g_y_xyzzzz_1[i] * pa_x[i];

        g_xy_xzzzzz_0[i] = g_x_xzzzzz_1[i] * pa_y[i];

        g_xy_yyyyyy_0[i] = g_y_yyyyyy_1[i] * pa_x[i];

        g_xy_yyyyyz_0[i] = g_y_yyyyyz_1[i] * pa_x[i];

        g_xy_yyyyzz_0[i] = g_y_yyyyzz_1[i] * pa_x[i];

        g_xy_yyyzzz_0[i] = g_y_yyyzzz_1[i] * pa_x[i];

        g_xy_yyzzzz_0[i] = g_y_yyzzzz_1[i] * pa_x[i];

        g_xy_yzzzzz_0[i] = g_y_yzzzzz_1[i] * pa_x[i];

        g_xy_zzzzzz_0[i] = g_y_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 56-84 components of targeted buffer : DI

    auto g_xz_xxxxxx_0 = pbuffer.data(idx_eri_0_di + 56);

    auto g_xz_xxxxxy_0 = pbuffer.data(idx_eri_0_di + 57);

    auto g_xz_xxxxxz_0 = pbuffer.data(idx_eri_0_di + 58);

    auto g_xz_xxxxyy_0 = pbuffer.data(idx_eri_0_di + 59);

    auto g_xz_xxxxyz_0 = pbuffer.data(idx_eri_0_di + 60);

    auto g_xz_xxxxzz_0 = pbuffer.data(idx_eri_0_di + 61);

    auto g_xz_xxxyyy_0 = pbuffer.data(idx_eri_0_di + 62);

    auto g_xz_xxxyyz_0 = pbuffer.data(idx_eri_0_di + 63);

    auto g_xz_xxxyzz_0 = pbuffer.data(idx_eri_0_di + 64);

    auto g_xz_xxxzzz_0 = pbuffer.data(idx_eri_0_di + 65);

    auto g_xz_xxyyyy_0 = pbuffer.data(idx_eri_0_di + 66);

    auto g_xz_xxyyyz_0 = pbuffer.data(idx_eri_0_di + 67);

    auto g_xz_xxyyzz_0 = pbuffer.data(idx_eri_0_di + 68);

    auto g_xz_xxyzzz_0 = pbuffer.data(idx_eri_0_di + 69);

    auto g_xz_xxzzzz_0 = pbuffer.data(idx_eri_0_di + 70);

    auto g_xz_xyyyyy_0 = pbuffer.data(idx_eri_0_di + 71);

    auto g_xz_xyyyyz_0 = pbuffer.data(idx_eri_0_di + 72);

    auto g_xz_xyyyzz_0 = pbuffer.data(idx_eri_0_di + 73);

    auto g_xz_xyyzzz_0 = pbuffer.data(idx_eri_0_di + 74);

    auto g_xz_xyzzzz_0 = pbuffer.data(idx_eri_0_di + 75);

    auto g_xz_xzzzzz_0 = pbuffer.data(idx_eri_0_di + 76);

    auto g_xz_yyyyyy_0 = pbuffer.data(idx_eri_0_di + 77);

    auto g_xz_yyyyyz_0 = pbuffer.data(idx_eri_0_di + 78);

    auto g_xz_yyyyzz_0 = pbuffer.data(idx_eri_0_di + 79);

    auto g_xz_yyyzzz_0 = pbuffer.data(idx_eri_0_di + 80);

    auto g_xz_yyzzzz_0 = pbuffer.data(idx_eri_0_di + 81);

    auto g_xz_yzzzzz_0 = pbuffer.data(idx_eri_0_di + 82);

    auto g_xz_zzzzzz_0 = pbuffer.data(idx_eri_0_di + 83);

    #pragma omp simd aligned(g_x_xxxxxx_1, g_x_xxxxxy_1, g_x_xxxxyy_1, g_x_xxxyyy_1, g_x_xxyyyy_1, g_x_xyyyyy_1, g_xz_xxxxxx_0, g_xz_xxxxxy_0, g_xz_xxxxxz_0, g_xz_xxxxyy_0, g_xz_xxxxyz_0, g_xz_xxxxzz_0, g_xz_xxxyyy_0, g_xz_xxxyyz_0, g_xz_xxxyzz_0, g_xz_xxxzzz_0, g_xz_xxyyyy_0, g_xz_xxyyyz_0, g_xz_xxyyzz_0, g_xz_xxyzzz_0, g_xz_xxzzzz_0, g_xz_xyyyyy_0, g_xz_xyyyyz_0, g_xz_xyyyzz_0, g_xz_xyyzzz_0, g_xz_xyzzzz_0, g_xz_xzzzzz_0, g_xz_yyyyyy_0, g_xz_yyyyyz_0, g_xz_yyyyzz_0, g_xz_yyyzzz_0, g_xz_yyzzzz_0, g_xz_yzzzzz_0, g_xz_zzzzzz_0, g_z_xxxxxz_1, g_z_xxxxyz_1, g_z_xxxxz_1, g_z_xxxxzz_1, g_z_xxxyyz_1, g_z_xxxyz_1, g_z_xxxyzz_1, g_z_xxxzz_1, g_z_xxxzzz_1, g_z_xxyyyz_1, g_z_xxyyz_1, g_z_xxyyzz_1, g_z_xxyzz_1, g_z_xxyzzz_1, g_z_xxzzz_1, g_z_xxzzzz_1, g_z_xyyyyz_1, g_z_xyyyz_1, g_z_xyyyzz_1, g_z_xyyzz_1, g_z_xyyzzz_1, g_z_xyzzz_1, g_z_xyzzzz_1, g_z_xzzzz_1, g_z_xzzzzz_1, g_z_yyyyyy_1, g_z_yyyyyz_1, g_z_yyyyz_1, g_z_yyyyzz_1, g_z_yyyzz_1, g_z_yyyzzz_1, g_z_yyzzz_1, g_z_yyzzzz_1, g_z_yzzzz_1, g_z_yzzzzz_1, g_z_zzzzz_1, g_z_zzzzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xz_xxxxxx_0[i] = g_x_xxxxxx_1[i] * pa_z[i];

        g_xz_xxxxxy_0[i] = g_x_xxxxxy_1[i] * pa_z[i];

        g_xz_xxxxxz_0[i] = 5.0 * g_z_xxxxz_1[i] * fe_0 + g_z_xxxxxz_1[i] * pa_x[i];

        g_xz_xxxxyy_0[i] = g_x_xxxxyy_1[i] * pa_z[i];

        g_xz_xxxxyz_0[i] = 4.0 * g_z_xxxyz_1[i] * fe_0 + g_z_xxxxyz_1[i] * pa_x[i];

        g_xz_xxxxzz_0[i] = 4.0 * g_z_xxxzz_1[i] * fe_0 + g_z_xxxxzz_1[i] * pa_x[i];

        g_xz_xxxyyy_0[i] = g_x_xxxyyy_1[i] * pa_z[i];

        g_xz_xxxyyz_0[i] = 3.0 * g_z_xxyyz_1[i] * fe_0 + g_z_xxxyyz_1[i] * pa_x[i];

        g_xz_xxxyzz_0[i] = 3.0 * g_z_xxyzz_1[i] * fe_0 + g_z_xxxyzz_1[i] * pa_x[i];

        g_xz_xxxzzz_0[i] = 3.0 * g_z_xxzzz_1[i] * fe_0 + g_z_xxxzzz_1[i] * pa_x[i];

        g_xz_xxyyyy_0[i] = g_x_xxyyyy_1[i] * pa_z[i];

        g_xz_xxyyyz_0[i] = 2.0 * g_z_xyyyz_1[i] * fe_0 + g_z_xxyyyz_1[i] * pa_x[i];

        g_xz_xxyyzz_0[i] = 2.0 * g_z_xyyzz_1[i] * fe_0 + g_z_xxyyzz_1[i] * pa_x[i];

        g_xz_xxyzzz_0[i] = 2.0 * g_z_xyzzz_1[i] * fe_0 + g_z_xxyzzz_1[i] * pa_x[i];

        g_xz_xxzzzz_0[i] = 2.0 * g_z_xzzzz_1[i] * fe_0 + g_z_xxzzzz_1[i] * pa_x[i];

        g_xz_xyyyyy_0[i] = g_x_xyyyyy_1[i] * pa_z[i];

        g_xz_xyyyyz_0[i] = g_z_yyyyz_1[i] * fe_0 + g_z_xyyyyz_1[i] * pa_x[i];

        g_xz_xyyyzz_0[i] = g_z_yyyzz_1[i] * fe_0 + g_z_xyyyzz_1[i] * pa_x[i];

        g_xz_xyyzzz_0[i] = g_z_yyzzz_1[i] * fe_0 + g_z_xyyzzz_1[i] * pa_x[i];

        g_xz_xyzzzz_0[i] = g_z_yzzzz_1[i] * fe_0 + g_z_xyzzzz_1[i] * pa_x[i];

        g_xz_xzzzzz_0[i] = g_z_zzzzz_1[i] * fe_0 + g_z_xzzzzz_1[i] * pa_x[i];

        g_xz_yyyyyy_0[i] = g_z_yyyyyy_1[i] * pa_x[i];

        g_xz_yyyyyz_0[i] = g_z_yyyyyz_1[i] * pa_x[i];

        g_xz_yyyyzz_0[i] = g_z_yyyyzz_1[i] * pa_x[i];

        g_xz_yyyzzz_0[i] = g_z_yyyzzz_1[i] * pa_x[i];

        g_xz_yyzzzz_0[i] = g_z_yyzzzz_1[i] * pa_x[i];

        g_xz_yzzzzz_0[i] = g_z_yzzzzz_1[i] * pa_x[i];

        g_xz_zzzzzz_0[i] = g_z_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 84-112 components of targeted buffer : DI

    auto g_yy_xxxxxx_0 = pbuffer.data(idx_eri_0_di + 84);

    auto g_yy_xxxxxy_0 = pbuffer.data(idx_eri_0_di + 85);

    auto g_yy_xxxxxz_0 = pbuffer.data(idx_eri_0_di + 86);

    auto g_yy_xxxxyy_0 = pbuffer.data(idx_eri_0_di + 87);

    auto g_yy_xxxxyz_0 = pbuffer.data(idx_eri_0_di + 88);

    auto g_yy_xxxxzz_0 = pbuffer.data(idx_eri_0_di + 89);

    auto g_yy_xxxyyy_0 = pbuffer.data(idx_eri_0_di + 90);

    auto g_yy_xxxyyz_0 = pbuffer.data(idx_eri_0_di + 91);

    auto g_yy_xxxyzz_0 = pbuffer.data(idx_eri_0_di + 92);

    auto g_yy_xxxzzz_0 = pbuffer.data(idx_eri_0_di + 93);

    auto g_yy_xxyyyy_0 = pbuffer.data(idx_eri_0_di + 94);

    auto g_yy_xxyyyz_0 = pbuffer.data(idx_eri_0_di + 95);

    auto g_yy_xxyyzz_0 = pbuffer.data(idx_eri_0_di + 96);

    auto g_yy_xxyzzz_0 = pbuffer.data(idx_eri_0_di + 97);

    auto g_yy_xxzzzz_0 = pbuffer.data(idx_eri_0_di + 98);

    auto g_yy_xyyyyy_0 = pbuffer.data(idx_eri_0_di + 99);

    auto g_yy_xyyyyz_0 = pbuffer.data(idx_eri_0_di + 100);

    auto g_yy_xyyyzz_0 = pbuffer.data(idx_eri_0_di + 101);

    auto g_yy_xyyzzz_0 = pbuffer.data(idx_eri_0_di + 102);

    auto g_yy_xyzzzz_0 = pbuffer.data(idx_eri_0_di + 103);

    auto g_yy_xzzzzz_0 = pbuffer.data(idx_eri_0_di + 104);

    auto g_yy_yyyyyy_0 = pbuffer.data(idx_eri_0_di + 105);

    auto g_yy_yyyyyz_0 = pbuffer.data(idx_eri_0_di + 106);

    auto g_yy_yyyyzz_0 = pbuffer.data(idx_eri_0_di + 107);

    auto g_yy_yyyzzz_0 = pbuffer.data(idx_eri_0_di + 108);

    auto g_yy_yyzzzz_0 = pbuffer.data(idx_eri_0_di + 109);

    auto g_yy_yzzzzz_0 = pbuffer.data(idx_eri_0_di + 110);

    auto g_yy_zzzzzz_0 = pbuffer.data(idx_eri_0_di + 111);

    #pragma omp simd aligned(g_0_xxxxxx_0, g_0_xxxxxx_1, g_0_xxxxxy_0, g_0_xxxxxy_1, g_0_xxxxxz_0, g_0_xxxxxz_1, g_0_xxxxyy_0, g_0_xxxxyy_1, g_0_xxxxyz_0, g_0_xxxxyz_1, g_0_xxxxzz_0, g_0_xxxxzz_1, g_0_xxxyyy_0, g_0_xxxyyy_1, g_0_xxxyyz_0, g_0_xxxyyz_1, g_0_xxxyzz_0, g_0_xxxyzz_1, g_0_xxxzzz_0, g_0_xxxzzz_1, g_0_xxyyyy_0, g_0_xxyyyy_1, g_0_xxyyyz_0, g_0_xxyyyz_1, g_0_xxyyzz_0, g_0_xxyyzz_1, g_0_xxyzzz_0, g_0_xxyzzz_1, g_0_xxzzzz_0, g_0_xxzzzz_1, g_0_xyyyyy_0, g_0_xyyyyy_1, g_0_xyyyyz_0, g_0_xyyyyz_1, g_0_xyyyzz_0, g_0_xyyyzz_1, g_0_xyyzzz_0, g_0_xyyzzz_1, g_0_xyzzzz_0, g_0_xyzzzz_1, g_0_xzzzzz_0, g_0_xzzzzz_1, g_0_yyyyyy_0, g_0_yyyyyy_1, g_0_yyyyyz_0, g_0_yyyyyz_1, g_0_yyyyzz_0, g_0_yyyyzz_1, g_0_yyyzzz_0, g_0_yyyzzz_1, g_0_yyzzzz_0, g_0_yyzzzz_1, g_0_yzzzzz_0, g_0_yzzzzz_1, g_0_zzzzzz_0, g_0_zzzzzz_1, g_y_xxxxx_1, g_y_xxxxxx_1, g_y_xxxxxy_1, g_y_xxxxxz_1, g_y_xxxxy_1, g_y_xxxxyy_1, g_y_xxxxyz_1, g_y_xxxxz_1, g_y_xxxxzz_1, g_y_xxxyy_1, g_y_xxxyyy_1, g_y_xxxyyz_1, g_y_xxxyz_1, g_y_xxxyzz_1, g_y_xxxzz_1, g_y_xxxzzz_1, g_y_xxyyy_1, g_y_xxyyyy_1, g_y_xxyyyz_1, g_y_xxyyz_1, g_y_xxyyzz_1, g_y_xxyzz_1, g_y_xxyzzz_1, g_y_xxzzz_1, g_y_xxzzzz_1, g_y_xyyyy_1, g_y_xyyyyy_1, g_y_xyyyyz_1, g_y_xyyyz_1, g_y_xyyyzz_1, g_y_xyyzz_1, g_y_xyyzzz_1, g_y_xyzzz_1, g_y_xyzzzz_1, g_y_xzzzz_1, g_y_xzzzzz_1, g_y_yyyyy_1, g_y_yyyyyy_1, g_y_yyyyyz_1, g_y_yyyyz_1, g_y_yyyyzz_1, g_y_yyyzz_1, g_y_yyyzzz_1, g_y_yyzzz_1, g_y_yyzzzz_1, g_y_yzzzz_1, g_y_yzzzzz_1, g_y_zzzzz_1, g_y_zzzzzz_1, g_yy_xxxxxx_0, g_yy_xxxxxy_0, g_yy_xxxxxz_0, g_yy_xxxxyy_0, g_yy_xxxxyz_0, g_yy_xxxxzz_0, g_yy_xxxyyy_0, g_yy_xxxyyz_0, g_yy_xxxyzz_0, g_yy_xxxzzz_0, g_yy_xxyyyy_0, g_yy_xxyyyz_0, g_yy_xxyyzz_0, g_yy_xxyzzz_0, g_yy_xxzzzz_0, g_yy_xyyyyy_0, g_yy_xyyyyz_0, g_yy_xyyyzz_0, g_yy_xyyzzz_0, g_yy_xyzzzz_0, g_yy_xzzzzz_0, g_yy_yyyyyy_0, g_yy_yyyyyz_0, g_yy_yyyyzz_0, g_yy_yyyzzz_0, g_yy_yyzzzz_0, g_yy_yzzzzz_0, g_yy_zzzzzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yy_xxxxxx_0[i] = g_0_xxxxxx_0[i] * fbe_0 - g_0_xxxxxx_1[i] * fz_be_0 + g_y_xxxxxx_1[i] * pa_y[i];

        g_yy_xxxxxy_0[i] = g_0_xxxxxy_0[i] * fbe_0 - g_0_xxxxxy_1[i] * fz_be_0 + g_y_xxxxx_1[i] * fe_0 + g_y_xxxxxy_1[i] * pa_y[i];

        g_yy_xxxxxz_0[i] = g_0_xxxxxz_0[i] * fbe_0 - g_0_xxxxxz_1[i] * fz_be_0 + g_y_xxxxxz_1[i] * pa_y[i];

        g_yy_xxxxyy_0[i] = g_0_xxxxyy_0[i] * fbe_0 - g_0_xxxxyy_1[i] * fz_be_0 + 2.0 * g_y_xxxxy_1[i] * fe_0 + g_y_xxxxyy_1[i] * pa_y[i];

        g_yy_xxxxyz_0[i] = g_0_xxxxyz_0[i] * fbe_0 - g_0_xxxxyz_1[i] * fz_be_0 + g_y_xxxxz_1[i] * fe_0 + g_y_xxxxyz_1[i] * pa_y[i];

        g_yy_xxxxzz_0[i] = g_0_xxxxzz_0[i] * fbe_0 - g_0_xxxxzz_1[i] * fz_be_0 + g_y_xxxxzz_1[i] * pa_y[i];

        g_yy_xxxyyy_0[i] = g_0_xxxyyy_0[i] * fbe_0 - g_0_xxxyyy_1[i] * fz_be_0 + 3.0 * g_y_xxxyy_1[i] * fe_0 + g_y_xxxyyy_1[i] * pa_y[i];

        g_yy_xxxyyz_0[i] = g_0_xxxyyz_0[i] * fbe_0 - g_0_xxxyyz_1[i] * fz_be_0 + 2.0 * g_y_xxxyz_1[i] * fe_0 + g_y_xxxyyz_1[i] * pa_y[i];

        g_yy_xxxyzz_0[i] = g_0_xxxyzz_0[i] * fbe_0 - g_0_xxxyzz_1[i] * fz_be_0 + g_y_xxxzz_1[i] * fe_0 + g_y_xxxyzz_1[i] * pa_y[i];

        g_yy_xxxzzz_0[i] = g_0_xxxzzz_0[i] * fbe_0 - g_0_xxxzzz_1[i] * fz_be_0 + g_y_xxxzzz_1[i] * pa_y[i];

        g_yy_xxyyyy_0[i] = g_0_xxyyyy_0[i] * fbe_0 - g_0_xxyyyy_1[i] * fz_be_0 + 4.0 * g_y_xxyyy_1[i] * fe_0 + g_y_xxyyyy_1[i] * pa_y[i];

        g_yy_xxyyyz_0[i] = g_0_xxyyyz_0[i] * fbe_0 - g_0_xxyyyz_1[i] * fz_be_0 + 3.0 * g_y_xxyyz_1[i] * fe_0 + g_y_xxyyyz_1[i] * pa_y[i];

        g_yy_xxyyzz_0[i] = g_0_xxyyzz_0[i] * fbe_0 - g_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_y_xxyzz_1[i] * fe_0 + g_y_xxyyzz_1[i] * pa_y[i];

        g_yy_xxyzzz_0[i] = g_0_xxyzzz_0[i] * fbe_0 - g_0_xxyzzz_1[i] * fz_be_0 + g_y_xxzzz_1[i] * fe_0 + g_y_xxyzzz_1[i] * pa_y[i];

        g_yy_xxzzzz_0[i] = g_0_xxzzzz_0[i] * fbe_0 - g_0_xxzzzz_1[i] * fz_be_0 + g_y_xxzzzz_1[i] * pa_y[i];

        g_yy_xyyyyy_0[i] = g_0_xyyyyy_0[i] * fbe_0 - g_0_xyyyyy_1[i] * fz_be_0 + 5.0 * g_y_xyyyy_1[i] * fe_0 + g_y_xyyyyy_1[i] * pa_y[i];

        g_yy_xyyyyz_0[i] = g_0_xyyyyz_0[i] * fbe_0 - g_0_xyyyyz_1[i] * fz_be_0 + 4.0 * g_y_xyyyz_1[i] * fe_0 + g_y_xyyyyz_1[i] * pa_y[i];

        g_yy_xyyyzz_0[i] = g_0_xyyyzz_0[i] * fbe_0 - g_0_xyyyzz_1[i] * fz_be_0 + 3.0 * g_y_xyyzz_1[i] * fe_0 + g_y_xyyyzz_1[i] * pa_y[i];

        g_yy_xyyzzz_0[i] = g_0_xyyzzz_0[i] * fbe_0 - g_0_xyyzzz_1[i] * fz_be_0 + 2.0 * g_y_xyzzz_1[i] * fe_0 + g_y_xyyzzz_1[i] * pa_y[i];

        g_yy_xyzzzz_0[i] = g_0_xyzzzz_0[i] * fbe_0 - g_0_xyzzzz_1[i] * fz_be_0 + g_y_xzzzz_1[i] * fe_0 + g_y_xyzzzz_1[i] * pa_y[i];

        g_yy_xzzzzz_0[i] = g_0_xzzzzz_0[i] * fbe_0 - g_0_xzzzzz_1[i] * fz_be_0 + g_y_xzzzzz_1[i] * pa_y[i];

        g_yy_yyyyyy_0[i] = g_0_yyyyyy_0[i] * fbe_0 - g_0_yyyyyy_1[i] * fz_be_0 + 6.0 * g_y_yyyyy_1[i] * fe_0 + g_y_yyyyyy_1[i] * pa_y[i];

        g_yy_yyyyyz_0[i] = g_0_yyyyyz_0[i] * fbe_0 - g_0_yyyyyz_1[i] * fz_be_0 + 5.0 * g_y_yyyyz_1[i] * fe_0 + g_y_yyyyyz_1[i] * pa_y[i];

        g_yy_yyyyzz_0[i] = g_0_yyyyzz_0[i] * fbe_0 - g_0_yyyyzz_1[i] * fz_be_0 + 4.0 * g_y_yyyzz_1[i] * fe_0 + g_y_yyyyzz_1[i] * pa_y[i];

        g_yy_yyyzzz_0[i] = g_0_yyyzzz_0[i] * fbe_0 - g_0_yyyzzz_1[i] * fz_be_0 + 3.0 * g_y_yyzzz_1[i] * fe_0 + g_y_yyyzzz_1[i] * pa_y[i];

        g_yy_yyzzzz_0[i] = g_0_yyzzzz_0[i] * fbe_0 - g_0_yyzzzz_1[i] * fz_be_0 + 2.0 * g_y_yzzzz_1[i] * fe_0 + g_y_yyzzzz_1[i] * pa_y[i];

        g_yy_yzzzzz_0[i] = g_0_yzzzzz_0[i] * fbe_0 - g_0_yzzzzz_1[i] * fz_be_0 + g_y_zzzzz_1[i] * fe_0 + g_y_yzzzzz_1[i] * pa_y[i];

        g_yy_zzzzzz_0[i] = g_0_zzzzzz_0[i] * fbe_0 - g_0_zzzzzz_1[i] * fz_be_0 + g_y_zzzzzz_1[i] * pa_y[i];
    }

    // Set up 112-140 components of targeted buffer : DI

    auto g_yz_xxxxxx_0 = pbuffer.data(idx_eri_0_di + 112);

    auto g_yz_xxxxxy_0 = pbuffer.data(idx_eri_0_di + 113);

    auto g_yz_xxxxxz_0 = pbuffer.data(idx_eri_0_di + 114);

    auto g_yz_xxxxyy_0 = pbuffer.data(idx_eri_0_di + 115);

    auto g_yz_xxxxyz_0 = pbuffer.data(idx_eri_0_di + 116);

    auto g_yz_xxxxzz_0 = pbuffer.data(idx_eri_0_di + 117);

    auto g_yz_xxxyyy_0 = pbuffer.data(idx_eri_0_di + 118);

    auto g_yz_xxxyyz_0 = pbuffer.data(idx_eri_0_di + 119);

    auto g_yz_xxxyzz_0 = pbuffer.data(idx_eri_0_di + 120);

    auto g_yz_xxxzzz_0 = pbuffer.data(idx_eri_0_di + 121);

    auto g_yz_xxyyyy_0 = pbuffer.data(idx_eri_0_di + 122);

    auto g_yz_xxyyyz_0 = pbuffer.data(idx_eri_0_di + 123);

    auto g_yz_xxyyzz_0 = pbuffer.data(idx_eri_0_di + 124);

    auto g_yz_xxyzzz_0 = pbuffer.data(idx_eri_0_di + 125);

    auto g_yz_xxzzzz_0 = pbuffer.data(idx_eri_0_di + 126);

    auto g_yz_xyyyyy_0 = pbuffer.data(idx_eri_0_di + 127);

    auto g_yz_xyyyyz_0 = pbuffer.data(idx_eri_0_di + 128);

    auto g_yz_xyyyzz_0 = pbuffer.data(idx_eri_0_di + 129);

    auto g_yz_xyyzzz_0 = pbuffer.data(idx_eri_0_di + 130);

    auto g_yz_xyzzzz_0 = pbuffer.data(idx_eri_0_di + 131);

    auto g_yz_xzzzzz_0 = pbuffer.data(idx_eri_0_di + 132);

    auto g_yz_yyyyyy_0 = pbuffer.data(idx_eri_0_di + 133);

    auto g_yz_yyyyyz_0 = pbuffer.data(idx_eri_0_di + 134);

    auto g_yz_yyyyzz_0 = pbuffer.data(idx_eri_0_di + 135);

    auto g_yz_yyyzzz_0 = pbuffer.data(idx_eri_0_di + 136);

    auto g_yz_yyzzzz_0 = pbuffer.data(idx_eri_0_di + 137);

    auto g_yz_yzzzzz_0 = pbuffer.data(idx_eri_0_di + 138);

    auto g_yz_zzzzzz_0 = pbuffer.data(idx_eri_0_di + 139);

    #pragma omp simd aligned(g_y_xxxxxy_1, g_y_xxxxyy_1, g_y_xxxyyy_1, g_y_xxyyyy_1, g_y_xyyyyy_1, g_y_yyyyyy_1, g_yz_xxxxxx_0, g_yz_xxxxxy_0, g_yz_xxxxxz_0, g_yz_xxxxyy_0, g_yz_xxxxyz_0, g_yz_xxxxzz_0, g_yz_xxxyyy_0, g_yz_xxxyyz_0, g_yz_xxxyzz_0, g_yz_xxxzzz_0, g_yz_xxyyyy_0, g_yz_xxyyyz_0, g_yz_xxyyzz_0, g_yz_xxyzzz_0, g_yz_xxzzzz_0, g_yz_xyyyyy_0, g_yz_xyyyyz_0, g_yz_xyyyzz_0, g_yz_xyyzzz_0, g_yz_xyzzzz_0, g_yz_xzzzzz_0, g_yz_yyyyyy_0, g_yz_yyyyyz_0, g_yz_yyyyzz_0, g_yz_yyyzzz_0, g_yz_yyzzzz_0, g_yz_yzzzzz_0, g_yz_zzzzzz_0, g_z_xxxxxx_1, g_z_xxxxxz_1, g_z_xxxxyz_1, g_z_xxxxz_1, g_z_xxxxzz_1, g_z_xxxyyz_1, g_z_xxxyz_1, g_z_xxxyzz_1, g_z_xxxzz_1, g_z_xxxzzz_1, g_z_xxyyyz_1, g_z_xxyyz_1, g_z_xxyyzz_1, g_z_xxyzz_1, g_z_xxyzzz_1, g_z_xxzzz_1, g_z_xxzzzz_1, g_z_xyyyyz_1, g_z_xyyyz_1, g_z_xyyyzz_1, g_z_xyyzz_1, g_z_xyyzzz_1, g_z_xyzzz_1, g_z_xyzzzz_1, g_z_xzzzz_1, g_z_xzzzzz_1, g_z_yyyyyz_1, g_z_yyyyz_1, g_z_yyyyzz_1, g_z_yyyzz_1, g_z_yyyzzz_1, g_z_yyzzz_1, g_z_yyzzzz_1, g_z_yzzzz_1, g_z_yzzzzz_1, g_z_zzzzz_1, g_z_zzzzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yz_xxxxxx_0[i] = g_z_xxxxxx_1[i] * pa_y[i];

        g_yz_xxxxxy_0[i] = g_y_xxxxxy_1[i] * pa_z[i];

        g_yz_xxxxxz_0[i] = g_z_xxxxxz_1[i] * pa_y[i];

        g_yz_xxxxyy_0[i] = g_y_xxxxyy_1[i] * pa_z[i];

        g_yz_xxxxyz_0[i] = g_z_xxxxz_1[i] * fe_0 + g_z_xxxxyz_1[i] * pa_y[i];

        g_yz_xxxxzz_0[i] = g_z_xxxxzz_1[i] * pa_y[i];

        g_yz_xxxyyy_0[i] = g_y_xxxyyy_1[i] * pa_z[i];

        g_yz_xxxyyz_0[i] = 2.0 * g_z_xxxyz_1[i] * fe_0 + g_z_xxxyyz_1[i] * pa_y[i];

        g_yz_xxxyzz_0[i] = g_z_xxxzz_1[i] * fe_0 + g_z_xxxyzz_1[i] * pa_y[i];

        g_yz_xxxzzz_0[i] = g_z_xxxzzz_1[i] * pa_y[i];

        g_yz_xxyyyy_0[i] = g_y_xxyyyy_1[i] * pa_z[i];

        g_yz_xxyyyz_0[i] = 3.0 * g_z_xxyyz_1[i] * fe_0 + g_z_xxyyyz_1[i] * pa_y[i];

        g_yz_xxyyzz_0[i] = 2.0 * g_z_xxyzz_1[i] * fe_0 + g_z_xxyyzz_1[i] * pa_y[i];

        g_yz_xxyzzz_0[i] = g_z_xxzzz_1[i] * fe_0 + g_z_xxyzzz_1[i] * pa_y[i];

        g_yz_xxzzzz_0[i] = g_z_xxzzzz_1[i] * pa_y[i];

        g_yz_xyyyyy_0[i] = g_y_xyyyyy_1[i] * pa_z[i];

        g_yz_xyyyyz_0[i] = 4.0 * g_z_xyyyz_1[i] * fe_0 + g_z_xyyyyz_1[i] * pa_y[i];

        g_yz_xyyyzz_0[i] = 3.0 * g_z_xyyzz_1[i] * fe_0 + g_z_xyyyzz_1[i] * pa_y[i];

        g_yz_xyyzzz_0[i] = 2.0 * g_z_xyzzz_1[i] * fe_0 + g_z_xyyzzz_1[i] * pa_y[i];

        g_yz_xyzzzz_0[i] = g_z_xzzzz_1[i] * fe_0 + g_z_xyzzzz_1[i] * pa_y[i];

        g_yz_xzzzzz_0[i] = g_z_xzzzzz_1[i] * pa_y[i];

        g_yz_yyyyyy_0[i] = g_y_yyyyyy_1[i] * pa_z[i];

        g_yz_yyyyyz_0[i] = 5.0 * g_z_yyyyz_1[i] * fe_0 + g_z_yyyyyz_1[i] * pa_y[i];

        g_yz_yyyyzz_0[i] = 4.0 * g_z_yyyzz_1[i] * fe_0 + g_z_yyyyzz_1[i] * pa_y[i];

        g_yz_yyyzzz_0[i] = 3.0 * g_z_yyzzz_1[i] * fe_0 + g_z_yyyzzz_1[i] * pa_y[i];

        g_yz_yyzzzz_0[i] = 2.0 * g_z_yzzzz_1[i] * fe_0 + g_z_yyzzzz_1[i] * pa_y[i];

        g_yz_yzzzzz_0[i] = g_z_zzzzz_1[i] * fe_0 + g_z_yzzzzz_1[i] * pa_y[i];

        g_yz_zzzzzz_0[i] = g_z_zzzzzz_1[i] * pa_y[i];
    }

    // Set up 140-168 components of targeted buffer : DI

    auto g_zz_xxxxxx_0 = pbuffer.data(idx_eri_0_di + 140);

    auto g_zz_xxxxxy_0 = pbuffer.data(idx_eri_0_di + 141);

    auto g_zz_xxxxxz_0 = pbuffer.data(idx_eri_0_di + 142);

    auto g_zz_xxxxyy_0 = pbuffer.data(idx_eri_0_di + 143);

    auto g_zz_xxxxyz_0 = pbuffer.data(idx_eri_0_di + 144);

    auto g_zz_xxxxzz_0 = pbuffer.data(idx_eri_0_di + 145);

    auto g_zz_xxxyyy_0 = pbuffer.data(idx_eri_0_di + 146);

    auto g_zz_xxxyyz_0 = pbuffer.data(idx_eri_0_di + 147);

    auto g_zz_xxxyzz_0 = pbuffer.data(idx_eri_0_di + 148);

    auto g_zz_xxxzzz_0 = pbuffer.data(idx_eri_0_di + 149);

    auto g_zz_xxyyyy_0 = pbuffer.data(idx_eri_0_di + 150);

    auto g_zz_xxyyyz_0 = pbuffer.data(idx_eri_0_di + 151);

    auto g_zz_xxyyzz_0 = pbuffer.data(idx_eri_0_di + 152);

    auto g_zz_xxyzzz_0 = pbuffer.data(idx_eri_0_di + 153);

    auto g_zz_xxzzzz_0 = pbuffer.data(idx_eri_0_di + 154);

    auto g_zz_xyyyyy_0 = pbuffer.data(idx_eri_0_di + 155);

    auto g_zz_xyyyyz_0 = pbuffer.data(idx_eri_0_di + 156);

    auto g_zz_xyyyzz_0 = pbuffer.data(idx_eri_0_di + 157);

    auto g_zz_xyyzzz_0 = pbuffer.data(idx_eri_0_di + 158);

    auto g_zz_xyzzzz_0 = pbuffer.data(idx_eri_0_di + 159);

    auto g_zz_xzzzzz_0 = pbuffer.data(idx_eri_0_di + 160);

    auto g_zz_yyyyyy_0 = pbuffer.data(idx_eri_0_di + 161);

    auto g_zz_yyyyyz_0 = pbuffer.data(idx_eri_0_di + 162);

    auto g_zz_yyyyzz_0 = pbuffer.data(idx_eri_0_di + 163);

    auto g_zz_yyyzzz_0 = pbuffer.data(idx_eri_0_di + 164);

    auto g_zz_yyzzzz_0 = pbuffer.data(idx_eri_0_di + 165);

    auto g_zz_yzzzzz_0 = pbuffer.data(idx_eri_0_di + 166);

    auto g_zz_zzzzzz_0 = pbuffer.data(idx_eri_0_di + 167);

    #pragma omp simd aligned(g_0_xxxxxx_0, g_0_xxxxxx_1, g_0_xxxxxy_0, g_0_xxxxxy_1, g_0_xxxxxz_0, g_0_xxxxxz_1, g_0_xxxxyy_0, g_0_xxxxyy_1, g_0_xxxxyz_0, g_0_xxxxyz_1, g_0_xxxxzz_0, g_0_xxxxzz_1, g_0_xxxyyy_0, g_0_xxxyyy_1, g_0_xxxyyz_0, g_0_xxxyyz_1, g_0_xxxyzz_0, g_0_xxxyzz_1, g_0_xxxzzz_0, g_0_xxxzzz_1, g_0_xxyyyy_0, g_0_xxyyyy_1, g_0_xxyyyz_0, g_0_xxyyyz_1, g_0_xxyyzz_0, g_0_xxyyzz_1, g_0_xxyzzz_0, g_0_xxyzzz_1, g_0_xxzzzz_0, g_0_xxzzzz_1, g_0_xyyyyy_0, g_0_xyyyyy_1, g_0_xyyyyz_0, g_0_xyyyyz_1, g_0_xyyyzz_0, g_0_xyyyzz_1, g_0_xyyzzz_0, g_0_xyyzzz_1, g_0_xyzzzz_0, g_0_xyzzzz_1, g_0_xzzzzz_0, g_0_xzzzzz_1, g_0_yyyyyy_0, g_0_yyyyyy_1, g_0_yyyyyz_0, g_0_yyyyyz_1, g_0_yyyyzz_0, g_0_yyyyzz_1, g_0_yyyzzz_0, g_0_yyyzzz_1, g_0_yyzzzz_0, g_0_yyzzzz_1, g_0_yzzzzz_0, g_0_yzzzzz_1, g_0_zzzzzz_0, g_0_zzzzzz_1, g_z_xxxxx_1, g_z_xxxxxx_1, g_z_xxxxxy_1, g_z_xxxxxz_1, g_z_xxxxy_1, g_z_xxxxyy_1, g_z_xxxxyz_1, g_z_xxxxz_1, g_z_xxxxzz_1, g_z_xxxyy_1, g_z_xxxyyy_1, g_z_xxxyyz_1, g_z_xxxyz_1, g_z_xxxyzz_1, g_z_xxxzz_1, g_z_xxxzzz_1, g_z_xxyyy_1, g_z_xxyyyy_1, g_z_xxyyyz_1, g_z_xxyyz_1, g_z_xxyyzz_1, g_z_xxyzz_1, g_z_xxyzzz_1, g_z_xxzzz_1, g_z_xxzzzz_1, g_z_xyyyy_1, g_z_xyyyyy_1, g_z_xyyyyz_1, g_z_xyyyz_1, g_z_xyyyzz_1, g_z_xyyzz_1, g_z_xyyzzz_1, g_z_xyzzz_1, g_z_xyzzzz_1, g_z_xzzzz_1, g_z_xzzzzz_1, g_z_yyyyy_1, g_z_yyyyyy_1, g_z_yyyyyz_1, g_z_yyyyz_1, g_z_yyyyzz_1, g_z_yyyzz_1, g_z_yyyzzz_1, g_z_yyzzz_1, g_z_yyzzzz_1, g_z_yzzzz_1, g_z_yzzzzz_1, g_z_zzzzz_1, g_z_zzzzzz_1, g_zz_xxxxxx_0, g_zz_xxxxxy_0, g_zz_xxxxxz_0, g_zz_xxxxyy_0, g_zz_xxxxyz_0, g_zz_xxxxzz_0, g_zz_xxxyyy_0, g_zz_xxxyyz_0, g_zz_xxxyzz_0, g_zz_xxxzzz_0, g_zz_xxyyyy_0, g_zz_xxyyyz_0, g_zz_xxyyzz_0, g_zz_xxyzzz_0, g_zz_xxzzzz_0, g_zz_xyyyyy_0, g_zz_xyyyyz_0, g_zz_xyyyzz_0, g_zz_xyyzzz_0, g_zz_xyzzzz_0, g_zz_xzzzzz_0, g_zz_yyyyyy_0, g_zz_yyyyyz_0, g_zz_yyyyzz_0, g_zz_yyyzzz_0, g_zz_yyzzzz_0, g_zz_yzzzzz_0, g_zz_zzzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_zz_xxxxxx_0[i] = g_0_xxxxxx_0[i] * fbe_0 - g_0_xxxxxx_1[i] * fz_be_0 + g_z_xxxxxx_1[i] * pa_z[i];

        g_zz_xxxxxy_0[i] = g_0_xxxxxy_0[i] * fbe_0 - g_0_xxxxxy_1[i] * fz_be_0 + g_z_xxxxxy_1[i] * pa_z[i];

        g_zz_xxxxxz_0[i] = g_0_xxxxxz_0[i] * fbe_0 - g_0_xxxxxz_1[i] * fz_be_0 + g_z_xxxxx_1[i] * fe_0 + g_z_xxxxxz_1[i] * pa_z[i];

        g_zz_xxxxyy_0[i] = g_0_xxxxyy_0[i] * fbe_0 - g_0_xxxxyy_1[i] * fz_be_0 + g_z_xxxxyy_1[i] * pa_z[i];

        g_zz_xxxxyz_0[i] = g_0_xxxxyz_0[i] * fbe_0 - g_0_xxxxyz_1[i] * fz_be_0 + g_z_xxxxy_1[i] * fe_0 + g_z_xxxxyz_1[i] * pa_z[i];

        g_zz_xxxxzz_0[i] = g_0_xxxxzz_0[i] * fbe_0 - g_0_xxxxzz_1[i] * fz_be_0 + 2.0 * g_z_xxxxz_1[i] * fe_0 + g_z_xxxxzz_1[i] * pa_z[i];

        g_zz_xxxyyy_0[i] = g_0_xxxyyy_0[i] * fbe_0 - g_0_xxxyyy_1[i] * fz_be_0 + g_z_xxxyyy_1[i] * pa_z[i];

        g_zz_xxxyyz_0[i] = g_0_xxxyyz_0[i] * fbe_0 - g_0_xxxyyz_1[i] * fz_be_0 + g_z_xxxyy_1[i] * fe_0 + g_z_xxxyyz_1[i] * pa_z[i];

        g_zz_xxxyzz_0[i] = g_0_xxxyzz_0[i] * fbe_0 - g_0_xxxyzz_1[i] * fz_be_0 + 2.0 * g_z_xxxyz_1[i] * fe_0 + g_z_xxxyzz_1[i] * pa_z[i];

        g_zz_xxxzzz_0[i] = g_0_xxxzzz_0[i] * fbe_0 - g_0_xxxzzz_1[i] * fz_be_0 + 3.0 * g_z_xxxzz_1[i] * fe_0 + g_z_xxxzzz_1[i] * pa_z[i];

        g_zz_xxyyyy_0[i] = g_0_xxyyyy_0[i] * fbe_0 - g_0_xxyyyy_1[i] * fz_be_0 + g_z_xxyyyy_1[i] * pa_z[i];

        g_zz_xxyyyz_0[i] = g_0_xxyyyz_0[i] * fbe_0 - g_0_xxyyyz_1[i] * fz_be_0 + g_z_xxyyy_1[i] * fe_0 + g_z_xxyyyz_1[i] * pa_z[i];

        g_zz_xxyyzz_0[i] = g_0_xxyyzz_0[i] * fbe_0 - g_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_z_xxyyz_1[i] * fe_0 + g_z_xxyyzz_1[i] * pa_z[i];

        g_zz_xxyzzz_0[i] = g_0_xxyzzz_0[i] * fbe_0 - g_0_xxyzzz_1[i] * fz_be_0 + 3.0 * g_z_xxyzz_1[i] * fe_0 + g_z_xxyzzz_1[i] * pa_z[i];

        g_zz_xxzzzz_0[i] = g_0_xxzzzz_0[i] * fbe_0 - g_0_xxzzzz_1[i] * fz_be_0 + 4.0 * g_z_xxzzz_1[i] * fe_0 + g_z_xxzzzz_1[i] * pa_z[i];

        g_zz_xyyyyy_0[i] = g_0_xyyyyy_0[i] * fbe_0 - g_0_xyyyyy_1[i] * fz_be_0 + g_z_xyyyyy_1[i] * pa_z[i];

        g_zz_xyyyyz_0[i] = g_0_xyyyyz_0[i] * fbe_0 - g_0_xyyyyz_1[i] * fz_be_0 + g_z_xyyyy_1[i] * fe_0 + g_z_xyyyyz_1[i] * pa_z[i];

        g_zz_xyyyzz_0[i] = g_0_xyyyzz_0[i] * fbe_0 - g_0_xyyyzz_1[i] * fz_be_0 + 2.0 * g_z_xyyyz_1[i] * fe_0 + g_z_xyyyzz_1[i] * pa_z[i];

        g_zz_xyyzzz_0[i] = g_0_xyyzzz_0[i] * fbe_0 - g_0_xyyzzz_1[i] * fz_be_0 + 3.0 * g_z_xyyzz_1[i] * fe_0 + g_z_xyyzzz_1[i] * pa_z[i];

        g_zz_xyzzzz_0[i] = g_0_xyzzzz_0[i] * fbe_0 - g_0_xyzzzz_1[i] * fz_be_0 + 4.0 * g_z_xyzzz_1[i] * fe_0 + g_z_xyzzzz_1[i] * pa_z[i];

        g_zz_xzzzzz_0[i] = g_0_xzzzzz_0[i] * fbe_0 - g_0_xzzzzz_1[i] * fz_be_0 + 5.0 * g_z_xzzzz_1[i] * fe_0 + g_z_xzzzzz_1[i] * pa_z[i];

        g_zz_yyyyyy_0[i] = g_0_yyyyyy_0[i] * fbe_0 - g_0_yyyyyy_1[i] * fz_be_0 + g_z_yyyyyy_1[i] * pa_z[i];

        g_zz_yyyyyz_0[i] = g_0_yyyyyz_0[i] * fbe_0 - g_0_yyyyyz_1[i] * fz_be_0 + g_z_yyyyy_1[i] * fe_0 + g_z_yyyyyz_1[i] * pa_z[i];

        g_zz_yyyyzz_0[i] = g_0_yyyyzz_0[i] * fbe_0 - g_0_yyyyzz_1[i] * fz_be_0 + 2.0 * g_z_yyyyz_1[i] * fe_0 + g_z_yyyyzz_1[i] * pa_z[i];

        g_zz_yyyzzz_0[i] = g_0_yyyzzz_0[i] * fbe_0 - g_0_yyyzzz_1[i] * fz_be_0 + 3.0 * g_z_yyyzz_1[i] * fe_0 + g_z_yyyzzz_1[i] * pa_z[i];

        g_zz_yyzzzz_0[i] = g_0_yyzzzz_0[i] * fbe_0 - g_0_yyzzzz_1[i] * fz_be_0 + 4.0 * g_z_yyzzz_1[i] * fe_0 + g_z_yyzzzz_1[i] * pa_z[i];

        g_zz_yzzzzz_0[i] = g_0_yzzzzz_0[i] * fbe_0 - g_0_yzzzzz_1[i] * fz_be_0 + 5.0 * g_z_yzzzz_1[i] * fe_0 + g_z_yzzzzz_1[i] * pa_z[i];

        g_zz_zzzzzz_0[i] = g_0_zzzzzz_0[i] * fbe_0 - g_0_zzzzzz_1[i] * fz_be_0 + 6.0 * g_z_zzzzz_1[i] * fe_0 + g_z_zzzzzz_1[i] * pa_z[i];
    }

}

} // t2ceri namespace


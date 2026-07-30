#include "ProjectedCorePotentialPrimRecDIForS.hpp"

namespace t2pecp { // t2pecp namespace

auto
comp_prim_projected_core_potential_di_s(CSimdArray<double>& pbuffer, 
                                        const size_t idx_di_s_0_0_0,
                                        const size_t idx_si_s_0_0_0,
                                        const size_t idx_pi_s_0_0_0,
                                        const size_t idx_si_s_1_0_0,
                                        const size_t idx_pi_s_1_0_0,
                                        const int p,
                                        const size_t idx_si_s_0_0_1,
                                        const size_t idx_pi_s_0_0_1,
                                        const CSimdArray<double>& factors,
                                        const TPoint<double>& r_a,
                                        const double a_exp,
                                        const double c_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents on ket side

    auto b_exps = factors.data(0);

    // Set up A center coordinates

    const auto xyz = r_a.coordinates();

    const auto a_x = xyz[0];

    const auto a_y = xyz[1];

    const auto a_z = xyz[2];

    // Set up components of auxiliary buffer : SI

    auto tg_0_xxxxxx_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0);

    auto tg_0_xxxxxy_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 1);

    auto tg_0_xxxxxz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 2);

    auto tg_0_xxxxyy_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 3);

    auto tg_0_xxxxyz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 4);

    auto tg_0_xxxxzz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 5);

    auto tg_0_xxxyyy_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 6);

    auto tg_0_xxxyyz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 7);

    auto tg_0_xxxyzz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 8);

    auto tg_0_xxxzzz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 9);

    auto tg_0_xxyyyy_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 10);

    auto tg_0_xxyyyz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 11);

    auto tg_0_xxyyzz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 12);

    auto tg_0_xxyzzz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 13);

    auto tg_0_xxzzzz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 14);

    auto tg_0_xyyyyy_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 15);

    auto tg_0_xyyyyz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 16);

    auto tg_0_xyyyzz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 17);

    auto tg_0_xyyzzz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 18);

    auto tg_0_xyzzzz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 19);

    auto tg_0_xzzzzz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 20);

    auto tg_0_yyyyyy_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 21);

    auto tg_0_yyyyyz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 22);

    auto tg_0_yyyyzz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 23);

    auto tg_0_yyyzzz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 24);

    auto tg_0_yyzzzz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 25);

    auto tg_0_yzzzzz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 26);

    auto tg_0_zzzzzz_s_0_0_0 = pbuffer.data(idx_si_s_0_0_0 + 27);

    // Set up components of auxiliary buffer : PI

    auto tg_x_xxxxxx_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0);

    auto tg_x_xxxxxy_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 1);

    auto tg_x_xxxxxz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 2);

    auto tg_x_xxxxyy_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 3);

    auto tg_x_xxxxyz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 4);

    auto tg_x_xxxxzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 5);

    auto tg_x_xxxyyy_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 6);

    auto tg_x_xxxyyz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 7);

    auto tg_x_xxxyzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 8);

    auto tg_x_xxxzzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 9);

    auto tg_x_xxyyyy_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 10);

    auto tg_x_xxyyyz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 11);

    auto tg_x_xxyyzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 12);

    auto tg_x_xxyzzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 13);

    auto tg_x_xxzzzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 14);

    auto tg_x_xyyyyy_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 15);

    auto tg_x_xyyyyz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 16);

    auto tg_x_xyyyzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 17);

    auto tg_x_xyyzzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 18);

    auto tg_x_xyzzzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 19);

    auto tg_x_xzzzzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 20);

    auto tg_x_yyyyyy_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 21);

    auto tg_x_yyyyyz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 22);

    auto tg_x_yyyyzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 23);

    auto tg_x_yyyzzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 24);

    auto tg_x_yyzzzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 25);

    auto tg_x_yzzzzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 26);

    auto tg_x_zzzzzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 27);

    auto tg_y_xxxxxx_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 28);

    auto tg_y_xxxxxy_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 29);

    auto tg_y_xxxxxz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 30);

    auto tg_y_xxxxyy_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 31);

    auto tg_y_xxxxyz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 32);

    auto tg_y_xxxxzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 33);

    auto tg_y_xxxyyy_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 34);

    auto tg_y_xxxyyz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 35);

    auto tg_y_xxxyzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 36);

    auto tg_y_xxxzzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 37);

    auto tg_y_xxyyyy_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 38);

    auto tg_y_xxyyyz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 39);

    auto tg_y_xxyyzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 40);

    auto tg_y_xxyzzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 41);

    auto tg_y_xxzzzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 42);

    auto tg_y_xyyyyy_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 43);

    auto tg_y_xyyyyz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 44);

    auto tg_y_xyyyzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 45);

    auto tg_y_xyyzzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 46);

    auto tg_y_xyzzzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 47);

    auto tg_y_xzzzzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 48);

    auto tg_y_yyyyyy_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 49);

    auto tg_y_yyyyyz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 50);

    auto tg_y_yyyyzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 51);

    auto tg_y_yyyzzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 52);

    auto tg_y_yyzzzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 53);

    auto tg_y_yzzzzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 54);

    auto tg_y_zzzzzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 55);

    auto tg_z_xxxxxx_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 56);

    auto tg_z_xxxxxy_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 57);

    auto tg_z_xxxxxz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 58);

    auto tg_z_xxxxyy_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 59);

    auto tg_z_xxxxyz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 60);

    auto tg_z_xxxxzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 61);

    auto tg_z_xxxyyy_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 62);

    auto tg_z_xxxyyz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 63);

    auto tg_z_xxxyzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 64);

    auto tg_z_xxxzzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 65);

    auto tg_z_xxyyyy_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 66);

    auto tg_z_xxyyyz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 67);

    auto tg_z_xxyyzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 68);

    auto tg_z_xxyzzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 69);

    auto tg_z_xxzzzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 70);

    auto tg_z_xyyyyy_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 71);

    auto tg_z_xyyyyz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 72);

    auto tg_z_xyyyzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 73);

    auto tg_z_xyyzzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 74);

    auto tg_z_xyzzzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 75);

    auto tg_z_xzzzzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 76);

    auto tg_z_yyyyyy_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 77);

    auto tg_z_yyyyyz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 78);

    auto tg_z_yyyyzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 79);

    auto tg_z_yyyzzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 80);

    auto tg_z_yyzzzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 81);

    auto tg_z_yzzzzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 82);

    auto tg_z_zzzzzz_s_0_0_0 = pbuffer.data(idx_pi_s_0_0_0 + 83);

    // Set up components of auxiliary buffer : SI

    auto tg_0_xxxxxx_s_1_0_0 = pbuffer.data(idx_si_s_1_0_0);

    auto tg_0_xxxxxy_s_1_0_0 = pbuffer.data(idx_si_s_1_0_0 + 1);

    auto tg_0_xxxxxz_s_1_0_0 = pbuffer.data(idx_si_s_1_0_0 + 2);

    auto tg_0_xxxxyy_s_1_0_0 = pbuffer.data(idx_si_s_1_0_0 + 3);

    auto tg_0_xxxxyz_s_1_0_0 = pbuffer.data(idx_si_s_1_0_0 + 4);

    auto tg_0_xxxxzz_s_1_0_0 = pbuffer.data(idx_si_s_1_0_0 + 5);

    auto tg_0_xxxyyy_s_1_0_0 = pbuffer.data(idx_si_s_1_0_0 + 6);

    auto tg_0_xxxyyz_s_1_0_0 = pbuffer.data(idx_si_s_1_0_0 + 7);

    auto tg_0_xxxyzz_s_1_0_0 = pbuffer.data(idx_si_s_1_0_0 + 8);

    auto tg_0_xxxzzz_s_1_0_0 = pbuffer.data(idx_si_s_1_0_0 + 9);

    auto tg_0_xxyyyy_s_1_0_0 = pbuffer.data(idx_si_s_1_0_0 + 10);

    auto tg_0_xxyyyz_s_1_0_0 = pbuffer.data(idx_si_s_1_0_0 + 11);

    auto tg_0_xxyyzz_s_1_0_0 = pbuffer.data(idx_si_s_1_0_0 + 12);

    auto tg_0_xxyzzz_s_1_0_0 = pbuffer.data(idx_si_s_1_0_0 + 13);

    auto tg_0_xxzzzz_s_1_0_0 = pbuffer.data(idx_si_s_1_0_0 + 14);

    auto tg_0_xyyyyy_s_1_0_0 = pbuffer.data(idx_si_s_1_0_0 + 15);

    auto tg_0_xyyyyz_s_1_0_0 = pbuffer.data(idx_si_s_1_0_0 + 16);

    auto tg_0_xyyyzz_s_1_0_0 = pbuffer.data(idx_si_s_1_0_0 + 17);

    auto tg_0_xyyzzz_s_1_0_0 = pbuffer.data(idx_si_s_1_0_0 + 18);

    auto tg_0_xyzzzz_s_1_0_0 = pbuffer.data(idx_si_s_1_0_0 + 19);

    auto tg_0_xzzzzz_s_1_0_0 = pbuffer.data(idx_si_s_1_0_0 + 20);

    auto tg_0_yyyyyy_s_1_0_0 = pbuffer.data(idx_si_s_1_0_0 + 21);

    auto tg_0_yyyyyz_s_1_0_0 = pbuffer.data(idx_si_s_1_0_0 + 22);

    auto tg_0_yyyyzz_s_1_0_0 = pbuffer.data(idx_si_s_1_0_0 + 23);

    auto tg_0_yyyzzz_s_1_0_0 = pbuffer.data(idx_si_s_1_0_0 + 24);

    auto tg_0_yyzzzz_s_1_0_0 = pbuffer.data(idx_si_s_1_0_0 + 25);

    auto tg_0_yzzzzz_s_1_0_0 = pbuffer.data(idx_si_s_1_0_0 + 26);

    auto tg_0_zzzzzz_s_1_0_0 = pbuffer.data(idx_si_s_1_0_0 + 27);

    // Set up components of auxiliary buffer : PI

    auto tg_x_xxxxxx_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0);

    auto tg_x_xxxxxy_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 1);

    auto tg_x_xxxxxz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 2);

    auto tg_x_xxxxyy_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 3);

    auto tg_x_xxxxyz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 4);

    auto tg_x_xxxxzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 5);

    auto tg_x_xxxyyy_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 6);

    auto tg_x_xxxyyz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 7);

    auto tg_x_xxxyzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 8);

    auto tg_x_xxxzzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 9);

    auto tg_x_xxyyyy_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 10);

    auto tg_x_xxyyyz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 11);

    auto tg_x_xxyyzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 12);

    auto tg_x_xxyzzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 13);

    auto tg_x_xxzzzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 14);

    auto tg_x_xyyyyy_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 15);

    auto tg_x_xyyyyz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 16);

    auto tg_x_xyyyzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 17);

    auto tg_x_xyyzzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 18);

    auto tg_x_xyzzzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 19);

    auto tg_x_xzzzzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 20);

    auto tg_x_yyyyyy_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 21);

    auto tg_x_yyyyyz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 22);

    auto tg_x_yyyyzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 23);

    auto tg_x_yyyzzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 24);

    auto tg_x_yyzzzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 25);

    auto tg_x_yzzzzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 26);

    auto tg_x_zzzzzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 27);

    auto tg_y_xxxxxx_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 28);

    auto tg_y_xxxxxy_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 29);

    auto tg_y_xxxxxz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 30);

    auto tg_y_xxxxyy_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 31);

    auto tg_y_xxxxyz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 32);

    auto tg_y_xxxxzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 33);

    auto tg_y_xxxyyy_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 34);

    auto tg_y_xxxyyz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 35);

    auto tg_y_xxxyzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 36);

    auto tg_y_xxxzzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 37);

    auto tg_y_xxyyyy_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 38);

    auto tg_y_xxyyyz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 39);

    auto tg_y_xxyyzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 40);

    auto tg_y_xxyzzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 41);

    auto tg_y_xxzzzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 42);

    auto tg_y_xyyyyy_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 43);

    auto tg_y_xyyyyz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 44);

    auto tg_y_xyyyzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 45);

    auto tg_y_xyyzzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 46);

    auto tg_y_xyzzzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 47);

    auto tg_y_xzzzzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 48);

    auto tg_y_yyyyyy_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 49);

    auto tg_y_yyyyyz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 50);

    auto tg_y_yyyyzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 51);

    auto tg_y_yyyzzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 52);

    auto tg_y_yyzzzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 53);

    auto tg_y_yzzzzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 54);

    auto tg_y_zzzzzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 55);

    auto tg_z_xxxxxx_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 56);

    auto tg_z_xxxxxy_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 57);

    auto tg_z_xxxxxz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 58);

    auto tg_z_xxxxyy_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 59);

    auto tg_z_xxxxyz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 60);

    auto tg_z_xxxxzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 61);

    auto tg_z_xxxyyy_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 62);

    auto tg_z_xxxyyz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 63);

    auto tg_z_xxxyzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 64);

    auto tg_z_xxxzzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 65);

    auto tg_z_xxyyyy_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 66);

    auto tg_z_xxyyyz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 67);

    auto tg_z_xxyyzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 68);

    auto tg_z_xxyzzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 69);

    auto tg_z_xxzzzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 70);

    auto tg_z_xyyyyy_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 71);

    auto tg_z_xyyyyz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 72);

    auto tg_z_xyyyzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 73);

    auto tg_z_xyyzzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 74);

    auto tg_z_xyzzzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 75);

    auto tg_z_xzzzzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 76);

    auto tg_z_yyyyyy_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 77);

    auto tg_z_yyyyyz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 78);

    auto tg_z_yyyyzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 79);

    auto tg_z_yyyzzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 80);

    auto tg_z_yyzzzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 81);

    auto tg_z_yzzzzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 82);

    auto tg_z_zzzzzz_s_1_0_0 = pbuffer.data(idx_pi_s_1_0_0 + 83);

    // Set up components of targeted buffer : DI

    auto tg_xx_xxxxxx_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0);

    auto tg_xx_xxxxxy_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 1);

    auto tg_xx_xxxxxz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 2);

    auto tg_xx_xxxxyy_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 3);

    auto tg_xx_xxxxyz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 4);

    auto tg_xx_xxxxzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 5);

    auto tg_xx_xxxyyy_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 6);

    auto tg_xx_xxxyyz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 7);

    auto tg_xx_xxxyzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 8);

    auto tg_xx_xxxzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 9);

    auto tg_xx_xxyyyy_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 10);

    auto tg_xx_xxyyyz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 11);

    auto tg_xx_xxyyzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 12);

    auto tg_xx_xxyzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 13);

    auto tg_xx_xxzzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 14);

    auto tg_xx_xyyyyy_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 15);

    auto tg_xx_xyyyyz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 16);

    auto tg_xx_xyyyzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 17);

    auto tg_xx_xyyzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 18);

    auto tg_xx_xyzzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 19);

    auto tg_xx_xzzzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 20);

    auto tg_xx_yyyyyy_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 21);

    auto tg_xx_yyyyyz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 22);

    auto tg_xx_yyyyzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 23);

    auto tg_xx_yyyzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 24);

    auto tg_xx_yyzzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 25);

    auto tg_xx_yzzzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 26);

    auto tg_xx_zzzzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 27);

    auto tg_xy_xxxxxx_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 28);

    auto tg_xy_xxxxxy_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 29);

    auto tg_xy_xxxxxz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 30);

    auto tg_xy_xxxxyy_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 31);

    auto tg_xy_xxxxyz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 32);

    auto tg_xy_xxxxzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 33);

    auto tg_xy_xxxyyy_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 34);

    auto tg_xy_xxxyyz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 35);

    auto tg_xy_xxxyzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 36);

    auto tg_xy_xxxzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 37);

    auto tg_xy_xxyyyy_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 38);

    auto tg_xy_xxyyyz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 39);

    auto tg_xy_xxyyzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 40);

    auto tg_xy_xxyzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 41);

    auto tg_xy_xxzzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 42);

    auto tg_xy_xyyyyy_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 43);

    auto tg_xy_xyyyyz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 44);

    auto tg_xy_xyyyzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 45);

    auto tg_xy_xyyzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 46);

    auto tg_xy_xyzzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 47);

    auto tg_xy_xzzzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 48);

    auto tg_xy_yyyyyy_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 49);

    auto tg_xy_yyyyyz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 50);

    auto tg_xy_yyyyzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 51);

    auto tg_xy_yyyzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 52);

    auto tg_xy_yyzzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 53);

    auto tg_xy_yzzzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 54);

    auto tg_xy_zzzzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 55);

    auto tg_xz_xxxxxx_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 56);

    auto tg_xz_xxxxxy_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 57);

    auto tg_xz_xxxxxz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 58);

    auto tg_xz_xxxxyy_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 59);

    auto tg_xz_xxxxyz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 60);

    auto tg_xz_xxxxzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 61);

    auto tg_xz_xxxyyy_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 62);

    auto tg_xz_xxxyyz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 63);

    auto tg_xz_xxxyzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 64);

    auto tg_xz_xxxzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 65);

    auto tg_xz_xxyyyy_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 66);

    auto tg_xz_xxyyyz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 67);

    auto tg_xz_xxyyzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 68);

    auto tg_xz_xxyzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 69);

    auto tg_xz_xxzzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 70);

    auto tg_xz_xyyyyy_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 71);

    auto tg_xz_xyyyyz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 72);

    auto tg_xz_xyyyzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 73);

    auto tg_xz_xyyzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 74);

    auto tg_xz_xyzzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 75);

    auto tg_xz_xzzzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 76);

    auto tg_xz_yyyyyy_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 77);

    auto tg_xz_yyyyyz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 78);

    auto tg_xz_yyyyzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 79);

    auto tg_xz_yyyzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 80);

    auto tg_xz_yyzzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 81);

    auto tg_xz_yzzzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 82);

    auto tg_xz_zzzzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 83);

    auto tg_yy_xxxxxx_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 84);

    auto tg_yy_xxxxxy_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 85);

    auto tg_yy_xxxxxz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 86);

    auto tg_yy_xxxxyy_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 87);

    auto tg_yy_xxxxyz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 88);

    auto tg_yy_xxxxzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 89);

    auto tg_yy_xxxyyy_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 90);

    auto tg_yy_xxxyyz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 91);

    auto tg_yy_xxxyzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 92);

    auto tg_yy_xxxzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 93);

    auto tg_yy_xxyyyy_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 94);

    auto tg_yy_xxyyyz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 95);

    auto tg_yy_xxyyzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 96);

    auto tg_yy_xxyzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 97);

    auto tg_yy_xxzzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 98);

    auto tg_yy_xyyyyy_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 99);

    auto tg_yy_xyyyyz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 100);

    auto tg_yy_xyyyzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 101);

    auto tg_yy_xyyzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 102);

    auto tg_yy_xyzzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 103);

    auto tg_yy_xzzzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 104);

    auto tg_yy_yyyyyy_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 105);

    auto tg_yy_yyyyyz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 106);

    auto tg_yy_yyyyzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 107);

    auto tg_yy_yyyzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 108);

    auto tg_yy_yyzzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 109);

    auto tg_yy_yzzzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 110);

    auto tg_yy_zzzzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 111);

    auto tg_yz_xxxxxx_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 112);

    auto tg_yz_xxxxxy_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 113);

    auto tg_yz_xxxxxz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 114);

    auto tg_yz_xxxxyy_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 115);

    auto tg_yz_xxxxyz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 116);

    auto tg_yz_xxxxzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 117);

    auto tg_yz_xxxyyy_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 118);

    auto tg_yz_xxxyyz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 119);

    auto tg_yz_xxxyzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 120);

    auto tg_yz_xxxzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 121);

    auto tg_yz_xxyyyy_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 122);

    auto tg_yz_xxyyyz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 123);

    auto tg_yz_xxyyzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 124);

    auto tg_yz_xxyzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 125);

    auto tg_yz_xxzzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 126);

    auto tg_yz_xyyyyy_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 127);

    auto tg_yz_xyyyyz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 128);

    auto tg_yz_xyyyzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 129);

    auto tg_yz_xyyzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 130);

    auto tg_yz_xyzzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 131);

    auto tg_yz_xzzzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 132);

    auto tg_yz_yyyyyy_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 133);

    auto tg_yz_yyyyyz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 134);

    auto tg_yz_yyyyzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 135);

    auto tg_yz_yyyzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 136);

    auto tg_yz_yyzzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 137);

    auto tg_yz_yzzzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 138);

    auto tg_yz_zzzzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 139);

    auto tg_zz_xxxxxx_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 140);

    auto tg_zz_xxxxxy_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 141);

    auto tg_zz_xxxxxz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 142);

    auto tg_zz_xxxxyy_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 143);

    auto tg_zz_xxxxyz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 144);

    auto tg_zz_xxxxzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 145);

    auto tg_zz_xxxyyy_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 146);

    auto tg_zz_xxxyyz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 147);

    auto tg_zz_xxxyzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 148);

    auto tg_zz_xxxzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 149);

    auto tg_zz_xxyyyy_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 150);

    auto tg_zz_xxyyyz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 151);

    auto tg_zz_xxyyzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 152);

    auto tg_zz_xxyzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 153);

    auto tg_zz_xxzzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 154);

    auto tg_zz_xyyyyy_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 155);

    auto tg_zz_xyyyyz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 156);

    auto tg_zz_xyyyzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 157);

    auto tg_zz_xyyzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 158);

    auto tg_zz_xyzzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 159);

    auto tg_zz_xzzzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 160);

    auto tg_zz_yyyyyy_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 161);

    auto tg_zz_yyyyyz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 162);

    auto tg_zz_yyyyzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 163);

    auto tg_zz_yyyzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 164);

    auto tg_zz_yyzzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 165);

    auto tg_zz_yzzzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 166);

    auto tg_zz_zzzzzz_s_0_0_0 = pbuffer.data(idx_di_s_0_0_0 + 167);

    #pragma omp simd aligned(b_exps, tg_0_xxxxxx_s_0_0_0, tg_0_xxxxxx_s_1_0_0, tg_0_xxxxxy_s_0_0_0, tg_0_xxxxxy_s_1_0_0, tg_0_xxxxxz_s_0_0_0, tg_0_xxxxxz_s_1_0_0, tg_0_xxxxyy_s_0_0_0, tg_0_xxxxyy_s_1_0_0, tg_0_xxxxyz_s_0_0_0, tg_0_xxxxyz_s_1_0_0, tg_0_xxxxzz_s_0_0_0, tg_0_xxxxzz_s_1_0_0, tg_0_xxxyyy_s_0_0_0, tg_0_xxxyyy_s_1_0_0, tg_0_xxxyyz_s_0_0_0, tg_0_xxxyyz_s_1_0_0, tg_0_xxxyzz_s_0_0_0, tg_0_xxxyzz_s_1_0_0, tg_0_xxxzzz_s_0_0_0, tg_0_xxxzzz_s_1_0_0, tg_0_xxyyyy_s_0_0_0, tg_0_xxyyyy_s_1_0_0, tg_0_xxyyyz_s_0_0_0, tg_0_xxyyyz_s_1_0_0, tg_0_xxyyzz_s_0_0_0, tg_0_xxyyzz_s_1_0_0, tg_0_xxyzzz_s_0_0_0, tg_0_xxyzzz_s_1_0_0, tg_0_xxzzzz_s_0_0_0, tg_0_xxzzzz_s_1_0_0, tg_0_xyyyyy_s_0_0_0, tg_0_xyyyyy_s_1_0_0, tg_0_xyyyyz_s_0_0_0, tg_0_xyyyyz_s_1_0_0, tg_0_xyyyzz_s_0_0_0, tg_0_xyyyzz_s_1_0_0, tg_0_xyyzzz_s_0_0_0, tg_0_xyyzzz_s_1_0_0, tg_0_xyzzzz_s_0_0_0, tg_0_xyzzzz_s_1_0_0, tg_0_xzzzzz_s_0_0_0, tg_0_xzzzzz_s_1_0_0, tg_0_yyyyyy_s_0_0_0, tg_0_yyyyyy_s_1_0_0, tg_0_yyyyyz_s_0_0_0, tg_0_yyyyyz_s_1_0_0, tg_0_yyyyzz_s_0_0_0, tg_0_yyyyzz_s_1_0_0, tg_0_yyyzzz_s_0_0_0, tg_0_yyyzzz_s_1_0_0, tg_0_yyzzzz_s_0_0_0, tg_0_yyzzzz_s_1_0_0, tg_0_yzzzzz_s_0_0_0, tg_0_yzzzzz_s_1_0_0, tg_0_zzzzzz_s_0_0_0, tg_0_zzzzzz_s_1_0_0, tg_x_xxxxxx_s_0_0_0, tg_x_xxxxxx_s_1_0_0, tg_x_xxxxxy_s_0_0_0, tg_x_xxxxxy_s_1_0_0, tg_x_xxxxxz_s_0_0_0, tg_x_xxxxxz_s_1_0_0, tg_x_xxxxyy_s_0_0_0, tg_x_xxxxyy_s_1_0_0, tg_x_xxxxyz_s_0_0_0, tg_x_xxxxyz_s_1_0_0, tg_x_xxxxzz_s_0_0_0, tg_x_xxxxzz_s_1_0_0, tg_x_xxxyyy_s_0_0_0, tg_x_xxxyyy_s_1_0_0, tg_x_xxxyyz_s_0_0_0, tg_x_xxxyyz_s_1_0_0, tg_x_xxxyzz_s_0_0_0, tg_x_xxxyzz_s_1_0_0, tg_x_xxxzzz_s_0_0_0, tg_x_xxxzzz_s_1_0_0, tg_x_xxyyyy_s_0_0_0, tg_x_xxyyyy_s_1_0_0, tg_x_xxyyyz_s_0_0_0, tg_x_xxyyyz_s_1_0_0, tg_x_xxyyzz_s_0_0_0, tg_x_xxyyzz_s_1_0_0, tg_x_xxyzzz_s_0_0_0, tg_x_xxyzzz_s_1_0_0, tg_x_xxzzzz_s_0_0_0, tg_x_xxzzzz_s_1_0_0, tg_x_xyyyyy_s_0_0_0, tg_x_xyyyyy_s_1_0_0, tg_x_xyyyyz_s_0_0_0, tg_x_xyyyyz_s_1_0_0, tg_x_xyyyzz_s_0_0_0, tg_x_xyyyzz_s_1_0_0, tg_x_xyyzzz_s_0_0_0, tg_x_xyyzzz_s_1_0_0, tg_x_xyzzzz_s_0_0_0, tg_x_xyzzzz_s_1_0_0, tg_x_xzzzzz_s_0_0_0, tg_x_xzzzzz_s_1_0_0, tg_x_yyyyyy_s_0_0_0, tg_x_yyyyyy_s_1_0_0, tg_x_yyyyyz_s_0_0_0, tg_x_yyyyyz_s_1_0_0, tg_x_yyyyzz_s_0_0_0, tg_x_yyyyzz_s_1_0_0, tg_x_yyyzzz_s_0_0_0, tg_x_yyyzzz_s_1_0_0, tg_x_yyzzzz_s_0_0_0, tg_x_yyzzzz_s_1_0_0, tg_x_yzzzzz_s_0_0_0, tg_x_yzzzzz_s_1_0_0, tg_x_zzzzzz_s_0_0_0, tg_x_zzzzzz_s_1_0_0, tg_xx_xxxxxx_s_0_0_0, tg_xx_xxxxxy_s_0_0_0, tg_xx_xxxxxz_s_0_0_0, tg_xx_xxxxyy_s_0_0_0, tg_xx_xxxxyz_s_0_0_0, tg_xx_xxxxzz_s_0_0_0, tg_xx_xxxyyy_s_0_0_0, tg_xx_xxxyyz_s_0_0_0, tg_xx_xxxyzz_s_0_0_0, tg_xx_xxxzzz_s_0_0_0, tg_xx_xxyyyy_s_0_0_0, tg_xx_xxyyyz_s_0_0_0, tg_xx_xxyyzz_s_0_0_0, tg_xx_xxyzzz_s_0_0_0, tg_xx_xxzzzz_s_0_0_0, tg_xx_xyyyyy_s_0_0_0, tg_xx_xyyyyz_s_0_0_0, tg_xx_xyyyzz_s_0_0_0, tg_xx_xyyzzz_s_0_0_0, tg_xx_xyzzzz_s_0_0_0, tg_xx_xzzzzz_s_0_0_0, tg_xx_yyyyyy_s_0_0_0, tg_xx_yyyyyz_s_0_0_0, tg_xx_yyyyzz_s_0_0_0, tg_xx_yyyzzz_s_0_0_0, tg_xx_yyzzzz_s_0_0_0, tg_xx_yzzzzz_s_0_0_0, tg_xx_zzzzzz_s_0_0_0, tg_xy_xxxxxx_s_0_0_0, tg_xy_xxxxxy_s_0_0_0, tg_xy_xxxxxz_s_0_0_0, tg_xy_xxxxyy_s_0_0_0, tg_xy_xxxxyz_s_0_0_0, tg_xy_xxxxzz_s_0_0_0, tg_xy_xxxyyy_s_0_0_0, tg_xy_xxxyyz_s_0_0_0, tg_xy_xxxyzz_s_0_0_0, tg_xy_xxxzzz_s_0_0_0, tg_xy_xxyyyy_s_0_0_0, tg_xy_xxyyyz_s_0_0_0, tg_xy_xxyyzz_s_0_0_0, tg_xy_xxyzzz_s_0_0_0, tg_xy_xxzzzz_s_0_0_0, tg_xy_xyyyyy_s_0_0_0, tg_xy_xyyyyz_s_0_0_0, tg_xy_xyyyzz_s_0_0_0, tg_xy_xyyzzz_s_0_0_0, tg_xy_xyzzzz_s_0_0_0, tg_xy_xzzzzz_s_0_0_0, tg_xy_yyyyyy_s_0_0_0, tg_xy_yyyyyz_s_0_0_0, tg_xy_yyyyzz_s_0_0_0, tg_xy_yyyzzz_s_0_0_0, tg_xy_yyzzzz_s_0_0_0, tg_xy_yzzzzz_s_0_0_0, tg_xy_zzzzzz_s_0_0_0, tg_xz_xxxxxx_s_0_0_0, tg_xz_xxxxxy_s_0_0_0, tg_xz_xxxxxz_s_0_0_0, tg_xz_xxxxyy_s_0_0_0, tg_xz_xxxxyz_s_0_0_0, tg_xz_xxxxzz_s_0_0_0, tg_xz_xxxyyy_s_0_0_0, tg_xz_xxxyyz_s_0_0_0, tg_xz_xxxyzz_s_0_0_0, tg_xz_xxxzzz_s_0_0_0, tg_xz_xxyyyy_s_0_0_0, tg_xz_xxyyyz_s_0_0_0, tg_xz_xxyyzz_s_0_0_0, tg_xz_xxyzzz_s_0_0_0, tg_xz_xxzzzz_s_0_0_0, tg_xz_xyyyyy_s_0_0_0, tg_xz_xyyyyz_s_0_0_0, tg_xz_xyyyzz_s_0_0_0, tg_xz_xyyzzz_s_0_0_0, tg_xz_xyzzzz_s_0_0_0, tg_xz_xzzzzz_s_0_0_0, tg_xz_yyyyyy_s_0_0_0, tg_xz_yyyyyz_s_0_0_0, tg_xz_yyyyzz_s_0_0_0, tg_xz_yyyzzz_s_0_0_0, tg_xz_yyzzzz_s_0_0_0, tg_xz_yzzzzz_s_0_0_0, tg_xz_zzzzzz_s_0_0_0, tg_y_xxxxxx_s_0_0_0, tg_y_xxxxxx_s_1_0_0, tg_y_xxxxxy_s_0_0_0, tg_y_xxxxxy_s_1_0_0, tg_y_xxxxxz_s_0_0_0, tg_y_xxxxxz_s_1_0_0, tg_y_xxxxyy_s_0_0_0, tg_y_xxxxyy_s_1_0_0, tg_y_xxxxyz_s_0_0_0, tg_y_xxxxyz_s_1_0_0, tg_y_xxxxzz_s_0_0_0, tg_y_xxxxzz_s_1_0_0, tg_y_xxxyyy_s_0_0_0, tg_y_xxxyyy_s_1_0_0, tg_y_xxxyyz_s_0_0_0, tg_y_xxxyyz_s_1_0_0, tg_y_xxxyzz_s_0_0_0, tg_y_xxxyzz_s_1_0_0, tg_y_xxxzzz_s_0_0_0, tg_y_xxxzzz_s_1_0_0, tg_y_xxyyyy_s_0_0_0, tg_y_xxyyyy_s_1_0_0, tg_y_xxyyyz_s_0_0_0, tg_y_xxyyyz_s_1_0_0, tg_y_xxyyzz_s_0_0_0, tg_y_xxyyzz_s_1_0_0, tg_y_xxyzzz_s_0_0_0, tg_y_xxyzzz_s_1_0_0, tg_y_xxzzzz_s_0_0_0, tg_y_xxzzzz_s_1_0_0, tg_y_xyyyyy_s_0_0_0, tg_y_xyyyyy_s_1_0_0, tg_y_xyyyyz_s_0_0_0, tg_y_xyyyyz_s_1_0_0, tg_y_xyyyzz_s_0_0_0, tg_y_xyyyzz_s_1_0_0, tg_y_xyyzzz_s_0_0_0, tg_y_xyyzzz_s_1_0_0, tg_y_xyzzzz_s_0_0_0, tg_y_xyzzzz_s_1_0_0, tg_y_xzzzzz_s_0_0_0, tg_y_xzzzzz_s_1_0_0, tg_y_yyyyyy_s_0_0_0, tg_y_yyyyyy_s_1_0_0, tg_y_yyyyyz_s_0_0_0, tg_y_yyyyyz_s_1_0_0, tg_y_yyyyzz_s_0_0_0, tg_y_yyyyzz_s_1_0_0, tg_y_yyyzzz_s_0_0_0, tg_y_yyyzzz_s_1_0_0, tg_y_yyzzzz_s_0_0_0, tg_y_yyzzzz_s_1_0_0, tg_y_yzzzzz_s_0_0_0, tg_y_yzzzzz_s_1_0_0, tg_y_zzzzzz_s_0_0_0, tg_y_zzzzzz_s_1_0_0, tg_yy_xxxxxx_s_0_0_0, tg_yy_xxxxxy_s_0_0_0, tg_yy_xxxxxz_s_0_0_0, tg_yy_xxxxyy_s_0_0_0, tg_yy_xxxxyz_s_0_0_0, tg_yy_xxxxzz_s_0_0_0, tg_yy_xxxyyy_s_0_0_0, tg_yy_xxxyyz_s_0_0_0, tg_yy_xxxyzz_s_0_0_0, tg_yy_xxxzzz_s_0_0_0, tg_yy_xxyyyy_s_0_0_0, tg_yy_xxyyyz_s_0_0_0, tg_yy_xxyyzz_s_0_0_0, tg_yy_xxyzzz_s_0_0_0, tg_yy_xxzzzz_s_0_0_0, tg_yy_xyyyyy_s_0_0_0, tg_yy_xyyyyz_s_0_0_0, tg_yy_xyyyzz_s_0_0_0, tg_yy_xyyzzz_s_0_0_0, tg_yy_xyzzzz_s_0_0_0, tg_yy_xzzzzz_s_0_0_0, tg_yy_yyyyyy_s_0_0_0, tg_yy_yyyyyz_s_0_0_0, tg_yy_yyyyzz_s_0_0_0, tg_yy_yyyzzz_s_0_0_0, tg_yy_yyzzzz_s_0_0_0, tg_yy_yzzzzz_s_0_0_0, tg_yy_zzzzzz_s_0_0_0, tg_yz_xxxxxx_s_0_0_0, tg_yz_xxxxxy_s_0_0_0, tg_yz_xxxxxz_s_0_0_0, tg_yz_xxxxyy_s_0_0_0, tg_yz_xxxxyz_s_0_0_0, tg_yz_xxxxzz_s_0_0_0, tg_yz_xxxyyy_s_0_0_0, tg_yz_xxxyyz_s_0_0_0, tg_yz_xxxyzz_s_0_0_0, tg_yz_xxxzzz_s_0_0_0, tg_yz_xxyyyy_s_0_0_0, tg_yz_xxyyyz_s_0_0_0, tg_yz_xxyyzz_s_0_0_0, tg_yz_xxyzzz_s_0_0_0, tg_yz_xxzzzz_s_0_0_0, tg_yz_xyyyyy_s_0_0_0, tg_yz_xyyyyz_s_0_0_0, tg_yz_xyyyzz_s_0_0_0, tg_yz_xyyzzz_s_0_0_0, tg_yz_xyzzzz_s_0_0_0, tg_yz_xzzzzz_s_0_0_0, tg_yz_yyyyyy_s_0_0_0, tg_yz_yyyyyz_s_0_0_0, tg_yz_yyyyzz_s_0_0_0, tg_yz_yyyzzz_s_0_0_0, tg_yz_yyzzzz_s_0_0_0, tg_yz_yzzzzz_s_0_0_0, tg_yz_zzzzzz_s_0_0_0, tg_z_xxxxxx_s_0_0_0, tg_z_xxxxxx_s_1_0_0, tg_z_xxxxxy_s_0_0_0, tg_z_xxxxxy_s_1_0_0, tg_z_xxxxxz_s_0_0_0, tg_z_xxxxxz_s_1_0_0, tg_z_xxxxyy_s_0_0_0, tg_z_xxxxyy_s_1_0_0, tg_z_xxxxyz_s_0_0_0, tg_z_xxxxyz_s_1_0_0, tg_z_xxxxzz_s_0_0_0, tg_z_xxxxzz_s_1_0_0, tg_z_xxxyyy_s_0_0_0, tg_z_xxxyyy_s_1_0_0, tg_z_xxxyyz_s_0_0_0, tg_z_xxxyyz_s_1_0_0, tg_z_xxxyzz_s_0_0_0, tg_z_xxxyzz_s_1_0_0, tg_z_xxxzzz_s_0_0_0, tg_z_xxxzzz_s_1_0_0, tg_z_xxyyyy_s_0_0_0, tg_z_xxyyyy_s_1_0_0, tg_z_xxyyyz_s_0_0_0, tg_z_xxyyyz_s_1_0_0, tg_z_xxyyzz_s_0_0_0, tg_z_xxyyzz_s_1_0_0, tg_z_xxyzzz_s_0_0_0, tg_z_xxyzzz_s_1_0_0, tg_z_xxzzzz_s_0_0_0, tg_z_xxzzzz_s_1_0_0, tg_z_xyyyyy_s_0_0_0, tg_z_xyyyyy_s_1_0_0, tg_z_xyyyyz_s_0_0_0, tg_z_xyyyyz_s_1_0_0, tg_z_xyyyzz_s_0_0_0, tg_z_xyyyzz_s_1_0_0, tg_z_xyyzzz_s_0_0_0, tg_z_xyyzzz_s_1_0_0, tg_z_xyzzzz_s_0_0_0, tg_z_xyzzzz_s_1_0_0, tg_z_xzzzzz_s_0_0_0, tg_z_xzzzzz_s_1_0_0, tg_z_yyyyyy_s_0_0_0, tg_z_yyyyyy_s_1_0_0, tg_z_yyyyyz_s_0_0_0, tg_z_yyyyyz_s_1_0_0, tg_z_yyyyzz_s_0_0_0, tg_z_yyyyzz_s_1_0_0, tg_z_yyyzzz_s_0_0_0, tg_z_yyyzzz_s_1_0_0, tg_z_yyzzzz_s_0_0_0, tg_z_yyzzzz_s_1_0_0, tg_z_yzzzzz_s_0_0_0, tg_z_yzzzzz_s_1_0_0, tg_z_zzzzzz_s_0_0_0, tg_z_zzzzzz_s_1_0_0, tg_zz_xxxxxx_s_0_0_0, tg_zz_xxxxxy_s_0_0_0, tg_zz_xxxxxz_s_0_0_0, tg_zz_xxxxyy_s_0_0_0, tg_zz_xxxxyz_s_0_0_0, tg_zz_xxxxzz_s_0_0_0, tg_zz_xxxyyy_s_0_0_0, tg_zz_xxxyyz_s_0_0_0, tg_zz_xxxyzz_s_0_0_0, tg_zz_xxxzzz_s_0_0_0, tg_zz_xxyyyy_s_0_0_0, tg_zz_xxyyyz_s_0_0_0, tg_zz_xxyyzz_s_0_0_0, tg_zz_xxyzzz_s_0_0_0, tg_zz_xxzzzz_s_0_0_0, tg_zz_xyyyyy_s_0_0_0, tg_zz_xyyyyz_s_0_0_0, tg_zz_xyyyzz_s_0_0_0, tg_zz_xyyzzz_s_0_0_0, tg_zz_xyzzzz_s_0_0_0, tg_zz_xzzzzz_s_0_0_0, tg_zz_yyyyyy_s_0_0_0, tg_zz_yyyyyz_s_0_0_0, tg_zz_yyyyzz_s_0_0_0, tg_zz_yyyzzz_s_0_0_0, tg_zz_yyzzzz_s_0_0_0, tg_zz_yzzzzz_s_0_0_0, tg_zz_zzzzzz_s_0_0_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fzi_0 = 1.0 / (a_exp + b_exps[i] + c_exp);

        const double faz_0 = -(b_exps[i] + c_exp) * fzi_0;

        const double fbzi_0 = b_exps[i] * fzi_0;

        tg_xx_xxxxxx_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxxxxx_s_0_0_0[i] * fzi_0 + tg_0_xxxxxx_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_xxxxxx_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xxxxxx_s_0_0_0[i] * a_x * faz_0;

        tg_xx_xxxxxy_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxxxxy_s_0_0_0[i] * fzi_0 + tg_0_xxxxxy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_xxxxxy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xxxxxy_s_0_0_0[i] * a_x * faz_0;

        tg_xx_xxxxxz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxxxxz_s_0_0_0[i] * fzi_0 + tg_0_xxxxxz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_xxxxxz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xxxxxz_s_0_0_0[i] * a_x * faz_0;

        tg_xx_xxxxyy_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxxxyy_s_0_0_0[i] * fzi_0 + tg_0_xxxxyy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_xxxxyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xxxxyy_s_0_0_0[i] * a_x * faz_0;

        tg_xx_xxxxyz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxxxyz_s_0_0_0[i] * fzi_0 + tg_0_xxxxyz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_xxxxyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xxxxyz_s_0_0_0[i] * a_x * faz_0;

        tg_xx_xxxxzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxxxzz_s_0_0_0[i] * fzi_0 + tg_0_xxxxzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_xxxxzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xxxxzz_s_0_0_0[i] * a_x * faz_0;

        tg_xx_xxxyyy_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxxyyy_s_0_0_0[i] * fzi_0 + tg_0_xxxyyy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_xxxyyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xxxyyy_s_0_0_0[i] * a_x * faz_0;

        tg_xx_xxxyyz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxxyyz_s_0_0_0[i] * fzi_0 + tg_0_xxxyyz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_xxxyyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xxxyyz_s_0_0_0[i] * a_x * faz_0;

        tg_xx_xxxyzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxxyzz_s_0_0_0[i] * fzi_0 + tg_0_xxxyzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_xxxyzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xxxyzz_s_0_0_0[i] * a_x * faz_0;

        tg_xx_xxxzzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxxzzz_s_0_0_0[i] * fzi_0 + tg_0_xxxzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_xxxzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xxxzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xx_xxyyyy_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxyyyy_s_0_0_0[i] * fzi_0 + tg_0_xxyyyy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_xxyyyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xxyyyy_s_0_0_0[i] * a_x * faz_0;

        tg_xx_xxyyyz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxyyyz_s_0_0_0[i] * fzi_0 + tg_0_xxyyyz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_xxyyyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xxyyyz_s_0_0_0[i] * a_x * faz_0;

        tg_xx_xxyyzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxyyzz_s_0_0_0[i] * fzi_0 + tg_0_xxyyzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_xxyyzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xxyyzz_s_0_0_0[i] * a_x * faz_0;

        tg_xx_xxyzzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxyzzz_s_0_0_0[i] * fzi_0 + tg_0_xxyzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_xxyzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xxyzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xx_xxzzzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxzzzz_s_0_0_0[i] * fzi_0 + tg_0_xxzzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_xxzzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xxzzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xx_xyyyyy_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xyyyyy_s_0_0_0[i] * fzi_0 + tg_0_xyyyyy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_xyyyyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xyyyyy_s_0_0_0[i] * a_x * faz_0;

        tg_xx_xyyyyz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xyyyyz_s_0_0_0[i] * fzi_0 + tg_0_xyyyyz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_xyyyyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xyyyyz_s_0_0_0[i] * a_x * faz_0;

        tg_xx_xyyyzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xyyyzz_s_0_0_0[i] * fzi_0 + tg_0_xyyyzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_xyyyzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xyyyzz_s_0_0_0[i] * a_x * faz_0;

        tg_xx_xyyzzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xyyzzz_s_0_0_0[i] * fzi_0 + tg_0_xyyzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_xyyzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xyyzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xx_xyzzzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xyzzzz_s_0_0_0[i] * fzi_0 + tg_0_xyzzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_xyzzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xyzzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xx_xzzzzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xzzzzz_s_0_0_0[i] * fzi_0 + tg_0_xzzzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_xzzzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xzzzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xx_yyyyyy_s_0_0_0[i] = 1.0 / 2.0 * tg_0_yyyyyy_s_0_0_0[i] * fzi_0 + tg_0_yyyyyy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_yyyyyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_yyyyyy_s_0_0_0[i] * a_x * faz_0;

        tg_xx_yyyyyz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_yyyyyz_s_0_0_0[i] * fzi_0 + tg_0_yyyyyz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_yyyyyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_yyyyyz_s_0_0_0[i] * a_x * faz_0;

        tg_xx_yyyyzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_yyyyzz_s_0_0_0[i] * fzi_0 + tg_0_yyyyzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_yyyyzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_yyyyzz_s_0_0_0[i] * a_x * faz_0;

        tg_xx_yyyzzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_yyyzzz_s_0_0_0[i] * fzi_0 + tg_0_yyyzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_yyyzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_yyyzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xx_yyzzzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_yyzzzz_s_0_0_0[i] * fzi_0 + tg_0_yyzzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_yyzzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_yyzzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xx_yzzzzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_yzzzzz_s_0_0_0[i] * fzi_0 + tg_0_yzzzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_yzzzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_yzzzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xx_zzzzzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_zzzzzz_s_0_0_0[i] * fzi_0 + tg_0_zzzzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_zzzzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_zzzzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xy_xxxxxx_s_0_0_0[i] = 2.0 * tg_y_xxxxxx_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xxxxxx_s_0_0_0[i] * a_x * faz_0;

        tg_xy_xxxxxy_s_0_0_0[i] = 2.0 * tg_y_xxxxxy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xxxxxy_s_0_0_0[i] * a_x * faz_0;

        tg_xy_xxxxxz_s_0_0_0[i] = 2.0 * tg_y_xxxxxz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xxxxxz_s_0_0_0[i] * a_x * faz_0;

        tg_xy_xxxxyy_s_0_0_0[i] = 2.0 * tg_y_xxxxyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xxxxyy_s_0_0_0[i] * a_x * faz_0;

        tg_xy_xxxxyz_s_0_0_0[i] = 2.0 * tg_y_xxxxyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xxxxyz_s_0_0_0[i] * a_x * faz_0;

        tg_xy_xxxxzz_s_0_0_0[i] = 2.0 * tg_y_xxxxzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xxxxzz_s_0_0_0[i] * a_x * faz_0;

        tg_xy_xxxyyy_s_0_0_0[i] = 2.0 * tg_y_xxxyyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xxxyyy_s_0_0_0[i] * a_x * faz_0;

        tg_xy_xxxyyz_s_0_0_0[i] = 2.0 * tg_y_xxxyyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xxxyyz_s_0_0_0[i] * a_x * faz_0;

        tg_xy_xxxyzz_s_0_0_0[i] = 2.0 * tg_y_xxxyzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xxxyzz_s_0_0_0[i] * a_x * faz_0;

        tg_xy_xxxzzz_s_0_0_0[i] = 2.0 * tg_y_xxxzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xxxzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xy_xxyyyy_s_0_0_0[i] = 2.0 * tg_y_xxyyyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xxyyyy_s_0_0_0[i] * a_x * faz_0;

        tg_xy_xxyyyz_s_0_0_0[i] = 2.0 * tg_y_xxyyyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xxyyyz_s_0_0_0[i] * a_x * faz_0;

        tg_xy_xxyyzz_s_0_0_0[i] = 2.0 * tg_y_xxyyzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xxyyzz_s_0_0_0[i] * a_x * faz_0;

        tg_xy_xxyzzz_s_0_0_0[i] = 2.0 * tg_y_xxyzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xxyzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xy_xxzzzz_s_0_0_0[i] = 2.0 * tg_y_xxzzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xxzzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xy_xyyyyy_s_0_0_0[i] = 2.0 * tg_y_xyyyyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xyyyyy_s_0_0_0[i] * a_x * faz_0;

        tg_xy_xyyyyz_s_0_0_0[i] = 2.0 * tg_y_xyyyyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xyyyyz_s_0_0_0[i] * a_x * faz_0;

        tg_xy_xyyyzz_s_0_0_0[i] = 2.0 * tg_y_xyyyzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xyyyzz_s_0_0_0[i] * a_x * faz_0;

        tg_xy_xyyzzz_s_0_0_0[i] = 2.0 * tg_y_xyyzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xyyzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xy_xyzzzz_s_0_0_0[i] = 2.0 * tg_y_xyzzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xyzzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xy_xzzzzz_s_0_0_0[i] = 2.0 * tg_y_xzzzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xzzzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xy_yyyyyy_s_0_0_0[i] = 2.0 * tg_y_yyyyyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_yyyyyy_s_0_0_0[i] * a_x * faz_0;

        tg_xy_yyyyyz_s_0_0_0[i] = 2.0 * tg_y_yyyyyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_yyyyyz_s_0_0_0[i] * a_x * faz_0;

        tg_xy_yyyyzz_s_0_0_0[i] = 2.0 * tg_y_yyyyzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_yyyyzz_s_0_0_0[i] * a_x * faz_0;

        tg_xy_yyyzzz_s_0_0_0[i] = 2.0 * tg_y_yyyzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_yyyzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xy_yyzzzz_s_0_0_0[i] = 2.0 * tg_y_yyzzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_yyzzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xy_yzzzzz_s_0_0_0[i] = 2.0 * tg_y_yzzzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_yzzzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xy_zzzzzz_s_0_0_0[i] = 2.0 * tg_y_zzzzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_zzzzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xz_xxxxxx_s_0_0_0[i] = 2.0 * tg_z_xxxxxx_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxxxx_s_0_0_0[i] * a_x * faz_0;

        tg_xz_xxxxxy_s_0_0_0[i] = 2.0 * tg_z_xxxxxy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxxxy_s_0_0_0[i] * a_x * faz_0;

        tg_xz_xxxxxz_s_0_0_0[i] = 2.0 * tg_z_xxxxxz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxxxz_s_0_0_0[i] * a_x * faz_0;

        tg_xz_xxxxyy_s_0_0_0[i] = 2.0 * tg_z_xxxxyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxxyy_s_0_0_0[i] * a_x * faz_0;

        tg_xz_xxxxyz_s_0_0_0[i] = 2.0 * tg_z_xxxxyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxxyz_s_0_0_0[i] * a_x * faz_0;

        tg_xz_xxxxzz_s_0_0_0[i] = 2.0 * tg_z_xxxxzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxxzz_s_0_0_0[i] * a_x * faz_0;

        tg_xz_xxxyyy_s_0_0_0[i] = 2.0 * tg_z_xxxyyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxyyy_s_0_0_0[i] * a_x * faz_0;

        tg_xz_xxxyyz_s_0_0_0[i] = 2.0 * tg_z_xxxyyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxyyz_s_0_0_0[i] * a_x * faz_0;

        tg_xz_xxxyzz_s_0_0_0[i] = 2.0 * tg_z_xxxyzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxyzz_s_0_0_0[i] * a_x * faz_0;

        tg_xz_xxxzzz_s_0_0_0[i] = 2.0 * tg_z_xxxzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xz_xxyyyy_s_0_0_0[i] = 2.0 * tg_z_xxyyyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xxyyyy_s_0_0_0[i] * a_x * faz_0;

        tg_xz_xxyyyz_s_0_0_0[i] = 2.0 * tg_z_xxyyyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xxyyyz_s_0_0_0[i] * a_x * faz_0;

        tg_xz_xxyyzz_s_0_0_0[i] = 2.0 * tg_z_xxyyzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xxyyzz_s_0_0_0[i] * a_x * faz_0;

        tg_xz_xxyzzz_s_0_0_0[i] = 2.0 * tg_z_xxyzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xxyzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xz_xxzzzz_s_0_0_0[i] = 2.0 * tg_z_xxzzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xxzzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xz_xyyyyy_s_0_0_0[i] = 2.0 * tg_z_xyyyyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xyyyyy_s_0_0_0[i] * a_x * faz_0;

        tg_xz_xyyyyz_s_0_0_0[i] = 2.0 * tg_z_xyyyyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xyyyyz_s_0_0_0[i] * a_x * faz_0;

        tg_xz_xyyyzz_s_0_0_0[i] = 2.0 * tg_z_xyyyzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xyyyzz_s_0_0_0[i] * a_x * faz_0;

        tg_xz_xyyzzz_s_0_0_0[i] = 2.0 * tg_z_xyyzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xyyzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xz_xyzzzz_s_0_0_0[i] = 2.0 * tg_z_xyzzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xyzzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xz_xzzzzz_s_0_0_0[i] = 2.0 * tg_z_xzzzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xzzzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xz_yyyyyy_s_0_0_0[i] = 2.0 * tg_z_yyyyyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_yyyyyy_s_0_0_0[i] * a_x * faz_0;

        tg_xz_yyyyyz_s_0_0_0[i] = 2.0 * tg_z_yyyyyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_yyyyyz_s_0_0_0[i] * a_x * faz_0;

        tg_xz_yyyyzz_s_0_0_0[i] = 2.0 * tg_z_yyyyzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_yyyyzz_s_0_0_0[i] * a_x * faz_0;

        tg_xz_yyyzzz_s_0_0_0[i] = 2.0 * tg_z_yyyzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_yyyzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xz_yyzzzz_s_0_0_0[i] = 2.0 * tg_z_yyzzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_yyzzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xz_yzzzzz_s_0_0_0[i] = 2.0 * tg_z_yzzzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_yzzzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xz_zzzzzz_s_0_0_0[i] = 2.0 * tg_z_zzzzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_zzzzzz_s_0_0_0[i] * a_x * faz_0;

        tg_yy_xxxxxx_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxxxxx_s_0_0_0[i] * fzi_0 + tg_0_xxxxxx_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_xxxxxx_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xxxxxx_s_0_0_0[i] * a_y * faz_0;

        tg_yy_xxxxxy_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxxxxy_s_0_0_0[i] * fzi_0 + tg_0_xxxxxy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_xxxxxy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xxxxxy_s_0_0_0[i] * a_y * faz_0;

        tg_yy_xxxxxz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxxxxz_s_0_0_0[i] * fzi_0 + tg_0_xxxxxz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_xxxxxz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xxxxxz_s_0_0_0[i] * a_y * faz_0;

        tg_yy_xxxxyy_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxxxyy_s_0_0_0[i] * fzi_0 + tg_0_xxxxyy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_xxxxyy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xxxxyy_s_0_0_0[i] * a_y * faz_0;

        tg_yy_xxxxyz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxxxyz_s_0_0_0[i] * fzi_0 + tg_0_xxxxyz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_xxxxyz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xxxxyz_s_0_0_0[i] * a_y * faz_0;

        tg_yy_xxxxzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxxxzz_s_0_0_0[i] * fzi_0 + tg_0_xxxxzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_xxxxzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xxxxzz_s_0_0_0[i] * a_y * faz_0;

        tg_yy_xxxyyy_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxxyyy_s_0_0_0[i] * fzi_0 + tg_0_xxxyyy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_xxxyyy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xxxyyy_s_0_0_0[i] * a_y * faz_0;

        tg_yy_xxxyyz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxxyyz_s_0_0_0[i] * fzi_0 + tg_0_xxxyyz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_xxxyyz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xxxyyz_s_0_0_0[i] * a_y * faz_0;

        tg_yy_xxxyzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxxyzz_s_0_0_0[i] * fzi_0 + tg_0_xxxyzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_xxxyzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xxxyzz_s_0_0_0[i] * a_y * faz_0;

        tg_yy_xxxzzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxxzzz_s_0_0_0[i] * fzi_0 + tg_0_xxxzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_xxxzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xxxzzz_s_0_0_0[i] * a_y * faz_0;

        tg_yy_xxyyyy_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxyyyy_s_0_0_0[i] * fzi_0 + tg_0_xxyyyy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_xxyyyy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xxyyyy_s_0_0_0[i] * a_y * faz_0;

        tg_yy_xxyyyz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxyyyz_s_0_0_0[i] * fzi_0 + tg_0_xxyyyz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_xxyyyz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xxyyyz_s_0_0_0[i] * a_y * faz_0;

        tg_yy_xxyyzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxyyzz_s_0_0_0[i] * fzi_0 + tg_0_xxyyzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_xxyyzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xxyyzz_s_0_0_0[i] * a_y * faz_0;

        tg_yy_xxyzzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxyzzz_s_0_0_0[i] * fzi_0 + tg_0_xxyzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_xxyzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xxyzzz_s_0_0_0[i] * a_y * faz_0;

        tg_yy_xxzzzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxzzzz_s_0_0_0[i] * fzi_0 + tg_0_xxzzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_xxzzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xxzzzz_s_0_0_0[i] * a_y * faz_0;

        tg_yy_xyyyyy_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xyyyyy_s_0_0_0[i] * fzi_0 + tg_0_xyyyyy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_xyyyyy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xyyyyy_s_0_0_0[i] * a_y * faz_0;

        tg_yy_xyyyyz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xyyyyz_s_0_0_0[i] * fzi_0 + tg_0_xyyyyz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_xyyyyz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xyyyyz_s_0_0_0[i] * a_y * faz_0;

        tg_yy_xyyyzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xyyyzz_s_0_0_0[i] * fzi_0 + tg_0_xyyyzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_xyyyzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xyyyzz_s_0_0_0[i] * a_y * faz_0;

        tg_yy_xyyzzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xyyzzz_s_0_0_0[i] * fzi_0 + tg_0_xyyzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_xyyzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xyyzzz_s_0_0_0[i] * a_y * faz_0;

        tg_yy_xyzzzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xyzzzz_s_0_0_0[i] * fzi_0 + tg_0_xyzzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_xyzzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xyzzzz_s_0_0_0[i] * a_y * faz_0;

        tg_yy_xzzzzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xzzzzz_s_0_0_0[i] * fzi_0 + tg_0_xzzzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_xzzzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xzzzzz_s_0_0_0[i] * a_y * faz_0;

        tg_yy_yyyyyy_s_0_0_0[i] = 1.0 / 2.0 * tg_0_yyyyyy_s_0_0_0[i] * fzi_0 + tg_0_yyyyyy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_yyyyyy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_yyyyyy_s_0_0_0[i] * a_y * faz_0;

        tg_yy_yyyyyz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_yyyyyz_s_0_0_0[i] * fzi_0 + tg_0_yyyyyz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_yyyyyz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_yyyyyz_s_0_0_0[i] * a_y * faz_0;

        tg_yy_yyyyzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_yyyyzz_s_0_0_0[i] * fzi_0 + tg_0_yyyyzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_yyyyzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_yyyyzz_s_0_0_0[i] * a_y * faz_0;

        tg_yy_yyyzzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_yyyzzz_s_0_0_0[i] * fzi_0 + tg_0_yyyzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_yyyzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_yyyzzz_s_0_0_0[i] * a_y * faz_0;

        tg_yy_yyzzzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_yyzzzz_s_0_0_0[i] * fzi_0 + tg_0_yyzzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_yyzzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_yyzzzz_s_0_0_0[i] * a_y * faz_0;

        tg_yy_yzzzzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_yzzzzz_s_0_0_0[i] * fzi_0 + tg_0_yzzzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_yzzzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_yzzzzz_s_0_0_0[i] * a_y * faz_0;

        tg_yy_zzzzzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_zzzzzz_s_0_0_0[i] * fzi_0 + tg_0_zzzzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_zzzzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_zzzzzz_s_0_0_0[i] * a_y * faz_0;

        tg_yz_xxxxxx_s_0_0_0[i] = 2.0 * tg_z_xxxxxx_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxxxx_s_0_0_0[i] * a_y * faz_0;

        tg_yz_xxxxxy_s_0_0_0[i] = 2.0 * tg_z_xxxxxy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxxxy_s_0_0_0[i] * a_y * faz_0;

        tg_yz_xxxxxz_s_0_0_0[i] = 2.0 * tg_z_xxxxxz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxxxz_s_0_0_0[i] * a_y * faz_0;

        tg_yz_xxxxyy_s_0_0_0[i] = 2.0 * tg_z_xxxxyy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxxyy_s_0_0_0[i] * a_y * faz_0;

        tg_yz_xxxxyz_s_0_0_0[i] = 2.0 * tg_z_xxxxyz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxxyz_s_0_0_0[i] * a_y * faz_0;

        tg_yz_xxxxzz_s_0_0_0[i] = 2.0 * tg_z_xxxxzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxxzz_s_0_0_0[i] * a_y * faz_0;

        tg_yz_xxxyyy_s_0_0_0[i] = 2.0 * tg_z_xxxyyy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxyyy_s_0_0_0[i] * a_y * faz_0;

        tg_yz_xxxyyz_s_0_0_0[i] = 2.0 * tg_z_xxxyyz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxyyz_s_0_0_0[i] * a_y * faz_0;

        tg_yz_xxxyzz_s_0_0_0[i] = 2.0 * tg_z_xxxyzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxyzz_s_0_0_0[i] * a_y * faz_0;

        tg_yz_xxxzzz_s_0_0_0[i] = 2.0 * tg_z_xxxzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxzzz_s_0_0_0[i] * a_y * faz_0;

        tg_yz_xxyyyy_s_0_0_0[i] = 2.0 * tg_z_xxyyyy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xxyyyy_s_0_0_0[i] * a_y * faz_0;

        tg_yz_xxyyyz_s_0_0_0[i] = 2.0 * tg_z_xxyyyz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xxyyyz_s_0_0_0[i] * a_y * faz_0;

        tg_yz_xxyyzz_s_0_0_0[i] = 2.0 * tg_z_xxyyzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xxyyzz_s_0_0_0[i] * a_y * faz_0;

        tg_yz_xxyzzz_s_0_0_0[i] = 2.0 * tg_z_xxyzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xxyzzz_s_0_0_0[i] * a_y * faz_0;

        tg_yz_xxzzzz_s_0_0_0[i] = 2.0 * tg_z_xxzzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xxzzzz_s_0_0_0[i] * a_y * faz_0;

        tg_yz_xyyyyy_s_0_0_0[i] = 2.0 * tg_z_xyyyyy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xyyyyy_s_0_0_0[i] * a_y * faz_0;

        tg_yz_xyyyyz_s_0_0_0[i] = 2.0 * tg_z_xyyyyz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xyyyyz_s_0_0_0[i] * a_y * faz_0;

        tg_yz_xyyyzz_s_0_0_0[i] = 2.0 * tg_z_xyyyzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xyyyzz_s_0_0_0[i] * a_y * faz_0;

        tg_yz_xyyzzz_s_0_0_0[i] = 2.0 * tg_z_xyyzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xyyzzz_s_0_0_0[i] * a_y * faz_0;

        tg_yz_xyzzzz_s_0_0_0[i] = 2.0 * tg_z_xyzzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xyzzzz_s_0_0_0[i] * a_y * faz_0;

        tg_yz_xzzzzz_s_0_0_0[i] = 2.0 * tg_z_xzzzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xzzzzz_s_0_0_0[i] * a_y * faz_0;

        tg_yz_yyyyyy_s_0_0_0[i] = 2.0 * tg_z_yyyyyy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_yyyyyy_s_0_0_0[i] * a_y * faz_0;

        tg_yz_yyyyyz_s_0_0_0[i] = 2.0 * tg_z_yyyyyz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_yyyyyz_s_0_0_0[i] * a_y * faz_0;

        tg_yz_yyyyzz_s_0_0_0[i] = 2.0 * tg_z_yyyyzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_yyyyzz_s_0_0_0[i] * a_y * faz_0;

        tg_yz_yyyzzz_s_0_0_0[i] = 2.0 * tg_z_yyyzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_yyyzzz_s_0_0_0[i] * a_y * faz_0;

        tg_yz_yyzzzz_s_0_0_0[i] = 2.0 * tg_z_yyzzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_yyzzzz_s_0_0_0[i] * a_y * faz_0;

        tg_yz_yzzzzz_s_0_0_0[i] = 2.0 * tg_z_yzzzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_yzzzzz_s_0_0_0[i] * a_y * faz_0;

        tg_yz_zzzzzz_s_0_0_0[i] = 2.0 * tg_z_zzzzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_zzzzzz_s_0_0_0[i] * a_y * faz_0;

        tg_zz_xxxxxx_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxxxxx_s_0_0_0[i] * fzi_0 + tg_0_xxxxxx_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_xxxxxx_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxxxx_s_0_0_0[i] * a_z * faz_0;

        tg_zz_xxxxxy_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxxxxy_s_0_0_0[i] * fzi_0 + tg_0_xxxxxy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_xxxxxy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxxxy_s_0_0_0[i] * a_z * faz_0;

        tg_zz_xxxxxz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxxxxz_s_0_0_0[i] * fzi_0 + tg_0_xxxxxz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_xxxxxz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxxxz_s_0_0_0[i] * a_z * faz_0;

        tg_zz_xxxxyy_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxxxyy_s_0_0_0[i] * fzi_0 + tg_0_xxxxyy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_xxxxyy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxxyy_s_0_0_0[i] * a_z * faz_0;

        tg_zz_xxxxyz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxxxyz_s_0_0_0[i] * fzi_0 + tg_0_xxxxyz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_xxxxyz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxxyz_s_0_0_0[i] * a_z * faz_0;

        tg_zz_xxxxzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxxxzz_s_0_0_0[i] * fzi_0 + tg_0_xxxxzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_xxxxzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxxzz_s_0_0_0[i] * a_z * faz_0;

        tg_zz_xxxyyy_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxxyyy_s_0_0_0[i] * fzi_0 + tg_0_xxxyyy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_xxxyyy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxyyy_s_0_0_0[i] * a_z * faz_0;

        tg_zz_xxxyyz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxxyyz_s_0_0_0[i] * fzi_0 + tg_0_xxxyyz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_xxxyyz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxyyz_s_0_0_0[i] * a_z * faz_0;

        tg_zz_xxxyzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxxyzz_s_0_0_0[i] * fzi_0 + tg_0_xxxyzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_xxxyzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxyzz_s_0_0_0[i] * a_z * faz_0;

        tg_zz_xxxzzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxxzzz_s_0_0_0[i] * fzi_0 + tg_0_xxxzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_xxxzzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxzzz_s_0_0_0[i] * a_z * faz_0;

        tg_zz_xxyyyy_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxyyyy_s_0_0_0[i] * fzi_0 + tg_0_xxyyyy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_xxyyyy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xxyyyy_s_0_0_0[i] * a_z * faz_0;

        tg_zz_xxyyyz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxyyyz_s_0_0_0[i] * fzi_0 + tg_0_xxyyyz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_xxyyyz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xxyyyz_s_0_0_0[i] * a_z * faz_0;

        tg_zz_xxyyzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxyyzz_s_0_0_0[i] * fzi_0 + tg_0_xxyyzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_xxyyzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xxyyzz_s_0_0_0[i] * a_z * faz_0;

        tg_zz_xxyzzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxyzzz_s_0_0_0[i] * fzi_0 + tg_0_xxyzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_xxyzzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xxyzzz_s_0_0_0[i] * a_z * faz_0;

        tg_zz_xxzzzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xxzzzz_s_0_0_0[i] * fzi_0 + tg_0_xxzzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_xxzzzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xxzzzz_s_0_0_0[i] * a_z * faz_0;

        tg_zz_xyyyyy_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xyyyyy_s_0_0_0[i] * fzi_0 + tg_0_xyyyyy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_xyyyyy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xyyyyy_s_0_0_0[i] * a_z * faz_0;

        tg_zz_xyyyyz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xyyyyz_s_0_0_0[i] * fzi_0 + tg_0_xyyyyz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_xyyyyz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xyyyyz_s_0_0_0[i] * a_z * faz_0;

        tg_zz_xyyyzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xyyyzz_s_0_0_0[i] * fzi_0 + tg_0_xyyyzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_xyyyzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xyyyzz_s_0_0_0[i] * a_z * faz_0;

        tg_zz_xyyzzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xyyzzz_s_0_0_0[i] * fzi_0 + tg_0_xyyzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_xyyzzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xyyzzz_s_0_0_0[i] * a_z * faz_0;

        tg_zz_xyzzzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xyzzzz_s_0_0_0[i] * fzi_0 + tg_0_xyzzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_xyzzzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xyzzzz_s_0_0_0[i] * a_z * faz_0;

        tg_zz_xzzzzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_xzzzzz_s_0_0_0[i] * fzi_0 + tg_0_xzzzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_xzzzzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xzzzzz_s_0_0_0[i] * a_z * faz_0;

        tg_zz_yyyyyy_s_0_0_0[i] = 1.0 / 2.0 * tg_0_yyyyyy_s_0_0_0[i] * fzi_0 + tg_0_yyyyyy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_yyyyyy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_yyyyyy_s_0_0_0[i] * a_z * faz_0;

        tg_zz_yyyyyz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_yyyyyz_s_0_0_0[i] * fzi_0 + tg_0_yyyyyz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_yyyyyz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_yyyyyz_s_0_0_0[i] * a_z * faz_0;

        tg_zz_yyyyzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_yyyyzz_s_0_0_0[i] * fzi_0 + tg_0_yyyyzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_yyyyzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_yyyyzz_s_0_0_0[i] * a_z * faz_0;

        tg_zz_yyyzzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_yyyzzz_s_0_0_0[i] * fzi_0 + tg_0_yyyzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_yyyzzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_yyyzzz_s_0_0_0[i] * a_z * faz_0;

        tg_zz_yyzzzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_yyzzzz_s_0_0_0[i] * fzi_0 + tg_0_yyzzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_yyzzzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_yyzzzz_s_0_0_0[i] * a_z * faz_0;

        tg_zz_yzzzzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_yzzzzz_s_0_0_0[i] * fzi_0 + tg_0_yzzzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_yzzzzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_yzzzzz_s_0_0_0[i] * a_z * faz_0;

        tg_zz_zzzzzz_s_0_0_0[i] = 1.0 / 2.0 * tg_0_zzzzzz_s_0_0_0[i] * fzi_0 + tg_0_zzzzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_zzzzzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_zzzzzz_s_0_0_0[i] * a_z * faz_0;
    }

    if (p > 0)
    {
        const double fp_0 = (double)p;

        // Set up components of auxiliary buffer : SI

        auto tg_0_xxxxxx_s_0_0_1 = pbuffer.data(idx_si_s_0_0_1);

        auto tg_0_xxxxxy_s_0_0_1 = pbuffer.data(idx_si_s_0_0_1 + 1);

        auto tg_0_xxxxxz_s_0_0_1 = pbuffer.data(idx_si_s_0_0_1 + 2);

        auto tg_0_xxxxyy_s_0_0_1 = pbuffer.data(idx_si_s_0_0_1 + 3);

        auto tg_0_xxxxyz_s_0_0_1 = pbuffer.data(idx_si_s_0_0_1 + 4);

        auto tg_0_xxxxzz_s_0_0_1 = pbuffer.data(idx_si_s_0_0_1 + 5);

        auto tg_0_xxxyyy_s_0_0_1 = pbuffer.data(idx_si_s_0_0_1 + 6);

        auto tg_0_xxxyyz_s_0_0_1 = pbuffer.data(idx_si_s_0_0_1 + 7);

        auto tg_0_xxxyzz_s_0_0_1 = pbuffer.data(idx_si_s_0_0_1 + 8);

        auto tg_0_xxxzzz_s_0_0_1 = pbuffer.data(idx_si_s_0_0_1 + 9);

        auto tg_0_xxyyyy_s_0_0_1 = pbuffer.data(idx_si_s_0_0_1 + 10);

        auto tg_0_xxyyyz_s_0_0_1 = pbuffer.data(idx_si_s_0_0_1 + 11);

        auto tg_0_xxyyzz_s_0_0_1 = pbuffer.data(idx_si_s_0_0_1 + 12);

        auto tg_0_xxyzzz_s_0_0_1 = pbuffer.data(idx_si_s_0_0_1 + 13);

        auto tg_0_xxzzzz_s_0_0_1 = pbuffer.data(idx_si_s_0_0_1 + 14);

        auto tg_0_xyyyyy_s_0_0_1 = pbuffer.data(idx_si_s_0_0_1 + 15);

        auto tg_0_xyyyyz_s_0_0_1 = pbuffer.data(idx_si_s_0_0_1 + 16);

        auto tg_0_xyyyzz_s_0_0_1 = pbuffer.data(idx_si_s_0_0_1 + 17);

        auto tg_0_xyyzzz_s_0_0_1 = pbuffer.data(idx_si_s_0_0_1 + 18);

        auto tg_0_xyzzzz_s_0_0_1 = pbuffer.data(idx_si_s_0_0_1 + 19);

        auto tg_0_xzzzzz_s_0_0_1 = pbuffer.data(idx_si_s_0_0_1 + 20);

        auto tg_0_yyyyyy_s_0_0_1 = pbuffer.data(idx_si_s_0_0_1 + 21);

        auto tg_0_yyyyyz_s_0_0_1 = pbuffer.data(idx_si_s_0_0_1 + 22);

        auto tg_0_yyyyzz_s_0_0_1 = pbuffer.data(idx_si_s_0_0_1 + 23);

        auto tg_0_yyyzzz_s_0_0_1 = pbuffer.data(idx_si_s_0_0_1 + 24);

        auto tg_0_yyzzzz_s_0_0_1 = pbuffer.data(idx_si_s_0_0_1 + 25);

        auto tg_0_yzzzzz_s_0_0_1 = pbuffer.data(idx_si_s_0_0_1 + 26);

        auto tg_0_zzzzzz_s_0_0_1 = pbuffer.data(idx_si_s_0_0_1 + 27);

        // Set up components of auxiliary buffer : PI

        auto tg_x_xxxxxx_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1);

        auto tg_x_xxxxxy_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 1);

        auto tg_x_xxxxxz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 2);

        auto tg_x_xxxxyy_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 3);

        auto tg_x_xxxxyz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 4);

        auto tg_x_xxxxzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 5);

        auto tg_x_xxxyyy_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 6);

        auto tg_x_xxxyyz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 7);

        auto tg_x_xxxyzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 8);

        auto tg_x_xxxzzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 9);

        auto tg_x_xxyyyy_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 10);

        auto tg_x_xxyyyz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 11);

        auto tg_x_xxyyzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 12);

        auto tg_x_xxyzzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 13);

        auto tg_x_xxzzzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 14);

        auto tg_x_xyyyyy_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 15);

        auto tg_x_xyyyyz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 16);

        auto tg_x_xyyyzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 17);

        auto tg_x_xyyzzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 18);

        auto tg_x_xyzzzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 19);

        auto tg_x_xzzzzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 20);

        auto tg_x_yyyyyy_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 21);

        auto tg_x_yyyyyz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 22);

        auto tg_x_yyyyzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 23);

        auto tg_x_yyyzzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 24);

        auto tg_x_yyzzzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 25);

        auto tg_x_yzzzzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 26);

        auto tg_x_zzzzzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 27);

        auto tg_y_xxxxxx_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 28);

        auto tg_y_xxxxxy_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 29);

        auto tg_y_xxxxxz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 30);

        auto tg_y_xxxxyy_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 31);

        auto tg_y_xxxxyz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 32);

        auto tg_y_xxxxzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 33);

        auto tg_y_xxxyyy_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 34);

        auto tg_y_xxxyyz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 35);

        auto tg_y_xxxyzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 36);

        auto tg_y_xxxzzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 37);

        auto tg_y_xxyyyy_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 38);

        auto tg_y_xxyyyz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 39);

        auto tg_y_xxyyzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 40);

        auto tg_y_xxyzzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 41);

        auto tg_y_xxzzzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 42);

        auto tg_y_xyyyyy_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 43);

        auto tg_y_xyyyyz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 44);

        auto tg_y_xyyyzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 45);

        auto tg_y_xyyzzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 46);

        auto tg_y_xyzzzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 47);

        auto tg_y_xzzzzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 48);

        auto tg_y_yyyyyy_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 49);

        auto tg_y_yyyyyz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 50);

        auto tg_y_yyyyzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 51);

        auto tg_y_yyyzzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 52);

        auto tg_y_yyzzzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 53);

        auto tg_y_yzzzzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 54);

        auto tg_y_zzzzzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 55);

        auto tg_z_xxxxxx_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 56);

        auto tg_z_xxxxxy_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 57);

        auto tg_z_xxxxxz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 58);

        auto tg_z_xxxxyy_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 59);

        auto tg_z_xxxxyz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 60);

        auto tg_z_xxxxzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 61);

        auto tg_z_xxxyyy_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 62);

        auto tg_z_xxxyyz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 63);

        auto tg_z_xxxyzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 64);

        auto tg_z_xxxzzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 65);

        auto tg_z_xxyyyy_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 66);

        auto tg_z_xxyyyz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 67);

        auto tg_z_xxyyzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 68);

        auto tg_z_xxyzzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 69);

        auto tg_z_xxzzzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 70);

        auto tg_z_xyyyyy_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 71);

        auto tg_z_xyyyyz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 72);

        auto tg_z_xyyyzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 73);

        auto tg_z_xyyzzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 74);

        auto tg_z_xyzzzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 75);

        auto tg_z_xzzzzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 76);

        auto tg_z_yyyyyy_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 77);

        auto tg_z_yyyyyz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 78);

        auto tg_z_yyyyzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 79);

        auto tg_z_yyyzzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 80);

        auto tg_z_yyzzzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 81);

        auto tg_z_yzzzzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 82);

        auto tg_z_zzzzzz_s_0_0_1 = pbuffer.data(idx_pi_s_0_0_1 + 83);

        #pragma omp simd aligned(b_exps, tg_0_xxxxxx_s_0_0_1, tg_0_xxxxxy_s_0_0_1, tg_0_xxxxxz_s_0_0_1, tg_0_xxxxyy_s_0_0_1, tg_0_xxxxyz_s_0_0_1, tg_0_xxxxzz_s_0_0_1, tg_0_xxxyyy_s_0_0_1, tg_0_xxxyyz_s_0_0_1, tg_0_xxxyzz_s_0_0_1, tg_0_xxxzzz_s_0_0_1, tg_0_xxyyyy_s_0_0_1, tg_0_xxyyyz_s_0_0_1, tg_0_xxyyzz_s_0_0_1, tg_0_xxyzzz_s_0_0_1, tg_0_xxzzzz_s_0_0_1, tg_0_xyyyyy_s_0_0_1, tg_0_xyyyyz_s_0_0_1, tg_0_xyyyzz_s_0_0_1, tg_0_xyyzzz_s_0_0_1, tg_0_xyzzzz_s_0_0_1, tg_0_xzzzzz_s_0_0_1, tg_0_yyyyyy_s_0_0_1, tg_0_yyyyyz_s_0_0_1, tg_0_yyyyzz_s_0_0_1, tg_0_yyyzzz_s_0_0_1, tg_0_yyzzzz_s_0_0_1, tg_0_yzzzzz_s_0_0_1, tg_0_zzzzzz_s_0_0_1, tg_x_xxxxxx_s_0_0_1, tg_x_xxxxxy_s_0_0_1, tg_x_xxxxxz_s_0_0_1, tg_x_xxxxyy_s_0_0_1, tg_x_xxxxyz_s_0_0_1, tg_x_xxxxzz_s_0_0_1, tg_x_xxxyyy_s_0_0_1, tg_x_xxxyyz_s_0_0_1, tg_x_xxxyzz_s_0_0_1, tg_x_xxxzzz_s_0_0_1, tg_x_xxyyyy_s_0_0_1, tg_x_xxyyyz_s_0_0_1, tg_x_xxyyzz_s_0_0_1, tg_x_xxyzzz_s_0_0_1, tg_x_xxzzzz_s_0_0_1, tg_x_xyyyyy_s_0_0_1, tg_x_xyyyyz_s_0_0_1, tg_x_xyyyzz_s_0_0_1, tg_x_xyyzzz_s_0_0_1, tg_x_xyzzzz_s_0_0_1, tg_x_xzzzzz_s_0_0_1, tg_x_yyyyyy_s_0_0_1, tg_x_yyyyyz_s_0_0_1, tg_x_yyyyzz_s_0_0_1, tg_x_yyyzzz_s_0_0_1, tg_x_yyzzzz_s_0_0_1, tg_x_yzzzzz_s_0_0_1, tg_x_zzzzzz_s_0_0_1, tg_xx_xxxxxx_s_0_0_0, tg_xx_xxxxxy_s_0_0_0, tg_xx_xxxxxz_s_0_0_0, tg_xx_xxxxyy_s_0_0_0, tg_xx_xxxxyz_s_0_0_0, tg_xx_xxxxzz_s_0_0_0, tg_xx_xxxyyy_s_0_0_0, tg_xx_xxxyyz_s_0_0_0, tg_xx_xxxyzz_s_0_0_0, tg_xx_xxxzzz_s_0_0_0, tg_xx_xxyyyy_s_0_0_0, tg_xx_xxyyyz_s_0_0_0, tg_xx_xxyyzz_s_0_0_0, tg_xx_xxyzzz_s_0_0_0, tg_xx_xxzzzz_s_0_0_0, tg_xx_xyyyyy_s_0_0_0, tg_xx_xyyyyz_s_0_0_0, tg_xx_xyyyzz_s_0_0_0, tg_xx_xyyzzz_s_0_0_0, tg_xx_xyzzzz_s_0_0_0, tg_xx_xzzzzz_s_0_0_0, tg_xx_yyyyyy_s_0_0_0, tg_xx_yyyyyz_s_0_0_0, tg_xx_yyyyzz_s_0_0_0, tg_xx_yyyzzz_s_0_0_0, tg_xx_yyzzzz_s_0_0_0, tg_xx_yzzzzz_s_0_0_0, tg_xx_zzzzzz_s_0_0_0, tg_xy_xxxxxx_s_0_0_0, tg_xy_xxxxxy_s_0_0_0, tg_xy_xxxxxz_s_0_0_0, tg_xy_xxxxyy_s_0_0_0, tg_xy_xxxxyz_s_0_0_0, tg_xy_xxxxzz_s_0_0_0, tg_xy_xxxyyy_s_0_0_0, tg_xy_xxxyyz_s_0_0_0, tg_xy_xxxyzz_s_0_0_0, tg_xy_xxxzzz_s_0_0_0, tg_xy_xxyyyy_s_0_0_0, tg_xy_xxyyyz_s_0_0_0, tg_xy_xxyyzz_s_0_0_0, tg_xy_xxyzzz_s_0_0_0, tg_xy_xxzzzz_s_0_0_0, tg_xy_xyyyyy_s_0_0_0, tg_xy_xyyyyz_s_0_0_0, tg_xy_xyyyzz_s_0_0_0, tg_xy_xyyzzz_s_0_0_0, tg_xy_xyzzzz_s_0_0_0, tg_xy_xzzzzz_s_0_0_0, tg_xy_yyyyyy_s_0_0_0, tg_xy_yyyyyz_s_0_0_0, tg_xy_yyyyzz_s_0_0_0, tg_xy_yyyzzz_s_0_0_0, tg_xy_yyzzzz_s_0_0_0, tg_xy_yzzzzz_s_0_0_0, tg_xy_zzzzzz_s_0_0_0, tg_xz_xxxxxx_s_0_0_0, tg_xz_xxxxxy_s_0_0_0, tg_xz_xxxxxz_s_0_0_0, tg_xz_xxxxyy_s_0_0_0, tg_xz_xxxxyz_s_0_0_0, tg_xz_xxxxzz_s_0_0_0, tg_xz_xxxyyy_s_0_0_0, tg_xz_xxxyyz_s_0_0_0, tg_xz_xxxyzz_s_0_0_0, tg_xz_xxxzzz_s_0_0_0, tg_xz_xxyyyy_s_0_0_0, tg_xz_xxyyyz_s_0_0_0, tg_xz_xxyyzz_s_0_0_0, tg_xz_xxyzzz_s_0_0_0, tg_xz_xxzzzz_s_0_0_0, tg_xz_xyyyyy_s_0_0_0, tg_xz_xyyyyz_s_0_0_0, tg_xz_xyyyzz_s_0_0_0, tg_xz_xyyzzz_s_0_0_0, tg_xz_xyzzzz_s_0_0_0, tg_xz_xzzzzz_s_0_0_0, tg_xz_yyyyyy_s_0_0_0, tg_xz_yyyyyz_s_0_0_0, tg_xz_yyyyzz_s_0_0_0, tg_xz_yyyzzz_s_0_0_0, tg_xz_yyzzzz_s_0_0_0, tg_xz_yzzzzz_s_0_0_0, tg_xz_zzzzzz_s_0_0_0, tg_y_xxxxxx_s_0_0_1, tg_y_xxxxxy_s_0_0_1, tg_y_xxxxxz_s_0_0_1, tg_y_xxxxyy_s_0_0_1, tg_y_xxxxyz_s_0_0_1, tg_y_xxxxzz_s_0_0_1, tg_y_xxxyyy_s_0_0_1, tg_y_xxxyyz_s_0_0_1, tg_y_xxxyzz_s_0_0_1, tg_y_xxxzzz_s_0_0_1, tg_y_xxyyyy_s_0_0_1, tg_y_xxyyyz_s_0_0_1, tg_y_xxyyzz_s_0_0_1, tg_y_xxyzzz_s_0_0_1, tg_y_xxzzzz_s_0_0_1, tg_y_xyyyyy_s_0_0_1, tg_y_xyyyyz_s_0_0_1, tg_y_xyyyzz_s_0_0_1, tg_y_xyyzzz_s_0_0_1, tg_y_xyzzzz_s_0_0_1, tg_y_xzzzzz_s_0_0_1, tg_y_yyyyyy_s_0_0_1, tg_y_yyyyyz_s_0_0_1, tg_y_yyyyzz_s_0_0_1, tg_y_yyyzzz_s_0_0_1, tg_y_yyzzzz_s_0_0_1, tg_y_yzzzzz_s_0_0_1, tg_y_zzzzzz_s_0_0_1, tg_yy_xxxxxx_s_0_0_0, tg_yy_xxxxxy_s_0_0_0, tg_yy_xxxxxz_s_0_0_0, tg_yy_xxxxyy_s_0_0_0, tg_yy_xxxxyz_s_0_0_0, tg_yy_xxxxzz_s_0_0_0, tg_yy_xxxyyy_s_0_0_0, tg_yy_xxxyyz_s_0_0_0, tg_yy_xxxyzz_s_0_0_0, tg_yy_xxxzzz_s_0_0_0, tg_yy_xxyyyy_s_0_0_0, tg_yy_xxyyyz_s_0_0_0, tg_yy_xxyyzz_s_0_0_0, tg_yy_xxyzzz_s_0_0_0, tg_yy_xxzzzz_s_0_0_0, tg_yy_xyyyyy_s_0_0_0, tg_yy_xyyyyz_s_0_0_0, tg_yy_xyyyzz_s_0_0_0, tg_yy_xyyzzz_s_0_0_0, tg_yy_xyzzzz_s_0_0_0, tg_yy_xzzzzz_s_0_0_0, tg_yy_yyyyyy_s_0_0_0, tg_yy_yyyyyz_s_0_0_0, tg_yy_yyyyzz_s_0_0_0, tg_yy_yyyzzz_s_0_0_0, tg_yy_yyzzzz_s_0_0_0, tg_yy_yzzzzz_s_0_0_0, tg_yy_zzzzzz_s_0_0_0, tg_yz_xxxxxx_s_0_0_0, tg_yz_xxxxxy_s_0_0_0, tg_yz_xxxxxz_s_0_0_0, tg_yz_xxxxyy_s_0_0_0, tg_yz_xxxxyz_s_0_0_0, tg_yz_xxxxzz_s_0_0_0, tg_yz_xxxyyy_s_0_0_0, tg_yz_xxxyyz_s_0_0_0, tg_yz_xxxyzz_s_0_0_0, tg_yz_xxxzzz_s_0_0_0, tg_yz_xxyyyy_s_0_0_0, tg_yz_xxyyyz_s_0_0_0, tg_yz_xxyyzz_s_0_0_0, tg_yz_xxyzzz_s_0_0_0, tg_yz_xxzzzz_s_0_0_0, tg_yz_xyyyyy_s_0_0_0, tg_yz_xyyyyz_s_0_0_0, tg_yz_xyyyzz_s_0_0_0, tg_yz_xyyzzz_s_0_0_0, tg_yz_xyzzzz_s_0_0_0, tg_yz_xzzzzz_s_0_0_0, tg_yz_yyyyyy_s_0_0_0, tg_yz_yyyyyz_s_0_0_0, tg_yz_yyyyzz_s_0_0_0, tg_yz_yyyzzz_s_0_0_0, tg_yz_yyzzzz_s_0_0_0, tg_yz_yzzzzz_s_0_0_0, tg_yz_zzzzzz_s_0_0_0, tg_z_xxxxxx_s_0_0_1, tg_z_xxxxxy_s_0_0_1, tg_z_xxxxxz_s_0_0_1, tg_z_xxxxyy_s_0_0_1, tg_z_xxxxyz_s_0_0_1, tg_z_xxxxzz_s_0_0_1, tg_z_xxxyyy_s_0_0_1, tg_z_xxxyyz_s_0_0_1, tg_z_xxxyzz_s_0_0_1, tg_z_xxxzzz_s_0_0_1, tg_z_xxyyyy_s_0_0_1, tg_z_xxyyyz_s_0_0_1, tg_z_xxyyzz_s_0_0_1, tg_z_xxyzzz_s_0_0_1, tg_z_xxzzzz_s_0_0_1, tg_z_xyyyyy_s_0_0_1, tg_z_xyyyyz_s_0_0_1, tg_z_xyyyzz_s_0_0_1, tg_z_xyyzzz_s_0_0_1, tg_z_xyzzzz_s_0_0_1, tg_z_xzzzzz_s_0_0_1, tg_z_yyyyyy_s_0_0_1, tg_z_yyyyyz_s_0_0_1, tg_z_yyyyzz_s_0_0_1, tg_z_yyyzzz_s_0_0_1, tg_z_yyzzzz_s_0_0_1, tg_z_yzzzzz_s_0_0_1, tg_z_zzzzzz_s_0_0_1, tg_zz_xxxxxx_s_0_0_0, tg_zz_xxxxxy_s_0_0_0, tg_zz_xxxxxz_s_0_0_0, tg_zz_xxxxyy_s_0_0_0, tg_zz_xxxxyz_s_0_0_0, tg_zz_xxxxzz_s_0_0_0, tg_zz_xxxyyy_s_0_0_0, tg_zz_xxxyyz_s_0_0_0, tg_zz_xxxyzz_s_0_0_0, tg_zz_xxxzzz_s_0_0_0, tg_zz_xxyyyy_s_0_0_0, tg_zz_xxyyyz_s_0_0_0, tg_zz_xxyyzz_s_0_0_0, tg_zz_xxyzzz_s_0_0_0, tg_zz_xxzzzz_s_0_0_0, tg_zz_xyyyyy_s_0_0_0, tg_zz_xyyyyz_s_0_0_0, tg_zz_xyyyzz_s_0_0_0, tg_zz_xyyzzz_s_0_0_0, tg_zz_xyzzzz_s_0_0_0, tg_zz_xzzzzz_s_0_0_0, tg_zz_yyyyyy_s_0_0_0, tg_zz_yyyyyz_s_0_0_0, tg_zz_yyyyzz_s_0_0_0, tg_zz_yyyzzz_s_0_0_0, tg_zz_yyzzzz_s_0_0_0, tg_zz_yzzzzz_s_0_0_0, tg_zz_zzzzzz_s_0_0_0  : 64)
        for (size_t i = 0; i < nelems; i++)
        {
            const double fai_0 = 1.0 / a_exp;

            tg_xx_xxxxxx_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxxxx_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xxxxxx_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xxxxxy_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxxxy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xxxxxy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xxxxxz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxxxz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xxxxxz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xxxxyy_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxxyy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xxxxyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xxxxyz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxxyz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xxxxyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xxxxzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxxzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xxxxzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xxxyyy_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxyyy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xxxyyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xxxyyz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxyyz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xxxyyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xxxyzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxyzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xxxyzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xxxzzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xxxzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xxyyyy_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxyyyy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xxyyyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xxyyyz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxyyyz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xxyyyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xxyyzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxyyzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xxyyzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xxyzzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxyzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xxyzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xxzzzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxzzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xxzzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xyyyyy_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xyyyyy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xyyyyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xyyyyz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xyyyyz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xyyyyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xyyyzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xyyyzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xyyyzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xyyzzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xyyzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xyyzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xyzzzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xyzzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xyzzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xzzzzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xzzzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xzzzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_yyyyyy_s_0_0_0[i] += 1.0 / 2.0 * tg_0_yyyyyy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_yyyyyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_yyyyyz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_yyyyyz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_yyyyyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_yyyyzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_yyyyzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_yyyyzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_yyyzzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_yyyzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_yyyzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_yyzzzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_yyzzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_yyzzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_yzzzzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_yzzzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_yzzzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_zzzzzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_zzzzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_zzzzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xxxxxx_s_0_0_0[i] += tg_y_xxxxxx_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xxxxxy_s_0_0_0[i] += tg_y_xxxxxy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xxxxxz_s_0_0_0[i] += tg_y_xxxxxz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xxxxyy_s_0_0_0[i] += tg_y_xxxxyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xxxxyz_s_0_0_0[i] += tg_y_xxxxyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xxxxzz_s_0_0_0[i] += tg_y_xxxxzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xxxyyy_s_0_0_0[i] += tg_y_xxxyyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xxxyyz_s_0_0_0[i] += tg_y_xxxyyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xxxyzz_s_0_0_0[i] += tg_y_xxxyzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xxxzzz_s_0_0_0[i] += tg_y_xxxzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xxyyyy_s_0_0_0[i] += tg_y_xxyyyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xxyyyz_s_0_0_0[i] += tg_y_xxyyyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xxyyzz_s_0_0_0[i] += tg_y_xxyyzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xxyzzz_s_0_0_0[i] += tg_y_xxyzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xxzzzz_s_0_0_0[i] += tg_y_xxzzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xyyyyy_s_0_0_0[i] += tg_y_xyyyyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xyyyyz_s_0_0_0[i] += tg_y_xyyyyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xyyyzz_s_0_0_0[i] += tg_y_xyyyzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xyyzzz_s_0_0_0[i] += tg_y_xyyzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xyzzzz_s_0_0_0[i] += tg_y_xyzzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xzzzzz_s_0_0_0[i] += tg_y_xzzzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_yyyyyy_s_0_0_0[i] += tg_y_yyyyyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_yyyyyz_s_0_0_0[i] += tg_y_yyyyyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_yyyyzz_s_0_0_0[i] += tg_y_yyyyzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_yyyzzz_s_0_0_0[i] += tg_y_yyyzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_yyzzzz_s_0_0_0[i] += tg_y_yyzzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_yzzzzz_s_0_0_0[i] += tg_y_yzzzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_zzzzzz_s_0_0_0[i] += tg_y_zzzzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xxxxxx_s_0_0_0[i] += tg_z_xxxxxx_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xxxxxy_s_0_0_0[i] += tg_z_xxxxxy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xxxxxz_s_0_0_0[i] += tg_z_xxxxxz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xxxxyy_s_0_0_0[i] += tg_z_xxxxyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xxxxyz_s_0_0_0[i] += tg_z_xxxxyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xxxxzz_s_0_0_0[i] += tg_z_xxxxzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xxxyyy_s_0_0_0[i] += tg_z_xxxyyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xxxyyz_s_0_0_0[i] += tg_z_xxxyyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xxxyzz_s_0_0_0[i] += tg_z_xxxyzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xxxzzz_s_0_0_0[i] += tg_z_xxxzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xxyyyy_s_0_0_0[i] += tg_z_xxyyyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xxyyyz_s_0_0_0[i] += tg_z_xxyyyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xxyyzz_s_0_0_0[i] += tg_z_xxyyzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xxyzzz_s_0_0_0[i] += tg_z_xxyzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xxzzzz_s_0_0_0[i] += tg_z_xxzzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xyyyyy_s_0_0_0[i] += tg_z_xyyyyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xyyyyz_s_0_0_0[i] += tg_z_xyyyyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xyyyzz_s_0_0_0[i] += tg_z_xyyyzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xyyzzz_s_0_0_0[i] += tg_z_xyyzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xyzzzz_s_0_0_0[i] += tg_z_xyzzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xzzzzz_s_0_0_0[i] += tg_z_xzzzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_yyyyyy_s_0_0_0[i] += tg_z_yyyyyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_yyyyyz_s_0_0_0[i] += tg_z_yyyyyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_yyyyzz_s_0_0_0[i] += tg_z_yyyyzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_yyyzzz_s_0_0_0[i] += tg_z_yyyzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_yyzzzz_s_0_0_0[i] += tg_z_yyzzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_yzzzzz_s_0_0_0[i] += tg_z_yzzzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_zzzzzz_s_0_0_0[i] += tg_z_zzzzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_yy_xxxxxx_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxxxx_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xxxxxx_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xxxxxy_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxxxy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xxxxxy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xxxxxz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxxxz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xxxxxz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xxxxyy_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxxyy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xxxxyy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xxxxyz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxxyz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xxxxyz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xxxxzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxxzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xxxxzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xxxyyy_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxyyy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xxxyyy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xxxyyz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxyyz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xxxyyz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xxxyzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxyzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xxxyzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xxxzzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xxxzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xxyyyy_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxyyyy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xxyyyy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xxyyyz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxyyyz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xxyyyz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xxyyzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxyyzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xxyyzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xxyzzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxyzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xxyzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xxzzzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxzzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xxzzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xyyyyy_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xyyyyy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xyyyyy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xyyyyz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xyyyyz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xyyyyz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xyyyzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xyyyzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xyyyzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xyyzzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xyyzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xyyzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xyzzzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xyzzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xyzzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xzzzzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xzzzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xzzzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_yyyyyy_s_0_0_0[i] += 1.0 / 2.0 * tg_0_yyyyyy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_yyyyyy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_yyyyyz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_yyyyyz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_yyyyyz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_yyyyzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_yyyyzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_yyyyzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_yyyzzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_yyyzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_yyyzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_yyzzzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_yyzzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_yyzzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_yzzzzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_yzzzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_yzzzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_zzzzzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_zzzzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_zzzzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xxxxxx_s_0_0_0[i] += tg_z_xxxxxx_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xxxxxy_s_0_0_0[i] += tg_z_xxxxxy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xxxxxz_s_0_0_0[i] += tg_z_xxxxxz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xxxxyy_s_0_0_0[i] += tg_z_xxxxyy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xxxxyz_s_0_0_0[i] += tg_z_xxxxyz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xxxxzz_s_0_0_0[i] += tg_z_xxxxzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xxxyyy_s_0_0_0[i] += tg_z_xxxyyy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xxxyyz_s_0_0_0[i] += tg_z_xxxyyz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xxxyzz_s_0_0_0[i] += tg_z_xxxyzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xxxzzz_s_0_0_0[i] += tg_z_xxxzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xxyyyy_s_0_0_0[i] += tg_z_xxyyyy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xxyyyz_s_0_0_0[i] += tg_z_xxyyyz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xxyyzz_s_0_0_0[i] += tg_z_xxyyzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xxyzzz_s_0_0_0[i] += tg_z_xxyzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xxzzzz_s_0_0_0[i] += tg_z_xxzzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xyyyyy_s_0_0_0[i] += tg_z_xyyyyy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xyyyyz_s_0_0_0[i] += tg_z_xyyyyz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xyyyzz_s_0_0_0[i] += tg_z_xyyyzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xyyzzz_s_0_0_0[i] += tg_z_xyyzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xyzzzz_s_0_0_0[i] += tg_z_xyzzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xzzzzz_s_0_0_0[i] += tg_z_xzzzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_yyyyyy_s_0_0_0[i] += tg_z_yyyyyy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_yyyyyz_s_0_0_0[i] += tg_z_yyyyyz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_yyyyzz_s_0_0_0[i] += tg_z_yyyyzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_yyyzzz_s_0_0_0[i] += tg_z_yyyzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_yyzzzz_s_0_0_0[i] += tg_z_yyzzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_yzzzzz_s_0_0_0[i] += tg_z_yzzzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_zzzzzz_s_0_0_0[i] += tg_z_zzzzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_zz_xxxxxx_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxxxx_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xxxxxx_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xxxxxy_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxxxy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xxxxxy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xxxxxz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxxxz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xxxxxz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xxxxyy_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxxyy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xxxxyy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xxxxyz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxxyz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xxxxyz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xxxxzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxxzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xxxxzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xxxyyy_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxyyy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xxxyyy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xxxyyz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxyyz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xxxyyz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xxxyzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxyzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xxxyzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xxxzzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xxxzzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xxyyyy_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxyyyy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xxyyyy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xxyyyz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxyyyz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xxyyyz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xxyyzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxyyzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xxyyzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xxyzzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxyzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xxyzzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xxzzzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxzzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xxzzzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xyyyyy_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xyyyyy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xyyyyy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xyyyyz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xyyyyz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xyyyyz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xyyyzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xyyyzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xyyyzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xyyzzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xyyzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xyyzzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xyzzzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xyzzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xyzzzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xzzzzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xzzzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xzzzzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_yyyyyy_s_0_0_0[i] += 1.0 / 2.0 * tg_0_yyyyyy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_yyyyyy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_yyyyyz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_yyyyyz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_yyyyyz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_yyyyzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_yyyyzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_yyyyzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_yyyzzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_yyyzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_yyyzzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_yyzzzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_yyzzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_yyzzzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_yzzzzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_yzzzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_yzzzzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_zzzzzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_zzzzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_zzzzzz_s_0_0_1[i] * fai_0 * a_z * fp_0;
        }
    }
}

} // t2pecp namespace


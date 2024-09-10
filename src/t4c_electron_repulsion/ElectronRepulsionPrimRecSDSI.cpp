#include "ElectronRepulsionPrimRecSDSI.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_sdsi(CSimdArray<double>& pbuffer,
                                  const size_t idx_eri_0_sdsi,
                                  size_t idx_eri_0_sssi,
                                  size_t idx_eri_1_sssi,
                                  size_t idx_eri_1_spsh,
                                  size_t idx_eri_0_spsi,
                                  size_t idx_eri_1_spsi,
                                  CSimdArray<double>& factors,
                                  const size_t idx_wp,
                                  const TPoint<double>& r_pb,
                                  const double a_exp,
                                  const double b_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto c_exps = factors.data(0);

    auto d_exps = factors.data(1);

    // Set up R(WP) distances

    auto wp_x = factors.data(idx_wp);

    auto wp_y = factors.data(idx_wp + 1);

    auto wp_z = factors.data(idx_wp + 2);

    // set up R(PB) distances

    const auto xyz = r_pb.coordinates();

    const auto pb_x = xyz[0];

    const auto pb_y = xyz[1];

    const auto pb_z = xyz[2];

    /// Set up components of auxilary buffer : SSSI

    auto g_0_0_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sssi);

    auto g_0_0_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sssi + 1);

    auto g_0_0_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sssi + 2);

    auto g_0_0_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sssi + 3);

    auto g_0_0_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sssi + 4);

    auto g_0_0_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sssi + 5);

    auto g_0_0_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sssi + 6);

    auto g_0_0_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sssi + 7);

    auto g_0_0_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sssi + 8);

    auto g_0_0_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sssi + 9);

    auto g_0_0_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sssi + 10);

    auto g_0_0_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sssi + 11);

    auto g_0_0_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sssi + 12);

    auto g_0_0_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sssi + 13);

    auto g_0_0_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sssi + 14);

    auto g_0_0_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sssi + 15);

    auto g_0_0_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sssi + 16);

    auto g_0_0_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sssi + 17);

    auto g_0_0_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sssi + 18);

    auto g_0_0_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sssi + 19);

    auto g_0_0_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sssi + 20);

    auto g_0_0_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sssi + 21);

    auto g_0_0_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sssi + 22);

    auto g_0_0_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sssi + 23);

    auto g_0_0_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sssi + 24);

    auto g_0_0_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sssi + 25);

    auto g_0_0_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sssi + 26);

    auto g_0_0_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sssi + 27);

    /// Set up components of auxilary buffer : SSSI

    auto g_0_0_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sssi);

    auto g_0_0_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sssi + 1);

    auto g_0_0_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sssi + 2);

    auto g_0_0_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sssi + 3);

    auto g_0_0_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sssi + 4);

    auto g_0_0_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sssi + 5);

    auto g_0_0_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sssi + 6);

    auto g_0_0_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sssi + 7);

    auto g_0_0_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sssi + 8);

    auto g_0_0_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sssi + 9);

    auto g_0_0_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sssi + 10);

    auto g_0_0_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sssi + 11);

    auto g_0_0_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sssi + 12);

    auto g_0_0_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sssi + 13);

    auto g_0_0_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sssi + 14);

    auto g_0_0_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sssi + 15);

    auto g_0_0_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sssi + 16);

    auto g_0_0_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sssi + 17);

    auto g_0_0_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sssi + 18);

    auto g_0_0_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sssi + 19);

    auto g_0_0_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sssi + 20);

    auto g_0_0_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sssi + 21);

    auto g_0_0_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sssi + 22);

    auto g_0_0_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sssi + 23);

    auto g_0_0_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sssi + 24);

    auto g_0_0_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sssi + 25);

    auto g_0_0_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sssi + 26);

    auto g_0_0_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sssi + 27);

    /// Set up components of auxilary buffer : SPSH

    auto g_0_x_0_xxxxx_1 = pbuffer.data(idx_eri_1_spsh);

    auto g_0_x_0_xxxxy_1 = pbuffer.data(idx_eri_1_spsh + 1);

    auto g_0_x_0_xxxxz_1 = pbuffer.data(idx_eri_1_spsh + 2);

    auto g_0_x_0_xxxyy_1 = pbuffer.data(idx_eri_1_spsh + 3);

    auto g_0_x_0_xxxyz_1 = pbuffer.data(idx_eri_1_spsh + 4);

    auto g_0_x_0_xxxzz_1 = pbuffer.data(idx_eri_1_spsh + 5);

    auto g_0_x_0_xxyyy_1 = pbuffer.data(idx_eri_1_spsh + 6);

    auto g_0_x_0_xxyyz_1 = pbuffer.data(idx_eri_1_spsh + 7);

    auto g_0_x_0_xxyzz_1 = pbuffer.data(idx_eri_1_spsh + 8);

    auto g_0_x_0_xxzzz_1 = pbuffer.data(idx_eri_1_spsh + 9);

    auto g_0_x_0_xyyyy_1 = pbuffer.data(idx_eri_1_spsh + 10);

    auto g_0_x_0_xyyyz_1 = pbuffer.data(idx_eri_1_spsh + 11);

    auto g_0_x_0_xyyzz_1 = pbuffer.data(idx_eri_1_spsh + 12);

    auto g_0_x_0_xyzzz_1 = pbuffer.data(idx_eri_1_spsh + 13);

    auto g_0_x_0_xzzzz_1 = pbuffer.data(idx_eri_1_spsh + 14);

    auto g_0_x_0_yyyyy_1 = pbuffer.data(idx_eri_1_spsh + 15);

    auto g_0_x_0_yyyyz_1 = pbuffer.data(idx_eri_1_spsh + 16);

    auto g_0_x_0_yyyzz_1 = pbuffer.data(idx_eri_1_spsh + 17);

    auto g_0_x_0_yyzzz_1 = pbuffer.data(idx_eri_1_spsh + 18);

    auto g_0_x_0_yzzzz_1 = pbuffer.data(idx_eri_1_spsh + 19);

    auto g_0_x_0_zzzzz_1 = pbuffer.data(idx_eri_1_spsh + 20);

    auto g_0_y_0_xxxxx_1 = pbuffer.data(idx_eri_1_spsh + 21);

    auto g_0_y_0_xxxxy_1 = pbuffer.data(idx_eri_1_spsh + 22);

    auto g_0_y_0_xxxxz_1 = pbuffer.data(idx_eri_1_spsh + 23);

    auto g_0_y_0_xxxyy_1 = pbuffer.data(idx_eri_1_spsh + 24);

    auto g_0_y_0_xxxyz_1 = pbuffer.data(idx_eri_1_spsh + 25);

    auto g_0_y_0_xxxzz_1 = pbuffer.data(idx_eri_1_spsh + 26);

    auto g_0_y_0_xxyyy_1 = pbuffer.data(idx_eri_1_spsh + 27);

    auto g_0_y_0_xxyyz_1 = pbuffer.data(idx_eri_1_spsh + 28);

    auto g_0_y_0_xxyzz_1 = pbuffer.data(idx_eri_1_spsh + 29);

    auto g_0_y_0_xxzzz_1 = pbuffer.data(idx_eri_1_spsh + 30);

    auto g_0_y_0_xyyyy_1 = pbuffer.data(idx_eri_1_spsh + 31);

    auto g_0_y_0_xyyyz_1 = pbuffer.data(idx_eri_1_spsh + 32);

    auto g_0_y_0_xyyzz_1 = pbuffer.data(idx_eri_1_spsh + 33);

    auto g_0_y_0_xyzzz_1 = pbuffer.data(idx_eri_1_spsh + 34);

    auto g_0_y_0_xzzzz_1 = pbuffer.data(idx_eri_1_spsh + 35);

    auto g_0_y_0_yyyyy_1 = pbuffer.data(idx_eri_1_spsh + 36);

    auto g_0_y_0_yyyyz_1 = pbuffer.data(idx_eri_1_spsh + 37);

    auto g_0_y_0_yyyzz_1 = pbuffer.data(idx_eri_1_spsh + 38);

    auto g_0_y_0_yyzzz_1 = pbuffer.data(idx_eri_1_spsh + 39);

    auto g_0_y_0_yzzzz_1 = pbuffer.data(idx_eri_1_spsh + 40);

    auto g_0_y_0_zzzzz_1 = pbuffer.data(idx_eri_1_spsh + 41);

    auto g_0_z_0_xxxxx_1 = pbuffer.data(idx_eri_1_spsh + 42);

    auto g_0_z_0_xxxxy_1 = pbuffer.data(idx_eri_1_spsh + 43);

    auto g_0_z_0_xxxxz_1 = pbuffer.data(idx_eri_1_spsh + 44);

    auto g_0_z_0_xxxyy_1 = pbuffer.data(idx_eri_1_spsh + 45);

    auto g_0_z_0_xxxyz_1 = pbuffer.data(idx_eri_1_spsh + 46);

    auto g_0_z_0_xxxzz_1 = pbuffer.data(idx_eri_1_spsh + 47);

    auto g_0_z_0_xxyyy_1 = pbuffer.data(idx_eri_1_spsh + 48);

    auto g_0_z_0_xxyyz_1 = pbuffer.data(idx_eri_1_spsh + 49);

    auto g_0_z_0_xxyzz_1 = pbuffer.data(idx_eri_1_spsh + 50);

    auto g_0_z_0_xxzzz_1 = pbuffer.data(idx_eri_1_spsh + 51);

    auto g_0_z_0_xyyyy_1 = pbuffer.data(idx_eri_1_spsh + 52);

    auto g_0_z_0_xyyyz_1 = pbuffer.data(idx_eri_1_spsh + 53);

    auto g_0_z_0_xyyzz_1 = pbuffer.data(idx_eri_1_spsh + 54);

    auto g_0_z_0_xyzzz_1 = pbuffer.data(idx_eri_1_spsh + 55);

    auto g_0_z_0_xzzzz_1 = pbuffer.data(idx_eri_1_spsh + 56);

    auto g_0_z_0_yyyyy_1 = pbuffer.data(idx_eri_1_spsh + 57);

    auto g_0_z_0_yyyyz_1 = pbuffer.data(idx_eri_1_spsh + 58);

    auto g_0_z_0_yyyzz_1 = pbuffer.data(idx_eri_1_spsh + 59);

    auto g_0_z_0_yyzzz_1 = pbuffer.data(idx_eri_1_spsh + 60);

    auto g_0_z_0_yzzzz_1 = pbuffer.data(idx_eri_1_spsh + 61);

    auto g_0_z_0_zzzzz_1 = pbuffer.data(idx_eri_1_spsh + 62);

    /// Set up components of auxilary buffer : SPSI

    auto g_0_x_0_xxxxxx_0 = pbuffer.data(idx_eri_0_spsi);

    auto g_0_x_0_xxxxxy_0 = pbuffer.data(idx_eri_0_spsi + 1);

    auto g_0_x_0_xxxxxz_0 = pbuffer.data(idx_eri_0_spsi + 2);

    auto g_0_x_0_xxxxyy_0 = pbuffer.data(idx_eri_0_spsi + 3);

    auto g_0_x_0_xxxxyz_0 = pbuffer.data(idx_eri_0_spsi + 4);

    auto g_0_x_0_xxxxzz_0 = pbuffer.data(idx_eri_0_spsi + 5);

    auto g_0_x_0_xxxyyy_0 = pbuffer.data(idx_eri_0_spsi + 6);

    auto g_0_x_0_xxxyyz_0 = pbuffer.data(idx_eri_0_spsi + 7);

    auto g_0_x_0_xxxyzz_0 = pbuffer.data(idx_eri_0_spsi + 8);

    auto g_0_x_0_xxxzzz_0 = pbuffer.data(idx_eri_0_spsi + 9);

    auto g_0_x_0_xxyyyy_0 = pbuffer.data(idx_eri_0_spsi + 10);

    auto g_0_x_0_xxyyyz_0 = pbuffer.data(idx_eri_0_spsi + 11);

    auto g_0_x_0_xxyyzz_0 = pbuffer.data(idx_eri_0_spsi + 12);

    auto g_0_x_0_xxyzzz_0 = pbuffer.data(idx_eri_0_spsi + 13);

    auto g_0_x_0_xxzzzz_0 = pbuffer.data(idx_eri_0_spsi + 14);

    auto g_0_x_0_xyyyyy_0 = pbuffer.data(idx_eri_0_spsi + 15);

    auto g_0_x_0_xyyyyz_0 = pbuffer.data(idx_eri_0_spsi + 16);

    auto g_0_x_0_xyyyzz_0 = pbuffer.data(idx_eri_0_spsi + 17);

    auto g_0_x_0_xyyzzz_0 = pbuffer.data(idx_eri_0_spsi + 18);

    auto g_0_x_0_xyzzzz_0 = pbuffer.data(idx_eri_0_spsi + 19);

    auto g_0_x_0_xzzzzz_0 = pbuffer.data(idx_eri_0_spsi + 20);

    auto g_0_x_0_yyyyyy_0 = pbuffer.data(idx_eri_0_spsi + 21);

    auto g_0_x_0_yyyyyz_0 = pbuffer.data(idx_eri_0_spsi + 22);

    auto g_0_x_0_yyyyzz_0 = pbuffer.data(idx_eri_0_spsi + 23);

    auto g_0_x_0_yyyzzz_0 = pbuffer.data(idx_eri_0_spsi + 24);

    auto g_0_x_0_yyzzzz_0 = pbuffer.data(idx_eri_0_spsi + 25);

    auto g_0_x_0_yzzzzz_0 = pbuffer.data(idx_eri_0_spsi + 26);

    auto g_0_x_0_zzzzzz_0 = pbuffer.data(idx_eri_0_spsi + 27);

    auto g_0_y_0_xxxxxx_0 = pbuffer.data(idx_eri_0_spsi + 28);

    auto g_0_y_0_xxxxxy_0 = pbuffer.data(idx_eri_0_spsi + 29);

    auto g_0_y_0_xxxxxz_0 = pbuffer.data(idx_eri_0_spsi + 30);

    auto g_0_y_0_xxxxyy_0 = pbuffer.data(idx_eri_0_spsi + 31);

    auto g_0_y_0_xxxxyz_0 = pbuffer.data(idx_eri_0_spsi + 32);

    auto g_0_y_0_xxxxzz_0 = pbuffer.data(idx_eri_0_spsi + 33);

    auto g_0_y_0_xxxyyy_0 = pbuffer.data(idx_eri_0_spsi + 34);

    auto g_0_y_0_xxxyyz_0 = pbuffer.data(idx_eri_0_spsi + 35);

    auto g_0_y_0_xxxyzz_0 = pbuffer.data(idx_eri_0_spsi + 36);

    auto g_0_y_0_xxxzzz_0 = pbuffer.data(idx_eri_0_spsi + 37);

    auto g_0_y_0_xxyyyy_0 = pbuffer.data(idx_eri_0_spsi + 38);

    auto g_0_y_0_xxyyyz_0 = pbuffer.data(idx_eri_0_spsi + 39);

    auto g_0_y_0_xxyyzz_0 = pbuffer.data(idx_eri_0_spsi + 40);

    auto g_0_y_0_xxyzzz_0 = pbuffer.data(idx_eri_0_spsi + 41);

    auto g_0_y_0_xxzzzz_0 = pbuffer.data(idx_eri_0_spsi + 42);

    auto g_0_y_0_xyyyyy_0 = pbuffer.data(idx_eri_0_spsi + 43);

    auto g_0_y_0_xyyyyz_0 = pbuffer.data(idx_eri_0_spsi + 44);

    auto g_0_y_0_xyyyzz_0 = pbuffer.data(idx_eri_0_spsi + 45);

    auto g_0_y_0_xyyzzz_0 = pbuffer.data(idx_eri_0_spsi + 46);

    auto g_0_y_0_xyzzzz_0 = pbuffer.data(idx_eri_0_spsi + 47);

    auto g_0_y_0_xzzzzz_0 = pbuffer.data(idx_eri_0_spsi + 48);

    auto g_0_y_0_yyyyyy_0 = pbuffer.data(idx_eri_0_spsi + 49);

    auto g_0_y_0_yyyyyz_0 = pbuffer.data(idx_eri_0_spsi + 50);

    auto g_0_y_0_yyyyzz_0 = pbuffer.data(idx_eri_0_spsi + 51);

    auto g_0_y_0_yyyzzz_0 = pbuffer.data(idx_eri_0_spsi + 52);

    auto g_0_y_0_yyzzzz_0 = pbuffer.data(idx_eri_0_spsi + 53);

    auto g_0_y_0_yzzzzz_0 = pbuffer.data(idx_eri_0_spsi + 54);

    auto g_0_y_0_zzzzzz_0 = pbuffer.data(idx_eri_0_spsi + 55);

    auto g_0_z_0_xxxxxx_0 = pbuffer.data(idx_eri_0_spsi + 56);

    auto g_0_z_0_xxxxxy_0 = pbuffer.data(idx_eri_0_spsi + 57);

    auto g_0_z_0_xxxxxz_0 = pbuffer.data(idx_eri_0_spsi + 58);

    auto g_0_z_0_xxxxyy_0 = pbuffer.data(idx_eri_0_spsi + 59);

    auto g_0_z_0_xxxxyz_0 = pbuffer.data(idx_eri_0_spsi + 60);

    auto g_0_z_0_xxxxzz_0 = pbuffer.data(idx_eri_0_spsi + 61);

    auto g_0_z_0_xxxyyy_0 = pbuffer.data(idx_eri_0_spsi + 62);

    auto g_0_z_0_xxxyyz_0 = pbuffer.data(idx_eri_0_spsi + 63);

    auto g_0_z_0_xxxyzz_0 = pbuffer.data(idx_eri_0_spsi + 64);

    auto g_0_z_0_xxxzzz_0 = pbuffer.data(idx_eri_0_spsi + 65);

    auto g_0_z_0_xxyyyy_0 = pbuffer.data(idx_eri_0_spsi + 66);

    auto g_0_z_0_xxyyyz_0 = pbuffer.data(idx_eri_0_spsi + 67);

    auto g_0_z_0_xxyyzz_0 = pbuffer.data(idx_eri_0_spsi + 68);

    auto g_0_z_0_xxyzzz_0 = pbuffer.data(idx_eri_0_spsi + 69);

    auto g_0_z_0_xxzzzz_0 = pbuffer.data(idx_eri_0_spsi + 70);

    auto g_0_z_0_xyyyyy_0 = pbuffer.data(idx_eri_0_spsi + 71);

    auto g_0_z_0_xyyyyz_0 = pbuffer.data(idx_eri_0_spsi + 72);

    auto g_0_z_0_xyyyzz_0 = pbuffer.data(idx_eri_0_spsi + 73);

    auto g_0_z_0_xyyzzz_0 = pbuffer.data(idx_eri_0_spsi + 74);

    auto g_0_z_0_xyzzzz_0 = pbuffer.data(idx_eri_0_spsi + 75);

    auto g_0_z_0_xzzzzz_0 = pbuffer.data(idx_eri_0_spsi + 76);

    auto g_0_z_0_yyyyyy_0 = pbuffer.data(idx_eri_0_spsi + 77);

    auto g_0_z_0_yyyyyz_0 = pbuffer.data(idx_eri_0_spsi + 78);

    auto g_0_z_0_yyyyzz_0 = pbuffer.data(idx_eri_0_spsi + 79);

    auto g_0_z_0_yyyzzz_0 = pbuffer.data(idx_eri_0_spsi + 80);

    auto g_0_z_0_yyzzzz_0 = pbuffer.data(idx_eri_0_spsi + 81);

    auto g_0_z_0_yzzzzz_0 = pbuffer.data(idx_eri_0_spsi + 82);

    auto g_0_z_0_zzzzzz_0 = pbuffer.data(idx_eri_0_spsi + 83);

    /// Set up components of auxilary buffer : SPSI

    auto g_0_x_0_xxxxxx_1 = pbuffer.data(idx_eri_1_spsi);

    auto g_0_x_0_xxxxxy_1 = pbuffer.data(idx_eri_1_spsi + 1);

    auto g_0_x_0_xxxxxz_1 = pbuffer.data(idx_eri_1_spsi + 2);

    auto g_0_x_0_xxxxyy_1 = pbuffer.data(idx_eri_1_spsi + 3);

    auto g_0_x_0_xxxxyz_1 = pbuffer.data(idx_eri_1_spsi + 4);

    auto g_0_x_0_xxxxzz_1 = pbuffer.data(idx_eri_1_spsi + 5);

    auto g_0_x_0_xxxyyy_1 = pbuffer.data(idx_eri_1_spsi + 6);

    auto g_0_x_0_xxxyyz_1 = pbuffer.data(idx_eri_1_spsi + 7);

    auto g_0_x_0_xxxyzz_1 = pbuffer.data(idx_eri_1_spsi + 8);

    auto g_0_x_0_xxxzzz_1 = pbuffer.data(idx_eri_1_spsi + 9);

    auto g_0_x_0_xxyyyy_1 = pbuffer.data(idx_eri_1_spsi + 10);

    auto g_0_x_0_xxyyyz_1 = pbuffer.data(idx_eri_1_spsi + 11);

    auto g_0_x_0_xxyyzz_1 = pbuffer.data(idx_eri_1_spsi + 12);

    auto g_0_x_0_xxyzzz_1 = pbuffer.data(idx_eri_1_spsi + 13);

    auto g_0_x_0_xxzzzz_1 = pbuffer.data(idx_eri_1_spsi + 14);

    auto g_0_x_0_xyyyyy_1 = pbuffer.data(idx_eri_1_spsi + 15);

    auto g_0_x_0_xyyyyz_1 = pbuffer.data(idx_eri_1_spsi + 16);

    auto g_0_x_0_xyyyzz_1 = pbuffer.data(idx_eri_1_spsi + 17);

    auto g_0_x_0_xyyzzz_1 = pbuffer.data(idx_eri_1_spsi + 18);

    auto g_0_x_0_xyzzzz_1 = pbuffer.data(idx_eri_1_spsi + 19);

    auto g_0_x_0_xzzzzz_1 = pbuffer.data(idx_eri_1_spsi + 20);

    auto g_0_x_0_yyyyyy_1 = pbuffer.data(idx_eri_1_spsi + 21);

    auto g_0_x_0_yyyyyz_1 = pbuffer.data(idx_eri_1_spsi + 22);

    auto g_0_x_0_yyyyzz_1 = pbuffer.data(idx_eri_1_spsi + 23);

    auto g_0_x_0_yyyzzz_1 = pbuffer.data(idx_eri_1_spsi + 24);

    auto g_0_x_0_yyzzzz_1 = pbuffer.data(idx_eri_1_spsi + 25);

    auto g_0_x_0_yzzzzz_1 = pbuffer.data(idx_eri_1_spsi + 26);

    auto g_0_x_0_zzzzzz_1 = pbuffer.data(idx_eri_1_spsi + 27);

    auto g_0_y_0_xxxxxx_1 = pbuffer.data(idx_eri_1_spsi + 28);

    auto g_0_y_0_xxxxxy_1 = pbuffer.data(idx_eri_1_spsi + 29);

    auto g_0_y_0_xxxxxz_1 = pbuffer.data(idx_eri_1_spsi + 30);

    auto g_0_y_0_xxxxyy_1 = pbuffer.data(idx_eri_1_spsi + 31);

    auto g_0_y_0_xxxxyz_1 = pbuffer.data(idx_eri_1_spsi + 32);

    auto g_0_y_0_xxxxzz_1 = pbuffer.data(idx_eri_1_spsi + 33);

    auto g_0_y_0_xxxyyy_1 = pbuffer.data(idx_eri_1_spsi + 34);

    auto g_0_y_0_xxxyyz_1 = pbuffer.data(idx_eri_1_spsi + 35);

    auto g_0_y_0_xxxyzz_1 = pbuffer.data(idx_eri_1_spsi + 36);

    auto g_0_y_0_xxxzzz_1 = pbuffer.data(idx_eri_1_spsi + 37);

    auto g_0_y_0_xxyyyy_1 = pbuffer.data(idx_eri_1_spsi + 38);

    auto g_0_y_0_xxyyyz_1 = pbuffer.data(idx_eri_1_spsi + 39);

    auto g_0_y_0_xxyyzz_1 = pbuffer.data(idx_eri_1_spsi + 40);

    auto g_0_y_0_xxyzzz_1 = pbuffer.data(idx_eri_1_spsi + 41);

    auto g_0_y_0_xxzzzz_1 = pbuffer.data(idx_eri_1_spsi + 42);

    auto g_0_y_0_xyyyyy_1 = pbuffer.data(idx_eri_1_spsi + 43);

    auto g_0_y_0_xyyyyz_1 = pbuffer.data(idx_eri_1_spsi + 44);

    auto g_0_y_0_xyyyzz_1 = pbuffer.data(idx_eri_1_spsi + 45);

    auto g_0_y_0_xyyzzz_1 = pbuffer.data(idx_eri_1_spsi + 46);

    auto g_0_y_0_xyzzzz_1 = pbuffer.data(idx_eri_1_spsi + 47);

    auto g_0_y_0_xzzzzz_1 = pbuffer.data(idx_eri_1_spsi + 48);

    auto g_0_y_0_yyyyyy_1 = pbuffer.data(idx_eri_1_spsi + 49);

    auto g_0_y_0_yyyyyz_1 = pbuffer.data(idx_eri_1_spsi + 50);

    auto g_0_y_0_yyyyzz_1 = pbuffer.data(idx_eri_1_spsi + 51);

    auto g_0_y_0_yyyzzz_1 = pbuffer.data(idx_eri_1_spsi + 52);

    auto g_0_y_0_yyzzzz_1 = pbuffer.data(idx_eri_1_spsi + 53);

    auto g_0_y_0_yzzzzz_1 = pbuffer.data(idx_eri_1_spsi + 54);

    auto g_0_y_0_zzzzzz_1 = pbuffer.data(idx_eri_1_spsi + 55);

    auto g_0_z_0_xxxxxx_1 = pbuffer.data(idx_eri_1_spsi + 56);

    auto g_0_z_0_xxxxxy_1 = pbuffer.data(idx_eri_1_spsi + 57);

    auto g_0_z_0_xxxxxz_1 = pbuffer.data(idx_eri_1_spsi + 58);

    auto g_0_z_0_xxxxyy_1 = pbuffer.data(idx_eri_1_spsi + 59);

    auto g_0_z_0_xxxxyz_1 = pbuffer.data(idx_eri_1_spsi + 60);

    auto g_0_z_0_xxxxzz_1 = pbuffer.data(idx_eri_1_spsi + 61);

    auto g_0_z_0_xxxyyy_1 = pbuffer.data(idx_eri_1_spsi + 62);

    auto g_0_z_0_xxxyyz_1 = pbuffer.data(idx_eri_1_spsi + 63);

    auto g_0_z_0_xxxyzz_1 = pbuffer.data(idx_eri_1_spsi + 64);

    auto g_0_z_0_xxxzzz_1 = pbuffer.data(idx_eri_1_spsi + 65);

    auto g_0_z_0_xxyyyy_1 = pbuffer.data(idx_eri_1_spsi + 66);

    auto g_0_z_0_xxyyyz_1 = pbuffer.data(idx_eri_1_spsi + 67);

    auto g_0_z_0_xxyyzz_1 = pbuffer.data(idx_eri_1_spsi + 68);

    auto g_0_z_0_xxyzzz_1 = pbuffer.data(idx_eri_1_spsi + 69);

    auto g_0_z_0_xxzzzz_1 = pbuffer.data(idx_eri_1_spsi + 70);

    auto g_0_z_0_xyyyyy_1 = pbuffer.data(idx_eri_1_spsi + 71);

    auto g_0_z_0_xyyyyz_1 = pbuffer.data(idx_eri_1_spsi + 72);

    auto g_0_z_0_xyyyzz_1 = pbuffer.data(idx_eri_1_spsi + 73);

    auto g_0_z_0_xyyzzz_1 = pbuffer.data(idx_eri_1_spsi + 74);

    auto g_0_z_0_xyzzzz_1 = pbuffer.data(idx_eri_1_spsi + 75);

    auto g_0_z_0_xzzzzz_1 = pbuffer.data(idx_eri_1_spsi + 76);

    auto g_0_z_0_yyyyyy_1 = pbuffer.data(idx_eri_1_spsi + 77);

    auto g_0_z_0_yyyyyz_1 = pbuffer.data(idx_eri_1_spsi + 78);

    auto g_0_z_0_yyyyzz_1 = pbuffer.data(idx_eri_1_spsi + 79);

    auto g_0_z_0_yyyzzz_1 = pbuffer.data(idx_eri_1_spsi + 80);

    auto g_0_z_0_yyzzzz_1 = pbuffer.data(idx_eri_1_spsi + 81);

    auto g_0_z_0_yzzzzz_1 = pbuffer.data(idx_eri_1_spsi + 82);

    auto g_0_z_0_zzzzzz_1 = pbuffer.data(idx_eri_1_spsi + 83);

    /// Set up 0-28 components of targeted buffer : SDSI

    auto g_0_xx_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sdsi);

    auto g_0_xx_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sdsi + 1);

    auto g_0_xx_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sdsi + 2);

    auto g_0_xx_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sdsi + 3);

    auto g_0_xx_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sdsi + 4);

    auto g_0_xx_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sdsi + 5);

    auto g_0_xx_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sdsi + 6);

    auto g_0_xx_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sdsi + 7);

    auto g_0_xx_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sdsi + 8);

    auto g_0_xx_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sdsi + 9);

    auto g_0_xx_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sdsi + 10);

    auto g_0_xx_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sdsi + 11);

    auto g_0_xx_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sdsi + 12);

    auto g_0_xx_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sdsi + 13);

    auto g_0_xx_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 14);

    auto g_0_xx_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sdsi + 15);

    auto g_0_xx_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sdsi + 16);

    auto g_0_xx_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sdsi + 17);

    auto g_0_xx_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sdsi + 18);

    auto g_0_xx_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 19);

    auto g_0_xx_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 20);

    auto g_0_xx_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sdsi + 21);

    auto g_0_xx_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sdsi + 22);

    auto g_0_xx_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sdsi + 23);

    auto g_0_xx_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sdsi + 24);

    auto g_0_xx_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 25);

    auto g_0_xx_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 26);

    auto g_0_xx_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 27);

    #pragma omp simd aligned(g_0_0_0_xxxxxx_0, g_0_0_0_xxxxxx_1, g_0_0_0_xxxxxy_0, g_0_0_0_xxxxxy_1, g_0_0_0_xxxxxz_0, g_0_0_0_xxxxxz_1, g_0_0_0_xxxxyy_0, g_0_0_0_xxxxyy_1, g_0_0_0_xxxxyz_0, g_0_0_0_xxxxyz_1, g_0_0_0_xxxxzz_0, g_0_0_0_xxxxzz_1, g_0_0_0_xxxyyy_0, g_0_0_0_xxxyyy_1, g_0_0_0_xxxyyz_0, g_0_0_0_xxxyyz_1, g_0_0_0_xxxyzz_0, g_0_0_0_xxxyzz_1, g_0_0_0_xxxzzz_0, g_0_0_0_xxxzzz_1, g_0_0_0_xxyyyy_0, g_0_0_0_xxyyyy_1, g_0_0_0_xxyyyz_0, g_0_0_0_xxyyyz_1, g_0_0_0_xxyyzz_0, g_0_0_0_xxyyzz_1, g_0_0_0_xxyzzz_0, g_0_0_0_xxyzzz_1, g_0_0_0_xxzzzz_0, g_0_0_0_xxzzzz_1, g_0_0_0_xyyyyy_0, g_0_0_0_xyyyyy_1, g_0_0_0_xyyyyz_0, g_0_0_0_xyyyyz_1, g_0_0_0_xyyyzz_0, g_0_0_0_xyyyzz_1, g_0_0_0_xyyzzz_0, g_0_0_0_xyyzzz_1, g_0_0_0_xyzzzz_0, g_0_0_0_xyzzzz_1, g_0_0_0_xzzzzz_0, g_0_0_0_xzzzzz_1, g_0_0_0_yyyyyy_0, g_0_0_0_yyyyyy_1, g_0_0_0_yyyyyz_0, g_0_0_0_yyyyyz_1, g_0_0_0_yyyyzz_0, g_0_0_0_yyyyzz_1, g_0_0_0_yyyzzz_0, g_0_0_0_yyyzzz_1, g_0_0_0_yyzzzz_0, g_0_0_0_yyzzzz_1, g_0_0_0_yzzzzz_0, g_0_0_0_yzzzzz_1, g_0_0_0_zzzzzz_0, g_0_0_0_zzzzzz_1, g_0_x_0_xxxxx_1, g_0_x_0_xxxxxx_0, g_0_x_0_xxxxxx_1, g_0_x_0_xxxxxy_0, g_0_x_0_xxxxxy_1, g_0_x_0_xxxxxz_0, g_0_x_0_xxxxxz_1, g_0_x_0_xxxxy_1, g_0_x_0_xxxxyy_0, g_0_x_0_xxxxyy_1, g_0_x_0_xxxxyz_0, g_0_x_0_xxxxyz_1, g_0_x_0_xxxxz_1, g_0_x_0_xxxxzz_0, g_0_x_0_xxxxzz_1, g_0_x_0_xxxyy_1, g_0_x_0_xxxyyy_0, g_0_x_0_xxxyyy_1, g_0_x_0_xxxyyz_0, g_0_x_0_xxxyyz_1, g_0_x_0_xxxyz_1, g_0_x_0_xxxyzz_0, g_0_x_0_xxxyzz_1, g_0_x_0_xxxzz_1, g_0_x_0_xxxzzz_0, g_0_x_0_xxxzzz_1, g_0_x_0_xxyyy_1, g_0_x_0_xxyyyy_0, g_0_x_0_xxyyyy_1, g_0_x_0_xxyyyz_0, g_0_x_0_xxyyyz_1, g_0_x_0_xxyyz_1, g_0_x_0_xxyyzz_0, g_0_x_0_xxyyzz_1, g_0_x_0_xxyzz_1, g_0_x_0_xxyzzz_0, g_0_x_0_xxyzzz_1, g_0_x_0_xxzzz_1, g_0_x_0_xxzzzz_0, g_0_x_0_xxzzzz_1, g_0_x_0_xyyyy_1, g_0_x_0_xyyyyy_0, g_0_x_0_xyyyyy_1, g_0_x_0_xyyyyz_0, g_0_x_0_xyyyyz_1, g_0_x_0_xyyyz_1, g_0_x_0_xyyyzz_0, g_0_x_0_xyyyzz_1, g_0_x_0_xyyzz_1, g_0_x_0_xyyzzz_0, g_0_x_0_xyyzzz_1, g_0_x_0_xyzzz_1, g_0_x_0_xyzzzz_0, g_0_x_0_xyzzzz_1, g_0_x_0_xzzzz_1, g_0_x_0_xzzzzz_0, g_0_x_0_xzzzzz_1, g_0_x_0_yyyyy_1, g_0_x_0_yyyyyy_0, g_0_x_0_yyyyyy_1, g_0_x_0_yyyyyz_0, g_0_x_0_yyyyyz_1, g_0_x_0_yyyyz_1, g_0_x_0_yyyyzz_0, g_0_x_0_yyyyzz_1, g_0_x_0_yyyzz_1, g_0_x_0_yyyzzz_0, g_0_x_0_yyyzzz_1, g_0_x_0_yyzzz_1, g_0_x_0_yyzzzz_0, g_0_x_0_yyzzzz_1, g_0_x_0_yzzzz_1, g_0_x_0_yzzzzz_0, g_0_x_0_yzzzzz_1, g_0_x_0_zzzzz_1, g_0_x_0_zzzzzz_0, g_0_x_0_zzzzzz_1, g_0_xx_0_xxxxxx_0, g_0_xx_0_xxxxxy_0, g_0_xx_0_xxxxxz_0, g_0_xx_0_xxxxyy_0, g_0_xx_0_xxxxyz_0, g_0_xx_0_xxxxzz_0, g_0_xx_0_xxxyyy_0, g_0_xx_0_xxxyyz_0, g_0_xx_0_xxxyzz_0, g_0_xx_0_xxxzzz_0, g_0_xx_0_xxyyyy_0, g_0_xx_0_xxyyyz_0, g_0_xx_0_xxyyzz_0, g_0_xx_0_xxyzzz_0, g_0_xx_0_xxzzzz_0, g_0_xx_0_xyyyyy_0, g_0_xx_0_xyyyyz_0, g_0_xx_0_xyyyzz_0, g_0_xx_0_xyyzzz_0, g_0_xx_0_xyzzzz_0, g_0_xx_0_xzzzzz_0, g_0_xx_0_yyyyyy_0, g_0_xx_0_yyyyyz_0, g_0_xx_0_yyyyzz_0, g_0_xx_0_yyyzzz_0, g_0_xx_0_yyzzzz_0, g_0_xx_0_yzzzzz_0, g_0_xx_0_zzzzzz_0, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xx_0_xxxxxx_0[i] = g_0_0_0_xxxxxx_0[i] * fi_ab_0 - g_0_0_0_xxxxxx_1[i] * fti_ab_0 + 6.0 * g_0_x_0_xxxxx_1[i] * fi_abcd_0 + g_0_x_0_xxxxxx_0[i] * pb_x + g_0_x_0_xxxxxx_1[i] * wp_x[i];

        g_0_xx_0_xxxxxy_0[i] = g_0_0_0_xxxxxy_0[i] * fi_ab_0 - g_0_0_0_xxxxxy_1[i] * fti_ab_0 + 5.0 * g_0_x_0_xxxxy_1[i] * fi_abcd_0 + g_0_x_0_xxxxxy_0[i] * pb_x + g_0_x_0_xxxxxy_1[i] * wp_x[i];

        g_0_xx_0_xxxxxz_0[i] = g_0_0_0_xxxxxz_0[i] * fi_ab_0 - g_0_0_0_xxxxxz_1[i] * fti_ab_0 + 5.0 * g_0_x_0_xxxxz_1[i] * fi_abcd_0 + g_0_x_0_xxxxxz_0[i] * pb_x + g_0_x_0_xxxxxz_1[i] * wp_x[i];

        g_0_xx_0_xxxxyy_0[i] = g_0_0_0_xxxxyy_0[i] * fi_ab_0 - g_0_0_0_xxxxyy_1[i] * fti_ab_0 + 4.0 * g_0_x_0_xxxyy_1[i] * fi_abcd_0 + g_0_x_0_xxxxyy_0[i] * pb_x + g_0_x_0_xxxxyy_1[i] * wp_x[i];

        g_0_xx_0_xxxxyz_0[i] = g_0_0_0_xxxxyz_0[i] * fi_ab_0 - g_0_0_0_xxxxyz_1[i] * fti_ab_0 + 4.0 * g_0_x_0_xxxyz_1[i] * fi_abcd_0 + g_0_x_0_xxxxyz_0[i] * pb_x + g_0_x_0_xxxxyz_1[i] * wp_x[i];

        g_0_xx_0_xxxxzz_0[i] = g_0_0_0_xxxxzz_0[i] * fi_ab_0 - g_0_0_0_xxxxzz_1[i] * fti_ab_0 + 4.0 * g_0_x_0_xxxzz_1[i] * fi_abcd_0 + g_0_x_0_xxxxzz_0[i] * pb_x + g_0_x_0_xxxxzz_1[i] * wp_x[i];

        g_0_xx_0_xxxyyy_0[i] = g_0_0_0_xxxyyy_0[i] * fi_ab_0 - g_0_0_0_xxxyyy_1[i] * fti_ab_0 + 3.0 * g_0_x_0_xxyyy_1[i] * fi_abcd_0 + g_0_x_0_xxxyyy_0[i] * pb_x + g_0_x_0_xxxyyy_1[i] * wp_x[i];

        g_0_xx_0_xxxyyz_0[i] = g_0_0_0_xxxyyz_0[i] * fi_ab_0 - g_0_0_0_xxxyyz_1[i] * fti_ab_0 + 3.0 * g_0_x_0_xxyyz_1[i] * fi_abcd_0 + g_0_x_0_xxxyyz_0[i] * pb_x + g_0_x_0_xxxyyz_1[i] * wp_x[i];

        g_0_xx_0_xxxyzz_0[i] = g_0_0_0_xxxyzz_0[i] * fi_ab_0 - g_0_0_0_xxxyzz_1[i] * fti_ab_0 + 3.0 * g_0_x_0_xxyzz_1[i] * fi_abcd_0 + g_0_x_0_xxxyzz_0[i] * pb_x + g_0_x_0_xxxyzz_1[i] * wp_x[i];

        g_0_xx_0_xxxzzz_0[i] = g_0_0_0_xxxzzz_0[i] * fi_ab_0 - g_0_0_0_xxxzzz_1[i] * fti_ab_0 + 3.0 * g_0_x_0_xxzzz_1[i] * fi_abcd_0 + g_0_x_0_xxxzzz_0[i] * pb_x + g_0_x_0_xxxzzz_1[i] * wp_x[i];

        g_0_xx_0_xxyyyy_0[i] = g_0_0_0_xxyyyy_0[i] * fi_ab_0 - g_0_0_0_xxyyyy_1[i] * fti_ab_0 + 2.0 * g_0_x_0_xyyyy_1[i] * fi_abcd_0 + g_0_x_0_xxyyyy_0[i] * pb_x + g_0_x_0_xxyyyy_1[i] * wp_x[i];

        g_0_xx_0_xxyyyz_0[i] = g_0_0_0_xxyyyz_0[i] * fi_ab_0 - g_0_0_0_xxyyyz_1[i] * fti_ab_0 + 2.0 * g_0_x_0_xyyyz_1[i] * fi_abcd_0 + g_0_x_0_xxyyyz_0[i] * pb_x + g_0_x_0_xxyyyz_1[i] * wp_x[i];

        g_0_xx_0_xxyyzz_0[i] = g_0_0_0_xxyyzz_0[i] * fi_ab_0 - g_0_0_0_xxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_x_0_xyyzz_1[i] * fi_abcd_0 + g_0_x_0_xxyyzz_0[i] * pb_x + g_0_x_0_xxyyzz_1[i] * wp_x[i];

        g_0_xx_0_xxyzzz_0[i] = g_0_0_0_xxyzzz_0[i] * fi_ab_0 - g_0_0_0_xxyzzz_1[i] * fti_ab_0 + 2.0 * g_0_x_0_xyzzz_1[i] * fi_abcd_0 + g_0_x_0_xxyzzz_0[i] * pb_x + g_0_x_0_xxyzzz_1[i] * wp_x[i];

        g_0_xx_0_xxzzzz_0[i] = g_0_0_0_xxzzzz_0[i] * fi_ab_0 - g_0_0_0_xxzzzz_1[i] * fti_ab_0 + 2.0 * g_0_x_0_xzzzz_1[i] * fi_abcd_0 + g_0_x_0_xxzzzz_0[i] * pb_x + g_0_x_0_xxzzzz_1[i] * wp_x[i];

        g_0_xx_0_xyyyyy_0[i] = g_0_0_0_xyyyyy_0[i] * fi_ab_0 - g_0_0_0_xyyyyy_1[i] * fti_ab_0 + g_0_x_0_yyyyy_1[i] * fi_abcd_0 + g_0_x_0_xyyyyy_0[i] * pb_x + g_0_x_0_xyyyyy_1[i] * wp_x[i];

        g_0_xx_0_xyyyyz_0[i] = g_0_0_0_xyyyyz_0[i] * fi_ab_0 - g_0_0_0_xyyyyz_1[i] * fti_ab_0 + g_0_x_0_yyyyz_1[i] * fi_abcd_0 + g_0_x_0_xyyyyz_0[i] * pb_x + g_0_x_0_xyyyyz_1[i] * wp_x[i];

        g_0_xx_0_xyyyzz_0[i] = g_0_0_0_xyyyzz_0[i] * fi_ab_0 - g_0_0_0_xyyyzz_1[i] * fti_ab_0 + g_0_x_0_yyyzz_1[i] * fi_abcd_0 + g_0_x_0_xyyyzz_0[i] * pb_x + g_0_x_0_xyyyzz_1[i] * wp_x[i];

        g_0_xx_0_xyyzzz_0[i] = g_0_0_0_xyyzzz_0[i] * fi_ab_0 - g_0_0_0_xyyzzz_1[i] * fti_ab_0 + g_0_x_0_yyzzz_1[i] * fi_abcd_0 + g_0_x_0_xyyzzz_0[i] * pb_x + g_0_x_0_xyyzzz_1[i] * wp_x[i];

        g_0_xx_0_xyzzzz_0[i] = g_0_0_0_xyzzzz_0[i] * fi_ab_0 - g_0_0_0_xyzzzz_1[i] * fti_ab_0 + g_0_x_0_yzzzz_1[i] * fi_abcd_0 + g_0_x_0_xyzzzz_0[i] * pb_x + g_0_x_0_xyzzzz_1[i] * wp_x[i];

        g_0_xx_0_xzzzzz_0[i] = g_0_0_0_xzzzzz_0[i] * fi_ab_0 - g_0_0_0_xzzzzz_1[i] * fti_ab_0 + g_0_x_0_zzzzz_1[i] * fi_abcd_0 + g_0_x_0_xzzzzz_0[i] * pb_x + g_0_x_0_xzzzzz_1[i] * wp_x[i];

        g_0_xx_0_yyyyyy_0[i] = g_0_0_0_yyyyyy_0[i] * fi_ab_0 - g_0_0_0_yyyyyy_1[i] * fti_ab_0 + g_0_x_0_yyyyyy_0[i] * pb_x + g_0_x_0_yyyyyy_1[i] * wp_x[i];

        g_0_xx_0_yyyyyz_0[i] = g_0_0_0_yyyyyz_0[i] * fi_ab_0 - g_0_0_0_yyyyyz_1[i] * fti_ab_0 + g_0_x_0_yyyyyz_0[i] * pb_x + g_0_x_0_yyyyyz_1[i] * wp_x[i];

        g_0_xx_0_yyyyzz_0[i] = g_0_0_0_yyyyzz_0[i] * fi_ab_0 - g_0_0_0_yyyyzz_1[i] * fti_ab_0 + g_0_x_0_yyyyzz_0[i] * pb_x + g_0_x_0_yyyyzz_1[i] * wp_x[i];

        g_0_xx_0_yyyzzz_0[i] = g_0_0_0_yyyzzz_0[i] * fi_ab_0 - g_0_0_0_yyyzzz_1[i] * fti_ab_0 + g_0_x_0_yyyzzz_0[i] * pb_x + g_0_x_0_yyyzzz_1[i] * wp_x[i];

        g_0_xx_0_yyzzzz_0[i] = g_0_0_0_yyzzzz_0[i] * fi_ab_0 - g_0_0_0_yyzzzz_1[i] * fti_ab_0 + g_0_x_0_yyzzzz_0[i] * pb_x + g_0_x_0_yyzzzz_1[i] * wp_x[i];

        g_0_xx_0_yzzzzz_0[i] = g_0_0_0_yzzzzz_0[i] * fi_ab_0 - g_0_0_0_yzzzzz_1[i] * fti_ab_0 + g_0_x_0_yzzzzz_0[i] * pb_x + g_0_x_0_yzzzzz_1[i] * wp_x[i];

        g_0_xx_0_zzzzzz_0[i] = g_0_0_0_zzzzzz_0[i] * fi_ab_0 - g_0_0_0_zzzzzz_1[i] * fti_ab_0 + g_0_x_0_zzzzzz_0[i] * pb_x + g_0_x_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 28-56 components of targeted buffer : SDSI

    auto g_0_xy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sdsi + 28);

    auto g_0_xy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sdsi + 29);

    auto g_0_xy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sdsi + 30);

    auto g_0_xy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sdsi + 31);

    auto g_0_xy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sdsi + 32);

    auto g_0_xy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sdsi + 33);

    auto g_0_xy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sdsi + 34);

    auto g_0_xy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sdsi + 35);

    auto g_0_xy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sdsi + 36);

    auto g_0_xy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sdsi + 37);

    auto g_0_xy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sdsi + 38);

    auto g_0_xy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sdsi + 39);

    auto g_0_xy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sdsi + 40);

    auto g_0_xy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sdsi + 41);

    auto g_0_xy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 42);

    auto g_0_xy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sdsi + 43);

    auto g_0_xy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sdsi + 44);

    auto g_0_xy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sdsi + 45);

    auto g_0_xy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sdsi + 46);

    auto g_0_xy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 47);

    auto g_0_xy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 48);

    auto g_0_xy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sdsi + 49);

    auto g_0_xy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sdsi + 50);

    auto g_0_xy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sdsi + 51);

    auto g_0_xy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sdsi + 52);

    auto g_0_xy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 53);

    auto g_0_xy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 54);

    auto g_0_xy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 55);

    #pragma omp simd aligned(g_0_x_0_xxxxxx_0, g_0_x_0_xxxxxx_1, g_0_x_0_xxxxxz_0, g_0_x_0_xxxxxz_1, g_0_x_0_xxxxzz_0, g_0_x_0_xxxxzz_1, g_0_x_0_xxxzzz_0, g_0_x_0_xxxzzz_1, g_0_x_0_xxzzzz_0, g_0_x_0_xxzzzz_1, g_0_x_0_xzzzzz_0, g_0_x_0_xzzzzz_1, g_0_xy_0_xxxxxx_0, g_0_xy_0_xxxxxy_0, g_0_xy_0_xxxxxz_0, g_0_xy_0_xxxxyy_0, g_0_xy_0_xxxxyz_0, g_0_xy_0_xxxxzz_0, g_0_xy_0_xxxyyy_0, g_0_xy_0_xxxyyz_0, g_0_xy_0_xxxyzz_0, g_0_xy_0_xxxzzz_0, g_0_xy_0_xxyyyy_0, g_0_xy_0_xxyyyz_0, g_0_xy_0_xxyyzz_0, g_0_xy_0_xxyzzz_0, g_0_xy_0_xxzzzz_0, g_0_xy_0_xyyyyy_0, g_0_xy_0_xyyyyz_0, g_0_xy_0_xyyyzz_0, g_0_xy_0_xyyzzz_0, g_0_xy_0_xyzzzz_0, g_0_xy_0_xzzzzz_0, g_0_xy_0_yyyyyy_0, g_0_xy_0_yyyyyz_0, g_0_xy_0_yyyyzz_0, g_0_xy_0_yyyzzz_0, g_0_xy_0_yyzzzz_0, g_0_xy_0_yzzzzz_0, g_0_xy_0_zzzzzz_0, g_0_y_0_xxxxxy_0, g_0_y_0_xxxxxy_1, g_0_y_0_xxxxy_1, g_0_y_0_xxxxyy_0, g_0_y_0_xxxxyy_1, g_0_y_0_xxxxyz_0, g_0_y_0_xxxxyz_1, g_0_y_0_xxxyy_1, g_0_y_0_xxxyyy_0, g_0_y_0_xxxyyy_1, g_0_y_0_xxxyyz_0, g_0_y_0_xxxyyz_1, g_0_y_0_xxxyz_1, g_0_y_0_xxxyzz_0, g_0_y_0_xxxyzz_1, g_0_y_0_xxyyy_1, g_0_y_0_xxyyyy_0, g_0_y_0_xxyyyy_1, g_0_y_0_xxyyyz_0, g_0_y_0_xxyyyz_1, g_0_y_0_xxyyz_1, g_0_y_0_xxyyzz_0, g_0_y_0_xxyyzz_1, g_0_y_0_xxyzz_1, g_0_y_0_xxyzzz_0, g_0_y_0_xxyzzz_1, g_0_y_0_xyyyy_1, g_0_y_0_xyyyyy_0, g_0_y_0_xyyyyy_1, g_0_y_0_xyyyyz_0, g_0_y_0_xyyyyz_1, g_0_y_0_xyyyz_1, g_0_y_0_xyyyzz_0, g_0_y_0_xyyyzz_1, g_0_y_0_xyyzz_1, g_0_y_0_xyyzzz_0, g_0_y_0_xyyzzz_1, g_0_y_0_xyzzz_1, g_0_y_0_xyzzzz_0, g_0_y_0_xyzzzz_1, g_0_y_0_yyyyy_1, g_0_y_0_yyyyyy_0, g_0_y_0_yyyyyy_1, g_0_y_0_yyyyyz_0, g_0_y_0_yyyyyz_1, g_0_y_0_yyyyz_1, g_0_y_0_yyyyzz_0, g_0_y_0_yyyyzz_1, g_0_y_0_yyyzz_1, g_0_y_0_yyyzzz_0, g_0_y_0_yyyzzz_1, g_0_y_0_yyzzz_1, g_0_y_0_yyzzzz_0, g_0_y_0_yyzzzz_1, g_0_y_0_yzzzz_1, g_0_y_0_yzzzzz_0, g_0_y_0_yzzzzz_1, g_0_y_0_zzzzzz_0, g_0_y_0_zzzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xy_0_xxxxxx_0[i] = g_0_x_0_xxxxxx_0[i] * pb_y + g_0_x_0_xxxxxx_1[i] * wp_y[i];

        g_0_xy_0_xxxxxy_0[i] = 5.0 * g_0_y_0_xxxxy_1[i] * fi_abcd_0 + g_0_y_0_xxxxxy_0[i] * pb_x + g_0_y_0_xxxxxy_1[i] * wp_x[i];

        g_0_xy_0_xxxxxz_0[i] = g_0_x_0_xxxxxz_0[i] * pb_y + g_0_x_0_xxxxxz_1[i] * wp_y[i];

        g_0_xy_0_xxxxyy_0[i] = 4.0 * g_0_y_0_xxxyy_1[i] * fi_abcd_0 + g_0_y_0_xxxxyy_0[i] * pb_x + g_0_y_0_xxxxyy_1[i] * wp_x[i];

        g_0_xy_0_xxxxyz_0[i] = 4.0 * g_0_y_0_xxxyz_1[i] * fi_abcd_0 + g_0_y_0_xxxxyz_0[i] * pb_x + g_0_y_0_xxxxyz_1[i] * wp_x[i];

        g_0_xy_0_xxxxzz_0[i] = g_0_x_0_xxxxzz_0[i] * pb_y + g_0_x_0_xxxxzz_1[i] * wp_y[i];

        g_0_xy_0_xxxyyy_0[i] = 3.0 * g_0_y_0_xxyyy_1[i] * fi_abcd_0 + g_0_y_0_xxxyyy_0[i] * pb_x + g_0_y_0_xxxyyy_1[i] * wp_x[i];

        g_0_xy_0_xxxyyz_0[i] = 3.0 * g_0_y_0_xxyyz_1[i] * fi_abcd_0 + g_0_y_0_xxxyyz_0[i] * pb_x + g_0_y_0_xxxyyz_1[i] * wp_x[i];

        g_0_xy_0_xxxyzz_0[i] = 3.0 * g_0_y_0_xxyzz_1[i] * fi_abcd_0 + g_0_y_0_xxxyzz_0[i] * pb_x + g_0_y_0_xxxyzz_1[i] * wp_x[i];

        g_0_xy_0_xxxzzz_0[i] = g_0_x_0_xxxzzz_0[i] * pb_y + g_0_x_0_xxxzzz_1[i] * wp_y[i];

        g_0_xy_0_xxyyyy_0[i] = 2.0 * g_0_y_0_xyyyy_1[i] * fi_abcd_0 + g_0_y_0_xxyyyy_0[i] * pb_x + g_0_y_0_xxyyyy_1[i] * wp_x[i];

        g_0_xy_0_xxyyyz_0[i] = 2.0 * g_0_y_0_xyyyz_1[i] * fi_abcd_0 + g_0_y_0_xxyyyz_0[i] * pb_x + g_0_y_0_xxyyyz_1[i] * wp_x[i];

        g_0_xy_0_xxyyzz_0[i] = 2.0 * g_0_y_0_xyyzz_1[i] * fi_abcd_0 + g_0_y_0_xxyyzz_0[i] * pb_x + g_0_y_0_xxyyzz_1[i] * wp_x[i];

        g_0_xy_0_xxyzzz_0[i] = 2.0 * g_0_y_0_xyzzz_1[i] * fi_abcd_0 + g_0_y_0_xxyzzz_0[i] * pb_x + g_0_y_0_xxyzzz_1[i] * wp_x[i];

        g_0_xy_0_xxzzzz_0[i] = g_0_x_0_xxzzzz_0[i] * pb_y + g_0_x_0_xxzzzz_1[i] * wp_y[i];

        g_0_xy_0_xyyyyy_0[i] = g_0_y_0_yyyyy_1[i] * fi_abcd_0 + g_0_y_0_xyyyyy_0[i] * pb_x + g_0_y_0_xyyyyy_1[i] * wp_x[i];

        g_0_xy_0_xyyyyz_0[i] = g_0_y_0_yyyyz_1[i] * fi_abcd_0 + g_0_y_0_xyyyyz_0[i] * pb_x + g_0_y_0_xyyyyz_1[i] * wp_x[i];

        g_0_xy_0_xyyyzz_0[i] = g_0_y_0_yyyzz_1[i] * fi_abcd_0 + g_0_y_0_xyyyzz_0[i] * pb_x + g_0_y_0_xyyyzz_1[i] * wp_x[i];

        g_0_xy_0_xyyzzz_0[i] = g_0_y_0_yyzzz_1[i] * fi_abcd_0 + g_0_y_0_xyyzzz_0[i] * pb_x + g_0_y_0_xyyzzz_1[i] * wp_x[i];

        g_0_xy_0_xyzzzz_0[i] = g_0_y_0_yzzzz_1[i] * fi_abcd_0 + g_0_y_0_xyzzzz_0[i] * pb_x + g_0_y_0_xyzzzz_1[i] * wp_x[i];

        g_0_xy_0_xzzzzz_0[i] = g_0_x_0_xzzzzz_0[i] * pb_y + g_0_x_0_xzzzzz_1[i] * wp_y[i];

        g_0_xy_0_yyyyyy_0[i] = g_0_y_0_yyyyyy_0[i] * pb_x + g_0_y_0_yyyyyy_1[i] * wp_x[i];

        g_0_xy_0_yyyyyz_0[i] = g_0_y_0_yyyyyz_0[i] * pb_x + g_0_y_0_yyyyyz_1[i] * wp_x[i];

        g_0_xy_0_yyyyzz_0[i] = g_0_y_0_yyyyzz_0[i] * pb_x + g_0_y_0_yyyyzz_1[i] * wp_x[i];

        g_0_xy_0_yyyzzz_0[i] = g_0_y_0_yyyzzz_0[i] * pb_x + g_0_y_0_yyyzzz_1[i] * wp_x[i];

        g_0_xy_0_yyzzzz_0[i] = g_0_y_0_yyzzzz_0[i] * pb_x + g_0_y_0_yyzzzz_1[i] * wp_x[i];

        g_0_xy_0_yzzzzz_0[i] = g_0_y_0_yzzzzz_0[i] * pb_x + g_0_y_0_yzzzzz_1[i] * wp_x[i];

        g_0_xy_0_zzzzzz_0[i] = g_0_y_0_zzzzzz_0[i] * pb_x + g_0_y_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 56-84 components of targeted buffer : SDSI

    auto g_0_xz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sdsi + 56);

    auto g_0_xz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sdsi + 57);

    auto g_0_xz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sdsi + 58);

    auto g_0_xz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sdsi + 59);

    auto g_0_xz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sdsi + 60);

    auto g_0_xz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sdsi + 61);

    auto g_0_xz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sdsi + 62);

    auto g_0_xz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sdsi + 63);

    auto g_0_xz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sdsi + 64);

    auto g_0_xz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sdsi + 65);

    auto g_0_xz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sdsi + 66);

    auto g_0_xz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sdsi + 67);

    auto g_0_xz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sdsi + 68);

    auto g_0_xz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sdsi + 69);

    auto g_0_xz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 70);

    auto g_0_xz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sdsi + 71);

    auto g_0_xz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sdsi + 72);

    auto g_0_xz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sdsi + 73);

    auto g_0_xz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sdsi + 74);

    auto g_0_xz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 75);

    auto g_0_xz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 76);

    auto g_0_xz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sdsi + 77);

    auto g_0_xz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sdsi + 78);

    auto g_0_xz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sdsi + 79);

    auto g_0_xz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sdsi + 80);

    auto g_0_xz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 81);

    auto g_0_xz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 82);

    auto g_0_xz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 83);

    #pragma omp simd aligned(g_0_x_0_xxxxxx_0, g_0_x_0_xxxxxx_1, g_0_x_0_xxxxxy_0, g_0_x_0_xxxxxy_1, g_0_x_0_xxxxyy_0, g_0_x_0_xxxxyy_1, g_0_x_0_xxxyyy_0, g_0_x_0_xxxyyy_1, g_0_x_0_xxyyyy_0, g_0_x_0_xxyyyy_1, g_0_x_0_xyyyyy_0, g_0_x_0_xyyyyy_1, g_0_xz_0_xxxxxx_0, g_0_xz_0_xxxxxy_0, g_0_xz_0_xxxxxz_0, g_0_xz_0_xxxxyy_0, g_0_xz_0_xxxxyz_0, g_0_xz_0_xxxxzz_0, g_0_xz_0_xxxyyy_0, g_0_xz_0_xxxyyz_0, g_0_xz_0_xxxyzz_0, g_0_xz_0_xxxzzz_0, g_0_xz_0_xxyyyy_0, g_0_xz_0_xxyyyz_0, g_0_xz_0_xxyyzz_0, g_0_xz_0_xxyzzz_0, g_0_xz_0_xxzzzz_0, g_0_xz_0_xyyyyy_0, g_0_xz_0_xyyyyz_0, g_0_xz_0_xyyyzz_0, g_0_xz_0_xyyzzz_0, g_0_xz_0_xyzzzz_0, g_0_xz_0_xzzzzz_0, g_0_xz_0_yyyyyy_0, g_0_xz_0_yyyyyz_0, g_0_xz_0_yyyyzz_0, g_0_xz_0_yyyzzz_0, g_0_xz_0_yyzzzz_0, g_0_xz_0_yzzzzz_0, g_0_xz_0_zzzzzz_0, g_0_z_0_xxxxxz_0, g_0_z_0_xxxxxz_1, g_0_z_0_xxxxyz_0, g_0_z_0_xxxxyz_1, g_0_z_0_xxxxz_1, g_0_z_0_xxxxzz_0, g_0_z_0_xxxxzz_1, g_0_z_0_xxxyyz_0, g_0_z_0_xxxyyz_1, g_0_z_0_xxxyz_1, g_0_z_0_xxxyzz_0, g_0_z_0_xxxyzz_1, g_0_z_0_xxxzz_1, g_0_z_0_xxxzzz_0, g_0_z_0_xxxzzz_1, g_0_z_0_xxyyyz_0, g_0_z_0_xxyyyz_1, g_0_z_0_xxyyz_1, g_0_z_0_xxyyzz_0, g_0_z_0_xxyyzz_1, g_0_z_0_xxyzz_1, g_0_z_0_xxyzzz_0, g_0_z_0_xxyzzz_1, g_0_z_0_xxzzz_1, g_0_z_0_xxzzzz_0, g_0_z_0_xxzzzz_1, g_0_z_0_xyyyyz_0, g_0_z_0_xyyyyz_1, g_0_z_0_xyyyz_1, g_0_z_0_xyyyzz_0, g_0_z_0_xyyyzz_1, g_0_z_0_xyyzz_1, g_0_z_0_xyyzzz_0, g_0_z_0_xyyzzz_1, g_0_z_0_xyzzz_1, g_0_z_0_xyzzzz_0, g_0_z_0_xyzzzz_1, g_0_z_0_xzzzz_1, g_0_z_0_xzzzzz_0, g_0_z_0_xzzzzz_1, g_0_z_0_yyyyyy_0, g_0_z_0_yyyyyy_1, g_0_z_0_yyyyyz_0, g_0_z_0_yyyyyz_1, g_0_z_0_yyyyz_1, g_0_z_0_yyyyzz_0, g_0_z_0_yyyyzz_1, g_0_z_0_yyyzz_1, g_0_z_0_yyyzzz_0, g_0_z_0_yyyzzz_1, g_0_z_0_yyzzz_1, g_0_z_0_yyzzzz_0, g_0_z_0_yyzzzz_1, g_0_z_0_yzzzz_1, g_0_z_0_yzzzzz_0, g_0_z_0_yzzzzz_1, g_0_z_0_zzzzz_1, g_0_z_0_zzzzzz_0, g_0_z_0_zzzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xz_0_xxxxxx_0[i] = g_0_x_0_xxxxxx_0[i] * pb_z + g_0_x_0_xxxxxx_1[i] * wp_z[i];

        g_0_xz_0_xxxxxy_0[i] = g_0_x_0_xxxxxy_0[i] * pb_z + g_0_x_0_xxxxxy_1[i] * wp_z[i];

        g_0_xz_0_xxxxxz_0[i] = 5.0 * g_0_z_0_xxxxz_1[i] * fi_abcd_0 + g_0_z_0_xxxxxz_0[i] * pb_x + g_0_z_0_xxxxxz_1[i] * wp_x[i];

        g_0_xz_0_xxxxyy_0[i] = g_0_x_0_xxxxyy_0[i] * pb_z + g_0_x_0_xxxxyy_1[i] * wp_z[i];

        g_0_xz_0_xxxxyz_0[i] = 4.0 * g_0_z_0_xxxyz_1[i] * fi_abcd_0 + g_0_z_0_xxxxyz_0[i] * pb_x + g_0_z_0_xxxxyz_1[i] * wp_x[i];

        g_0_xz_0_xxxxzz_0[i] = 4.0 * g_0_z_0_xxxzz_1[i] * fi_abcd_0 + g_0_z_0_xxxxzz_0[i] * pb_x + g_0_z_0_xxxxzz_1[i] * wp_x[i];

        g_0_xz_0_xxxyyy_0[i] = g_0_x_0_xxxyyy_0[i] * pb_z + g_0_x_0_xxxyyy_1[i] * wp_z[i];

        g_0_xz_0_xxxyyz_0[i] = 3.0 * g_0_z_0_xxyyz_1[i] * fi_abcd_0 + g_0_z_0_xxxyyz_0[i] * pb_x + g_0_z_0_xxxyyz_1[i] * wp_x[i];

        g_0_xz_0_xxxyzz_0[i] = 3.0 * g_0_z_0_xxyzz_1[i] * fi_abcd_0 + g_0_z_0_xxxyzz_0[i] * pb_x + g_0_z_0_xxxyzz_1[i] * wp_x[i];

        g_0_xz_0_xxxzzz_0[i] = 3.0 * g_0_z_0_xxzzz_1[i] * fi_abcd_0 + g_0_z_0_xxxzzz_0[i] * pb_x + g_0_z_0_xxxzzz_1[i] * wp_x[i];

        g_0_xz_0_xxyyyy_0[i] = g_0_x_0_xxyyyy_0[i] * pb_z + g_0_x_0_xxyyyy_1[i] * wp_z[i];

        g_0_xz_0_xxyyyz_0[i] = 2.0 * g_0_z_0_xyyyz_1[i] * fi_abcd_0 + g_0_z_0_xxyyyz_0[i] * pb_x + g_0_z_0_xxyyyz_1[i] * wp_x[i];

        g_0_xz_0_xxyyzz_0[i] = 2.0 * g_0_z_0_xyyzz_1[i] * fi_abcd_0 + g_0_z_0_xxyyzz_0[i] * pb_x + g_0_z_0_xxyyzz_1[i] * wp_x[i];

        g_0_xz_0_xxyzzz_0[i] = 2.0 * g_0_z_0_xyzzz_1[i] * fi_abcd_0 + g_0_z_0_xxyzzz_0[i] * pb_x + g_0_z_0_xxyzzz_1[i] * wp_x[i];

        g_0_xz_0_xxzzzz_0[i] = 2.0 * g_0_z_0_xzzzz_1[i] * fi_abcd_0 + g_0_z_0_xxzzzz_0[i] * pb_x + g_0_z_0_xxzzzz_1[i] * wp_x[i];

        g_0_xz_0_xyyyyy_0[i] = g_0_x_0_xyyyyy_0[i] * pb_z + g_0_x_0_xyyyyy_1[i] * wp_z[i];

        g_0_xz_0_xyyyyz_0[i] = g_0_z_0_yyyyz_1[i] * fi_abcd_0 + g_0_z_0_xyyyyz_0[i] * pb_x + g_0_z_0_xyyyyz_1[i] * wp_x[i];

        g_0_xz_0_xyyyzz_0[i] = g_0_z_0_yyyzz_1[i] * fi_abcd_0 + g_0_z_0_xyyyzz_0[i] * pb_x + g_0_z_0_xyyyzz_1[i] * wp_x[i];

        g_0_xz_0_xyyzzz_0[i] = g_0_z_0_yyzzz_1[i] * fi_abcd_0 + g_0_z_0_xyyzzz_0[i] * pb_x + g_0_z_0_xyyzzz_1[i] * wp_x[i];

        g_0_xz_0_xyzzzz_0[i] = g_0_z_0_yzzzz_1[i] * fi_abcd_0 + g_0_z_0_xyzzzz_0[i] * pb_x + g_0_z_0_xyzzzz_1[i] * wp_x[i];

        g_0_xz_0_xzzzzz_0[i] = g_0_z_0_zzzzz_1[i] * fi_abcd_0 + g_0_z_0_xzzzzz_0[i] * pb_x + g_0_z_0_xzzzzz_1[i] * wp_x[i];

        g_0_xz_0_yyyyyy_0[i] = g_0_z_0_yyyyyy_0[i] * pb_x + g_0_z_0_yyyyyy_1[i] * wp_x[i];

        g_0_xz_0_yyyyyz_0[i] = g_0_z_0_yyyyyz_0[i] * pb_x + g_0_z_0_yyyyyz_1[i] * wp_x[i];

        g_0_xz_0_yyyyzz_0[i] = g_0_z_0_yyyyzz_0[i] * pb_x + g_0_z_0_yyyyzz_1[i] * wp_x[i];

        g_0_xz_0_yyyzzz_0[i] = g_0_z_0_yyyzzz_0[i] * pb_x + g_0_z_0_yyyzzz_1[i] * wp_x[i];

        g_0_xz_0_yyzzzz_0[i] = g_0_z_0_yyzzzz_0[i] * pb_x + g_0_z_0_yyzzzz_1[i] * wp_x[i];

        g_0_xz_0_yzzzzz_0[i] = g_0_z_0_yzzzzz_0[i] * pb_x + g_0_z_0_yzzzzz_1[i] * wp_x[i];

        g_0_xz_0_zzzzzz_0[i] = g_0_z_0_zzzzzz_0[i] * pb_x + g_0_z_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 84-112 components of targeted buffer : SDSI

    auto g_0_yy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sdsi + 84);

    auto g_0_yy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sdsi + 85);

    auto g_0_yy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sdsi + 86);

    auto g_0_yy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sdsi + 87);

    auto g_0_yy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sdsi + 88);

    auto g_0_yy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sdsi + 89);

    auto g_0_yy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sdsi + 90);

    auto g_0_yy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sdsi + 91);

    auto g_0_yy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sdsi + 92);

    auto g_0_yy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sdsi + 93);

    auto g_0_yy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sdsi + 94);

    auto g_0_yy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sdsi + 95);

    auto g_0_yy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sdsi + 96);

    auto g_0_yy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sdsi + 97);

    auto g_0_yy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 98);

    auto g_0_yy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sdsi + 99);

    auto g_0_yy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sdsi + 100);

    auto g_0_yy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sdsi + 101);

    auto g_0_yy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sdsi + 102);

    auto g_0_yy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 103);

    auto g_0_yy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 104);

    auto g_0_yy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sdsi + 105);

    auto g_0_yy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sdsi + 106);

    auto g_0_yy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sdsi + 107);

    auto g_0_yy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sdsi + 108);

    auto g_0_yy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 109);

    auto g_0_yy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 110);

    auto g_0_yy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 111);

    #pragma omp simd aligned(g_0_0_0_xxxxxx_0, g_0_0_0_xxxxxx_1, g_0_0_0_xxxxxy_0, g_0_0_0_xxxxxy_1, g_0_0_0_xxxxxz_0, g_0_0_0_xxxxxz_1, g_0_0_0_xxxxyy_0, g_0_0_0_xxxxyy_1, g_0_0_0_xxxxyz_0, g_0_0_0_xxxxyz_1, g_0_0_0_xxxxzz_0, g_0_0_0_xxxxzz_1, g_0_0_0_xxxyyy_0, g_0_0_0_xxxyyy_1, g_0_0_0_xxxyyz_0, g_0_0_0_xxxyyz_1, g_0_0_0_xxxyzz_0, g_0_0_0_xxxyzz_1, g_0_0_0_xxxzzz_0, g_0_0_0_xxxzzz_1, g_0_0_0_xxyyyy_0, g_0_0_0_xxyyyy_1, g_0_0_0_xxyyyz_0, g_0_0_0_xxyyyz_1, g_0_0_0_xxyyzz_0, g_0_0_0_xxyyzz_1, g_0_0_0_xxyzzz_0, g_0_0_0_xxyzzz_1, g_0_0_0_xxzzzz_0, g_0_0_0_xxzzzz_1, g_0_0_0_xyyyyy_0, g_0_0_0_xyyyyy_1, g_0_0_0_xyyyyz_0, g_0_0_0_xyyyyz_1, g_0_0_0_xyyyzz_0, g_0_0_0_xyyyzz_1, g_0_0_0_xyyzzz_0, g_0_0_0_xyyzzz_1, g_0_0_0_xyzzzz_0, g_0_0_0_xyzzzz_1, g_0_0_0_xzzzzz_0, g_0_0_0_xzzzzz_1, g_0_0_0_yyyyyy_0, g_0_0_0_yyyyyy_1, g_0_0_0_yyyyyz_0, g_0_0_0_yyyyyz_1, g_0_0_0_yyyyzz_0, g_0_0_0_yyyyzz_1, g_0_0_0_yyyzzz_0, g_0_0_0_yyyzzz_1, g_0_0_0_yyzzzz_0, g_0_0_0_yyzzzz_1, g_0_0_0_yzzzzz_0, g_0_0_0_yzzzzz_1, g_0_0_0_zzzzzz_0, g_0_0_0_zzzzzz_1, g_0_y_0_xxxxx_1, g_0_y_0_xxxxxx_0, g_0_y_0_xxxxxx_1, g_0_y_0_xxxxxy_0, g_0_y_0_xxxxxy_1, g_0_y_0_xxxxxz_0, g_0_y_0_xxxxxz_1, g_0_y_0_xxxxy_1, g_0_y_0_xxxxyy_0, g_0_y_0_xxxxyy_1, g_0_y_0_xxxxyz_0, g_0_y_0_xxxxyz_1, g_0_y_0_xxxxz_1, g_0_y_0_xxxxzz_0, g_0_y_0_xxxxzz_1, g_0_y_0_xxxyy_1, g_0_y_0_xxxyyy_0, g_0_y_0_xxxyyy_1, g_0_y_0_xxxyyz_0, g_0_y_0_xxxyyz_1, g_0_y_0_xxxyz_1, g_0_y_0_xxxyzz_0, g_0_y_0_xxxyzz_1, g_0_y_0_xxxzz_1, g_0_y_0_xxxzzz_0, g_0_y_0_xxxzzz_1, g_0_y_0_xxyyy_1, g_0_y_0_xxyyyy_0, g_0_y_0_xxyyyy_1, g_0_y_0_xxyyyz_0, g_0_y_0_xxyyyz_1, g_0_y_0_xxyyz_1, g_0_y_0_xxyyzz_0, g_0_y_0_xxyyzz_1, g_0_y_0_xxyzz_1, g_0_y_0_xxyzzz_0, g_0_y_0_xxyzzz_1, g_0_y_0_xxzzz_1, g_0_y_0_xxzzzz_0, g_0_y_0_xxzzzz_1, g_0_y_0_xyyyy_1, g_0_y_0_xyyyyy_0, g_0_y_0_xyyyyy_1, g_0_y_0_xyyyyz_0, g_0_y_0_xyyyyz_1, g_0_y_0_xyyyz_1, g_0_y_0_xyyyzz_0, g_0_y_0_xyyyzz_1, g_0_y_0_xyyzz_1, g_0_y_0_xyyzzz_0, g_0_y_0_xyyzzz_1, g_0_y_0_xyzzz_1, g_0_y_0_xyzzzz_0, g_0_y_0_xyzzzz_1, g_0_y_0_xzzzz_1, g_0_y_0_xzzzzz_0, g_0_y_0_xzzzzz_1, g_0_y_0_yyyyy_1, g_0_y_0_yyyyyy_0, g_0_y_0_yyyyyy_1, g_0_y_0_yyyyyz_0, g_0_y_0_yyyyyz_1, g_0_y_0_yyyyz_1, g_0_y_0_yyyyzz_0, g_0_y_0_yyyyzz_1, g_0_y_0_yyyzz_1, g_0_y_0_yyyzzz_0, g_0_y_0_yyyzzz_1, g_0_y_0_yyzzz_1, g_0_y_0_yyzzzz_0, g_0_y_0_yyzzzz_1, g_0_y_0_yzzzz_1, g_0_y_0_yzzzzz_0, g_0_y_0_yzzzzz_1, g_0_y_0_zzzzz_1, g_0_y_0_zzzzzz_0, g_0_y_0_zzzzzz_1, g_0_yy_0_xxxxxx_0, g_0_yy_0_xxxxxy_0, g_0_yy_0_xxxxxz_0, g_0_yy_0_xxxxyy_0, g_0_yy_0_xxxxyz_0, g_0_yy_0_xxxxzz_0, g_0_yy_0_xxxyyy_0, g_0_yy_0_xxxyyz_0, g_0_yy_0_xxxyzz_0, g_0_yy_0_xxxzzz_0, g_0_yy_0_xxyyyy_0, g_0_yy_0_xxyyyz_0, g_0_yy_0_xxyyzz_0, g_0_yy_0_xxyzzz_0, g_0_yy_0_xxzzzz_0, g_0_yy_0_xyyyyy_0, g_0_yy_0_xyyyyz_0, g_0_yy_0_xyyyzz_0, g_0_yy_0_xyyzzz_0, g_0_yy_0_xyzzzz_0, g_0_yy_0_xzzzzz_0, g_0_yy_0_yyyyyy_0, g_0_yy_0_yyyyyz_0, g_0_yy_0_yyyyzz_0, g_0_yy_0_yyyzzz_0, g_0_yy_0_yyzzzz_0, g_0_yy_0_yzzzzz_0, g_0_yy_0_zzzzzz_0, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yy_0_xxxxxx_0[i] = g_0_0_0_xxxxxx_0[i] * fi_ab_0 - g_0_0_0_xxxxxx_1[i] * fti_ab_0 + g_0_y_0_xxxxxx_0[i] * pb_y + g_0_y_0_xxxxxx_1[i] * wp_y[i];

        g_0_yy_0_xxxxxy_0[i] = g_0_0_0_xxxxxy_0[i] * fi_ab_0 - g_0_0_0_xxxxxy_1[i] * fti_ab_0 + g_0_y_0_xxxxx_1[i] * fi_abcd_0 + g_0_y_0_xxxxxy_0[i] * pb_y + g_0_y_0_xxxxxy_1[i] * wp_y[i];

        g_0_yy_0_xxxxxz_0[i] = g_0_0_0_xxxxxz_0[i] * fi_ab_0 - g_0_0_0_xxxxxz_1[i] * fti_ab_0 + g_0_y_0_xxxxxz_0[i] * pb_y + g_0_y_0_xxxxxz_1[i] * wp_y[i];

        g_0_yy_0_xxxxyy_0[i] = g_0_0_0_xxxxyy_0[i] * fi_ab_0 - g_0_0_0_xxxxyy_1[i] * fti_ab_0 + 2.0 * g_0_y_0_xxxxy_1[i] * fi_abcd_0 + g_0_y_0_xxxxyy_0[i] * pb_y + g_0_y_0_xxxxyy_1[i] * wp_y[i];

        g_0_yy_0_xxxxyz_0[i] = g_0_0_0_xxxxyz_0[i] * fi_ab_0 - g_0_0_0_xxxxyz_1[i] * fti_ab_0 + g_0_y_0_xxxxz_1[i] * fi_abcd_0 + g_0_y_0_xxxxyz_0[i] * pb_y + g_0_y_0_xxxxyz_1[i] * wp_y[i];

        g_0_yy_0_xxxxzz_0[i] = g_0_0_0_xxxxzz_0[i] * fi_ab_0 - g_0_0_0_xxxxzz_1[i] * fti_ab_0 + g_0_y_0_xxxxzz_0[i] * pb_y + g_0_y_0_xxxxzz_1[i] * wp_y[i];

        g_0_yy_0_xxxyyy_0[i] = g_0_0_0_xxxyyy_0[i] * fi_ab_0 - g_0_0_0_xxxyyy_1[i] * fti_ab_0 + 3.0 * g_0_y_0_xxxyy_1[i] * fi_abcd_0 + g_0_y_0_xxxyyy_0[i] * pb_y + g_0_y_0_xxxyyy_1[i] * wp_y[i];

        g_0_yy_0_xxxyyz_0[i] = g_0_0_0_xxxyyz_0[i] * fi_ab_0 - g_0_0_0_xxxyyz_1[i] * fti_ab_0 + 2.0 * g_0_y_0_xxxyz_1[i] * fi_abcd_0 + g_0_y_0_xxxyyz_0[i] * pb_y + g_0_y_0_xxxyyz_1[i] * wp_y[i];

        g_0_yy_0_xxxyzz_0[i] = g_0_0_0_xxxyzz_0[i] * fi_ab_0 - g_0_0_0_xxxyzz_1[i] * fti_ab_0 + g_0_y_0_xxxzz_1[i] * fi_abcd_0 + g_0_y_0_xxxyzz_0[i] * pb_y + g_0_y_0_xxxyzz_1[i] * wp_y[i];

        g_0_yy_0_xxxzzz_0[i] = g_0_0_0_xxxzzz_0[i] * fi_ab_0 - g_0_0_0_xxxzzz_1[i] * fti_ab_0 + g_0_y_0_xxxzzz_0[i] * pb_y + g_0_y_0_xxxzzz_1[i] * wp_y[i];

        g_0_yy_0_xxyyyy_0[i] = g_0_0_0_xxyyyy_0[i] * fi_ab_0 - g_0_0_0_xxyyyy_1[i] * fti_ab_0 + 4.0 * g_0_y_0_xxyyy_1[i] * fi_abcd_0 + g_0_y_0_xxyyyy_0[i] * pb_y + g_0_y_0_xxyyyy_1[i] * wp_y[i];

        g_0_yy_0_xxyyyz_0[i] = g_0_0_0_xxyyyz_0[i] * fi_ab_0 - g_0_0_0_xxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_y_0_xxyyz_1[i] * fi_abcd_0 + g_0_y_0_xxyyyz_0[i] * pb_y + g_0_y_0_xxyyyz_1[i] * wp_y[i];

        g_0_yy_0_xxyyzz_0[i] = g_0_0_0_xxyyzz_0[i] * fi_ab_0 - g_0_0_0_xxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_y_0_xxyzz_1[i] * fi_abcd_0 + g_0_y_0_xxyyzz_0[i] * pb_y + g_0_y_0_xxyyzz_1[i] * wp_y[i];

        g_0_yy_0_xxyzzz_0[i] = g_0_0_0_xxyzzz_0[i] * fi_ab_0 - g_0_0_0_xxyzzz_1[i] * fti_ab_0 + g_0_y_0_xxzzz_1[i] * fi_abcd_0 + g_0_y_0_xxyzzz_0[i] * pb_y + g_0_y_0_xxyzzz_1[i] * wp_y[i];

        g_0_yy_0_xxzzzz_0[i] = g_0_0_0_xxzzzz_0[i] * fi_ab_0 - g_0_0_0_xxzzzz_1[i] * fti_ab_0 + g_0_y_0_xxzzzz_0[i] * pb_y + g_0_y_0_xxzzzz_1[i] * wp_y[i];

        g_0_yy_0_xyyyyy_0[i] = g_0_0_0_xyyyyy_0[i] * fi_ab_0 - g_0_0_0_xyyyyy_1[i] * fti_ab_0 + 5.0 * g_0_y_0_xyyyy_1[i] * fi_abcd_0 + g_0_y_0_xyyyyy_0[i] * pb_y + g_0_y_0_xyyyyy_1[i] * wp_y[i];

        g_0_yy_0_xyyyyz_0[i] = g_0_0_0_xyyyyz_0[i] * fi_ab_0 - g_0_0_0_xyyyyz_1[i] * fti_ab_0 + 4.0 * g_0_y_0_xyyyz_1[i] * fi_abcd_0 + g_0_y_0_xyyyyz_0[i] * pb_y + g_0_y_0_xyyyyz_1[i] * wp_y[i];

        g_0_yy_0_xyyyzz_0[i] = g_0_0_0_xyyyzz_0[i] * fi_ab_0 - g_0_0_0_xyyyzz_1[i] * fti_ab_0 + 3.0 * g_0_y_0_xyyzz_1[i] * fi_abcd_0 + g_0_y_0_xyyyzz_0[i] * pb_y + g_0_y_0_xyyyzz_1[i] * wp_y[i];

        g_0_yy_0_xyyzzz_0[i] = g_0_0_0_xyyzzz_0[i] * fi_ab_0 - g_0_0_0_xyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_y_0_xyzzz_1[i] * fi_abcd_0 + g_0_y_0_xyyzzz_0[i] * pb_y + g_0_y_0_xyyzzz_1[i] * wp_y[i];

        g_0_yy_0_xyzzzz_0[i] = g_0_0_0_xyzzzz_0[i] * fi_ab_0 - g_0_0_0_xyzzzz_1[i] * fti_ab_0 + g_0_y_0_xzzzz_1[i] * fi_abcd_0 + g_0_y_0_xyzzzz_0[i] * pb_y + g_0_y_0_xyzzzz_1[i] * wp_y[i];

        g_0_yy_0_xzzzzz_0[i] = g_0_0_0_xzzzzz_0[i] * fi_ab_0 - g_0_0_0_xzzzzz_1[i] * fti_ab_0 + g_0_y_0_xzzzzz_0[i] * pb_y + g_0_y_0_xzzzzz_1[i] * wp_y[i];

        g_0_yy_0_yyyyyy_0[i] = g_0_0_0_yyyyyy_0[i] * fi_ab_0 - g_0_0_0_yyyyyy_1[i] * fti_ab_0 + 6.0 * g_0_y_0_yyyyy_1[i] * fi_abcd_0 + g_0_y_0_yyyyyy_0[i] * pb_y + g_0_y_0_yyyyyy_1[i] * wp_y[i];

        g_0_yy_0_yyyyyz_0[i] = g_0_0_0_yyyyyz_0[i] * fi_ab_0 - g_0_0_0_yyyyyz_1[i] * fti_ab_0 + 5.0 * g_0_y_0_yyyyz_1[i] * fi_abcd_0 + g_0_y_0_yyyyyz_0[i] * pb_y + g_0_y_0_yyyyyz_1[i] * wp_y[i];

        g_0_yy_0_yyyyzz_0[i] = g_0_0_0_yyyyzz_0[i] * fi_ab_0 - g_0_0_0_yyyyzz_1[i] * fti_ab_0 + 4.0 * g_0_y_0_yyyzz_1[i] * fi_abcd_0 + g_0_y_0_yyyyzz_0[i] * pb_y + g_0_y_0_yyyyzz_1[i] * wp_y[i];

        g_0_yy_0_yyyzzz_0[i] = g_0_0_0_yyyzzz_0[i] * fi_ab_0 - g_0_0_0_yyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_y_0_yyzzz_1[i] * fi_abcd_0 + g_0_y_0_yyyzzz_0[i] * pb_y + g_0_y_0_yyyzzz_1[i] * wp_y[i];

        g_0_yy_0_yyzzzz_0[i] = g_0_0_0_yyzzzz_0[i] * fi_ab_0 - g_0_0_0_yyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_y_0_yzzzz_1[i] * fi_abcd_0 + g_0_y_0_yyzzzz_0[i] * pb_y + g_0_y_0_yyzzzz_1[i] * wp_y[i];

        g_0_yy_0_yzzzzz_0[i] = g_0_0_0_yzzzzz_0[i] * fi_ab_0 - g_0_0_0_yzzzzz_1[i] * fti_ab_0 + g_0_y_0_zzzzz_1[i] * fi_abcd_0 + g_0_y_0_yzzzzz_0[i] * pb_y + g_0_y_0_yzzzzz_1[i] * wp_y[i];

        g_0_yy_0_zzzzzz_0[i] = g_0_0_0_zzzzzz_0[i] * fi_ab_0 - g_0_0_0_zzzzzz_1[i] * fti_ab_0 + g_0_y_0_zzzzzz_0[i] * pb_y + g_0_y_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 112-140 components of targeted buffer : SDSI

    auto g_0_yz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sdsi + 112);

    auto g_0_yz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sdsi + 113);

    auto g_0_yz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sdsi + 114);

    auto g_0_yz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sdsi + 115);

    auto g_0_yz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sdsi + 116);

    auto g_0_yz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sdsi + 117);

    auto g_0_yz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sdsi + 118);

    auto g_0_yz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sdsi + 119);

    auto g_0_yz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sdsi + 120);

    auto g_0_yz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sdsi + 121);

    auto g_0_yz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sdsi + 122);

    auto g_0_yz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sdsi + 123);

    auto g_0_yz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sdsi + 124);

    auto g_0_yz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sdsi + 125);

    auto g_0_yz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 126);

    auto g_0_yz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sdsi + 127);

    auto g_0_yz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sdsi + 128);

    auto g_0_yz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sdsi + 129);

    auto g_0_yz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sdsi + 130);

    auto g_0_yz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 131);

    auto g_0_yz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 132);

    auto g_0_yz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sdsi + 133);

    auto g_0_yz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sdsi + 134);

    auto g_0_yz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sdsi + 135);

    auto g_0_yz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sdsi + 136);

    auto g_0_yz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 137);

    auto g_0_yz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 138);

    auto g_0_yz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 139);

    #pragma omp simd aligned(g_0_y_0_xxxxxy_0, g_0_y_0_xxxxxy_1, g_0_y_0_xxxxyy_0, g_0_y_0_xxxxyy_1, g_0_y_0_xxxyyy_0, g_0_y_0_xxxyyy_1, g_0_y_0_xxyyyy_0, g_0_y_0_xxyyyy_1, g_0_y_0_xyyyyy_0, g_0_y_0_xyyyyy_1, g_0_y_0_yyyyyy_0, g_0_y_0_yyyyyy_1, g_0_yz_0_xxxxxx_0, g_0_yz_0_xxxxxy_0, g_0_yz_0_xxxxxz_0, g_0_yz_0_xxxxyy_0, g_0_yz_0_xxxxyz_0, g_0_yz_0_xxxxzz_0, g_0_yz_0_xxxyyy_0, g_0_yz_0_xxxyyz_0, g_0_yz_0_xxxyzz_0, g_0_yz_0_xxxzzz_0, g_0_yz_0_xxyyyy_0, g_0_yz_0_xxyyyz_0, g_0_yz_0_xxyyzz_0, g_0_yz_0_xxyzzz_0, g_0_yz_0_xxzzzz_0, g_0_yz_0_xyyyyy_0, g_0_yz_0_xyyyyz_0, g_0_yz_0_xyyyzz_0, g_0_yz_0_xyyzzz_0, g_0_yz_0_xyzzzz_0, g_0_yz_0_xzzzzz_0, g_0_yz_0_yyyyyy_0, g_0_yz_0_yyyyyz_0, g_0_yz_0_yyyyzz_0, g_0_yz_0_yyyzzz_0, g_0_yz_0_yyzzzz_0, g_0_yz_0_yzzzzz_0, g_0_yz_0_zzzzzz_0, g_0_z_0_xxxxxx_0, g_0_z_0_xxxxxx_1, g_0_z_0_xxxxxz_0, g_0_z_0_xxxxxz_1, g_0_z_0_xxxxyz_0, g_0_z_0_xxxxyz_1, g_0_z_0_xxxxz_1, g_0_z_0_xxxxzz_0, g_0_z_0_xxxxzz_1, g_0_z_0_xxxyyz_0, g_0_z_0_xxxyyz_1, g_0_z_0_xxxyz_1, g_0_z_0_xxxyzz_0, g_0_z_0_xxxyzz_1, g_0_z_0_xxxzz_1, g_0_z_0_xxxzzz_0, g_0_z_0_xxxzzz_1, g_0_z_0_xxyyyz_0, g_0_z_0_xxyyyz_1, g_0_z_0_xxyyz_1, g_0_z_0_xxyyzz_0, g_0_z_0_xxyyzz_1, g_0_z_0_xxyzz_1, g_0_z_0_xxyzzz_0, g_0_z_0_xxyzzz_1, g_0_z_0_xxzzz_1, g_0_z_0_xxzzzz_0, g_0_z_0_xxzzzz_1, g_0_z_0_xyyyyz_0, g_0_z_0_xyyyyz_1, g_0_z_0_xyyyz_1, g_0_z_0_xyyyzz_0, g_0_z_0_xyyyzz_1, g_0_z_0_xyyzz_1, g_0_z_0_xyyzzz_0, g_0_z_0_xyyzzz_1, g_0_z_0_xyzzz_1, g_0_z_0_xyzzzz_0, g_0_z_0_xyzzzz_1, g_0_z_0_xzzzz_1, g_0_z_0_xzzzzz_0, g_0_z_0_xzzzzz_1, g_0_z_0_yyyyyz_0, g_0_z_0_yyyyyz_1, g_0_z_0_yyyyz_1, g_0_z_0_yyyyzz_0, g_0_z_0_yyyyzz_1, g_0_z_0_yyyzz_1, g_0_z_0_yyyzzz_0, g_0_z_0_yyyzzz_1, g_0_z_0_yyzzz_1, g_0_z_0_yyzzzz_0, g_0_z_0_yyzzzz_1, g_0_z_0_yzzzz_1, g_0_z_0_yzzzzz_0, g_0_z_0_yzzzzz_1, g_0_z_0_zzzzz_1, g_0_z_0_zzzzzz_0, g_0_z_0_zzzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yz_0_xxxxxx_0[i] = g_0_z_0_xxxxxx_0[i] * pb_y + g_0_z_0_xxxxxx_1[i] * wp_y[i];

        g_0_yz_0_xxxxxy_0[i] = g_0_y_0_xxxxxy_0[i] * pb_z + g_0_y_0_xxxxxy_1[i] * wp_z[i];

        g_0_yz_0_xxxxxz_0[i] = g_0_z_0_xxxxxz_0[i] * pb_y + g_0_z_0_xxxxxz_1[i] * wp_y[i];

        g_0_yz_0_xxxxyy_0[i] = g_0_y_0_xxxxyy_0[i] * pb_z + g_0_y_0_xxxxyy_1[i] * wp_z[i];

        g_0_yz_0_xxxxyz_0[i] = g_0_z_0_xxxxz_1[i] * fi_abcd_0 + g_0_z_0_xxxxyz_0[i] * pb_y + g_0_z_0_xxxxyz_1[i] * wp_y[i];

        g_0_yz_0_xxxxzz_0[i] = g_0_z_0_xxxxzz_0[i] * pb_y + g_0_z_0_xxxxzz_1[i] * wp_y[i];

        g_0_yz_0_xxxyyy_0[i] = g_0_y_0_xxxyyy_0[i] * pb_z + g_0_y_0_xxxyyy_1[i] * wp_z[i];

        g_0_yz_0_xxxyyz_0[i] = 2.0 * g_0_z_0_xxxyz_1[i] * fi_abcd_0 + g_0_z_0_xxxyyz_0[i] * pb_y + g_0_z_0_xxxyyz_1[i] * wp_y[i];

        g_0_yz_0_xxxyzz_0[i] = g_0_z_0_xxxzz_1[i] * fi_abcd_0 + g_0_z_0_xxxyzz_0[i] * pb_y + g_0_z_0_xxxyzz_1[i] * wp_y[i];

        g_0_yz_0_xxxzzz_0[i] = g_0_z_0_xxxzzz_0[i] * pb_y + g_0_z_0_xxxzzz_1[i] * wp_y[i];

        g_0_yz_0_xxyyyy_0[i] = g_0_y_0_xxyyyy_0[i] * pb_z + g_0_y_0_xxyyyy_1[i] * wp_z[i];

        g_0_yz_0_xxyyyz_0[i] = 3.0 * g_0_z_0_xxyyz_1[i] * fi_abcd_0 + g_0_z_0_xxyyyz_0[i] * pb_y + g_0_z_0_xxyyyz_1[i] * wp_y[i];

        g_0_yz_0_xxyyzz_0[i] = 2.0 * g_0_z_0_xxyzz_1[i] * fi_abcd_0 + g_0_z_0_xxyyzz_0[i] * pb_y + g_0_z_0_xxyyzz_1[i] * wp_y[i];

        g_0_yz_0_xxyzzz_0[i] = g_0_z_0_xxzzz_1[i] * fi_abcd_0 + g_0_z_0_xxyzzz_0[i] * pb_y + g_0_z_0_xxyzzz_1[i] * wp_y[i];

        g_0_yz_0_xxzzzz_0[i] = g_0_z_0_xxzzzz_0[i] * pb_y + g_0_z_0_xxzzzz_1[i] * wp_y[i];

        g_0_yz_0_xyyyyy_0[i] = g_0_y_0_xyyyyy_0[i] * pb_z + g_0_y_0_xyyyyy_1[i] * wp_z[i];

        g_0_yz_0_xyyyyz_0[i] = 4.0 * g_0_z_0_xyyyz_1[i] * fi_abcd_0 + g_0_z_0_xyyyyz_0[i] * pb_y + g_0_z_0_xyyyyz_1[i] * wp_y[i];

        g_0_yz_0_xyyyzz_0[i] = 3.0 * g_0_z_0_xyyzz_1[i] * fi_abcd_0 + g_0_z_0_xyyyzz_0[i] * pb_y + g_0_z_0_xyyyzz_1[i] * wp_y[i];

        g_0_yz_0_xyyzzz_0[i] = 2.0 * g_0_z_0_xyzzz_1[i] * fi_abcd_0 + g_0_z_0_xyyzzz_0[i] * pb_y + g_0_z_0_xyyzzz_1[i] * wp_y[i];

        g_0_yz_0_xyzzzz_0[i] = g_0_z_0_xzzzz_1[i] * fi_abcd_0 + g_0_z_0_xyzzzz_0[i] * pb_y + g_0_z_0_xyzzzz_1[i] * wp_y[i];

        g_0_yz_0_xzzzzz_0[i] = g_0_z_0_xzzzzz_0[i] * pb_y + g_0_z_0_xzzzzz_1[i] * wp_y[i];

        g_0_yz_0_yyyyyy_0[i] = g_0_y_0_yyyyyy_0[i] * pb_z + g_0_y_0_yyyyyy_1[i] * wp_z[i];

        g_0_yz_0_yyyyyz_0[i] = 5.0 * g_0_z_0_yyyyz_1[i] * fi_abcd_0 + g_0_z_0_yyyyyz_0[i] * pb_y + g_0_z_0_yyyyyz_1[i] * wp_y[i];

        g_0_yz_0_yyyyzz_0[i] = 4.0 * g_0_z_0_yyyzz_1[i] * fi_abcd_0 + g_0_z_0_yyyyzz_0[i] * pb_y + g_0_z_0_yyyyzz_1[i] * wp_y[i];

        g_0_yz_0_yyyzzz_0[i] = 3.0 * g_0_z_0_yyzzz_1[i] * fi_abcd_0 + g_0_z_0_yyyzzz_0[i] * pb_y + g_0_z_0_yyyzzz_1[i] * wp_y[i];

        g_0_yz_0_yyzzzz_0[i] = 2.0 * g_0_z_0_yzzzz_1[i] * fi_abcd_0 + g_0_z_0_yyzzzz_0[i] * pb_y + g_0_z_0_yyzzzz_1[i] * wp_y[i];

        g_0_yz_0_yzzzzz_0[i] = g_0_z_0_zzzzz_1[i] * fi_abcd_0 + g_0_z_0_yzzzzz_0[i] * pb_y + g_0_z_0_yzzzzz_1[i] * wp_y[i];

        g_0_yz_0_zzzzzz_0[i] = g_0_z_0_zzzzzz_0[i] * pb_y + g_0_z_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 140-168 components of targeted buffer : SDSI

    auto g_0_zz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sdsi + 140);

    auto g_0_zz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sdsi + 141);

    auto g_0_zz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sdsi + 142);

    auto g_0_zz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sdsi + 143);

    auto g_0_zz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sdsi + 144);

    auto g_0_zz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sdsi + 145);

    auto g_0_zz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sdsi + 146);

    auto g_0_zz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sdsi + 147);

    auto g_0_zz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sdsi + 148);

    auto g_0_zz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sdsi + 149);

    auto g_0_zz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sdsi + 150);

    auto g_0_zz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sdsi + 151);

    auto g_0_zz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sdsi + 152);

    auto g_0_zz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sdsi + 153);

    auto g_0_zz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 154);

    auto g_0_zz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sdsi + 155);

    auto g_0_zz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sdsi + 156);

    auto g_0_zz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sdsi + 157);

    auto g_0_zz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sdsi + 158);

    auto g_0_zz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 159);

    auto g_0_zz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 160);

    auto g_0_zz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sdsi + 161);

    auto g_0_zz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sdsi + 162);

    auto g_0_zz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sdsi + 163);

    auto g_0_zz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sdsi + 164);

    auto g_0_zz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 165);

    auto g_0_zz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 166);

    auto g_0_zz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 167);

    #pragma omp simd aligned(g_0_0_0_xxxxxx_0, g_0_0_0_xxxxxx_1, g_0_0_0_xxxxxy_0, g_0_0_0_xxxxxy_1, g_0_0_0_xxxxxz_0, g_0_0_0_xxxxxz_1, g_0_0_0_xxxxyy_0, g_0_0_0_xxxxyy_1, g_0_0_0_xxxxyz_0, g_0_0_0_xxxxyz_1, g_0_0_0_xxxxzz_0, g_0_0_0_xxxxzz_1, g_0_0_0_xxxyyy_0, g_0_0_0_xxxyyy_1, g_0_0_0_xxxyyz_0, g_0_0_0_xxxyyz_1, g_0_0_0_xxxyzz_0, g_0_0_0_xxxyzz_1, g_0_0_0_xxxzzz_0, g_0_0_0_xxxzzz_1, g_0_0_0_xxyyyy_0, g_0_0_0_xxyyyy_1, g_0_0_0_xxyyyz_0, g_0_0_0_xxyyyz_1, g_0_0_0_xxyyzz_0, g_0_0_0_xxyyzz_1, g_0_0_0_xxyzzz_0, g_0_0_0_xxyzzz_1, g_0_0_0_xxzzzz_0, g_0_0_0_xxzzzz_1, g_0_0_0_xyyyyy_0, g_0_0_0_xyyyyy_1, g_0_0_0_xyyyyz_0, g_0_0_0_xyyyyz_1, g_0_0_0_xyyyzz_0, g_0_0_0_xyyyzz_1, g_0_0_0_xyyzzz_0, g_0_0_0_xyyzzz_1, g_0_0_0_xyzzzz_0, g_0_0_0_xyzzzz_1, g_0_0_0_xzzzzz_0, g_0_0_0_xzzzzz_1, g_0_0_0_yyyyyy_0, g_0_0_0_yyyyyy_1, g_0_0_0_yyyyyz_0, g_0_0_0_yyyyyz_1, g_0_0_0_yyyyzz_0, g_0_0_0_yyyyzz_1, g_0_0_0_yyyzzz_0, g_0_0_0_yyyzzz_1, g_0_0_0_yyzzzz_0, g_0_0_0_yyzzzz_1, g_0_0_0_yzzzzz_0, g_0_0_0_yzzzzz_1, g_0_0_0_zzzzzz_0, g_0_0_0_zzzzzz_1, g_0_z_0_xxxxx_1, g_0_z_0_xxxxxx_0, g_0_z_0_xxxxxx_1, g_0_z_0_xxxxxy_0, g_0_z_0_xxxxxy_1, g_0_z_0_xxxxxz_0, g_0_z_0_xxxxxz_1, g_0_z_0_xxxxy_1, g_0_z_0_xxxxyy_0, g_0_z_0_xxxxyy_1, g_0_z_0_xxxxyz_0, g_0_z_0_xxxxyz_1, g_0_z_0_xxxxz_1, g_0_z_0_xxxxzz_0, g_0_z_0_xxxxzz_1, g_0_z_0_xxxyy_1, g_0_z_0_xxxyyy_0, g_0_z_0_xxxyyy_1, g_0_z_0_xxxyyz_0, g_0_z_0_xxxyyz_1, g_0_z_0_xxxyz_1, g_0_z_0_xxxyzz_0, g_0_z_0_xxxyzz_1, g_0_z_0_xxxzz_1, g_0_z_0_xxxzzz_0, g_0_z_0_xxxzzz_1, g_0_z_0_xxyyy_1, g_0_z_0_xxyyyy_0, g_0_z_0_xxyyyy_1, g_0_z_0_xxyyyz_0, g_0_z_0_xxyyyz_1, g_0_z_0_xxyyz_1, g_0_z_0_xxyyzz_0, g_0_z_0_xxyyzz_1, g_0_z_0_xxyzz_1, g_0_z_0_xxyzzz_0, g_0_z_0_xxyzzz_1, g_0_z_0_xxzzz_1, g_0_z_0_xxzzzz_0, g_0_z_0_xxzzzz_1, g_0_z_0_xyyyy_1, g_0_z_0_xyyyyy_0, g_0_z_0_xyyyyy_1, g_0_z_0_xyyyyz_0, g_0_z_0_xyyyyz_1, g_0_z_0_xyyyz_1, g_0_z_0_xyyyzz_0, g_0_z_0_xyyyzz_1, g_0_z_0_xyyzz_1, g_0_z_0_xyyzzz_0, g_0_z_0_xyyzzz_1, g_0_z_0_xyzzz_1, g_0_z_0_xyzzzz_0, g_0_z_0_xyzzzz_1, g_0_z_0_xzzzz_1, g_0_z_0_xzzzzz_0, g_0_z_0_xzzzzz_1, g_0_z_0_yyyyy_1, g_0_z_0_yyyyyy_0, g_0_z_0_yyyyyy_1, g_0_z_0_yyyyyz_0, g_0_z_0_yyyyyz_1, g_0_z_0_yyyyz_1, g_0_z_0_yyyyzz_0, g_0_z_0_yyyyzz_1, g_0_z_0_yyyzz_1, g_0_z_0_yyyzzz_0, g_0_z_0_yyyzzz_1, g_0_z_0_yyzzz_1, g_0_z_0_yyzzzz_0, g_0_z_0_yyzzzz_1, g_0_z_0_yzzzz_1, g_0_z_0_yzzzzz_0, g_0_z_0_yzzzzz_1, g_0_z_0_zzzzz_1, g_0_z_0_zzzzzz_0, g_0_z_0_zzzzzz_1, g_0_zz_0_xxxxxx_0, g_0_zz_0_xxxxxy_0, g_0_zz_0_xxxxxz_0, g_0_zz_0_xxxxyy_0, g_0_zz_0_xxxxyz_0, g_0_zz_0_xxxxzz_0, g_0_zz_0_xxxyyy_0, g_0_zz_0_xxxyyz_0, g_0_zz_0_xxxyzz_0, g_0_zz_0_xxxzzz_0, g_0_zz_0_xxyyyy_0, g_0_zz_0_xxyyyz_0, g_0_zz_0_xxyyzz_0, g_0_zz_0_xxyzzz_0, g_0_zz_0_xxzzzz_0, g_0_zz_0_xyyyyy_0, g_0_zz_0_xyyyyz_0, g_0_zz_0_xyyyzz_0, g_0_zz_0_xyyzzz_0, g_0_zz_0_xyzzzz_0, g_0_zz_0_xzzzzz_0, g_0_zz_0_yyyyyy_0, g_0_zz_0_yyyyyz_0, g_0_zz_0_yyyyzz_0, g_0_zz_0_yyyzzz_0, g_0_zz_0_yyzzzz_0, g_0_zz_0_yzzzzz_0, g_0_zz_0_zzzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zz_0_xxxxxx_0[i] = g_0_0_0_xxxxxx_0[i] * fi_ab_0 - g_0_0_0_xxxxxx_1[i] * fti_ab_0 + g_0_z_0_xxxxxx_0[i] * pb_z + g_0_z_0_xxxxxx_1[i] * wp_z[i];

        g_0_zz_0_xxxxxy_0[i] = g_0_0_0_xxxxxy_0[i] * fi_ab_0 - g_0_0_0_xxxxxy_1[i] * fti_ab_0 + g_0_z_0_xxxxxy_0[i] * pb_z + g_0_z_0_xxxxxy_1[i] * wp_z[i];

        g_0_zz_0_xxxxxz_0[i] = g_0_0_0_xxxxxz_0[i] * fi_ab_0 - g_0_0_0_xxxxxz_1[i] * fti_ab_0 + g_0_z_0_xxxxx_1[i] * fi_abcd_0 + g_0_z_0_xxxxxz_0[i] * pb_z + g_0_z_0_xxxxxz_1[i] * wp_z[i];

        g_0_zz_0_xxxxyy_0[i] = g_0_0_0_xxxxyy_0[i] * fi_ab_0 - g_0_0_0_xxxxyy_1[i] * fti_ab_0 + g_0_z_0_xxxxyy_0[i] * pb_z + g_0_z_0_xxxxyy_1[i] * wp_z[i];

        g_0_zz_0_xxxxyz_0[i] = g_0_0_0_xxxxyz_0[i] * fi_ab_0 - g_0_0_0_xxxxyz_1[i] * fti_ab_0 + g_0_z_0_xxxxy_1[i] * fi_abcd_0 + g_0_z_0_xxxxyz_0[i] * pb_z + g_0_z_0_xxxxyz_1[i] * wp_z[i];

        g_0_zz_0_xxxxzz_0[i] = g_0_0_0_xxxxzz_0[i] * fi_ab_0 - g_0_0_0_xxxxzz_1[i] * fti_ab_0 + 2.0 * g_0_z_0_xxxxz_1[i] * fi_abcd_0 + g_0_z_0_xxxxzz_0[i] * pb_z + g_0_z_0_xxxxzz_1[i] * wp_z[i];

        g_0_zz_0_xxxyyy_0[i] = g_0_0_0_xxxyyy_0[i] * fi_ab_0 - g_0_0_0_xxxyyy_1[i] * fti_ab_0 + g_0_z_0_xxxyyy_0[i] * pb_z + g_0_z_0_xxxyyy_1[i] * wp_z[i];

        g_0_zz_0_xxxyyz_0[i] = g_0_0_0_xxxyyz_0[i] * fi_ab_0 - g_0_0_0_xxxyyz_1[i] * fti_ab_0 + g_0_z_0_xxxyy_1[i] * fi_abcd_0 + g_0_z_0_xxxyyz_0[i] * pb_z + g_0_z_0_xxxyyz_1[i] * wp_z[i];

        g_0_zz_0_xxxyzz_0[i] = g_0_0_0_xxxyzz_0[i] * fi_ab_0 - g_0_0_0_xxxyzz_1[i] * fti_ab_0 + 2.0 * g_0_z_0_xxxyz_1[i] * fi_abcd_0 + g_0_z_0_xxxyzz_0[i] * pb_z + g_0_z_0_xxxyzz_1[i] * wp_z[i];

        g_0_zz_0_xxxzzz_0[i] = g_0_0_0_xxxzzz_0[i] * fi_ab_0 - g_0_0_0_xxxzzz_1[i] * fti_ab_0 + 3.0 * g_0_z_0_xxxzz_1[i] * fi_abcd_0 + g_0_z_0_xxxzzz_0[i] * pb_z + g_0_z_0_xxxzzz_1[i] * wp_z[i];

        g_0_zz_0_xxyyyy_0[i] = g_0_0_0_xxyyyy_0[i] * fi_ab_0 - g_0_0_0_xxyyyy_1[i] * fti_ab_0 + g_0_z_0_xxyyyy_0[i] * pb_z + g_0_z_0_xxyyyy_1[i] * wp_z[i];

        g_0_zz_0_xxyyyz_0[i] = g_0_0_0_xxyyyz_0[i] * fi_ab_0 - g_0_0_0_xxyyyz_1[i] * fti_ab_0 + g_0_z_0_xxyyy_1[i] * fi_abcd_0 + g_0_z_0_xxyyyz_0[i] * pb_z + g_0_z_0_xxyyyz_1[i] * wp_z[i];

        g_0_zz_0_xxyyzz_0[i] = g_0_0_0_xxyyzz_0[i] * fi_ab_0 - g_0_0_0_xxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_z_0_xxyyz_1[i] * fi_abcd_0 + g_0_z_0_xxyyzz_0[i] * pb_z + g_0_z_0_xxyyzz_1[i] * wp_z[i];

        g_0_zz_0_xxyzzz_0[i] = g_0_0_0_xxyzzz_0[i] * fi_ab_0 - g_0_0_0_xxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_z_0_xxyzz_1[i] * fi_abcd_0 + g_0_z_0_xxyzzz_0[i] * pb_z + g_0_z_0_xxyzzz_1[i] * wp_z[i];

        g_0_zz_0_xxzzzz_0[i] = g_0_0_0_xxzzzz_0[i] * fi_ab_0 - g_0_0_0_xxzzzz_1[i] * fti_ab_0 + 4.0 * g_0_z_0_xxzzz_1[i] * fi_abcd_0 + g_0_z_0_xxzzzz_0[i] * pb_z + g_0_z_0_xxzzzz_1[i] * wp_z[i];

        g_0_zz_0_xyyyyy_0[i] = g_0_0_0_xyyyyy_0[i] * fi_ab_0 - g_0_0_0_xyyyyy_1[i] * fti_ab_0 + g_0_z_0_xyyyyy_0[i] * pb_z + g_0_z_0_xyyyyy_1[i] * wp_z[i];

        g_0_zz_0_xyyyyz_0[i] = g_0_0_0_xyyyyz_0[i] * fi_ab_0 - g_0_0_0_xyyyyz_1[i] * fti_ab_0 + g_0_z_0_xyyyy_1[i] * fi_abcd_0 + g_0_z_0_xyyyyz_0[i] * pb_z + g_0_z_0_xyyyyz_1[i] * wp_z[i];

        g_0_zz_0_xyyyzz_0[i] = g_0_0_0_xyyyzz_0[i] * fi_ab_0 - g_0_0_0_xyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_z_0_xyyyz_1[i] * fi_abcd_0 + g_0_z_0_xyyyzz_0[i] * pb_z + g_0_z_0_xyyyzz_1[i] * wp_z[i];

        g_0_zz_0_xyyzzz_0[i] = g_0_0_0_xyyzzz_0[i] * fi_ab_0 - g_0_0_0_xyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_z_0_xyyzz_1[i] * fi_abcd_0 + g_0_z_0_xyyzzz_0[i] * pb_z + g_0_z_0_xyyzzz_1[i] * wp_z[i];

        g_0_zz_0_xyzzzz_0[i] = g_0_0_0_xyzzzz_0[i] * fi_ab_0 - g_0_0_0_xyzzzz_1[i] * fti_ab_0 + 4.0 * g_0_z_0_xyzzz_1[i] * fi_abcd_0 + g_0_z_0_xyzzzz_0[i] * pb_z + g_0_z_0_xyzzzz_1[i] * wp_z[i];

        g_0_zz_0_xzzzzz_0[i] = g_0_0_0_xzzzzz_0[i] * fi_ab_0 - g_0_0_0_xzzzzz_1[i] * fti_ab_0 + 5.0 * g_0_z_0_xzzzz_1[i] * fi_abcd_0 + g_0_z_0_xzzzzz_0[i] * pb_z + g_0_z_0_xzzzzz_1[i] * wp_z[i];

        g_0_zz_0_yyyyyy_0[i] = g_0_0_0_yyyyyy_0[i] * fi_ab_0 - g_0_0_0_yyyyyy_1[i] * fti_ab_0 + g_0_z_0_yyyyyy_0[i] * pb_z + g_0_z_0_yyyyyy_1[i] * wp_z[i];

        g_0_zz_0_yyyyyz_0[i] = g_0_0_0_yyyyyz_0[i] * fi_ab_0 - g_0_0_0_yyyyyz_1[i] * fti_ab_0 + g_0_z_0_yyyyy_1[i] * fi_abcd_0 + g_0_z_0_yyyyyz_0[i] * pb_z + g_0_z_0_yyyyyz_1[i] * wp_z[i];

        g_0_zz_0_yyyyzz_0[i] = g_0_0_0_yyyyzz_0[i] * fi_ab_0 - g_0_0_0_yyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_z_0_yyyyz_1[i] * fi_abcd_0 + g_0_z_0_yyyyzz_0[i] * pb_z + g_0_z_0_yyyyzz_1[i] * wp_z[i];

        g_0_zz_0_yyyzzz_0[i] = g_0_0_0_yyyzzz_0[i] * fi_ab_0 - g_0_0_0_yyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_z_0_yyyzz_1[i] * fi_abcd_0 + g_0_z_0_yyyzzz_0[i] * pb_z + g_0_z_0_yyyzzz_1[i] * wp_z[i];

        g_0_zz_0_yyzzzz_0[i] = g_0_0_0_yyzzzz_0[i] * fi_ab_0 - g_0_0_0_yyzzzz_1[i] * fti_ab_0 + 4.0 * g_0_z_0_yyzzz_1[i] * fi_abcd_0 + g_0_z_0_yyzzzz_0[i] * pb_z + g_0_z_0_yyzzzz_1[i] * wp_z[i];

        g_0_zz_0_yzzzzz_0[i] = g_0_0_0_yzzzzz_0[i] * fi_ab_0 - g_0_0_0_yzzzzz_1[i] * fti_ab_0 + 5.0 * g_0_z_0_yzzzz_1[i] * fi_abcd_0 + g_0_z_0_yzzzzz_0[i] * pb_z + g_0_z_0_yzzzzz_1[i] * wp_z[i];

        g_0_zz_0_zzzzzz_0[i] = g_0_0_0_zzzzzz_0[i] * fi_ab_0 - g_0_0_0_zzzzzz_1[i] * fti_ab_0 + 6.0 * g_0_z_0_zzzzz_1[i] * fi_abcd_0 + g_0_z_0_zzzzzz_0[i] * pb_z + g_0_z_0_zzzzzz_1[i] * wp_z[i];
    }
}

} // erirec namespace


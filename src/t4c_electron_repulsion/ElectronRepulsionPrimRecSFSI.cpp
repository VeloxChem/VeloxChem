#include "ElectronRepulsionPrimRecSFSI.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_sfsi(CSimdArray<double>& pbuffer,
                                  const size_t idx_eri_0_sfsi,
                                  size_t idx_eri_0_spsi,
                                  size_t idx_eri_1_spsi,
                                  size_t idx_eri_1_sdsh,
                                  size_t idx_eri_0_sdsi,
                                  size_t idx_eri_1_sdsi,
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

    /// Set up components of auxilary buffer : SDSH

    auto g_0_xx_0_xxxxx_1 = pbuffer.data(idx_eri_1_sdsh);

    auto g_0_xx_0_xxxxy_1 = pbuffer.data(idx_eri_1_sdsh + 1);

    auto g_0_xx_0_xxxxz_1 = pbuffer.data(idx_eri_1_sdsh + 2);

    auto g_0_xx_0_xxxyy_1 = pbuffer.data(idx_eri_1_sdsh + 3);

    auto g_0_xx_0_xxxyz_1 = pbuffer.data(idx_eri_1_sdsh + 4);

    auto g_0_xx_0_xxxzz_1 = pbuffer.data(idx_eri_1_sdsh + 5);

    auto g_0_xx_0_xxyyy_1 = pbuffer.data(idx_eri_1_sdsh + 6);

    auto g_0_xx_0_xxyyz_1 = pbuffer.data(idx_eri_1_sdsh + 7);

    auto g_0_xx_0_xxyzz_1 = pbuffer.data(idx_eri_1_sdsh + 8);

    auto g_0_xx_0_xxzzz_1 = pbuffer.data(idx_eri_1_sdsh + 9);

    auto g_0_xx_0_xyyyy_1 = pbuffer.data(idx_eri_1_sdsh + 10);

    auto g_0_xx_0_xyyyz_1 = pbuffer.data(idx_eri_1_sdsh + 11);

    auto g_0_xx_0_xyyzz_1 = pbuffer.data(idx_eri_1_sdsh + 12);

    auto g_0_xx_0_xyzzz_1 = pbuffer.data(idx_eri_1_sdsh + 13);

    auto g_0_xx_0_xzzzz_1 = pbuffer.data(idx_eri_1_sdsh + 14);

    auto g_0_xx_0_yyyyy_1 = pbuffer.data(idx_eri_1_sdsh + 15);

    auto g_0_xx_0_yyyyz_1 = pbuffer.data(idx_eri_1_sdsh + 16);

    auto g_0_xx_0_yyyzz_1 = pbuffer.data(idx_eri_1_sdsh + 17);

    auto g_0_xx_0_yyzzz_1 = pbuffer.data(idx_eri_1_sdsh + 18);

    auto g_0_xx_0_yzzzz_1 = pbuffer.data(idx_eri_1_sdsh + 19);

    auto g_0_xx_0_zzzzz_1 = pbuffer.data(idx_eri_1_sdsh + 20);

    auto g_0_yy_0_xxxxx_1 = pbuffer.data(idx_eri_1_sdsh + 63);

    auto g_0_yy_0_xxxxy_1 = pbuffer.data(idx_eri_1_sdsh + 64);

    auto g_0_yy_0_xxxxz_1 = pbuffer.data(idx_eri_1_sdsh + 65);

    auto g_0_yy_0_xxxyy_1 = pbuffer.data(idx_eri_1_sdsh + 66);

    auto g_0_yy_0_xxxyz_1 = pbuffer.data(idx_eri_1_sdsh + 67);

    auto g_0_yy_0_xxxzz_1 = pbuffer.data(idx_eri_1_sdsh + 68);

    auto g_0_yy_0_xxyyy_1 = pbuffer.data(idx_eri_1_sdsh + 69);

    auto g_0_yy_0_xxyyz_1 = pbuffer.data(idx_eri_1_sdsh + 70);

    auto g_0_yy_0_xxyzz_1 = pbuffer.data(idx_eri_1_sdsh + 71);

    auto g_0_yy_0_xxzzz_1 = pbuffer.data(idx_eri_1_sdsh + 72);

    auto g_0_yy_0_xyyyy_1 = pbuffer.data(idx_eri_1_sdsh + 73);

    auto g_0_yy_0_xyyyz_1 = pbuffer.data(idx_eri_1_sdsh + 74);

    auto g_0_yy_0_xyyzz_1 = pbuffer.data(idx_eri_1_sdsh + 75);

    auto g_0_yy_0_xyzzz_1 = pbuffer.data(idx_eri_1_sdsh + 76);

    auto g_0_yy_0_xzzzz_1 = pbuffer.data(idx_eri_1_sdsh + 77);

    auto g_0_yy_0_yyyyy_1 = pbuffer.data(idx_eri_1_sdsh + 78);

    auto g_0_yy_0_yyyyz_1 = pbuffer.data(idx_eri_1_sdsh + 79);

    auto g_0_yy_0_yyyzz_1 = pbuffer.data(idx_eri_1_sdsh + 80);

    auto g_0_yy_0_yyzzz_1 = pbuffer.data(idx_eri_1_sdsh + 81);

    auto g_0_yy_0_yzzzz_1 = pbuffer.data(idx_eri_1_sdsh + 82);

    auto g_0_yy_0_zzzzz_1 = pbuffer.data(idx_eri_1_sdsh + 83);

    auto g_0_yz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sdsh + 88);

    auto g_0_yz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sdsh + 91);

    auto g_0_yz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sdsh + 92);

    auto g_0_yz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sdsh + 95);

    auto g_0_yz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sdsh + 96);

    auto g_0_yz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sdsh + 97);

    auto g_0_yz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sdsh + 100);

    auto g_0_yz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sdsh + 101);

    auto g_0_yz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sdsh + 102);

    auto g_0_yz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sdsh + 103);

    auto g_0_zz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sdsh + 105);

    auto g_0_zz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sdsh + 106);

    auto g_0_zz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sdsh + 107);

    auto g_0_zz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sdsh + 108);

    auto g_0_zz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sdsh + 109);

    auto g_0_zz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sdsh + 110);

    auto g_0_zz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sdsh + 111);

    auto g_0_zz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sdsh + 112);

    auto g_0_zz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sdsh + 113);

    auto g_0_zz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sdsh + 114);

    auto g_0_zz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sdsh + 115);

    auto g_0_zz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sdsh + 116);

    auto g_0_zz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sdsh + 117);

    auto g_0_zz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sdsh + 118);

    auto g_0_zz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sdsh + 119);

    auto g_0_zz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sdsh + 120);

    auto g_0_zz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sdsh + 121);

    auto g_0_zz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sdsh + 122);

    auto g_0_zz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sdsh + 123);

    auto g_0_zz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sdsh + 124);

    auto g_0_zz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sdsh + 125);

    /// Set up components of auxilary buffer : SDSI

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

    auto g_0_xy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sdsi + 29);

    auto g_0_xy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sdsi + 31);

    auto g_0_xy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sdsi + 34);

    auto g_0_xy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sdsi + 38);

    auto g_0_xy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sdsi + 43);

    auto g_0_xz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sdsi + 56);

    auto g_0_xz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sdsi + 58);

    auto g_0_xz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sdsi + 61);

    auto g_0_xz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sdsi + 65);

    auto g_0_xz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 70);

    auto g_0_xz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 76);

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

    auto g_0_yz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sdsi + 116);

    auto g_0_yz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sdsi + 119);

    auto g_0_yz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sdsi + 120);

    auto g_0_yz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sdsi + 123);

    auto g_0_yz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sdsi + 124);

    auto g_0_yz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sdsi + 125);

    auto g_0_yz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sdsi + 128);

    auto g_0_yz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sdsi + 129);

    auto g_0_yz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sdsi + 130);

    auto g_0_yz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 131);

    auto g_0_yz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sdsi + 133);

    auto g_0_yz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sdsi + 134);

    auto g_0_yz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sdsi + 135);

    auto g_0_yz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sdsi + 136);

    auto g_0_yz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 137);

    auto g_0_yz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 138);

    auto g_0_yz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sdsi + 139);

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

    /// Set up components of auxilary buffer : SDSI

    auto g_0_xx_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sdsi);

    auto g_0_xx_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sdsi + 1);

    auto g_0_xx_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sdsi + 2);

    auto g_0_xx_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sdsi + 3);

    auto g_0_xx_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sdsi + 4);

    auto g_0_xx_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sdsi + 5);

    auto g_0_xx_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sdsi + 6);

    auto g_0_xx_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sdsi + 7);

    auto g_0_xx_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sdsi + 8);

    auto g_0_xx_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sdsi + 9);

    auto g_0_xx_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sdsi + 10);

    auto g_0_xx_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sdsi + 11);

    auto g_0_xx_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sdsi + 12);

    auto g_0_xx_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sdsi + 13);

    auto g_0_xx_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 14);

    auto g_0_xx_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sdsi + 15);

    auto g_0_xx_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sdsi + 16);

    auto g_0_xx_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sdsi + 17);

    auto g_0_xx_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sdsi + 18);

    auto g_0_xx_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 19);

    auto g_0_xx_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 20);

    auto g_0_xx_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sdsi + 21);

    auto g_0_xx_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sdsi + 22);

    auto g_0_xx_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sdsi + 23);

    auto g_0_xx_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sdsi + 24);

    auto g_0_xx_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 25);

    auto g_0_xx_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 26);

    auto g_0_xx_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 27);

    auto g_0_xy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sdsi + 29);

    auto g_0_xy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sdsi + 31);

    auto g_0_xy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sdsi + 34);

    auto g_0_xy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sdsi + 38);

    auto g_0_xy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sdsi + 43);

    auto g_0_xz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sdsi + 56);

    auto g_0_xz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sdsi + 58);

    auto g_0_xz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sdsi + 61);

    auto g_0_xz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sdsi + 65);

    auto g_0_xz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 70);

    auto g_0_xz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 76);

    auto g_0_yy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sdsi + 84);

    auto g_0_yy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sdsi + 85);

    auto g_0_yy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sdsi + 86);

    auto g_0_yy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sdsi + 87);

    auto g_0_yy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sdsi + 88);

    auto g_0_yy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sdsi + 89);

    auto g_0_yy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sdsi + 90);

    auto g_0_yy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sdsi + 91);

    auto g_0_yy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sdsi + 92);

    auto g_0_yy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sdsi + 93);

    auto g_0_yy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sdsi + 94);

    auto g_0_yy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sdsi + 95);

    auto g_0_yy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sdsi + 96);

    auto g_0_yy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sdsi + 97);

    auto g_0_yy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 98);

    auto g_0_yy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sdsi + 99);

    auto g_0_yy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sdsi + 100);

    auto g_0_yy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sdsi + 101);

    auto g_0_yy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sdsi + 102);

    auto g_0_yy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 103);

    auto g_0_yy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 104);

    auto g_0_yy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sdsi + 105);

    auto g_0_yy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sdsi + 106);

    auto g_0_yy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sdsi + 107);

    auto g_0_yy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sdsi + 108);

    auto g_0_yy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 109);

    auto g_0_yy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 110);

    auto g_0_yy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 111);

    auto g_0_yz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sdsi + 116);

    auto g_0_yz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sdsi + 119);

    auto g_0_yz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sdsi + 120);

    auto g_0_yz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sdsi + 123);

    auto g_0_yz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sdsi + 124);

    auto g_0_yz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sdsi + 125);

    auto g_0_yz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sdsi + 128);

    auto g_0_yz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sdsi + 129);

    auto g_0_yz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sdsi + 130);

    auto g_0_yz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 131);

    auto g_0_yz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sdsi + 133);

    auto g_0_yz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sdsi + 134);

    auto g_0_yz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sdsi + 135);

    auto g_0_yz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sdsi + 136);

    auto g_0_yz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 137);

    auto g_0_yz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 138);

    auto g_0_yz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 139);

    auto g_0_zz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sdsi + 140);

    auto g_0_zz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sdsi + 141);

    auto g_0_zz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sdsi + 142);

    auto g_0_zz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sdsi + 143);

    auto g_0_zz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sdsi + 144);

    auto g_0_zz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sdsi + 145);

    auto g_0_zz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sdsi + 146);

    auto g_0_zz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sdsi + 147);

    auto g_0_zz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sdsi + 148);

    auto g_0_zz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sdsi + 149);

    auto g_0_zz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sdsi + 150);

    auto g_0_zz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sdsi + 151);

    auto g_0_zz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sdsi + 152);

    auto g_0_zz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sdsi + 153);

    auto g_0_zz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 154);

    auto g_0_zz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sdsi + 155);

    auto g_0_zz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sdsi + 156);

    auto g_0_zz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sdsi + 157);

    auto g_0_zz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sdsi + 158);

    auto g_0_zz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 159);

    auto g_0_zz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 160);

    auto g_0_zz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sdsi + 161);

    auto g_0_zz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sdsi + 162);

    auto g_0_zz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sdsi + 163);

    auto g_0_zz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sdsi + 164);

    auto g_0_zz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 165);

    auto g_0_zz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 166);

    auto g_0_zz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sdsi + 167);

    /// Set up 0-28 components of targeted buffer : SFSI

    auto g_0_xxx_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sfsi);

    auto g_0_xxx_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sfsi + 1);

    auto g_0_xxx_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sfsi + 2);

    auto g_0_xxx_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sfsi + 3);

    auto g_0_xxx_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sfsi + 4);

    auto g_0_xxx_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sfsi + 5);

    auto g_0_xxx_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sfsi + 6);

    auto g_0_xxx_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sfsi + 7);

    auto g_0_xxx_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sfsi + 8);

    auto g_0_xxx_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sfsi + 9);

    auto g_0_xxx_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 10);

    auto g_0_xxx_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 11);

    auto g_0_xxx_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 12);

    auto g_0_xxx_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 13);

    auto g_0_xxx_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 14);

    auto g_0_xxx_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 15);

    auto g_0_xxx_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 16);

    auto g_0_xxx_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 17);

    auto g_0_xxx_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 18);

    auto g_0_xxx_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 19);

    auto g_0_xxx_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 20);

    auto g_0_xxx_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 21);

    auto g_0_xxx_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 22);

    auto g_0_xxx_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 23);

    auto g_0_xxx_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 24);

    auto g_0_xxx_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 25);

    auto g_0_xxx_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 26);

    auto g_0_xxx_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 27);

    #pragma omp simd aligned(g_0_x_0_xxxxxx_0, g_0_x_0_xxxxxx_1, g_0_x_0_xxxxxy_0, g_0_x_0_xxxxxy_1, g_0_x_0_xxxxxz_0, g_0_x_0_xxxxxz_1, g_0_x_0_xxxxyy_0, g_0_x_0_xxxxyy_1, g_0_x_0_xxxxyz_0, g_0_x_0_xxxxyz_1, g_0_x_0_xxxxzz_0, g_0_x_0_xxxxzz_1, g_0_x_0_xxxyyy_0, g_0_x_0_xxxyyy_1, g_0_x_0_xxxyyz_0, g_0_x_0_xxxyyz_1, g_0_x_0_xxxyzz_0, g_0_x_0_xxxyzz_1, g_0_x_0_xxxzzz_0, g_0_x_0_xxxzzz_1, g_0_x_0_xxyyyy_0, g_0_x_0_xxyyyy_1, g_0_x_0_xxyyyz_0, g_0_x_0_xxyyyz_1, g_0_x_0_xxyyzz_0, g_0_x_0_xxyyzz_1, g_0_x_0_xxyzzz_0, g_0_x_0_xxyzzz_1, g_0_x_0_xxzzzz_0, g_0_x_0_xxzzzz_1, g_0_x_0_xyyyyy_0, g_0_x_0_xyyyyy_1, g_0_x_0_xyyyyz_0, g_0_x_0_xyyyyz_1, g_0_x_0_xyyyzz_0, g_0_x_0_xyyyzz_1, g_0_x_0_xyyzzz_0, g_0_x_0_xyyzzz_1, g_0_x_0_xyzzzz_0, g_0_x_0_xyzzzz_1, g_0_x_0_xzzzzz_0, g_0_x_0_xzzzzz_1, g_0_x_0_yyyyyy_0, g_0_x_0_yyyyyy_1, g_0_x_0_yyyyyz_0, g_0_x_0_yyyyyz_1, g_0_x_0_yyyyzz_0, g_0_x_0_yyyyzz_1, g_0_x_0_yyyzzz_0, g_0_x_0_yyyzzz_1, g_0_x_0_yyzzzz_0, g_0_x_0_yyzzzz_1, g_0_x_0_yzzzzz_0, g_0_x_0_yzzzzz_1, g_0_x_0_zzzzzz_0, g_0_x_0_zzzzzz_1, g_0_xx_0_xxxxx_1, g_0_xx_0_xxxxxx_0, g_0_xx_0_xxxxxx_1, g_0_xx_0_xxxxxy_0, g_0_xx_0_xxxxxy_1, g_0_xx_0_xxxxxz_0, g_0_xx_0_xxxxxz_1, g_0_xx_0_xxxxy_1, g_0_xx_0_xxxxyy_0, g_0_xx_0_xxxxyy_1, g_0_xx_0_xxxxyz_0, g_0_xx_0_xxxxyz_1, g_0_xx_0_xxxxz_1, g_0_xx_0_xxxxzz_0, g_0_xx_0_xxxxzz_1, g_0_xx_0_xxxyy_1, g_0_xx_0_xxxyyy_0, g_0_xx_0_xxxyyy_1, g_0_xx_0_xxxyyz_0, g_0_xx_0_xxxyyz_1, g_0_xx_0_xxxyz_1, g_0_xx_0_xxxyzz_0, g_0_xx_0_xxxyzz_1, g_0_xx_0_xxxzz_1, g_0_xx_0_xxxzzz_0, g_0_xx_0_xxxzzz_1, g_0_xx_0_xxyyy_1, g_0_xx_0_xxyyyy_0, g_0_xx_0_xxyyyy_1, g_0_xx_0_xxyyyz_0, g_0_xx_0_xxyyyz_1, g_0_xx_0_xxyyz_1, g_0_xx_0_xxyyzz_0, g_0_xx_0_xxyyzz_1, g_0_xx_0_xxyzz_1, g_0_xx_0_xxyzzz_0, g_0_xx_0_xxyzzz_1, g_0_xx_0_xxzzz_1, g_0_xx_0_xxzzzz_0, g_0_xx_0_xxzzzz_1, g_0_xx_0_xyyyy_1, g_0_xx_0_xyyyyy_0, g_0_xx_0_xyyyyy_1, g_0_xx_0_xyyyyz_0, g_0_xx_0_xyyyyz_1, g_0_xx_0_xyyyz_1, g_0_xx_0_xyyyzz_0, g_0_xx_0_xyyyzz_1, g_0_xx_0_xyyzz_1, g_0_xx_0_xyyzzz_0, g_0_xx_0_xyyzzz_1, g_0_xx_0_xyzzz_1, g_0_xx_0_xyzzzz_0, g_0_xx_0_xyzzzz_1, g_0_xx_0_xzzzz_1, g_0_xx_0_xzzzzz_0, g_0_xx_0_xzzzzz_1, g_0_xx_0_yyyyy_1, g_0_xx_0_yyyyyy_0, g_0_xx_0_yyyyyy_1, g_0_xx_0_yyyyyz_0, g_0_xx_0_yyyyyz_1, g_0_xx_0_yyyyz_1, g_0_xx_0_yyyyzz_0, g_0_xx_0_yyyyzz_1, g_0_xx_0_yyyzz_1, g_0_xx_0_yyyzzz_0, g_0_xx_0_yyyzzz_1, g_0_xx_0_yyzzz_1, g_0_xx_0_yyzzzz_0, g_0_xx_0_yyzzzz_1, g_0_xx_0_yzzzz_1, g_0_xx_0_yzzzzz_0, g_0_xx_0_yzzzzz_1, g_0_xx_0_zzzzz_1, g_0_xx_0_zzzzzz_0, g_0_xx_0_zzzzzz_1, g_0_xxx_0_xxxxxx_0, g_0_xxx_0_xxxxxy_0, g_0_xxx_0_xxxxxz_0, g_0_xxx_0_xxxxyy_0, g_0_xxx_0_xxxxyz_0, g_0_xxx_0_xxxxzz_0, g_0_xxx_0_xxxyyy_0, g_0_xxx_0_xxxyyz_0, g_0_xxx_0_xxxyzz_0, g_0_xxx_0_xxxzzz_0, g_0_xxx_0_xxyyyy_0, g_0_xxx_0_xxyyyz_0, g_0_xxx_0_xxyyzz_0, g_0_xxx_0_xxyzzz_0, g_0_xxx_0_xxzzzz_0, g_0_xxx_0_xyyyyy_0, g_0_xxx_0_xyyyyz_0, g_0_xxx_0_xyyyzz_0, g_0_xxx_0_xyyzzz_0, g_0_xxx_0_xyzzzz_0, g_0_xxx_0_xzzzzz_0, g_0_xxx_0_yyyyyy_0, g_0_xxx_0_yyyyyz_0, g_0_xxx_0_yyyyzz_0, g_0_xxx_0_yyyzzz_0, g_0_xxx_0_yyzzzz_0, g_0_xxx_0_yzzzzz_0, g_0_xxx_0_zzzzzz_0, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxx_0_xxxxxx_0[i] = 2.0 * g_0_x_0_xxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxxxx_1[i] * fti_ab_0 + 6.0 * g_0_xx_0_xxxxx_1[i] * fi_abcd_0 + g_0_xx_0_xxxxxx_0[i] * pb_x + g_0_xx_0_xxxxxx_1[i] * wp_x[i];

        g_0_xxx_0_xxxxxy_0[i] = 2.0 * g_0_x_0_xxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxxxy_1[i] * fti_ab_0 + 5.0 * g_0_xx_0_xxxxy_1[i] * fi_abcd_0 + g_0_xx_0_xxxxxy_0[i] * pb_x + g_0_xx_0_xxxxxy_1[i] * wp_x[i];

        g_0_xxx_0_xxxxxz_0[i] = 2.0 * g_0_x_0_xxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxxxz_1[i] * fti_ab_0 + 5.0 * g_0_xx_0_xxxxz_1[i] * fi_abcd_0 + g_0_xx_0_xxxxxz_0[i] * pb_x + g_0_xx_0_xxxxxz_1[i] * wp_x[i];

        g_0_xxx_0_xxxxyy_0[i] = 2.0 * g_0_x_0_xxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxxyy_1[i] * fti_ab_0 + 4.0 * g_0_xx_0_xxxyy_1[i] * fi_abcd_0 + g_0_xx_0_xxxxyy_0[i] * pb_x + g_0_xx_0_xxxxyy_1[i] * wp_x[i];

        g_0_xxx_0_xxxxyz_0[i] = 2.0 * g_0_x_0_xxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxxyz_1[i] * fti_ab_0 + 4.0 * g_0_xx_0_xxxyz_1[i] * fi_abcd_0 + g_0_xx_0_xxxxyz_0[i] * pb_x + g_0_xx_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxx_0_xxxxzz_0[i] = 2.0 * g_0_x_0_xxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxxzz_1[i] * fti_ab_0 + 4.0 * g_0_xx_0_xxxzz_1[i] * fi_abcd_0 + g_0_xx_0_xxxxzz_0[i] * pb_x + g_0_xx_0_xxxxzz_1[i] * wp_x[i];

        g_0_xxx_0_xxxyyy_0[i] = 2.0 * g_0_x_0_xxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxyyy_1[i] * fti_ab_0 + 3.0 * g_0_xx_0_xxyyy_1[i] * fi_abcd_0 + g_0_xx_0_xxxyyy_0[i] * pb_x + g_0_xx_0_xxxyyy_1[i] * wp_x[i];

        g_0_xxx_0_xxxyyz_0[i] = 2.0 * g_0_x_0_xxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxyyz_1[i] * fti_ab_0 + 3.0 * g_0_xx_0_xxyyz_1[i] * fi_abcd_0 + g_0_xx_0_xxxyyz_0[i] * pb_x + g_0_xx_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxx_0_xxxyzz_0[i] = 2.0 * g_0_x_0_xxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxyzz_1[i] * fti_ab_0 + 3.0 * g_0_xx_0_xxyzz_1[i] * fi_abcd_0 + g_0_xx_0_xxxyzz_0[i] * pb_x + g_0_xx_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxx_0_xxxzzz_0[i] = 2.0 * g_0_x_0_xxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxzzz_1[i] * fti_ab_0 + 3.0 * g_0_xx_0_xxzzz_1[i] * fi_abcd_0 + g_0_xx_0_xxxzzz_0[i] * pb_x + g_0_xx_0_xxxzzz_1[i] * wp_x[i];

        g_0_xxx_0_xxyyyy_0[i] = 2.0 * g_0_x_0_xxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxyyyy_1[i] * fti_ab_0 + 2.0 * g_0_xx_0_xyyyy_1[i] * fi_abcd_0 + g_0_xx_0_xxyyyy_0[i] * pb_x + g_0_xx_0_xxyyyy_1[i] * wp_x[i];

        g_0_xxx_0_xxyyyz_0[i] = 2.0 * g_0_x_0_xxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xx_0_xyyyz_1[i] * fi_abcd_0 + g_0_xx_0_xxyyyz_0[i] * pb_x + g_0_xx_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxx_0_xxyyzz_0[i] = 2.0 * g_0_x_0_xxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xx_0_xyyzz_1[i] * fi_abcd_0 + g_0_xx_0_xxyyzz_0[i] * pb_x + g_0_xx_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxx_0_xxyzzz_0[i] = 2.0 * g_0_x_0_xxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xx_0_xyzzz_1[i] * fi_abcd_0 + g_0_xx_0_xxyzzz_0[i] * pb_x + g_0_xx_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxx_0_xxzzzz_0[i] = 2.0 * g_0_x_0_xxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xx_0_xzzzz_1[i] * fi_abcd_0 + g_0_xx_0_xxzzzz_0[i] * pb_x + g_0_xx_0_xxzzzz_1[i] * wp_x[i];

        g_0_xxx_0_xyyyyy_0[i] = 2.0 * g_0_x_0_xyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xyyyyy_1[i] * fti_ab_0 + g_0_xx_0_yyyyy_1[i] * fi_abcd_0 + g_0_xx_0_xyyyyy_0[i] * pb_x + g_0_xx_0_xyyyyy_1[i] * wp_x[i];

        g_0_xxx_0_xyyyyz_0[i] = 2.0 * g_0_x_0_xyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xyyyyz_1[i] * fti_ab_0 + g_0_xx_0_yyyyz_1[i] * fi_abcd_0 + g_0_xx_0_xyyyyz_0[i] * pb_x + g_0_xx_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxx_0_xyyyzz_0[i] = 2.0 * g_0_x_0_xyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xyyyzz_1[i] * fti_ab_0 + g_0_xx_0_yyyzz_1[i] * fi_abcd_0 + g_0_xx_0_xyyyzz_0[i] * pb_x + g_0_xx_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxx_0_xyyzzz_0[i] = 2.0 * g_0_x_0_xyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xyyzzz_1[i] * fti_ab_0 + g_0_xx_0_yyzzz_1[i] * fi_abcd_0 + g_0_xx_0_xyyzzz_0[i] * pb_x + g_0_xx_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxx_0_xyzzzz_0[i] = 2.0 * g_0_x_0_xyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xyzzzz_1[i] * fti_ab_0 + g_0_xx_0_yzzzz_1[i] * fi_abcd_0 + g_0_xx_0_xyzzzz_0[i] * pb_x + g_0_xx_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxx_0_xzzzzz_0[i] = 2.0 * g_0_x_0_xzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xzzzzz_1[i] * fti_ab_0 + g_0_xx_0_zzzzz_1[i] * fi_abcd_0 + g_0_xx_0_xzzzzz_0[i] * pb_x + g_0_xx_0_xzzzzz_1[i] * wp_x[i];

        g_0_xxx_0_yyyyyy_0[i] = 2.0 * g_0_x_0_yyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_x_0_yyyyyy_1[i] * fti_ab_0 + g_0_xx_0_yyyyyy_0[i] * pb_x + g_0_xx_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxx_0_yyyyyz_0[i] = 2.0 * g_0_x_0_yyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_yyyyyz_1[i] * fti_ab_0 + g_0_xx_0_yyyyyz_0[i] * pb_x + g_0_xx_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxx_0_yyyyzz_0[i] = 2.0 * g_0_x_0_yyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_yyyyzz_1[i] * fti_ab_0 + g_0_xx_0_yyyyzz_0[i] * pb_x + g_0_xx_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxx_0_yyyzzz_0[i] = 2.0 * g_0_x_0_yyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_yyyzzz_1[i] * fti_ab_0 + g_0_xx_0_yyyzzz_0[i] * pb_x + g_0_xx_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxx_0_yyzzzz_0[i] = 2.0 * g_0_x_0_yyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_yyzzzz_1[i] * fti_ab_0 + g_0_xx_0_yyzzzz_0[i] * pb_x + g_0_xx_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxx_0_yzzzzz_0[i] = 2.0 * g_0_x_0_yzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_yzzzzz_1[i] * fti_ab_0 + g_0_xx_0_yzzzzz_0[i] * pb_x + g_0_xx_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxx_0_zzzzzz_0[i] = 2.0 * g_0_x_0_zzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_zzzzzz_1[i] * fti_ab_0 + g_0_xx_0_zzzzzz_0[i] * pb_x + g_0_xx_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 28-56 components of targeted buffer : SFSI

    auto g_0_xxy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sfsi + 28);

    auto g_0_xxy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sfsi + 29);

    auto g_0_xxy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sfsi + 30);

    auto g_0_xxy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sfsi + 31);

    auto g_0_xxy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sfsi + 32);

    auto g_0_xxy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sfsi + 33);

    auto g_0_xxy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sfsi + 34);

    auto g_0_xxy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sfsi + 35);

    auto g_0_xxy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sfsi + 36);

    auto g_0_xxy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sfsi + 37);

    auto g_0_xxy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 38);

    auto g_0_xxy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 39);

    auto g_0_xxy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 40);

    auto g_0_xxy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 41);

    auto g_0_xxy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 42);

    auto g_0_xxy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 43);

    auto g_0_xxy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 44);

    auto g_0_xxy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 45);

    auto g_0_xxy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 46);

    auto g_0_xxy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 47);

    auto g_0_xxy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 48);

    auto g_0_xxy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 49);

    auto g_0_xxy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 50);

    auto g_0_xxy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 51);

    auto g_0_xxy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 52);

    auto g_0_xxy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 53);

    auto g_0_xxy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 54);

    auto g_0_xxy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 55);

    #pragma omp simd aligned(g_0_xx_0_xxxxx_1, g_0_xx_0_xxxxxx_0, g_0_xx_0_xxxxxx_1, g_0_xx_0_xxxxxy_0, g_0_xx_0_xxxxxy_1, g_0_xx_0_xxxxxz_0, g_0_xx_0_xxxxxz_1, g_0_xx_0_xxxxy_1, g_0_xx_0_xxxxyy_0, g_0_xx_0_xxxxyy_1, g_0_xx_0_xxxxyz_0, g_0_xx_0_xxxxyz_1, g_0_xx_0_xxxxz_1, g_0_xx_0_xxxxzz_0, g_0_xx_0_xxxxzz_1, g_0_xx_0_xxxyy_1, g_0_xx_0_xxxyyy_0, g_0_xx_0_xxxyyy_1, g_0_xx_0_xxxyyz_0, g_0_xx_0_xxxyyz_1, g_0_xx_0_xxxyz_1, g_0_xx_0_xxxyzz_0, g_0_xx_0_xxxyzz_1, g_0_xx_0_xxxzz_1, g_0_xx_0_xxxzzz_0, g_0_xx_0_xxxzzz_1, g_0_xx_0_xxyyy_1, g_0_xx_0_xxyyyy_0, g_0_xx_0_xxyyyy_1, g_0_xx_0_xxyyyz_0, g_0_xx_0_xxyyyz_1, g_0_xx_0_xxyyz_1, g_0_xx_0_xxyyzz_0, g_0_xx_0_xxyyzz_1, g_0_xx_0_xxyzz_1, g_0_xx_0_xxyzzz_0, g_0_xx_0_xxyzzz_1, g_0_xx_0_xxzzz_1, g_0_xx_0_xxzzzz_0, g_0_xx_0_xxzzzz_1, g_0_xx_0_xyyyy_1, g_0_xx_0_xyyyyy_0, g_0_xx_0_xyyyyy_1, g_0_xx_0_xyyyyz_0, g_0_xx_0_xyyyyz_1, g_0_xx_0_xyyyz_1, g_0_xx_0_xyyyzz_0, g_0_xx_0_xyyyzz_1, g_0_xx_0_xyyzz_1, g_0_xx_0_xyyzzz_0, g_0_xx_0_xyyzzz_1, g_0_xx_0_xyzzz_1, g_0_xx_0_xyzzzz_0, g_0_xx_0_xyzzzz_1, g_0_xx_0_xzzzz_1, g_0_xx_0_xzzzzz_0, g_0_xx_0_xzzzzz_1, g_0_xx_0_yyyyy_1, g_0_xx_0_yyyyyy_0, g_0_xx_0_yyyyyy_1, g_0_xx_0_yyyyyz_0, g_0_xx_0_yyyyyz_1, g_0_xx_0_yyyyz_1, g_0_xx_0_yyyyzz_0, g_0_xx_0_yyyyzz_1, g_0_xx_0_yyyzz_1, g_0_xx_0_yyyzzz_0, g_0_xx_0_yyyzzz_1, g_0_xx_0_yyzzz_1, g_0_xx_0_yyzzzz_0, g_0_xx_0_yyzzzz_1, g_0_xx_0_yzzzz_1, g_0_xx_0_yzzzzz_0, g_0_xx_0_yzzzzz_1, g_0_xx_0_zzzzz_1, g_0_xx_0_zzzzzz_0, g_0_xx_0_zzzzzz_1, g_0_xxy_0_xxxxxx_0, g_0_xxy_0_xxxxxy_0, g_0_xxy_0_xxxxxz_0, g_0_xxy_0_xxxxyy_0, g_0_xxy_0_xxxxyz_0, g_0_xxy_0_xxxxzz_0, g_0_xxy_0_xxxyyy_0, g_0_xxy_0_xxxyyz_0, g_0_xxy_0_xxxyzz_0, g_0_xxy_0_xxxzzz_0, g_0_xxy_0_xxyyyy_0, g_0_xxy_0_xxyyyz_0, g_0_xxy_0_xxyyzz_0, g_0_xxy_0_xxyzzz_0, g_0_xxy_0_xxzzzz_0, g_0_xxy_0_xyyyyy_0, g_0_xxy_0_xyyyyz_0, g_0_xxy_0_xyyyzz_0, g_0_xxy_0_xyyzzz_0, g_0_xxy_0_xyzzzz_0, g_0_xxy_0_xzzzzz_0, g_0_xxy_0_yyyyyy_0, g_0_xxy_0_yyyyyz_0, g_0_xxy_0_yyyyzz_0, g_0_xxy_0_yyyzzz_0, g_0_xxy_0_yyzzzz_0, g_0_xxy_0_yzzzzz_0, g_0_xxy_0_zzzzzz_0, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxy_0_xxxxxx_0[i] = g_0_xx_0_xxxxxx_0[i] * pb_y + g_0_xx_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxy_0_xxxxxy_0[i] = g_0_xx_0_xxxxx_1[i] * fi_abcd_0 + g_0_xx_0_xxxxxy_0[i] * pb_y + g_0_xx_0_xxxxxy_1[i] * wp_y[i];

        g_0_xxy_0_xxxxxz_0[i] = g_0_xx_0_xxxxxz_0[i] * pb_y + g_0_xx_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxy_0_xxxxyy_0[i] = 2.0 * g_0_xx_0_xxxxy_1[i] * fi_abcd_0 + g_0_xx_0_xxxxyy_0[i] * pb_y + g_0_xx_0_xxxxyy_1[i] * wp_y[i];

        g_0_xxy_0_xxxxyz_0[i] = g_0_xx_0_xxxxz_1[i] * fi_abcd_0 + g_0_xx_0_xxxxyz_0[i] * pb_y + g_0_xx_0_xxxxyz_1[i] * wp_y[i];

        g_0_xxy_0_xxxxzz_0[i] = g_0_xx_0_xxxxzz_0[i] * pb_y + g_0_xx_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxy_0_xxxyyy_0[i] = 3.0 * g_0_xx_0_xxxyy_1[i] * fi_abcd_0 + g_0_xx_0_xxxyyy_0[i] * pb_y + g_0_xx_0_xxxyyy_1[i] * wp_y[i];

        g_0_xxy_0_xxxyyz_0[i] = 2.0 * g_0_xx_0_xxxyz_1[i] * fi_abcd_0 + g_0_xx_0_xxxyyz_0[i] * pb_y + g_0_xx_0_xxxyyz_1[i] * wp_y[i];

        g_0_xxy_0_xxxyzz_0[i] = g_0_xx_0_xxxzz_1[i] * fi_abcd_0 + g_0_xx_0_xxxyzz_0[i] * pb_y + g_0_xx_0_xxxyzz_1[i] * wp_y[i];

        g_0_xxy_0_xxxzzz_0[i] = g_0_xx_0_xxxzzz_0[i] * pb_y + g_0_xx_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxy_0_xxyyyy_0[i] = 4.0 * g_0_xx_0_xxyyy_1[i] * fi_abcd_0 + g_0_xx_0_xxyyyy_0[i] * pb_y + g_0_xx_0_xxyyyy_1[i] * wp_y[i];

        g_0_xxy_0_xxyyyz_0[i] = 3.0 * g_0_xx_0_xxyyz_1[i] * fi_abcd_0 + g_0_xx_0_xxyyyz_0[i] * pb_y + g_0_xx_0_xxyyyz_1[i] * wp_y[i];

        g_0_xxy_0_xxyyzz_0[i] = 2.0 * g_0_xx_0_xxyzz_1[i] * fi_abcd_0 + g_0_xx_0_xxyyzz_0[i] * pb_y + g_0_xx_0_xxyyzz_1[i] * wp_y[i];

        g_0_xxy_0_xxyzzz_0[i] = g_0_xx_0_xxzzz_1[i] * fi_abcd_0 + g_0_xx_0_xxyzzz_0[i] * pb_y + g_0_xx_0_xxyzzz_1[i] * wp_y[i];

        g_0_xxy_0_xxzzzz_0[i] = g_0_xx_0_xxzzzz_0[i] * pb_y + g_0_xx_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxy_0_xyyyyy_0[i] = 5.0 * g_0_xx_0_xyyyy_1[i] * fi_abcd_0 + g_0_xx_0_xyyyyy_0[i] * pb_y + g_0_xx_0_xyyyyy_1[i] * wp_y[i];

        g_0_xxy_0_xyyyyz_0[i] = 4.0 * g_0_xx_0_xyyyz_1[i] * fi_abcd_0 + g_0_xx_0_xyyyyz_0[i] * pb_y + g_0_xx_0_xyyyyz_1[i] * wp_y[i];

        g_0_xxy_0_xyyyzz_0[i] = 3.0 * g_0_xx_0_xyyzz_1[i] * fi_abcd_0 + g_0_xx_0_xyyyzz_0[i] * pb_y + g_0_xx_0_xyyyzz_1[i] * wp_y[i];

        g_0_xxy_0_xyyzzz_0[i] = 2.0 * g_0_xx_0_xyzzz_1[i] * fi_abcd_0 + g_0_xx_0_xyyzzz_0[i] * pb_y + g_0_xx_0_xyyzzz_1[i] * wp_y[i];

        g_0_xxy_0_xyzzzz_0[i] = g_0_xx_0_xzzzz_1[i] * fi_abcd_0 + g_0_xx_0_xyzzzz_0[i] * pb_y + g_0_xx_0_xyzzzz_1[i] * wp_y[i];

        g_0_xxy_0_xzzzzz_0[i] = g_0_xx_0_xzzzzz_0[i] * pb_y + g_0_xx_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxy_0_yyyyyy_0[i] = 6.0 * g_0_xx_0_yyyyy_1[i] * fi_abcd_0 + g_0_xx_0_yyyyyy_0[i] * pb_y + g_0_xx_0_yyyyyy_1[i] * wp_y[i];

        g_0_xxy_0_yyyyyz_0[i] = 5.0 * g_0_xx_0_yyyyz_1[i] * fi_abcd_0 + g_0_xx_0_yyyyyz_0[i] * pb_y + g_0_xx_0_yyyyyz_1[i] * wp_y[i];

        g_0_xxy_0_yyyyzz_0[i] = 4.0 * g_0_xx_0_yyyzz_1[i] * fi_abcd_0 + g_0_xx_0_yyyyzz_0[i] * pb_y + g_0_xx_0_yyyyzz_1[i] * wp_y[i];

        g_0_xxy_0_yyyzzz_0[i] = 3.0 * g_0_xx_0_yyzzz_1[i] * fi_abcd_0 + g_0_xx_0_yyyzzz_0[i] * pb_y + g_0_xx_0_yyyzzz_1[i] * wp_y[i];

        g_0_xxy_0_yyzzzz_0[i] = 2.0 * g_0_xx_0_yzzzz_1[i] * fi_abcd_0 + g_0_xx_0_yyzzzz_0[i] * pb_y + g_0_xx_0_yyzzzz_1[i] * wp_y[i];

        g_0_xxy_0_yzzzzz_0[i] = g_0_xx_0_zzzzz_1[i] * fi_abcd_0 + g_0_xx_0_yzzzzz_0[i] * pb_y + g_0_xx_0_yzzzzz_1[i] * wp_y[i];

        g_0_xxy_0_zzzzzz_0[i] = g_0_xx_0_zzzzzz_0[i] * pb_y + g_0_xx_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 56-84 components of targeted buffer : SFSI

    auto g_0_xxz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sfsi + 56);

    auto g_0_xxz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sfsi + 57);

    auto g_0_xxz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sfsi + 58);

    auto g_0_xxz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sfsi + 59);

    auto g_0_xxz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sfsi + 60);

    auto g_0_xxz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sfsi + 61);

    auto g_0_xxz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sfsi + 62);

    auto g_0_xxz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sfsi + 63);

    auto g_0_xxz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sfsi + 64);

    auto g_0_xxz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sfsi + 65);

    auto g_0_xxz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 66);

    auto g_0_xxz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 67);

    auto g_0_xxz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 68);

    auto g_0_xxz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 69);

    auto g_0_xxz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 70);

    auto g_0_xxz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 71);

    auto g_0_xxz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 72);

    auto g_0_xxz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 73);

    auto g_0_xxz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 74);

    auto g_0_xxz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 75);

    auto g_0_xxz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 76);

    auto g_0_xxz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 77);

    auto g_0_xxz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 78);

    auto g_0_xxz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 79);

    auto g_0_xxz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 80);

    auto g_0_xxz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 81);

    auto g_0_xxz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 82);

    auto g_0_xxz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 83);

    #pragma omp simd aligned(g_0_xx_0_xxxxx_1, g_0_xx_0_xxxxxx_0, g_0_xx_0_xxxxxx_1, g_0_xx_0_xxxxxy_0, g_0_xx_0_xxxxxy_1, g_0_xx_0_xxxxxz_0, g_0_xx_0_xxxxxz_1, g_0_xx_0_xxxxy_1, g_0_xx_0_xxxxyy_0, g_0_xx_0_xxxxyy_1, g_0_xx_0_xxxxyz_0, g_0_xx_0_xxxxyz_1, g_0_xx_0_xxxxz_1, g_0_xx_0_xxxxzz_0, g_0_xx_0_xxxxzz_1, g_0_xx_0_xxxyy_1, g_0_xx_0_xxxyyy_0, g_0_xx_0_xxxyyy_1, g_0_xx_0_xxxyyz_0, g_0_xx_0_xxxyyz_1, g_0_xx_0_xxxyz_1, g_0_xx_0_xxxyzz_0, g_0_xx_0_xxxyzz_1, g_0_xx_0_xxxzz_1, g_0_xx_0_xxxzzz_0, g_0_xx_0_xxxzzz_1, g_0_xx_0_xxyyy_1, g_0_xx_0_xxyyyy_0, g_0_xx_0_xxyyyy_1, g_0_xx_0_xxyyyz_0, g_0_xx_0_xxyyyz_1, g_0_xx_0_xxyyz_1, g_0_xx_0_xxyyzz_0, g_0_xx_0_xxyyzz_1, g_0_xx_0_xxyzz_1, g_0_xx_0_xxyzzz_0, g_0_xx_0_xxyzzz_1, g_0_xx_0_xxzzz_1, g_0_xx_0_xxzzzz_0, g_0_xx_0_xxzzzz_1, g_0_xx_0_xyyyy_1, g_0_xx_0_xyyyyy_0, g_0_xx_0_xyyyyy_1, g_0_xx_0_xyyyyz_0, g_0_xx_0_xyyyyz_1, g_0_xx_0_xyyyz_1, g_0_xx_0_xyyyzz_0, g_0_xx_0_xyyyzz_1, g_0_xx_0_xyyzz_1, g_0_xx_0_xyyzzz_0, g_0_xx_0_xyyzzz_1, g_0_xx_0_xyzzz_1, g_0_xx_0_xyzzzz_0, g_0_xx_0_xyzzzz_1, g_0_xx_0_xzzzz_1, g_0_xx_0_xzzzzz_0, g_0_xx_0_xzzzzz_1, g_0_xx_0_yyyyy_1, g_0_xx_0_yyyyyy_0, g_0_xx_0_yyyyyy_1, g_0_xx_0_yyyyyz_0, g_0_xx_0_yyyyyz_1, g_0_xx_0_yyyyz_1, g_0_xx_0_yyyyzz_0, g_0_xx_0_yyyyzz_1, g_0_xx_0_yyyzz_1, g_0_xx_0_yyyzzz_0, g_0_xx_0_yyyzzz_1, g_0_xx_0_yyzzz_1, g_0_xx_0_yyzzzz_0, g_0_xx_0_yyzzzz_1, g_0_xx_0_yzzzz_1, g_0_xx_0_yzzzzz_0, g_0_xx_0_yzzzzz_1, g_0_xx_0_zzzzz_1, g_0_xx_0_zzzzzz_0, g_0_xx_0_zzzzzz_1, g_0_xxz_0_xxxxxx_0, g_0_xxz_0_xxxxxy_0, g_0_xxz_0_xxxxxz_0, g_0_xxz_0_xxxxyy_0, g_0_xxz_0_xxxxyz_0, g_0_xxz_0_xxxxzz_0, g_0_xxz_0_xxxyyy_0, g_0_xxz_0_xxxyyz_0, g_0_xxz_0_xxxyzz_0, g_0_xxz_0_xxxzzz_0, g_0_xxz_0_xxyyyy_0, g_0_xxz_0_xxyyyz_0, g_0_xxz_0_xxyyzz_0, g_0_xxz_0_xxyzzz_0, g_0_xxz_0_xxzzzz_0, g_0_xxz_0_xyyyyy_0, g_0_xxz_0_xyyyyz_0, g_0_xxz_0_xyyyzz_0, g_0_xxz_0_xyyzzz_0, g_0_xxz_0_xyzzzz_0, g_0_xxz_0_xzzzzz_0, g_0_xxz_0_yyyyyy_0, g_0_xxz_0_yyyyyz_0, g_0_xxz_0_yyyyzz_0, g_0_xxz_0_yyyzzz_0, g_0_xxz_0_yyzzzz_0, g_0_xxz_0_yzzzzz_0, g_0_xxz_0_zzzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxz_0_xxxxxx_0[i] = g_0_xx_0_xxxxxx_0[i] * pb_z + g_0_xx_0_xxxxxx_1[i] * wp_z[i];

        g_0_xxz_0_xxxxxy_0[i] = g_0_xx_0_xxxxxy_0[i] * pb_z + g_0_xx_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxz_0_xxxxxz_0[i] = g_0_xx_0_xxxxx_1[i] * fi_abcd_0 + g_0_xx_0_xxxxxz_0[i] * pb_z + g_0_xx_0_xxxxxz_1[i] * wp_z[i];

        g_0_xxz_0_xxxxyy_0[i] = g_0_xx_0_xxxxyy_0[i] * pb_z + g_0_xx_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxz_0_xxxxyz_0[i] = g_0_xx_0_xxxxy_1[i] * fi_abcd_0 + g_0_xx_0_xxxxyz_0[i] * pb_z + g_0_xx_0_xxxxyz_1[i] * wp_z[i];

        g_0_xxz_0_xxxxzz_0[i] = 2.0 * g_0_xx_0_xxxxz_1[i] * fi_abcd_0 + g_0_xx_0_xxxxzz_0[i] * pb_z + g_0_xx_0_xxxxzz_1[i] * wp_z[i];

        g_0_xxz_0_xxxyyy_0[i] = g_0_xx_0_xxxyyy_0[i] * pb_z + g_0_xx_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxz_0_xxxyyz_0[i] = g_0_xx_0_xxxyy_1[i] * fi_abcd_0 + g_0_xx_0_xxxyyz_0[i] * pb_z + g_0_xx_0_xxxyyz_1[i] * wp_z[i];

        g_0_xxz_0_xxxyzz_0[i] = 2.0 * g_0_xx_0_xxxyz_1[i] * fi_abcd_0 + g_0_xx_0_xxxyzz_0[i] * pb_z + g_0_xx_0_xxxyzz_1[i] * wp_z[i];

        g_0_xxz_0_xxxzzz_0[i] = 3.0 * g_0_xx_0_xxxzz_1[i] * fi_abcd_0 + g_0_xx_0_xxxzzz_0[i] * pb_z + g_0_xx_0_xxxzzz_1[i] * wp_z[i];

        g_0_xxz_0_xxyyyy_0[i] = g_0_xx_0_xxyyyy_0[i] * pb_z + g_0_xx_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxz_0_xxyyyz_0[i] = g_0_xx_0_xxyyy_1[i] * fi_abcd_0 + g_0_xx_0_xxyyyz_0[i] * pb_z + g_0_xx_0_xxyyyz_1[i] * wp_z[i];

        g_0_xxz_0_xxyyzz_0[i] = 2.0 * g_0_xx_0_xxyyz_1[i] * fi_abcd_0 + g_0_xx_0_xxyyzz_0[i] * pb_z + g_0_xx_0_xxyyzz_1[i] * wp_z[i];

        g_0_xxz_0_xxyzzz_0[i] = 3.0 * g_0_xx_0_xxyzz_1[i] * fi_abcd_0 + g_0_xx_0_xxyzzz_0[i] * pb_z + g_0_xx_0_xxyzzz_1[i] * wp_z[i];

        g_0_xxz_0_xxzzzz_0[i] = 4.0 * g_0_xx_0_xxzzz_1[i] * fi_abcd_0 + g_0_xx_0_xxzzzz_0[i] * pb_z + g_0_xx_0_xxzzzz_1[i] * wp_z[i];

        g_0_xxz_0_xyyyyy_0[i] = g_0_xx_0_xyyyyy_0[i] * pb_z + g_0_xx_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxz_0_xyyyyz_0[i] = g_0_xx_0_xyyyy_1[i] * fi_abcd_0 + g_0_xx_0_xyyyyz_0[i] * pb_z + g_0_xx_0_xyyyyz_1[i] * wp_z[i];

        g_0_xxz_0_xyyyzz_0[i] = 2.0 * g_0_xx_0_xyyyz_1[i] * fi_abcd_0 + g_0_xx_0_xyyyzz_0[i] * pb_z + g_0_xx_0_xyyyzz_1[i] * wp_z[i];

        g_0_xxz_0_xyyzzz_0[i] = 3.0 * g_0_xx_0_xyyzz_1[i] * fi_abcd_0 + g_0_xx_0_xyyzzz_0[i] * pb_z + g_0_xx_0_xyyzzz_1[i] * wp_z[i];

        g_0_xxz_0_xyzzzz_0[i] = 4.0 * g_0_xx_0_xyzzz_1[i] * fi_abcd_0 + g_0_xx_0_xyzzzz_0[i] * pb_z + g_0_xx_0_xyzzzz_1[i] * wp_z[i];

        g_0_xxz_0_xzzzzz_0[i] = 5.0 * g_0_xx_0_xzzzz_1[i] * fi_abcd_0 + g_0_xx_0_xzzzzz_0[i] * pb_z + g_0_xx_0_xzzzzz_1[i] * wp_z[i];

        g_0_xxz_0_yyyyyy_0[i] = g_0_xx_0_yyyyyy_0[i] * pb_z + g_0_xx_0_yyyyyy_1[i] * wp_z[i];

        g_0_xxz_0_yyyyyz_0[i] = g_0_xx_0_yyyyy_1[i] * fi_abcd_0 + g_0_xx_0_yyyyyz_0[i] * pb_z + g_0_xx_0_yyyyyz_1[i] * wp_z[i];

        g_0_xxz_0_yyyyzz_0[i] = 2.0 * g_0_xx_0_yyyyz_1[i] * fi_abcd_0 + g_0_xx_0_yyyyzz_0[i] * pb_z + g_0_xx_0_yyyyzz_1[i] * wp_z[i];

        g_0_xxz_0_yyyzzz_0[i] = 3.0 * g_0_xx_0_yyyzz_1[i] * fi_abcd_0 + g_0_xx_0_yyyzzz_0[i] * pb_z + g_0_xx_0_yyyzzz_1[i] * wp_z[i];

        g_0_xxz_0_yyzzzz_0[i] = 4.0 * g_0_xx_0_yyzzz_1[i] * fi_abcd_0 + g_0_xx_0_yyzzzz_0[i] * pb_z + g_0_xx_0_yyzzzz_1[i] * wp_z[i];

        g_0_xxz_0_yzzzzz_0[i] = 5.0 * g_0_xx_0_yzzzz_1[i] * fi_abcd_0 + g_0_xx_0_yzzzzz_0[i] * pb_z + g_0_xx_0_yzzzzz_1[i] * wp_z[i];

        g_0_xxz_0_zzzzzz_0[i] = 6.0 * g_0_xx_0_zzzzz_1[i] * fi_abcd_0 + g_0_xx_0_zzzzzz_0[i] * pb_z + g_0_xx_0_zzzzzz_1[i] * wp_z[i];
    }

    /// Set up 84-112 components of targeted buffer : SFSI

    auto g_0_xyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sfsi + 84);

    auto g_0_xyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sfsi + 85);

    auto g_0_xyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sfsi + 86);

    auto g_0_xyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sfsi + 87);

    auto g_0_xyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sfsi + 88);

    auto g_0_xyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sfsi + 89);

    auto g_0_xyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sfsi + 90);

    auto g_0_xyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sfsi + 91);

    auto g_0_xyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sfsi + 92);

    auto g_0_xyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sfsi + 93);

    auto g_0_xyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 94);

    auto g_0_xyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 95);

    auto g_0_xyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 96);

    auto g_0_xyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 97);

    auto g_0_xyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 98);

    auto g_0_xyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 99);

    auto g_0_xyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 100);

    auto g_0_xyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 101);

    auto g_0_xyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 102);

    auto g_0_xyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 103);

    auto g_0_xyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 104);

    auto g_0_xyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 105);

    auto g_0_xyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 106);

    auto g_0_xyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 107);

    auto g_0_xyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 108);

    auto g_0_xyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 109);

    auto g_0_xyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 110);

    auto g_0_xyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 111);

    #pragma omp simd aligned(g_0_xyy_0_xxxxxx_0, g_0_xyy_0_xxxxxy_0, g_0_xyy_0_xxxxxz_0, g_0_xyy_0_xxxxyy_0, g_0_xyy_0_xxxxyz_0, g_0_xyy_0_xxxxzz_0, g_0_xyy_0_xxxyyy_0, g_0_xyy_0_xxxyyz_0, g_0_xyy_0_xxxyzz_0, g_0_xyy_0_xxxzzz_0, g_0_xyy_0_xxyyyy_0, g_0_xyy_0_xxyyyz_0, g_0_xyy_0_xxyyzz_0, g_0_xyy_0_xxyzzz_0, g_0_xyy_0_xxzzzz_0, g_0_xyy_0_xyyyyy_0, g_0_xyy_0_xyyyyz_0, g_0_xyy_0_xyyyzz_0, g_0_xyy_0_xyyzzz_0, g_0_xyy_0_xyzzzz_0, g_0_xyy_0_xzzzzz_0, g_0_xyy_0_yyyyyy_0, g_0_xyy_0_yyyyyz_0, g_0_xyy_0_yyyyzz_0, g_0_xyy_0_yyyzzz_0, g_0_xyy_0_yyzzzz_0, g_0_xyy_0_yzzzzz_0, g_0_xyy_0_zzzzzz_0, g_0_yy_0_xxxxx_1, g_0_yy_0_xxxxxx_0, g_0_yy_0_xxxxxx_1, g_0_yy_0_xxxxxy_0, g_0_yy_0_xxxxxy_1, g_0_yy_0_xxxxxz_0, g_0_yy_0_xxxxxz_1, g_0_yy_0_xxxxy_1, g_0_yy_0_xxxxyy_0, g_0_yy_0_xxxxyy_1, g_0_yy_0_xxxxyz_0, g_0_yy_0_xxxxyz_1, g_0_yy_0_xxxxz_1, g_0_yy_0_xxxxzz_0, g_0_yy_0_xxxxzz_1, g_0_yy_0_xxxyy_1, g_0_yy_0_xxxyyy_0, g_0_yy_0_xxxyyy_1, g_0_yy_0_xxxyyz_0, g_0_yy_0_xxxyyz_1, g_0_yy_0_xxxyz_1, g_0_yy_0_xxxyzz_0, g_0_yy_0_xxxyzz_1, g_0_yy_0_xxxzz_1, g_0_yy_0_xxxzzz_0, g_0_yy_0_xxxzzz_1, g_0_yy_0_xxyyy_1, g_0_yy_0_xxyyyy_0, g_0_yy_0_xxyyyy_1, g_0_yy_0_xxyyyz_0, g_0_yy_0_xxyyyz_1, g_0_yy_0_xxyyz_1, g_0_yy_0_xxyyzz_0, g_0_yy_0_xxyyzz_1, g_0_yy_0_xxyzz_1, g_0_yy_0_xxyzzz_0, g_0_yy_0_xxyzzz_1, g_0_yy_0_xxzzz_1, g_0_yy_0_xxzzzz_0, g_0_yy_0_xxzzzz_1, g_0_yy_0_xyyyy_1, g_0_yy_0_xyyyyy_0, g_0_yy_0_xyyyyy_1, g_0_yy_0_xyyyyz_0, g_0_yy_0_xyyyyz_1, g_0_yy_0_xyyyz_1, g_0_yy_0_xyyyzz_0, g_0_yy_0_xyyyzz_1, g_0_yy_0_xyyzz_1, g_0_yy_0_xyyzzz_0, g_0_yy_0_xyyzzz_1, g_0_yy_0_xyzzz_1, g_0_yy_0_xyzzzz_0, g_0_yy_0_xyzzzz_1, g_0_yy_0_xzzzz_1, g_0_yy_0_xzzzzz_0, g_0_yy_0_xzzzzz_1, g_0_yy_0_yyyyy_1, g_0_yy_0_yyyyyy_0, g_0_yy_0_yyyyyy_1, g_0_yy_0_yyyyyz_0, g_0_yy_0_yyyyyz_1, g_0_yy_0_yyyyz_1, g_0_yy_0_yyyyzz_0, g_0_yy_0_yyyyzz_1, g_0_yy_0_yyyzz_1, g_0_yy_0_yyyzzz_0, g_0_yy_0_yyyzzz_1, g_0_yy_0_yyzzz_1, g_0_yy_0_yyzzzz_0, g_0_yy_0_yyzzzz_1, g_0_yy_0_yzzzz_1, g_0_yy_0_yzzzzz_0, g_0_yy_0_yzzzzz_1, g_0_yy_0_zzzzz_1, g_0_yy_0_zzzzzz_0, g_0_yy_0_zzzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyy_0_xxxxxx_0[i] = 6.0 * g_0_yy_0_xxxxx_1[i] * fi_abcd_0 + g_0_yy_0_xxxxxx_0[i] * pb_x + g_0_yy_0_xxxxxx_1[i] * wp_x[i];

        g_0_xyy_0_xxxxxy_0[i] = 5.0 * g_0_yy_0_xxxxy_1[i] * fi_abcd_0 + g_0_yy_0_xxxxxy_0[i] * pb_x + g_0_yy_0_xxxxxy_1[i] * wp_x[i];

        g_0_xyy_0_xxxxxz_0[i] = 5.0 * g_0_yy_0_xxxxz_1[i] * fi_abcd_0 + g_0_yy_0_xxxxxz_0[i] * pb_x + g_0_yy_0_xxxxxz_1[i] * wp_x[i];

        g_0_xyy_0_xxxxyy_0[i] = 4.0 * g_0_yy_0_xxxyy_1[i] * fi_abcd_0 + g_0_yy_0_xxxxyy_0[i] * pb_x + g_0_yy_0_xxxxyy_1[i] * wp_x[i];

        g_0_xyy_0_xxxxyz_0[i] = 4.0 * g_0_yy_0_xxxyz_1[i] * fi_abcd_0 + g_0_yy_0_xxxxyz_0[i] * pb_x + g_0_yy_0_xxxxyz_1[i] * wp_x[i];

        g_0_xyy_0_xxxxzz_0[i] = 4.0 * g_0_yy_0_xxxzz_1[i] * fi_abcd_0 + g_0_yy_0_xxxxzz_0[i] * pb_x + g_0_yy_0_xxxxzz_1[i] * wp_x[i];

        g_0_xyy_0_xxxyyy_0[i] = 3.0 * g_0_yy_0_xxyyy_1[i] * fi_abcd_0 + g_0_yy_0_xxxyyy_0[i] * pb_x + g_0_yy_0_xxxyyy_1[i] * wp_x[i];

        g_0_xyy_0_xxxyyz_0[i] = 3.0 * g_0_yy_0_xxyyz_1[i] * fi_abcd_0 + g_0_yy_0_xxxyyz_0[i] * pb_x + g_0_yy_0_xxxyyz_1[i] * wp_x[i];

        g_0_xyy_0_xxxyzz_0[i] = 3.0 * g_0_yy_0_xxyzz_1[i] * fi_abcd_0 + g_0_yy_0_xxxyzz_0[i] * pb_x + g_0_yy_0_xxxyzz_1[i] * wp_x[i];

        g_0_xyy_0_xxxzzz_0[i] = 3.0 * g_0_yy_0_xxzzz_1[i] * fi_abcd_0 + g_0_yy_0_xxxzzz_0[i] * pb_x + g_0_yy_0_xxxzzz_1[i] * wp_x[i];

        g_0_xyy_0_xxyyyy_0[i] = 2.0 * g_0_yy_0_xyyyy_1[i] * fi_abcd_0 + g_0_yy_0_xxyyyy_0[i] * pb_x + g_0_yy_0_xxyyyy_1[i] * wp_x[i];

        g_0_xyy_0_xxyyyz_0[i] = 2.0 * g_0_yy_0_xyyyz_1[i] * fi_abcd_0 + g_0_yy_0_xxyyyz_0[i] * pb_x + g_0_yy_0_xxyyyz_1[i] * wp_x[i];

        g_0_xyy_0_xxyyzz_0[i] = 2.0 * g_0_yy_0_xyyzz_1[i] * fi_abcd_0 + g_0_yy_0_xxyyzz_0[i] * pb_x + g_0_yy_0_xxyyzz_1[i] * wp_x[i];

        g_0_xyy_0_xxyzzz_0[i] = 2.0 * g_0_yy_0_xyzzz_1[i] * fi_abcd_0 + g_0_yy_0_xxyzzz_0[i] * pb_x + g_0_yy_0_xxyzzz_1[i] * wp_x[i];

        g_0_xyy_0_xxzzzz_0[i] = 2.0 * g_0_yy_0_xzzzz_1[i] * fi_abcd_0 + g_0_yy_0_xxzzzz_0[i] * pb_x + g_0_yy_0_xxzzzz_1[i] * wp_x[i];

        g_0_xyy_0_xyyyyy_0[i] = g_0_yy_0_yyyyy_1[i] * fi_abcd_0 + g_0_yy_0_xyyyyy_0[i] * pb_x + g_0_yy_0_xyyyyy_1[i] * wp_x[i];

        g_0_xyy_0_xyyyyz_0[i] = g_0_yy_0_yyyyz_1[i] * fi_abcd_0 + g_0_yy_0_xyyyyz_0[i] * pb_x + g_0_yy_0_xyyyyz_1[i] * wp_x[i];

        g_0_xyy_0_xyyyzz_0[i] = g_0_yy_0_yyyzz_1[i] * fi_abcd_0 + g_0_yy_0_xyyyzz_0[i] * pb_x + g_0_yy_0_xyyyzz_1[i] * wp_x[i];

        g_0_xyy_0_xyyzzz_0[i] = g_0_yy_0_yyzzz_1[i] * fi_abcd_0 + g_0_yy_0_xyyzzz_0[i] * pb_x + g_0_yy_0_xyyzzz_1[i] * wp_x[i];

        g_0_xyy_0_xyzzzz_0[i] = g_0_yy_0_yzzzz_1[i] * fi_abcd_0 + g_0_yy_0_xyzzzz_0[i] * pb_x + g_0_yy_0_xyzzzz_1[i] * wp_x[i];

        g_0_xyy_0_xzzzzz_0[i] = g_0_yy_0_zzzzz_1[i] * fi_abcd_0 + g_0_yy_0_xzzzzz_0[i] * pb_x + g_0_yy_0_xzzzzz_1[i] * wp_x[i];

        g_0_xyy_0_yyyyyy_0[i] = g_0_yy_0_yyyyyy_0[i] * pb_x + g_0_yy_0_yyyyyy_1[i] * wp_x[i];

        g_0_xyy_0_yyyyyz_0[i] = g_0_yy_0_yyyyyz_0[i] * pb_x + g_0_yy_0_yyyyyz_1[i] * wp_x[i];

        g_0_xyy_0_yyyyzz_0[i] = g_0_yy_0_yyyyzz_0[i] * pb_x + g_0_yy_0_yyyyzz_1[i] * wp_x[i];

        g_0_xyy_0_yyyzzz_0[i] = g_0_yy_0_yyyzzz_0[i] * pb_x + g_0_yy_0_yyyzzz_1[i] * wp_x[i];

        g_0_xyy_0_yyzzzz_0[i] = g_0_yy_0_yyzzzz_0[i] * pb_x + g_0_yy_0_yyzzzz_1[i] * wp_x[i];

        g_0_xyy_0_yzzzzz_0[i] = g_0_yy_0_yzzzzz_0[i] * pb_x + g_0_yy_0_yzzzzz_1[i] * wp_x[i];

        g_0_xyy_0_zzzzzz_0[i] = g_0_yy_0_zzzzzz_0[i] * pb_x + g_0_yy_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 112-140 components of targeted buffer : SFSI

    auto g_0_xyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sfsi + 112);

    auto g_0_xyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sfsi + 113);

    auto g_0_xyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sfsi + 114);

    auto g_0_xyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sfsi + 115);

    auto g_0_xyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sfsi + 116);

    auto g_0_xyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sfsi + 117);

    auto g_0_xyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sfsi + 118);

    auto g_0_xyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sfsi + 119);

    auto g_0_xyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sfsi + 120);

    auto g_0_xyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sfsi + 121);

    auto g_0_xyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 122);

    auto g_0_xyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 123);

    auto g_0_xyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 124);

    auto g_0_xyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 125);

    auto g_0_xyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 126);

    auto g_0_xyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 127);

    auto g_0_xyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 128);

    auto g_0_xyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 129);

    auto g_0_xyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 130);

    auto g_0_xyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 131);

    auto g_0_xyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 132);

    auto g_0_xyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 133);

    auto g_0_xyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 134);

    auto g_0_xyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 135);

    auto g_0_xyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 136);

    auto g_0_xyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 137);

    auto g_0_xyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 138);

    auto g_0_xyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 139);

    #pragma omp simd aligned(g_0_xy_0_xxxxxy_0, g_0_xy_0_xxxxxy_1, g_0_xy_0_xxxxyy_0, g_0_xy_0_xxxxyy_1, g_0_xy_0_xxxyyy_0, g_0_xy_0_xxxyyy_1, g_0_xy_0_xxyyyy_0, g_0_xy_0_xxyyyy_1, g_0_xy_0_xyyyyy_0, g_0_xy_0_xyyyyy_1, g_0_xyz_0_xxxxxx_0, g_0_xyz_0_xxxxxy_0, g_0_xyz_0_xxxxxz_0, g_0_xyz_0_xxxxyy_0, g_0_xyz_0_xxxxyz_0, g_0_xyz_0_xxxxzz_0, g_0_xyz_0_xxxyyy_0, g_0_xyz_0_xxxyyz_0, g_0_xyz_0_xxxyzz_0, g_0_xyz_0_xxxzzz_0, g_0_xyz_0_xxyyyy_0, g_0_xyz_0_xxyyyz_0, g_0_xyz_0_xxyyzz_0, g_0_xyz_0_xxyzzz_0, g_0_xyz_0_xxzzzz_0, g_0_xyz_0_xyyyyy_0, g_0_xyz_0_xyyyyz_0, g_0_xyz_0_xyyyzz_0, g_0_xyz_0_xyyzzz_0, g_0_xyz_0_xyzzzz_0, g_0_xyz_0_xzzzzz_0, g_0_xyz_0_yyyyyy_0, g_0_xyz_0_yyyyyz_0, g_0_xyz_0_yyyyzz_0, g_0_xyz_0_yyyzzz_0, g_0_xyz_0_yyzzzz_0, g_0_xyz_0_yzzzzz_0, g_0_xyz_0_zzzzzz_0, g_0_xz_0_xxxxxx_0, g_0_xz_0_xxxxxx_1, g_0_xz_0_xxxxxz_0, g_0_xz_0_xxxxxz_1, g_0_xz_0_xxxxzz_0, g_0_xz_0_xxxxzz_1, g_0_xz_0_xxxzzz_0, g_0_xz_0_xxxzzz_1, g_0_xz_0_xxzzzz_0, g_0_xz_0_xxzzzz_1, g_0_xz_0_xzzzzz_0, g_0_xz_0_xzzzzz_1, g_0_yz_0_xxxxyz_0, g_0_yz_0_xxxxyz_1, g_0_yz_0_xxxyyz_0, g_0_yz_0_xxxyyz_1, g_0_yz_0_xxxyz_1, g_0_yz_0_xxxyzz_0, g_0_yz_0_xxxyzz_1, g_0_yz_0_xxyyyz_0, g_0_yz_0_xxyyyz_1, g_0_yz_0_xxyyz_1, g_0_yz_0_xxyyzz_0, g_0_yz_0_xxyyzz_1, g_0_yz_0_xxyzz_1, g_0_yz_0_xxyzzz_0, g_0_yz_0_xxyzzz_1, g_0_yz_0_xyyyyz_0, g_0_yz_0_xyyyyz_1, g_0_yz_0_xyyyz_1, g_0_yz_0_xyyyzz_0, g_0_yz_0_xyyyzz_1, g_0_yz_0_xyyzz_1, g_0_yz_0_xyyzzz_0, g_0_yz_0_xyyzzz_1, g_0_yz_0_xyzzz_1, g_0_yz_0_xyzzzz_0, g_0_yz_0_xyzzzz_1, g_0_yz_0_yyyyyy_0, g_0_yz_0_yyyyyy_1, g_0_yz_0_yyyyyz_0, g_0_yz_0_yyyyyz_1, g_0_yz_0_yyyyz_1, g_0_yz_0_yyyyzz_0, g_0_yz_0_yyyyzz_1, g_0_yz_0_yyyzz_1, g_0_yz_0_yyyzzz_0, g_0_yz_0_yyyzzz_1, g_0_yz_0_yyzzz_1, g_0_yz_0_yyzzzz_0, g_0_yz_0_yyzzzz_1, g_0_yz_0_yzzzz_1, g_0_yz_0_yzzzzz_0, g_0_yz_0_yzzzzz_1, g_0_yz_0_zzzzzz_0, g_0_yz_0_zzzzzz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyz_0_xxxxxx_0[i] = g_0_xz_0_xxxxxx_0[i] * pb_y + g_0_xz_0_xxxxxx_1[i] * wp_y[i];

        g_0_xyz_0_xxxxxy_0[i] = g_0_xy_0_xxxxxy_0[i] * pb_z + g_0_xy_0_xxxxxy_1[i] * wp_z[i];

        g_0_xyz_0_xxxxxz_0[i] = g_0_xz_0_xxxxxz_0[i] * pb_y + g_0_xz_0_xxxxxz_1[i] * wp_y[i];

        g_0_xyz_0_xxxxyy_0[i] = g_0_xy_0_xxxxyy_0[i] * pb_z + g_0_xy_0_xxxxyy_1[i] * wp_z[i];

        g_0_xyz_0_xxxxyz_0[i] = 4.0 * g_0_yz_0_xxxyz_1[i] * fi_abcd_0 + g_0_yz_0_xxxxyz_0[i] * pb_x + g_0_yz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xyz_0_xxxxzz_0[i] = g_0_xz_0_xxxxzz_0[i] * pb_y + g_0_xz_0_xxxxzz_1[i] * wp_y[i];

        g_0_xyz_0_xxxyyy_0[i] = g_0_xy_0_xxxyyy_0[i] * pb_z + g_0_xy_0_xxxyyy_1[i] * wp_z[i];

        g_0_xyz_0_xxxyyz_0[i] = 3.0 * g_0_yz_0_xxyyz_1[i] * fi_abcd_0 + g_0_yz_0_xxxyyz_0[i] * pb_x + g_0_yz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xyz_0_xxxyzz_0[i] = 3.0 * g_0_yz_0_xxyzz_1[i] * fi_abcd_0 + g_0_yz_0_xxxyzz_0[i] * pb_x + g_0_yz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xyz_0_xxxzzz_0[i] = g_0_xz_0_xxxzzz_0[i] * pb_y + g_0_xz_0_xxxzzz_1[i] * wp_y[i];

        g_0_xyz_0_xxyyyy_0[i] = g_0_xy_0_xxyyyy_0[i] * pb_z + g_0_xy_0_xxyyyy_1[i] * wp_z[i];

        g_0_xyz_0_xxyyyz_0[i] = 2.0 * g_0_yz_0_xyyyz_1[i] * fi_abcd_0 + g_0_yz_0_xxyyyz_0[i] * pb_x + g_0_yz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xyz_0_xxyyzz_0[i] = 2.0 * g_0_yz_0_xyyzz_1[i] * fi_abcd_0 + g_0_yz_0_xxyyzz_0[i] * pb_x + g_0_yz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xyz_0_xxyzzz_0[i] = 2.0 * g_0_yz_0_xyzzz_1[i] * fi_abcd_0 + g_0_yz_0_xxyzzz_0[i] * pb_x + g_0_yz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xyz_0_xxzzzz_0[i] = g_0_xz_0_xxzzzz_0[i] * pb_y + g_0_xz_0_xxzzzz_1[i] * wp_y[i];

        g_0_xyz_0_xyyyyy_0[i] = g_0_xy_0_xyyyyy_0[i] * pb_z + g_0_xy_0_xyyyyy_1[i] * wp_z[i];

        g_0_xyz_0_xyyyyz_0[i] = g_0_yz_0_yyyyz_1[i] * fi_abcd_0 + g_0_yz_0_xyyyyz_0[i] * pb_x + g_0_yz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xyz_0_xyyyzz_0[i] = g_0_yz_0_yyyzz_1[i] * fi_abcd_0 + g_0_yz_0_xyyyzz_0[i] * pb_x + g_0_yz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xyz_0_xyyzzz_0[i] = g_0_yz_0_yyzzz_1[i] * fi_abcd_0 + g_0_yz_0_xyyzzz_0[i] * pb_x + g_0_yz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xyz_0_xyzzzz_0[i] = g_0_yz_0_yzzzz_1[i] * fi_abcd_0 + g_0_yz_0_xyzzzz_0[i] * pb_x + g_0_yz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xyz_0_xzzzzz_0[i] = g_0_xz_0_xzzzzz_0[i] * pb_y + g_0_xz_0_xzzzzz_1[i] * wp_y[i];

        g_0_xyz_0_yyyyyy_0[i] = g_0_yz_0_yyyyyy_0[i] * pb_x + g_0_yz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xyz_0_yyyyyz_0[i] = g_0_yz_0_yyyyyz_0[i] * pb_x + g_0_yz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xyz_0_yyyyzz_0[i] = g_0_yz_0_yyyyzz_0[i] * pb_x + g_0_yz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xyz_0_yyyzzz_0[i] = g_0_yz_0_yyyzzz_0[i] * pb_x + g_0_yz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xyz_0_yyzzzz_0[i] = g_0_yz_0_yyzzzz_0[i] * pb_x + g_0_yz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xyz_0_yzzzzz_0[i] = g_0_yz_0_yzzzzz_0[i] * pb_x + g_0_yz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xyz_0_zzzzzz_0[i] = g_0_yz_0_zzzzzz_0[i] * pb_x + g_0_yz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 140-168 components of targeted buffer : SFSI

    auto g_0_xzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sfsi + 140);

    auto g_0_xzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sfsi + 141);

    auto g_0_xzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sfsi + 142);

    auto g_0_xzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sfsi + 143);

    auto g_0_xzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sfsi + 144);

    auto g_0_xzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sfsi + 145);

    auto g_0_xzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sfsi + 146);

    auto g_0_xzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sfsi + 147);

    auto g_0_xzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sfsi + 148);

    auto g_0_xzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sfsi + 149);

    auto g_0_xzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 150);

    auto g_0_xzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 151);

    auto g_0_xzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 152);

    auto g_0_xzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 153);

    auto g_0_xzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 154);

    auto g_0_xzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 155);

    auto g_0_xzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 156);

    auto g_0_xzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 157);

    auto g_0_xzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 158);

    auto g_0_xzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 159);

    auto g_0_xzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 160);

    auto g_0_xzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 161);

    auto g_0_xzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 162);

    auto g_0_xzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 163);

    auto g_0_xzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 164);

    auto g_0_xzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 165);

    auto g_0_xzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 166);

    auto g_0_xzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 167);

    #pragma omp simd aligned(g_0_xzz_0_xxxxxx_0, g_0_xzz_0_xxxxxy_0, g_0_xzz_0_xxxxxz_0, g_0_xzz_0_xxxxyy_0, g_0_xzz_0_xxxxyz_0, g_0_xzz_0_xxxxzz_0, g_0_xzz_0_xxxyyy_0, g_0_xzz_0_xxxyyz_0, g_0_xzz_0_xxxyzz_0, g_0_xzz_0_xxxzzz_0, g_0_xzz_0_xxyyyy_0, g_0_xzz_0_xxyyyz_0, g_0_xzz_0_xxyyzz_0, g_0_xzz_0_xxyzzz_0, g_0_xzz_0_xxzzzz_0, g_0_xzz_0_xyyyyy_0, g_0_xzz_0_xyyyyz_0, g_0_xzz_0_xyyyzz_0, g_0_xzz_0_xyyzzz_0, g_0_xzz_0_xyzzzz_0, g_0_xzz_0_xzzzzz_0, g_0_xzz_0_yyyyyy_0, g_0_xzz_0_yyyyyz_0, g_0_xzz_0_yyyyzz_0, g_0_xzz_0_yyyzzz_0, g_0_xzz_0_yyzzzz_0, g_0_xzz_0_yzzzzz_0, g_0_xzz_0_zzzzzz_0, g_0_zz_0_xxxxx_1, g_0_zz_0_xxxxxx_0, g_0_zz_0_xxxxxx_1, g_0_zz_0_xxxxxy_0, g_0_zz_0_xxxxxy_1, g_0_zz_0_xxxxxz_0, g_0_zz_0_xxxxxz_1, g_0_zz_0_xxxxy_1, g_0_zz_0_xxxxyy_0, g_0_zz_0_xxxxyy_1, g_0_zz_0_xxxxyz_0, g_0_zz_0_xxxxyz_1, g_0_zz_0_xxxxz_1, g_0_zz_0_xxxxzz_0, g_0_zz_0_xxxxzz_1, g_0_zz_0_xxxyy_1, g_0_zz_0_xxxyyy_0, g_0_zz_0_xxxyyy_1, g_0_zz_0_xxxyyz_0, g_0_zz_0_xxxyyz_1, g_0_zz_0_xxxyz_1, g_0_zz_0_xxxyzz_0, g_0_zz_0_xxxyzz_1, g_0_zz_0_xxxzz_1, g_0_zz_0_xxxzzz_0, g_0_zz_0_xxxzzz_1, g_0_zz_0_xxyyy_1, g_0_zz_0_xxyyyy_0, g_0_zz_0_xxyyyy_1, g_0_zz_0_xxyyyz_0, g_0_zz_0_xxyyyz_1, g_0_zz_0_xxyyz_1, g_0_zz_0_xxyyzz_0, g_0_zz_0_xxyyzz_1, g_0_zz_0_xxyzz_1, g_0_zz_0_xxyzzz_0, g_0_zz_0_xxyzzz_1, g_0_zz_0_xxzzz_1, g_0_zz_0_xxzzzz_0, g_0_zz_0_xxzzzz_1, g_0_zz_0_xyyyy_1, g_0_zz_0_xyyyyy_0, g_0_zz_0_xyyyyy_1, g_0_zz_0_xyyyyz_0, g_0_zz_0_xyyyyz_1, g_0_zz_0_xyyyz_1, g_0_zz_0_xyyyzz_0, g_0_zz_0_xyyyzz_1, g_0_zz_0_xyyzz_1, g_0_zz_0_xyyzzz_0, g_0_zz_0_xyyzzz_1, g_0_zz_0_xyzzz_1, g_0_zz_0_xyzzzz_0, g_0_zz_0_xyzzzz_1, g_0_zz_0_xzzzz_1, g_0_zz_0_xzzzzz_0, g_0_zz_0_xzzzzz_1, g_0_zz_0_yyyyy_1, g_0_zz_0_yyyyyy_0, g_0_zz_0_yyyyyy_1, g_0_zz_0_yyyyyz_0, g_0_zz_0_yyyyyz_1, g_0_zz_0_yyyyz_1, g_0_zz_0_yyyyzz_0, g_0_zz_0_yyyyzz_1, g_0_zz_0_yyyzz_1, g_0_zz_0_yyyzzz_0, g_0_zz_0_yyyzzz_1, g_0_zz_0_yyzzz_1, g_0_zz_0_yyzzzz_0, g_0_zz_0_yyzzzz_1, g_0_zz_0_yzzzz_1, g_0_zz_0_yzzzzz_0, g_0_zz_0_yzzzzz_1, g_0_zz_0_zzzzz_1, g_0_zz_0_zzzzzz_0, g_0_zz_0_zzzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzz_0_xxxxxx_0[i] = 6.0 * g_0_zz_0_xxxxx_1[i] * fi_abcd_0 + g_0_zz_0_xxxxxx_0[i] * pb_x + g_0_zz_0_xxxxxx_1[i] * wp_x[i];

        g_0_xzz_0_xxxxxy_0[i] = 5.0 * g_0_zz_0_xxxxy_1[i] * fi_abcd_0 + g_0_zz_0_xxxxxy_0[i] * pb_x + g_0_zz_0_xxxxxy_1[i] * wp_x[i];

        g_0_xzz_0_xxxxxz_0[i] = 5.0 * g_0_zz_0_xxxxz_1[i] * fi_abcd_0 + g_0_zz_0_xxxxxz_0[i] * pb_x + g_0_zz_0_xxxxxz_1[i] * wp_x[i];

        g_0_xzz_0_xxxxyy_0[i] = 4.0 * g_0_zz_0_xxxyy_1[i] * fi_abcd_0 + g_0_zz_0_xxxxyy_0[i] * pb_x + g_0_zz_0_xxxxyy_1[i] * wp_x[i];

        g_0_xzz_0_xxxxyz_0[i] = 4.0 * g_0_zz_0_xxxyz_1[i] * fi_abcd_0 + g_0_zz_0_xxxxyz_0[i] * pb_x + g_0_zz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xzz_0_xxxxzz_0[i] = 4.0 * g_0_zz_0_xxxzz_1[i] * fi_abcd_0 + g_0_zz_0_xxxxzz_0[i] * pb_x + g_0_zz_0_xxxxzz_1[i] * wp_x[i];

        g_0_xzz_0_xxxyyy_0[i] = 3.0 * g_0_zz_0_xxyyy_1[i] * fi_abcd_0 + g_0_zz_0_xxxyyy_0[i] * pb_x + g_0_zz_0_xxxyyy_1[i] * wp_x[i];

        g_0_xzz_0_xxxyyz_0[i] = 3.0 * g_0_zz_0_xxyyz_1[i] * fi_abcd_0 + g_0_zz_0_xxxyyz_0[i] * pb_x + g_0_zz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xzz_0_xxxyzz_0[i] = 3.0 * g_0_zz_0_xxyzz_1[i] * fi_abcd_0 + g_0_zz_0_xxxyzz_0[i] * pb_x + g_0_zz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xzz_0_xxxzzz_0[i] = 3.0 * g_0_zz_0_xxzzz_1[i] * fi_abcd_0 + g_0_zz_0_xxxzzz_0[i] * pb_x + g_0_zz_0_xxxzzz_1[i] * wp_x[i];

        g_0_xzz_0_xxyyyy_0[i] = 2.0 * g_0_zz_0_xyyyy_1[i] * fi_abcd_0 + g_0_zz_0_xxyyyy_0[i] * pb_x + g_0_zz_0_xxyyyy_1[i] * wp_x[i];

        g_0_xzz_0_xxyyyz_0[i] = 2.0 * g_0_zz_0_xyyyz_1[i] * fi_abcd_0 + g_0_zz_0_xxyyyz_0[i] * pb_x + g_0_zz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xzz_0_xxyyzz_0[i] = 2.0 * g_0_zz_0_xyyzz_1[i] * fi_abcd_0 + g_0_zz_0_xxyyzz_0[i] * pb_x + g_0_zz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xzz_0_xxyzzz_0[i] = 2.0 * g_0_zz_0_xyzzz_1[i] * fi_abcd_0 + g_0_zz_0_xxyzzz_0[i] * pb_x + g_0_zz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xzz_0_xxzzzz_0[i] = 2.0 * g_0_zz_0_xzzzz_1[i] * fi_abcd_0 + g_0_zz_0_xxzzzz_0[i] * pb_x + g_0_zz_0_xxzzzz_1[i] * wp_x[i];

        g_0_xzz_0_xyyyyy_0[i] = g_0_zz_0_yyyyy_1[i] * fi_abcd_0 + g_0_zz_0_xyyyyy_0[i] * pb_x + g_0_zz_0_xyyyyy_1[i] * wp_x[i];

        g_0_xzz_0_xyyyyz_0[i] = g_0_zz_0_yyyyz_1[i] * fi_abcd_0 + g_0_zz_0_xyyyyz_0[i] * pb_x + g_0_zz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xzz_0_xyyyzz_0[i] = g_0_zz_0_yyyzz_1[i] * fi_abcd_0 + g_0_zz_0_xyyyzz_0[i] * pb_x + g_0_zz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xzz_0_xyyzzz_0[i] = g_0_zz_0_yyzzz_1[i] * fi_abcd_0 + g_0_zz_0_xyyzzz_0[i] * pb_x + g_0_zz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xzz_0_xyzzzz_0[i] = g_0_zz_0_yzzzz_1[i] * fi_abcd_0 + g_0_zz_0_xyzzzz_0[i] * pb_x + g_0_zz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xzz_0_xzzzzz_0[i] = g_0_zz_0_zzzzz_1[i] * fi_abcd_0 + g_0_zz_0_xzzzzz_0[i] * pb_x + g_0_zz_0_xzzzzz_1[i] * wp_x[i];

        g_0_xzz_0_yyyyyy_0[i] = g_0_zz_0_yyyyyy_0[i] * pb_x + g_0_zz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xzz_0_yyyyyz_0[i] = g_0_zz_0_yyyyyz_0[i] * pb_x + g_0_zz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xzz_0_yyyyzz_0[i] = g_0_zz_0_yyyyzz_0[i] * pb_x + g_0_zz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xzz_0_yyyzzz_0[i] = g_0_zz_0_yyyzzz_0[i] * pb_x + g_0_zz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xzz_0_yyzzzz_0[i] = g_0_zz_0_yyzzzz_0[i] * pb_x + g_0_zz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xzz_0_yzzzzz_0[i] = g_0_zz_0_yzzzzz_0[i] * pb_x + g_0_zz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xzz_0_zzzzzz_0[i] = g_0_zz_0_zzzzzz_0[i] * pb_x + g_0_zz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 168-196 components of targeted buffer : SFSI

    auto g_0_yyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sfsi + 168);

    auto g_0_yyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sfsi + 169);

    auto g_0_yyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sfsi + 170);

    auto g_0_yyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sfsi + 171);

    auto g_0_yyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sfsi + 172);

    auto g_0_yyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sfsi + 173);

    auto g_0_yyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sfsi + 174);

    auto g_0_yyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sfsi + 175);

    auto g_0_yyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sfsi + 176);

    auto g_0_yyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sfsi + 177);

    auto g_0_yyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 178);

    auto g_0_yyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 179);

    auto g_0_yyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 180);

    auto g_0_yyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 181);

    auto g_0_yyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 182);

    auto g_0_yyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 183);

    auto g_0_yyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 184);

    auto g_0_yyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 185);

    auto g_0_yyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 186);

    auto g_0_yyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 187);

    auto g_0_yyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 188);

    auto g_0_yyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 189);

    auto g_0_yyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 190);

    auto g_0_yyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 191);

    auto g_0_yyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 192);

    auto g_0_yyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 193);

    auto g_0_yyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 194);

    auto g_0_yyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 195);

    #pragma omp simd aligned(g_0_y_0_xxxxxx_0, g_0_y_0_xxxxxx_1, g_0_y_0_xxxxxy_0, g_0_y_0_xxxxxy_1, g_0_y_0_xxxxxz_0, g_0_y_0_xxxxxz_1, g_0_y_0_xxxxyy_0, g_0_y_0_xxxxyy_1, g_0_y_0_xxxxyz_0, g_0_y_0_xxxxyz_1, g_0_y_0_xxxxzz_0, g_0_y_0_xxxxzz_1, g_0_y_0_xxxyyy_0, g_0_y_0_xxxyyy_1, g_0_y_0_xxxyyz_0, g_0_y_0_xxxyyz_1, g_0_y_0_xxxyzz_0, g_0_y_0_xxxyzz_1, g_0_y_0_xxxzzz_0, g_0_y_0_xxxzzz_1, g_0_y_0_xxyyyy_0, g_0_y_0_xxyyyy_1, g_0_y_0_xxyyyz_0, g_0_y_0_xxyyyz_1, g_0_y_0_xxyyzz_0, g_0_y_0_xxyyzz_1, g_0_y_0_xxyzzz_0, g_0_y_0_xxyzzz_1, g_0_y_0_xxzzzz_0, g_0_y_0_xxzzzz_1, g_0_y_0_xyyyyy_0, g_0_y_0_xyyyyy_1, g_0_y_0_xyyyyz_0, g_0_y_0_xyyyyz_1, g_0_y_0_xyyyzz_0, g_0_y_0_xyyyzz_1, g_0_y_0_xyyzzz_0, g_0_y_0_xyyzzz_1, g_0_y_0_xyzzzz_0, g_0_y_0_xyzzzz_1, g_0_y_0_xzzzzz_0, g_0_y_0_xzzzzz_1, g_0_y_0_yyyyyy_0, g_0_y_0_yyyyyy_1, g_0_y_0_yyyyyz_0, g_0_y_0_yyyyyz_1, g_0_y_0_yyyyzz_0, g_0_y_0_yyyyzz_1, g_0_y_0_yyyzzz_0, g_0_y_0_yyyzzz_1, g_0_y_0_yyzzzz_0, g_0_y_0_yyzzzz_1, g_0_y_0_yzzzzz_0, g_0_y_0_yzzzzz_1, g_0_y_0_zzzzzz_0, g_0_y_0_zzzzzz_1, g_0_yy_0_xxxxx_1, g_0_yy_0_xxxxxx_0, g_0_yy_0_xxxxxx_1, g_0_yy_0_xxxxxy_0, g_0_yy_0_xxxxxy_1, g_0_yy_0_xxxxxz_0, g_0_yy_0_xxxxxz_1, g_0_yy_0_xxxxy_1, g_0_yy_0_xxxxyy_0, g_0_yy_0_xxxxyy_1, g_0_yy_0_xxxxyz_0, g_0_yy_0_xxxxyz_1, g_0_yy_0_xxxxz_1, g_0_yy_0_xxxxzz_0, g_0_yy_0_xxxxzz_1, g_0_yy_0_xxxyy_1, g_0_yy_0_xxxyyy_0, g_0_yy_0_xxxyyy_1, g_0_yy_0_xxxyyz_0, g_0_yy_0_xxxyyz_1, g_0_yy_0_xxxyz_1, g_0_yy_0_xxxyzz_0, g_0_yy_0_xxxyzz_1, g_0_yy_0_xxxzz_1, g_0_yy_0_xxxzzz_0, g_0_yy_0_xxxzzz_1, g_0_yy_0_xxyyy_1, g_0_yy_0_xxyyyy_0, g_0_yy_0_xxyyyy_1, g_0_yy_0_xxyyyz_0, g_0_yy_0_xxyyyz_1, g_0_yy_0_xxyyz_1, g_0_yy_0_xxyyzz_0, g_0_yy_0_xxyyzz_1, g_0_yy_0_xxyzz_1, g_0_yy_0_xxyzzz_0, g_0_yy_0_xxyzzz_1, g_0_yy_0_xxzzz_1, g_0_yy_0_xxzzzz_0, g_0_yy_0_xxzzzz_1, g_0_yy_0_xyyyy_1, g_0_yy_0_xyyyyy_0, g_0_yy_0_xyyyyy_1, g_0_yy_0_xyyyyz_0, g_0_yy_0_xyyyyz_1, g_0_yy_0_xyyyz_1, g_0_yy_0_xyyyzz_0, g_0_yy_0_xyyyzz_1, g_0_yy_0_xyyzz_1, g_0_yy_0_xyyzzz_0, g_0_yy_0_xyyzzz_1, g_0_yy_0_xyzzz_1, g_0_yy_0_xyzzzz_0, g_0_yy_0_xyzzzz_1, g_0_yy_0_xzzzz_1, g_0_yy_0_xzzzzz_0, g_0_yy_0_xzzzzz_1, g_0_yy_0_yyyyy_1, g_0_yy_0_yyyyyy_0, g_0_yy_0_yyyyyy_1, g_0_yy_0_yyyyyz_0, g_0_yy_0_yyyyyz_1, g_0_yy_0_yyyyz_1, g_0_yy_0_yyyyzz_0, g_0_yy_0_yyyyzz_1, g_0_yy_0_yyyzz_1, g_0_yy_0_yyyzzz_0, g_0_yy_0_yyyzzz_1, g_0_yy_0_yyzzz_1, g_0_yy_0_yyzzzz_0, g_0_yy_0_yyzzzz_1, g_0_yy_0_yzzzz_1, g_0_yy_0_yzzzzz_0, g_0_yy_0_yzzzzz_1, g_0_yy_0_zzzzz_1, g_0_yy_0_zzzzzz_0, g_0_yy_0_zzzzzz_1, g_0_yyy_0_xxxxxx_0, g_0_yyy_0_xxxxxy_0, g_0_yyy_0_xxxxxz_0, g_0_yyy_0_xxxxyy_0, g_0_yyy_0_xxxxyz_0, g_0_yyy_0_xxxxzz_0, g_0_yyy_0_xxxyyy_0, g_0_yyy_0_xxxyyz_0, g_0_yyy_0_xxxyzz_0, g_0_yyy_0_xxxzzz_0, g_0_yyy_0_xxyyyy_0, g_0_yyy_0_xxyyyz_0, g_0_yyy_0_xxyyzz_0, g_0_yyy_0_xxyzzz_0, g_0_yyy_0_xxzzzz_0, g_0_yyy_0_xyyyyy_0, g_0_yyy_0_xyyyyz_0, g_0_yyy_0_xyyyzz_0, g_0_yyy_0_xyyzzz_0, g_0_yyy_0_xyzzzz_0, g_0_yyy_0_xzzzzz_0, g_0_yyy_0_yyyyyy_0, g_0_yyy_0_yyyyyz_0, g_0_yyy_0_yyyyzz_0, g_0_yyy_0_yyyzzz_0, g_0_yyy_0_yyzzzz_0, g_0_yyy_0_yzzzzz_0, g_0_yyy_0_zzzzzz_0, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyy_0_xxxxxx_0[i] = 2.0 * g_0_y_0_xxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxxxx_1[i] * fti_ab_0 + g_0_yy_0_xxxxxx_0[i] * pb_y + g_0_yy_0_xxxxxx_1[i] * wp_y[i];

        g_0_yyy_0_xxxxxy_0[i] = 2.0 * g_0_y_0_xxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxxxy_1[i] * fti_ab_0 + g_0_yy_0_xxxxx_1[i] * fi_abcd_0 + g_0_yy_0_xxxxxy_0[i] * pb_y + g_0_yy_0_xxxxxy_1[i] * wp_y[i];

        g_0_yyy_0_xxxxxz_0[i] = 2.0 * g_0_y_0_xxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxxxz_1[i] * fti_ab_0 + g_0_yy_0_xxxxxz_0[i] * pb_y + g_0_yy_0_xxxxxz_1[i] * wp_y[i];

        g_0_yyy_0_xxxxyy_0[i] = 2.0 * g_0_y_0_xxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxxyy_1[i] * fti_ab_0 + 2.0 * g_0_yy_0_xxxxy_1[i] * fi_abcd_0 + g_0_yy_0_xxxxyy_0[i] * pb_y + g_0_yy_0_xxxxyy_1[i] * wp_y[i];

        g_0_yyy_0_xxxxyz_0[i] = 2.0 * g_0_y_0_xxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxxyz_1[i] * fti_ab_0 + g_0_yy_0_xxxxz_1[i] * fi_abcd_0 + g_0_yy_0_xxxxyz_0[i] * pb_y + g_0_yy_0_xxxxyz_1[i] * wp_y[i];

        g_0_yyy_0_xxxxzz_0[i] = 2.0 * g_0_y_0_xxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxxzz_1[i] * fti_ab_0 + g_0_yy_0_xxxxzz_0[i] * pb_y + g_0_yy_0_xxxxzz_1[i] * wp_y[i];

        g_0_yyy_0_xxxyyy_0[i] = 2.0 * g_0_y_0_xxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxyyy_1[i] * fti_ab_0 + 3.0 * g_0_yy_0_xxxyy_1[i] * fi_abcd_0 + g_0_yy_0_xxxyyy_0[i] * pb_y + g_0_yy_0_xxxyyy_1[i] * wp_y[i];

        g_0_yyy_0_xxxyyz_0[i] = 2.0 * g_0_y_0_xxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yy_0_xxxyz_1[i] * fi_abcd_0 + g_0_yy_0_xxxyyz_0[i] * pb_y + g_0_yy_0_xxxyyz_1[i] * wp_y[i];

        g_0_yyy_0_xxxyzz_0[i] = 2.0 * g_0_y_0_xxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxyzz_1[i] * fti_ab_0 + g_0_yy_0_xxxzz_1[i] * fi_abcd_0 + g_0_yy_0_xxxyzz_0[i] * pb_y + g_0_yy_0_xxxyzz_1[i] * wp_y[i];

        g_0_yyy_0_xxxzzz_0[i] = 2.0 * g_0_y_0_xxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxzzz_1[i] * fti_ab_0 + g_0_yy_0_xxxzzz_0[i] * pb_y + g_0_yy_0_xxxzzz_1[i] * wp_y[i];

        g_0_yyy_0_xxyyyy_0[i] = 2.0 * g_0_y_0_xxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxyyyy_1[i] * fti_ab_0 + 4.0 * g_0_yy_0_xxyyy_1[i] * fi_abcd_0 + g_0_yy_0_xxyyyy_0[i] * pb_y + g_0_yy_0_xxyyyy_1[i] * wp_y[i];

        g_0_yyy_0_xxyyyz_0[i] = 2.0 * g_0_y_0_xxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yy_0_xxyyz_1[i] * fi_abcd_0 + g_0_yy_0_xxyyyz_0[i] * pb_y + g_0_yy_0_xxyyyz_1[i] * wp_y[i];

        g_0_yyy_0_xxyyzz_0[i] = 2.0 * g_0_y_0_xxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yy_0_xxyzz_1[i] * fi_abcd_0 + g_0_yy_0_xxyyzz_0[i] * pb_y + g_0_yy_0_xxyyzz_1[i] * wp_y[i];

        g_0_yyy_0_xxyzzz_0[i] = 2.0 * g_0_y_0_xxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxyzzz_1[i] * fti_ab_0 + g_0_yy_0_xxzzz_1[i] * fi_abcd_0 + g_0_yy_0_xxyzzz_0[i] * pb_y + g_0_yy_0_xxyzzz_1[i] * wp_y[i];

        g_0_yyy_0_xxzzzz_0[i] = 2.0 * g_0_y_0_xxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxzzzz_1[i] * fti_ab_0 + g_0_yy_0_xxzzzz_0[i] * pb_y + g_0_yy_0_xxzzzz_1[i] * wp_y[i];

        g_0_yyy_0_xyyyyy_0[i] = 2.0 * g_0_y_0_xyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xyyyyy_1[i] * fti_ab_0 + 5.0 * g_0_yy_0_xyyyy_1[i] * fi_abcd_0 + g_0_yy_0_xyyyyy_0[i] * pb_y + g_0_yy_0_xyyyyy_1[i] * wp_y[i];

        g_0_yyy_0_xyyyyz_0[i] = 2.0 * g_0_y_0_xyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xyyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yy_0_xyyyz_1[i] * fi_abcd_0 + g_0_yy_0_xyyyyz_0[i] * pb_y + g_0_yy_0_xyyyyz_1[i] * wp_y[i];

        g_0_yyy_0_xyyyzz_0[i] = 2.0 * g_0_y_0_xyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xyyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yy_0_xyyzz_1[i] * fi_abcd_0 + g_0_yy_0_xyyyzz_0[i] * pb_y + g_0_yy_0_xyyyzz_1[i] * wp_y[i];

        g_0_yyy_0_xyyzzz_0[i] = 2.0 * g_0_y_0_xyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yy_0_xyzzz_1[i] * fi_abcd_0 + g_0_yy_0_xyyzzz_0[i] * pb_y + g_0_yy_0_xyyzzz_1[i] * wp_y[i];

        g_0_yyy_0_xyzzzz_0[i] = 2.0 * g_0_y_0_xyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xyzzzz_1[i] * fti_ab_0 + g_0_yy_0_xzzzz_1[i] * fi_abcd_0 + g_0_yy_0_xyzzzz_0[i] * pb_y + g_0_yy_0_xyzzzz_1[i] * wp_y[i];

        g_0_yyy_0_xzzzzz_0[i] = 2.0 * g_0_y_0_xzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xzzzzz_1[i] * fti_ab_0 + g_0_yy_0_xzzzzz_0[i] * pb_y + g_0_yy_0_xzzzzz_1[i] * wp_y[i];

        g_0_yyy_0_yyyyyy_0[i] = 2.0 * g_0_y_0_yyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_y_0_yyyyyy_1[i] * fti_ab_0 + 6.0 * g_0_yy_0_yyyyy_1[i] * fi_abcd_0 + g_0_yy_0_yyyyyy_0[i] * pb_y + g_0_yy_0_yyyyyy_1[i] * wp_y[i];

        g_0_yyy_0_yyyyyz_0[i] = 2.0 * g_0_y_0_yyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_yyyyyz_1[i] * fti_ab_0 + 5.0 * g_0_yy_0_yyyyz_1[i] * fi_abcd_0 + g_0_yy_0_yyyyyz_0[i] * pb_y + g_0_yy_0_yyyyyz_1[i] * wp_y[i];

        g_0_yyy_0_yyyyzz_0[i] = 2.0 * g_0_y_0_yyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_yyyyzz_1[i] * fti_ab_0 + 4.0 * g_0_yy_0_yyyzz_1[i] * fi_abcd_0 + g_0_yy_0_yyyyzz_0[i] * pb_y + g_0_yy_0_yyyyzz_1[i] * wp_y[i];

        g_0_yyy_0_yyyzzz_0[i] = 2.0 * g_0_y_0_yyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_yyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_yy_0_yyzzz_1[i] * fi_abcd_0 + g_0_yy_0_yyyzzz_0[i] * pb_y + g_0_yy_0_yyyzzz_1[i] * wp_y[i];

        g_0_yyy_0_yyzzzz_0[i] = 2.0 * g_0_y_0_yyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_yyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yy_0_yzzzz_1[i] * fi_abcd_0 + g_0_yy_0_yyzzzz_0[i] * pb_y + g_0_yy_0_yyzzzz_1[i] * wp_y[i];

        g_0_yyy_0_yzzzzz_0[i] = 2.0 * g_0_y_0_yzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_yzzzzz_1[i] * fti_ab_0 + g_0_yy_0_zzzzz_1[i] * fi_abcd_0 + g_0_yy_0_yzzzzz_0[i] * pb_y + g_0_yy_0_yzzzzz_1[i] * wp_y[i];

        g_0_yyy_0_zzzzzz_0[i] = 2.0 * g_0_y_0_zzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_zzzzzz_1[i] * fti_ab_0 + g_0_yy_0_zzzzzz_0[i] * pb_y + g_0_yy_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 196-224 components of targeted buffer : SFSI

    auto g_0_yyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sfsi + 196);

    auto g_0_yyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sfsi + 197);

    auto g_0_yyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sfsi + 198);

    auto g_0_yyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sfsi + 199);

    auto g_0_yyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sfsi + 200);

    auto g_0_yyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sfsi + 201);

    auto g_0_yyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sfsi + 202);

    auto g_0_yyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sfsi + 203);

    auto g_0_yyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sfsi + 204);

    auto g_0_yyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sfsi + 205);

    auto g_0_yyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 206);

    auto g_0_yyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 207);

    auto g_0_yyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 208);

    auto g_0_yyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 209);

    auto g_0_yyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 210);

    auto g_0_yyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 211);

    auto g_0_yyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 212);

    auto g_0_yyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 213);

    auto g_0_yyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 214);

    auto g_0_yyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 215);

    auto g_0_yyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 216);

    auto g_0_yyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 217);

    auto g_0_yyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 218);

    auto g_0_yyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 219);

    auto g_0_yyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 220);

    auto g_0_yyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 221);

    auto g_0_yyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 222);

    auto g_0_yyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 223);

    #pragma omp simd aligned(g_0_yy_0_xxxxx_1, g_0_yy_0_xxxxxx_0, g_0_yy_0_xxxxxx_1, g_0_yy_0_xxxxxy_0, g_0_yy_0_xxxxxy_1, g_0_yy_0_xxxxxz_0, g_0_yy_0_xxxxxz_1, g_0_yy_0_xxxxy_1, g_0_yy_0_xxxxyy_0, g_0_yy_0_xxxxyy_1, g_0_yy_0_xxxxyz_0, g_0_yy_0_xxxxyz_1, g_0_yy_0_xxxxz_1, g_0_yy_0_xxxxzz_0, g_0_yy_0_xxxxzz_1, g_0_yy_0_xxxyy_1, g_0_yy_0_xxxyyy_0, g_0_yy_0_xxxyyy_1, g_0_yy_0_xxxyyz_0, g_0_yy_0_xxxyyz_1, g_0_yy_0_xxxyz_1, g_0_yy_0_xxxyzz_0, g_0_yy_0_xxxyzz_1, g_0_yy_0_xxxzz_1, g_0_yy_0_xxxzzz_0, g_0_yy_0_xxxzzz_1, g_0_yy_0_xxyyy_1, g_0_yy_0_xxyyyy_0, g_0_yy_0_xxyyyy_1, g_0_yy_0_xxyyyz_0, g_0_yy_0_xxyyyz_1, g_0_yy_0_xxyyz_1, g_0_yy_0_xxyyzz_0, g_0_yy_0_xxyyzz_1, g_0_yy_0_xxyzz_1, g_0_yy_0_xxyzzz_0, g_0_yy_0_xxyzzz_1, g_0_yy_0_xxzzz_1, g_0_yy_0_xxzzzz_0, g_0_yy_0_xxzzzz_1, g_0_yy_0_xyyyy_1, g_0_yy_0_xyyyyy_0, g_0_yy_0_xyyyyy_1, g_0_yy_0_xyyyyz_0, g_0_yy_0_xyyyyz_1, g_0_yy_0_xyyyz_1, g_0_yy_0_xyyyzz_0, g_0_yy_0_xyyyzz_1, g_0_yy_0_xyyzz_1, g_0_yy_0_xyyzzz_0, g_0_yy_0_xyyzzz_1, g_0_yy_0_xyzzz_1, g_0_yy_0_xyzzzz_0, g_0_yy_0_xyzzzz_1, g_0_yy_0_xzzzz_1, g_0_yy_0_xzzzzz_0, g_0_yy_0_xzzzzz_1, g_0_yy_0_yyyyy_1, g_0_yy_0_yyyyyy_0, g_0_yy_0_yyyyyy_1, g_0_yy_0_yyyyyz_0, g_0_yy_0_yyyyyz_1, g_0_yy_0_yyyyz_1, g_0_yy_0_yyyyzz_0, g_0_yy_0_yyyyzz_1, g_0_yy_0_yyyzz_1, g_0_yy_0_yyyzzz_0, g_0_yy_0_yyyzzz_1, g_0_yy_0_yyzzz_1, g_0_yy_0_yyzzzz_0, g_0_yy_0_yyzzzz_1, g_0_yy_0_yzzzz_1, g_0_yy_0_yzzzzz_0, g_0_yy_0_yzzzzz_1, g_0_yy_0_zzzzz_1, g_0_yy_0_zzzzzz_0, g_0_yy_0_zzzzzz_1, g_0_yyz_0_xxxxxx_0, g_0_yyz_0_xxxxxy_0, g_0_yyz_0_xxxxxz_0, g_0_yyz_0_xxxxyy_0, g_0_yyz_0_xxxxyz_0, g_0_yyz_0_xxxxzz_0, g_0_yyz_0_xxxyyy_0, g_0_yyz_0_xxxyyz_0, g_0_yyz_0_xxxyzz_0, g_0_yyz_0_xxxzzz_0, g_0_yyz_0_xxyyyy_0, g_0_yyz_0_xxyyyz_0, g_0_yyz_0_xxyyzz_0, g_0_yyz_0_xxyzzz_0, g_0_yyz_0_xxzzzz_0, g_0_yyz_0_xyyyyy_0, g_0_yyz_0_xyyyyz_0, g_0_yyz_0_xyyyzz_0, g_0_yyz_0_xyyzzz_0, g_0_yyz_0_xyzzzz_0, g_0_yyz_0_xzzzzz_0, g_0_yyz_0_yyyyyy_0, g_0_yyz_0_yyyyyz_0, g_0_yyz_0_yyyyzz_0, g_0_yyz_0_yyyzzz_0, g_0_yyz_0_yyzzzz_0, g_0_yyz_0_yzzzzz_0, g_0_yyz_0_zzzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyz_0_xxxxxx_0[i] = g_0_yy_0_xxxxxx_0[i] * pb_z + g_0_yy_0_xxxxxx_1[i] * wp_z[i];

        g_0_yyz_0_xxxxxy_0[i] = g_0_yy_0_xxxxxy_0[i] * pb_z + g_0_yy_0_xxxxxy_1[i] * wp_z[i];

        g_0_yyz_0_xxxxxz_0[i] = g_0_yy_0_xxxxx_1[i] * fi_abcd_0 + g_0_yy_0_xxxxxz_0[i] * pb_z + g_0_yy_0_xxxxxz_1[i] * wp_z[i];

        g_0_yyz_0_xxxxyy_0[i] = g_0_yy_0_xxxxyy_0[i] * pb_z + g_0_yy_0_xxxxyy_1[i] * wp_z[i];

        g_0_yyz_0_xxxxyz_0[i] = g_0_yy_0_xxxxy_1[i] * fi_abcd_0 + g_0_yy_0_xxxxyz_0[i] * pb_z + g_0_yy_0_xxxxyz_1[i] * wp_z[i];

        g_0_yyz_0_xxxxzz_0[i] = 2.0 * g_0_yy_0_xxxxz_1[i] * fi_abcd_0 + g_0_yy_0_xxxxzz_0[i] * pb_z + g_0_yy_0_xxxxzz_1[i] * wp_z[i];

        g_0_yyz_0_xxxyyy_0[i] = g_0_yy_0_xxxyyy_0[i] * pb_z + g_0_yy_0_xxxyyy_1[i] * wp_z[i];

        g_0_yyz_0_xxxyyz_0[i] = g_0_yy_0_xxxyy_1[i] * fi_abcd_0 + g_0_yy_0_xxxyyz_0[i] * pb_z + g_0_yy_0_xxxyyz_1[i] * wp_z[i];

        g_0_yyz_0_xxxyzz_0[i] = 2.0 * g_0_yy_0_xxxyz_1[i] * fi_abcd_0 + g_0_yy_0_xxxyzz_0[i] * pb_z + g_0_yy_0_xxxyzz_1[i] * wp_z[i];

        g_0_yyz_0_xxxzzz_0[i] = 3.0 * g_0_yy_0_xxxzz_1[i] * fi_abcd_0 + g_0_yy_0_xxxzzz_0[i] * pb_z + g_0_yy_0_xxxzzz_1[i] * wp_z[i];

        g_0_yyz_0_xxyyyy_0[i] = g_0_yy_0_xxyyyy_0[i] * pb_z + g_0_yy_0_xxyyyy_1[i] * wp_z[i];

        g_0_yyz_0_xxyyyz_0[i] = g_0_yy_0_xxyyy_1[i] * fi_abcd_0 + g_0_yy_0_xxyyyz_0[i] * pb_z + g_0_yy_0_xxyyyz_1[i] * wp_z[i];

        g_0_yyz_0_xxyyzz_0[i] = 2.0 * g_0_yy_0_xxyyz_1[i] * fi_abcd_0 + g_0_yy_0_xxyyzz_0[i] * pb_z + g_0_yy_0_xxyyzz_1[i] * wp_z[i];

        g_0_yyz_0_xxyzzz_0[i] = 3.0 * g_0_yy_0_xxyzz_1[i] * fi_abcd_0 + g_0_yy_0_xxyzzz_0[i] * pb_z + g_0_yy_0_xxyzzz_1[i] * wp_z[i];

        g_0_yyz_0_xxzzzz_0[i] = 4.0 * g_0_yy_0_xxzzz_1[i] * fi_abcd_0 + g_0_yy_0_xxzzzz_0[i] * pb_z + g_0_yy_0_xxzzzz_1[i] * wp_z[i];

        g_0_yyz_0_xyyyyy_0[i] = g_0_yy_0_xyyyyy_0[i] * pb_z + g_0_yy_0_xyyyyy_1[i] * wp_z[i];

        g_0_yyz_0_xyyyyz_0[i] = g_0_yy_0_xyyyy_1[i] * fi_abcd_0 + g_0_yy_0_xyyyyz_0[i] * pb_z + g_0_yy_0_xyyyyz_1[i] * wp_z[i];

        g_0_yyz_0_xyyyzz_0[i] = 2.0 * g_0_yy_0_xyyyz_1[i] * fi_abcd_0 + g_0_yy_0_xyyyzz_0[i] * pb_z + g_0_yy_0_xyyyzz_1[i] * wp_z[i];

        g_0_yyz_0_xyyzzz_0[i] = 3.0 * g_0_yy_0_xyyzz_1[i] * fi_abcd_0 + g_0_yy_0_xyyzzz_0[i] * pb_z + g_0_yy_0_xyyzzz_1[i] * wp_z[i];

        g_0_yyz_0_xyzzzz_0[i] = 4.0 * g_0_yy_0_xyzzz_1[i] * fi_abcd_0 + g_0_yy_0_xyzzzz_0[i] * pb_z + g_0_yy_0_xyzzzz_1[i] * wp_z[i];

        g_0_yyz_0_xzzzzz_0[i] = 5.0 * g_0_yy_0_xzzzz_1[i] * fi_abcd_0 + g_0_yy_0_xzzzzz_0[i] * pb_z + g_0_yy_0_xzzzzz_1[i] * wp_z[i];

        g_0_yyz_0_yyyyyy_0[i] = g_0_yy_0_yyyyyy_0[i] * pb_z + g_0_yy_0_yyyyyy_1[i] * wp_z[i];

        g_0_yyz_0_yyyyyz_0[i] = g_0_yy_0_yyyyy_1[i] * fi_abcd_0 + g_0_yy_0_yyyyyz_0[i] * pb_z + g_0_yy_0_yyyyyz_1[i] * wp_z[i];

        g_0_yyz_0_yyyyzz_0[i] = 2.0 * g_0_yy_0_yyyyz_1[i] * fi_abcd_0 + g_0_yy_0_yyyyzz_0[i] * pb_z + g_0_yy_0_yyyyzz_1[i] * wp_z[i];

        g_0_yyz_0_yyyzzz_0[i] = 3.0 * g_0_yy_0_yyyzz_1[i] * fi_abcd_0 + g_0_yy_0_yyyzzz_0[i] * pb_z + g_0_yy_0_yyyzzz_1[i] * wp_z[i];

        g_0_yyz_0_yyzzzz_0[i] = 4.0 * g_0_yy_0_yyzzz_1[i] * fi_abcd_0 + g_0_yy_0_yyzzzz_0[i] * pb_z + g_0_yy_0_yyzzzz_1[i] * wp_z[i];

        g_0_yyz_0_yzzzzz_0[i] = 5.0 * g_0_yy_0_yzzzz_1[i] * fi_abcd_0 + g_0_yy_0_yzzzzz_0[i] * pb_z + g_0_yy_0_yzzzzz_1[i] * wp_z[i];

        g_0_yyz_0_zzzzzz_0[i] = 6.0 * g_0_yy_0_zzzzz_1[i] * fi_abcd_0 + g_0_yy_0_zzzzzz_0[i] * pb_z + g_0_yy_0_zzzzzz_1[i] * wp_z[i];
    }

    /// Set up 224-252 components of targeted buffer : SFSI

    auto g_0_yzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sfsi + 224);

    auto g_0_yzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sfsi + 225);

    auto g_0_yzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sfsi + 226);

    auto g_0_yzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sfsi + 227);

    auto g_0_yzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sfsi + 228);

    auto g_0_yzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sfsi + 229);

    auto g_0_yzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sfsi + 230);

    auto g_0_yzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sfsi + 231);

    auto g_0_yzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sfsi + 232);

    auto g_0_yzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sfsi + 233);

    auto g_0_yzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 234);

    auto g_0_yzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 235);

    auto g_0_yzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 236);

    auto g_0_yzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 237);

    auto g_0_yzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 238);

    auto g_0_yzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 239);

    auto g_0_yzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 240);

    auto g_0_yzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 241);

    auto g_0_yzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 242);

    auto g_0_yzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 243);

    auto g_0_yzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 244);

    auto g_0_yzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 245);

    auto g_0_yzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 246);

    auto g_0_yzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 247);

    auto g_0_yzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 248);

    auto g_0_yzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 249);

    auto g_0_yzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 250);

    auto g_0_yzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 251);

    #pragma omp simd aligned(g_0_yzz_0_xxxxxx_0, g_0_yzz_0_xxxxxy_0, g_0_yzz_0_xxxxxz_0, g_0_yzz_0_xxxxyy_0, g_0_yzz_0_xxxxyz_0, g_0_yzz_0_xxxxzz_0, g_0_yzz_0_xxxyyy_0, g_0_yzz_0_xxxyyz_0, g_0_yzz_0_xxxyzz_0, g_0_yzz_0_xxxzzz_0, g_0_yzz_0_xxyyyy_0, g_0_yzz_0_xxyyyz_0, g_0_yzz_0_xxyyzz_0, g_0_yzz_0_xxyzzz_0, g_0_yzz_0_xxzzzz_0, g_0_yzz_0_xyyyyy_0, g_0_yzz_0_xyyyyz_0, g_0_yzz_0_xyyyzz_0, g_0_yzz_0_xyyzzz_0, g_0_yzz_0_xyzzzz_0, g_0_yzz_0_xzzzzz_0, g_0_yzz_0_yyyyyy_0, g_0_yzz_0_yyyyyz_0, g_0_yzz_0_yyyyzz_0, g_0_yzz_0_yyyzzz_0, g_0_yzz_0_yyzzzz_0, g_0_yzz_0_yzzzzz_0, g_0_yzz_0_zzzzzz_0, g_0_zz_0_xxxxx_1, g_0_zz_0_xxxxxx_0, g_0_zz_0_xxxxxx_1, g_0_zz_0_xxxxxy_0, g_0_zz_0_xxxxxy_1, g_0_zz_0_xxxxxz_0, g_0_zz_0_xxxxxz_1, g_0_zz_0_xxxxy_1, g_0_zz_0_xxxxyy_0, g_0_zz_0_xxxxyy_1, g_0_zz_0_xxxxyz_0, g_0_zz_0_xxxxyz_1, g_0_zz_0_xxxxz_1, g_0_zz_0_xxxxzz_0, g_0_zz_0_xxxxzz_1, g_0_zz_0_xxxyy_1, g_0_zz_0_xxxyyy_0, g_0_zz_0_xxxyyy_1, g_0_zz_0_xxxyyz_0, g_0_zz_0_xxxyyz_1, g_0_zz_0_xxxyz_1, g_0_zz_0_xxxyzz_0, g_0_zz_0_xxxyzz_1, g_0_zz_0_xxxzz_1, g_0_zz_0_xxxzzz_0, g_0_zz_0_xxxzzz_1, g_0_zz_0_xxyyy_1, g_0_zz_0_xxyyyy_0, g_0_zz_0_xxyyyy_1, g_0_zz_0_xxyyyz_0, g_0_zz_0_xxyyyz_1, g_0_zz_0_xxyyz_1, g_0_zz_0_xxyyzz_0, g_0_zz_0_xxyyzz_1, g_0_zz_0_xxyzz_1, g_0_zz_0_xxyzzz_0, g_0_zz_0_xxyzzz_1, g_0_zz_0_xxzzz_1, g_0_zz_0_xxzzzz_0, g_0_zz_0_xxzzzz_1, g_0_zz_0_xyyyy_1, g_0_zz_0_xyyyyy_0, g_0_zz_0_xyyyyy_1, g_0_zz_0_xyyyyz_0, g_0_zz_0_xyyyyz_1, g_0_zz_0_xyyyz_1, g_0_zz_0_xyyyzz_0, g_0_zz_0_xyyyzz_1, g_0_zz_0_xyyzz_1, g_0_zz_0_xyyzzz_0, g_0_zz_0_xyyzzz_1, g_0_zz_0_xyzzz_1, g_0_zz_0_xyzzzz_0, g_0_zz_0_xyzzzz_1, g_0_zz_0_xzzzz_1, g_0_zz_0_xzzzzz_0, g_0_zz_0_xzzzzz_1, g_0_zz_0_yyyyy_1, g_0_zz_0_yyyyyy_0, g_0_zz_0_yyyyyy_1, g_0_zz_0_yyyyyz_0, g_0_zz_0_yyyyyz_1, g_0_zz_0_yyyyz_1, g_0_zz_0_yyyyzz_0, g_0_zz_0_yyyyzz_1, g_0_zz_0_yyyzz_1, g_0_zz_0_yyyzzz_0, g_0_zz_0_yyyzzz_1, g_0_zz_0_yyzzz_1, g_0_zz_0_yyzzzz_0, g_0_zz_0_yyzzzz_1, g_0_zz_0_yzzzz_1, g_0_zz_0_yzzzzz_0, g_0_zz_0_yzzzzz_1, g_0_zz_0_zzzzz_1, g_0_zz_0_zzzzzz_0, g_0_zz_0_zzzzzz_1, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzz_0_xxxxxx_0[i] = g_0_zz_0_xxxxxx_0[i] * pb_y + g_0_zz_0_xxxxxx_1[i] * wp_y[i];

        g_0_yzz_0_xxxxxy_0[i] = g_0_zz_0_xxxxx_1[i] * fi_abcd_0 + g_0_zz_0_xxxxxy_0[i] * pb_y + g_0_zz_0_xxxxxy_1[i] * wp_y[i];

        g_0_yzz_0_xxxxxz_0[i] = g_0_zz_0_xxxxxz_0[i] * pb_y + g_0_zz_0_xxxxxz_1[i] * wp_y[i];

        g_0_yzz_0_xxxxyy_0[i] = 2.0 * g_0_zz_0_xxxxy_1[i] * fi_abcd_0 + g_0_zz_0_xxxxyy_0[i] * pb_y + g_0_zz_0_xxxxyy_1[i] * wp_y[i];

        g_0_yzz_0_xxxxyz_0[i] = g_0_zz_0_xxxxz_1[i] * fi_abcd_0 + g_0_zz_0_xxxxyz_0[i] * pb_y + g_0_zz_0_xxxxyz_1[i] * wp_y[i];

        g_0_yzz_0_xxxxzz_0[i] = g_0_zz_0_xxxxzz_0[i] * pb_y + g_0_zz_0_xxxxzz_1[i] * wp_y[i];

        g_0_yzz_0_xxxyyy_0[i] = 3.0 * g_0_zz_0_xxxyy_1[i] * fi_abcd_0 + g_0_zz_0_xxxyyy_0[i] * pb_y + g_0_zz_0_xxxyyy_1[i] * wp_y[i];

        g_0_yzz_0_xxxyyz_0[i] = 2.0 * g_0_zz_0_xxxyz_1[i] * fi_abcd_0 + g_0_zz_0_xxxyyz_0[i] * pb_y + g_0_zz_0_xxxyyz_1[i] * wp_y[i];

        g_0_yzz_0_xxxyzz_0[i] = g_0_zz_0_xxxzz_1[i] * fi_abcd_0 + g_0_zz_0_xxxyzz_0[i] * pb_y + g_0_zz_0_xxxyzz_1[i] * wp_y[i];

        g_0_yzz_0_xxxzzz_0[i] = g_0_zz_0_xxxzzz_0[i] * pb_y + g_0_zz_0_xxxzzz_1[i] * wp_y[i];

        g_0_yzz_0_xxyyyy_0[i] = 4.0 * g_0_zz_0_xxyyy_1[i] * fi_abcd_0 + g_0_zz_0_xxyyyy_0[i] * pb_y + g_0_zz_0_xxyyyy_1[i] * wp_y[i];

        g_0_yzz_0_xxyyyz_0[i] = 3.0 * g_0_zz_0_xxyyz_1[i] * fi_abcd_0 + g_0_zz_0_xxyyyz_0[i] * pb_y + g_0_zz_0_xxyyyz_1[i] * wp_y[i];

        g_0_yzz_0_xxyyzz_0[i] = 2.0 * g_0_zz_0_xxyzz_1[i] * fi_abcd_0 + g_0_zz_0_xxyyzz_0[i] * pb_y + g_0_zz_0_xxyyzz_1[i] * wp_y[i];

        g_0_yzz_0_xxyzzz_0[i] = g_0_zz_0_xxzzz_1[i] * fi_abcd_0 + g_0_zz_0_xxyzzz_0[i] * pb_y + g_0_zz_0_xxyzzz_1[i] * wp_y[i];

        g_0_yzz_0_xxzzzz_0[i] = g_0_zz_0_xxzzzz_0[i] * pb_y + g_0_zz_0_xxzzzz_1[i] * wp_y[i];

        g_0_yzz_0_xyyyyy_0[i] = 5.0 * g_0_zz_0_xyyyy_1[i] * fi_abcd_0 + g_0_zz_0_xyyyyy_0[i] * pb_y + g_0_zz_0_xyyyyy_1[i] * wp_y[i];

        g_0_yzz_0_xyyyyz_0[i] = 4.0 * g_0_zz_0_xyyyz_1[i] * fi_abcd_0 + g_0_zz_0_xyyyyz_0[i] * pb_y + g_0_zz_0_xyyyyz_1[i] * wp_y[i];

        g_0_yzz_0_xyyyzz_0[i] = 3.0 * g_0_zz_0_xyyzz_1[i] * fi_abcd_0 + g_0_zz_0_xyyyzz_0[i] * pb_y + g_0_zz_0_xyyyzz_1[i] * wp_y[i];

        g_0_yzz_0_xyyzzz_0[i] = 2.0 * g_0_zz_0_xyzzz_1[i] * fi_abcd_0 + g_0_zz_0_xyyzzz_0[i] * pb_y + g_0_zz_0_xyyzzz_1[i] * wp_y[i];

        g_0_yzz_0_xyzzzz_0[i] = g_0_zz_0_xzzzz_1[i] * fi_abcd_0 + g_0_zz_0_xyzzzz_0[i] * pb_y + g_0_zz_0_xyzzzz_1[i] * wp_y[i];

        g_0_yzz_0_xzzzzz_0[i] = g_0_zz_0_xzzzzz_0[i] * pb_y + g_0_zz_0_xzzzzz_1[i] * wp_y[i];

        g_0_yzz_0_yyyyyy_0[i] = 6.0 * g_0_zz_0_yyyyy_1[i] * fi_abcd_0 + g_0_zz_0_yyyyyy_0[i] * pb_y + g_0_zz_0_yyyyyy_1[i] * wp_y[i];

        g_0_yzz_0_yyyyyz_0[i] = 5.0 * g_0_zz_0_yyyyz_1[i] * fi_abcd_0 + g_0_zz_0_yyyyyz_0[i] * pb_y + g_0_zz_0_yyyyyz_1[i] * wp_y[i];

        g_0_yzz_0_yyyyzz_0[i] = 4.0 * g_0_zz_0_yyyzz_1[i] * fi_abcd_0 + g_0_zz_0_yyyyzz_0[i] * pb_y + g_0_zz_0_yyyyzz_1[i] * wp_y[i];

        g_0_yzz_0_yyyzzz_0[i] = 3.0 * g_0_zz_0_yyzzz_1[i] * fi_abcd_0 + g_0_zz_0_yyyzzz_0[i] * pb_y + g_0_zz_0_yyyzzz_1[i] * wp_y[i];

        g_0_yzz_0_yyzzzz_0[i] = 2.0 * g_0_zz_0_yzzzz_1[i] * fi_abcd_0 + g_0_zz_0_yyzzzz_0[i] * pb_y + g_0_zz_0_yyzzzz_1[i] * wp_y[i];

        g_0_yzz_0_yzzzzz_0[i] = g_0_zz_0_zzzzz_1[i] * fi_abcd_0 + g_0_zz_0_yzzzzz_0[i] * pb_y + g_0_zz_0_yzzzzz_1[i] * wp_y[i];

        g_0_yzz_0_zzzzzz_0[i] = g_0_zz_0_zzzzzz_0[i] * pb_y + g_0_zz_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 252-280 components of targeted buffer : SFSI

    auto g_0_zzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sfsi + 252);

    auto g_0_zzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sfsi + 253);

    auto g_0_zzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sfsi + 254);

    auto g_0_zzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sfsi + 255);

    auto g_0_zzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sfsi + 256);

    auto g_0_zzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sfsi + 257);

    auto g_0_zzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sfsi + 258);

    auto g_0_zzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sfsi + 259);

    auto g_0_zzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sfsi + 260);

    auto g_0_zzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sfsi + 261);

    auto g_0_zzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 262);

    auto g_0_zzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 263);

    auto g_0_zzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 264);

    auto g_0_zzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 265);

    auto g_0_zzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 266);

    auto g_0_zzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 267);

    auto g_0_zzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 268);

    auto g_0_zzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 269);

    auto g_0_zzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 270);

    auto g_0_zzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 271);

    auto g_0_zzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 272);

    auto g_0_zzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sfsi + 273);

    auto g_0_zzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sfsi + 274);

    auto g_0_zzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sfsi + 275);

    auto g_0_zzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sfsi + 276);

    auto g_0_zzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 277);

    auto g_0_zzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 278);

    auto g_0_zzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sfsi + 279);

    #pragma omp simd aligned(g_0_z_0_xxxxxx_0, g_0_z_0_xxxxxx_1, g_0_z_0_xxxxxy_0, g_0_z_0_xxxxxy_1, g_0_z_0_xxxxxz_0, g_0_z_0_xxxxxz_1, g_0_z_0_xxxxyy_0, g_0_z_0_xxxxyy_1, g_0_z_0_xxxxyz_0, g_0_z_0_xxxxyz_1, g_0_z_0_xxxxzz_0, g_0_z_0_xxxxzz_1, g_0_z_0_xxxyyy_0, g_0_z_0_xxxyyy_1, g_0_z_0_xxxyyz_0, g_0_z_0_xxxyyz_1, g_0_z_0_xxxyzz_0, g_0_z_0_xxxyzz_1, g_0_z_0_xxxzzz_0, g_0_z_0_xxxzzz_1, g_0_z_0_xxyyyy_0, g_0_z_0_xxyyyy_1, g_0_z_0_xxyyyz_0, g_0_z_0_xxyyyz_1, g_0_z_0_xxyyzz_0, g_0_z_0_xxyyzz_1, g_0_z_0_xxyzzz_0, g_0_z_0_xxyzzz_1, g_0_z_0_xxzzzz_0, g_0_z_0_xxzzzz_1, g_0_z_0_xyyyyy_0, g_0_z_0_xyyyyy_1, g_0_z_0_xyyyyz_0, g_0_z_0_xyyyyz_1, g_0_z_0_xyyyzz_0, g_0_z_0_xyyyzz_1, g_0_z_0_xyyzzz_0, g_0_z_0_xyyzzz_1, g_0_z_0_xyzzzz_0, g_0_z_0_xyzzzz_1, g_0_z_0_xzzzzz_0, g_0_z_0_xzzzzz_1, g_0_z_0_yyyyyy_0, g_0_z_0_yyyyyy_1, g_0_z_0_yyyyyz_0, g_0_z_0_yyyyyz_1, g_0_z_0_yyyyzz_0, g_0_z_0_yyyyzz_1, g_0_z_0_yyyzzz_0, g_0_z_0_yyyzzz_1, g_0_z_0_yyzzzz_0, g_0_z_0_yyzzzz_1, g_0_z_0_yzzzzz_0, g_0_z_0_yzzzzz_1, g_0_z_0_zzzzzz_0, g_0_z_0_zzzzzz_1, g_0_zz_0_xxxxx_1, g_0_zz_0_xxxxxx_0, g_0_zz_0_xxxxxx_1, g_0_zz_0_xxxxxy_0, g_0_zz_0_xxxxxy_1, g_0_zz_0_xxxxxz_0, g_0_zz_0_xxxxxz_1, g_0_zz_0_xxxxy_1, g_0_zz_0_xxxxyy_0, g_0_zz_0_xxxxyy_1, g_0_zz_0_xxxxyz_0, g_0_zz_0_xxxxyz_1, g_0_zz_0_xxxxz_1, g_0_zz_0_xxxxzz_0, g_0_zz_0_xxxxzz_1, g_0_zz_0_xxxyy_1, g_0_zz_0_xxxyyy_0, g_0_zz_0_xxxyyy_1, g_0_zz_0_xxxyyz_0, g_0_zz_0_xxxyyz_1, g_0_zz_0_xxxyz_1, g_0_zz_0_xxxyzz_0, g_0_zz_0_xxxyzz_1, g_0_zz_0_xxxzz_1, g_0_zz_0_xxxzzz_0, g_0_zz_0_xxxzzz_1, g_0_zz_0_xxyyy_1, g_0_zz_0_xxyyyy_0, g_0_zz_0_xxyyyy_1, g_0_zz_0_xxyyyz_0, g_0_zz_0_xxyyyz_1, g_0_zz_0_xxyyz_1, g_0_zz_0_xxyyzz_0, g_0_zz_0_xxyyzz_1, g_0_zz_0_xxyzz_1, g_0_zz_0_xxyzzz_0, g_0_zz_0_xxyzzz_1, g_0_zz_0_xxzzz_1, g_0_zz_0_xxzzzz_0, g_0_zz_0_xxzzzz_1, g_0_zz_0_xyyyy_1, g_0_zz_0_xyyyyy_0, g_0_zz_0_xyyyyy_1, g_0_zz_0_xyyyyz_0, g_0_zz_0_xyyyyz_1, g_0_zz_0_xyyyz_1, g_0_zz_0_xyyyzz_0, g_0_zz_0_xyyyzz_1, g_0_zz_0_xyyzz_1, g_0_zz_0_xyyzzz_0, g_0_zz_0_xyyzzz_1, g_0_zz_0_xyzzz_1, g_0_zz_0_xyzzzz_0, g_0_zz_0_xyzzzz_1, g_0_zz_0_xzzzz_1, g_0_zz_0_xzzzzz_0, g_0_zz_0_xzzzzz_1, g_0_zz_0_yyyyy_1, g_0_zz_0_yyyyyy_0, g_0_zz_0_yyyyyy_1, g_0_zz_0_yyyyyz_0, g_0_zz_0_yyyyyz_1, g_0_zz_0_yyyyz_1, g_0_zz_0_yyyyzz_0, g_0_zz_0_yyyyzz_1, g_0_zz_0_yyyzz_1, g_0_zz_0_yyyzzz_0, g_0_zz_0_yyyzzz_1, g_0_zz_0_yyzzz_1, g_0_zz_0_yyzzzz_0, g_0_zz_0_yyzzzz_1, g_0_zz_0_yzzzz_1, g_0_zz_0_yzzzzz_0, g_0_zz_0_yzzzzz_1, g_0_zz_0_zzzzz_1, g_0_zz_0_zzzzzz_0, g_0_zz_0_zzzzzz_1, g_0_zzz_0_xxxxxx_0, g_0_zzz_0_xxxxxy_0, g_0_zzz_0_xxxxxz_0, g_0_zzz_0_xxxxyy_0, g_0_zzz_0_xxxxyz_0, g_0_zzz_0_xxxxzz_0, g_0_zzz_0_xxxyyy_0, g_0_zzz_0_xxxyyz_0, g_0_zzz_0_xxxyzz_0, g_0_zzz_0_xxxzzz_0, g_0_zzz_0_xxyyyy_0, g_0_zzz_0_xxyyyz_0, g_0_zzz_0_xxyyzz_0, g_0_zzz_0_xxyzzz_0, g_0_zzz_0_xxzzzz_0, g_0_zzz_0_xyyyyy_0, g_0_zzz_0_xyyyyz_0, g_0_zzz_0_xyyyzz_0, g_0_zzz_0_xyyzzz_0, g_0_zzz_0_xyzzzz_0, g_0_zzz_0_xzzzzz_0, g_0_zzz_0_yyyyyy_0, g_0_zzz_0_yyyyyz_0, g_0_zzz_0_yyyyzz_0, g_0_zzz_0_yyyzzz_0, g_0_zzz_0_yyzzzz_0, g_0_zzz_0_yzzzzz_0, g_0_zzz_0_zzzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzz_0_xxxxxx_0[i] = 2.0 * g_0_z_0_xxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxxxx_1[i] * fti_ab_0 + g_0_zz_0_xxxxxx_0[i] * pb_z + g_0_zz_0_xxxxxx_1[i] * wp_z[i];

        g_0_zzz_0_xxxxxy_0[i] = 2.0 * g_0_z_0_xxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxxxy_1[i] * fti_ab_0 + g_0_zz_0_xxxxxy_0[i] * pb_z + g_0_zz_0_xxxxxy_1[i] * wp_z[i];

        g_0_zzz_0_xxxxxz_0[i] = 2.0 * g_0_z_0_xxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxxxz_1[i] * fti_ab_0 + g_0_zz_0_xxxxx_1[i] * fi_abcd_0 + g_0_zz_0_xxxxxz_0[i] * pb_z + g_0_zz_0_xxxxxz_1[i] * wp_z[i];

        g_0_zzz_0_xxxxyy_0[i] = 2.0 * g_0_z_0_xxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxxyy_1[i] * fti_ab_0 + g_0_zz_0_xxxxyy_0[i] * pb_z + g_0_zz_0_xxxxyy_1[i] * wp_z[i];

        g_0_zzz_0_xxxxyz_0[i] = 2.0 * g_0_z_0_xxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxxyz_1[i] * fti_ab_0 + g_0_zz_0_xxxxy_1[i] * fi_abcd_0 + g_0_zz_0_xxxxyz_0[i] * pb_z + g_0_zz_0_xxxxyz_1[i] * wp_z[i];

        g_0_zzz_0_xxxxzz_0[i] = 2.0 * g_0_z_0_xxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxxzz_1[i] * fti_ab_0 + 2.0 * g_0_zz_0_xxxxz_1[i] * fi_abcd_0 + g_0_zz_0_xxxxzz_0[i] * pb_z + g_0_zz_0_xxxxzz_1[i] * wp_z[i];

        g_0_zzz_0_xxxyyy_0[i] = 2.0 * g_0_z_0_xxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxyyy_1[i] * fti_ab_0 + g_0_zz_0_xxxyyy_0[i] * pb_z + g_0_zz_0_xxxyyy_1[i] * wp_z[i];

        g_0_zzz_0_xxxyyz_0[i] = 2.0 * g_0_z_0_xxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxyyz_1[i] * fti_ab_0 + g_0_zz_0_xxxyy_1[i] * fi_abcd_0 + g_0_zz_0_xxxyyz_0[i] * pb_z + g_0_zz_0_xxxyyz_1[i] * wp_z[i];

        g_0_zzz_0_xxxyzz_0[i] = 2.0 * g_0_z_0_xxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxyzz_1[i] * fti_ab_0 + 2.0 * g_0_zz_0_xxxyz_1[i] * fi_abcd_0 + g_0_zz_0_xxxyzz_0[i] * pb_z + g_0_zz_0_xxxyzz_1[i] * wp_z[i];

        g_0_zzz_0_xxxzzz_0[i] = 2.0 * g_0_z_0_xxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxzzz_1[i] * fti_ab_0 + 3.0 * g_0_zz_0_xxxzz_1[i] * fi_abcd_0 + g_0_zz_0_xxxzzz_0[i] * pb_z + g_0_zz_0_xxxzzz_1[i] * wp_z[i];

        g_0_zzz_0_xxyyyy_0[i] = 2.0 * g_0_z_0_xxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxyyyy_1[i] * fti_ab_0 + g_0_zz_0_xxyyyy_0[i] * pb_z + g_0_zz_0_xxyyyy_1[i] * wp_z[i];

        g_0_zzz_0_xxyyyz_0[i] = 2.0 * g_0_z_0_xxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxyyyz_1[i] * fti_ab_0 + g_0_zz_0_xxyyy_1[i] * fi_abcd_0 + g_0_zz_0_xxyyyz_0[i] * pb_z + g_0_zz_0_xxyyyz_1[i] * wp_z[i];

        g_0_zzz_0_xxyyzz_0[i] = 2.0 * g_0_z_0_xxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zz_0_xxyyz_1[i] * fi_abcd_0 + g_0_zz_0_xxyyzz_0[i] * pb_z + g_0_zz_0_xxyyzz_1[i] * wp_z[i];

        g_0_zzz_0_xxyzzz_0[i] = 2.0 * g_0_z_0_xxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zz_0_xxyzz_1[i] * fi_abcd_0 + g_0_zz_0_xxyzzz_0[i] * pb_z + g_0_zz_0_xxyzzz_1[i] * wp_z[i];

        g_0_zzz_0_xxzzzz_0[i] = 2.0 * g_0_z_0_xxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zz_0_xxzzz_1[i] * fi_abcd_0 + g_0_zz_0_xxzzzz_0[i] * pb_z + g_0_zz_0_xxzzzz_1[i] * wp_z[i];

        g_0_zzz_0_xyyyyy_0[i] = 2.0 * g_0_z_0_xyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xyyyyy_1[i] * fti_ab_0 + g_0_zz_0_xyyyyy_0[i] * pb_z + g_0_zz_0_xyyyyy_1[i] * wp_z[i];

        g_0_zzz_0_xyyyyz_0[i] = 2.0 * g_0_z_0_xyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xyyyyz_1[i] * fti_ab_0 + g_0_zz_0_xyyyy_1[i] * fi_abcd_0 + g_0_zz_0_xyyyyz_0[i] * pb_z + g_0_zz_0_xyyyyz_1[i] * wp_z[i];

        g_0_zzz_0_xyyyzz_0[i] = 2.0 * g_0_z_0_xyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zz_0_xyyyz_1[i] * fi_abcd_0 + g_0_zz_0_xyyyzz_0[i] * pb_z + g_0_zz_0_xyyyzz_1[i] * wp_z[i];

        g_0_zzz_0_xyyzzz_0[i] = 2.0 * g_0_z_0_xyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zz_0_xyyzz_1[i] * fi_abcd_0 + g_0_zz_0_xyyzzz_0[i] * pb_z + g_0_zz_0_xyyzzz_1[i] * wp_z[i];

        g_0_zzz_0_xyzzzz_0[i] = 2.0 * g_0_z_0_xyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xyzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zz_0_xyzzz_1[i] * fi_abcd_0 + g_0_zz_0_xyzzzz_0[i] * pb_z + g_0_zz_0_xyzzzz_1[i] * wp_z[i];

        g_0_zzz_0_xzzzzz_0[i] = 2.0 * g_0_z_0_xzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xzzzzz_1[i] * fti_ab_0 + 5.0 * g_0_zz_0_xzzzz_1[i] * fi_abcd_0 + g_0_zz_0_xzzzzz_0[i] * pb_z + g_0_zz_0_xzzzzz_1[i] * wp_z[i];

        g_0_zzz_0_yyyyyy_0[i] = 2.0 * g_0_z_0_yyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_z_0_yyyyyy_1[i] * fti_ab_0 + g_0_zz_0_yyyyyy_0[i] * pb_z + g_0_zz_0_yyyyyy_1[i] * wp_z[i];

        g_0_zzz_0_yyyyyz_0[i] = 2.0 * g_0_z_0_yyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_yyyyyz_1[i] * fti_ab_0 + g_0_zz_0_yyyyy_1[i] * fi_abcd_0 + g_0_zz_0_yyyyyz_0[i] * pb_z + g_0_zz_0_yyyyyz_1[i] * wp_z[i];

        g_0_zzz_0_yyyyzz_0[i] = 2.0 * g_0_z_0_yyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_yyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zz_0_yyyyz_1[i] * fi_abcd_0 + g_0_zz_0_yyyyzz_0[i] * pb_z + g_0_zz_0_yyyyzz_1[i] * wp_z[i];

        g_0_zzz_0_yyyzzz_0[i] = 2.0 * g_0_z_0_yyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_yyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zz_0_yyyzz_1[i] * fi_abcd_0 + g_0_zz_0_yyyzzz_0[i] * pb_z + g_0_zz_0_yyyzzz_1[i] * wp_z[i];

        g_0_zzz_0_yyzzzz_0[i] = 2.0 * g_0_z_0_yyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_yyzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zz_0_yyzzz_1[i] * fi_abcd_0 + g_0_zz_0_yyzzzz_0[i] * pb_z + g_0_zz_0_yyzzzz_1[i] * wp_z[i];

        g_0_zzz_0_yzzzzz_0[i] = 2.0 * g_0_z_0_yzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_yzzzzz_1[i] * fti_ab_0 + 5.0 * g_0_zz_0_yzzzz_1[i] * fi_abcd_0 + g_0_zz_0_yzzzzz_0[i] * pb_z + g_0_zz_0_yzzzzz_1[i] * wp_z[i];

        g_0_zzz_0_zzzzzz_0[i] = 2.0 * g_0_z_0_zzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_zzzzzz_1[i] * fti_ab_0 + 6.0 * g_0_zz_0_zzzzz_1[i] * fi_abcd_0 + g_0_zz_0_zzzzzz_0[i] * pb_z + g_0_zz_0_zzzzzz_1[i] * wp_z[i];
    }
}

} // erirec namespace


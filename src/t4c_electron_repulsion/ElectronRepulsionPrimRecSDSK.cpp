#include "ElectronRepulsionPrimRecSDSK.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_sdsk(CSimdArray<double>& pbuffer,
                                  const size_t idx_eri_0_sdsk,
                                  size_t idx_eri_0_sssk,
                                  size_t idx_eri_1_sssk,
                                  size_t idx_eri_1_spsi,
                                  size_t idx_eri_0_spsk,
                                  size_t idx_eri_1_spsk,
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

    /// Set up components of auxilary buffer : SSSK

    auto g_0_0_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sssk);

    auto g_0_0_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sssk + 1);

    auto g_0_0_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sssk + 2);

    auto g_0_0_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sssk + 3);

    auto g_0_0_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sssk + 4);

    auto g_0_0_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sssk + 5);

    auto g_0_0_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sssk + 6);

    auto g_0_0_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sssk + 7);

    auto g_0_0_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sssk + 8);

    auto g_0_0_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sssk + 9);

    auto g_0_0_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sssk + 10);

    auto g_0_0_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sssk + 11);

    auto g_0_0_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sssk + 12);

    auto g_0_0_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sssk + 13);

    auto g_0_0_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sssk + 14);

    auto g_0_0_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sssk + 15);

    auto g_0_0_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sssk + 16);

    auto g_0_0_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sssk + 17);

    auto g_0_0_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sssk + 18);

    auto g_0_0_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sssk + 19);

    auto g_0_0_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sssk + 20);

    auto g_0_0_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sssk + 21);

    auto g_0_0_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sssk + 22);

    auto g_0_0_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sssk + 23);

    auto g_0_0_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sssk + 24);

    auto g_0_0_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sssk + 25);

    auto g_0_0_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sssk + 26);

    auto g_0_0_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sssk + 27);

    auto g_0_0_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sssk + 28);

    auto g_0_0_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sssk + 29);

    auto g_0_0_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sssk + 30);

    auto g_0_0_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sssk + 31);

    auto g_0_0_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sssk + 32);

    auto g_0_0_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sssk + 33);

    auto g_0_0_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sssk + 34);

    auto g_0_0_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sssk + 35);

    /// Set up components of auxilary buffer : SSSK

    auto g_0_0_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sssk);

    auto g_0_0_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sssk + 1);

    auto g_0_0_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sssk + 2);

    auto g_0_0_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sssk + 3);

    auto g_0_0_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sssk + 4);

    auto g_0_0_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sssk + 5);

    auto g_0_0_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sssk + 6);

    auto g_0_0_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sssk + 7);

    auto g_0_0_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sssk + 8);

    auto g_0_0_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sssk + 9);

    auto g_0_0_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sssk + 10);

    auto g_0_0_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sssk + 11);

    auto g_0_0_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sssk + 12);

    auto g_0_0_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sssk + 13);

    auto g_0_0_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sssk + 14);

    auto g_0_0_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sssk + 15);

    auto g_0_0_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sssk + 16);

    auto g_0_0_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sssk + 17);

    auto g_0_0_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sssk + 18);

    auto g_0_0_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sssk + 19);

    auto g_0_0_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sssk + 20);

    auto g_0_0_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sssk + 21);

    auto g_0_0_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sssk + 22);

    auto g_0_0_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sssk + 23);

    auto g_0_0_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sssk + 24);

    auto g_0_0_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sssk + 25);

    auto g_0_0_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sssk + 26);

    auto g_0_0_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sssk + 27);

    auto g_0_0_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sssk + 28);

    auto g_0_0_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sssk + 29);

    auto g_0_0_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sssk + 30);

    auto g_0_0_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sssk + 31);

    auto g_0_0_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sssk + 32);

    auto g_0_0_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sssk + 33);

    auto g_0_0_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sssk + 34);

    auto g_0_0_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sssk + 35);

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

    /// Set up components of auxilary buffer : SPSK

    auto g_0_x_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_spsk);

    auto g_0_x_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_spsk + 1);

    auto g_0_x_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_spsk + 2);

    auto g_0_x_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_spsk + 3);

    auto g_0_x_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_spsk + 4);

    auto g_0_x_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_spsk + 5);

    auto g_0_x_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_spsk + 6);

    auto g_0_x_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_spsk + 7);

    auto g_0_x_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_spsk + 8);

    auto g_0_x_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_spsk + 9);

    auto g_0_x_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_spsk + 10);

    auto g_0_x_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_spsk + 11);

    auto g_0_x_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_spsk + 12);

    auto g_0_x_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_spsk + 13);

    auto g_0_x_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_spsk + 14);

    auto g_0_x_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_spsk + 15);

    auto g_0_x_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_spsk + 16);

    auto g_0_x_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_spsk + 17);

    auto g_0_x_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_spsk + 18);

    auto g_0_x_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_spsk + 19);

    auto g_0_x_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_spsk + 20);

    auto g_0_x_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_spsk + 21);

    auto g_0_x_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_spsk + 22);

    auto g_0_x_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_spsk + 23);

    auto g_0_x_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_spsk + 24);

    auto g_0_x_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_spsk + 25);

    auto g_0_x_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_spsk + 26);

    auto g_0_x_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_spsk + 27);

    auto g_0_x_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_spsk + 28);

    auto g_0_x_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_spsk + 29);

    auto g_0_x_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_spsk + 30);

    auto g_0_x_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_spsk + 31);

    auto g_0_x_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_spsk + 32);

    auto g_0_x_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_spsk + 33);

    auto g_0_x_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_spsk + 34);

    auto g_0_x_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_spsk + 35);

    auto g_0_y_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_spsk + 36);

    auto g_0_y_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_spsk + 37);

    auto g_0_y_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_spsk + 38);

    auto g_0_y_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_spsk + 39);

    auto g_0_y_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_spsk + 40);

    auto g_0_y_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_spsk + 41);

    auto g_0_y_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_spsk + 42);

    auto g_0_y_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_spsk + 43);

    auto g_0_y_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_spsk + 44);

    auto g_0_y_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_spsk + 45);

    auto g_0_y_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_spsk + 46);

    auto g_0_y_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_spsk + 47);

    auto g_0_y_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_spsk + 48);

    auto g_0_y_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_spsk + 49);

    auto g_0_y_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_spsk + 50);

    auto g_0_y_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_spsk + 51);

    auto g_0_y_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_spsk + 52);

    auto g_0_y_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_spsk + 53);

    auto g_0_y_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_spsk + 54);

    auto g_0_y_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_spsk + 55);

    auto g_0_y_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_spsk + 56);

    auto g_0_y_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_spsk + 57);

    auto g_0_y_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_spsk + 58);

    auto g_0_y_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_spsk + 59);

    auto g_0_y_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_spsk + 60);

    auto g_0_y_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_spsk + 61);

    auto g_0_y_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_spsk + 62);

    auto g_0_y_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_spsk + 63);

    auto g_0_y_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_spsk + 64);

    auto g_0_y_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_spsk + 65);

    auto g_0_y_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_spsk + 66);

    auto g_0_y_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_spsk + 67);

    auto g_0_y_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_spsk + 68);

    auto g_0_y_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_spsk + 69);

    auto g_0_y_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_spsk + 70);

    auto g_0_y_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_spsk + 71);

    auto g_0_z_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_spsk + 72);

    auto g_0_z_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_spsk + 73);

    auto g_0_z_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_spsk + 74);

    auto g_0_z_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_spsk + 75);

    auto g_0_z_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_spsk + 76);

    auto g_0_z_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_spsk + 77);

    auto g_0_z_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_spsk + 78);

    auto g_0_z_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_spsk + 79);

    auto g_0_z_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_spsk + 80);

    auto g_0_z_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_spsk + 81);

    auto g_0_z_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_spsk + 82);

    auto g_0_z_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_spsk + 83);

    auto g_0_z_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_spsk + 84);

    auto g_0_z_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_spsk + 85);

    auto g_0_z_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_spsk + 86);

    auto g_0_z_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_spsk + 87);

    auto g_0_z_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_spsk + 88);

    auto g_0_z_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_spsk + 89);

    auto g_0_z_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_spsk + 90);

    auto g_0_z_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_spsk + 91);

    auto g_0_z_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_spsk + 92);

    auto g_0_z_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_spsk + 93);

    auto g_0_z_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_spsk + 94);

    auto g_0_z_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_spsk + 95);

    auto g_0_z_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_spsk + 96);

    auto g_0_z_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_spsk + 97);

    auto g_0_z_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_spsk + 98);

    auto g_0_z_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_spsk + 99);

    auto g_0_z_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_spsk + 100);

    auto g_0_z_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_spsk + 101);

    auto g_0_z_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_spsk + 102);

    auto g_0_z_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_spsk + 103);

    auto g_0_z_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_spsk + 104);

    auto g_0_z_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_spsk + 105);

    auto g_0_z_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_spsk + 106);

    auto g_0_z_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_spsk + 107);

    /// Set up components of auxilary buffer : SPSK

    auto g_0_x_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_spsk);

    auto g_0_x_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_spsk + 1);

    auto g_0_x_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_spsk + 2);

    auto g_0_x_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_spsk + 3);

    auto g_0_x_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_spsk + 4);

    auto g_0_x_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_spsk + 5);

    auto g_0_x_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_spsk + 6);

    auto g_0_x_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_spsk + 7);

    auto g_0_x_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_spsk + 8);

    auto g_0_x_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_spsk + 9);

    auto g_0_x_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_spsk + 10);

    auto g_0_x_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_spsk + 11);

    auto g_0_x_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_spsk + 12);

    auto g_0_x_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_spsk + 13);

    auto g_0_x_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_spsk + 14);

    auto g_0_x_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_spsk + 15);

    auto g_0_x_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_spsk + 16);

    auto g_0_x_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_spsk + 17);

    auto g_0_x_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_spsk + 18);

    auto g_0_x_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_spsk + 19);

    auto g_0_x_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_spsk + 20);

    auto g_0_x_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_spsk + 21);

    auto g_0_x_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_spsk + 22);

    auto g_0_x_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_spsk + 23);

    auto g_0_x_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_spsk + 24);

    auto g_0_x_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_spsk + 25);

    auto g_0_x_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_spsk + 26);

    auto g_0_x_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_spsk + 27);

    auto g_0_x_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_spsk + 28);

    auto g_0_x_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_spsk + 29);

    auto g_0_x_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_spsk + 30);

    auto g_0_x_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_spsk + 31);

    auto g_0_x_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_spsk + 32);

    auto g_0_x_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_spsk + 33);

    auto g_0_x_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_spsk + 34);

    auto g_0_x_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_spsk + 35);

    auto g_0_y_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_spsk + 36);

    auto g_0_y_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_spsk + 37);

    auto g_0_y_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_spsk + 38);

    auto g_0_y_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_spsk + 39);

    auto g_0_y_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_spsk + 40);

    auto g_0_y_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_spsk + 41);

    auto g_0_y_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_spsk + 42);

    auto g_0_y_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_spsk + 43);

    auto g_0_y_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_spsk + 44);

    auto g_0_y_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_spsk + 45);

    auto g_0_y_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_spsk + 46);

    auto g_0_y_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_spsk + 47);

    auto g_0_y_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_spsk + 48);

    auto g_0_y_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_spsk + 49);

    auto g_0_y_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_spsk + 50);

    auto g_0_y_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_spsk + 51);

    auto g_0_y_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_spsk + 52);

    auto g_0_y_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_spsk + 53);

    auto g_0_y_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_spsk + 54);

    auto g_0_y_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_spsk + 55);

    auto g_0_y_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_spsk + 56);

    auto g_0_y_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_spsk + 57);

    auto g_0_y_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_spsk + 58);

    auto g_0_y_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_spsk + 59);

    auto g_0_y_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_spsk + 60);

    auto g_0_y_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_spsk + 61);

    auto g_0_y_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_spsk + 62);

    auto g_0_y_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_spsk + 63);

    auto g_0_y_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_spsk + 64);

    auto g_0_y_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_spsk + 65);

    auto g_0_y_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_spsk + 66);

    auto g_0_y_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_spsk + 67);

    auto g_0_y_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_spsk + 68);

    auto g_0_y_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_spsk + 69);

    auto g_0_y_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_spsk + 70);

    auto g_0_y_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_spsk + 71);

    auto g_0_z_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_spsk + 72);

    auto g_0_z_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_spsk + 73);

    auto g_0_z_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_spsk + 74);

    auto g_0_z_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_spsk + 75);

    auto g_0_z_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_spsk + 76);

    auto g_0_z_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_spsk + 77);

    auto g_0_z_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_spsk + 78);

    auto g_0_z_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_spsk + 79);

    auto g_0_z_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_spsk + 80);

    auto g_0_z_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_spsk + 81);

    auto g_0_z_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_spsk + 82);

    auto g_0_z_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_spsk + 83);

    auto g_0_z_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_spsk + 84);

    auto g_0_z_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_spsk + 85);

    auto g_0_z_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_spsk + 86);

    auto g_0_z_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_spsk + 87);

    auto g_0_z_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_spsk + 88);

    auto g_0_z_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_spsk + 89);

    auto g_0_z_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_spsk + 90);

    auto g_0_z_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_spsk + 91);

    auto g_0_z_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_spsk + 92);

    auto g_0_z_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_spsk + 93);

    auto g_0_z_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_spsk + 94);

    auto g_0_z_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_spsk + 95);

    auto g_0_z_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_spsk + 96);

    auto g_0_z_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_spsk + 97);

    auto g_0_z_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_spsk + 98);

    auto g_0_z_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_spsk + 99);

    auto g_0_z_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_spsk + 100);

    auto g_0_z_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_spsk + 101);

    auto g_0_z_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_spsk + 102);

    auto g_0_z_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_spsk + 103);

    auto g_0_z_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_spsk + 104);

    auto g_0_z_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_spsk + 105);

    auto g_0_z_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_spsk + 106);

    auto g_0_z_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_spsk + 107);

    /// Set up 0-36 components of targeted buffer : SDSK

    auto g_0_xx_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sdsk);

    auto g_0_xx_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sdsk + 1);

    auto g_0_xx_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sdsk + 2);

    auto g_0_xx_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sdsk + 3);

    auto g_0_xx_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sdsk + 4);

    auto g_0_xx_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sdsk + 5);

    auto g_0_xx_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sdsk + 6);

    auto g_0_xx_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sdsk + 7);

    auto g_0_xx_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sdsk + 8);

    auto g_0_xx_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sdsk + 9);

    auto g_0_xx_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sdsk + 10);

    auto g_0_xx_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sdsk + 11);

    auto g_0_xx_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sdsk + 12);

    auto g_0_xx_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sdsk + 13);

    auto g_0_xx_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 14);

    auto g_0_xx_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sdsk + 15);

    auto g_0_xx_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sdsk + 16);

    auto g_0_xx_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sdsk + 17);

    auto g_0_xx_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sdsk + 18);

    auto g_0_xx_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 19);

    auto g_0_xx_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 20);

    auto g_0_xx_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sdsk + 21);

    auto g_0_xx_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sdsk + 22);

    auto g_0_xx_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sdsk + 23);

    auto g_0_xx_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sdsk + 24);

    auto g_0_xx_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 25);

    auto g_0_xx_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 26);

    auto g_0_xx_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 27);

    auto g_0_xx_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sdsk + 28);

    auto g_0_xx_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sdsk + 29);

    auto g_0_xx_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sdsk + 30);

    auto g_0_xx_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sdsk + 31);

    auto g_0_xx_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 32);

    auto g_0_xx_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 33);

    auto g_0_xx_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 34);

    auto g_0_xx_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 35);

    #pragma omp simd aligned(g_0_0_0_xxxxxxx_0, g_0_0_0_xxxxxxx_1, g_0_0_0_xxxxxxy_0, g_0_0_0_xxxxxxy_1, g_0_0_0_xxxxxxz_0, g_0_0_0_xxxxxxz_1, g_0_0_0_xxxxxyy_0, g_0_0_0_xxxxxyy_1, g_0_0_0_xxxxxyz_0, g_0_0_0_xxxxxyz_1, g_0_0_0_xxxxxzz_0, g_0_0_0_xxxxxzz_1, g_0_0_0_xxxxyyy_0, g_0_0_0_xxxxyyy_1, g_0_0_0_xxxxyyz_0, g_0_0_0_xxxxyyz_1, g_0_0_0_xxxxyzz_0, g_0_0_0_xxxxyzz_1, g_0_0_0_xxxxzzz_0, g_0_0_0_xxxxzzz_1, g_0_0_0_xxxyyyy_0, g_0_0_0_xxxyyyy_1, g_0_0_0_xxxyyyz_0, g_0_0_0_xxxyyyz_1, g_0_0_0_xxxyyzz_0, g_0_0_0_xxxyyzz_1, g_0_0_0_xxxyzzz_0, g_0_0_0_xxxyzzz_1, g_0_0_0_xxxzzzz_0, g_0_0_0_xxxzzzz_1, g_0_0_0_xxyyyyy_0, g_0_0_0_xxyyyyy_1, g_0_0_0_xxyyyyz_0, g_0_0_0_xxyyyyz_1, g_0_0_0_xxyyyzz_0, g_0_0_0_xxyyyzz_1, g_0_0_0_xxyyzzz_0, g_0_0_0_xxyyzzz_1, g_0_0_0_xxyzzzz_0, g_0_0_0_xxyzzzz_1, g_0_0_0_xxzzzzz_0, g_0_0_0_xxzzzzz_1, g_0_0_0_xyyyyyy_0, g_0_0_0_xyyyyyy_1, g_0_0_0_xyyyyyz_0, g_0_0_0_xyyyyyz_1, g_0_0_0_xyyyyzz_0, g_0_0_0_xyyyyzz_1, g_0_0_0_xyyyzzz_0, g_0_0_0_xyyyzzz_1, g_0_0_0_xyyzzzz_0, g_0_0_0_xyyzzzz_1, g_0_0_0_xyzzzzz_0, g_0_0_0_xyzzzzz_1, g_0_0_0_xzzzzzz_0, g_0_0_0_xzzzzzz_1, g_0_0_0_yyyyyyy_0, g_0_0_0_yyyyyyy_1, g_0_0_0_yyyyyyz_0, g_0_0_0_yyyyyyz_1, g_0_0_0_yyyyyzz_0, g_0_0_0_yyyyyzz_1, g_0_0_0_yyyyzzz_0, g_0_0_0_yyyyzzz_1, g_0_0_0_yyyzzzz_0, g_0_0_0_yyyzzzz_1, g_0_0_0_yyzzzzz_0, g_0_0_0_yyzzzzz_1, g_0_0_0_yzzzzzz_0, g_0_0_0_yzzzzzz_1, g_0_0_0_zzzzzzz_0, g_0_0_0_zzzzzzz_1, g_0_x_0_xxxxxx_1, g_0_x_0_xxxxxxx_0, g_0_x_0_xxxxxxx_1, g_0_x_0_xxxxxxy_0, g_0_x_0_xxxxxxy_1, g_0_x_0_xxxxxxz_0, g_0_x_0_xxxxxxz_1, g_0_x_0_xxxxxy_1, g_0_x_0_xxxxxyy_0, g_0_x_0_xxxxxyy_1, g_0_x_0_xxxxxyz_0, g_0_x_0_xxxxxyz_1, g_0_x_0_xxxxxz_1, g_0_x_0_xxxxxzz_0, g_0_x_0_xxxxxzz_1, g_0_x_0_xxxxyy_1, g_0_x_0_xxxxyyy_0, g_0_x_0_xxxxyyy_1, g_0_x_0_xxxxyyz_0, g_0_x_0_xxxxyyz_1, g_0_x_0_xxxxyz_1, g_0_x_0_xxxxyzz_0, g_0_x_0_xxxxyzz_1, g_0_x_0_xxxxzz_1, g_0_x_0_xxxxzzz_0, g_0_x_0_xxxxzzz_1, g_0_x_0_xxxyyy_1, g_0_x_0_xxxyyyy_0, g_0_x_0_xxxyyyy_1, g_0_x_0_xxxyyyz_0, g_0_x_0_xxxyyyz_1, g_0_x_0_xxxyyz_1, g_0_x_0_xxxyyzz_0, g_0_x_0_xxxyyzz_1, g_0_x_0_xxxyzz_1, g_0_x_0_xxxyzzz_0, g_0_x_0_xxxyzzz_1, g_0_x_0_xxxzzz_1, g_0_x_0_xxxzzzz_0, g_0_x_0_xxxzzzz_1, g_0_x_0_xxyyyy_1, g_0_x_0_xxyyyyy_0, g_0_x_0_xxyyyyy_1, g_0_x_0_xxyyyyz_0, g_0_x_0_xxyyyyz_1, g_0_x_0_xxyyyz_1, g_0_x_0_xxyyyzz_0, g_0_x_0_xxyyyzz_1, g_0_x_0_xxyyzz_1, g_0_x_0_xxyyzzz_0, g_0_x_0_xxyyzzz_1, g_0_x_0_xxyzzz_1, g_0_x_0_xxyzzzz_0, g_0_x_0_xxyzzzz_1, g_0_x_0_xxzzzz_1, g_0_x_0_xxzzzzz_0, g_0_x_0_xxzzzzz_1, g_0_x_0_xyyyyy_1, g_0_x_0_xyyyyyy_0, g_0_x_0_xyyyyyy_1, g_0_x_0_xyyyyyz_0, g_0_x_0_xyyyyyz_1, g_0_x_0_xyyyyz_1, g_0_x_0_xyyyyzz_0, g_0_x_0_xyyyyzz_1, g_0_x_0_xyyyzz_1, g_0_x_0_xyyyzzz_0, g_0_x_0_xyyyzzz_1, g_0_x_0_xyyzzz_1, g_0_x_0_xyyzzzz_0, g_0_x_0_xyyzzzz_1, g_0_x_0_xyzzzz_1, g_0_x_0_xyzzzzz_0, g_0_x_0_xyzzzzz_1, g_0_x_0_xzzzzz_1, g_0_x_0_xzzzzzz_0, g_0_x_0_xzzzzzz_1, g_0_x_0_yyyyyy_1, g_0_x_0_yyyyyyy_0, g_0_x_0_yyyyyyy_1, g_0_x_0_yyyyyyz_0, g_0_x_0_yyyyyyz_1, g_0_x_0_yyyyyz_1, g_0_x_0_yyyyyzz_0, g_0_x_0_yyyyyzz_1, g_0_x_0_yyyyzz_1, g_0_x_0_yyyyzzz_0, g_0_x_0_yyyyzzz_1, g_0_x_0_yyyzzz_1, g_0_x_0_yyyzzzz_0, g_0_x_0_yyyzzzz_1, g_0_x_0_yyzzzz_1, g_0_x_0_yyzzzzz_0, g_0_x_0_yyzzzzz_1, g_0_x_0_yzzzzz_1, g_0_x_0_yzzzzzz_0, g_0_x_0_yzzzzzz_1, g_0_x_0_zzzzzz_1, g_0_x_0_zzzzzzz_0, g_0_x_0_zzzzzzz_1, g_0_xx_0_xxxxxxx_0, g_0_xx_0_xxxxxxy_0, g_0_xx_0_xxxxxxz_0, g_0_xx_0_xxxxxyy_0, g_0_xx_0_xxxxxyz_0, g_0_xx_0_xxxxxzz_0, g_0_xx_0_xxxxyyy_0, g_0_xx_0_xxxxyyz_0, g_0_xx_0_xxxxyzz_0, g_0_xx_0_xxxxzzz_0, g_0_xx_0_xxxyyyy_0, g_0_xx_0_xxxyyyz_0, g_0_xx_0_xxxyyzz_0, g_0_xx_0_xxxyzzz_0, g_0_xx_0_xxxzzzz_0, g_0_xx_0_xxyyyyy_0, g_0_xx_0_xxyyyyz_0, g_0_xx_0_xxyyyzz_0, g_0_xx_0_xxyyzzz_0, g_0_xx_0_xxyzzzz_0, g_0_xx_0_xxzzzzz_0, g_0_xx_0_xyyyyyy_0, g_0_xx_0_xyyyyyz_0, g_0_xx_0_xyyyyzz_0, g_0_xx_0_xyyyzzz_0, g_0_xx_0_xyyzzzz_0, g_0_xx_0_xyzzzzz_0, g_0_xx_0_xzzzzzz_0, g_0_xx_0_yyyyyyy_0, g_0_xx_0_yyyyyyz_0, g_0_xx_0_yyyyyzz_0, g_0_xx_0_yyyyzzz_0, g_0_xx_0_yyyzzzz_0, g_0_xx_0_yyzzzzz_0, g_0_xx_0_yzzzzzz_0, g_0_xx_0_zzzzzzz_0, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xx_0_xxxxxxx_0[i] = g_0_0_0_xxxxxxx_0[i] * fi_ab_0 - g_0_0_0_xxxxxxx_1[i] * fti_ab_0 + 7.0 * g_0_x_0_xxxxxx_1[i] * fi_abcd_0 + g_0_x_0_xxxxxxx_0[i] * pb_x + g_0_x_0_xxxxxxx_1[i] * wp_x[i];

        g_0_xx_0_xxxxxxy_0[i] = g_0_0_0_xxxxxxy_0[i] * fi_ab_0 - g_0_0_0_xxxxxxy_1[i] * fti_ab_0 + 6.0 * g_0_x_0_xxxxxy_1[i] * fi_abcd_0 + g_0_x_0_xxxxxxy_0[i] * pb_x + g_0_x_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xx_0_xxxxxxz_0[i] = g_0_0_0_xxxxxxz_0[i] * fi_ab_0 - g_0_0_0_xxxxxxz_1[i] * fti_ab_0 + 6.0 * g_0_x_0_xxxxxz_1[i] * fi_abcd_0 + g_0_x_0_xxxxxxz_0[i] * pb_x + g_0_x_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xx_0_xxxxxyy_0[i] = g_0_0_0_xxxxxyy_0[i] * fi_ab_0 - g_0_0_0_xxxxxyy_1[i] * fti_ab_0 + 5.0 * g_0_x_0_xxxxyy_1[i] * fi_abcd_0 + g_0_x_0_xxxxxyy_0[i] * pb_x + g_0_x_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xx_0_xxxxxyz_0[i] = g_0_0_0_xxxxxyz_0[i] * fi_ab_0 - g_0_0_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_x_0_xxxxyz_1[i] * fi_abcd_0 + g_0_x_0_xxxxxyz_0[i] * pb_x + g_0_x_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xx_0_xxxxxzz_0[i] = g_0_0_0_xxxxxzz_0[i] * fi_ab_0 - g_0_0_0_xxxxxzz_1[i] * fti_ab_0 + 5.0 * g_0_x_0_xxxxzz_1[i] * fi_abcd_0 + g_0_x_0_xxxxxzz_0[i] * pb_x + g_0_x_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xx_0_xxxxyyy_0[i] = g_0_0_0_xxxxyyy_0[i] * fi_ab_0 - g_0_0_0_xxxxyyy_1[i] * fti_ab_0 + 4.0 * g_0_x_0_xxxyyy_1[i] * fi_abcd_0 + g_0_x_0_xxxxyyy_0[i] * pb_x + g_0_x_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xx_0_xxxxyyz_0[i] = g_0_0_0_xxxxyyz_0[i] * fi_ab_0 - g_0_0_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_x_0_xxxyyz_1[i] * fi_abcd_0 + g_0_x_0_xxxxyyz_0[i] * pb_x + g_0_x_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xx_0_xxxxyzz_0[i] = g_0_0_0_xxxxyzz_0[i] * fi_ab_0 - g_0_0_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_x_0_xxxyzz_1[i] * fi_abcd_0 + g_0_x_0_xxxxyzz_0[i] * pb_x + g_0_x_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xx_0_xxxxzzz_0[i] = g_0_0_0_xxxxzzz_0[i] * fi_ab_0 - g_0_0_0_xxxxzzz_1[i] * fti_ab_0 + 4.0 * g_0_x_0_xxxzzz_1[i] * fi_abcd_0 + g_0_x_0_xxxxzzz_0[i] * pb_x + g_0_x_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xx_0_xxxyyyy_0[i] = g_0_0_0_xxxyyyy_0[i] * fi_ab_0 - g_0_0_0_xxxyyyy_1[i] * fti_ab_0 + 3.0 * g_0_x_0_xxyyyy_1[i] * fi_abcd_0 + g_0_x_0_xxxyyyy_0[i] * pb_x + g_0_x_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xx_0_xxxyyyz_0[i] = g_0_0_0_xxxyyyz_0[i] * fi_ab_0 - g_0_0_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_x_0_xxyyyz_1[i] * fi_abcd_0 + g_0_x_0_xxxyyyz_0[i] * pb_x + g_0_x_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xx_0_xxxyyzz_0[i] = g_0_0_0_xxxyyzz_0[i] * fi_ab_0 - g_0_0_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_x_0_xxyyzz_1[i] * fi_abcd_0 + g_0_x_0_xxxyyzz_0[i] * pb_x + g_0_x_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xx_0_xxxyzzz_0[i] = g_0_0_0_xxxyzzz_0[i] * fi_ab_0 - g_0_0_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_x_0_xxyzzz_1[i] * fi_abcd_0 + g_0_x_0_xxxyzzz_0[i] * pb_x + g_0_x_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xx_0_xxxzzzz_0[i] = g_0_0_0_xxxzzzz_0[i] * fi_ab_0 - g_0_0_0_xxxzzzz_1[i] * fti_ab_0 + 3.0 * g_0_x_0_xxzzzz_1[i] * fi_abcd_0 + g_0_x_0_xxxzzzz_0[i] * pb_x + g_0_x_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xx_0_xxyyyyy_0[i] = g_0_0_0_xxyyyyy_0[i] * fi_ab_0 - g_0_0_0_xxyyyyy_1[i] * fti_ab_0 + 2.0 * g_0_x_0_xyyyyy_1[i] * fi_abcd_0 + g_0_x_0_xxyyyyy_0[i] * pb_x + g_0_x_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xx_0_xxyyyyz_0[i] = g_0_0_0_xxyyyyz_0[i] * fi_ab_0 - g_0_0_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_x_0_xyyyyz_1[i] * fi_abcd_0 + g_0_x_0_xxyyyyz_0[i] * pb_x + g_0_x_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xx_0_xxyyyzz_0[i] = g_0_0_0_xxyyyzz_0[i] * fi_ab_0 - g_0_0_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_x_0_xyyyzz_1[i] * fi_abcd_0 + g_0_x_0_xxyyyzz_0[i] * pb_x + g_0_x_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xx_0_xxyyzzz_0[i] = g_0_0_0_xxyyzzz_0[i] * fi_ab_0 - g_0_0_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_x_0_xyyzzz_1[i] * fi_abcd_0 + g_0_x_0_xxyyzzz_0[i] * pb_x + g_0_x_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xx_0_xxyzzzz_0[i] = g_0_0_0_xxyzzzz_0[i] * fi_ab_0 - g_0_0_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_x_0_xyzzzz_1[i] * fi_abcd_0 + g_0_x_0_xxyzzzz_0[i] * pb_x + g_0_x_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xx_0_xxzzzzz_0[i] = g_0_0_0_xxzzzzz_0[i] * fi_ab_0 - g_0_0_0_xxzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_x_0_xzzzzz_1[i] * fi_abcd_0 + g_0_x_0_xxzzzzz_0[i] * pb_x + g_0_x_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xx_0_xyyyyyy_0[i] = g_0_0_0_xyyyyyy_0[i] * fi_ab_0 - g_0_0_0_xyyyyyy_1[i] * fti_ab_0 + g_0_x_0_yyyyyy_1[i] * fi_abcd_0 + g_0_x_0_xyyyyyy_0[i] * pb_x + g_0_x_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xx_0_xyyyyyz_0[i] = g_0_0_0_xyyyyyz_0[i] * fi_ab_0 - g_0_0_0_xyyyyyz_1[i] * fti_ab_0 + g_0_x_0_yyyyyz_1[i] * fi_abcd_0 + g_0_x_0_xyyyyyz_0[i] * pb_x + g_0_x_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xx_0_xyyyyzz_0[i] = g_0_0_0_xyyyyzz_0[i] * fi_ab_0 - g_0_0_0_xyyyyzz_1[i] * fti_ab_0 + g_0_x_0_yyyyzz_1[i] * fi_abcd_0 + g_0_x_0_xyyyyzz_0[i] * pb_x + g_0_x_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xx_0_xyyyzzz_0[i] = g_0_0_0_xyyyzzz_0[i] * fi_ab_0 - g_0_0_0_xyyyzzz_1[i] * fti_ab_0 + g_0_x_0_yyyzzz_1[i] * fi_abcd_0 + g_0_x_0_xyyyzzz_0[i] * pb_x + g_0_x_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xx_0_xyyzzzz_0[i] = g_0_0_0_xyyzzzz_0[i] * fi_ab_0 - g_0_0_0_xyyzzzz_1[i] * fti_ab_0 + g_0_x_0_yyzzzz_1[i] * fi_abcd_0 + g_0_x_0_xyyzzzz_0[i] * pb_x + g_0_x_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xx_0_xyzzzzz_0[i] = g_0_0_0_xyzzzzz_0[i] * fi_ab_0 - g_0_0_0_xyzzzzz_1[i] * fti_ab_0 + g_0_x_0_yzzzzz_1[i] * fi_abcd_0 + g_0_x_0_xyzzzzz_0[i] * pb_x + g_0_x_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xx_0_xzzzzzz_0[i] = g_0_0_0_xzzzzzz_0[i] * fi_ab_0 - g_0_0_0_xzzzzzz_1[i] * fti_ab_0 + g_0_x_0_zzzzzz_1[i] * fi_abcd_0 + g_0_x_0_xzzzzzz_0[i] * pb_x + g_0_x_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xx_0_yyyyyyy_0[i] = g_0_0_0_yyyyyyy_0[i] * fi_ab_0 - g_0_0_0_yyyyyyy_1[i] * fti_ab_0 + g_0_x_0_yyyyyyy_0[i] * pb_x + g_0_x_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xx_0_yyyyyyz_0[i] = g_0_0_0_yyyyyyz_0[i] * fi_ab_0 - g_0_0_0_yyyyyyz_1[i] * fti_ab_0 + g_0_x_0_yyyyyyz_0[i] * pb_x + g_0_x_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xx_0_yyyyyzz_0[i] = g_0_0_0_yyyyyzz_0[i] * fi_ab_0 - g_0_0_0_yyyyyzz_1[i] * fti_ab_0 + g_0_x_0_yyyyyzz_0[i] * pb_x + g_0_x_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xx_0_yyyyzzz_0[i] = g_0_0_0_yyyyzzz_0[i] * fi_ab_0 - g_0_0_0_yyyyzzz_1[i] * fti_ab_0 + g_0_x_0_yyyyzzz_0[i] * pb_x + g_0_x_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xx_0_yyyzzzz_0[i] = g_0_0_0_yyyzzzz_0[i] * fi_ab_0 - g_0_0_0_yyyzzzz_1[i] * fti_ab_0 + g_0_x_0_yyyzzzz_0[i] * pb_x + g_0_x_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xx_0_yyzzzzz_0[i] = g_0_0_0_yyzzzzz_0[i] * fi_ab_0 - g_0_0_0_yyzzzzz_1[i] * fti_ab_0 + g_0_x_0_yyzzzzz_0[i] * pb_x + g_0_x_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xx_0_yzzzzzz_0[i] = g_0_0_0_yzzzzzz_0[i] * fi_ab_0 - g_0_0_0_yzzzzzz_1[i] * fti_ab_0 + g_0_x_0_yzzzzzz_0[i] * pb_x + g_0_x_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xx_0_zzzzzzz_0[i] = g_0_0_0_zzzzzzz_0[i] * fi_ab_0 - g_0_0_0_zzzzzzz_1[i] * fti_ab_0 + g_0_x_0_zzzzzzz_0[i] * pb_x + g_0_x_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 36-72 components of targeted buffer : SDSK

    auto g_0_xy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sdsk + 36);

    auto g_0_xy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sdsk + 37);

    auto g_0_xy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sdsk + 38);

    auto g_0_xy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sdsk + 39);

    auto g_0_xy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sdsk + 40);

    auto g_0_xy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sdsk + 41);

    auto g_0_xy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sdsk + 42);

    auto g_0_xy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sdsk + 43);

    auto g_0_xy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sdsk + 44);

    auto g_0_xy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sdsk + 45);

    auto g_0_xy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sdsk + 46);

    auto g_0_xy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sdsk + 47);

    auto g_0_xy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sdsk + 48);

    auto g_0_xy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sdsk + 49);

    auto g_0_xy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 50);

    auto g_0_xy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sdsk + 51);

    auto g_0_xy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sdsk + 52);

    auto g_0_xy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sdsk + 53);

    auto g_0_xy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sdsk + 54);

    auto g_0_xy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 55);

    auto g_0_xy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 56);

    auto g_0_xy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sdsk + 57);

    auto g_0_xy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sdsk + 58);

    auto g_0_xy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sdsk + 59);

    auto g_0_xy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sdsk + 60);

    auto g_0_xy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 61);

    auto g_0_xy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 62);

    auto g_0_xy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 63);

    auto g_0_xy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sdsk + 64);

    auto g_0_xy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sdsk + 65);

    auto g_0_xy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sdsk + 66);

    auto g_0_xy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sdsk + 67);

    auto g_0_xy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 68);

    auto g_0_xy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 69);

    auto g_0_xy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 70);

    auto g_0_xy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 71);

    #pragma omp simd aligned(g_0_x_0_xxxxxxx_0, g_0_x_0_xxxxxxx_1, g_0_x_0_xxxxxxz_0, g_0_x_0_xxxxxxz_1, g_0_x_0_xxxxxzz_0, g_0_x_0_xxxxxzz_1, g_0_x_0_xxxxzzz_0, g_0_x_0_xxxxzzz_1, g_0_x_0_xxxzzzz_0, g_0_x_0_xxxzzzz_1, g_0_x_0_xxzzzzz_0, g_0_x_0_xxzzzzz_1, g_0_x_0_xzzzzzz_0, g_0_x_0_xzzzzzz_1, g_0_xy_0_xxxxxxx_0, g_0_xy_0_xxxxxxy_0, g_0_xy_0_xxxxxxz_0, g_0_xy_0_xxxxxyy_0, g_0_xy_0_xxxxxyz_0, g_0_xy_0_xxxxxzz_0, g_0_xy_0_xxxxyyy_0, g_0_xy_0_xxxxyyz_0, g_0_xy_0_xxxxyzz_0, g_0_xy_0_xxxxzzz_0, g_0_xy_0_xxxyyyy_0, g_0_xy_0_xxxyyyz_0, g_0_xy_0_xxxyyzz_0, g_0_xy_0_xxxyzzz_0, g_0_xy_0_xxxzzzz_0, g_0_xy_0_xxyyyyy_0, g_0_xy_0_xxyyyyz_0, g_0_xy_0_xxyyyzz_0, g_0_xy_0_xxyyzzz_0, g_0_xy_0_xxyzzzz_0, g_0_xy_0_xxzzzzz_0, g_0_xy_0_xyyyyyy_0, g_0_xy_0_xyyyyyz_0, g_0_xy_0_xyyyyzz_0, g_0_xy_0_xyyyzzz_0, g_0_xy_0_xyyzzzz_0, g_0_xy_0_xyzzzzz_0, g_0_xy_0_xzzzzzz_0, g_0_xy_0_yyyyyyy_0, g_0_xy_0_yyyyyyz_0, g_0_xy_0_yyyyyzz_0, g_0_xy_0_yyyyzzz_0, g_0_xy_0_yyyzzzz_0, g_0_xy_0_yyzzzzz_0, g_0_xy_0_yzzzzzz_0, g_0_xy_0_zzzzzzz_0, g_0_y_0_xxxxxxy_0, g_0_y_0_xxxxxxy_1, g_0_y_0_xxxxxy_1, g_0_y_0_xxxxxyy_0, g_0_y_0_xxxxxyy_1, g_0_y_0_xxxxxyz_0, g_0_y_0_xxxxxyz_1, g_0_y_0_xxxxyy_1, g_0_y_0_xxxxyyy_0, g_0_y_0_xxxxyyy_1, g_0_y_0_xxxxyyz_0, g_0_y_0_xxxxyyz_1, g_0_y_0_xxxxyz_1, g_0_y_0_xxxxyzz_0, g_0_y_0_xxxxyzz_1, g_0_y_0_xxxyyy_1, g_0_y_0_xxxyyyy_0, g_0_y_0_xxxyyyy_1, g_0_y_0_xxxyyyz_0, g_0_y_0_xxxyyyz_1, g_0_y_0_xxxyyz_1, g_0_y_0_xxxyyzz_0, g_0_y_0_xxxyyzz_1, g_0_y_0_xxxyzz_1, g_0_y_0_xxxyzzz_0, g_0_y_0_xxxyzzz_1, g_0_y_0_xxyyyy_1, g_0_y_0_xxyyyyy_0, g_0_y_0_xxyyyyy_1, g_0_y_0_xxyyyyz_0, g_0_y_0_xxyyyyz_1, g_0_y_0_xxyyyz_1, g_0_y_0_xxyyyzz_0, g_0_y_0_xxyyyzz_1, g_0_y_0_xxyyzz_1, g_0_y_0_xxyyzzz_0, g_0_y_0_xxyyzzz_1, g_0_y_0_xxyzzz_1, g_0_y_0_xxyzzzz_0, g_0_y_0_xxyzzzz_1, g_0_y_0_xyyyyy_1, g_0_y_0_xyyyyyy_0, g_0_y_0_xyyyyyy_1, g_0_y_0_xyyyyyz_0, g_0_y_0_xyyyyyz_1, g_0_y_0_xyyyyz_1, g_0_y_0_xyyyyzz_0, g_0_y_0_xyyyyzz_1, g_0_y_0_xyyyzz_1, g_0_y_0_xyyyzzz_0, g_0_y_0_xyyyzzz_1, g_0_y_0_xyyzzz_1, g_0_y_0_xyyzzzz_0, g_0_y_0_xyyzzzz_1, g_0_y_0_xyzzzz_1, g_0_y_0_xyzzzzz_0, g_0_y_0_xyzzzzz_1, g_0_y_0_yyyyyy_1, g_0_y_0_yyyyyyy_0, g_0_y_0_yyyyyyy_1, g_0_y_0_yyyyyyz_0, g_0_y_0_yyyyyyz_1, g_0_y_0_yyyyyz_1, g_0_y_0_yyyyyzz_0, g_0_y_0_yyyyyzz_1, g_0_y_0_yyyyzz_1, g_0_y_0_yyyyzzz_0, g_0_y_0_yyyyzzz_1, g_0_y_0_yyyzzz_1, g_0_y_0_yyyzzzz_0, g_0_y_0_yyyzzzz_1, g_0_y_0_yyzzzz_1, g_0_y_0_yyzzzzz_0, g_0_y_0_yyzzzzz_1, g_0_y_0_yzzzzz_1, g_0_y_0_yzzzzzz_0, g_0_y_0_yzzzzzz_1, g_0_y_0_zzzzzzz_0, g_0_y_0_zzzzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xy_0_xxxxxxx_0[i] = g_0_x_0_xxxxxxx_0[i] * pb_y + g_0_x_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xy_0_xxxxxxy_0[i] = 6.0 * g_0_y_0_xxxxxy_1[i] * fi_abcd_0 + g_0_y_0_xxxxxxy_0[i] * pb_x + g_0_y_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xy_0_xxxxxxz_0[i] = g_0_x_0_xxxxxxz_0[i] * pb_y + g_0_x_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xy_0_xxxxxyy_0[i] = 5.0 * g_0_y_0_xxxxyy_1[i] * fi_abcd_0 + g_0_y_0_xxxxxyy_0[i] * pb_x + g_0_y_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xy_0_xxxxxyz_0[i] = 5.0 * g_0_y_0_xxxxyz_1[i] * fi_abcd_0 + g_0_y_0_xxxxxyz_0[i] * pb_x + g_0_y_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xy_0_xxxxxzz_0[i] = g_0_x_0_xxxxxzz_0[i] * pb_y + g_0_x_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xy_0_xxxxyyy_0[i] = 4.0 * g_0_y_0_xxxyyy_1[i] * fi_abcd_0 + g_0_y_0_xxxxyyy_0[i] * pb_x + g_0_y_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xy_0_xxxxyyz_0[i] = 4.0 * g_0_y_0_xxxyyz_1[i] * fi_abcd_0 + g_0_y_0_xxxxyyz_0[i] * pb_x + g_0_y_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xy_0_xxxxyzz_0[i] = 4.0 * g_0_y_0_xxxyzz_1[i] * fi_abcd_0 + g_0_y_0_xxxxyzz_0[i] * pb_x + g_0_y_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xy_0_xxxxzzz_0[i] = g_0_x_0_xxxxzzz_0[i] * pb_y + g_0_x_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xy_0_xxxyyyy_0[i] = 3.0 * g_0_y_0_xxyyyy_1[i] * fi_abcd_0 + g_0_y_0_xxxyyyy_0[i] * pb_x + g_0_y_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xy_0_xxxyyyz_0[i] = 3.0 * g_0_y_0_xxyyyz_1[i] * fi_abcd_0 + g_0_y_0_xxxyyyz_0[i] * pb_x + g_0_y_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xy_0_xxxyyzz_0[i] = 3.0 * g_0_y_0_xxyyzz_1[i] * fi_abcd_0 + g_0_y_0_xxxyyzz_0[i] * pb_x + g_0_y_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xy_0_xxxyzzz_0[i] = 3.0 * g_0_y_0_xxyzzz_1[i] * fi_abcd_0 + g_0_y_0_xxxyzzz_0[i] * pb_x + g_0_y_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xy_0_xxxzzzz_0[i] = g_0_x_0_xxxzzzz_0[i] * pb_y + g_0_x_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xy_0_xxyyyyy_0[i] = 2.0 * g_0_y_0_xyyyyy_1[i] * fi_abcd_0 + g_0_y_0_xxyyyyy_0[i] * pb_x + g_0_y_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xy_0_xxyyyyz_0[i] = 2.0 * g_0_y_0_xyyyyz_1[i] * fi_abcd_0 + g_0_y_0_xxyyyyz_0[i] * pb_x + g_0_y_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xy_0_xxyyyzz_0[i] = 2.0 * g_0_y_0_xyyyzz_1[i] * fi_abcd_0 + g_0_y_0_xxyyyzz_0[i] * pb_x + g_0_y_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xy_0_xxyyzzz_0[i] = 2.0 * g_0_y_0_xyyzzz_1[i] * fi_abcd_0 + g_0_y_0_xxyyzzz_0[i] * pb_x + g_0_y_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xy_0_xxyzzzz_0[i] = 2.0 * g_0_y_0_xyzzzz_1[i] * fi_abcd_0 + g_0_y_0_xxyzzzz_0[i] * pb_x + g_0_y_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xy_0_xxzzzzz_0[i] = g_0_x_0_xxzzzzz_0[i] * pb_y + g_0_x_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xy_0_xyyyyyy_0[i] = g_0_y_0_yyyyyy_1[i] * fi_abcd_0 + g_0_y_0_xyyyyyy_0[i] * pb_x + g_0_y_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xy_0_xyyyyyz_0[i] = g_0_y_0_yyyyyz_1[i] * fi_abcd_0 + g_0_y_0_xyyyyyz_0[i] * pb_x + g_0_y_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xy_0_xyyyyzz_0[i] = g_0_y_0_yyyyzz_1[i] * fi_abcd_0 + g_0_y_0_xyyyyzz_0[i] * pb_x + g_0_y_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xy_0_xyyyzzz_0[i] = g_0_y_0_yyyzzz_1[i] * fi_abcd_0 + g_0_y_0_xyyyzzz_0[i] * pb_x + g_0_y_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xy_0_xyyzzzz_0[i] = g_0_y_0_yyzzzz_1[i] * fi_abcd_0 + g_0_y_0_xyyzzzz_0[i] * pb_x + g_0_y_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xy_0_xyzzzzz_0[i] = g_0_y_0_yzzzzz_1[i] * fi_abcd_0 + g_0_y_0_xyzzzzz_0[i] * pb_x + g_0_y_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xy_0_xzzzzzz_0[i] = g_0_x_0_xzzzzzz_0[i] * pb_y + g_0_x_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xy_0_yyyyyyy_0[i] = g_0_y_0_yyyyyyy_0[i] * pb_x + g_0_y_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xy_0_yyyyyyz_0[i] = g_0_y_0_yyyyyyz_0[i] * pb_x + g_0_y_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xy_0_yyyyyzz_0[i] = g_0_y_0_yyyyyzz_0[i] * pb_x + g_0_y_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xy_0_yyyyzzz_0[i] = g_0_y_0_yyyyzzz_0[i] * pb_x + g_0_y_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xy_0_yyyzzzz_0[i] = g_0_y_0_yyyzzzz_0[i] * pb_x + g_0_y_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xy_0_yyzzzzz_0[i] = g_0_y_0_yyzzzzz_0[i] * pb_x + g_0_y_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xy_0_yzzzzzz_0[i] = g_0_y_0_yzzzzzz_0[i] * pb_x + g_0_y_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xy_0_zzzzzzz_0[i] = g_0_y_0_zzzzzzz_0[i] * pb_x + g_0_y_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 72-108 components of targeted buffer : SDSK

    auto g_0_xz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sdsk + 72);

    auto g_0_xz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sdsk + 73);

    auto g_0_xz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sdsk + 74);

    auto g_0_xz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sdsk + 75);

    auto g_0_xz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sdsk + 76);

    auto g_0_xz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sdsk + 77);

    auto g_0_xz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sdsk + 78);

    auto g_0_xz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sdsk + 79);

    auto g_0_xz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sdsk + 80);

    auto g_0_xz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sdsk + 81);

    auto g_0_xz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sdsk + 82);

    auto g_0_xz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sdsk + 83);

    auto g_0_xz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sdsk + 84);

    auto g_0_xz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sdsk + 85);

    auto g_0_xz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 86);

    auto g_0_xz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sdsk + 87);

    auto g_0_xz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sdsk + 88);

    auto g_0_xz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sdsk + 89);

    auto g_0_xz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sdsk + 90);

    auto g_0_xz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 91);

    auto g_0_xz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 92);

    auto g_0_xz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sdsk + 93);

    auto g_0_xz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sdsk + 94);

    auto g_0_xz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sdsk + 95);

    auto g_0_xz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sdsk + 96);

    auto g_0_xz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 97);

    auto g_0_xz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 98);

    auto g_0_xz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 99);

    auto g_0_xz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sdsk + 100);

    auto g_0_xz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sdsk + 101);

    auto g_0_xz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sdsk + 102);

    auto g_0_xz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sdsk + 103);

    auto g_0_xz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 104);

    auto g_0_xz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 105);

    auto g_0_xz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 106);

    auto g_0_xz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 107);

    #pragma omp simd aligned(g_0_x_0_xxxxxxx_0, g_0_x_0_xxxxxxx_1, g_0_x_0_xxxxxxy_0, g_0_x_0_xxxxxxy_1, g_0_x_0_xxxxxyy_0, g_0_x_0_xxxxxyy_1, g_0_x_0_xxxxyyy_0, g_0_x_0_xxxxyyy_1, g_0_x_0_xxxyyyy_0, g_0_x_0_xxxyyyy_1, g_0_x_0_xxyyyyy_0, g_0_x_0_xxyyyyy_1, g_0_x_0_xyyyyyy_0, g_0_x_0_xyyyyyy_1, g_0_xz_0_xxxxxxx_0, g_0_xz_0_xxxxxxy_0, g_0_xz_0_xxxxxxz_0, g_0_xz_0_xxxxxyy_0, g_0_xz_0_xxxxxyz_0, g_0_xz_0_xxxxxzz_0, g_0_xz_0_xxxxyyy_0, g_0_xz_0_xxxxyyz_0, g_0_xz_0_xxxxyzz_0, g_0_xz_0_xxxxzzz_0, g_0_xz_0_xxxyyyy_0, g_0_xz_0_xxxyyyz_0, g_0_xz_0_xxxyyzz_0, g_0_xz_0_xxxyzzz_0, g_0_xz_0_xxxzzzz_0, g_0_xz_0_xxyyyyy_0, g_0_xz_0_xxyyyyz_0, g_0_xz_0_xxyyyzz_0, g_0_xz_0_xxyyzzz_0, g_0_xz_0_xxyzzzz_0, g_0_xz_0_xxzzzzz_0, g_0_xz_0_xyyyyyy_0, g_0_xz_0_xyyyyyz_0, g_0_xz_0_xyyyyzz_0, g_0_xz_0_xyyyzzz_0, g_0_xz_0_xyyzzzz_0, g_0_xz_0_xyzzzzz_0, g_0_xz_0_xzzzzzz_0, g_0_xz_0_yyyyyyy_0, g_0_xz_0_yyyyyyz_0, g_0_xz_0_yyyyyzz_0, g_0_xz_0_yyyyzzz_0, g_0_xz_0_yyyzzzz_0, g_0_xz_0_yyzzzzz_0, g_0_xz_0_yzzzzzz_0, g_0_xz_0_zzzzzzz_0, g_0_z_0_xxxxxxz_0, g_0_z_0_xxxxxxz_1, g_0_z_0_xxxxxyz_0, g_0_z_0_xxxxxyz_1, g_0_z_0_xxxxxz_1, g_0_z_0_xxxxxzz_0, g_0_z_0_xxxxxzz_1, g_0_z_0_xxxxyyz_0, g_0_z_0_xxxxyyz_1, g_0_z_0_xxxxyz_1, g_0_z_0_xxxxyzz_0, g_0_z_0_xxxxyzz_1, g_0_z_0_xxxxzz_1, g_0_z_0_xxxxzzz_0, g_0_z_0_xxxxzzz_1, g_0_z_0_xxxyyyz_0, g_0_z_0_xxxyyyz_1, g_0_z_0_xxxyyz_1, g_0_z_0_xxxyyzz_0, g_0_z_0_xxxyyzz_1, g_0_z_0_xxxyzz_1, g_0_z_0_xxxyzzz_0, g_0_z_0_xxxyzzz_1, g_0_z_0_xxxzzz_1, g_0_z_0_xxxzzzz_0, g_0_z_0_xxxzzzz_1, g_0_z_0_xxyyyyz_0, g_0_z_0_xxyyyyz_1, g_0_z_0_xxyyyz_1, g_0_z_0_xxyyyzz_0, g_0_z_0_xxyyyzz_1, g_0_z_0_xxyyzz_1, g_0_z_0_xxyyzzz_0, g_0_z_0_xxyyzzz_1, g_0_z_0_xxyzzz_1, g_0_z_0_xxyzzzz_0, g_0_z_0_xxyzzzz_1, g_0_z_0_xxzzzz_1, g_0_z_0_xxzzzzz_0, g_0_z_0_xxzzzzz_1, g_0_z_0_xyyyyyz_0, g_0_z_0_xyyyyyz_1, g_0_z_0_xyyyyz_1, g_0_z_0_xyyyyzz_0, g_0_z_0_xyyyyzz_1, g_0_z_0_xyyyzz_1, g_0_z_0_xyyyzzz_0, g_0_z_0_xyyyzzz_1, g_0_z_0_xyyzzz_1, g_0_z_0_xyyzzzz_0, g_0_z_0_xyyzzzz_1, g_0_z_0_xyzzzz_1, g_0_z_0_xyzzzzz_0, g_0_z_0_xyzzzzz_1, g_0_z_0_xzzzzz_1, g_0_z_0_xzzzzzz_0, g_0_z_0_xzzzzzz_1, g_0_z_0_yyyyyyy_0, g_0_z_0_yyyyyyy_1, g_0_z_0_yyyyyyz_0, g_0_z_0_yyyyyyz_1, g_0_z_0_yyyyyz_1, g_0_z_0_yyyyyzz_0, g_0_z_0_yyyyyzz_1, g_0_z_0_yyyyzz_1, g_0_z_0_yyyyzzz_0, g_0_z_0_yyyyzzz_1, g_0_z_0_yyyzzz_1, g_0_z_0_yyyzzzz_0, g_0_z_0_yyyzzzz_1, g_0_z_0_yyzzzz_1, g_0_z_0_yyzzzzz_0, g_0_z_0_yyzzzzz_1, g_0_z_0_yzzzzz_1, g_0_z_0_yzzzzzz_0, g_0_z_0_yzzzzzz_1, g_0_z_0_zzzzzz_1, g_0_z_0_zzzzzzz_0, g_0_z_0_zzzzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xz_0_xxxxxxx_0[i] = g_0_x_0_xxxxxxx_0[i] * pb_z + g_0_x_0_xxxxxxx_1[i] * wp_z[i];

        g_0_xz_0_xxxxxxy_0[i] = g_0_x_0_xxxxxxy_0[i] * pb_z + g_0_x_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xz_0_xxxxxxz_0[i] = 6.0 * g_0_z_0_xxxxxz_1[i] * fi_abcd_0 + g_0_z_0_xxxxxxz_0[i] * pb_x + g_0_z_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xz_0_xxxxxyy_0[i] = g_0_x_0_xxxxxyy_0[i] * pb_z + g_0_x_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xz_0_xxxxxyz_0[i] = 5.0 * g_0_z_0_xxxxyz_1[i] * fi_abcd_0 + g_0_z_0_xxxxxyz_0[i] * pb_x + g_0_z_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xz_0_xxxxxzz_0[i] = 5.0 * g_0_z_0_xxxxzz_1[i] * fi_abcd_0 + g_0_z_0_xxxxxzz_0[i] * pb_x + g_0_z_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xz_0_xxxxyyy_0[i] = g_0_x_0_xxxxyyy_0[i] * pb_z + g_0_x_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xz_0_xxxxyyz_0[i] = 4.0 * g_0_z_0_xxxyyz_1[i] * fi_abcd_0 + g_0_z_0_xxxxyyz_0[i] * pb_x + g_0_z_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xz_0_xxxxyzz_0[i] = 4.0 * g_0_z_0_xxxyzz_1[i] * fi_abcd_0 + g_0_z_0_xxxxyzz_0[i] * pb_x + g_0_z_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xz_0_xxxxzzz_0[i] = 4.0 * g_0_z_0_xxxzzz_1[i] * fi_abcd_0 + g_0_z_0_xxxxzzz_0[i] * pb_x + g_0_z_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xz_0_xxxyyyy_0[i] = g_0_x_0_xxxyyyy_0[i] * pb_z + g_0_x_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xz_0_xxxyyyz_0[i] = 3.0 * g_0_z_0_xxyyyz_1[i] * fi_abcd_0 + g_0_z_0_xxxyyyz_0[i] * pb_x + g_0_z_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xz_0_xxxyyzz_0[i] = 3.0 * g_0_z_0_xxyyzz_1[i] * fi_abcd_0 + g_0_z_0_xxxyyzz_0[i] * pb_x + g_0_z_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xz_0_xxxyzzz_0[i] = 3.0 * g_0_z_0_xxyzzz_1[i] * fi_abcd_0 + g_0_z_0_xxxyzzz_0[i] * pb_x + g_0_z_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xz_0_xxxzzzz_0[i] = 3.0 * g_0_z_0_xxzzzz_1[i] * fi_abcd_0 + g_0_z_0_xxxzzzz_0[i] * pb_x + g_0_z_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xz_0_xxyyyyy_0[i] = g_0_x_0_xxyyyyy_0[i] * pb_z + g_0_x_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xz_0_xxyyyyz_0[i] = 2.0 * g_0_z_0_xyyyyz_1[i] * fi_abcd_0 + g_0_z_0_xxyyyyz_0[i] * pb_x + g_0_z_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xz_0_xxyyyzz_0[i] = 2.0 * g_0_z_0_xyyyzz_1[i] * fi_abcd_0 + g_0_z_0_xxyyyzz_0[i] * pb_x + g_0_z_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xz_0_xxyyzzz_0[i] = 2.0 * g_0_z_0_xyyzzz_1[i] * fi_abcd_0 + g_0_z_0_xxyyzzz_0[i] * pb_x + g_0_z_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xz_0_xxyzzzz_0[i] = 2.0 * g_0_z_0_xyzzzz_1[i] * fi_abcd_0 + g_0_z_0_xxyzzzz_0[i] * pb_x + g_0_z_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xz_0_xxzzzzz_0[i] = 2.0 * g_0_z_0_xzzzzz_1[i] * fi_abcd_0 + g_0_z_0_xxzzzzz_0[i] * pb_x + g_0_z_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xz_0_xyyyyyy_0[i] = g_0_x_0_xyyyyyy_0[i] * pb_z + g_0_x_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xz_0_xyyyyyz_0[i] = g_0_z_0_yyyyyz_1[i] * fi_abcd_0 + g_0_z_0_xyyyyyz_0[i] * pb_x + g_0_z_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xz_0_xyyyyzz_0[i] = g_0_z_0_yyyyzz_1[i] * fi_abcd_0 + g_0_z_0_xyyyyzz_0[i] * pb_x + g_0_z_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xz_0_xyyyzzz_0[i] = g_0_z_0_yyyzzz_1[i] * fi_abcd_0 + g_0_z_0_xyyyzzz_0[i] * pb_x + g_0_z_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xz_0_xyyzzzz_0[i] = g_0_z_0_yyzzzz_1[i] * fi_abcd_0 + g_0_z_0_xyyzzzz_0[i] * pb_x + g_0_z_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xz_0_xyzzzzz_0[i] = g_0_z_0_yzzzzz_1[i] * fi_abcd_0 + g_0_z_0_xyzzzzz_0[i] * pb_x + g_0_z_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xz_0_xzzzzzz_0[i] = g_0_z_0_zzzzzz_1[i] * fi_abcd_0 + g_0_z_0_xzzzzzz_0[i] * pb_x + g_0_z_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xz_0_yyyyyyy_0[i] = g_0_z_0_yyyyyyy_0[i] * pb_x + g_0_z_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xz_0_yyyyyyz_0[i] = g_0_z_0_yyyyyyz_0[i] * pb_x + g_0_z_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xz_0_yyyyyzz_0[i] = g_0_z_0_yyyyyzz_0[i] * pb_x + g_0_z_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xz_0_yyyyzzz_0[i] = g_0_z_0_yyyyzzz_0[i] * pb_x + g_0_z_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xz_0_yyyzzzz_0[i] = g_0_z_0_yyyzzzz_0[i] * pb_x + g_0_z_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xz_0_yyzzzzz_0[i] = g_0_z_0_yyzzzzz_0[i] * pb_x + g_0_z_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xz_0_yzzzzzz_0[i] = g_0_z_0_yzzzzzz_0[i] * pb_x + g_0_z_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xz_0_zzzzzzz_0[i] = g_0_z_0_zzzzzzz_0[i] * pb_x + g_0_z_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 108-144 components of targeted buffer : SDSK

    auto g_0_yy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sdsk + 108);

    auto g_0_yy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sdsk + 109);

    auto g_0_yy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sdsk + 110);

    auto g_0_yy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sdsk + 111);

    auto g_0_yy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sdsk + 112);

    auto g_0_yy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sdsk + 113);

    auto g_0_yy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sdsk + 114);

    auto g_0_yy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sdsk + 115);

    auto g_0_yy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sdsk + 116);

    auto g_0_yy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sdsk + 117);

    auto g_0_yy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sdsk + 118);

    auto g_0_yy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sdsk + 119);

    auto g_0_yy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sdsk + 120);

    auto g_0_yy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sdsk + 121);

    auto g_0_yy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 122);

    auto g_0_yy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sdsk + 123);

    auto g_0_yy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sdsk + 124);

    auto g_0_yy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sdsk + 125);

    auto g_0_yy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sdsk + 126);

    auto g_0_yy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 127);

    auto g_0_yy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 128);

    auto g_0_yy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sdsk + 129);

    auto g_0_yy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sdsk + 130);

    auto g_0_yy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sdsk + 131);

    auto g_0_yy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sdsk + 132);

    auto g_0_yy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 133);

    auto g_0_yy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 134);

    auto g_0_yy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 135);

    auto g_0_yy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sdsk + 136);

    auto g_0_yy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sdsk + 137);

    auto g_0_yy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sdsk + 138);

    auto g_0_yy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sdsk + 139);

    auto g_0_yy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 140);

    auto g_0_yy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 141);

    auto g_0_yy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 142);

    auto g_0_yy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 143);

    #pragma omp simd aligned(g_0_0_0_xxxxxxx_0, g_0_0_0_xxxxxxx_1, g_0_0_0_xxxxxxy_0, g_0_0_0_xxxxxxy_1, g_0_0_0_xxxxxxz_0, g_0_0_0_xxxxxxz_1, g_0_0_0_xxxxxyy_0, g_0_0_0_xxxxxyy_1, g_0_0_0_xxxxxyz_0, g_0_0_0_xxxxxyz_1, g_0_0_0_xxxxxzz_0, g_0_0_0_xxxxxzz_1, g_0_0_0_xxxxyyy_0, g_0_0_0_xxxxyyy_1, g_0_0_0_xxxxyyz_0, g_0_0_0_xxxxyyz_1, g_0_0_0_xxxxyzz_0, g_0_0_0_xxxxyzz_1, g_0_0_0_xxxxzzz_0, g_0_0_0_xxxxzzz_1, g_0_0_0_xxxyyyy_0, g_0_0_0_xxxyyyy_1, g_0_0_0_xxxyyyz_0, g_0_0_0_xxxyyyz_1, g_0_0_0_xxxyyzz_0, g_0_0_0_xxxyyzz_1, g_0_0_0_xxxyzzz_0, g_0_0_0_xxxyzzz_1, g_0_0_0_xxxzzzz_0, g_0_0_0_xxxzzzz_1, g_0_0_0_xxyyyyy_0, g_0_0_0_xxyyyyy_1, g_0_0_0_xxyyyyz_0, g_0_0_0_xxyyyyz_1, g_0_0_0_xxyyyzz_0, g_0_0_0_xxyyyzz_1, g_0_0_0_xxyyzzz_0, g_0_0_0_xxyyzzz_1, g_0_0_0_xxyzzzz_0, g_0_0_0_xxyzzzz_1, g_0_0_0_xxzzzzz_0, g_0_0_0_xxzzzzz_1, g_0_0_0_xyyyyyy_0, g_0_0_0_xyyyyyy_1, g_0_0_0_xyyyyyz_0, g_0_0_0_xyyyyyz_1, g_0_0_0_xyyyyzz_0, g_0_0_0_xyyyyzz_1, g_0_0_0_xyyyzzz_0, g_0_0_0_xyyyzzz_1, g_0_0_0_xyyzzzz_0, g_0_0_0_xyyzzzz_1, g_0_0_0_xyzzzzz_0, g_0_0_0_xyzzzzz_1, g_0_0_0_xzzzzzz_0, g_0_0_0_xzzzzzz_1, g_0_0_0_yyyyyyy_0, g_0_0_0_yyyyyyy_1, g_0_0_0_yyyyyyz_0, g_0_0_0_yyyyyyz_1, g_0_0_0_yyyyyzz_0, g_0_0_0_yyyyyzz_1, g_0_0_0_yyyyzzz_0, g_0_0_0_yyyyzzz_1, g_0_0_0_yyyzzzz_0, g_0_0_0_yyyzzzz_1, g_0_0_0_yyzzzzz_0, g_0_0_0_yyzzzzz_1, g_0_0_0_yzzzzzz_0, g_0_0_0_yzzzzzz_1, g_0_0_0_zzzzzzz_0, g_0_0_0_zzzzzzz_1, g_0_y_0_xxxxxx_1, g_0_y_0_xxxxxxx_0, g_0_y_0_xxxxxxx_1, g_0_y_0_xxxxxxy_0, g_0_y_0_xxxxxxy_1, g_0_y_0_xxxxxxz_0, g_0_y_0_xxxxxxz_1, g_0_y_0_xxxxxy_1, g_0_y_0_xxxxxyy_0, g_0_y_0_xxxxxyy_1, g_0_y_0_xxxxxyz_0, g_0_y_0_xxxxxyz_1, g_0_y_0_xxxxxz_1, g_0_y_0_xxxxxzz_0, g_0_y_0_xxxxxzz_1, g_0_y_0_xxxxyy_1, g_0_y_0_xxxxyyy_0, g_0_y_0_xxxxyyy_1, g_0_y_0_xxxxyyz_0, g_0_y_0_xxxxyyz_1, g_0_y_0_xxxxyz_1, g_0_y_0_xxxxyzz_0, g_0_y_0_xxxxyzz_1, g_0_y_0_xxxxzz_1, g_0_y_0_xxxxzzz_0, g_0_y_0_xxxxzzz_1, g_0_y_0_xxxyyy_1, g_0_y_0_xxxyyyy_0, g_0_y_0_xxxyyyy_1, g_0_y_0_xxxyyyz_0, g_0_y_0_xxxyyyz_1, g_0_y_0_xxxyyz_1, g_0_y_0_xxxyyzz_0, g_0_y_0_xxxyyzz_1, g_0_y_0_xxxyzz_1, g_0_y_0_xxxyzzz_0, g_0_y_0_xxxyzzz_1, g_0_y_0_xxxzzz_1, g_0_y_0_xxxzzzz_0, g_0_y_0_xxxzzzz_1, g_0_y_0_xxyyyy_1, g_0_y_0_xxyyyyy_0, g_0_y_0_xxyyyyy_1, g_0_y_0_xxyyyyz_0, g_0_y_0_xxyyyyz_1, g_0_y_0_xxyyyz_1, g_0_y_0_xxyyyzz_0, g_0_y_0_xxyyyzz_1, g_0_y_0_xxyyzz_1, g_0_y_0_xxyyzzz_0, g_0_y_0_xxyyzzz_1, g_0_y_0_xxyzzz_1, g_0_y_0_xxyzzzz_0, g_0_y_0_xxyzzzz_1, g_0_y_0_xxzzzz_1, g_0_y_0_xxzzzzz_0, g_0_y_0_xxzzzzz_1, g_0_y_0_xyyyyy_1, g_0_y_0_xyyyyyy_0, g_0_y_0_xyyyyyy_1, g_0_y_0_xyyyyyz_0, g_0_y_0_xyyyyyz_1, g_0_y_0_xyyyyz_1, g_0_y_0_xyyyyzz_0, g_0_y_0_xyyyyzz_1, g_0_y_0_xyyyzz_1, g_0_y_0_xyyyzzz_0, g_0_y_0_xyyyzzz_1, g_0_y_0_xyyzzz_1, g_0_y_0_xyyzzzz_0, g_0_y_0_xyyzzzz_1, g_0_y_0_xyzzzz_1, g_0_y_0_xyzzzzz_0, g_0_y_0_xyzzzzz_1, g_0_y_0_xzzzzz_1, g_0_y_0_xzzzzzz_0, g_0_y_0_xzzzzzz_1, g_0_y_0_yyyyyy_1, g_0_y_0_yyyyyyy_0, g_0_y_0_yyyyyyy_1, g_0_y_0_yyyyyyz_0, g_0_y_0_yyyyyyz_1, g_0_y_0_yyyyyz_1, g_0_y_0_yyyyyzz_0, g_0_y_0_yyyyyzz_1, g_0_y_0_yyyyzz_1, g_0_y_0_yyyyzzz_0, g_0_y_0_yyyyzzz_1, g_0_y_0_yyyzzz_1, g_0_y_0_yyyzzzz_0, g_0_y_0_yyyzzzz_1, g_0_y_0_yyzzzz_1, g_0_y_0_yyzzzzz_0, g_0_y_0_yyzzzzz_1, g_0_y_0_yzzzzz_1, g_0_y_0_yzzzzzz_0, g_0_y_0_yzzzzzz_1, g_0_y_0_zzzzzz_1, g_0_y_0_zzzzzzz_0, g_0_y_0_zzzzzzz_1, g_0_yy_0_xxxxxxx_0, g_0_yy_0_xxxxxxy_0, g_0_yy_0_xxxxxxz_0, g_0_yy_0_xxxxxyy_0, g_0_yy_0_xxxxxyz_0, g_0_yy_0_xxxxxzz_0, g_0_yy_0_xxxxyyy_0, g_0_yy_0_xxxxyyz_0, g_0_yy_0_xxxxyzz_0, g_0_yy_0_xxxxzzz_0, g_0_yy_0_xxxyyyy_0, g_0_yy_0_xxxyyyz_0, g_0_yy_0_xxxyyzz_0, g_0_yy_0_xxxyzzz_0, g_0_yy_0_xxxzzzz_0, g_0_yy_0_xxyyyyy_0, g_0_yy_0_xxyyyyz_0, g_0_yy_0_xxyyyzz_0, g_0_yy_0_xxyyzzz_0, g_0_yy_0_xxyzzzz_0, g_0_yy_0_xxzzzzz_0, g_0_yy_0_xyyyyyy_0, g_0_yy_0_xyyyyyz_0, g_0_yy_0_xyyyyzz_0, g_0_yy_0_xyyyzzz_0, g_0_yy_0_xyyzzzz_0, g_0_yy_0_xyzzzzz_0, g_0_yy_0_xzzzzzz_0, g_0_yy_0_yyyyyyy_0, g_0_yy_0_yyyyyyz_0, g_0_yy_0_yyyyyzz_0, g_0_yy_0_yyyyzzz_0, g_0_yy_0_yyyzzzz_0, g_0_yy_0_yyzzzzz_0, g_0_yy_0_yzzzzzz_0, g_0_yy_0_zzzzzzz_0, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yy_0_xxxxxxx_0[i] = g_0_0_0_xxxxxxx_0[i] * fi_ab_0 - g_0_0_0_xxxxxxx_1[i] * fti_ab_0 + g_0_y_0_xxxxxxx_0[i] * pb_y + g_0_y_0_xxxxxxx_1[i] * wp_y[i];

        g_0_yy_0_xxxxxxy_0[i] = g_0_0_0_xxxxxxy_0[i] * fi_ab_0 - g_0_0_0_xxxxxxy_1[i] * fti_ab_0 + g_0_y_0_xxxxxx_1[i] * fi_abcd_0 + g_0_y_0_xxxxxxy_0[i] * pb_y + g_0_y_0_xxxxxxy_1[i] * wp_y[i];

        g_0_yy_0_xxxxxxz_0[i] = g_0_0_0_xxxxxxz_0[i] * fi_ab_0 - g_0_0_0_xxxxxxz_1[i] * fti_ab_0 + g_0_y_0_xxxxxxz_0[i] * pb_y + g_0_y_0_xxxxxxz_1[i] * wp_y[i];

        g_0_yy_0_xxxxxyy_0[i] = g_0_0_0_xxxxxyy_0[i] * fi_ab_0 - g_0_0_0_xxxxxyy_1[i] * fti_ab_0 + 2.0 * g_0_y_0_xxxxxy_1[i] * fi_abcd_0 + g_0_y_0_xxxxxyy_0[i] * pb_y + g_0_y_0_xxxxxyy_1[i] * wp_y[i];

        g_0_yy_0_xxxxxyz_0[i] = g_0_0_0_xxxxxyz_0[i] * fi_ab_0 - g_0_0_0_xxxxxyz_1[i] * fti_ab_0 + g_0_y_0_xxxxxz_1[i] * fi_abcd_0 + g_0_y_0_xxxxxyz_0[i] * pb_y + g_0_y_0_xxxxxyz_1[i] * wp_y[i];

        g_0_yy_0_xxxxxzz_0[i] = g_0_0_0_xxxxxzz_0[i] * fi_ab_0 - g_0_0_0_xxxxxzz_1[i] * fti_ab_0 + g_0_y_0_xxxxxzz_0[i] * pb_y + g_0_y_0_xxxxxzz_1[i] * wp_y[i];

        g_0_yy_0_xxxxyyy_0[i] = g_0_0_0_xxxxyyy_0[i] * fi_ab_0 - g_0_0_0_xxxxyyy_1[i] * fti_ab_0 + 3.0 * g_0_y_0_xxxxyy_1[i] * fi_abcd_0 + g_0_y_0_xxxxyyy_0[i] * pb_y + g_0_y_0_xxxxyyy_1[i] * wp_y[i];

        g_0_yy_0_xxxxyyz_0[i] = g_0_0_0_xxxxyyz_0[i] * fi_ab_0 - g_0_0_0_xxxxyyz_1[i] * fti_ab_0 + 2.0 * g_0_y_0_xxxxyz_1[i] * fi_abcd_0 + g_0_y_0_xxxxyyz_0[i] * pb_y + g_0_y_0_xxxxyyz_1[i] * wp_y[i];

        g_0_yy_0_xxxxyzz_0[i] = g_0_0_0_xxxxyzz_0[i] * fi_ab_0 - g_0_0_0_xxxxyzz_1[i] * fti_ab_0 + g_0_y_0_xxxxzz_1[i] * fi_abcd_0 + g_0_y_0_xxxxyzz_0[i] * pb_y + g_0_y_0_xxxxyzz_1[i] * wp_y[i];

        g_0_yy_0_xxxxzzz_0[i] = g_0_0_0_xxxxzzz_0[i] * fi_ab_0 - g_0_0_0_xxxxzzz_1[i] * fti_ab_0 + g_0_y_0_xxxxzzz_0[i] * pb_y + g_0_y_0_xxxxzzz_1[i] * wp_y[i];

        g_0_yy_0_xxxyyyy_0[i] = g_0_0_0_xxxyyyy_0[i] * fi_ab_0 - g_0_0_0_xxxyyyy_1[i] * fti_ab_0 + 4.0 * g_0_y_0_xxxyyy_1[i] * fi_abcd_0 + g_0_y_0_xxxyyyy_0[i] * pb_y + g_0_y_0_xxxyyyy_1[i] * wp_y[i];

        g_0_yy_0_xxxyyyz_0[i] = g_0_0_0_xxxyyyz_0[i] * fi_ab_0 - g_0_0_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_y_0_xxxyyz_1[i] * fi_abcd_0 + g_0_y_0_xxxyyyz_0[i] * pb_y + g_0_y_0_xxxyyyz_1[i] * wp_y[i];

        g_0_yy_0_xxxyyzz_0[i] = g_0_0_0_xxxyyzz_0[i] * fi_ab_0 - g_0_0_0_xxxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_y_0_xxxyzz_1[i] * fi_abcd_0 + g_0_y_0_xxxyyzz_0[i] * pb_y + g_0_y_0_xxxyyzz_1[i] * wp_y[i];

        g_0_yy_0_xxxyzzz_0[i] = g_0_0_0_xxxyzzz_0[i] * fi_ab_0 - g_0_0_0_xxxyzzz_1[i] * fti_ab_0 + g_0_y_0_xxxzzz_1[i] * fi_abcd_0 + g_0_y_0_xxxyzzz_0[i] * pb_y + g_0_y_0_xxxyzzz_1[i] * wp_y[i];

        g_0_yy_0_xxxzzzz_0[i] = g_0_0_0_xxxzzzz_0[i] * fi_ab_0 - g_0_0_0_xxxzzzz_1[i] * fti_ab_0 + g_0_y_0_xxxzzzz_0[i] * pb_y + g_0_y_0_xxxzzzz_1[i] * wp_y[i];

        g_0_yy_0_xxyyyyy_0[i] = g_0_0_0_xxyyyyy_0[i] * fi_ab_0 - g_0_0_0_xxyyyyy_1[i] * fti_ab_0 + 5.0 * g_0_y_0_xxyyyy_1[i] * fi_abcd_0 + g_0_y_0_xxyyyyy_0[i] * pb_y + g_0_y_0_xxyyyyy_1[i] * wp_y[i];

        g_0_yy_0_xxyyyyz_0[i] = g_0_0_0_xxyyyyz_0[i] * fi_ab_0 - g_0_0_0_xxyyyyz_1[i] * fti_ab_0 + 4.0 * g_0_y_0_xxyyyz_1[i] * fi_abcd_0 + g_0_y_0_xxyyyyz_0[i] * pb_y + g_0_y_0_xxyyyyz_1[i] * wp_y[i];

        g_0_yy_0_xxyyyzz_0[i] = g_0_0_0_xxyyyzz_0[i] * fi_ab_0 - g_0_0_0_xxyyyzz_1[i] * fti_ab_0 + 3.0 * g_0_y_0_xxyyzz_1[i] * fi_abcd_0 + g_0_y_0_xxyyyzz_0[i] * pb_y + g_0_y_0_xxyyyzz_1[i] * wp_y[i];

        g_0_yy_0_xxyyzzz_0[i] = g_0_0_0_xxyyzzz_0[i] * fi_ab_0 - g_0_0_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_y_0_xxyzzz_1[i] * fi_abcd_0 + g_0_y_0_xxyyzzz_0[i] * pb_y + g_0_y_0_xxyyzzz_1[i] * wp_y[i];

        g_0_yy_0_xxyzzzz_0[i] = g_0_0_0_xxyzzzz_0[i] * fi_ab_0 - g_0_0_0_xxyzzzz_1[i] * fti_ab_0 + g_0_y_0_xxzzzz_1[i] * fi_abcd_0 + g_0_y_0_xxyzzzz_0[i] * pb_y + g_0_y_0_xxyzzzz_1[i] * wp_y[i];

        g_0_yy_0_xxzzzzz_0[i] = g_0_0_0_xxzzzzz_0[i] * fi_ab_0 - g_0_0_0_xxzzzzz_1[i] * fti_ab_0 + g_0_y_0_xxzzzzz_0[i] * pb_y + g_0_y_0_xxzzzzz_1[i] * wp_y[i];

        g_0_yy_0_xyyyyyy_0[i] = g_0_0_0_xyyyyyy_0[i] * fi_ab_0 - g_0_0_0_xyyyyyy_1[i] * fti_ab_0 + 6.0 * g_0_y_0_xyyyyy_1[i] * fi_abcd_0 + g_0_y_0_xyyyyyy_0[i] * pb_y + g_0_y_0_xyyyyyy_1[i] * wp_y[i];

        g_0_yy_0_xyyyyyz_0[i] = g_0_0_0_xyyyyyz_0[i] * fi_ab_0 - g_0_0_0_xyyyyyz_1[i] * fti_ab_0 + 5.0 * g_0_y_0_xyyyyz_1[i] * fi_abcd_0 + g_0_y_0_xyyyyyz_0[i] * pb_y + g_0_y_0_xyyyyyz_1[i] * wp_y[i];

        g_0_yy_0_xyyyyzz_0[i] = g_0_0_0_xyyyyzz_0[i] * fi_ab_0 - g_0_0_0_xyyyyzz_1[i] * fti_ab_0 + 4.0 * g_0_y_0_xyyyzz_1[i] * fi_abcd_0 + g_0_y_0_xyyyyzz_0[i] * pb_y + g_0_y_0_xyyyyzz_1[i] * wp_y[i];

        g_0_yy_0_xyyyzzz_0[i] = g_0_0_0_xyyyzzz_0[i] * fi_ab_0 - g_0_0_0_xyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_y_0_xyyzzz_1[i] * fi_abcd_0 + g_0_y_0_xyyyzzz_0[i] * pb_y + g_0_y_0_xyyyzzz_1[i] * wp_y[i];

        g_0_yy_0_xyyzzzz_0[i] = g_0_0_0_xyyzzzz_0[i] * fi_ab_0 - g_0_0_0_xyyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_y_0_xyzzzz_1[i] * fi_abcd_0 + g_0_y_0_xyyzzzz_0[i] * pb_y + g_0_y_0_xyyzzzz_1[i] * wp_y[i];

        g_0_yy_0_xyzzzzz_0[i] = g_0_0_0_xyzzzzz_0[i] * fi_ab_0 - g_0_0_0_xyzzzzz_1[i] * fti_ab_0 + g_0_y_0_xzzzzz_1[i] * fi_abcd_0 + g_0_y_0_xyzzzzz_0[i] * pb_y + g_0_y_0_xyzzzzz_1[i] * wp_y[i];

        g_0_yy_0_xzzzzzz_0[i] = g_0_0_0_xzzzzzz_0[i] * fi_ab_0 - g_0_0_0_xzzzzzz_1[i] * fti_ab_0 + g_0_y_0_xzzzzzz_0[i] * pb_y + g_0_y_0_xzzzzzz_1[i] * wp_y[i];

        g_0_yy_0_yyyyyyy_0[i] = g_0_0_0_yyyyyyy_0[i] * fi_ab_0 - g_0_0_0_yyyyyyy_1[i] * fti_ab_0 + 7.0 * g_0_y_0_yyyyyy_1[i] * fi_abcd_0 + g_0_y_0_yyyyyyy_0[i] * pb_y + g_0_y_0_yyyyyyy_1[i] * wp_y[i];

        g_0_yy_0_yyyyyyz_0[i] = g_0_0_0_yyyyyyz_0[i] * fi_ab_0 - g_0_0_0_yyyyyyz_1[i] * fti_ab_0 + 6.0 * g_0_y_0_yyyyyz_1[i] * fi_abcd_0 + g_0_y_0_yyyyyyz_0[i] * pb_y + g_0_y_0_yyyyyyz_1[i] * wp_y[i];

        g_0_yy_0_yyyyyzz_0[i] = g_0_0_0_yyyyyzz_0[i] * fi_ab_0 - g_0_0_0_yyyyyzz_1[i] * fti_ab_0 + 5.0 * g_0_y_0_yyyyzz_1[i] * fi_abcd_0 + g_0_y_0_yyyyyzz_0[i] * pb_y + g_0_y_0_yyyyyzz_1[i] * wp_y[i];

        g_0_yy_0_yyyyzzz_0[i] = g_0_0_0_yyyyzzz_0[i] * fi_ab_0 - g_0_0_0_yyyyzzz_1[i] * fti_ab_0 + 4.0 * g_0_y_0_yyyzzz_1[i] * fi_abcd_0 + g_0_y_0_yyyyzzz_0[i] * pb_y + g_0_y_0_yyyyzzz_1[i] * wp_y[i];

        g_0_yy_0_yyyzzzz_0[i] = g_0_0_0_yyyzzzz_0[i] * fi_ab_0 - g_0_0_0_yyyzzzz_1[i] * fti_ab_0 + 3.0 * g_0_y_0_yyzzzz_1[i] * fi_abcd_0 + g_0_y_0_yyyzzzz_0[i] * pb_y + g_0_y_0_yyyzzzz_1[i] * wp_y[i];

        g_0_yy_0_yyzzzzz_0[i] = g_0_0_0_yyzzzzz_0[i] * fi_ab_0 - g_0_0_0_yyzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_y_0_yzzzzz_1[i] * fi_abcd_0 + g_0_y_0_yyzzzzz_0[i] * pb_y + g_0_y_0_yyzzzzz_1[i] * wp_y[i];

        g_0_yy_0_yzzzzzz_0[i] = g_0_0_0_yzzzzzz_0[i] * fi_ab_0 - g_0_0_0_yzzzzzz_1[i] * fti_ab_0 + g_0_y_0_zzzzzz_1[i] * fi_abcd_0 + g_0_y_0_yzzzzzz_0[i] * pb_y + g_0_y_0_yzzzzzz_1[i] * wp_y[i];

        g_0_yy_0_zzzzzzz_0[i] = g_0_0_0_zzzzzzz_0[i] * fi_ab_0 - g_0_0_0_zzzzzzz_1[i] * fti_ab_0 + g_0_y_0_zzzzzzz_0[i] * pb_y + g_0_y_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 144-180 components of targeted buffer : SDSK

    auto g_0_yz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sdsk + 144);

    auto g_0_yz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sdsk + 145);

    auto g_0_yz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sdsk + 146);

    auto g_0_yz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sdsk + 147);

    auto g_0_yz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sdsk + 148);

    auto g_0_yz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sdsk + 149);

    auto g_0_yz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sdsk + 150);

    auto g_0_yz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sdsk + 151);

    auto g_0_yz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sdsk + 152);

    auto g_0_yz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sdsk + 153);

    auto g_0_yz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sdsk + 154);

    auto g_0_yz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sdsk + 155);

    auto g_0_yz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sdsk + 156);

    auto g_0_yz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sdsk + 157);

    auto g_0_yz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 158);

    auto g_0_yz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sdsk + 159);

    auto g_0_yz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sdsk + 160);

    auto g_0_yz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sdsk + 161);

    auto g_0_yz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sdsk + 162);

    auto g_0_yz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 163);

    auto g_0_yz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 164);

    auto g_0_yz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sdsk + 165);

    auto g_0_yz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sdsk + 166);

    auto g_0_yz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sdsk + 167);

    auto g_0_yz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sdsk + 168);

    auto g_0_yz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 169);

    auto g_0_yz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 170);

    auto g_0_yz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 171);

    auto g_0_yz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sdsk + 172);

    auto g_0_yz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sdsk + 173);

    auto g_0_yz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sdsk + 174);

    auto g_0_yz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sdsk + 175);

    auto g_0_yz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 176);

    auto g_0_yz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 177);

    auto g_0_yz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 178);

    auto g_0_yz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 179);

    #pragma omp simd aligned(g_0_y_0_xxxxxxy_0, g_0_y_0_xxxxxxy_1, g_0_y_0_xxxxxyy_0, g_0_y_0_xxxxxyy_1, g_0_y_0_xxxxyyy_0, g_0_y_0_xxxxyyy_1, g_0_y_0_xxxyyyy_0, g_0_y_0_xxxyyyy_1, g_0_y_0_xxyyyyy_0, g_0_y_0_xxyyyyy_1, g_0_y_0_xyyyyyy_0, g_0_y_0_xyyyyyy_1, g_0_y_0_yyyyyyy_0, g_0_y_0_yyyyyyy_1, g_0_yz_0_xxxxxxx_0, g_0_yz_0_xxxxxxy_0, g_0_yz_0_xxxxxxz_0, g_0_yz_0_xxxxxyy_0, g_0_yz_0_xxxxxyz_0, g_0_yz_0_xxxxxzz_0, g_0_yz_0_xxxxyyy_0, g_0_yz_0_xxxxyyz_0, g_0_yz_0_xxxxyzz_0, g_0_yz_0_xxxxzzz_0, g_0_yz_0_xxxyyyy_0, g_0_yz_0_xxxyyyz_0, g_0_yz_0_xxxyyzz_0, g_0_yz_0_xxxyzzz_0, g_0_yz_0_xxxzzzz_0, g_0_yz_0_xxyyyyy_0, g_0_yz_0_xxyyyyz_0, g_0_yz_0_xxyyyzz_0, g_0_yz_0_xxyyzzz_0, g_0_yz_0_xxyzzzz_0, g_0_yz_0_xxzzzzz_0, g_0_yz_0_xyyyyyy_0, g_0_yz_0_xyyyyyz_0, g_0_yz_0_xyyyyzz_0, g_0_yz_0_xyyyzzz_0, g_0_yz_0_xyyzzzz_0, g_0_yz_0_xyzzzzz_0, g_0_yz_0_xzzzzzz_0, g_0_yz_0_yyyyyyy_0, g_0_yz_0_yyyyyyz_0, g_0_yz_0_yyyyyzz_0, g_0_yz_0_yyyyzzz_0, g_0_yz_0_yyyzzzz_0, g_0_yz_0_yyzzzzz_0, g_0_yz_0_yzzzzzz_0, g_0_yz_0_zzzzzzz_0, g_0_z_0_xxxxxxx_0, g_0_z_0_xxxxxxx_1, g_0_z_0_xxxxxxz_0, g_0_z_0_xxxxxxz_1, g_0_z_0_xxxxxyz_0, g_0_z_0_xxxxxyz_1, g_0_z_0_xxxxxz_1, g_0_z_0_xxxxxzz_0, g_0_z_0_xxxxxzz_1, g_0_z_0_xxxxyyz_0, g_0_z_0_xxxxyyz_1, g_0_z_0_xxxxyz_1, g_0_z_0_xxxxyzz_0, g_0_z_0_xxxxyzz_1, g_0_z_0_xxxxzz_1, g_0_z_0_xxxxzzz_0, g_0_z_0_xxxxzzz_1, g_0_z_0_xxxyyyz_0, g_0_z_0_xxxyyyz_1, g_0_z_0_xxxyyz_1, g_0_z_0_xxxyyzz_0, g_0_z_0_xxxyyzz_1, g_0_z_0_xxxyzz_1, g_0_z_0_xxxyzzz_0, g_0_z_0_xxxyzzz_1, g_0_z_0_xxxzzz_1, g_0_z_0_xxxzzzz_0, g_0_z_0_xxxzzzz_1, g_0_z_0_xxyyyyz_0, g_0_z_0_xxyyyyz_1, g_0_z_0_xxyyyz_1, g_0_z_0_xxyyyzz_0, g_0_z_0_xxyyyzz_1, g_0_z_0_xxyyzz_1, g_0_z_0_xxyyzzz_0, g_0_z_0_xxyyzzz_1, g_0_z_0_xxyzzz_1, g_0_z_0_xxyzzzz_0, g_0_z_0_xxyzzzz_1, g_0_z_0_xxzzzz_1, g_0_z_0_xxzzzzz_0, g_0_z_0_xxzzzzz_1, g_0_z_0_xyyyyyz_0, g_0_z_0_xyyyyyz_1, g_0_z_0_xyyyyz_1, g_0_z_0_xyyyyzz_0, g_0_z_0_xyyyyzz_1, g_0_z_0_xyyyzz_1, g_0_z_0_xyyyzzz_0, g_0_z_0_xyyyzzz_1, g_0_z_0_xyyzzz_1, g_0_z_0_xyyzzzz_0, g_0_z_0_xyyzzzz_1, g_0_z_0_xyzzzz_1, g_0_z_0_xyzzzzz_0, g_0_z_0_xyzzzzz_1, g_0_z_0_xzzzzz_1, g_0_z_0_xzzzzzz_0, g_0_z_0_xzzzzzz_1, g_0_z_0_yyyyyyz_0, g_0_z_0_yyyyyyz_1, g_0_z_0_yyyyyz_1, g_0_z_0_yyyyyzz_0, g_0_z_0_yyyyyzz_1, g_0_z_0_yyyyzz_1, g_0_z_0_yyyyzzz_0, g_0_z_0_yyyyzzz_1, g_0_z_0_yyyzzz_1, g_0_z_0_yyyzzzz_0, g_0_z_0_yyyzzzz_1, g_0_z_0_yyzzzz_1, g_0_z_0_yyzzzzz_0, g_0_z_0_yyzzzzz_1, g_0_z_0_yzzzzz_1, g_0_z_0_yzzzzzz_0, g_0_z_0_yzzzzzz_1, g_0_z_0_zzzzzz_1, g_0_z_0_zzzzzzz_0, g_0_z_0_zzzzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yz_0_xxxxxxx_0[i] = g_0_z_0_xxxxxxx_0[i] * pb_y + g_0_z_0_xxxxxxx_1[i] * wp_y[i];

        g_0_yz_0_xxxxxxy_0[i] = g_0_y_0_xxxxxxy_0[i] * pb_z + g_0_y_0_xxxxxxy_1[i] * wp_z[i];

        g_0_yz_0_xxxxxxz_0[i] = g_0_z_0_xxxxxxz_0[i] * pb_y + g_0_z_0_xxxxxxz_1[i] * wp_y[i];

        g_0_yz_0_xxxxxyy_0[i] = g_0_y_0_xxxxxyy_0[i] * pb_z + g_0_y_0_xxxxxyy_1[i] * wp_z[i];

        g_0_yz_0_xxxxxyz_0[i] = g_0_z_0_xxxxxz_1[i] * fi_abcd_0 + g_0_z_0_xxxxxyz_0[i] * pb_y + g_0_z_0_xxxxxyz_1[i] * wp_y[i];

        g_0_yz_0_xxxxxzz_0[i] = g_0_z_0_xxxxxzz_0[i] * pb_y + g_0_z_0_xxxxxzz_1[i] * wp_y[i];

        g_0_yz_0_xxxxyyy_0[i] = g_0_y_0_xxxxyyy_0[i] * pb_z + g_0_y_0_xxxxyyy_1[i] * wp_z[i];

        g_0_yz_0_xxxxyyz_0[i] = 2.0 * g_0_z_0_xxxxyz_1[i] * fi_abcd_0 + g_0_z_0_xxxxyyz_0[i] * pb_y + g_0_z_0_xxxxyyz_1[i] * wp_y[i];

        g_0_yz_0_xxxxyzz_0[i] = g_0_z_0_xxxxzz_1[i] * fi_abcd_0 + g_0_z_0_xxxxyzz_0[i] * pb_y + g_0_z_0_xxxxyzz_1[i] * wp_y[i];

        g_0_yz_0_xxxxzzz_0[i] = g_0_z_0_xxxxzzz_0[i] * pb_y + g_0_z_0_xxxxzzz_1[i] * wp_y[i];

        g_0_yz_0_xxxyyyy_0[i] = g_0_y_0_xxxyyyy_0[i] * pb_z + g_0_y_0_xxxyyyy_1[i] * wp_z[i];

        g_0_yz_0_xxxyyyz_0[i] = 3.0 * g_0_z_0_xxxyyz_1[i] * fi_abcd_0 + g_0_z_0_xxxyyyz_0[i] * pb_y + g_0_z_0_xxxyyyz_1[i] * wp_y[i];

        g_0_yz_0_xxxyyzz_0[i] = 2.0 * g_0_z_0_xxxyzz_1[i] * fi_abcd_0 + g_0_z_0_xxxyyzz_0[i] * pb_y + g_0_z_0_xxxyyzz_1[i] * wp_y[i];

        g_0_yz_0_xxxyzzz_0[i] = g_0_z_0_xxxzzz_1[i] * fi_abcd_0 + g_0_z_0_xxxyzzz_0[i] * pb_y + g_0_z_0_xxxyzzz_1[i] * wp_y[i];

        g_0_yz_0_xxxzzzz_0[i] = g_0_z_0_xxxzzzz_0[i] * pb_y + g_0_z_0_xxxzzzz_1[i] * wp_y[i];

        g_0_yz_0_xxyyyyy_0[i] = g_0_y_0_xxyyyyy_0[i] * pb_z + g_0_y_0_xxyyyyy_1[i] * wp_z[i];

        g_0_yz_0_xxyyyyz_0[i] = 4.0 * g_0_z_0_xxyyyz_1[i] * fi_abcd_0 + g_0_z_0_xxyyyyz_0[i] * pb_y + g_0_z_0_xxyyyyz_1[i] * wp_y[i];

        g_0_yz_0_xxyyyzz_0[i] = 3.0 * g_0_z_0_xxyyzz_1[i] * fi_abcd_0 + g_0_z_0_xxyyyzz_0[i] * pb_y + g_0_z_0_xxyyyzz_1[i] * wp_y[i];

        g_0_yz_0_xxyyzzz_0[i] = 2.0 * g_0_z_0_xxyzzz_1[i] * fi_abcd_0 + g_0_z_0_xxyyzzz_0[i] * pb_y + g_0_z_0_xxyyzzz_1[i] * wp_y[i];

        g_0_yz_0_xxyzzzz_0[i] = g_0_z_0_xxzzzz_1[i] * fi_abcd_0 + g_0_z_0_xxyzzzz_0[i] * pb_y + g_0_z_0_xxyzzzz_1[i] * wp_y[i];

        g_0_yz_0_xxzzzzz_0[i] = g_0_z_0_xxzzzzz_0[i] * pb_y + g_0_z_0_xxzzzzz_1[i] * wp_y[i];

        g_0_yz_0_xyyyyyy_0[i] = g_0_y_0_xyyyyyy_0[i] * pb_z + g_0_y_0_xyyyyyy_1[i] * wp_z[i];

        g_0_yz_0_xyyyyyz_0[i] = 5.0 * g_0_z_0_xyyyyz_1[i] * fi_abcd_0 + g_0_z_0_xyyyyyz_0[i] * pb_y + g_0_z_0_xyyyyyz_1[i] * wp_y[i];

        g_0_yz_0_xyyyyzz_0[i] = 4.0 * g_0_z_0_xyyyzz_1[i] * fi_abcd_0 + g_0_z_0_xyyyyzz_0[i] * pb_y + g_0_z_0_xyyyyzz_1[i] * wp_y[i];

        g_0_yz_0_xyyyzzz_0[i] = 3.0 * g_0_z_0_xyyzzz_1[i] * fi_abcd_0 + g_0_z_0_xyyyzzz_0[i] * pb_y + g_0_z_0_xyyyzzz_1[i] * wp_y[i];

        g_0_yz_0_xyyzzzz_0[i] = 2.0 * g_0_z_0_xyzzzz_1[i] * fi_abcd_0 + g_0_z_0_xyyzzzz_0[i] * pb_y + g_0_z_0_xyyzzzz_1[i] * wp_y[i];

        g_0_yz_0_xyzzzzz_0[i] = g_0_z_0_xzzzzz_1[i] * fi_abcd_0 + g_0_z_0_xyzzzzz_0[i] * pb_y + g_0_z_0_xyzzzzz_1[i] * wp_y[i];

        g_0_yz_0_xzzzzzz_0[i] = g_0_z_0_xzzzzzz_0[i] * pb_y + g_0_z_0_xzzzzzz_1[i] * wp_y[i];

        g_0_yz_0_yyyyyyy_0[i] = g_0_y_0_yyyyyyy_0[i] * pb_z + g_0_y_0_yyyyyyy_1[i] * wp_z[i];

        g_0_yz_0_yyyyyyz_0[i] = 6.0 * g_0_z_0_yyyyyz_1[i] * fi_abcd_0 + g_0_z_0_yyyyyyz_0[i] * pb_y + g_0_z_0_yyyyyyz_1[i] * wp_y[i];

        g_0_yz_0_yyyyyzz_0[i] = 5.0 * g_0_z_0_yyyyzz_1[i] * fi_abcd_0 + g_0_z_0_yyyyyzz_0[i] * pb_y + g_0_z_0_yyyyyzz_1[i] * wp_y[i];

        g_0_yz_0_yyyyzzz_0[i] = 4.0 * g_0_z_0_yyyzzz_1[i] * fi_abcd_0 + g_0_z_0_yyyyzzz_0[i] * pb_y + g_0_z_0_yyyyzzz_1[i] * wp_y[i];

        g_0_yz_0_yyyzzzz_0[i] = 3.0 * g_0_z_0_yyzzzz_1[i] * fi_abcd_0 + g_0_z_0_yyyzzzz_0[i] * pb_y + g_0_z_0_yyyzzzz_1[i] * wp_y[i];

        g_0_yz_0_yyzzzzz_0[i] = 2.0 * g_0_z_0_yzzzzz_1[i] * fi_abcd_0 + g_0_z_0_yyzzzzz_0[i] * pb_y + g_0_z_0_yyzzzzz_1[i] * wp_y[i];

        g_0_yz_0_yzzzzzz_0[i] = g_0_z_0_zzzzzz_1[i] * fi_abcd_0 + g_0_z_0_yzzzzzz_0[i] * pb_y + g_0_z_0_yzzzzzz_1[i] * wp_y[i];

        g_0_yz_0_zzzzzzz_0[i] = g_0_z_0_zzzzzzz_0[i] * pb_y + g_0_z_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 180-216 components of targeted buffer : SDSK

    auto g_0_zz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sdsk + 180);

    auto g_0_zz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sdsk + 181);

    auto g_0_zz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sdsk + 182);

    auto g_0_zz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sdsk + 183);

    auto g_0_zz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sdsk + 184);

    auto g_0_zz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sdsk + 185);

    auto g_0_zz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sdsk + 186);

    auto g_0_zz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sdsk + 187);

    auto g_0_zz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sdsk + 188);

    auto g_0_zz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sdsk + 189);

    auto g_0_zz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sdsk + 190);

    auto g_0_zz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sdsk + 191);

    auto g_0_zz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sdsk + 192);

    auto g_0_zz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sdsk + 193);

    auto g_0_zz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 194);

    auto g_0_zz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sdsk + 195);

    auto g_0_zz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sdsk + 196);

    auto g_0_zz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sdsk + 197);

    auto g_0_zz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sdsk + 198);

    auto g_0_zz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 199);

    auto g_0_zz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 200);

    auto g_0_zz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sdsk + 201);

    auto g_0_zz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sdsk + 202);

    auto g_0_zz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sdsk + 203);

    auto g_0_zz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sdsk + 204);

    auto g_0_zz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 205);

    auto g_0_zz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 206);

    auto g_0_zz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 207);

    auto g_0_zz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sdsk + 208);

    auto g_0_zz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sdsk + 209);

    auto g_0_zz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sdsk + 210);

    auto g_0_zz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sdsk + 211);

    auto g_0_zz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 212);

    auto g_0_zz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 213);

    auto g_0_zz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 214);

    auto g_0_zz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sdsk + 215);

    #pragma omp simd aligned(g_0_0_0_xxxxxxx_0, g_0_0_0_xxxxxxx_1, g_0_0_0_xxxxxxy_0, g_0_0_0_xxxxxxy_1, g_0_0_0_xxxxxxz_0, g_0_0_0_xxxxxxz_1, g_0_0_0_xxxxxyy_0, g_0_0_0_xxxxxyy_1, g_0_0_0_xxxxxyz_0, g_0_0_0_xxxxxyz_1, g_0_0_0_xxxxxzz_0, g_0_0_0_xxxxxzz_1, g_0_0_0_xxxxyyy_0, g_0_0_0_xxxxyyy_1, g_0_0_0_xxxxyyz_0, g_0_0_0_xxxxyyz_1, g_0_0_0_xxxxyzz_0, g_0_0_0_xxxxyzz_1, g_0_0_0_xxxxzzz_0, g_0_0_0_xxxxzzz_1, g_0_0_0_xxxyyyy_0, g_0_0_0_xxxyyyy_1, g_0_0_0_xxxyyyz_0, g_0_0_0_xxxyyyz_1, g_0_0_0_xxxyyzz_0, g_0_0_0_xxxyyzz_1, g_0_0_0_xxxyzzz_0, g_0_0_0_xxxyzzz_1, g_0_0_0_xxxzzzz_0, g_0_0_0_xxxzzzz_1, g_0_0_0_xxyyyyy_0, g_0_0_0_xxyyyyy_1, g_0_0_0_xxyyyyz_0, g_0_0_0_xxyyyyz_1, g_0_0_0_xxyyyzz_0, g_0_0_0_xxyyyzz_1, g_0_0_0_xxyyzzz_0, g_0_0_0_xxyyzzz_1, g_0_0_0_xxyzzzz_0, g_0_0_0_xxyzzzz_1, g_0_0_0_xxzzzzz_0, g_0_0_0_xxzzzzz_1, g_0_0_0_xyyyyyy_0, g_0_0_0_xyyyyyy_1, g_0_0_0_xyyyyyz_0, g_0_0_0_xyyyyyz_1, g_0_0_0_xyyyyzz_0, g_0_0_0_xyyyyzz_1, g_0_0_0_xyyyzzz_0, g_0_0_0_xyyyzzz_1, g_0_0_0_xyyzzzz_0, g_0_0_0_xyyzzzz_1, g_0_0_0_xyzzzzz_0, g_0_0_0_xyzzzzz_1, g_0_0_0_xzzzzzz_0, g_0_0_0_xzzzzzz_1, g_0_0_0_yyyyyyy_0, g_0_0_0_yyyyyyy_1, g_0_0_0_yyyyyyz_0, g_0_0_0_yyyyyyz_1, g_0_0_0_yyyyyzz_0, g_0_0_0_yyyyyzz_1, g_0_0_0_yyyyzzz_0, g_0_0_0_yyyyzzz_1, g_0_0_0_yyyzzzz_0, g_0_0_0_yyyzzzz_1, g_0_0_0_yyzzzzz_0, g_0_0_0_yyzzzzz_1, g_0_0_0_yzzzzzz_0, g_0_0_0_yzzzzzz_1, g_0_0_0_zzzzzzz_0, g_0_0_0_zzzzzzz_1, g_0_z_0_xxxxxx_1, g_0_z_0_xxxxxxx_0, g_0_z_0_xxxxxxx_1, g_0_z_0_xxxxxxy_0, g_0_z_0_xxxxxxy_1, g_0_z_0_xxxxxxz_0, g_0_z_0_xxxxxxz_1, g_0_z_0_xxxxxy_1, g_0_z_0_xxxxxyy_0, g_0_z_0_xxxxxyy_1, g_0_z_0_xxxxxyz_0, g_0_z_0_xxxxxyz_1, g_0_z_0_xxxxxz_1, g_0_z_0_xxxxxzz_0, g_0_z_0_xxxxxzz_1, g_0_z_0_xxxxyy_1, g_0_z_0_xxxxyyy_0, g_0_z_0_xxxxyyy_1, g_0_z_0_xxxxyyz_0, g_0_z_0_xxxxyyz_1, g_0_z_0_xxxxyz_1, g_0_z_0_xxxxyzz_0, g_0_z_0_xxxxyzz_1, g_0_z_0_xxxxzz_1, g_0_z_0_xxxxzzz_0, g_0_z_0_xxxxzzz_1, g_0_z_0_xxxyyy_1, g_0_z_0_xxxyyyy_0, g_0_z_0_xxxyyyy_1, g_0_z_0_xxxyyyz_0, g_0_z_0_xxxyyyz_1, g_0_z_0_xxxyyz_1, g_0_z_0_xxxyyzz_0, g_0_z_0_xxxyyzz_1, g_0_z_0_xxxyzz_1, g_0_z_0_xxxyzzz_0, g_0_z_0_xxxyzzz_1, g_0_z_0_xxxzzz_1, g_0_z_0_xxxzzzz_0, g_0_z_0_xxxzzzz_1, g_0_z_0_xxyyyy_1, g_0_z_0_xxyyyyy_0, g_0_z_0_xxyyyyy_1, g_0_z_0_xxyyyyz_0, g_0_z_0_xxyyyyz_1, g_0_z_0_xxyyyz_1, g_0_z_0_xxyyyzz_0, g_0_z_0_xxyyyzz_1, g_0_z_0_xxyyzz_1, g_0_z_0_xxyyzzz_0, g_0_z_0_xxyyzzz_1, g_0_z_0_xxyzzz_1, g_0_z_0_xxyzzzz_0, g_0_z_0_xxyzzzz_1, g_0_z_0_xxzzzz_1, g_0_z_0_xxzzzzz_0, g_0_z_0_xxzzzzz_1, g_0_z_0_xyyyyy_1, g_0_z_0_xyyyyyy_0, g_0_z_0_xyyyyyy_1, g_0_z_0_xyyyyyz_0, g_0_z_0_xyyyyyz_1, g_0_z_0_xyyyyz_1, g_0_z_0_xyyyyzz_0, g_0_z_0_xyyyyzz_1, g_0_z_0_xyyyzz_1, g_0_z_0_xyyyzzz_0, g_0_z_0_xyyyzzz_1, g_0_z_0_xyyzzz_1, g_0_z_0_xyyzzzz_0, g_0_z_0_xyyzzzz_1, g_0_z_0_xyzzzz_1, g_0_z_0_xyzzzzz_0, g_0_z_0_xyzzzzz_1, g_0_z_0_xzzzzz_1, g_0_z_0_xzzzzzz_0, g_0_z_0_xzzzzzz_1, g_0_z_0_yyyyyy_1, g_0_z_0_yyyyyyy_0, g_0_z_0_yyyyyyy_1, g_0_z_0_yyyyyyz_0, g_0_z_0_yyyyyyz_1, g_0_z_0_yyyyyz_1, g_0_z_0_yyyyyzz_0, g_0_z_0_yyyyyzz_1, g_0_z_0_yyyyzz_1, g_0_z_0_yyyyzzz_0, g_0_z_0_yyyyzzz_1, g_0_z_0_yyyzzz_1, g_0_z_0_yyyzzzz_0, g_0_z_0_yyyzzzz_1, g_0_z_0_yyzzzz_1, g_0_z_0_yyzzzzz_0, g_0_z_0_yyzzzzz_1, g_0_z_0_yzzzzz_1, g_0_z_0_yzzzzzz_0, g_0_z_0_yzzzzzz_1, g_0_z_0_zzzzzz_1, g_0_z_0_zzzzzzz_0, g_0_z_0_zzzzzzz_1, g_0_zz_0_xxxxxxx_0, g_0_zz_0_xxxxxxy_0, g_0_zz_0_xxxxxxz_0, g_0_zz_0_xxxxxyy_0, g_0_zz_0_xxxxxyz_0, g_0_zz_0_xxxxxzz_0, g_0_zz_0_xxxxyyy_0, g_0_zz_0_xxxxyyz_0, g_0_zz_0_xxxxyzz_0, g_0_zz_0_xxxxzzz_0, g_0_zz_0_xxxyyyy_0, g_0_zz_0_xxxyyyz_0, g_0_zz_0_xxxyyzz_0, g_0_zz_0_xxxyzzz_0, g_0_zz_0_xxxzzzz_0, g_0_zz_0_xxyyyyy_0, g_0_zz_0_xxyyyyz_0, g_0_zz_0_xxyyyzz_0, g_0_zz_0_xxyyzzz_0, g_0_zz_0_xxyzzzz_0, g_0_zz_0_xxzzzzz_0, g_0_zz_0_xyyyyyy_0, g_0_zz_0_xyyyyyz_0, g_0_zz_0_xyyyyzz_0, g_0_zz_0_xyyyzzz_0, g_0_zz_0_xyyzzzz_0, g_0_zz_0_xyzzzzz_0, g_0_zz_0_xzzzzzz_0, g_0_zz_0_yyyyyyy_0, g_0_zz_0_yyyyyyz_0, g_0_zz_0_yyyyyzz_0, g_0_zz_0_yyyyzzz_0, g_0_zz_0_yyyzzzz_0, g_0_zz_0_yyzzzzz_0, g_0_zz_0_yzzzzzz_0, g_0_zz_0_zzzzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zz_0_xxxxxxx_0[i] = g_0_0_0_xxxxxxx_0[i] * fi_ab_0 - g_0_0_0_xxxxxxx_1[i] * fti_ab_0 + g_0_z_0_xxxxxxx_0[i] * pb_z + g_0_z_0_xxxxxxx_1[i] * wp_z[i];

        g_0_zz_0_xxxxxxy_0[i] = g_0_0_0_xxxxxxy_0[i] * fi_ab_0 - g_0_0_0_xxxxxxy_1[i] * fti_ab_0 + g_0_z_0_xxxxxxy_0[i] * pb_z + g_0_z_0_xxxxxxy_1[i] * wp_z[i];

        g_0_zz_0_xxxxxxz_0[i] = g_0_0_0_xxxxxxz_0[i] * fi_ab_0 - g_0_0_0_xxxxxxz_1[i] * fti_ab_0 + g_0_z_0_xxxxxx_1[i] * fi_abcd_0 + g_0_z_0_xxxxxxz_0[i] * pb_z + g_0_z_0_xxxxxxz_1[i] * wp_z[i];

        g_0_zz_0_xxxxxyy_0[i] = g_0_0_0_xxxxxyy_0[i] * fi_ab_0 - g_0_0_0_xxxxxyy_1[i] * fti_ab_0 + g_0_z_0_xxxxxyy_0[i] * pb_z + g_0_z_0_xxxxxyy_1[i] * wp_z[i];

        g_0_zz_0_xxxxxyz_0[i] = g_0_0_0_xxxxxyz_0[i] * fi_ab_0 - g_0_0_0_xxxxxyz_1[i] * fti_ab_0 + g_0_z_0_xxxxxy_1[i] * fi_abcd_0 + g_0_z_0_xxxxxyz_0[i] * pb_z + g_0_z_0_xxxxxyz_1[i] * wp_z[i];

        g_0_zz_0_xxxxxzz_0[i] = g_0_0_0_xxxxxzz_0[i] * fi_ab_0 - g_0_0_0_xxxxxzz_1[i] * fti_ab_0 + 2.0 * g_0_z_0_xxxxxz_1[i] * fi_abcd_0 + g_0_z_0_xxxxxzz_0[i] * pb_z + g_0_z_0_xxxxxzz_1[i] * wp_z[i];

        g_0_zz_0_xxxxyyy_0[i] = g_0_0_0_xxxxyyy_0[i] * fi_ab_0 - g_0_0_0_xxxxyyy_1[i] * fti_ab_0 + g_0_z_0_xxxxyyy_0[i] * pb_z + g_0_z_0_xxxxyyy_1[i] * wp_z[i];

        g_0_zz_0_xxxxyyz_0[i] = g_0_0_0_xxxxyyz_0[i] * fi_ab_0 - g_0_0_0_xxxxyyz_1[i] * fti_ab_0 + g_0_z_0_xxxxyy_1[i] * fi_abcd_0 + g_0_z_0_xxxxyyz_0[i] * pb_z + g_0_z_0_xxxxyyz_1[i] * wp_z[i];

        g_0_zz_0_xxxxyzz_0[i] = g_0_0_0_xxxxyzz_0[i] * fi_ab_0 - g_0_0_0_xxxxyzz_1[i] * fti_ab_0 + 2.0 * g_0_z_0_xxxxyz_1[i] * fi_abcd_0 + g_0_z_0_xxxxyzz_0[i] * pb_z + g_0_z_0_xxxxyzz_1[i] * wp_z[i];

        g_0_zz_0_xxxxzzz_0[i] = g_0_0_0_xxxxzzz_0[i] * fi_ab_0 - g_0_0_0_xxxxzzz_1[i] * fti_ab_0 + 3.0 * g_0_z_0_xxxxzz_1[i] * fi_abcd_0 + g_0_z_0_xxxxzzz_0[i] * pb_z + g_0_z_0_xxxxzzz_1[i] * wp_z[i];

        g_0_zz_0_xxxyyyy_0[i] = g_0_0_0_xxxyyyy_0[i] * fi_ab_0 - g_0_0_0_xxxyyyy_1[i] * fti_ab_0 + g_0_z_0_xxxyyyy_0[i] * pb_z + g_0_z_0_xxxyyyy_1[i] * wp_z[i];

        g_0_zz_0_xxxyyyz_0[i] = g_0_0_0_xxxyyyz_0[i] * fi_ab_0 - g_0_0_0_xxxyyyz_1[i] * fti_ab_0 + g_0_z_0_xxxyyy_1[i] * fi_abcd_0 + g_0_z_0_xxxyyyz_0[i] * pb_z + g_0_z_0_xxxyyyz_1[i] * wp_z[i];

        g_0_zz_0_xxxyyzz_0[i] = g_0_0_0_xxxyyzz_0[i] * fi_ab_0 - g_0_0_0_xxxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_z_0_xxxyyz_1[i] * fi_abcd_0 + g_0_z_0_xxxyyzz_0[i] * pb_z + g_0_z_0_xxxyyzz_1[i] * wp_z[i];

        g_0_zz_0_xxxyzzz_0[i] = g_0_0_0_xxxyzzz_0[i] * fi_ab_0 - g_0_0_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_z_0_xxxyzz_1[i] * fi_abcd_0 + g_0_z_0_xxxyzzz_0[i] * pb_z + g_0_z_0_xxxyzzz_1[i] * wp_z[i];

        g_0_zz_0_xxxzzzz_0[i] = g_0_0_0_xxxzzzz_0[i] * fi_ab_0 - g_0_0_0_xxxzzzz_1[i] * fti_ab_0 + 4.0 * g_0_z_0_xxxzzz_1[i] * fi_abcd_0 + g_0_z_0_xxxzzzz_0[i] * pb_z + g_0_z_0_xxxzzzz_1[i] * wp_z[i];

        g_0_zz_0_xxyyyyy_0[i] = g_0_0_0_xxyyyyy_0[i] * fi_ab_0 - g_0_0_0_xxyyyyy_1[i] * fti_ab_0 + g_0_z_0_xxyyyyy_0[i] * pb_z + g_0_z_0_xxyyyyy_1[i] * wp_z[i];

        g_0_zz_0_xxyyyyz_0[i] = g_0_0_0_xxyyyyz_0[i] * fi_ab_0 - g_0_0_0_xxyyyyz_1[i] * fti_ab_0 + g_0_z_0_xxyyyy_1[i] * fi_abcd_0 + g_0_z_0_xxyyyyz_0[i] * pb_z + g_0_z_0_xxyyyyz_1[i] * wp_z[i];

        g_0_zz_0_xxyyyzz_0[i] = g_0_0_0_xxyyyzz_0[i] * fi_ab_0 - g_0_0_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_z_0_xxyyyz_1[i] * fi_abcd_0 + g_0_z_0_xxyyyzz_0[i] * pb_z + g_0_z_0_xxyyyzz_1[i] * wp_z[i];

        g_0_zz_0_xxyyzzz_0[i] = g_0_0_0_xxyyzzz_0[i] * fi_ab_0 - g_0_0_0_xxyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_z_0_xxyyzz_1[i] * fi_abcd_0 + g_0_z_0_xxyyzzz_0[i] * pb_z + g_0_z_0_xxyyzzz_1[i] * wp_z[i];

        g_0_zz_0_xxyzzzz_0[i] = g_0_0_0_xxyzzzz_0[i] * fi_ab_0 - g_0_0_0_xxyzzzz_1[i] * fti_ab_0 + 4.0 * g_0_z_0_xxyzzz_1[i] * fi_abcd_0 + g_0_z_0_xxyzzzz_0[i] * pb_z + g_0_z_0_xxyzzzz_1[i] * wp_z[i];

        g_0_zz_0_xxzzzzz_0[i] = g_0_0_0_xxzzzzz_0[i] * fi_ab_0 - g_0_0_0_xxzzzzz_1[i] * fti_ab_0 + 5.0 * g_0_z_0_xxzzzz_1[i] * fi_abcd_0 + g_0_z_0_xxzzzzz_0[i] * pb_z + g_0_z_0_xxzzzzz_1[i] * wp_z[i];

        g_0_zz_0_xyyyyyy_0[i] = g_0_0_0_xyyyyyy_0[i] * fi_ab_0 - g_0_0_0_xyyyyyy_1[i] * fti_ab_0 + g_0_z_0_xyyyyyy_0[i] * pb_z + g_0_z_0_xyyyyyy_1[i] * wp_z[i];

        g_0_zz_0_xyyyyyz_0[i] = g_0_0_0_xyyyyyz_0[i] * fi_ab_0 - g_0_0_0_xyyyyyz_1[i] * fti_ab_0 + g_0_z_0_xyyyyy_1[i] * fi_abcd_0 + g_0_z_0_xyyyyyz_0[i] * pb_z + g_0_z_0_xyyyyyz_1[i] * wp_z[i];

        g_0_zz_0_xyyyyzz_0[i] = g_0_0_0_xyyyyzz_0[i] * fi_ab_0 - g_0_0_0_xyyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_z_0_xyyyyz_1[i] * fi_abcd_0 + g_0_z_0_xyyyyzz_0[i] * pb_z + g_0_z_0_xyyyyzz_1[i] * wp_z[i];

        g_0_zz_0_xyyyzzz_0[i] = g_0_0_0_xyyyzzz_0[i] * fi_ab_0 - g_0_0_0_xyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_z_0_xyyyzz_1[i] * fi_abcd_0 + g_0_z_0_xyyyzzz_0[i] * pb_z + g_0_z_0_xyyyzzz_1[i] * wp_z[i];

        g_0_zz_0_xyyzzzz_0[i] = g_0_0_0_xyyzzzz_0[i] * fi_ab_0 - g_0_0_0_xyyzzzz_1[i] * fti_ab_0 + 4.0 * g_0_z_0_xyyzzz_1[i] * fi_abcd_0 + g_0_z_0_xyyzzzz_0[i] * pb_z + g_0_z_0_xyyzzzz_1[i] * wp_z[i];

        g_0_zz_0_xyzzzzz_0[i] = g_0_0_0_xyzzzzz_0[i] * fi_ab_0 - g_0_0_0_xyzzzzz_1[i] * fti_ab_0 + 5.0 * g_0_z_0_xyzzzz_1[i] * fi_abcd_0 + g_0_z_0_xyzzzzz_0[i] * pb_z + g_0_z_0_xyzzzzz_1[i] * wp_z[i];

        g_0_zz_0_xzzzzzz_0[i] = g_0_0_0_xzzzzzz_0[i] * fi_ab_0 - g_0_0_0_xzzzzzz_1[i] * fti_ab_0 + 6.0 * g_0_z_0_xzzzzz_1[i] * fi_abcd_0 + g_0_z_0_xzzzzzz_0[i] * pb_z + g_0_z_0_xzzzzzz_1[i] * wp_z[i];

        g_0_zz_0_yyyyyyy_0[i] = g_0_0_0_yyyyyyy_0[i] * fi_ab_0 - g_0_0_0_yyyyyyy_1[i] * fti_ab_0 + g_0_z_0_yyyyyyy_0[i] * pb_z + g_0_z_0_yyyyyyy_1[i] * wp_z[i];

        g_0_zz_0_yyyyyyz_0[i] = g_0_0_0_yyyyyyz_0[i] * fi_ab_0 - g_0_0_0_yyyyyyz_1[i] * fti_ab_0 + g_0_z_0_yyyyyy_1[i] * fi_abcd_0 + g_0_z_0_yyyyyyz_0[i] * pb_z + g_0_z_0_yyyyyyz_1[i] * wp_z[i];

        g_0_zz_0_yyyyyzz_0[i] = g_0_0_0_yyyyyzz_0[i] * fi_ab_0 - g_0_0_0_yyyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_z_0_yyyyyz_1[i] * fi_abcd_0 + g_0_z_0_yyyyyzz_0[i] * pb_z + g_0_z_0_yyyyyzz_1[i] * wp_z[i];

        g_0_zz_0_yyyyzzz_0[i] = g_0_0_0_yyyyzzz_0[i] * fi_ab_0 - g_0_0_0_yyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_z_0_yyyyzz_1[i] * fi_abcd_0 + g_0_z_0_yyyyzzz_0[i] * pb_z + g_0_z_0_yyyyzzz_1[i] * wp_z[i];

        g_0_zz_0_yyyzzzz_0[i] = g_0_0_0_yyyzzzz_0[i] * fi_ab_0 - g_0_0_0_yyyzzzz_1[i] * fti_ab_0 + 4.0 * g_0_z_0_yyyzzz_1[i] * fi_abcd_0 + g_0_z_0_yyyzzzz_0[i] * pb_z + g_0_z_0_yyyzzzz_1[i] * wp_z[i];

        g_0_zz_0_yyzzzzz_0[i] = g_0_0_0_yyzzzzz_0[i] * fi_ab_0 - g_0_0_0_yyzzzzz_1[i] * fti_ab_0 + 5.0 * g_0_z_0_yyzzzz_1[i] * fi_abcd_0 + g_0_z_0_yyzzzzz_0[i] * pb_z + g_0_z_0_yyzzzzz_1[i] * wp_z[i];

        g_0_zz_0_yzzzzzz_0[i] = g_0_0_0_yzzzzzz_0[i] * fi_ab_0 - g_0_0_0_yzzzzzz_1[i] * fti_ab_0 + 6.0 * g_0_z_0_yzzzzz_1[i] * fi_abcd_0 + g_0_z_0_yzzzzzz_0[i] * pb_z + g_0_z_0_yzzzzzz_1[i] * wp_z[i];

        g_0_zz_0_zzzzzzz_0[i] = g_0_0_0_zzzzzzz_0[i] * fi_ab_0 - g_0_0_0_zzzzzzz_1[i] * fti_ab_0 + 7.0 * g_0_z_0_zzzzzz_1[i] * fi_abcd_0 + g_0_z_0_zzzzzzz_0[i] * pb_z + g_0_z_0_zzzzzzz_1[i] * wp_z[i];
    }
}

} // erirec namespace


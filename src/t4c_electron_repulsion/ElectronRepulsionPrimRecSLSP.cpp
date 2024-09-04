#include "ElectronRepulsionPrimRecSLSP.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_slsp(CSimdArray<double>& pbuffer,
                                  const size_t idx_eri_0_slsp,
                                  size_t idx_eri_0_sisp,
                                  size_t idx_eri_1_sisp,
                                  size_t idx_eri_1_skss,
                                  size_t idx_eri_0_sksp,
                                  size_t idx_eri_1_sksp,
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

    /// Set up components of auxilary buffer : SISP

    auto g_0_xxxxxx_0_x_0 = pbuffer.data(idx_eri_0_sisp);

    auto g_0_xxxxxx_0_y_0 = pbuffer.data(idx_eri_0_sisp + 1);

    auto g_0_xxxxxx_0_z_0 = pbuffer.data(idx_eri_0_sisp + 2);

    auto g_0_xxxxxy_0_x_0 = pbuffer.data(idx_eri_0_sisp + 3);

    auto g_0_xxxxxz_0_x_0 = pbuffer.data(idx_eri_0_sisp + 6);

    auto g_0_xxxxyy_0_x_0 = pbuffer.data(idx_eri_0_sisp + 9);

    auto g_0_xxxxyy_0_y_0 = pbuffer.data(idx_eri_0_sisp + 10);

    auto g_0_xxxxyy_0_z_0 = pbuffer.data(idx_eri_0_sisp + 11);

    auto g_0_xxxxzz_0_x_0 = pbuffer.data(idx_eri_0_sisp + 15);

    auto g_0_xxxxzz_0_y_0 = pbuffer.data(idx_eri_0_sisp + 16);

    auto g_0_xxxxzz_0_z_0 = pbuffer.data(idx_eri_0_sisp + 17);

    auto g_0_xxxyyy_0_x_0 = pbuffer.data(idx_eri_0_sisp + 18);

    auto g_0_xxxyyy_0_y_0 = pbuffer.data(idx_eri_0_sisp + 19);

    auto g_0_xxxyyy_0_z_0 = pbuffer.data(idx_eri_0_sisp + 20);

    auto g_0_xxxyzz_0_x_0 = pbuffer.data(idx_eri_0_sisp + 24);

    auto g_0_xxxzzz_0_x_0 = pbuffer.data(idx_eri_0_sisp + 27);

    auto g_0_xxxzzz_0_y_0 = pbuffer.data(idx_eri_0_sisp + 28);

    auto g_0_xxxzzz_0_z_0 = pbuffer.data(idx_eri_0_sisp + 29);

    auto g_0_xxyyyy_0_x_0 = pbuffer.data(idx_eri_0_sisp + 30);

    auto g_0_xxyyyy_0_y_0 = pbuffer.data(idx_eri_0_sisp + 31);

    auto g_0_xxyyyy_0_z_0 = pbuffer.data(idx_eri_0_sisp + 32);

    auto g_0_xxyyzz_0_x_0 = pbuffer.data(idx_eri_0_sisp + 36);

    auto g_0_xxyyzz_0_y_0 = pbuffer.data(idx_eri_0_sisp + 37);

    auto g_0_xxyyzz_0_z_0 = pbuffer.data(idx_eri_0_sisp + 38);

    auto g_0_xxyzzz_0_x_0 = pbuffer.data(idx_eri_0_sisp + 39);

    auto g_0_xxzzzz_0_x_0 = pbuffer.data(idx_eri_0_sisp + 42);

    auto g_0_xxzzzz_0_y_0 = pbuffer.data(idx_eri_0_sisp + 43);

    auto g_0_xxzzzz_0_z_0 = pbuffer.data(idx_eri_0_sisp + 44);

    auto g_0_xyyyyy_0_y_0 = pbuffer.data(idx_eri_0_sisp + 46);

    auto g_0_xyyyyy_0_z_0 = pbuffer.data(idx_eri_0_sisp + 47);

    auto g_0_xyyyzz_0_y_0 = pbuffer.data(idx_eri_0_sisp + 52);

    auto g_0_xyyyzz_0_z_0 = pbuffer.data(idx_eri_0_sisp + 53);

    auto g_0_xyyzzz_0_y_0 = pbuffer.data(idx_eri_0_sisp + 55);

    auto g_0_xyyzzz_0_z_0 = pbuffer.data(idx_eri_0_sisp + 56);

    auto g_0_xzzzzz_0_y_0 = pbuffer.data(idx_eri_0_sisp + 61);

    auto g_0_xzzzzz_0_z_0 = pbuffer.data(idx_eri_0_sisp + 62);

    auto g_0_yyyyyy_0_x_0 = pbuffer.data(idx_eri_0_sisp + 63);

    auto g_0_yyyyyy_0_y_0 = pbuffer.data(idx_eri_0_sisp + 64);

    auto g_0_yyyyyy_0_z_0 = pbuffer.data(idx_eri_0_sisp + 65);

    auto g_0_yyyyyz_0_y_0 = pbuffer.data(idx_eri_0_sisp + 67);

    auto g_0_yyyyzz_0_x_0 = pbuffer.data(idx_eri_0_sisp + 69);

    auto g_0_yyyyzz_0_y_0 = pbuffer.data(idx_eri_0_sisp + 70);

    auto g_0_yyyyzz_0_z_0 = pbuffer.data(idx_eri_0_sisp + 71);

    auto g_0_yyyzzz_0_x_0 = pbuffer.data(idx_eri_0_sisp + 72);

    auto g_0_yyyzzz_0_y_0 = pbuffer.data(idx_eri_0_sisp + 73);

    auto g_0_yyyzzz_0_z_0 = pbuffer.data(idx_eri_0_sisp + 74);

    auto g_0_yyzzzz_0_x_0 = pbuffer.data(idx_eri_0_sisp + 75);

    auto g_0_yyzzzz_0_y_0 = pbuffer.data(idx_eri_0_sisp + 76);

    auto g_0_yyzzzz_0_z_0 = pbuffer.data(idx_eri_0_sisp + 77);

    auto g_0_yzzzzz_0_x_0 = pbuffer.data(idx_eri_0_sisp + 78);

    auto g_0_yzzzzz_0_z_0 = pbuffer.data(idx_eri_0_sisp + 80);

    auto g_0_zzzzzz_0_x_0 = pbuffer.data(idx_eri_0_sisp + 81);

    auto g_0_zzzzzz_0_y_0 = pbuffer.data(idx_eri_0_sisp + 82);

    auto g_0_zzzzzz_0_z_0 = pbuffer.data(idx_eri_0_sisp + 83);

    /// Set up components of auxilary buffer : SISP

    auto g_0_xxxxxx_0_x_1 = pbuffer.data(idx_eri_1_sisp);

    auto g_0_xxxxxx_0_y_1 = pbuffer.data(idx_eri_1_sisp + 1);

    auto g_0_xxxxxx_0_z_1 = pbuffer.data(idx_eri_1_sisp + 2);

    auto g_0_xxxxxy_0_x_1 = pbuffer.data(idx_eri_1_sisp + 3);

    auto g_0_xxxxxz_0_x_1 = pbuffer.data(idx_eri_1_sisp + 6);

    auto g_0_xxxxyy_0_x_1 = pbuffer.data(idx_eri_1_sisp + 9);

    auto g_0_xxxxyy_0_y_1 = pbuffer.data(idx_eri_1_sisp + 10);

    auto g_0_xxxxyy_0_z_1 = pbuffer.data(idx_eri_1_sisp + 11);

    auto g_0_xxxxzz_0_x_1 = pbuffer.data(idx_eri_1_sisp + 15);

    auto g_0_xxxxzz_0_y_1 = pbuffer.data(idx_eri_1_sisp + 16);

    auto g_0_xxxxzz_0_z_1 = pbuffer.data(idx_eri_1_sisp + 17);

    auto g_0_xxxyyy_0_x_1 = pbuffer.data(idx_eri_1_sisp + 18);

    auto g_0_xxxyyy_0_y_1 = pbuffer.data(idx_eri_1_sisp + 19);

    auto g_0_xxxyyy_0_z_1 = pbuffer.data(idx_eri_1_sisp + 20);

    auto g_0_xxxyzz_0_x_1 = pbuffer.data(idx_eri_1_sisp + 24);

    auto g_0_xxxzzz_0_x_1 = pbuffer.data(idx_eri_1_sisp + 27);

    auto g_0_xxxzzz_0_y_1 = pbuffer.data(idx_eri_1_sisp + 28);

    auto g_0_xxxzzz_0_z_1 = pbuffer.data(idx_eri_1_sisp + 29);

    auto g_0_xxyyyy_0_x_1 = pbuffer.data(idx_eri_1_sisp + 30);

    auto g_0_xxyyyy_0_y_1 = pbuffer.data(idx_eri_1_sisp + 31);

    auto g_0_xxyyyy_0_z_1 = pbuffer.data(idx_eri_1_sisp + 32);

    auto g_0_xxyyzz_0_x_1 = pbuffer.data(idx_eri_1_sisp + 36);

    auto g_0_xxyyzz_0_y_1 = pbuffer.data(idx_eri_1_sisp + 37);

    auto g_0_xxyyzz_0_z_1 = pbuffer.data(idx_eri_1_sisp + 38);

    auto g_0_xxyzzz_0_x_1 = pbuffer.data(idx_eri_1_sisp + 39);

    auto g_0_xxzzzz_0_x_1 = pbuffer.data(idx_eri_1_sisp + 42);

    auto g_0_xxzzzz_0_y_1 = pbuffer.data(idx_eri_1_sisp + 43);

    auto g_0_xxzzzz_0_z_1 = pbuffer.data(idx_eri_1_sisp + 44);

    auto g_0_xyyyyy_0_y_1 = pbuffer.data(idx_eri_1_sisp + 46);

    auto g_0_xyyyyy_0_z_1 = pbuffer.data(idx_eri_1_sisp + 47);

    auto g_0_xyyyzz_0_y_1 = pbuffer.data(idx_eri_1_sisp + 52);

    auto g_0_xyyyzz_0_z_1 = pbuffer.data(idx_eri_1_sisp + 53);

    auto g_0_xyyzzz_0_y_1 = pbuffer.data(idx_eri_1_sisp + 55);

    auto g_0_xyyzzz_0_z_1 = pbuffer.data(idx_eri_1_sisp + 56);

    auto g_0_xzzzzz_0_y_1 = pbuffer.data(idx_eri_1_sisp + 61);

    auto g_0_xzzzzz_0_z_1 = pbuffer.data(idx_eri_1_sisp + 62);

    auto g_0_yyyyyy_0_x_1 = pbuffer.data(idx_eri_1_sisp + 63);

    auto g_0_yyyyyy_0_y_1 = pbuffer.data(idx_eri_1_sisp + 64);

    auto g_0_yyyyyy_0_z_1 = pbuffer.data(idx_eri_1_sisp + 65);

    auto g_0_yyyyyz_0_y_1 = pbuffer.data(idx_eri_1_sisp + 67);

    auto g_0_yyyyzz_0_x_1 = pbuffer.data(idx_eri_1_sisp + 69);

    auto g_0_yyyyzz_0_y_1 = pbuffer.data(idx_eri_1_sisp + 70);

    auto g_0_yyyyzz_0_z_1 = pbuffer.data(idx_eri_1_sisp + 71);

    auto g_0_yyyzzz_0_x_1 = pbuffer.data(idx_eri_1_sisp + 72);

    auto g_0_yyyzzz_0_y_1 = pbuffer.data(idx_eri_1_sisp + 73);

    auto g_0_yyyzzz_0_z_1 = pbuffer.data(idx_eri_1_sisp + 74);

    auto g_0_yyzzzz_0_x_1 = pbuffer.data(idx_eri_1_sisp + 75);

    auto g_0_yyzzzz_0_y_1 = pbuffer.data(idx_eri_1_sisp + 76);

    auto g_0_yyzzzz_0_z_1 = pbuffer.data(idx_eri_1_sisp + 77);

    auto g_0_yzzzzz_0_x_1 = pbuffer.data(idx_eri_1_sisp + 78);

    auto g_0_yzzzzz_0_z_1 = pbuffer.data(idx_eri_1_sisp + 80);

    auto g_0_zzzzzz_0_x_1 = pbuffer.data(idx_eri_1_sisp + 81);

    auto g_0_zzzzzz_0_y_1 = pbuffer.data(idx_eri_1_sisp + 82);

    auto g_0_zzzzzz_0_z_1 = pbuffer.data(idx_eri_1_sisp + 83);

    /// Set up components of auxilary buffer : SKSS

    auto g_0_xxxxxxx_0_0_1 = pbuffer.data(idx_eri_1_skss);

    auto g_0_xxxxxyy_0_0_1 = pbuffer.data(idx_eri_1_skss + 3);

    auto g_0_xxxxxzz_0_0_1 = pbuffer.data(idx_eri_1_skss + 5);

    auto g_0_xxxxyyy_0_0_1 = pbuffer.data(idx_eri_1_skss + 6);

    auto g_0_xxxxzzz_0_0_1 = pbuffer.data(idx_eri_1_skss + 9);

    auto g_0_xxxyyyy_0_0_1 = pbuffer.data(idx_eri_1_skss + 10);

    auto g_0_xxxzzzz_0_0_1 = pbuffer.data(idx_eri_1_skss + 14);

    auto g_0_xxyyyyy_0_0_1 = pbuffer.data(idx_eri_1_skss + 15);

    auto g_0_xxzzzzz_0_0_1 = pbuffer.data(idx_eri_1_skss + 20);

    auto g_0_yyyyyyy_0_0_1 = pbuffer.data(idx_eri_1_skss + 28);

    auto g_0_yyyyyzz_0_0_1 = pbuffer.data(idx_eri_1_skss + 30);

    auto g_0_yyyyzzz_0_0_1 = pbuffer.data(idx_eri_1_skss + 31);

    auto g_0_yyyzzzz_0_0_1 = pbuffer.data(idx_eri_1_skss + 32);

    auto g_0_yyzzzzz_0_0_1 = pbuffer.data(idx_eri_1_skss + 33);

    auto g_0_zzzzzzz_0_0_1 = pbuffer.data(idx_eri_1_skss + 35);

    /// Set up components of auxilary buffer : SKSP

    auto g_0_xxxxxxx_0_x_0 = pbuffer.data(idx_eri_0_sksp);

    auto g_0_xxxxxxx_0_y_0 = pbuffer.data(idx_eri_0_sksp + 1);

    auto g_0_xxxxxxx_0_z_0 = pbuffer.data(idx_eri_0_sksp + 2);

    auto g_0_xxxxxxy_0_x_0 = pbuffer.data(idx_eri_0_sksp + 3);

    auto g_0_xxxxxxy_0_y_0 = pbuffer.data(idx_eri_0_sksp + 4);

    auto g_0_xxxxxxz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 6);

    auto g_0_xxxxxxz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 8);

    auto g_0_xxxxxyy_0_x_0 = pbuffer.data(idx_eri_0_sksp + 9);

    auto g_0_xxxxxyy_0_y_0 = pbuffer.data(idx_eri_0_sksp + 10);

    auto g_0_xxxxxyy_0_z_0 = pbuffer.data(idx_eri_0_sksp + 11);

    auto g_0_xxxxxzz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 15);

    auto g_0_xxxxxzz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 16);

    auto g_0_xxxxxzz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 17);

    auto g_0_xxxxyyy_0_x_0 = pbuffer.data(idx_eri_0_sksp + 18);

    auto g_0_xxxxyyy_0_y_0 = pbuffer.data(idx_eri_0_sksp + 19);

    auto g_0_xxxxyyy_0_z_0 = pbuffer.data(idx_eri_0_sksp + 20);

    auto g_0_xxxxyzz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 24);

    auto g_0_xxxxzzz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 27);

    auto g_0_xxxxzzz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 28);

    auto g_0_xxxxzzz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 29);

    auto g_0_xxxyyyy_0_x_0 = pbuffer.data(idx_eri_0_sksp + 30);

    auto g_0_xxxyyyy_0_y_0 = pbuffer.data(idx_eri_0_sksp + 31);

    auto g_0_xxxyyyy_0_z_0 = pbuffer.data(idx_eri_0_sksp + 32);

    auto g_0_xxxyyzz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 36);

    auto g_0_xxxyyzz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 37);

    auto g_0_xxxyyzz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 38);

    auto g_0_xxxyzzz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 39);

    auto g_0_xxxzzzz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 42);

    auto g_0_xxxzzzz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 43);

    auto g_0_xxxzzzz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 44);

    auto g_0_xxyyyyy_0_x_0 = pbuffer.data(idx_eri_0_sksp + 45);

    auto g_0_xxyyyyy_0_y_0 = pbuffer.data(idx_eri_0_sksp + 46);

    auto g_0_xxyyyyy_0_z_0 = pbuffer.data(idx_eri_0_sksp + 47);

    auto g_0_xxyyyzz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 51);

    auto g_0_xxyyyzz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 52);

    auto g_0_xxyyyzz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 53);

    auto g_0_xxyyzzz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 54);

    auto g_0_xxyyzzz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 55);

    auto g_0_xxyyzzz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 56);

    auto g_0_xxyzzzz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 57);

    auto g_0_xxzzzzz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 60);

    auto g_0_xxzzzzz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 61);

    auto g_0_xxzzzzz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 62);

    auto g_0_xyyyyyy_0_x_0 = pbuffer.data(idx_eri_0_sksp + 63);

    auto g_0_xyyyyyy_0_y_0 = pbuffer.data(idx_eri_0_sksp + 64);

    auto g_0_xyyyyyy_0_z_0 = pbuffer.data(idx_eri_0_sksp + 65);

    auto g_0_xyyyyzz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 70);

    auto g_0_xyyyyzz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 71);

    auto g_0_xyyyzzz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 73);

    auto g_0_xyyyzzz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 74);

    auto g_0_xyyzzzz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 76);

    auto g_0_xyyzzzz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 77);

    auto g_0_xzzzzzz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 81);

    auto g_0_xzzzzzz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 82);

    auto g_0_xzzzzzz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 83);

    auto g_0_yyyyyyy_0_x_0 = pbuffer.data(idx_eri_0_sksp + 84);

    auto g_0_yyyyyyy_0_y_0 = pbuffer.data(idx_eri_0_sksp + 85);

    auto g_0_yyyyyyy_0_z_0 = pbuffer.data(idx_eri_0_sksp + 86);

    auto g_0_yyyyyyz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 88);

    auto g_0_yyyyyyz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 89);

    auto g_0_yyyyyzz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 90);

    auto g_0_yyyyyzz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 91);

    auto g_0_yyyyyzz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 92);

    auto g_0_yyyyzzz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 93);

    auto g_0_yyyyzzz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 94);

    auto g_0_yyyyzzz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 95);

    auto g_0_yyyzzzz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 96);

    auto g_0_yyyzzzz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 97);

    auto g_0_yyyzzzz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 98);

    auto g_0_yyzzzzz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 99);

    auto g_0_yyzzzzz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 100);

    auto g_0_yyzzzzz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 101);

    auto g_0_yzzzzzz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 102);

    auto g_0_yzzzzzz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 103);

    auto g_0_yzzzzzz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 104);

    auto g_0_zzzzzzz_0_x_0 = pbuffer.data(idx_eri_0_sksp + 105);

    auto g_0_zzzzzzz_0_y_0 = pbuffer.data(idx_eri_0_sksp + 106);

    auto g_0_zzzzzzz_0_z_0 = pbuffer.data(idx_eri_0_sksp + 107);

    /// Set up components of auxilary buffer : SKSP

    auto g_0_xxxxxxx_0_x_1 = pbuffer.data(idx_eri_1_sksp);

    auto g_0_xxxxxxx_0_y_1 = pbuffer.data(idx_eri_1_sksp + 1);

    auto g_0_xxxxxxx_0_z_1 = pbuffer.data(idx_eri_1_sksp + 2);

    auto g_0_xxxxxxy_0_x_1 = pbuffer.data(idx_eri_1_sksp + 3);

    auto g_0_xxxxxxy_0_y_1 = pbuffer.data(idx_eri_1_sksp + 4);

    auto g_0_xxxxxxz_0_x_1 = pbuffer.data(idx_eri_1_sksp + 6);

    auto g_0_xxxxxxz_0_z_1 = pbuffer.data(idx_eri_1_sksp + 8);

    auto g_0_xxxxxyy_0_x_1 = pbuffer.data(idx_eri_1_sksp + 9);

    auto g_0_xxxxxyy_0_y_1 = pbuffer.data(idx_eri_1_sksp + 10);

    auto g_0_xxxxxyy_0_z_1 = pbuffer.data(idx_eri_1_sksp + 11);

    auto g_0_xxxxxzz_0_x_1 = pbuffer.data(idx_eri_1_sksp + 15);

    auto g_0_xxxxxzz_0_y_1 = pbuffer.data(idx_eri_1_sksp + 16);

    auto g_0_xxxxxzz_0_z_1 = pbuffer.data(idx_eri_1_sksp + 17);

    auto g_0_xxxxyyy_0_x_1 = pbuffer.data(idx_eri_1_sksp + 18);

    auto g_0_xxxxyyy_0_y_1 = pbuffer.data(idx_eri_1_sksp + 19);

    auto g_0_xxxxyyy_0_z_1 = pbuffer.data(idx_eri_1_sksp + 20);

    auto g_0_xxxxyzz_0_x_1 = pbuffer.data(idx_eri_1_sksp + 24);

    auto g_0_xxxxzzz_0_x_1 = pbuffer.data(idx_eri_1_sksp + 27);

    auto g_0_xxxxzzz_0_y_1 = pbuffer.data(idx_eri_1_sksp + 28);

    auto g_0_xxxxzzz_0_z_1 = pbuffer.data(idx_eri_1_sksp + 29);

    auto g_0_xxxyyyy_0_x_1 = pbuffer.data(idx_eri_1_sksp + 30);

    auto g_0_xxxyyyy_0_y_1 = pbuffer.data(idx_eri_1_sksp + 31);

    auto g_0_xxxyyyy_0_z_1 = pbuffer.data(idx_eri_1_sksp + 32);

    auto g_0_xxxyyzz_0_x_1 = pbuffer.data(idx_eri_1_sksp + 36);

    auto g_0_xxxyyzz_0_y_1 = pbuffer.data(idx_eri_1_sksp + 37);

    auto g_0_xxxyyzz_0_z_1 = pbuffer.data(idx_eri_1_sksp + 38);

    auto g_0_xxxyzzz_0_x_1 = pbuffer.data(idx_eri_1_sksp + 39);

    auto g_0_xxxzzzz_0_x_1 = pbuffer.data(idx_eri_1_sksp + 42);

    auto g_0_xxxzzzz_0_y_1 = pbuffer.data(idx_eri_1_sksp + 43);

    auto g_0_xxxzzzz_0_z_1 = pbuffer.data(idx_eri_1_sksp + 44);

    auto g_0_xxyyyyy_0_x_1 = pbuffer.data(idx_eri_1_sksp + 45);

    auto g_0_xxyyyyy_0_y_1 = pbuffer.data(idx_eri_1_sksp + 46);

    auto g_0_xxyyyyy_0_z_1 = pbuffer.data(idx_eri_1_sksp + 47);

    auto g_0_xxyyyzz_0_x_1 = pbuffer.data(idx_eri_1_sksp + 51);

    auto g_0_xxyyyzz_0_y_1 = pbuffer.data(idx_eri_1_sksp + 52);

    auto g_0_xxyyyzz_0_z_1 = pbuffer.data(idx_eri_1_sksp + 53);

    auto g_0_xxyyzzz_0_x_1 = pbuffer.data(idx_eri_1_sksp + 54);

    auto g_0_xxyyzzz_0_y_1 = pbuffer.data(idx_eri_1_sksp + 55);

    auto g_0_xxyyzzz_0_z_1 = pbuffer.data(idx_eri_1_sksp + 56);

    auto g_0_xxyzzzz_0_x_1 = pbuffer.data(idx_eri_1_sksp + 57);

    auto g_0_xxzzzzz_0_x_1 = pbuffer.data(idx_eri_1_sksp + 60);

    auto g_0_xxzzzzz_0_y_1 = pbuffer.data(idx_eri_1_sksp + 61);

    auto g_0_xxzzzzz_0_z_1 = pbuffer.data(idx_eri_1_sksp + 62);

    auto g_0_xyyyyyy_0_x_1 = pbuffer.data(idx_eri_1_sksp + 63);

    auto g_0_xyyyyyy_0_y_1 = pbuffer.data(idx_eri_1_sksp + 64);

    auto g_0_xyyyyyy_0_z_1 = pbuffer.data(idx_eri_1_sksp + 65);

    auto g_0_xyyyyzz_0_y_1 = pbuffer.data(idx_eri_1_sksp + 70);

    auto g_0_xyyyyzz_0_z_1 = pbuffer.data(idx_eri_1_sksp + 71);

    auto g_0_xyyyzzz_0_y_1 = pbuffer.data(idx_eri_1_sksp + 73);

    auto g_0_xyyyzzz_0_z_1 = pbuffer.data(idx_eri_1_sksp + 74);

    auto g_0_xyyzzzz_0_y_1 = pbuffer.data(idx_eri_1_sksp + 76);

    auto g_0_xyyzzzz_0_z_1 = pbuffer.data(idx_eri_1_sksp + 77);

    auto g_0_xzzzzzz_0_x_1 = pbuffer.data(idx_eri_1_sksp + 81);

    auto g_0_xzzzzzz_0_y_1 = pbuffer.data(idx_eri_1_sksp + 82);

    auto g_0_xzzzzzz_0_z_1 = pbuffer.data(idx_eri_1_sksp + 83);

    auto g_0_yyyyyyy_0_x_1 = pbuffer.data(idx_eri_1_sksp + 84);

    auto g_0_yyyyyyy_0_y_1 = pbuffer.data(idx_eri_1_sksp + 85);

    auto g_0_yyyyyyy_0_z_1 = pbuffer.data(idx_eri_1_sksp + 86);

    auto g_0_yyyyyyz_0_y_1 = pbuffer.data(idx_eri_1_sksp + 88);

    auto g_0_yyyyyyz_0_z_1 = pbuffer.data(idx_eri_1_sksp + 89);

    auto g_0_yyyyyzz_0_x_1 = pbuffer.data(idx_eri_1_sksp + 90);

    auto g_0_yyyyyzz_0_y_1 = pbuffer.data(idx_eri_1_sksp + 91);

    auto g_0_yyyyyzz_0_z_1 = pbuffer.data(idx_eri_1_sksp + 92);

    auto g_0_yyyyzzz_0_x_1 = pbuffer.data(idx_eri_1_sksp + 93);

    auto g_0_yyyyzzz_0_y_1 = pbuffer.data(idx_eri_1_sksp + 94);

    auto g_0_yyyyzzz_0_z_1 = pbuffer.data(idx_eri_1_sksp + 95);

    auto g_0_yyyzzzz_0_x_1 = pbuffer.data(idx_eri_1_sksp + 96);

    auto g_0_yyyzzzz_0_y_1 = pbuffer.data(idx_eri_1_sksp + 97);

    auto g_0_yyyzzzz_0_z_1 = pbuffer.data(idx_eri_1_sksp + 98);

    auto g_0_yyzzzzz_0_x_1 = pbuffer.data(idx_eri_1_sksp + 99);

    auto g_0_yyzzzzz_0_y_1 = pbuffer.data(idx_eri_1_sksp + 100);

    auto g_0_yyzzzzz_0_z_1 = pbuffer.data(idx_eri_1_sksp + 101);

    auto g_0_yzzzzzz_0_x_1 = pbuffer.data(idx_eri_1_sksp + 102);

    auto g_0_yzzzzzz_0_y_1 = pbuffer.data(idx_eri_1_sksp + 103);

    auto g_0_yzzzzzz_0_z_1 = pbuffer.data(idx_eri_1_sksp + 104);

    auto g_0_zzzzzzz_0_x_1 = pbuffer.data(idx_eri_1_sksp + 105);

    auto g_0_zzzzzzz_0_y_1 = pbuffer.data(idx_eri_1_sksp + 106);

    auto g_0_zzzzzzz_0_z_1 = pbuffer.data(idx_eri_1_sksp + 107);

    /// Set up 0-3 components of targeted buffer : SLSP

    auto g_0_xxxxxxxx_0_x_0 = pbuffer.data(idx_eri_0_slsp);

    auto g_0_xxxxxxxx_0_y_0 = pbuffer.data(idx_eri_0_slsp + 1);

    auto g_0_xxxxxxxx_0_z_0 = pbuffer.data(idx_eri_0_slsp + 2);

    #pragma omp simd aligned(g_0_xxxxxx_0_x_0, g_0_xxxxxx_0_x_1, g_0_xxxxxx_0_y_0, g_0_xxxxxx_0_y_1, g_0_xxxxxx_0_z_0, g_0_xxxxxx_0_z_1, g_0_xxxxxxx_0_0_1, g_0_xxxxxxx_0_x_0, g_0_xxxxxxx_0_x_1, g_0_xxxxxxx_0_y_0, g_0_xxxxxxx_0_y_1, g_0_xxxxxxx_0_z_0, g_0_xxxxxxx_0_z_1, g_0_xxxxxxxx_0_x_0, g_0_xxxxxxxx_0_y_0, g_0_xxxxxxxx_0_z_0, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxxxx_0_x_0[i] = 7.0 * g_0_xxxxxx_0_x_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_x_1[i] * fti_ab_0 + g_0_xxxxxxx_0_0_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_x_0[i] * pb_x + g_0_xxxxxxx_0_x_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_y_0[i] = 7.0 * g_0_xxxxxx_0_y_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_y_1[i] * fti_ab_0 + g_0_xxxxxxx_0_y_0[i] * pb_x + g_0_xxxxxxx_0_y_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_z_0[i] = 7.0 * g_0_xxxxxx_0_z_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_z_1[i] * fti_ab_0 + g_0_xxxxxxx_0_z_0[i] * pb_x + g_0_xxxxxxx_0_z_1[i] * wp_x[i];
    }

    /// Set up 3-6 components of targeted buffer : SLSP

    auto g_0_xxxxxxxy_0_x_0 = pbuffer.data(idx_eri_0_slsp + 3);

    auto g_0_xxxxxxxy_0_y_0 = pbuffer.data(idx_eri_0_slsp + 4);

    auto g_0_xxxxxxxy_0_z_0 = pbuffer.data(idx_eri_0_slsp + 5);

    #pragma omp simd aligned(g_0_xxxxxxx_0_0_1, g_0_xxxxxxx_0_x_0, g_0_xxxxxxx_0_x_1, g_0_xxxxxxx_0_y_0, g_0_xxxxxxx_0_y_1, g_0_xxxxxxx_0_z_0, g_0_xxxxxxx_0_z_1, g_0_xxxxxxxy_0_x_0, g_0_xxxxxxxy_0_y_0, g_0_xxxxxxxy_0_z_0, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxxxy_0_x_0[i] = g_0_xxxxxxx_0_x_0[i] * pb_y + g_0_xxxxxxx_0_x_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_y_0[i] = g_0_xxxxxxx_0_0_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_y_0[i] * pb_y + g_0_xxxxxxx_0_y_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_z_0[i] = g_0_xxxxxxx_0_z_0[i] * pb_y + g_0_xxxxxxx_0_z_1[i] * wp_y[i];
    }

    /// Set up 6-9 components of targeted buffer : SLSP

    auto g_0_xxxxxxxz_0_x_0 = pbuffer.data(idx_eri_0_slsp + 6);

    auto g_0_xxxxxxxz_0_y_0 = pbuffer.data(idx_eri_0_slsp + 7);

    auto g_0_xxxxxxxz_0_z_0 = pbuffer.data(idx_eri_0_slsp + 8);

    #pragma omp simd aligned(g_0_xxxxxxx_0_0_1, g_0_xxxxxxx_0_x_0, g_0_xxxxxxx_0_x_1, g_0_xxxxxxx_0_y_0, g_0_xxxxxxx_0_y_1, g_0_xxxxxxx_0_z_0, g_0_xxxxxxx_0_z_1, g_0_xxxxxxxz_0_x_0, g_0_xxxxxxxz_0_y_0, g_0_xxxxxxxz_0_z_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxxxz_0_x_0[i] = g_0_xxxxxxx_0_x_0[i] * pb_z + g_0_xxxxxxx_0_x_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_y_0[i] = g_0_xxxxxxx_0_y_0[i] * pb_z + g_0_xxxxxxx_0_y_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_z_0[i] = g_0_xxxxxxx_0_0_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_z_0[i] * pb_z + g_0_xxxxxxx_0_z_1[i] * wp_z[i];
    }

    /// Set up 9-12 components of targeted buffer : SLSP

    auto g_0_xxxxxxyy_0_x_0 = pbuffer.data(idx_eri_0_slsp + 9);

    auto g_0_xxxxxxyy_0_y_0 = pbuffer.data(idx_eri_0_slsp + 10);

    auto g_0_xxxxxxyy_0_z_0 = pbuffer.data(idx_eri_0_slsp + 11);

    #pragma omp simd aligned(g_0_xxxxxx_0_x_0, g_0_xxxxxx_0_x_1, g_0_xxxxxxy_0_x_0, g_0_xxxxxxy_0_x_1, g_0_xxxxxxyy_0_x_0, g_0_xxxxxxyy_0_y_0, g_0_xxxxxxyy_0_z_0, g_0_xxxxxyy_0_y_0, g_0_xxxxxyy_0_y_1, g_0_xxxxxyy_0_z_0, g_0_xxxxxyy_0_z_1, g_0_xxxxyy_0_y_0, g_0_xxxxyy_0_y_1, g_0_xxxxyy_0_z_0, g_0_xxxxyy_0_z_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxxyy_0_x_0[i] = g_0_xxxxxx_0_x_0[i] * fi_ab_0 - g_0_xxxxxx_0_x_1[i] * fti_ab_0 + g_0_xxxxxxy_0_x_0[i] * pb_y + g_0_xxxxxxy_0_x_1[i] * wp_y[i];

        g_0_xxxxxxyy_0_y_0[i] = 5.0 * g_0_xxxxyy_0_y_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_y_1[i] * fti_ab_0 + g_0_xxxxxyy_0_y_0[i] * pb_x + g_0_xxxxxyy_0_y_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_z_0[i] = 5.0 * g_0_xxxxyy_0_z_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_z_1[i] * fti_ab_0 + g_0_xxxxxyy_0_z_0[i] * pb_x + g_0_xxxxxyy_0_z_1[i] * wp_x[i];
    }

    /// Set up 12-15 components of targeted buffer : SLSP

    auto g_0_xxxxxxyz_0_x_0 = pbuffer.data(idx_eri_0_slsp + 12);

    auto g_0_xxxxxxyz_0_y_0 = pbuffer.data(idx_eri_0_slsp + 13);

    auto g_0_xxxxxxyz_0_z_0 = pbuffer.data(idx_eri_0_slsp + 14);

    #pragma omp simd aligned(g_0_xxxxxxy_0_y_0, g_0_xxxxxxy_0_y_1, g_0_xxxxxxyz_0_x_0, g_0_xxxxxxyz_0_y_0, g_0_xxxxxxyz_0_z_0, g_0_xxxxxxz_0_x_0, g_0_xxxxxxz_0_x_1, g_0_xxxxxxz_0_z_0, g_0_xxxxxxz_0_z_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_0_xxxxxxyz_0_x_0[i] = g_0_xxxxxxz_0_x_0[i] * pb_y + g_0_xxxxxxz_0_x_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_y_0[i] = g_0_xxxxxxy_0_y_0[i] * pb_z + g_0_xxxxxxy_0_y_1[i] * wp_z[i];

        g_0_xxxxxxyz_0_z_0[i] = g_0_xxxxxxz_0_z_0[i] * pb_y + g_0_xxxxxxz_0_z_1[i] * wp_y[i];
    }

    /// Set up 15-18 components of targeted buffer : SLSP

    auto g_0_xxxxxxzz_0_x_0 = pbuffer.data(idx_eri_0_slsp + 15);

    auto g_0_xxxxxxzz_0_y_0 = pbuffer.data(idx_eri_0_slsp + 16);

    auto g_0_xxxxxxzz_0_z_0 = pbuffer.data(idx_eri_0_slsp + 17);

    #pragma omp simd aligned(g_0_xxxxxx_0_x_0, g_0_xxxxxx_0_x_1, g_0_xxxxxxz_0_x_0, g_0_xxxxxxz_0_x_1, g_0_xxxxxxzz_0_x_0, g_0_xxxxxxzz_0_y_0, g_0_xxxxxxzz_0_z_0, g_0_xxxxxzz_0_y_0, g_0_xxxxxzz_0_y_1, g_0_xxxxxzz_0_z_0, g_0_xxxxxzz_0_z_1, g_0_xxxxzz_0_y_0, g_0_xxxxzz_0_y_1, g_0_xxxxzz_0_z_0, g_0_xxxxzz_0_z_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxxzz_0_x_0[i] = g_0_xxxxxx_0_x_0[i] * fi_ab_0 - g_0_xxxxxx_0_x_1[i] * fti_ab_0 + g_0_xxxxxxz_0_x_0[i] * pb_z + g_0_xxxxxxz_0_x_1[i] * wp_z[i];

        g_0_xxxxxxzz_0_y_0[i] = 5.0 * g_0_xxxxzz_0_y_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_y_1[i] * fti_ab_0 + g_0_xxxxxzz_0_y_0[i] * pb_x + g_0_xxxxxzz_0_y_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_z_0[i] = 5.0 * g_0_xxxxzz_0_z_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_z_1[i] * fti_ab_0 + g_0_xxxxxzz_0_z_0[i] * pb_x + g_0_xxxxxzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 18-21 components of targeted buffer : SLSP

    auto g_0_xxxxxyyy_0_x_0 = pbuffer.data(idx_eri_0_slsp + 18);

    auto g_0_xxxxxyyy_0_y_0 = pbuffer.data(idx_eri_0_slsp + 19);

    auto g_0_xxxxxyyy_0_z_0 = pbuffer.data(idx_eri_0_slsp + 20);

    #pragma omp simd aligned(g_0_xxxxxy_0_x_0, g_0_xxxxxy_0_x_1, g_0_xxxxxyy_0_x_0, g_0_xxxxxyy_0_x_1, g_0_xxxxxyyy_0_x_0, g_0_xxxxxyyy_0_y_0, g_0_xxxxxyyy_0_z_0, g_0_xxxxyyy_0_y_0, g_0_xxxxyyy_0_y_1, g_0_xxxxyyy_0_z_0, g_0_xxxxyyy_0_z_1, g_0_xxxyyy_0_y_0, g_0_xxxyyy_0_y_1, g_0_xxxyyy_0_z_0, g_0_xxxyyy_0_z_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxyyy_0_x_0[i] = 2.0 * g_0_xxxxxy_0_x_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxy_0_x_1[i] * fti_ab_0 + g_0_xxxxxyy_0_x_0[i] * pb_y + g_0_xxxxxyy_0_x_1[i] * wp_y[i];

        g_0_xxxxxyyy_0_y_0[i] = 4.0 * g_0_xxxyyy_0_y_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_y_1[i] * fti_ab_0 + g_0_xxxxyyy_0_y_0[i] * pb_x + g_0_xxxxyyy_0_y_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_z_0[i] = 4.0 * g_0_xxxyyy_0_z_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_z_1[i] * fti_ab_0 + g_0_xxxxyyy_0_z_0[i] * pb_x + g_0_xxxxyyy_0_z_1[i] * wp_x[i];
    }

    /// Set up 21-24 components of targeted buffer : SLSP

    auto g_0_xxxxxyyz_0_x_0 = pbuffer.data(idx_eri_0_slsp + 21);

    auto g_0_xxxxxyyz_0_y_0 = pbuffer.data(idx_eri_0_slsp + 22);

    auto g_0_xxxxxyyz_0_z_0 = pbuffer.data(idx_eri_0_slsp + 23);

    #pragma omp simd aligned(g_0_xxxxxyy_0_0_1, g_0_xxxxxyy_0_x_0, g_0_xxxxxyy_0_x_1, g_0_xxxxxyy_0_y_0, g_0_xxxxxyy_0_y_1, g_0_xxxxxyy_0_z_0, g_0_xxxxxyy_0_z_1, g_0_xxxxxyyz_0_x_0, g_0_xxxxxyyz_0_y_0, g_0_xxxxxyyz_0_z_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxyyz_0_x_0[i] = g_0_xxxxxyy_0_x_0[i] * pb_z + g_0_xxxxxyy_0_x_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_y_0[i] = g_0_xxxxxyy_0_y_0[i] * pb_z + g_0_xxxxxyy_0_y_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_z_0[i] = g_0_xxxxxyy_0_0_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_z_0[i] * pb_z + g_0_xxxxxyy_0_z_1[i] * wp_z[i];
    }

    /// Set up 24-27 components of targeted buffer : SLSP

    auto g_0_xxxxxyzz_0_x_0 = pbuffer.data(idx_eri_0_slsp + 24);

    auto g_0_xxxxxyzz_0_y_0 = pbuffer.data(idx_eri_0_slsp + 25);

    auto g_0_xxxxxyzz_0_z_0 = pbuffer.data(idx_eri_0_slsp + 26);

    #pragma omp simd aligned(g_0_xxxxxyzz_0_x_0, g_0_xxxxxyzz_0_y_0, g_0_xxxxxyzz_0_z_0, g_0_xxxxxzz_0_0_1, g_0_xxxxxzz_0_x_0, g_0_xxxxxzz_0_x_1, g_0_xxxxxzz_0_y_0, g_0_xxxxxzz_0_y_1, g_0_xxxxxzz_0_z_0, g_0_xxxxxzz_0_z_1, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxyzz_0_x_0[i] = g_0_xxxxxzz_0_x_0[i] * pb_y + g_0_xxxxxzz_0_x_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_y_0[i] = g_0_xxxxxzz_0_0_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_y_0[i] * pb_y + g_0_xxxxxzz_0_y_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_z_0[i] = g_0_xxxxxzz_0_z_0[i] * pb_y + g_0_xxxxxzz_0_z_1[i] * wp_y[i];
    }

    /// Set up 27-30 components of targeted buffer : SLSP

    auto g_0_xxxxxzzz_0_x_0 = pbuffer.data(idx_eri_0_slsp + 27);

    auto g_0_xxxxxzzz_0_y_0 = pbuffer.data(idx_eri_0_slsp + 28);

    auto g_0_xxxxxzzz_0_z_0 = pbuffer.data(idx_eri_0_slsp + 29);

    #pragma omp simd aligned(g_0_xxxxxz_0_x_0, g_0_xxxxxz_0_x_1, g_0_xxxxxzz_0_x_0, g_0_xxxxxzz_0_x_1, g_0_xxxxxzzz_0_x_0, g_0_xxxxxzzz_0_y_0, g_0_xxxxxzzz_0_z_0, g_0_xxxxzzz_0_y_0, g_0_xxxxzzz_0_y_1, g_0_xxxxzzz_0_z_0, g_0_xxxxzzz_0_z_1, g_0_xxxzzz_0_y_0, g_0_xxxzzz_0_y_1, g_0_xxxzzz_0_z_0, g_0_xxxzzz_0_z_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxzzz_0_x_0[i] = 2.0 * g_0_xxxxxz_0_x_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxz_0_x_1[i] * fti_ab_0 + g_0_xxxxxzz_0_x_0[i] * pb_z + g_0_xxxxxzz_0_x_1[i] * wp_z[i];

        g_0_xxxxxzzz_0_y_0[i] = 4.0 * g_0_xxxzzz_0_y_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_y_1[i] * fti_ab_0 + g_0_xxxxzzz_0_y_0[i] * pb_x + g_0_xxxxzzz_0_y_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_z_0[i] = 4.0 * g_0_xxxzzz_0_z_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_z_1[i] * fti_ab_0 + g_0_xxxxzzz_0_z_0[i] * pb_x + g_0_xxxxzzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 30-33 components of targeted buffer : SLSP

    auto g_0_xxxxyyyy_0_x_0 = pbuffer.data(idx_eri_0_slsp + 30);

    auto g_0_xxxxyyyy_0_y_0 = pbuffer.data(idx_eri_0_slsp + 31);

    auto g_0_xxxxyyyy_0_z_0 = pbuffer.data(idx_eri_0_slsp + 32);

    #pragma omp simd aligned(g_0_xxxxyy_0_x_0, g_0_xxxxyy_0_x_1, g_0_xxxxyyy_0_x_0, g_0_xxxxyyy_0_x_1, g_0_xxxxyyyy_0_x_0, g_0_xxxxyyyy_0_y_0, g_0_xxxxyyyy_0_z_0, g_0_xxxyyyy_0_y_0, g_0_xxxyyyy_0_y_1, g_0_xxxyyyy_0_z_0, g_0_xxxyyyy_0_z_1, g_0_xxyyyy_0_y_0, g_0_xxyyyy_0_y_1, g_0_xxyyyy_0_z_0, g_0_xxyyyy_0_z_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxyyyy_0_x_0[i] = 3.0 * g_0_xxxxyy_0_x_0[i] * fi_ab_0 - 3.0 * g_0_xxxxyy_0_x_1[i] * fti_ab_0 + g_0_xxxxyyy_0_x_0[i] * pb_y + g_0_xxxxyyy_0_x_1[i] * wp_y[i];

        g_0_xxxxyyyy_0_y_0[i] = 3.0 * g_0_xxyyyy_0_y_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_y_1[i] * fti_ab_0 + g_0_xxxyyyy_0_y_0[i] * pb_x + g_0_xxxyyyy_0_y_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_z_0[i] = 3.0 * g_0_xxyyyy_0_z_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_z_1[i] * fti_ab_0 + g_0_xxxyyyy_0_z_0[i] * pb_x + g_0_xxxyyyy_0_z_1[i] * wp_x[i];
    }

    /// Set up 33-36 components of targeted buffer : SLSP

    auto g_0_xxxxyyyz_0_x_0 = pbuffer.data(idx_eri_0_slsp + 33);

    auto g_0_xxxxyyyz_0_y_0 = pbuffer.data(idx_eri_0_slsp + 34);

    auto g_0_xxxxyyyz_0_z_0 = pbuffer.data(idx_eri_0_slsp + 35);

    #pragma omp simd aligned(g_0_xxxxyyy_0_0_1, g_0_xxxxyyy_0_x_0, g_0_xxxxyyy_0_x_1, g_0_xxxxyyy_0_y_0, g_0_xxxxyyy_0_y_1, g_0_xxxxyyy_0_z_0, g_0_xxxxyyy_0_z_1, g_0_xxxxyyyz_0_x_0, g_0_xxxxyyyz_0_y_0, g_0_xxxxyyyz_0_z_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxyyyz_0_x_0[i] = g_0_xxxxyyy_0_x_0[i] * pb_z + g_0_xxxxyyy_0_x_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_y_0[i] = g_0_xxxxyyy_0_y_0[i] * pb_z + g_0_xxxxyyy_0_y_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_z_0[i] = g_0_xxxxyyy_0_0_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_z_0[i] * pb_z + g_0_xxxxyyy_0_z_1[i] * wp_z[i];
    }

    /// Set up 36-39 components of targeted buffer : SLSP

    auto g_0_xxxxyyzz_0_x_0 = pbuffer.data(idx_eri_0_slsp + 36);

    auto g_0_xxxxyyzz_0_y_0 = pbuffer.data(idx_eri_0_slsp + 37);

    auto g_0_xxxxyyzz_0_z_0 = pbuffer.data(idx_eri_0_slsp + 38);

    #pragma omp simd aligned(g_0_xxxxyyzz_0_x_0, g_0_xxxxyyzz_0_y_0, g_0_xxxxyyzz_0_z_0, g_0_xxxxyzz_0_x_0, g_0_xxxxyzz_0_x_1, g_0_xxxxzz_0_x_0, g_0_xxxxzz_0_x_1, g_0_xxxyyzz_0_y_0, g_0_xxxyyzz_0_y_1, g_0_xxxyyzz_0_z_0, g_0_xxxyyzz_0_z_1, g_0_xxyyzz_0_y_0, g_0_xxyyzz_0_y_1, g_0_xxyyzz_0_z_0, g_0_xxyyzz_0_z_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxyyzz_0_x_0[i] = g_0_xxxxzz_0_x_0[i] * fi_ab_0 - g_0_xxxxzz_0_x_1[i] * fti_ab_0 + g_0_xxxxyzz_0_x_0[i] * pb_y + g_0_xxxxyzz_0_x_1[i] * wp_y[i];

        g_0_xxxxyyzz_0_y_0[i] = 3.0 * g_0_xxyyzz_0_y_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_y_1[i] * fti_ab_0 + g_0_xxxyyzz_0_y_0[i] * pb_x + g_0_xxxyyzz_0_y_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_z_0[i] = 3.0 * g_0_xxyyzz_0_z_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_z_1[i] * fti_ab_0 + g_0_xxxyyzz_0_z_0[i] * pb_x + g_0_xxxyyzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 39-42 components of targeted buffer : SLSP

    auto g_0_xxxxyzzz_0_x_0 = pbuffer.data(idx_eri_0_slsp + 39);

    auto g_0_xxxxyzzz_0_y_0 = pbuffer.data(idx_eri_0_slsp + 40);

    auto g_0_xxxxyzzz_0_z_0 = pbuffer.data(idx_eri_0_slsp + 41);

    #pragma omp simd aligned(g_0_xxxxyzzz_0_x_0, g_0_xxxxyzzz_0_y_0, g_0_xxxxyzzz_0_z_0, g_0_xxxxzzz_0_0_1, g_0_xxxxzzz_0_x_0, g_0_xxxxzzz_0_x_1, g_0_xxxxzzz_0_y_0, g_0_xxxxzzz_0_y_1, g_0_xxxxzzz_0_z_0, g_0_xxxxzzz_0_z_1, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxyzzz_0_x_0[i] = g_0_xxxxzzz_0_x_0[i] * pb_y + g_0_xxxxzzz_0_x_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_y_0[i] = g_0_xxxxzzz_0_0_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_y_0[i] * pb_y + g_0_xxxxzzz_0_y_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_z_0[i] = g_0_xxxxzzz_0_z_0[i] * pb_y + g_0_xxxxzzz_0_z_1[i] * wp_y[i];
    }

    /// Set up 42-45 components of targeted buffer : SLSP

    auto g_0_xxxxzzzz_0_x_0 = pbuffer.data(idx_eri_0_slsp + 42);

    auto g_0_xxxxzzzz_0_y_0 = pbuffer.data(idx_eri_0_slsp + 43);

    auto g_0_xxxxzzzz_0_z_0 = pbuffer.data(idx_eri_0_slsp + 44);

    #pragma omp simd aligned(g_0_xxxxzz_0_x_0, g_0_xxxxzz_0_x_1, g_0_xxxxzzz_0_x_0, g_0_xxxxzzz_0_x_1, g_0_xxxxzzzz_0_x_0, g_0_xxxxzzzz_0_y_0, g_0_xxxxzzzz_0_z_0, g_0_xxxzzzz_0_y_0, g_0_xxxzzzz_0_y_1, g_0_xxxzzzz_0_z_0, g_0_xxxzzzz_0_z_1, g_0_xxzzzz_0_y_0, g_0_xxzzzz_0_y_1, g_0_xxzzzz_0_z_0, g_0_xxzzzz_0_z_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxzzzz_0_x_0[i] = 3.0 * g_0_xxxxzz_0_x_0[i] * fi_ab_0 - 3.0 * g_0_xxxxzz_0_x_1[i] * fti_ab_0 + g_0_xxxxzzz_0_x_0[i] * pb_z + g_0_xxxxzzz_0_x_1[i] * wp_z[i];

        g_0_xxxxzzzz_0_y_0[i] = 3.0 * g_0_xxzzzz_0_y_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_y_1[i] * fti_ab_0 + g_0_xxxzzzz_0_y_0[i] * pb_x + g_0_xxxzzzz_0_y_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_z_0[i] = 3.0 * g_0_xxzzzz_0_z_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_z_1[i] * fti_ab_0 + g_0_xxxzzzz_0_z_0[i] * pb_x + g_0_xxxzzzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 45-48 components of targeted buffer : SLSP

    auto g_0_xxxyyyyy_0_x_0 = pbuffer.data(idx_eri_0_slsp + 45);

    auto g_0_xxxyyyyy_0_y_0 = pbuffer.data(idx_eri_0_slsp + 46);

    auto g_0_xxxyyyyy_0_z_0 = pbuffer.data(idx_eri_0_slsp + 47);

    #pragma omp simd aligned(g_0_xxxyyy_0_x_0, g_0_xxxyyy_0_x_1, g_0_xxxyyyy_0_x_0, g_0_xxxyyyy_0_x_1, g_0_xxxyyyyy_0_x_0, g_0_xxxyyyyy_0_y_0, g_0_xxxyyyyy_0_z_0, g_0_xxyyyyy_0_y_0, g_0_xxyyyyy_0_y_1, g_0_xxyyyyy_0_z_0, g_0_xxyyyyy_0_z_1, g_0_xyyyyy_0_y_0, g_0_xyyyyy_0_y_1, g_0_xyyyyy_0_z_0, g_0_xyyyyy_0_z_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyyyyy_0_x_0[i] = 4.0 * g_0_xxxyyy_0_x_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_x_1[i] * fti_ab_0 + g_0_xxxyyyy_0_x_0[i] * pb_y + g_0_xxxyyyy_0_x_1[i] * wp_y[i];

        g_0_xxxyyyyy_0_y_0[i] = 2.0 * g_0_xyyyyy_0_y_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_y_1[i] * fti_ab_0 + g_0_xxyyyyy_0_y_0[i] * pb_x + g_0_xxyyyyy_0_y_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_z_0[i] = 2.0 * g_0_xyyyyy_0_z_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_z_1[i] * fti_ab_0 + g_0_xxyyyyy_0_z_0[i] * pb_x + g_0_xxyyyyy_0_z_1[i] * wp_x[i];
    }

    /// Set up 48-51 components of targeted buffer : SLSP

    auto g_0_xxxyyyyz_0_x_0 = pbuffer.data(idx_eri_0_slsp + 48);

    auto g_0_xxxyyyyz_0_y_0 = pbuffer.data(idx_eri_0_slsp + 49);

    auto g_0_xxxyyyyz_0_z_0 = pbuffer.data(idx_eri_0_slsp + 50);

    #pragma omp simd aligned(g_0_xxxyyyy_0_0_1, g_0_xxxyyyy_0_x_0, g_0_xxxyyyy_0_x_1, g_0_xxxyyyy_0_y_0, g_0_xxxyyyy_0_y_1, g_0_xxxyyyy_0_z_0, g_0_xxxyyyy_0_z_1, g_0_xxxyyyyz_0_x_0, g_0_xxxyyyyz_0_y_0, g_0_xxxyyyyz_0_z_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyyyyz_0_x_0[i] = g_0_xxxyyyy_0_x_0[i] * pb_z + g_0_xxxyyyy_0_x_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_y_0[i] = g_0_xxxyyyy_0_y_0[i] * pb_z + g_0_xxxyyyy_0_y_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_z_0[i] = g_0_xxxyyyy_0_0_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_z_0[i] * pb_z + g_0_xxxyyyy_0_z_1[i] * wp_z[i];
    }

    /// Set up 51-54 components of targeted buffer : SLSP

    auto g_0_xxxyyyzz_0_x_0 = pbuffer.data(idx_eri_0_slsp + 51);

    auto g_0_xxxyyyzz_0_y_0 = pbuffer.data(idx_eri_0_slsp + 52);

    auto g_0_xxxyyyzz_0_z_0 = pbuffer.data(idx_eri_0_slsp + 53);

    #pragma omp simd aligned(g_0_xxxyyyzz_0_x_0, g_0_xxxyyyzz_0_y_0, g_0_xxxyyyzz_0_z_0, g_0_xxxyyzz_0_x_0, g_0_xxxyyzz_0_x_1, g_0_xxxyzz_0_x_0, g_0_xxxyzz_0_x_1, g_0_xxyyyzz_0_y_0, g_0_xxyyyzz_0_y_1, g_0_xxyyyzz_0_z_0, g_0_xxyyyzz_0_z_1, g_0_xyyyzz_0_y_0, g_0_xyyyzz_0_y_1, g_0_xyyyzz_0_z_0, g_0_xyyyzz_0_z_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyyyzz_0_x_0[i] = 2.0 * g_0_xxxyzz_0_x_0[i] * fi_ab_0 - 2.0 * g_0_xxxyzz_0_x_1[i] * fti_ab_0 + g_0_xxxyyzz_0_x_0[i] * pb_y + g_0_xxxyyzz_0_x_1[i] * wp_y[i];

        g_0_xxxyyyzz_0_y_0[i] = 2.0 * g_0_xyyyzz_0_y_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_y_1[i] * fti_ab_0 + g_0_xxyyyzz_0_y_0[i] * pb_x + g_0_xxyyyzz_0_y_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_z_0[i] = 2.0 * g_0_xyyyzz_0_z_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_z_1[i] * fti_ab_0 + g_0_xxyyyzz_0_z_0[i] * pb_x + g_0_xxyyyzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 54-57 components of targeted buffer : SLSP

    auto g_0_xxxyyzzz_0_x_0 = pbuffer.data(idx_eri_0_slsp + 54);

    auto g_0_xxxyyzzz_0_y_0 = pbuffer.data(idx_eri_0_slsp + 55);

    auto g_0_xxxyyzzz_0_z_0 = pbuffer.data(idx_eri_0_slsp + 56);

    #pragma omp simd aligned(g_0_xxxyyzzz_0_x_0, g_0_xxxyyzzz_0_y_0, g_0_xxxyyzzz_0_z_0, g_0_xxxyzzz_0_x_0, g_0_xxxyzzz_0_x_1, g_0_xxxzzz_0_x_0, g_0_xxxzzz_0_x_1, g_0_xxyyzzz_0_y_0, g_0_xxyyzzz_0_y_1, g_0_xxyyzzz_0_z_0, g_0_xxyyzzz_0_z_1, g_0_xyyzzz_0_y_0, g_0_xyyzzz_0_y_1, g_0_xyyzzz_0_z_0, g_0_xyyzzz_0_z_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyyzzz_0_x_0[i] = g_0_xxxzzz_0_x_0[i] * fi_ab_0 - g_0_xxxzzz_0_x_1[i] * fti_ab_0 + g_0_xxxyzzz_0_x_0[i] * pb_y + g_0_xxxyzzz_0_x_1[i] * wp_y[i];

        g_0_xxxyyzzz_0_y_0[i] = 2.0 * g_0_xyyzzz_0_y_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_y_1[i] * fti_ab_0 + g_0_xxyyzzz_0_y_0[i] * pb_x + g_0_xxyyzzz_0_y_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_z_0[i] = 2.0 * g_0_xyyzzz_0_z_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_z_1[i] * fti_ab_0 + g_0_xxyyzzz_0_z_0[i] * pb_x + g_0_xxyyzzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 57-60 components of targeted buffer : SLSP

    auto g_0_xxxyzzzz_0_x_0 = pbuffer.data(idx_eri_0_slsp + 57);

    auto g_0_xxxyzzzz_0_y_0 = pbuffer.data(idx_eri_0_slsp + 58);

    auto g_0_xxxyzzzz_0_z_0 = pbuffer.data(idx_eri_0_slsp + 59);

    #pragma omp simd aligned(g_0_xxxyzzzz_0_x_0, g_0_xxxyzzzz_0_y_0, g_0_xxxyzzzz_0_z_0, g_0_xxxzzzz_0_0_1, g_0_xxxzzzz_0_x_0, g_0_xxxzzzz_0_x_1, g_0_xxxzzzz_0_y_0, g_0_xxxzzzz_0_y_1, g_0_xxxzzzz_0_z_0, g_0_xxxzzzz_0_z_1, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyzzzz_0_x_0[i] = g_0_xxxzzzz_0_x_0[i] * pb_y + g_0_xxxzzzz_0_x_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_y_0[i] = g_0_xxxzzzz_0_0_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_y_0[i] * pb_y + g_0_xxxzzzz_0_y_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_z_0[i] = g_0_xxxzzzz_0_z_0[i] * pb_y + g_0_xxxzzzz_0_z_1[i] * wp_y[i];
    }

    /// Set up 60-63 components of targeted buffer : SLSP

    auto g_0_xxxzzzzz_0_x_0 = pbuffer.data(idx_eri_0_slsp + 60);

    auto g_0_xxxzzzzz_0_y_0 = pbuffer.data(idx_eri_0_slsp + 61);

    auto g_0_xxxzzzzz_0_z_0 = pbuffer.data(idx_eri_0_slsp + 62);

    #pragma omp simd aligned(g_0_xxxzzz_0_x_0, g_0_xxxzzz_0_x_1, g_0_xxxzzzz_0_x_0, g_0_xxxzzzz_0_x_1, g_0_xxxzzzzz_0_x_0, g_0_xxxzzzzz_0_y_0, g_0_xxxzzzzz_0_z_0, g_0_xxzzzzz_0_y_0, g_0_xxzzzzz_0_y_1, g_0_xxzzzzz_0_z_0, g_0_xxzzzzz_0_z_1, g_0_xzzzzz_0_y_0, g_0_xzzzzz_0_y_1, g_0_xzzzzz_0_z_0, g_0_xzzzzz_0_z_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxzzzzz_0_x_0[i] = 4.0 * g_0_xxxzzz_0_x_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_x_1[i] * fti_ab_0 + g_0_xxxzzzz_0_x_0[i] * pb_z + g_0_xxxzzzz_0_x_1[i] * wp_z[i];

        g_0_xxxzzzzz_0_y_0[i] = 2.0 * g_0_xzzzzz_0_y_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_y_1[i] * fti_ab_0 + g_0_xxzzzzz_0_y_0[i] * pb_x + g_0_xxzzzzz_0_y_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_z_0[i] = 2.0 * g_0_xzzzzz_0_z_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_z_1[i] * fti_ab_0 + g_0_xxzzzzz_0_z_0[i] * pb_x + g_0_xxzzzzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 63-66 components of targeted buffer : SLSP

    auto g_0_xxyyyyyy_0_x_0 = pbuffer.data(idx_eri_0_slsp + 63);

    auto g_0_xxyyyyyy_0_y_0 = pbuffer.data(idx_eri_0_slsp + 64);

    auto g_0_xxyyyyyy_0_z_0 = pbuffer.data(idx_eri_0_slsp + 65);

    #pragma omp simd aligned(g_0_xxyyyy_0_x_0, g_0_xxyyyy_0_x_1, g_0_xxyyyyy_0_x_0, g_0_xxyyyyy_0_x_1, g_0_xxyyyyyy_0_x_0, g_0_xxyyyyyy_0_y_0, g_0_xxyyyyyy_0_z_0, g_0_xyyyyyy_0_y_0, g_0_xyyyyyy_0_y_1, g_0_xyyyyyy_0_z_0, g_0_xyyyyyy_0_z_1, g_0_yyyyyy_0_y_0, g_0_yyyyyy_0_y_1, g_0_yyyyyy_0_z_0, g_0_yyyyyy_0_z_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyyyyy_0_x_0[i] = 5.0 * g_0_xxyyyy_0_x_0[i] * fi_ab_0 - 5.0 * g_0_xxyyyy_0_x_1[i] * fti_ab_0 + g_0_xxyyyyy_0_x_0[i] * pb_y + g_0_xxyyyyy_0_x_1[i] * wp_y[i];

        g_0_xxyyyyyy_0_y_0[i] = g_0_yyyyyy_0_y_0[i] * fi_ab_0 - g_0_yyyyyy_0_y_1[i] * fti_ab_0 + g_0_xyyyyyy_0_y_0[i] * pb_x + g_0_xyyyyyy_0_y_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_z_0[i] = g_0_yyyyyy_0_z_0[i] * fi_ab_0 - g_0_yyyyyy_0_z_1[i] * fti_ab_0 + g_0_xyyyyyy_0_z_0[i] * pb_x + g_0_xyyyyyy_0_z_1[i] * wp_x[i];
    }

    /// Set up 66-69 components of targeted buffer : SLSP

    auto g_0_xxyyyyyz_0_x_0 = pbuffer.data(idx_eri_0_slsp + 66);

    auto g_0_xxyyyyyz_0_y_0 = pbuffer.data(idx_eri_0_slsp + 67);

    auto g_0_xxyyyyyz_0_z_0 = pbuffer.data(idx_eri_0_slsp + 68);

    #pragma omp simd aligned(g_0_xxyyyyy_0_0_1, g_0_xxyyyyy_0_x_0, g_0_xxyyyyy_0_x_1, g_0_xxyyyyy_0_y_0, g_0_xxyyyyy_0_y_1, g_0_xxyyyyy_0_z_0, g_0_xxyyyyy_0_z_1, g_0_xxyyyyyz_0_x_0, g_0_xxyyyyyz_0_y_0, g_0_xxyyyyyz_0_z_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyyyyz_0_x_0[i] = g_0_xxyyyyy_0_x_0[i] * pb_z + g_0_xxyyyyy_0_x_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_y_0[i] = g_0_xxyyyyy_0_y_0[i] * pb_z + g_0_xxyyyyy_0_y_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_z_0[i] = g_0_xxyyyyy_0_0_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_z_0[i] * pb_z + g_0_xxyyyyy_0_z_1[i] * wp_z[i];
    }

    /// Set up 69-72 components of targeted buffer : SLSP

    auto g_0_xxyyyyzz_0_x_0 = pbuffer.data(idx_eri_0_slsp + 69);

    auto g_0_xxyyyyzz_0_y_0 = pbuffer.data(idx_eri_0_slsp + 70);

    auto g_0_xxyyyyzz_0_z_0 = pbuffer.data(idx_eri_0_slsp + 71);

    #pragma omp simd aligned(g_0_xxyyyyzz_0_x_0, g_0_xxyyyyzz_0_y_0, g_0_xxyyyyzz_0_z_0, g_0_xxyyyzz_0_x_0, g_0_xxyyyzz_0_x_1, g_0_xxyyzz_0_x_0, g_0_xxyyzz_0_x_1, g_0_xyyyyzz_0_y_0, g_0_xyyyyzz_0_y_1, g_0_xyyyyzz_0_z_0, g_0_xyyyyzz_0_z_1, g_0_yyyyzz_0_y_0, g_0_yyyyzz_0_y_1, g_0_yyyyzz_0_z_0, g_0_yyyyzz_0_z_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyyyzz_0_x_0[i] = 3.0 * g_0_xxyyzz_0_x_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_x_1[i] * fti_ab_0 + g_0_xxyyyzz_0_x_0[i] * pb_y + g_0_xxyyyzz_0_x_1[i] * wp_y[i];

        g_0_xxyyyyzz_0_y_0[i] = g_0_yyyyzz_0_y_0[i] * fi_ab_0 - g_0_yyyyzz_0_y_1[i] * fti_ab_0 + g_0_xyyyyzz_0_y_0[i] * pb_x + g_0_xyyyyzz_0_y_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_z_0[i] = g_0_yyyyzz_0_z_0[i] * fi_ab_0 - g_0_yyyyzz_0_z_1[i] * fti_ab_0 + g_0_xyyyyzz_0_z_0[i] * pb_x + g_0_xyyyyzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 72-75 components of targeted buffer : SLSP

    auto g_0_xxyyyzzz_0_x_0 = pbuffer.data(idx_eri_0_slsp + 72);

    auto g_0_xxyyyzzz_0_y_0 = pbuffer.data(idx_eri_0_slsp + 73);

    auto g_0_xxyyyzzz_0_z_0 = pbuffer.data(idx_eri_0_slsp + 74);

    #pragma omp simd aligned(g_0_xxyyyzzz_0_x_0, g_0_xxyyyzzz_0_y_0, g_0_xxyyyzzz_0_z_0, g_0_xxyyzzz_0_x_0, g_0_xxyyzzz_0_x_1, g_0_xxyzzz_0_x_0, g_0_xxyzzz_0_x_1, g_0_xyyyzzz_0_y_0, g_0_xyyyzzz_0_y_1, g_0_xyyyzzz_0_z_0, g_0_xyyyzzz_0_z_1, g_0_yyyzzz_0_y_0, g_0_yyyzzz_0_y_1, g_0_yyyzzz_0_z_0, g_0_yyyzzz_0_z_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyyzzz_0_x_0[i] = 2.0 * g_0_xxyzzz_0_x_0[i] * fi_ab_0 - 2.0 * g_0_xxyzzz_0_x_1[i] * fti_ab_0 + g_0_xxyyzzz_0_x_0[i] * pb_y + g_0_xxyyzzz_0_x_1[i] * wp_y[i];

        g_0_xxyyyzzz_0_y_0[i] = g_0_yyyzzz_0_y_0[i] * fi_ab_0 - g_0_yyyzzz_0_y_1[i] * fti_ab_0 + g_0_xyyyzzz_0_y_0[i] * pb_x + g_0_xyyyzzz_0_y_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_z_0[i] = g_0_yyyzzz_0_z_0[i] * fi_ab_0 - g_0_yyyzzz_0_z_1[i] * fti_ab_0 + g_0_xyyyzzz_0_z_0[i] * pb_x + g_0_xyyyzzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 75-78 components of targeted buffer : SLSP

    auto g_0_xxyyzzzz_0_x_0 = pbuffer.data(idx_eri_0_slsp + 75);

    auto g_0_xxyyzzzz_0_y_0 = pbuffer.data(idx_eri_0_slsp + 76);

    auto g_0_xxyyzzzz_0_z_0 = pbuffer.data(idx_eri_0_slsp + 77);

    #pragma omp simd aligned(g_0_xxyyzzzz_0_x_0, g_0_xxyyzzzz_0_y_0, g_0_xxyyzzzz_0_z_0, g_0_xxyzzzz_0_x_0, g_0_xxyzzzz_0_x_1, g_0_xxzzzz_0_x_0, g_0_xxzzzz_0_x_1, g_0_xyyzzzz_0_y_0, g_0_xyyzzzz_0_y_1, g_0_xyyzzzz_0_z_0, g_0_xyyzzzz_0_z_1, g_0_yyzzzz_0_y_0, g_0_yyzzzz_0_y_1, g_0_yyzzzz_0_z_0, g_0_yyzzzz_0_z_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyzzzz_0_x_0[i] = g_0_xxzzzz_0_x_0[i] * fi_ab_0 - g_0_xxzzzz_0_x_1[i] * fti_ab_0 + g_0_xxyzzzz_0_x_0[i] * pb_y + g_0_xxyzzzz_0_x_1[i] * wp_y[i];

        g_0_xxyyzzzz_0_y_0[i] = g_0_yyzzzz_0_y_0[i] * fi_ab_0 - g_0_yyzzzz_0_y_1[i] * fti_ab_0 + g_0_xyyzzzz_0_y_0[i] * pb_x + g_0_xyyzzzz_0_y_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_z_0[i] = g_0_yyzzzz_0_z_0[i] * fi_ab_0 - g_0_yyzzzz_0_z_1[i] * fti_ab_0 + g_0_xyyzzzz_0_z_0[i] * pb_x + g_0_xyyzzzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 78-81 components of targeted buffer : SLSP

    auto g_0_xxyzzzzz_0_x_0 = pbuffer.data(idx_eri_0_slsp + 78);

    auto g_0_xxyzzzzz_0_y_0 = pbuffer.data(idx_eri_0_slsp + 79);

    auto g_0_xxyzzzzz_0_z_0 = pbuffer.data(idx_eri_0_slsp + 80);

    #pragma omp simd aligned(g_0_xxyzzzzz_0_x_0, g_0_xxyzzzzz_0_y_0, g_0_xxyzzzzz_0_z_0, g_0_xxzzzzz_0_0_1, g_0_xxzzzzz_0_x_0, g_0_xxzzzzz_0_x_1, g_0_xxzzzzz_0_y_0, g_0_xxzzzzz_0_y_1, g_0_xxzzzzz_0_z_0, g_0_xxzzzzz_0_z_1, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyzzzzz_0_x_0[i] = g_0_xxzzzzz_0_x_0[i] * pb_y + g_0_xxzzzzz_0_x_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_y_0[i] = g_0_xxzzzzz_0_0_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_y_0[i] * pb_y + g_0_xxzzzzz_0_y_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_z_0[i] = g_0_xxzzzzz_0_z_0[i] * pb_y + g_0_xxzzzzz_0_z_1[i] * wp_y[i];
    }

    /// Set up 81-84 components of targeted buffer : SLSP

    auto g_0_xxzzzzzz_0_x_0 = pbuffer.data(idx_eri_0_slsp + 81);

    auto g_0_xxzzzzzz_0_y_0 = pbuffer.data(idx_eri_0_slsp + 82);

    auto g_0_xxzzzzzz_0_z_0 = pbuffer.data(idx_eri_0_slsp + 83);

    #pragma omp simd aligned(g_0_xxzzzz_0_x_0, g_0_xxzzzz_0_x_1, g_0_xxzzzzz_0_x_0, g_0_xxzzzzz_0_x_1, g_0_xxzzzzzz_0_x_0, g_0_xxzzzzzz_0_y_0, g_0_xxzzzzzz_0_z_0, g_0_xzzzzzz_0_y_0, g_0_xzzzzzz_0_y_1, g_0_xzzzzzz_0_z_0, g_0_xzzzzzz_0_z_1, g_0_zzzzzz_0_y_0, g_0_zzzzzz_0_y_1, g_0_zzzzzz_0_z_0, g_0_zzzzzz_0_z_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxzzzzzz_0_x_0[i] = 5.0 * g_0_xxzzzz_0_x_0[i] * fi_ab_0 - 5.0 * g_0_xxzzzz_0_x_1[i] * fti_ab_0 + g_0_xxzzzzz_0_x_0[i] * pb_z + g_0_xxzzzzz_0_x_1[i] * wp_z[i];

        g_0_xxzzzzzz_0_y_0[i] = g_0_zzzzzz_0_y_0[i] * fi_ab_0 - g_0_zzzzzz_0_y_1[i] * fti_ab_0 + g_0_xzzzzzz_0_y_0[i] * pb_x + g_0_xzzzzzz_0_y_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_z_0[i] = g_0_zzzzzz_0_z_0[i] * fi_ab_0 - g_0_zzzzzz_0_z_1[i] * fti_ab_0 + g_0_xzzzzzz_0_z_0[i] * pb_x + g_0_xzzzzzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 84-87 components of targeted buffer : SLSP

    auto g_0_xyyyyyyy_0_x_0 = pbuffer.data(idx_eri_0_slsp + 84);

    auto g_0_xyyyyyyy_0_y_0 = pbuffer.data(idx_eri_0_slsp + 85);

    auto g_0_xyyyyyyy_0_z_0 = pbuffer.data(idx_eri_0_slsp + 86);

    #pragma omp simd aligned(g_0_xyyyyyyy_0_x_0, g_0_xyyyyyyy_0_y_0, g_0_xyyyyyyy_0_z_0, g_0_yyyyyyy_0_0_1, g_0_yyyyyyy_0_x_0, g_0_yyyyyyy_0_x_1, g_0_yyyyyyy_0_y_0, g_0_yyyyyyy_0_y_1, g_0_yyyyyyy_0_z_0, g_0_yyyyyyy_0_z_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyyyy_0_x_0[i] = g_0_yyyyyyy_0_0_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_x_0[i] * pb_x + g_0_yyyyyyy_0_x_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_y_0[i] = g_0_yyyyyyy_0_y_0[i] * pb_x + g_0_yyyyyyy_0_y_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_z_0[i] = g_0_yyyyyyy_0_z_0[i] * pb_x + g_0_yyyyyyy_0_z_1[i] * wp_x[i];
    }

    /// Set up 87-90 components of targeted buffer : SLSP

    auto g_0_xyyyyyyz_0_x_0 = pbuffer.data(idx_eri_0_slsp + 87);

    auto g_0_xyyyyyyz_0_y_0 = pbuffer.data(idx_eri_0_slsp + 88);

    auto g_0_xyyyyyyz_0_z_0 = pbuffer.data(idx_eri_0_slsp + 89);

    #pragma omp simd aligned(g_0_xyyyyyy_0_x_0, g_0_xyyyyyy_0_x_1, g_0_xyyyyyyz_0_x_0, g_0_xyyyyyyz_0_y_0, g_0_xyyyyyyz_0_z_0, g_0_yyyyyyz_0_y_0, g_0_yyyyyyz_0_y_1, g_0_yyyyyyz_0_z_0, g_0_yyyyyyz_0_z_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_0_xyyyyyyz_0_x_0[i] = g_0_xyyyyyy_0_x_0[i] * pb_z + g_0_xyyyyyy_0_x_1[i] * wp_z[i];

        g_0_xyyyyyyz_0_y_0[i] = g_0_yyyyyyz_0_y_0[i] * pb_x + g_0_yyyyyyz_0_y_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_z_0[i] = g_0_yyyyyyz_0_z_0[i] * pb_x + g_0_yyyyyyz_0_z_1[i] * wp_x[i];
    }

    /// Set up 90-93 components of targeted buffer : SLSP

    auto g_0_xyyyyyzz_0_x_0 = pbuffer.data(idx_eri_0_slsp + 90);

    auto g_0_xyyyyyzz_0_y_0 = pbuffer.data(idx_eri_0_slsp + 91);

    auto g_0_xyyyyyzz_0_z_0 = pbuffer.data(idx_eri_0_slsp + 92);

    #pragma omp simd aligned(g_0_xyyyyyzz_0_x_0, g_0_xyyyyyzz_0_y_0, g_0_xyyyyyzz_0_z_0, g_0_yyyyyzz_0_0_1, g_0_yyyyyzz_0_x_0, g_0_yyyyyzz_0_x_1, g_0_yyyyyzz_0_y_0, g_0_yyyyyzz_0_y_1, g_0_yyyyyzz_0_z_0, g_0_yyyyyzz_0_z_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyyzz_0_x_0[i] = g_0_yyyyyzz_0_0_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_x_0[i] * pb_x + g_0_yyyyyzz_0_x_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_y_0[i] = g_0_yyyyyzz_0_y_0[i] * pb_x + g_0_yyyyyzz_0_y_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_z_0[i] = g_0_yyyyyzz_0_z_0[i] * pb_x + g_0_yyyyyzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 93-96 components of targeted buffer : SLSP

    auto g_0_xyyyyzzz_0_x_0 = pbuffer.data(idx_eri_0_slsp + 93);

    auto g_0_xyyyyzzz_0_y_0 = pbuffer.data(idx_eri_0_slsp + 94);

    auto g_0_xyyyyzzz_0_z_0 = pbuffer.data(idx_eri_0_slsp + 95);

    #pragma omp simd aligned(g_0_xyyyyzzz_0_x_0, g_0_xyyyyzzz_0_y_0, g_0_xyyyyzzz_0_z_0, g_0_yyyyzzz_0_0_1, g_0_yyyyzzz_0_x_0, g_0_yyyyzzz_0_x_1, g_0_yyyyzzz_0_y_0, g_0_yyyyzzz_0_y_1, g_0_yyyyzzz_0_z_0, g_0_yyyyzzz_0_z_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyzzz_0_x_0[i] = g_0_yyyyzzz_0_0_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_x_0[i] * pb_x + g_0_yyyyzzz_0_x_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_y_0[i] = g_0_yyyyzzz_0_y_0[i] * pb_x + g_0_yyyyzzz_0_y_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_z_0[i] = g_0_yyyyzzz_0_z_0[i] * pb_x + g_0_yyyyzzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 96-99 components of targeted buffer : SLSP

    auto g_0_xyyyzzzz_0_x_0 = pbuffer.data(idx_eri_0_slsp + 96);

    auto g_0_xyyyzzzz_0_y_0 = pbuffer.data(idx_eri_0_slsp + 97);

    auto g_0_xyyyzzzz_0_z_0 = pbuffer.data(idx_eri_0_slsp + 98);

    #pragma omp simd aligned(g_0_xyyyzzzz_0_x_0, g_0_xyyyzzzz_0_y_0, g_0_xyyyzzzz_0_z_0, g_0_yyyzzzz_0_0_1, g_0_yyyzzzz_0_x_0, g_0_yyyzzzz_0_x_1, g_0_yyyzzzz_0_y_0, g_0_yyyzzzz_0_y_1, g_0_yyyzzzz_0_z_0, g_0_yyyzzzz_0_z_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyzzzz_0_x_0[i] = g_0_yyyzzzz_0_0_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_x_0[i] * pb_x + g_0_yyyzzzz_0_x_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_y_0[i] = g_0_yyyzzzz_0_y_0[i] * pb_x + g_0_yyyzzzz_0_y_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_z_0[i] = g_0_yyyzzzz_0_z_0[i] * pb_x + g_0_yyyzzzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 99-102 components of targeted buffer : SLSP

    auto g_0_xyyzzzzz_0_x_0 = pbuffer.data(idx_eri_0_slsp + 99);

    auto g_0_xyyzzzzz_0_y_0 = pbuffer.data(idx_eri_0_slsp + 100);

    auto g_0_xyyzzzzz_0_z_0 = pbuffer.data(idx_eri_0_slsp + 101);

    #pragma omp simd aligned(g_0_xyyzzzzz_0_x_0, g_0_xyyzzzzz_0_y_0, g_0_xyyzzzzz_0_z_0, g_0_yyzzzzz_0_0_1, g_0_yyzzzzz_0_x_0, g_0_yyzzzzz_0_x_1, g_0_yyzzzzz_0_y_0, g_0_yyzzzzz_0_y_1, g_0_yyzzzzz_0_z_0, g_0_yyzzzzz_0_z_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyzzzzz_0_x_0[i] = g_0_yyzzzzz_0_0_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_x_0[i] * pb_x + g_0_yyzzzzz_0_x_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_y_0[i] = g_0_yyzzzzz_0_y_0[i] * pb_x + g_0_yyzzzzz_0_y_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_z_0[i] = g_0_yyzzzzz_0_z_0[i] * pb_x + g_0_yyzzzzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 102-105 components of targeted buffer : SLSP

    auto g_0_xyzzzzzz_0_x_0 = pbuffer.data(idx_eri_0_slsp + 102);

    auto g_0_xyzzzzzz_0_y_0 = pbuffer.data(idx_eri_0_slsp + 103);

    auto g_0_xyzzzzzz_0_z_0 = pbuffer.data(idx_eri_0_slsp + 104);

    #pragma omp simd aligned(g_0_xyzzzzzz_0_x_0, g_0_xyzzzzzz_0_y_0, g_0_xyzzzzzz_0_z_0, g_0_xzzzzzz_0_x_0, g_0_xzzzzzz_0_x_1, g_0_yzzzzzz_0_y_0, g_0_yzzzzzz_0_y_1, g_0_yzzzzzz_0_z_0, g_0_yzzzzzz_0_z_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_0_xyzzzzzz_0_x_0[i] = g_0_xzzzzzz_0_x_0[i] * pb_y + g_0_xzzzzzz_0_x_1[i] * wp_y[i];

        g_0_xyzzzzzz_0_y_0[i] = g_0_yzzzzzz_0_y_0[i] * pb_x + g_0_yzzzzzz_0_y_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_z_0[i] = g_0_yzzzzzz_0_z_0[i] * pb_x + g_0_yzzzzzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 105-108 components of targeted buffer : SLSP

    auto g_0_xzzzzzzz_0_x_0 = pbuffer.data(idx_eri_0_slsp + 105);

    auto g_0_xzzzzzzz_0_y_0 = pbuffer.data(idx_eri_0_slsp + 106);

    auto g_0_xzzzzzzz_0_z_0 = pbuffer.data(idx_eri_0_slsp + 107);

    #pragma omp simd aligned(g_0_xzzzzzzz_0_x_0, g_0_xzzzzzzz_0_y_0, g_0_xzzzzzzz_0_z_0, g_0_zzzzzzz_0_0_1, g_0_zzzzzzz_0_x_0, g_0_zzzzzzz_0_x_1, g_0_zzzzzzz_0_y_0, g_0_zzzzzzz_0_y_1, g_0_zzzzzzz_0_z_0, g_0_zzzzzzz_0_z_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzzzzzzz_0_x_0[i] = g_0_zzzzzzz_0_0_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_x_0[i] * pb_x + g_0_zzzzzzz_0_x_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_y_0[i] = g_0_zzzzzzz_0_y_0[i] * pb_x + g_0_zzzzzzz_0_y_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_z_0[i] = g_0_zzzzzzz_0_z_0[i] * pb_x + g_0_zzzzzzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 108-111 components of targeted buffer : SLSP

    auto g_0_yyyyyyyy_0_x_0 = pbuffer.data(idx_eri_0_slsp + 108);

    auto g_0_yyyyyyyy_0_y_0 = pbuffer.data(idx_eri_0_slsp + 109);

    auto g_0_yyyyyyyy_0_z_0 = pbuffer.data(idx_eri_0_slsp + 110);

    #pragma omp simd aligned(g_0_yyyyyy_0_x_0, g_0_yyyyyy_0_x_1, g_0_yyyyyy_0_y_0, g_0_yyyyyy_0_y_1, g_0_yyyyyy_0_z_0, g_0_yyyyyy_0_z_1, g_0_yyyyyyy_0_0_1, g_0_yyyyyyy_0_x_0, g_0_yyyyyyy_0_x_1, g_0_yyyyyyy_0_y_0, g_0_yyyyyyy_0_y_1, g_0_yyyyyyy_0_z_0, g_0_yyyyyyy_0_z_1, g_0_yyyyyyyy_0_x_0, g_0_yyyyyyyy_0_y_0, g_0_yyyyyyyy_0_z_0, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyyyyy_0_x_0[i] = 7.0 * g_0_yyyyyy_0_x_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_x_1[i] * fti_ab_0 + g_0_yyyyyyy_0_x_0[i] * pb_y + g_0_yyyyyyy_0_x_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_y_0[i] = 7.0 * g_0_yyyyyy_0_y_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_y_1[i] * fti_ab_0 + g_0_yyyyyyy_0_0_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_y_0[i] * pb_y + g_0_yyyyyyy_0_y_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_z_0[i] = 7.0 * g_0_yyyyyy_0_z_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_z_1[i] * fti_ab_0 + g_0_yyyyyyy_0_z_0[i] * pb_y + g_0_yyyyyyy_0_z_1[i] * wp_y[i];
    }

    /// Set up 111-114 components of targeted buffer : SLSP

    auto g_0_yyyyyyyz_0_x_0 = pbuffer.data(idx_eri_0_slsp + 111);

    auto g_0_yyyyyyyz_0_y_0 = pbuffer.data(idx_eri_0_slsp + 112);

    auto g_0_yyyyyyyz_0_z_0 = pbuffer.data(idx_eri_0_slsp + 113);

    #pragma omp simd aligned(g_0_yyyyyyy_0_0_1, g_0_yyyyyyy_0_x_0, g_0_yyyyyyy_0_x_1, g_0_yyyyyyy_0_y_0, g_0_yyyyyyy_0_y_1, g_0_yyyyyyy_0_z_0, g_0_yyyyyyy_0_z_1, g_0_yyyyyyyz_0_x_0, g_0_yyyyyyyz_0_y_0, g_0_yyyyyyyz_0_z_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyyyyyz_0_x_0[i] = g_0_yyyyyyy_0_x_0[i] * pb_z + g_0_yyyyyyy_0_x_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_y_0[i] = g_0_yyyyyyy_0_y_0[i] * pb_z + g_0_yyyyyyy_0_y_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_z_0[i] = g_0_yyyyyyy_0_0_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_z_0[i] * pb_z + g_0_yyyyyyy_0_z_1[i] * wp_z[i];
    }

    /// Set up 114-117 components of targeted buffer : SLSP

    auto g_0_yyyyyyzz_0_x_0 = pbuffer.data(idx_eri_0_slsp + 114);

    auto g_0_yyyyyyzz_0_y_0 = pbuffer.data(idx_eri_0_slsp + 115);

    auto g_0_yyyyyyzz_0_z_0 = pbuffer.data(idx_eri_0_slsp + 116);

    #pragma omp simd aligned(g_0_yyyyyy_0_y_0, g_0_yyyyyy_0_y_1, g_0_yyyyyyz_0_y_0, g_0_yyyyyyz_0_y_1, g_0_yyyyyyzz_0_x_0, g_0_yyyyyyzz_0_y_0, g_0_yyyyyyzz_0_z_0, g_0_yyyyyzz_0_x_0, g_0_yyyyyzz_0_x_1, g_0_yyyyyzz_0_z_0, g_0_yyyyyzz_0_z_1, g_0_yyyyzz_0_x_0, g_0_yyyyzz_0_x_1, g_0_yyyyzz_0_z_0, g_0_yyyyzz_0_z_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyyyyzz_0_x_0[i] = 5.0 * g_0_yyyyzz_0_x_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_x_1[i] * fti_ab_0 + g_0_yyyyyzz_0_x_0[i] * pb_y + g_0_yyyyyzz_0_x_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_y_0[i] = g_0_yyyyyy_0_y_0[i] * fi_ab_0 - g_0_yyyyyy_0_y_1[i] * fti_ab_0 + g_0_yyyyyyz_0_y_0[i] * pb_z + g_0_yyyyyyz_0_y_1[i] * wp_z[i];

        g_0_yyyyyyzz_0_z_0[i] = 5.0 * g_0_yyyyzz_0_z_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_z_1[i] * fti_ab_0 + g_0_yyyyyzz_0_z_0[i] * pb_y + g_0_yyyyyzz_0_z_1[i] * wp_y[i];
    }

    /// Set up 117-120 components of targeted buffer : SLSP

    auto g_0_yyyyyzzz_0_x_0 = pbuffer.data(idx_eri_0_slsp + 117);

    auto g_0_yyyyyzzz_0_y_0 = pbuffer.data(idx_eri_0_slsp + 118);

    auto g_0_yyyyyzzz_0_z_0 = pbuffer.data(idx_eri_0_slsp + 119);

    #pragma omp simd aligned(g_0_yyyyyz_0_y_0, g_0_yyyyyz_0_y_1, g_0_yyyyyzz_0_y_0, g_0_yyyyyzz_0_y_1, g_0_yyyyyzzz_0_x_0, g_0_yyyyyzzz_0_y_0, g_0_yyyyyzzz_0_z_0, g_0_yyyyzzz_0_x_0, g_0_yyyyzzz_0_x_1, g_0_yyyyzzz_0_z_0, g_0_yyyyzzz_0_z_1, g_0_yyyzzz_0_x_0, g_0_yyyzzz_0_x_1, g_0_yyyzzz_0_z_0, g_0_yyyzzz_0_z_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyyyzzz_0_x_0[i] = 4.0 * g_0_yyyzzz_0_x_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_x_1[i] * fti_ab_0 + g_0_yyyyzzz_0_x_0[i] * pb_y + g_0_yyyyzzz_0_x_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_y_0[i] = 2.0 * g_0_yyyyyz_0_y_0[i] * fi_ab_0 - 2.0 * g_0_yyyyyz_0_y_1[i] * fti_ab_0 + g_0_yyyyyzz_0_y_0[i] * pb_z + g_0_yyyyyzz_0_y_1[i] * wp_z[i];

        g_0_yyyyyzzz_0_z_0[i] = 4.0 * g_0_yyyzzz_0_z_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_z_1[i] * fti_ab_0 + g_0_yyyyzzz_0_z_0[i] * pb_y + g_0_yyyyzzz_0_z_1[i] * wp_y[i];
    }

    /// Set up 120-123 components of targeted buffer : SLSP

    auto g_0_yyyyzzzz_0_x_0 = pbuffer.data(idx_eri_0_slsp + 120);

    auto g_0_yyyyzzzz_0_y_0 = pbuffer.data(idx_eri_0_slsp + 121);

    auto g_0_yyyyzzzz_0_z_0 = pbuffer.data(idx_eri_0_slsp + 122);

    #pragma omp simd aligned(g_0_yyyyzz_0_y_0, g_0_yyyyzz_0_y_1, g_0_yyyyzzz_0_y_0, g_0_yyyyzzz_0_y_1, g_0_yyyyzzzz_0_x_0, g_0_yyyyzzzz_0_y_0, g_0_yyyyzzzz_0_z_0, g_0_yyyzzzz_0_x_0, g_0_yyyzzzz_0_x_1, g_0_yyyzzzz_0_z_0, g_0_yyyzzzz_0_z_1, g_0_yyzzzz_0_x_0, g_0_yyzzzz_0_x_1, g_0_yyzzzz_0_z_0, g_0_yyzzzz_0_z_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyyzzzz_0_x_0[i] = 3.0 * g_0_yyzzzz_0_x_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_x_1[i] * fti_ab_0 + g_0_yyyzzzz_0_x_0[i] * pb_y + g_0_yyyzzzz_0_x_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_y_0[i] = 3.0 * g_0_yyyyzz_0_y_0[i] * fi_ab_0 - 3.0 * g_0_yyyyzz_0_y_1[i] * fti_ab_0 + g_0_yyyyzzz_0_y_0[i] * pb_z + g_0_yyyyzzz_0_y_1[i] * wp_z[i];

        g_0_yyyyzzzz_0_z_0[i] = 3.0 * g_0_yyzzzz_0_z_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_z_1[i] * fti_ab_0 + g_0_yyyzzzz_0_z_0[i] * pb_y + g_0_yyyzzzz_0_z_1[i] * wp_y[i];
    }

    /// Set up 123-126 components of targeted buffer : SLSP

    auto g_0_yyyzzzzz_0_x_0 = pbuffer.data(idx_eri_0_slsp + 123);

    auto g_0_yyyzzzzz_0_y_0 = pbuffer.data(idx_eri_0_slsp + 124);

    auto g_0_yyyzzzzz_0_z_0 = pbuffer.data(idx_eri_0_slsp + 125);

    #pragma omp simd aligned(g_0_yyyzzz_0_y_0, g_0_yyyzzz_0_y_1, g_0_yyyzzzz_0_y_0, g_0_yyyzzzz_0_y_1, g_0_yyyzzzzz_0_x_0, g_0_yyyzzzzz_0_y_0, g_0_yyyzzzzz_0_z_0, g_0_yyzzzzz_0_x_0, g_0_yyzzzzz_0_x_1, g_0_yyzzzzz_0_z_0, g_0_yyzzzzz_0_z_1, g_0_yzzzzz_0_x_0, g_0_yzzzzz_0_x_1, g_0_yzzzzz_0_z_0, g_0_yzzzzz_0_z_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyzzzzz_0_x_0[i] = 2.0 * g_0_yzzzzz_0_x_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_x_1[i] * fti_ab_0 + g_0_yyzzzzz_0_x_0[i] * pb_y + g_0_yyzzzzz_0_x_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_y_0[i] = 4.0 * g_0_yyyzzz_0_y_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_y_1[i] * fti_ab_0 + g_0_yyyzzzz_0_y_0[i] * pb_z + g_0_yyyzzzz_0_y_1[i] * wp_z[i];

        g_0_yyyzzzzz_0_z_0[i] = 2.0 * g_0_yzzzzz_0_z_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_z_1[i] * fti_ab_0 + g_0_yyzzzzz_0_z_0[i] * pb_y + g_0_yyzzzzz_0_z_1[i] * wp_y[i];
    }

    /// Set up 126-129 components of targeted buffer : SLSP

    auto g_0_yyzzzzzz_0_x_0 = pbuffer.data(idx_eri_0_slsp + 126);

    auto g_0_yyzzzzzz_0_y_0 = pbuffer.data(idx_eri_0_slsp + 127);

    auto g_0_yyzzzzzz_0_z_0 = pbuffer.data(idx_eri_0_slsp + 128);

    #pragma omp simd aligned(g_0_yyzzzz_0_y_0, g_0_yyzzzz_0_y_1, g_0_yyzzzzz_0_y_0, g_0_yyzzzzz_0_y_1, g_0_yyzzzzzz_0_x_0, g_0_yyzzzzzz_0_y_0, g_0_yyzzzzzz_0_z_0, g_0_yzzzzzz_0_x_0, g_0_yzzzzzz_0_x_1, g_0_yzzzzzz_0_z_0, g_0_yzzzzzz_0_z_1, g_0_zzzzzz_0_x_0, g_0_zzzzzz_0_x_1, g_0_zzzzzz_0_z_0, g_0_zzzzzz_0_z_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyzzzzzz_0_x_0[i] = g_0_zzzzzz_0_x_0[i] * fi_ab_0 - g_0_zzzzzz_0_x_1[i] * fti_ab_0 + g_0_yzzzzzz_0_x_0[i] * pb_y + g_0_yzzzzzz_0_x_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_y_0[i] = 5.0 * g_0_yyzzzz_0_y_0[i] * fi_ab_0 - 5.0 * g_0_yyzzzz_0_y_1[i] * fti_ab_0 + g_0_yyzzzzz_0_y_0[i] * pb_z + g_0_yyzzzzz_0_y_1[i] * wp_z[i];

        g_0_yyzzzzzz_0_z_0[i] = g_0_zzzzzz_0_z_0[i] * fi_ab_0 - g_0_zzzzzz_0_z_1[i] * fti_ab_0 + g_0_yzzzzzz_0_z_0[i] * pb_y + g_0_yzzzzzz_0_z_1[i] * wp_y[i];
    }

    /// Set up 129-132 components of targeted buffer : SLSP

    auto g_0_yzzzzzzz_0_x_0 = pbuffer.data(idx_eri_0_slsp + 129);

    auto g_0_yzzzzzzz_0_y_0 = pbuffer.data(idx_eri_0_slsp + 130);

    auto g_0_yzzzzzzz_0_z_0 = pbuffer.data(idx_eri_0_slsp + 131);

    #pragma omp simd aligned(g_0_yzzzzzzz_0_x_0, g_0_yzzzzzzz_0_y_0, g_0_yzzzzzzz_0_z_0, g_0_zzzzzzz_0_0_1, g_0_zzzzzzz_0_x_0, g_0_zzzzzzz_0_x_1, g_0_zzzzzzz_0_y_0, g_0_zzzzzzz_0_y_1, g_0_zzzzzzz_0_z_0, g_0_zzzzzzz_0_z_1, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzzzzzzz_0_x_0[i] = g_0_zzzzzzz_0_x_0[i] * pb_y + g_0_zzzzzzz_0_x_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_y_0[i] = g_0_zzzzzzz_0_0_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_y_0[i] * pb_y + g_0_zzzzzzz_0_y_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_z_0[i] = g_0_zzzzzzz_0_z_0[i] * pb_y + g_0_zzzzzzz_0_z_1[i] * wp_y[i];
    }

    /// Set up 132-135 components of targeted buffer : SLSP

    auto g_0_zzzzzzzz_0_x_0 = pbuffer.data(idx_eri_0_slsp + 132);

    auto g_0_zzzzzzzz_0_y_0 = pbuffer.data(idx_eri_0_slsp + 133);

    auto g_0_zzzzzzzz_0_z_0 = pbuffer.data(idx_eri_0_slsp + 134);

    #pragma omp simd aligned(g_0_zzzzzz_0_x_0, g_0_zzzzzz_0_x_1, g_0_zzzzzz_0_y_0, g_0_zzzzzz_0_y_1, g_0_zzzzzz_0_z_0, g_0_zzzzzz_0_z_1, g_0_zzzzzzz_0_0_1, g_0_zzzzzzz_0_x_0, g_0_zzzzzzz_0_x_1, g_0_zzzzzzz_0_y_0, g_0_zzzzzzz_0_y_1, g_0_zzzzzzz_0_z_0, g_0_zzzzzzz_0_z_1, g_0_zzzzzzzz_0_x_0, g_0_zzzzzzzz_0_y_0, g_0_zzzzzzzz_0_z_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzzzzzzz_0_x_0[i] = 7.0 * g_0_zzzzzz_0_x_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_x_1[i] * fti_ab_0 + g_0_zzzzzzz_0_x_0[i] * pb_z + g_0_zzzzzzz_0_x_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_y_0[i] = 7.0 * g_0_zzzzzz_0_y_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_y_1[i] * fti_ab_0 + g_0_zzzzzzz_0_y_0[i] * pb_z + g_0_zzzzzzz_0_y_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_z_0[i] = 7.0 * g_0_zzzzzz_0_z_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_z_1[i] * fti_ab_0 + g_0_zzzzzzz_0_0_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_z_0[i] * pb_z + g_0_zzzzzzz_0_z_1[i] * wp_z[i];
    }
}

} // erirec namespace


#include "ElectronRepulsionPrimRecSFSH.hpp"

namespace erirec {  // erirec namespace

auto
comp_prim_electron_repulsion_sfsh(CSimdArray<double>&   pbuffer,
                                  const size_t          idx_eri_0_sfsh,
                                  size_t                idx_eri_0_spsh,
                                  size_t                idx_eri_1_spsh,
                                  size_t                idx_eri_1_sdsg,
                                  size_t                idx_eri_0_sdsh,
                                  size_t                idx_eri_1_sdsh,
                                  CSimdArray<double>&   factors,
                                  const size_t          idx_wp,
                                  const TPoint<double>& r_pb,
                                  const double          a_exp,
                                  const double          b_exp) -> void
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

    /// Set up components of auxilary buffer : SPSH

    auto g_0_x_0_xxxxx_0 = pbuffer.data(idx_eri_0_spsh);

    auto g_0_x_0_xxxxy_0 = pbuffer.data(idx_eri_0_spsh + 1);

    auto g_0_x_0_xxxxz_0 = pbuffer.data(idx_eri_0_spsh + 2);

    auto g_0_x_0_xxxyy_0 = pbuffer.data(idx_eri_0_spsh + 3);

    auto g_0_x_0_xxxyz_0 = pbuffer.data(idx_eri_0_spsh + 4);

    auto g_0_x_0_xxxzz_0 = pbuffer.data(idx_eri_0_spsh + 5);

    auto g_0_x_0_xxyyy_0 = pbuffer.data(idx_eri_0_spsh + 6);

    auto g_0_x_0_xxyyz_0 = pbuffer.data(idx_eri_0_spsh + 7);

    auto g_0_x_0_xxyzz_0 = pbuffer.data(idx_eri_0_spsh + 8);

    auto g_0_x_0_xxzzz_0 = pbuffer.data(idx_eri_0_spsh + 9);

    auto g_0_x_0_xyyyy_0 = pbuffer.data(idx_eri_0_spsh + 10);

    auto g_0_x_0_xyyyz_0 = pbuffer.data(idx_eri_0_spsh + 11);

    auto g_0_x_0_xyyzz_0 = pbuffer.data(idx_eri_0_spsh + 12);

    auto g_0_x_0_xyzzz_0 = pbuffer.data(idx_eri_0_spsh + 13);

    auto g_0_x_0_xzzzz_0 = pbuffer.data(idx_eri_0_spsh + 14);

    auto g_0_x_0_yyyyy_0 = pbuffer.data(idx_eri_0_spsh + 15);

    auto g_0_x_0_yyyyz_0 = pbuffer.data(idx_eri_0_spsh + 16);

    auto g_0_x_0_yyyzz_0 = pbuffer.data(idx_eri_0_spsh + 17);

    auto g_0_x_0_yyzzz_0 = pbuffer.data(idx_eri_0_spsh + 18);

    auto g_0_x_0_yzzzz_0 = pbuffer.data(idx_eri_0_spsh + 19);

    auto g_0_x_0_zzzzz_0 = pbuffer.data(idx_eri_0_spsh + 20);

    auto g_0_y_0_xxxxx_0 = pbuffer.data(idx_eri_0_spsh + 21);

    auto g_0_y_0_xxxxy_0 = pbuffer.data(idx_eri_0_spsh + 22);

    auto g_0_y_0_xxxxz_0 = pbuffer.data(idx_eri_0_spsh + 23);

    auto g_0_y_0_xxxyy_0 = pbuffer.data(idx_eri_0_spsh + 24);

    auto g_0_y_0_xxxyz_0 = pbuffer.data(idx_eri_0_spsh + 25);

    auto g_0_y_0_xxxzz_0 = pbuffer.data(idx_eri_0_spsh + 26);

    auto g_0_y_0_xxyyy_0 = pbuffer.data(idx_eri_0_spsh + 27);

    auto g_0_y_0_xxyyz_0 = pbuffer.data(idx_eri_0_spsh + 28);

    auto g_0_y_0_xxyzz_0 = pbuffer.data(idx_eri_0_spsh + 29);

    auto g_0_y_0_xxzzz_0 = pbuffer.data(idx_eri_0_spsh + 30);

    auto g_0_y_0_xyyyy_0 = pbuffer.data(idx_eri_0_spsh + 31);

    auto g_0_y_0_xyyyz_0 = pbuffer.data(idx_eri_0_spsh + 32);

    auto g_0_y_0_xyyzz_0 = pbuffer.data(idx_eri_0_spsh + 33);

    auto g_0_y_0_xyzzz_0 = pbuffer.data(idx_eri_0_spsh + 34);

    auto g_0_y_0_xzzzz_0 = pbuffer.data(idx_eri_0_spsh + 35);

    auto g_0_y_0_yyyyy_0 = pbuffer.data(idx_eri_0_spsh + 36);

    auto g_0_y_0_yyyyz_0 = pbuffer.data(idx_eri_0_spsh + 37);

    auto g_0_y_0_yyyzz_0 = pbuffer.data(idx_eri_0_spsh + 38);

    auto g_0_y_0_yyzzz_0 = pbuffer.data(idx_eri_0_spsh + 39);

    auto g_0_y_0_yzzzz_0 = pbuffer.data(idx_eri_0_spsh + 40);

    auto g_0_y_0_zzzzz_0 = pbuffer.data(idx_eri_0_spsh + 41);

    auto g_0_z_0_xxxxx_0 = pbuffer.data(idx_eri_0_spsh + 42);

    auto g_0_z_0_xxxxy_0 = pbuffer.data(idx_eri_0_spsh + 43);

    auto g_0_z_0_xxxxz_0 = pbuffer.data(idx_eri_0_spsh + 44);

    auto g_0_z_0_xxxyy_0 = pbuffer.data(idx_eri_0_spsh + 45);

    auto g_0_z_0_xxxyz_0 = pbuffer.data(idx_eri_0_spsh + 46);

    auto g_0_z_0_xxxzz_0 = pbuffer.data(idx_eri_0_spsh + 47);

    auto g_0_z_0_xxyyy_0 = pbuffer.data(idx_eri_0_spsh + 48);

    auto g_0_z_0_xxyyz_0 = pbuffer.data(idx_eri_0_spsh + 49);

    auto g_0_z_0_xxyzz_0 = pbuffer.data(idx_eri_0_spsh + 50);

    auto g_0_z_0_xxzzz_0 = pbuffer.data(idx_eri_0_spsh + 51);

    auto g_0_z_0_xyyyy_0 = pbuffer.data(idx_eri_0_spsh + 52);

    auto g_0_z_0_xyyyz_0 = pbuffer.data(idx_eri_0_spsh + 53);

    auto g_0_z_0_xyyzz_0 = pbuffer.data(idx_eri_0_spsh + 54);

    auto g_0_z_0_xyzzz_0 = pbuffer.data(idx_eri_0_spsh + 55);

    auto g_0_z_0_xzzzz_0 = pbuffer.data(idx_eri_0_spsh + 56);

    auto g_0_z_0_yyyyy_0 = pbuffer.data(idx_eri_0_spsh + 57);

    auto g_0_z_0_yyyyz_0 = pbuffer.data(idx_eri_0_spsh + 58);

    auto g_0_z_0_yyyzz_0 = pbuffer.data(idx_eri_0_spsh + 59);

    auto g_0_z_0_yyzzz_0 = pbuffer.data(idx_eri_0_spsh + 60);

    auto g_0_z_0_yzzzz_0 = pbuffer.data(idx_eri_0_spsh + 61);

    auto g_0_z_0_zzzzz_0 = pbuffer.data(idx_eri_0_spsh + 62);

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

    /// Set up components of auxilary buffer : SDSG

    auto g_0_xx_0_xxxx_1 = pbuffer.data(idx_eri_1_sdsg);

    auto g_0_xx_0_xxxy_1 = pbuffer.data(idx_eri_1_sdsg + 1);

    auto g_0_xx_0_xxxz_1 = pbuffer.data(idx_eri_1_sdsg + 2);

    auto g_0_xx_0_xxyy_1 = pbuffer.data(idx_eri_1_sdsg + 3);

    auto g_0_xx_0_xxyz_1 = pbuffer.data(idx_eri_1_sdsg + 4);

    auto g_0_xx_0_xxzz_1 = pbuffer.data(idx_eri_1_sdsg + 5);

    auto g_0_xx_0_xyyy_1 = pbuffer.data(idx_eri_1_sdsg + 6);

    auto g_0_xx_0_xyyz_1 = pbuffer.data(idx_eri_1_sdsg + 7);

    auto g_0_xx_0_xyzz_1 = pbuffer.data(idx_eri_1_sdsg + 8);

    auto g_0_xx_0_xzzz_1 = pbuffer.data(idx_eri_1_sdsg + 9);

    auto g_0_xx_0_yyyy_1 = pbuffer.data(idx_eri_1_sdsg + 10);

    auto g_0_xx_0_yyyz_1 = pbuffer.data(idx_eri_1_sdsg + 11);

    auto g_0_xx_0_yyzz_1 = pbuffer.data(idx_eri_1_sdsg + 12);

    auto g_0_xx_0_yzzz_1 = pbuffer.data(idx_eri_1_sdsg + 13);

    auto g_0_xx_0_zzzz_1 = pbuffer.data(idx_eri_1_sdsg + 14);

    auto g_0_yy_0_xxxx_1 = pbuffer.data(idx_eri_1_sdsg + 45);

    auto g_0_yy_0_xxxy_1 = pbuffer.data(idx_eri_1_sdsg + 46);

    auto g_0_yy_0_xxxz_1 = pbuffer.data(idx_eri_1_sdsg + 47);

    auto g_0_yy_0_xxyy_1 = pbuffer.data(idx_eri_1_sdsg + 48);

    auto g_0_yy_0_xxyz_1 = pbuffer.data(idx_eri_1_sdsg + 49);

    auto g_0_yy_0_xxzz_1 = pbuffer.data(idx_eri_1_sdsg + 50);

    auto g_0_yy_0_xyyy_1 = pbuffer.data(idx_eri_1_sdsg + 51);

    auto g_0_yy_0_xyyz_1 = pbuffer.data(idx_eri_1_sdsg + 52);

    auto g_0_yy_0_xyzz_1 = pbuffer.data(idx_eri_1_sdsg + 53);

    auto g_0_yy_0_xzzz_1 = pbuffer.data(idx_eri_1_sdsg + 54);

    auto g_0_yy_0_yyyy_1 = pbuffer.data(idx_eri_1_sdsg + 55);

    auto g_0_yy_0_yyyz_1 = pbuffer.data(idx_eri_1_sdsg + 56);

    auto g_0_yy_0_yyzz_1 = pbuffer.data(idx_eri_1_sdsg + 57);

    auto g_0_yy_0_yzzz_1 = pbuffer.data(idx_eri_1_sdsg + 58);

    auto g_0_yy_0_zzzz_1 = pbuffer.data(idx_eri_1_sdsg + 59);

    auto g_0_yz_0_xxyz_1 = pbuffer.data(idx_eri_1_sdsg + 64);

    auto g_0_yz_0_xyyz_1 = pbuffer.data(idx_eri_1_sdsg + 67);

    auto g_0_yz_0_xyzz_1 = pbuffer.data(idx_eri_1_sdsg + 68);

    auto g_0_yz_0_yyyz_1 = pbuffer.data(idx_eri_1_sdsg + 71);

    auto g_0_yz_0_yyzz_1 = pbuffer.data(idx_eri_1_sdsg + 72);

    auto g_0_yz_0_yzzz_1 = pbuffer.data(idx_eri_1_sdsg + 73);

    auto g_0_zz_0_xxxx_1 = pbuffer.data(idx_eri_1_sdsg + 75);

    auto g_0_zz_0_xxxy_1 = pbuffer.data(idx_eri_1_sdsg + 76);

    auto g_0_zz_0_xxxz_1 = pbuffer.data(idx_eri_1_sdsg + 77);

    auto g_0_zz_0_xxyy_1 = pbuffer.data(idx_eri_1_sdsg + 78);

    auto g_0_zz_0_xxyz_1 = pbuffer.data(idx_eri_1_sdsg + 79);

    auto g_0_zz_0_xxzz_1 = pbuffer.data(idx_eri_1_sdsg + 80);

    auto g_0_zz_0_xyyy_1 = pbuffer.data(idx_eri_1_sdsg + 81);

    auto g_0_zz_0_xyyz_1 = pbuffer.data(idx_eri_1_sdsg + 82);

    auto g_0_zz_0_xyzz_1 = pbuffer.data(idx_eri_1_sdsg + 83);

    auto g_0_zz_0_xzzz_1 = pbuffer.data(idx_eri_1_sdsg + 84);

    auto g_0_zz_0_yyyy_1 = pbuffer.data(idx_eri_1_sdsg + 85);

    auto g_0_zz_0_yyyz_1 = pbuffer.data(idx_eri_1_sdsg + 86);

    auto g_0_zz_0_yyzz_1 = pbuffer.data(idx_eri_1_sdsg + 87);

    auto g_0_zz_0_yzzz_1 = pbuffer.data(idx_eri_1_sdsg + 88);

    auto g_0_zz_0_zzzz_1 = pbuffer.data(idx_eri_1_sdsg + 89);

    /// Set up components of auxilary buffer : SDSH

    auto g_0_xx_0_xxxxx_0 = pbuffer.data(idx_eri_0_sdsh);

    auto g_0_xx_0_xxxxy_0 = pbuffer.data(idx_eri_0_sdsh + 1);

    auto g_0_xx_0_xxxxz_0 = pbuffer.data(idx_eri_0_sdsh + 2);

    auto g_0_xx_0_xxxyy_0 = pbuffer.data(idx_eri_0_sdsh + 3);

    auto g_0_xx_0_xxxyz_0 = pbuffer.data(idx_eri_0_sdsh + 4);

    auto g_0_xx_0_xxxzz_0 = pbuffer.data(idx_eri_0_sdsh + 5);

    auto g_0_xx_0_xxyyy_0 = pbuffer.data(idx_eri_0_sdsh + 6);

    auto g_0_xx_0_xxyyz_0 = pbuffer.data(idx_eri_0_sdsh + 7);

    auto g_0_xx_0_xxyzz_0 = pbuffer.data(idx_eri_0_sdsh + 8);

    auto g_0_xx_0_xxzzz_0 = pbuffer.data(idx_eri_0_sdsh + 9);

    auto g_0_xx_0_xyyyy_0 = pbuffer.data(idx_eri_0_sdsh + 10);

    auto g_0_xx_0_xyyyz_0 = pbuffer.data(idx_eri_0_sdsh + 11);

    auto g_0_xx_0_xyyzz_0 = pbuffer.data(idx_eri_0_sdsh + 12);

    auto g_0_xx_0_xyzzz_0 = pbuffer.data(idx_eri_0_sdsh + 13);

    auto g_0_xx_0_xzzzz_0 = pbuffer.data(idx_eri_0_sdsh + 14);

    auto g_0_xx_0_yyyyy_0 = pbuffer.data(idx_eri_0_sdsh + 15);

    auto g_0_xx_0_yyyyz_0 = pbuffer.data(idx_eri_0_sdsh + 16);

    auto g_0_xx_0_yyyzz_0 = pbuffer.data(idx_eri_0_sdsh + 17);

    auto g_0_xx_0_yyzzz_0 = pbuffer.data(idx_eri_0_sdsh + 18);

    auto g_0_xx_0_yzzzz_0 = pbuffer.data(idx_eri_0_sdsh + 19);

    auto g_0_xx_0_zzzzz_0 = pbuffer.data(idx_eri_0_sdsh + 20);

    auto g_0_xy_0_xxxxy_0 = pbuffer.data(idx_eri_0_sdsh + 22);

    auto g_0_xy_0_xxxyy_0 = pbuffer.data(idx_eri_0_sdsh + 24);

    auto g_0_xy_0_xxyyy_0 = pbuffer.data(idx_eri_0_sdsh + 27);

    auto g_0_xy_0_xyyyy_0 = pbuffer.data(idx_eri_0_sdsh + 31);

    auto g_0_xz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sdsh + 42);

    auto g_0_xz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sdsh + 44);

    auto g_0_xz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sdsh + 47);

    auto g_0_xz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sdsh + 51);

    auto g_0_xz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sdsh + 56);

    auto g_0_yy_0_xxxxx_0 = pbuffer.data(idx_eri_0_sdsh + 63);

    auto g_0_yy_0_xxxxy_0 = pbuffer.data(idx_eri_0_sdsh + 64);

    auto g_0_yy_0_xxxxz_0 = pbuffer.data(idx_eri_0_sdsh + 65);

    auto g_0_yy_0_xxxyy_0 = pbuffer.data(idx_eri_0_sdsh + 66);

    auto g_0_yy_0_xxxyz_0 = pbuffer.data(idx_eri_0_sdsh + 67);

    auto g_0_yy_0_xxxzz_0 = pbuffer.data(idx_eri_0_sdsh + 68);

    auto g_0_yy_0_xxyyy_0 = pbuffer.data(idx_eri_0_sdsh + 69);

    auto g_0_yy_0_xxyyz_0 = pbuffer.data(idx_eri_0_sdsh + 70);

    auto g_0_yy_0_xxyzz_0 = pbuffer.data(idx_eri_0_sdsh + 71);

    auto g_0_yy_0_xxzzz_0 = pbuffer.data(idx_eri_0_sdsh + 72);

    auto g_0_yy_0_xyyyy_0 = pbuffer.data(idx_eri_0_sdsh + 73);

    auto g_0_yy_0_xyyyz_0 = pbuffer.data(idx_eri_0_sdsh + 74);

    auto g_0_yy_0_xyyzz_0 = pbuffer.data(idx_eri_0_sdsh + 75);

    auto g_0_yy_0_xyzzz_0 = pbuffer.data(idx_eri_0_sdsh + 76);

    auto g_0_yy_0_xzzzz_0 = pbuffer.data(idx_eri_0_sdsh + 77);

    auto g_0_yy_0_yyyyy_0 = pbuffer.data(idx_eri_0_sdsh + 78);

    auto g_0_yy_0_yyyyz_0 = pbuffer.data(idx_eri_0_sdsh + 79);

    auto g_0_yy_0_yyyzz_0 = pbuffer.data(idx_eri_0_sdsh + 80);

    auto g_0_yy_0_yyzzz_0 = pbuffer.data(idx_eri_0_sdsh + 81);

    auto g_0_yy_0_yzzzz_0 = pbuffer.data(idx_eri_0_sdsh + 82);

    auto g_0_yy_0_zzzzz_0 = pbuffer.data(idx_eri_0_sdsh + 83);

    auto g_0_yz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sdsh + 88);

    auto g_0_yz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sdsh + 91);

    auto g_0_yz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sdsh + 92);

    auto g_0_yz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sdsh + 95);

    auto g_0_yz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sdsh + 96);

    auto g_0_yz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sdsh + 97);

    auto g_0_yz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sdsh + 99);

    auto g_0_yz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sdsh + 100);

    auto g_0_yz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sdsh + 101);

    auto g_0_yz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sdsh + 102);

    auto g_0_yz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sdsh + 103);

    auto g_0_yz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sdsh + 104);

    auto g_0_zz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sdsh + 105);

    auto g_0_zz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sdsh + 106);

    auto g_0_zz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sdsh + 107);

    auto g_0_zz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sdsh + 108);

    auto g_0_zz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sdsh + 109);

    auto g_0_zz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sdsh + 110);

    auto g_0_zz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sdsh + 111);

    auto g_0_zz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sdsh + 112);

    auto g_0_zz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sdsh + 113);

    auto g_0_zz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sdsh + 114);

    auto g_0_zz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sdsh + 115);

    auto g_0_zz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sdsh + 116);

    auto g_0_zz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sdsh + 117);

    auto g_0_zz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sdsh + 118);

    auto g_0_zz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sdsh + 119);

    auto g_0_zz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sdsh + 120);

    auto g_0_zz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sdsh + 121);

    auto g_0_zz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sdsh + 122);

    auto g_0_zz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sdsh + 123);

    auto g_0_zz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sdsh + 124);

    auto g_0_zz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sdsh + 125);

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

    auto g_0_xy_0_xxxxy_1 = pbuffer.data(idx_eri_1_sdsh + 22);

    auto g_0_xy_0_xxxyy_1 = pbuffer.data(idx_eri_1_sdsh + 24);

    auto g_0_xy_0_xxyyy_1 = pbuffer.data(idx_eri_1_sdsh + 27);

    auto g_0_xy_0_xyyyy_1 = pbuffer.data(idx_eri_1_sdsh + 31);

    auto g_0_xz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sdsh + 42);

    auto g_0_xz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sdsh + 44);

    auto g_0_xz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sdsh + 47);

    auto g_0_xz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sdsh + 51);

    auto g_0_xz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sdsh + 56);

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

    auto g_0_yz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sdsh + 99);

    auto g_0_yz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sdsh + 100);

    auto g_0_yz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sdsh + 101);

    auto g_0_yz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sdsh + 102);

    auto g_0_yz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sdsh + 103);

    auto g_0_yz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sdsh + 104);

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

    /// Set up 0-21 components of targeted buffer : SFSH

    auto g_0_xxx_0_xxxxx_0 = pbuffer.data(idx_eri_0_sfsh);

    auto g_0_xxx_0_xxxxy_0 = pbuffer.data(idx_eri_0_sfsh + 1);

    auto g_0_xxx_0_xxxxz_0 = pbuffer.data(idx_eri_0_sfsh + 2);

    auto g_0_xxx_0_xxxyy_0 = pbuffer.data(idx_eri_0_sfsh + 3);

    auto g_0_xxx_0_xxxyz_0 = pbuffer.data(idx_eri_0_sfsh + 4);

    auto g_0_xxx_0_xxxzz_0 = pbuffer.data(idx_eri_0_sfsh + 5);

    auto g_0_xxx_0_xxyyy_0 = pbuffer.data(idx_eri_0_sfsh + 6);

    auto g_0_xxx_0_xxyyz_0 = pbuffer.data(idx_eri_0_sfsh + 7);

    auto g_0_xxx_0_xxyzz_0 = pbuffer.data(idx_eri_0_sfsh + 8);

    auto g_0_xxx_0_xxzzz_0 = pbuffer.data(idx_eri_0_sfsh + 9);

    auto g_0_xxx_0_xyyyy_0 = pbuffer.data(idx_eri_0_sfsh + 10);

    auto g_0_xxx_0_xyyyz_0 = pbuffer.data(idx_eri_0_sfsh + 11);

    auto g_0_xxx_0_xyyzz_0 = pbuffer.data(idx_eri_0_sfsh + 12);

    auto g_0_xxx_0_xyzzz_0 = pbuffer.data(idx_eri_0_sfsh + 13);

    auto g_0_xxx_0_xzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 14);

    auto g_0_xxx_0_yyyyy_0 = pbuffer.data(idx_eri_0_sfsh + 15);

    auto g_0_xxx_0_yyyyz_0 = pbuffer.data(idx_eri_0_sfsh + 16);

    auto g_0_xxx_0_yyyzz_0 = pbuffer.data(idx_eri_0_sfsh + 17);

    auto g_0_xxx_0_yyzzz_0 = pbuffer.data(idx_eri_0_sfsh + 18);

    auto g_0_xxx_0_yzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 19);

    auto g_0_xxx_0_zzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 20);

#pragma omp simd aligned(g_0_x_0_xxxxx_0,       \
                             g_0_x_0_xxxxx_1,   \
                             g_0_x_0_xxxxy_0,   \
                             g_0_x_0_xxxxy_1,   \
                             g_0_x_0_xxxxz_0,   \
                             g_0_x_0_xxxxz_1,   \
                             g_0_x_0_xxxyy_0,   \
                             g_0_x_0_xxxyy_1,   \
                             g_0_x_0_xxxyz_0,   \
                             g_0_x_0_xxxyz_1,   \
                             g_0_x_0_xxxzz_0,   \
                             g_0_x_0_xxxzz_1,   \
                             g_0_x_0_xxyyy_0,   \
                             g_0_x_0_xxyyy_1,   \
                             g_0_x_0_xxyyz_0,   \
                             g_0_x_0_xxyyz_1,   \
                             g_0_x_0_xxyzz_0,   \
                             g_0_x_0_xxyzz_1,   \
                             g_0_x_0_xxzzz_0,   \
                             g_0_x_0_xxzzz_1,   \
                             g_0_x_0_xyyyy_0,   \
                             g_0_x_0_xyyyy_1,   \
                             g_0_x_0_xyyyz_0,   \
                             g_0_x_0_xyyyz_1,   \
                             g_0_x_0_xyyzz_0,   \
                             g_0_x_0_xyyzz_1,   \
                             g_0_x_0_xyzzz_0,   \
                             g_0_x_0_xyzzz_1,   \
                             g_0_x_0_xzzzz_0,   \
                             g_0_x_0_xzzzz_1,   \
                             g_0_x_0_yyyyy_0,   \
                             g_0_x_0_yyyyy_1,   \
                             g_0_x_0_yyyyz_0,   \
                             g_0_x_0_yyyyz_1,   \
                             g_0_x_0_yyyzz_0,   \
                             g_0_x_0_yyyzz_1,   \
                             g_0_x_0_yyzzz_0,   \
                             g_0_x_0_yyzzz_1,   \
                             g_0_x_0_yzzzz_0,   \
                             g_0_x_0_yzzzz_1,   \
                             g_0_x_0_zzzzz_0,   \
                             g_0_x_0_zzzzz_1,   \
                             g_0_xx_0_xxxx_1,   \
                             g_0_xx_0_xxxxx_0,  \
                             g_0_xx_0_xxxxx_1,  \
                             g_0_xx_0_xxxxy_0,  \
                             g_0_xx_0_xxxxy_1,  \
                             g_0_xx_0_xxxxz_0,  \
                             g_0_xx_0_xxxxz_1,  \
                             g_0_xx_0_xxxy_1,   \
                             g_0_xx_0_xxxyy_0,  \
                             g_0_xx_0_xxxyy_1,  \
                             g_0_xx_0_xxxyz_0,  \
                             g_0_xx_0_xxxyz_1,  \
                             g_0_xx_0_xxxz_1,   \
                             g_0_xx_0_xxxzz_0,  \
                             g_0_xx_0_xxxzz_1,  \
                             g_0_xx_0_xxyy_1,   \
                             g_0_xx_0_xxyyy_0,  \
                             g_0_xx_0_xxyyy_1,  \
                             g_0_xx_0_xxyyz_0,  \
                             g_0_xx_0_xxyyz_1,  \
                             g_0_xx_0_xxyz_1,   \
                             g_0_xx_0_xxyzz_0,  \
                             g_0_xx_0_xxyzz_1,  \
                             g_0_xx_0_xxzz_1,   \
                             g_0_xx_0_xxzzz_0,  \
                             g_0_xx_0_xxzzz_1,  \
                             g_0_xx_0_xyyy_1,   \
                             g_0_xx_0_xyyyy_0,  \
                             g_0_xx_0_xyyyy_1,  \
                             g_0_xx_0_xyyyz_0,  \
                             g_0_xx_0_xyyyz_1,  \
                             g_0_xx_0_xyyz_1,   \
                             g_0_xx_0_xyyzz_0,  \
                             g_0_xx_0_xyyzz_1,  \
                             g_0_xx_0_xyzz_1,   \
                             g_0_xx_0_xyzzz_0,  \
                             g_0_xx_0_xyzzz_1,  \
                             g_0_xx_0_xzzz_1,   \
                             g_0_xx_0_xzzzz_0,  \
                             g_0_xx_0_xzzzz_1,  \
                             g_0_xx_0_yyyy_1,   \
                             g_0_xx_0_yyyyy_0,  \
                             g_0_xx_0_yyyyy_1,  \
                             g_0_xx_0_yyyyz_0,  \
                             g_0_xx_0_yyyyz_1,  \
                             g_0_xx_0_yyyz_1,   \
                             g_0_xx_0_yyyzz_0,  \
                             g_0_xx_0_yyyzz_1,  \
                             g_0_xx_0_yyzz_1,   \
                             g_0_xx_0_yyzzz_0,  \
                             g_0_xx_0_yyzzz_1,  \
                             g_0_xx_0_yzzz_1,   \
                             g_0_xx_0_yzzzz_0,  \
                             g_0_xx_0_yzzzz_1,  \
                             g_0_xx_0_zzzz_1,   \
                             g_0_xx_0_zzzzz_0,  \
                             g_0_xx_0_zzzzz_1,  \
                             g_0_xxx_0_xxxxx_0, \
                             g_0_xxx_0_xxxxy_0, \
                             g_0_xxx_0_xxxxz_0, \
                             g_0_xxx_0_xxxyy_0, \
                             g_0_xxx_0_xxxyz_0, \
                             g_0_xxx_0_xxxzz_0, \
                             g_0_xxx_0_xxyyy_0, \
                             g_0_xxx_0_xxyyz_0, \
                             g_0_xxx_0_xxyzz_0, \
                             g_0_xxx_0_xxzzz_0, \
                             g_0_xxx_0_xyyyy_0, \
                             g_0_xxx_0_xyyyz_0, \
                             g_0_xxx_0_xyyzz_0, \
                             g_0_xxx_0_xyzzz_0, \
                             g_0_xxx_0_xzzzz_0, \
                             g_0_xxx_0_yyyyy_0, \
                             g_0_xxx_0_yyyyz_0, \
                             g_0_xxx_0_yyyzz_0, \
                             g_0_xxx_0_yyzzz_0, \
                             g_0_xxx_0_yzzzz_0, \
                             g_0_xxx_0_zzzzz_0, \
                             wp_x,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxx_0_xxxxx_0[i] = 2.0 * g_0_x_0_xxxxx_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxxx_1[i] * fti_ab_0 + 5.0 * g_0_xx_0_xxxx_1[i] * fi_abcd_0 +
                               g_0_xx_0_xxxxx_0[i] * pb_x + g_0_xx_0_xxxxx_1[i] * wp_x[i];

        g_0_xxx_0_xxxxy_0[i] = 2.0 * g_0_x_0_xxxxy_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxxy_1[i] * fti_ab_0 + 4.0 * g_0_xx_0_xxxy_1[i] * fi_abcd_0 +
                               g_0_xx_0_xxxxy_0[i] * pb_x + g_0_xx_0_xxxxy_1[i] * wp_x[i];

        g_0_xxx_0_xxxxz_0[i] = 2.0 * g_0_x_0_xxxxz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxxz_1[i] * fti_ab_0 + 4.0 * g_0_xx_0_xxxz_1[i] * fi_abcd_0 +
                               g_0_xx_0_xxxxz_0[i] * pb_x + g_0_xx_0_xxxxz_1[i] * wp_x[i];

        g_0_xxx_0_xxxyy_0[i] = 2.0 * g_0_x_0_xxxyy_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxyy_1[i] * fti_ab_0 + 3.0 * g_0_xx_0_xxyy_1[i] * fi_abcd_0 +
                               g_0_xx_0_xxxyy_0[i] * pb_x + g_0_xx_0_xxxyy_1[i] * wp_x[i];

        g_0_xxx_0_xxxyz_0[i] = 2.0 * g_0_x_0_xxxyz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xx_0_xxyz_1[i] * fi_abcd_0 +
                               g_0_xx_0_xxxyz_0[i] * pb_x + g_0_xx_0_xxxyz_1[i] * wp_x[i];

        g_0_xxx_0_xxxzz_0[i] = 2.0 * g_0_x_0_xxxzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxxzz_1[i] * fti_ab_0 + 3.0 * g_0_xx_0_xxzz_1[i] * fi_abcd_0 +
                               g_0_xx_0_xxxzz_0[i] * pb_x + g_0_xx_0_xxxzz_1[i] * wp_x[i];

        g_0_xxx_0_xxyyy_0[i] = 2.0 * g_0_x_0_xxyyy_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxyyy_1[i] * fti_ab_0 + 2.0 * g_0_xx_0_xyyy_1[i] * fi_abcd_0 +
                               g_0_xx_0_xxyyy_0[i] * pb_x + g_0_xx_0_xxyyy_1[i] * wp_x[i];

        g_0_xxx_0_xxyyz_0[i] = 2.0 * g_0_x_0_xxyyz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xx_0_xyyz_1[i] * fi_abcd_0 +
                               g_0_xx_0_xxyyz_0[i] * pb_x + g_0_xx_0_xxyyz_1[i] * wp_x[i];

        g_0_xxx_0_xxyzz_0[i] = 2.0 * g_0_x_0_xxyzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xx_0_xyzz_1[i] * fi_abcd_0 +
                               g_0_xx_0_xxyzz_0[i] * pb_x + g_0_xx_0_xxyzz_1[i] * wp_x[i];

        g_0_xxx_0_xxzzz_0[i] = 2.0 * g_0_x_0_xxzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxzzz_1[i] * fti_ab_0 + 2.0 * g_0_xx_0_xzzz_1[i] * fi_abcd_0 +
                               g_0_xx_0_xxzzz_0[i] * pb_x + g_0_xx_0_xxzzz_1[i] * wp_x[i];

        g_0_xxx_0_xyyyy_0[i] = 2.0 * g_0_x_0_xyyyy_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xyyyy_1[i] * fti_ab_0 + g_0_xx_0_yyyy_1[i] * fi_abcd_0 +
                               g_0_xx_0_xyyyy_0[i] * pb_x + g_0_xx_0_xyyyy_1[i] * wp_x[i];

        g_0_xxx_0_xyyyz_0[i] = 2.0 * g_0_x_0_xyyyz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xyyyz_1[i] * fti_ab_0 + g_0_xx_0_yyyz_1[i] * fi_abcd_0 +
                               g_0_xx_0_xyyyz_0[i] * pb_x + g_0_xx_0_xyyyz_1[i] * wp_x[i];

        g_0_xxx_0_xyyzz_0[i] = 2.0 * g_0_x_0_xyyzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xyyzz_1[i] * fti_ab_0 + g_0_xx_0_yyzz_1[i] * fi_abcd_0 +
                               g_0_xx_0_xyyzz_0[i] * pb_x + g_0_xx_0_xyyzz_1[i] * wp_x[i];

        g_0_xxx_0_xyzzz_0[i] = 2.0 * g_0_x_0_xyzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xyzzz_1[i] * fti_ab_0 + g_0_xx_0_yzzz_1[i] * fi_abcd_0 +
                               g_0_xx_0_xyzzz_0[i] * pb_x + g_0_xx_0_xyzzz_1[i] * wp_x[i];

        g_0_xxx_0_xzzzz_0[i] = 2.0 * g_0_x_0_xzzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xzzzz_1[i] * fti_ab_0 + g_0_xx_0_zzzz_1[i] * fi_abcd_0 +
                               g_0_xx_0_xzzzz_0[i] * pb_x + g_0_xx_0_xzzzz_1[i] * wp_x[i];

        g_0_xxx_0_yyyyy_0[i] =
            2.0 * g_0_x_0_yyyyy_0[i] * fi_ab_0 - 2.0 * g_0_x_0_yyyyy_1[i] * fti_ab_0 + g_0_xx_0_yyyyy_0[i] * pb_x + g_0_xx_0_yyyyy_1[i] * wp_x[i];

        g_0_xxx_0_yyyyz_0[i] =
            2.0 * g_0_x_0_yyyyz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_yyyyz_1[i] * fti_ab_0 + g_0_xx_0_yyyyz_0[i] * pb_x + g_0_xx_0_yyyyz_1[i] * wp_x[i];

        g_0_xxx_0_yyyzz_0[i] =
            2.0 * g_0_x_0_yyyzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_yyyzz_1[i] * fti_ab_0 + g_0_xx_0_yyyzz_0[i] * pb_x + g_0_xx_0_yyyzz_1[i] * wp_x[i];

        g_0_xxx_0_yyzzz_0[i] =
            2.0 * g_0_x_0_yyzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_yyzzz_1[i] * fti_ab_0 + g_0_xx_0_yyzzz_0[i] * pb_x + g_0_xx_0_yyzzz_1[i] * wp_x[i];

        g_0_xxx_0_yzzzz_0[i] =
            2.0 * g_0_x_0_yzzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_yzzzz_1[i] * fti_ab_0 + g_0_xx_0_yzzzz_0[i] * pb_x + g_0_xx_0_yzzzz_1[i] * wp_x[i];

        g_0_xxx_0_zzzzz_0[i] =
            2.0 * g_0_x_0_zzzzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_zzzzz_1[i] * fti_ab_0 + g_0_xx_0_zzzzz_0[i] * pb_x + g_0_xx_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 21-42 components of targeted buffer : SFSH

    auto g_0_xxy_0_xxxxx_0 = pbuffer.data(idx_eri_0_sfsh + 21);

    auto g_0_xxy_0_xxxxy_0 = pbuffer.data(idx_eri_0_sfsh + 22);

    auto g_0_xxy_0_xxxxz_0 = pbuffer.data(idx_eri_0_sfsh + 23);

    auto g_0_xxy_0_xxxyy_0 = pbuffer.data(idx_eri_0_sfsh + 24);

    auto g_0_xxy_0_xxxyz_0 = pbuffer.data(idx_eri_0_sfsh + 25);

    auto g_0_xxy_0_xxxzz_0 = pbuffer.data(idx_eri_0_sfsh + 26);

    auto g_0_xxy_0_xxyyy_0 = pbuffer.data(idx_eri_0_sfsh + 27);

    auto g_0_xxy_0_xxyyz_0 = pbuffer.data(idx_eri_0_sfsh + 28);

    auto g_0_xxy_0_xxyzz_0 = pbuffer.data(idx_eri_0_sfsh + 29);

    auto g_0_xxy_0_xxzzz_0 = pbuffer.data(idx_eri_0_sfsh + 30);

    auto g_0_xxy_0_xyyyy_0 = pbuffer.data(idx_eri_0_sfsh + 31);

    auto g_0_xxy_0_xyyyz_0 = pbuffer.data(idx_eri_0_sfsh + 32);

    auto g_0_xxy_0_xyyzz_0 = pbuffer.data(idx_eri_0_sfsh + 33);

    auto g_0_xxy_0_xyzzz_0 = pbuffer.data(idx_eri_0_sfsh + 34);

    auto g_0_xxy_0_xzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 35);

    auto g_0_xxy_0_yyyyy_0 = pbuffer.data(idx_eri_0_sfsh + 36);

    auto g_0_xxy_0_yyyyz_0 = pbuffer.data(idx_eri_0_sfsh + 37);

    auto g_0_xxy_0_yyyzz_0 = pbuffer.data(idx_eri_0_sfsh + 38);

    auto g_0_xxy_0_yyzzz_0 = pbuffer.data(idx_eri_0_sfsh + 39);

    auto g_0_xxy_0_yzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 40);

    auto g_0_xxy_0_zzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 41);

#pragma omp simd aligned(g_0_xx_0_xxxx_1,       \
                             g_0_xx_0_xxxxx_0,  \
                             g_0_xx_0_xxxxx_1,  \
                             g_0_xx_0_xxxxy_0,  \
                             g_0_xx_0_xxxxy_1,  \
                             g_0_xx_0_xxxxz_0,  \
                             g_0_xx_0_xxxxz_1,  \
                             g_0_xx_0_xxxy_1,   \
                             g_0_xx_0_xxxyy_0,  \
                             g_0_xx_0_xxxyy_1,  \
                             g_0_xx_0_xxxyz_0,  \
                             g_0_xx_0_xxxyz_1,  \
                             g_0_xx_0_xxxz_1,   \
                             g_0_xx_0_xxxzz_0,  \
                             g_0_xx_0_xxxzz_1,  \
                             g_0_xx_0_xxyy_1,   \
                             g_0_xx_0_xxyyy_0,  \
                             g_0_xx_0_xxyyy_1,  \
                             g_0_xx_0_xxyyz_0,  \
                             g_0_xx_0_xxyyz_1,  \
                             g_0_xx_0_xxyz_1,   \
                             g_0_xx_0_xxyzz_0,  \
                             g_0_xx_0_xxyzz_1,  \
                             g_0_xx_0_xxzz_1,   \
                             g_0_xx_0_xxzzz_0,  \
                             g_0_xx_0_xxzzz_1,  \
                             g_0_xx_0_xyyy_1,   \
                             g_0_xx_0_xyyyy_0,  \
                             g_0_xx_0_xyyyy_1,  \
                             g_0_xx_0_xyyyz_0,  \
                             g_0_xx_0_xyyyz_1,  \
                             g_0_xx_0_xyyz_1,   \
                             g_0_xx_0_xyyzz_0,  \
                             g_0_xx_0_xyyzz_1,  \
                             g_0_xx_0_xyzz_1,   \
                             g_0_xx_0_xyzzz_0,  \
                             g_0_xx_0_xyzzz_1,  \
                             g_0_xx_0_xzzz_1,   \
                             g_0_xx_0_xzzzz_0,  \
                             g_0_xx_0_xzzzz_1,  \
                             g_0_xx_0_yyyy_1,   \
                             g_0_xx_0_yyyyy_0,  \
                             g_0_xx_0_yyyyy_1,  \
                             g_0_xx_0_yyyyz_0,  \
                             g_0_xx_0_yyyyz_1,  \
                             g_0_xx_0_yyyz_1,   \
                             g_0_xx_0_yyyzz_0,  \
                             g_0_xx_0_yyyzz_1,  \
                             g_0_xx_0_yyzz_1,   \
                             g_0_xx_0_yyzzz_0,  \
                             g_0_xx_0_yyzzz_1,  \
                             g_0_xx_0_yzzz_1,   \
                             g_0_xx_0_yzzzz_0,  \
                             g_0_xx_0_yzzzz_1,  \
                             g_0_xx_0_zzzz_1,   \
                             g_0_xx_0_zzzzz_0,  \
                             g_0_xx_0_zzzzz_1,  \
                             g_0_xxy_0_xxxxx_0, \
                             g_0_xxy_0_xxxxy_0, \
                             g_0_xxy_0_xxxxz_0, \
                             g_0_xxy_0_xxxyy_0, \
                             g_0_xxy_0_xxxyz_0, \
                             g_0_xxy_0_xxxzz_0, \
                             g_0_xxy_0_xxyyy_0, \
                             g_0_xxy_0_xxyyz_0, \
                             g_0_xxy_0_xxyzz_0, \
                             g_0_xxy_0_xxzzz_0, \
                             g_0_xxy_0_xyyyy_0, \
                             g_0_xxy_0_xyyyz_0, \
                             g_0_xxy_0_xyyzz_0, \
                             g_0_xxy_0_xyzzz_0, \
                             g_0_xxy_0_xzzzz_0, \
                             g_0_xxy_0_yyyyy_0, \
                             g_0_xxy_0_yyyyz_0, \
                             g_0_xxy_0_yyyzz_0, \
                             g_0_xxy_0_yyzzz_0, \
                             g_0_xxy_0_yzzzz_0, \
                             g_0_xxy_0_zzzzz_0, \
                             wp_y,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxy_0_xxxxx_0[i] = g_0_xx_0_xxxxx_0[i] * pb_y + g_0_xx_0_xxxxx_1[i] * wp_y[i];

        g_0_xxy_0_xxxxy_0[i] = g_0_xx_0_xxxx_1[i] * fi_abcd_0 + g_0_xx_0_xxxxy_0[i] * pb_y + g_0_xx_0_xxxxy_1[i] * wp_y[i];

        g_0_xxy_0_xxxxz_0[i] = g_0_xx_0_xxxxz_0[i] * pb_y + g_0_xx_0_xxxxz_1[i] * wp_y[i];

        g_0_xxy_0_xxxyy_0[i] = 2.0 * g_0_xx_0_xxxy_1[i] * fi_abcd_0 + g_0_xx_0_xxxyy_0[i] * pb_y + g_0_xx_0_xxxyy_1[i] * wp_y[i];

        g_0_xxy_0_xxxyz_0[i] = g_0_xx_0_xxxz_1[i] * fi_abcd_0 + g_0_xx_0_xxxyz_0[i] * pb_y + g_0_xx_0_xxxyz_1[i] * wp_y[i];

        g_0_xxy_0_xxxzz_0[i] = g_0_xx_0_xxxzz_0[i] * pb_y + g_0_xx_0_xxxzz_1[i] * wp_y[i];

        g_0_xxy_0_xxyyy_0[i] = 3.0 * g_0_xx_0_xxyy_1[i] * fi_abcd_0 + g_0_xx_0_xxyyy_0[i] * pb_y + g_0_xx_0_xxyyy_1[i] * wp_y[i];

        g_0_xxy_0_xxyyz_0[i] = 2.0 * g_0_xx_0_xxyz_1[i] * fi_abcd_0 + g_0_xx_0_xxyyz_0[i] * pb_y + g_0_xx_0_xxyyz_1[i] * wp_y[i];

        g_0_xxy_0_xxyzz_0[i] = g_0_xx_0_xxzz_1[i] * fi_abcd_0 + g_0_xx_0_xxyzz_0[i] * pb_y + g_0_xx_0_xxyzz_1[i] * wp_y[i];

        g_0_xxy_0_xxzzz_0[i] = g_0_xx_0_xxzzz_0[i] * pb_y + g_0_xx_0_xxzzz_1[i] * wp_y[i];

        g_0_xxy_0_xyyyy_0[i] = 4.0 * g_0_xx_0_xyyy_1[i] * fi_abcd_0 + g_0_xx_0_xyyyy_0[i] * pb_y + g_0_xx_0_xyyyy_1[i] * wp_y[i];

        g_0_xxy_0_xyyyz_0[i] = 3.0 * g_0_xx_0_xyyz_1[i] * fi_abcd_0 + g_0_xx_0_xyyyz_0[i] * pb_y + g_0_xx_0_xyyyz_1[i] * wp_y[i];

        g_0_xxy_0_xyyzz_0[i] = 2.0 * g_0_xx_0_xyzz_1[i] * fi_abcd_0 + g_0_xx_0_xyyzz_0[i] * pb_y + g_0_xx_0_xyyzz_1[i] * wp_y[i];

        g_0_xxy_0_xyzzz_0[i] = g_0_xx_0_xzzz_1[i] * fi_abcd_0 + g_0_xx_0_xyzzz_0[i] * pb_y + g_0_xx_0_xyzzz_1[i] * wp_y[i];

        g_0_xxy_0_xzzzz_0[i] = g_0_xx_0_xzzzz_0[i] * pb_y + g_0_xx_0_xzzzz_1[i] * wp_y[i];

        g_0_xxy_0_yyyyy_0[i] = 5.0 * g_0_xx_0_yyyy_1[i] * fi_abcd_0 + g_0_xx_0_yyyyy_0[i] * pb_y + g_0_xx_0_yyyyy_1[i] * wp_y[i];

        g_0_xxy_0_yyyyz_0[i] = 4.0 * g_0_xx_0_yyyz_1[i] * fi_abcd_0 + g_0_xx_0_yyyyz_0[i] * pb_y + g_0_xx_0_yyyyz_1[i] * wp_y[i];

        g_0_xxy_0_yyyzz_0[i] = 3.0 * g_0_xx_0_yyzz_1[i] * fi_abcd_0 + g_0_xx_0_yyyzz_0[i] * pb_y + g_0_xx_0_yyyzz_1[i] * wp_y[i];

        g_0_xxy_0_yyzzz_0[i] = 2.0 * g_0_xx_0_yzzz_1[i] * fi_abcd_0 + g_0_xx_0_yyzzz_0[i] * pb_y + g_0_xx_0_yyzzz_1[i] * wp_y[i];

        g_0_xxy_0_yzzzz_0[i] = g_0_xx_0_zzzz_1[i] * fi_abcd_0 + g_0_xx_0_yzzzz_0[i] * pb_y + g_0_xx_0_yzzzz_1[i] * wp_y[i];

        g_0_xxy_0_zzzzz_0[i] = g_0_xx_0_zzzzz_0[i] * pb_y + g_0_xx_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 42-63 components of targeted buffer : SFSH

    auto g_0_xxz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sfsh + 42);

    auto g_0_xxz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sfsh + 43);

    auto g_0_xxz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sfsh + 44);

    auto g_0_xxz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sfsh + 45);

    auto g_0_xxz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sfsh + 46);

    auto g_0_xxz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sfsh + 47);

    auto g_0_xxz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sfsh + 48);

    auto g_0_xxz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sfsh + 49);

    auto g_0_xxz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sfsh + 50);

    auto g_0_xxz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sfsh + 51);

    auto g_0_xxz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sfsh + 52);

    auto g_0_xxz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sfsh + 53);

    auto g_0_xxz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sfsh + 54);

    auto g_0_xxz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sfsh + 55);

    auto g_0_xxz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 56);

    auto g_0_xxz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sfsh + 57);

    auto g_0_xxz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sfsh + 58);

    auto g_0_xxz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sfsh + 59);

    auto g_0_xxz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sfsh + 60);

    auto g_0_xxz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 61);

    auto g_0_xxz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 62);

#pragma omp simd aligned(g_0_xx_0_xxxx_1,       \
                             g_0_xx_0_xxxxx_0,  \
                             g_0_xx_0_xxxxx_1,  \
                             g_0_xx_0_xxxxy_0,  \
                             g_0_xx_0_xxxxy_1,  \
                             g_0_xx_0_xxxxz_0,  \
                             g_0_xx_0_xxxxz_1,  \
                             g_0_xx_0_xxxy_1,   \
                             g_0_xx_0_xxxyy_0,  \
                             g_0_xx_0_xxxyy_1,  \
                             g_0_xx_0_xxxyz_0,  \
                             g_0_xx_0_xxxyz_1,  \
                             g_0_xx_0_xxxz_1,   \
                             g_0_xx_0_xxxzz_0,  \
                             g_0_xx_0_xxxzz_1,  \
                             g_0_xx_0_xxyy_1,   \
                             g_0_xx_0_xxyyy_0,  \
                             g_0_xx_0_xxyyy_1,  \
                             g_0_xx_0_xxyyz_0,  \
                             g_0_xx_0_xxyyz_1,  \
                             g_0_xx_0_xxyz_1,   \
                             g_0_xx_0_xxyzz_0,  \
                             g_0_xx_0_xxyzz_1,  \
                             g_0_xx_0_xxzz_1,   \
                             g_0_xx_0_xxzzz_0,  \
                             g_0_xx_0_xxzzz_1,  \
                             g_0_xx_0_xyyy_1,   \
                             g_0_xx_0_xyyyy_0,  \
                             g_0_xx_0_xyyyy_1,  \
                             g_0_xx_0_xyyyz_0,  \
                             g_0_xx_0_xyyyz_1,  \
                             g_0_xx_0_xyyz_1,   \
                             g_0_xx_0_xyyzz_0,  \
                             g_0_xx_0_xyyzz_1,  \
                             g_0_xx_0_xyzz_1,   \
                             g_0_xx_0_xyzzz_0,  \
                             g_0_xx_0_xyzzz_1,  \
                             g_0_xx_0_xzzz_1,   \
                             g_0_xx_0_xzzzz_0,  \
                             g_0_xx_0_xzzzz_1,  \
                             g_0_xx_0_yyyy_1,   \
                             g_0_xx_0_yyyyy_0,  \
                             g_0_xx_0_yyyyy_1,  \
                             g_0_xx_0_yyyyz_0,  \
                             g_0_xx_0_yyyyz_1,  \
                             g_0_xx_0_yyyz_1,   \
                             g_0_xx_0_yyyzz_0,  \
                             g_0_xx_0_yyyzz_1,  \
                             g_0_xx_0_yyzz_1,   \
                             g_0_xx_0_yyzzz_0,  \
                             g_0_xx_0_yyzzz_1,  \
                             g_0_xx_0_yzzz_1,   \
                             g_0_xx_0_yzzzz_0,  \
                             g_0_xx_0_yzzzz_1,  \
                             g_0_xx_0_zzzz_1,   \
                             g_0_xx_0_zzzzz_0,  \
                             g_0_xx_0_zzzzz_1,  \
                             g_0_xxz_0_xxxxx_0, \
                             g_0_xxz_0_xxxxy_0, \
                             g_0_xxz_0_xxxxz_0, \
                             g_0_xxz_0_xxxyy_0, \
                             g_0_xxz_0_xxxyz_0, \
                             g_0_xxz_0_xxxzz_0, \
                             g_0_xxz_0_xxyyy_0, \
                             g_0_xxz_0_xxyyz_0, \
                             g_0_xxz_0_xxyzz_0, \
                             g_0_xxz_0_xxzzz_0, \
                             g_0_xxz_0_xyyyy_0, \
                             g_0_xxz_0_xyyyz_0, \
                             g_0_xxz_0_xyyzz_0, \
                             g_0_xxz_0_xyzzz_0, \
                             g_0_xxz_0_xzzzz_0, \
                             g_0_xxz_0_yyyyy_0, \
                             g_0_xxz_0_yyyyz_0, \
                             g_0_xxz_0_yyyzz_0, \
                             g_0_xxz_0_yyzzz_0, \
                             g_0_xxz_0_yzzzz_0, \
                             g_0_xxz_0_zzzzz_0, \
                             wp_z,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxz_0_xxxxx_0[i] = g_0_xx_0_xxxxx_0[i] * pb_z + g_0_xx_0_xxxxx_1[i] * wp_z[i];

        g_0_xxz_0_xxxxy_0[i] = g_0_xx_0_xxxxy_0[i] * pb_z + g_0_xx_0_xxxxy_1[i] * wp_z[i];

        g_0_xxz_0_xxxxz_0[i] = g_0_xx_0_xxxx_1[i] * fi_abcd_0 + g_0_xx_0_xxxxz_0[i] * pb_z + g_0_xx_0_xxxxz_1[i] * wp_z[i];

        g_0_xxz_0_xxxyy_0[i] = g_0_xx_0_xxxyy_0[i] * pb_z + g_0_xx_0_xxxyy_1[i] * wp_z[i];

        g_0_xxz_0_xxxyz_0[i] = g_0_xx_0_xxxy_1[i] * fi_abcd_0 + g_0_xx_0_xxxyz_0[i] * pb_z + g_0_xx_0_xxxyz_1[i] * wp_z[i];

        g_0_xxz_0_xxxzz_0[i] = 2.0 * g_0_xx_0_xxxz_1[i] * fi_abcd_0 + g_0_xx_0_xxxzz_0[i] * pb_z + g_0_xx_0_xxxzz_1[i] * wp_z[i];

        g_0_xxz_0_xxyyy_0[i] = g_0_xx_0_xxyyy_0[i] * pb_z + g_0_xx_0_xxyyy_1[i] * wp_z[i];

        g_0_xxz_0_xxyyz_0[i] = g_0_xx_0_xxyy_1[i] * fi_abcd_0 + g_0_xx_0_xxyyz_0[i] * pb_z + g_0_xx_0_xxyyz_1[i] * wp_z[i];

        g_0_xxz_0_xxyzz_0[i] = 2.0 * g_0_xx_0_xxyz_1[i] * fi_abcd_0 + g_0_xx_0_xxyzz_0[i] * pb_z + g_0_xx_0_xxyzz_1[i] * wp_z[i];

        g_0_xxz_0_xxzzz_0[i] = 3.0 * g_0_xx_0_xxzz_1[i] * fi_abcd_0 + g_0_xx_0_xxzzz_0[i] * pb_z + g_0_xx_0_xxzzz_1[i] * wp_z[i];

        g_0_xxz_0_xyyyy_0[i] = g_0_xx_0_xyyyy_0[i] * pb_z + g_0_xx_0_xyyyy_1[i] * wp_z[i];

        g_0_xxz_0_xyyyz_0[i] = g_0_xx_0_xyyy_1[i] * fi_abcd_0 + g_0_xx_0_xyyyz_0[i] * pb_z + g_0_xx_0_xyyyz_1[i] * wp_z[i];

        g_0_xxz_0_xyyzz_0[i] = 2.0 * g_0_xx_0_xyyz_1[i] * fi_abcd_0 + g_0_xx_0_xyyzz_0[i] * pb_z + g_0_xx_0_xyyzz_1[i] * wp_z[i];

        g_0_xxz_0_xyzzz_0[i] = 3.0 * g_0_xx_0_xyzz_1[i] * fi_abcd_0 + g_0_xx_0_xyzzz_0[i] * pb_z + g_0_xx_0_xyzzz_1[i] * wp_z[i];

        g_0_xxz_0_xzzzz_0[i] = 4.0 * g_0_xx_0_xzzz_1[i] * fi_abcd_0 + g_0_xx_0_xzzzz_0[i] * pb_z + g_0_xx_0_xzzzz_1[i] * wp_z[i];

        g_0_xxz_0_yyyyy_0[i] = g_0_xx_0_yyyyy_0[i] * pb_z + g_0_xx_0_yyyyy_1[i] * wp_z[i];

        g_0_xxz_0_yyyyz_0[i] = g_0_xx_0_yyyy_1[i] * fi_abcd_0 + g_0_xx_0_yyyyz_0[i] * pb_z + g_0_xx_0_yyyyz_1[i] * wp_z[i];

        g_0_xxz_0_yyyzz_0[i] = 2.0 * g_0_xx_0_yyyz_1[i] * fi_abcd_0 + g_0_xx_0_yyyzz_0[i] * pb_z + g_0_xx_0_yyyzz_1[i] * wp_z[i];

        g_0_xxz_0_yyzzz_0[i] = 3.0 * g_0_xx_0_yyzz_1[i] * fi_abcd_0 + g_0_xx_0_yyzzz_0[i] * pb_z + g_0_xx_0_yyzzz_1[i] * wp_z[i];

        g_0_xxz_0_yzzzz_0[i] = 4.0 * g_0_xx_0_yzzz_1[i] * fi_abcd_0 + g_0_xx_0_yzzzz_0[i] * pb_z + g_0_xx_0_yzzzz_1[i] * wp_z[i];

        g_0_xxz_0_zzzzz_0[i] = 5.0 * g_0_xx_0_zzzz_1[i] * fi_abcd_0 + g_0_xx_0_zzzzz_0[i] * pb_z + g_0_xx_0_zzzzz_1[i] * wp_z[i];
    }

    /// Set up 63-84 components of targeted buffer : SFSH

    auto g_0_xyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_sfsh + 63);

    auto g_0_xyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_sfsh + 64);

    auto g_0_xyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_sfsh + 65);

    auto g_0_xyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_sfsh + 66);

    auto g_0_xyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_sfsh + 67);

    auto g_0_xyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_sfsh + 68);

    auto g_0_xyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_sfsh + 69);

    auto g_0_xyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_sfsh + 70);

    auto g_0_xyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_sfsh + 71);

    auto g_0_xyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_sfsh + 72);

    auto g_0_xyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_sfsh + 73);

    auto g_0_xyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_sfsh + 74);

    auto g_0_xyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_sfsh + 75);

    auto g_0_xyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_sfsh + 76);

    auto g_0_xyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 77);

    auto g_0_xyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_sfsh + 78);

    auto g_0_xyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_sfsh + 79);

    auto g_0_xyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_sfsh + 80);

    auto g_0_xyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_sfsh + 81);

    auto g_0_xyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 82);

    auto g_0_xyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 83);

#pragma omp simd aligned(g_0_xyy_0_xxxxx_0,     \
                             g_0_xyy_0_xxxxy_0, \
                             g_0_xyy_0_xxxxz_0, \
                             g_0_xyy_0_xxxyy_0, \
                             g_0_xyy_0_xxxyz_0, \
                             g_0_xyy_0_xxxzz_0, \
                             g_0_xyy_0_xxyyy_0, \
                             g_0_xyy_0_xxyyz_0, \
                             g_0_xyy_0_xxyzz_0, \
                             g_0_xyy_0_xxzzz_0, \
                             g_0_xyy_0_xyyyy_0, \
                             g_0_xyy_0_xyyyz_0, \
                             g_0_xyy_0_xyyzz_0, \
                             g_0_xyy_0_xyzzz_0, \
                             g_0_xyy_0_xzzzz_0, \
                             g_0_xyy_0_yyyyy_0, \
                             g_0_xyy_0_yyyyz_0, \
                             g_0_xyy_0_yyyzz_0, \
                             g_0_xyy_0_yyzzz_0, \
                             g_0_xyy_0_yzzzz_0, \
                             g_0_xyy_0_zzzzz_0, \
                             g_0_yy_0_xxxx_1,   \
                             g_0_yy_0_xxxxx_0,  \
                             g_0_yy_0_xxxxx_1,  \
                             g_0_yy_0_xxxxy_0,  \
                             g_0_yy_0_xxxxy_1,  \
                             g_0_yy_0_xxxxz_0,  \
                             g_0_yy_0_xxxxz_1,  \
                             g_0_yy_0_xxxy_1,   \
                             g_0_yy_0_xxxyy_0,  \
                             g_0_yy_0_xxxyy_1,  \
                             g_0_yy_0_xxxyz_0,  \
                             g_0_yy_0_xxxyz_1,  \
                             g_0_yy_0_xxxz_1,   \
                             g_0_yy_0_xxxzz_0,  \
                             g_0_yy_0_xxxzz_1,  \
                             g_0_yy_0_xxyy_1,   \
                             g_0_yy_0_xxyyy_0,  \
                             g_0_yy_0_xxyyy_1,  \
                             g_0_yy_0_xxyyz_0,  \
                             g_0_yy_0_xxyyz_1,  \
                             g_0_yy_0_xxyz_1,   \
                             g_0_yy_0_xxyzz_0,  \
                             g_0_yy_0_xxyzz_1,  \
                             g_0_yy_0_xxzz_1,   \
                             g_0_yy_0_xxzzz_0,  \
                             g_0_yy_0_xxzzz_1,  \
                             g_0_yy_0_xyyy_1,   \
                             g_0_yy_0_xyyyy_0,  \
                             g_0_yy_0_xyyyy_1,  \
                             g_0_yy_0_xyyyz_0,  \
                             g_0_yy_0_xyyyz_1,  \
                             g_0_yy_0_xyyz_1,   \
                             g_0_yy_0_xyyzz_0,  \
                             g_0_yy_0_xyyzz_1,  \
                             g_0_yy_0_xyzz_1,   \
                             g_0_yy_0_xyzzz_0,  \
                             g_0_yy_0_xyzzz_1,  \
                             g_0_yy_0_xzzz_1,   \
                             g_0_yy_0_xzzzz_0,  \
                             g_0_yy_0_xzzzz_1,  \
                             g_0_yy_0_yyyy_1,   \
                             g_0_yy_0_yyyyy_0,  \
                             g_0_yy_0_yyyyy_1,  \
                             g_0_yy_0_yyyyz_0,  \
                             g_0_yy_0_yyyyz_1,  \
                             g_0_yy_0_yyyz_1,   \
                             g_0_yy_0_yyyzz_0,  \
                             g_0_yy_0_yyyzz_1,  \
                             g_0_yy_0_yyzz_1,   \
                             g_0_yy_0_yyzzz_0,  \
                             g_0_yy_0_yyzzz_1,  \
                             g_0_yy_0_yzzz_1,   \
                             g_0_yy_0_yzzzz_0,  \
                             g_0_yy_0_yzzzz_1,  \
                             g_0_yy_0_zzzz_1,   \
                             g_0_yy_0_zzzzz_0,  \
                             g_0_yy_0_zzzzz_1,  \
                             wp_x,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyy_0_xxxxx_0[i] = 5.0 * g_0_yy_0_xxxx_1[i] * fi_abcd_0 + g_0_yy_0_xxxxx_0[i] * pb_x + g_0_yy_0_xxxxx_1[i] * wp_x[i];

        g_0_xyy_0_xxxxy_0[i] = 4.0 * g_0_yy_0_xxxy_1[i] * fi_abcd_0 + g_0_yy_0_xxxxy_0[i] * pb_x + g_0_yy_0_xxxxy_1[i] * wp_x[i];

        g_0_xyy_0_xxxxz_0[i] = 4.0 * g_0_yy_0_xxxz_1[i] * fi_abcd_0 + g_0_yy_0_xxxxz_0[i] * pb_x + g_0_yy_0_xxxxz_1[i] * wp_x[i];

        g_0_xyy_0_xxxyy_0[i] = 3.0 * g_0_yy_0_xxyy_1[i] * fi_abcd_0 + g_0_yy_0_xxxyy_0[i] * pb_x + g_0_yy_0_xxxyy_1[i] * wp_x[i];

        g_0_xyy_0_xxxyz_0[i] = 3.0 * g_0_yy_0_xxyz_1[i] * fi_abcd_0 + g_0_yy_0_xxxyz_0[i] * pb_x + g_0_yy_0_xxxyz_1[i] * wp_x[i];

        g_0_xyy_0_xxxzz_0[i] = 3.0 * g_0_yy_0_xxzz_1[i] * fi_abcd_0 + g_0_yy_0_xxxzz_0[i] * pb_x + g_0_yy_0_xxxzz_1[i] * wp_x[i];

        g_0_xyy_0_xxyyy_0[i] = 2.0 * g_0_yy_0_xyyy_1[i] * fi_abcd_0 + g_0_yy_0_xxyyy_0[i] * pb_x + g_0_yy_0_xxyyy_1[i] * wp_x[i];

        g_0_xyy_0_xxyyz_0[i] = 2.0 * g_0_yy_0_xyyz_1[i] * fi_abcd_0 + g_0_yy_0_xxyyz_0[i] * pb_x + g_0_yy_0_xxyyz_1[i] * wp_x[i];

        g_0_xyy_0_xxyzz_0[i] = 2.0 * g_0_yy_0_xyzz_1[i] * fi_abcd_0 + g_0_yy_0_xxyzz_0[i] * pb_x + g_0_yy_0_xxyzz_1[i] * wp_x[i];

        g_0_xyy_0_xxzzz_0[i] = 2.0 * g_0_yy_0_xzzz_1[i] * fi_abcd_0 + g_0_yy_0_xxzzz_0[i] * pb_x + g_0_yy_0_xxzzz_1[i] * wp_x[i];

        g_0_xyy_0_xyyyy_0[i] = g_0_yy_0_yyyy_1[i] * fi_abcd_0 + g_0_yy_0_xyyyy_0[i] * pb_x + g_0_yy_0_xyyyy_1[i] * wp_x[i];

        g_0_xyy_0_xyyyz_0[i] = g_0_yy_0_yyyz_1[i] * fi_abcd_0 + g_0_yy_0_xyyyz_0[i] * pb_x + g_0_yy_0_xyyyz_1[i] * wp_x[i];

        g_0_xyy_0_xyyzz_0[i] = g_0_yy_0_yyzz_1[i] * fi_abcd_0 + g_0_yy_0_xyyzz_0[i] * pb_x + g_0_yy_0_xyyzz_1[i] * wp_x[i];

        g_0_xyy_0_xyzzz_0[i] = g_0_yy_0_yzzz_1[i] * fi_abcd_0 + g_0_yy_0_xyzzz_0[i] * pb_x + g_0_yy_0_xyzzz_1[i] * wp_x[i];

        g_0_xyy_0_xzzzz_0[i] = g_0_yy_0_zzzz_1[i] * fi_abcd_0 + g_0_yy_0_xzzzz_0[i] * pb_x + g_0_yy_0_xzzzz_1[i] * wp_x[i];

        g_0_xyy_0_yyyyy_0[i] = g_0_yy_0_yyyyy_0[i] * pb_x + g_0_yy_0_yyyyy_1[i] * wp_x[i];

        g_0_xyy_0_yyyyz_0[i] = g_0_yy_0_yyyyz_0[i] * pb_x + g_0_yy_0_yyyyz_1[i] * wp_x[i];

        g_0_xyy_0_yyyzz_0[i] = g_0_yy_0_yyyzz_0[i] * pb_x + g_0_yy_0_yyyzz_1[i] * wp_x[i];

        g_0_xyy_0_yyzzz_0[i] = g_0_yy_0_yyzzz_0[i] * pb_x + g_0_yy_0_yyzzz_1[i] * wp_x[i];

        g_0_xyy_0_yzzzz_0[i] = g_0_yy_0_yzzzz_0[i] * pb_x + g_0_yy_0_yzzzz_1[i] * wp_x[i];

        g_0_xyy_0_zzzzz_0[i] = g_0_yy_0_zzzzz_0[i] * pb_x + g_0_yy_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 84-105 components of targeted buffer : SFSH

    auto g_0_xyz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sfsh + 84);

    auto g_0_xyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sfsh + 85);

    auto g_0_xyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sfsh + 86);

    auto g_0_xyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sfsh + 87);

    auto g_0_xyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sfsh + 88);

    auto g_0_xyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sfsh + 89);

    auto g_0_xyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sfsh + 90);

    auto g_0_xyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sfsh + 91);

    auto g_0_xyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sfsh + 92);

    auto g_0_xyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sfsh + 93);

    auto g_0_xyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sfsh + 94);

    auto g_0_xyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sfsh + 95);

    auto g_0_xyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sfsh + 96);

    auto g_0_xyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sfsh + 97);

    auto g_0_xyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 98);

    auto g_0_xyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sfsh + 99);

    auto g_0_xyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sfsh + 100);

    auto g_0_xyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sfsh + 101);

    auto g_0_xyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sfsh + 102);

    auto g_0_xyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 103);

    auto g_0_xyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 104);

#pragma omp simd aligned(g_0_xy_0_xxxxy_0,      \
                             g_0_xy_0_xxxxy_1,  \
                             g_0_xy_0_xxxyy_0,  \
                             g_0_xy_0_xxxyy_1,  \
                             g_0_xy_0_xxyyy_0,  \
                             g_0_xy_0_xxyyy_1,  \
                             g_0_xy_0_xyyyy_0,  \
                             g_0_xy_0_xyyyy_1,  \
                             g_0_xyz_0_xxxxx_0, \
                             g_0_xyz_0_xxxxy_0, \
                             g_0_xyz_0_xxxxz_0, \
                             g_0_xyz_0_xxxyy_0, \
                             g_0_xyz_0_xxxyz_0, \
                             g_0_xyz_0_xxxzz_0, \
                             g_0_xyz_0_xxyyy_0, \
                             g_0_xyz_0_xxyyz_0, \
                             g_0_xyz_0_xxyzz_0, \
                             g_0_xyz_0_xxzzz_0, \
                             g_0_xyz_0_xyyyy_0, \
                             g_0_xyz_0_xyyyz_0, \
                             g_0_xyz_0_xyyzz_0, \
                             g_0_xyz_0_xyzzz_0, \
                             g_0_xyz_0_xzzzz_0, \
                             g_0_xyz_0_yyyyy_0, \
                             g_0_xyz_0_yyyyz_0, \
                             g_0_xyz_0_yyyzz_0, \
                             g_0_xyz_0_yyzzz_0, \
                             g_0_xyz_0_yzzzz_0, \
                             g_0_xyz_0_zzzzz_0, \
                             g_0_xz_0_xxxxx_0,  \
                             g_0_xz_0_xxxxx_1,  \
                             g_0_xz_0_xxxxz_0,  \
                             g_0_xz_0_xxxxz_1,  \
                             g_0_xz_0_xxxzz_0,  \
                             g_0_xz_0_xxxzz_1,  \
                             g_0_xz_0_xxzzz_0,  \
                             g_0_xz_0_xxzzz_1,  \
                             g_0_xz_0_xzzzz_0,  \
                             g_0_xz_0_xzzzz_1,  \
                             g_0_yz_0_xxxyz_0,  \
                             g_0_yz_0_xxxyz_1,  \
                             g_0_yz_0_xxyyz_0,  \
                             g_0_yz_0_xxyyz_1,  \
                             g_0_yz_0_xxyz_1,   \
                             g_0_yz_0_xxyzz_0,  \
                             g_0_yz_0_xxyzz_1,  \
                             g_0_yz_0_xyyyz_0,  \
                             g_0_yz_0_xyyyz_1,  \
                             g_0_yz_0_xyyz_1,   \
                             g_0_yz_0_xyyzz_0,  \
                             g_0_yz_0_xyyzz_1,  \
                             g_0_yz_0_xyzz_1,   \
                             g_0_yz_0_xyzzz_0,  \
                             g_0_yz_0_xyzzz_1,  \
                             g_0_yz_0_yyyyy_0,  \
                             g_0_yz_0_yyyyy_1,  \
                             g_0_yz_0_yyyyz_0,  \
                             g_0_yz_0_yyyyz_1,  \
                             g_0_yz_0_yyyz_1,   \
                             g_0_yz_0_yyyzz_0,  \
                             g_0_yz_0_yyyzz_1,  \
                             g_0_yz_0_yyzz_1,   \
                             g_0_yz_0_yyzzz_0,  \
                             g_0_yz_0_yyzzz_1,  \
                             g_0_yz_0_yzzz_1,   \
                             g_0_yz_0_yzzzz_0,  \
                             g_0_yz_0_yzzzz_1,  \
                             g_0_yz_0_zzzzz_0,  \
                             g_0_yz_0_zzzzz_1,  \
                             wp_x,              \
                             wp_y,              \
                             wp_z,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyz_0_xxxxx_0[i] = g_0_xz_0_xxxxx_0[i] * pb_y + g_0_xz_0_xxxxx_1[i] * wp_y[i];

        g_0_xyz_0_xxxxy_0[i] = g_0_xy_0_xxxxy_0[i] * pb_z + g_0_xy_0_xxxxy_1[i] * wp_z[i];

        g_0_xyz_0_xxxxz_0[i] = g_0_xz_0_xxxxz_0[i] * pb_y + g_0_xz_0_xxxxz_1[i] * wp_y[i];

        g_0_xyz_0_xxxyy_0[i] = g_0_xy_0_xxxyy_0[i] * pb_z + g_0_xy_0_xxxyy_1[i] * wp_z[i];

        g_0_xyz_0_xxxyz_0[i] = 3.0 * g_0_yz_0_xxyz_1[i] * fi_abcd_0 + g_0_yz_0_xxxyz_0[i] * pb_x + g_0_yz_0_xxxyz_1[i] * wp_x[i];

        g_0_xyz_0_xxxzz_0[i] = g_0_xz_0_xxxzz_0[i] * pb_y + g_0_xz_0_xxxzz_1[i] * wp_y[i];

        g_0_xyz_0_xxyyy_0[i] = g_0_xy_0_xxyyy_0[i] * pb_z + g_0_xy_0_xxyyy_1[i] * wp_z[i];

        g_0_xyz_0_xxyyz_0[i] = 2.0 * g_0_yz_0_xyyz_1[i] * fi_abcd_0 + g_0_yz_0_xxyyz_0[i] * pb_x + g_0_yz_0_xxyyz_1[i] * wp_x[i];

        g_0_xyz_0_xxyzz_0[i] = 2.0 * g_0_yz_0_xyzz_1[i] * fi_abcd_0 + g_0_yz_0_xxyzz_0[i] * pb_x + g_0_yz_0_xxyzz_1[i] * wp_x[i];

        g_0_xyz_0_xxzzz_0[i] = g_0_xz_0_xxzzz_0[i] * pb_y + g_0_xz_0_xxzzz_1[i] * wp_y[i];

        g_0_xyz_0_xyyyy_0[i] = g_0_xy_0_xyyyy_0[i] * pb_z + g_0_xy_0_xyyyy_1[i] * wp_z[i];

        g_0_xyz_0_xyyyz_0[i] = g_0_yz_0_yyyz_1[i] * fi_abcd_0 + g_0_yz_0_xyyyz_0[i] * pb_x + g_0_yz_0_xyyyz_1[i] * wp_x[i];

        g_0_xyz_0_xyyzz_0[i] = g_0_yz_0_yyzz_1[i] * fi_abcd_0 + g_0_yz_0_xyyzz_0[i] * pb_x + g_0_yz_0_xyyzz_1[i] * wp_x[i];

        g_0_xyz_0_xyzzz_0[i] = g_0_yz_0_yzzz_1[i] * fi_abcd_0 + g_0_yz_0_xyzzz_0[i] * pb_x + g_0_yz_0_xyzzz_1[i] * wp_x[i];

        g_0_xyz_0_xzzzz_0[i] = g_0_xz_0_xzzzz_0[i] * pb_y + g_0_xz_0_xzzzz_1[i] * wp_y[i];

        g_0_xyz_0_yyyyy_0[i] = g_0_yz_0_yyyyy_0[i] * pb_x + g_0_yz_0_yyyyy_1[i] * wp_x[i];

        g_0_xyz_0_yyyyz_0[i] = g_0_yz_0_yyyyz_0[i] * pb_x + g_0_yz_0_yyyyz_1[i] * wp_x[i];

        g_0_xyz_0_yyyzz_0[i] = g_0_yz_0_yyyzz_0[i] * pb_x + g_0_yz_0_yyyzz_1[i] * wp_x[i];

        g_0_xyz_0_yyzzz_0[i] = g_0_yz_0_yyzzz_0[i] * pb_x + g_0_yz_0_yyzzz_1[i] * wp_x[i];

        g_0_xyz_0_yzzzz_0[i] = g_0_yz_0_yzzzz_0[i] * pb_x + g_0_yz_0_yzzzz_1[i] * wp_x[i];

        g_0_xyz_0_zzzzz_0[i] = g_0_yz_0_zzzzz_0[i] * pb_x + g_0_yz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 105-126 components of targeted buffer : SFSH

    auto g_0_xzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sfsh + 105);

    auto g_0_xzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sfsh + 106);

    auto g_0_xzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sfsh + 107);

    auto g_0_xzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sfsh + 108);

    auto g_0_xzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sfsh + 109);

    auto g_0_xzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sfsh + 110);

    auto g_0_xzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sfsh + 111);

    auto g_0_xzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sfsh + 112);

    auto g_0_xzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sfsh + 113);

    auto g_0_xzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sfsh + 114);

    auto g_0_xzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sfsh + 115);

    auto g_0_xzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sfsh + 116);

    auto g_0_xzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sfsh + 117);

    auto g_0_xzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sfsh + 118);

    auto g_0_xzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 119);

    auto g_0_xzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sfsh + 120);

    auto g_0_xzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sfsh + 121);

    auto g_0_xzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sfsh + 122);

    auto g_0_xzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sfsh + 123);

    auto g_0_xzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 124);

    auto g_0_xzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 125);

#pragma omp simd aligned(g_0_xzz_0_xxxxx_0,     \
                             g_0_xzz_0_xxxxy_0, \
                             g_0_xzz_0_xxxxz_0, \
                             g_0_xzz_0_xxxyy_0, \
                             g_0_xzz_0_xxxyz_0, \
                             g_0_xzz_0_xxxzz_0, \
                             g_0_xzz_0_xxyyy_0, \
                             g_0_xzz_0_xxyyz_0, \
                             g_0_xzz_0_xxyzz_0, \
                             g_0_xzz_0_xxzzz_0, \
                             g_0_xzz_0_xyyyy_0, \
                             g_0_xzz_0_xyyyz_0, \
                             g_0_xzz_0_xyyzz_0, \
                             g_0_xzz_0_xyzzz_0, \
                             g_0_xzz_0_xzzzz_0, \
                             g_0_xzz_0_yyyyy_0, \
                             g_0_xzz_0_yyyyz_0, \
                             g_0_xzz_0_yyyzz_0, \
                             g_0_xzz_0_yyzzz_0, \
                             g_0_xzz_0_yzzzz_0, \
                             g_0_xzz_0_zzzzz_0, \
                             g_0_zz_0_xxxx_1,   \
                             g_0_zz_0_xxxxx_0,  \
                             g_0_zz_0_xxxxx_1,  \
                             g_0_zz_0_xxxxy_0,  \
                             g_0_zz_0_xxxxy_1,  \
                             g_0_zz_0_xxxxz_0,  \
                             g_0_zz_0_xxxxz_1,  \
                             g_0_zz_0_xxxy_1,   \
                             g_0_zz_0_xxxyy_0,  \
                             g_0_zz_0_xxxyy_1,  \
                             g_0_zz_0_xxxyz_0,  \
                             g_0_zz_0_xxxyz_1,  \
                             g_0_zz_0_xxxz_1,   \
                             g_0_zz_0_xxxzz_0,  \
                             g_0_zz_0_xxxzz_1,  \
                             g_0_zz_0_xxyy_1,   \
                             g_0_zz_0_xxyyy_0,  \
                             g_0_zz_0_xxyyy_1,  \
                             g_0_zz_0_xxyyz_0,  \
                             g_0_zz_0_xxyyz_1,  \
                             g_0_zz_0_xxyz_1,   \
                             g_0_zz_0_xxyzz_0,  \
                             g_0_zz_0_xxyzz_1,  \
                             g_0_zz_0_xxzz_1,   \
                             g_0_zz_0_xxzzz_0,  \
                             g_0_zz_0_xxzzz_1,  \
                             g_0_zz_0_xyyy_1,   \
                             g_0_zz_0_xyyyy_0,  \
                             g_0_zz_0_xyyyy_1,  \
                             g_0_zz_0_xyyyz_0,  \
                             g_0_zz_0_xyyyz_1,  \
                             g_0_zz_0_xyyz_1,   \
                             g_0_zz_0_xyyzz_0,  \
                             g_0_zz_0_xyyzz_1,  \
                             g_0_zz_0_xyzz_1,   \
                             g_0_zz_0_xyzzz_0,  \
                             g_0_zz_0_xyzzz_1,  \
                             g_0_zz_0_xzzz_1,   \
                             g_0_zz_0_xzzzz_0,  \
                             g_0_zz_0_xzzzz_1,  \
                             g_0_zz_0_yyyy_1,   \
                             g_0_zz_0_yyyyy_0,  \
                             g_0_zz_0_yyyyy_1,  \
                             g_0_zz_0_yyyyz_0,  \
                             g_0_zz_0_yyyyz_1,  \
                             g_0_zz_0_yyyz_1,   \
                             g_0_zz_0_yyyzz_0,  \
                             g_0_zz_0_yyyzz_1,  \
                             g_0_zz_0_yyzz_1,   \
                             g_0_zz_0_yyzzz_0,  \
                             g_0_zz_0_yyzzz_1,  \
                             g_0_zz_0_yzzz_1,   \
                             g_0_zz_0_yzzzz_0,  \
                             g_0_zz_0_yzzzz_1,  \
                             g_0_zz_0_zzzz_1,   \
                             g_0_zz_0_zzzzz_0,  \
                             g_0_zz_0_zzzzz_1,  \
                             wp_x,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzz_0_xxxxx_0[i] = 5.0 * g_0_zz_0_xxxx_1[i] * fi_abcd_0 + g_0_zz_0_xxxxx_0[i] * pb_x + g_0_zz_0_xxxxx_1[i] * wp_x[i];

        g_0_xzz_0_xxxxy_0[i] = 4.0 * g_0_zz_0_xxxy_1[i] * fi_abcd_0 + g_0_zz_0_xxxxy_0[i] * pb_x + g_0_zz_0_xxxxy_1[i] * wp_x[i];

        g_0_xzz_0_xxxxz_0[i] = 4.0 * g_0_zz_0_xxxz_1[i] * fi_abcd_0 + g_0_zz_0_xxxxz_0[i] * pb_x + g_0_zz_0_xxxxz_1[i] * wp_x[i];

        g_0_xzz_0_xxxyy_0[i] = 3.0 * g_0_zz_0_xxyy_1[i] * fi_abcd_0 + g_0_zz_0_xxxyy_0[i] * pb_x + g_0_zz_0_xxxyy_1[i] * wp_x[i];

        g_0_xzz_0_xxxyz_0[i] = 3.0 * g_0_zz_0_xxyz_1[i] * fi_abcd_0 + g_0_zz_0_xxxyz_0[i] * pb_x + g_0_zz_0_xxxyz_1[i] * wp_x[i];

        g_0_xzz_0_xxxzz_0[i] = 3.0 * g_0_zz_0_xxzz_1[i] * fi_abcd_0 + g_0_zz_0_xxxzz_0[i] * pb_x + g_0_zz_0_xxxzz_1[i] * wp_x[i];

        g_0_xzz_0_xxyyy_0[i] = 2.0 * g_0_zz_0_xyyy_1[i] * fi_abcd_0 + g_0_zz_0_xxyyy_0[i] * pb_x + g_0_zz_0_xxyyy_1[i] * wp_x[i];

        g_0_xzz_0_xxyyz_0[i] = 2.0 * g_0_zz_0_xyyz_1[i] * fi_abcd_0 + g_0_zz_0_xxyyz_0[i] * pb_x + g_0_zz_0_xxyyz_1[i] * wp_x[i];

        g_0_xzz_0_xxyzz_0[i] = 2.0 * g_0_zz_0_xyzz_1[i] * fi_abcd_0 + g_0_zz_0_xxyzz_0[i] * pb_x + g_0_zz_0_xxyzz_1[i] * wp_x[i];

        g_0_xzz_0_xxzzz_0[i] = 2.0 * g_0_zz_0_xzzz_1[i] * fi_abcd_0 + g_0_zz_0_xxzzz_0[i] * pb_x + g_0_zz_0_xxzzz_1[i] * wp_x[i];

        g_0_xzz_0_xyyyy_0[i] = g_0_zz_0_yyyy_1[i] * fi_abcd_0 + g_0_zz_0_xyyyy_0[i] * pb_x + g_0_zz_0_xyyyy_1[i] * wp_x[i];

        g_0_xzz_0_xyyyz_0[i] = g_0_zz_0_yyyz_1[i] * fi_abcd_0 + g_0_zz_0_xyyyz_0[i] * pb_x + g_0_zz_0_xyyyz_1[i] * wp_x[i];

        g_0_xzz_0_xyyzz_0[i] = g_0_zz_0_yyzz_1[i] * fi_abcd_0 + g_0_zz_0_xyyzz_0[i] * pb_x + g_0_zz_0_xyyzz_1[i] * wp_x[i];

        g_0_xzz_0_xyzzz_0[i] = g_0_zz_0_yzzz_1[i] * fi_abcd_0 + g_0_zz_0_xyzzz_0[i] * pb_x + g_0_zz_0_xyzzz_1[i] * wp_x[i];

        g_0_xzz_0_xzzzz_0[i] = g_0_zz_0_zzzz_1[i] * fi_abcd_0 + g_0_zz_0_xzzzz_0[i] * pb_x + g_0_zz_0_xzzzz_1[i] * wp_x[i];

        g_0_xzz_0_yyyyy_0[i] = g_0_zz_0_yyyyy_0[i] * pb_x + g_0_zz_0_yyyyy_1[i] * wp_x[i];

        g_0_xzz_0_yyyyz_0[i] = g_0_zz_0_yyyyz_0[i] * pb_x + g_0_zz_0_yyyyz_1[i] * wp_x[i];

        g_0_xzz_0_yyyzz_0[i] = g_0_zz_0_yyyzz_0[i] * pb_x + g_0_zz_0_yyyzz_1[i] * wp_x[i];

        g_0_xzz_0_yyzzz_0[i] = g_0_zz_0_yyzzz_0[i] * pb_x + g_0_zz_0_yyzzz_1[i] * wp_x[i];

        g_0_xzz_0_yzzzz_0[i] = g_0_zz_0_yzzzz_0[i] * pb_x + g_0_zz_0_yzzzz_1[i] * wp_x[i];

        g_0_xzz_0_zzzzz_0[i] = g_0_zz_0_zzzzz_0[i] * pb_x + g_0_zz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 126-147 components of targeted buffer : SFSH

    auto g_0_yyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_sfsh + 126);

    auto g_0_yyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_sfsh + 127);

    auto g_0_yyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_sfsh + 128);

    auto g_0_yyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_sfsh + 129);

    auto g_0_yyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_sfsh + 130);

    auto g_0_yyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_sfsh + 131);

    auto g_0_yyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_sfsh + 132);

    auto g_0_yyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_sfsh + 133);

    auto g_0_yyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_sfsh + 134);

    auto g_0_yyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_sfsh + 135);

    auto g_0_yyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_sfsh + 136);

    auto g_0_yyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_sfsh + 137);

    auto g_0_yyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_sfsh + 138);

    auto g_0_yyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_sfsh + 139);

    auto g_0_yyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 140);

    auto g_0_yyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_sfsh + 141);

    auto g_0_yyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_sfsh + 142);

    auto g_0_yyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_sfsh + 143);

    auto g_0_yyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_sfsh + 144);

    auto g_0_yyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 145);

    auto g_0_yyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 146);

#pragma omp simd aligned(g_0_y_0_xxxxx_0,       \
                             g_0_y_0_xxxxx_1,   \
                             g_0_y_0_xxxxy_0,   \
                             g_0_y_0_xxxxy_1,   \
                             g_0_y_0_xxxxz_0,   \
                             g_0_y_0_xxxxz_1,   \
                             g_0_y_0_xxxyy_0,   \
                             g_0_y_0_xxxyy_1,   \
                             g_0_y_0_xxxyz_0,   \
                             g_0_y_0_xxxyz_1,   \
                             g_0_y_0_xxxzz_0,   \
                             g_0_y_0_xxxzz_1,   \
                             g_0_y_0_xxyyy_0,   \
                             g_0_y_0_xxyyy_1,   \
                             g_0_y_0_xxyyz_0,   \
                             g_0_y_0_xxyyz_1,   \
                             g_0_y_0_xxyzz_0,   \
                             g_0_y_0_xxyzz_1,   \
                             g_0_y_0_xxzzz_0,   \
                             g_0_y_0_xxzzz_1,   \
                             g_0_y_0_xyyyy_0,   \
                             g_0_y_0_xyyyy_1,   \
                             g_0_y_0_xyyyz_0,   \
                             g_0_y_0_xyyyz_1,   \
                             g_0_y_0_xyyzz_0,   \
                             g_0_y_0_xyyzz_1,   \
                             g_0_y_0_xyzzz_0,   \
                             g_0_y_0_xyzzz_1,   \
                             g_0_y_0_xzzzz_0,   \
                             g_0_y_0_xzzzz_1,   \
                             g_0_y_0_yyyyy_0,   \
                             g_0_y_0_yyyyy_1,   \
                             g_0_y_0_yyyyz_0,   \
                             g_0_y_0_yyyyz_1,   \
                             g_0_y_0_yyyzz_0,   \
                             g_0_y_0_yyyzz_1,   \
                             g_0_y_0_yyzzz_0,   \
                             g_0_y_0_yyzzz_1,   \
                             g_0_y_0_yzzzz_0,   \
                             g_0_y_0_yzzzz_1,   \
                             g_0_y_0_zzzzz_0,   \
                             g_0_y_0_zzzzz_1,   \
                             g_0_yy_0_xxxx_1,   \
                             g_0_yy_0_xxxxx_0,  \
                             g_0_yy_0_xxxxx_1,  \
                             g_0_yy_0_xxxxy_0,  \
                             g_0_yy_0_xxxxy_1,  \
                             g_0_yy_0_xxxxz_0,  \
                             g_0_yy_0_xxxxz_1,  \
                             g_0_yy_0_xxxy_1,   \
                             g_0_yy_0_xxxyy_0,  \
                             g_0_yy_0_xxxyy_1,  \
                             g_0_yy_0_xxxyz_0,  \
                             g_0_yy_0_xxxyz_1,  \
                             g_0_yy_0_xxxz_1,   \
                             g_0_yy_0_xxxzz_0,  \
                             g_0_yy_0_xxxzz_1,  \
                             g_0_yy_0_xxyy_1,   \
                             g_0_yy_0_xxyyy_0,  \
                             g_0_yy_0_xxyyy_1,  \
                             g_0_yy_0_xxyyz_0,  \
                             g_0_yy_0_xxyyz_1,  \
                             g_0_yy_0_xxyz_1,   \
                             g_0_yy_0_xxyzz_0,  \
                             g_0_yy_0_xxyzz_1,  \
                             g_0_yy_0_xxzz_1,   \
                             g_0_yy_0_xxzzz_0,  \
                             g_0_yy_0_xxzzz_1,  \
                             g_0_yy_0_xyyy_1,   \
                             g_0_yy_0_xyyyy_0,  \
                             g_0_yy_0_xyyyy_1,  \
                             g_0_yy_0_xyyyz_0,  \
                             g_0_yy_0_xyyyz_1,  \
                             g_0_yy_0_xyyz_1,   \
                             g_0_yy_0_xyyzz_0,  \
                             g_0_yy_0_xyyzz_1,  \
                             g_0_yy_0_xyzz_1,   \
                             g_0_yy_0_xyzzz_0,  \
                             g_0_yy_0_xyzzz_1,  \
                             g_0_yy_0_xzzz_1,   \
                             g_0_yy_0_xzzzz_0,  \
                             g_0_yy_0_xzzzz_1,  \
                             g_0_yy_0_yyyy_1,   \
                             g_0_yy_0_yyyyy_0,  \
                             g_0_yy_0_yyyyy_1,  \
                             g_0_yy_0_yyyyz_0,  \
                             g_0_yy_0_yyyyz_1,  \
                             g_0_yy_0_yyyz_1,   \
                             g_0_yy_0_yyyzz_0,  \
                             g_0_yy_0_yyyzz_1,  \
                             g_0_yy_0_yyzz_1,   \
                             g_0_yy_0_yyzzz_0,  \
                             g_0_yy_0_yyzzz_1,  \
                             g_0_yy_0_yzzz_1,   \
                             g_0_yy_0_yzzzz_0,  \
                             g_0_yy_0_yzzzz_1,  \
                             g_0_yy_0_zzzz_1,   \
                             g_0_yy_0_zzzzz_0,  \
                             g_0_yy_0_zzzzz_1,  \
                             g_0_yyy_0_xxxxx_0, \
                             g_0_yyy_0_xxxxy_0, \
                             g_0_yyy_0_xxxxz_0, \
                             g_0_yyy_0_xxxyy_0, \
                             g_0_yyy_0_xxxyz_0, \
                             g_0_yyy_0_xxxzz_0, \
                             g_0_yyy_0_xxyyy_0, \
                             g_0_yyy_0_xxyyz_0, \
                             g_0_yyy_0_xxyzz_0, \
                             g_0_yyy_0_xxzzz_0, \
                             g_0_yyy_0_xyyyy_0, \
                             g_0_yyy_0_xyyyz_0, \
                             g_0_yyy_0_xyyzz_0, \
                             g_0_yyy_0_xyzzz_0, \
                             g_0_yyy_0_xzzzz_0, \
                             g_0_yyy_0_yyyyy_0, \
                             g_0_yyy_0_yyyyz_0, \
                             g_0_yyy_0_yyyzz_0, \
                             g_0_yyy_0_yyzzz_0, \
                             g_0_yyy_0_yzzzz_0, \
                             g_0_yyy_0_zzzzz_0, \
                             wp_y,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyy_0_xxxxx_0[i] =
            2.0 * g_0_y_0_xxxxx_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxxx_1[i] * fti_ab_0 + g_0_yy_0_xxxxx_0[i] * pb_y + g_0_yy_0_xxxxx_1[i] * wp_y[i];

        g_0_yyy_0_xxxxy_0[i] = 2.0 * g_0_y_0_xxxxy_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxxy_1[i] * fti_ab_0 + g_0_yy_0_xxxx_1[i] * fi_abcd_0 +
                               g_0_yy_0_xxxxy_0[i] * pb_y + g_0_yy_0_xxxxy_1[i] * wp_y[i];

        g_0_yyy_0_xxxxz_0[i] =
            2.0 * g_0_y_0_xxxxz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxxz_1[i] * fti_ab_0 + g_0_yy_0_xxxxz_0[i] * pb_y + g_0_yy_0_xxxxz_1[i] * wp_y[i];

        g_0_yyy_0_xxxyy_0[i] = 2.0 * g_0_y_0_xxxyy_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxyy_1[i] * fti_ab_0 + 2.0 * g_0_yy_0_xxxy_1[i] * fi_abcd_0 +
                               g_0_yy_0_xxxyy_0[i] * pb_y + g_0_yy_0_xxxyy_1[i] * wp_y[i];

        g_0_yyy_0_xxxyz_0[i] = 2.0 * g_0_y_0_xxxyz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxyz_1[i] * fti_ab_0 + g_0_yy_0_xxxz_1[i] * fi_abcd_0 +
                               g_0_yy_0_xxxyz_0[i] * pb_y + g_0_yy_0_xxxyz_1[i] * wp_y[i];

        g_0_yyy_0_xxxzz_0[i] =
            2.0 * g_0_y_0_xxxzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxxzz_1[i] * fti_ab_0 + g_0_yy_0_xxxzz_0[i] * pb_y + g_0_yy_0_xxxzz_1[i] * wp_y[i];

        g_0_yyy_0_xxyyy_0[i] = 2.0 * g_0_y_0_xxyyy_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxyyy_1[i] * fti_ab_0 + 3.0 * g_0_yy_0_xxyy_1[i] * fi_abcd_0 +
                               g_0_yy_0_xxyyy_0[i] * pb_y + g_0_yy_0_xxyyy_1[i] * wp_y[i];

        g_0_yyy_0_xxyyz_0[i] = 2.0 * g_0_y_0_xxyyz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yy_0_xxyz_1[i] * fi_abcd_0 +
                               g_0_yy_0_xxyyz_0[i] * pb_y + g_0_yy_0_xxyyz_1[i] * wp_y[i];

        g_0_yyy_0_xxyzz_0[i] = 2.0 * g_0_y_0_xxyzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxyzz_1[i] * fti_ab_0 + g_0_yy_0_xxzz_1[i] * fi_abcd_0 +
                               g_0_yy_0_xxyzz_0[i] * pb_y + g_0_yy_0_xxyzz_1[i] * wp_y[i];

        g_0_yyy_0_xxzzz_0[i] =
            2.0 * g_0_y_0_xxzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxzzz_1[i] * fti_ab_0 + g_0_yy_0_xxzzz_0[i] * pb_y + g_0_yy_0_xxzzz_1[i] * wp_y[i];

        g_0_yyy_0_xyyyy_0[i] = 2.0 * g_0_y_0_xyyyy_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xyyyy_1[i] * fti_ab_0 + 4.0 * g_0_yy_0_xyyy_1[i] * fi_abcd_0 +
                               g_0_yy_0_xyyyy_0[i] * pb_y + g_0_yy_0_xyyyy_1[i] * wp_y[i];

        g_0_yyy_0_xyyyz_0[i] = 2.0 * g_0_y_0_xyyyz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yy_0_xyyz_1[i] * fi_abcd_0 +
                               g_0_yy_0_xyyyz_0[i] * pb_y + g_0_yy_0_xyyyz_1[i] * wp_y[i];

        g_0_yyy_0_xyyzz_0[i] = 2.0 * g_0_y_0_xyyzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yy_0_xyzz_1[i] * fi_abcd_0 +
                               g_0_yy_0_xyyzz_0[i] * pb_y + g_0_yy_0_xyyzz_1[i] * wp_y[i];

        g_0_yyy_0_xyzzz_0[i] = 2.0 * g_0_y_0_xyzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xyzzz_1[i] * fti_ab_0 + g_0_yy_0_xzzz_1[i] * fi_abcd_0 +
                               g_0_yy_0_xyzzz_0[i] * pb_y + g_0_yy_0_xyzzz_1[i] * wp_y[i];

        g_0_yyy_0_xzzzz_0[i] =
            2.0 * g_0_y_0_xzzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xzzzz_1[i] * fti_ab_0 + g_0_yy_0_xzzzz_0[i] * pb_y + g_0_yy_0_xzzzz_1[i] * wp_y[i];

        g_0_yyy_0_yyyyy_0[i] = 2.0 * g_0_y_0_yyyyy_0[i] * fi_ab_0 - 2.0 * g_0_y_0_yyyyy_1[i] * fti_ab_0 + 5.0 * g_0_yy_0_yyyy_1[i] * fi_abcd_0 +
                               g_0_yy_0_yyyyy_0[i] * pb_y + g_0_yy_0_yyyyy_1[i] * wp_y[i];

        g_0_yyy_0_yyyyz_0[i] = 2.0 * g_0_y_0_yyyyz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_yyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yy_0_yyyz_1[i] * fi_abcd_0 +
                               g_0_yy_0_yyyyz_0[i] * pb_y + g_0_yy_0_yyyyz_1[i] * wp_y[i];

        g_0_yyy_0_yyyzz_0[i] = 2.0 * g_0_y_0_yyyzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_yyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yy_0_yyzz_1[i] * fi_abcd_0 +
                               g_0_yy_0_yyyzz_0[i] * pb_y + g_0_yy_0_yyyzz_1[i] * wp_y[i];

        g_0_yyy_0_yyzzz_0[i] = 2.0 * g_0_y_0_yyzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_yyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yy_0_yzzz_1[i] * fi_abcd_0 +
                               g_0_yy_0_yyzzz_0[i] * pb_y + g_0_yy_0_yyzzz_1[i] * wp_y[i];

        g_0_yyy_0_yzzzz_0[i] = 2.0 * g_0_y_0_yzzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_yzzzz_1[i] * fti_ab_0 + g_0_yy_0_zzzz_1[i] * fi_abcd_0 +
                               g_0_yy_0_yzzzz_0[i] * pb_y + g_0_yy_0_yzzzz_1[i] * wp_y[i];

        g_0_yyy_0_zzzzz_0[i] =
            2.0 * g_0_y_0_zzzzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_zzzzz_1[i] * fti_ab_0 + g_0_yy_0_zzzzz_0[i] * pb_y + g_0_yy_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 147-168 components of targeted buffer : SFSH

    auto g_0_yyz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sfsh + 147);

    auto g_0_yyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sfsh + 148);

    auto g_0_yyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sfsh + 149);

    auto g_0_yyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sfsh + 150);

    auto g_0_yyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sfsh + 151);

    auto g_0_yyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sfsh + 152);

    auto g_0_yyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sfsh + 153);

    auto g_0_yyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sfsh + 154);

    auto g_0_yyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sfsh + 155);

    auto g_0_yyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sfsh + 156);

    auto g_0_yyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sfsh + 157);

    auto g_0_yyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sfsh + 158);

    auto g_0_yyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sfsh + 159);

    auto g_0_yyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sfsh + 160);

    auto g_0_yyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 161);

    auto g_0_yyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sfsh + 162);

    auto g_0_yyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sfsh + 163);

    auto g_0_yyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sfsh + 164);

    auto g_0_yyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sfsh + 165);

    auto g_0_yyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 166);

    auto g_0_yyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 167);

#pragma omp simd aligned(g_0_yy_0_xxxx_1,       \
                             g_0_yy_0_xxxxx_0,  \
                             g_0_yy_0_xxxxx_1,  \
                             g_0_yy_0_xxxxy_0,  \
                             g_0_yy_0_xxxxy_1,  \
                             g_0_yy_0_xxxxz_0,  \
                             g_0_yy_0_xxxxz_1,  \
                             g_0_yy_0_xxxy_1,   \
                             g_0_yy_0_xxxyy_0,  \
                             g_0_yy_0_xxxyy_1,  \
                             g_0_yy_0_xxxyz_0,  \
                             g_0_yy_0_xxxyz_1,  \
                             g_0_yy_0_xxxz_1,   \
                             g_0_yy_0_xxxzz_0,  \
                             g_0_yy_0_xxxzz_1,  \
                             g_0_yy_0_xxyy_1,   \
                             g_0_yy_0_xxyyy_0,  \
                             g_0_yy_0_xxyyy_1,  \
                             g_0_yy_0_xxyyz_0,  \
                             g_0_yy_0_xxyyz_1,  \
                             g_0_yy_0_xxyz_1,   \
                             g_0_yy_0_xxyzz_0,  \
                             g_0_yy_0_xxyzz_1,  \
                             g_0_yy_0_xxzz_1,   \
                             g_0_yy_0_xxzzz_0,  \
                             g_0_yy_0_xxzzz_1,  \
                             g_0_yy_0_xyyy_1,   \
                             g_0_yy_0_xyyyy_0,  \
                             g_0_yy_0_xyyyy_1,  \
                             g_0_yy_0_xyyyz_0,  \
                             g_0_yy_0_xyyyz_1,  \
                             g_0_yy_0_xyyz_1,   \
                             g_0_yy_0_xyyzz_0,  \
                             g_0_yy_0_xyyzz_1,  \
                             g_0_yy_0_xyzz_1,   \
                             g_0_yy_0_xyzzz_0,  \
                             g_0_yy_0_xyzzz_1,  \
                             g_0_yy_0_xzzz_1,   \
                             g_0_yy_0_xzzzz_0,  \
                             g_0_yy_0_xzzzz_1,  \
                             g_0_yy_0_yyyy_1,   \
                             g_0_yy_0_yyyyy_0,  \
                             g_0_yy_0_yyyyy_1,  \
                             g_0_yy_0_yyyyz_0,  \
                             g_0_yy_0_yyyyz_1,  \
                             g_0_yy_0_yyyz_1,   \
                             g_0_yy_0_yyyzz_0,  \
                             g_0_yy_0_yyyzz_1,  \
                             g_0_yy_0_yyzz_1,   \
                             g_0_yy_0_yyzzz_0,  \
                             g_0_yy_0_yyzzz_1,  \
                             g_0_yy_0_yzzz_1,   \
                             g_0_yy_0_yzzzz_0,  \
                             g_0_yy_0_yzzzz_1,  \
                             g_0_yy_0_zzzz_1,   \
                             g_0_yy_0_zzzzz_0,  \
                             g_0_yy_0_zzzzz_1,  \
                             g_0_yyz_0_xxxxx_0, \
                             g_0_yyz_0_xxxxy_0, \
                             g_0_yyz_0_xxxxz_0, \
                             g_0_yyz_0_xxxyy_0, \
                             g_0_yyz_0_xxxyz_0, \
                             g_0_yyz_0_xxxzz_0, \
                             g_0_yyz_0_xxyyy_0, \
                             g_0_yyz_0_xxyyz_0, \
                             g_0_yyz_0_xxyzz_0, \
                             g_0_yyz_0_xxzzz_0, \
                             g_0_yyz_0_xyyyy_0, \
                             g_0_yyz_0_xyyyz_0, \
                             g_0_yyz_0_xyyzz_0, \
                             g_0_yyz_0_xyzzz_0, \
                             g_0_yyz_0_xzzzz_0, \
                             g_0_yyz_0_yyyyy_0, \
                             g_0_yyz_0_yyyyz_0, \
                             g_0_yyz_0_yyyzz_0, \
                             g_0_yyz_0_yyzzz_0, \
                             g_0_yyz_0_yzzzz_0, \
                             g_0_yyz_0_zzzzz_0, \
                             wp_z,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyz_0_xxxxx_0[i] = g_0_yy_0_xxxxx_0[i] * pb_z + g_0_yy_0_xxxxx_1[i] * wp_z[i];

        g_0_yyz_0_xxxxy_0[i] = g_0_yy_0_xxxxy_0[i] * pb_z + g_0_yy_0_xxxxy_1[i] * wp_z[i];

        g_0_yyz_0_xxxxz_0[i] = g_0_yy_0_xxxx_1[i] * fi_abcd_0 + g_0_yy_0_xxxxz_0[i] * pb_z + g_0_yy_0_xxxxz_1[i] * wp_z[i];

        g_0_yyz_0_xxxyy_0[i] = g_0_yy_0_xxxyy_0[i] * pb_z + g_0_yy_0_xxxyy_1[i] * wp_z[i];

        g_0_yyz_0_xxxyz_0[i] = g_0_yy_0_xxxy_1[i] * fi_abcd_0 + g_0_yy_0_xxxyz_0[i] * pb_z + g_0_yy_0_xxxyz_1[i] * wp_z[i];

        g_0_yyz_0_xxxzz_0[i] = 2.0 * g_0_yy_0_xxxz_1[i] * fi_abcd_0 + g_0_yy_0_xxxzz_0[i] * pb_z + g_0_yy_0_xxxzz_1[i] * wp_z[i];

        g_0_yyz_0_xxyyy_0[i] = g_0_yy_0_xxyyy_0[i] * pb_z + g_0_yy_0_xxyyy_1[i] * wp_z[i];

        g_0_yyz_0_xxyyz_0[i] = g_0_yy_0_xxyy_1[i] * fi_abcd_0 + g_0_yy_0_xxyyz_0[i] * pb_z + g_0_yy_0_xxyyz_1[i] * wp_z[i];

        g_0_yyz_0_xxyzz_0[i] = 2.0 * g_0_yy_0_xxyz_1[i] * fi_abcd_0 + g_0_yy_0_xxyzz_0[i] * pb_z + g_0_yy_0_xxyzz_1[i] * wp_z[i];

        g_0_yyz_0_xxzzz_0[i] = 3.0 * g_0_yy_0_xxzz_1[i] * fi_abcd_0 + g_0_yy_0_xxzzz_0[i] * pb_z + g_0_yy_0_xxzzz_1[i] * wp_z[i];

        g_0_yyz_0_xyyyy_0[i] = g_0_yy_0_xyyyy_0[i] * pb_z + g_0_yy_0_xyyyy_1[i] * wp_z[i];

        g_0_yyz_0_xyyyz_0[i] = g_0_yy_0_xyyy_1[i] * fi_abcd_0 + g_0_yy_0_xyyyz_0[i] * pb_z + g_0_yy_0_xyyyz_1[i] * wp_z[i];

        g_0_yyz_0_xyyzz_0[i] = 2.0 * g_0_yy_0_xyyz_1[i] * fi_abcd_0 + g_0_yy_0_xyyzz_0[i] * pb_z + g_0_yy_0_xyyzz_1[i] * wp_z[i];

        g_0_yyz_0_xyzzz_0[i] = 3.0 * g_0_yy_0_xyzz_1[i] * fi_abcd_0 + g_0_yy_0_xyzzz_0[i] * pb_z + g_0_yy_0_xyzzz_1[i] * wp_z[i];

        g_0_yyz_0_xzzzz_0[i] = 4.0 * g_0_yy_0_xzzz_1[i] * fi_abcd_0 + g_0_yy_0_xzzzz_0[i] * pb_z + g_0_yy_0_xzzzz_1[i] * wp_z[i];

        g_0_yyz_0_yyyyy_0[i] = g_0_yy_0_yyyyy_0[i] * pb_z + g_0_yy_0_yyyyy_1[i] * wp_z[i];

        g_0_yyz_0_yyyyz_0[i] = g_0_yy_0_yyyy_1[i] * fi_abcd_0 + g_0_yy_0_yyyyz_0[i] * pb_z + g_0_yy_0_yyyyz_1[i] * wp_z[i];

        g_0_yyz_0_yyyzz_0[i] = 2.0 * g_0_yy_0_yyyz_1[i] * fi_abcd_0 + g_0_yy_0_yyyzz_0[i] * pb_z + g_0_yy_0_yyyzz_1[i] * wp_z[i];

        g_0_yyz_0_yyzzz_0[i] = 3.0 * g_0_yy_0_yyzz_1[i] * fi_abcd_0 + g_0_yy_0_yyzzz_0[i] * pb_z + g_0_yy_0_yyzzz_1[i] * wp_z[i];

        g_0_yyz_0_yzzzz_0[i] = 4.0 * g_0_yy_0_yzzz_1[i] * fi_abcd_0 + g_0_yy_0_yzzzz_0[i] * pb_z + g_0_yy_0_yzzzz_1[i] * wp_z[i];

        g_0_yyz_0_zzzzz_0[i] = 5.0 * g_0_yy_0_zzzz_1[i] * fi_abcd_0 + g_0_yy_0_zzzzz_0[i] * pb_z + g_0_yy_0_zzzzz_1[i] * wp_z[i];
    }

    /// Set up 168-189 components of targeted buffer : SFSH

    auto g_0_yzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sfsh + 168);

    auto g_0_yzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sfsh + 169);

    auto g_0_yzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sfsh + 170);

    auto g_0_yzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sfsh + 171);

    auto g_0_yzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sfsh + 172);

    auto g_0_yzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sfsh + 173);

    auto g_0_yzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sfsh + 174);

    auto g_0_yzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sfsh + 175);

    auto g_0_yzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sfsh + 176);

    auto g_0_yzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sfsh + 177);

    auto g_0_yzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sfsh + 178);

    auto g_0_yzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sfsh + 179);

    auto g_0_yzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sfsh + 180);

    auto g_0_yzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sfsh + 181);

    auto g_0_yzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 182);

    auto g_0_yzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sfsh + 183);

    auto g_0_yzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sfsh + 184);

    auto g_0_yzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sfsh + 185);

    auto g_0_yzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sfsh + 186);

    auto g_0_yzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 187);

    auto g_0_yzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 188);

#pragma omp simd aligned(g_0_yzz_0_xxxxx_0,     \
                             g_0_yzz_0_xxxxy_0, \
                             g_0_yzz_0_xxxxz_0, \
                             g_0_yzz_0_xxxyy_0, \
                             g_0_yzz_0_xxxyz_0, \
                             g_0_yzz_0_xxxzz_0, \
                             g_0_yzz_0_xxyyy_0, \
                             g_0_yzz_0_xxyyz_0, \
                             g_0_yzz_0_xxyzz_0, \
                             g_0_yzz_0_xxzzz_0, \
                             g_0_yzz_0_xyyyy_0, \
                             g_0_yzz_0_xyyyz_0, \
                             g_0_yzz_0_xyyzz_0, \
                             g_0_yzz_0_xyzzz_0, \
                             g_0_yzz_0_xzzzz_0, \
                             g_0_yzz_0_yyyyy_0, \
                             g_0_yzz_0_yyyyz_0, \
                             g_0_yzz_0_yyyzz_0, \
                             g_0_yzz_0_yyzzz_0, \
                             g_0_yzz_0_yzzzz_0, \
                             g_0_yzz_0_zzzzz_0, \
                             g_0_zz_0_xxxx_1,   \
                             g_0_zz_0_xxxxx_0,  \
                             g_0_zz_0_xxxxx_1,  \
                             g_0_zz_0_xxxxy_0,  \
                             g_0_zz_0_xxxxy_1,  \
                             g_0_zz_0_xxxxz_0,  \
                             g_0_zz_0_xxxxz_1,  \
                             g_0_zz_0_xxxy_1,   \
                             g_0_zz_0_xxxyy_0,  \
                             g_0_zz_0_xxxyy_1,  \
                             g_0_zz_0_xxxyz_0,  \
                             g_0_zz_0_xxxyz_1,  \
                             g_0_zz_0_xxxz_1,   \
                             g_0_zz_0_xxxzz_0,  \
                             g_0_zz_0_xxxzz_1,  \
                             g_0_zz_0_xxyy_1,   \
                             g_0_zz_0_xxyyy_0,  \
                             g_0_zz_0_xxyyy_1,  \
                             g_0_zz_0_xxyyz_0,  \
                             g_0_zz_0_xxyyz_1,  \
                             g_0_zz_0_xxyz_1,   \
                             g_0_zz_0_xxyzz_0,  \
                             g_0_zz_0_xxyzz_1,  \
                             g_0_zz_0_xxzz_1,   \
                             g_0_zz_0_xxzzz_0,  \
                             g_0_zz_0_xxzzz_1,  \
                             g_0_zz_0_xyyy_1,   \
                             g_0_zz_0_xyyyy_0,  \
                             g_0_zz_0_xyyyy_1,  \
                             g_0_zz_0_xyyyz_0,  \
                             g_0_zz_0_xyyyz_1,  \
                             g_0_zz_0_xyyz_1,   \
                             g_0_zz_0_xyyzz_0,  \
                             g_0_zz_0_xyyzz_1,  \
                             g_0_zz_0_xyzz_1,   \
                             g_0_zz_0_xyzzz_0,  \
                             g_0_zz_0_xyzzz_1,  \
                             g_0_zz_0_xzzz_1,   \
                             g_0_zz_0_xzzzz_0,  \
                             g_0_zz_0_xzzzz_1,  \
                             g_0_zz_0_yyyy_1,   \
                             g_0_zz_0_yyyyy_0,  \
                             g_0_zz_0_yyyyy_1,  \
                             g_0_zz_0_yyyyz_0,  \
                             g_0_zz_0_yyyyz_1,  \
                             g_0_zz_0_yyyz_1,   \
                             g_0_zz_0_yyyzz_0,  \
                             g_0_zz_0_yyyzz_1,  \
                             g_0_zz_0_yyzz_1,   \
                             g_0_zz_0_yyzzz_0,  \
                             g_0_zz_0_yyzzz_1,  \
                             g_0_zz_0_yzzz_1,   \
                             g_0_zz_0_yzzzz_0,  \
                             g_0_zz_0_yzzzz_1,  \
                             g_0_zz_0_zzzz_1,   \
                             g_0_zz_0_zzzzz_0,  \
                             g_0_zz_0_zzzzz_1,  \
                             wp_y,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzz_0_xxxxx_0[i] = g_0_zz_0_xxxxx_0[i] * pb_y + g_0_zz_0_xxxxx_1[i] * wp_y[i];

        g_0_yzz_0_xxxxy_0[i] = g_0_zz_0_xxxx_1[i] * fi_abcd_0 + g_0_zz_0_xxxxy_0[i] * pb_y + g_0_zz_0_xxxxy_1[i] * wp_y[i];

        g_0_yzz_0_xxxxz_0[i] = g_0_zz_0_xxxxz_0[i] * pb_y + g_0_zz_0_xxxxz_1[i] * wp_y[i];

        g_0_yzz_0_xxxyy_0[i] = 2.0 * g_0_zz_0_xxxy_1[i] * fi_abcd_0 + g_0_zz_0_xxxyy_0[i] * pb_y + g_0_zz_0_xxxyy_1[i] * wp_y[i];

        g_0_yzz_0_xxxyz_0[i] = g_0_zz_0_xxxz_1[i] * fi_abcd_0 + g_0_zz_0_xxxyz_0[i] * pb_y + g_0_zz_0_xxxyz_1[i] * wp_y[i];

        g_0_yzz_0_xxxzz_0[i] = g_0_zz_0_xxxzz_0[i] * pb_y + g_0_zz_0_xxxzz_1[i] * wp_y[i];

        g_0_yzz_0_xxyyy_0[i] = 3.0 * g_0_zz_0_xxyy_1[i] * fi_abcd_0 + g_0_zz_0_xxyyy_0[i] * pb_y + g_0_zz_0_xxyyy_1[i] * wp_y[i];

        g_0_yzz_0_xxyyz_0[i] = 2.0 * g_0_zz_0_xxyz_1[i] * fi_abcd_0 + g_0_zz_0_xxyyz_0[i] * pb_y + g_0_zz_0_xxyyz_1[i] * wp_y[i];

        g_0_yzz_0_xxyzz_0[i] = g_0_zz_0_xxzz_1[i] * fi_abcd_0 + g_0_zz_0_xxyzz_0[i] * pb_y + g_0_zz_0_xxyzz_1[i] * wp_y[i];

        g_0_yzz_0_xxzzz_0[i] = g_0_zz_0_xxzzz_0[i] * pb_y + g_0_zz_0_xxzzz_1[i] * wp_y[i];

        g_0_yzz_0_xyyyy_0[i] = 4.0 * g_0_zz_0_xyyy_1[i] * fi_abcd_0 + g_0_zz_0_xyyyy_0[i] * pb_y + g_0_zz_0_xyyyy_1[i] * wp_y[i];

        g_0_yzz_0_xyyyz_0[i] = 3.0 * g_0_zz_0_xyyz_1[i] * fi_abcd_0 + g_0_zz_0_xyyyz_0[i] * pb_y + g_0_zz_0_xyyyz_1[i] * wp_y[i];

        g_0_yzz_0_xyyzz_0[i] = 2.0 * g_0_zz_0_xyzz_1[i] * fi_abcd_0 + g_0_zz_0_xyyzz_0[i] * pb_y + g_0_zz_0_xyyzz_1[i] * wp_y[i];

        g_0_yzz_0_xyzzz_0[i] = g_0_zz_0_xzzz_1[i] * fi_abcd_0 + g_0_zz_0_xyzzz_0[i] * pb_y + g_0_zz_0_xyzzz_1[i] * wp_y[i];

        g_0_yzz_0_xzzzz_0[i] = g_0_zz_0_xzzzz_0[i] * pb_y + g_0_zz_0_xzzzz_1[i] * wp_y[i];

        g_0_yzz_0_yyyyy_0[i] = 5.0 * g_0_zz_0_yyyy_1[i] * fi_abcd_0 + g_0_zz_0_yyyyy_0[i] * pb_y + g_0_zz_0_yyyyy_1[i] * wp_y[i];

        g_0_yzz_0_yyyyz_0[i] = 4.0 * g_0_zz_0_yyyz_1[i] * fi_abcd_0 + g_0_zz_0_yyyyz_0[i] * pb_y + g_0_zz_0_yyyyz_1[i] * wp_y[i];

        g_0_yzz_0_yyyzz_0[i] = 3.0 * g_0_zz_0_yyzz_1[i] * fi_abcd_0 + g_0_zz_0_yyyzz_0[i] * pb_y + g_0_zz_0_yyyzz_1[i] * wp_y[i];

        g_0_yzz_0_yyzzz_0[i] = 2.0 * g_0_zz_0_yzzz_1[i] * fi_abcd_0 + g_0_zz_0_yyzzz_0[i] * pb_y + g_0_zz_0_yyzzz_1[i] * wp_y[i];

        g_0_yzz_0_yzzzz_0[i] = g_0_zz_0_zzzz_1[i] * fi_abcd_0 + g_0_zz_0_yzzzz_0[i] * pb_y + g_0_zz_0_yzzzz_1[i] * wp_y[i];

        g_0_yzz_0_zzzzz_0[i] = g_0_zz_0_zzzzz_0[i] * pb_y + g_0_zz_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 189-210 components of targeted buffer : SFSH

    auto g_0_zzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sfsh + 189);

    auto g_0_zzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sfsh + 190);

    auto g_0_zzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sfsh + 191);

    auto g_0_zzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sfsh + 192);

    auto g_0_zzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sfsh + 193);

    auto g_0_zzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sfsh + 194);

    auto g_0_zzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sfsh + 195);

    auto g_0_zzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sfsh + 196);

    auto g_0_zzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sfsh + 197);

    auto g_0_zzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sfsh + 198);

    auto g_0_zzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sfsh + 199);

    auto g_0_zzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sfsh + 200);

    auto g_0_zzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sfsh + 201);

    auto g_0_zzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sfsh + 202);

    auto g_0_zzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 203);

    auto g_0_zzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sfsh + 204);

    auto g_0_zzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sfsh + 205);

    auto g_0_zzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sfsh + 206);

    auto g_0_zzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sfsh + 207);

    auto g_0_zzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 208);

    auto g_0_zzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 209);

#pragma omp simd aligned(g_0_z_0_xxxxx_0,       \
                             g_0_z_0_xxxxx_1,   \
                             g_0_z_0_xxxxy_0,   \
                             g_0_z_0_xxxxy_1,   \
                             g_0_z_0_xxxxz_0,   \
                             g_0_z_0_xxxxz_1,   \
                             g_0_z_0_xxxyy_0,   \
                             g_0_z_0_xxxyy_1,   \
                             g_0_z_0_xxxyz_0,   \
                             g_0_z_0_xxxyz_1,   \
                             g_0_z_0_xxxzz_0,   \
                             g_0_z_0_xxxzz_1,   \
                             g_0_z_0_xxyyy_0,   \
                             g_0_z_0_xxyyy_1,   \
                             g_0_z_0_xxyyz_0,   \
                             g_0_z_0_xxyyz_1,   \
                             g_0_z_0_xxyzz_0,   \
                             g_0_z_0_xxyzz_1,   \
                             g_0_z_0_xxzzz_0,   \
                             g_0_z_0_xxzzz_1,   \
                             g_0_z_0_xyyyy_0,   \
                             g_0_z_0_xyyyy_1,   \
                             g_0_z_0_xyyyz_0,   \
                             g_0_z_0_xyyyz_1,   \
                             g_0_z_0_xyyzz_0,   \
                             g_0_z_0_xyyzz_1,   \
                             g_0_z_0_xyzzz_0,   \
                             g_0_z_0_xyzzz_1,   \
                             g_0_z_0_xzzzz_0,   \
                             g_0_z_0_xzzzz_1,   \
                             g_0_z_0_yyyyy_0,   \
                             g_0_z_0_yyyyy_1,   \
                             g_0_z_0_yyyyz_0,   \
                             g_0_z_0_yyyyz_1,   \
                             g_0_z_0_yyyzz_0,   \
                             g_0_z_0_yyyzz_1,   \
                             g_0_z_0_yyzzz_0,   \
                             g_0_z_0_yyzzz_1,   \
                             g_0_z_0_yzzzz_0,   \
                             g_0_z_0_yzzzz_1,   \
                             g_0_z_0_zzzzz_0,   \
                             g_0_z_0_zzzzz_1,   \
                             g_0_zz_0_xxxx_1,   \
                             g_0_zz_0_xxxxx_0,  \
                             g_0_zz_0_xxxxx_1,  \
                             g_0_zz_0_xxxxy_0,  \
                             g_0_zz_0_xxxxy_1,  \
                             g_0_zz_0_xxxxz_0,  \
                             g_0_zz_0_xxxxz_1,  \
                             g_0_zz_0_xxxy_1,   \
                             g_0_zz_0_xxxyy_0,  \
                             g_0_zz_0_xxxyy_1,  \
                             g_0_zz_0_xxxyz_0,  \
                             g_0_zz_0_xxxyz_1,  \
                             g_0_zz_0_xxxz_1,   \
                             g_0_zz_0_xxxzz_0,  \
                             g_0_zz_0_xxxzz_1,  \
                             g_0_zz_0_xxyy_1,   \
                             g_0_zz_0_xxyyy_0,  \
                             g_0_zz_0_xxyyy_1,  \
                             g_0_zz_0_xxyyz_0,  \
                             g_0_zz_0_xxyyz_1,  \
                             g_0_zz_0_xxyz_1,   \
                             g_0_zz_0_xxyzz_0,  \
                             g_0_zz_0_xxyzz_1,  \
                             g_0_zz_0_xxzz_1,   \
                             g_0_zz_0_xxzzz_0,  \
                             g_0_zz_0_xxzzz_1,  \
                             g_0_zz_0_xyyy_1,   \
                             g_0_zz_0_xyyyy_0,  \
                             g_0_zz_0_xyyyy_1,  \
                             g_0_zz_0_xyyyz_0,  \
                             g_0_zz_0_xyyyz_1,  \
                             g_0_zz_0_xyyz_1,   \
                             g_0_zz_0_xyyzz_0,  \
                             g_0_zz_0_xyyzz_1,  \
                             g_0_zz_0_xyzz_1,   \
                             g_0_zz_0_xyzzz_0,  \
                             g_0_zz_0_xyzzz_1,  \
                             g_0_zz_0_xzzz_1,   \
                             g_0_zz_0_xzzzz_0,  \
                             g_0_zz_0_xzzzz_1,  \
                             g_0_zz_0_yyyy_1,   \
                             g_0_zz_0_yyyyy_0,  \
                             g_0_zz_0_yyyyy_1,  \
                             g_0_zz_0_yyyyz_0,  \
                             g_0_zz_0_yyyyz_1,  \
                             g_0_zz_0_yyyz_1,   \
                             g_0_zz_0_yyyzz_0,  \
                             g_0_zz_0_yyyzz_1,  \
                             g_0_zz_0_yyzz_1,   \
                             g_0_zz_0_yyzzz_0,  \
                             g_0_zz_0_yyzzz_1,  \
                             g_0_zz_0_yzzz_1,   \
                             g_0_zz_0_yzzzz_0,  \
                             g_0_zz_0_yzzzz_1,  \
                             g_0_zz_0_zzzz_1,   \
                             g_0_zz_0_zzzzz_0,  \
                             g_0_zz_0_zzzzz_1,  \
                             g_0_zzz_0_xxxxx_0, \
                             g_0_zzz_0_xxxxy_0, \
                             g_0_zzz_0_xxxxz_0, \
                             g_0_zzz_0_xxxyy_0, \
                             g_0_zzz_0_xxxyz_0, \
                             g_0_zzz_0_xxxzz_0, \
                             g_0_zzz_0_xxyyy_0, \
                             g_0_zzz_0_xxyyz_0, \
                             g_0_zzz_0_xxyzz_0, \
                             g_0_zzz_0_xxzzz_0, \
                             g_0_zzz_0_xyyyy_0, \
                             g_0_zzz_0_xyyyz_0, \
                             g_0_zzz_0_xyyzz_0, \
                             g_0_zzz_0_xyzzz_0, \
                             g_0_zzz_0_xzzzz_0, \
                             g_0_zzz_0_yyyyy_0, \
                             g_0_zzz_0_yyyyz_0, \
                             g_0_zzz_0_yyyzz_0, \
                             g_0_zzz_0_yyzzz_0, \
                             g_0_zzz_0_yzzzz_0, \
                             g_0_zzz_0_zzzzz_0, \
                             wp_z,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzz_0_xxxxx_0[i] =
            2.0 * g_0_z_0_xxxxx_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxxx_1[i] * fti_ab_0 + g_0_zz_0_xxxxx_0[i] * pb_z + g_0_zz_0_xxxxx_1[i] * wp_z[i];

        g_0_zzz_0_xxxxy_0[i] =
            2.0 * g_0_z_0_xxxxy_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxxy_1[i] * fti_ab_0 + g_0_zz_0_xxxxy_0[i] * pb_z + g_0_zz_0_xxxxy_1[i] * wp_z[i];

        g_0_zzz_0_xxxxz_0[i] = 2.0 * g_0_z_0_xxxxz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxxz_1[i] * fti_ab_0 + g_0_zz_0_xxxx_1[i] * fi_abcd_0 +
                               g_0_zz_0_xxxxz_0[i] * pb_z + g_0_zz_0_xxxxz_1[i] * wp_z[i];

        g_0_zzz_0_xxxyy_0[i] =
            2.0 * g_0_z_0_xxxyy_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxyy_1[i] * fti_ab_0 + g_0_zz_0_xxxyy_0[i] * pb_z + g_0_zz_0_xxxyy_1[i] * wp_z[i];

        g_0_zzz_0_xxxyz_0[i] = 2.0 * g_0_z_0_xxxyz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxyz_1[i] * fti_ab_0 + g_0_zz_0_xxxy_1[i] * fi_abcd_0 +
                               g_0_zz_0_xxxyz_0[i] * pb_z + g_0_zz_0_xxxyz_1[i] * wp_z[i];

        g_0_zzz_0_xxxzz_0[i] = 2.0 * g_0_z_0_xxxzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxxzz_1[i] * fti_ab_0 + 2.0 * g_0_zz_0_xxxz_1[i] * fi_abcd_0 +
                               g_0_zz_0_xxxzz_0[i] * pb_z + g_0_zz_0_xxxzz_1[i] * wp_z[i];

        g_0_zzz_0_xxyyy_0[i] =
            2.0 * g_0_z_0_xxyyy_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxyyy_1[i] * fti_ab_0 + g_0_zz_0_xxyyy_0[i] * pb_z + g_0_zz_0_xxyyy_1[i] * wp_z[i];

        g_0_zzz_0_xxyyz_0[i] = 2.0 * g_0_z_0_xxyyz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxyyz_1[i] * fti_ab_0 + g_0_zz_0_xxyy_1[i] * fi_abcd_0 +
                               g_0_zz_0_xxyyz_0[i] * pb_z + g_0_zz_0_xxyyz_1[i] * wp_z[i];

        g_0_zzz_0_xxyzz_0[i] = 2.0 * g_0_z_0_xxyzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_zz_0_xxyz_1[i] * fi_abcd_0 +
                               g_0_zz_0_xxyzz_0[i] * pb_z + g_0_zz_0_xxyzz_1[i] * wp_z[i];

        g_0_zzz_0_xxzzz_0[i] = 2.0 * g_0_z_0_xxzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxzzz_1[i] * fti_ab_0 + 3.0 * g_0_zz_0_xxzz_1[i] * fi_abcd_0 +
                               g_0_zz_0_xxzzz_0[i] * pb_z + g_0_zz_0_xxzzz_1[i] * wp_z[i];

        g_0_zzz_0_xyyyy_0[i] =
            2.0 * g_0_z_0_xyyyy_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xyyyy_1[i] * fti_ab_0 + g_0_zz_0_xyyyy_0[i] * pb_z + g_0_zz_0_xyyyy_1[i] * wp_z[i];

        g_0_zzz_0_xyyyz_0[i] = 2.0 * g_0_z_0_xyyyz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xyyyz_1[i] * fti_ab_0 + g_0_zz_0_xyyy_1[i] * fi_abcd_0 +
                               g_0_zz_0_xyyyz_0[i] * pb_z + g_0_zz_0_xyyyz_1[i] * wp_z[i];

        g_0_zzz_0_xyyzz_0[i] = 2.0 * g_0_z_0_xyyzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zz_0_xyyz_1[i] * fi_abcd_0 +
                               g_0_zz_0_xyyzz_0[i] * pb_z + g_0_zz_0_xyyzz_1[i] * wp_z[i];

        g_0_zzz_0_xyzzz_0[i] = 2.0 * g_0_z_0_xyzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zz_0_xyzz_1[i] * fi_abcd_0 +
                               g_0_zz_0_xyzzz_0[i] * pb_z + g_0_zz_0_xyzzz_1[i] * wp_z[i];

        g_0_zzz_0_xzzzz_0[i] = 2.0 * g_0_z_0_xzzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zz_0_xzzz_1[i] * fi_abcd_0 +
                               g_0_zz_0_xzzzz_0[i] * pb_z + g_0_zz_0_xzzzz_1[i] * wp_z[i];

        g_0_zzz_0_yyyyy_0[i] =
            2.0 * g_0_z_0_yyyyy_0[i] * fi_ab_0 - 2.0 * g_0_z_0_yyyyy_1[i] * fti_ab_0 + g_0_zz_0_yyyyy_0[i] * pb_z + g_0_zz_0_yyyyy_1[i] * wp_z[i];

        g_0_zzz_0_yyyyz_0[i] = 2.0 * g_0_z_0_yyyyz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_yyyyz_1[i] * fti_ab_0 + g_0_zz_0_yyyy_1[i] * fi_abcd_0 +
                               g_0_zz_0_yyyyz_0[i] * pb_z + g_0_zz_0_yyyyz_1[i] * wp_z[i];

        g_0_zzz_0_yyyzz_0[i] = 2.0 * g_0_z_0_yyyzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_yyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zz_0_yyyz_1[i] * fi_abcd_0 +
                               g_0_zz_0_yyyzz_0[i] * pb_z + g_0_zz_0_yyyzz_1[i] * wp_z[i];

        g_0_zzz_0_yyzzz_0[i] = 2.0 * g_0_z_0_yyzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_yyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zz_0_yyzz_1[i] * fi_abcd_0 +
                               g_0_zz_0_yyzzz_0[i] * pb_z + g_0_zz_0_yyzzz_1[i] * wp_z[i];

        g_0_zzz_0_yzzzz_0[i] = 2.0 * g_0_z_0_yzzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_yzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zz_0_yzzz_1[i] * fi_abcd_0 +
                               g_0_zz_0_yzzzz_0[i] * pb_z + g_0_zz_0_yzzzz_1[i] * wp_z[i];

        g_0_zzz_0_zzzzz_0[i] = 2.0 * g_0_z_0_zzzzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_zzzzz_1[i] * fti_ab_0 + 5.0 * g_0_zz_0_zzzz_1[i] * fi_abcd_0 +
                               g_0_zz_0_zzzzz_0[i] * pb_z + g_0_zz_0_zzzzz_1[i] * wp_z[i];
    }
}

}  // namespace erirec

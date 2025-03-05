#include "ThreeCenterElectronRepulsionPrimRecFSH.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_fsh(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_fsh,
                                 size_t idx_eri_0_psh,
                                 size_t idx_eri_1_psh,
                                 size_t idx_eri_1_dsg,
                                 size_t idx_eri_1_dsh,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto c_exps = factors.data(0);

    auto d_exps = factors.data(1);

    // Set up R(WA) distances

    auto wa_x = factors.data(idx_wa);

    auto wa_y = factors.data(idx_wa + 1);

    auto wa_z = factors.data(idx_wa + 2);

    /// Set up components of auxilary buffer : PSH

    auto g_x_0_xxxxx_0 = pbuffer.data(idx_eri_0_psh);

    auto g_x_0_xxxxy_0 = pbuffer.data(idx_eri_0_psh + 1);

    auto g_x_0_xxxxz_0 = pbuffer.data(idx_eri_0_psh + 2);

    auto g_x_0_xxxyy_0 = pbuffer.data(idx_eri_0_psh + 3);

    auto g_x_0_xxxyz_0 = pbuffer.data(idx_eri_0_psh + 4);

    auto g_x_0_xxxzz_0 = pbuffer.data(idx_eri_0_psh + 5);

    auto g_x_0_xxyyy_0 = pbuffer.data(idx_eri_0_psh + 6);

    auto g_x_0_xxyyz_0 = pbuffer.data(idx_eri_0_psh + 7);

    auto g_x_0_xxyzz_0 = pbuffer.data(idx_eri_0_psh + 8);

    auto g_x_0_xxzzz_0 = pbuffer.data(idx_eri_0_psh + 9);

    auto g_x_0_xyyyy_0 = pbuffer.data(idx_eri_0_psh + 10);

    auto g_x_0_xyyyz_0 = pbuffer.data(idx_eri_0_psh + 11);

    auto g_x_0_xyyzz_0 = pbuffer.data(idx_eri_0_psh + 12);

    auto g_x_0_xyzzz_0 = pbuffer.data(idx_eri_0_psh + 13);

    auto g_x_0_xzzzz_0 = pbuffer.data(idx_eri_0_psh + 14);

    auto g_x_0_yyyyy_0 = pbuffer.data(idx_eri_0_psh + 15);

    auto g_x_0_yyyyz_0 = pbuffer.data(idx_eri_0_psh + 16);

    auto g_x_0_yyyzz_0 = pbuffer.data(idx_eri_0_psh + 17);

    auto g_x_0_yyzzz_0 = pbuffer.data(idx_eri_0_psh + 18);

    auto g_x_0_yzzzz_0 = pbuffer.data(idx_eri_0_psh + 19);

    auto g_x_0_zzzzz_0 = pbuffer.data(idx_eri_0_psh + 20);

    auto g_y_0_xxxxx_0 = pbuffer.data(idx_eri_0_psh + 21);

    auto g_y_0_xxxxy_0 = pbuffer.data(idx_eri_0_psh + 22);

    auto g_y_0_xxxxz_0 = pbuffer.data(idx_eri_0_psh + 23);

    auto g_y_0_xxxyy_0 = pbuffer.data(idx_eri_0_psh + 24);

    auto g_y_0_xxxyz_0 = pbuffer.data(idx_eri_0_psh + 25);

    auto g_y_0_xxxzz_0 = pbuffer.data(idx_eri_0_psh + 26);

    auto g_y_0_xxyyy_0 = pbuffer.data(idx_eri_0_psh + 27);

    auto g_y_0_xxyyz_0 = pbuffer.data(idx_eri_0_psh + 28);

    auto g_y_0_xxyzz_0 = pbuffer.data(idx_eri_0_psh + 29);

    auto g_y_0_xxzzz_0 = pbuffer.data(idx_eri_0_psh + 30);

    auto g_y_0_xyyyy_0 = pbuffer.data(idx_eri_0_psh + 31);

    auto g_y_0_xyyyz_0 = pbuffer.data(idx_eri_0_psh + 32);

    auto g_y_0_xyyzz_0 = pbuffer.data(idx_eri_0_psh + 33);

    auto g_y_0_xyzzz_0 = pbuffer.data(idx_eri_0_psh + 34);

    auto g_y_0_xzzzz_0 = pbuffer.data(idx_eri_0_psh + 35);

    auto g_y_0_yyyyy_0 = pbuffer.data(idx_eri_0_psh + 36);

    auto g_y_0_yyyyz_0 = pbuffer.data(idx_eri_0_psh + 37);

    auto g_y_0_yyyzz_0 = pbuffer.data(idx_eri_0_psh + 38);

    auto g_y_0_yyzzz_0 = pbuffer.data(idx_eri_0_psh + 39);

    auto g_y_0_yzzzz_0 = pbuffer.data(idx_eri_0_psh + 40);

    auto g_y_0_zzzzz_0 = pbuffer.data(idx_eri_0_psh + 41);

    auto g_z_0_xxxxx_0 = pbuffer.data(idx_eri_0_psh + 42);

    auto g_z_0_xxxxy_0 = pbuffer.data(idx_eri_0_psh + 43);

    auto g_z_0_xxxxz_0 = pbuffer.data(idx_eri_0_psh + 44);

    auto g_z_0_xxxyy_0 = pbuffer.data(idx_eri_0_psh + 45);

    auto g_z_0_xxxyz_0 = pbuffer.data(idx_eri_0_psh + 46);

    auto g_z_0_xxxzz_0 = pbuffer.data(idx_eri_0_psh + 47);

    auto g_z_0_xxyyy_0 = pbuffer.data(idx_eri_0_psh + 48);

    auto g_z_0_xxyyz_0 = pbuffer.data(idx_eri_0_psh + 49);

    auto g_z_0_xxyzz_0 = pbuffer.data(idx_eri_0_psh + 50);

    auto g_z_0_xxzzz_0 = pbuffer.data(idx_eri_0_psh + 51);

    auto g_z_0_xyyyy_0 = pbuffer.data(idx_eri_0_psh + 52);

    auto g_z_0_xyyyz_0 = pbuffer.data(idx_eri_0_psh + 53);

    auto g_z_0_xyyzz_0 = pbuffer.data(idx_eri_0_psh + 54);

    auto g_z_0_xyzzz_0 = pbuffer.data(idx_eri_0_psh + 55);

    auto g_z_0_xzzzz_0 = pbuffer.data(idx_eri_0_psh + 56);

    auto g_z_0_yyyyy_0 = pbuffer.data(idx_eri_0_psh + 57);

    auto g_z_0_yyyyz_0 = pbuffer.data(idx_eri_0_psh + 58);

    auto g_z_0_yyyzz_0 = pbuffer.data(idx_eri_0_psh + 59);

    auto g_z_0_yyzzz_0 = pbuffer.data(idx_eri_0_psh + 60);

    auto g_z_0_yzzzz_0 = pbuffer.data(idx_eri_0_psh + 61);

    auto g_z_0_zzzzz_0 = pbuffer.data(idx_eri_0_psh + 62);

    /// Set up components of auxilary buffer : PSH

    auto g_x_0_xxxxx_1 = pbuffer.data(idx_eri_1_psh);

    auto g_x_0_xxxxy_1 = pbuffer.data(idx_eri_1_psh + 1);

    auto g_x_0_xxxxz_1 = pbuffer.data(idx_eri_1_psh + 2);

    auto g_x_0_xxxyy_1 = pbuffer.data(idx_eri_1_psh + 3);

    auto g_x_0_xxxyz_1 = pbuffer.data(idx_eri_1_psh + 4);

    auto g_x_0_xxxzz_1 = pbuffer.data(idx_eri_1_psh + 5);

    auto g_x_0_xxyyy_1 = pbuffer.data(idx_eri_1_psh + 6);

    auto g_x_0_xxyyz_1 = pbuffer.data(idx_eri_1_psh + 7);

    auto g_x_0_xxyzz_1 = pbuffer.data(idx_eri_1_psh + 8);

    auto g_x_0_xxzzz_1 = pbuffer.data(idx_eri_1_psh + 9);

    auto g_x_0_xyyyy_1 = pbuffer.data(idx_eri_1_psh + 10);

    auto g_x_0_xyyyz_1 = pbuffer.data(idx_eri_1_psh + 11);

    auto g_x_0_xyyzz_1 = pbuffer.data(idx_eri_1_psh + 12);

    auto g_x_0_xyzzz_1 = pbuffer.data(idx_eri_1_psh + 13);

    auto g_x_0_xzzzz_1 = pbuffer.data(idx_eri_1_psh + 14);

    auto g_x_0_yyyyy_1 = pbuffer.data(idx_eri_1_psh + 15);

    auto g_x_0_yyyyz_1 = pbuffer.data(idx_eri_1_psh + 16);

    auto g_x_0_yyyzz_1 = pbuffer.data(idx_eri_1_psh + 17);

    auto g_x_0_yyzzz_1 = pbuffer.data(idx_eri_1_psh + 18);

    auto g_x_0_yzzzz_1 = pbuffer.data(idx_eri_1_psh + 19);

    auto g_x_0_zzzzz_1 = pbuffer.data(idx_eri_1_psh + 20);

    auto g_y_0_xxxxx_1 = pbuffer.data(idx_eri_1_psh + 21);

    auto g_y_0_xxxxy_1 = pbuffer.data(idx_eri_1_psh + 22);

    auto g_y_0_xxxxz_1 = pbuffer.data(idx_eri_1_psh + 23);

    auto g_y_0_xxxyy_1 = pbuffer.data(idx_eri_1_psh + 24);

    auto g_y_0_xxxyz_1 = pbuffer.data(idx_eri_1_psh + 25);

    auto g_y_0_xxxzz_1 = pbuffer.data(idx_eri_1_psh + 26);

    auto g_y_0_xxyyy_1 = pbuffer.data(idx_eri_1_psh + 27);

    auto g_y_0_xxyyz_1 = pbuffer.data(idx_eri_1_psh + 28);

    auto g_y_0_xxyzz_1 = pbuffer.data(idx_eri_1_psh + 29);

    auto g_y_0_xxzzz_1 = pbuffer.data(idx_eri_1_psh + 30);

    auto g_y_0_xyyyy_1 = pbuffer.data(idx_eri_1_psh + 31);

    auto g_y_0_xyyyz_1 = pbuffer.data(idx_eri_1_psh + 32);

    auto g_y_0_xyyzz_1 = pbuffer.data(idx_eri_1_psh + 33);

    auto g_y_0_xyzzz_1 = pbuffer.data(idx_eri_1_psh + 34);

    auto g_y_0_xzzzz_1 = pbuffer.data(idx_eri_1_psh + 35);

    auto g_y_0_yyyyy_1 = pbuffer.data(idx_eri_1_psh + 36);

    auto g_y_0_yyyyz_1 = pbuffer.data(idx_eri_1_psh + 37);

    auto g_y_0_yyyzz_1 = pbuffer.data(idx_eri_1_psh + 38);

    auto g_y_0_yyzzz_1 = pbuffer.data(idx_eri_1_psh + 39);

    auto g_y_0_yzzzz_1 = pbuffer.data(idx_eri_1_psh + 40);

    auto g_y_0_zzzzz_1 = pbuffer.data(idx_eri_1_psh + 41);

    auto g_z_0_xxxxx_1 = pbuffer.data(idx_eri_1_psh + 42);

    auto g_z_0_xxxxy_1 = pbuffer.data(idx_eri_1_psh + 43);

    auto g_z_0_xxxxz_1 = pbuffer.data(idx_eri_1_psh + 44);

    auto g_z_0_xxxyy_1 = pbuffer.data(idx_eri_1_psh + 45);

    auto g_z_0_xxxyz_1 = pbuffer.data(idx_eri_1_psh + 46);

    auto g_z_0_xxxzz_1 = pbuffer.data(idx_eri_1_psh + 47);

    auto g_z_0_xxyyy_1 = pbuffer.data(idx_eri_1_psh + 48);

    auto g_z_0_xxyyz_1 = pbuffer.data(idx_eri_1_psh + 49);

    auto g_z_0_xxyzz_1 = pbuffer.data(idx_eri_1_psh + 50);

    auto g_z_0_xxzzz_1 = pbuffer.data(idx_eri_1_psh + 51);

    auto g_z_0_xyyyy_1 = pbuffer.data(idx_eri_1_psh + 52);

    auto g_z_0_xyyyz_1 = pbuffer.data(idx_eri_1_psh + 53);

    auto g_z_0_xyyzz_1 = pbuffer.data(idx_eri_1_psh + 54);

    auto g_z_0_xyzzz_1 = pbuffer.data(idx_eri_1_psh + 55);

    auto g_z_0_xzzzz_1 = pbuffer.data(idx_eri_1_psh + 56);

    auto g_z_0_yyyyy_1 = pbuffer.data(idx_eri_1_psh + 57);

    auto g_z_0_yyyyz_1 = pbuffer.data(idx_eri_1_psh + 58);

    auto g_z_0_yyyzz_1 = pbuffer.data(idx_eri_1_psh + 59);

    auto g_z_0_yyzzz_1 = pbuffer.data(idx_eri_1_psh + 60);

    auto g_z_0_yzzzz_1 = pbuffer.data(idx_eri_1_psh + 61);

    auto g_z_0_zzzzz_1 = pbuffer.data(idx_eri_1_psh + 62);

    /// Set up components of auxilary buffer : DSG

    auto g_xx_0_xxxx_1 = pbuffer.data(idx_eri_1_dsg);

    auto g_xx_0_xxxy_1 = pbuffer.data(idx_eri_1_dsg + 1);

    auto g_xx_0_xxxz_1 = pbuffer.data(idx_eri_1_dsg + 2);

    auto g_xx_0_xxyy_1 = pbuffer.data(idx_eri_1_dsg + 3);

    auto g_xx_0_xxyz_1 = pbuffer.data(idx_eri_1_dsg + 4);

    auto g_xx_0_xxzz_1 = pbuffer.data(idx_eri_1_dsg + 5);

    auto g_xx_0_xyyy_1 = pbuffer.data(idx_eri_1_dsg + 6);

    auto g_xx_0_xyyz_1 = pbuffer.data(idx_eri_1_dsg + 7);

    auto g_xx_0_xyzz_1 = pbuffer.data(idx_eri_1_dsg + 8);

    auto g_xx_0_xzzz_1 = pbuffer.data(idx_eri_1_dsg + 9);

    auto g_xx_0_yyyy_1 = pbuffer.data(idx_eri_1_dsg + 10);

    auto g_xx_0_yyyz_1 = pbuffer.data(idx_eri_1_dsg + 11);

    auto g_xx_0_yyzz_1 = pbuffer.data(idx_eri_1_dsg + 12);

    auto g_xx_0_yzzz_1 = pbuffer.data(idx_eri_1_dsg + 13);

    auto g_xx_0_zzzz_1 = pbuffer.data(idx_eri_1_dsg + 14);

    auto g_yy_0_xxxx_1 = pbuffer.data(idx_eri_1_dsg + 45);

    auto g_yy_0_xxxy_1 = pbuffer.data(idx_eri_1_dsg + 46);

    auto g_yy_0_xxxz_1 = pbuffer.data(idx_eri_1_dsg + 47);

    auto g_yy_0_xxyy_1 = pbuffer.data(idx_eri_1_dsg + 48);

    auto g_yy_0_xxyz_1 = pbuffer.data(idx_eri_1_dsg + 49);

    auto g_yy_0_xxzz_1 = pbuffer.data(idx_eri_1_dsg + 50);

    auto g_yy_0_xyyy_1 = pbuffer.data(idx_eri_1_dsg + 51);

    auto g_yy_0_xyyz_1 = pbuffer.data(idx_eri_1_dsg + 52);

    auto g_yy_0_xyzz_1 = pbuffer.data(idx_eri_1_dsg + 53);

    auto g_yy_0_xzzz_1 = pbuffer.data(idx_eri_1_dsg + 54);

    auto g_yy_0_yyyy_1 = pbuffer.data(idx_eri_1_dsg + 55);

    auto g_yy_0_yyyz_1 = pbuffer.data(idx_eri_1_dsg + 56);

    auto g_yy_0_yyzz_1 = pbuffer.data(idx_eri_1_dsg + 57);

    auto g_yy_0_yzzz_1 = pbuffer.data(idx_eri_1_dsg + 58);

    auto g_yy_0_zzzz_1 = pbuffer.data(idx_eri_1_dsg + 59);

    auto g_yz_0_xxyz_1 = pbuffer.data(idx_eri_1_dsg + 64);

    auto g_yz_0_xyyz_1 = pbuffer.data(idx_eri_1_dsg + 67);

    auto g_yz_0_xyzz_1 = pbuffer.data(idx_eri_1_dsg + 68);

    auto g_yz_0_yyyz_1 = pbuffer.data(idx_eri_1_dsg + 71);

    auto g_yz_0_yyzz_1 = pbuffer.data(idx_eri_1_dsg + 72);

    auto g_yz_0_yzzz_1 = pbuffer.data(idx_eri_1_dsg + 73);

    auto g_zz_0_xxxx_1 = pbuffer.data(idx_eri_1_dsg + 75);

    auto g_zz_0_xxxy_1 = pbuffer.data(idx_eri_1_dsg + 76);

    auto g_zz_0_xxxz_1 = pbuffer.data(idx_eri_1_dsg + 77);

    auto g_zz_0_xxyy_1 = pbuffer.data(idx_eri_1_dsg + 78);

    auto g_zz_0_xxyz_1 = pbuffer.data(idx_eri_1_dsg + 79);

    auto g_zz_0_xxzz_1 = pbuffer.data(idx_eri_1_dsg + 80);

    auto g_zz_0_xyyy_1 = pbuffer.data(idx_eri_1_dsg + 81);

    auto g_zz_0_xyyz_1 = pbuffer.data(idx_eri_1_dsg + 82);

    auto g_zz_0_xyzz_1 = pbuffer.data(idx_eri_1_dsg + 83);

    auto g_zz_0_xzzz_1 = pbuffer.data(idx_eri_1_dsg + 84);

    auto g_zz_0_yyyy_1 = pbuffer.data(idx_eri_1_dsg + 85);

    auto g_zz_0_yyyz_1 = pbuffer.data(idx_eri_1_dsg + 86);

    auto g_zz_0_yyzz_1 = pbuffer.data(idx_eri_1_dsg + 87);

    auto g_zz_0_yzzz_1 = pbuffer.data(idx_eri_1_dsg + 88);

    auto g_zz_0_zzzz_1 = pbuffer.data(idx_eri_1_dsg + 89);

    /// Set up components of auxilary buffer : DSH

    auto g_xx_0_xxxxx_1 = pbuffer.data(idx_eri_1_dsh);

    auto g_xx_0_xxxxy_1 = pbuffer.data(idx_eri_1_dsh + 1);

    auto g_xx_0_xxxxz_1 = pbuffer.data(idx_eri_1_dsh + 2);

    auto g_xx_0_xxxyy_1 = pbuffer.data(idx_eri_1_dsh + 3);

    auto g_xx_0_xxxyz_1 = pbuffer.data(idx_eri_1_dsh + 4);

    auto g_xx_0_xxxzz_1 = pbuffer.data(idx_eri_1_dsh + 5);

    auto g_xx_0_xxyyy_1 = pbuffer.data(idx_eri_1_dsh + 6);

    auto g_xx_0_xxyyz_1 = pbuffer.data(idx_eri_1_dsh + 7);

    auto g_xx_0_xxyzz_1 = pbuffer.data(idx_eri_1_dsh + 8);

    auto g_xx_0_xxzzz_1 = pbuffer.data(idx_eri_1_dsh + 9);

    auto g_xx_0_xyyyy_1 = pbuffer.data(idx_eri_1_dsh + 10);

    auto g_xx_0_xyyyz_1 = pbuffer.data(idx_eri_1_dsh + 11);

    auto g_xx_0_xyyzz_1 = pbuffer.data(idx_eri_1_dsh + 12);

    auto g_xx_0_xyzzz_1 = pbuffer.data(idx_eri_1_dsh + 13);

    auto g_xx_0_xzzzz_1 = pbuffer.data(idx_eri_1_dsh + 14);

    auto g_xx_0_yyyyy_1 = pbuffer.data(idx_eri_1_dsh + 15);

    auto g_xx_0_yyyyz_1 = pbuffer.data(idx_eri_1_dsh + 16);

    auto g_xx_0_yyyzz_1 = pbuffer.data(idx_eri_1_dsh + 17);

    auto g_xx_0_yyzzz_1 = pbuffer.data(idx_eri_1_dsh + 18);

    auto g_xx_0_yzzzz_1 = pbuffer.data(idx_eri_1_dsh + 19);

    auto g_xx_0_zzzzz_1 = pbuffer.data(idx_eri_1_dsh + 20);

    auto g_xy_0_xxxxy_1 = pbuffer.data(idx_eri_1_dsh + 22);

    auto g_xy_0_xxxyy_1 = pbuffer.data(idx_eri_1_dsh + 24);

    auto g_xy_0_xxyyy_1 = pbuffer.data(idx_eri_1_dsh + 27);

    auto g_xy_0_xyyyy_1 = pbuffer.data(idx_eri_1_dsh + 31);

    auto g_xz_0_xxxxx_1 = pbuffer.data(idx_eri_1_dsh + 42);

    auto g_xz_0_xxxxz_1 = pbuffer.data(idx_eri_1_dsh + 44);

    auto g_xz_0_xxxzz_1 = pbuffer.data(idx_eri_1_dsh + 47);

    auto g_xz_0_xxzzz_1 = pbuffer.data(idx_eri_1_dsh + 51);

    auto g_xz_0_xzzzz_1 = pbuffer.data(idx_eri_1_dsh + 56);

    auto g_yy_0_xxxxx_1 = pbuffer.data(idx_eri_1_dsh + 63);

    auto g_yy_0_xxxxy_1 = pbuffer.data(idx_eri_1_dsh + 64);

    auto g_yy_0_xxxxz_1 = pbuffer.data(idx_eri_1_dsh + 65);

    auto g_yy_0_xxxyy_1 = pbuffer.data(idx_eri_1_dsh + 66);

    auto g_yy_0_xxxyz_1 = pbuffer.data(idx_eri_1_dsh + 67);

    auto g_yy_0_xxxzz_1 = pbuffer.data(idx_eri_1_dsh + 68);

    auto g_yy_0_xxyyy_1 = pbuffer.data(idx_eri_1_dsh + 69);

    auto g_yy_0_xxyyz_1 = pbuffer.data(idx_eri_1_dsh + 70);

    auto g_yy_0_xxyzz_1 = pbuffer.data(idx_eri_1_dsh + 71);

    auto g_yy_0_xxzzz_1 = pbuffer.data(idx_eri_1_dsh + 72);

    auto g_yy_0_xyyyy_1 = pbuffer.data(idx_eri_1_dsh + 73);

    auto g_yy_0_xyyyz_1 = pbuffer.data(idx_eri_1_dsh + 74);

    auto g_yy_0_xyyzz_1 = pbuffer.data(idx_eri_1_dsh + 75);

    auto g_yy_0_xyzzz_1 = pbuffer.data(idx_eri_1_dsh + 76);

    auto g_yy_0_xzzzz_1 = pbuffer.data(idx_eri_1_dsh + 77);

    auto g_yy_0_yyyyy_1 = pbuffer.data(idx_eri_1_dsh + 78);

    auto g_yy_0_yyyyz_1 = pbuffer.data(idx_eri_1_dsh + 79);

    auto g_yy_0_yyyzz_1 = pbuffer.data(idx_eri_1_dsh + 80);

    auto g_yy_0_yyzzz_1 = pbuffer.data(idx_eri_1_dsh + 81);

    auto g_yy_0_yzzzz_1 = pbuffer.data(idx_eri_1_dsh + 82);

    auto g_yy_0_zzzzz_1 = pbuffer.data(idx_eri_1_dsh + 83);

    auto g_yz_0_xxxyz_1 = pbuffer.data(idx_eri_1_dsh + 88);

    auto g_yz_0_xxyyz_1 = pbuffer.data(idx_eri_1_dsh + 91);

    auto g_yz_0_xxyzz_1 = pbuffer.data(idx_eri_1_dsh + 92);

    auto g_yz_0_xyyyz_1 = pbuffer.data(idx_eri_1_dsh + 95);

    auto g_yz_0_xyyzz_1 = pbuffer.data(idx_eri_1_dsh + 96);

    auto g_yz_0_xyzzz_1 = pbuffer.data(idx_eri_1_dsh + 97);

    auto g_yz_0_yyyyy_1 = pbuffer.data(idx_eri_1_dsh + 99);

    auto g_yz_0_yyyyz_1 = pbuffer.data(idx_eri_1_dsh + 100);

    auto g_yz_0_yyyzz_1 = pbuffer.data(idx_eri_1_dsh + 101);

    auto g_yz_0_yyzzz_1 = pbuffer.data(idx_eri_1_dsh + 102);

    auto g_yz_0_yzzzz_1 = pbuffer.data(idx_eri_1_dsh + 103);

    auto g_yz_0_zzzzz_1 = pbuffer.data(idx_eri_1_dsh + 104);

    auto g_zz_0_xxxxx_1 = pbuffer.data(idx_eri_1_dsh + 105);

    auto g_zz_0_xxxxy_1 = pbuffer.data(idx_eri_1_dsh + 106);

    auto g_zz_0_xxxxz_1 = pbuffer.data(idx_eri_1_dsh + 107);

    auto g_zz_0_xxxyy_1 = pbuffer.data(idx_eri_1_dsh + 108);

    auto g_zz_0_xxxyz_1 = pbuffer.data(idx_eri_1_dsh + 109);

    auto g_zz_0_xxxzz_1 = pbuffer.data(idx_eri_1_dsh + 110);

    auto g_zz_0_xxyyy_1 = pbuffer.data(idx_eri_1_dsh + 111);

    auto g_zz_0_xxyyz_1 = pbuffer.data(idx_eri_1_dsh + 112);

    auto g_zz_0_xxyzz_1 = pbuffer.data(idx_eri_1_dsh + 113);

    auto g_zz_0_xxzzz_1 = pbuffer.data(idx_eri_1_dsh + 114);

    auto g_zz_0_xyyyy_1 = pbuffer.data(idx_eri_1_dsh + 115);

    auto g_zz_0_xyyyz_1 = pbuffer.data(idx_eri_1_dsh + 116);

    auto g_zz_0_xyyzz_1 = pbuffer.data(idx_eri_1_dsh + 117);

    auto g_zz_0_xyzzz_1 = pbuffer.data(idx_eri_1_dsh + 118);

    auto g_zz_0_xzzzz_1 = pbuffer.data(idx_eri_1_dsh + 119);

    auto g_zz_0_yyyyy_1 = pbuffer.data(idx_eri_1_dsh + 120);

    auto g_zz_0_yyyyz_1 = pbuffer.data(idx_eri_1_dsh + 121);

    auto g_zz_0_yyyzz_1 = pbuffer.data(idx_eri_1_dsh + 122);

    auto g_zz_0_yyzzz_1 = pbuffer.data(idx_eri_1_dsh + 123);

    auto g_zz_0_yzzzz_1 = pbuffer.data(idx_eri_1_dsh + 124);

    auto g_zz_0_zzzzz_1 = pbuffer.data(idx_eri_1_dsh + 125);

    /// Set up 0-21 components of targeted buffer : FSH

    auto g_xxx_0_xxxxx_0 = pbuffer.data(idx_eri_0_fsh);

    auto g_xxx_0_xxxxy_0 = pbuffer.data(idx_eri_0_fsh + 1);

    auto g_xxx_0_xxxxz_0 = pbuffer.data(idx_eri_0_fsh + 2);

    auto g_xxx_0_xxxyy_0 = pbuffer.data(idx_eri_0_fsh + 3);

    auto g_xxx_0_xxxyz_0 = pbuffer.data(idx_eri_0_fsh + 4);

    auto g_xxx_0_xxxzz_0 = pbuffer.data(idx_eri_0_fsh + 5);

    auto g_xxx_0_xxyyy_0 = pbuffer.data(idx_eri_0_fsh + 6);

    auto g_xxx_0_xxyyz_0 = pbuffer.data(idx_eri_0_fsh + 7);

    auto g_xxx_0_xxyzz_0 = pbuffer.data(idx_eri_0_fsh + 8);

    auto g_xxx_0_xxzzz_0 = pbuffer.data(idx_eri_0_fsh + 9);

    auto g_xxx_0_xyyyy_0 = pbuffer.data(idx_eri_0_fsh + 10);

    auto g_xxx_0_xyyyz_0 = pbuffer.data(idx_eri_0_fsh + 11);

    auto g_xxx_0_xyyzz_0 = pbuffer.data(idx_eri_0_fsh + 12);

    auto g_xxx_0_xyzzz_0 = pbuffer.data(idx_eri_0_fsh + 13);

    auto g_xxx_0_xzzzz_0 = pbuffer.data(idx_eri_0_fsh + 14);

    auto g_xxx_0_yyyyy_0 = pbuffer.data(idx_eri_0_fsh + 15);

    auto g_xxx_0_yyyyz_0 = pbuffer.data(idx_eri_0_fsh + 16);

    auto g_xxx_0_yyyzz_0 = pbuffer.data(idx_eri_0_fsh + 17);

    auto g_xxx_0_yyzzz_0 = pbuffer.data(idx_eri_0_fsh + 18);

    auto g_xxx_0_yzzzz_0 = pbuffer.data(idx_eri_0_fsh + 19);

    auto g_xxx_0_zzzzz_0 = pbuffer.data(idx_eri_0_fsh + 20);

    #pragma omp simd aligned(g_x_0_xxxxx_0, g_x_0_xxxxx_1, g_x_0_xxxxy_0, g_x_0_xxxxy_1, g_x_0_xxxxz_0, g_x_0_xxxxz_1, g_x_0_xxxyy_0, g_x_0_xxxyy_1, g_x_0_xxxyz_0, g_x_0_xxxyz_1, g_x_0_xxxzz_0, g_x_0_xxxzz_1, g_x_0_xxyyy_0, g_x_0_xxyyy_1, g_x_0_xxyyz_0, g_x_0_xxyyz_1, g_x_0_xxyzz_0, g_x_0_xxyzz_1, g_x_0_xxzzz_0, g_x_0_xxzzz_1, g_x_0_xyyyy_0, g_x_0_xyyyy_1, g_x_0_xyyyz_0, g_x_0_xyyyz_1, g_x_0_xyyzz_0, g_x_0_xyyzz_1, g_x_0_xyzzz_0, g_x_0_xyzzz_1, g_x_0_xzzzz_0, g_x_0_xzzzz_1, g_x_0_yyyyy_0, g_x_0_yyyyy_1, g_x_0_yyyyz_0, g_x_0_yyyyz_1, g_x_0_yyyzz_0, g_x_0_yyyzz_1, g_x_0_yyzzz_0, g_x_0_yyzzz_1, g_x_0_yzzzz_0, g_x_0_yzzzz_1, g_x_0_zzzzz_0, g_x_0_zzzzz_1, g_xx_0_xxxx_1, g_xx_0_xxxxx_1, g_xx_0_xxxxy_1, g_xx_0_xxxxz_1, g_xx_0_xxxy_1, g_xx_0_xxxyy_1, g_xx_0_xxxyz_1, g_xx_0_xxxz_1, g_xx_0_xxxzz_1, g_xx_0_xxyy_1, g_xx_0_xxyyy_1, g_xx_0_xxyyz_1, g_xx_0_xxyz_1, g_xx_0_xxyzz_1, g_xx_0_xxzz_1, g_xx_0_xxzzz_1, g_xx_0_xyyy_1, g_xx_0_xyyyy_1, g_xx_0_xyyyz_1, g_xx_0_xyyz_1, g_xx_0_xyyzz_1, g_xx_0_xyzz_1, g_xx_0_xyzzz_1, g_xx_0_xzzz_1, g_xx_0_xzzzz_1, g_xx_0_yyyy_1, g_xx_0_yyyyy_1, g_xx_0_yyyyz_1, g_xx_0_yyyz_1, g_xx_0_yyyzz_1, g_xx_0_yyzz_1, g_xx_0_yyzzz_1, g_xx_0_yzzz_1, g_xx_0_yzzzz_1, g_xx_0_zzzz_1, g_xx_0_zzzzz_1, g_xxx_0_xxxxx_0, g_xxx_0_xxxxy_0, g_xxx_0_xxxxz_0, g_xxx_0_xxxyy_0, g_xxx_0_xxxyz_0, g_xxx_0_xxxzz_0, g_xxx_0_xxyyy_0, g_xxx_0_xxyyz_0, g_xxx_0_xxyzz_0, g_xxx_0_xxzzz_0, g_xxx_0_xyyyy_0, g_xxx_0_xyyyz_0, g_xxx_0_xyyzz_0, g_xxx_0_xyzzz_0, g_xxx_0_xzzzz_0, g_xxx_0_yyyyy_0, g_xxx_0_yyyyz_0, g_xxx_0_yyyzz_0, g_xxx_0_yyzzz_0, g_xxx_0_yzzzz_0, g_xxx_0_zzzzz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxx_0_xxxxx_0[i] = 2.0 * g_x_0_xxxxx_0[i] * fbe_0 - 2.0 * g_x_0_xxxxx_1[i] * fz_be_0 + 5.0 * g_xx_0_xxxx_1[i] * fi_acd_0 + g_xx_0_xxxxx_1[i] * wa_x[i];

        g_xxx_0_xxxxy_0[i] = 2.0 * g_x_0_xxxxy_0[i] * fbe_0 - 2.0 * g_x_0_xxxxy_1[i] * fz_be_0 + 4.0 * g_xx_0_xxxy_1[i] * fi_acd_0 + g_xx_0_xxxxy_1[i] * wa_x[i];

        g_xxx_0_xxxxz_0[i] = 2.0 * g_x_0_xxxxz_0[i] * fbe_0 - 2.0 * g_x_0_xxxxz_1[i] * fz_be_0 + 4.0 * g_xx_0_xxxz_1[i] * fi_acd_0 + g_xx_0_xxxxz_1[i] * wa_x[i];

        g_xxx_0_xxxyy_0[i] = 2.0 * g_x_0_xxxyy_0[i] * fbe_0 - 2.0 * g_x_0_xxxyy_1[i] * fz_be_0 + 3.0 * g_xx_0_xxyy_1[i] * fi_acd_0 + g_xx_0_xxxyy_1[i] * wa_x[i];

        g_xxx_0_xxxyz_0[i] = 2.0 * g_x_0_xxxyz_0[i] * fbe_0 - 2.0 * g_x_0_xxxyz_1[i] * fz_be_0 + 3.0 * g_xx_0_xxyz_1[i] * fi_acd_0 + g_xx_0_xxxyz_1[i] * wa_x[i];

        g_xxx_0_xxxzz_0[i] = 2.0 * g_x_0_xxxzz_0[i] * fbe_0 - 2.0 * g_x_0_xxxzz_1[i] * fz_be_0 + 3.0 * g_xx_0_xxzz_1[i] * fi_acd_0 + g_xx_0_xxxzz_1[i] * wa_x[i];

        g_xxx_0_xxyyy_0[i] = 2.0 * g_x_0_xxyyy_0[i] * fbe_0 - 2.0 * g_x_0_xxyyy_1[i] * fz_be_0 + 2.0 * g_xx_0_xyyy_1[i] * fi_acd_0 + g_xx_0_xxyyy_1[i] * wa_x[i];

        g_xxx_0_xxyyz_0[i] = 2.0 * g_x_0_xxyyz_0[i] * fbe_0 - 2.0 * g_x_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_xx_0_xyyz_1[i] * fi_acd_0 + g_xx_0_xxyyz_1[i] * wa_x[i];

        g_xxx_0_xxyzz_0[i] = 2.0 * g_x_0_xxyzz_0[i] * fbe_0 - 2.0 * g_x_0_xxyzz_1[i] * fz_be_0 + 2.0 * g_xx_0_xyzz_1[i] * fi_acd_0 + g_xx_0_xxyzz_1[i] * wa_x[i];

        g_xxx_0_xxzzz_0[i] = 2.0 * g_x_0_xxzzz_0[i] * fbe_0 - 2.0 * g_x_0_xxzzz_1[i] * fz_be_0 + 2.0 * g_xx_0_xzzz_1[i] * fi_acd_0 + g_xx_0_xxzzz_1[i] * wa_x[i];

        g_xxx_0_xyyyy_0[i] = 2.0 * g_x_0_xyyyy_0[i] * fbe_0 - 2.0 * g_x_0_xyyyy_1[i] * fz_be_0 + g_xx_0_yyyy_1[i] * fi_acd_0 + g_xx_0_xyyyy_1[i] * wa_x[i];

        g_xxx_0_xyyyz_0[i] = 2.0 * g_x_0_xyyyz_0[i] * fbe_0 - 2.0 * g_x_0_xyyyz_1[i] * fz_be_0 + g_xx_0_yyyz_1[i] * fi_acd_0 + g_xx_0_xyyyz_1[i] * wa_x[i];

        g_xxx_0_xyyzz_0[i] = 2.0 * g_x_0_xyyzz_0[i] * fbe_0 - 2.0 * g_x_0_xyyzz_1[i] * fz_be_0 + g_xx_0_yyzz_1[i] * fi_acd_0 + g_xx_0_xyyzz_1[i] * wa_x[i];

        g_xxx_0_xyzzz_0[i] = 2.0 * g_x_0_xyzzz_0[i] * fbe_0 - 2.0 * g_x_0_xyzzz_1[i] * fz_be_0 + g_xx_0_yzzz_1[i] * fi_acd_0 + g_xx_0_xyzzz_1[i] * wa_x[i];

        g_xxx_0_xzzzz_0[i] = 2.0 * g_x_0_xzzzz_0[i] * fbe_0 - 2.0 * g_x_0_xzzzz_1[i] * fz_be_0 + g_xx_0_zzzz_1[i] * fi_acd_0 + g_xx_0_xzzzz_1[i] * wa_x[i];

        g_xxx_0_yyyyy_0[i] = 2.0 * g_x_0_yyyyy_0[i] * fbe_0 - 2.0 * g_x_0_yyyyy_1[i] * fz_be_0 + g_xx_0_yyyyy_1[i] * wa_x[i];

        g_xxx_0_yyyyz_0[i] = 2.0 * g_x_0_yyyyz_0[i] * fbe_0 - 2.0 * g_x_0_yyyyz_1[i] * fz_be_0 + g_xx_0_yyyyz_1[i] * wa_x[i];

        g_xxx_0_yyyzz_0[i] = 2.0 * g_x_0_yyyzz_0[i] * fbe_0 - 2.0 * g_x_0_yyyzz_1[i] * fz_be_0 + g_xx_0_yyyzz_1[i] * wa_x[i];

        g_xxx_0_yyzzz_0[i] = 2.0 * g_x_0_yyzzz_0[i] * fbe_0 - 2.0 * g_x_0_yyzzz_1[i] * fz_be_0 + g_xx_0_yyzzz_1[i] * wa_x[i];

        g_xxx_0_yzzzz_0[i] = 2.0 * g_x_0_yzzzz_0[i] * fbe_0 - 2.0 * g_x_0_yzzzz_1[i] * fz_be_0 + g_xx_0_yzzzz_1[i] * wa_x[i];

        g_xxx_0_zzzzz_0[i] = 2.0 * g_x_0_zzzzz_0[i] * fbe_0 - 2.0 * g_x_0_zzzzz_1[i] * fz_be_0 + g_xx_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 21-42 components of targeted buffer : FSH

    auto g_xxy_0_xxxxx_0 = pbuffer.data(idx_eri_0_fsh + 21);

    auto g_xxy_0_xxxxy_0 = pbuffer.data(idx_eri_0_fsh + 22);

    auto g_xxy_0_xxxxz_0 = pbuffer.data(idx_eri_0_fsh + 23);

    auto g_xxy_0_xxxyy_0 = pbuffer.data(idx_eri_0_fsh + 24);

    auto g_xxy_0_xxxyz_0 = pbuffer.data(idx_eri_0_fsh + 25);

    auto g_xxy_0_xxxzz_0 = pbuffer.data(idx_eri_0_fsh + 26);

    auto g_xxy_0_xxyyy_0 = pbuffer.data(idx_eri_0_fsh + 27);

    auto g_xxy_0_xxyyz_0 = pbuffer.data(idx_eri_0_fsh + 28);

    auto g_xxy_0_xxyzz_0 = pbuffer.data(idx_eri_0_fsh + 29);

    auto g_xxy_0_xxzzz_0 = pbuffer.data(idx_eri_0_fsh + 30);

    auto g_xxy_0_xyyyy_0 = pbuffer.data(idx_eri_0_fsh + 31);

    auto g_xxy_0_xyyyz_0 = pbuffer.data(idx_eri_0_fsh + 32);

    auto g_xxy_0_xyyzz_0 = pbuffer.data(idx_eri_0_fsh + 33);

    auto g_xxy_0_xyzzz_0 = pbuffer.data(idx_eri_0_fsh + 34);

    auto g_xxy_0_xzzzz_0 = pbuffer.data(idx_eri_0_fsh + 35);

    auto g_xxy_0_yyyyy_0 = pbuffer.data(idx_eri_0_fsh + 36);

    auto g_xxy_0_yyyyz_0 = pbuffer.data(idx_eri_0_fsh + 37);

    auto g_xxy_0_yyyzz_0 = pbuffer.data(idx_eri_0_fsh + 38);

    auto g_xxy_0_yyzzz_0 = pbuffer.data(idx_eri_0_fsh + 39);

    auto g_xxy_0_yzzzz_0 = pbuffer.data(idx_eri_0_fsh + 40);

    auto g_xxy_0_zzzzz_0 = pbuffer.data(idx_eri_0_fsh + 41);

    #pragma omp simd aligned(g_xx_0_xxxx_1, g_xx_0_xxxxx_1, g_xx_0_xxxxy_1, g_xx_0_xxxxz_1, g_xx_0_xxxy_1, g_xx_0_xxxyy_1, g_xx_0_xxxyz_1, g_xx_0_xxxz_1, g_xx_0_xxxzz_1, g_xx_0_xxyy_1, g_xx_0_xxyyy_1, g_xx_0_xxyyz_1, g_xx_0_xxyz_1, g_xx_0_xxyzz_1, g_xx_0_xxzz_1, g_xx_0_xxzzz_1, g_xx_0_xyyy_1, g_xx_0_xyyyy_1, g_xx_0_xyyyz_1, g_xx_0_xyyz_1, g_xx_0_xyyzz_1, g_xx_0_xyzz_1, g_xx_0_xyzzz_1, g_xx_0_xzzz_1, g_xx_0_xzzzz_1, g_xx_0_yyyy_1, g_xx_0_yyyyy_1, g_xx_0_yyyyz_1, g_xx_0_yyyz_1, g_xx_0_yyyzz_1, g_xx_0_yyzz_1, g_xx_0_yyzzz_1, g_xx_0_yzzz_1, g_xx_0_yzzzz_1, g_xx_0_zzzz_1, g_xx_0_zzzzz_1, g_xxy_0_xxxxx_0, g_xxy_0_xxxxy_0, g_xxy_0_xxxxz_0, g_xxy_0_xxxyy_0, g_xxy_0_xxxyz_0, g_xxy_0_xxxzz_0, g_xxy_0_xxyyy_0, g_xxy_0_xxyyz_0, g_xxy_0_xxyzz_0, g_xxy_0_xxzzz_0, g_xxy_0_xyyyy_0, g_xxy_0_xyyyz_0, g_xxy_0_xyyzz_0, g_xxy_0_xyzzz_0, g_xxy_0_xzzzz_0, g_xxy_0_yyyyy_0, g_xxy_0_yyyyz_0, g_xxy_0_yyyzz_0, g_xxy_0_yyzzz_0, g_xxy_0_yzzzz_0, g_xxy_0_zzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxy_0_xxxxx_0[i] = g_xx_0_xxxxx_1[i] * wa_y[i];

        g_xxy_0_xxxxy_0[i] = g_xx_0_xxxx_1[i] * fi_acd_0 + g_xx_0_xxxxy_1[i] * wa_y[i];

        g_xxy_0_xxxxz_0[i] = g_xx_0_xxxxz_1[i] * wa_y[i];

        g_xxy_0_xxxyy_0[i] = 2.0 * g_xx_0_xxxy_1[i] * fi_acd_0 + g_xx_0_xxxyy_1[i] * wa_y[i];

        g_xxy_0_xxxyz_0[i] = g_xx_0_xxxz_1[i] * fi_acd_0 + g_xx_0_xxxyz_1[i] * wa_y[i];

        g_xxy_0_xxxzz_0[i] = g_xx_0_xxxzz_1[i] * wa_y[i];

        g_xxy_0_xxyyy_0[i] = 3.0 * g_xx_0_xxyy_1[i] * fi_acd_0 + g_xx_0_xxyyy_1[i] * wa_y[i];

        g_xxy_0_xxyyz_0[i] = 2.0 * g_xx_0_xxyz_1[i] * fi_acd_0 + g_xx_0_xxyyz_1[i] * wa_y[i];

        g_xxy_0_xxyzz_0[i] = g_xx_0_xxzz_1[i] * fi_acd_0 + g_xx_0_xxyzz_1[i] * wa_y[i];

        g_xxy_0_xxzzz_0[i] = g_xx_0_xxzzz_1[i] * wa_y[i];

        g_xxy_0_xyyyy_0[i] = 4.0 * g_xx_0_xyyy_1[i] * fi_acd_0 + g_xx_0_xyyyy_1[i] * wa_y[i];

        g_xxy_0_xyyyz_0[i] = 3.0 * g_xx_0_xyyz_1[i] * fi_acd_0 + g_xx_0_xyyyz_1[i] * wa_y[i];

        g_xxy_0_xyyzz_0[i] = 2.0 * g_xx_0_xyzz_1[i] * fi_acd_0 + g_xx_0_xyyzz_1[i] * wa_y[i];

        g_xxy_0_xyzzz_0[i] = g_xx_0_xzzz_1[i] * fi_acd_0 + g_xx_0_xyzzz_1[i] * wa_y[i];

        g_xxy_0_xzzzz_0[i] = g_xx_0_xzzzz_1[i] * wa_y[i];

        g_xxy_0_yyyyy_0[i] = 5.0 * g_xx_0_yyyy_1[i] * fi_acd_0 + g_xx_0_yyyyy_1[i] * wa_y[i];

        g_xxy_0_yyyyz_0[i] = 4.0 * g_xx_0_yyyz_1[i] * fi_acd_0 + g_xx_0_yyyyz_1[i] * wa_y[i];

        g_xxy_0_yyyzz_0[i] = 3.0 * g_xx_0_yyzz_1[i] * fi_acd_0 + g_xx_0_yyyzz_1[i] * wa_y[i];

        g_xxy_0_yyzzz_0[i] = 2.0 * g_xx_0_yzzz_1[i] * fi_acd_0 + g_xx_0_yyzzz_1[i] * wa_y[i];

        g_xxy_0_yzzzz_0[i] = g_xx_0_zzzz_1[i] * fi_acd_0 + g_xx_0_yzzzz_1[i] * wa_y[i];

        g_xxy_0_zzzzz_0[i] = g_xx_0_zzzzz_1[i] * wa_y[i];
    }

    /// Set up 42-63 components of targeted buffer : FSH

    auto g_xxz_0_xxxxx_0 = pbuffer.data(idx_eri_0_fsh + 42);

    auto g_xxz_0_xxxxy_0 = pbuffer.data(idx_eri_0_fsh + 43);

    auto g_xxz_0_xxxxz_0 = pbuffer.data(idx_eri_0_fsh + 44);

    auto g_xxz_0_xxxyy_0 = pbuffer.data(idx_eri_0_fsh + 45);

    auto g_xxz_0_xxxyz_0 = pbuffer.data(idx_eri_0_fsh + 46);

    auto g_xxz_0_xxxzz_0 = pbuffer.data(idx_eri_0_fsh + 47);

    auto g_xxz_0_xxyyy_0 = pbuffer.data(idx_eri_0_fsh + 48);

    auto g_xxz_0_xxyyz_0 = pbuffer.data(idx_eri_0_fsh + 49);

    auto g_xxz_0_xxyzz_0 = pbuffer.data(idx_eri_0_fsh + 50);

    auto g_xxz_0_xxzzz_0 = pbuffer.data(idx_eri_0_fsh + 51);

    auto g_xxz_0_xyyyy_0 = pbuffer.data(idx_eri_0_fsh + 52);

    auto g_xxz_0_xyyyz_0 = pbuffer.data(idx_eri_0_fsh + 53);

    auto g_xxz_0_xyyzz_0 = pbuffer.data(idx_eri_0_fsh + 54);

    auto g_xxz_0_xyzzz_0 = pbuffer.data(idx_eri_0_fsh + 55);

    auto g_xxz_0_xzzzz_0 = pbuffer.data(idx_eri_0_fsh + 56);

    auto g_xxz_0_yyyyy_0 = pbuffer.data(idx_eri_0_fsh + 57);

    auto g_xxz_0_yyyyz_0 = pbuffer.data(idx_eri_0_fsh + 58);

    auto g_xxz_0_yyyzz_0 = pbuffer.data(idx_eri_0_fsh + 59);

    auto g_xxz_0_yyzzz_0 = pbuffer.data(idx_eri_0_fsh + 60);

    auto g_xxz_0_yzzzz_0 = pbuffer.data(idx_eri_0_fsh + 61);

    auto g_xxz_0_zzzzz_0 = pbuffer.data(idx_eri_0_fsh + 62);

    #pragma omp simd aligned(g_xx_0_xxxx_1, g_xx_0_xxxxx_1, g_xx_0_xxxxy_1, g_xx_0_xxxxz_1, g_xx_0_xxxy_1, g_xx_0_xxxyy_1, g_xx_0_xxxyz_1, g_xx_0_xxxz_1, g_xx_0_xxxzz_1, g_xx_0_xxyy_1, g_xx_0_xxyyy_1, g_xx_0_xxyyz_1, g_xx_0_xxyz_1, g_xx_0_xxyzz_1, g_xx_0_xxzz_1, g_xx_0_xxzzz_1, g_xx_0_xyyy_1, g_xx_0_xyyyy_1, g_xx_0_xyyyz_1, g_xx_0_xyyz_1, g_xx_0_xyyzz_1, g_xx_0_xyzz_1, g_xx_0_xyzzz_1, g_xx_0_xzzz_1, g_xx_0_xzzzz_1, g_xx_0_yyyy_1, g_xx_0_yyyyy_1, g_xx_0_yyyyz_1, g_xx_0_yyyz_1, g_xx_0_yyyzz_1, g_xx_0_yyzz_1, g_xx_0_yyzzz_1, g_xx_0_yzzz_1, g_xx_0_yzzzz_1, g_xx_0_zzzz_1, g_xx_0_zzzzz_1, g_xxz_0_xxxxx_0, g_xxz_0_xxxxy_0, g_xxz_0_xxxxz_0, g_xxz_0_xxxyy_0, g_xxz_0_xxxyz_0, g_xxz_0_xxxzz_0, g_xxz_0_xxyyy_0, g_xxz_0_xxyyz_0, g_xxz_0_xxyzz_0, g_xxz_0_xxzzz_0, g_xxz_0_xyyyy_0, g_xxz_0_xyyyz_0, g_xxz_0_xyyzz_0, g_xxz_0_xyzzz_0, g_xxz_0_xzzzz_0, g_xxz_0_yyyyy_0, g_xxz_0_yyyyz_0, g_xxz_0_yyyzz_0, g_xxz_0_yyzzz_0, g_xxz_0_yzzzz_0, g_xxz_0_zzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxz_0_xxxxx_0[i] = g_xx_0_xxxxx_1[i] * wa_z[i];

        g_xxz_0_xxxxy_0[i] = g_xx_0_xxxxy_1[i] * wa_z[i];

        g_xxz_0_xxxxz_0[i] = g_xx_0_xxxx_1[i] * fi_acd_0 + g_xx_0_xxxxz_1[i] * wa_z[i];

        g_xxz_0_xxxyy_0[i] = g_xx_0_xxxyy_1[i] * wa_z[i];

        g_xxz_0_xxxyz_0[i] = g_xx_0_xxxy_1[i] * fi_acd_0 + g_xx_0_xxxyz_1[i] * wa_z[i];

        g_xxz_0_xxxzz_0[i] = 2.0 * g_xx_0_xxxz_1[i] * fi_acd_0 + g_xx_0_xxxzz_1[i] * wa_z[i];

        g_xxz_0_xxyyy_0[i] = g_xx_0_xxyyy_1[i] * wa_z[i];

        g_xxz_0_xxyyz_0[i] = g_xx_0_xxyy_1[i] * fi_acd_0 + g_xx_0_xxyyz_1[i] * wa_z[i];

        g_xxz_0_xxyzz_0[i] = 2.0 * g_xx_0_xxyz_1[i] * fi_acd_0 + g_xx_0_xxyzz_1[i] * wa_z[i];

        g_xxz_0_xxzzz_0[i] = 3.0 * g_xx_0_xxzz_1[i] * fi_acd_0 + g_xx_0_xxzzz_1[i] * wa_z[i];

        g_xxz_0_xyyyy_0[i] = g_xx_0_xyyyy_1[i] * wa_z[i];

        g_xxz_0_xyyyz_0[i] = g_xx_0_xyyy_1[i] * fi_acd_0 + g_xx_0_xyyyz_1[i] * wa_z[i];

        g_xxz_0_xyyzz_0[i] = 2.0 * g_xx_0_xyyz_1[i] * fi_acd_0 + g_xx_0_xyyzz_1[i] * wa_z[i];

        g_xxz_0_xyzzz_0[i] = 3.0 * g_xx_0_xyzz_1[i] * fi_acd_0 + g_xx_0_xyzzz_1[i] * wa_z[i];

        g_xxz_0_xzzzz_0[i] = 4.0 * g_xx_0_xzzz_1[i] * fi_acd_0 + g_xx_0_xzzzz_1[i] * wa_z[i];

        g_xxz_0_yyyyy_0[i] = g_xx_0_yyyyy_1[i] * wa_z[i];

        g_xxz_0_yyyyz_0[i] = g_xx_0_yyyy_1[i] * fi_acd_0 + g_xx_0_yyyyz_1[i] * wa_z[i];

        g_xxz_0_yyyzz_0[i] = 2.0 * g_xx_0_yyyz_1[i] * fi_acd_0 + g_xx_0_yyyzz_1[i] * wa_z[i];

        g_xxz_0_yyzzz_0[i] = 3.0 * g_xx_0_yyzz_1[i] * fi_acd_0 + g_xx_0_yyzzz_1[i] * wa_z[i];

        g_xxz_0_yzzzz_0[i] = 4.0 * g_xx_0_yzzz_1[i] * fi_acd_0 + g_xx_0_yzzzz_1[i] * wa_z[i];

        g_xxz_0_zzzzz_0[i] = 5.0 * g_xx_0_zzzz_1[i] * fi_acd_0 + g_xx_0_zzzzz_1[i] * wa_z[i];
    }

    /// Set up 63-84 components of targeted buffer : FSH

    auto g_xyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_fsh + 63);

    auto g_xyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_fsh + 64);

    auto g_xyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_fsh + 65);

    auto g_xyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_fsh + 66);

    auto g_xyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_fsh + 67);

    auto g_xyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_fsh + 68);

    auto g_xyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_fsh + 69);

    auto g_xyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_fsh + 70);

    auto g_xyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_fsh + 71);

    auto g_xyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_fsh + 72);

    auto g_xyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_fsh + 73);

    auto g_xyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_fsh + 74);

    auto g_xyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_fsh + 75);

    auto g_xyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_fsh + 76);

    auto g_xyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_fsh + 77);

    auto g_xyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_fsh + 78);

    auto g_xyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_fsh + 79);

    auto g_xyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_fsh + 80);

    auto g_xyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_fsh + 81);

    auto g_xyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_fsh + 82);

    auto g_xyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_fsh + 83);

    #pragma omp simd aligned(g_xyy_0_xxxxx_0, g_xyy_0_xxxxy_0, g_xyy_0_xxxxz_0, g_xyy_0_xxxyy_0, g_xyy_0_xxxyz_0, g_xyy_0_xxxzz_0, g_xyy_0_xxyyy_0, g_xyy_0_xxyyz_0, g_xyy_0_xxyzz_0, g_xyy_0_xxzzz_0, g_xyy_0_xyyyy_0, g_xyy_0_xyyyz_0, g_xyy_0_xyyzz_0, g_xyy_0_xyzzz_0, g_xyy_0_xzzzz_0, g_xyy_0_yyyyy_0, g_xyy_0_yyyyz_0, g_xyy_0_yyyzz_0, g_xyy_0_yyzzz_0, g_xyy_0_yzzzz_0, g_xyy_0_zzzzz_0, g_yy_0_xxxx_1, g_yy_0_xxxxx_1, g_yy_0_xxxxy_1, g_yy_0_xxxxz_1, g_yy_0_xxxy_1, g_yy_0_xxxyy_1, g_yy_0_xxxyz_1, g_yy_0_xxxz_1, g_yy_0_xxxzz_1, g_yy_0_xxyy_1, g_yy_0_xxyyy_1, g_yy_0_xxyyz_1, g_yy_0_xxyz_1, g_yy_0_xxyzz_1, g_yy_0_xxzz_1, g_yy_0_xxzzz_1, g_yy_0_xyyy_1, g_yy_0_xyyyy_1, g_yy_0_xyyyz_1, g_yy_0_xyyz_1, g_yy_0_xyyzz_1, g_yy_0_xyzz_1, g_yy_0_xyzzz_1, g_yy_0_xzzz_1, g_yy_0_xzzzz_1, g_yy_0_yyyy_1, g_yy_0_yyyyy_1, g_yy_0_yyyyz_1, g_yy_0_yyyz_1, g_yy_0_yyyzz_1, g_yy_0_yyzz_1, g_yy_0_yyzzz_1, g_yy_0_yzzz_1, g_yy_0_yzzzz_1, g_yy_0_zzzz_1, g_yy_0_zzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyy_0_xxxxx_0[i] = 5.0 * g_yy_0_xxxx_1[i] * fi_acd_0 + g_yy_0_xxxxx_1[i] * wa_x[i];

        g_xyy_0_xxxxy_0[i] = 4.0 * g_yy_0_xxxy_1[i] * fi_acd_0 + g_yy_0_xxxxy_1[i] * wa_x[i];

        g_xyy_0_xxxxz_0[i] = 4.0 * g_yy_0_xxxz_1[i] * fi_acd_0 + g_yy_0_xxxxz_1[i] * wa_x[i];

        g_xyy_0_xxxyy_0[i] = 3.0 * g_yy_0_xxyy_1[i] * fi_acd_0 + g_yy_0_xxxyy_1[i] * wa_x[i];

        g_xyy_0_xxxyz_0[i] = 3.0 * g_yy_0_xxyz_1[i] * fi_acd_0 + g_yy_0_xxxyz_1[i] * wa_x[i];

        g_xyy_0_xxxzz_0[i] = 3.0 * g_yy_0_xxzz_1[i] * fi_acd_0 + g_yy_0_xxxzz_1[i] * wa_x[i];

        g_xyy_0_xxyyy_0[i] = 2.0 * g_yy_0_xyyy_1[i] * fi_acd_0 + g_yy_0_xxyyy_1[i] * wa_x[i];

        g_xyy_0_xxyyz_0[i] = 2.0 * g_yy_0_xyyz_1[i] * fi_acd_0 + g_yy_0_xxyyz_1[i] * wa_x[i];

        g_xyy_0_xxyzz_0[i] = 2.0 * g_yy_0_xyzz_1[i] * fi_acd_0 + g_yy_0_xxyzz_1[i] * wa_x[i];

        g_xyy_0_xxzzz_0[i] = 2.0 * g_yy_0_xzzz_1[i] * fi_acd_0 + g_yy_0_xxzzz_1[i] * wa_x[i];

        g_xyy_0_xyyyy_0[i] = g_yy_0_yyyy_1[i] * fi_acd_0 + g_yy_0_xyyyy_1[i] * wa_x[i];

        g_xyy_0_xyyyz_0[i] = g_yy_0_yyyz_1[i] * fi_acd_0 + g_yy_0_xyyyz_1[i] * wa_x[i];

        g_xyy_0_xyyzz_0[i] = g_yy_0_yyzz_1[i] * fi_acd_0 + g_yy_0_xyyzz_1[i] * wa_x[i];

        g_xyy_0_xyzzz_0[i] = g_yy_0_yzzz_1[i] * fi_acd_0 + g_yy_0_xyzzz_1[i] * wa_x[i];

        g_xyy_0_xzzzz_0[i] = g_yy_0_zzzz_1[i] * fi_acd_0 + g_yy_0_xzzzz_1[i] * wa_x[i];

        g_xyy_0_yyyyy_0[i] = g_yy_0_yyyyy_1[i] * wa_x[i];

        g_xyy_0_yyyyz_0[i] = g_yy_0_yyyyz_1[i] * wa_x[i];

        g_xyy_0_yyyzz_0[i] = g_yy_0_yyyzz_1[i] * wa_x[i];

        g_xyy_0_yyzzz_0[i] = g_yy_0_yyzzz_1[i] * wa_x[i];

        g_xyy_0_yzzzz_0[i] = g_yy_0_yzzzz_1[i] * wa_x[i];

        g_xyy_0_zzzzz_0[i] = g_yy_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 84-105 components of targeted buffer : FSH

    auto g_xyz_0_xxxxx_0 = pbuffer.data(idx_eri_0_fsh + 84);

    auto g_xyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_fsh + 85);

    auto g_xyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_fsh + 86);

    auto g_xyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_fsh + 87);

    auto g_xyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_fsh + 88);

    auto g_xyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_fsh + 89);

    auto g_xyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_fsh + 90);

    auto g_xyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_fsh + 91);

    auto g_xyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_fsh + 92);

    auto g_xyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_fsh + 93);

    auto g_xyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_fsh + 94);

    auto g_xyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_fsh + 95);

    auto g_xyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_fsh + 96);

    auto g_xyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_fsh + 97);

    auto g_xyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_fsh + 98);

    auto g_xyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_fsh + 99);

    auto g_xyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_fsh + 100);

    auto g_xyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_fsh + 101);

    auto g_xyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_fsh + 102);

    auto g_xyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_fsh + 103);

    auto g_xyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_fsh + 104);

    #pragma omp simd aligned(g_xy_0_xxxxy_1, g_xy_0_xxxyy_1, g_xy_0_xxyyy_1, g_xy_0_xyyyy_1, g_xyz_0_xxxxx_0, g_xyz_0_xxxxy_0, g_xyz_0_xxxxz_0, g_xyz_0_xxxyy_0, g_xyz_0_xxxyz_0, g_xyz_0_xxxzz_0, g_xyz_0_xxyyy_0, g_xyz_0_xxyyz_0, g_xyz_0_xxyzz_0, g_xyz_0_xxzzz_0, g_xyz_0_xyyyy_0, g_xyz_0_xyyyz_0, g_xyz_0_xyyzz_0, g_xyz_0_xyzzz_0, g_xyz_0_xzzzz_0, g_xyz_0_yyyyy_0, g_xyz_0_yyyyz_0, g_xyz_0_yyyzz_0, g_xyz_0_yyzzz_0, g_xyz_0_yzzzz_0, g_xyz_0_zzzzz_0, g_xz_0_xxxxx_1, g_xz_0_xxxxz_1, g_xz_0_xxxzz_1, g_xz_0_xxzzz_1, g_xz_0_xzzzz_1, g_yz_0_xxxyz_1, g_yz_0_xxyyz_1, g_yz_0_xxyz_1, g_yz_0_xxyzz_1, g_yz_0_xyyyz_1, g_yz_0_xyyz_1, g_yz_0_xyyzz_1, g_yz_0_xyzz_1, g_yz_0_xyzzz_1, g_yz_0_yyyyy_1, g_yz_0_yyyyz_1, g_yz_0_yyyz_1, g_yz_0_yyyzz_1, g_yz_0_yyzz_1, g_yz_0_yyzzz_1, g_yz_0_yzzz_1, g_yz_0_yzzzz_1, g_yz_0_zzzzz_1, wa_x, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyz_0_xxxxx_0[i] = g_xz_0_xxxxx_1[i] * wa_y[i];

        g_xyz_0_xxxxy_0[i] = g_xy_0_xxxxy_1[i] * wa_z[i];

        g_xyz_0_xxxxz_0[i] = g_xz_0_xxxxz_1[i] * wa_y[i];

        g_xyz_0_xxxyy_0[i] = g_xy_0_xxxyy_1[i] * wa_z[i];

        g_xyz_0_xxxyz_0[i] = 3.0 * g_yz_0_xxyz_1[i] * fi_acd_0 + g_yz_0_xxxyz_1[i] * wa_x[i];

        g_xyz_0_xxxzz_0[i] = g_xz_0_xxxzz_1[i] * wa_y[i];

        g_xyz_0_xxyyy_0[i] = g_xy_0_xxyyy_1[i] * wa_z[i];

        g_xyz_0_xxyyz_0[i] = 2.0 * g_yz_0_xyyz_1[i] * fi_acd_0 + g_yz_0_xxyyz_1[i] * wa_x[i];

        g_xyz_0_xxyzz_0[i] = 2.0 * g_yz_0_xyzz_1[i] * fi_acd_0 + g_yz_0_xxyzz_1[i] * wa_x[i];

        g_xyz_0_xxzzz_0[i] = g_xz_0_xxzzz_1[i] * wa_y[i];

        g_xyz_0_xyyyy_0[i] = g_xy_0_xyyyy_1[i] * wa_z[i];

        g_xyz_0_xyyyz_0[i] = g_yz_0_yyyz_1[i] * fi_acd_0 + g_yz_0_xyyyz_1[i] * wa_x[i];

        g_xyz_0_xyyzz_0[i] = g_yz_0_yyzz_1[i] * fi_acd_0 + g_yz_0_xyyzz_1[i] * wa_x[i];

        g_xyz_0_xyzzz_0[i] = g_yz_0_yzzz_1[i] * fi_acd_0 + g_yz_0_xyzzz_1[i] * wa_x[i];

        g_xyz_0_xzzzz_0[i] = g_xz_0_xzzzz_1[i] * wa_y[i];

        g_xyz_0_yyyyy_0[i] = g_yz_0_yyyyy_1[i] * wa_x[i];

        g_xyz_0_yyyyz_0[i] = g_yz_0_yyyyz_1[i] * wa_x[i];

        g_xyz_0_yyyzz_0[i] = g_yz_0_yyyzz_1[i] * wa_x[i];

        g_xyz_0_yyzzz_0[i] = g_yz_0_yyzzz_1[i] * wa_x[i];

        g_xyz_0_yzzzz_0[i] = g_yz_0_yzzzz_1[i] * wa_x[i];

        g_xyz_0_zzzzz_0[i] = g_yz_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 105-126 components of targeted buffer : FSH

    auto g_xzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_fsh + 105);

    auto g_xzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_fsh + 106);

    auto g_xzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_fsh + 107);

    auto g_xzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_fsh + 108);

    auto g_xzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_fsh + 109);

    auto g_xzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_fsh + 110);

    auto g_xzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_fsh + 111);

    auto g_xzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_fsh + 112);

    auto g_xzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_fsh + 113);

    auto g_xzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_fsh + 114);

    auto g_xzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_fsh + 115);

    auto g_xzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_fsh + 116);

    auto g_xzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_fsh + 117);

    auto g_xzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_fsh + 118);

    auto g_xzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_fsh + 119);

    auto g_xzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_fsh + 120);

    auto g_xzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_fsh + 121);

    auto g_xzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_fsh + 122);

    auto g_xzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_fsh + 123);

    auto g_xzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_fsh + 124);

    auto g_xzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_fsh + 125);

    #pragma omp simd aligned(g_xzz_0_xxxxx_0, g_xzz_0_xxxxy_0, g_xzz_0_xxxxz_0, g_xzz_0_xxxyy_0, g_xzz_0_xxxyz_0, g_xzz_0_xxxzz_0, g_xzz_0_xxyyy_0, g_xzz_0_xxyyz_0, g_xzz_0_xxyzz_0, g_xzz_0_xxzzz_0, g_xzz_0_xyyyy_0, g_xzz_0_xyyyz_0, g_xzz_0_xyyzz_0, g_xzz_0_xyzzz_0, g_xzz_0_xzzzz_0, g_xzz_0_yyyyy_0, g_xzz_0_yyyyz_0, g_xzz_0_yyyzz_0, g_xzz_0_yyzzz_0, g_xzz_0_yzzzz_0, g_xzz_0_zzzzz_0, g_zz_0_xxxx_1, g_zz_0_xxxxx_1, g_zz_0_xxxxy_1, g_zz_0_xxxxz_1, g_zz_0_xxxy_1, g_zz_0_xxxyy_1, g_zz_0_xxxyz_1, g_zz_0_xxxz_1, g_zz_0_xxxzz_1, g_zz_0_xxyy_1, g_zz_0_xxyyy_1, g_zz_0_xxyyz_1, g_zz_0_xxyz_1, g_zz_0_xxyzz_1, g_zz_0_xxzz_1, g_zz_0_xxzzz_1, g_zz_0_xyyy_1, g_zz_0_xyyyy_1, g_zz_0_xyyyz_1, g_zz_0_xyyz_1, g_zz_0_xyyzz_1, g_zz_0_xyzz_1, g_zz_0_xyzzz_1, g_zz_0_xzzz_1, g_zz_0_xzzzz_1, g_zz_0_yyyy_1, g_zz_0_yyyyy_1, g_zz_0_yyyyz_1, g_zz_0_yyyz_1, g_zz_0_yyyzz_1, g_zz_0_yyzz_1, g_zz_0_yyzzz_1, g_zz_0_yzzz_1, g_zz_0_yzzzz_1, g_zz_0_zzzz_1, g_zz_0_zzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzz_0_xxxxx_0[i] = 5.0 * g_zz_0_xxxx_1[i] * fi_acd_0 + g_zz_0_xxxxx_1[i] * wa_x[i];

        g_xzz_0_xxxxy_0[i] = 4.0 * g_zz_0_xxxy_1[i] * fi_acd_0 + g_zz_0_xxxxy_1[i] * wa_x[i];

        g_xzz_0_xxxxz_0[i] = 4.0 * g_zz_0_xxxz_1[i] * fi_acd_0 + g_zz_0_xxxxz_1[i] * wa_x[i];

        g_xzz_0_xxxyy_0[i] = 3.0 * g_zz_0_xxyy_1[i] * fi_acd_0 + g_zz_0_xxxyy_1[i] * wa_x[i];

        g_xzz_0_xxxyz_0[i] = 3.0 * g_zz_0_xxyz_1[i] * fi_acd_0 + g_zz_0_xxxyz_1[i] * wa_x[i];

        g_xzz_0_xxxzz_0[i] = 3.0 * g_zz_0_xxzz_1[i] * fi_acd_0 + g_zz_0_xxxzz_1[i] * wa_x[i];

        g_xzz_0_xxyyy_0[i] = 2.0 * g_zz_0_xyyy_1[i] * fi_acd_0 + g_zz_0_xxyyy_1[i] * wa_x[i];

        g_xzz_0_xxyyz_0[i] = 2.0 * g_zz_0_xyyz_1[i] * fi_acd_0 + g_zz_0_xxyyz_1[i] * wa_x[i];

        g_xzz_0_xxyzz_0[i] = 2.0 * g_zz_0_xyzz_1[i] * fi_acd_0 + g_zz_0_xxyzz_1[i] * wa_x[i];

        g_xzz_0_xxzzz_0[i] = 2.0 * g_zz_0_xzzz_1[i] * fi_acd_0 + g_zz_0_xxzzz_1[i] * wa_x[i];

        g_xzz_0_xyyyy_0[i] = g_zz_0_yyyy_1[i] * fi_acd_0 + g_zz_0_xyyyy_1[i] * wa_x[i];

        g_xzz_0_xyyyz_0[i] = g_zz_0_yyyz_1[i] * fi_acd_0 + g_zz_0_xyyyz_1[i] * wa_x[i];

        g_xzz_0_xyyzz_0[i] = g_zz_0_yyzz_1[i] * fi_acd_0 + g_zz_0_xyyzz_1[i] * wa_x[i];

        g_xzz_0_xyzzz_0[i] = g_zz_0_yzzz_1[i] * fi_acd_0 + g_zz_0_xyzzz_1[i] * wa_x[i];

        g_xzz_0_xzzzz_0[i] = g_zz_0_zzzz_1[i] * fi_acd_0 + g_zz_0_xzzzz_1[i] * wa_x[i];

        g_xzz_0_yyyyy_0[i] = g_zz_0_yyyyy_1[i] * wa_x[i];

        g_xzz_0_yyyyz_0[i] = g_zz_0_yyyyz_1[i] * wa_x[i];

        g_xzz_0_yyyzz_0[i] = g_zz_0_yyyzz_1[i] * wa_x[i];

        g_xzz_0_yyzzz_0[i] = g_zz_0_yyzzz_1[i] * wa_x[i];

        g_xzz_0_yzzzz_0[i] = g_zz_0_yzzzz_1[i] * wa_x[i];

        g_xzz_0_zzzzz_0[i] = g_zz_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 126-147 components of targeted buffer : FSH

    auto g_yyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_fsh + 126);

    auto g_yyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_fsh + 127);

    auto g_yyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_fsh + 128);

    auto g_yyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_fsh + 129);

    auto g_yyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_fsh + 130);

    auto g_yyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_fsh + 131);

    auto g_yyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_fsh + 132);

    auto g_yyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_fsh + 133);

    auto g_yyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_fsh + 134);

    auto g_yyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_fsh + 135);

    auto g_yyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_fsh + 136);

    auto g_yyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_fsh + 137);

    auto g_yyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_fsh + 138);

    auto g_yyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_fsh + 139);

    auto g_yyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_fsh + 140);

    auto g_yyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_fsh + 141);

    auto g_yyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_fsh + 142);

    auto g_yyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_fsh + 143);

    auto g_yyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_fsh + 144);

    auto g_yyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_fsh + 145);

    auto g_yyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_fsh + 146);

    #pragma omp simd aligned(g_y_0_xxxxx_0, g_y_0_xxxxx_1, g_y_0_xxxxy_0, g_y_0_xxxxy_1, g_y_0_xxxxz_0, g_y_0_xxxxz_1, g_y_0_xxxyy_0, g_y_0_xxxyy_1, g_y_0_xxxyz_0, g_y_0_xxxyz_1, g_y_0_xxxzz_0, g_y_0_xxxzz_1, g_y_0_xxyyy_0, g_y_0_xxyyy_1, g_y_0_xxyyz_0, g_y_0_xxyyz_1, g_y_0_xxyzz_0, g_y_0_xxyzz_1, g_y_0_xxzzz_0, g_y_0_xxzzz_1, g_y_0_xyyyy_0, g_y_0_xyyyy_1, g_y_0_xyyyz_0, g_y_0_xyyyz_1, g_y_0_xyyzz_0, g_y_0_xyyzz_1, g_y_0_xyzzz_0, g_y_0_xyzzz_1, g_y_0_xzzzz_0, g_y_0_xzzzz_1, g_y_0_yyyyy_0, g_y_0_yyyyy_1, g_y_0_yyyyz_0, g_y_0_yyyyz_1, g_y_0_yyyzz_0, g_y_0_yyyzz_1, g_y_0_yyzzz_0, g_y_0_yyzzz_1, g_y_0_yzzzz_0, g_y_0_yzzzz_1, g_y_0_zzzzz_0, g_y_0_zzzzz_1, g_yy_0_xxxx_1, g_yy_0_xxxxx_1, g_yy_0_xxxxy_1, g_yy_0_xxxxz_1, g_yy_0_xxxy_1, g_yy_0_xxxyy_1, g_yy_0_xxxyz_1, g_yy_0_xxxz_1, g_yy_0_xxxzz_1, g_yy_0_xxyy_1, g_yy_0_xxyyy_1, g_yy_0_xxyyz_1, g_yy_0_xxyz_1, g_yy_0_xxyzz_1, g_yy_0_xxzz_1, g_yy_0_xxzzz_1, g_yy_0_xyyy_1, g_yy_0_xyyyy_1, g_yy_0_xyyyz_1, g_yy_0_xyyz_1, g_yy_0_xyyzz_1, g_yy_0_xyzz_1, g_yy_0_xyzzz_1, g_yy_0_xzzz_1, g_yy_0_xzzzz_1, g_yy_0_yyyy_1, g_yy_0_yyyyy_1, g_yy_0_yyyyz_1, g_yy_0_yyyz_1, g_yy_0_yyyzz_1, g_yy_0_yyzz_1, g_yy_0_yyzzz_1, g_yy_0_yzzz_1, g_yy_0_yzzzz_1, g_yy_0_zzzz_1, g_yy_0_zzzzz_1, g_yyy_0_xxxxx_0, g_yyy_0_xxxxy_0, g_yyy_0_xxxxz_0, g_yyy_0_xxxyy_0, g_yyy_0_xxxyz_0, g_yyy_0_xxxzz_0, g_yyy_0_xxyyy_0, g_yyy_0_xxyyz_0, g_yyy_0_xxyzz_0, g_yyy_0_xxzzz_0, g_yyy_0_xyyyy_0, g_yyy_0_xyyyz_0, g_yyy_0_xyyzz_0, g_yyy_0_xyzzz_0, g_yyy_0_xzzzz_0, g_yyy_0_yyyyy_0, g_yyy_0_yyyyz_0, g_yyy_0_yyyzz_0, g_yyy_0_yyzzz_0, g_yyy_0_yzzzz_0, g_yyy_0_zzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyy_0_xxxxx_0[i] = 2.0 * g_y_0_xxxxx_0[i] * fbe_0 - 2.0 * g_y_0_xxxxx_1[i] * fz_be_0 + g_yy_0_xxxxx_1[i] * wa_y[i];

        g_yyy_0_xxxxy_0[i] = 2.0 * g_y_0_xxxxy_0[i] * fbe_0 - 2.0 * g_y_0_xxxxy_1[i] * fz_be_0 + g_yy_0_xxxx_1[i] * fi_acd_0 + g_yy_0_xxxxy_1[i] * wa_y[i];

        g_yyy_0_xxxxz_0[i] = 2.0 * g_y_0_xxxxz_0[i] * fbe_0 - 2.0 * g_y_0_xxxxz_1[i] * fz_be_0 + g_yy_0_xxxxz_1[i] * wa_y[i];

        g_yyy_0_xxxyy_0[i] = 2.0 * g_y_0_xxxyy_0[i] * fbe_0 - 2.0 * g_y_0_xxxyy_1[i] * fz_be_0 + 2.0 * g_yy_0_xxxy_1[i] * fi_acd_0 + g_yy_0_xxxyy_1[i] * wa_y[i];

        g_yyy_0_xxxyz_0[i] = 2.0 * g_y_0_xxxyz_0[i] * fbe_0 - 2.0 * g_y_0_xxxyz_1[i] * fz_be_0 + g_yy_0_xxxz_1[i] * fi_acd_0 + g_yy_0_xxxyz_1[i] * wa_y[i];

        g_yyy_0_xxxzz_0[i] = 2.0 * g_y_0_xxxzz_0[i] * fbe_0 - 2.0 * g_y_0_xxxzz_1[i] * fz_be_0 + g_yy_0_xxxzz_1[i] * wa_y[i];

        g_yyy_0_xxyyy_0[i] = 2.0 * g_y_0_xxyyy_0[i] * fbe_0 - 2.0 * g_y_0_xxyyy_1[i] * fz_be_0 + 3.0 * g_yy_0_xxyy_1[i] * fi_acd_0 + g_yy_0_xxyyy_1[i] * wa_y[i];

        g_yyy_0_xxyyz_0[i] = 2.0 * g_y_0_xxyyz_0[i] * fbe_0 - 2.0 * g_y_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_yy_0_xxyz_1[i] * fi_acd_0 + g_yy_0_xxyyz_1[i] * wa_y[i];

        g_yyy_0_xxyzz_0[i] = 2.0 * g_y_0_xxyzz_0[i] * fbe_0 - 2.0 * g_y_0_xxyzz_1[i] * fz_be_0 + g_yy_0_xxzz_1[i] * fi_acd_0 + g_yy_0_xxyzz_1[i] * wa_y[i];

        g_yyy_0_xxzzz_0[i] = 2.0 * g_y_0_xxzzz_0[i] * fbe_0 - 2.0 * g_y_0_xxzzz_1[i] * fz_be_0 + g_yy_0_xxzzz_1[i] * wa_y[i];

        g_yyy_0_xyyyy_0[i] = 2.0 * g_y_0_xyyyy_0[i] * fbe_0 - 2.0 * g_y_0_xyyyy_1[i] * fz_be_0 + 4.0 * g_yy_0_xyyy_1[i] * fi_acd_0 + g_yy_0_xyyyy_1[i] * wa_y[i];

        g_yyy_0_xyyyz_0[i] = 2.0 * g_y_0_xyyyz_0[i] * fbe_0 - 2.0 * g_y_0_xyyyz_1[i] * fz_be_0 + 3.0 * g_yy_0_xyyz_1[i] * fi_acd_0 + g_yy_0_xyyyz_1[i] * wa_y[i];

        g_yyy_0_xyyzz_0[i] = 2.0 * g_y_0_xyyzz_0[i] * fbe_0 - 2.0 * g_y_0_xyyzz_1[i] * fz_be_0 + 2.0 * g_yy_0_xyzz_1[i] * fi_acd_0 + g_yy_0_xyyzz_1[i] * wa_y[i];

        g_yyy_0_xyzzz_0[i] = 2.0 * g_y_0_xyzzz_0[i] * fbe_0 - 2.0 * g_y_0_xyzzz_1[i] * fz_be_0 + g_yy_0_xzzz_1[i] * fi_acd_0 + g_yy_0_xyzzz_1[i] * wa_y[i];

        g_yyy_0_xzzzz_0[i] = 2.0 * g_y_0_xzzzz_0[i] * fbe_0 - 2.0 * g_y_0_xzzzz_1[i] * fz_be_0 + g_yy_0_xzzzz_1[i] * wa_y[i];

        g_yyy_0_yyyyy_0[i] = 2.0 * g_y_0_yyyyy_0[i] * fbe_0 - 2.0 * g_y_0_yyyyy_1[i] * fz_be_0 + 5.0 * g_yy_0_yyyy_1[i] * fi_acd_0 + g_yy_0_yyyyy_1[i] * wa_y[i];

        g_yyy_0_yyyyz_0[i] = 2.0 * g_y_0_yyyyz_0[i] * fbe_0 - 2.0 * g_y_0_yyyyz_1[i] * fz_be_0 + 4.0 * g_yy_0_yyyz_1[i] * fi_acd_0 + g_yy_0_yyyyz_1[i] * wa_y[i];

        g_yyy_0_yyyzz_0[i] = 2.0 * g_y_0_yyyzz_0[i] * fbe_0 - 2.0 * g_y_0_yyyzz_1[i] * fz_be_0 + 3.0 * g_yy_0_yyzz_1[i] * fi_acd_0 + g_yy_0_yyyzz_1[i] * wa_y[i];

        g_yyy_0_yyzzz_0[i] = 2.0 * g_y_0_yyzzz_0[i] * fbe_0 - 2.0 * g_y_0_yyzzz_1[i] * fz_be_0 + 2.0 * g_yy_0_yzzz_1[i] * fi_acd_0 + g_yy_0_yyzzz_1[i] * wa_y[i];

        g_yyy_0_yzzzz_0[i] = 2.0 * g_y_0_yzzzz_0[i] * fbe_0 - 2.0 * g_y_0_yzzzz_1[i] * fz_be_0 + g_yy_0_zzzz_1[i] * fi_acd_0 + g_yy_0_yzzzz_1[i] * wa_y[i];

        g_yyy_0_zzzzz_0[i] = 2.0 * g_y_0_zzzzz_0[i] * fbe_0 - 2.0 * g_y_0_zzzzz_1[i] * fz_be_0 + g_yy_0_zzzzz_1[i] * wa_y[i];
    }

    /// Set up 147-168 components of targeted buffer : FSH

    auto g_yyz_0_xxxxx_0 = pbuffer.data(idx_eri_0_fsh + 147);

    auto g_yyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_fsh + 148);

    auto g_yyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_fsh + 149);

    auto g_yyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_fsh + 150);

    auto g_yyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_fsh + 151);

    auto g_yyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_fsh + 152);

    auto g_yyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_fsh + 153);

    auto g_yyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_fsh + 154);

    auto g_yyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_fsh + 155);

    auto g_yyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_fsh + 156);

    auto g_yyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_fsh + 157);

    auto g_yyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_fsh + 158);

    auto g_yyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_fsh + 159);

    auto g_yyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_fsh + 160);

    auto g_yyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_fsh + 161);

    auto g_yyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_fsh + 162);

    auto g_yyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_fsh + 163);

    auto g_yyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_fsh + 164);

    auto g_yyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_fsh + 165);

    auto g_yyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_fsh + 166);

    auto g_yyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_fsh + 167);

    #pragma omp simd aligned(g_yy_0_xxxx_1, g_yy_0_xxxxx_1, g_yy_0_xxxxy_1, g_yy_0_xxxxz_1, g_yy_0_xxxy_1, g_yy_0_xxxyy_1, g_yy_0_xxxyz_1, g_yy_0_xxxz_1, g_yy_0_xxxzz_1, g_yy_0_xxyy_1, g_yy_0_xxyyy_1, g_yy_0_xxyyz_1, g_yy_0_xxyz_1, g_yy_0_xxyzz_1, g_yy_0_xxzz_1, g_yy_0_xxzzz_1, g_yy_0_xyyy_1, g_yy_0_xyyyy_1, g_yy_0_xyyyz_1, g_yy_0_xyyz_1, g_yy_0_xyyzz_1, g_yy_0_xyzz_1, g_yy_0_xyzzz_1, g_yy_0_xzzz_1, g_yy_0_xzzzz_1, g_yy_0_yyyy_1, g_yy_0_yyyyy_1, g_yy_0_yyyyz_1, g_yy_0_yyyz_1, g_yy_0_yyyzz_1, g_yy_0_yyzz_1, g_yy_0_yyzzz_1, g_yy_0_yzzz_1, g_yy_0_yzzzz_1, g_yy_0_zzzz_1, g_yy_0_zzzzz_1, g_yyz_0_xxxxx_0, g_yyz_0_xxxxy_0, g_yyz_0_xxxxz_0, g_yyz_0_xxxyy_0, g_yyz_0_xxxyz_0, g_yyz_0_xxxzz_0, g_yyz_0_xxyyy_0, g_yyz_0_xxyyz_0, g_yyz_0_xxyzz_0, g_yyz_0_xxzzz_0, g_yyz_0_xyyyy_0, g_yyz_0_xyyyz_0, g_yyz_0_xyyzz_0, g_yyz_0_xyzzz_0, g_yyz_0_xzzzz_0, g_yyz_0_yyyyy_0, g_yyz_0_yyyyz_0, g_yyz_0_yyyzz_0, g_yyz_0_yyzzz_0, g_yyz_0_yzzzz_0, g_yyz_0_zzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyz_0_xxxxx_0[i] = g_yy_0_xxxxx_1[i] * wa_z[i];

        g_yyz_0_xxxxy_0[i] = g_yy_0_xxxxy_1[i] * wa_z[i];

        g_yyz_0_xxxxz_0[i] = g_yy_0_xxxx_1[i] * fi_acd_0 + g_yy_0_xxxxz_1[i] * wa_z[i];

        g_yyz_0_xxxyy_0[i] = g_yy_0_xxxyy_1[i] * wa_z[i];

        g_yyz_0_xxxyz_0[i] = g_yy_0_xxxy_1[i] * fi_acd_0 + g_yy_0_xxxyz_1[i] * wa_z[i];

        g_yyz_0_xxxzz_0[i] = 2.0 * g_yy_0_xxxz_1[i] * fi_acd_0 + g_yy_0_xxxzz_1[i] * wa_z[i];

        g_yyz_0_xxyyy_0[i] = g_yy_0_xxyyy_1[i] * wa_z[i];

        g_yyz_0_xxyyz_0[i] = g_yy_0_xxyy_1[i] * fi_acd_0 + g_yy_0_xxyyz_1[i] * wa_z[i];

        g_yyz_0_xxyzz_0[i] = 2.0 * g_yy_0_xxyz_1[i] * fi_acd_0 + g_yy_0_xxyzz_1[i] * wa_z[i];

        g_yyz_0_xxzzz_0[i] = 3.0 * g_yy_0_xxzz_1[i] * fi_acd_0 + g_yy_0_xxzzz_1[i] * wa_z[i];

        g_yyz_0_xyyyy_0[i] = g_yy_0_xyyyy_1[i] * wa_z[i];

        g_yyz_0_xyyyz_0[i] = g_yy_0_xyyy_1[i] * fi_acd_0 + g_yy_0_xyyyz_1[i] * wa_z[i];

        g_yyz_0_xyyzz_0[i] = 2.0 * g_yy_0_xyyz_1[i] * fi_acd_0 + g_yy_0_xyyzz_1[i] * wa_z[i];

        g_yyz_0_xyzzz_0[i] = 3.0 * g_yy_0_xyzz_1[i] * fi_acd_0 + g_yy_0_xyzzz_1[i] * wa_z[i];

        g_yyz_0_xzzzz_0[i] = 4.0 * g_yy_0_xzzz_1[i] * fi_acd_0 + g_yy_0_xzzzz_1[i] * wa_z[i];

        g_yyz_0_yyyyy_0[i] = g_yy_0_yyyyy_1[i] * wa_z[i];

        g_yyz_0_yyyyz_0[i] = g_yy_0_yyyy_1[i] * fi_acd_0 + g_yy_0_yyyyz_1[i] * wa_z[i];

        g_yyz_0_yyyzz_0[i] = 2.0 * g_yy_0_yyyz_1[i] * fi_acd_0 + g_yy_0_yyyzz_1[i] * wa_z[i];

        g_yyz_0_yyzzz_0[i] = 3.0 * g_yy_0_yyzz_1[i] * fi_acd_0 + g_yy_0_yyzzz_1[i] * wa_z[i];

        g_yyz_0_yzzzz_0[i] = 4.0 * g_yy_0_yzzz_1[i] * fi_acd_0 + g_yy_0_yzzzz_1[i] * wa_z[i];

        g_yyz_0_zzzzz_0[i] = 5.0 * g_yy_0_zzzz_1[i] * fi_acd_0 + g_yy_0_zzzzz_1[i] * wa_z[i];
    }

    /// Set up 168-189 components of targeted buffer : FSH

    auto g_yzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_fsh + 168);

    auto g_yzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_fsh + 169);

    auto g_yzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_fsh + 170);

    auto g_yzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_fsh + 171);

    auto g_yzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_fsh + 172);

    auto g_yzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_fsh + 173);

    auto g_yzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_fsh + 174);

    auto g_yzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_fsh + 175);

    auto g_yzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_fsh + 176);

    auto g_yzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_fsh + 177);

    auto g_yzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_fsh + 178);

    auto g_yzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_fsh + 179);

    auto g_yzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_fsh + 180);

    auto g_yzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_fsh + 181);

    auto g_yzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_fsh + 182);

    auto g_yzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_fsh + 183);

    auto g_yzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_fsh + 184);

    auto g_yzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_fsh + 185);

    auto g_yzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_fsh + 186);

    auto g_yzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_fsh + 187);

    auto g_yzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_fsh + 188);

    #pragma omp simd aligned(g_yzz_0_xxxxx_0, g_yzz_0_xxxxy_0, g_yzz_0_xxxxz_0, g_yzz_0_xxxyy_0, g_yzz_0_xxxyz_0, g_yzz_0_xxxzz_0, g_yzz_0_xxyyy_0, g_yzz_0_xxyyz_0, g_yzz_0_xxyzz_0, g_yzz_0_xxzzz_0, g_yzz_0_xyyyy_0, g_yzz_0_xyyyz_0, g_yzz_0_xyyzz_0, g_yzz_0_xyzzz_0, g_yzz_0_xzzzz_0, g_yzz_0_yyyyy_0, g_yzz_0_yyyyz_0, g_yzz_0_yyyzz_0, g_yzz_0_yyzzz_0, g_yzz_0_yzzzz_0, g_yzz_0_zzzzz_0, g_zz_0_xxxx_1, g_zz_0_xxxxx_1, g_zz_0_xxxxy_1, g_zz_0_xxxxz_1, g_zz_0_xxxy_1, g_zz_0_xxxyy_1, g_zz_0_xxxyz_1, g_zz_0_xxxz_1, g_zz_0_xxxzz_1, g_zz_0_xxyy_1, g_zz_0_xxyyy_1, g_zz_0_xxyyz_1, g_zz_0_xxyz_1, g_zz_0_xxyzz_1, g_zz_0_xxzz_1, g_zz_0_xxzzz_1, g_zz_0_xyyy_1, g_zz_0_xyyyy_1, g_zz_0_xyyyz_1, g_zz_0_xyyz_1, g_zz_0_xyyzz_1, g_zz_0_xyzz_1, g_zz_0_xyzzz_1, g_zz_0_xzzz_1, g_zz_0_xzzzz_1, g_zz_0_yyyy_1, g_zz_0_yyyyy_1, g_zz_0_yyyyz_1, g_zz_0_yyyz_1, g_zz_0_yyyzz_1, g_zz_0_yyzz_1, g_zz_0_yyzzz_1, g_zz_0_yzzz_1, g_zz_0_yzzzz_1, g_zz_0_zzzz_1, g_zz_0_zzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzz_0_xxxxx_0[i] = g_zz_0_xxxxx_1[i] * wa_y[i];

        g_yzz_0_xxxxy_0[i] = g_zz_0_xxxx_1[i] * fi_acd_0 + g_zz_0_xxxxy_1[i] * wa_y[i];

        g_yzz_0_xxxxz_0[i] = g_zz_0_xxxxz_1[i] * wa_y[i];

        g_yzz_0_xxxyy_0[i] = 2.0 * g_zz_0_xxxy_1[i] * fi_acd_0 + g_zz_0_xxxyy_1[i] * wa_y[i];

        g_yzz_0_xxxyz_0[i] = g_zz_0_xxxz_1[i] * fi_acd_0 + g_zz_0_xxxyz_1[i] * wa_y[i];

        g_yzz_0_xxxzz_0[i] = g_zz_0_xxxzz_1[i] * wa_y[i];

        g_yzz_0_xxyyy_0[i] = 3.0 * g_zz_0_xxyy_1[i] * fi_acd_0 + g_zz_0_xxyyy_1[i] * wa_y[i];

        g_yzz_0_xxyyz_0[i] = 2.0 * g_zz_0_xxyz_1[i] * fi_acd_0 + g_zz_0_xxyyz_1[i] * wa_y[i];

        g_yzz_0_xxyzz_0[i] = g_zz_0_xxzz_1[i] * fi_acd_0 + g_zz_0_xxyzz_1[i] * wa_y[i];

        g_yzz_0_xxzzz_0[i] = g_zz_0_xxzzz_1[i] * wa_y[i];

        g_yzz_0_xyyyy_0[i] = 4.0 * g_zz_0_xyyy_1[i] * fi_acd_0 + g_zz_0_xyyyy_1[i] * wa_y[i];

        g_yzz_0_xyyyz_0[i] = 3.0 * g_zz_0_xyyz_1[i] * fi_acd_0 + g_zz_0_xyyyz_1[i] * wa_y[i];

        g_yzz_0_xyyzz_0[i] = 2.0 * g_zz_0_xyzz_1[i] * fi_acd_0 + g_zz_0_xyyzz_1[i] * wa_y[i];

        g_yzz_0_xyzzz_0[i] = g_zz_0_xzzz_1[i] * fi_acd_0 + g_zz_0_xyzzz_1[i] * wa_y[i];

        g_yzz_0_xzzzz_0[i] = g_zz_0_xzzzz_1[i] * wa_y[i];

        g_yzz_0_yyyyy_0[i] = 5.0 * g_zz_0_yyyy_1[i] * fi_acd_0 + g_zz_0_yyyyy_1[i] * wa_y[i];

        g_yzz_0_yyyyz_0[i] = 4.0 * g_zz_0_yyyz_1[i] * fi_acd_0 + g_zz_0_yyyyz_1[i] * wa_y[i];

        g_yzz_0_yyyzz_0[i] = 3.0 * g_zz_0_yyzz_1[i] * fi_acd_0 + g_zz_0_yyyzz_1[i] * wa_y[i];

        g_yzz_0_yyzzz_0[i] = 2.0 * g_zz_0_yzzz_1[i] * fi_acd_0 + g_zz_0_yyzzz_1[i] * wa_y[i];

        g_yzz_0_yzzzz_0[i] = g_zz_0_zzzz_1[i] * fi_acd_0 + g_zz_0_yzzzz_1[i] * wa_y[i];

        g_yzz_0_zzzzz_0[i] = g_zz_0_zzzzz_1[i] * wa_y[i];
    }

    /// Set up 189-210 components of targeted buffer : FSH

    auto g_zzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_fsh + 189);

    auto g_zzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_fsh + 190);

    auto g_zzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_fsh + 191);

    auto g_zzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_fsh + 192);

    auto g_zzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_fsh + 193);

    auto g_zzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_fsh + 194);

    auto g_zzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_fsh + 195);

    auto g_zzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_fsh + 196);

    auto g_zzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_fsh + 197);

    auto g_zzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_fsh + 198);

    auto g_zzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_fsh + 199);

    auto g_zzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_fsh + 200);

    auto g_zzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_fsh + 201);

    auto g_zzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_fsh + 202);

    auto g_zzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_fsh + 203);

    auto g_zzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_fsh + 204);

    auto g_zzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_fsh + 205);

    auto g_zzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_fsh + 206);

    auto g_zzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_fsh + 207);

    auto g_zzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_fsh + 208);

    auto g_zzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_fsh + 209);

    #pragma omp simd aligned(g_z_0_xxxxx_0, g_z_0_xxxxx_1, g_z_0_xxxxy_0, g_z_0_xxxxy_1, g_z_0_xxxxz_0, g_z_0_xxxxz_1, g_z_0_xxxyy_0, g_z_0_xxxyy_1, g_z_0_xxxyz_0, g_z_0_xxxyz_1, g_z_0_xxxzz_0, g_z_0_xxxzz_1, g_z_0_xxyyy_0, g_z_0_xxyyy_1, g_z_0_xxyyz_0, g_z_0_xxyyz_1, g_z_0_xxyzz_0, g_z_0_xxyzz_1, g_z_0_xxzzz_0, g_z_0_xxzzz_1, g_z_0_xyyyy_0, g_z_0_xyyyy_1, g_z_0_xyyyz_0, g_z_0_xyyyz_1, g_z_0_xyyzz_0, g_z_0_xyyzz_1, g_z_0_xyzzz_0, g_z_0_xyzzz_1, g_z_0_xzzzz_0, g_z_0_xzzzz_1, g_z_0_yyyyy_0, g_z_0_yyyyy_1, g_z_0_yyyyz_0, g_z_0_yyyyz_1, g_z_0_yyyzz_0, g_z_0_yyyzz_1, g_z_0_yyzzz_0, g_z_0_yyzzz_1, g_z_0_yzzzz_0, g_z_0_yzzzz_1, g_z_0_zzzzz_0, g_z_0_zzzzz_1, g_zz_0_xxxx_1, g_zz_0_xxxxx_1, g_zz_0_xxxxy_1, g_zz_0_xxxxz_1, g_zz_0_xxxy_1, g_zz_0_xxxyy_1, g_zz_0_xxxyz_1, g_zz_0_xxxz_1, g_zz_0_xxxzz_1, g_zz_0_xxyy_1, g_zz_0_xxyyy_1, g_zz_0_xxyyz_1, g_zz_0_xxyz_1, g_zz_0_xxyzz_1, g_zz_0_xxzz_1, g_zz_0_xxzzz_1, g_zz_0_xyyy_1, g_zz_0_xyyyy_1, g_zz_0_xyyyz_1, g_zz_0_xyyz_1, g_zz_0_xyyzz_1, g_zz_0_xyzz_1, g_zz_0_xyzzz_1, g_zz_0_xzzz_1, g_zz_0_xzzzz_1, g_zz_0_yyyy_1, g_zz_0_yyyyy_1, g_zz_0_yyyyz_1, g_zz_0_yyyz_1, g_zz_0_yyyzz_1, g_zz_0_yyzz_1, g_zz_0_yyzzz_1, g_zz_0_yzzz_1, g_zz_0_yzzzz_1, g_zz_0_zzzz_1, g_zz_0_zzzzz_1, g_zzz_0_xxxxx_0, g_zzz_0_xxxxy_0, g_zzz_0_xxxxz_0, g_zzz_0_xxxyy_0, g_zzz_0_xxxyz_0, g_zzz_0_xxxzz_0, g_zzz_0_xxyyy_0, g_zzz_0_xxyyz_0, g_zzz_0_xxyzz_0, g_zzz_0_xxzzz_0, g_zzz_0_xyyyy_0, g_zzz_0_xyyyz_0, g_zzz_0_xyyzz_0, g_zzz_0_xyzzz_0, g_zzz_0_xzzzz_0, g_zzz_0_yyyyy_0, g_zzz_0_yyyyz_0, g_zzz_0_yyyzz_0, g_zzz_0_yyzzz_0, g_zzz_0_yzzzz_0, g_zzz_0_zzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzz_0_xxxxx_0[i] = 2.0 * g_z_0_xxxxx_0[i] * fbe_0 - 2.0 * g_z_0_xxxxx_1[i] * fz_be_0 + g_zz_0_xxxxx_1[i] * wa_z[i];

        g_zzz_0_xxxxy_0[i] = 2.0 * g_z_0_xxxxy_0[i] * fbe_0 - 2.0 * g_z_0_xxxxy_1[i] * fz_be_0 + g_zz_0_xxxxy_1[i] * wa_z[i];

        g_zzz_0_xxxxz_0[i] = 2.0 * g_z_0_xxxxz_0[i] * fbe_0 - 2.0 * g_z_0_xxxxz_1[i] * fz_be_0 + g_zz_0_xxxx_1[i] * fi_acd_0 + g_zz_0_xxxxz_1[i] * wa_z[i];

        g_zzz_0_xxxyy_0[i] = 2.0 * g_z_0_xxxyy_0[i] * fbe_0 - 2.0 * g_z_0_xxxyy_1[i] * fz_be_0 + g_zz_0_xxxyy_1[i] * wa_z[i];

        g_zzz_0_xxxyz_0[i] = 2.0 * g_z_0_xxxyz_0[i] * fbe_0 - 2.0 * g_z_0_xxxyz_1[i] * fz_be_0 + g_zz_0_xxxy_1[i] * fi_acd_0 + g_zz_0_xxxyz_1[i] * wa_z[i];

        g_zzz_0_xxxzz_0[i] = 2.0 * g_z_0_xxxzz_0[i] * fbe_0 - 2.0 * g_z_0_xxxzz_1[i] * fz_be_0 + 2.0 * g_zz_0_xxxz_1[i] * fi_acd_0 + g_zz_0_xxxzz_1[i] * wa_z[i];

        g_zzz_0_xxyyy_0[i] = 2.0 * g_z_0_xxyyy_0[i] * fbe_0 - 2.0 * g_z_0_xxyyy_1[i] * fz_be_0 + g_zz_0_xxyyy_1[i] * wa_z[i];

        g_zzz_0_xxyyz_0[i] = 2.0 * g_z_0_xxyyz_0[i] * fbe_0 - 2.0 * g_z_0_xxyyz_1[i] * fz_be_0 + g_zz_0_xxyy_1[i] * fi_acd_0 + g_zz_0_xxyyz_1[i] * wa_z[i];

        g_zzz_0_xxyzz_0[i] = 2.0 * g_z_0_xxyzz_0[i] * fbe_0 - 2.0 * g_z_0_xxyzz_1[i] * fz_be_0 + 2.0 * g_zz_0_xxyz_1[i] * fi_acd_0 + g_zz_0_xxyzz_1[i] * wa_z[i];

        g_zzz_0_xxzzz_0[i] = 2.0 * g_z_0_xxzzz_0[i] * fbe_0 - 2.0 * g_z_0_xxzzz_1[i] * fz_be_0 + 3.0 * g_zz_0_xxzz_1[i] * fi_acd_0 + g_zz_0_xxzzz_1[i] * wa_z[i];

        g_zzz_0_xyyyy_0[i] = 2.0 * g_z_0_xyyyy_0[i] * fbe_0 - 2.0 * g_z_0_xyyyy_1[i] * fz_be_0 + g_zz_0_xyyyy_1[i] * wa_z[i];

        g_zzz_0_xyyyz_0[i] = 2.0 * g_z_0_xyyyz_0[i] * fbe_0 - 2.0 * g_z_0_xyyyz_1[i] * fz_be_0 + g_zz_0_xyyy_1[i] * fi_acd_0 + g_zz_0_xyyyz_1[i] * wa_z[i];

        g_zzz_0_xyyzz_0[i] = 2.0 * g_z_0_xyyzz_0[i] * fbe_0 - 2.0 * g_z_0_xyyzz_1[i] * fz_be_0 + 2.0 * g_zz_0_xyyz_1[i] * fi_acd_0 + g_zz_0_xyyzz_1[i] * wa_z[i];

        g_zzz_0_xyzzz_0[i] = 2.0 * g_z_0_xyzzz_0[i] * fbe_0 - 2.0 * g_z_0_xyzzz_1[i] * fz_be_0 + 3.0 * g_zz_0_xyzz_1[i] * fi_acd_0 + g_zz_0_xyzzz_1[i] * wa_z[i];

        g_zzz_0_xzzzz_0[i] = 2.0 * g_z_0_xzzzz_0[i] * fbe_0 - 2.0 * g_z_0_xzzzz_1[i] * fz_be_0 + 4.0 * g_zz_0_xzzz_1[i] * fi_acd_0 + g_zz_0_xzzzz_1[i] * wa_z[i];

        g_zzz_0_yyyyy_0[i] = 2.0 * g_z_0_yyyyy_0[i] * fbe_0 - 2.0 * g_z_0_yyyyy_1[i] * fz_be_0 + g_zz_0_yyyyy_1[i] * wa_z[i];

        g_zzz_0_yyyyz_0[i] = 2.0 * g_z_0_yyyyz_0[i] * fbe_0 - 2.0 * g_z_0_yyyyz_1[i] * fz_be_0 + g_zz_0_yyyy_1[i] * fi_acd_0 + g_zz_0_yyyyz_1[i] * wa_z[i];

        g_zzz_0_yyyzz_0[i] = 2.0 * g_z_0_yyyzz_0[i] * fbe_0 - 2.0 * g_z_0_yyyzz_1[i] * fz_be_0 + 2.0 * g_zz_0_yyyz_1[i] * fi_acd_0 + g_zz_0_yyyzz_1[i] * wa_z[i];

        g_zzz_0_yyzzz_0[i] = 2.0 * g_z_0_yyzzz_0[i] * fbe_0 - 2.0 * g_z_0_yyzzz_1[i] * fz_be_0 + 3.0 * g_zz_0_yyzz_1[i] * fi_acd_0 + g_zz_0_yyzzz_1[i] * wa_z[i];

        g_zzz_0_yzzzz_0[i] = 2.0 * g_z_0_yzzzz_0[i] * fbe_0 - 2.0 * g_z_0_yzzzz_1[i] * fz_be_0 + 4.0 * g_zz_0_yzzz_1[i] * fi_acd_0 + g_zz_0_yzzzz_1[i] * wa_z[i];

        g_zzz_0_zzzzz_0[i] = 2.0 * g_z_0_zzzzz_0[i] * fbe_0 - 2.0 * g_z_0_zzzzz_1[i] * fz_be_0 + 5.0 * g_zz_0_zzzz_1[i] * fi_acd_0 + g_zz_0_zzzzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace


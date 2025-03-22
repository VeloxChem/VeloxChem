#include "ThreeCenterElectronRepulsionPrimRecISD.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_isd(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_isd,
                                 size_t idx_eri_0_gsd,
                                 size_t idx_eri_1_gsd,
                                 size_t idx_eri_1_hsp,
                                 size_t idx_eri_1_hsd,
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

    /// Set up components of auxilary buffer : GSD

    auto g_xxxx_0_xx_0 = pbuffer.data(idx_eri_0_gsd);

    auto g_xxxx_0_xy_0 = pbuffer.data(idx_eri_0_gsd + 1);

    auto g_xxxx_0_xz_0 = pbuffer.data(idx_eri_0_gsd + 2);

    auto g_xxxx_0_yy_0 = pbuffer.data(idx_eri_0_gsd + 3);

    auto g_xxxx_0_yz_0 = pbuffer.data(idx_eri_0_gsd + 4);

    auto g_xxxx_0_zz_0 = pbuffer.data(idx_eri_0_gsd + 5);

    auto g_xxxy_0_xx_0 = pbuffer.data(idx_eri_0_gsd + 6);

    auto g_xxxy_0_xz_0 = pbuffer.data(idx_eri_0_gsd + 8);

    auto g_xxxz_0_xx_0 = pbuffer.data(idx_eri_0_gsd + 12);

    auto g_xxxz_0_xy_0 = pbuffer.data(idx_eri_0_gsd + 13);

    auto g_xxyy_0_xx_0 = pbuffer.data(idx_eri_0_gsd + 18);

    auto g_xxyy_0_xy_0 = pbuffer.data(idx_eri_0_gsd + 19);

    auto g_xxyy_0_xz_0 = pbuffer.data(idx_eri_0_gsd + 20);

    auto g_xxyy_0_yy_0 = pbuffer.data(idx_eri_0_gsd + 21);

    auto g_xxyy_0_yz_0 = pbuffer.data(idx_eri_0_gsd + 22);

    auto g_xxyy_0_zz_0 = pbuffer.data(idx_eri_0_gsd + 23);

    auto g_xxzz_0_xx_0 = pbuffer.data(idx_eri_0_gsd + 30);

    auto g_xxzz_0_xy_0 = pbuffer.data(idx_eri_0_gsd + 31);

    auto g_xxzz_0_xz_0 = pbuffer.data(idx_eri_0_gsd + 32);

    auto g_xxzz_0_yy_0 = pbuffer.data(idx_eri_0_gsd + 33);

    auto g_xxzz_0_yz_0 = pbuffer.data(idx_eri_0_gsd + 34);

    auto g_xxzz_0_zz_0 = pbuffer.data(idx_eri_0_gsd + 35);

    auto g_xyyy_0_xy_0 = pbuffer.data(idx_eri_0_gsd + 37);

    auto g_xyyy_0_yy_0 = pbuffer.data(idx_eri_0_gsd + 39);

    auto g_xyyy_0_yz_0 = pbuffer.data(idx_eri_0_gsd + 40);

    auto g_xyyy_0_zz_0 = pbuffer.data(idx_eri_0_gsd + 41);

    auto g_xzzz_0_xz_0 = pbuffer.data(idx_eri_0_gsd + 56);

    auto g_xzzz_0_yy_0 = pbuffer.data(idx_eri_0_gsd + 57);

    auto g_xzzz_0_yz_0 = pbuffer.data(idx_eri_0_gsd + 58);

    auto g_xzzz_0_zz_0 = pbuffer.data(idx_eri_0_gsd + 59);

    auto g_yyyy_0_xx_0 = pbuffer.data(idx_eri_0_gsd + 60);

    auto g_yyyy_0_xy_0 = pbuffer.data(idx_eri_0_gsd + 61);

    auto g_yyyy_0_xz_0 = pbuffer.data(idx_eri_0_gsd + 62);

    auto g_yyyy_0_yy_0 = pbuffer.data(idx_eri_0_gsd + 63);

    auto g_yyyy_0_yz_0 = pbuffer.data(idx_eri_0_gsd + 64);

    auto g_yyyy_0_zz_0 = pbuffer.data(idx_eri_0_gsd + 65);

    auto g_yyyz_0_xy_0 = pbuffer.data(idx_eri_0_gsd + 67);

    auto g_yyyz_0_yy_0 = pbuffer.data(idx_eri_0_gsd + 69);

    auto g_yyzz_0_xx_0 = pbuffer.data(idx_eri_0_gsd + 72);

    auto g_yyzz_0_xy_0 = pbuffer.data(idx_eri_0_gsd + 73);

    auto g_yyzz_0_xz_0 = pbuffer.data(idx_eri_0_gsd + 74);

    auto g_yyzz_0_yy_0 = pbuffer.data(idx_eri_0_gsd + 75);

    auto g_yyzz_0_yz_0 = pbuffer.data(idx_eri_0_gsd + 76);

    auto g_yyzz_0_zz_0 = pbuffer.data(idx_eri_0_gsd + 77);

    auto g_yzzz_0_xx_0 = pbuffer.data(idx_eri_0_gsd + 78);

    auto g_yzzz_0_xz_0 = pbuffer.data(idx_eri_0_gsd + 80);

    auto g_yzzz_0_yz_0 = pbuffer.data(idx_eri_0_gsd + 82);

    auto g_yzzz_0_zz_0 = pbuffer.data(idx_eri_0_gsd + 83);

    auto g_zzzz_0_xx_0 = pbuffer.data(idx_eri_0_gsd + 84);

    auto g_zzzz_0_xy_0 = pbuffer.data(idx_eri_0_gsd + 85);

    auto g_zzzz_0_xz_0 = pbuffer.data(idx_eri_0_gsd + 86);

    auto g_zzzz_0_yy_0 = pbuffer.data(idx_eri_0_gsd + 87);

    auto g_zzzz_0_yz_0 = pbuffer.data(idx_eri_0_gsd + 88);

    auto g_zzzz_0_zz_0 = pbuffer.data(idx_eri_0_gsd + 89);

    /// Set up components of auxilary buffer : GSD

    auto g_xxxx_0_xx_1 = pbuffer.data(idx_eri_1_gsd);

    auto g_xxxx_0_xy_1 = pbuffer.data(idx_eri_1_gsd + 1);

    auto g_xxxx_0_xz_1 = pbuffer.data(idx_eri_1_gsd + 2);

    auto g_xxxx_0_yy_1 = pbuffer.data(idx_eri_1_gsd + 3);

    auto g_xxxx_0_yz_1 = pbuffer.data(idx_eri_1_gsd + 4);

    auto g_xxxx_0_zz_1 = pbuffer.data(idx_eri_1_gsd + 5);

    auto g_xxxy_0_xx_1 = pbuffer.data(idx_eri_1_gsd + 6);

    auto g_xxxy_0_xz_1 = pbuffer.data(idx_eri_1_gsd + 8);

    auto g_xxxz_0_xx_1 = pbuffer.data(idx_eri_1_gsd + 12);

    auto g_xxxz_0_xy_1 = pbuffer.data(idx_eri_1_gsd + 13);

    auto g_xxyy_0_xx_1 = pbuffer.data(idx_eri_1_gsd + 18);

    auto g_xxyy_0_xy_1 = pbuffer.data(idx_eri_1_gsd + 19);

    auto g_xxyy_0_xz_1 = pbuffer.data(idx_eri_1_gsd + 20);

    auto g_xxyy_0_yy_1 = pbuffer.data(idx_eri_1_gsd + 21);

    auto g_xxyy_0_yz_1 = pbuffer.data(idx_eri_1_gsd + 22);

    auto g_xxyy_0_zz_1 = pbuffer.data(idx_eri_1_gsd + 23);

    auto g_xxzz_0_xx_1 = pbuffer.data(idx_eri_1_gsd + 30);

    auto g_xxzz_0_xy_1 = pbuffer.data(idx_eri_1_gsd + 31);

    auto g_xxzz_0_xz_1 = pbuffer.data(idx_eri_1_gsd + 32);

    auto g_xxzz_0_yy_1 = pbuffer.data(idx_eri_1_gsd + 33);

    auto g_xxzz_0_yz_1 = pbuffer.data(idx_eri_1_gsd + 34);

    auto g_xxzz_0_zz_1 = pbuffer.data(idx_eri_1_gsd + 35);

    auto g_xyyy_0_xy_1 = pbuffer.data(idx_eri_1_gsd + 37);

    auto g_xyyy_0_yy_1 = pbuffer.data(idx_eri_1_gsd + 39);

    auto g_xyyy_0_yz_1 = pbuffer.data(idx_eri_1_gsd + 40);

    auto g_xyyy_0_zz_1 = pbuffer.data(idx_eri_1_gsd + 41);

    auto g_xzzz_0_xz_1 = pbuffer.data(idx_eri_1_gsd + 56);

    auto g_xzzz_0_yy_1 = pbuffer.data(idx_eri_1_gsd + 57);

    auto g_xzzz_0_yz_1 = pbuffer.data(idx_eri_1_gsd + 58);

    auto g_xzzz_0_zz_1 = pbuffer.data(idx_eri_1_gsd + 59);

    auto g_yyyy_0_xx_1 = pbuffer.data(idx_eri_1_gsd + 60);

    auto g_yyyy_0_xy_1 = pbuffer.data(idx_eri_1_gsd + 61);

    auto g_yyyy_0_xz_1 = pbuffer.data(idx_eri_1_gsd + 62);

    auto g_yyyy_0_yy_1 = pbuffer.data(idx_eri_1_gsd + 63);

    auto g_yyyy_0_yz_1 = pbuffer.data(idx_eri_1_gsd + 64);

    auto g_yyyy_0_zz_1 = pbuffer.data(idx_eri_1_gsd + 65);

    auto g_yyyz_0_xy_1 = pbuffer.data(idx_eri_1_gsd + 67);

    auto g_yyyz_0_yy_1 = pbuffer.data(idx_eri_1_gsd + 69);

    auto g_yyzz_0_xx_1 = pbuffer.data(idx_eri_1_gsd + 72);

    auto g_yyzz_0_xy_1 = pbuffer.data(idx_eri_1_gsd + 73);

    auto g_yyzz_0_xz_1 = pbuffer.data(idx_eri_1_gsd + 74);

    auto g_yyzz_0_yy_1 = pbuffer.data(idx_eri_1_gsd + 75);

    auto g_yyzz_0_yz_1 = pbuffer.data(idx_eri_1_gsd + 76);

    auto g_yyzz_0_zz_1 = pbuffer.data(idx_eri_1_gsd + 77);

    auto g_yzzz_0_xx_1 = pbuffer.data(idx_eri_1_gsd + 78);

    auto g_yzzz_0_xz_1 = pbuffer.data(idx_eri_1_gsd + 80);

    auto g_yzzz_0_yz_1 = pbuffer.data(idx_eri_1_gsd + 82);

    auto g_yzzz_0_zz_1 = pbuffer.data(idx_eri_1_gsd + 83);

    auto g_zzzz_0_xx_1 = pbuffer.data(idx_eri_1_gsd + 84);

    auto g_zzzz_0_xy_1 = pbuffer.data(idx_eri_1_gsd + 85);

    auto g_zzzz_0_xz_1 = pbuffer.data(idx_eri_1_gsd + 86);

    auto g_zzzz_0_yy_1 = pbuffer.data(idx_eri_1_gsd + 87);

    auto g_zzzz_0_yz_1 = pbuffer.data(idx_eri_1_gsd + 88);

    auto g_zzzz_0_zz_1 = pbuffer.data(idx_eri_1_gsd + 89);

    /// Set up components of auxilary buffer : HSP

    auto g_xxxxx_0_x_1 = pbuffer.data(idx_eri_1_hsp);

    auto g_xxxxx_0_y_1 = pbuffer.data(idx_eri_1_hsp + 1);

    auto g_xxxxx_0_z_1 = pbuffer.data(idx_eri_1_hsp + 2);

    auto g_xxxxz_0_z_1 = pbuffer.data(idx_eri_1_hsp + 8);

    auto g_xxxyy_0_x_1 = pbuffer.data(idx_eri_1_hsp + 9);

    auto g_xxxyy_0_y_1 = pbuffer.data(idx_eri_1_hsp + 10);

    auto g_xxxyy_0_z_1 = pbuffer.data(idx_eri_1_hsp + 11);

    auto g_xxxzz_0_x_1 = pbuffer.data(idx_eri_1_hsp + 15);

    auto g_xxxzz_0_y_1 = pbuffer.data(idx_eri_1_hsp + 16);

    auto g_xxxzz_0_z_1 = pbuffer.data(idx_eri_1_hsp + 17);

    auto g_xxyyy_0_x_1 = pbuffer.data(idx_eri_1_hsp + 18);

    auto g_xxyyy_0_y_1 = pbuffer.data(idx_eri_1_hsp + 19);

    auto g_xxyyy_0_z_1 = pbuffer.data(idx_eri_1_hsp + 20);

    auto g_xxzzz_0_x_1 = pbuffer.data(idx_eri_1_hsp + 27);

    auto g_xxzzz_0_y_1 = pbuffer.data(idx_eri_1_hsp + 28);

    auto g_xxzzz_0_z_1 = pbuffer.data(idx_eri_1_hsp + 29);

    auto g_xyyyy_0_y_1 = pbuffer.data(idx_eri_1_hsp + 31);

    auto g_xzzzz_0_z_1 = pbuffer.data(idx_eri_1_hsp + 44);

    auto g_yyyyy_0_x_1 = pbuffer.data(idx_eri_1_hsp + 45);

    auto g_yyyyy_0_y_1 = pbuffer.data(idx_eri_1_hsp + 46);

    auto g_yyyyy_0_z_1 = pbuffer.data(idx_eri_1_hsp + 47);

    auto g_yyyyz_0_z_1 = pbuffer.data(idx_eri_1_hsp + 50);

    auto g_yyyzz_0_x_1 = pbuffer.data(idx_eri_1_hsp + 51);

    auto g_yyyzz_0_y_1 = pbuffer.data(idx_eri_1_hsp + 52);

    auto g_yyyzz_0_z_1 = pbuffer.data(idx_eri_1_hsp + 53);

    auto g_yyzzz_0_x_1 = pbuffer.data(idx_eri_1_hsp + 54);

    auto g_yyzzz_0_y_1 = pbuffer.data(idx_eri_1_hsp + 55);

    auto g_yyzzz_0_z_1 = pbuffer.data(idx_eri_1_hsp + 56);

    auto g_yzzzz_0_y_1 = pbuffer.data(idx_eri_1_hsp + 58);

    auto g_yzzzz_0_z_1 = pbuffer.data(idx_eri_1_hsp + 59);

    auto g_zzzzz_0_x_1 = pbuffer.data(idx_eri_1_hsp + 60);

    auto g_zzzzz_0_y_1 = pbuffer.data(idx_eri_1_hsp + 61);

    auto g_zzzzz_0_z_1 = pbuffer.data(idx_eri_1_hsp + 62);

    /// Set up components of auxilary buffer : HSD

    auto g_xxxxx_0_xx_1 = pbuffer.data(idx_eri_1_hsd);

    auto g_xxxxx_0_xy_1 = pbuffer.data(idx_eri_1_hsd + 1);

    auto g_xxxxx_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 2);

    auto g_xxxxx_0_yy_1 = pbuffer.data(idx_eri_1_hsd + 3);

    auto g_xxxxx_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 4);

    auto g_xxxxx_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 5);

    auto g_xxxxy_0_xx_1 = pbuffer.data(idx_eri_1_hsd + 6);

    auto g_xxxxy_0_xy_1 = pbuffer.data(idx_eri_1_hsd + 7);

    auto g_xxxxy_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 8);

    auto g_xxxxy_0_yy_1 = pbuffer.data(idx_eri_1_hsd + 9);

    auto g_xxxxz_0_xx_1 = pbuffer.data(idx_eri_1_hsd + 12);

    auto g_xxxxz_0_xy_1 = pbuffer.data(idx_eri_1_hsd + 13);

    auto g_xxxxz_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 14);

    auto g_xxxxz_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 16);

    auto g_xxxxz_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 17);

    auto g_xxxyy_0_xx_1 = pbuffer.data(idx_eri_1_hsd + 18);

    auto g_xxxyy_0_xy_1 = pbuffer.data(idx_eri_1_hsd + 19);

    auto g_xxxyy_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 20);

    auto g_xxxyy_0_yy_1 = pbuffer.data(idx_eri_1_hsd + 21);

    auto g_xxxyy_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 22);

    auto g_xxxyy_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 23);

    auto g_xxxzz_0_xx_1 = pbuffer.data(idx_eri_1_hsd + 30);

    auto g_xxxzz_0_xy_1 = pbuffer.data(idx_eri_1_hsd + 31);

    auto g_xxxzz_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 32);

    auto g_xxxzz_0_yy_1 = pbuffer.data(idx_eri_1_hsd + 33);

    auto g_xxxzz_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 34);

    auto g_xxxzz_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 35);

    auto g_xxyyy_0_xx_1 = pbuffer.data(idx_eri_1_hsd + 36);

    auto g_xxyyy_0_xy_1 = pbuffer.data(idx_eri_1_hsd + 37);

    auto g_xxyyy_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 38);

    auto g_xxyyy_0_yy_1 = pbuffer.data(idx_eri_1_hsd + 39);

    auto g_xxyyy_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 40);

    auto g_xxyyy_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 41);

    auto g_xxyyz_0_xy_1 = pbuffer.data(idx_eri_1_hsd + 43);

    auto g_xxyzz_0_xx_1 = pbuffer.data(idx_eri_1_hsd + 48);

    auto g_xxyzz_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 50);

    auto g_xxzzz_0_xx_1 = pbuffer.data(idx_eri_1_hsd + 54);

    auto g_xxzzz_0_xy_1 = pbuffer.data(idx_eri_1_hsd + 55);

    auto g_xxzzz_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 56);

    auto g_xxzzz_0_yy_1 = pbuffer.data(idx_eri_1_hsd + 57);

    auto g_xxzzz_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 58);

    auto g_xxzzz_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 59);

    auto g_xyyyy_0_xx_1 = pbuffer.data(idx_eri_1_hsd + 60);

    auto g_xyyyy_0_xy_1 = pbuffer.data(idx_eri_1_hsd + 61);

    auto g_xyyyy_0_yy_1 = pbuffer.data(idx_eri_1_hsd + 63);

    auto g_xyyyy_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 64);

    auto g_xyyyy_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 65);

    auto g_xyyzz_0_yy_1 = pbuffer.data(idx_eri_1_hsd + 75);

    auto g_xyyzz_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 76);

    auto g_xyyzz_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 77);

    auto g_xzzzz_0_xx_1 = pbuffer.data(idx_eri_1_hsd + 84);

    auto g_xzzzz_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 86);

    auto g_xzzzz_0_yy_1 = pbuffer.data(idx_eri_1_hsd + 87);

    auto g_xzzzz_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 88);

    auto g_xzzzz_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 89);

    auto g_yyyyy_0_xx_1 = pbuffer.data(idx_eri_1_hsd + 90);

    auto g_yyyyy_0_xy_1 = pbuffer.data(idx_eri_1_hsd + 91);

    auto g_yyyyy_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 92);

    auto g_yyyyy_0_yy_1 = pbuffer.data(idx_eri_1_hsd + 93);

    auto g_yyyyy_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 94);

    auto g_yyyyy_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 95);

    auto g_yyyyz_0_xy_1 = pbuffer.data(idx_eri_1_hsd + 97);

    auto g_yyyyz_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 98);

    auto g_yyyyz_0_yy_1 = pbuffer.data(idx_eri_1_hsd + 99);

    auto g_yyyyz_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 100);

    auto g_yyyyz_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 101);

    auto g_yyyzz_0_xx_1 = pbuffer.data(idx_eri_1_hsd + 102);

    auto g_yyyzz_0_xy_1 = pbuffer.data(idx_eri_1_hsd + 103);

    auto g_yyyzz_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 104);

    auto g_yyyzz_0_yy_1 = pbuffer.data(idx_eri_1_hsd + 105);

    auto g_yyyzz_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 106);

    auto g_yyyzz_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 107);

    auto g_yyzzz_0_xx_1 = pbuffer.data(idx_eri_1_hsd + 108);

    auto g_yyzzz_0_xy_1 = pbuffer.data(idx_eri_1_hsd + 109);

    auto g_yyzzz_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 110);

    auto g_yyzzz_0_yy_1 = pbuffer.data(idx_eri_1_hsd + 111);

    auto g_yyzzz_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 112);

    auto g_yyzzz_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 113);

    auto g_yzzzz_0_xx_1 = pbuffer.data(idx_eri_1_hsd + 114);

    auto g_yzzzz_0_xy_1 = pbuffer.data(idx_eri_1_hsd + 115);

    auto g_yzzzz_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 116);

    auto g_yzzzz_0_yy_1 = pbuffer.data(idx_eri_1_hsd + 117);

    auto g_yzzzz_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 118);

    auto g_yzzzz_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 119);

    auto g_zzzzz_0_xx_1 = pbuffer.data(idx_eri_1_hsd + 120);

    auto g_zzzzz_0_xy_1 = pbuffer.data(idx_eri_1_hsd + 121);

    auto g_zzzzz_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 122);

    auto g_zzzzz_0_yy_1 = pbuffer.data(idx_eri_1_hsd + 123);

    auto g_zzzzz_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 124);

    auto g_zzzzz_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 125);

    /// Set up 0-6 components of targeted buffer : ISD

    auto g_xxxxxx_0_xx_0 = pbuffer.data(idx_eri_0_isd);

    auto g_xxxxxx_0_xy_0 = pbuffer.data(idx_eri_0_isd + 1);

    auto g_xxxxxx_0_xz_0 = pbuffer.data(idx_eri_0_isd + 2);

    auto g_xxxxxx_0_yy_0 = pbuffer.data(idx_eri_0_isd + 3);

    auto g_xxxxxx_0_yz_0 = pbuffer.data(idx_eri_0_isd + 4);

    auto g_xxxxxx_0_zz_0 = pbuffer.data(idx_eri_0_isd + 5);

    #pragma omp simd aligned(g_xxxx_0_xx_0, g_xxxx_0_xx_1, g_xxxx_0_xy_0, g_xxxx_0_xy_1, g_xxxx_0_xz_0, g_xxxx_0_xz_1, g_xxxx_0_yy_0, g_xxxx_0_yy_1, g_xxxx_0_yz_0, g_xxxx_0_yz_1, g_xxxx_0_zz_0, g_xxxx_0_zz_1, g_xxxxx_0_x_1, g_xxxxx_0_xx_1, g_xxxxx_0_xy_1, g_xxxxx_0_xz_1, g_xxxxx_0_y_1, g_xxxxx_0_yy_1, g_xxxxx_0_yz_1, g_xxxxx_0_z_1, g_xxxxx_0_zz_1, g_xxxxxx_0_xx_0, g_xxxxxx_0_xy_0, g_xxxxxx_0_xz_0, g_xxxxxx_0_yy_0, g_xxxxxx_0_yz_0, g_xxxxxx_0_zz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxxx_0_xx_0[i] = 5.0 * g_xxxx_0_xx_0[i] * fbe_0 - 5.0 * g_xxxx_0_xx_1[i] * fz_be_0 + 2.0 * g_xxxxx_0_x_1[i] * fi_acd_0 + g_xxxxx_0_xx_1[i] * wa_x[i];

        g_xxxxxx_0_xy_0[i] = 5.0 * g_xxxx_0_xy_0[i] * fbe_0 - 5.0 * g_xxxx_0_xy_1[i] * fz_be_0 + g_xxxxx_0_y_1[i] * fi_acd_0 + g_xxxxx_0_xy_1[i] * wa_x[i];

        g_xxxxxx_0_xz_0[i] = 5.0 * g_xxxx_0_xz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xz_1[i] * fz_be_0 + g_xxxxx_0_z_1[i] * fi_acd_0 + g_xxxxx_0_xz_1[i] * wa_x[i];

        g_xxxxxx_0_yy_0[i] = 5.0 * g_xxxx_0_yy_0[i] * fbe_0 - 5.0 * g_xxxx_0_yy_1[i] * fz_be_0 + g_xxxxx_0_yy_1[i] * wa_x[i];

        g_xxxxxx_0_yz_0[i] = 5.0 * g_xxxx_0_yz_0[i] * fbe_0 - 5.0 * g_xxxx_0_yz_1[i] * fz_be_0 + g_xxxxx_0_yz_1[i] * wa_x[i];

        g_xxxxxx_0_zz_0[i] = 5.0 * g_xxxx_0_zz_0[i] * fbe_0 - 5.0 * g_xxxx_0_zz_1[i] * fz_be_0 + g_xxxxx_0_zz_1[i] * wa_x[i];
    }

    /// Set up 6-12 components of targeted buffer : ISD

    auto g_xxxxxy_0_xx_0 = pbuffer.data(idx_eri_0_isd + 6);

    auto g_xxxxxy_0_xy_0 = pbuffer.data(idx_eri_0_isd + 7);

    auto g_xxxxxy_0_xz_0 = pbuffer.data(idx_eri_0_isd + 8);

    auto g_xxxxxy_0_yy_0 = pbuffer.data(idx_eri_0_isd + 9);

    auto g_xxxxxy_0_yz_0 = pbuffer.data(idx_eri_0_isd + 10);

    auto g_xxxxxy_0_zz_0 = pbuffer.data(idx_eri_0_isd + 11);

    #pragma omp simd aligned(g_xxxxx_0_x_1, g_xxxxx_0_xx_1, g_xxxxx_0_xy_1, g_xxxxx_0_xz_1, g_xxxxx_0_y_1, g_xxxxx_0_yy_1, g_xxxxx_0_yz_1, g_xxxxx_0_z_1, g_xxxxx_0_zz_1, g_xxxxxy_0_xx_0, g_xxxxxy_0_xy_0, g_xxxxxy_0_xz_0, g_xxxxxy_0_yy_0, g_xxxxxy_0_yz_0, g_xxxxxy_0_zz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxy_0_xx_0[i] = g_xxxxx_0_xx_1[i] * wa_y[i];

        g_xxxxxy_0_xy_0[i] = g_xxxxx_0_x_1[i] * fi_acd_0 + g_xxxxx_0_xy_1[i] * wa_y[i];

        g_xxxxxy_0_xz_0[i] = g_xxxxx_0_xz_1[i] * wa_y[i];

        g_xxxxxy_0_yy_0[i] = 2.0 * g_xxxxx_0_y_1[i] * fi_acd_0 + g_xxxxx_0_yy_1[i] * wa_y[i];

        g_xxxxxy_0_yz_0[i] = g_xxxxx_0_z_1[i] * fi_acd_0 + g_xxxxx_0_yz_1[i] * wa_y[i];

        g_xxxxxy_0_zz_0[i] = g_xxxxx_0_zz_1[i] * wa_y[i];
    }

    /// Set up 12-18 components of targeted buffer : ISD

    auto g_xxxxxz_0_xx_0 = pbuffer.data(idx_eri_0_isd + 12);

    auto g_xxxxxz_0_xy_0 = pbuffer.data(idx_eri_0_isd + 13);

    auto g_xxxxxz_0_xz_0 = pbuffer.data(idx_eri_0_isd + 14);

    auto g_xxxxxz_0_yy_0 = pbuffer.data(idx_eri_0_isd + 15);

    auto g_xxxxxz_0_yz_0 = pbuffer.data(idx_eri_0_isd + 16);

    auto g_xxxxxz_0_zz_0 = pbuffer.data(idx_eri_0_isd + 17);

    #pragma omp simd aligned(g_xxxxx_0_x_1, g_xxxxx_0_xx_1, g_xxxxx_0_xy_1, g_xxxxx_0_xz_1, g_xxxxx_0_y_1, g_xxxxx_0_yy_1, g_xxxxx_0_yz_1, g_xxxxx_0_z_1, g_xxxxx_0_zz_1, g_xxxxxz_0_xx_0, g_xxxxxz_0_xy_0, g_xxxxxz_0_xz_0, g_xxxxxz_0_yy_0, g_xxxxxz_0_yz_0, g_xxxxxz_0_zz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxz_0_xx_0[i] = g_xxxxx_0_xx_1[i] * wa_z[i];

        g_xxxxxz_0_xy_0[i] = g_xxxxx_0_xy_1[i] * wa_z[i];

        g_xxxxxz_0_xz_0[i] = g_xxxxx_0_x_1[i] * fi_acd_0 + g_xxxxx_0_xz_1[i] * wa_z[i];

        g_xxxxxz_0_yy_0[i] = g_xxxxx_0_yy_1[i] * wa_z[i];

        g_xxxxxz_0_yz_0[i] = g_xxxxx_0_y_1[i] * fi_acd_0 + g_xxxxx_0_yz_1[i] * wa_z[i];

        g_xxxxxz_0_zz_0[i] = 2.0 * g_xxxxx_0_z_1[i] * fi_acd_0 + g_xxxxx_0_zz_1[i] * wa_z[i];
    }

    /// Set up 18-24 components of targeted buffer : ISD

    auto g_xxxxyy_0_xx_0 = pbuffer.data(idx_eri_0_isd + 18);

    auto g_xxxxyy_0_xy_0 = pbuffer.data(idx_eri_0_isd + 19);

    auto g_xxxxyy_0_xz_0 = pbuffer.data(idx_eri_0_isd + 20);

    auto g_xxxxyy_0_yy_0 = pbuffer.data(idx_eri_0_isd + 21);

    auto g_xxxxyy_0_yz_0 = pbuffer.data(idx_eri_0_isd + 22);

    auto g_xxxxyy_0_zz_0 = pbuffer.data(idx_eri_0_isd + 23);

    #pragma omp simd aligned(g_xxxx_0_xx_0, g_xxxx_0_xx_1, g_xxxx_0_xz_0, g_xxxx_0_xz_1, g_xxxxy_0_xx_1, g_xxxxy_0_xz_1, g_xxxxyy_0_xx_0, g_xxxxyy_0_xy_0, g_xxxxyy_0_xz_0, g_xxxxyy_0_yy_0, g_xxxxyy_0_yz_0, g_xxxxyy_0_zz_0, g_xxxyy_0_xy_1, g_xxxyy_0_y_1, g_xxxyy_0_yy_1, g_xxxyy_0_yz_1, g_xxxyy_0_zz_1, g_xxyy_0_xy_0, g_xxyy_0_xy_1, g_xxyy_0_yy_0, g_xxyy_0_yy_1, g_xxyy_0_yz_0, g_xxyy_0_yz_1, g_xxyy_0_zz_0, g_xxyy_0_zz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxyy_0_xx_0[i] = g_xxxx_0_xx_0[i] * fbe_0 - g_xxxx_0_xx_1[i] * fz_be_0 + g_xxxxy_0_xx_1[i] * wa_y[i];

        g_xxxxyy_0_xy_0[i] = 3.0 * g_xxyy_0_xy_0[i] * fbe_0 - 3.0 * g_xxyy_0_xy_1[i] * fz_be_0 + g_xxxyy_0_y_1[i] * fi_acd_0 + g_xxxyy_0_xy_1[i] * wa_x[i];

        g_xxxxyy_0_xz_0[i] = g_xxxx_0_xz_0[i] * fbe_0 - g_xxxx_0_xz_1[i] * fz_be_0 + g_xxxxy_0_xz_1[i] * wa_y[i];

        g_xxxxyy_0_yy_0[i] = 3.0 * g_xxyy_0_yy_0[i] * fbe_0 - 3.0 * g_xxyy_0_yy_1[i] * fz_be_0 + g_xxxyy_0_yy_1[i] * wa_x[i];

        g_xxxxyy_0_yz_0[i] = 3.0 * g_xxyy_0_yz_0[i] * fbe_0 - 3.0 * g_xxyy_0_yz_1[i] * fz_be_0 + g_xxxyy_0_yz_1[i] * wa_x[i];

        g_xxxxyy_0_zz_0[i] = 3.0 * g_xxyy_0_zz_0[i] * fbe_0 - 3.0 * g_xxyy_0_zz_1[i] * fz_be_0 + g_xxxyy_0_zz_1[i] * wa_x[i];
    }

    /// Set up 24-30 components of targeted buffer : ISD

    auto g_xxxxyz_0_xx_0 = pbuffer.data(idx_eri_0_isd + 24);

    auto g_xxxxyz_0_xy_0 = pbuffer.data(idx_eri_0_isd + 25);

    auto g_xxxxyz_0_xz_0 = pbuffer.data(idx_eri_0_isd + 26);

    auto g_xxxxyz_0_yy_0 = pbuffer.data(idx_eri_0_isd + 27);

    auto g_xxxxyz_0_yz_0 = pbuffer.data(idx_eri_0_isd + 28);

    auto g_xxxxyz_0_zz_0 = pbuffer.data(idx_eri_0_isd + 29);

    #pragma omp simd aligned(g_xxxxy_0_xy_1, g_xxxxy_0_yy_1, g_xxxxyz_0_xx_0, g_xxxxyz_0_xy_0, g_xxxxyz_0_xz_0, g_xxxxyz_0_yy_0, g_xxxxyz_0_yz_0, g_xxxxyz_0_zz_0, g_xxxxz_0_xx_1, g_xxxxz_0_xz_1, g_xxxxz_0_yz_1, g_xxxxz_0_z_1, g_xxxxz_0_zz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxyz_0_xx_0[i] = g_xxxxz_0_xx_1[i] * wa_y[i];

        g_xxxxyz_0_xy_0[i] = g_xxxxy_0_xy_1[i] * wa_z[i];

        g_xxxxyz_0_xz_0[i] = g_xxxxz_0_xz_1[i] * wa_y[i];

        g_xxxxyz_0_yy_0[i] = g_xxxxy_0_yy_1[i] * wa_z[i];

        g_xxxxyz_0_yz_0[i] = g_xxxxz_0_z_1[i] * fi_acd_0 + g_xxxxz_0_yz_1[i] * wa_y[i];

        g_xxxxyz_0_zz_0[i] = g_xxxxz_0_zz_1[i] * wa_y[i];
    }

    /// Set up 30-36 components of targeted buffer : ISD

    auto g_xxxxzz_0_xx_0 = pbuffer.data(idx_eri_0_isd + 30);

    auto g_xxxxzz_0_xy_0 = pbuffer.data(idx_eri_0_isd + 31);

    auto g_xxxxzz_0_xz_0 = pbuffer.data(idx_eri_0_isd + 32);

    auto g_xxxxzz_0_yy_0 = pbuffer.data(idx_eri_0_isd + 33);

    auto g_xxxxzz_0_yz_0 = pbuffer.data(idx_eri_0_isd + 34);

    auto g_xxxxzz_0_zz_0 = pbuffer.data(idx_eri_0_isd + 35);

    #pragma omp simd aligned(g_xxxx_0_xx_0, g_xxxx_0_xx_1, g_xxxx_0_xy_0, g_xxxx_0_xy_1, g_xxxxz_0_xx_1, g_xxxxz_0_xy_1, g_xxxxzz_0_xx_0, g_xxxxzz_0_xy_0, g_xxxxzz_0_xz_0, g_xxxxzz_0_yy_0, g_xxxxzz_0_yz_0, g_xxxxzz_0_zz_0, g_xxxzz_0_xz_1, g_xxxzz_0_yy_1, g_xxxzz_0_yz_1, g_xxxzz_0_z_1, g_xxxzz_0_zz_1, g_xxzz_0_xz_0, g_xxzz_0_xz_1, g_xxzz_0_yy_0, g_xxzz_0_yy_1, g_xxzz_0_yz_0, g_xxzz_0_yz_1, g_xxzz_0_zz_0, g_xxzz_0_zz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxzz_0_xx_0[i] = g_xxxx_0_xx_0[i] * fbe_0 - g_xxxx_0_xx_1[i] * fz_be_0 + g_xxxxz_0_xx_1[i] * wa_z[i];

        g_xxxxzz_0_xy_0[i] = g_xxxx_0_xy_0[i] * fbe_0 - g_xxxx_0_xy_1[i] * fz_be_0 + g_xxxxz_0_xy_1[i] * wa_z[i];

        g_xxxxzz_0_xz_0[i] = 3.0 * g_xxzz_0_xz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xz_1[i] * fz_be_0 + g_xxxzz_0_z_1[i] * fi_acd_0 + g_xxxzz_0_xz_1[i] * wa_x[i];

        g_xxxxzz_0_yy_0[i] = 3.0 * g_xxzz_0_yy_0[i] * fbe_0 - 3.0 * g_xxzz_0_yy_1[i] * fz_be_0 + g_xxxzz_0_yy_1[i] * wa_x[i];

        g_xxxxzz_0_yz_0[i] = 3.0 * g_xxzz_0_yz_0[i] * fbe_0 - 3.0 * g_xxzz_0_yz_1[i] * fz_be_0 + g_xxxzz_0_yz_1[i] * wa_x[i];

        g_xxxxzz_0_zz_0[i] = 3.0 * g_xxzz_0_zz_0[i] * fbe_0 - 3.0 * g_xxzz_0_zz_1[i] * fz_be_0 + g_xxxzz_0_zz_1[i] * wa_x[i];
    }

    /// Set up 36-42 components of targeted buffer : ISD

    auto g_xxxyyy_0_xx_0 = pbuffer.data(idx_eri_0_isd + 36);

    auto g_xxxyyy_0_xy_0 = pbuffer.data(idx_eri_0_isd + 37);

    auto g_xxxyyy_0_xz_0 = pbuffer.data(idx_eri_0_isd + 38);

    auto g_xxxyyy_0_yy_0 = pbuffer.data(idx_eri_0_isd + 39);

    auto g_xxxyyy_0_yz_0 = pbuffer.data(idx_eri_0_isd + 40);

    auto g_xxxyyy_0_zz_0 = pbuffer.data(idx_eri_0_isd + 41);

    #pragma omp simd aligned(g_xxxy_0_xx_0, g_xxxy_0_xx_1, g_xxxy_0_xz_0, g_xxxy_0_xz_1, g_xxxyy_0_xx_1, g_xxxyy_0_xz_1, g_xxxyyy_0_xx_0, g_xxxyyy_0_xy_0, g_xxxyyy_0_xz_0, g_xxxyyy_0_yy_0, g_xxxyyy_0_yz_0, g_xxxyyy_0_zz_0, g_xxyyy_0_xy_1, g_xxyyy_0_y_1, g_xxyyy_0_yy_1, g_xxyyy_0_yz_1, g_xxyyy_0_zz_1, g_xyyy_0_xy_0, g_xyyy_0_xy_1, g_xyyy_0_yy_0, g_xyyy_0_yy_1, g_xyyy_0_yz_0, g_xyyy_0_yz_1, g_xyyy_0_zz_0, g_xyyy_0_zz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxyyy_0_xx_0[i] = 2.0 * g_xxxy_0_xx_0[i] * fbe_0 - 2.0 * g_xxxy_0_xx_1[i] * fz_be_0 + g_xxxyy_0_xx_1[i] * wa_y[i];

        g_xxxyyy_0_xy_0[i] = 2.0 * g_xyyy_0_xy_0[i] * fbe_0 - 2.0 * g_xyyy_0_xy_1[i] * fz_be_0 + g_xxyyy_0_y_1[i] * fi_acd_0 + g_xxyyy_0_xy_1[i] * wa_x[i];

        g_xxxyyy_0_xz_0[i] = 2.0 * g_xxxy_0_xz_0[i] * fbe_0 - 2.0 * g_xxxy_0_xz_1[i] * fz_be_0 + g_xxxyy_0_xz_1[i] * wa_y[i];

        g_xxxyyy_0_yy_0[i] = 2.0 * g_xyyy_0_yy_0[i] * fbe_0 - 2.0 * g_xyyy_0_yy_1[i] * fz_be_0 + g_xxyyy_0_yy_1[i] * wa_x[i];

        g_xxxyyy_0_yz_0[i] = 2.0 * g_xyyy_0_yz_0[i] * fbe_0 - 2.0 * g_xyyy_0_yz_1[i] * fz_be_0 + g_xxyyy_0_yz_1[i] * wa_x[i];

        g_xxxyyy_0_zz_0[i] = 2.0 * g_xyyy_0_zz_0[i] * fbe_0 - 2.0 * g_xyyy_0_zz_1[i] * fz_be_0 + g_xxyyy_0_zz_1[i] * wa_x[i];
    }

    /// Set up 42-48 components of targeted buffer : ISD

    auto g_xxxyyz_0_xx_0 = pbuffer.data(idx_eri_0_isd + 42);

    auto g_xxxyyz_0_xy_0 = pbuffer.data(idx_eri_0_isd + 43);

    auto g_xxxyyz_0_xz_0 = pbuffer.data(idx_eri_0_isd + 44);

    auto g_xxxyyz_0_yy_0 = pbuffer.data(idx_eri_0_isd + 45);

    auto g_xxxyyz_0_yz_0 = pbuffer.data(idx_eri_0_isd + 46);

    auto g_xxxyyz_0_zz_0 = pbuffer.data(idx_eri_0_isd + 47);

    #pragma omp simd aligned(g_xxxyy_0_x_1, g_xxxyy_0_xx_1, g_xxxyy_0_xy_1, g_xxxyy_0_xz_1, g_xxxyy_0_y_1, g_xxxyy_0_yy_1, g_xxxyy_0_yz_1, g_xxxyy_0_z_1, g_xxxyy_0_zz_1, g_xxxyyz_0_xx_0, g_xxxyyz_0_xy_0, g_xxxyyz_0_xz_0, g_xxxyyz_0_yy_0, g_xxxyyz_0_yz_0, g_xxxyyz_0_zz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyyz_0_xx_0[i] = g_xxxyy_0_xx_1[i] * wa_z[i];

        g_xxxyyz_0_xy_0[i] = g_xxxyy_0_xy_1[i] * wa_z[i];

        g_xxxyyz_0_xz_0[i] = g_xxxyy_0_x_1[i] * fi_acd_0 + g_xxxyy_0_xz_1[i] * wa_z[i];

        g_xxxyyz_0_yy_0[i] = g_xxxyy_0_yy_1[i] * wa_z[i];

        g_xxxyyz_0_yz_0[i] = g_xxxyy_0_y_1[i] * fi_acd_0 + g_xxxyy_0_yz_1[i] * wa_z[i];

        g_xxxyyz_0_zz_0[i] = 2.0 * g_xxxyy_0_z_1[i] * fi_acd_0 + g_xxxyy_0_zz_1[i] * wa_z[i];
    }

    /// Set up 48-54 components of targeted buffer : ISD

    auto g_xxxyzz_0_xx_0 = pbuffer.data(idx_eri_0_isd + 48);

    auto g_xxxyzz_0_xy_0 = pbuffer.data(idx_eri_0_isd + 49);

    auto g_xxxyzz_0_xz_0 = pbuffer.data(idx_eri_0_isd + 50);

    auto g_xxxyzz_0_yy_0 = pbuffer.data(idx_eri_0_isd + 51);

    auto g_xxxyzz_0_yz_0 = pbuffer.data(idx_eri_0_isd + 52);

    auto g_xxxyzz_0_zz_0 = pbuffer.data(idx_eri_0_isd + 53);

    #pragma omp simd aligned(g_xxxyzz_0_xx_0, g_xxxyzz_0_xy_0, g_xxxyzz_0_xz_0, g_xxxyzz_0_yy_0, g_xxxyzz_0_yz_0, g_xxxyzz_0_zz_0, g_xxxzz_0_x_1, g_xxxzz_0_xx_1, g_xxxzz_0_xy_1, g_xxxzz_0_xz_1, g_xxxzz_0_y_1, g_xxxzz_0_yy_1, g_xxxzz_0_yz_1, g_xxxzz_0_z_1, g_xxxzz_0_zz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyzz_0_xx_0[i] = g_xxxzz_0_xx_1[i] * wa_y[i];

        g_xxxyzz_0_xy_0[i] = g_xxxzz_0_x_1[i] * fi_acd_0 + g_xxxzz_0_xy_1[i] * wa_y[i];

        g_xxxyzz_0_xz_0[i] = g_xxxzz_0_xz_1[i] * wa_y[i];

        g_xxxyzz_0_yy_0[i] = 2.0 * g_xxxzz_0_y_1[i] * fi_acd_0 + g_xxxzz_0_yy_1[i] * wa_y[i];

        g_xxxyzz_0_yz_0[i] = g_xxxzz_0_z_1[i] * fi_acd_0 + g_xxxzz_0_yz_1[i] * wa_y[i];

        g_xxxyzz_0_zz_0[i] = g_xxxzz_0_zz_1[i] * wa_y[i];
    }

    /// Set up 54-60 components of targeted buffer : ISD

    auto g_xxxzzz_0_xx_0 = pbuffer.data(idx_eri_0_isd + 54);

    auto g_xxxzzz_0_xy_0 = pbuffer.data(idx_eri_0_isd + 55);

    auto g_xxxzzz_0_xz_0 = pbuffer.data(idx_eri_0_isd + 56);

    auto g_xxxzzz_0_yy_0 = pbuffer.data(idx_eri_0_isd + 57);

    auto g_xxxzzz_0_yz_0 = pbuffer.data(idx_eri_0_isd + 58);

    auto g_xxxzzz_0_zz_0 = pbuffer.data(idx_eri_0_isd + 59);

    #pragma omp simd aligned(g_xxxz_0_xx_0, g_xxxz_0_xx_1, g_xxxz_0_xy_0, g_xxxz_0_xy_1, g_xxxzz_0_xx_1, g_xxxzz_0_xy_1, g_xxxzzz_0_xx_0, g_xxxzzz_0_xy_0, g_xxxzzz_0_xz_0, g_xxxzzz_0_yy_0, g_xxxzzz_0_yz_0, g_xxxzzz_0_zz_0, g_xxzzz_0_xz_1, g_xxzzz_0_yy_1, g_xxzzz_0_yz_1, g_xxzzz_0_z_1, g_xxzzz_0_zz_1, g_xzzz_0_xz_0, g_xzzz_0_xz_1, g_xzzz_0_yy_0, g_xzzz_0_yy_1, g_xzzz_0_yz_0, g_xzzz_0_yz_1, g_xzzz_0_zz_0, g_xzzz_0_zz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxzzz_0_xx_0[i] = 2.0 * g_xxxz_0_xx_0[i] * fbe_0 - 2.0 * g_xxxz_0_xx_1[i] * fz_be_0 + g_xxxzz_0_xx_1[i] * wa_z[i];

        g_xxxzzz_0_xy_0[i] = 2.0 * g_xxxz_0_xy_0[i] * fbe_0 - 2.0 * g_xxxz_0_xy_1[i] * fz_be_0 + g_xxxzz_0_xy_1[i] * wa_z[i];

        g_xxxzzz_0_xz_0[i] = 2.0 * g_xzzz_0_xz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xz_1[i] * fz_be_0 + g_xxzzz_0_z_1[i] * fi_acd_0 + g_xxzzz_0_xz_1[i] * wa_x[i];

        g_xxxzzz_0_yy_0[i] = 2.0 * g_xzzz_0_yy_0[i] * fbe_0 - 2.0 * g_xzzz_0_yy_1[i] * fz_be_0 + g_xxzzz_0_yy_1[i] * wa_x[i];

        g_xxxzzz_0_yz_0[i] = 2.0 * g_xzzz_0_yz_0[i] * fbe_0 - 2.0 * g_xzzz_0_yz_1[i] * fz_be_0 + g_xxzzz_0_yz_1[i] * wa_x[i];

        g_xxxzzz_0_zz_0[i] = 2.0 * g_xzzz_0_zz_0[i] * fbe_0 - 2.0 * g_xzzz_0_zz_1[i] * fz_be_0 + g_xxzzz_0_zz_1[i] * wa_x[i];
    }

    /// Set up 60-66 components of targeted buffer : ISD

    auto g_xxyyyy_0_xx_0 = pbuffer.data(idx_eri_0_isd + 60);

    auto g_xxyyyy_0_xy_0 = pbuffer.data(idx_eri_0_isd + 61);

    auto g_xxyyyy_0_xz_0 = pbuffer.data(idx_eri_0_isd + 62);

    auto g_xxyyyy_0_yy_0 = pbuffer.data(idx_eri_0_isd + 63);

    auto g_xxyyyy_0_yz_0 = pbuffer.data(idx_eri_0_isd + 64);

    auto g_xxyyyy_0_zz_0 = pbuffer.data(idx_eri_0_isd + 65);

    #pragma omp simd aligned(g_xxyy_0_xx_0, g_xxyy_0_xx_1, g_xxyy_0_xz_0, g_xxyy_0_xz_1, g_xxyyy_0_xx_1, g_xxyyy_0_xz_1, g_xxyyyy_0_xx_0, g_xxyyyy_0_xy_0, g_xxyyyy_0_xz_0, g_xxyyyy_0_yy_0, g_xxyyyy_0_yz_0, g_xxyyyy_0_zz_0, g_xyyyy_0_xy_1, g_xyyyy_0_y_1, g_xyyyy_0_yy_1, g_xyyyy_0_yz_1, g_xyyyy_0_zz_1, g_yyyy_0_xy_0, g_yyyy_0_xy_1, g_yyyy_0_yy_0, g_yyyy_0_yy_1, g_yyyy_0_yz_0, g_yyyy_0_yz_1, g_yyyy_0_zz_0, g_yyyy_0_zz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyyyy_0_xx_0[i] = 3.0 * g_xxyy_0_xx_0[i] * fbe_0 - 3.0 * g_xxyy_0_xx_1[i] * fz_be_0 + g_xxyyy_0_xx_1[i] * wa_y[i];

        g_xxyyyy_0_xy_0[i] = g_yyyy_0_xy_0[i] * fbe_0 - g_yyyy_0_xy_1[i] * fz_be_0 + g_xyyyy_0_y_1[i] * fi_acd_0 + g_xyyyy_0_xy_1[i] * wa_x[i];

        g_xxyyyy_0_xz_0[i] = 3.0 * g_xxyy_0_xz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xz_1[i] * fz_be_0 + g_xxyyy_0_xz_1[i] * wa_y[i];

        g_xxyyyy_0_yy_0[i] = g_yyyy_0_yy_0[i] * fbe_0 - g_yyyy_0_yy_1[i] * fz_be_0 + g_xyyyy_0_yy_1[i] * wa_x[i];

        g_xxyyyy_0_yz_0[i] = g_yyyy_0_yz_0[i] * fbe_0 - g_yyyy_0_yz_1[i] * fz_be_0 + g_xyyyy_0_yz_1[i] * wa_x[i];

        g_xxyyyy_0_zz_0[i] = g_yyyy_0_zz_0[i] * fbe_0 - g_yyyy_0_zz_1[i] * fz_be_0 + g_xyyyy_0_zz_1[i] * wa_x[i];
    }

    /// Set up 66-72 components of targeted buffer : ISD

    auto g_xxyyyz_0_xx_0 = pbuffer.data(idx_eri_0_isd + 66);

    auto g_xxyyyz_0_xy_0 = pbuffer.data(idx_eri_0_isd + 67);

    auto g_xxyyyz_0_xz_0 = pbuffer.data(idx_eri_0_isd + 68);

    auto g_xxyyyz_0_yy_0 = pbuffer.data(idx_eri_0_isd + 69);

    auto g_xxyyyz_0_yz_0 = pbuffer.data(idx_eri_0_isd + 70);

    auto g_xxyyyz_0_zz_0 = pbuffer.data(idx_eri_0_isd + 71);

    #pragma omp simd aligned(g_xxyyy_0_x_1, g_xxyyy_0_xx_1, g_xxyyy_0_xy_1, g_xxyyy_0_xz_1, g_xxyyy_0_y_1, g_xxyyy_0_yy_1, g_xxyyy_0_yz_1, g_xxyyy_0_z_1, g_xxyyy_0_zz_1, g_xxyyyz_0_xx_0, g_xxyyyz_0_xy_0, g_xxyyyz_0_xz_0, g_xxyyyz_0_yy_0, g_xxyyyz_0_yz_0, g_xxyyyz_0_zz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyyyz_0_xx_0[i] = g_xxyyy_0_xx_1[i] * wa_z[i];

        g_xxyyyz_0_xy_0[i] = g_xxyyy_0_xy_1[i] * wa_z[i];

        g_xxyyyz_0_xz_0[i] = g_xxyyy_0_x_1[i] * fi_acd_0 + g_xxyyy_0_xz_1[i] * wa_z[i];

        g_xxyyyz_0_yy_0[i] = g_xxyyy_0_yy_1[i] * wa_z[i];

        g_xxyyyz_0_yz_0[i] = g_xxyyy_0_y_1[i] * fi_acd_0 + g_xxyyy_0_yz_1[i] * wa_z[i];

        g_xxyyyz_0_zz_0[i] = 2.0 * g_xxyyy_0_z_1[i] * fi_acd_0 + g_xxyyy_0_zz_1[i] * wa_z[i];
    }

    /// Set up 72-78 components of targeted buffer : ISD

    auto g_xxyyzz_0_xx_0 = pbuffer.data(idx_eri_0_isd + 72);

    auto g_xxyyzz_0_xy_0 = pbuffer.data(idx_eri_0_isd + 73);

    auto g_xxyyzz_0_xz_0 = pbuffer.data(idx_eri_0_isd + 74);

    auto g_xxyyzz_0_yy_0 = pbuffer.data(idx_eri_0_isd + 75);

    auto g_xxyyzz_0_yz_0 = pbuffer.data(idx_eri_0_isd + 76);

    auto g_xxyyzz_0_zz_0 = pbuffer.data(idx_eri_0_isd + 77);

    #pragma omp simd aligned(g_xxyy_0_xy_0, g_xxyy_0_xy_1, g_xxyyz_0_xy_1, g_xxyyzz_0_xx_0, g_xxyyzz_0_xy_0, g_xxyyzz_0_xz_0, g_xxyyzz_0_yy_0, g_xxyyzz_0_yz_0, g_xxyyzz_0_zz_0, g_xxyzz_0_xx_1, g_xxyzz_0_xz_1, g_xxzz_0_xx_0, g_xxzz_0_xx_1, g_xxzz_0_xz_0, g_xxzz_0_xz_1, g_xyyzz_0_yy_1, g_xyyzz_0_yz_1, g_xyyzz_0_zz_1, g_yyzz_0_yy_0, g_yyzz_0_yy_1, g_yyzz_0_yz_0, g_yyzz_0_yz_1, g_yyzz_0_zz_0, g_yyzz_0_zz_1, wa_x, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyyzz_0_xx_0[i] = g_xxzz_0_xx_0[i] * fbe_0 - g_xxzz_0_xx_1[i] * fz_be_0 + g_xxyzz_0_xx_1[i] * wa_y[i];

        g_xxyyzz_0_xy_0[i] = g_xxyy_0_xy_0[i] * fbe_0 - g_xxyy_0_xy_1[i] * fz_be_0 + g_xxyyz_0_xy_1[i] * wa_z[i];

        g_xxyyzz_0_xz_0[i] = g_xxzz_0_xz_0[i] * fbe_0 - g_xxzz_0_xz_1[i] * fz_be_0 + g_xxyzz_0_xz_1[i] * wa_y[i];

        g_xxyyzz_0_yy_0[i] = g_yyzz_0_yy_0[i] * fbe_0 - g_yyzz_0_yy_1[i] * fz_be_0 + g_xyyzz_0_yy_1[i] * wa_x[i];

        g_xxyyzz_0_yz_0[i] = g_yyzz_0_yz_0[i] * fbe_0 - g_yyzz_0_yz_1[i] * fz_be_0 + g_xyyzz_0_yz_1[i] * wa_x[i];

        g_xxyyzz_0_zz_0[i] = g_yyzz_0_zz_0[i] * fbe_0 - g_yyzz_0_zz_1[i] * fz_be_0 + g_xyyzz_0_zz_1[i] * wa_x[i];
    }

    /// Set up 78-84 components of targeted buffer : ISD

    auto g_xxyzzz_0_xx_0 = pbuffer.data(idx_eri_0_isd + 78);

    auto g_xxyzzz_0_xy_0 = pbuffer.data(idx_eri_0_isd + 79);

    auto g_xxyzzz_0_xz_0 = pbuffer.data(idx_eri_0_isd + 80);

    auto g_xxyzzz_0_yy_0 = pbuffer.data(idx_eri_0_isd + 81);

    auto g_xxyzzz_0_yz_0 = pbuffer.data(idx_eri_0_isd + 82);

    auto g_xxyzzz_0_zz_0 = pbuffer.data(idx_eri_0_isd + 83);

    #pragma omp simd aligned(g_xxyzzz_0_xx_0, g_xxyzzz_0_xy_0, g_xxyzzz_0_xz_0, g_xxyzzz_0_yy_0, g_xxyzzz_0_yz_0, g_xxyzzz_0_zz_0, g_xxzzz_0_x_1, g_xxzzz_0_xx_1, g_xxzzz_0_xy_1, g_xxzzz_0_xz_1, g_xxzzz_0_y_1, g_xxzzz_0_yy_1, g_xxzzz_0_yz_1, g_xxzzz_0_z_1, g_xxzzz_0_zz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyzzz_0_xx_0[i] = g_xxzzz_0_xx_1[i] * wa_y[i];

        g_xxyzzz_0_xy_0[i] = g_xxzzz_0_x_1[i] * fi_acd_0 + g_xxzzz_0_xy_1[i] * wa_y[i];

        g_xxyzzz_0_xz_0[i] = g_xxzzz_0_xz_1[i] * wa_y[i];

        g_xxyzzz_0_yy_0[i] = 2.0 * g_xxzzz_0_y_1[i] * fi_acd_0 + g_xxzzz_0_yy_1[i] * wa_y[i];

        g_xxyzzz_0_yz_0[i] = g_xxzzz_0_z_1[i] * fi_acd_0 + g_xxzzz_0_yz_1[i] * wa_y[i];

        g_xxyzzz_0_zz_0[i] = g_xxzzz_0_zz_1[i] * wa_y[i];
    }

    /// Set up 84-90 components of targeted buffer : ISD

    auto g_xxzzzz_0_xx_0 = pbuffer.data(idx_eri_0_isd + 84);

    auto g_xxzzzz_0_xy_0 = pbuffer.data(idx_eri_0_isd + 85);

    auto g_xxzzzz_0_xz_0 = pbuffer.data(idx_eri_0_isd + 86);

    auto g_xxzzzz_0_yy_0 = pbuffer.data(idx_eri_0_isd + 87);

    auto g_xxzzzz_0_yz_0 = pbuffer.data(idx_eri_0_isd + 88);

    auto g_xxzzzz_0_zz_0 = pbuffer.data(idx_eri_0_isd + 89);

    #pragma omp simd aligned(g_xxzz_0_xx_0, g_xxzz_0_xx_1, g_xxzz_0_xy_0, g_xxzz_0_xy_1, g_xxzzz_0_xx_1, g_xxzzz_0_xy_1, g_xxzzzz_0_xx_0, g_xxzzzz_0_xy_0, g_xxzzzz_0_xz_0, g_xxzzzz_0_yy_0, g_xxzzzz_0_yz_0, g_xxzzzz_0_zz_0, g_xzzzz_0_xz_1, g_xzzzz_0_yy_1, g_xzzzz_0_yz_1, g_xzzzz_0_z_1, g_xzzzz_0_zz_1, g_zzzz_0_xz_0, g_zzzz_0_xz_1, g_zzzz_0_yy_0, g_zzzz_0_yy_1, g_zzzz_0_yz_0, g_zzzz_0_yz_1, g_zzzz_0_zz_0, g_zzzz_0_zz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxzzzz_0_xx_0[i] = 3.0 * g_xxzz_0_xx_0[i] * fbe_0 - 3.0 * g_xxzz_0_xx_1[i] * fz_be_0 + g_xxzzz_0_xx_1[i] * wa_z[i];

        g_xxzzzz_0_xy_0[i] = 3.0 * g_xxzz_0_xy_0[i] * fbe_0 - 3.0 * g_xxzz_0_xy_1[i] * fz_be_0 + g_xxzzz_0_xy_1[i] * wa_z[i];

        g_xxzzzz_0_xz_0[i] = g_zzzz_0_xz_0[i] * fbe_0 - g_zzzz_0_xz_1[i] * fz_be_0 + g_xzzzz_0_z_1[i] * fi_acd_0 + g_xzzzz_0_xz_1[i] * wa_x[i];

        g_xxzzzz_0_yy_0[i] = g_zzzz_0_yy_0[i] * fbe_0 - g_zzzz_0_yy_1[i] * fz_be_0 + g_xzzzz_0_yy_1[i] * wa_x[i];

        g_xxzzzz_0_yz_0[i] = g_zzzz_0_yz_0[i] * fbe_0 - g_zzzz_0_yz_1[i] * fz_be_0 + g_xzzzz_0_yz_1[i] * wa_x[i];

        g_xxzzzz_0_zz_0[i] = g_zzzz_0_zz_0[i] * fbe_0 - g_zzzz_0_zz_1[i] * fz_be_0 + g_xzzzz_0_zz_1[i] * wa_x[i];
    }

    /// Set up 90-96 components of targeted buffer : ISD

    auto g_xyyyyy_0_xx_0 = pbuffer.data(idx_eri_0_isd + 90);

    auto g_xyyyyy_0_xy_0 = pbuffer.data(idx_eri_0_isd + 91);

    auto g_xyyyyy_0_xz_0 = pbuffer.data(idx_eri_0_isd + 92);

    auto g_xyyyyy_0_yy_0 = pbuffer.data(idx_eri_0_isd + 93);

    auto g_xyyyyy_0_yz_0 = pbuffer.data(idx_eri_0_isd + 94);

    auto g_xyyyyy_0_zz_0 = pbuffer.data(idx_eri_0_isd + 95);

    #pragma omp simd aligned(g_xyyyyy_0_xx_0, g_xyyyyy_0_xy_0, g_xyyyyy_0_xz_0, g_xyyyyy_0_yy_0, g_xyyyyy_0_yz_0, g_xyyyyy_0_zz_0, g_yyyyy_0_x_1, g_yyyyy_0_xx_1, g_yyyyy_0_xy_1, g_yyyyy_0_xz_1, g_yyyyy_0_y_1, g_yyyyy_0_yy_1, g_yyyyy_0_yz_1, g_yyyyy_0_z_1, g_yyyyy_0_zz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyyy_0_xx_0[i] = 2.0 * g_yyyyy_0_x_1[i] * fi_acd_0 + g_yyyyy_0_xx_1[i] * wa_x[i];

        g_xyyyyy_0_xy_0[i] = g_yyyyy_0_y_1[i] * fi_acd_0 + g_yyyyy_0_xy_1[i] * wa_x[i];

        g_xyyyyy_0_xz_0[i] = g_yyyyy_0_z_1[i] * fi_acd_0 + g_yyyyy_0_xz_1[i] * wa_x[i];

        g_xyyyyy_0_yy_0[i] = g_yyyyy_0_yy_1[i] * wa_x[i];

        g_xyyyyy_0_yz_0[i] = g_yyyyy_0_yz_1[i] * wa_x[i];

        g_xyyyyy_0_zz_0[i] = g_yyyyy_0_zz_1[i] * wa_x[i];
    }

    /// Set up 96-102 components of targeted buffer : ISD

    auto g_xyyyyz_0_xx_0 = pbuffer.data(idx_eri_0_isd + 96);

    auto g_xyyyyz_0_xy_0 = pbuffer.data(idx_eri_0_isd + 97);

    auto g_xyyyyz_0_xz_0 = pbuffer.data(idx_eri_0_isd + 98);

    auto g_xyyyyz_0_yy_0 = pbuffer.data(idx_eri_0_isd + 99);

    auto g_xyyyyz_0_yz_0 = pbuffer.data(idx_eri_0_isd + 100);

    auto g_xyyyyz_0_zz_0 = pbuffer.data(idx_eri_0_isd + 101);

    #pragma omp simd aligned(g_xyyyy_0_xx_1, g_xyyyy_0_xy_1, g_xyyyyz_0_xx_0, g_xyyyyz_0_xy_0, g_xyyyyz_0_xz_0, g_xyyyyz_0_yy_0, g_xyyyyz_0_yz_0, g_xyyyyz_0_zz_0, g_yyyyz_0_xz_1, g_yyyyz_0_yy_1, g_yyyyz_0_yz_1, g_yyyyz_0_z_1, g_yyyyz_0_zz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyyz_0_xx_0[i] = g_xyyyy_0_xx_1[i] * wa_z[i];

        g_xyyyyz_0_xy_0[i] = g_xyyyy_0_xy_1[i] * wa_z[i];

        g_xyyyyz_0_xz_0[i] = g_yyyyz_0_z_1[i] * fi_acd_0 + g_yyyyz_0_xz_1[i] * wa_x[i];

        g_xyyyyz_0_yy_0[i] = g_yyyyz_0_yy_1[i] * wa_x[i];

        g_xyyyyz_0_yz_0[i] = g_yyyyz_0_yz_1[i] * wa_x[i];

        g_xyyyyz_0_zz_0[i] = g_yyyyz_0_zz_1[i] * wa_x[i];
    }

    /// Set up 102-108 components of targeted buffer : ISD

    auto g_xyyyzz_0_xx_0 = pbuffer.data(idx_eri_0_isd + 102);

    auto g_xyyyzz_0_xy_0 = pbuffer.data(idx_eri_0_isd + 103);

    auto g_xyyyzz_0_xz_0 = pbuffer.data(idx_eri_0_isd + 104);

    auto g_xyyyzz_0_yy_0 = pbuffer.data(idx_eri_0_isd + 105);

    auto g_xyyyzz_0_yz_0 = pbuffer.data(idx_eri_0_isd + 106);

    auto g_xyyyzz_0_zz_0 = pbuffer.data(idx_eri_0_isd + 107);

    #pragma omp simd aligned(g_xyyyzz_0_xx_0, g_xyyyzz_0_xy_0, g_xyyyzz_0_xz_0, g_xyyyzz_0_yy_0, g_xyyyzz_0_yz_0, g_xyyyzz_0_zz_0, g_yyyzz_0_x_1, g_yyyzz_0_xx_1, g_yyyzz_0_xy_1, g_yyyzz_0_xz_1, g_yyyzz_0_y_1, g_yyyzz_0_yy_1, g_yyyzz_0_yz_1, g_yyyzz_0_z_1, g_yyyzz_0_zz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyzz_0_xx_0[i] = 2.0 * g_yyyzz_0_x_1[i] * fi_acd_0 + g_yyyzz_0_xx_1[i] * wa_x[i];

        g_xyyyzz_0_xy_0[i] = g_yyyzz_0_y_1[i] * fi_acd_0 + g_yyyzz_0_xy_1[i] * wa_x[i];

        g_xyyyzz_0_xz_0[i] = g_yyyzz_0_z_1[i] * fi_acd_0 + g_yyyzz_0_xz_1[i] * wa_x[i];

        g_xyyyzz_0_yy_0[i] = g_yyyzz_0_yy_1[i] * wa_x[i];

        g_xyyyzz_0_yz_0[i] = g_yyyzz_0_yz_1[i] * wa_x[i];

        g_xyyyzz_0_zz_0[i] = g_yyyzz_0_zz_1[i] * wa_x[i];
    }

    /// Set up 108-114 components of targeted buffer : ISD

    auto g_xyyzzz_0_xx_0 = pbuffer.data(idx_eri_0_isd + 108);

    auto g_xyyzzz_0_xy_0 = pbuffer.data(idx_eri_0_isd + 109);

    auto g_xyyzzz_0_xz_0 = pbuffer.data(idx_eri_0_isd + 110);

    auto g_xyyzzz_0_yy_0 = pbuffer.data(idx_eri_0_isd + 111);

    auto g_xyyzzz_0_yz_0 = pbuffer.data(idx_eri_0_isd + 112);

    auto g_xyyzzz_0_zz_0 = pbuffer.data(idx_eri_0_isd + 113);

    #pragma omp simd aligned(g_xyyzzz_0_xx_0, g_xyyzzz_0_xy_0, g_xyyzzz_0_xz_0, g_xyyzzz_0_yy_0, g_xyyzzz_0_yz_0, g_xyyzzz_0_zz_0, g_yyzzz_0_x_1, g_yyzzz_0_xx_1, g_yyzzz_0_xy_1, g_yyzzz_0_xz_1, g_yyzzz_0_y_1, g_yyzzz_0_yy_1, g_yyzzz_0_yz_1, g_yyzzz_0_z_1, g_yyzzz_0_zz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyzzz_0_xx_0[i] = 2.0 * g_yyzzz_0_x_1[i] * fi_acd_0 + g_yyzzz_0_xx_1[i] * wa_x[i];

        g_xyyzzz_0_xy_0[i] = g_yyzzz_0_y_1[i] * fi_acd_0 + g_yyzzz_0_xy_1[i] * wa_x[i];

        g_xyyzzz_0_xz_0[i] = g_yyzzz_0_z_1[i] * fi_acd_0 + g_yyzzz_0_xz_1[i] * wa_x[i];

        g_xyyzzz_0_yy_0[i] = g_yyzzz_0_yy_1[i] * wa_x[i];

        g_xyyzzz_0_yz_0[i] = g_yyzzz_0_yz_1[i] * wa_x[i];

        g_xyyzzz_0_zz_0[i] = g_yyzzz_0_zz_1[i] * wa_x[i];
    }

    /// Set up 114-120 components of targeted buffer : ISD

    auto g_xyzzzz_0_xx_0 = pbuffer.data(idx_eri_0_isd + 114);

    auto g_xyzzzz_0_xy_0 = pbuffer.data(idx_eri_0_isd + 115);

    auto g_xyzzzz_0_xz_0 = pbuffer.data(idx_eri_0_isd + 116);

    auto g_xyzzzz_0_yy_0 = pbuffer.data(idx_eri_0_isd + 117);

    auto g_xyzzzz_0_yz_0 = pbuffer.data(idx_eri_0_isd + 118);

    auto g_xyzzzz_0_zz_0 = pbuffer.data(idx_eri_0_isd + 119);

    #pragma omp simd aligned(g_xyzzzz_0_xx_0, g_xyzzzz_0_xy_0, g_xyzzzz_0_xz_0, g_xyzzzz_0_yy_0, g_xyzzzz_0_yz_0, g_xyzzzz_0_zz_0, g_xzzzz_0_xx_1, g_xzzzz_0_xz_1, g_yzzzz_0_xy_1, g_yzzzz_0_y_1, g_yzzzz_0_yy_1, g_yzzzz_0_yz_1, g_yzzzz_0_zz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyzzzz_0_xx_0[i] = g_xzzzz_0_xx_1[i] * wa_y[i];

        g_xyzzzz_0_xy_0[i] = g_yzzzz_0_y_1[i] * fi_acd_0 + g_yzzzz_0_xy_1[i] * wa_x[i];

        g_xyzzzz_0_xz_0[i] = g_xzzzz_0_xz_1[i] * wa_y[i];

        g_xyzzzz_0_yy_0[i] = g_yzzzz_0_yy_1[i] * wa_x[i];

        g_xyzzzz_0_yz_0[i] = g_yzzzz_0_yz_1[i] * wa_x[i];

        g_xyzzzz_0_zz_0[i] = g_yzzzz_0_zz_1[i] * wa_x[i];
    }

    /// Set up 120-126 components of targeted buffer : ISD

    auto g_xzzzzz_0_xx_0 = pbuffer.data(idx_eri_0_isd + 120);

    auto g_xzzzzz_0_xy_0 = pbuffer.data(idx_eri_0_isd + 121);

    auto g_xzzzzz_0_xz_0 = pbuffer.data(idx_eri_0_isd + 122);

    auto g_xzzzzz_0_yy_0 = pbuffer.data(idx_eri_0_isd + 123);

    auto g_xzzzzz_0_yz_0 = pbuffer.data(idx_eri_0_isd + 124);

    auto g_xzzzzz_0_zz_0 = pbuffer.data(idx_eri_0_isd + 125);

    #pragma omp simd aligned(g_xzzzzz_0_xx_0, g_xzzzzz_0_xy_0, g_xzzzzz_0_xz_0, g_xzzzzz_0_yy_0, g_xzzzzz_0_yz_0, g_xzzzzz_0_zz_0, g_zzzzz_0_x_1, g_zzzzz_0_xx_1, g_zzzzz_0_xy_1, g_zzzzz_0_xz_1, g_zzzzz_0_y_1, g_zzzzz_0_yy_1, g_zzzzz_0_yz_1, g_zzzzz_0_z_1, g_zzzzz_0_zz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzzzzz_0_xx_0[i] = 2.0 * g_zzzzz_0_x_1[i] * fi_acd_0 + g_zzzzz_0_xx_1[i] * wa_x[i];

        g_xzzzzz_0_xy_0[i] = g_zzzzz_0_y_1[i] * fi_acd_0 + g_zzzzz_0_xy_1[i] * wa_x[i];

        g_xzzzzz_0_xz_0[i] = g_zzzzz_0_z_1[i] * fi_acd_0 + g_zzzzz_0_xz_1[i] * wa_x[i];

        g_xzzzzz_0_yy_0[i] = g_zzzzz_0_yy_1[i] * wa_x[i];

        g_xzzzzz_0_yz_0[i] = g_zzzzz_0_yz_1[i] * wa_x[i];

        g_xzzzzz_0_zz_0[i] = g_zzzzz_0_zz_1[i] * wa_x[i];
    }

    /// Set up 126-132 components of targeted buffer : ISD

    auto g_yyyyyy_0_xx_0 = pbuffer.data(idx_eri_0_isd + 126);

    auto g_yyyyyy_0_xy_0 = pbuffer.data(idx_eri_0_isd + 127);

    auto g_yyyyyy_0_xz_0 = pbuffer.data(idx_eri_0_isd + 128);

    auto g_yyyyyy_0_yy_0 = pbuffer.data(idx_eri_0_isd + 129);

    auto g_yyyyyy_0_yz_0 = pbuffer.data(idx_eri_0_isd + 130);

    auto g_yyyyyy_0_zz_0 = pbuffer.data(idx_eri_0_isd + 131);

    #pragma omp simd aligned(g_yyyy_0_xx_0, g_yyyy_0_xx_1, g_yyyy_0_xy_0, g_yyyy_0_xy_1, g_yyyy_0_xz_0, g_yyyy_0_xz_1, g_yyyy_0_yy_0, g_yyyy_0_yy_1, g_yyyy_0_yz_0, g_yyyy_0_yz_1, g_yyyy_0_zz_0, g_yyyy_0_zz_1, g_yyyyy_0_x_1, g_yyyyy_0_xx_1, g_yyyyy_0_xy_1, g_yyyyy_0_xz_1, g_yyyyy_0_y_1, g_yyyyy_0_yy_1, g_yyyyy_0_yz_1, g_yyyyy_0_z_1, g_yyyyy_0_zz_1, g_yyyyyy_0_xx_0, g_yyyyyy_0_xy_0, g_yyyyyy_0_xz_0, g_yyyyyy_0_yy_0, g_yyyyyy_0_yz_0, g_yyyyyy_0_zz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyyy_0_xx_0[i] = 5.0 * g_yyyy_0_xx_0[i] * fbe_0 - 5.0 * g_yyyy_0_xx_1[i] * fz_be_0 + g_yyyyy_0_xx_1[i] * wa_y[i];

        g_yyyyyy_0_xy_0[i] = 5.0 * g_yyyy_0_xy_0[i] * fbe_0 - 5.0 * g_yyyy_0_xy_1[i] * fz_be_0 + g_yyyyy_0_x_1[i] * fi_acd_0 + g_yyyyy_0_xy_1[i] * wa_y[i];

        g_yyyyyy_0_xz_0[i] = 5.0 * g_yyyy_0_xz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xz_1[i] * fz_be_0 + g_yyyyy_0_xz_1[i] * wa_y[i];

        g_yyyyyy_0_yy_0[i] = 5.0 * g_yyyy_0_yy_0[i] * fbe_0 - 5.0 * g_yyyy_0_yy_1[i] * fz_be_0 + 2.0 * g_yyyyy_0_y_1[i] * fi_acd_0 + g_yyyyy_0_yy_1[i] * wa_y[i];

        g_yyyyyy_0_yz_0[i] = 5.0 * g_yyyy_0_yz_0[i] * fbe_0 - 5.0 * g_yyyy_0_yz_1[i] * fz_be_0 + g_yyyyy_0_z_1[i] * fi_acd_0 + g_yyyyy_0_yz_1[i] * wa_y[i];

        g_yyyyyy_0_zz_0[i] = 5.0 * g_yyyy_0_zz_0[i] * fbe_0 - 5.0 * g_yyyy_0_zz_1[i] * fz_be_0 + g_yyyyy_0_zz_1[i] * wa_y[i];
    }

    /// Set up 132-138 components of targeted buffer : ISD

    auto g_yyyyyz_0_xx_0 = pbuffer.data(idx_eri_0_isd + 132);

    auto g_yyyyyz_0_xy_0 = pbuffer.data(idx_eri_0_isd + 133);

    auto g_yyyyyz_0_xz_0 = pbuffer.data(idx_eri_0_isd + 134);

    auto g_yyyyyz_0_yy_0 = pbuffer.data(idx_eri_0_isd + 135);

    auto g_yyyyyz_0_yz_0 = pbuffer.data(idx_eri_0_isd + 136);

    auto g_yyyyyz_0_zz_0 = pbuffer.data(idx_eri_0_isd + 137);

    #pragma omp simd aligned(g_yyyyy_0_x_1, g_yyyyy_0_xx_1, g_yyyyy_0_xy_1, g_yyyyy_0_xz_1, g_yyyyy_0_y_1, g_yyyyy_0_yy_1, g_yyyyy_0_yz_1, g_yyyyy_0_z_1, g_yyyyy_0_zz_1, g_yyyyyz_0_xx_0, g_yyyyyz_0_xy_0, g_yyyyyz_0_xz_0, g_yyyyyz_0_yy_0, g_yyyyyz_0_yz_0, g_yyyyyz_0_zz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyyyyz_0_xx_0[i] = g_yyyyy_0_xx_1[i] * wa_z[i];

        g_yyyyyz_0_xy_0[i] = g_yyyyy_0_xy_1[i] * wa_z[i];

        g_yyyyyz_0_xz_0[i] = g_yyyyy_0_x_1[i] * fi_acd_0 + g_yyyyy_0_xz_1[i] * wa_z[i];

        g_yyyyyz_0_yy_0[i] = g_yyyyy_0_yy_1[i] * wa_z[i];

        g_yyyyyz_0_yz_0[i] = g_yyyyy_0_y_1[i] * fi_acd_0 + g_yyyyy_0_yz_1[i] * wa_z[i];

        g_yyyyyz_0_zz_0[i] = 2.0 * g_yyyyy_0_z_1[i] * fi_acd_0 + g_yyyyy_0_zz_1[i] * wa_z[i];
    }

    /// Set up 138-144 components of targeted buffer : ISD

    auto g_yyyyzz_0_xx_0 = pbuffer.data(idx_eri_0_isd + 138);

    auto g_yyyyzz_0_xy_0 = pbuffer.data(idx_eri_0_isd + 139);

    auto g_yyyyzz_0_xz_0 = pbuffer.data(idx_eri_0_isd + 140);

    auto g_yyyyzz_0_yy_0 = pbuffer.data(idx_eri_0_isd + 141);

    auto g_yyyyzz_0_yz_0 = pbuffer.data(idx_eri_0_isd + 142);

    auto g_yyyyzz_0_zz_0 = pbuffer.data(idx_eri_0_isd + 143);

    #pragma omp simd aligned(g_yyyy_0_xy_0, g_yyyy_0_xy_1, g_yyyy_0_yy_0, g_yyyy_0_yy_1, g_yyyyz_0_xy_1, g_yyyyz_0_yy_1, g_yyyyzz_0_xx_0, g_yyyyzz_0_xy_0, g_yyyyzz_0_xz_0, g_yyyyzz_0_yy_0, g_yyyyzz_0_yz_0, g_yyyyzz_0_zz_0, g_yyyzz_0_xx_1, g_yyyzz_0_xz_1, g_yyyzz_0_yz_1, g_yyyzz_0_z_1, g_yyyzz_0_zz_1, g_yyzz_0_xx_0, g_yyzz_0_xx_1, g_yyzz_0_xz_0, g_yyzz_0_xz_1, g_yyzz_0_yz_0, g_yyzz_0_yz_1, g_yyzz_0_zz_0, g_yyzz_0_zz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyzz_0_xx_0[i] = 3.0 * g_yyzz_0_xx_0[i] * fbe_0 - 3.0 * g_yyzz_0_xx_1[i] * fz_be_0 + g_yyyzz_0_xx_1[i] * wa_y[i];

        g_yyyyzz_0_xy_0[i] = g_yyyy_0_xy_0[i] * fbe_0 - g_yyyy_0_xy_1[i] * fz_be_0 + g_yyyyz_0_xy_1[i] * wa_z[i];

        g_yyyyzz_0_xz_0[i] = 3.0 * g_yyzz_0_xz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xz_1[i] * fz_be_0 + g_yyyzz_0_xz_1[i] * wa_y[i];

        g_yyyyzz_0_yy_0[i] = g_yyyy_0_yy_0[i] * fbe_0 - g_yyyy_0_yy_1[i] * fz_be_0 + g_yyyyz_0_yy_1[i] * wa_z[i];

        g_yyyyzz_0_yz_0[i] = 3.0 * g_yyzz_0_yz_0[i] * fbe_0 - 3.0 * g_yyzz_0_yz_1[i] * fz_be_0 + g_yyyzz_0_z_1[i] * fi_acd_0 + g_yyyzz_0_yz_1[i] * wa_y[i];

        g_yyyyzz_0_zz_0[i] = 3.0 * g_yyzz_0_zz_0[i] * fbe_0 - 3.0 * g_yyzz_0_zz_1[i] * fz_be_0 + g_yyyzz_0_zz_1[i] * wa_y[i];
    }

    /// Set up 144-150 components of targeted buffer : ISD

    auto g_yyyzzz_0_xx_0 = pbuffer.data(idx_eri_0_isd + 144);

    auto g_yyyzzz_0_xy_0 = pbuffer.data(idx_eri_0_isd + 145);

    auto g_yyyzzz_0_xz_0 = pbuffer.data(idx_eri_0_isd + 146);

    auto g_yyyzzz_0_yy_0 = pbuffer.data(idx_eri_0_isd + 147);

    auto g_yyyzzz_0_yz_0 = pbuffer.data(idx_eri_0_isd + 148);

    auto g_yyyzzz_0_zz_0 = pbuffer.data(idx_eri_0_isd + 149);

    #pragma omp simd aligned(g_yyyz_0_xy_0, g_yyyz_0_xy_1, g_yyyz_0_yy_0, g_yyyz_0_yy_1, g_yyyzz_0_xy_1, g_yyyzz_0_yy_1, g_yyyzzz_0_xx_0, g_yyyzzz_0_xy_0, g_yyyzzz_0_xz_0, g_yyyzzz_0_yy_0, g_yyyzzz_0_yz_0, g_yyyzzz_0_zz_0, g_yyzzz_0_xx_1, g_yyzzz_0_xz_1, g_yyzzz_0_yz_1, g_yyzzz_0_z_1, g_yyzzz_0_zz_1, g_yzzz_0_xx_0, g_yzzz_0_xx_1, g_yzzz_0_xz_0, g_yzzz_0_xz_1, g_yzzz_0_yz_0, g_yzzz_0_yz_1, g_yzzz_0_zz_0, g_yzzz_0_zz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyzzz_0_xx_0[i] = 2.0 * g_yzzz_0_xx_0[i] * fbe_0 - 2.0 * g_yzzz_0_xx_1[i] * fz_be_0 + g_yyzzz_0_xx_1[i] * wa_y[i];

        g_yyyzzz_0_xy_0[i] = 2.0 * g_yyyz_0_xy_0[i] * fbe_0 - 2.0 * g_yyyz_0_xy_1[i] * fz_be_0 + g_yyyzz_0_xy_1[i] * wa_z[i];

        g_yyyzzz_0_xz_0[i] = 2.0 * g_yzzz_0_xz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xz_1[i] * fz_be_0 + g_yyzzz_0_xz_1[i] * wa_y[i];

        g_yyyzzz_0_yy_0[i] = 2.0 * g_yyyz_0_yy_0[i] * fbe_0 - 2.0 * g_yyyz_0_yy_1[i] * fz_be_0 + g_yyyzz_0_yy_1[i] * wa_z[i];

        g_yyyzzz_0_yz_0[i] = 2.0 * g_yzzz_0_yz_0[i] * fbe_0 - 2.0 * g_yzzz_0_yz_1[i] * fz_be_0 + g_yyzzz_0_z_1[i] * fi_acd_0 + g_yyzzz_0_yz_1[i] * wa_y[i];

        g_yyyzzz_0_zz_0[i] = 2.0 * g_yzzz_0_zz_0[i] * fbe_0 - 2.0 * g_yzzz_0_zz_1[i] * fz_be_0 + g_yyzzz_0_zz_1[i] * wa_y[i];
    }

    /// Set up 150-156 components of targeted buffer : ISD

    auto g_yyzzzz_0_xx_0 = pbuffer.data(idx_eri_0_isd + 150);

    auto g_yyzzzz_0_xy_0 = pbuffer.data(idx_eri_0_isd + 151);

    auto g_yyzzzz_0_xz_0 = pbuffer.data(idx_eri_0_isd + 152);

    auto g_yyzzzz_0_yy_0 = pbuffer.data(idx_eri_0_isd + 153);

    auto g_yyzzzz_0_yz_0 = pbuffer.data(idx_eri_0_isd + 154);

    auto g_yyzzzz_0_zz_0 = pbuffer.data(idx_eri_0_isd + 155);

    #pragma omp simd aligned(g_yyzz_0_xy_0, g_yyzz_0_xy_1, g_yyzz_0_yy_0, g_yyzz_0_yy_1, g_yyzzz_0_xy_1, g_yyzzz_0_yy_1, g_yyzzzz_0_xx_0, g_yyzzzz_0_xy_0, g_yyzzzz_0_xz_0, g_yyzzzz_0_yy_0, g_yyzzzz_0_yz_0, g_yyzzzz_0_zz_0, g_yzzzz_0_xx_1, g_yzzzz_0_xz_1, g_yzzzz_0_yz_1, g_yzzzz_0_z_1, g_yzzzz_0_zz_1, g_zzzz_0_xx_0, g_zzzz_0_xx_1, g_zzzz_0_xz_0, g_zzzz_0_xz_1, g_zzzz_0_yz_0, g_zzzz_0_yz_1, g_zzzz_0_zz_0, g_zzzz_0_zz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyzzzz_0_xx_0[i] = g_zzzz_0_xx_0[i] * fbe_0 - g_zzzz_0_xx_1[i] * fz_be_0 + g_yzzzz_0_xx_1[i] * wa_y[i];

        g_yyzzzz_0_xy_0[i] = 3.0 * g_yyzz_0_xy_0[i] * fbe_0 - 3.0 * g_yyzz_0_xy_1[i] * fz_be_0 + g_yyzzz_0_xy_1[i] * wa_z[i];

        g_yyzzzz_0_xz_0[i] = g_zzzz_0_xz_0[i] * fbe_0 - g_zzzz_0_xz_1[i] * fz_be_0 + g_yzzzz_0_xz_1[i] * wa_y[i];

        g_yyzzzz_0_yy_0[i] = 3.0 * g_yyzz_0_yy_0[i] * fbe_0 - 3.0 * g_yyzz_0_yy_1[i] * fz_be_0 + g_yyzzz_0_yy_1[i] * wa_z[i];

        g_yyzzzz_0_yz_0[i] = g_zzzz_0_yz_0[i] * fbe_0 - g_zzzz_0_yz_1[i] * fz_be_0 + g_yzzzz_0_z_1[i] * fi_acd_0 + g_yzzzz_0_yz_1[i] * wa_y[i];

        g_yyzzzz_0_zz_0[i] = g_zzzz_0_zz_0[i] * fbe_0 - g_zzzz_0_zz_1[i] * fz_be_0 + g_yzzzz_0_zz_1[i] * wa_y[i];
    }

    /// Set up 156-162 components of targeted buffer : ISD

    auto g_yzzzzz_0_xx_0 = pbuffer.data(idx_eri_0_isd + 156);

    auto g_yzzzzz_0_xy_0 = pbuffer.data(idx_eri_0_isd + 157);

    auto g_yzzzzz_0_xz_0 = pbuffer.data(idx_eri_0_isd + 158);

    auto g_yzzzzz_0_yy_0 = pbuffer.data(idx_eri_0_isd + 159);

    auto g_yzzzzz_0_yz_0 = pbuffer.data(idx_eri_0_isd + 160);

    auto g_yzzzzz_0_zz_0 = pbuffer.data(idx_eri_0_isd + 161);

    #pragma omp simd aligned(g_yzzzzz_0_xx_0, g_yzzzzz_0_xy_0, g_yzzzzz_0_xz_0, g_yzzzzz_0_yy_0, g_yzzzzz_0_yz_0, g_yzzzzz_0_zz_0, g_zzzzz_0_x_1, g_zzzzz_0_xx_1, g_zzzzz_0_xy_1, g_zzzzz_0_xz_1, g_zzzzz_0_y_1, g_zzzzz_0_yy_1, g_zzzzz_0_yz_1, g_zzzzz_0_z_1, g_zzzzz_0_zz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzzzzz_0_xx_0[i] = g_zzzzz_0_xx_1[i] * wa_y[i];

        g_yzzzzz_0_xy_0[i] = g_zzzzz_0_x_1[i] * fi_acd_0 + g_zzzzz_0_xy_1[i] * wa_y[i];

        g_yzzzzz_0_xz_0[i] = g_zzzzz_0_xz_1[i] * wa_y[i];

        g_yzzzzz_0_yy_0[i] = 2.0 * g_zzzzz_0_y_1[i] * fi_acd_0 + g_zzzzz_0_yy_1[i] * wa_y[i];

        g_yzzzzz_0_yz_0[i] = g_zzzzz_0_z_1[i] * fi_acd_0 + g_zzzzz_0_yz_1[i] * wa_y[i];

        g_yzzzzz_0_zz_0[i] = g_zzzzz_0_zz_1[i] * wa_y[i];
    }

    /// Set up 162-168 components of targeted buffer : ISD

    auto g_zzzzzz_0_xx_0 = pbuffer.data(idx_eri_0_isd + 162);

    auto g_zzzzzz_0_xy_0 = pbuffer.data(idx_eri_0_isd + 163);

    auto g_zzzzzz_0_xz_0 = pbuffer.data(idx_eri_0_isd + 164);

    auto g_zzzzzz_0_yy_0 = pbuffer.data(idx_eri_0_isd + 165);

    auto g_zzzzzz_0_yz_0 = pbuffer.data(idx_eri_0_isd + 166);

    auto g_zzzzzz_0_zz_0 = pbuffer.data(idx_eri_0_isd + 167);

    #pragma omp simd aligned(g_zzzz_0_xx_0, g_zzzz_0_xx_1, g_zzzz_0_xy_0, g_zzzz_0_xy_1, g_zzzz_0_xz_0, g_zzzz_0_xz_1, g_zzzz_0_yy_0, g_zzzz_0_yy_1, g_zzzz_0_yz_0, g_zzzz_0_yz_1, g_zzzz_0_zz_0, g_zzzz_0_zz_1, g_zzzzz_0_x_1, g_zzzzz_0_xx_1, g_zzzzz_0_xy_1, g_zzzzz_0_xz_1, g_zzzzz_0_y_1, g_zzzzz_0_yy_1, g_zzzzz_0_yz_1, g_zzzzz_0_z_1, g_zzzzz_0_zz_1, g_zzzzzz_0_xx_0, g_zzzzzz_0_xy_0, g_zzzzzz_0_xz_0, g_zzzzzz_0_yy_0, g_zzzzzz_0_yz_0, g_zzzzzz_0_zz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzzzzz_0_xx_0[i] = 5.0 * g_zzzz_0_xx_0[i] * fbe_0 - 5.0 * g_zzzz_0_xx_1[i] * fz_be_0 + g_zzzzz_0_xx_1[i] * wa_z[i];

        g_zzzzzz_0_xy_0[i] = 5.0 * g_zzzz_0_xy_0[i] * fbe_0 - 5.0 * g_zzzz_0_xy_1[i] * fz_be_0 + g_zzzzz_0_xy_1[i] * wa_z[i];

        g_zzzzzz_0_xz_0[i] = 5.0 * g_zzzz_0_xz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xz_1[i] * fz_be_0 + g_zzzzz_0_x_1[i] * fi_acd_0 + g_zzzzz_0_xz_1[i] * wa_z[i];

        g_zzzzzz_0_yy_0[i] = 5.0 * g_zzzz_0_yy_0[i] * fbe_0 - 5.0 * g_zzzz_0_yy_1[i] * fz_be_0 + g_zzzzz_0_yy_1[i] * wa_z[i];

        g_zzzzzz_0_yz_0[i] = 5.0 * g_zzzz_0_yz_0[i] * fbe_0 - 5.0 * g_zzzz_0_yz_1[i] * fz_be_0 + g_zzzzz_0_y_1[i] * fi_acd_0 + g_zzzzz_0_yz_1[i] * wa_z[i];

        g_zzzzzz_0_zz_0[i] = 5.0 * g_zzzz_0_zz_0[i] * fbe_0 - 5.0 * g_zzzz_0_zz_1[i] * fz_be_0 + 2.0 * g_zzzzz_0_z_1[i] * fi_acd_0 + g_zzzzz_0_zz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace


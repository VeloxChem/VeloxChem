#include "TwoCenterElectronRepulsionPrimRecGH.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_gh(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_gh,
                                const size_t idx_eri_0_dh,
                                const size_t idx_eri_1_dh,
                                const size_t idx_eri_1_fg,
                                const size_t idx_eri_1_fh,
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

    // Set up components of auxiliary buffer : DH

    auto g_xx_xxxxx_0 = pbuffer.data(idx_eri_0_dh);

    auto g_xx_xxxxy_0 = pbuffer.data(idx_eri_0_dh + 1);

    auto g_xx_xxxxz_0 = pbuffer.data(idx_eri_0_dh + 2);

    auto g_xx_xxxyy_0 = pbuffer.data(idx_eri_0_dh + 3);

    auto g_xx_xxxyz_0 = pbuffer.data(idx_eri_0_dh + 4);

    auto g_xx_xxxzz_0 = pbuffer.data(idx_eri_0_dh + 5);

    auto g_xx_xxyyy_0 = pbuffer.data(idx_eri_0_dh + 6);

    auto g_xx_xxyyz_0 = pbuffer.data(idx_eri_0_dh + 7);

    auto g_xx_xxyzz_0 = pbuffer.data(idx_eri_0_dh + 8);

    auto g_xx_xxzzz_0 = pbuffer.data(idx_eri_0_dh + 9);

    auto g_xx_xyyyy_0 = pbuffer.data(idx_eri_0_dh + 10);

    auto g_xx_xyyyz_0 = pbuffer.data(idx_eri_0_dh + 11);

    auto g_xx_xyyzz_0 = pbuffer.data(idx_eri_0_dh + 12);

    auto g_xx_xyzzz_0 = pbuffer.data(idx_eri_0_dh + 13);

    auto g_xx_xzzzz_0 = pbuffer.data(idx_eri_0_dh + 14);

    auto g_xx_yyyyy_0 = pbuffer.data(idx_eri_0_dh + 15);

    auto g_xx_yyyyz_0 = pbuffer.data(idx_eri_0_dh + 16);

    auto g_xx_yyyzz_0 = pbuffer.data(idx_eri_0_dh + 17);

    auto g_xx_yyzzz_0 = pbuffer.data(idx_eri_0_dh + 18);

    auto g_xx_yzzzz_0 = pbuffer.data(idx_eri_0_dh + 19);

    auto g_xx_zzzzz_0 = pbuffer.data(idx_eri_0_dh + 20);

    auto g_yy_xxxxx_0 = pbuffer.data(idx_eri_0_dh + 63);

    auto g_yy_xxxxy_0 = pbuffer.data(idx_eri_0_dh + 64);

    auto g_yy_xxxxz_0 = pbuffer.data(idx_eri_0_dh + 65);

    auto g_yy_xxxyy_0 = pbuffer.data(idx_eri_0_dh + 66);

    auto g_yy_xxxyz_0 = pbuffer.data(idx_eri_0_dh + 67);

    auto g_yy_xxxzz_0 = pbuffer.data(idx_eri_0_dh + 68);

    auto g_yy_xxyyy_0 = pbuffer.data(idx_eri_0_dh + 69);

    auto g_yy_xxyyz_0 = pbuffer.data(idx_eri_0_dh + 70);

    auto g_yy_xxyzz_0 = pbuffer.data(idx_eri_0_dh + 71);

    auto g_yy_xxzzz_0 = pbuffer.data(idx_eri_0_dh + 72);

    auto g_yy_xyyyy_0 = pbuffer.data(idx_eri_0_dh + 73);

    auto g_yy_xyyyz_0 = pbuffer.data(idx_eri_0_dh + 74);

    auto g_yy_xyyzz_0 = pbuffer.data(idx_eri_0_dh + 75);

    auto g_yy_xyzzz_0 = pbuffer.data(idx_eri_0_dh + 76);

    auto g_yy_xzzzz_0 = pbuffer.data(idx_eri_0_dh + 77);

    auto g_yy_yyyyy_0 = pbuffer.data(idx_eri_0_dh + 78);

    auto g_yy_yyyyz_0 = pbuffer.data(idx_eri_0_dh + 79);

    auto g_yy_yyyzz_0 = pbuffer.data(idx_eri_0_dh + 80);

    auto g_yy_yyzzz_0 = pbuffer.data(idx_eri_0_dh + 81);

    auto g_yy_yzzzz_0 = pbuffer.data(idx_eri_0_dh + 82);

    auto g_yy_zzzzz_0 = pbuffer.data(idx_eri_0_dh + 83);

    auto g_zz_xxxxx_0 = pbuffer.data(idx_eri_0_dh + 105);

    auto g_zz_xxxxy_0 = pbuffer.data(idx_eri_0_dh + 106);

    auto g_zz_xxxxz_0 = pbuffer.data(idx_eri_0_dh + 107);

    auto g_zz_xxxyy_0 = pbuffer.data(idx_eri_0_dh + 108);

    auto g_zz_xxxyz_0 = pbuffer.data(idx_eri_0_dh + 109);

    auto g_zz_xxxzz_0 = pbuffer.data(idx_eri_0_dh + 110);

    auto g_zz_xxyyy_0 = pbuffer.data(idx_eri_0_dh + 111);

    auto g_zz_xxyyz_0 = pbuffer.data(idx_eri_0_dh + 112);

    auto g_zz_xxyzz_0 = pbuffer.data(idx_eri_0_dh + 113);

    auto g_zz_xxzzz_0 = pbuffer.data(idx_eri_0_dh + 114);

    auto g_zz_xyyyy_0 = pbuffer.data(idx_eri_0_dh + 115);

    auto g_zz_xyyyz_0 = pbuffer.data(idx_eri_0_dh + 116);

    auto g_zz_xyyzz_0 = pbuffer.data(idx_eri_0_dh + 117);

    auto g_zz_xyzzz_0 = pbuffer.data(idx_eri_0_dh + 118);

    auto g_zz_xzzzz_0 = pbuffer.data(idx_eri_0_dh + 119);

    auto g_zz_yyyyy_0 = pbuffer.data(idx_eri_0_dh + 120);

    auto g_zz_yyyyz_0 = pbuffer.data(idx_eri_0_dh + 121);

    auto g_zz_yyyzz_0 = pbuffer.data(idx_eri_0_dh + 122);

    auto g_zz_yyzzz_0 = pbuffer.data(idx_eri_0_dh + 123);

    auto g_zz_yzzzz_0 = pbuffer.data(idx_eri_0_dh + 124);

    auto g_zz_zzzzz_0 = pbuffer.data(idx_eri_0_dh + 125);

    // Set up components of auxiliary buffer : DH

    auto g_xx_xxxxx_1 = pbuffer.data(idx_eri_1_dh);

    auto g_xx_xxxxy_1 = pbuffer.data(idx_eri_1_dh + 1);

    auto g_xx_xxxxz_1 = pbuffer.data(idx_eri_1_dh + 2);

    auto g_xx_xxxyy_1 = pbuffer.data(idx_eri_1_dh + 3);

    auto g_xx_xxxyz_1 = pbuffer.data(idx_eri_1_dh + 4);

    auto g_xx_xxxzz_1 = pbuffer.data(idx_eri_1_dh + 5);

    auto g_xx_xxyyy_1 = pbuffer.data(idx_eri_1_dh + 6);

    auto g_xx_xxyyz_1 = pbuffer.data(idx_eri_1_dh + 7);

    auto g_xx_xxyzz_1 = pbuffer.data(idx_eri_1_dh + 8);

    auto g_xx_xxzzz_1 = pbuffer.data(idx_eri_1_dh + 9);

    auto g_xx_xyyyy_1 = pbuffer.data(idx_eri_1_dh + 10);

    auto g_xx_xyyyz_1 = pbuffer.data(idx_eri_1_dh + 11);

    auto g_xx_xyyzz_1 = pbuffer.data(idx_eri_1_dh + 12);

    auto g_xx_xyzzz_1 = pbuffer.data(idx_eri_1_dh + 13);

    auto g_xx_xzzzz_1 = pbuffer.data(idx_eri_1_dh + 14);

    auto g_xx_yyyyy_1 = pbuffer.data(idx_eri_1_dh + 15);

    auto g_xx_yyyyz_1 = pbuffer.data(idx_eri_1_dh + 16);

    auto g_xx_yyyzz_1 = pbuffer.data(idx_eri_1_dh + 17);

    auto g_xx_yyzzz_1 = pbuffer.data(idx_eri_1_dh + 18);

    auto g_xx_yzzzz_1 = pbuffer.data(idx_eri_1_dh + 19);

    auto g_xx_zzzzz_1 = pbuffer.data(idx_eri_1_dh + 20);

    auto g_yy_xxxxx_1 = pbuffer.data(idx_eri_1_dh + 63);

    auto g_yy_xxxxy_1 = pbuffer.data(idx_eri_1_dh + 64);

    auto g_yy_xxxxz_1 = pbuffer.data(idx_eri_1_dh + 65);

    auto g_yy_xxxyy_1 = pbuffer.data(idx_eri_1_dh + 66);

    auto g_yy_xxxyz_1 = pbuffer.data(idx_eri_1_dh + 67);

    auto g_yy_xxxzz_1 = pbuffer.data(idx_eri_1_dh + 68);

    auto g_yy_xxyyy_1 = pbuffer.data(idx_eri_1_dh + 69);

    auto g_yy_xxyyz_1 = pbuffer.data(idx_eri_1_dh + 70);

    auto g_yy_xxyzz_1 = pbuffer.data(idx_eri_1_dh + 71);

    auto g_yy_xxzzz_1 = pbuffer.data(idx_eri_1_dh + 72);

    auto g_yy_xyyyy_1 = pbuffer.data(idx_eri_1_dh + 73);

    auto g_yy_xyyyz_1 = pbuffer.data(idx_eri_1_dh + 74);

    auto g_yy_xyyzz_1 = pbuffer.data(idx_eri_1_dh + 75);

    auto g_yy_xyzzz_1 = pbuffer.data(idx_eri_1_dh + 76);

    auto g_yy_xzzzz_1 = pbuffer.data(idx_eri_1_dh + 77);

    auto g_yy_yyyyy_1 = pbuffer.data(idx_eri_1_dh + 78);

    auto g_yy_yyyyz_1 = pbuffer.data(idx_eri_1_dh + 79);

    auto g_yy_yyyzz_1 = pbuffer.data(idx_eri_1_dh + 80);

    auto g_yy_yyzzz_1 = pbuffer.data(idx_eri_1_dh + 81);

    auto g_yy_yzzzz_1 = pbuffer.data(idx_eri_1_dh + 82);

    auto g_yy_zzzzz_1 = pbuffer.data(idx_eri_1_dh + 83);

    auto g_zz_xxxxx_1 = pbuffer.data(idx_eri_1_dh + 105);

    auto g_zz_xxxxy_1 = pbuffer.data(idx_eri_1_dh + 106);

    auto g_zz_xxxxz_1 = pbuffer.data(idx_eri_1_dh + 107);

    auto g_zz_xxxyy_1 = pbuffer.data(idx_eri_1_dh + 108);

    auto g_zz_xxxyz_1 = pbuffer.data(idx_eri_1_dh + 109);

    auto g_zz_xxxzz_1 = pbuffer.data(idx_eri_1_dh + 110);

    auto g_zz_xxyyy_1 = pbuffer.data(idx_eri_1_dh + 111);

    auto g_zz_xxyyz_1 = pbuffer.data(idx_eri_1_dh + 112);

    auto g_zz_xxyzz_1 = pbuffer.data(idx_eri_1_dh + 113);

    auto g_zz_xxzzz_1 = pbuffer.data(idx_eri_1_dh + 114);

    auto g_zz_xyyyy_1 = pbuffer.data(idx_eri_1_dh + 115);

    auto g_zz_xyyyz_1 = pbuffer.data(idx_eri_1_dh + 116);

    auto g_zz_xyyzz_1 = pbuffer.data(idx_eri_1_dh + 117);

    auto g_zz_xyzzz_1 = pbuffer.data(idx_eri_1_dh + 118);

    auto g_zz_xzzzz_1 = pbuffer.data(idx_eri_1_dh + 119);

    auto g_zz_yyyyy_1 = pbuffer.data(idx_eri_1_dh + 120);

    auto g_zz_yyyyz_1 = pbuffer.data(idx_eri_1_dh + 121);

    auto g_zz_yyyzz_1 = pbuffer.data(idx_eri_1_dh + 122);

    auto g_zz_yyzzz_1 = pbuffer.data(idx_eri_1_dh + 123);

    auto g_zz_yzzzz_1 = pbuffer.data(idx_eri_1_dh + 124);

    auto g_zz_zzzzz_1 = pbuffer.data(idx_eri_1_dh + 125);

    // Set up components of auxiliary buffer : FG

    auto g_xxx_xxxx_1 = pbuffer.data(idx_eri_1_fg);

    auto g_xxx_xxxy_1 = pbuffer.data(idx_eri_1_fg + 1);

    auto g_xxx_xxxz_1 = pbuffer.data(idx_eri_1_fg + 2);

    auto g_xxx_xxyy_1 = pbuffer.data(idx_eri_1_fg + 3);

    auto g_xxx_xxyz_1 = pbuffer.data(idx_eri_1_fg + 4);

    auto g_xxx_xxzz_1 = pbuffer.data(idx_eri_1_fg + 5);

    auto g_xxx_xyyy_1 = pbuffer.data(idx_eri_1_fg + 6);

    auto g_xxx_xyyz_1 = pbuffer.data(idx_eri_1_fg + 7);

    auto g_xxx_xyzz_1 = pbuffer.data(idx_eri_1_fg + 8);

    auto g_xxx_xzzz_1 = pbuffer.data(idx_eri_1_fg + 9);

    auto g_xxx_yyyy_1 = pbuffer.data(idx_eri_1_fg + 10);

    auto g_xxx_yyyz_1 = pbuffer.data(idx_eri_1_fg + 11);

    auto g_xxx_yyzz_1 = pbuffer.data(idx_eri_1_fg + 12);

    auto g_xxx_yzzz_1 = pbuffer.data(idx_eri_1_fg + 13);

    auto g_xxx_zzzz_1 = pbuffer.data(idx_eri_1_fg + 14);

    auto g_xxz_xxxz_1 = pbuffer.data(idx_eri_1_fg + 32);

    auto g_xxz_xxyz_1 = pbuffer.data(idx_eri_1_fg + 34);

    auto g_xxz_xxzz_1 = pbuffer.data(idx_eri_1_fg + 35);

    auto g_xxz_xyyz_1 = pbuffer.data(idx_eri_1_fg + 37);

    auto g_xxz_xyzz_1 = pbuffer.data(idx_eri_1_fg + 38);

    auto g_xxz_xzzz_1 = pbuffer.data(idx_eri_1_fg + 39);

    auto g_xxz_yyyz_1 = pbuffer.data(idx_eri_1_fg + 41);

    auto g_xxz_yyzz_1 = pbuffer.data(idx_eri_1_fg + 42);

    auto g_xxz_yzzz_1 = pbuffer.data(idx_eri_1_fg + 43);

    auto g_xxz_zzzz_1 = pbuffer.data(idx_eri_1_fg + 44);

    auto g_xyy_xxxy_1 = pbuffer.data(idx_eri_1_fg + 46);

    auto g_xyy_xxyy_1 = pbuffer.data(idx_eri_1_fg + 48);

    auto g_xyy_xxyz_1 = pbuffer.data(idx_eri_1_fg + 49);

    auto g_xyy_xyyy_1 = pbuffer.data(idx_eri_1_fg + 51);

    auto g_xyy_xyyz_1 = pbuffer.data(idx_eri_1_fg + 52);

    auto g_xyy_xyzz_1 = pbuffer.data(idx_eri_1_fg + 53);

    auto g_xyy_yyyy_1 = pbuffer.data(idx_eri_1_fg + 55);

    auto g_xyy_yyyz_1 = pbuffer.data(idx_eri_1_fg + 56);

    auto g_xyy_yyzz_1 = pbuffer.data(idx_eri_1_fg + 57);

    auto g_xyy_yzzz_1 = pbuffer.data(idx_eri_1_fg + 58);

    auto g_xzz_xxxz_1 = pbuffer.data(idx_eri_1_fg + 77);

    auto g_xzz_xxyz_1 = pbuffer.data(idx_eri_1_fg + 79);

    auto g_xzz_xxzz_1 = pbuffer.data(idx_eri_1_fg + 80);

    auto g_xzz_xyyz_1 = pbuffer.data(idx_eri_1_fg + 82);

    auto g_xzz_xyzz_1 = pbuffer.data(idx_eri_1_fg + 83);

    auto g_xzz_xzzz_1 = pbuffer.data(idx_eri_1_fg + 84);

    auto g_xzz_yyyz_1 = pbuffer.data(idx_eri_1_fg + 86);

    auto g_xzz_yyzz_1 = pbuffer.data(idx_eri_1_fg + 87);

    auto g_xzz_yzzz_1 = pbuffer.data(idx_eri_1_fg + 88);

    auto g_xzz_zzzz_1 = pbuffer.data(idx_eri_1_fg + 89);

    auto g_yyy_xxxx_1 = pbuffer.data(idx_eri_1_fg + 90);

    auto g_yyy_xxxy_1 = pbuffer.data(idx_eri_1_fg + 91);

    auto g_yyy_xxxz_1 = pbuffer.data(idx_eri_1_fg + 92);

    auto g_yyy_xxyy_1 = pbuffer.data(idx_eri_1_fg + 93);

    auto g_yyy_xxyz_1 = pbuffer.data(idx_eri_1_fg + 94);

    auto g_yyy_xxzz_1 = pbuffer.data(idx_eri_1_fg + 95);

    auto g_yyy_xyyy_1 = pbuffer.data(idx_eri_1_fg + 96);

    auto g_yyy_xyyz_1 = pbuffer.data(idx_eri_1_fg + 97);

    auto g_yyy_xyzz_1 = pbuffer.data(idx_eri_1_fg + 98);

    auto g_yyy_xzzz_1 = pbuffer.data(idx_eri_1_fg + 99);

    auto g_yyy_yyyy_1 = pbuffer.data(idx_eri_1_fg + 100);

    auto g_yyy_yyyz_1 = pbuffer.data(idx_eri_1_fg + 101);

    auto g_yyy_yyzz_1 = pbuffer.data(idx_eri_1_fg + 102);

    auto g_yyy_yzzz_1 = pbuffer.data(idx_eri_1_fg + 103);

    auto g_yyy_zzzz_1 = pbuffer.data(idx_eri_1_fg + 104);

    auto g_yyz_xxxz_1 = pbuffer.data(idx_eri_1_fg + 107);

    auto g_yyz_xxyz_1 = pbuffer.data(idx_eri_1_fg + 109);

    auto g_yyz_xxzz_1 = pbuffer.data(idx_eri_1_fg + 110);

    auto g_yyz_xyyz_1 = pbuffer.data(idx_eri_1_fg + 112);

    auto g_yyz_xyzz_1 = pbuffer.data(idx_eri_1_fg + 113);

    auto g_yyz_xzzz_1 = pbuffer.data(idx_eri_1_fg + 114);

    auto g_yyz_yyyz_1 = pbuffer.data(idx_eri_1_fg + 116);

    auto g_yyz_yyzz_1 = pbuffer.data(idx_eri_1_fg + 117);

    auto g_yyz_yzzz_1 = pbuffer.data(idx_eri_1_fg + 118);

    auto g_yyz_zzzz_1 = pbuffer.data(idx_eri_1_fg + 119);

    auto g_yzz_xxxy_1 = pbuffer.data(idx_eri_1_fg + 121);

    auto g_yzz_xxxz_1 = pbuffer.data(idx_eri_1_fg + 122);

    auto g_yzz_xxyy_1 = pbuffer.data(idx_eri_1_fg + 123);

    auto g_yzz_xxyz_1 = pbuffer.data(idx_eri_1_fg + 124);

    auto g_yzz_xxzz_1 = pbuffer.data(idx_eri_1_fg + 125);

    auto g_yzz_xyyy_1 = pbuffer.data(idx_eri_1_fg + 126);

    auto g_yzz_xyyz_1 = pbuffer.data(idx_eri_1_fg + 127);

    auto g_yzz_xyzz_1 = pbuffer.data(idx_eri_1_fg + 128);

    auto g_yzz_xzzz_1 = pbuffer.data(idx_eri_1_fg + 129);

    auto g_yzz_yyyy_1 = pbuffer.data(idx_eri_1_fg + 130);

    auto g_yzz_yyyz_1 = pbuffer.data(idx_eri_1_fg + 131);

    auto g_yzz_yyzz_1 = pbuffer.data(idx_eri_1_fg + 132);

    auto g_yzz_yzzz_1 = pbuffer.data(idx_eri_1_fg + 133);

    auto g_yzz_zzzz_1 = pbuffer.data(idx_eri_1_fg + 134);

    auto g_zzz_xxxx_1 = pbuffer.data(idx_eri_1_fg + 135);

    auto g_zzz_xxxy_1 = pbuffer.data(idx_eri_1_fg + 136);

    auto g_zzz_xxxz_1 = pbuffer.data(idx_eri_1_fg + 137);

    auto g_zzz_xxyy_1 = pbuffer.data(idx_eri_1_fg + 138);

    auto g_zzz_xxyz_1 = pbuffer.data(idx_eri_1_fg + 139);

    auto g_zzz_xxzz_1 = pbuffer.data(idx_eri_1_fg + 140);

    auto g_zzz_xyyy_1 = pbuffer.data(idx_eri_1_fg + 141);

    auto g_zzz_xyyz_1 = pbuffer.data(idx_eri_1_fg + 142);

    auto g_zzz_xyzz_1 = pbuffer.data(idx_eri_1_fg + 143);

    auto g_zzz_xzzz_1 = pbuffer.data(idx_eri_1_fg + 144);

    auto g_zzz_yyyy_1 = pbuffer.data(idx_eri_1_fg + 145);

    auto g_zzz_yyyz_1 = pbuffer.data(idx_eri_1_fg + 146);

    auto g_zzz_yyzz_1 = pbuffer.data(idx_eri_1_fg + 147);

    auto g_zzz_yzzz_1 = pbuffer.data(idx_eri_1_fg + 148);

    auto g_zzz_zzzz_1 = pbuffer.data(idx_eri_1_fg + 149);

    // Set up components of auxiliary buffer : FH

    auto g_xxx_xxxxx_1 = pbuffer.data(idx_eri_1_fh);

    auto g_xxx_xxxxy_1 = pbuffer.data(idx_eri_1_fh + 1);

    auto g_xxx_xxxxz_1 = pbuffer.data(idx_eri_1_fh + 2);

    auto g_xxx_xxxyy_1 = pbuffer.data(idx_eri_1_fh + 3);

    auto g_xxx_xxxyz_1 = pbuffer.data(idx_eri_1_fh + 4);

    auto g_xxx_xxxzz_1 = pbuffer.data(idx_eri_1_fh + 5);

    auto g_xxx_xxyyy_1 = pbuffer.data(idx_eri_1_fh + 6);

    auto g_xxx_xxyyz_1 = pbuffer.data(idx_eri_1_fh + 7);

    auto g_xxx_xxyzz_1 = pbuffer.data(idx_eri_1_fh + 8);

    auto g_xxx_xxzzz_1 = pbuffer.data(idx_eri_1_fh + 9);

    auto g_xxx_xyyyy_1 = pbuffer.data(idx_eri_1_fh + 10);

    auto g_xxx_xyyyz_1 = pbuffer.data(idx_eri_1_fh + 11);

    auto g_xxx_xyyzz_1 = pbuffer.data(idx_eri_1_fh + 12);

    auto g_xxx_xyzzz_1 = pbuffer.data(idx_eri_1_fh + 13);

    auto g_xxx_xzzzz_1 = pbuffer.data(idx_eri_1_fh + 14);

    auto g_xxx_yyyyy_1 = pbuffer.data(idx_eri_1_fh + 15);

    auto g_xxx_yyyyz_1 = pbuffer.data(idx_eri_1_fh + 16);

    auto g_xxx_yyyzz_1 = pbuffer.data(idx_eri_1_fh + 17);

    auto g_xxx_yyzzz_1 = pbuffer.data(idx_eri_1_fh + 18);

    auto g_xxx_yzzzz_1 = pbuffer.data(idx_eri_1_fh + 19);

    auto g_xxx_zzzzz_1 = pbuffer.data(idx_eri_1_fh + 20);

    auto g_xxy_xxxxx_1 = pbuffer.data(idx_eri_1_fh + 21);

    auto g_xxy_xxxxy_1 = pbuffer.data(idx_eri_1_fh + 22);

    auto g_xxy_xxxxz_1 = pbuffer.data(idx_eri_1_fh + 23);

    auto g_xxy_xxxyy_1 = pbuffer.data(idx_eri_1_fh + 24);

    auto g_xxy_xxxzz_1 = pbuffer.data(idx_eri_1_fh + 26);

    auto g_xxy_xxyyy_1 = pbuffer.data(idx_eri_1_fh + 27);

    auto g_xxy_xxzzz_1 = pbuffer.data(idx_eri_1_fh + 30);

    auto g_xxy_xyyyy_1 = pbuffer.data(idx_eri_1_fh + 31);

    auto g_xxy_xzzzz_1 = pbuffer.data(idx_eri_1_fh + 35);

    auto g_xxy_yyyyy_1 = pbuffer.data(idx_eri_1_fh + 36);

    auto g_xxz_xxxxx_1 = pbuffer.data(idx_eri_1_fh + 42);

    auto g_xxz_xxxxy_1 = pbuffer.data(idx_eri_1_fh + 43);

    auto g_xxz_xxxxz_1 = pbuffer.data(idx_eri_1_fh + 44);

    auto g_xxz_xxxyy_1 = pbuffer.data(idx_eri_1_fh + 45);

    auto g_xxz_xxxyz_1 = pbuffer.data(idx_eri_1_fh + 46);

    auto g_xxz_xxxzz_1 = pbuffer.data(idx_eri_1_fh + 47);

    auto g_xxz_xxyyy_1 = pbuffer.data(idx_eri_1_fh + 48);

    auto g_xxz_xxyyz_1 = pbuffer.data(idx_eri_1_fh + 49);

    auto g_xxz_xxyzz_1 = pbuffer.data(idx_eri_1_fh + 50);

    auto g_xxz_xxzzz_1 = pbuffer.data(idx_eri_1_fh + 51);

    auto g_xxz_xyyyy_1 = pbuffer.data(idx_eri_1_fh + 52);

    auto g_xxz_xyyyz_1 = pbuffer.data(idx_eri_1_fh + 53);

    auto g_xxz_xyyzz_1 = pbuffer.data(idx_eri_1_fh + 54);

    auto g_xxz_xyzzz_1 = pbuffer.data(idx_eri_1_fh + 55);

    auto g_xxz_xzzzz_1 = pbuffer.data(idx_eri_1_fh + 56);

    auto g_xxz_yyyyz_1 = pbuffer.data(idx_eri_1_fh + 58);

    auto g_xxz_yyyzz_1 = pbuffer.data(idx_eri_1_fh + 59);

    auto g_xxz_yyzzz_1 = pbuffer.data(idx_eri_1_fh + 60);

    auto g_xxz_yzzzz_1 = pbuffer.data(idx_eri_1_fh + 61);

    auto g_xxz_zzzzz_1 = pbuffer.data(idx_eri_1_fh + 62);

    auto g_xyy_xxxxx_1 = pbuffer.data(idx_eri_1_fh + 63);

    auto g_xyy_xxxxy_1 = pbuffer.data(idx_eri_1_fh + 64);

    auto g_xyy_xxxyy_1 = pbuffer.data(idx_eri_1_fh + 66);

    auto g_xyy_xxxyz_1 = pbuffer.data(idx_eri_1_fh + 67);

    auto g_xyy_xxyyy_1 = pbuffer.data(idx_eri_1_fh + 69);

    auto g_xyy_xxyyz_1 = pbuffer.data(idx_eri_1_fh + 70);

    auto g_xyy_xxyzz_1 = pbuffer.data(idx_eri_1_fh + 71);

    auto g_xyy_xyyyy_1 = pbuffer.data(idx_eri_1_fh + 73);

    auto g_xyy_xyyyz_1 = pbuffer.data(idx_eri_1_fh + 74);

    auto g_xyy_xyyzz_1 = pbuffer.data(idx_eri_1_fh + 75);

    auto g_xyy_xyzzz_1 = pbuffer.data(idx_eri_1_fh + 76);

    auto g_xyy_yyyyy_1 = pbuffer.data(idx_eri_1_fh + 78);

    auto g_xyy_yyyyz_1 = pbuffer.data(idx_eri_1_fh + 79);

    auto g_xyy_yyyzz_1 = pbuffer.data(idx_eri_1_fh + 80);

    auto g_xyy_yyzzz_1 = pbuffer.data(idx_eri_1_fh + 81);

    auto g_xyy_yzzzz_1 = pbuffer.data(idx_eri_1_fh + 82);

    auto g_xyy_zzzzz_1 = pbuffer.data(idx_eri_1_fh + 83);

    auto g_xzz_xxxxx_1 = pbuffer.data(idx_eri_1_fh + 105);

    auto g_xzz_xxxxz_1 = pbuffer.data(idx_eri_1_fh + 107);

    auto g_xzz_xxxyz_1 = pbuffer.data(idx_eri_1_fh + 109);

    auto g_xzz_xxxzz_1 = pbuffer.data(idx_eri_1_fh + 110);

    auto g_xzz_xxyyz_1 = pbuffer.data(idx_eri_1_fh + 112);

    auto g_xzz_xxyzz_1 = pbuffer.data(idx_eri_1_fh + 113);

    auto g_xzz_xxzzz_1 = pbuffer.data(idx_eri_1_fh + 114);

    auto g_xzz_xyyyz_1 = pbuffer.data(idx_eri_1_fh + 116);

    auto g_xzz_xyyzz_1 = pbuffer.data(idx_eri_1_fh + 117);

    auto g_xzz_xyzzz_1 = pbuffer.data(idx_eri_1_fh + 118);

    auto g_xzz_xzzzz_1 = pbuffer.data(idx_eri_1_fh + 119);

    auto g_xzz_yyyyy_1 = pbuffer.data(idx_eri_1_fh + 120);

    auto g_xzz_yyyyz_1 = pbuffer.data(idx_eri_1_fh + 121);

    auto g_xzz_yyyzz_1 = pbuffer.data(idx_eri_1_fh + 122);

    auto g_xzz_yyzzz_1 = pbuffer.data(idx_eri_1_fh + 123);

    auto g_xzz_yzzzz_1 = pbuffer.data(idx_eri_1_fh + 124);

    auto g_xzz_zzzzz_1 = pbuffer.data(idx_eri_1_fh + 125);

    auto g_yyy_xxxxx_1 = pbuffer.data(idx_eri_1_fh + 126);

    auto g_yyy_xxxxy_1 = pbuffer.data(idx_eri_1_fh + 127);

    auto g_yyy_xxxxz_1 = pbuffer.data(idx_eri_1_fh + 128);

    auto g_yyy_xxxyy_1 = pbuffer.data(idx_eri_1_fh + 129);

    auto g_yyy_xxxyz_1 = pbuffer.data(idx_eri_1_fh + 130);

    auto g_yyy_xxxzz_1 = pbuffer.data(idx_eri_1_fh + 131);

    auto g_yyy_xxyyy_1 = pbuffer.data(idx_eri_1_fh + 132);

    auto g_yyy_xxyyz_1 = pbuffer.data(idx_eri_1_fh + 133);

    auto g_yyy_xxyzz_1 = pbuffer.data(idx_eri_1_fh + 134);

    auto g_yyy_xxzzz_1 = pbuffer.data(idx_eri_1_fh + 135);

    auto g_yyy_xyyyy_1 = pbuffer.data(idx_eri_1_fh + 136);

    auto g_yyy_xyyyz_1 = pbuffer.data(idx_eri_1_fh + 137);

    auto g_yyy_xyyzz_1 = pbuffer.data(idx_eri_1_fh + 138);

    auto g_yyy_xyzzz_1 = pbuffer.data(idx_eri_1_fh + 139);

    auto g_yyy_xzzzz_1 = pbuffer.data(idx_eri_1_fh + 140);

    auto g_yyy_yyyyy_1 = pbuffer.data(idx_eri_1_fh + 141);

    auto g_yyy_yyyyz_1 = pbuffer.data(idx_eri_1_fh + 142);

    auto g_yyy_yyyzz_1 = pbuffer.data(idx_eri_1_fh + 143);

    auto g_yyy_yyzzz_1 = pbuffer.data(idx_eri_1_fh + 144);

    auto g_yyy_yzzzz_1 = pbuffer.data(idx_eri_1_fh + 145);

    auto g_yyy_zzzzz_1 = pbuffer.data(idx_eri_1_fh + 146);

    auto g_yyz_xxxxy_1 = pbuffer.data(idx_eri_1_fh + 148);

    auto g_yyz_xxxxz_1 = pbuffer.data(idx_eri_1_fh + 149);

    auto g_yyz_xxxyy_1 = pbuffer.data(idx_eri_1_fh + 150);

    auto g_yyz_xxxyz_1 = pbuffer.data(idx_eri_1_fh + 151);

    auto g_yyz_xxxzz_1 = pbuffer.data(idx_eri_1_fh + 152);

    auto g_yyz_xxyyy_1 = pbuffer.data(idx_eri_1_fh + 153);

    auto g_yyz_xxyyz_1 = pbuffer.data(idx_eri_1_fh + 154);

    auto g_yyz_xxyzz_1 = pbuffer.data(idx_eri_1_fh + 155);

    auto g_yyz_xxzzz_1 = pbuffer.data(idx_eri_1_fh + 156);

    auto g_yyz_xyyyy_1 = pbuffer.data(idx_eri_1_fh + 157);

    auto g_yyz_xyyyz_1 = pbuffer.data(idx_eri_1_fh + 158);

    auto g_yyz_xyyzz_1 = pbuffer.data(idx_eri_1_fh + 159);

    auto g_yyz_xyzzz_1 = pbuffer.data(idx_eri_1_fh + 160);

    auto g_yyz_xzzzz_1 = pbuffer.data(idx_eri_1_fh + 161);

    auto g_yyz_yyyyy_1 = pbuffer.data(idx_eri_1_fh + 162);

    auto g_yyz_yyyyz_1 = pbuffer.data(idx_eri_1_fh + 163);

    auto g_yyz_yyyzz_1 = pbuffer.data(idx_eri_1_fh + 164);

    auto g_yyz_yyzzz_1 = pbuffer.data(idx_eri_1_fh + 165);

    auto g_yyz_yzzzz_1 = pbuffer.data(idx_eri_1_fh + 166);

    auto g_yyz_zzzzz_1 = pbuffer.data(idx_eri_1_fh + 167);

    auto g_yzz_xxxxx_1 = pbuffer.data(idx_eri_1_fh + 168);

    auto g_yzz_xxxxy_1 = pbuffer.data(idx_eri_1_fh + 169);

    auto g_yzz_xxxxz_1 = pbuffer.data(idx_eri_1_fh + 170);

    auto g_yzz_xxxyy_1 = pbuffer.data(idx_eri_1_fh + 171);

    auto g_yzz_xxxyz_1 = pbuffer.data(idx_eri_1_fh + 172);

    auto g_yzz_xxxzz_1 = pbuffer.data(idx_eri_1_fh + 173);

    auto g_yzz_xxyyy_1 = pbuffer.data(idx_eri_1_fh + 174);

    auto g_yzz_xxyyz_1 = pbuffer.data(idx_eri_1_fh + 175);

    auto g_yzz_xxyzz_1 = pbuffer.data(idx_eri_1_fh + 176);

    auto g_yzz_xxzzz_1 = pbuffer.data(idx_eri_1_fh + 177);

    auto g_yzz_xyyyy_1 = pbuffer.data(idx_eri_1_fh + 178);

    auto g_yzz_xyyyz_1 = pbuffer.data(idx_eri_1_fh + 179);

    auto g_yzz_xyyzz_1 = pbuffer.data(idx_eri_1_fh + 180);

    auto g_yzz_xyzzz_1 = pbuffer.data(idx_eri_1_fh + 181);

    auto g_yzz_xzzzz_1 = pbuffer.data(idx_eri_1_fh + 182);

    auto g_yzz_yyyyy_1 = pbuffer.data(idx_eri_1_fh + 183);

    auto g_yzz_yyyyz_1 = pbuffer.data(idx_eri_1_fh + 184);

    auto g_yzz_yyyzz_1 = pbuffer.data(idx_eri_1_fh + 185);

    auto g_yzz_yyzzz_1 = pbuffer.data(idx_eri_1_fh + 186);

    auto g_yzz_yzzzz_1 = pbuffer.data(idx_eri_1_fh + 187);

    auto g_yzz_zzzzz_1 = pbuffer.data(idx_eri_1_fh + 188);

    auto g_zzz_xxxxx_1 = pbuffer.data(idx_eri_1_fh + 189);

    auto g_zzz_xxxxy_1 = pbuffer.data(idx_eri_1_fh + 190);

    auto g_zzz_xxxxz_1 = pbuffer.data(idx_eri_1_fh + 191);

    auto g_zzz_xxxyy_1 = pbuffer.data(idx_eri_1_fh + 192);

    auto g_zzz_xxxyz_1 = pbuffer.data(idx_eri_1_fh + 193);

    auto g_zzz_xxxzz_1 = pbuffer.data(idx_eri_1_fh + 194);

    auto g_zzz_xxyyy_1 = pbuffer.data(idx_eri_1_fh + 195);

    auto g_zzz_xxyyz_1 = pbuffer.data(idx_eri_1_fh + 196);

    auto g_zzz_xxyzz_1 = pbuffer.data(idx_eri_1_fh + 197);

    auto g_zzz_xxzzz_1 = pbuffer.data(idx_eri_1_fh + 198);

    auto g_zzz_xyyyy_1 = pbuffer.data(idx_eri_1_fh + 199);

    auto g_zzz_xyyyz_1 = pbuffer.data(idx_eri_1_fh + 200);

    auto g_zzz_xyyzz_1 = pbuffer.data(idx_eri_1_fh + 201);

    auto g_zzz_xyzzz_1 = pbuffer.data(idx_eri_1_fh + 202);

    auto g_zzz_xzzzz_1 = pbuffer.data(idx_eri_1_fh + 203);

    auto g_zzz_yyyyy_1 = pbuffer.data(idx_eri_1_fh + 204);

    auto g_zzz_yyyyz_1 = pbuffer.data(idx_eri_1_fh + 205);

    auto g_zzz_yyyzz_1 = pbuffer.data(idx_eri_1_fh + 206);

    auto g_zzz_yyzzz_1 = pbuffer.data(idx_eri_1_fh + 207);

    auto g_zzz_yzzzz_1 = pbuffer.data(idx_eri_1_fh + 208);

    auto g_zzz_zzzzz_1 = pbuffer.data(idx_eri_1_fh + 209);

    // Set up 0-21 components of targeted buffer : GH

    auto g_xxxx_xxxxx_0 = pbuffer.data(idx_eri_0_gh);

    auto g_xxxx_xxxxy_0 = pbuffer.data(idx_eri_0_gh + 1);

    auto g_xxxx_xxxxz_0 = pbuffer.data(idx_eri_0_gh + 2);

    auto g_xxxx_xxxyy_0 = pbuffer.data(idx_eri_0_gh + 3);

    auto g_xxxx_xxxyz_0 = pbuffer.data(idx_eri_0_gh + 4);

    auto g_xxxx_xxxzz_0 = pbuffer.data(idx_eri_0_gh + 5);

    auto g_xxxx_xxyyy_0 = pbuffer.data(idx_eri_0_gh + 6);

    auto g_xxxx_xxyyz_0 = pbuffer.data(idx_eri_0_gh + 7);

    auto g_xxxx_xxyzz_0 = pbuffer.data(idx_eri_0_gh + 8);

    auto g_xxxx_xxzzz_0 = pbuffer.data(idx_eri_0_gh + 9);

    auto g_xxxx_xyyyy_0 = pbuffer.data(idx_eri_0_gh + 10);

    auto g_xxxx_xyyyz_0 = pbuffer.data(idx_eri_0_gh + 11);

    auto g_xxxx_xyyzz_0 = pbuffer.data(idx_eri_0_gh + 12);

    auto g_xxxx_xyzzz_0 = pbuffer.data(idx_eri_0_gh + 13);

    auto g_xxxx_xzzzz_0 = pbuffer.data(idx_eri_0_gh + 14);

    auto g_xxxx_yyyyy_0 = pbuffer.data(idx_eri_0_gh + 15);

    auto g_xxxx_yyyyz_0 = pbuffer.data(idx_eri_0_gh + 16);

    auto g_xxxx_yyyzz_0 = pbuffer.data(idx_eri_0_gh + 17);

    auto g_xxxx_yyzzz_0 = pbuffer.data(idx_eri_0_gh + 18);

    auto g_xxxx_yzzzz_0 = pbuffer.data(idx_eri_0_gh + 19);

    auto g_xxxx_zzzzz_0 = pbuffer.data(idx_eri_0_gh + 20);

    #pragma omp simd aligned(g_xx_xxxxx_0, g_xx_xxxxx_1, g_xx_xxxxy_0, g_xx_xxxxy_1, g_xx_xxxxz_0, g_xx_xxxxz_1, g_xx_xxxyy_0, g_xx_xxxyy_1, g_xx_xxxyz_0, g_xx_xxxyz_1, g_xx_xxxzz_0, g_xx_xxxzz_1, g_xx_xxyyy_0, g_xx_xxyyy_1, g_xx_xxyyz_0, g_xx_xxyyz_1, g_xx_xxyzz_0, g_xx_xxyzz_1, g_xx_xxzzz_0, g_xx_xxzzz_1, g_xx_xyyyy_0, g_xx_xyyyy_1, g_xx_xyyyz_0, g_xx_xyyyz_1, g_xx_xyyzz_0, g_xx_xyyzz_1, g_xx_xyzzz_0, g_xx_xyzzz_1, g_xx_xzzzz_0, g_xx_xzzzz_1, g_xx_yyyyy_0, g_xx_yyyyy_1, g_xx_yyyyz_0, g_xx_yyyyz_1, g_xx_yyyzz_0, g_xx_yyyzz_1, g_xx_yyzzz_0, g_xx_yyzzz_1, g_xx_yzzzz_0, g_xx_yzzzz_1, g_xx_zzzzz_0, g_xx_zzzzz_1, g_xxx_xxxx_1, g_xxx_xxxxx_1, g_xxx_xxxxy_1, g_xxx_xxxxz_1, g_xxx_xxxy_1, g_xxx_xxxyy_1, g_xxx_xxxyz_1, g_xxx_xxxz_1, g_xxx_xxxzz_1, g_xxx_xxyy_1, g_xxx_xxyyy_1, g_xxx_xxyyz_1, g_xxx_xxyz_1, g_xxx_xxyzz_1, g_xxx_xxzz_1, g_xxx_xxzzz_1, g_xxx_xyyy_1, g_xxx_xyyyy_1, g_xxx_xyyyz_1, g_xxx_xyyz_1, g_xxx_xyyzz_1, g_xxx_xyzz_1, g_xxx_xyzzz_1, g_xxx_xzzz_1, g_xxx_xzzzz_1, g_xxx_yyyy_1, g_xxx_yyyyy_1, g_xxx_yyyyz_1, g_xxx_yyyz_1, g_xxx_yyyzz_1, g_xxx_yyzz_1, g_xxx_yyzzz_1, g_xxx_yzzz_1, g_xxx_yzzzz_1, g_xxx_zzzz_1, g_xxx_zzzzz_1, g_xxxx_xxxxx_0, g_xxxx_xxxxy_0, g_xxxx_xxxxz_0, g_xxxx_xxxyy_0, g_xxxx_xxxyz_0, g_xxxx_xxxzz_0, g_xxxx_xxyyy_0, g_xxxx_xxyyz_0, g_xxxx_xxyzz_0, g_xxxx_xxzzz_0, g_xxxx_xyyyy_0, g_xxxx_xyyyz_0, g_xxxx_xyyzz_0, g_xxxx_xyzzz_0, g_xxxx_xzzzz_0, g_xxxx_yyyyy_0, g_xxxx_yyyyz_0, g_xxxx_yyyzz_0, g_xxxx_yyzzz_0, g_xxxx_yzzzz_0, g_xxxx_zzzzz_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxx_xxxxx_0[i] = 3.0 * g_xx_xxxxx_0[i] * fbe_0 - 3.0 * g_xx_xxxxx_1[i] * fz_be_0 + 5.0 * g_xxx_xxxx_1[i] * fe_0 + g_xxx_xxxxx_1[i] * pa_x[i];

        g_xxxx_xxxxy_0[i] = 3.0 * g_xx_xxxxy_0[i] * fbe_0 - 3.0 * g_xx_xxxxy_1[i] * fz_be_0 + 4.0 * g_xxx_xxxy_1[i] * fe_0 + g_xxx_xxxxy_1[i] * pa_x[i];

        g_xxxx_xxxxz_0[i] = 3.0 * g_xx_xxxxz_0[i] * fbe_0 - 3.0 * g_xx_xxxxz_1[i] * fz_be_0 + 4.0 * g_xxx_xxxz_1[i] * fe_0 + g_xxx_xxxxz_1[i] * pa_x[i];

        g_xxxx_xxxyy_0[i] = 3.0 * g_xx_xxxyy_0[i] * fbe_0 - 3.0 * g_xx_xxxyy_1[i] * fz_be_0 + 3.0 * g_xxx_xxyy_1[i] * fe_0 + g_xxx_xxxyy_1[i] * pa_x[i];

        g_xxxx_xxxyz_0[i] = 3.0 * g_xx_xxxyz_0[i] * fbe_0 - 3.0 * g_xx_xxxyz_1[i] * fz_be_0 + 3.0 * g_xxx_xxyz_1[i] * fe_0 + g_xxx_xxxyz_1[i] * pa_x[i];

        g_xxxx_xxxzz_0[i] = 3.0 * g_xx_xxxzz_0[i] * fbe_0 - 3.0 * g_xx_xxxzz_1[i] * fz_be_0 + 3.0 * g_xxx_xxzz_1[i] * fe_0 + g_xxx_xxxzz_1[i] * pa_x[i];

        g_xxxx_xxyyy_0[i] = 3.0 * g_xx_xxyyy_0[i] * fbe_0 - 3.0 * g_xx_xxyyy_1[i] * fz_be_0 + 2.0 * g_xxx_xyyy_1[i] * fe_0 + g_xxx_xxyyy_1[i] * pa_x[i];

        g_xxxx_xxyyz_0[i] = 3.0 * g_xx_xxyyz_0[i] * fbe_0 - 3.0 * g_xx_xxyyz_1[i] * fz_be_0 + 2.0 * g_xxx_xyyz_1[i] * fe_0 + g_xxx_xxyyz_1[i] * pa_x[i];

        g_xxxx_xxyzz_0[i] = 3.0 * g_xx_xxyzz_0[i] * fbe_0 - 3.0 * g_xx_xxyzz_1[i] * fz_be_0 + 2.0 * g_xxx_xyzz_1[i] * fe_0 + g_xxx_xxyzz_1[i] * pa_x[i];

        g_xxxx_xxzzz_0[i] = 3.0 * g_xx_xxzzz_0[i] * fbe_0 - 3.0 * g_xx_xxzzz_1[i] * fz_be_0 + 2.0 * g_xxx_xzzz_1[i] * fe_0 + g_xxx_xxzzz_1[i] * pa_x[i];

        g_xxxx_xyyyy_0[i] = 3.0 * g_xx_xyyyy_0[i] * fbe_0 - 3.0 * g_xx_xyyyy_1[i] * fz_be_0 + g_xxx_yyyy_1[i] * fe_0 + g_xxx_xyyyy_1[i] * pa_x[i];

        g_xxxx_xyyyz_0[i] = 3.0 * g_xx_xyyyz_0[i] * fbe_0 - 3.0 * g_xx_xyyyz_1[i] * fz_be_0 + g_xxx_yyyz_1[i] * fe_0 + g_xxx_xyyyz_1[i] * pa_x[i];

        g_xxxx_xyyzz_0[i] = 3.0 * g_xx_xyyzz_0[i] * fbe_0 - 3.0 * g_xx_xyyzz_1[i] * fz_be_0 + g_xxx_yyzz_1[i] * fe_0 + g_xxx_xyyzz_1[i] * pa_x[i];

        g_xxxx_xyzzz_0[i] = 3.0 * g_xx_xyzzz_0[i] * fbe_0 - 3.0 * g_xx_xyzzz_1[i] * fz_be_0 + g_xxx_yzzz_1[i] * fe_0 + g_xxx_xyzzz_1[i] * pa_x[i];

        g_xxxx_xzzzz_0[i] = 3.0 * g_xx_xzzzz_0[i] * fbe_0 - 3.0 * g_xx_xzzzz_1[i] * fz_be_0 + g_xxx_zzzz_1[i] * fe_0 + g_xxx_xzzzz_1[i] * pa_x[i];

        g_xxxx_yyyyy_0[i] = 3.0 * g_xx_yyyyy_0[i] * fbe_0 - 3.0 * g_xx_yyyyy_1[i] * fz_be_0 + g_xxx_yyyyy_1[i] * pa_x[i];

        g_xxxx_yyyyz_0[i] = 3.0 * g_xx_yyyyz_0[i] * fbe_0 - 3.0 * g_xx_yyyyz_1[i] * fz_be_0 + g_xxx_yyyyz_1[i] * pa_x[i];

        g_xxxx_yyyzz_0[i] = 3.0 * g_xx_yyyzz_0[i] * fbe_0 - 3.0 * g_xx_yyyzz_1[i] * fz_be_0 + g_xxx_yyyzz_1[i] * pa_x[i];

        g_xxxx_yyzzz_0[i] = 3.0 * g_xx_yyzzz_0[i] * fbe_0 - 3.0 * g_xx_yyzzz_1[i] * fz_be_0 + g_xxx_yyzzz_1[i] * pa_x[i];

        g_xxxx_yzzzz_0[i] = 3.0 * g_xx_yzzzz_0[i] * fbe_0 - 3.0 * g_xx_yzzzz_1[i] * fz_be_0 + g_xxx_yzzzz_1[i] * pa_x[i];

        g_xxxx_zzzzz_0[i] = 3.0 * g_xx_zzzzz_0[i] * fbe_0 - 3.0 * g_xx_zzzzz_1[i] * fz_be_0 + g_xxx_zzzzz_1[i] * pa_x[i];
    }

    // Set up 21-42 components of targeted buffer : GH

    auto g_xxxy_xxxxx_0 = pbuffer.data(idx_eri_0_gh + 21);

    auto g_xxxy_xxxxy_0 = pbuffer.data(idx_eri_0_gh + 22);

    auto g_xxxy_xxxxz_0 = pbuffer.data(idx_eri_0_gh + 23);

    auto g_xxxy_xxxyy_0 = pbuffer.data(idx_eri_0_gh + 24);

    auto g_xxxy_xxxyz_0 = pbuffer.data(idx_eri_0_gh + 25);

    auto g_xxxy_xxxzz_0 = pbuffer.data(idx_eri_0_gh + 26);

    auto g_xxxy_xxyyy_0 = pbuffer.data(idx_eri_0_gh + 27);

    auto g_xxxy_xxyyz_0 = pbuffer.data(idx_eri_0_gh + 28);

    auto g_xxxy_xxyzz_0 = pbuffer.data(idx_eri_0_gh + 29);

    auto g_xxxy_xxzzz_0 = pbuffer.data(idx_eri_0_gh + 30);

    auto g_xxxy_xyyyy_0 = pbuffer.data(idx_eri_0_gh + 31);

    auto g_xxxy_xyyyz_0 = pbuffer.data(idx_eri_0_gh + 32);

    auto g_xxxy_xyyzz_0 = pbuffer.data(idx_eri_0_gh + 33);

    auto g_xxxy_xyzzz_0 = pbuffer.data(idx_eri_0_gh + 34);

    auto g_xxxy_xzzzz_0 = pbuffer.data(idx_eri_0_gh + 35);

    auto g_xxxy_yyyyy_0 = pbuffer.data(idx_eri_0_gh + 36);

    auto g_xxxy_yyyyz_0 = pbuffer.data(idx_eri_0_gh + 37);

    auto g_xxxy_yyyzz_0 = pbuffer.data(idx_eri_0_gh + 38);

    auto g_xxxy_yyzzz_0 = pbuffer.data(idx_eri_0_gh + 39);

    auto g_xxxy_yzzzz_0 = pbuffer.data(idx_eri_0_gh + 40);

    auto g_xxxy_zzzzz_0 = pbuffer.data(idx_eri_0_gh + 41);

    #pragma omp simd aligned(g_xxx_xxxx_1, g_xxx_xxxxx_1, g_xxx_xxxxy_1, g_xxx_xxxxz_1, g_xxx_xxxy_1, g_xxx_xxxyy_1, g_xxx_xxxyz_1, g_xxx_xxxz_1, g_xxx_xxxzz_1, g_xxx_xxyy_1, g_xxx_xxyyy_1, g_xxx_xxyyz_1, g_xxx_xxyz_1, g_xxx_xxyzz_1, g_xxx_xxzz_1, g_xxx_xxzzz_1, g_xxx_xyyy_1, g_xxx_xyyyy_1, g_xxx_xyyyz_1, g_xxx_xyyz_1, g_xxx_xyyzz_1, g_xxx_xyzz_1, g_xxx_xyzzz_1, g_xxx_xzzz_1, g_xxx_xzzzz_1, g_xxx_yyyy_1, g_xxx_yyyyy_1, g_xxx_yyyyz_1, g_xxx_yyyz_1, g_xxx_yyyzz_1, g_xxx_yyzz_1, g_xxx_yyzzz_1, g_xxx_yzzz_1, g_xxx_yzzzz_1, g_xxx_zzzz_1, g_xxx_zzzzz_1, g_xxxy_xxxxx_0, g_xxxy_xxxxy_0, g_xxxy_xxxxz_0, g_xxxy_xxxyy_0, g_xxxy_xxxyz_0, g_xxxy_xxxzz_0, g_xxxy_xxyyy_0, g_xxxy_xxyyz_0, g_xxxy_xxyzz_0, g_xxxy_xxzzz_0, g_xxxy_xyyyy_0, g_xxxy_xyyyz_0, g_xxxy_xyyzz_0, g_xxxy_xyzzz_0, g_xxxy_xzzzz_0, g_xxxy_yyyyy_0, g_xxxy_yyyyz_0, g_xxxy_yyyzz_0, g_xxxy_yyzzz_0, g_xxxy_yzzzz_0, g_xxxy_zzzzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxy_xxxxx_0[i] = g_xxx_xxxxx_1[i] * pa_y[i];

        g_xxxy_xxxxy_0[i] = g_xxx_xxxx_1[i] * fe_0 + g_xxx_xxxxy_1[i] * pa_y[i];

        g_xxxy_xxxxz_0[i] = g_xxx_xxxxz_1[i] * pa_y[i];

        g_xxxy_xxxyy_0[i] = 2.0 * g_xxx_xxxy_1[i] * fe_0 + g_xxx_xxxyy_1[i] * pa_y[i];

        g_xxxy_xxxyz_0[i] = g_xxx_xxxz_1[i] * fe_0 + g_xxx_xxxyz_1[i] * pa_y[i];

        g_xxxy_xxxzz_0[i] = g_xxx_xxxzz_1[i] * pa_y[i];

        g_xxxy_xxyyy_0[i] = 3.0 * g_xxx_xxyy_1[i] * fe_0 + g_xxx_xxyyy_1[i] * pa_y[i];

        g_xxxy_xxyyz_0[i] = 2.0 * g_xxx_xxyz_1[i] * fe_0 + g_xxx_xxyyz_1[i] * pa_y[i];

        g_xxxy_xxyzz_0[i] = g_xxx_xxzz_1[i] * fe_0 + g_xxx_xxyzz_1[i] * pa_y[i];

        g_xxxy_xxzzz_0[i] = g_xxx_xxzzz_1[i] * pa_y[i];

        g_xxxy_xyyyy_0[i] = 4.0 * g_xxx_xyyy_1[i] * fe_0 + g_xxx_xyyyy_1[i] * pa_y[i];

        g_xxxy_xyyyz_0[i] = 3.0 * g_xxx_xyyz_1[i] * fe_0 + g_xxx_xyyyz_1[i] * pa_y[i];

        g_xxxy_xyyzz_0[i] = 2.0 * g_xxx_xyzz_1[i] * fe_0 + g_xxx_xyyzz_1[i] * pa_y[i];

        g_xxxy_xyzzz_0[i] = g_xxx_xzzz_1[i] * fe_0 + g_xxx_xyzzz_1[i] * pa_y[i];

        g_xxxy_xzzzz_0[i] = g_xxx_xzzzz_1[i] * pa_y[i];

        g_xxxy_yyyyy_0[i] = 5.0 * g_xxx_yyyy_1[i] * fe_0 + g_xxx_yyyyy_1[i] * pa_y[i];

        g_xxxy_yyyyz_0[i] = 4.0 * g_xxx_yyyz_1[i] * fe_0 + g_xxx_yyyyz_1[i] * pa_y[i];

        g_xxxy_yyyzz_0[i] = 3.0 * g_xxx_yyzz_1[i] * fe_0 + g_xxx_yyyzz_1[i] * pa_y[i];

        g_xxxy_yyzzz_0[i] = 2.0 * g_xxx_yzzz_1[i] * fe_0 + g_xxx_yyzzz_1[i] * pa_y[i];

        g_xxxy_yzzzz_0[i] = g_xxx_zzzz_1[i] * fe_0 + g_xxx_yzzzz_1[i] * pa_y[i];

        g_xxxy_zzzzz_0[i] = g_xxx_zzzzz_1[i] * pa_y[i];
    }

    // Set up 42-63 components of targeted buffer : GH

    auto g_xxxz_xxxxx_0 = pbuffer.data(idx_eri_0_gh + 42);

    auto g_xxxz_xxxxy_0 = pbuffer.data(idx_eri_0_gh + 43);

    auto g_xxxz_xxxxz_0 = pbuffer.data(idx_eri_0_gh + 44);

    auto g_xxxz_xxxyy_0 = pbuffer.data(idx_eri_0_gh + 45);

    auto g_xxxz_xxxyz_0 = pbuffer.data(idx_eri_0_gh + 46);

    auto g_xxxz_xxxzz_0 = pbuffer.data(idx_eri_0_gh + 47);

    auto g_xxxz_xxyyy_0 = pbuffer.data(idx_eri_0_gh + 48);

    auto g_xxxz_xxyyz_0 = pbuffer.data(idx_eri_0_gh + 49);

    auto g_xxxz_xxyzz_0 = pbuffer.data(idx_eri_0_gh + 50);

    auto g_xxxz_xxzzz_0 = pbuffer.data(idx_eri_0_gh + 51);

    auto g_xxxz_xyyyy_0 = pbuffer.data(idx_eri_0_gh + 52);

    auto g_xxxz_xyyyz_0 = pbuffer.data(idx_eri_0_gh + 53);

    auto g_xxxz_xyyzz_0 = pbuffer.data(idx_eri_0_gh + 54);

    auto g_xxxz_xyzzz_0 = pbuffer.data(idx_eri_0_gh + 55);

    auto g_xxxz_xzzzz_0 = pbuffer.data(idx_eri_0_gh + 56);

    auto g_xxxz_yyyyy_0 = pbuffer.data(idx_eri_0_gh + 57);

    auto g_xxxz_yyyyz_0 = pbuffer.data(idx_eri_0_gh + 58);

    auto g_xxxz_yyyzz_0 = pbuffer.data(idx_eri_0_gh + 59);

    auto g_xxxz_yyzzz_0 = pbuffer.data(idx_eri_0_gh + 60);

    auto g_xxxz_yzzzz_0 = pbuffer.data(idx_eri_0_gh + 61);

    auto g_xxxz_zzzzz_0 = pbuffer.data(idx_eri_0_gh + 62);

    #pragma omp simd aligned(g_xxx_xxxx_1, g_xxx_xxxxx_1, g_xxx_xxxxy_1, g_xxx_xxxxz_1, g_xxx_xxxy_1, g_xxx_xxxyy_1, g_xxx_xxxyz_1, g_xxx_xxxz_1, g_xxx_xxxzz_1, g_xxx_xxyy_1, g_xxx_xxyyy_1, g_xxx_xxyyz_1, g_xxx_xxyz_1, g_xxx_xxyzz_1, g_xxx_xxzz_1, g_xxx_xxzzz_1, g_xxx_xyyy_1, g_xxx_xyyyy_1, g_xxx_xyyyz_1, g_xxx_xyyz_1, g_xxx_xyyzz_1, g_xxx_xyzz_1, g_xxx_xyzzz_1, g_xxx_xzzz_1, g_xxx_xzzzz_1, g_xxx_yyyy_1, g_xxx_yyyyy_1, g_xxx_yyyyz_1, g_xxx_yyyz_1, g_xxx_yyyzz_1, g_xxx_yyzz_1, g_xxx_yyzzz_1, g_xxx_yzzz_1, g_xxx_yzzzz_1, g_xxx_zzzz_1, g_xxx_zzzzz_1, g_xxxz_xxxxx_0, g_xxxz_xxxxy_0, g_xxxz_xxxxz_0, g_xxxz_xxxyy_0, g_xxxz_xxxyz_0, g_xxxz_xxxzz_0, g_xxxz_xxyyy_0, g_xxxz_xxyyz_0, g_xxxz_xxyzz_0, g_xxxz_xxzzz_0, g_xxxz_xyyyy_0, g_xxxz_xyyyz_0, g_xxxz_xyyzz_0, g_xxxz_xyzzz_0, g_xxxz_xzzzz_0, g_xxxz_yyyyy_0, g_xxxz_yyyyz_0, g_xxxz_yyyzz_0, g_xxxz_yyzzz_0, g_xxxz_yzzzz_0, g_xxxz_zzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxz_xxxxx_0[i] = g_xxx_xxxxx_1[i] * pa_z[i];

        g_xxxz_xxxxy_0[i] = g_xxx_xxxxy_1[i] * pa_z[i];

        g_xxxz_xxxxz_0[i] = g_xxx_xxxx_1[i] * fe_0 + g_xxx_xxxxz_1[i] * pa_z[i];

        g_xxxz_xxxyy_0[i] = g_xxx_xxxyy_1[i] * pa_z[i];

        g_xxxz_xxxyz_0[i] = g_xxx_xxxy_1[i] * fe_0 + g_xxx_xxxyz_1[i] * pa_z[i];

        g_xxxz_xxxzz_0[i] = 2.0 * g_xxx_xxxz_1[i] * fe_0 + g_xxx_xxxzz_1[i] * pa_z[i];

        g_xxxz_xxyyy_0[i] = g_xxx_xxyyy_1[i] * pa_z[i];

        g_xxxz_xxyyz_0[i] = g_xxx_xxyy_1[i] * fe_0 + g_xxx_xxyyz_1[i] * pa_z[i];

        g_xxxz_xxyzz_0[i] = 2.0 * g_xxx_xxyz_1[i] * fe_0 + g_xxx_xxyzz_1[i] * pa_z[i];

        g_xxxz_xxzzz_0[i] = 3.0 * g_xxx_xxzz_1[i] * fe_0 + g_xxx_xxzzz_1[i] * pa_z[i];

        g_xxxz_xyyyy_0[i] = g_xxx_xyyyy_1[i] * pa_z[i];

        g_xxxz_xyyyz_0[i] = g_xxx_xyyy_1[i] * fe_0 + g_xxx_xyyyz_1[i] * pa_z[i];

        g_xxxz_xyyzz_0[i] = 2.0 * g_xxx_xyyz_1[i] * fe_0 + g_xxx_xyyzz_1[i] * pa_z[i];

        g_xxxz_xyzzz_0[i] = 3.0 * g_xxx_xyzz_1[i] * fe_0 + g_xxx_xyzzz_1[i] * pa_z[i];

        g_xxxz_xzzzz_0[i] = 4.0 * g_xxx_xzzz_1[i] * fe_0 + g_xxx_xzzzz_1[i] * pa_z[i];

        g_xxxz_yyyyy_0[i] = g_xxx_yyyyy_1[i] * pa_z[i];

        g_xxxz_yyyyz_0[i] = g_xxx_yyyy_1[i] * fe_0 + g_xxx_yyyyz_1[i] * pa_z[i];

        g_xxxz_yyyzz_0[i] = 2.0 * g_xxx_yyyz_1[i] * fe_0 + g_xxx_yyyzz_1[i] * pa_z[i];

        g_xxxz_yyzzz_0[i] = 3.0 * g_xxx_yyzz_1[i] * fe_0 + g_xxx_yyzzz_1[i] * pa_z[i];

        g_xxxz_yzzzz_0[i] = 4.0 * g_xxx_yzzz_1[i] * fe_0 + g_xxx_yzzzz_1[i] * pa_z[i];

        g_xxxz_zzzzz_0[i] = 5.0 * g_xxx_zzzz_1[i] * fe_0 + g_xxx_zzzzz_1[i] * pa_z[i];
    }

    // Set up 63-84 components of targeted buffer : GH

    auto g_xxyy_xxxxx_0 = pbuffer.data(idx_eri_0_gh + 63);

    auto g_xxyy_xxxxy_0 = pbuffer.data(idx_eri_0_gh + 64);

    auto g_xxyy_xxxxz_0 = pbuffer.data(idx_eri_0_gh + 65);

    auto g_xxyy_xxxyy_0 = pbuffer.data(idx_eri_0_gh + 66);

    auto g_xxyy_xxxyz_0 = pbuffer.data(idx_eri_0_gh + 67);

    auto g_xxyy_xxxzz_0 = pbuffer.data(idx_eri_0_gh + 68);

    auto g_xxyy_xxyyy_0 = pbuffer.data(idx_eri_0_gh + 69);

    auto g_xxyy_xxyyz_0 = pbuffer.data(idx_eri_0_gh + 70);

    auto g_xxyy_xxyzz_0 = pbuffer.data(idx_eri_0_gh + 71);

    auto g_xxyy_xxzzz_0 = pbuffer.data(idx_eri_0_gh + 72);

    auto g_xxyy_xyyyy_0 = pbuffer.data(idx_eri_0_gh + 73);

    auto g_xxyy_xyyyz_0 = pbuffer.data(idx_eri_0_gh + 74);

    auto g_xxyy_xyyzz_0 = pbuffer.data(idx_eri_0_gh + 75);

    auto g_xxyy_xyzzz_0 = pbuffer.data(idx_eri_0_gh + 76);

    auto g_xxyy_xzzzz_0 = pbuffer.data(idx_eri_0_gh + 77);

    auto g_xxyy_yyyyy_0 = pbuffer.data(idx_eri_0_gh + 78);

    auto g_xxyy_yyyyz_0 = pbuffer.data(idx_eri_0_gh + 79);

    auto g_xxyy_yyyzz_0 = pbuffer.data(idx_eri_0_gh + 80);

    auto g_xxyy_yyzzz_0 = pbuffer.data(idx_eri_0_gh + 81);

    auto g_xxyy_yzzzz_0 = pbuffer.data(idx_eri_0_gh + 82);

    auto g_xxyy_zzzzz_0 = pbuffer.data(idx_eri_0_gh + 83);

    #pragma omp simd aligned(g_xx_xxxxx_0, g_xx_xxxxx_1, g_xx_xxxxz_0, g_xx_xxxxz_1, g_xx_xxxzz_0, g_xx_xxxzz_1, g_xx_xxzzz_0, g_xx_xxzzz_1, g_xx_xzzzz_0, g_xx_xzzzz_1, g_xxy_xxxxx_1, g_xxy_xxxxz_1, g_xxy_xxxzz_1, g_xxy_xxzzz_1, g_xxy_xzzzz_1, g_xxyy_xxxxx_0, g_xxyy_xxxxy_0, g_xxyy_xxxxz_0, g_xxyy_xxxyy_0, g_xxyy_xxxyz_0, g_xxyy_xxxzz_0, g_xxyy_xxyyy_0, g_xxyy_xxyyz_0, g_xxyy_xxyzz_0, g_xxyy_xxzzz_0, g_xxyy_xyyyy_0, g_xxyy_xyyyz_0, g_xxyy_xyyzz_0, g_xxyy_xyzzz_0, g_xxyy_xzzzz_0, g_xxyy_yyyyy_0, g_xxyy_yyyyz_0, g_xxyy_yyyzz_0, g_xxyy_yyzzz_0, g_xxyy_yzzzz_0, g_xxyy_zzzzz_0, g_xyy_xxxxy_1, g_xyy_xxxy_1, g_xyy_xxxyy_1, g_xyy_xxxyz_1, g_xyy_xxyy_1, g_xyy_xxyyy_1, g_xyy_xxyyz_1, g_xyy_xxyz_1, g_xyy_xxyzz_1, g_xyy_xyyy_1, g_xyy_xyyyy_1, g_xyy_xyyyz_1, g_xyy_xyyz_1, g_xyy_xyyzz_1, g_xyy_xyzz_1, g_xyy_xyzzz_1, g_xyy_yyyy_1, g_xyy_yyyyy_1, g_xyy_yyyyz_1, g_xyy_yyyz_1, g_xyy_yyyzz_1, g_xyy_yyzz_1, g_xyy_yyzzz_1, g_xyy_yzzz_1, g_xyy_yzzzz_1, g_xyy_zzzzz_1, g_yy_xxxxy_0, g_yy_xxxxy_1, g_yy_xxxyy_0, g_yy_xxxyy_1, g_yy_xxxyz_0, g_yy_xxxyz_1, g_yy_xxyyy_0, g_yy_xxyyy_1, g_yy_xxyyz_0, g_yy_xxyyz_1, g_yy_xxyzz_0, g_yy_xxyzz_1, g_yy_xyyyy_0, g_yy_xyyyy_1, g_yy_xyyyz_0, g_yy_xyyyz_1, g_yy_xyyzz_0, g_yy_xyyzz_1, g_yy_xyzzz_0, g_yy_xyzzz_1, g_yy_yyyyy_0, g_yy_yyyyy_1, g_yy_yyyyz_0, g_yy_yyyyz_1, g_yy_yyyzz_0, g_yy_yyyzz_1, g_yy_yyzzz_0, g_yy_yyzzz_1, g_yy_yzzzz_0, g_yy_yzzzz_1, g_yy_zzzzz_0, g_yy_zzzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxyy_xxxxx_0[i] = g_xx_xxxxx_0[i] * fbe_0 - g_xx_xxxxx_1[i] * fz_be_0 + g_xxy_xxxxx_1[i] * pa_y[i];

        g_xxyy_xxxxy_0[i] = g_yy_xxxxy_0[i] * fbe_0 - g_yy_xxxxy_1[i] * fz_be_0 + 4.0 * g_xyy_xxxy_1[i] * fe_0 + g_xyy_xxxxy_1[i] * pa_x[i];

        g_xxyy_xxxxz_0[i] = g_xx_xxxxz_0[i] * fbe_0 - g_xx_xxxxz_1[i] * fz_be_0 + g_xxy_xxxxz_1[i] * pa_y[i];

        g_xxyy_xxxyy_0[i] = g_yy_xxxyy_0[i] * fbe_0 - g_yy_xxxyy_1[i] * fz_be_0 + 3.0 * g_xyy_xxyy_1[i] * fe_0 + g_xyy_xxxyy_1[i] * pa_x[i];

        g_xxyy_xxxyz_0[i] = g_yy_xxxyz_0[i] * fbe_0 - g_yy_xxxyz_1[i] * fz_be_0 + 3.0 * g_xyy_xxyz_1[i] * fe_0 + g_xyy_xxxyz_1[i] * pa_x[i];

        g_xxyy_xxxzz_0[i] = g_xx_xxxzz_0[i] * fbe_0 - g_xx_xxxzz_1[i] * fz_be_0 + g_xxy_xxxzz_1[i] * pa_y[i];

        g_xxyy_xxyyy_0[i] = g_yy_xxyyy_0[i] * fbe_0 - g_yy_xxyyy_1[i] * fz_be_0 + 2.0 * g_xyy_xyyy_1[i] * fe_0 + g_xyy_xxyyy_1[i] * pa_x[i];

        g_xxyy_xxyyz_0[i] = g_yy_xxyyz_0[i] * fbe_0 - g_yy_xxyyz_1[i] * fz_be_0 + 2.0 * g_xyy_xyyz_1[i] * fe_0 + g_xyy_xxyyz_1[i] * pa_x[i];

        g_xxyy_xxyzz_0[i] = g_yy_xxyzz_0[i] * fbe_0 - g_yy_xxyzz_1[i] * fz_be_0 + 2.0 * g_xyy_xyzz_1[i] * fe_0 + g_xyy_xxyzz_1[i] * pa_x[i];

        g_xxyy_xxzzz_0[i] = g_xx_xxzzz_0[i] * fbe_0 - g_xx_xxzzz_1[i] * fz_be_0 + g_xxy_xxzzz_1[i] * pa_y[i];

        g_xxyy_xyyyy_0[i] = g_yy_xyyyy_0[i] * fbe_0 - g_yy_xyyyy_1[i] * fz_be_0 + g_xyy_yyyy_1[i] * fe_0 + g_xyy_xyyyy_1[i] * pa_x[i];

        g_xxyy_xyyyz_0[i] = g_yy_xyyyz_0[i] * fbe_0 - g_yy_xyyyz_1[i] * fz_be_0 + g_xyy_yyyz_1[i] * fe_0 + g_xyy_xyyyz_1[i] * pa_x[i];

        g_xxyy_xyyzz_0[i] = g_yy_xyyzz_0[i] * fbe_0 - g_yy_xyyzz_1[i] * fz_be_0 + g_xyy_yyzz_1[i] * fe_0 + g_xyy_xyyzz_1[i] * pa_x[i];

        g_xxyy_xyzzz_0[i] = g_yy_xyzzz_0[i] * fbe_0 - g_yy_xyzzz_1[i] * fz_be_0 + g_xyy_yzzz_1[i] * fe_0 + g_xyy_xyzzz_1[i] * pa_x[i];

        g_xxyy_xzzzz_0[i] = g_xx_xzzzz_0[i] * fbe_0 - g_xx_xzzzz_1[i] * fz_be_0 + g_xxy_xzzzz_1[i] * pa_y[i];

        g_xxyy_yyyyy_0[i] = g_yy_yyyyy_0[i] * fbe_0 - g_yy_yyyyy_1[i] * fz_be_0 + g_xyy_yyyyy_1[i] * pa_x[i];

        g_xxyy_yyyyz_0[i] = g_yy_yyyyz_0[i] * fbe_0 - g_yy_yyyyz_1[i] * fz_be_0 + g_xyy_yyyyz_1[i] * pa_x[i];

        g_xxyy_yyyzz_0[i] = g_yy_yyyzz_0[i] * fbe_0 - g_yy_yyyzz_1[i] * fz_be_0 + g_xyy_yyyzz_1[i] * pa_x[i];

        g_xxyy_yyzzz_0[i] = g_yy_yyzzz_0[i] * fbe_0 - g_yy_yyzzz_1[i] * fz_be_0 + g_xyy_yyzzz_1[i] * pa_x[i];

        g_xxyy_yzzzz_0[i] = g_yy_yzzzz_0[i] * fbe_0 - g_yy_yzzzz_1[i] * fz_be_0 + g_xyy_yzzzz_1[i] * pa_x[i];

        g_xxyy_zzzzz_0[i] = g_yy_zzzzz_0[i] * fbe_0 - g_yy_zzzzz_1[i] * fz_be_0 + g_xyy_zzzzz_1[i] * pa_x[i];
    }

    // Set up 84-105 components of targeted buffer : GH

    auto g_xxyz_xxxxx_0 = pbuffer.data(idx_eri_0_gh + 84);

    auto g_xxyz_xxxxy_0 = pbuffer.data(idx_eri_0_gh + 85);

    auto g_xxyz_xxxxz_0 = pbuffer.data(idx_eri_0_gh + 86);

    auto g_xxyz_xxxyy_0 = pbuffer.data(idx_eri_0_gh + 87);

    auto g_xxyz_xxxyz_0 = pbuffer.data(idx_eri_0_gh + 88);

    auto g_xxyz_xxxzz_0 = pbuffer.data(idx_eri_0_gh + 89);

    auto g_xxyz_xxyyy_0 = pbuffer.data(idx_eri_0_gh + 90);

    auto g_xxyz_xxyyz_0 = pbuffer.data(idx_eri_0_gh + 91);

    auto g_xxyz_xxyzz_0 = pbuffer.data(idx_eri_0_gh + 92);

    auto g_xxyz_xxzzz_0 = pbuffer.data(idx_eri_0_gh + 93);

    auto g_xxyz_xyyyy_0 = pbuffer.data(idx_eri_0_gh + 94);

    auto g_xxyz_xyyyz_0 = pbuffer.data(idx_eri_0_gh + 95);

    auto g_xxyz_xyyzz_0 = pbuffer.data(idx_eri_0_gh + 96);

    auto g_xxyz_xyzzz_0 = pbuffer.data(idx_eri_0_gh + 97);

    auto g_xxyz_xzzzz_0 = pbuffer.data(idx_eri_0_gh + 98);

    auto g_xxyz_yyyyy_0 = pbuffer.data(idx_eri_0_gh + 99);

    auto g_xxyz_yyyyz_0 = pbuffer.data(idx_eri_0_gh + 100);

    auto g_xxyz_yyyzz_0 = pbuffer.data(idx_eri_0_gh + 101);

    auto g_xxyz_yyzzz_0 = pbuffer.data(idx_eri_0_gh + 102);

    auto g_xxyz_yzzzz_0 = pbuffer.data(idx_eri_0_gh + 103);

    auto g_xxyz_zzzzz_0 = pbuffer.data(idx_eri_0_gh + 104);

    #pragma omp simd aligned(g_xxy_xxxxy_1, g_xxy_xxxyy_1, g_xxy_xxyyy_1, g_xxy_xyyyy_1, g_xxy_yyyyy_1, g_xxyz_xxxxx_0, g_xxyz_xxxxy_0, g_xxyz_xxxxz_0, g_xxyz_xxxyy_0, g_xxyz_xxxyz_0, g_xxyz_xxxzz_0, g_xxyz_xxyyy_0, g_xxyz_xxyyz_0, g_xxyz_xxyzz_0, g_xxyz_xxzzz_0, g_xxyz_xyyyy_0, g_xxyz_xyyyz_0, g_xxyz_xyyzz_0, g_xxyz_xyzzz_0, g_xxyz_xzzzz_0, g_xxyz_yyyyy_0, g_xxyz_yyyyz_0, g_xxyz_yyyzz_0, g_xxyz_yyzzz_0, g_xxyz_yzzzz_0, g_xxyz_zzzzz_0, g_xxz_xxxxx_1, g_xxz_xxxxz_1, g_xxz_xxxyz_1, g_xxz_xxxz_1, g_xxz_xxxzz_1, g_xxz_xxyyz_1, g_xxz_xxyz_1, g_xxz_xxyzz_1, g_xxz_xxzz_1, g_xxz_xxzzz_1, g_xxz_xyyyz_1, g_xxz_xyyz_1, g_xxz_xyyzz_1, g_xxz_xyzz_1, g_xxz_xyzzz_1, g_xxz_xzzz_1, g_xxz_xzzzz_1, g_xxz_yyyyz_1, g_xxz_yyyz_1, g_xxz_yyyzz_1, g_xxz_yyzz_1, g_xxz_yyzzz_1, g_xxz_yzzz_1, g_xxz_yzzzz_1, g_xxz_zzzz_1, g_xxz_zzzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyz_xxxxx_0[i] = g_xxz_xxxxx_1[i] * pa_y[i];

        g_xxyz_xxxxy_0[i] = g_xxy_xxxxy_1[i] * pa_z[i];

        g_xxyz_xxxxz_0[i] = g_xxz_xxxxz_1[i] * pa_y[i];

        g_xxyz_xxxyy_0[i] = g_xxy_xxxyy_1[i] * pa_z[i];

        g_xxyz_xxxyz_0[i] = g_xxz_xxxz_1[i] * fe_0 + g_xxz_xxxyz_1[i] * pa_y[i];

        g_xxyz_xxxzz_0[i] = g_xxz_xxxzz_1[i] * pa_y[i];

        g_xxyz_xxyyy_0[i] = g_xxy_xxyyy_1[i] * pa_z[i];

        g_xxyz_xxyyz_0[i] = 2.0 * g_xxz_xxyz_1[i] * fe_0 + g_xxz_xxyyz_1[i] * pa_y[i];

        g_xxyz_xxyzz_0[i] = g_xxz_xxzz_1[i] * fe_0 + g_xxz_xxyzz_1[i] * pa_y[i];

        g_xxyz_xxzzz_0[i] = g_xxz_xxzzz_1[i] * pa_y[i];

        g_xxyz_xyyyy_0[i] = g_xxy_xyyyy_1[i] * pa_z[i];

        g_xxyz_xyyyz_0[i] = 3.0 * g_xxz_xyyz_1[i] * fe_0 + g_xxz_xyyyz_1[i] * pa_y[i];

        g_xxyz_xyyzz_0[i] = 2.0 * g_xxz_xyzz_1[i] * fe_0 + g_xxz_xyyzz_1[i] * pa_y[i];

        g_xxyz_xyzzz_0[i] = g_xxz_xzzz_1[i] * fe_0 + g_xxz_xyzzz_1[i] * pa_y[i];

        g_xxyz_xzzzz_0[i] = g_xxz_xzzzz_1[i] * pa_y[i];

        g_xxyz_yyyyy_0[i] = g_xxy_yyyyy_1[i] * pa_z[i];

        g_xxyz_yyyyz_0[i] = 4.0 * g_xxz_yyyz_1[i] * fe_0 + g_xxz_yyyyz_1[i] * pa_y[i];

        g_xxyz_yyyzz_0[i] = 3.0 * g_xxz_yyzz_1[i] * fe_0 + g_xxz_yyyzz_1[i] * pa_y[i];

        g_xxyz_yyzzz_0[i] = 2.0 * g_xxz_yzzz_1[i] * fe_0 + g_xxz_yyzzz_1[i] * pa_y[i];

        g_xxyz_yzzzz_0[i] = g_xxz_zzzz_1[i] * fe_0 + g_xxz_yzzzz_1[i] * pa_y[i];

        g_xxyz_zzzzz_0[i] = g_xxz_zzzzz_1[i] * pa_y[i];
    }

    // Set up 105-126 components of targeted buffer : GH

    auto g_xxzz_xxxxx_0 = pbuffer.data(idx_eri_0_gh + 105);

    auto g_xxzz_xxxxy_0 = pbuffer.data(idx_eri_0_gh + 106);

    auto g_xxzz_xxxxz_0 = pbuffer.data(idx_eri_0_gh + 107);

    auto g_xxzz_xxxyy_0 = pbuffer.data(idx_eri_0_gh + 108);

    auto g_xxzz_xxxyz_0 = pbuffer.data(idx_eri_0_gh + 109);

    auto g_xxzz_xxxzz_0 = pbuffer.data(idx_eri_0_gh + 110);

    auto g_xxzz_xxyyy_0 = pbuffer.data(idx_eri_0_gh + 111);

    auto g_xxzz_xxyyz_0 = pbuffer.data(idx_eri_0_gh + 112);

    auto g_xxzz_xxyzz_0 = pbuffer.data(idx_eri_0_gh + 113);

    auto g_xxzz_xxzzz_0 = pbuffer.data(idx_eri_0_gh + 114);

    auto g_xxzz_xyyyy_0 = pbuffer.data(idx_eri_0_gh + 115);

    auto g_xxzz_xyyyz_0 = pbuffer.data(idx_eri_0_gh + 116);

    auto g_xxzz_xyyzz_0 = pbuffer.data(idx_eri_0_gh + 117);

    auto g_xxzz_xyzzz_0 = pbuffer.data(idx_eri_0_gh + 118);

    auto g_xxzz_xzzzz_0 = pbuffer.data(idx_eri_0_gh + 119);

    auto g_xxzz_yyyyy_0 = pbuffer.data(idx_eri_0_gh + 120);

    auto g_xxzz_yyyyz_0 = pbuffer.data(idx_eri_0_gh + 121);

    auto g_xxzz_yyyzz_0 = pbuffer.data(idx_eri_0_gh + 122);

    auto g_xxzz_yyzzz_0 = pbuffer.data(idx_eri_0_gh + 123);

    auto g_xxzz_yzzzz_0 = pbuffer.data(idx_eri_0_gh + 124);

    auto g_xxzz_zzzzz_0 = pbuffer.data(idx_eri_0_gh + 125);

    #pragma omp simd aligned(g_xx_xxxxx_0, g_xx_xxxxx_1, g_xx_xxxxy_0, g_xx_xxxxy_1, g_xx_xxxyy_0, g_xx_xxxyy_1, g_xx_xxyyy_0, g_xx_xxyyy_1, g_xx_xyyyy_0, g_xx_xyyyy_1, g_xxz_xxxxx_1, g_xxz_xxxxy_1, g_xxz_xxxyy_1, g_xxz_xxyyy_1, g_xxz_xyyyy_1, g_xxzz_xxxxx_0, g_xxzz_xxxxy_0, g_xxzz_xxxxz_0, g_xxzz_xxxyy_0, g_xxzz_xxxyz_0, g_xxzz_xxxzz_0, g_xxzz_xxyyy_0, g_xxzz_xxyyz_0, g_xxzz_xxyzz_0, g_xxzz_xxzzz_0, g_xxzz_xyyyy_0, g_xxzz_xyyyz_0, g_xxzz_xyyzz_0, g_xxzz_xyzzz_0, g_xxzz_xzzzz_0, g_xxzz_yyyyy_0, g_xxzz_yyyyz_0, g_xxzz_yyyzz_0, g_xxzz_yyzzz_0, g_xxzz_yzzzz_0, g_xxzz_zzzzz_0, g_xzz_xxxxz_1, g_xzz_xxxyz_1, g_xzz_xxxz_1, g_xzz_xxxzz_1, g_xzz_xxyyz_1, g_xzz_xxyz_1, g_xzz_xxyzz_1, g_xzz_xxzz_1, g_xzz_xxzzz_1, g_xzz_xyyyz_1, g_xzz_xyyz_1, g_xzz_xyyzz_1, g_xzz_xyzz_1, g_xzz_xyzzz_1, g_xzz_xzzz_1, g_xzz_xzzzz_1, g_xzz_yyyyy_1, g_xzz_yyyyz_1, g_xzz_yyyz_1, g_xzz_yyyzz_1, g_xzz_yyzz_1, g_xzz_yyzzz_1, g_xzz_yzzz_1, g_xzz_yzzzz_1, g_xzz_zzzz_1, g_xzz_zzzzz_1, g_zz_xxxxz_0, g_zz_xxxxz_1, g_zz_xxxyz_0, g_zz_xxxyz_1, g_zz_xxxzz_0, g_zz_xxxzz_1, g_zz_xxyyz_0, g_zz_xxyyz_1, g_zz_xxyzz_0, g_zz_xxyzz_1, g_zz_xxzzz_0, g_zz_xxzzz_1, g_zz_xyyyz_0, g_zz_xyyyz_1, g_zz_xyyzz_0, g_zz_xyyzz_1, g_zz_xyzzz_0, g_zz_xyzzz_1, g_zz_xzzzz_0, g_zz_xzzzz_1, g_zz_yyyyy_0, g_zz_yyyyy_1, g_zz_yyyyz_0, g_zz_yyyyz_1, g_zz_yyyzz_0, g_zz_yyyzz_1, g_zz_yyzzz_0, g_zz_yyzzz_1, g_zz_yzzzz_0, g_zz_yzzzz_1, g_zz_zzzzz_0, g_zz_zzzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxzz_xxxxx_0[i] = g_xx_xxxxx_0[i] * fbe_0 - g_xx_xxxxx_1[i] * fz_be_0 + g_xxz_xxxxx_1[i] * pa_z[i];

        g_xxzz_xxxxy_0[i] = g_xx_xxxxy_0[i] * fbe_0 - g_xx_xxxxy_1[i] * fz_be_0 + g_xxz_xxxxy_1[i] * pa_z[i];

        g_xxzz_xxxxz_0[i] = g_zz_xxxxz_0[i] * fbe_0 - g_zz_xxxxz_1[i] * fz_be_0 + 4.0 * g_xzz_xxxz_1[i] * fe_0 + g_xzz_xxxxz_1[i] * pa_x[i];

        g_xxzz_xxxyy_0[i] = g_xx_xxxyy_0[i] * fbe_0 - g_xx_xxxyy_1[i] * fz_be_0 + g_xxz_xxxyy_1[i] * pa_z[i];

        g_xxzz_xxxyz_0[i] = g_zz_xxxyz_0[i] * fbe_0 - g_zz_xxxyz_1[i] * fz_be_0 + 3.0 * g_xzz_xxyz_1[i] * fe_0 + g_xzz_xxxyz_1[i] * pa_x[i];

        g_xxzz_xxxzz_0[i] = g_zz_xxxzz_0[i] * fbe_0 - g_zz_xxxzz_1[i] * fz_be_0 + 3.0 * g_xzz_xxzz_1[i] * fe_0 + g_xzz_xxxzz_1[i] * pa_x[i];

        g_xxzz_xxyyy_0[i] = g_xx_xxyyy_0[i] * fbe_0 - g_xx_xxyyy_1[i] * fz_be_0 + g_xxz_xxyyy_1[i] * pa_z[i];

        g_xxzz_xxyyz_0[i] = g_zz_xxyyz_0[i] * fbe_0 - g_zz_xxyyz_1[i] * fz_be_0 + 2.0 * g_xzz_xyyz_1[i] * fe_0 + g_xzz_xxyyz_1[i] * pa_x[i];

        g_xxzz_xxyzz_0[i] = g_zz_xxyzz_0[i] * fbe_0 - g_zz_xxyzz_1[i] * fz_be_0 + 2.0 * g_xzz_xyzz_1[i] * fe_0 + g_xzz_xxyzz_1[i] * pa_x[i];

        g_xxzz_xxzzz_0[i] = g_zz_xxzzz_0[i] * fbe_0 - g_zz_xxzzz_1[i] * fz_be_0 + 2.0 * g_xzz_xzzz_1[i] * fe_0 + g_xzz_xxzzz_1[i] * pa_x[i];

        g_xxzz_xyyyy_0[i] = g_xx_xyyyy_0[i] * fbe_0 - g_xx_xyyyy_1[i] * fz_be_0 + g_xxz_xyyyy_1[i] * pa_z[i];

        g_xxzz_xyyyz_0[i] = g_zz_xyyyz_0[i] * fbe_0 - g_zz_xyyyz_1[i] * fz_be_0 + g_xzz_yyyz_1[i] * fe_0 + g_xzz_xyyyz_1[i] * pa_x[i];

        g_xxzz_xyyzz_0[i] = g_zz_xyyzz_0[i] * fbe_0 - g_zz_xyyzz_1[i] * fz_be_0 + g_xzz_yyzz_1[i] * fe_0 + g_xzz_xyyzz_1[i] * pa_x[i];

        g_xxzz_xyzzz_0[i] = g_zz_xyzzz_0[i] * fbe_0 - g_zz_xyzzz_1[i] * fz_be_0 + g_xzz_yzzz_1[i] * fe_0 + g_xzz_xyzzz_1[i] * pa_x[i];

        g_xxzz_xzzzz_0[i] = g_zz_xzzzz_0[i] * fbe_0 - g_zz_xzzzz_1[i] * fz_be_0 + g_xzz_zzzz_1[i] * fe_0 + g_xzz_xzzzz_1[i] * pa_x[i];

        g_xxzz_yyyyy_0[i] = g_zz_yyyyy_0[i] * fbe_0 - g_zz_yyyyy_1[i] * fz_be_0 + g_xzz_yyyyy_1[i] * pa_x[i];

        g_xxzz_yyyyz_0[i] = g_zz_yyyyz_0[i] * fbe_0 - g_zz_yyyyz_1[i] * fz_be_0 + g_xzz_yyyyz_1[i] * pa_x[i];

        g_xxzz_yyyzz_0[i] = g_zz_yyyzz_0[i] * fbe_0 - g_zz_yyyzz_1[i] * fz_be_0 + g_xzz_yyyzz_1[i] * pa_x[i];

        g_xxzz_yyzzz_0[i] = g_zz_yyzzz_0[i] * fbe_0 - g_zz_yyzzz_1[i] * fz_be_0 + g_xzz_yyzzz_1[i] * pa_x[i];

        g_xxzz_yzzzz_0[i] = g_zz_yzzzz_0[i] * fbe_0 - g_zz_yzzzz_1[i] * fz_be_0 + g_xzz_yzzzz_1[i] * pa_x[i];

        g_xxzz_zzzzz_0[i] = g_zz_zzzzz_0[i] * fbe_0 - g_zz_zzzzz_1[i] * fz_be_0 + g_xzz_zzzzz_1[i] * pa_x[i];
    }

    // Set up 126-147 components of targeted buffer : GH

    auto g_xyyy_xxxxx_0 = pbuffer.data(idx_eri_0_gh + 126);

    auto g_xyyy_xxxxy_0 = pbuffer.data(idx_eri_0_gh + 127);

    auto g_xyyy_xxxxz_0 = pbuffer.data(idx_eri_0_gh + 128);

    auto g_xyyy_xxxyy_0 = pbuffer.data(idx_eri_0_gh + 129);

    auto g_xyyy_xxxyz_0 = pbuffer.data(idx_eri_0_gh + 130);

    auto g_xyyy_xxxzz_0 = pbuffer.data(idx_eri_0_gh + 131);

    auto g_xyyy_xxyyy_0 = pbuffer.data(idx_eri_0_gh + 132);

    auto g_xyyy_xxyyz_0 = pbuffer.data(idx_eri_0_gh + 133);

    auto g_xyyy_xxyzz_0 = pbuffer.data(idx_eri_0_gh + 134);

    auto g_xyyy_xxzzz_0 = pbuffer.data(idx_eri_0_gh + 135);

    auto g_xyyy_xyyyy_0 = pbuffer.data(idx_eri_0_gh + 136);

    auto g_xyyy_xyyyz_0 = pbuffer.data(idx_eri_0_gh + 137);

    auto g_xyyy_xyyzz_0 = pbuffer.data(idx_eri_0_gh + 138);

    auto g_xyyy_xyzzz_0 = pbuffer.data(idx_eri_0_gh + 139);

    auto g_xyyy_xzzzz_0 = pbuffer.data(idx_eri_0_gh + 140);

    auto g_xyyy_yyyyy_0 = pbuffer.data(idx_eri_0_gh + 141);

    auto g_xyyy_yyyyz_0 = pbuffer.data(idx_eri_0_gh + 142);

    auto g_xyyy_yyyzz_0 = pbuffer.data(idx_eri_0_gh + 143);

    auto g_xyyy_yyzzz_0 = pbuffer.data(idx_eri_0_gh + 144);

    auto g_xyyy_yzzzz_0 = pbuffer.data(idx_eri_0_gh + 145);

    auto g_xyyy_zzzzz_0 = pbuffer.data(idx_eri_0_gh + 146);

    #pragma omp simd aligned(g_xyyy_xxxxx_0, g_xyyy_xxxxy_0, g_xyyy_xxxxz_0, g_xyyy_xxxyy_0, g_xyyy_xxxyz_0, g_xyyy_xxxzz_0, g_xyyy_xxyyy_0, g_xyyy_xxyyz_0, g_xyyy_xxyzz_0, g_xyyy_xxzzz_0, g_xyyy_xyyyy_0, g_xyyy_xyyyz_0, g_xyyy_xyyzz_0, g_xyyy_xyzzz_0, g_xyyy_xzzzz_0, g_xyyy_yyyyy_0, g_xyyy_yyyyz_0, g_xyyy_yyyzz_0, g_xyyy_yyzzz_0, g_xyyy_yzzzz_0, g_xyyy_zzzzz_0, g_yyy_xxxx_1, g_yyy_xxxxx_1, g_yyy_xxxxy_1, g_yyy_xxxxz_1, g_yyy_xxxy_1, g_yyy_xxxyy_1, g_yyy_xxxyz_1, g_yyy_xxxz_1, g_yyy_xxxzz_1, g_yyy_xxyy_1, g_yyy_xxyyy_1, g_yyy_xxyyz_1, g_yyy_xxyz_1, g_yyy_xxyzz_1, g_yyy_xxzz_1, g_yyy_xxzzz_1, g_yyy_xyyy_1, g_yyy_xyyyy_1, g_yyy_xyyyz_1, g_yyy_xyyz_1, g_yyy_xyyzz_1, g_yyy_xyzz_1, g_yyy_xyzzz_1, g_yyy_xzzz_1, g_yyy_xzzzz_1, g_yyy_yyyy_1, g_yyy_yyyyy_1, g_yyy_yyyyz_1, g_yyy_yyyz_1, g_yyy_yyyzz_1, g_yyy_yyzz_1, g_yyy_yyzzz_1, g_yyy_yzzz_1, g_yyy_yzzzz_1, g_yyy_zzzz_1, g_yyy_zzzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyy_xxxxx_0[i] = 5.0 * g_yyy_xxxx_1[i] * fe_0 + g_yyy_xxxxx_1[i] * pa_x[i];

        g_xyyy_xxxxy_0[i] = 4.0 * g_yyy_xxxy_1[i] * fe_0 + g_yyy_xxxxy_1[i] * pa_x[i];

        g_xyyy_xxxxz_0[i] = 4.0 * g_yyy_xxxz_1[i] * fe_0 + g_yyy_xxxxz_1[i] * pa_x[i];

        g_xyyy_xxxyy_0[i] = 3.0 * g_yyy_xxyy_1[i] * fe_0 + g_yyy_xxxyy_1[i] * pa_x[i];

        g_xyyy_xxxyz_0[i] = 3.0 * g_yyy_xxyz_1[i] * fe_0 + g_yyy_xxxyz_1[i] * pa_x[i];

        g_xyyy_xxxzz_0[i] = 3.0 * g_yyy_xxzz_1[i] * fe_0 + g_yyy_xxxzz_1[i] * pa_x[i];

        g_xyyy_xxyyy_0[i] = 2.0 * g_yyy_xyyy_1[i] * fe_0 + g_yyy_xxyyy_1[i] * pa_x[i];

        g_xyyy_xxyyz_0[i] = 2.0 * g_yyy_xyyz_1[i] * fe_0 + g_yyy_xxyyz_1[i] * pa_x[i];

        g_xyyy_xxyzz_0[i] = 2.0 * g_yyy_xyzz_1[i] * fe_0 + g_yyy_xxyzz_1[i] * pa_x[i];

        g_xyyy_xxzzz_0[i] = 2.0 * g_yyy_xzzz_1[i] * fe_0 + g_yyy_xxzzz_1[i] * pa_x[i];

        g_xyyy_xyyyy_0[i] = g_yyy_yyyy_1[i] * fe_0 + g_yyy_xyyyy_1[i] * pa_x[i];

        g_xyyy_xyyyz_0[i] = g_yyy_yyyz_1[i] * fe_0 + g_yyy_xyyyz_1[i] * pa_x[i];

        g_xyyy_xyyzz_0[i] = g_yyy_yyzz_1[i] * fe_0 + g_yyy_xyyzz_1[i] * pa_x[i];

        g_xyyy_xyzzz_0[i] = g_yyy_yzzz_1[i] * fe_0 + g_yyy_xyzzz_1[i] * pa_x[i];

        g_xyyy_xzzzz_0[i] = g_yyy_zzzz_1[i] * fe_0 + g_yyy_xzzzz_1[i] * pa_x[i];

        g_xyyy_yyyyy_0[i] = g_yyy_yyyyy_1[i] * pa_x[i];

        g_xyyy_yyyyz_0[i] = g_yyy_yyyyz_1[i] * pa_x[i];

        g_xyyy_yyyzz_0[i] = g_yyy_yyyzz_1[i] * pa_x[i];

        g_xyyy_yyzzz_0[i] = g_yyy_yyzzz_1[i] * pa_x[i];

        g_xyyy_yzzzz_0[i] = g_yyy_yzzzz_1[i] * pa_x[i];

        g_xyyy_zzzzz_0[i] = g_yyy_zzzzz_1[i] * pa_x[i];
    }

    // Set up 147-168 components of targeted buffer : GH

    auto g_xyyz_xxxxx_0 = pbuffer.data(idx_eri_0_gh + 147);

    auto g_xyyz_xxxxy_0 = pbuffer.data(idx_eri_0_gh + 148);

    auto g_xyyz_xxxxz_0 = pbuffer.data(idx_eri_0_gh + 149);

    auto g_xyyz_xxxyy_0 = pbuffer.data(idx_eri_0_gh + 150);

    auto g_xyyz_xxxyz_0 = pbuffer.data(idx_eri_0_gh + 151);

    auto g_xyyz_xxxzz_0 = pbuffer.data(idx_eri_0_gh + 152);

    auto g_xyyz_xxyyy_0 = pbuffer.data(idx_eri_0_gh + 153);

    auto g_xyyz_xxyyz_0 = pbuffer.data(idx_eri_0_gh + 154);

    auto g_xyyz_xxyzz_0 = pbuffer.data(idx_eri_0_gh + 155);

    auto g_xyyz_xxzzz_0 = pbuffer.data(idx_eri_0_gh + 156);

    auto g_xyyz_xyyyy_0 = pbuffer.data(idx_eri_0_gh + 157);

    auto g_xyyz_xyyyz_0 = pbuffer.data(idx_eri_0_gh + 158);

    auto g_xyyz_xyyzz_0 = pbuffer.data(idx_eri_0_gh + 159);

    auto g_xyyz_xyzzz_0 = pbuffer.data(idx_eri_0_gh + 160);

    auto g_xyyz_xzzzz_0 = pbuffer.data(idx_eri_0_gh + 161);

    auto g_xyyz_yyyyy_0 = pbuffer.data(idx_eri_0_gh + 162);

    auto g_xyyz_yyyyz_0 = pbuffer.data(idx_eri_0_gh + 163);

    auto g_xyyz_yyyzz_0 = pbuffer.data(idx_eri_0_gh + 164);

    auto g_xyyz_yyzzz_0 = pbuffer.data(idx_eri_0_gh + 165);

    auto g_xyyz_yzzzz_0 = pbuffer.data(idx_eri_0_gh + 166);

    auto g_xyyz_zzzzz_0 = pbuffer.data(idx_eri_0_gh + 167);

    #pragma omp simd aligned(g_xyy_xxxxx_1, g_xyy_xxxxy_1, g_xyy_xxxyy_1, g_xyy_xxyyy_1, g_xyy_xyyyy_1, g_xyyz_xxxxx_0, g_xyyz_xxxxy_0, g_xyyz_xxxxz_0, g_xyyz_xxxyy_0, g_xyyz_xxxyz_0, g_xyyz_xxxzz_0, g_xyyz_xxyyy_0, g_xyyz_xxyyz_0, g_xyyz_xxyzz_0, g_xyyz_xxzzz_0, g_xyyz_xyyyy_0, g_xyyz_xyyyz_0, g_xyyz_xyyzz_0, g_xyyz_xyzzz_0, g_xyyz_xzzzz_0, g_xyyz_yyyyy_0, g_xyyz_yyyyz_0, g_xyyz_yyyzz_0, g_xyyz_yyzzz_0, g_xyyz_yzzzz_0, g_xyyz_zzzzz_0, g_yyz_xxxxz_1, g_yyz_xxxyz_1, g_yyz_xxxz_1, g_yyz_xxxzz_1, g_yyz_xxyyz_1, g_yyz_xxyz_1, g_yyz_xxyzz_1, g_yyz_xxzz_1, g_yyz_xxzzz_1, g_yyz_xyyyz_1, g_yyz_xyyz_1, g_yyz_xyyzz_1, g_yyz_xyzz_1, g_yyz_xyzzz_1, g_yyz_xzzz_1, g_yyz_xzzzz_1, g_yyz_yyyyy_1, g_yyz_yyyyz_1, g_yyz_yyyz_1, g_yyz_yyyzz_1, g_yyz_yyzz_1, g_yyz_yyzzz_1, g_yyz_yzzz_1, g_yyz_yzzzz_1, g_yyz_zzzz_1, g_yyz_zzzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyz_xxxxx_0[i] = g_xyy_xxxxx_1[i] * pa_z[i];

        g_xyyz_xxxxy_0[i] = g_xyy_xxxxy_1[i] * pa_z[i];

        g_xyyz_xxxxz_0[i] = 4.0 * g_yyz_xxxz_1[i] * fe_0 + g_yyz_xxxxz_1[i] * pa_x[i];

        g_xyyz_xxxyy_0[i] = g_xyy_xxxyy_1[i] * pa_z[i];

        g_xyyz_xxxyz_0[i] = 3.0 * g_yyz_xxyz_1[i] * fe_0 + g_yyz_xxxyz_1[i] * pa_x[i];

        g_xyyz_xxxzz_0[i] = 3.0 * g_yyz_xxzz_1[i] * fe_0 + g_yyz_xxxzz_1[i] * pa_x[i];

        g_xyyz_xxyyy_0[i] = g_xyy_xxyyy_1[i] * pa_z[i];

        g_xyyz_xxyyz_0[i] = 2.0 * g_yyz_xyyz_1[i] * fe_0 + g_yyz_xxyyz_1[i] * pa_x[i];

        g_xyyz_xxyzz_0[i] = 2.0 * g_yyz_xyzz_1[i] * fe_0 + g_yyz_xxyzz_1[i] * pa_x[i];

        g_xyyz_xxzzz_0[i] = 2.0 * g_yyz_xzzz_1[i] * fe_0 + g_yyz_xxzzz_1[i] * pa_x[i];

        g_xyyz_xyyyy_0[i] = g_xyy_xyyyy_1[i] * pa_z[i];

        g_xyyz_xyyyz_0[i] = g_yyz_yyyz_1[i] * fe_0 + g_yyz_xyyyz_1[i] * pa_x[i];

        g_xyyz_xyyzz_0[i] = g_yyz_yyzz_1[i] * fe_0 + g_yyz_xyyzz_1[i] * pa_x[i];

        g_xyyz_xyzzz_0[i] = g_yyz_yzzz_1[i] * fe_0 + g_yyz_xyzzz_1[i] * pa_x[i];

        g_xyyz_xzzzz_0[i] = g_yyz_zzzz_1[i] * fe_0 + g_yyz_xzzzz_1[i] * pa_x[i];

        g_xyyz_yyyyy_0[i] = g_yyz_yyyyy_1[i] * pa_x[i];

        g_xyyz_yyyyz_0[i] = g_yyz_yyyyz_1[i] * pa_x[i];

        g_xyyz_yyyzz_0[i] = g_yyz_yyyzz_1[i] * pa_x[i];

        g_xyyz_yyzzz_0[i] = g_yyz_yyzzz_1[i] * pa_x[i];

        g_xyyz_yzzzz_0[i] = g_yyz_yzzzz_1[i] * pa_x[i];

        g_xyyz_zzzzz_0[i] = g_yyz_zzzzz_1[i] * pa_x[i];
    }

    // Set up 168-189 components of targeted buffer : GH

    auto g_xyzz_xxxxx_0 = pbuffer.data(idx_eri_0_gh + 168);

    auto g_xyzz_xxxxy_0 = pbuffer.data(idx_eri_0_gh + 169);

    auto g_xyzz_xxxxz_0 = pbuffer.data(idx_eri_0_gh + 170);

    auto g_xyzz_xxxyy_0 = pbuffer.data(idx_eri_0_gh + 171);

    auto g_xyzz_xxxyz_0 = pbuffer.data(idx_eri_0_gh + 172);

    auto g_xyzz_xxxzz_0 = pbuffer.data(idx_eri_0_gh + 173);

    auto g_xyzz_xxyyy_0 = pbuffer.data(idx_eri_0_gh + 174);

    auto g_xyzz_xxyyz_0 = pbuffer.data(idx_eri_0_gh + 175);

    auto g_xyzz_xxyzz_0 = pbuffer.data(idx_eri_0_gh + 176);

    auto g_xyzz_xxzzz_0 = pbuffer.data(idx_eri_0_gh + 177);

    auto g_xyzz_xyyyy_0 = pbuffer.data(idx_eri_0_gh + 178);

    auto g_xyzz_xyyyz_0 = pbuffer.data(idx_eri_0_gh + 179);

    auto g_xyzz_xyyzz_0 = pbuffer.data(idx_eri_0_gh + 180);

    auto g_xyzz_xyzzz_0 = pbuffer.data(idx_eri_0_gh + 181);

    auto g_xyzz_xzzzz_0 = pbuffer.data(idx_eri_0_gh + 182);

    auto g_xyzz_yyyyy_0 = pbuffer.data(idx_eri_0_gh + 183);

    auto g_xyzz_yyyyz_0 = pbuffer.data(idx_eri_0_gh + 184);

    auto g_xyzz_yyyzz_0 = pbuffer.data(idx_eri_0_gh + 185);

    auto g_xyzz_yyzzz_0 = pbuffer.data(idx_eri_0_gh + 186);

    auto g_xyzz_yzzzz_0 = pbuffer.data(idx_eri_0_gh + 187);

    auto g_xyzz_zzzzz_0 = pbuffer.data(idx_eri_0_gh + 188);

    #pragma omp simd aligned(g_xyzz_xxxxx_0, g_xyzz_xxxxy_0, g_xyzz_xxxxz_0, g_xyzz_xxxyy_0, g_xyzz_xxxyz_0, g_xyzz_xxxzz_0, g_xyzz_xxyyy_0, g_xyzz_xxyyz_0, g_xyzz_xxyzz_0, g_xyzz_xxzzz_0, g_xyzz_xyyyy_0, g_xyzz_xyyyz_0, g_xyzz_xyyzz_0, g_xyzz_xyzzz_0, g_xyzz_xzzzz_0, g_xyzz_yyyyy_0, g_xyzz_yyyyz_0, g_xyzz_yyyzz_0, g_xyzz_yyzzz_0, g_xyzz_yzzzz_0, g_xyzz_zzzzz_0, g_xzz_xxxxx_1, g_xzz_xxxxz_1, g_xzz_xxxzz_1, g_xzz_xxzzz_1, g_xzz_xzzzz_1, g_yzz_xxxxy_1, g_yzz_xxxy_1, g_yzz_xxxyy_1, g_yzz_xxxyz_1, g_yzz_xxyy_1, g_yzz_xxyyy_1, g_yzz_xxyyz_1, g_yzz_xxyz_1, g_yzz_xxyzz_1, g_yzz_xyyy_1, g_yzz_xyyyy_1, g_yzz_xyyyz_1, g_yzz_xyyz_1, g_yzz_xyyzz_1, g_yzz_xyzz_1, g_yzz_xyzzz_1, g_yzz_yyyy_1, g_yzz_yyyyy_1, g_yzz_yyyyz_1, g_yzz_yyyz_1, g_yzz_yyyzz_1, g_yzz_yyzz_1, g_yzz_yyzzz_1, g_yzz_yzzz_1, g_yzz_yzzzz_1, g_yzz_zzzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyzz_xxxxx_0[i] = g_xzz_xxxxx_1[i] * pa_y[i];

        g_xyzz_xxxxy_0[i] = 4.0 * g_yzz_xxxy_1[i] * fe_0 + g_yzz_xxxxy_1[i] * pa_x[i];

        g_xyzz_xxxxz_0[i] = g_xzz_xxxxz_1[i] * pa_y[i];

        g_xyzz_xxxyy_0[i] = 3.0 * g_yzz_xxyy_1[i] * fe_0 + g_yzz_xxxyy_1[i] * pa_x[i];

        g_xyzz_xxxyz_0[i] = 3.0 * g_yzz_xxyz_1[i] * fe_0 + g_yzz_xxxyz_1[i] * pa_x[i];

        g_xyzz_xxxzz_0[i] = g_xzz_xxxzz_1[i] * pa_y[i];

        g_xyzz_xxyyy_0[i] = 2.0 * g_yzz_xyyy_1[i] * fe_0 + g_yzz_xxyyy_1[i] * pa_x[i];

        g_xyzz_xxyyz_0[i] = 2.0 * g_yzz_xyyz_1[i] * fe_0 + g_yzz_xxyyz_1[i] * pa_x[i];

        g_xyzz_xxyzz_0[i] = 2.0 * g_yzz_xyzz_1[i] * fe_0 + g_yzz_xxyzz_1[i] * pa_x[i];

        g_xyzz_xxzzz_0[i] = g_xzz_xxzzz_1[i] * pa_y[i];

        g_xyzz_xyyyy_0[i] = g_yzz_yyyy_1[i] * fe_0 + g_yzz_xyyyy_1[i] * pa_x[i];

        g_xyzz_xyyyz_0[i] = g_yzz_yyyz_1[i] * fe_0 + g_yzz_xyyyz_1[i] * pa_x[i];

        g_xyzz_xyyzz_0[i] = g_yzz_yyzz_1[i] * fe_0 + g_yzz_xyyzz_1[i] * pa_x[i];

        g_xyzz_xyzzz_0[i] = g_yzz_yzzz_1[i] * fe_0 + g_yzz_xyzzz_1[i] * pa_x[i];

        g_xyzz_xzzzz_0[i] = g_xzz_xzzzz_1[i] * pa_y[i];

        g_xyzz_yyyyy_0[i] = g_yzz_yyyyy_1[i] * pa_x[i];

        g_xyzz_yyyyz_0[i] = g_yzz_yyyyz_1[i] * pa_x[i];

        g_xyzz_yyyzz_0[i] = g_yzz_yyyzz_1[i] * pa_x[i];

        g_xyzz_yyzzz_0[i] = g_yzz_yyzzz_1[i] * pa_x[i];

        g_xyzz_yzzzz_0[i] = g_yzz_yzzzz_1[i] * pa_x[i];

        g_xyzz_zzzzz_0[i] = g_yzz_zzzzz_1[i] * pa_x[i];
    }

    // Set up 189-210 components of targeted buffer : GH

    auto g_xzzz_xxxxx_0 = pbuffer.data(idx_eri_0_gh + 189);

    auto g_xzzz_xxxxy_0 = pbuffer.data(idx_eri_0_gh + 190);

    auto g_xzzz_xxxxz_0 = pbuffer.data(idx_eri_0_gh + 191);

    auto g_xzzz_xxxyy_0 = pbuffer.data(idx_eri_0_gh + 192);

    auto g_xzzz_xxxyz_0 = pbuffer.data(idx_eri_0_gh + 193);

    auto g_xzzz_xxxzz_0 = pbuffer.data(idx_eri_0_gh + 194);

    auto g_xzzz_xxyyy_0 = pbuffer.data(idx_eri_0_gh + 195);

    auto g_xzzz_xxyyz_0 = pbuffer.data(idx_eri_0_gh + 196);

    auto g_xzzz_xxyzz_0 = pbuffer.data(idx_eri_0_gh + 197);

    auto g_xzzz_xxzzz_0 = pbuffer.data(idx_eri_0_gh + 198);

    auto g_xzzz_xyyyy_0 = pbuffer.data(idx_eri_0_gh + 199);

    auto g_xzzz_xyyyz_0 = pbuffer.data(idx_eri_0_gh + 200);

    auto g_xzzz_xyyzz_0 = pbuffer.data(idx_eri_0_gh + 201);

    auto g_xzzz_xyzzz_0 = pbuffer.data(idx_eri_0_gh + 202);

    auto g_xzzz_xzzzz_0 = pbuffer.data(idx_eri_0_gh + 203);

    auto g_xzzz_yyyyy_0 = pbuffer.data(idx_eri_0_gh + 204);

    auto g_xzzz_yyyyz_0 = pbuffer.data(idx_eri_0_gh + 205);

    auto g_xzzz_yyyzz_0 = pbuffer.data(idx_eri_0_gh + 206);

    auto g_xzzz_yyzzz_0 = pbuffer.data(idx_eri_0_gh + 207);

    auto g_xzzz_yzzzz_0 = pbuffer.data(idx_eri_0_gh + 208);

    auto g_xzzz_zzzzz_0 = pbuffer.data(idx_eri_0_gh + 209);

    #pragma omp simd aligned(g_xzzz_xxxxx_0, g_xzzz_xxxxy_0, g_xzzz_xxxxz_0, g_xzzz_xxxyy_0, g_xzzz_xxxyz_0, g_xzzz_xxxzz_0, g_xzzz_xxyyy_0, g_xzzz_xxyyz_0, g_xzzz_xxyzz_0, g_xzzz_xxzzz_0, g_xzzz_xyyyy_0, g_xzzz_xyyyz_0, g_xzzz_xyyzz_0, g_xzzz_xyzzz_0, g_xzzz_xzzzz_0, g_xzzz_yyyyy_0, g_xzzz_yyyyz_0, g_xzzz_yyyzz_0, g_xzzz_yyzzz_0, g_xzzz_yzzzz_0, g_xzzz_zzzzz_0, g_zzz_xxxx_1, g_zzz_xxxxx_1, g_zzz_xxxxy_1, g_zzz_xxxxz_1, g_zzz_xxxy_1, g_zzz_xxxyy_1, g_zzz_xxxyz_1, g_zzz_xxxz_1, g_zzz_xxxzz_1, g_zzz_xxyy_1, g_zzz_xxyyy_1, g_zzz_xxyyz_1, g_zzz_xxyz_1, g_zzz_xxyzz_1, g_zzz_xxzz_1, g_zzz_xxzzz_1, g_zzz_xyyy_1, g_zzz_xyyyy_1, g_zzz_xyyyz_1, g_zzz_xyyz_1, g_zzz_xyyzz_1, g_zzz_xyzz_1, g_zzz_xyzzz_1, g_zzz_xzzz_1, g_zzz_xzzzz_1, g_zzz_yyyy_1, g_zzz_yyyyy_1, g_zzz_yyyyz_1, g_zzz_yyyz_1, g_zzz_yyyzz_1, g_zzz_yyzz_1, g_zzz_yyzzz_1, g_zzz_yzzz_1, g_zzz_yzzzz_1, g_zzz_zzzz_1, g_zzz_zzzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xzzz_xxxxx_0[i] = 5.0 * g_zzz_xxxx_1[i] * fe_0 + g_zzz_xxxxx_1[i] * pa_x[i];

        g_xzzz_xxxxy_0[i] = 4.0 * g_zzz_xxxy_1[i] * fe_0 + g_zzz_xxxxy_1[i] * pa_x[i];

        g_xzzz_xxxxz_0[i] = 4.0 * g_zzz_xxxz_1[i] * fe_0 + g_zzz_xxxxz_1[i] * pa_x[i];

        g_xzzz_xxxyy_0[i] = 3.0 * g_zzz_xxyy_1[i] * fe_0 + g_zzz_xxxyy_1[i] * pa_x[i];

        g_xzzz_xxxyz_0[i] = 3.0 * g_zzz_xxyz_1[i] * fe_0 + g_zzz_xxxyz_1[i] * pa_x[i];

        g_xzzz_xxxzz_0[i] = 3.0 * g_zzz_xxzz_1[i] * fe_0 + g_zzz_xxxzz_1[i] * pa_x[i];

        g_xzzz_xxyyy_0[i] = 2.0 * g_zzz_xyyy_1[i] * fe_0 + g_zzz_xxyyy_1[i] * pa_x[i];

        g_xzzz_xxyyz_0[i] = 2.0 * g_zzz_xyyz_1[i] * fe_0 + g_zzz_xxyyz_1[i] * pa_x[i];

        g_xzzz_xxyzz_0[i] = 2.0 * g_zzz_xyzz_1[i] * fe_0 + g_zzz_xxyzz_1[i] * pa_x[i];

        g_xzzz_xxzzz_0[i] = 2.0 * g_zzz_xzzz_1[i] * fe_0 + g_zzz_xxzzz_1[i] * pa_x[i];

        g_xzzz_xyyyy_0[i] = g_zzz_yyyy_1[i] * fe_0 + g_zzz_xyyyy_1[i] * pa_x[i];

        g_xzzz_xyyyz_0[i] = g_zzz_yyyz_1[i] * fe_0 + g_zzz_xyyyz_1[i] * pa_x[i];

        g_xzzz_xyyzz_0[i] = g_zzz_yyzz_1[i] * fe_0 + g_zzz_xyyzz_1[i] * pa_x[i];

        g_xzzz_xyzzz_0[i] = g_zzz_yzzz_1[i] * fe_0 + g_zzz_xyzzz_1[i] * pa_x[i];

        g_xzzz_xzzzz_0[i] = g_zzz_zzzz_1[i] * fe_0 + g_zzz_xzzzz_1[i] * pa_x[i];

        g_xzzz_yyyyy_0[i] = g_zzz_yyyyy_1[i] * pa_x[i];

        g_xzzz_yyyyz_0[i] = g_zzz_yyyyz_1[i] * pa_x[i];

        g_xzzz_yyyzz_0[i] = g_zzz_yyyzz_1[i] * pa_x[i];

        g_xzzz_yyzzz_0[i] = g_zzz_yyzzz_1[i] * pa_x[i];

        g_xzzz_yzzzz_0[i] = g_zzz_yzzzz_1[i] * pa_x[i];

        g_xzzz_zzzzz_0[i] = g_zzz_zzzzz_1[i] * pa_x[i];
    }

    // Set up 210-231 components of targeted buffer : GH

    auto g_yyyy_xxxxx_0 = pbuffer.data(idx_eri_0_gh + 210);

    auto g_yyyy_xxxxy_0 = pbuffer.data(idx_eri_0_gh + 211);

    auto g_yyyy_xxxxz_0 = pbuffer.data(idx_eri_0_gh + 212);

    auto g_yyyy_xxxyy_0 = pbuffer.data(idx_eri_0_gh + 213);

    auto g_yyyy_xxxyz_0 = pbuffer.data(idx_eri_0_gh + 214);

    auto g_yyyy_xxxzz_0 = pbuffer.data(idx_eri_0_gh + 215);

    auto g_yyyy_xxyyy_0 = pbuffer.data(idx_eri_0_gh + 216);

    auto g_yyyy_xxyyz_0 = pbuffer.data(idx_eri_0_gh + 217);

    auto g_yyyy_xxyzz_0 = pbuffer.data(idx_eri_0_gh + 218);

    auto g_yyyy_xxzzz_0 = pbuffer.data(idx_eri_0_gh + 219);

    auto g_yyyy_xyyyy_0 = pbuffer.data(idx_eri_0_gh + 220);

    auto g_yyyy_xyyyz_0 = pbuffer.data(idx_eri_0_gh + 221);

    auto g_yyyy_xyyzz_0 = pbuffer.data(idx_eri_0_gh + 222);

    auto g_yyyy_xyzzz_0 = pbuffer.data(idx_eri_0_gh + 223);

    auto g_yyyy_xzzzz_0 = pbuffer.data(idx_eri_0_gh + 224);

    auto g_yyyy_yyyyy_0 = pbuffer.data(idx_eri_0_gh + 225);

    auto g_yyyy_yyyyz_0 = pbuffer.data(idx_eri_0_gh + 226);

    auto g_yyyy_yyyzz_0 = pbuffer.data(idx_eri_0_gh + 227);

    auto g_yyyy_yyzzz_0 = pbuffer.data(idx_eri_0_gh + 228);

    auto g_yyyy_yzzzz_0 = pbuffer.data(idx_eri_0_gh + 229);

    auto g_yyyy_zzzzz_0 = pbuffer.data(idx_eri_0_gh + 230);

    #pragma omp simd aligned(g_yy_xxxxx_0, g_yy_xxxxx_1, g_yy_xxxxy_0, g_yy_xxxxy_1, g_yy_xxxxz_0, g_yy_xxxxz_1, g_yy_xxxyy_0, g_yy_xxxyy_1, g_yy_xxxyz_0, g_yy_xxxyz_1, g_yy_xxxzz_0, g_yy_xxxzz_1, g_yy_xxyyy_0, g_yy_xxyyy_1, g_yy_xxyyz_0, g_yy_xxyyz_1, g_yy_xxyzz_0, g_yy_xxyzz_1, g_yy_xxzzz_0, g_yy_xxzzz_1, g_yy_xyyyy_0, g_yy_xyyyy_1, g_yy_xyyyz_0, g_yy_xyyyz_1, g_yy_xyyzz_0, g_yy_xyyzz_1, g_yy_xyzzz_0, g_yy_xyzzz_1, g_yy_xzzzz_0, g_yy_xzzzz_1, g_yy_yyyyy_0, g_yy_yyyyy_1, g_yy_yyyyz_0, g_yy_yyyyz_1, g_yy_yyyzz_0, g_yy_yyyzz_1, g_yy_yyzzz_0, g_yy_yyzzz_1, g_yy_yzzzz_0, g_yy_yzzzz_1, g_yy_zzzzz_0, g_yy_zzzzz_1, g_yyy_xxxx_1, g_yyy_xxxxx_1, g_yyy_xxxxy_1, g_yyy_xxxxz_1, g_yyy_xxxy_1, g_yyy_xxxyy_1, g_yyy_xxxyz_1, g_yyy_xxxz_1, g_yyy_xxxzz_1, g_yyy_xxyy_1, g_yyy_xxyyy_1, g_yyy_xxyyz_1, g_yyy_xxyz_1, g_yyy_xxyzz_1, g_yyy_xxzz_1, g_yyy_xxzzz_1, g_yyy_xyyy_1, g_yyy_xyyyy_1, g_yyy_xyyyz_1, g_yyy_xyyz_1, g_yyy_xyyzz_1, g_yyy_xyzz_1, g_yyy_xyzzz_1, g_yyy_xzzz_1, g_yyy_xzzzz_1, g_yyy_yyyy_1, g_yyy_yyyyy_1, g_yyy_yyyyz_1, g_yyy_yyyz_1, g_yyy_yyyzz_1, g_yyy_yyzz_1, g_yyy_yyzzz_1, g_yyy_yzzz_1, g_yyy_yzzzz_1, g_yyy_zzzz_1, g_yyy_zzzzz_1, g_yyyy_xxxxx_0, g_yyyy_xxxxy_0, g_yyyy_xxxxz_0, g_yyyy_xxxyy_0, g_yyyy_xxxyz_0, g_yyyy_xxxzz_0, g_yyyy_xxyyy_0, g_yyyy_xxyyz_0, g_yyyy_xxyzz_0, g_yyyy_xxzzz_0, g_yyyy_xyyyy_0, g_yyyy_xyyyz_0, g_yyyy_xyyzz_0, g_yyyy_xyzzz_0, g_yyyy_xzzzz_0, g_yyyy_yyyyy_0, g_yyyy_yyyyz_0, g_yyyy_yyyzz_0, g_yyyy_yyzzz_0, g_yyyy_yzzzz_0, g_yyyy_zzzzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyy_xxxxx_0[i] = 3.0 * g_yy_xxxxx_0[i] * fbe_0 - 3.0 * g_yy_xxxxx_1[i] * fz_be_0 + g_yyy_xxxxx_1[i] * pa_y[i];

        g_yyyy_xxxxy_0[i] = 3.0 * g_yy_xxxxy_0[i] * fbe_0 - 3.0 * g_yy_xxxxy_1[i] * fz_be_0 + g_yyy_xxxx_1[i] * fe_0 + g_yyy_xxxxy_1[i] * pa_y[i];

        g_yyyy_xxxxz_0[i] = 3.0 * g_yy_xxxxz_0[i] * fbe_0 - 3.0 * g_yy_xxxxz_1[i] * fz_be_0 + g_yyy_xxxxz_1[i] * pa_y[i];

        g_yyyy_xxxyy_0[i] = 3.0 * g_yy_xxxyy_0[i] * fbe_0 - 3.0 * g_yy_xxxyy_1[i] * fz_be_0 + 2.0 * g_yyy_xxxy_1[i] * fe_0 + g_yyy_xxxyy_1[i] * pa_y[i];

        g_yyyy_xxxyz_0[i] = 3.0 * g_yy_xxxyz_0[i] * fbe_0 - 3.0 * g_yy_xxxyz_1[i] * fz_be_0 + g_yyy_xxxz_1[i] * fe_0 + g_yyy_xxxyz_1[i] * pa_y[i];

        g_yyyy_xxxzz_0[i] = 3.0 * g_yy_xxxzz_0[i] * fbe_0 - 3.0 * g_yy_xxxzz_1[i] * fz_be_0 + g_yyy_xxxzz_1[i] * pa_y[i];

        g_yyyy_xxyyy_0[i] = 3.0 * g_yy_xxyyy_0[i] * fbe_0 - 3.0 * g_yy_xxyyy_1[i] * fz_be_0 + 3.0 * g_yyy_xxyy_1[i] * fe_0 + g_yyy_xxyyy_1[i] * pa_y[i];

        g_yyyy_xxyyz_0[i] = 3.0 * g_yy_xxyyz_0[i] * fbe_0 - 3.0 * g_yy_xxyyz_1[i] * fz_be_0 + 2.0 * g_yyy_xxyz_1[i] * fe_0 + g_yyy_xxyyz_1[i] * pa_y[i];

        g_yyyy_xxyzz_0[i] = 3.0 * g_yy_xxyzz_0[i] * fbe_0 - 3.0 * g_yy_xxyzz_1[i] * fz_be_0 + g_yyy_xxzz_1[i] * fe_0 + g_yyy_xxyzz_1[i] * pa_y[i];

        g_yyyy_xxzzz_0[i] = 3.0 * g_yy_xxzzz_0[i] * fbe_0 - 3.0 * g_yy_xxzzz_1[i] * fz_be_0 + g_yyy_xxzzz_1[i] * pa_y[i];

        g_yyyy_xyyyy_0[i] = 3.0 * g_yy_xyyyy_0[i] * fbe_0 - 3.0 * g_yy_xyyyy_1[i] * fz_be_0 + 4.0 * g_yyy_xyyy_1[i] * fe_0 + g_yyy_xyyyy_1[i] * pa_y[i];

        g_yyyy_xyyyz_0[i] = 3.0 * g_yy_xyyyz_0[i] * fbe_0 - 3.0 * g_yy_xyyyz_1[i] * fz_be_0 + 3.0 * g_yyy_xyyz_1[i] * fe_0 + g_yyy_xyyyz_1[i] * pa_y[i];

        g_yyyy_xyyzz_0[i] = 3.0 * g_yy_xyyzz_0[i] * fbe_0 - 3.0 * g_yy_xyyzz_1[i] * fz_be_0 + 2.0 * g_yyy_xyzz_1[i] * fe_0 + g_yyy_xyyzz_1[i] * pa_y[i];

        g_yyyy_xyzzz_0[i] = 3.0 * g_yy_xyzzz_0[i] * fbe_0 - 3.0 * g_yy_xyzzz_1[i] * fz_be_0 + g_yyy_xzzz_1[i] * fe_0 + g_yyy_xyzzz_1[i] * pa_y[i];

        g_yyyy_xzzzz_0[i] = 3.0 * g_yy_xzzzz_0[i] * fbe_0 - 3.0 * g_yy_xzzzz_1[i] * fz_be_0 + g_yyy_xzzzz_1[i] * pa_y[i];

        g_yyyy_yyyyy_0[i] = 3.0 * g_yy_yyyyy_0[i] * fbe_0 - 3.0 * g_yy_yyyyy_1[i] * fz_be_0 + 5.0 * g_yyy_yyyy_1[i] * fe_0 + g_yyy_yyyyy_1[i] * pa_y[i];

        g_yyyy_yyyyz_0[i] = 3.0 * g_yy_yyyyz_0[i] * fbe_0 - 3.0 * g_yy_yyyyz_1[i] * fz_be_0 + 4.0 * g_yyy_yyyz_1[i] * fe_0 + g_yyy_yyyyz_1[i] * pa_y[i];

        g_yyyy_yyyzz_0[i] = 3.0 * g_yy_yyyzz_0[i] * fbe_0 - 3.0 * g_yy_yyyzz_1[i] * fz_be_0 + 3.0 * g_yyy_yyzz_1[i] * fe_0 + g_yyy_yyyzz_1[i] * pa_y[i];

        g_yyyy_yyzzz_0[i] = 3.0 * g_yy_yyzzz_0[i] * fbe_0 - 3.0 * g_yy_yyzzz_1[i] * fz_be_0 + 2.0 * g_yyy_yzzz_1[i] * fe_0 + g_yyy_yyzzz_1[i] * pa_y[i];

        g_yyyy_yzzzz_0[i] = 3.0 * g_yy_yzzzz_0[i] * fbe_0 - 3.0 * g_yy_yzzzz_1[i] * fz_be_0 + g_yyy_zzzz_1[i] * fe_0 + g_yyy_yzzzz_1[i] * pa_y[i];

        g_yyyy_zzzzz_0[i] = 3.0 * g_yy_zzzzz_0[i] * fbe_0 - 3.0 * g_yy_zzzzz_1[i] * fz_be_0 + g_yyy_zzzzz_1[i] * pa_y[i];
    }

    // Set up 231-252 components of targeted buffer : GH

    auto g_yyyz_xxxxx_0 = pbuffer.data(idx_eri_0_gh + 231);

    auto g_yyyz_xxxxy_0 = pbuffer.data(idx_eri_0_gh + 232);

    auto g_yyyz_xxxxz_0 = pbuffer.data(idx_eri_0_gh + 233);

    auto g_yyyz_xxxyy_0 = pbuffer.data(idx_eri_0_gh + 234);

    auto g_yyyz_xxxyz_0 = pbuffer.data(idx_eri_0_gh + 235);

    auto g_yyyz_xxxzz_0 = pbuffer.data(idx_eri_0_gh + 236);

    auto g_yyyz_xxyyy_0 = pbuffer.data(idx_eri_0_gh + 237);

    auto g_yyyz_xxyyz_0 = pbuffer.data(idx_eri_0_gh + 238);

    auto g_yyyz_xxyzz_0 = pbuffer.data(idx_eri_0_gh + 239);

    auto g_yyyz_xxzzz_0 = pbuffer.data(idx_eri_0_gh + 240);

    auto g_yyyz_xyyyy_0 = pbuffer.data(idx_eri_0_gh + 241);

    auto g_yyyz_xyyyz_0 = pbuffer.data(idx_eri_0_gh + 242);

    auto g_yyyz_xyyzz_0 = pbuffer.data(idx_eri_0_gh + 243);

    auto g_yyyz_xyzzz_0 = pbuffer.data(idx_eri_0_gh + 244);

    auto g_yyyz_xzzzz_0 = pbuffer.data(idx_eri_0_gh + 245);

    auto g_yyyz_yyyyy_0 = pbuffer.data(idx_eri_0_gh + 246);

    auto g_yyyz_yyyyz_0 = pbuffer.data(idx_eri_0_gh + 247);

    auto g_yyyz_yyyzz_0 = pbuffer.data(idx_eri_0_gh + 248);

    auto g_yyyz_yyzzz_0 = pbuffer.data(idx_eri_0_gh + 249);

    auto g_yyyz_yzzzz_0 = pbuffer.data(idx_eri_0_gh + 250);

    auto g_yyyz_zzzzz_0 = pbuffer.data(idx_eri_0_gh + 251);

    #pragma omp simd aligned(g_yyy_xxxx_1, g_yyy_xxxxx_1, g_yyy_xxxxy_1, g_yyy_xxxxz_1, g_yyy_xxxy_1, g_yyy_xxxyy_1, g_yyy_xxxyz_1, g_yyy_xxxz_1, g_yyy_xxxzz_1, g_yyy_xxyy_1, g_yyy_xxyyy_1, g_yyy_xxyyz_1, g_yyy_xxyz_1, g_yyy_xxyzz_1, g_yyy_xxzz_1, g_yyy_xxzzz_1, g_yyy_xyyy_1, g_yyy_xyyyy_1, g_yyy_xyyyz_1, g_yyy_xyyz_1, g_yyy_xyyzz_1, g_yyy_xyzz_1, g_yyy_xyzzz_1, g_yyy_xzzz_1, g_yyy_xzzzz_1, g_yyy_yyyy_1, g_yyy_yyyyy_1, g_yyy_yyyyz_1, g_yyy_yyyz_1, g_yyy_yyyzz_1, g_yyy_yyzz_1, g_yyy_yyzzz_1, g_yyy_yzzz_1, g_yyy_yzzzz_1, g_yyy_zzzz_1, g_yyy_zzzzz_1, g_yyyz_xxxxx_0, g_yyyz_xxxxy_0, g_yyyz_xxxxz_0, g_yyyz_xxxyy_0, g_yyyz_xxxyz_0, g_yyyz_xxxzz_0, g_yyyz_xxyyy_0, g_yyyz_xxyyz_0, g_yyyz_xxyzz_0, g_yyyz_xxzzz_0, g_yyyz_xyyyy_0, g_yyyz_xyyyz_0, g_yyyz_xyyzz_0, g_yyyz_xyzzz_0, g_yyyz_xzzzz_0, g_yyyz_yyyyy_0, g_yyyz_yyyyz_0, g_yyyz_yyyzz_0, g_yyyz_yyzzz_0, g_yyyz_yzzzz_0, g_yyyz_zzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yyyz_xxxxx_0[i] = g_yyy_xxxxx_1[i] * pa_z[i];

        g_yyyz_xxxxy_0[i] = g_yyy_xxxxy_1[i] * pa_z[i];

        g_yyyz_xxxxz_0[i] = g_yyy_xxxx_1[i] * fe_0 + g_yyy_xxxxz_1[i] * pa_z[i];

        g_yyyz_xxxyy_0[i] = g_yyy_xxxyy_1[i] * pa_z[i];

        g_yyyz_xxxyz_0[i] = g_yyy_xxxy_1[i] * fe_0 + g_yyy_xxxyz_1[i] * pa_z[i];

        g_yyyz_xxxzz_0[i] = 2.0 * g_yyy_xxxz_1[i] * fe_0 + g_yyy_xxxzz_1[i] * pa_z[i];

        g_yyyz_xxyyy_0[i] = g_yyy_xxyyy_1[i] * pa_z[i];

        g_yyyz_xxyyz_0[i] = g_yyy_xxyy_1[i] * fe_0 + g_yyy_xxyyz_1[i] * pa_z[i];

        g_yyyz_xxyzz_0[i] = 2.0 * g_yyy_xxyz_1[i] * fe_0 + g_yyy_xxyzz_1[i] * pa_z[i];

        g_yyyz_xxzzz_0[i] = 3.0 * g_yyy_xxzz_1[i] * fe_0 + g_yyy_xxzzz_1[i] * pa_z[i];

        g_yyyz_xyyyy_0[i] = g_yyy_xyyyy_1[i] * pa_z[i];

        g_yyyz_xyyyz_0[i] = g_yyy_xyyy_1[i] * fe_0 + g_yyy_xyyyz_1[i] * pa_z[i];

        g_yyyz_xyyzz_0[i] = 2.0 * g_yyy_xyyz_1[i] * fe_0 + g_yyy_xyyzz_1[i] * pa_z[i];

        g_yyyz_xyzzz_0[i] = 3.0 * g_yyy_xyzz_1[i] * fe_0 + g_yyy_xyzzz_1[i] * pa_z[i];

        g_yyyz_xzzzz_0[i] = 4.0 * g_yyy_xzzz_1[i] * fe_0 + g_yyy_xzzzz_1[i] * pa_z[i];

        g_yyyz_yyyyy_0[i] = g_yyy_yyyyy_1[i] * pa_z[i];

        g_yyyz_yyyyz_0[i] = g_yyy_yyyy_1[i] * fe_0 + g_yyy_yyyyz_1[i] * pa_z[i];

        g_yyyz_yyyzz_0[i] = 2.0 * g_yyy_yyyz_1[i] * fe_0 + g_yyy_yyyzz_1[i] * pa_z[i];

        g_yyyz_yyzzz_0[i] = 3.0 * g_yyy_yyzz_1[i] * fe_0 + g_yyy_yyzzz_1[i] * pa_z[i];

        g_yyyz_yzzzz_0[i] = 4.0 * g_yyy_yzzz_1[i] * fe_0 + g_yyy_yzzzz_1[i] * pa_z[i];

        g_yyyz_zzzzz_0[i] = 5.0 * g_yyy_zzzz_1[i] * fe_0 + g_yyy_zzzzz_1[i] * pa_z[i];
    }

    // Set up 252-273 components of targeted buffer : GH

    auto g_yyzz_xxxxx_0 = pbuffer.data(idx_eri_0_gh + 252);

    auto g_yyzz_xxxxy_0 = pbuffer.data(idx_eri_0_gh + 253);

    auto g_yyzz_xxxxz_0 = pbuffer.data(idx_eri_0_gh + 254);

    auto g_yyzz_xxxyy_0 = pbuffer.data(idx_eri_0_gh + 255);

    auto g_yyzz_xxxyz_0 = pbuffer.data(idx_eri_0_gh + 256);

    auto g_yyzz_xxxzz_0 = pbuffer.data(idx_eri_0_gh + 257);

    auto g_yyzz_xxyyy_0 = pbuffer.data(idx_eri_0_gh + 258);

    auto g_yyzz_xxyyz_0 = pbuffer.data(idx_eri_0_gh + 259);

    auto g_yyzz_xxyzz_0 = pbuffer.data(idx_eri_0_gh + 260);

    auto g_yyzz_xxzzz_0 = pbuffer.data(idx_eri_0_gh + 261);

    auto g_yyzz_xyyyy_0 = pbuffer.data(idx_eri_0_gh + 262);

    auto g_yyzz_xyyyz_0 = pbuffer.data(idx_eri_0_gh + 263);

    auto g_yyzz_xyyzz_0 = pbuffer.data(idx_eri_0_gh + 264);

    auto g_yyzz_xyzzz_0 = pbuffer.data(idx_eri_0_gh + 265);

    auto g_yyzz_xzzzz_0 = pbuffer.data(idx_eri_0_gh + 266);

    auto g_yyzz_yyyyy_0 = pbuffer.data(idx_eri_0_gh + 267);

    auto g_yyzz_yyyyz_0 = pbuffer.data(idx_eri_0_gh + 268);

    auto g_yyzz_yyyzz_0 = pbuffer.data(idx_eri_0_gh + 269);

    auto g_yyzz_yyzzz_0 = pbuffer.data(idx_eri_0_gh + 270);

    auto g_yyzz_yzzzz_0 = pbuffer.data(idx_eri_0_gh + 271);

    auto g_yyzz_zzzzz_0 = pbuffer.data(idx_eri_0_gh + 272);

    #pragma omp simd aligned(g_yy_xxxxy_0, g_yy_xxxxy_1, g_yy_xxxyy_0, g_yy_xxxyy_1, g_yy_xxyyy_0, g_yy_xxyyy_1, g_yy_xyyyy_0, g_yy_xyyyy_1, g_yy_yyyyy_0, g_yy_yyyyy_1, g_yyz_xxxxy_1, g_yyz_xxxyy_1, g_yyz_xxyyy_1, g_yyz_xyyyy_1, g_yyz_yyyyy_1, g_yyzz_xxxxx_0, g_yyzz_xxxxy_0, g_yyzz_xxxxz_0, g_yyzz_xxxyy_0, g_yyzz_xxxyz_0, g_yyzz_xxxzz_0, g_yyzz_xxyyy_0, g_yyzz_xxyyz_0, g_yyzz_xxyzz_0, g_yyzz_xxzzz_0, g_yyzz_xyyyy_0, g_yyzz_xyyyz_0, g_yyzz_xyyzz_0, g_yyzz_xyzzz_0, g_yyzz_xzzzz_0, g_yyzz_yyyyy_0, g_yyzz_yyyyz_0, g_yyzz_yyyzz_0, g_yyzz_yyzzz_0, g_yyzz_yzzzz_0, g_yyzz_zzzzz_0, g_yzz_xxxxx_1, g_yzz_xxxxz_1, g_yzz_xxxyz_1, g_yzz_xxxz_1, g_yzz_xxxzz_1, g_yzz_xxyyz_1, g_yzz_xxyz_1, g_yzz_xxyzz_1, g_yzz_xxzz_1, g_yzz_xxzzz_1, g_yzz_xyyyz_1, g_yzz_xyyz_1, g_yzz_xyyzz_1, g_yzz_xyzz_1, g_yzz_xyzzz_1, g_yzz_xzzz_1, g_yzz_xzzzz_1, g_yzz_yyyyz_1, g_yzz_yyyz_1, g_yzz_yyyzz_1, g_yzz_yyzz_1, g_yzz_yyzzz_1, g_yzz_yzzz_1, g_yzz_yzzzz_1, g_yzz_zzzz_1, g_yzz_zzzzz_1, g_zz_xxxxx_0, g_zz_xxxxx_1, g_zz_xxxxz_0, g_zz_xxxxz_1, g_zz_xxxyz_0, g_zz_xxxyz_1, g_zz_xxxzz_0, g_zz_xxxzz_1, g_zz_xxyyz_0, g_zz_xxyyz_1, g_zz_xxyzz_0, g_zz_xxyzz_1, g_zz_xxzzz_0, g_zz_xxzzz_1, g_zz_xyyyz_0, g_zz_xyyyz_1, g_zz_xyyzz_0, g_zz_xyyzz_1, g_zz_xyzzz_0, g_zz_xyzzz_1, g_zz_xzzzz_0, g_zz_xzzzz_1, g_zz_yyyyz_0, g_zz_yyyyz_1, g_zz_yyyzz_0, g_zz_yyyzz_1, g_zz_yyzzz_0, g_zz_yyzzz_1, g_zz_yzzzz_0, g_zz_yzzzz_1, g_zz_zzzzz_0, g_zz_zzzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyzz_xxxxx_0[i] = g_zz_xxxxx_0[i] * fbe_0 - g_zz_xxxxx_1[i] * fz_be_0 + g_yzz_xxxxx_1[i] * pa_y[i];

        g_yyzz_xxxxy_0[i] = g_yy_xxxxy_0[i] * fbe_0 - g_yy_xxxxy_1[i] * fz_be_0 + g_yyz_xxxxy_1[i] * pa_z[i];

        g_yyzz_xxxxz_0[i] = g_zz_xxxxz_0[i] * fbe_0 - g_zz_xxxxz_1[i] * fz_be_0 + g_yzz_xxxxz_1[i] * pa_y[i];

        g_yyzz_xxxyy_0[i] = g_yy_xxxyy_0[i] * fbe_0 - g_yy_xxxyy_1[i] * fz_be_0 + g_yyz_xxxyy_1[i] * pa_z[i];

        g_yyzz_xxxyz_0[i] = g_zz_xxxyz_0[i] * fbe_0 - g_zz_xxxyz_1[i] * fz_be_0 + g_yzz_xxxz_1[i] * fe_0 + g_yzz_xxxyz_1[i] * pa_y[i];

        g_yyzz_xxxzz_0[i] = g_zz_xxxzz_0[i] * fbe_0 - g_zz_xxxzz_1[i] * fz_be_0 + g_yzz_xxxzz_1[i] * pa_y[i];

        g_yyzz_xxyyy_0[i] = g_yy_xxyyy_0[i] * fbe_0 - g_yy_xxyyy_1[i] * fz_be_0 + g_yyz_xxyyy_1[i] * pa_z[i];

        g_yyzz_xxyyz_0[i] = g_zz_xxyyz_0[i] * fbe_0 - g_zz_xxyyz_1[i] * fz_be_0 + 2.0 * g_yzz_xxyz_1[i] * fe_0 + g_yzz_xxyyz_1[i] * pa_y[i];

        g_yyzz_xxyzz_0[i] = g_zz_xxyzz_0[i] * fbe_0 - g_zz_xxyzz_1[i] * fz_be_0 + g_yzz_xxzz_1[i] * fe_0 + g_yzz_xxyzz_1[i] * pa_y[i];

        g_yyzz_xxzzz_0[i] = g_zz_xxzzz_0[i] * fbe_0 - g_zz_xxzzz_1[i] * fz_be_0 + g_yzz_xxzzz_1[i] * pa_y[i];

        g_yyzz_xyyyy_0[i] = g_yy_xyyyy_0[i] * fbe_0 - g_yy_xyyyy_1[i] * fz_be_0 + g_yyz_xyyyy_1[i] * pa_z[i];

        g_yyzz_xyyyz_0[i] = g_zz_xyyyz_0[i] * fbe_0 - g_zz_xyyyz_1[i] * fz_be_0 + 3.0 * g_yzz_xyyz_1[i] * fe_0 + g_yzz_xyyyz_1[i] * pa_y[i];

        g_yyzz_xyyzz_0[i] = g_zz_xyyzz_0[i] * fbe_0 - g_zz_xyyzz_1[i] * fz_be_0 + 2.0 * g_yzz_xyzz_1[i] * fe_0 + g_yzz_xyyzz_1[i] * pa_y[i];

        g_yyzz_xyzzz_0[i] = g_zz_xyzzz_0[i] * fbe_0 - g_zz_xyzzz_1[i] * fz_be_0 + g_yzz_xzzz_1[i] * fe_0 + g_yzz_xyzzz_1[i] * pa_y[i];

        g_yyzz_xzzzz_0[i] = g_zz_xzzzz_0[i] * fbe_0 - g_zz_xzzzz_1[i] * fz_be_0 + g_yzz_xzzzz_1[i] * pa_y[i];

        g_yyzz_yyyyy_0[i] = g_yy_yyyyy_0[i] * fbe_0 - g_yy_yyyyy_1[i] * fz_be_0 + g_yyz_yyyyy_1[i] * pa_z[i];

        g_yyzz_yyyyz_0[i] = g_zz_yyyyz_0[i] * fbe_0 - g_zz_yyyyz_1[i] * fz_be_0 + 4.0 * g_yzz_yyyz_1[i] * fe_0 + g_yzz_yyyyz_1[i] * pa_y[i];

        g_yyzz_yyyzz_0[i] = g_zz_yyyzz_0[i] * fbe_0 - g_zz_yyyzz_1[i] * fz_be_0 + 3.0 * g_yzz_yyzz_1[i] * fe_0 + g_yzz_yyyzz_1[i] * pa_y[i];

        g_yyzz_yyzzz_0[i] = g_zz_yyzzz_0[i] * fbe_0 - g_zz_yyzzz_1[i] * fz_be_0 + 2.0 * g_yzz_yzzz_1[i] * fe_0 + g_yzz_yyzzz_1[i] * pa_y[i];

        g_yyzz_yzzzz_0[i] = g_zz_yzzzz_0[i] * fbe_0 - g_zz_yzzzz_1[i] * fz_be_0 + g_yzz_zzzz_1[i] * fe_0 + g_yzz_yzzzz_1[i] * pa_y[i];

        g_yyzz_zzzzz_0[i] = g_zz_zzzzz_0[i] * fbe_0 - g_zz_zzzzz_1[i] * fz_be_0 + g_yzz_zzzzz_1[i] * pa_y[i];
    }

    // Set up 273-294 components of targeted buffer : GH

    auto g_yzzz_xxxxx_0 = pbuffer.data(idx_eri_0_gh + 273);

    auto g_yzzz_xxxxy_0 = pbuffer.data(idx_eri_0_gh + 274);

    auto g_yzzz_xxxxz_0 = pbuffer.data(idx_eri_0_gh + 275);

    auto g_yzzz_xxxyy_0 = pbuffer.data(idx_eri_0_gh + 276);

    auto g_yzzz_xxxyz_0 = pbuffer.data(idx_eri_0_gh + 277);

    auto g_yzzz_xxxzz_0 = pbuffer.data(idx_eri_0_gh + 278);

    auto g_yzzz_xxyyy_0 = pbuffer.data(idx_eri_0_gh + 279);

    auto g_yzzz_xxyyz_0 = pbuffer.data(idx_eri_0_gh + 280);

    auto g_yzzz_xxyzz_0 = pbuffer.data(idx_eri_0_gh + 281);

    auto g_yzzz_xxzzz_0 = pbuffer.data(idx_eri_0_gh + 282);

    auto g_yzzz_xyyyy_0 = pbuffer.data(idx_eri_0_gh + 283);

    auto g_yzzz_xyyyz_0 = pbuffer.data(idx_eri_0_gh + 284);

    auto g_yzzz_xyyzz_0 = pbuffer.data(idx_eri_0_gh + 285);

    auto g_yzzz_xyzzz_0 = pbuffer.data(idx_eri_0_gh + 286);

    auto g_yzzz_xzzzz_0 = pbuffer.data(idx_eri_0_gh + 287);

    auto g_yzzz_yyyyy_0 = pbuffer.data(idx_eri_0_gh + 288);

    auto g_yzzz_yyyyz_0 = pbuffer.data(idx_eri_0_gh + 289);

    auto g_yzzz_yyyzz_0 = pbuffer.data(idx_eri_0_gh + 290);

    auto g_yzzz_yyzzz_0 = pbuffer.data(idx_eri_0_gh + 291);

    auto g_yzzz_yzzzz_0 = pbuffer.data(idx_eri_0_gh + 292);

    auto g_yzzz_zzzzz_0 = pbuffer.data(idx_eri_0_gh + 293);

    #pragma omp simd aligned(g_yzzz_xxxxx_0, g_yzzz_xxxxy_0, g_yzzz_xxxxz_0, g_yzzz_xxxyy_0, g_yzzz_xxxyz_0, g_yzzz_xxxzz_0, g_yzzz_xxyyy_0, g_yzzz_xxyyz_0, g_yzzz_xxyzz_0, g_yzzz_xxzzz_0, g_yzzz_xyyyy_0, g_yzzz_xyyyz_0, g_yzzz_xyyzz_0, g_yzzz_xyzzz_0, g_yzzz_xzzzz_0, g_yzzz_yyyyy_0, g_yzzz_yyyyz_0, g_yzzz_yyyzz_0, g_yzzz_yyzzz_0, g_yzzz_yzzzz_0, g_yzzz_zzzzz_0, g_zzz_xxxx_1, g_zzz_xxxxx_1, g_zzz_xxxxy_1, g_zzz_xxxxz_1, g_zzz_xxxy_1, g_zzz_xxxyy_1, g_zzz_xxxyz_1, g_zzz_xxxz_1, g_zzz_xxxzz_1, g_zzz_xxyy_1, g_zzz_xxyyy_1, g_zzz_xxyyz_1, g_zzz_xxyz_1, g_zzz_xxyzz_1, g_zzz_xxzz_1, g_zzz_xxzzz_1, g_zzz_xyyy_1, g_zzz_xyyyy_1, g_zzz_xyyyz_1, g_zzz_xyyz_1, g_zzz_xyyzz_1, g_zzz_xyzz_1, g_zzz_xyzzz_1, g_zzz_xzzz_1, g_zzz_xzzzz_1, g_zzz_yyyy_1, g_zzz_yyyyy_1, g_zzz_yyyyz_1, g_zzz_yyyz_1, g_zzz_yyyzz_1, g_zzz_yyzz_1, g_zzz_yyzzz_1, g_zzz_yzzz_1, g_zzz_yzzzz_1, g_zzz_zzzz_1, g_zzz_zzzzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yzzz_xxxxx_0[i] = g_zzz_xxxxx_1[i] * pa_y[i];

        g_yzzz_xxxxy_0[i] = g_zzz_xxxx_1[i] * fe_0 + g_zzz_xxxxy_1[i] * pa_y[i];

        g_yzzz_xxxxz_0[i] = g_zzz_xxxxz_1[i] * pa_y[i];

        g_yzzz_xxxyy_0[i] = 2.0 * g_zzz_xxxy_1[i] * fe_0 + g_zzz_xxxyy_1[i] * pa_y[i];

        g_yzzz_xxxyz_0[i] = g_zzz_xxxz_1[i] * fe_0 + g_zzz_xxxyz_1[i] * pa_y[i];

        g_yzzz_xxxzz_0[i] = g_zzz_xxxzz_1[i] * pa_y[i];

        g_yzzz_xxyyy_0[i] = 3.0 * g_zzz_xxyy_1[i] * fe_0 + g_zzz_xxyyy_1[i] * pa_y[i];

        g_yzzz_xxyyz_0[i] = 2.0 * g_zzz_xxyz_1[i] * fe_0 + g_zzz_xxyyz_1[i] * pa_y[i];

        g_yzzz_xxyzz_0[i] = g_zzz_xxzz_1[i] * fe_0 + g_zzz_xxyzz_1[i] * pa_y[i];

        g_yzzz_xxzzz_0[i] = g_zzz_xxzzz_1[i] * pa_y[i];

        g_yzzz_xyyyy_0[i] = 4.0 * g_zzz_xyyy_1[i] * fe_0 + g_zzz_xyyyy_1[i] * pa_y[i];

        g_yzzz_xyyyz_0[i] = 3.0 * g_zzz_xyyz_1[i] * fe_0 + g_zzz_xyyyz_1[i] * pa_y[i];

        g_yzzz_xyyzz_0[i] = 2.0 * g_zzz_xyzz_1[i] * fe_0 + g_zzz_xyyzz_1[i] * pa_y[i];

        g_yzzz_xyzzz_0[i] = g_zzz_xzzz_1[i] * fe_0 + g_zzz_xyzzz_1[i] * pa_y[i];

        g_yzzz_xzzzz_0[i] = g_zzz_xzzzz_1[i] * pa_y[i];

        g_yzzz_yyyyy_0[i] = 5.0 * g_zzz_yyyy_1[i] * fe_0 + g_zzz_yyyyy_1[i] * pa_y[i];

        g_yzzz_yyyyz_0[i] = 4.0 * g_zzz_yyyz_1[i] * fe_0 + g_zzz_yyyyz_1[i] * pa_y[i];

        g_yzzz_yyyzz_0[i] = 3.0 * g_zzz_yyzz_1[i] * fe_0 + g_zzz_yyyzz_1[i] * pa_y[i];

        g_yzzz_yyzzz_0[i] = 2.0 * g_zzz_yzzz_1[i] * fe_0 + g_zzz_yyzzz_1[i] * pa_y[i];

        g_yzzz_yzzzz_0[i] = g_zzz_zzzz_1[i] * fe_0 + g_zzz_yzzzz_1[i] * pa_y[i];

        g_yzzz_zzzzz_0[i] = g_zzz_zzzzz_1[i] * pa_y[i];
    }

    // Set up 294-315 components of targeted buffer : GH

    auto g_zzzz_xxxxx_0 = pbuffer.data(idx_eri_0_gh + 294);

    auto g_zzzz_xxxxy_0 = pbuffer.data(idx_eri_0_gh + 295);

    auto g_zzzz_xxxxz_0 = pbuffer.data(idx_eri_0_gh + 296);

    auto g_zzzz_xxxyy_0 = pbuffer.data(idx_eri_0_gh + 297);

    auto g_zzzz_xxxyz_0 = pbuffer.data(idx_eri_0_gh + 298);

    auto g_zzzz_xxxzz_0 = pbuffer.data(idx_eri_0_gh + 299);

    auto g_zzzz_xxyyy_0 = pbuffer.data(idx_eri_0_gh + 300);

    auto g_zzzz_xxyyz_0 = pbuffer.data(idx_eri_0_gh + 301);

    auto g_zzzz_xxyzz_0 = pbuffer.data(idx_eri_0_gh + 302);

    auto g_zzzz_xxzzz_0 = pbuffer.data(idx_eri_0_gh + 303);

    auto g_zzzz_xyyyy_0 = pbuffer.data(idx_eri_0_gh + 304);

    auto g_zzzz_xyyyz_0 = pbuffer.data(idx_eri_0_gh + 305);

    auto g_zzzz_xyyzz_0 = pbuffer.data(idx_eri_0_gh + 306);

    auto g_zzzz_xyzzz_0 = pbuffer.data(idx_eri_0_gh + 307);

    auto g_zzzz_xzzzz_0 = pbuffer.data(idx_eri_0_gh + 308);

    auto g_zzzz_yyyyy_0 = pbuffer.data(idx_eri_0_gh + 309);

    auto g_zzzz_yyyyz_0 = pbuffer.data(idx_eri_0_gh + 310);

    auto g_zzzz_yyyzz_0 = pbuffer.data(idx_eri_0_gh + 311);

    auto g_zzzz_yyzzz_0 = pbuffer.data(idx_eri_0_gh + 312);

    auto g_zzzz_yzzzz_0 = pbuffer.data(idx_eri_0_gh + 313);

    auto g_zzzz_zzzzz_0 = pbuffer.data(idx_eri_0_gh + 314);

    #pragma omp simd aligned(g_zz_xxxxx_0, g_zz_xxxxx_1, g_zz_xxxxy_0, g_zz_xxxxy_1, g_zz_xxxxz_0, g_zz_xxxxz_1, g_zz_xxxyy_0, g_zz_xxxyy_1, g_zz_xxxyz_0, g_zz_xxxyz_1, g_zz_xxxzz_0, g_zz_xxxzz_1, g_zz_xxyyy_0, g_zz_xxyyy_1, g_zz_xxyyz_0, g_zz_xxyyz_1, g_zz_xxyzz_0, g_zz_xxyzz_1, g_zz_xxzzz_0, g_zz_xxzzz_1, g_zz_xyyyy_0, g_zz_xyyyy_1, g_zz_xyyyz_0, g_zz_xyyyz_1, g_zz_xyyzz_0, g_zz_xyyzz_1, g_zz_xyzzz_0, g_zz_xyzzz_1, g_zz_xzzzz_0, g_zz_xzzzz_1, g_zz_yyyyy_0, g_zz_yyyyy_1, g_zz_yyyyz_0, g_zz_yyyyz_1, g_zz_yyyzz_0, g_zz_yyyzz_1, g_zz_yyzzz_0, g_zz_yyzzz_1, g_zz_yzzzz_0, g_zz_yzzzz_1, g_zz_zzzzz_0, g_zz_zzzzz_1, g_zzz_xxxx_1, g_zzz_xxxxx_1, g_zzz_xxxxy_1, g_zzz_xxxxz_1, g_zzz_xxxy_1, g_zzz_xxxyy_1, g_zzz_xxxyz_1, g_zzz_xxxz_1, g_zzz_xxxzz_1, g_zzz_xxyy_1, g_zzz_xxyyy_1, g_zzz_xxyyz_1, g_zzz_xxyz_1, g_zzz_xxyzz_1, g_zzz_xxzz_1, g_zzz_xxzzz_1, g_zzz_xyyy_1, g_zzz_xyyyy_1, g_zzz_xyyyz_1, g_zzz_xyyz_1, g_zzz_xyyzz_1, g_zzz_xyzz_1, g_zzz_xyzzz_1, g_zzz_xzzz_1, g_zzz_xzzzz_1, g_zzz_yyyy_1, g_zzz_yyyyy_1, g_zzz_yyyyz_1, g_zzz_yyyz_1, g_zzz_yyyzz_1, g_zzz_yyzz_1, g_zzz_yyzzz_1, g_zzz_yzzz_1, g_zzz_yzzzz_1, g_zzz_zzzz_1, g_zzz_zzzzz_1, g_zzzz_xxxxx_0, g_zzzz_xxxxy_0, g_zzzz_xxxxz_0, g_zzzz_xxxyy_0, g_zzzz_xxxyz_0, g_zzzz_xxxzz_0, g_zzzz_xxyyy_0, g_zzzz_xxyyz_0, g_zzzz_xxyzz_0, g_zzzz_xxzzz_0, g_zzzz_xyyyy_0, g_zzzz_xyyyz_0, g_zzzz_xyyzz_0, g_zzzz_xyzzz_0, g_zzzz_xzzzz_0, g_zzzz_yyyyy_0, g_zzzz_yyyyz_0, g_zzzz_yyyzz_0, g_zzzz_yyzzz_0, g_zzzz_yzzzz_0, g_zzzz_zzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_zzzz_xxxxx_0[i] = 3.0 * g_zz_xxxxx_0[i] * fbe_0 - 3.0 * g_zz_xxxxx_1[i] * fz_be_0 + g_zzz_xxxxx_1[i] * pa_z[i];

        g_zzzz_xxxxy_0[i] = 3.0 * g_zz_xxxxy_0[i] * fbe_0 - 3.0 * g_zz_xxxxy_1[i] * fz_be_0 + g_zzz_xxxxy_1[i] * pa_z[i];

        g_zzzz_xxxxz_0[i] = 3.0 * g_zz_xxxxz_0[i] * fbe_0 - 3.0 * g_zz_xxxxz_1[i] * fz_be_0 + g_zzz_xxxx_1[i] * fe_0 + g_zzz_xxxxz_1[i] * pa_z[i];

        g_zzzz_xxxyy_0[i] = 3.0 * g_zz_xxxyy_0[i] * fbe_0 - 3.0 * g_zz_xxxyy_1[i] * fz_be_0 + g_zzz_xxxyy_1[i] * pa_z[i];

        g_zzzz_xxxyz_0[i] = 3.0 * g_zz_xxxyz_0[i] * fbe_0 - 3.0 * g_zz_xxxyz_1[i] * fz_be_0 + g_zzz_xxxy_1[i] * fe_0 + g_zzz_xxxyz_1[i] * pa_z[i];

        g_zzzz_xxxzz_0[i] = 3.0 * g_zz_xxxzz_0[i] * fbe_0 - 3.0 * g_zz_xxxzz_1[i] * fz_be_0 + 2.0 * g_zzz_xxxz_1[i] * fe_0 + g_zzz_xxxzz_1[i] * pa_z[i];

        g_zzzz_xxyyy_0[i] = 3.0 * g_zz_xxyyy_0[i] * fbe_0 - 3.0 * g_zz_xxyyy_1[i] * fz_be_0 + g_zzz_xxyyy_1[i] * pa_z[i];

        g_zzzz_xxyyz_0[i] = 3.0 * g_zz_xxyyz_0[i] * fbe_0 - 3.0 * g_zz_xxyyz_1[i] * fz_be_0 + g_zzz_xxyy_1[i] * fe_0 + g_zzz_xxyyz_1[i] * pa_z[i];

        g_zzzz_xxyzz_0[i] = 3.0 * g_zz_xxyzz_0[i] * fbe_0 - 3.0 * g_zz_xxyzz_1[i] * fz_be_0 + 2.0 * g_zzz_xxyz_1[i] * fe_0 + g_zzz_xxyzz_1[i] * pa_z[i];

        g_zzzz_xxzzz_0[i] = 3.0 * g_zz_xxzzz_0[i] * fbe_0 - 3.0 * g_zz_xxzzz_1[i] * fz_be_0 + 3.0 * g_zzz_xxzz_1[i] * fe_0 + g_zzz_xxzzz_1[i] * pa_z[i];

        g_zzzz_xyyyy_0[i] = 3.0 * g_zz_xyyyy_0[i] * fbe_0 - 3.0 * g_zz_xyyyy_1[i] * fz_be_0 + g_zzz_xyyyy_1[i] * pa_z[i];

        g_zzzz_xyyyz_0[i] = 3.0 * g_zz_xyyyz_0[i] * fbe_0 - 3.0 * g_zz_xyyyz_1[i] * fz_be_0 + g_zzz_xyyy_1[i] * fe_0 + g_zzz_xyyyz_1[i] * pa_z[i];

        g_zzzz_xyyzz_0[i] = 3.0 * g_zz_xyyzz_0[i] * fbe_0 - 3.0 * g_zz_xyyzz_1[i] * fz_be_0 + 2.0 * g_zzz_xyyz_1[i] * fe_0 + g_zzz_xyyzz_1[i] * pa_z[i];

        g_zzzz_xyzzz_0[i] = 3.0 * g_zz_xyzzz_0[i] * fbe_0 - 3.0 * g_zz_xyzzz_1[i] * fz_be_0 + 3.0 * g_zzz_xyzz_1[i] * fe_0 + g_zzz_xyzzz_1[i] * pa_z[i];

        g_zzzz_xzzzz_0[i] = 3.0 * g_zz_xzzzz_0[i] * fbe_0 - 3.0 * g_zz_xzzzz_1[i] * fz_be_0 + 4.0 * g_zzz_xzzz_1[i] * fe_0 + g_zzz_xzzzz_1[i] * pa_z[i];

        g_zzzz_yyyyy_0[i] = 3.0 * g_zz_yyyyy_0[i] * fbe_0 - 3.0 * g_zz_yyyyy_1[i] * fz_be_0 + g_zzz_yyyyy_1[i] * pa_z[i];

        g_zzzz_yyyyz_0[i] = 3.0 * g_zz_yyyyz_0[i] * fbe_0 - 3.0 * g_zz_yyyyz_1[i] * fz_be_0 + g_zzz_yyyy_1[i] * fe_0 + g_zzz_yyyyz_1[i] * pa_z[i];

        g_zzzz_yyyzz_0[i] = 3.0 * g_zz_yyyzz_0[i] * fbe_0 - 3.0 * g_zz_yyyzz_1[i] * fz_be_0 + 2.0 * g_zzz_yyyz_1[i] * fe_0 + g_zzz_yyyzz_1[i] * pa_z[i];

        g_zzzz_yyzzz_0[i] = 3.0 * g_zz_yyzzz_0[i] * fbe_0 - 3.0 * g_zz_yyzzz_1[i] * fz_be_0 + 3.0 * g_zzz_yyzz_1[i] * fe_0 + g_zzz_yyzzz_1[i] * pa_z[i];

        g_zzzz_yzzzz_0[i] = 3.0 * g_zz_yzzzz_0[i] * fbe_0 - 3.0 * g_zz_yzzzz_1[i] * fz_be_0 + 4.0 * g_zzz_yzzz_1[i] * fe_0 + g_zzz_yzzzz_1[i] * pa_z[i];

        g_zzzz_zzzzz_0[i] = 3.0 * g_zz_zzzzz_0[i] * fbe_0 - 3.0 * g_zz_zzzzz_1[i] * fz_be_0 + 5.0 * g_zzz_zzzz_1[i] * fe_0 + g_zzz_zzzzz_1[i] * pa_z[i];
    }

}

} // t2ceri namespace


#include "ThreeCenterElectronRepulsionPrimRecISG.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_isg(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_isg,
                                 size_t idx_eri_0_gsg,
                                 size_t idx_eri_1_gsg,
                                 size_t idx_eri_1_hsf,
                                 size_t idx_eri_1_hsg,
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

    /// Set up components of auxilary buffer : GSG

    auto g_xxxx_0_xxxx_0 = pbuffer.data(idx_eri_0_gsg);

    auto g_xxxx_0_xxxy_0 = pbuffer.data(idx_eri_0_gsg + 1);

    auto g_xxxx_0_xxxz_0 = pbuffer.data(idx_eri_0_gsg + 2);

    auto g_xxxx_0_xxyy_0 = pbuffer.data(idx_eri_0_gsg + 3);

    auto g_xxxx_0_xxyz_0 = pbuffer.data(idx_eri_0_gsg + 4);

    auto g_xxxx_0_xxzz_0 = pbuffer.data(idx_eri_0_gsg + 5);

    auto g_xxxx_0_xyyy_0 = pbuffer.data(idx_eri_0_gsg + 6);

    auto g_xxxx_0_xyyz_0 = pbuffer.data(idx_eri_0_gsg + 7);

    auto g_xxxx_0_xyzz_0 = pbuffer.data(idx_eri_0_gsg + 8);

    auto g_xxxx_0_xzzz_0 = pbuffer.data(idx_eri_0_gsg + 9);

    auto g_xxxx_0_yyyy_0 = pbuffer.data(idx_eri_0_gsg + 10);

    auto g_xxxx_0_yyyz_0 = pbuffer.data(idx_eri_0_gsg + 11);

    auto g_xxxx_0_yyzz_0 = pbuffer.data(idx_eri_0_gsg + 12);

    auto g_xxxx_0_yzzz_0 = pbuffer.data(idx_eri_0_gsg + 13);

    auto g_xxxx_0_zzzz_0 = pbuffer.data(idx_eri_0_gsg + 14);

    auto g_xxxy_0_xxxx_0 = pbuffer.data(idx_eri_0_gsg + 15);

    auto g_xxxy_0_xxxz_0 = pbuffer.data(idx_eri_0_gsg + 17);

    auto g_xxxy_0_xxzz_0 = pbuffer.data(idx_eri_0_gsg + 20);

    auto g_xxxy_0_xzzz_0 = pbuffer.data(idx_eri_0_gsg + 24);

    auto g_xxxz_0_xxxx_0 = pbuffer.data(idx_eri_0_gsg + 30);

    auto g_xxxz_0_xxxy_0 = pbuffer.data(idx_eri_0_gsg + 31);

    auto g_xxxz_0_xxyy_0 = pbuffer.data(idx_eri_0_gsg + 33);

    auto g_xxxz_0_xyyy_0 = pbuffer.data(idx_eri_0_gsg + 36);

    auto g_xxyy_0_xxxx_0 = pbuffer.data(idx_eri_0_gsg + 45);

    auto g_xxyy_0_xxxy_0 = pbuffer.data(idx_eri_0_gsg + 46);

    auto g_xxyy_0_xxxz_0 = pbuffer.data(idx_eri_0_gsg + 47);

    auto g_xxyy_0_xxyy_0 = pbuffer.data(idx_eri_0_gsg + 48);

    auto g_xxyy_0_xxyz_0 = pbuffer.data(idx_eri_0_gsg + 49);

    auto g_xxyy_0_xxzz_0 = pbuffer.data(idx_eri_0_gsg + 50);

    auto g_xxyy_0_xyyy_0 = pbuffer.data(idx_eri_0_gsg + 51);

    auto g_xxyy_0_xyyz_0 = pbuffer.data(idx_eri_0_gsg + 52);

    auto g_xxyy_0_xyzz_0 = pbuffer.data(idx_eri_0_gsg + 53);

    auto g_xxyy_0_xzzz_0 = pbuffer.data(idx_eri_0_gsg + 54);

    auto g_xxyy_0_yyyy_0 = pbuffer.data(idx_eri_0_gsg + 55);

    auto g_xxyy_0_yyyz_0 = pbuffer.data(idx_eri_0_gsg + 56);

    auto g_xxyy_0_yyzz_0 = pbuffer.data(idx_eri_0_gsg + 57);

    auto g_xxyy_0_yzzz_0 = pbuffer.data(idx_eri_0_gsg + 58);

    auto g_xxyy_0_zzzz_0 = pbuffer.data(idx_eri_0_gsg + 59);

    auto g_xxzz_0_xxxx_0 = pbuffer.data(idx_eri_0_gsg + 75);

    auto g_xxzz_0_xxxy_0 = pbuffer.data(idx_eri_0_gsg + 76);

    auto g_xxzz_0_xxxz_0 = pbuffer.data(idx_eri_0_gsg + 77);

    auto g_xxzz_0_xxyy_0 = pbuffer.data(idx_eri_0_gsg + 78);

    auto g_xxzz_0_xxyz_0 = pbuffer.data(idx_eri_0_gsg + 79);

    auto g_xxzz_0_xxzz_0 = pbuffer.data(idx_eri_0_gsg + 80);

    auto g_xxzz_0_xyyy_0 = pbuffer.data(idx_eri_0_gsg + 81);

    auto g_xxzz_0_xyyz_0 = pbuffer.data(idx_eri_0_gsg + 82);

    auto g_xxzz_0_xyzz_0 = pbuffer.data(idx_eri_0_gsg + 83);

    auto g_xxzz_0_xzzz_0 = pbuffer.data(idx_eri_0_gsg + 84);

    auto g_xxzz_0_yyyy_0 = pbuffer.data(idx_eri_0_gsg + 85);

    auto g_xxzz_0_yyyz_0 = pbuffer.data(idx_eri_0_gsg + 86);

    auto g_xxzz_0_yyzz_0 = pbuffer.data(idx_eri_0_gsg + 87);

    auto g_xxzz_0_yzzz_0 = pbuffer.data(idx_eri_0_gsg + 88);

    auto g_xxzz_0_zzzz_0 = pbuffer.data(idx_eri_0_gsg + 89);

    auto g_xyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_gsg + 91);

    auto g_xyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_gsg + 93);

    auto g_xyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_gsg + 94);

    auto g_xyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_gsg + 96);

    auto g_xyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_gsg + 97);

    auto g_xyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_gsg + 98);

    auto g_xyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_gsg + 100);

    auto g_xyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_gsg + 101);

    auto g_xyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_gsg + 102);

    auto g_xyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_gsg + 103);

    auto g_xyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_gsg + 104);

    auto g_xzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_gsg + 137);

    auto g_xzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_gsg + 139);

    auto g_xzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_gsg + 140);

    auto g_xzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_gsg + 142);

    auto g_xzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_gsg + 143);

    auto g_xzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_gsg + 144);

    auto g_xzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_gsg + 145);

    auto g_xzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_gsg + 146);

    auto g_xzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_gsg + 147);

    auto g_xzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_gsg + 148);

    auto g_xzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_gsg + 149);

    auto g_yyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_gsg + 150);

    auto g_yyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_gsg + 151);

    auto g_yyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_gsg + 152);

    auto g_yyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_gsg + 153);

    auto g_yyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_gsg + 154);

    auto g_yyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_gsg + 155);

    auto g_yyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_gsg + 156);

    auto g_yyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_gsg + 157);

    auto g_yyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_gsg + 158);

    auto g_yyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_gsg + 159);

    auto g_yyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_gsg + 160);

    auto g_yyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_gsg + 161);

    auto g_yyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_gsg + 162);

    auto g_yyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_gsg + 163);

    auto g_yyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_gsg + 164);

    auto g_yyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_gsg + 166);

    auto g_yyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_gsg + 168);

    auto g_yyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_gsg + 171);

    auto g_yyyz_0_yyyy_0 = pbuffer.data(idx_eri_0_gsg + 175);

    auto g_yyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_gsg + 180);

    auto g_yyzz_0_xxxy_0 = pbuffer.data(idx_eri_0_gsg + 181);

    auto g_yyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_gsg + 182);

    auto g_yyzz_0_xxyy_0 = pbuffer.data(idx_eri_0_gsg + 183);

    auto g_yyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_gsg + 184);

    auto g_yyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_gsg + 185);

    auto g_yyzz_0_xyyy_0 = pbuffer.data(idx_eri_0_gsg + 186);

    auto g_yyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_gsg + 187);

    auto g_yyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_gsg + 188);

    auto g_yyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_gsg + 189);

    auto g_yyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_gsg + 190);

    auto g_yyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_gsg + 191);

    auto g_yyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_gsg + 192);

    auto g_yyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_gsg + 193);

    auto g_yyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_gsg + 194);

    auto g_yzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_gsg + 195);

    auto g_yzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_gsg + 197);

    auto g_yzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_gsg + 199);

    auto g_yzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_gsg + 200);

    auto g_yzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_gsg + 202);

    auto g_yzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_gsg + 203);

    auto g_yzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_gsg + 204);

    auto g_yzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_gsg + 206);

    auto g_yzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_gsg + 207);

    auto g_yzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_gsg + 208);

    auto g_yzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_gsg + 209);

    auto g_zzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_gsg + 210);

    auto g_zzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_gsg + 211);

    auto g_zzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_gsg + 212);

    auto g_zzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_gsg + 213);

    auto g_zzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_gsg + 214);

    auto g_zzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_gsg + 215);

    auto g_zzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_gsg + 216);

    auto g_zzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_gsg + 217);

    auto g_zzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_gsg + 218);

    auto g_zzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_gsg + 219);

    auto g_zzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_gsg + 220);

    auto g_zzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_gsg + 221);

    auto g_zzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_gsg + 222);

    auto g_zzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_gsg + 223);

    auto g_zzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_gsg + 224);

    /// Set up components of auxilary buffer : GSG

    auto g_xxxx_0_xxxx_1 = pbuffer.data(idx_eri_1_gsg);

    auto g_xxxx_0_xxxy_1 = pbuffer.data(idx_eri_1_gsg + 1);

    auto g_xxxx_0_xxxz_1 = pbuffer.data(idx_eri_1_gsg + 2);

    auto g_xxxx_0_xxyy_1 = pbuffer.data(idx_eri_1_gsg + 3);

    auto g_xxxx_0_xxyz_1 = pbuffer.data(idx_eri_1_gsg + 4);

    auto g_xxxx_0_xxzz_1 = pbuffer.data(idx_eri_1_gsg + 5);

    auto g_xxxx_0_xyyy_1 = pbuffer.data(idx_eri_1_gsg + 6);

    auto g_xxxx_0_xyyz_1 = pbuffer.data(idx_eri_1_gsg + 7);

    auto g_xxxx_0_xyzz_1 = pbuffer.data(idx_eri_1_gsg + 8);

    auto g_xxxx_0_xzzz_1 = pbuffer.data(idx_eri_1_gsg + 9);

    auto g_xxxx_0_yyyy_1 = pbuffer.data(idx_eri_1_gsg + 10);

    auto g_xxxx_0_yyyz_1 = pbuffer.data(idx_eri_1_gsg + 11);

    auto g_xxxx_0_yyzz_1 = pbuffer.data(idx_eri_1_gsg + 12);

    auto g_xxxx_0_yzzz_1 = pbuffer.data(idx_eri_1_gsg + 13);

    auto g_xxxx_0_zzzz_1 = pbuffer.data(idx_eri_1_gsg + 14);

    auto g_xxxy_0_xxxx_1 = pbuffer.data(idx_eri_1_gsg + 15);

    auto g_xxxy_0_xxxz_1 = pbuffer.data(idx_eri_1_gsg + 17);

    auto g_xxxy_0_xxzz_1 = pbuffer.data(idx_eri_1_gsg + 20);

    auto g_xxxy_0_xzzz_1 = pbuffer.data(idx_eri_1_gsg + 24);

    auto g_xxxz_0_xxxx_1 = pbuffer.data(idx_eri_1_gsg + 30);

    auto g_xxxz_0_xxxy_1 = pbuffer.data(idx_eri_1_gsg + 31);

    auto g_xxxz_0_xxyy_1 = pbuffer.data(idx_eri_1_gsg + 33);

    auto g_xxxz_0_xyyy_1 = pbuffer.data(idx_eri_1_gsg + 36);

    auto g_xxyy_0_xxxx_1 = pbuffer.data(idx_eri_1_gsg + 45);

    auto g_xxyy_0_xxxy_1 = pbuffer.data(idx_eri_1_gsg + 46);

    auto g_xxyy_0_xxxz_1 = pbuffer.data(idx_eri_1_gsg + 47);

    auto g_xxyy_0_xxyy_1 = pbuffer.data(idx_eri_1_gsg + 48);

    auto g_xxyy_0_xxyz_1 = pbuffer.data(idx_eri_1_gsg + 49);

    auto g_xxyy_0_xxzz_1 = pbuffer.data(idx_eri_1_gsg + 50);

    auto g_xxyy_0_xyyy_1 = pbuffer.data(idx_eri_1_gsg + 51);

    auto g_xxyy_0_xyyz_1 = pbuffer.data(idx_eri_1_gsg + 52);

    auto g_xxyy_0_xyzz_1 = pbuffer.data(idx_eri_1_gsg + 53);

    auto g_xxyy_0_xzzz_1 = pbuffer.data(idx_eri_1_gsg + 54);

    auto g_xxyy_0_yyyy_1 = pbuffer.data(idx_eri_1_gsg + 55);

    auto g_xxyy_0_yyyz_1 = pbuffer.data(idx_eri_1_gsg + 56);

    auto g_xxyy_0_yyzz_1 = pbuffer.data(idx_eri_1_gsg + 57);

    auto g_xxyy_0_yzzz_1 = pbuffer.data(idx_eri_1_gsg + 58);

    auto g_xxyy_0_zzzz_1 = pbuffer.data(idx_eri_1_gsg + 59);

    auto g_xxzz_0_xxxx_1 = pbuffer.data(idx_eri_1_gsg + 75);

    auto g_xxzz_0_xxxy_1 = pbuffer.data(idx_eri_1_gsg + 76);

    auto g_xxzz_0_xxxz_1 = pbuffer.data(idx_eri_1_gsg + 77);

    auto g_xxzz_0_xxyy_1 = pbuffer.data(idx_eri_1_gsg + 78);

    auto g_xxzz_0_xxyz_1 = pbuffer.data(idx_eri_1_gsg + 79);

    auto g_xxzz_0_xxzz_1 = pbuffer.data(idx_eri_1_gsg + 80);

    auto g_xxzz_0_xyyy_1 = pbuffer.data(idx_eri_1_gsg + 81);

    auto g_xxzz_0_xyyz_1 = pbuffer.data(idx_eri_1_gsg + 82);

    auto g_xxzz_0_xyzz_1 = pbuffer.data(idx_eri_1_gsg + 83);

    auto g_xxzz_0_xzzz_1 = pbuffer.data(idx_eri_1_gsg + 84);

    auto g_xxzz_0_yyyy_1 = pbuffer.data(idx_eri_1_gsg + 85);

    auto g_xxzz_0_yyyz_1 = pbuffer.data(idx_eri_1_gsg + 86);

    auto g_xxzz_0_yyzz_1 = pbuffer.data(idx_eri_1_gsg + 87);

    auto g_xxzz_0_yzzz_1 = pbuffer.data(idx_eri_1_gsg + 88);

    auto g_xxzz_0_zzzz_1 = pbuffer.data(idx_eri_1_gsg + 89);

    auto g_xyyy_0_xxxy_1 = pbuffer.data(idx_eri_1_gsg + 91);

    auto g_xyyy_0_xxyy_1 = pbuffer.data(idx_eri_1_gsg + 93);

    auto g_xyyy_0_xxyz_1 = pbuffer.data(idx_eri_1_gsg + 94);

    auto g_xyyy_0_xyyy_1 = pbuffer.data(idx_eri_1_gsg + 96);

    auto g_xyyy_0_xyyz_1 = pbuffer.data(idx_eri_1_gsg + 97);

    auto g_xyyy_0_xyzz_1 = pbuffer.data(idx_eri_1_gsg + 98);

    auto g_xyyy_0_yyyy_1 = pbuffer.data(idx_eri_1_gsg + 100);

    auto g_xyyy_0_yyyz_1 = pbuffer.data(idx_eri_1_gsg + 101);

    auto g_xyyy_0_yyzz_1 = pbuffer.data(idx_eri_1_gsg + 102);

    auto g_xyyy_0_yzzz_1 = pbuffer.data(idx_eri_1_gsg + 103);

    auto g_xyyy_0_zzzz_1 = pbuffer.data(idx_eri_1_gsg + 104);

    auto g_xzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_gsg + 137);

    auto g_xzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_gsg + 139);

    auto g_xzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_gsg + 140);

    auto g_xzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_gsg + 142);

    auto g_xzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_gsg + 143);

    auto g_xzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_gsg + 144);

    auto g_xzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_gsg + 145);

    auto g_xzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_gsg + 146);

    auto g_xzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_gsg + 147);

    auto g_xzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_gsg + 148);

    auto g_xzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_gsg + 149);

    auto g_yyyy_0_xxxx_1 = pbuffer.data(idx_eri_1_gsg + 150);

    auto g_yyyy_0_xxxy_1 = pbuffer.data(idx_eri_1_gsg + 151);

    auto g_yyyy_0_xxxz_1 = pbuffer.data(idx_eri_1_gsg + 152);

    auto g_yyyy_0_xxyy_1 = pbuffer.data(idx_eri_1_gsg + 153);

    auto g_yyyy_0_xxyz_1 = pbuffer.data(idx_eri_1_gsg + 154);

    auto g_yyyy_0_xxzz_1 = pbuffer.data(idx_eri_1_gsg + 155);

    auto g_yyyy_0_xyyy_1 = pbuffer.data(idx_eri_1_gsg + 156);

    auto g_yyyy_0_xyyz_1 = pbuffer.data(idx_eri_1_gsg + 157);

    auto g_yyyy_0_xyzz_1 = pbuffer.data(idx_eri_1_gsg + 158);

    auto g_yyyy_0_xzzz_1 = pbuffer.data(idx_eri_1_gsg + 159);

    auto g_yyyy_0_yyyy_1 = pbuffer.data(idx_eri_1_gsg + 160);

    auto g_yyyy_0_yyyz_1 = pbuffer.data(idx_eri_1_gsg + 161);

    auto g_yyyy_0_yyzz_1 = pbuffer.data(idx_eri_1_gsg + 162);

    auto g_yyyy_0_yzzz_1 = pbuffer.data(idx_eri_1_gsg + 163);

    auto g_yyyy_0_zzzz_1 = pbuffer.data(idx_eri_1_gsg + 164);

    auto g_yyyz_0_xxxy_1 = pbuffer.data(idx_eri_1_gsg + 166);

    auto g_yyyz_0_xxyy_1 = pbuffer.data(idx_eri_1_gsg + 168);

    auto g_yyyz_0_xyyy_1 = pbuffer.data(idx_eri_1_gsg + 171);

    auto g_yyyz_0_yyyy_1 = pbuffer.data(idx_eri_1_gsg + 175);

    auto g_yyzz_0_xxxx_1 = pbuffer.data(idx_eri_1_gsg + 180);

    auto g_yyzz_0_xxxy_1 = pbuffer.data(idx_eri_1_gsg + 181);

    auto g_yyzz_0_xxxz_1 = pbuffer.data(idx_eri_1_gsg + 182);

    auto g_yyzz_0_xxyy_1 = pbuffer.data(idx_eri_1_gsg + 183);

    auto g_yyzz_0_xxyz_1 = pbuffer.data(idx_eri_1_gsg + 184);

    auto g_yyzz_0_xxzz_1 = pbuffer.data(idx_eri_1_gsg + 185);

    auto g_yyzz_0_xyyy_1 = pbuffer.data(idx_eri_1_gsg + 186);

    auto g_yyzz_0_xyyz_1 = pbuffer.data(idx_eri_1_gsg + 187);

    auto g_yyzz_0_xyzz_1 = pbuffer.data(idx_eri_1_gsg + 188);

    auto g_yyzz_0_xzzz_1 = pbuffer.data(idx_eri_1_gsg + 189);

    auto g_yyzz_0_yyyy_1 = pbuffer.data(idx_eri_1_gsg + 190);

    auto g_yyzz_0_yyyz_1 = pbuffer.data(idx_eri_1_gsg + 191);

    auto g_yyzz_0_yyzz_1 = pbuffer.data(idx_eri_1_gsg + 192);

    auto g_yyzz_0_yzzz_1 = pbuffer.data(idx_eri_1_gsg + 193);

    auto g_yyzz_0_zzzz_1 = pbuffer.data(idx_eri_1_gsg + 194);

    auto g_yzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_gsg + 195);

    auto g_yzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_gsg + 197);

    auto g_yzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_gsg + 199);

    auto g_yzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_gsg + 200);

    auto g_yzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_gsg + 202);

    auto g_yzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_gsg + 203);

    auto g_yzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_gsg + 204);

    auto g_yzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_gsg + 206);

    auto g_yzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_gsg + 207);

    auto g_yzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_gsg + 208);

    auto g_yzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_gsg + 209);

    auto g_zzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_gsg + 210);

    auto g_zzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_gsg + 211);

    auto g_zzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_gsg + 212);

    auto g_zzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_gsg + 213);

    auto g_zzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_gsg + 214);

    auto g_zzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_gsg + 215);

    auto g_zzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_gsg + 216);

    auto g_zzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_gsg + 217);

    auto g_zzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_gsg + 218);

    auto g_zzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_gsg + 219);

    auto g_zzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_gsg + 220);

    auto g_zzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_gsg + 221);

    auto g_zzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_gsg + 222);

    auto g_zzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_gsg + 223);

    auto g_zzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_gsg + 224);

    /// Set up components of auxilary buffer : HSF

    auto g_xxxxx_0_xxx_1 = pbuffer.data(idx_eri_1_hsf);

    auto g_xxxxx_0_xxy_1 = pbuffer.data(idx_eri_1_hsf + 1);

    auto g_xxxxx_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 2);

    auto g_xxxxx_0_xyy_1 = pbuffer.data(idx_eri_1_hsf + 3);

    auto g_xxxxx_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 4);

    auto g_xxxxx_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 5);

    auto g_xxxxx_0_yyy_1 = pbuffer.data(idx_eri_1_hsf + 6);

    auto g_xxxxx_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 7);

    auto g_xxxxx_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 8);

    auto g_xxxxx_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 9);

    auto g_xxxxz_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 22);

    auto g_xxxxz_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 24);

    auto g_xxxxz_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 25);

    auto g_xxxxz_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 27);

    auto g_xxxxz_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 28);

    auto g_xxxxz_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 29);

    auto g_xxxyy_0_xxx_1 = pbuffer.data(idx_eri_1_hsf + 30);

    auto g_xxxyy_0_xxy_1 = pbuffer.data(idx_eri_1_hsf + 31);

    auto g_xxxyy_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 32);

    auto g_xxxyy_0_xyy_1 = pbuffer.data(idx_eri_1_hsf + 33);

    auto g_xxxyy_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 34);

    auto g_xxxyy_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 35);

    auto g_xxxyy_0_yyy_1 = pbuffer.data(idx_eri_1_hsf + 36);

    auto g_xxxyy_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 37);

    auto g_xxxyy_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 38);

    auto g_xxxyy_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 39);

    auto g_xxxzz_0_xxx_1 = pbuffer.data(idx_eri_1_hsf + 50);

    auto g_xxxzz_0_xxy_1 = pbuffer.data(idx_eri_1_hsf + 51);

    auto g_xxxzz_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 52);

    auto g_xxxzz_0_xyy_1 = pbuffer.data(idx_eri_1_hsf + 53);

    auto g_xxxzz_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 54);

    auto g_xxxzz_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 55);

    auto g_xxxzz_0_yyy_1 = pbuffer.data(idx_eri_1_hsf + 56);

    auto g_xxxzz_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 57);

    auto g_xxxzz_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 58);

    auto g_xxxzz_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 59);

    auto g_xxyyy_0_xxx_1 = pbuffer.data(idx_eri_1_hsf + 60);

    auto g_xxyyy_0_xxy_1 = pbuffer.data(idx_eri_1_hsf + 61);

    auto g_xxyyy_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 62);

    auto g_xxyyy_0_xyy_1 = pbuffer.data(idx_eri_1_hsf + 63);

    auto g_xxyyy_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 64);

    auto g_xxyyy_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 65);

    auto g_xxyyy_0_yyy_1 = pbuffer.data(idx_eri_1_hsf + 66);

    auto g_xxyyy_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 67);

    auto g_xxyyy_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 68);

    auto g_xxyyy_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 69);

    auto g_xxzzz_0_xxx_1 = pbuffer.data(idx_eri_1_hsf + 90);

    auto g_xxzzz_0_xxy_1 = pbuffer.data(idx_eri_1_hsf + 91);

    auto g_xxzzz_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 92);

    auto g_xxzzz_0_xyy_1 = pbuffer.data(idx_eri_1_hsf + 93);

    auto g_xxzzz_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 94);

    auto g_xxzzz_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 95);

    auto g_xxzzz_0_yyy_1 = pbuffer.data(idx_eri_1_hsf + 96);

    auto g_xxzzz_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 97);

    auto g_xxzzz_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 98);

    auto g_xxzzz_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 99);

    auto g_xyyyy_0_xxy_1 = pbuffer.data(idx_eri_1_hsf + 101);

    auto g_xyyyy_0_xyy_1 = pbuffer.data(idx_eri_1_hsf + 103);

    auto g_xyyyy_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 104);

    auto g_xyyyy_0_yyy_1 = pbuffer.data(idx_eri_1_hsf + 106);

    auto g_xyyyy_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 107);

    auto g_xyyyy_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 108);

    auto g_xyyzz_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 124);

    auto g_xyyzz_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 127);

    auto g_xyyzz_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 128);

    auto g_xzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 142);

    auto g_xzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 144);

    auto g_xzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 145);

    auto g_xzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 147);

    auto g_xzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 148);

    auto g_xzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 149);

    auto g_yyyyy_0_xxx_1 = pbuffer.data(idx_eri_1_hsf + 150);

    auto g_yyyyy_0_xxy_1 = pbuffer.data(idx_eri_1_hsf + 151);

    auto g_yyyyy_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 152);

    auto g_yyyyy_0_xyy_1 = pbuffer.data(idx_eri_1_hsf + 153);

    auto g_yyyyy_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 154);

    auto g_yyyyy_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 155);

    auto g_yyyyy_0_yyy_1 = pbuffer.data(idx_eri_1_hsf + 156);

    auto g_yyyyy_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 157);

    auto g_yyyyy_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 158);

    auto g_yyyyy_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 159);

    auto g_yyyyz_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 162);

    auto g_yyyyz_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 164);

    auto g_yyyyz_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 165);

    auto g_yyyyz_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 167);

    auto g_yyyyz_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 168);

    auto g_yyyyz_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 169);

    auto g_yyyzz_0_xxx_1 = pbuffer.data(idx_eri_1_hsf + 170);

    auto g_yyyzz_0_xxy_1 = pbuffer.data(idx_eri_1_hsf + 171);

    auto g_yyyzz_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 172);

    auto g_yyyzz_0_xyy_1 = pbuffer.data(idx_eri_1_hsf + 173);

    auto g_yyyzz_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 174);

    auto g_yyyzz_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 175);

    auto g_yyyzz_0_yyy_1 = pbuffer.data(idx_eri_1_hsf + 176);

    auto g_yyyzz_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 177);

    auto g_yyyzz_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 178);

    auto g_yyyzz_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 179);

    auto g_yyzzz_0_xxx_1 = pbuffer.data(idx_eri_1_hsf + 180);

    auto g_yyzzz_0_xxy_1 = pbuffer.data(idx_eri_1_hsf + 181);

    auto g_yyzzz_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 182);

    auto g_yyzzz_0_xyy_1 = pbuffer.data(idx_eri_1_hsf + 183);

    auto g_yyzzz_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 184);

    auto g_yyzzz_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 185);

    auto g_yyzzz_0_yyy_1 = pbuffer.data(idx_eri_1_hsf + 186);

    auto g_yyzzz_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 187);

    auto g_yyzzz_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 188);

    auto g_yyzzz_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 189);

    auto g_yzzzz_0_xxy_1 = pbuffer.data(idx_eri_1_hsf + 191);

    auto g_yzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 192);

    auto g_yzzzz_0_xyy_1 = pbuffer.data(idx_eri_1_hsf + 193);

    auto g_yzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 194);

    auto g_yzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 195);

    auto g_yzzzz_0_yyy_1 = pbuffer.data(idx_eri_1_hsf + 196);

    auto g_yzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 197);

    auto g_yzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 198);

    auto g_yzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 199);

    auto g_zzzzz_0_xxx_1 = pbuffer.data(idx_eri_1_hsf + 200);

    auto g_zzzzz_0_xxy_1 = pbuffer.data(idx_eri_1_hsf + 201);

    auto g_zzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 202);

    auto g_zzzzz_0_xyy_1 = pbuffer.data(idx_eri_1_hsf + 203);

    auto g_zzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 204);

    auto g_zzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 205);

    auto g_zzzzz_0_yyy_1 = pbuffer.data(idx_eri_1_hsf + 206);

    auto g_zzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 207);

    auto g_zzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 208);

    auto g_zzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 209);

    /// Set up components of auxilary buffer : HSG

    auto g_xxxxx_0_xxxx_1 = pbuffer.data(idx_eri_1_hsg);

    auto g_xxxxx_0_xxxy_1 = pbuffer.data(idx_eri_1_hsg + 1);

    auto g_xxxxx_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 2);

    auto g_xxxxx_0_xxyy_1 = pbuffer.data(idx_eri_1_hsg + 3);

    auto g_xxxxx_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 4);

    auto g_xxxxx_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 5);

    auto g_xxxxx_0_xyyy_1 = pbuffer.data(idx_eri_1_hsg + 6);

    auto g_xxxxx_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 7);

    auto g_xxxxx_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 8);

    auto g_xxxxx_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 9);

    auto g_xxxxx_0_yyyy_1 = pbuffer.data(idx_eri_1_hsg + 10);

    auto g_xxxxx_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 11);

    auto g_xxxxx_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 12);

    auto g_xxxxx_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 13);

    auto g_xxxxx_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 14);

    auto g_xxxxy_0_xxxx_1 = pbuffer.data(idx_eri_1_hsg + 15);

    auto g_xxxxy_0_xxxy_1 = pbuffer.data(idx_eri_1_hsg + 16);

    auto g_xxxxy_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 17);

    auto g_xxxxy_0_xxyy_1 = pbuffer.data(idx_eri_1_hsg + 18);

    auto g_xxxxy_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 20);

    auto g_xxxxy_0_xyyy_1 = pbuffer.data(idx_eri_1_hsg + 21);

    auto g_xxxxy_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 24);

    auto g_xxxxy_0_yyyy_1 = pbuffer.data(idx_eri_1_hsg + 25);

    auto g_xxxxz_0_xxxx_1 = pbuffer.data(idx_eri_1_hsg + 30);

    auto g_xxxxz_0_xxxy_1 = pbuffer.data(idx_eri_1_hsg + 31);

    auto g_xxxxz_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 32);

    auto g_xxxxz_0_xxyy_1 = pbuffer.data(idx_eri_1_hsg + 33);

    auto g_xxxxz_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 34);

    auto g_xxxxz_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 35);

    auto g_xxxxz_0_xyyy_1 = pbuffer.data(idx_eri_1_hsg + 36);

    auto g_xxxxz_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 37);

    auto g_xxxxz_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 38);

    auto g_xxxxz_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 39);

    auto g_xxxxz_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 41);

    auto g_xxxxz_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 42);

    auto g_xxxxz_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 43);

    auto g_xxxxz_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 44);

    auto g_xxxyy_0_xxxx_1 = pbuffer.data(idx_eri_1_hsg + 45);

    auto g_xxxyy_0_xxxy_1 = pbuffer.data(idx_eri_1_hsg + 46);

    auto g_xxxyy_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 47);

    auto g_xxxyy_0_xxyy_1 = pbuffer.data(idx_eri_1_hsg + 48);

    auto g_xxxyy_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 49);

    auto g_xxxyy_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 50);

    auto g_xxxyy_0_xyyy_1 = pbuffer.data(idx_eri_1_hsg + 51);

    auto g_xxxyy_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 52);

    auto g_xxxyy_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 53);

    auto g_xxxyy_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 54);

    auto g_xxxyy_0_yyyy_1 = pbuffer.data(idx_eri_1_hsg + 55);

    auto g_xxxyy_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 56);

    auto g_xxxyy_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 57);

    auto g_xxxyy_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 58);

    auto g_xxxyy_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 59);

    auto g_xxxzz_0_xxxx_1 = pbuffer.data(idx_eri_1_hsg + 75);

    auto g_xxxzz_0_xxxy_1 = pbuffer.data(idx_eri_1_hsg + 76);

    auto g_xxxzz_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 77);

    auto g_xxxzz_0_xxyy_1 = pbuffer.data(idx_eri_1_hsg + 78);

    auto g_xxxzz_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 79);

    auto g_xxxzz_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 80);

    auto g_xxxzz_0_xyyy_1 = pbuffer.data(idx_eri_1_hsg + 81);

    auto g_xxxzz_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 82);

    auto g_xxxzz_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 83);

    auto g_xxxzz_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 84);

    auto g_xxxzz_0_yyyy_1 = pbuffer.data(idx_eri_1_hsg + 85);

    auto g_xxxzz_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 86);

    auto g_xxxzz_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 87);

    auto g_xxxzz_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 88);

    auto g_xxxzz_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 89);

    auto g_xxyyy_0_xxxx_1 = pbuffer.data(idx_eri_1_hsg + 90);

    auto g_xxyyy_0_xxxy_1 = pbuffer.data(idx_eri_1_hsg + 91);

    auto g_xxyyy_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 92);

    auto g_xxyyy_0_xxyy_1 = pbuffer.data(idx_eri_1_hsg + 93);

    auto g_xxyyy_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 94);

    auto g_xxyyy_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 95);

    auto g_xxyyy_0_xyyy_1 = pbuffer.data(idx_eri_1_hsg + 96);

    auto g_xxyyy_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 97);

    auto g_xxyyy_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 98);

    auto g_xxyyy_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 99);

    auto g_xxyyy_0_yyyy_1 = pbuffer.data(idx_eri_1_hsg + 100);

    auto g_xxyyy_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 101);

    auto g_xxyyy_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 102);

    auto g_xxyyy_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 103);

    auto g_xxyyy_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 104);

    auto g_xxyyz_0_xxxy_1 = pbuffer.data(idx_eri_1_hsg + 106);

    auto g_xxyyz_0_xxyy_1 = pbuffer.data(idx_eri_1_hsg + 108);

    auto g_xxyyz_0_xyyy_1 = pbuffer.data(idx_eri_1_hsg + 111);

    auto g_xxyzz_0_xxxx_1 = pbuffer.data(idx_eri_1_hsg + 120);

    auto g_xxyzz_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 122);

    auto g_xxyzz_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 125);

    auto g_xxyzz_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 129);

    auto g_xxzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_hsg + 135);

    auto g_xxzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_hsg + 136);

    auto g_xxzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 137);

    auto g_xxzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_hsg + 138);

    auto g_xxzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 139);

    auto g_xxzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 140);

    auto g_xxzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_hsg + 141);

    auto g_xxzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 142);

    auto g_xxzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 143);

    auto g_xxzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 144);

    auto g_xxzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_hsg + 145);

    auto g_xxzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 146);

    auto g_xxzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 147);

    auto g_xxzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 148);

    auto g_xxzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 149);

    auto g_xyyyy_0_xxxx_1 = pbuffer.data(idx_eri_1_hsg + 150);

    auto g_xyyyy_0_xxxy_1 = pbuffer.data(idx_eri_1_hsg + 151);

    auto g_xyyyy_0_xxyy_1 = pbuffer.data(idx_eri_1_hsg + 153);

    auto g_xyyyy_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 154);

    auto g_xyyyy_0_xyyy_1 = pbuffer.data(idx_eri_1_hsg + 156);

    auto g_xyyyy_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 157);

    auto g_xyyyy_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 158);

    auto g_xyyyy_0_yyyy_1 = pbuffer.data(idx_eri_1_hsg + 160);

    auto g_xyyyy_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 161);

    auto g_xyyyy_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 162);

    auto g_xyyyy_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 163);

    auto g_xyyyy_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 164);

    auto g_xyyzz_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 184);

    auto g_xyyzz_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 187);

    auto g_xyyzz_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 188);

    auto g_xyyzz_0_yyyy_1 = pbuffer.data(idx_eri_1_hsg + 190);

    auto g_xyyzz_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 191);

    auto g_xyyzz_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 192);

    auto g_xyyzz_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 193);

    auto g_xyyzz_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 194);

    auto g_xzzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_hsg + 210);

    auto g_xzzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 212);

    auto g_xzzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 214);

    auto g_xzzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 215);

    auto g_xzzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 217);

    auto g_xzzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 218);

    auto g_xzzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 219);

    auto g_xzzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_hsg + 220);

    auto g_xzzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 221);

    auto g_xzzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 222);

    auto g_xzzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 223);

    auto g_xzzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 224);

    auto g_yyyyy_0_xxxx_1 = pbuffer.data(idx_eri_1_hsg + 225);

    auto g_yyyyy_0_xxxy_1 = pbuffer.data(idx_eri_1_hsg + 226);

    auto g_yyyyy_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 227);

    auto g_yyyyy_0_xxyy_1 = pbuffer.data(idx_eri_1_hsg + 228);

    auto g_yyyyy_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 229);

    auto g_yyyyy_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 230);

    auto g_yyyyy_0_xyyy_1 = pbuffer.data(idx_eri_1_hsg + 231);

    auto g_yyyyy_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 232);

    auto g_yyyyy_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 233);

    auto g_yyyyy_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 234);

    auto g_yyyyy_0_yyyy_1 = pbuffer.data(idx_eri_1_hsg + 235);

    auto g_yyyyy_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 236);

    auto g_yyyyy_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 237);

    auto g_yyyyy_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 238);

    auto g_yyyyy_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 239);

    auto g_yyyyz_0_xxxy_1 = pbuffer.data(idx_eri_1_hsg + 241);

    auto g_yyyyz_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 242);

    auto g_yyyyz_0_xxyy_1 = pbuffer.data(idx_eri_1_hsg + 243);

    auto g_yyyyz_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 244);

    auto g_yyyyz_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 245);

    auto g_yyyyz_0_xyyy_1 = pbuffer.data(idx_eri_1_hsg + 246);

    auto g_yyyyz_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 247);

    auto g_yyyyz_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 248);

    auto g_yyyyz_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 249);

    auto g_yyyyz_0_yyyy_1 = pbuffer.data(idx_eri_1_hsg + 250);

    auto g_yyyyz_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 251);

    auto g_yyyyz_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 252);

    auto g_yyyyz_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 253);

    auto g_yyyyz_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 254);

    auto g_yyyzz_0_xxxx_1 = pbuffer.data(idx_eri_1_hsg + 255);

    auto g_yyyzz_0_xxxy_1 = pbuffer.data(idx_eri_1_hsg + 256);

    auto g_yyyzz_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 257);

    auto g_yyyzz_0_xxyy_1 = pbuffer.data(idx_eri_1_hsg + 258);

    auto g_yyyzz_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 259);

    auto g_yyyzz_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 260);

    auto g_yyyzz_0_xyyy_1 = pbuffer.data(idx_eri_1_hsg + 261);

    auto g_yyyzz_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 262);

    auto g_yyyzz_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 263);

    auto g_yyyzz_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 264);

    auto g_yyyzz_0_yyyy_1 = pbuffer.data(idx_eri_1_hsg + 265);

    auto g_yyyzz_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 266);

    auto g_yyyzz_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 267);

    auto g_yyyzz_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 268);

    auto g_yyyzz_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 269);

    auto g_yyzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_hsg + 270);

    auto g_yyzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_hsg + 271);

    auto g_yyzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 272);

    auto g_yyzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_hsg + 273);

    auto g_yyzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 274);

    auto g_yyzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 275);

    auto g_yyzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_hsg + 276);

    auto g_yyzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 277);

    auto g_yyzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 278);

    auto g_yyzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 279);

    auto g_yyzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_hsg + 280);

    auto g_yyzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 281);

    auto g_yyzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 282);

    auto g_yyzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 283);

    auto g_yyzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 284);

    auto g_yzzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_hsg + 285);

    auto g_yzzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_hsg + 286);

    auto g_yzzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 287);

    auto g_yzzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_hsg + 288);

    auto g_yzzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 289);

    auto g_yzzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 290);

    auto g_yzzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_hsg + 291);

    auto g_yzzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 292);

    auto g_yzzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 293);

    auto g_yzzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 294);

    auto g_yzzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_hsg + 295);

    auto g_yzzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 296);

    auto g_yzzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 297);

    auto g_yzzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 298);

    auto g_yzzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 299);

    auto g_zzzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_hsg + 300);

    auto g_zzzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_hsg + 301);

    auto g_zzzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 302);

    auto g_zzzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_hsg + 303);

    auto g_zzzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 304);

    auto g_zzzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 305);

    auto g_zzzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_hsg + 306);

    auto g_zzzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 307);

    auto g_zzzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 308);

    auto g_zzzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 309);

    auto g_zzzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_hsg + 310);

    auto g_zzzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 311);

    auto g_zzzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 312);

    auto g_zzzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 313);

    auto g_zzzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 314);

    /// Set up 0-15 components of targeted buffer : ISG

    auto g_xxxxxx_0_xxxx_0 = pbuffer.data(idx_eri_0_isg);

    auto g_xxxxxx_0_xxxy_0 = pbuffer.data(idx_eri_0_isg + 1);

    auto g_xxxxxx_0_xxxz_0 = pbuffer.data(idx_eri_0_isg + 2);

    auto g_xxxxxx_0_xxyy_0 = pbuffer.data(idx_eri_0_isg + 3);

    auto g_xxxxxx_0_xxyz_0 = pbuffer.data(idx_eri_0_isg + 4);

    auto g_xxxxxx_0_xxzz_0 = pbuffer.data(idx_eri_0_isg + 5);

    auto g_xxxxxx_0_xyyy_0 = pbuffer.data(idx_eri_0_isg + 6);

    auto g_xxxxxx_0_xyyz_0 = pbuffer.data(idx_eri_0_isg + 7);

    auto g_xxxxxx_0_xyzz_0 = pbuffer.data(idx_eri_0_isg + 8);

    auto g_xxxxxx_0_xzzz_0 = pbuffer.data(idx_eri_0_isg + 9);

    auto g_xxxxxx_0_yyyy_0 = pbuffer.data(idx_eri_0_isg + 10);

    auto g_xxxxxx_0_yyyz_0 = pbuffer.data(idx_eri_0_isg + 11);

    auto g_xxxxxx_0_yyzz_0 = pbuffer.data(idx_eri_0_isg + 12);

    auto g_xxxxxx_0_yzzz_0 = pbuffer.data(idx_eri_0_isg + 13);

    auto g_xxxxxx_0_zzzz_0 = pbuffer.data(idx_eri_0_isg + 14);

    #pragma omp simd aligned(g_xxxx_0_xxxx_0, g_xxxx_0_xxxx_1, g_xxxx_0_xxxy_0, g_xxxx_0_xxxy_1, g_xxxx_0_xxxz_0, g_xxxx_0_xxxz_1, g_xxxx_0_xxyy_0, g_xxxx_0_xxyy_1, g_xxxx_0_xxyz_0, g_xxxx_0_xxyz_1, g_xxxx_0_xxzz_0, g_xxxx_0_xxzz_1, g_xxxx_0_xyyy_0, g_xxxx_0_xyyy_1, g_xxxx_0_xyyz_0, g_xxxx_0_xyyz_1, g_xxxx_0_xyzz_0, g_xxxx_0_xyzz_1, g_xxxx_0_xzzz_0, g_xxxx_0_xzzz_1, g_xxxx_0_yyyy_0, g_xxxx_0_yyyy_1, g_xxxx_0_yyyz_0, g_xxxx_0_yyyz_1, g_xxxx_0_yyzz_0, g_xxxx_0_yyzz_1, g_xxxx_0_yzzz_0, g_xxxx_0_yzzz_1, g_xxxx_0_zzzz_0, g_xxxx_0_zzzz_1, g_xxxxx_0_xxx_1, g_xxxxx_0_xxxx_1, g_xxxxx_0_xxxy_1, g_xxxxx_0_xxxz_1, g_xxxxx_0_xxy_1, g_xxxxx_0_xxyy_1, g_xxxxx_0_xxyz_1, g_xxxxx_0_xxz_1, g_xxxxx_0_xxzz_1, g_xxxxx_0_xyy_1, g_xxxxx_0_xyyy_1, g_xxxxx_0_xyyz_1, g_xxxxx_0_xyz_1, g_xxxxx_0_xyzz_1, g_xxxxx_0_xzz_1, g_xxxxx_0_xzzz_1, g_xxxxx_0_yyy_1, g_xxxxx_0_yyyy_1, g_xxxxx_0_yyyz_1, g_xxxxx_0_yyz_1, g_xxxxx_0_yyzz_1, g_xxxxx_0_yzz_1, g_xxxxx_0_yzzz_1, g_xxxxx_0_zzz_1, g_xxxxx_0_zzzz_1, g_xxxxxx_0_xxxx_0, g_xxxxxx_0_xxxy_0, g_xxxxxx_0_xxxz_0, g_xxxxxx_0_xxyy_0, g_xxxxxx_0_xxyz_0, g_xxxxxx_0_xxzz_0, g_xxxxxx_0_xyyy_0, g_xxxxxx_0_xyyz_0, g_xxxxxx_0_xyzz_0, g_xxxxxx_0_xzzz_0, g_xxxxxx_0_yyyy_0, g_xxxxxx_0_yyyz_0, g_xxxxxx_0_yyzz_0, g_xxxxxx_0_yzzz_0, g_xxxxxx_0_zzzz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxxx_0_xxxx_0[i] = 5.0 * g_xxxx_0_xxxx_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxx_1[i] * fz_be_0 + 4.0 * g_xxxxx_0_xxx_1[i] * fi_acd_0 + g_xxxxx_0_xxxx_1[i] * wa_x[i];

        g_xxxxxx_0_xxxy_0[i] = 5.0 * g_xxxx_0_xxxy_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxy_1[i] * fz_be_0 + 3.0 * g_xxxxx_0_xxy_1[i] * fi_acd_0 + g_xxxxx_0_xxxy_1[i] * wa_x[i];

        g_xxxxxx_0_xxxz_0[i] = 5.0 * g_xxxx_0_xxxz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxz_1[i] * fz_be_0 + 3.0 * g_xxxxx_0_xxz_1[i] * fi_acd_0 + g_xxxxx_0_xxxz_1[i] * wa_x[i];

        g_xxxxxx_0_xxyy_0[i] = 5.0 * g_xxxx_0_xxyy_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxyy_1[i] * fz_be_0 + 2.0 * g_xxxxx_0_xyy_1[i] * fi_acd_0 + g_xxxxx_0_xxyy_1[i] * wa_x[i];

        g_xxxxxx_0_xxyz_0[i] = 5.0 * g_xxxx_0_xxyz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxyz_1[i] * fz_be_0 + 2.0 * g_xxxxx_0_xyz_1[i] * fi_acd_0 + g_xxxxx_0_xxyz_1[i] * wa_x[i];

        g_xxxxxx_0_xxzz_0[i] = 5.0 * g_xxxx_0_xxzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxzz_1[i] * fz_be_0 + 2.0 * g_xxxxx_0_xzz_1[i] * fi_acd_0 + g_xxxxx_0_xxzz_1[i] * wa_x[i];

        g_xxxxxx_0_xyyy_0[i] = 5.0 * g_xxxx_0_xyyy_0[i] * fbe_0 - 5.0 * g_xxxx_0_xyyy_1[i] * fz_be_0 + g_xxxxx_0_yyy_1[i] * fi_acd_0 + g_xxxxx_0_xyyy_1[i] * wa_x[i];

        g_xxxxxx_0_xyyz_0[i] = 5.0 * g_xxxx_0_xyyz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xyyz_1[i] * fz_be_0 + g_xxxxx_0_yyz_1[i] * fi_acd_0 + g_xxxxx_0_xyyz_1[i] * wa_x[i];

        g_xxxxxx_0_xyzz_0[i] = 5.0 * g_xxxx_0_xyzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xyzz_1[i] * fz_be_0 + g_xxxxx_0_yzz_1[i] * fi_acd_0 + g_xxxxx_0_xyzz_1[i] * wa_x[i];

        g_xxxxxx_0_xzzz_0[i] = 5.0 * g_xxxx_0_xzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xzzz_1[i] * fz_be_0 + g_xxxxx_0_zzz_1[i] * fi_acd_0 + g_xxxxx_0_xzzz_1[i] * wa_x[i];

        g_xxxxxx_0_yyyy_0[i] = 5.0 * g_xxxx_0_yyyy_0[i] * fbe_0 - 5.0 * g_xxxx_0_yyyy_1[i] * fz_be_0 + g_xxxxx_0_yyyy_1[i] * wa_x[i];

        g_xxxxxx_0_yyyz_0[i] = 5.0 * g_xxxx_0_yyyz_0[i] * fbe_0 - 5.0 * g_xxxx_0_yyyz_1[i] * fz_be_0 + g_xxxxx_0_yyyz_1[i] * wa_x[i];

        g_xxxxxx_0_yyzz_0[i] = 5.0 * g_xxxx_0_yyzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_yyzz_1[i] * fz_be_0 + g_xxxxx_0_yyzz_1[i] * wa_x[i];

        g_xxxxxx_0_yzzz_0[i] = 5.0 * g_xxxx_0_yzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_yzzz_1[i] * fz_be_0 + g_xxxxx_0_yzzz_1[i] * wa_x[i];

        g_xxxxxx_0_zzzz_0[i] = 5.0 * g_xxxx_0_zzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_zzzz_1[i] * fz_be_0 + g_xxxxx_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 15-30 components of targeted buffer : ISG

    auto g_xxxxxy_0_xxxx_0 = pbuffer.data(idx_eri_0_isg + 15);

    auto g_xxxxxy_0_xxxy_0 = pbuffer.data(idx_eri_0_isg + 16);

    auto g_xxxxxy_0_xxxz_0 = pbuffer.data(idx_eri_0_isg + 17);

    auto g_xxxxxy_0_xxyy_0 = pbuffer.data(idx_eri_0_isg + 18);

    auto g_xxxxxy_0_xxyz_0 = pbuffer.data(idx_eri_0_isg + 19);

    auto g_xxxxxy_0_xxzz_0 = pbuffer.data(idx_eri_0_isg + 20);

    auto g_xxxxxy_0_xyyy_0 = pbuffer.data(idx_eri_0_isg + 21);

    auto g_xxxxxy_0_xyyz_0 = pbuffer.data(idx_eri_0_isg + 22);

    auto g_xxxxxy_0_xyzz_0 = pbuffer.data(idx_eri_0_isg + 23);

    auto g_xxxxxy_0_xzzz_0 = pbuffer.data(idx_eri_0_isg + 24);

    auto g_xxxxxy_0_yyyy_0 = pbuffer.data(idx_eri_0_isg + 25);

    auto g_xxxxxy_0_yyyz_0 = pbuffer.data(idx_eri_0_isg + 26);

    auto g_xxxxxy_0_yyzz_0 = pbuffer.data(idx_eri_0_isg + 27);

    auto g_xxxxxy_0_yzzz_0 = pbuffer.data(idx_eri_0_isg + 28);

    auto g_xxxxxy_0_zzzz_0 = pbuffer.data(idx_eri_0_isg + 29);

    #pragma omp simd aligned(g_xxxxx_0_xxx_1, g_xxxxx_0_xxxx_1, g_xxxxx_0_xxxy_1, g_xxxxx_0_xxxz_1, g_xxxxx_0_xxy_1, g_xxxxx_0_xxyy_1, g_xxxxx_0_xxyz_1, g_xxxxx_0_xxz_1, g_xxxxx_0_xxzz_1, g_xxxxx_0_xyy_1, g_xxxxx_0_xyyy_1, g_xxxxx_0_xyyz_1, g_xxxxx_0_xyz_1, g_xxxxx_0_xyzz_1, g_xxxxx_0_xzz_1, g_xxxxx_0_xzzz_1, g_xxxxx_0_yyy_1, g_xxxxx_0_yyyy_1, g_xxxxx_0_yyyz_1, g_xxxxx_0_yyz_1, g_xxxxx_0_yyzz_1, g_xxxxx_0_yzz_1, g_xxxxx_0_yzzz_1, g_xxxxx_0_zzz_1, g_xxxxx_0_zzzz_1, g_xxxxxy_0_xxxx_0, g_xxxxxy_0_xxxy_0, g_xxxxxy_0_xxxz_0, g_xxxxxy_0_xxyy_0, g_xxxxxy_0_xxyz_0, g_xxxxxy_0_xxzz_0, g_xxxxxy_0_xyyy_0, g_xxxxxy_0_xyyz_0, g_xxxxxy_0_xyzz_0, g_xxxxxy_0_xzzz_0, g_xxxxxy_0_yyyy_0, g_xxxxxy_0_yyyz_0, g_xxxxxy_0_yyzz_0, g_xxxxxy_0_yzzz_0, g_xxxxxy_0_zzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxy_0_xxxx_0[i] = g_xxxxx_0_xxxx_1[i] * wa_y[i];

        g_xxxxxy_0_xxxy_0[i] = g_xxxxx_0_xxx_1[i] * fi_acd_0 + g_xxxxx_0_xxxy_1[i] * wa_y[i];

        g_xxxxxy_0_xxxz_0[i] = g_xxxxx_0_xxxz_1[i] * wa_y[i];

        g_xxxxxy_0_xxyy_0[i] = 2.0 * g_xxxxx_0_xxy_1[i] * fi_acd_0 + g_xxxxx_0_xxyy_1[i] * wa_y[i];

        g_xxxxxy_0_xxyz_0[i] = g_xxxxx_0_xxz_1[i] * fi_acd_0 + g_xxxxx_0_xxyz_1[i] * wa_y[i];

        g_xxxxxy_0_xxzz_0[i] = g_xxxxx_0_xxzz_1[i] * wa_y[i];

        g_xxxxxy_0_xyyy_0[i] = 3.0 * g_xxxxx_0_xyy_1[i] * fi_acd_0 + g_xxxxx_0_xyyy_1[i] * wa_y[i];

        g_xxxxxy_0_xyyz_0[i] = 2.0 * g_xxxxx_0_xyz_1[i] * fi_acd_0 + g_xxxxx_0_xyyz_1[i] * wa_y[i];

        g_xxxxxy_0_xyzz_0[i] = g_xxxxx_0_xzz_1[i] * fi_acd_0 + g_xxxxx_0_xyzz_1[i] * wa_y[i];

        g_xxxxxy_0_xzzz_0[i] = g_xxxxx_0_xzzz_1[i] * wa_y[i];

        g_xxxxxy_0_yyyy_0[i] = 4.0 * g_xxxxx_0_yyy_1[i] * fi_acd_0 + g_xxxxx_0_yyyy_1[i] * wa_y[i];

        g_xxxxxy_0_yyyz_0[i] = 3.0 * g_xxxxx_0_yyz_1[i] * fi_acd_0 + g_xxxxx_0_yyyz_1[i] * wa_y[i];

        g_xxxxxy_0_yyzz_0[i] = 2.0 * g_xxxxx_0_yzz_1[i] * fi_acd_0 + g_xxxxx_0_yyzz_1[i] * wa_y[i];

        g_xxxxxy_0_yzzz_0[i] = g_xxxxx_0_zzz_1[i] * fi_acd_0 + g_xxxxx_0_yzzz_1[i] * wa_y[i];

        g_xxxxxy_0_zzzz_0[i] = g_xxxxx_0_zzzz_1[i] * wa_y[i];
    }

    /// Set up 30-45 components of targeted buffer : ISG

    auto g_xxxxxz_0_xxxx_0 = pbuffer.data(idx_eri_0_isg + 30);

    auto g_xxxxxz_0_xxxy_0 = pbuffer.data(idx_eri_0_isg + 31);

    auto g_xxxxxz_0_xxxz_0 = pbuffer.data(idx_eri_0_isg + 32);

    auto g_xxxxxz_0_xxyy_0 = pbuffer.data(idx_eri_0_isg + 33);

    auto g_xxxxxz_0_xxyz_0 = pbuffer.data(idx_eri_0_isg + 34);

    auto g_xxxxxz_0_xxzz_0 = pbuffer.data(idx_eri_0_isg + 35);

    auto g_xxxxxz_0_xyyy_0 = pbuffer.data(idx_eri_0_isg + 36);

    auto g_xxxxxz_0_xyyz_0 = pbuffer.data(idx_eri_0_isg + 37);

    auto g_xxxxxz_0_xyzz_0 = pbuffer.data(idx_eri_0_isg + 38);

    auto g_xxxxxz_0_xzzz_0 = pbuffer.data(idx_eri_0_isg + 39);

    auto g_xxxxxz_0_yyyy_0 = pbuffer.data(idx_eri_0_isg + 40);

    auto g_xxxxxz_0_yyyz_0 = pbuffer.data(idx_eri_0_isg + 41);

    auto g_xxxxxz_0_yyzz_0 = pbuffer.data(idx_eri_0_isg + 42);

    auto g_xxxxxz_0_yzzz_0 = pbuffer.data(idx_eri_0_isg + 43);

    auto g_xxxxxz_0_zzzz_0 = pbuffer.data(idx_eri_0_isg + 44);

    #pragma omp simd aligned(g_xxxxx_0_xxx_1, g_xxxxx_0_xxxx_1, g_xxxxx_0_xxxy_1, g_xxxxx_0_xxxz_1, g_xxxxx_0_xxy_1, g_xxxxx_0_xxyy_1, g_xxxxx_0_xxyz_1, g_xxxxx_0_xxz_1, g_xxxxx_0_xxzz_1, g_xxxxx_0_xyy_1, g_xxxxx_0_xyyy_1, g_xxxxx_0_xyyz_1, g_xxxxx_0_xyz_1, g_xxxxx_0_xyzz_1, g_xxxxx_0_xzz_1, g_xxxxx_0_xzzz_1, g_xxxxx_0_yyy_1, g_xxxxx_0_yyyy_1, g_xxxxx_0_yyyz_1, g_xxxxx_0_yyz_1, g_xxxxx_0_yyzz_1, g_xxxxx_0_yzz_1, g_xxxxx_0_yzzz_1, g_xxxxx_0_zzz_1, g_xxxxx_0_zzzz_1, g_xxxxxz_0_xxxx_0, g_xxxxxz_0_xxxy_0, g_xxxxxz_0_xxxz_0, g_xxxxxz_0_xxyy_0, g_xxxxxz_0_xxyz_0, g_xxxxxz_0_xxzz_0, g_xxxxxz_0_xyyy_0, g_xxxxxz_0_xyyz_0, g_xxxxxz_0_xyzz_0, g_xxxxxz_0_xzzz_0, g_xxxxxz_0_yyyy_0, g_xxxxxz_0_yyyz_0, g_xxxxxz_0_yyzz_0, g_xxxxxz_0_yzzz_0, g_xxxxxz_0_zzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxz_0_xxxx_0[i] = g_xxxxx_0_xxxx_1[i] * wa_z[i];

        g_xxxxxz_0_xxxy_0[i] = g_xxxxx_0_xxxy_1[i] * wa_z[i];

        g_xxxxxz_0_xxxz_0[i] = g_xxxxx_0_xxx_1[i] * fi_acd_0 + g_xxxxx_0_xxxz_1[i] * wa_z[i];

        g_xxxxxz_0_xxyy_0[i] = g_xxxxx_0_xxyy_1[i] * wa_z[i];

        g_xxxxxz_0_xxyz_0[i] = g_xxxxx_0_xxy_1[i] * fi_acd_0 + g_xxxxx_0_xxyz_1[i] * wa_z[i];

        g_xxxxxz_0_xxzz_0[i] = 2.0 * g_xxxxx_0_xxz_1[i] * fi_acd_0 + g_xxxxx_0_xxzz_1[i] * wa_z[i];

        g_xxxxxz_0_xyyy_0[i] = g_xxxxx_0_xyyy_1[i] * wa_z[i];

        g_xxxxxz_0_xyyz_0[i] = g_xxxxx_0_xyy_1[i] * fi_acd_0 + g_xxxxx_0_xyyz_1[i] * wa_z[i];

        g_xxxxxz_0_xyzz_0[i] = 2.0 * g_xxxxx_0_xyz_1[i] * fi_acd_0 + g_xxxxx_0_xyzz_1[i] * wa_z[i];

        g_xxxxxz_0_xzzz_0[i] = 3.0 * g_xxxxx_0_xzz_1[i] * fi_acd_0 + g_xxxxx_0_xzzz_1[i] * wa_z[i];

        g_xxxxxz_0_yyyy_0[i] = g_xxxxx_0_yyyy_1[i] * wa_z[i];

        g_xxxxxz_0_yyyz_0[i] = g_xxxxx_0_yyy_1[i] * fi_acd_0 + g_xxxxx_0_yyyz_1[i] * wa_z[i];

        g_xxxxxz_0_yyzz_0[i] = 2.0 * g_xxxxx_0_yyz_1[i] * fi_acd_0 + g_xxxxx_0_yyzz_1[i] * wa_z[i];

        g_xxxxxz_0_yzzz_0[i] = 3.0 * g_xxxxx_0_yzz_1[i] * fi_acd_0 + g_xxxxx_0_yzzz_1[i] * wa_z[i];

        g_xxxxxz_0_zzzz_0[i] = 4.0 * g_xxxxx_0_zzz_1[i] * fi_acd_0 + g_xxxxx_0_zzzz_1[i] * wa_z[i];
    }

    /// Set up 45-60 components of targeted buffer : ISG

    auto g_xxxxyy_0_xxxx_0 = pbuffer.data(idx_eri_0_isg + 45);

    auto g_xxxxyy_0_xxxy_0 = pbuffer.data(idx_eri_0_isg + 46);

    auto g_xxxxyy_0_xxxz_0 = pbuffer.data(idx_eri_0_isg + 47);

    auto g_xxxxyy_0_xxyy_0 = pbuffer.data(idx_eri_0_isg + 48);

    auto g_xxxxyy_0_xxyz_0 = pbuffer.data(idx_eri_0_isg + 49);

    auto g_xxxxyy_0_xxzz_0 = pbuffer.data(idx_eri_0_isg + 50);

    auto g_xxxxyy_0_xyyy_0 = pbuffer.data(idx_eri_0_isg + 51);

    auto g_xxxxyy_0_xyyz_0 = pbuffer.data(idx_eri_0_isg + 52);

    auto g_xxxxyy_0_xyzz_0 = pbuffer.data(idx_eri_0_isg + 53);

    auto g_xxxxyy_0_xzzz_0 = pbuffer.data(idx_eri_0_isg + 54);

    auto g_xxxxyy_0_yyyy_0 = pbuffer.data(idx_eri_0_isg + 55);

    auto g_xxxxyy_0_yyyz_0 = pbuffer.data(idx_eri_0_isg + 56);

    auto g_xxxxyy_0_yyzz_0 = pbuffer.data(idx_eri_0_isg + 57);

    auto g_xxxxyy_0_yzzz_0 = pbuffer.data(idx_eri_0_isg + 58);

    auto g_xxxxyy_0_zzzz_0 = pbuffer.data(idx_eri_0_isg + 59);

    #pragma omp simd aligned(g_xxxx_0_xxxx_0, g_xxxx_0_xxxx_1, g_xxxx_0_xxxz_0, g_xxxx_0_xxxz_1, g_xxxx_0_xxzz_0, g_xxxx_0_xxzz_1, g_xxxx_0_xzzz_0, g_xxxx_0_xzzz_1, g_xxxxy_0_xxxx_1, g_xxxxy_0_xxxz_1, g_xxxxy_0_xxzz_1, g_xxxxy_0_xzzz_1, g_xxxxyy_0_xxxx_0, g_xxxxyy_0_xxxy_0, g_xxxxyy_0_xxxz_0, g_xxxxyy_0_xxyy_0, g_xxxxyy_0_xxyz_0, g_xxxxyy_0_xxzz_0, g_xxxxyy_0_xyyy_0, g_xxxxyy_0_xyyz_0, g_xxxxyy_0_xyzz_0, g_xxxxyy_0_xzzz_0, g_xxxxyy_0_yyyy_0, g_xxxxyy_0_yyyz_0, g_xxxxyy_0_yyzz_0, g_xxxxyy_0_yzzz_0, g_xxxxyy_0_zzzz_0, g_xxxyy_0_xxxy_1, g_xxxyy_0_xxy_1, g_xxxyy_0_xxyy_1, g_xxxyy_0_xxyz_1, g_xxxyy_0_xyy_1, g_xxxyy_0_xyyy_1, g_xxxyy_0_xyyz_1, g_xxxyy_0_xyz_1, g_xxxyy_0_xyzz_1, g_xxxyy_0_yyy_1, g_xxxyy_0_yyyy_1, g_xxxyy_0_yyyz_1, g_xxxyy_0_yyz_1, g_xxxyy_0_yyzz_1, g_xxxyy_0_yzz_1, g_xxxyy_0_yzzz_1, g_xxxyy_0_zzzz_1, g_xxyy_0_xxxy_0, g_xxyy_0_xxxy_1, g_xxyy_0_xxyy_0, g_xxyy_0_xxyy_1, g_xxyy_0_xxyz_0, g_xxyy_0_xxyz_1, g_xxyy_0_xyyy_0, g_xxyy_0_xyyy_1, g_xxyy_0_xyyz_0, g_xxyy_0_xyyz_1, g_xxyy_0_xyzz_0, g_xxyy_0_xyzz_1, g_xxyy_0_yyyy_0, g_xxyy_0_yyyy_1, g_xxyy_0_yyyz_0, g_xxyy_0_yyyz_1, g_xxyy_0_yyzz_0, g_xxyy_0_yyzz_1, g_xxyy_0_yzzz_0, g_xxyy_0_yzzz_1, g_xxyy_0_zzzz_0, g_xxyy_0_zzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxyy_0_xxxx_0[i] = g_xxxx_0_xxxx_0[i] * fbe_0 - g_xxxx_0_xxxx_1[i] * fz_be_0 + g_xxxxy_0_xxxx_1[i] * wa_y[i];

        g_xxxxyy_0_xxxy_0[i] = 3.0 * g_xxyy_0_xxxy_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxy_1[i] * fz_be_0 + 3.0 * g_xxxyy_0_xxy_1[i] * fi_acd_0 + g_xxxyy_0_xxxy_1[i] * wa_x[i];

        g_xxxxyy_0_xxxz_0[i] = g_xxxx_0_xxxz_0[i] * fbe_0 - g_xxxx_0_xxxz_1[i] * fz_be_0 + g_xxxxy_0_xxxz_1[i] * wa_y[i];

        g_xxxxyy_0_xxyy_0[i] = 3.0 * g_xxyy_0_xxyy_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxyy_1[i] * fz_be_0 + 2.0 * g_xxxyy_0_xyy_1[i] * fi_acd_0 + g_xxxyy_0_xxyy_1[i] * wa_x[i];

        g_xxxxyy_0_xxyz_0[i] = 3.0 * g_xxyy_0_xxyz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxyz_1[i] * fz_be_0 + 2.0 * g_xxxyy_0_xyz_1[i] * fi_acd_0 + g_xxxyy_0_xxyz_1[i] * wa_x[i];

        g_xxxxyy_0_xxzz_0[i] = g_xxxx_0_xxzz_0[i] * fbe_0 - g_xxxx_0_xxzz_1[i] * fz_be_0 + g_xxxxy_0_xxzz_1[i] * wa_y[i];

        g_xxxxyy_0_xyyy_0[i] = 3.0 * g_xxyy_0_xyyy_0[i] * fbe_0 - 3.0 * g_xxyy_0_xyyy_1[i] * fz_be_0 + g_xxxyy_0_yyy_1[i] * fi_acd_0 + g_xxxyy_0_xyyy_1[i] * wa_x[i];

        g_xxxxyy_0_xyyz_0[i] = 3.0 * g_xxyy_0_xyyz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xyyz_1[i] * fz_be_0 + g_xxxyy_0_yyz_1[i] * fi_acd_0 + g_xxxyy_0_xyyz_1[i] * wa_x[i];

        g_xxxxyy_0_xyzz_0[i] = 3.0 * g_xxyy_0_xyzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xyzz_1[i] * fz_be_0 + g_xxxyy_0_yzz_1[i] * fi_acd_0 + g_xxxyy_0_xyzz_1[i] * wa_x[i];

        g_xxxxyy_0_xzzz_0[i] = g_xxxx_0_xzzz_0[i] * fbe_0 - g_xxxx_0_xzzz_1[i] * fz_be_0 + g_xxxxy_0_xzzz_1[i] * wa_y[i];

        g_xxxxyy_0_yyyy_0[i] = 3.0 * g_xxyy_0_yyyy_0[i] * fbe_0 - 3.0 * g_xxyy_0_yyyy_1[i] * fz_be_0 + g_xxxyy_0_yyyy_1[i] * wa_x[i];

        g_xxxxyy_0_yyyz_0[i] = 3.0 * g_xxyy_0_yyyz_0[i] * fbe_0 - 3.0 * g_xxyy_0_yyyz_1[i] * fz_be_0 + g_xxxyy_0_yyyz_1[i] * wa_x[i];

        g_xxxxyy_0_yyzz_0[i] = 3.0 * g_xxyy_0_yyzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_yyzz_1[i] * fz_be_0 + g_xxxyy_0_yyzz_1[i] * wa_x[i];

        g_xxxxyy_0_yzzz_0[i] = 3.0 * g_xxyy_0_yzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_yzzz_1[i] * fz_be_0 + g_xxxyy_0_yzzz_1[i] * wa_x[i];

        g_xxxxyy_0_zzzz_0[i] = 3.0 * g_xxyy_0_zzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_zzzz_1[i] * fz_be_0 + g_xxxyy_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 60-75 components of targeted buffer : ISG

    auto g_xxxxyz_0_xxxx_0 = pbuffer.data(idx_eri_0_isg + 60);

    auto g_xxxxyz_0_xxxy_0 = pbuffer.data(idx_eri_0_isg + 61);

    auto g_xxxxyz_0_xxxz_0 = pbuffer.data(idx_eri_0_isg + 62);

    auto g_xxxxyz_0_xxyy_0 = pbuffer.data(idx_eri_0_isg + 63);

    auto g_xxxxyz_0_xxyz_0 = pbuffer.data(idx_eri_0_isg + 64);

    auto g_xxxxyz_0_xxzz_0 = pbuffer.data(idx_eri_0_isg + 65);

    auto g_xxxxyz_0_xyyy_0 = pbuffer.data(idx_eri_0_isg + 66);

    auto g_xxxxyz_0_xyyz_0 = pbuffer.data(idx_eri_0_isg + 67);

    auto g_xxxxyz_0_xyzz_0 = pbuffer.data(idx_eri_0_isg + 68);

    auto g_xxxxyz_0_xzzz_0 = pbuffer.data(idx_eri_0_isg + 69);

    auto g_xxxxyz_0_yyyy_0 = pbuffer.data(idx_eri_0_isg + 70);

    auto g_xxxxyz_0_yyyz_0 = pbuffer.data(idx_eri_0_isg + 71);

    auto g_xxxxyz_0_yyzz_0 = pbuffer.data(idx_eri_0_isg + 72);

    auto g_xxxxyz_0_yzzz_0 = pbuffer.data(idx_eri_0_isg + 73);

    auto g_xxxxyz_0_zzzz_0 = pbuffer.data(idx_eri_0_isg + 74);

    #pragma omp simd aligned(g_xxxxy_0_xxxy_1, g_xxxxy_0_xxyy_1, g_xxxxy_0_xyyy_1, g_xxxxy_0_yyyy_1, g_xxxxyz_0_xxxx_0, g_xxxxyz_0_xxxy_0, g_xxxxyz_0_xxxz_0, g_xxxxyz_0_xxyy_0, g_xxxxyz_0_xxyz_0, g_xxxxyz_0_xxzz_0, g_xxxxyz_0_xyyy_0, g_xxxxyz_0_xyyz_0, g_xxxxyz_0_xyzz_0, g_xxxxyz_0_xzzz_0, g_xxxxyz_0_yyyy_0, g_xxxxyz_0_yyyz_0, g_xxxxyz_0_yyzz_0, g_xxxxyz_0_yzzz_0, g_xxxxyz_0_zzzz_0, g_xxxxz_0_xxxx_1, g_xxxxz_0_xxxz_1, g_xxxxz_0_xxyz_1, g_xxxxz_0_xxz_1, g_xxxxz_0_xxzz_1, g_xxxxz_0_xyyz_1, g_xxxxz_0_xyz_1, g_xxxxz_0_xyzz_1, g_xxxxz_0_xzz_1, g_xxxxz_0_xzzz_1, g_xxxxz_0_yyyz_1, g_xxxxz_0_yyz_1, g_xxxxz_0_yyzz_1, g_xxxxz_0_yzz_1, g_xxxxz_0_yzzz_1, g_xxxxz_0_zzz_1, g_xxxxz_0_zzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxyz_0_xxxx_0[i] = g_xxxxz_0_xxxx_1[i] * wa_y[i];

        g_xxxxyz_0_xxxy_0[i] = g_xxxxy_0_xxxy_1[i] * wa_z[i];

        g_xxxxyz_0_xxxz_0[i] = g_xxxxz_0_xxxz_1[i] * wa_y[i];

        g_xxxxyz_0_xxyy_0[i] = g_xxxxy_0_xxyy_1[i] * wa_z[i];

        g_xxxxyz_0_xxyz_0[i] = g_xxxxz_0_xxz_1[i] * fi_acd_0 + g_xxxxz_0_xxyz_1[i] * wa_y[i];

        g_xxxxyz_0_xxzz_0[i] = g_xxxxz_0_xxzz_1[i] * wa_y[i];

        g_xxxxyz_0_xyyy_0[i] = g_xxxxy_0_xyyy_1[i] * wa_z[i];

        g_xxxxyz_0_xyyz_0[i] = 2.0 * g_xxxxz_0_xyz_1[i] * fi_acd_0 + g_xxxxz_0_xyyz_1[i] * wa_y[i];

        g_xxxxyz_0_xyzz_0[i] = g_xxxxz_0_xzz_1[i] * fi_acd_0 + g_xxxxz_0_xyzz_1[i] * wa_y[i];

        g_xxxxyz_0_xzzz_0[i] = g_xxxxz_0_xzzz_1[i] * wa_y[i];

        g_xxxxyz_0_yyyy_0[i] = g_xxxxy_0_yyyy_1[i] * wa_z[i];

        g_xxxxyz_0_yyyz_0[i] = 3.0 * g_xxxxz_0_yyz_1[i] * fi_acd_0 + g_xxxxz_0_yyyz_1[i] * wa_y[i];

        g_xxxxyz_0_yyzz_0[i] = 2.0 * g_xxxxz_0_yzz_1[i] * fi_acd_0 + g_xxxxz_0_yyzz_1[i] * wa_y[i];

        g_xxxxyz_0_yzzz_0[i] = g_xxxxz_0_zzz_1[i] * fi_acd_0 + g_xxxxz_0_yzzz_1[i] * wa_y[i];

        g_xxxxyz_0_zzzz_0[i] = g_xxxxz_0_zzzz_1[i] * wa_y[i];
    }

    /// Set up 75-90 components of targeted buffer : ISG

    auto g_xxxxzz_0_xxxx_0 = pbuffer.data(idx_eri_0_isg + 75);

    auto g_xxxxzz_0_xxxy_0 = pbuffer.data(idx_eri_0_isg + 76);

    auto g_xxxxzz_0_xxxz_0 = pbuffer.data(idx_eri_0_isg + 77);

    auto g_xxxxzz_0_xxyy_0 = pbuffer.data(idx_eri_0_isg + 78);

    auto g_xxxxzz_0_xxyz_0 = pbuffer.data(idx_eri_0_isg + 79);

    auto g_xxxxzz_0_xxzz_0 = pbuffer.data(idx_eri_0_isg + 80);

    auto g_xxxxzz_0_xyyy_0 = pbuffer.data(idx_eri_0_isg + 81);

    auto g_xxxxzz_0_xyyz_0 = pbuffer.data(idx_eri_0_isg + 82);

    auto g_xxxxzz_0_xyzz_0 = pbuffer.data(idx_eri_0_isg + 83);

    auto g_xxxxzz_0_xzzz_0 = pbuffer.data(idx_eri_0_isg + 84);

    auto g_xxxxzz_0_yyyy_0 = pbuffer.data(idx_eri_0_isg + 85);

    auto g_xxxxzz_0_yyyz_0 = pbuffer.data(idx_eri_0_isg + 86);

    auto g_xxxxzz_0_yyzz_0 = pbuffer.data(idx_eri_0_isg + 87);

    auto g_xxxxzz_0_yzzz_0 = pbuffer.data(idx_eri_0_isg + 88);

    auto g_xxxxzz_0_zzzz_0 = pbuffer.data(idx_eri_0_isg + 89);

    #pragma omp simd aligned(g_xxxx_0_xxxx_0, g_xxxx_0_xxxx_1, g_xxxx_0_xxxy_0, g_xxxx_0_xxxy_1, g_xxxx_0_xxyy_0, g_xxxx_0_xxyy_1, g_xxxx_0_xyyy_0, g_xxxx_0_xyyy_1, g_xxxxz_0_xxxx_1, g_xxxxz_0_xxxy_1, g_xxxxz_0_xxyy_1, g_xxxxz_0_xyyy_1, g_xxxxzz_0_xxxx_0, g_xxxxzz_0_xxxy_0, g_xxxxzz_0_xxxz_0, g_xxxxzz_0_xxyy_0, g_xxxxzz_0_xxyz_0, g_xxxxzz_0_xxzz_0, g_xxxxzz_0_xyyy_0, g_xxxxzz_0_xyyz_0, g_xxxxzz_0_xyzz_0, g_xxxxzz_0_xzzz_0, g_xxxxzz_0_yyyy_0, g_xxxxzz_0_yyyz_0, g_xxxxzz_0_yyzz_0, g_xxxxzz_0_yzzz_0, g_xxxxzz_0_zzzz_0, g_xxxzz_0_xxxz_1, g_xxxzz_0_xxyz_1, g_xxxzz_0_xxz_1, g_xxxzz_0_xxzz_1, g_xxxzz_0_xyyz_1, g_xxxzz_0_xyz_1, g_xxxzz_0_xyzz_1, g_xxxzz_0_xzz_1, g_xxxzz_0_xzzz_1, g_xxxzz_0_yyyy_1, g_xxxzz_0_yyyz_1, g_xxxzz_0_yyz_1, g_xxxzz_0_yyzz_1, g_xxxzz_0_yzz_1, g_xxxzz_0_yzzz_1, g_xxxzz_0_zzz_1, g_xxxzz_0_zzzz_1, g_xxzz_0_xxxz_0, g_xxzz_0_xxxz_1, g_xxzz_0_xxyz_0, g_xxzz_0_xxyz_1, g_xxzz_0_xxzz_0, g_xxzz_0_xxzz_1, g_xxzz_0_xyyz_0, g_xxzz_0_xyyz_1, g_xxzz_0_xyzz_0, g_xxzz_0_xyzz_1, g_xxzz_0_xzzz_0, g_xxzz_0_xzzz_1, g_xxzz_0_yyyy_0, g_xxzz_0_yyyy_1, g_xxzz_0_yyyz_0, g_xxzz_0_yyyz_1, g_xxzz_0_yyzz_0, g_xxzz_0_yyzz_1, g_xxzz_0_yzzz_0, g_xxzz_0_yzzz_1, g_xxzz_0_zzzz_0, g_xxzz_0_zzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxzz_0_xxxx_0[i] = g_xxxx_0_xxxx_0[i] * fbe_0 - g_xxxx_0_xxxx_1[i] * fz_be_0 + g_xxxxz_0_xxxx_1[i] * wa_z[i];

        g_xxxxzz_0_xxxy_0[i] = g_xxxx_0_xxxy_0[i] * fbe_0 - g_xxxx_0_xxxy_1[i] * fz_be_0 + g_xxxxz_0_xxxy_1[i] * wa_z[i];

        g_xxxxzz_0_xxxz_0[i] = 3.0 * g_xxzz_0_xxxz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxz_1[i] * fz_be_0 + 3.0 * g_xxxzz_0_xxz_1[i] * fi_acd_0 + g_xxxzz_0_xxxz_1[i] * wa_x[i];

        g_xxxxzz_0_xxyy_0[i] = g_xxxx_0_xxyy_0[i] * fbe_0 - g_xxxx_0_xxyy_1[i] * fz_be_0 + g_xxxxz_0_xxyy_1[i] * wa_z[i];

        g_xxxxzz_0_xxyz_0[i] = 3.0 * g_xxzz_0_xxyz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxyz_1[i] * fz_be_0 + 2.0 * g_xxxzz_0_xyz_1[i] * fi_acd_0 + g_xxxzz_0_xxyz_1[i] * wa_x[i];

        g_xxxxzz_0_xxzz_0[i] = 3.0 * g_xxzz_0_xxzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxzz_1[i] * fz_be_0 + 2.0 * g_xxxzz_0_xzz_1[i] * fi_acd_0 + g_xxxzz_0_xxzz_1[i] * wa_x[i];

        g_xxxxzz_0_xyyy_0[i] = g_xxxx_0_xyyy_0[i] * fbe_0 - g_xxxx_0_xyyy_1[i] * fz_be_0 + g_xxxxz_0_xyyy_1[i] * wa_z[i];

        g_xxxxzz_0_xyyz_0[i] = 3.0 * g_xxzz_0_xyyz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xyyz_1[i] * fz_be_0 + g_xxxzz_0_yyz_1[i] * fi_acd_0 + g_xxxzz_0_xyyz_1[i] * wa_x[i];

        g_xxxxzz_0_xyzz_0[i] = 3.0 * g_xxzz_0_xyzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xyzz_1[i] * fz_be_0 + g_xxxzz_0_yzz_1[i] * fi_acd_0 + g_xxxzz_0_xyzz_1[i] * wa_x[i];

        g_xxxxzz_0_xzzz_0[i] = 3.0 * g_xxzz_0_xzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xzzz_1[i] * fz_be_0 + g_xxxzz_0_zzz_1[i] * fi_acd_0 + g_xxxzz_0_xzzz_1[i] * wa_x[i];

        g_xxxxzz_0_yyyy_0[i] = 3.0 * g_xxzz_0_yyyy_0[i] * fbe_0 - 3.0 * g_xxzz_0_yyyy_1[i] * fz_be_0 + g_xxxzz_0_yyyy_1[i] * wa_x[i];

        g_xxxxzz_0_yyyz_0[i] = 3.0 * g_xxzz_0_yyyz_0[i] * fbe_0 - 3.0 * g_xxzz_0_yyyz_1[i] * fz_be_0 + g_xxxzz_0_yyyz_1[i] * wa_x[i];

        g_xxxxzz_0_yyzz_0[i] = 3.0 * g_xxzz_0_yyzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_yyzz_1[i] * fz_be_0 + g_xxxzz_0_yyzz_1[i] * wa_x[i];

        g_xxxxzz_0_yzzz_0[i] = 3.0 * g_xxzz_0_yzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_yzzz_1[i] * fz_be_0 + g_xxxzz_0_yzzz_1[i] * wa_x[i];

        g_xxxxzz_0_zzzz_0[i] = 3.0 * g_xxzz_0_zzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_zzzz_1[i] * fz_be_0 + g_xxxzz_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 90-105 components of targeted buffer : ISG

    auto g_xxxyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_isg + 90);

    auto g_xxxyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_isg + 91);

    auto g_xxxyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_isg + 92);

    auto g_xxxyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_isg + 93);

    auto g_xxxyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_isg + 94);

    auto g_xxxyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_isg + 95);

    auto g_xxxyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_isg + 96);

    auto g_xxxyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_isg + 97);

    auto g_xxxyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_isg + 98);

    auto g_xxxyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_isg + 99);

    auto g_xxxyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_isg + 100);

    auto g_xxxyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_isg + 101);

    auto g_xxxyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_isg + 102);

    auto g_xxxyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_isg + 103);

    auto g_xxxyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_isg + 104);

    #pragma omp simd aligned(g_xxxy_0_xxxx_0, g_xxxy_0_xxxx_1, g_xxxy_0_xxxz_0, g_xxxy_0_xxxz_1, g_xxxy_0_xxzz_0, g_xxxy_0_xxzz_1, g_xxxy_0_xzzz_0, g_xxxy_0_xzzz_1, g_xxxyy_0_xxxx_1, g_xxxyy_0_xxxz_1, g_xxxyy_0_xxzz_1, g_xxxyy_0_xzzz_1, g_xxxyyy_0_xxxx_0, g_xxxyyy_0_xxxy_0, g_xxxyyy_0_xxxz_0, g_xxxyyy_0_xxyy_0, g_xxxyyy_0_xxyz_0, g_xxxyyy_0_xxzz_0, g_xxxyyy_0_xyyy_0, g_xxxyyy_0_xyyz_0, g_xxxyyy_0_xyzz_0, g_xxxyyy_0_xzzz_0, g_xxxyyy_0_yyyy_0, g_xxxyyy_0_yyyz_0, g_xxxyyy_0_yyzz_0, g_xxxyyy_0_yzzz_0, g_xxxyyy_0_zzzz_0, g_xxyyy_0_xxxy_1, g_xxyyy_0_xxy_1, g_xxyyy_0_xxyy_1, g_xxyyy_0_xxyz_1, g_xxyyy_0_xyy_1, g_xxyyy_0_xyyy_1, g_xxyyy_0_xyyz_1, g_xxyyy_0_xyz_1, g_xxyyy_0_xyzz_1, g_xxyyy_0_yyy_1, g_xxyyy_0_yyyy_1, g_xxyyy_0_yyyz_1, g_xxyyy_0_yyz_1, g_xxyyy_0_yyzz_1, g_xxyyy_0_yzz_1, g_xxyyy_0_yzzz_1, g_xxyyy_0_zzzz_1, g_xyyy_0_xxxy_0, g_xyyy_0_xxxy_1, g_xyyy_0_xxyy_0, g_xyyy_0_xxyy_1, g_xyyy_0_xxyz_0, g_xyyy_0_xxyz_1, g_xyyy_0_xyyy_0, g_xyyy_0_xyyy_1, g_xyyy_0_xyyz_0, g_xyyy_0_xyyz_1, g_xyyy_0_xyzz_0, g_xyyy_0_xyzz_1, g_xyyy_0_yyyy_0, g_xyyy_0_yyyy_1, g_xyyy_0_yyyz_0, g_xyyy_0_yyyz_1, g_xyyy_0_yyzz_0, g_xyyy_0_yyzz_1, g_xyyy_0_yzzz_0, g_xyyy_0_yzzz_1, g_xyyy_0_zzzz_0, g_xyyy_0_zzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxyyy_0_xxxx_0[i] = 2.0 * g_xxxy_0_xxxx_0[i] * fbe_0 - 2.0 * g_xxxy_0_xxxx_1[i] * fz_be_0 + g_xxxyy_0_xxxx_1[i] * wa_y[i];

        g_xxxyyy_0_xxxy_0[i] = 2.0 * g_xyyy_0_xxxy_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxxy_1[i] * fz_be_0 + 3.0 * g_xxyyy_0_xxy_1[i] * fi_acd_0 + g_xxyyy_0_xxxy_1[i] * wa_x[i];

        g_xxxyyy_0_xxxz_0[i] = 2.0 * g_xxxy_0_xxxz_0[i] * fbe_0 - 2.0 * g_xxxy_0_xxxz_1[i] * fz_be_0 + g_xxxyy_0_xxxz_1[i] * wa_y[i];

        g_xxxyyy_0_xxyy_0[i] = 2.0 * g_xyyy_0_xxyy_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxyy_1[i] * fz_be_0 + 2.0 * g_xxyyy_0_xyy_1[i] * fi_acd_0 + g_xxyyy_0_xxyy_1[i] * wa_x[i];

        g_xxxyyy_0_xxyz_0[i] = 2.0 * g_xyyy_0_xxyz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxyz_1[i] * fz_be_0 + 2.0 * g_xxyyy_0_xyz_1[i] * fi_acd_0 + g_xxyyy_0_xxyz_1[i] * wa_x[i];

        g_xxxyyy_0_xxzz_0[i] = 2.0 * g_xxxy_0_xxzz_0[i] * fbe_0 - 2.0 * g_xxxy_0_xxzz_1[i] * fz_be_0 + g_xxxyy_0_xxzz_1[i] * wa_y[i];

        g_xxxyyy_0_xyyy_0[i] = 2.0 * g_xyyy_0_xyyy_0[i] * fbe_0 - 2.0 * g_xyyy_0_xyyy_1[i] * fz_be_0 + g_xxyyy_0_yyy_1[i] * fi_acd_0 + g_xxyyy_0_xyyy_1[i] * wa_x[i];

        g_xxxyyy_0_xyyz_0[i] = 2.0 * g_xyyy_0_xyyz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xyyz_1[i] * fz_be_0 + g_xxyyy_0_yyz_1[i] * fi_acd_0 + g_xxyyy_0_xyyz_1[i] * wa_x[i];

        g_xxxyyy_0_xyzz_0[i] = 2.0 * g_xyyy_0_xyzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xyzz_1[i] * fz_be_0 + g_xxyyy_0_yzz_1[i] * fi_acd_0 + g_xxyyy_0_xyzz_1[i] * wa_x[i];

        g_xxxyyy_0_xzzz_0[i] = 2.0 * g_xxxy_0_xzzz_0[i] * fbe_0 - 2.0 * g_xxxy_0_xzzz_1[i] * fz_be_0 + g_xxxyy_0_xzzz_1[i] * wa_y[i];

        g_xxxyyy_0_yyyy_0[i] = 2.0 * g_xyyy_0_yyyy_0[i] * fbe_0 - 2.0 * g_xyyy_0_yyyy_1[i] * fz_be_0 + g_xxyyy_0_yyyy_1[i] * wa_x[i];

        g_xxxyyy_0_yyyz_0[i] = 2.0 * g_xyyy_0_yyyz_0[i] * fbe_0 - 2.0 * g_xyyy_0_yyyz_1[i] * fz_be_0 + g_xxyyy_0_yyyz_1[i] * wa_x[i];

        g_xxxyyy_0_yyzz_0[i] = 2.0 * g_xyyy_0_yyzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_yyzz_1[i] * fz_be_0 + g_xxyyy_0_yyzz_1[i] * wa_x[i];

        g_xxxyyy_0_yzzz_0[i] = 2.0 * g_xyyy_0_yzzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_yzzz_1[i] * fz_be_0 + g_xxyyy_0_yzzz_1[i] * wa_x[i];

        g_xxxyyy_0_zzzz_0[i] = 2.0 * g_xyyy_0_zzzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_zzzz_1[i] * fz_be_0 + g_xxyyy_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 105-120 components of targeted buffer : ISG

    auto g_xxxyyz_0_xxxx_0 = pbuffer.data(idx_eri_0_isg + 105);

    auto g_xxxyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_isg + 106);

    auto g_xxxyyz_0_xxxz_0 = pbuffer.data(idx_eri_0_isg + 107);

    auto g_xxxyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_isg + 108);

    auto g_xxxyyz_0_xxyz_0 = pbuffer.data(idx_eri_0_isg + 109);

    auto g_xxxyyz_0_xxzz_0 = pbuffer.data(idx_eri_0_isg + 110);

    auto g_xxxyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_isg + 111);

    auto g_xxxyyz_0_xyyz_0 = pbuffer.data(idx_eri_0_isg + 112);

    auto g_xxxyyz_0_xyzz_0 = pbuffer.data(idx_eri_0_isg + 113);

    auto g_xxxyyz_0_xzzz_0 = pbuffer.data(idx_eri_0_isg + 114);

    auto g_xxxyyz_0_yyyy_0 = pbuffer.data(idx_eri_0_isg + 115);

    auto g_xxxyyz_0_yyyz_0 = pbuffer.data(idx_eri_0_isg + 116);

    auto g_xxxyyz_0_yyzz_0 = pbuffer.data(idx_eri_0_isg + 117);

    auto g_xxxyyz_0_yzzz_0 = pbuffer.data(idx_eri_0_isg + 118);

    auto g_xxxyyz_0_zzzz_0 = pbuffer.data(idx_eri_0_isg + 119);

    #pragma omp simd aligned(g_xxxyy_0_xxx_1, g_xxxyy_0_xxxx_1, g_xxxyy_0_xxxy_1, g_xxxyy_0_xxxz_1, g_xxxyy_0_xxy_1, g_xxxyy_0_xxyy_1, g_xxxyy_0_xxyz_1, g_xxxyy_0_xxz_1, g_xxxyy_0_xxzz_1, g_xxxyy_0_xyy_1, g_xxxyy_0_xyyy_1, g_xxxyy_0_xyyz_1, g_xxxyy_0_xyz_1, g_xxxyy_0_xyzz_1, g_xxxyy_0_xzz_1, g_xxxyy_0_xzzz_1, g_xxxyy_0_yyy_1, g_xxxyy_0_yyyy_1, g_xxxyy_0_yyyz_1, g_xxxyy_0_yyz_1, g_xxxyy_0_yyzz_1, g_xxxyy_0_yzz_1, g_xxxyy_0_yzzz_1, g_xxxyy_0_zzz_1, g_xxxyy_0_zzzz_1, g_xxxyyz_0_xxxx_0, g_xxxyyz_0_xxxy_0, g_xxxyyz_0_xxxz_0, g_xxxyyz_0_xxyy_0, g_xxxyyz_0_xxyz_0, g_xxxyyz_0_xxzz_0, g_xxxyyz_0_xyyy_0, g_xxxyyz_0_xyyz_0, g_xxxyyz_0_xyzz_0, g_xxxyyz_0_xzzz_0, g_xxxyyz_0_yyyy_0, g_xxxyyz_0_yyyz_0, g_xxxyyz_0_yyzz_0, g_xxxyyz_0_yzzz_0, g_xxxyyz_0_zzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyyz_0_xxxx_0[i] = g_xxxyy_0_xxxx_1[i] * wa_z[i];

        g_xxxyyz_0_xxxy_0[i] = g_xxxyy_0_xxxy_1[i] * wa_z[i];

        g_xxxyyz_0_xxxz_0[i] = g_xxxyy_0_xxx_1[i] * fi_acd_0 + g_xxxyy_0_xxxz_1[i] * wa_z[i];

        g_xxxyyz_0_xxyy_0[i] = g_xxxyy_0_xxyy_1[i] * wa_z[i];

        g_xxxyyz_0_xxyz_0[i] = g_xxxyy_0_xxy_1[i] * fi_acd_0 + g_xxxyy_0_xxyz_1[i] * wa_z[i];

        g_xxxyyz_0_xxzz_0[i] = 2.0 * g_xxxyy_0_xxz_1[i] * fi_acd_0 + g_xxxyy_0_xxzz_1[i] * wa_z[i];

        g_xxxyyz_0_xyyy_0[i] = g_xxxyy_0_xyyy_1[i] * wa_z[i];

        g_xxxyyz_0_xyyz_0[i] = g_xxxyy_0_xyy_1[i] * fi_acd_0 + g_xxxyy_0_xyyz_1[i] * wa_z[i];

        g_xxxyyz_0_xyzz_0[i] = 2.0 * g_xxxyy_0_xyz_1[i] * fi_acd_0 + g_xxxyy_0_xyzz_1[i] * wa_z[i];

        g_xxxyyz_0_xzzz_0[i] = 3.0 * g_xxxyy_0_xzz_1[i] * fi_acd_0 + g_xxxyy_0_xzzz_1[i] * wa_z[i];

        g_xxxyyz_0_yyyy_0[i] = g_xxxyy_0_yyyy_1[i] * wa_z[i];

        g_xxxyyz_0_yyyz_0[i] = g_xxxyy_0_yyy_1[i] * fi_acd_0 + g_xxxyy_0_yyyz_1[i] * wa_z[i];

        g_xxxyyz_0_yyzz_0[i] = 2.0 * g_xxxyy_0_yyz_1[i] * fi_acd_0 + g_xxxyy_0_yyzz_1[i] * wa_z[i];

        g_xxxyyz_0_yzzz_0[i] = 3.0 * g_xxxyy_0_yzz_1[i] * fi_acd_0 + g_xxxyy_0_yzzz_1[i] * wa_z[i];

        g_xxxyyz_0_zzzz_0[i] = 4.0 * g_xxxyy_0_zzz_1[i] * fi_acd_0 + g_xxxyy_0_zzzz_1[i] * wa_z[i];
    }

    /// Set up 120-135 components of targeted buffer : ISG

    auto g_xxxyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_isg + 120);

    auto g_xxxyzz_0_xxxy_0 = pbuffer.data(idx_eri_0_isg + 121);

    auto g_xxxyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_isg + 122);

    auto g_xxxyzz_0_xxyy_0 = pbuffer.data(idx_eri_0_isg + 123);

    auto g_xxxyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_isg + 124);

    auto g_xxxyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_isg + 125);

    auto g_xxxyzz_0_xyyy_0 = pbuffer.data(idx_eri_0_isg + 126);

    auto g_xxxyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_isg + 127);

    auto g_xxxyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_isg + 128);

    auto g_xxxyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_isg + 129);

    auto g_xxxyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_isg + 130);

    auto g_xxxyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_isg + 131);

    auto g_xxxyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_isg + 132);

    auto g_xxxyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_isg + 133);

    auto g_xxxyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_isg + 134);

    #pragma omp simd aligned(g_xxxyzz_0_xxxx_0, g_xxxyzz_0_xxxy_0, g_xxxyzz_0_xxxz_0, g_xxxyzz_0_xxyy_0, g_xxxyzz_0_xxyz_0, g_xxxyzz_0_xxzz_0, g_xxxyzz_0_xyyy_0, g_xxxyzz_0_xyyz_0, g_xxxyzz_0_xyzz_0, g_xxxyzz_0_xzzz_0, g_xxxyzz_0_yyyy_0, g_xxxyzz_0_yyyz_0, g_xxxyzz_0_yyzz_0, g_xxxyzz_0_yzzz_0, g_xxxyzz_0_zzzz_0, g_xxxzz_0_xxx_1, g_xxxzz_0_xxxx_1, g_xxxzz_0_xxxy_1, g_xxxzz_0_xxxz_1, g_xxxzz_0_xxy_1, g_xxxzz_0_xxyy_1, g_xxxzz_0_xxyz_1, g_xxxzz_0_xxz_1, g_xxxzz_0_xxzz_1, g_xxxzz_0_xyy_1, g_xxxzz_0_xyyy_1, g_xxxzz_0_xyyz_1, g_xxxzz_0_xyz_1, g_xxxzz_0_xyzz_1, g_xxxzz_0_xzz_1, g_xxxzz_0_xzzz_1, g_xxxzz_0_yyy_1, g_xxxzz_0_yyyy_1, g_xxxzz_0_yyyz_1, g_xxxzz_0_yyz_1, g_xxxzz_0_yyzz_1, g_xxxzz_0_yzz_1, g_xxxzz_0_yzzz_1, g_xxxzz_0_zzz_1, g_xxxzz_0_zzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyzz_0_xxxx_0[i] = g_xxxzz_0_xxxx_1[i] * wa_y[i];

        g_xxxyzz_0_xxxy_0[i] = g_xxxzz_0_xxx_1[i] * fi_acd_0 + g_xxxzz_0_xxxy_1[i] * wa_y[i];

        g_xxxyzz_0_xxxz_0[i] = g_xxxzz_0_xxxz_1[i] * wa_y[i];

        g_xxxyzz_0_xxyy_0[i] = 2.0 * g_xxxzz_0_xxy_1[i] * fi_acd_0 + g_xxxzz_0_xxyy_1[i] * wa_y[i];

        g_xxxyzz_0_xxyz_0[i] = g_xxxzz_0_xxz_1[i] * fi_acd_0 + g_xxxzz_0_xxyz_1[i] * wa_y[i];

        g_xxxyzz_0_xxzz_0[i] = g_xxxzz_0_xxzz_1[i] * wa_y[i];

        g_xxxyzz_0_xyyy_0[i] = 3.0 * g_xxxzz_0_xyy_1[i] * fi_acd_0 + g_xxxzz_0_xyyy_1[i] * wa_y[i];

        g_xxxyzz_0_xyyz_0[i] = 2.0 * g_xxxzz_0_xyz_1[i] * fi_acd_0 + g_xxxzz_0_xyyz_1[i] * wa_y[i];

        g_xxxyzz_0_xyzz_0[i] = g_xxxzz_0_xzz_1[i] * fi_acd_0 + g_xxxzz_0_xyzz_1[i] * wa_y[i];

        g_xxxyzz_0_xzzz_0[i] = g_xxxzz_0_xzzz_1[i] * wa_y[i];

        g_xxxyzz_0_yyyy_0[i] = 4.0 * g_xxxzz_0_yyy_1[i] * fi_acd_0 + g_xxxzz_0_yyyy_1[i] * wa_y[i];

        g_xxxyzz_0_yyyz_0[i] = 3.0 * g_xxxzz_0_yyz_1[i] * fi_acd_0 + g_xxxzz_0_yyyz_1[i] * wa_y[i];

        g_xxxyzz_0_yyzz_0[i] = 2.0 * g_xxxzz_0_yzz_1[i] * fi_acd_0 + g_xxxzz_0_yyzz_1[i] * wa_y[i];

        g_xxxyzz_0_yzzz_0[i] = g_xxxzz_0_zzz_1[i] * fi_acd_0 + g_xxxzz_0_yzzz_1[i] * wa_y[i];

        g_xxxyzz_0_zzzz_0[i] = g_xxxzz_0_zzzz_1[i] * wa_y[i];
    }

    /// Set up 135-150 components of targeted buffer : ISG

    auto g_xxxzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_isg + 135);

    auto g_xxxzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_isg + 136);

    auto g_xxxzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_isg + 137);

    auto g_xxxzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_isg + 138);

    auto g_xxxzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_isg + 139);

    auto g_xxxzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_isg + 140);

    auto g_xxxzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_isg + 141);

    auto g_xxxzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_isg + 142);

    auto g_xxxzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_isg + 143);

    auto g_xxxzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_isg + 144);

    auto g_xxxzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_isg + 145);

    auto g_xxxzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_isg + 146);

    auto g_xxxzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_isg + 147);

    auto g_xxxzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_isg + 148);

    auto g_xxxzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_isg + 149);

    #pragma omp simd aligned(g_xxxz_0_xxxx_0, g_xxxz_0_xxxx_1, g_xxxz_0_xxxy_0, g_xxxz_0_xxxy_1, g_xxxz_0_xxyy_0, g_xxxz_0_xxyy_1, g_xxxz_0_xyyy_0, g_xxxz_0_xyyy_1, g_xxxzz_0_xxxx_1, g_xxxzz_0_xxxy_1, g_xxxzz_0_xxyy_1, g_xxxzz_0_xyyy_1, g_xxxzzz_0_xxxx_0, g_xxxzzz_0_xxxy_0, g_xxxzzz_0_xxxz_0, g_xxxzzz_0_xxyy_0, g_xxxzzz_0_xxyz_0, g_xxxzzz_0_xxzz_0, g_xxxzzz_0_xyyy_0, g_xxxzzz_0_xyyz_0, g_xxxzzz_0_xyzz_0, g_xxxzzz_0_xzzz_0, g_xxxzzz_0_yyyy_0, g_xxxzzz_0_yyyz_0, g_xxxzzz_0_yyzz_0, g_xxxzzz_0_yzzz_0, g_xxxzzz_0_zzzz_0, g_xxzzz_0_xxxz_1, g_xxzzz_0_xxyz_1, g_xxzzz_0_xxz_1, g_xxzzz_0_xxzz_1, g_xxzzz_0_xyyz_1, g_xxzzz_0_xyz_1, g_xxzzz_0_xyzz_1, g_xxzzz_0_xzz_1, g_xxzzz_0_xzzz_1, g_xxzzz_0_yyyy_1, g_xxzzz_0_yyyz_1, g_xxzzz_0_yyz_1, g_xxzzz_0_yyzz_1, g_xxzzz_0_yzz_1, g_xxzzz_0_yzzz_1, g_xxzzz_0_zzz_1, g_xxzzz_0_zzzz_1, g_xzzz_0_xxxz_0, g_xzzz_0_xxxz_1, g_xzzz_0_xxyz_0, g_xzzz_0_xxyz_1, g_xzzz_0_xxzz_0, g_xzzz_0_xxzz_1, g_xzzz_0_xyyz_0, g_xzzz_0_xyyz_1, g_xzzz_0_xyzz_0, g_xzzz_0_xyzz_1, g_xzzz_0_xzzz_0, g_xzzz_0_xzzz_1, g_xzzz_0_yyyy_0, g_xzzz_0_yyyy_1, g_xzzz_0_yyyz_0, g_xzzz_0_yyyz_1, g_xzzz_0_yyzz_0, g_xzzz_0_yyzz_1, g_xzzz_0_yzzz_0, g_xzzz_0_yzzz_1, g_xzzz_0_zzzz_0, g_xzzz_0_zzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxzzz_0_xxxx_0[i] = 2.0 * g_xxxz_0_xxxx_0[i] * fbe_0 - 2.0 * g_xxxz_0_xxxx_1[i] * fz_be_0 + g_xxxzz_0_xxxx_1[i] * wa_z[i];

        g_xxxzzz_0_xxxy_0[i] = 2.0 * g_xxxz_0_xxxy_0[i] * fbe_0 - 2.0 * g_xxxz_0_xxxy_1[i] * fz_be_0 + g_xxxzz_0_xxxy_1[i] * wa_z[i];

        g_xxxzzz_0_xxxz_0[i] = 2.0 * g_xzzz_0_xxxz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxxz_1[i] * fz_be_0 + 3.0 * g_xxzzz_0_xxz_1[i] * fi_acd_0 + g_xxzzz_0_xxxz_1[i] * wa_x[i];

        g_xxxzzz_0_xxyy_0[i] = 2.0 * g_xxxz_0_xxyy_0[i] * fbe_0 - 2.0 * g_xxxz_0_xxyy_1[i] * fz_be_0 + g_xxxzz_0_xxyy_1[i] * wa_z[i];

        g_xxxzzz_0_xxyz_0[i] = 2.0 * g_xzzz_0_xxyz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxyz_1[i] * fz_be_0 + 2.0 * g_xxzzz_0_xyz_1[i] * fi_acd_0 + g_xxzzz_0_xxyz_1[i] * wa_x[i];

        g_xxxzzz_0_xxzz_0[i] = 2.0 * g_xzzz_0_xxzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxzz_1[i] * fz_be_0 + 2.0 * g_xxzzz_0_xzz_1[i] * fi_acd_0 + g_xxzzz_0_xxzz_1[i] * wa_x[i];

        g_xxxzzz_0_xyyy_0[i] = 2.0 * g_xxxz_0_xyyy_0[i] * fbe_0 - 2.0 * g_xxxz_0_xyyy_1[i] * fz_be_0 + g_xxxzz_0_xyyy_1[i] * wa_z[i];

        g_xxxzzz_0_xyyz_0[i] = 2.0 * g_xzzz_0_xyyz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xyyz_1[i] * fz_be_0 + g_xxzzz_0_yyz_1[i] * fi_acd_0 + g_xxzzz_0_xyyz_1[i] * wa_x[i];

        g_xxxzzz_0_xyzz_0[i] = 2.0 * g_xzzz_0_xyzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xyzz_1[i] * fz_be_0 + g_xxzzz_0_yzz_1[i] * fi_acd_0 + g_xxzzz_0_xyzz_1[i] * wa_x[i];

        g_xxxzzz_0_xzzz_0[i] = 2.0 * g_xzzz_0_xzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xzzz_1[i] * fz_be_0 + g_xxzzz_0_zzz_1[i] * fi_acd_0 + g_xxzzz_0_xzzz_1[i] * wa_x[i];

        g_xxxzzz_0_yyyy_0[i] = 2.0 * g_xzzz_0_yyyy_0[i] * fbe_0 - 2.0 * g_xzzz_0_yyyy_1[i] * fz_be_0 + g_xxzzz_0_yyyy_1[i] * wa_x[i];

        g_xxxzzz_0_yyyz_0[i] = 2.0 * g_xzzz_0_yyyz_0[i] * fbe_0 - 2.0 * g_xzzz_0_yyyz_1[i] * fz_be_0 + g_xxzzz_0_yyyz_1[i] * wa_x[i];

        g_xxxzzz_0_yyzz_0[i] = 2.0 * g_xzzz_0_yyzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_yyzz_1[i] * fz_be_0 + g_xxzzz_0_yyzz_1[i] * wa_x[i];

        g_xxxzzz_0_yzzz_0[i] = 2.0 * g_xzzz_0_yzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_yzzz_1[i] * fz_be_0 + g_xxzzz_0_yzzz_1[i] * wa_x[i];

        g_xxxzzz_0_zzzz_0[i] = 2.0 * g_xzzz_0_zzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_zzzz_1[i] * fz_be_0 + g_xxzzz_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 150-165 components of targeted buffer : ISG

    auto g_xxyyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_isg + 150);

    auto g_xxyyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_isg + 151);

    auto g_xxyyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_isg + 152);

    auto g_xxyyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_isg + 153);

    auto g_xxyyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_isg + 154);

    auto g_xxyyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_isg + 155);

    auto g_xxyyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_isg + 156);

    auto g_xxyyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_isg + 157);

    auto g_xxyyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_isg + 158);

    auto g_xxyyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_isg + 159);

    auto g_xxyyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_isg + 160);

    auto g_xxyyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_isg + 161);

    auto g_xxyyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_isg + 162);

    auto g_xxyyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_isg + 163);

    auto g_xxyyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_isg + 164);

    #pragma omp simd aligned(g_xxyy_0_xxxx_0, g_xxyy_0_xxxx_1, g_xxyy_0_xxxz_0, g_xxyy_0_xxxz_1, g_xxyy_0_xxzz_0, g_xxyy_0_xxzz_1, g_xxyy_0_xzzz_0, g_xxyy_0_xzzz_1, g_xxyyy_0_xxxx_1, g_xxyyy_0_xxxz_1, g_xxyyy_0_xxzz_1, g_xxyyy_0_xzzz_1, g_xxyyyy_0_xxxx_0, g_xxyyyy_0_xxxy_0, g_xxyyyy_0_xxxz_0, g_xxyyyy_0_xxyy_0, g_xxyyyy_0_xxyz_0, g_xxyyyy_0_xxzz_0, g_xxyyyy_0_xyyy_0, g_xxyyyy_0_xyyz_0, g_xxyyyy_0_xyzz_0, g_xxyyyy_0_xzzz_0, g_xxyyyy_0_yyyy_0, g_xxyyyy_0_yyyz_0, g_xxyyyy_0_yyzz_0, g_xxyyyy_0_yzzz_0, g_xxyyyy_0_zzzz_0, g_xyyyy_0_xxxy_1, g_xyyyy_0_xxy_1, g_xyyyy_0_xxyy_1, g_xyyyy_0_xxyz_1, g_xyyyy_0_xyy_1, g_xyyyy_0_xyyy_1, g_xyyyy_0_xyyz_1, g_xyyyy_0_xyz_1, g_xyyyy_0_xyzz_1, g_xyyyy_0_yyy_1, g_xyyyy_0_yyyy_1, g_xyyyy_0_yyyz_1, g_xyyyy_0_yyz_1, g_xyyyy_0_yyzz_1, g_xyyyy_0_yzz_1, g_xyyyy_0_yzzz_1, g_xyyyy_0_zzzz_1, g_yyyy_0_xxxy_0, g_yyyy_0_xxxy_1, g_yyyy_0_xxyy_0, g_yyyy_0_xxyy_1, g_yyyy_0_xxyz_0, g_yyyy_0_xxyz_1, g_yyyy_0_xyyy_0, g_yyyy_0_xyyy_1, g_yyyy_0_xyyz_0, g_yyyy_0_xyyz_1, g_yyyy_0_xyzz_0, g_yyyy_0_xyzz_1, g_yyyy_0_yyyy_0, g_yyyy_0_yyyy_1, g_yyyy_0_yyyz_0, g_yyyy_0_yyyz_1, g_yyyy_0_yyzz_0, g_yyyy_0_yyzz_1, g_yyyy_0_yzzz_0, g_yyyy_0_yzzz_1, g_yyyy_0_zzzz_0, g_yyyy_0_zzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyyyy_0_xxxx_0[i] = 3.0 * g_xxyy_0_xxxx_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxx_1[i] * fz_be_0 + g_xxyyy_0_xxxx_1[i] * wa_y[i];

        g_xxyyyy_0_xxxy_0[i] = g_yyyy_0_xxxy_0[i] * fbe_0 - g_yyyy_0_xxxy_1[i] * fz_be_0 + 3.0 * g_xyyyy_0_xxy_1[i] * fi_acd_0 + g_xyyyy_0_xxxy_1[i] * wa_x[i];

        g_xxyyyy_0_xxxz_0[i] = 3.0 * g_xxyy_0_xxxz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxz_1[i] * fz_be_0 + g_xxyyy_0_xxxz_1[i] * wa_y[i];

        g_xxyyyy_0_xxyy_0[i] = g_yyyy_0_xxyy_0[i] * fbe_0 - g_yyyy_0_xxyy_1[i] * fz_be_0 + 2.0 * g_xyyyy_0_xyy_1[i] * fi_acd_0 + g_xyyyy_0_xxyy_1[i] * wa_x[i];

        g_xxyyyy_0_xxyz_0[i] = g_yyyy_0_xxyz_0[i] * fbe_0 - g_yyyy_0_xxyz_1[i] * fz_be_0 + 2.0 * g_xyyyy_0_xyz_1[i] * fi_acd_0 + g_xyyyy_0_xxyz_1[i] * wa_x[i];

        g_xxyyyy_0_xxzz_0[i] = 3.0 * g_xxyy_0_xxzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxzz_1[i] * fz_be_0 + g_xxyyy_0_xxzz_1[i] * wa_y[i];

        g_xxyyyy_0_xyyy_0[i] = g_yyyy_0_xyyy_0[i] * fbe_0 - g_yyyy_0_xyyy_1[i] * fz_be_0 + g_xyyyy_0_yyy_1[i] * fi_acd_0 + g_xyyyy_0_xyyy_1[i] * wa_x[i];

        g_xxyyyy_0_xyyz_0[i] = g_yyyy_0_xyyz_0[i] * fbe_0 - g_yyyy_0_xyyz_1[i] * fz_be_0 + g_xyyyy_0_yyz_1[i] * fi_acd_0 + g_xyyyy_0_xyyz_1[i] * wa_x[i];

        g_xxyyyy_0_xyzz_0[i] = g_yyyy_0_xyzz_0[i] * fbe_0 - g_yyyy_0_xyzz_1[i] * fz_be_0 + g_xyyyy_0_yzz_1[i] * fi_acd_0 + g_xyyyy_0_xyzz_1[i] * wa_x[i];

        g_xxyyyy_0_xzzz_0[i] = 3.0 * g_xxyy_0_xzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xzzz_1[i] * fz_be_0 + g_xxyyy_0_xzzz_1[i] * wa_y[i];

        g_xxyyyy_0_yyyy_0[i] = g_yyyy_0_yyyy_0[i] * fbe_0 - g_yyyy_0_yyyy_1[i] * fz_be_0 + g_xyyyy_0_yyyy_1[i] * wa_x[i];

        g_xxyyyy_0_yyyz_0[i] = g_yyyy_0_yyyz_0[i] * fbe_0 - g_yyyy_0_yyyz_1[i] * fz_be_0 + g_xyyyy_0_yyyz_1[i] * wa_x[i];

        g_xxyyyy_0_yyzz_0[i] = g_yyyy_0_yyzz_0[i] * fbe_0 - g_yyyy_0_yyzz_1[i] * fz_be_0 + g_xyyyy_0_yyzz_1[i] * wa_x[i];

        g_xxyyyy_0_yzzz_0[i] = g_yyyy_0_yzzz_0[i] * fbe_0 - g_yyyy_0_yzzz_1[i] * fz_be_0 + g_xyyyy_0_yzzz_1[i] * wa_x[i];

        g_xxyyyy_0_zzzz_0[i] = g_yyyy_0_zzzz_0[i] * fbe_0 - g_yyyy_0_zzzz_1[i] * fz_be_0 + g_xyyyy_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 165-180 components of targeted buffer : ISG

    auto g_xxyyyz_0_xxxx_0 = pbuffer.data(idx_eri_0_isg + 165);

    auto g_xxyyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_isg + 166);

    auto g_xxyyyz_0_xxxz_0 = pbuffer.data(idx_eri_0_isg + 167);

    auto g_xxyyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_isg + 168);

    auto g_xxyyyz_0_xxyz_0 = pbuffer.data(idx_eri_0_isg + 169);

    auto g_xxyyyz_0_xxzz_0 = pbuffer.data(idx_eri_0_isg + 170);

    auto g_xxyyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_isg + 171);

    auto g_xxyyyz_0_xyyz_0 = pbuffer.data(idx_eri_0_isg + 172);

    auto g_xxyyyz_0_xyzz_0 = pbuffer.data(idx_eri_0_isg + 173);

    auto g_xxyyyz_0_xzzz_0 = pbuffer.data(idx_eri_0_isg + 174);

    auto g_xxyyyz_0_yyyy_0 = pbuffer.data(idx_eri_0_isg + 175);

    auto g_xxyyyz_0_yyyz_0 = pbuffer.data(idx_eri_0_isg + 176);

    auto g_xxyyyz_0_yyzz_0 = pbuffer.data(idx_eri_0_isg + 177);

    auto g_xxyyyz_0_yzzz_0 = pbuffer.data(idx_eri_0_isg + 178);

    auto g_xxyyyz_0_zzzz_0 = pbuffer.data(idx_eri_0_isg + 179);

    #pragma omp simd aligned(g_xxyyy_0_xxx_1, g_xxyyy_0_xxxx_1, g_xxyyy_0_xxxy_1, g_xxyyy_0_xxxz_1, g_xxyyy_0_xxy_1, g_xxyyy_0_xxyy_1, g_xxyyy_0_xxyz_1, g_xxyyy_0_xxz_1, g_xxyyy_0_xxzz_1, g_xxyyy_0_xyy_1, g_xxyyy_0_xyyy_1, g_xxyyy_0_xyyz_1, g_xxyyy_0_xyz_1, g_xxyyy_0_xyzz_1, g_xxyyy_0_xzz_1, g_xxyyy_0_xzzz_1, g_xxyyy_0_yyy_1, g_xxyyy_0_yyyy_1, g_xxyyy_0_yyyz_1, g_xxyyy_0_yyz_1, g_xxyyy_0_yyzz_1, g_xxyyy_0_yzz_1, g_xxyyy_0_yzzz_1, g_xxyyy_0_zzz_1, g_xxyyy_0_zzzz_1, g_xxyyyz_0_xxxx_0, g_xxyyyz_0_xxxy_0, g_xxyyyz_0_xxxz_0, g_xxyyyz_0_xxyy_0, g_xxyyyz_0_xxyz_0, g_xxyyyz_0_xxzz_0, g_xxyyyz_0_xyyy_0, g_xxyyyz_0_xyyz_0, g_xxyyyz_0_xyzz_0, g_xxyyyz_0_xzzz_0, g_xxyyyz_0_yyyy_0, g_xxyyyz_0_yyyz_0, g_xxyyyz_0_yyzz_0, g_xxyyyz_0_yzzz_0, g_xxyyyz_0_zzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyyyz_0_xxxx_0[i] = g_xxyyy_0_xxxx_1[i] * wa_z[i];

        g_xxyyyz_0_xxxy_0[i] = g_xxyyy_0_xxxy_1[i] * wa_z[i];

        g_xxyyyz_0_xxxz_0[i] = g_xxyyy_0_xxx_1[i] * fi_acd_0 + g_xxyyy_0_xxxz_1[i] * wa_z[i];

        g_xxyyyz_0_xxyy_0[i] = g_xxyyy_0_xxyy_1[i] * wa_z[i];

        g_xxyyyz_0_xxyz_0[i] = g_xxyyy_0_xxy_1[i] * fi_acd_0 + g_xxyyy_0_xxyz_1[i] * wa_z[i];

        g_xxyyyz_0_xxzz_0[i] = 2.0 * g_xxyyy_0_xxz_1[i] * fi_acd_0 + g_xxyyy_0_xxzz_1[i] * wa_z[i];

        g_xxyyyz_0_xyyy_0[i] = g_xxyyy_0_xyyy_1[i] * wa_z[i];

        g_xxyyyz_0_xyyz_0[i] = g_xxyyy_0_xyy_1[i] * fi_acd_0 + g_xxyyy_0_xyyz_1[i] * wa_z[i];

        g_xxyyyz_0_xyzz_0[i] = 2.0 * g_xxyyy_0_xyz_1[i] * fi_acd_0 + g_xxyyy_0_xyzz_1[i] * wa_z[i];

        g_xxyyyz_0_xzzz_0[i] = 3.0 * g_xxyyy_0_xzz_1[i] * fi_acd_0 + g_xxyyy_0_xzzz_1[i] * wa_z[i];

        g_xxyyyz_0_yyyy_0[i] = g_xxyyy_0_yyyy_1[i] * wa_z[i];

        g_xxyyyz_0_yyyz_0[i] = g_xxyyy_0_yyy_1[i] * fi_acd_0 + g_xxyyy_0_yyyz_1[i] * wa_z[i];

        g_xxyyyz_0_yyzz_0[i] = 2.0 * g_xxyyy_0_yyz_1[i] * fi_acd_0 + g_xxyyy_0_yyzz_1[i] * wa_z[i];

        g_xxyyyz_0_yzzz_0[i] = 3.0 * g_xxyyy_0_yzz_1[i] * fi_acd_0 + g_xxyyy_0_yzzz_1[i] * wa_z[i];

        g_xxyyyz_0_zzzz_0[i] = 4.0 * g_xxyyy_0_zzz_1[i] * fi_acd_0 + g_xxyyy_0_zzzz_1[i] * wa_z[i];
    }

    /// Set up 180-195 components of targeted buffer : ISG

    auto g_xxyyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_isg + 180);

    auto g_xxyyzz_0_xxxy_0 = pbuffer.data(idx_eri_0_isg + 181);

    auto g_xxyyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_isg + 182);

    auto g_xxyyzz_0_xxyy_0 = pbuffer.data(idx_eri_0_isg + 183);

    auto g_xxyyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_isg + 184);

    auto g_xxyyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_isg + 185);

    auto g_xxyyzz_0_xyyy_0 = pbuffer.data(idx_eri_0_isg + 186);

    auto g_xxyyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_isg + 187);

    auto g_xxyyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_isg + 188);

    auto g_xxyyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_isg + 189);

    auto g_xxyyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_isg + 190);

    auto g_xxyyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_isg + 191);

    auto g_xxyyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_isg + 192);

    auto g_xxyyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_isg + 193);

    auto g_xxyyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_isg + 194);

    #pragma omp simd aligned(g_xxyy_0_xxxy_0, g_xxyy_0_xxxy_1, g_xxyy_0_xxyy_0, g_xxyy_0_xxyy_1, g_xxyy_0_xyyy_0, g_xxyy_0_xyyy_1, g_xxyyz_0_xxxy_1, g_xxyyz_0_xxyy_1, g_xxyyz_0_xyyy_1, g_xxyyzz_0_xxxx_0, g_xxyyzz_0_xxxy_0, g_xxyyzz_0_xxxz_0, g_xxyyzz_0_xxyy_0, g_xxyyzz_0_xxyz_0, g_xxyyzz_0_xxzz_0, g_xxyyzz_0_xyyy_0, g_xxyyzz_0_xyyz_0, g_xxyyzz_0_xyzz_0, g_xxyyzz_0_xzzz_0, g_xxyyzz_0_yyyy_0, g_xxyyzz_0_yyyz_0, g_xxyyzz_0_yyzz_0, g_xxyyzz_0_yzzz_0, g_xxyyzz_0_zzzz_0, g_xxyzz_0_xxxx_1, g_xxyzz_0_xxxz_1, g_xxyzz_0_xxzz_1, g_xxyzz_0_xzzz_1, g_xxzz_0_xxxx_0, g_xxzz_0_xxxx_1, g_xxzz_0_xxxz_0, g_xxzz_0_xxxz_1, g_xxzz_0_xxzz_0, g_xxzz_0_xxzz_1, g_xxzz_0_xzzz_0, g_xxzz_0_xzzz_1, g_xyyzz_0_xxyz_1, g_xyyzz_0_xyyz_1, g_xyyzz_0_xyz_1, g_xyyzz_0_xyzz_1, g_xyyzz_0_yyyy_1, g_xyyzz_0_yyyz_1, g_xyyzz_0_yyz_1, g_xyyzz_0_yyzz_1, g_xyyzz_0_yzz_1, g_xyyzz_0_yzzz_1, g_xyyzz_0_zzzz_1, g_yyzz_0_xxyz_0, g_yyzz_0_xxyz_1, g_yyzz_0_xyyz_0, g_yyzz_0_xyyz_1, g_yyzz_0_xyzz_0, g_yyzz_0_xyzz_1, g_yyzz_0_yyyy_0, g_yyzz_0_yyyy_1, g_yyzz_0_yyyz_0, g_yyzz_0_yyyz_1, g_yyzz_0_yyzz_0, g_yyzz_0_yyzz_1, g_yyzz_0_yzzz_0, g_yyzz_0_yzzz_1, g_yyzz_0_zzzz_0, g_yyzz_0_zzzz_1, wa_x, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyyzz_0_xxxx_0[i] = g_xxzz_0_xxxx_0[i] * fbe_0 - g_xxzz_0_xxxx_1[i] * fz_be_0 + g_xxyzz_0_xxxx_1[i] * wa_y[i];

        g_xxyyzz_0_xxxy_0[i] = g_xxyy_0_xxxy_0[i] * fbe_0 - g_xxyy_0_xxxy_1[i] * fz_be_0 + g_xxyyz_0_xxxy_1[i] * wa_z[i];

        g_xxyyzz_0_xxxz_0[i] = g_xxzz_0_xxxz_0[i] * fbe_0 - g_xxzz_0_xxxz_1[i] * fz_be_0 + g_xxyzz_0_xxxz_1[i] * wa_y[i];

        g_xxyyzz_0_xxyy_0[i] = g_xxyy_0_xxyy_0[i] * fbe_0 - g_xxyy_0_xxyy_1[i] * fz_be_0 + g_xxyyz_0_xxyy_1[i] * wa_z[i];

        g_xxyyzz_0_xxyz_0[i] = g_yyzz_0_xxyz_0[i] * fbe_0 - g_yyzz_0_xxyz_1[i] * fz_be_0 + 2.0 * g_xyyzz_0_xyz_1[i] * fi_acd_0 + g_xyyzz_0_xxyz_1[i] * wa_x[i];

        g_xxyyzz_0_xxzz_0[i] = g_xxzz_0_xxzz_0[i] * fbe_0 - g_xxzz_0_xxzz_1[i] * fz_be_0 + g_xxyzz_0_xxzz_1[i] * wa_y[i];

        g_xxyyzz_0_xyyy_0[i] = g_xxyy_0_xyyy_0[i] * fbe_0 - g_xxyy_0_xyyy_1[i] * fz_be_0 + g_xxyyz_0_xyyy_1[i] * wa_z[i];

        g_xxyyzz_0_xyyz_0[i] = g_yyzz_0_xyyz_0[i] * fbe_0 - g_yyzz_0_xyyz_1[i] * fz_be_0 + g_xyyzz_0_yyz_1[i] * fi_acd_0 + g_xyyzz_0_xyyz_1[i] * wa_x[i];

        g_xxyyzz_0_xyzz_0[i] = g_yyzz_0_xyzz_0[i] * fbe_0 - g_yyzz_0_xyzz_1[i] * fz_be_0 + g_xyyzz_0_yzz_1[i] * fi_acd_0 + g_xyyzz_0_xyzz_1[i] * wa_x[i];

        g_xxyyzz_0_xzzz_0[i] = g_xxzz_0_xzzz_0[i] * fbe_0 - g_xxzz_0_xzzz_1[i] * fz_be_0 + g_xxyzz_0_xzzz_1[i] * wa_y[i];

        g_xxyyzz_0_yyyy_0[i] = g_yyzz_0_yyyy_0[i] * fbe_0 - g_yyzz_0_yyyy_1[i] * fz_be_0 + g_xyyzz_0_yyyy_1[i] * wa_x[i];

        g_xxyyzz_0_yyyz_0[i] = g_yyzz_0_yyyz_0[i] * fbe_0 - g_yyzz_0_yyyz_1[i] * fz_be_0 + g_xyyzz_0_yyyz_1[i] * wa_x[i];

        g_xxyyzz_0_yyzz_0[i] = g_yyzz_0_yyzz_0[i] * fbe_0 - g_yyzz_0_yyzz_1[i] * fz_be_0 + g_xyyzz_0_yyzz_1[i] * wa_x[i];

        g_xxyyzz_0_yzzz_0[i] = g_yyzz_0_yzzz_0[i] * fbe_0 - g_yyzz_0_yzzz_1[i] * fz_be_0 + g_xyyzz_0_yzzz_1[i] * wa_x[i];

        g_xxyyzz_0_zzzz_0[i] = g_yyzz_0_zzzz_0[i] * fbe_0 - g_yyzz_0_zzzz_1[i] * fz_be_0 + g_xyyzz_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 195-210 components of targeted buffer : ISG

    auto g_xxyzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_isg + 195);

    auto g_xxyzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_isg + 196);

    auto g_xxyzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_isg + 197);

    auto g_xxyzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_isg + 198);

    auto g_xxyzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_isg + 199);

    auto g_xxyzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_isg + 200);

    auto g_xxyzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_isg + 201);

    auto g_xxyzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_isg + 202);

    auto g_xxyzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_isg + 203);

    auto g_xxyzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_isg + 204);

    auto g_xxyzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_isg + 205);

    auto g_xxyzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_isg + 206);

    auto g_xxyzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_isg + 207);

    auto g_xxyzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_isg + 208);

    auto g_xxyzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_isg + 209);

    #pragma omp simd aligned(g_xxyzzz_0_xxxx_0, g_xxyzzz_0_xxxy_0, g_xxyzzz_0_xxxz_0, g_xxyzzz_0_xxyy_0, g_xxyzzz_0_xxyz_0, g_xxyzzz_0_xxzz_0, g_xxyzzz_0_xyyy_0, g_xxyzzz_0_xyyz_0, g_xxyzzz_0_xyzz_0, g_xxyzzz_0_xzzz_0, g_xxyzzz_0_yyyy_0, g_xxyzzz_0_yyyz_0, g_xxyzzz_0_yyzz_0, g_xxyzzz_0_yzzz_0, g_xxyzzz_0_zzzz_0, g_xxzzz_0_xxx_1, g_xxzzz_0_xxxx_1, g_xxzzz_0_xxxy_1, g_xxzzz_0_xxxz_1, g_xxzzz_0_xxy_1, g_xxzzz_0_xxyy_1, g_xxzzz_0_xxyz_1, g_xxzzz_0_xxz_1, g_xxzzz_0_xxzz_1, g_xxzzz_0_xyy_1, g_xxzzz_0_xyyy_1, g_xxzzz_0_xyyz_1, g_xxzzz_0_xyz_1, g_xxzzz_0_xyzz_1, g_xxzzz_0_xzz_1, g_xxzzz_0_xzzz_1, g_xxzzz_0_yyy_1, g_xxzzz_0_yyyy_1, g_xxzzz_0_yyyz_1, g_xxzzz_0_yyz_1, g_xxzzz_0_yyzz_1, g_xxzzz_0_yzz_1, g_xxzzz_0_yzzz_1, g_xxzzz_0_zzz_1, g_xxzzz_0_zzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyzzz_0_xxxx_0[i] = g_xxzzz_0_xxxx_1[i] * wa_y[i];

        g_xxyzzz_0_xxxy_0[i] = g_xxzzz_0_xxx_1[i] * fi_acd_0 + g_xxzzz_0_xxxy_1[i] * wa_y[i];

        g_xxyzzz_0_xxxz_0[i] = g_xxzzz_0_xxxz_1[i] * wa_y[i];

        g_xxyzzz_0_xxyy_0[i] = 2.0 * g_xxzzz_0_xxy_1[i] * fi_acd_0 + g_xxzzz_0_xxyy_1[i] * wa_y[i];

        g_xxyzzz_0_xxyz_0[i] = g_xxzzz_0_xxz_1[i] * fi_acd_0 + g_xxzzz_0_xxyz_1[i] * wa_y[i];

        g_xxyzzz_0_xxzz_0[i] = g_xxzzz_0_xxzz_1[i] * wa_y[i];

        g_xxyzzz_0_xyyy_0[i] = 3.0 * g_xxzzz_0_xyy_1[i] * fi_acd_0 + g_xxzzz_0_xyyy_1[i] * wa_y[i];

        g_xxyzzz_0_xyyz_0[i] = 2.0 * g_xxzzz_0_xyz_1[i] * fi_acd_0 + g_xxzzz_0_xyyz_1[i] * wa_y[i];

        g_xxyzzz_0_xyzz_0[i] = g_xxzzz_0_xzz_1[i] * fi_acd_0 + g_xxzzz_0_xyzz_1[i] * wa_y[i];

        g_xxyzzz_0_xzzz_0[i] = g_xxzzz_0_xzzz_1[i] * wa_y[i];

        g_xxyzzz_0_yyyy_0[i] = 4.0 * g_xxzzz_0_yyy_1[i] * fi_acd_0 + g_xxzzz_0_yyyy_1[i] * wa_y[i];

        g_xxyzzz_0_yyyz_0[i] = 3.0 * g_xxzzz_0_yyz_1[i] * fi_acd_0 + g_xxzzz_0_yyyz_1[i] * wa_y[i];

        g_xxyzzz_0_yyzz_0[i] = 2.0 * g_xxzzz_0_yzz_1[i] * fi_acd_0 + g_xxzzz_0_yyzz_1[i] * wa_y[i];

        g_xxyzzz_0_yzzz_0[i] = g_xxzzz_0_zzz_1[i] * fi_acd_0 + g_xxzzz_0_yzzz_1[i] * wa_y[i];

        g_xxyzzz_0_zzzz_0[i] = g_xxzzz_0_zzzz_1[i] * wa_y[i];
    }

    /// Set up 210-225 components of targeted buffer : ISG

    auto g_xxzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_isg + 210);

    auto g_xxzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_isg + 211);

    auto g_xxzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_isg + 212);

    auto g_xxzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_isg + 213);

    auto g_xxzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_isg + 214);

    auto g_xxzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_isg + 215);

    auto g_xxzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_isg + 216);

    auto g_xxzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_isg + 217);

    auto g_xxzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_isg + 218);

    auto g_xxzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_isg + 219);

    auto g_xxzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_isg + 220);

    auto g_xxzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_isg + 221);

    auto g_xxzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_isg + 222);

    auto g_xxzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_isg + 223);

    auto g_xxzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_isg + 224);

    #pragma omp simd aligned(g_xxzz_0_xxxx_0, g_xxzz_0_xxxx_1, g_xxzz_0_xxxy_0, g_xxzz_0_xxxy_1, g_xxzz_0_xxyy_0, g_xxzz_0_xxyy_1, g_xxzz_0_xyyy_0, g_xxzz_0_xyyy_1, g_xxzzz_0_xxxx_1, g_xxzzz_0_xxxy_1, g_xxzzz_0_xxyy_1, g_xxzzz_0_xyyy_1, g_xxzzzz_0_xxxx_0, g_xxzzzz_0_xxxy_0, g_xxzzzz_0_xxxz_0, g_xxzzzz_0_xxyy_0, g_xxzzzz_0_xxyz_0, g_xxzzzz_0_xxzz_0, g_xxzzzz_0_xyyy_0, g_xxzzzz_0_xyyz_0, g_xxzzzz_0_xyzz_0, g_xxzzzz_0_xzzz_0, g_xxzzzz_0_yyyy_0, g_xxzzzz_0_yyyz_0, g_xxzzzz_0_yyzz_0, g_xxzzzz_0_yzzz_0, g_xxzzzz_0_zzzz_0, g_xzzzz_0_xxxz_1, g_xzzzz_0_xxyz_1, g_xzzzz_0_xxz_1, g_xzzzz_0_xxzz_1, g_xzzzz_0_xyyz_1, g_xzzzz_0_xyz_1, g_xzzzz_0_xyzz_1, g_xzzzz_0_xzz_1, g_xzzzz_0_xzzz_1, g_xzzzz_0_yyyy_1, g_xzzzz_0_yyyz_1, g_xzzzz_0_yyz_1, g_xzzzz_0_yyzz_1, g_xzzzz_0_yzz_1, g_xzzzz_0_yzzz_1, g_xzzzz_0_zzz_1, g_xzzzz_0_zzzz_1, g_zzzz_0_xxxz_0, g_zzzz_0_xxxz_1, g_zzzz_0_xxyz_0, g_zzzz_0_xxyz_1, g_zzzz_0_xxzz_0, g_zzzz_0_xxzz_1, g_zzzz_0_xyyz_0, g_zzzz_0_xyyz_1, g_zzzz_0_xyzz_0, g_zzzz_0_xyzz_1, g_zzzz_0_xzzz_0, g_zzzz_0_xzzz_1, g_zzzz_0_yyyy_0, g_zzzz_0_yyyy_1, g_zzzz_0_yyyz_0, g_zzzz_0_yyyz_1, g_zzzz_0_yyzz_0, g_zzzz_0_yyzz_1, g_zzzz_0_yzzz_0, g_zzzz_0_yzzz_1, g_zzzz_0_zzzz_0, g_zzzz_0_zzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxzzzz_0_xxxx_0[i] = 3.0 * g_xxzz_0_xxxx_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxx_1[i] * fz_be_0 + g_xxzzz_0_xxxx_1[i] * wa_z[i];

        g_xxzzzz_0_xxxy_0[i] = 3.0 * g_xxzz_0_xxxy_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxy_1[i] * fz_be_0 + g_xxzzz_0_xxxy_1[i] * wa_z[i];

        g_xxzzzz_0_xxxz_0[i] = g_zzzz_0_xxxz_0[i] * fbe_0 - g_zzzz_0_xxxz_1[i] * fz_be_0 + 3.0 * g_xzzzz_0_xxz_1[i] * fi_acd_0 + g_xzzzz_0_xxxz_1[i] * wa_x[i];

        g_xxzzzz_0_xxyy_0[i] = 3.0 * g_xxzz_0_xxyy_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxyy_1[i] * fz_be_0 + g_xxzzz_0_xxyy_1[i] * wa_z[i];

        g_xxzzzz_0_xxyz_0[i] = g_zzzz_0_xxyz_0[i] * fbe_0 - g_zzzz_0_xxyz_1[i] * fz_be_0 + 2.0 * g_xzzzz_0_xyz_1[i] * fi_acd_0 + g_xzzzz_0_xxyz_1[i] * wa_x[i];

        g_xxzzzz_0_xxzz_0[i] = g_zzzz_0_xxzz_0[i] * fbe_0 - g_zzzz_0_xxzz_1[i] * fz_be_0 + 2.0 * g_xzzzz_0_xzz_1[i] * fi_acd_0 + g_xzzzz_0_xxzz_1[i] * wa_x[i];

        g_xxzzzz_0_xyyy_0[i] = 3.0 * g_xxzz_0_xyyy_0[i] * fbe_0 - 3.0 * g_xxzz_0_xyyy_1[i] * fz_be_0 + g_xxzzz_0_xyyy_1[i] * wa_z[i];

        g_xxzzzz_0_xyyz_0[i] = g_zzzz_0_xyyz_0[i] * fbe_0 - g_zzzz_0_xyyz_1[i] * fz_be_0 + g_xzzzz_0_yyz_1[i] * fi_acd_0 + g_xzzzz_0_xyyz_1[i] * wa_x[i];

        g_xxzzzz_0_xyzz_0[i] = g_zzzz_0_xyzz_0[i] * fbe_0 - g_zzzz_0_xyzz_1[i] * fz_be_0 + g_xzzzz_0_yzz_1[i] * fi_acd_0 + g_xzzzz_0_xyzz_1[i] * wa_x[i];

        g_xxzzzz_0_xzzz_0[i] = g_zzzz_0_xzzz_0[i] * fbe_0 - g_zzzz_0_xzzz_1[i] * fz_be_0 + g_xzzzz_0_zzz_1[i] * fi_acd_0 + g_xzzzz_0_xzzz_1[i] * wa_x[i];

        g_xxzzzz_0_yyyy_0[i] = g_zzzz_0_yyyy_0[i] * fbe_0 - g_zzzz_0_yyyy_1[i] * fz_be_0 + g_xzzzz_0_yyyy_1[i] * wa_x[i];

        g_xxzzzz_0_yyyz_0[i] = g_zzzz_0_yyyz_0[i] * fbe_0 - g_zzzz_0_yyyz_1[i] * fz_be_0 + g_xzzzz_0_yyyz_1[i] * wa_x[i];

        g_xxzzzz_0_yyzz_0[i] = g_zzzz_0_yyzz_0[i] * fbe_0 - g_zzzz_0_yyzz_1[i] * fz_be_0 + g_xzzzz_0_yyzz_1[i] * wa_x[i];

        g_xxzzzz_0_yzzz_0[i] = g_zzzz_0_yzzz_0[i] * fbe_0 - g_zzzz_0_yzzz_1[i] * fz_be_0 + g_xzzzz_0_yzzz_1[i] * wa_x[i];

        g_xxzzzz_0_zzzz_0[i] = g_zzzz_0_zzzz_0[i] * fbe_0 - g_zzzz_0_zzzz_1[i] * fz_be_0 + g_xzzzz_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 225-240 components of targeted buffer : ISG

    auto g_xyyyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_isg + 225);

    auto g_xyyyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_isg + 226);

    auto g_xyyyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_isg + 227);

    auto g_xyyyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_isg + 228);

    auto g_xyyyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_isg + 229);

    auto g_xyyyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_isg + 230);

    auto g_xyyyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_isg + 231);

    auto g_xyyyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_isg + 232);

    auto g_xyyyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_isg + 233);

    auto g_xyyyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_isg + 234);

    auto g_xyyyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_isg + 235);

    auto g_xyyyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_isg + 236);

    auto g_xyyyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_isg + 237);

    auto g_xyyyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_isg + 238);

    auto g_xyyyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_isg + 239);

    #pragma omp simd aligned(g_xyyyyy_0_xxxx_0, g_xyyyyy_0_xxxy_0, g_xyyyyy_0_xxxz_0, g_xyyyyy_0_xxyy_0, g_xyyyyy_0_xxyz_0, g_xyyyyy_0_xxzz_0, g_xyyyyy_0_xyyy_0, g_xyyyyy_0_xyyz_0, g_xyyyyy_0_xyzz_0, g_xyyyyy_0_xzzz_0, g_xyyyyy_0_yyyy_0, g_xyyyyy_0_yyyz_0, g_xyyyyy_0_yyzz_0, g_xyyyyy_0_yzzz_0, g_xyyyyy_0_zzzz_0, g_yyyyy_0_xxx_1, g_yyyyy_0_xxxx_1, g_yyyyy_0_xxxy_1, g_yyyyy_0_xxxz_1, g_yyyyy_0_xxy_1, g_yyyyy_0_xxyy_1, g_yyyyy_0_xxyz_1, g_yyyyy_0_xxz_1, g_yyyyy_0_xxzz_1, g_yyyyy_0_xyy_1, g_yyyyy_0_xyyy_1, g_yyyyy_0_xyyz_1, g_yyyyy_0_xyz_1, g_yyyyy_0_xyzz_1, g_yyyyy_0_xzz_1, g_yyyyy_0_xzzz_1, g_yyyyy_0_yyy_1, g_yyyyy_0_yyyy_1, g_yyyyy_0_yyyz_1, g_yyyyy_0_yyz_1, g_yyyyy_0_yyzz_1, g_yyyyy_0_yzz_1, g_yyyyy_0_yzzz_1, g_yyyyy_0_zzz_1, g_yyyyy_0_zzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyyy_0_xxxx_0[i] = 4.0 * g_yyyyy_0_xxx_1[i] * fi_acd_0 + g_yyyyy_0_xxxx_1[i] * wa_x[i];

        g_xyyyyy_0_xxxy_0[i] = 3.0 * g_yyyyy_0_xxy_1[i] * fi_acd_0 + g_yyyyy_0_xxxy_1[i] * wa_x[i];

        g_xyyyyy_0_xxxz_0[i] = 3.0 * g_yyyyy_0_xxz_1[i] * fi_acd_0 + g_yyyyy_0_xxxz_1[i] * wa_x[i];

        g_xyyyyy_0_xxyy_0[i] = 2.0 * g_yyyyy_0_xyy_1[i] * fi_acd_0 + g_yyyyy_0_xxyy_1[i] * wa_x[i];

        g_xyyyyy_0_xxyz_0[i] = 2.0 * g_yyyyy_0_xyz_1[i] * fi_acd_0 + g_yyyyy_0_xxyz_1[i] * wa_x[i];

        g_xyyyyy_0_xxzz_0[i] = 2.0 * g_yyyyy_0_xzz_1[i] * fi_acd_0 + g_yyyyy_0_xxzz_1[i] * wa_x[i];

        g_xyyyyy_0_xyyy_0[i] = g_yyyyy_0_yyy_1[i] * fi_acd_0 + g_yyyyy_0_xyyy_1[i] * wa_x[i];

        g_xyyyyy_0_xyyz_0[i] = g_yyyyy_0_yyz_1[i] * fi_acd_0 + g_yyyyy_0_xyyz_1[i] * wa_x[i];

        g_xyyyyy_0_xyzz_0[i] = g_yyyyy_0_yzz_1[i] * fi_acd_0 + g_yyyyy_0_xyzz_1[i] * wa_x[i];

        g_xyyyyy_0_xzzz_0[i] = g_yyyyy_0_zzz_1[i] * fi_acd_0 + g_yyyyy_0_xzzz_1[i] * wa_x[i];

        g_xyyyyy_0_yyyy_0[i] = g_yyyyy_0_yyyy_1[i] * wa_x[i];

        g_xyyyyy_0_yyyz_0[i] = g_yyyyy_0_yyyz_1[i] * wa_x[i];

        g_xyyyyy_0_yyzz_0[i] = g_yyyyy_0_yyzz_1[i] * wa_x[i];

        g_xyyyyy_0_yzzz_0[i] = g_yyyyy_0_yzzz_1[i] * wa_x[i];

        g_xyyyyy_0_zzzz_0[i] = g_yyyyy_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 240-255 components of targeted buffer : ISG

    auto g_xyyyyz_0_xxxx_0 = pbuffer.data(idx_eri_0_isg + 240);

    auto g_xyyyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_isg + 241);

    auto g_xyyyyz_0_xxxz_0 = pbuffer.data(idx_eri_0_isg + 242);

    auto g_xyyyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_isg + 243);

    auto g_xyyyyz_0_xxyz_0 = pbuffer.data(idx_eri_0_isg + 244);

    auto g_xyyyyz_0_xxzz_0 = pbuffer.data(idx_eri_0_isg + 245);

    auto g_xyyyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_isg + 246);

    auto g_xyyyyz_0_xyyz_0 = pbuffer.data(idx_eri_0_isg + 247);

    auto g_xyyyyz_0_xyzz_0 = pbuffer.data(idx_eri_0_isg + 248);

    auto g_xyyyyz_0_xzzz_0 = pbuffer.data(idx_eri_0_isg + 249);

    auto g_xyyyyz_0_yyyy_0 = pbuffer.data(idx_eri_0_isg + 250);

    auto g_xyyyyz_0_yyyz_0 = pbuffer.data(idx_eri_0_isg + 251);

    auto g_xyyyyz_0_yyzz_0 = pbuffer.data(idx_eri_0_isg + 252);

    auto g_xyyyyz_0_yzzz_0 = pbuffer.data(idx_eri_0_isg + 253);

    auto g_xyyyyz_0_zzzz_0 = pbuffer.data(idx_eri_0_isg + 254);

    #pragma omp simd aligned(g_xyyyy_0_xxxx_1, g_xyyyy_0_xxxy_1, g_xyyyy_0_xxyy_1, g_xyyyy_0_xyyy_1, g_xyyyyz_0_xxxx_0, g_xyyyyz_0_xxxy_0, g_xyyyyz_0_xxxz_0, g_xyyyyz_0_xxyy_0, g_xyyyyz_0_xxyz_0, g_xyyyyz_0_xxzz_0, g_xyyyyz_0_xyyy_0, g_xyyyyz_0_xyyz_0, g_xyyyyz_0_xyzz_0, g_xyyyyz_0_xzzz_0, g_xyyyyz_0_yyyy_0, g_xyyyyz_0_yyyz_0, g_xyyyyz_0_yyzz_0, g_xyyyyz_0_yzzz_0, g_xyyyyz_0_zzzz_0, g_yyyyz_0_xxxz_1, g_yyyyz_0_xxyz_1, g_yyyyz_0_xxz_1, g_yyyyz_0_xxzz_1, g_yyyyz_0_xyyz_1, g_yyyyz_0_xyz_1, g_yyyyz_0_xyzz_1, g_yyyyz_0_xzz_1, g_yyyyz_0_xzzz_1, g_yyyyz_0_yyyy_1, g_yyyyz_0_yyyz_1, g_yyyyz_0_yyz_1, g_yyyyz_0_yyzz_1, g_yyyyz_0_yzz_1, g_yyyyz_0_yzzz_1, g_yyyyz_0_zzz_1, g_yyyyz_0_zzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyyz_0_xxxx_0[i] = g_xyyyy_0_xxxx_1[i] * wa_z[i];

        g_xyyyyz_0_xxxy_0[i] = g_xyyyy_0_xxxy_1[i] * wa_z[i];

        g_xyyyyz_0_xxxz_0[i] = 3.0 * g_yyyyz_0_xxz_1[i] * fi_acd_0 + g_yyyyz_0_xxxz_1[i] * wa_x[i];

        g_xyyyyz_0_xxyy_0[i] = g_xyyyy_0_xxyy_1[i] * wa_z[i];

        g_xyyyyz_0_xxyz_0[i] = 2.0 * g_yyyyz_0_xyz_1[i] * fi_acd_0 + g_yyyyz_0_xxyz_1[i] * wa_x[i];

        g_xyyyyz_0_xxzz_0[i] = 2.0 * g_yyyyz_0_xzz_1[i] * fi_acd_0 + g_yyyyz_0_xxzz_1[i] * wa_x[i];

        g_xyyyyz_0_xyyy_0[i] = g_xyyyy_0_xyyy_1[i] * wa_z[i];

        g_xyyyyz_0_xyyz_0[i] = g_yyyyz_0_yyz_1[i] * fi_acd_0 + g_yyyyz_0_xyyz_1[i] * wa_x[i];

        g_xyyyyz_0_xyzz_0[i] = g_yyyyz_0_yzz_1[i] * fi_acd_0 + g_yyyyz_0_xyzz_1[i] * wa_x[i];

        g_xyyyyz_0_xzzz_0[i] = g_yyyyz_0_zzz_1[i] * fi_acd_0 + g_yyyyz_0_xzzz_1[i] * wa_x[i];

        g_xyyyyz_0_yyyy_0[i] = g_yyyyz_0_yyyy_1[i] * wa_x[i];

        g_xyyyyz_0_yyyz_0[i] = g_yyyyz_0_yyyz_1[i] * wa_x[i];

        g_xyyyyz_0_yyzz_0[i] = g_yyyyz_0_yyzz_1[i] * wa_x[i];

        g_xyyyyz_0_yzzz_0[i] = g_yyyyz_0_yzzz_1[i] * wa_x[i];

        g_xyyyyz_0_zzzz_0[i] = g_yyyyz_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 255-270 components of targeted buffer : ISG

    auto g_xyyyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_isg + 255);

    auto g_xyyyzz_0_xxxy_0 = pbuffer.data(idx_eri_0_isg + 256);

    auto g_xyyyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_isg + 257);

    auto g_xyyyzz_0_xxyy_0 = pbuffer.data(idx_eri_0_isg + 258);

    auto g_xyyyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_isg + 259);

    auto g_xyyyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_isg + 260);

    auto g_xyyyzz_0_xyyy_0 = pbuffer.data(idx_eri_0_isg + 261);

    auto g_xyyyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_isg + 262);

    auto g_xyyyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_isg + 263);

    auto g_xyyyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_isg + 264);

    auto g_xyyyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_isg + 265);

    auto g_xyyyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_isg + 266);

    auto g_xyyyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_isg + 267);

    auto g_xyyyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_isg + 268);

    auto g_xyyyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_isg + 269);

    #pragma omp simd aligned(g_xyyyzz_0_xxxx_0, g_xyyyzz_0_xxxy_0, g_xyyyzz_0_xxxz_0, g_xyyyzz_0_xxyy_0, g_xyyyzz_0_xxyz_0, g_xyyyzz_0_xxzz_0, g_xyyyzz_0_xyyy_0, g_xyyyzz_0_xyyz_0, g_xyyyzz_0_xyzz_0, g_xyyyzz_0_xzzz_0, g_xyyyzz_0_yyyy_0, g_xyyyzz_0_yyyz_0, g_xyyyzz_0_yyzz_0, g_xyyyzz_0_yzzz_0, g_xyyyzz_0_zzzz_0, g_yyyzz_0_xxx_1, g_yyyzz_0_xxxx_1, g_yyyzz_0_xxxy_1, g_yyyzz_0_xxxz_1, g_yyyzz_0_xxy_1, g_yyyzz_0_xxyy_1, g_yyyzz_0_xxyz_1, g_yyyzz_0_xxz_1, g_yyyzz_0_xxzz_1, g_yyyzz_0_xyy_1, g_yyyzz_0_xyyy_1, g_yyyzz_0_xyyz_1, g_yyyzz_0_xyz_1, g_yyyzz_0_xyzz_1, g_yyyzz_0_xzz_1, g_yyyzz_0_xzzz_1, g_yyyzz_0_yyy_1, g_yyyzz_0_yyyy_1, g_yyyzz_0_yyyz_1, g_yyyzz_0_yyz_1, g_yyyzz_0_yyzz_1, g_yyyzz_0_yzz_1, g_yyyzz_0_yzzz_1, g_yyyzz_0_zzz_1, g_yyyzz_0_zzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyzz_0_xxxx_0[i] = 4.0 * g_yyyzz_0_xxx_1[i] * fi_acd_0 + g_yyyzz_0_xxxx_1[i] * wa_x[i];

        g_xyyyzz_0_xxxy_0[i] = 3.0 * g_yyyzz_0_xxy_1[i] * fi_acd_0 + g_yyyzz_0_xxxy_1[i] * wa_x[i];

        g_xyyyzz_0_xxxz_0[i] = 3.0 * g_yyyzz_0_xxz_1[i] * fi_acd_0 + g_yyyzz_0_xxxz_1[i] * wa_x[i];

        g_xyyyzz_0_xxyy_0[i] = 2.0 * g_yyyzz_0_xyy_1[i] * fi_acd_0 + g_yyyzz_0_xxyy_1[i] * wa_x[i];

        g_xyyyzz_0_xxyz_0[i] = 2.0 * g_yyyzz_0_xyz_1[i] * fi_acd_0 + g_yyyzz_0_xxyz_1[i] * wa_x[i];

        g_xyyyzz_0_xxzz_0[i] = 2.0 * g_yyyzz_0_xzz_1[i] * fi_acd_0 + g_yyyzz_0_xxzz_1[i] * wa_x[i];

        g_xyyyzz_0_xyyy_0[i] = g_yyyzz_0_yyy_1[i] * fi_acd_0 + g_yyyzz_0_xyyy_1[i] * wa_x[i];

        g_xyyyzz_0_xyyz_0[i] = g_yyyzz_0_yyz_1[i] * fi_acd_0 + g_yyyzz_0_xyyz_1[i] * wa_x[i];

        g_xyyyzz_0_xyzz_0[i] = g_yyyzz_0_yzz_1[i] * fi_acd_0 + g_yyyzz_0_xyzz_1[i] * wa_x[i];

        g_xyyyzz_0_xzzz_0[i] = g_yyyzz_0_zzz_1[i] * fi_acd_0 + g_yyyzz_0_xzzz_1[i] * wa_x[i];

        g_xyyyzz_0_yyyy_0[i] = g_yyyzz_0_yyyy_1[i] * wa_x[i];

        g_xyyyzz_0_yyyz_0[i] = g_yyyzz_0_yyyz_1[i] * wa_x[i];

        g_xyyyzz_0_yyzz_0[i] = g_yyyzz_0_yyzz_1[i] * wa_x[i];

        g_xyyyzz_0_yzzz_0[i] = g_yyyzz_0_yzzz_1[i] * wa_x[i];

        g_xyyyzz_0_zzzz_0[i] = g_yyyzz_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 270-285 components of targeted buffer : ISG

    auto g_xyyzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_isg + 270);

    auto g_xyyzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_isg + 271);

    auto g_xyyzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_isg + 272);

    auto g_xyyzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_isg + 273);

    auto g_xyyzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_isg + 274);

    auto g_xyyzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_isg + 275);

    auto g_xyyzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_isg + 276);

    auto g_xyyzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_isg + 277);

    auto g_xyyzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_isg + 278);

    auto g_xyyzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_isg + 279);

    auto g_xyyzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_isg + 280);

    auto g_xyyzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_isg + 281);

    auto g_xyyzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_isg + 282);

    auto g_xyyzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_isg + 283);

    auto g_xyyzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_isg + 284);

    #pragma omp simd aligned(g_xyyzzz_0_xxxx_0, g_xyyzzz_0_xxxy_0, g_xyyzzz_0_xxxz_0, g_xyyzzz_0_xxyy_0, g_xyyzzz_0_xxyz_0, g_xyyzzz_0_xxzz_0, g_xyyzzz_0_xyyy_0, g_xyyzzz_0_xyyz_0, g_xyyzzz_0_xyzz_0, g_xyyzzz_0_xzzz_0, g_xyyzzz_0_yyyy_0, g_xyyzzz_0_yyyz_0, g_xyyzzz_0_yyzz_0, g_xyyzzz_0_yzzz_0, g_xyyzzz_0_zzzz_0, g_yyzzz_0_xxx_1, g_yyzzz_0_xxxx_1, g_yyzzz_0_xxxy_1, g_yyzzz_0_xxxz_1, g_yyzzz_0_xxy_1, g_yyzzz_0_xxyy_1, g_yyzzz_0_xxyz_1, g_yyzzz_0_xxz_1, g_yyzzz_0_xxzz_1, g_yyzzz_0_xyy_1, g_yyzzz_0_xyyy_1, g_yyzzz_0_xyyz_1, g_yyzzz_0_xyz_1, g_yyzzz_0_xyzz_1, g_yyzzz_0_xzz_1, g_yyzzz_0_xzzz_1, g_yyzzz_0_yyy_1, g_yyzzz_0_yyyy_1, g_yyzzz_0_yyyz_1, g_yyzzz_0_yyz_1, g_yyzzz_0_yyzz_1, g_yyzzz_0_yzz_1, g_yyzzz_0_yzzz_1, g_yyzzz_0_zzz_1, g_yyzzz_0_zzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyzzz_0_xxxx_0[i] = 4.0 * g_yyzzz_0_xxx_1[i] * fi_acd_0 + g_yyzzz_0_xxxx_1[i] * wa_x[i];

        g_xyyzzz_0_xxxy_0[i] = 3.0 * g_yyzzz_0_xxy_1[i] * fi_acd_0 + g_yyzzz_0_xxxy_1[i] * wa_x[i];

        g_xyyzzz_0_xxxz_0[i] = 3.0 * g_yyzzz_0_xxz_1[i] * fi_acd_0 + g_yyzzz_0_xxxz_1[i] * wa_x[i];

        g_xyyzzz_0_xxyy_0[i] = 2.0 * g_yyzzz_0_xyy_1[i] * fi_acd_0 + g_yyzzz_0_xxyy_1[i] * wa_x[i];

        g_xyyzzz_0_xxyz_0[i] = 2.0 * g_yyzzz_0_xyz_1[i] * fi_acd_0 + g_yyzzz_0_xxyz_1[i] * wa_x[i];

        g_xyyzzz_0_xxzz_0[i] = 2.0 * g_yyzzz_0_xzz_1[i] * fi_acd_0 + g_yyzzz_0_xxzz_1[i] * wa_x[i];

        g_xyyzzz_0_xyyy_0[i] = g_yyzzz_0_yyy_1[i] * fi_acd_0 + g_yyzzz_0_xyyy_1[i] * wa_x[i];

        g_xyyzzz_0_xyyz_0[i] = g_yyzzz_0_yyz_1[i] * fi_acd_0 + g_yyzzz_0_xyyz_1[i] * wa_x[i];

        g_xyyzzz_0_xyzz_0[i] = g_yyzzz_0_yzz_1[i] * fi_acd_0 + g_yyzzz_0_xyzz_1[i] * wa_x[i];

        g_xyyzzz_0_xzzz_0[i] = g_yyzzz_0_zzz_1[i] * fi_acd_0 + g_yyzzz_0_xzzz_1[i] * wa_x[i];

        g_xyyzzz_0_yyyy_0[i] = g_yyzzz_0_yyyy_1[i] * wa_x[i];

        g_xyyzzz_0_yyyz_0[i] = g_yyzzz_0_yyyz_1[i] * wa_x[i];

        g_xyyzzz_0_yyzz_0[i] = g_yyzzz_0_yyzz_1[i] * wa_x[i];

        g_xyyzzz_0_yzzz_0[i] = g_yyzzz_0_yzzz_1[i] * wa_x[i];

        g_xyyzzz_0_zzzz_0[i] = g_yyzzz_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 285-300 components of targeted buffer : ISG

    auto g_xyzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_isg + 285);

    auto g_xyzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_isg + 286);

    auto g_xyzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_isg + 287);

    auto g_xyzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_isg + 288);

    auto g_xyzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_isg + 289);

    auto g_xyzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_isg + 290);

    auto g_xyzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_isg + 291);

    auto g_xyzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_isg + 292);

    auto g_xyzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_isg + 293);

    auto g_xyzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_isg + 294);

    auto g_xyzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_isg + 295);

    auto g_xyzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_isg + 296);

    auto g_xyzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_isg + 297);

    auto g_xyzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_isg + 298);

    auto g_xyzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_isg + 299);

    #pragma omp simd aligned(g_xyzzzz_0_xxxx_0, g_xyzzzz_0_xxxy_0, g_xyzzzz_0_xxxz_0, g_xyzzzz_0_xxyy_0, g_xyzzzz_0_xxyz_0, g_xyzzzz_0_xxzz_0, g_xyzzzz_0_xyyy_0, g_xyzzzz_0_xyyz_0, g_xyzzzz_0_xyzz_0, g_xyzzzz_0_xzzz_0, g_xyzzzz_0_yyyy_0, g_xyzzzz_0_yyyz_0, g_xyzzzz_0_yyzz_0, g_xyzzzz_0_yzzz_0, g_xyzzzz_0_zzzz_0, g_xzzzz_0_xxxx_1, g_xzzzz_0_xxxz_1, g_xzzzz_0_xxzz_1, g_xzzzz_0_xzzz_1, g_yzzzz_0_xxxy_1, g_yzzzz_0_xxy_1, g_yzzzz_0_xxyy_1, g_yzzzz_0_xxyz_1, g_yzzzz_0_xyy_1, g_yzzzz_0_xyyy_1, g_yzzzz_0_xyyz_1, g_yzzzz_0_xyz_1, g_yzzzz_0_xyzz_1, g_yzzzz_0_yyy_1, g_yzzzz_0_yyyy_1, g_yzzzz_0_yyyz_1, g_yzzzz_0_yyz_1, g_yzzzz_0_yyzz_1, g_yzzzz_0_yzz_1, g_yzzzz_0_yzzz_1, g_yzzzz_0_zzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyzzzz_0_xxxx_0[i] = g_xzzzz_0_xxxx_1[i] * wa_y[i];

        g_xyzzzz_0_xxxy_0[i] = 3.0 * g_yzzzz_0_xxy_1[i] * fi_acd_0 + g_yzzzz_0_xxxy_1[i] * wa_x[i];

        g_xyzzzz_0_xxxz_0[i] = g_xzzzz_0_xxxz_1[i] * wa_y[i];

        g_xyzzzz_0_xxyy_0[i] = 2.0 * g_yzzzz_0_xyy_1[i] * fi_acd_0 + g_yzzzz_0_xxyy_1[i] * wa_x[i];

        g_xyzzzz_0_xxyz_0[i] = 2.0 * g_yzzzz_0_xyz_1[i] * fi_acd_0 + g_yzzzz_0_xxyz_1[i] * wa_x[i];

        g_xyzzzz_0_xxzz_0[i] = g_xzzzz_0_xxzz_1[i] * wa_y[i];

        g_xyzzzz_0_xyyy_0[i] = g_yzzzz_0_yyy_1[i] * fi_acd_0 + g_yzzzz_0_xyyy_1[i] * wa_x[i];

        g_xyzzzz_0_xyyz_0[i] = g_yzzzz_0_yyz_1[i] * fi_acd_0 + g_yzzzz_0_xyyz_1[i] * wa_x[i];

        g_xyzzzz_0_xyzz_0[i] = g_yzzzz_0_yzz_1[i] * fi_acd_0 + g_yzzzz_0_xyzz_1[i] * wa_x[i];

        g_xyzzzz_0_xzzz_0[i] = g_xzzzz_0_xzzz_1[i] * wa_y[i];

        g_xyzzzz_0_yyyy_0[i] = g_yzzzz_0_yyyy_1[i] * wa_x[i];

        g_xyzzzz_0_yyyz_0[i] = g_yzzzz_0_yyyz_1[i] * wa_x[i];

        g_xyzzzz_0_yyzz_0[i] = g_yzzzz_0_yyzz_1[i] * wa_x[i];

        g_xyzzzz_0_yzzz_0[i] = g_yzzzz_0_yzzz_1[i] * wa_x[i];

        g_xyzzzz_0_zzzz_0[i] = g_yzzzz_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 300-315 components of targeted buffer : ISG

    auto g_xzzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_isg + 300);

    auto g_xzzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_isg + 301);

    auto g_xzzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_isg + 302);

    auto g_xzzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_isg + 303);

    auto g_xzzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_isg + 304);

    auto g_xzzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_isg + 305);

    auto g_xzzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_isg + 306);

    auto g_xzzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_isg + 307);

    auto g_xzzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_isg + 308);

    auto g_xzzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_isg + 309);

    auto g_xzzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_isg + 310);

    auto g_xzzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_isg + 311);

    auto g_xzzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_isg + 312);

    auto g_xzzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_isg + 313);

    auto g_xzzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_isg + 314);

    #pragma omp simd aligned(g_xzzzzz_0_xxxx_0, g_xzzzzz_0_xxxy_0, g_xzzzzz_0_xxxz_0, g_xzzzzz_0_xxyy_0, g_xzzzzz_0_xxyz_0, g_xzzzzz_0_xxzz_0, g_xzzzzz_0_xyyy_0, g_xzzzzz_0_xyyz_0, g_xzzzzz_0_xyzz_0, g_xzzzzz_0_xzzz_0, g_xzzzzz_0_yyyy_0, g_xzzzzz_0_yyyz_0, g_xzzzzz_0_yyzz_0, g_xzzzzz_0_yzzz_0, g_xzzzzz_0_zzzz_0, g_zzzzz_0_xxx_1, g_zzzzz_0_xxxx_1, g_zzzzz_0_xxxy_1, g_zzzzz_0_xxxz_1, g_zzzzz_0_xxy_1, g_zzzzz_0_xxyy_1, g_zzzzz_0_xxyz_1, g_zzzzz_0_xxz_1, g_zzzzz_0_xxzz_1, g_zzzzz_0_xyy_1, g_zzzzz_0_xyyy_1, g_zzzzz_0_xyyz_1, g_zzzzz_0_xyz_1, g_zzzzz_0_xyzz_1, g_zzzzz_0_xzz_1, g_zzzzz_0_xzzz_1, g_zzzzz_0_yyy_1, g_zzzzz_0_yyyy_1, g_zzzzz_0_yyyz_1, g_zzzzz_0_yyz_1, g_zzzzz_0_yyzz_1, g_zzzzz_0_yzz_1, g_zzzzz_0_yzzz_1, g_zzzzz_0_zzz_1, g_zzzzz_0_zzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzzzzz_0_xxxx_0[i] = 4.0 * g_zzzzz_0_xxx_1[i] * fi_acd_0 + g_zzzzz_0_xxxx_1[i] * wa_x[i];

        g_xzzzzz_0_xxxy_0[i] = 3.0 * g_zzzzz_0_xxy_1[i] * fi_acd_0 + g_zzzzz_0_xxxy_1[i] * wa_x[i];

        g_xzzzzz_0_xxxz_0[i] = 3.0 * g_zzzzz_0_xxz_1[i] * fi_acd_0 + g_zzzzz_0_xxxz_1[i] * wa_x[i];

        g_xzzzzz_0_xxyy_0[i] = 2.0 * g_zzzzz_0_xyy_1[i] * fi_acd_0 + g_zzzzz_0_xxyy_1[i] * wa_x[i];

        g_xzzzzz_0_xxyz_0[i] = 2.0 * g_zzzzz_0_xyz_1[i] * fi_acd_0 + g_zzzzz_0_xxyz_1[i] * wa_x[i];

        g_xzzzzz_0_xxzz_0[i] = 2.0 * g_zzzzz_0_xzz_1[i] * fi_acd_0 + g_zzzzz_0_xxzz_1[i] * wa_x[i];

        g_xzzzzz_0_xyyy_0[i] = g_zzzzz_0_yyy_1[i] * fi_acd_0 + g_zzzzz_0_xyyy_1[i] * wa_x[i];

        g_xzzzzz_0_xyyz_0[i] = g_zzzzz_0_yyz_1[i] * fi_acd_0 + g_zzzzz_0_xyyz_1[i] * wa_x[i];

        g_xzzzzz_0_xyzz_0[i] = g_zzzzz_0_yzz_1[i] * fi_acd_0 + g_zzzzz_0_xyzz_1[i] * wa_x[i];

        g_xzzzzz_0_xzzz_0[i] = g_zzzzz_0_zzz_1[i] * fi_acd_0 + g_zzzzz_0_xzzz_1[i] * wa_x[i];

        g_xzzzzz_0_yyyy_0[i] = g_zzzzz_0_yyyy_1[i] * wa_x[i];

        g_xzzzzz_0_yyyz_0[i] = g_zzzzz_0_yyyz_1[i] * wa_x[i];

        g_xzzzzz_0_yyzz_0[i] = g_zzzzz_0_yyzz_1[i] * wa_x[i];

        g_xzzzzz_0_yzzz_0[i] = g_zzzzz_0_yzzz_1[i] * wa_x[i];

        g_xzzzzz_0_zzzz_0[i] = g_zzzzz_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 315-330 components of targeted buffer : ISG

    auto g_yyyyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_isg + 315);

    auto g_yyyyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_isg + 316);

    auto g_yyyyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_isg + 317);

    auto g_yyyyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_isg + 318);

    auto g_yyyyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_isg + 319);

    auto g_yyyyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_isg + 320);

    auto g_yyyyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_isg + 321);

    auto g_yyyyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_isg + 322);

    auto g_yyyyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_isg + 323);

    auto g_yyyyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_isg + 324);

    auto g_yyyyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_isg + 325);

    auto g_yyyyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_isg + 326);

    auto g_yyyyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_isg + 327);

    auto g_yyyyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_isg + 328);

    auto g_yyyyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_isg + 329);

    #pragma omp simd aligned(g_yyyy_0_xxxx_0, g_yyyy_0_xxxx_1, g_yyyy_0_xxxy_0, g_yyyy_0_xxxy_1, g_yyyy_0_xxxz_0, g_yyyy_0_xxxz_1, g_yyyy_0_xxyy_0, g_yyyy_0_xxyy_1, g_yyyy_0_xxyz_0, g_yyyy_0_xxyz_1, g_yyyy_0_xxzz_0, g_yyyy_0_xxzz_1, g_yyyy_0_xyyy_0, g_yyyy_0_xyyy_1, g_yyyy_0_xyyz_0, g_yyyy_0_xyyz_1, g_yyyy_0_xyzz_0, g_yyyy_0_xyzz_1, g_yyyy_0_xzzz_0, g_yyyy_0_xzzz_1, g_yyyy_0_yyyy_0, g_yyyy_0_yyyy_1, g_yyyy_0_yyyz_0, g_yyyy_0_yyyz_1, g_yyyy_0_yyzz_0, g_yyyy_0_yyzz_1, g_yyyy_0_yzzz_0, g_yyyy_0_yzzz_1, g_yyyy_0_zzzz_0, g_yyyy_0_zzzz_1, g_yyyyy_0_xxx_1, g_yyyyy_0_xxxx_1, g_yyyyy_0_xxxy_1, g_yyyyy_0_xxxz_1, g_yyyyy_0_xxy_1, g_yyyyy_0_xxyy_1, g_yyyyy_0_xxyz_1, g_yyyyy_0_xxz_1, g_yyyyy_0_xxzz_1, g_yyyyy_0_xyy_1, g_yyyyy_0_xyyy_1, g_yyyyy_0_xyyz_1, g_yyyyy_0_xyz_1, g_yyyyy_0_xyzz_1, g_yyyyy_0_xzz_1, g_yyyyy_0_xzzz_1, g_yyyyy_0_yyy_1, g_yyyyy_0_yyyy_1, g_yyyyy_0_yyyz_1, g_yyyyy_0_yyz_1, g_yyyyy_0_yyzz_1, g_yyyyy_0_yzz_1, g_yyyyy_0_yzzz_1, g_yyyyy_0_zzz_1, g_yyyyy_0_zzzz_1, g_yyyyyy_0_xxxx_0, g_yyyyyy_0_xxxy_0, g_yyyyyy_0_xxxz_0, g_yyyyyy_0_xxyy_0, g_yyyyyy_0_xxyz_0, g_yyyyyy_0_xxzz_0, g_yyyyyy_0_xyyy_0, g_yyyyyy_0_xyyz_0, g_yyyyyy_0_xyzz_0, g_yyyyyy_0_xzzz_0, g_yyyyyy_0_yyyy_0, g_yyyyyy_0_yyyz_0, g_yyyyyy_0_yyzz_0, g_yyyyyy_0_yzzz_0, g_yyyyyy_0_zzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyyy_0_xxxx_0[i] = 5.0 * g_yyyy_0_xxxx_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxx_1[i] * fz_be_0 + g_yyyyy_0_xxxx_1[i] * wa_y[i];

        g_yyyyyy_0_xxxy_0[i] = 5.0 * g_yyyy_0_xxxy_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxy_1[i] * fz_be_0 + g_yyyyy_0_xxx_1[i] * fi_acd_0 + g_yyyyy_0_xxxy_1[i] * wa_y[i];

        g_yyyyyy_0_xxxz_0[i] = 5.0 * g_yyyy_0_xxxz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxz_1[i] * fz_be_0 + g_yyyyy_0_xxxz_1[i] * wa_y[i];

        g_yyyyyy_0_xxyy_0[i] = 5.0 * g_yyyy_0_xxyy_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxyy_1[i] * fz_be_0 + 2.0 * g_yyyyy_0_xxy_1[i] * fi_acd_0 + g_yyyyy_0_xxyy_1[i] * wa_y[i];

        g_yyyyyy_0_xxyz_0[i] = 5.0 * g_yyyy_0_xxyz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxyz_1[i] * fz_be_0 + g_yyyyy_0_xxz_1[i] * fi_acd_0 + g_yyyyy_0_xxyz_1[i] * wa_y[i];

        g_yyyyyy_0_xxzz_0[i] = 5.0 * g_yyyy_0_xxzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxzz_1[i] * fz_be_0 + g_yyyyy_0_xxzz_1[i] * wa_y[i];

        g_yyyyyy_0_xyyy_0[i] = 5.0 * g_yyyy_0_xyyy_0[i] * fbe_0 - 5.0 * g_yyyy_0_xyyy_1[i] * fz_be_0 + 3.0 * g_yyyyy_0_xyy_1[i] * fi_acd_0 + g_yyyyy_0_xyyy_1[i] * wa_y[i];

        g_yyyyyy_0_xyyz_0[i] = 5.0 * g_yyyy_0_xyyz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xyyz_1[i] * fz_be_0 + 2.0 * g_yyyyy_0_xyz_1[i] * fi_acd_0 + g_yyyyy_0_xyyz_1[i] * wa_y[i];

        g_yyyyyy_0_xyzz_0[i] = 5.0 * g_yyyy_0_xyzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xyzz_1[i] * fz_be_0 + g_yyyyy_0_xzz_1[i] * fi_acd_0 + g_yyyyy_0_xyzz_1[i] * wa_y[i];

        g_yyyyyy_0_xzzz_0[i] = 5.0 * g_yyyy_0_xzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xzzz_1[i] * fz_be_0 + g_yyyyy_0_xzzz_1[i] * wa_y[i];

        g_yyyyyy_0_yyyy_0[i] = 5.0 * g_yyyy_0_yyyy_0[i] * fbe_0 - 5.0 * g_yyyy_0_yyyy_1[i] * fz_be_0 + 4.0 * g_yyyyy_0_yyy_1[i] * fi_acd_0 + g_yyyyy_0_yyyy_1[i] * wa_y[i];

        g_yyyyyy_0_yyyz_0[i] = 5.0 * g_yyyy_0_yyyz_0[i] * fbe_0 - 5.0 * g_yyyy_0_yyyz_1[i] * fz_be_0 + 3.0 * g_yyyyy_0_yyz_1[i] * fi_acd_0 + g_yyyyy_0_yyyz_1[i] * wa_y[i];

        g_yyyyyy_0_yyzz_0[i] = 5.0 * g_yyyy_0_yyzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_yyzz_1[i] * fz_be_0 + 2.0 * g_yyyyy_0_yzz_1[i] * fi_acd_0 + g_yyyyy_0_yyzz_1[i] * wa_y[i];

        g_yyyyyy_0_yzzz_0[i] = 5.0 * g_yyyy_0_yzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_yzzz_1[i] * fz_be_0 + g_yyyyy_0_zzz_1[i] * fi_acd_0 + g_yyyyy_0_yzzz_1[i] * wa_y[i];

        g_yyyyyy_0_zzzz_0[i] = 5.0 * g_yyyy_0_zzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_zzzz_1[i] * fz_be_0 + g_yyyyy_0_zzzz_1[i] * wa_y[i];
    }

    /// Set up 330-345 components of targeted buffer : ISG

    auto g_yyyyyz_0_xxxx_0 = pbuffer.data(idx_eri_0_isg + 330);

    auto g_yyyyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_isg + 331);

    auto g_yyyyyz_0_xxxz_0 = pbuffer.data(idx_eri_0_isg + 332);

    auto g_yyyyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_isg + 333);

    auto g_yyyyyz_0_xxyz_0 = pbuffer.data(idx_eri_0_isg + 334);

    auto g_yyyyyz_0_xxzz_0 = pbuffer.data(idx_eri_0_isg + 335);

    auto g_yyyyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_isg + 336);

    auto g_yyyyyz_0_xyyz_0 = pbuffer.data(idx_eri_0_isg + 337);

    auto g_yyyyyz_0_xyzz_0 = pbuffer.data(idx_eri_0_isg + 338);

    auto g_yyyyyz_0_xzzz_0 = pbuffer.data(idx_eri_0_isg + 339);

    auto g_yyyyyz_0_yyyy_0 = pbuffer.data(idx_eri_0_isg + 340);

    auto g_yyyyyz_0_yyyz_0 = pbuffer.data(idx_eri_0_isg + 341);

    auto g_yyyyyz_0_yyzz_0 = pbuffer.data(idx_eri_0_isg + 342);

    auto g_yyyyyz_0_yzzz_0 = pbuffer.data(idx_eri_0_isg + 343);

    auto g_yyyyyz_0_zzzz_0 = pbuffer.data(idx_eri_0_isg + 344);

    #pragma omp simd aligned(g_yyyyy_0_xxx_1, g_yyyyy_0_xxxx_1, g_yyyyy_0_xxxy_1, g_yyyyy_0_xxxz_1, g_yyyyy_0_xxy_1, g_yyyyy_0_xxyy_1, g_yyyyy_0_xxyz_1, g_yyyyy_0_xxz_1, g_yyyyy_0_xxzz_1, g_yyyyy_0_xyy_1, g_yyyyy_0_xyyy_1, g_yyyyy_0_xyyz_1, g_yyyyy_0_xyz_1, g_yyyyy_0_xyzz_1, g_yyyyy_0_xzz_1, g_yyyyy_0_xzzz_1, g_yyyyy_0_yyy_1, g_yyyyy_0_yyyy_1, g_yyyyy_0_yyyz_1, g_yyyyy_0_yyz_1, g_yyyyy_0_yyzz_1, g_yyyyy_0_yzz_1, g_yyyyy_0_yzzz_1, g_yyyyy_0_zzz_1, g_yyyyy_0_zzzz_1, g_yyyyyz_0_xxxx_0, g_yyyyyz_0_xxxy_0, g_yyyyyz_0_xxxz_0, g_yyyyyz_0_xxyy_0, g_yyyyyz_0_xxyz_0, g_yyyyyz_0_xxzz_0, g_yyyyyz_0_xyyy_0, g_yyyyyz_0_xyyz_0, g_yyyyyz_0_xyzz_0, g_yyyyyz_0_xzzz_0, g_yyyyyz_0_yyyy_0, g_yyyyyz_0_yyyz_0, g_yyyyyz_0_yyzz_0, g_yyyyyz_0_yzzz_0, g_yyyyyz_0_zzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyyyyz_0_xxxx_0[i] = g_yyyyy_0_xxxx_1[i] * wa_z[i];

        g_yyyyyz_0_xxxy_0[i] = g_yyyyy_0_xxxy_1[i] * wa_z[i];

        g_yyyyyz_0_xxxz_0[i] = g_yyyyy_0_xxx_1[i] * fi_acd_0 + g_yyyyy_0_xxxz_1[i] * wa_z[i];

        g_yyyyyz_0_xxyy_0[i] = g_yyyyy_0_xxyy_1[i] * wa_z[i];

        g_yyyyyz_0_xxyz_0[i] = g_yyyyy_0_xxy_1[i] * fi_acd_0 + g_yyyyy_0_xxyz_1[i] * wa_z[i];

        g_yyyyyz_0_xxzz_0[i] = 2.0 * g_yyyyy_0_xxz_1[i] * fi_acd_0 + g_yyyyy_0_xxzz_1[i] * wa_z[i];

        g_yyyyyz_0_xyyy_0[i] = g_yyyyy_0_xyyy_1[i] * wa_z[i];

        g_yyyyyz_0_xyyz_0[i] = g_yyyyy_0_xyy_1[i] * fi_acd_0 + g_yyyyy_0_xyyz_1[i] * wa_z[i];

        g_yyyyyz_0_xyzz_0[i] = 2.0 * g_yyyyy_0_xyz_1[i] * fi_acd_0 + g_yyyyy_0_xyzz_1[i] * wa_z[i];

        g_yyyyyz_0_xzzz_0[i] = 3.0 * g_yyyyy_0_xzz_1[i] * fi_acd_0 + g_yyyyy_0_xzzz_1[i] * wa_z[i];

        g_yyyyyz_0_yyyy_0[i] = g_yyyyy_0_yyyy_1[i] * wa_z[i];

        g_yyyyyz_0_yyyz_0[i] = g_yyyyy_0_yyy_1[i] * fi_acd_0 + g_yyyyy_0_yyyz_1[i] * wa_z[i];

        g_yyyyyz_0_yyzz_0[i] = 2.0 * g_yyyyy_0_yyz_1[i] * fi_acd_0 + g_yyyyy_0_yyzz_1[i] * wa_z[i];

        g_yyyyyz_0_yzzz_0[i] = 3.0 * g_yyyyy_0_yzz_1[i] * fi_acd_0 + g_yyyyy_0_yzzz_1[i] * wa_z[i];

        g_yyyyyz_0_zzzz_0[i] = 4.0 * g_yyyyy_0_zzz_1[i] * fi_acd_0 + g_yyyyy_0_zzzz_1[i] * wa_z[i];
    }

    /// Set up 345-360 components of targeted buffer : ISG

    auto g_yyyyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_isg + 345);

    auto g_yyyyzz_0_xxxy_0 = pbuffer.data(idx_eri_0_isg + 346);

    auto g_yyyyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_isg + 347);

    auto g_yyyyzz_0_xxyy_0 = pbuffer.data(idx_eri_0_isg + 348);

    auto g_yyyyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_isg + 349);

    auto g_yyyyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_isg + 350);

    auto g_yyyyzz_0_xyyy_0 = pbuffer.data(idx_eri_0_isg + 351);

    auto g_yyyyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_isg + 352);

    auto g_yyyyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_isg + 353);

    auto g_yyyyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_isg + 354);

    auto g_yyyyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_isg + 355);

    auto g_yyyyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_isg + 356);

    auto g_yyyyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_isg + 357);

    auto g_yyyyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_isg + 358);

    auto g_yyyyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_isg + 359);

    #pragma omp simd aligned(g_yyyy_0_xxxy_0, g_yyyy_0_xxxy_1, g_yyyy_0_xxyy_0, g_yyyy_0_xxyy_1, g_yyyy_0_xyyy_0, g_yyyy_0_xyyy_1, g_yyyy_0_yyyy_0, g_yyyy_0_yyyy_1, g_yyyyz_0_xxxy_1, g_yyyyz_0_xxyy_1, g_yyyyz_0_xyyy_1, g_yyyyz_0_yyyy_1, g_yyyyzz_0_xxxx_0, g_yyyyzz_0_xxxy_0, g_yyyyzz_0_xxxz_0, g_yyyyzz_0_xxyy_0, g_yyyyzz_0_xxyz_0, g_yyyyzz_0_xxzz_0, g_yyyyzz_0_xyyy_0, g_yyyyzz_0_xyyz_0, g_yyyyzz_0_xyzz_0, g_yyyyzz_0_xzzz_0, g_yyyyzz_0_yyyy_0, g_yyyyzz_0_yyyz_0, g_yyyyzz_0_yyzz_0, g_yyyyzz_0_yzzz_0, g_yyyyzz_0_zzzz_0, g_yyyzz_0_xxxx_1, g_yyyzz_0_xxxz_1, g_yyyzz_0_xxyz_1, g_yyyzz_0_xxz_1, g_yyyzz_0_xxzz_1, g_yyyzz_0_xyyz_1, g_yyyzz_0_xyz_1, g_yyyzz_0_xyzz_1, g_yyyzz_0_xzz_1, g_yyyzz_0_xzzz_1, g_yyyzz_0_yyyz_1, g_yyyzz_0_yyz_1, g_yyyzz_0_yyzz_1, g_yyyzz_0_yzz_1, g_yyyzz_0_yzzz_1, g_yyyzz_0_zzz_1, g_yyyzz_0_zzzz_1, g_yyzz_0_xxxx_0, g_yyzz_0_xxxx_1, g_yyzz_0_xxxz_0, g_yyzz_0_xxxz_1, g_yyzz_0_xxyz_0, g_yyzz_0_xxyz_1, g_yyzz_0_xxzz_0, g_yyzz_0_xxzz_1, g_yyzz_0_xyyz_0, g_yyzz_0_xyyz_1, g_yyzz_0_xyzz_0, g_yyzz_0_xyzz_1, g_yyzz_0_xzzz_0, g_yyzz_0_xzzz_1, g_yyzz_0_yyyz_0, g_yyzz_0_yyyz_1, g_yyzz_0_yyzz_0, g_yyzz_0_yyzz_1, g_yyzz_0_yzzz_0, g_yyzz_0_yzzz_1, g_yyzz_0_zzzz_0, g_yyzz_0_zzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyzz_0_xxxx_0[i] = 3.0 * g_yyzz_0_xxxx_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxx_1[i] * fz_be_0 + g_yyyzz_0_xxxx_1[i] * wa_y[i];

        g_yyyyzz_0_xxxy_0[i] = g_yyyy_0_xxxy_0[i] * fbe_0 - g_yyyy_0_xxxy_1[i] * fz_be_0 + g_yyyyz_0_xxxy_1[i] * wa_z[i];

        g_yyyyzz_0_xxxz_0[i] = 3.0 * g_yyzz_0_xxxz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxz_1[i] * fz_be_0 + g_yyyzz_0_xxxz_1[i] * wa_y[i];

        g_yyyyzz_0_xxyy_0[i] = g_yyyy_0_xxyy_0[i] * fbe_0 - g_yyyy_0_xxyy_1[i] * fz_be_0 + g_yyyyz_0_xxyy_1[i] * wa_z[i];

        g_yyyyzz_0_xxyz_0[i] = 3.0 * g_yyzz_0_xxyz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxyz_1[i] * fz_be_0 + g_yyyzz_0_xxz_1[i] * fi_acd_0 + g_yyyzz_0_xxyz_1[i] * wa_y[i];

        g_yyyyzz_0_xxzz_0[i] = 3.0 * g_yyzz_0_xxzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxzz_1[i] * fz_be_0 + g_yyyzz_0_xxzz_1[i] * wa_y[i];

        g_yyyyzz_0_xyyy_0[i] = g_yyyy_0_xyyy_0[i] * fbe_0 - g_yyyy_0_xyyy_1[i] * fz_be_0 + g_yyyyz_0_xyyy_1[i] * wa_z[i];

        g_yyyyzz_0_xyyz_0[i] = 3.0 * g_yyzz_0_xyyz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xyyz_1[i] * fz_be_0 + 2.0 * g_yyyzz_0_xyz_1[i] * fi_acd_0 + g_yyyzz_0_xyyz_1[i] * wa_y[i];

        g_yyyyzz_0_xyzz_0[i] = 3.0 * g_yyzz_0_xyzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xyzz_1[i] * fz_be_0 + g_yyyzz_0_xzz_1[i] * fi_acd_0 + g_yyyzz_0_xyzz_1[i] * wa_y[i];

        g_yyyyzz_0_xzzz_0[i] = 3.0 * g_yyzz_0_xzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xzzz_1[i] * fz_be_0 + g_yyyzz_0_xzzz_1[i] * wa_y[i];

        g_yyyyzz_0_yyyy_0[i] = g_yyyy_0_yyyy_0[i] * fbe_0 - g_yyyy_0_yyyy_1[i] * fz_be_0 + g_yyyyz_0_yyyy_1[i] * wa_z[i];

        g_yyyyzz_0_yyyz_0[i] = 3.0 * g_yyzz_0_yyyz_0[i] * fbe_0 - 3.0 * g_yyzz_0_yyyz_1[i] * fz_be_0 + 3.0 * g_yyyzz_0_yyz_1[i] * fi_acd_0 + g_yyyzz_0_yyyz_1[i] * wa_y[i];

        g_yyyyzz_0_yyzz_0[i] = 3.0 * g_yyzz_0_yyzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_yyzz_1[i] * fz_be_0 + 2.0 * g_yyyzz_0_yzz_1[i] * fi_acd_0 + g_yyyzz_0_yyzz_1[i] * wa_y[i];

        g_yyyyzz_0_yzzz_0[i] = 3.0 * g_yyzz_0_yzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_yzzz_1[i] * fz_be_0 + g_yyyzz_0_zzz_1[i] * fi_acd_0 + g_yyyzz_0_yzzz_1[i] * wa_y[i];

        g_yyyyzz_0_zzzz_0[i] = 3.0 * g_yyzz_0_zzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_zzzz_1[i] * fz_be_0 + g_yyyzz_0_zzzz_1[i] * wa_y[i];
    }

    /// Set up 360-375 components of targeted buffer : ISG

    auto g_yyyzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_isg + 360);

    auto g_yyyzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_isg + 361);

    auto g_yyyzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_isg + 362);

    auto g_yyyzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_isg + 363);

    auto g_yyyzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_isg + 364);

    auto g_yyyzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_isg + 365);

    auto g_yyyzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_isg + 366);

    auto g_yyyzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_isg + 367);

    auto g_yyyzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_isg + 368);

    auto g_yyyzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_isg + 369);

    auto g_yyyzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_isg + 370);

    auto g_yyyzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_isg + 371);

    auto g_yyyzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_isg + 372);

    auto g_yyyzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_isg + 373);

    auto g_yyyzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_isg + 374);

    #pragma omp simd aligned(g_yyyz_0_xxxy_0, g_yyyz_0_xxxy_1, g_yyyz_0_xxyy_0, g_yyyz_0_xxyy_1, g_yyyz_0_xyyy_0, g_yyyz_0_xyyy_1, g_yyyz_0_yyyy_0, g_yyyz_0_yyyy_1, g_yyyzz_0_xxxy_1, g_yyyzz_0_xxyy_1, g_yyyzz_0_xyyy_1, g_yyyzz_0_yyyy_1, g_yyyzzz_0_xxxx_0, g_yyyzzz_0_xxxy_0, g_yyyzzz_0_xxxz_0, g_yyyzzz_0_xxyy_0, g_yyyzzz_0_xxyz_0, g_yyyzzz_0_xxzz_0, g_yyyzzz_0_xyyy_0, g_yyyzzz_0_xyyz_0, g_yyyzzz_0_xyzz_0, g_yyyzzz_0_xzzz_0, g_yyyzzz_0_yyyy_0, g_yyyzzz_0_yyyz_0, g_yyyzzz_0_yyzz_0, g_yyyzzz_0_yzzz_0, g_yyyzzz_0_zzzz_0, g_yyzzz_0_xxxx_1, g_yyzzz_0_xxxz_1, g_yyzzz_0_xxyz_1, g_yyzzz_0_xxz_1, g_yyzzz_0_xxzz_1, g_yyzzz_0_xyyz_1, g_yyzzz_0_xyz_1, g_yyzzz_0_xyzz_1, g_yyzzz_0_xzz_1, g_yyzzz_0_xzzz_1, g_yyzzz_0_yyyz_1, g_yyzzz_0_yyz_1, g_yyzzz_0_yyzz_1, g_yyzzz_0_yzz_1, g_yyzzz_0_yzzz_1, g_yyzzz_0_zzz_1, g_yyzzz_0_zzzz_1, g_yzzz_0_xxxx_0, g_yzzz_0_xxxx_1, g_yzzz_0_xxxz_0, g_yzzz_0_xxxz_1, g_yzzz_0_xxyz_0, g_yzzz_0_xxyz_1, g_yzzz_0_xxzz_0, g_yzzz_0_xxzz_1, g_yzzz_0_xyyz_0, g_yzzz_0_xyyz_1, g_yzzz_0_xyzz_0, g_yzzz_0_xyzz_1, g_yzzz_0_xzzz_0, g_yzzz_0_xzzz_1, g_yzzz_0_yyyz_0, g_yzzz_0_yyyz_1, g_yzzz_0_yyzz_0, g_yzzz_0_yyzz_1, g_yzzz_0_yzzz_0, g_yzzz_0_yzzz_1, g_yzzz_0_zzzz_0, g_yzzz_0_zzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyzzz_0_xxxx_0[i] = 2.0 * g_yzzz_0_xxxx_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxxx_1[i] * fz_be_0 + g_yyzzz_0_xxxx_1[i] * wa_y[i];

        g_yyyzzz_0_xxxy_0[i] = 2.0 * g_yyyz_0_xxxy_0[i] * fbe_0 - 2.0 * g_yyyz_0_xxxy_1[i] * fz_be_0 + g_yyyzz_0_xxxy_1[i] * wa_z[i];

        g_yyyzzz_0_xxxz_0[i] = 2.0 * g_yzzz_0_xxxz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxxz_1[i] * fz_be_0 + g_yyzzz_0_xxxz_1[i] * wa_y[i];

        g_yyyzzz_0_xxyy_0[i] = 2.0 * g_yyyz_0_xxyy_0[i] * fbe_0 - 2.0 * g_yyyz_0_xxyy_1[i] * fz_be_0 + g_yyyzz_0_xxyy_1[i] * wa_z[i];

        g_yyyzzz_0_xxyz_0[i] = 2.0 * g_yzzz_0_xxyz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxyz_1[i] * fz_be_0 + g_yyzzz_0_xxz_1[i] * fi_acd_0 + g_yyzzz_0_xxyz_1[i] * wa_y[i];

        g_yyyzzz_0_xxzz_0[i] = 2.0 * g_yzzz_0_xxzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxzz_1[i] * fz_be_0 + g_yyzzz_0_xxzz_1[i] * wa_y[i];

        g_yyyzzz_0_xyyy_0[i] = 2.0 * g_yyyz_0_xyyy_0[i] * fbe_0 - 2.0 * g_yyyz_0_xyyy_1[i] * fz_be_0 + g_yyyzz_0_xyyy_1[i] * wa_z[i];

        g_yyyzzz_0_xyyz_0[i] = 2.0 * g_yzzz_0_xyyz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xyyz_1[i] * fz_be_0 + 2.0 * g_yyzzz_0_xyz_1[i] * fi_acd_0 + g_yyzzz_0_xyyz_1[i] * wa_y[i];

        g_yyyzzz_0_xyzz_0[i] = 2.0 * g_yzzz_0_xyzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xyzz_1[i] * fz_be_0 + g_yyzzz_0_xzz_1[i] * fi_acd_0 + g_yyzzz_0_xyzz_1[i] * wa_y[i];

        g_yyyzzz_0_xzzz_0[i] = 2.0 * g_yzzz_0_xzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xzzz_1[i] * fz_be_0 + g_yyzzz_0_xzzz_1[i] * wa_y[i];

        g_yyyzzz_0_yyyy_0[i] = 2.0 * g_yyyz_0_yyyy_0[i] * fbe_0 - 2.0 * g_yyyz_0_yyyy_1[i] * fz_be_0 + g_yyyzz_0_yyyy_1[i] * wa_z[i];

        g_yyyzzz_0_yyyz_0[i] = 2.0 * g_yzzz_0_yyyz_0[i] * fbe_0 - 2.0 * g_yzzz_0_yyyz_1[i] * fz_be_0 + 3.0 * g_yyzzz_0_yyz_1[i] * fi_acd_0 + g_yyzzz_0_yyyz_1[i] * wa_y[i];

        g_yyyzzz_0_yyzz_0[i] = 2.0 * g_yzzz_0_yyzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_yyzz_1[i] * fz_be_0 + 2.0 * g_yyzzz_0_yzz_1[i] * fi_acd_0 + g_yyzzz_0_yyzz_1[i] * wa_y[i];

        g_yyyzzz_0_yzzz_0[i] = 2.0 * g_yzzz_0_yzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_yzzz_1[i] * fz_be_0 + g_yyzzz_0_zzz_1[i] * fi_acd_0 + g_yyzzz_0_yzzz_1[i] * wa_y[i];

        g_yyyzzz_0_zzzz_0[i] = 2.0 * g_yzzz_0_zzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_zzzz_1[i] * fz_be_0 + g_yyzzz_0_zzzz_1[i] * wa_y[i];
    }

    /// Set up 375-390 components of targeted buffer : ISG

    auto g_yyzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_isg + 375);

    auto g_yyzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_isg + 376);

    auto g_yyzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_isg + 377);

    auto g_yyzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_isg + 378);

    auto g_yyzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_isg + 379);

    auto g_yyzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_isg + 380);

    auto g_yyzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_isg + 381);

    auto g_yyzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_isg + 382);

    auto g_yyzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_isg + 383);

    auto g_yyzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_isg + 384);

    auto g_yyzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_isg + 385);

    auto g_yyzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_isg + 386);

    auto g_yyzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_isg + 387);

    auto g_yyzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_isg + 388);

    auto g_yyzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_isg + 389);

    #pragma omp simd aligned(g_yyzz_0_xxxy_0, g_yyzz_0_xxxy_1, g_yyzz_0_xxyy_0, g_yyzz_0_xxyy_1, g_yyzz_0_xyyy_0, g_yyzz_0_xyyy_1, g_yyzz_0_yyyy_0, g_yyzz_0_yyyy_1, g_yyzzz_0_xxxy_1, g_yyzzz_0_xxyy_1, g_yyzzz_0_xyyy_1, g_yyzzz_0_yyyy_1, g_yyzzzz_0_xxxx_0, g_yyzzzz_0_xxxy_0, g_yyzzzz_0_xxxz_0, g_yyzzzz_0_xxyy_0, g_yyzzzz_0_xxyz_0, g_yyzzzz_0_xxzz_0, g_yyzzzz_0_xyyy_0, g_yyzzzz_0_xyyz_0, g_yyzzzz_0_xyzz_0, g_yyzzzz_0_xzzz_0, g_yyzzzz_0_yyyy_0, g_yyzzzz_0_yyyz_0, g_yyzzzz_0_yyzz_0, g_yyzzzz_0_yzzz_0, g_yyzzzz_0_zzzz_0, g_yzzzz_0_xxxx_1, g_yzzzz_0_xxxz_1, g_yzzzz_0_xxyz_1, g_yzzzz_0_xxz_1, g_yzzzz_0_xxzz_1, g_yzzzz_0_xyyz_1, g_yzzzz_0_xyz_1, g_yzzzz_0_xyzz_1, g_yzzzz_0_xzz_1, g_yzzzz_0_xzzz_1, g_yzzzz_0_yyyz_1, g_yzzzz_0_yyz_1, g_yzzzz_0_yyzz_1, g_yzzzz_0_yzz_1, g_yzzzz_0_yzzz_1, g_yzzzz_0_zzz_1, g_yzzzz_0_zzzz_1, g_zzzz_0_xxxx_0, g_zzzz_0_xxxx_1, g_zzzz_0_xxxz_0, g_zzzz_0_xxxz_1, g_zzzz_0_xxyz_0, g_zzzz_0_xxyz_1, g_zzzz_0_xxzz_0, g_zzzz_0_xxzz_1, g_zzzz_0_xyyz_0, g_zzzz_0_xyyz_1, g_zzzz_0_xyzz_0, g_zzzz_0_xyzz_1, g_zzzz_0_xzzz_0, g_zzzz_0_xzzz_1, g_zzzz_0_yyyz_0, g_zzzz_0_yyyz_1, g_zzzz_0_yyzz_0, g_zzzz_0_yyzz_1, g_zzzz_0_yzzz_0, g_zzzz_0_yzzz_1, g_zzzz_0_zzzz_0, g_zzzz_0_zzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyzzzz_0_xxxx_0[i] = g_zzzz_0_xxxx_0[i] * fbe_0 - g_zzzz_0_xxxx_1[i] * fz_be_0 + g_yzzzz_0_xxxx_1[i] * wa_y[i];

        g_yyzzzz_0_xxxy_0[i] = 3.0 * g_yyzz_0_xxxy_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxy_1[i] * fz_be_0 + g_yyzzz_0_xxxy_1[i] * wa_z[i];

        g_yyzzzz_0_xxxz_0[i] = g_zzzz_0_xxxz_0[i] * fbe_0 - g_zzzz_0_xxxz_1[i] * fz_be_0 + g_yzzzz_0_xxxz_1[i] * wa_y[i];

        g_yyzzzz_0_xxyy_0[i] = 3.0 * g_yyzz_0_xxyy_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxyy_1[i] * fz_be_0 + g_yyzzz_0_xxyy_1[i] * wa_z[i];

        g_yyzzzz_0_xxyz_0[i] = g_zzzz_0_xxyz_0[i] * fbe_0 - g_zzzz_0_xxyz_1[i] * fz_be_0 + g_yzzzz_0_xxz_1[i] * fi_acd_0 + g_yzzzz_0_xxyz_1[i] * wa_y[i];

        g_yyzzzz_0_xxzz_0[i] = g_zzzz_0_xxzz_0[i] * fbe_0 - g_zzzz_0_xxzz_1[i] * fz_be_0 + g_yzzzz_0_xxzz_1[i] * wa_y[i];

        g_yyzzzz_0_xyyy_0[i] = 3.0 * g_yyzz_0_xyyy_0[i] * fbe_0 - 3.0 * g_yyzz_0_xyyy_1[i] * fz_be_0 + g_yyzzz_0_xyyy_1[i] * wa_z[i];

        g_yyzzzz_0_xyyz_0[i] = g_zzzz_0_xyyz_0[i] * fbe_0 - g_zzzz_0_xyyz_1[i] * fz_be_0 + 2.0 * g_yzzzz_0_xyz_1[i] * fi_acd_0 + g_yzzzz_0_xyyz_1[i] * wa_y[i];

        g_yyzzzz_0_xyzz_0[i] = g_zzzz_0_xyzz_0[i] * fbe_0 - g_zzzz_0_xyzz_1[i] * fz_be_0 + g_yzzzz_0_xzz_1[i] * fi_acd_0 + g_yzzzz_0_xyzz_1[i] * wa_y[i];

        g_yyzzzz_0_xzzz_0[i] = g_zzzz_0_xzzz_0[i] * fbe_0 - g_zzzz_0_xzzz_1[i] * fz_be_0 + g_yzzzz_0_xzzz_1[i] * wa_y[i];

        g_yyzzzz_0_yyyy_0[i] = 3.0 * g_yyzz_0_yyyy_0[i] * fbe_0 - 3.0 * g_yyzz_0_yyyy_1[i] * fz_be_0 + g_yyzzz_0_yyyy_1[i] * wa_z[i];

        g_yyzzzz_0_yyyz_0[i] = g_zzzz_0_yyyz_0[i] * fbe_0 - g_zzzz_0_yyyz_1[i] * fz_be_0 + 3.0 * g_yzzzz_0_yyz_1[i] * fi_acd_0 + g_yzzzz_0_yyyz_1[i] * wa_y[i];

        g_yyzzzz_0_yyzz_0[i] = g_zzzz_0_yyzz_0[i] * fbe_0 - g_zzzz_0_yyzz_1[i] * fz_be_0 + 2.0 * g_yzzzz_0_yzz_1[i] * fi_acd_0 + g_yzzzz_0_yyzz_1[i] * wa_y[i];

        g_yyzzzz_0_yzzz_0[i] = g_zzzz_0_yzzz_0[i] * fbe_0 - g_zzzz_0_yzzz_1[i] * fz_be_0 + g_yzzzz_0_zzz_1[i] * fi_acd_0 + g_yzzzz_0_yzzz_1[i] * wa_y[i];

        g_yyzzzz_0_zzzz_0[i] = g_zzzz_0_zzzz_0[i] * fbe_0 - g_zzzz_0_zzzz_1[i] * fz_be_0 + g_yzzzz_0_zzzz_1[i] * wa_y[i];
    }

    /// Set up 390-405 components of targeted buffer : ISG

    auto g_yzzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_isg + 390);

    auto g_yzzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_isg + 391);

    auto g_yzzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_isg + 392);

    auto g_yzzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_isg + 393);

    auto g_yzzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_isg + 394);

    auto g_yzzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_isg + 395);

    auto g_yzzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_isg + 396);

    auto g_yzzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_isg + 397);

    auto g_yzzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_isg + 398);

    auto g_yzzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_isg + 399);

    auto g_yzzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_isg + 400);

    auto g_yzzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_isg + 401);

    auto g_yzzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_isg + 402);

    auto g_yzzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_isg + 403);

    auto g_yzzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_isg + 404);

    #pragma omp simd aligned(g_yzzzzz_0_xxxx_0, g_yzzzzz_0_xxxy_0, g_yzzzzz_0_xxxz_0, g_yzzzzz_0_xxyy_0, g_yzzzzz_0_xxyz_0, g_yzzzzz_0_xxzz_0, g_yzzzzz_0_xyyy_0, g_yzzzzz_0_xyyz_0, g_yzzzzz_0_xyzz_0, g_yzzzzz_0_xzzz_0, g_yzzzzz_0_yyyy_0, g_yzzzzz_0_yyyz_0, g_yzzzzz_0_yyzz_0, g_yzzzzz_0_yzzz_0, g_yzzzzz_0_zzzz_0, g_zzzzz_0_xxx_1, g_zzzzz_0_xxxx_1, g_zzzzz_0_xxxy_1, g_zzzzz_0_xxxz_1, g_zzzzz_0_xxy_1, g_zzzzz_0_xxyy_1, g_zzzzz_0_xxyz_1, g_zzzzz_0_xxz_1, g_zzzzz_0_xxzz_1, g_zzzzz_0_xyy_1, g_zzzzz_0_xyyy_1, g_zzzzz_0_xyyz_1, g_zzzzz_0_xyz_1, g_zzzzz_0_xyzz_1, g_zzzzz_0_xzz_1, g_zzzzz_0_xzzz_1, g_zzzzz_0_yyy_1, g_zzzzz_0_yyyy_1, g_zzzzz_0_yyyz_1, g_zzzzz_0_yyz_1, g_zzzzz_0_yyzz_1, g_zzzzz_0_yzz_1, g_zzzzz_0_yzzz_1, g_zzzzz_0_zzz_1, g_zzzzz_0_zzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzzzzz_0_xxxx_0[i] = g_zzzzz_0_xxxx_1[i] * wa_y[i];

        g_yzzzzz_0_xxxy_0[i] = g_zzzzz_0_xxx_1[i] * fi_acd_0 + g_zzzzz_0_xxxy_1[i] * wa_y[i];

        g_yzzzzz_0_xxxz_0[i] = g_zzzzz_0_xxxz_1[i] * wa_y[i];

        g_yzzzzz_0_xxyy_0[i] = 2.0 * g_zzzzz_0_xxy_1[i] * fi_acd_0 + g_zzzzz_0_xxyy_1[i] * wa_y[i];

        g_yzzzzz_0_xxyz_0[i] = g_zzzzz_0_xxz_1[i] * fi_acd_0 + g_zzzzz_0_xxyz_1[i] * wa_y[i];

        g_yzzzzz_0_xxzz_0[i] = g_zzzzz_0_xxzz_1[i] * wa_y[i];

        g_yzzzzz_0_xyyy_0[i] = 3.0 * g_zzzzz_0_xyy_1[i] * fi_acd_0 + g_zzzzz_0_xyyy_1[i] * wa_y[i];

        g_yzzzzz_0_xyyz_0[i] = 2.0 * g_zzzzz_0_xyz_1[i] * fi_acd_0 + g_zzzzz_0_xyyz_1[i] * wa_y[i];

        g_yzzzzz_0_xyzz_0[i] = g_zzzzz_0_xzz_1[i] * fi_acd_0 + g_zzzzz_0_xyzz_1[i] * wa_y[i];

        g_yzzzzz_0_xzzz_0[i] = g_zzzzz_0_xzzz_1[i] * wa_y[i];

        g_yzzzzz_0_yyyy_0[i] = 4.0 * g_zzzzz_0_yyy_1[i] * fi_acd_0 + g_zzzzz_0_yyyy_1[i] * wa_y[i];

        g_yzzzzz_0_yyyz_0[i] = 3.0 * g_zzzzz_0_yyz_1[i] * fi_acd_0 + g_zzzzz_0_yyyz_1[i] * wa_y[i];

        g_yzzzzz_0_yyzz_0[i] = 2.0 * g_zzzzz_0_yzz_1[i] * fi_acd_0 + g_zzzzz_0_yyzz_1[i] * wa_y[i];

        g_yzzzzz_0_yzzz_0[i] = g_zzzzz_0_zzz_1[i] * fi_acd_0 + g_zzzzz_0_yzzz_1[i] * wa_y[i];

        g_yzzzzz_0_zzzz_0[i] = g_zzzzz_0_zzzz_1[i] * wa_y[i];
    }

    /// Set up 405-420 components of targeted buffer : ISG

    auto g_zzzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_isg + 405);

    auto g_zzzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_isg + 406);

    auto g_zzzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_isg + 407);

    auto g_zzzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_isg + 408);

    auto g_zzzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_isg + 409);

    auto g_zzzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_isg + 410);

    auto g_zzzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_isg + 411);

    auto g_zzzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_isg + 412);

    auto g_zzzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_isg + 413);

    auto g_zzzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_isg + 414);

    auto g_zzzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_isg + 415);

    auto g_zzzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_isg + 416);

    auto g_zzzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_isg + 417);

    auto g_zzzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_isg + 418);

    auto g_zzzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_isg + 419);

    #pragma omp simd aligned(g_zzzz_0_xxxx_0, g_zzzz_0_xxxx_1, g_zzzz_0_xxxy_0, g_zzzz_0_xxxy_1, g_zzzz_0_xxxz_0, g_zzzz_0_xxxz_1, g_zzzz_0_xxyy_0, g_zzzz_0_xxyy_1, g_zzzz_0_xxyz_0, g_zzzz_0_xxyz_1, g_zzzz_0_xxzz_0, g_zzzz_0_xxzz_1, g_zzzz_0_xyyy_0, g_zzzz_0_xyyy_1, g_zzzz_0_xyyz_0, g_zzzz_0_xyyz_1, g_zzzz_0_xyzz_0, g_zzzz_0_xyzz_1, g_zzzz_0_xzzz_0, g_zzzz_0_xzzz_1, g_zzzz_0_yyyy_0, g_zzzz_0_yyyy_1, g_zzzz_0_yyyz_0, g_zzzz_0_yyyz_1, g_zzzz_0_yyzz_0, g_zzzz_0_yyzz_1, g_zzzz_0_yzzz_0, g_zzzz_0_yzzz_1, g_zzzz_0_zzzz_0, g_zzzz_0_zzzz_1, g_zzzzz_0_xxx_1, g_zzzzz_0_xxxx_1, g_zzzzz_0_xxxy_1, g_zzzzz_0_xxxz_1, g_zzzzz_0_xxy_1, g_zzzzz_0_xxyy_1, g_zzzzz_0_xxyz_1, g_zzzzz_0_xxz_1, g_zzzzz_0_xxzz_1, g_zzzzz_0_xyy_1, g_zzzzz_0_xyyy_1, g_zzzzz_0_xyyz_1, g_zzzzz_0_xyz_1, g_zzzzz_0_xyzz_1, g_zzzzz_0_xzz_1, g_zzzzz_0_xzzz_1, g_zzzzz_0_yyy_1, g_zzzzz_0_yyyy_1, g_zzzzz_0_yyyz_1, g_zzzzz_0_yyz_1, g_zzzzz_0_yyzz_1, g_zzzzz_0_yzz_1, g_zzzzz_0_yzzz_1, g_zzzzz_0_zzz_1, g_zzzzz_0_zzzz_1, g_zzzzzz_0_xxxx_0, g_zzzzzz_0_xxxy_0, g_zzzzzz_0_xxxz_0, g_zzzzzz_0_xxyy_0, g_zzzzzz_0_xxyz_0, g_zzzzzz_0_xxzz_0, g_zzzzzz_0_xyyy_0, g_zzzzzz_0_xyyz_0, g_zzzzzz_0_xyzz_0, g_zzzzzz_0_xzzz_0, g_zzzzzz_0_yyyy_0, g_zzzzzz_0_yyyz_0, g_zzzzzz_0_yyzz_0, g_zzzzzz_0_yzzz_0, g_zzzzzz_0_zzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzzzzz_0_xxxx_0[i] = 5.0 * g_zzzz_0_xxxx_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxx_1[i] * fz_be_0 + g_zzzzz_0_xxxx_1[i] * wa_z[i];

        g_zzzzzz_0_xxxy_0[i] = 5.0 * g_zzzz_0_xxxy_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxy_1[i] * fz_be_0 + g_zzzzz_0_xxxy_1[i] * wa_z[i];

        g_zzzzzz_0_xxxz_0[i] = 5.0 * g_zzzz_0_xxxz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxz_1[i] * fz_be_0 + g_zzzzz_0_xxx_1[i] * fi_acd_0 + g_zzzzz_0_xxxz_1[i] * wa_z[i];

        g_zzzzzz_0_xxyy_0[i] = 5.0 * g_zzzz_0_xxyy_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxyy_1[i] * fz_be_0 + g_zzzzz_0_xxyy_1[i] * wa_z[i];

        g_zzzzzz_0_xxyz_0[i] = 5.0 * g_zzzz_0_xxyz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxyz_1[i] * fz_be_0 + g_zzzzz_0_xxy_1[i] * fi_acd_0 + g_zzzzz_0_xxyz_1[i] * wa_z[i];

        g_zzzzzz_0_xxzz_0[i] = 5.0 * g_zzzz_0_xxzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_0_xxz_1[i] * fi_acd_0 + g_zzzzz_0_xxzz_1[i] * wa_z[i];

        g_zzzzzz_0_xyyy_0[i] = 5.0 * g_zzzz_0_xyyy_0[i] * fbe_0 - 5.0 * g_zzzz_0_xyyy_1[i] * fz_be_0 + g_zzzzz_0_xyyy_1[i] * wa_z[i];

        g_zzzzzz_0_xyyz_0[i] = 5.0 * g_zzzz_0_xyyz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xyyz_1[i] * fz_be_0 + g_zzzzz_0_xyy_1[i] * fi_acd_0 + g_zzzzz_0_xyyz_1[i] * wa_z[i];

        g_zzzzzz_0_xyzz_0[i] = 5.0 * g_zzzz_0_xyzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xyzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_0_xyz_1[i] * fi_acd_0 + g_zzzzz_0_xyzz_1[i] * wa_z[i];

        g_zzzzzz_0_xzzz_0[i] = 5.0 * g_zzzz_0_xzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xzzz_1[i] * fz_be_0 + 3.0 * g_zzzzz_0_xzz_1[i] * fi_acd_0 + g_zzzzz_0_xzzz_1[i] * wa_z[i];

        g_zzzzzz_0_yyyy_0[i] = 5.0 * g_zzzz_0_yyyy_0[i] * fbe_0 - 5.0 * g_zzzz_0_yyyy_1[i] * fz_be_0 + g_zzzzz_0_yyyy_1[i] * wa_z[i];

        g_zzzzzz_0_yyyz_0[i] = 5.0 * g_zzzz_0_yyyz_0[i] * fbe_0 - 5.0 * g_zzzz_0_yyyz_1[i] * fz_be_0 + g_zzzzz_0_yyy_1[i] * fi_acd_0 + g_zzzzz_0_yyyz_1[i] * wa_z[i];

        g_zzzzzz_0_yyzz_0[i] = 5.0 * g_zzzz_0_yyzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_yyzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_0_yyz_1[i] * fi_acd_0 + g_zzzzz_0_yyzz_1[i] * wa_z[i];

        g_zzzzzz_0_yzzz_0[i] = 5.0 * g_zzzz_0_yzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_yzzz_1[i] * fz_be_0 + 3.0 * g_zzzzz_0_yzz_1[i] * fi_acd_0 + g_zzzzz_0_yzzz_1[i] * wa_z[i];

        g_zzzzzz_0_zzzz_0[i] = 5.0 * g_zzzz_0_zzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_zzzz_1[i] * fz_be_0 + 4.0 * g_zzzzz_0_zzz_1[i] * fi_acd_0 + g_zzzzz_0_zzzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace


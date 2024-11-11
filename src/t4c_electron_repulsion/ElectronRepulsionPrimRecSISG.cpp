#include "ElectronRepulsionPrimRecSISG.hpp"

namespace erirec {  // erirec namespace

auto
comp_prim_electron_repulsion_sisg(CSimdArray<double>&   pbuffer,
                                  const size_t          idx_eri_0_sisg,
                                  size_t                idx_eri_0_sgsg,
                                  size_t                idx_eri_1_sgsg,
                                  size_t                idx_eri_1_shsf,
                                  size_t                idx_eri_0_shsg,
                                  size_t                idx_eri_1_shsg,
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

    /// Set up components of auxilary buffer : SGSG

    auto g_0_xxxx_0_xxxx_0 = pbuffer.data(idx_eri_0_sgsg);

    auto g_0_xxxx_0_xxxy_0 = pbuffer.data(idx_eri_0_sgsg + 1);

    auto g_0_xxxx_0_xxxz_0 = pbuffer.data(idx_eri_0_sgsg + 2);

    auto g_0_xxxx_0_xxyy_0 = pbuffer.data(idx_eri_0_sgsg + 3);

    auto g_0_xxxx_0_xxyz_0 = pbuffer.data(idx_eri_0_sgsg + 4);

    auto g_0_xxxx_0_xxzz_0 = pbuffer.data(idx_eri_0_sgsg + 5);

    auto g_0_xxxx_0_xyyy_0 = pbuffer.data(idx_eri_0_sgsg + 6);

    auto g_0_xxxx_0_xyyz_0 = pbuffer.data(idx_eri_0_sgsg + 7);

    auto g_0_xxxx_0_xyzz_0 = pbuffer.data(idx_eri_0_sgsg + 8);

    auto g_0_xxxx_0_xzzz_0 = pbuffer.data(idx_eri_0_sgsg + 9);

    auto g_0_xxxx_0_yyyy_0 = pbuffer.data(idx_eri_0_sgsg + 10);

    auto g_0_xxxx_0_yyyz_0 = pbuffer.data(idx_eri_0_sgsg + 11);

    auto g_0_xxxx_0_yyzz_0 = pbuffer.data(idx_eri_0_sgsg + 12);

    auto g_0_xxxx_0_yzzz_0 = pbuffer.data(idx_eri_0_sgsg + 13);

    auto g_0_xxxx_0_zzzz_0 = pbuffer.data(idx_eri_0_sgsg + 14);

    auto g_0_xxxy_0_xxxx_0 = pbuffer.data(idx_eri_0_sgsg + 15);

    auto g_0_xxxy_0_xxxz_0 = pbuffer.data(idx_eri_0_sgsg + 17);

    auto g_0_xxxy_0_xxzz_0 = pbuffer.data(idx_eri_0_sgsg + 20);

    auto g_0_xxxy_0_xzzz_0 = pbuffer.data(idx_eri_0_sgsg + 24);

    auto g_0_xxxz_0_xxxx_0 = pbuffer.data(idx_eri_0_sgsg + 30);

    auto g_0_xxxz_0_xxxy_0 = pbuffer.data(idx_eri_0_sgsg + 31);

    auto g_0_xxxz_0_xxyy_0 = pbuffer.data(idx_eri_0_sgsg + 33);

    auto g_0_xxxz_0_xyyy_0 = pbuffer.data(idx_eri_0_sgsg + 36);

    auto g_0_xxyy_0_xxxx_0 = pbuffer.data(idx_eri_0_sgsg + 45);

    auto g_0_xxyy_0_xxxy_0 = pbuffer.data(idx_eri_0_sgsg + 46);

    auto g_0_xxyy_0_xxxz_0 = pbuffer.data(idx_eri_0_sgsg + 47);

    auto g_0_xxyy_0_xxyy_0 = pbuffer.data(idx_eri_0_sgsg + 48);

    auto g_0_xxyy_0_xxyz_0 = pbuffer.data(idx_eri_0_sgsg + 49);

    auto g_0_xxyy_0_xxzz_0 = pbuffer.data(idx_eri_0_sgsg + 50);

    auto g_0_xxyy_0_xyyy_0 = pbuffer.data(idx_eri_0_sgsg + 51);

    auto g_0_xxyy_0_xyyz_0 = pbuffer.data(idx_eri_0_sgsg + 52);

    auto g_0_xxyy_0_xyzz_0 = pbuffer.data(idx_eri_0_sgsg + 53);

    auto g_0_xxyy_0_xzzz_0 = pbuffer.data(idx_eri_0_sgsg + 54);

    auto g_0_xxyy_0_yyyy_0 = pbuffer.data(idx_eri_0_sgsg + 55);

    auto g_0_xxyy_0_yyyz_0 = pbuffer.data(idx_eri_0_sgsg + 56);

    auto g_0_xxyy_0_yyzz_0 = pbuffer.data(idx_eri_0_sgsg + 57);

    auto g_0_xxyy_0_yzzz_0 = pbuffer.data(idx_eri_0_sgsg + 58);

    auto g_0_xxyy_0_zzzz_0 = pbuffer.data(idx_eri_0_sgsg + 59);

    auto g_0_xxzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sgsg + 75);

    auto g_0_xxzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sgsg + 76);

    auto g_0_xxzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sgsg + 77);

    auto g_0_xxzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sgsg + 78);

    auto g_0_xxzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sgsg + 79);

    auto g_0_xxzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sgsg + 80);

    auto g_0_xxzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sgsg + 81);

    auto g_0_xxzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sgsg + 82);

    auto g_0_xxzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sgsg + 83);

    auto g_0_xxzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sgsg + 84);

    auto g_0_xxzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sgsg + 85);

    auto g_0_xxzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sgsg + 86);

    auto g_0_xxzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sgsg + 87);

    auto g_0_xxzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sgsg + 88);

    auto g_0_xxzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sgsg + 89);

    auto g_0_xyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_sgsg + 91);

    auto g_0_xyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_sgsg + 93);

    auto g_0_xyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_sgsg + 94);

    auto g_0_xyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_sgsg + 96);

    auto g_0_xyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_sgsg + 97);

    auto g_0_xyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_sgsg + 98);

    auto g_0_xyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_sgsg + 100);

    auto g_0_xyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_sgsg + 101);

    auto g_0_xyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_sgsg + 102);

    auto g_0_xyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_sgsg + 103);

    auto g_0_xyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_sgsg + 104);

    auto g_0_xzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sgsg + 137);

    auto g_0_xzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sgsg + 139);

    auto g_0_xzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sgsg + 140);

    auto g_0_xzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sgsg + 142);

    auto g_0_xzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sgsg + 143);

    auto g_0_xzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sgsg + 144);

    auto g_0_xzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sgsg + 145);

    auto g_0_xzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sgsg + 146);

    auto g_0_xzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sgsg + 147);

    auto g_0_xzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sgsg + 148);

    auto g_0_xzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sgsg + 149);

    auto g_0_yyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_sgsg + 150);

    auto g_0_yyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_sgsg + 151);

    auto g_0_yyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_sgsg + 152);

    auto g_0_yyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_sgsg + 153);

    auto g_0_yyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_sgsg + 154);

    auto g_0_yyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_sgsg + 155);

    auto g_0_yyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_sgsg + 156);

    auto g_0_yyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_sgsg + 157);

    auto g_0_yyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_sgsg + 158);

    auto g_0_yyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_sgsg + 159);

    auto g_0_yyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_sgsg + 160);

    auto g_0_yyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_sgsg + 161);

    auto g_0_yyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_sgsg + 162);

    auto g_0_yyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_sgsg + 163);

    auto g_0_yyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_sgsg + 164);

    auto g_0_yyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_sgsg + 166);

    auto g_0_yyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_sgsg + 168);

    auto g_0_yyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_sgsg + 171);

    auto g_0_yyyz_0_yyyy_0 = pbuffer.data(idx_eri_0_sgsg + 175);

    auto g_0_yyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sgsg + 180);

    auto g_0_yyzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sgsg + 181);

    auto g_0_yyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sgsg + 182);

    auto g_0_yyzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sgsg + 183);

    auto g_0_yyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sgsg + 184);

    auto g_0_yyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sgsg + 185);

    auto g_0_yyzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sgsg + 186);

    auto g_0_yyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sgsg + 187);

    auto g_0_yyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sgsg + 188);

    auto g_0_yyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sgsg + 189);

    auto g_0_yyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sgsg + 190);

    auto g_0_yyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sgsg + 191);

    auto g_0_yyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sgsg + 192);

    auto g_0_yyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sgsg + 193);

    auto g_0_yyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sgsg + 194);

    auto g_0_yzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sgsg + 195);

    auto g_0_yzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sgsg + 197);

    auto g_0_yzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sgsg + 199);

    auto g_0_yzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sgsg + 200);

    auto g_0_yzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sgsg + 202);

    auto g_0_yzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sgsg + 203);

    auto g_0_yzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sgsg + 204);

    auto g_0_yzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sgsg + 206);

    auto g_0_yzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sgsg + 207);

    auto g_0_yzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sgsg + 208);

    auto g_0_yzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sgsg + 209);

    auto g_0_zzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sgsg + 210);

    auto g_0_zzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sgsg + 211);

    auto g_0_zzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sgsg + 212);

    auto g_0_zzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sgsg + 213);

    auto g_0_zzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sgsg + 214);

    auto g_0_zzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sgsg + 215);

    auto g_0_zzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sgsg + 216);

    auto g_0_zzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sgsg + 217);

    auto g_0_zzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sgsg + 218);

    auto g_0_zzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sgsg + 219);

    auto g_0_zzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sgsg + 220);

    auto g_0_zzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sgsg + 221);

    auto g_0_zzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sgsg + 222);

    auto g_0_zzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sgsg + 223);

    auto g_0_zzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sgsg + 224);

    /// Set up components of auxilary buffer : SGSG

    auto g_0_xxxx_0_xxxx_1 = pbuffer.data(idx_eri_1_sgsg);

    auto g_0_xxxx_0_xxxy_1 = pbuffer.data(idx_eri_1_sgsg + 1);

    auto g_0_xxxx_0_xxxz_1 = pbuffer.data(idx_eri_1_sgsg + 2);

    auto g_0_xxxx_0_xxyy_1 = pbuffer.data(idx_eri_1_sgsg + 3);

    auto g_0_xxxx_0_xxyz_1 = pbuffer.data(idx_eri_1_sgsg + 4);

    auto g_0_xxxx_0_xxzz_1 = pbuffer.data(idx_eri_1_sgsg + 5);

    auto g_0_xxxx_0_xyyy_1 = pbuffer.data(idx_eri_1_sgsg + 6);

    auto g_0_xxxx_0_xyyz_1 = pbuffer.data(idx_eri_1_sgsg + 7);

    auto g_0_xxxx_0_xyzz_1 = pbuffer.data(idx_eri_1_sgsg + 8);

    auto g_0_xxxx_0_xzzz_1 = pbuffer.data(idx_eri_1_sgsg + 9);

    auto g_0_xxxx_0_yyyy_1 = pbuffer.data(idx_eri_1_sgsg + 10);

    auto g_0_xxxx_0_yyyz_1 = pbuffer.data(idx_eri_1_sgsg + 11);

    auto g_0_xxxx_0_yyzz_1 = pbuffer.data(idx_eri_1_sgsg + 12);

    auto g_0_xxxx_0_yzzz_1 = pbuffer.data(idx_eri_1_sgsg + 13);

    auto g_0_xxxx_0_zzzz_1 = pbuffer.data(idx_eri_1_sgsg + 14);

    auto g_0_xxxy_0_xxxx_1 = pbuffer.data(idx_eri_1_sgsg + 15);

    auto g_0_xxxy_0_xxxz_1 = pbuffer.data(idx_eri_1_sgsg + 17);

    auto g_0_xxxy_0_xxzz_1 = pbuffer.data(idx_eri_1_sgsg + 20);

    auto g_0_xxxy_0_xzzz_1 = pbuffer.data(idx_eri_1_sgsg + 24);

    auto g_0_xxxz_0_xxxx_1 = pbuffer.data(idx_eri_1_sgsg + 30);

    auto g_0_xxxz_0_xxxy_1 = pbuffer.data(idx_eri_1_sgsg + 31);

    auto g_0_xxxz_0_xxyy_1 = pbuffer.data(idx_eri_1_sgsg + 33);

    auto g_0_xxxz_0_xyyy_1 = pbuffer.data(idx_eri_1_sgsg + 36);

    auto g_0_xxyy_0_xxxx_1 = pbuffer.data(idx_eri_1_sgsg + 45);

    auto g_0_xxyy_0_xxxy_1 = pbuffer.data(idx_eri_1_sgsg + 46);

    auto g_0_xxyy_0_xxxz_1 = pbuffer.data(idx_eri_1_sgsg + 47);

    auto g_0_xxyy_0_xxyy_1 = pbuffer.data(idx_eri_1_sgsg + 48);

    auto g_0_xxyy_0_xxyz_1 = pbuffer.data(idx_eri_1_sgsg + 49);

    auto g_0_xxyy_0_xxzz_1 = pbuffer.data(idx_eri_1_sgsg + 50);

    auto g_0_xxyy_0_xyyy_1 = pbuffer.data(idx_eri_1_sgsg + 51);

    auto g_0_xxyy_0_xyyz_1 = pbuffer.data(idx_eri_1_sgsg + 52);

    auto g_0_xxyy_0_xyzz_1 = pbuffer.data(idx_eri_1_sgsg + 53);

    auto g_0_xxyy_0_xzzz_1 = pbuffer.data(idx_eri_1_sgsg + 54);

    auto g_0_xxyy_0_yyyy_1 = pbuffer.data(idx_eri_1_sgsg + 55);

    auto g_0_xxyy_0_yyyz_1 = pbuffer.data(idx_eri_1_sgsg + 56);

    auto g_0_xxyy_0_yyzz_1 = pbuffer.data(idx_eri_1_sgsg + 57);

    auto g_0_xxyy_0_yzzz_1 = pbuffer.data(idx_eri_1_sgsg + 58);

    auto g_0_xxyy_0_zzzz_1 = pbuffer.data(idx_eri_1_sgsg + 59);

    auto g_0_xxzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sgsg + 75);

    auto g_0_xxzz_0_xxxy_1 = pbuffer.data(idx_eri_1_sgsg + 76);

    auto g_0_xxzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sgsg + 77);

    auto g_0_xxzz_0_xxyy_1 = pbuffer.data(idx_eri_1_sgsg + 78);

    auto g_0_xxzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sgsg + 79);

    auto g_0_xxzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sgsg + 80);

    auto g_0_xxzz_0_xyyy_1 = pbuffer.data(idx_eri_1_sgsg + 81);

    auto g_0_xxzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sgsg + 82);

    auto g_0_xxzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sgsg + 83);

    auto g_0_xxzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sgsg + 84);

    auto g_0_xxzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sgsg + 85);

    auto g_0_xxzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sgsg + 86);

    auto g_0_xxzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sgsg + 87);

    auto g_0_xxzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sgsg + 88);

    auto g_0_xxzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sgsg + 89);

    auto g_0_xyyy_0_xxxy_1 = pbuffer.data(idx_eri_1_sgsg + 91);

    auto g_0_xyyy_0_xxyy_1 = pbuffer.data(idx_eri_1_sgsg + 93);

    auto g_0_xyyy_0_xxyz_1 = pbuffer.data(idx_eri_1_sgsg + 94);

    auto g_0_xyyy_0_xyyy_1 = pbuffer.data(idx_eri_1_sgsg + 96);

    auto g_0_xyyy_0_xyyz_1 = pbuffer.data(idx_eri_1_sgsg + 97);

    auto g_0_xyyy_0_xyzz_1 = pbuffer.data(idx_eri_1_sgsg + 98);

    auto g_0_xyyy_0_yyyy_1 = pbuffer.data(idx_eri_1_sgsg + 100);

    auto g_0_xyyy_0_yyyz_1 = pbuffer.data(idx_eri_1_sgsg + 101);

    auto g_0_xyyy_0_yyzz_1 = pbuffer.data(idx_eri_1_sgsg + 102);

    auto g_0_xyyy_0_yzzz_1 = pbuffer.data(idx_eri_1_sgsg + 103);

    auto g_0_xyyy_0_zzzz_1 = pbuffer.data(idx_eri_1_sgsg + 104);

    auto g_0_xzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sgsg + 137);

    auto g_0_xzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sgsg + 139);

    auto g_0_xzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sgsg + 140);

    auto g_0_xzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sgsg + 142);

    auto g_0_xzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sgsg + 143);

    auto g_0_xzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sgsg + 144);

    auto g_0_xzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sgsg + 145);

    auto g_0_xzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sgsg + 146);

    auto g_0_xzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sgsg + 147);

    auto g_0_xzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sgsg + 148);

    auto g_0_xzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sgsg + 149);

    auto g_0_yyyy_0_xxxx_1 = pbuffer.data(idx_eri_1_sgsg + 150);

    auto g_0_yyyy_0_xxxy_1 = pbuffer.data(idx_eri_1_sgsg + 151);

    auto g_0_yyyy_0_xxxz_1 = pbuffer.data(idx_eri_1_sgsg + 152);

    auto g_0_yyyy_0_xxyy_1 = pbuffer.data(idx_eri_1_sgsg + 153);

    auto g_0_yyyy_0_xxyz_1 = pbuffer.data(idx_eri_1_sgsg + 154);

    auto g_0_yyyy_0_xxzz_1 = pbuffer.data(idx_eri_1_sgsg + 155);

    auto g_0_yyyy_0_xyyy_1 = pbuffer.data(idx_eri_1_sgsg + 156);

    auto g_0_yyyy_0_xyyz_1 = pbuffer.data(idx_eri_1_sgsg + 157);

    auto g_0_yyyy_0_xyzz_1 = pbuffer.data(idx_eri_1_sgsg + 158);

    auto g_0_yyyy_0_xzzz_1 = pbuffer.data(idx_eri_1_sgsg + 159);

    auto g_0_yyyy_0_yyyy_1 = pbuffer.data(idx_eri_1_sgsg + 160);

    auto g_0_yyyy_0_yyyz_1 = pbuffer.data(idx_eri_1_sgsg + 161);

    auto g_0_yyyy_0_yyzz_1 = pbuffer.data(idx_eri_1_sgsg + 162);

    auto g_0_yyyy_0_yzzz_1 = pbuffer.data(idx_eri_1_sgsg + 163);

    auto g_0_yyyy_0_zzzz_1 = pbuffer.data(idx_eri_1_sgsg + 164);

    auto g_0_yyyz_0_xxxy_1 = pbuffer.data(idx_eri_1_sgsg + 166);

    auto g_0_yyyz_0_xxyy_1 = pbuffer.data(idx_eri_1_sgsg + 168);

    auto g_0_yyyz_0_xyyy_1 = pbuffer.data(idx_eri_1_sgsg + 171);

    auto g_0_yyyz_0_yyyy_1 = pbuffer.data(idx_eri_1_sgsg + 175);

    auto g_0_yyzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sgsg + 180);

    auto g_0_yyzz_0_xxxy_1 = pbuffer.data(idx_eri_1_sgsg + 181);

    auto g_0_yyzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sgsg + 182);

    auto g_0_yyzz_0_xxyy_1 = pbuffer.data(idx_eri_1_sgsg + 183);

    auto g_0_yyzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sgsg + 184);

    auto g_0_yyzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sgsg + 185);

    auto g_0_yyzz_0_xyyy_1 = pbuffer.data(idx_eri_1_sgsg + 186);

    auto g_0_yyzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sgsg + 187);

    auto g_0_yyzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sgsg + 188);

    auto g_0_yyzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sgsg + 189);

    auto g_0_yyzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sgsg + 190);

    auto g_0_yyzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sgsg + 191);

    auto g_0_yyzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sgsg + 192);

    auto g_0_yyzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sgsg + 193);

    auto g_0_yyzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sgsg + 194);

    auto g_0_yzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sgsg + 195);

    auto g_0_yzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sgsg + 197);

    auto g_0_yzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sgsg + 199);

    auto g_0_yzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sgsg + 200);

    auto g_0_yzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sgsg + 202);

    auto g_0_yzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sgsg + 203);

    auto g_0_yzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sgsg + 204);

    auto g_0_yzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sgsg + 206);

    auto g_0_yzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sgsg + 207);

    auto g_0_yzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sgsg + 208);

    auto g_0_yzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sgsg + 209);

    auto g_0_zzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sgsg + 210);

    auto g_0_zzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_sgsg + 211);

    auto g_0_zzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sgsg + 212);

    auto g_0_zzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_sgsg + 213);

    auto g_0_zzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sgsg + 214);

    auto g_0_zzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sgsg + 215);

    auto g_0_zzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_sgsg + 216);

    auto g_0_zzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sgsg + 217);

    auto g_0_zzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sgsg + 218);

    auto g_0_zzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sgsg + 219);

    auto g_0_zzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sgsg + 220);

    auto g_0_zzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sgsg + 221);

    auto g_0_zzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sgsg + 222);

    auto g_0_zzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sgsg + 223);

    auto g_0_zzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sgsg + 224);

    /// Set up components of auxilary buffer : SHSF

    auto g_0_xxxxx_0_xxx_1 = pbuffer.data(idx_eri_1_shsf);

    auto g_0_xxxxx_0_xxy_1 = pbuffer.data(idx_eri_1_shsf + 1);

    auto g_0_xxxxx_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 2);

    auto g_0_xxxxx_0_xyy_1 = pbuffer.data(idx_eri_1_shsf + 3);

    auto g_0_xxxxx_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 4);

    auto g_0_xxxxx_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 5);

    auto g_0_xxxxx_0_yyy_1 = pbuffer.data(idx_eri_1_shsf + 6);

    auto g_0_xxxxx_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 7);

    auto g_0_xxxxx_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 8);

    auto g_0_xxxxx_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 9);

    auto g_0_xxxxz_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 22);

    auto g_0_xxxxz_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 24);

    auto g_0_xxxxz_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 25);

    auto g_0_xxxxz_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 27);

    auto g_0_xxxxz_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 28);

    auto g_0_xxxxz_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 29);

    auto g_0_xxxyy_0_xxx_1 = pbuffer.data(idx_eri_1_shsf + 30);

    auto g_0_xxxyy_0_xxy_1 = pbuffer.data(idx_eri_1_shsf + 31);

    auto g_0_xxxyy_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 32);

    auto g_0_xxxyy_0_xyy_1 = pbuffer.data(idx_eri_1_shsf + 33);

    auto g_0_xxxyy_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 34);

    auto g_0_xxxyy_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 35);

    auto g_0_xxxyy_0_yyy_1 = pbuffer.data(idx_eri_1_shsf + 36);

    auto g_0_xxxyy_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 37);

    auto g_0_xxxyy_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 38);

    auto g_0_xxxyy_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 39);

    auto g_0_xxxzz_0_xxx_1 = pbuffer.data(idx_eri_1_shsf + 50);

    auto g_0_xxxzz_0_xxy_1 = pbuffer.data(idx_eri_1_shsf + 51);

    auto g_0_xxxzz_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 52);

    auto g_0_xxxzz_0_xyy_1 = pbuffer.data(idx_eri_1_shsf + 53);

    auto g_0_xxxzz_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 54);

    auto g_0_xxxzz_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 55);

    auto g_0_xxxzz_0_yyy_1 = pbuffer.data(idx_eri_1_shsf + 56);

    auto g_0_xxxzz_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 57);

    auto g_0_xxxzz_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 58);

    auto g_0_xxxzz_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 59);

    auto g_0_xxyyy_0_xxx_1 = pbuffer.data(idx_eri_1_shsf + 60);

    auto g_0_xxyyy_0_xxy_1 = pbuffer.data(idx_eri_1_shsf + 61);

    auto g_0_xxyyy_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 62);

    auto g_0_xxyyy_0_xyy_1 = pbuffer.data(idx_eri_1_shsf + 63);

    auto g_0_xxyyy_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 64);

    auto g_0_xxyyy_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 65);

    auto g_0_xxyyy_0_yyy_1 = pbuffer.data(idx_eri_1_shsf + 66);

    auto g_0_xxyyy_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 67);

    auto g_0_xxyyy_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 68);

    auto g_0_xxyyy_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 69);

    auto g_0_xxzzz_0_xxx_1 = pbuffer.data(idx_eri_1_shsf + 90);

    auto g_0_xxzzz_0_xxy_1 = pbuffer.data(idx_eri_1_shsf + 91);

    auto g_0_xxzzz_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 92);

    auto g_0_xxzzz_0_xyy_1 = pbuffer.data(idx_eri_1_shsf + 93);

    auto g_0_xxzzz_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 94);

    auto g_0_xxzzz_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 95);

    auto g_0_xxzzz_0_yyy_1 = pbuffer.data(idx_eri_1_shsf + 96);

    auto g_0_xxzzz_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 97);

    auto g_0_xxzzz_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 98);

    auto g_0_xxzzz_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 99);

    auto g_0_xyyyy_0_xxy_1 = pbuffer.data(idx_eri_1_shsf + 101);

    auto g_0_xyyyy_0_xyy_1 = pbuffer.data(idx_eri_1_shsf + 103);

    auto g_0_xyyyy_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 104);

    auto g_0_xyyyy_0_yyy_1 = pbuffer.data(idx_eri_1_shsf + 106);

    auto g_0_xyyyy_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 107);

    auto g_0_xyyyy_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 108);

    auto g_0_xyyzz_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 124);

    auto g_0_xyyzz_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 127);

    auto g_0_xyyzz_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 128);

    auto g_0_xzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 142);

    auto g_0_xzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 144);

    auto g_0_xzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 145);

    auto g_0_xzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 147);

    auto g_0_xzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 148);

    auto g_0_xzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 149);

    auto g_0_yyyyy_0_xxx_1 = pbuffer.data(idx_eri_1_shsf + 150);

    auto g_0_yyyyy_0_xxy_1 = pbuffer.data(idx_eri_1_shsf + 151);

    auto g_0_yyyyy_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 152);

    auto g_0_yyyyy_0_xyy_1 = pbuffer.data(idx_eri_1_shsf + 153);

    auto g_0_yyyyy_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 154);

    auto g_0_yyyyy_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 155);

    auto g_0_yyyyy_0_yyy_1 = pbuffer.data(idx_eri_1_shsf + 156);

    auto g_0_yyyyy_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 157);

    auto g_0_yyyyy_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 158);

    auto g_0_yyyyy_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 159);

    auto g_0_yyyyz_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 162);

    auto g_0_yyyyz_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 164);

    auto g_0_yyyyz_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 165);

    auto g_0_yyyyz_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 167);

    auto g_0_yyyyz_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 168);

    auto g_0_yyyyz_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 169);

    auto g_0_yyyzz_0_xxx_1 = pbuffer.data(idx_eri_1_shsf + 170);

    auto g_0_yyyzz_0_xxy_1 = pbuffer.data(idx_eri_1_shsf + 171);

    auto g_0_yyyzz_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 172);

    auto g_0_yyyzz_0_xyy_1 = pbuffer.data(idx_eri_1_shsf + 173);

    auto g_0_yyyzz_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 174);

    auto g_0_yyyzz_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 175);

    auto g_0_yyyzz_0_yyy_1 = pbuffer.data(idx_eri_1_shsf + 176);

    auto g_0_yyyzz_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 177);

    auto g_0_yyyzz_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 178);

    auto g_0_yyyzz_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 179);

    auto g_0_yyzzz_0_xxx_1 = pbuffer.data(idx_eri_1_shsf + 180);

    auto g_0_yyzzz_0_xxy_1 = pbuffer.data(idx_eri_1_shsf + 181);

    auto g_0_yyzzz_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 182);

    auto g_0_yyzzz_0_xyy_1 = pbuffer.data(idx_eri_1_shsf + 183);

    auto g_0_yyzzz_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 184);

    auto g_0_yyzzz_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 185);

    auto g_0_yyzzz_0_yyy_1 = pbuffer.data(idx_eri_1_shsf + 186);

    auto g_0_yyzzz_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 187);

    auto g_0_yyzzz_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 188);

    auto g_0_yyzzz_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 189);

    auto g_0_yzzzz_0_xxy_1 = pbuffer.data(idx_eri_1_shsf + 191);

    auto g_0_yzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 192);

    auto g_0_yzzzz_0_xyy_1 = pbuffer.data(idx_eri_1_shsf + 193);

    auto g_0_yzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 194);

    auto g_0_yzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 195);

    auto g_0_yzzzz_0_yyy_1 = pbuffer.data(idx_eri_1_shsf + 196);

    auto g_0_yzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 197);

    auto g_0_yzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 198);

    auto g_0_yzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 199);

    auto g_0_zzzzz_0_xxx_1 = pbuffer.data(idx_eri_1_shsf + 200);

    auto g_0_zzzzz_0_xxy_1 = pbuffer.data(idx_eri_1_shsf + 201);

    auto g_0_zzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 202);

    auto g_0_zzzzz_0_xyy_1 = pbuffer.data(idx_eri_1_shsf + 203);

    auto g_0_zzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 204);

    auto g_0_zzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 205);

    auto g_0_zzzzz_0_yyy_1 = pbuffer.data(idx_eri_1_shsf + 206);

    auto g_0_zzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 207);

    auto g_0_zzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 208);

    auto g_0_zzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 209);

    /// Set up components of auxilary buffer : SHSG

    auto g_0_xxxxx_0_xxxx_0 = pbuffer.data(idx_eri_0_shsg);

    auto g_0_xxxxx_0_xxxy_0 = pbuffer.data(idx_eri_0_shsg + 1);

    auto g_0_xxxxx_0_xxxz_0 = pbuffer.data(idx_eri_0_shsg + 2);

    auto g_0_xxxxx_0_xxyy_0 = pbuffer.data(idx_eri_0_shsg + 3);

    auto g_0_xxxxx_0_xxyz_0 = pbuffer.data(idx_eri_0_shsg + 4);

    auto g_0_xxxxx_0_xxzz_0 = pbuffer.data(idx_eri_0_shsg + 5);

    auto g_0_xxxxx_0_xyyy_0 = pbuffer.data(idx_eri_0_shsg + 6);

    auto g_0_xxxxx_0_xyyz_0 = pbuffer.data(idx_eri_0_shsg + 7);

    auto g_0_xxxxx_0_xyzz_0 = pbuffer.data(idx_eri_0_shsg + 8);

    auto g_0_xxxxx_0_xzzz_0 = pbuffer.data(idx_eri_0_shsg + 9);

    auto g_0_xxxxx_0_yyyy_0 = pbuffer.data(idx_eri_0_shsg + 10);

    auto g_0_xxxxx_0_yyyz_0 = pbuffer.data(idx_eri_0_shsg + 11);

    auto g_0_xxxxx_0_yyzz_0 = pbuffer.data(idx_eri_0_shsg + 12);

    auto g_0_xxxxx_0_yzzz_0 = pbuffer.data(idx_eri_0_shsg + 13);

    auto g_0_xxxxx_0_zzzz_0 = pbuffer.data(idx_eri_0_shsg + 14);

    auto g_0_xxxxy_0_xxxx_0 = pbuffer.data(idx_eri_0_shsg + 15);

    auto g_0_xxxxy_0_xxxy_0 = pbuffer.data(idx_eri_0_shsg + 16);

    auto g_0_xxxxy_0_xxxz_0 = pbuffer.data(idx_eri_0_shsg + 17);

    auto g_0_xxxxy_0_xxyy_0 = pbuffer.data(idx_eri_0_shsg + 18);

    auto g_0_xxxxy_0_xxzz_0 = pbuffer.data(idx_eri_0_shsg + 20);

    auto g_0_xxxxy_0_xyyy_0 = pbuffer.data(idx_eri_0_shsg + 21);

    auto g_0_xxxxy_0_xzzz_0 = pbuffer.data(idx_eri_0_shsg + 24);

    auto g_0_xxxxy_0_yyyy_0 = pbuffer.data(idx_eri_0_shsg + 25);

    auto g_0_xxxxz_0_xxxx_0 = pbuffer.data(idx_eri_0_shsg + 30);

    auto g_0_xxxxz_0_xxxy_0 = pbuffer.data(idx_eri_0_shsg + 31);

    auto g_0_xxxxz_0_xxxz_0 = pbuffer.data(idx_eri_0_shsg + 32);

    auto g_0_xxxxz_0_xxyy_0 = pbuffer.data(idx_eri_0_shsg + 33);

    auto g_0_xxxxz_0_xxyz_0 = pbuffer.data(idx_eri_0_shsg + 34);

    auto g_0_xxxxz_0_xxzz_0 = pbuffer.data(idx_eri_0_shsg + 35);

    auto g_0_xxxxz_0_xyyy_0 = pbuffer.data(idx_eri_0_shsg + 36);

    auto g_0_xxxxz_0_xyyz_0 = pbuffer.data(idx_eri_0_shsg + 37);

    auto g_0_xxxxz_0_xyzz_0 = pbuffer.data(idx_eri_0_shsg + 38);

    auto g_0_xxxxz_0_xzzz_0 = pbuffer.data(idx_eri_0_shsg + 39);

    auto g_0_xxxxz_0_yyyz_0 = pbuffer.data(idx_eri_0_shsg + 41);

    auto g_0_xxxxz_0_yyzz_0 = pbuffer.data(idx_eri_0_shsg + 42);

    auto g_0_xxxxz_0_yzzz_0 = pbuffer.data(idx_eri_0_shsg + 43);

    auto g_0_xxxxz_0_zzzz_0 = pbuffer.data(idx_eri_0_shsg + 44);

    auto g_0_xxxyy_0_xxxx_0 = pbuffer.data(idx_eri_0_shsg + 45);

    auto g_0_xxxyy_0_xxxy_0 = pbuffer.data(idx_eri_0_shsg + 46);

    auto g_0_xxxyy_0_xxxz_0 = pbuffer.data(idx_eri_0_shsg + 47);

    auto g_0_xxxyy_0_xxyy_0 = pbuffer.data(idx_eri_0_shsg + 48);

    auto g_0_xxxyy_0_xxyz_0 = pbuffer.data(idx_eri_0_shsg + 49);

    auto g_0_xxxyy_0_xxzz_0 = pbuffer.data(idx_eri_0_shsg + 50);

    auto g_0_xxxyy_0_xyyy_0 = pbuffer.data(idx_eri_0_shsg + 51);

    auto g_0_xxxyy_0_xyyz_0 = pbuffer.data(idx_eri_0_shsg + 52);

    auto g_0_xxxyy_0_xyzz_0 = pbuffer.data(idx_eri_0_shsg + 53);

    auto g_0_xxxyy_0_xzzz_0 = pbuffer.data(idx_eri_0_shsg + 54);

    auto g_0_xxxyy_0_yyyy_0 = pbuffer.data(idx_eri_0_shsg + 55);

    auto g_0_xxxyy_0_yyyz_0 = pbuffer.data(idx_eri_0_shsg + 56);

    auto g_0_xxxyy_0_yyzz_0 = pbuffer.data(idx_eri_0_shsg + 57);

    auto g_0_xxxyy_0_yzzz_0 = pbuffer.data(idx_eri_0_shsg + 58);

    auto g_0_xxxyy_0_zzzz_0 = pbuffer.data(idx_eri_0_shsg + 59);

    auto g_0_xxxzz_0_xxxx_0 = pbuffer.data(idx_eri_0_shsg + 75);

    auto g_0_xxxzz_0_xxxy_0 = pbuffer.data(idx_eri_0_shsg + 76);

    auto g_0_xxxzz_0_xxxz_0 = pbuffer.data(idx_eri_0_shsg + 77);

    auto g_0_xxxzz_0_xxyy_0 = pbuffer.data(idx_eri_0_shsg + 78);

    auto g_0_xxxzz_0_xxyz_0 = pbuffer.data(idx_eri_0_shsg + 79);

    auto g_0_xxxzz_0_xxzz_0 = pbuffer.data(idx_eri_0_shsg + 80);

    auto g_0_xxxzz_0_xyyy_0 = pbuffer.data(idx_eri_0_shsg + 81);

    auto g_0_xxxzz_0_xyyz_0 = pbuffer.data(idx_eri_0_shsg + 82);

    auto g_0_xxxzz_0_xyzz_0 = pbuffer.data(idx_eri_0_shsg + 83);

    auto g_0_xxxzz_0_xzzz_0 = pbuffer.data(idx_eri_0_shsg + 84);

    auto g_0_xxxzz_0_yyyy_0 = pbuffer.data(idx_eri_0_shsg + 85);

    auto g_0_xxxzz_0_yyyz_0 = pbuffer.data(idx_eri_0_shsg + 86);

    auto g_0_xxxzz_0_yyzz_0 = pbuffer.data(idx_eri_0_shsg + 87);

    auto g_0_xxxzz_0_yzzz_0 = pbuffer.data(idx_eri_0_shsg + 88);

    auto g_0_xxxzz_0_zzzz_0 = pbuffer.data(idx_eri_0_shsg + 89);

    auto g_0_xxyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_shsg + 90);

    auto g_0_xxyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_shsg + 91);

    auto g_0_xxyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_shsg + 92);

    auto g_0_xxyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_shsg + 93);

    auto g_0_xxyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_shsg + 94);

    auto g_0_xxyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_shsg + 95);

    auto g_0_xxyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_shsg + 96);

    auto g_0_xxyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_shsg + 97);

    auto g_0_xxyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_shsg + 98);

    auto g_0_xxyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_shsg + 99);

    auto g_0_xxyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_shsg + 100);

    auto g_0_xxyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_shsg + 101);

    auto g_0_xxyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_shsg + 102);

    auto g_0_xxyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_shsg + 103);

    auto g_0_xxyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_shsg + 104);

    auto g_0_xxyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_shsg + 106);

    auto g_0_xxyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_shsg + 108);

    auto g_0_xxyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_shsg + 111);

    auto g_0_xxyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_shsg + 120);

    auto g_0_xxyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_shsg + 122);

    auto g_0_xxyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_shsg + 125);

    auto g_0_xxyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_shsg + 129);

    auto g_0_xxzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_shsg + 135);

    auto g_0_xxzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_shsg + 136);

    auto g_0_xxzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_shsg + 137);

    auto g_0_xxzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_shsg + 138);

    auto g_0_xxzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_shsg + 139);

    auto g_0_xxzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_shsg + 140);

    auto g_0_xxzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_shsg + 141);

    auto g_0_xxzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_shsg + 142);

    auto g_0_xxzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_shsg + 143);

    auto g_0_xxzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_shsg + 144);

    auto g_0_xxzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_shsg + 145);

    auto g_0_xxzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_shsg + 146);

    auto g_0_xxzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_shsg + 147);

    auto g_0_xxzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_shsg + 148);

    auto g_0_xxzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_shsg + 149);

    auto g_0_xyyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_shsg + 150);

    auto g_0_xyyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_shsg + 151);

    auto g_0_xyyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_shsg + 153);

    auto g_0_xyyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_shsg + 154);

    auto g_0_xyyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_shsg + 156);

    auto g_0_xyyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_shsg + 157);

    auto g_0_xyyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_shsg + 158);

    auto g_0_xyyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_shsg + 160);

    auto g_0_xyyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_shsg + 161);

    auto g_0_xyyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_shsg + 162);

    auto g_0_xyyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_shsg + 163);

    auto g_0_xyyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_shsg + 164);

    auto g_0_xyyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_shsg + 184);

    auto g_0_xyyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_shsg + 187);

    auto g_0_xyyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_shsg + 188);

    auto g_0_xyyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_shsg + 190);

    auto g_0_xyyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_shsg + 191);

    auto g_0_xyyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_shsg + 192);

    auto g_0_xyyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_shsg + 193);

    auto g_0_xyyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_shsg + 194);

    auto g_0_xzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_shsg + 210);

    auto g_0_xzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_shsg + 212);

    auto g_0_xzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_shsg + 214);

    auto g_0_xzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_shsg + 215);

    auto g_0_xzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_shsg + 217);

    auto g_0_xzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_shsg + 218);

    auto g_0_xzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_shsg + 219);

    auto g_0_xzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_shsg + 220);

    auto g_0_xzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_shsg + 221);

    auto g_0_xzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_shsg + 222);

    auto g_0_xzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_shsg + 223);

    auto g_0_xzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_shsg + 224);

    auto g_0_yyyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_shsg + 225);

    auto g_0_yyyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_shsg + 226);

    auto g_0_yyyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_shsg + 227);

    auto g_0_yyyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_shsg + 228);

    auto g_0_yyyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_shsg + 229);

    auto g_0_yyyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_shsg + 230);

    auto g_0_yyyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_shsg + 231);

    auto g_0_yyyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_shsg + 232);

    auto g_0_yyyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_shsg + 233);

    auto g_0_yyyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_shsg + 234);

    auto g_0_yyyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_shsg + 235);

    auto g_0_yyyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_shsg + 236);

    auto g_0_yyyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_shsg + 237);

    auto g_0_yyyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_shsg + 238);

    auto g_0_yyyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_shsg + 239);

    auto g_0_yyyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_shsg + 241);

    auto g_0_yyyyz_0_xxxz_0 = pbuffer.data(idx_eri_0_shsg + 242);

    auto g_0_yyyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_shsg + 243);

    auto g_0_yyyyz_0_xxyz_0 = pbuffer.data(idx_eri_0_shsg + 244);

    auto g_0_yyyyz_0_xxzz_0 = pbuffer.data(idx_eri_0_shsg + 245);

    auto g_0_yyyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_shsg + 246);

    auto g_0_yyyyz_0_xyyz_0 = pbuffer.data(idx_eri_0_shsg + 247);

    auto g_0_yyyyz_0_xyzz_0 = pbuffer.data(idx_eri_0_shsg + 248);

    auto g_0_yyyyz_0_xzzz_0 = pbuffer.data(idx_eri_0_shsg + 249);

    auto g_0_yyyyz_0_yyyy_0 = pbuffer.data(idx_eri_0_shsg + 250);

    auto g_0_yyyyz_0_yyyz_0 = pbuffer.data(idx_eri_0_shsg + 251);

    auto g_0_yyyyz_0_yyzz_0 = pbuffer.data(idx_eri_0_shsg + 252);

    auto g_0_yyyyz_0_yzzz_0 = pbuffer.data(idx_eri_0_shsg + 253);

    auto g_0_yyyyz_0_zzzz_0 = pbuffer.data(idx_eri_0_shsg + 254);

    auto g_0_yyyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_shsg + 255);

    auto g_0_yyyzz_0_xxxy_0 = pbuffer.data(idx_eri_0_shsg + 256);

    auto g_0_yyyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_shsg + 257);

    auto g_0_yyyzz_0_xxyy_0 = pbuffer.data(idx_eri_0_shsg + 258);

    auto g_0_yyyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_shsg + 259);

    auto g_0_yyyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_shsg + 260);

    auto g_0_yyyzz_0_xyyy_0 = pbuffer.data(idx_eri_0_shsg + 261);

    auto g_0_yyyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_shsg + 262);

    auto g_0_yyyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_shsg + 263);

    auto g_0_yyyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_shsg + 264);

    auto g_0_yyyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_shsg + 265);

    auto g_0_yyyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_shsg + 266);

    auto g_0_yyyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_shsg + 267);

    auto g_0_yyyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_shsg + 268);

    auto g_0_yyyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_shsg + 269);

    auto g_0_yyzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_shsg + 270);

    auto g_0_yyzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_shsg + 271);

    auto g_0_yyzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_shsg + 272);

    auto g_0_yyzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_shsg + 273);

    auto g_0_yyzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_shsg + 274);

    auto g_0_yyzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_shsg + 275);

    auto g_0_yyzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_shsg + 276);

    auto g_0_yyzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_shsg + 277);

    auto g_0_yyzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_shsg + 278);

    auto g_0_yyzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_shsg + 279);

    auto g_0_yyzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_shsg + 280);

    auto g_0_yyzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_shsg + 281);

    auto g_0_yyzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_shsg + 282);

    auto g_0_yyzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_shsg + 283);

    auto g_0_yyzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_shsg + 284);

    auto g_0_yzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_shsg + 285);

    auto g_0_yzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_shsg + 286);

    auto g_0_yzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_shsg + 287);

    auto g_0_yzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_shsg + 288);

    auto g_0_yzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_shsg + 289);

    auto g_0_yzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_shsg + 290);

    auto g_0_yzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_shsg + 291);

    auto g_0_yzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_shsg + 292);

    auto g_0_yzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_shsg + 293);

    auto g_0_yzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_shsg + 294);

    auto g_0_yzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_shsg + 295);

    auto g_0_yzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_shsg + 296);

    auto g_0_yzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_shsg + 297);

    auto g_0_yzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_shsg + 298);

    auto g_0_yzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_shsg + 299);

    auto g_0_zzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_shsg + 300);

    auto g_0_zzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_shsg + 301);

    auto g_0_zzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_shsg + 302);

    auto g_0_zzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_shsg + 303);

    auto g_0_zzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_shsg + 304);

    auto g_0_zzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_shsg + 305);

    auto g_0_zzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_shsg + 306);

    auto g_0_zzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_shsg + 307);

    auto g_0_zzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_shsg + 308);

    auto g_0_zzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_shsg + 309);

    auto g_0_zzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_shsg + 310);

    auto g_0_zzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_shsg + 311);

    auto g_0_zzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_shsg + 312);

    auto g_0_zzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_shsg + 313);

    auto g_0_zzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_shsg + 314);

    /// Set up components of auxilary buffer : SHSG

    auto g_0_xxxxx_0_xxxx_1 = pbuffer.data(idx_eri_1_shsg);

    auto g_0_xxxxx_0_xxxy_1 = pbuffer.data(idx_eri_1_shsg + 1);

    auto g_0_xxxxx_0_xxxz_1 = pbuffer.data(idx_eri_1_shsg + 2);

    auto g_0_xxxxx_0_xxyy_1 = pbuffer.data(idx_eri_1_shsg + 3);

    auto g_0_xxxxx_0_xxyz_1 = pbuffer.data(idx_eri_1_shsg + 4);

    auto g_0_xxxxx_0_xxzz_1 = pbuffer.data(idx_eri_1_shsg + 5);

    auto g_0_xxxxx_0_xyyy_1 = pbuffer.data(idx_eri_1_shsg + 6);

    auto g_0_xxxxx_0_xyyz_1 = pbuffer.data(idx_eri_1_shsg + 7);

    auto g_0_xxxxx_0_xyzz_1 = pbuffer.data(idx_eri_1_shsg + 8);

    auto g_0_xxxxx_0_xzzz_1 = pbuffer.data(idx_eri_1_shsg + 9);

    auto g_0_xxxxx_0_yyyy_1 = pbuffer.data(idx_eri_1_shsg + 10);

    auto g_0_xxxxx_0_yyyz_1 = pbuffer.data(idx_eri_1_shsg + 11);

    auto g_0_xxxxx_0_yyzz_1 = pbuffer.data(idx_eri_1_shsg + 12);

    auto g_0_xxxxx_0_yzzz_1 = pbuffer.data(idx_eri_1_shsg + 13);

    auto g_0_xxxxx_0_zzzz_1 = pbuffer.data(idx_eri_1_shsg + 14);

    auto g_0_xxxxy_0_xxxx_1 = pbuffer.data(idx_eri_1_shsg + 15);

    auto g_0_xxxxy_0_xxxy_1 = pbuffer.data(idx_eri_1_shsg + 16);

    auto g_0_xxxxy_0_xxxz_1 = pbuffer.data(idx_eri_1_shsg + 17);

    auto g_0_xxxxy_0_xxyy_1 = pbuffer.data(idx_eri_1_shsg + 18);

    auto g_0_xxxxy_0_xxzz_1 = pbuffer.data(idx_eri_1_shsg + 20);

    auto g_0_xxxxy_0_xyyy_1 = pbuffer.data(idx_eri_1_shsg + 21);

    auto g_0_xxxxy_0_xzzz_1 = pbuffer.data(idx_eri_1_shsg + 24);

    auto g_0_xxxxy_0_yyyy_1 = pbuffer.data(idx_eri_1_shsg + 25);

    auto g_0_xxxxz_0_xxxx_1 = pbuffer.data(idx_eri_1_shsg + 30);

    auto g_0_xxxxz_0_xxxy_1 = pbuffer.data(idx_eri_1_shsg + 31);

    auto g_0_xxxxz_0_xxxz_1 = pbuffer.data(idx_eri_1_shsg + 32);

    auto g_0_xxxxz_0_xxyy_1 = pbuffer.data(idx_eri_1_shsg + 33);

    auto g_0_xxxxz_0_xxyz_1 = pbuffer.data(idx_eri_1_shsg + 34);

    auto g_0_xxxxz_0_xxzz_1 = pbuffer.data(idx_eri_1_shsg + 35);

    auto g_0_xxxxz_0_xyyy_1 = pbuffer.data(idx_eri_1_shsg + 36);

    auto g_0_xxxxz_0_xyyz_1 = pbuffer.data(idx_eri_1_shsg + 37);

    auto g_0_xxxxz_0_xyzz_1 = pbuffer.data(idx_eri_1_shsg + 38);

    auto g_0_xxxxz_0_xzzz_1 = pbuffer.data(idx_eri_1_shsg + 39);

    auto g_0_xxxxz_0_yyyz_1 = pbuffer.data(idx_eri_1_shsg + 41);

    auto g_0_xxxxz_0_yyzz_1 = pbuffer.data(idx_eri_1_shsg + 42);

    auto g_0_xxxxz_0_yzzz_1 = pbuffer.data(idx_eri_1_shsg + 43);

    auto g_0_xxxxz_0_zzzz_1 = pbuffer.data(idx_eri_1_shsg + 44);

    auto g_0_xxxyy_0_xxxx_1 = pbuffer.data(idx_eri_1_shsg + 45);

    auto g_0_xxxyy_0_xxxy_1 = pbuffer.data(idx_eri_1_shsg + 46);

    auto g_0_xxxyy_0_xxxz_1 = pbuffer.data(idx_eri_1_shsg + 47);

    auto g_0_xxxyy_0_xxyy_1 = pbuffer.data(idx_eri_1_shsg + 48);

    auto g_0_xxxyy_0_xxyz_1 = pbuffer.data(idx_eri_1_shsg + 49);

    auto g_0_xxxyy_0_xxzz_1 = pbuffer.data(idx_eri_1_shsg + 50);

    auto g_0_xxxyy_0_xyyy_1 = pbuffer.data(idx_eri_1_shsg + 51);

    auto g_0_xxxyy_0_xyyz_1 = pbuffer.data(idx_eri_1_shsg + 52);

    auto g_0_xxxyy_0_xyzz_1 = pbuffer.data(idx_eri_1_shsg + 53);

    auto g_0_xxxyy_0_xzzz_1 = pbuffer.data(idx_eri_1_shsg + 54);

    auto g_0_xxxyy_0_yyyy_1 = pbuffer.data(idx_eri_1_shsg + 55);

    auto g_0_xxxyy_0_yyyz_1 = pbuffer.data(idx_eri_1_shsg + 56);

    auto g_0_xxxyy_0_yyzz_1 = pbuffer.data(idx_eri_1_shsg + 57);

    auto g_0_xxxyy_0_yzzz_1 = pbuffer.data(idx_eri_1_shsg + 58);

    auto g_0_xxxyy_0_zzzz_1 = pbuffer.data(idx_eri_1_shsg + 59);

    auto g_0_xxxzz_0_xxxx_1 = pbuffer.data(idx_eri_1_shsg + 75);

    auto g_0_xxxzz_0_xxxy_1 = pbuffer.data(idx_eri_1_shsg + 76);

    auto g_0_xxxzz_0_xxxz_1 = pbuffer.data(idx_eri_1_shsg + 77);

    auto g_0_xxxzz_0_xxyy_1 = pbuffer.data(idx_eri_1_shsg + 78);

    auto g_0_xxxzz_0_xxyz_1 = pbuffer.data(idx_eri_1_shsg + 79);

    auto g_0_xxxzz_0_xxzz_1 = pbuffer.data(idx_eri_1_shsg + 80);

    auto g_0_xxxzz_0_xyyy_1 = pbuffer.data(idx_eri_1_shsg + 81);

    auto g_0_xxxzz_0_xyyz_1 = pbuffer.data(idx_eri_1_shsg + 82);

    auto g_0_xxxzz_0_xyzz_1 = pbuffer.data(idx_eri_1_shsg + 83);

    auto g_0_xxxzz_0_xzzz_1 = pbuffer.data(idx_eri_1_shsg + 84);

    auto g_0_xxxzz_0_yyyy_1 = pbuffer.data(idx_eri_1_shsg + 85);

    auto g_0_xxxzz_0_yyyz_1 = pbuffer.data(idx_eri_1_shsg + 86);

    auto g_0_xxxzz_0_yyzz_1 = pbuffer.data(idx_eri_1_shsg + 87);

    auto g_0_xxxzz_0_yzzz_1 = pbuffer.data(idx_eri_1_shsg + 88);

    auto g_0_xxxzz_0_zzzz_1 = pbuffer.data(idx_eri_1_shsg + 89);

    auto g_0_xxyyy_0_xxxx_1 = pbuffer.data(idx_eri_1_shsg + 90);

    auto g_0_xxyyy_0_xxxy_1 = pbuffer.data(idx_eri_1_shsg + 91);

    auto g_0_xxyyy_0_xxxz_1 = pbuffer.data(idx_eri_1_shsg + 92);

    auto g_0_xxyyy_0_xxyy_1 = pbuffer.data(idx_eri_1_shsg + 93);

    auto g_0_xxyyy_0_xxyz_1 = pbuffer.data(idx_eri_1_shsg + 94);

    auto g_0_xxyyy_0_xxzz_1 = pbuffer.data(idx_eri_1_shsg + 95);

    auto g_0_xxyyy_0_xyyy_1 = pbuffer.data(idx_eri_1_shsg + 96);

    auto g_0_xxyyy_0_xyyz_1 = pbuffer.data(idx_eri_1_shsg + 97);

    auto g_0_xxyyy_0_xyzz_1 = pbuffer.data(idx_eri_1_shsg + 98);

    auto g_0_xxyyy_0_xzzz_1 = pbuffer.data(idx_eri_1_shsg + 99);

    auto g_0_xxyyy_0_yyyy_1 = pbuffer.data(idx_eri_1_shsg + 100);

    auto g_0_xxyyy_0_yyyz_1 = pbuffer.data(idx_eri_1_shsg + 101);

    auto g_0_xxyyy_0_yyzz_1 = pbuffer.data(idx_eri_1_shsg + 102);

    auto g_0_xxyyy_0_yzzz_1 = pbuffer.data(idx_eri_1_shsg + 103);

    auto g_0_xxyyy_0_zzzz_1 = pbuffer.data(idx_eri_1_shsg + 104);

    auto g_0_xxyyz_0_xxxy_1 = pbuffer.data(idx_eri_1_shsg + 106);

    auto g_0_xxyyz_0_xxyy_1 = pbuffer.data(idx_eri_1_shsg + 108);

    auto g_0_xxyyz_0_xyyy_1 = pbuffer.data(idx_eri_1_shsg + 111);

    auto g_0_xxyzz_0_xxxx_1 = pbuffer.data(idx_eri_1_shsg + 120);

    auto g_0_xxyzz_0_xxxz_1 = pbuffer.data(idx_eri_1_shsg + 122);

    auto g_0_xxyzz_0_xxzz_1 = pbuffer.data(idx_eri_1_shsg + 125);

    auto g_0_xxyzz_0_xzzz_1 = pbuffer.data(idx_eri_1_shsg + 129);

    auto g_0_xxzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_shsg + 135);

    auto g_0_xxzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_shsg + 136);

    auto g_0_xxzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_shsg + 137);

    auto g_0_xxzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_shsg + 138);

    auto g_0_xxzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_shsg + 139);

    auto g_0_xxzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_shsg + 140);

    auto g_0_xxzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_shsg + 141);

    auto g_0_xxzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_shsg + 142);

    auto g_0_xxzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_shsg + 143);

    auto g_0_xxzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_shsg + 144);

    auto g_0_xxzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_shsg + 145);

    auto g_0_xxzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_shsg + 146);

    auto g_0_xxzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_shsg + 147);

    auto g_0_xxzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_shsg + 148);

    auto g_0_xxzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_shsg + 149);

    auto g_0_xyyyy_0_xxxx_1 = pbuffer.data(idx_eri_1_shsg + 150);

    auto g_0_xyyyy_0_xxxy_1 = pbuffer.data(idx_eri_1_shsg + 151);

    auto g_0_xyyyy_0_xxyy_1 = pbuffer.data(idx_eri_1_shsg + 153);

    auto g_0_xyyyy_0_xxyz_1 = pbuffer.data(idx_eri_1_shsg + 154);

    auto g_0_xyyyy_0_xyyy_1 = pbuffer.data(idx_eri_1_shsg + 156);

    auto g_0_xyyyy_0_xyyz_1 = pbuffer.data(idx_eri_1_shsg + 157);

    auto g_0_xyyyy_0_xyzz_1 = pbuffer.data(idx_eri_1_shsg + 158);

    auto g_0_xyyyy_0_yyyy_1 = pbuffer.data(idx_eri_1_shsg + 160);

    auto g_0_xyyyy_0_yyyz_1 = pbuffer.data(idx_eri_1_shsg + 161);

    auto g_0_xyyyy_0_yyzz_1 = pbuffer.data(idx_eri_1_shsg + 162);

    auto g_0_xyyyy_0_yzzz_1 = pbuffer.data(idx_eri_1_shsg + 163);

    auto g_0_xyyyy_0_zzzz_1 = pbuffer.data(idx_eri_1_shsg + 164);

    auto g_0_xyyzz_0_xxyz_1 = pbuffer.data(idx_eri_1_shsg + 184);

    auto g_0_xyyzz_0_xyyz_1 = pbuffer.data(idx_eri_1_shsg + 187);

    auto g_0_xyyzz_0_xyzz_1 = pbuffer.data(idx_eri_1_shsg + 188);

    auto g_0_xyyzz_0_yyyy_1 = pbuffer.data(idx_eri_1_shsg + 190);

    auto g_0_xyyzz_0_yyyz_1 = pbuffer.data(idx_eri_1_shsg + 191);

    auto g_0_xyyzz_0_yyzz_1 = pbuffer.data(idx_eri_1_shsg + 192);

    auto g_0_xyyzz_0_yzzz_1 = pbuffer.data(idx_eri_1_shsg + 193);

    auto g_0_xyyzz_0_zzzz_1 = pbuffer.data(idx_eri_1_shsg + 194);

    auto g_0_xzzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_shsg + 210);

    auto g_0_xzzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_shsg + 212);

    auto g_0_xzzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_shsg + 214);

    auto g_0_xzzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_shsg + 215);

    auto g_0_xzzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_shsg + 217);

    auto g_0_xzzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_shsg + 218);

    auto g_0_xzzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_shsg + 219);

    auto g_0_xzzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_shsg + 220);

    auto g_0_xzzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_shsg + 221);

    auto g_0_xzzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_shsg + 222);

    auto g_0_xzzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_shsg + 223);

    auto g_0_xzzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_shsg + 224);

    auto g_0_yyyyy_0_xxxx_1 = pbuffer.data(idx_eri_1_shsg + 225);

    auto g_0_yyyyy_0_xxxy_1 = pbuffer.data(idx_eri_1_shsg + 226);

    auto g_0_yyyyy_0_xxxz_1 = pbuffer.data(idx_eri_1_shsg + 227);

    auto g_0_yyyyy_0_xxyy_1 = pbuffer.data(idx_eri_1_shsg + 228);

    auto g_0_yyyyy_0_xxyz_1 = pbuffer.data(idx_eri_1_shsg + 229);

    auto g_0_yyyyy_0_xxzz_1 = pbuffer.data(idx_eri_1_shsg + 230);

    auto g_0_yyyyy_0_xyyy_1 = pbuffer.data(idx_eri_1_shsg + 231);

    auto g_0_yyyyy_0_xyyz_1 = pbuffer.data(idx_eri_1_shsg + 232);

    auto g_0_yyyyy_0_xyzz_1 = pbuffer.data(idx_eri_1_shsg + 233);

    auto g_0_yyyyy_0_xzzz_1 = pbuffer.data(idx_eri_1_shsg + 234);

    auto g_0_yyyyy_0_yyyy_1 = pbuffer.data(idx_eri_1_shsg + 235);

    auto g_0_yyyyy_0_yyyz_1 = pbuffer.data(idx_eri_1_shsg + 236);

    auto g_0_yyyyy_0_yyzz_1 = pbuffer.data(idx_eri_1_shsg + 237);

    auto g_0_yyyyy_0_yzzz_1 = pbuffer.data(idx_eri_1_shsg + 238);

    auto g_0_yyyyy_0_zzzz_1 = pbuffer.data(idx_eri_1_shsg + 239);

    auto g_0_yyyyz_0_xxxy_1 = pbuffer.data(idx_eri_1_shsg + 241);

    auto g_0_yyyyz_0_xxxz_1 = pbuffer.data(idx_eri_1_shsg + 242);

    auto g_0_yyyyz_0_xxyy_1 = pbuffer.data(idx_eri_1_shsg + 243);

    auto g_0_yyyyz_0_xxyz_1 = pbuffer.data(idx_eri_1_shsg + 244);

    auto g_0_yyyyz_0_xxzz_1 = pbuffer.data(idx_eri_1_shsg + 245);

    auto g_0_yyyyz_0_xyyy_1 = pbuffer.data(idx_eri_1_shsg + 246);

    auto g_0_yyyyz_0_xyyz_1 = pbuffer.data(idx_eri_1_shsg + 247);

    auto g_0_yyyyz_0_xyzz_1 = pbuffer.data(idx_eri_1_shsg + 248);

    auto g_0_yyyyz_0_xzzz_1 = pbuffer.data(idx_eri_1_shsg + 249);

    auto g_0_yyyyz_0_yyyy_1 = pbuffer.data(idx_eri_1_shsg + 250);

    auto g_0_yyyyz_0_yyyz_1 = pbuffer.data(idx_eri_1_shsg + 251);

    auto g_0_yyyyz_0_yyzz_1 = pbuffer.data(idx_eri_1_shsg + 252);

    auto g_0_yyyyz_0_yzzz_1 = pbuffer.data(idx_eri_1_shsg + 253);

    auto g_0_yyyyz_0_zzzz_1 = pbuffer.data(idx_eri_1_shsg + 254);

    auto g_0_yyyzz_0_xxxx_1 = pbuffer.data(idx_eri_1_shsg + 255);

    auto g_0_yyyzz_0_xxxy_1 = pbuffer.data(idx_eri_1_shsg + 256);

    auto g_0_yyyzz_0_xxxz_1 = pbuffer.data(idx_eri_1_shsg + 257);

    auto g_0_yyyzz_0_xxyy_1 = pbuffer.data(idx_eri_1_shsg + 258);

    auto g_0_yyyzz_0_xxyz_1 = pbuffer.data(idx_eri_1_shsg + 259);

    auto g_0_yyyzz_0_xxzz_1 = pbuffer.data(idx_eri_1_shsg + 260);

    auto g_0_yyyzz_0_xyyy_1 = pbuffer.data(idx_eri_1_shsg + 261);

    auto g_0_yyyzz_0_xyyz_1 = pbuffer.data(idx_eri_1_shsg + 262);

    auto g_0_yyyzz_0_xyzz_1 = pbuffer.data(idx_eri_1_shsg + 263);

    auto g_0_yyyzz_0_xzzz_1 = pbuffer.data(idx_eri_1_shsg + 264);

    auto g_0_yyyzz_0_yyyy_1 = pbuffer.data(idx_eri_1_shsg + 265);

    auto g_0_yyyzz_0_yyyz_1 = pbuffer.data(idx_eri_1_shsg + 266);

    auto g_0_yyyzz_0_yyzz_1 = pbuffer.data(idx_eri_1_shsg + 267);

    auto g_0_yyyzz_0_yzzz_1 = pbuffer.data(idx_eri_1_shsg + 268);

    auto g_0_yyyzz_0_zzzz_1 = pbuffer.data(idx_eri_1_shsg + 269);

    auto g_0_yyzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_shsg + 270);

    auto g_0_yyzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_shsg + 271);

    auto g_0_yyzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_shsg + 272);

    auto g_0_yyzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_shsg + 273);

    auto g_0_yyzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_shsg + 274);

    auto g_0_yyzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_shsg + 275);

    auto g_0_yyzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_shsg + 276);

    auto g_0_yyzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_shsg + 277);

    auto g_0_yyzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_shsg + 278);

    auto g_0_yyzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_shsg + 279);

    auto g_0_yyzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_shsg + 280);

    auto g_0_yyzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_shsg + 281);

    auto g_0_yyzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_shsg + 282);

    auto g_0_yyzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_shsg + 283);

    auto g_0_yyzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_shsg + 284);

    auto g_0_yzzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_shsg + 285);

    auto g_0_yzzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_shsg + 286);

    auto g_0_yzzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_shsg + 287);

    auto g_0_yzzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_shsg + 288);

    auto g_0_yzzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_shsg + 289);

    auto g_0_yzzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_shsg + 290);

    auto g_0_yzzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_shsg + 291);

    auto g_0_yzzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_shsg + 292);

    auto g_0_yzzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_shsg + 293);

    auto g_0_yzzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_shsg + 294);

    auto g_0_yzzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_shsg + 295);

    auto g_0_yzzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_shsg + 296);

    auto g_0_yzzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_shsg + 297);

    auto g_0_yzzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_shsg + 298);

    auto g_0_yzzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_shsg + 299);

    auto g_0_zzzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_shsg + 300);

    auto g_0_zzzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_shsg + 301);

    auto g_0_zzzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_shsg + 302);

    auto g_0_zzzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_shsg + 303);

    auto g_0_zzzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_shsg + 304);

    auto g_0_zzzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_shsg + 305);

    auto g_0_zzzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_shsg + 306);

    auto g_0_zzzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_shsg + 307);

    auto g_0_zzzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_shsg + 308);

    auto g_0_zzzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_shsg + 309);

    auto g_0_zzzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_shsg + 310);

    auto g_0_zzzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_shsg + 311);

    auto g_0_zzzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_shsg + 312);

    auto g_0_zzzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_shsg + 313);

    auto g_0_zzzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_shsg + 314);

    /// Set up 0-15 components of targeted buffer : SISG

    auto g_0_xxxxxx_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg);

    auto g_0_xxxxxx_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 1);

    auto g_0_xxxxxx_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 2);

    auto g_0_xxxxxx_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 3);

    auto g_0_xxxxxx_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 4);

    auto g_0_xxxxxx_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 5);

    auto g_0_xxxxxx_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 6);

    auto g_0_xxxxxx_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 7);

    auto g_0_xxxxxx_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 8);

    auto g_0_xxxxxx_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 9);

    auto g_0_xxxxxx_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 10);

    auto g_0_xxxxxx_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 11);

    auto g_0_xxxxxx_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 12);

    auto g_0_xxxxxx_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 13);

    auto g_0_xxxxxx_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 14);

#pragma omp simd aligned(g_0_xxxx_0_xxxx_0,       \
                             g_0_xxxx_0_xxxx_1,   \
                             g_0_xxxx_0_xxxy_0,   \
                             g_0_xxxx_0_xxxy_1,   \
                             g_0_xxxx_0_xxxz_0,   \
                             g_0_xxxx_0_xxxz_1,   \
                             g_0_xxxx_0_xxyy_0,   \
                             g_0_xxxx_0_xxyy_1,   \
                             g_0_xxxx_0_xxyz_0,   \
                             g_0_xxxx_0_xxyz_1,   \
                             g_0_xxxx_0_xxzz_0,   \
                             g_0_xxxx_0_xxzz_1,   \
                             g_0_xxxx_0_xyyy_0,   \
                             g_0_xxxx_0_xyyy_1,   \
                             g_0_xxxx_0_xyyz_0,   \
                             g_0_xxxx_0_xyyz_1,   \
                             g_0_xxxx_0_xyzz_0,   \
                             g_0_xxxx_0_xyzz_1,   \
                             g_0_xxxx_0_xzzz_0,   \
                             g_0_xxxx_0_xzzz_1,   \
                             g_0_xxxx_0_yyyy_0,   \
                             g_0_xxxx_0_yyyy_1,   \
                             g_0_xxxx_0_yyyz_0,   \
                             g_0_xxxx_0_yyyz_1,   \
                             g_0_xxxx_0_yyzz_0,   \
                             g_0_xxxx_0_yyzz_1,   \
                             g_0_xxxx_0_yzzz_0,   \
                             g_0_xxxx_0_yzzz_1,   \
                             g_0_xxxx_0_zzzz_0,   \
                             g_0_xxxx_0_zzzz_1,   \
                             g_0_xxxxx_0_xxx_1,   \
                             g_0_xxxxx_0_xxxx_0,  \
                             g_0_xxxxx_0_xxxx_1,  \
                             g_0_xxxxx_0_xxxy_0,  \
                             g_0_xxxxx_0_xxxy_1,  \
                             g_0_xxxxx_0_xxxz_0,  \
                             g_0_xxxxx_0_xxxz_1,  \
                             g_0_xxxxx_0_xxy_1,   \
                             g_0_xxxxx_0_xxyy_0,  \
                             g_0_xxxxx_0_xxyy_1,  \
                             g_0_xxxxx_0_xxyz_0,  \
                             g_0_xxxxx_0_xxyz_1,  \
                             g_0_xxxxx_0_xxz_1,   \
                             g_0_xxxxx_0_xxzz_0,  \
                             g_0_xxxxx_0_xxzz_1,  \
                             g_0_xxxxx_0_xyy_1,   \
                             g_0_xxxxx_0_xyyy_0,  \
                             g_0_xxxxx_0_xyyy_1,  \
                             g_0_xxxxx_0_xyyz_0,  \
                             g_0_xxxxx_0_xyyz_1,  \
                             g_0_xxxxx_0_xyz_1,   \
                             g_0_xxxxx_0_xyzz_0,  \
                             g_0_xxxxx_0_xyzz_1,  \
                             g_0_xxxxx_0_xzz_1,   \
                             g_0_xxxxx_0_xzzz_0,  \
                             g_0_xxxxx_0_xzzz_1,  \
                             g_0_xxxxx_0_yyy_1,   \
                             g_0_xxxxx_0_yyyy_0,  \
                             g_0_xxxxx_0_yyyy_1,  \
                             g_0_xxxxx_0_yyyz_0,  \
                             g_0_xxxxx_0_yyyz_1,  \
                             g_0_xxxxx_0_yyz_1,   \
                             g_0_xxxxx_0_yyzz_0,  \
                             g_0_xxxxx_0_yyzz_1,  \
                             g_0_xxxxx_0_yzz_1,   \
                             g_0_xxxxx_0_yzzz_0,  \
                             g_0_xxxxx_0_yzzz_1,  \
                             g_0_xxxxx_0_zzz_1,   \
                             g_0_xxxxx_0_zzzz_0,  \
                             g_0_xxxxx_0_zzzz_1,  \
                             g_0_xxxxxx_0_xxxx_0, \
                             g_0_xxxxxx_0_xxxy_0, \
                             g_0_xxxxxx_0_xxxz_0, \
                             g_0_xxxxxx_0_xxyy_0, \
                             g_0_xxxxxx_0_xxyz_0, \
                             g_0_xxxxxx_0_xxzz_0, \
                             g_0_xxxxxx_0_xyyy_0, \
                             g_0_xxxxxx_0_xyyz_0, \
                             g_0_xxxxxx_0_xyzz_0, \
                             g_0_xxxxxx_0_xzzz_0, \
                             g_0_xxxxxx_0_yyyy_0, \
                             g_0_xxxxxx_0_yyyz_0, \
                             g_0_xxxxxx_0_yyzz_0, \
                             g_0_xxxxxx_0_yzzz_0, \
                             g_0_xxxxxx_0_zzzz_0, \
                             wp_x,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxx_0_xxxx_0[i] = 5.0 * g_0_xxxx_0_xxxx_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxx_1[i] * fti_ab_0 +
                                 4.0 * g_0_xxxxx_0_xxx_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxx_0[i] * pb_x + g_0_xxxxx_0_xxxx_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxy_0[i] = 5.0 * g_0_xxxx_0_xxxy_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxy_1[i] * fti_ab_0 +
                                 3.0 * g_0_xxxxx_0_xxy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxy_0[i] * pb_x + g_0_xxxxx_0_xxxy_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxz_0[i] = 5.0 * g_0_xxxx_0_xxxz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxz_1[i] * fti_ab_0 +
                                 3.0 * g_0_xxxxx_0_xxz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxz_0[i] * pb_x + g_0_xxxxx_0_xxxz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxyy_0[i] = 5.0 * g_0_xxxx_0_xxyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxyy_1[i] * fti_ab_0 +
                                 2.0 * g_0_xxxxx_0_xyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyy_0[i] * pb_x + g_0_xxxxx_0_xxyy_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxyz_0[i] = 5.0 * g_0_xxxx_0_xxyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxyz_1[i] * fti_ab_0 +
                                 2.0 * g_0_xxxxx_0_xyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyz_0[i] * pb_x + g_0_xxxxx_0_xxyz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxzz_0[i] = 5.0 * g_0_xxxx_0_xxzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxzz_1[i] * fti_ab_0 +
                                 2.0 * g_0_xxxxx_0_xzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxzz_0[i] * pb_x + g_0_xxxxx_0_xxzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xyyy_0[i] = 5.0 * g_0_xxxx_0_xyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xyyy_1[i] * fti_ab_0 + g_0_xxxxx_0_yyy_1[i] * fi_abcd_0 +
                                 g_0_xxxxx_0_xyyy_0[i] * pb_x + g_0_xxxxx_0_xyyy_1[i] * wp_x[i];

        g_0_xxxxxx_0_xyyz_0[i] = 5.0 * g_0_xxxx_0_xyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xyyz_1[i] * fti_ab_0 + g_0_xxxxx_0_yyz_1[i] * fi_abcd_0 +
                                 g_0_xxxxx_0_xyyz_0[i] * pb_x + g_0_xxxxx_0_xyyz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xyzz_0[i] = 5.0 * g_0_xxxx_0_xyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xyzz_1[i] * fti_ab_0 + g_0_xxxxx_0_yzz_1[i] * fi_abcd_0 +
                                 g_0_xxxxx_0_xyzz_0[i] * pb_x + g_0_xxxxx_0_xyzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xzzz_0[i] = 5.0 * g_0_xxxx_0_xzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xzzz_1[i] * fti_ab_0 + g_0_xxxxx_0_zzz_1[i] * fi_abcd_0 +
                                 g_0_xxxxx_0_xzzz_0[i] * pb_x + g_0_xxxxx_0_xzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_yyyy_0[i] = 5.0 * g_0_xxxx_0_yyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_yyyy_1[i] * fti_ab_0 + g_0_xxxxx_0_yyyy_0[i] * pb_x +
                                 g_0_xxxxx_0_yyyy_1[i] * wp_x[i];

        g_0_xxxxxx_0_yyyz_0[i] = 5.0 * g_0_xxxx_0_yyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_yyyz_1[i] * fti_ab_0 + g_0_xxxxx_0_yyyz_0[i] * pb_x +
                                 g_0_xxxxx_0_yyyz_1[i] * wp_x[i];

        g_0_xxxxxx_0_yyzz_0[i] = 5.0 * g_0_xxxx_0_yyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_yyzz_1[i] * fti_ab_0 + g_0_xxxxx_0_yyzz_0[i] * pb_x +
                                 g_0_xxxxx_0_yyzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_yzzz_0[i] = 5.0 * g_0_xxxx_0_yzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_yzzz_1[i] * fti_ab_0 + g_0_xxxxx_0_yzzz_0[i] * pb_x +
                                 g_0_xxxxx_0_yzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_zzzz_0[i] = 5.0 * g_0_xxxx_0_zzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_zzzz_1[i] * fti_ab_0 + g_0_xxxxx_0_zzzz_0[i] * pb_x +
                                 g_0_xxxxx_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 15-30 components of targeted buffer : SISG

    auto g_0_xxxxxy_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 15);

    auto g_0_xxxxxy_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 16);

    auto g_0_xxxxxy_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 17);

    auto g_0_xxxxxy_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 18);

    auto g_0_xxxxxy_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 19);

    auto g_0_xxxxxy_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 20);

    auto g_0_xxxxxy_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 21);

    auto g_0_xxxxxy_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 22);

    auto g_0_xxxxxy_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 23);

    auto g_0_xxxxxy_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 24);

    auto g_0_xxxxxy_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 25);

    auto g_0_xxxxxy_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 26);

    auto g_0_xxxxxy_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 27);

    auto g_0_xxxxxy_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 28);

    auto g_0_xxxxxy_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 29);

#pragma omp simd aligned(g_0_xxxxx_0_xxx_1,       \
                             g_0_xxxxx_0_xxxx_0,  \
                             g_0_xxxxx_0_xxxx_1,  \
                             g_0_xxxxx_0_xxxy_0,  \
                             g_0_xxxxx_0_xxxy_1,  \
                             g_0_xxxxx_0_xxxz_0,  \
                             g_0_xxxxx_0_xxxz_1,  \
                             g_0_xxxxx_0_xxy_1,   \
                             g_0_xxxxx_0_xxyy_0,  \
                             g_0_xxxxx_0_xxyy_1,  \
                             g_0_xxxxx_0_xxyz_0,  \
                             g_0_xxxxx_0_xxyz_1,  \
                             g_0_xxxxx_0_xxz_1,   \
                             g_0_xxxxx_0_xxzz_0,  \
                             g_0_xxxxx_0_xxzz_1,  \
                             g_0_xxxxx_0_xyy_1,   \
                             g_0_xxxxx_0_xyyy_0,  \
                             g_0_xxxxx_0_xyyy_1,  \
                             g_0_xxxxx_0_xyyz_0,  \
                             g_0_xxxxx_0_xyyz_1,  \
                             g_0_xxxxx_0_xyz_1,   \
                             g_0_xxxxx_0_xyzz_0,  \
                             g_0_xxxxx_0_xyzz_1,  \
                             g_0_xxxxx_0_xzz_1,   \
                             g_0_xxxxx_0_xzzz_0,  \
                             g_0_xxxxx_0_xzzz_1,  \
                             g_0_xxxxx_0_yyy_1,   \
                             g_0_xxxxx_0_yyyy_0,  \
                             g_0_xxxxx_0_yyyy_1,  \
                             g_0_xxxxx_0_yyyz_0,  \
                             g_0_xxxxx_0_yyyz_1,  \
                             g_0_xxxxx_0_yyz_1,   \
                             g_0_xxxxx_0_yyzz_0,  \
                             g_0_xxxxx_0_yyzz_1,  \
                             g_0_xxxxx_0_yzz_1,   \
                             g_0_xxxxx_0_yzzz_0,  \
                             g_0_xxxxx_0_yzzz_1,  \
                             g_0_xxxxx_0_zzz_1,   \
                             g_0_xxxxx_0_zzzz_0,  \
                             g_0_xxxxx_0_zzzz_1,  \
                             g_0_xxxxxy_0_xxxx_0, \
                             g_0_xxxxxy_0_xxxy_0, \
                             g_0_xxxxxy_0_xxxz_0, \
                             g_0_xxxxxy_0_xxyy_0, \
                             g_0_xxxxxy_0_xxyz_0, \
                             g_0_xxxxxy_0_xxzz_0, \
                             g_0_xxxxxy_0_xyyy_0, \
                             g_0_xxxxxy_0_xyyz_0, \
                             g_0_xxxxxy_0_xyzz_0, \
                             g_0_xxxxxy_0_xzzz_0, \
                             g_0_xxxxxy_0_yyyy_0, \
                             g_0_xxxxxy_0_yyyz_0, \
                             g_0_xxxxxy_0_yyzz_0, \
                             g_0_xxxxxy_0_yzzz_0, \
                             g_0_xxxxxy_0_zzzz_0, \
                             wp_y,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxy_0_xxxx_0[i] = g_0_xxxxx_0_xxxx_0[i] * pb_y + g_0_xxxxx_0_xxxx_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxy_0[i] = g_0_xxxxx_0_xxx_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxy_0[i] * pb_y + g_0_xxxxx_0_xxxy_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxz_0[i] = g_0_xxxxx_0_xxxz_0[i] * pb_y + g_0_xxxxx_0_xxxz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxyy_0[i] = 2.0 * g_0_xxxxx_0_xxy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyy_0[i] * pb_y + g_0_xxxxx_0_xxyy_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxyz_0[i] = g_0_xxxxx_0_xxz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyz_0[i] * pb_y + g_0_xxxxx_0_xxyz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxzz_0[i] = g_0_xxxxx_0_xxzz_0[i] * pb_y + g_0_xxxxx_0_xxzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xyyy_0[i] = 3.0 * g_0_xxxxx_0_xyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyy_0[i] * pb_y + g_0_xxxxx_0_xyyy_1[i] * wp_y[i];

        g_0_xxxxxy_0_xyyz_0[i] = 2.0 * g_0_xxxxx_0_xyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyz_0[i] * pb_y + g_0_xxxxx_0_xyyz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xyzz_0[i] = g_0_xxxxx_0_xzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyzz_0[i] * pb_y + g_0_xxxxx_0_xyzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xzzz_0[i] = g_0_xxxxx_0_xzzz_0[i] * pb_y + g_0_xxxxx_0_xzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_yyyy_0[i] = 4.0 * g_0_xxxxx_0_yyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyyy_0[i] * pb_y + g_0_xxxxx_0_yyyy_1[i] * wp_y[i];

        g_0_xxxxxy_0_yyyz_0[i] = 3.0 * g_0_xxxxx_0_yyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyyz_0[i] * pb_y + g_0_xxxxx_0_yyyz_1[i] * wp_y[i];

        g_0_xxxxxy_0_yyzz_0[i] = 2.0 * g_0_xxxxx_0_yzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyzz_0[i] * pb_y + g_0_xxxxx_0_yyzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_yzzz_0[i] = g_0_xxxxx_0_zzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yzzz_0[i] * pb_y + g_0_xxxxx_0_yzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_zzzz_0[i] = g_0_xxxxx_0_zzzz_0[i] * pb_y + g_0_xxxxx_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 30-45 components of targeted buffer : SISG

    auto g_0_xxxxxz_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 30);

    auto g_0_xxxxxz_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 31);

    auto g_0_xxxxxz_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 32);

    auto g_0_xxxxxz_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 33);

    auto g_0_xxxxxz_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 34);

    auto g_0_xxxxxz_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 35);

    auto g_0_xxxxxz_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 36);

    auto g_0_xxxxxz_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 37);

    auto g_0_xxxxxz_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 38);

    auto g_0_xxxxxz_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 39);

    auto g_0_xxxxxz_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 40);

    auto g_0_xxxxxz_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 41);

    auto g_0_xxxxxz_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 42);

    auto g_0_xxxxxz_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 43);

    auto g_0_xxxxxz_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 44);

#pragma omp simd aligned(g_0_xxxxx_0_xxx_1,       \
                             g_0_xxxxx_0_xxxx_0,  \
                             g_0_xxxxx_0_xxxx_1,  \
                             g_0_xxxxx_0_xxxy_0,  \
                             g_0_xxxxx_0_xxxy_1,  \
                             g_0_xxxxx_0_xxxz_0,  \
                             g_0_xxxxx_0_xxxz_1,  \
                             g_0_xxxxx_0_xxy_1,   \
                             g_0_xxxxx_0_xxyy_0,  \
                             g_0_xxxxx_0_xxyy_1,  \
                             g_0_xxxxx_0_xxyz_0,  \
                             g_0_xxxxx_0_xxyz_1,  \
                             g_0_xxxxx_0_xxz_1,   \
                             g_0_xxxxx_0_xxzz_0,  \
                             g_0_xxxxx_0_xxzz_1,  \
                             g_0_xxxxx_0_xyy_1,   \
                             g_0_xxxxx_0_xyyy_0,  \
                             g_0_xxxxx_0_xyyy_1,  \
                             g_0_xxxxx_0_xyyz_0,  \
                             g_0_xxxxx_0_xyyz_1,  \
                             g_0_xxxxx_0_xyz_1,   \
                             g_0_xxxxx_0_xyzz_0,  \
                             g_0_xxxxx_0_xyzz_1,  \
                             g_0_xxxxx_0_xzz_1,   \
                             g_0_xxxxx_0_xzzz_0,  \
                             g_0_xxxxx_0_xzzz_1,  \
                             g_0_xxxxx_0_yyy_1,   \
                             g_0_xxxxx_0_yyyy_0,  \
                             g_0_xxxxx_0_yyyy_1,  \
                             g_0_xxxxx_0_yyyz_0,  \
                             g_0_xxxxx_0_yyyz_1,  \
                             g_0_xxxxx_0_yyz_1,   \
                             g_0_xxxxx_0_yyzz_0,  \
                             g_0_xxxxx_0_yyzz_1,  \
                             g_0_xxxxx_0_yzz_1,   \
                             g_0_xxxxx_0_yzzz_0,  \
                             g_0_xxxxx_0_yzzz_1,  \
                             g_0_xxxxx_0_zzz_1,   \
                             g_0_xxxxx_0_zzzz_0,  \
                             g_0_xxxxx_0_zzzz_1,  \
                             g_0_xxxxxz_0_xxxx_0, \
                             g_0_xxxxxz_0_xxxy_0, \
                             g_0_xxxxxz_0_xxxz_0, \
                             g_0_xxxxxz_0_xxyy_0, \
                             g_0_xxxxxz_0_xxyz_0, \
                             g_0_xxxxxz_0_xxzz_0, \
                             g_0_xxxxxz_0_xyyy_0, \
                             g_0_xxxxxz_0_xyyz_0, \
                             g_0_xxxxxz_0_xyzz_0, \
                             g_0_xxxxxz_0_xzzz_0, \
                             g_0_xxxxxz_0_yyyy_0, \
                             g_0_xxxxxz_0_yyyz_0, \
                             g_0_xxxxxz_0_yyzz_0, \
                             g_0_xxxxxz_0_yzzz_0, \
                             g_0_xxxxxz_0_zzzz_0, \
                             wp_z,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxz_0_xxxx_0[i] = g_0_xxxxx_0_xxxx_0[i] * pb_z + g_0_xxxxx_0_xxxx_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxy_0[i] = g_0_xxxxx_0_xxxy_0[i] * pb_z + g_0_xxxxx_0_xxxy_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxz_0[i] = g_0_xxxxx_0_xxx_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxz_0[i] * pb_z + g_0_xxxxx_0_xxxz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxyy_0[i] = g_0_xxxxx_0_xxyy_0[i] * pb_z + g_0_xxxxx_0_xxyy_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxyz_0[i] = g_0_xxxxx_0_xxy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyz_0[i] * pb_z + g_0_xxxxx_0_xxyz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxzz_0[i] = 2.0 * g_0_xxxxx_0_xxz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxzz_0[i] * pb_z + g_0_xxxxx_0_xxzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xyyy_0[i] = g_0_xxxxx_0_xyyy_0[i] * pb_z + g_0_xxxxx_0_xyyy_1[i] * wp_z[i];

        g_0_xxxxxz_0_xyyz_0[i] = g_0_xxxxx_0_xyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyz_0[i] * pb_z + g_0_xxxxx_0_xyyz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xyzz_0[i] = 2.0 * g_0_xxxxx_0_xyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyzz_0[i] * pb_z + g_0_xxxxx_0_xyzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xzzz_0[i] = 3.0 * g_0_xxxxx_0_xzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xzzz_0[i] * pb_z + g_0_xxxxx_0_xzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_yyyy_0[i] = g_0_xxxxx_0_yyyy_0[i] * pb_z + g_0_xxxxx_0_yyyy_1[i] * wp_z[i];

        g_0_xxxxxz_0_yyyz_0[i] = g_0_xxxxx_0_yyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyyz_0[i] * pb_z + g_0_xxxxx_0_yyyz_1[i] * wp_z[i];

        g_0_xxxxxz_0_yyzz_0[i] = 2.0 * g_0_xxxxx_0_yyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyzz_0[i] * pb_z + g_0_xxxxx_0_yyzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_yzzz_0[i] = 3.0 * g_0_xxxxx_0_yzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yzzz_0[i] * pb_z + g_0_xxxxx_0_yzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_zzzz_0[i] = 4.0 * g_0_xxxxx_0_zzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_zzzz_0[i] * pb_z + g_0_xxxxx_0_zzzz_1[i] * wp_z[i];
    }

    /// Set up 45-60 components of targeted buffer : SISG

    auto g_0_xxxxyy_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 45);

    auto g_0_xxxxyy_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 46);

    auto g_0_xxxxyy_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 47);

    auto g_0_xxxxyy_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 48);

    auto g_0_xxxxyy_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 49);

    auto g_0_xxxxyy_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 50);

    auto g_0_xxxxyy_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 51);

    auto g_0_xxxxyy_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 52);

    auto g_0_xxxxyy_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 53);

    auto g_0_xxxxyy_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 54);

    auto g_0_xxxxyy_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 55);

    auto g_0_xxxxyy_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 56);

    auto g_0_xxxxyy_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 57);

    auto g_0_xxxxyy_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 58);

    auto g_0_xxxxyy_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 59);

#pragma omp simd aligned(g_0_xxxx_0_xxxx_0,       \
                             g_0_xxxx_0_xxxx_1,   \
                             g_0_xxxx_0_xxxz_0,   \
                             g_0_xxxx_0_xxxz_1,   \
                             g_0_xxxx_0_xxzz_0,   \
                             g_0_xxxx_0_xxzz_1,   \
                             g_0_xxxx_0_xzzz_0,   \
                             g_0_xxxx_0_xzzz_1,   \
                             g_0_xxxxy_0_xxxx_0,  \
                             g_0_xxxxy_0_xxxx_1,  \
                             g_0_xxxxy_0_xxxz_0,  \
                             g_0_xxxxy_0_xxxz_1,  \
                             g_0_xxxxy_0_xxzz_0,  \
                             g_0_xxxxy_0_xxzz_1,  \
                             g_0_xxxxy_0_xzzz_0,  \
                             g_0_xxxxy_0_xzzz_1,  \
                             g_0_xxxxyy_0_xxxx_0, \
                             g_0_xxxxyy_0_xxxy_0, \
                             g_0_xxxxyy_0_xxxz_0, \
                             g_0_xxxxyy_0_xxyy_0, \
                             g_0_xxxxyy_0_xxyz_0, \
                             g_0_xxxxyy_0_xxzz_0, \
                             g_0_xxxxyy_0_xyyy_0, \
                             g_0_xxxxyy_0_xyyz_0, \
                             g_0_xxxxyy_0_xyzz_0, \
                             g_0_xxxxyy_0_xzzz_0, \
                             g_0_xxxxyy_0_yyyy_0, \
                             g_0_xxxxyy_0_yyyz_0, \
                             g_0_xxxxyy_0_yyzz_0, \
                             g_0_xxxxyy_0_yzzz_0, \
                             g_0_xxxxyy_0_zzzz_0, \
                             g_0_xxxyy_0_xxxy_0,  \
                             g_0_xxxyy_0_xxxy_1,  \
                             g_0_xxxyy_0_xxy_1,   \
                             g_0_xxxyy_0_xxyy_0,  \
                             g_0_xxxyy_0_xxyy_1,  \
                             g_0_xxxyy_0_xxyz_0,  \
                             g_0_xxxyy_0_xxyz_1,  \
                             g_0_xxxyy_0_xyy_1,   \
                             g_0_xxxyy_0_xyyy_0,  \
                             g_0_xxxyy_0_xyyy_1,  \
                             g_0_xxxyy_0_xyyz_0,  \
                             g_0_xxxyy_0_xyyz_1,  \
                             g_0_xxxyy_0_xyz_1,   \
                             g_0_xxxyy_0_xyzz_0,  \
                             g_0_xxxyy_0_xyzz_1,  \
                             g_0_xxxyy_0_yyy_1,   \
                             g_0_xxxyy_0_yyyy_0,  \
                             g_0_xxxyy_0_yyyy_1,  \
                             g_0_xxxyy_0_yyyz_0,  \
                             g_0_xxxyy_0_yyyz_1,  \
                             g_0_xxxyy_0_yyz_1,   \
                             g_0_xxxyy_0_yyzz_0,  \
                             g_0_xxxyy_0_yyzz_1,  \
                             g_0_xxxyy_0_yzz_1,   \
                             g_0_xxxyy_0_yzzz_0,  \
                             g_0_xxxyy_0_yzzz_1,  \
                             g_0_xxxyy_0_zzzz_0,  \
                             g_0_xxxyy_0_zzzz_1,  \
                             g_0_xxyy_0_xxxy_0,   \
                             g_0_xxyy_0_xxxy_1,   \
                             g_0_xxyy_0_xxyy_0,   \
                             g_0_xxyy_0_xxyy_1,   \
                             g_0_xxyy_0_xxyz_0,   \
                             g_0_xxyy_0_xxyz_1,   \
                             g_0_xxyy_0_xyyy_0,   \
                             g_0_xxyy_0_xyyy_1,   \
                             g_0_xxyy_0_xyyz_0,   \
                             g_0_xxyy_0_xyyz_1,   \
                             g_0_xxyy_0_xyzz_0,   \
                             g_0_xxyy_0_xyzz_1,   \
                             g_0_xxyy_0_yyyy_0,   \
                             g_0_xxyy_0_yyyy_1,   \
                             g_0_xxyy_0_yyyz_0,   \
                             g_0_xxyy_0_yyyz_1,   \
                             g_0_xxyy_0_yyzz_0,   \
                             g_0_xxyy_0_yyzz_1,   \
                             g_0_xxyy_0_yzzz_0,   \
                             g_0_xxyy_0_yzzz_1,   \
                             g_0_xxyy_0_zzzz_0,   \
                             g_0_xxyy_0_zzzz_1,   \
                             wp_x,                \
                             wp_y,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxyy_0_xxxx_0[i] =
            g_0_xxxx_0_xxxx_0[i] * fi_ab_0 - g_0_xxxx_0_xxxx_1[i] * fti_ab_0 + g_0_xxxxy_0_xxxx_0[i] * pb_y + g_0_xxxxy_0_xxxx_1[i] * wp_y[i];

        g_0_xxxxyy_0_xxxy_0[i] = 3.0 * g_0_xxyy_0_xxxy_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxy_1[i] * fti_ab_0 +
                                 3.0 * g_0_xxxyy_0_xxy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxy_0[i] * pb_x + g_0_xxxyy_0_xxxy_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxxz_0[i] =
            g_0_xxxx_0_xxxz_0[i] * fi_ab_0 - g_0_xxxx_0_xxxz_1[i] * fti_ab_0 + g_0_xxxxy_0_xxxz_0[i] * pb_y + g_0_xxxxy_0_xxxz_1[i] * wp_y[i];

        g_0_xxxxyy_0_xxyy_0[i] = 3.0 * g_0_xxyy_0_xxyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxyy_1[i] * fti_ab_0 +
                                 2.0 * g_0_xxxyy_0_xyy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxyy_0[i] * pb_x + g_0_xxxyy_0_xxyy_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxyz_0[i] = 3.0 * g_0_xxyy_0_xxyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxyz_1[i] * fti_ab_0 +
                                 2.0 * g_0_xxxyy_0_xyz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxyz_0[i] * pb_x + g_0_xxxyy_0_xxyz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxzz_0[i] =
            g_0_xxxx_0_xxzz_0[i] * fi_ab_0 - g_0_xxxx_0_xxzz_1[i] * fti_ab_0 + g_0_xxxxy_0_xxzz_0[i] * pb_y + g_0_xxxxy_0_xxzz_1[i] * wp_y[i];

        g_0_xxxxyy_0_xyyy_0[i] = 3.0 * g_0_xxyy_0_xyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xyyy_1[i] * fti_ab_0 + g_0_xxxyy_0_yyy_1[i] * fi_abcd_0 +
                                 g_0_xxxyy_0_xyyy_0[i] * pb_x + g_0_xxxyy_0_xyyy_1[i] * wp_x[i];

        g_0_xxxxyy_0_xyyz_0[i] = 3.0 * g_0_xxyy_0_xyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xyyz_1[i] * fti_ab_0 + g_0_xxxyy_0_yyz_1[i] * fi_abcd_0 +
                                 g_0_xxxyy_0_xyyz_0[i] * pb_x + g_0_xxxyy_0_xyyz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xyzz_0[i] = 3.0 * g_0_xxyy_0_xyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xyzz_1[i] * fti_ab_0 + g_0_xxxyy_0_yzz_1[i] * fi_abcd_0 +
                                 g_0_xxxyy_0_xyzz_0[i] * pb_x + g_0_xxxyy_0_xyzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xzzz_0[i] =
            g_0_xxxx_0_xzzz_0[i] * fi_ab_0 - g_0_xxxx_0_xzzz_1[i] * fti_ab_0 + g_0_xxxxy_0_xzzz_0[i] * pb_y + g_0_xxxxy_0_xzzz_1[i] * wp_y[i];

        g_0_xxxxyy_0_yyyy_0[i] = 3.0 * g_0_xxyy_0_yyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_yyyy_1[i] * fti_ab_0 + g_0_xxxyy_0_yyyy_0[i] * pb_x +
                                 g_0_xxxyy_0_yyyy_1[i] * wp_x[i];

        g_0_xxxxyy_0_yyyz_0[i] = 3.0 * g_0_xxyy_0_yyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_yyyz_1[i] * fti_ab_0 + g_0_xxxyy_0_yyyz_0[i] * pb_x +
                                 g_0_xxxyy_0_yyyz_1[i] * wp_x[i];

        g_0_xxxxyy_0_yyzz_0[i] = 3.0 * g_0_xxyy_0_yyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_yyzz_1[i] * fti_ab_0 + g_0_xxxyy_0_yyzz_0[i] * pb_x +
                                 g_0_xxxyy_0_yyzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_yzzz_0[i] = 3.0 * g_0_xxyy_0_yzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_yzzz_1[i] * fti_ab_0 + g_0_xxxyy_0_yzzz_0[i] * pb_x +
                                 g_0_xxxyy_0_yzzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_zzzz_0[i] = 3.0 * g_0_xxyy_0_zzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_zzzz_1[i] * fti_ab_0 + g_0_xxxyy_0_zzzz_0[i] * pb_x +
                                 g_0_xxxyy_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 60-75 components of targeted buffer : SISG

    auto g_0_xxxxyz_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 60);

    auto g_0_xxxxyz_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 61);

    auto g_0_xxxxyz_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 62);

    auto g_0_xxxxyz_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 63);

    auto g_0_xxxxyz_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 64);

    auto g_0_xxxxyz_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 65);

    auto g_0_xxxxyz_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 66);

    auto g_0_xxxxyz_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 67);

    auto g_0_xxxxyz_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 68);

    auto g_0_xxxxyz_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 69);

    auto g_0_xxxxyz_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 70);

    auto g_0_xxxxyz_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 71);

    auto g_0_xxxxyz_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 72);

    auto g_0_xxxxyz_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 73);

    auto g_0_xxxxyz_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 74);

#pragma omp simd aligned(g_0_xxxxy_0_xxxy_0,      \
                             g_0_xxxxy_0_xxxy_1,  \
                             g_0_xxxxy_0_xxyy_0,  \
                             g_0_xxxxy_0_xxyy_1,  \
                             g_0_xxxxy_0_xyyy_0,  \
                             g_0_xxxxy_0_xyyy_1,  \
                             g_0_xxxxy_0_yyyy_0,  \
                             g_0_xxxxy_0_yyyy_1,  \
                             g_0_xxxxyz_0_xxxx_0, \
                             g_0_xxxxyz_0_xxxy_0, \
                             g_0_xxxxyz_0_xxxz_0, \
                             g_0_xxxxyz_0_xxyy_0, \
                             g_0_xxxxyz_0_xxyz_0, \
                             g_0_xxxxyz_0_xxzz_0, \
                             g_0_xxxxyz_0_xyyy_0, \
                             g_0_xxxxyz_0_xyyz_0, \
                             g_0_xxxxyz_0_xyzz_0, \
                             g_0_xxxxyz_0_xzzz_0, \
                             g_0_xxxxyz_0_yyyy_0, \
                             g_0_xxxxyz_0_yyyz_0, \
                             g_0_xxxxyz_0_yyzz_0, \
                             g_0_xxxxyz_0_yzzz_0, \
                             g_0_xxxxyz_0_zzzz_0, \
                             g_0_xxxxz_0_xxxx_0,  \
                             g_0_xxxxz_0_xxxx_1,  \
                             g_0_xxxxz_0_xxxz_0,  \
                             g_0_xxxxz_0_xxxz_1,  \
                             g_0_xxxxz_0_xxyz_0,  \
                             g_0_xxxxz_0_xxyz_1,  \
                             g_0_xxxxz_0_xxz_1,   \
                             g_0_xxxxz_0_xxzz_0,  \
                             g_0_xxxxz_0_xxzz_1,  \
                             g_0_xxxxz_0_xyyz_0,  \
                             g_0_xxxxz_0_xyyz_1,  \
                             g_0_xxxxz_0_xyz_1,   \
                             g_0_xxxxz_0_xyzz_0,  \
                             g_0_xxxxz_0_xyzz_1,  \
                             g_0_xxxxz_0_xzz_1,   \
                             g_0_xxxxz_0_xzzz_0,  \
                             g_0_xxxxz_0_xzzz_1,  \
                             g_0_xxxxz_0_yyyz_0,  \
                             g_0_xxxxz_0_yyyz_1,  \
                             g_0_xxxxz_0_yyz_1,   \
                             g_0_xxxxz_0_yyzz_0,  \
                             g_0_xxxxz_0_yyzz_1,  \
                             g_0_xxxxz_0_yzz_1,   \
                             g_0_xxxxz_0_yzzz_0,  \
                             g_0_xxxxz_0_yzzz_1,  \
                             g_0_xxxxz_0_zzz_1,   \
                             g_0_xxxxz_0_zzzz_0,  \
                             g_0_xxxxz_0_zzzz_1,  \
                             wp_y,                \
                             wp_z,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxyz_0_xxxx_0[i] = g_0_xxxxz_0_xxxx_0[i] * pb_y + g_0_xxxxz_0_xxxx_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxxy_0[i] = g_0_xxxxy_0_xxxy_0[i] * pb_z + g_0_xxxxy_0_xxxy_1[i] * wp_z[i];

        g_0_xxxxyz_0_xxxz_0[i] = g_0_xxxxz_0_xxxz_0[i] * pb_y + g_0_xxxxz_0_xxxz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxyy_0[i] = g_0_xxxxy_0_xxyy_0[i] * pb_z + g_0_xxxxy_0_xxyy_1[i] * wp_z[i];

        g_0_xxxxyz_0_xxyz_0[i] = g_0_xxxxz_0_xxz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xxyz_0[i] * pb_y + g_0_xxxxz_0_xxyz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxzz_0[i] = g_0_xxxxz_0_xxzz_0[i] * pb_y + g_0_xxxxz_0_xxzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xyyy_0[i] = g_0_xxxxy_0_xyyy_0[i] * pb_z + g_0_xxxxy_0_xyyy_1[i] * wp_z[i];

        g_0_xxxxyz_0_xyyz_0[i] = 2.0 * g_0_xxxxz_0_xyz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xyyz_0[i] * pb_y + g_0_xxxxz_0_xyyz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xyzz_0[i] = g_0_xxxxz_0_xzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xyzz_0[i] * pb_y + g_0_xxxxz_0_xyzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xzzz_0[i] = g_0_xxxxz_0_xzzz_0[i] * pb_y + g_0_xxxxz_0_xzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_yyyy_0[i] = g_0_xxxxy_0_yyyy_0[i] * pb_z + g_0_xxxxy_0_yyyy_1[i] * wp_z[i];

        g_0_xxxxyz_0_yyyz_0[i] = 3.0 * g_0_xxxxz_0_yyz_1[i] * fi_abcd_0 + g_0_xxxxz_0_yyyz_0[i] * pb_y + g_0_xxxxz_0_yyyz_1[i] * wp_y[i];

        g_0_xxxxyz_0_yyzz_0[i] = 2.0 * g_0_xxxxz_0_yzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_yyzz_0[i] * pb_y + g_0_xxxxz_0_yyzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_yzzz_0[i] = g_0_xxxxz_0_zzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_yzzz_0[i] * pb_y + g_0_xxxxz_0_yzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_zzzz_0[i] = g_0_xxxxz_0_zzzz_0[i] * pb_y + g_0_xxxxz_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 75-90 components of targeted buffer : SISG

    auto g_0_xxxxzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 75);

    auto g_0_xxxxzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 76);

    auto g_0_xxxxzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 77);

    auto g_0_xxxxzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 78);

    auto g_0_xxxxzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 79);

    auto g_0_xxxxzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 80);

    auto g_0_xxxxzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 81);

    auto g_0_xxxxzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 82);

    auto g_0_xxxxzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 83);

    auto g_0_xxxxzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 84);

    auto g_0_xxxxzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 85);

    auto g_0_xxxxzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 86);

    auto g_0_xxxxzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 87);

    auto g_0_xxxxzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 88);

    auto g_0_xxxxzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 89);

#pragma omp simd aligned(g_0_xxxx_0_xxxx_0,       \
                             g_0_xxxx_0_xxxx_1,   \
                             g_0_xxxx_0_xxxy_0,   \
                             g_0_xxxx_0_xxxy_1,   \
                             g_0_xxxx_0_xxyy_0,   \
                             g_0_xxxx_0_xxyy_1,   \
                             g_0_xxxx_0_xyyy_0,   \
                             g_0_xxxx_0_xyyy_1,   \
                             g_0_xxxxz_0_xxxx_0,  \
                             g_0_xxxxz_0_xxxx_1,  \
                             g_0_xxxxz_0_xxxy_0,  \
                             g_0_xxxxz_0_xxxy_1,  \
                             g_0_xxxxz_0_xxyy_0,  \
                             g_0_xxxxz_0_xxyy_1,  \
                             g_0_xxxxz_0_xyyy_0,  \
                             g_0_xxxxz_0_xyyy_1,  \
                             g_0_xxxxzz_0_xxxx_0, \
                             g_0_xxxxzz_0_xxxy_0, \
                             g_0_xxxxzz_0_xxxz_0, \
                             g_0_xxxxzz_0_xxyy_0, \
                             g_0_xxxxzz_0_xxyz_0, \
                             g_0_xxxxzz_0_xxzz_0, \
                             g_0_xxxxzz_0_xyyy_0, \
                             g_0_xxxxzz_0_xyyz_0, \
                             g_0_xxxxzz_0_xyzz_0, \
                             g_0_xxxxzz_0_xzzz_0, \
                             g_0_xxxxzz_0_yyyy_0, \
                             g_0_xxxxzz_0_yyyz_0, \
                             g_0_xxxxzz_0_yyzz_0, \
                             g_0_xxxxzz_0_yzzz_0, \
                             g_0_xxxxzz_0_zzzz_0, \
                             g_0_xxxzz_0_xxxz_0,  \
                             g_0_xxxzz_0_xxxz_1,  \
                             g_0_xxxzz_0_xxyz_0,  \
                             g_0_xxxzz_0_xxyz_1,  \
                             g_0_xxxzz_0_xxz_1,   \
                             g_0_xxxzz_0_xxzz_0,  \
                             g_0_xxxzz_0_xxzz_1,  \
                             g_0_xxxzz_0_xyyz_0,  \
                             g_0_xxxzz_0_xyyz_1,  \
                             g_0_xxxzz_0_xyz_1,   \
                             g_0_xxxzz_0_xyzz_0,  \
                             g_0_xxxzz_0_xyzz_1,  \
                             g_0_xxxzz_0_xzz_1,   \
                             g_0_xxxzz_0_xzzz_0,  \
                             g_0_xxxzz_0_xzzz_1,  \
                             g_0_xxxzz_0_yyyy_0,  \
                             g_0_xxxzz_0_yyyy_1,  \
                             g_0_xxxzz_0_yyyz_0,  \
                             g_0_xxxzz_0_yyyz_1,  \
                             g_0_xxxzz_0_yyz_1,   \
                             g_0_xxxzz_0_yyzz_0,  \
                             g_0_xxxzz_0_yyzz_1,  \
                             g_0_xxxzz_0_yzz_1,   \
                             g_0_xxxzz_0_yzzz_0,  \
                             g_0_xxxzz_0_yzzz_1,  \
                             g_0_xxxzz_0_zzz_1,   \
                             g_0_xxxzz_0_zzzz_0,  \
                             g_0_xxxzz_0_zzzz_1,  \
                             g_0_xxzz_0_xxxz_0,   \
                             g_0_xxzz_0_xxxz_1,   \
                             g_0_xxzz_0_xxyz_0,   \
                             g_0_xxzz_0_xxyz_1,   \
                             g_0_xxzz_0_xxzz_0,   \
                             g_0_xxzz_0_xxzz_1,   \
                             g_0_xxzz_0_xyyz_0,   \
                             g_0_xxzz_0_xyyz_1,   \
                             g_0_xxzz_0_xyzz_0,   \
                             g_0_xxzz_0_xyzz_1,   \
                             g_0_xxzz_0_xzzz_0,   \
                             g_0_xxzz_0_xzzz_1,   \
                             g_0_xxzz_0_yyyy_0,   \
                             g_0_xxzz_0_yyyy_1,   \
                             g_0_xxzz_0_yyyz_0,   \
                             g_0_xxzz_0_yyyz_1,   \
                             g_0_xxzz_0_yyzz_0,   \
                             g_0_xxzz_0_yyzz_1,   \
                             g_0_xxzz_0_yzzz_0,   \
                             g_0_xxzz_0_yzzz_1,   \
                             g_0_xxzz_0_zzzz_0,   \
                             g_0_xxzz_0_zzzz_1,   \
                             wp_x,                \
                             wp_z,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxzz_0_xxxx_0[i] =
            g_0_xxxx_0_xxxx_0[i] * fi_ab_0 - g_0_xxxx_0_xxxx_1[i] * fti_ab_0 + g_0_xxxxz_0_xxxx_0[i] * pb_z + g_0_xxxxz_0_xxxx_1[i] * wp_z[i];

        g_0_xxxxzz_0_xxxy_0[i] =
            g_0_xxxx_0_xxxy_0[i] * fi_ab_0 - g_0_xxxx_0_xxxy_1[i] * fti_ab_0 + g_0_xxxxz_0_xxxy_0[i] * pb_z + g_0_xxxxz_0_xxxy_1[i] * wp_z[i];

        g_0_xxxxzz_0_xxxz_0[i] = 3.0 * g_0_xxzz_0_xxxz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxz_1[i] * fti_ab_0 +
                                 3.0 * g_0_xxxzz_0_xxz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxz_0[i] * pb_x + g_0_xxxzz_0_xxxz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxyy_0[i] =
            g_0_xxxx_0_xxyy_0[i] * fi_ab_0 - g_0_xxxx_0_xxyy_1[i] * fti_ab_0 + g_0_xxxxz_0_xxyy_0[i] * pb_z + g_0_xxxxz_0_xxyy_1[i] * wp_z[i];

        g_0_xxxxzz_0_xxyz_0[i] = 3.0 * g_0_xxzz_0_xxyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxyz_1[i] * fti_ab_0 +
                                 2.0 * g_0_xxxzz_0_xyz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxyz_0[i] * pb_x + g_0_xxxzz_0_xxyz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxzz_0[i] = 3.0 * g_0_xxzz_0_xxzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxzz_1[i] * fti_ab_0 +
                                 2.0 * g_0_xxxzz_0_xzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxzz_0[i] * pb_x + g_0_xxxzz_0_xxzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xyyy_0[i] =
            g_0_xxxx_0_xyyy_0[i] * fi_ab_0 - g_0_xxxx_0_xyyy_1[i] * fti_ab_0 + g_0_xxxxz_0_xyyy_0[i] * pb_z + g_0_xxxxz_0_xyyy_1[i] * wp_z[i];

        g_0_xxxxzz_0_xyyz_0[i] = 3.0 * g_0_xxzz_0_xyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xyyz_1[i] * fti_ab_0 + g_0_xxxzz_0_yyz_1[i] * fi_abcd_0 +
                                 g_0_xxxzz_0_xyyz_0[i] * pb_x + g_0_xxxzz_0_xyyz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xyzz_0[i] = 3.0 * g_0_xxzz_0_xyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xyzz_1[i] * fti_ab_0 + g_0_xxxzz_0_yzz_1[i] * fi_abcd_0 +
                                 g_0_xxxzz_0_xyzz_0[i] * pb_x + g_0_xxxzz_0_xyzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xzzz_0[i] = 3.0 * g_0_xxzz_0_xzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xzzz_1[i] * fti_ab_0 + g_0_xxxzz_0_zzz_1[i] * fi_abcd_0 +
                                 g_0_xxxzz_0_xzzz_0[i] * pb_x + g_0_xxxzz_0_xzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_yyyy_0[i] = 3.0 * g_0_xxzz_0_yyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_yyyy_1[i] * fti_ab_0 + g_0_xxxzz_0_yyyy_0[i] * pb_x +
                                 g_0_xxxzz_0_yyyy_1[i] * wp_x[i];

        g_0_xxxxzz_0_yyyz_0[i] = 3.0 * g_0_xxzz_0_yyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_yyyz_1[i] * fti_ab_0 + g_0_xxxzz_0_yyyz_0[i] * pb_x +
                                 g_0_xxxzz_0_yyyz_1[i] * wp_x[i];

        g_0_xxxxzz_0_yyzz_0[i] = 3.0 * g_0_xxzz_0_yyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_yyzz_1[i] * fti_ab_0 + g_0_xxxzz_0_yyzz_0[i] * pb_x +
                                 g_0_xxxzz_0_yyzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_yzzz_0[i] = 3.0 * g_0_xxzz_0_yzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_yzzz_1[i] * fti_ab_0 + g_0_xxxzz_0_yzzz_0[i] * pb_x +
                                 g_0_xxxzz_0_yzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_zzzz_0[i] = 3.0 * g_0_xxzz_0_zzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_zzzz_1[i] * fti_ab_0 + g_0_xxxzz_0_zzzz_0[i] * pb_x +
                                 g_0_xxxzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 90-105 components of targeted buffer : SISG

    auto g_0_xxxyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 90);

    auto g_0_xxxyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 91);

    auto g_0_xxxyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 92);

    auto g_0_xxxyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 93);

    auto g_0_xxxyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 94);

    auto g_0_xxxyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 95);

    auto g_0_xxxyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 96);

    auto g_0_xxxyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 97);

    auto g_0_xxxyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 98);

    auto g_0_xxxyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 99);

    auto g_0_xxxyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 100);

    auto g_0_xxxyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 101);

    auto g_0_xxxyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 102);

    auto g_0_xxxyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 103);

    auto g_0_xxxyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 104);

#pragma omp simd aligned(g_0_xxxy_0_xxxx_0,       \
                             g_0_xxxy_0_xxxx_1,   \
                             g_0_xxxy_0_xxxz_0,   \
                             g_0_xxxy_0_xxxz_1,   \
                             g_0_xxxy_0_xxzz_0,   \
                             g_0_xxxy_0_xxzz_1,   \
                             g_0_xxxy_0_xzzz_0,   \
                             g_0_xxxy_0_xzzz_1,   \
                             g_0_xxxyy_0_xxxx_0,  \
                             g_0_xxxyy_0_xxxx_1,  \
                             g_0_xxxyy_0_xxxz_0,  \
                             g_0_xxxyy_0_xxxz_1,  \
                             g_0_xxxyy_0_xxzz_0,  \
                             g_0_xxxyy_0_xxzz_1,  \
                             g_0_xxxyy_0_xzzz_0,  \
                             g_0_xxxyy_0_xzzz_1,  \
                             g_0_xxxyyy_0_xxxx_0, \
                             g_0_xxxyyy_0_xxxy_0, \
                             g_0_xxxyyy_0_xxxz_0, \
                             g_0_xxxyyy_0_xxyy_0, \
                             g_0_xxxyyy_0_xxyz_0, \
                             g_0_xxxyyy_0_xxzz_0, \
                             g_0_xxxyyy_0_xyyy_0, \
                             g_0_xxxyyy_0_xyyz_0, \
                             g_0_xxxyyy_0_xyzz_0, \
                             g_0_xxxyyy_0_xzzz_0, \
                             g_0_xxxyyy_0_yyyy_0, \
                             g_0_xxxyyy_0_yyyz_0, \
                             g_0_xxxyyy_0_yyzz_0, \
                             g_0_xxxyyy_0_yzzz_0, \
                             g_0_xxxyyy_0_zzzz_0, \
                             g_0_xxyyy_0_xxxy_0,  \
                             g_0_xxyyy_0_xxxy_1,  \
                             g_0_xxyyy_0_xxy_1,   \
                             g_0_xxyyy_0_xxyy_0,  \
                             g_0_xxyyy_0_xxyy_1,  \
                             g_0_xxyyy_0_xxyz_0,  \
                             g_0_xxyyy_0_xxyz_1,  \
                             g_0_xxyyy_0_xyy_1,   \
                             g_0_xxyyy_0_xyyy_0,  \
                             g_0_xxyyy_0_xyyy_1,  \
                             g_0_xxyyy_0_xyyz_0,  \
                             g_0_xxyyy_0_xyyz_1,  \
                             g_0_xxyyy_0_xyz_1,   \
                             g_0_xxyyy_0_xyzz_0,  \
                             g_0_xxyyy_0_xyzz_1,  \
                             g_0_xxyyy_0_yyy_1,   \
                             g_0_xxyyy_0_yyyy_0,  \
                             g_0_xxyyy_0_yyyy_1,  \
                             g_0_xxyyy_0_yyyz_0,  \
                             g_0_xxyyy_0_yyyz_1,  \
                             g_0_xxyyy_0_yyz_1,   \
                             g_0_xxyyy_0_yyzz_0,  \
                             g_0_xxyyy_0_yyzz_1,  \
                             g_0_xxyyy_0_yzz_1,   \
                             g_0_xxyyy_0_yzzz_0,  \
                             g_0_xxyyy_0_yzzz_1,  \
                             g_0_xxyyy_0_zzzz_0,  \
                             g_0_xxyyy_0_zzzz_1,  \
                             g_0_xyyy_0_xxxy_0,   \
                             g_0_xyyy_0_xxxy_1,   \
                             g_0_xyyy_0_xxyy_0,   \
                             g_0_xyyy_0_xxyy_1,   \
                             g_0_xyyy_0_xxyz_0,   \
                             g_0_xyyy_0_xxyz_1,   \
                             g_0_xyyy_0_xyyy_0,   \
                             g_0_xyyy_0_xyyy_1,   \
                             g_0_xyyy_0_xyyz_0,   \
                             g_0_xyyy_0_xyyz_1,   \
                             g_0_xyyy_0_xyzz_0,   \
                             g_0_xyyy_0_xyzz_1,   \
                             g_0_xyyy_0_yyyy_0,   \
                             g_0_xyyy_0_yyyy_1,   \
                             g_0_xyyy_0_yyyz_0,   \
                             g_0_xyyy_0_yyyz_1,   \
                             g_0_xyyy_0_yyzz_0,   \
                             g_0_xyyy_0_yyzz_1,   \
                             g_0_xyyy_0_yzzz_0,   \
                             g_0_xyyy_0_yzzz_1,   \
                             g_0_xyyy_0_zzzz_0,   \
                             g_0_xyyy_0_zzzz_1,   \
                             wp_x,                \
                             wp_y,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxyyy_0_xxxx_0[i] = 2.0 * g_0_xxxy_0_xxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxxy_0_xxxx_1[i] * fti_ab_0 + g_0_xxxyy_0_xxxx_0[i] * pb_y +
                                 g_0_xxxyy_0_xxxx_1[i] * wp_y[i];

        g_0_xxxyyy_0_xxxy_0[i] = 2.0 * g_0_xyyy_0_xxxy_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxxy_1[i] * fti_ab_0 +
                                 3.0 * g_0_xxyyy_0_xxy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxy_0[i] * pb_x + g_0_xxyyy_0_xxxy_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxxz_0[i] = 2.0 * g_0_xxxy_0_xxxz_0[i] * fi_ab_0 - 2.0 * g_0_xxxy_0_xxxz_1[i] * fti_ab_0 + g_0_xxxyy_0_xxxz_0[i] * pb_y +
                                 g_0_xxxyy_0_xxxz_1[i] * wp_y[i];

        g_0_xxxyyy_0_xxyy_0[i] = 2.0 * g_0_xyyy_0_xxyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxyy_1[i] * fti_ab_0 +
                                 2.0 * g_0_xxyyy_0_xyy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxyy_0[i] * pb_x + g_0_xxyyy_0_xxyy_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxyz_0[i] = 2.0 * g_0_xyyy_0_xxyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxyz_1[i] * fti_ab_0 +
                                 2.0 * g_0_xxyyy_0_xyz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxyz_0[i] * pb_x + g_0_xxyyy_0_xxyz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxzz_0[i] = 2.0 * g_0_xxxy_0_xxzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxy_0_xxzz_1[i] * fti_ab_0 + g_0_xxxyy_0_xxzz_0[i] * pb_y +
                                 g_0_xxxyy_0_xxzz_1[i] * wp_y[i];

        g_0_xxxyyy_0_xyyy_0[i] = 2.0 * g_0_xyyy_0_xyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xyyy_1[i] * fti_ab_0 + g_0_xxyyy_0_yyy_1[i] * fi_abcd_0 +
                                 g_0_xxyyy_0_xyyy_0[i] * pb_x + g_0_xxyyy_0_xyyy_1[i] * wp_x[i];

        g_0_xxxyyy_0_xyyz_0[i] = 2.0 * g_0_xyyy_0_xyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xyyz_1[i] * fti_ab_0 + g_0_xxyyy_0_yyz_1[i] * fi_abcd_0 +
                                 g_0_xxyyy_0_xyyz_0[i] * pb_x + g_0_xxyyy_0_xyyz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xyzz_0[i] = 2.0 * g_0_xyyy_0_xyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xyzz_1[i] * fti_ab_0 + g_0_xxyyy_0_yzz_1[i] * fi_abcd_0 +
                                 g_0_xxyyy_0_xyzz_0[i] * pb_x + g_0_xxyyy_0_xyzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xzzz_0[i] = 2.0 * g_0_xxxy_0_xzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxy_0_xzzz_1[i] * fti_ab_0 + g_0_xxxyy_0_xzzz_0[i] * pb_y +
                                 g_0_xxxyy_0_xzzz_1[i] * wp_y[i];

        g_0_xxxyyy_0_yyyy_0[i] = 2.0 * g_0_xyyy_0_yyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_yyyy_1[i] * fti_ab_0 + g_0_xxyyy_0_yyyy_0[i] * pb_x +
                                 g_0_xxyyy_0_yyyy_1[i] * wp_x[i];

        g_0_xxxyyy_0_yyyz_0[i] = 2.0 * g_0_xyyy_0_yyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_yyyz_1[i] * fti_ab_0 + g_0_xxyyy_0_yyyz_0[i] * pb_x +
                                 g_0_xxyyy_0_yyyz_1[i] * wp_x[i];

        g_0_xxxyyy_0_yyzz_0[i] = 2.0 * g_0_xyyy_0_yyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_yyzz_1[i] * fti_ab_0 + g_0_xxyyy_0_yyzz_0[i] * pb_x +
                                 g_0_xxyyy_0_yyzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_yzzz_0[i] = 2.0 * g_0_xyyy_0_yzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_yzzz_1[i] * fti_ab_0 + g_0_xxyyy_0_yzzz_0[i] * pb_x +
                                 g_0_xxyyy_0_yzzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_zzzz_0[i] = 2.0 * g_0_xyyy_0_zzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_zzzz_1[i] * fti_ab_0 + g_0_xxyyy_0_zzzz_0[i] * pb_x +
                                 g_0_xxyyy_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 105-120 components of targeted buffer : SISG

    auto g_0_xxxyyz_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 105);

    auto g_0_xxxyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 106);

    auto g_0_xxxyyz_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 107);

    auto g_0_xxxyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 108);

    auto g_0_xxxyyz_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 109);

    auto g_0_xxxyyz_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 110);

    auto g_0_xxxyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 111);

    auto g_0_xxxyyz_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 112);

    auto g_0_xxxyyz_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 113);

    auto g_0_xxxyyz_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 114);

    auto g_0_xxxyyz_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 115);

    auto g_0_xxxyyz_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 116);

    auto g_0_xxxyyz_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 117);

    auto g_0_xxxyyz_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 118);

    auto g_0_xxxyyz_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 119);

#pragma omp simd aligned(g_0_xxxyy_0_xxx_1,       \
                             g_0_xxxyy_0_xxxx_0,  \
                             g_0_xxxyy_0_xxxx_1,  \
                             g_0_xxxyy_0_xxxy_0,  \
                             g_0_xxxyy_0_xxxy_1,  \
                             g_0_xxxyy_0_xxxz_0,  \
                             g_0_xxxyy_0_xxxz_1,  \
                             g_0_xxxyy_0_xxy_1,   \
                             g_0_xxxyy_0_xxyy_0,  \
                             g_0_xxxyy_0_xxyy_1,  \
                             g_0_xxxyy_0_xxyz_0,  \
                             g_0_xxxyy_0_xxyz_1,  \
                             g_0_xxxyy_0_xxz_1,   \
                             g_0_xxxyy_0_xxzz_0,  \
                             g_0_xxxyy_0_xxzz_1,  \
                             g_0_xxxyy_0_xyy_1,   \
                             g_0_xxxyy_0_xyyy_0,  \
                             g_0_xxxyy_0_xyyy_1,  \
                             g_0_xxxyy_0_xyyz_0,  \
                             g_0_xxxyy_0_xyyz_1,  \
                             g_0_xxxyy_0_xyz_1,   \
                             g_0_xxxyy_0_xyzz_0,  \
                             g_0_xxxyy_0_xyzz_1,  \
                             g_0_xxxyy_0_xzz_1,   \
                             g_0_xxxyy_0_xzzz_0,  \
                             g_0_xxxyy_0_xzzz_1,  \
                             g_0_xxxyy_0_yyy_1,   \
                             g_0_xxxyy_0_yyyy_0,  \
                             g_0_xxxyy_0_yyyy_1,  \
                             g_0_xxxyy_0_yyyz_0,  \
                             g_0_xxxyy_0_yyyz_1,  \
                             g_0_xxxyy_0_yyz_1,   \
                             g_0_xxxyy_0_yyzz_0,  \
                             g_0_xxxyy_0_yyzz_1,  \
                             g_0_xxxyy_0_yzz_1,   \
                             g_0_xxxyy_0_yzzz_0,  \
                             g_0_xxxyy_0_yzzz_1,  \
                             g_0_xxxyy_0_zzz_1,   \
                             g_0_xxxyy_0_zzzz_0,  \
                             g_0_xxxyy_0_zzzz_1,  \
                             g_0_xxxyyz_0_xxxx_0, \
                             g_0_xxxyyz_0_xxxy_0, \
                             g_0_xxxyyz_0_xxxz_0, \
                             g_0_xxxyyz_0_xxyy_0, \
                             g_0_xxxyyz_0_xxyz_0, \
                             g_0_xxxyyz_0_xxzz_0, \
                             g_0_xxxyyz_0_xyyy_0, \
                             g_0_xxxyyz_0_xyyz_0, \
                             g_0_xxxyyz_0_xyzz_0, \
                             g_0_xxxyyz_0_xzzz_0, \
                             g_0_xxxyyz_0_yyyy_0, \
                             g_0_xxxyyz_0_yyyz_0, \
                             g_0_xxxyyz_0_yyzz_0, \
                             g_0_xxxyyz_0_yzzz_0, \
                             g_0_xxxyyz_0_zzzz_0, \
                             wp_z,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyyz_0_xxxx_0[i] = g_0_xxxyy_0_xxxx_0[i] * pb_z + g_0_xxxyy_0_xxxx_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxy_0[i] = g_0_xxxyy_0_xxxy_0[i] * pb_z + g_0_xxxyy_0_xxxy_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxz_0[i] = g_0_xxxyy_0_xxx_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxz_0[i] * pb_z + g_0_xxxyy_0_xxxz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxyy_0[i] = g_0_xxxyy_0_xxyy_0[i] * pb_z + g_0_xxxyy_0_xxyy_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxyz_0[i] = g_0_xxxyy_0_xxy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxyz_0[i] * pb_z + g_0_xxxyy_0_xxyz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxzz_0[i] = 2.0 * g_0_xxxyy_0_xxz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxzz_0[i] * pb_z + g_0_xxxyy_0_xxzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xyyy_0[i] = g_0_xxxyy_0_xyyy_0[i] * pb_z + g_0_xxxyy_0_xyyy_1[i] * wp_z[i];

        g_0_xxxyyz_0_xyyz_0[i] = g_0_xxxyy_0_xyy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyyz_0[i] * pb_z + g_0_xxxyy_0_xyyz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xyzz_0[i] = 2.0 * g_0_xxxyy_0_xyz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyzz_0[i] * pb_z + g_0_xxxyy_0_xyzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xzzz_0[i] = 3.0 * g_0_xxxyy_0_xzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xzzz_0[i] * pb_z + g_0_xxxyy_0_xzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_yyyy_0[i] = g_0_xxxyy_0_yyyy_0[i] * pb_z + g_0_xxxyy_0_yyyy_1[i] * wp_z[i];

        g_0_xxxyyz_0_yyyz_0[i] = g_0_xxxyy_0_yyy_1[i] * fi_abcd_0 + g_0_xxxyy_0_yyyz_0[i] * pb_z + g_0_xxxyy_0_yyyz_1[i] * wp_z[i];

        g_0_xxxyyz_0_yyzz_0[i] = 2.0 * g_0_xxxyy_0_yyz_1[i] * fi_abcd_0 + g_0_xxxyy_0_yyzz_0[i] * pb_z + g_0_xxxyy_0_yyzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_yzzz_0[i] = 3.0 * g_0_xxxyy_0_yzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_yzzz_0[i] * pb_z + g_0_xxxyy_0_yzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_zzzz_0[i] = 4.0 * g_0_xxxyy_0_zzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_zzzz_0[i] * pb_z + g_0_xxxyy_0_zzzz_1[i] * wp_z[i];
    }

    /// Set up 120-135 components of targeted buffer : SISG

    auto g_0_xxxyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 120);

    auto g_0_xxxyzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 121);

    auto g_0_xxxyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 122);

    auto g_0_xxxyzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 123);

    auto g_0_xxxyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 124);

    auto g_0_xxxyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 125);

    auto g_0_xxxyzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 126);

    auto g_0_xxxyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 127);

    auto g_0_xxxyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 128);

    auto g_0_xxxyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 129);

    auto g_0_xxxyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 130);

    auto g_0_xxxyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 131);

    auto g_0_xxxyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 132);

    auto g_0_xxxyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 133);

    auto g_0_xxxyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 134);

#pragma omp simd aligned(g_0_xxxyzz_0_xxxx_0,     \
                             g_0_xxxyzz_0_xxxy_0, \
                             g_0_xxxyzz_0_xxxz_0, \
                             g_0_xxxyzz_0_xxyy_0, \
                             g_0_xxxyzz_0_xxyz_0, \
                             g_0_xxxyzz_0_xxzz_0, \
                             g_0_xxxyzz_0_xyyy_0, \
                             g_0_xxxyzz_0_xyyz_0, \
                             g_0_xxxyzz_0_xyzz_0, \
                             g_0_xxxyzz_0_xzzz_0, \
                             g_0_xxxyzz_0_yyyy_0, \
                             g_0_xxxyzz_0_yyyz_0, \
                             g_0_xxxyzz_0_yyzz_0, \
                             g_0_xxxyzz_0_yzzz_0, \
                             g_0_xxxyzz_0_zzzz_0, \
                             g_0_xxxzz_0_xxx_1,   \
                             g_0_xxxzz_0_xxxx_0,  \
                             g_0_xxxzz_0_xxxx_1,  \
                             g_0_xxxzz_0_xxxy_0,  \
                             g_0_xxxzz_0_xxxy_1,  \
                             g_0_xxxzz_0_xxxz_0,  \
                             g_0_xxxzz_0_xxxz_1,  \
                             g_0_xxxzz_0_xxy_1,   \
                             g_0_xxxzz_0_xxyy_0,  \
                             g_0_xxxzz_0_xxyy_1,  \
                             g_0_xxxzz_0_xxyz_0,  \
                             g_0_xxxzz_0_xxyz_1,  \
                             g_0_xxxzz_0_xxz_1,   \
                             g_0_xxxzz_0_xxzz_0,  \
                             g_0_xxxzz_0_xxzz_1,  \
                             g_0_xxxzz_0_xyy_1,   \
                             g_0_xxxzz_0_xyyy_0,  \
                             g_0_xxxzz_0_xyyy_1,  \
                             g_0_xxxzz_0_xyyz_0,  \
                             g_0_xxxzz_0_xyyz_1,  \
                             g_0_xxxzz_0_xyz_1,   \
                             g_0_xxxzz_0_xyzz_0,  \
                             g_0_xxxzz_0_xyzz_1,  \
                             g_0_xxxzz_0_xzz_1,   \
                             g_0_xxxzz_0_xzzz_0,  \
                             g_0_xxxzz_0_xzzz_1,  \
                             g_0_xxxzz_0_yyy_1,   \
                             g_0_xxxzz_0_yyyy_0,  \
                             g_0_xxxzz_0_yyyy_1,  \
                             g_0_xxxzz_0_yyyz_0,  \
                             g_0_xxxzz_0_yyyz_1,  \
                             g_0_xxxzz_0_yyz_1,   \
                             g_0_xxxzz_0_yyzz_0,  \
                             g_0_xxxzz_0_yyzz_1,  \
                             g_0_xxxzz_0_yzz_1,   \
                             g_0_xxxzz_0_yzzz_0,  \
                             g_0_xxxzz_0_yzzz_1,  \
                             g_0_xxxzz_0_zzz_1,   \
                             g_0_xxxzz_0_zzzz_0,  \
                             g_0_xxxzz_0_zzzz_1,  \
                             wp_y,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyzz_0_xxxx_0[i] = g_0_xxxzz_0_xxxx_0[i] * pb_y + g_0_xxxzz_0_xxxx_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxy_0[i] = g_0_xxxzz_0_xxx_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxy_0[i] * pb_y + g_0_xxxzz_0_xxxy_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxz_0[i] = g_0_xxxzz_0_xxxz_0[i] * pb_y + g_0_xxxzz_0_xxxz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxyy_0[i] = 2.0 * g_0_xxxzz_0_xxy_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxyy_0[i] * pb_y + g_0_xxxzz_0_xxyy_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxyz_0[i] = g_0_xxxzz_0_xxz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxyz_0[i] * pb_y + g_0_xxxzz_0_xxyz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxzz_0[i] = g_0_xxxzz_0_xxzz_0[i] * pb_y + g_0_xxxzz_0_xxzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xyyy_0[i] = 3.0 * g_0_xxxzz_0_xyy_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyyy_0[i] * pb_y + g_0_xxxzz_0_xyyy_1[i] * wp_y[i];

        g_0_xxxyzz_0_xyyz_0[i] = 2.0 * g_0_xxxzz_0_xyz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyyz_0[i] * pb_y + g_0_xxxzz_0_xyyz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xyzz_0[i] = g_0_xxxzz_0_xzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyzz_0[i] * pb_y + g_0_xxxzz_0_xyzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xzzz_0[i] = g_0_xxxzz_0_xzzz_0[i] * pb_y + g_0_xxxzz_0_xzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_yyyy_0[i] = 4.0 * g_0_xxxzz_0_yyy_1[i] * fi_abcd_0 + g_0_xxxzz_0_yyyy_0[i] * pb_y + g_0_xxxzz_0_yyyy_1[i] * wp_y[i];

        g_0_xxxyzz_0_yyyz_0[i] = 3.0 * g_0_xxxzz_0_yyz_1[i] * fi_abcd_0 + g_0_xxxzz_0_yyyz_0[i] * pb_y + g_0_xxxzz_0_yyyz_1[i] * wp_y[i];

        g_0_xxxyzz_0_yyzz_0[i] = 2.0 * g_0_xxxzz_0_yzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_yyzz_0[i] * pb_y + g_0_xxxzz_0_yyzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_yzzz_0[i] = g_0_xxxzz_0_zzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_yzzz_0[i] * pb_y + g_0_xxxzz_0_yzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_zzzz_0[i] = g_0_xxxzz_0_zzzz_0[i] * pb_y + g_0_xxxzz_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 135-150 components of targeted buffer : SISG

    auto g_0_xxxzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 135);

    auto g_0_xxxzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 136);

    auto g_0_xxxzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 137);

    auto g_0_xxxzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 138);

    auto g_0_xxxzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 139);

    auto g_0_xxxzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 140);

    auto g_0_xxxzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 141);

    auto g_0_xxxzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 142);

    auto g_0_xxxzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 143);

    auto g_0_xxxzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 144);

    auto g_0_xxxzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 145);

    auto g_0_xxxzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 146);

    auto g_0_xxxzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 147);

    auto g_0_xxxzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 148);

    auto g_0_xxxzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 149);

#pragma omp simd aligned(g_0_xxxz_0_xxxx_0,       \
                             g_0_xxxz_0_xxxx_1,   \
                             g_0_xxxz_0_xxxy_0,   \
                             g_0_xxxz_0_xxxy_1,   \
                             g_0_xxxz_0_xxyy_0,   \
                             g_0_xxxz_0_xxyy_1,   \
                             g_0_xxxz_0_xyyy_0,   \
                             g_0_xxxz_0_xyyy_1,   \
                             g_0_xxxzz_0_xxxx_0,  \
                             g_0_xxxzz_0_xxxx_1,  \
                             g_0_xxxzz_0_xxxy_0,  \
                             g_0_xxxzz_0_xxxy_1,  \
                             g_0_xxxzz_0_xxyy_0,  \
                             g_0_xxxzz_0_xxyy_1,  \
                             g_0_xxxzz_0_xyyy_0,  \
                             g_0_xxxzz_0_xyyy_1,  \
                             g_0_xxxzzz_0_xxxx_0, \
                             g_0_xxxzzz_0_xxxy_0, \
                             g_0_xxxzzz_0_xxxz_0, \
                             g_0_xxxzzz_0_xxyy_0, \
                             g_0_xxxzzz_0_xxyz_0, \
                             g_0_xxxzzz_0_xxzz_0, \
                             g_0_xxxzzz_0_xyyy_0, \
                             g_0_xxxzzz_0_xyyz_0, \
                             g_0_xxxzzz_0_xyzz_0, \
                             g_0_xxxzzz_0_xzzz_0, \
                             g_0_xxxzzz_0_yyyy_0, \
                             g_0_xxxzzz_0_yyyz_0, \
                             g_0_xxxzzz_0_yyzz_0, \
                             g_0_xxxzzz_0_yzzz_0, \
                             g_0_xxxzzz_0_zzzz_0, \
                             g_0_xxzzz_0_xxxz_0,  \
                             g_0_xxzzz_0_xxxz_1,  \
                             g_0_xxzzz_0_xxyz_0,  \
                             g_0_xxzzz_0_xxyz_1,  \
                             g_0_xxzzz_0_xxz_1,   \
                             g_0_xxzzz_0_xxzz_0,  \
                             g_0_xxzzz_0_xxzz_1,  \
                             g_0_xxzzz_0_xyyz_0,  \
                             g_0_xxzzz_0_xyyz_1,  \
                             g_0_xxzzz_0_xyz_1,   \
                             g_0_xxzzz_0_xyzz_0,  \
                             g_0_xxzzz_0_xyzz_1,  \
                             g_0_xxzzz_0_xzz_1,   \
                             g_0_xxzzz_0_xzzz_0,  \
                             g_0_xxzzz_0_xzzz_1,  \
                             g_0_xxzzz_0_yyyy_0,  \
                             g_0_xxzzz_0_yyyy_1,  \
                             g_0_xxzzz_0_yyyz_0,  \
                             g_0_xxzzz_0_yyyz_1,  \
                             g_0_xxzzz_0_yyz_1,   \
                             g_0_xxzzz_0_yyzz_0,  \
                             g_0_xxzzz_0_yyzz_1,  \
                             g_0_xxzzz_0_yzz_1,   \
                             g_0_xxzzz_0_yzzz_0,  \
                             g_0_xxzzz_0_yzzz_1,  \
                             g_0_xxzzz_0_zzz_1,   \
                             g_0_xxzzz_0_zzzz_0,  \
                             g_0_xxzzz_0_zzzz_1,  \
                             g_0_xzzz_0_xxxz_0,   \
                             g_0_xzzz_0_xxxz_1,   \
                             g_0_xzzz_0_xxyz_0,   \
                             g_0_xzzz_0_xxyz_1,   \
                             g_0_xzzz_0_xxzz_0,   \
                             g_0_xzzz_0_xxzz_1,   \
                             g_0_xzzz_0_xyyz_0,   \
                             g_0_xzzz_0_xyyz_1,   \
                             g_0_xzzz_0_xyzz_0,   \
                             g_0_xzzz_0_xyzz_1,   \
                             g_0_xzzz_0_xzzz_0,   \
                             g_0_xzzz_0_xzzz_1,   \
                             g_0_xzzz_0_yyyy_0,   \
                             g_0_xzzz_0_yyyy_1,   \
                             g_0_xzzz_0_yyyz_0,   \
                             g_0_xzzz_0_yyyz_1,   \
                             g_0_xzzz_0_yyzz_0,   \
                             g_0_xzzz_0_yyzz_1,   \
                             g_0_xzzz_0_yzzz_0,   \
                             g_0_xzzz_0_yzzz_1,   \
                             g_0_xzzz_0_zzzz_0,   \
                             g_0_xzzz_0_zzzz_1,   \
                             wp_x,                \
                             wp_z,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxzzz_0_xxxx_0[i] = 2.0 * g_0_xxxz_0_xxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxxz_0_xxxx_1[i] * fti_ab_0 + g_0_xxxzz_0_xxxx_0[i] * pb_z +
                                 g_0_xxxzz_0_xxxx_1[i] * wp_z[i];

        g_0_xxxzzz_0_xxxy_0[i] = 2.0 * g_0_xxxz_0_xxxy_0[i] * fi_ab_0 - 2.0 * g_0_xxxz_0_xxxy_1[i] * fti_ab_0 + g_0_xxxzz_0_xxxy_0[i] * pb_z +
                                 g_0_xxxzz_0_xxxy_1[i] * wp_z[i];

        g_0_xxxzzz_0_xxxz_0[i] = 2.0 * g_0_xzzz_0_xxxz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxxz_1[i] * fti_ab_0 +
                                 3.0 * g_0_xxzzz_0_xxz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxz_0[i] * pb_x + g_0_xxzzz_0_xxxz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxyy_0[i] = 2.0 * g_0_xxxz_0_xxyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxz_0_xxyy_1[i] * fti_ab_0 + g_0_xxxzz_0_xxyy_0[i] * pb_z +
                                 g_0_xxxzz_0_xxyy_1[i] * wp_z[i];

        g_0_xxxzzz_0_xxyz_0[i] = 2.0 * g_0_xzzz_0_xxyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxyz_1[i] * fti_ab_0 +
                                 2.0 * g_0_xxzzz_0_xyz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxyz_0[i] * pb_x + g_0_xxzzz_0_xxyz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxzz_0[i] = 2.0 * g_0_xzzz_0_xxzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxzz_1[i] * fti_ab_0 +
                                 2.0 * g_0_xxzzz_0_xzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxzz_0[i] * pb_x + g_0_xxzzz_0_xxzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xyyy_0[i] = 2.0 * g_0_xxxz_0_xyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxz_0_xyyy_1[i] * fti_ab_0 + g_0_xxxzz_0_xyyy_0[i] * pb_z +
                                 g_0_xxxzz_0_xyyy_1[i] * wp_z[i];

        g_0_xxxzzz_0_xyyz_0[i] = 2.0 * g_0_xzzz_0_xyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xyyz_1[i] * fti_ab_0 + g_0_xxzzz_0_yyz_1[i] * fi_abcd_0 +
                                 g_0_xxzzz_0_xyyz_0[i] * pb_x + g_0_xxzzz_0_xyyz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xyzz_0[i] = 2.0 * g_0_xzzz_0_xyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xyzz_1[i] * fti_ab_0 + g_0_xxzzz_0_yzz_1[i] * fi_abcd_0 +
                                 g_0_xxzzz_0_xyzz_0[i] * pb_x + g_0_xxzzz_0_xyzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xzzz_0[i] = 2.0 * g_0_xzzz_0_xzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xzzz_1[i] * fti_ab_0 + g_0_xxzzz_0_zzz_1[i] * fi_abcd_0 +
                                 g_0_xxzzz_0_xzzz_0[i] * pb_x + g_0_xxzzz_0_xzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_yyyy_0[i] = 2.0 * g_0_xzzz_0_yyyy_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_yyyy_1[i] * fti_ab_0 + g_0_xxzzz_0_yyyy_0[i] * pb_x +
                                 g_0_xxzzz_0_yyyy_1[i] * wp_x[i];

        g_0_xxxzzz_0_yyyz_0[i] = 2.0 * g_0_xzzz_0_yyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_yyyz_1[i] * fti_ab_0 + g_0_xxzzz_0_yyyz_0[i] * pb_x +
                                 g_0_xxzzz_0_yyyz_1[i] * wp_x[i];

        g_0_xxxzzz_0_yyzz_0[i] = 2.0 * g_0_xzzz_0_yyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_yyzz_1[i] * fti_ab_0 + g_0_xxzzz_0_yyzz_0[i] * pb_x +
                                 g_0_xxzzz_0_yyzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_yzzz_0[i] = 2.0 * g_0_xzzz_0_yzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_yzzz_1[i] * fti_ab_0 + g_0_xxzzz_0_yzzz_0[i] * pb_x +
                                 g_0_xxzzz_0_yzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_zzzz_0[i] = 2.0 * g_0_xzzz_0_zzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_zzzz_1[i] * fti_ab_0 + g_0_xxzzz_0_zzzz_0[i] * pb_x +
                                 g_0_xxzzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 150-165 components of targeted buffer : SISG

    auto g_0_xxyyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 150);

    auto g_0_xxyyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 151);

    auto g_0_xxyyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 152);

    auto g_0_xxyyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 153);

    auto g_0_xxyyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 154);

    auto g_0_xxyyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 155);

    auto g_0_xxyyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 156);

    auto g_0_xxyyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 157);

    auto g_0_xxyyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 158);

    auto g_0_xxyyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 159);

    auto g_0_xxyyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 160);

    auto g_0_xxyyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 161);

    auto g_0_xxyyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 162);

    auto g_0_xxyyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 163);

    auto g_0_xxyyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 164);

#pragma omp simd aligned(g_0_xxyy_0_xxxx_0,       \
                             g_0_xxyy_0_xxxx_1,   \
                             g_0_xxyy_0_xxxz_0,   \
                             g_0_xxyy_0_xxxz_1,   \
                             g_0_xxyy_0_xxzz_0,   \
                             g_0_xxyy_0_xxzz_1,   \
                             g_0_xxyy_0_xzzz_0,   \
                             g_0_xxyy_0_xzzz_1,   \
                             g_0_xxyyy_0_xxxx_0,  \
                             g_0_xxyyy_0_xxxx_1,  \
                             g_0_xxyyy_0_xxxz_0,  \
                             g_0_xxyyy_0_xxxz_1,  \
                             g_0_xxyyy_0_xxzz_0,  \
                             g_0_xxyyy_0_xxzz_1,  \
                             g_0_xxyyy_0_xzzz_0,  \
                             g_0_xxyyy_0_xzzz_1,  \
                             g_0_xxyyyy_0_xxxx_0, \
                             g_0_xxyyyy_0_xxxy_0, \
                             g_0_xxyyyy_0_xxxz_0, \
                             g_0_xxyyyy_0_xxyy_0, \
                             g_0_xxyyyy_0_xxyz_0, \
                             g_0_xxyyyy_0_xxzz_0, \
                             g_0_xxyyyy_0_xyyy_0, \
                             g_0_xxyyyy_0_xyyz_0, \
                             g_0_xxyyyy_0_xyzz_0, \
                             g_0_xxyyyy_0_xzzz_0, \
                             g_0_xxyyyy_0_yyyy_0, \
                             g_0_xxyyyy_0_yyyz_0, \
                             g_0_xxyyyy_0_yyzz_0, \
                             g_0_xxyyyy_0_yzzz_0, \
                             g_0_xxyyyy_0_zzzz_0, \
                             g_0_xyyyy_0_xxxy_0,  \
                             g_0_xyyyy_0_xxxy_1,  \
                             g_0_xyyyy_0_xxy_1,   \
                             g_0_xyyyy_0_xxyy_0,  \
                             g_0_xyyyy_0_xxyy_1,  \
                             g_0_xyyyy_0_xxyz_0,  \
                             g_0_xyyyy_0_xxyz_1,  \
                             g_0_xyyyy_0_xyy_1,   \
                             g_0_xyyyy_0_xyyy_0,  \
                             g_0_xyyyy_0_xyyy_1,  \
                             g_0_xyyyy_0_xyyz_0,  \
                             g_0_xyyyy_0_xyyz_1,  \
                             g_0_xyyyy_0_xyz_1,   \
                             g_0_xyyyy_0_xyzz_0,  \
                             g_0_xyyyy_0_xyzz_1,  \
                             g_0_xyyyy_0_yyy_1,   \
                             g_0_xyyyy_0_yyyy_0,  \
                             g_0_xyyyy_0_yyyy_1,  \
                             g_0_xyyyy_0_yyyz_0,  \
                             g_0_xyyyy_0_yyyz_1,  \
                             g_0_xyyyy_0_yyz_1,   \
                             g_0_xyyyy_0_yyzz_0,  \
                             g_0_xyyyy_0_yyzz_1,  \
                             g_0_xyyyy_0_yzz_1,   \
                             g_0_xyyyy_0_yzzz_0,  \
                             g_0_xyyyy_0_yzzz_1,  \
                             g_0_xyyyy_0_zzzz_0,  \
                             g_0_xyyyy_0_zzzz_1,  \
                             g_0_yyyy_0_xxxy_0,   \
                             g_0_yyyy_0_xxxy_1,   \
                             g_0_yyyy_0_xxyy_0,   \
                             g_0_yyyy_0_xxyy_1,   \
                             g_0_yyyy_0_xxyz_0,   \
                             g_0_yyyy_0_xxyz_1,   \
                             g_0_yyyy_0_xyyy_0,   \
                             g_0_yyyy_0_xyyy_1,   \
                             g_0_yyyy_0_xyyz_0,   \
                             g_0_yyyy_0_xyyz_1,   \
                             g_0_yyyy_0_xyzz_0,   \
                             g_0_yyyy_0_xyzz_1,   \
                             g_0_yyyy_0_yyyy_0,   \
                             g_0_yyyy_0_yyyy_1,   \
                             g_0_yyyy_0_yyyz_0,   \
                             g_0_yyyy_0_yyyz_1,   \
                             g_0_yyyy_0_yyzz_0,   \
                             g_0_yyyy_0_yyzz_1,   \
                             g_0_yyyy_0_yzzz_0,   \
                             g_0_yyyy_0_yzzz_1,   \
                             g_0_yyyy_0_zzzz_0,   \
                             g_0_yyyy_0_zzzz_1,   \
                             wp_x,                \
                             wp_y,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyyy_0_xxxx_0[i] = 3.0 * g_0_xxyy_0_xxxx_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxx_1[i] * fti_ab_0 + g_0_xxyyy_0_xxxx_0[i] * pb_y +
                                 g_0_xxyyy_0_xxxx_1[i] * wp_y[i];

        g_0_xxyyyy_0_xxxy_0[i] = g_0_yyyy_0_xxxy_0[i] * fi_ab_0 - g_0_yyyy_0_xxxy_1[i] * fti_ab_0 + 3.0 * g_0_xyyyy_0_xxy_1[i] * fi_abcd_0 +
                                 g_0_xyyyy_0_xxxy_0[i] * pb_x + g_0_xyyyy_0_xxxy_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxxz_0[i] = 3.0 * g_0_xxyy_0_xxxz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxz_1[i] * fti_ab_0 + g_0_xxyyy_0_xxxz_0[i] * pb_y +
                                 g_0_xxyyy_0_xxxz_1[i] * wp_y[i];

        g_0_xxyyyy_0_xxyy_0[i] = g_0_yyyy_0_xxyy_0[i] * fi_ab_0 - g_0_yyyy_0_xxyy_1[i] * fti_ab_0 + 2.0 * g_0_xyyyy_0_xyy_1[i] * fi_abcd_0 +
                                 g_0_xyyyy_0_xxyy_0[i] * pb_x + g_0_xyyyy_0_xxyy_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxyz_0[i] = g_0_yyyy_0_xxyz_0[i] * fi_ab_0 - g_0_yyyy_0_xxyz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyy_0_xyz_1[i] * fi_abcd_0 +
                                 g_0_xyyyy_0_xxyz_0[i] * pb_x + g_0_xyyyy_0_xxyz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxzz_0[i] = 3.0 * g_0_xxyy_0_xxzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxzz_1[i] * fti_ab_0 + g_0_xxyyy_0_xxzz_0[i] * pb_y +
                                 g_0_xxyyy_0_xxzz_1[i] * wp_y[i];

        g_0_xxyyyy_0_xyyy_0[i] = g_0_yyyy_0_xyyy_0[i] * fi_ab_0 - g_0_yyyy_0_xyyy_1[i] * fti_ab_0 + g_0_xyyyy_0_yyy_1[i] * fi_abcd_0 +
                                 g_0_xyyyy_0_xyyy_0[i] * pb_x + g_0_xyyyy_0_xyyy_1[i] * wp_x[i];

        g_0_xxyyyy_0_xyyz_0[i] = g_0_yyyy_0_xyyz_0[i] * fi_ab_0 - g_0_yyyy_0_xyyz_1[i] * fti_ab_0 + g_0_xyyyy_0_yyz_1[i] * fi_abcd_0 +
                                 g_0_xyyyy_0_xyyz_0[i] * pb_x + g_0_xyyyy_0_xyyz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xyzz_0[i] = g_0_yyyy_0_xyzz_0[i] * fi_ab_0 - g_0_yyyy_0_xyzz_1[i] * fti_ab_0 + g_0_xyyyy_0_yzz_1[i] * fi_abcd_0 +
                                 g_0_xyyyy_0_xyzz_0[i] * pb_x + g_0_xyyyy_0_xyzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xzzz_0[i] = 3.0 * g_0_xxyy_0_xzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xzzz_1[i] * fti_ab_0 + g_0_xxyyy_0_xzzz_0[i] * pb_y +
                                 g_0_xxyyy_0_xzzz_1[i] * wp_y[i];

        g_0_xxyyyy_0_yyyy_0[i] =
            g_0_yyyy_0_yyyy_0[i] * fi_ab_0 - g_0_yyyy_0_yyyy_1[i] * fti_ab_0 + g_0_xyyyy_0_yyyy_0[i] * pb_x + g_0_xyyyy_0_yyyy_1[i] * wp_x[i];

        g_0_xxyyyy_0_yyyz_0[i] =
            g_0_yyyy_0_yyyz_0[i] * fi_ab_0 - g_0_yyyy_0_yyyz_1[i] * fti_ab_0 + g_0_xyyyy_0_yyyz_0[i] * pb_x + g_0_xyyyy_0_yyyz_1[i] * wp_x[i];

        g_0_xxyyyy_0_yyzz_0[i] =
            g_0_yyyy_0_yyzz_0[i] * fi_ab_0 - g_0_yyyy_0_yyzz_1[i] * fti_ab_0 + g_0_xyyyy_0_yyzz_0[i] * pb_x + g_0_xyyyy_0_yyzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_yzzz_0[i] =
            g_0_yyyy_0_yzzz_0[i] * fi_ab_0 - g_0_yyyy_0_yzzz_1[i] * fti_ab_0 + g_0_xyyyy_0_yzzz_0[i] * pb_x + g_0_xyyyy_0_yzzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_zzzz_0[i] =
            g_0_yyyy_0_zzzz_0[i] * fi_ab_0 - g_0_yyyy_0_zzzz_1[i] * fti_ab_0 + g_0_xyyyy_0_zzzz_0[i] * pb_x + g_0_xyyyy_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 165-180 components of targeted buffer : SISG

    auto g_0_xxyyyz_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 165);

    auto g_0_xxyyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 166);

    auto g_0_xxyyyz_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 167);

    auto g_0_xxyyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 168);

    auto g_0_xxyyyz_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 169);

    auto g_0_xxyyyz_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 170);

    auto g_0_xxyyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 171);

    auto g_0_xxyyyz_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 172);

    auto g_0_xxyyyz_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 173);

    auto g_0_xxyyyz_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 174);

    auto g_0_xxyyyz_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 175);

    auto g_0_xxyyyz_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 176);

    auto g_0_xxyyyz_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 177);

    auto g_0_xxyyyz_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 178);

    auto g_0_xxyyyz_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 179);

#pragma omp simd aligned(g_0_xxyyy_0_xxx_1,       \
                             g_0_xxyyy_0_xxxx_0,  \
                             g_0_xxyyy_0_xxxx_1,  \
                             g_0_xxyyy_0_xxxy_0,  \
                             g_0_xxyyy_0_xxxy_1,  \
                             g_0_xxyyy_0_xxxz_0,  \
                             g_0_xxyyy_0_xxxz_1,  \
                             g_0_xxyyy_0_xxy_1,   \
                             g_0_xxyyy_0_xxyy_0,  \
                             g_0_xxyyy_0_xxyy_1,  \
                             g_0_xxyyy_0_xxyz_0,  \
                             g_0_xxyyy_0_xxyz_1,  \
                             g_0_xxyyy_0_xxz_1,   \
                             g_0_xxyyy_0_xxzz_0,  \
                             g_0_xxyyy_0_xxzz_1,  \
                             g_0_xxyyy_0_xyy_1,   \
                             g_0_xxyyy_0_xyyy_0,  \
                             g_0_xxyyy_0_xyyy_1,  \
                             g_0_xxyyy_0_xyyz_0,  \
                             g_0_xxyyy_0_xyyz_1,  \
                             g_0_xxyyy_0_xyz_1,   \
                             g_0_xxyyy_0_xyzz_0,  \
                             g_0_xxyyy_0_xyzz_1,  \
                             g_0_xxyyy_0_xzz_1,   \
                             g_0_xxyyy_0_xzzz_0,  \
                             g_0_xxyyy_0_xzzz_1,  \
                             g_0_xxyyy_0_yyy_1,   \
                             g_0_xxyyy_0_yyyy_0,  \
                             g_0_xxyyy_0_yyyy_1,  \
                             g_0_xxyyy_0_yyyz_0,  \
                             g_0_xxyyy_0_yyyz_1,  \
                             g_0_xxyyy_0_yyz_1,   \
                             g_0_xxyyy_0_yyzz_0,  \
                             g_0_xxyyy_0_yyzz_1,  \
                             g_0_xxyyy_0_yzz_1,   \
                             g_0_xxyyy_0_yzzz_0,  \
                             g_0_xxyyy_0_yzzz_1,  \
                             g_0_xxyyy_0_zzz_1,   \
                             g_0_xxyyy_0_zzzz_0,  \
                             g_0_xxyyy_0_zzzz_1,  \
                             g_0_xxyyyz_0_xxxx_0, \
                             g_0_xxyyyz_0_xxxy_0, \
                             g_0_xxyyyz_0_xxxz_0, \
                             g_0_xxyyyz_0_xxyy_0, \
                             g_0_xxyyyz_0_xxyz_0, \
                             g_0_xxyyyz_0_xxzz_0, \
                             g_0_xxyyyz_0_xyyy_0, \
                             g_0_xxyyyz_0_xyyz_0, \
                             g_0_xxyyyz_0_xyzz_0, \
                             g_0_xxyyyz_0_xzzz_0, \
                             g_0_xxyyyz_0_yyyy_0, \
                             g_0_xxyyyz_0_yyyz_0, \
                             g_0_xxyyyz_0_yyzz_0, \
                             g_0_xxyyyz_0_yzzz_0, \
                             g_0_xxyyyz_0_zzzz_0, \
                             wp_z,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyyz_0_xxxx_0[i] = g_0_xxyyy_0_xxxx_0[i] * pb_z + g_0_xxyyy_0_xxxx_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxy_0[i] = g_0_xxyyy_0_xxxy_0[i] * pb_z + g_0_xxyyy_0_xxxy_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxz_0[i] = g_0_xxyyy_0_xxx_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxz_0[i] * pb_z + g_0_xxyyy_0_xxxz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxyy_0[i] = g_0_xxyyy_0_xxyy_0[i] * pb_z + g_0_xxyyy_0_xxyy_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxyz_0[i] = g_0_xxyyy_0_xxy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxyz_0[i] * pb_z + g_0_xxyyy_0_xxyz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxzz_0[i] = 2.0 * g_0_xxyyy_0_xxz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxzz_0[i] * pb_z + g_0_xxyyy_0_xxzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xyyy_0[i] = g_0_xxyyy_0_xyyy_0[i] * pb_z + g_0_xxyyy_0_xyyy_1[i] * wp_z[i];

        g_0_xxyyyz_0_xyyz_0[i] = g_0_xxyyy_0_xyy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyyz_0[i] * pb_z + g_0_xxyyy_0_xyyz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xyzz_0[i] = 2.0 * g_0_xxyyy_0_xyz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyzz_0[i] * pb_z + g_0_xxyyy_0_xyzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xzzz_0[i] = 3.0 * g_0_xxyyy_0_xzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xzzz_0[i] * pb_z + g_0_xxyyy_0_xzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_yyyy_0[i] = g_0_xxyyy_0_yyyy_0[i] * pb_z + g_0_xxyyy_0_yyyy_1[i] * wp_z[i];

        g_0_xxyyyz_0_yyyz_0[i] = g_0_xxyyy_0_yyy_1[i] * fi_abcd_0 + g_0_xxyyy_0_yyyz_0[i] * pb_z + g_0_xxyyy_0_yyyz_1[i] * wp_z[i];

        g_0_xxyyyz_0_yyzz_0[i] = 2.0 * g_0_xxyyy_0_yyz_1[i] * fi_abcd_0 + g_0_xxyyy_0_yyzz_0[i] * pb_z + g_0_xxyyy_0_yyzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_yzzz_0[i] = 3.0 * g_0_xxyyy_0_yzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_yzzz_0[i] * pb_z + g_0_xxyyy_0_yzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_zzzz_0[i] = 4.0 * g_0_xxyyy_0_zzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_zzzz_0[i] * pb_z + g_0_xxyyy_0_zzzz_1[i] * wp_z[i];
    }

    /// Set up 180-195 components of targeted buffer : SISG

    auto g_0_xxyyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 180);

    auto g_0_xxyyzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 181);

    auto g_0_xxyyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 182);

    auto g_0_xxyyzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 183);

    auto g_0_xxyyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 184);

    auto g_0_xxyyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 185);

    auto g_0_xxyyzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 186);

    auto g_0_xxyyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 187);

    auto g_0_xxyyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 188);

    auto g_0_xxyyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 189);

    auto g_0_xxyyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 190);

    auto g_0_xxyyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 191);

    auto g_0_xxyyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 192);

    auto g_0_xxyyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 193);

    auto g_0_xxyyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 194);

#pragma omp simd aligned(g_0_xxyy_0_xxxy_0,       \
                             g_0_xxyy_0_xxxy_1,   \
                             g_0_xxyy_0_xxyy_0,   \
                             g_0_xxyy_0_xxyy_1,   \
                             g_0_xxyy_0_xyyy_0,   \
                             g_0_xxyy_0_xyyy_1,   \
                             g_0_xxyyz_0_xxxy_0,  \
                             g_0_xxyyz_0_xxxy_1,  \
                             g_0_xxyyz_0_xxyy_0,  \
                             g_0_xxyyz_0_xxyy_1,  \
                             g_0_xxyyz_0_xyyy_0,  \
                             g_0_xxyyz_0_xyyy_1,  \
                             g_0_xxyyzz_0_xxxx_0, \
                             g_0_xxyyzz_0_xxxy_0, \
                             g_0_xxyyzz_0_xxxz_0, \
                             g_0_xxyyzz_0_xxyy_0, \
                             g_0_xxyyzz_0_xxyz_0, \
                             g_0_xxyyzz_0_xxzz_0, \
                             g_0_xxyyzz_0_xyyy_0, \
                             g_0_xxyyzz_0_xyyz_0, \
                             g_0_xxyyzz_0_xyzz_0, \
                             g_0_xxyyzz_0_xzzz_0, \
                             g_0_xxyyzz_0_yyyy_0, \
                             g_0_xxyyzz_0_yyyz_0, \
                             g_0_xxyyzz_0_yyzz_0, \
                             g_0_xxyyzz_0_yzzz_0, \
                             g_0_xxyyzz_0_zzzz_0, \
                             g_0_xxyzz_0_xxxx_0,  \
                             g_0_xxyzz_0_xxxx_1,  \
                             g_0_xxyzz_0_xxxz_0,  \
                             g_0_xxyzz_0_xxxz_1,  \
                             g_0_xxyzz_0_xxzz_0,  \
                             g_0_xxyzz_0_xxzz_1,  \
                             g_0_xxyzz_0_xzzz_0,  \
                             g_0_xxyzz_0_xzzz_1,  \
                             g_0_xxzz_0_xxxx_0,   \
                             g_0_xxzz_0_xxxx_1,   \
                             g_0_xxzz_0_xxxz_0,   \
                             g_0_xxzz_0_xxxz_1,   \
                             g_0_xxzz_0_xxzz_0,   \
                             g_0_xxzz_0_xxzz_1,   \
                             g_0_xxzz_0_xzzz_0,   \
                             g_0_xxzz_0_xzzz_1,   \
                             g_0_xyyzz_0_xxyz_0,  \
                             g_0_xyyzz_0_xxyz_1,  \
                             g_0_xyyzz_0_xyyz_0,  \
                             g_0_xyyzz_0_xyyz_1,  \
                             g_0_xyyzz_0_xyz_1,   \
                             g_0_xyyzz_0_xyzz_0,  \
                             g_0_xyyzz_0_xyzz_1,  \
                             g_0_xyyzz_0_yyyy_0,  \
                             g_0_xyyzz_0_yyyy_1,  \
                             g_0_xyyzz_0_yyyz_0,  \
                             g_0_xyyzz_0_yyyz_1,  \
                             g_0_xyyzz_0_yyz_1,   \
                             g_0_xyyzz_0_yyzz_0,  \
                             g_0_xyyzz_0_yyzz_1,  \
                             g_0_xyyzz_0_yzz_1,   \
                             g_0_xyyzz_0_yzzz_0,  \
                             g_0_xyyzz_0_yzzz_1,  \
                             g_0_xyyzz_0_zzzz_0,  \
                             g_0_xyyzz_0_zzzz_1,  \
                             g_0_yyzz_0_xxyz_0,   \
                             g_0_yyzz_0_xxyz_1,   \
                             g_0_yyzz_0_xyyz_0,   \
                             g_0_yyzz_0_xyyz_1,   \
                             g_0_yyzz_0_xyzz_0,   \
                             g_0_yyzz_0_xyzz_1,   \
                             g_0_yyzz_0_yyyy_0,   \
                             g_0_yyzz_0_yyyy_1,   \
                             g_0_yyzz_0_yyyz_0,   \
                             g_0_yyzz_0_yyyz_1,   \
                             g_0_yyzz_0_yyzz_0,   \
                             g_0_yyzz_0_yyzz_1,   \
                             g_0_yyzz_0_yzzz_0,   \
                             g_0_yyzz_0_yzzz_1,   \
                             g_0_yyzz_0_zzzz_0,   \
                             g_0_yyzz_0_zzzz_1,   \
                             wp_x,                \
                             wp_y,                \
                             wp_z,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyzz_0_xxxx_0[i] =
            g_0_xxzz_0_xxxx_0[i] * fi_ab_0 - g_0_xxzz_0_xxxx_1[i] * fti_ab_0 + g_0_xxyzz_0_xxxx_0[i] * pb_y + g_0_xxyzz_0_xxxx_1[i] * wp_y[i];

        g_0_xxyyzz_0_xxxy_0[i] =
            g_0_xxyy_0_xxxy_0[i] * fi_ab_0 - g_0_xxyy_0_xxxy_1[i] * fti_ab_0 + g_0_xxyyz_0_xxxy_0[i] * pb_z + g_0_xxyyz_0_xxxy_1[i] * wp_z[i];

        g_0_xxyyzz_0_xxxz_0[i] =
            g_0_xxzz_0_xxxz_0[i] * fi_ab_0 - g_0_xxzz_0_xxxz_1[i] * fti_ab_0 + g_0_xxyzz_0_xxxz_0[i] * pb_y + g_0_xxyzz_0_xxxz_1[i] * wp_y[i];

        g_0_xxyyzz_0_xxyy_0[i] =
            g_0_xxyy_0_xxyy_0[i] * fi_ab_0 - g_0_xxyy_0_xxyy_1[i] * fti_ab_0 + g_0_xxyyz_0_xxyy_0[i] * pb_z + g_0_xxyyz_0_xxyy_1[i] * wp_z[i];

        g_0_xxyyzz_0_xxyz_0[i] = g_0_yyzz_0_xxyz_0[i] * fi_ab_0 - g_0_yyzz_0_xxyz_1[i] * fti_ab_0 + 2.0 * g_0_xyyzz_0_xyz_1[i] * fi_abcd_0 +
                                 g_0_xyyzz_0_xxyz_0[i] * pb_x + g_0_xyyzz_0_xxyz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xxzz_0[i] =
            g_0_xxzz_0_xxzz_0[i] * fi_ab_0 - g_0_xxzz_0_xxzz_1[i] * fti_ab_0 + g_0_xxyzz_0_xxzz_0[i] * pb_y + g_0_xxyzz_0_xxzz_1[i] * wp_y[i];

        g_0_xxyyzz_0_xyyy_0[i] =
            g_0_xxyy_0_xyyy_0[i] * fi_ab_0 - g_0_xxyy_0_xyyy_1[i] * fti_ab_0 + g_0_xxyyz_0_xyyy_0[i] * pb_z + g_0_xxyyz_0_xyyy_1[i] * wp_z[i];

        g_0_xxyyzz_0_xyyz_0[i] = g_0_yyzz_0_xyyz_0[i] * fi_ab_0 - g_0_yyzz_0_xyyz_1[i] * fti_ab_0 + g_0_xyyzz_0_yyz_1[i] * fi_abcd_0 +
                                 g_0_xyyzz_0_xyyz_0[i] * pb_x + g_0_xyyzz_0_xyyz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xyzz_0[i] = g_0_yyzz_0_xyzz_0[i] * fi_ab_0 - g_0_yyzz_0_xyzz_1[i] * fti_ab_0 + g_0_xyyzz_0_yzz_1[i] * fi_abcd_0 +
                                 g_0_xyyzz_0_xyzz_0[i] * pb_x + g_0_xyyzz_0_xyzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xzzz_0[i] =
            g_0_xxzz_0_xzzz_0[i] * fi_ab_0 - g_0_xxzz_0_xzzz_1[i] * fti_ab_0 + g_0_xxyzz_0_xzzz_0[i] * pb_y + g_0_xxyzz_0_xzzz_1[i] * wp_y[i];

        g_0_xxyyzz_0_yyyy_0[i] =
            g_0_yyzz_0_yyyy_0[i] * fi_ab_0 - g_0_yyzz_0_yyyy_1[i] * fti_ab_0 + g_0_xyyzz_0_yyyy_0[i] * pb_x + g_0_xyyzz_0_yyyy_1[i] * wp_x[i];

        g_0_xxyyzz_0_yyyz_0[i] =
            g_0_yyzz_0_yyyz_0[i] * fi_ab_0 - g_0_yyzz_0_yyyz_1[i] * fti_ab_0 + g_0_xyyzz_0_yyyz_0[i] * pb_x + g_0_xyyzz_0_yyyz_1[i] * wp_x[i];

        g_0_xxyyzz_0_yyzz_0[i] =
            g_0_yyzz_0_yyzz_0[i] * fi_ab_0 - g_0_yyzz_0_yyzz_1[i] * fti_ab_0 + g_0_xyyzz_0_yyzz_0[i] * pb_x + g_0_xyyzz_0_yyzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_yzzz_0[i] =
            g_0_yyzz_0_yzzz_0[i] * fi_ab_0 - g_0_yyzz_0_yzzz_1[i] * fti_ab_0 + g_0_xyyzz_0_yzzz_0[i] * pb_x + g_0_xyyzz_0_yzzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_zzzz_0[i] =
            g_0_yyzz_0_zzzz_0[i] * fi_ab_0 - g_0_yyzz_0_zzzz_1[i] * fti_ab_0 + g_0_xyyzz_0_zzzz_0[i] * pb_x + g_0_xyyzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 195-210 components of targeted buffer : SISG

    auto g_0_xxyzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 195);

    auto g_0_xxyzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 196);

    auto g_0_xxyzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 197);

    auto g_0_xxyzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 198);

    auto g_0_xxyzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 199);

    auto g_0_xxyzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 200);

    auto g_0_xxyzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 201);

    auto g_0_xxyzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 202);

    auto g_0_xxyzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 203);

    auto g_0_xxyzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 204);

    auto g_0_xxyzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 205);

    auto g_0_xxyzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 206);

    auto g_0_xxyzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 207);

    auto g_0_xxyzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 208);

    auto g_0_xxyzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 209);

#pragma omp simd aligned(g_0_xxyzzz_0_xxxx_0,     \
                             g_0_xxyzzz_0_xxxy_0, \
                             g_0_xxyzzz_0_xxxz_0, \
                             g_0_xxyzzz_0_xxyy_0, \
                             g_0_xxyzzz_0_xxyz_0, \
                             g_0_xxyzzz_0_xxzz_0, \
                             g_0_xxyzzz_0_xyyy_0, \
                             g_0_xxyzzz_0_xyyz_0, \
                             g_0_xxyzzz_0_xyzz_0, \
                             g_0_xxyzzz_0_xzzz_0, \
                             g_0_xxyzzz_0_yyyy_0, \
                             g_0_xxyzzz_0_yyyz_0, \
                             g_0_xxyzzz_0_yyzz_0, \
                             g_0_xxyzzz_0_yzzz_0, \
                             g_0_xxyzzz_0_zzzz_0, \
                             g_0_xxzzz_0_xxx_1,   \
                             g_0_xxzzz_0_xxxx_0,  \
                             g_0_xxzzz_0_xxxx_1,  \
                             g_0_xxzzz_0_xxxy_0,  \
                             g_0_xxzzz_0_xxxy_1,  \
                             g_0_xxzzz_0_xxxz_0,  \
                             g_0_xxzzz_0_xxxz_1,  \
                             g_0_xxzzz_0_xxy_1,   \
                             g_0_xxzzz_0_xxyy_0,  \
                             g_0_xxzzz_0_xxyy_1,  \
                             g_0_xxzzz_0_xxyz_0,  \
                             g_0_xxzzz_0_xxyz_1,  \
                             g_0_xxzzz_0_xxz_1,   \
                             g_0_xxzzz_0_xxzz_0,  \
                             g_0_xxzzz_0_xxzz_1,  \
                             g_0_xxzzz_0_xyy_1,   \
                             g_0_xxzzz_0_xyyy_0,  \
                             g_0_xxzzz_0_xyyy_1,  \
                             g_0_xxzzz_0_xyyz_0,  \
                             g_0_xxzzz_0_xyyz_1,  \
                             g_0_xxzzz_0_xyz_1,   \
                             g_0_xxzzz_0_xyzz_0,  \
                             g_0_xxzzz_0_xyzz_1,  \
                             g_0_xxzzz_0_xzz_1,   \
                             g_0_xxzzz_0_xzzz_0,  \
                             g_0_xxzzz_0_xzzz_1,  \
                             g_0_xxzzz_0_yyy_1,   \
                             g_0_xxzzz_0_yyyy_0,  \
                             g_0_xxzzz_0_yyyy_1,  \
                             g_0_xxzzz_0_yyyz_0,  \
                             g_0_xxzzz_0_yyyz_1,  \
                             g_0_xxzzz_0_yyz_1,   \
                             g_0_xxzzz_0_yyzz_0,  \
                             g_0_xxzzz_0_yyzz_1,  \
                             g_0_xxzzz_0_yzz_1,   \
                             g_0_xxzzz_0_yzzz_0,  \
                             g_0_xxzzz_0_yzzz_1,  \
                             g_0_xxzzz_0_zzz_1,   \
                             g_0_xxzzz_0_zzzz_0,  \
                             g_0_xxzzz_0_zzzz_1,  \
                             wp_y,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyzzz_0_xxxx_0[i] = g_0_xxzzz_0_xxxx_0[i] * pb_y + g_0_xxzzz_0_xxxx_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxy_0[i] = g_0_xxzzz_0_xxx_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxy_0[i] * pb_y + g_0_xxzzz_0_xxxy_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxz_0[i] = g_0_xxzzz_0_xxxz_0[i] * pb_y + g_0_xxzzz_0_xxxz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxyy_0[i] = 2.0 * g_0_xxzzz_0_xxy_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxyy_0[i] * pb_y + g_0_xxzzz_0_xxyy_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxyz_0[i] = g_0_xxzzz_0_xxz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxyz_0[i] * pb_y + g_0_xxzzz_0_xxyz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxzz_0[i] = g_0_xxzzz_0_xxzz_0[i] * pb_y + g_0_xxzzz_0_xxzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xyyy_0[i] = 3.0 * g_0_xxzzz_0_xyy_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyyy_0[i] * pb_y + g_0_xxzzz_0_xyyy_1[i] * wp_y[i];

        g_0_xxyzzz_0_xyyz_0[i] = 2.0 * g_0_xxzzz_0_xyz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyyz_0[i] * pb_y + g_0_xxzzz_0_xyyz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xyzz_0[i] = g_0_xxzzz_0_xzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyzz_0[i] * pb_y + g_0_xxzzz_0_xyzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xzzz_0[i] = g_0_xxzzz_0_xzzz_0[i] * pb_y + g_0_xxzzz_0_xzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_yyyy_0[i] = 4.0 * g_0_xxzzz_0_yyy_1[i] * fi_abcd_0 + g_0_xxzzz_0_yyyy_0[i] * pb_y + g_0_xxzzz_0_yyyy_1[i] * wp_y[i];

        g_0_xxyzzz_0_yyyz_0[i] = 3.0 * g_0_xxzzz_0_yyz_1[i] * fi_abcd_0 + g_0_xxzzz_0_yyyz_0[i] * pb_y + g_0_xxzzz_0_yyyz_1[i] * wp_y[i];

        g_0_xxyzzz_0_yyzz_0[i] = 2.0 * g_0_xxzzz_0_yzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_yyzz_0[i] * pb_y + g_0_xxzzz_0_yyzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_yzzz_0[i] = g_0_xxzzz_0_zzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_yzzz_0[i] * pb_y + g_0_xxzzz_0_yzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_zzzz_0[i] = g_0_xxzzz_0_zzzz_0[i] * pb_y + g_0_xxzzz_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 210-225 components of targeted buffer : SISG

    auto g_0_xxzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 210);

    auto g_0_xxzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 211);

    auto g_0_xxzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 212);

    auto g_0_xxzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 213);

    auto g_0_xxzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 214);

    auto g_0_xxzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 215);

    auto g_0_xxzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 216);

    auto g_0_xxzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 217);

    auto g_0_xxzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 218);

    auto g_0_xxzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 219);

    auto g_0_xxzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 220);

    auto g_0_xxzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 221);

    auto g_0_xxzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 222);

    auto g_0_xxzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 223);

    auto g_0_xxzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 224);

#pragma omp simd aligned(g_0_xxzz_0_xxxx_0,       \
                             g_0_xxzz_0_xxxx_1,   \
                             g_0_xxzz_0_xxxy_0,   \
                             g_0_xxzz_0_xxxy_1,   \
                             g_0_xxzz_0_xxyy_0,   \
                             g_0_xxzz_0_xxyy_1,   \
                             g_0_xxzz_0_xyyy_0,   \
                             g_0_xxzz_0_xyyy_1,   \
                             g_0_xxzzz_0_xxxx_0,  \
                             g_0_xxzzz_0_xxxx_1,  \
                             g_0_xxzzz_0_xxxy_0,  \
                             g_0_xxzzz_0_xxxy_1,  \
                             g_0_xxzzz_0_xxyy_0,  \
                             g_0_xxzzz_0_xxyy_1,  \
                             g_0_xxzzz_0_xyyy_0,  \
                             g_0_xxzzz_0_xyyy_1,  \
                             g_0_xxzzzz_0_xxxx_0, \
                             g_0_xxzzzz_0_xxxy_0, \
                             g_0_xxzzzz_0_xxxz_0, \
                             g_0_xxzzzz_0_xxyy_0, \
                             g_0_xxzzzz_0_xxyz_0, \
                             g_0_xxzzzz_0_xxzz_0, \
                             g_0_xxzzzz_0_xyyy_0, \
                             g_0_xxzzzz_0_xyyz_0, \
                             g_0_xxzzzz_0_xyzz_0, \
                             g_0_xxzzzz_0_xzzz_0, \
                             g_0_xxzzzz_0_yyyy_0, \
                             g_0_xxzzzz_0_yyyz_0, \
                             g_0_xxzzzz_0_yyzz_0, \
                             g_0_xxzzzz_0_yzzz_0, \
                             g_0_xxzzzz_0_zzzz_0, \
                             g_0_xzzzz_0_xxxz_0,  \
                             g_0_xzzzz_0_xxxz_1,  \
                             g_0_xzzzz_0_xxyz_0,  \
                             g_0_xzzzz_0_xxyz_1,  \
                             g_0_xzzzz_0_xxz_1,   \
                             g_0_xzzzz_0_xxzz_0,  \
                             g_0_xzzzz_0_xxzz_1,  \
                             g_0_xzzzz_0_xyyz_0,  \
                             g_0_xzzzz_0_xyyz_1,  \
                             g_0_xzzzz_0_xyz_1,   \
                             g_0_xzzzz_0_xyzz_0,  \
                             g_0_xzzzz_0_xyzz_1,  \
                             g_0_xzzzz_0_xzz_1,   \
                             g_0_xzzzz_0_xzzz_0,  \
                             g_0_xzzzz_0_xzzz_1,  \
                             g_0_xzzzz_0_yyyy_0,  \
                             g_0_xzzzz_0_yyyy_1,  \
                             g_0_xzzzz_0_yyyz_0,  \
                             g_0_xzzzz_0_yyyz_1,  \
                             g_0_xzzzz_0_yyz_1,   \
                             g_0_xzzzz_0_yyzz_0,  \
                             g_0_xzzzz_0_yyzz_1,  \
                             g_0_xzzzz_0_yzz_1,   \
                             g_0_xzzzz_0_yzzz_0,  \
                             g_0_xzzzz_0_yzzz_1,  \
                             g_0_xzzzz_0_zzz_1,   \
                             g_0_xzzzz_0_zzzz_0,  \
                             g_0_xzzzz_0_zzzz_1,  \
                             g_0_zzzz_0_xxxz_0,   \
                             g_0_zzzz_0_xxxz_1,   \
                             g_0_zzzz_0_xxyz_0,   \
                             g_0_zzzz_0_xxyz_1,   \
                             g_0_zzzz_0_xxzz_0,   \
                             g_0_zzzz_0_xxzz_1,   \
                             g_0_zzzz_0_xyyz_0,   \
                             g_0_zzzz_0_xyyz_1,   \
                             g_0_zzzz_0_xyzz_0,   \
                             g_0_zzzz_0_xyzz_1,   \
                             g_0_zzzz_0_xzzz_0,   \
                             g_0_zzzz_0_xzzz_1,   \
                             g_0_zzzz_0_yyyy_0,   \
                             g_0_zzzz_0_yyyy_1,   \
                             g_0_zzzz_0_yyyz_0,   \
                             g_0_zzzz_0_yyyz_1,   \
                             g_0_zzzz_0_yyzz_0,   \
                             g_0_zzzz_0_yyzz_1,   \
                             g_0_zzzz_0_yzzz_0,   \
                             g_0_zzzz_0_yzzz_1,   \
                             g_0_zzzz_0_zzzz_0,   \
                             g_0_zzzz_0_zzzz_1,   \
                             wp_x,                \
                             wp_z,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxzzzz_0_xxxx_0[i] = 3.0 * g_0_xxzz_0_xxxx_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxx_1[i] * fti_ab_0 + g_0_xxzzz_0_xxxx_0[i] * pb_z +
                                 g_0_xxzzz_0_xxxx_1[i] * wp_z[i];

        g_0_xxzzzz_0_xxxy_0[i] = 3.0 * g_0_xxzz_0_xxxy_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxy_1[i] * fti_ab_0 + g_0_xxzzz_0_xxxy_0[i] * pb_z +
                                 g_0_xxzzz_0_xxxy_1[i] * wp_z[i];

        g_0_xxzzzz_0_xxxz_0[i] = g_0_zzzz_0_xxxz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxz_1[i] * fti_ab_0 + 3.0 * g_0_xzzzz_0_xxz_1[i] * fi_abcd_0 +
                                 g_0_xzzzz_0_xxxz_0[i] * pb_x + g_0_xzzzz_0_xxxz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxyy_0[i] = 3.0 * g_0_xxzz_0_xxyy_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxyy_1[i] * fti_ab_0 + g_0_xxzzz_0_xxyy_0[i] * pb_z +
                                 g_0_xxzzz_0_xxyy_1[i] * wp_z[i];

        g_0_xxzzzz_0_xxyz_0[i] = g_0_zzzz_0_xxyz_0[i] * fi_ab_0 - g_0_zzzz_0_xxyz_1[i] * fti_ab_0 + 2.0 * g_0_xzzzz_0_xyz_1[i] * fi_abcd_0 +
                                 g_0_xzzzz_0_xxyz_0[i] * pb_x + g_0_xzzzz_0_xxyz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxzz_0[i] = g_0_zzzz_0_xxzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxzz_1[i] * fti_ab_0 + 2.0 * g_0_xzzzz_0_xzz_1[i] * fi_abcd_0 +
                                 g_0_xzzzz_0_xxzz_0[i] * pb_x + g_0_xzzzz_0_xxzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xyyy_0[i] = 3.0 * g_0_xxzz_0_xyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xyyy_1[i] * fti_ab_0 + g_0_xxzzz_0_xyyy_0[i] * pb_z +
                                 g_0_xxzzz_0_xyyy_1[i] * wp_z[i];

        g_0_xxzzzz_0_xyyz_0[i] = g_0_zzzz_0_xyyz_0[i] * fi_ab_0 - g_0_zzzz_0_xyyz_1[i] * fti_ab_0 + g_0_xzzzz_0_yyz_1[i] * fi_abcd_0 +
                                 g_0_xzzzz_0_xyyz_0[i] * pb_x + g_0_xzzzz_0_xyyz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xyzz_0[i] = g_0_zzzz_0_xyzz_0[i] * fi_ab_0 - g_0_zzzz_0_xyzz_1[i] * fti_ab_0 + g_0_xzzzz_0_yzz_1[i] * fi_abcd_0 +
                                 g_0_xzzzz_0_xyzz_0[i] * pb_x + g_0_xzzzz_0_xyzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xzzz_0[i] = g_0_zzzz_0_xzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xzzz_1[i] * fti_ab_0 + g_0_xzzzz_0_zzz_1[i] * fi_abcd_0 +
                                 g_0_xzzzz_0_xzzz_0[i] * pb_x + g_0_xzzzz_0_xzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_yyyy_0[i] =
            g_0_zzzz_0_yyyy_0[i] * fi_ab_0 - g_0_zzzz_0_yyyy_1[i] * fti_ab_0 + g_0_xzzzz_0_yyyy_0[i] * pb_x + g_0_xzzzz_0_yyyy_1[i] * wp_x[i];

        g_0_xxzzzz_0_yyyz_0[i] =
            g_0_zzzz_0_yyyz_0[i] * fi_ab_0 - g_0_zzzz_0_yyyz_1[i] * fti_ab_0 + g_0_xzzzz_0_yyyz_0[i] * pb_x + g_0_xzzzz_0_yyyz_1[i] * wp_x[i];

        g_0_xxzzzz_0_yyzz_0[i] =
            g_0_zzzz_0_yyzz_0[i] * fi_ab_0 - g_0_zzzz_0_yyzz_1[i] * fti_ab_0 + g_0_xzzzz_0_yyzz_0[i] * pb_x + g_0_xzzzz_0_yyzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_yzzz_0[i] =
            g_0_zzzz_0_yzzz_0[i] * fi_ab_0 - g_0_zzzz_0_yzzz_1[i] * fti_ab_0 + g_0_xzzzz_0_yzzz_0[i] * pb_x + g_0_xzzzz_0_yzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_zzzz_0[i] =
            g_0_zzzz_0_zzzz_0[i] * fi_ab_0 - g_0_zzzz_0_zzzz_1[i] * fti_ab_0 + g_0_xzzzz_0_zzzz_0[i] * pb_x + g_0_xzzzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 225-240 components of targeted buffer : SISG

    auto g_0_xyyyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 225);

    auto g_0_xyyyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 226);

    auto g_0_xyyyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 227);

    auto g_0_xyyyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 228);

    auto g_0_xyyyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 229);

    auto g_0_xyyyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 230);

    auto g_0_xyyyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 231);

    auto g_0_xyyyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 232);

    auto g_0_xyyyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 233);

    auto g_0_xyyyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 234);

    auto g_0_xyyyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 235);

    auto g_0_xyyyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 236);

    auto g_0_xyyyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 237);

    auto g_0_xyyyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 238);

    auto g_0_xyyyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 239);

#pragma omp simd aligned(g_0_xyyyyy_0_xxxx_0,     \
                             g_0_xyyyyy_0_xxxy_0, \
                             g_0_xyyyyy_0_xxxz_0, \
                             g_0_xyyyyy_0_xxyy_0, \
                             g_0_xyyyyy_0_xxyz_0, \
                             g_0_xyyyyy_0_xxzz_0, \
                             g_0_xyyyyy_0_xyyy_0, \
                             g_0_xyyyyy_0_xyyz_0, \
                             g_0_xyyyyy_0_xyzz_0, \
                             g_0_xyyyyy_0_xzzz_0, \
                             g_0_xyyyyy_0_yyyy_0, \
                             g_0_xyyyyy_0_yyyz_0, \
                             g_0_xyyyyy_0_yyzz_0, \
                             g_0_xyyyyy_0_yzzz_0, \
                             g_0_xyyyyy_0_zzzz_0, \
                             g_0_yyyyy_0_xxx_1,   \
                             g_0_yyyyy_0_xxxx_0,  \
                             g_0_yyyyy_0_xxxx_1,  \
                             g_0_yyyyy_0_xxxy_0,  \
                             g_0_yyyyy_0_xxxy_1,  \
                             g_0_yyyyy_0_xxxz_0,  \
                             g_0_yyyyy_0_xxxz_1,  \
                             g_0_yyyyy_0_xxy_1,   \
                             g_0_yyyyy_0_xxyy_0,  \
                             g_0_yyyyy_0_xxyy_1,  \
                             g_0_yyyyy_0_xxyz_0,  \
                             g_0_yyyyy_0_xxyz_1,  \
                             g_0_yyyyy_0_xxz_1,   \
                             g_0_yyyyy_0_xxzz_0,  \
                             g_0_yyyyy_0_xxzz_1,  \
                             g_0_yyyyy_0_xyy_1,   \
                             g_0_yyyyy_0_xyyy_0,  \
                             g_0_yyyyy_0_xyyy_1,  \
                             g_0_yyyyy_0_xyyz_0,  \
                             g_0_yyyyy_0_xyyz_1,  \
                             g_0_yyyyy_0_xyz_1,   \
                             g_0_yyyyy_0_xyzz_0,  \
                             g_0_yyyyy_0_xyzz_1,  \
                             g_0_yyyyy_0_xzz_1,   \
                             g_0_yyyyy_0_xzzz_0,  \
                             g_0_yyyyy_0_xzzz_1,  \
                             g_0_yyyyy_0_yyy_1,   \
                             g_0_yyyyy_0_yyyy_0,  \
                             g_0_yyyyy_0_yyyy_1,  \
                             g_0_yyyyy_0_yyyz_0,  \
                             g_0_yyyyy_0_yyyz_1,  \
                             g_0_yyyyy_0_yyz_1,   \
                             g_0_yyyyy_0_yyzz_0,  \
                             g_0_yyyyy_0_yyzz_1,  \
                             g_0_yyyyy_0_yzz_1,   \
                             g_0_yyyyy_0_yzzz_0,  \
                             g_0_yyyyy_0_yzzz_1,  \
                             g_0_yyyyy_0_zzz_1,   \
                             g_0_yyyyy_0_zzzz_0,  \
                             g_0_yyyyy_0_zzzz_1,  \
                             wp_x,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyy_0_xxxx_0[i] = 4.0 * g_0_yyyyy_0_xxx_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxx_0[i] * pb_x + g_0_yyyyy_0_xxxx_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxy_0[i] = 3.0 * g_0_yyyyy_0_xxy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxy_0[i] * pb_x + g_0_yyyyy_0_xxxy_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxz_0[i] = 3.0 * g_0_yyyyy_0_xxz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxz_0[i] * pb_x + g_0_yyyyy_0_xxxz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxyy_0[i] = 2.0 * g_0_yyyyy_0_xyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyy_0[i] * pb_x + g_0_yyyyy_0_xxyy_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxyz_0[i] = 2.0 * g_0_yyyyy_0_xyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyz_0[i] * pb_x + g_0_yyyyy_0_xxyz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxzz_0[i] = 2.0 * g_0_yyyyy_0_xzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxzz_0[i] * pb_x + g_0_yyyyy_0_xxzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xyyy_0[i] = g_0_yyyyy_0_yyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyy_0[i] * pb_x + g_0_yyyyy_0_xyyy_1[i] * wp_x[i];

        g_0_xyyyyy_0_xyyz_0[i] = g_0_yyyyy_0_yyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyz_0[i] * pb_x + g_0_yyyyy_0_xyyz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xyzz_0[i] = g_0_yyyyy_0_yzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyzz_0[i] * pb_x + g_0_yyyyy_0_xyzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xzzz_0[i] = g_0_yyyyy_0_zzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xzzz_0[i] * pb_x + g_0_yyyyy_0_xzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_yyyy_0[i] = g_0_yyyyy_0_yyyy_0[i] * pb_x + g_0_yyyyy_0_yyyy_1[i] * wp_x[i];

        g_0_xyyyyy_0_yyyz_0[i] = g_0_yyyyy_0_yyyz_0[i] * pb_x + g_0_yyyyy_0_yyyz_1[i] * wp_x[i];

        g_0_xyyyyy_0_yyzz_0[i] = g_0_yyyyy_0_yyzz_0[i] * pb_x + g_0_yyyyy_0_yyzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_yzzz_0[i] = g_0_yyyyy_0_yzzz_0[i] * pb_x + g_0_yyyyy_0_yzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_zzzz_0[i] = g_0_yyyyy_0_zzzz_0[i] * pb_x + g_0_yyyyy_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 240-255 components of targeted buffer : SISG

    auto g_0_xyyyyz_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 240);

    auto g_0_xyyyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 241);

    auto g_0_xyyyyz_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 242);

    auto g_0_xyyyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 243);

    auto g_0_xyyyyz_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 244);

    auto g_0_xyyyyz_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 245);

    auto g_0_xyyyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 246);

    auto g_0_xyyyyz_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 247);

    auto g_0_xyyyyz_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 248);

    auto g_0_xyyyyz_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 249);

    auto g_0_xyyyyz_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 250);

    auto g_0_xyyyyz_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 251);

    auto g_0_xyyyyz_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 252);

    auto g_0_xyyyyz_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 253);

    auto g_0_xyyyyz_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 254);

#pragma omp simd aligned(g_0_xyyyy_0_xxxx_0,      \
                             g_0_xyyyy_0_xxxx_1,  \
                             g_0_xyyyy_0_xxxy_0,  \
                             g_0_xyyyy_0_xxxy_1,  \
                             g_0_xyyyy_0_xxyy_0,  \
                             g_0_xyyyy_0_xxyy_1,  \
                             g_0_xyyyy_0_xyyy_0,  \
                             g_0_xyyyy_0_xyyy_1,  \
                             g_0_xyyyyz_0_xxxx_0, \
                             g_0_xyyyyz_0_xxxy_0, \
                             g_0_xyyyyz_0_xxxz_0, \
                             g_0_xyyyyz_0_xxyy_0, \
                             g_0_xyyyyz_0_xxyz_0, \
                             g_0_xyyyyz_0_xxzz_0, \
                             g_0_xyyyyz_0_xyyy_0, \
                             g_0_xyyyyz_0_xyyz_0, \
                             g_0_xyyyyz_0_xyzz_0, \
                             g_0_xyyyyz_0_xzzz_0, \
                             g_0_xyyyyz_0_yyyy_0, \
                             g_0_xyyyyz_0_yyyz_0, \
                             g_0_xyyyyz_0_yyzz_0, \
                             g_0_xyyyyz_0_yzzz_0, \
                             g_0_xyyyyz_0_zzzz_0, \
                             g_0_yyyyz_0_xxxz_0,  \
                             g_0_yyyyz_0_xxxz_1,  \
                             g_0_yyyyz_0_xxyz_0,  \
                             g_0_yyyyz_0_xxyz_1,  \
                             g_0_yyyyz_0_xxz_1,   \
                             g_0_yyyyz_0_xxzz_0,  \
                             g_0_yyyyz_0_xxzz_1,  \
                             g_0_yyyyz_0_xyyz_0,  \
                             g_0_yyyyz_0_xyyz_1,  \
                             g_0_yyyyz_0_xyz_1,   \
                             g_0_yyyyz_0_xyzz_0,  \
                             g_0_yyyyz_0_xyzz_1,  \
                             g_0_yyyyz_0_xzz_1,   \
                             g_0_yyyyz_0_xzzz_0,  \
                             g_0_yyyyz_0_xzzz_1,  \
                             g_0_yyyyz_0_yyyy_0,  \
                             g_0_yyyyz_0_yyyy_1,  \
                             g_0_yyyyz_0_yyyz_0,  \
                             g_0_yyyyz_0_yyyz_1,  \
                             g_0_yyyyz_0_yyz_1,   \
                             g_0_yyyyz_0_yyzz_0,  \
                             g_0_yyyyz_0_yyzz_1,  \
                             g_0_yyyyz_0_yzz_1,   \
                             g_0_yyyyz_0_yzzz_0,  \
                             g_0_yyyyz_0_yzzz_1,  \
                             g_0_yyyyz_0_zzz_1,   \
                             g_0_yyyyz_0_zzzz_0,  \
                             g_0_yyyyz_0_zzzz_1,  \
                             wp_x,                \
                             wp_z,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyz_0_xxxx_0[i] = g_0_xyyyy_0_xxxx_0[i] * pb_z + g_0_xyyyy_0_xxxx_1[i] * wp_z[i];

        g_0_xyyyyz_0_xxxy_0[i] = g_0_xyyyy_0_xxxy_0[i] * pb_z + g_0_xyyyy_0_xxxy_1[i] * wp_z[i];

        g_0_xyyyyz_0_xxxz_0[i] = 3.0 * g_0_yyyyz_0_xxz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxxz_0[i] * pb_x + g_0_yyyyz_0_xxxz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxyy_0[i] = g_0_xyyyy_0_xxyy_0[i] * pb_z + g_0_xyyyy_0_xxyy_1[i] * wp_z[i];

        g_0_xyyyyz_0_xxyz_0[i] = 2.0 * g_0_yyyyz_0_xyz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxyz_0[i] * pb_x + g_0_yyyyz_0_xxyz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxzz_0[i] = 2.0 * g_0_yyyyz_0_xzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxzz_0[i] * pb_x + g_0_yyyyz_0_xxzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xyyy_0[i] = g_0_xyyyy_0_xyyy_0[i] * pb_z + g_0_xyyyy_0_xyyy_1[i] * wp_z[i];

        g_0_xyyyyz_0_xyyz_0[i] = g_0_yyyyz_0_yyz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xyyz_0[i] * pb_x + g_0_yyyyz_0_xyyz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xyzz_0[i] = g_0_yyyyz_0_yzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xyzz_0[i] * pb_x + g_0_yyyyz_0_xyzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xzzz_0[i] = g_0_yyyyz_0_zzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xzzz_0[i] * pb_x + g_0_yyyyz_0_xzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_yyyy_0[i] = g_0_yyyyz_0_yyyy_0[i] * pb_x + g_0_yyyyz_0_yyyy_1[i] * wp_x[i];

        g_0_xyyyyz_0_yyyz_0[i] = g_0_yyyyz_0_yyyz_0[i] * pb_x + g_0_yyyyz_0_yyyz_1[i] * wp_x[i];

        g_0_xyyyyz_0_yyzz_0[i] = g_0_yyyyz_0_yyzz_0[i] * pb_x + g_0_yyyyz_0_yyzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_yzzz_0[i] = g_0_yyyyz_0_yzzz_0[i] * pb_x + g_0_yyyyz_0_yzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_zzzz_0[i] = g_0_yyyyz_0_zzzz_0[i] * pb_x + g_0_yyyyz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 255-270 components of targeted buffer : SISG

    auto g_0_xyyyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 255);

    auto g_0_xyyyzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 256);

    auto g_0_xyyyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 257);

    auto g_0_xyyyzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 258);

    auto g_0_xyyyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 259);

    auto g_0_xyyyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 260);

    auto g_0_xyyyzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 261);

    auto g_0_xyyyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 262);

    auto g_0_xyyyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 263);

    auto g_0_xyyyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 264);

    auto g_0_xyyyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 265);

    auto g_0_xyyyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 266);

    auto g_0_xyyyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 267);

    auto g_0_xyyyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 268);

    auto g_0_xyyyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 269);

#pragma omp simd aligned(g_0_xyyyzz_0_xxxx_0,     \
                             g_0_xyyyzz_0_xxxy_0, \
                             g_0_xyyyzz_0_xxxz_0, \
                             g_0_xyyyzz_0_xxyy_0, \
                             g_0_xyyyzz_0_xxyz_0, \
                             g_0_xyyyzz_0_xxzz_0, \
                             g_0_xyyyzz_0_xyyy_0, \
                             g_0_xyyyzz_0_xyyz_0, \
                             g_0_xyyyzz_0_xyzz_0, \
                             g_0_xyyyzz_0_xzzz_0, \
                             g_0_xyyyzz_0_yyyy_0, \
                             g_0_xyyyzz_0_yyyz_0, \
                             g_0_xyyyzz_0_yyzz_0, \
                             g_0_xyyyzz_0_yzzz_0, \
                             g_0_xyyyzz_0_zzzz_0, \
                             g_0_yyyzz_0_xxx_1,   \
                             g_0_yyyzz_0_xxxx_0,  \
                             g_0_yyyzz_0_xxxx_1,  \
                             g_0_yyyzz_0_xxxy_0,  \
                             g_0_yyyzz_0_xxxy_1,  \
                             g_0_yyyzz_0_xxxz_0,  \
                             g_0_yyyzz_0_xxxz_1,  \
                             g_0_yyyzz_0_xxy_1,   \
                             g_0_yyyzz_0_xxyy_0,  \
                             g_0_yyyzz_0_xxyy_1,  \
                             g_0_yyyzz_0_xxyz_0,  \
                             g_0_yyyzz_0_xxyz_1,  \
                             g_0_yyyzz_0_xxz_1,   \
                             g_0_yyyzz_0_xxzz_0,  \
                             g_0_yyyzz_0_xxzz_1,  \
                             g_0_yyyzz_0_xyy_1,   \
                             g_0_yyyzz_0_xyyy_0,  \
                             g_0_yyyzz_0_xyyy_1,  \
                             g_0_yyyzz_0_xyyz_0,  \
                             g_0_yyyzz_0_xyyz_1,  \
                             g_0_yyyzz_0_xyz_1,   \
                             g_0_yyyzz_0_xyzz_0,  \
                             g_0_yyyzz_0_xyzz_1,  \
                             g_0_yyyzz_0_xzz_1,   \
                             g_0_yyyzz_0_xzzz_0,  \
                             g_0_yyyzz_0_xzzz_1,  \
                             g_0_yyyzz_0_yyy_1,   \
                             g_0_yyyzz_0_yyyy_0,  \
                             g_0_yyyzz_0_yyyy_1,  \
                             g_0_yyyzz_0_yyyz_0,  \
                             g_0_yyyzz_0_yyyz_1,  \
                             g_0_yyyzz_0_yyz_1,   \
                             g_0_yyyzz_0_yyzz_0,  \
                             g_0_yyyzz_0_yyzz_1,  \
                             g_0_yyyzz_0_yzz_1,   \
                             g_0_yyyzz_0_yzzz_0,  \
                             g_0_yyyzz_0_yzzz_1,  \
                             g_0_yyyzz_0_zzz_1,   \
                             g_0_yyyzz_0_zzzz_0,  \
                             g_0_yyyzz_0_zzzz_1,  \
                             wp_x,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyzz_0_xxxx_0[i] = 4.0 * g_0_yyyzz_0_xxx_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxx_0[i] * pb_x + g_0_yyyzz_0_xxxx_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxy_0[i] = 3.0 * g_0_yyyzz_0_xxy_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxy_0[i] * pb_x + g_0_yyyzz_0_xxxy_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxz_0[i] = 3.0 * g_0_yyyzz_0_xxz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxz_0[i] * pb_x + g_0_yyyzz_0_xxxz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxyy_0[i] = 2.0 * g_0_yyyzz_0_xyy_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxyy_0[i] * pb_x + g_0_yyyzz_0_xxyy_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxyz_0[i] = 2.0 * g_0_yyyzz_0_xyz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxyz_0[i] * pb_x + g_0_yyyzz_0_xxyz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxzz_0[i] = 2.0 * g_0_yyyzz_0_xzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxzz_0[i] * pb_x + g_0_yyyzz_0_xxzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xyyy_0[i] = g_0_yyyzz_0_yyy_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyyy_0[i] * pb_x + g_0_yyyzz_0_xyyy_1[i] * wp_x[i];

        g_0_xyyyzz_0_xyyz_0[i] = g_0_yyyzz_0_yyz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyyz_0[i] * pb_x + g_0_yyyzz_0_xyyz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xyzz_0[i] = g_0_yyyzz_0_yzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyzz_0[i] * pb_x + g_0_yyyzz_0_xyzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xzzz_0[i] = g_0_yyyzz_0_zzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xzzz_0[i] * pb_x + g_0_yyyzz_0_xzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_yyyy_0[i] = g_0_yyyzz_0_yyyy_0[i] * pb_x + g_0_yyyzz_0_yyyy_1[i] * wp_x[i];

        g_0_xyyyzz_0_yyyz_0[i] = g_0_yyyzz_0_yyyz_0[i] * pb_x + g_0_yyyzz_0_yyyz_1[i] * wp_x[i];

        g_0_xyyyzz_0_yyzz_0[i] = g_0_yyyzz_0_yyzz_0[i] * pb_x + g_0_yyyzz_0_yyzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_yzzz_0[i] = g_0_yyyzz_0_yzzz_0[i] * pb_x + g_0_yyyzz_0_yzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_zzzz_0[i] = g_0_yyyzz_0_zzzz_0[i] * pb_x + g_0_yyyzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 270-285 components of targeted buffer : SISG

    auto g_0_xyyzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 270);

    auto g_0_xyyzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 271);

    auto g_0_xyyzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 272);

    auto g_0_xyyzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 273);

    auto g_0_xyyzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 274);

    auto g_0_xyyzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 275);

    auto g_0_xyyzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 276);

    auto g_0_xyyzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 277);

    auto g_0_xyyzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 278);

    auto g_0_xyyzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 279);

    auto g_0_xyyzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 280);

    auto g_0_xyyzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 281);

    auto g_0_xyyzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 282);

    auto g_0_xyyzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 283);

    auto g_0_xyyzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 284);

#pragma omp simd aligned(g_0_xyyzzz_0_xxxx_0,     \
                             g_0_xyyzzz_0_xxxy_0, \
                             g_0_xyyzzz_0_xxxz_0, \
                             g_0_xyyzzz_0_xxyy_0, \
                             g_0_xyyzzz_0_xxyz_0, \
                             g_0_xyyzzz_0_xxzz_0, \
                             g_0_xyyzzz_0_xyyy_0, \
                             g_0_xyyzzz_0_xyyz_0, \
                             g_0_xyyzzz_0_xyzz_0, \
                             g_0_xyyzzz_0_xzzz_0, \
                             g_0_xyyzzz_0_yyyy_0, \
                             g_0_xyyzzz_0_yyyz_0, \
                             g_0_xyyzzz_0_yyzz_0, \
                             g_0_xyyzzz_0_yzzz_0, \
                             g_0_xyyzzz_0_zzzz_0, \
                             g_0_yyzzz_0_xxx_1,   \
                             g_0_yyzzz_0_xxxx_0,  \
                             g_0_yyzzz_0_xxxx_1,  \
                             g_0_yyzzz_0_xxxy_0,  \
                             g_0_yyzzz_0_xxxy_1,  \
                             g_0_yyzzz_0_xxxz_0,  \
                             g_0_yyzzz_0_xxxz_1,  \
                             g_0_yyzzz_0_xxy_1,   \
                             g_0_yyzzz_0_xxyy_0,  \
                             g_0_yyzzz_0_xxyy_1,  \
                             g_0_yyzzz_0_xxyz_0,  \
                             g_0_yyzzz_0_xxyz_1,  \
                             g_0_yyzzz_0_xxz_1,   \
                             g_0_yyzzz_0_xxzz_0,  \
                             g_0_yyzzz_0_xxzz_1,  \
                             g_0_yyzzz_0_xyy_1,   \
                             g_0_yyzzz_0_xyyy_0,  \
                             g_0_yyzzz_0_xyyy_1,  \
                             g_0_yyzzz_0_xyyz_0,  \
                             g_0_yyzzz_0_xyyz_1,  \
                             g_0_yyzzz_0_xyz_1,   \
                             g_0_yyzzz_0_xyzz_0,  \
                             g_0_yyzzz_0_xyzz_1,  \
                             g_0_yyzzz_0_xzz_1,   \
                             g_0_yyzzz_0_xzzz_0,  \
                             g_0_yyzzz_0_xzzz_1,  \
                             g_0_yyzzz_0_yyy_1,   \
                             g_0_yyzzz_0_yyyy_0,  \
                             g_0_yyzzz_0_yyyy_1,  \
                             g_0_yyzzz_0_yyyz_0,  \
                             g_0_yyzzz_0_yyyz_1,  \
                             g_0_yyzzz_0_yyz_1,   \
                             g_0_yyzzz_0_yyzz_0,  \
                             g_0_yyzzz_0_yyzz_1,  \
                             g_0_yyzzz_0_yzz_1,   \
                             g_0_yyzzz_0_yzzz_0,  \
                             g_0_yyzzz_0_yzzz_1,  \
                             g_0_yyzzz_0_zzz_1,   \
                             g_0_yyzzz_0_zzzz_0,  \
                             g_0_yyzzz_0_zzzz_1,  \
                             wp_x,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyzzz_0_xxxx_0[i] = 4.0 * g_0_yyzzz_0_xxx_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxx_0[i] * pb_x + g_0_yyzzz_0_xxxx_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxy_0[i] = 3.0 * g_0_yyzzz_0_xxy_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxy_0[i] * pb_x + g_0_yyzzz_0_xxxy_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxz_0[i] = 3.0 * g_0_yyzzz_0_xxz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxz_0[i] * pb_x + g_0_yyzzz_0_xxxz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxyy_0[i] = 2.0 * g_0_yyzzz_0_xyy_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxyy_0[i] * pb_x + g_0_yyzzz_0_xxyy_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxyz_0[i] = 2.0 * g_0_yyzzz_0_xyz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxyz_0[i] * pb_x + g_0_yyzzz_0_xxyz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxzz_0[i] = 2.0 * g_0_yyzzz_0_xzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxzz_0[i] * pb_x + g_0_yyzzz_0_xxzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xyyy_0[i] = g_0_yyzzz_0_yyy_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyyy_0[i] * pb_x + g_0_yyzzz_0_xyyy_1[i] * wp_x[i];

        g_0_xyyzzz_0_xyyz_0[i] = g_0_yyzzz_0_yyz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyyz_0[i] * pb_x + g_0_yyzzz_0_xyyz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xyzz_0[i] = g_0_yyzzz_0_yzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyzz_0[i] * pb_x + g_0_yyzzz_0_xyzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xzzz_0[i] = g_0_yyzzz_0_zzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xzzz_0[i] * pb_x + g_0_yyzzz_0_xzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_yyyy_0[i] = g_0_yyzzz_0_yyyy_0[i] * pb_x + g_0_yyzzz_0_yyyy_1[i] * wp_x[i];

        g_0_xyyzzz_0_yyyz_0[i] = g_0_yyzzz_0_yyyz_0[i] * pb_x + g_0_yyzzz_0_yyyz_1[i] * wp_x[i];

        g_0_xyyzzz_0_yyzz_0[i] = g_0_yyzzz_0_yyzz_0[i] * pb_x + g_0_yyzzz_0_yyzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_yzzz_0[i] = g_0_yyzzz_0_yzzz_0[i] * pb_x + g_0_yyzzz_0_yzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_zzzz_0[i] = g_0_yyzzz_0_zzzz_0[i] * pb_x + g_0_yyzzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 285-300 components of targeted buffer : SISG

    auto g_0_xyzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 285);

    auto g_0_xyzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 286);

    auto g_0_xyzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 287);

    auto g_0_xyzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 288);

    auto g_0_xyzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 289);

    auto g_0_xyzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 290);

    auto g_0_xyzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 291);

    auto g_0_xyzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 292);

    auto g_0_xyzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 293);

    auto g_0_xyzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 294);

    auto g_0_xyzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 295);

    auto g_0_xyzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 296);

    auto g_0_xyzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 297);

    auto g_0_xyzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 298);

    auto g_0_xyzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 299);

#pragma omp simd aligned(g_0_xyzzzz_0_xxxx_0,     \
                             g_0_xyzzzz_0_xxxy_0, \
                             g_0_xyzzzz_0_xxxz_0, \
                             g_0_xyzzzz_0_xxyy_0, \
                             g_0_xyzzzz_0_xxyz_0, \
                             g_0_xyzzzz_0_xxzz_0, \
                             g_0_xyzzzz_0_xyyy_0, \
                             g_0_xyzzzz_0_xyyz_0, \
                             g_0_xyzzzz_0_xyzz_0, \
                             g_0_xyzzzz_0_xzzz_0, \
                             g_0_xyzzzz_0_yyyy_0, \
                             g_0_xyzzzz_0_yyyz_0, \
                             g_0_xyzzzz_0_yyzz_0, \
                             g_0_xyzzzz_0_yzzz_0, \
                             g_0_xyzzzz_0_zzzz_0, \
                             g_0_xzzzz_0_xxxx_0,  \
                             g_0_xzzzz_0_xxxx_1,  \
                             g_0_xzzzz_0_xxxz_0,  \
                             g_0_xzzzz_0_xxxz_1,  \
                             g_0_xzzzz_0_xxzz_0,  \
                             g_0_xzzzz_0_xxzz_1,  \
                             g_0_xzzzz_0_xzzz_0,  \
                             g_0_xzzzz_0_xzzz_1,  \
                             g_0_yzzzz_0_xxxy_0,  \
                             g_0_yzzzz_0_xxxy_1,  \
                             g_0_yzzzz_0_xxy_1,   \
                             g_0_yzzzz_0_xxyy_0,  \
                             g_0_yzzzz_0_xxyy_1,  \
                             g_0_yzzzz_0_xxyz_0,  \
                             g_0_yzzzz_0_xxyz_1,  \
                             g_0_yzzzz_0_xyy_1,   \
                             g_0_yzzzz_0_xyyy_0,  \
                             g_0_yzzzz_0_xyyy_1,  \
                             g_0_yzzzz_0_xyyz_0,  \
                             g_0_yzzzz_0_xyyz_1,  \
                             g_0_yzzzz_0_xyz_1,   \
                             g_0_yzzzz_0_xyzz_0,  \
                             g_0_yzzzz_0_xyzz_1,  \
                             g_0_yzzzz_0_yyy_1,   \
                             g_0_yzzzz_0_yyyy_0,  \
                             g_0_yzzzz_0_yyyy_1,  \
                             g_0_yzzzz_0_yyyz_0,  \
                             g_0_yzzzz_0_yyyz_1,  \
                             g_0_yzzzz_0_yyz_1,   \
                             g_0_yzzzz_0_yyzz_0,  \
                             g_0_yzzzz_0_yyzz_1,  \
                             g_0_yzzzz_0_yzz_1,   \
                             g_0_yzzzz_0_yzzz_0,  \
                             g_0_yzzzz_0_yzzz_1,  \
                             g_0_yzzzz_0_zzzz_0,  \
                             g_0_yzzzz_0_zzzz_1,  \
                             wp_x,                \
                             wp_y,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyzzzz_0_xxxx_0[i] = g_0_xzzzz_0_xxxx_0[i] * pb_y + g_0_xzzzz_0_xxxx_1[i] * wp_y[i];

        g_0_xyzzzz_0_xxxy_0[i] = 3.0 * g_0_yzzzz_0_xxy_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxy_0[i] * pb_x + g_0_yzzzz_0_xxxy_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxxz_0[i] = g_0_xzzzz_0_xxxz_0[i] * pb_y + g_0_xzzzz_0_xxxz_1[i] * wp_y[i];

        g_0_xyzzzz_0_xxyy_0[i] = 2.0 * g_0_yzzzz_0_xyy_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxyy_0[i] * pb_x + g_0_yzzzz_0_xxyy_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxyz_0[i] = 2.0 * g_0_yzzzz_0_xyz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxyz_0[i] * pb_x + g_0_yzzzz_0_xxyz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxzz_0[i] = g_0_xzzzz_0_xxzz_0[i] * pb_y + g_0_xzzzz_0_xxzz_1[i] * wp_y[i];

        g_0_xyzzzz_0_xyyy_0[i] = g_0_yzzzz_0_yyy_1[i] * fi_abcd_0 + g_0_yzzzz_0_xyyy_0[i] * pb_x + g_0_yzzzz_0_xyyy_1[i] * wp_x[i];

        g_0_xyzzzz_0_xyyz_0[i] = g_0_yzzzz_0_yyz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xyyz_0[i] * pb_x + g_0_yzzzz_0_xyyz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xyzz_0[i] = g_0_yzzzz_0_yzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xyzz_0[i] * pb_x + g_0_yzzzz_0_xyzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xzzz_0[i] = g_0_xzzzz_0_xzzz_0[i] * pb_y + g_0_xzzzz_0_xzzz_1[i] * wp_y[i];

        g_0_xyzzzz_0_yyyy_0[i] = g_0_yzzzz_0_yyyy_0[i] * pb_x + g_0_yzzzz_0_yyyy_1[i] * wp_x[i];

        g_0_xyzzzz_0_yyyz_0[i] = g_0_yzzzz_0_yyyz_0[i] * pb_x + g_0_yzzzz_0_yyyz_1[i] * wp_x[i];

        g_0_xyzzzz_0_yyzz_0[i] = g_0_yzzzz_0_yyzz_0[i] * pb_x + g_0_yzzzz_0_yyzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_yzzz_0[i] = g_0_yzzzz_0_yzzz_0[i] * pb_x + g_0_yzzzz_0_yzzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_zzzz_0[i] = g_0_yzzzz_0_zzzz_0[i] * pb_x + g_0_yzzzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 300-315 components of targeted buffer : SISG

    auto g_0_xzzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 300);

    auto g_0_xzzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 301);

    auto g_0_xzzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 302);

    auto g_0_xzzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 303);

    auto g_0_xzzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 304);

    auto g_0_xzzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 305);

    auto g_0_xzzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 306);

    auto g_0_xzzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 307);

    auto g_0_xzzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 308);

    auto g_0_xzzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 309);

    auto g_0_xzzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 310);

    auto g_0_xzzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 311);

    auto g_0_xzzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 312);

    auto g_0_xzzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 313);

    auto g_0_xzzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 314);

#pragma omp simd aligned(g_0_xzzzzz_0_xxxx_0,     \
                             g_0_xzzzzz_0_xxxy_0, \
                             g_0_xzzzzz_0_xxxz_0, \
                             g_0_xzzzzz_0_xxyy_0, \
                             g_0_xzzzzz_0_xxyz_0, \
                             g_0_xzzzzz_0_xxzz_0, \
                             g_0_xzzzzz_0_xyyy_0, \
                             g_0_xzzzzz_0_xyyz_0, \
                             g_0_xzzzzz_0_xyzz_0, \
                             g_0_xzzzzz_0_xzzz_0, \
                             g_0_xzzzzz_0_yyyy_0, \
                             g_0_xzzzzz_0_yyyz_0, \
                             g_0_xzzzzz_0_yyzz_0, \
                             g_0_xzzzzz_0_yzzz_0, \
                             g_0_xzzzzz_0_zzzz_0, \
                             g_0_zzzzz_0_xxx_1,   \
                             g_0_zzzzz_0_xxxx_0,  \
                             g_0_zzzzz_0_xxxx_1,  \
                             g_0_zzzzz_0_xxxy_0,  \
                             g_0_zzzzz_0_xxxy_1,  \
                             g_0_zzzzz_0_xxxz_0,  \
                             g_0_zzzzz_0_xxxz_1,  \
                             g_0_zzzzz_0_xxy_1,   \
                             g_0_zzzzz_0_xxyy_0,  \
                             g_0_zzzzz_0_xxyy_1,  \
                             g_0_zzzzz_0_xxyz_0,  \
                             g_0_zzzzz_0_xxyz_1,  \
                             g_0_zzzzz_0_xxz_1,   \
                             g_0_zzzzz_0_xxzz_0,  \
                             g_0_zzzzz_0_xxzz_1,  \
                             g_0_zzzzz_0_xyy_1,   \
                             g_0_zzzzz_0_xyyy_0,  \
                             g_0_zzzzz_0_xyyy_1,  \
                             g_0_zzzzz_0_xyyz_0,  \
                             g_0_zzzzz_0_xyyz_1,  \
                             g_0_zzzzz_0_xyz_1,   \
                             g_0_zzzzz_0_xyzz_0,  \
                             g_0_zzzzz_0_xyzz_1,  \
                             g_0_zzzzz_0_xzz_1,   \
                             g_0_zzzzz_0_xzzz_0,  \
                             g_0_zzzzz_0_xzzz_1,  \
                             g_0_zzzzz_0_yyy_1,   \
                             g_0_zzzzz_0_yyyy_0,  \
                             g_0_zzzzz_0_yyyy_1,  \
                             g_0_zzzzz_0_yyyz_0,  \
                             g_0_zzzzz_0_yyyz_1,  \
                             g_0_zzzzz_0_yyz_1,   \
                             g_0_zzzzz_0_yyzz_0,  \
                             g_0_zzzzz_0_yyzz_1,  \
                             g_0_zzzzz_0_yzz_1,   \
                             g_0_zzzzz_0_yzzz_0,  \
                             g_0_zzzzz_0_yzzz_1,  \
                             g_0_zzzzz_0_zzz_1,   \
                             g_0_zzzzz_0_zzzz_0,  \
                             g_0_zzzzz_0_zzzz_1,  \
                             wp_x,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzzzzz_0_xxxx_0[i] = 4.0 * g_0_zzzzz_0_xxx_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxx_0[i] * pb_x + g_0_zzzzz_0_xxxx_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxy_0[i] = 3.0 * g_0_zzzzz_0_xxy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxy_0[i] * pb_x + g_0_zzzzz_0_xxxy_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxz_0[i] = 3.0 * g_0_zzzzz_0_xxz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxz_0[i] * pb_x + g_0_zzzzz_0_xxxz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxyy_0[i] = 2.0 * g_0_zzzzz_0_xyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyy_0[i] * pb_x + g_0_zzzzz_0_xxyy_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxyz_0[i] = 2.0 * g_0_zzzzz_0_xyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyz_0[i] * pb_x + g_0_zzzzz_0_xxyz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxzz_0[i] = 2.0 * g_0_zzzzz_0_xzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxzz_0[i] * pb_x + g_0_zzzzz_0_xxzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xyyy_0[i] = g_0_zzzzz_0_yyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyy_0[i] * pb_x + g_0_zzzzz_0_xyyy_1[i] * wp_x[i];

        g_0_xzzzzz_0_xyyz_0[i] = g_0_zzzzz_0_yyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyz_0[i] * pb_x + g_0_zzzzz_0_xyyz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xyzz_0[i] = g_0_zzzzz_0_yzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyzz_0[i] * pb_x + g_0_zzzzz_0_xyzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xzzz_0[i] = g_0_zzzzz_0_zzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xzzz_0[i] * pb_x + g_0_zzzzz_0_xzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_yyyy_0[i] = g_0_zzzzz_0_yyyy_0[i] * pb_x + g_0_zzzzz_0_yyyy_1[i] * wp_x[i];

        g_0_xzzzzz_0_yyyz_0[i] = g_0_zzzzz_0_yyyz_0[i] * pb_x + g_0_zzzzz_0_yyyz_1[i] * wp_x[i];

        g_0_xzzzzz_0_yyzz_0[i] = g_0_zzzzz_0_yyzz_0[i] * pb_x + g_0_zzzzz_0_yyzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_yzzz_0[i] = g_0_zzzzz_0_yzzz_0[i] * pb_x + g_0_zzzzz_0_yzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_zzzz_0[i] = g_0_zzzzz_0_zzzz_0[i] * pb_x + g_0_zzzzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 315-330 components of targeted buffer : SISG

    auto g_0_yyyyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 315);

    auto g_0_yyyyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 316);

    auto g_0_yyyyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 317);

    auto g_0_yyyyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 318);

    auto g_0_yyyyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 319);

    auto g_0_yyyyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 320);

    auto g_0_yyyyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 321);

    auto g_0_yyyyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 322);

    auto g_0_yyyyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 323);

    auto g_0_yyyyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 324);

    auto g_0_yyyyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 325);

    auto g_0_yyyyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 326);

    auto g_0_yyyyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 327);

    auto g_0_yyyyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 328);

    auto g_0_yyyyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 329);

#pragma omp simd aligned(g_0_yyyy_0_xxxx_0,       \
                             g_0_yyyy_0_xxxx_1,   \
                             g_0_yyyy_0_xxxy_0,   \
                             g_0_yyyy_0_xxxy_1,   \
                             g_0_yyyy_0_xxxz_0,   \
                             g_0_yyyy_0_xxxz_1,   \
                             g_0_yyyy_0_xxyy_0,   \
                             g_0_yyyy_0_xxyy_1,   \
                             g_0_yyyy_0_xxyz_0,   \
                             g_0_yyyy_0_xxyz_1,   \
                             g_0_yyyy_0_xxzz_0,   \
                             g_0_yyyy_0_xxzz_1,   \
                             g_0_yyyy_0_xyyy_0,   \
                             g_0_yyyy_0_xyyy_1,   \
                             g_0_yyyy_0_xyyz_0,   \
                             g_0_yyyy_0_xyyz_1,   \
                             g_0_yyyy_0_xyzz_0,   \
                             g_0_yyyy_0_xyzz_1,   \
                             g_0_yyyy_0_xzzz_0,   \
                             g_0_yyyy_0_xzzz_1,   \
                             g_0_yyyy_0_yyyy_0,   \
                             g_0_yyyy_0_yyyy_1,   \
                             g_0_yyyy_0_yyyz_0,   \
                             g_0_yyyy_0_yyyz_1,   \
                             g_0_yyyy_0_yyzz_0,   \
                             g_0_yyyy_0_yyzz_1,   \
                             g_0_yyyy_0_yzzz_0,   \
                             g_0_yyyy_0_yzzz_1,   \
                             g_0_yyyy_0_zzzz_0,   \
                             g_0_yyyy_0_zzzz_1,   \
                             g_0_yyyyy_0_xxx_1,   \
                             g_0_yyyyy_0_xxxx_0,  \
                             g_0_yyyyy_0_xxxx_1,  \
                             g_0_yyyyy_0_xxxy_0,  \
                             g_0_yyyyy_0_xxxy_1,  \
                             g_0_yyyyy_0_xxxz_0,  \
                             g_0_yyyyy_0_xxxz_1,  \
                             g_0_yyyyy_0_xxy_1,   \
                             g_0_yyyyy_0_xxyy_0,  \
                             g_0_yyyyy_0_xxyy_1,  \
                             g_0_yyyyy_0_xxyz_0,  \
                             g_0_yyyyy_0_xxyz_1,  \
                             g_0_yyyyy_0_xxz_1,   \
                             g_0_yyyyy_0_xxzz_0,  \
                             g_0_yyyyy_0_xxzz_1,  \
                             g_0_yyyyy_0_xyy_1,   \
                             g_0_yyyyy_0_xyyy_0,  \
                             g_0_yyyyy_0_xyyy_1,  \
                             g_0_yyyyy_0_xyyz_0,  \
                             g_0_yyyyy_0_xyyz_1,  \
                             g_0_yyyyy_0_xyz_1,   \
                             g_0_yyyyy_0_xyzz_0,  \
                             g_0_yyyyy_0_xyzz_1,  \
                             g_0_yyyyy_0_xzz_1,   \
                             g_0_yyyyy_0_xzzz_0,  \
                             g_0_yyyyy_0_xzzz_1,  \
                             g_0_yyyyy_0_yyy_1,   \
                             g_0_yyyyy_0_yyyy_0,  \
                             g_0_yyyyy_0_yyyy_1,  \
                             g_0_yyyyy_0_yyyz_0,  \
                             g_0_yyyyy_0_yyyz_1,  \
                             g_0_yyyyy_0_yyz_1,   \
                             g_0_yyyyy_0_yyzz_0,  \
                             g_0_yyyyy_0_yyzz_1,  \
                             g_0_yyyyy_0_yzz_1,   \
                             g_0_yyyyy_0_yzzz_0,  \
                             g_0_yyyyy_0_yzzz_1,  \
                             g_0_yyyyy_0_zzz_1,   \
                             g_0_yyyyy_0_zzzz_0,  \
                             g_0_yyyyy_0_zzzz_1,  \
                             g_0_yyyyyy_0_xxxx_0, \
                             g_0_yyyyyy_0_xxxy_0, \
                             g_0_yyyyyy_0_xxxz_0, \
                             g_0_yyyyyy_0_xxyy_0, \
                             g_0_yyyyyy_0_xxyz_0, \
                             g_0_yyyyyy_0_xxzz_0, \
                             g_0_yyyyyy_0_xyyy_0, \
                             g_0_yyyyyy_0_xyyz_0, \
                             g_0_yyyyyy_0_xyzz_0, \
                             g_0_yyyyyy_0_xzzz_0, \
                             g_0_yyyyyy_0_yyyy_0, \
                             g_0_yyyyyy_0_yyyz_0, \
                             g_0_yyyyyy_0_yyzz_0, \
                             g_0_yyyyyy_0_yzzz_0, \
                             g_0_yyyyyy_0_zzzz_0, \
                             wp_y,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyyy_0_xxxx_0[i] = 5.0 * g_0_yyyy_0_xxxx_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxx_1[i] * fti_ab_0 + g_0_yyyyy_0_xxxx_0[i] * pb_y +
                                 g_0_yyyyy_0_xxxx_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxy_0[i] = 5.0 * g_0_yyyy_0_xxxy_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxy_1[i] * fti_ab_0 + g_0_yyyyy_0_xxx_1[i] * fi_abcd_0 +
                                 g_0_yyyyy_0_xxxy_0[i] * pb_y + g_0_yyyyy_0_xxxy_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxz_0[i] = 5.0 * g_0_yyyy_0_xxxz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxz_1[i] * fti_ab_0 + g_0_yyyyy_0_xxxz_0[i] * pb_y +
                                 g_0_yyyyy_0_xxxz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxyy_0[i] = 5.0 * g_0_yyyy_0_xxyy_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxyy_1[i] * fti_ab_0 +
                                 2.0 * g_0_yyyyy_0_xxy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyy_0[i] * pb_y + g_0_yyyyy_0_xxyy_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxyz_0[i] = 5.0 * g_0_yyyy_0_xxyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxyz_1[i] * fti_ab_0 + g_0_yyyyy_0_xxz_1[i] * fi_abcd_0 +
                                 g_0_yyyyy_0_xxyz_0[i] * pb_y + g_0_yyyyy_0_xxyz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxzz_0[i] = 5.0 * g_0_yyyy_0_xxzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxzz_1[i] * fti_ab_0 + g_0_yyyyy_0_xxzz_0[i] * pb_y +
                                 g_0_yyyyy_0_xxzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xyyy_0[i] = 5.0 * g_0_yyyy_0_xyyy_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xyyy_1[i] * fti_ab_0 +
                                 3.0 * g_0_yyyyy_0_xyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyy_0[i] * pb_y + g_0_yyyyy_0_xyyy_1[i] * wp_y[i];

        g_0_yyyyyy_0_xyyz_0[i] = 5.0 * g_0_yyyy_0_xyyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xyyz_1[i] * fti_ab_0 +
                                 2.0 * g_0_yyyyy_0_xyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyz_0[i] * pb_y + g_0_yyyyy_0_xyyz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xyzz_0[i] = 5.0 * g_0_yyyy_0_xyzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xyzz_1[i] * fti_ab_0 + g_0_yyyyy_0_xzz_1[i] * fi_abcd_0 +
                                 g_0_yyyyy_0_xyzz_0[i] * pb_y + g_0_yyyyy_0_xyzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xzzz_0[i] = 5.0 * g_0_yyyy_0_xzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xzzz_1[i] * fti_ab_0 + g_0_yyyyy_0_xzzz_0[i] * pb_y +
                                 g_0_yyyyy_0_xzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_yyyy_0[i] = 5.0 * g_0_yyyy_0_yyyy_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_yyyy_1[i] * fti_ab_0 +
                                 4.0 * g_0_yyyyy_0_yyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyyy_0[i] * pb_y + g_0_yyyyy_0_yyyy_1[i] * wp_y[i];

        g_0_yyyyyy_0_yyyz_0[i] = 5.0 * g_0_yyyy_0_yyyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_yyyz_1[i] * fti_ab_0 +
                                 3.0 * g_0_yyyyy_0_yyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyyz_0[i] * pb_y + g_0_yyyyy_0_yyyz_1[i] * wp_y[i];

        g_0_yyyyyy_0_yyzz_0[i] = 5.0 * g_0_yyyy_0_yyzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_yyzz_1[i] * fti_ab_0 +
                                 2.0 * g_0_yyyyy_0_yzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyzz_0[i] * pb_y + g_0_yyyyy_0_yyzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_yzzz_0[i] = 5.0 * g_0_yyyy_0_yzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_yzzz_1[i] * fti_ab_0 + g_0_yyyyy_0_zzz_1[i] * fi_abcd_0 +
                                 g_0_yyyyy_0_yzzz_0[i] * pb_y + g_0_yyyyy_0_yzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_zzzz_0[i] = 5.0 * g_0_yyyy_0_zzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_zzzz_1[i] * fti_ab_0 + g_0_yyyyy_0_zzzz_0[i] * pb_y +
                                 g_0_yyyyy_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 330-345 components of targeted buffer : SISG

    auto g_0_yyyyyz_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 330);

    auto g_0_yyyyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 331);

    auto g_0_yyyyyz_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 332);

    auto g_0_yyyyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 333);

    auto g_0_yyyyyz_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 334);

    auto g_0_yyyyyz_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 335);

    auto g_0_yyyyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 336);

    auto g_0_yyyyyz_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 337);

    auto g_0_yyyyyz_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 338);

    auto g_0_yyyyyz_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 339);

    auto g_0_yyyyyz_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 340);

    auto g_0_yyyyyz_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 341);

    auto g_0_yyyyyz_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 342);

    auto g_0_yyyyyz_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 343);

    auto g_0_yyyyyz_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 344);

#pragma omp simd aligned(g_0_yyyyy_0_xxx_1,       \
                             g_0_yyyyy_0_xxxx_0,  \
                             g_0_yyyyy_0_xxxx_1,  \
                             g_0_yyyyy_0_xxxy_0,  \
                             g_0_yyyyy_0_xxxy_1,  \
                             g_0_yyyyy_0_xxxz_0,  \
                             g_0_yyyyy_0_xxxz_1,  \
                             g_0_yyyyy_0_xxy_1,   \
                             g_0_yyyyy_0_xxyy_0,  \
                             g_0_yyyyy_0_xxyy_1,  \
                             g_0_yyyyy_0_xxyz_0,  \
                             g_0_yyyyy_0_xxyz_1,  \
                             g_0_yyyyy_0_xxz_1,   \
                             g_0_yyyyy_0_xxzz_0,  \
                             g_0_yyyyy_0_xxzz_1,  \
                             g_0_yyyyy_0_xyy_1,   \
                             g_0_yyyyy_0_xyyy_0,  \
                             g_0_yyyyy_0_xyyy_1,  \
                             g_0_yyyyy_0_xyyz_0,  \
                             g_0_yyyyy_0_xyyz_1,  \
                             g_0_yyyyy_0_xyz_1,   \
                             g_0_yyyyy_0_xyzz_0,  \
                             g_0_yyyyy_0_xyzz_1,  \
                             g_0_yyyyy_0_xzz_1,   \
                             g_0_yyyyy_0_xzzz_0,  \
                             g_0_yyyyy_0_xzzz_1,  \
                             g_0_yyyyy_0_yyy_1,   \
                             g_0_yyyyy_0_yyyy_0,  \
                             g_0_yyyyy_0_yyyy_1,  \
                             g_0_yyyyy_0_yyyz_0,  \
                             g_0_yyyyy_0_yyyz_1,  \
                             g_0_yyyyy_0_yyz_1,   \
                             g_0_yyyyy_0_yyzz_0,  \
                             g_0_yyyyy_0_yyzz_1,  \
                             g_0_yyyyy_0_yzz_1,   \
                             g_0_yyyyy_0_yzzz_0,  \
                             g_0_yyyyy_0_yzzz_1,  \
                             g_0_yyyyy_0_zzz_1,   \
                             g_0_yyyyy_0_zzzz_0,  \
                             g_0_yyyyy_0_zzzz_1,  \
                             g_0_yyyyyz_0_xxxx_0, \
                             g_0_yyyyyz_0_xxxy_0, \
                             g_0_yyyyyz_0_xxxz_0, \
                             g_0_yyyyyz_0_xxyy_0, \
                             g_0_yyyyyz_0_xxyz_0, \
                             g_0_yyyyyz_0_xxzz_0, \
                             g_0_yyyyyz_0_xyyy_0, \
                             g_0_yyyyyz_0_xyyz_0, \
                             g_0_yyyyyz_0_xyzz_0, \
                             g_0_yyyyyz_0_xzzz_0, \
                             g_0_yyyyyz_0_yyyy_0, \
                             g_0_yyyyyz_0_yyyz_0, \
                             g_0_yyyyyz_0_yyzz_0, \
                             g_0_yyyyyz_0_yzzz_0, \
                             g_0_yyyyyz_0_zzzz_0, \
                             wp_z,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyyyz_0_xxxx_0[i] = g_0_yyyyy_0_xxxx_0[i] * pb_z + g_0_yyyyy_0_xxxx_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxy_0[i] = g_0_yyyyy_0_xxxy_0[i] * pb_z + g_0_yyyyy_0_xxxy_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxz_0[i] = g_0_yyyyy_0_xxx_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxz_0[i] * pb_z + g_0_yyyyy_0_xxxz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxyy_0[i] = g_0_yyyyy_0_xxyy_0[i] * pb_z + g_0_yyyyy_0_xxyy_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxyz_0[i] = g_0_yyyyy_0_xxy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyz_0[i] * pb_z + g_0_yyyyy_0_xxyz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxzz_0[i] = 2.0 * g_0_yyyyy_0_xxz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxzz_0[i] * pb_z + g_0_yyyyy_0_xxzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xyyy_0[i] = g_0_yyyyy_0_xyyy_0[i] * pb_z + g_0_yyyyy_0_xyyy_1[i] * wp_z[i];

        g_0_yyyyyz_0_xyyz_0[i] = g_0_yyyyy_0_xyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyz_0[i] * pb_z + g_0_yyyyy_0_xyyz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xyzz_0[i] = 2.0 * g_0_yyyyy_0_xyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyzz_0[i] * pb_z + g_0_yyyyy_0_xyzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xzzz_0[i] = 3.0 * g_0_yyyyy_0_xzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xzzz_0[i] * pb_z + g_0_yyyyy_0_xzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_yyyy_0[i] = g_0_yyyyy_0_yyyy_0[i] * pb_z + g_0_yyyyy_0_yyyy_1[i] * wp_z[i];

        g_0_yyyyyz_0_yyyz_0[i] = g_0_yyyyy_0_yyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyyz_0[i] * pb_z + g_0_yyyyy_0_yyyz_1[i] * wp_z[i];

        g_0_yyyyyz_0_yyzz_0[i] = 2.0 * g_0_yyyyy_0_yyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyzz_0[i] * pb_z + g_0_yyyyy_0_yyzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_yzzz_0[i] = 3.0 * g_0_yyyyy_0_yzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yzzz_0[i] * pb_z + g_0_yyyyy_0_yzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_zzzz_0[i] = 4.0 * g_0_yyyyy_0_zzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_zzzz_0[i] * pb_z + g_0_yyyyy_0_zzzz_1[i] * wp_z[i];
    }

    /// Set up 345-360 components of targeted buffer : SISG

    auto g_0_yyyyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 345);

    auto g_0_yyyyzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 346);

    auto g_0_yyyyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 347);

    auto g_0_yyyyzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 348);

    auto g_0_yyyyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 349);

    auto g_0_yyyyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 350);

    auto g_0_yyyyzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 351);

    auto g_0_yyyyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 352);

    auto g_0_yyyyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 353);

    auto g_0_yyyyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 354);

    auto g_0_yyyyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 355);

    auto g_0_yyyyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 356);

    auto g_0_yyyyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 357);

    auto g_0_yyyyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 358);

    auto g_0_yyyyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 359);

#pragma omp simd aligned(g_0_yyyy_0_xxxy_0,       \
                             g_0_yyyy_0_xxxy_1,   \
                             g_0_yyyy_0_xxyy_0,   \
                             g_0_yyyy_0_xxyy_1,   \
                             g_0_yyyy_0_xyyy_0,   \
                             g_0_yyyy_0_xyyy_1,   \
                             g_0_yyyy_0_yyyy_0,   \
                             g_0_yyyy_0_yyyy_1,   \
                             g_0_yyyyz_0_xxxy_0,  \
                             g_0_yyyyz_0_xxxy_1,  \
                             g_0_yyyyz_0_xxyy_0,  \
                             g_0_yyyyz_0_xxyy_1,  \
                             g_0_yyyyz_0_xyyy_0,  \
                             g_0_yyyyz_0_xyyy_1,  \
                             g_0_yyyyz_0_yyyy_0,  \
                             g_0_yyyyz_0_yyyy_1,  \
                             g_0_yyyyzz_0_xxxx_0, \
                             g_0_yyyyzz_0_xxxy_0, \
                             g_0_yyyyzz_0_xxxz_0, \
                             g_0_yyyyzz_0_xxyy_0, \
                             g_0_yyyyzz_0_xxyz_0, \
                             g_0_yyyyzz_0_xxzz_0, \
                             g_0_yyyyzz_0_xyyy_0, \
                             g_0_yyyyzz_0_xyyz_0, \
                             g_0_yyyyzz_0_xyzz_0, \
                             g_0_yyyyzz_0_xzzz_0, \
                             g_0_yyyyzz_0_yyyy_0, \
                             g_0_yyyyzz_0_yyyz_0, \
                             g_0_yyyyzz_0_yyzz_0, \
                             g_0_yyyyzz_0_yzzz_0, \
                             g_0_yyyyzz_0_zzzz_0, \
                             g_0_yyyzz_0_xxxx_0,  \
                             g_0_yyyzz_0_xxxx_1,  \
                             g_0_yyyzz_0_xxxz_0,  \
                             g_0_yyyzz_0_xxxz_1,  \
                             g_0_yyyzz_0_xxyz_0,  \
                             g_0_yyyzz_0_xxyz_1,  \
                             g_0_yyyzz_0_xxz_1,   \
                             g_0_yyyzz_0_xxzz_0,  \
                             g_0_yyyzz_0_xxzz_1,  \
                             g_0_yyyzz_0_xyyz_0,  \
                             g_0_yyyzz_0_xyyz_1,  \
                             g_0_yyyzz_0_xyz_1,   \
                             g_0_yyyzz_0_xyzz_0,  \
                             g_0_yyyzz_0_xyzz_1,  \
                             g_0_yyyzz_0_xzz_1,   \
                             g_0_yyyzz_0_xzzz_0,  \
                             g_0_yyyzz_0_xzzz_1,  \
                             g_0_yyyzz_0_yyyz_0,  \
                             g_0_yyyzz_0_yyyz_1,  \
                             g_0_yyyzz_0_yyz_1,   \
                             g_0_yyyzz_0_yyzz_0,  \
                             g_0_yyyzz_0_yyzz_1,  \
                             g_0_yyyzz_0_yzz_1,   \
                             g_0_yyyzz_0_yzzz_0,  \
                             g_0_yyyzz_0_yzzz_1,  \
                             g_0_yyyzz_0_zzz_1,   \
                             g_0_yyyzz_0_zzzz_0,  \
                             g_0_yyyzz_0_zzzz_1,  \
                             g_0_yyzz_0_xxxx_0,   \
                             g_0_yyzz_0_xxxx_1,   \
                             g_0_yyzz_0_xxxz_0,   \
                             g_0_yyzz_0_xxxz_1,   \
                             g_0_yyzz_0_xxyz_0,   \
                             g_0_yyzz_0_xxyz_1,   \
                             g_0_yyzz_0_xxzz_0,   \
                             g_0_yyzz_0_xxzz_1,   \
                             g_0_yyzz_0_xyyz_0,   \
                             g_0_yyzz_0_xyyz_1,   \
                             g_0_yyzz_0_xyzz_0,   \
                             g_0_yyzz_0_xyzz_1,   \
                             g_0_yyzz_0_xzzz_0,   \
                             g_0_yyzz_0_xzzz_1,   \
                             g_0_yyzz_0_yyyz_0,   \
                             g_0_yyzz_0_yyyz_1,   \
                             g_0_yyzz_0_yyzz_0,   \
                             g_0_yyzz_0_yyzz_1,   \
                             g_0_yyzz_0_yzzz_0,   \
                             g_0_yyzz_0_yzzz_1,   \
                             g_0_yyzz_0_zzzz_0,   \
                             g_0_yyzz_0_zzzz_1,   \
                             wp_y,                \
                             wp_z,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyzz_0_xxxx_0[i] = 3.0 * g_0_yyzz_0_xxxx_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxx_1[i] * fti_ab_0 + g_0_yyyzz_0_xxxx_0[i] * pb_y +
                                 g_0_yyyzz_0_xxxx_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxxy_0[i] =
            g_0_yyyy_0_xxxy_0[i] * fi_ab_0 - g_0_yyyy_0_xxxy_1[i] * fti_ab_0 + g_0_yyyyz_0_xxxy_0[i] * pb_z + g_0_yyyyz_0_xxxy_1[i] * wp_z[i];

        g_0_yyyyzz_0_xxxz_0[i] = 3.0 * g_0_yyzz_0_xxxz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxz_1[i] * fti_ab_0 + g_0_yyyzz_0_xxxz_0[i] * pb_y +
                                 g_0_yyyzz_0_xxxz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxyy_0[i] =
            g_0_yyyy_0_xxyy_0[i] * fi_ab_0 - g_0_yyyy_0_xxyy_1[i] * fti_ab_0 + g_0_yyyyz_0_xxyy_0[i] * pb_z + g_0_yyyyz_0_xxyy_1[i] * wp_z[i];

        g_0_yyyyzz_0_xxyz_0[i] = 3.0 * g_0_yyzz_0_xxyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxyz_1[i] * fti_ab_0 + g_0_yyyzz_0_xxz_1[i] * fi_abcd_0 +
                                 g_0_yyyzz_0_xxyz_0[i] * pb_y + g_0_yyyzz_0_xxyz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxzz_0[i] = 3.0 * g_0_yyzz_0_xxzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxzz_1[i] * fti_ab_0 + g_0_yyyzz_0_xxzz_0[i] * pb_y +
                                 g_0_yyyzz_0_xxzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xyyy_0[i] =
            g_0_yyyy_0_xyyy_0[i] * fi_ab_0 - g_0_yyyy_0_xyyy_1[i] * fti_ab_0 + g_0_yyyyz_0_xyyy_0[i] * pb_z + g_0_yyyyz_0_xyyy_1[i] * wp_z[i];

        g_0_yyyyzz_0_xyyz_0[i] = 3.0 * g_0_yyzz_0_xyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xyyz_1[i] * fti_ab_0 +
                                 2.0 * g_0_yyyzz_0_xyz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyyz_0[i] * pb_y + g_0_yyyzz_0_xyyz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xyzz_0[i] = 3.0 * g_0_yyzz_0_xyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xyzz_1[i] * fti_ab_0 + g_0_yyyzz_0_xzz_1[i] * fi_abcd_0 +
                                 g_0_yyyzz_0_xyzz_0[i] * pb_y + g_0_yyyzz_0_xyzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xzzz_0[i] = 3.0 * g_0_yyzz_0_xzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xzzz_1[i] * fti_ab_0 + g_0_yyyzz_0_xzzz_0[i] * pb_y +
                                 g_0_yyyzz_0_xzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_yyyy_0[i] =
            g_0_yyyy_0_yyyy_0[i] * fi_ab_0 - g_0_yyyy_0_yyyy_1[i] * fti_ab_0 + g_0_yyyyz_0_yyyy_0[i] * pb_z + g_0_yyyyz_0_yyyy_1[i] * wp_z[i];

        g_0_yyyyzz_0_yyyz_0[i] = 3.0 * g_0_yyzz_0_yyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_yyyz_1[i] * fti_ab_0 +
                                 3.0 * g_0_yyyzz_0_yyz_1[i] * fi_abcd_0 + g_0_yyyzz_0_yyyz_0[i] * pb_y + g_0_yyyzz_0_yyyz_1[i] * wp_y[i];

        g_0_yyyyzz_0_yyzz_0[i] = 3.0 * g_0_yyzz_0_yyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_yyzz_1[i] * fti_ab_0 +
                                 2.0 * g_0_yyyzz_0_yzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_yyzz_0[i] * pb_y + g_0_yyyzz_0_yyzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_yzzz_0[i] = 3.0 * g_0_yyzz_0_yzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_yzzz_1[i] * fti_ab_0 + g_0_yyyzz_0_zzz_1[i] * fi_abcd_0 +
                                 g_0_yyyzz_0_yzzz_0[i] * pb_y + g_0_yyyzz_0_yzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_zzzz_0[i] = 3.0 * g_0_yyzz_0_zzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_zzzz_1[i] * fti_ab_0 + g_0_yyyzz_0_zzzz_0[i] * pb_y +
                                 g_0_yyyzz_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 360-375 components of targeted buffer : SISG

    auto g_0_yyyzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 360);

    auto g_0_yyyzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 361);

    auto g_0_yyyzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 362);

    auto g_0_yyyzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 363);

    auto g_0_yyyzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 364);

    auto g_0_yyyzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 365);

    auto g_0_yyyzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 366);

    auto g_0_yyyzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 367);

    auto g_0_yyyzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 368);

    auto g_0_yyyzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 369);

    auto g_0_yyyzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 370);

    auto g_0_yyyzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 371);

    auto g_0_yyyzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 372);

    auto g_0_yyyzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 373);

    auto g_0_yyyzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 374);

#pragma omp simd aligned(g_0_yyyz_0_xxxy_0,       \
                             g_0_yyyz_0_xxxy_1,   \
                             g_0_yyyz_0_xxyy_0,   \
                             g_0_yyyz_0_xxyy_1,   \
                             g_0_yyyz_0_xyyy_0,   \
                             g_0_yyyz_0_xyyy_1,   \
                             g_0_yyyz_0_yyyy_0,   \
                             g_0_yyyz_0_yyyy_1,   \
                             g_0_yyyzz_0_xxxy_0,  \
                             g_0_yyyzz_0_xxxy_1,  \
                             g_0_yyyzz_0_xxyy_0,  \
                             g_0_yyyzz_0_xxyy_1,  \
                             g_0_yyyzz_0_xyyy_0,  \
                             g_0_yyyzz_0_xyyy_1,  \
                             g_0_yyyzz_0_yyyy_0,  \
                             g_0_yyyzz_0_yyyy_1,  \
                             g_0_yyyzzz_0_xxxx_0, \
                             g_0_yyyzzz_0_xxxy_0, \
                             g_0_yyyzzz_0_xxxz_0, \
                             g_0_yyyzzz_0_xxyy_0, \
                             g_0_yyyzzz_0_xxyz_0, \
                             g_0_yyyzzz_0_xxzz_0, \
                             g_0_yyyzzz_0_xyyy_0, \
                             g_0_yyyzzz_0_xyyz_0, \
                             g_0_yyyzzz_0_xyzz_0, \
                             g_0_yyyzzz_0_xzzz_0, \
                             g_0_yyyzzz_0_yyyy_0, \
                             g_0_yyyzzz_0_yyyz_0, \
                             g_0_yyyzzz_0_yyzz_0, \
                             g_0_yyyzzz_0_yzzz_0, \
                             g_0_yyyzzz_0_zzzz_0, \
                             g_0_yyzzz_0_xxxx_0,  \
                             g_0_yyzzz_0_xxxx_1,  \
                             g_0_yyzzz_0_xxxz_0,  \
                             g_0_yyzzz_0_xxxz_1,  \
                             g_0_yyzzz_0_xxyz_0,  \
                             g_0_yyzzz_0_xxyz_1,  \
                             g_0_yyzzz_0_xxz_1,   \
                             g_0_yyzzz_0_xxzz_0,  \
                             g_0_yyzzz_0_xxzz_1,  \
                             g_0_yyzzz_0_xyyz_0,  \
                             g_0_yyzzz_0_xyyz_1,  \
                             g_0_yyzzz_0_xyz_1,   \
                             g_0_yyzzz_0_xyzz_0,  \
                             g_0_yyzzz_0_xyzz_1,  \
                             g_0_yyzzz_0_xzz_1,   \
                             g_0_yyzzz_0_xzzz_0,  \
                             g_0_yyzzz_0_xzzz_1,  \
                             g_0_yyzzz_0_yyyz_0,  \
                             g_0_yyzzz_0_yyyz_1,  \
                             g_0_yyzzz_0_yyz_1,   \
                             g_0_yyzzz_0_yyzz_0,  \
                             g_0_yyzzz_0_yyzz_1,  \
                             g_0_yyzzz_0_yzz_1,   \
                             g_0_yyzzz_0_yzzz_0,  \
                             g_0_yyzzz_0_yzzz_1,  \
                             g_0_yyzzz_0_zzz_1,   \
                             g_0_yyzzz_0_zzzz_0,  \
                             g_0_yyzzz_0_zzzz_1,  \
                             g_0_yzzz_0_xxxx_0,   \
                             g_0_yzzz_0_xxxx_1,   \
                             g_0_yzzz_0_xxxz_0,   \
                             g_0_yzzz_0_xxxz_1,   \
                             g_0_yzzz_0_xxyz_0,   \
                             g_0_yzzz_0_xxyz_1,   \
                             g_0_yzzz_0_xxzz_0,   \
                             g_0_yzzz_0_xxzz_1,   \
                             g_0_yzzz_0_xyyz_0,   \
                             g_0_yzzz_0_xyyz_1,   \
                             g_0_yzzz_0_xyzz_0,   \
                             g_0_yzzz_0_xyzz_1,   \
                             g_0_yzzz_0_xzzz_0,   \
                             g_0_yzzz_0_xzzz_1,   \
                             g_0_yzzz_0_yyyz_0,   \
                             g_0_yzzz_0_yyyz_1,   \
                             g_0_yzzz_0_yyzz_0,   \
                             g_0_yzzz_0_yyzz_1,   \
                             g_0_yzzz_0_yzzz_0,   \
                             g_0_yzzz_0_yzzz_1,   \
                             g_0_yzzz_0_zzzz_0,   \
                             g_0_yzzz_0_zzzz_1,   \
                             wp_y,                \
                             wp_z,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyzzz_0_xxxx_0[i] = 2.0 * g_0_yzzz_0_xxxx_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxxx_1[i] * fti_ab_0 + g_0_yyzzz_0_xxxx_0[i] * pb_y +
                                 g_0_yyzzz_0_xxxx_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxxy_0[i] = 2.0 * g_0_yyyz_0_xxxy_0[i] * fi_ab_0 - 2.0 * g_0_yyyz_0_xxxy_1[i] * fti_ab_0 + g_0_yyyzz_0_xxxy_0[i] * pb_z +
                                 g_0_yyyzz_0_xxxy_1[i] * wp_z[i];

        g_0_yyyzzz_0_xxxz_0[i] = 2.0 * g_0_yzzz_0_xxxz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxxz_1[i] * fti_ab_0 + g_0_yyzzz_0_xxxz_0[i] * pb_y +
                                 g_0_yyzzz_0_xxxz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxyy_0[i] = 2.0 * g_0_yyyz_0_xxyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyz_0_xxyy_1[i] * fti_ab_0 + g_0_yyyzz_0_xxyy_0[i] * pb_z +
                                 g_0_yyyzz_0_xxyy_1[i] * wp_z[i];

        g_0_yyyzzz_0_xxyz_0[i] = 2.0 * g_0_yzzz_0_xxyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxyz_1[i] * fti_ab_0 + g_0_yyzzz_0_xxz_1[i] * fi_abcd_0 +
                                 g_0_yyzzz_0_xxyz_0[i] * pb_y + g_0_yyzzz_0_xxyz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxzz_0[i] = 2.0 * g_0_yzzz_0_xxzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxzz_1[i] * fti_ab_0 + g_0_yyzzz_0_xxzz_0[i] * pb_y +
                                 g_0_yyzzz_0_xxzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xyyy_0[i] = 2.0 * g_0_yyyz_0_xyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyz_0_xyyy_1[i] * fti_ab_0 + g_0_yyyzz_0_xyyy_0[i] * pb_z +
                                 g_0_yyyzz_0_xyyy_1[i] * wp_z[i];

        g_0_yyyzzz_0_xyyz_0[i] = 2.0 * g_0_yzzz_0_xyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xyyz_1[i] * fti_ab_0 +
                                 2.0 * g_0_yyzzz_0_xyz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyyz_0[i] * pb_y + g_0_yyzzz_0_xyyz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xyzz_0[i] = 2.0 * g_0_yzzz_0_xyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xyzz_1[i] * fti_ab_0 + g_0_yyzzz_0_xzz_1[i] * fi_abcd_0 +
                                 g_0_yyzzz_0_xyzz_0[i] * pb_y + g_0_yyzzz_0_xyzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xzzz_0[i] = 2.0 * g_0_yzzz_0_xzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xzzz_1[i] * fti_ab_0 + g_0_yyzzz_0_xzzz_0[i] * pb_y +
                                 g_0_yyzzz_0_xzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_yyyy_0[i] = 2.0 * g_0_yyyz_0_yyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyz_0_yyyy_1[i] * fti_ab_0 + g_0_yyyzz_0_yyyy_0[i] * pb_z +
                                 g_0_yyyzz_0_yyyy_1[i] * wp_z[i];

        g_0_yyyzzz_0_yyyz_0[i] = 2.0 * g_0_yzzz_0_yyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_yyyz_1[i] * fti_ab_0 +
                                 3.0 * g_0_yyzzz_0_yyz_1[i] * fi_abcd_0 + g_0_yyzzz_0_yyyz_0[i] * pb_y + g_0_yyzzz_0_yyyz_1[i] * wp_y[i];

        g_0_yyyzzz_0_yyzz_0[i] = 2.0 * g_0_yzzz_0_yyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_yyzz_1[i] * fti_ab_0 +
                                 2.0 * g_0_yyzzz_0_yzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_yyzz_0[i] * pb_y + g_0_yyzzz_0_yyzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_yzzz_0[i] = 2.0 * g_0_yzzz_0_yzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_yzzz_1[i] * fti_ab_0 + g_0_yyzzz_0_zzz_1[i] * fi_abcd_0 +
                                 g_0_yyzzz_0_yzzz_0[i] * pb_y + g_0_yyzzz_0_yzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_zzzz_0[i] = 2.0 * g_0_yzzz_0_zzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_zzzz_1[i] * fti_ab_0 + g_0_yyzzz_0_zzzz_0[i] * pb_y +
                                 g_0_yyzzz_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 375-390 components of targeted buffer : SISG

    auto g_0_yyzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 375);

    auto g_0_yyzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 376);

    auto g_0_yyzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 377);

    auto g_0_yyzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 378);

    auto g_0_yyzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 379);

    auto g_0_yyzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 380);

    auto g_0_yyzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 381);

    auto g_0_yyzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 382);

    auto g_0_yyzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 383);

    auto g_0_yyzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 384);

    auto g_0_yyzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 385);

    auto g_0_yyzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 386);

    auto g_0_yyzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 387);

    auto g_0_yyzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 388);

    auto g_0_yyzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 389);

#pragma omp simd aligned(g_0_yyzz_0_xxxy_0,       \
                             g_0_yyzz_0_xxxy_1,   \
                             g_0_yyzz_0_xxyy_0,   \
                             g_0_yyzz_0_xxyy_1,   \
                             g_0_yyzz_0_xyyy_0,   \
                             g_0_yyzz_0_xyyy_1,   \
                             g_0_yyzz_0_yyyy_0,   \
                             g_0_yyzz_0_yyyy_1,   \
                             g_0_yyzzz_0_xxxy_0,  \
                             g_0_yyzzz_0_xxxy_1,  \
                             g_0_yyzzz_0_xxyy_0,  \
                             g_0_yyzzz_0_xxyy_1,  \
                             g_0_yyzzz_0_xyyy_0,  \
                             g_0_yyzzz_0_xyyy_1,  \
                             g_0_yyzzz_0_yyyy_0,  \
                             g_0_yyzzz_0_yyyy_1,  \
                             g_0_yyzzzz_0_xxxx_0, \
                             g_0_yyzzzz_0_xxxy_0, \
                             g_0_yyzzzz_0_xxxz_0, \
                             g_0_yyzzzz_0_xxyy_0, \
                             g_0_yyzzzz_0_xxyz_0, \
                             g_0_yyzzzz_0_xxzz_0, \
                             g_0_yyzzzz_0_xyyy_0, \
                             g_0_yyzzzz_0_xyyz_0, \
                             g_0_yyzzzz_0_xyzz_0, \
                             g_0_yyzzzz_0_xzzz_0, \
                             g_0_yyzzzz_0_yyyy_0, \
                             g_0_yyzzzz_0_yyyz_0, \
                             g_0_yyzzzz_0_yyzz_0, \
                             g_0_yyzzzz_0_yzzz_0, \
                             g_0_yyzzzz_0_zzzz_0, \
                             g_0_yzzzz_0_xxxx_0,  \
                             g_0_yzzzz_0_xxxx_1,  \
                             g_0_yzzzz_0_xxxz_0,  \
                             g_0_yzzzz_0_xxxz_1,  \
                             g_0_yzzzz_0_xxyz_0,  \
                             g_0_yzzzz_0_xxyz_1,  \
                             g_0_yzzzz_0_xxz_1,   \
                             g_0_yzzzz_0_xxzz_0,  \
                             g_0_yzzzz_0_xxzz_1,  \
                             g_0_yzzzz_0_xyyz_0,  \
                             g_0_yzzzz_0_xyyz_1,  \
                             g_0_yzzzz_0_xyz_1,   \
                             g_0_yzzzz_0_xyzz_0,  \
                             g_0_yzzzz_0_xyzz_1,  \
                             g_0_yzzzz_0_xzz_1,   \
                             g_0_yzzzz_0_xzzz_0,  \
                             g_0_yzzzz_0_xzzz_1,  \
                             g_0_yzzzz_0_yyyz_0,  \
                             g_0_yzzzz_0_yyyz_1,  \
                             g_0_yzzzz_0_yyz_1,   \
                             g_0_yzzzz_0_yyzz_0,  \
                             g_0_yzzzz_0_yyzz_1,  \
                             g_0_yzzzz_0_yzz_1,   \
                             g_0_yzzzz_0_yzzz_0,  \
                             g_0_yzzzz_0_yzzz_1,  \
                             g_0_yzzzz_0_zzz_1,   \
                             g_0_yzzzz_0_zzzz_0,  \
                             g_0_yzzzz_0_zzzz_1,  \
                             g_0_zzzz_0_xxxx_0,   \
                             g_0_zzzz_0_xxxx_1,   \
                             g_0_zzzz_0_xxxz_0,   \
                             g_0_zzzz_0_xxxz_1,   \
                             g_0_zzzz_0_xxyz_0,   \
                             g_0_zzzz_0_xxyz_1,   \
                             g_0_zzzz_0_xxzz_0,   \
                             g_0_zzzz_0_xxzz_1,   \
                             g_0_zzzz_0_xyyz_0,   \
                             g_0_zzzz_0_xyyz_1,   \
                             g_0_zzzz_0_xyzz_0,   \
                             g_0_zzzz_0_xyzz_1,   \
                             g_0_zzzz_0_xzzz_0,   \
                             g_0_zzzz_0_xzzz_1,   \
                             g_0_zzzz_0_yyyz_0,   \
                             g_0_zzzz_0_yyyz_1,   \
                             g_0_zzzz_0_yyzz_0,   \
                             g_0_zzzz_0_yyzz_1,   \
                             g_0_zzzz_0_yzzz_0,   \
                             g_0_zzzz_0_yzzz_1,   \
                             g_0_zzzz_0_zzzz_0,   \
                             g_0_zzzz_0_zzzz_1,   \
                             wp_y,                \
                             wp_z,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyzzzz_0_xxxx_0[i] =
            g_0_zzzz_0_xxxx_0[i] * fi_ab_0 - g_0_zzzz_0_xxxx_1[i] * fti_ab_0 + g_0_yzzzz_0_xxxx_0[i] * pb_y + g_0_yzzzz_0_xxxx_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxxy_0[i] = 3.0 * g_0_yyzz_0_xxxy_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxy_1[i] * fti_ab_0 + g_0_yyzzz_0_xxxy_0[i] * pb_z +
                                 g_0_yyzzz_0_xxxy_1[i] * wp_z[i];

        g_0_yyzzzz_0_xxxz_0[i] =
            g_0_zzzz_0_xxxz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxz_1[i] * fti_ab_0 + g_0_yzzzz_0_xxxz_0[i] * pb_y + g_0_yzzzz_0_xxxz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxyy_0[i] = 3.0 * g_0_yyzz_0_xxyy_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxyy_1[i] * fti_ab_0 + g_0_yyzzz_0_xxyy_0[i] * pb_z +
                                 g_0_yyzzz_0_xxyy_1[i] * wp_z[i];

        g_0_yyzzzz_0_xxyz_0[i] = g_0_zzzz_0_xxyz_0[i] * fi_ab_0 - g_0_zzzz_0_xxyz_1[i] * fti_ab_0 + g_0_yzzzz_0_xxz_1[i] * fi_abcd_0 +
                                 g_0_yzzzz_0_xxyz_0[i] * pb_y + g_0_yzzzz_0_xxyz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxzz_0[i] =
            g_0_zzzz_0_xxzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxzz_1[i] * fti_ab_0 + g_0_yzzzz_0_xxzz_0[i] * pb_y + g_0_yzzzz_0_xxzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xyyy_0[i] = 3.0 * g_0_yyzz_0_xyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xyyy_1[i] * fti_ab_0 + g_0_yyzzz_0_xyyy_0[i] * pb_z +
                                 g_0_yyzzz_0_xyyy_1[i] * wp_z[i];

        g_0_yyzzzz_0_xyyz_0[i] = g_0_zzzz_0_xyyz_0[i] * fi_ab_0 - g_0_zzzz_0_xyyz_1[i] * fti_ab_0 + 2.0 * g_0_yzzzz_0_xyz_1[i] * fi_abcd_0 +
                                 g_0_yzzzz_0_xyyz_0[i] * pb_y + g_0_yzzzz_0_xyyz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xyzz_0[i] = g_0_zzzz_0_xyzz_0[i] * fi_ab_0 - g_0_zzzz_0_xyzz_1[i] * fti_ab_0 + g_0_yzzzz_0_xzz_1[i] * fi_abcd_0 +
                                 g_0_yzzzz_0_xyzz_0[i] * pb_y + g_0_yzzzz_0_xyzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xzzz_0[i] =
            g_0_zzzz_0_xzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xzzz_1[i] * fti_ab_0 + g_0_yzzzz_0_xzzz_0[i] * pb_y + g_0_yzzzz_0_xzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_yyyy_0[i] = 3.0 * g_0_yyzz_0_yyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_yyyy_1[i] * fti_ab_0 + g_0_yyzzz_0_yyyy_0[i] * pb_z +
                                 g_0_yyzzz_0_yyyy_1[i] * wp_z[i];

        g_0_yyzzzz_0_yyyz_0[i] = g_0_zzzz_0_yyyz_0[i] * fi_ab_0 - g_0_zzzz_0_yyyz_1[i] * fti_ab_0 + 3.0 * g_0_yzzzz_0_yyz_1[i] * fi_abcd_0 +
                                 g_0_yzzzz_0_yyyz_0[i] * pb_y + g_0_yzzzz_0_yyyz_1[i] * wp_y[i];

        g_0_yyzzzz_0_yyzz_0[i] = g_0_zzzz_0_yyzz_0[i] * fi_ab_0 - g_0_zzzz_0_yyzz_1[i] * fti_ab_0 + 2.0 * g_0_yzzzz_0_yzz_1[i] * fi_abcd_0 +
                                 g_0_yzzzz_0_yyzz_0[i] * pb_y + g_0_yzzzz_0_yyzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_yzzz_0[i] = g_0_zzzz_0_yzzz_0[i] * fi_ab_0 - g_0_zzzz_0_yzzz_1[i] * fti_ab_0 + g_0_yzzzz_0_zzz_1[i] * fi_abcd_0 +
                                 g_0_yzzzz_0_yzzz_0[i] * pb_y + g_0_yzzzz_0_yzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_zzzz_0[i] =
            g_0_zzzz_0_zzzz_0[i] * fi_ab_0 - g_0_zzzz_0_zzzz_1[i] * fti_ab_0 + g_0_yzzzz_0_zzzz_0[i] * pb_y + g_0_yzzzz_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 390-405 components of targeted buffer : SISG

    auto g_0_yzzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 390);

    auto g_0_yzzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 391);

    auto g_0_yzzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 392);

    auto g_0_yzzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 393);

    auto g_0_yzzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 394);

    auto g_0_yzzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 395);

    auto g_0_yzzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 396);

    auto g_0_yzzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 397);

    auto g_0_yzzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 398);

    auto g_0_yzzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 399);

    auto g_0_yzzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 400);

    auto g_0_yzzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 401);

    auto g_0_yzzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 402);

    auto g_0_yzzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 403);

    auto g_0_yzzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 404);

#pragma omp simd aligned(g_0_yzzzzz_0_xxxx_0,     \
                             g_0_yzzzzz_0_xxxy_0, \
                             g_0_yzzzzz_0_xxxz_0, \
                             g_0_yzzzzz_0_xxyy_0, \
                             g_0_yzzzzz_0_xxyz_0, \
                             g_0_yzzzzz_0_xxzz_0, \
                             g_0_yzzzzz_0_xyyy_0, \
                             g_0_yzzzzz_0_xyyz_0, \
                             g_0_yzzzzz_0_xyzz_0, \
                             g_0_yzzzzz_0_xzzz_0, \
                             g_0_yzzzzz_0_yyyy_0, \
                             g_0_yzzzzz_0_yyyz_0, \
                             g_0_yzzzzz_0_yyzz_0, \
                             g_0_yzzzzz_0_yzzz_0, \
                             g_0_yzzzzz_0_zzzz_0, \
                             g_0_zzzzz_0_xxx_1,   \
                             g_0_zzzzz_0_xxxx_0,  \
                             g_0_zzzzz_0_xxxx_1,  \
                             g_0_zzzzz_0_xxxy_0,  \
                             g_0_zzzzz_0_xxxy_1,  \
                             g_0_zzzzz_0_xxxz_0,  \
                             g_0_zzzzz_0_xxxz_1,  \
                             g_0_zzzzz_0_xxy_1,   \
                             g_0_zzzzz_0_xxyy_0,  \
                             g_0_zzzzz_0_xxyy_1,  \
                             g_0_zzzzz_0_xxyz_0,  \
                             g_0_zzzzz_0_xxyz_1,  \
                             g_0_zzzzz_0_xxz_1,   \
                             g_0_zzzzz_0_xxzz_0,  \
                             g_0_zzzzz_0_xxzz_1,  \
                             g_0_zzzzz_0_xyy_1,   \
                             g_0_zzzzz_0_xyyy_0,  \
                             g_0_zzzzz_0_xyyy_1,  \
                             g_0_zzzzz_0_xyyz_0,  \
                             g_0_zzzzz_0_xyyz_1,  \
                             g_0_zzzzz_0_xyz_1,   \
                             g_0_zzzzz_0_xyzz_0,  \
                             g_0_zzzzz_0_xyzz_1,  \
                             g_0_zzzzz_0_xzz_1,   \
                             g_0_zzzzz_0_xzzz_0,  \
                             g_0_zzzzz_0_xzzz_1,  \
                             g_0_zzzzz_0_yyy_1,   \
                             g_0_zzzzz_0_yyyy_0,  \
                             g_0_zzzzz_0_yyyy_1,  \
                             g_0_zzzzz_0_yyyz_0,  \
                             g_0_zzzzz_0_yyyz_1,  \
                             g_0_zzzzz_0_yyz_1,   \
                             g_0_zzzzz_0_yyzz_0,  \
                             g_0_zzzzz_0_yyzz_1,  \
                             g_0_zzzzz_0_yzz_1,   \
                             g_0_zzzzz_0_yzzz_0,  \
                             g_0_zzzzz_0_yzzz_1,  \
                             g_0_zzzzz_0_zzz_1,   \
                             g_0_zzzzz_0_zzzz_0,  \
                             g_0_zzzzz_0_zzzz_1,  \
                             wp_y,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzzzzz_0_xxxx_0[i] = g_0_zzzzz_0_xxxx_0[i] * pb_y + g_0_zzzzz_0_xxxx_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxy_0[i] = g_0_zzzzz_0_xxx_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxy_0[i] * pb_y + g_0_zzzzz_0_xxxy_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxz_0[i] = g_0_zzzzz_0_xxxz_0[i] * pb_y + g_0_zzzzz_0_xxxz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxyy_0[i] = 2.0 * g_0_zzzzz_0_xxy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyy_0[i] * pb_y + g_0_zzzzz_0_xxyy_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxyz_0[i] = g_0_zzzzz_0_xxz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyz_0[i] * pb_y + g_0_zzzzz_0_xxyz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxzz_0[i] = g_0_zzzzz_0_xxzz_0[i] * pb_y + g_0_zzzzz_0_xxzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xyyy_0[i] = 3.0 * g_0_zzzzz_0_xyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyy_0[i] * pb_y + g_0_zzzzz_0_xyyy_1[i] * wp_y[i];

        g_0_yzzzzz_0_xyyz_0[i] = 2.0 * g_0_zzzzz_0_xyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyz_0[i] * pb_y + g_0_zzzzz_0_xyyz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xyzz_0[i] = g_0_zzzzz_0_xzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyzz_0[i] * pb_y + g_0_zzzzz_0_xyzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xzzz_0[i] = g_0_zzzzz_0_xzzz_0[i] * pb_y + g_0_zzzzz_0_xzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_yyyy_0[i] = 4.0 * g_0_zzzzz_0_yyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyyy_0[i] * pb_y + g_0_zzzzz_0_yyyy_1[i] * wp_y[i];

        g_0_yzzzzz_0_yyyz_0[i] = 3.0 * g_0_zzzzz_0_yyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyyz_0[i] * pb_y + g_0_zzzzz_0_yyyz_1[i] * wp_y[i];

        g_0_yzzzzz_0_yyzz_0[i] = 2.0 * g_0_zzzzz_0_yzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyzz_0[i] * pb_y + g_0_zzzzz_0_yyzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_yzzz_0[i] = g_0_zzzzz_0_zzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yzzz_0[i] * pb_y + g_0_zzzzz_0_yzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_zzzz_0[i] = g_0_zzzzz_0_zzzz_0[i] * pb_y + g_0_zzzzz_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 405-420 components of targeted buffer : SISG

    auto g_0_zzzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sisg + 405);

    auto g_0_zzzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sisg + 406);

    auto g_0_zzzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sisg + 407);

    auto g_0_zzzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sisg + 408);

    auto g_0_zzzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sisg + 409);

    auto g_0_zzzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sisg + 410);

    auto g_0_zzzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sisg + 411);

    auto g_0_zzzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sisg + 412);

    auto g_0_zzzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sisg + 413);

    auto g_0_zzzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sisg + 414);

    auto g_0_zzzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sisg + 415);

    auto g_0_zzzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sisg + 416);

    auto g_0_zzzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sisg + 417);

    auto g_0_zzzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sisg + 418);

    auto g_0_zzzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sisg + 419);

#pragma omp simd aligned(g_0_zzzz_0_xxxx_0,       \
                             g_0_zzzz_0_xxxx_1,   \
                             g_0_zzzz_0_xxxy_0,   \
                             g_0_zzzz_0_xxxy_1,   \
                             g_0_zzzz_0_xxxz_0,   \
                             g_0_zzzz_0_xxxz_1,   \
                             g_0_zzzz_0_xxyy_0,   \
                             g_0_zzzz_0_xxyy_1,   \
                             g_0_zzzz_0_xxyz_0,   \
                             g_0_zzzz_0_xxyz_1,   \
                             g_0_zzzz_0_xxzz_0,   \
                             g_0_zzzz_0_xxzz_1,   \
                             g_0_zzzz_0_xyyy_0,   \
                             g_0_zzzz_0_xyyy_1,   \
                             g_0_zzzz_0_xyyz_0,   \
                             g_0_zzzz_0_xyyz_1,   \
                             g_0_zzzz_0_xyzz_0,   \
                             g_0_zzzz_0_xyzz_1,   \
                             g_0_zzzz_0_xzzz_0,   \
                             g_0_zzzz_0_xzzz_1,   \
                             g_0_zzzz_0_yyyy_0,   \
                             g_0_zzzz_0_yyyy_1,   \
                             g_0_zzzz_0_yyyz_0,   \
                             g_0_zzzz_0_yyyz_1,   \
                             g_0_zzzz_0_yyzz_0,   \
                             g_0_zzzz_0_yyzz_1,   \
                             g_0_zzzz_0_yzzz_0,   \
                             g_0_zzzz_0_yzzz_1,   \
                             g_0_zzzz_0_zzzz_0,   \
                             g_0_zzzz_0_zzzz_1,   \
                             g_0_zzzzz_0_xxx_1,   \
                             g_0_zzzzz_0_xxxx_0,  \
                             g_0_zzzzz_0_xxxx_1,  \
                             g_0_zzzzz_0_xxxy_0,  \
                             g_0_zzzzz_0_xxxy_1,  \
                             g_0_zzzzz_0_xxxz_0,  \
                             g_0_zzzzz_0_xxxz_1,  \
                             g_0_zzzzz_0_xxy_1,   \
                             g_0_zzzzz_0_xxyy_0,  \
                             g_0_zzzzz_0_xxyy_1,  \
                             g_0_zzzzz_0_xxyz_0,  \
                             g_0_zzzzz_0_xxyz_1,  \
                             g_0_zzzzz_0_xxz_1,   \
                             g_0_zzzzz_0_xxzz_0,  \
                             g_0_zzzzz_0_xxzz_1,  \
                             g_0_zzzzz_0_xyy_1,   \
                             g_0_zzzzz_0_xyyy_0,  \
                             g_0_zzzzz_0_xyyy_1,  \
                             g_0_zzzzz_0_xyyz_0,  \
                             g_0_zzzzz_0_xyyz_1,  \
                             g_0_zzzzz_0_xyz_1,   \
                             g_0_zzzzz_0_xyzz_0,  \
                             g_0_zzzzz_0_xyzz_1,  \
                             g_0_zzzzz_0_xzz_1,   \
                             g_0_zzzzz_0_xzzz_0,  \
                             g_0_zzzzz_0_xzzz_1,  \
                             g_0_zzzzz_0_yyy_1,   \
                             g_0_zzzzz_0_yyyy_0,  \
                             g_0_zzzzz_0_yyyy_1,  \
                             g_0_zzzzz_0_yyyz_0,  \
                             g_0_zzzzz_0_yyyz_1,  \
                             g_0_zzzzz_0_yyz_1,   \
                             g_0_zzzzz_0_yyzz_0,  \
                             g_0_zzzzz_0_yyzz_1,  \
                             g_0_zzzzz_0_yzz_1,   \
                             g_0_zzzzz_0_yzzz_0,  \
                             g_0_zzzzz_0_yzzz_1,  \
                             g_0_zzzzz_0_zzz_1,   \
                             g_0_zzzzz_0_zzzz_0,  \
                             g_0_zzzzz_0_zzzz_1,  \
                             g_0_zzzzzz_0_xxxx_0, \
                             g_0_zzzzzz_0_xxxy_0, \
                             g_0_zzzzzz_0_xxxz_0, \
                             g_0_zzzzzz_0_xxyy_0, \
                             g_0_zzzzzz_0_xxyz_0, \
                             g_0_zzzzzz_0_xxzz_0, \
                             g_0_zzzzzz_0_xyyy_0, \
                             g_0_zzzzzz_0_xyyz_0, \
                             g_0_zzzzzz_0_xyzz_0, \
                             g_0_zzzzzz_0_xzzz_0, \
                             g_0_zzzzzz_0_yyyy_0, \
                             g_0_zzzzzz_0_yyyz_0, \
                             g_0_zzzzzz_0_yyzz_0, \
                             g_0_zzzzzz_0_yzzz_0, \
                             g_0_zzzzzz_0_zzzz_0, \
                             wp_z,                \
                             c_exps,              \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzzzzz_0_xxxx_0[i] = 5.0 * g_0_zzzz_0_xxxx_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxx_1[i] * fti_ab_0 + g_0_zzzzz_0_xxxx_0[i] * pb_z +
                                 g_0_zzzzz_0_xxxx_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxy_0[i] = 5.0 * g_0_zzzz_0_xxxy_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxy_1[i] * fti_ab_0 + g_0_zzzzz_0_xxxy_0[i] * pb_z +
                                 g_0_zzzzz_0_xxxy_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxz_0[i] = 5.0 * g_0_zzzz_0_xxxz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxz_1[i] * fti_ab_0 + g_0_zzzzz_0_xxx_1[i] * fi_abcd_0 +
                                 g_0_zzzzz_0_xxxz_0[i] * pb_z + g_0_zzzzz_0_xxxz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxyy_0[i] = 5.0 * g_0_zzzz_0_xxyy_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxyy_1[i] * fti_ab_0 + g_0_zzzzz_0_xxyy_0[i] * pb_z +
                                 g_0_zzzzz_0_xxyy_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxyz_0[i] = 5.0 * g_0_zzzz_0_xxyz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxyz_1[i] * fti_ab_0 + g_0_zzzzz_0_xxy_1[i] * fi_abcd_0 +
                                 g_0_zzzzz_0_xxyz_0[i] * pb_z + g_0_zzzzz_0_xxyz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxzz_0[i] = 5.0 * g_0_zzzz_0_xxzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxzz_1[i] * fti_ab_0 +
                                 2.0 * g_0_zzzzz_0_xxz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxzz_0[i] * pb_z + g_0_zzzzz_0_xxzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xyyy_0[i] = 5.0 * g_0_zzzz_0_xyyy_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xyyy_1[i] * fti_ab_0 + g_0_zzzzz_0_xyyy_0[i] * pb_z +
                                 g_0_zzzzz_0_xyyy_1[i] * wp_z[i];

        g_0_zzzzzz_0_xyyz_0[i] = 5.0 * g_0_zzzz_0_xyyz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xyyz_1[i] * fti_ab_0 + g_0_zzzzz_0_xyy_1[i] * fi_abcd_0 +
                                 g_0_zzzzz_0_xyyz_0[i] * pb_z + g_0_zzzzz_0_xyyz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xyzz_0[i] = 5.0 * g_0_zzzz_0_xyzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xyzz_1[i] * fti_ab_0 +
                                 2.0 * g_0_zzzzz_0_xyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyzz_0[i] * pb_z + g_0_zzzzz_0_xyzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xzzz_0[i] = 5.0 * g_0_zzzz_0_xzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xzzz_1[i] * fti_ab_0 +
                                 3.0 * g_0_zzzzz_0_xzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xzzz_0[i] * pb_z + g_0_zzzzz_0_xzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_yyyy_0[i] = 5.0 * g_0_zzzz_0_yyyy_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_yyyy_1[i] * fti_ab_0 + g_0_zzzzz_0_yyyy_0[i] * pb_z +
                                 g_0_zzzzz_0_yyyy_1[i] * wp_z[i];

        g_0_zzzzzz_0_yyyz_0[i] = 5.0 * g_0_zzzz_0_yyyz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_yyyz_1[i] * fti_ab_0 + g_0_zzzzz_0_yyy_1[i] * fi_abcd_0 +
                                 g_0_zzzzz_0_yyyz_0[i] * pb_z + g_0_zzzzz_0_yyyz_1[i] * wp_z[i];

        g_0_zzzzzz_0_yyzz_0[i] = 5.0 * g_0_zzzz_0_yyzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_yyzz_1[i] * fti_ab_0 +
                                 2.0 * g_0_zzzzz_0_yyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyzz_0[i] * pb_z + g_0_zzzzz_0_yyzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_yzzz_0[i] = 5.0 * g_0_zzzz_0_yzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_yzzz_1[i] * fti_ab_0 +
                                 3.0 * g_0_zzzzz_0_yzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yzzz_0[i] * pb_z + g_0_zzzzz_0_yzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_zzzz_0[i] = 5.0 * g_0_zzzz_0_zzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_zzzz_1[i] * fti_ab_0 +
                                 4.0 * g_0_zzzzz_0_zzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_zzzz_0[i] * pb_z + g_0_zzzzz_0_zzzz_1[i] * wp_z[i];
    }
}

}  // namespace erirec

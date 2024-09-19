#include "ElectronRepulsionPrimRecSGSH.hpp"

namespace erirec {  // erirec namespace

auto
comp_prim_electron_repulsion_sgsh(CSimdArray<double>&   pbuffer,
                                  const size_t          idx_eri_0_sgsh,
                                  size_t                idx_eri_0_sdsh,
                                  size_t                idx_eri_1_sdsh,
                                  size_t                idx_eri_1_sfsg,
                                  size_t                idx_eri_0_sfsh,
                                  size_t                idx_eri_1_sfsh,
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

    /// Set up components of auxilary buffer : SFSG

    auto g_0_xxx_0_xxxx_1 = pbuffer.data(idx_eri_1_sfsg);

    auto g_0_xxx_0_xxxy_1 = pbuffer.data(idx_eri_1_sfsg + 1);

    auto g_0_xxx_0_xxxz_1 = pbuffer.data(idx_eri_1_sfsg + 2);

    auto g_0_xxx_0_xxyy_1 = pbuffer.data(idx_eri_1_sfsg + 3);

    auto g_0_xxx_0_xxyz_1 = pbuffer.data(idx_eri_1_sfsg + 4);

    auto g_0_xxx_0_xxzz_1 = pbuffer.data(idx_eri_1_sfsg + 5);

    auto g_0_xxx_0_xyyy_1 = pbuffer.data(idx_eri_1_sfsg + 6);

    auto g_0_xxx_0_xyyz_1 = pbuffer.data(idx_eri_1_sfsg + 7);

    auto g_0_xxx_0_xyzz_1 = pbuffer.data(idx_eri_1_sfsg + 8);

    auto g_0_xxx_0_xzzz_1 = pbuffer.data(idx_eri_1_sfsg + 9);

    auto g_0_xxx_0_yyyy_1 = pbuffer.data(idx_eri_1_sfsg + 10);

    auto g_0_xxx_0_yyyz_1 = pbuffer.data(idx_eri_1_sfsg + 11);

    auto g_0_xxx_0_yyzz_1 = pbuffer.data(idx_eri_1_sfsg + 12);

    auto g_0_xxx_0_yzzz_1 = pbuffer.data(idx_eri_1_sfsg + 13);

    auto g_0_xxx_0_zzzz_1 = pbuffer.data(idx_eri_1_sfsg + 14);

    auto g_0_xxz_0_xxxz_1 = pbuffer.data(idx_eri_1_sfsg + 32);

    auto g_0_xxz_0_xxyz_1 = pbuffer.data(idx_eri_1_sfsg + 34);

    auto g_0_xxz_0_xxzz_1 = pbuffer.data(idx_eri_1_sfsg + 35);

    auto g_0_xxz_0_xyyz_1 = pbuffer.data(idx_eri_1_sfsg + 37);

    auto g_0_xxz_0_xyzz_1 = pbuffer.data(idx_eri_1_sfsg + 38);

    auto g_0_xxz_0_xzzz_1 = pbuffer.data(idx_eri_1_sfsg + 39);

    auto g_0_xxz_0_yyyz_1 = pbuffer.data(idx_eri_1_sfsg + 41);

    auto g_0_xxz_0_yyzz_1 = pbuffer.data(idx_eri_1_sfsg + 42);

    auto g_0_xxz_0_yzzz_1 = pbuffer.data(idx_eri_1_sfsg + 43);

    auto g_0_xxz_0_zzzz_1 = pbuffer.data(idx_eri_1_sfsg + 44);

    auto g_0_xyy_0_xxxy_1 = pbuffer.data(idx_eri_1_sfsg + 46);

    auto g_0_xyy_0_xxyy_1 = pbuffer.data(idx_eri_1_sfsg + 48);

    auto g_0_xyy_0_xxyz_1 = pbuffer.data(idx_eri_1_sfsg + 49);

    auto g_0_xyy_0_xyyy_1 = pbuffer.data(idx_eri_1_sfsg + 51);

    auto g_0_xyy_0_xyyz_1 = pbuffer.data(idx_eri_1_sfsg + 52);

    auto g_0_xyy_0_xyzz_1 = pbuffer.data(idx_eri_1_sfsg + 53);

    auto g_0_xyy_0_yyyy_1 = pbuffer.data(idx_eri_1_sfsg + 55);

    auto g_0_xyy_0_yyyz_1 = pbuffer.data(idx_eri_1_sfsg + 56);

    auto g_0_xyy_0_yyzz_1 = pbuffer.data(idx_eri_1_sfsg + 57);

    auto g_0_xyy_0_yzzz_1 = pbuffer.data(idx_eri_1_sfsg + 58);

    auto g_0_xzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sfsg + 77);

    auto g_0_xzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sfsg + 79);

    auto g_0_xzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sfsg + 80);

    auto g_0_xzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sfsg + 82);

    auto g_0_xzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sfsg + 83);

    auto g_0_xzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sfsg + 84);

    auto g_0_xzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sfsg + 86);

    auto g_0_xzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sfsg + 87);

    auto g_0_xzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sfsg + 88);

    auto g_0_xzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sfsg + 89);

    auto g_0_yyy_0_xxxx_1 = pbuffer.data(idx_eri_1_sfsg + 90);

    auto g_0_yyy_0_xxxy_1 = pbuffer.data(idx_eri_1_sfsg + 91);

    auto g_0_yyy_0_xxxz_1 = pbuffer.data(idx_eri_1_sfsg + 92);

    auto g_0_yyy_0_xxyy_1 = pbuffer.data(idx_eri_1_sfsg + 93);

    auto g_0_yyy_0_xxyz_1 = pbuffer.data(idx_eri_1_sfsg + 94);

    auto g_0_yyy_0_xxzz_1 = pbuffer.data(idx_eri_1_sfsg + 95);

    auto g_0_yyy_0_xyyy_1 = pbuffer.data(idx_eri_1_sfsg + 96);

    auto g_0_yyy_0_xyyz_1 = pbuffer.data(idx_eri_1_sfsg + 97);

    auto g_0_yyy_0_xyzz_1 = pbuffer.data(idx_eri_1_sfsg + 98);

    auto g_0_yyy_0_xzzz_1 = pbuffer.data(idx_eri_1_sfsg + 99);

    auto g_0_yyy_0_yyyy_1 = pbuffer.data(idx_eri_1_sfsg + 100);

    auto g_0_yyy_0_yyyz_1 = pbuffer.data(idx_eri_1_sfsg + 101);

    auto g_0_yyy_0_yyzz_1 = pbuffer.data(idx_eri_1_sfsg + 102);

    auto g_0_yyy_0_yzzz_1 = pbuffer.data(idx_eri_1_sfsg + 103);

    auto g_0_yyy_0_zzzz_1 = pbuffer.data(idx_eri_1_sfsg + 104);

    auto g_0_yyz_0_xxxz_1 = pbuffer.data(idx_eri_1_sfsg + 107);

    auto g_0_yyz_0_xxyz_1 = pbuffer.data(idx_eri_1_sfsg + 109);

    auto g_0_yyz_0_xxzz_1 = pbuffer.data(idx_eri_1_sfsg + 110);

    auto g_0_yyz_0_xyyz_1 = pbuffer.data(idx_eri_1_sfsg + 112);

    auto g_0_yyz_0_xyzz_1 = pbuffer.data(idx_eri_1_sfsg + 113);

    auto g_0_yyz_0_xzzz_1 = pbuffer.data(idx_eri_1_sfsg + 114);

    auto g_0_yyz_0_yyyz_1 = pbuffer.data(idx_eri_1_sfsg + 116);

    auto g_0_yyz_0_yyzz_1 = pbuffer.data(idx_eri_1_sfsg + 117);

    auto g_0_yyz_0_yzzz_1 = pbuffer.data(idx_eri_1_sfsg + 118);

    auto g_0_yyz_0_zzzz_1 = pbuffer.data(idx_eri_1_sfsg + 119);

    auto g_0_yzz_0_xxxy_1 = pbuffer.data(idx_eri_1_sfsg + 121);

    auto g_0_yzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sfsg + 122);

    auto g_0_yzz_0_xxyy_1 = pbuffer.data(idx_eri_1_sfsg + 123);

    auto g_0_yzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sfsg + 124);

    auto g_0_yzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sfsg + 125);

    auto g_0_yzz_0_xyyy_1 = pbuffer.data(idx_eri_1_sfsg + 126);

    auto g_0_yzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sfsg + 127);

    auto g_0_yzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sfsg + 128);

    auto g_0_yzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sfsg + 129);

    auto g_0_yzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sfsg + 130);

    auto g_0_yzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sfsg + 131);

    auto g_0_yzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sfsg + 132);

    auto g_0_yzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sfsg + 133);

    auto g_0_yzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sfsg + 134);

    auto g_0_zzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sfsg + 135);

    auto g_0_zzz_0_xxxy_1 = pbuffer.data(idx_eri_1_sfsg + 136);

    auto g_0_zzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sfsg + 137);

    auto g_0_zzz_0_xxyy_1 = pbuffer.data(idx_eri_1_sfsg + 138);

    auto g_0_zzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sfsg + 139);

    auto g_0_zzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sfsg + 140);

    auto g_0_zzz_0_xyyy_1 = pbuffer.data(idx_eri_1_sfsg + 141);

    auto g_0_zzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sfsg + 142);

    auto g_0_zzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sfsg + 143);

    auto g_0_zzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sfsg + 144);

    auto g_0_zzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sfsg + 145);

    auto g_0_zzz_0_yyyz_1 = pbuffer.data(idx_eri_1_sfsg + 146);

    auto g_0_zzz_0_yyzz_1 = pbuffer.data(idx_eri_1_sfsg + 147);

    auto g_0_zzz_0_yzzz_1 = pbuffer.data(idx_eri_1_sfsg + 148);

    auto g_0_zzz_0_zzzz_1 = pbuffer.data(idx_eri_1_sfsg + 149);

    /// Set up components of auxilary buffer : SFSH

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

    auto g_0_xxy_0_xxxxx_0 = pbuffer.data(idx_eri_0_sfsh + 21);

    auto g_0_xxy_0_xxxxy_0 = pbuffer.data(idx_eri_0_sfsh + 22);

    auto g_0_xxy_0_xxxxz_0 = pbuffer.data(idx_eri_0_sfsh + 23);

    auto g_0_xxy_0_xxxyy_0 = pbuffer.data(idx_eri_0_sfsh + 24);

    auto g_0_xxy_0_xxxzz_0 = pbuffer.data(idx_eri_0_sfsh + 26);

    auto g_0_xxy_0_xxyyy_0 = pbuffer.data(idx_eri_0_sfsh + 27);

    auto g_0_xxy_0_xxzzz_0 = pbuffer.data(idx_eri_0_sfsh + 30);

    auto g_0_xxy_0_xyyyy_0 = pbuffer.data(idx_eri_0_sfsh + 31);

    auto g_0_xxy_0_xzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 35);

    auto g_0_xxy_0_yyyyy_0 = pbuffer.data(idx_eri_0_sfsh + 36);

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

    auto g_0_xxz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sfsh + 58);

    auto g_0_xxz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sfsh + 59);

    auto g_0_xxz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sfsh + 60);

    auto g_0_xxz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 61);

    auto g_0_xxz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 62);

    auto g_0_xyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_sfsh + 63);

    auto g_0_xyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_sfsh + 64);

    auto g_0_xyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_sfsh + 66);

    auto g_0_xyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_sfsh + 67);

    auto g_0_xyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_sfsh + 69);

    auto g_0_xyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_sfsh + 70);

    auto g_0_xyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_sfsh + 71);

    auto g_0_xyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_sfsh + 73);

    auto g_0_xyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_sfsh + 74);

    auto g_0_xyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_sfsh + 75);

    auto g_0_xyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_sfsh + 76);

    auto g_0_xyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_sfsh + 78);

    auto g_0_xyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_sfsh + 79);

    auto g_0_xyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_sfsh + 80);

    auto g_0_xyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_sfsh + 81);

    auto g_0_xyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 82);

    auto g_0_xyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_sfsh + 83);

    auto g_0_xzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sfsh + 105);

    auto g_0_xzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sfsh + 107);

    auto g_0_xzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sfsh + 109);

    auto g_0_xzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sfsh + 110);

    auto g_0_xzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sfsh + 112);

    auto g_0_xzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sfsh + 113);

    auto g_0_xzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sfsh + 114);

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

    /// Set up components of auxilary buffer : SFSH

    auto g_0_xxx_0_xxxxx_1 = pbuffer.data(idx_eri_1_sfsh);

    auto g_0_xxx_0_xxxxy_1 = pbuffer.data(idx_eri_1_sfsh + 1);

    auto g_0_xxx_0_xxxxz_1 = pbuffer.data(idx_eri_1_sfsh + 2);

    auto g_0_xxx_0_xxxyy_1 = pbuffer.data(idx_eri_1_sfsh + 3);

    auto g_0_xxx_0_xxxyz_1 = pbuffer.data(idx_eri_1_sfsh + 4);

    auto g_0_xxx_0_xxxzz_1 = pbuffer.data(idx_eri_1_sfsh + 5);

    auto g_0_xxx_0_xxyyy_1 = pbuffer.data(idx_eri_1_sfsh + 6);

    auto g_0_xxx_0_xxyyz_1 = pbuffer.data(idx_eri_1_sfsh + 7);

    auto g_0_xxx_0_xxyzz_1 = pbuffer.data(idx_eri_1_sfsh + 8);

    auto g_0_xxx_0_xxzzz_1 = pbuffer.data(idx_eri_1_sfsh + 9);

    auto g_0_xxx_0_xyyyy_1 = pbuffer.data(idx_eri_1_sfsh + 10);

    auto g_0_xxx_0_xyyyz_1 = pbuffer.data(idx_eri_1_sfsh + 11);

    auto g_0_xxx_0_xyyzz_1 = pbuffer.data(idx_eri_1_sfsh + 12);

    auto g_0_xxx_0_xyzzz_1 = pbuffer.data(idx_eri_1_sfsh + 13);

    auto g_0_xxx_0_xzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 14);

    auto g_0_xxx_0_yyyyy_1 = pbuffer.data(idx_eri_1_sfsh + 15);

    auto g_0_xxx_0_yyyyz_1 = pbuffer.data(idx_eri_1_sfsh + 16);

    auto g_0_xxx_0_yyyzz_1 = pbuffer.data(idx_eri_1_sfsh + 17);

    auto g_0_xxx_0_yyzzz_1 = pbuffer.data(idx_eri_1_sfsh + 18);

    auto g_0_xxx_0_yzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 19);

    auto g_0_xxx_0_zzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 20);

    auto g_0_xxy_0_xxxxx_1 = pbuffer.data(idx_eri_1_sfsh + 21);

    auto g_0_xxy_0_xxxxy_1 = pbuffer.data(idx_eri_1_sfsh + 22);

    auto g_0_xxy_0_xxxxz_1 = pbuffer.data(idx_eri_1_sfsh + 23);

    auto g_0_xxy_0_xxxyy_1 = pbuffer.data(idx_eri_1_sfsh + 24);

    auto g_0_xxy_0_xxxzz_1 = pbuffer.data(idx_eri_1_sfsh + 26);

    auto g_0_xxy_0_xxyyy_1 = pbuffer.data(idx_eri_1_sfsh + 27);

    auto g_0_xxy_0_xxzzz_1 = pbuffer.data(idx_eri_1_sfsh + 30);

    auto g_0_xxy_0_xyyyy_1 = pbuffer.data(idx_eri_1_sfsh + 31);

    auto g_0_xxy_0_xzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 35);

    auto g_0_xxy_0_yyyyy_1 = pbuffer.data(idx_eri_1_sfsh + 36);

    auto g_0_xxz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sfsh + 42);

    auto g_0_xxz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sfsh + 43);

    auto g_0_xxz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sfsh + 44);

    auto g_0_xxz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sfsh + 45);

    auto g_0_xxz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sfsh + 46);

    auto g_0_xxz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sfsh + 47);

    auto g_0_xxz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sfsh + 48);

    auto g_0_xxz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sfsh + 49);

    auto g_0_xxz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sfsh + 50);

    auto g_0_xxz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sfsh + 51);

    auto g_0_xxz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sfsh + 52);

    auto g_0_xxz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sfsh + 53);

    auto g_0_xxz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sfsh + 54);

    auto g_0_xxz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sfsh + 55);

    auto g_0_xxz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 56);

    auto g_0_xxz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sfsh + 58);

    auto g_0_xxz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sfsh + 59);

    auto g_0_xxz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sfsh + 60);

    auto g_0_xxz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 61);

    auto g_0_xxz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 62);

    auto g_0_xyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_sfsh + 63);

    auto g_0_xyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_sfsh + 64);

    auto g_0_xyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_sfsh + 66);

    auto g_0_xyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_sfsh + 67);

    auto g_0_xyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_sfsh + 69);

    auto g_0_xyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_sfsh + 70);

    auto g_0_xyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_sfsh + 71);

    auto g_0_xyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_sfsh + 73);

    auto g_0_xyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_sfsh + 74);

    auto g_0_xyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_sfsh + 75);

    auto g_0_xyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_sfsh + 76);

    auto g_0_xyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_sfsh + 78);

    auto g_0_xyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_sfsh + 79);

    auto g_0_xyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_sfsh + 80);

    auto g_0_xyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_sfsh + 81);

    auto g_0_xyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 82);

    auto g_0_xyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 83);

    auto g_0_xzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sfsh + 105);

    auto g_0_xzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sfsh + 107);

    auto g_0_xzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sfsh + 109);

    auto g_0_xzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sfsh + 110);

    auto g_0_xzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sfsh + 112);

    auto g_0_xzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sfsh + 113);

    auto g_0_xzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sfsh + 114);

    auto g_0_xzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sfsh + 116);

    auto g_0_xzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sfsh + 117);

    auto g_0_xzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sfsh + 118);

    auto g_0_xzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 119);

    auto g_0_xzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sfsh + 120);

    auto g_0_xzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sfsh + 121);

    auto g_0_xzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sfsh + 122);

    auto g_0_xzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sfsh + 123);

    auto g_0_xzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 124);

    auto g_0_xzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 125);

    auto g_0_yyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_sfsh + 126);

    auto g_0_yyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_sfsh + 127);

    auto g_0_yyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_sfsh + 128);

    auto g_0_yyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_sfsh + 129);

    auto g_0_yyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_sfsh + 130);

    auto g_0_yyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_sfsh + 131);

    auto g_0_yyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_sfsh + 132);

    auto g_0_yyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_sfsh + 133);

    auto g_0_yyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_sfsh + 134);

    auto g_0_yyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_sfsh + 135);

    auto g_0_yyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_sfsh + 136);

    auto g_0_yyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_sfsh + 137);

    auto g_0_yyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_sfsh + 138);

    auto g_0_yyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_sfsh + 139);

    auto g_0_yyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 140);

    auto g_0_yyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_sfsh + 141);

    auto g_0_yyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_sfsh + 142);

    auto g_0_yyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_sfsh + 143);

    auto g_0_yyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_sfsh + 144);

    auto g_0_yyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 145);

    auto g_0_yyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 146);

    auto g_0_yyz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sfsh + 148);

    auto g_0_yyz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sfsh + 149);

    auto g_0_yyz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sfsh + 150);

    auto g_0_yyz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sfsh + 151);

    auto g_0_yyz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sfsh + 152);

    auto g_0_yyz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sfsh + 153);

    auto g_0_yyz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sfsh + 154);

    auto g_0_yyz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sfsh + 155);

    auto g_0_yyz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sfsh + 156);

    auto g_0_yyz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sfsh + 157);

    auto g_0_yyz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sfsh + 158);

    auto g_0_yyz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sfsh + 159);

    auto g_0_yyz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sfsh + 160);

    auto g_0_yyz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 161);

    auto g_0_yyz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sfsh + 162);

    auto g_0_yyz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sfsh + 163);

    auto g_0_yyz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sfsh + 164);

    auto g_0_yyz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sfsh + 165);

    auto g_0_yyz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 166);

    auto g_0_yyz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 167);

    auto g_0_yzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sfsh + 168);

    auto g_0_yzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sfsh + 169);

    auto g_0_yzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sfsh + 170);

    auto g_0_yzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sfsh + 171);

    auto g_0_yzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sfsh + 172);

    auto g_0_yzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sfsh + 173);

    auto g_0_yzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sfsh + 174);

    auto g_0_yzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sfsh + 175);

    auto g_0_yzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sfsh + 176);

    auto g_0_yzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sfsh + 177);

    auto g_0_yzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sfsh + 178);

    auto g_0_yzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sfsh + 179);

    auto g_0_yzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sfsh + 180);

    auto g_0_yzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sfsh + 181);

    auto g_0_yzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 182);

    auto g_0_yzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sfsh + 183);

    auto g_0_yzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sfsh + 184);

    auto g_0_yzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sfsh + 185);

    auto g_0_yzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sfsh + 186);

    auto g_0_yzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 187);

    auto g_0_yzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 188);

    auto g_0_zzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sfsh + 189);

    auto g_0_zzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sfsh + 190);

    auto g_0_zzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sfsh + 191);

    auto g_0_zzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sfsh + 192);

    auto g_0_zzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sfsh + 193);

    auto g_0_zzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sfsh + 194);

    auto g_0_zzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sfsh + 195);

    auto g_0_zzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sfsh + 196);

    auto g_0_zzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sfsh + 197);

    auto g_0_zzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sfsh + 198);

    auto g_0_zzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sfsh + 199);

    auto g_0_zzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sfsh + 200);

    auto g_0_zzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sfsh + 201);

    auto g_0_zzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sfsh + 202);

    auto g_0_zzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 203);

    auto g_0_zzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sfsh + 204);

    auto g_0_zzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sfsh + 205);

    auto g_0_zzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sfsh + 206);

    auto g_0_zzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sfsh + 207);

    auto g_0_zzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 208);

    auto g_0_zzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sfsh + 209);

    /// Set up 0-21 components of targeted buffer : SGSH

    auto g_0_xxxx_0_xxxxx_0 = pbuffer.data(idx_eri_0_sgsh);

    auto g_0_xxxx_0_xxxxy_0 = pbuffer.data(idx_eri_0_sgsh + 1);

    auto g_0_xxxx_0_xxxxz_0 = pbuffer.data(idx_eri_0_sgsh + 2);

    auto g_0_xxxx_0_xxxyy_0 = pbuffer.data(idx_eri_0_sgsh + 3);

    auto g_0_xxxx_0_xxxyz_0 = pbuffer.data(idx_eri_0_sgsh + 4);

    auto g_0_xxxx_0_xxxzz_0 = pbuffer.data(idx_eri_0_sgsh + 5);

    auto g_0_xxxx_0_xxyyy_0 = pbuffer.data(idx_eri_0_sgsh + 6);

    auto g_0_xxxx_0_xxyyz_0 = pbuffer.data(idx_eri_0_sgsh + 7);

    auto g_0_xxxx_0_xxyzz_0 = pbuffer.data(idx_eri_0_sgsh + 8);

    auto g_0_xxxx_0_xxzzz_0 = pbuffer.data(idx_eri_0_sgsh + 9);

    auto g_0_xxxx_0_xyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 10);

    auto g_0_xxxx_0_xyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 11);

    auto g_0_xxxx_0_xyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 12);

    auto g_0_xxxx_0_xyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 13);

    auto g_0_xxxx_0_xzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 14);

    auto g_0_xxxx_0_yyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 15);

    auto g_0_xxxx_0_yyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 16);

    auto g_0_xxxx_0_yyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 17);

    auto g_0_xxxx_0_yyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 18);

    auto g_0_xxxx_0_yzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 19);

    auto g_0_xxxx_0_zzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 20);

#pragma omp simd aligned(g_0_xx_0_xxxxx_0,       \
                             g_0_xx_0_xxxxx_1,   \
                             g_0_xx_0_xxxxy_0,   \
                             g_0_xx_0_xxxxy_1,   \
                             g_0_xx_0_xxxxz_0,   \
                             g_0_xx_0_xxxxz_1,   \
                             g_0_xx_0_xxxyy_0,   \
                             g_0_xx_0_xxxyy_1,   \
                             g_0_xx_0_xxxyz_0,   \
                             g_0_xx_0_xxxyz_1,   \
                             g_0_xx_0_xxxzz_0,   \
                             g_0_xx_0_xxxzz_1,   \
                             g_0_xx_0_xxyyy_0,   \
                             g_0_xx_0_xxyyy_1,   \
                             g_0_xx_0_xxyyz_0,   \
                             g_0_xx_0_xxyyz_1,   \
                             g_0_xx_0_xxyzz_0,   \
                             g_0_xx_0_xxyzz_1,   \
                             g_0_xx_0_xxzzz_0,   \
                             g_0_xx_0_xxzzz_1,   \
                             g_0_xx_0_xyyyy_0,   \
                             g_0_xx_0_xyyyy_1,   \
                             g_0_xx_0_xyyyz_0,   \
                             g_0_xx_0_xyyyz_1,   \
                             g_0_xx_0_xyyzz_0,   \
                             g_0_xx_0_xyyzz_1,   \
                             g_0_xx_0_xyzzz_0,   \
                             g_0_xx_0_xyzzz_1,   \
                             g_0_xx_0_xzzzz_0,   \
                             g_0_xx_0_xzzzz_1,   \
                             g_0_xx_0_yyyyy_0,   \
                             g_0_xx_0_yyyyy_1,   \
                             g_0_xx_0_yyyyz_0,   \
                             g_0_xx_0_yyyyz_1,   \
                             g_0_xx_0_yyyzz_0,   \
                             g_0_xx_0_yyyzz_1,   \
                             g_0_xx_0_yyzzz_0,   \
                             g_0_xx_0_yyzzz_1,   \
                             g_0_xx_0_yzzzz_0,   \
                             g_0_xx_0_yzzzz_1,   \
                             g_0_xx_0_zzzzz_0,   \
                             g_0_xx_0_zzzzz_1,   \
                             g_0_xxx_0_xxxx_1,   \
                             g_0_xxx_0_xxxxx_0,  \
                             g_0_xxx_0_xxxxx_1,  \
                             g_0_xxx_0_xxxxy_0,  \
                             g_0_xxx_0_xxxxy_1,  \
                             g_0_xxx_0_xxxxz_0,  \
                             g_0_xxx_0_xxxxz_1,  \
                             g_0_xxx_0_xxxy_1,   \
                             g_0_xxx_0_xxxyy_0,  \
                             g_0_xxx_0_xxxyy_1,  \
                             g_0_xxx_0_xxxyz_0,  \
                             g_0_xxx_0_xxxyz_1,  \
                             g_0_xxx_0_xxxz_1,   \
                             g_0_xxx_0_xxxzz_0,  \
                             g_0_xxx_0_xxxzz_1,  \
                             g_0_xxx_0_xxyy_1,   \
                             g_0_xxx_0_xxyyy_0,  \
                             g_0_xxx_0_xxyyy_1,  \
                             g_0_xxx_0_xxyyz_0,  \
                             g_0_xxx_0_xxyyz_1,  \
                             g_0_xxx_0_xxyz_1,   \
                             g_0_xxx_0_xxyzz_0,  \
                             g_0_xxx_0_xxyzz_1,  \
                             g_0_xxx_0_xxzz_1,   \
                             g_0_xxx_0_xxzzz_0,  \
                             g_0_xxx_0_xxzzz_1,  \
                             g_0_xxx_0_xyyy_1,   \
                             g_0_xxx_0_xyyyy_0,  \
                             g_0_xxx_0_xyyyy_1,  \
                             g_0_xxx_0_xyyyz_0,  \
                             g_0_xxx_0_xyyyz_1,  \
                             g_0_xxx_0_xyyz_1,   \
                             g_0_xxx_0_xyyzz_0,  \
                             g_0_xxx_0_xyyzz_1,  \
                             g_0_xxx_0_xyzz_1,   \
                             g_0_xxx_0_xyzzz_0,  \
                             g_0_xxx_0_xyzzz_1,  \
                             g_0_xxx_0_xzzz_1,   \
                             g_0_xxx_0_xzzzz_0,  \
                             g_0_xxx_0_xzzzz_1,  \
                             g_0_xxx_0_yyyy_1,   \
                             g_0_xxx_0_yyyyy_0,  \
                             g_0_xxx_0_yyyyy_1,  \
                             g_0_xxx_0_yyyyz_0,  \
                             g_0_xxx_0_yyyyz_1,  \
                             g_0_xxx_0_yyyz_1,   \
                             g_0_xxx_0_yyyzz_0,  \
                             g_0_xxx_0_yyyzz_1,  \
                             g_0_xxx_0_yyzz_1,   \
                             g_0_xxx_0_yyzzz_0,  \
                             g_0_xxx_0_yyzzz_1,  \
                             g_0_xxx_0_yzzz_1,   \
                             g_0_xxx_0_yzzzz_0,  \
                             g_0_xxx_0_yzzzz_1,  \
                             g_0_xxx_0_zzzz_1,   \
                             g_0_xxx_0_zzzzz_0,  \
                             g_0_xxx_0_zzzzz_1,  \
                             g_0_xxxx_0_xxxxx_0, \
                             g_0_xxxx_0_xxxxy_0, \
                             g_0_xxxx_0_xxxxz_0, \
                             g_0_xxxx_0_xxxyy_0, \
                             g_0_xxxx_0_xxxyz_0, \
                             g_0_xxxx_0_xxxzz_0, \
                             g_0_xxxx_0_xxyyy_0, \
                             g_0_xxxx_0_xxyyz_0, \
                             g_0_xxxx_0_xxyzz_0, \
                             g_0_xxxx_0_xxzzz_0, \
                             g_0_xxxx_0_xyyyy_0, \
                             g_0_xxxx_0_xyyyz_0, \
                             g_0_xxxx_0_xyyzz_0, \
                             g_0_xxxx_0_xyzzz_0, \
                             g_0_xxxx_0_xzzzz_0, \
                             g_0_xxxx_0_yyyyy_0, \
                             g_0_xxxx_0_yyyyz_0, \
                             g_0_xxxx_0_yyyzz_0, \
                             g_0_xxxx_0_yyzzz_0, \
                             g_0_xxxx_0_yzzzz_0, \
                             g_0_xxxx_0_zzzzz_0, \
                             wp_x,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxx_0_xxxxx_0[i] = 3.0 * g_0_xx_0_xxxxx_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxxx_1[i] * fti_ab_0 + 5.0 * g_0_xxx_0_xxxx_1[i] * fi_abcd_0 +
                                g_0_xxx_0_xxxxx_0[i] * pb_x + g_0_xxx_0_xxxxx_1[i] * wp_x[i];

        g_0_xxxx_0_xxxxy_0[i] = 3.0 * g_0_xx_0_xxxxy_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxxy_1[i] * fti_ab_0 + 4.0 * g_0_xxx_0_xxxy_1[i] * fi_abcd_0 +
                                g_0_xxx_0_xxxxy_0[i] * pb_x + g_0_xxx_0_xxxxy_1[i] * wp_x[i];

        g_0_xxxx_0_xxxxz_0[i] = 3.0 * g_0_xx_0_xxxxz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxxz_1[i] * fti_ab_0 + 4.0 * g_0_xxx_0_xxxz_1[i] * fi_abcd_0 +
                                g_0_xxx_0_xxxxz_0[i] * pb_x + g_0_xxx_0_xxxxz_1[i] * wp_x[i];

        g_0_xxxx_0_xxxyy_0[i] = 3.0 * g_0_xx_0_xxxyy_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxyy_1[i] * fti_ab_0 + 3.0 * g_0_xxx_0_xxyy_1[i] * fi_abcd_0 +
                                g_0_xxx_0_xxxyy_0[i] * pb_x + g_0_xxx_0_xxxyy_1[i] * wp_x[i];

        g_0_xxxx_0_xxxyz_0[i] = 3.0 * g_0_xx_0_xxxyz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xxx_0_xxyz_1[i] * fi_abcd_0 +
                                g_0_xxx_0_xxxyz_0[i] * pb_x + g_0_xxx_0_xxxyz_1[i] * wp_x[i];

        g_0_xxxx_0_xxxzz_0[i] = 3.0 * g_0_xx_0_xxxzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxzz_1[i] * fti_ab_0 + 3.0 * g_0_xxx_0_xxzz_1[i] * fi_abcd_0 +
                                g_0_xxx_0_xxxzz_0[i] * pb_x + g_0_xxx_0_xxxzz_1[i] * wp_x[i];

        g_0_xxxx_0_xxyyy_0[i] = 3.0 * g_0_xx_0_xxyyy_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxyyy_1[i] * fti_ab_0 + 2.0 * g_0_xxx_0_xyyy_1[i] * fi_abcd_0 +
                                g_0_xxx_0_xxyyy_0[i] * pb_x + g_0_xxx_0_xxyyy_1[i] * wp_x[i];

        g_0_xxxx_0_xxyyz_0[i] = 3.0 * g_0_xx_0_xxyyz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxx_0_xyyz_1[i] * fi_abcd_0 +
                                g_0_xxx_0_xxyyz_0[i] * pb_x + g_0_xxx_0_xxyyz_1[i] * wp_x[i];

        g_0_xxxx_0_xxyzz_0[i] = 3.0 * g_0_xx_0_xxyzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxx_0_xyzz_1[i] * fi_abcd_0 +
                                g_0_xxx_0_xxyzz_0[i] * pb_x + g_0_xxx_0_xxyzz_1[i] * wp_x[i];

        g_0_xxxx_0_xxzzz_0[i] = 3.0 * g_0_xx_0_xxzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxx_0_xzzz_1[i] * fi_abcd_0 +
                                g_0_xxx_0_xxzzz_0[i] * pb_x + g_0_xxx_0_xxzzz_1[i] * wp_x[i];

        g_0_xxxx_0_xyyyy_0[i] = 3.0 * g_0_xx_0_xyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xyyyy_1[i] * fti_ab_0 + g_0_xxx_0_yyyy_1[i] * fi_abcd_0 +
                                g_0_xxx_0_xyyyy_0[i] * pb_x + g_0_xxx_0_xyyyy_1[i] * wp_x[i];

        g_0_xxxx_0_xyyyz_0[i] = 3.0 * g_0_xx_0_xyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xyyyz_1[i] * fti_ab_0 + g_0_xxx_0_yyyz_1[i] * fi_abcd_0 +
                                g_0_xxx_0_xyyyz_0[i] * pb_x + g_0_xxx_0_xyyyz_1[i] * wp_x[i];

        g_0_xxxx_0_xyyzz_0[i] = 3.0 * g_0_xx_0_xyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xyyzz_1[i] * fti_ab_0 + g_0_xxx_0_yyzz_1[i] * fi_abcd_0 +
                                g_0_xxx_0_xyyzz_0[i] * pb_x + g_0_xxx_0_xyyzz_1[i] * wp_x[i];

        g_0_xxxx_0_xyzzz_0[i] = 3.0 * g_0_xx_0_xyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xyzzz_1[i] * fti_ab_0 + g_0_xxx_0_yzzz_1[i] * fi_abcd_0 +
                                g_0_xxx_0_xyzzz_0[i] * pb_x + g_0_xxx_0_xyzzz_1[i] * wp_x[i];

        g_0_xxxx_0_xzzzz_0[i] = 3.0 * g_0_xx_0_xzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xzzzz_1[i] * fti_ab_0 + g_0_xxx_0_zzzz_1[i] * fi_abcd_0 +
                                g_0_xxx_0_xzzzz_0[i] * pb_x + g_0_xxx_0_xzzzz_1[i] * wp_x[i];

        g_0_xxxx_0_yyyyy_0[i] =
            3.0 * g_0_xx_0_yyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_yyyyy_1[i] * fti_ab_0 + g_0_xxx_0_yyyyy_0[i] * pb_x + g_0_xxx_0_yyyyy_1[i] * wp_x[i];

        g_0_xxxx_0_yyyyz_0[i] =
            3.0 * g_0_xx_0_yyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_yyyyz_1[i] * fti_ab_0 + g_0_xxx_0_yyyyz_0[i] * pb_x + g_0_xxx_0_yyyyz_1[i] * wp_x[i];

        g_0_xxxx_0_yyyzz_0[i] =
            3.0 * g_0_xx_0_yyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_yyyzz_1[i] * fti_ab_0 + g_0_xxx_0_yyyzz_0[i] * pb_x + g_0_xxx_0_yyyzz_1[i] * wp_x[i];

        g_0_xxxx_0_yyzzz_0[i] =
            3.0 * g_0_xx_0_yyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_yyzzz_1[i] * fti_ab_0 + g_0_xxx_0_yyzzz_0[i] * pb_x + g_0_xxx_0_yyzzz_1[i] * wp_x[i];

        g_0_xxxx_0_yzzzz_0[i] =
            3.0 * g_0_xx_0_yzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_yzzzz_1[i] * fti_ab_0 + g_0_xxx_0_yzzzz_0[i] * pb_x + g_0_xxx_0_yzzzz_1[i] * wp_x[i];

        g_0_xxxx_0_zzzzz_0[i] =
            3.0 * g_0_xx_0_zzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_zzzzz_1[i] * fti_ab_0 + g_0_xxx_0_zzzzz_0[i] * pb_x + g_0_xxx_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 21-42 components of targeted buffer : SGSH

    auto g_0_xxxy_0_xxxxx_0 = pbuffer.data(idx_eri_0_sgsh + 21);

    auto g_0_xxxy_0_xxxxy_0 = pbuffer.data(idx_eri_0_sgsh + 22);

    auto g_0_xxxy_0_xxxxz_0 = pbuffer.data(idx_eri_0_sgsh + 23);

    auto g_0_xxxy_0_xxxyy_0 = pbuffer.data(idx_eri_0_sgsh + 24);

    auto g_0_xxxy_0_xxxyz_0 = pbuffer.data(idx_eri_0_sgsh + 25);

    auto g_0_xxxy_0_xxxzz_0 = pbuffer.data(idx_eri_0_sgsh + 26);

    auto g_0_xxxy_0_xxyyy_0 = pbuffer.data(idx_eri_0_sgsh + 27);

    auto g_0_xxxy_0_xxyyz_0 = pbuffer.data(idx_eri_0_sgsh + 28);

    auto g_0_xxxy_0_xxyzz_0 = pbuffer.data(idx_eri_0_sgsh + 29);

    auto g_0_xxxy_0_xxzzz_0 = pbuffer.data(idx_eri_0_sgsh + 30);

    auto g_0_xxxy_0_xyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 31);

    auto g_0_xxxy_0_xyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 32);

    auto g_0_xxxy_0_xyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 33);

    auto g_0_xxxy_0_xyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 34);

    auto g_0_xxxy_0_xzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 35);

    auto g_0_xxxy_0_yyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 36);

    auto g_0_xxxy_0_yyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 37);

    auto g_0_xxxy_0_yyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 38);

    auto g_0_xxxy_0_yyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 39);

    auto g_0_xxxy_0_yzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 40);

    auto g_0_xxxy_0_zzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 41);

#pragma omp simd aligned(g_0_xxx_0_xxxx_1,       \
                             g_0_xxx_0_xxxxx_0,  \
                             g_0_xxx_0_xxxxx_1,  \
                             g_0_xxx_0_xxxxy_0,  \
                             g_0_xxx_0_xxxxy_1,  \
                             g_0_xxx_0_xxxxz_0,  \
                             g_0_xxx_0_xxxxz_1,  \
                             g_0_xxx_0_xxxy_1,   \
                             g_0_xxx_0_xxxyy_0,  \
                             g_0_xxx_0_xxxyy_1,  \
                             g_0_xxx_0_xxxyz_0,  \
                             g_0_xxx_0_xxxyz_1,  \
                             g_0_xxx_0_xxxz_1,   \
                             g_0_xxx_0_xxxzz_0,  \
                             g_0_xxx_0_xxxzz_1,  \
                             g_0_xxx_0_xxyy_1,   \
                             g_0_xxx_0_xxyyy_0,  \
                             g_0_xxx_0_xxyyy_1,  \
                             g_0_xxx_0_xxyyz_0,  \
                             g_0_xxx_0_xxyyz_1,  \
                             g_0_xxx_0_xxyz_1,   \
                             g_0_xxx_0_xxyzz_0,  \
                             g_0_xxx_0_xxyzz_1,  \
                             g_0_xxx_0_xxzz_1,   \
                             g_0_xxx_0_xxzzz_0,  \
                             g_0_xxx_0_xxzzz_1,  \
                             g_0_xxx_0_xyyy_1,   \
                             g_0_xxx_0_xyyyy_0,  \
                             g_0_xxx_0_xyyyy_1,  \
                             g_0_xxx_0_xyyyz_0,  \
                             g_0_xxx_0_xyyyz_1,  \
                             g_0_xxx_0_xyyz_1,   \
                             g_0_xxx_0_xyyzz_0,  \
                             g_0_xxx_0_xyyzz_1,  \
                             g_0_xxx_0_xyzz_1,   \
                             g_0_xxx_0_xyzzz_0,  \
                             g_0_xxx_0_xyzzz_1,  \
                             g_0_xxx_0_xzzz_1,   \
                             g_0_xxx_0_xzzzz_0,  \
                             g_0_xxx_0_xzzzz_1,  \
                             g_0_xxx_0_yyyy_1,   \
                             g_0_xxx_0_yyyyy_0,  \
                             g_0_xxx_0_yyyyy_1,  \
                             g_0_xxx_0_yyyyz_0,  \
                             g_0_xxx_0_yyyyz_1,  \
                             g_0_xxx_0_yyyz_1,   \
                             g_0_xxx_0_yyyzz_0,  \
                             g_0_xxx_0_yyyzz_1,  \
                             g_0_xxx_0_yyzz_1,   \
                             g_0_xxx_0_yyzzz_0,  \
                             g_0_xxx_0_yyzzz_1,  \
                             g_0_xxx_0_yzzz_1,   \
                             g_0_xxx_0_yzzzz_0,  \
                             g_0_xxx_0_yzzzz_1,  \
                             g_0_xxx_0_zzzz_1,   \
                             g_0_xxx_0_zzzzz_0,  \
                             g_0_xxx_0_zzzzz_1,  \
                             g_0_xxxy_0_xxxxx_0, \
                             g_0_xxxy_0_xxxxy_0, \
                             g_0_xxxy_0_xxxxz_0, \
                             g_0_xxxy_0_xxxyy_0, \
                             g_0_xxxy_0_xxxyz_0, \
                             g_0_xxxy_0_xxxzz_0, \
                             g_0_xxxy_0_xxyyy_0, \
                             g_0_xxxy_0_xxyyz_0, \
                             g_0_xxxy_0_xxyzz_0, \
                             g_0_xxxy_0_xxzzz_0, \
                             g_0_xxxy_0_xyyyy_0, \
                             g_0_xxxy_0_xyyyz_0, \
                             g_0_xxxy_0_xyyzz_0, \
                             g_0_xxxy_0_xyzzz_0, \
                             g_0_xxxy_0_xzzzz_0, \
                             g_0_xxxy_0_yyyyy_0, \
                             g_0_xxxy_0_yyyyz_0, \
                             g_0_xxxy_0_yyyzz_0, \
                             g_0_xxxy_0_yyzzz_0, \
                             g_0_xxxy_0_yzzzz_0, \
                             g_0_xxxy_0_zzzzz_0, \
                             wp_y,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxy_0_xxxxx_0[i] = g_0_xxx_0_xxxxx_0[i] * pb_y + g_0_xxx_0_xxxxx_1[i] * wp_y[i];

        g_0_xxxy_0_xxxxy_0[i] = g_0_xxx_0_xxxx_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxy_0[i] * pb_y + g_0_xxx_0_xxxxy_1[i] * wp_y[i];

        g_0_xxxy_0_xxxxz_0[i] = g_0_xxx_0_xxxxz_0[i] * pb_y + g_0_xxx_0_xxxxz_1[i] * wp_y[i];

        g_0_xxxy_0_xxxyy_0[i] = 2.0 * g_0_xxx_0_xxxy_1[i] * fi_abcd_0 + g_0_xxx_0_xxxyy_0[i] * pb_y + g_0_xxx_0_xxxyy_1[i] * wp_y[i];

        g_0_xxxy_0_xxxyz_0[i] = g_0_xxx_0_xxxz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxyz_0[i] * pb_y + g_0_xxx_0_xxxyz_1[i] * wp_y[i];

        g_0_xxxy_0_xxxzz_0[i] = g_0_xxx_0_xxxzz_0[i] * pb_y + g_0_xxx_0_xxxzz_1[i] * wp_y[i];

        g_0_xxxy_0_xxyyy_0[i] = 3.0 * g_0_xxx_0_xxyy_1[i] * fi_abcd_0 + g_0_xxx_0_xxyyy_0[i] * pb_y + g_0_xxx_0_xxyyy_1[i] * wp_y[i];

        g_0_xxxy_0_xxyyz_0[i] = 2.0 * g_0_xxx_0_xxyz_1[i] * fi_abcd_0 + g_0_xxx_0_xxyyz_0[i] * pb_y + g_0_xxx_0_xxyyz_1[i] * wp_y[i];

        g_0_xxxy_0_xxyzz_0[i] = g_0_xxx_0_xxzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxyzz_0[i] * pb_y + g_0_xxx_0_xxyzz_1[i] * wp_y[i];

        g_0_xxxy_0_xxzzz_0[i] = g_0_xxx_0_xxzzz_0[i] * pb_y + g_0_xxx_0_xxzzz_1[i] * wp_y[i];

        g_0_xxxy_0_xyyyy_0[i] = 4.0 * g_0_xxx_0_xyyy_1[i] * fi_abcd_0 + g_0_xxx_0_xyyyy_0[i] * pb_y + g_0_xxx_0_xyyyy_1[i] * wp_y[i];

        g_0_xxxy_0_xyyyz_0[i] = 3.0 * g_0_xxx_0_xyyz_1[i] * fi_abcd_0 + g_0_xxx_0_xyyyz_0[i] * pb_y + g_0_xxx_0_xyyyz_1[i] * wp_y[i];

        g_0_xxxy_0_xyyzz_0[i] = 2.0 * g_0_xxx_0_xyzz_1[i] * fi_abcd_0 + g_0_xxx_0_xyyzz_0[i] * pb_y + g_0_xxx_0_xyyzz_1[i] * wp_y[i];

        g_0_xxxy_0_xyzzz_0[i] = g_0_xxx_0_xzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xyzzz_0[i] * pb_y + g_0_xxx_0_xyzzz_1[i] * wp_y[i];

        g_0_xxxy_0_xzzzz_0[i] = g_0_xxx_0_xzzzz_0[i] * pb_y + g_0_xxx_0_xzzzz_1[i] * wp_y[i];

        g_0_xxxy_0_yyyyy_0[i] = 5.0 * g_0_xxx_0_yyyy_1[i] * fi_abcd_0 + g_0_xxx_0_yyyyy_0[i] * pb_y + g_0_xxx_0_yyyyy_1[i] * wp_y[i];

        g_0_xxxy_0_yyyyz_0[i] = 4.0 * g_0_xxx_0_yyyz_1[i] * fi_abcd_0 + g_0_xxx_0_yyyyz_0[i] * pb_y + g_0_xxx_0_yyyyz_1[i] * wp_y[i];

        g_0_xxxy_0_yyyzz_0[i] = 3.0 * g_0_xxx_0_yyzz_1[i] * fi_abcd_0 + g_0_xxx_0_yyyzz_0[i] * pb_y + g_0_xxx_0_yyyzz_1[i] * wp_y[i];

        g_0_xxxy_0_yyzzz_0[i] = 2.0 * g_0_xxx_0_yzzz_1[i] * fi_abcd_0 + g_0_xxx_0_yyzzz_0[i] * pb_y + g_0_xxx_0_yyzzz_1[i] * wp_y[i];

        g_0_xxxy_0_yzzzz_0[i] = g_0_xxx_0_zzzz_1[i] * fi_abcd_0 + g_0_xxx_0_yzzzz_0[i] * pb_y + g_0_xxx_0_yzzzz_1[i] * wp_y[i];

        g_0_xxxy_0_zzzzz_0[i] = g_0_xxx_0_zzzzz_0[i] * pb_y + g_0_xxx_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 42-63 components of targeted buffer : SGSH

    auto g_0_xxxz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sgsh + 42);

    auto g_0_xxxz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sgsh + 43);

    auto g_0_xxxz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sgsh + 44);

    auto g_0_xxxz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sgsh + 45);

    auto g_0_xxxz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sgsh + 46);

    auto g_0_xxxz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sgsh + 47);

    auto g_0_xxxz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sgsh + 48);

    auto g_0_xxxz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sgsh + 49);

    auto g_0_xxxz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sgsh + 50);

    auto g_0_xxxz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sgsh + 51);

    auto g_0_xxxz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 52);

    auto g_0_xxxz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 53);

    auto g_0_xxxz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 54);

    auto g_0_xxxz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 55);

    auto g_0_xxxz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 56);

    auto g_0_xxxz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 57);

    auto g_0_xxxz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 58);

    auto g_0_xxxz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 59);

    auto g_0_xxxz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 60);

    auto g_0_xxxz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 61);

    auto g_0_xxxz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 62);

#pragma omp simd aligned(g_0_xxx_0_xxxx_1,       \
                             g_0_xxx_0_xxxxx_0,  \
                             g_0_xxx_0_xxxxx_1,  \
                             g_0_xxx_0_xxxxy_0,  \
                             g_0_xxx_0_xxxxy_1,  \
                             g_0_xxx_0_xxxxz_0,  \
                             g_0_xxx_0_xxxxz_1,  \
                             g_0_xxx_0_xxxy_1,   \
                             g_0_xxx_0_xxxyy_0,  \
                             g_0_xxx_0_xxxyy_1,  \
                             g_0_xxx_0_xxxyz_0,  \
                             g_0_xxx_0_xxxyz_1,  \
                             g_0_xxx_0_xxxz_1,   \
                             g_0_xxx_0_xxxzz_0,  \
                             g_0_xxx_0_xxxzz_1,  \
                             g_0_xxx_0_xxyy_1,   \
                             g_0_xxx_0_xxyyy_0,  \
                             g_0_xxx_0_xxyyy_1,  \
                             g_0_xxx_0_xxyyz_0,  \
                             g_0_xxx_0_xxyyz_1,  \
                             g_0_xxx_0_xxyz_1,   \
                             g_0_xxx_0_xxyzz_0,  \
                             g_0_xxx_0_xxyzz_1,  \
                             g_0_xxx_0_xxzz_1,   \
                             g_0_xxx_0_xxzzz_0,  \
                             g_0_xxx_0_xxzzz_1,  \
                             g_0_xxx_0_xyyy_1,   \
                             g_0_xxx_0_xyyyy_0,  \
                             g_0_xxx_0_xyyyy_1,  \
                             g_0_xxx_0_xyyyz_0,  \
                             g_0_xxx_0_xyyyz_1,  \
                             g_0_xxx_0_xyyz_1,   \
                             g_0_xxx_0_xyyzz_0,  \
                             g_0_xxx_0_xyyzz_1,  \
                             g_0_xxx_0_xyzz_1,   \
                             g_0_xxx_0_xyzzz_0,  \
                             g_0_xxx_0_xyzzz_1,  \
                             g_0_xxx_0_xzzz_1,   \
                             g_0_xxx_0_xzzzz_0,  \
                             g_0_xxx_0_xzzzz_1,  \
                             g_0_xxx_0_yyyy_1,   \
                             g_0_xxx_0_yyyyy_0,  \
                             g_0_xxx_0_yyyyy_1,  \
                             g_0_xxx_0_yyyyz_0,  \
                             g_0_xxx_0_yyyyz_1,  \
                             g_0_xxx_0_yyyz_1,   \
                             g_0_xxx_0_yyyzz_0,  \
                             g_0_xxx_0_yyyzz_1,  \
                             g_0_xxx_0_yyzz_1,   \
                             g_0_xxx_0_yyzzz_0,  \
                             g_0_xxx_0_yyzzz_1,  \
                             g_0_xxx_0_yzzz_1,   \
                             g_0_xxx_0_yzzzz_0,  \
                             g_0_xxx_0_yzzzz_1,  \
                             g_0_xxx_0_zzzz_1,   \
                             g_0_xxx_0_zzzzz_0,  \
                             g_0_xxx_0_zzzzz_1,  \
                             g_0_xxxz_0_xxxxx_0, \
                             g_0_xxxz_0_xxxxy_0, \
                             g_0_xxxz_0_xxxxz_0, \
                             g_0_xxxz_0_xxxyy_0, \
                             g_0_xxxz_0_xxxyz_0, \
                             g_0_xxxz_0_xxxzz_0, \
                             g_0_xxxz_0_xxyyy_0, \
                             g_0_xxxz_0_xxyyz_0, \
                             g_0_xxxz_0_xxyzz_0, \
                             g_0_xxxz_0_xxzzz_0, \
                             g_0_xxxz_0_xyyyy_0, \
                             g_0_xxxz_0_xyyyz_0, \
                             g_0_xxxz_0_xyyzz_0, \
                             g_0_xxxz_0_xyzzz_0, \
                             g_0_xxxz_0_xzzzz_0, \
                             g_0_xxxz_0_yyyyy_0, \
                             g_0_xxxz_0_yyyyz_0, \
                             g_0_xxxz_0_yyyzz_0, \
                             g_0_xxxz_0_yyzzz_0, \
                             g_0_xxxz_0_yzzzz_0, \
                             g_0_xxxz_0_zzzzz_0, \
                             wp_z,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxz_0_xxxxx_0[i] = g_0_xxx_0_xxxxx_0[i] * pb_z + g_0_xxx_0_xxxxx_1[i] * wp_z[i];

        g_0_xxxz_0_xxxxy_0[i] = g_0_xxx_0_xxxxy_0[i] * pb_z + g_0_xxx_0_xxxxy_1[i] * wp_z[i];

        g_0_xxxz_0_xxxxz_0[i] = g_0_xxx_0_xxxx_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxz_0[i] * pb_z + g_0_xxx_0_xxxxz_1[i] * wp_z[i];

        g_0_xxxz_0_xxxyy_0[i] = g_0_xxx_0_xxxyy_0[i] * pb_z + g_0_xxx_0_xxxyy_1[i] * wp_z[i];

        g_0_xxxz_0_xxxyz_0[i] = g_0_xxx_0_xxxy_1[i] * fi_abcd_0 + g_0_xxx_0_xxxyz_0[i] * pb_z + g_0_xxx_0_xxxyz_1[i] * wp_z[i];

        g_0_xxxz_0_xxxzz_0[i] = 2.0 * g_0_xxx_0_xxxz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxzz_0[i] * pb_z + g_0_xxx_0_xxxzz_1[i] * wp_z[i];

        g_0_xxxz_0_xxyyy_0[i] = g_0_xxx_0_xxyyy_0[i] * pb_z + g_0_xxx_0_xxyyy_1[i] * wp_z[i];

        g_0_xxxz_0_xxyyz_0[i] = g_0_xxx_0_xxyy_1[i] * fi_abcd_0 + g_0_xxx_0_xxyyz_0[i] * pb_z + g_0_xxx_0_xxyyz_1[i] * wp_z[i];

        g_0_xxxz_0_xxyzz_0[i] = 2.0 * g_0_xxx_0_xxyz_1[i] * fi_abcd_0 + g_0_xxx_0_xxyzz_0[i] * pb_z + g_0_xxx_0_xxyzz_1[i] * wp_z[i];

        g_0_xxxz_0_xxzzz_0[i] = 3.0 * g_0_xxx_0_xxzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxzzz_0[i] * pb_z + g_0_xxx_0_xxzzz_1[i] * wp_z[i];

        g_0_xxxz_0_xyyyy_0[i] = g_0_xxx_0_xyyyy_0[i] * pb_z + g_0_xxx_0_xyyyy_1[i] * wp_z[i];

        g_0_xxxz_0_xyyyz_0[i] = g_0_xxx_0_xyyy_1[i] * fi_abcd_0 + g_0_xxx_0_xyyyz_0[i] * pb_z + g_0_xxx_0_xyyyz_1[i] * wp_z[i];

        g_0_xxxz_0_xyyzz_0[i] = 2.0 * g_0_xxx_0_xyyz_1[i] * fi_abcd_0 + g_0_xxx_0_xyyzz_0[i] * pb_z + g_0_xxx_0_xyyzz_1[i] * wp_z[i];

        g_0_xxxz_0_xyzzz_0[i] = 3.0 * g_0_xxx_0_xyzz_1[i] * fi_abcd_0 + g_0_xxx_0_xyzzz_0[i] * pb_z + g_0_xxx_0_xyzzz_1[i] * wp_z[i];

        g_0_xxxz_0_xzzzz_0[i] = 4.0 * g_0_xxx_0_xzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xzzzz_0[i] * pb_z + g_0_xxx_0_xzzzz_1[i] * wp_z[i];

        g_0_xxxz_0_yyyyy_0[i] = g_0_xxx_0_yyyyy_0[i] * pb_z + g_0_xxx_0_yyyyy_1[i] * wp_z[i];

        g_0_xxxz_0_yyyyz_0[i] = g_0_xxx_0_yyyy_1[i] * fi_abcd_0 + g_0_xxx_0_yyyyz_0[i] * pb_z + g_0_xxx_0_yyyyz_1[i] * wp_z[i];

        g_0_xxxz_0_yyyzz_0[i] = 2.0 * g_0_xxx_0_yyyz_1[i] * fi_abcd_0 + g_0_xxx_0_yyyzz_0[i] * pb_z + g_0_xxx_0_yyyzz_1[i] * wp_z[i];

        g_0_xxxz_0_yyzzz_0[i] = 3.0 * g_0_xxx_0_yyzz_1[i] * fi_abcd_0 + g_0_xxx_0_yyzzz_0[i] * pb_z + g_0_xxx_0_yyzzz_1[i] * wp_z[i];

        g_0_xxxz_0_yzzzz_0[i] = 4.0 * g_0_xxx_0_yzzz_1[i] * fi_abcd_0 + g_0_xxx_0_yzzzz_0[i] * pb_z + g_0_xxx_0_yzzzz_1[i] * wp_z[i];

        g_0_xxxz_0_zzzzz_0[i] = 5.0 * g_0_xxx_0_zzzz_1[i] * fi_abcd_0 + g_0_xxx_0_zzzzz_0[i] * pb_z + g_0_xxx_0_zzzzz_1[i] * wp_z[i];
    }

    /// Set up 63-84 components of targeted buffer : SGSH

    auto g_0_xxyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_sgsh + 63);

    auto g_0_xxyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_sgsh + 64);

    auto g_0_xxyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_sgsh + 65);

    auto g_0_xxyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_sgsh + 66);

    auto g_0_xxyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_sgsh + 67);

    auto g_0_xxyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_sgsh + 68);

    auto g_0_xxyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_sgsh + 69);

    auto g_0_xxyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_sgsh + 70);

    auto g_0_xxyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_sgsh + 71);

    auto g_0_xxyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_sgsh + 72);

    auto g_0_xxyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 73);

    auto g_0_xxyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 74);

    auto g_0_xxyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 75);

    auto g_0_xxyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 76);

    auto g_0_xxyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 77);

    auto g_0_xxyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 78);

    auto g_0_xxyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 79);

    auto g_0_xxyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 80);

    auto g_0_xxyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 81);

    auto g_0_xxyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 82);

    auto g_0_xxyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 83);

#pragma omp simd aligned(g_0_xx_0_xxxxx_0,       \
                             g_0_xx_0_xxxxx_1,   \
                             g_0_xx_0_xxxxz_0,   \
                             g_0_xx_0_xxxxz_1,   \
                             g_0_xx_0_xxxzz_0,   \
                             g_0_xx_0_xxxzz_1,   \
                             g_0_xx_0_xxzzz_0,   \
                             g_0_xx_0_xxzzz_1,   \
                             g_0_xx_0_xzzzz_0,   \
                             g_0_xx_0_xzzzz_1,   \
                             g_0_xxy_0_xxxxx_0,  \
                             g_0_xxy_0_xxxxx_1,  \
                             g_0_xxy_0_xxxxz_0,  \
                             g_0_xxy_0_xxxxz_1,  \
                             g_0_xxy_0_xxxzz_0,  \
                             g_0_xxy_0_xxxzz_1,  \
                             g_0_xxy_0_xxzzz_0,  \
                             g_0_xxy_0_xxzzz_1,  \
                             g_0_xxy_0_xzzzz_0,  \
                             g_0_xxy_0_xzzzz_1,  \
                             g_0_xxyy_0_xxxxx_0, \
                             g_0_xxyy_0_xxxxy_0, \
                             g_0_xxyy_0_xxxxz_0, \
                             g_0_xxyy_0_xxxyy_0, \
                             g_0_xxyy_0_xxxyz_0, \
                             g_0_xxyy_0_xxxzz_0, \
                             g_0_xxyy_0_xxyyy_0, \
                             g_0_xxyy_0_xxyyz_0, \
                             g_0_xxyy_0_xxyzz_0, \
                             g_0_xxyy_0_xxzzz_0, \
                             g_0_xxyy_0_xyyyy_0, \
                             g_0_xxyy_0_xyyyz_0, \
                             g_0_xxyy_0_xyyzz_0, \
                             g_0_xxyy_0_xyzzz_0, \
                             g_0_xxyy_0_xzzzz_0, \
                             g_0_xxyy_0_yyyyy_0, \
                             g_0_xxyy_0_yyyyz_0, \
                             g_0_xxyy_0_yyyzz_0, \
                             g_0_xxyy_0_yyzzz_0, \
                             g_0_xxyy_0_yzzzz_0, \
                             g_0_xxyy_0_zzzzz_0, \
                             g_0_xyy_0_xxxxy_0,  \
                             g_0_xyy_0_xxxxy_1,  \
                             g_0_xyy_0_xxxy_1,   \
                             g_0_xyy_0_xxxyy_0,  \
                             g_0_xyy_0_xxxyy_1,  \
                             g_0_xyy_0_xxxyz_0,  \
                             g_0_xyy_0_xxxyz_1,  \
                             g_0_xyy_0_xxyy_1,   \
                             g_0_xyy_0_xxyyy_0,  \
                             g_0_xyy_0_xxyyy_1,  \
                             g_0_xyy_0_xxyyz_0,  \
                             g_0_xyy_0_xxyyz_1,  \
                             g_0_xyy_0_xxyz_1,   \
                             g_0_xyy_0_xxyzz_0,  \
                             g_0_xyy_0_xxyzz_1,  \
                             g_0_xyy_0_xyyy_1,   \
                             g_0_xyy_0_xyyyy_0,  \
                             g_0_xyy_0_xyyyy_1,  \
                             g_0_xyy_0_xyyyz_0,  \
                             g_0_xyy_0_xyyyz_1,  \
                             g_0_xyy_0_xyyz_1,   \
                             g_0_xyy_0_xyyzz_0,  \
                             g_0_xyy_0_xyyzz_1,  \
                             g_0_xyy_0_xyzz_1,   \
                             g_0_xyy_0_xyzzz_0,  \
                             g_0_xyy_0_xyzzz_1,  \
                             g_0_xyy_0_yyyy_1,   \
                             g_0_xyy_0_yyyyy_0,  \
                             g_0_xyy_0_yyyyy_1,  \
                             g_0_xyy_0_yyyyz_0,  \
                             g_0_xyy_0_yyyyz_1,  \
                             g_0_xyy_0_yyyz_1,   \
                             g_0_xyy_0_yyyzz_0,  \
                             g_0_xyy_0_yyyzz_1,  \
                             g_0_xyy_0_yyzz_1,   \
                             g_0_xyy_0_yyzzz_0,  \
                             g_0_xyy_0_yyzzz_1,  \
                             g_0_xyy_0_yzzz_1,   \
                             g_0_xyy_0_yzzzz_0,  \
                             g_0_xyy_0_yzzzz_1,  \
                             g_0_xyy_0_zzzzz_0,  \
                             g_0_xyy_0_zzzzz_1,  \
                             g_0_yy_0_xxxxy_0,   \
                             g_0_yy_0_xxxxy_1,   \
                             g_0_yy_0_xxxyy_0,   \
                             g_0_yy_0_xxxyy_1,   \
                             g_0_yy_0_xxxyz_0,   \
                             g_0_yy_0_xxxyz_1,   \
                             g_0_yy_0_xxyyy_0,   \
                             g_0_yy_0_xxyyy_1,   \
                             g_0_yy_0_xxyyz_0,   \
                             g_0_yy_0_xxyyz_1,   \
                             g_0_yy_0_xxyzz_0,   \
                             g_0_yy_0_xxyzz_1,   \
                             g_0_yy_0_xyyyy_0,   \
                             g_0_yy_0_xyyyy_1,   \
                             g_0_yy_0_xyyyz_0,   \
                             g_0_yy_0_xyyyz_1,   \
                             g_0_yy_0_xyyzz_0,   \
                             g_0_yy_0_xyyzz_1,   \
                             g_0_yy_0_xyzzz_0,   \
                             g_0_yy_0_xyzzz_1,   \
                             g_0_yy_0_yyyyy_0,   \
                             g_0_yy_0_yyyyy_1,   \
                             g_0_yy_0_yyyyz_0,   \
                             g_0_yy_0_yyyyz_1,   \
                             g_0_yy_0_yyyzz_0,   \
                             g_0_yy_0_yyyzz_1,   \
                             g_0_yy_0_yyzzz_0,   \
                             g_0_yy_0_yyzzz_1,   \
                             g_0_yy_0_yzzzz_0,   \
                             g_0_yy_0_yzzzz_1,   \
                             g_0_yy_0_zzzzz_0,   \
                             g_0_yy_0_zzzzz_1,   \
                             wp_x,               \
                             wp_y,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyy_0_xxxxx_0[i] =
            g_0_xx_0_xxxxx_0[i] * fi_ab_0 - g_0_xx_0_xxxxx_1[i] * fti_ab_0 + g_0_xxy_0_xxxxx_0[i] * pb_y + g_0_xxy_0_xxxxx_1[i] * wp_y[i];

        g_0_xxyy_0_xxxxy_0[i] = g_0_yy_0_xxxxy_0[i] * fi_ab_0 - g_0_yy_0_xxxxy_1[i] * fti_ab_0 + 4.0 * g_0_xyy_0_xxxy_1[i] * fi_abcd_0 +
                                g_0_xyy_0_xxxxy_0[i] * pb_x + g_0_xyy_0_xxxxy_1[i] * wp_x[i];

        g_0_xxyy_0_xxxxz_0[i] =
            g_0_xx_0_xxxxz_0[i] * fi_ab_0 - g_0_xx_0_xxxxz_1[i] * fti_ab_0 + g_0_xxy_0_xxxxz_0[i] * pb_y + g_0_xxy_0_xxxxz_1[i] * wp_y[i];

        g_0_xxyy_0_xxxyy_0[i] = g_0_yy_0_xxxyy_0[i] * fi_ab_0 - g_0_yy_0_xxxyy_1[i] * fti_ab_0 + 3.0 * g_0_xyy_0_xxyy_1[i] * fi_abcd_0 +
                                g_0_xyy_0_xxxyy_0[i] * pb_x + g_0_xyy_0_xxxyy_1[i] * wp_x[i];

        g_0_xxyy_0_xxxyz_0[i] = g_0_yy_0_xxxyz_0[i] * fi_ab_0 - g_0_yy_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xyy_0_xxyz_1[i] * fi_abcd_0 +
                                g_0_xyy_0_xxxyz_0[i] * pb_x + g_0_xyy_0_xxxyz_1[i] * wp_x[i];

        g_0_xxyy_0_xxxzz_0[i] =
            g_0_xx_0_xxxzz_0[i] * fi_ab_0 - g_0_xx_0_xxxzz_1[i] * fti_ab_0 + g_0_xxy_0_xxxzz_0[i] * pb_y + g_0_xxy_0_xxxzz_1[i] * wp_y[i];

        g_0_xxyy_0_xxyyy_0[i] = g_0_yy_0_xxyyy_0[i] * fi_ab_0 - g_0_yy_0_xxyyy_1[i] * fti_ab_0 + 2.0 * g_0_xyy_0_xyyy_1[i] * fi_abcd_0 +
                                g_0_xyy_0_xxyyy_0[i] * pb_x + g_0_xyy_0_xxyyy_1[i] * wp_x[i];

        g_0_xxyy_0_xxyyz_0[i] = g_0_yy_0_xxyyz_0[i] * fi_ab_0 - g_0_yy_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xyy_0_xyyz_1[i] * fi_abcd_0 +
                                g_0_xyy_0_xxyyz_0[i] * pb_x + g_0_xyy_0_xxyyz_1[i] * wp_x[i];

        g_0_xxyy_0_xxyzz_0[i] = g_0_yy_0_xxyzz_0[i] * fi_ab_0 - g_0_yy_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xyy_0_xyzz_1[i] * fi_abcd_0 +
                                g_0_xyy_0_xxyzz_0[i] * pb_x + g_0_xyy_0_xxyzz_1[i] * wp_x[i];

        g_0_xxyy_0_xxzzz_0[i] =
            g_0_xx_0_xxzzz_0[i] * fi_ab_0 - g_0_xx_0_xxzzz_1[i] * fti_ab_0 + g_0_xxy_0_xxzzz_0[i] * pb_y + g_0_xxy_0_xxzzz_1[i] * wp_y[i];

        g_0_xxyy_0_xyyyy_0[i] = g_0_yy_0_xyyyy_0[i] * fi_ab_0 - g_0_yy_0_xyyyy_1[i] * fti_ab_0 + g_0_xyy_0_yyyy_1[i] * fi_abcd_0 +
                                g_0_xyy_0_xyyyy_0[i] * pb_x + g_0_xyy_0_xyyyy_1[i] * wp_x[i];

        g_0_xxyy_0_xyyyz_0[i] = g_0_yy_0_xyyyz_0[i] * fi_ab_0 - g_0_yy_0_xyyyz_1[i] * fti_ab_0 + g_0_xyy_0_yyyz_1[i] * fi_abcd_0 +
                                g_0_xyy_0_xyyyz_0[i] * pb_x + g_0_xyy_0_xyyyz_1[i] * wp_x[i];

        g_0_xxyy_0_xyyzz_0[i] = g_0_yy_0_xyyzz_0[i] * fi_ab_0 - g_0_yy_0_xyyzz_1[i] * fti_ab_0 + g_0_xyy_0_yyzz_1[i] * fi_abcd_0 +
                                g_0_xyy_0_xyyzz_0[i] * pb_x + g_0_xyy_0_xyyzz_1[i] * wp_x[i];

        g_0_xxyy_0_xyzzz_0[i] = g_0_yy_0_xyzzz_0[i] * fi_ab_0 - g_0_yy_0_xyzzz_1[i] * fti_ab_0 + g_0_xyy_0_yzzz_1[i] * fi_abcd_0 +
                                g_0_xyy_0_xyzzz_0[i] * pb_x + g_0_xyy_0_xyzzz_1[i] * wp_x[i];

        g_0_xxyy_0_xzzzz_0[i] =
            g_0_xx_0_xzzzz_0[i] * fi_ab_0 - g_0_xx_0_xzzzz_1[i] * fti_ab_0 + g_0_xxy_0_xzzzz_0[i] * pb_y + g_0_xxy_0_xzzzz_1[i] * wp_y[i];

        g_0_xxyy_0_yyyyy_0[i] =
            g_0_yy_0_yyyyy_0[i] * fi_ab_0 - g_0_yy_0_yyyyy_1[i] * fti_ab_0 + g_0_xyy_0_yyyyy_0[i] * pb_x + g_0_xyy_0_yyyyy_1[i] * wp_x[i];

        g_0_xxyy_0_yyyyz_0[i] =
            g_0_yy_0_yyyyz_0[i] * fi_ab_0 - g_0_yy_0_yyyyz_1[i] * fti_ab_0 + g_0_xyy_0_yyyyz_0[i] * pb_x + g_0_xyy_0_yyyyz_1[i] * wp_x[i];

        g_0_xxyy_0_yyyzz_0[i] =
            g_0_yy_0_yyyzz_0[i] * fi_ab_0 - g_0_yy_0_yyyzz_1[i] * fti_ab_0 + g_0_xyy_0_yyyzz_0[i] * pb_x + g_0_xyy_0_yyyzz_1[i] * wp_x[i];

        g_0_xxyy_0_yyzzz_0[i] =
            g_0_yy_0_yyzzz_0[i] * fi_ab_0 - g_0_yy_0_yyzzz_1[i] * fti_ab_0 + g_0_xyy_0_yyzzz_0[i] * pb_x + g_0_xyy_0_yyzzz_1[i] * wp_x[i];

        g_0_xxyy_0_yzzzz_0[i] =
            g_0_yy_0_yzzzz_0[i] * fi_ab_0 - g_0_yy_0_yzzzz_1[i] * fti_ab_0 + g_0_xyy_0_yzzzz_0[i] * pb_x + g_0_xyy_0_yzzzz_1[i] * wp_x[i];

        g_0_xxyy_0_zzzzz_0[i] =
            g_0_yy_0_zzzzz_0[i] * fi_ab_0 - g_0_yy_0_zzzzz_1[i] * fti_ab_0 + g_0_xyy_0_zzzzz_0[i] * pb_x + g_0_xyy_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 84-105 components of targeted buffer : SGSH

    auto g_0_xxyz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sgsh + 84);

    auto g_0_xxyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sgsh + 85);

    auto g_0_xxyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sgsh + 86);

    auto g_0_xxyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sgsh + 87);

    auto g_0_xxyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sgsh + 88);

    auto g_0_xxyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sgsh + 89);

    auto g_0_xxyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sgsh + 90);

    auto g_0_xxyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sgsh + 91);

    auto g_0_xxyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sgsh + 92);

    auto g_0_xxyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sgsh + 93);

    auto g_0_xxyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 94);

    auto g_0_xxyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 95);

    auto g_0_xxyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 96);

    auto g_0_xxyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 97);

    auto g_0_xxyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 98);

    auto g_0_xxyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 99);

    auto g_0_xxyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 100);

    auto g_0_xxyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 101);

    auto g_0_xxyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 102);

    auto g_0_xxyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 103);

    auto g_0_xxyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 104);

#pragma omp simd aligned(g_0_xxy_0_xxxxy_0,      \
                             g_0_xxy_0_xxxxy_1,  \
                             g_0_xxy_0_xxxyy_0,  \
                             g_0_xxy_0_xxxyy_1,  \
                             g_0_xxy_0_xxyyy_0,  \
                             g_0_xxy_0_xxyyy_1,  \
                             g_0_xxy_0_xyyyy_0,  \
                             g_0_xxy_0_xyyyy_1,  \
                             g_0_xxy_0_yyyyy_0,  \
                             g_0_xxy_0_yyyyy_1,  \
                             g_0_xxyz_0_xxxxx_0, \
                             g_0_xxyz_0_xxxxy_0, \
                             g_0_xxyz_0_xxxxz_0, \
                             g_0_xxyz_0_xxxyy_0, \
                             g_0_xxyz_0_xxxyz_0, \
                             g_0_xxyz_0_xxxzz_0, \
                             g_0_xxyz_0_xxyyy_0, \
                             g_0_xxyz_0_xxyyz_0, \
                             g_0_xxyz_0_xxyzz_0, \
                             g_0_xxyz_0_xxzzz_0, \
                             g_0_xxyz_0_xyyyy_0, \
                             g_0_xxyz_0_xyyyz_0, \
                             g_0_xxyz_0_xyyzz_0, \
                             g_0_xxyz_0_xyzzz_0, \
                             g_0_xxyz_0_xzzzz_0, \
                             g_0_xxyz_0_yyyyy_0, \
                             g_0_xxyz_0_yyyyz_0, \
                             g_0_xxyz_0_yyyzz_0, \
                             g_0_xxyz_0_yyzzz_0, \
                             g_0_xxyz_0_yzzzz_0, \
                             g_0_xxyz_0_zzzzz_0, \
                             g_0_xxz_0_xxxxx_0,  \
                             g_0_xxz_0_xxxxx_1,  \
                             g_0_xxz_0_xxxxz_0,  \
                             g_0_xxz_0_xxxxz_1,  \
                             g_0_xxz_0_xxxyz_0,  \
                             g_0_xxz_0_xxxyz_1,  \
                             g_0_xxz_0_xxxz_1,   \
                             g_0_xxz_0_xxxzz_0,  \
                             g_0_xxz_0_xxxzz_1,  \
                             g_0_xxz_0_xxyyz_0,  \
                             g_0_xxz_0_xxyyz_1,  \
                             g_0_xxz_0_xxyz_1,   \
                             g_0_xxz_0_xxyzz_0,  \
                             g_0_xxz_0_xxyzz_1,  \
                             g_0_xxz_0_xxzz_1,   \
                             g_0_xxz_0_xxzzz_0,  \
                             g_0_xxz_0_xxzzz_1,  \
                             g_0_xxz_0_xyyyz_0,  \
                             g_0_xxz_0_xyyyz_1,  \
                             g_0_xxz_0_xyyz_1,   \
                             g_0_xxz_0_xyyzz_0,  \
                             g_0_xxz_0_xyyzz_1,  \
                             g_0_xxz_0_xyzz_1,   \
                             g_0_xxz_0_xyzzz_0,  \
                             g_0_xxz_0_xyzzz_1,  \
                             g_0_xxz_0_xzzz_1,   \
                             g_0_xxz_0_xzzzz_0,  \
                             g_0_xxz_0_xzzzz_1,  \
                             g_0_xxz_0_yyyyz_0,  \
                             g_0_xxz_0_yyyyz_1,  \
                             g_0_xxz_0_yyyz_1,   \
                             g_0_xxz_0_yyyzz_0,  \
                             g_0_xxz_0_yyyzz_1,  \
                             g_0_xxz_0_yyzz_1,   \
                             g_0_xxz_0_yyzzz_0,  \
                             g_0_xxz_0_yyzzz_1,  \
                             g_0_xxz_0_yzzz_1,   \
                             g_0_xxz_0_yzzzz_0,  \
                             g_0_xxz_0_yzzzz_1,  \
                             g_0_xxz_0_zzzz_1,   \
                             g_0_xxz_0_zzzzz_0,  \
                             g_0_xxz_0_zzzzz_1,  \
                             wp_y,               \
                             wp_z,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyz_0_xxxxx_0[i] = g_0_xxz_0_xxxxx_0[i] * pb_y + g_0_xxz_0_xxxxx_1[i] * wp_y[i];

        g_0_xxyz_0_xxxxy_0[i] = g_0_xxy_0_xxxxy_0[i] * pb_z + g_0_xxy_0_xxxxy_1[i] * wp_z[i];

        g_0_xxyz_0_xxxxz_0[i] = g_0_xxz_0_xxxxz_0[i] * pb_y + g_0_xxz_0_xxxxz_1[i] * wp_y[i];

        g_0_xxyz_0_xxxyy_0[i] = g_0_xxy_0_xxxyy_0[i] * pb_z + g_0_xxy_0_xxxyy_1[i] * wp_z[i];

        g_0_xxyz_0_xxxyz_0[i] = g_0_xxz_0_xxxz_1[i] * fi_abcd_0 + g_0_xxz_0_xxxyz_0[i] * pb_y + g_0_xxz_0_xxxyz_1[i] * wp_y[i];

        g_0_xxyz_0_xxxzz_0[i] = g_0_xxz_0_xxxzz_0[i] * pb_y + g_0_xxz_0_xxxzz_1[i] * wp_y[i];

        g_0_xxyz_0_xxyyy_0[i] = g_0_xxy_0_xxyyy_0[i] * pb_z + g_0_xxy_0_xxyyy_1[i] * wp_z[i];

        g_0_xxyz_0_xxyyz_0[i] = 2.0 * g_0_xxz_0_xxyz_1[i] * fi_abcd_0 + g_0_xxz_0_xxyyz_0[i] * pb_y + g_0_xxz_0_xxyyz_1[i] * wp_y[i];

        g_0_xxyz_0_xxyzz_0[i] = g_0_xxz_0_xxzz_1[i] * fi_abcd_0 + g_0_xxz_0_xxyzz_0[i] * pb_y + g_0_xxz_0_xxyzz_1[i] * wp_y[i];

        g_0_xxyz_0_xxzzz_0[i] = g_0_xxz_0_xxzzz_0[i] * pb_y + g_0_xxz_0_xxzzz_1[i] * wp_y[i];

        g_0_xxyz_0_xyyyy_0[i] = g_0_xxy_0_xyyyy_0[i] * pb_z + g_0_xxy_0_xyyyy_1[i] * wp_z[i];

        g_0_xxyz_0_xyyyz_0[i] = 3.0 * g_0_xxz_0_xyyz_1[i] * fi_abcd_0 + g_0_xxz_0_xyyyz_0[i] * pb_y + g_0_xxz_0_xyyyz_1[i] * wp_y[i];

        g_0_xxyz_0_xyyzz_0[i] = 2.0 * g_0_xxz_0_xyzz_1[i] * fi_abcd_0 + g_0_xxz_0_xyyzz_0[i] * pb_y + g_0_xxz_0_xyyzz_1[i] * wp_y[i];

        g_0_xxyz_0_xyzzz_0[i] = g_0_xxz_0_xzzz_1[i] * fi_abcd_0 + g_0_xxz_0_xyzzz_0[i] * pb_y + g_0_xxz_0_xyzzz_1[i] * wp_y[i];

        g_0_xxyz_0_xzzzz_0[i] = g_0_xxz_0_xzzzz_0[i] * pb_y + g_0_xxz_0_xzzzz_1[i] * wp_y[i];

        g_0_xxyz_0_yyyyy_0[i] = g_0_xxy_0_yyyyy_0[i] * pb_z + g_0_xxy_0_yyyyy_1[i] * wp_z[i];

        g_0_xxyz_0_yyyyz_0[i] = 4.0 * g_0_xxz_0_yyyz_1[i] * fi_abcd_0 + g_0_xxz_0_yyyyz_0[i] * pb_y + g_0_xxz_0_yyyyz_1[i] * wp_y[i];

        g_0_xxyz_0_yyyzz_0[i] = 3.0 * g_0_xxz_0_yyzz_1[i] * fi_abcd_0 + g_0_xxz_0_yyyzz_0[i] * pb_y + g_0_xxz_0_yyyzz_1[i] * wp_y[i];

        g_0_xxyz_0_yyzzz_0[i] = 2.0 * g_0_xxz_0_yzzz_1[i] * fi_abcd_0 + g_0_xxz_0_yyzzz_0[i] * pb_y + g_0_xxz_0_yyzzz_1[i] * wp_y[i];

        g_0_xxyz_0_yzzzz_0[i] = g_0_xxz_0_zzzz_1[i] * fi_abcd_0 + g_0_xxz_0_yzzzz_0[i] * pb_y + g_0_xxz_0_yzzzz_1[i] * wp_y[i];

        g_0_xxyz_0_zzzzz_0[i] = g_0_xxz_0_zzzzz_0[i] * pb_y + g_0_xxz_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 105-126 components of targeted buffer : SGSH

    auto g_0_xxzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sgsh + 105);

    auto g_0_xxzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sgsh + 106);

    auto g_0_xxzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sgsh + 107);

    auto g_0_xxzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sgsh + 108);

    auto g_0_xxzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sgsh + 109);

    auto g_0_xxzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sgsh + 110);

    auto g_0_xxzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sgsh + 111);

    auto g_0_xxzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sgsh + 112);

    auto g_0_xxzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sgsh + 113);

    auto g_0_xxzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sgsh + 114);

    auto g_0_xxzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 115);

    auto g_0_xxzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 116);

    auto g_0_xxzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 117);

    auto g_0_xxzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 118);

    auto g_0_xxzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 119);

    auto g_0_xxzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 120);

    auto g_0_xxzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 121);

    auto g_0_xxzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 122);

    auto g_0_xxzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 123);

    auto g_0_xxzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 124);

    auto g_0_xxzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 125);

#pragma omp simd aligned(g_0_xx_0_xxxxx_0,       \
                             g_0_xx_0_xxxxx_1,   \
                             g_0_xx_0_xxxxy_0,   \
                             g_0_xx_0_xxxxy_1,   \
                             g_0_xx_0_xxxyy_0,   \
                             g_0_xx_0_xxxyy_1,   \
                             g_0_xx_0_xxyyy_0,   \
                             g_0_xx_0_xxyyy_1,   \
                             g_0_xx_0_xyyyy_0,   \
                             g_0_xx_0_xyyyy_1,   \
                             g_0_xxz_0_xxxxx_0,  \
                             g_0_xxz_0_xxxxx_1,  \
                             g_0_xxz_0_xxxxy_0,  \
                             g_0_xxz_0_xxxxy_1,  \
                             g_0_xxz_0_xxxyy_0,  \
                             g_0_xxz_0_xxxyy_1,  \
                             g_0_xxz_0_xxyyy_0,  \
                             g_0_xxz_0_xxyyy_1,  \
                             g_0_xxz_0_xyyyy_0,  \
                             g_0_xxz_0_xyyyy_1,  \
                             g_0_xxzz_0_xxxxx_0, \
                             g_0_xxzz_0_xxxxy_0, \
                             g_0_xxzz_0_xxxxz_0, \
                             g_0_xxzz_0_xxxyy_0, \
                             g_0_xxzz_0_xxxyz_0, \
                             g_0_xxzz_0_xxxzz_0, \
                             g_0_xxzz_0_xxyyy_0, \
                             g_0_xxzz_0_xxyyz_0, \
                             g_0_xxzz_0_xxyzz_0, \
                             g_0_xxzz_0_xxzzz_0, \
                             g_0_xxzz_0_xyyyy_0, \
                             g_0_xxzz_0_xyyyz_0, \
                             g_0_xxzz_0_xyyzz_0, \
                             g_0_xxzz_0_xyzzz_0, \
                             g_0_xxzz_0_xzzzz_0, \
                             g_0_xxzz_0_yyyyy_0, \
                             g_0_xxzz_0_yyyyz_0, \
                             g_0_xxzz_0_yyyzz_0, \
                             g_0_xxzz_0_yyzzz_0, \
                             g_0_xxzz_0_yzzzz_0, \
                             g_0_xxzz_0_zzzzz_0, \
                             g_0_xzz_0_xxxxz_0,  \
                             g_0_xzz_0_xxxxz_1,  \
                             g_0_xzz_0_xxxyz_0,  \
                             g_0_xzz_0_xxxyz_1,  \
                             g_0_xzz_0_xxxz_1,   \
                             g_0_xzz_0_xxxzz_0,  \
                             g_0_xzz_0_xxxzz_1,  \
                             g_0_xzz_0_xxyyz_0,  \
                             g_0_xzz_0_xxyyz_1,  \
                             g_0_xzz_0_xxyz_1,   \
                             g_0_xzz_0_xxyzz_0,  \
                             g_0_xzz_0_xxyzz_1,  \
                             g_0_xzz_0_xxzz_1,   \
                             g_0_xzz_0_xxzzz_0,  \
                             g_0_xzz_0_xxzzz_1,  \
                             g_0_xzz_0_xyyyz_0,  \
                             g_0_xzz_0_xyyyz_1,  \
                             g_0_xzz_0_xyyz_1,   \
                             g_0_xzz_0_xyyzz_0,  \
                             g_0_xzz_0_xyyzz_1,  \
                             g_0_xzz_0_xyzz_1,   \
                             g_0_xzz_0_xyzzz_0,  \
                             g_0_xzz_0_xyzzz_1,  \
                             g_0_xzz_0_xzzz_1,   \
                             g_0_xzz_0_xzzzz_0,  \
                             g_0_xzz_0_xzzzz_1,  \
                             g_0_xzz_0_yyyyy_0,  \
                             g_0_xzz_0_yyyyy_1,  \
                             g_0_xzz_0_yyyyz_0,  \
                             g_0_xzz_0_yyyyz_1,  \
                             g_0_xzz_0_yyyz_1,   \
                             g_0_xzz_0_yyyzz_0,  \
                             g_0_xzz_0_yyyzz_1,  \
                             g_0_xzz_0_yyzz_1,   \
                             g_0_xzz_0_yyzzz_0,  \
                             g_0_xzz_0_yyzzz_1,  \
                             g_0_xzz_0_yzzz_1,   \
                             g_0_xzz_0_yzzzz_0,  \
                             g_0_xzz_0_yzzzz_1,  \
                             g_0_xzz_0_zzzz_1,   \
                             g_0_xzz_0_zzzzz_0,  \
                             g_0_xzz_0_zzzzz_1,  \
                             g_0_zz_0_xxxxz_0,   \
                             g_0_zz_0_xxxxz_1,   \
                             g_0_zz_0_xxxyz_0,   \
                             g_0_zz_0_xxxyz_1,   \
                             g_0_zz_0_xxxzz_0,   \
                             g_0_zz_0_xxxzz_1,   \
                             g_0_zz_0_xxyyz_0,   \
                             g_0_zz_0_xxyyz_1,   \
                             g_0_zz_0_xxyzz_0,   \
                             g_0_zz_0_xxyzz_1,   \
                             g_0_zz_0_xxzzz_0,   \
                             g_0_zz_0_xxzzz_1,   \
                             g_0_zz_0_xyyyz_0,   \
                             g_0_zz_0_xyyyz_1,   \
                             g_0_zz_0_xyyzz_0,   \
                             g_0_zz_0_xyyzz_1,   \
                             g_0_zz_0_xyzzz_0,   \
                             g_0_zz_0_xyzzz_1,   \
                             g_0_zz_0_xzzzz_0,   \
                             g_0_zz_0_xzzzz_1,   \
                             g_0_zz_0_yyyyy_0,   \
                             g_0_zz_0_yyyyy_1,   \
                             g_0_zz_0_yyyyz_0,   \
                             g_0_zz_0_yyyyz_1,   \
                             g_0_zz_0_yyyzz_0,   \
                             g_0_zz_0_yyyzz_1,   \
                             g_0_zz_0_yyzzz_0,   \
                             g_0_zz_0_yyzzz_1,   \
                             g_0_zz_0_yzzzz_0,   \
                             g_0_zz_0_yzzzz_1,   \
                             g_0_zz_0_zzzzz_0,   \
                             g_0_zz_0_zzzzz_1,   \
                             wp_x,               \
                             wp_z,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxzz_0_xxxxx_0[i] =
            g_0_xx_0_xxxxx_0[i] * fi_ab_0 - g_0_xx_0_xxxxx_1[i] * fti_ab_0 + g_0_xxz_0_xxxxx_0[i] * pb_z + g_0_xxz_0_xxxxx_1[i] * wp_z[i];

        g_0_xxzz_0_xxxxy_0[i] =
            g_0_xx_0_xxxxy_0[i] * fi_ab_0 - g_0_xx_0_xxxxy_1[i] * fti_ab_0 + g_0_xxz_0_xxxxy_0[i] * pb_z + g_0_xxz_0_xxxxy_1[i] * wp_z[i];

        g_0_xxzz_0_xxxxz_0[i] = g_0_zz_0_xxxxz_0[i] * fi_ab_0 - g_0_zz_0_xxxxz_1[i] * fti_ab_0 + 4.0 * g_0_xzz_0_xxxz_1[i] * fi_abcd_0 +
                                g_0_xzz_0_xxxxz_0[i] * pb_x + g_0_xzz_0_xxxxz_1[i] * wp_x[i];

        g_0_xxzz_0_xxxyy_0[i] =
            g_0_xx_0_xxxyy_0[i] * fi_ab_0 - g_0_xx_0_xxxyy_1[i] * fti_ab_0 + g_0_xxz_0_xxxyy_0[i] * pb_z + g_0_xxz_0_xxxyy_1[i] * wp_z[i];

        g_0_xxzz_0_xxxyz_0[i] = g_0_zz_0_xxxyz_0[i] * fi_ab_0 - g_0_zz_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xzz_0_xxyz_1[i] * fi_abcd_0 +
                                g_0_xzz_0_xxxyz_0[i] * pb_x + g_0_xzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xxzz_0_xxxzz_0[i] = g_0_zz_0_xxxzz_0[i] * fi_ab_0 - g_0_zz_0_xxxzz_1[i] * fti_ab_0 + 3.0 * g_0_xzz_0_xxzz_1[i] * fi_abcd_0 +
                                g_0_xzz_0_xxxzz_0[i] * pb_x + g_0_xzz_0_xxxzz_1[i] * wp_x[i];

        g_0_xxzz_0_xxyyy_0[i] =
            g_0_xx_0_xxyyy_0[i] * fi_ab_0 - g_0_xx_0_xxyyy_1[i] * fti_ab_0 + g_0_xxz_0_xxyyy_0[i] * pb_z + g_0_xxz_0_xxyyy_1[i] * wp_z[i];

        g_0_xxzz_0_xxyyz_0[i] = g_0_zz_0_xxyyz_0[i] * fi_ab_0 - g_0_zz_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xzz_0_xyyz_1[i] * fi_abcd_0 +
                                g_0_xzz_0_xxyyz_0[i] * pb_x + g_0_xzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xxzz_0_xxyzz_0[i] = g_0_zz_0_xxyzz_0[i] * fi_ab_0 - g_0_zz_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xzz_0_xyzz_1[i] * fi_abcd_0 +
                                g_0_xzz_0_xxyzz_0[i] * pb_x + g_0_xzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xxzz_0_xxzzz_0[i] = g_0_zz_0_xxzzz_0[i] * fi_ab_0 - g_0_zz_0_xxzzz_1[i] * fti_ab_0 + 2.0 * g_0_xzz_0_xzzz_1[i] * fi_abcd_0 +
                                g_0_xzz_0_xxzzz_0[i] * pb_x + g_0_xzz_0_xxzzz_1[i] * wp_x[i];

        g_0_xxzz_0_xyyyy_0[i] =
            g_0_xx_0_xyyyy_0[i] * fi_ab_0 - g_0_xx_0_xyyyy_1[i] * fti_ab_0 + g_0_xxz_0_xyyyy_0[i] * pb_z + g_0_xxz_0_xyyyy_1[i] * wp_z[i];

        g_0_xxzz_0_xyyyz_0[i] = g_0_zz_0_xyyyz_0[i] * fi_ab_0 - g_0_zz_0_xyyyz_1[i] * fti_ab_0 + g_0_xzz_0_yyyz_1[i] * fi_abcd_0 +
                                g_0_xzz_0_xyyyz_0[i] * pb_x + g_0_xzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xxzz_0_xyyzz_0[i] = g_0_zz_0_xyyzz_0[i] * fi_ab_0 - g_0_zz_0_xyyzz_1[i] * fti_ab_0 + g_0_xzz_0_yyzz_1[i] * fi_abcd_0 +
                                g_0_xzz_0_xyyzz_0[i] * pb_x + g_0_xzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xxzz_0_xyzzz_0[i] = g_0_zz_0_xyzzz_0[i] * fi_ab_0 - g_0_zz_0_xyzzz_1[i] * fti_ab_0 + g_0_xzz_0_yzzz_1[i] * fi_abcd_0 +
                                g_0_xzz_0_xyzzz_0[i] * pb_x + g_0_xzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xxzz_0_xzzzz_0[i] = g_0_zz_0_xzzzz_0[i] * fi_ab_0 - g_0_zz_0_xzzzz_1[i] * fti_ab_0 + g_0_xzz_0_zzzz_1[i] * fi_abcd_0 +
                                g_0_xzz_0_xzzzz_0[i] * pb_x + g_0_xzz_0_xzzzz_1[i] * wp_x[i];

        g_0_xxzz_0_yyyyy_0[i] =
            g_0_zz_0_yyyyy_0[i] * fi_ab_0 - g_0_zz_0_yyyyy_1[i] * fti_ab_0 + g_0_xzz_0_yyyyy_0[i] * pb_x + g_0_xzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xxzz_0_yyyyz_0[i] =
            g_0_zz_0_yyyyz_0[i] * fi_ab_0 - g_0_zz_0_yyyyz_1[i] * fti_ab_0 + g_0_xzz_0_yyyyz_0[i] * pb_x + g_0_xzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xxzz_0_yyyzz_0[i] =
            g_0_zz_0_yyyzz_0[i] * fi_ab_0 - g_0_zz_0_yyyzz_1[i] * fti_ab_0 + g_0_xzz_0_yyyzz_0[i] * pb_x + g_0_xzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xxzz_0_yyzzz_0[i] =
            g_0_zz_0_yyzzz_0[i] * fi_ab_0 - g_0_zz_0_yyzzz_1[i] * fti_ab_0 + g_0_xzz_0_yyzzz_0[i] * pb_x + g_0_xzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xxzz_0_yzzzz_0[i] =
            g_0_zz_0_yzzzz_0[i] * fi_ab_0 - g_0_zz_0_yzzzz_1[i] * fti_ab_0 + g_0_xzz_0_yzzzz_0[i] * pb_x + g_0_xzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xxzz_0_zzzzz_0[i] =
            g_0_zz_0_zzzzz_0[i] * fi_ab_0 - g_0_zz_0_zzzzz_1[i] * fti_ab_0 + g_0_xzz_0_zzzzz_0[i] * pb_x + g_0_xzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 126-147 components of targeted buffer : SGSH

    auto g_0_xyyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_sgsh + 126);

    auto g_0_xyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_sgsh + 127);

    auto g_0_xyyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_sgsh + 128);

    auto g_0_xyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_sgsh + 129);

    auto g_0_xyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_sgsh + 130);

    auto g_0_xyyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_sgsh + 131);

    auto g_0_xyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_sgsh + 132);

    auto g_0_xyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_sgsh + 133);

    auto g_0_xyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_sgsh + 134);

    auto g_0_xyyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_sgsh + 135);

    auto g_0_xyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 136);

    auto g_0_xyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 137);

    auto g_0_xyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 138);

    auto g_0_xyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 139);

    auto g_0_xyyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 140);

    auto g_0_xyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 141);

    auto g_0_xyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 142);

    auto g_0_xyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 143);

    auto g_0_xyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 144);

    auto g_0_xyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 145);

    auto g_0_xyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 146);

#pragma omp simd aligned(g_0_xyyy_0_xxxxx_0,     \
                             g_0_xyyy_0_xxxxy_0, \
                             g_0_xyyy_0_xxxxz_0, \
                             g_0_xyyy_0_xxxyy_0, \
                             g_0_xyyy_0_xxxyz_0, \
                             g_0_xyyy_0_xxxzz_0, \
                             g_0_xyyy_0_xxyyy_0, \
                             g_0_xyyy_0_xxyyz_0, \
                             g_0_xyyy_0_xxyzz_0, \
                             g_0_xyyy_0_xxzzz_0, \
                             g_0_xyyy_0_xyyyy_0, \
                             g_0_xyyy_0_xyyyz_0, \
                             g_0_xyyy_0_xyyzz_0, \
                             g_0_xyyy_0_xyzzz_0, \
                             g_0_xyyy_0_xzzzz_0, \
                             g_0_xyyy_0_yyyyy_0, \
                             g_0_xyyy_0_yyyyz_0, \
                             g_0_xyyy_0_yyyzz_0, \
                             g_0_xyyy_0_yyzzz_0, \
                             g_0_xyyy_0_yzzzz_0, \
                             g_0_xyyy_0_zzzzz_0, \
                             g_0_yyy_0_xxxx_1,   \
                             g_0_yyy_0_xxxxx_0,  \
                             g_0_yyy_0_xxxxx_1,  \
                             g_0_yyy_0_xxxxy_0,  \
                             g_0_yyy_0_xxxxy_1,  \
                             g_0_yyy_0_xxxxz_0,  \
                             g_0_yyy_0_xxxxz_1,  \
                             g_0_yyy_0_xxxy_1,   \
                             g_0_yyy_0_xxxyy_0,  \
                             g_0_yyy_0_xxxyy_1,  \
                             g_0_yyy_0_xxxyz_0,  \
                             g_0_yyy_0_xxxyz_1,  \
                             g_0_yyy_0_xxxz_1,   \
                             g_0_yyy_0_xxxzz_0,  \
                             g_0_yyy_0_xxxzz_1,  \
                             g_0_yyy_0_xxyy_1,   \
                             g_0_yyy_0_xxyyy_0,  \
                             g_0_yyy_0_xxyyy_1,  \
                             g_0_yyy_0_xxyyz_0,  \
                             g_0_yyy_0_xxyyz_1,  \
                             g_0_yyy_0_xxyz_1,   \
                             g_0_yyy_0_xxyzz_0,  \
                             g_0_yyy_0_xxyzz_1,  \
                             g_0_yyy_0_xxzz_1,   \
                             g_0_yyy_0_xxzzz_0,  \
                             g_0_yyy_0_xxzzz_1,  \
                             g_0_yyy_0_xyyy_1,   \
                             g_0_yyy_0_xyyyy_0,  \
                             g_0_yyy_0_xyyyy_1,  \
                             g_0_yyy_0_xyyyz_0,  \
                             g_0_yyy_0_xyyyz_1,  \
                             g_0_yyy_0_xyyz_1,   \
                             g_0_yyy_0_xyyzz_0,  \
                             g_0_yyy_0_xyyzz_1,  \
                             g_0_yyy_0_xyzz_1,   \
                             g_0_yyy_0_xyzzz_0,  \
                             g_0_yyy_0_xyzzz_1,  \
                             g_0_yyy_0_xzzz_1,   \
                             g_0_yyy_0_xzzzz_0,  \
                             g_0_yyy_0_xzzzz_1,  \
                             g_0_yyy_0_yyyy_1,   \
                             g_0_yyy_0_yyyyy_0,  \
                             g_0_yyy_0_yyyyy_1,  \
                             g_0_yyy_0_yyyyz_0,  \
                             g_0_yyy_0_yyyyz_1,  \
                             g_0_yyy_0_yyyz_1,   \
                             g_0_yyy_0_yyyzz_0,  \
                             g_0_yyy_0_yyyzz_1,  \
                             g_0_yyy_0_yyzz_1,   \
                             g_0_yyy_0_yyzzz_0,  \
                             g_0_yyy_0_yyzzz_1,  \
                             g_0_yyy_0_yzzz_1,   \
                             g_0_yyy_0_yzzzz_0,  \
                             g_0_yyy_0_yzzzz_1,  \
                             g_0_yyy_0_zzzz_1,   \
                             g_0_yyy_0_zzzzz_0,  \
                             g_0_yyy_0_zzzzz_1,  \
                             wp_x,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyy_0_xxxxx_0[i] = 5.0 * g_0_yyy_0_xxxx_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxx_0[i] * pb_x + g_0_yyy_0_xxxxx_1[i] * wp_x[i];

        g_0_xyyy_0_xxxxy_0[i] = 4.0 * g_0_yyy_0_xxxy_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxy_0[i] * pb_x + g_0_yyy_0_xxxxy_1[i] * wp_x[i];

        g_0_xyyy_0_xxxxz_0[i] = 4.0 * g_0_yyy_0_xxxz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxz_0[i] * pb_x + g_0_yyy_0_xxxxz_1[i] * wp_x[i];

        g_0_xyyy_0_xxxyy_0[i] = 3.0 * g_0_yyy_0_xxyy_1[i] * fi_abcd_0 + g_0_yyy_0_xxxyy_0[i] * pb_x + g_0_yyy_0_xxxyy_1[i] * wp_x[i];

        g_0_xyyy_0_xxxyz_0[i] = 3.0 * g_0_yyy_0_xxyz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxyz_0[i] * pb_x + g_0_yyy_0_xxxyz_1[i] * wp_x[i];

        g_0_xyyy_0_xxxzz_0[i] = 3.0 * g_0_yyy_0_xxzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxzz_0[i] * pb_x + g_0_yyy_0_xxxzz_1[i] * wp_x[i];

        g_0_xyyy_0_xxyyy_0[i] = 2.0 * g_0_yyy_0_xyyy_1[i] * fi_abcd_0 + g_0_yyy_0_xxyyy_0[i] * pb_x + g_0_yyy_0_xxyyy_1[i] * wp_x[i];

        g_0_xyyy_0_xxyyz_0[i] = 2.0 * g_0_yyy_0_xyyz_1[i] * fi_abcd_0 + g_0_yyy_0_xxyyz_0[i] * pb_x + g_0_yyy_0_xxyyz_1[i] * wp_x[i];

        g_0_xyyy_0_xxyzz_0[i] = 2.0 * g_0_yyy_0_xyzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxyzz_0[i] * pb_x + g_0_yyy_0_xxyzz_1[i] * wp_x[i];

        g_0_xyyy_0_xxzzz_0[i] = 2.0 * g_0_yyy_0_xzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxzzz_0[i] * pb_x + g_0_yyy_0_xxzzz_1[i] * wp_x[i];

        g_0_xyyy_0_xyyyy_0[i] = g_0_yyy_0_yyyy_1[i] * fi_abcd_0 + g_0_yyy_0_xyyyy_0[i] * pb_x + g_0_yyy_0_xyyyy_1[i] * wp_x[i];

        g_0_xyyy_0_xyyyz_0[i] = g_0_yyy_0_yyyz_1[i] * fi_abcd_0 + g_0_yyy_0_xyyyz_0[i] * pb_x + g_0_yyy_0_xyyyz_1[i] * wp_x[i];

        g_0_xyyy_0_xyyzz_0[i] = g_0_yyy_0_yyzz_1[i] * fi_abcd_0 + g_0_yyy_0_xyyzz_0[i] * pb_x + g_0_yyy_0_xyyzz_1[i] * wp_x[i];

        g_0_xyyy_0_xyzzz_0[i] = g_0_yyy_0_yzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xyzzz_0[i] * pb_x + g_0_yyy_0_xyzzz_1[i] * wp_x[i];

        g_0_xyyy_0_xzzzz_0[i] = g_0_yyy_0_zzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xzzzz_0[i] * pb_x + g_0_yyy_0_xzzzz_1[i] * wp_x[i];

        g_0_xyyy_0_yyyyy_0[i] = g_0_yyy_0_yyyyy_0[i] * pb_x + g_0_yyy_0_yyyyy_1[i] * wp_x[i];

        g_0_xyyy_0_yyyyz_0[i] = g_0_yyy_0_yyyyz_0[i] * pb_x + g_0_yyy_0_yyyyz_1[i] * wp_x[i];

        g_0_xyyy_0_yyyzz_0[i] = g_0_yyy_0_yyyzz_0[i] * pb_x + g_0_yyy_0_yyyzz_1[i] * wp_x[i];

        g_0_xyyy_0_yyzzz_0[i] = g_0_yyy_0_yyzzz_0[i] * pb_x + g_0_yyy_0_yyzzz_1[i] * wp_x[i];

        g_0_xyyy_0_yzzzz_0[i] = g_0_yyy_0_yzzzz_0[i] * pb_x + g_0_yyy_0_yzzzz_1[i] * wp_x[i];

        g_0_xyyy_0_zzzzz_0[i] = g_0_yyy_0_zzzzz_0[i] * pb_x + g_0_yyy_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 147-168 components of targeted buffer : SGSH

    auto g_0_xyyz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sgsh + 147);

    auto g_0_xyyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sgsh + 148);

    auto g_0_xyyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sgsh + 149);

    auto g_0_xyyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sgsh + 150);

    auto g_0_xyyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sgsh + 151);

    auto g_0_xyyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sgsh + 152);

    auto g_0_xyyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sgsh + 153);

    auto g_0_xyyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sgsh + 154);

    auto g_0_xyyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sgsh + 155);

    auto g_0_xyyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sgsh + 156);

    auto g_0_xyyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 157);

    auto g_0_xyyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 158);

    auto g_0_xyyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 159);

    auto g_0_xyyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 160);

    auto g_0_xyyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 161);

    auto g_0_xyyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 162);

    auto g_0_xyyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 163);

    auto g_0_xyyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 164);

    auto g_0_xyyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 165);

    auto g_0_xyyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 166);

    auto g_0_xyyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 167);

#pragma omp simd aligned(g_0_xyy_0_xxxxx_0,      \
                             g_0_xyy_0_xxxxx_1,  \
                             g_0_xyy_0_xxxxy_0,  \
                             g_0_xyy_0_xxxxy_1,  \
                             g_0_xyy_0_xxxyy_0,  \
                             g_0_xyy_0_xxxyy_1,  \
                             g_0_xyy_0_xxyyy_0,  \
                             g_0_xyy_0_xxyyy_1,  \
                             g_0_xyy_0_xyyyy_0,  \
                             g_0_xyy_0_xyyyy_1,  \
                             g_0_xyyz_0_xxxxx_0, \
                             g_0_xyyz_0_xxxxy_0, \
                             g_0_xyyz_0_xxxxz_0, \
                             g_0_xyyz_0_xxxyy_0, \
                             g_0_xyyz_0_xxxyz_0, \
                             g_0_xyyz_0_xxxzz_0, \
                             g_0_xyyz_0_xxyyy_0, \
                             g_0_xyyz_0_xxyyz_0, \
                             g_0_xyyz_0_xxyzz_0, \
                             g_0_xyyz_0_xxzzz_0, \
                             g_0_xyyz_0_xyyyy_0, \
                             g_0_xyyz_0_xyyyz_0, \
                             g_0_xyyz_0_xyyzz_0, \
                             g_0_xyyz_0_xyzzz_0, \
                             g_0_xyyz_0_xzzzz_0, \
                             g_0_xyyz_0_yyyyy_0, \
                             g_0_xyyz_0_yyyyz_0, \
                             g_0_xyyz_0_yyyzz_0, \
                             g_0_xyyz_0_yyzzz_0, \
                             g_0_xyyz_0_yzzzz_0, \
                             g_0_xyyz_0_zzzzz_0, \
                             g_0_yyz_0_xxxxz_0,  \
                             g_0_yyz_0_xxxxz_1,  \
                             g_0_yyz_0_xxxyz_0,  \
                             g_0_yyz_0_xxxyz_1,  \
                             g_0_yyz_0_xxxz_1,   \
                             g_0_yyz_0_xxxzz_0,  \
                             g_0_yyz_0_xxxzz_1,  \
                             g_0_yyz_0_xxyyz_0,  \
                             g_0_yyz_0_xxyyz_1,  \
                             g_0_yyz_0_xxyz_1,   \
                             g_0_yyz_0_xxyzz_0,  \
                             g_0_yyz_0_xxyzz_1,  \
                             g_0_yyz_0_xxzz_1,   \
                             g_0_yyz_0_xxzzz_0,  \
                             g_0_yyz_0_xxzzz_1,  \
                             g_0_yyz_0_xyyyz_0,  \
                             g_0_yyz_0_xyyyz_1,  \
                             g_0_yyz_0_xyyz_1,   \
                             g_0_yyz_0_xyyzz_0,  \
                             g_0_yyz_0_xyyzz_1,  \
                             g_0_yyz_0_xyzz_1,   \
                             g_0_yyz_0_xyzzz_0,  \
                             g_0_yyz_0_xyzzz_1,  \
                             g_0_yyz_0_xzzz_1,   \
                             g_0_yyz_0_xzzzz_0,  \
                             g_0_yyz_0_xzzzz_1,  \
                             g_0_yyz_0_yyyyy_0,  \
                             g_0_yyz_0_yyyyy_1,  \
                             g_0_yyz_0_yyyyz_0,  \
                             g_0_yyz_0_yyyyz_1,  \
                             g_0_yyz_0_yyyz_1,   \
                             g_0_yyz_0_yyyzz_0,  \
                             g_0_yyz_0_yyyzz_1,  \
                             g_0_yyz_0_yyzz_1,   \
                             g_0_yyz_0_yyzzz_0,  \
                             g_0_yyz_0_yyzzz_1,  \
                             g_0_yyz_0_yzzz_1,   \
                             g_0_yyz_0_yzzzz_0,  \
                             g_0_yyz_0_yzzzz_1,  \
                             g_0_yyz_0_zzzz_1,   \
                             g_0_yyz_0_zzzzz_0,  \
                             g_0_yyz_0_zzzzz_1,  \
                             wp_x,               \
                             wp_z,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyz_0_xxxxx_0[i] = g_0_xyy_0_xxxxx_0[i] * pb_z + g_0_xyy_0_xxxxx_1[i] * wp_z[i];

        g_0_xyyz_0_xxxxy_0[i] = g_0_xyy_0_xxxxy_0[i] * pb_z + g_0_xyy_0_xxxxy_1[i] * wp_z[i];

        g_0_xyyz_0_xxxxz_0[i] = 4.0 * g_0_yyz_0_xxxz_1[i] * fi_abcd_0 + g_0_yyz_0_xxxxz_0[i] * pb_x + g_0_yyz_0_xxxxz_1[i] * wp_x[i];

        g_0_xyyz_0_xxxyy_0[i] = g_0_xyy_0_xxxyy_0[i] * pb_z + g_0_xyy_0_xxxyy_1[i] * wp_z[i];

        g_0_xyyz_0_xxxyz_0[i] = 3.0 * g_0_yyz_0_xxyz_1[i] * fi_abcd_0 + g_0_yyz_0_xxxyz_0[i] * pb_x + g_0_yyz_0_xxxyz_1[i] * wp_x[i];

        g_0_xyyz_0_xxxzz_0[i] = 3.0 * g_0_yyz_0_xxzz_1[i] * fi_abcd_0 + g_0_yyz_0_xxxzz_0[i] * pb_x + g_0_yyz_0_xxxzz_1[i] * wp_x[i];

        g_0_xyyz_0_xxyyy_0[i] = g_0_xyy_0_xxyyy_0[i] * pb_z + g_0_xyy_0_xxyyy_1[i] * wp_z[i];

        g_0_xyyz_0_xxyyz_0[i] = 2.0 * g_0_yyz_0_xyyz_1[i] * fi_abcd_0 + g_0_yyz_0_xxyyz_0[i] * pb_x + g_0_yyz_0_xxyyz_1[i] * wp_x[i];

        g_0_xyyz_0_xxyzz_0[i] = 2.0 * g_0_yyz_0_xyzz_1[i] * fi_abcd_0 + g_0_yyz_0_xxyzz_0[i] * pb_x + g_0_yyz_0_xxyzz_1[i] * wp_x[i];

        g_0_xyyz_0_xxzzz_0[i] = 2.0 * g_0_yyz_0_xzzz_1[i] * fi_abcd_0 + g_0_yyz_0_xxzzz_0[i] * pb_x + g_0_yyz_0_xxzzz_1[i] * wp_x[i];

        g_0_xyyz_0_xyyyy_0[i] = g_0_xyy_0_xyyyy_0[i] * pb_z + g_0_xyy_0_xyyyy_1[i] * wp_z[i];

        g_0_xyyz_0_xyyyz_0[i] = g_0_yyz_0_yyyz_1[i] * fi_abcd_0 + g_0_yyz_0_xyyyz_0[i] * pb_x + g_0_yyz_0_xyyyz_1[i] * wp_x[i];

        g_0_xyyz_0_xyyzz_0[i] = g_0_yyz_0_yyzz_1[i] * fi_abcd_0 + g_0_yyz_0_xyyzz_0[i] * pb_x + g_0_yyz_0_xyyzz_1[i] * wp_x[i];

        g_0_xyyz_0_xyzzz_0[i] = g_0_yyz_0_yzzz_1[i] * fi_abcd_0 + g_0_yyz_0_xyzzz_0[i] * pb_x + g_0_yyz_0_xyzzz_1[i] * wp_x[i];

        g_0_xyyz_0_xzzzz_0[i] = g_0_yyz_0_zzzz_1[i] * fi_abcd_0 + g_0_yyz_0_xzzzz_0[i] * pb_x + g_0_yyz_0_xzzzz_1[i] * wp_x[i];

        g_0_xyyz_0_yyyyy_0[i] = g_0_yyz_0_yyyyy_0[i] * pb_x + g_0_yyz_0_yyyyy_1[i] * wp_x[i];

        g_0_xyyz_0_yyyyz_0[i] = g_0_yyz_0_yyyyz_0[i] * pb_x + g_0_yyz_0_yyyyz_1[i] * wp_x[i];

        g_0_xyyz_0_yyyzz_0[i] = g_0_yyz_0_yyyzz_0[i] * pb_x + g_0_yyz_0_yyyzz_1[i] * wp_x[i];

        g_0_xyyz_0_yyzzz_0[i] = g_0_yyz_0_yyzzz_0[i] * pb_x + g_0_yyz_0_yyzzz_1[i] * wp_x[i];

        g_0_xyyz_0_yzzzz_0[i] = g_0_yyz_0_yzzzz_0[i] * pb_x + g_0_yyz_0_yzzzz_1[i] * wp_x[i];

        g_0_xyyz_0_zzzzz_0[i] = g_0_yyz_0_zzzzz_0[i] * pb_x + g_0_yyz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 168-189 components of targeted buffer : SGSH

    auto g_0_xyzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sgsh + 168);

    auto g_0_xyzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sgsh + 169);

    auto g_0_xyzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sgsh + 170);

    auto g_0_xyzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sgsh + 171);

    auto g_0_xyzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sgsh + 172);

    auto g_0_xyzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sgsh + 173);

    auto g_0_xyzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sgsh + 174);

    auto g_0_xyzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sgsh + 175);

    auto g_0_xyzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sgsh + 176);

    auto g_0_xyzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sgsh + 177);

    auto g_0_xyzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 178);

    auto g_0_xyzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 179);

    auto g_0_xyzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 180);

    auto g_0_xyzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 181);

    auto g_0_xyzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 182);

    auto g_0_xyzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 183);

    auto g_0_xyzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 184);

    auto g_0_xyzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 185);

    auto g_0_xyzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 186);

    auto g_0_xyzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 187);

    auto g_0_xyzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 188);

#pragma omp simd aligned(g_0_xyzz_0_xxxxx_0,     \
                             g_0_xyzz_0_xxxxy_0, \
                             g_0_xyzz_0_xxxxz_0, \
                             g_0_xyzz_0_xxxyy_0, \
                             g_0_xyzz_0_xxxyz_0, \
                             g_0_xyzz_0_xxxzz_0, \
                             g_0_xyzz_0_xxyyy_0, \
                             g_0_xyzz_0_xxyyz_0, \
                             g_0_xyzz_0_xxyzz_0, \
                             g_0_xyzz_0_xxzzz_0, \
                             g_0_xyzz_0_xyyyy_0, \
                             g_0_xyzz_0_xyyyz_0, \
                             g_0_xyzz_0_xyyzz_0, \
                             g_0_xyzz_0_xyzzz_0, \
                             g_0_xyzz_0_xzzzz_0, \
                             g_0_xyzz_0_yyyyy_0, \
                             g_0_xyzz_0_yyyyz_0, \
                             g_0_xyzz_0_yyyzz_0, \
                             g_0_xyzz_0_yyzzz_0, \
                             g_0_xyzz_0_yzzzz_0, \
                             g_0_xyzz_0_zzzzz_0, \
                             g_0_xzz_0_xxxxx_0,  \
                             g_0_xzz_0_xxxxx_1,  \
                             g_0_xzz_0_xxxxz_0,  \
                             g_0_xzz_0_xxxxz_1,  \
                             g_0_xzz_0_xxxzz_0,  \
                             g_0_xzz_0_xxxzz_1,  \
                             g_0_xzz_0_xxzzz_0,  \
                             g_0_xzz_0_xxzzz_1,  \
                             g_0_xzz_0_xzzzz_0,  \
                             g_0_xzz_0_xzzzz_1,  \
                             g_0_yzz_0_xxxxy_0,  \
                             g_0_yzz_0_xxxxy_1,  \
                             g_0_yzz_0_xxxy_1,   \
                             g_0_yzz_0_xxxyy_0,  \
                             g_0_yzz_0_xxxyy_1,  \
                             g_0_yzz_0_xxxyz_0,  \
                             g_0_yzz_0_xxxyz_1,  \
                             g_0_yzz_0_xxyy_1,   \
                             g_0_yzz_0_xxyyy_0,  \
                             g_0_yzz_0_xxyyy_1,  \
                             g_0_yzz_0_xxyyz_0,  \
                             g_0_yzz_0_xxyyz_1,  \
                             g_0_yzz_0_xxyz_1,   \
                             g_0_yzz_0_xxyzz_0,  \
                             g_0_yzz_0_xxyzz_1,  \
                             g_0_yzz_0_xyyy_1,   \
                             g_0_yzz_0_xyyyy_0,  \
                             g_0_yzz_0_xyyyy_1,  \
                             g_0_yzz_0_xyyyz_0,  \
                             g_0_yzz_0_xyyyz_1,  \
                             g_0_yzz_0_xyyz_1,   \
                             g_0_yzz_0_xyyzz_0,  \
                             g_0_yzz_0_xyyzz_1,  \
                             g_0_yzz_0_xyzz_1,   \
                             g_0_yzz_0_xyzzz_0,  \
                             g_0_yzz_0_xyzzz_1,  \
                             g_0_yzz_0_yyyy_1,   \
                             g_0_yzz_0_yyyyy_0,  \
                             g_0_yzz_0_yyyyy_1,  \
                             g_0_yzz_0_yyyyz_0,  \
                             g_0_yzz_0_yyyyz_1,  \
                             g_0_yzz_0_yyyz_1,   \
                             g_0_yzz_0_yyyzz_0,  \
                             g_0_yzz_0_yyyzz_1,  \
                             g_0_yzz_0_yyzz_1,   \
                             g_0_yzz_0_yyzzz_0,  \
                             g_0_yzz_0_yyzzz_1,  \
                             g_0_yzz_0_yzzz_1,   \
                             g_0_yzz_0_yzzzz_0,  \
                             g_0_yzz_0_yzzzz_1,  \
                             g_0_yzz_0_zzzzz_0,  \
                             g_0_yzz_0_zzzzz_1,  \
                             wp_x,               \
                             wp_y,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyzz_0_xxxxx_0[i] = g_0_xzz_0_xxxxx_0[i] * pb_y + g_0_xzz_0_xxxxx_1[i] * wp_y[i];

        g_0_xyzz_0_xxxxy_0[i] = 4.0 * g_0_yzz_0_xxxy_1[i] * fi_abcd_0 + g_0_yzz_0_xxxxy_0[i] * pb_x + g_0_yzz_0_xxxxy_1[i] * wp_x[i];

        g_0_xyzz_0_xxxxz_0[i] = g_0_xzz_0_xxxxz_0[i] * pb_y + g_0_xzz_0_xxxxz_1[i] * wp_y[i];

        g_0_xyzz_0_xxxyy_0[i] = 3.0 * g_0_yzz_0_xxyy_1[i] * fi_abcd_0 + g_0_yzz_0_xxxyy_0[i] * pb_x + g_0_yzz_0_xxxyy_1[i] * wp_x[i];

        g_0_xyzz_0_xxxyz_0[i] = 3.0 * g_0_yzz_0_xxyz_1[i] * fi_abcd_0 + g_0_yzz_0_xxxyz_0[i] * pb_x + g_0_yzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xyzz_0_xxxzz_0[i] = g_0_xzz_0_xxxzz_0[i] * pb_y + g_0_xzz_0_xxxzz_1[i] * wp_y[i];

        g_0_xyzz_0_xxyyy_0[i] = 2.0 * g_0_yzz_0_xyyy_1[i] * fi_abcd_0 + g_0_yzz_0_xxyyy_0[i] * pb_x + g_0_yzz_0_xxyyy_1[i] * wp_x[i];

        g_0_xyzz_0_xxyyz_0[i] = 2.0 * g_0_yzz_0_xyyz_1[i] * fi_abcd_0 + g_0_yzz_0_xxyyz_0[i] * pb_x + g_0_yzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xyzz_0_xxyzz_0[i] = 2.0 * g_0_yzz_0_xyzz_1[i] * fi_abcd_0 + g_0_yzz_0_xxyzz_0[i] * pb_x + g_0_yzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xyzz_0_xxzzz_0[i] = g_0_xzz_0_xxzzz_0[i] * pb_y + g_0_xzz_0_xxzzz_1[i] * wp_y[i];

        g_0_xyzz_0_xyyyy_0[i] = g_0_yzz_0_yyyy_1[i] * fi_abcd_0 + g_0_yzz_0_xyyyy_0[i] * pb_x + g_0_yzz_0_xyyyy_1[i] * wp_x[i];

        g_0_xyzz_0_xyyyz_0[i] = g_0_yzz_0_yyyz_1[i] * fi_abcd_0 + g_0_yzz_0_xyyyz_0[i] * pb_x + g_0_yzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xyzz_0_xyyzz_0[i] = g_0_yzz_0_yyzz_1[i] * fi_abcd_0 + g_0_yzz_0_xyyzz_0[i] * pb_x + g_0_yzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xyzz_0_xyzzz_0[i] = g_0_yzz_0_yzzz_1[i] * fi_abcd_0 + g_0_yzz_0_xyzzz_0[i] * pb_x + g_0_yzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xyzz_0_xzzzz_0[i] = g_0_xzz_0_xzzzz_0[i] * pb_y + g_0_xzz_0_xzzzz_1[i] * wp_y[i];

        g_0_xyzz_0_yyyyy_0[i] = g_0_yzz_0_yyyyy_0[i] * pb_x + g_0_yzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xyzz_0_yyyyz_0[i] = g_0_yzz_0_yyyyz_0[i] * pb_x + g_0_yzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xyzz_0_yyyzz_0[i] = g_0_yzz_0_yyyzz_0[i] * pb_x + g_0_yzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xyzz_0_yyzzz_0[i] = g_0_yzz_0_yyzzz_0[i] * pb_x + g_0_yzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xyzz_0_yzzzz_0[i] = g_0_yzz_0_yzzzz_0[i] * pb_x + g_0_yzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xyzz_0_zzzzz_0[i] = g_0_yzz_0_zzzzz_0[i] * pb_x + g_0_yzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 189-210 components of targeted buffer : SGSH

    auto g_0_xzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sgsh + 189);

    auto g_0_xzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sgsh + 190);

    auto g_0_xzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sgsh + 191);

    auto g_0_xzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sgsh + 192);

    auto g_0_xzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sgsh + 193);

    auto g_0_xzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sgsh + 194);

    auto g_0_xzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sgsh + 195);

    auto g_0_xzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sgsh + 196);

    auto g_0_xzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sgsh + 197);

    auto g_0_xzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sgsh + 198);

    auto g_0_xzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 199);

    auto g_0_xzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 200);

    auto g_0_xzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 201);

    auto g_0_xzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 202);

    auto g_0_xzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 203);

    auto g_0_xzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 204);

    auto g_0_xzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 205);

    auto g_0_xzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 206);

    auto g_0_xzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 207);

    auto g_0_xzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 208);

    auto g_0_xzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 209);

#pragma omp simd aligned(g_0_xzzz_0_xxxxx_0,     \
                             g_0_xzzz_0_xxxxy_0, \
                             g_0_xzzz_0_xxxxz_0, \
                             g_0_xzzz_0_xxxyy_0, \
                             g_0_xzzz_0_xxxyz_0, \
                             g_0_xzzz_0_xxxzz_0, \
                             g_0_xzzz_0_xxyyy_0, \
                             g_0_xzzz_0_xxyyz_0, \
                             g_0_xzzz_0_xxyzz_0, \
                             g_0_xzzz_0_xxzzz_0, \
                             g_0_xzzz_0_xyyyy_0, \
                             g_0_xzzz_0_xyyyz_0, \
                             g_0_xzzz_0_xyyzz_0, \
                             g_0_xzzz_0_xyzzz_0, \
                             g_0_xzzz_0_xzzzz_0, \
                             g_0_xzzz_0_yyyyy_0, \
                             g_0_xzzz_0_yyyyz_0, \
                             g_0_xzzz_0_yyyzz_0, \
                             g_0_xzzz_0_yyzzz_0, \
                             g_0_xzzz_0_yzzzz_0, \
                             g_0_xzzz_0_zzzzz_0, \
                             g_0_zzz_0_xxxx_1,   \
                             g_0_zzz_0_xxxxx_0,  \
                             g_0_zzz_0_xxxxx_1,  \
                             g_0_zzz_0_xxxxy_0,  \
                             g_0_zzz_0_xxxxy_1,  \
                             g_0_zzz_0_xxxxz_0,  \
                             g_0_zzz_0_xxxxz_1,  \
                             g_0_zzz_0_xxxy_1,   \
                             g_0_zzz_0_xxxyy_0,  \
                             g_0_zzz_0_xxxyy_1,  \
                             g_0_zzz_0_xxxyz_0,  \
                             g_0_zzz_0_xxxyz_1,  \
                             g_0_zzz_0_xxxz_1,   \
                             g_0_zzz_0_xxxzz_0,  \
                             g_0_zzz_0_xxxzz_1,  \
                             g_0_zzz_0_xxyy_1,   \
                             g_0_zzz_0_xxyyy_0,  \
                             g_0_zzz_0_xxyyy_1,  \
                             g_0_zzz_0_xxyyz_0,  \
                             g_0_zzz_0_xxyyz_1,  \
                             g_0_zzz_0_xxyz_1,   \
                             g_0_zzz_0_xxyzz_0,  \
                             g_0_zzz_0_xxyzz_1,  \
                             g_0_zzz_0_xxzz_1,   \
                             g_0_zzz_0_xxzzz_0,  \
                             g_0_zzz_0_xxzzz_1,  \
                             g_0_zzz_0_xyyy_1,   \
                             g_0_zzz_0_xyyyy_0,  \
                             g_0_zzz_0_xyyyy_1,  \
                             g_0_zzz_0_xyyyz_0,  \
                             g_0_zzz_0_xyyyz_1,  \
                             g_0_zzz_0_xyyz_1,   \
                             g_0_zzz_0_xyyzz_0,  \
                             g_0_zzz_0_xyyzz_1,  \
                             g_0_zzz_0_xyzz_1,   \
                             g_0_zzz_0_xyzzz_0,  \
                             g_0_zzz_0_xyzzz_1,  \
                             g_0_zzz_0_xzzz_1,   \
                             g_0_zzz_0_xzzzz_0,  \
                             g_0_zzz_0_xzzzz_1,  \
                             g_0_zzz_0_yyyy_1,   \
                             g_0_zzz_0_yyyyy_0,  \
                             g_0_zzz_0_yyyyy_1,  \
                             g_0_zzz_0_yyyyz_0,  \
                             g_0_zzz_0_yyyyz_1,  \
                             g_0_zzz_0_yyyz_1,   \
                             g_0_zzz_0_yyyzz_0,  \
                             g_0_zzz_0_yyyzz_1,  \
                             g_0_zzz_0_yyzz_1,   \
                             g_0_zzz_0_yyzzz_0,  \
                             g_0_zzz_0_yyzzz_1,  \
                             g_0_zzz_0_yzzz_1,   \
                             g_0_zzz_0_yzzzz_0,  \
                             g_0_zzz_0_yzzzz_1,  \
                             g_0_zzz_0_zzzz_1,   \
                             g_0_zzz_0_zzzzz_0,  \
                             g_0_zzz_0_zzzzz_1,  \
                             wp_x,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzzz_0_xxxxx_0[i] = 5.0 * g_0_zzz_0_xxxx_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxx_0[i] * pb_x + g_0_zzz_0_xxxxx_1[i] * wp_x[i];

        g_0_xzzz_0_xxxxy_0[i] = 4.0 * g_0_zzz_0_xxxy_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxy_0[i] * pb_x + g_0_zzz_0_xxxxy_1[i] * wp_x[i];

        g_0_xzzz_0_xxxxz_0[i] = 4.0 * g_0_zzz_0_xxxz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxz_0[i] * pb_x + g_0_zzz_0_xxxxz_1[i] * wp_x[i];

        g_0_xzzz_0_xxxyy_0[i] = 3.0 * g_0_zzz_0_xxyy_1[i] * fi_abcd_0 + g_0_zzz_0_xxxyy_0[i] * pb_x + g_0_zzz_0_xxxyy_1[i] * wp_x[i];

        g_0_xzzz_0_xxxyz_0[i] = 3.0 * g_0_zzz_0_xxyz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxyz_0[i] * pb_x + g_0_zzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xzzz_0_xxxzz_0[i] = 3.0 * g_0_zzz_0_xxzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxzz_0[i] * pb_x + g_0_zzz_0_xxxzz_1[i] * wp_x[i];

        g_0_xzzz_0_xxyyy_0[i] = 2.0 * g_0_zzz_0_xyyy_1[i] * fi_abcd_0 + g_0_zzz_0_xxyyy_0[i] * pb_x + g_0_zzz_0_xxyyy_1[i] * wp_x[i];

        g_0_xzzz_0_xxyyz_0[i] = 2.0 * g_0_zzz_0_xyyz_1[i] * fi_abcd_0 + g_0_zzz_0_xxyyz_0[i] * pb_x + g_0_zzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xzzz_0_xxyzz_0[i] = 2.0 * g_0_zzz_0_xyzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxyzz_0[i] * pb_x + g_0_zzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xzzz_0_xxzzz_0[i] = 2.0 * g_0_zzz_0_xzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxzzz_0[i] * pb_x + g_0_zzz_0_xxzzz_1[i] * wp_x[i];

        g_0_xzzz_0_xyyyy_0[i] = g_0_zzz_0_yyyy_1[i] * fi_abcd_0 + g_0_zzz_0_xyyyy_0[i] * pb_x + g_0_zzz_0_xyyyy_1[i] * wp_x[i];

        g_0_xzzz_0_xyyyz_0[i] = g_0_zzz_0_yyyz_1[i] * fi_abcd_0 + g_0_zzz_0_xyyyz_0[i] * pb_x + g_0_zzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xzzz_0_xyyzz_0[i] = g_0_zzz_0_yyzz_1[i] * fi_abcd_0 + g_0_zzz_0_xyyzz_0[i] * pb_x + g_0_zzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xzzz_0_xyzzz_0[i] = g_0_zzz_0_yzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xyzzz_0[i] * pb_x + g_0_zzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xzzz_0_xzzzz_0[i] = g_0_zzz_0_zzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xzzzz_0[i] * pb_x + g_0_zzz_0_xzzzz_1[i] * wp_x[i];

        g_0_xzzz_0_yyyyy_0[i] = g_0_zzz_0_yyyyy_0[i] * pb_x + g_0_zzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xzzz_0_yyyyz_0[i] = g_0_zzz_0_yyyyz_0[i] * pb_x + g_0_zzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xzzz_0_yyyzz_0[i] = g_0_zzz_0_yyyzz_0[i] * pb_x + g_0_zzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xzzz_0_yyzzz_0[i] = g_0_zzz_0_yyzzz_0[i] * pb_x + g_0_zzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xzzz_0_yzzzz_0[i] = g_0_zzz_0_yzzzz_0[i] * pb_x + g_0_zzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xzzz_0_zzzzz_0[i] = g_0_zzz_0_zzzzz_0[i] * pb_x + g_0_zzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 210-231 components of targeted buffer : SGSH

    auto g_0_yyyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_sgsh + 210);

    auto g_0_yyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_sgsh + 211);

    auto g_0_yyyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_sgsh + 212);

    auto g_0_yyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_sgsh + 213);

    auto g_0_yyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_sgsh + 214);

    auto g_0_yyyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_sgsh + 215);

    auto g_0_yyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_sgsh + 216);

    auto g_0_yyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_sgsh + 217);

    auto g_0_yyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_sgsh + 218);

    auto g_0_yyyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_sgsh + 219);

    auto g_0_yyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 220);

    auto g_0_yyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 221);

    auto g_0_yyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 222);

    auto g_0_yyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 223);

    auto g_0_yyyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 224);

    auto g_0_yyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 225);

    auto g_0_yyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 226);

    auto g_0_yyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 227);

    auto g_0_yyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 228);

    auto g_0_yyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 229);

    auto g_0_yyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 230);

#pragma omp simd aligned(g_0_yy_0_xxxxx_0,       \
                             g_0_yy_0_xxxxx_1,   \
                             g_0_yy_0_xxxxy_0,   \
                             g_0_yy_0_xxxxy_1,   \
                             g_0_yy_0_xxxxz_0,   \
                             g_0_yy_0_xxxxz_1,   \
                             g_0_yy_0_xxxyy_0,   \
                             g_0_yy_0_xxxyy_1,   \
                             g_0_yy_0_xxxyz_0,   \
                             g_0_yy_0_xxxyz_1,   \
                             g_0_yy_0_xxxzz_0,   \
                             g_0_yy_0_xxxzz_1,   \
                             g_0_yy_0_xxyyy_0,   \
                             g_0_yy_0_xxyyy_1,   \
                             g_0_yy_0_xxyyz_0,   \
                             g_0_yy_0_xxyyz_1,   \
                             g_0_yy_0_xxyzz_0,   \
                             g_0_yy_0_xxyzz_1,   \
                             g_0_yy_0_xxzzz_0,   \
                             g_0_yy_0_xxzzz_1,   \
                             g_0_yy_0_xyyyy_0,   \
                             g_0_yy_0_xyyyy_1,   \
                             g_0_yy_0_xyyyz_0,   \
                             g_0_yy_0_xyyyz_1,   \
                             g_0_yy_0_xyyzz_0,   \
                             g_0_yy_0_xyyzz_1,   \
                             g_0_yy_0_xyzzz_0,   \
                             g_0_yy_0_xyzzz_1,   \
                             g_0_yy_0_xzzzz_0,   \
                             g_0_yy_0_xzzzz_1,   \
                             g_0_yy_0_yyyyy_0,   \
                             g_0_yy_0_yyyyy_1,   \
                             g_0_yy_0_yyyyz_0,   \
                             g_0_yy_0_yyyyz_1,   \
                             g_0_yy_0_yyyzz_0,   \
                             g_0_yy_0_yyyzz_1,   \
                             g_0_yy_0_yyzzz_0,   \
                             g_0_yy_0_yyzzz_1,   \
                             g_0_yy_0_yzzzz_0,   \
                             g_0_yy_0_yzzzz_1,   \
                             g_0_yy_0_zzzzz_0,   \
                             g_0_yy_0_zzzzz_1,   \
                             g_0_yyy_0_xxxx_1,   \
                             g_0_yyy_0_xxxxx_0,  \
                             g_0_yyy_0_xxxxx_1,  \
                             g_0_yyy_0_xxxxy_0,  \
                             g_0_yyy_0_xxxxy_1,  \
                             g_0_yyy_0_xxxxz_0,  \
                             g_0_yyy_0_xxxxz_1,  \
                             g_0_yyy_0_xxxy_1,   \
                             g_0_yyy_0_xxxyy_0,  \
                             g_0_yyy_0_xxxyy_1,  \
                             g_0_yyy_0_xxxyz_0,  \
                             g_0_yyy_0_xxxyz_1,  \
                             g_0_yyy_0_xxxz_1,   \
                             g_0_yyy_0_xxxzz_0,  \
                             g_0_yyy_0_xxxzz_1,  \
                             g_0_yyy_0_xxyy_1,   \
                             g_0_yyy_0_xxyyy_0,  \
                             g_0_yyy_0_xxyyy_1,  \
                             g_0_yyy_0_xxyyz_0,  \
                             g_0_yyy_0_xxyyz_1,  \
                             g_0_yyy_0_xxyz_1,   \
                             g_0_yyy_0_xxyzz_0,  \
                             g_0_yyy_0_xxyzz_1,  \
                             g_0_yyy_0_xxzz_1,   \
                             g_0_yyy_0_xxzzz_0,  \
                             g_0_yyy_0_xxzzz_1,  \
                             g_0_yyy_0_xyyy_1,   \
                             g_0_yyy_0_xyyyy_0,  \
                             g_0_yyy_0_xyyyy_1,  \
                             g_0_yyy_0_xyyyz_0,  \
                             g_0_yyy_0_xyyyz_1,  \
                             g_0_yyy_0_xyyz_1,   \
                             g_0_yyy_0_xyyzz_0,  \
                             g_0_yyy_0_xyyzz_1,  \
                             g_0_yyy_0_xyzz_1,   \
                             g_0_yyy_0_xyzzz_0,  \
                             g_0_yyy_0_xyzzz_1,  \
                             g_0_yyy_0_xzzz_1,   \
                             g_0_yyy_0_xzzzz_0,  \
                             g_0_yyy_0_xzzzz_1,  \
                             g_0_yyy_0_yyyy_1,   \
                             g_0_yyy_0_yyyyy_0,  \
                             g_0_yyy_0_yyyyy_1,  \
                             g_0_yyy_0_yyyyz_0,  \
                             g_0_yyy_0_yyyyz_1,  \
                             g_0_yyy_0_yyyz_1,   \
                             g_0_yyy_0_yyyzz_0,  \
                             g_0_yyy_0_yyyzz_1,  \
                             g_0_yyy_0_yyzz_1,   \
                             g_0_yyy_0_yyzzz_0,  \
                             g_0_yyy_0_yyzzz_1,  \
                             g_0_yyy_0_yzzz_1,   \
                             g_0_yyy_0_yzzzz_0,  \
                             g_0_yyy_0_yzzzz_1,  \
                             g_0_yyy_0_zzzz_1,   \
                             g_0_yyy_0_zzzzz_0,  \
                             g_0_yyy_0_zzzzz_1,  \
                             g_0_yyyy_0_xxxxx_0, \
                             g_0_yyyy_0_xxxxy_0, \
                             g_0_yyyy_0_xxxxz_0, \
                             g_0_yyyy_0_xxxyy_0, \
                             g_0_yyyy_0_xxxyz_0, \
                             g_0_yyyy_0_xxxzz_0, \
                             g_0_yyyy_0_xxyyy_0, \
                             g_0_yyyy_0_xxyyz_0, \
                             g_0_yyyy_0_xxyzz_0, \
                             g_0_yyyy_0_xxzzz_0, \
                             g_0_yyyy_0_xyyyy_0, \
                             g_0_yyyy_0_xyyyz_0, \
                             g_0_yyyy_0_xyyzz_0, \
                             g_0_yyyy_0_xyzzz_0, \
                             g_0_yyyy_0_xzzzz_0, \
                             g_0_yyyy_0_yyyyy_0, \
                             g_0_yyyy_0_yyyyz_0, \
                             g_0_yyyy_0_yyyzz_0, \
                             g_0_yyyy_0_yyzzz_0, \
                             g_0_yyyy_0_yzzzz_0, \
                             g_0_yyyy_0_zzzzz_0, \
                             wp_y,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyy_0_xxxxx_0[i] =
            3.0 * g_0_yy_0_xxxxx_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxxx_1[i] * fti_ab_0 + g_0_yyy_0_xxxxx_0[i] * pb_y + g_0_yyy_0_xxxxx_1[i] * wp_y[i];

        g_0_yyyy_0_xxxxy_0[i] = 3.0 * g_0_yy_0_xxxxy_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxxy_1[i] * fti_ab_0 + g_0_yyy_0_xxxx_1[i] * fi_abcd_0 +
                                g_0_yyy_0_xxxxy_0[i] * pb_y + g_0_yyy_0_xxxxy_1[i] * wp_y[i];

        g_0_yyyy_0_xxxxz_0[i] =
            3.0 * g_0_yy_0_xxxxz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxxz_1[i] * fti_ab_0 + g_0_yyy_0_xxxxz_0[i] * pb_y + g_0_yyy_0_xxxxz_1[i] * wp_y[i];

        g_0_yyyy_0_xxxyy_0[i] = 3.0 * g_0_yy_0_xxxyy_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxyy_1[i] * fti_ab_0 + 2.0 * g_0_yyy_0_xxxy_1[i] * fi_abcd_0 +
                                g_0_yyy_0_xxxyy_0[i] * pb_y + g_0_yyy_0_xxxyy_1[i] * wp_y[i];

        g_0_yyyy_0_xxxyz_0[i] = 3.0 * g_0_yy_0_xxxyz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxyz_1[i] * fti_ab_0 + g_0_yyy_0_xxxz_1[i] * fi_abcd_0 +
                                g_0_yyy_0_xxxyz_0[i] * pb_y + g_0_yyy_0_xxxyz_1[i] * wp_y[i];

        g_0_yyyy_0_xxxzz_0[i] =
            3.0 * g_0_yy_0_xxxzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxzz_1[i] * fti_ab_0 + g_0_yyy_0_xxxzz_0[i] * pb_y + g_0_yyy_0_xxxzz_1[i] * wp_y[i];

        g_0_yyyy_0_xxyyy_0[i] = 3.0 * g_0_yy_0_xxyyy_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxyyy_1[i] * fti_ab_0 + 3.0 * g_0_yyy_0_xxyy_1[i] * fi_abcd_0 +
                                g_0_yyy_0_xxyyy_0[i] * pb_y + g_0_yyy_0_xxyyy_1[i] * wp_y[i];

        g_0_yyyy_0_xxyyz_0[i] = 3.0 * g_0_yy_0_xxyyz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yyy_0_xxyz_1[i] * fi_abcd_0 +
                                g_0_yyy_0_xxyyz_0[i] * pb_y + g_0_yyy_0_xxyyz_1[i] * wp_y[i];

        g_0_yyyy_0_xxyzz_0[i] = 3.0 * g_0_yy_0_xxyzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxyzz_1[i] * fti_ab_0 + g_0_yyy_0_xxzz_1[i] * fi_abcd_0 +
                                g_0_yyy_0_xxyzz_0[i] * pb_y + g_0_yyy_0_xxyzz_1[i] * wp_y[i];

        g_0_yyyy_0_xxzzz_0[i] =
            3.0 * g_0_yy_0_xxzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxzzz_1[i] * fti_ab_0 + g_0_yyy_0_xxzzz_0[i] * pb_y + g_0_yyy_0_xxzzz_1[i] * wp_y[i];

        g_0_yyyy_0_xyyyy_0[i] = 3.0 * g_0_yy_0_xyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xyyyy_1[i] * fti_ab_0 + 4.0 * g_0_yyy_0_xyyy_1[i] * fi_abcd_0 +
                                g_0_yyy_0_xyyyy_0[i] * pb_y + g_0_yyy_0_xyyyy_1[i] * wp_y[i];

        g_0_yyyy_0_xyyyz_0[i] = 3.0 * g_0_yy_0_xyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yyy_0_xyyz_1[i] * fi_abcd_0 +
                                g_0_yyy_0_xyyyz_0[i] * pb_y + g_0_yyy_0_xyyyz_1[i] * wp_y[i];

        g_0_yyyy_0_xyyzz_0[i] = 3.0 * g_0_yy_0_xyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yyy_0_xyzz_1[i] * fi_abcd_0 +
                                g_0_yyy_0_xyyzz_0[i] * pb_y + g_0_yyy_0_xyyzz_1[i] * wp_y[i];

        g_0_yyyy_0_xyzzz_0[i] = 3.0 * g_0_yy_0_xyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xyzzz_1[i] * fti_ab_0 + g_0_yyy_0_xzzz_1[i] * fi_abcd_0 +
                                g_0_yyy_0_xyzzz_0[i] * pb_y + g_0_yyy_0_xyzzz_1[i] * wp_y[i];

        g_0_yyyy_0_xzzzz_0[i] =
            3.0 * g_0_yy_0_xzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xzzzz_1[i] * fti_ab_0 + g_0_yyy_0_xzzzz_0[i] * pb_y + g_0_yyy_0_xzzzz_1[i] * wp_y[i];

        g_0_yyyy_0_yyyyy_0[i] = 3.0 * g_0_yy_0_yyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_yyyyy_1[i] * fti_ab_0 + 5.0 * g_0_yyy_0_yyyy_1[i] * fi_abcd_0 +
                                g_0_yyy_0_yyyyy_0[i] * pb_y + g_0_yyy_0_yyyyy_1[i] * wp_y[i];

        g_0_yyyy_0_yyyyz_0[i] = 3.0 * g_0_yy_0_yyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_yyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yyy_0_yyyz_1[i] * fi_abcd_0 +
                                g_0_yyy_0_yyyyz_0[i] * pb_y + g_0_yyy_0_yyyyz_1[i] * wp_y[i];

        g_0_yyyy_0_yyyzz_0[i] = 3.0 * g_0_yy_0_yyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_yyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yyy_0_yyzz_1[i] * fi_abcd_0 +
                                g_0_yyy_0_yyyzz_0[i] * pb_y + g_0_yyy_0_yyyzz_1[i] * wp_y[i];

        g_0_yyyy_0_yyzzz_0[i] = 3.0 * g_0_yy_0_yyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_yyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyy_0_yzzz_1[i] * fi_abcd_0 +
                                g_0_yyy_0_yyzzz_0[i] * pb_y + g_0_yyy_0_yyzzz_1[i] * wp_y[i];

        g_0_yyyy_0_yzzzz_0[i] = 3.0 * g_0_yy_0_yzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_yzzzz_1[i] * fti_ab_0 + g_0_yyy_0_zzzz_1[i] * fi_abcd_0 +
                                g_0_yyy_0_yzzzz_0[i] * pb_y + g_0_yyy_0_yzzzz_1[i] * wp_y[i];

        g_0_yyyy_0_zzzzz_0[i] =
            3.0 * g_0_yy_0_zzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_zzzzz_1[i] * fti_ab_0 + g_0_yyy_0_zzzzz_0[i] * pb_y + g_0_yyy_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 231-252 components of targeted buffer : SGSH

    auto g_0_yyyz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sgsh + 231);

    auto g_0_yyyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sgsh + 232);

    auto g_0_yyyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sgsh + 233);

    auto g_0_yyyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sgsh + 234);

    auto g_0_yyyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sgsh + 235);

    auto g_0_yyyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sgsh + 236);

    auto g_0_yyyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sgsh + 237);

    auto g_0_yyyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sgsh + 238);

    auto g_0_yyyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sgsh + 239);

    auto g_0_yyyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sgsh + 240);

    auto g_0_yyyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 241);

    auto g_0_yyyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 242);

    auto g_0_yyyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 243);

    auto g_0_yyyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 244);

    auto g_0_yyyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 245);

    auto g_0_yyyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 246);

    auto g_0_yyyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 247);

    auto g_0_yyyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 248);

    auto g_0_yyyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 249);

    auto g_0_yyyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 250);

    auto g_0_yyyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 251);

#pragma omp simd aligned(g_0_yyy_0_xxxx_1,       \
                             g_0_yyy_0_xxxxx_0,  \
                             g_0_yyy_0_xxxxx_1,  \
                             g_0_yyy_0_xxxxy_0,  \
                             g_0_yyy_0_xxxxy_1,  \
                             g_0_yyy_0_xxxxz_0,  \
                             g_0_yyy_0_xxxxz_1,  \
                             g_0_yyy_0_xxxy_1,   \
                             g_0_yyy_0_xxxyy_0,  \
                             g_0_yyy_0_xxxyy_1,  \
                             g_0_yyy_0_xxxyz_0,  \
                             g_0_yyy_0_xxxyz_1,  \
                             g_0_yyy_0_xxxz_1,   \
                             g_0_yyy_0_xxxzz_0,  \
                             g_0_yyy_0_xxxzz_1,  \
                             g_0_yyy_0_xxyy_1,   \
                             g_0_yyy_0_xxyyy_0,  \
                             g_0_yyy_0_xxyyy_1,  \
                             g_0_yyy_0_xxyyz_0,  \
                             g_0_yyy_0_xxyyz_1,  \
                             g_0_yyy_0_xxyz_1,   \
                             g_0_yyy_0_xxyzz_0,  \
                             g_0_yyy_0_xxyzz_1,  \
                             g_0_yyy_0_xxzz_1,   \
                             g_0_yyy_0_xxzzz_0,  \
                             g_0_yyy_0_xxzzz_1,  \
                             g_0_yyy_0_xyyy_1,   \
                             g_0_yyy_0_xyyyy_0,  \
                             g_0_yyy_0_xyyyy_1,  \
                             g_0_yyy_0_xyyyz_0,  \
                             g_0_yyy_0_xyyyz_1,  \
                             g_0_yyy_0_xyyz_1,   \
                             g_0_yyy_0_xyyzz_0,  \
                             g_0_yyy_0_xyyzz_1,  \
                             g_0_yyy_0_xyzz_1,   \
                             g_0_yyy_0_xyzzz_0,  \
                             g_0_yyy_0_xyzzz_1,  \
                             g_0_yyy_0_xzzz_1,   \
                             g_0_yyy_0_xzzzz_0,  \
                             g_0_yyy_0_xzzzz_1,  \
                             g_0_yyy_0_yyyy_1,   \
                             g_0_yyy_0_yyyyy_0,  \
                             g_0_yyy_0_yyyyy_1,  \
                             g_0_yyy_0_yyyyz_0,  \
                             g_0_yyy_0_yyyyz_1,  \
                             g_0_yyy_0_yyyz_1,   \
                             g_0_yyy_0_yyyzz_0,  \
                             g_0_yyy_0_yyyzz_1,  \
                             g_0_yyy_0_yyzz_1,   \
                             g_0_yyy_0_yyzzz_0,  \
                             g_0_yyy_0_yyzzz_1,  \
                             g_0_yyy_0_yzzz_1,   \
                             g_0_yyy_0_yzzzz_0,  \
                             g_0_yyy_0_yzzzz_1,  \
                             g_0_yyy_0_zzzz_1,   \
                             g_0_yyy_0_zzzzz_0,  \
                             g_0_yyy_0_zzzzz_1,  \
                             g_0_yyyz_0_xxxxx_0, \
                             g_0_yyyz_0_xxxxy_0, \
                             g_0_yyyz_0_xxxxz_0, \
                             g_0_yyyz_0_xxxyy_0, \
                             g_0_yyyz_0_xxxyz_0, \
                             g_0_yyyz_0_xxxzz_0, \
                             g_0_yyyz_0_xxyyy_0, \
                             g_0_yyyz_0_xxyyz_0, \
                             g_0_yyyz_0_xxyzz_0, \
                             g_0_yyyz_0_xxzzz_0, \
                             g_0_yyyz_0_xyyyy_0, \
                             g_0_yyyz_0_xyyyz_0, \
                             g_0_yyyz_0_xyyzz_0, \
                             g_0_yyyz_0_xyzzz_0, \
                             g_0_yyyz_0_xzzzz_0, \
                             g_0_yyyz_0_yyyyy_0, \
                             g_0_yyyz_0_yyyyz_0, \
                             g_0_yyyz_0_yyyzz_0, \
                             g_0_yyyz_0_yyzzz_0, \
                             g_0_yyyz_0_yzzzz_0, \
                             g_0_yyyz_0_zzzzz_0, \
                             wp_z,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyz_0_xxxxx_0[i] = g_0_yyy_0_xxxxx_0[i] * pb_z + g_0_yyy_0_xxxxx_1[i] * wp_z[i];

        g_0_yyyz_0_xxxxy_0[i] = g_0_yyy_0_xxxxy_0[i] * pb_z + g_0_yyy_0_xxxxy_1[i] * wp_z[i];

        g_0_yyyz_0_xxxxz_0[i] = g_0_yyy_0_xxxx_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxz_0[i] * pb_z + g_0_yyy_0_xxxxz_1[i] * wp_z[i];

        g_0_yyyz_0_xxxyy_0[i] = g_0_yyy_0_xxxyy_0[i] * pb_z + g_0_yyy_0_xxxyy_1[i] * wp_z[i];

        g_0_yyyz_0_xxxyz_0[i] = g_0_yyy_0_xxxy_1[i] * fi_abcd_0 + g_0_yyy_0_xxxyz_0[i] * pb_z + g_0_yyy_0_xxxyz_1[i] * wp_z[i];

        g_0_yyyz_0_xxxzz_0[i] = 2.0 * g_0_yyy_0_xxxz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxzz_0[i] * pb_z + g_0_yyy_0_xxxzz_1[i] * wp_z[i];

        g_0_yyyz_0_xxyyy_0[i] = g_0_yyy_0_xxyyy_0[i] * pb_z + g_0_yyy_0_xxyyy_1[i] * wp_z[i];

        g_0_yyyz_0_xxyyz_0[i] = g_0_yyy_0_xxyy_1[i] * fi_abcd_0 + g_0_yyy_0_xxyyz_0[i] * pb_z + g_0_yyy_0_xxyyz_1[i] * wp_z[i];

        g_0_yyyz_0_xxyzz_0[i] = 2.0 * g_0_yyy_0_xxyz_1[i] * fi_abcd_0 + g_0_yyy_0_xxyzz_0[i] * pb_z + g_0_yyy_0_xxyzz_1[i] * wp_z[i];

        g_0_yyyz_0_xxzzz_0[i] = 3.0 * g_0_yyy_0_xxzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxzzz_0[i] * pb_z + g_0_yyy_0_xxzzz_1[i] * wp_z[i];

        g_0_yyyz_0_xyyyy_0[i] = g_0_yyy_0_xyyyy_0[i] * pb_z + g_0_yyy_0_xyyyy_1[i] * wp_z[i];

        g_0_yyyz_0_xyyyz_0[i] = g_0_yyy_0_xyyy_1[i] * fi_abcd_0 + g_0_yyy_0_xyyyz_0[i] * pb_z + g_0_yyy_0_xyyyz_1[i] * wp_z[i];

        g_0_yyyz_0_xyyzz_0[i] = 2.0 * g_0_yyy_0_xyyz_1[i] * fi_abcd_0 + g_0_yyy_0_xyyzz_0[i] * pb_z + g_0_yyy_0_xyyzz_1[i] * wp_z[i];

        g_0_yyyz_0_xyzzz_0[i] = 3.0 * g_0_yyy_0_xyzz_1[i] * fi_abcd_0 + g_0_yyy_0_xyzzz_0[i] * pb_z + g_0_yyy_0_xyzzz_1[i] * wp_z[i];

        g_0_yyyz_0_xzzzz_0[i] = 4.0 * g_0_yyy_0_xzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xzzzz_0[i] * pb_z + g_0_yyy_0_xzzzz_1[i] * wp_z[i];

        g_0_yyyz_0_yyyyy_0[i] = g_0_yyy_0_yyyyy_0[i] * pb_z + g_0_yyy_0_yyyyy_1[i] * wp_z[i];

        g_0_yyyz_0_yyyyz_0[i] = g_0_yyy_0_yyyy_1[i] * fi_abcd_0 + g_0_yyy_0_yyyyz_0[i] * pb_z + g_0_yyy_0_yyyyz_1[i] * wp_z[i];

        g_0_yyyz_0_yyyzz_0[i] = 2.0 * g_0_yyy_0_yyyz_1[i] * fi_abcd_0 + g_0_yyy_0_yyyzz_0[i] * pb_z + g_0_yyy_0_yyyzz_1[i] * wp_z[i];

        g_0_yyyz_0_yyzzz_0[i] = 3.0 * g_0_yyy_0_yyzz_1[i] * fi_abcd_0 + g_0_yyy_0_yyzzz_0[i] * pb_z + g_0_yyy_0_yyzzz_1[i] * wp_z[i];

        g_0_yyyz_0_yzzzz_0[i] = 4.0 * g_0_yyy_0_yzzz_1[i] * fi_abcd_0 + g_0_yyy_0_yzzzz_0[i] * pb_z + g_0_yyy_0_yzzzz_1[i] * wp_z[i];

        g_0_yyyz_0_zzzzz_0[i] = 5.0 * g_0_yyy_0_zzzz_1[i] * fi_abcd_0 + g_0_yyy_0_zzzzz_0[i] * pb_z + g_0_yyy_0_zzzzz_1[i] * wp_z[i];
    }

    /// Set up 252-273 components of targeted buffer : SGSH

    auto g_0_yyzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sgsh + 252);

    auto g_0_yyzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sgsh + 253);

    auto g_0_yyzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sgsh + 254);

    auto g_0_yyzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sgsh + 255);

    auto g_0_yyzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sgsh + 256);

    auto g_0_yyzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sgsh + 257);

    auto g_0_yyzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sgsh + 258);

    auto g_0_yyzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sgsh + 259);

    auto g_0_yyzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sgsh + 260);

    auto g_0_yyzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sgsh + 261);

    auto g_0_yyzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 262);

    auto g_0_yyzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 263);

    auto g_0_yyzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 264);

    auto g_0_yyzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 265);

    auto g_0_yyzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 266);

    auto g_0_yyzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 267);

    auto g_0_yyzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 268);

    auto g_0_yyzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 269);

    auto g_0_yyzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 270);

    auto g_0_yyzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 271);

    auto g_0_yyzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 272);

#pragma omp simd aligned(g_0_yy_0_xxxxy_0,       \
                             g_0_yy_0_xxxxy_1,   \
                             g_0_yy_0_xxxyy_0,   \
                             g_0_yy_0_xxxyy_1,   \
                             g_0_yy_0_xxyyy_0,   \
                             g_0_yy_0_xxyyy_1,   \
                             g_0_yy_0_xyyyy_0,   \
                             g_0_yy_0_xyyyy_1,   \
                             g_0_yy_0_yyyyy_0,   \
                             g_0_yy_0_yyyyy_1,   \
                             g_0_yyz_0_xxxxy_0,  \
                             g_0_yyz_0_xxxxy_1,  \
                             g_0_yyz_0_xxxyy_0,  \
                             g_0_yyz_0_xxxyy_1,  \
                             g_0_yyz_0_xxyyy_0,  \
                             g_0_yyz_0_xxyyy_1,  \
                             g_0_yyz_0_xyyyy_0,  \
                             g_0_yyz_0_xyyyy_1,  \
                             g_0_yyz_0_yyyyy_0,  \
                             g_0_yyz_0_yyyyy_1,  \
                             g_0_yyzz_0_xxxxx_0, \
                             g_0_yyzz_0_xxxxy_0, \
                             g_0_yyzz_0_xxxxz_0, \
                             g_0_yyzz_0_xxxyy_0, \
                             g_0_yyzz_0_xxxyz_0, \
                             g_0_yyzz_0_xxxzz_0, \
                             g_0_yyzz_0_xxyyy_0, \
                             g_0_yyzz_0_xxyyz_0, \
                             g_0_yyzz_0_xxyzz_0, \
                             g_0_yyzz_0_xxzzz_0, \
                             g_0_yyzz_0_xyyyy_0, \
                             g_0_yyzz_0_xyyyz_0, \
                             g_0_yyzz_0_xyyzz_0, \
                             g_0_yyzz_0_xyzzz_0, \
                             g_0_yyzz_0_xzzzz_0, \
                             g_0_yyzz_0_yyyyy_0, \
                             g_0_yyzz_0_yyyyz_0, \
                             g_0_yyzz_0_yyyzz_0, \
                             g_0_yyzz_0_yyzzz_0, \
                             g_0_yyzz_0_yzzzz_0, \
                             g_0_yyzz_0_zzzzz_0, \
                             g_0_yzz_0_xxxxx_0,  \
                             g_0_yzz_0_xxxxx_1,  \
                             g_0_yzz_0_xxxxz_0,  \
                             g_0_yzz_0_xxxxz_1,  \
                             g_0_yzz_0_xxxyz_0,  \
                             g_0_yzz_0_xxxyz_1,  \
                             g_0_yzz_0_xxxz_1,   \
                             g_0_yzz_0_xxxzz_0,  \
                             g_0_yzz_0_xxxzz_1,  \
                             g_0_yzz_0_xxyyz_0,  \
                             g_0_yzz_0_xxyyz_1,  \
                             g_0_yzz_0_xxyz_1,   \
                             g_0_yzz_0_xxyzz_0,  \
                             g_0_yzz_0_xxyzz_1,  \
                             g_0_yzz_0_xxzz_1,   \
                             g_0_yzz_0_xxzzz_0,  \
                             g_0_yzz_0_xxzzz_1,  \
                             g_0_yzz_0_xyyyz_0,  \
                             g_0_yzz_0_xyyyz_1,  \
                             g_0_yzz_0_xyyz_1,   \
                             g_0_yzz_0_xyyzz_0,  \
                             g_0_yzz_0_xyyzz_1,  \
                             g_0_yzz_0_xyzz_1,   \
                             g_0_yzz_0_xyzzz_0,  \
                             g_0_yzz_0_xyzzz_1,  \
                             g_0_yzz_0_xzzz_1,   \
                             g_0_yzz_0_xzzzz_0,  \
                             g_0_yzz_0_xzzzz_1,  \
                             g_0_yzz_0_yyyyz_0,  \
                             g_0_yzz_0_yyyyz_1,  \
                             g_0_yzz_0_yyyz_1,   \
                             g_0_yzz_0_yyyzz_0,  \
                             g_0_yzz_0_yyyzz_1,  \
                             g_0_yzz_0_yyzz_1,   \
                             g_0_yzz_0_yyzzz_0,  \
                             g_0_yzz_0_yyzzz_1,  \
                             g_0_yzz_0_yzzz_1,   \
                             g_0_yzz_0_yzzzz_0,  \
                             g_0_yzz_0_yzzzz_1,  \
                             g_0_yzz_0_zzzz_1,   \
                             g_0_yzz_0_zzzzz_0,  \
                             g_0_yzz_0_zzzzz_1,  \
                             g_0_zz_0_xxxxx_0,   \
                             g_0_zz_0_xxxxx_1,   \
                             g_0_zz_0_xxxxz_0,   \
                             g_0_zz_0_xxxxz_1,   \
                             g_0_zz_0_xxxyz_0,   \
                             g_0_zz_0_xxxyz_1,   \
                             g_0_zz_0_xxxzz_0,   \
                             g_0_zz_0_xxxzz_1,   \
                             g_0_zz_0_xxyyz_0,   \
                             g_0_zz_0_xxyyz_1,   \
                             g_0_zz_0_xxyzz_0,   \
                             g_0_zz_0_xxyzz_1,   \
                             g_0_zz_0_xxzzz_0,   \
                             g_0_zz_0_xxzzz_1,   \
                             g_0_zz_0_xyyyz_0,   \
                             g_0_zz_0_xyyyz_1,   \
                             g_0_zz_0_xyyzz_0,   \
                             g_0_zz_0_xyyzz_1,   \
                             g_0_zz_0_xyzzz_0,   \
                             g_0_zz_0_xyzzz_1,   \
                             g_0_zz_0_xzzzz_0,   \
                             g_0_zz_0_xzzzz_1,   \
                             g_0_zz_0_yyyyz_0,   \
                             g_0_zz_0_yyyyz_1,   \
                             g_0_zz_0_yyyzz_0,   \
                             g_0_zz_0_yyyzz_1,   \
                             g_0_zz_0_yyzzz_0,   \
                             g_0_zz_0_yyzzz_1,   \
                             g_0_zz_0_yzzzz_0,   \
                             g_0_zz_0_yzzzz_1,   \
                             g_0_zz_0_zzzzz_0,   \
                             g_0_zz_0_zzzzz_1,   \
                             wp_y,               \
                             wp_z,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyzz_0_xxxxx_0[i] =
            g_0_zz_0_xxxxx_0[i] * fi_ab_0 - g_0_zz_0_xxxxx_1[i] * fti_ab_0 + g_0_yzz_0_xxxxx_0[i] * pb_y + g_0_yzz_0_xxxxx_1[i] * wp_y[i];

        g_0_yyzz_0_xxxxy_0[i] =
            g_0_yy_0_xxxxy_0[i] * fi_ab_0 - g_0_yy_0_xxxxy_1[i] * fti_ab_0 + g_0_yyz_0_xxxxy_0[i] * pb_z + g_0_yyz_0_xxxxy_1[i] * wp_z[i];

        g_0_yyzz_0_xxxxz_0[i] =
            g_0_zz_0_xxxxz_0[i] * fi_ab_0 - g_0_zz_0_xxxxz_1[i] * fti_ab_0 + g_0_yzz_0_xxxxz_0[i] * pb_y + g_0_yzz_0_xxxxz_1[i] * wp_y[i];

        g_0_yyzz_0_xxxyy_0[i] =
            g_0_yy_0_xxxyy_0[i] * fi_ab_0 - g_0_yy_0_xxxyy_1[i] * fti_ab_0 + g_0_yyz_0_xxxyy_0[i] * pb_z + g_0_yyz_0_xxxyy_1[i] * wp_z[i];

        g_0_yyzz_0_xxxyz_0[i] = g_0_zz_0_xxxyz_0[i] * fi_ab_0 - g_0_zz_0_xxxyz_1[i] * fti_ab_0 + g_0_yzz_0_xxxz_1[i] * fi_abcd_0 +
                                g_0_yzz_0_xxxyz_0[i] * pb_y + g_0_yzz_0_xxxyz_1[i] * wp_y[i];

        g_0_yyzz_0_xxxzz_0[i] =
            g_0_zz_0_xxxzz_0[i] * fi_ab_0 - g_0_zz_0_xxxzz_1[i] * fti_ab_0 + g_0_yzz_0_xxxzz_0[i] * pb_y + g_0_yzz_0_xxxzz_1[i] * wp_y[i];

        g_0_yyzz_0_xxyyy_0[i] =
            g_0_yy_0_xxyyy_0[i] * fi_ab_0 - g_0_yy_0_xxyyy_1[i] * fti_ab_0 + g_0_yyz_0_xxyyy_0[i] * pb_z + g_0_yyz_0_xxyyy_1[i] * wp_z[i];

        g_0_yyzz_0_xxyyz_0[i] = g_0_zz_0_xxyyz_0[i] * fi_ab_0 - g_0_zz_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yzz_0_xxyz_1[i] * fi_abcd_0 +
                                g_0_yzz_0_xxyyz_0[i] * pb_y + g_0_yzz_0_xxyyz_1[i] * wp_y[i];

        g_0_yyzz_0_xxyzz_0[i] = g_0_zz_0_xxyzz_0[i] * fi_ab_0 - g_0_zz_0_xxyzz_1[i] * fti_ab_0 + g_0_yzz_0_xxzz_1[i] * fi_abcd_0 +
                                g_0_yzz_0_xxyzz_0[i] * pb_y + g_0_yzz_0_xxyzz_1[i] * wp_y[i];

        g_0_yyzz_0_xxzzz_0[i] =
            g_0_zz_0_xxzzz_0[i] * fi_ab_0 - g_0_zz_0_xxzzz_1[i] * fti_ab_0 + g_0_yzz_0_xxzzz_0[i] * pb_y + g_0_yzz_0_xxzzz_1[i] * wp_y[i];

        g_0_yyzz_0_xyyyy_0[i] =
            g_0_yy_0_xyyyy_0[i] * fi_ab_0 - g_0_yy_0_xyyyy_1[i] * fti_ab_0 + g_0_yyz_0_xyyyy_0[i] * pb_z + g_0_yyz_0_xyyyy_1[i] * wp_z[i];

        g_0_yyzz_0_xyyyz_0[i] = g_0_zz_0_xyyyz_0[i] * fi_ab_0 - g_0_zz_0_xyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yzz_0_xyyz_1[i] * fi_abcd_0 +
                                g_0_yzz_0_xyyyz_0[i] * pb_y + g_0_yzz_0_xyyyz_1[i] * wp_y[i];

        g_0_yyzz_0_xyyzz_0[i] = g_0_zz_0_xyyzz_0[i] * fi_ab_0 - g_0_zz_0_xyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yzz_0_xyzz_1[i] * fi_abcd_0 +
                                g_0_yzz_0_xyyzz_0[i] * pb_y + g_0_yzz_0_xyyzz_1[i] * wp_y[i];

        g_0_yyzz_0_xyzzz_0[i] = g_0_zz_0_xyzzz_0[i] * fi_ab_0 - g_0_zz_0_xyzzz_1[i] * fti_ab_0 + g_0_yzz_0_xzzz_1[i] * fi_abcd_0 +
                                g_0_yzz_0_xyzzz_0[i] * pb_y + g_0_yzz_0_xyzzz_1[i] * wp_y[i];

        g_0_yyzz_0_xzzzz_0[i] =
            g_0_zz_0_xzzzz_0[i] * fi_ab_0 - g_0_zz_0_xzzzz_1[i] * fti_ab_0 + g_0_yzz_0_xzzzz_0[i] * pb_y + g_0_yzz_0_xzzzz_1[i] * wp_y[i];

        g_0_yyzz_0_yyyyy_0[i] =
            g_0_yy_0_yyyyy_0[i] * fi_ab_0 - g_0_yy_0_yyyyy_1[i] * fti_ab_0 + g_0_yyz_0_yyyyy_0[i] * pb_z + g_0_yyz_0_yyyyy_1[i] * wp_z[i];

        g_0_yyzz_0_yyyyz_0[i] = g_0_zz_0_yyyyz_0[i] * fi_ab_0 - g_0_zz_0_yyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yzz_0_yyyz_1[i] * fi_abcd_0 +
                                g_0_yzz_0_yyyyz_0[i] * pb_y + g_0_yzz_0_yyyyz_1[i] * wp_y[i];

        g_0_yyzz_0_yyyzz_0[i] = g_0_zz_0_yyyzz_0[i] * fi_ab_0 - g_0_zz_0_yyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yzz_0_yyzz_1[i] * fi_abcd_0 +
                                g_0_yzz_0_yyyzz_0[i] * pb_y + g_0_yzz_0_yyyzz_1[i] * wp_y[i];

        g_0_yyzz_0_yyzzz_0[i] = g_0_zz_0_yyzzz_0[i] * fi_ab_0 - g_0_zz_0_yyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yzz_0_yzzz_1[i] * fi_abcd_0 +
                                g_0_yzz_0_yyzzz_0[i] * pb_y + g_0_yzz_0_yyzzz_1[i] * wp_y[i];

        g_0_yyzz_0_yzzzz_0[i] = g_0_zz_0_yzzzz_0[i] * fi_ab_0 - g_0_zz_0_yzzzz_1[i] * fti_ab_0 + g_0_yzz_0_zzzz_1[i] * fi_abcd_0 +
                                g_0_yzz_0_yzzzz_0[i] * pb_y + g_0_yzz_0_yzzzz_1[i] * wp_y[i];

        g_0_yyzz_0_zzzzz_0[i] =
            g_0_zz_0_zzzzz_0[i] * fi_ab_0 - g_0_zz_0_zzzzz_1[i] * fti_ab_0 + g_0_yzz_0_zzzzz_0[i] * pb_y + g_0_yzz_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 273-294 components of targeted buffer : SGSH

    auto g_0_yzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sgsh + 273);

    auto g_0_yzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sgsh + 274);

    auto g_0_yzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sgsh + 275);

    auto g_0_yzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sgsh + 276);

    auto g_0_yzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sgsh + 277);

    auto g_0_yzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sgsh + 278);

    auto g_0_yzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sgsh + 279);

    auto g_0_yzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sgsh + 280);

    auto g_0_yzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sgsh + 281);

    auto g_0_yzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sgsh + 282);

    auto g_0_yzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 283);

    auto g_0_yzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 284);

    auto g_0_yzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 285);

    auto g_0_yzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 286);

    auto g_0_yzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 287);

    auto g_0_yzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 288);

    auto g_0_yzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 289);

    auto g_0_yzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 290);

    auto g_0_yzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 291);

    auto g_0_yzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 292);

    auto g_0_yzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 293);

#pragma omp simd aligned(g_0_yzzz_0_xxxxx_0,     \
                             g_0_yzzz_0_xxxxy_0, \
                             g_0_yzzz_0_xxxxz_0, \
                             g_0_yzzz_0_xxxyy_0, \
                             g_0_yzzz_0_xxxyz_0, \
                             g_0_yzzz_0_xxxzz_0, \
                             g_0_yzzz_0_xxyyy_0, \
                             g_0_yzzz_0_xxyyz_0, \
                             g_0_yzzz_0_xxyzz_0, \
                             g_0_yzzz_0_xxzzz_0, \
                             g_0_yzzz_0_xyyyy_0, \
                             g_0_yzzz_0_xyyyz_0, \
                             g_0_yzzz_0_xyyzz_0, \
                             g_0_yzzz_0_xyzzz_0, \
                             g_0_yzzz_0_xzzzz_0, \
                             g_0_yzzz_0_yyyyy_0, \
                             g_0_yzzz_0_yyyyz_0, \
                             g_0_yzzz_0_yyyzz_0, \
                             g_0_yzzz_0_yyzzz_0, \
                             g_0_yzzz_0_yzzzz_0, \
                             g_0_yzzz_0_zzzzz_0, \
                             g_0_zzz_0_xxxx_1,   \
                             g_0_zzz_0_xxxxx_0,  \
                             g_0_zzz_0_xxxxx_1,  \
                             g_0_zzz_0_xxxxy_0,  \
                             g_0_zzz_0_xxxxy_1,  \
                             g_0_zzz_0_xxxxz_0,  \
                             g_0_zzz_0_xxxxz_1,  \
                             g_0_zzz_0_xxxy_1,   \
                             g_0_zzz_0_xxxyy_0,  \
                             g_0_zzz_0_xxxyy_1,  \
                             g_0_zzz_0_xxxyz_0,  \
                             g_0_zzz_0_xxxyz_1,  \
                             g_0_zzz_0_xxxz_1,   \
                             g_0_zzz_0_xxxzz_0,  \
                             g_0_zzz_0_xxxzz_1,  \
                             g_0_zzz_0_xxyy_1,   \
                             g_0_zzz_0_xxyyy_0,  \
                             g_0_zzz_0_xxyyy_1,  \
                             g_0_zzz_0_xxyyz_0,  \
                             g_0_zzz_0_xxyyz_1,  \
                             g_0_zzz_0_xxyz_1,   \
                             g_0_zzz_0_xxyzz_0,  \
                             g_0_zzz_0_xxyzz_1,  \
                             g_0_zzz_0_xxzz_1,   \
                             g_0_zzz_0_xxzzz_0,  \
                             g_0_zzz_0_xxzzz_1,  \
                             g_0_zzz_0_xyyy_1,   \
                             g_0_zzz_0_xyyyy_0,  \
                             g_0_zzz_0_xyyyy_1,  \
                             g_0_zzz_0_xyyyz_0,  \
                             g_0_zzz_0_xyyyz_1,  \
                             g_0_zzz_0_xyyz_1,   \
                             g_0_zzz_0_xyyzz_0,  \
                             g_0_zzz_0_xyyzz_1,  \
                             g_0_zzz_0_xyzz_1,   \
                             g_0_zzz_0_xyzzz_0,  \
                             g_0_zzz_0_xyzzz_1,  \
                             g_0_zzz_0_xzzz_1,   \
                             g_0_zzz_0_xzzzz_0,  \
                             g_0_zzz_0_xzzzz_1,  \
                             g_0_zzz_0_yyyy_1,   \
                             g_0_zzz_0_yyyyy_0,  \
                             g_0_zzz_0_yyyyy_1,  \
                             g_0_zzz_0_yyyyz_0,  \
                             g_0_zzz_0_yyyyz_1,  \
                             g_0_zzz_0_yyyz_1,   \
                             g_0_zzz_0_yyyzz_0,  \
                             g_0_zzz_0_yyyzz_1,  \
                             g_0_zzz_0_yyzz_1,   \
                             g_0_zzz_0_yyzzz_0,  \
                             g_0_zzz_0_yyzzz_1,  \
                             g_0_zzz_0_yzzz_1,   \
                             g_0_zzz_0_yzzzz_0,  \
                             g_0_zzz_0_yzzzz_1,  \
                             g_0_zzz_0_zzzz_1,   \
                             g_0_zzz_0_zzzzz_0,  \
                             g_0_zzz_0_zzzzz_1,  \
                             wp_y,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzzz_0_xxxxx_0[i] = g_0_zzz_0_xxxxx_0[i] * pb_y + g_0_zzz_0_xxxxx_1[i] * wp_y[i];

        g_0_yzzz_0_xxxxy_0[i] = g_0_zzz_0_xxxx_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxy_0[i] * pb_y + g_0_zzz_0_xxxxy_1[i] * wp_y[i];

        g_0_yzzz_0_xxxxz_0[i] = g_0_zzz_0_xxxxz_0[i] * pb_y + g_0_zzz_0_xxxxz_1[i] * wp_y[i];

        g_0_yzzz_0_xxxyy_0[i] = 2.0 * g_0_zzz_0_xxxy_1[i] * fi_abcd_0 + g_0_zzz_0_xxxyy_0[i] * pb_y + g_0_zzz_0_xxxyy_1[i] * wp_y[i];

        g_0_yzzz_0_xxxyz_0[i] = g_0_zzz_0_xxxz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxyz_0[i] * pb_y + g_0_zzz_0_xxxyz_1[i] * wp_y[i];

        g_0_yzzz_0_xxxzz_0[i] = g_0_zzz_0_xxxzz_0[i] * pb_y + g_0_zzz_0_xxxzz_1[i] * wp_y[i];

        g_0_yzzz_0_xxyyy_0[i] = 3.0 * g_0_zzz_0_xxyy_1[i] * fi_abcd_0 + g_0_zzz_0_xxyyy_0[i] * pb_y + g_0_zzz_0_xxyyy_1[i] * wp_y[i];

        g_0_yzzz_0_xxyyz_0[i] = 2.0 * g_0_zzz_0_xxyz_1[i] * fi_abcd_0 + g_0_zzz_0_xxyyz_0[i] * pb_y + g_0_zzz_0_xxyyz_1[i] * wp_y[i];

        g_0_yzzz_0_xxyzz_0[i] = g_0_zzz_0_xxzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxyzz_0[i] * pb_y + g_0_zzz_0_xxyzz_1[i] * wp_y[i];

        g_0_yzzz_0_xxzzz_0[i] = g_0_zzz_0_xxzzz_0[i] * pb_y + g_0_zzz_0_xxzzz_1[i] * wp_y[i];

        g_0_yzzz_0_xyyyy_0[i] = 4.0 * g_0_zzz_0_xyyy_1[i] * fi_abcd_0 + g_0_zzz_0_xyyyy_0[i] * pb_y + g_0_zzz_0_xyyyy_1[i] * wp_y[i];

        g_0_yzzz_0_xyyyz_0[i] = 3.0 * g_0_zzz_0_xyyz_1[i] * fi_abcd_0 + g_0_zzz_0_xyyyz_0[i] * pb_y + g_0_zzz_0_xyyyz_1[i] * wp_y[i];

        g_0_yzzz_0_xyyzz_0[i] = 2.0 * g_0_zzz_0_xyzz_1[i] * fi_abcd_0 + g_0_zzz_0_xyyzz_0[i] * pb_y + g_0_zzz_0_xyyzz_1[i] * wp_y[i];

        g_0_yzzz_0_xyzzz_0[i] = g_0_zzz_0_xzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xyzzz_0[i] * pb_y + g_0_zzz_0_xyzzz_1[i] * wp_y[i];

        g_0_yzzz_0_xzzzz_0[i] = g_0_zzz_0_xzzzz_0[i] * pb_y + g_0_zzz_0_xzzzz_1[i] * wp_y[i];

        g_0_yzzz_0_yyyyy_0[i] = 5.0 * g_0_zzz_0_yyyy_1[i] * fi_abcd_0 + g_0_zzz_0_yyyyy_0[i] * pb_y + g_0_zzz_0_yyyyy_1[i] * wp_y[i];

        g_0_yzzz_0_yyyyz_0[i] = 4.0 * g_0_zzz_0_yyyz_1[i] * fi_abcd_0 + g_0_zzz_0_yyyyz_0[i] * pb_y + g_0_zzz_0_yyyyz_1[i] * wp_y[i];

        g_0_yzzz_0_yyyzz_0[i] = 3.0 * g_0_zzz_0_yyzz_1[i] * fi_abcd_0 + g_0_zzz_0_yyyzz_0[i] * pb_y + g_0_zzz_0_yyyzz_1[i] * wp_y[i];

        g_0_yzzz_0_yyzzz_0[i] = 2.0 * g_0_zzz_0_yzzz_1[i] * fi_abcd_0 + g_0_zzz_0_yyzzz_0[i] * pb_y + g_0_zzz_0_yyzzz_1[i] * wp_y[i];

        g_0_yzzz_0_yzzzz_0[i] = g_0_zzz_0_zzzz_1[i] * fi_abcd_0 + g_0_zzz_0_yzzzz_0[i] * pb_y + g_0_zzz_0_yzzzz_1[i] * wp_y[i];

        g_0_yzzz_0_zzzzz_0[i] = g_0_zzz_0_zzzzz_0[i] * pb_y + g_0_zzz_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 294-315 components of targeted buffer : SGSH

    auto g_0_zzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sgsh + 294);

    auto g_0_zzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sgsh + 295);

    auto g_0_zzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sgsh + 296);

    auto g_0_zzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sgsh + 297);

    auto g_0_zzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sgsh + 298);

    auto g_0_zzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sgsh + 299);

    auto g_0_zzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sgsh + 300);

    auto g_0_zzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sgsh + 301);

    auto g_0_zzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sgsh + 302);

    auto g_0_zzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sgsh + 303);

    auto g_0_zzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 304);

    auto g_0_zzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 305);

    auto g_0_zzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 306);

    auto g_0_zzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 307);

    auto g_0_zzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 308);

    auto g_0_zzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sgsh + 309);

    auto g_0_zzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sgsh + 310);

    auto g_0_zzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sgsh + 311);

    auto g_0_zzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sgsh + 312);

    auto g_0_zzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 313);

    auto g_0_zzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sgsh + 314);

#pragma omp simd aligned(g_0_zz_0_xxxxx_0,       \
                             g_0_zz_0_xxxxx_1,   \
                             g_0_zz_0_xxxxy_0,   \
                             g_0_zz_0_xxxxy_1,   \
                             g_0_zz_0_xxxxz_0,   \
                             g_0_zz_0_xxxxz_1,   \
                             g_0_zz_0_xxxyy_0,   \
                             g_0_zz_0_xxxyy_1,   \
                             g_0_zz_0_xxxyz_0,   \
                             g_0_zz_0_xxxyz_1,   \
                             g_0_zz_0_xxxzz_0,   \
                             g_0_zz_0_xxxzz_1,   \
                             g_0_zz_0_xxyyy_0,   \
                             g_0_zz_0_xxyyy_1,   \
                             g_0_zz_0_xxyyz_0,   \
                             g_0_zz_0_xxyyz_1,   \
                             g_0_zz_0_xxyzz_0,   \
                             g_0_zz_0_xxyzz_1,   \
                             g_0_zz_0_xxzzz_0,   \
                             g_0_zz_0_xxzzz_1,   \
                             g_0_zz_0_xyyyy_0,   \
                             g_0_zz_0_xyyyy_1,   \
                             g_0_zz_0_xyyyz_0,   \
                             g_0_zz_0_xyyyz_1,   \
                             g_0_zz_0_xyyzz_0,   \
                             g_0_zz_0_xyyzz_1,   \
                             g_0_zz_0_xyzzz_0,   \
                             g_0_zz_0_xyzzz_1,   \
                             g_0_zz_0_xzzzz_0,   \
                             g_0_zz_0_xzzzz_1,   \
                             g_0_zz_0_yyyyy_0,   \
                             g_0_zz_0_yyyyy_1,   \
                             g_0_zz_0_yyyyz_0,   \
                             g_0_zz_0_yyyyz_1,   \
                             g_0_zz_0_yyyzz_0,   \
                             g_0_zz_0_yyyzz_1,   \
                             g_0_zz_0_yyzzz_0,   \
                             g_0_zz_0_yyzzz_1,   \
                             g_0_zz_0_yzzzz_0,   \
                             g_0_zz_0_yzzzz_1,   \
                             g_0_zz_0_zzzzz_0,   \
                             g_0_zz_0_zzzzz_1,   \
                             g_0_zzz_0_xxxx_1,   \
                             g_0_zzz_0_xxxxx_0,  \
                             g_0_zzz_0_xxxxx_1,  \
                             g_0_zzz_0_xxxxy_0,  \
                             g_0_zzz_0_xxxxy_1,  \
                             g_0_zzz_0_xxxxz_0,  \
                             g_0_zzz_0_xxxxz_1,  \
                             g_0_zzz_0_xxxy_1,   \
                             g_0_zzz_0_xxxyy_0,  \
                             g_0_zzz_0_xxxyy_1,  \
                             g_0_zzz_0_xxxyz_0,  \
                             g_0_zzz_0_xxxyz_1,  \
                             g_0_zzz_0_xxxz_1,   \
                             g_0_zzz_0_xxxzz_0,  \
                             g_0_zzz_0_xxxzz_1,  \
                             g_0_zzz_0_xxyy_1,   \
                             g_0_zzz_0_xxyyy_0,  \
                             g_0_zzz_0_xxyyy_1,  \
                             g_0_zzz_0_xxyyz_0,  \
                             g_0_zzz_0_xxyyz_1,  \
                             g_0_zzz_0_xxyz_1,   \
                             g_0_zzz_0_xxyzz_0,  \
                             g_0_zzz_0_xxyzz_1,  \
                             g_0_zzz_0_xxzz_1,   \
                             g_0_zzz_0_xxzzz_0,  \
                             g_0_zzz_0_xxzzz_1,  \
                             g_0_zzz_0_xyyy_1,   \
                             g_0_zzz_0_xyyyy_0,  \
                             g_0_zzz_0_xyyyy_1,  \
                             g_0_zzz_0_xyyyz_0,  \
                             g_0_zzz_0_xyyyz_1,  \
                             g_0_zzz_0_xyyz_1,   \
                             g_0_zzz_0_xyyzz_0,  \
                             g_0_zzz_0_xyyzz_1,  \
                             g_0_zzz_0_xyzz_1,   \
                             g_0_zzz_0_xyzzz_0,  \
                             g_0_zzz_0_xyzzz_1,  \
                             g_0_zzz_0_xzzz_1,   \
                             g_0_zzz_0_xzzzz_0,  \
                             g_0_zzz_0_xzzzz_1,  \
                             g_0_zzz_0_yyyy_1,   \
                             g_0_zzz_0_yyyyy_0,  \
                             g_0_zzz_0_yyyyy_1,  \
                             g_0_zzz_0_yyyyz_0,  \
                             g_0_zzz_0_yyyyz_1,  \
                             g_0_zzz_0_yyyz_1,   \
                             g_0_zzz_0_yyyzz_0,  \
                             g_0_zzz_0_yyyzz_1,  \
                             g_0_zzz_0_yyzz_1,   \
                             g_0_zzz_0_yyzzz_0,  \
                             g_0_zzz_0_yyzzz_1,  \
                             g_0_zzz_0_yzzz_1,   \
                             g_0_zzz_0_yzzzz_0,  \
                             g_0_zzz_0_yzzzz_1,  \
                             g_0_zzz_0_zzzz_1,   \
                             g_0_zzz_0_zzzzz_0,  \
                             g_0_zzz_0_zzzzz_1,  \
                             g_0_zzzz_0_xxxxx_0, \
                             g_0_zzzz_0_xxxxy_0, \
                             g_0_zzzz_0_xxxxz_0, \
                             g_0_zzzz_0_xxxyy_0, \
                             g_0_zzzz_0_xxxyz_0, \
                             g_0_zzzz_0_xxxzz_0, \
                             g_0_zzzz_0_xxyyy_0, \
                             g_0_zzzz_0_xxyyz_0, \
                             g_0_zzzz_0_xxyzz_0, \
                             g_0_zzzz_0_xxzzz_0, \
                             g_0_zzzz_0_xyyyy_0, \
                             g_0_zzzz_0_xyyyz_0, \
                             g_0_zzzz_0_xyyzz_0, \
                             g_0_zzzz_0_xyzzz_0, \
                             g_0_zzzz_0_xzzzz_0, \
                             g_0_zzzz_0_yyyyy_0, \
                             g_0_zzzz_0_yyyyz_0, \
                             g_0_zzzz_0_yyyzz_0, \
                             g_0_zzzz_0_yyzzz_0, \
                             g_0_zzzz_0_yzzzz_0, \
                             g_0_zzzz_0_zzzzz_0, \
                             wp_z,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzzz_0_xxxxx_0[i] =
            3.0 * g_0_zz_0_xxxxx_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxxx_1[i] * fti_ab_0 + g_0_zzz_0_xxxxx_0[i] * pb_z + g_0_zzz_0_xxxxx_1[i] * wp_z[i];

        g_0_zzzz_0_xxxxy_0[i] =
            3.0 * g_0_zz_0_xxxxy_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxxy_1[i] * fti_ab_0 + g_0_zzz_0_xxxxy_0[i] * pb_z + g_0_zzz_0_xxxxy_1[i] * wp_z[i];

        g_0_zzzz_0_xxxxz_0[i] = 3.0 * g_0_zz_0_xxxxz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxxz_1[i] * fti_ab_0 + g_0_zzz_0_xxxx_1[i] * fi_abcd_0 +
                                g_0_zzz_0_xxxxz_0[i] * pb_z + g_0_zzz_0_xxxxz_1[i] * wp_z[i];

        g_0_zzzz_0_xxxyy_0[i] =
            3.0 * g_0_zz_0_xxxyy_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxyy_1[i] * fti_ab_0 + g_0_zzz_0_xxxyy_0[i] * pb_z + g_0_zzz_0_xxxyy_1[i] * wp_z[i];

        g_0_zzzz_0_xxxyz_0[i] = 3.0 * g_0_zz_0_xxxyz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxyz_1[i] * fti_ab_0 + g_0_zzz_0_xxxy_1[i] * fi_abcd_0 +
                                g_0_zzz_0_xxxyz_0[i] * pb_z + g_0_zzz_0_xxxyz_1[i] * wp_z[i];

        g_0_zzzz_0_xxxzz_0[i] = 3.0 * g_0_zz_0_xxxzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxzz_1[i] * fti_ab_0 + 2.0 * g_0_zzz_0_xxxz_1[i] * fi_abcd_0 +
                                g_0_zzz_0_xxxzz_0[i] * pb_z + g_0_zzz_0_xxxzz_1[i] * wp_z[i];

        g_0_zzzz_0_xxyyy_0[i] =
            3.0 * g_0_zz_0_xxyyy_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxyyy_1[i] * fti_ab_0 + g_0_zzz_0_xxyyy_0[i] * pb_z + g_0_zzz_0_xxyyy_1[i] * wp_z[i];

        g_0_zzzz_0_xxyyz_0[i] = 3.0 * g_0_zz_0_xxyyz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxyyz_1[i] * fti_ab_0 + g_0_zzz_0_xxyy_1[i] * fi_abcd_0 +
                                g_0_zzz_0_xxyyz_0[i] * pb_z + g_0_zzz_0_xxyyz_1[i] * wp_z[i];

        g_0_zzzz_0_xxyzz_0[i] = 3.0 * g_0_zz_0_xxyzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzz_0_xxyz_1[i] * fi_abcd_0 +
                                g_0_zzz_0_xxyzz_0[i] * pb_z + g_0_zzz_0_xxyzz_1[i] * wp_z[i];

        g_0_zzzz_0_xxzzz_0[i] = 3.0 * g_0_zz_0_xxzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzz_0_xxzz_1[i] * fi_abcd_0 +
                                g_0_zzz_0_xxzzz_0[i] * pb_z + g_0_zzz_0_xxzzz_1[i] * wp_z[i];

        g_0_zzzz_0_xyyyy_0[i] =
            3.0 * g_0_zz_0_xyyyy_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xyyyy_1[i] * fti_ab_0 + g_0_zzz_0_xyyyy_0[i] * pb_z + g_0_zzz_0_xyyyy_1[i] * wp_z[i];

        g_0_zzzz_0_xyyyz_0[i] = 3.0 * g_0_zz_0_xyyyz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xyyyz_1[i] * fti_ab_0 + g_0_zzz_0_xyyy_1[i] * fi_abcd_0 +
                                g_0_zzz_0_xyyyz_0[i] * pb_z + g_0_zzz_0_xyyyz_1[i] * wp_z[i];

        g_0_zzzz_0_xyyzz_0[i] = 3.0 * g_0_zz_0_xyyzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzz_0_xyyz_1[i] * fi_abcd_0 +
                                g_0_zzz_0_xyyzz_0[i] * pb_z + g_0_zzz_0_xyyzz_1[i] * wp_z[i];

        g_0_zzzz_0_xyzzz_0[i] = 3.0 * g_0_zz_0_xyzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzz_0_xyzz_1[i] * fi_abcd_0 +
                                g_0_zzz_0_xyzzz_0[i] * pb_z + g_0_zzz_0_xyzzz_1[i] * wp_z[i];

        g_0_zzzz_0_xzzzz_0[i] = 3.0 * g_0_zz_0_xzzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zzz_0_xzzz_1[i] * fi_abcd_0 +
                                g_0_zzz_0_xzzzz_0[i] * pb_z + g_0_zzz_0_xzzzz_1[i] * wp_z[i];

        g_0_zzzz_0_yyyyy_0[i] =
            3.0 * g_0_zz_0_yyyyy_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_yyyyy_1[i] * fti_ab_0 + g_0_zzz_0_yyyyy_0[i] * pb_z + g_0_zzz_0_yyyyy_1[i] * wp_z[i];

        g_0_zzzz_0_yyyyz_0[i] = 3.0 * g_0_zz_0_yyyyz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_yyyyz_1[i] * fti_ab_0 + g_0_zzz_0_yyyy_1[i] * fi_abcd_0 +
                                g_0_zzz_0_yyyyz_0[i] * pb_z + g_0_zzz_0_yyyyz_1[i] * wp_z[i];

        g_0_zzzz_0_yyyzz_0[i] = 3.0 * g_0_zz_0_yyyzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_yyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzz_0_yyyz_1[i] * fi_abcd_0 +
                                g_0_zzz_0_yyyzz_0[i] * pb_z + g_0_zzz_0_yyyzz_1[i] * wp_z[i];

        g_0_zzzz_0_yyzzz_0[i] = 3.0 * g_0_zz_0_yyzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_yyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzz_0_yyzz_1[i] * fi_abcd_0 +
                                g_0_zzz_0_yyzzz_0[i] * pb_z + g_0_zzz_0_yyzzz_1[i] * wp_z[i];

        g_0_zzzz_0_yzzzz_0[i] = 3.0 * g_0_zz_0_yzzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_yzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zzz_0_yzzz_1[i] * fi_abcd_0 +
                                g_0_zzz_0_yzzzz_0[i] * pb_z + g_0_zzz_0_yzzzz_1[i] * wp_z[i];

        g_0_zzzz_0_zzzzz_0[i] = 3.0 * g_0_zz_0_zzzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_zzzzz_1[i] * fti_ab_0 + 5.0 * g_0_zzz_0_zzzz_1[i] * fi_abcd_0 +
                                g_0_zzz_0_zzzzz_0[i] * pb_z + g_0_zzz_0_zzzzz_1[i] * wp_z[i];
    }
}

}  // namespace erirec

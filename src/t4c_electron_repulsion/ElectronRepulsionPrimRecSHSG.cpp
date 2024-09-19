#include "ElectronRepulsionPrimRecSHSG.hpp"

namespace erirec {  // erirec namespace

auto
comp_prim_electron_repulsion_shsg(CSimdArray<double>&   pbuffer,
                                  const size_t          idx_eri_0_shsg,
                                  size_t                idx_eri_0_sfsg,
                                  size_t                idx_eri_1_sfsg,
                                  size_t                idx_eri_1_sgsf,
                                  size_t                idx_eri_0_sgsg,
                                  size_t                idx_eri_1_sgsg,
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

    /// Set up components of auxilary buffer : SFSG

    auto g_0_xxx_0_xxxx_0 = pbuffer.data(idx_eri_0_sfsg);

    auto g_0_xxx_0_xxxy_0 = pbuffer.data(idx_eri_0_sfsg + 1);

    auto g_0_xxx_0_xxxz_0 = pbuffer.data(idx_eri_0_sfsg + 2);

    auto g_0_xxx_0_xxyy_0 = pbuffer.data(idx_eri_0_sfsg + 3);

    auto g_0_xxx_0_xxyz_0 = pbuffer.data(idx_eri_0_sfsg + 4);

    auto g_0_xxx_0_xxzz_0 = pbuffer.data(idx_eri_0_sfsg + 5);

    auto g_0_xxx_0_xyyy_0 = pbuffer.data(idx_eri_0_sfsg + 6);

    auto g_0_xxx_0_xyyz_0 = pbuffer.data(idx_eri_0_sfsg + 7);

    auto g_0_xxx_0_xyzz_0 = pbuffer.data(idx_eri_0_sfsg + 8);

    auto g_0_xxx_0_xzzz_0 = pbuffer.data(idx_eri_0_sfsg + 9);

    auto g_0_xxx_0_yyyy_0 = pbuffer.data(idx_eri_0_sfsg + 10);

    auto g_0_xxx_0_yyyz_0 = pbuffer.data(idx_eri_0_sfsg + 11);

    auto g_0_xxx_0_yyzz_0 = pbuffer.data(idx_eri_0_sfsg + 12);

    auto g_0_xxx_0_yzzz_0 = pbuffer.data(idx_eri_0_sfsg + 13);

    auto g_0_xxx_0_zzzz_0 = pbuffer.data(idx_eri_0_sfsg + 14);

    auto g_0_xxy_0_xxxx_0 = pbuffer.data(idx_eri_0_sfsg + 15);

    auto g_0_xxy_0_xxxz_0 = pbuffer.data(idx_eri_0_sfsg + 17);

    auto g_0_xxy_0_xxzz_0 = pbuffer.data(idx_eri_0_sfsg + 20);

    auto g_0_xxy_0_xzzz_0 = pbuffer.data(idx_eri_0_sfsg + 24);

    auto g_0_xxz_0_xxxx_0 = pbuffer.data(idx_eri_0_sfsg + 30);

    auto g_0_xxz_0_xxxy_0 = pbuffer.data(idx_eri_0_sfsg + 31);

    auto g_0_xxz_0_xxyy_0 = pbuffer.data(idx_eri_0_sfsg + 33);

    auto g_0_xxz_0_xyyy_0 = pbuffer.data(idx_eri_0_sfsg + 36);

    auto g_0_xyy_0_xxxy_0 = pbuffer.data(idx_eri_0_sfsg + 46);

    auto g_0_xyy_0_xxyy_0 = pbuffer.data(idx_eri_0_sfsg + 48);

    auto g_0_xyy_0_xxyz_0 = pbuffer.data(idx_eri_0_sfsg + 49);

    auto g_0_xyy_0_xyyy_0 = pbuffer.data(idx_eri_0_sfsg + 51);

    auto g_0_xyy_0_xyyz_0 = pbuffer.data(idx_eri_0_sfsg + 52);

    auto g_0_xyy_0_xyzz_0 = pbuffer.data(idx_eri_0_sfsg + 53);

    auto g_0_xyy_0_yyyy_0 = pbuffer.data(idx_eri_0_sfsg + 55);

    auto g_0_xyy_0_yyyz_0 = pbuffer.data(idx_eri_0_sfsg + 56);

    auto g_0_xyy_0_yyzz_0 = pbuffer.data(idx_eri_0_sfsg + 57);

    auto g_0_xyy_0_yzzz_0 = pbuffer.data(idx_eri_0_sfsg + 58);

    auto g_0_xyy_0_zzzz_0 = pbuffer.data(idx_eri_0_sfsg + 59);

    auto g_0_xzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sfsg + 77);

    auto g_0_xzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sfsg + 79);

    auto g_0_xzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sfsg + 80);

    auto g_0_xzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sfsg + 82);

    auto g_0_xzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sfsg + 83);

    auto g_0_xzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sfsg + 84);

    auto g_0_xzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sfsg + 85);

    auto g_0_xzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sfsg + 86);

    auto g_0_xzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sfsg + 87);

    auto g_0_xzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sfsg + 88);

    auto g_0_xzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sfsg + 89);

    auto g_0_yyy_0_xxxx_0 = pbuffer.data(idx_eri_0_sfsg + 90);

    auto g_0_yyy_0_xxxy_0 = pbuffer.data(idx_eri_0_sfsg + 91);

    auto g_0_yyy_0_xxxz_0 = pbuffer.data(idx_eri_0_sfsg + 92);

    auto g_0_yyy_0_xxyy_0 = pbuffer.data(idx_eri_0_sfsg + 93);

    auto g_0_yyy_0_xxyz_0 = pbuffer.data(idx_eri_0_sfsg + 94);

    auto g_0_yyy_0_xxzz_0 = pbuffer.data(idx_eri_0_sfsg + 95);

    auto g_0_yyy_0_xyyy_0 = pbuffer.data(idx_eri_0_sfsg + 96);

    auto g_0_yyy_0_xyyz_0 = pbuffer.data(idx_eri_0_sfsg + 97);

    auto g_0_yyy_0_xyzz_0 = pbuffer.data(idx_eri_0_sfsg + 98);

    auto g_0_yyy_0_xzzz_0 = pbuffer.data(idx_eri_0_sfsg + 99);

    auto g_0_yyy_0_yyyy_0 = pbuffer.data(idx_eri_0_sfsg + 100);

    auto g_0_yyy_0_yyyz_0 = pbuffer.data(idx_eri_0_sfsg + 101);

    auto g_0_yyy_0_yyzz_0 = pbuffer.data(idx_eri_0_sfsg + 102);

    auto g_0_yyy_0_yzzz_0 = pbuffer.data(idx_eri_0_sfsg + 103);

    auto g_0_yyy_0_zzzz_0 = pbuffer.data(idx_eri_0_sfsg + 104);

    auto g_0_yyz_0_xxxy_0 = pbuffer.data(idx_eri_0_sfsg + 106);

    auto g_0_yyz_0_xxyy_0 = pbuffer.data(idx_eri_0_sfsg + 108);

    auto g_0_yyz_0_xyyy_0 = pbuffer.data(idx_eri_0_sfsg + 111);

    auto g_0_yyz_0_yyyy_0 = pbuffer.data(idx_eri_0_sfsg + 115);

    auto g_0_yzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sfsg + 120);

    auto g_0_yzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sfsg + 122);

    auto g_0_yzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sfsg + 124);

    auto g_0_yzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sfsg + 125);

    auto g_0_yzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sfsg + 127);

    auto g_0_yzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sfsg + 128);

    auto g_0_yzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sfsg + 129);

    auto g_0_yzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sfsg + 131);

    auto g_0_yzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sfsg + 132);

    auto g_0_yzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sfsg + 133);

    auto g_0_yzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sfsg + 134);

    auto g_0_zzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sfsg + 135);

    auto g_0_zzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sfsg + 136);

    auto g_0_zzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sfsg + 137);

    auto g_0_zzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sfsg + 138);

    auto g_0_zzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sfsg + 139);

    auto g_0_zzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sfsg + 140);

    auto g_0_zzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sfsg + 141);

    auto g_0_zzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sfsg + 142);

    auto g_0_zzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sfsg + 143);

    auto g_0_zzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sfsg + 144);

    auto g_0_zzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sfsg + 145);

    auto g_0_zzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sfsg + 146);

    auto g_0_zzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sfsg + 147);

    auto g_0_zzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sfsg + 148);

    auto g_0_zzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sfsg + 149);

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

    auto g_0_xxy_0_xxxx_1 = pbuffer.data(idx_eri_1_sfsg + 15);

    auto g_0_xxy_0_xxxz_1 = pbuffer.data(idx_eri_1_sfsg + 17);

    auto g_0_xxy_0_xxzz_1 = pbuffer.data(idx_eri_1_sfsg + 20);

    auto g_0_xxy_0_xzzz_1 = pbuffer.data(idx_eri_1_sfsg + 24);

    auto g_0_xxz_0_xxxx_1 = pbuffer.data(idx_eri_1_sfsg + 30);

    auto g_0_xxz_0_xxxy_1 = pbuffer.data(idx_eri_1_sfsg + 31);

    auto g_0_xxz_0_xxyy_1 = pbuffer.data(idx_eri_1_sfsg + 33);

    auto g_0_xxz_0_xyyy_1 = pbuffer.data(idx_eri_1_sfsg + 36);

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

    auto g_0_xyy_0_zzzz_1 = pbuffer.data(idx_eri_1_sfsg + 59);

    auto g_0_xzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sfsg + 77);

    auto g_0_xzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sfsg + 79);

    auto g_0_xzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sfsg + 80);

    auto g_0_xzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sfsg + 82);

    auto g_0_xzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sfsg + 83);

    auto g_0_xzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sfsg + 84);

    auto g_0_xzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sfsg + 85);

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

    auto g_0_yyz_0_xxxy_1 = pbuffer.data(idx_eri_1_sfsg + 106);

    auto g_0_yyz_0_xxyy_1 = pbuffer.data(idx_eri_1_sfsg + 108);

    auto g_0_yyz_0_xyyy_1 = pbuffer.data(idx_eri_1_sfsg + 111);

    auto g_0_yyz_0_yyyy_1 = pbuffer.data(idx_eri_1_sfsg + 115);

    auto g_0_yzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sfsg + 120);

    auto g_0_yzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sfsg + 122);

    auto g_0_yzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sfsg + 124);

    auto g_0_yzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sfsg + 125);

    auto g_0_yzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sfsg + 127);

    auto g_0_yzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sfsg + 128);

    auto g_0_yzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sfsg + 129);

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

    /// Set up components of auxilary buffer : SGSF

    auto g_0_xxxx_0_xxx_1 = pbuffer.data(idx_eri_1_sgsf);

    auto g_0_xxxx_0_xxy_1 = pbuffer.data(idx_eri_1_sgsf + 1);

    auto g_0_xxxx_0_xxz_1 = pbuffer.data(idx_eri_1_sgsf + 2);

    auto g_0_xxxx_0_xyy_1 = pbuffer.data(idx_eri_1_sgsf + 3);

    auto g_0_xxxx_0_xyz_1 = pbuffer.data(idx_eri_1_sgsf + 4);

    auto g_0_xxxx_0_xzz_1 = pbuffer.data(idx_eri_1_sgsf + 5);

    auto g_0_xxxx_0_yyy_1 = pbuffer.data(idx_eri_1_sgsf + 6);

    auto g_0_xxxx_0_yyz_1 = pbuffer.data(idx_eri_1_sgsf + 7);

    auto g_0_xxxx_0_yzz_1 = pbuffer.data(idx_eri_1_sgsf + 8);

    auto g_0_xxxx_0_zzz_1 = pbuffer.data(idx_eri_1_sgsf + 9);

    auto g_0_xxxz_0_xxz_1 = pbuffer.data(idx_eri_1_sgsf + 22);

    auto g_0_xxxz_0_xyz_1 = pbuffer.data(idx_eri_1_sgsf + 24);

    auto g_0_xxxz_0_xzz_1 = pbuffer.data(idx_eri_1_sgsf + 25);

    auto g_0_xxxz_0_yyz_1 = pbuffer.data(idx_eri_1_sgsf + 27);

    auto g_0_xxxz_0_yzz_1 = pbuffer.data(idx_eri_1_sgsf + 28);

    auto g_0_xxxz_0_zzz_1 = pbuffer.data(idx_eri_1_sgsf + 29);

    auto g_0_xxyy_0_xxx_1 = pbuffer.data(idx_eri_1_sgsf + 30);

    auto g_0_xxyy_0_xxy_1 = pbuffer.data(idx_eri_1_sgsf + 31);

    auto g_0_xxyy_0_xxz_1 = pbuffer.data(idx_eri_1_sgsf + 32);

    auto g_0_xxyy_0_xyy_1 = pbuffer.data(idx_eri_1_sgsf + 33);

    auto g_0_xxyy_0_xyz_1 = pbuffer.data(idx_eri_1_sgsf + 34);

    auto g_0_xxyy_0_xzz_1 = pbuffer.data(idx_eri_1_sgsf + 35);

    auto g_0_xxyy_0_yyy_1 = pbuffer.data(idx_eri_1_sgsf + 36);

    auto g_0_xxyy_0_yyz_1 = pbuffer.data(idx_eri_1_sgsf + 37);

    auto g_0_xxyy_0_yzz_1 = pbuffer.data(idx_eri_1_sgsf + 38);

    auto g_0_xxyy_0_zzz_1 = pbuffer.data(idx_eri_1_sgsf + 39);

    auto g_0_xxzz_0_xxx_1 = pbuffer.data(idx_eri_1_sgsf + 50);

    auto g_0_xxzz_0_xxy_1 = pbuffer.data(idx_eri_1_sgsf + 51);

    auto g_0_xxzz_0_xxz_1 = pbuffer.data(idx_eri_1_sgsf + 52);

    auto g_0_xxzz_0_xyy_1 = pbuffer.data(idx_eri_1_sgsf + 53);

    auto g_0_xxzz_0_xyz_1 = pbuffer.data(idx_eri_1_sgsf + 54);

    auto g_0_xxzz_0_xzz_1 = pbuffer.data(idx_eri_1_sgsf + 55);

    auto g_0_xxzz_0_yyy_1 = pbuffer.data(idx_eri_1_sgsf + 56);

    auto g_0_xxzz_0_yyz_1 = pbuffer.data(idx_eri_1_sgsf + 57);

    auto g_0_xxzz_0_yzz_1 = pbuffer.data(idx_eri_1_sgsf + 58);

    auto g_0_xxzz_0_zzz_1 = pbuffer.data(idx_eri_1_sgsf + 59);

    auto g_0_xyyy_0_xxy_1 = pbuffer.data(idx_eri_1_sgsf + 61);

    auto g_0_xyyy_0_xyy_1 = pbuffer.data(idx_eri_1_sgsf + 63);

    auto g_0_xyyy_0_xyz_1 = pbuffer.data(idx_eri_1_sgsf + 64);

    auto g_0_xyyy_0_yyy_1 = pbuffer.data(idx_eri_1_sgsf + 66);

    auto g_0_xyyy_0_yyz_1 = pbuffer.data(idx_eri_1_sgsf + 67);

    auto g_0_xyyy_0_yzz_1 = pbuffer.data(idx_eri_1_sgsf + 68);

    auto g_0_xzzz_0_xxz_1 = pbuffer.data(idx_eri_1_sgsf + 92);

    auto g_0_xzzz_0_xyz_1 = pbuffer.data(idx_eri_1_sgsf + 94);

    auto g_0_xzzz_0_xzz_1 = pbuffer.data(idx_eri_1_sgsf + 95);

    auto g_0_xzzz_0_yyz_1 = pbuffer.data(idx_eri_1_sgsf + 97);

    auto g_0_xzzz_0_yzz_1 = pbuffer.data(idx_eri_1_sgsf + 98);

    auto g_0_xzzz_0_zzz_1 = pbuffer.data(idx_eri_1_sgsf + 99);

    auto g_0_yyyy_0_xxx_1 = pbuffer.data(idx_eri_1_sgsf + 100);

    auto g_0_yyyy_0_xxy_1 = pbuffer.data(idx_eri_1_sgsf + 101);

    auto g_0_yyyy_0_xxz_1 = pbuffer.data(idx_eri_1_sgsf + 102);

    auto g_0_yyyy_0_xyy_1 = pbuffer.data(idx_eri_1_sgsf + 103);

    auto g_0_yyyy_0_xyz_1 = pbuffer.data(idx_eri_1_sgsf + 104);

    auto g_0_yyyy_0_xzz_1 = pbuffer.data(idx_eri_1_sgsf + 105);

    auto g_0_yyyy_0_yyy_1 = pbuffer.data(idx_eri_1_sgsf + 106);

    auto g_0_yyyy_0_yyz_1 = pbuffer.data(idx_eri_1_sgsf + 107);

    auto g_0_yyyy_0_yzz_1 = pbuffer.data(idx_eri_1_sgsf + 108);

    auto g_0_yyyy_0_zzz_1 = pbuffer.data(idx_eri_1_sgsf + 109);

    auto g_0_yyyz_0_xxz_1 = pbuffer.data(idx_eri_1_sgsf + 112);

    auto g_0_yyyz_0_xyz_1 = pbuffer.data(idx_eri_1_sgsf + 114);

    auto g_0_yyyz_0_xzz_1 = pbuffer.data(idx_eri_1_sgsf + 115);

    auto g_0_yyyz_0_yyz_1 = pbuffer.data(idx_eri_1_sgsf + 117);

    auto g_0_yyyz_0_yzz_1 = pbuffer.data(idx_eri_1_sgsf + 118);

    auto g_0_yyyz_0_zzz_1 = pbuffer.data(idx_eri_1_sgsf + 119);

    auto g_0_yyzz_0_xxx_1 = pbuffer.data(idx_eri_1_sgsf + 120);

    auto g_0_yyzz_0_xxy_1 = pbuffer.data(idx_eri_1_sgsf + 121);

    auto g_0_yyzz_0_xxz_1 = pbuffer.data(idx_eri_1_sgsf + 122);

    auto g_0_yyzz_0_xyy_1 = pbuffer.data(idx_eri_1_sgsf + 123);

    auto g_0_yyzz_0_xyz_1 = pbuffer.data(idx_eri_1_sgsf + 124);

    auto g_0_yyzz_0_xzz_1 = pbuffer.data(idx_eri_1_sgsf + 125);

    auto g_0_yyzz_0_yyy_1 = pbuffer.data(idx_eri_1_sgsf + 126);

    auto g_0_yyzz_0_yyz_1 = pbuffer.data(idx_eri_1_sgsf + 127);

    auto g_0_yyzz_0_yzz_1 = pbuffer.data(idx_eri_1_sgsf + 128);

    auto g_0_yyzz_0_zzz_1 = pbuffer.data(idx_eri_1_sgsf + 129);

    auto g_0_yzzz_0_xxy_1 = pbuffer.data(idx_eri_1_sgsf + 131);

    auto g_0_yzzz_0_xxz_1 = pbuffer.data(idx_eri_1_sgsf + 132);

    auto g_0_yzzz_0_xyy_1 = pbuffer.data(idx_eri_1_sgsf + 133);

    auto g_0_yzzz_0_xyz_1 = pbuffer.data(idx_eri_1_sgsf + 134);

    auto g_0_yzzz_0_xzz_1 = pbuffer.data(idx_eri_1_sgsf + 135);

    auto g_0_yzzz_0_yyy_1 = pbuffer.data(idx_eri_1_sgsf + 136);

    auto g_0_yzzz_0_yyz_1 = pbuffer.data(idx_eri_1_sgsf + 137);

    auto g_0_yzzz_0_yzz_1 = pbuffer.data(idx_eri_1_sgsf + 138);

    auto g_0_yzzz_0_zzz_1 = pbuffer.data(idx_eri_1_sgsf + 139);

    auto g_0_zzzz_0_xxx_1 = pbuffer.data(idx_eri_1_sgsf + 140);

    auto g_0_zzzz_0_xxy_1 = pbuffer.data(idx_eri_1_sgsf + 141);

    auto g_0_zzzz_0_xxz_1 = pbuffer.data(idx_eri_1_sgsf + 142);

    auto g_0_zzzz_0_xyy_1 = pbuffer.data(idx_eri_1_sgsf + 143);

    auto g_0_zzzz_0_xyz_1 = pbuffer.data(idx_eri_1_sgsf + 144);

    auto g_0_zzzz_0_xzz_1 = pbuffer.data(idx_eri_1_sgsf + 145);

    auto g_0_zzzz_0_yyy_1 = pbuffer.data(idx_eri_1_sgsf + 146);

    auto g_0_zzzz_0_yyz_1 = pbuffer.data(idx_eri_1_sgsf + 147);

    auto g_0_zzzz_0_yzz_1 = pbuffer.data(idx_eri_1_sgsf + 148);

    auto g_0_zzzz_0_zzz_1 = pbuffer.data(idx_eri_1_sgsf + 149);

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

    auto g_0_xxxy_0_xxxy_0 = pbuffer.data(idx_eri_0_sgsg + 16);

    auto g_0_xxxy_0_xxxz_0 = pbuffer.data(idx_eri_0_sgsg + 17);

    auto g_0_xxxy_0_xxyy_0 = pbuffer.data(idx_eri_0_sgsg + 18);

    auto g_0_xxxy_0_xxzz_0 = pbuffer.data(idx_eri_0_sgsg + 20);

    auto g_0_xxxy_0_xyyy_0 = pbuffer.data(idx_eri_0_sgsg + 21);

    auto g_0_xxxy_0_xzzz_0 = pbuffer.data(idx_eri_0_sgsg + 24);

    auto g_0_xxxy_0_yyyy_0 = pbuffer.data(idx_eri_0_sgsg + 25);

    auto g_0_xxxz_0_xxxx_0 = pbuffer.data(idx_eri_0_sgsg + 30);

    auto g_0_xxxz_0_xxxy_0 = pbuffer.data(idx_eri_0_sgsg + 31);

    auto g_0_xxxz_0_xxxz_0 = pbuffer.data(idx_eri_0_sgsg + 32);

    auto g_0_xxxz_0_xxyy_0 = pbuffer.data(idx_eri_0_sgsg + 33);

    auto g_0_xxxz_0_xxyz_0 = pbuffer.data(idx_eri_0_sgsg + 34);

    auto g_0_xxxz_0_xxzz_0 = pbuffer.data(idx_eri_0_sgsg + 35);

    auto g_0_xxxz_0_xyyy_0 = pbuffer.data(idx_eri_0_sgsg + 36);

    auto g_0_xxxz_0_xyyz_0 = pbuffer.data(idx_eri_0_sgsg + 37);

    auto g_0_xxxz_0_xyzz_0 = pbuffer.data(idx_eri_0_sgsg + 38);

    auto g_0_xxxz_0_xzzz_0 = pbuffer.data(idx_eri_0_sgsg + 39);

    auto g_0_xxxz_0_yyyz_0 = pbuffer.data(idx_eri_0_sgsg + 41);

    auto g_0_xxxz_0_yyzz_0 = pbuffer.data(idx_eri_0_sgsg + 42);

    auto g_0_xxxz_0_yzzz_0 = pbuffer.data(idx_eri_0_sgsg + 43);

    auto g_0_xxxz_0_zzzz_0 = pbuffer.data(idx_eri_0_sgsg + 44);

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

    auto g_0_xyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_sgsg + 90);

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

    auto g_0_xzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sgsg + 135);

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

    auto g_0_yyyz_0_xxxz_0 = pbuffer.data(idx_eri_0_sgsg + 167);

    auto g_0_yyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_sgsg + 168);

    auto g_0_yyyz_0_xxyz_0 = pbuffer.data(idx_eri_0_sgsg + 169);

    auto g_0_yyyz_0_xxzz_0 = pbuffer.data(idx_eri_0_sgsg + 170);

    auto g_0_yyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_sgsg + 171);

    auto g_0_yyyz_0_xyyz_0 = pbuffer.data(idx_eri_0_sgsg + 172);

    auto g_0_yyyz_0_xyzz_0 = pbuffer.data(idx_eri_0_sgsg + 173);

    auto g_0_yyyz_0_xzzz_0 = pbuffer.data(idx_eri_0_sgsg + 174);

    auto g_0_yyyz_0_yyyy_0 = pbuffer.data(idx_eri_0_sgsg + 175);

    auto g_0_yyyz_0_yyyz_0 = pbuffer.data(idx_eri_0_sgsg + 176);

    auto g_0_yyyz_0_yyzz_0 = pbuffer.data(idx_eri_0_sgsg + 177);

    auto g_0_yyyz_0_yzzz_0 = pbuffer.data(idx_eri_0_sgsg + 178);

    auto g_0_yyyz_0_zzzz_0 = pbuffer.data(idx_eri_0_sgsg + 179);

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

    auto g_0_yzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sgsg + 196);

    auto g_0_yzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sgsg + 197);

    auto g_0_yzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sgsg + 198);

    auto g_0_yzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sgsg + 199);

    auto g_0_yzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sgsg + 200);

    auto g_0_yzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sgsg + 201);

    auto g_0_yzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sgsg + 202);

    auto g_0_yzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sgsg + 203);

    auto g_0_yzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sgsg + 204);

    auto g_0_yzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sgsg + 205);

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

    auto g_0_xxxy_0_xxxy_1 = pbuffer.data(idx_eri_1_sgsg + 16);

    auto g_0_xxxy_0_xxxz_1 = pbuffer.data(idx_eri_1_sgsg + 17);

    auto g_0_xxxy_0_xxyy_1 = pbuffer.data(idx_eri_1_sgsg + 18);

    auto g_0_xxxy_0_xxzz_1 = pbuffer.data(idx_eri_1_sgsg + 20);

    auto g_0_xxxy_0_xyyy_1 = pbuffer.data(idx_eri_1_sgsg + 21);

    auto g_0_xxxy_0_xzzz_1 = pbuffer.data(idx_eri_1_sgsg + 24);

    auto g_0_xxxy_0_yyyy_1 = pbuffer.data(idx_eri_1_sgsg + 25);

    auto g_0_xxxz_0_xxxx_1 = pbuffer.data(idx_eri_1_sgsg + 30);

    auto g_0_xxxz_0_xxxy_1 = pbuffer.data(idx_eri_1_sgsg + 31);

    auto g_0_xxxz_0_xxxz_1 = pbuffer.data(idx_eri_1_sgsg + 32);

    auto g_0_xxxz_0_xxyy_1 = pbuffer.data(idx_eri_1_sgsg + 33);

    auto g_0_xxxz_0_xxyz_1 = pbuffer.data(idx_eri_1_sgsg + 34);

    auto g_0_xxxz_0_xxzz_1 = pbuffer.data(idx_eri_1_sgsg + 35);

    auto g_0_xxxz_0_xyyy_1 = pbuffer.data(idx_eri_1_sgsg + 36);

    auto g_0_xxxz_0_xyyz_1 = pbuffer.data(idx_eri_1_sgsg + 37);

    auto g_0_xxxz_0_xyzz_1 = pbuffer.data(idx_eri_1_sgsg + 38);

    auto g_0_xxxz_0_xzzz_1 = pbuffer.data(idx_eri_1_sgsg + 39);

    auto g_0_xxxz_0_yyyz_1 = pbuffer.data(idx_eri_1_sgsg + 41);

    auto g_0_xxxz_0_yyzz_1 = pbuffer.data(idx_eri_1_sgsg + 42);

    auto g_0_xxxz_0_yzzz_1 = pbuffer.data(idx_eri_1_sgsg + 43);

    auto g_0_xxxz_0_zzzz_1 = pbuffer.data(idx_eri_1_sgsg + 44);

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

    auto g_0_xyyy_0_xxxx_1 = pbuffer.data(idx_eri_1_sgsg + 90);

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

    auto g_0_xzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sgsg + 135);

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

    auto g_0_yyyz_0_xxxz_1 = pbuffer.data(idx_eri_1_sgsg + 167);

    auto g_0_yyyz_0_xxyy_1 = pbuffer.data(idx_eri_1_sgsg + 168);

    auto g_0_yyyz_0_xxyz_1 = pbuffer.data(idx_eri_1_sgsg + 169);

    auto g_0_yyyz_0_xxzz_1 = pbuffer.data(idx_eri_1_sgsg + 170);

    auto g_0_yyyz_0_xyyy_1 = pbuffer.data(idx_eri_1_sgsg + 171);

    auto g_0_yyyz_0_xyyz_1 = pbuffer.data(idx_eri_1_sgsg + 172);

    auto g_0_yyyz_0_xyzz_1 = pbuffer.data(idx_eri_1_sgsg + 173);

    auto g_0_yyyz_0_xzzz_1 = pbuffer.data(idx_eri_1_sgsg + 174);

    auto g_0_yyyz_0_yyyy_1 = pbuffer.data(idx_eri_1_sgsg + 175);

    auto g_0_yyyz_0_yyyz_1 = pbuffer.data(idx_eri_1_sgsg + 176);

    auto g_0_yyyz_0_yyzz_1 = pbuffer.data(idx_eri_1_sgsg + 177);

    auto g_0_yyyz_0_yzzz_1 = pbuffer.data(idx_eri_1_sgsg + 178);

    auto g_0_yyyz_0_zzzz_1 = pbuffer.data(idx_eri_1_sgsg + 179);

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

    auto g_0_yzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_sgsg + 196);

    auto g_0_yzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_sgsg + 197);

    auto g_0_yzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_sgsg + 198);

    auto g_0_yzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_sgsg + 199);

    auto g_0_yzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_sgsg + 200);

    auto g_0_yzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_sgsg + 201);

    auto g_0_yzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_sgsg + 202);

    auto g_0_yzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_sgsg + 203);

    auto g_0_yzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_sgsg + 204);

    auto g_0_yzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_sgsg + 205);

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

    /// Set up 0-15 components of targeted buffer : SHSG

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

#pragma omp simd aligned(g_0_xxx_0_xxxx_0,       \
                             g_0_xxx_0_xxxx_1,   \
                             g_0_xxx_0_xxxy_0,   \
                             g_0_xxx_0_xxxy_1,   \
                             g_0_xxx_0_xxxz_0,   \
                             g_0_xxx_0_xxxz_1,   \
                             g_0_xxx_0_xxyy_0,   \
                             g_0_xxx_0_xxyy_1,   \
                             g_0_xxx_0_xxyz_0,   \
                             g_0_xxx_0_xxyz_1,   \
                             g_0_xxx_0_xxzz_0,   \
                             g_0_xxx_0_xxzz_1,   \
                             g_0_xxx_0_xyyy_0,   \
                             g_0_xxx_0_xyyy_1,   \
                             g_0_xxx_0_xyyz_0,   \
                             g_0_xxx_0_xyyz_1,   \
                             g_0_xxx_0_xyzz_0,   \
                             g_0_xxx_0_xyzz_1,   \
                             g_0_xxx_0_xzzz_0,   \
                             g_0_xxx_0_xzzz_1,   \
                             g_0_xxx_0_yyyy_0,   \
                             g_0_xxx_0_yyyy_1,   \
                             g_0_xxx_0_yyyz_0,   \
                             g_0_xxx_0_yyyz_1,   \
                             g_0_xxx_0_yyzz_0,   \
                             g_0_xxx_0_yyzz_1,   \
                             g_0_xxx_0_yzzz_0,   \
                             g_0_xxx_0_yzzz_1,   \
                             g_0_xxx_0_zzzz_0,   \
                             g_0_xxx_0_zzzz_1,   \
                             g_0_xxxx_0_xxx_1,   \
                             g_0_xxxx_0_xxxx_0,  \
                             g_0_xxxx_0_xxxx_1,  \
                             g_0_xxxx_0_xxxy_0,  \
                             g_0_xxxx_0_xxxy_1,  \
                             g_0_xxxx_0_xxxz_0,  \
                             g_0_xxxx_0_xxxz_1,  \
                             g_0_xxxx_0_xxy_1,   \
                             g_0_xxxx_0_xxyy_0,  \
                             g_0_xxxx_0_xxyy_1,  \
                             g_0_xxxx_0_xxyz_0,  \
                             g_0_xxxx_0_xxyz_1,  \
                             g_0_xxxx_0_xxz_1,   \
                             g_0_xxxx_0_xxzz_0,  \
                             g_0_xxxx_0_xxzz_1,  \
                             g_0_xxxx_0_xyy_1,   \
                             g_0_xxxx_0_xyyy_0,  \
                             g_0_xxxx_0_xyyy_1,  \
                             g_0_xxxx_0_xyyz_0,  \
                             g_0_xxxx_0_xyyz_1,  \
                             g_0_xxxx_0_xyz_1,   \
                             g_0_xxxx_0_xyzz_0,  \
                             g_0_xxxx_0_xyzz_1,  \
                             g_0_xxxx_0_xzz_1,   \
                             g_0_xxxx_0_xzzz_0,  \
                             g_0_xxxx_0_xzzz_1,  \
                             g_0_xxxx_0_yyy_1,   \
                             g_0_xxxx_0_yyyy_0,  \
                             g_0_xxxx_0_yyyy_1,  \
                             g_0_xxxx_0_yyyz_0,  \
                             g_0_xxxx_0_yyyz_1,  \
                             g_0_xxxx_0_yyz_1,   \
                             g_0_xxxx_0_yyzz_0,  \
                             g_0_xxxx_0_yyzz_1,  \
                             g_0_xxxx_0_yzz_1,   \
                             g_0_xxxx_0_yzzz_0,  \
                             g_0_xxxx_0_yzzz_1,  \
                             g_0_xxxx_0_zzz_1,   \
                             g_0_xxxx_0_zzzz_0,  \
                             g_0_xxxx_0_zzzz_1,  \
                             g_0_xxxxx_0_xxxx_0, \
                             g_0_xxxxx_0_xxxy_0, \
                             g_0_xxxxx_0_xxxz_0, \
                             g_0_xxxxx_0_xxyy_0, \
                             g_0_xxxxx_0_xxyz_0, \
                             g_0_xxxxx_0_xxzz_0, \
                             g_0_xxxxx_0_xyyy_0, \
                             g_0_xxxxx_0_xyyz_0, \
                             g_0_xxxxx_0_xyzz_0, \
                             g_0_xxxxx_0_xzzz_0, \
                             g_0_xxxxx_0_yyyy_0, \
                             g_0_xxxxx_0_yyyz_0, \
                             g_0_xxxxx_0_yyzz_0, \
                             g_0_xxxxx_0_yzzz_0, \
                             g_0_xxxxx_0_zzzz_0, \
                             wp_x,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxx_0_xxxx_0[i] = 4.0 * g_0_xxx_0_xxxx_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxx_1[i] * fti_ab_0 + 4.0 * g_0_xxxx_0_xxx_1[i] * fi_abcd_0 +
                                g_0_xxxx_0_xxxx_0[i] * pb_x + g_0_xxxx_0_xxxx_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxy_0[i] = 4.0 * g_0_xxx_0_xxxy_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxy_1[i] * fti_ab_0 + 3.0 * g_0_xxxx_0_xxy_1[i] * fi_abcd_0 +
                                g_0_xxxx_0_xxxy_0[i] * pb_x + g_0_xxxx_0_xxxy_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxz_0[i] = 4.0 * g_0_xxx_0_xxxz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxz_1[i] * fti_ab_0 + 3.0 * g_0_xxxx_0_xxz_1[i] * fi_abcd_0 +
                                g_0_xxxx_0_xxxz_0[i] * pb_x + g_0_xxxx_0_xxxz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxyy_0[i] = 4.0 * g_0_xxx_0_xxyy_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxyy_1[i] * fti_ab_0 + 2.0 * g_0_xxxx_0_xyy_1[i] * fi_abcd_0 +
                                g_0_xxxx_0_xxyy_0[i] * pb_x + g_0_xxxx_0_xxyy_1[i] * wp_x[i];

        g_0_xxxxx_0_xxyz_0[i] = 4.0 * g_0_xxx_0_xxyz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxx_0_xyz_1[i] * fi_abcd_0 +
                                g_0_xxxx_0_xxyz_0[i] * pb_x + g_0_xxxx_0_xxyz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxzz_0[i] = 4.0 * g_0_xxx_0_xxzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxx_0_xzz_1[i] * fi_abcd_0 +
                                g_0_xxxx_0_xxzz_0[i] * pb_x + g_0_xxxx_0_xxzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xyyy_0[i] = 4.0 * g_0_xxx_0_xyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xyyy_1[i] * fti_ab_0 + g_0_xxxx_0_yyy_1[i] * fi_abcd_0 +
                                g_0_xxxx_0_xyyy_0[i] * pb_x + g_0_xxxx_0_xyyy_1[i] * wp_x[i];

        g_0_xxxxx_0_xyyz_0[i] = 4.0 * g_0_xxx_0_xyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xyyz_1[i] * fti_ab_0 + g_0_xxxx_0_yyz_1[i] * fi_abcd_0 +
                                g_0_xxxx_0_xyyz_0[i] * pb_x + g_0_xxxx_0_xyyz_1[i] * wp_x[i];

        g_0_xxxxx_0_xyzz_0[i] = 4.0 * g_0_xxx_0_xyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xyzz_1[i] * fti_ab_0 + g_0_xxxx_0_yzz_1[i] * fi_abcd_0 +
                                g_0_xxxx_0_xyzz_0[i] * pb_x + g_0_xxxx_0_xyzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xzzz_0[i] = 4.0 * g_0_xxx_0_xzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xzzz_1[i] * fti_ab_0 + g_0_xxxx_0_zzz_1[i] * fi_abcd_0 +
                                g_0_xxxx_0_xzzz_0[i] * pb_x + g_0_xxxx_0_xzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_yyyy_0[i] =
            4.0 * g_0_xxx_0_yyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_yyyy_1[i] * fti_ab_0 + g_0_xxxx_0_yyyy_0[i] * pb_x + g_0_xxxx_0_yyyy_1[i] * wp_x[i];

        g_0_xxxxx_0_yyyz_0[i] =
            4.0 * g_0_xxx_0_yyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_yyyz_1[i] * fti_ab_0 + g_0_xxxx_0_yyyz_0[i] * pb_x + g_0_xxxx_0_yyyz_1[i] * wp_x[i];

        g_0_xxxxx_0_yyzz_0[i] =
            4.0 * g_0_xxx_0_yyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_yyzz_1[i] * fti_ab_0 + g_0_xxxx_0_yyzz_0[i] * pb_x + g_0_xxxx_0_yyzz_1[i] * wp_x[i];

        g_0_xxxxx_0_yzzz_0[i] =
            4.0 * g_0_xxx_0_yzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_yzzz_1[i] * fti_ab_0 + g_0_xxxx_0_yzzz_0[i] * pb_x + g_0_xxxx_0_yzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_zzzz_0[i] =
            4.0 * g_0_xxx_0_zzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_zzzz_1[i] * fti_ab_0 + g_0_xxxx_0_zzzz_0[i] * pb_x + g_0_xxxx_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 15-30 components of targeted buffer : SHSG

    auto g_0_xxxxy_0_xxxx_0 = pbuffer.data(idx_eri_0_shsg + 15);

    auto g_0_xxxxy_0_xxxy_0 = pbuffer.data(idx_eri_0_shsg + 16);

    auto g_0_xxxxy_0_xxxz_0 = pbuffer.data(idx_eri_0_shsg + 17);

    auto g_0_xxxxy_0_xxyy_0 = pbuffer.data(idx_eri_0_shsg + 18);

    auto g_0_xxxxy_0_xxyz_0 = pbuffer.data(idx_eri_0_shsg + 19);

    auto g_0_xxxxy_0_xxzz_0 = pbuffer.data(idx_eri_0_shsg + 20);

    auto g_0_xxxxy_0_xyyy_0 = pbuffer.data(idx_eri_0_shsg + 21);

    auto g_0_xxxxy_0_xyyz_0 = pbuffer.data(idx_eri_0_shsg + 22);

    auto g_0_xxxxy_0_xyzz_0 = pbuffer.data(idx_eri_0_shsg + 23);

    auto g_0_xxxxy_0_xzzz_0 = pbuffer.data(idx_eri_0_shsg + 24);

    auto g_0_xxxxy_0_yyyy_0 = pbuffer.data(idx_eri_0_shsg + 25);

    auto g_0_xxxxy_0_yyyz_0 = pbuffer.data(idx_eri_0_shsg + 26);

    auto g_0_xxxxy_0_yyzz_0 = pbuffer.data(idx_eri_0_shsg + 27);

    auto g_0_xxxxy_0_yzzz_0 = pbuffer.data(idx_eri_0_shsg + 28);

    auto g_0_xxxxy_0_zzzz_0 = pbuffer.data(idx_eri_0_shsg + 29);

#pragma omp simd aligned(g_0_xxxx_0_xxx_1,       \
                             g_0_xxxx_0_xxxx_0,  \
                             g_0_xxxx_0_xxxx_1,  \
                             g_0_xxxx_0_xxxy_0,  \
                             g_0_xxxx_0_xxxy_1,  \
                             g_0_xxxx_0_xxxz_0,  \
                             g_0_xxxx_0_xxxz_1,  \
                             g_0_xxxx_0_xxy_1,   \
                             g_0_xxxx_0_xxyy_0,  \
                             g_0_xxxx_0_xxyy_1,  \
                             g_0_xxxx_0_xxyz_0,  \
                             g_0_xxxx_0_xxyz_1,  \
                             g_0_xxxx_0_xxz_1,   \
                             g_0_xxxx_0_xxzz_0,  \
                             g_0_xxxx_0_xxzz_1,  \
                             g_0_xxxx_0_xyy_1,   \
                             g_0_xxxx_0_xyyy_0,  \
                             g_0_xxxx_0_xyyy_1,  \
                             g_0_xxxx_0_xyyz_0,  \
                             g_0_xxxx_0_xyyz_1,  \
                             g_0_xxxx_0_xyz_1,   \
                             g_0_xxxx_0_xyzz_0,  \
                             g_0_xxxx_0_xyzz_1,  \
                             g_0_xxxx_0_xzz_1,   \
                             g_0_xxxx_0_xzzz_0,  \
                             g_0_xxxx_0_xzzz_1,  \
                             g_0_xxxx_0_yyy_1,   \
                             g_0_xxxx_0_yyyy_0,  \
                             g_0_xxxx_0_yyyy_1,  \
                             g_0_xxxx_0_yyyz_0,  \
                             g_0_xxxx_0_yyyz_1,  \
                             g_0_xxxx_0_yyz_1,   \
                             g_0_xxxx_0_yyzz_0,  \
                             g_0_xxxx_0_yyzz_1,  \
                             g_0_xxxx_0_yzz_1,   \
                             g_0_xxxx_0_yzzz_0,  \
                             g_0_xxxx_0_yzzz_1,  \
                             g_0_xxxx_0_zzz_1,   \
                             g_0_xxxx_0_zzzz_0,  \
                             g_0_xxxx_0_zzzz_1,  \
                             g_0_xxxxy_0_xxxx_0, \
                             g_0_xxxxy_0_xxxy_0, \
                             g_0_xxxxy_0_xxxz_0, \
                             g_0_xxxxy_0_xxyy_0, \
                             g_0_xxxxy_0_xxyz_0, \
                             g_0_xxxxy_0_xxzz_0, \
                             g_0_xxxxy_0_xyyy_0, \
                             g_0_xxxxy_0_xyyz_0, \
                             g_0_xxxxy_0_xyzz_0, \
                             g_0_xxxxy_0_xzzz_0, \
                             g_0_xxxxy_0_yyyy_0, \
                             g_0_xxxxy_0_yyyz_0, \
                             g_0_xxxxy_0_yyzz_0, \
                             g_0_xxxxy_0_yzzz_0, \
                             g_0_xxxxy_0_zzzz_0, \
                             wp_y,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxy_0_xxxx_0[i] = g_0_xxxx_0_xxxx_0[i] * pb_y + g_0_xxxx_0_xxxx_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxy_0[i] = g_0_xxxx_0_xxx_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxy_0[i] * pb_y + g_0_xxxx_0_xxxy_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxz_0[i] = g_0_xxxx_0_xxxz_0[i] * pb_y + g_0_xxxx_0_xxxz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxyy_0[i] = 2.0 * g_0_xxxx_0_xxy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyy_0[i] * pb_y + g_0_xxxx_0_xxyy_1[i] * wp_y[i];

        g_0_xxxxy_0_xxyz_0[i] = g_0_xxxx_0_xxz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyz_0[i] * pb_y + g_0_xxxx_0_xxyz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxzz_0[i] = g_0_xxxx_0_xxzz_0[i] * pb_y + g_0_xxxx_0_xxzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xyyy_0[i] = 3.0 * g_0_xxxx_0_xyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyy_0[i] * pb_y + g_0_xxxx_0_xyyy_1[i] * wp_y[i];

        g_0_xxxxy_0_xyyz_0[i] = 2.0 * g_0_xxxx_0_xyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyz_0[i] * pb_y + g_0_xxxx_0_xyyz_1[i] * wp_y[i];

        g_0_xxxxy_0_xyzz_0[i] = g_0_xxxx_0_xzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyzz_0[i] * pb_y + g_0_xxxx_0_xyzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xzzz_0[i] = g_0_xxxx_0_xzzz_0[i] * pb_y + g_0_xxxx_0_xzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_yyyy_0[i] = 4.0 * g_0_xxxx_0_yyy_1[i] * fi_abcd_0 + g_0_xxxx_0_yyyy_0[i] * pb_y + g_0_xxxx_0_yyyy_1[i] * wp_y[i];

        g_0_xxxxy_0_yyyz_0[i] = 3.0 * g_0_xxxx_0_yyz_1[i] * fi_abcd_0 + g_0_xxxx_0_yyyz_0[i] * pb_y + g_0_xxxx_0_yyyz_1[i] * wp_y[i];

        g_0_xxxxy_0_yyzz_0[i] = 2.0 * g_0_xxxx_0_yzz_1[i] * fi_abcd_0 + g_0_xxxx_0_yyzz_0[i] * pb_y + g_0_xxxx_0_yyzz_1[i] * wp_y[i];

        g_0_xxxxy_0_yzzz_0[i] = g_0_xxxx_0_zzz_1[i] * fi_abcd_0 + g_0_xxxx_0_yzzz_0[i] * pb_y + g_0_xxxx_0_yzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_zzzz_0[i] = g_0_xxxx_0_zzzz_0[i] * pb_y + g_0_xxxx_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 30-45 components of targeted buffer : SHSG

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

    auto g_0_xxxxz_0_yyyy_0 = pbuffer.data(idx_eri_0_shsg + 40);

    auto g_0_xxxxz_0_yyyz_0 = pbuffer.data(idx_eri_0_shsg + 41);

    auto g_0_xxxxz_0_yyzz_0 = pbuffer.data(idx_eri_0_shsg + 42);

    auto g_0_xxxxz_0_yzzz_0 = pbuffer.data(idx_eri_0_shsg + 43);

    auto g_0_xxxxz_0_zzzz_0 = pbuffer.data(idx_eri_0_shsg + 44);

#pragma omp simd aligned(g_0_xxxx_0_xxx_1,       \
                             g_0_xxxx_0_xxxx_0,  \
                             g_0_xxxx_0_xxxx_1,  \
                             g_0_xxxx_0_xxxy_0,  \
                             g_0_xxxx_0_xxxy_1,  \
                             g_0_xxxx_0_xxxz_0,  \
                             g_0_xxxx_0_xxxz_1,  \
                             g_0_xxxx_0_xxy_1,   \
                             g_0_xxxx_0_xxyy_0,  \
                             g_0_xxxx_0_xxyy_1,  \
                             g_0_xxxx_0_xxyz_0,  \
                             g_0_xxxx_0_xxyz_1,  \
                             g_0_xxxx_0_xxz_1,   \
                             g_0_xxxx_0_xxzz_0,  \
                             g_0_xxxx_0_xxzz_1,  \
                             g_0_xxxx_0_xyy_1,   \
                             g_0_xxxx_0_xyyy_0,  \
                             g_0_xxxx_0_xyyy_1,  \
                             g_0_xxxx_0_xyyz_0,  \
                             g_0_xxxx_0_xyyz_1,  \
                             g_0_xxxx_0_xyz_1,   \
                             g_0_xxxx_0_xyzz_0,  \
                             g_0_xxxx_0_xyzz_1,  \
                             g_0_xxxx_0_xzz_1,   \
                             g_0_xxxx_0_xzzz_0,  \
                             g_0_xxxx_0_xzzz_1,  \
                             g_0_xxxx_0_yyy_1,   \
                             g_0_xxxx_0_yyyy_0,  \
                             g_0_xxxx_0_yyyy_1,  \
                             g_0_xxxx_0_yyyz_0,  \
                             g_0_xxxx_0_yyyz_1,  \
                             g_0_xxxx_0_yyz_1,   \
                             g_0_xxxx_0_yyzz_0,  \
                             g_0_xxxx_0_yyzz_1,  \
                             g_0_xxxx_0_yzz_1,   \
                             g_0_xxxx_0_yzzz_0,  \
                             g_0_xxxx_0_yzzz_1,  \
                             g_0_xxxx_0_zzz_1,   \
                             g_0_xxxx_0_zzzz_0,  \
                             g_0_xxxx_0_zzzz_1,  \
                             g_0_xxxxz_0_xxxx_0, \
                             g_0_xxxxz_0_xxxy_0, \
                             g_0_xxxxz_0_xxxz_0, \
                             g_0_xxxxz_0_xxyy_0, \
                             g_0_xxxxz_0_xxyz_0, \
                             g_0_xxxxz_0_xxzz_0, \
                             g_0_xxxxz_0_xyyy_0, \
                             g_0_xxxxz_0_xyyz_0, \
                             g_0_xxxxz_0_xyzz_0, \
                             g_0_xxxxz_0_xzzz_0, \
                             g_0_xxxxz_0_yyyy_0, \
                             g_0_xxxxz_0_yyyz_0, \
                             g_0_xxxxz_0_yyzz_0, \
                             g_0_xxxxz_0_yzzz_0, \
                             g_0_xxxxz_0_zzzz_0, \
                             wp_z,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxz_0_xxxx_0[i] = g_0_xxxx_0_xxxx_0[i] * pb_z + g_0_xxxx_0_xxxx_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxy_0[i] = g_0_xxxx_0_xxxy_0[i] * pb_z + g_0_xxxx_0_xxxy_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxz_0[i] = g_0_xxxx_0_xxx_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxz_0[i] * pb_z + g_0_xxxx_0_xxxz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxyy_0[i] = g_0_xxxx_0_xxyy_0[i] * pb_z + g_0_xxxx_0_xxyy_1[i] * wp_z[i];

        g_0_xxxxz_0_xxyz_0[i] = g_0_xxxx_0_xxy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyz_0[i] * pb_z + g_0_xxxx_0_xxyz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxzz_0[i] = 2.0 * g_0_xxxx_0_xxz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxzz_0[i] * pb_z + g_0_xxxx_0_xxzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xyyy_0[i] = g_0_xxxx_0_xyyy_0[i] * pb_z + g_0_xxxx_0_xyyy_1[i] * wp_z[i];

        g_0_xxxxz_0_xyyz_0[i] = g_0_xxxx_0_xyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyz_0[i] * pb_z + g_0_xxxx_0_xyyz_1[i] * wp_z[i];

        g_0_xxxxz_0_xyzz_0[i] = 2.0 * g_0_xxxx_0_xyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyzz_0[i] * pb_z + g_0_xxxx_0_xyzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xzzz_0[i] = 3.0 * g_0_xxxx_0_xzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xzzz_0[i] * pb_z + g_0_xxxx_0_xzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_yyyy_0[i] = g_0_xxxx_0_yyyy_0[i] * pb_z + g_0_xxxx_0_yyyy_1[i] * wp_z[i];

        g_0_xxxxz_0_yyyz_0[i] = g_0_xxxx_0_yyy_1[i] * fi_abcd_0 + g_0_xxxx_0_yyyz_0[i] * pb_z + g_0_xxxx_0_yyyz_1[i] * wp_z[i];

        g_0_xxxxz_0_yyzz_0[i] = 2.0 * g_0_xxxx_0_yyz_1[i] * fi_abcd_0 + g_0_xxxx_0_yyzz_0[i] * pb_z + g_0_xxxx_0_yyzz_1[i] * wp_z[i];

        g_0_xxxxz_0_yzzz_0[i] = 3.0 * g_0_xxxx_0_yzz_1[i] * fi_abcd_0 + g_0_xxxx_0_yzzz_0[i] * pb_z + g_0_xxxx_0_yzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_zzzz_0[i] = 4.0 * g_0_xxxx_0_zzz_1[i] * fi_abcd_0 + g_0_xxxx_0_zzzz_0[i] * pb_z + g_0_xxxx_0_zzzz_1[i] * wp_z[i];
    }

    /// Set up 45-60 components of targeted buffer : SHSG

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

#pragma omp simd aligned(g_0_xxx_0_xxxx_0,       \
                             g_0_xxx_0_xxxx_1,   \
                             g_0_xxx_0_xxxz_0,   \
                             g_0_xxx_0_xxxz_1,   \
                             g_0_xxx_0_xxzz_0,   \
                             g_0_xxx_0_xxzz_1,   \
                             g_0_xxx_0_xzzz_0,   \
                             g_0_xxx_0_xzzz_1,   \
                             g_0_xxxy_0_xxxx_0,  \
                             g_0_xxxy_0_xxxx_1,  \
                             g_0_xxxy_0_xxxz_0,  \
                             g_0_xxxy_0_xxxz_1,  \
                             g_0_xxxy_0_xxzz_0,  \
                             g_0_xxxy_0_xxzz_1,  \
                             g_0_xxxy_0_xzzz_0,  \
                             g_0_xxxy_0_xzzz_1,  \
                             g_0_xxxyy_0_xxxx_0, \
                             g_0_xxxyy_0_xxxy_0, \
                             g_0_xxxyy_0_xxxz_0, \
                             g_0_xxxyy_0_xxyy_0, \
                             g_0_xxxyy_0_xxyz_0, \
                             g_0_xxxyy_0_xxzz_0, \
                             g_0_xxxyy_0_xyyy_0, \
                             g_0_xxxyy_0_xyyz_0, \
                             g_0_xxxyy_0_xyzz_0, \
                             g_0_xxxyy_0_xzzz_0, \
                             g_0_xxxyy_0_yyyy_0, \
                             g_0_xxxyy_0_yyyz_0, \
                             g_0_xxxyy_0_yyzz_0, \
                             g_0_xxxyy_0_yzzz_0, \
                             g_0_xxxyy_0_zzzz_0, \
                             g_0_xxyy_0_xxxy_0,  \
                             g_0_xxyy_0_xxxy_1,  \
                             g_0_xxyy_0_xxy_1,   \
                             g_0_xxyy_0_xxyy_0,  \
                             g_0_xxyy_0_xxyy_1,  \
                             g_0_xxyy_0_xxyz_0,  \
                             g_0_xxyy_0_xxyz_1,  \
                             g_0_xxyy_0_xyy_1,   \
                             g_0_xxyy_0_xyyy_0,  \
                             g_0_xxyy_0_xyyy_1,  \
                             g_0_xxyy_0_xyyz_0,  \
                             g_0_xxyy_0_xyyz_1,  \
                             g_0_xxyy_0_xyz_1,   \
                             g_0_xxyy_0_xyzz_0,  \
                             g_0_xxyy_0_xyzz_1,  \
                             g_0_xxyy_0_yyy_1,   \
                             g_0_xxyy_0_yyyy_0,  \
                             g_0_xxyy_0_yyyy_1,  \
                             g_0_xxyy_0_yyyz_0,  \
                             g_0_xxyy_0_yyyz_1,  \
                             g_0_xxyy_0_yyz_1,   \
                             g_0_xxyy_0_yyzz_0,  \
                             g_0_xxyy_0_yyzz_1,  \
                             g_0_xxyy_0_yzz_1,   \
                             g_0_xxyy_0_yzzz_0,  \
                             g_0_xxyy_0_yzzz_1,  \
                             g_0_xxyy_0_zzzz_0,  \
                             g_0_xxyy_0_zzzz_1,  \
                             g_0_xyy_0_xxxy_0,   \
                             g_0_xyy_0_xxxy_1,   \
                             g_0_xyy_0_xxyy_0,   \
                             g_0_xyy_0_xxyy_1,   \
                             g_0_xyy_0_xxyz_0,   \
                             g_0_xyy_0_xxyz_1,   \
                             g_0_xyy_0_xyyy_0,   \
                             g_0_xyy_0_xyyy_1,   \
                             g_0_xyy_0_xyyz_0,   \
                             g_0_xyy_0_xyyz_1,   \
                             g_0_xyy_0_xyzz_0,   \
                             g_0_xyy_0_xyzz_1,   \
                             g_0_xyy_0_yyyy_0,   \
                             g_0_xyy_0_yyyy_1,   \
                             g_0_xyy_0_yyyz_0,   \
                             g_0_xyy_0_yyyz_1,   \
                             g_0_xyy_0_yyzz_0,   \
                             g_0_xyy_0_yyzz_1,   \
                             g_0_xyy_0_yzzz_0,   \
                             g_0_xyy_0_yzzz_1,   \
                             g_0_xyy_0_zzzz_0,   \
                             g_0_xyy_0_zzzz_1,   \
                             wp_x,               \
                             wp_y,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxyy_0_xxxx_0[i] =
            g_0_xxx_0_xxxx_0[i] * fi_ab_0 - g_0_xxx_0_xxxx_1[i] * fti_ab_0 + g_0_xxxy_0_xxxx_0[i] * pb_y + g_0_xxxy_0_xxxx_1[i] * wp_y[i];

        g_0_xxxyy_0_xxxy_0[i] = 2.0 * g_0_xyy_0_xxxy_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxxy_1[i] * fti_ab_0 + 3.0 * g_0_xxyy_0_xxy_1[i] * fi_abcd_0 +
                                g_0_xxyy_0_xxxy_0[i] * pb_x + g_0_xxyy_0_xxxy_1[i] * wp_x[i];

        g_0_xxxyy_0_xxxz_0[i] =
            g_0_xxx_0_xxxz_0[i] * fi_ab_0 - g_0_xxx_0_xxxz_1[i] * fti_ab_0 + g_0_xxxy_0_xxxz_0[i] * pb_y + g_0_xxxy_0_xxxz_1[i] * wp_y[i];

        g_0_xxxyy_0_xxyy_0[i] = 2.0 * g_0_xyy_0_xxyy_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxyy_1[i] * fti_ab_0 + 2.0 * g_0_xxyy_0_xyy_1[i] * fi_abcd_0 +
                                g_0_xxyy_0_xxyy_0[i] * pb_x + g_0_xxyy_0_xxyy_1[i] * wp_x[i];

        g_0_xxxyy_0_xxyz_0[i] = 2.0 * g_0_xyy_0_xxyz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxyz_1[i] * fti_ab_0 + 2.0 * g_0_xxyy_0_xyz_1[i] * fi_abcd_0 +
                                g_0_xxyy_0_xxyz_0[i] * pb_x + g_0_xxyy_0_xxyz_1[i] * wp_x[i];

        g_0_xxxyy_0_xxzz_0[i] =
            g_0_xxx_0_xxzz_0[i] * fi_ab_0 - g_0_xxx_0_xxzz_1[i] * fti_ab_0 + g_0_xxxy_0_xxzz_0[i] * pb_y + g_0_xxxy_0_xxzz_1[i] * wp_y[i];

        g_0_xxxyy_0_xyyy_0[i] = 2.0 * g_0_xyy_0_xyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xyyy_1[i] * fti_ab_0 + g_0_xxyy_0_yyy_1[i] * fi_abcd_0 +
                                g_0_xxyy_0_xyyy_0[i] * pb_x + g_0_xxyy_0_xyyy_1[i] * wp_x[i];

        g_0_xxxyy_0_xyyz_0[i] = 2.0 * g_0_xyy_0_xyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xyyz_1[i] * fti_ab_0 + g_0_xxyy_0_yyz_1[i] * fi_abcd_0 +
                                g_0_xxyy_0_xyyz_0[i] * pb_x + g_0_xxyy_0_xyyz_1[i] * wp_x[i];

        g_0_xxxyy_0_xyzz_0[i] = 2.0 * g_0_xyy_0_xyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xyzz_1[i] * fti_ab_0 + g_0_xxyy_0_yzz_1[i] * fi_abcd_0 +
                                g_0_xxyy_0_xyzz_0[i] * pb_x + g_0_xxyy_0_xyzz_1[i] * wp_x[i];

        g_0_xxxyy_0_xzzz_0[i] =
            g_0_xxx_0_xzzz_0[i] * fi_ab_0 - g_0_xxx_0_xzzz_1[i] * fti_ab_0 + g_0_xxxy_0_xzzz_0[i] * pb_y + g_0_xxxy_0_xzzz_1[i] * wp_y[i];

        g_0_xxxyy_0_yyyy_0[i] =
            2.0 * g_0_xyy_0_yyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_yyyy_1[i] * fti_ab_0 + g_0_xxyy_0_yyyy_0[i] * pb_x + g_0_xxyy_0_yyyy_1[i] * wp_x[i];

        g_0_xxxyy_0_yyyz_0[i] =
            2.0 * g_0_xyy_0_yyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_yyyz_1[i] * fti_ab_0 + g_0_xxyy_0_yyyz_0[i] * pb_x + g_0_xxyy_0_yyyz_1[i] * wp_x[i];

        g_0_xxxyy_0_yyzz_0[i] =
            2.0 * g_0_xyy_0_yyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_yyzz_1[i] * fti_ab_0 + g_0_xxyy_0_yyzz_0[i] * pb_x + g_0_xxyy_0_yyzz_1[i] * wp_x[i];

        g_0_xxxyy_0_yzzz_0[i] =
            2.0 * g_0_xyy_0_yzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_yzzz_1[i] * fti_ab_0 + g_0_xxyy_0_yzzz_0[i] * pb_x + g_0_xxyy_0_yzzz_1[i] * wp_x[i];

        g_0_xxxyy_0_zzzz_0[i] =
            2.0 * g_0_xyy_0_zzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_zzzz_1[i] * fti_ab_0 + g_0_xxyy_0_zzzz_0[i] * pb_x + g_0_xxyy_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 60-75 components of targeted buffer : SHSG

    auto g_0_xxxyz_0_xxxx_0 = pbuffer.data(idx_eri_0_shsg + 60);

    auto g_0_xxxyz_0_xxxy_0 = pbuffer.data(idx_eri_0_shsg + 61);

    auto g_0_xxxyz_0_xxxz_0 = pbuffer.data(idx_eri_0_shsg + 62);

    auto g_0_xxxyz_0_xxyy_0 = pbuffer.data(idx_eri_0_shsg + 63);

    auto g_0_xxxyz_0_xxyz_0 = pbuffer.data(idx_eri_0_shsg + 64);

    auto g_0_xxxyz_0_xxzz_0 = pbuffer.data(idx_eri_0_shsg + 65);

    auto g_0_xxxyz_0_xyyy_0 = pbuffer.data(idx_eri_0_shsg + 66);

    auto g_0_xxxyz_0_xyyz_0 = pbuffer.data(idx_eri_0_shsg + 67);

    auto g_0_xxxyz_0_xyzz_0 = pbuffer.data(idx_eri_0_shsg + 68);

    auto g_0_xxxyz_0_xzzz_0 = pbuffer.data(idx_eri_0_shsg + 69);

    auto g_0_xxxyz_0_yyyy_0 = pbuffer.data(idx_eri_0_shsg + 70);

    auto g_0_xxxyz_0_yyyz_0 = pbuffer.data(idx_eri_0_shsg + 71);

    auto g_0_xxxyz_0_yyzz_0 = pbuffer.data(idx_eri_0_shsg + 72);

    auto g_0_xxxyz_0_yzzz_0 = pbuffer.data(idx_eri_0_shsg + 73);

    auto g_0_xxxyz_0_zzzz_0 = pbuffer.data(idx_eri_0_shsg + 74);

#pragma omp simd aligned(g_0_xxxy_0_xxxy_0,      \
                             g_0_xxxy_0_xxxy_1,  \
                             g_0_xxxy_0_xxyy_0,  \
                             g_0_xxxy_0_xxyy_1,  \
                             g_0_xxxy_0_xyyy_0,  \
                             g_0_xxxy_0_xyyy_1,  \
                             g_0_xxxy_0_yyyy_0,  \
                             g_0_xxxy_0_yyyy_1,  \
                             g_0_xxxyz_0_xxxx_0, \
                             g_0_xxxyz_0_xxxy_0, \
                             g_0_xxxyz_0_xxxz_0, \
                             g_0_xxxyz_0_xxyy_0, \
                             g_0_xxxyz_0_xxyz_0, \
                             g_0_xxxyz_0_xxzz_0, \
                             g_0_xxxyz_0_xyyy_0, \
                             g_0_xxxyz_0_xyyz_0, \
                             g_0_xxxyz_0_xyzz_0, \
                             g_0_xxxyz_0_xzzz_0, \
                             g_0_xxxyz_0_yyyy_0, \
                             g_0_xxxyz_0_yyyz_0, \
                             g_0_xxxyz_0_yyzz_0, \
                             g_0_xxxyz_0_yzzz_0, \
                             g_0_xxxyz_0_zzzz_0, \
                             g_0_xxxz_0_xxxx_0,  \
                             g_0_xxxz_0_xxxx_1,  \
                             g_0_xxxz_0_xxxz_0,  \
                             g_0_xxxz_0_xxxz_1,  \
                             g_0_xxxz_0_xxyz_0,  \
                             g_0_xxxz_0_xxyz_1,  \
                             g_0_xxxz_0_xxz_1,   \
                             g_0_xxxz_0_xxzz_0,  \
                             g_0_xxxz_0_xxzz_1,  \
                             g_0_xxxz_0_xyyz_0,  \
                             g_0_xxxz_0_xyyz_1,  \
                             g_0_xxxz_0_xyz_1,   \
                             g_0_xxxz_0_xyzz_0,  \
                             g_0_xxxz_0_xyzz_1,  \
                             g_0_xxxz_0_xzz_1,   \
                             g_0_xxxz_0_xzzz_0,  \
                             g_0_xxxz_0_xzzz_1,  \
                             g_0_xxxz_0_yyyz_0,  \
                             g_0_xxxz_0_yyyz_1,  \
                             g_0_xxxz_0_yyz_1,   \
                             g_0_xxxz_0_yyzz_0,  \
                             g_0_xxxz_0_yyzz_1,  \
                             g_0_xxxz_0_yzz_1,   \
                             g_0_xxxz_0_yzzz_0,  \
                             g_0_xxxz_0_yzzz_1,  \
                             g_0_xxxz_0_zzz_1,   \
                             g_0_xxxz_0_zzzz_0,  \
                             g_0_xxxz_0_zzzz_1,  \
                             wp_y,               \
                             wp_z,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyz_0_xxxx_0[i] = g_0_xxxz_0_xxxx_0[i] * pb_y + g_0_xxxz_0_xxxx_1[i] * wp_y[i];

        g_0_xxxyz_0_xxxy_0[i] = g_0_xxxy_0_xxxy_0[i] * pb_z + g_0_xxxy_0_xxxy_1[i] * wp_z[i];

        g_0_xxxyz_0_xxxz_0[i] = g_0_xxxz_0_xxxz_0[i] * pb_y + g_0_xxxz_0_xxxz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxyy_0[i] = g_0_xxxy_0_xxyy_0[i] * pb_z + g_0_xxxy_0_xxyy_1[i] * wp_z[i];

        g_0_xxxyz_0_xxyz_0[i] = g_0_xxxz_0_xxz_1[i] * fi_abcd_0 + g_0_xxxz_0_xxyz_0[i] * pb_y + g_0_xxxz_0_xxyz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxzz_0[i] = g_0_xxxz_0_xxzz_0[i] * pb_y + g_0_xxxz_0_xxzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xyyy_0[i] = g_0_xxxy_0_xyyy_0[i] * pb_z + g_0_xxxy_0_xyyy_1[i] * wp_z[i];

        g_0_xxxyz_0_xyyz_0[i] = 2.0 * g_0_xxxz_0_xyz_1[i] * fi_abcd_0 + g_0_xxxz_0_xyyz_0[i] * pb_y + g_0_xxxz_0_xyyz_1[i] * wp_y[i];

        g_0_xxxyz_0_xyzz_0[i] = g_0_xxxz_0_xzz_1[i] * fi_abcd_0 + g_0_xxxz_0_xyzz_0[i] * pb_y + g_0_xxxz_0_xyzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xzzz_0[i] = g_0_xxxz_0_xzzz_0[i] * pb_y + g_0_xxxz_0_xzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_yyyy_0[i] = g_0_xxxy_0_yyyy_0[i] * pb_z + g_0_xxxy_0_yyyy_1[i] * wp_z[i];

        g_0_xxxyz_0_yyyz_0[i] = 3.0 * g_0_xxxz_0_yyz_1[i] * fi_abcd_0 + g_0_xxxz_0_yyyz_0[i] * pb_y + g_0_xxxz_0_yyyz_1[i] * wp_y[i];

        g_0_xxxyz_0_yyzz_0[i] = 2.0 * g_0_xxxz_0_yzz_1[i] * fi_abcd_0 + g_0_xxxz_0_yyzz_0[i] * pb_y + g_0_xxxz_0_yyzz_1[i] * wp_y[i];

        g_0_xxxyz_0_yzzz_0[i] = g_0_xxxz_0_zzz_1[i] * fi_abcd_0 + g_0_xxxz_0_yzzz_0[i] * pb_y + g_0_xxxz_0_yzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_zzzz_0[i] = g_0_xxxz_0_zzzz_0[i] * pb_y + g_0_xxxz_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 75-90 components of targeted buffer : SHSG

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

#pragma omp simd aligned(g_0_xxx_0_xxxx_0,       \
                             g_0_xxx_0_xxxx_1,   \
                             g_0_xxx_0_xxxy_0,   \
                             g_0_xxx_0_xxxy_1,   \
                             g_0_xxx_0_xxyy_0,   \
                             g_0_xxx_0_xxyy_1,   \
                             g_0_xxx_0_xyyy_0,   \
                             g_0_xxx_0_xyyy_1,   \
                             g_0_xxxz_0_xxxx_0,  \
                             g_0_xxxz_0_xxxx_1,  \
                             g_0_xxxz_0_xxxy_0,  \
                             g_0_xxxz_0_xxxy_1,  \
                             g_0_xxxz_0_xxyy_0,  \
                             g_0_xxxz_0_xxyy_1,  \
                             g_0_xxxz_0_xyyy_0,  \
                             g_0_xxxz_0_xyyy_1,  \
                             g_0_xxxzz_0_xxxx_0, \
                             g_0_xxxzz_0_xxxy_0, \
                             g_0_xxxzz_0_xxxz_0, \
                             g_0_xxxzz_0_xxyy_0, \
                             g_0_xxxzz_0_xxyz_0, \
                             g_0_xxxzz_0_xxzz_0, \
                             g_0_xxxzz_0_xyyy_0, \
                             g_0_xxxzz_0_xyyz_0, \
                             g_0_xxxzz_0_xyzz_0, \
                             g_0_xxxzz_0_xzzz_0, \
                             g_0_xxxzz_0_yyyy_0, \
                             g_0_xxxzz_0_yyyz_0, \
                             g_0_xxxzz_0_yyzz_0, \
                             g_0_xxxzz_0_yzzz_0, \
                             g_0_xxxzz_0_zzzz_0, \
                             g_0_xxzz_0_xxxz_0,  \
                             g_0_xxzz_0_xxxz_1,  \
                             g_0_xxzz_0_xxyz_0,  \
                             g_0_xxzz_0_xxyz_1,  \
                             g_0_xxzz_0_xxz_1,   \
                             g_0_xxzz_0_xxzz_0,  \
                             g_0_xxzz_0_xxzz_1,  \
                             g_0_xxzz_0_xyyz_0,  \
                             g_0_xxzz_0_xyyz_1,  \
                             g_0_xxzz_0_xyz_1,   \
                             g_0_xxzz_0_xyzz_0,  \
                             g_0_xxzz_0_xyzz_1,  \
                             g_0_xxzz_0_xzz_1,   \
                             g_0_xxzz_0_xzzz_0,  \
                             g_0_xxzz_0_xzzz_1,  \
                             g_0_xxzz_0_yyyy_0,  \
                             g_0_xxzz_0_yyyy_1,  \
                             g_0_xxzz_0_yyyz_0,  \
                             g_0_xxzz_0_yyyz_1,  \
                             g_0_xxzz_0_yyz_1,   \
                             g_0_xxzz_0_yyzz_0,  \
                             g_0_xxzz_0_yyzz_1,  \
                             g_0_xxzz_0_yzz_1,   \
                             g_0_xxzz_0_yzzz_0,  \
                             g_0_xxzz_0_yzzz_1,  \
                             g_0_xxzz_0_zzz_1,   \
                             g_0_xxzz_0_zzzz_0,  \
                             g_0_xxzz_0_zzzz_1,  \
                             g_0_xzz_0_xxxz_0,   \
                             g_0_xzz_0_xxxz_1,   \
                             g_0_xzz_0_xxyz_0,   \
                             g_0_xzz_0_xxyz_1,   \
                             g_0_xzz_0_xxzz_0,   \
                             g_0_xzz_0_xxzz_1,   \
                             g_0_xzz_0_xyyz_0,   \
                             g_0_xzz_0_xyyz_1,   \
                             g_0_xzz_0_xyzz_0,   \
                             g_0_xzz_0_xyzz_1,   \
                             g_0_xzz_0_xzzz_0,   \
                             g_0_xzz_0_xzzz_1,   \
                             g_0_xzz_0_yyyy_0,   \
                             g_0_xzz_0_yyyy_1,   \
                             g_0_xzz_0_yyyz_0,   \
                             g_0_xzz_0_yyyz_1,   \
                             g_0_xzz_0_yyzz_0,   \
                             g_0_xzz_0_yyzz_1,   \
                             g_0_xzz_0_yzzz_0,   \
                             g_0_xzz_0_yzzz_1,   \
                             g_0_xzz_0_zzzz_0,   \
                             g_0_xzz_0_zzzz_1,   \
                             wp_x,               \
                             wp_z,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxzz_0_xxxx_0[i] =
            g_0_xxx_0_xxxx_0[i] * fi_ab_0 - g_0_xxx_0_xxxx_1[i] * fti_ab_0 + g_0_xxxz_0_xxxx_0[i] * pb_z + g_0_xxxz_0_xxxx_1[i] * wp_z[i];

        g_0_xxxzz_0_xxxy_0[i] =
            g_0_xxx_0_xxxy_0[i] * fi_ab_0 - g_0_xxx_0_xxxy_1[i] * fti_ab_0 + g_0_xxxz_0_xxxy_0[i] * pb_z + g_0_xxxz_0_xxxy_1[i] * wp_z[i];

        g_0_xxxzz_0_xxxz_0[i] = 2.0 * g_0_xzz_0_xxxz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxxz_1[i] * fti_ab_0 + 3.0 * g_0_xxzz_0_xxz_1[i] * fi_abcd_0 +
                                g_0_xxzz_0_xxxz_0[i] * pb_x + g_0_xxzz_0_xxxz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxyy_0[i] =
            g_0_xxx_0_xxyy_0[i] * fi_ab_0 - g_0_xxx_0_xxyy_1[i] * fti_ab_0 + g_0_xxxz_0_xxyy_0[i] * pb_z + g_0_xxxz_0_xxyy_1[i] * wp_z[i];

        g_0_xxxzz_0_xxyz_0[i] = 2.0 * g_0_xzz_0_xxyz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxyz_1[i] * fti_ab_0 + 2.0 * g_0_xxzz_0_xyz_1[i] * fi_abcd_0 +
                                g_0_xxzz_0_xxyz_0[i] * pb_x + g_0_xxzz_0_xxyz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxzz_0[i] = 2.0 * g_0_xzz_0_xxzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxzz_1[i] * fti_ab_0 + 2.0 * g_0_xxzz_0_xzz_1[i] * fi_abcd_0 +
                                g_0_xxzz_0_xxzz_0[i] * pb_x + g_0_xxzz_0_xxzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xyyy_0[i] =
            g_0_xxx_0_xyyy_0[i] * fi_ab_0 - g_0_xxx_0_xyyy_1[i] * fti_ab_0 + g_0_xxxz_0_xyyy_0[i] * pb_z + g_0_xxxz_0_xyyy_1[i] * wp_z[i];

        g_0_xxxzz_0_xyyz_0[i] = 2.0 * g_0_xzz_0_xyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xyyz_1[i] * fti_ab_0 + g_0_xxzz_0_yyz_1[i] * fi_abcd_0 +
                                g_0_xxzz_0_xyyz_0[i] * pb_x + g_0_xxzz_0_xyyz_1[i] * wp_x[i];

        g_0_xxxzz_0_xyzz_0[i] = 2.0 * g_0_xzz_0_xyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xyzz_1[i] * fti_ab_0 + g_0_xxzz_0_yzz_1[i] * fi_abcd_0 +
                                g_0_xxzz_0_xyzz_0[i] * pb_x + g_0_xxzz_0_xyzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xzzz_0[i] = 2.0 * g_0_xzz_0_xzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xzzz_1[i] * fti_ab_0 + g_0_xxzz_0_zzz_1[i] * fi_abcd_0 +
                                g_0_xxzz_0_xzzz_0[i] * pb_x + g_0_xxzz_0_xzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_yyyy_0[i] =
            2.0 * g_0_xzz_0_yyyy_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_yyyy_1[i] * fti_ab_0 + g_0_xxzz_0_yyyy_0[i] * pb_x + g_0_xxzz_0_yyyy_1[i] * wp_x[i];

        g_0_xxxzz_0_yyyz_0[i] =
            2.0 * g_0_xzz_0_yyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_yyyz_1[i] * fti_ab_0 + g_0_xxzz_0_yyyz_0[i] * pb_x + g_0_xxzz_0_yyyz_1[i] * wp_x[i];

        g_0_xxxzz_0_yyzz_0[i] =
            2.0 * g_0_xzz_0_yyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_yyzz_1[i] * fti_ab_0 + g_0_xxzz_0_yyzz_0[i] * pb_x + g_0_xxzz_0_yyzz_1[i] * wp_x[i];

        g_0_xxxzz_0_yzzz_0[i] =
            2.0 * g_0_xzz_0_yzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_yzzz_1[i] * fti_ab_0 + g_0_xxzz_0_yzzz_0[i] * pb_x + g_0_xxzz_0_yzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_zzzz_0[i] =
            2.0 * g_0_xzz_0_zzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_zzzz_1[i] * fti_ab_0 + g_0_xxzz_0_zzzz_0[i] * pb_x + g_0_xxzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 90-105 components of targeted buffer : SHSG

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

#pragma omp simd aligned(g_0_xxy_0_xxxx_0,       \
                             g_0_xxy_0_xxxx_1,   \
                             g_0_xxy_0_xxxz_0,   \
                             g_0_xxy_0_xxxz_1,   \
                             g_0_xxy_0_xxzz_0,   \
                             g_0_xxy_0_xxzz_1,   \
                             g_0_xxy_0_xzzz_0,   \
                             g_0_xxy_0_xzzz_1,   \
                             g_0_xxyy_0_xxxx_0,  \
                             g_0_xxyy_0_xxxx_1,  \
                             g_0_xxyy_0_xxxz_0,  \
                             g_0_xxyy_0_xxxz_1,  \
                             g_0_xxyy_0_xxzz_0,  \
                             g_0_xxyy_0_xxzz_1,  \
                             g_0_xxyy_0_xzzz_0,  \
                             g_0_xxyy_0_xzzz_1,  \
                             g_0_xxyyy_0_xxxx_0, \
                             g_0_xxyyy_0_xxxy_0, \
                             g_0_xxyyy_0_xxxz_0, \
                             g_0_xxyyy_0_xxyy_0, \
                             g_0_xxyyy_0_xxyz_0, \
                             g_0_xxyyy_0_xxzz_0, \
                             g_0_xxyyy_0_xyyy_0, \
                             g_0_xxyyy_0_xyyz_0, \
                             g_0_xxyyy_0_xyzz_0, \
                             g_0_xxyyy_0_xzzz_0, \
                             g_0_xxyyy_0_yyyy_0, \
                             g_0_xxyyy_0_yyyz_0, \
                             g_0_xxyyy_0_yyzz_0, \
                             g_0_xxyyy_0_yzzz_0, \
                             g_0_xxyyy_0_zzzz_0, \
                             g_0_xyyy_0_xxxy_0,  \
                             g_0_xyyy_0_xxxy_1,  \
                             g_0_xyyy_0_xxy_1,   \
                             g_0_xyyy_0_xxyy_0,  \
                             g_0_xyyy_0_xxyy_1,  \
                             g_0_xyyy_0_xxyz_0,  \
                             g_0_xyyy_0_xxyz_1,  \
                             g_0_xyyy_0_xyy_1,   \
                             g_0_xyyy_0_xyyy_0,  \
                             g_0_xyyy_0_xyyy_1,  \
                             g_0_xyyy_0_xyyz_0,  \
                             g_0_xyyy_0_xyyz_1,  \
                             g_0_xyyy_0_xyz_1,   \
                             g_0_xyyy_0_xyzz_0,  \
                             g_0_xyyy_0_xyzz_1,  \
                             g_0_xyyy_0_yyy_1,   \
                             g_0_xyyy_0_yyyy_0,  \
                             g_0_xyyy_0_yyyy_1,  \
                             g_0_xyyy_0_yyyz_0,  \
                             g_0_xyyy_0_yyyz_1,  \
                             g_0_xyyy_0_yyz_1,   \
                             g_0_xyyy_0_yyzz_0,  \
                             g_0_xyyy_0_yyzz_1,  \
                             g_0_xyyy_0_yzz_1,   \
                             g_0_xyyy_0_yzzz_0,  \
                             g_0_xyyy_0_yzzz_1,  \
                             g_0_xyyy_0_zzzz_0,  \
                             g_0_xyyy_0_zzzz_1,  \
                             g_0_yyy_0_xxxy_0,   \
                             g_0_yyy_0_xxxy_1,   \
                             g_0_yyy_0_xxyy_0,   \
                             g_0_yyy_0_xxyy_1,   \
                             g_0_yyy_0_xxyz_0,   \
                             g_0_yyy_0_xxyz_1,   \
                             g_0_yyy_0_xyyy_0,   \
                             g_0_yyy_0_xyyy_1,   \
                             g_0_yyy_0_xyyz_0,   \
                             g_0_yyy_0_xyyz_1,   \
                             g_0_yyy_0_xyzz_0,   \
                             g_0_yyy_0_xyzz_1,   \
                             g_0_yyy_0_yyyy_0,   \
                             g_0_yyy_0_yyyy_1,   \
                             g_0_yyy_0_yyyz_0,   \
                             g_0_yyy_0_yyyz_1,   \
                             g_0_yyy_0_yyzz_0,   \
                             g_0_yyy_0_yyzz_1,   \
                             g_0_yyy_0_yzzz_0,   \
                             g_0_yyy_0_yzzz_1,   \
                             g_0_yyy_0_zzzz_0,   \
                             g_0_yyy_0_zzzz_1,   \
                             wp_x,               \
                             wp_y,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyy_0_xxxx_0[i] =
            2.0 * g_0_xxy_0_xxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxy_0_xxxx_1[i] * fti_ab_0 + g_0_xxyy_0_xxxx_0[i] * pb_y + g_0_xxyy_0_xxxx_1[i] * wp_y[i];

        g_0_xxyyy_0_xxxy_0[i] = g_0_yyy_0_xxxy_0[i] * fi_ab_0 - g_0_yyy_0_xxxy_1[i] * fti_ab_0 + 3.0 * g_0_xyyy_0_xxy_1[i] * fi_abcd_0 +
                                g_0_xyyy_0_xxxy_0[i] * pb_x + g_0_xyyy_0_xxxy_1[i] * wp_x[i];

        g_0_xxyyy_0_xxxz_0[i] =
            2.0 * g_0_xxy_0_xxxz_0[i] * fi_ab_0 - 2.0 * g_0_xxy_0_xxxz_1[i] * fti_ab_0 + g_0_xxyy_0_xxxz_0[i] * pb_y + g_0_xxyy_0_xxxz_1[i] * wp_y[i];

        g_0_xxyyy_0_xxyy_0[i] = g_0_yyy_0_xxyy_0[i] * fi_ab_0 - g_0_yyy_0_xxyy_1[i] * fti_ab_0 + 2.0 * g_0_xyyy_0_xyy_1[i] * fi_abcd_0 +
                                g_0_xyyy_0_xxyy_0[i] * pb_x + g_0_xyyy_0_xxyy_1[i] * wp_x[i];

        g_0_xxyyy_0_xxyz_0[i] = g_0_yyy_0_xxyz_0[i] * fi_ab_0 - g_0_yyy_0_xxyz_1[i] * fti_ab_0 + 2.0 * g_0_xyyy_0_xyz_1[i] * fi_abcd_0 +
                                g_0_xyyy_0_xxyz_0[i] * pb_x + g_0_xyyy_0_xxyz_1[i] * wp_x[i];

        g_0_xxyyy_0_xxzz_0[i] =
            2.0 * g_0_xxy_0_xxzz_0[i] * fi_ab_0 - 2.0 * g_0_xxy_0_xxzz_1[i] * fti_ab_0 + g_0_xxyy_0_xxzz_0[i] * pb_y + g_0_xxyy_0_xxzz_1[i] * wp_y[i];

        g_0_xxyyy_0_xyyy_0[i] = g_0_yyy_0_xyyy_0[i] * fi_ab_0 - g_0_yyy_0_xyyy_1[i] * fti_ab_0 + g_0_xyyy_0_yyy_1[i] * fi_abcd_0 +
                                g_0_xyyy_0_xyyy_0[i] * pb_x + g_0_xyyy_0_xyyy_1[i] * wp_x[i];

        g_0_xxyyy_0_xyyz_0[i] = g_0_yyy_0_xyyz_0[i] * fi_ab_0 - g_0_yyy_0_xyyz_1[i] * fti_ab_0 + g_0_xyyy_0_yyz_1[i] * fi_abcd_0 +
                                g_0_xyyy_0_xyyz_0[i] * pb_x + g_0_xyyy_0_xyyz_1[i] * wp_x[i];

        g_0_xxyyy_0_xyzz_0[i] = g_0_yyy_0_xyzz_0[i] * fi_ab_0 - g_0_yyy_0_xyzz_1[i] * fti_ab_0 + g_0_xyyy_0_yzz_1[i] * fi_abcd_0 +
                                g_0_xyyy_0_xyzz_0[i] * pb_x + g_0_xyyy_0_xyzz_1[i] * wp_x[i];

        g_0_xxyyy_0_xzzz_0[i] =
            2.0 * g_0_xxy_0_xzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxy_0_xzzz_1[i] * fti_ab_0 + g_0_xxyy_0_xzzz_0[i] * pb_y + g_0_xxyy_0_xzzz_1[i] * wp_y[i];

        g_0_xxyyy_0_yyyy_0[i] =
            g_0_yyy_0_yyyy_0[i] * fi_ab_0 - g_0_yyy_0_yyyy_1[i] * fti_ab_0 + g_0_xyyy_0_yyyy_0[i] * pb_x + g_0_xyyy_0_yyyy_1[i] * wp_x[i];

        g_0_xxyyy_0_yyyz_0[i] =
            g_0_yyy_0_yyyz_0[i] * fi_ab_0 - g_0_yyy_0_yyyz_1[i] * fti_ab_0 + g_0_xyyy_0_yyyz_0[i] * pb_x + g_0_xyyy_0_yyyz_1[i] * wp_x[i];

        g_0_xxyyy_0_yyzz_0[i] =
            g_0_yyy_0_yyzz_0[i] * fi_ab_0 - g_0_yyy_0_yyzz_1[i] * fti_ab_0 + g_0_xyyy_0_yyzz_0[i] * pb_x + g_0_xyyy_0_yyzz_1[i] * wp_x[i];

        g_0_xxyyy_0_yzzz_0[i] =
            g_0_yyy_0_yzzz_0[i] * fi_ab_0 - g_0_yyy_0_yzzz_1[i] * fti_ab_0 + g_0_xyyy_0_yzzz_0[i] * pb_x + g_0_xyyy_0_yzzz_1[i] * wp_x[i];

        g_0_xxyyy_0_zzzz_0[i] =
            g_0_yyy_0_zzzz_0[i] * fi_ab_0 - g_0_yyy_0_zzzz_1[i] * fti_ab_0 + g_0_xyyy_0_zzzz_0[i] * pb_x + g_0_xyyy_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 105-120 components of targeted buffer : SHSG

    auto g_0_xxyyz_0_xxxx_0 = pbuffer.data(idx_eri_0_shsg + 105);

    auto g_0_xxyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_shsg + 106);

    auto g_0_xxyyz_0_xxxz_0 = pbuffer.data(idx_eri_0_shsg + 107);

    auto g_0_xxyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_shsg + 108);

    auto g_0_xxyyz_0_xxyz_0 = pbuffer.data(idx_eri_0_shsg + 109);

    auto g_0_xxyyz_0_xxzz_0 = pbuffer.data(idx_eri_0_shsg + 110);

    auto g_0_xxyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_shsg + 111);

    auto g_0_xxyyz_0_xyyz_0 = pbuffer.data(idx_eri_0_shsg + 112);

    auto g_0_xxyyz_0_xyzz_0 = pbuffer.data(idx_eri_0_shsg + 113);

    auto g_0_xxyyz_0_xzzz_0 = pbuffer.data(idx_eri_0_shsg + 114);

    auto g_0_xxyyz_0_yyyy_0 = pbuffer.data(idx_eri_0_shsg + 115);

    auto g_0_xxyyz_0_yyyz_0 = pbuffer.data(idx_eri_0_shsg + 116);

    auto g_0_xxyyz_0_yyzz_0 = pbuffer.data(idx_eri_0_shsg + 117);

    auto g_0_xxyyz_0_yzzz_0 = pbuffer.data(idx_eri_0_shsg + 118);

    auto g_0_xxyyz_0_zzzz_0 = pbuffer.data(idx_eri_0_shsg + 119);

#pragma omp simd aligned(g_0_xxyy_0_xxx_1,       \
                             g_0_xxyy_0_xxxx_0,  \
                             g_0_xxyy_0_xxxx_1,  \
                             g_0_xxyy_0_xxxy_0,  \
                             g_0_xxyy_0_xxxy_1,  \
                             g_0_xxyy_0_xxxz_0,  \
                             g_0_xxyy_0_xxxz_1,  \
                             g_0_xxyy_0_xxy_1,   \
                             g_0_xxyy_0_xxyy_0,  \
                             g_0_xxyy_0_xxyy_1,  \
                             g_0_xxyy_0_xxyz_0,  \
                             g_0_xxyy_0_xxyz_1,  \
                             g_0_xxyy_0_xxz_1,   \
                             g_0_xxyy_0_xxzz_0,  \
                             g_0_xxyy_0_xxzz_1,  \
                             g_0_xxyy_0_xyy_1,   \
                             g_0_xxyy_0_xyyy_0,  \
                             g_0_xxyy_0_xyyy_1,  \
                             g_0_xxyy_0_xyyz_0,  \
                             g_0_xxyy_0_xyyz_1,  \
                             g_0_xxyy_0_xyz_1,   \
                             g_0_xxyy_0_xyzz_0,  \
                             g_0_xxyy_0_xyzz_1,  \
                             g_0_xxyy_0_xzz_1,   \
                             g_0_xxyy_0_xzzz_0,  \
                             g_0_xxyy_0_xzzz_1,  \
                             g_0_xxyy_0_yyy_1,   \
                             g_0_xxyy_0_yyyy_0,  \
                             g_0_xxyy_0_yyyy_1,  \
                             g_0_xxyy_0_yyyz_0,  \
                             g_0_xxyy_0_yyyz_1,  \
                             g_0_xxyy_0_yyz_1,   \
                             g_0_xxyy_0_yyzz_0,  \
                             g_0_xxyy_0_yyzz_1,  \
                             g_0_xxyy_0_yzz_1,   \
                             g_0_xxyy_0_yzzz_0,  \
                             g_0_xxyy_0_yzzz_1,  \
                             g_0_xxyy_0_zzz_1,   \
                             g_0_xxyy_0_zzzz_0,  \
                             g_0_xxyy_0_zzzz_1,  \
                             g_0_xxyyz_0_xxxx_0, \
                             g_0_xxyyz_0_xxxy_0, \
                             g_0_xxyyz_0_xxxz_0, \
                             g_0_xxyyz_0_xxyy_0, \
                             g_0_xxyyz_0_xxyz_0, \
                             g_0_xxyyz_0_xxzz_0, \
                             g_0_xxyyz_0_xyyy_0, \
                             g_0_xxyyz_0_xyyz_0, \
                             g_0_xxyyz_0_xyzz_0, \
                             g_0_xxyyz_0_xzzz_0, \
                             g_0_xxyyz_0_yyyy_0, \
                             g_0_xxyyz_0_yyyz_0, \
                             g_0_xxyyz_0_yyzz_0, \
                             g_0_xxyyz_0_yzzz_0, \
                             g_0_xxyyz_0_zzzz_0, \
                             wp_z,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyz_0_xxxx_0[i] = g_0_xxyy_0_xxxx_0[i] * pb_z + g_0_xxyy_0_xxxx_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxy_0[i] = g_0_xxyy_0_xxxy_0[i] * pb_z + g_0_xxyy_0_xxxy_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxz_0[i] = g_0_xxyy_0_xxx_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxz_0[i] * pb_z + g_0_xxyy_0_xxxz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxyy_0[i] = g_0_xxyy_0_xxyy_0[i] * pb_z + g_0_xxyy_0_xxyy_1[i] * wp_z[i];

        g_0_xxyyz_0_xxyz_0[i] = g_0_xxyy_0_xxy_1[i] * fi_abcd_0 + g_0_xxyy_0_xxyz_0[i] * pb_z + g_0_xxyy_0_xxyz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxzz_0[i] = 2.0 * g_0_xxyy_0_xxz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxzz_0[i] * pb_z + g_0_xxyy_0_xxzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xyyy_0[i] = g_0_xxyy_0_xyyy_0[i] * pb_z + g_0_xxyy_0_xyyy_1[i] * wp_z[i];

        g_0_xxyyz_0_xyyz_0[i] = g_0_xxyy_0_xyy_1[i] * fi_abcd_0 + g_0_xxyy_0_xyyz_0[i] * pb_z + g_0_xxyy_0_xyyz_1[i] * wp_z[i];

        g_0_xxyyz_0_xyzz_0[i] = 2.0 * g_0_xxyy_0_xyz_1[i] * fi_abcd_0 + g_0_xxyy_0_xyzz_0[i] * pb_z + g_0_xxyy_0_xyzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xzzz_0[i] = 3.0 * g_0_xxyy_0_xzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xzzz_0[i] * pb_z + g_0_xxyy_0_xzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_yyyy_0[i] = g_0_xxyy_0_yyyy_0[i] * pb_z + g_0_xxyy_0_yyyy_1[i] * wp_z[i];

        g_0_xxyyz_0_yyyz_0[i] = g_0_xxyy_0_yyy_1[i] * fi_abcd_0 + g_0_xxyy_0_yyyz_0[i] * pb_z + g_0_xxyy_0_yyyz_1[i] * wp_z[i];

        g_0_xxyyz_0_yyzz_0[i] = 2.0 * g_0_xxyy_0_yyz_1[i] * fi_abcd_0 + g_0_xxyy_0_yyzz_0[i] * pb_z + g_0_xxyy_0_yyzz_1[i] * wp_z[i];

        g_0_xxyyz_0_yzzz_0[i] = 3.0 * g_0_xxyy_0_yzz_1[i] * fi_abcd_0 + g_0_xxyy_0_yzzz_0[i] * pb_z + g_0_xxyy_0_yzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_zzzz_0[i] = 4.0 * g_0_xxyy_0_zzz_1[i] * fi_abcd_0 + g_0_xxyy_0_zzzz_0[i] * pb_z + g_0_xxyy_0_zzzz_1[i] * wp_z[i];
    }

    /// Set up 120-135 components of targeted buffer : SHSG

    auto g_0_xxyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_shsg + 120);

    auto g_0_xxyzz_0_xxxy_0 = pbuffer.data(idx_eri_0_shsg + 121);

    auto g_0_xxyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_shsg + 122);

    auto g_0_xxyzz_0_xxyy_0 = pbuffer.data(idx_eri_0_shsg + 123);

    auto g_0_xxyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_shsg + 124);

    auto g_0_xxyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_shsg + 125);

    auto g_0_xxyzz_0_xyyy_0 = pbuffer.data(idx_eri_0_shsg + 126);

    auto g_0_xxyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_shsg + 127);

    auto g_0_xxyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_shsg + 128);

    auto g_0_xxyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_shsg + 129);

    auto g_0_xxyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_shsg + 130);

    auto g_0_xxyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_shsg + 131);

    auto g_0_xxyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_shsg + 132);

    auto g_0_xxyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_shsg + 133);

    auto g_0_xxyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_shsg + 134);

#pragma omp simd aligned(g_0_xxyzz_0_xxxx_0,     \
                             g_0_xxyzz_0_xxxy_0, \
                             g_0_xxyzz_0_xxxz_0, \
                             g_0_xxyzz_0_xxyy_0, \
                             g_0_xxyzz_0_xxyz_0, \
                             g_0_xxyzz_0_xxzz_0, \
                             g_0_xxyzz_0_xyyy_0, \
                             g_0_xxyzz_0_xyyz_0, \
                             g_0_xxyzz_0_xyzz_0, \
                             g_0_xxyzz_0_xzzz_0, \
                             g_0_xxyzz_0_yyyy_0, \
                             g_0_xxyzz_0_yyyz_0, \
                             g_0_xxyzz_0_yyzz_0, \
                             g_0_xxyzz_0_yzzz_0, \
                             g_0_xxyzz_0_zzzz_0, \
                             g_0_xxzz_0_xxx_1,   \
                             g_0_xxzz_0_xxxx_0,  \
                             g_0_xxzz_0_xxxx_1,  \
                             g_0_xxzz_0_xxxy_0,  \
                             g_0_xxzz_0_xxxy_1,  \
                             g_0_xxzz_0_xxxz_0,  \
                             g_0_xxzz_0_xxxz_1,  \
                             g_0_xxzz_0_xxy_1,   \
                             g_0_xxzz_0_xxyy_0,  \
                             g_0_xxzz_0_xxyy_1,  \
                             g_0_xxzz_0_xxyz_0,  \
                             g_0_xxzz_0_xxyz_1,  \
                             g_0_xxzz_0_xxz_1,   \
                             g_0_xxzz_0_xxzz_0,  \
                             g_0_xxzz_0_xxzz_1,  \
                             g_0_xxzz_0_xyy_1,   \
                             g_0_xxzz_0_xyyy_0,  \
                             g_0_xxzz_0_xyyy_1,  \
                             g_0_xxzz_0_xyyz_0,  \
                             g_0_xxzz_0_xyyz_1,  \
                             g_0_xxzz_0_xyz_1,   \
                             g_0_xxzz_0_xyzz_0,  \
                             g_0_xxzz_0_xyzz_1,  \
                             g_0_xxzz_0_xzz_1,   \
                             g_0_xxzz_0_xzzz_0,  \
                             g_0_xxzz_0_xzzz_1,  \
                             g_0_xxzz_0_yyy_1,   \
                             g_0_xxzz_0_yyyy_0,  \
                             g_0_xxzz_0_yyyy_1,  \
                             g_0_xxzz_0_yyyz_0,  \
                             g_0_xxzz_0_yyyz_1,  \
                             g_0_xxzz_0_yyz_1,   \
                             g_0_xxzz_0_yyzz_0,  \
                             g_0_xxzz_0_yyzz_1,  \
                             g_0_xxzz_0_yzz_1,   \
                             g_0_xxzz_0_yzzz_0,  \
                             g_0_xxzz_0_yzzz_1,  \
                             g_0_xxzz_0_zzz_1,   \
                             g_0_xxzz_0_zzzz_0,  \
                             g_0_xxzz_0_zzzz_1,  \
                             wp_y,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyzz_0_xxxx_0[i] = g_0_xxzz_0_xxxx_0[i] * pb_y + g_0_xxzz_0_xxxx_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxy_0[i] = g_0_xxzz_0_xxx_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxy_0[i] * pb_y + g_0_xxzz_0_xxxy_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxz_0[i] = g_0_xxzz_0_xxxz_0[i] * pb_y + g_0_xxzz_0_xxxz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxyy_0[i] = 2.0 * g_0_xxzz_0_xxy_1[i] * fi_abcd_0 + g_0_xxzz_0_xxyy_0[i] * pb_y + g_0_xxzz_0_xxyy_1[i] * wp_y[i];

        g_0_xxyzz_0_xxyz_0[i] = g_0_xxzz_0_xxz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxyz_0[i] * pb_y + g_0_xxzz_0_xxyz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxzz_0[i] = g_0_xxzz_0_xxzz_0[i] * pb_y + g_0_xxzz_0_xxzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xyyy_0[i] = 3.0 * g_0_xxzz_0_xyy_1[i] * fi_abcd_0 + g_0_xxzz_0_xyyy_0[i] * pb_y + g_0_xxzz_0_xyyy_1[i] * wp_y[i];

        g_0_xxyzz_0_xyyz_0[i] = 2.0 * g_0_xxzz_0_xyz_1[i] * fi_abcd_0 + g_0_xxzz_0_xyyz_0[i] * pb_y + g_0_xxzz_0_xyyz_1[i] * wp_y[i];

        g_0_xxyzz_0_xyzz_0[i] = g_0_xxzz_0_xzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xyzz_0[i] * pb_y + g_0_xxzz_0_xyzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xzzz_0[i] = g_0_xxzz_0_xzzz_0[i] * pb_y + g_0_xxzz_0_xzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_yyyy_0[i] = 4.0 * g_0_xxzz_0_yyy_1[i] * fi_abcd_0 + g_0_xxzz_0_yyyy_0[i] * pb_y + g_0_xxzz_0_yyyy_1[i] * wp_y[i];

        g_0_xxyzz_0_yyyz_0[i] = 3.0 * g_0_xxzz_0_yyz_1[i] * fi_abcd_0 + g_0_xxzz_0_yyyz_0[i] * pb_y + g_0_xxzz_0_yyyz_1[i] * wp_y[i];

        g_0_xxyzz_0_yyzz_0[i] = 2.0 * g_0_xxzz_0_yzz_1[i] * fi_abcd_0 + g_0_xxzz_0_yyzz_0[i] * pb_y + g_0_xxzz_0_yyzz_1[i] * wp_y[i];

        g_0_xxyzz_0_yzzz_0[i] = g_0_xxzz_0_zzz_1[i] * fi_abcd_0 + g_0_xxzz_0_yzzz_0[i] * pb_y + g_0_xxzz_0_yzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_zzzz_0[i] = g_0_xxzz_0_zzzz_0[i] * pb_y + g_0_xxzz_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 135-150 components of targeted buffer : SHSG

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

#pragma omp simd aligned(g_0_xxz_0_xxxx_0,       \
                             g_0_xxz_0_xxxx_1,   \
                             g_0_xxz_0_xxxy_0,   \
                             g_0_xxz_0_xxxy_1,   \
                             g_0_xxz_0_xxyy_0,   \
                             g_0_xxz_0_xxyy_1,   \
                             g_0_xxz_0_xyyy_0,   \
                             g_0_xxz_0_xyyy_1,   \
                             g_0_xxzz_0_xxxx_0,  \
                             g_0_xxzz_0_xxxx_1,  \
                             g_0_xxzz_0_xxxy_0,  \
                             g_0_xxzz_0_xxxy_1,  \
                             g_0_xxzz_0_xxyy_0,  \
                             g_0_xxzz_0_xxyy_1,  \
                             g_0_xxzz_0_xyyy_0,  \
                             g_0_xxzz_0_xyyy_1,  \
                             g_0_xxzzz_0_xxxx_0, \
                             g_0_xxzzz_0_xxxy_0, \
                             g_0_xxzzz_0_xxxz_0, \
                             g_0_xxzzz_0_xxyy_0, \
                             g_0_xxzzz_0_xxyz_0, \
                             g_0_xxzzz_0_xxzz_0, \
                             g_0_xxzzz_0_xyyy_0, \
                             g_0_xxzzz_0_xyyz_0, \
                             g_0_xxzzz_0_xyzz_0, \
                             g_0_xxzzz_0_xzzz_0, \
                             g_0_xxzzz_0_yyyy_0, \
                             g_0_xxzzz_0_yyyz_0, \
                             g_0_xxzzz_0_yyzz_0, \
                             g_0_xxzzz_0_yzzz_0, \
                             g_0_xxzzz_0_zzzz_0, \
                             g_0_xzzz_0_xxxz_0,  \
                             g_0_xzzz_0_xxxz_1,  \
                             g_0_xzzz_0_xxyz_0,  \
                             g_0_xzzz_0_xxyz_1,  \
                             g_0_xzzz_0_xxz_1,   \
                             g_0_xzzz_0_xxzz_0,  \
                             g_0_xzzz_0_xxzz_1,  \
                             g_0_xzzz_0_xyyz_0,  \
                             g_0_xzzz_0_xyyz_1,  \
                             g_0_xzzz_0_xyz_1,   \
                             g_0_xzzz_0_xyzz_0,  \
                             g_0_xzzz_0_xyzz_1,  \
                             g_0_xzzz_0_xzz_1,   \
                             g_0_xzzz_0_xzzz_0,  \
                             g_0_xzzz_0_xzzz_1,  \
                             g_0_xzzz_0_yyyy_0,  \
                             g_0_xzzz_0_yyyy_1,  \
                             g_0_xzzz_0_yyyz_0,  \
                             g_0_xzzz_0_yyyz_1,  \
                             g_0_xzzz_0_yyz_1,   \
                             g_0_xzzz_0_yyzz_0,  \
                             g_0_xzzz_0_yyzz_1,  \
                             g_0_xzzz_0_yzz_1,   \
                             g_0_xzzz_0_yzzz_0,  \
                             g_0_xzzz_0_yzzz_1,  \
                             g_0_xzzz_0_zzz_1,   \
                             g_0_xzzz_0_zzzz_0,  \
                             g_0_xzzz_0_zzzz_1,  \
                             g_0_zzz_0_xxxz_0,   \
                             g_0_zzz_0_xxxz_1,   \
                             g_0_zzz_0_xxyz_0,   \
                             g_0_zzz_0_xxyz_1,   \
                             g_0_zzz_0_xxzz_0,   \
                             g_0_zzz_0_xxzz_1,   \
                             g_0_zzz_0_xyyz_0,   \
                             g_0_zzz_0_xyyz_1,   \
                             g_0_zzz_0_xyzz_0,   \
                             g_0_zzz_0_xyzz_1,   \
                             g_0_zzz_0_xzzz_0,   \
                             g_0_zzz_0_xzzz_1,   \
                             g_0_zzz_0_yyyy_0,   \
                             g_0_zzz_0_yyyy_1,   \
                             g_0_zzz_0_yyyz_0,   \
                             g_0_zzz_0_yyyz_1,   \
                             g_0_zzz_0_yyzz_0,   \
                             g_0_zzz_0_yyzz_1,   \
                             g_0_zzz_0_yzzz_0,   \
                             g_0_zzz_0_yzzz_1,   \
                             g_0_zzz_0_zzzz_0,   \
                             g_0_zzz_0_zzzz_1,   \
                             wp_x,               \
                             wp_z,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxzzz_0_xxxx_0[i] =
            2.0 * g_0_xxz_0_xxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxz_0_xxxx_1[i] * fti_ab_0 + g_0_xxzz_0_xxxx_0[i] * pb_z + g_0_xxzz_0_xxxx_1[i] * wp_z[i];

        g_0_xxzzz_0_xxxy_0[i] =
            2.0 * g_0_xxz_0_xxxy_0[i] * fi_ab_0 - 2.0 * g_0_xxz_0_xxxy_1[i] * fti_ab_0 + g_0_xxzz_0_xxxy_0[i] * pb_z + g_0_xxzz_0_xxxy_1[i] * wp_z[i];

        g_0_xxzzz_0_xxxz_0[i] = g_0_zzz_0_xxxz_0[i] * fi_ab_0 - g_0_zzz_0_xxxz_1[i] * fti_ab_0 + 3.0 * g_0_xzzz_0_xxz_1[i] * fi_abcd_0 +
                                g_0_xzzz_0_xxxz_0[i] * pb_x + g_0_xzzz_0_xxxz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxyy_0[i] =
            2.0 * g_0_xxz_0_xxyy_0[i] * fi_ab_0 - 2.0 * g_0_xxz_0_xxyy_1[i] * fti_ab_0 + g_0_xxzz_0_xxyy_0[i] * pb_z + g_0_xxzz_0_xxyy_1[i] * wp_z[i];

        g_0_xxzzz_0_xxyz_0[i] = g_0_zzz_0_xxyz_0[i] * fi_ab_0 - g_0_zzz_0_xxyz_1[i] * fti_ab_0 + 2.0 * g_0_xzzz_0_xyz_1[i] * fi_abcd_0 +
                                g_0_xzzz_0_xxyz_0[i] * pb_x + g_0_xzzz_0_xxyz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxzz_0[i] = g_0_zzz_0_xxzz_0[i] * fi_ab_0 - g_0_zzz_0_xxzz_1[i] * fti_ab_0 + 2.0 * g_0_xzzz_0_xzz_1[i] * fi_abcd_0 +
                                g_0_xzzz_0_xxzz_0[i] * pb_x + g_0_xzzz_0_xxzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xyyy_0[i] =
            2.0 * g_0_xxz_0_xyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxz_0_xyyy_1[i] * fti_ab_0 + g_0_xxzz_0_xyyy_0[i] * pb_z + g_0_xxzz_0_xyyy_1[i] * wp_z[i];

        g_0_xxzzz_0_xyyz_0[i] = g_0_zzz_0_xyyz_0[i] * fi_ab_0 - g_0_zzz_0_xyyz_1[i] * fti_ab_0 + g_0_xzzz_0_yyz_1[i] * fi_abcd_0 +
                                g_0_xzzz_0_xyyz_0[i] * pb_x + g_0_xzzz_0_xyyz_1[i] * wp_x[i];

        g_0_xxzzz_0_xyzz_0[i] = g_0_zzz_0_xyzz_0[i] * fi_ab_0 - g_0_zzz_0_xyzz_1[i] * fti_ab_0 + g_0_xzzz_0_yzz_1[i] * fi_abcd_0 +
                                g_0_xzzz_0_xyzz_0[i] * pb_x + g_0_xzzz_0_xyzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xzzz_0[i] = g_0_zzz_0_xzzz_0[i] * fi_ab_0 - g_0_zzz_0_xzzz_1[i] * fti_ab_0 + g_0_xzzz_0_zzz_1[i] * fi_abcd_0 +
                                g_0_xzzz_0_xzzz_0[i] * pb_x + g_0_xzzz_0_xzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_yyyy_0[i] =
            g_0_zzz_0_yyyy_0[i] * fi_ab_0 - g_0_zzz_0_yyyy_1[i] * fti_ab_0 + g_0_xzzz_0_yyyy_0[i] * pb_x + g_0_xzzz_0_yyyy_1[i] * wp_x[i];

        g_0_xxzzz_0_yyyz_0[i] =
            g_0_zzz_0_yyyz_0[i] * fi_ab_0 - g_0_zzz_0_yyyz_1[i] * fti_ab_0 + g_0_xzzz_0_yyyz_0[i] * pb_x + g_0_xzzz_0_yyyz_1[i] * wp_x[i];

        g_0_xxzzz_0_yyzz_0[i] =
            g_0_zzz_0_yyzz_0[i] * fi_ab_0 - g_0_zzz_0_yyzz_1[i] * fti_ab_0 + g_0_xzzz_0_yyzz_0[i] * pb_x + g_0_xzzz_0_yyzz_1[i] * wp_x[i];

        g_0_xxzzz_0_yzzz_0[i] =
            g_0_zzz_0_yzzz_0[i] * fi_ab_0 - g_0_zzz_0_yzzz_1[i] * fti_ab_0 + g_0_xzzz_0_yzzz_0[i] * pb_x + g_0_xzzz_0_yzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_zzzz_0[i] =
            g_0_zzz_0_zzzz_0[i] * fi_ab_0 - g_0_zzz_0_zzzz_1[i] * fti_ab_0 + g_0_xzzz_0_zzzz_0[i] * pb_x + g_0_xzzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 150-165 components of targeted buffer : SHSG

    auto g_0_xyyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_shsg + 150);

    auto g_0_xyyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_shsg + 151);

    auto g_0_xyyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_shsg + 152);

    auto g_0_xyyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_shsg + 153);

    auto g_0_xyyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_shsg + 154);

    auto g_0_xyyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_shsg + 155);

    auto g_0_xyyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_shsg + 156);

    auto g_0_xyyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_shsg + 157);

    auto g_0_xyyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_shsg + 158);

    auto g_0_xyyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_shsg + 159);

    auto g_0_xyyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_shsg + 160);

    auto g_0_xyyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_shsg + 161);

    auto g_0_xyyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_shsg + 162);

    auto g_0_xyyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_shsg + 163);

    auto g_0_xyyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_shsg + 164);

#pragma omp simd aligned(g_0_xyyyy_0_xxxx_0,     \
                             g_0_xyyyy_0_xxxy_0, \
                             g_0_xyyyy_0_xxxz_0, \
                             g_0_xyyyy_0_xxyy_0, \
                             g_0_xyyyy_0_xxyz_0, \
                             g_0_xyyyy_0_xxzz_0, \
                             g_0_xyyyy_0_xyyy_0, \
                             g_0_xyyyy_0_xyyz_0, \
                             g_0_xyyyy_0_xyzz_0, \
                             g_0_xyyyy_0_xzzz_0, \
                             g_0_xyyyy_0_yyyy_0, \
                             g_0_xyyyy_0_yyyz_0, \
                             g_0_xyyyy_0_yyzz_0, \
                             g_0_xyyyy_0_yzzz_0, \
                             g_0_xyyyy_0_zzzz_0, \
                             g_0_yyyy_0_xxx_1,   \
                             g_0_yyyy_0_xxxx_0,  \
                             g_0_yyyy_0_xxxx_1,  \
                             g_0_yyyy_0_xxxy_0,  \
                             g_0_yyyy_0_xxxy_1,  \
                             g_0_yyyy_0_xxxz_0,  \
                             g_0_yyyy_0_xxxz_1,  \
                             g_0_yyyy_0_xxy_1,   \
                             g_0_yyyy_0_xxyy_0,  \
                             g_0_yyyy_0_xxyy_1,  \
                             g_0_yyyy_0_xxyz_0,  \
                             g_0_yyyy_0_xxyz_1,  \
                             g_0_yyyy_0_xxz_1,   \
                             g_0_yyyy_0_xxzz_0,  \
                             g_0_yyyy_0_xxzz_1,  \
                             g_0_yyyy_0_xyy_1,   \
                             g_0_yyyy_0_xyyy_0,  \
                             g_0_yyyy_0_xyyy_1,  \
                             g_0_yyyy_0_xyyz_0,  \
                             g_0_yyyy_0_xyyz_1,  \
                             g_0_yyyy_0_xyz_1,   \
                             g_0_yyyy_0_xyzz_0,  \
                             g_0_yyyy_0_xyzz_1,  \
                             g_0_yyyy_0_xzz_1,   \
                             g_0_yyyy_0_xzzz_0,  \
                             g_0_yyyy_0_xzzz_1,  \
                             g_0_yyyy_0_yyy_1,   \
                             g_0_yyyy_0_yyyy_0,  \
                             g_0_yyyy_0_yyyy_1,  \
                             g_0_yyyy_0_yyyz_0,  \
                             g_0_yyyy_0_yyyz_1,  \
                             g_0_yyyy_0_yyz_1,   \
                             g_0_yyyy_0_yyzz_0,  \
                             g_0_yyyy_0_yyzz_1,  \
                             g_0_yyyy_0_yzz_1,   \
                             g_0_yyyy_0_yzzz_0,  \
                             g_0_yyyy_0_yzzz_1,  \
                             g_0_yyyy_0_zzz_1,   \
                             g_0_yyyy_0_zzzz_0,  \
                             g_0_yyyy_0_zzzz_1,  \
                             wp_x,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyy_0_xxxx_0[i] = 4.0 * g_0_yyyy_0_xxx_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxx_0[i] * pb_x + g_0_yyyy_0_xxxx_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxy_0[i] = 3.0 * g_0_yyyy_0_xxy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxy_0[i] * pb_x + g_0_yyyy_0_xxxy_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxz_0[i] = 3.0 * g_0_yyyy_0_xxz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxz_0[i] * pb_x + g_0_yyyy_0_xxxz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxyy_0[i] = 2.0 * g_0_yyyy_0_xyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyy_0[i] * pb_x + g_0_yyyy_0_xxyy_1[i] * wp_x[i];

        g_0_xyyyy_0_xxyz_0[i] = 2.0 * g_0_yyyy_0_xyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyz_0[i] * pb_x + g_0_yyyy_0_xxyz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxzz_0[i] = 2.0 * g_0_yyyy_0_xzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxzz_0[i] * pb_x + g_0_yyyy_0_xxzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xyyy_0[i] = g_0_yyyy_0_yyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyy_0[i] * pb_x + g_0_yyyy_0_xyyy_1[i] * wp_x[i];

        g_0_xyyyy_0_xyyz_0[i] = g_0_yyyy_0_yyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyz_0[i] * pb_x + g_0_yyyy_0_xyyz_1[i] * wp_x[i];

        g_0_xyyyy_0_xyzz_0[i] = g_0_yyyy_0_yzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyzz_0[i] * pb_x + g_0_yyyy_0_xyzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xzzz_0[i] = g_0_yyyy_0_zzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xzzz_0[i] * pb_x + g_0_yyyy_0_xzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_yyyy_0[i] = g_0_yyyy_0_yyyy_0[i] * pb_x + g_0_yyyy_0_yyyy_1[i] * wp_x[i];

        g_0_xyyyy_0_yyyz_0[i] = g_0_yyyy_0_yyyz_0[i] * pb_x + g_0_yyyy_0_yyyz_1[i] * wp_x[i];

        g_0_xyyyy_0_yyzz_0[i] = g_0_yyyy_0_yyzz_0[i] * pb_x + g_0_yyyy_0_yyzz_1[i] * wp_x[i];

        g_0_xyyyy_0_yzzz_0[i] = g_0_yyyy_0_yzzz_0[i] * pb_x + g_0_yyyy_0_yzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_zzzz_0[i] = g_0_yyyy_0_zzzz_0[i] * pb_x + g_0_yyyy_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 165-180 components of targeted buffer : SHSG

    auto g_0_xyyyz_0_xxxx_0 = pbuffer.data(idx_eri_0_shsg + 165);

    auto g_0_xyyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_shsg + 166);

    auto g_0_xyyyz_0_xxxz_0 = pbuffer.data(idx_eri_0_shsg + 167);

    auto g_0_xyyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_shsg + 168);

    auto g_0_xyyyz_0_xxyz_0 = pbuffer.data(idx_eri_0_shsg + 169);

    auto g_0_xyyyz_0_xxzz_0 = pbuffer.data(idx_eri_0_shsg + 170);

    auto g_0_xyyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_shsg + 171);

    auto g_0_xyyyz_0_xyyz_0 = pbuffer.data(idx_eri_0_shsg + 172);

    auto g_0_xyyyz_0_xyzz_0 = pbuffer.data(idx_eri_0_shsg + 173);

    auto g_0_xyyyz_0_xzzz_0 = pbuffer.data(idx_eri_0_shsg + 174);

    auto g_0_xyyyz_0_yyyy_0 = pbuffer.data(idx_eri_0_shsg + 175);

    auto g_0_xyyyz_0_yyyz_0 = pbuffer.data(idx_eri_0_shsg + 176);

    auto g_0_xyyyz_0_yyzz_0 = pbuffer.data(idx_eri_0_shsg + 177);

    auto g_0_xyyyz_0_yzzz_0 = pbuffer.data(idx_eri_0_shsg + 178);

    auto g_0_xyyyz_0_zzzz_0 = pbuffer.data(idx_eri_0_shsg + 179);

#pragma omp simd aligned(g_0_xyyy_0_xxxx_0,      \
                             g_0_xyyy_0_xxxx_1,  \
                             g_0_xyyy_0_xxxy_0,  \
                             g_0_xyyy_0_xxxy_1,  \
                             g_0_xyyy_0_xxyy_0,  \
                             g_0_xyyy_0_xxyy_1,  \
                             g_0_xyyy_0_xyyy_0,  \
                             g_0_xyyy_0_xyyy_1,  \
                             g_0_xyyyz_0_xxxx_0, \
                             g_0_xyyyz_0_xxxy_0, \
                             g_0_xyyyz_0_xxxz_0, \
                             g_0_xyyyz_0_xxyy_0, \
                             g_0_xyyyz_0_xxyz_0, \
                             g_0_xyyyz_0_xxzz_0, \
                             g_0_xyyyz_0_xyyy_0, \
                             g_0_xyyyz_0_xyyz_0, \
                             g_0_xyyyz_0_xyzz_0, \
                             g_0_xyyyz_0_xzzz_0, \
                             g_0_xyyyz_0_yyyy_0, \
                             g_0_xyyyz_0_yyyz_0, \
                             g_0_xyyyz_0_yyzz_0, \
                             g_0_xyyyz_0_yzzz_0, \
                             g_0_xyyyz_0_zzzz_0, \
                             g_0_yyyz_0_xxxz_0,  \
                             g_0_yyyz_0_xxxz_1,  \
                             g_0_yyyz_0_xxyz_0,  \
                             g_0_yyyz_0_xxyz_1,  \
                             g_0_yyyz_0_xxz_1,   \
                             g_0_yyyz_0_xxzz_0,  \
                             g_0_yyyz_0_xxzz_1,  \
                             g_0_yyyz_0_xyyz_0,  \
                             g_0_yyyz_0_xyyz_1,  \
                             g_0_yyyz_0_xyz_1,   \
                             g_0_yyyz_0_xyzz_0,  \
                             g_0_yyyz_0_xyzz_1,  \
                             g_0_yyyz_0_xzz_1,   \
                             g_0_yyyz_0_xzzz_0,  \
                             g_0_yyyz_0_xzzz_1,  \
                             g_0_yyyz_0_yyyy_0,  \
                             g_0_yyyz_0_yyyy_1,  \
                             g_0_yyyz_0_yyyz_0,  \
                             g_0_yyyz_0_yyyz_1,  \
                             g_0_yyyz_0_yyz_1,   \
                             g_0_yyyz_0_yyzz_0,  \
                             g_0_yyyz_0_yyzz_1,  \
                             g_0_yyyz_0_yzz_1,   \
                             g_0_yyyz_0_yzzz_0,  \
                             g_0_yyyz_0_yzzz_1,  \
                             g_0_yyyz_0_zzz_1,   \
                             g_0_yyyz_0_zzzz_0,  \
                             g_0_yyyz_0_zzzz_1,  \
                             wp_x,               \
                             wp_z,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyz_0_xxxx_0[i] = g_0_xyyy_0_xxxx_0[i] * pb_z + g_0_xyyy_0_xxxx_1[i] * wp_z[i];

        g_0_xyyyz_0_xxxy_0[i] = g_0_xyyy_0_xxxy_0[i] * pb_z + g_0_xyyy_0_xxxy_1[i] * wp_z[i];

        g_0_xyyyz_0_xxxz_0[i] = 3.0 * g_0_yyyz_0_xxz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxxz_0[i] * pb_x + g_0_yyyz_0_xxxz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxyy_0[i] = g_0_xyyy_0_xxyy_0[i] * pb_z + g_0_xyyy_0_xxyy_1[i] * wp_z[i];

        g_0_xyyyz_0_xxyz_0[i] = 2.0 * g_0_yyyz_0_xyz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxyz_0[i] * pb_x + g_0_yyyz_0_xxyz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxzz_0[i] = 2.0 * g_0_yyyz_0_xzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxzz_0[i] * pb_x + g_0_yyyz_0_xxzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xyyy_0[i] = g_0_xyyy_0_xyyy_0[i] * pb_z + g_0_xyyy_0_xyyy_1[i] * wp_z[i];

        g_0_xyyyz_0_xyyz_0[i] = g_0_yyyz_0_yyz_1[i] * fi_abcd_0 + g_0_yyyz_0_xyyz_0[i] * pb_x + g_0_yyyz_0_xyyz_1[i] * wp_x[i];

        g_0_xyyyz_0_xyzz_0[i] = g_0_yyyz_0_yzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xyzz_0[i] * pb_x + g_0_yyyz_0_xyzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xzzz_0[i] = g_0_yyyz_0_zzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xzzz_0[i] * pb_x + g_0_yyyz_0_xzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_yyyy_0[i] = g_0_yyyz_0_yyyy_0[i] * pb_x + g_0_yyyz_0_yyyy_1[i] * wp_x[i];

        g_0_xyyyz_0_yyyz_0[i] = g_0_yyyz_0_yyyz_0[i] * pb_x + g_0_yyyz_0_yyyz_1[i] * wp_x[i];

        g_0_xyyyz_0_yyzz_0[i] = g_0_yyyz_0_yyzz_0[i] * pb_x + g_0_yyyz_0_yyzz_1[i] * wp_x[i];

        g_0_xyyyz_0_yzzz_0[i] = g_0_yyyz_0_yzzz_0[i] * pb_x + g_0_yyyz_0_yzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_zzzz_0[i] = g_0_yyyz_0_zzzz_0[i] * pb_x + g_0_yyyz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 180-195 components of targeted buffer : SHSG

    auto g_0_xyyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_shsg + 180);

    auto g_0_xyyzz_0_xxxy_0 = pbuffer.data(idx_eri_0_shsg + 181);

    auto g_0_xyyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_shsg + 182);

    auto g_0_xyyzz_0_xxyy_0 = pbuffer.data(idx_eri_0_shsg + 183);

    auto g_0_xyyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_shsg + 184);

    auto g_0_xyyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_shsg + 185);

    auto g_0_xyyzz_0_xyyy_0 = pbuffer.data(idx_eri_0_shsg + 186);

    auto g_0_xyyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_shsg + 187);

    auto g_0_xyyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_shsg + 188);

    auto g_0_xyyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_shsg + 189);

    auto g_0_xyyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_shsg + 190);

    auto g_0_xyyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_shsg + 191);

    auto g_0_xyyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_shsg + 192);

    auto g_0_xyyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_shsg + 193);

    auto g_0_xyyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_shsg + 194);

#pragma omp simd aligned(g_0_xyyzz_0_xxxx_0,     \
                             g_0_xyyzz_0_xxxy_0, \
                             g_0_xyyzz_0_xxxz_0, \
                             g_0_xyyzz_0_xxyy_0, \
                             g_0_xyyzz_0_xxyz_0, \
                             g_0_xyyzz_0_xxzz_0, \
                             g_0_xyyzz_0_xyyy_0, \
                             g_0_xyyzz_0_xyyz_0, \
                             g_0_xyyzz_0_xyzz_0, \
                             g_0_xyyzz_0_xzzz_0, \
                             g_0_xyyzz_0_yyyy_0, \
                             g_0_xyyzz_0_yyyz_0, \
                             g_0_xyyzz_0_yyzz_0, \
                             g_0_xyyzz_0_yzzz_0, \
                             g_0_xyyzz_0_zzzz_0, \
                             g_0_yyzz_0_xxx_1,   \
                             g_0_yyzz_0_xxxx_0,  \
                             g_0_yyzz_0_xxxx_1,  \
                             g_0_yyzz_0_xxxy_0,  \
                             g_0_yyzz_0_xxxy_1,  \
                             g_0_yyzz_0_xxxz_0,  \
                             g_0_yyzz_0_xxxz_1,  \
                             g_0_yyzz_0_xxy_1,   \
                             g_0_yyzz_0_xxyy_0,  \
                             g_0_yyzz_0_xxyy_1,  \
                             g_0_yyzz_0_xxyz_0,  \
                             g_0_yyzz_0_xxyz_1,  \
                             g_0_yyzz_0_xxz_1,   \
                             g_0_yyzz_0_xxzz_0,  \
                             g_0_yyzz_0_xxzz_1,  \
                             g_0_yyzz_0_xyy_1,   \
                             g_0_yyzz_0_xyyy_0,  \
                             g_0_yyzz_0_xyyy_1,  \
                             g_0_yyzz_0_xyyz_0,  \
                             g_0_yyzz_0_xyyz_1,  \
                             g_0_yyzz_0_xyz_1,   \
                             g_0_yyzz_0_xyzz_0,  \
                             g_0_yyzz_0_xyzz_1,  \
                             g_0_yyzz_0_xzz_1,   \
                             g_0_yyzz_0_xzzz_0,  \
                             g_0_yyzz_0_xzzz_1,  \
                             g_0_yyzz_0_yyy_1,   \
                             g_0_yyzz_0_yyyy_0,  \
                             g_0_yyzz_0_yyyy_1,  \
                             g_0_yyzz_0_yyyz_0,  \
                             g_0_yyzz_0_yyyz_1,  \
                             g_0_yyzz_0_yyz_1,   \
                             g_0_yyzz_0_yyzz_0,  \
                             g_0_yyzz_0_yyzz_1,  \
                             g_0_yyzz_0_yzz_1,   \
                             g_0_yyzz_0_yzzz_0,  \
                             g_0_yyzz_0_yzzz_1,  \
                             g_0_yyzz_0_zzz_1,   \
                             g_0_yyzz_0_zzzz_0,  \
                             g_0_yyzz_0_zzzz_1,  \
                             wp_x,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyzz_0_xxxx_0[i] = 4.0 * g_0_yyzz_0_xxx_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxx_0[i] * pb_x + g_0_yyzz_0_xxxx_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxy_0[i] = 3.0 * g_0_yyzz_0_xxy_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxy_0[i] * pb_x + g_0_yyzz_0_xxxy_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxz_0[i] = 3.0 * g_0_yyzz_0_xxz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxz_0[i] * pb_x + g_0_yyzz_0_xxxz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxyy_0[i] = 2.0 * g_0_yyzz_0_xyy_1[i] * fi_abcd_0 + g_0_yyzz_0_xxyy_0[i] * pb_x + g_0_yyzz_0_xxyy_1[i] * wp_x[i];

        g_0_xyyzz_0_xxyz_0[i] = 2.0 * g_0_yyzz_0_xyz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxyz_0[i] * pb_x + g_0_yyzz_0_xxyz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxzz_0[i] = 2.0 * g_0_yyzz_0_xzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxzz_0[i] * pb_x + g_0_yyzz_0_xxzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xyyy_0[i] = g_0_yyzz_0_yyy_1[i] * fi_abcd_0 + g_0_yyzz_0_xyyy_0[i] * pb_x + g_0_yyzz_0_xyyy_1[i] * wp_x[i];

        g_0_xyyzz_0_xyyz_0[i] = g_0_yyzz_0_yyz_1[i] * fi_abcd_0 + g_0_yyzz_0_xyyz_0[i] * pb_x + g_0_yyzz_0_xyyz_1[i] * wp_x[i];

        g_0_xyyzz_0_xyzz_0[i] = g_0_yyzz_0_yzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xyzz_0[i] * pb_x + g_0_yyzz_0_xyzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xzzz_0[i] = g_0_yyzz_0_zzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xzzz_0[i] * pb_x + g_0_yyzz_0_xzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_yyyy_0[i] = g_0_yyzz_0_yyyy_0[i] * pb_x + g_0_yyzz_0_yyyy_1[i] * wp_x[i];

        g_0_xyyzz_0_yyyz_0[i] = g_0_yyzz_0_yyyz_0[i] * pb_x + g_0_yyzz_0_yyyz_1[i] * wp_x[i];

        g_0_xyyzz_0_yyzz_0[i] = g_0_yyzz_0_yyzz_0[i] * pb_x + g_0_yyzz_0_yyzz_1[i] * wp_x[i];

        g_0_xyyzz_0_yzzz_0[i] = g_0_yyzz_0_yzzz_0[i] * pb_x + g_0_yyzz_0_yzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_zzzz_0[i] = g_0_yyzz_0_zzzz_0[i] * pb_x + g_0_yyzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 195-210 components of targeted buffer : SHSG

    auto g_0_xyzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_shsg + 195);

    auto g_0_xyzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_shsg + 196);

    auto g_0_xyzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_shsg + 197);

    auto g_0_xyzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_shsg + 198);

    auto g_0_xyzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_shsg + 199);

    auto g_0_xyzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_shsg + 200);

    auto g_0_xyzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_shsg + 201);

    auto g_0_xyzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_shsg + 202);

    auto g_0_xyzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_shsg + 203);

    auto g_0_xyzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_shsg + 204);

    auto g_0_xyzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_shsg + 205);

    auto g_0_xyzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_shsg + 206);

    auto g_0_xyzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_shsg + 207);

    auto g_0_xyzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_shsg + 208);

    auto g_0_xyzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_shsg + 209);

#pragma omp simd aligned(g_0_xyzzz_0_xxxx_0,     \
                             g_0_xyzzz_0_xxxy_0, \
                             g_0_xyzzz_0_xxxz_0, \
                             g_0_xyzzz_0_xxyy_0, \
                             g_0_xyzzz_0_xxyz_0, \
                             g_0_xyzzz_0_xxzz_0, \
                             g_0_xyzzz_0_xyyy_0, \
                             g_0_xyzzz_0_xyyz_0, \
                             g_0_xyzzz_0_xyzz_0, \
                             g_0_xyzzz_0_xzzz_0, \
                             g_0_xyzzz_0_yyyy_0, \
                             g_0_xyzzz_0_yyyz_0, \
                             g_0_xyzzz_0_yyzz_0, \
                             g_0_xyzzz_0_yzzz_0, \
                             g_0_xyzzz_0_zzzz_0, \
                             g_0_xzzz_0_xxxx_0,  \
                             g_0_xzzz_0_xxxx_1,  \
                             g_0_xzzz_0_xxxz_0,  \
                             g_0_xzzz_0_xxxz_1,  \
                             g_0_xzzz_0_xxzz_0,  \
                             g_0_xzzz_0_xxzz_1,  \
                             g_0_xzzz_0_xzzz_0,  \
                             g_0_xzzz_0_xzzz_1,  \
                             g_0_yzzz_0_xxxy_0,  \
                             g_0_yzzz_0_xxxy_1,  \
                             g_0_yzzz_0_xxy_1,   \
                             g_0_yzzz_0_xxyy_0,  \
                             g_0_yzzz_0_xxyy_1,  \
                             g_0_yzzz_0_xxyz_0,  \
                             g_0_yzzz_0_xxyz_1,  \
                             g_0_yzzz_0_xyy_1,   \
                             g_0_yzzz_0_xyyy_0,  \
                             g_0_yzzz_0_xyyy_1,  \
                             g_0_yzzz_0_xyyz_0,  \
                             g_0_yzzz_0_xyyz_1,  \
                             g_0_yzzz_0_xyz_1,   \
                             g_0_yzzz_0_xyzz_0,  \
                             g_0_yzzz_0_xyzz_1,  \
                             g_0_yzzz_0_yyy_1,   \
                             g_0_yzzz_0_yyyy_0,  \
                             g_0_yzzz_0_yyyy_1,  \
                             g_0_yzzz_0_yyyz_0,  \
                             g_0_yzzz_0_yyyz_1,  \
                             g_0_yzzz_0_yyz_1,   \
                             g_0_yzzz_0_yyzz_0,  \
                             g_0_yzzz_0_yyzz_1,  \
                             g_0_yzzz_0_yzz_1,   \
                             g_0_yzzz_0_yzzz_0,  \
                             g_0_yzzz_0_yzzz_1,  \
                             g_0_yzzz_0_zzzz_0,  \
                             g_0_yzzz_0_zzzz_1,  \
                             wp_x,               \
                             wp_y,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyzzz_0_xxxx_0[i] = g_0_xzzz_0_xxxx_0[i] * pb_y + g_0_xzzz_0_xxxx_1[i] * wp_y[i];

        g_0_xyzzz_0_xxxy_0[i] = 3.0 * g_0_yzzz_0_xxy_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxy_0[i] * pb_x + g_0_yzzz_0_xxxy_1[i] * wp_x[i];

        g_0_xyzzz_0_xxxz_0[i] = g_0_xzzz_0_xxxz_0[i] * pb_y + g_0_xzzz_0_xxxz_1[i] * wp_y[i];

        g_0_xyzzz_0_xxyy_0[i] = 2.0 * g_0_yzzz_0_xyy_1[i] * fi_abcd_0 + g_0_yzzz_0_xxyy_0[i] * pb_x + g_0_yzzz_0_xxyy_1[i] * wp_x[i];

        g_0_xyzzz_0_xxyz_0[i] = 2.0 * g_0_yzzz_0_xyz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxyz_0[i] * pb_x + g_0_yzzz_0_xxyz_1[i] * wp_x[i];

        g_0_xyzzz_0_xxzz_0[i] = g_0_xzzz_0_xxzz_0[i] * pb_y + g_0_xzzz_0_xxzz_1[i] * wp_y[i];

        g_0_xyzzz_0_xyyy_0[i] = g_0_yzzz_0_yyy_1[i] * fi_abcd_0 + g_0_yzzz_0_xyyy_0[i] * pb_x + g_0_yzzz_0_xyyy_1[i] * wp_x[i];

        g_0_xyzzz_0_xyyz_0[i] = g_0_yzzz_0_yyz_1[i] * fi_abcd_0 + g_0_yzzz_0_xyyz_0[i] * pb_x + g_0_yzzz_0_xyyz_1[i] * wp_x[i];

        g_0_xyzzz_0_xyzz_0[i] = g_0_yzzz_0_yzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xyzz_0[i] * pb_x + g_0_yzzz_0_xyzz_1[i] * wp_x[i];

        g_0_xyzzz_0_xzzz_0[i] = g_0_xzzz_0_xzzz_0[i] * pb_y + g_0_xzzz_0_xzzz_1[i] * wp_y[i];

        g_0_xyzzz_0_yyyy_0[i] = g_0_yzzz_0_yyyy_0[i] * pb_x + g_0_yzzz_0_yyyy_1[i] * wp_x[i];

        g_0_xyzzz_0_yyyz_0[i] = g_0_yzzz_0_yyyz_0[i] * pb_x + g_0_yzzz_0_yyyz_1[i] * wp_x[i];

        g_0_xyzzz_0_yyzz_0[i] = g_0_yzzz_0_yyzz_0[i] * pb_x + g_0_yzzz_0_yyzz_1[i] * wp_x[i];

        g_0_xyzzz_0_yzzz_0[i] = g_0_yzzz_0_yzzz_0[i] * pb_x + g_0_yzzz_0_yzzz_1[i] * wp_x[i];

        g_0_xyzzz_0_zzzz_0[i] = g_0_yzzz_0_zzzz_0[i] * pb_x + g_0_yzzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 210-225 components of targeted buffer : SHSG

    auto g_0_xzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_shsg + 210);

    auto g_0_xzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_shsg + 211);

    auto g_0_xzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_shsg + 212);

    auto g_0_xzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_shsg + 213);

    auto g_0_xzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_shsg + 214);

    auto g_0_xzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_shsg + 215);

    auto g_0_xzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_shsg + 216);

    auto g_0_xzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_shsg + 217);

    auto g_0_xzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_shsg + 218);

    auto g_0_xzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_shsg + 219);

    auto g_0_xzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_shsg + 220);

    auto g_0_xzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_shsg + 221);

    auto g_0_xzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_shsg + 222);

    auto g_0_xzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_shsg + 223);

    auto g_0_xzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_shsg + 224);

#pragma omp simd aligned(g_0_xzzzz_0_xxxx_0,     \
                             g_0_xzzzz_0_xxxy_0, \
                             g_0_xzzzz_0_xxxz_0, \
                             g_0_xzzzz_0_xxyy_0, \
                             g_0_xzzzz_0_xxyz_0, \
                             g_0_xzzzz_0_xxzz_0, \
                             g_0_xzzzz_0_xyyy_0, \
                             g_0_xzzzz_0_xyyz_0, \
                             g_0_xzzzz_0_xyzz_0, \
                             g_0_xzzzz_0_xzzz_0, \
                             g_0_xzzzz_0_yyyy_0, \
                             g_0_xzzzz_0_yyyz_0, \
                             g_0_xzzzz_0_yyzz_0, \
                             g_0_xzzzz_0_yzzz_0, \
                             g_0_xzzzz_0_zzzz_0, \
                             g_0_zzzz_0_xxx_1,   \
                             g_0_zzzz_0_xxxx_0,  \
                             g_0_zzzz_0_xxxx_1,  \
                             g_0_zzzz_0_xxxy_0,  \
                             g_0_zzzz_0_xxxy_1,  \
                             g_0_zzzz_0_xxxz_0,  \
                             g_0_zzzz_0_xxxz_1,  \
                             g_0_zzzz_0_xxy_1,   \
                             g_0_zzzz_0_xxyy_0,  \
                             g_0_zzzz_0_xxyy_1,  \
                             g_0_zzzz_0_xxyz_0,  \
                             g_0_zzzz_0_xxyz_1,  \
                             g_0_zzzz_0_xxz_1,   \
                             g_0_zzzz_0_xxzz_0,  \
                             g_0_zzzz_0_xxzz_1,  \
                             g_0_zzzz_0_xyy_1,   \
                             g_0_zzzz_0_xyyy_0,  \
                             g_0_zzzz_0_xyyy_1,  \
                             g_0_zzzz_0_xyyz_0,  \
                             g_0_zzzz_0_xyyz_1,  \
                             g_0_zzzz_0_xyz_1,   \
                             g_0_zzzz_0_xyzz_0,  \
                             g_0_zzzz_0_xyzz_1,  \
                             g_0_zzzz_0_xzz_1,   \
                             g_0_zzzz_0_xzzz_0,  \
                             g_0_zzzz_0_xzzz_1,  \
                             g_0_zzzz_0_yyy_1,   \
                             g_0_zzzz_0_yyyy_0,  \
                             g_0_zzzz_0_yyyy_1,  \
                             g_0_zzzz_0_yyyz_0,  \
                             g_0_zzzz_0_yyyz_1,  \
                             g_0_zzzz_0_yyz_1,   \
                             g_0_zzzz_0_yyzz_0,  \
                             g_0_zzzz_0_yyzz_1,  \
                             g_0_zzzz_0_yzz_1,   \
                             g_0_zzzz_0_yzzz_0,  \
                             g_0_zzzz_0_yzzz_1,  \
                             g_0_zzzz_0_zzz_1,   \
                             g_0_zzzz_0_zzzz_0,  \
                             g_0_zzzz_0_zzzz_1,  \
                             wp_x,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzzzz_0_xxxx_0[i] = 4.0 * g_0_zzzz_0_xxx_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxx_0[i] * pb_x + g_0_zzzz_0_xxxx_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxy_0[i] = 3.0 * g_0_zzzz_0_xxy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxy_0[i] * pb_x + g_0_zzzz_0_xxxy_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxz_0[i] = 3.0 * g_0_zzzz_0_xxz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxz_0[i] * pb_x + g_0_zzzz_0_xxxz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxyy_0[i] = 2.0 * g_0_zzzz_0_xyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyy_0[i] * pb_x + g_0_zzzz_0_xxyy_1[i] * wp_x[i];

        g_0_xzzzz_0_xxyz_0[i] = 2.0 * g_0_zzzz_0_xyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyz_0[i] * pb_x + g_0_zzzz_0_xxyz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxzz_0[i] = 2.0 * g_0_zzzz_0_xzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxzz_0[i] * pb_x + g_0_zzzz_0_xxzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xyyy_0[i] = g_0_zzzz_0_yyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyy_0[i] * pb_x + g_0_zzzz_0_xyyy_1[i] * wp_x[i];

        g_0_xzzzz_0_xyyz_0[i] = g_0_zzzz_0_yyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyz_0[i] * pb_x + g_0_zzzz_0_xyyz_1[i] * wp_x[i];

        g_0_xzzzz_0_xyzz_0[i] = g_0_zzzz_0_yzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyzz_0[i] * pb_x + g_0_zzzz_0_xyzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xzzz_0[i] = g_0_zzzz_0_zzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xzzz_0[i] * pb_x + g_0_zzzz_0_xzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_yyyy_0[i] = g_0_zzzz_0_yyyy_0[i] * pb_x + g_0_zzzz_0_yyyy_1[i] * wp_x[i];

        g_0_xzzzz_0_yyyz_0[i] = g_0_zzzz_0_yyyz_0[i] * pb_x + g_0_zzzz_0_yyyz_1[i] * wp_x[i];

        g_0_xzzzz_0_yyzz_0[i] = g_0_zzzz_0_yyzz_0[i] * pb_x + g_0_zzzz_0_yyzz_1[i] * wp_x[i];

        g_0_xzzzz_0_yzzz_0[i] = g_0_zzzz_0_yzzz_0[i] * pb_x + g_0_zzzz_0_yzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_zzzz_0[i] = g_0_zzzz_0_zzzz_0[i] * pb_x + g_0_zzzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 225-240 components of targeted buffer : SHSG

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

#pragma omp simd aligned(g_0_yyy_0_xxxx_0,       \
                             g_0_yyy_0_xxxx_1,   \
                             g_0_yyy_0_xxxy_0,   \
                             g_0_yyy_0_xxxy_1,   \
                             g_0_yyy_0_xxxz_0,   \
                             g_0_yyy_0_xxxz_1,   \
                             g_0_yyy_0_xxyy_0,   \
                             g_0_yyy_0_xxyy_1,   \
                             g_0_yyy_0_xxyz_0,   \
                             g_0_yyy_0_xxyz_1,   \
                             g_0_yyy_0_xxzz_0,   \
                             g_0_yyy_0_xxzz_1,   \
                             g_0_yyy_0_xyyy_0,   \
                             g_0_yyy_0_xyyy_1,   \
                             g_0_yyy_0_xyyz_0,   \
                             g_0_yyy_0_xyyz_1,   \
                             g_0_yyy_0_xyzz_0,   \
                             g_0_yyy_0_xyzz_1,   \
                             g_0_yyy_0_xzzz_0,   \
                             g_0_yyy_0_xzzz_1,   \
                             g_0_yyy_0_yyyy_0,   \
                             g_0_yyy_0_yyyy_1,   \
                             g_0_yyy_0_yyyz_0,   \
                             g_0_yyy_0_yyyz_1,   \
                             g_0_yyy_0_yyzz_0,   \
                             g_0_yyy_0_yyzz_1,   \
                             g_0_yyy_0_yzzz_0,   \
                             g_0_yyy_0_yzzz_1,   \
                             g_0_yyy_0_zzzz_0,   \
                             g_0_yyy_0_zzzz_1,   \
                             g_0_yyyy_0_xxx_1,   \
                             g_0_yyyy_0_xxxx_0,  \
                             g_0_yyyy_0_xxxx_1,  \
                             g_0_yyyy_0_xxxy_0,  \
                             g_0_yyyy_0_xxxy_1,  \
                             g_0_yyyy_0_xxxz_0,  \
                             g_0_yyyy_0_xxxz_1,  \
                             g_0_yyyy_0_xxy_1,   \
                             g_0_yyyy_0_xxyy_0,  \
                             g_0_yyyy_0_xxyy_1,  \
                             g_0_yyyy_0_xxyz_0,  \
                             g_0_yyyy_0_xxyz_1,  \
                             g_0_yyyy_0_xxz_1,   \
                             g_0_yyyy_0_xxzz_0,  \
                             g_0_yyyy_0_xxzz_1,  \
                             g_0_yyyy_0_xyy_1,   \
                             g_0_yyyy_0_xyyy_0,  \
                             g_0_yyyy_0_xyyy_1,  \
                             g_0_yyyy_0_xyyz_0,  \
                             g_0_yyyy_0_xyyz_1,  \
                             g_0_yyyy_0_xyz_1,   \
                             g_0_yyyy_0_xyzz_0,  \
                             g_0_yyyy_0_xyzz_1,  \
                             g_0_yyyy_0_xzz_1,   \
                             g_0_yyyy_0_xzzz_0,  \
                             g_0_yyyy_0_xzzz_1,  \
                             g_0_yyyy_0_yyy_1,   \
                             g_0_yyyy_0_yyyy_0,  \
                             g_0_yyyy_0_yyyy_1,  \
                             g_0_yyyy_0_yyyz_0,  \
                             g_0_yyyy_0_yyyz_1,  \
                             g_0_yyyy_0_yyz_1,   \
                             g_0_yyyy_0_yyzz_0,  \
                             g_0_yyyy_0_yyzz_1,  \
                             g_0_yyyy_0_yzz_1,   \
                             g_0_yyyy_0_yzzz_0,  \
                             g_0_yyyy_0_yzzz_1,  \
                             g_0_yyyy_0_zzz_1,   \
                             g_0_yyyy_0_zzzz_0,  \
                             g_0_yyyy_0_zzzz_1,  \
                             g_0_yyyyy_0_xxxx_0, \
                             g_0_yyyyy_0_xxxy_0, \
                             g_0_yyyyy_0_xxxz_0, \
                             g_0_yyyyy_0_xxyy_0, \
                             g_0_yyyyy_0_xxyz_0, \
                             g_0_yyyyy_0_xxzz_0, \
                             g_0_yyyyy_0_xyyy_0, \
                             g_0_yyyyy_0_xyyz_0, \
                             g_0_yyyyy_0_xyzz_0, \
                             g_0_yyyyy_0_xzzz_0, \
                             g_0_yyyyy_0_yyyy_0, \
                             g_0_yyyyy_0_yyyz_0, \
                             g_0_yyyyy_0_yyzz_0, \
                             g_0_yyyyy_0_yzzz_0, \
                             g_0_yyyyy_0_zzzz_0, \
                             wp_y,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyy_0_xxxx_0[i] =
            4.0 * g_0_yyy_0_xxxx_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxx_1[i] * fti_ab_0 + g_0_yyyy_0_xxxx_0[i] * pb_y + g_0_yyyy_0_xxxx_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxy_0[i] = 4.0 * g_0_yyy_0_xxxy_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxy_1[i] * fti_ab_0 + g_0_yyyy_0_xxx_1[i] * fi_abcd_0 +
                                g_0_yyyy_0_xxxy_0[i] * pb_y + g_0_yyyy_0_xxxy_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxz_0[i] =
            4.0 * g_0_yyy_0_xxxz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxz_1[i] * fti_ab_0 + g_0_yyyy_0_xxxz_0[i] * pb_y + g_0_yyyy_0_xxxz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxyy_0[i] = 4.0 * g_0_yyy_0_xxyy_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxyy_1[i] * fti_ab_0 + 2.0 * g_0_yyyy_0_xxy_1[i] * fi_abcd_0 +
                                g_0_yyyy_0_xxyy_0[i] * pb_y + g_0_yyyy_0_xxyy_1[i] * wp_y[i];

        g_0_yyyyy_0_xxyz_0[i] = 4.0 * g_0_yyy_0_xxyz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxyz_1[i] * fti_ab_0 + g_0_yyyy_0_xxz_1[i] * fi_abcd_0 +
                                g_0_yyyy_0_xxyz_0[i] * pb_y + g_0_yyyy_0_xxyz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxzz_0[i] =
            4.0 * g_0_yyy_0_xxzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxzz_1[i] * fti_ab_0 + g_0_yyyy_0_xxzz_0[i] * pb_y + g_0_yyyy_0_xxzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xyyy_0[i] = 4.0 * g_0_yyy_0_xyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xyyy_1[i] * fti_ab_0 + 3.0 * g_0_yyyy_0_xyy_1[i] * fi_abcd_0 +
                                g_0_yyyy_0_xyyy_0[i] * pb_y + g_0_yyyy_0_xyyy_1[i] * wp_y[i];

        g_0_yyyyy_0_xyyz_0[i] = 4.0 * g_0_yyy_0_xyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xyyz_1[i] * fti_ab_0 + 2.0 * g_0_yyyy_0_xyz_1[i] * fi_abcd_0 +
                                g_0_yyyy_0_xyyz_0[i] * pb_y + g_0_yyyy_0_xyyz_1[i] * wp_y[i];

        g_0_yyyyy_0_xyzz_0[i] = 4.0 * g_0_yyy_0_xyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xyzz_1[i] * fti_ab_0 + g_0_yyyy_0_xzz_1[i] * fi_abcd_0 +
                                g_0_yyyy_0_xyzz_0[i] * pb_y + g_0_yyyy_0_xyzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xzzz_0[i] =
            4.0 * g_0_yyy_0_xzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xzzz_1[i] * fti_ab_0 + g_0_yyyy_0_xzzz_0[i] * pb_y + g_0_yyyy_0_xzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_yyyy_0[i] = 4.0 * g_0_yyy_0_yyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_yyyy_1[i] * fti_ab_0 + 4.0 * g_0_yyyy_0_yyy_1[i] * fi_abcd_0 +
                                g_0_yyyy_0_yyyy_0[i] * pb_y + g_0_yyyy_0_yyyy_1[i] * wp_y[i];

        g_0_yyyyy_0_yyyz_0[i] = 4.0 * g_0_yyy_0_yyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_yyyz_1[i] * fti_ab_0 + 3.0 * g_0_yyyy_0_yyz_1[i] * fi_abcd_0 +
                                g_0_yyyy_0_yyyz_0[i] * pb_y + g_0_yyyy_0_yyyz_1[i] * wp_y[i];

        g_0_yyyyy_0_yyzz_0[i] = 4.0 * g_0_yyy_0_yyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_yyzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyy_0_yzz_1[i] * fi_abcd_0 +
                                g_0_yyyy_0_yyzz_0[i] * pb_y + g_0_yyyy_0_yyzz_1[i] * wp_y[i];

        g_0_yyyyy_0_yzzz_0[i] = 4.0 * g_0_yyy_0_yzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_yzzz_1[i] * fti_ab_0 + g_0_yyyy_0_zzz_1[i] * fi_abcd_0 +
                                g_0_yyyy_0_yzzz_0[i] * pb_y + g_0_yyyy_0_yzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_zzzz_0[i] =
            4.0 * g_0_yyy_0_zzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_zzzz_1[i] * fti_ab_0 + g_0_yyyy_0_zzzz_0[i] * pb_y + g_0_yyyy_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 240-255 components of targeted buffer : SHSG

    auto g_0_yyyyz_0_xxxx_0 = pbuffer.data(idx_eri_0_shsg + 240);

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

#pragma omp simd aligned(g_0_yyyy_0_xxx_1,       \
                             g_0_yyyy_0_xxxx_0,  \
                             g_0_yyyy_0_xxxx_1,  \
                             g_0_yyyy_0_xxxy_0,  \
                             g_0_yyyy_0_xxxy_1,  \
                             g_0_yyyy_0_xxxz_0,  \
                             g_0_yyyy_0_xxxz_1,  \
                             g_0_yyyy_0_xxy_1,   \
                             g_0_yyyy_0_xxyy_0,  \
                             g_0_yyyy_0_xxyy_1,  \
                             g_0_yyyy_0_xxyz_0,  \
                             g_0_yyyy_0_xxyz_1,  \
                             g_0_yyyy_0_xxz_1,   \
                             g_0_yyyy_0_xxzz_0,  \
                             g_0_yyyy_0_xxzz_1,  \
                             g_0_yyyy_0_xyy_1,   \
                             g_0_yyyy_0_xyyy_0,  \
                             g_0_yyyy_0_xyyy_1,  \
                             g_0_yyyy_0_xyyz_0,  \
                             g_0_yyyy_0_xyyz_1,  \
                             g_0_yyyy_0_xyz_1,   \
                             g_0_yyyy_0_xyzz_0,  \
                             g_0_yyyy_0_xyzz_1,  \
                             g_0_yyyy_0_xzz_1,   \
                             g_0_yyyy_0_xzzz_0,  \
                             g_0_yyyy_0_xzzz_1,  \
                             g_0_yyyy_0_yyy_1,   \
                             g_0_yyyy_0_yyyy_0,  \
                             g_0_yyyy_0_yyyy_1,  \
                             g_0_yyyy_0_yyyz_0,  \
                             g_0_yyyy_0_yyyz_1,  \
                             g_0_yyyy_0_yyz_1,   \
                             g_0_yyyy_0_yyzz_0,  \
                             g_0_yyyy_0_yyzz_1,  \
                             g_0_yyyy_0_yzz_1,   \
                             g_0_yyyy_0_yzzz_0,  \
                             g_0_yyyy_0_yzzz_1,  \
                             g_0_yyyy_0_zzz_1,   \
                             g_0_yyyy_0_zzzz_0,  \
                             g_0_yyyy_0_zzzz_1,  \
                             g_0_yyyyz_0_xxxx_0, \
                             g_0_yyyyz_0_xxxy_0, \
                             g_0_yyyyz_0_xxxz_0, \
                             g_0_yyyyz_0_xxyy_0, \
                             g_0_yyyyz_0_xxyz_0, \
                             g_0_yyyyz_0_xxzz_0, \
                             g_0_yyyyz_0_xyyy_0, \
                             g_0_yyyyz_0_xyyz_0, \
                             g_0_yyyyz_0_xyzz_0, \
                             g_0_yyyyz_0_xzzz_0, \
                             g_0_yyyyz_0_yyyy_0, \
                             g_0_yyyyz_0_yyyz_0, \
                             g_0_yyyyz_0_yyzz_0, \
                             g_0_yyyyz_0_yzzz_0, \
                             g_0_yyyyz_0_zzzz_0, \
                             wp_z,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyyz_0_xxxx_0[i] = g_0_yyyy_0_xxxx_0[i] * pb_z + g_0_yyyy_0_xxxx_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxy_0[i] = g_0_yyyy_0_xxxy_0[i] * pb_z + g_0_yyyy_0_xxxy_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxz_0[i] = g_0_yyyy_0_xxx_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxz_0[i] * pb_z + g_0_yyyy_0_xxxz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxyy_0[i] = g_0_yyyy_0_xxyy_0[i] * pb_z + g_0_yyyy_0_xxyy_1[i] * wp_z[i];

        g_0_yyyyz_0_xxyz_0[i] = g_0_yyyy_0_xxy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyz_0[i] * pb_z + g_0_yyyy_0_xxyz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxzz_0[i] = 2.0 * g_0_yyyy_0_xxz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxzz_0[i] * pb_z + g_0_yyyy_0_xxzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xyyy_0[i] = g_0_yyyy_0_xyyy_0[i] * pb_z + g_0_yyyy_0_xyyy_1[i] * wp_z[i];

        g_0_yyyyz_0_xyyz_0[i] = g_0_yyyy_0_xyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyz_0[i] * pb_z + g_0_yyyy_0_xyyz_1[i] * wp_z[i];

        g_0_yyyyz_0_xyzz_0[i] = 2.0 * g_0_yyyy_0_xyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyzz_0[i] * pb_z + g_0_yyyy_0_xyzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xzzz_0[i] = 3.0 * g_0_yyyy_0_xzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xzzz_0[i] * pb_z + g_0_yyyy_0_xzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_yyyy_0[i] = g_0_yyyy_0_yyyy_0[i] * pb_z + g_0_yyyy_0_yyyy_1[i] * wp_z[i];

        g_0_yyyyz_0_yyyz_0[i] = g_0_yyyy_0_yyy_1[i] * fi_abcd_0 + g_0_yyyy_0_yyyz_0[i] * pb_z + g_0_yyyy_0_yyyz_1[i] * wp_z[i];

        g_0_yyyyz_0_yyzz_0[i] = 2.0 * g_0_yyyy_0_yyz_1[i] * fi_abcd_0 + g_0_yyyy_0_yyzz_0[i] * pb_z + g_0_yyyy_0_yyzz_1[i] * wp_z[i];

        g_0_yyyyz_0_yzzz_0[i] = 3.0 * g_0_yyyy_0_yzz_1[i] * fi_abcd_0 + g_0_yyyy_0_yzzz_0[i] * pb_z + g_0_yyyy_0_yzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_zzzz_0[i] = 4.0 * g_0_yyyy_0_zzz_1[i] * fi_abcd_0 + g_0_yyyy_0_zzzz_0[i] * pb_z + g_0_yyyy_0_zzzz_1[i] * wp_z[i];
    }

    /// Set up 255-270 components of targeted buffer : SHSG

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

#pragma omp simd aligned(g_0_yyy_0_xxxy_0,       \
                             g_0_yyy_0_xxxy_1,   \
                             g_0_yyy_0_xxyy_0,   \
                             g_0_yyy_0_xxyy_1,   \
                             g_0_yyy_0_xyyy_0,   \
                             g_0_yyy_0_xyyy_1,   \
                             g_0_yyy_0_yyyy_0,   \
                             g_0_yyy_0_yyyy_1,   \
                             g_0_yyyz_0_xxxy_0,  \
                             g_0_yyyz_0_xxxy_1,  \
                             g_0_yyyz_0_xxyy_0,  \
                             g_0_yyyz_0_xxyy_1,  \
                             g_0_yyyz_0_xyyy_0,  \
                             g_0_yyyz_0_xyyy_1,  \
                             g_0_yyyz_0_yyyy_0,  \
                             g_0_yyyz_0_yyyy_1,  \
                             g_0_yyyzz_0_xxxx_0, \
                             g_0_yyyzz_0_xxxy_0, \
                             g_0_yyyzz_0_xxxz_0, \
                             g_0_yyyzz_0_xxyy_0, \
                             g_0_yyyzz_0_xxyz_0, \
                             g_0_yyyzz_0_xxzz_0, \
                             g_0_yyyzz_0_xyyy_0, \
                             g_0_yyyzz_0_xyyz_0, \
                             g_0_yyyzz_0_xyzz_0, \
                             g_0_yyyzz_0_xzzz_0, \
                             g_0_yyyzz_0_yyyy_0, \
                             g_0_yyyzz_0_yyyz_0, \
                             g_0_yyyzz_0_yyzz_0, \
                             g_0_yyyzz_0_yzzz_0, \
                             g_0_yyyzz_0_zzzz_0, \
                             g_0_yyzz_0_xxxx_0,  \
                             g_0_yyzz_0_xxxx_1,  \
                             g_0_yyzz_0_xxxz_0,  \
                             g_0_yyzz_0_xxxz_1,  \
                             g_0_yyzz_0_xxyz_0,  \
                             g_0_yyzz_0_xxyz_1,  \
                             g_0_yyzz_0_xxz_1,   \
                             g_0_yyzz_0_xxzz_0,  \
                             g_0_yyzz_0_xxzz_1,  \
                             g_0_yyzz_0_xyyz_0,  \
                             g_0_yyzz_0_xyyz_1,  \
                             g_0_yyzz_0_xyz_1,   \
                             g_0_yyzz_0_xyzz_0,  \
                             g_0_yyzz_0_xyzz_1,  \
                             g_0_yyzz_0_xzz_1,   \
                             g_0_yyzz_0_xzzz_0,  \
                             g_0_yyzz_0_xzzz_1,  \
                             g_0_yyzz_0_yyyz_0,  \
                             g_0_yyzz_0_yyyz_1,  \
                             g_0_yyzz_0_yyz_1,   \
                             g_0_yyzz_0_yyzz_0,  \
                             g_0_yyzz_0_yyzz_1,  \
                             g_0_yyzz_0_yzz_1,   \
                             g_0_yyzz_0_yzzz_0,  \
                             g_0_yyzz_0_yzzz_1,  \
                             g_0_yyzz_0_zzz_1,   \
                             g_0_yyzz_0_zzzz_0,  \
                             g_0_yyzz_0_zzzz_1,  \
                             g_0_yzz_0_xxxx_0,   \
                             g_0_yzz_0_xxxx_1,   \
                             g_0_yzz_0_xxxz_0,   \
                             g_0_yzz_0_xxxz_1,   \
                             g_0_yzz_0_xxyz_0,   \
                             g_0_yzz_0_xxyz_1,   \
                             g_0_yzz_0_xxzz_0,   \
                             g_0_yzz_0_xxzz_1,   \
                             g_0_yzz_0_xyyz_0,   \
                             g_0_yzz_0_xyyz_1,   \
                             g_0_yzz_0_xyzz_0,   \
                             g_0_yzz_0_xyzz_1,   \
                             g_0_yzz_0_xzzz_0,   \
                             g_0_yzz_0_xzzz_1,   \
                             g_0_yzz_0_yyyz_0,   \
                             g_0_yzz_0_yyyz_1,   \
                             g_0_yzz_0_yyzz_0,   \
                             g_0_yzz_0_yyzz_1,   \
                             g_0_yzz_0_yzzz_0,   \
                             g_0_yzz_0_yzzz_1,   \
                             g_0_yzz_0_zzzz_0,   \
                             g_0_yzz_0_zzzz_1,   \
                             wp_y,               \
                             wp_z,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyzz_0_xxxx_0[i] =
            2.0 * g_0_yzz_0_xxxx_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxxx_1[i] * fti_ab_0 + g_0_yyzz_0_xxxx_0[i] * pb_y + g_0_yyzz_0_xxxx_1[i] * wp_y[i];

        g_0_yyyzz_0_xxxy_0[i] =
            g_0_yyy_0_xxxy_0[i] * fi_ab_0 - g_0_yyy_0_xxxy_1[i] * fti_ab_0 + g_0_yyyz_0_xxxy_0[i] * pb_z + g_0_yyyz_0_xxxy_1[i] * wp_z[i];

        g_0_yyyzz_0_xxxz_0[i] =
            2.0 * g_0_yzz_0_xxxz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxxz_1[i] * fti_ab_0 + g_0_yyzz_0_xxxz_0[i] * pb_y + g_0_yyzz_0_xxxz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxyy_0[i] =
            g_0_yyy_0_xxyy_0[i] * fi_ab_0 - g_0_yyy_0_xxyy_1[i] * fti_ab_0 + g_0_yyyz_0_xxyy_0[i] * pb_z + g_0_yyyz_0_xxyy_1[i] * wp_z[i];

        g_0_yyyzz_0_xxyz_0[i] = 2.0 * g_0_yzz_0_xxyz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxyz_1[i] * fti_ab_0 + g_0_yyzz_0_xxz_1[i] * fi_abcd_0 +
                                g_0_yyzz_0_xxyz_0[i] * pb_y + g_0_yyzz_0_xxyz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxzz_0[i] =
            2.0 * g_0_yzz_0_xxzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxzz_1[i] * fti_ab_0 + g_0_yyzz_0_xxzz_0[i] * pb_y + g_0_yyzz_0_xxzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xyyy_0[i] =
            g_0_yyy_0_xyyy_0[i] * fi_ab_0 - g_0_yyy_0_xyyy_1[i] * fti_ab_0 + g_0_yyyz_0_xyyy_0[i] * pb_z + g_0_yyyz_0_xyyy_1[i] * wp_z[i];

        g_0_yyyzz_0_xyyz_0[i] = 2.0 * g_0_yzz_0_xyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xyyz_1[i] * fti_ab_0 + 2.0 * g_0_yyzz_0_xyz_1[i] * fi_abcd_0 +
                                g_0_yyzz_0_xyyz_0[i] * pb_y + g_0_yyzz_0_xyyz_1[i] * wp_y[i];

        g_0_yyyzz_0_xyzz_0[i] = 2.0 * g_0_yzz_0_xyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xyzz_1[i] * fti_ab_0 + g_0_yyzz_0_xzz_1[i] * fi_abcd_0 +
                                g_0_yyzz_0_xyzz_0[i] * pb_y + g_0_yyzz_0_xyzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xzzz_0[i] =
            2.0 * g_0_yzz_0_xzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xzzz_1[i] * fti_ab_0 + g_0_yyzz_0_xzzz_0[i] * pb_y + g_0_yyzz_0_xzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_yyyy_0[i] =
            g_0_yyy_0_yyyy_0[i] * fi_ab_0 - g_0_yyy_0_yyyy_1[i] * fti_ab_0 + g_0_yyyz_0_yyyy_0[i] * pb_z + g_0_yyyz_0_yyyy_1[i] * wp_z[i];

        g_0_yyyzz_0_yyyz_0[i] = 2.0 * g_0_yzz_0_yyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_yyyz_1[i] * fti_ab_0 + 3.0 * g_0_yyzz_0_yyz_1[i] * fi_abcd_0 +
                                g_0_yyzz_0_yyyz_0[i] * pb_y + g_0_yyzz_0_yyyz_1[i] * wp_y[i];

        g_0_yyyzz_0_yyzz_0[i] = 2.0 * g_0_yzz_0_yyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_yyzz_1[i] * fti_ab_0 + 2.0 * g_0_yyzz_0_yzz_1[i] * fi_abcd_0 +
                                g_0_yyzz_0_yyzz_0[i] * pb_y + g_0_yyzz_0_yyzz_1[i] * wp_y[i];

        g_0_yyyzz_0_yzzz_0[i] = 2.0 * g_0_yzz_0_yzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_yzzz_1[i] * fti_ab_0 + g_0_yyzz_0_zzz_1[i] * fi_abcd_0 +
                                g_0_yyzz_0_yzzz_0[i] * pb_y + g_0_yyzz_0_yzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_zzzz_0[i] =
            2.0 * g_0_yzz_0_zzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_zzzz_1[i] * fti_ab_0 + g_0_yyzz_0_zzzz_0[i] * pb_y + g_0_yyzz_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 270-285 components of targeted buffer : SHSG

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

#pragma omp simd aligned(g_0_yyz_0_xxxy_0,       \
                             g_0_yyz_0_xxxy_1,   \
                             g_0_yyz_0_xxyy_0,   \
                             g_0_yyz_0_xxyy_1,   \
                             g_0_yyz_0_xyyy_0,   \
                             g_0_yyz_0_xyyy_1,   \
                             g_0_yyz_0_yyyy_0,   \
                             g_0_yyz_0_yyyy_1,   \
                             g_0_yyzz_0_xxxy_0,  \
                             g_0_yyzz_0_xxxy_1,  \
                             g_0_yyzz_0_xxyy_0,  \
                             g_0_yyzz_0_xxyy_1,  \
                             g_0_yyzz_0_xyyy_0,  \
                             g_0_yyzz_0_xyyy_1,  \
                             g_0_yyzz_0_yyyy_0,  \
                             g_0_yyzz_0_yyyy_1,  \
                             g_0_yyzzz_0_xxxx_0, \
                             g_0_yyzzz_0_xxxy_0, \
                             g_0_yyzzz_0_xxxz_0, \
                             g_0_yyzzz_0_xxyy_0, \
                             g_0_yyzzz_0_xxyz_0, \
                             g_0_yyzzz_0_xxzz_0, \
                             g_0_yyzzz_0_xyyy_0, \
                             g_0_yyzzz_0_xyyz_0, \
                             g_0_yyzzz_0_xyzz_0, \
                             g_0_yyzzz_0_xzzz_0, \
                             g_0_yyzzz_0_yyyy_0, \
                             g_0_yyzzz_0_yyyz_0, \
                             g_0_yyzzz_0_yyzz_0, \
                             g_0_yyzzz_0_yzzz_0, \
                             g_0_yyzzz_0_zzzz_0, \
                             g_0_yzzz_0_xxxx_0,  \
                             g_0_yzzz_0_xxxx_1,  \
                             g_0_yzzz_0_xxxz_0,  \
                             g_0_yzzz_0_xxxz_1,  \
                             g_0_yzzz_0_xxyz_0,  \
                             g_0_yzzz_0_xxyz_1,  \
                             g_0_yzzz_0_xxz_1,   \
                             g_0_yzzz_0_xxzz_0,  \
                             g_0_yzzz_0_xxzz_1,  \
                             g_0_yzzz_0_xyyz_0,  \
                             g_0_yzzz_0_xyyz_1,  \
                             g_0_yzzz_0_xyz_1,   \
                             g_0_yzzz_0_xyzz_0,  \
                             g_0_yzzz_0_xyzz_1,  \
                             g_0_yzzz_0_xzz_1,   \
                             g_0_yzzz_0_xzzz_0,  \
                             g_0_yzzz_0_xzzz_1,  \
                             g_0_yzzz_0_yyyz_0,  \
                             g_0_yzzz_0_yyyz_1,  \
                             g_0_yzzz_0_yyz_1,   \
                             g_0_yzzz_0_yyzz_0,  \
                             g_0_yzzz_0_yyzz_1,  \
                             g_0_yzzz_0_yzz_1,   \
                             g_0_yzzz_0_yzzz_0,  \
                             g_0_yzzz_0_yzzz_1,  \
                             g_0_yzzz_0_zzz_1,   \
                             g_0_yzzz_0_zzzz_0,  \
                             g_0_yzzz_0_zzzz_1,  \
                             g_0_zzz_0_xxxx_0,   \
                             g_0_zzz_0_xxxx_1,   \
                             g_0_zzz_0_xxxz_0,   \
                             g_0_zzz_0_xxxz_1,   \
                             g_0_zzz_0_xxyz_0,   \
                             g_0_zzz_0_xxyz_1,   \
                             g_0_zzz_0_xxzz_0,   \
                             g_0_zzz_0_xxzz_1,   \
                             g_0_zzz_0_xyyz_0,   \
                             g_0_zzz_0_xyyz_1,   \
                             g_0_zzz_0_xyzz_0,   \
                             g_0_zzz_0_xyzz_1,   \
                             g_0_zzz_0_xzzz_0,   \
                             g_0_zzz_0_xzzz_1,   \
                             g_0_zzz_0_yyyz_0,   \
                             g_0_zzz_0_yyyz_1,   \
                             g_0_zzz_0_yyzz_0,   \
                             g_0_zzz_0_yyzz_1,   \
                             g_0_zzz_0_yzzz_0,   \
                             g_0_zzz_0_yzzz_1,   \
                             g_0_zzz_0_zzzz_0,   \
                             g_0_zzz_0_zzzz_1,   \
                             wp_y,               \
                             wp_z,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyzzz_0_xxxx_0[i] =
            g_0_zzz_0_xxxx_0[i] * fi_ab_0 - g_0_zzz_0_xxxx_1[i] * fti_ab_0 + g_0_yzzz_0_xxxx_0[i] * pb_y + g_0_yzzz_0_xxxx_1[i] * wp_y[i];

        g_0_yyzzz_0_xxxy_0[i] =
            2.0 * g_0_yyz_0_xxxy_0[i] * fi_ab_0 - 2.0 * g_0_yyz_0_xxxy_1[i] * fti_ab_0 + g_0_yyzz_0_xxxy_0[i] * pb_z + g_0_yyzz_0_xxxy_1[i] * wp_z[i];

        g_0_yyzzz_0_xxxz_0[i] =
            g_0_zzz_0_xxxz_0[i] * fi_ab_0 - g_0_zzz_0_xxxz_1[i] * fti_ab_0 + g_0_yzzz_0_xxxz_0[i] * pb_y + g_0_yzzz_0_xxxz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxyy_0[i] =
            2.0 * g_0_yyz_0_xxyy_0[i] * fi_ab_0 - 2.0 * g_0_yyz_0_xxyy_1[i] * fti_ab_0 + g_0_yyzz_0_xxyy_0[i] * pb_z + g_0_yyzz_0_xxyy_1[i] * wp_z[i];

        g_0_yyzzz_0_xxyz_0[i] = g_0_zzz_0_xxyz_0[i] * fi_ab_0 - g_0_zzz_0_xxyz_1[i] * fti_ab_0 + g_0_yzzz_0_xxz_1[i] * fi_abcd_0 +
                                g_0_yzzz_0_xxyz_0[i] * pb_y + g_0_yzzz_0_xxyz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxzz_0[i] =
            g_0_zzz_0_xxzz_0[i] * fi_ab_0 - g_0_zzz_0_xxzz_1[i] * fti_ab_0 + g_0_yzzz_0_xxzz_0[i] * pb_y + g_0_yzzz_0_xxzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xyyy_0[i] =
            2.0 * g_0_yyz_0_xyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyz_0_xyyy_1[i] * fti_ab_0 + g_0_yyzz_0_xyyy_0[i] * pb_z + g_0_yyzz_0_xyyy_1[i] * wp_z[i];

        g_0_yyzzz_0_xyyz_0[i] = g_0_zzz_0_xyyz_0[i] * fi_ab_0 - g_0_zzz_0_xyyz_1[i] * fti_ab_0 + 2.0 * g_0_yzzz_0_xyz_1[i] * fi_abcd_0 +
                                g_0_yzzz_0_xyyz_0[i] * pb_y + g_0_yzzz_0_xyyz_1[i] * wp_y[i];

        g_0_yyzzz_0_xyzz_0[i] = g_0_zzz_0_xyzz_0[i] * fi_ab_0 - g_0_zzz_0_xyzz_1[i] * fti_ab_0 + g_0_yzzz_0_xzz_1[i] * fi_abcd_0 +
                                g_0_yzzz_0_xyzz_0[i] * pb_y + g_0_yzzz_0_xyzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xzzz_0[i] =
            g_0_zzz_0_xzzz_0[i] * fi_ab_0 - g_0_zzz_0_xzzz_1[i] * fti_ab_0 + g_0_yzzz_0_xzzz_0[i] * pb_y + g_0_yzzz_0_xzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_yyyy_0[i] =
            2.0 * g_0_yyz_0_yyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyz_0_yyyy_1[i] * fti_ab_0 + g_0_yyzz_0_yyyy_0[i] * pb_z + g_0_yyzz_0_yyyy_1[i] * wp_z[i];

        g_0_yyzzz_0_yyyz_0[i] = g_0_zzz_0_yyyz_0[i] * fi_ab_0 - g_0_zzz_0_yyyz_1[i] * fti_ab_0 + 3.0 * g_0_yzzz_0_yyz_1[i] * fi_abcd_0 +
                                g_0_yzzz_0_yyyz_0[i] * pb_y + g_0_yzzz_0_yyyz_1[i] * wp_y[i];

        g_0_yyzzz_0_yyzz_0[i] = g_0_zzz_0_yyzz_0[i] * fi_ab_0 - g_0_zzz_0_yyzz_1[i] * fti_ab_0 + 2.0 * g_0_yzzz_0_yzz_1[i] * fi_abcd_0 +
                                g_0_yzzz_0_yyzz_0[i] * pb_y + g_0_yzzz_0_yyzz_1[i] * wp_y[i];

        g_0_yyzzz_0_yzzz_0[i] = g_0_zzz_0_yzzz_0[i] * fi_ab_0 - g_0_zzz_0_yzzz_1[i] * fti_ab_0 + g_0_yzzz_0_zzz_1[i] * fi_abcd_0 +
                                g_0_yzzz_0_yzzz_0[i] * pb_y + g_0_yzzz_0_yzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_zzzz_0[i] =
            g_0_zzz_0_zzzz_0[i] * fi_ab_0 - g_0_zzz_0_zzzz_1[i] * fti_ab_0 + g_0_yzzz_0_zzzz_0[i] * pb_y + g_0_yzzz_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 285-300 components of targeted buffer : SHSG

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

#pragma omp simd aligned(g_0_yzzzz_0_xxxx_0,     \
                             g_0_yzzzz_0_xxxy_0, \
                             g_0_yzzzz_0_xxxz_0, \
                             g_0_yzzzz_0_xxyy_0, \
                             g_0_yzzzz_0_xxyz_0, \
                             g_0_yzzzz_0_xxzz_0, \
                             g_0_yzzzz_0_xyyy_0, \
                             g_0_yzzzz_0_xyyz_0, \
                             g_0_yzzzz_0_xyzz_0, \
                             g_0_yzzzz_0_xzzz_0, \
                             g_0_yzzzz_0_yyyy_0, \
                             g_0_yzzzz_0_yyyz_0, \
                             g_0_yzzzz_0_yyzz_0, \
                             g_0_yzzzz_0_yzzz_0, \
                             g_0_yzzzz_0_zzzz_0, \
                             g_0_zzzz_0_xxx_1,   \
                             g_0_zzzz_0_xxxx_0,  \
                             g_0_zzzz_0_xxxx_1,  \
                             g_0_zzzz_0_xxxy_0,  \
                             g_0_zzzz_0_xxxy_1,  \
                             g_0_zzzz_0_xxxz_0,  \
                             g_0_zzzz_0_xxxz_1,  \
                             g_0_zzzz_0_xxy_1,   \
                             g_0_zzzz_0_xxyy_0,  \
                             g_0_zzzz_0_xxyy_1,  \
                             g_0_zzzz_0_xxyz_0,  \
                             g_0_zzzz_0_xxyz_1,  \
                             g_0_zzzz_0_xxz_1,   \
                             g_0_zzzz_0_xxzz_0,  \
                             g_0_zzzz_0_xxzz_1,  \
                             g_0_zzzz_0_xyy_1,   \
                             g_0_zzzz_0_xyyy_0,  \
                             g_0_zzzz_0_xyyy_1,  \
                             g_0_zzzz_0_xyyz_0,  \
                             g_0_zzzz_0_xyyz_1,  \
                             g_0_zzzz_0_xyz_1,   \
                             g_0_zzzz_0_xyzz_0,  \
                             g_0_zzzz_0_xyzz_1,  \
                             g_0_zzzz_0_xzz_1,   \
                             g_0_zzzz_0_xzzz_0,  \
                             g_0_zzzz_0_xzzz_1,  \
                             g_0_zzzz_0_yyy_1,   \
                             g_0_zzzz_0_yyyy_0,  \
                             g_0_zzzz_0_yyyy_1,  \
                             g_0_zzzz_0_yyyz_0,  \
                             g_0_zzzz_0_yyyz_1,  \
                             g_0_zzzz_0_yyz_1,   \
                             g_0_zzzz_0_yyzz_0,  \
                             g_0_zzzz_0_yyzz_1,  \
                             g_0_zzzz_0_yzz_1,   \
                             g_0_zzzz_0_yzzz_0,  \
                             g_0_zzzz_0_yzzz_1,  \
                             g_0_zzzz_0_zzz_1,   \
                             g_0_zzzz_0_zzzz_0,  \
                             g_0_zzzz_0_zzzz_1,  \
                             wp_y,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzzzz_0_xxxx_0[i] = g_0_zzzz_0_xxxx_0[i] * pb_y + g_0_zzzz_0_xxxx_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxy_0[i] = g_0_zzzz_0_xxx_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxy_0[i] * pb_y + g_0_zzzz_0_xxxy_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxz_0[i] = g_0_zzzz_0_xxxz_0[i] * pb_y + g_0_zzzz_0_xxxz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxyy_0[i] = 2.0 * g_0_zzzz_0_xxy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyy_0[i] * pb_y + g_0_zzzz_0_xxyy_1[i] * wp_y[i];

        g_0_yzzzz_0_xxyz_0[i] = g_0_zzzz_0_xxz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyz_0[i] * pb_y + g_0_zzzz_0_xxyz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxzz_0[i] = g_0_zzzz_0_xxzz_0[i] * pb_y + g_0_zzzz_0_xxzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xyyy_0[i] = 3.0 * g_0_zzzz_0_xyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyy_0[i] * pb_y + g_0_zzzz_0_xyyy_1[i] * wp_y[i];

        g_0_yzzzz_0_xyyz_0[i] = 2.0 * g_0_zzzz_0_xyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyz_0[i] * pb_y + g_0_zzzz_0_xyyz_1[i] * wp_y[i];

        g_0_yzzzz_0_xyzz_0[i] = g_0_zzzz_0_xzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyzz_0[i] * pb_y + g_0_zzzz_0_xyzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xzzz_0[i] = g_0_zzzz_0_xzzz_0[i] * pb_y + g_0_zzzz_0_xzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_yyyy_0[i] = 4.0 * g_0_zzzz_0_yyy_1[i] * fi_abcd_0 + g_0_zzzz_0_yyyy_0[i] * pb_y + g_0_zzzz_0_yyyy_1[i] * wp_y[i];

        g_0_yzzzz_0_yyyz_0[i] = 3.0 * g_0_zzzz_0_yyz_1[i] * fi_abcd_0 + g_0_zzzz_0_yyyz_0[i] * pb_y + g_0_zzzz_0_yyyz_1[i] * wp_y[i];

        g_0_yzzzz_0_yyzz_0[i] = 2.0 * g_0_zzzz_0_yzz_1[i] * fi_abcd_0 + g_0_zzzz_0_yyzz_0[i] * pb_y + g_0_zzzz_0_yyzz_1[i] * wp_y[i];

        g_0_yzzzz_0_yzzz_0[i] = g_0_zzzz_0_zzz_1[i] * fi_abcd_0 + g_0_zzzz_0_yzzz_0[i] * pb_y + g_0_zzzz_0_yzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_zzzz_0[i] = g_0_zzzz_0_zzzz_0[i] * pb_y + g_0_zzzz_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 300-315 components of targeted buffer : SHSG

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

#pragma omp simd aligned(g_0_zzz_0_xxxx_0,       \
                             g_0_zzz_0_xxxx_1,   \
                             g_0_zzz_0_xxxy_0,   \
                             g_0_zzz_0_xxxy_1,   \
                             g_0_zzz_0_xxxz_0,   \
                             g_0_zzz_0_xxxz_1,   \
                             g_0_zzz_0_xxyy_0,   \
                             g_0_zzz_0_xxyy_1,   \
                             g_0_zzz_0_xxyz_0,   \
                             g_0_zzz_0_xxyz_1,   \
                             g_0_zzz_0_xxzz_0,   \
                             g_0_zzz_0_xxzz_1,   \
                             g_0_zzz_0_xyyy_0,   \
                             g_0_zzz_0_xyyy_1,   \
                             g_0_zzz_0_xyyz_0,   \
                             g_0_zzz_0_xyyz_1,   \
                             g_0_zzz_0_xyzz_0,   \
                             g_0_zzz_0_xyzz_1,   \
                             g_0_zzz_0_xzzz_0,   \
                             g_0_zzz_0_xzzz_1,   \
                             g_0_zzz_0_yyyy_0,   \
                             g_0_zzz_0_yyyy_1,   \
                             g_0_zzz_0_yyyz_0,   \
                             g_0_zzz_0_yyyz_1,   \
                             g_0_zzz_0_yyzz_0,   \
                             g_0_zzz_0_yyzz_1,   \
                             g_0_zzz_0_yzzz_0,   \
                             g_0_zzz_0_yzzz_1,   \
                             g_0_zzz_0_zzzz_0,   \
                             g_0_zzz_0_zzzz_1,   \
                             g_0_zzzz_0_xxx_1,   \
                             g_0_zzzz_0_xxxx_0,  \
                             g_0_zzzz_0_xxxx_1,  \
                             g_0_zzzz_0_xxxy_0,  \
                             g_0_zzzz_0_xxxy_1,  \
                             g_0_zzzz_0_xxxz_0,  \
                             g_0_zzzz_0_xxxz_1,  \
                             g_0_zzzz_0_xxy_1,   \
                             g_0_zzzz_0_xxyy_0,  \
                             g_0_zzzz_0_xxyy_1,  \
                             g_0_zzzz_0_xxyz_0,  \
                             g_0_zzzz_0_xxyz_1,  \
                             g_0_zzzz_0_xxz_1,   \
                             g_0_zzzz_0_xxzz_0,  \
                             g_0_zzzz_0_xxzz_1,  \
                             g_0_zzzz_0_xyy_1,   \
                             g_0_zzzz_0_xyyy_0,  \
                             g_0_zzzz_0_xyyy_1,  \
                             g_0_zzzz_0_xyyz_0,  \
                             g_0_zzzz_0_xyyz_1,  \
                             g_0_zzzz_0_xyz_1,   \
                             g_0_zzzz_0_xyzz_0,  \
                             g_0_zzzz_0_xyzz_1,  \
                             g_0_zzzz_0_xzz_1,   \
                             g_0_zzzz_0_xzzz_0,  \
                             g_0_zzzz_0_xzzz_1,  \
                             g_0_zzzz_0_yyy_1,   \
                             g_0_zzzz_0_yyyy_0,  \
                             g_0_zzzz_0_yyyy_1,  \
                             g_0_zzzz_0_yyyz_0,  \
                             g_0_zzzz_0_yyyz_1,  \
                             g_0_zzzz_0_yyz_1,   \
                             g_0_zzzz_0_yyzz_0,  \
                             g_0_zzzz_0_yyzz_1,  \
                             g_0_zzzz_0_yzz_1,   \
                             g_0_zzzz_0_yzzz_0,  \
                             g_0_zzzz_0_yzzz_1,  \
                             g_0_zzzz_0_zzz_1,   \
                             g_0_zzzz_0_zzzz_0,  \
                             g_0_zzzz_0_zzzz_1,  \
                             g_0_zzzzz_0_xxxx_0, \
                             g_0_zzzzz_0_xxxy_0, \
                             g_0_zzzzz_0_xxxz_0, \
                             g_0_zzzzz_0_xxyy_0, \
                             g_0_zzzzz_0_xxyz_0, \
                             g_0_zzzzz_0_xxzz_0, \
                             g_0_zzzzz_0_xyyy_0, \
                             g_0_zzzzz_0_xyyz_0, \
                             g_0_zzzzz_0_xyzz_0, \
                             g_0_zzzzz_0_xzzz_0, \
                             g_0_zzzzz_0_yyyy_0, \
                             g_0_zzzzz_0_yyyz_0, \
                             g_0_zzzzz_0_yyzz_0, \
                             g_0_zzzzz_0_yzzz_0, \
                             g_0_zzzzz_0_zzzz_0, \
                             wp_z,               \
                             c_exps,             \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzzzz_0_xxxx_0[i] =
            4.0 * g_0_zzz_0_xxxx_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxx_1[i] * fti_ab_0 + g_0_zzzz_0_xxxx_0[i] * pb_z + g_0_zzzz_0_xxxx_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxy_0[i] =
            4.0 * g_0_zzz_0_xxxy_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxy_1[i] * fti_ab_0 + g_0_zzzz_0_xxxy_0[i] * pb_z + g_0_zzzz_0_xxxy_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxz_0[i] = 4.0 * g_0_zzz_0_xxxz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxz_1[i] * fti_ab_0 + g_0_zzzz_0_xxx_1[i] * fi_abcd_0 +
                                g_0_zzzz_0_xxxz_0[i] * pb_z + g_0_zzzz_0_xxxz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxyy_0[i] =
            4.0 * g_0_zzz_0_xxyy_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxyy_1[i] * fti_ab_0 + g_0_zzzz_0_xxyy_0[i] * pb_z + g_0_zzzz_0_xxyy_1[i] * wp_z[i];

        g_0_zzzzz_0_xxyz_0[i] = 4.0 * g_0_zzz_0_xxyz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxyz_1[i] * fti_ab_0 + g_0_zzzz_0_xxy_1[i] * fi_abcd_0 +
                                g_0_zzzz_0_xxyz_0[i] * pb_z + g_0_zzzz_0_xxyz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxzz_0[i] = 4.0 * g_0_zzz_0_xxzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzz_0_xxz_1[i] * fi_abcd_0 +
                                g_0_zzzz_0_xxzz_0[i] * pb_z + g_0_zzzz_0_xxzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xyyy_0[i] =
            4.0 * g_0_zzz_0_xyyy_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xyyy_1[i] * fti_ab_0 + g_0_zzzz_0_xyyy_0[i] * pb_z + g_0_zzzz_0_xyyy_1[i] * wp_z[i];

        g_0_zzzzz_0_xyyz_0[i] = 4.0 * g_0_zzz_0_xyyz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xyyz_1[i] * fti_ab_0 + g_0_zzzz_0_xyy_1[i] * fi_abcd_0 +
                                g_0_zzzz_0_xyyz_0[i] * pb_z + g_0_zzzz_0_xyyz_1[i] * wp_z[i];

        g_0_zzzzz_0_xyzz_0[i] = 4.0 * g_0_zzz_0_xyzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzz_0_xyz_1[i] * fi_abcd_0 +
                                g_0_zzzz_0_xyzz_0[i] * pb_z + g_0_zzzz_0_xyzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xzzz_0[i] = 4.0 * g_0_zzz_0_xzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzz_0_xzz_1[i] * fi_abcd_0 +
                                g_0_zzzz_0_xzzz_0[i] * pb_z + g_0_zzzz_0_xzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_yyyy_0[i] =
            4.0 * g_0_zzz_0_yyyy_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_yyyy_1[i] * fti_ab_0 + g_0_zzzz_0_yyyy_0[i] * pb_z + g_0_zzzz_0_yyyy_1[i] * wp_z[i];

        g_0_zzzzz_0_yyyz_0[i] = 4.0 * g_0_zzz_0_yyyz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_yyyz_1[i] * fti_ab_0 + g_0_zzzz_0_yyy_1[i] * fi_abcd_0 +
                                g_0_zzzz_0_yyyz_0[i] * pb_z + g_0_zzzz_0_yyyz_1[i] * wp_z[i];

        g_0_zzzzz_0_yyzz_0[i] = 4.0 * g_0_zzz_0_yyzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_yyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzz_0_yyz_1[i] * fi_abcd_0 +
                                g_0_zzzz_0_yyzz_0[i] * pb_z + g_0_zzzz_0_yyzz_1[i] * wp_z[i];

        g_0_zzzzz_0_yzzz_0[i] = 4.0 * g_0_zzz_0_yzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_yzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzz_0_yzz_1[i] * fi_abcd_0 +
                                g_0_zzzz_0_yzzz_0[i] * pb_z + g_0_zzzz_0_yzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_zzzz_0[i] = 4.0 * g_0_zzz_0_zzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_zzzz_1[i] * fti_ab_0 + 4.0 * g_0_zzzz_0_zzz_1[i] * fi_abcd_0 +
                                g_0_zzzz_0_zzzz_0[i] * pb_z + g_0_zzzz_0_zzzz_1[i] * wp_z[i];
    }
}

}  // namespace erirec

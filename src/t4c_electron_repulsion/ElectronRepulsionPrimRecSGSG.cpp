#include "ElectronRepulsionPrimRecSGSG.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_sgsg(CSimdArray<double>& pbuffer,
                                  const size_t idx_eri_0_sgsg,
                                  size_t idx_eri_0_sdsg,
                                  size_t idx_eri_1_sdsg,
                                  size_t idx_eri_1_sfsf,
                                  size_t idx_eri_0_sfsg,
                                  size_t idx_eri_1_sfsg,
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

    /// Set up components of auxilary buffer : SDSG

    auto g_0_xx_0_xxxx_0 = pbuffer.data(idx_eri_0_sdsg);

    auto g_0_xx_0_xxxy_0 = pbuffer.data(idx_eri_0_sdsg + 1);

    auto g_0_xx_0_xxxz_0 = pbuffer.data(idx_eri_0_sdsg + 2);

    auto g_0_xx_0_xxyy_0 = pbuffer.data(idx_eri_0_sdsg + 3);

    auto g_0_xx_0_xxyz_0 = pbuffer.data(idx_eri_0_sdsg + 4);

    auto g_0_xx_0_xxzz_0 = pbuffer.data(idx_eri_0_sdsg + 5);

    auto g_0_xx_0_xyyy_0 = pbuffer.data(idx_eri_0_sdsg + 6);

    auto g_0_xx_0_xyyz_0 = pbuffer.data(idx_eri_0_sdsg + 7);

    auto g_0_xx_0_xyzz_0 = pbuffer.data(idx_eri_0_sdsg + 8);

    auto g_0_xx_0_xzzz_0 = pbuffer.data(idx_eri_0_sdsg + 9);

    auto g_0_xx_0_yyyy_0 = pbuffer.data(idx_eri_0_sdsg + 10);

    auto g_0_xx_0_yyyz_0 = pbuffer.data(idx_eri_0_sdsg + 11);

    auto g_0_xx_0_yyzz_0 = pbuffer.data(idx_eri_0_sdsg + 12);

    auto g_0_xx_0_yzzz_0 = pbuffer.data(idx_eri_0_sdsg + 13);

    auto g_0_xx_0_zzzz_0 = pbuffer.data(idx_eri_0_sdsg + 14);

    auto g_0_yy_0_xxxx_0 = pbuffer.data(idx_eri_0_sdsg + 45);

    auto g_0_yy_0_xxxy_0 = pbuffer.data(idx_eri_0_sdsg + 46);

    auto g_0_yy_0_xxxz_0 = pbuffer.data(idx_eri_0_sdsg + 47);

    auto g_0_yy_0_xxyy_0 = pbuffer.data(idx_eri_0_sdsg + 48);

    auto g_0_yy_0_xxyz_0 = pbuffer.data(idx_eri_0_sdsg + 49);

    auto g_0_yy_0_xxzz_0 = pbuffer.data(idx_eri_0_sdsg + 50);

    auto g_0_yy_0_xyyy_0 = pbuffer.data(idx_eri_0_sdsg + 51);

    auto g_0_yy_0_xyyz_0 = pbuffer.data(idx_eri_0_sdsg + 52);

    auto g_0_yy_0_xyzz_0 = pbuffer.data(idx_eri_0_sdsg + 53);

    auto g_0_yy_0_xzzz_0 = pbuffer.data(idx_eri_0_sdsg + 54);

    auto g_0_yy_0_yyyy_0 = pbuffer.data(idx_eri_0_sdsg + 55);

    auto g_0_yy_0_yyyz_0 = pbuffer.data(idx_eri_0_sdsg + 56);

    auto g_0_yy_0_yyzz_0 = pbuffer.data(idx_eri_0_sdsg + 57);

    auto g_0_yy_0_yzzz_0 = pbuffer.data(idx_eri_0_sdsg + 58);

    auto g_0_yy_0_zzzz_0 = pbuffer.data(idx_eri_0_sdsg + 59);

    auto g_0_zz_0_xxxx_0 = pbuffer.data(idx_eri_0_sdsg + 75);

    auto g_0_zz_0_xxxy_0 = pbuffer.data(idx_eri_0_sdsg + 76);

    auto g_0_zz_0_xxxz_0 = pbuffer.data(idx_eri_0_sdsg + 77);

    auto g_0_zz_0_xxyy_0 = pbuffer.data(idx_eri_0_sdsg + 78);

    auto g_0_zz_0_xxyz_0 = pbuffer.data(idx_eri_0_sdsg + 79);

    auto g_0_zz_0_xxzz_0 = pbuffer.data(idx_eri_0_sdsg + 80);

    auto g_0_zz_0_xyyy_0 = pbuffer.data(idx_eri_0_sdsg + 81);

    auto g_0_zz_0_xyyz_0 = pbuffer.data(idx_eri_0_sdsg + 82);

    auto g_0_zz_0_xyzz_0 = pbuffer.data(idx_eri_0_sdsg + 83);

    auto g_0_zz_0_xzzz_0 = pbuffer.data(idx_eri_0_sdsg + 84);

    auto g_0_zz_0_yyyy_0 = pbuffer.data(idx_eri_0_sdsg + 85);

    auto g_0_zz_0_yyyz_0 = pbuffer.data(idx_eri_0_sdsg + 86);

    auto g_0_zz_0_yyzz_0 = pbuffer.data(idx_eri_0_sdsg + 87);

    auto g_0_zz_0_yzzz_0 = pbuffer.data(idx_eri_0_sdsg + 88);

    auto g_0_zz_0_zzzz_0 = pbuffer.data(idx_eri_0_sdsg + 89);

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

    /// Set up components of auxilary buffer : SFSF

    auto g_0_xxx_0_xxx_1 = pbuffer.data(idx_eri_1_sfsf);

    auto g_0_xxx_0_xxy_1 = pbuffer.data(idx_eri_1_sfsf + 1);

    auto g_0_xxx_0_xxz_1 = pbuffer.data(idx_eri_1_sfsf + 2);

    auto g_0_xxx_0_xyy_1 = pbuffer.data(idx_eri_1_sfsf + 3);

    auto g_0_xxx_0_xyz_1 = pbuffer.data(idx_eri_1_sfsf + 4);

    auto g_0_xxx_0_xzz_1 = pbuffer.data(idx_eri_1_sfsf + 5);

    auto g_0_xxx_0_yyy_1 = pbuffer.data(idx_eri_1_sfsf + 6);

    auto g_0_xxx_0_yyz_1 = pbuffer.data(idx_eri_1_sfsf + 7);

    auto g_0_xxx_0_yzz_1 = pbuffer.data(idx_eri_1_sfsf + 8);

    auto g_0_xxx_0_zzz_1 = pbuffer.data(idx_eri_1_sfsf + 9);

    auto g_0_xxz_0_xxz_1 = pbuffer.data(idx_eri_1_sfsf + 22);

    auto g_0_xxz_0_xyz_1 = pbuffer.data(idx_eri_1_sfsf + 24);

    auto g_0_xxz_0_xzz_1 = pbuffer.data(idx_eri_1_sfsf + 25);

    auto g_0_xxz_0_yyz_1 = pbuffer.data(idx_eri_1_sfsf + 27);

    auto g_0_xxz_0_yzz_1 = pbuffer.data(idx_eri_1_sfsf + 28);

    auto g_0_xxz_0_zzz_1 = pbuffer.data(idx_eri_1_sfsf + 29);

    auto g_0_xyy_0_xxy_1 = pbuffer.data(idx_eri_1_sfsf + 31);

    auto g_0_xyy_0_xyy_1 = pbuffer.data(idx_eri_1_sfsf + 33);

    auto g_0_xyy_0_xyz_1 = pbuffer.data(idx_eri_1_sfsf + 34);

    auto g_0_xyy_0_yyy_1 = pbuffer.data(idx_eri_1_sfsf + 36);

    auto g_0_xyy_0_yyz_1 = pbuffer.data(idx_eri_1_sfsf + 37);

    auto g_0_xyy_0_yzz_1 = pbuffer.data(idx_eri_1_sfsf + 38);

    auto g_0_xzz_0_xxz_1 = pbuffer.data(idx_eri_1_sfsf + 52);

    auto g_0_xzz_0_xyz_1 = pbuffer.data(idx_eri_1_sfsf + 54);

    auto g_0_xzz_0_xzz_1 = pbuffer.data(idx_eri_1_sfsf + 55);

    auto g_0_xzz_0_yyz_1 = pbuffer.data(idx_eri_1_sfsf + 57);

    auto g_0_xzz_0_yzz_1 = pbuffer.data(idx_eri_1_sfsf + 58);

    auto g_0_xzz_0_zzz_1 = pbuffer.data(idx_eri_1_sfsf + 59);

    auto g_0_yyy_0_xxx_1 = pbuffer.data(idx_eri_1_sfsf + 60);

    auto g_0_yyy_0_xxy_1 = pbuffer.data(idx_eri_1_sfsf + 61);

    auto g_0_yyy_0_xxz_1 = pbuffer.data(idx_eri_1_sfsf + 62);

    auto g_0_yyy_0_xyy_1 = pbuffer.data(idx_eri_1_sfsf + 63);

    auto g_0_yyy_0_xyz_1 = pbuffer.data(idx_eri_1_sfsf + 64);

    auto g_0_yyy_0_xzz_1 = pbuffer.data(idx_eri_1_sfsf + 65);

    auto g_0_yyy_0_yyy_1 = pbuffer.data(idx_eri_1_sfsf + 66);

    auto g_0_yyy_0_yyz_1 = pbuffer.data(idx_eri_1_sfsf + 67);

    auto g_0_yyy_0_yzz_1 = pbuffer.data(idx_eri_1_sfsf + 68);

    auto g_0_yyy_0_zzz_1 = pbuffer.data(idx_eri_1_sfsf + 69);

    auto g_0_yyz_0_xxz_1 = pbuffer.data(idx_eri_1_sfsf + 72);

    auto g_0_yyz_0_xyz_1 = pbuffer.data(idx_eri_1_sfsf + 74);

    auto g_0_yyz_0_xzz_1 = pbuffer.data(idx_eri_1_sfsf + 75);

    auto g_0_yyz_0_yyz_1 = pbuffer.data(idx_eri_1_sfsf + 77);

    auto g_0_yyz_0_yzz_1 = pbuffer.data(idx_eri_1_sfsf + 78);

    auto g_0_yyz_0_zzz_1 = pbuffer.data(idx_eri_1_sfsf + 79);

    auto g_0_yzz_0_xxy_1 = pbuffer.data(idx_eri_1_sfsf + 81);

    auto g_0_yzz_0_xxz_1 = pbuffer.data(idx_eri_1_sfsf + 82);

    auto g_0_yzz_0_xyy_1 = pbuffer.data(idx_eri_1_sfsf + 83);

    auto g_0_yzz_0_xyz_1 = pbuffer.data(idx_eri_1_sfsf + 84);

    auto g_0_yzz_0_xzz_1 = pbuffer.data(idx_eri_1_sfsf + 85);

    auto g_0_yzz_0_yyy_1 = pbuffer.data(idx_eri_1_sfsf + 86);

    auto g_0_yzz_0_yyz_1 = pbuffer.data(idx_eri_1_sfsf + 87);

    auto g_0_yzz_0_yzz_1 = pbuffer.data(idx_eri_1_sfsf + 88);

    auto g_0_yzz_0_zzz_1 = pbuffer.data(idx_eri_1_sfsf + 89);

    auto g_0_zzz_0_xxx_1 = pbuffer.data(idx_eri_1_sfsf + 90);

    auto g_0_zzz_0_xxy_1 = pbuffer.data(idx_eri_1_sfsf + 91);

    auto g_0_zzz_0_xxz_1 = pbuffer.data(idx_eri_1_sfsf + 92);

    auto g_0_zzz_0_xyy_1 = pbuffer.data(idx_eri_1_sfsf + 93);

    auto g_0_zzz_0_xyz_1 = pbuffer.data(idx_eri_1_sfsf + 94);

    auto g_0_zzz_0_xzz_1 = pbuffer.data(idx_eri_1_sfsf + 95);

    auto g_0_zzz_0_yyy_1 = pbuffer.data(idx_eri_1_sfsf + 96);

    auto g_0_zzz_0_yyz_1 = pbuffer.data(idx_eri_1_sfsf + 97);

    auto g_0_zzz_0_yzz_1 = pbuffer.data(idx_eri_1_sfsf + 98);

    auto g_0_zzz_0_zzz_1 = pbuffer.data(idx_eri_1_sfsf + 99);

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

    auto g_0_xxy_0_xxxy_0 = pbuffer.data(idx_eri_0_sfsg + 16);

    auto g_0_xxy_0_xxxz_0 = pbuffer.data(idx_eri_0_sfsg + 17);

    auto g_0_xxy_0_xxyy_0 = pbuffer.data(idx_eri_0_sfsg + 18);

    auto g_0_xxy_0_xxzz_0 = pbuffer.data(idx_eri_0_sfsg + 20);

    auto g_0_xxy_0_xyyy_0 = pbuffer.data(idx_eri_0_sfsg + 21);

    auto g_0_xxy_0_xzzz_0 = pbuffer.data(idx_eri_0_sfsg + 24);

    auto g_0_xxy_0_yyyy_0 = pbuffer.data(idx_eri_0_sfsg + 25);

    auto g_0_xxz_0_xxxx_0 = pbuffer.data(idx_eri_0_sfsg + 30);

    auto g_0_xxz_0_xxxy_0 = pbuffer.data(idx_eri_0_sfsg + 31);

    auto g_0_xxz_0_xxxz_0 = pbuffer.data(idx_eri_0_sfsg + 32);

    auto g_0_xxz_0_xxyy_0 = pbuffer.data(idx_eri_0_sfsg + 33);

    auto g_0_xxz_0_xxyz_0 = pbuffer.data(idx_eri_0_sfsg + 34);

    auto g_0_xxz_0_xxzz_0 = pbuffer.data(idx_eri_0_sfsg + 35);

    auto g_0_xxz_0_xyyy_0 = pbuffer.data(idx_eri_0_sfsg + 36);

    auto g_0_xxz_0_xyyz_0 = pbuffer.data(idx_eri_0_sfsg + 37);

    auto g_0_xxz_0_xyzz_0 = pbuffer.data(idx_eri_0_sfsg + 38);

    auto g_0_xxz_0_xzzz_0 = pbuffer.data(idx_eri_0_sfsg + 39);

    auto g_0_xxz_0_yyyz_0 = pbuffer.data(idx_eri_0_sfsg + 41);

    auto g_0_xxz_0_yyzz_0 = pbuffer.data(idx_eri_0_sfsg + 42);

    auto g_0_xxz_0_yzzz_0 = pbuffer.data(idx_eri_0_sfsg + 43);

    auto g_0_xxz_0_zzzz_0 = pbuffer.data(idx_eri_0_sfsg + 44);

    auto g_0_xyy_0_xxxx_0 = pbuffer.data(idx_eri_0_sfsg + 45);

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

    auto g_0_xzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sfsg + 75);

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

    auto g_0_yyz_0_xxxz_0 = pbuffer.data(idx_eri_0_sfsg + 107);

    auto g_0_yyz_0_xxyy_0 = pbuffer.data(idx_eri_0_sfsg + 108);

    auto g_0_yyz_0_xxyz_0 = pbuffer.data(idx_eri_0_sfsg + 109);

    auto g_0_yyz_0_xxzz_0 = pbuffer.data(idx_eri_0_sfsg + 110);

    auto g_0_yyz_0_xyyy_0 = pbuffer.data(idx_eri_0_sfsg + 111);

    auto g_0_yyz_0_xyyz_0 = pbuffer.data(idx_eri_0_sfsg + 112);

    auto g_0_yyz_0_xyzz_0 = pbuffer.data(idx_eri_0_sfsg + 113);

    auto g_0_yyz_0_xzzz_0 = pbuffer.data(idx_eri_0_sfsg + 114);

    auto g_0_yyz_0_yyyy_0 = pbuffer.data(idx_eri_0_sfsg + 115);

    auto g_0_yyz_0_yyyz_0 = pbuffer.data(idx_eri_0_sfsg + 116);

    auto g_0_yyz_0_yyzz_0 = pbuffer.data(idx_eri_0_sfsg + 117);

    auto g_0_yyz_0_yzzz_0 = pbuffer.data(idx_eri_0_sfsg + 118);

    auto g_0_yyz_0_zzzz_0 = pbuffer.data(idx_eri_0_sfsg + 119);

    auto g_0_yzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sfsg + 120);

    auto g_0_yzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sfsg + 121);

    auto g_0_yzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sfsg + 122);

    auto g_0_yzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sfsg + 123);

    auto g_0_yzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sfsg + 124);

    auto g_0_yzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sfsg + 125);

    auto g_0_yzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sfsg + 126);

    auto g_0_yzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sfsg + 127);

    auto g_0_yzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sfsg + 128);

    auto g_0_yzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sfsg + 129);

    auto g_0_yzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sfsg + 130);

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

    auto g_0_xxy_0_xxxy_1 = pbuffer.data(idx_eri_1_sfsg + 16);

    auto g_0_xxy_0_xxxz_1 = pbuffer.data(idx_eri_1_sfsg + 17);

    auto g_0_xxy_0_xxyy_1 = pbuffer.data(idx_eri_1_sfsg + 18);

    auto g_0_xxy_0_xxzz_1 = pbuffer.data(idx_eri_1_sfsg + 20);

    auto g_0_xxy_0_xyyy_1 = pbuffer.data(idx_eri_1_sfsg + 21);

    auto g_0_xxy_0_xzzz_1 = pbuffer.data(idx_eri_1_sfsg + 24);

    auto g_0_xxy_0_yyyy_1 = pbuffer.data(idx_eri_1_sfsg + 25);

    auto g_0_xxz_0_xxxx_1 = pbuffer.data(idx_eri_1_sfsg + 30);

    auto g_0_xxz_0_xxxy_1 = pbuffer.data(idx_eri_1_sfsg + 31);

    auto g_0_xxz_0_xxxz_1 = pbuffer.data(idx_eri_1_sfsg + 32);

    auto g_0_xxz_0_xxyy_1 = pbuffer.data(idx_eri_1_sfsg + 33);

    auto g_0_xxz_0_xxyz_1 = pbuffer.data(idx_eri_1_sfsg + 34);

    auto g_0_xxz_0_xxzz_1 = pbuffer.data(idx_eri_1_sfsg + 35);

    auto g_0_xxz_0_xyyy_1 = pbuffer.data(idx_eri_1_sfsg + 36);

    auto g_0_xxz_0_xyyz_1 = pbuffer.data(idx_eri_1_sfsg + 37);

    auto g_0_xxz_0_xyzz_1 = pbuffer.data(idx_eri_1_sfsg + 38);

    auto g_0_xxz_0_xzzz_1 = pbuffer.data(idx_eri_1_sfsg + 39);

    auto g_0_xxz_0_yyyz_1 = pbuffer.data(idx_eri_1_sfsg + 41);

    auto g_0_xxz_0_yyzz_1 = pbuffer.data(idx_eri_1_sfsg + 42);

    auto g_0_xxz_0_yzzz_1 = pbuffer.data(idx_eri_1_sfsg + 43);

    auto g_0_xxz_0_zzzz_1 = pbuffer.data(idx_eri_1_sfsg + 44);

    auto g_0_xyy_0_xxxx_1 = pbuffer.data(idx_eri_1_sfsg + 45);

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

    auto g_0_xzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sfsg + 75);

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

    auto g_0_yyz_0_xxxz_1 = pbuffer.data(idx_eri_1_sfsg + 107);

    auto g_0_yyz_0_xxyy_1 = pbuffer.data(idx_eri_1_sfsg + 108);

    auto g_0_yyz_0_xxyz_1 = pbuffer.data(idx_eri_1_sfsg + 109);

    auto g_0_yyz_0_xxzz_1 = pbuffer.data(idx_eri_1_sfsg + 110);

    auto g_0_yyz_0_xyyy_1 = pbuffer.data(idx_eri_1_sfsg + 111);

    auto g_0_yyz_0_xyyz_1 = pbuffer.data(idx_eri_1_sfsg + 112);

    auto g_0_yyz_0_xyzz_1 = pbuffer.data(idx_eri_1_sfsg + 113);

    auto g_0_yyz_0_xzzz_1 = pbuffer.data(idx_eri_1_sfsg + 114);

    auto g_0_yyz_0_yyyy_1 = pbuffer.data(idx_eri_1_sfsg + 115);

    auto g_0_yyz_0_yyyz_1 = pbuffer.data(idx_eri_1_sfsg + 116);

    auto g_0_yyz_0_yyzz_1 = pbuffer.data(idx_eri_1_sfsg + 117);

    auto g_0_yyz_0_yzzz_1 = pbuffer.data(idx_eri_1_sfsg + 118);

    auto g_0_yyz_0_zzzz_1 = pbuffer.data(idx_eri_1_sfsg + 119);

    auto g_0_yzz_0_xxxx_1 = pbuffer.data(idx_eri_1_sfsg + 120);

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

    /// Set up 0-15 components of targeted buffer : SGSG

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

    #pragma omp simd aligned(g_0_xx_0_xxxx_0, g_0_xx_0_xxxx_1, g_0_xx_0_xxxy_0, g_0_xx_0_xxxy_1, g_0_xx_0_xxxz_0, g_0_xx_0_xxxz_1, g_0_xx_0_xxyy_0, g_0_xx_0_xxyy_1, g_0_xx_0_xxyz_0, g_0_xx_0_xxyz_1, g_0_xx_0_xxzz_0, g_0_xx_0_xxzz_1, g_0_xx_0_xyyy_0, g_0_xx_0_xyyy_1, g_0_xx_0_xyyz_0, g_0_xx_0_xyyz_1, g_0_xx_0_xyzz_0, g_0_xx_0_xyzz_1, g_0_xx_0_xzzz_0, g_0_xx_0_xzzz_1, g_0_xx_0_yyyy_0, g_0_xx_0_yyyy_1, g_0_xx_0_yyyz_0, g_0_xx_0_yyyz_1, g_0_xx_0_yyzz_0, g_0_xx_0_yyzz_1, g_0_xx_0_yzzz_0, g_0_xx_0_yzzz_1, g_0_xx_0_zzzz_0, g_0_xx_0_zzzz_1, g_0_xxx_0_xxx_1, g_0_xxx_0_xxxx_0, g_0_xxx_0_xxxx_1, g_0_xxx_0_xxxy_0, g_0_xxx_0_xxxy_1, g_0_xxx_0_xxxz_0, g_0_xxx_0_xxxz_1, g_0_xxx_0_xxy_1, g_0_xxx_0_xxyy_0, g_0_xxx_0_xxyy_1, g_0_xxx_0_xxyz_0, g_0_xxx_0_xxyz_1, g_0_xxx_0_xxz_1, g_0_xxx_0_xxzz_0, g_0_xxx_0_xxzz_1, g_0_xxx_0_xyy_1, g_0_xxx_0_xyyy_0, g_0_xxx_0_xyyy_1, g_0_xxx_0_xyyz_0, g_0_xxx_0_xyyz_1, g_0_xxx_0_xyz_1, g_0_xxx_0_xyzz_0, g_0_xxx_0_xyzz_1, g_0_xxx_0_xzz_1, g_0_xxx_0_xzzz_0, g_0_xxx_0_xzzz_1, g_0_xxx_0_yyy_1, g_0_xxx_0_yyyy_0, g_0_xxx_0_yyyy_1, g_0_xxx_0_yyyz_0, g_0_xxx_0_yyyz_1, g_0_xxx_0_yyz_1, g_0_xxx_0_yyzz_0, g_0_xxx_0_yyzz_1, g_0_xxx_0_yzz_1, g_0_xxx_0_yzzz_0, g_0_xxx_0_yzzz_1, g_0_xxx_0_zzz_1, g_0_xxx_0_zzzz_0, g_0_xxx_0_zzzz_1, g_0_xxxx_0_xxxx_0, g_0_xxxx_0_xxxy_0, g_0_xxxx_0_xxxz_0, g_0_xxxx_0_xxyy_0, g_0_xxxx_0_xxyz_0, g_0_xxxx_0_xxzz_0, g_0_xxxx_0_xyyy_0, g_0_xxxx_0_xyyz_0, g_0_xxxx_0_xyzz_0, g_0_xxxx_0_xzzz_0, g_0_xxxx_0_yyyy_0, g_0_xxxx_0_yyyz_0, g_0_xxxx_0_yyzz_0, g_0_xxxx_0_yzzz_0, g_0_xxxx_0_zzzz_0, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxx_0_xxxx_0[i] = 3.0 * g_0_xx_0_xxxx_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxx_1[i] * fti_ab_0 + 4.0 * g_0_xxx_0_xxx_1[i] * fi_abcd_0 + g_0_xxx_0_xxxx_0[i] * pb_x + g_0_xxx_0_xxxx_1[i] * wp_x[i];

        g_0_xxxx_0_xxxy_0[i] = 3.0 * g_0_xx_0_xxxy_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxy_1[i] * fti_ab_0 + 3.0 * g_0_xxx_0_xxy_1[i] * fi_abcd_0 + g_0_xxx_0_xxxy_0[i] * pb_x + g_0_xxx_0_xxxy_1[i] * wp_x[i];

        g_0_xxxx_0_xxxz_0[i] = 3.0 * g_0_xx_0_xxxz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxz_1[i] * fti_ab_0 + 3.0 * g_0_xxx_0_xxz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxz_0[i] * pb_x + g_0_xxx_0_xxxz_1[i] * wp_x[i];

        g_0_xxxx_0_xxyy_0[i] = 3.0 * g_0_xx_0_xxyy_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxyy_1[i] * fti_ab_0 + 2.0 * g_0_xxx_0_xyy_1[i] * fi_abcd_0 + g_0_xxx_0_xxyy_0[i] * pb_x + g_0_xxx_0_xxyy_1[i] * wp_x[i];

        g_0_xxxx_0_xxyz_0[i] = 3.0 * g_0_xx_0_xxyz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxyz_1[i] * fti_ab_0 + 2.0 * g_0_xxx_0_xyz_1[i] * fi_abcd_0 + g_0_xxx_0_xxyz_0[i] * pb_x + g_0_xxx_0_xxyz_1[i] * wp_x[i];

        g_0_xxxx_0_xxzz_0[i] = 3.0 * g_0_xx_0_xxzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxzz_1[i] * fti_ab_0 + 2.0 * g_0_xxx_0_xzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxzz_0[i] * pb_x + g_0_xxx_0_xxzz_1[i] * wp_x[i];

        g_0_xxxx_0_xyyy_0[i] = 3.0 * g_0_xx_0_xyyy_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xyyy_1[i] * fti_ab_0 + g_0_xxx_0_yyy_1[i] * fi_abcd_0 + g_0_xxx_0_xyyy_0[i] * pb_x + g_0_xxx_0_xyyy_1[i] * wp_x[i];

        g_0_xxxx_0_xyyz_0[i] = 3.0 * g_0_xx_0_xyyz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xyyz_1[i] * fti_ab_0 + g_0_xxx_0_yyz_1[i] * fi_abcd_0 + g_0_xxx_0_xyyz_0[i] * pb_x + g_0_xxx_0_xyyz_1[i] * wp_x[i];

        g_0_xxxx_0_xyzz_0[i] = 3.0 * g_0_xx_0_xyzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xyzz_1[i] * fti_ab_0 + g_0_xxx_0_yzz_1[i] * fi_abcd_0 + g_0_xxx_0_xyzz_0[i] * pb_x + g_0_xxx_0_xyzz_1[i] * wp_x[i];

        g_0_xxxx_0_xzzz_0[i] = 3.0 * g_0_xx_0_xzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xzzz_1[i] * fti_ab_0 + g_0_xxx_0_zzz_1[i] * fi_abcd_0 + g_0_xxx_0_xzzz_0[i] * pb_x + g_0_xxx_0_xzzz_1[i] * wp_x[i];

        g_0_xxxx_0_yyyy_0[i] = 3.0 * g_0_xx_0_yyyy_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_yyyy_1[i] * fti_ab_0 + g_0_xxx_0_yyyy_0[i] * pb_x + g_0_xxx_0_yyyy_1[i] * wp_x[i];

        g_0_xxxx_0_yyyz_0[i] = 3.0 * g_0_xx_0_yyyz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_yyyz_1[i] * fti_ab_0 + g_0_xxx_0_yyyz_0[i] * pb_x + g_0_xxx_0_yyyz_1[i] * wp_x[i];

        g_0_xxxx_0_yyzz_0[i] = 3.0 * g_0_xx_0_yyzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_yyzz_1[i] * fti_ab_0 + g_0_xxx_0_yyzz_0[i] * pb_x + g_0_xxx_0_yyzz_1[i] * wp_x[i];

        g_0_xxxx_0_yzzz_0[i] = 3.0 * g_0_xx_0_yzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_yzzz_1[i] * fti_ab_0 + g_0_xxx_0_yzzz_0[i] * pb_x + g_0_xxx_0_yzzz_1[i] * wp_x[i];

        g_0_xxxx_0_zzzz_0[i] = 3.0 * g_0_xx_0_zzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_zzzz_1[i] * fti_ab_0 + g_0_xxx_0_zzzz_0[i] * pb_x + g_0_xxx_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 15-30 components of targeted buffer : SGSG

    auto g_0_xxxy_0_xxxx_0 = pbuffer.data(idx_eri_0_sgsg + 15);

    auto g_0_xxxy_0_xxxy_0 = pbuffer.data(idx_eri_0_sgsg + 16);

    auto g_0_xxxy_0_xxxz_0 = pbuffer.data(idx_eri_0_sgsg + 17);

    auto g_0_xxxy_0_xxyy_0 = pbuffer.data(idx_eri_0_sgsg + 18);

    auto g_0_xxxy_0_xxyz_0 = pbuffer.data(idx_eri_0_sgsg + 19);

    auto g_0_xxxy_0_xxzz_0 = pbuffer.data(idx_eri_0_sgsg + 20);

    auto g_0_xxxy_0_xyyy_0 = pbuffer.data(idx_eri_0_sgsg + 21);

    auto g_0_xxxy_0_xyyz_0 = pbuffer.data(idx_eri_0_sgsg + 22);

    auto g_0_xxxy_0_xyzz_0 = pbuffer.data(idx_eri_0_sgsg + 23);

    auto g_0_xxxy_0_xzzz_0 = pbuffer.data(idx_eri_0_sgsg + 24);

    auto g_0_xxxy_0_yyyy_0 = pbuffer.data(idx_eri_0_sgsg + 25);

    auto g_0_xxxy_0_yyyz_0 = pbuffer.data(idx_eri_0_sgsg + 26);

    auto g_0_xxxy_0_yyzz_0 = pbuffer.data(idx_eri_0_sgsg + 27);

    auto g_0_xxxy_0_yzzz_0 = pbuffer.data(idx_eri_0_sgsg + 28);

    auto g_0_xxxy_0_zzzz_0 = pbuffer.data(idx_eri_0_sgsg + 29);

    #pragma omp simd aligned(g_0_xxx_0_xxx_1, g_0_xxx_0_xxxx_0, g_0_xxx_0_xxxx_1, g_0_xxx_0_xxxy_0, g_0_xxx_0_xxxy_1, g_0_xxx_0_xxxz_0, g_0_xxx_0_xxxz_1, g_0_xxx_0_xxy_1, g_0_xxx_0_xxyy_0, g_0_xxx_0_xxyy_1, g_0_xxx_0_xxyz_0, g_0_xxx_0_xxyz_1, g_0_xxx_0_xxz_1, g_0_xxx_0_xxzz_0, g_0_xxx_0_xxzz_1, g_0_xxx_0_xyy_1, g_0_xxx_0_xyyy_0, g_0_xxx_0_xyyy_1, g_0_xxx_0_xyyz_0, g_0_xxx_0_xyyz_1, g_0_xxx_0_xyz_1, g_0_xxx_0_xyzz_0, g_0_xxx_0_xyzz_1, g_0_xxx_0_xzz_1, g_0_xxx_0_xzzz_0, g_0_xxx_0_xzzz_1, g_0_xxx_0_yyy_1, g_0_xxx_0_yyyy_0, g_0_xxx_0_yyyy_1, g_0_xxx_0_yyyz_0, g_0_xxx_0_yyyz_1, g_0_xxx_0_yyz_1, g_0_xxx_0_yyzz_0, g_0_xxx_0_yyzz_1, g_0_xxx_0_yzz_1, g_0_xxx_0_yzzz_0, g_0_xxx_0_yzzz_1, g_0_xxx_0_zzz_1, g_0_xxx_0_zzzz_0, g_0_xxx_0_zzzz_1, g_0_xxxy_0_xxxx_0, g_0_xxxy_0_xxxy_0, g_0_xxxy_0_xxxz_0, g_0_xxxy_0_xxyy_0, g_0_xxxy_0_xxyz_0, g_0_xxxy_0_xxzz_0, g_0_xxxy_0_xyyy_0, g_0_xxxy_0_xyyz_0, g_0_xxxy_0_xyzz_0, g_0_xxxy_0_xzzz_0, g_0_xxxy_0_yyyy_0, g_0_xxxy_0_yyyz_0, g_0_xxxy_0_yyzz_0, g_0_xxxy_0_yzzz_0, g_0_xxxy_0_zzzz_0, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxy_0_xxxx_0[i] = g_0_xxx_0_xxxx_0[i] * pb_y + g_0_xxx_0_xxxx_1[i] * wp_y[i];

        g_0_xxxy_0_xxxy_0[i] = g_0_xxx_0_xxx_1[i] * fi_abcd_0 + g_0_xxx_0_xxxy_0[i] * pb_y + g_0_xxx_0_xxxy_1[i] * wp_y[i];

        g_0_xxxy_0_xxxz_0[i] = g_0_xxx_0_xxxz_0[i] * pb_y + g_0_xxx_0_xxxz_1[i] * wp_y[i];

        g_0_xxxy_0_xxyy_0[i] = 2.0 * g_0_xxx_0_xxy_1[i] * fi_abcd_0 + g_0_xxx_0_xxyy_0[i] * pb_y + g_0_xxx_0_xxyy_1[i] * wp_y[i];

        g_0_xxxy_0_xxyz_0[i] = g_0_xxx_0_xxz_1[i] * fi_abcd_0 + g_0_xxx_0_xxyz_0[i] * pb_y + g_0_xxx_0_xxyz_1[i] * wp_y[i];

        g_0_xxxy_0_xxzz_0[i] = g_0_xxx_0_xxzz_0[i] * pb_y + g_0_xxx_0_xxzz_1[i] * wp_y[i];

        g_0_xxxy_0_xyyy_0[i] = 3.0 * g_0_xxx_0_xyy_1[i] * fi_abcd_0 + g_0_xxx_0_xyyy_0[i] * pb_y + g_0_xxx_0_xyyy_1[i] * wp_y[i];

        g_0_xxxy_0_xyyz_0[i] = 2.0 * g_0_xxx_0_xyz_1[i] * fi_abcd_0 + g_0_xxx_0_xyyz_0[i] * pb_y + g_0_xxx_0_xyyz_1[i] * wp_y[i];

        g_0_xxxy_0_xyzz_0[i] = g_0_xxx_0_xzz_1[i] * fi_abcd_0 + g_0_xxx_0_xyzz_0[i] * pb_y + g_0_xxx_0_xyzz_1[i] * wp_y[i];

        g_0_xxxy_0_xzzz_0[i] = g_0_xxx_0_xzzz_0[i] * pb_y + g_0_xxx_0_xzzz_1[i] * wp_y[i];

        g_0_xxxy_0_yyyy_0[i] = 4.0 * g_0_xxx_0_yyy_1[i] * fi_abcd_0 + g_0_xxx_0_yyyy_0[i] * pb_y + g_0_xxx_0_yyyy_1[i] * wp_y[i];

        g_0_xxxy_0_yyyz_0[i] = 3.0 * g_0_xxx_0_yyz_1[i] * fi_abcd_0 + g_0_xxx_0_yyyz_0[i] * pb_y + g_0_xxx_0_yyyz_1[i] * wp_y[i];

        g_0_xxxy_0_yyzz_0[i] = 2.0 * g_0_xxx_0_yzz_1[i] * fi_abcd_0 + g_0_xxx_0_yyzz_0[i] * pb_y + g_0_xxx_0_yyzz_1[i] * wp_y[i];

        g_0_xxxy_0_yzzz_0[i] = g_0_xxx_0_zzz_1[i] * fi_abcd_0 + g_0_xxx_0_yzzz_0[i] * pb_y + g_0_xxx_0_yzzz_1[i] * wp_y[i];

        g_0_xxxy_0_zzzz_0[i] = g_0_xxx_0_zzzz_0[i] * pb_y + g_0_xxx_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 30-45 components of targeted buffer : SGSG

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

    auto g_0_xxxz_0_yyyy_0 = pbuffer.data(idx_eri_0_sgsg + 40);

    auto g_0_xxxz_0_yyyz_0 = pbuffer.data(idx_eri_0_sgsg + 41);

    auto g_0_xxxz_0_yyzz_0 = pbuffer.data(idx_eri_0_sgsg + 42);

    auto g_0_xxxz_0_yzzz_0 = pbuffer.data(idx_eri_0_sgsg + 43);

    auto g_0_xxxz_0_zzzz_0 = pbuffer.data(idx_eri_0_sgsg + 44);

    #pragma omp simd aligned(g_0_xxx_0_xxx_1, g_0_xxx_0_xxxx_0, g_0_xxx_0_xxxx_1, g_0_xxx_0_xxxy_0, g_0_xxx_0_xxxy_1, g_0_xxx_0_xxxz_0, g_0_xxx_0_xxxz_1, g_0_xxx_0_xxy_1, g_0_xxx_0_xxyy_0, g_0_xxx_0_xxyy_1, g_0_xxx_0_xxyz_0, g_0_xxx_0_xxyz_1, g_0_xxx_0_xxz_1, g_0_xxx_0_xxzz_0, g_0_xxx_0_xxzz_1, g_0_xxx_0_xyy_1, g_0_xxx_0_xyyy_0, g_0_xxx_0_xyyy_1, g_0_xxx_0_xyyz_0, g_0_xxx_0_xyyz_1, g_0_xxx_0_xyz_1, g_0_xxx_0_xyzz_0, g_0_xxx_0_xyzz_1, g_0_xxx_0_xzz_1, g_0_xxx_0_xzzz_0, g_0_xxx_0_xzzz_1, g_0_xxx_0_yyy_1, g_0_xxx_0_yyyy_0, g_0_xxx_0_yyyy_1, g_0_xxx_0_yyyz_0, g_0_xxx_0_yyyz_1, g_0_xxx_0_yyz_1, g_0_xxx_0_yyzz_0, g_0_xxx_0_yyzz_1, g_0_xxx_0_yzz_1, g_0_xxx_0_yzzz_0, g_0_xxx_0_yzzz_1, g_0_xxx_0_zzz_1, g_0_xxx_0_zzzz_0, g_0_xxx_0_zzzz_1, g_0_xxxz_0_xxxx_0, g_0_xxxz_0_xxxy_0, g_0_xxxz_0_xxxz_0, g_0_xxxz_0_xxyy_0, g_0_xxxz_0_xxyz_0, g_0_xxxz_0_xxzz_0, g_0_xxxz_0_xyyy_0, g_0_xxxz_0_xyyz_0, g_0_xxxz_0_xyzz_0, g_0_xxxz_0_xzzz_0, g_0_xxxz_0_yyyy_0, g_0_xxxz_0_yyyz_0, g_0_xxxz_0_yyzz_0, g_0_xxxz_0_yzzz_0, g_0_xxxz_0_zzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxz_0_xxxx_0[i] = g_0_xxx_0_xxxx_0[i] * pb_z + g_0_xxx_0_xxxx_1[i] * wp_z[i];

        g_0_xxxz_0_xxxy_0[i] = g_0_xxx_0_xxxy_0[i] * pb_z + g_0_xxx_0_xxxy_1[i] * wp_z[i];

        g_0_xxxz_0_xxxz_0[i] = g_0_xxx_0_xxx_1[i] * fi_abcd_0 + g_0_xxx_0_xxxz_0[i] * pb_z + g_0_xxx_0_xxxz_1[i] * wp_z[i];

        g_0_xxxz_0_xxyy_0[i] = g_0_xxx_0_xxyy_0[i] * pb_z + g_0_xxx_0_xxyy_1[i] * wp_z[i];

        g_0_xxxz_0_xxyz_0[i] = g_0_xxx_0_xxy_1[i] * fi_abcd_0 + g_0_xxx_0_xxyz_0[i] * pb_z + g_0_xxx_0_xxyz_1[i] * wp_z[i];

        g_0_xxxz_0_xxzz_0[i] = 2.0 * g_0_xxx_0_xxz_1[i] * fi_abcd_0 + g_0_xxx_0_xxzz_0[i] * pb_z + g_0_xxx_0_xxzz_1[i] * wp_z[i];

        g_0_xxxz_0_xyyy_0[i] = g_0_xxx_0_xyyy_0[i] * pb_z + g_0_xxx_0_xyyy_1[i] * wp_z[i];

        g_0_xxxz_0_xyyz_0[i] = g_0_xxx_0_xyy_1[i] * fi_abcd_0 + g_0_xxx_0_xyyz_0[i] * pb_z + g_0_xxx_0_xyyz_1[i] * wp_z[i];

        g_0_xxxz_0_xyzz_0[i] = 2.0 * g_0_xxx_0_xyz_1[i] * fi_abcd_0 + g_0_xxx_0_xyzz_0[i] * pb_z + g_0_xxx_0_xyzz_1[i] * wp_z[i];

        g_0_xxxz_0_xzzz_0[i] = 3.0 * g_0_xxx_0_xzz_1[i] * fi_abcd_0 + g_0_xxx_0_xzzz_0[i] * pb_z + g_0_xxx_0_xzzz_1[i] * wp_z[i];

        g_0_xxxz_0_yyyy_0[i] = g_0_xxx_0_yyyy_0[i] * pb_z + g_0_xxx_0_yyyy_1[i] * wp_z[i];

        g_0_xxxz_0_yyyz_0[i] = g_0_xxx_0_yyy_1[i] * fi_abcd_0 + g_0_xxx_0_yyyz_0[i] * pb_z + g_0_xxx_0_yyyz_1[i] * wp_z[i];

        g_0_xxxz_0_yyzz_0[i] = 2.0 * g_0_xxx_0_yyz_1[i] * fi_abcd_0 + g_0_xxx_0_yyzz_0[i] * pb_z + g_0_xxx_0_yyzz_1[i] * wp_z[i];

        g_0_xxxz_0_yzzz_0[i] = 3.0 * g_0_xxx_0_yzz_1[i] * fi_abcd_0 + g_0_xxx_0_yzzz_0[i] * pb_z + g_0_xxx_0_yzzz_1[i] * wp_z[i];

        g_0_xxxz_0_zzzz_0[i] = 4.0 * g_0_xxx_0_zzz_1[i] * fi_abcd_0 + g_0_xxx_0_zzzz_0[i] * pb_z + g_0_xxx_0_zzzz_1[i] * wp_z[i];
    }

    /// Set up 45-60 components of targeted buffer : SGSG

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

    #pragma omp simd aligned(g_0_xx_0_xxxx_0, g_0_xx_0_xxxx_1, g_0_xx_0_xxxz_0, g_0_xx_0_xxxz_1, g_0_xx_0_xxzz_0, g_0_xx_0_xxzz_1, g_0_xx_0_xzzz_0, g_0_xx_0_xzzz_1, g_0_xxy_0_xxxx_0, g_0_xxy_0_xxxx_1, g_0_xxy_0_xxxz_0, g_0_xxy_0_xxxz_1, g_0_xxy_0_xxzz_0, g_0_xxy_0_xxzz_1, g_0_xxy_0_xzzz_0, g_0_xxy_0_xzzz_1, g_0_xxyy_0_xxxx_0, g_0_xxyy_0_xxxy_0, g_0_xxyy_0_xxxz_0, g_0_xxyy_0_xxyy_0, g_0_xxyy_0_xxyz_0, g_0_xxyy_0_xxzz_0, g_0_xxyy_0_xyyy_0, g_0_xxyy_0_xyyz_0, g_0_xxyy_0_xyzz_0, g_0_xxyy_0_xzzz_0, g_0_xxyy_0_yyyy_0, g_0_xxyy_0_yyyz_0, g_0_xxyy_0_yyzz_0, g_0_xxyy_0_yzzz_0, g_0_xxyy_0_zzzz_0, g_0_xyy_0_xxxy_0, g_0_xyy_0_xxxy_1, g_0_xyy_0_xxy_1, g_0_xyy_0_xxyy_0, g_0_xyy_0_xxyy_1, g_0_xyy_0_xxyz_0, g_0_xyy_0_xxyz_1, g_0_xyy_0_xyy_1, g_0_xyy_0_xyyy_0, g_0_xyy_0_xyyy_1, g_0_xyy_0_xyyz_0, g_0_xyy_0_xyyz_1, g_0_xyy_0_xyz_1, g_0_xyy_0_xyzz_0, g_0_xyy_0_xyzz_1, g_0_xyy_0_yyy_1, g_0_xyy_0_yyyy_0, g_0_xyy_0_yyyy_1, g_0_xyy_0_yyyz_0, g_0_xyy_0_yyyz_1, g_0_xyy_0_yyz_1, g_0_xyy_0_yyzz_0, g_0_xyy_0_yyzz_1, g_0_xyy_0_yzz_1, g_0_xyy_0_yzzz_0, g_0_xyy_0_yzzz_1, g_0_xyy_0_zzzz_0, g_0_xyy_0_zzzz_1, g_0_yy_0_xxxy_0, g_0_yy_0_xxxy_1, g_0_yy_0_xxyy_0, g_0_yy_0_xxyy_1, g_0_yy_0_xxyz_0, g_0_yy_0_xxyz_1, g_0_yy_0_xyyy_0, g_0_yy_0_xyyy_1, g_0_yy_0_xyyz_0, g_0_yy_0_xyyz_1, g_0_yy_0_xyzz_0, g_0_yy_0_xyzz_1, g_0_yy_0_yyyy_0, g_0_yy_0_yyyy_1, g_0_yy_0_yyyz_0, g_0_yy_0_yyyz_1, g_0_yy_0_yyzz_0, g_0_yy_0_yyzz_1, g_0_yy_0_yzzz_0, g_0_yy_0_yzzz_1, g_0_yy_0_zzzz_0, g_0_yy_0_zzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyy_0_xxxx_0[i] = g_0_xx_0_xxxx_0[i] * fi_ab_0 - g_0_xx_0_xxxx_1[i] * fti_ab_0 + g_0_xxy_0_xxxx_0[i] * pb_y + g_0_xxy_0_xxxx_1[i] * wp_y[i];

        g_0_xxyy_0_xxxy_0[i] = g_0_yy_0_xxxy_0[i] * fi_ab_0 - g_0_yy_0_xxxy_1[i] * fti_ab_0 + 3.0 * g_0_xyy_0_xxy_1[i] * fi_abcd_0 + g_0_xyy_0_xxxy_0[i] * pb_x + g_0_xyy_0_xxxy_1[i] * wp_x[i];

        g_0_xxyy_0_xxxz_0[i] = g_0_xx_0_xxxz_0[i] * fi_ab_0 - g_0_xx_0_xxxz_1[i] * fti_ab_0 + g_0_xxy_0_xxxz_0[i] * pb_y + g_0_xxy_0_xxxz_1[i] * wp_y[i];

        g_0_xxyy_0_xxyy_0[i] = g_0_yy_0_xxyy_0[i] * fi_ab_0 - g_0_yy_0_xxyy_1[i] * fti_ab_0 + 2.0 * g_0_xyy_0_xyy_1[i] * fi_abcd_0 + g_0_xyy_0_xxyy_0[i] * pb_x + g_0_xyy_0_xxyy_1[i] * wp_x[i];

        g_0_xxyy_0_xxyz_0[i] = g_0_yy_0_xxyz_0[i] * fi_ab_0 - g_0_yy_0_xxyz_1[i] * fti_ab_0 + 2.0 * g_0_xyy_0_xyz_1[i] * fi_abcd_0 + g_0_xyy_0_xxyz_0[i] * pb_x + g_0_xyy_0_xxyz_1[i] * wp_x[i];

        g_0_xxyy_0_xxzz_0[i] = g_0_xx_0_xxzz_0[i] * fi_ab_0 - g_0_xx_0_xxzz_1[i] * fti_ab_0 + g_0_xxy_0_xxzz_0[i] * pb_y + g_0_xxy_0_xxzz_1[i] * wp_y[i];

        g_0_xxyy_0_xyyy_0[i] = g_0_yy_0_xyyy_0[i] * fi_ab_0 - g_0_yy_0_xyyy_1[i] * fti_ab_0 + g_0_xyy_0_yyy_1[i] * fi_abcd_0 + g_0_xyy_0_xyyy_0[i] * pb_x + g_0_xyy_0_xyyy_1[i] * wp_x[i];

        g_0_xxyy_0_xyyz_0[i] = g_0_yy_0_xyyz_0[i] * fi_ab_0 - g_0_yy_0_xyyz_1[i] * fti_ab_0 + g_0_xyy_0_yyz_1[i] * fi_abcd_0 + g_0_xyy_0_xyyz_0[i] * pb_x + g_0_xyy_0_xyyz_1[i] * wp_x[i];

        g_0_xxyy_0_xyzz_0[i] = g_0_yy_0_xyzz_0[i] * fi_ab_0 - g_0_yy_0_xyzz_1[i] * fti_ab_0 + g_0_xyy_0_yzz_1[i] * fi_abcd_0 + g_0_xyy_0_xyzz_0[i] * pb_x + g_0_xyy_0_xyzz_1[i] * wp_x[i];

        g_0_xxyy_0_xzzz_0[i] = g_0_xx_0_xzzz_0[i] * fi_ab_0 - g_0_xx_0_xzzz_1[i] * fti_ab_0 + g_0_xxy_0_xzzz_0[i] * pb_y + g_0_xxy_0_xzzz_1[i] * wp_y[i];

        g_0_xxyy_0_yyyy_0[i] = g_0_yy_0_yyyy_0[i] * fi_ab_0 - g_0_yy_0_yyyy_1[i] * fti_ab_0 + g_0_xyy_0_yyyy_0[i] * pb_x + g_0_xyy_0_yyyy_1[i] * wp_x[i];

        g_0_xxyy_0_yyyz_0[i] = g_0_yy_0_yyyz_0[i] * fi_ab_0 - g_0_yy_0_yyyz_1[i] * fti_ab_0 + g_0_xyy_0_yyyz_0[i] * pb_x + g_0_xyy_0_yyyz_1[i] * wp_x[i];

        g_0_xxyy_0_yyzz_0[i] = g_0_yy_0_yyzz_0[i] * fi_ab_0 - g_0_yy_0_yyzz_1[i] * fti_ab_0 + g_0_xyy_0_yyzz_0[i] * pb_x + g_0_xyy_0_yyzz_1[i] * wp_x[i];

        g_0_xxyy_0_yzzz_0[i] = g_0_yy_0_yzzz_0[i] * fi_ab_0 - g_0_yy_0_yzzz_1[i] * fti_ab_0 + g_0_xyy_0_yzzz_0[i] * pb_x + g_0_xyy_0_yzzz_1[i] * wp_x[i];

        g_0_xxyy_0_zzzz_0[i] = g_0_yy_0_zzzz_0[i] * fi_ab_0 - g_0_yy_0_zzzz_1[i] * fti_ab_0 + g_0_xyy_0_zzzz_0[i] * pb_x + g_0_xyy_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 60-75 components of targeted buffer : SGSG

    auto g_0_xxyz_0_xxxx_0 = pbuffer.data(idx_eri_0_sgsg + 60);

    auto g_0_xxyz_0_xxxy_0 = pbuffer.data(idx_eri_0_sgsg + 61);

    auto g_0_xxyz_0_xxxz_0 = pbuffer.data(idx_eri_0_sgsg + 62);

    auto g_0_xxyz_0_xxyy_0 = pbuffer.data(idx_eri_0_sgsg + 63);

    auto g_0_xxyz_0_xxyz_0 = pbuffer.data(idx_eri_0_sgsg + 64);

    auto g_0_xxyz_0_xxzz_0 = pbuffer.data(idx_eri_0_sgsg + 65);

    auto g_0_xxyz_0_xyyy_0 = pbuffer.data(idx_eri_0_sgsg + 66);

    auto g_0_xxyz_0_xyyz_0 = pbuffer.data(idx_eri_0_sgsg + 67);

    auto g_0_xxyz_0_xyzz_0 = pbuffer.data(idx_eri_0_sgsg + 68);

    auto g_0_xxyz_0_xzzz_0 = pbuffer.data(idx_eri_0_sgsg + 69);

    auto g_0_xxyz_0_yyyy_0 = pbuffer.data(idx_eri_0_sgsg + 70);

    auto g_0_xxyz_0_yyyz_0 = pbuffer.data(idx_eri_0_sgsg + 71);

    auto g_0_xxyz_0_yyzz_0 = pbuffer.data(idx_eri_0_sgsg + 72);

    auto g_0_xxyz_0_yzzz_0 = pbuffer.data(idx_eri_0_sgsg + 73);

    auto g_0_xxyz_0_zzzz_0 = pbuffer.data(idx_eri_0_sgsg + 74);

    #pragma omp simd aligned(g_0_xxy_0_xxxy_0, g_0_xxy_0_xxxy_1, g_0_xxy_0_xxyy_0, g_0_xxy_0_xxyy_1, g_0_xxy_0_xyyy_0, g_0_xxy_0_xyyy_1, g_0_xxy_0_yyyy_0, g_0_xxy_0_yyyy_1, g_0_xxyz_0_xxxx_0, g_0_xxyz_0_xxxy_0, g_0_xxyz_0_xxxz_0, g_0_xxyz_0_xxyy_0, g_0_xxyz_0_xxyz_0, g_0_xxyz_0_xxzz_0, g_0_xxyz_0_xyyy_0, g_0_xxyz_0_xyyz_0, g_0_xxyz_0_xyzz_0, g_0_xxyz_0_xzzz_0, g_0_xxyz_0_yyyy_0, g_0_xxyz_0_yyyz_0, g_0_xxyz_0_yyzz_0, g_0_xxyz_0_yzzz_0, g_0_xxyz_0_zzzz_0, g_0_xxz_0_xxxx_0, g_0_xxz_0_xxxx_1, g_0_xxz_0_xxxz_0, g_0_xxz_0_xxxz_1, g_0_xxz_0_xxyz_0, g_0_xxz_0_xxyz_1, g_0_xxz_0_xxz_1, g_0_xxz_0_xxzz_0, g_0_xxz_0_xxzz_1, g_0_xxz_0_xyyz_0, g_0_xxz_0_xyyz_1, g_0_xxz_0_xyz_1, g_0_xxz_0_xyzz_0, g_0_xxz_0_xyzz_1, g_0_xxz_0_xzz_1, g_0_xxz_0_xzzz_0, g_0_xxz_0_xzzz_1, g_0_xxz_0_yyyz_0, g_0_xxz_0_yyyz_1, g_0_xxz_0_yyz_1, g_0_xxz_0_yyzz_0, g_0_xxz_0_yyzz_1, g_0_xxz_0_yzz_1, g_0_xxz_0_yzzz_0, g_0_xxz_0_yzzz_1, g_0_xxz_0_zzz_1, g_0_xxz_0_zzzz_0, g_0_xxz_0_zzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyz_0_xxxx_0[i] = g_0_xxz_0_xxxx_0[i] * pb_y + g_0_xxz_0_xxxx_1[i] * wp_y[i];

        g_0_xxyz_0_xxxy_0[i] = g_0_xxy_0_xxxy_0[i] * pb_z + g_0_xxy_0_xxxy_1[i] * wp_z[i];

        g_0_xxyz_0_xxxz_0[i] = g_0_xxz_0_xxxz_0[i] * pb_y + g_0_xxz_0_xxxz_1[i] * wp_y[i];

        g_0_xxyz_0_xxyy_0[i] = g_0_xxy_0_xxyy_0[i] * pb_z + g_0_xxy_0_xxyy_1[i] * wp_z[i];

        g_0_xxyz_0_xxyz_0[i] = g_0_xxz_0_xxz_1[i] * fi_abcd_0 + g_0_xxz_0_xxyz_0[i] * pb_y + g_0_xxz_0_xxyz_1[i] * wp_y[i];

        g_0_xxyz_0_xxzz_0[i] = g_0_xxz_0_xxzz_0[i] * pb_y + g_0_xxz_0_xxzz_1[i] * wp_y[i];

        g_0_xxyz_0_xyyy_0[i] = g_0_xxy_0_xyyy_0[i] * pb_z + g_0_xxy_0_xyyy_1[i] * wp_z[i];

        g_0_xxyz_0_xyyz_0[i] = 2.0 * g_0_xxz_0_xyz_1[i] * fi_abcd_0 + g_0_xxz_0_xyyz_0[i] * pb_y + g_0_xxz_0_xyyz_1[i] * wp_y[i];

        g_0_xxyz_0_xyzz_0[i] = g_0_xxz_0_xzz_1[i] * fi_abcd_0 + g_0_xxz_0_xyzz_0[i] * pb_y + g_0_xxz_0_xyzz_1[i] * wp_y[i];

        g_0_xxyz_0_xzzz_0[i] = g_0_xxz_0_xzzz_0[i] * pb_y + g_0_xxz_0_xzzz_1[i] * wp_y[i];

        g_0_xxyz_0_yyyy_0[i] = g_0_xxy_0_yyyy_0[i] * pb_z + g_0_xxy_0_yyyy_1[i] * wp_z[i];

        g_0_xxyz_0_yyyz_0[i] = 3.0 * g_0_xxz_0_yyz_1[i] * fi_abcd_0 + g_0_xxz_0_yyyz_0[i] * pb_y + g_0_xxz_0_yyyz_1[i] * wp_y[i];

        g_0_xxyz_0_yyzz_0[i] = 2.0 * g_0_xxz_0_yzz_1[i] * fi_abcd_0 + g_0_xxz_0_yyzz_0[i] * pb_y + g_0_xxz_0_yyzz_1[i] * wp_y[i];

        g_0_xxyz_0_yzzz_0[i] = g_0_xxz_0_zzz_1[i] * fi_abcd_0 + g_0_xxz_0_yzzz_0[i] * pb_y + g_0_xxz_0_yzzz_1[i] * wp_y[i];

        g_0_xxyz_0_zzzz_0[i] = g_0_xxz_0_zzzz_0[i] * pb_y + g_0_xxz_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 75-90 components of targeted buffer : SGSG

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

    #pragma omp simd aligned(g_0_xx_0_xxxx_0, g_0_xx_0_xxxx_1, g_0_xx_0_xxxy_0, g_0_xx_0_xxxy_1, g_0_xx_0_xxyy_0, g_0_xx_0_xxyy_1, g_0_xx_0_xyyy_0, g_0_xx_0_xyyy_1, g_0_xxz_0_xxxx_0, g_0_xxz_0_xxxx_1, g_0_xxz_0_xxxy_0, g_0_xxz_0_xxxy_1, g_0_xxz_0_xxyy_0, g_0_xxz_0_xxyy_1, g_0_xxz_0_xyyy_0, g_0_xxz_0_xyyy_1, g_0_xxzz_0_xxxx_0, g_0_xxzz_0_xxxy_0, g_0_xxzz_0_xxxz_0, g_0_xxzz_0_xxyy_0, g_0_xxzz_0_xxyz_0, g_0_xxzz_0_xxzz_0, g_0_xxzz_0_xyyy_0, g_0_xxzz_0_xyyz_0, g_0_xxzz_0_xyzz_0, g_0_xxzz_0_xzzz_0, g_0_xxzz_0_yyyy_0, g_0_xxzz_0_yyyz_0, g_0_xxzz_0_yyzz_0, g_0_xxzz_0_yzzz_0, g_0_xxzz_0_zzzz_0, g_0_xzz_0_xxxz_0, g_0_xzz_0_xxxz_1, g_0_xzz_0_xxyz_0, g_0_xzz_0_xxyz_1, g_0_xzz_0_xxz_1, g_0_xzz_0_xxzz_0, g_0_xzz_0_xxzz_1, g_0_xzz_0_xyyz_0, g_0_xzz_0_xyyz_1, g_0_xzz_0_xyz_1, g_0_xzz_0_xyzz_0, g_0_xzz_0_xyzz_1, g_0_xzz_0_xzz_1, g_0_xzz_0_xzzz_0, g_0_xzz_0_xzzz_1, g_0_xzz_0_yyyy_0, g_0_xzz_0_yyyy_1, g_0_xzz_0_yyyz_0, g_0_xzz_0_yyyz_1, g_0_xzz_0_yyz_1, g_0_xzz_0_yyzz_0, g_0_xzz_0_yyzz_1, g_0_xzz_0_yzz_1, g_0_xzz_0_yzzz_0, g_0_xzz_0_yzzz_1, g_0_xzz_0_zzz_1, g_0_xzz_0_zzzz_0, g_0_xzz_0_zzzz_1, g_0_zz_0_xxxz_0, g_0_zz_0_xxxz_1, g_0_zz_0_xxyz_0, g_0_zz_0_xxyz_1, g_0_zz_0_xxzz_0, g_0_zz_0_xxzz_1, g_0_zz_0_xyyz_0, g_0_zz_0_xyyz_1, g_0_zz_0_xyzz_0, g_0_zz_0_xyzz_1, g_0_zz_0_xzzz_0, g_0_zz_0_xzzz_1, g_0_zz_0_yyyy_0, g_0_zz_0_yyyy_1, g_0_zz_0_yyyz_0, g_0_zz_0_yyyz_1, g_0_zz_0_yyzz_0, g_0_zz_0_yyzz_1, g_0_zz_0_yzzz_0, g_0_zz_0_yzzz_1, g_0_zz_0_zzzz_0, g_0_zz_0_zzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxzz_0_xxxx_0[i] = g_0_xx_0_xxxx_0[i] * fi_ab_0 - g_0_xx_0_xxxx_1[i] * fti_ab_0 + g_0_xxz_0_xxxx_0[i] * pb_z + g_0_xxz_0_xxxx_1[i] * wp_z[i];

        g_0_xxzz_0_xxxy_0[i] = g_0_xx_0_xxxy_0[i] * fi_ab_0 - g_0_xx_0_xxxy_1[i] * fti_ab_0 + g_0_xxz_0_xxxy_0[i] * pb_z + g_0_xxz_0_xxxy_1[i] * wp_z[i];

        g_0_xxzz_0_xxxz_0[i] = g_0_zz_0_xxxz_0[i] * fi_ab_0 - g_0_zz_0_xxxz_1[i] * fti_ab_0 + 3.0 * g_0_xzz_0_xxz_1[i] * fi_abcd_0 + g_0_xzz_0_xxxz_0[i] * pb_x + g_0_xzz_0_xxxz_1[i] * wp_x[i];

        g_0_xxzz_0_xxyy_0[i] = g_0_xx_0_xxyy_0[i] * fi_ab_0 - g_0_xx_0_xxyy_1[i] * fti_ab_0 + g_0_xxz_0_xxyy_0[i] * pb_z + g_0_xxz_0_xxyy_1[i] * wp_z[i];

        g_0_xxzz_0_xxyz_0[i] = g_0_zz_0_xxyz_0[i] * fi_ab_0 - g_0_zz_0_xxyz_1[i] * fti_ab_0 + 2.0 * g_0_xzz_0_xyz_1[i] * fi_abcd_0 + g_0_xzz_0_xxyz_0[i] * pb_x + g_0_xzz_0_xxyz_1[i] * wp_x[i];

        g_0_xxzz_0_xxzz_0[i] = g_0_zz_0_xxzz_0[i] * fi_ab_0 - g_0_zz_0_xxzz_1[i] * fti_ab_0 + 2.0 * g_0_xzz_0_xzz_1[i] * fi_abcd_0 + g_0_xzz_0_xxzz_0[i] * pb_x + g_0_xzz_0_xxzz_1[i] * wp_x[i];

        g_0_xxzz_0_xyyy_0[i] = g_0_xx_0_xyyy_0[i] * fi_ab_0 - g_0_xx_0_xyyy_1[i] * fti_ab_0 + g_0_xxz_0_xyyy_0[i] * pb_z + g_0_xxz_0_xyyy_1[i] * wp_z[i];

        g_0_xxzz_0_xyyz_0[i] = g_0_zz_0_xyyz_0[i] * fi_ab_0 - g_0_zz_0_xyyz_1[i] * fti_ab_0 + g_0_xzz_0_yyz_1[i] * fi_abcd_0 + g_0_xzz_0_xyyz_0[i] * pb_x + g_0_xzz_0_xyyz_1[i] * wp_x[i];

        g_0_xxzz_0_xyzz_0[i] = g_0_zz_0_xyzz_0[i] * fi_ab_0 - g_0_zz_0_xyzz_1[i] * fti_ab_0 + g_0_xzz_0_yzz_1[i] * fi_abcd_0 + g_0_xzz_0_xyzz_0[i] * pb_x + g_0_xzz_0_xyzz_1[i] * wp_x[i];

        g_0_xxzz_0_xzzz_0[i] = g_0_zz_0_xzzz_0[i] * fi_ab_0 - g_0_zz_0_xzzz_1[i] * fti_ab_0 + g_0_xzz_0_zzz_1[i] * fi_abcd_0 + g_0_xzz_0_xzzz_0[i] * pb_x + g_0_xzz_0_xzzz_1[i] * wp_x[i];

        g_0_xxzz_0_yyyy_0[i] = g_0_zz_0_yyyy_0[i] * fi_ab_0 - g_0_zz_0_yyyy_1[i] * fti_ab_0 + g_0_xzz_0_yyyy_0[i] * pb_x + g_0_xzz_0_yyyy_1[i] * wp_x[i];

        g_0_xxzz_0_yyyz_0[i] = g_0_zz_0_yyyz_0[i] * fi_ab_0 - g_0_zz_0_yyyz_1[i] * fti_ab_0 + g_0_xzz_0_yyyz_0[i] * pb_x + g_0_xzz_0_yyyz_1[i] * wp_x[i];

        g_0_xxzz_0_yyzz_0[i] = g_0_zz_0_yyzz_0[i] * fi_ab_0 - g_0_zz_0_yyzz_1[i] * fti_ab_0 + g_0_xzz_0_yyzz_0[i] * pb_x + g_0_xzz_0_yyzz_1[i] * wp_x[i];

        g_0_xxzz_0_yzzz_0[i] = g_0_zz_0_yzzz_0[i] * fi_ab_0 - g_0_zz_0_yzzz_1[i] * fti_ab_0 + g_0_xzz_0_yzzz_0[i] * pb_x + g_0_xzz_0_yzzz_1[i] * wp_x[i];

        g_0_xxzz_0_zzzz_0[i] = g_0_zz_0_zzzz_0[i] * fi_ab_0 - g_0_zz_0_zzzz_1[i] * fti_ab_0 + g_0_xzz_0_zzzz_0[i] * pb_x + g_0_xzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 90-105 components of targeted buffer : SGSG

    auto g_0_xyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_sgsg + 90);

    auto g_0_xyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_sgsg + 91);

    auto g_0_xyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_sgsg + 92);

    auto g_0_xyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_sgsg + 93);

    auto g_0_xyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_sgsg + 94);

    auto g_0_xyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_sgsg + 95);

    auto g_0_xyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_sgsg + 96);

    auto g_0_xyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_sgsg + 97);

    auto g_0_xyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_sgsg + 98);

    auto g_0_xyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_sgsg + 99);

    auto g_0_xyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_sgsg + 100);

    auto g_0_xyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_sgsg + 101);

    auto g_0_xyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_sgsg + 102);

    auto g_0_xyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_sgsg + 103);

    auto g_0_xyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_sgsg + 104);

    #pragma omp simd aligned(g_0_xyyy_0_xxxx_0, g_0_xyyy_0_xxxy_0, g_0_xyyy_0_xxxz_0, g_0_xyyy_0_xxyy_0, g_0_xyyy_0_xxyz_0, g_0_xyyy_0_xxzz_0, g_0_xyyy_0_xyyy_0, g_0_xyyy_0_xyyz_0, g_0_xyyy_0_xyzz_0, g_0_xyyy_0_xzzz_0, g_0_xyyy_0_yyyy_0, g_0_xyyy_0_yyyz_0, g_0_xyyy_0_yyzz_0, g_0_xyyy_0_yzzz_0, g_0_xyyy_0_zzzz_0, g_0_yyy_0_xxx_1, g_0_yyy_0_xxxx_0, g_0_yyy_0_xxxx_1, g_0_yyy_0_xxxy_0, g_0_yyy_0_xxxy_1, g_0_yyy_0_xxxz_0, g_0_yyy_0_xxxz_1, g_0_yyy_0_xxy_1, g_0_yyy_0_xxyy_0, g_0_yyy_0_xxyy_1, g_0_yyy_0_xxyz_0, g_0_yyy_0_xxyz_1, g_0_yyy_0_xxz_1, g_0_yyy_0_xxzz_0, g_0_yyy_0_xxzz_1, g_0_yyy_0_xyy_1, g_0_yyy_0_xyyy_0, g_0_yyy_0_xyyy_1, g_0_yyy_0_xyyz_0, g_0_yyy_0_xyyz_1, g_0_yyy_0_xyz_1, g_0_yyy_0_xyzz_0, g_0_yyy_0_xyzz_1, g_0_yyy_0_xzz_1, g_0_yyy_0_xzzz_0, g_0_yyy_0_xzzz_1, g_0_yyy_0_yyy_1, g_0_yyy_0_yyyy_0, g_0_yyy_0_yyyy_1, g_0_yyy_0_yyyz_0, g_0_yyy_0_yyyz_1, g_0_yyy_0_yyz_1, g_0_yyy_0_yyzz_0, g_0_yyy_0_yyzz_1, g_0_yyy_0_yzz_1, g_0_yyy_0_yzzz_0, g_0_yyy_0_yzzz_1, g_0_yyy_0_zzz_1, g_0_yyy_0_zzzz_0, g_0_yyy_0_zzzz_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyy_0_xxxx_0[i] = 4.0 * g_0_yyy_0_xxx_1[i] * fi_abcd_0 + g_0_yyy_0_xxxx_0[i] * pb_x + g_0_yyy_0_xxxx_1[i] * wp_x[i];

        g_0_xyyy_0_xxxy_0[i] = 3.0 * g_0_yyy_0_xxy_1[i] * fi_abcd_0 + g_0_yyy_0_xxxy_0[i] * pb_x + g_0_yyy_0_xxxy_1[i] * wp_x[i];

        g_0_xyyy_0_xxxz_0[i] = 3.0 * g_0_yyy_0_xxz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxz_0[i] * pb_x + g_0_yyy_0_xxxz_1[i] * wp_x[i];

        g_0_xyyy_0_xxyy_0[i] = 2.0 * g_0_yyy_0_xyy_1[i] * fi_abcd_0 + g_0_yyy_0_xxyy_0[i] * pb_x + g_0_yyy_0_xxyy_1[i] * wp_x[i];

        g_0_xyyy_0_xxyz_0[i] = 2.0 * g_0_yyy_0_xyz_1[i] * fi_abcd_0 + g_0_yyy_0_xxyz_0[i] * pb_x + g_0_yyy_0_xxyz_1[i] * wp_x[i];

        g_0_xyyy_0_xxzz_0[i] = 2.0 * g_0_yyy_0_xzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxzz_0[i] * pb_x + g_0_yyy_0_xxzz_1[i] * wp_x[i];

        g_0_xyyy_0_xyyy_0[i] = g_0_yyy_0_yyy_1[i] * fi_abcd_0 + g_0_yyy_0_xyyy_0[i] * pb_x + g_0_yyy_0_xyyy_1[i] * wp_x[i];

        g_0_xyyy_0_xyyz_0[i] = g_0_yyy_0_yyz_1[i] * fi_abcd_0 + g_0_yyy_0_xyyz_0[i] * pb_x + g_0_yyy_0_xyyz_1[i] * wp_x[i];

        g_0_xyyy_0_xyzz_0[i] = g_0_yyy_0_yzz_1[i] * fi_abcd_0 + g_0_yyy_0_xyzz_0[i] * pb_x + g_0_yyy_0_xyzz_1[i] * wp_x[i];

        g_0_xyyy_0_xzzz_0[i] = g_0_yyy_0_zzz_1[i] * fi_abcd_0 + g_0_yyy_0_xzzz_0[i] * pb_x + g_0_yyy_0_xzzz_1[i] * wp_x[i];

        g_0_xyyy_0_yyyy_0[i] = g_0_yyy_0_yyyy_0[i] * pb_x + g_0_yyy_0_yyyy_1[i] * wp_x[i];

        g_0_xyyy_0_yyyz_0[i] = g_0_yyy_0_yyyz_0[i] * pb_x + g_0_yyy_0_yyyz_1[i] * wp_x[i];

        g_0_xyyy_0_yyzz_0[i] = g_0_yyy_0_yyzz_0[i] * pb_x + g_0_yyy_0_yyzz_1[i] * wp_x[i];

        g_0_xyyy_0_yzzz_0[i] = g_0_yyy_0_yzzz_0[i] * pb_x + g_0_yyy_0_yzzz_1[i] * wp_x[i];

        g_0_xyyy_0_zzzz_0[i] = g_0_yyy_0_zzzz_0[i] * pb_x + g_0_yyy_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 105-120 components of targeted buffer : SGSG

    auto g_0_xyyz_0_xxxx_0 = pbuffer.data(idx_eri_0_sgsg + 105);

    auto g_0_xyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_sgsg + 106);

    auto g_0_xyyz_0_xxxz_0 = pbuffer.data(idx_eri_0_sgsg + 107);

    auto g_0_xyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_sgsg + 108);

    auto g_0_xyyz_0_xxyz_0 = pbuffer.data(idx_eri_0_sgsg + 109);

    auto g_0_xyyz_0_xxzz_0 = pbuffer.data(idx_eri_0_sgsg + 110);

    auto g_0_xyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_sgsg + 111);

    auto g_0_xyyz_0_xyyz_0 = pbuffer.data(idx_eri_0_sgsg + 112);

    auto g_0_xyyz_0_xyzz_0 = pbuffer.data(idx_eri_0_sgsg + 113);

    auto g_0_xyyz_0_xzzz_0 = pbuffer.data(idx_eri_0_sgsg + 114);

    auto g_0_xyyz_0_yyyy_0 = pbuffer.data(idx_eri_0_sgsg + 115);

    auto g_0_xyyz_0_yyyz_0 = pbuffer.data(idx_eri_0_sgsg + 116);

    auto g_0_xyyz_0_yyzz_0 = pbuffer.data(idx_eri_0_sgsg + 117);

    auto g_0_xyyz_0_yzzz_0 = pbuffer.data(idx_eri_0_sgsg + 118);

    auto g_0_xyyz_0_zzzz_0 = pbuffer.data(idx_eri_0_sgsg + 119);

    #pragma omp simd aligned(g_0_xyy_0_xxxx_0, g_0_xyy_0_xxxx_1, g_0_xyy_0_xxxy_0, g_0_xyy_0_xxxy_1, g_0_xyy_0_xxyy_0, g_0_xyy_0_xxyy_1, g_0_xyy_0_xyyy_0, g_0_xyy_0_xyyy_1, g_0_xyyz_0_xxxx_0, g_0_xyyz_0_xxxy_0, g_0_xyyz_0_xxxz_0, g_0_xyyz_0_xxyy_0, g_0_xyyz_0_xxyz_0, g_0_xyyz_0_xxzz_0, g_0_xyyz_0_xyyy_0, g_0_xyyz_0_xyyz_0, g_0_xyyz_0_xyzz_0, g_0_xyyz_0_xzzz_0, g_0_xyyz_0_yyyy_0, g_0_xyyz_0_yyyz_0, g_0_xyyz_0_yyzz_0, g_0_xyyz_0_yzzz_0, g_0_xyyz_0_zzzz_0, g_0_yyz_0_xxxz_0, g_0_yyz_0_xxxz_1, g_0_yyz_0_xxyz_0, g_0_yyz_0_xxyz_1, g_0_yyz_0_xxz_1, g_0_yyz_0_xxzz_0, g_0_yyz_0_xxzz_1, g_0_yyz_0_xyyz_0, g_0_yyz_0_xyyz_1, g_0_yyz_0_xyz_1, g_0_yyz_0_xyzz_0, g_0_yyz_0_xyzz_1, g_0_yyz_0_xzz_1, g_0_yyz_0_xzzz_0, g_0_yyz_0_xzzz_1, g_0_yyz_0_yyyy_0, g_0_yyz_0_yyyy_1, g_0_yyz_0_yyyz_0, g_0_yyz_0_yyyz_1, g_0_yyz_0_yyz_1, g_0_yyz_0_yyzz_0, g_0_yyz_0_yyzz_1, g_0_yyz_0_yzz_1, g_0_yyz_0_yzzz_0, g_0_yyz_0_yzzz_1, g_0_yyz_0_zzz_1, g_0_yyz_0_zzzz_0, g_0_yyz_0_zzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyz_0_xxxx_0[i] = g_0_xyy_0_xxxx_0[i] * pb_z + g_0_xyy_0_xxxx_1[i] * wp_z[i];

        g_0_xyyz_0_xxxy_0[i] = g_0_xyy_0_xxxy_0[i] * pb_z + g_0_xyy_0_xxxy_1[i] * wp_z[i];

        g_0_xyyz_0_xxxz_0[i] = 3.0 * g_0_yyz_0_xxz_1[i] * fi_abcd_0 + g_0_yyz_0_xxxz_0[i] * pb_x + g_0_yyz_0_xxxz_1[i] * wp_x[i];

        g_0_xyyz_0_xxyy_0[i] = g_0_xyy_0_xxyy_0[i] * pb_z + g_0_xyy_0_xxyy_1[i] * wp_z[i];

        g_0_xyyz_0_xxyz_0[i] = 2.0 * g_0_yyz_0_xyz_1[i] * fi_abcd_0 + g_0_yyz_0_xxyz_0[i] * pb_x + g_0_yyz_0_xxyz_1[i] * wp_x[i];

        g_0_xyyz_0_xxzz_0[i] = 2.0 * g_0_yyz_0_xzz_1[i] * fi_abcd_0 + g_0_yyz_0_xxzz_0[i] * pb_x + g_0_yyz_0_xxzz_1[i] * wp_x[i];

        g_0_xyyz_0_xyyy_0[i] = g_0_xyy_0_xyyy_0[i] * pb_z + g_0_xyy_0_xyyy_1[i] * wp_z[i];

        g_0_xyyz_0_xyyz_0[i] = g_0_yyz_0_yyz_1[i] * fi_abcd_0 + g_0_yyz_0_xyyz_0[i] * pb_x + g_0_yyz_0_xyyz_1[i] * wp_x[i];

        g_0_xyyz_0_xyzz_0[i] = g_0_yyz_0_yzz_1[i] * fi_abcd_0 + g_0_yyz_0_xyzz_0[i] * pb_x + g_0_yyz_0_xyzz_1[i] * wp_x[i];

        g_0_xyyz_0_xzzz_0[i] = g_0_yyz_0_zzz_1[i] * fi_abcd_0 + g_0_yyz_0_xzzz_0[i] * pb_x + g_0_yyz_0_xzzz_1[i] * wp_x[i];

        g_0_xyyz_0_yyyy_0[i] = g_0_yyz_0_yyyy_0[i] * pb_x + g_0_yyz_0_yyyy_1[i] * wp_x[i];

        g_0_xyyz_0_yyyz_0[i] = g_0_yyz_0_yyyz_0[i] * pb_x + g_0_yyz_0_yyyz_1[i] * wp_x[i];

        g_0_xyyz_0_yyzz_0[i] = g_0_yyz_0_yyzz_0[i] * pb_x + g_0_yyz_0_yyzz_1[i] * wp_x[i];

        g_0_xyyz_0_yzzz_0[i] = g_0_yyz_0_yzzz_0[i] * pb_x + g_0_yyz_0_yzzz_1[i] * wp_x[i];

        g_0_xyyz_0_zzzz_0[i] = g_0_yyz_0_zzzz_0[i] * pb_x + g_0_yyz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 120-135 components of targeted buffer : SGSG

    auto g_0_xyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sgsg + 120);

    auto g_0_xyzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sgsg + 121);

    auto g_0_xyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sgsg + 122);

    auto g_0_xyzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sgsg + 123);

    auto g_0_xyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sgsg + 124);

    auto g_0_xyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sgsg + 125);

    auto g_0_xyzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sgsg + 126);

    auto g_0_xyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sgsg + 127);

    auto g_0_xyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sgsg + 128);

    auto g_0_xyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sgsg + 129);

    auto g_0_xyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sgsg + 130);

    auto g_0_xyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sgsg + 131);

    auto g_0_xyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sgsg + 132);

    auto g_0_xyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sgsg + 133);

    auto g_0_xyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sgsg + 134);

    #pragma omp simd aligned(g_0_xyzz_0_xxxx_0, g_0_xyzz_0_xxxy_0, g_0_xyzz_0_xxxz_0, g_0_xyzz_0_xxyy_0, g_0_xyzz_0_xxyz_0, g_0_xyzz_0_xxzz_0, g_0_xyzz_0_xyyy_0, g_0_xyzz_0_xyyz_0, g_0_xyzz_0_xyzz_0, g_0_xyzz_0_xzzz_0, g_0_xyzz_0_yyyy_0, g_0_xyzz_0_yyyz_0, g_0_xyzz_0_yyzz_0, g_0_xyzz_0_yzzz_0, g_0_xyzz_0_zzzz_0, g_0_xzz_0_xxxx_0, g_0_xzz_0_xxxx_1, g_0_xzz_0_xxxz_0, g_0_xzz_0_xxxz_1, g_0_xzz_0_xxzz_0, g_0_xzz_0_xxzz_1, g_0_xzz_0_xzzz_0, g_0_xzz_0_xzzz_1, g_0_yzz_0_xxxy_0, g_0_yzz_0_xxxy_1, g_0_yzz_0_xxy_1, g_0_yzz_0_xxyy_0, g_0_yzz_0_xxyy_1, g_0_yzz_0_xxyz_0, g_0_yzz_0_xxyz_1, g_0_yzz_0_xyy_1, g_0_yzz_0_xyyy_0, g_0_yzz_0_xyyy_1, g_0_yzz_0_xyyz_0, g_0_yzz_0_xyyz_1, g_0_yzz_0_xyz_1, g_0_yzz_0_xyzz_0, g_0_yzz_0_xyzz_1, g_0_yzz_0_yyy_1, g_0_yzz_0_yyyy_0, g_0_yzz_0_yyyy_1, g_0_yzz_0_yyyz_0, g_0_yzz_0_yyyz_1, g_0_yzz_0_yyz_1, g_0_yzz_0_yyzz_0, g_0_yzz_0_yyzz_1, g_0_yzz_0_yzz_1, g_0_yzz_0_yzzz_0, g_0_yzz_0_yzzz_1, g_0_yzz_0_zzzz_0, g_0_yzz_0_zzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyzz_0_xxxx_0[i] = g_0_xzz_0_xxxx_0[i] * pb_y + g_0_xzz_0_xxxx_1[i] * wp_y[i];

        g_0_xyzz_0_xxxy_0[i] = 3.0 * g_0_yzz_0_xxy_1[i] * fi_abcd_0 + g_0_yzz_0_xxxy_0[i] * pb_x + g_0_yzz_0_xxxy_1[i] * wp_x[i];

        g_0_xyzz_0_xxxz_0[i] = g_0_xzz_0_xxxz_0[i] * pb_y + g_0_xzz_0_xxxz_1[i] * wp_y[i];

        g_0_xyzz_0_xxyy_0[i] = 2.0 * g_0_yzz_0_xyy_1[i] * fi_abcd_0 + g_0_yzz_0_xxyy_0[i] * pb_x + g_0_yzz_0_xxyy_1[i] * wp_x[i];

        g_0_xyzz_0_xxyz_0[i] = 2.0 * g_0_yzz_0_xyz_1[i] * fi_abcd_0 + g_0_yzz_0_xxyz_0[i] * pb_x + g_0_yzz_0_xxyz_1[i] * wp_x[i];

        g_0_xyzz_0_xxzz_0[i] = g_0_xzz_0_xxzz_0[i] * pb_y + g_0_xzz_0_xxzz_1[i] * wp_y[i];

        g_0_xyzz_0_xyyy_0[i] = g_0_yzz_0_yyy_1[i] * fi_abcd_0 + g_0_yzz_0_xyyy_0[i] * pb_x + g_0_yzz_0_xyyy_1[i] * wp_x[i];

        g_0_xyzz_0_xyyz_0[i] = g_0_yzz_0_yyz_1[i] * fi_abcd_0 + g_0_yzz_0_xyyz_0[i] * pb_x + g_0_yzz_0_xyyz_1[i] * wp_x[i];

        g_0_xyzz_0_xyzz_0[i] = g_0_yzz_0_yzz_1[i] * fi_abcd_0 + g_0_yzz_0_xyzz_0[i] * pb_x + g_0_yzz_0_xyzz_1[i] * wp_x[i];

        g_0_xyzz_0_xzzz_0[i] = g_0_xzz_0_xzzz_0[i] * pb_y + g_0_xzz_0_xzzz_1[i] * wp_y[i];

        g_0_xyzz_0_yyyy_0[i] = g_0_yzz_0_yyyy_0[i] * pb_x + g_0_yzz_0_yyyy_1[i] * wp_x[i];

        g_0_xyzz_0_yyyz_0[i] = g_0_yzz_0_yyyz_0[i] * pb_x + g_0_yzz_0_yyyz_1[i] * wp_x[i];

        g_0_xyzz_0_yyzz_0[i] = g_0_yzz_0_yyzz_0[i] * pb_x + g_0_yzz_0_yyzz_1[i] * wp_x[i];

        g_0_xyzz_0_yzzz_0[i] = g_0_yzz_0_yzzz_0[i] * pb_x + g_0_yzz_0_yzzz_1[i] * wp_x[i];

        g_0_xyzz_0_zzzz_0[i] = g_0_yzz_0_zzzz_0[i] * pb_x + g_0_yzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 135-150 components of targeted buffer : SGSG

    auto g_0_xzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_sgsg + 135);

    auto g_0_xzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_sgsg + 136);

    auto g_0_xzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_sgsg + 137);

    auto g_0_xzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_sgsg + 138);

    auto g_0_xzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_sgsg + 139);

    auto g_0_xzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_sgsg + 140);

    auto g_0_xzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_sgsg + 141);

    auto g_0_xzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_sgsg + 142);

    auto g_0_xzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_sgsg + 143);

    auto g_0_xzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_sgsg + 144);

    auto g_0_xzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_sgsg + 145);

    auto g_0_xzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_sgsg + 146);

    auto g_0_xzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_sgsg + 147);

    auto g_0_xzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_sgsg + 148);

    auto g_0_xzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_sgsg + 149);

    #pragma omp simd aligned(g_0_xzzz_0_xxxx_0, g_0_xzzz_0_xxxy_0, g_0_xzzz_0_xxxz_0, g_0_xzzz_0_xxyy_0, g_0_xzzz_0_xxyz_0, g_0_xzzz_0_xxzz_0, g_0_xzzz_0_xyyy_0, g_0_xzzz_0_xyyz_0, g_0_xzzz_0_xyzz_0, g_0_xzzz_0_xzzz_0, g_0_xzzz_0_yyyy_0, g_0_xzzz_0_yyyz_0, g_0_xzzz_0_yyzz_0, g_0_xzzz_0_yzzz_0, g_0_xzzz_0_zzzz_0, g_0_zzz_0_xxx_1, g_0_zzz_0_xxxx_0, g_0_zzz_0_xxxx_1, g_0_zzz_0_xxxy_0, g_0_zzz_0_xxxy_1, g_0_zzz_0_xxxz_0, g_0_zzz_0_xxxz_1, g_0_zzz_0_xxy_1, g_0_zzz_0_xxyy_0, g_0_zzz_0_xxyy_1, g_0_zzz_0_xxyz_0, g_0_zzz_0_xxyz_1, g_0_zzz_0_xxz_1, g_0_zzz_0_xxzz_0, g_0_zzz_0_xxzz_1, g_0_zzz_0_xyy_1, g_0_zzz_0_xyyy_0, g_0_zzz_0_xyyy_1, g_0_zzz_0_xyyz_0, g_0_zzz_0_xyyz_1, g_0_zzz_0_xyz_1, g_0_zzz_0_xyzz_0, g_0_zzz_0_xyzz_1, g_0_zzz_0_xzz_1, g_0_zzz_0_xzzz_0, g_0_zzz_0_xzzz_1, g_0_zzz_0_yyy_1, g_0_zzz_0_yyyy_0, g_0_zzz_0_yyyy_1, g_0_zzz_0_yyyz_0, g_0_zzz_0_yyyz_1, g_0_zzz_0_yyz_1, g_0_zzz_0_yyzz_0, g_0_zzz_0_yyzz_1, g_0_zzz_0_yzz_1, g_0_zzz_0_yzzz_0, g_0_zzz_0_yzzz_1, g_0_zzz_0_zzz_1, g_0_zzz_0_zzzz_0, g_0_zzz_0_zzzz_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzzz_0_xxxx_0[i] = 4.0 * g_0_zzz_0_xxx_1[i] * fi_abcd_0 + g_0_zzz_0_xxxx_0[i] * pb_x + g_0_zzz_0_xxxx_1[i] * wp_x[i];

        g_0_xzzz_0_xxxy_0[i] = 3.0 * g_0_zzz_0_xxy_1[i] * fi_abcd_0 + g_0_zzz_0_xxxy_0[i] * pb_x + g_0_zzz_0_xxxy_1[i] * wp_x[i];

        g_0_xzzz_0_xxxz_0[i] = 3.0 * g_0_zzz_0_xxz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxz_0[i] * pb_x + g_0_zzz_0_xxxz_1[i] * wp_x[i];

        g_0_xzzz_0_xxyy_0[i] = 2.0 * g_0_zzz_0_xyy_1[i] * fi_abcd_0 + g_0_zzz_0_xxyy_0[i] * pb_x + g_0_zzz_0_xxyy_1[i] * wp_x[i];

        g_0_xzzz_0_xxyz_0[i] = 2.0 * g_0_zzz_0_xyz_1[i] * fi_abcd_0 + g_0_zzz_0_xxyz_0[i] * pb_x + g_0_zzz_0_xxyz_1[i] * wp_x[i];

        g_0_xzzz_0_xxzz_0[i] = 2.0 * g_0_zzz_0_xzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxzz_0[i] * pb_x + g_0_zzz_0_xxzz_1[i] * wp_x[i];

        g_0_xzzz_0_xyyy_0[i] = g_0_zzz_0_yyy_1[i] * fi_abcd_0 + g_0_zzz_0_xyyy_0[i] * pb_x + g_0_zzz_0_xyyy_1[i] * wp_x[i];

        g_0_xzzz_0_xyyz_0[i] = g_0_zzz_0_yyz_1[i] * fi_abcd_0 + g_0_zzz_0_xyyz_0[i] * pb_x + g_0_zzz_0_xyyz_1[i] * wp_x[i];

        g_0_xzzz_0_xyzz_0[i] = g_0_zzz_0_yzz_1[i] * fi_abcd_0 + g_0_zzz_0_xyzz_0[i] * pb_x + g_0_zzz_0_xyzz_1[i] * wp_x[i];

        g_0_xzzz_0_xzzz_0[i] = g_0_zzz_0_zzz_1[i] * fi_abcd_0 + g_0_zzz_0_xzzz_0[i] * pb_x + g_0_zzz_0_xzzz_1[i] * wp_x[i];

        g_0_xzzz_0_yyyy_0[i] = g_0_zzz_0_yyyy_0[i] * pb_x + g_0_zzz_0_yyyy_1[i] * wp_x[i];

        g_0_xzzz_0_yyyz_0[i] = g_0_zzz_0_yyyz_0[i] * pb_x + g_0_zzz_0_yyyz_1[i] * wp_x[i];

        g_0_xzzz_0_yyzz_0[i] = g_0_zzz_0_yyzz_0[i] * pb_x + g_0_zzz_0_yyzz_1[i] * wp_x[i];

        g_0_xzzz_0_yzzz_0[i] = g_0_zzz_0_yzzz_0[i] * pb_x + g_0_zzz_0_yzzz_1[i] * wp_x[i];

        g_0_xzzz_0_zzzz_0[i] = g_0_zzz_0_zzzz_0[i] * pb_x + g_0_zzz_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 150-165 components of targeted buffer : SGSG

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

    #pragma omp simd aligned(g_0_yy_0_xxxx_0, g_0_yy_0_xxxx_1, g_0_yy_0_xxxy_0, g_0_yy_0_xxxy_1, g_0_yy_0_xxxz_0, g_0_yy_0_xxxz_1, g_0_yy_0_xxyy_0, g_0_yy_0_xxyy_1, g_0_yy_0_xxyz_0, g_0_yy_0_xxyz_1, g_0_yy_0_xxzz_0, g_0_yy_0_xxzz_1, g_0_yy_0_xyyy_0, g_0_yy_0_xyyy_1, g_0_yy_0_xyyz_0, g_0_yy_0_xyyz_1, g_0_yy_0_xyzz_0, g_0_yy_0_xyzz_1, g_0_yy_0_xzzz_0, g_0_yy_0_xzzz_1, g_0_yy_0_yyyy_0, g_0_yy_0_yyyy_1, g_0_yy_0_yyyz_0, g_0_yy_0_yyyz_1, g_0_yy_0_yyzz_0, g_0_yy_0_yyzz_1, g_0_yy_0_yzzz_0, g_0_yy_0_yzzz_1, g_0_yy_0_zzzz_0, g_0_yy_0_zzzz_1, g_0_yyy_0_xxx_1, g_0_yyy_0_xxxx_0, g_0_yyy_0_xxxx_1, g_0_yyy_0_xxxy_0, g_0_yyy_0_xxxy_1, g_0_yyy_0_xxxz_0, g_0_yyy_0_xxxz_1, g_0_yyy_0_xxy_1, g_0_yyy_0_xxyy_0, g_0_yyy_0_xxyy_1, g_0_yyy_0_xxyz_0, g_0_yyy_0_xxyz_1, g_0_yyy_0_xxz_1, g_0_yyy_0_xxzz_0, g_0_yyy_0_xxzz_1, g_0_yyy_0_xyy_1, g_0_yyy_0_xyyy_0, g_0_yyy_0_xyyy_1, g_0_yyy_0_xyyz_0, g_0_yyy_0_xyyz_1, g_0_yyy_0_xyz_1, g_0_yyy_0_xyzz_0, g_0_yyy_0_xyzz_1, g_0_yyy_0_xzz_1, g_0_yyy_0_xzzz_0, g_0_yyy_0_xzzz_1, g_0_yyy_0_yyy_1, g_0_yyy_0_yyyy_0, g_0_yyy_0_yyyy_1, g_0_yyy_0_yyyz_0, g_0_yyy_0_yyyz_1, g_0_yyy_0_yyz_1, g_0_yyy_0_yyzz_0, g_0_yyy_0_yyzz_1, g_0_yyy_0_yzz_1, g_0_yyy_0_yzzz_0, g_0_yyy_0_yzzz_1, g_0_yyy_0_zzz_1, g_0_yyy_0_zzzz_0, g_0_yyy_0_zzzz_1, g_0_yyyy_0_xxxx_0, g_0_yyyy_0_xxxy_0, g_0_yyyy_0_xxxz_0, g_0_yyyy_0_xxyy_0, g_0_yyyy_0_xxyz_0, g_0_yyyy_0_xxzz_0, g_0_yyyy_0_xyyy_0, g_0_yyyy_0_xyyz_0, g_0_yyyy_0_xyzz_0, g_0_yyyy_0_xzzz_0, g_0_yyyy_0_yyyy_0, g_0_yyyy_0_yyyz_0, g_0_yyyy_0_yyzz_0, g_0_yyyy_0_yzzz_0, g_0_yyyy_0_zzzz_0, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyy_0_xxxx_0[i] = 3.0 * g_0_yy_0_xxxx_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxx_1[i] * fti_ab_0 + g_0_yyy_0_xxxx_0[i] * pb_y + g_0_yyy_0_xxxx_1[i] * wp_y[i];

        g_0_yyyy_0_xxxy_0[i] = 3.0 * g_0_yy_0_xxxy_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxy_1[i] * fti_ab_0 + g_0_yyy_0_xxx_1[i] * fi_abcd_0 + g_0_yyy_0_xxxy_0[i] * pb_y + g_0_yyy_0_xxxy_1[i] * wp_y[i];

        g_0_yyyy_0_xxxz_0[i] = 3.0 * g_0_yy_0_xxxz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxz_1[i] * fti_ab_0 + g_0_yyy_0_xxxz_0[i] * pb_y + g_0_yyy_0_xxxz_1[i] * wp_y[i];

        g_0_yyyy_0_xxyy_0[i] = 3.0 * g_0_yy_0_xxyy_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxyy_1[i] * fti_ab_0 + 2.0 * g_0_yyy_0_xxy_1[i] * fi_abcd_0 + g_0_yyy_0_xxyy_0[i] * pb_y + g_0_yyy_0_xxyy_1[i] * wp_y[i];

        g_0_yyyy_0_xxyz_0[i] = 3.0 * g_0_yy_0_xxyz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxyz_1[i] * fti_ab_0 + g_0_yyy_0_xxz_1[i] * fi_abcd_0 + g_0_yyy_0_xxyz_0[i] * pb_y + g_0_yyy_0_xxyz_1[i] * wp_y[i];

        g_0_yyyy_0_xxzz_0[i] = 3.0 * g_0_yy_0_xxzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxzz_1[i] * fti_ab_0 + g_0_yyy_0_xxzz_0[i] * pb_y + g_0_yyy_0_xxzz_1[i] * wp_y[i];

        g_0_yyyy_0_xyyy_0[i] = 3.0 * g_0_yy_0_xyyy_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xyyy_1[i] * fti_ab_0 + 3.0 * g_0_yyy_0_xyy_1[i] * fi_abcd_0 + g_0_yyy_0_xyyy_0[i] * pb_y + g_0_yyy_0_xyyy_1[i] * wp_y[i];

        g_0_yyyy_0_xyyz_0[i] = 3.0 * g_0_yy_0_xyyz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xyyz_1[i] * fti_ab_0 + 2.0 * g_0_yyy_0_xyz_1[i] * fi_abcd_0 + g_0_yyy_0_xyyz_0[i] * pb_y + g_0_yyy_0_xyyz_1[i] * wp_y[i];

        g_0_yyyy_0_xyzz_0[i] = 3.0 * g_0_yy_0_xyzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xyzz_1[i] * fti_ab_0 + g_0_yyy_0_xzz_1[i] * fi_abcd_0 + g_0_yyy_0_xyzz_0[i] * pb_y + g_0_yyy_0_xyzz_1[i] * wp_y[i];

        g_0_yyyy_0_xzzz_0[i] = 3.0 * g_0_yy_0_xzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xzzz_1[i] * fti_ab_0 + g_0_yyy_0_xzzz_0[i] * pb_y + g_0_yyy_0_xzzz_1[i] * wp_y[i];

        g_0_yyyy_0_yyyy_0[i] = 3.0 * g_0_yy_0_yyyy_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_yyyy_1[i] * fti_ab_0 + 4.0 * g_0_yyy_0_yyy_1[i] * fi_abcd_0 + g_0_yyy_0_yyyy_0[i] * pb_y + g_0_yyy_0_yyyy_1[i] * wp_y[i];

        g_0_yyyy_0_yyyz_0[i] = 3.0 * g_0_yy_0_yyyz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_yyyz_1[i] * fti_ab_0 + 3.0 * g_0_yyy_0_yyz_1[i] * fi_abcd_0 + g_0_yyy_0_yyyz_0[i] * pb_y + g_0_yyy_0_yyyz_1[i] * wp_y[i];

        g_0_yyyy_0_yyzz_0[i] = 3.0 * g_0_yy_0_yyzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_yyzz_1[i] * fti_ab_0 + 2.0 * g_0_yyy_0_yzz_1[i] * fi_abcd_0 + g_0_yyy_0_yyzz_0[i] * pb_y + g_0_yyy_0_yyzz_1[i] * wp_y[i];

        g_0_yyyy_0_yzzz_0[i] = 3.0 * g_0_yy_0_yzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_yzzz_1[i] * fti_ab_0 + g_0_yyy_0_zzz_1[i] * fi_abcd_0 + g_0_yyy_0_yzzz_0[i] * pb_y + g_0_yyy_0_yzzz_1[i] * wp_y[i];

        g_0_yyyy_0_zzzz_0[i] = 3.0 * g_0_yy_0_zzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_zzzz_1[i] * fti_ab_0 + g_0_yyy_0_zzzz_0[i] * pb_y + g_0_yyy_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 165-180 components of targeted buffer : SGSG

    auto g_0_yyyz_0_xxxx_0 = pbuffer.data(idx_eri_0_sgsg + 165);

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

    #pragma omp simd aligned(g_0_yyy_0_xxx_1, g_0_yyy_0_xxxx_0, g_0_yyy_0_xxxx_1, g_0_yyy_0_xxxy_0, g_0_yyy_0_xxxy_1, g_0_yyy_0_xxxz_0, g_0_yyy_0_xxxz_1, g_0_yyy_0_xxy_1, g_0_yyy_0_xxyy_0, g_0_yyy_0_xxyy_1, g_0_yyy_0_xxyz_0, g_0_yyy_0_xxyz_1, g_0_yyy_0_xxz_1, g_0_yyy_0_xxzz_0, g_0_yyy_0_xxzz_1, g_0_yyy_0_xyy_1, g_0_yyy_0_xyyy_0, g_0_yyy_0_xyyy_1, g_0_yyy_0_xyyz_0, g_0_yyy_0_xyyz_1, g_0_yyy_0_xyz_1, g_0_yyy_0_xyzz_0, g_0_yyy_0_xyzz_1, g_0_yyy_0_xzz_1, g_0_yyy_0_xzzz_0, g_0_yyy_0_xzzz_1, g_0_yyy_0_yyy_1, g_0_yyy_0_yyyy_0, g_0_yyy_0_yyyy_1, g_0_yyy_0_yyyz_0, g_0_yyy_0_yyyz_1, g_0_yyy_0_yyz_1, g_0_yyy_0_yyzz_0, g_0_yyy_0_yyzz_1, g_0_yyy_0_yzz_1, g_0_yyy_0_yzzz_0, g_0_yyy_0_yzzz_1, g_0_yyy_0_zzz_1, g_0_yyy_0_zzzz_0, g_0_yyy_0_zzzz_1, g_0_yyyz_0_xxxx_0, g_0_yyyz_0_xxxy_0, g_0_yyyz_0_xxxz_0, g_0_yyyz_0_xxyy_0, g_0_yyyz_0_xxyz_0, g_0_yyyz_0_xxzz_0, g_0_yyyz_0_xyyy_0, g_0_yyyz_0_xyyz_0, g_0_yyyz_0_xyzz_0, g_0_yyyz_0_xzzz_0, g_0_yyyz_0_yyyy_0, g_0_yyyz_0_yyyz_0, g_0_yyyz_0_yyzz_0, g_0_yyyz_0_yzzz_0, g_0_yyyz_0_zzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyz_0_xxxx_0[i] = g_0_yyy_0_xxxx_0[i] * pb_z + g_0_yyy_0_xxxx_1[i] * wp_z[i];

        g_0_yyyz_0_xxxy_0[i] = g_0_yyy_0_xxxy_0[i] * pb_z + g_0_yyy_0_xxxy_1[i] * wp_z[i];

        g_0_yyyz_0_xxxz_0[i] = g_0_yyy_0_xxx_1[i] * fi_abcd_0 + g_0_yyy_0_xxxz_0[i] * pb_z + g_0_yyy_0_xxxz_1[i] * wp_z[i];

        g_0_yyyz_0_xxyy_0[i] = g_0_yyy_0_xxyy_0[i] * pb_z + g_0_yyy_0_xxyy_1[i] * wp_z[i];

        g_0_yyyz_0_xxyz_0[i] = g_0_yyy_0_xxy_1[i] * fi_abcd_0 + g_0_yyy_0_xxyz_0[i] * pb_z + g_0_yyy_0_xxyz_1[i] * wp_z[i];

        g_0_yyyz_0_xxzz_0[i] = 2.0 * g_0_yyy_0_xxz_1[i] * fi_abcd_0 + g_0_yyy_0_xxzz_0[i] * pb_z + g_0_yyy_0_xxzz_1[i] * wp_z[i];

        g_0_yyyz_0_xyyy_0[i] = g_0_yyy_0_xyyy_0[i] * pb_z + g_0_yyy_0_xyyy_1[i] * wp_z[i];

        g_0_yyyz_0_xyyz_0[i] = g_0_yyy_0_xyy_1[i] * fi_abcd_0 + g_0_yyy_0_xyyz_0[i] * pb_z + g_0_yyy_0_xyyz_1[i] * wp_z[i];

        g_0_yyyz_0_xyzz_0[i] = 2.0 * g_0_yyy_0_xyz_1[i] * fi_abcd_0 + g_0_yyy_0_xyzz_0[i] * pb_z + g_0_yyy_0_xyzz_1[i] * wp_z[i];

        g_0_yyyz_0_xzzz_0[i] = 3.0 * g_0_yyy_0_xzz_1[i] * fi_abcd_0 + g_0_yyy_0_xzzz_0[i] * pb_z + g_0_yyy_0_xzzz_1[i] * wp_z[i];

        g_0_yyyz_0_yyyy_0[i] = g_0_yyy_0_yyyy_0[i] * pb_z + g_0_yyy_0_yyyy_1[i] * wp_z[i];

        g_0_yyyz_0_yyyz_0[i] = g_0_yyy_0_yyy_1[i] * fi_abcd_0 + g_0_yyy_0_yyyz_0[i] * pb_z + g_0_yyy_0_yyyz_1[i] * wp_z[i];

        g_0_yyyz_0_yyzz_0[i] = 2.0 * g_0_yyy_0_yyz_1[i] * fi_abcd_0 + g_0_yyy_0_yyzz_0[i] * pb_z + g_0_yyy_0_yyzz_1[i] * wp_z[i];

        g_0_yyyz_0_yzzz_0[i] = 3.0 * g_0_yyy_0_yzz_1[i] * fi_abcd_0 + g_0_yyy_0_yzzz_0[i] * pb_z + g_0_yyy_0_yzzz_1[i] * wp_z[i];

        g_0_yyyz_0_zzzz_0[i] = 4.0 * g_0_yyy_0_zzz_1[i] * fi_abcd_0 + g_0_yyy_0_zzzz_0[i] * pb_z + g_0_yyy_0_zzzz_1[i] * wp_z[i];
    }

    /// Set up 180-195 components of targeted buffer : SGSG

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

    #pragma omp simd aligned(g_0_yy_0_xxxy_0, g_0_yy_0_xxxy_1, g_0_yy_0_xxyy_0, g_0_yy_0_xxyy_1, g_0_yy_0_xyyy_0, g_0_yy_0_xyyy_1, g_0_yy_0_yyyy_0, g_0_yy_0_yyyy_1, g_0_yyz_0_xxxy_0, g_0_yyz_0_xxxy_1, g_0_yyz_0_xxyy_0, g_0_yyz_0_xxyy_1, g_0_yyz_0_xyyy_0, g_0_yyz_0_xyyy_1, g_0_yyz_0_yyyy_0, g_0_yyz_0_yyyy_1, g_0_yyzz_0_xxxx_0, g_0_yyzz_0_xxxy_0, g_0_yyzz_0_xxxz_0, g_0_yyzz_0_xxyy_0, g_0_yyzz_0_xxyz_0, g_0_yyzz_0_xxzz_0, g_0_yyzz_0_xyyy_0, g_0_yyzz_0_xyyz_0, g_0_yyzz_0_xyzz_0, g_0_yyzz_0_xzzz_0, g_0_yyzz_0_yyyy_0, g_0_yyzz_0_yyyz_0, g_0_yyzz_0_yyzz_0, g_0_yyzz_0_yzzz_0, g_0_yyzz_0_zzzz_0, g_0_yzz_0_xxxx_0, g_0_yzz_0_xxxx_1, g_0_yzz_0_xxxz_0, g_0_yzz_0_xxxz_1, g_0_yzz_0_xxyz_0, g_0_yzz_0_xxyz_1, g_0_yzz_0_xxz_1, g_0_yzz_0_xxzz_0, g_0_yzz_0_xxzz_1, g_0_yzz_0_xyyz_0, g_0_yzz_0_xyyz_1, g_0_yzz_0_xyz_1, g_0_yzz_0_xyzz_0, g_0_yzz_0_xyzz_1, g_0_yzz_0_xzz_1, g_0_yzz_0_xzzz_0, g_0_yzz_0_xzzz_1, g_0_yzz_0_yyyz_0, g_0_yzz_0_yyyz_1, g_0_yzz_0_yyz_1, g_0_yzz_0_yyzz_0, g_0_yzz_0_yyzz_1, g_0_yzz_0_yzz_1, g_0_yzz_0_yzzz_0, g_0_yzz_0_yzzz_1, g_0_yzz_0_zzz_1, g_0_yzz_0_zzzz_0, g_0_yzz_0_zzzz_1, g_0_zz_0_xxxx_0, g_0_zz_0_xxxx_1, g_0_zz_0_xxxz_0, g_0_zz_0_xxxz_1, g_0_zz_0_xxyz_0, g_0_zz_0_xxyz_1, g_0_zz_0_xxzz_0, g_0_zz_0_xxzz_1, g_0_zz_0_xyyz_0, g_0_zz_0_xyyz_1, g_0_zz_0_xyzz_0, g_0_zz_0_xyzz_1, g_0_zz_0_xzzz_0, g_0_zz_0_xzzz_1, g_0_zz_0_yyyz_0, g_0_zz_0_yyyz_1, g_0_zz_0_yyzz_0, g_0_zz_0_yyzz_1, g_0_zz_0_yzzz_0, g_0_zz_0_yzzz_1, g_0_zz_0_zzzz_0, g_0_zz_0_zzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyzz_0_xxxx_0[i] = g_0_zz_0_xxxx_0[i] * fi_ab_0 - g_0_zz_0_xxxx_1[i] * fti_ab_0 + g_0_yzz_0_xxxx_0[i] * pb_y + g_0_yzz_0_xxxx_1[i] * wp_y[i];

        g_0_yyzz_0_xxxy_0[i] = g_0_yy_0_xxxy_0[i] * fi_ab_0 - g_0_yy_0_xxxy_1[i] * fti_ab_0 + g_0_yyz_0_xxxy_0[i] * pb_z + g_0_yyz_0_xxxy_1[i] * wp_z[i];

        g_0_yyzz_0_xxxz_0[i] = g_0_zz_0_xxxz_0[i] * fi_ab_0 - g_0_zz_0_xxxz_1[i] * fti_ab_0 + g_0_yzz_0_xxxz_0[i] * pb_y + g_0_yzz_0_xxxz_1[i] * wp_y[i];

        g_0_yyzz_0_xxyy_0[i] = g_0_yy_0_xxyy_0[i] * fi_ab_0 - g_0_yy_0_xxyy_1[i] * fti_ab_0 + g_0_yyz_0_xxyy_0[i] * pb_z + g_0_yyz_0_xxyy_1[i] * wp_z[i];

        g_0_yyzz_0_xxyz_0[i] = g_0_zz_0_xxyz_0[i] * fi_ab_0 - g_0_zz_0_xxyz_1[i] * fti_ab_0 + g_0_yzz_0_xxz_1[i] * fi_abcd_0 + g_0_yzz_0_xxyz_0[i] * pb_y + g_0_yzz_0_xxyz_1[i] * wp_y[i];

        g_0_yyzz_0_xxzz_0[i] = g_0_zz_0_xxzz_0[i] * fi_ab_0 - g_0_zz_0_xxzz_1[i] * fti_ab_0 + g_0_yzz_0_xxzz_0[i] * pb_y + g_0_yzz_0_xxzz_1[i] * wp_y[i];

        g_0_yyzz_0_xyyy_0[i] = g_0_yy_0_xyyy_0[i] * fi_ab_0 - g_0_yy_0_xyyy_1[i] * fti_ab_0 + g_0_yyz_0_xyyy_0[i] * pb_z + g_0_yyz_0_xyyy_1[i] * wp_z[i];

        g_0_yyzz_0_xyyz_0[i] = g_0_zz_0_xyyz_0[i] * fi_ab_0 - g_0_zz_0_xyyz_1[i] * fti_ab_0 + 2.0 * g_0_yzz_0_xyz_1[i] * fi_abcd_0 + g_0_yzz_0_xyyz_0[i] * pb_y + g_0_yzz_0_xyyz_1[i] * wp_y[i];

        g_0_yyzz_0_xyzz_0[i] = g_0_zz_0_xyzz_0[i] * fi_ab_0 - g_0_zz_0_xyzz_1[i] * fti_ab_0 + g_0_yzz_0_xzz_1[i] * fi_abcd_0 + g_0_yzz_0_xyzz_0[i] * pb_y + g_0_yzz_0_xyzz_1[i] * wp_y[i];

        g_0_yyzz_0_xzzz_0[i] = g_0_zz_0_xzzz_0[i] * fi_ab_0 - g_0_zz_0_xzzz_1[i] * fti_ab_0 + g_0_yzz_0_xzzz_0[i] * pb_y + g_0_yzz_0_xzzz_1[i] * wp_y[i];

        g_0_yyzz_0_yyyy_0[i] = g_0_yy_0_yyyy_0[i] * fi_ab_0 - g_0_yy_0_yyyy_1[i] * fti_ab_0 + g_0_yyz_0_yyyy_0[i] * pb_z + g_0_yyz_0_yyyy_1[i] * wp_z[i];

        g_0_yyzz_0_yyyz_0[i] = g_0_zz_0_yyyz_0[i] * fi_ab_0 - g_0_zz_0_yyyz_1[i] * fti_ab_0 + 3.0 * g_0_yzz_0_yyz_1[i] * fi_abcd_0 + g_0_yzz_0_yyyz_0[i] * pb_y + g_0_yzz_0_yyyz_1[i] * wp_y[i];

        g_0_yyzz_0_yyzz_0[i] = g_0_zz_0_yyzz_0[i] * fi_ab_0 - g_0_zz_0_yyzz_1[i] * fti_ab_0 + 2.0 * g_0_yzz_0_yzz_1[i] * fi_abcd_0 + g_0_yzz_0_yyzz_0[i] * pb_y + g_0_yzz_0_yyzz_1[i] * wp_y[i];

        g_0_yyzz_0_yzzz_0[i] = g_0_zz_0_yzzz_0[i] * fi_ab_0 - g_0_zz_0_yzzz_1[i] * fti_ab_0 + g_0_yzz_0_zzz_1[i] * fi_abcd_0 + g_0_yzz_0_yzzz_0[i] * pb_y + g_0_yzz_0_yzzz_1[i] * wp_y[i];

        g_0_yyzz_0_zzzz_0[i] = g_0_zz_0_zzzz_0[i] * fi_ab_0 - g_0_zz_0_zzzz_1[i] * fti_ab_0 + g_0_yzz_0_zzzz_0[i] * pb_y + g_0_yzz_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 195-210 components of targeted buffer : SGSG

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

    #pragma omp simd aligned(g_0_yzzz_0_xxxx_0, g_0_yzzz_0_xxxy_0, g_0_yzzz_0_xxxz_0, g_0_yzzz_0_xxyy_0, g_0_yzzz_0_xxyz_0, g_0_yzzz_0_xxzz_0, g_0_yzzz_0_xyyy_0, g_0_yzzz_0_xyyz_0, g_0_yzzz_0_xyzz_0, g_0_yzzz_0_xzzz_0, g_0_yzzz_0_yyyy_0, g_0_yzzz_0_yyyz_0, g_0_yzzz_0_yyzz_0, g_0_yzzz_0_yzzz_0, g_0_yzzz_0_zzzz_0, g_0_zzz_0_xxx_1, g_0_zzz_0_xxxx_0, g_0_zzz_0_xxxx_1, g_0_zzz_0_xxxy_0, g_0_zzz_0_xxxy_1, g_0_zzz_0_xxxz_0, g_0_zzz_0_xxxz_1, g_0_zzz_0_xxy_1, g_0_zzz_0_xxyy_0, g_0_zzz_0_xxyy_1, g_0_zzz_0_xxyz_0, g_0_zzz_0_xxyz_1, g_0_zzz_0_xxz_1, g_0_zzz_0_xxzz_0, g_0_zzz_0_xxzz_1, g_0_zzz_0_xyy_1, g_0_zzz_0_xyyy_0, g_0_zzz_0_xyyy_1, g_0_zzz_0_xyyz_0, g_0_zzz_0_xyyz_1, g_0_zzz_0_xyz_1, g_0_zzz_0_xyzz_0, g_0_zzz_0_xyzz_1, g_0_zzz_0_xzz_1, g_0_zzz_0_xzzz_0, g_0_zzz_0_xzzz_1, g_0_zzz_0_yyy_1, g_0_zzz_0_yyyy_0, g_0_zzz_0_yyyy_1, g_0_zzz_0_yyyz_0, g_0_zzz_0_yyyz_1, g_0_zzz_0_yyz_1, g_0_zzz_0_yyzz_0, g_0_zzz_0_yyzz_1, g_0_zzz_0_yzz_1, g_0_zzz_0_yzzz_0, g_0_zzz_0_yzzz_1, g_0_zzz_0_zzz_1, g_0_zzz_0_zzzz_0, g_0_zzz_0_zzzz_1, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzzz_0_xxxx_0[i] = g_0_zzz_0_xxxx_0[i] * pb_y + g_0_zzz_0_xxxx_1[i] * wp_y[i];

        g_0_yzzz_0_xxxy_0[i] = g_0_zzz_0_xxx_1[i] * fi_abcd_0 + g_0_zzz_0_xxxy_0[i] * pb_y + g_0_zzz_0_xxxy_1[i] * wp_y[i];

        g_0_yzzz_0_xxxz_0[i] = g_0_zzz_0_xxxz_0[i] * pb_y + g_0_zzz_0_xxxz_1[i] * wp_y[i];

        g_0_yzzz_0_xxyy_0[i] = 2.0 * g_0_zzz_0_xxy_1[i] * fi_abcd_0 + g_0_zzz_0_xxyy_0[i] * pb_y + g_0_zzz_0_xxyy_1[i] * wp_y[i];

        g_0_yzzz_0_xxyz_0[i] = g_0_zzz_0_xxz_1[i] * fi_abcd_0 + g_0_zzz_0_xxyz_0[i] * pb_y + g_0_zzz_0_xxyz_1[i] * wp_y[i];

        g_0_yzzz_0_xxzz_0[i] = g_0_zzz_0_xxzz_0[i] * pb_y + g_0_zzz_0_xxzz_1[i] * wp_y[i];

        g_0_yzzz_0_xyyy_0[i] = 3.0 * g_0_zzz_0_xyy_1[i] * fi_abcd_0 + g_0_zzz_0_xyyy_0[i] * pb_y + g_0_zzz_0_xyyy_1[i] * wp_y[i];

        g_0_yzzz_0_xyyz_0[i] = 2.0 * g_0_zzz_0_xyz_1[i] * fi_abcd_0 + g_0_zzz_0_xyyz_0[i] * pb_y + g_0_zzz_0_xyyz_1[i] * wp_y[i];

        g_0_yzzz_0_xyzz_0[i] = g_0_zzz_0_xzz_1[i] * fi_abcd_0 + g_0_zzz_0_xyzz_0[i] * pb_y + g_0_zzz_0_xyzz_1[i] * wp_y[i];

        g_0_yzzz_0_xzzz_0[i] = g_0_zzz_0_xzzz_0[i] * pb_y + g_0_zzz_0_xzzz_1[i] * wp_y[i];

        g_0_yzzz_0_yyyy_0[i] = 4.0 * g_0_zzz_0_yyy_1[i] * fi_abcd_0 + g_0_zzz_0_yyyy_0[i] * pb_y + g_0_zzz_0_yyyy_1[i] * wp_y[i];

        g_0_yzzz_0_yyyz_0[i] = 3.0 * g_0_zzz_0_yyz_1[i] * fi_abcd_0 + g_0_zzz_0_yyyz_0[i] * pb_y + g_0_zzz_0_yyyz_1[i] * wp_y[i];

        g_0_yzzz_0_yyzz_0[i] = 2.0 * g_0_zzz_0_yzz_1[i] * fi_abcd_0 + g_0_zzz_0_yyzz_0[i] * pb_y + g_0_zzz_0_yyzz_1[i] * wp_y[i];

        g_0_yzzz_0_yzzz_0[i] = g_0_zzz_0_zzz_1[i] * fi_abcd_0 + g_0_zzz_0_yzzz_0[i] * pb_y + g_0_zzz_0_yzzz_1[i] * wp_y[i];

        g_0_yzzz_0_zzzz_0[i] = g_0_zzz_0_zzzz_0[i] * pb_y + g_0_zzz_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 210-225 components of targeted buffer : SGSG

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

    #pragma omp simd aligned(g_0_zz_0_xxxx_0, g_0_zz_0_xxxx_1, g_0_zz_0_xxxy_0, g_0_zz_0_xxxy_1, g_0_zz_0_xxxz_0, g_0_zz_0_xxxz_1, g_0_zz_0_xxyy_0, g_0_zz_0_xxyy_1, g_0_zz_0_xxyz_0, g_0_zz_0_xxyz_1, g_0_zz_0_xxzz_0, g_0_zz_0_xxzz_1, g_0_zz_0_xyyy_0, g_0_zz_0_xyyy_1, g_0_zz_0_xyyz_0, g_0_zz_0_xyyz_1, g_0_zz_0_xyzz_0, g_0_zz_0_xyzz_1, g_0_zz_0_xzzz_0, g_0_zz_0_xzzz_1, g_0_zz_0_yyyy_0, g_0_zz_0_yyyy_1, g_0_zz_0_yyyz_0, g_0_zz_0_yyyz_1, g_0_zz_0_yyzz_0, g_0_zz_0_yyzz_1, g_0_zz_0_yzzz_0, g_0_zz_0_yzzz_1, g_0_zz_0_zzzz_0, g_0_zz_0_zzzz_1, g_0_zzz_0_xxx_1, g_0_zzz_0_xxxx_0, g_0_zzz_0_xxxx_1, g_0_zzz_0_xxxy_0, g_0_zzz_0_xxxy_1, g_0_zzz_0_xxxz_0, g_0_zzz_0_xxxz_1, g_0_zzz_0_xxy_1, g_0_zzz_0_xxyy_0, g_0_zzz_0_xxyy_1, g_0_zzz_0_xxyz_0, g_0_zzz_0_xxyz_1, g_0_zzz_0_xxz_1, g_0_zzz_0_xxzz_0, g_0_zzz_0_xxzz_1, g_0_zzz_0_xyy_1, g_0_zzz_0_xyyy_0, g_0_zzz_0_xyyy_1, g_0_zzz_0_xyyz_0, g_0_zzz_0_xyyz_1, g_0_zzz_0_xyz_1, g_0_zzz_0_xyzz_0, g_0_zzz_0_xyzz_1, g_0_zzz_0_xzz_1, g_0_zzz_0_xzzz_0, g_0_zzz_0_xzzz_1, g_0_zzz_0_yyy_1, g_0_zzz_0_yyyy_0, g_0_zzz_0_yyyy_1, g_0_zzz_0_yyyz_0, g_0_zzz_0_yyyz_1, g_0_zzz_0_yyz_1, g_0_zzz_0_yyzz_0, g_0_zzz_0_yyzz_1, g_0_zzz_0_yzz_1, g_0_zzz_0_yzzz_0, g_0_zzz_0_yzzz_1, g_0_zzz_0_zzz_1, g_0_zzz_0_zzzz_0, g_0_zzz_0_zzzz_1, g_0_zzzz_0_xxxx_0, g_0_zzzz_0_xxxy_0, g_0_zzzz_0_xxxz_0, g_0_zzzz_0_xxyy_0, g_0_zzzz_0_xxyz_0, g_0_zzzz_0_xxzz_0, g_0_zzzz_0_xyyy_0, g_0_zzzz_0_xyyz_0, g_0_zzzz_0_xyzz_0, g_0_zzzz_0_xzzz_0, g_0_zzzz_0_yyyy_0, g_0_zzzz_0_yyyz_0, g_0_zzzz_0_yyzz_0, g_0_zzzz_0_yzzz_0, g_0_zzzz_0_zzzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzzz_0_xxxx_0[i] = 3.0 * g_0_zz_0_xxxx_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxx_1[i] * fti_ab_0 + g_0_zzz_0_xxxx_0[i] * pb_z + g_0_zzz_0_xxxx_1[i] * wp_z[i];

        g_0_zzzz_0_xxxy_0[i] = 3.0 * g_0_zz_0_xxxy_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxy_1[i] * fti_ab_0 + g_0_zzz_0_xxxy_0[i] * pb_z + g_0_zzz_0_xxxy_1[i] * wp_z[i];

        g_0_zzzz_0_xxxz_0[i] = 3.0 * g_0_zz_0_xxxz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxz_1[i] * fti_ab_0 + g_0_zzz_0_xxx_1[i] * fi_abcd_0 + g_0_zzz_0_xxxz_0[i] * pb_z + g_0_zzz_0_xxxz_1[i] * wp_z[i];

        g_0_zzzz_0_xxyy_0[i] = 3.0 * g_0_zz_0_xxyy_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxyy_1[i] * fti_ab_0 + g_0_zzz_0_xxyy_0[i] * pb_z + g_0_zzz_0_xxyy_1[i] * wp_z[i];

        g_0_zzzz_0_xxyz_0[i] = 3.0 * g_0_zz_0_xxyz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxyz_1[i] * fti_ab_0 + g_0_zzz_0_xxy_1[i] * fi_abcd_0 + g_0_zzz_0_xxyz_0[i] * pb_z + g_0_zzz_0_xxyz_1[i] * wp_z[i];

        g_0_zzzz_0_xxzz_0[i] = 3.0 * g_0_zz_0_xxzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxzz_1[i] * fti_ab_0 + 2.0 * g_0_zzz_0_xxz_1[i] * fi_abcd_0 + g_0_zzz_0_xxzz_0[i] * pb_z + g_0_zzz_0_xxzz_1[i] * wp_z[i];

        g_0_zzzz_0_xyyy_0[i] = 3.0 * g_0_zz_0_xyyy_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xyyy_1[i] * fti_ab_0 + g_0_zzz_0_xyyy_0[i] * pb_z + g_0_zzz_0_xyyy_1[i] * wp_z[i];

        g_0_zzzz_0_xyyz_0[i] = 3.0 * g_0_zz_0_xyyz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xyyz_1[i] * fti_ab_0 + g_0_zzz_0_xyy_1[i] * fi_abcd_0 + g_0_zzz_0_xyyz_0[i] * pb_z + g_0_zzz_0_xyyz_1[i] * wp_z[i];

        g_0_zzzz_0_xyzz_0[i] = 3.0 * g_0_zz_0_xyzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzz_0_xyz_1[i] * fi_abcd_0 + g_0_zzz_0_xyzz_0[i] * pb_z + g_0_zzz_0_xyzz_1[i] * wp_z[i];

        g_0_zzzz_0_xzzz_0[i] = 3.0 * g_0_zz_0_xzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzz_0_xzz_1[i] * fi_abcd_0 + g_0_zzz_0_xzzz_0[i] * pb_z + g_0_zzz_0_xzzz_1[i] * wp_z[i];

        g_0_zzzz_0_yyyy_0[i] = 3.0 * g_0_zz_0_yyyy_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_yyyy_1[i] * fti_ab_0 + g_0_zzz_0_yyyy_0[i] * pb_z + g_0_zzz_0_yyyy_1[i] * wp_z[i];

        g_0_zzzz_0_yyyz_0[i] = 3.0 * g_0_zz_0_yyyz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_yyyz_1[i] * fti_ab_0 + g_0_zzz_0_yyy_1[i] * fi_abcd_0 + g_0_zzz_0_yyyz_0[i] * pb_z + g_0_zzz_0_yyyz_1[i] * wp_z[i];

        g_0_zzzz_0_yyzz_0[i] = 3.0 * g_0_zz_0_yyzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_yyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzz_0_yyz_1[i] * fi_abcd_0 + g_0_zzz_0_yyzz_0[i] * pb_z + g_0_zzz_0_yyzz_1[i] * wp_z[i];

        g_0_zzzz_0_yzzz_0[i] = 3.0 * g_0_zz_0_yzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_yzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzz_0_yzz_1[i] * fi_abcd_0 + g_0_zzz_0_yzzz_0[i] * pb_z + g_0_zzz_0_yzzz_1[i] * wp_z[i];

        g_0_zzzz_0_zzzz_0[i] = 3.0 * g_0_zz_0_zzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_zzzz_1[i] * fti_ab_0 + 4.0 * g_0_zzz_0_zzz_1[i] * fi_abcd_0 + g_0_zzz_0_zzzz_0[i] * pb_z + g_0_zzz_0_zzzz_1[i] * wp_z[i];
    }
}

} // erirec namespace


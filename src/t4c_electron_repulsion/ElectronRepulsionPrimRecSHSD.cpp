#include "ElectronRepulsionPrimRecSHSD.hpp"

namespace erirec {  // erirec namespace

auto
comp_prim_electron_repulsion_shsd(CSimdArray<double>&   pbuffer,
                                  const size_t          idx_eri_0_shsd,
                                  size_t                idx_eri_0_sfsd,
                                  size_t                idx_eri_1_sfsd,
                                  size_t                idx_eri_1_sgsp,
                                  size_t                idx_eri_0_sgsd,
                                  size_t                idx_eri_1_sgsd,
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

    /// Set up components of auxilary buffer : SFSD

    auto g_0_xxx_0_xx_0 = pbuffer.data(idx_eri_0_sfsd);

    auto g_0_xxx_0_xy_0 = pbuffer.data(idx_eri_0_sfsd + 1);

    auto g_0_xxx_0_xz_0 = pbuffer.data(idx_eri_0_sfsd + 2);

    auto g_0_xxx_0_yy_0 = pbuffer.data(idx_eri_0_sfsd + 3);

    auto g_0_xxx_0_yz_0 = pbuffer.data(idx_eri_0_sfsd + 4);

    auto g_0_xxx_0_zz_0 = pbuffer.data(idx_eri_0_sfsd + 5);

    auto g_0_xxy_0_xx_0 = pbuffer.data(idx_eri_0_sfsd + 6);

    auto g_0_xxy_0_xz_0 = pbuffer.data(idx_eri_0_sfsd + 8);

    auto g_0_xxz_0_xx_0 = pbuffer.data(idx_eri_0_sfsd + 12);

    auto g_0_xxz_0_xy_0 = pbuffer.data(idx_eri_0_sfsd + 13);

    auto g_0_xyy_0_xy_0 = pbuffer.data(idx_eri_0_sfsd + 19);

    auto g_0_xyy_0_yy_0 = pbuffer.data(idx_eri_0_sfsd + 21);

    auto g_0_xyy_0_yz_0 = pbuffer.data(idx_eri_0_sfsd + 22);

    auto g_0_xyy_0_zz_0 = pbuffer.data(idx_eri_0_sfsd + 23);

    auto g_0_xzz_0_xz_0 = pbuffer.data(idx_eri_0_sfsd + 32);

    auto g_0_xzz_0_yy_0 = pbuffer.data(idx_eri_0_sfsd + 33);

    auto g_0_xzz_0_yz_0 = pbuffer.data(idx_eri_0_sfsd + 34);

    auto g_0_xzz_0_zz_0 = pbuffer.data(idx_eri_0_sfsd + 35);

    auto g_0_yyy_0_xx_0 = pbuffer.data(idx_eri_0_sfsd + 36);

    auto g_0_yyy_0_xy_0 = pbuffer.data(idx_eri_0_sfsd + 37);

    auto g_0_yyy_0_xz_0 = pbuffer.data(idx_eri_0_sfsd + 38);

    auto g_0_yyy_0_yy_0 = pbuffer.data(idx_eri_0_sfsd + 39);

    auto g_0_yyy_0_yz_0 = pbuffer.data(idx_eri_0_sfsd + 40);

    auto g_0_yyy_0_zz_0 = pbuffer.data(idx_eri_0_sfsd + 41);

    auto g_0_yyz_0_xy_0 = pbuffer.data(idx_eri_0_sfsd + 43);

    auto g_0_yyz_0_yy_0 = pbuffer.data(idx_eri_0_sfsd + 45);

    auto g_0_yzz_0_xx_0 = pbuffer.data(idx_eri_0_sfsd + 48);

    auto g_0_yzz_0_xz_0 = pbuffer.data(idx_eri_0_sfsd + 50);

    auto g_0_yzz_0_yz_0 = pbuffer.data(idx_eri_0_sfsd + 52);

    auto g_0_yzz_0_zz_0 = pbuffer.data(idx_eri_0_sfsd + 53);

    auto g_0_zzz_0_xx_0 = pbuffer.data(idx_eri_0_sfsd + 54);

    auto g_0_zzz_0_xy_0 = pbuffer.data(idx_eri_0_sfsd + 55);

    auto g_0_zzz_0_xz_0 = pbuffer.data(idx_eri_0_sfsd + 56);

    auto g_0_zzz_0_yy_0 = pbuffer.data(idx_eri_0_sfsd + 57);

    auto g_0_zzz_0_yz_0 = pbuffer.data(idx_eri_0_sfsd + 58);

    auto g_0_zzz_0_zz_0 = pbuffer.data(idx_eri_0_sfsd + 59);

    /// Set up components of auxilary buffer : SFSD

    auto g_0_xxx_0_xx_1 = pbuffer.data(idx_eri_1_sfsd);

    auto g_0_xxx_0_xy_1 = pbuffer.data(idx_eri_1_sfsd + 1);

    auto g_0_xxx_0_xz_1 = pbuffer.data(idx_eri_1_sfsd + 2);

    auto g_0_xxx_0_yy_1 = pbuffer.data(idx_eri_1_sfsd + 3);

    auto g_0_xxx_0_yz_1 = pbuffer.data(idx_eri_1_sfsd + 4);

    auto g_0_xxx_0_zz_1 = pbuffer.data(idx_eri_1_sfsd + 5);

    auto g_0_xxy_0_xx_1 = pbuffer.data(idx_eri_1_sfsd + 6);

    auto g_0_xxy_0_xz_1 = pbuffer.data(idx_eri_1_sfsd + 8);

    auto g_0_xxz_0_xx_1 = pbuffer.data(idx_eri_1_sfsd + 12);

    auto g_0_xxz_0_xy_1 = pbuffer.data(idx_eri_1_sfsd + 13);

    auto g_0_xyy_0_xy_1 = pbuffer.data(idx_eri_1_sfsd + 19);

    auto g_0_xyy_0_yy_1 = pbuffer.data(idx_eri_1_sfsd + 21);

    auto g_0_xyy_0_yz_1 = pbuffer.data(idx_eri_1_sfsd + 22);

    auto g_0_xyy_0_zz_1 = pbuffer.data(idx_eri_1_sfsd + 23);

    auto g_0_xzz_0_xz_1 = pbuffer.data(idx_eri_1_sfsd + 32);

    auto g_0_xzz_0_yy_1 = pbuffer.data(idx_eri_1_sfsd + 33);

    auto g_0_xzz_0_yz_1 = pbuffer.data(idx_eri_1_sfsd + 34);

    auto g_0_xzz_0_zz_1 = pbuffer.data(idx_eri_1_sfsd + 35);

    auto g_0_yyy_0_xx_1 = pbuffer.data(idx_eri_1_sfsd + 36);

    auto g_0_yyy_0_xy_1 = pbuffer.data(idx_eri_1_sfsd + 37);

    auto g_0_yyy_0_xz_1 = pbuffer.data(idx_eri_1_sfsd + 38);

    auto g_0_yyy_0_yy_1 = pbuffer.data(idx_eri_1_sfsd + 39);

    auto g_0_yyy_0_yz_1 = pbuffer.data(idx_eri_1_sfsd + 40);

    auto g_0_yyy_0_zz_1 = pbuffer.data(idx_eri_1_sfsd + 41);

    auto g_0_yyz_0_xy_1 = pbuffer.data(idx_eri_1_sfsd + 43);

    auto g_0_yyz_0_yy_1 = pbuffer.data(idx_eri_1_sfsd + 45);

    auto g_0_yzz_0_xx_1 = pbuffer.data(idx_eri_1_sfsd + 48);

    auto g_0_yzz_0_xz_1 = pbuffer.data(idx_eri_1_sfsd + 50);

    auto g_0_yzz_0_yz_1 = pbuffer.data(idx_eri_1_sfsd + 52);

    auto g_0_yzz_0_zz_1 = pbuffer.data(idx_eri_1_sfsd + 53);

    auto g_0_zzz_0_xx_1 = pbuffer.data(idx_eri_1_sfsd + 54);

    auto g_0_zzz_0_xy_1 = pbuffer.data(idx_eri_1_sfsd + 55);

    auto g_0_zzz_0_xz_1 = pbuffer.data(idx_eri_1_sfsd + 56);

    auto g_0_zzz_0_yy_1 = pbuffer.data(idx_eri_1_sfsd + 57);

    auto g_0_zzz_0_yz_1 = pbuffer.data(idx_eri_1_sfsd + 58);

    auto g_0_zzz_0_zz_1 = pbuffer.data(idx_eri_1_sfsd + 59);

    /// Set up components of auxilary buffer : SGSP

    auto g_0_xxxx_0_x_1 = pbuffer.data(idx_eri_1_sgsp);

    auto g_0_xxxx_0_y_1 = pbuffer.data(idx_eri_1_sgsp + 1);

    auto g_0_xxxx_0_z_1 = pbuffer.data(idx_eri_1_sgsp + 2);

    auto g_0_xxxz_0_z_1 = pbuffer.data(idx_eri_1_sgsp + 8);

    auto g_0_xxyy_0_x_1 = pbuffer.data(idx_eri_1_sgsp + 9);

    auto g_0_xxyy_0_y_1 = pbuffer.data(idx_eri_1_sgsp + 10);

    auto g_0_xxyy_0_z_1 = pbuffer.data(idx_eri_1_sgsp + 11);

    auto g_0_xxzz_0_x_1 = pbuffer.data(idx_eri_1_sgsp + 15);

    auto g_0_xxzz_0_y_1 = pbuffer.data(idx_eri_1_sgsp + 16);

    auto g_0_xxzz_0_z_1 = pbuffer.data(idx_eri_1_sgsp + 17);

    auto g_0_xyyy_0_y_1 = pbuffer.data(idx_eri_1_sgsp + 19);

    auto g_0_xzzz_0_z_1 = pbuffer.data(idx_eri_1_sgsp + 29);

    auto g_0_yyyy_0_x_1 = pbuffer.data(idx_eri_1_sgsp + 30);

    auto g_0_yyyy_0_y_1 = pbuffer.data(idx_eri_1_sgsp + 31);

    auto g_0_yyyy_0_z_1 = pbuffer.data(idx_eri_1_sgsp + 32);

    auto g_0_yyyz_0_z_1 = pbuffer.data(idx_eri_1_sgsp + 35);

    auto g_0_yyzz_0_x_1 = pbuffer.data(idx_eri_1_sgsp + 36);

    auto g_0_yyzz_0_y_1 = pbuffer.data(idx_eri_1_sgsp + 37);

    auto g_0_yyzz_0_z_1 = pbuffer.data(idx_eri_1_sgsp + 38);

    auto g_0_yzzz_0_y_1 = pbuffer.data(idx_eri_1_sgsp + 40);

    auto g_0_yzzz_0_z_1 = pbuffer.data(idx_eri_1_sgsp + 41);

    auto g_0_zzzz_0_x_1 = pbuffer.data(idx_eri_1_sgsp + 42);

    auto g_0_zzzz_0_y_1 = pbuffer.data(idx_eri_1_sgsp + 43);

    auto g_0_zzzz_0_z_1 = pbuffer.data(idx_eri_1_sgsp + 44);

    /// Set up components of auxilary buffer : SGSD

    auto g_0_xxxx_0_xx_0 = pbuffer.data(idx_eri_0_sgsd);

    auto g_0_xxxx_0_xy_0 = pbuffer.data(idx_eri_0_sgsd + 1);

    auto g_0_xxxx_0_xz_0 = pbuffer.data(idx_eri_0_sgsd + 2);

    auto g_0_xxxx_0_yy_0 = pbuffer.data(idx_eri_0_sgsd + 3);

    auto g_0_xxxx_0_yz_0 = pbuffer.data(idx_eri_0_sgsd + 4);

    auto g_0_xxxx_0_zz_0 = pbuffer.data(idx_eri_0_sgsd + 5);

    auto g_0_xxxy_0_xx_0 = pbuffer.data(idx_eri_0_sgsd + 6);

    auto g_0_xxxy_0_xy_0 = pbuffer.data(idx_eri_0_sgsd + 7);

    auto g_0_xxxy_0_xz_0 = pbuffer.data(idx_eri_0_sgsd + 8);

    auto g_0_xxxy_0_yy_0 = pbuffer.data(idx_eri_0_sgsd + 9);

    auto g_0_xxxz_0_xx_0 = pbuffer.data(idx_eri_0_sgsd + 12);

    auto g_0_xxxz_0_xy_0 = pbuffer.data(idx_eri_0_sgsd + 13);

    auto g_0_xxxz_0_xz_0 = pbuffer.data(idx_eri_0_sgsd + 14);

    auto g_0_xxxz_0_yz_0 = pbuffer.data(idx_eri_0_sgsd + 16);

    auto g_0_xxxz_0_zz_0 = pbuffer.data(idx_eri_0_sgsd + 17);

    auto g_0_xxyy_0_xx_0 = pbuffer.data(idx_eri_0_sgsd + 18);

    auto g_0_xxyy_0_xy_0 = pbuffer.data(idx_eri_0_sgsd + 19);

    auto g_0_xxyy_0_xz_0 = pbuffer.data(idx_eri_0_sgsd + 20);

    auto g_0_xxyy_0_yy_0 = pbuffer.data(idx_eri_0_sgsd + 21);

    auto g_0_xxyy_0_yz_0 = pbuffer.data(idx_eri_0_sgsd + 22);

    auto g_0_xxyy_0_zz_0 = pbuffer.data(idx_eri_0_sgsd + 23);

    auto g_0_xxzz_0_xx_0 = pbuffer.data(idx_eri_0_sgsd + 30);

    auto g_0_xxzz_0_xy_0 = pbuffer.data(idx_eri_0_sgsd + 31);

    auto g_0_xxzz_0_xz_0 = pbuffer.data(idx_eri_0_sgsd + 32);

    auto g_0_xxzz_0_yy_0 = pbuffer.data(idx_eri_0_sgsd + 33);

    auto g_0_xxzz_0_yz_0 = pbuffer.data(idx_eri_0_sgsd + 34);

    auto g_0_xxzz_0_zz_0 = pbuffer.data(idx_eri_0_sgsd + 35);

    auto g_0_xyyy_0_xx_0 = pbuffer.data(idx_eri_0_sgsd + 36);

    auto g_0_xyyy_0_xy_0 = pbuffer.data(idx_eri_0_sgsd + 37);

    auto g_0_xyyy_0_yy_0 = pbuffer.data(idx_eri_0_sgsd + 39);

    auto g_0_xyyy_0_yz_0 = pbuffer.data(idx_eri_0_sgsd + 40);

    auto g_0_xyyy_0_zz_0 = pbuffer.data(idx_eri_0_sgsd + 41);

    auto g_0_xzzz_0_xx_0 = pbuffer.data(idx_eri_0_sgsd + 54);

    auto g_0_xzzz_0_xz_0 = pbuffer.data(idx_eri_0_sgsd + 56);

    auto g_0_xzzz_0_yy_0 = pbuffer.data(idx_eri_0_sgsd + 57);

    auto g_0_xzzz_0_yz_0 = pbuffer.data(idx_eri_0_sgsd + 58);

    auto g_0_xzzz_0_zz_0 = pbuffer.data(idx_eri_0_sgsd + 59);

    auto g_0_yyyy_0_xx_0 = pbuffer.data(idx_eri_0_sgsd + 60);

    auto g_0_yyyy_0_xy_0 = pbuffer.data(idx_eri_0_sgsd + 61);

    auto g_0_yyyy_0_xz_0 = pbuffer.data(idx_eri_0_sgsd + 62);

    auto g_0_yyyy_0_yy_0 = pbuffer.data(idx_eri_0_sgsd + 63);

    auto g_0_yyyy_0_yz_0 = pbuffer.data(idx_eri_0_sgsd + 64);

    auto g_0_yyyy_0_zz_0 = pbuffer.data(idx_eri_0_sgsd + 65);

    auto g_0_yyyz_0_xy_0 = pbuffer.data(idx_eri_0_sgsd + 67);

    auto g_0_yyyz_0_xz_0 = pbuffer.data(idx_eri_0_sgsd + 68);

    auto g_0_yyyz_0_yy_0 = pbuffer.data(idx_eri_0_sgsd + 69);

    auto g_0_yyyz_0_yz_0 = pbuffer.data(idx_eri_0_sgsd + 70);

    auto g_0_yyyz_0_zz_0 = pbuffer.data(idx_eri_0_sgsd + 71);

    auto g_0_yyzz_0_xx_0 = pbuffer.data(idx_eri_0_sgsd + 72);

    auto g_0_yyzz_0_xy_0 = pbuffer.data(idx_eri_0_sgsd + 73);

    auto g_0_yyzz_0_xz_0 = pbuffer.data(idx_eri_0_sgsd + 74);

    auto g_0_yyzz_0_yy_0 = pbuffer.data(idx_eri_0_sgsd + 75);

    auto g_0_yyzz_0_yz_0 = pbuffer.data(idx_eri_0_sgsd + 76);

    auto g_0_yyzz_0_zz_0 = pbuffer.data(idx_eri_0_sgsd + 77);

    auto g_0_yzzz_0_xx_0 = pbuffer.data(idx_eri_0_sgsd + 78);

    auto g_0_yzzz_0_xy_0 = pbuffer.data(idx_eri_0_sgsd + 79);

    auto g_0_yzzz_0_xz_0 = pbuffer.data(idx_eri_0_sgsd + 80);

    auto g_0_yzzz_0_yy_0 = pbuffer.data(idx_eri_0_sgsd + 81);

    auto g_0_yzzz_0_yz_0 = pbuffer.data(idx_eri_0_sgsd + 82);

    auto g_0_yzzz_0_zz_0 = pbuffer.data(idx_eri_0_sgsd + 83);

    auto g_0_zzzz_0_xx_0 = pbuffer.data(idx_eri_0_sgsd + 84);

    auto g_0_zzzz_0_xy_0 = pbuffer.data(idx_eri_0_sgsd + 85);

    auto g_0_zzzz_0_xz_0 = pbuffer.data(idx_eri_0_sgsd + 86);

    auto g_0_zzzz_0_yy_0 = pbuffer.data(idx_eri_0_sgsd + 87);

    auto g_0_zzzz_0_yz_0 = pbuffer.data(idx_eri_0_sgsd + 88);

    auto g_0_zzzz_0_zz_0 = pbuffer.data(idx_eri_0_sgsd + 89);

    /// Set up components of auxilary buffer : SGSD

    auto g_0_xxxx_0_xx_1 = pbuffer.data(idx_eri_1_sgsd);

    auto g_0_xxxx_0_xy_1 = pbuffer.data(idx_eri_1_sgsd + 1);

    auto g_0_xxxx_0_xz_1 = pbuffer.data(idx_eri_1_sgsd + 2);

    auto g_0_xxxx_0_yy_1 = pbuffer.data(idx_eri_1_sgsd + 3);

    auto g_0_xxxx_0_yz_1 = pbuffer.data(idx_eri_1_sgsd + 4);

    auto g_0_xxxx_0_zz_1 = pbuffer.data(idx_eri_1_sgsd + 5);

    auto g_0_xxxy_0_xx_1 = pbuffer.data(idx_eri_1_sgsd + 6);

    auto g_0_xxxy_0_xy_1 = pbuffer.data(idx_eri_1_sgsd + 7);

    auto g_0_xxxy_0_xz_1 = pbuffer.data(idx_eri_1_sgsd + 8);

    auto g_0_xxxy_0_yy_1 = pbuffer.data(idx_eri_1_sgsd + 9);

    auto g_0_xxxz_0_xx_1 = pbuffer.data(idx_eri_1_sgsd + 12);

    auto g_0_xxxz_0_xy_1 = pbuffer.data(idx_eri_1_sgsd + 13);

    auto g_0_xxxz_0_xz_1 = pbuffer.data(idx_eri_1_sgsd + 14);

    auto g_0_xxxz_0_yz_1 = pbuffer.data(idx_eri_1_sgsd + 16);

    auto g_0_xxxz_0_zz_1 = pbuffer.data(idx_eri_1_sgsd + 17);

    auto g_0_xxyy_0_xx_1 = pbuffer.data(idx_eri_1_sgsd + 18);

    auto g_0_xxyy_0_xy_1 = pbuffer.data(idx_eri_1_sgsd + 19);

    auto g_0_xxyy_0_xz_1 = pbuffer.data(idx_eri_1_sgsd + 20);

    auto g_0_xxyy_0_yy_1 = pbuffer.data(idx_eri_1_sgsd + 21);

    auto g_0_xxyy_0_yz_1 = pbuffer.data(idx_eri_1_sgsd + 22);

    auto g_0_xxyy_0_zz_1 = pbuffer.data(idx_eri_1_sgsd + 23);

    auto g_0_xxzz_0_xx_1 = pbuffer.data(idx_eri_1_sgsd + 30);

    auto g_0_xxzz_0_xy_1 = pbuffer.data(idx_eri_1_sgsd + 31);

    auto g_0_xxzz_0_xz_1 = pbuffer.data(idx_eri_1_sgsd + 32);

    auto g_0_xxzz_0_yy_1 = pbuffer.data(idx_eri_1_sgsd + 33);

    auto g_0_xxzz_0_yz_1 = pbuffer.data(idx_eri_1_sgsd + 34);

    auto g_0_xxzz_0_zz_1 = pbuffer.data(idx_eri_1_sgsd + 35);

    auto g_0_xyyy_0_xx_1 = pbuffer.data(idx_eri_1_sgsd + 36);

    auto g_0_xyyy_0_xy_1 = pbuffer.data(idx_eri_1_sgsd + 37);

    auto g_0_xyyy_0_yy_1 = pbuffer.data(idx_eri_1_sgsd + 39);

    auto g_0_xyyy_0_yz_1 = pbuffer.data(idx_eri_1_sgsd + 40);

    auto g_0_xyyy_0_zz_1 = pbuffer.data(idx_eri_1_sgsd + 41);

    auto g_0_xzzz_0_xx_1 = pbuffer.data(idx_eri_1_sgsd + 54);

    auto g_0_xzzz_0_xz_1 = pbuffer.data(idx_eri_1_sgsd + 56);

    auto g_0_xzzz_0_yy_1 = pbuffer.data(idx_eri_1_sgsd + 57);

    auto g_0_xzzz_0_yz_1 = pbuffer.data(idx_eri_1_sgsd + 58);

    auto g_0_xzzz_0_zz_1 = pbuffer.data(idx_eri_1_sgsd + 59);

    auto g_0_yyyy_0_xx_1 = pbuffer.data(idx_eri_1_sgsd + 60);

    auto g_0_yyyy_0_xy_1 = pbuffer.data(idx_eri_1_sgsd + 61);

    auto g_0_yyyy_0_xz_1 = pbuffer.data(idx_eri_1_sgsd + 62);

    auto g_0_yyyy_0_yy_1 = pbuffer.data(idx_eri_1_sgsd + 63);

    auto g_0_yyyy_0_yz_1 = pbuffer.data(idx_eri_1_sgsd + 64);

    auto g_0_yyyy_0_zz_1 = pbuffer.data(idx_eri_1_sgsd + 65);

    auto g_0_yyyz_0_xy_1 = pbuffer.data(idx_eri_1_sgsd + 67);

    auto g_0_yyyz_0_xz_1 = pbuffer.data(idx_eri_1_sgsd + 68);

    auto g_0_yyyz_0_yy_1 = pbuffer.data(idx_eri_1_sgsd + 69);

    auto g_0_yyyz_0_yz_1 = pbuffer.data(idx_eri_1_sgsd + 70);

    auto g_0_yyyz_0_zz_1 = pbuffer.data(idx_eri_1_sgsd + 71);

    auto g_0_yyzz_0_xx_1 = pbuffer.data(idx_eri_1_sgsd + 72);

    auto g_0_yyzz_0_xy_1 = pbuffer.data(idx_eri_1_sgsd + 73);

    auto g_0_yyzz_0_xz_1 = pbuffer.data(idx_eri_1_sgsd + 74);

    auto g_0_yyzz_0_yy_1 = pbuffer.data(idx_eri_1_sgsd + 75);

    auto g_0_yyzz_0_yz_1 = pbuffer.data(idx_eri_1_sgsd + 76);

    auto g_0_yyzz_0_zz_1 = pbuffer.data(idx_eri_1_sgsd + 77);

    auto g_0_yzzz_0_xx_1 = pbuffer.data(idx_eri_1_sgsd + 78);

    auto g_0_yzzz_0_xy_1 = pbuffer.data(idx_eri_1_sgsd + 79);

    auto g_0_yzzz_0_xz_1 = pbuffer.data(idx_eri_1_sgsd + 80);

    auto g_0_yzzz_0_yy_1 = pbuffer.data(idx_eri_1_sgsd + 81);

    auto g_0_yzzz_0_yz_1 = pbuffer.data(idx_eri_1_sgsd + 82);

    auto g_0_yzzz_0_zz_1 = pbuffer.data(idx_eri_1_sgsd + 83);

    auto g_0_zzzz_0_xx_1 = pbuffer.data(idx_eri_1_sgsd + 84);

    auto g_0_zzzz_0_xy_1 = pbuffer.data(idx_eri_1_sgsd + 85);

    auto g_0_zzzz_0_xz_1 = pbuffer.data(idx_eri_1_sgsd + 86);

    auto g_0_zzzz_0_yy_1 = pbuffer.data(idx_eri_1_sgsd + 87);

    auto g_0_zzzz_0_yz_1 = pbuffer.data(idx_eri_1_sgsd + 88);

    auto g_0_zzzz_0_zz_1 = pbuffer.data(idx_eri_1_sgsd + 89);

    /// Set up 0-6 components of targeted buffer : SHSD

    auto g_0_xxxxx_0_xx_0 = pbuffer.data(idx_eri_0_shsd);

    auto g_0_xxxxx_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 1);

    auto g_0_xxxxx_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 2);

    auto g_0_xxxxx_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 3);

    auto g_0_xxxxx_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 4);

    auto g_0_xxxxx_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 5);

#pragma omp simd aligned(g_0_xxx_0_xx_0,       \
                             g_0_xxx_0_xx_1,   \
                             g_0_xxx_0_xy_0,   \
                             g_0_xxx_0_xy_1,   \
                             g_0_xxx_0_xz_0,   \
                             g_0_xxx_0_xz_1,   \
                             g_0_xxx_0_yy_0,   \
                             g_0_xxx_0_yy_1,   \
                             g_0_xxx_0_yz_0,   \
                             g_0_xxx_0_yz_1,   \
                             g_0_xxx_0_zz_0,   \
                             g_0_xxx_0_zz_1,   \
                             g_0_xxxx_0_x_1,   \
                             g_0_xxxx_0_xx_0,  \
                             g_0_xxxx_0_xx_1,  \
                             g_0_xxxx_0_xy_0,  \
                             g_0_xxxx_0_xy_1,  \
                             g_0_xxxx_0_xz_0,  \
                             g_0_xxxx_0_xz_1,  \
                             g_0_xxxx_0_y_1,   \
                             g_0_xxxx_0_yy_0,  \
                             g_0_xxxx_0_yy_1,  \
                             g_0_xxxx_0_yz_0,  \
                             g_0_xxxx_0_yz_1,  \
                             g_0_xxxx_0_z_1,   \
                             g_0_xxxx_0_zz_0,  \
                             g_0_xxxx_0_zz_1,  \
                             g_0_xxxxx_0_xx_0, \
                             g_0_xxxxx_0_xy_0, \
                             g_0_xxxxx_0_xz_0, \
                             g_0_xxxxx_0_yy_0, \
                             g_0_xxxxx_0_yz_0, \
                             g_0_xxxxx_0_zz_0, \
                             wp_x,             \
                             c_exps,           \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxx_0_xx_0[i] = 4.0 * g_0_xxx_0_xx_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xx_1[i] * fti_ab_0 + 2.0 * g_0_xxxx_0_x_1[i] * fi_abcd_0 +
                              g_0_xxxx_0_xx_0[i] * pb_x + g_0_xxxx_0_xx_1[i] * wp_x[i];

        g_0_xxxxx_0_xy_0[i] = 4.0 * g_0_xxx_0_xy_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xy_1[i] * fti_ab_0 + g_0_xxxx_0_y_1[i] * fi_abcd_0 +
                              g_0_xxxx_0_xy_0[i] * pb_x + g_0_xxxx_0_xy_1[i] * wp_x[i];

        g_0_xxxxx_0_xz_0[i] = 4.0 * g_0_xxx_0_xz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xz_1[i] * fti_ab_0 + g_0_xxxx_0_z_1[i] * fi_abcd_0 +
                              g_0_xxxx_0_xz_0[i] * pb_x + g_0_xxxx_0_xz_1[i] * wp_x[i];

        g_0_xxxxx_0_yy_0[i] =
            4.0 * g_0_xxx_0_yy_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_yy_1[i] * fti_ab_0 + g_0_xxxx_0_yy_0[i] * pb_x + g_0_xxxx_0_yy_1[i] * wp_x[i];

        g_0_xxxxx_0_yz_0[i] =
            4.0 * g_0_xxx_0_yz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_yz_1[i] * fti_ab_0 + g_0_xxxx_0_yz_0[i] * pb_x + g_0_xxxx_0_yz_1[i] * wp_x[i];

        g_0_xxxxx_0_zz_0[i] =
            4.0 * g_0_xxx_0_zz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_zz_1[i] * fti_ab_0 + g_0_xxxx_0_zz_0[i] * pb_x + g_0_xxxx_0_zz_1[i] * wp_x[i];
    }

    /// Set up 6-12 components of targeted buffer : SHSD

    auto g_0_xxxxy_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 6);

    auto g_0_xxxxy_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 7);

    auto g_0_xxxxy_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 8);

    auto g_0_xxxxy_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 9);

    auto g_0_xxxxy_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 10);

    auto g_0_xxxxy_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 11);

#pragma omp simd aligned(g_0_xxxx_0_x_1,       \
                             g_0_xxxx_0_xx_0,  \
                             g_0_xxxx_0_xx_1,  \
                             g_0_xxxx_0_xy_0,  \
                             g_0_xxxx_0_xy_1,  \
                             g_0_xxxx_0_xz_0,  \
                             g_0_xxxx_0_xz_1,  \
                             g_0_xxxx_0_y_1,   \
                             g_0_xxxx_0_yy_0,  \
                             g_0_xxxx_0_yy_1,  \
                             g_0_xxxx_0_yz_0,  \
                             g_0_xxxx_0_yz_1,  \
                             g_0_xxxx_0_z_1,   \
                             g_0_xxxx_0_zz_0,  \
                             g_0_xxxx_0_zz_1,  \
                             g_0_xxxxy_0_xx_0, \
                             g_0_xxxxy_0_xy_0, \
                             g_0_xxxxy_0_xz_0, \
                             g_0_xxxxy_0_yy_0, \
                             g_0_xxxxy_0_yz_0, \
                             g_0_xxxxy_0_zz_0, \
                             wp_y,             \
                             c_exps,           \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxy_0_xx_0[i] = g_0_xxxx_0_xx_0[i] * pb_y + g_0_xxxx_0_xx_1[i] * wp_y[i];

        g_0_xxxxy_0_xy_0[i] = g_0_xxxx_0_x_1[i] * fi_abcd_0 + g_0_xxxx_0_xy_0[i] * pb_y + g_0_xxxx_0_xy_1[i] * wp_y[i];

        g_0_xxxxy_0_xz_0[i] = g_0_xxxx_0_xz_0[i] * pb_y + g_0_xxxx_0_xz_1[i] * wp_y[i];

        g_0_xxxxy_0_yy_0[i] = 2.0 * g_0_xxxx_0_y_1[i] * fi_abcd_0 + g_0_xxxx_0_yy_0[i] * pb_y + g_0_xxxx_0_yy_1[i] * wp_y[i];

        g_0_xxxxy_0_yz_0[i] = g_0_xxxx_0_z_1[i] * fi_abcd_0 + g_0_xxxx_0_yz_0[i] * pb_y + g_0_xxxx_0_yz_1[i] * wp_y[i];

        g_0_xxxxy_0_zz_0[i] = g_0_xxxx_0_zz_0[i] * pb_y + g_0_xxxx_0_zz_1[i] * wp_y[i];
    }

    /// Set up 12-18 components of targeted buffer : SHSD

    auto g_0_xxxxz_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 12);

    auto g_0_xxxxz_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 13);

    auto g_0_xxxxz_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 14);

    auto g_0_xxxxz_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 15);

    auto g_0_xxxxz_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 16);

    auto g_0_xxxxz_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 17);

#pragma omp simd aligned(g_0_xxxx_0_x_1,       \
                             g_0_xxxx_0_xx_0,  \
                             g_0_xxxx_0_xx_1,  \
                             g_0_xxxx_0_xy_0,  \
                             g_0_xxxx_0_xy_1,  \
                             g_0_xxxx_0_xz_0,  \
                             g_0_xxxx_0_xz_1,  \
                             g_0_xxxx_0_y_1,   \
                             g_0_xxxx_0_yy_0,  \
                             g_0_xxxx_0_yy_1,  \
                             g_0_xxxx_0_yz_0,  \
                             g_0_xxxx_0_yz_1,  \
                             g_0_xxxx_0_z_1,   \
                             g_0_xxxx_0_zz_0,  \
                             g_0_xxxx_0_zz_1,  \
                             g_0_xxxxz_0_xx_0, \
                             g_0_xxxxz_0_xy_0, \
                             g_0_xxxxz_0_xz_0, \
                             g_0_xxxxz_0_yy_0, \
                             g_0_xxxxz_0_yz_0, \
                             g_0_xxxxz_0_zz_0, \
                             wp_z,             \
                             c_exps,           \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxz_0_xx_0[i] = g_0_xxxx_0_xx_0[i] * pb_z + g_0_xxxx_0_xx_1[i] * wp_z[i];

        g_0_xxxxz_0_xy_0[i] = g_0_xxxx_0_xy_0[i] * pb_z + g_0_xxxx_0_xy_1[i] * wp_z[i];

        g_0_xxxxz_0_xz_0[i] = g_0_xxxx_0_x_1[i] * fi_abcd_0 + g_0_xxxx_0_xz_0[i] * pb_z + g_0_xxxx_0_xz_1[i] * wp_z[i];

        g_0_xxxxz_0_yy_0[i] = g_0_xxxx_0_yy_0[i] * pb_z + g_0_xxxx_0_yy_1[i] * wp_z[i];

        g_0_xxxxz_0_yz_0[i] = g_0_xxxx_0_y_1[i] * fi_abcd_0 + g_0_xxxx_0_yz_0[i] * pb_z + g_0_xxxx_0_yz_1[i] * wp_z[i];

        g_0_xxxxz_0_zz_0[i] = 2.0 * g_0_xxxx_0_z_1[i] * fi_abcd_0 + g_0_xxxx_0_zz_0[i] * pb_z + g_0_xxxx_0_zz_1[i] * wp_z[i];
    }

    /// Set up 18-24 components of targeted buffer : SHSD

    auto g_0_xxxyy_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 18);

    auto g_0_xxxyy_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 19);

    auto g_0_xxxyy_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 20);

    auto g_0_xxxyy_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 21);

    auto g_0_xxxyy_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 22);

    auto g_0_xxxyy_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 23);

#pragma omp simd aligned(g_0_xxx_0_xx_0,       \
                             g_0_xxx_0_xx_1,   \
                             g_0_xxx_0_xz_0,   \
                             g_0_xxx_0_xz_1,   \
                             g_0_xxxy_0_xx_0,  \
                             g_0_xxxy_0_xx_1,  \
                             g_0_xxxy_0_xz_0,  \
                             g_0_xxxy_0_xz_1,  \
                             g_0_xxxyy_0_xx_0, \
                             g_0_xxxyy_0_xy_0, \
                             g_0_xxxyy_0_xz_0, \
                             g_0_xxxyy_0_yy_0, \
                             g_0_xxxyy_0_yz_0, \
                             g_0_xxxyy_0_zz_0, \
                             g_0_xxyy_0_xy_0,  \
                             g_0_xxyy_0_xy_1,  \
                             g_0_xxyy_0_y_1,   \
                             g_0_xxyy_0_yy_0,  \
                             g_0_xxyy_0_yy_1,  \
                             g_0_xxyy_0_yz_0,  \
                             g_0_xxyy_0_yz_1,  \
                             g_0_xxyy_0_zz_0,  \
                             g_0_xxyy_0_zz_1,  \
                             g_0_xyy_0_xy_0,   \
                             g_0_xyy_0_xy_1,   \
                             g_0_xyy_0_yy_0,   \
                             g_0_xyy_0_yy_1,   \
                             g_0_xyy_0_yz_0,   \
                             g_0_xyy_0_yz_1,   \
                             g_0_xyy_0_zz_0,   \
                             g_0_xyy_0_zz_1,   \
                             wp_x,             \
                             wp_y,             \
                             c_exps,           \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxyy_0_xx_0[i] = g_0_xxx_0_xx_0[i] * fi_ab_0 - g_0_xxx_0_xx_1[i] * fti_ab_0 + g_0_xxxy_0_xx_0[i] * pb_y + g_0_xxxy_0_xx_1[i] * wp_y[i];

        g_0_xxxyy_0_xy_0[i] = 2.0 * g_0_xyy_0_xy_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xy_1[i] * fti_ab_0 + g_0_xxyy_0_y_1[i] * fi_abcd_0 +
                              g_0_xxyy_0_xy_0[i] * pb_x + g_0_xxyy_0_xy_1[i] * wp_x[i];

        g_0_xxxyy_0_xz_0[i] = g_0_xxx_0_xz_0[i] * fi_ab_0 - g_0_xxx_0_xz_1[i] * fti_ab_0 + g_0_xxxy_0_xz_0[i] * pb_y + g_0_xxxy_0_xz_1[i] * wp_y[i];

        g_0_xxxyy_0_yy_0[i] =
            2.0 * g_0_xyy_0_yy_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_yy_1[i] * fti_ab_0 + g_0_xxyy_0_yy_0[i] * pb_x + g_0_xxyy_0_yy_1[i] * wp_x[i];

        g_0_xxxyy_0_yz_0[i] =
            2.0 * g_0_xyy_0_yz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_yz_1[i] * fti_ab_0 + g_0_xxyy_0_yz_0[i] * pb_x + g_0_xxyy_0_yz_1[i] * wp_x[i];

        g_0_xxxyy_0_zz_0[i] =
            2.0 * g_0_xyy_0_zz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_zz_1[i] * fti_ab_0 + g_0_xxyy_0_zz_0[i] * pb_x + g_0_xxyy_0_zz_1[i] * wp_x[i];
    }

    /// Set up 24-30 components of targeted buffer : SHSD

    auto g_0_xxxyz_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 24);

    auto g_0_xxxyz_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 25);

    auto g_0_xxxyz_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 26);

    auto g_0_xxxyz_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 27);

    auto g_0_xxxyz_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 28);

    auto g_0_xxxyz_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 29);

#pragma omp simd aligned(g_0_xxxy_0_xy_0,      \
                             g_0_xxxy_0_xy_1,  \
                             g_0_xxxy_0_yy_0,  \
                             g_0_xxxy_0_yy_1,  \
                             g_0_xxxyz_0_xx_0, \
                             g_0_xxxyz_0_xy_0, \
                             g_0_xxxyz_0_xz_0, \
                             g_0_xxxyz_0_yy_0, \
                             g_0_xxxyz_0_yz_0, \
                             g_0_xxxyz_0_zz_0, \
                             g_0_xxxz_0_xx_0,  \
                             g_0_xxxz_0_xx_1,  \
                             g_0_xxxz_0_xz_0,  \
                             g_0_xxxz_0_xz_1,  \
                             g_0_xxxz_0_yz_0,  \
                             g_0_xxxz_0_yz_1,  \
                             g_0_xxxz_0_z_1,   \
                             g_0_xxxz_0_zz_0,  \
                             g_0_xxxz_0_zz_1,  \
                             wp_y,             \
                             wp_z,             \
                             c_exps,           \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyz_0_xx_0[i] = g_0_xxxz_0_xx_0[i] * pb_y + g_0_xxxz_0_xx_1[i] * wp_y[i];

        g_0_xxxyz_0_xy_0[i] = g_0_xxxy_0_xy_0[i] * pb_z + g_0_xxxy_0_xy_1[i] * wp_z[i];

        g_0_xxxyz_0_xz_0[i] = g_0_xxxz_0_xz_0[i] * pb_y + g_0_xxxz_0_xz_1[i] * wp_y[i];

        g_0_xxxyz_0_yy_0[i] = g_0_xxxy_0_yy_0[i] * pb_z + g_0_xxxy_0_yy_1[i] * wp_z[i];

        g_0_xxxyz_0_yz_0[i] = g_0_xxxz_0_z_1[i] * fi_abcd_0 + g_0_xxxz_0_yz_0[i] * pb_y + g_0_xxxz_0_yz_1[i] * wp_y[i];

        g_0_xxxyz_0_zz_0[i] = g_0_xxxz_0_zz_0[i] * pb_y + g_0_xxxz_0_zz_1[i] * wp_y[i];
    }

    /// Set up 30-36 components of targeted buffer : SHSD

    auto g_0_xxxzz_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 30);

    auto g_0_xxxzz_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 31);

    auto g_0_xxxzz_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 32);

    auto g_0_xxxzz_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 33);

    auto g_0_xxxzz_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 34);

    auto g_0_xxxzz_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 35);

#pragma omp simd aligned(g_0_xxx_0_xx_0,       \
                             g_0_xxx_0_xx_1,   \
                             g_0_xxx_0_xy_0,   \
                             g_0_xxx_0_xy_1,   \
                             g_0_xxxz_0_xx_0,  \
                             g_0_xxxz_0_xx_1,  \
                             g_0_xxxz_0_xy_0,  \
                             g_0_xxxz_0_xy_1,  \
                             g_0_xxxzz_0_xx_0, \
                             g_0_xxxzz_0_xy_0, \
                             g_0_xxxzz_0_xz_0, \
                             g_0_xxxzz_0_yy_0, \
                             g_0_xxxzz_0_yz_0, \
                             g_0_xxxzz_0_zz_0, \
                             g_0_xxzz_0_xz_0,  \
                             g_0_xxzz_0_xz_1,  \
                             g_0_xxzz_0_yy_0,  \
                             g_0_xxzz_0_yy_1,  \
                             g_0_xxzz_0_yz_0,  \
                             g_0_xxzz_0_yz_1,  \
                             g_0_xxzz_0_z_1,   \
                             g_0_xxzz_0_zz_0,  \
                             g_0_xxzz_0_zz_1,  \
                             g_0_xzz_0_xz_0,   \
                             g_0_xzz_0_xz_1,   \
                             g_0_xzz_0_yy_0,   \
                             g_0_xzz_0_yy_1,   \
                             g_0_xzz_0_yz_0,   \
                             g_0_xzz_0_yz_1,   \
                             g_0_xzz_0_zz_0,   \
                             g_0_xzz_0_zz_1,   \
                             wp_x,             \
                             wp_z,             \
                             c_exps,           \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxzz_0_xx_0[i] = g_0_xxx_0_xx_0[i] * fi_ab_0 - g_0_xxx_0_xx_1[i] * fti_ab_0 + g_0_xxxz_0_xx_0[i] * pb_z + g_0_xxxz_0_xx_1[i] * wp_z[i];

        g_0_xxxzz_0_xy_0[i] = g_0_xxx_0_xy_0[i] * fi_ab_0 - g_0_xxx_0_xy_1[i] * fti_ab_0 + g_0_xxxz_0_xy_0[i] * pb_z + g_0_xxxz_0_xy_1[i] * wp_z[i];

        g_0_xxxzz_0_xz_0[i] = 2.0 * g_0_xzz_0_xz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xz_1[i] * fti_ab_0 + g_0_xxzz_0_z_1[i] * fi_abcd_0 +
                              g_0_xxzz_0_xz_0[i] * pb_x + g_0_xxzz_0_xz_1[i] * wp_x[i];

        g_0_xxxzz_0_yy_0[i] =
            2.0 * g_0_xzz_0_yy_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_yy_1[i] * fti_ab_0 + g_0_xxzz_0_yy_0[i] * pb_x + g_0_xxzz_0_yy_1[i] * wp_x[i];

        g_0_xxxzz_0_yz_0[i] =
            2.0 * g_0_xzz_0_yz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_yz_1[i] * fti_ab_0 + g_0_xxzz_0_yz_0[i] * pb_x + g_0_xxzz_0_yz_1[i] * wp_x[i];

        g_0_xxxzz_0_zz_0[i] =
            2.0 * g_0_xzz_0_zz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_zz_1[i] * fti_ab_0 + g_0_xxzz_0_zz_0[i] * pb_x + g_0_xxzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 36-42 components of targeted buffer : SHSD

    auto g_0_xxyyy_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 36);

    auto g_0_xxyyy_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 37);

    auto g_0_xxyyy_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 38);

    auto g_0_xxyyy_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 39);

    auto g_0_xxyyy_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 40);

    auto g_0_xxyyy_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 41);

#pragma omp simd aligned(g_0_xxy_0_xx_0,       \
                             g_0_xxy_0_xx_1,   \
                             g_0_xxy_0_xz_0,   \
                             g_0_xxy_0_xz_1,   \
                             g_0_xxyy_0_xx_0,  \
                             g_0_xxyy_0_xx_1,  \
                             g_0_xxyy_0_xz_0,  \
                             g_0_xxyy_0_xz_1,  \
                             g_0_xxyyy_0_xx_0, \
                             g_0_xxyyy_0_xy_0, \
                             g_0_xxyyy_0_xz_0, \
                             g_0_xxyyy_0_yy_0, \
                             g_0_xxyyy_0_yz_0, \
                             g_0_xxyyy_0_zz_0, \
                             g_0_xyyy_0_xy_0,  \
                             g_0_xyyy_0_xy_1,  \
                             g_0_xyyy_0_y_1,   \
                             g_0_xyyy_0_yy_0,  \
                             g_0_xyyy_0_yy_1,  \
                             g_0_xyyy_0_yz_0,  \
                             g_0_xyyy_0_yz_1,  \
                             g_0_xyyy_0_zz_0,  \
                             g_0_xyyy_0_zz_1,  \
                             g_0_yyy_0_xy_0,   \
                             g_0_yyy_0_xy_1,   \
                             g_0_yyy_0_yy_0,   \
                             g_0_yyy_0_yy_1,   \
                             g_0_yyy_0_yz_0,   \
                             g_0_yyy_0_yz_1,   \
                             g_0_yyy_0_zz_0,   \
                             g_0_yyy_0_zz_1,   \
                             wp_x,             \
                             wp_y,             \
                             c_exps,           \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyy_0_xx_0[i] =
            2.0 * g_0_xxy_0_xx_0[i] * fi_ab_0 - 2.0 * g_0_xxy_0_xx_1[i] * fti_ab_0 + g_0_xxyy_0_xx_0[i] * pb_y + g_0_xxyy_0_xx_1[i] * wp_y[i];

        g_0_xxyyy_0_xy_0[i] = g_0_yyy_0_xy_0[i] * fi_ab_0 - g_0_yyy_0_xy_1[i] * fti_ab_0 + g_0_xyyy_0_y_1[i] * fi_abcd_0 + g_0_xyyy_0_xy_0[i] * pb_x +
                              g_0_xyyy_0_xy_1[i] * wp_x[i];

        g_0_xxyyy_0_xz_0[i] =
            2.0 * g_0_xxy_0_xz_0[i] * fi_ab_0 - 2.0 * g_0_xxy_0_xz_1[i] * fti_ab_0 + g_0_xxyy_0_xz_0[i] * pb_y + g_0_xxyy_0_xz_1[i] * wp_y[i];

        g_0_xxyyy_0_yy_0[i] = g_0_yyy_0_yy_0[i] * fi_ab_0 - g_0_yyy_0_yy_1[i] * fti_ab_0 + g_0_xyyy_0_yy_0[i] * pb_x + g_0_xyyy_0_yy_1[i] * wp_x[i];

        g_0_xxyyy_0_yz_0[i] = g_0_yyy_0_yz_0[i] * fi_ab_0 - g_0_yyy_0_yz_1[i] * fti_ab_0 + g_0_xyyy_0_yz_0[i] * pb_x + g_0_xyyy_0_yz_1[i] * wp_x[i];

        g_0_xxyyy_0_zz_0[i] = g_0_yyy_0_zz_0[i] * fi_ab_0 - g_0_yyy_0_zz_1[i] * fti_ab_0 + g_0_xyyy_0_zz_0[i] * pb_x + g_0_xyyy_0_zz_1[i] * wp_x[i];
    }

    /// Set up 42-48 components of targeted buffer : SHSD

    auto g_0_xxyyz_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 42);

    auto g_0_xxyyz_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 43);

    auto g_0_xxyyz_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 44);

    auto g_0_xxyyz_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 45);

    auto g_0_xxyyz_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 46);

    auto g_0_xxyyz_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 47);

#pragma omp simd aligned(g_0_xxyy_0_x_1,       \
                             g_0_xxyy_0_xx_0,  \
                             g_0_xxyy_0_xx_1,  \
                             g_0_xxyy_0_xy_0,  \
                             g_0_xxyy_0_xy_1,  \
                             g_0_xxyy_0_xz_0,  \
                             g_0_xxyy_0_xz_1,  \
                             g_0_xxyy_0_y_1,   \
                             g_0_xxyy_0_yy_0,  \
                             g_0_xxyy_0_yy_1,  \
                             g_0_xxyy_0_yz_0,  \
                             g_0_xxyy_0_yz_1,  \
                             g_0_xxyy_0_z_1,   \
                             g_0_xxyy_0_zz_0,  \
                             g_0_xxyy_0_zz_1,  \
                             g_0_xxyyz_0_xx_0, \
                             g_0_xxyyz_0_xy_0, \
                             g_0_xxyyz_0_xz_0, \
                             g_0_xxyyz_0_yy_0, \
                             g_0_xxyyz_0_yz_0, \
                             g_0_xxyyz_0_zz_0, \
                             wp_z,             \
                             c_exps,           \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyz_0_xx_0[i] = g_0_xxyy_0_xx_0[i] * pb_z + g_0_xxyy_0_xx_1[i] * wp_z[i];

        g_0_xxyyz_0_xy_0[i] = g_0_xxyy_0_xy_0[i] * pb_z + g_0_xxyy_0_xy_1[i] * wp_z[i];

        g_0_xxyyz_0_xz_0[i] = g_0_xxyy_0_x_1[i] * fi_abcd_0 + g_0_xxyy_0_xz_0[i] * pb_z + g_0_xxyy_0_xz_1[i] * wp_z[i];

        g_0_xxyyz_0_yy_0[i] = g_0_xxyy_0_yy_0[i] * pb_z + g_0_xxyy_0_yy_1[i] * wp_z[i];

        g_0_xxyyz_0_yz_0[i] = g_0_xxyy_0_y_1[i] * fi_abcd_0 + g_0_xxyy_0_yz_0[i] * pb_z + g_0_xxyy_0_yz_1[i] * wp_z[i];

        g_0_xxyyz_0_zz_0[i] = 2.0 * g_0_xxyy_0_z_1[i] * fi_abcd_0 + g_0_xxyy_0_zz_0[i] * pb_z + g_0_xxyy_0_zz_1[i] * wp_z[i];
    }

    /// Set up 48-54 components of targeted buffer : SHSD

    auto g_0_xxyzz_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 48);

    auto g_0_xxyzz_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 49);

    auto g_0_xxyzz_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 50);

    auto g_0_xxyzz_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 51);

    auto g_0_xxyzz_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 52);

    auto g_0_xxyzz_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 53);

#pragma omp simd aligned(g_0_xxyzz_0_xx_0,     \
                             g_0_xxyzz_0_xy_0, \
                             g_0_xxyzz_0_xz_0, \
                             g_0_xxyzz_0_yy_0, \
                             g_0_xxyzz_0_yz_0, \
                             g_0_xxyzz_0_zz_0, \
                             g_0_xxzz_0_x_1,   \
                             g_0_xxzz_0_xx_0,  \
                             g_0_xxzz_0_xx_1,  \
                             g_0_xxzz_0_xy_0,  \
                             g_0_xxzz_0_xy_1,  \
                             g_0_xxzz_0_xz_0,  \
                             g_0_xxzz_0_xz_1,  \
                             g_0_xxzz_0_y_1,   \
                             g_0_xxzz_0_yy_0,  \
                             g_0_xxzz_0_yy_1,  \
                             g_0_xxzz_0_yz_0,  \
                             g_0_xxzz_0_yz_1,  \
                             g_0_xxzz_0_z_1,   \
                             g_0_xxzz_0_zz_0,  \
                             g_0_xxzz_0_zz_1,  \
                             wp_y,             \
                             c_exps,           \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyzz_0_xx_0[i] = g_0_xxzz_0_xx_0[i] * pb_y + g_0_xxzz_0_xx_1[i] * wp_y[i];

        g_0_xxyzz_0_xy_0[i] = g_0_xxzz_0_x_1[i] * fi_abcd_0 + g_0_xxzz_0_xy_0[i] * pb_y + g_0_xxzz_0_xy_1[i] * wp_y[i];

        g_0_xxyzz_0_xz_0[i] = g_0_xxzz_0_xz_0[i] * pb_y + g_0_xxzz_0_xz_1[i] * wp_y[i];

        g_0_xxyzz_0_yy_0[i] = 2.0 * g_0_xxzz_0_y_1[i] * fi_abcd_0 + g_0_xxzz_0_yy_0[i] * pb_y + g_0_xxzz_0_yy_1[i] * wp_y[i];

        g_0_xxyzz_0_yz_0[i] = g_0_xxzz_0_z_1[i] * fi_abcd_0 + g_0_xxzz_0_yz_0[i] * pb_y + g_0_xxzz_0_yz_1[i] * wp_y[i];

        g_0_xxyzz_0_zz_0[i] = g_0_xxzz_0_zz_0[i] * pb_y + g_0_xxzz_0_zz_1[i] * wp_y[i];
    }

    /// Set up 54-60 components of targeted buffer : SHSD

    auto g_0_xxzzz_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 54);

    auto g_0_xxzzz_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 55);

    auto g_0_xxzzz_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 56);

    auto g_0_xxzzz_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 57);

    auto g_0_xxzzz_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 58);

    auto g_0_xxzzz_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 59);

#pragma omp simd aligned(g_0_xxz_0_xx_0,       \
                             g_0_xxz_0_xx_1,   \
                             g_0_xxz_0_xy_0,   \
                             g_0_xxz_0_xy_1,   \
                             g_0_xxzz_0_xx_0,  \
                             g_0_xxzz_0_xx_1,  \
                             g_0_xxzz_0_xy_0,  \
                             g_0_xxzz_0_xy_1,  \
                             g_0_xxzzz_0_xx_0, \
                             g_0_xxzzz_0_xy_0, \
                             g_0_xxzzz_0_xz_0, \
                             g_0_xxzzz_0_yy_0, \
                             g_0_xxzzz_0_yz_0, \
                             g_0_xxzzz_0_zz_0, \
                             g_0_xzzz_0_xz_0,  \
                             g_0_xzzz_0_xz_1,  \
                             g_0_xzzz_0_yy_0,  \
                             g_0_xzzz_0_yy_1,  \
                             g_0_xzzz_0_yz_0,  \
                             g_0_xzzz_0_yz_1,  \
                             g_0_xzzz_0_z_1,   \
                             g_0_xzzz_0_zz_0,  \
                             g_0_xzzz_0_zz_1,  \
                             g_0_zzz_0_xz_0,   \
                             g_0_zzz_0_xz_1,   \
                             g_0_zzz_0_yy_0,   \
                             g_0_zzz_0_yy_1,   \
                             g_0_zzz_0_yz_0,   \
                             g_0_zzz_0_yz_1,   \
                             g_0_zzz_0_zz_0,   \
                             g_0_zzz_0_zz_1,   \
                             wp_x,             \
                             wp_z,             \
                             c_exps,           \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxzzz_0_xx_0[i] =
            2.0 * g_0_xxz_0_xx_0[i] * fi_ab_0 - 2.0 * g_0_xxz_0_xx_1[i] * fti_ab_0 + g_0_xxzz_0_xx_0[i] * pb_z + g_0_xxzz_0_xx_1[i] * wp_z[i];

        g_0_xxzzz_0_xy_0[i] =
            2.0 * g_0_xxz_0_xy_0[i] * fi_ab_0 - 2.0 * g_0_xxz_0_xy_1[i] * fti_ab_0 + g_0_xxzz_0_xy_0[i] * pb_z + g_0_xxzz_0_xy_1[i] * wp_z[i];

        g_0_xxzzz_0_xz_0[i] = g_0_zzz_0_xz_0[i] * fi_ab_0 - g_0_zzz_0_xz_1[i] * fti_ab_0 + g_0_xzzz_0_z_1[i] * fi_abcd_0 + g_0_xzzz_0_xz_0[i] * pb_x +
                              g_0_xzzz_0_xz_1[i] * wp_x[i];

        g_0_xxzzz_0_yy_0[i] = g_0_zzz_0_yy_0[i] * fi_ab_0 - g_0_zzz_0_yy_1[i] * fti_ab_0 + g_0_xzzz_0_yy_0[i] * pb_x + g_0_xzzz_0_yy_1[i] * wp_x[i];

        g_0_xxzzz_0_yz_0[i] = g_0_zzz_0_yz_0[i] * fi_ab_0 - g_0_zzz_0_yz_1[i] * fti_ab_0 + g_0_xzzz_0_yz_0[i] * pb_x + g_0_xzzz_0_yz_1[i] * wp_x[i];

        g_0_xxzzz_0_zz_0[i] = g_0_zzz_0_zz_0[i] * fi_ab_0 - g_0_zzz_0_zz_1[i] * fti_ab_0 + g_0_xzzz_0_zz_0[i] * pb_x + g_0_xzzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 60-66 components of targeted buffer : SHSD

    auto g_0_xyyyy_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 60);

    auto g_0_xyyyy_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 61);

    auto g_0_xyyyy_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 62);

    auto g_0_xyyyy_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 63);

    auto g_0_xyyyy_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 64);

    auto g_0_xyyyy_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 65);

#pragma omp simd aligned(g_0_xyyyy_0_xx_0,     \
                             g_0_xyyyy_0_xy_0, \
                             g_0_xyyyy_0_xz_0, \
                             g_0_xyyyy_0_yy_0, \
                             g_0_xyyyy_0_yz_0, \
                             g_0_xyyyy_0_zz_0, \
                             g_0_yyyy_0_x_1,   \
                             g_0_yyyy_0_xx_0,  \
                             g_0_yyyy_0_xx_1,  \
                             g_0_yyyy_0_xy_0,  \
                             g_0_yyyy_0_xy_1,  \
                             g_0_yyyy_0_xz_0,  \
                             g_0_yyyy_0_xz_1,  \
                             g_0_yyyy_0_y_1,   \
                             g_0_yyyy_0_yy_0,  \
                             g_0_yyyy_0_yy_1,  \
                             g_0_yyyy_0_yz_0,  \
                             g_0_yyyy_0_yz_1,  \
                             g_0_yyyy_0_z_1,   \
                             g_0_yyyy_0_zz_0,  \
                             g_0_yyyy_0_zz_1,  \
                             wp_x,             \
                             c_exps,           \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyy_0_xx_0[i] = 2.0 * g_0_yyyy_0_x_1[i] * fi_abcd_0 + g_0_yyyy_0_xx_0[i] * pb_x + g_0_yyyy_0_xx_1[i] * wp_x[i];

        g_0_xyyyy_0_xy_0[i] = g_0_yyyy_0_y_1[i] * fi_abcd_0 + g_0_yyyy_0_xy_0[i] * pb_x + g_0_yyyy_0_xy_1[i] * wp_x[i];

        g_0_xyyyy_0_xz_0[i] = g_0_yyyy_0_z_1[i] * fi_abcd_0 + g_0_yyyy_0_xz_0[i] * pb_x + g_0_yyyy_0_xz_1[i] * wp_x[i];

        g_0_xyyyy_0_yy_0[i] = g_0_yyyy_0_yy_0[i] * pb_x + g_0_yyyy_0_yy_1[i] * wp_x[i];

        g_0_xyyyy_0_yz_0[i] = g_0_yyyy_0_yz_0[i] * pb_x + g_0_yyyy_0_yz_1[i] * wp_x[i];

        g_0_xyyyy_0_zz_0[i] = g_0_yyyy_0_zz_0[i] * pb_x + g_0_yyyy_0_zz_1[i] * wp_x[i];
    }

    /// Set up 66-72 components of targeted buffer : SHSD

    auto g_0_xyyyz_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 66);

    auto g_0_xyyyz_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 67);

    auto g_0_xyyyz_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 68);

    auto g_0_xyyyz_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 69);

    auto g_0_xyyyz_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 70);

    auto g_0_xyyyz_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 71);

#pragma omp simd aligned(g_0_xyyy_0_xx_0,      \
                             g_0_xyyy_0_xx_1,  \
                             g_0_xyyy_0_xy_0,  \
                             g_0_xyyy_0_xy_1,  \
                             g_0_xyyyz_0_xx_0, \
                             g_0_xyyyz_0_xy_0, \
                             g_0_xyyyz_0_xz_0, \
                             g_0_xyyyz_0_yy_0, \
                             g_0_xyyyz_0_yz_0, \
                             g_0_xyyyz_0_zz_0, \
                             g_0_yyyz_0_xz_0,  \
                             g_0_yyyz_0_xz_1,  \
                             g_0_yyyz_0_yy_0,  \
                             g_0_yyyz_0_yy_1,  \
                             g_0_yyyz_0_yz_0,  \
                             g_0_yyyz_0_yz_1,  \
                             g_0_yyyz_0_z_1,   \
                             g_0_yyyz_0_zz_0,  \
                             g_0_yyyz_0_zz_1,  \
                             wp_x,             \
                             wp_z,             \
                             c_exps,           \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyz_0_xx_0[i] = g_0_xyyy_0_xx_0[i] * pb_z + g_0_xyyy_0_xx_1[i] * wp_z[i];

        g_0_xyyyz_0_xy_0[i] = g_0_xyyy_0_xy_0[i] * pb_z + g_0_xyyy_0_xy_1[i] * wp_z[i];

        g_0_xyyyz_0_xz_0[i] = g_0_yyyz_0_z_1[i] * fi_abcd_0 + g_0_yyyz_0_xz_0[i] * pb_x + g_0_yyyz_0_xz_1[i] * wp_x[i];

        g_0_xyyyz_0_yy_0[i] = g_0_yyyz_0_yy_0[i] * pb_x + g_0_yyyz_0_yy_1[i] * wp_x[i];

        g_0_xyyyz_0_yz_0[i] = g_0_yyyz_0_yz_0[i] * pb_x + g_0_yyyz_0_yz_1[i] * wp_x[i];

        g_0_xyyyz_0_zz_0[i] = g_0_yyyz_0_zz_0[i] * pb_x + g_0_yyyz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 72-78 components of targeted buffer : SHSD

    auto g_0_xyyzz_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 72);

    auto g_0_xyyzz_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 73);

    auto g_0_xyyzz_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 74);

    auto g_0_xyyzz_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 75);

    auto g_0_xyyzz_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 76);

    auto g_0_xyyzz_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 77);

#pragma omp simd aligned(g_0_xyyzz_0_xx_0,     \
                             g_0_xyyzz_0_xy_0, \
                             g_0_xyyzz_0_xz_0, \
                             g_0_xyyzz_0_yy_0, \
                             g_0_xyyzz_0_yz_0, \
                             g_0_xyyzz_0_zz_0, \
                             g_0_yyzz_0_x_1,   \
                             g_0_yyzz_0_xx_0,  \
                             g_0_yyzz_0_xx_1,  \
                             g_0_yyzz_0_xy_0,  \
                             g_0_yyzz_0_xy_1,  \
                             g_0_yyzz_0_xz_0,  \
                             g_0_yyzz_0_xz_1,  \
                             g_0_yyzz_0_y_1,   \
                             g_0_yyzz_0_yy_0,  \
                             g_0_yyzz_0_yy_1,  \
                             g_0_yyzz_0_yz_0,  \
                             g_0_yyzz_0_yz_1,  \
                             g_0_yyzz_0_z_1,   \
                             g_0_yyzz_0_zz_0,  \
                             g_0_yyzz_0_zz_1,  \
                             wp_x,             \
                             c_exps,           \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyzz_0_xx_0[i] = 2.0 * g_0_yyzz_0_x_1[i] * fi_abcd_0 + g_0_yyzz_0_xx_0[i] * pb_x + g_0_yyzz_0_xx_1[i] * wp_x[i];

        g_0_xyyzz_0_xy_0[i] = g_0_yyzz_0_y_1[i] * fi_abcd_0 + g_0_yyzz_0_xy_0[i] * pb_x + g_0_yyzz_0_xy_1[i] * wp_x[i];

        g_0_xyyzz_0_xz_0[i] = g_0_yyzz_0_z_1[i] * fi_abcd_0 + g_0_yyzz_0_xz_0[i] * pb_x + g_0_yyzz_0_xz_1[i] * wp_x[i];

        g_0_xyyzz_0_yy_0[i] = g_0_yyzz_0_yy_0[i] * pb_x + g_0_yyzz_0_yy_1[i] * wp_x[i];

        g_0_xyyzz_0_yz_0[i] = g_0_yyzz_0_yz_0[i] * pb_x + g_0_yyzz_0_yz_1[i] * wp_x[i];

        g_0_xyyzz_0_zz_0[i] = g_0_yyzz_0_zz_0[i] * pb_x + g_0_yyzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 78-84 components of targeted buffer : SHSD

    auto g_0_xyzzz_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 78);

    auto g_0_xyzzz_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 79);

    auto g_0_xyzzz_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 80);

    auto g_0_xyzzz_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 81);

    auto g_0_xyzzz_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 82);

    auto g_0_xyzzz_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 83);

#pragma omp simd aligned(g_0_xyzzz_0_xx_0,     \
                             g_0_xyzzz_0_xy_0, \
                             g_0_xyzzz_0_xz_0, \
                             g_0_xyzzz_0_yy_0, \
                             g_0_xyzzz_0_yz_0, \
                             g_0_xyzzz_0_zz_0, \
                             g_0_xzzz_0_xx_0,  \
                             g_0_xzzz_0_xx_1,  \
                             g_0_xzzz_0_xz_0,  \
                             g_0_xzzz_0_xz_1,  \
                             g_0_yzzz_0_xy_0,  \
                             g_0_yzzz_0_xy_1,  \
                             g_0_yzzz_0_y_1,   \
                             g_0_yzzz_0_yy_0,  \
                             g_0_yzzz_0_yy_1,  \
                             g_0_yzzz_0_yz_0,  \
                             g_0_yzzz_0_yz_1,  \
                             g_0_yzzz_0_zz_0,  \
                             g_0_yzzz_0_zz_1,  \
                             wp_x,             \
                             wp_y,             \
                             c_exps,           \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyzzz_0_xx_0[i] = g_0_xzzz_0_xx_0[i] * pb_y + g_0_xzzz_0_xx_1[i] * wp_y[i];

        g_0_xyzzz_0_xy_0[i] = g_0_yzzz_0_y_1[i] * fi_abcd_0 + g_0_yzzz_0_xy_0[i] * pb_x + g_0_yzzz_0_xy_1[i] * wp_x[i];

        g_0_xyzzz_0_xz_0[i] = g_0_xzzz_0_xz_0[i] * pb_y + g_0_xzzz_0_xz_1[i] * wp_y[i];

        g_0_xyzzz_0_yy_0[i] = g_0_yzzz_0_yy_0[i] * pb_x + g_0_yzzz_0_yy_1[i] * wp_x[i];

        g_0_xyzzz_0_yz_0[i] = g_0_yzzz_0_yz_0[i] * pb_x + g_0_yzzz_0_yz_1[i] * wp_x[i];

        g_0_xyzzz_0_zz_0[i] = g_0_yzzz_0_zz_0[i] * pb_x + g_0_yzzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 84-90 components of targeted buffer : SHSD

    auto g_0_xzzzz_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 84);

    auto g_0_xzzzz_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 85);

    auto g_0_xzzzz_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 86);

    auto g_0_xzzzz_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 87);

    auto g_0_xzzzz_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 88);

    auto g_0_xzzzz_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 89);

#pragma omp simd aligned(g_0_xzzzz_0_xx_0,     \
                             g_0_xzzzz_0_xy_0, \
                             g_0_xzzzz_0_xz_0, \
                             g_0_xzzzz_0_yy_0, \
                             g_0_xzzzz_0_yz_0, \
                             g_0_xzzzz_0_zz_0, \
                             g_0_zzzz_0_x_1,   \
                             g_0_zzzz_0_xx_0,  \
                             g_0_zzzz_0_xx_1,  \
                             g_0_zzzz_0_xy_0,  \
                             g_0_zzzz_0_xy_1,  \
                             g_0_zzzz_0_xz_0,  \
                             g_0_zzzz_0_xz_1,  \
                             g_0_zzzz_0_y_1,   \
                             g_0_zzzz_0_yy_0,  \
                             g_0_zzzz_0_yy_1,  \
                             g_0_zzzz_0_yz_0,  \
                             g_0_zzzz_0_yz_1,  \
                             g_0_zzzz_0_z_1,   \
                             g_0_zzzz_0_zz_0,  \
                             g_0_zzzz_0_zz_1,  \
                             wp_x,             \
                             c_exps,           \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzzzz_0_xx_0[i] = 2.0 * g_0_zzzz_0_x_1[i] * fi_abcd_0 + g_0_zzzz_0_xx_0[i] * pb_x + g_0_zzzz_0_xx_1[i] * wp_x[i];

        g_0_xzzzz_0_xy_0[i] = g_0_zzzz_0_y_1[i] * fi_abcd_0 + g_0_zzzz_0_xy_0[i] * pb_x + g_0_zzzz_0_xy_1[i] * wp_x[i];

        g_0_xzzzz_0_xz_0[i] = g_0_zzzz_0_z_1[i] * fi_abcd_0 + g_0_zzzz_0_xz_0[i] * pb_x + g_0_zzzz_0_xz_1[i] * wp_x[i];

        g_0_xzzzz_0_yy_0[i] = g_0_zzzz_0_yy_0[i] * pb_x + g_0_zzzz_0_yy_1[i] * wp_x[i];

        g_0_xzzzz_0_yz_0[i] = g_0_zzzz_0_yz_0[i] * pb_x + g_0_zzzz_0_yz_1[i] * wp_x[i];

        g_0_xzzzz_0_zz_0[i] = g_0_zzzz_0_zz_0[i] * pb_x + g_0_zzzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 90-96 components of targeted buffer : SHSD

    auto g_0_yyyyy_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 90);

    auto g_0_yyyyy_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 91);

    auto g_0_yyyyy_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 92);

    auto g_0_yyyyy_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 93);

    auto g_0_yyyyy_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 94);

    auto g_0_yyyyy_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 95);

#pragma omp simd aligned(g_0_yyy_0_xx_0,       \
                             g_0_yyy_0_xx_1,   \
                             g_0_yyy_0_xy_0,   \
                             g_0_yyy_0_xy_1,   \
                             g_0_yyy_0_xz_0,   \
                             g_0_yyy_0_xz_1,   \
                             g_0_yyy_0_yy_0,   \
                             g_0_yyy_0_yy_1,   \
                             g_0_yyy_0_yz_0,   \
                             g_0_yyy_0_yz_1,   \
                             g_0_yyy_0_zz_0,   \
                             g_0_yyy_0_zz_1,   \
                             g_0_yyyy_0_x_1,   \
                             g_0_yyyy_0_xx_0,  \
                             g_0_yyyy_0_xx_1,  \
                             g_0_yyyy_0_xy_0,  \
                             g_0_yyyy_0_xy_1,  \
                             g_0_yyyy_0_xz_0,  \
                             g_0_yyyy_0_xz_1,  \
                             g_0_yyyy_0_y_1,   \
                             g_0_yyyy_0_yy_0,  \
                             g_0_yyyy_0_yy_1,  \
                             g_0_yyyy_0_yz_0,  \
                             g_0_yyyy_0_yz_1,  \
                             g_0_yyyy_0_z_1,   \
                             g_0_yyyy_0_zz_0,  \
                             g_0_yyyy_0_zz_1,  \
                             g_0_yyyyy_0_xx_0, \
                             g_0_yyyyy_0_xy_0, \
                             g_0_yyyyy_0_xz_0, \
                             g_0_yyyyy_0_yy_0, \
                             g_0_yyyyy_0_yz_0, \
                             g_0_yyyyy_0_zz_0, \
                             wp_y,             \
                             c_exps,           \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyy_0_xx_0[i] =
            4.0 * g_0_yyy_0_xx_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xx_1[i] * fti_ab_0 + g_0_yyyy_0_xx_0[i] * pb_y + g_0_yyyy_0_xx_1[i] * wp_y[i];

        g_0_yyyyy_0_xy_0[i] = 4.0 * g_0_yyy_0_xy_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xy_1[i] * fti_ab_0 + g_0_yyyy_0_x_1[i] * fi_abcd_0 +
                              g_0_yyyy_0_xy_0[i] * pb_y + g_0_yyyy_0_xy_1[i] * wp_y[i];

        g_0_yyyyy_0_xz_0[i] =
            4.0 * g_0_yyy_0_xz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xz_1[i] * fti_ab_0 + g_0_yyyy_0_xz_0[i] * pb_y + g_0_yyyy_0_xz_1[i] * wp_y[i];

        g_0_yyyyy_0_yy_0[i] = 4.0 * g_0_yyy_0_yy_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_yy_1[i] * fti_ab_0 + 2.0 * g_0_yyyy_0_y_1[i] * fi_abcd_0 +
                              g_0_yyyy_0_yy_0[i] * pb_y + g_0_yyyy_0_yy_1[i] * wp_y[i];

        g_0_yyyyy_0_yz_0[i] = 4.0 * g_0_yyy_0_yz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_yz_1[i] * fti_ab_0 + g_0_yyyy_0_z_1[i] * fi_abcd_0 +
                              g_0_yyyy_0_yz_0[i] * pb_y + g_0_yyyy_0_yz_1[i] * wp_y[i];

        g_0_yyyyy_0_zz_0[i] =
            4.0 * g_0_yyy_0_zz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_zz_1[i] * fti_ab_0 + g_0_yyyy_0_zz_0[i] * pb_y + g_0_yyyy_0_zz_1[i] * wp_y[i];
    }

    /// Set up 96-102 components of targeted buffer : SHSD

    auto g_0_yyyyz_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 96);

    auto g_0_yyyyz_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 97);

    auto g_0_yyyyz_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 98);

    auto g_0_yyyyz_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 99);

    auto g_0_yyyyz_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 100);

    auto g_0_yyyyz_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 101);

#pragma omp simd aligned(g_0_yyyy_0_x_1,       \
                             g_0_yyyy_0_xx_0,  \
                             g_0_yyyy_0_xx_1,  \
                             g_0_yyyy_0_xy_0,  \
                             g_0_yyyy_0_xy_1,  \
                             g_0_yyyy_0_xz_0,  \
                             g_0_yyyy_0_xz_1,  \
                             g_0_yyyy_0_y_1,   \
                             g_0_yyyy_0_yy_0,  \
                             g_0_yyyy_0_yy_1,  \
                             g_0_yyyy_0_yz_0,  \
                             g_0_yyyy_0_yz_1,  \
                             g_0_yyyy_0_z_1,   \
                             g_0_yyyy_0_zz_0,  \
                             g_0_yyyy_0_zz_1,  \
                             g_0_yyyyz_0_xx_0, \
                             g_0_yyyyz_0_xy_0, \
                             g_0_yyyyz_0_xz_0, \
                             g_0_yyyyz_0_yy_0, \
                             g_0_yyyyz_0_yz_0, \
                             g_0_yyyyz_0_zz_0, \
                             wp_z,             \
                             c_exps,           \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyyz_0_xx_0[i] = g_0_yyyy_0_xx_0[i] * pb_z + g_0_yyyy_0_xx_1[i] * wp_z[i];

        g_0_yyyyz_0_xy_0[i] = g_0_yyyy_0_xy_0[i] * pb_z + g_0_yyyy_0_xy_1[i] * wp_z[i];

        g_0_yyyyz_0_xz_0[i] = g_0_yyyy_0_x_1[i] * fi_abcd_0 + g_0_yyyy_0_xz_0[i] * pb_z + g_0_yyyy_0_xz_1[i] * wp_z[i];

        g_0_yyyyz_0_yy_0[i] = g_0_yyyy_0_yy_0[i] * pb_z + g_0_yyyy_0_yy_1[i] * wp_z[i];

        g_0_yyyyz_0_yz_0[i] = g_0_yyyy_0_y_1[i] * fi_abcd_0 + g_0_yyyy_0_yz_0[i] * pb_z + g_0_yyyy_0_yz_1[i] * wp_z[i];

        g_0_yyyyz_0_zz_0[i] = 2.0 * g_0_yyyy_0_z_1[i] * fi_abcd_0 + g_0_yyyy_0_zz_0[i] * pb_z + g_0_yyyy_0_zz_1[i] * wp_z[i];
    }

    /// Set up 102-108 components of targeted buffer : SHSD

    auto g_0_yyyzz_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 102);

    auto g_0_yyyzz_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 103);

    auto g_0_yyyzz_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 104);

    auto g_0_yyyzz_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 105);

    auto g_0_yyyzz_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 106);

    auto g_0_yyyzz_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 107);

#pragma omp simd aligned(g_0_yyy_0_xy_0,       \
                             g_0_yyy_0_xy_1,   \
                             g_0_yyy_0_yy_0,   \
                             g_0_yyy_0_yy_1,   \
                             g_0_yyyz_0_xy_0,  \
                             g_0_yyyz_0_xy_1,  \
                             g_0_yyyz_0_yy_0,  \
                             g_0_yyyz_0_yy_1,  \
                             g_0_yyyzz_0_xx_0, \
                             g_0_yyyzz_0_xy_0, \
                             g_0_yyyzz_0_xz_0, \
                             g_0_yyyzz_0_yy_0, \
                             g_0_yyyzz_0_yz_0, \
                             g_0_yyyzz_0_zz_0, \
                             g_0_yyzz_0_xx_0,  \
                             g_0_yyzz_0_xx_1,  \
                             g_0_yyzz_0_xz_0,  \
                             g_0_yyzz_0_xz_1,  \
                             g_0_yyzz_0_yz_0,  \
                             g_0_yyzz_0_yz_1,  \
                             g_0_yyzz_0_z_1,   \
                             g_0_yyzz_0_zz_0,  \
                             g_0_yyzz_0_zz_1,  \
                             g_0_yzz_0_xx_0,   \
                             g_0_yzz_0_xx_1,   \
                             g_0_yzz_0_xz_0,   \
                             g_0_yzz_0_xz_1,   \
                             g_0_yzz_0_yz_0,   \
                             g_0_yzz_0_yz_1,   \
                             g_0_yzz_0_zz_0,   \
                             g_0_yzz_0_zz_1,   \
                             wp_y,             \
                             wp_z,             \
                             c_exps,           \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyzz_0_xx_0[i] =
            2.0 * g_0_yzz_0_xx_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xx_1[i] * fti_ab_0 + g_0_yyzz_0_xx_0[i] * pb_y + g_0_yyzz_0_xx_1[i] * wp_y[i];

        g_0_yyyzz_0_xy_0[i] = g_0_yyy_0_xy_0[i] * fi_ab_0 - g_0_yyy_0_xy_1[i] * fti_ab_0 + g_0_yyyz_0_xy_0[i] * pb_z + g_0_yyyz_0_xy_1[i] * wp_z[i];

        g_0_yyyzz_0_xz_0[i] =
            2.0 * g_0_yzz_0_xz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xz_1[i] * fti_ab_0 + g_0_yyzz_0_xz_0[i] * pb_y + g_0_yyzz_0_xz_1[i] * wp_y[i];

        g_0_yyyzz_0_yy_0[i] = g_0_yyy_0_yy_0[i] * fi_ab_0 - g_0_yyy_0_yy_1[i] * fti_ab_0 + g_0_yyyz_0_yy_0[i] * pb_z + g_0_yyyz_0_yy_1[i] * wp_z[i];

        g_0_yyyzz_0_yz_0[i] = 2.0 * g_0_yzz_0_yz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_yz_1[i] * fti_ab_0 + g_0_yyzz_0_z_1[i] * fi_abcd_0 +
                              g_0_yyzz_0_yz_0[i] * pb_y + g_0_yyzz_0_yz_1[i] * wp_y[i];

        g_0_yyyzz_0_zz_0[i] =
            2.0 * g_0_yzz_0_zz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_zz_1[i] * fti_ab_0 + g_0_yyzz_0_zz_0[i] * pb_y + g_0_yyzz_0_zz_1[i] * wp_y[i];
    }

    /// Set up 108-114 components of targeted buffer : SHSD

    auto g_0_yyzzz_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 108);

    auto g_0_yyzzz_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 109);

    auto g_0_yyzzz_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 110);

    auto g_0_yyzzz_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 111);

    auto g_0_yyzzz_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 112);

    auto g_0_yyzzz_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 113);

#pragma omp simd aligned(g_0_yyz_0_xy_0,       \
                             g_0_yyz_0_xy_1,   \
                             g_0_yyz_0_yy_0,   \
                             g_0_yyz_0_yy_1,   \
                             g_0_yyzz_0_xy_0,  \
                             g_0_yyzz_0_xy_1,  \
                             g_0_yyzz_0_yy_0,  \
                             g_0_yyzz_0_yy_1,  \
                             g_0_yyzzz_0_xx_0, \
                             g_0_yyzzz_0_xy_0, \
                             g_0_yyzzz_0_xz_0, \
                             g_0_yyzzz_0_yy_0, \
                             g_0_yyzzz_0_yz_0, \
                             g_0_yyzzz_0_zz_0, \
                             g_0_yzzz_0_xx_0,  \
                             g_0_yzzz_0_xx_1,  \
                             g_0_yzzz_0_xz_0,  \
                             g_0_yzzz_0_xz_1,  \
                             g_0_yzzz_0_yz_0,  \
                             g_0_yzzz_0_yz_1,  \
                             g_0_yzzz_0_z_1,   \
                             g_0_yzzz_0_zz_0,  \
                             g_0_yzzz_0_zz_1,  \
                             g_0_zzz_0_xx_0,   \
                             g_0_zzz_0_xx_1,   \
                             g_0_zzz_0_xz_0,   \
                             g_0_zzz_0_xz_1,   \
                             g_0_zzz_0_yz_0,   \
                             g_0_zzz_0_yz_1,   \
                             g_0_zzz_0_zz_0,   \
                             g_0_zzz_0_zz_1,   \
                             wp_y,             \
                             wp_z,             \
                             c_exps,           \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyzzz_0_xx_0[i] = g_0_zzz_0_xx_0[i] * fi_ab_0 - g_0_zzz_0_xx_1[i] * fti_ab_0 + g_0_yzzz_0_xx_0[i] * pb_y + g_0_yzzz_0_xx_1[i] * wp_y[i];

        g_0_yyzzz_0_xy_0[i] =
            2.0 * g_0_yyz_0_xy_0[i] * fi_ab_0 - 2.0 * g_0_yyz_0_xy_1[i] * fti_ab_0 + g_0_yyzz_0_xy_0[i] * pb_z + g_0_yyzz_0_xy_1[i] * wp_z[i];

        g_0_yyzzz_0_xz_0[i] = g_0_zzz_0_xz_0[i] * fi_ab_0 - g_0_zzz_0_xz_1[i] * fti_ab_0 + g_0_yzzz_0_xz_0[i] * pb_y + g_0_yzzz_0_xz_1[i] * wp_y[i];

        g_0_yyzzz_0_yy_0[i] =
            2.0 * g_0_yyz_0_yy_0[i] * fi_ab_0 - 2.0 * g_0_yyz_0_yy_1[i] * fti_ab_0 + g_0_yyzz_0_yy_0[i] * pb_z + g_0_yyzz_0_yy_1[i] * wp_z[i];

        g_0_yyzzz_0_yz_0[i] = g_0_zzz_0_yz_0[i] * fi_ab_0 - g_0_zzz_0_yz_1[i] * fti_ab_0 + g_0_yzzz_0_z_1[i] * fi_abcd_0 + g_0_yzzz_0_yz_0[i] * pb_y +
                              g_0_yzzz_0_yz_1[i] * wp_y[i];

        g_0_yyzzz_0_zz_0[i] = g_0_zzz_0_zz_0[i] * fi_ab_0 - g_0_zzz_0_zz_1[i] * fti_ab_0 + g_0_yzzz_0_zz_0[i] * pb_y + g_0_yzzz_0_zz_1[i] * wp_y[i];
    }

    /// Set up 114-120 components of targeted buffer : SHSD

    auto g_0_yzzzz_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 114);

    auto g_0_yzzzz_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 115);

    auto g_0_yzzzz_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 116);

    auto g_0_yzzzz_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 117);

    auto g_0_yzzzz_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 118);

    auto g_0_yzzzz_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 119);

#pragma omp simd aligned(g_0_yzzzz_0_xx_0,     \
                             g_0_yzzzz_0_xy_0, \
                             g_0_yzzzz_0_xz_0, \
                             g_0_yzzzz_0_yy_0, \
                             g_0_yzzzz_0_yz_0, \
                             g_0_yzzzz_0_zz_0, \
                             g_0_zzzz_0_x_1,   \
                             g_0_zzzz_0_xx_0,  \
                             g_0_zzzz_0_xx_1,  \
                             g_0_zzzz_0_xy_0,  \
                             g_0_zzzz_0_xy_1,  \
                             g_0_zzzz_0_xz_0,  \
                             g_0_zzzz_0_xz_1,  \
                             g_0_zzzz_0_y_1,   \
                             g_0_zzzz_0_yy_0,  \
                             g_0_zzzz_0_yy_1,  \
                             g_0_zzzz_0_yz_0,  \
                             g_0_zzzz_0_yz_1,  \
                             g_0_zzzz_0_z_1,   \
                             g_0_zzzz_0_zz_0,  \
                             g_0_zzzz_0_zz_1,  \
                             wp_y,             \
                             c_exps,           \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzzzz_0_xx_0[i] = g_0_zzzz_0_xx_0[i] * pb_y + g_0_zzzz_0_xx_1[i] * wp_y[i];

        g_0_yzzzz_0_xy_0[i] = g_0_zzzz_0_x_1[i] * fi_abcd_0 + g_0_zzzz_0_xy_0[i] * pb_y + g_0_zzzz_0_xy_1[i] * wp_y[i];

        g_0_yzzzz_0_xz_0[i] = g_0_zzzz_0_xz_0[i] * pb_y + g_0_zzzz_0_xz_1[i] * wp_y[i];

        g_0_yzzzz_0_yy_0[i] = 2.0 * g_0_zzzz_0_y_1[i] * fi_abcd_0 + g_0_zzzz_0_yy_0[i] * pb_y + g_0_zzzz_0_yy_1[i] * wp_y[i];

        g_0_yzzzz_0_yz_0[i] = g_0_zzzz_0_z_1[i] * fi_abcd_0 + g_0_zzzz_0_yz_0[i] * pb_y + g_0_zzzz_0_yz_1[i] * wp_y[i];

        g_0_yzzzz_0_zz_0[i] = g_0_zzzz_0_zz_0[i] * pb_y + g_0_zzzz_0_zz_1[i] * wp_y[i];
    }

    /// Set up 120-126 components of targeted buffer : SHSD

    auto g_0_zzzzz_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 120);

    auto g_0_zzzzz_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 121);

    auto g_0_zzzzz_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 122);

    auto g_0_zzzzz_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 123);

    auto g_0_zzzzz_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 124);

    auto g_0_zzzzz_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 125);

#pragma omp simd aligned(g_0_zzz_0_xx_0,       \
                             g_0_zzz_0_xx_1,   \
                             g_0_zzz_0_xy_0,   \
                             g_0_zzz_0_xy_1,   \
                             g_0_zzz_0_xz_0,   \
                             g_0_zzz_0_xz_1,   \
                             g_0_zzz_0_yy_0,   \
                             g_0_zzz_0_yy_1,   \
                             g_0_zzz_0_yz_0,   \
                             g_0_zzz_0_yz_1,   \
                             g_0_zzz_0_zz_0,   \
                             g_0_zzz_0_zz_1,   \
                             g_0_zzzz_0_x_1,   \
                             g_0_zzzz_0_xx_0,  \
                             g_0_zzzz_0_xx_1,  \
                             g_0_zzzz_0_xy_0,  \
                             g_0_zzzz_0_xy_1,  \
                             g_0_zzzz_0_xz_0,  \
                             g_0_zzzz_0_xz_1,  \
                             g_0_zzzz_0_y_1,   \
                             g_0_zzzz_0_yy_0,  \
                             g_0_zzzz_0_yy_1,  \
                             g_0_zzzz_0_yz_0,  \
                             g_0_zzzz_0_yz_1,  \
                             g_0_zzzz_0_z_1,   \
                             g_0_zzzz_0_zz_0,  \
                             g_0_zzzz_0_zz_1,  \
                             g_0_zzzzz_0_xx_0, \
                             g_0_zzzzz_0_xy_0, \
                             g_0_zzzzz_0_xz_0, \
                             g_0_zzzzz_0_yy_0, \
                             g_0_zzzzz_0_yz_0, \
                             g_0_zzzzz_0_zz_0, \
                             wp_z,             \
                             c_exps,           \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzzzz_0_xx_0[i] =
            4.0 * g_0_zzz_0_xx_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xx_1[i] * fti_ab_0 + g_0_zzzz_0_xx_0[i] * pb_z + g_0_zzzz_0_xx_1[i] * wp_z[i];

        g_0_zzzzz_0_xy_0[i] =
            4.0 * g_0_zzz_0_xy_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xy_1[i] * fti_ab_0 + g_0_zzzz_0_xy_0[i] * pb_z + g_0_zzzz_0_xy_1[i] * wp_z[i];

        g_0_zzzzz_0_xz_0[i] = 4.0 * g_0_zzz_0_xz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xz_1[i] * fti_ab_0 + g_0_zzzz_0_x_1[i] * fi_abcd_0 +
                              g_0_zzzz_0_xz_0[i] * pb_z + g_0_zzzz_0_xz_1[i] * wp_z[i];

        g_0_zzzzz_0_yy_0[i] =
            4.0 * g_0_zzz_0_yy_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_yy_1[i] * fti_ab_0 + g_0_zzzz_0_yy_0[i] * pb_z + g_0_zzzz_0_yy_1[i] * wp_z[i];

        g_0_zzzzz_0_yz_0[i] = 4.0 * g_0_zzz_0_yz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_yz_1[i] * fti_ab_0 + g_0_zzzz_0_y_1[i] * fi_abcd_0 +
                              g_0_zzzz_0_yz_0[i] * pb_z + g_0_zzzz_0_yz_1[i] * wp_z[i];

        g_0_zzzzz_0_zz_0[i] = 4.0 * g_0_zzz_0_zz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_zz_1[i] * fti_ab_0 + 2.0 * g_0_zzzz_0_z_1[i] * fi_abcd_0 +
                              g_0_zzzz_0_zz_0[i] * pb_z + g_0_zzzz_0_zz_1[i] * wp_z[i];
    }
}

}  // namespace erirec

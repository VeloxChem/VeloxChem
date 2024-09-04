#include "ElectronRepulsionPrimRecSISD.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_sisd(CSimdArray<double>& pbuffer,
                                  const size_t idx_eri_0_sisd,
                                  size_t idx_eri_0_sgsd,
                                  size_t idx_eri_1_sgsd,
                                  size_t idx_eri_1_shsp,
                                  size_t idx_eri_0_shsd,
                                  size_t idx_eri_1_shsd,
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

    /// Set up components of auxilary buffer : SGSD

    auto g_0_xxxx_0_xx_0 = pbuffer.data(idx_eri_0_sgsd);

    auto g_0_xxxx_0_xy_0 = pbuffer.data(idx_eri_0_sgsd + 1);

    auto g_0_xxxx_0_xz_0 = pbuffer.data(idx_eri_0_sgsd + 2);

    auto g_0_xxxx_0_yy_0 = pbuffer.data(idx_eri_0_sgsd + 3);

    auto g_0_xxxx_0_yz_0 = pbuffer.data(idx_eri_0_sgsd + 4);

    auto g_0_xxxx_0_zz_0 = pbuffer.data(idx_eri_0_sgsd + 5);

    auto g_0_xxxy_0_xx_0 = pbuffer.data(idx_eri_0_sgsd + 6);

    auto g_0_xxxy_0_xz_0 = pbuffer.data(idx_eri_0_sgsd + 8);

    auto g_0_xxxz_0_xx_0 = pbuffer.data(idx_eri_0_sgsd + 12);

    auto g_0_xxxz_0_xy_0 = pbuffer.data(idx_eri_0_sgsd + 13);

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

    auto g_0_xyyy_0_xy_0 = pbuffer.data(idx_eri_0_sgsd + 37);

    auto g_0_xyyy_0_yy_0 = pbuffer.data(idx_eri_0_sgsd + 39);

    auto g_0_xyyy_0_yz_0 = pbuffer.data(idx_eri_0_sgsd + 40);

    auto g_0_xyyy_0_zz_0 = pbuffer.data(idx_eri_0_sgsd + 41);

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

    auto g_0_yyyz_0_yy_0 = pbuffer.data(idx_eri_0_sgsd + 69);

    auto g_0_yyzz_0_xx_0 = pbuffer.data(idx_eri_0_sgsd + 72);

    auto g_0_yyzz_0_xy_0 = pbuffer.data(idx_eri_0_sgsd + 73);

    auto g_0_yyzz_0_xz_0 = pbuffer.data(idx_eri_0_sgsd + 74);

    auto g_0_yyzz_0_yy_0 = pbuffer.data(idx_eri_0_sgsd + 75);

    auto g_0_yyzz_0_yz_0 = pbuffer.data(idx_eri_0_sgsd + 76);

    auto g_0_yyzz_0_zz_0 = pbuffer.data(idx_eri_0_sgsd + 77);

    auto g_0_yzzz_0_xx_0 = pbuffer.data(idx_eri_0_sgsd + 78);

    auto g_0_yzzz_0_xz_0 = pbuffer.data(idx_eri_0_sgsd + 80);

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

    auto g_0_xxxy_0_xz_1 = pbuffer.data(idx_eri_1_sgsd + 8);

    auto g_0_xxxz_0_xx_1 = pbuffer.data(idx_eri_1_sgsd + 12);

    auto g_0_xxxz_0_xy_1 = pbuffer.data(idx_eri_1_sgsd + 13);

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

    auto g_0_xyyy_0_xy_1 = pbuffer.data(idx_eri_1_sgsd + 37);

    auto g_0_xyyy_0_yy_1 = pbuffer.data(idx_eri_1_sgsd + 39);

    auto g_0_xyyy_0_yz_1 = pbuffer.data(idx_eri_1_sgsd + 40);

    auto g_0_xyyy_0_zz_1 = pbuffer.data(idx_eri_1_sgsd + 41);

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

    auto g_0_yyyz_0_yy_1 = pbuffer.data(idx_eri_1_sgsd + 69);

    auto g_0_yyzz_0_xx_1 = pbuffer.data(idx_eri_1_sgsd + 72);

    auto g_0_yyzz_0_xy_1 = pbuffer.data(idx_eri_1_sgsd + 73);

    auto g_0_yyzz_0_xz_1 = pbuffer.data(idx_eri_1_sgsd + 74);

    auto g_0_yyzz_0_yy_1 = pbuffer.data(idx_eri_1_sgsd + 75);

    auto g_0_yyzz_0_yz_1 = pbuffer.data(idx_eri_1_sgsd + 76);

    auto g_0_yyzz_0_zz_1 = pbuffer.data(idx_eri_1_sgsd + 77);

    auto g_0_yzzz_0_xx_1 = pbuffer.data(idx_eri_1_sgsd + 78);

    auto g_0_yzzz_0_xz_1 = pbuffer.data(idx_eri_1_sgsd + 80);

    auto g_0_yzzz_0_yz_1 = pbuffer.data(idx_eri_1_sgsd + 82);

    auto g_0_yzzz_0_zz_1 = pbuffer.data(idx_eri_1_sgsd + 83);

    auto g_0_zzzz_0_xx_1 = pbuffer.data(idx_eri_1_sgsd + 84);

    auto g_0_zzzz_0_xy_1 = pbuffer.data(idx_eri_1_sgsd + 85);

    auto g_0_zzzz_0_xz_1 = pbuffer.data(idx_eri_1_sgsd + 86);

    auto g_0_zzzz_0_yy_1 = pbuffer.data(idx_eri_1_sgsd + 87);

    auto g_0_zzzz_0_yz_1 = pbuffer.data(idx_eri_1_sgsd + 88);

    auto g_0_zzzz_0_zz_1 = pbuffer.data(idx_eri_1_sgsd + 89);

    /// Set up components of auxilary buffer : SHSP

    auto g_0_xxxxx_0_x_1 = pbuffer.data(idx_eri_1_shsp);

    auto g_0_xxxxx_0_y_1 = pbuffer.data(idx_eri_1_shsp + 1);

    auto g_0_xxxxx_0_z_1 = pbuffer.data(idx_eri_1_shsp + 2);

    auto g_0_xxxxz_0_z_1 = pbuffer.data(idx_eri_1_shsp + 8);

    auto g_0_xxxyy_0_x_1 = pbuffer.data(idx_eri_1_shsp + 9);

    auto g_0_xxxyy_0_y_1 = pbuffer.data(idx_eri_1_shsp + 10);

    auto g_0_xxxyy_0_z_1 = pbuffer.data(idx_eri_1_shsp + 11);

    auto g_0_xxxzz_0_x_1 = pbuffer.data(idx_eri_1_shsp + 15);

    auto g_0_xxxzz_0_y_1 = pbuffer.data(idx_eri_1_shsp + 16);

    auto g_0_xxxzz_0_z_1 = pbuffer.data(idx_eri_1_shsp + 17);

    auto g_0_xxyyy_0_x_1 = pbuffer.data(idx_eri_1_shsp + 18);

    auto g_0_xxyyy_0_y_1 = pbuffer.data(idx_eri_1_shsp + 19);

    auto g_0_xxyyy_0_z_1 = pbuffer.data(idx_eri_1_shsp + 20);

    auto g_0_xxzzz_0_x_1 = pbuffer.data(idx_eri_1_shsp + 27);

    auto g_0_xxzzz_0_y_1 = pbuffer.data(idx_eri_1_shsp + 28);

    auto g_0_xxzzz_0_z_1 = pbuffer.data(idx_eri_1_shsp + 29);

    auto g_0_xyyyy_0_y_1 = pbuffer.data(idx_eri_1_shsp + 31);

    auto g_0_xzzzz_0_z_1 = pbuffer.data(idx_eri_1_shsp + 44);

    auto g_0_yyyyy_0_x_1 = pbuffer.data(idx_eri_1_shsp + 45);

    auto g_0_yyyyy_0_y_1 = pbuffer.data(idx_eri_1_shsp + 46);

    auto g_0_yyyyy_0_z_1 = pbuffer.data(idx_eri_1_shsp + 47);

    auto g_0_yyyyz_0_z_1 = pbuffer.data(idx_eri_1_shsp + 50);

    auto g_0_yyyzz_0_x_1 = pbuffer.data(idx_eri_1_shsp + 51);

    auto g_0_yyyzz_0_y_1 = pbuffer.data(idx_eri_1_shsp + 52);

    auto g_0_yyyzz_0_z_1 = pbuffer.data(idx_eri_1_shsp + 53);

    auto g_0_yyzzz_0_x_1 = pbuffer.data(idx_eri_1_shsp + 54);

    auto g_0_yyzzz_0_y_1 = pbuffer.data(idx_eri_1_shsp + 55);

    auto g_0_yyzzz_0_z_1 = pbuffer.data(idx_eri_1_shsp + 56);

    auto g_0_yzzzz_0_y_1 = pbuffer.data(idx_eri_1_shsp + 58);

    auto g_0_yzzzz_0_z_1 = pbuffer.data(idx_eri_1_shsp + 59);

    auto g_0_zzzzz_0_x_1 = pbuffer.data(idx_eri_1_shsp + 60);

    auto g_0_zzzzz_0_y_1 = pbuffer.data(idx_eri_1_shsp + 61);

    auto g_0_zzzzz_0_z_1 = pbuffer.data(idx_eri_1_shsp + 62);

    /// Set up components of auxilary buffer : SHSD

    auto g_0_xxxxx_0_xx_0 = pbuffer.data(idx_eri_0_shsd);

    auto g_0_xxxxx_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 1);

    auto g_0_xxxxx_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 2);

    auto g_0_xxxxx_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 3);

    auto g_0_xxxxx_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 4);

    auto g_0_xxxxx_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 5);

    auto g_0_xxxxy_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 6);

    auto g_0_xxxxy_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 7);

    auto g_0_xxxxy_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 8);

    auto g_0_xxxxy_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 9);

    auto g_0_xxxxz_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 12);

    auto g_0_xxxxz_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 13);

    auto g_0_xxxxz_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 14);

    auto g_0_xxxxz_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 16);

    auto g_0_xxxxz_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 17);

    auto g_0_xxxyy_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 18);

    auto g_0_xxxyy_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 19);

    auto g_0_xxxyy_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 20);

    auto g_0_xxxyy_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 21);

    auto g_0_xxxyy_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 22);

    auto g_0_xxxyy_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 23);

    auto g_0_xxxzz_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 30);

    auto g_0_xxxzz_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 31);

    auto g_0_xxxzz_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 32);

    auto g_0_xxxzz_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 33);

    auto g_0_xxxzz_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 34);

    auto g_0_xxxzz_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 35);

    auto g_0_xxyyy_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 36);

    auto g_0_xxyyy_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 37);

    auto g_0_xxyyy_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 38);

    auto g_0_xxyyy_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 39);

    auto g_0_xxyyy_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 40);

    auto g_0_xxyyy_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 41);

    auto g_0_xxyyz_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 43);

    auto g_0_xxyzz_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 48);

    auto g_0_xxyzz_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 50);

    auto g_0_xxzzz_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 54);

    auto g_0_xxzzz_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 55);

    auto g_0_xxzzz_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 56);

    auto g_0_xxzzz_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 57);

    auto g_0_xxzzz_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 58);

    auto g_0_xxzzz_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 59);

    auto g_0_xyyyy_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 60);

    auto g_0_xyyyy_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 61);

    auto g_0_xyyyy_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 63);

    auto g_0_xyyyy_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 64);

    auto g_0_xyyyy_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 65);

    auto g_0_xyyzz_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 75);

    auto g_0_xyyzz_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 76);

    auto g_0_xyyzz_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 77);

    auto g_0_xzzzz_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 84);

    auto g_0_xzzzz_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 86);

    auto g_0_xzzzz_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 87);

    auto g_0_xzzzz_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 88);

    auto g_0_xzzzz_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 89);

    auto g_0_yyyyy_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 90);

    auto g_0_yyyyy_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 91);

    auto g_0_yyyyy_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 92);

    auto g_0_yyyyy_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 93);

    auto g_0_yyyyy_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 94);

    auto g_0_yyyyy_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 95);

    auto g_0_yyyyz_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 97);

    auto g_0_yyyyz_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 98);

    auto g_0_yyyyz_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 99);

    auto g_0_yyyyz_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 100);

    auto g_0_yyyyz_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 101);

    auto g_0_yyyzz_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 102);

    auto g_0_yyyzz_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 103);

    auto g_0_yyyzz_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 104);

    auto g_0_yyyzz_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 105);

    auto g_0_yyyzz_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 106);

    auto g_0_yyyzz_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 107);

    auto g_0_yyzzz_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 108);

    auto g_0_yyzzz_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 109);

    auto g_0_yyzzz_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 110);

    auto g_0_yyzzz_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 111);

    auto g_0_yyzzz_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 112);

    auto g_0_yyzzz_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 113);

    auto g_0_yzzzz_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 114);

    auto g_0_yzzzz_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 115);

    auto g_0_yzzzz_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 116);

    auto g_0_yzzzz_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 117);

    auto g_0_yzzzz_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 118);

    auto g_0_yzzzz_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 119);

    auto g_0_zzzzz_0_xx_0 = pbuffer.data(idx_eri_0_shsd + 120);

    auto g_0_zzzzz_0_xy_0 = pbuffer.data(idx_eri_0_shsd + 121);

    auto g_0_zzzzz_0_xz_0 = pbuffer.data(idx_eri_0_shsd + 122);

    auto g_0_zzzzz_0_yy_0 = pbuffer.data(idx_eri_0_shsd + 123);

    auto g_0_zzzzz_0_yz_0 = pbuffer.data(idx_eri_0_shsd + 124);

    auto g_0_zzzzz_0_zz_0 = pbuffer.data(idx_eri_0_shsd + 125);

    /// Set up components of auxilary buffer : SHSD

    auto g_0_xxxxx_0_xx_1 = pbuffer.data(idx_eri_1_shsd);

    auto g_0_xxxxx_0_xy_1 = pbuffer.data(idx_eri_1_shsd + 1);

    auto g_0_xxxxx_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 2);

    auto g_0_xxxxx_0_yy_1 = pbuffer.data(idx_eri_1_shsd + 3);

    auto g_0_xxxxx_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 4);

    auto g_0_xxxxx_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 5);

    auto g_0_xxxxy_0_xx_1 = pbuffer.data(idx_eri_1_shsd + 6);

    auto g_0_xxxxy_0_xy_1 = pbuffer.data(idx_eri_1_shsd + 7);

    auto g_0_xxxxy_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 8);

    auto g_0_xxxxy_0_yy_1 = pbuffer.data(idx_eri_1_shsd + 9);

    auto g_0_xxxxz_0_xx_1 = pbuffer.data(idx_eri_1_shsd + 12);

    auto g_0_xxxxz_0_xy_1 = pbuffer.data(idx_eri_1_shsd + 13);

    auto g_0_xxxxz_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 14);

    auto g_0_xxxxz_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 16);

    auto g_0_xxxxz_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 17);

    auto g_0_xxxyy_0_xx_1 = pbuffer.data(idx_eri_1_shsd + 18);

    auto g_0_xxxyy_0_xy_1 = pbuffer.data(idx_eri_1_shsd + 19);

    auto g_0_xxxyy_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 20);

    auto g_0_xxxyy_0_yy_1 = pbuffer.data(idx_eri_1_shsd + 21);

    auto g_0_xxxyy_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 22);

    auto g_0_xxxyy_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 23);

    auto g_0_xxxzz_0_xx_1 = pbuffer.data(idx_eri_1_shsd + 30);

    auto g_0_xxxzz_0_xy_1 = pbuffer.data(idx_eri_1_shsd + 31);

    auto g_0_xxxzz_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 32);

    auto g_0_xxxzz_0_yy_1 = pbuffer.data(idx_eri_1_shsd + 33);

    auto g_0_xxxzz_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 34);

    auto g_0_xxxzz_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 35);

    auto g_0_xxyyy_0_xx_1 = pbuffer.data(idx_eri_1_shsd + 36);

    auto g_0_xxyyy_0_xy_1 = pbuffer.data(idx_eri_1_shsd + 37);

    auto g_0_xxyyy_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 38);

    auto g_0_xxyyy_0_yy_1 = pbuffer.data(idx_eri_1_shsd + 39);

    auto g_0_xxyyy_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 40);

    auto g_0_xxyyy_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 41);

    auto g_0_xxyyz_0_xy_1 = pbuffer.data(idx_eri_1_shsd + 43);

    auto g_0_xxyzz_0_xx_1 = pbuffer.data(idx_eri_1_shsd + 48);

    auto g_0_xxyzz_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 50);

    auto g_0_xxzzz_0_xx_1 = pbuffer.data(idx_eri_1_shsd + 54);

    auto g_0_xxzzz_0_xy_1 = pbuffer.data(idx_eri_1_shsd + 55);

    auto g_0_xxzzz_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 56);

    auto g_0_xxzzz_0_yy_1 = pbuffer.data(idx_eri_1_shsd + 57);

    auto g_0_xxzzz_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 58);

    auto g_0_xxzzz_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 59);

    auto g_0_xyyyy_0_xx_1 = pbuffer.data(idx_eri_1_shsd + 60);

    auto g_0_xyyyy_0_xy_1 = pbuffer.data(idx_eri_1_shsd + 61);

    auto g_0_xyyyy_0_yy_1 = pbuffer.data(idx_eri_1_shsd + 63);

    auto g_0_xyyyy_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 64);

    auto g_0_xyyyy_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 65);

    auto g_0_xyyzz_0_yy_1 = pbuffer.data(idx_eri_1_shsd + 75);

    auto g_0_xyyzz_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 76);

    auto g_0_xyyzz_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 77);

    auto g_0_xzzzz_0_xx_1 = pbuffer.data(idx_eri_1_shsd + 84);

    auto g_0_xzzzz_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 86);

    auto g_0_xzzzz_0_yy_1 = pbuffer.data(idx_eri_1_shsd + 87);

    auto g_0_xzzzz_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 88);

    auto g_0_xzzzz_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 89);

    auto g_0_yyyyy_0_xx_1 = pbuffer.data(idx_eri_1_shsd + 90);

    auto g_0_yyyyy_0_xy_1 = pbuffer.data(idx_eri_1_shsd + 91);

    auto g_0_yyyyy_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 92);

    auto g_0_yyyyy_0_yy_1 = pbuffer.data(idx_eri_1_shsd + 93);

    auto g_0_yyyyy_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 94);

    auto g_0_yyyyy_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 95);

    auto g_0_yyyyz_0_xy_1 = pbuffer.data(idx_eri_1_shsd + 97);

    auto g_0_yyyyz_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 98);

    auto g_0_yyyyz_0_yy_1 = pbuffer.data(idx_eri_1_shsd + 99);

    auto g_0_yyyyz_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 100);

    auto g_0_yyyyz_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 101);

    auto g_0_yyyzz_0_xx_1 = pbuffer.data(idx_eri_1_shsd + 102);

    auto g_0_yyyzz_0_xy_1 = pbuffer.data(idx_eri_1_shsd + 103);

    auto g_0_yyyzz_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 104);

    auto g_0_yyyzz_0_yy_1 = pbuffer.data(idx_eri_1_shsd + 105);

    auto g_0_yyyzz_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 106);

    auto g_0_yyyzz_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 107);

    auto g_0_yyzzz_0_xx_1 = pbuffer.data(idx_eri_1_shsd + 108);

    auto g_0_yyzzz_0_xy_1 = pbuffer.data(idx_eri_1_shsd + 109);

    auto g_0_yyzzz_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 110);

    auto g_0_yyzzz_0_yy_1 = pbuffer.data(idx_eri_1_shsd + 111);

    auto g_0_yyzzz_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 112);

    auto g_0_yyzzz_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 113);

    auto g_0_yzzzz_0_xx_1 = pbuffer.data(idx_eri_1_shsd + 114);

    auto g_0_yzzzz_0_xy_1 = pbuffer.data(idx_eri_1_shsd + 115);

    auto g_0_yzzzz_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 116);

    auto g_0_yzzzz_0_yy_1 = pbuffer.data(idx_eri_1_shsd + 117);

    auto g_0_yzzzz_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 118);

    auto g_0_yzzzz_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 119);

    auto g_0_zzzzz_0_xx_1 = pbuffer.data(idx_eri_1_shsd + 120);

    auto g_0_zzzzz_0_xy_1 = pbuffer.data(idx_eri_1_shsd + 121);

    auto g_0_zzzzz_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 122);

    auto g_0_zzzzz_0_yy_1 = pbuffer.data(idx_eri_1_shsd + 123);

    auto g_0_zzzzz_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 124);

    auto g_0_zzzzz_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 125);

    /// Set up 0-6 components of targeted buffer : SISD

    auto g_0_xxxxxx_0_xx_0 = pbuffer.data(idx_eri_0_sisd);

    auto g_0_xxxxxx_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 1);

    auto g_0_xxxxxx_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 2);

    auto g_0_xxxxxx_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 3);

    auto g_0_xxxxxx_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 4);

    auto g_0_xxxxxx_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 5);

    #pragma omp simd aligned(g_0_xxxx_0_xx_0, g_0_xxxx_0_xx_1, g_0_xxxx_0_xy_0, g_0_xxxx_0_xy_1, g_0_xxxx_0_xz_0, g_0_xxxx_0_xz_1, g_0_xxxx_0_yy_0, g_0_xxxx_0_yy_1, g_0_xxxx_0_yz_0, g_0_xxxx_0_yz_1, g_0_xxxx_0_zz_0, g_0_xxxx_0_zz_1, g_0_xxxxx_0_x_1, g_0_xxxxx_0_xx_0, g_0_xxxxx_0_xx_1, g_0_xxxxx_0_xy_0, g_0_xxxxx_0_xy_1, g_0_xxxxx_0_xz_0, g_0_xxxxx_0_xz_1, g_0_xxxxx_0_y_1, g_0_xxxxx_0_yy_0, g_0_xxxxx_0_yy_1, g_0_xxxxx_0_yz_0, g_0_xxxxx_0_yz_1, g_0_xxxxx_0_z_1, g_0_xxxxx_0_zz_0, g_0_xxxxx_0_zz_1, g_0_xxxxxx_0_xx_0, g_0_xxxxxx_0_xy_0, g_0_xxxxxx_0_xz_0, g_0_xxxxxx_0_yy_0, g_0_xxxxxx_0_yz_0, g_0_xxxxxx_0_zz_0, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxx_0_xx_0[i] = 5.0 * g_0_xxxx_0_xx_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xx_1[i] * fti_ab_0 + 2.0 * g_0_xxxxx_0_x_1[i] * fi_abcd_0 + g_0_xxxxx_0_xx_0[i] * pb_x + g_0_xxxxx_0_xx_1[i] * wp_x[i];

        g_0_xxxxxx_0_xy_0[i] = 5.0 * g_0_xxxx_0_xy_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xy_1[i] * fti_ab_0 + g_0_xxxxx_0_y_1[i] * fi_abcd_0 + g_0_xxxxx_0_xy_0[i] * pb_x + g_0_xxxxx_0_xy_1[i] * wp_x[i];

        g_0_xxxxxx_0_xz_0[i] = 5.0 * g_0_xxxx_0_xz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xz_1[i] * fti_ab_0 + g_0_xxxxx_0_z_1[i] * fi_abcd_0 + g_0_xxxxx_0_xz_0[i] * pb_x + g_0_xxxxx_0_xz_1[i] * wp_x[i];

        g_0_xxxxxx_0_yy_0[i] = 5.0 * g_0_xxxx_0_yy_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_yy_1[i] * fti_ab_0 + g_0_xxxxx_0_yy_0[i] * pb_x + g_0_xxxxx_0_yy_1[i] * wp_x[i];

        g_0_xxxxxx_0_yz_0[i] = 5.0 * g_0_xxxx_0_yz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_yz_1[i] * fti_ab_0 + g_0_xxxxx_0_yz_0[i] * pb_x + g_0_xxxxx_0_yz_1[i] * wp_x[i];

        g_0_xxxxxx_0_zz_0[i] = 5.0 * g_0_xxxx_0_zz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_zz_1[i] * fti_ab_0 + g_0_xxxxx_0_zz_0[i] * pb_x + g_0_xxxxx_0_zz_1[i] * wp_x[i];
    }

    /// Set up 6-12 components of targeted buffer : SISD

    auto g_0_xxxxxy_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 6);

    auto g_0_xxxxxy_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 7);

    auto g_0_xxxxxy_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 8);

    auto g_0_xxxxxy_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 9);

    auto g_0_xxxxxy_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 10);

    auto g_0_xxxxxy_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 11);

    #pragma omp simd aligned(g_0_xxxxx_0_x_1, g_0_xxxxx_0_xx_0, g_0_xxxxx_0_xx_1, g_0_xxxxx_0_xy_0, g_0_xxxxx_0_xy_1, g_0_xxxxx_0_xz_0, g_0_xxxxx_0_xz_1, g_0_xxxxx_0_y_1, g_0_xxxxx_0_yy_0, g_0_xxxxx_0_yy_1, g_0_xxxxx_0_yz_0, g_0_xxxxx_0_yz_1, g_0_xxxxx_0_z_1, g_0_xxxxx_0_zz_0, g_0_xxxxx_0_zz_1, g_0_xxxxxy_0_xx_0, g_0_xxxxxy_0_xy_0, g_0_xxxxxy_0_xz_0, g_0_xxxxxy_0_yy_0, g_0_xxxxxy_0_yz_0, g_0_xxxxxy_0_zz_0, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxy_0_xx_0[i] = g_0_xxxxx_0_xx_0[i] * pb_y + g_0_xxxxx_0_xx_1[i] * wp_y[i];

        g_0_xxxxxy_0_xy_0[i] = g_0_xxxxx_0_x_1[i] * fi_abcd_0 + g_0_xxxxx_0_xy_0[i] * pb_y + g_0_xxxxx_0_xy_1[i] * wp_y[i];

        g_0_xxxxxy_0_xz_0[i] = g_0_xxxxx_0_xz_0[i] * pb_y + g_0_xxxxx_0_xz_1[i] * wp_y[i];

        g_0_xxxxxy_0_yy_0[i] = 2.0 * g_0_xxxxx_0_y_1[i] * fi_abcd_0 + g_0_xxxxx_0_yy_0[i] * pb_y + g_0_xxxxx_0_yy_1[i] * wp_y[i];

        g_0_xxxxxy_0_yz_0[i] = g_0_xxxxx_0_z_1[i] * fi_abcd_0 + g_0_xxxxx_0_yz_0[i] * pb_y + g_0_xxxxx_0_yz_1[i] * wp_y[i];

        g_0_xxxxxy_0_zz_0[i] = g_0_xxxxx_0_zz_0[i] * pb_y + g_0_xxxxx_0_zz_1[i] * wp_y[i];
    }

    /// Set up 12-18 components of targeted buffer : SISD

    auto g_0_xxxxxz_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 12);

    auto g_0_xxxxxz_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 13);

    auto g_0_xxxxxz_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 14);

    auto g_0_xxxxxz_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 15);

    auto g_0_xxxxxz_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 16);

    auto g_0_xxxxxz_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 17);

    #pragma omp simd aligned(g_0_xxxxx_0_x_1, g_0_xxxxx_0_xx_0, g_0_xxxxx_0_xx_1, g_0_xxxxx_0_xy_0, g_0_xxxxx_0_xy_1, g_0_xxxxx_0_xz_0, g_0_xxxxx_0_xz_1, g_0_xxxxx_0_y_1, g_0_xxxxx_0_yy_0, g_0_xxxxx_0_yy_1, g_0_xxxxx_0_yz_0, g_0_xxxxx_0_yz_1, g_0_xxxxx_0_z_1, g_0_xxxxx_0_zz_0, g_0_xxxxx_0_zz_1, g_0_xxxxxz_0_xx_0, g_0_xxxxxz_0_xy_0, g_0_xxxxxz_0_xz_0, g_0_xxxxxz_0_yy_0, g_0_xxxxxz_0_yz_0, g_0_xxxxxz_0_zz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxz_0_xx_0[i] = g_0_xxxxx_0_xx_0[i] * pb_z + g_0_xxxxx_0_xx_1[i] * wp_z[i];

        g_0_xxxxxz_0_xy_0[i] = g_0_xxxxx_0_xy_0[i] * pb_z + g_0_xxxxx_0_xy_1[i] * wp_z[i];

        g_0_xxxxxz_0_xz_0[i] = g_0_xxxxx_0_x_1[i] * fi_abcd_0 + g_0_xxxxx_0_xz_0[i] * pb_z + g_0_xxxxx_0_xz_1[i] * wp_z[i];

        g_0_xxxxxz_0_yy_0[i] = g_0_xxxxx_0_yy_0[i] * pb_z + g_0_xxxxx_0_yy_1[i] * wp_z[i];

        g_0_xxxxxz_0_yz_0[i] = g_0_xxxxx_0_y_1[i] * fi_abcd_0 + g_0_xxxxx_0_yz_0[i] * pb_z + g_0_xxxxx_0_yz_1[i] * wp_z[i];

        g_0_xxxxxz_0_zz_0[i] = 2.0 * g_0_xxxxx_0_z_1[i] * fi_abcd_0 + g_0_xxxxx_0_zz_0[i] * pb_z + g_0_xxxxx_0_zz_1[i] * wp_z[i];
    }

    /// Set up 18-24 components of targeted buffer : SISD

    auto g_0_xxxxyy_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 18);

    auto g_0_xxxxyy_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 19);

    auto g_0_xxxxyy_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 20);

    auto g_0_xxxxyy_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 21);

    auto g_0_xxxxyy_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 22);

    auto g_0_xxxxyy_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 23);

    #pragma omp simd aligned(g_0_xxxx_0_xx_0, g_0_xxxx_0_xx_1, g_0_xxxx_0_xz_0, g_0_xxxx_0_xz_1, g_0_xxxxy_0_xx_0, g_0_xxxxy_0_xx_1, g_0_xxxxy_0_xz_0, g_0_xxxxy_0_xz_1, g_0_xxxxyy_0_xx_0, g_0_xxxxyy_0_xy_0, g_0_xxxxyy_0_xz_0, g_0_xxxxyy_0_yy_0, g_0_xxxxyy_0_yz_0, g_0_xxxxyy_0_zz_0, g_0_xxxyy_0_xy_0, g_0_xxxyy_0_xy_1, g_0_xxxyy_0_y_1, g_0_xxxyy_0_yy_0, g_0_xxxyy_0_yy_1, g_0_xxxyy_0_yz_0, g_0_xxxyy_0_yz_1, g_0_xxxyy_0_zz_0, g_0_xxxyy_0_zz_1, g_0_xxyy_0_xy_0, g_0_xxyy_0_xy_1, g_0_xxyy_0_yy_0, g_0_xxyy_0_yy_1, g_0_xxyy_0_yz_0, g_0_xxyy_0_yz_1, g_0_xxyy_0_zz_0, g_0_xxyy_0_zz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxyy_0_xx_0[i] = g_0_xxxx_0_xx_0[i] * fi_ab_0 - g_0_xxxx_0_xx_1[i] * fti_ab_0 + g_0_xxxxy_0_xx_0[i] * pb_y + g_0_xxxxy_0_xx_1[i] * wp_y[i];

        g_0_xxxxyy_0_xy_0[i] = 3.0 * g_0_xxyy_0_xy_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xy_1[i] * fti_ab_0 + g_0_xxxyy_0_y_1[i] * fi_abcd_0 + g_0_xxxyy_0_xy_0[i] * pb_x + g_0_xxxyy_0_xy_1[i] * wp_x[i];

        g_0_xxxxyy_0_xz_0[i] = g_0_xxxx_0_xz_0[i] * fi_ab_0 - g_0_xxxx_0_xz_1[i] * fti_ab_0 + g_0_xxxxy_0_xz_0[i] * pb_y + g_0_xxxxy_0_xz_1[i] * wp_y[i];

        g_0_xxxxyy_0_yy_0[i] = 3.0 * g_0_xxyy_0_yy_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_yy_1[i] * fti_ab_0 + g_0_xxxyy_0_yy_0[i] * pb_x + g_0_xxxyy_0_yy_1[i] * wp_x[i];

        g_0_xxxxyy_0_yz_0[i] = 3.0 * g_0_xxyy_0_yz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_yz_1[i] * fti_ab_0 + g_0_xxxyy_0_yz_0[i] * pb_x + g_0_xxxyy_0_yz_1[i] * wp_x[i];

        g_0_xxxxyy_0_zz_0[i] = 3.0 * g_0_xxyy_0_zz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_zz_1[i] * fti_ab_0 + g_0_xxxyy_0_zz_0[i] * pb_x + g_0_xxxyy_0_zz_1[i] * wp_x[i];
    }

    /// Set up 24-30 components of targeted buffer : SISD

    auto g_0_xxxxyz_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 24);

    auto g_0_xxxxyz_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 25);

    auto g_0_xxxxyz_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 26);

    auto g_0_xxxxyz_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 27);

    auto g_0_xxxxyz_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 28);

    auto g_0_xxxxyz_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 29);

    #pragma omp simd aligned(g_0_xxxxy_0_xy_0, g_0_xxxxy_0_xy_1, g_0_xxxxy_0_yy_0, g_0_xxxxy_0_yy_1, g_0_xxxxyz_0_xx_0, g_0_xxxxyz_0_xy_0, g_0_xxxxyz_0_xz_0, g_0_xxxxyz_0_yy_0, g_0_xxxxyz_0_yz_0, g_0_xxxxyz_0_zz_0, g_0_xxxxz_0_xx_0, g_0_xxxxz_0_xx_1, g_0_xxxxz_0_xz_0, g_0_xxxxz_0_xz_1, g_0_xxxxz_0_yz_0, g_0_xxxxz_0_yz_1, g_0_xxxxz_0_z_1, g_0_xxxxz_0_zz_0, g_0_xxxxz_0_zz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxyz_0_xx_0[i] = g_0_xxxxz_0_xx_0[i] * pb_y + g_0_xxxxz_0_xx_1[i] * wp_y[i];

        g_0_xxxxyz_0_xy_0[i] = g_0_xxxxy_0_xy_0[i] * pb_z + g_0_xxxxy_0_xy_1[i] * wp_z[i];

        g_0_xxxxyz_0_xz_0[i] = g_0_xxxxz_0_xz_0[i] * pb_y + g_0_xxxxz_0_xz_1[i] * wp_y[i];

        g_0_xxxxyz_0_yy_0[i] = g_0_xxxxy_0_yy_0[i] * pb_z + g_0_xxxxy_0_yy_1[i] * wp_z[i];

        g_0_xxxxyz_0_yz_0[i] = g_0_xxxxz_0_z_1[i] * fi_abcd_0 + g_0_xxxxz_0_yz_0[i] * pb_y + g_0_xxxxz_0_yz_1[i] * wp_y[i];

        g_0_xxxxyz_0_zz_0[i] = g_0_xxxxz_0_zz_0[i] * pb_y + g_0_xxxxz_0_zz_1[i] * wp_y[i];
    }

    /// Set up 30-36 components of targeted buffer : SISD

    auto g_0_xxxxzz_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 30);

    auto g_0_xxxxzz_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 31);

    auto g_0_xxxxzz_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 32);

    auto g_0_xxxxzz_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 33);

    auto g_0_xxxxzz_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 34);

    auto g_0_xxxxzz_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 35);

    #pragma omp simd aligned(g_0_xxxx_0_xx_0, g_0_xxxx_0_xx_1, g_0_xxxx_0_xy_0, g_0_xxxx_0_xy_1, g_0_xxxxz_0_xx_0, g_0_xxxxz_0_xx_1, g_0_xxxxz_0_xy_0, g_0_xxxxz_0_xy_1, g_0_xxxxzz_0_xx_0, g_0_xxxxzz_0_xy_0, g_0_xxxxzz_0_xz_0, g_0_xxxxzz_0_yy_0, g_0_xxxxzz_0_yz_0, g_0_xxxxzz_0_zz_0, g_0_xxxzz_0_xz_0, g_0_xxxzz_0_xz_1, g_0_xxxzz_0_yy_0, g_0_xxxzz_0_yy_1, g_0_xxxzz_0_yz_0, g_0_xxxzz_0_yz_1, g_0_xxxzz_0_z_1, g_0_xxxzz_0_zz_0, g_0_xxxzz_0_zz_1, g_0_xxzz_0_xz_0, g_0_xxzz_0_xz_1, g_0_xxzz_0_yy_0, g_0_xxzz_0_yy_1, g_0_xxzz_0_yz_0, g_0_xxzz_0_yz_1, g_0_xxzz_0_zz_0, g_0_xxzz_0_zz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxzz_0_xx_0[i] = g_0_xxxx_0_xx_0[i] * fi_ab_0 - g_0_xxxx_0_xx_1[i] * fti_ab_0 + g_0_xxxxz_0_xx_0[i] * pb_z + g_0_xxxxz_0_xx_1[i] * wp_z[i];

        g_0_xxxxzz_0_xy_0[i] = g_0_xxxx_0_xy_0[i] * fi_ab_0 - g_0_xxxx_0_xy_1[i] * fti_ab_0 + g_0_xxxxz_0_xy_0[i] * pb_z + g_0_xxxxz_0_xy_1[i] * wp_z[i];

        g_0_xxxxzz_0_xz_0[i] = 3.0 * g_0_xxzz_0_xz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xz_1[i] * fti_ab_0 + g_0_xxxzz_0_z_1[i] * fi_abcd_0 + g_0_xxxzz_0_xz_0[i] * pb_x + g_0_xxxzz_0_xz_1[i] * wp_x[i];

        g_0_xxxxzz_0_yy_0[i] = 3.0 * g_0_xxzz_0_yy_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_yy_1[i] * fti_ab_0 + g_0_xxxzz_0_yy_0[i] * pb_x + g_0_xxxzz_0_yy_1[i] * wp_x[i];

        g_0_xxxxzz_0_yz_0[i] = 3.0 * g_0_xxzz_0_yz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_yz_1[i] * fti_ab_0 + g_0_xxxzz_0_yz_0[i] * pb_x + g_0_xxxzz_0_yz_1[i] * wp_x[i];

        g_0_xxxxzz_0_zz_0[i] = 3.0 * g_0_xxzz_0_zz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_zz_1[i] * fti_ab_0 + g_0_xxxzz_0_zz_0[i] * pb_x + g_0_xxxzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 36-42 components of targeted buffer : SISD

    auto g_0_xxxyyy_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 36);

    auto g_0_xxxyyy_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 37);

    auto g_0_xxxyyy_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 38);

    auto g_0_xxxyyy_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 39);

    auto g_0_xxxyyy_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 40);

    auto g_0_xxxyyy_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 41);

    #pragma omp simd aligned(g_0_xxxy_0_xx_0, g_0_xxxy_0_xx_1, g_0_xxxy_0_xz_0, g_0_xxxy_0_xz_1, g_0_xxxyy_0_xx_0, g_0_xxxyy_0_xx_1, g_0_xxxyy_0_xz_0, g_0_xxxyy_0_xz_1, g_0_xxxyyy_0_xx_0, g_0_xxxyyy_0_xy_0, g_0_xxxyyy_0_xz_0, g_0_xxxyyy_0_yy_0, g_0_xxxyyy_0_yz_0, g_0_xxxyyy_0_zz_0, g_0_xxyyy_0_xy_0, g_0_xxyyy_0_xy_1, g_0_xxyyy_0_y_1, g_0_xxyyy_0_yy_0, g_0_xxyyy_0_yy_1, g_0_xxyyy_0_yz_0, g_0_xxyyy_0_yz_1, g_0_xxyyy_0_zz_0, g_0_xxyyy_0_zz_1, g_0_xyyy_0_xy_0, g_0_xyyy_0_xy_1, g_0_xyyy_0_yy_0, g_0_xyyy_0_yy_1, g_0_xyyy_0_yz_0, g_0_xyyy_0_yz_1, g_0_xyyy_0_zz_0, g_0_xyyy_0_zz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxyyy_0_xx_0[i] = 2.0 * g_0_xxxy_0_xx_0[i] * fi_ab_0 - 2.0 * g_0_xxxy_0_xx_1[i] * fti_ab_0 + g_0_xxxyy_0_xx_0[i] * pb_y + g_0_xxxyy_0_xx_1[i] * wp_y[i];

        g_0_xxxyyy_0_xy_0[i] = 2.0 * g_0_xyyy_0_xy_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xy_1[i] * fti_ab_0 + g_0_xxyyy_0_y_1[i] * fi_abcd_0 + g_0_xxyyy_0_xy_0[i] * pb_x + g_0_xxyyy_0_xy_1[i] * wp_x[i];

        g_0_xxxyyy_0_xz_0[i] = 2.0 * g_0_xxxy_0_xz_0[i] * fi_ab_0 - 2.0 * g_0_xxxy_0_xz_1[i] * fti_ab_0 + g_0_xxxyy_0_xz_0[i] * pb_y + g_0_xxxyy_0_xz_1[i] * wp_y[i];

        g_0_xxxyyy_0_yy_0[i] = 2.0 * g_0_xyyy_0_yy_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_yy_1[i] * fti_ab_0 + g_0_xxyyy_0_yy_0[i] * pb_x + g_0_xxyyy_0_yy_1[i] * wp_x[i];

        g_0_xxxyyy_0_yz_0[i] = 2.0 * g_0_xyyy_0_yz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_yz_1[i] * fti_ab_0 + g_0_xxyyy_0_yz_0[i] * pb_x + g_0_xxyyy_0_yz_1[i] * wp_x[i];

        g_0_xxxyyy_0_zz_0[i] = 2.0 * g_0_xyyy_0_zz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_zz_1[i] * fti_ab_0 + g_0_xxyyy_0_zz_0[i] * pb_x + g_0_xxyyy_0_zz_1[i] * wp_x[i];
    }

    /// Set up 42-48 components of targeted buffer : SISD

    auto g_0_xxxyyz_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 42);

    auto g_0_xxxyyz_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 43);

    auto g_0_xxxyyz_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 44);

    auto g_0_xxxyyz_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 45);

    auto g_0_xxxyyz_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 46);

    auto g_0_xxxyyz_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 47);

    #pragma omp simd aligned(g_0_xxxyy_0_x_1, g_0_xxxyy_0_xx_0, g_0_xxxyy_0_xx_1, g_0_xxxyy_0_xy_0, g_0_xxxyy_0_xy_1, g_0_xxxyy_0_xz_0, g_0_xxxyy_0_xz_1, g_0_xxxyy_0_y_1, g_0_xxxyy_0_yy_0, g_0_xxxyy_0_yy_1, g_0_xxxyy_0_yz_0, g_0_xxxyy_0_yz_1, g_0_xxxyy_0_z_1, g_0_xxxyy_0_zz_0, g_0_xxxyy_0_zz_1, g_0_xxxyyz_0_xx_0, g_0_xxxyyz_0_xy_0, g_0_xxxyyz_0_xz_0, g_0_xxxyyz_0_yy_0, g_0_xxxyyz_0_yz_0, g_0_xxxyyz_0_zz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyyz_0_xx_0[i] = g_0_xxxyy_0_xx_0[i] * pb_z + g_0_xxxyy_0_xx_1[i] * wp_z[i];

        g_0_xxxyyz_0_xy_0[i] = g_0_xxxyy_0_xy_0[i] * pb_z + g_0_xxxyy_0_xy_1[i] * wp_z[i];

        g_0_xxxyyz_0_xz_0[i] = g_0_xxxyy_0_x_1[i] * fi_abcd_0 + g_0_xxxyy_0_xz_0[i] * pb_z + g_0_xxxyy_0_xz_1[i] * wp_z[i];

        g_0_xxxyyz_0_yy_0[i] = g_0_xxxyy_0_yy_0[i] * pb_z + g_0_xxxyy_0_yy_1[i] * wp_z[i];

        g_0_xxxyyz_0_yz_0[i] = g_0_xxxyy_0_y_1[i] * fi_abcd_0 + g_0_xxxyy_0_yz_0[i] * pb_z + g_0_xxxyy_0_yz_1[i] * wp_z[i];

        g_0_xxxyyz_0_zz_0[i] = 2.0 * g_0_xxxyy_0_z_1[i] * fi_abcd_0 + g_0_xxxyy_0_zz_0[i] * pb_z + g_0_xxxyy_0_zz_1[i] * wp_z[i];
    }

    /// Set up 48-54 components of targeted buffer : SISD

    auto g_0_xxxyzz_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 48);

    auto g_0_xxxyzz_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 49);

    auto g_0_xxxyzz_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 50);

    auto g_0_xxxyzz_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 51);

    auto g_0_xxxyzz_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 52);

    auto g_0_xxxyzz_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 53);

    #pragma omp simd aligned(g_0_xxxyzz_0_xx_0, g_0_xxxyzz_0_xy_0, g_0_xxxyzz_0_xz_0, g_0_xxxyzz_0_yy_0, g_0_xxxyzz_0_yz_0, g_0_xxxyzz_0_zz_0, g_0_xxxzz_0_x_1, g_0_xxxzz_0_xx_0, g_0_xxxzz_0_xx_1, g_0_xxxzz_0_xy_0, g_0_xxxzz_0_xy_1, g_0_xxxzz_0_xz_0, g_0_xxxzz_0_xz_1, g_0_xxxzz_0_y_1, g_0_xxxzz_0_yy_0, g_0_xxxzz_0_yy_1, g_0_xxxzz_0_yz_0, g_0_xxxzz_0_yz_1, g_0_xxxzz_0_z_1, g_0_xxxzz_0_zz_0, g_0_xxxzz_0_zz_1, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyzz_0_xx_0[i] = g_0_xxxzz_0_xx_0[i] * pb_y + g_0_xxxzz_0_xx_1[i] * wp_y[i];

        g_0_xxxyzz_0_xy_0[i] = g_0_xxxzz_0_x_1[i] * fi_abcd_0 + g_0_xxxzz_0_xy_0[i] * pb_y + g_0_xxxzz_0_xy_1[i] * wp_y[i];

        g_0_xxxyzz_0_xz_0[i] = g_0_xxxzz_0_xz_0[i] * pb_y + g_0_xxxzz_0_xz_1[i] * wp_y[i];

        g_0_xxxyzz_0_yy_0[i] = 2.0 * g_0_xxxzz_0_y_1[i] * fi_abcd_0 + g_0_xxxzz_0_yy_0[i] * pb_y + g_0_xxxzz_0_yy_1[i] * wp_y[i];

        g_0_xxxyzz_0_yz_0[i] = g_0_xxxzz_0_z_1[i] * fi_abcd_0 + g_0_xxxzz_0_yz_0[i] * pb_y + g_0_xxxzz_0_yz_1[i] * wp_y[i];

        g_0_xxxyzz_0_zz_0[i] = g_0_xxxzz_0_zz_0[i] * pb_y + g_0_xxxzz_0_zz_1[i] * wp_y[i];
    }

    /// Set up 54-60 components of targeted buffer : SISD

    auto g_0_xxxzzz_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 54);

    auto g_0_xxxzzz_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 55);

    auto g_0_xxxzzz_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 56);

    auto g_0_xxxzzz_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 57);

    auto g_0_xxxzzz_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 58);

    auto g_0_xxxzzz_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 59);

    #pragma omp simd aligned(g_0_xxxz_0_xx_0, g_0_xxxz_0_xx_1, g_0_xxxz_0_xy_0, g_0_xxxz_0_xy_1, g_0_xxxzz_0_xx_0, g_0_xxxzz_0_xx_1, g_0_xxxzz_0_xy_0, g_0_xxxzz_0_xy_1, g_0_xxxzzz_0_xx_0, g_0_xxxzzz_0_xy_0, g_0_xxxzzz_0_xz_0, g_0_xxxzzz_0_yy_0, g_0_xxxzzz_0_yz_0, g_0_xxxzzz_0_zz_0, g_0_xxzzz_0_xz_0, g_0_xxzzz_0_xz_1, g_0_xxzzz_0_yy_0, g_0_xxzzz_0_yy_1, g_0_xxzzz_0_yz_0, g_0_xxzzz_0_yz_1, g_0_xxzzz_0_z_1, g_0_xxzzz_0_zz_0, g_0_xxzzz_0_zz_1, g_0_xzzz_0_xz_0, g_0_xzzz_0_xz_1, g_0_xzzz_0_yy_0, g_0_xzzz_0_yy_1, g_0_xzzz_0_yz_0, g_0_xzzz_0_yz_1, g_0_xzzz_0_zz_0, g_0_xzzz_0_zz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxzzz_0_xx_0[i] = 2.0 * g_0_xxxz_0_xx_0[i] * fi_ab_0 - 2.0 * g_0_xxxz_0_xx_1[i] * fti_ab_0 + g_0_xxxzz_0_xx_0[i] * pb_z + g_0_xxxzz_0_xx_1[i] * wp_z[i];

        g_0_xxxzzz_0_xy_0[i] = 2.0 * g_0_xxxz_0_xy_0[i] * fi_ab_0 - 2.0 * g_0_xxxz_0_xy_1[i] * fti_ab_0 + g_0_xxxzz_0_xy_0[i] * pb_z + g_0_xxxzz_0_xy_1[i] * wp_z[i];

        g_0_xxxzzz_0_xz_0[i] = 2.0 * g_0_xzzz_0_xz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xz_1[i] * fti_ab_0 + g_0_xxzzz_0_z_1[i] * fi_abcd_0 + g_0_xxzzz_0_xz_0[i] * pb_x + g_0_xxzzz_0_xz_1[i] * wp_x[i];

        g_0_xxxzzz_0_yy_0[i] = 2.0 * g_0_xzzz_0_yy_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_yy_1[i] * fti_ab_0 + g_0_xxzzz_0_yy_0[i] * pb_x + g_0_xxzzz_0_yy_1[i] * wp_x[i];

        g_0_xxxzzz_0_yz_0[i] = 2.0 * g_0_xzzz_0_yz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_yz_1[i] * fti_ab_0 + g_0_xxzzz_0_yz_0[i] * pb_x + g_0_xxzzz_0_yz_1[i] * wp_x[i];

        g_0_xxxzzz_0_zz_0[i] = 2.0 * g_0_xzzz_0_zz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_zz_1[i] * fti_ab_0 + g_0_xxzzz_0_zz_0[i] * pb_x + g_0_xxzzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 60-66 components of targeted buffer : SISD

    auto g_0_xxyyyy_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 60);

    auto g_0_xxyyyy_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 61);

    auto g_0_xxyyyy_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 62);

    auto g_0_xxyyyy_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 63);

    auto g_0_xxyyyy_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 64);

    auto g_0_xxyyyy_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 65);

    #pragma omp simd aligned(g_0_xxyy_0_xx_0, g_0_xxyy_0_xx_1, g_0_xxyy_0_xz_0, g_0_xxyy_0_xz_1, g_0_xxyyy_0_xx_0, g_0_xxyyy_0_xx_1, g_0_xxyyy_0_xz_0, g_0_xxyyy_0_xz_1, g_0_xxyyyy_0_xx_0, g_0_xxyyyy_0_xy_0, g_0_xxyyyy_0_xz_0, g_0_xxyyyy_0_yy_0, g_0_xxyyyy_0_yz_0, g_0_xxyyyy_0_zz_0, g_0_xyyyy_0_xy_0, g_0_xyyyy_0_xy_1, g_0_xyyyy_0_y_1, g_0_xyyyy_0_yy_0, g_0_xyyyy_0_yy_1, g_0_xyyyy_0_yz_0, g_0_xyyyy_0_yz_1, g_0_xyyyy_0_zz_0, g_0_xyyyy_0_zz_1, g_0_yyyy_0_xy_0, g_0_yyyy_0_xy_1, g_0_yyyy_0_yy_0, g_0_yyyy_0_yy_1, g_0_yyyy_0_yz_0, g_0_yyyy_0_yz_1, g_0_yyyy_0_zz_0, g_0_yyyy_0_zz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyyy_0_xx_0[i] = 3.0 * g_0_xxyy_0_xx_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xx_1[i] * fti_ab_0 + g_0_xxyyy_0_xx_0[i] * pb_y + g_0_xxyyy_0_xx_1[i] * wp_y[i];

        g_0_xxyyyy_0_xy_0[i] = g_0_yyyy_0_xy_0[i] * fi_ab_0 - g_0_yyyy_0_xy_1[i] * fti_ab_0 + g_0_xyyyy_0_y_1[i] * fi_abcd_0 + g_0_xyyyy_0_xy_0[i] * pb_x + g_0_xyyyy_0_xy_1[i] * wp_x[i];

        g_0_xxyyyy_0_xz_0[i] = 3.0 * g_0_xxyy_0_xz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xz_1[i] * fti_ab_0 + g_0_xxyyy_0_xz_0[i] * pb_y + g_0_xxyyy_0_xz_1[i] * wp_y[i];

        g_0_xxyyyy_0_yy_0[i] = g_0_yyyy_0_yy_0[i] * fi_ab_0 - g_0_yyyy_0_yy_1[i] * fti_ab_0 + g_0_xyyyy_0_yy_0[i] * pb_x + g_0_xyyyy_0_yy_1[i] * wp_x[i];

        g_0_xxyyyy_0_yz_0[i] = g_0_yyyy_0_yz_0[i] * fi_ab_0 - g_0_yyyy_0_yz_1[i] * fti_ab_0 + g_0_xyyyy_0_yz_0[i] * pb_x + g_0_xyyyy_0_yz_1[i] * wp_x[i];

        g_0_xxyyyy_0_zz_0[i] = g_0_yyyy_0_zz_0[i] * fi_ab_0 - g_0_yyyy_0_zz_1[i] * fti_ab_0 + g_0_xyyyy_0_zz_0[i] * pb_x + g_0_xyyyy_0_zz_1[i] * wp_x[i];
    }

    /// Set up 66-72 components of targeted buffer : SISD

    auto g_0_xxyyyz_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 66);

    auto g_0_xxyyyz_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 67);

    auto g_0_xxyyyz_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 68);

    auto g_0_xxyyyz_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 69);

    auto g_0_xxyyyz_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 70);

    auto g_0_xxyyyz_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 71);

    #pragma omp simd aligned(g_0_xxyyy_0_x_1, g_0_xxyyy_0_xx_0, g_0_xxyyy_0_xx_1, g_0_xxyyy_0_xy_0, g_0_xxyyy_0_xy_1, g_0_xxyyy_0_xz_0, g_0_xxyyy_0_xz_1, g_0_xxyyy_0_y_1, g_0_xxyyy_0_yy_0, g_0_xxyyy_0_yy_1, g_0_xxyyy_0_yz_0, g_0_xxyyy_0_yz_1, g_0_xxyyy_0_z_1, g_0_xxyyy_0_zz_0, g_0_xxyyy_0_zz_1, g_0_xxyyyz_0_xx_0, g_0_xxyyyz_0_xy_0, g_0_xxyyyz_0_xz_0, g_0_xxyyyz_0_yy_0, g_0_xxyyyz_0_yz_0, g_0_xxyyyz_0_zz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyyz_0_xx_0[i] = g_0_xxyyy_0_xx_0[i] * pb_z + g_0_xxyyy_0_xx_1[i] * wp_z[i];

        g_0_xxyyyz_0_xy_0[i] = g_0_xxyyy_0_xy_0[i] * pb_z + g_0_xxyyy_0_xy_1[i] * wp_z[i];

        g_0_xxyyyz_0_xz_0[i] = g_0_xxyyy_0_x_1[i] * fi_abcd_0 + g_0_xxyyy_0_xz_0[i] * pb_z + g_0_xxyyy_0_xz_1[i] * wp_z[i];

        g_0_xxyyyz_0_yy_0[i] = g_0_xxyyy_0_yy_0[i] * pb_z + g_0_xxyyy_0_yy_1[i] * wp_z[i];

        g_0_xxyyyz_0_yz_0[i] = g_0_xxyyy_0_y_1[i] * fi_abcd_0 + g_0_xxyyy_0_yz_0[i] * pb_z + g_0_xxyyy_0_yz_1[i] * wp_z[i];

        g_0_xxyyyz_0_zz_0[i] = 2.0 * g_0_xxyyy_0_z_1[i] * fi_abcd_0 + g_0_xxyyy_0_zz_0[i] * pb_z + g_0_xxyyy_0_zz_1[i] * wp_z[i];
    }

    /// Set up 72-78 components of targeted buffer : SISD

    auto g_0_xxyyzz_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 72);

    auto g_0_xxyyzz_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 73);

    auto g_0_xxyyzz_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 74);

    auto g_0_xxyyzz_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 75);

    auto g_0_xxyyzz_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 76);

    auto g_0_xxyyzz_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 77);

    #pragma omp simd aligned(g_0_xxyy_0_xy_0, g_0_xxyy_0_xy_1, g_0_xxyyz_0_xy_0, g_0_xxyyz_0_xy_1, g_0_xxyyzz_0_xx_0, g_0_xxyyzz_0_xy_0, g_0_xxyyzz_0_xz_0, g_0_xxyyzz_0_yy_0, g_0_xxyyzz_0_yz_0, g_0_xxyyzz_0_zz_0, g_0_xxyzz_0_xx_0, g_0_xxyzz_0_xx_1, g_0_xxyzz_0_xz_0, g_0_xxyzz_0_xz_1, g_0_xxzz_0_xx_0, g_0_xxzz_0_xx_1, g_0_xxzz_0_xz_0, g_0_xxzz_0_xz_1, g_0_xyyzz_0_yy_0, g_0_xyyzz_0_yy_1, g_0_xyyzz_0_yz_0, g_0_xyyzz_0_yz_1, g_0_xyyzz_0_zz_0, g_0_xyyzz_0_zz_1, g_0_yyzz_0_yy_0, g_0_yyzz_0_yy_1, g_0_yyzz_0_yz_0, g_0_yyzz_0_yz_1, g_0_yyzz_0_zz_0, g_0_yyzz_0_zz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyzz_0_xx_0[i] = g_0_xxzz_0_xx_0[i] * fi_ab_0 - g_0_xxzz_0_xx_1[i] * fti_ab_0 + g_0_xxyzz_0_xx_0[i] * pb_y + g_0_xxyzz_0_xx_1[i] * wp_y[i];

        g_0_xxyyzz_0_xy_0[i] = g_0_xxyy_0_xy_0[i] * fi_ab_0 - g_0_xxyy_0_xy_1[i] * fti_ab_0 + g_0_xxyyz_0_xy_0[i] * pb_z + g_0_xxyyz_0_xy_1[i] * wp_z[i];

        g_0_xxyyzz_0_xz_0[i] = g_0_xxzz_0_xz_0[i] * fi_ab_0 - g_0_xxzz_0_xz_1[i] * fti_ab_0 + g_0_xxyzz_0_xz_0[i] * pb_y + g_0_xxyzz_0_xz_1[i] * wp_y[i];

        g_0_xxyyzz_0_yy_0[i] = g_0_yyzz_0_yy_0[i] * fi_ab_0 - g_0_yyzz_0_yy_1[i] * fti_ab_0 + g_0_xyyzz_0_yy_0[i] * pb_x + g_0_xyyzz_0_yy_1[i] * wp_x[i];

        g_0_xxyyzz_0_yz_0[i] = g_0_yyzz_0_yz_0[i] * fi_ab_0 - g_0_yyzz_0_yz_1[i] * fti_ab_0 + g_0_xyyzz_0_yz_0[i] * pb_x + g_0_xyyzz_0_yz_1[i] * wp_x[i];

        g_0_xxyyzz_0_zz_0[i] = g_0_yyzz_0_zz_0[i] * fi_ab_0 - g_0_yyzz_0_zz_1[i] * fti_ab_0 + g_0_xyyzz_0_zz_0[i] * pb_x + g_0_xyyzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 78-84 components of targeted buffer : SISD

    auto g_0_xxyzzz_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 78);

    auto g_0_xxyzzz_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 79);

    auto g_0_xxyzzz_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 80);

    auto g_0_xxyzzz_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 81);

    auto g_0_xxyzzz_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 82);

    auto g_0_xxyzzz_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 83);

    #pragma omp simd aligned(g_0_xxyzzz_0_xx_0, g_0_xxyzzz_0_xy_0, g_0_xxyzzz_0_xz_0, g_0_xxyzzz_0_yy_0, g_0_xxyzzz_0_yz_0, g_0_xxyzzz_0_zz_0, g_0_xxzzz_0_x_1, g_0_xxzzz_0_xx_0, g_0_xxzzz_0_xx_1, g_0_xxzzz_0_xy_0, g_0_xxzzz_0_xy_1, g_0_xxzzz_0_xz_0, g_0_xxzzz_0_xz_1, g_0_xxzzz_0_y_1, g_0_xxzzz_0_yy_0, g_0_xxzzz_0_yy_1, g_0_xxzzz_0_yz_0, g_0_xxzzz_0_yz_1, g_0_xxzzz_0_z_1, g_0_xxzzz_0_zz_0, g_0_xxzzz_0_zz_1, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyzzz_0_xx_0[i] = g_0_xxzzz_0_xx_0[i] * pb_y + g_0_xxzzz_0_xx_1[i] * wp_y[i];

        g_0_xxyzzz_0_xy_0[i] = g_0_xxzzz_0_x_1[i] * fi_abcd_0 + g_0_xxzzz_0_xy_0[i] * pb_y + g_0_xxzzz_0_xy_1[i] * wp_y[i];

        g_0_xxyzzz_0_xz_0[i] = g_0_xxzzz_0_xz_0[i] * pb_y + g_0_xxzzz_0_xz_1[i] * wp_y[i];

        g_0_xxyzzz_0_yy_0[i] = 2.0 * g_0_xxzzz_0_y_1[i] * fi_abcd_0 + g_0_xxzzz_0_yy_0[i] * pb_y + g_0_xxzzz_0_yy_1[i] * wp_y[i];

        g_0_xxyzzz_0_yz_0[i] = g_0_xxzzz_0_z_1[i] * fi_abcd_0 + g_0_xxzzz_0_yz_0[i] * pb_y + g_0_xxzzz_0_yz_1[i] * wp_y[i];

        g_0_xxyzzz_0_zz_0[i] = g_0_xxzzz_0_zz_0[i] * pb_y + g_0_xxzzz_0_zz_1[i] * wp_y[i];
    }

    /// Set up 84-90 components of targeted buffer : SISD

    auto g_0_xxzzzz_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 84);

    auto g_0_xxzzzz_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 85);

    auto g_0_xxzzzz_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 86);

    auto g_0_xxzzzz_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 87);

    auto g_0_xxzzzz_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 88);

    auto g_0_xxzzzz_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 89);

    #pragma omp simd aligned(g_0_xxzz_0_xx_0, g_0_xxzz_0_xx_1, g_0_xxzz_0_xy_0, g_0_xxzz_0_xy_1, g_0_xxzzz_0_xx_0, g_0_xxzzz_0_xx_1, g_0_xxzzz_0_xy_0, g_0_xxzzz_0_xy_1, g_0_xxzzzz_0_xx_0, g_0_xxzzzz_0_xy_0, g_0_xxzzzz_0_xz_0, g_0_xxzzzz_0_yy_0, g_0_xxzzzz_0_yz_0, g_0_xxzzzz_0_zz_0, g_0_xzzzz_0_xz_0, g_0_xzzzz_0_xz_1, g_0_xzzzz_0_yy_0, g_0_xzzzz_0_yy_1, g_0_xzzzz_0_yz_0, g_0_xzzzz_0_yz_1, g_0_xzzzz_0_z_1, g_0_xzzzz_0_zz_0, g_0_xzzzz_0_zz_1, g_0_zzzz_0_xz_0, g_0_zzzz_0_xz_1, g_0_zzzz_0_yy_0, g_0_zzzz_0_yy_1, g_0_zzzz_0_yz_0, g_0_zzzz_0_yz_1, g_0_zzzz_0_zz_0, g_0_zzzz_0_zz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxzzzz_0_xx_0[i] = 3.0 * g_0_xxzz_0_xx_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xx_1[i] * fti_ab_0 + g_0_xxzzz_0_xx_0[i] * pb_z + g_0_xxzzz_0_xx_1[i] * wp_z[i];

        g_0_xxzzzz_0_xy_0[i] = 3.0 * g_0_xxzz_0_xy_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xy_1[i] * fti_ab_0 + g_0_xxzzz_0_xy_0[i] * pb_z + g_0_xxzzz_0_xy_1[i] * wp_z[i];

        g_0_xxzzzz_0_xz_0[i] = g_0_zzzz_0_xz_0[i] * fi_ab_0 - g_0_zzzz_0_xz_1[i] * fti_ab_0 + g_0_xzzzz_0_z_1[i] * fi_abcd_0 + g_0_xzzzz_0_xz_0[i] * pb_x + g_0_xzzzz_0_xz_1[i] * wp_x[i];

        g_0_xxzzzz_0_yy_0[i] = g_0_zzzz_0_yy_0[i] * fi_ab_0 - g_0_zzzz_0_yy_1[i] * fti_ab_0 + g_0_xzzzz_0_yy_0[i] * pb_x + g_0_xzzzz_0_yy_1[i] * wp_x[i];

        g_0_xxzzzz_0_yz_0[i] = g_0_zzzz_0_yz_0[i] * fi_ab_0 - g_0_zzzz_0_yz_1[i] * fti_ab_0 + g_0_xzzzz_0_yz_0[i] * pb_x + g_0_xzzzz_0_yz_1[i] * wp_x[i];

        g_0_xxzzzz_0_zz_0[i] = g_0_zzzz_0_zz_0[i] * fi_ab_0 - g_0_zzzz_0_zz_1[i] * fti_ab_0 + g_0_xzzzz_0_zz_0[i] * pb_x + g_0_xzzzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 90-96 components of targeted buffer : SISD

    auto g_0_xyyyyy_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 90);

    auto g_0_xyyyyy_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 91);

    auto g_0_xyyyyy_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 92);

    auto g_0_xyyyyy_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 93);

    auto g_0_xyyyyy_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 94);

    auto g_0_xyyyyy_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 95);

    #pragma omp simd aligned(g_0_xyyyyy_0_xx_0, g_0_xyyyyy_0_xy_0, g_0_xyyyyy_0_xz_0, g_0_xyyyyy_0_yy_0, g_0_xyyyyy_0_yz_0, g_0_xyyyyy_0_zz_0, g_0_yyyyy_0_x_1, g_0_yyyyy_0_xx_0, g_0_yyyyy_0_xx_1, g_0_yyyyy_0_xy_0, g_0_yyyyy_0_xy_1, g_0_yyyyy_0_xz_0, g_0_yyyyy_0_xz_1, g_0_yyyyy_0_y_1, g_0_yyyyy_0_yy_0, g_0_yyyyy_0_yy_1, g_0_yyyyy_0_yz_0, g_0_yyyyy_0_yz_1, g_0_yyyyy_0_z_1, g_0_yyyyy_0_zz_0, g_0_yyyyy_0_zz_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyy_0_xx_0[i] = 2.0 * g_0_yyyyy_0_x_1[i] * fi_abcd_0 + g_0_yyyyy_0_xx_0[i] * pb_x + g_0_yyyyy_0_xx_1[i] * wp_x[i];

        g_0_xyyyyy_0_xy_0[i] = g_0_yyyyy_0_y_1[i] * fi_abcd_0 + g_0_yyyyy_0_xy_0[i] * pb_x + g_0_yyyyy_0_xy_1[i] * wp_x[i];

        g_0_xyyyyy_0_xz_0[i] = g_0_yyyyy_0_z_1[i] * fi_abcd_0 + g_0_yyyyy_0_xz_0[i] * pb_x + g_0_yyyyy_0_xz_1[i] * wp_x[i];

        g_0_xyyyyy_0_yy_0[i] = g_0_yyyyy_0_yy_0[i] * pb_x + g_0_yyyyy_0_yy_1[i] * wp_x[i];

        g_0_xyyyyy_0_yz_0[i] = g_0_yyyyy_0_yz_0[i] * pb_x + g_0_yyyyy_0_yz_1[i] * wp_x[i];

        g_0_xyyyyy_0_zz_0[i] = g_0_yyyyy_0_zz_0[i] * pb_x + g_0_yyyyy_0_zz_1[i] * wp_x[i];
    }

    /// Set up 96-102 components of targeted buffer : SISD

    auto g_0_xyyyyz_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 96);

    auto g_0_xyyyyz_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 97);

    auto g_0_xyyyyz_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 98);

    auto g_0_xyyyyz_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 99);

    auto g_0_xyyyyz_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 100);

    auto g_0_xyyyyz_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 101);

    #pragma omp simd aligned(g_0_xyyyy_0_xx_0, g_0_xyyyy_0_xx_1, g_0_xyyyy_0_xy_0, g_0_xyyyy_0_xy_1, g_0_xyyyyz_0_xx_0, g_0_xyyyyz_0_xy_0, g_0_xyyyyz_0_xz_0, g_0_xyyyyz_0_yy_0, g_0_xyyyyz_0_yz_0, g_0_xyyyyz_0_zz_0, g_0_yyyyz_0_xz_0, g_0_yyyyz_0_xz_1, g_0_yyyyz_0_yy_0, g_0_yyyyz_0_yy_1, g_0_yyyyz_0_yz_0, g_0_yyyyz_0_yz_1, g_0_yyyyz_0_z_1, g_0_yyyyz_0_zz_0, g_0_yyyyz_0_zz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyz_0_xx_0[i] = g_0_xyyyy_0_xx_0[i] * pb_z + g_0_xyyyy_0_xx_1[i] * wp_z[i];

        g_0_xyyyyz_0_xy_0[i] = g_0_xyyyy_0_xy_0[i] * pb_z + g_0_xyyyy_0_xy_1[i] * wp_z[i];

        g_0_xyyyyz_0_xz_0[i] = g_0_yyyyz_0_z_1[i] * fi_abcd_0 + g_0_yyyyz_0_xz_0[i] * pb_x + g_0_yyyyz_0_xz_1[i] * wp_x[i];

        g_0_xyyyyz_0_yy_0[i] = g_0_yyyyz_0_yy_0[i] * pb_x + g_0_yyyyz_0_yy_1[i] * wp_x[i];

        g_0_xyyyyz_0_yz_0[i] = g_0_yyyyz_0_yz_0[i] * pb_x + g_0_yyyyz_0_yz_1[i] * wp_x[i];

        g_0_xyyyyz_0_zz_0[i] = g_0_yyyyz_0_zz_0[i] * pb_x + g_0_yyyyz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 102-108 components of targeted buffer : SISD

    auto g_0_xyyyzz_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 102);

    auto g_0_xyyyzz_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 103);

    auto g_0_xyyyzz_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 104);

    auto g_0_xyyyzz_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 105);

    auto g_0_xyyyzz_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 106);

    auto g_0_xyyyzz_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 107);

    #pragma omp simd aligned(g_0_xyyyzz_0_xx_0, g_0_xyyyzz_0_xy_0, g_0_xyyyzz_0_xz_0, g_0_xyyyzz_0_yy_0, g_0_xyyyzz_0_yz_0, g_0_xyyyzz_0_zz_0, g_0_yyyzz_0_x_1, g_0_yyyzz_0_xx_0, g_0_yyyzz_0_xx_1, g_0_yyyzz_0_xy_0, g_0_yyyzz_0_xy_1, g_0_yyyzz_0_xz_0, g_0_yyyzz_0_xz_1, g_0_yyyzz_0_y_1, g_0_yyyzz_0_yy_0, g_0_yyyzz_0_yy_1, g_0_yyyzz_0_yz_0, g_0_yyyzz_0_yz_1, g_0_yyyzz_0_z_1, g_0_yyyzz_0_zz_0, g_0_yyyzz_0_zz_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyzz_0_xx_0[i] = 2.0 * g_0_yyyzz_0_x_1[i] * fi_abcd_0 + g_0_yyyzz_0_xx_0[i] * pb_x + g_0_yyyzz_0_xx_1[i] * wp_x[i];

        g_0_xyyyzz_0_xy_0[i] = g_0_yyyzz_0_y_1[i] * fi_abcd_0 + g_0_yyyzz_0_xy_0[i] * pb_x + g_0_yyyzz_0_xy_1[i] * wp_x[i];

        g_0_xyyyzz_0_xz_0[i] = g_0_yyyzz_0_z_1[i] * fi_abcd_0 + g_0_yyyzz_0_xz_0[i] * pb_x + g_0_yyyzz_0_xz_1[i] * wp_x[i];

        g_0_xyyyzz_0_yy_0[i] = g_0_yyyzz_0_yy_0[i] * pb_x + g_0_yyyzz_0_yy_1[i] * wp_x[i];

        g_0_xyyyzz_0_yz_0[i] = g_0_yyyzz_0_yz_0[i] * pb_x + g_0_yyyzz_0_yz_1[i] * wp_x[i];

        g_0_xyyyzz_0_zz_0[i] = g_0_yyyzz_0_zz_0[i] * pb_x + g_0_yyyzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 108-114 components of targeted buffer : SISD

    auto g_0_xyyzzz_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 108);

    auto g_0_xyyzzz_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 109);

    auto g_0_xyyzzz_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 110);

    auto g_0_xyyzzz_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 111);

    auto g_0_xyyzzz_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 112);

    auto g_0_xyyzzz_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 113);

    #pragma omp simd aligned(g_0_xyyzzz_0_xx_0, g_0_xyyzzz_0_xy_0, g_0_xyyzzz_0_xz_0, g_0_xyyzzz_0_yy_0, g_0_xyyzzz_0_yz_0, g_0_xyyzzz_0_zz_0, g_0_yyzzz_0_x_1, g_0_yyzzz_0_xx_0, g_0_yyzzz_0_xx_1, g_0_yyzzz_0_xy_0, g_0_yyzzz_0_xy_1, g_0_yyzzz_0_xz_0, g_0_yyzzz_0_xz_1, g_0_yyzzz_0_y_1, g_0_yyzzz_0_yy_0, g_0_yyzzz_0_yy_1, g_0_yyzzz_0_yz_0, g_0_yyzzz_0_yz_1, g_0_yyzzz_0_z_1, g_0_yyzzz_0_zz_0, g_0_yyzzz_0_zz_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyzzz_0_xx_0[i] = 2.0 * g_0_yyzzz_0_x_1[i] * fi_abcd_0 + g_0_yyzzz_0_xx_0[i] * pb_x + g_0_yyzzz_0_xx_1[i] * wp_x[i];

        g_0_xyyzzz_0_xy_0[i] = g_0_yyzzz_0_y_1[i] * fi_abcd_0 + g_0_yyzzz_0_xy_0[i] * pb_x + g_0_yyzzz_0_xy_1[i] * wp_x[i];

        g_0_xyyzzz_0_xz_0[i] = g_0_yyzzz_0_z_1[i] * fi_abcd_0 + g_0_yyzzz_0_xz_0[i] * pb_x + g_0_yyzzz_0_xz_1[i] * wp_x[i];

        g_0_xyyzzz_0_yy_0[i] = g_0_yyzzz_0_yy_0[i] * pb_x + g_0_yyzzz_0_yy_1[i] * wp_x[i];

        g_0_xyyzzz_0_yz_0[i] = g_0_yyzzz_0_yz_0[i] * pb_x + g_0_yyzzz_0_yz_1[i] * wp_x[i];

        g_0_xyyzzz_0_zz_0[i] = g_0_yyzzz_0_zz_0[i] * pb_x + g_0_yyzzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 114-120 components of targeted buffer : SISD

    auto g_0_xyzzzz_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 114);

    auto g_0_xyzzzz_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 115);

    auto g_0_xyzzzz_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 116);

    auto g_0_xyzzzz_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 117);

    auto g_0_xyzzzz_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 118);

    auto g_0_xyzzzz_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 119);

    #pragma omp simd aligned(g_0_xyzzzz_0_xx_0, g_0_xyzzzz_0_xy_0, g_0_xyzzzz_0_xz_0, g_0_xyzzzz_0_yy_0, g_0_xyzzzz_0_yz_0, g_0_xyzzzz_0_zz_0, g_0_xzzzz_0_xx_0, g_0_xzzzz_0_xx_1, g_0_xzzzz_0_xz_0, g_0_xzzzz_0_xz_1, g_0_yzzzz_0_xy_0, g_0_yzzzz_0_xy_1, g_0_yzzzz_0_y_1, g_0_yzzzz_0_yy_0, g_0_yzzzz_0_yy_1, g_0_yzzzz_0_yz_0, g_0_yzzzz_0_yz_1, g_0_yzzzz_0_zz_0, g_0_yzzzz_0_zz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyzzzz_0_xx_0[i] = g_0_xzzzz_0_xx_0[i] * pb_y + g_0_xzzzz_0_xx_1[i] * wp_y[i];

        g_0_xyzzzz_0_xy_0[i] = g_0_yzzzz_0_y_1[i] * fi_abcd_0 + g_0_yzzzz_0_xy_0[i] * pb_x + g_0_yzzzz_0_xy_1[i] * wp_x[i];

        g_0_xyzzzz_0_xz_0[i] = g_0_xzzzz_0_xz_0[i] * pb_y + g_0_xzzzz_0_xz_1[i] * wp_y[i];

        g_0_xyzzzz_0_yy_0[i] = g_0_yzzzz_0_yy_0[i] * pb_x + g_0_yzzzz_0_yy_1[i] * wp_x[i];

        g_0_xyzzzz_0_yz_0[i] = g_0_yzzzz_0_yz_0[i] * pb_x + g_0_yzzzz_0_yz_1[i] * wp_x[i];

        g_0_xyzzzz_0_zz_0[i] = g_0_yzzzz_0_zz_0[i] * pb_x + g_0_yzzzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 120-126 components of targeted buffer : SISD

    auto g_0_xzzzzz_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 120);

    auto g_0_xzzzzz_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 121);

    auto g_0_xzzzzz_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 122);

    auto g_0_xzzzzz_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 123);

    auto g_0_xzzzzz_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 124);

    auto g_0_xzzzzz_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 125);

    #pragma omp simd aligned(g_0_xzzzzz_0_xx_0, g_0_xzzzzz_0_xy_0, g_0_xzzzzz_0_xz_0, g_0_xzzzzz_0_yy_0, g_0_xzzzzz_0_yz_0, g_0_xzzzzz_0_zz_0, g_0_zzzzz_0_x_1, g_0_zzzzz_0_xx_0, g_0_zzzzz_0_xx_1, g_0_zzzzz_0_xy_0, g_0_zzzzz_0_xy_1, g_0_zzzzz_0_xz_0, g_0_zzzzz_0_xz_1, g_0_zzzzz_0_y_1, g_0_zzzzz_0_yy_0, g_0_zzzzz_0_yy_1, g_0_zzzzz_0_yz_0, g_0_zzzzz_0_yz_1, g_0_zzzzz_0_z_1, g_0_zzzzz_0_zz_0, g_0_zzzzz_0_zz_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzzzzz_0_xx_0[i] = 2.0 * g_0_zzzzz_0_x_1[i] * fi_abcd_0 + g_0_zzzzz_0_xx_0[i] * pb_x + g_0_zzzzz_0_xx_1[i] * wp_x[i];

        g_0_xzzzzz_0_xy_0[i] = g_0_zzzzz_0_y_1[i] * fi_abcd_0 + g_0_zzzzz_0_xy_0[i] * pb_x + g_0_zzzzz_0_xy_1[i] * wp_x[i];

        g_0_xzzzzz_0_xz_0[i] = g_0_zzzzz_0_z_1[i] * fi_abcd_0 + g_0_zzzzz_0_xz_0[i] * pb_x + g_0_zzzzz_0_xz_1[i] * wp_x[i];

        g_0_xzzzzz_0_yy_0[i] = g_0_zzzzz_0_yy_0[i] * pb_x + g_0_zzzzz_0_yy_1[i] * wp_x[i];

        g_0_xzzzzz_0_yz_0[i] = g_0_zzzzz_0_yz_0[i] * pb_x + g_0_zzzzz_0_yz_1[i] * wp_x[i];

        g_0_xzzzzz_0_zz_0[i] = g_0_zzzzz_0_zz_0[i] * pb_x + g_0_zzzzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 126-132 components of targeted buffer : SISD

    auto g_0_yyyyyy_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 126);

    auto g_0_yyyyyy_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 127);

    auto g_0_yyyyyy_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 128);

    auto g_0_yyyyyy_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 129);

    auto g_0_yyyyyy_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 130);

    auto g_0_yyyyyy_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 131);

    #pragma omp simd aligned(g_0_yyyy_0_xx_0, g_0_yyyy_0_xx_1, g_0_yyyy_0_xy_0, g_0_yyyy_0_xy_1, g_0_yyyy_0_xz_0, g_0_yyyy_0_xz_1, g_0_yyyy_0_yy_0, g_0_yyyy_0_yy_1, g_0_yyyy_0_yz_0, g_0_yyyy_0_yz_1, g_0_yyyy_0_zz_0, g_0_yyyy_0_zz_1, g_0_yyyyy_0_x_1, g_0_yyyyy_0_xx_0, g_0_yyyyy_0_xx_1, g_0_yyyyy_0_xy_0, g_0_yyyyy_0_xy_1, g_0_yyyyy_0_xz_0, g_0_yyyyy_0_xz_1, g_0_yyyyy_0_y_1, g_0_yyyyy_0_yy_0, g_0_yyyyy_0_yy_1, g_0_yyyyy_0_yz_0, g_0_yyyyy_0_yz_1, g_0_yyyyy_0_z_1, g_0_yyyyy_0_zz_0, g_0_yyyyy_0_zz_1, g_0_yyyyyy_0_xx_0, g_0_yyyyyy_0_xy_0, g_0_yyyyyy_0_xz_0, g_0_yyyyyy_0_yy_0, g_0_yyyyyy_0_yz_0, g_0_yyyyyy_0_zz_0, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyyy_0_xx_0[i] = 5.0 * g_0_yyyy_0_xx_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xx_1[i] * fti_ab_0 + g_0_yyyyy_0_xx_0[i] * pb_y + g_0_yyyyy_0_xx_1[i] * wp_y[i];

        g_0_yyyyyy_0_xy_0[i] = 5.0 * g_0_yyyy_0_xy_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xy_1[i] * fti_ab_0 + g_0_yyyyy_0_x_1[i] * fi_abcd_0 + g_0_yyyyy_0_xy_0[i] * pb_y + g_0_yyyyy_0_xy_1[i] * wp_y[i];

        g_0_yyyyyy_0_xz_0[i] = 5.0 * g_0_yyyy_0_xz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xz_1[i] * fti_ab_0 + g_0_yyyyy_0_xz_0[i] * pb_y + g_0_yyyyy_0_xz_1[i] * wp_y[i];

        g_0_yyyyyy_0_yy_0[i] = 5.0 * g_0_yyyy_0_yy_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_yy_1[i] * fti_ab_0 + 2.0 * g_0_yyyyy_0_y_1[i] * fi_abcd_0 + g_0_yyyyy_0_yy_0[i] * pb_y + g_0_yyyyy_0_yy_1[i] * wp_y[i];

        g_0_yyyyyy_0_yz_0[i] = 5.0 * g_0_yyyy_0_yz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_yz_1[i] * fti_ab_0 + g_0_yyyyy_0_z_1[i] * fi_abcd_0 + g_0_yyyyy_0_yz_0[i] * pb_y + g_0_yyyyy_0_yz_1[i] * wp_y[i];

        g_0_yyyyyy_0_zz_0[i] = 5.0 * g_0_yyyy_0_zz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_zz_1[i] * fti_ab_0 + g_0_yyyyy_0_zz_0[i] * pb_y + g_0_yyyyy_0_zz_1[i] * wp_y[i];
    }

    /// Set up 132-138 components of targeted buffer : SISD

    auto g_0_yyyyyz_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 132);

    auto g_0_yyyyyz_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 133);

    auto g_0_yyyyyz_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 134);

    auto g_0_yyyyyz_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 135);

    auto g_0_yyyyyz_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 136);

    auto g_0_yyyyyz_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 137);

    #pragma omp simd aligned(g_0_yyyyy_0_x_1, g_0_yyyyy_0_xx_0, g_0_yyyyy_0_xx_1, g_0_yyyyy_0_xy_0, g_0_yyyyy_0_xy_1, g_0_yyyyy_0_xz_0, g_0_yyyyy_0_xz_1, g_0_yyyyy_0_y_1, g_0_yyyyy_0_yy_0, g_0_yyyyy_0_yy_1, g_0_yyyyy_0_yz_0, g_0_yyyyy_0_yz_1, g_0_yyyyy_0_z_1, g_0_yyyyy_0_zz_0, g_0_yyyyy_0_zz_1, g_0_yyyyyz_0_xx_0, g_0_yyyyyz_0_xy_0, g_0_yyyyyz_0_xz_0, g_0_yyyyyz_0_yy_0, g_0_yyyyyz_0_yz_0, g_0_yyyyyz_0_zz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyyyz_0_xx_0[i] = g_0_yyyyy_0_xx_0[i] * pb_z + g_0_yyyyy_0_xx_1[i] * wp_z[i];

        g_0_yyyyyz_0_xy_0[i] = g_0_yyyyy_0_xy_0[i] * pb_z + g_0_yyyyy_0_xy_1[i] * wp_z[i];

        g_0_yyyyyz_0_xz_0[i] = g_0_yyyyy_0_x_1[i] * fi_abcd_0 + g_0_yyyyy_0_xz_0[i] * pb_z + g_0_yyyyy_0_xz_1[i] * wp_z[i];

        g_0_yyyyyz_0_yy_0[i] = g_0_yyyyy_0_yy_0[i] * pb_z + g_0_yyyyy_0_yy_1[i] * wp_z[i];

        g_0_yyyyyz_0_yz_0[i] = g_0_yyyyy_0_y_1[i] * fi_abcd_0 + g_0_yyyyy_0_yz_0[i] * pb_z + g_0_yyyyy_0_yz_1[i] * wp_z[i];

        g_0_yyyyyz_0_zz_0[i] = 2.0 * g_0_yyyyy_0_z_1[i] * fi_abcd_0 + g_0_yyyyy_0_zz_0[i] * pb_z + g_0_yyyyy_0_zz_1[i] * wp_z[i];
    }

    /// Set up 138-144 components of targeted buffer : SISD

    auto g_0_yyyyzz_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 138);

    auto g_0_yyyyzz_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 139);

    auto g_0_yyyyzz_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 140);

    auto g_0_yyyyzz_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 141);

    auto g_0_yyyyzz_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 142);

    auto g_0_yyyyzz_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 143);

    #pragma omp simd aligned(g_0_yyyy_0_xy_0, g_0_yyyy_0_xy_1, g_0_yyyy_0_yy_0, g_0_yyyy_0_yy_1, g_0_yyyyz_0_xy_0, g_0_yyyyz_0_xy_1, g_0_yyyyz_0_yy_0, g_0_yyyyz_0_yy_1, g_0_yyyyzz_0_xx_0, g_0_yyyyzz_0_xy_0, g_0_yyyyzz_0_xz_0, g_0_yyyyzz_0_yy_0, g_0_yyyyzz_0_yz_0, g_0_yyyyzz_0_zz_0, g_0_yyyzz_0_xx_0, g_0_yyyzz_0_xx_1, g_0_yyyzz_0_xz_0, g_0_yyyzz_0_xz_1, g_0_yyyzz_0_yz_0, g_0_yyyzz_0_yz_1, g_0_yyyzz_0_z_1, g_0_yyyzz_0_zz_0, g_0_yyyzz_0_zz_1, g_0_yyzz_0_xx_0, g_0_yyzz_0_xx_1, g_0_yyzz_0_xz_0, g_0_yyzz_0_xz_1, g_0_yyzz_0_yz_0, g_0_yyzz_0_yz_1, g_0_yyzz_0_zz_0, g_0_yyzz_0_zz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyzz_0_xx_0[i] = 3.0 * g_0_yyzz_0_xx_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xx_1[i] * fti_ab_0 + g_0_yyyzz_0_xx_0[i] * pb_y + g_0_yyyzz_0_xx_1[i] * wp_y[i];

        g_0_yyyyzz_0_xy_0[i] = g_0_yyyy_0_xy_0[i] * fi_ab_0 - g_0_yyyy_0_xy_1[i] * fti_ab_0 + g_0_yyyyz_0_xy_0[i] * pb_z + g_0_yyyyz_0_xy_1[i] * wp_z[i];

        g_0_yyyyzz_0_xz_0[i] = 3.0 * g_0_yyzz_0_xz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xz_1[i] * fti_ab_0 + g_0_yyyzz_0_xz_0[i] * pb_y + g_0_yyyzz_0_xz_1[i] * wp_y[i];

        g_0_yyyyzz_0_yy_0[i] = g_0_yyyy_0_yy_0[i] * fi_ab_0 - g_0_yyyy_0_yy_1[i] * fti_ab_0 + g_0_yyyyz_0_yy_0[i] * pb_z + g_0_yyyyz_0_yy_1[i] * wp_z[i];

        g_0_yyyyzz_0_yz_0[i] = 3.0 * g_0_yyzz_0_yz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_yz_1[i] * fti_ab_0 + g_0_yyyzz_0_z_1[i] * fi_abcd_0 + g_0_yyyzz_0_yz_0[i] * pb_y + g_0_yyyzz_0_yz_1[i] * wp_y[i];

        g_0_yyyyzz_0_zz_0[i] = 3.0 * g_0_yyzz_0_zz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_zz_1[i] * fti_ab_0 + g_0_yyyzz_0_zz_0[i] * pb_y + g_0_yyyzz_0_zz_1[i] * wp_y[i];
    }

    /// Set up 144-150 components of targeted buffer : SISD

    auto g_0_yyyzzz_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 144);

    auto g_0_yyyzzz_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 145);

    auto g_0_yyyzzz_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 146);

    auto g_0_yyyzzz_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 147);

    auto g_0_yyyzzz_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 148);

    auto g_0_yyyzzz_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 149);

    #pragma omp simd aligned(g_0_yyyz_0_xy_0, g_0_yyyz_0_xy_1, g_0_yyyz_0_yy_0, g_0_yyyz_0_yy_1, g_0_yyyzz_0_xy_0, g_0_yyyzz_0_xy_1, g_0_yyyzz_0_yy_0, g_0_yyyzz_0_yy_1, g_0_yyyzzz_0_xx_0, g_0_yyyzzz_0_xy_0, g_0_yyyzzz_0_xz_0, g_0_yyyzzz_0_yy_0, g_0_yyyzzz_0_yz_0, g_0_yyyzzz_0_zz_0, g_0_yyzzz_0_xx_0, g_0_yyzzz_0_xx_1, g_0_yyzzz_0_xz_0, g_0_yyzzz_0_xz_1, g_0_yyzzz_0_yz_0, g_0_yyzzz_0_yz_1, g_0_yyzzz_0_z_1, g_0_yyzzz_0_zz_0, g_0_yyzzz_0_zz_1, g_0_yzzz_0_xx_0, g_0_yzzz_0_xx_1, g_0_yzzz_0_xz_0, g_0_yzzz_0_xz_1, g_0_yzzz_0_yz_0, g_0_yzzz_0_yz_1, g_0_yzzz_0_zz_0, g_0_yzzz_0_zz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyzzz_0_xx_0[i] = 2.0 * g_0_yzzz_0_xx_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xx_1[i] * fti_ab_0 + g_0_yyzzz_0_xx_0[i] * pb_y + g_0_yyzzz_0_xx_1[i] * wp_y[i];

        g_0_yyyzzz_0_xy_0[i] = 2.0 * g_0_yyyz_0_xy_0[i] * fi_ab_0 - 2.0 * g_0_yyyz_0_xy_1[i] * fti_ab_0 + g_0_yyyzz_0_xy_0[i] * pb_z + g_0_yyyzz_0_xy_1[i] * wp_z[i];

        g_0_yyyzzz_0_xz_0[i] = 2.0 * g_0_yzzz_0_xz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xz_1[i] * fti_ab_0 + g_0_yyzzz_0_xz_0[i] * pb_y + g_0_yyzzz_0_xz_1[i] * wp_y[i];

        g_0_yyyzzz_0_yy_0[i] = 2.0 * g_0_yyyz_0_yy_0[i] * fi_ab_0 - 2.0 * g_0_yyyz_0_yy_1[i] * fti_ab_0 + g_0_yyyzz_0_yy_0[i] * pb_z + g_0_yyyzz_0_yy_1[i] * wp_z[i];

        g_0_yyyzzz_0_yz_0[i] = 2.0 * g_0_yzzz_0_yz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_yz_1[i] * fti_ab_0 + g_0_yyzzz_0_z_1[i] * fi_abcd_0 + g_0_yyzzz_0_yz_0[i] * pb_y + g_0_yyzzz_0_yz_1[i] * wp_y[i];

        g_0_yyyzzz_0_zz_0[i] = 2.0 * g_0_yzzz_0_zz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_zz_1[i] * fti_ab_0 + g_0_yyzzz_0_zz_0[i] * pb_y + g_0_yyzzz_0_zz_1[i] * wp_y[i];
    }

    /// Set up 150-156 components of targeted buffer : SISD

    auto g_0_yyzzzz_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 150);

    auto g_0_yyzzzz_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 151);

    auto g_0_yyzzzz_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 152);

    auto g_0_yyzzzz_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 153);

    auto g_0_yyzzzz_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 154);

    auto g_0_yyzzzz_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 155);

    #pragma omp simd aligned(g_0_yyzz_0_xy_0, g_0_yyzz_0_xy_1, g_0_yyzz_0_yy_0, g_0_yyzz_0_yy_1, g_0_yyzzz_0_xy_0, g_0_yyzzz_0_xy_1, g_0_yyzzz_0_yy_0, g_0_yyzzz_0_yy_1, g_0_yyzzzz_0_xx_0, g_0_yyzzzz_0_xy_0, g_0_yyzzzz_0_xz_0, g_0_yyzzzz_0_yy_0, g_0_yyzzzz_0_yz_0, g_0_yyzzzz_0_zz_0, g_0_yzzzz_0_xx_0, g_0_yzzzz_0_xx_1, g_0_yzzzz_0_xz_0, g_0_yzzzz_0_xz_1, g_0_yzzzz_0_yz_0, g_0_yzzzz_0_yz_1, g_0_yzzzz_0_z_1, g_0_yzzzz_0_zz_0, g_0_yzzzz_0_zz_1, g_0_zzzz_0_xx_0, g_0_zzzz_0_xx_1, g_0_zzzz_0_xz_0, g_0_zzzz_0_xz_1, g_0_zzzz_0_yz_0, g_0_zzzz_0_yz_1, g_0_zzzz_0_zz_0, g_0_zzzz_0_zz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyzzzz_0_xx_0[i] = g_0_zzzz_0_xx_0[i] * fi_ab_0 - g_0_zzzz_0_xx_1[i] * fti_ab_0 + g_0_yzzzz_0_xx_0[i] * pb_y + g_0_yzzzz_0_xx_1[i] * wp_y[i];

        g_0_yyzzzz_0_xy_0[i] = 3.0 * g_0_yyzz_0_xy_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xy_1[i] * fti_ab_0 + g_0_yyzzz_0_xy_0[i] * pb_z + g_0_yyzzz_0_xy_1[i] * wp_z[i];

        g_0_yyzzzz_0_xz_0[i] = g_0_zzzz_0_xz_0[i] * fi_ab_0 - g_0_zzzz_0_xz_1[i] * fti_ab_0 + g_0_yzzzz_0_xz_0[i] * pb_y + g_0_yzzzz_0_xz_1[i] * wp_y[i];

        g_0_yyzzzz_0_yy_0[i] = 3.0 * g_0_yyzz_0_yy_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_yy_1[i] * fti_ab_0 + g_0_yyzzz_0_yy_0[i] * pb_z + g_0_yyzzz_0_yy_1[i] * wp_z[i];

        g_0_yyzzzz_0_yz_0[i] = g_0_zzzz_0_yz_0[i] * fi_ab_0 - g_0_zzzz_0_yz_1[i] * fti_ab_0 + g_0_yzzzz_0_z_1[i] * fi_abcd_0 + g_0_yzzzz_0_yz_0[i] * pb_y + g_0_yzzzz_0_yz_1[i] * wp_y[i];

        g_0_yyzzzz_0_zz_0[i] = g_0_zzzz_0_zz_0[i] * fi_ab_0 - g_0_zzzz_0_zz_1[i] * fti_ab_0 + g_0_yzzzz_0_zz_0[i] * pb_y + g_0_yzzzz_0_zz_1[i] * wp_y[i];
    }

    /// Set up 156-162 components of targeted buffer : SISD

    auto g_0_yzzzzz_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 156);

    auto g_0_yzzzzz_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 157);

    auto g_0_yzzzzz_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 158);

    auto g_0_yzzzzz_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 159);

    auto g_0_yzzzzz_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 160);

    auto g_0_yzzzzz_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 161);

    #pragma omp simd aligned(g_0_yzzzzz_0_xx_0, g_0_yzzzzz_0_xy_0, g_0_yzzzzz_0_xz_0, g_0_yzzzzz_0_yy_0, g_0_yzzzzz_0_yz_0, g_0_yzzzzz_0_zz_0, g_0_zzzzz_0_x_1, g_0_zzzzz_0_xx_0, g_0_zzzzz_0_xx_1, g_0_zzzzz_0_xy_0, g_0_zzzzz_0_xy_1, g_0_zzzzz_0_xz_0, g_0_zzzzz_0_xz_1, g_0_zzzzz_0_y_1, g_0_zzzzz_0_yy_0, g_0_zzzzz_0_yy_1, g_0_zzzzz_0_yz_0, g_0_zzzzz_0_yz_1, g_0_zzzzz_0_z_1, g_0_zzzzz_0_zz_0, g_0_zzzzz_0_zz_1, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzzzzz_0_xx_0[i] = g_0_zzzzz_0_xx_0[i] * pb_y + g_0_zzzzz_0_xx_1[i] * wp_y[i];

        g_0_yzzzzz_0_xy_0[i] = g_0_zzzzz_0_x_1[i] * fi_abcd_0 + g_0_zzzzz_0_xy_0[i] * pb_y + g_0_zzzzz_0_xy_1[i] * wp_y[i];

        g_0_yzzzzz_0_xz_0[i] = g_0_zzzzz_0_xz_0[i] * pb_y + g_0_zzzzz_0_xz_1[i] * wp_y[i];

        g_0_yzzzzz_0_yy_0[i] = 2.0 * g_0_zzzzz_0_y_1[i] * fi_abcd_0 + g_0_zzzzz_0_yy_0[i] * pb_y + g_0_zzzzz_0_yy_1[i] * wp_y[i];

        g_0_yzzzzz_0_yz_0[i] = g_0_zzzzz_0_z_1[i] * fi_abcd_0 + g_0_zzzzz_0_yz_0[i] * pb_y + g_0_zzzzz_0_yz_1[i] * wp_y[i];

        g_0_yzzzzz_0_zz_0[i] = g_0_zzzzz_0_zz_0[i] * pb_y + g_0_zzzzz_0_zz_1[i] * wp_y[i];
    }

    /// Set up 162-168 components of targeted buffer : SISD

    auto g_0_zzzzzz_0_xx_0 = pbuffer.data(idx_eri_0_sisd + 162);

    auto g_0_zzzzzz_0_xy_0 = pbuffer.data(idx_eri_0_sisd + 163);

    auto g_0_zzzzzz_0_xz_0 = pbuffer.data(idx_eri_0_sisd + 164);

    auto g_0_zzzzzz_0_yy_0 = pbuffer.data(idx_eri_0_sisd + 165);

    auto g_0_zzzzzz_0_yz_0 = pbuffer.data(idx_eri_0_sisd + 166);

    auto g_0_zzzzzz_0_zz_0 = pbuffer.data(idx_eri_0_sisd + 167);

    #pragma omp simd aligned(g_0_zzzz_0_xx_0, g_0_zzzz_0_xx_1, g_0_zzzz_0_xy_0, g_0_zzzz_0_xy_1, g_0_zzzz_0_xz_0, g_0_zzzz_0_xz_1, g_0_zzzz_0_yy_0, g_0_zzzz_0_yy_1, g_0_zzzz_0_yz_0, g_0_zzzz_0_yz_1, g_0_zzzz_0_zz_0, g_0_zzzz_0_zz_1, g_0_zzzzz_0_x_1, g_0_zzzzz_0_xx_0, g_0_zzzzz_0_xx_1, g_0_zzzzz_0_xy_0, g_0_zzzzz_0_xy_1, g_0_zzzzz_0_xz_0, g_0_zzzzz_0_xz_1, g_0_zzzzz_0_y_1, g_0_zzzzz_0_yy_0, g_0_zzzzz_0_yy_1, g_0_zzzzz_0_yz_0, g_0_zzzzz_0_yz_1, g_0_zzzzz_0_z_1, g_0_zzzzz_0_zz_0, g_0_zzzzz_0_zz_1, g_0_zzzzzz_0_xx_0, g_0_zzzzzz_0_xy_0, g_0_zzzzzz_0_xz_0, g_0_zzzzzz_0_yy_0, g_0_zzzzzz_0_yz_0, g_0_zzzzzz_0_zz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzzzzz_0_xx_0[i] = 5.0 * g_0_zzzz_0_xx_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xx_1[i] * fti_ab_0 + g_0_zzzzz_0_xx_0[i] * pb_z + g_0_zzzzz_0_xx_1[i] * wp_z[i];

        g_0_zzzzzz_0_xy_0[i] = 5.0 * g_0_zzzz_0_xy_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xy_1[i] * fti_ab_0 + g_0_zzzzz_0_xy_0[i] * pb_z + g_0_zzzzz_0_xy_1[i] * wp_z[i];

        g_0_zzzzzz_0_xz_0[i] = 5.0 * g_0_zzzz_0_xz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xz_1[i] * fti_ab_0 + g_0_zzzzz_0_x_1[i] * fi_abcd_0 + g_0_zzzzz_0_xz_0[i] * pb_z + g_0_zzzzz_0_xz_1[i] * wp_z[i];

        g_0_zzzzzz_0_yy_0[i] = 5.0 * g_0_zzzz_0_yy_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_yy_1[i] * fti_ab_0 + g_0_zzzzz_0_yy_0[i] * pb_z + g_0_zzzzz_0_yy_1[i] * wp_z[i];

        g_0_zzzzzz_0_yz_0[i] = 5.0 * g_0_zzzz_0_yz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_yz_1[i] * fti_ab_0 + g_0_zzzzz_0_y_1[i] * fi_abcd_0 + g_0_zzzzz_0_yz_0[i] * pb_z + g_0_zzzzz_0_yz_1[i] * wp_z[i];

        g_0_zzzzzz_0_zz_0[i] = 5.0 * g_0_zzzz_0_zz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_zz_1[i] * fti_ab_0 + 2.0 * g_0_zzzzz_0_z_1[i] * fi_abcd_0 + g_0_zzzzz_0_zz_0[i] * pb_z + g_0_zzzzz_0_zz_1[i] * wp_z[i];
    }
}

} // erirec namespace


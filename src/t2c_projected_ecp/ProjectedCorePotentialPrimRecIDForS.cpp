#include "ProjectedCorePotentialPrimRecIDForS.hpp"

namespace t2pecp { // t2pecp namespace

auto
comp_prim_projected_core_potential_id_s(CSimdArray<double>& pbuffer, 
                                        const size_t idx_id_s_0_0_0,
                                        const size_t idx_gd_s_0_0_0,
                                        const size_t idx_hd_s_0_0_0,
                                        const size_t idx_gd_s_1_0_0,
                                        const size_t idx_hd_s_1_0_0,
                                        const int p,
                                        const size_t idx_gd_s_0_0_1,
                                        const size_t idx_hd_s_0_0_1,
                                        const CSimdArray<double>& factors,
                                        const TPoint<double>& r_a,
                                        const double a_exp,
                                        const double c_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents on ket side

    auto b_exps = factors.data(0);

    // Set up A center coordinates

    const auto xyz = r_a.coordinates();

    const auto a_x = xyz[0];

    const auto a_y = xyz[1];

    const auto a_z = xyz[2];

    // Set up components of auxiliary buffer : GD

    auto tg_xxxx_xx_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0);

    auto tg_xxxx_xy_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 1);

    auto tg_xxxx_xz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 2);

    auto tg_xxxx_yy_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 3);

    auto tg_xxxx_yz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 4);

    auto tg_xxxx_zz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 5);

    auto tg_xxxy_xx_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 6);

    auto tg_xxxy_xy_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 7);

    auto tg_xxxy_xz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 8);

    auto tg_xxxy_yy_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 9);

    auto tg_xxxy_yz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 10);

    auto tg_xxxy_zz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 11);

    auto tg_xxxz_xx_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 12);

    auto tg_xxxz_xy_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 13);

    auto tg_xxxz_xz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 14);

    auto tg_xxxz_yy_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 15);

    auto tg_xxxz_yz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 16);

    auto tg_xxxz_zz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 17);

    auto tg_xxyy_xx_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 18);

    auto tg_xxyy_xy_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 19);

    auto tg_xxyy_xz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 20);

    auto tg_xxyy_yy_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 21);

    auto tg_xxyy_yz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 22);

    auto tg_xxyy_zz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 23);

    auto tg_xxyz_xx_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 24);

    auto tg_xxyz_xy_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 25);

    auto tg_xxyz_xz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 26);

    auto tg_xxyz_yy_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 27);

    auto tg_xxyz_yz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 28);

    auto tg_xxyz_zz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 29);

    auto tg_xxzz_xx_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 30);

    auto tg_xxzz_xy_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 31);

    auto tg_xxzz_xz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 32);

    auto tg_xxzz_yy_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 33);

    auto tg_xxzz_yz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 34);

    auto tg_xxzz_zz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 35);

    auto tg_xyyy_xx_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 36);

    auto tg_xyyy_xy_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 37);

    auto tg_xyyy_xz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 38);

    auto tg_xyyy_yy_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 39);

    auto tg_xyyy_yz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 40);

    auto tg_xyyy_zz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 41);

    auto tg_xyyz_xx_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 42);

    auto tg_xyyz_xy_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 43);

    auto tg_xyyz_xz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 44);

    auto tg_xyyz_yy_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 45);

    auto tg_xyyz_yz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 46);

    auto tg_xyyz_zz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 47);

    auto tg_xyzz_xx_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 48);

    auto tg_xyzz_xy_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 49);

    auto tg_xyzz_xz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 50);

    auto tg_xyzz_yy_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 51);

    auto tg_xyzz_yz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 52);

    auto tg_xyzz_zz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 53);

    auto tg_xzzz_xx_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 54);

    auto tg_xzzz_xy_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 55);

    auto tg_xzzz_xz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 56);

    auto tg_xzzz_yy_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 57);

    auto tg_xzzz_yz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 58);

    auto tg_xzzz_zz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 59);

    auto tg_yyyy_xx_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 60);

    auto tg_yyyy_xy_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 61);

    auto tg_yyyy_xz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 62);

    auto tg_yyyy_yy_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 63);

    auto tg_yyyy_yz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 64);

    auto tg_yyyy_zz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 65);

    auto tg_yyyz_xx_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 66);

    auto tg_yyyz_xy_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 67);

    auto tg_yyyz_xz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 68);

    auto tg_yyyz_yy_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 69);

    auto tg_yyyz_yz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 70);

    auto tg_yyyz_zz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 71);

    auto tg_yyzz_xx_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 72);

    auto tg_yyzz_xy_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 73);

    auto tg_yyzz_xz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 74);

    auto tg_yyzz_yy_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 75);

    auto tg_yyzz_yz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 76);

    auto tg_yyzz_zz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 77);

    auto tg_yzzz_xx_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 78);

    auto tg_yzzz_xy_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 79);

    auto tg_yzzz_xz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 80);

    auto tg_yzzz_yy_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 81);

    auto tg_yzzz_yz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 82);

    auto tg_yzzz_zz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 83);

    auto tg_zzzz_xx_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 84);

    auto tg_zzzz_xy_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 85);

    auto tg_zzzz_xz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 86);

    auto tg_zzzz_yy_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 87);

    auto tg_zzzz_yz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 88);

    auto tg_zzzz_zz_s_0_0_0 = pbuffer.data(idx_gd_s_0_0_0 + 89);

    // Set up components of auxiliary buffer : HD

    auto tg_xxxxx_xx_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0);

    auto tg_xxxxx_xy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 1);

    auto tg_xxxxx_xz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 2);

    auto tg_xxxxx_yy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 3);

    auto tg_xxxxx_yz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 4);

    auto tg_xxxxx_zz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 5);

    auto tg_xxxxy_xx_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 6);

    auto tg_xxxxy_xy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 7);

    auto tg_xxxxy_xz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 8);

    auto tg_xxxxy_yy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 9);

    auto tg_xxxxy_yz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 10);

    auto tg_xxxxy_zz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 11);

    auto tg_xxxxz_xx_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 12);

    auto tg_xxxxz_xy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 13);

    auto tg_xxxxz_xz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 14);

    auto tg_xxxxz_yy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 15);

    auto tg_xxxxz_yz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 16);

    auto tg_xxxxz_zz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 17);

    auto tg_xxxyy_xx_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 18);

    auto tg_xxxyy_xy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 19);

    auto tg_xxxyy_xz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 20);

    auto tg_xxxyy_yy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 21);

    auto tg_xxxyy_yz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 22);

    auto tg_xxxyy_zz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 23);

    auto tg_xxxyz_xx_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 24);

    auto tg_xxxyz_xy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 25);

    auto tg_xxxyz_xz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 26);

    auto tg_xxxyz_yy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 27);

    auto tg_xxxyz_yz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 28);

    auto tg_xxxyz_zz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 29);

    auto tg_xxxzz_xx_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 30);

    auto tg_xxxzz_xy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 31);

    auto tg_xxxzz_xz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 32);

    auto tg_xxxzz_yy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 33);

    auto tg_xxxzz_yz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 34);

    auto tg_xxxzz_zz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 35);

    auto tg_xxyyy_xx_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 36);

    auto tg_xxyyy_xy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 37);

    auto tg_xxyyy_xz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 38);

    auto tg_xxyyy_yy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 39);

    auto tg_xxyyy_yz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 40);

    auto tg_xxyyy_zz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 41);

    auto tg_xxyyz_xx_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 42);

    auto tg_xxyyz_xy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 43);

    auto tg_xxyyz_xz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 44);

    auto tg_xxyyz_yy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 45);

    auto tg_xxyyz_yz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 46);

    auto tg_xxyyz_zz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 47);

    auto tg_xxyzz_xx_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 48);

    auto tg_xxyzz_xy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 49);

    auto tg_xxyzz_xz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 50);

    auto tg_xxyzz_yy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 51);

    auto tg_xxyzz_yz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 52);

    auto tg_xxyzz_zz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 53);

    auto tg_xxzzz_xx_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 54);

    auto tg_xxzzz_xy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 55);

    auto tg_xxzzz_xz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 56);

    auto tg_xxzzz_yy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 57);

    auto tg_xxzzz_yz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 58);

    auto tg_xxzzz_zz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 59);

    auto tg_xyyyy_xx_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 60);

    auto tg_xyyyy_xy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 61);

    auto tg_xyyyy_xz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 62);

    auto tg_xyyyy_yy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 63);

    auto tg_xyyyy_yz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 64);

    auto tg_xyyyy_zz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 65);

    auto tg_xyyyz_xx_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 66);

    auto tg_xyyyz_xy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 67);

    auto tg_xyyyz_xz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 68);

    auto tg_xyyyz_yy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 69);

    auto tg_xyyyz_yz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 70);

    auto tg_xyyyz_zz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 71);

    auto tg_xyyzz_xx_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 72);

    auto tg_xyyzz_xy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 73);

    auto tg_xyyzz_xz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 74);

    auto tg_xyyzz_yy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 75);

    auto tg_xyyzz_yz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 76);

    auto tg_xyyzz_zz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 77);

    auto tg_xyzzz_xx_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 78);

    auto tg_xyzzz_xy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 79);

    auto tg_xyzzz_xz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 80);

    auto tg_xyzzz_yy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 81);

    auto tg_xyzzz_yz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 82);

    auto tg_xyzzz_zz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 83);

    auto tg_xzzzz_xx_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 84);

    auto tg_xzzzz_xy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 85);

    auto tg_xzzzz_xz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 86);

    auto tg_xzzzz_yy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 87);

    auto tg_xzzzz_yz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 88);

    auto tg_xzzzz_zz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 89);

    auto tg_yyyyy_xx_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 90);

    auto tg_yyyyy_xy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 91);

    auto tg_yyyyy_xz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 92);

    auto tg_yyyyy_yy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 93);

    auto tg_yyyyy_yz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 94);

    auto tg_yyyyy_zz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 95);

    auto tg_yyyyz_xx_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 96);

    auto tg_yyyyz_xy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 97);

    auto tg_yyyyz_xz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 98);

    auto tg_yyyyz_yy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 99);

    auto tg_yyyyz_yz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 100);

    auto tg_yyyyz_zz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 101);

    auto tg_yyyzz_xx_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 102);

    auto tg_yyyzz_xy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 103);

    auto tg_yyyzz_xz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 104);

    auto tg_yyyzz_yy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 105);

    auto tg_yyyzz_yz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 106);

    auto tg_yyyzz_zz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 107);

    auto tg_yyzzz_xx_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 108);

    auto tg_yyzzz_xy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 109);

    auto tg_yyzzz_xz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 110);

    auto tg_yyzzz_yy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 111);

    auto tg_yyzzz_yz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 112);

    auto tg_yyzzz_zz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 113);

    auto tg_yzzzz_xx_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 114);

    auto tg_yzzzz_xy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 115);

    auto tg_yzzzz_xz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 116);

    auto tg_yzzzz_yy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 117);

    auto tg_yzzzz_yz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 118);

    auto tg_yzzzz_zz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 119);

    auto tg_zzzzz_xx_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 120);

    auto tg_zzzzz_xy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 121);

    auto tg_zzzzz_xz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 122);

    auto tg_zzzzz_yy_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 123);

    auto tg_zzzzz_yz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 124);

    auto tg_zzzzz_zz_s_0_0_0 = pbuffer.data(idx_hd_s_0_0_0 + 125);

    // Set up components of auxiliary buffer : GD

    auto tg_xxxx_xx_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0);

    auto tg_xxxx_xy_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 1);

    auto tg_xxxx_xz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 2);

    auto tg_xxxx_yy_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 3);

    auto tg_xxxx_yz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 4);

    auto tg_xxxx_zz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 5);

    auto tg_xxxy_xx_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 6);

    auto tg_xxxy_xy_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 7);

    auto tg_xxxy_xz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 8);

    auto tg_xxxy_yy_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 9);

    auto tg_xxxy_yz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 10);

    auto tg_xxxy_zz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 11);

    auto tg_xxxz_xx_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 12);

    auto tg_xxxz_xy_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 13);

    auto tg_xxxz_xz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 14);

    auto tg_xxxz_yy_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 15);

    auto tg_xxxz_yz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 16);

    auto tg_xxxz_zz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 17);

    auto tg_xxyy_xx_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 18);

    auto tg_xxyy_xy_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 19);

    auto tg_xxyy_xz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 20);

    auto tg_xxyy_yy_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 21);

    auto tg_xxyy_yz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 22);

    auto tg_xxyy_zz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 23);

    auto tg_xxyz_xx_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 24);

    auto tg_xxyz_xy_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 25);

    auto tg_xxyz_xz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 26);

    auto tg_xxyz_yy_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 27);

    auto tg_xxyz_yz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 28);

    auto tg_xxyz_zz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 29);

    auto tg_xxzz_xx_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 30);

    auto tg_xxzz_xy_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 31);

    auto tg_xxzz_xz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 32);

    auto tg_xxzz_yy_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 33);

    auto tg_xxzz_yz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 34);

    auto tg_xxzz_zz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 35);

    auto tg_xyyy_xx_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 36);

    auto tg_xyyy_xy_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 37);

    auto tg_xyyy_xz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 38);

    auto tg_xyyy_yy_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 39);

    auto tg_xyyy_yz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 40);

    auto tg_xyyy_zz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 41);

    auto tg_xyyz_xx_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 42);

    auto tg_xyyz_xy_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 43);

    auto tg_xyyz_xz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 44);

    auto tg_xyyz_yy_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 45);

    auto tg_xyyz_yz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 46);

    auto tg_xyyz_zz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 47);

    auto tg_xyzz_xx_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 48);

    auto tg_xyzz_xy_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 49);

    auto tg_xyzz_xz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 50);

    auto tg_xyzz_yy_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 51);

    auto tg_xyzz_yz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 52);

    auto tg_xyzz_zz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 53);

    auto tg_xzzz_xx_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 54);

    auto tg_xzzz_xy_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 55);

    auto tg_xzzz_xz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 56);

    auto tg_xzzz_yy_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 57);

    auto tg_xzzz_yz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 58);

    auto tg_xzzz_zz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 59);

    auto tg_yyyy_xx_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 60);

    auto tg_yyyy_xy_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 61);

    auto tg_yyyy_xz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 62);

    auto tg_yyyy_yy_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 63);

    auto tg_yyyy_yz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 64);

    auto tg_yyyy_zz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 65);

    auto tg_yyyz_xx_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 66);

    auto tg_yyyz_xy_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 67);

    auto tg_yyyz_xz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 68);

    auto tg_yyyz_yy_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 69);

    auto tg_yyyz_yz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 70);

    auto tg_yyyz_zz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 71);

    auto tg_yyzz_xx_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 72);

    auto tg_yyzz_xy_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 73);

    auto tg_yyzz_xz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 74);

    auto tg_yyzz_yy_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 75);

    auto tg_yyzz_yz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 76);

    auto tg_yyzz_zz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 77);

    auto tg_yzzz_xx_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 78);

    auto tg_yzzz_xy_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 79);

    auto tg_yzzz_xz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 80);

    auto tg_yzzz_yy_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 81);

    auto tg_yzzz_yz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 82);

    auto tg_yzzz_zz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 83);

    auto tg_zzzz_xx_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 84);

    auto tg_zzzz_xy_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 85);

    auto tg_zzzz_xz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 86);

    auto tg_zzzz_yy_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 87);

    auto tg_zzzz_yz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 88);

    auto tg_zzzz_zz_s_1_0_0 = pbuffer.data(idx_gd_s_1_0_0 + 89);

    // Set up components of auxiliary buffer : HD

    auto tg_xxxxx_xx_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0);

    auto tg_xxxxx_xy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 1);

    auto tg_xxxxx_xz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 2);

    auto tg_xxxxx_yy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 3);

    auto tg_xxxxx_yz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 4);

    auto tg_xxxxx_zz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 5);

    auto tg_xxxxy_xx_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 6);

    auto tg_xxxxy_xy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 7);

    auto tg_xxxxy_xz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 8);

    auto tg_xxxxy_yy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 9);

    auto tg_xxxxy_yz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 10);

    auto tg_xxxxy_zz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 11);

    auto tg_xxxxz_xx_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 12);

    auto tg_xxxxz_xy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 13);

    auto tg_xxxxz_xz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 14);

    auto tg_xxxxz_yy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 15);

    auto tg_xxxxz_yz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 16);

    auto tg_xxxxz_zz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 17);

    auto tg_xxxyy_xx_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 18);

    auto tg_xxxyy_xy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 19);

    auto tg_xxxyy_xz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 20);

    auto tg_xxxyy_yy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 21);

    auto tg_xxxyy_yz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 22);

    auto tg_xxxyy_zz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 23);

    auto tg_xxxyz_xx_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 24);

    auto tg_xxxyz_xy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 25);

    auto tg_xxxyz_xz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 26);

    auto tg_xxxyz_yy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 27);

    auto tg_xxxyz_yz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 28);

    auto tg_xxxyz_zz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 29);

    auto tg_xxxzz_xx_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 30);

    auto tg_xxxzz_xy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 31);

    auto tg_xxxzz_xz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 32);

    auto tg_xxxzz_yy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 33);

    auto tg_xxxzz_yz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 34);

    auto tg_xxxzz_zz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 35);

    auto tg_xxyyy_xx_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 36);

    auto tg_xxyyy_xy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 37);

    auto tg_xxyyy_xz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 38);

    auto tg_xxyyy_yy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 39);

    auto tg_xxyyy_yz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 40);

    auto tg_xxyyy_zz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 41);

    auto tg_xxyyz_xx_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 42);

    auto tg_xxyyz_xy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 43);

    auto tg_xxyyz_xz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 44);

    auto tg_xxyyz_yy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 45);

    auto tg_xxyyz_yz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 46);

    auto tg_xxyyz_zz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 47);

    auto tg_xxyzz_xx_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 48);

    auto tg_xxyzz_xy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 49);

    auto tg_xxyzz_xz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 50);

    auto tg_xxyzz_yy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 51);

    auto tg_xxyzz_yz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 52);

    auto tg_xxyzz_zz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 53);

    auto tg_xxzzz_xx_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 54);

    auto tg_xxzzz_xy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 55);

    auto tg_xxzzz_xz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 56);

    auto tg_xxzzz_yy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 57);

    auto tg_xxzzz_yz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 58);

    auto tg_xxzzz_zz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 59);

    auto tg_xyyyy_xx_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 60);

    auto tg_xyyyy_xy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 61);

    auto tg_xyyyy_xz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 62);

    auto tg_xyyyy_yy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 63);

    auto tg_xyyyy_yz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 64);

    auto tg_xyyyy_zz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 65);

    auto tg_xyyyz_xx_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 66);

    auto tg_xyyyz_xy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 67);

    auto tg_xyyyz_xz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 68);

    auto tg_xyyyz_yy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 69);

    auto tg_xyyyz_yz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 70);

    auto tg_xyyyz_zz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 71);

    auto tg_xyyzz_xx_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 72);

    auto tg_xyyzz_xy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 73);

    auto tg_xyyzz_xz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 74);

    auto tg_xyyzz_yy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 75);

    auto tg_xyyzz_yz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 76);

    auto tg_xyyzz_zz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 77);

    auto tg_xyzzz_xx_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 78);

    auto tg_xyzzz_xy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 79);

    auto tg_xyzzz_xz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 80);

    auto tg_xyzzz_yy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 81);

    auto tg_xyzzz_yz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 82);

    auto tg_xyzzz_zz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 83);

    auto tg_xzzzz_xx_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 84);

    auto tg_xzzzz_xy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 85);

    auto tg_xzzzz_xz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 86);

    auto tg_xzzzz_yy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 87);

    auto tg_xzzzz_yz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 88);

    auto tg_xzzzz_zz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 89);

    auto tg_yyyyy_xx_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 90);

    auto tg_yyyyy_xy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 91);

    auto tg_yyyyy_xz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 92);

    auto tg_yyyyy_yy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 93);

    auto tg_yyyyy_yz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 94);

    auto tg_yyyyy_zz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 95);

    auto tg_yyyyz_xx_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 96);

    auto tg_yyyyz_xy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 97);

    auto tg_yyyyz_xz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 98);

    auto tg_yyyyz_yy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 99);

    auto tg_yyyyz_yz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 100);

    auto tg_yyyyz_zz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 101);

    auto tg_yyyzz_xx_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 102);

    auto tg_yyyzz_xy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 103);

    auto tg_yyyzz_xz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 104);

    auto tg_yyyzz_yy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 105);

    auto tg_yyyzz_yz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 106);

    auto tg_yyyzz_zz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 107);

    auto tg_yyzzz_xx_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 108);

    auto tg_yyzzz_xy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 109);

    auto tg_yyzzz_xz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 110);

    auto tg_yyzzz_yy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 111);

    auto tg_yyzzz_yz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 112);

    auto tg_yyzzz_zz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 113);

    auto tg_yzzzz_xx_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 114);

    auto tg_yzzzz_xy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 115);

    auto tg_yzzzz_xz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 116);

    auto tg_yzzzz_yy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 117);

    auto tg_yzzzz_yz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 118);

    auto tg_yzzzz_zz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 119);

    auto tg_zzzzz_xx_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 120);

    auto tg_zzzzz_xy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 121);

    auto tg_zzzzz_xz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 122);

    auto tg_zzzzz_yy_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 123);

    auto tg_zzzzz_yz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 124);

    auto tg_zzzzz_zz_s_1_0_0 = pbuffer.data(idx_hd_s_1_0_0 + 125);

    // Set up components of targeted buffer : ID

    auto tg_xxxxxx_xx_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0);

    auto tg_xxxxxx_xy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 1);

    auto tg_xxxxxx_xz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 2);

    auto tg_xxxxxx_yy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 3);

    auto tg_xxxxxx_yz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 4);

    auto tg_xxxxxx_zz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 5);

    auto tg_xxxxxy_xx_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 6);

    auto tg_xxxxxy_xy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 7);

    auto tg_xxxxxy_xz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 8);

    auto tg_xxxxxy_yy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 9);

    auto tg_xxxxxy_yz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 10);

    auto tg_xxxxxy_zz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 11);

    auto tg_xxxxxz_xx_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 12);

    auto tg_xxxxxz_xy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 13);

    auto tg_xxxxxz_xz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 14);

    auto tg_xxxxxz_yy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 15);

    auto tg_xxxxxz_yz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 16);

    auto tg_xxxxxz_zz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 17);

    auto tg_xxxxyy_xx_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 18);

    auto tg_xxxxyy_xy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 19);

    auto tg_xxxxyy_xz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 20);

    auto tg_xxxxyy_yy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 21);

    auto tg_xxxxyy_yz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 22);

    auto tg_xxxxyy_zz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 23);

    auto tg_xxxxyz_xx_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 24);

    auto tg_xxxxyz_xy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 25);

    auto tg_xxxxyz_xz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 26);

    auto tg_xxxxyz_yy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 27);

    auto tg_xxxxyz_yz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 28);

    auto tg_xxxxyz_zz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 29);

    auto tg_xxxxzz_xx_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 30);

    auto tg_xxxxzz_xy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 31);

    auto tg_xxxxzz_xz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 32);

    auto tg_xxxxzz_yy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 33);

    auto tg_xxxxzz_yz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 34);

    auto tg_xxxxzz_zz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 35);

    auto tg_xxxyyy_xx_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 36);

    auto tg_xxxyyy_xy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 37);

    auto tg_xxxyyy_xz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 38);

    auto tg_xxxyyy_yy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 39);

    auto tg_xxxyyy_yz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 40);

    auto tg_xxxyyy_zz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 41);

    auto tg_xxxyyz_xx_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 42);

    auto tg_xxxyyz_xy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 43);

    auto tg_xxxyyz_xz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 44);

    auto tg_xxxyyz_yy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 45);

    auto tg_xxxyyz_yz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 46);

    auto tg_xxxyyz_zz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 47);

    auto tg_xxxyzz_xx_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 48);

    auto tg_xxxyzz_xy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 49);

    auto tg_xxxyzz_xz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 50);

    auto tg_xxxyzz_yy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 51);

    auto tg_xxxyzz_yz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 52);

    auto tg_xxxyzz_zz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 53);

    auto tg_xxxzzz_xx_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 54);

    auto tg_xxxzzz_xy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 55);

    auto tg_xxxzzz_xz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 56);

    auto tg_xxxzzz_yy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 57);

    auto tg_xxxzzz_yz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 58);

    auto tg_xxxzzz_zz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 59);

    auto tg_xxyyyy_xx_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 60);

    auto tg_xxyyyy_xy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 61);

    auto tg_xxyyyy_xz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 62);

    auto tg_xxyyyy_yy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 63);

    auto tg_xxyyyy_yz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 64);

    auto tg_xxyyyy_zz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 65);

    auto tg_xxyyyz_xx_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 66);

    auto tg_xxyyyz_xy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 67);

    auto tg_xxyyyz_xz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 68);

    auto tg_xxyyyz_yy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 69);

    auto tg_xxyyyz_yz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 70);

    auto tg_xxyyyz_zz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 71);

    auto tg_xxyyzz_xx_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 72);

    auto tg_xxyyzz_xy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 73);

    auto tg_xxyyzz_xz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 74);

    auto tg_xxyyzz_yy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 75);

    auto tg_xxyyzz_yz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 76);

    auto tg_xxyyzz_zz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 77);

    auto tg_xxyzzz_xx_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 78);

    auto tg_xxyzzz_xy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 79);

    auto tg_xxyzzz_xz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 80);

    auto tg_xxyzzz_yy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 81);

    auto tg_xxyzzz_yz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 82);

    auto tg_xxyzzz_zz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 83);

    auto tg_xxzzzz_xx_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 84);

    auto tg_xxzzzz_xy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 85);

    auto tg_xxzzzz_xz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 86);

    auto tg_xxzzzz_yy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 87);

    auto tg_xxzzzz_yz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 88);

    auto tg_xxzzzz_zz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 89);

    auto tg_xyyyyy_xx_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 90);

    auto tg_xyyyyy_xy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 91);

    auto tg_xyyyyy_xz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 92);

    auto tg_xyyyyy_yy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 93);

    auto tg_xyyyyy_yz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 94);

    auto tg_xyyyyy_zz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 95);

    auto tg_xyyyyz_xx_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 96);

    auto tg_xyyyyz_xy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 97);

    auto tg_xyyyyz_xz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 98);

    auto tg_xyyyyz_yy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 99);

    auto tg_xyyyyz_yz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 100);

    auto tg_xyyyyz_zz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 101);

    auto tg_xyyyzz_xx_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 102);

    auto tg_xyyyzz_xy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 103);

    auto tg_xyyyzz_xz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 104);

    auto tg_xyyyzz_yy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 105);

    auto tg_xyyyzz_yz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 106);

    auto tg_xyyyzz_zz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 107);

    auto tg_xyyzzz_xx_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 108);

    auto tg_xyyzzz_xy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 109);

    auto tg_xyyzzz_xz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 110);

    auto tg_xyyzzz_yy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 111);

    auto tg_xyyzzz_yz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 112);

    auto tg_xyyzzz_zz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 113);

    auto tg_xyzzzz_xx_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 114);

    auto tg_xyzzzz_xy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 115);

    auto tg_xyzzzz_xz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 116);

    auto tg_xyzzzz_yy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 117);

    auto tg_xyzzzz_yz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 118);

    auto tg_xyzzzz_zz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 119);

    auto tg_xzzzzz_xx_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 120);

    auto tg_xzzzzz_xy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 121);

    auto tg_xzzzzz_xz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 122);

    auto tg_xzzzzz_yy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 123);

    auto tg_xzzzzz_yz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 124);

    auto tg_xzzzzz_zz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 125);

    auto tg_yyyyyy_xx_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 126);

    auto tg_yyyyyy_xy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 127);

    auto tg_yyyyyy_xz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 128);

    auto tg_yyyyyy_yy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 129);

    auto tg_yyyyyy_yz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 130);

    auto tg_yyyyyy_zz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 131);

    auto tg_yyyyyz_xx_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 132);

    auto tg_yyyyyz_xy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 133);

    auto tg_yyyyyz_xz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 134);

    auto tg_yyyyyz_yy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 135);

    auto tg_yyyyyz_yz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 136);

    auto tg_yyyyyz_zz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 137);

    auto tg_yyyyzz_xx_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 138);

    auto tg_yyyyzz_xy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 139);

    auto tg_yyyyzz_xz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 140);

    auto tg_yyyyzz_yy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 141);

    auto tg_yyyyzz_yz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 142);

    auto tg_yyyyzz_zz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 143);

    auto tg_yyyzzz_xx_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 144);

    auto tg_yyyzzz_xy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 145);

    auto tg_yyyzzz_xz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 146);

    auto tg_yyyzzz_yy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 147);

    auto tg_yyyzzz_yz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 148);

    auto tg_yyyzzz_zz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 149);

    auto tg_yyzzzz_xx_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 150);

    auto tg_yyzzzz_xy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 151);

    auto tg_yyzzzz_xz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 152);

    auto tg_yyzzzz_yy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 153);

    auto tg_yyzzzz_yz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 154);

    auto tg_yyzzzz_zz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 155);

    auto tg_yzzzzz_xx_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 156);

    auto tg_yzzzzz_xy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 157);

    auto tg_yzzzzz_xz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 158);

    auto tg_yzzzzz_yy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 159);

    auto tg_yzzzzz_yz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 160);

    auto tg_yzzzzz_zz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 161);

    auto tg_zzzzzz_xx_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 162);

    auto tg_zzzzzz_xy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 163);

    auto tg_zzzzzz_xz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 164);

    auto tg_zzzzzz_yy_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 165);

    auto tg_zzzzzz_yz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 166);

    auto tg_zzzzzz_zz_s_0_0_0 = pbuffer.data(idx_id_s_0_0_0 + 167);

    #pragma omp simd aligned(b_exps, tg_xxxx_xx_s_0_0_0, tg_xxxx_xx_s_1_0_0, tg_xxxx_xy_s_0_0_0, tg_xxxx_xy_s_1_0_0, tg_xxxx_xz_s_0_0_0, tg_xxxx_xz_s_1_0_0, tg_xxxx_yy_s_0_0_0, tg_xxxx_yy_s_1_0_0, tg_xxxx_yz_s_0_0_0, tg_xxxx_yz_s_1_0_0, tg_xxxx_zz_s_0_0_0, tg_xxxx_zz_s_1_0_0, tg_xxxxx_xx_s_0_0_0, tg_xxxxx_xx_s_1_0_0, tg_xxxxx_xy_s_0_0_0, tg_xxxxx_xy_s_1_0_0, tg_xxxxx_xz_s_0_0_0, tg_xxxxx_xz_s_1_0_0, tg_xxxxx_yy_s_0_0_0, tg_xxxxx_yy_s_1_0_0, tg_xxxxx_yz_s_0_0_0, tg_xxxxx_yz_s_1_0_0, tg_xxxxx_zz_s_0_0_0, tg_xxxxx_zz_s_1_0_0, tg_xxxxxx_xx_s_0_0_0, tg_xxxxxx_xy_s_0_0_0, tg_xxxxxx_xz_s_0_0_0, tg_xxxxxx_yy_s_0_0_0, tg_xxxxxx_yz_s_0_0_0, tg_xxxxxx_zz_s_0_0_0, tg_xxxxxy_xx_s_0_0_0, tg_xxxxxy_xy_s_0_0_0, tg_xxxxxy_xz_s_0_0_0, tg_xxxxxy_yy_s_0_0_0, tg_xxxxxy_yz_s_0_0_0, tg_xxxxxy_zz_s_0_0_0, tg_xxxxxz_xx_s_0_0_0, tg_xxxxxz_xy_s_0_0_0, tg_xxxxxz_xz_s_0_0_0, tg_xxxxxz_yy_s_0_0_0, tg_xxxxxz_yz_s_0_0_0, tg_xxxxxz_zz_s_0_0_0, tg_xxxxyy_xx_s_0_0_0, tg_xxxxyy_xy_s_0_0_0, tg_xxxxyy_xz_s_0_0_0, tg_xxxxyy_yy_s_0_0_0, tg_xxxxyy_yz_s_0_0_0, tg_xxxxyy_zz_s_0_0_0, tg_xxxxyz_xx_s_0_0_0, tg_xxxxyz_xy_s_0_0_0, tg_xxxxyz_xz_s_0_0_0, tg_xxxxyz_yy_s_0_0_0, tg_xxxxyz_yz_s_0_0_0, tg_xxxxyz_zz_s_0_0_0, tg_xxxxz_xx_s_0_0_0, tg_xxxxz_xx_s_1_0_0, tg_xxxxz_xy_s_0_0_0, tg_xxxxz_xy_s_1_0_0, tg_xxxxz_xz_s_0_0_0, tg_xxxxz_xz_s_1_0_0, tg_xxxxz_yy_s_0_0_0, tg_xxxxz_yy_s_1_0_0, tg_xxxxz_yz_s_0_0_0, tg_xxxxz_yz_s_1_0_0, tg_xxxxz_zz_s_0_0_0, tg_xxxxz_zz_s_1_0_0, tg_xxxxzz_xx_s_0_0_0, tg_xxxxzz_xy_s_0_0_0, tg_xxxxzz_xz_s_0_0_0, tg_xxxxzz_yy_s_0_0_0, tg_xxxxzz_yz_s_0_0_0, tg_xxxxzz_zz_s_0_0_0, tg_xxxyy_xx_s_0_0_0, tg_xxxyy_xx_s_1_0_0, tg_xxxyy_xy_s_0_0_0, tg_xxxyy_xy_s_1_0_0, tg_xxxyy_xz_s_0_0_0, tg_xxxyy_xz_s_1_0_0, tg_xxxyy_yy_s_0_0_0, tg_xxxyy_yy_s_1_0_0, tg_xxxyy_yz_s_0_0_0, tg_xxxyy_yz_s_1_0_0, tg_xxxyy_zz_s_0_0_0, tg_xxxyy_zz_s_1_0_0, tg_xxxyyy_xx_s_0_0_0, tg_xxxyyy_xy_s_0_0_0, tg_xxxyyy_xz_s_0_0_0, tg_xxxyyy_yy_s_0_0_0, tg_xxxyyy_yz_s_0_0_0, tg_xxxyyy_zz_s_0_0_0, tg_xxxyyz_xx_s_0_0_0, tg_xxxyyz_xy_s_0_0_0, tg_xxxyyz_xz_s_0_0_0, tg_xxxyyz_yy_s_0_0_0, tg_xxxyyz_yz_s_0_0_0, tg_xxxyyz_zz_s_0_0_0, tg_xxxyzz_xx_s_0_0_0, tg_xxxyzz_xy_s_0_0_0, tg_xxxyzz_xz_s_0_0_0, tg_xxxyzz_yy_s_0_0_0, tg_xxxyzz_yz_s_0_0_0, tg_xxxyzz_zz_s_0_0_0, tg_xxxzz_xx_s_0_0_0, tg_xxxzz_xx_s_1_0_0, tg_xxxzz_xy_s_0_0_0, tg_xxxzz_xy_s_1_0_0, tg_xxxzz_xz_s_0_0_0, tg_xxxzz_xz_s_1_0_0, tg_xxxzz_yy_s_0_0_0, tg_xxxzz_yy_s_1_0_0, tg_xxxzz_yz_s_0_0_0, tg_xxxzz_yz_s_1_0_0, tg_xxxzz_zz_s_0_0_0, tg_xxxzz_zz_s_1_0_0, tg_xxxzzz_xx_s_0_0_0, tg_xxxzzz_xy_s_0_0_0, tg_xxxzzz_xz_s_0_0_0, tg_xxxzzz_yy_s_0_0_0, tg_xxxzzz_yz_s_0_0_0, tg_xxxzzz_zz_s_0_0_0, tg_xxyy_xx_s_0_0_0, tg_xxyy_xx_s_1_0_0, tg_xxyy_xy_s_0_0_0, tg_xxyy_xy_s_1_0_0, tg_xxyy_xz_s_0_0_0, tg_xxyy_xz_s_1_0_0, tg_xxyy_yy_s_0_0_0, tg_xxyy_yy_s_1_0_0, tg_xxyy_yz_s_0_0_0, tg_xxyy_yz_s_1_0_0, tg_xxyy_zz_s_0_0_0, tg_xxyy_zz_s_1_0_0, tg_xxyyy_xx_s_0_0_0, tg_xxyyy_xx_s_1_0_0, tg_xxyyy_xy_s_0_0_0, tg_xxyyy_xy_s_1_0_0, tg_xxyyy_xz_s_0_0_0, tg_xxyyy_xz_s_1_0_0, tg_xxyyy_yy_s_0_0_0, tg_xxyyy_yy_s_1_0_0, tg_xxyyy_yz_s_0_0_0, tg_xxyyy_yz_s_1_0_0, tg_xxyyy_zz_s_0_0_0, tg_xxyyy_zz_s_1_0_0, tg_xxyyyy_xx_s_0_0_0, tg_xxyyyy_xy_s_0_0_0, tg_xxyyyy_xz_s_0_0_0, tg_xxyyyy_yy_s_0_0_0, tg_xxyyyy_yz_s_0_0_0, tg_xxyyyy_zz_s_0_0_0, tg_xxyyyz_xx_s_0_0_0, tg_xxyyyz_xy_s_0_0_0, tg_xxyyyz_xz_s_0_0_0, tg_xxyyyz_yy_s_0_0_0, tg_xxyyyz_yz_s_0_0_0, tg_xxyyyz_zz_s_0_0_0, tg_xxyyzz_xx_s_0_0_0, tg_xxyyzz_xy_s_0_0_0, tg_xxyyzz_xz_s_0_0_0, tg_xxyyzz_yy_s_0_0_0, tg_xxyyzz_yz_s_0_0_0, tg_xxyyzz_zz_s_0_0_0, tg_xxyzzz_xx_s_0_0_0, tg_xxyzzz_xy_s_0_0_0, tg_xxyzzz_xz_s_0_0_0, tg_xxyzzz_yy_s_0_0_0, tg_xxyzzz_yz_s_0_0_0, tg_xxyzzz_zz_s_0_0_0, tg_xxzz_xx_s_0_0_0, tg_xxzz_xx_s_1_0_0, tg_xxzz_xy_s_0_0_0, tg_xxzz_xy_s_1_0_0, tg_xxzz_xz_s_0_0_0, tg_xxzz_xz_s_1_0_0, tg_xxzz_yy_s_0_0_0, tg_xxzz_yy_s_1_0_0, tg_xxzz_yz_s_0_0_0, tg_xxzz_yz_s_1_0_0, tg_xxzz_zz_s_0_0_0, tg_xxzz_zz_s_1_0_0, tg_xxzzz_xx_s_0_0_0, tg_xxzzz_xx_s_1_0_0, tg_xxzzz_xy_s_0_0_0, tg_xxzzz_xy_s_1_0_0, tg_xxzzz_xz_s_0_0_0, tg_xxzzz_xz_s_1_0_0, tg_xxzzz_yy_s_0_0_0, tg_xxzzz_yy_s_1_0_0, tg_xxzzz_yz_s_0_0_0, tg_xxzzz_yz_s_1_0_0, tg_xxzzz_zz_s_0_0_0, tg_xxzzz_zz_s_1_0_0, tg_xxzzzz_xx_s_0_0_0, tg_xxzzzz_xy_s_0_0_0, tg_xxzzzz_xz_s_0_0_0, tg_xxzzzz_yy_s_0_0_0, tg_xxzzzz_yz_s_0_0_0, tg_xxzzzz_zz_s_0_0_0, tg_xyyy_xx_s_0_0_0, tg_xyyy_xx_s_1_0_0, tg_xyyy_xy_s_0_0_0, tg_xyyy_xy_s_1_0_0, tg_xyyy_xz_s_0_0_0, tg_xyyy_xz_s_1_0_0, tg_xyyy_yy_s_0_0_0, tg_xyyy_yy_s_1_0_0, tg_xyyy_yz_s_0_0_0, tg_xyyy_yz_s_1_0_0, tg_xyyy_zz_s_0_0_0, tg_xyyy_zz_s_1_0_0, tg_xyyyy_xx_s_0_0_0, tg_xyyyy_xx_s_1_0_0, tg_xyyyy_xy_s_0_0_0, tg_xyyyy_xy_s_1_0_0, tg_xyyyy_xz_s_0_0_0, tg_xyyyy_xz_s_1_0_0, tg_xyyyy_yy_s_0_0_0, tg_xyyyy_yy_s_1_0_0, tg_xyyyy_yz_s_0_0_0, tg_xyyyy_yz_s_1_0_0, tg_xyyyy_zz_s_0_0_0, tg_xyyyy_zz_s_1_0_0, tg_xyyyyy_xx_s_0_0_0, tg_xyyyyy_xy_s_0_0_0, tg_xyyyyy_xz_s_0_0_0, tg_xyyyyy_yy_s_0_0_0, tg_xyyyyy_yz_s_0_0_0, tg_xyyyyy_zz_s_0_0_0, tg_xyyyyz_xx_s_0_0_0, tg_xyyyyz_xy_s_0_0_0, tg_xyyyyz_xz_s_0_0_0, tg_xyyyyz_yy_s_0_0_0, tg_xyyyyz_yz_s_0_0_0, tg_xyyyyz_zz_s_0_0_0, tg_xyyyzz_xx_s_0_0_0, tg_xyyyzz_xy_s_0_0_0, tg_xyyyzz_xz_s_0_0_0, tg_xyyyzz_yy_s_0_0_0, tg_xyyyzz_yz_s_0_0_0, tg_xyyyzz_zz_s_0_0_0, tg_xyyzz_xx_s_0_0_0, tg_xyyzz_xx_s_1_0_0, tg_xyyzz_xy_s_0_0_0, tg_xyyzz_xy_s_1_0_0, tg_xyyzz_xz_s_0_0_0, tg_xyyzz_xz_s_1_0_0, tg_xyyzz_yy_s_0_0_0, tg_xyyzz_yy_s_1_0_0, tg_xyyzz_yz_s_0_0_0, tg_xyyzz_yz_s_1_0_0, tg_xyyzz_zz_s_0_0_0, tg_xyyzz_zz_s_1_0_0, tg_xyyzzz_xx_s_0_0_0, tg_xyyzzz_xy_s_0_0_0, tg_xyyzzz_xz_s_0_0_0, tg_xyyzzz_yy_s_0_0_0, tg_xyyzzz_yz_s_0_0_0, tg_xyyzzz_zz_s_0_0_0, tg_xyzzzz_xx_s_0_0_0, tg_xyzzzz_xy_s_0_0_0, tg_xyzzzz_xz_s_0_0_0, tg_xyzzzz_yy_s_0_0_0, tg_xyzzzz_yz_s_0_0_0, tg_xyzzzz_zz_s_0_0_0, tg_xzzz_xx_s_0_0_0, tg_xzzz_xx_s_1_0_0, tg_xzzz_xy_s_0_0_0, tg_xzzz_xy_s_1_0_0, tg_xzzz_xz_s_0_0_0, tg_xzzz_xz_s_1_0_0, tg_xzzz_yy_s_0_0_0, tg_xzzz_yy_s_1_0_0, tg_xzzz_yz_s_0_0_0, tg_xzzz_yz_s_1_0_0, tg_xzzz_zz_s_0_0_0, tg_xzzz_zz_s_1_0_0, tg_xzzzz_xx_s_0_0_0, tg_xzzzz_xx_s_1_0_0, tg_xzzzz_xy_s_0_0_0, tg_xzzzz_xy_s_1_0_0, tg_xzzzz_xz_s_0_0_0, tg_xzzzz_xz_s_1_0_0, tg_xzzzz_yy_s_0_0_0, tg_xzzzz_yy_s_1_0_0, tg_xzzzz_yz_s_0_0_0, tg_xzzzz_yz_s_1_0_0, tg_xzzzz_zz_s_0_0_0, tg_xzzzz_zz_s_1_0_0, tg_xzzzzz_xx_s_0_0_0, tg_xzzzzz_xy_s_0_0_0, tg_xzzzzz_xz_s_0_0_0, tg_xzzzzz_yy_s_0_0_0, tg_xzzzzz_yz_s_0_0_0, tg_xzzzzz_zz_s_0_0_0, tg_yyyy_xx_s_0_0_0, tg_yyyy_xx_s_1_0_0, tg_yyyy_xy_s_0_0_0, tg_yyyy_xy_s_1_0_0, tg_yyyy_xz_s_0_0_0, tg_yyyy_xz_s_1_0_0, tg_yyyy_yy_s_0_0_0, tg_yyyy_yy_s_1_0_0, tg_yyyy_yz_s_0_0_0, tg_yyyy_yz_s_1_0_0, tg_yyyy_zz_s_0_0_0, tg_yyyy_zz_s_1_0_0, tg_yyyyy_xx_s_0_0_0, tg_yyyyy_xx_s_1_0_0, tg_yyyyy_xy_s_0_0_0, tg_yyyyy_xy_s_1_0_0, tg_yyyyy_xz_s_0_0_0, tg_yyyyy_xz_s_1_0_0, tg_yyyyy_yy_s_0_0_0, tg_yyyyy_yy_s_1_0_0, tg_yyyyy_yz_s_0_0_0, tg_yyyyy_yz_s_1_0_0, tg_yyyyy_zz_s_0_0_0, tg_yyyyy_zz_s_1_0_0, tg_yyyyyy_xx_s_0_0_0, tg_yyyyyy_xy_s_0_0_0, tg_yyyyyy_xz_s_0_0_0, tg_yyyyyy_yy_s_0_0_0, tg_yyyyyy_yz_s_0_0_0, tg_yyyyyy_zz_s_0_0_0, tg_yyyyyz_xx_s_0_0_0, tg_yyyyyz_xy_s_0_0_0, tg_yyyyyz_xz_s_0_0_0, tg_yyyyyz_yy_s_0_0_0, tg_yyyyyz_yz_s_0_0_0, tg_yyyyyz_zz_s_0_0_0, tg_yyyyz_xx_s_0_0_0, tg_yyyyz_xx_s_1_0_0, tg_yyyyz_xy_s_0_0_0, tg_yyyyz_xy_s_1_0_0, tg_yyyyz_xz_s_0_0_0, tg_yyyyz_xz_s_1_0_0, tg_yyyyz_yy_s_0_0_0, tg_yyyyz_yy_s_1_0_0, tg_yyyyz_yz_s_0_0_0, tg_yyyyz_yz_s_1_0_0, tg_yyyyz_zz_s_0_0_0, tg_yyyyz_zz_s_1_0_0, tg_yyyyzz_xx_s_0_0_0, tg_yyyyzz_xy_s_0_0_0, tg_yyyyzz_xz_s_0_0_0, tg_yyyyzz_yy_s_0_0_0, tg_yyyyzz_yz_s_0_0_0, tg_yyyyzz_zz_s_0_0_0, tg_yyyzz_xx_s_0_0_0, tg_yyyzz_xx_s_1_0_0, tg_yyyzz_xy_s_0_0_0, tg_yyyzz_xy_s_1_0_0, tg_yyyzz_xz_s_0_0_0, tg_yyyzz_xz_s_1_0_0, tg_yyyzz_yy_s_0_0_0, tg_yyyzz_yy_s_1_0_0, tg_yyyzz_yz_s_0_0_0, tg_yyyzz_yz_s_1_0_0, tg_yyyzz_zz_s_0_0_0, tg_yyyzz_zz_s_1_0_0, tg_yyyzzz_xx_s_0_0_0, tg_yyyzzz_xy_s_0_0_0, tg_yyyzzz_xz_s_0_0_0, tg_yyyzzz_yy_s_0_0_0, tg_yyyzzz_yz_s_0_0_0, tg_yyyzzz_zz_s_0_0_0, tg_yyzz_xx_s_0_0_0, tg_yyzz_xx_s_1_0_0, tg_yyzz_xy_s_0_0_0, tg_yyzz_xy_s_1_0_0, tg_yyzz_xz_s_0_0_0, tg_yyzz_xz_s_1_0_0, tg_yyzz_yy_s_0_0_0, tg_yyzz_yy_s_1_0_0, tg_yyzz_yz_s_0_0_0, tg_yyzz_yz_s_1_0_0, tg_yyzz_zz_s_0_0_0, tg_yyzz_zz_s_1_0_0, tg_yyzzz_xx_s_0_0_0, tg_yyzzz_xx_s_1_0_0, tg_yyzzz_xy_s_0_0_0, tg_yyzzz_xy_s_1_0_0, tg_yyzzz_xz_s_0_0_0, tg_yyzzz_xz_s_1_0_0, tg_yyzzz_yy_s_0_0_0, tg_yyzzz_yy_s_1_0_0, tg_yyzzz_yz_s_0_0_0, tg_yyzzz_yz_s_1_0_0, tg_yyzzz_zz_s_0_0_0, tg_yyzzz_zz_s_1_0_0, tg_yyzzzz_xx_s_0_0_0, tg_yyzzzz_xy_s_0_0_0, tg_yyzzzz_xz_s_0_0_0, tg_yyzzzz_yy_s_0_0_0, tg_yyzzzz_yz_s_0_0_0, tg_yyzzzz_zz_s_0_0_0, tg_yzzz_xx_s_0_0_0, tg_yzzz_xx_s_1_0_0, tg_yzzz_xy_s_0_0_0, tg_yzzz_xy_s_1_0_0, tg_yzzz_xz_s_0_0_0, tg_yzzz_xz_s_1_0_0, tg_yzzz_yy_s_0_0_0, tg_yzzz_yy_s_1_0_0, tg_yzzz_yz_s_0_0_0, tg_yzzz_yz_s_1_0_0, tg_yzzz_zz_s_0_0_0, tg_yzzz_zz_s_1_0_0, tg_yzzzz_xx_s_0_0_0, tg_yzzzz_xx_s_1_0_0, tg_yzzzz_xy_s_0_0_0, tg_yzzzz_xy_s_1_0_0, tg_yzzzz_xz_s_0_0_0, tg_yzzzz_xz_s_1_0_0, tg_yzzzz_yy_s_0_0_0, tg_yzzzz_yy_s_1_0_0, tg_yzzzz_yz_s_0_0_0, tg_yzzzz_yz_s_1_0_0, tg_yzzzz_zz_s_0_0_0, tg_yzzzz_zz_s_1_0_0, tg_yzzzzz_xx_s_0_0_0, tg_yzzzzz_xy_s_0_0_0, tg_yzzzzz_xz_s_0_0_0, tg_yzzzzz_yy_s_0_0_0, tg_yzzzzz_yz_s_0_0_0, tg_yzzzzz_zz_s_0_0_0, tg_zzzz_xx_s_0_0_0, tg_zzzz_xx_s_1_0_0, tg_zzzz_xy_s_0_0_0, tg_zzzz_xy_s_1_0_0, tg_zzzz_xz_s_0_0_0, tg_zzzz_xz_s_1_0_0, tg_zzzz_yy_s_0_0_0, tg_zzzz_yy_s_1_0_0, tg_zzzz_yz_s_0_0_0, tg_zzzz_yz_s_1_0_0, tg_zzzz_zz_s_0_0_0, tg_zzzz_zz_s_1_0_0, tg_zzzzz_xx_s_0_0_0, tg_zzzzz_xx_s_1_0_0, tg_zzzzz_xy_s_0_0_0, tg_zzzzz_xy_s_1_0_0, tg_zzzzz_xz_s_0_0_0, tg_zzzzz_xz_s_1_0_0, tg_zzzzz_yy_s_0_0_0, tg_zzzzz_yy_s_1_0_0, tg_zzzzz_yz_s_0_0_0, tg_zzzzz_yz_s_1_0_0, tg_zzzzz_zz_s_0_0_0, tg_zzzzz_zz_s_1_0_0, tg_zzzzzz_xx_s_0_0_0, tg_zzzzzz_xy_s_0_0_0, tg_zzzzzz_xz_s_0_0_0, tg_zzzzzz_yy_s_0_0_0, tg_zzzzzz_yz_s_0_0_0, tg_zzzzzz_zz_s_0_0_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fzi_0 = 1.0 / (a_exp + b_exps[i] + c_exp);

        const double faz_0 = -(b_exps[i] + c_exp) * fzi_0;

        const double fbzi_0 = b_exps[i] * fzi_0;

        tg_xxxxxx_xx_s_0_0_0[i] = 5.0 / 2.0 * tg_xxxx_xx_s_0_0_0[i] * fzi_0 + 5.0 * tg_xxxx_xx_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xxxxx_xx_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxxxx_xx_s_0_0_0[i] * a_x * faz_0;

        tg_xxxxxx_xy_s_0_0_0[i] = 5.0 / 2.0 * tg_xxxx_xy_s_0_0_0[i] * fzi_0 + 5.0 * tg_xxxx_xy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xxxxx_xy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxxxx_xy_s_0_0_0[i] * a_x * faz_0;

        tg_xxxxxx_xz_s_0_0_0[i] = 5.0 / 2.0 * tg_xxxx_xz_s_0_0_0[i] * fzi_0 + 5.0 * tg_xxxx_xz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xxxxx_xz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxxxx_xz_s_0_0_0[i] * a_x * faz_0;

        tg_xxxxxx_yy_s_0_0_0[i] = 5.0 / 2.0 * tg_xxxx_yy_s_0_0_0[i] * fzi_0 + 5.0 * tg_xxxx_yy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xxxxx_yy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxxxx_yy_s_0_0_0[i] * a_x * faz_0;

        tg_xxxxxx_yz_s_0_0_0[i] = 5.0 / 2.0 * tg_xxxx_yz_s_0_0_0[i] * fzi_0 + 5.0 * tg_xxxx_yz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xxxxx_yz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxxxx_yz_s_0_0_0[i] * a_x * faz_0;

        tg_xxxxxx_zz_s_0_0_0[i] = 5.0 / 2.0 * tg_xxxx_zz_s_0_0_0[i] * fzi_0 + 5.0 * tg_xxxx_zz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xxxxx_zz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxxxx_zz_s_0_0_0[i] * a_x * faz_0;

        tg_xxxxxy_xx_s_0_0_0[i] = 2.0 * tg_xxxxx_xx_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxxx_xx_s_0_0_0[i] * a_y * faz_0;

        tg_xxxxxy_xy_s_0_0_0[i] = 2.0 * tg_xxxxx_xy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxxx_xy_s_0_0_0[i] * a_y * faz_0;

        tg_xxxxxy_xz_s_0_0_0[i] = 2.0 * tg_xxxxx_xz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxxx_xz_s_0_0_0[i] * a_y * faz_0;

        tg_xxxxxy_yy_s_0_0_0[i] = 2.0 * tg_xxxxx_yy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxxx_yy_s_0_0_0[i] * a_y * faz_0;

        tg_xxxxxy_yz_s_0_0_0[i] = 2.0 * tg_xxxxx_yz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxxx_yz_s_0_0_0[i] * a_y * faz_0;

        tg_xxxxxy_zz_s_0_0_0[i] = 2.0 * tg_xxxxx_zz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxxx_zz_s_0_0_0[i] * a_y * faz_0;

        tg_xxxxxz_xx_s_0_0_0[i] = 2.0 * tg_xxxxx_xx_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxxxx_xx_s_0_0_0[i] * a_z * faz_0;

        tg_xxxxxz_xy_s_0_0_0[i] = 2.0 * tg_xxxxx_xy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxxxx_xy_s_0_0_0[i] * a_z * faz_0;

        tg_xxxxxz_xz_s_0_0_0[i] = 2.0 * tg_xxxxx_xz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxxxx_xz_s_0_0_0[i] * a_z * faz_0;

        tg_xxxxxz_yy_s_0_0_0[i] = 2.0 * tg_xxxxx_yy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxxxx_yy_s_0_0_0[i] * a_z * faz_0;

        tg_xxxxxz_yz_s_0_0_0[i] = 2.0 * tg_xxxxx_yz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxxxx_yz_s_0_0_0[i] * a_z * faz_0;

        tg_xxxxxz_zz_s_0_0_0[i] = 2.0 * tg_xxxxx_zz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxxxx_zz_s_0_0_0[i] * a_z * faz_0;

        tg_xxxxyy_xx_s_0_0_0[i] = 3.0 / 2.0 * tg_xxyy_xx_s_0_0_0[i] * fzi_0 + 3.0 * tg_xxyy_xx_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xxxyy_xx_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxxyy_xx_s_0_0_0[i] * a_x * faz_0;

        tg_xxxxyy_xy_s_0_0_0[i] = 3.0 / 2.0 * tg_xxyy_xy_s_0_0_0[i] * fzi_0 + 3.0 * tg_xxyy_xy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xxxyy_xy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxxyy_xy_s_0_0_0[i] * a_x * faz_0;

        tg_xxxxyy_xz_s_0_0_0[i] = 3.0 / 2.0 * tg_xxyy_xz_s_0_0_0[i] * fzi_0 + 3.0 * tg_xxyy_xz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xxxyy_xz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxxyy_xz_s_0_0_0[i] * a_x * faz_0;

        tg_xxxxyy_yy_s_0_0_0[i] = 3.0 / 2.0 * tg_xxyy_yy_s_0_0_0[i] * fzi_0 + 3.0 * tg_xxyy_yy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xxxyy_yy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxxyy_yy_s_0_0_0[i] * a_x * faz_0;

        tg_xxxxyy_yz_s_0_0_0[i] = 3.0 / 2.0 * tg_xxyy_yz_s_0_0_0[i] * fzi_0 + 3.0 * tg_xxyy_yz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xxxyy_yz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxxyy_yz_s_0_0_0[i] * a_x * faz_0;

        tg_xxxxyy_zz_s_0_0_0[i] = 3.0 / 2.0 * tg_xxyy_zz_s_0_0_0[i] * fzi_0 + 3.0 * tg_xxyy_zz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xxxyy_zz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxxyy_zz_s_0_0_0[i] * a_x * faz_0;

        tg_xxxxyz_xx_s_0_0_0[i] = 2.0 * tg_xxxxz_xx_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxxz_xx_s_0_0_0[i] * a_y * faz_0;

        tg_xxxxyz_xy_s_0_0_0[i] = 2.0 * tg_xxxxz_xy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxxz_xy_s_0_0_0[i] * a_y * faz_0;

        tg_xxxxyz_xz_s_0_0_0[i] = 2.0 * tg_xxxxz_xz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxxz_xz_s_0_0_0[i] * a_y * faz_0;

        tg_xxxxyz_yy_s_0_0_0[i] = 2.0 * tg_xxxxz_yy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxxz_yy_s_0_0_0[i] * a_y * faz_0;

        tg_xxxxyz_yz_s_0_0_0[i] = 2.0 * tg_xxxxz_yz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxxz_yz_s_0_0_0[i] * a_y * faz_0;

        tg_xxxxyz_zz_s_0_0_0[i] = 2.0 * tg_xxxxz_zz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxxz_zz_s_0_0_0[i] * a_y * faz_0;

        tg_xxxxzz_xx_s_0_0_0[i] = 3.0 / 2.0 * tg_xxzz_xx_s_0_0_0[i] * fzi_0 + 3.0 * tg_xxzz_xx_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xxxzz_xx_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxxzz_xx_s_0_0_0[i] * a_x * faz_0;

        tg_xxxxzz_xy_s_0_0_0[i] = 3.0 / 2.0 * tg_xxzz_xy_s_0_0_0[i] * fzi_0 + 3.0 * tg_xxzz_xy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xxxzz_xy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxxzz_xy_s_0_0_0[i] * a_x * faz_0;

        tg_xxxxzz_xz_s_0_0_0[i] = 3.0 / 2.0 * tg_xxzz_xz_s_0_0_0[i] * fzi_0 + 3.0 * tg_xxzz_xz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xxxzz_xz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxxzz_xz_s_0_0_0[i] * a_x * faz_0;

        tg_xxxxzz_yy_s_0_0_0[i] = 3.0 / 2.0 * tg_xxzz_yy_s_0_0_0[i] * fzi_0 + 3.0 * tg_xxzz_yy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xxxzz_yy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxxzz_yy_s_0_0_0[i] * a_x * faz_0;

        tg_xxxxzz_yz_s_0_0_0[i] = 3.0 / 2.0 * tg_xxzz_yz_s_0_0_0[i] * fzi_0 + 3.0 * tg_xxzz_yz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xxxzz_yz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxxzz_yz_s_0_0_0[i] * a_x * faz_0;

        tg_xxxxzz_zz_s_0_0_0[i] = 3.0 / 2.0 * tg_xxzz_zz_s_0_0_0[i] * fzi_0 + 3.0 * tg_xxzz_zz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xxxzz_zz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxxzz_zz_s_0_0_0[i] * a_x * faz_0;

        tg_xxxyyy_xx_s_0_0_0[i] = tg_xyyy_xx_s_0_0_0[i] * fzi_0 + 2.0 * tg_xyyy_xx_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xxyyy_xx_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxyyy_xx_s_0_0_0[i] * a_x * faz_0;

        tg_xxxyyy_xy_s_0_0_0[i] = tg_xyyy_xy_s_0_0_0[i] * fzi_0 + 2.0 * tg_xyyy_xy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xxyyy_xy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxyyy_xy_s_0_0_0[i] * a_x * faz_0;

        tg_xxxyyy_xz_s_0_0_0[i] = tg_xyyy_xz_s_0_0_0[i] * fzi_0 + 2.0 * tg_xyyy_xz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xxyyy_xz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxyyy_xz_s_0_0_0[i] * a_x * faz_0;

        tg_xxxyyy_yy_s_0_0_0[i] = tg_xyyy_yy_s_0_0_0[i] * fzi_0 + 2.0 * tg_xyyy_yy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xxyyy_yy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxyyy_yy_s_0_0_0[i] * a_x * faz_0;

        tg_xxxyyy_yz_s_0_0_0[i] = tg_xyyy_yz_s_0_0_0[i] * fzi_0 + 2.0 * tg_xyyy_yz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xxyyy_yz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxyyy_yz_s_0_0_0[i] * a_x * faz_0;

        tg_xxxyyy_zz_s_0_0_0[i] = tg_xyyy_zz_s_0_0_0[i] * fzi_0 + 2.0 * tg_xyyy_zz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xxyyy_zz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxyyy_zz_s_0_0_0[i] * a_x * faz_0;

        tg_xxxyyz_xx_s_0_0_0[i] = 2.0 * tg_xxxyy_xx_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxxyy_xx_s_0_0_0[i] * a_z * faz_0;

        tg_xxxyyz_xy_s_0_0_0[i] = 2.0 * tg_xxxyy_xy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxxyy_xy_s_0_0_0[i] * a_z * faz_0;

        tg_xxxyyz_xz_s_0_0_0[i] = 2.0 * tg_xxxyy_xz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxxyy_xz_s_0_0_0[i] * a_z * faz_0;

        tg_xxxyyz_yy_s_0_0_0[i] = 2.0 * tg_xxxyy_yy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxxyy_yy_s_0_0_0[i] * a_z * faz_0;

        tg_xxxyyz_yz_s_0_0_0[i] = 2.0 * tg_xxxyy_yz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxxyy_yz_s_0_0_0[i] * a_z * faz_0;

        tg_xxxyyz_zz_s_0_0_0[i] = 2.0 * tg_xxxyy_zz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxxyy_zz_s_0_0_0[i] * a_z * faz_0;

        tg_xxxyzz_xx_s_0_0_0[i] = 2.0 * tg_xxxzz_xx_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxzz_xx_s_0_0_0[i] * a_y * faz_0;

        tg_xxxyzz_xy_s_0_0_0[i] = 2.0 * tg_xxxzz_xy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxzz_xy_s_0_0_0[i] * a_y * faz_0;

        tg_xxxyzz_xz_s_0_0_0[i] = 2.0 * tg_xxxzz_xz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxzz_xz_s_0_0_0[i] * a_y * faz_0;

        tg_xxxyzz_yy_s_0_0_0[i] = 2.0 * tg_xxxzz_yy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxzz_yy_s_0_0_0[i] * a_y * faz_0;

        tg_xxxyzz_yz_s_0_0_0[i] = 2.0 * tg_xxxzz_yz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxzz_yz_s_0_0_0[i] * a_y * faz_0;

        tg_xxxyzz_zz_s_0_0_0[i] = 2.0 * tg_xxxzz_zz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxzz_zz_s_0_0_0[i] * a_y * faz_0;

        tg_xxxzzz_xx_s_0_0_0[i] = tg_xzzz_xx_s_0_0_0[i] * fzi_0 + 2.0 * tg_xzzz_xx_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xxzzz_xx_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxzzz_xx_s_0_0_0[i] * a_x * faz_0;

        tg_xxxzzz_xy_s_0_0_0[i] = tg_xzzz_xy_s_0_0_0[i] * fzi_0 + 2.0 * tg_xzzz_xy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xxzzz_xy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxzzz_xy_s_0_0_0[i] * a_x * faz_0;

        tg_xxxzzz_xz_s_0_0_0[i] = tg_xzzz_xz_s_0_0_0[i] * fzi_0 + 2.0 * tg_xzzz_xz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xxzzz_xz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxzzz_xz_s_0_0_0[i] * a_x * faz_0;

        tg_xxxzzz_yy_s_0_0_0[i] = tg_xzzz_yy_s_0_0_0[i] * fzi_0 + 2.0 * tg_xzzz_yy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xxzzz_yy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxzzz_yy_s_0_0_0[i] * a_x * faz_0;

        tg_xxxzzz_yz_s_0_0_0[i] = tg_xzzz_yz_s_0_0_0[i] * fzi_0 + 2.0 * tg_xzzz_yz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xxzzz_yz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxzzz_yz_s_0_0_0[i] * a_x * faz_0;

        tg_xxxzzz_zz_s_0_0_0[i] = tg_xzzz_zz_s_0_0_0[i] * fzi_0 + 2.0 * tg_xzzz_zz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xxzzz_zz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxzzz_zz_s_0_0_0[i] * a_x * faz_0;

        tg_xxyyyy_xx_s_0_0_0[i] = 1.0 / 2.0 * tg_yyyy_xx_s_0_0_0[i] * fzi_0 + tg_yyyy_xx_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xyyyy_xx_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyyyy_xx_s_0_0_0[i] * a_x * faz_0;

        tg_xxyyyy_xy_s_0_0_0[i] = 1.0 / 2.0 * tg_yyyy_xy_s_0_0_0[i] * fzi_0 + tg_yyyy_xy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xyyyy_xy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyyyy_xy_s_0_0_0[i] * a_x * faz_0;

        tg_xxyyyy_xz_s_0_0_0[i] = 1.0 / 2.0 * tg_yyyy_xz_s_0_0_0[i] * fzi_0 + tg_yyyy_xz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xyyyy_xz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyyyy_xz_s_0_0_0[i] * a_x * faz_0;

        tg_xxyyyy_yy_s_0_0_0[i] = 1.0 / 2.0 * tg_yyyy_yy_s_0_0_0[i] * fzi_0 + tg_yyyy_yy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xyyyy_yy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyyyy_yy_s_0_0_0[i] * a_x * faz_0;

        tg_xxyyyy_yz_s_0_0_0[i] = 1.0 / 2.0 * tg_yyyy_yz_s_0_0_0[i] * fzi_0 + tg_yyyy_yz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xyyyy_yz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyyyy_yz_s_0_0_0[i] * a_x * faz_0;

        tg_xxyyyy_zz_s_0_0_0[i] = 1.0 / 2.0 * tg_yyyy_zz_s_0_0_0[i] * fzi_0 + tg_yyyy_zz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xyyyy_zz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyyyy_zz_s_0_0_0[i] * a_x * faz_0;

        tg_xxyyyz_xx_s_0_0_0[i] = 2.0 * tg_xxyyy_xx_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxyyy_xx_s_0_0_0[i] * a_z * faz_0;

        tg_xxyyyz_xy_s_0_0_0[i] = 2.0 * tg_xxyyy_xy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxyyy_xy_s_0_0_0[i] * a_z * faz_0;

        tg_xxyyyz_xz_s_0_0_0[i] = 2.0 * tg_xxyyy_xz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxyyy_xz_s_0_0_0[i] * a_z * faz_0;

        tg_xxyyyz_yy_s_0_0_0[i] = 2.0 * tg_xxyyy_yy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxyyy_yy_s_0_0_0[i] * a_z * faz_0;

        tg_xxyyyz_yz_s_0_0_0[i] = 2.0 * tg_xxyyy_yz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxyyy_yz_s_0_0_0[i] * a_z * faz_0;

        tg_xxyyyz_zz_s_0_0_0[i] = 2.0 * tg_xxyyy_zz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxyyy_zz_s_0_0_0[i] * a_z * faz_0;

        tg_xxyyzz_xx_s_0_0_0[i] = 1.0 / 2.0 * tg_yyzz_xx_s_0_0_0[i] * fzi_0 + tg_yyzz_xx_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xyyzz_xx_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyyzz_xx_s_0_0_0[i] * a_x * faz_0;

        tg_xxyyzz_xy_s_0_0_0[i] = 1.0 / 2.0 * tg_yyzz_xy_s_0_0_0[i] * fzi_0 + tg_yyzz_xy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xyyzz_xy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyyzz_xy_s_0_0_0[i] * a_x * faz_0;

        tg_xxyyzz_xz_s_0_0_0[i] = 1.0 / 2.0 * tg_yyzz_xz_s_0_0_0[i] * fzi_0 + tg_yyzz_xz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xyyzz_xz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyyzz_xz_s_0_0_0[i] * a_x * faz_0;

        tg_xxyyzz_yy_s_0_0_0[i] = 1.0 / 2.0 * tg_yyzz_yy_s_0_0_0[i] * fzi_0 + tg_yyzz_yy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xyyzz_yy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyyzz_yy_s_0_0_0[i] * a_x * faz_0;

        tg_xxyyzz_yz_s_0_0_0[i] = 1.0 / 2.0 * tg_yyzz_yz_s_0_0_0[i] * fzi_0 + tg_yyzz_yz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xyyzz_yz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyyzz_yz_s_0_0_0[i] * a_x * faz_0;

        tg_xxyyzz_zz_s_0_0_0[i] = 1.0 / 2.0 * tg_yyzz_zz_s_0_0_0[i] * fzi_0 + tg_yyzz_zz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xyyzz_zz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyyzz_zz_s_0_0_0[i] * a_x * faz_0;

        tg_xxyzzz_xx_s_0_0_0[i] = 2.0 * tg_xxzzz_xx_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxzzz_xx_s_0_0_0[i] * a_y * faz_0;

        tg_xxyzzz_xy_s_0_0_0[i] = 2.0 * tg_xxzzz_xy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxzzz_xy_s_0_0_0[i] * a_y * faz_0;

        tg_xxyzzz_xz_s_0_0_0[i] = 2.0 * tg_xxzzz_xz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxzzz_xz_s_0_0_0[i] * a_y * faz_0;

        tg_xxyzzz_yy_s_0_0_0[i] = 2.0 * tg_xxzzz_yy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxzzz_yy_s_0_0_0[i] * a_y * faz_0;

        tg_xxyzzz_yz_s_0_0_0[i] = 2.0 * tg_xxzzz_yz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxzzz_yz_s_0_0_0[i] * a_y * faz_0;

        tg_xxyzzz_zz_s_0_0_0[i] = 2.0 * tg_xxzzz_zz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxzzz_zz_s_0_0_0[i] * a_y * faz_0;

        tg_xxzzzz_xx_s_0_0_0[i] = 1.0 / 2.0 * tg_zzzz_xx_s_0_0_0[i] * fzi_0 + tg_zzzz_xx_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xzzzz_xx_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xzzzz_xx_s_0_0_0[i] * a_x * faz_0;

        tg_xxzzzz_xy_s_0_0_0[i] = 1.0 / 2.0 * tg_zzzz_xy_s_0_0_0[i] * fzi_0 + tg_zzzz_xy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xzzzz_xy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xzzzz_xy_s_0_0_0[i] * a_x * faz_0;

        tg_xxzzzz_xz_s_0_0_0[i] = 1.0 / 2.0 * tg_zzzz_xz_s_0_0_0[i] * fzi_0 + tg_zzzz_xz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xzzzz_xz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xzzzz_xz_s_0_0_0[i] * a_x * faz_0;

        tg_xxzzzz_yy_s_0_0_0[i] = 1.0 / 2.0 * tg_zzzz_yy_s_0_0_0[i] * fzi_0 + tg_zzzz_yy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xzzzz_yy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xzzzz_yy_s_0_0_0[i] * a_x * faz_0;

        tg_xxzzzz_yz_s_0_0_0[i] = 1.0 / 2.0 * tg_zzzz_yz_s_0_0_0[i] * fzi_0 + tg_zzzz_yz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xzzzz_yz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xzzzz_yz_s_0_0_0[i] * a_x * faz_0;

        tg_xxzzzz_zz_s_0_0_0[i] = 1.0 / 2.0 * tg_zzzz_zz_s_0_0_0[i] * fzi_0 + tg_zzzz_zz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xzzzz_zz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xzzzz_zz_s_0_0_0[i] * a_x * faz_0;

        tg_xyyyyy_xx_s_0_0_0[i] = 2.0 * tg_yyyyy_xx_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyyy_xx_s_0_0_0[i] * a_x * faz_0;

        tg_xyyyyy_xy_s_0_0_0[i] = 2.0 * tg_yyyyy_xy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyyy_xy_s_0_0_0[i] * a_x * faz_0;

        tg_xyyyyy_xz_s_0_0_0[i] = 2.0 * tg_yyyyy_xz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyyy_xz_s_0_0_0[i] * a_x * faz_0;

        tg_xyyyyy_yy_s_0_0_0[i] = 2.0 * tg_yyyyy_yy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyyy_yy_s_0_0_0[i] * a_x * faz_0;

        tg_xyyyyy_yz_s_0_0_0[i] = 2.0 * tg_yyyyy_yz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyyy_yz_s_0_0_0[i] * a_x * faz_0;

        tg_xyyyyy_zz_s_0_0_0[i] = 2.0 * tg_yyyyy_zz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyyy_zz_s_0_0_0[i] * a_x * faz_0;

        tg_xyyyyz_xx_s_0_0_0[i] = 2.0 * tg_yyyyz_xx_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyyz_xx_s_0_0_0[i] * a_x * faz_0;

        tg_xyyyyz_xy_s_0_0_0[i] = 2.0 * tg_yyyyz_xy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyyz_xy_s_0_0_0[i] * a_x * faz_0;

        tg_xyyyyz_xz_s_0_0_0[i] = 2.0 * tg_yyyyz_xz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyyz_xz_s_0_0_0[i] * a_x * faz_0;

        tg_xyyyyz_yy_s_0_0_0[i] = 2.0 * tg_yyyyz_yy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyyz_yy_s_0_0_0[i] * a_x * faz_0;

        tg_xyyyyz_yz_s_0_0_0[i] = 2.0 * tg_yyyyz_yz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyyz_yz_s_0_0_0[i] * a_x * faz_0;

        tg_xyyyyz_zz_s_0_0_0[i] = 2.0 * tg_yyyyz_zz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyyz_zz_s_0_0_0[i] * a_x * faz_0;

        tg_xyyyzz_xx_s_0_0_0[i] = 2.0 * tg_yyyzz_xx_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyzz_xx_s_0_0_0[i] * a_x * faz_0;

        tg_xyyyzz_xy_s_0_0_0[i] = 2.0 * tg_yyyzz_xy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyzz_xy_s_0_0_0[i] * a_x * faz_0;

        tg_xyyyzz_xz_s_0_0_0[i] = 2.0 * tg_yyyzz_xz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyzz_xz_s_0_0_0[i] * a_x * faz_0;

        tg_xyyyzz_yy_s_0_0_0[i] = 2.0 * tg_yyyzz_yy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyzz_yy_s_0_0_0[i] * a_x * faz_0;

        tg_xyyyzz_yz_s_0_0_0[i] = 2.0 * tg_yyyzz_yz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyzz_yz_s_0_0_0[i] * a_x * faz_0;

        tg_xyyyzz_zz_s_0_0_0[i] = 2.0 * tg_yyyzz_zz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyzz_zz_s_0_0_0[i] * a_x * faz_0;

        tg_xyyzzz_xx_s_0_0_0[i] = 2.0 * tg_yyzzz_xx_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyzzz_xx_s_0_0_0[i] * a_x * faz_0;

        tg_xyyzzz_xy_s_0_0_0[i] = 2.0 * tg_yyzzz_xy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyzzz_xy_s_0_0_0[i] * a_x * faz_0;

        tg_xyyzzz_xz_s_0_0_0[i] = 2.0 * tg_yyzzz_xz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyzzz_xz_s_0_0_0[i] * a_x * faz_0;

        tg_xyyzzz_yy_s_0_0_0[i] = 2.0 * tg_yyzzz_yy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyzzz_yy_s_0_0_0[i] * a_x * faz_0;

        tg_xyyzzz_yz_s_0_0_0[i] = 2.0 * tg_yyzzz_yz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyzzz_yz_s_0_0_0[i] * a_x * faz_0;

        tg_xyyzzz_zz_s_0_0_0[i] = 2.0 * tg_yyzzz_zz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyzzz_zz_s_0_0_0[i] * a_x * faz_0;

        tg_xyzzzz_xx_s_0_0_0[i] = 2.0 * tg_yzzzz_xx_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yzzzz_xx_s_0_0_0[i] * a_x * faz_0;

        tg_xyzzzz_xy_s_0_0_0[i] = 2.0 * tg_yzzzz_xy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yzzzz_xy_s_0_0_0[i] * a_x * faz_0;

        tg_xyzzzz_xz_s_0_0_0[i] = 2.0 * tg_yzzzz_xz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yzzzz_xz_s_0_0_0[i] * a_x * faz_0;

        tg_xyzzzz_yy_s_0_0_0[i] = 2.0 * tg_yzzzz_yy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yzzzz_yy_s_0_0_0[i] * a_x * faz_0;

        tg_xyzzzz_yz_s_0_0_0[i] = 2.0 * tg_yzzzz_yz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yzzzz_yz_s_0_0_0[i] * a_x * faz_0;

        tg_xyzzzz_zz_s_0_0_0[i] = 2.0 * tg_yzzzz_zz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yzzzz_zz_s_0_0_0[i] * a_x * faz_0;

        tg_xzzzzz_xx_s_0_0_0[i] = 2.0 * tg_zzzzz_xx_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzzzz_xx_s_0_0_0[i] * a_x * faz_0;

        tg_xzzzzz_xy_s_0_0_0[i] = 2.0 * tg_zzzzz_xy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzzzz_xy_s_0_0_0[i] * a_x * faz_0;

        tg_xzzzzz_xz_s_0_0_0[i] = 2.0 * tg_zzzzz_xz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzzzz_xz_s_0_0_0[i] * a_x * faz_0;

        tg_xzzzzz_yy_s_0_0_0[i] = 2.0 * tg_zzzzz_yy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzzzz_yy_s_0_0_0[i] * a_x * faz_0;

        tg_xzzzzz_yz_s_0_0_0[i] = 2.0 * tg_zzzzz_yz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzzzz_yz_s_0_0_0[i] * a_x * faz_0;

        tg_xzzzzz_zz_s_0_0_0[i] = 2.0 * tg_zzzzz_zz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzzzz_zz_s_0_0_0[i] * a_x * faz_0;

        tg_yyyyyy_xx_s_0_0_0[i] = 5.0 / 2.0 * tg_yyyy_xx_s_0_0_0[i] * fzi_0 + 5.0 * tg_yyyy_xx_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yyyyy_xx_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyyyy_xx_s_0_0_0[i] * a_y * faz_0;

        tg_yyyyyy_xy_s_0_0_0[i] = 5.0 / 2.0 * tg_yyyy_xy_s_0_0_0[i] * fzi_0 + 5.0 * tg_yyyy_xy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yyyyy_xy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyyyy_xy_s_0_0_0[i] * a_y * faz_0;

        tg_yyyyyy_xz_s_0_0_0[i] = 5.0 / 2.0 * tg_yyyy_xz_s_0_0_0[i] * fzi_0 + 5.0 * tg_yyyy_xz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yyyyy_xz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyyyy_xz_s_0_0_0[i] * a_y * faz_0;

        tg_yyyyyy_yy_s_0_0_0[i] = 5.0 / 2.0 * tg_yyyy_yy_s_0_0_0[i] * fzi_0 + 5.0 * tg_yyyy_yy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yyyyy_yy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyyyy_yy_s_0_0_0[i] * a_y * faz_0;

        tg_yyyyyy_yz_s_0_0_0[i] = 5.0 / 2.0 * tg_yyyy_yz_s_0_0_0[i] * fzi_0 + 5.0 * tg_yyyy_yz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yyyyy_yz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyyyy_yz_s_0_0_0[i] * a_y * faz_0;

        tg_yyyyyy_zz_s_0_0_0[i] = 5.0 / 2.0 * tg_yyyy_zz_s_0_0_0[i] * fzi_0 + 5.0 * tg_yyyy_zz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yyyyy_zz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyyyy_zz_s_0_0_0[i] * a_y * faz_0;

        tg_yyyyyz_xx_s_0_0_0[i] = 2.0 * tg_yyyyy_xx_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyyyy_xx_s_0_0_0[i] * a_z * faz_0;

        tg_yyyyyz_xy_s_0_0_0[i] = 2.0 * tg_yyyyy_xy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyyyy_xy_s_0_0_0[i] * a_z * faz_0;

        tg_yyyyyz_xz_s_0_0_0[i] = 2.0 * tg_yyyyy_xz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyyyy_xz_s_0_0_0[i] * a_z * faz_0;

        tg_yyyyyz_yy_s_0_0_0[i] = 2.0 * tg_yyyyy_yy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyyyy_yy_s_0_0_0[i] * a_z * faz_0;

        tg_yyyyyz_yz_s_0_0_0[i] = 2.0 * tg_yyyyy_yz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyyyy_yz_s_0_0_0[i] * a_z * faz_0;

        tg_yyyyyz_zz_s_0_0_0[i] = 2.0 * tg_yyyyy_zz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyyyy_zz_s_0_0_0[i] * a_z * faz_0;

        tg_yyyyzz_xx_s_0_0_0[i] = 3.0 / 2.0 * tg_yyzz_xx_s_0_0_0[i] * fzi_0 + 3.0 * tg_yyzz_xx_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yyyzz_xx_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyyzz_xx_s_0_0_0[i] * a_y * faz_0;

        tg_yyyyzz_xy_s_0_0_0[i] = 3.0 / 2.0 * tg_yyzz_xy_s_0_0_0[i] * fzi_0 + 3.0 * tg_yyzz_xy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yyyzz_xy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyyzz_xy_s_0_0_0[i] * a_y * faz_0;

        tg_yyyyzz_xz_s_0_0_0[i] = 3.0 / 2.0 * tg_yyzz_xz_s_0_0_0[i] * fzi_0 + 3.0 * tg_yyzz_xz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yyyzz_xz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyyzz_xz_s_0_0_0[i] * a_y * faz_0;

        tg_yyyyzz_yy_s_0_0_0[i] = 3.0 / 2.0 * tg_yyzz_yy_s_0_0_0[i] * fzi_0 + 3.0 * tg_yyzz_yy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yyyzz_yy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyyzz_yy_s_0_0_0[i] * a_y * faz_0;

        tg_yyyyzz_yz_s_0_0_0[i] = 3.0 / 2.0 * tg_yyzz_yz_s_0_0_0[i] * fzi_0 + 3.0 * tg_yyzz_yz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yyyzz_yz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyyzz_yz_s_0_0_0[i] * a_y * faz_0;

        tg_yyyyzz_zz_s_0_0_0[i] = 3.0 / 2.0 * tg_yyzz_zz_s_0_0_0[i] * fzi_0 + 3.0 * tg_yyzz_zz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yyyzz_zz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyyzz_zz_s_0_0_0[i] * a_y * faz_0;

        tg_yyyzzz_xx_s_0_0_0[i] = tg_yzzz_xx_s_0_0_0[i] * fzi_0 + 2.0 * tg_yzzz_xx_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yyzzz_xx_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyzzz_xx_s_0_0_0[i] * a_y * faz_0;

        tg_yyyzzz_xy_s_0_0_0[i] = tg_yzzz_xy_s_0_0_0[i] * fzi_0 + 2.0 * tg_yzzz_xy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yyzzz_xy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyzzz_xy_s_0_0_0[i] * a_y * faz_0;

        tg_yyyzzz_xz_s_0_0_0[i] = tg_yzzz_xz_s_0_0_0[i] * fzi_0 + 2.0 * tg_yzzz_xz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yyzzz_xz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyzzz_xz_s_0_0_0[i] * a_y * faz_0;

        tg_yyyzzz_yy_s_0_0_0[i] = tg_yzzz_yy_s_0_0_0[i] * fzi_0 + 2.0 * tg_yzzz_yy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yyzzz_yy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyzzz_yy_s_0_0_0[i] * a_y * faz_0;

        tg_yyyzzz_yz_s_0_0_0[i] = tg_yzzz_yz_s_0_0_0[i] * fzi_0 + 2.0 * tg_yzzz_yz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yyzzz_yz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyzzz_yz_s_0_0_0[i] * a_y * faz_0;

        tg_yyyzzz_zz_s_0_0_0[i] = tg_yzzz_zz_s_0_0_0[i] * fzi_0 + 2.0 * tg_yzzz_zz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yyzzz_zz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyzzz_zz_s_0_0_0[i] * a_y * faz_0;

        tg_yyzzzz_xx_s_0_0_0[i] = 1.0 / 2.0 * tg_zzzz_xx_s_0_0_0[i] * fzi_0 + tg_zzzz_xx_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yzzzz_xx_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yzzzz_xx_s_0_0_0[i] * a_y * faz_0;

        tg_yyzzzz_xy_s_0_0_0[i] = 1.0 / 2.0 * tg_zzzz_xy_s_0_0_0[i] * fzi_0 + tg_zzzz_xy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yzzzz_xy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yzzzz_xy_s_0_0_0[i] * a_y * faz_0;

        tg_yyzzzz_xz_s_0_0_0[i] = 1.0 / 2.0 * tg_zzzz_xz_s_0_0_0[i] * fzi_0 + tg_zzzz_xz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yzzzz_xz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yzzzz_xz_s_0_0_0[i] * a_y * faz_0;

        tg_yyzzzz_yy_s_0_0_0[i] = 1.0 / 2.0 * tg_zzzz_yy_s_0_0_0[i] * fzi_0 + tg_zzzz_yy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yzzzz_yy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yzzzz_yy_s_0_0_0[i] * a_y * faz_0;

        tg_yyzzzz_yz_s_0_0_0[i] = 1.0 / 2.0 * tg_zzzz_yz_s_0_0_0[i] * fzi_0 + tg_zzzz_yz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yzzzz_yz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yzzzz_yz_s_0_0_0[i] * a_y * faz_0;

        tg_yyzzzz_zz_s_0_0_0[i] = 1.0 / 2.0 * tg_zzzz_zz_s_0_0_0[i] * fzi_0 + tg_zzzz_zz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yzzzz_zz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yzzzz_zz_s_0_0_0[i] * a_y * faz_0;

        tg_yzzzzz_xx_s_0_0_0[i] = 2.0 * tg_zzzzz_xx_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzzzz_xx_s_0_0_0[i] * a_y * faz_0;

        tg_yzzzzz_xy_s_0_0_0[i] = 2.0 * tg_zzzzz_xy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzzzz_xy_s_0_0_0[i] * a_y * faz_0;

        tg_yzzzzz_xz_s_0_0_0[i] = 2.0 * tg_zzzzz_xz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzzzz_xz_s_0_0_0[i] * a_y * faz_0;

        tg_yzzzzz_yy_s_0_0_0[i] = 2.0 * tg_zzzzz_yy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzzzz_yy_s_0_0_0[i] * a_y * faz_0;

        tg_yzzzzz_yz_s_0_0_0[i] = 2.0 * tg_zzzzz_yz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzzzz_yz_s_0_0_0[i] * a_y * faz_0;

        tg_yzzzzz_zz_s_0_0_0[i] = 2.0 * tg_zzzzz_zz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzzzz_zz_s_0_0_0[i] * a_y * faz_0;

        tg_zzzzzz_xx_s_0_0_0[i] = 5.0 / 2.0 * tg_zzzz_xx_s_0_0_0[i] * fzi_0 + 5.0 * tg_zzzz_xx_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_zzzzz_xx_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzzzz_xx_s_0_0_0[i] * a_z * faz_0;

        tg_zzzzzz_xy_s_0_0_0[i] = 5.0 / 2.0 * tg_zzzz_xy_s_0_0_0[i] * fzi_0 + 5.0 * tg_zzzz_xy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_zzzzz_xy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzzzz_xy_s_0_0_0[i] * a_z * faz_0;

        tg_zzzzzz_xz_s_0_0_0[i] = 5.0 / 2.0 * tg_zzzz_xz_s_0_0_0[i] * fzi_0 + 5.0 * tg_zzzz_xz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_zzzzz_xz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzzzz_xz_s_0_0_0[i] * a_z * faz_0;

        tg_zzzzzz_yy_s_0_0_0[i] = 5.0 / 2.0 * tg_zzzz_yy_s_0_0_0[i] * fzi_0 + 5.0 * tg_zzzz_yy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_zzzzz_yy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzzzz_yy_s_0_0_0[i] * a_z * faz_0;

        tg_zzzzzz_yz_s_0_0_0[i] = 5.0 / 2.0 * tg_zzzz_yz_s_0_0_0[i] * fzi_0 + 5.0 * tg_zzzz_yz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_zzzzz_yz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzzzz_yz_s_0_0_0[i] * a_z * faz_0;

        tg_zzzzzz_zz_s_0_0_0[i] = 5.0 / 2.0 * tg_zzzz_zz_s_0_0_0[i] * fzi_0 + 5.0 * tg_zzzz_zz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_zzzzz_zz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzzzz_zz_s_0_0_0[i] * a_z * faz_0;
    }

    if (p > 0)
    {
        const double fp_0 = (double)p;

        // Set up components of auxiliary buffer : GD

        auto tg_xxxx_xx_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1);

        auto tg_xxxx_xy_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 1);

        auto tg_xxxx_xz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 2);

        auto tg_xxxx_yy_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 3);

        auto tg_xxxx_yz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 4);

        auto tg_xxxx_zz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 5);

        auto tg_xxxy_xx_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 6);

        auto tg_xxxy_xy_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 7);

        auto tg_xxxy_xz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 8);

        auto tg_xxxy_yy_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 9);

        auto tg_xxxy_yz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 10);

        auto tg_xxxy_zz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 11);

        auto tg_xxxz_xx_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 12);

        auto tg_xxxz_xy_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 13);

        auto tg_xxxz_xz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 14);

        auto tg_xxxz_yy_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 15);

        auto tg_xxxz_yz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 16);

        auto tg_xxxz_zz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 17);

        auto tg_xxyy_xx_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 18);

        auto tg_xxyy_xy_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 19);

        auto tg_xxyy_xz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 20);

        auto tg_xxyy_yy_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 21);

        auto tg_xxyy_yz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 22);

        auto tg_xxyy_zz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 23);

        auto tg_xxyz_xx_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 24);

        auto tg_xxyz_xy_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 25);

        auto tg_xxyz_xz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 26);

        auto tg_xxyz_yy_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 27);

        auto tg_xxyz_yz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 28);

        auto tg_xxyz_zz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 29);

        auto tg_xxzz_xx_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 30);

        auto tg_xxzz_xy_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 31);

        auto tg_xxzz_xz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 32);

        auto tg_xxzz_yy_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 33);

        auto tg_xxzz_yz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 34);

        auto tg_xxzz_zz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 35);

        auto tg_xyyy_xx_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 36);

        auto tg_xyyy_xy_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 37);

        auto tg_xyyy_xz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 38);

        auto tg_xyyy_yy_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 39);

        auto tg_xyyy_yz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 40);

        auto tg_xyyy_zz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 41);

        auto tg_xyyz_xx_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 42);

        auto tg_xyyz_xy_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 43);

        auto tg_xyyz_xz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 44);

        auto tg_xyyz_yy_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 45);

        auto tg_xyyz_yz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 46);

        auto tg_xyyz_zz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 47);

        auto tg_xyzz_xx_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 48);

        auto tg_xyzz_xy_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 49);

        auto tg_xyzz_xz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 50);

        auto tg_xyzz_yy_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 51);

        auto tg_xyzz_yz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 52);

        auto tg_xyzz_zz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 53);

        auto tg_xzzz_xx_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 54);

        auto tg_xzzz_xy_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 55);

        auto tg_xzzz_xz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 56);

        auto tg_xzzz_yy_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 57);

        auto tg_xzzz_yz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 58);

        auto tg_xzzz_zz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 59);

        auto tg_yyyy_xx_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 60);

        auto tg_yyyy_xy_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 61);

        auto tg_yyyy_xz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 62);

        auto tg_yyyy_yy_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 63);

        auto tg_yyyy_yz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 64);

        auto tg_yyyy_zz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 65);

        auto tg_yyyz_xx_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 66);

        auto tg_yyyz_xy_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 67);

        auto tg_yyyz_xz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 68);

        auto tg_yyyz_yy_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 69);

        auto tg_yyyz_yz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 70);

        auto tg_yyyz_zz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 71);

        auto tg_yyzz_xx_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 72);

        auto tg_yyzz_xy_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 73);

        auto tg_yyzz_xz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 74);

        auto tg_yyzz_yy_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 75);

        auto tg_yyzz_yz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 76);

        auto tg_yyzz_zz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 77);

        auto tg_yzzz_xx_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 78);

        auto tg_yzzz_xy_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 79);

        auto tg_yzzz_xz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 80);

        auto tg_yzzz_yy_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 81);

        auto tg_yzzz_yz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 82);

        auto tg_yzzz_zz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 83);

        auto tg_zzzz_xx_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 84);

        auto tg_zzzz_xy_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 85);

        auto tg_zzzz_xz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 86);

        auto tg_zzzz_yy_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 87);

        auto tg_zzzz_yz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 88);

        auto tg_zzzz_zz_s_0_0_1 = pbuffer.data(idx_gd_s_0_0_1 + 89);

        // Set up components of auxiliary buffer : HD

        auto tg_xxxxx_xx_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1);

        auto tg_xxxxx_xy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 1);

        auto tg_xxxxx_xz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 2);

        auto tg_xxxxx_yy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 3);

        auto tg_xxxxx_yz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 4);

        auto tg_xxxxx_zz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 5);

        auto tg_xxxxy_xx_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 6);

        auto tg_xxxxy_xy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 7);

        auto tg_xxxxy_xz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 8);

        auto tg_xxxxy_yy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 9);

        auto tg_xxxxy_yz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 10);

        auto tg_xxxxy_zz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 11);

        auto tg_xxxxz_xx_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 12);

        auto tg_xxxxz_xy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 13);

        auto tg_xxxxz_xz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 14);

        auto tg_xxxxz_yy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 15);

        auto tg_xxxxz_yz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 16);

        auto tg_xxxxz_zz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 17);

        auto tg_xxxyy_xx_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 18);

        auto tg_xxxyy_xy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 19);

        auto tg_xxxyy_xz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 20);

        auto tg_xxxyy_yy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 21);

        auto tg_xxxyy_yz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 22);

        auto tg_xxxyy_zz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 23);

        auto tg_xxxyz_xx_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 24);

        auto tg_xxxyz_xy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 25);

        auto tg_xxxyz_xz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 26);

        auto tg_xxxyz_yy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 27);

        auto tg_xxxyz_yz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 28);

        auto tg_xxxyz_zz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 29);

        auto tg_xxxzz_xx_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 30);

        auto tg_xxxzz_xy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 31);

        auto tg_xxxzz_xz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 32);

        auto tg_xxxzz_yy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 33);

        auto tg_xxxzz_yz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 34);

        auto tg_xxxzz_zz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 35);

        auto tg_xxyyy_xx_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 36);

        auto tg_xxyyy_xy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 37);

        auto tg_xxyyy_xz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 38);

        auto tg_xxyyy_yy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 39);

        auto tg_xxyyy_yz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 40);

        auto tg_xxyyy_zz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 41);

        auto tg_xxyyz_xx_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 42);

        auto tg_xxyyz_xy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 43);

        auto tg_xxyyz_xz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 44);

        auto tg_xxyyz_yy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 45);

        auto tg_xxyyz_yz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 46);

        auto tg_xxyyz_zz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 47);

        auto tg_xxyzz_xx_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 48);

        auto tg_xxyzz_xy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 49);

        auto tg_xxyzz_xz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 50);

        auto tg_xxyzz_yy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 51);

        auto tg_xxyzz_yz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 52);

        auto tg_xxyzz_zz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 53);

        auto tg_xxzzz_xx_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 54);

        auto tg_xxzzz_xy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 55);

        auto tg_xxzzz_xz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 56);

        auto tg_xxzzz_yy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 57);

        auto tg_xxzzz_yz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 58);

        auto tg_xxzzz_zz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 59);

        auto tg_xyyyy_xx_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 60);

        auto tg_xyyyy_xy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 61);

        auto tg_xyyyy_xz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 62);

        auto tg_xyyyy_yy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 63);

        auto tg_xyyyy_yz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 64);

        auto tg_xyyyy_zz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 65);

        auto tg_xyyyz_xx_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 66);

        auto tg_xyyyz_xy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 67);

        auto tg_xyyyz_xz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 68);

        auto tg_xyyyz_yy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 69);

        auto tg_xyyyz_yz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 70);

        auto tg_xyyyz_zz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 71);

        auto tg_xyyzz_xx_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 72);

        auto tg_xyyzz_xy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 73);

        auto tg_xyyzz_xz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 74);

        auto tg_xyyzz_yy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 75);

        auto tg_xyyzz_yz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 76);

        auto tg_xyyzz_zz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 77);

        auto tg_xyzzz_xx_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 78);

        auto tg_xyzzz_xy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 79);

        auto tg_xyzzz_xz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 80);

        auto tg_xyzzz_yy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 81);

        auto tg_xyzzz_yz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 82);

        auto tg_xyzzz_zz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 83);

        auto tg_xzzzz_xx_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 84);

        auto tg_xzzzz_xy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 85);

        auto tg_xzzzz_xz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 86);

        auto tg_xzzzz_yy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 87);

        auto tg_xzzzz_yz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 88);

        auto tg_xzzzz_zz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 89);

        auto tg_yyyyy_xx_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 90);

        auto tg_yyyyy_xy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 91);

        auto tg_yyyyy_xz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 92);

        auto tg_yyyyy_yy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 93);

        auto tg_yyyyy_yz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 94);

        auto tg_yyyyy_zz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 95);

        auto tg_yyyyz_xx_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 96);

        auto tg_yyyyz_xy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 97);

        auto tg_yyyyz_xz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 98);

        auto tg_yyyyz_yy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 99);

        auto tg_yyyyz_yz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 100);

        auto tg_yyyyz_zz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 101);

        auto tg_yyyzz_xx_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 102);

        auto tg_yyyzz_xy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 103);

        auto tg_yyyzz_xz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 104);

        auto tg_yyyzz_yy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 105);

        auto tg_yyyzz_yz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 106);

        auto tg_yyyzz_zz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 107);

        auto tg_yyzzz_xx_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 108);

        auto tg_yyzzz_xy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 109);

        auto tg_yyzzz_xz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 110);

        auto tg_yyzzz_yy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 111);

        auto tg_yyzzz_yz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 112);

        auto tg_yyzzz_zz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 113);

        auto tg_yzzzz_xx_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 114);

        auto tg_yzzzz_xy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 115);

        auto tg_yzzzz_xz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 116);

        auto tg_yzzzz_yy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 117);

        auto tg_yzzzz_yz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 118);

        auto tg_yzzzz_zz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 119);

        auto tg_zzzzz_xx_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 120);

        auto tg_zzzzz_xy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 121);

        auto tg_zzzzz_xz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 122);

        auto tg_zzzzz_yy_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 123);

        auto tg_zzzzz_yz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 124);

        auto tg_zzzzz_zz_s_0_0_1 = pbuffer.data(idx_hd_s_0_0_1 + 125);

        #pragma omp simd aligned(b_exps, tg_xxxx_xx_s_0_0_1, tg_xxxx_xy_s_0_0_1, tg_xxxx_xz_s_0_0_1, tg_xxxx_yy_s_0_0_1, tg_xxxx_yz_s_0_0_1, tg_xxxx_zz_s_0_0_1, tg_xxxxx_xx_s_0_0_1, tg_xxxxx_xy_s_0_0_1, tg_xxxxx_xz_s_0_0_1, tg_xxxxx_yy_s_0_0_1, tg_xxxxx_yz_s_0_0_1, tg_xxxxx_zz_s_0_0_1, tg_xxxxxx_xx_s_0_0_0, tg_xxxxxx_xy_s_0_0_0, tg_xxxxxx_xz_s_0_0_0, tg_xxxxxx_yy_s_0_0_0, tg_xxxxxx_yz_s_0_0_0, tg_xxxxxx_zz_s_0_0_0, tg_xxxxxy_xx_s_0_0_0, tg_xxxxxy_xy_s_0_0_0, tg_xxxxxy_xz_s_0_0_0, tg_xxxxxy_yy_s_0_0_0, tg_xxxxxy_yz_s_0_0_0, tg_xxxxxy_zz_s_0_0_0, tg_xxxxxz_xx_s_0_0_0, tg_xxxxxz_xy_s_0_0_0, tg_xxxxxz_xz_s_0_0_0, tg_xxxxxz_yy_s_0_0_0, tg_xxxxxz_yz_s_0_0_0, tg_xxxxxz_zz_s_0_0_0, tg_xxxxyy_xx_s_0_0_0, tg_xxxxyy_xy_s_0_0_0, tg_xxxxyy_xz_s_0_0_0, tg_xxxxyy_yy_s_0_0_0, tg_xxxxyy_yz_s_0_0_0, tg_xxxxyy_zz_s_0_0_0, tg_xxxxyz_xx_s_0_0_0, tg_xxxxyz_xy_s_0_0_0, tg_xxxxyz_xz_s_0_0_0, tg_xxxxyz_yy_s_0_0_0, tg_xxxxyz_yz_s_0_0_0, tg_xxxxyz_zz_s_0_0_0, tg_xxxxz_xx_s_0_0_1, tg_xxxxz_xy_s_0_0_1, tg_xxxxz_xz_s_0_0_1, tg_xxxxz_yy_s_0_0_1, tg_xxxxz_yz_s_0_0_1, tg_xxxxz_zz_s_0_0_1, tg_xxxxzz_xx_s_0_0_0, tg_xxxxzz_xy_s_0_0_0, tg_xxxxzz_xz_s_0_0_0, tg_xxxxzz_yy_s_0_0_0, tg_xxxxzz_yz_s_0_0_0, tg_xxxxzz_zz_s_0_0_0, tg_xxxyy_xx_s_0_0_1, tg_xxxyy_xy_s_0_0_1, tg_xxxyy_xz_s_0_0_1, tg_xxxyy_yy_s_0_0_1, tg_xxxyy_yz_s_0_0_1, tg_xxxyy_zz_s_0_0_1, tg_xxxyyy_xx_s_0_0_0, tg_xxxyyy_xy_s_0_0_0, tg_xxxyyy_xz_s_0_0_0, tg_xxxyyy_yy_s_0_0_0, tg_xxxyyy_yz_s_0_0_0, tg_xxxyyy_zz_s_0_0_0, tg_xxxyyz_xx_s_0_0_0, tg_xxxyyz_xy_s_0_0_0, tg_xxxyyz_xz_s_0_0_0, tg_xxxyyz_yy_s_0_0_0, tg_xxxyyz_yz_s_0_0_0, tg_xxxyyz_zz_s_0_0_0, tg_xxxyzz_xx_s_0_0_0, tg_xxxyzz_xy_s_0_0_0, tg_xxxyzz_xz_s_0_0_0, tg_xxxyzz_yy_s_0_0_0, tg_xxxyzz_yz_s_0_0_0, tg_xxxyzz_zz_s_0_0_0, tg_xxxzz_xx_s_0_0_1, tg_xxxzz_xy_s_0_0_1, tg_xxxzz_xz_s_0_0_1, tg_xxxzz_yy_s_0_0_1, tg_xxxzz_yz_s_0_0_1, tg_xxxzz_zz_s_0_0_1, tg_xxxzzz_xx_s_0_0_0, tg_xxxzzz_xy_s_0_0_0, tg_xxxzzz_xz_s_0_0_0, tg_xxxzzz_yy_s_0_0_0, tg_xxxzzz_yz_s_0_0_0, tg_xxxzzz_zz_s_0_0_0, tg_xxyy_xx_s_0_0_1, tg_xxyy_xy_s_0_0_1, tg_xxyy_xz_s_0_0_1, tg_xxyy_yy_s_0_0_1, tg_xxyy_yz_s_0_0_1, tg_xxyy_zz_s_0_0_1, tg_xxyyy_xx_s_0_0_1, tg_xxyyy_xy_s_0_0_1, tg_xxyyy_xz_s_0_0_1, tg_xxyyy_yy_s_0_0_1, tg_xxyyy_yz_s_0_0_1, tg_xxyyy_zz_s_0_0_1, tg_xxyyyy_xx_s_0_0_0, tg_xxyyyy_xy_s_0_0_0, tg_xxyyyy_xz_s_0_0_0, tg_xxyyyy_yy_s_0_0_0, tg_xxyyyy_yz_s_0_0_0, tg_xxyyyy_zz_s_0_0_0, tg_xxyyyz_xx_s_0_0_0, tg_xxyyyz_xy_s_0_0_0, tg_xxyyyz_xz_s_0_0_0, tg_xxyyyz_yy_s_0_0_0, tg_xxyyyz_yz_s_0_0_0, tg_xxyyyz_zz_s_0_0_0, tg_xxyyzz_xx_s_0_0_0, tg_xxyyzz_xy_s_0_0_0, tg_xxyyzz_xz_s_0_0_0, tg_xxyyzz_yy_s_0_0_0, tg_xxyyzz_yz_s_0_0_0, tg_xxyyzz_zz_s_0_0_0, tg_xxyzzz_xx_s_0_0_0, tg_xxyzzz_xy_s_0_0_0, tg_xxyzzz_xz_s_0_0_0, tg_xxyzzz_yy_s_0_0_0, tg_xxyzzz_yz_s_0_0_0, tg_xxyzzz_zz_s_0_0_0, tg_xxzz_xx_s_0_0_1, tg_xxzz_xy_s_0_0_1, tg_xxzz_xz_s_0_0_1, tg_xxzz_yy_s_0_0_1, tg_xxzz_yz_s_0_0_1, tg_xxzz_zz_s_0_0_1, tg_xxzzz_xx_s_0_0_1, tg_xxzzz_xy_s_0_0_1, tg_xxzzz_xz_s_0_0_1, tg_xxzzz_yy_s_0_0_1, tg_xxzzz_yz_s_0_0_1, tg_xxzzz_zz_s_0_0_1, tg_xxzzzz_xx_s_0_0_0, tg_xxzzzz_xy_s_0_0_0, tg_xxzzzz_xz_s_0_0_0, tg_xxzzzz_yy_s_0_0_0, tg_xxzzzz_yz_s_0_0_0, tg_xxzzzz_zz_s_0_0_0, tg_xyyy_xx_s_0_0_1, tg_xyyy_xy_s_0_0_1, tg_xyyy_xz_s_0_0_1, tg_xyyy_yy_s_0_0_1, tg_xyyy_yz_s_0_0_1, tg_xyyy_zz_s_0_0_1, tg_xyyyy_xx_s_0_0_1, tg_xyyyy_xy_s_0_0_1, tg_xyyyy_xz_s_0_0_1, tg_xyyyy_yy_s_0_0_1, tg_xyyyy_yz_s_0_0_1, tg_xyyyy_zz_s_0_0_1, tg_xyyyyy_xx_s_0_0_0, tg_xyyyyy_xy_s_0_0_0, tg_xyyyyy_xz_s_0_0_0, tg_xyyyyy_yy_s_0_0_0, tg_xyyyyy_yz_s_0_0_0, tg_xyyyyy_zz_s_0_0_0, tg_xyyyyz_xx_s_0_0_0, tg_xyyyyz_xy_s_0_0_0, tg_xyyyyz_xz_s_0_0_0, tg_xyyyyz_yy_s_0_0_0, tg_xyyyyz_yz_s_0_0_0, tg_xyyyyz_zz_s_0_0_0, tg_xyyyzz_xx_s_0_0_0, tg_xyyyzz_xy_s_0_0_0, tg_xyyyzz_xz_s_0_0_0, tg_xyyyzz_yy_s_0_0_0, tg_xyyyzz_yz_s_0_0_0, tg_xyyyzz_zz_s_0_0_0, tg_xyyzz_xx_s_0_0_1, tg_xyyzz_xy_s_0_0_1, tg_xyyzz_xz_s_0_0_1, tg_xyyzz_yy_s_0_0_1, tg_xyyzz_yz_s_0_0_1, tg_xyyzz_zz_s_0_0_1, tg_xyyzzz_xx_s_0_0_0, tg_xyyzzz_xy_s_0_0_0, tg_xyyzzz_xz_s_0_0_0, tg_xyyzzz_yy_s_0_0_0, tg_xyyzzz_yz_s_0_0_0, tg_xyyzzz_zz_s_0_0_0, tg_xyzzzz_xx_s_0_0_0, tg_xyzzzz_xy_s_0_0_0, tg_xyzzzz_xz_s_0_0_0, tg_xyzzzz_yy_s_0_0_0, tg_xyzzzz_yz_s_0_0_0, tg_xyzzzz_zz_s_0_0_0, tg_xzzz_xx_s_0_0_1, tg_xzzz_xy_s_0_0_1, tg_xzzz_xz_s_0_0_1, tg_xzzz_yy_s_0_0_1, tg_xzzz_yz_s_0_0_1, tg_xzzz_zz_s_0_0_1, tg_xzzzz_xx_s_0_0_1, tg_xzzzz_xy_s_0_0_1, tg_xzzzz_xz_s_0_0_1, tg_xzzzz_yy_s_0_0_1, tg_xzzzz_yz_s_0_0_1, tg_xzzzz_zz_s_0_0_1, tg_xzzzzz_xx_s_0_0_0, tg_xzzzzz_xy_s_0_0_0, tg_xzzzzz_xz_s_0_0_0, tg_xzzzzz_yy_s_0_0_0, tg_xzzzzz_yz_s_0_0_0, tg_xzzzzz_zz_s_0_0_0, tg_yyyy_xx_s_0_0_1, tg_yyyy_xy_s_0_0_1, tg_yyyy_xz_s_0_0_1, tg_yyyy_yy_s_0_0_1, tg_yyyy_yz_s_0_0_1, tg_yyyy_zz_s_0_0_1, tg_yyyyy_xx_s_0_0_1, tg_yyyyy_xy_s_0_0_1, tg_yyyyy_xz_s_0_0_1, tg_yyyyy_yy_s_0_0_1, tg_yyyyy_yz_s_0_0_1, tg_yyyyy_zz_s_0_0_1, tg_yyyyyy_xx_s_0_0_0, tg_yyyyyy_xy_s_0_0_0, tg_yyyyyy_xz_s_0_0_0, tg_yyyyyy_yy_s_0_0_0, tg_yyyyyy_yz_s_0_0_0, tg_yyyyyy_zz_s_0_0_0, tg_yyyyyz_xx_s_0_0_0, tg_yyyyyz_xy_s_0_0_0, tg_yyyyyz_xz_s_0_0_0, tg_yyyyyz_yy_s_0_0_0, tg_yyyyyz_yz_s_0_0_0, tg_yyyyyz_zz_s_0_0_0, tg_yyyyz_xx_s_0_0_1, tg_yyyyz_xy_s_0_0_1, tg_yyyyz_xz_s_0_0_1, tg_yyyyz_yy_s_0_0_1, tg_yyyyz_yz_s_0_0_1, tg_yyyyz_zz_s_0_0_1, tg_yyyyzz_xx_s_0_0_0, tg_yyyyzz_xy_s_0_0_0, tg_yyyyzz_xz_s_0_0_0, tg_yyyyzz_yy_s_0_0_0, tg_yyyyzz_yz_s_0_0_0, tg_yyyyzz_zz_s_0_0_0, tg_yyyzz_xx_s_0_0_1, tg_yyyzz_xy_s_0_0_1, tg_yyyzz_xz_s_0_0_1, tg_yyyzz_yy_s_0_0_1, tg_yyyzz_yz_s_0_0_1, tg_yyyzz_zz_s_0_0_1, tg_yyyzzz_xx_s_0_0_0, tg_yyyzzz_xy_s_0_0_0, tg_yyyzzz_xz_s_0_0_0, tg_yyyzzz_yy_s_0_0_0, tg_yyyzzz_yz_s_0_0_0, tg_yyyzzz_zz_s_0_0_0, tg_yyzz_xx_s_0_0_1, tg_yyzz_xy_s_0_0_1, tg_yyzz_xz_s_0_0_1, tg_yyzz_yy_s_0_0_1, tg_yyzz_yz_s_0_0_1, tg_yyzz_zz_s_0_0_1, tg_yyzzz_xx_s_0_0_1, tg_yyzzz_xy_s_0_0_1, tg_yyzzz_xz_s_0_0_1, tg_yyzzz_yy_s_0_0_1, tg_yyzzz_yz_s_0_0_1, tg_yyzzz_zz_s_0_0_1, tg_yyzzzz_xx_s_0_0_0, tg_yyzzzz_xy_s_0_0_0, tg_yyzzzz_xz_s_0_0_0, tg_yyzzzz_yy_s_0_0_0, tg_yyzzzz_yz_s_0_0_0, tg_yyzzzz_zz_s_0_0_0, tg_yzzz_xx_s_0_0_1, tg_yzzz_xy_s_0_0_1, tg_yzzz_xz_s_0_0_1, tg_yzzz_yy_s_0_0_1, tg_yzzz_yz_s_0_0_1, tg_yzzz_zz_s_0_0_1, tg_yzzzz_xx_s_0_0_1, tg_yzzzz_xy_s_0_0_1, tg_yzzzz_xz_s_0_0_1, tg_yzzzz_yy_s_0_0_1, tg_yzzzz_yz_s_0_0_1, tg_yzzzz_zz_s_0_0_1, tg_yzzzzz_xx_s_0_0_0, tg_yzzzzz_xy_s_0_0_0, tg_yzzzzz_xz_s_0_0_0, tg_yzzzzz_yy_s_0_0_0, tg_yzzzzz_yz_s_0_0_0, tg_yzzzzz_zz_s_0_0_0, tg_zzzz_xx_s_0_0_1, tg_zzzz_xy_s_0_0_1, tg_zzzz_xz_s_0_0_1, tg_zzzz_yy_s_0_0_1, tg_zzzz_yz_s_0_0_1, tg_zzzz_zz_s_0_0_1, tg_zzzzz_xx_s_0_0_1, tg_zzzzz_xy_s_0_0_1, tg_zzzzz_xz_s_0_0_1, tg_zzzzz_yy_s_0_0_1, tg_zzzzz_yz_s_0_0_1, tg_zzzzz_zz_s_0_0_1, tg_zzzzzz_xx_s_0_0_0, tg_zzzzzz_xy_s_0_0_0, tg_zzzzzz_xz_s_0_0_0, tg_zzzzzz_yy_s_0_0_0, tg_zzzzzz_yz_s_0_0_0, tg_zzzzzz_zz_s_0_0_0  : 64)
        for (size_t i = 0; i < nelems; i++)
        {
            const double fai_0 = 1.0 / a_exp;

            tg_xxxxxx_xx_s_0_0_0[i] += 5.0 / 2.0 * tg_xxxx_xx_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxxxx_xx_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxxxx_xy_s_0_0_0[i] += 5.0 / 2.0 * tg_xxxx_xy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxxxx_xy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxxxx_xz_s_0_0_0[i] += 5.0 / 2.0 * tg_xxxx_xz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxxxx_xz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxxxx_yy_s_0_0_0[i] += 5.0 / 2.0 * tg_xxxx_yy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxxxx_yy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxxxx_yz_s_0_0_0[i] += 5.0 / 2.0 * tg_xxxx_yz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxxxx_yz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxxxx_zz_s_0_0_0[i] += 5.0 / 2.0 * tg_xxxx_zz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxxxx_zz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxxxy_xx_s_0_0_0[i] += tg_xxxxx_xx_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxxxy_xy_s_0_0_0[i] += tg_xxxxx_xy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxxxy_xz_s_0_0_0[i] += tg_xxxxx_xz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxxxy_yy_s_0_0_0[i] += tg_xxxxx_yy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxxxy_yz_s_0_0_0[i] += tg_xxxxx_yz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxxxy_zz_s_0_0_0[i] += tg_xxxxx_zz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxxxz_xx_s_0_0_0[i] += tg_xxxxx_xx_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxxxz_xy_s_0_0_0[i] += tg_xxxxx_xy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxxxz_xz_s_0_0_0[i] += tg_xxxxx_xz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxxxz_yy_s_0_0_0[i] += tg_xxxxx_yy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxxxz_yz_s_0_0_0[i] += tg_xxxxx_yz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxxxz_zz_s_0_0_0[i] += tg_xxxxx_zz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxxyy_xx_s_0_0_0[i] += 3.0 / 2.0 * tg_xxyy_xx_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxxyy_xx_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxxyy_xy_s_0_0_0[i] += 3.0 / 2.0 * tg_xxyy_xy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxxyy_xy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxxyy_xz_s_0_0_0[i] += 3.0 / 2.0 * tg_xxyy_xz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxxyy_xz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxxyy_yy_s_0_0_0[i] += 3.0 / 2.0 * tg_xxyy_yy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxxyy_yy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxxyy_yz_s_0_0_0[i] += 3.0 / 2.0 * tg_xxyy_yz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxxyy_yz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxxyy_zz_s_0_0_0[i] += 3.0 / 2.0 * tg_xxyy_zz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxxyy_zz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxxyz_xx_s_0_0_0[i] += tg_xxxxz_xx_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxxyz_xy_s_0_0_0[i] += tg_xxxxz_xy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxxyz_xz_s_0_0_0[i] += tg_xxxxz_xz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxxyz_yy_s_0_0_0[i] += tg_xxxxz_yy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxxyz_yz_s_0_0_0[i] += tg_xxxxz_yz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxxyz_zz_s_0_0_0[i] += tg_xxxxz_zz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxxzz_xx_s_0_0_0[i] += 3.0 / 2.0 * tg_xxzz_xx_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxxzz_xx_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxxzz_xy_s_0_0_0[i] += 3.0 / 2.0 * tg_xxzz_xy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxxzz_xy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxxzz_xz_s_0_0_0[i] += 3.0 / 2.0 * tg_xxzz_xz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxxzz_xz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxxzz_yy_s_0_0_0[i] += 3.0 / 2.0 * tg_xxzz_yy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxxzz_yy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxxzz_yz_s_0_0_0[i] += 3.0 / 2.0 * tg_xxzz_yz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxxzz_yz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxxzz_zz_s_0_0_0[i] += 3.0 / 2.0 * tg_xxzz_zz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxxzz_zz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxyyy_xx_s_0_0_0[i] += tg_xyyy_xx_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxyyy_xx_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxyyy_xy_s_0_0_0[i] += tg_xyyy_xy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxyyy_xy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxyyy_xz_s_0_0_0[i] += tg_xyyy_xz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxyyy_xz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxyyy_yy_s_0_0_0[i] += tg_xyyy_yy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxyyy_yy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxyyy_yz_s_0_0_0[i] += tg_xyyy_yz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxyyy_yz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxyyy_zz_s_0_0_0[i] += tg_xyyy_zz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxyyy_zz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxyyz_xx_s_0_0_0[i] += tg_xxxyy_xx_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxyyz_xy_s_0_0_0[i] += tg_xxxyy_xy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxyyz_xz_s_0_0_0[i] += tg_xxxyy_xz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxyyz_yy_s_0_0_0[i] += tg_xxxyy_yy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxyyz_yz_s_0_0_0[i] += tg_xxxyy_yz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxyyz_zz_s_0_0_0[i] += tg_xxxyy_zz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxyzz_xx_s_0_0_0[i] += tg_xxxzz_xx_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxyzz_xy_s_0_0_0[i] += tg_xxxzz_xy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxyzz_xz_s_0_0_0[i] += tg_xxxzz_xz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxyzz_yy_s_0_0_0[i] += tg_xxxzz_yy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxyzz_yz_s_0_0_0[i] += tg_xxxzz_yz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxyzz_zz_s_0_0_0[i] += tg_xxxzz_zz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxzzz_xx_s_0_0_0[i] += tg_xzzz_xx_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxzzz_xx_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxzzz_xy_s_0_0_0[i] += tg_xzzz_xy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxzzz_xy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxzzz_xz_s_0_0_0[i] += tg_xzzz_xz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxzzz_xz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxzzz_yy_s_0_0_0[i] += tg_xzzz_yy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxzzz_yy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxzzz_yz_s_0_0_0[i] += tg_xzzz_yz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxzzz_yz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxzzz_zz_s_0_0_0[i] += tg_xzzz_zz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxzzz_zz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyyyy_xx_s_0_0_0[i] += 1.0 / 2.0 * tg_yyyy_xx_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyyyy_xx_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyyyy_xy_s_0_0_0[i] += 1.0 / 2.0 * tg_yyyy_xy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyyyy_xy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyyyy_xz_s_0_0_0[i] += 1.0 / 2.0 * tg_yyyy_xz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyyyy_xz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyyyy_yy_s_0_0_0[i] += 1.0 / 2.0 * tg_yyyy_yy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyyyy_yy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyyyy_yz_s_0_0_0[i] += 1.0 / 2.0 * tg_yyyy_yz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyyyy_yz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyyyy_zz_s_0_0_0[i] += 1.0 / 2.0 * tg_yyyy_zz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyyyy_zz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyyyz_xx_s_0_0_0[i] += tg_xxyyy_xx_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxyyyz_xy_s_0_0_0[i] += tg_xxyyy_xy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxyyyz_xz_s_0_0_0[i] += tg_xxyyy_xz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxyyyz_yy_s_0_0_0[i] += tg_xxyyy_yy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxyyyz_yz_s_0_0_0[i] += tg_xxyyy_yz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxyyyz_zz_s_0_0_0[i] += tg_xxyyy_zz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxyyzz_xx_s_0_0_0[i] += 1.0 / 2.0 * tg_yyzz_xx_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyyzz_xx_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyyzz_xy_s_0_0_0[i] += 1.0 / 2.0 * tg_yyzz_xy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyyzz_xy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyyzz_xz_s_0_0_0[i] += 1.0 / 2.0 * tg_yyzz_xz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyyzz_xz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyyzz_yy_s_0_0_0[i] += 1.0 / 2.0 * tg_yyzz_yy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyyzz_yy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyyzz_yz_s_0_0_0[i] += 1.0 / 2.0 * tg_yyzz_yz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyyzz_yz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyyzz_zz_s_0_0_0[i] += 1.0 / 2.0 * tg_yyzz_zz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyyzz_zz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyzzz_xx_s_0_0_0[i] += tg_xxzzz_xx_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxyzzz_xy_s_0_0_0[i] += tg_xxzzz_xy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxyzzz_xz_s_0_0_0[i] += tg_xxzzz_xz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxyzzz_yy_s_0_0_0[i] += tg_xxzzz_yy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxyzzz_yz_s_0_0_0[i] += tg_xxzzz_yz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxyzzz_zz_s_0_0_0[i] += tg_xxzzz_zz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxzzzz_xx_s_0_0_0[i] += 1.0 / 2.0 * tg_zzzz_xx_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzzzz_xx_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxzzzz_xy_s_0_0_0[i] += 1.0 / 2.0 * tg_zzzz_xy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzzzz_xy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxzzzz_xz_s_0_0_0[i] += 1.0 / 2.0 * tg_zzzz_xz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzzzz_xz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxzzzz_yy_s_0_0_0[i] += 1.0 / 2.0 * tg_zzzz_yy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzzzz_yy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxzzzz_yz_s_0_0_0[i] += 1.0 / 2.0 * tg_zzzz_yz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzzzz_yz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxzzzz_zz_s_0_0_0[i] += 1.0 / 2.0 * tg_zzzz_zz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzzzz_zz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyyy_xx_s_0_0_0[i] += tg_yyyyy_xx_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyyy_xy_s_0_0_0[i] += tg_yyyyy_xy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyyy_xz_s_0_0_0[i] += tg_yyyyy_xz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyyy_yy_s_0_0_0[i] += tg_yyyyy_yy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyyy_yz_s_0_0_0[i] += tg_yyyyy_yz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyyy_zz_s_0_0_0[i] += tg_yyyyy_zz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyyz_xx_s_0_0_0[i] += tg_yyyyz_xx_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyyz_xy_s_0_0_0[i] += tg_yyyyz_xy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyyz_xz_s_0_0_0[i] += tg_yyyyz_xz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyyz_yy_s_0_0_0[i] += tg_yyyyz_yy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyyz_yz_s_0_0_0[i] += tg_yyyyz_yz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyyz_zz_s_0_0_0[i] += tg_yyyyz_zz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyzz_xx_s_0_0_0[i] += tg_yyyzz_xx_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyzz_xy_s_0_0_0[i] += tg_yyyzz_xy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyzz_xz_s_0_0_0[i] += tg_yyyzz_xz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyzz_yy_s_0_0_0[i] += tg_yyyzz_yy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyzz_yz_s_0_0_0[i] += tg_yyyzz_yz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyzz_zz_s_0_0_0[i] += tg_yyyzz_zz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyzzz_xx_s_0_0_0[i] += tg_yyzzz_xx_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyzzz_xy_s_0_0_0[i] += tg_yyzzz_xy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyzzz_xz_s_0_0_0[i] += tg_yyzzz_xz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyzzz_yy_s_0_0_0[i] += tg_yyzzz_yy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyzzz_yz_s_0_0_0[i] += tg_yyzzz_yz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyzzz_zz_s_0_0_0[i] += tg_yyzzz_zz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzzzz_xx_s_0_0_0[i] += tg_yzzzz_xx_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzzzz_xy_s_0_0_0[i] += tg_yzzzz_xy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzzzz_xz_s_0_0_0[i] += tg_yzzzz_xz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzzzz_yy_s_0_0_0[i] += tg_yzzzz_yy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzzzz_yz_s_0_0_0[i] += tg_yzzzz_yz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzzzz_zz_s_0_0_0[i] += tg_yzzzz_zz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzzzz_xx_s_0_0_0[i] += tg_zzzzz_xx_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzzzz_xy_s_0_0_0[i] += tg_zzzzz_xy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzzzz_xz_s_0_0_0[i] += tg_zzzzz_xz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzzzz_yy_s_0_0_0[i] += tg_zzzzz_yy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzzzz_yz_s_0_0_0[i] += tg_zzzzz_yz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzzzz_zz_s_0_0_0[i] += tg_zzzzz_zz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_yyyyyy_xx_s_0_0_0[i] += 5.0 / 2.0 * tg_yyyy_xx_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyyyy_xx_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyyyy_xy_s_0_0_0[i] += 5.0 / 2.0 * tg_yyyy_xy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyyyy_xy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyyyy_xz_s_0_0_0[i] += 5.0 / 2.0 * tg_yyyy_xz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyyyy_xz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyyyy_yy_s_0_0_0[i] += 5.0 / 2.0 * tg_yyyy_yy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyyyy_yy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyyyy_yz_s_0_0_0[i] += 5.0 / 2.0 * tg_yyyy_yz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyyyy_yz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyyyy_zz_s_0_0_0[i] += 5.0 / 2.0 * tg_yyyy_zz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyyyy_zz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyyyz_xx_s_0_0_0[i] += tg_yyyyy_xx_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyyyz_xy_s_0_0_0[i] += tg_yyyyy_xy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyyyz_xz_s_0_0_0[i] += tg_yyyyy_xz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyyyz_yy_s_0_0_0[i] += tg_yyyyy_yy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyyyz_yz_s_0_0_0[i] += tg_yyyyy_yz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyyyz_zz_s_0_0_0[i] += tg_yyyyy_zz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyyzz_xx_s_0_0_0[i] += 3.0 / 2.0 * tg_yyzz_xx_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyyzz_xx_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyyzz_xy_s_0_0_0[i] += 3.0 / 2.0 * tg_yyzz_xy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyyzz_xy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyyzz_xz_s_0_0_0[i] += 3.0 / 2.0 * tg_yyzz_xz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyyzz_xz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyyzz_yy_s_0_0_0[i] += 3.0 / 2.0 * tg_yyzz_yy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyyzz_yy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyyzz_yz_s_0_0_0[i] += 3.0 / 2.0 * tg_yyzz_yz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyyzz_yz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyyzz_zz_s_0_0_0[i] += 3.0 / 2.0 * tg_yyzz_zz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyyzz_zz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyzzz_xx_s_0_0_0[i] += tg_yzzz_xx_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyzzz_xx_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyzzz_xy_s_0_0_0[i] += tg_yzzz_xy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyzzz_xy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyzzz_xz_s_0_0_0[i] += tg_yzzz_xz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyzzz_xz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyzzz_yy_s_0_0_0[i] += tg_yzzz_yy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyzzz_yy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyzzz_yz_s_0_0_0[i] += tg_yzzz_yz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyzzz_yz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyzzz_zz_s_0_0_0[i] += tg_yzzz_zz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyzzz_zz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzzzz_xx_s_0_0_0[i] += 1.0 / 2.0 * tg_zzzz_xx_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzzzz_xx_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzzzz_xy_s_0_0_0[i] += 1.0 / 2.0 * tg_zzzz_xy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzzzz_xy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzzzz_xz_s_0_0_0[i] += 1.0 / 2.0 * tg_zzzz_xz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzzzz_xz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzzzz_yy_s_0_0_0[i] += 1.0 / 2.0 * tg_zzzz_yy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzzzz_yy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzzzz_yz_s_0_0_0[i] += 1.0 / 2.0 * tg_zzzz_yz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzzzz_yz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzzzz_zz_s_0_0_0[i] += 1.0 / 2.0 * tg_zzzz_zz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzzzz_zz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzzzz_xx_s_0_0_0[i] += tg_zzzzz_xx_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzzzz_xy_s_0_0_0[i] += tg_zzzzz_xy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzzzz_xz_s_0_0_0[i] += tg_zzzzz_xz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzzzz_yy_s_0_0_0[i] += tg_zzzzz_yy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzzzz_yz_s_0_0_0[i] += tg_zzzzz_yz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzzzz_zz_s_0_0_0[i] += tg_zzzzz_zz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_zzzzzz_xx_s_0_0_0[i] += 5.0 / 2.0 * tg_zzzz_xx_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzzzz_xx_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzzzzz_xy_s_0_0_0[i] += 5.0 / 2.0 * tg_zzzz_xy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzzzz_xy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzzzzz_xz_s_0_0_0[i] += 5.0 / 2.0 * tg_zzzz_xz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzzzz_xz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzzzzz_yy_s_0_0_0[i] += 5.0 / 2.0 * tg_zzzz_yy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzzzz_yy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzzzzz_yz_s_0_0_0[i] += 5.0 / 2.0 * tg_zzzz_yz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzzzz_yz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzzzzz_zz_s_0_0_0[i] += 5.0 / 2.0 * tg_zzzz_zz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzzzz_zz_s_0_0_1[i] * fai_0 * a_z * fp_0;
        }
    }
}

} // t2pecp namespace

